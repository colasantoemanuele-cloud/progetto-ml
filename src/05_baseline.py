#!/usr/bin/env python3
"""
05_baseline.py - Modelli di riferimento e divario fra split casuale e split per paziente.

Scopo
-----
Addestra due modelli di riferimento (DummyClassifier stratificato e LogisticRegression
L2 con pesi di classe bilanciati) in una Pipeline che contiene tutte le
trasformazioni dipendenti dai dati: `LogCPM` (log1p di CP10k, non CPM, nonostante il
nome; non stima parametri) e `HVGSelector` (geni piu' dispersi) da `scpipe.py`, e
`StandardScaler(with_mean=False)` di scikit-learn. Le trasformazioni sono stimate solo sul
training di ciascun fold.

Ogni task e' valutato due volte: con lo split per paziente (`fold_group_5`, oppure
`fold_lopo` per task1a_paired e rare) e con lo split casuale per cellula
(`fold_random_5`, controllo negativo). La differenza misura quanto la correlazione
fra cellule dello stesso paziente gonfia la stima. Il macro-F1 di ogni fold e'
calcolato sulle sole classi presenti nel fold di test; l'F1 per classe e' NaN
quando la classe manca dal fold di test.

Input
-----
  --h5ad    data/interim/sc_raw.h5ad (conteggi grezzi)
  --splits  data/interim/splits_sc.csv (da 04_splits.py)

Output
------
  results/tables/baseline_{tag}.csv            metriche per task x schema x modello x fold (e F1 per classe)
  results/tables/baseline_summary_{tag}.csv    media e deviazione standard sui fold di macro-F1, BA, MCC
  results/tables/baseline_gap_{tag}.csv        macro-F1 medio per schema, divario casuale - paziente,
                                               rapporto fra le deviazioni standard (logreg_l2)
  results/figures/gap_{tag}_{task}.png         F1 per classe, split casuale contro split per paziente

Argomenti
---------
  --h5ad, --splits (obbligatori), --tag (default sc), --tasks (default i 6 task del README),
  --n-hvg (default 2000), --quick (sottocampiona il training a 20.000 cellule)

Uso
---
    python src/05_baseline.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv \\
        --tasks task2_hier1 task2_hier3 task2_v1 task2_v2 task1a_paired task1a_all
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import HVGSelector, LogCPM, macro_f1  # noqa: E402

from sklearn.dummy import DummyClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import (  # noqa: E402
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
)
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0

# Schema di partizionamento primario associato a ciascun task
PRIMARY = {"task1a_paired": "fold_lopo", "rare": "fold_lopo"}
DEFAULT_PRIMARY = "fold_group_5"
ALL_TASKS = [
    "task2_hier1",
    "task2_hier3",
    "task2_v1",
    "task2_v2",
    "task1a_paired",
    "task1a_all",
]


def make_pipeline(model: str, n_classes: int) -> Pipeline:
    """Costruisce la Pipeline end-to-end (LogCPM -> HVGSelector -> StandardScaler -> Classificatore)."""
    if model == "dummy":
        clf = DummyClassifier(strategy="stratified", random_state=SEED)
    elif model == "logreg_l2":
        # Regolarizzazione L2 standard con bilanciamento delle classi tramite pesi inversi
        clf = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=2000,
            class_weight="balanced",
            random_state=SEED,
        )
    else:
        raise ValueError(model)
    return Pipeline(
        [
            ("logcpm", LogCPM()),
            ("hvg", HVGSelector(n_top=2000)),
            # StandardScaler con with_mean=False evita la sottrazione della media preservando la struttura CSR
            ("scale", StandardScaler(with_mean=False)),
            ("clf", clf),
        ]
    )


def evaluate(y_true, y_pred, classes) -> list[dict]:
    """Balanced accuracy, macro-F1 (classi presenti nel test), MCC e F1 per classe (NaN se assente dal test)."""
    rows = [
        {
            "metric": "balanced_accuracy",
            "class": "_overall",
            "value": balanced_accuracy_score(y_true, y_pred),
        },
        {
            "metric": "macro_f1",
            "class": "_overall",
            "value": macro_f1(y_true, y_pred),
        },
        {
            "metric": "mcc",
            "class": "_overall",
            "value": matthews_corrcoef(y_true, y_pred),
        },
    ]
    per = f1_score(y_true, y_pred, average=None, labels=classes, zero_division=0)
    present = np.isin(classes, np.unique(y_true))
    per = np.where(present, per, np.nan)
    rows += [{"metric": "f1", "class": c, "value": v} for c, v in zip(classes, per)]
    return rows


def run_task(
    adata, sp_df: pd.DataFrame, task: str, quick: bool, n_hvg: int
) -> pd.DataFrame:
    """Esegue la cross-validation completa per un singolo task su partizione onesta e casuale."""
    sub = sp_df[sp_df.task == task]
    cells = sub.cell_id.values
    idx = adata.obs_names.get_indexer(cells)
    if len(idx) == 0 or (idx < 0).any():
        raise SystemExit(
            f"[ERRORE] task {task}: cell_id assenti dall'h5ad o task vuoto"
        )

    # Estrazione della sottomatrice sparsa delle cellule appartenenti al task corrente
    X = adata.X[idx]
    y = sub.y.values
    classes = np.array(sorted(np.unique(y)))
    primary = PRIMARY.get(task, DEFAULT_PRIMARY)

    print(
        f"\n{'=' * 74}\n{task}: {X.shape[0]:,} cellule x {X.shape[1]:,} geni | "
        f"{len(classes)} classi | primaria: {primary}\n{'=' * 74}"
    )

    rows = []
    for scheme in (primary, "fold_random_5"):
        folds = sub[scheme].values
        for model in ("dummy", "logreg_l2"):
            t0 = time.time()
            for f in np.unique(folds):
                tr, te = folds != f, folds == f
                # Sottocampionamento opzionale del solo training set per esecuzioni diagnostiche rapide
                if quick and tr.sum() > 20000:
                    rng = np.random.default_rng(SEED + f)
                    keep = rng.choice(np.where(tr)[0], 20000, replace=False)
                    tr = np.zeros_like(tr)
                    tr[keep] = True

                if len(np.unique(y[tr])) < 2:
                    raise SystemExit(
                        f"[ERRORE] {task}, {scheme}, fold {f}: il training ha una sola classe"
                    )
                pipe = make_pipeline(model, len(classes))
                pipe.set_params(hvg__n_top=n_hvg)
                pipe.fit(X[tr], y[tr])
                pred = pipe.predict(X[te])

                for r in evaluate(y[te], pred, classes):
                    rows.append(
                        {
                            "task": task,
                            "scheme": scheme,
                            "model": model,
                            "fold": int(f),
                            **r,
                        }
                    )

            df = pd.DataFrame(rows)
            m = df[
                (df.scheme == scheme) & (df.model == model) & (df.metric == "macro_f1")
            ].value
            print(
                f"  {scheme:14s} {model:10s} macro-F1 = {m.mean():.3f} ± {m.std():.3f}"
                f"   [{time.time() - t0:.0f}s]"
            )

    return pd.DataFrame(rows)


def plot_gap(res: pd.DataFrame, tag: str) -> None:
    """Genera i grafici a barre orizzontali del divario di F1 per classe (split casuale vs per paziente)."""
    d = res[
        (res.model == "logreg_l2") & (res.metric == "f1") & (res["class"] != "_overall")
    ]
    if d.empty:
        return
    for task, g in d.groupby("task"):
        piv = g.pivot_table(
            index="class", columns="scheme", values="value", aggfunc="mean"
        )
        err = g.pivot_table(
            index="class", columns="scheme", values="value", aggfunc="std"
        )
        if piv.shape[1] < 2:
            continue
        rnd = "fold_random_5"
        grp = [c for c in piv.columns if c != rnd][0]
        piv = piv[[rnd, grp]].sort_values(grp)

        fig, ax = plt.subplots(figsize=(8, 0.55 * len(piv) + 2))
        ypos = np.arange(len(piv))
        ax.barh(
            ypos - 0.2,
            piv[rnd],
            0.38,
            xerr=err[rnd].reindex(piv.index),
            label="split casuale (leakage)",
            color="#d95f02",
            alpha=0.85,
        )
        ax.barh(
            ypos + 0.2,
            piv[grp],
            0.38,
            xerr=err[grp].reindex(piv.index),
            label="split per paziente (onesto)",
            color="#1b9e77",
            alpha=0.85,
        )
        ax.set_yticks(ypos, piv.index)
        ax.set_xlabel("F1 per classe")
        ax.set_xlim(0, 1)
        ax.set_title(f"{task} - divario tra split casuale e split per paziente")
        ax.legend(loc="lower right")
        for i, c in enumerate(piv.index):
            ax.text(
                1.01,
                i,
                f"Δ={piv.loc[c, rnd] - piv.loc[c, grp]:+.3f}",
                va="center",
                fontsize=8,
            )
        fig.tight_layout()
        fig.savefig(FIG / f"gap_{tag}_{task}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"[ok] figura: results/figures/gap_{tag}_{task}.png")


def main() -> None:
    """Punto di ingresso principale per l'addestramento e la valutazione dei modelli baseline."""
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Modelli di riferimento: split casuale contro split per paziente.",
    )
    ap.add_argument("--h5ad", type=Path, required=True, help="AnnData scRNA-seq")
    ap.add_argument(
        "--splits", type=Path, required=True, help="split prodotti da 04_splits.py"
    )
    ap.add_argument("--tag", default="sc", help="suffisso dei file di output")
    ap.add_argument("--tasks", nargs="+", default=ALL_TASKS, help="task da valutare")
    ap.add_argument(
        "--n-hvg",
        type=int,
        default=2000,
        help="numero di geni piu' dispersi (HVGSelector)",
    )
    ap.add_argument(
        "--quick", action="store_true", help="sottocampiona il TRAIN a 20k cellule"
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    # Caricamento in memoria dell'oggetto AnnData e conversione esplicita in formato CSR
    print("[..] carico h5ad in memoria (serve X) ...")
    adata = ad.read_h5ad(args.h5ad)
    adata.X = sp.csr_matrix(adata.X)
    print(f"[ok] {adata.shape}, {adata.X.nnz:,} nnz")

    sp_df = pd.read_csv(args.splits)
    missing = set(args.tasks) - set(sp_df.task.unique())
    if missing:
        raise SystemExit(f"task assenti da {args.splits}: {missing}")

    res = pd.concat(
        [run_task(adata, sp_df, t, args.quick, args.n_hvg) for t in args.tasks]
    )
    out = TAB / f"baseline_{args.tag}.csv"
    res.to_csv(out, index=False)
    print(f"\n[ok] scritto {out}")

    summ = (
        res[res.metric.isin(["macro_f1", "balanced_accuracy", "mcc"])]
        .groupby(["task", "scheme", "model", "metric"])
        .value.agg(["mean", "std"])
    )
    print("\n" + summ.round(3).to_string())
    summ.reset_index().to_csv(TAB / f"baseline_summary_{args.tag}.csv", index=False)

    # Calcolo del differenziale medio di Macro-F1 tra valutazione con leakage e valutazione per paziente
    gap = (
        res[(res.model == "logreg_l2") & (res.metric == "macro_f1")]
        .groupby(["task", "scheme"])
        .value.mean()
        .unstack()
    )
    gap["gap_random_minus_honest"] = gap.get("fold_random_5") - gap.drop(
        columns="fold_random_5"
    ).max(axis=1)
    print("\n>>> RIEPILOGO DIVARIO DI GENERALIZZAZIONE:\n" + gap.round(3).to_string())
    # rapporto fra le deviazioni standard sui fold (split per paziente / split casuale)
    sd = (
        res[(res.model == "logreg_l2") & (res.metric == "macro_f1")]
        .groupby(["task", "scheme"])
        .value.std()
        .unstack()
    )
    gap["sd_ratio_honest_over_random"] = sd.drop(columns="fold_random_5").max(
        axis=1
    ) / sd.get("fold_random_5")
    gap.reset_index().to_csv(TAB / f"baseline_gap_{args.tag}.csv", index=False)

    plot_gap(res, args.tag)


if __name__ == "__main__":
    main()

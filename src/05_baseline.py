#!/usr/bin/env python3
"""
05_baseline.py - Modelli lineari di riferimento e quantificazione dell'inflazione da leakage.

Questo modulo addestra i modelli di riferimento (DummyClassifier stratificato e
LogisticRegression regolarizzata L2 con pesatura inversamente proporzionale alle
frequenze di classe) incapsulando tutte le trasformazioni dipendenti dai dati
(normalizzazione LogCPM, selezione dei geni ad alta variabilita' HVGSelector e
standardizzazione StandardScaler senza centramento) all'interno di una Pipeline
Scikit-learn. In questo modo i parametri di dispersione genica e varianza vengono
stimati esclusivamente sul fold di training, prevenendo il data leakage da
feature selection globale.

Per ogni task, i modelli vengono valutati sia sulla partizione stratificata per
paziente (fold_group_5 / fold_lopo, stima di generalizzazione inter-paziente) sia
sulla partizione casuale a livello di cellula (fold_random_5, controllo negativo),
quantificando il differenziale di performance indotto dalla correlazione
intra-paziente.

Uso:
    python src/05_baseline.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0

# Schema di partizionamento primario associato a ciascun task
PRIMARY = {"task1a_paired": "fold_lopo", "rare": "fold_lopo"}
DEFAULT_PRIMARY = "fold_group_5"


class LogCPM(BaseEstimator, TransformerMixin):
    """Normalizzazione per-cellula log1p(CPM/100): operatore stateless inserito nella Pipeline."""

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        # Calcolo della profondita' di sequenziamento totale per cellula e riscalamento CSR in-place
        counts = np.asarray(X.sum(axis=1)).ravel()
        counts[counts == 0] = 1.0
        X = sp.csr_matrix(X.multiply(self.target_sum / counts[:, None]))
        X.data = np.log1p(X.data)
        return X


class HVGSelector(BaseEstimator, TransformerMixin):
    """Selezione dei top-k geni per indice di dispersione (varianza/media) stimato sul solo training set."""

    def __init__(self, n_top: int = 2000, min_mean: float = 1e-3):
        self.n_top = n_top
        self.min_mean = min_mean

    def fit(self, X, y=None):  # noqa: ARG002
        # Calcolo vettorizzato di media e momento secondo su matrice sparsa CSR senza densificazione
        mean = np.asarray(X.mean(axis=0)).ravel()
        sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = np.maximum(sq - mean**2, 0.0)
        disp = np.where(mean > self.min_mean, var / (mean + 1e-12), -np.inf)
        k = min(self.n_top, int(np.isfinite(disp).sum()))
        self.idx_ = np.sort(np.argsort(-disp)[:k])
        return self

    def transform(self, X):
        return X[:, self.idx_]


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
            n_jobs=-1 if n_classes == 2 else None,
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
    """Calcola le metriche globali (Balanced Accuracy, Macro-F1, MCC) e gli score F1 classe-specifici."""
    rows = [
        {
            "metric": "balanced_accuracy",
            "class": "_overall",
            "value": balanced_accuracy_score(y_true, y_pred),
        },
        {
            "metric": "macro_f1",
            "class": "_overall",
            "value": f1_score(
                y_true, y_pred, average="macro", labels=classes, zero_division=0
            ),
        },
        {
            "metric": "mcc",
            "class": "_overall",
            "value": matthews_corrcoef(y_true, y_pred),
        },
    ]
    per = f1_score(y_true, y_pred, average=None, labels=classes, zero_division=0)
    rows += [{"metric": "f1", "class": c, "value": v} for c, v in zip(classes, per)]
    return rows


def run_task(
    adata, sp_df: pd.DataFrame, task: str, quick: bool, n_hvg: int
) -> pd.DataFrame:
    """Esegue la cross-validation completa per un singolo task su partizione onesta e casuale."""
    sub = sp_df[sp_df.task == task]
    cells = sub.cell_id.values
    idx = adata.obs_names.get_indexer(cells)
    assert (idx >= 0).all(), "cell_id assenti dall'h5ad"

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
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument("--tag", default="sc")
    ap.add_argument(
        "--tasks", nargs="*", default=["task2_hier1", "task2_hier3", "task2_v1"]
    )
    ap.add_argument("--n-hvg", type=int, default=2000)
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
        .round(3)
    )
    print("\n" + summ.to_string())

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

    plot_gap(res, args.tag)


if __name__ == "__main__":
    main()

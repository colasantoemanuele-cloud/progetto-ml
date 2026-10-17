#!/usr/bin/env python3
"""
11_model_zoo.py - Confronto di 13 classificatori con validazione per paziente.

Scopo
-----
Confronta tredici classificatori:
    dummy, naive_bayes, knn, logreg_l2, logreg_l1, elasticnet,
    decision_tree, bagging, random_forest, lightgbm,
    svm_linear, svm_rbf, mlp

Scelte del disegno:
1. Stesso numero di cellule di training per tutti i modelli: ogni fold di training
   e' sottocampionato in modo casuale uniforme (non stratificato) a `--max-train`
   cellule (default 20.000, seme SEED + fold), dopo lo split per paziente. Il
   confronto misura quindi gli algoritmi a parita' di dati, non il loro costo.
2. Rappresentazione: tutti i modelli partono da `LogCPM` (log1p CP10k) e dai 2.000
   geni piu' dispersi (`HVGSelector`), stimati sul training del fold.
   - `knn`, `svm_rbf` e `mlp` lavorano su `TruncatedSVD` a 50 componenti seguita da
     `StandardScaler` (stimati sul training): per knn e svm_rbf, basati su distanze
     e kernel, si riduce la concentrazione delle distanze; per l'MLP, che non e'
     basato su distanze, la proiezione densa serve a ridurre dimensione e costo.
   - `naive_bayes` (MultinomialNB) usa direttamente i valori log-normalizzati, non
     negativi, senza scaling (lo scaling non e' necessario e non cambia il segno).
   - gli altri modelli usano la matrice CSR scalata senza centratura.
3. Confronto con il modello migliore: t appaiato corretto di Nadeau-Bengio
   (fattore n_test / n_train) e test di Wilcoxon appaiato sui fold. Con pochi fold la
   potenza e' bassa: un p non significativo NON dimostra che due modelli siano
   equivalenti. La calibrazione di Nadeau-Bengio e' verificata in
   `src/sim/sim_a_calibration.py` con StratifiedKFold a 5 fold di dimensione uguale;
   qui e' applicato a fold LOPO di dimensione diversa con training sottocampionato,
   caso non coperto da quella verifica.
Il macro-F1 di ogni fold e' calcolato sulle sole classi presenti nel fold di test.
Se un modello fallisce in un fold, tutti i suoi fold vengono scartati.
Riproducibilita': `knn`, `svm_rbf` e `mlp` (proiezione densa) sono addestrati con BLAS a
un solo thread: con piu' thread TruncatedSVD e l'MLP cambiano risultato secondo il numero di
thread (l'MLP fino a 0,014 di macro-F1 per fold, knn e svm_rbf meno di 0,001).
Convergenza: per i modelli con solutore iterativo (logreg_l2, logreg_l1, elasticnet,
svm_linear, mlp) ogni fold registra se l'addestramento ha emesso un
`ConvergenceWarning` (colonna `converged`; vuota per gli altri modelli).

Input
-----
  --h5ad    data/interim/sc_raw.h5ad
  --splits  data/interim/splits_sc.csv (da 04_splits.py)

Output
------
  results/tables/model_zoo_{tag}.csv            metriche e convergenza per task x modello x fold
  results/tables/model_zoo_summary_{tag}.csv    media e sd sui fold, differenza da logreg_l2 e
                                                fold non convergenti per task x modello
  results/tables/model_comparison_{task}.csv    differenze dal migliore e p (Nadeau-Bengio, Wilcoxon)
  results/figures/model_zoo_{task}.png

Argomenti
---------
  --h5ad, --splits (obbligatori), --tasks (default task2_hier3 task2_v1; con --from-table tutti i task della tabella),
  --models (default tutti), --max-train (default 20000; 0 = tutte le cellule),
  --scheme (forza fold_group_5 / fold_random_5 / fold_lopo), --tag (default sc),
  --update (sostituisce nella tabella esistente solo i modelli e i task rieseguiti),
  --from-table (non addestra: rigenera riepilogo, confronti e figure dalla tabella per fold)

Uso
---
    python src/11_model_zoo.py --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv --tasks task2_v1 task2_hier3 --scheme fold_lopo
    python src/11_model_zoo.py --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv --tasks task1a_paired --tag task1 --scheme fold_lopo
    # rieseguire un solo modello senza perdere gli altri:
    python src/11_model_zoo.py ... --models lightgbm --update
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import HVGSelector, LogCPM, macro_f1, nadeau_bengio_ttest  # noqa: E402

from sklearn.decomposition import TruncatedSVD  # noqa: E402
from sklearn.dummy import DummyClassifier  # noqa: E402
from sklearn.exceptions import ConvergenceWarning  # noqa: E402
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression, SGDClassifier  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef  # noqa: E402
from sklearn.naive_bayes import MultinomialNB  # noqa: E402
from sklearn.neighbors import KNeighborsClassifier  # noqa: E402
from sklearn.neural_network import MLPClassifier  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from sklearn.svm import SVC, LinearSVC  # noqa: E402
from sklearn.tree import DecisionTreeClassifier  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402
from threadpoolctl import threadpool_limits  # noqa: E402

# Solo l'avviso ripetitivo sui nomi delle feature (LightGBM); i ConvergenceWarning, che in
# scikit-learn derivano da UserWarning, restano visibili e sono registrati per fold
warnings.filterwarnings("ignore", message=".*feature names.*", category=UserWarning)

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_HVG = 2000
N_SVD = 50

# modelli con solutore iterativo, per cui si registra la convergenza
ITERATIVE = {"logreg_l2", "logreg_l1", "elasticnet", "svm_linear", "mlp"}

PRIMARY = {"task1a_paired": "fold_lopo", "rare": "fold_lopo"}
DEFAULT_PRIMARY = "fold_group_5"


def _models() -> dict:
    """Restituisce il dizionario dei costruttori dei modelli e il flag di proiezione SVD densa."""
    return {
        "dummy": (
            lambda: DummyClassifier(strategy="stratified", random_state=SEED),
            False,
        ),
        "naive_bayes": (lambda: MultinomialNB(), False),
        "knn": (lambda: KNeighborsClassifier(n_neighbors=15, n_jobs=-1), True),
        "logreg_l2": (
            lambda: LogisticRegression(
                C=1.0,
                solver="lbfgs",
                max_iter=2000,
                class_weight="balanced",
                random_state=SEED,
            ),
            False,
        ),
        "logreg_l1": (
            lambda: SGDClassifier(
                loss="log_loss",
                penalty="l1",
                alpha=1e-4,
                max_iter=50,
                class_weight="balanced",
                n_jobs=-1,
                random_state=SEED,
            ),
            False,
        ),
        "elasticnet": (
            lambda: SGDClassifier(
                loss="log_loss",
                penalty="elasticnet",
                l1_ratio=0.5,
                alpha=1e-4,
                max_iter=50,
                class_weight="balanced",
                n_jobs=-1,
                random_state=SEED,
            ),
            False,
        ),
        "decision_tree": (
            lambda: DecisionTreeClassifier(
                max_depth=12, class_weight="balanced", random_state=SEED
            ),
            False,
        ),
        "bagging": (
            lambda: BaggingClassifier(
                DecisionTreeClassifier(max_depth=12, random_state=SEED),
                n_estimators=50,
                n_jobs=-1,
                random_state=SEED,
            ),
            False,
        ),
        "random_forest": (
            lambda: RandomForestClassifier(
                n_estimators=300,
                n_jobs=-1,
                class_weight="balanced_subsample",
                random_state=SEED,
            ),
            False,
        ),
        "lightgbm": (_lgbm, False),
        "svm_linear": (
            lambda: LinearSVC(
                C=0.1, class_weight="balanced", max_iter=5000, random_state=SEED
            ),
            False,
        ),
        "svm_rbf": (
            lambda: SVC(
                C=1.0,
                gamma="scale",
                class_weight="balanced",
                cache_size=1000,
                random_state=SEED,
            ),
            True,
        ),
        "mlp": (
            lambda: MLPClassifier(
                hidden_layer_sizes=(128,),
                max_iter=300,
                early_stopping=True,
                random_state=SEED,
            ),
            True,
        ),
    }


def _lgbm():
    """Istanzia il classificatore gradient boosting LightGBM con bilanciamento delle classi."""
    import lightgbm as lgb

    return lgb.LGBMClassifier(
        n_estimators=300,
        learning_rate=0.1,
        num_leaves=31,
        class_weight="balanced",
        n_jobs=-1,
        verbose=-1,
        random_state=SEED,
    )


def build(model: str, dense: bool) -> Pipeline:
    """Costruisce la Pipeline Scikit-learn adattando lo stadio di pre-processing al modello."""
    ctor, _ = _models()[model]
    steps = [("logcpm", LogCPM()), ("hvg", HVGSelector(n_top=N_HVG))]
    if dense:
        # TruncatedSVD(50) + StandardScaler denso per knn, svm_rbf e mlp (vedi docstring)
        steps += [
            ("svd", TruncatedSVD(n_components=N_SVD, random_state=SEED)),
            ("scale", StandardScaler()),
        ]
    elif model != "naive_bayes":
        steps += [("scale", StandardScaler(with_mean=False))]
    # (naive_bayes: nessuno scaling, MultinomialNB lavora sui valori log-normalizzati non negativi)
    steps.append(("clf", ctor()))
    return Pipeline(steps)


def evaluate(y_true, y_pred, classes) -> dict:  # noqa: ARG001
    """Balanced accuracy, macro-F1 (classi presenti nel test) e MCC."""
    return {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": macro_f1(y_true, y_pred),
        "mcc": matthews_corrcoef(y_true, y_pred),
    }


def fit_counting_convergence(pipe: Pipeline, X, y) -> int:
    """Addestra la pipeline e restituisce il numero di ConvergenceWarning emessi.

    Gli altri avvisi vengono riemessi normalmente.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        pipe.fit(X, y)
    n_conv = 0
    for w in caught:
        if issubclass(w.category, ConvergenceWarning):
            n_conv += 1
        else:
            warnings.showwarning(w.message, w.category, w.filename, w.lineno)
    return n_conv


def run_task(adata, sp_df, task, models, max_train, scheme_override) -> pd.DataFrame:
    """Addestra e valuta la batteria di modelli su tutti i fold del task specificato."""
    sub = sp_df[sp_df.task == task]
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
    if (idx < 0).any():
        raise SystemExit(
            f"[ERRORE] {task}: {int((idx < 0).sum())} cell_id dello split assenti dall'h5ad "
            "(split e h5ad devono venire dalla stessa esecuzione di 02/04)"
        )
    X, y = adata.X[idx], sub.y.values
    classes = np.array(sorted(np.unique(y)))
    scheme = scheme_override or PRIMARY.get(task, DEFAULT_PRIMARY)
    folds = sub[scheme].values

    print(
        f"\n{'=' * 78}\n{task}: {X.shape[0]:,} cellule | {len(classes)} classi | "
        f"split: {scheme}\n{'=' * 78}"
    )

    rows = []
    for model in models:
        _, dense = _models()[model]
        t0 = time.time()
        model_rows = []
        try:
            for f in np.unique(folds):
                tr, te = np.where(folds != f)[0], np.where(folds == f)[0]
                if max_train and len(tr) > max_train:
                    rng = np.random.default_rng(SEED + int(f))
                    tr = np.sort(rng.choice(tr, max_train, replace=False))
                pipe = build(model, dense)
                # Modelli densi (SVD + scaling): BLAS a un thread, perche' TruncatedSVD e
                # l'MLP danno risultati diversi con un numero diverso di thread
                with threadpool_limits(limits=1 if dense else None, user_api="blas"):
                    n_conv = fit_counting_convergence(pipe, X[tr], y[tr])
                    pred = pipe.predict(X[te])
                model_rows.append(
                    {
                        "task": task,
                        "scheme": scheme,
                        "model": model,
                        "fold": int(f),
                        "n_train": len(tr),
                        "n_test": len(te),
                        **evaluate(y[te], pred, classes),
                        "converged": (n_conv == 0) if model in ITERATIVE else np.nan,
                    }
                )
        except Exception as e:  # noqa: BLE001
            # Un modello fallito in un fold viene escluso del tutto: niente medie su fold parziali
            print(
                f"  {model:15s} FALLITO: {type(e).__name__}: {e} "
                f"(scartati {len(model_rows)} fold gia' calcolati)"
            )
            continue
        rows += model_rows
        s = pd.DataFrame(model_rows).macro_f1
        print(
            f"  {model:15s} macro-F1 = {s.mean():.4f} ± {s.std():.4f}   [{time.time() - t0:6.0f}s]"
        )
    return pd.DataFrame(rows)


def compare(res: pd.DataFrame, task: str) -> None:
    """Esegue il confronto statistico (Nadeau-Bengio e Wilcoxon) rispetto al modello con Macro-F1 massimo."""
    d = res[res.task == task]
    if d.empty:
        return
    piv = d.pivot_table(index="fold", columns="model", values="macro_f1")
    order = piv.mean().sort_values(ascending=False)
    best = order.index[0]
    n_tr = int(d.n_train.mean())
    n_te = int(d.n_test.mean())

    print(f"\n--- {task}: confronto col migliore ({best}, {order.iloc[0]:.4f}) ---")
    print(
        "p-value: Nadeau-Bengio (corretto per fold sovrapposti) e Wilcoxon appaiato\n"
    )
    rows = []
    for m in order.index:
        if m == best:
            continue
        a, b = piv[best].values, piv[m].values
        p_nb = nadeau_bengio_ttest(a, b, n_tr, n_te)
        try:
            p_w = wilcoxon(a, b).pvalue
        except ValueError:
            p_w = np.nan
        rows.append(
            {
                "model": m,
                "macro_f1": order[m],
                "delta_vs_best": order[m] - order[best],
                "p_nadeau_bengio": p_nb,
                "p_wilcoxon": p_w,
            }
        )
        print(
            f"  {m:15s} {order[m]:.4f}  Δ={order[m] - order[best]:+.4f}  "
            f"p_NB={p_nb:.3f}  p_W={p_w:.3f}"
        )
    pd.DataFrame(rows).to_csv(TAB / f"model_comparison_{task}.csv", index=False)

    print(
        f"\nNOTA: con {len(piv)} fold i gradi di liberta' sono {len(piv) - 1} e la potenza e' bassa."
    )
    print("Un p non significativo non dimostra che due modelli siano equivalenti.")


def main() -> None:
    """Punto di ingresso principale per l'esecuzione del benchmark multi-modello."""
    ap = argparse.ArgumentParser(
        description="Confronto di 13 classificatori con validazione per paziente."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument(
        "--tasks",
        nargs="*",
        default=None,
        help="default task2_hier3 task2_v1; con --from-table tutti i task della tabella",
    )
    ap.add_argument("--models", nargs="*", default=list(_models()))
    ap.add_argument("--max-train", type=int, default=20000, help="0 = tutte le cellule")
    ap.add_argument(
        "--scheme", default=None, help="forza fold_group_5 / fold_random_5 / fold_lopo"
    )
    ap.add_argument("--tag", default="sc")
    ap.add_argument(
        "--from-table",
        action="store_true",
        help="non addestra: rigenera riepilogo, confronti e figure da model_zoo_{tag}.csv",
    )
    ap.add_argument(
        "--update",
        action="store_true",
        help="sostituisce nella tabella esistente solo i (task, modello) rieseguiti",
    )
    args = ap.parse_args()
    unknown = set(args.models) - set(_models())
    if unknown:
        raise SystemExit(f"[ERRORE] modelli sconosciuti: {sorted(unknown)}")

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    out = TAB / f"model_zoo_{args.tag}.csv"
    if args.from_table:
        # Nessun addestramento: ricalcola riepilogo, confronti e figure dai risultati per fold
        if not out.exists():
            raise SystemExit(f"[ERRORE] --from-table: {out} assente")
        res = pd.read_csv(out)
        if args.tasks is None:
            args.tasks = list(dict.fromkeys(res.task))
        missing = [t for t in args.tasks if t not in set(res.task)]
        if missing:
            print(f"[!!] --from-table: task assenti da {out}: {missing} (saltati)")
        args.tasks = [t for t in args.tasks if t in set(res.task)]
    else:
        if args.tasks is None:
            args.tasks = ["task2_hier3", "task2_v1"]
        print("[..] carico h5ad ...")
        adata = ad.read_h5ad(args.h5ad)
        adata.X = sp.csr_matrix(adata.X)
        sp_df = pd.read_csv(args.splits)
        missing = set(args.tasks) - set(sp_df.task.unique())
        if missing:
            raise SystemExit(
                f"[ERRORE] task assenti da {args.splits}: {sorted(missing)}"
            )
        res = pd.concat(
            [
                run_task(adata, sp_df, t, args.models, args.max_train, args.scheme)
                for t in args.tasks
            ],
            ignore_index=True,
        )
    if args.update and out.exists() and not args.from_table:
        old = pd.read_csv(out)
        drop = old.task.isin(args.tasks) & old.model.isin(args.models)
        res = pd.concat([old[~drop], res], ignore_index=True)
        print(f"[ok] --update: sostituite {int(drop.sum())} righe di {out}")
    if res.empty:
        raise SystemExit("[ERRORE] nessun modello completato: nessuna tabella scritta.")
    mixed = res.groupby("task").scheme.nunique()
    if (mixed > 1).any():
        raise SystemExit(
            f"[ERRORE] task con piu' schemi di split nella stessa tabella: "
            f"{sorted(mixed[mixed > 1].index)}; usa un --tag diverso per ogni schema"
        )
    if "converged" not in res.columns:
        res["converged"] = np.nan
    res.to_csv(out, index=False)
    print(f"\n[ok] scritto {out}")
    summ = (
        res.groupby(["task", "scheme", "model"])
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            mcc_mean=("mcc", "mean"),
            n_folds=("fold", "nunique"),
            n_train_mean=("n_train", "mean"),
            # fold con ConvergenceWarning (NaN per i modelli senza solutore iterativo)
            n_folds_not_converged=(
                "converged",
                lambda c: (c == False).sum() if c.notna().any() else np.nan,  # noqa: E712
            ),
        )
        .reset_index()
    )
    # differenza dalla logistica L2 dello stesso task (riferimento lineare)
    ref = summ.loc[summ.model == "logreg_l2", ["task", "scheme", "macro_f1_mean"]]
    summ = summ.merge(
        ref.rename(columns={"macro_f1_mean": "_ref"}), on=["task", "scheme"], how="left"
    )
    summ["delta_vs_logreg_l2"] = summ.macro_f1_mean - summ.pop("_ref")
    summ.to_csv(TAB / f"model_zoo_summary_{args.tag}.csv", index=False)

    for t in args.tasks:
        compare(res, t)

    for t in args.tasks:
        d = res[res.task == t]
        if d.empty:
            continue
        g = d.groupby("model").macro_f1.agg(["mean", "std"]).sort_values("mean")
        fig, ax = plt.subplots(figsize=(7, 0.42 * len(g) + 1.5))
        ax.barh(
            g.index, g["mean"], xerr=g["std"], capsize=3, color="#1b9e77", alpha=0.85
        )
        ax.set_xlabel("macro-F1 (split per paziente)")
        ax.set_title(f"{t} - {int(d.n_train.mean()):,} cellule di training per fold")
        fig.tight_layout()
        fig.savefig(FIG / f"model_zoo_{t}.png", dpi=140)
        plt.close(fig)
        print(f"[ok] figura: results/figures/model_zoo_{t}.png")


if __name__ == "__main__":
    main()

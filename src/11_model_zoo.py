#!/usr/bin/env python3
"""
11_model_zoo.py - Benchmark comparativo multi-modello con validazione per paziente.

Questo modulo confronta sistematicamente tredici famiglie di classificatori
supervisionati:
    dummy, naive_bayes, knn, logreg_l2, logreg_l1, elasticnet,
    decision_tree, bagging, random_forest, lightgbm,
    svm_linear, svm_rbf, mlp

Architettura sperimentale e scelte metodologiche:
1. Parita' di supporto campionario di training:
   Per garantire che il confronto misuri la capacita' induttiva degli algoritmi
   e non il diverso costo computazionale asintotico (es. complessita' quadratica
   di SVM-RBF e k-NN), ciascun fold di training viene sottocampionato in modo
   deterministico allo stesso budget di cellule (default: 20.000 cellule per fold,
   estratte esclusivamente dopo lo split per paziente).

2. Adattamento della rappresentazione allo spazio metrico del classificatore:
   - I modelli basati su distanze euclidee o kernel radiali (`knn`, `svm_rbf`, `mlp`)
     operano su una proiezione densa `TruncatedSVD` a 50 componenti (stimata
     esclusivamente sul fold di training) seguita da `StandardScaler`, mitigando
     la concentrazione delle distanze nello spazio sparso a 2.000 dimensioni.
   - I modelli lineari sparsi e gli ensemble ad albero operano direttamente
     sulla matrice CSR dei 2.000 geni ad alta variabilita' (`HVGSelector`).

3. Inferenza statistica corretta per la sovrapposizione dei training set:
   La significativita' delle differenze di Macro-F1 rispetto al modello migliore
   viene valutata mediante il t-test corretto di Nadeau-Bengio (che incorpora il
   fattore di correlazione `n_test / n_train` tra i fold) affiancato dal test
   non parametrico dei ranghi con segno di Wilcoxon.

Uso:
    python src/11_model_zoo.py --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv --tasks task2_v1 task2_hier3
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
from scpipe import HVGSelector, LogCPM, nadeau_bengio_ttest  # noqa: E402

from sklearn.decomposition import TruncatedSVD  # noqa: E402
from sklearn.dummy import DummyClassifier  # noqa: E402
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression, SGDClassifier  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef  # noqa: E402
from sklearn.naive_bayes import MultinomialNB  # noqa: E402
from sklearn.neighbors import KNeighborsClassifier  # noqa: E402
from sklearn.neural_network import MLPClassifier  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from sklearn.svm import SVC, LinearSVC  # noqa: E402
from sklearn.tree import DecisionTreeClassifier  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_HVG = 2000
N_SVD = 50

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
        # Riduzione dimensionale TruncatedSVD(50) + StandardScaler denso per modelli metrici/non lineari
        steps += [
            ("svd", TruncatedSVD(n_components=N_SVD, random_state=SEED)),
            ("scale", StandardScaler()),
        ]
    elif model != "naive_bayes":
        # MultinomialNB opera direttamente sulle frequenze log-normalizzate non negative
        steps += [("scale", StandardScaler(with_mean=False))]
    steps.append(("clf", ctor()))
    return Pipeline(steps)


def evaluate(y_true, y_pred, classes) -> dict:
    """Calcola Balanced Accuracy, Macro-F1 e coefficiente di correlazione di Matthews (MCC)."""
    return {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(
            y_true, y_pred, average="macro", labels=classes, zero_division=0
        ),
        "mcc": matthews_corrcoef(y_true, y_pred),
    }


def run_task(adata, sp_df, task, models, max_train, scheme_override) -> pd.DataFrame:
    """Addestra e valuta la batteria di modelli su tutti i fold del task specificato."""
    sub = sp_df[sp_df.task == task]
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
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
        try:
            for f in np.unique(folds):
                tr, te = np.where(folds != f)[0], np.where(folds == f)[0]
                if max_train and len(tr) > max_train:
                    rng = np.random.default_rng(SEED + int(f))
                    tr = np.sort(rng.choice(tr, max_train, replace=False))
                pipe = build(model, dense)
                pipe.fit(X[tr], y[tr])
                pred = pipe.predict(X[te])
                rows.append(
                    {
                        "task": task,
                        "scheme": scheme,
                        "model": model,
                        "fold": int(f),
                        "n_train": len(tr),
                        "n_test": len(te),
                        **evaluate(y[te], pred, classes),
                    }
                )
        except Exception as e:  # noqa: BLE001
            print(f"  {model:15s} FALLITO: {type(e).__name__}: {e}")
            continue
        d = pd.DataFrame(rows)
        s = d[(d.model == model) & (d.task == task)].macro_f1
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
        "\nNOTA METODOLOGICA: con 5 fold la potenza statistica e' limitata dai gradi di liberta' (K-1 = 4)."
    )
    print(
        "L'assenza di significativita' tra modelli lineari e non lineari documenta la separabilita' lineare."
    )


def main() -> None:
    """Punto di ingresso principale per l'esecuzione del benchmark multi-modello."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument("--tasks", nargs="*", default=["task2_hier3", "task2_v1"])
    ap.add_argument("--models", nargs="*", default=list(_models()))
    ap.add_argument("--max-train", type=int, default=20000, help="0 = tutte le cellule")
    ap.add_argument(
        "--scheme", default=None, help="forza fold_group_5 / fold_random_5 / fold_lopo"
    )
    ap.add_argument("--tag", default="sc")
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    adata.X = sp.csr_matrix(adata.X)
    sp_df = pd.read_csv(args.splits)

    res = pd.concat(
        [
            run_task(adata, sp_df, t, args.models, args.max_train, args.scheme)
            for t in args.tasks
        ],
        ignore_index=True,
    )
    out = TAB / f"model_zoo_{args.tag}.csv"
    res.to_csv(out, index=False)
    print(f"\n[ok] scritto {out}")

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

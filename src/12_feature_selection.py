#!/usr/bin/env python3
"""
12_feature_selection.py - Confronto di sei criteri di selezione dei geni: prestazione e stabilita' (Jaccard).

Scopo
-----
Valuta sei criteri di selezione dei geni a tre cardinalita' (`k in {50, 200, 1000}`).
Normalizzazione, prefiltro e selezione sono stimati solo sul training di ciascun fold.

Metodi:
  - `random`    : k geni estratti a caso (controllo negativo).
  - `hvg`       : rapporto varianza / media (non supervisionato).
  - `anova`     : statistica F univariata (`f_classif`).
  - `chi2`      : chi-quadrato su profili log-normalizzati non negativi.
  - `l1_logreg` : massimo |coefficiente| di una logistica L1 (solver SAGA, C = 0.05,
                  max_iter = 400). Se i coefficienti non nulli sono meno di k, i geni
                  restanti hanno punteggio zero e l'ordine fra loro e' arbitrario:
                  lo script lo dichiara (`n_nonzero`, `k_gt_nonzero`).
  - `rf_import` : importanza di Gini media di una Random Forest bilanciata (200 alberi).

Parita' del confronto: in ogni fold il training viene sottocampionato una sola volta
(campionamento casuale uniforme, senza stratificazione, `MAX_TRAIN` = 20.000 cellule,
seme SEED + fold) e TUTTI i selettori e il classificatore a valle usano esattamente
queste stesse cellule. (Nella versione precedente `l1_logreg` e `rf_import` usavano un
sotto-campione di 8.000 cellule e gli altri 20.000.)

Valutazione:
  1. Prestazione: macro-F1 e balanced accuracy di una `LogisticRegression` L2
     bilanciata fissa, addestrata sui k geni scelti.
  2. Stabilita': Jaccard medio a coppie fra gli insiemi di geni dei fold,
     J(A, B) = |A ∩ B| / |A ∪ B|. Con `--scheme fold_group_5` (default) le coppie sono
     fra i 5 fold per paziente, non fra fold LOPO.

Input
-----
  --h5ad    data/interim/sc_raw.h5ad (conteggi grezzi)
  --splits  data/interim/splits_sc.csv (da 04_splits.py)

Output
------
  results/tables/feature_selection_{task}.csv       macro-F1 e BA per metodo x k x fold
  results/tables/feature_selection_summary_{task}.csv  per metodo x k: macro-F1 medio e sd,
                                                    balanced accuracy media, n fold e
                                                    any_k_gt_nonzero (k > geni con punteggio
                                                    > 0 in almeno un fold: F1 e Jaccard a quel
                                                    k dipendono da pareggi arbitrari)
  results/tables/feature_stability_{task}.csv       Jaccard medio fra fold per metodo x k
  results/tables/feature_selection_diag_{task}.csv  per metodo x fold: cellule usate, geni con
                                                    punteggio > 0, convergenza di SAGA, iterazioni
  results/tables/feature_core_genes_{task}.csv      geni scelti in tutti i fold (k = 50) per metodo
  results/figures/feature_selection_{task}.png
  (i file dipendono solo da --task: un'esecuzione con un altro --scheme o con un
  sottoinsieme di --methods li sovrascrive)

Argomenti
---------
  --h5ad, --splits (obbligatori), --task (default task2_v1),
  --scheme (default fold_group_5), --methods (default tutti)

Uso
---
    python src/12_feature_selection.py --h5ad data/interim/sc_raw.h5ad \\
        --splits data/interim/splits_sc.csv --task task2_v1
    python src/12_feature_selection.py --h5ad data/interim/sc_raw.h5ad \\
        --splits data/interim/splits_sc.csv --task task1a_paired --scheme fold_lopo
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
from scpipe import LogCPM, macro_f1, mean_pairwise_jaccard  # noqa: E402

from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.exceptions import ConvergenceWarning  # noqa: E402
from sklearn.feature_selection import chi2, f_classif  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
K_VALUES = (50, 200, 1000)
# Soglia minima di presenza (numero di cellule di training con espressione > 0)
MIN_CELLS_GENE = 20
# Cellule di training per fold, condivise da tutti i selettori e dal classificatore
MAX_TRAIN = 20000
L1_MAX_ITER = 400


def sel_hvg(X, y, rng):  # noqa: ARG001
    """Rapporto di dispersione (varianza / media) sulle cellule di training."""
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = np.maximum(sq - mean**2, 0.0)
    return np.where(mean > 1e-3, var / (mean + 1e-12), -np.inf), {}


def sel_anova(X, y, rng):  # noqa: ARG001
    """Statistica F univariata fra le classi."""
    with warnings.catch_warnings():
        # geni costanti nel training: F non definita (divisione per zero), punteggio -inf
        warnings.simplefilter("ignore", RuntimeWarning)
        f, _ = f_classif(X, y)
    return np.nan_to_num(f, nan=-np.inf), {}


def sel_chi2(X, y, rng):  # noqa: ARG001
    """Statistica chi-quadrato univariata."""
    c, _ = chi2(X, y)
    return np.nan_to_num(c, nan=-np.inf), {}


def sel_l1(X, y, rng):  # noqa: ARG001
    """Massimo |coefficiente| di una logistica L1 (SAGA); registra la convergenza."""
    m = LogisticRegression(
        l1_ratio=1.0,  # penalita' L1 pura (equivale a penalty="l1", deprecato in sklearn 1.8)
        C=0.05,
        solver="saga",
        max_iter=L1_MAX_ITER,
        class_weight="balanced",
        random_state=SEED,
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", ConvergenceWarning)
        m.fit(X, y)
    # gli avvisi diversi dalla mancata convergenza vengono riemessi, non scartati
    for x in w:
        if not issubclass(x.category, ConvergenceWarning):
            warnings.showwarning(x.message, x.category, x.filename, x.lineno)
    n_iter = int(np.max(m.n_iter_))
    converged = not any(issubclass(x.category, ConvergenceWarning) for x in w)
    return np.abs(m.coef_).max(axis=0), {"n_iter": n_iter, "converged": converged}


def sel_rf(X, y, rng):  # noqa: ARG001
    """Importanza di Gini (Mean Decrease in Impurity) di una Random Forest bilanciata."""
    m = RandomForestClassifier(
        n_estimators=200,
        n_jobs=-1,
        class_weight="balanced_subsample",
        random_state=SEED,
    ).fit(X, y)
    return m.feature_importances_, {}


def sel_random(X, y, rng):  # noqa: ARG001
    """Controllo negativo: punteggi casuali uniformi."""
    return rng.random(X.shape[1]), {}


SELECTORS = {
    "random": sel_random,
    "hvg": sel_hvg,
    "anova": sel_anova,
    "chi2": sel_chi2,
    "l1_logreg": sel_l1,
    "rf_import": sel_rf,
}


def run_fold(
    Xtr_raw, ytr, Xte_raw, yte, genes, classes, methods, rng, k_values=K_VALUES
) -> dict:
    """Esegue tutti i selettori e il classificatore a valle su un fold.

    `classes` non e' usato (resta nella firma per compatibilita' con i test).
    `Xtr_raw` / `Xte_raw` sono conteggi grezzi delle cellule di training (gia'
    sottocampionate) e di test. Normalizzazione, prefiltro dei geni e selezione
    usano solo il training; il test serve solo per la predizione finale.
    Restituisce un dizionario con `rows` (metriche per metodo x k), `chosen`
    ((metodo, k) -> insieme di geni) e `diag` (diagnostica per metodo).
    """
    lc = LogCPM().fit(Xtr_raw)
    Xtr, Xte = lc.transform(Xtr_raw), lc.transform(Xte_raw)
    keep = np.asarray((Xtr > 0).sum(axis=0)).ravel() >= MIN_CELLS_GENE
    Xtr, Xte = Xtr[:, keep], Xte[:, keep]
    gsub = genes[keep]
    rows, chosen, diag = [], {}, []
    for meth in methods:
        t0 = time.time()
        score, info = SELECTORS[meth](Xtr, ytr, rng)
        n_pos = int(np.sum(np.isfinite(score) & (score > 0)))
        rank = np.argsort(-score)
        diag.append(
            {
                "method": meth,
                "n_train": Xtr.shape[0],
                "n_genes_prefilter": int(keep.sum()),
                "n_nonzero": n_pos,
                "converged": info.get("converged", np.nan),
                "n_iter": info.get("n_iter", np.nan),
            }
        )
        for k in k_values:
            sel = np.sort(rank[:k])
            chosen[(meth, k)] = set(gsub[sel])

            # Standardizzazione senza centramento e classificatore lineare fisso
            sc = StandardScaler(with_mean=False).fit(Xtr[:, sel])
            clf = LogisticRegression(
                C=1.0,
                solver="lbfgs",
                max_iter=2000,
                class_weight="balanced",
                random_state=SEED,
            )
            clf.fit(sc.transform(Xtr[:, sel]), ytr)
            pred = clf.predict(sc.transform(Xte[:, sel]))
            rows.append(
                {
                    "method": meth,
                    "k": k,
                    "macro_f1": macro_f1(yte, pred),
                    "balanced_accuracy": balanced_accuracy_score(yte, pred),
                    "n_nonzero": n_pos,
                    "k_gt_nonzero": k > n_pos,
                }
            )
        extra = ""
        if meth == "l1_logreg":
            extra = f" | iter {info['n_iter']} | convergenza {info['converged']}"
        warn = (
            f"  [!] solo {n_pos} geni con punteggio > 0"
            if n_pos < max(k_values)
            else ""
        )
        print(
            f"  {meth:11s} [{time.time() - t0:5.0f}s] non nulli {n_pos:,}{extra}{warn}"
        )
    return {"rows": rows, "chosen": chosen, "diag": diag}


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Mostra il default in --help solo quando esiste (non per argomenti obbligatori o None)."""

    def _get_help_string(self, action):
        if action.default is None or action.required:
            return action.help
        return super()._get_help_string(action)


def main() -> None:
    """Esegue la cross-validation dei metodi di selezione e calcola la stabilita' di Jaccard."""
    ap = argparse.ArgumentParser(
        formatter_class=_HelpFormatter,
        description="Confronto di sei criteri di selezione dei geni (prestazione e Jaccard).",
    )
    ap.add_argument("--h5ad", type=Path, required=True, help="AnnData scRNA-seq")
    ap.add_argument(
        "--splits", type=Path, required=True, help="split prodotti da 04_splits.py"
    )
    ap.add_argument("--task", default="task2_v1", help="task da valutare")
    ap.add_argument(
        "--scheme",
        default="fold_group_5",
        choices=["fold_group_5", "fold_lopo", "fold_random_5"],
        help="colonna dei fold (decide anche i fold del Jaccard)",
    )
    ap.add_argument(
        "--methods",
        nargs="+",
        default=list(SELECTORS),
        choices=list(SELECTORS),
        help="criteri di selezione dei geni",
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    adata.X = sp.csr_matrix(adata.X)
    genes = adata.var_names.values

    sp_df = pd.read_csv(args.splits)
    sub = sp_df[sp_df.task == args.task]
    if sub.empty:
        raise SystemExit(f"[ERRORE] task {args.task} assente da {args.splits}")
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
    if (idx < 0).any():
        raise SystemExit(
            f"[ERRORE] {args.task}: {int((idx < 0).sum())} cell_id dello split assenti "
            "dall'h5ad (split e h5ad devono venire dalla stessa esecuzione di 02/04)"
        )
    X, y = adata.X[idx], sub.y.values
    folds = sub[args.scheme].values
    classes = np.array(sorted(np.unique(y)))
    print(
        f"[ok] {args.task}: {X.shape[0]:,} cellule | {len(classes)} classi | {args.scheme}"
    )

    rows, diag = [], []
    chosen = {(m, k): [] for m in args.methods for k in K_VALUES}
    for f in np.unique(folds):
        rng = np.random.default_rng(SEED + int(f))
        tr, te = np.where(folds != f)[0], np.where(folds == f)[0]
        if len(np.unique(y[tr])) < 2:
            raise SystemExit(
                f"[ERRORE] {args.task}, fold {f}: il training ha una sola classe"
            )
        if len(tr) > MAX_TRAIN:
            # Sottocampione casuale uniforme del solo training, condiviso da tutti i selettori
            tr = np.sort(rng.choice(tr, MAX_TRAIN, replace=False))
        print(f"\nfold {f}: train {len(tr):,} | test {len(te):,}")
        out = run_fold(X[tr], y[tr], X[te], y[te], genes, classes, args.methods, rng)
        rows += [{**r, "fold": int(f)} for r in out["rows"]]
        diag += [{**d, "fold": int(f)} for d in out["diag"]]
        for key, s in out["chosen"].items():
            chosen[key].append(s)

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"feature_selection_{args.task}.csv", index=False)
    pd.DataFrame(diag).to_csv(
        TAB / f"feature_selection_diag_{args.task}.csv", index=False
    )

    print("\n" + "=" * 78)
    print("(a) PRESTAZIONE: macro-F1 (media ± sd sui fold)")
    print("=" * 78)
    perf = res.pivot_table(
        index="method", columns="k", values="macro_f1", aggfunc=["mean", "std"]
    )
    print(perf.round(4).to_string())
    res.groupby(["method", "k"]).agg(
        macro_f1_mean=("macro_f1", "mean"),
        macro_f1_std=("macro_f1", "std"),
        balanced_accuracy_mean=("balanced_accuracy", "mean"),
        n_folds=("fold", "nunique"),
        any_k_gt_nonzero=("k_gt_nonzero", "any"),
    ).reset_index().to_csv(
        TAB / f"feature_selection_summary_{args.task}.csv", index=False
    )
    short = res[res.k_gt_nonzero]
    if len(short):
        print("\n[!] (metodo, k, fold) con meno di k geni a punteggio > 0:")
        print(short[["method", "k", "fold", "n_nonzero"]].to_string(index=False))

    print("\n" + "=" * 78)
    print(f"(b) STABILITA': Jaccard medio a coppie fra i fold di {args.scheme}")
    print("=" * 78)
    stab = pd.DataFrame(
        {
            k: {m: mean_pairwise_jaccard(chosen[(m, k)]) for m in args.methods}
            for k in K_VALUES
        }
    )
    print(stab.round(3).to_string())
    stab.to_csv(TAB / f"feature_stability_{args.task}.csv")
    print("\nJaccard 1.0 = stessi geni in tutti i fold; 0.0 = insiemi disgiunti.")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for meth in args.methods:
        d = res[res.method == meth].groupby("k").macro_f1.agg(["mean", "std"])
        axes[0].errorbar(
            d.index, d["mean"], yerr=d["std"], marker="o", capsize=3, label=meth
        )
        axes[1].plot(stab.columns, stab.loc[meth], marker="s", label=meth)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("numero di geni (k)")
    axes[0].set_ylabel("macro-F1")
    axes[0].set_title("prestazione")
    axes[0].legend(fontsize=8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("numero di geni (k)")
    axes[1].set_ylabel("Jaccard medio fra fold")
    axes[1].set_ylim(0, 1)
    axes[1].set_title("stabilita' della selezione")
    axes[1].legend(fontsize=8)
    fig.suptitle(f"{args.task} - selezione stimata solo sul training di ogni fold")
    fig.tight_layout()
    fig.savefig(FIG / f"feature_selection_{args.task}.png", dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: results/figures/feature_selection_{args.task}.png")

    print("\n>>> geni scelti in tutti i fold (k=50), per metodo:")
    core_rows = []
    for meth in args.methods:
        core = set.intersection(*chosen[(meth, 50)]) if chosen[(meth, 50)] else set()
        sym = (
            adata.var.gene_symbol.reindex(sorted(core)).dropna().astype(str).tolist()
            if "gene_symbol" in adata.var
            else sorted(core)
        )
        core_rows += [{"method": meth, "k": 50, "gene": g} for g in sym]
        print(f"  {meth:11s} {len(core):2d}/50 stabili: {', '.join(sym[:12])}")
    pd.DataFrame(core_rows, columns=["method", "k", "gene"]).to_csv(
        TAB / f"feature_core_genes_{args.task}.csv", index=False
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
12_feature_selection.py - Benchmark comparativo di selezione delle feature e analisi di stabilita' di Jaccard.

Questo modulo valuta sei criteri di selezione dei geni a tre livelli di
cardinalita' (`k in {50, 200, 1000}`), applicando la selezione strettamente
all'interno del fold di training per prevenire il data leakage da feature
selection globale:

Metodi confrontati:
  - `random`   : estrazione casuale uniforme di `k` geni (controllo negativo per
                 quantificare la ridondanza informativa del trascrittoma).
  - `hvg`      : indice di dispersione non supervisionato (varianza / media).
  - `anova`    : statistica F di Fisher univariata (`f_classif`).
  - `chi2`     : statistica chi-quadrato su profili log-normalizzati non negativi.
  - `l1_logreg`: selezione embedded tramite norma massima dei coefficienti di
                 una regressione logistica con penalizzazione L1 (Lasso, solver SAGA).
  - `rf_import`: selezione embedded tramite riduzione media di impurita' di Gini
                 in una Random Forest bilanciata.

La valutazione avviene congiuntamente lungo due assi:
  1. Capacita' predittiva (Macro-F1 e Balanced Accuracy) a parita' di classificatore
     a valle (`LogisticRegression` L2 bilanciata).
  2. Stabilita' strutturale della firma genica attraverso i fold di pazienti
     distinti, misurata mediante l'indice di Jaccard medio a coppie:
         J(A, B) = |A intersezione B| / |A unione B|

Uso:
    python src/12_feature_selection.py --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv --task task2_v1
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
from scpipe import LogCPM, mean_pairwise_jaccard  # noqa: E402

from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.feature_selection import chi2, f_classif  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, f1_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

warnings.filterwarnings("ignore")

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
K_VALUES = (50, 200, 1000)
# Soglia minima di presenza (numero di cellule con espressione > 0 nel training set)
MIN_CELLS_GENE = 20
MAX_TRAIN = 20000
# Sottocampionamento dedicato per i selettori embedded ad alto costo computazionale (L1 SAGA e RF)
MAX_SELECT = 8000


def sel_hvg(X, y, rng):  # noqa: ARG001
    """Calcola il rapporto di dispersione (varianza / media) sui soli campioni di training."""
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = np.maximum(sq - mean**2, 0.0)
    return np.where(mean > 1e-3, var / (mean + 1e-12), -np.inf)


def sel_anova(X, y, rng):  # noqa: ARG001
    """Calcola la statistica F univariata tra le classi sul training set."""
    f, _ = f_classif(X, y)
    return np.nan_to_num(f, nan=-np.inf)


def sel_chi2(X, y, rng):  # noqa: ARG001
    """Calcola la statistica chi-quadrato univariata sul training set."""
    c, _ = chi2(X, y)
    return np.nan_to_num(c, nan=-np.inf)


def _sub(X, y, rng, n=MAX_SELECT):
    """Estrae un sottocampione casuale riproducibile del training set per i selettori embedded."""
    if X.shape[0] <= n:
        return X, y
    i = np.sort(rng.choice(X.shape[0], n, replace=False))
    return X[i], y[i]


def sel_l1(X, y, rng):
    """Stima l'importanza genica come massimo valore assoluto dei pesi di una regressione logistica L1."""
    Xs, ys = _sub(X, y, rng)
    m = LogisticRegression(
        penalty="l1",
        C=0.05,
        solver="saga",
        max_iter=400,
        class_weight="balanced",
        random_state=SEED,
    ).fit(Xs, ys)
    return np.abs(m.coef_).max(axis=0)


def sel_rf(X, y, rng):
    """Stima l'importanza genica tramite Mean Decrease in Impurity di una Random Forest bilanciata."""
    Xs, ys = _sub(X, y, rng)
    m = RandomForestClassifier(
        n_estimators=200,
        n_jobs=-1,
        class_weight="balanced_subsample",
        random_state=SEED,
    ).fit(Xs, ys)
    return m.feature_importances_


def sel_random(X, y, rng):  # noqa: ARG001
    """Controllo negativo: assegna punteggi casuali uniformi per selezionare k geni non informativi."""
    return rng.random(X.shape[1])


SELECTORS = {
    "random": sel_random,
    "hvg": sel_hvg,
    "anova": sel_anova,
    "chi2": sel_chi2,
    "l1_logreg": sel_l1,
    "rf_import": sel_rf,
}


def main() -> None:
    """Esegue la cross-validation dei metodi di feature selection e calcola la stabilita' di Jaccard."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument("--task", default="task2_v1")
    ap.add_argument("--scheme", default="fold_group_5")
    ap.add_argument("--methods", nargs="*", default=list(SELECTORS))
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    adata.X = sp.csr_matrix(adata.X)
    genes = adata.var_names.values

    sp_df = pd.read_csv(args.splits)
    sub = sp_df[sp_df.task == args.task]
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
    X, y = adata.X[idx], sub.y.values
    folds = sub[args.scheme].values
    classes = np.array(sorted(np.unique(y)))
    print(
        f"[ok] {args.task}: {X.shape[0]:,} cellule | {len(classes)} classi | {args.scheme}"
    )

    rng = np.random.default_rng(SEED)
    rows, chosen = [], {(m, k): [] for m in args.methods for k in K_VALUES}

    for f in np.unique(folds):
        tr, te = np.where(folds != f)[0], np.where(folds == f)[0]
        if len(tr) > MAX_TRAIN:
            tr = np.sort(rng.choice(tr, MAX_TRAIN, replace=False))

        # Normalizzazione LogCPM e prefiltraggio di espressione minima stimati esclusivamente sul training set
        lc = LogCPM().fit(X[tr])
        Xtr, Xte = lc.transform(X[tr]), lc.transform(X[te])
        keep = np.asarray((Xtr > 0).sum(axis=0)).ravel() >= MIN_CELLS_GENE
        Xtr, Xte = Xtr[:, keep], Xte[:, keep]
        gsub = genes[keep]
        print(
            f"\nfold {f}: train {len(tr):,} | test {len(te):,} | geni dopo prefiltro {keep.sum():,}"
        )

        for meth in args.methods:
            t0 = time.time()
            score = SELECTORS[meth](Xtr, y[tr], rng)
            rank = np.argsort(-score)
            for k in K_VALUES:
                sel = np.sort(rank[:k])
                chosen[(meth, k)].append(set(gsub[sel]))

                # Standardizzazione senza centramento e addestramento del classificatore lineare fisso
                sc = StandardScaler(with_mean=False).fit(Xtr[:, sel])
                clf = LogisticRegression(
                    C=1.0,
                    solver="lbfgs",
                    max_iter=2000,
                    class_weight="balanced",
                    random_state=SEED,
                )
                clf.fit(sc.transform(Xtr[:, sel]), y[tr])
                pred = clf.predict(sc.transform(Xte[:, sel]))
                rows.append(
                    {
                        "method": meth,
                        "k": k,
                        "fold": int(f),
                        "macro_f1": f1_score(
                            y[te],
                            pred,
                            average="macro",
                            labels=classes,
                            zero_division=0,
                        ),
                        "balanced_accuracy": balanced_accuracy_score(y[te], pred),
                    }
                )
            print(f"  {meth:11s} [{time.time() - t0:5.0f}s]")

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"feature_selection_{args.task}.csv", index=False)

    print("\n" + "=" * 78)
    print("(a) PERFORMANCE: macro-F1 (media ± sd sui fold)")
    print("=" * 78)
    perf = res.pivot_table(
        index="method", columns="k", values="macro_f1", aggfunc=["mean", "std"]
    )
    print(perf.round(4).to_string())

    print("\n" + "=" * 78)
    print("(b) STABILITA': Jaccard medio a coppie fra i fold")
    print("=" * 78)
    stab = pd.DataFrame(
        {
            k: {m: mean_pairwise_jaccard(chosen[(m, k)]) for m in args.methods}
            for k in K_VALUES
        }
    )
    print(stab.round(3).to_string())
    stab.to_csv(TAB / f"feature_stability_{args.task}.csv")

    print(
        "\nIndice di Jaccard: 1.0 indica invarianza perfetta dei geni selezionati attraverso i fold;"
    )
    print(
        "0.0 indica disgiunzione totale. La stabilita' inter-paziente e' prerequisito per l'interpretabilita'."
    )

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
    axes[0].set_title("performance")
    axes[0].legend(fontsize=8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("numero di geni (k)")
    axes[1].set_ylabel("Jaccard medio fra fold")
    axes[1].set_ylim(0, 1)
    axes[1].set_title("stabilita' della selezione")
    axes[1].legend(fontsize=8)
    fig.suptitle(
        f"{args.task} - selezione stimata esclusivamente sul train di ogni fold"
    )
    fig.tight_layout()
    fig.savefig(FIG / f"feature_selection_{args.task}.png", dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: results/figures/feature_selection_{args.task}.png")

    print("\n>>> geni selezionati da tutti i fold (k=50), per metodo:")
    for meth in args.methods:
        core = set.intersection(*chosen[(meth, 50)]) if chosen[(meth, 50)] else set()
        sym = (
            adata.var.gene_symbol.reindex(sorted(core)).dropna().tolist()
            if "gene_symbol" in adata.var
            else sorted(core)
        )
        print(f"  {meth:11s} {len(core):2d}/50 stabili: {', '.join(sym[:12])}")
    print(
        "\nVerifica di coerenza biologica: presenza dei marcatori di lignaggio (CD3D, CD8A, MS4A1, LYZ, COL1A1, EPCAM)"
    )
    print("nell'intersezione dei fold.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
20_xai_per_patient.py: Stabilita' inter-paziente delle spiegazioni XAI in validazione Leave-One-Patient-Out.

Inquadramento Metodologico
--------------------------
Quando il numero di osservazioni per soggetto e' elevato ma il numero di soggetti
indipendenti e' ridotto, l'analisi di spiegabilita' (Explainable AI, XAI) attraverso
i fold di una cross-validation Leave-One-Patient-Out (LOPO) consente di valutare se
un classificatore lineare apprenda una regola trascrizionale invariante rispetto al
paziente o se adatti la propria frontiera decisionale alle specificita' dei soggetti
inclusi nel training set.

Se le 67.662 cellule della coorte appaiata costituissero osservazioni indipendenti,
l'esclusione di circa 13.000 cellule (un paziente su cinque) lascerebbe pressoche'
inalterato il vettore dei coefficienti appresi. La divergenza tra le spiegazioni
dei fold LOPO rispetto a un controllo nullo (in cui viene rimosso il medesimo numero
di cellule campionate casualmente ignorando l'identita' del paziente) misura pertanto
direttamente il grado di ridondanza intra-paziente (pseudoreplicazione) a livello
dello spazio delle feature.

Metriche di Stabilita' e Metodo XAI
-----------------------------------
Lo script utilizza l'interpretabilita' intrinseca di un modello di Regressione
Logistica con penalizzazione L1 (Lasso, `C = 0.1`, pesi bilanciati per classe),
preferita ad approssimazioni locali post-hoc (quali SHAP) per l'esattezza analitica
dei pesi lineari e l'efficienza computazionale su matrici sparse ad alta dimensionalita'.
Per ciascuno dei 5 fold LOPO vengono calcolati:
  1. L'indice di similarita' di Jaccard tra gli insiemi dei `TOP_K = 50` geni a
     maggiore valore assoluto del coefficiente massimo multi-classe;
  2. La correlazione di rango di Spearman tra i vettori completi dei coefficienti;
  3. L'intersezione globale dei geni condivisi da tutti i 5 fold rispetto ai geni
     esclusivi di un singolo fold (spiegazioni paziente-specifiche);
  4. Il confronto diretto con la stabilita' di Jaccard ottenuta su 5 repliche di
     controllo prive di struttura per paziente (`N_CONTROL = 5`).

Esecuzione:
    python src/20_xai_per_patient.py \
        --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv \
        --task task1a_paired
"""

from __future__ import annotations

import argparse
import sys
import warnings
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import LogCPM  # noqa: E402

from scipy.stats import spearmanr  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

warnings.filterwarnings("ignore")

# Percorsi di output e iperparametri del protocollo di stabilita' XAI
FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
MAX_TRAIN = 20000  # Tetto massimo di cellule di training per fold (allineato a 11 e 12)
MIN_CELLS_GENE = (
    10  # Soglia minima di espressione sul training set per il pre-filtraggio genico
)
TOP_K = 50  # Cardinalita' del sottoinsieme di geni principali per la metrica di Jaccard
N_CONTROL = (
    5  # Numero di repliche casuali per il controllo privo di struttura di gruppo
)


def jaccard(a: set, b: set) -> float:
    """Calcola l'indice di similarita' di Jaccard tra due insiemi finiti di feature."""
    return len(a & b) / max(len(a | b), 1)


def fit_explanation(X, y, genes, rng):
    """Addestra un classificatore di Regressione Logistica L1 sul fold fornito
    (previa normalizzazione LogCPM, filtro di espressione minima e scaling senza
    centratura densa) e restituisce il vettore completo dei coefficienti massimi
    in valore assoluto insieme all'insieme dei TOP_K geni piu' rilevanti."""
    # Normalizzazione per-cellula LogCPM e filtro dei geni espressi in almeno MIN_CELLS_GENE
    lc = LogCPM().fit(X)
    Xn = lc.transform(X)
    keep = np.asarray((Xn > 0).sum(axis=0)).ravel() >= MIN_CELLS_GENE
    Xn, gsub = Xn[:, keep], genes[keep]
    # Standardizzazione della varianza preservando la sparsita' (with_mean=False)
    Xs = StandardScaler(with_mean=False).fit_transform(Xn)

    # Addestramento del modello lineare sparso L1 con bilanciamento inverso delle frequenze di classe
    clf = LogisticRegression(
        penalty="l1",
        C=0.1,
        solver="liblinear",
        max_iter=1000,
        class_weight="balanced",
        random_state=SEED,
    )
    clf.fit(Xs, y)
    # Aggregazione multi-classe tramite massimo valore assoluto del peso per ciascun gene
    coef = np.abs(clf.coef_).max(axis=0)

    top = set(gsub[np.argsort(-coef)[:TOP_K]])
    # Ricostruzione del vettore denso sull'intero spazio genico originale per il calcolo di Spearman
    full = pd.Series(0.0, index=genes)
    full[gsub] = coef
    return full, top


def main() -> None:
    """Orchestra l'estrazione delle spiegazioni L1 sui 5 fold LOPO, il calcolo delle
    matrici di similarita' di Jaccard e Spearman e il confronto con il controllo casuale."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument("--task", default="task1a_paired")
    ap.add_argument("--scheme", default="fold_lopo")
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    # Caricamento della matrice di conteggio CSR e indicizzazione delle cellule del task
    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    adata.X = sp.csr_matrix(adata.X)
    genes = adata.var_names.values

    sp_df = pd.read_csv(args.splits)
    sub = sp_df[sp_df.task == args.task]
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
    X, y = adata.X[idx], sub.y.values
    folds = sub[args.scheme].values
    patients = sub.patient.values
    rng = np.random.default_rng(SEED)

    print(
        f"[ok] {args.task}: {X.shape[0]:,} cellule | {len(np.unique(patients))} pazienti"
    )
    print("\ncellule per paziente:")
    for p, n in pd.Series(patients).value_counts().sort_index().items():
        print(f"  {p}: {n:,}")

    # Estrazione dei profili di coefficienti L1 e dei TOP_K geni per ciascun fold LOPO
    print(
        f"\n{'=' * 74}\nSPIEGAZIONI: una per fold LOPO (= una per paziente escluso)\n{'=' * 74}"
    )
    coefs, tops, labels = {}, {}, {}
    for f in np.unique(folds):
        tr = np.where(folds != f)[0]
        held_out = pd.Series(patients[folds == f]).mode()[0]
        if len(tr) > MAX_TRAIN:
            tr = np.sort(rng.choice(tr, MAX_TRAIN, replace=False))
        full, top = fit_explanation(X[tr], y[tr], genes, rng)
        coefs[f], tops[f], labels[f] = full, top, held_out
        nz = int((full > 0).sum())
        print(
            f"  fold {f} (escluso {held_out}): train {len(tr):,} cellule | "
            f"{nz:,} geni con coefficiente non nullo"
        )

    # Calcolo delle matrici simmetriche di similarita' (Jaccard sui top-K e Spearman globale)
    print(
        f"\n{'=' * 74}\n(1) QUANTO SI ASSOMIGLIANO LE SPIEGAZIONI DEI 5 FOLD\n{'=' * 74}"
    )
    fs = sorted(coefs)
    names = [labels[f] for f in fs]
    J = pd.DataFrame(np.eye(len(fs)), index=names, columns=names)
    R = pd.DataFrame(np.eye(len(fs)), index=names, columns=names)
    for a, b in combinations(fs, 2):
        na, nb = labels[a], labels[b]
        J.loc[na, nb] = J.loc[nb, na] = jaccard(tops[a], tops[b])
        r = spearmanr(coefs[a].values, coefs[b].values).correlation
        R.loc[na, nb] = R.loc[nb, na] = r

    off = ~np.eye(len(fs), dtype=bool)
    j_mean, r_mean = J.values[off].mean(), R.values[off].mean()
    print(f"\nJaccard sui top-{TOP_K} geni (0 = spiegazioni disgiunte, 1 = identiche):")
    print(J.round(3).to_string())
    print(f"\n  media fuori diagonale: {j_mean:.3f}")
    print(f"\nSpearman fra i vettori di coefficienti: media {r_mean:.3f}")
    J.to_csv(TAB / f"xai_jaccard_{args.task}.csv")
    R.to_csv(TAB / f"xai_spearman_{args.task}.csv")

    # Identificazione del nucleo di geni invarianti condivisi da tutti i 5 fold rispetto ai geni esclusivi
    print(f"\n{'=' * 74}\n(2) GENI CONDIVISI DA TUTTI I FOLD\n{'=' * 74}")
    common = set.intersection(*tops.values())
    union = set.union(*tops.values())
    print(f"geni nell'unione delle 5 spiegazioni : {len(union)}")
    print(
        f"geni presenti in TUTTE e 5           : {len(common)} ({len(common) / TOP_K:.0%} di top-{TOP_K})"
    )
    if common:
        sym = (
            adata.var.gene_symbol.reindex(sorted(common)).dropna().tolist()
            if "gene_symbol" in adata.var
            else sorted(common)
        )
        print(f"  -> {', '.join(map(str, sym[:20]))}")
    print("\ngeni esclusivi di un solo fold (spiegazione paziente-specifica):")
    for f in fs:
        others = set.union(*[tops[g] for g in fs if g != f])
        excl = tops[f] - others
        print(
            f"  {labels[f]}: {len(excl)}/{TOP_K} geni non condivisi con nessun altro fold"
        )

    # Esperimento di controllo: ricampionamento casuale di uguale taglia ignorando la struttura per paziente
    print(
        f"\n{'=' * 74}\n(3) CONTROLLO: la variabilita' e' effetto-paziente o rumore?\n{'=' * 74}"
    )
    print("Ripeto lo stesso esperimento ignorando il paziente: sottoinsiemi casuali")
    print("della stessa dimensione, campionati da tutte le cellule mescolate.\n")
    ctrl_tops = []
    n_train = min(MAX_TRAIN, int(np.median([np.sum(folds != f) for f in fs])))
    for i in range(N_CONTROL):
        sel = rng.choice(len(y), n_train, replace=False)
        _, top = fit_explanation(X[sel], y[sel], genes, rng)
        ctrl_tops.append(top)
        print(f"  replica {i + 1}/{N_CONTROL} completata")
    j_ctrl = np.mean([jaccard(a, b) for a, b in combinations(ctrl_tops, 2)])

    print(f"\n  Jaccard medio FRA PAZIENTI (LOPO)        : {j_mean:.3f}")
    print(f"  Jaccard medio SENZA struttura (controllo): {j_ctrl:.3f}")
    print(f"  differenza                                : {j_ctrl - j_mean:+.3f}")
    print()
    # Confronto tra la stabilita' osservata in LOPO e la baseline di campionamento casuale
    if j_ctrl - j_mean > 0.05:
        print(
            "  LETTURA: le spiegazioni sono piu' simili quando si ignora il paziente."
        )
        print(
            "  La variabilita' osservata in LOPO eccede quella di puro campionamento:"
        )
        print(
            "  rimuovere UN paziente cambia cio' che il modello impara piu' di quanto"
        )
        print(
            "  lo cambi rimuovere lo stesso numero di cellule prese a caso. Le cellule"
        )
        print("  dello stesso paziente sono quindi ridondanti fra loro, e l'unita'")
        print("  informativa e' il paziente, non la cellula.")
        print(
            "\n  CAUTELA: su dati simulati la differenza fra 'effetto forte' e 'nessun"
        )
        print(
            "  effetto' era di soli 0.02 punti. Il segno e' informativo, la magnitudine"
        )
        print("  va riportata come indicativa, non come stima puntuale dell'effetto.")
    elif j_ctrl - j_mean > 0.02:
        print("  LETTURA: differenza nella direzione attesa ma di entita' contenuta.")
        print("  Compatibile sia con un effetto-paziente moderato sia con rumore di")
        print("  campionamento: da riportare come indicazione, non come dimostrazione.")
    else:
        print(
            "  LETTURA: nessuna differenza apprezzabile, o differenza di segno opposto."
        )
        print("  Le spiegazioni variano fra i fold quanto varierebbero ricampionando a")
        print(
            "  caso: la variabilita' e' rumore, non effetto-paziente. Questo risultato"
        )
        print(
            "  NON supporta l'argomento sulla pseudoreplicazione per questa via, e va"
        )
        print("  riportato come tale. L'argomento resta comunque sostenuto dal divario")
        print("  fra split casuale e split per paziente misurato al paragrafo 3.2.")

    pd.DataFrame(
        [
            {
                "jaccard_lopo": j_mean,
                "jaccard_controllo": j_ctrl,
                "differenza": j_ctrl - j_mean,
                "spearman_lopo": r_mean,
                "geni_comuni": len(common),
                "top_k": TOP_K,
            }
        ]
    ).to_csv(TAB / f"xai_summary_{args.task}.csv", index=False)

    # Generazione della figura bicolonna (heatmap di Jaccard LOPO e confronto con il controllo)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    im = axes[0].imshow(J.values, cmap="Blues", vmin=0, vmax=1)
    axes[0].set_xticks(range(len(names)), names, rotation=45, ha="right")
    axes[0].set_yticks(range(len(names)), names)
    for i in range(len(names)):
        for j in range(len(names)):
            axes[0].text(
                j,
                i,
                f"{J.values[i, j]:.2f}",
                ha="center",
                va="center",
                color="white" if J.values[i, j] > 0.5 else "black",
                fontsize=10,
            )
    axes[0].set_title(
        f"Similarita' fra spiegazioni (Jaccard, top-{TOP_K} geni)\n"
        "ogni cella: due fold LOPO a confronto",
        fontsize=11,
    )
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    axes[1].bar(
        ["fra pazienti\n(LOPO)", "senza struttura\n(controllo)"],
        [j_mean, j_ctrl],
        color=["#C1502E", "#1C7293"],
        alpha=0.85,
        width=0.55,
    )
    axes[1].set_ylabel(f"Jaccard medio fra spiegazioni (top-{TOP_K})")
    axes[1].set_ylim(0, max(j_ctrl, j_mean) * 1.35 + 0.05)
    for i, v in enumerate([j_mean, j_ctrl]):
        axes[1].text(
            i, v + 0.01, f"{v:.3f}", ha="center", fontsize=12, fontweight="bold"
        )
    axes[1].set_title(
        "Le spiegazioni divergono per effetto-paziente\no per rumore di campionamento?",
        fontsize=11,
    )
    fig.tight_layout()
    out = FIG / f"xai_per_patient_{args.task}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: {out}")
    print(f"[ok] tabelle: {TAB}/xai_*_{args.task}.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
20_xai_per_patient.py - Stabilita' fra pazienti delle spiegazioni di una logistica L1 (fold LOPO contro controllo casuale).

Scopo
-----
Se le cellule della coorte appaiata (task1a_paired: 67.662 cellule, 5 pazienti)
fossero osservazioni indipendenti, togliere dal training un paziente (circa 13.000
cellule) cambierebbe poco i coefficienti appresi. Lo script confronta quanto si
somigliano le spiegazioni dei fold leave-one-patient-out (LOPO) con quanto si
somigliano spiegazioni ottenute da sottoinsiemi casuali di cellule della stessa
dimensione dei training LOPO sottocampionati (20.000 cellule): e' uguale la dimensione
del training, non il numero di cellule tolte.

Metodo: logistica L1 (C = 0.1, liblinear, classi bilanciate) su log1p(CP10k) dei geni
espressi in almeno 10 cellule di training, scalati senza centratura; importanza di
un gene = |coefficiente|. Training di ogni fold sottocampionato a 20.000 cellule.
  1. Jaccard fra i top-50 geni di ogni coppia di fold LOPO;
  2. Spearman fra i vettori completi dei coefficienti (tutti i geni di `var`: i geni
     esclusi dal prefiltro valgono 0 in ogni fold, quindi il valore dipende anche
     dall'insieme dei geni considerati);
  3. geni presenti nei top-50 di tutti i fold (salvati con il simbolo) e geni
     esclusivi di un fold;
  4. controllo: 5 sottoinsiemi casuali di cellule della stessa dimensione, senza
     struttura per paziente; Jaccard medio fra le loro spiegazioni.
La differenza controllo - LOPO e' una misura indicativa: la simulazione
`src/sim/sim_d_xai_stability.py` (risultati in `results/tables/sim_d_xai_stability.csv`)
riporta la differenza con e senza un effetto paziente costruito, su pochi pazienti e
con una sola forma di effetto. Il controllo non e' perfettamente appaiato: due
sottoinsiemi casuali condividono in media piu' cellule di due training LOPO, quindi
la differenza e' positiva anche senza effetto paziente (in simulazione circa +0,03).
I geni comuni a tutti i fold includono geni della firma di stress da dissociazione
(FOS, BTG2, HSPH1 nella versione pubblicata), geni di risposta immediata che marcano
anche attivazione o residenza tissutale delle T (CD69, DUSP4, CREM) e CTRB1, marcatore
acinare (probabile RNA ambientale): non sono necessariamente biologia tumore/normale.

Input
-----
  --h5ad    data/interim/sc_raw.h5ad
  --splits  data/interim/splits_sc.csv (da 04_splits.py)
  results/tables/sim_d_xai_stability.csv (facoltativo: se presente, il messaggio finale
            confronta la differenza osservata con gli scenari simulati; nessuna soglia fissa)

Output
------
  results/tables/xai_jaccard_{task}.csv        matrice di Jaccard fra fold (top-50); indice = paziente
                                               escluso dal training del fold
  results/tables/xai_spearman_{task}.csv       matrice di Spearman fra fold (stesso indice)
  results/tables/xai_summary_{task}.csv        Jaccard LOPO, controllo, differenza, Spearman, n geni comuni
  results/tables/xai_common_genes_{task}.csv   geni nei top-50 di tutti i fold (id e simbolo)
  results/figures/xai_per_patient_{task}.png
  (i file dipendono solo da --task)

Argomenti
---------
  --h5ad, --splits (obbligatori), --task (default task1a_paired), --scheme (solo fold_lopo: ogni fold deve escludere un solo paziente)

Uso
---
    python src/20_xai_per_patient.py --h5ad data/interim/sc_raw.h5ad \\
        --splits data/interim/splits_sc.csv --task task1a_paired
"""

from __future__ import annotations

import argparse
import sys
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
from scpipe import LogCPM, jaccard  # noqa: E402

from scipy.stats import spearmanr  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

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
# Firma dello stress da dissociazione (geni precoci e da shock termico), geni di risposta
# immediata che marcano anche attivazione o residenza tissutale delle T, e marcatori acinari,
# per il commento sui geni comuni a tutti i fold
DISSOCIATION = {"FOS", "FOSB", "JUN", "JUNB", "EGR1", "BTG2", "DUSP1", "IER2", "NR4A1",
                "HSPA1A", "HSPA1B", "HSPH1"}  # fmt: skip
T_ACTIVATION = {"CD69", "DUSP4", "CREM"}
ACINAR = {"CTRB1", "CTRB2", "PRSS1", "CPA1", "CELA3A", "PNLIP", "CLPS"}


def fit_explanation(X, y, genes):
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
        l1_ratio=1.0,  # penalita' L1 pura (equivale a penalty="l1", deprecato in sklearn 1.8)
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


def lopo_explanations(X, y, folds, patients, genes, rng):
    """Una spiegazione per fold: logistica L1 addestrata sul solo training del fold.

    Il training di ogni fold (al massimo MAX_TRAIN cellule estratte con `rng`) e' l'unica
    parte dei dati passata a `fit_explanation`: le cellule di test non entrano.
    Restituisce (coefficienti, top-K, paziente escluso), dizionari indicizzati per fold.
    """
    coefs, tops, labels = {}, {}, {}
    for f in np.unique(folds):
        tr = np.where(folds != f)[0]
        if len(np.unique(y[tr])) < 2:
            raise SystemExit(f"[ERRORE] fold {f}: il training ha una sola classe")
        held_out = pd.Series(patients[folds == f]).mode()[0]
        if len(tr) > MAX_TRAIN:
            tr = np.sort(rng.choice(tr, MAX_TRAIN, replace=False))
        full, top = fit_explanation(X[tr], y[tr], genes)
        coefs[f], tops[f], labels[f] = full, top, held_out
        nz = int((full > 0).sum())
        print(
            f"  fold {f} (escluso {held_out}): train {len(tr):,} cellule | "
            f"{nz:,} geni con coefficiente non nullo"
        )

    return coefs, tops, labels


def control_explanations(X, y, n_train, genes, rng, n_control=N_CONTROL):
    """Controllo senza struttura per paziente: `n_control` sottoinsiemi casuali di
    `n_train` cellule (estratte da tutte, ignorando il paziente), una spiegazione ciascuno.
    Restituisce la lista degli insiemi top-K."""
    ctrl_tops = []
    for i in range(n_control):
        sel = rng.choice(len(y), n_train, replace=False)
        _, top = fit_explanation(X[sel], y[sel], genes)
        ctrl_tops.append(top)
        print(f"  replica {i + 1}/{n_control} completata")
    return ctrl_tops


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Mostra il default in --help solo quando esiste (non per argomenti obbligatori o None)."""

    def _get_help_string(self, action):
        if action.default is None or action.required:
            return action.help
        return super()._get_help_string(action)


def main() -> None:
    """Orchestra l'estrazione delle spiegazioni L1 sui fold LOPO, il calcolo delle
    matrici di similarita' di Jaccard e Spearman e il confronto con il controllo casuale."""
    ap = argparse.ArgumentParser(
        formatter_class=_HelpFormatter,
        description="Stabilita' fra pazienti delle spiegazioni L1: fold LOPO contro controllo.",
    )
    ap.add_argument("--h5ad", type=Path, required=True, help="AnnData scRNA-seq")
    ap.add_argument(
        "--splits", type=Path, required=True, help="split prodotti da 04_splits.py"
    )
    ap.add_argument("--task", default="task1a_paired", help="task da spiegare")
    ap.add_argument(
        "--scheme",
        default="fold_lopo",
        choices=["fold_lopo"],
        help="colonna dei fold (un fold per paziente)",
    )
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
    if sub.empty:
        raise SystemExit(f"[ERRORE] task {args.task} assente da {args.splits}")
    # ogni fold deve escludere un solo paziente: le spiegazioni sono etichettate per paziente
    per_fold = sub.groupby(args.scheme).patient.nunique()
    if (per_fold != 1).any():
        raise SystemExit(
            f"[ERRORE] {args.scheme}: fold con piu' di un paziente; serve uno split LOPO"
        )
    idx = adata.obs_names.get_indexer(sub.cell_id.values)
    if (idx < 0).any():
        raise SystemExit(
            f"[ERRORE] {args.task}: {int((idx < 0).sum())} cell_id dello split assenti "
            "dall'h5ad (split e h5ad devono venire dalla stessa esecuzione di 02/04)"
        )
    X, y = adata.X[idx], sub.y.values
    folds = sub[args.scheme].values
    patients = sub.patient.values
    rng = np.random.default_rng(SEED)
    if sub.empty or len(np.unique(folds)) < 2:
        raise SystemExit(
            f"[ERRORE] task {args.task}: servono almeno 2 fold in {args.scheme}"
        )

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
    coefs, tops, labels = lopo_explanations(X, y, folds, patients, genes, rng)

    # Calcolo delle matrici simmetriche di similarita' (Jaccard sui top-K e Spearman globale)
    print(
        f"\n{'=' * 74}\n(1) QUANTO SI ASSOMIGLIANO LE SPIEGAZIONI DEI {len(coefs)} FOLD\n{'=' * 74}"
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

    # Identificazione del nucleo di geni invarianti condivisi da tutti i fold rispetto ai geni esclusivi
    print(f"\n{'=' * 74}\n(2) GENI CONDIVISI DA TUTTI I FOLD\n{'=' * 74}")
    common = set.intersection(*tops.values())
    union = set.union(*tops.values())
    print(f"geni nell'unione delle {len(fs)} spiegazioni : {len(union)}")
    print(
        f"geni presenti in tutte le {len(fs)}       : {len(common)} ({len(common) / TOP_K:.0%} di top-{TOP_K})"
    )
    ids = sorted(common)
    sym = (
        adata.var.gene_symbol.reindex(ids).astype(str).tolist()
        if "gene_symbol" in adata.var
        else ids
    )
    pd.DataFrame({"gene_id": ids, "gene_symbol": sym}).sort_values(
        "gene_symbol"
    ).to_csv(TAB / f"xai_common_genes_{args.task}.csv", index=False)
    if common:
        print(f"  -> {', '.join(sorted(map(str, sym)))}")
        # firma di dissociazione, attivazione/residenza delle T, acinari (RNA ambientale)
        diss = sorted(set(sym) & DISSOCIATION)
        tact = sorted(set(sym) & T_ACTIVATION)
        acinar = sorted(set(sym) & ACINAR)
        if diss:
            print(f"  firma dello stress da dissociazione: {', '.join(diss)}")
        if tact:
            print(
                "  risposta immediata che marca anche attivazione o residenza tissutale "
                f"delle T: {', '.join(tact)}"
            )
        if acinar:
            print(
                f"  marcatori acinari (probabile RNA ambientale): {', '.join(acinar)}"
            )
    print("\ngeni presenti solo nella spiegazione addestrata senza quel paziente:")
    for f in fs:
        others = set.union(*[tops[g] for g in fs if g != f])
        excl = tops[f] - others
        print(
            f"  senza {labels[f]}: {len(excl)}/{TOP_K} geni non condivisi con nessun altro fold"
        )

    # Esperimento di controllo: ricampionamento casuale di uguale taglia ignorando la struttura per paziente
    print(
        f"\n{'=' * 74}\n(3) CONTROLLO: la variabilita' e' effetto-paziente o rumore?\n{'=' * 74}"
    )
    print("Ripeto lo stesso esperimento ignorando il paziente: sottoinsiemi casuali")
    print("della stessa dimensione, campionati da tutte le cellule mescolate.\n")
    n_train = min(MAX_TRAIN, int(np.median([np.sum(folds != f) for f in fs])))
    ctrl_tops = control_explanations(X, y, n_train, genes, rng)
    j_ctrl = np.mean([jaccard(a, b) for a, b in combinations(ctrl_tops, 2)])

    print(f"\n  Jaccard medio FRA PAZIENTI (LOPO)        : {j_mean:.3f}")
    print(f"  Jaccard medio SENZA struttura (controllo): {j_ctrl:.3f}")
    print(f"  differenza                                : {j_ctrl - j_mean:+.3f}")
    print()
    # Confronto tra la stabilita' osservata in LOPO e la baseline di campionamento casuale
    # Lettura senza soglie fisse: confronto con gli scenari simulati di sim_d (stessa procedura)
    diff = j_ctrl - j_mean
    sim_path = TAB / "sim_d_xai_stability.csv"
    if sim_path.exists():
        sim = pd.read_csv(sim_path)
        sim = sim[sim["rep"].astype(str) != "mean"]
        rng_sim = {
            c: (g.difference.min(), g.difference.max())
            for c, g in sim.groupby("condition")
        }
        for c, (lo_s, hi_s) in rng_sim.items():
            print(f"  simulazione {c}: differenza fra {lo_s:+.3f} e {hi_s:+.3f}")
        lo0, hi0 = rng_sim.get("no_patient_effect", (np.nan, np.nan))
        if np.isnan(lo0) or np.isnan(hi0):
            print(
                "  scenario senza effetto paziente assente dalla simulazione: nessun confronto."
            )
        elif diff > hi0:
            print(
                "  LETTURA: la differenza osservata e' sopra tutte le repliche simulate senza"
            )
            print(
                "  effetto paziente: togliere un paziente cambia la regola appresa piu' che"
            )
            print(
                "  addestrare su un sottoinsieme casuale della stessa dimensione, come atteso"
            )
            print("  se le cellule dello stesso paziente")
            print(
                "  sono ridondanti. La simulazione ha pochi pazienti e una sola forma di"
            )
            print("  effetto: e' un'indicazione, non una stima calibrata.")
        elif diff >= lo0:
            print(
                "  LETTURA: la differenza e' nell'intervallo delle repliche senza effetto"
            )
            print(
                "  paziente: non distingue un effetto paziente dal rumore di campionamento."
            )
        else:
            print(
                "  LETTURA: la differenza e' sotto le repliche senza effetto paziente."
            )
    else:
        print(
            f"  {sim_path} assente: esegui src/sim/sim_d_xai_stability.py per una lettura"
        )
        print(
            "  confrontata con dati simulati; nessuna conclusione da questa sola differenza."
        )
    print(
        "  Le repliche del controllo condividono fra loro piu' cellule di due training"
    )
    print("  LOPO: parte della differenza c'e' anche senza effetto paziente.")

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

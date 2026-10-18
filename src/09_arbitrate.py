#!/usr/bin/env python3
"""
09_arbitrate.py - Identita' CD4/CD8 dei cloni dai marcatori, tasso d'errore per compartimento e misure fisiche dei doppietti.

Scopo
-----
1. Arbitrato (sezione A). Il lignaggio CD4/CD8 di un clone T maturo non cambia fra
   tessuti. Per ogni (clone, compartimento) con almeno 3 cellule si calcolano la
   frazione di cellule CD8+ (CD8A > 0 o CD8B > 0) e CD4+ (CD4 > 0 senza CD8A/CD8B);
   l'identita' del clone e' CD8T o CD4T secondo la media, a peso uguale, di queste
   frazioni sui compartimenti in cui il clone compare (scartata se il margine
   |f8 - f4| < 0.20). Entrano nell'arbitrato e nel tasso d'errore solo i cloni con
   almeno 3 cellule in almeno 2 compartimenti (`n_comp >= 2`). L'errore e'
   l'etichetta di maggioranza del compartimento diversa dall'identita'. Questo schema include il compartimento giudicato nel riferimento
   (circolare): gli schemi senza circolarita' sono in 10_loco.py. IC 95% con cluster
   bootstrap sui pazienti (`scpipe.cluster_bootstrap`).
   Etichetta di maggioranza: moda, pareggi risolti dalla prima comparsa
   (`scpipe.majority_label`).

2. Misure fisiche dei doppietti (sezione B): UMI totali per goccia (`n_counts`; in un
   doppietto T+NK attesi circa n_counts(NK TCR-) + n_counts(CD8T)) e UMI del TCR per
   1000 UMI totali (diluiti dall'mRNA della seconda cellula). Le misure sono fisiche,
   ma i gruppi confrontati (NK TCR-/TCR+, CD8T TCR+, ...) vengono dalle classi
   dell'annotazione. Mediane descrittive, senza test; la stima della frazione di
   doppietti con IC sui pazienti e' in 10_loco.py.

Input
-----
  --h5ad    data/interim/sc_raw.h5ad
  --tcrdir  cartella con i 30 file VDJ

Output
------
  results/tables/arbitration.csv                 una riga per (clone, compartimento) risolto
  results/tables/error_rate_by_compartment.csv   errore per compartimento con IC; `ic_calibrato` e' falso
                                                 con <= 5 pazienti (copertura sotto il 95% in
                                                 `sim_a`), e allora si riportano anche le
                                                 ricampionature distinte (`n_distinct_resamples`)
  results/tables/doublet_physical.csv            mediane di n_counts, n_genes, UMI TCR per 1000, per classe x TCR
  results/figures/error_rate_by_compartment.png (IC con <= 5 pazienti in grigio tratteggiato),
  results/figures/doublet_physical.png
  Se la sezione A non ha dati, le sue tabelle e la figura non vengono scritte e quelle di
  un'esecuzione precedente vengono cancellate; se nessun compartimento ha almeno 3
  pazienti, error_rate_by_compartment.csv contiene la sola intestazione e la figura
  viene cancellata. Se la sezione B non ha classi con almeno 20 cellule,
  doublet_physical.csv contiene la sola intestazione e la figura viene cancellata.

Argomenti
---------
  --h5ad (obbligatorio), --tcrdir (default data/raw/tcr)

Uso
---
    python src/09_arbitrate.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import (  # noqa: E402
    SMALL_K,
    cluster_bootstrap,
    majority_label,
    small_k_diagnostics,
)

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_BOOT = 2000

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}
COMPS = ["PBMC", "Adjacent_normal", "Tumor"]
MARKERS = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1"]
FNAME_RE = re.compile(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz")

# Soglie per la risoluzione dell'identita' clonale (supporto minimo e margine tra frazioni CD8+/CD4+)
MIN_CELLS = 3
MIN_MARGIN = 0.20


def load_tcr(tcrdir: Path) -> pd.DataFrame:
    """Carica i contig VDJ e calcola il totale di UMI recettoriali (TRA + TRB) per cellula."""
    frames = []
    for p in sorted(tcrdir.glob("*_filtered_contig_annotations.csv.gz")):
        m = FNAME_RE.match(p.name)
        if not m:
            continue
        pnum, tis = m.groups()
        df = pd.read_csv(p, compression="gzip", low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df["patient"] = f"PA{int(pnum):02d}"
        df["tissue"] = TISSUE_MAP[tis]
        df["barcode16"] = df.barcode.astype(str).str.extract(
            r"([ACGT]{16})", expand=False
        )
        frames.append(
            df[["patient", "tissue", "barcode16", "chain", "cdr3_nt", "umis"]]
        )
    if not frames:
        raise SystemExit(f"[ERRORE] nessun file VDJ riconosciuto in {tcrdir}")
    c = pd.concat(frames, ignore_index=True).dropna(subset=["barcode16"])
    c = c.sort_values("umis", ascending=False).drop_duplicates(
        ["patient", "tissue", "barcode16", "chain"]
    )
    seq = c.pivot_table(
        index=["patient", "tissue", "barcode16"],
        columns="chain",
        values="cdr3_nt",
        aggfunc="first",
    )
    umi = c.pivot_table(
        index=["patient", "tissue", "barcode16"],
        columns="chain",
        values="umis",
        aggfunc="first",
    )
    for ch in ("TRA", "TRB"):
        if ch not in seq.columns:
            seq[ch], umi[ch] = np.nan, np.nan
    w = (
        seq.rename(columns={"TRA": "cdr3_TRA", "TRB": "cdr3_TRB"})
        .join(umi.rename(columns={"TRA": "umi_TRA", "TRB": "umi_TRB"}))
        .reset_index()
    )
    w["clone_b"] = w.cdr3_TRB
    w["umi_tcr"] = w[["umi_TRA", "umi_TRB"]].sum(axis=1, min_count=1)
    return w[["patient", "tissue", "barcode16", "clone_b", "umi_TRB", "umi_tcr"]]


def cluster_boot(by_pat: dict[str, np.ndarray], n: int = N_BOOT) -> tuple:
    """Media pooled e IC 95% con cluster bootstrap sui pazienti (`scpipe.cluster_bootstrap`)."""
    return cluster_bootstrap(by_pat, n, SEED)


def clone_identity(lym: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Determina l'identita' CD4T/CD8T del clone mediando con peso uniforme le frazioni compartimentali."""
    lym = lym.assign(
        cd8_pos=(lym.g_CD8A > 0) | (lym.g_CD8B > 0),
        cd4_pos=(lym.g_CD4 > 0) & (lym.g_CD8A == 0) & (lym.g_CD8B == 0),
    )
    per = (
        lym.groupby(["patient", "clone_b", "tissue"], observed=True)
        .agg(
            n=("cd8_pos", "size"),
            f8=("cd8_pos", "mean"),
            f4=("cd4_pos", "mean"),
            # moda; pareggi risolti dalla prima comparsa nell'ordine delle cellule
            lab=("all_celltype", majority_label),
        )
        .reset_index()
    )
    per = per[per.n >= MIN_CELLS]

    ident = (
        per.groupby(["patient", "clone_b"], observed=True)
        .agg(
            n_comp=("tissue", "nunique"),
            f8=("f8", "mean"),
            f4=("f4", "mean"),
            n_cells=("n", "sum"),
        )
        .reset_index()
    )
    ident["margin"] = (ident.f8 - ident.f4).abs()
    ident["identity"] = np.where(ident.f8 > ident.f4, "CD8T", "CD4T")
    ident.loc[ident.margin < MIN_MARGIN, "identity"] = "ambiguous"
    return per, ident


def partA(lym: pd.DataFrame) -> None:
    """Stima il tasso di discordanza rispetto all'identita' clonale di consenso per compartimento."""
    print("\n" + "=" * 78)
    print("(A) ARBITRATO: tasso d'errore per compartimento anatomico")
    print("=" * 78)
    per, ident = clone_identity(lym)
    conf = ident[(ident.identity != "ambiguous") & (ident.n_comp >= 2)]
    print(
        f"cloni con >={MIN_CELLS} cellule in almeno un compartimento: {len(ident):,} "
        f"(di cui in >=2 compartimenti: {(ident.n_comp >= 2).sum():,})"
    )
    print(
        f"  in >=2 compartimenti con identita' risolta (margine >= {MIN_MARGIN}): {len(conf):,} "
        f"({conf.identity.value_counts().to_dict()})"
    )
    print(f"  ambigui (scartati): {(ident.identity == 'ambiguous').sum():,}")
    if conf.empty:
        print("nessun clone risolvibile.")
        drop_stale_partA()
        return

    j = per.merge(conf[["patient", "clone_b", "identity"]], on=["patient", "clone_b"])
    j = j[j.lab.isin(["CD4T", "CD8T", "NK"])]
    j["error"] = j.lab != j.identity
    j.to_csv(TAB / "arbitration.csv", index=False)

    print("\ntasso d'errore dell'etichetta degli autori, per compartimento")
    print(
        "(riferimento = identita' clonale dai marcatori, pesata ugualmente fra compartimenti)\n"
    )
    rows = []
    for tis in COMPS:
        g = j[j.tissue == tis]
        if g.patient.nunique() < 3:
            print(f"  {tis:16s} solo {g.patient.nunique()} pazienti, salto")
            continue
        by = {p: gg.error.values.astype(float) for p, gg in g.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        k = g.patient.nunique()
        # con <= SMALL_K pazienti l'IC non e' calibrato: si riportano le ricampionature distinte
        n_dist = (
            small_k_diagnostics(by, N_BOOT, SEED)["n_distinct_resamples"]
            if k <= SMALL_K
            else np.nan
        )
        rows.append(
            {
                "tissue": tis,
                "n_clone_comp": len(g),
                "n_patients": k,
                "error_rate": m,
                "lo": lo,
                "hi": hi,
                "ic_calibrato": k > SMALL_K,
                "n_distinct_resamples": n_dist,
            }
        )
        print(
            f"  {tis:16s} errore = {m:.3f}  IC95% [{lo:.3f}, {hi:.3f}]  "
            f"(n={len(g)}, {g.patient.nunique()} pazienti)"
        )
    res = pd.DataFrame(
        rows,
        columns=[
            "tissue",
            "n_clone_comp",
            "n_patients",
            "error_rate",
            "lo",
            "hi",
            "ic_calibrato",
            "n_distinct_resamples",
        ],
    )
    res.to_csv(TAB / "error_rate_by_compartment.csv", index=False)
    if res.empty:
        print("nessun compartimento con almeno 3 pazienti: niente figura.")
        (FIG / "error_rate_by_compartment.png").unlink(missing_ok=True)
        return

    print("\ndifferenza appaiata Tumor - PBMC, sugli STESSI cloni:")
    piv = j.pivot_table(
        index=["patient", "clone_b"], columns="tissue", values="error", aggfunc="max"
    )
    if {"Tumor", "PBMC"}.issubset(piv.columns):
        d = piv.dropna(subset=["Tumor", "PBMC"]).reset_index()
        d["diff"] = d.Tumor - d.PBMC
        by = {p: g["diff"].values for p, g in d.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        k = d.patient.nunique()
        if k <= SMALL_K:
            v = f"IC non calibrato ({k} pazienti)"
        else:
            v = "IC esclude 0" if lo > 0 or hi < 0 else "IC include 0"
        print(
            f"  {m:+.3f}  IC95% [{lo:+.3f}, {hi:+.3f}]  {v}   (n={len(d)} cloni, {k} pz)"
        )

    print("\ndirezione dell'errore (etichetta autori -> identita' clonale):")
    print(
        j[j.error]
        .groupby(["tissue", "lab", "identity"])
        .size()
        .sort_values(ascending=False)
        .head(12)
        .to_string()
    )

    fig, ax = plt.subplots(figsize=(6.5, 4))
    # IC con <= SMALL_K pazienti: non calibrati, barra grigia tratteggiata
    cal = (res.n_patients > SMALL_K).values
    bars = ax.bar(
        res.tissue,
        res.error_rate,
        yerr=[res.error_rate - res.lo, res.hi - res.error_rate],
        capsize=5,
        color=["#d95f02" if c else "#bbbbbb" for c in cal],
        alpha=0.85,
    )
    for b, c in zip(bars, cal):
        if not c:
            b.set_hatch("//")
    if not cal.all():
        from matplotlib.patches import Patch

        ax.legend(
            handles=[
                Patch(
                    color="#d95f02", alpha=0.85, label=f"IC 95% (> {SMALL_K} pazienti)"
                ),
                Patch(
                    facecolor="#bbbbbb",
                    hatch="//",
                    label=f"IC non calibrato (<= {SMALL_K} pazienti)",
                ),
            ],
            fontsize=8,
        )
    ax.set_ylabel("tasso d'errore dell'etichetta")
    ax.set_title("chi sbaglia? riferimento = identita' clonale (marcatori)")
    fig.tight_layout()
    fig.savefig(FIG / "error_rate_by_compartment.png", dpi=140)
    plt.close(fig)


def drop_stale_partA() -> None:
    """Cancella le tabelle e la figura della sezione A di un'esecuzione precedente."""
    for path in (
        TAB / "arbitration.csv",
        TAB / "error_rate_by_compartment.csv",
        FIG / "error_rate_by_compartment.png",
    ):
        path.unlink(missing_ok=True)


def partB(gex: pd.DataFrame) -> None:
    """Verifica se le cellule non-T con TCR presentano massa UMI e diluizione da doppietto fisico."""
    print("\n" + "=" * 78)
    print(
        "(B) DOPPIETTO FISICO o CELLULA DI CONFINE?  (misure fisiche per classe annotata)"
    )
    print("=" * 78)
    g = gex.copy()
    # Rapporto tra UMI del recettore T e conteggio UMI totale per 1000 molecole
    g["ratio"] = g.umi_tcr / np.maximum(g.n_counts, 1) * 1000

    rows = []
    for cls in [
        "CD8T",
        "CD4T",
        "NK",
        "Platelet",
        "Endothelial",
        "Myeloid",
        "B",
        "Ductal",
    ]:
        for st, sub in (
            ("TCR+", g[(g.all_celltype == cls) & g.has_tcr]),
            ("TCR-", g[(g.all_celltype == cls) & ~g.has_tcr]),
        ):
            if len(sub) < 20:
                continue
            rows.append(
                {
                    "class": cls,
                    "tcr": st,
                    "n": len(sub),
                    "median_n_counts": sub.n_counts.median(),
                    "median_n_genes": sub.n_genes.median(),
                    "median_tcr_per_1k": sub.ratio.median() if st == "TCR+" else np.nan,
                }
            )
    out = pd.DataFrame(
        rows,
        columns=[
            "class",
            "tcr",
            "n",
            "median_n_counts",
            "median_n_genes",
            "median_tcr_per_1k",
        ],
    )
    out.to_csv(TAB / "doublet_physical.csv", index=False)
    if out.empty:
        print("nessuna classe con almeno 20 cellule: niente misure dei doppietti.")
        (FIG / "doublet_physical.png").unlink(missing_ok=True)
        return
    print(out.round(2).to_string(index=False))

    def med(cls, tcr):
        s = out[(out["class"] == cls) & (out.tcr == tcr)]
        return float(s.median_n_counts.iloc[0]) if len(s) else np.nan

    print("\n>>> PREDIZIONE quantitativa, se le NK-TCR+ sono doppietti T+NK")
    nk_neg, cd8 = med("NK", "TCR-"), med("CD8T", "TCR+")
    nk_pos = med("NK", "TCR+")
    print(
        f"  n_counts atteso  = n_counts(NK TCR-) + n_counts(CD8T) = "
        f"{nk_neg:.0f} + {cd8:.0f} = {nk_neg + cd8:.0f}"
    )
    print(f"  n_counts osservato (NK TCR+)                          = {nk_pos:.0f}")
    print(
        f"  rapporto osservato/atteso: {nk_pos / max(nk_neg + cd8, 1):.2f}   "
        "(1.0 = doppietto perfetto; ~0.5 = singoletto)"
    )

    r_cd8 = out.loc[(out["class"] == "CD8T") & (out.tcr == "TCR+"), "median_tcr_per_1k"]
    r_nk = out.loc[(out["class"] == "NK") & (out.tcr == "TCR+"), "median_tcr_per_1k"]
    if len(r_cd8) and len(r_nk):
        print(
            f"\n  UMI TCR per 1000 conte:  CD8T = {float(r_cd8.iloc[0]):.2f}   "
            f"NK-TCR+ = {float(r_nk.iloc[0]):.2f}"
        )
        # doppietto NK + CD8T: gli UMI TCR della CD8T diluiti su n_counts(NK TCR-) + n_counts(CD8T)
        print(
            f"  atteso per doppietto: {float(r_cd8.iloc[0]) * cd8 / max(nk_neg + cd8, 1):.2f}"
        )
        print(
            "  Un rapporto nelle NK piu' basso che nelle CD8T e' compatibile con la diluizione"
        )
        print(
            "  del trascritto TCR da parte dell'mRNA di una seconda cellula co-incapsulata."
        )

    print("\nCONTROLLO: le CD8T-TCR+ costituiscono il riferimento per i singoletti T.")
    print(
        "Valori di n_counts simili a quelli dei singoletti sono compatibili con singoletti;"
    )
    print("valori vicini alla somma di due cellule sono compatibili con doppietti.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    cls_order = [
        "CD8T",
        "CD4T",
        "NK",
        "Platelet",
        "Endothelial",
        "Myeloid",
        "B",
        "Ductal",
    ]
    for ax, col, lab in zip(
        axes, ["n_counts", "ratio"], ["n_counts (mediana)", "UMI TCR / 1000 conte"]
    ):
        pos = [
            g.loc[(g.all_celltype == c) & g.has_tcr, col].median() for c in cls_order
        ]
        neg = [
            g.loc[(g.all_celltype == c) & ~g.has_tcr, col].median()
            if col == "n_counts"
            else np.nan
            for c in cls_order
        ]
        x = np.arange(len(cls_order))
        if col == "n_counts":
            ax.bar(x - 0.2, neg, 0.38, label="TCR−", color="#7570b3", alpha=0.85)
        ax.bar(x + 0.2, pos, 0.38, label="TCR+", color="#d95f02", alpha=0.85)
        ax.set_xticks(x, cls_order, rotation=45, ha="right")
        ax.set_ylabel(lab)
        ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "doublet_physical.png", dpi=140)
    plt.close(fig)


def main() -> None:
    """Punto di ingresso principale per l'arbitrato clonale e il test biofisico dei doppietti."""
    ap = argparse.ArgumentParser(
        description="Identita' clonale dai marcatori, errore per compartimento, doppietti.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--h5ad",
        type=Path,
        required=True,
        help="AnnData con conteggi grezzi (da 02_build_h5ad.py)",
    )
    ap.add_argument(
        "--tcrdir",
        type=Path,
        default=Path("data/raw/tcr"),
        help="cartella con i file VDJ filtered_contig_annotations",
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype", "barcode"):
        obs[c] = obs[c].astype(str)
    obs = obs.rename(columns={"barcode": "barcode16"})
    sym = adata.var.gene_symbol.astype(str).values
    for gname in MARKERS:
        idx = np.where(sym == gname)[0]
        if len(idx):
            x = np.asarray(adata.X[:, idx[0]].todense()).ravel()
            obs[f"g_{gname}"] = np.log1p(x / np.maximum(obs.n_counts.values, 1) * 1e4)

    tcr = load_tcr(args.tcrdir)
    gex = obs.merge(
        tcr,
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    print(f"[ok] {len(gex):,} cellule | TCR: {gex.has_tcr.sum():,}")

    missing = [g for g in ("CD4", "CD8A", "CD8B") if f"g_{g}" not in obs]
    if missing:
        raise SystemExit(f"[ERRORE] marcatori assenti dall'h5ad: {missing}")
    lym = gex[gex.has_tcr & gex.all_celltype.isin(["CD4T", "CD8T", "NK"])].copy()
    if lym.empty:
        print(
            "[!!] nessuna cellula CD4T/CD8T/NK con TCR: salto l'arbitrato (sezione A)."
        )
        drop_stale_partA()
    else:
        partA(lym)
    partB(gex)


if __name__ == "__main__":
    main()

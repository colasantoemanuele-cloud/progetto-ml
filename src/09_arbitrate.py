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
   |f8 - f4| < 0.20). L'errore e' l'etichetta di maggioranza del compartimento diversa
   dall'identita'. Questo schema include il compartimento giudicato nel riferimento
   (circolare): gli schemi senza circolarita' sono in 10_loco.py. IC 95% con cluster
   bootstrap sui pazienti (`scpipe.cluster_bootstrap`).
   Etichetta di maggioranza: moda, pareggi risolti dalla prima comparsa
   (`scpipe.majority_label`).

2. Misure fisiche dei doppietti (sezione B), che non usano il clustering: UMI totali per
   goccia (`n_counts`; in un doppietto T+NK attesi circa n_counts(NK TCR-) + n_counts(CD8T))
   e UMI del TCR per 1000 UMI totali (diluiti dall'mRNA della seconda cellula).

Input
-----
  --h5ad    data/interim/sc_raw.h5ad
  --tcrdir  cartella con i 30 file VDJ

Output
------
  results/tables/arbitration.csv                 una riga per (clone, compartimento) risolto
  results/tables/error_rate_by_compartment.csv   errore per compartimento con IC
  results/tables/doublet_physical.csv            mediane di n_counts, n_genes, UMI TCR per 1000, per classe x TCR
  results/figures/error_rate_by_compartment.png, results/figures/doublet_physical.png

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
from scpipe import cluster_bootstrap, majority_label  # noqa: E402

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


def clone_identity(lym: pd.DataFrame) -> pd.DataFrame:
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
        rows.append(
            {
                "tissue": tis,
                "n_clone_comp": len(g),
                "n_patients": g.patient.nunique(),
                "error_rate": m,
                "lo": lo,
                "hi": hi,
            }
        )
        print(
            f"  {tis:16s} errore = {m:.3f}  IC95% [{lo:.3f}, {hi:.3f}]  "
            f"(n={len(g)}, {g.patient.nunique()} pazienti)"
        )
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "error_rate_by_compartment.csv", index=False)

    print("\ndifferenza appaiata Tumor - PBMC, sugli STESSI cloni:")
    piv = j.pivot_table(
        index=["patient", "clone_b"], columns="tissue", values="error", aggfunc="max"
    )
    if {"Tumor", "PBMC"}.issubset(piv.columns):
        d = piv.dropna(subset=["Tumor", "PBMC"]).reset_index()
        d["diff"] = d.Tumor - d.PBMC
        by = {p: g["diff"].values for p, g in d.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        v = "IC esclude 0" if lo > 0 or hi < 0 else "IC include 0"
        print(f"  {m:+.3f}  IC95% [{lo:+.3f}, {hi:+.3f}]  {v}   (n={len(d)} cloni)")

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
    ax.bar(
        res.tissue,
        res.error_rate,
        yerr=[res.error_rate - res.lo, res.hi - res.error_rate],
        capsize=5,
        color="#d95f02",
        alpha=0.85,
    )
    ax.set_ylabel("tasso d'errore dell'etichetta")
    ax.set_title("chi sbaglia? riferimento = identita' clonale (marcatori)")
    fig.tight_layout()
    fig.savefig(FIG / "error_rate_by_compartment.png", dpi=140)
    plt.close(fig)


def partB(gex: pd.DataFrame) -> None:
    """Verifica se le cellule non-T con TCR presentano massa UMI e diluizione da doppietto fisico."""
    print("\n" + "=" * 78)
    print(
        "(B) DOPPIETTO FISICO o CELLULA DI CONFINE?  (misure fisiche non trascrittomiche)"
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
    out = pd.DataFrame(rows)
    print(out.round(2).to_string(index=False))
    out.to_csv(TAB / "doublet_physical.csv", index=False)

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
        print(
            f"  atteso per doppietto: {float(r_cd8.iloc[0]) * cd8 / max(nk_pos, 1):.2f}"
        )
        print(
            "  Un rapporto osservato nelle NK ridotto rispetto alle CD8T indica diluizione"
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
        description="Identita' clonale dai marcatori, errore per compartimento, doppietti."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tcrdir", type=Path, default=Path("data/raw/tcr"))
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

    lym = gex[gex.has_tcr & gex.all_celltype.isin(["CD4T", "CD8T", "NK"])].copy()
    partA(lym)
    partB(gex)


if __name__ == "__main__":
    main()

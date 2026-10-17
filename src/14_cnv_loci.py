#!/usr/bin/env python3
"""
14_cnv_loci.py - Score CNV sui loci ricorrenti del PDAC: separazione fra duttali tumorali e adiacenti.

Scopo
-----
Lo score CNV globale (`mean |X_cnv|`) aggrega l'intero genoma e risente dei
programmi trascrizionali tessuto-specifici. Questo modulo restringe il calcolo ai
loci con alterazioni ricorrenti nell'adenocarcinoma duttale pancreatico e misura
l'AUC tumore-vs-adiacente delle duttali, sia pooled sia entro paziente:
  - perdite attese (`direction = -1`): 9p21 (`CDKN2A`), 17p13 (`TP53`), 18q21 (`SMAD4`), 1p, 6q
  - guadagni attesi (`direction = +1`): 8q24 (`MYC`), 12p12 (`KRAS`), 20q
Calcola uno score direzionale per locus e un indice composito (somma degli z-score,
`PDAC_COMPOSITE`).

Come in `13_infercnv.py`, inferCNV gira una sola volta con il riferimento formato
dalle immunitarie di TUTTI i pazienti. Le coordinate delle finestre di `X_cnv` non
sono restituite da infercnvpy: sono interpolate in modo uniforme fra il primo e
l'ultimo gene annotato di ciascun cromosoma, quindi i confini dei loci sono
approssimati (una finestra puo' cadere dentro o fuori un locus per qualche Mb).
L'AUC pooled mescola differenze fra pazienti e fra tessuti; l'AUC entro paziente
(pazienti con almeno 10 duttali per tessuto) isola il contrasto tissutale.

Input
-----
  --h5ad  data/interim/sc_raw.h5ad
  cache   data/interim/cnv_matrix.npz (letta se esiste, altrimenti calcolata e scritta)

Output
------
  data/interim/cnv_matrix.npz                    matrice X_cnv e coordinate delle finestre (cache)
  results/tables/cnv_loci_auc.csv                AUC pooled per segnale (score medio, loci, composito)
  results/tables/cnv_loci_auc_by_patient.csv     AUC entro paziente di PDAC_COMPOSITE e dello score medio
                                                 (ultima riga `_mean`: media sui pazienti valutabili)
  results/figures/cnv_loci_separation.png, results/figures/cnv_loci_auc.png

Argomenti
---------
  --h5ad (obbligatorio), --method {infercnv,copykat,auto} (default auto)

Uso
---
    python src/14_cnv_loci.py --h5ad data/interim/sc_raw.h5ad
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
import scanpy as sc
from scipy.sparse import issparse
from sklearn.metrics import roc_auc_score

FIG = Path("results/figures")
TAB = Path("results/tables")
CACHE = Path("data/interim/cnv_matrix.npz")
SEED = 0

IMMUNE_REF = ["CD8T", "CD4T", "NK", "B", "Myeloid", "Plasma"]

# Coordinate genomiche (GRCh38/hg38, in bp) e direzionalita' attesa dei loci citogenetici PDAC
PDAC_LOCI = [
    {"name": "CDKN2A_9p21", "chrom": "9", "start": 20e6, "end": 24e6, "direction": -1},
    {"name": "TP53_17p13", "chrom": "17", "start": 0, "end": 10e6, "direction": -1},
    {"name": "SMAD4_18q21", "chrom": "18", "start": 45e6, "end": 60e6, "direction": -1},
    {"name": "1p_loss", "chrom": "1", "start": 0, "end": 60e6, "direction": -1},
    {"name": "6q_loss", "chrom": "6", "start": 90e6, "end": 170e6, "direction": -1},
    {"name": "MYC_8q24", "chrom": "8", "start": 120e6, "end": 145e6, "direction": +1},
    {"name": "KRAS_12p12", "chrom": "12", "start": 20e6, "end": 30e6, "direction": +1},
    {"name": "20q_gain", "chrom": "20", "start": 30e6, "end": 65e6, "direction": +1},
]


def add_positions(adata: ad.AnnData) -> bool:
    """Annota le coordinate cromosomiche dei geni tramite Ensembl BioMart."""
    import infercnvpy as cnv

    if "gene_symbol" in adata.var:
        adata.var_names = adata.var["gene_symbol"].astype(str).values
    adata.var_names_make_unique()
    try:
        cnv.io.genomic_position_from_biomart(
            adata, biomart_gene_id="hgnc_symbol", species="hsapiens"
        )
        return adata.var["chromosome"].notna().sum() > 5000
    except Exception as e:  # noqa: BLE001
        print(f"[!!] biomart fallito: {e}")
        return False


def compute_cnv(sub: ad.AnnData, method: str):
    """Calcola la matrice completa `X_cnv` e ricostruisce la mappa genomica (cromosoma, posizione) delle finestre."""
    import infercnvpy as cnv

    a = sub.copy()
    if method == "infercnv":
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        a = a[:, a.var["chromosome"].notna()].copy()
        ref = [c for c in IMMUNE_REF if c in a.obs["ref_group"].unique()]
        print(f"[..] infercnv, riferimento = {ref} di tutti i pazienti (lento) ...")
        cnv.tl.infercnv(
            a, reference_key="ref_group", reference_cat=ref, window_size=100
        )
    else:
        print("[..] copykat (lento) ...")
        cnv.tl.copykat(a, n_jobs=-1)

    X = a.obsm["X_cnv"]
    X = np.asarray(X.todense()) if issparse(X) else np.asarray(X)

    # Ricostruzione del cromosoma e della coordinata di ciascuna finestra di X_cnv.
    # infercnvpy restituisce solo la colonna d'inizio di ogni cromosoma (chr_pos), non le
    # coordinate delle finestre: la posizione e' interpolata in modo UNIFORME fra il primo
    # e l'ultimo gene annotato del cromosoma. I geni non sono equispaziati, quindi i
    # confini dei loci in PDAC_LOCI sono applicati in modo approssimato.
    chr_pos = a.uns["cnv"]["chr_pos"]
    n_win = X.shape[1]
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    win_chrom = np.empty(n_win, dtype=object)
    for i, (chrom, start_col) in enumerate(order):
        end_col = order[i + 1][1] if i + 1 < len(order) else n_win
        win_chrom[start_col:end_col] = chrom.replace("chr", "")
    var = a.var.copy()
    var["chromosome"] = var["chromosome"].astype(str).str.replace("chr", "")
    win_pos = np.full(n_win, np.nan)
    for chrom, start_col in order:
        cc = chrom.replace("chr", "")
        end_col = next((sc_ for ch2, sc_ in order if sc_ > start_col), n_win)
        genes = var[var.chromosome == cc].sort_values("start")
        if len(genes) == 0:
            continue
        pos = np.linspace(genes.start.min(), genes.start.max(), end_col - start_col)
        win_pos[start_col:end_col] = pos
    windows = pd.DataFrame({"chrom": win_chrom, "pos": win_pos})
    return X, windows, a.obs_names


def locus_score(X, windows, locus) -> np.ndarray:
    """Calcola lo score CNV medio moltiplicato per la direzione citogenetica attesa (+1 guadagno, -1 perdita)."""
    m = (
        (windows.chrom == locus["chrom"])
        & (windows.pos >= locus["start"])
        & (windows.pos <= locus["end"])
    )
    if m.sum() == 0:
        return np.full(X.shape[0], np.nan)
    return X[:, m.values].mean(axis=1) * locus["direction"]


def main() -> None:
    """Punto di ingresso principale per la valutazione citogenetica locus-specifica."""
    ap = argparse.ArgumentParser(
        description="Score CNV sui loci PDAC: AUC tumore vs adiacente, pooled ed entro paziente."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--method", choices=["infercnv", "copykat", "auto"], default="auto")
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    for c in ("tissue", "patients", "all_celltype"):
        adata.obs[c] = adata.obs[c].astype(str)

    is_duct = adata.obs.all_celltype == "Ductal"
    is_imm = adata.obs.all_celltype.isin(IMMUNE_REF)
    keep = (is_duct | is_imm) & adata.obs.tissue.isin(["Tumor", "Adjacent_normal"])
    sub = adata[keep].copy()
    sub.obs["ref_group"] = np.where(
        sub.obs.all_celltype.isin(IMMUNE_REF), sub.obs.all_celltype, "Ductal"
    )

    # Caricamento o calcolo con salvataggio compresso (.npz) della matrice delle finestre genomiche X_cnv
    if CACHE.exists():
        print(f"[ok] X_cnv da cache {CACHE}")
        npz = np.load(CACHE, allow_pickle=True)
        X, obs_names = npz["X"], npz["obs_names"]
        windows = pd.DataFrame({"chrom": npz["win_chrom"], "pos": npz["win_pos"]})
        sub = sub[obs_names].copy()
    else:
        method = args.method
        if method == "auto":
            method = "infercnv" if add_positions(sub) else "copykat"
            print(f"[ok] metodo: {method}")
        elif method == "infercnv" and not add_positions(sub):
            sys.exit("[ERRORE] posizioni assenti; usa --method copykat")
        t0 = time.time()
        X, windows, obs_names = compute_cnv(sub, method)
        print(f"[ok] X_cnv {X.shape} in {time.time() - t0:.0f}s")
        np.savez_compressed(
            CACHE,
            X=X,
            obs_names=np.array(obs_names),
            win_chrom=windows.chrom.values,
            win_pos=windows.pos.values,
        )
        sub = sub[obs_names].copy()

    obs = sub.obs
    duct = obs.all_celltype.values == "Ductal"
    tumor = duct & (obs.tissue.values == "Tumor")
    normal = duct & (obs.tissue.values == "Adjacent_normal")
    print(f"\nduttali: {tumor.sum():,} tumorali | {normal.sum():,} normali adiacenti")
    if tumor.sum() == 0 or normal.sum() == 0:
        raise SystemExit(
            "[ERRORE] servono duttali sia in Tumor sia in Adjacent_normal."
        )

    print("\n" + "=" * 74)
    print("SEPARAZIONE TUMORE vs NORMALE: score medio vs loci PDAC specifici")
    print("=" * 74)
    print("AUC = capacita' discriminativa tra duttali TUMORALI e NORMALI adiacenti.")
    print("0.5 = distribuzioni sovrapposte | >0.65 = separazione locus-specifica\n")

    def auc_tn(score):
        y = np.r_[np.ones(tumor.sum()), np.zeros(normal.sum())]
        idx = np.r_[np.where(tumor)[0], np.where(normal)[0]]
        s = score[idx]
        ok = ~np.isnan(s)
        return (
            roc_auc_score(y[ok], s[ok])
            if ok.sum() > 10 and len(np.unique(y[ok])) > 1
            else np.nan
        )

    mean_score = np.abs(X).mean(axis=1)
    rows = [
        {
            "signal": "MEAN_|CNV|",
            "auc_tumor_vs_normal": auc_tn(mean_score),
            "n_windows": X.shape[1],
        }
    ]
    loci_scores = {}
    for L in PDAC_LOCI:
        s = locus_score(X, windows, L)
        loci_scores[L["name"]] = s
        nwin = int(
            (
                (windows.chrom == L["chrom"])
                & (windows.pos >= L["start"])
                & (windows.pos <= L["end"])
            ).sum()
        )
        rows.append(
            {"signal": L["name"], "auc_tumor_vs_normal": auc_tn(s), "n_windows": nwin}
        )

    # Somma standardizzata (z-score) dei contributi direzionali dei singoli loci PDAC
    valid = [s for s in loci_scores.values() if not np.all(np.isnan(s))]
    if not valid:
        raise SystemExit("[ERRORE] nessun locus PDAC coperto dalle finestre di X_cnv.")
    comp = np.nansum(
        [(s - np.nanmean(s)) / (np.nanstd(s) + 1e-9) for s in valid], axis=0
    )
    rows.append(
        {
            "signal": "PDAC_COMPOSITE",
            "auc_tumor_vs_normal": auc_tn(comp),
            "n_windows": len(valid),
        }
    )

    res = pd.DataFrame(rows)
    print(res.round(3).to_string(index=False))
    res.to_csv(TAB / "cnv_loci_auc.csv", index=False)

    mean_auc = res.loc[res.signal == "MEAN_|CNV|", "auc_tumor_vs_normal"].iloc[0]
    comp_auc = res.loc[res.signal == "PDAC_COMPOSITE", "auc_tumor_vs_normal"].iloc[0]
    print("\n" + "=" * 74)
    if comp_auc > mean_auc + 0.05 and comp_auc > 0.65:
        print(f"ESITO: i loci PDAC separano (AUC {comp_auc:.2f}) meglio dello score")
        print(
            f"  medio (AUC {mean_auc:.2f}). Il segnale citogenetico e' locus-specifico."
        )
    elif comp_auc < 0.6:
        print(f"ESITO: sovrapposizione anche sui loci PDAC (AUC {comp_auc:.2f}).")
        print(
            "  Le cellule duttali adiacenti condividono il profilo delle duttali tumorali"
        )
        print(
            "  (compatibile con effetto di campo peri-tumorale o risoluzione limitata dell'mRNA)."
        )
    else:
        print(
            f"ESITO: separazione parziale (loci {comp_auc:.2f} vs medio {mean_auc:.2f})."
        )
        print("  Segnale intermedio lungo il continuum epiteliale.")
    print("=" * 74)

    print(
        "\nAUC tumore-vs-normale PDAC_COMPOSITE, per paziente con abbastanza normali:"
    )
    pat = obs.patients.values
    by_pat = []
    for p in np.unique(pat[normal]):
        t = tumor & (pat == p)
        n = normal & (pat == p)
        if n.sum() < 10 or t.sum() < 10:
            print(
                f"  {p}: non valutabile (tumor n={t.sum()}, normal n={n.sum()}, minimo 10)"
            )
            continue
        idx = np.r_[np.where(t)[0], np.where(n)[0]]
        y = np.r_[np.ones(t.sum()), np.zeros(n.sum())]
        row = {"patient": p, "n_tumor": int(t.sum()), "n_normal": int(n.sum())}
        for name, score in (("PDAC_COMPOSITE", comp), ("MEAN_|CNV|", mean_score)):
            s = score[idx]
            ok = ~np.isnan(s)
            row[f"auc_{name}"] = (
                roc_auc_score(y[ok], s[ok]) if len(np.unique(y[ok])) > 1 else np.nan
            )
        by_pat.append(row)
        print(
            f"  {p}: AUC={row['auc_PDAC_COMPOSITE']:.3f}  "
            f"(tumor n={t.sum()}, normal n={n.sum()})"
        )
    by_pat = pd.DataFrame(by_pat)
    if not by_pat.empty:
        # riga di riepilogo: media delle AUC entro paziente
        mean_row = by_pat.drop(columns="patient").mean(numeric_only=True).to_dict()
        by_pat = pd.concat(
            [by_pat, pd.DataFrame([{"patient": "_mean", **mean_row}])],
            ignore_index=True,
        )
    by_pat.to_csv(TAB / "cnv_loci_auc_by_patient.csv", index=False)
    by_pat = by_pat[by_pat.patient != "_mean"]
    if not by_pat.empty:
        print(
            f"  AUC entro paziente PDAC_COMPOSITE: media {by_pat.auc_PDAC_COMPOSITE.mean():.3f}, "
            f"range [{by_pat.auc_PDAC_COMPOSITE.min():.3f}, {by_pat.auc_PDAC_COMPOSITE.max():.3f}] "
            f"su {len(by_pat)} pazienti"
        )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].hist(
        mean_score[tumor],
        bins=40,
        alpha=0.55,
        density=True,
        label="Tumor",
        color="#d95f02",
    )
    axes[0].hist(
        mean_score[normal],
        bins=40,
        alpha=0.55,
        density=True,
        label="Normal adj",
        color="#1b9e77",
    )
    axes[0].set_title(f"score MEDIO |CNV|  (AUC {mean_auc:.2f})")
    axes[0].set_xlabel("mean |CNV|")
    axes[0].legend(fontsize=8)

    axes[1].hist(
        comp[tumor], bins=40, alpha=0.55, density=True, label="Tumor", color="#d95f02"
    )
    axes[1].hist(
        comp[normal],
        bins=40,
        alpha=0.55,
        density=True,
        label="Normal adj",
        color="#1b9e77",
    )
    axes[1].set_title(f"score LOCI PDAC  (AUC {comp_auc:.2f})")
    axes[1].set_xlabel("score composito loci")
    axes[1].legend(fontsize=8)
    fig.suptitle("Separabilita' sui loci citogenetici PDAC vs score medio")
    fig.tight_layout()
    fig.savefig(FIG / "cnv_loci_separation.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    d = (
        res[res.signal != "MEAN_|CNV|"]
        .set_index("signal")
        .auc_tumor_vs_normal.sort_values()
    )
    colors = ["#d95f02" if v > 0.6 else "#999999" for v in d.values]
    ax.barh(d.index, d.values, color=colors, alpha=0.85)
    ax.axvline(0.5, ls="--", c="k", lw=0.8)
    ax.axvline(
        mean_auc, ls=":", c="#d95f02", lw=1.2, label=f"score medio ({mean_auc:.2f})"
    )
    ax.set_xlabel("AUC tumore vs normale adiacente")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "cnv_loci_auc.png", dpi=140)
    plt.close(fig)
    print("\n[ok] figure: results/figures/cnv_loci_separation.png, cnv_loci_auc.png")


if __name__ == "__main__":
    main()

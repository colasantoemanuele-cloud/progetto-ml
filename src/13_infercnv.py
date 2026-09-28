#!/usr/bin/env python3
"""
13_infercnv.py - Inferenza delle variazioni del numero di copie (CNV) e classificazione LOPO delle cellule duttali.

Poiche' l'annotazione originale del dataset GSE278694 raggruppa tutte le cellule
epiteliali duttali in un'unica classe (`Ductal`) senza distinguere l'epitelio
neoplastico da quello residuo non trasformato, questo modulo implementa un
protocollo in tre stadi:

1. Stima del profilo di aberrazione cromosomica tramite `infercnvpy` (o `copykat`):
   Calcola lo scostamento di espressione lungo finestre genomiche di 100 geni
   utilizzando come riferimento diploide interno le popolazioni immunitarie
   (`CD8T`, `CD4T`, `NK`, `B`, `Myeloid`, `Plasma`) del medesimo paziente.
   La soglia individuale di aneuploidia e' definita come:
       soglia(paziente) = mediana(CNV_immuni) + 3 * MAD_normalizzata(CNV_immuni)

2. Classificazione supervisionata Leave-One-Patient-Out (LOPO):
   Addestra un classificatore Random Forest bilanciato sulle cellule duttali
   del tessuto tumorale dei pazienti con almeno 150 cellule duttali, valutando
   la predizione dello stato di aneuploidia sui pazienti esclusi dal training.

3. Validazione ortogonale sul tessuto adiacente normale (`Adjacent_normal`):
   Applica il modello addestrato esclusivamente sul compartimento tumorale alle
   cellule duttali del tessuto adiacente non neoplastico per verificare il tasso
   di specificita' su un compartimento anatomico indipendente da inferCNV.

Uso:
    python src/13_infercnv.py --h5ad data/interim/sc_raw.h5ad \
        --splits data/interim/splits_sc.csv
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
import scipy.sparse as sp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0

IMMUNE_REF = ["CD8T", "CD4T", "NK", "B", "Myeloid", "Plasma"]
# Supporto minimo di cellule duttali tumorali per l'inclusione del paziente nella validazione LOPO
MIN_DUCTAL = 150


def add_positions(adata: ad.AnnData) -> bool:
    """Interroga Ensembl BioMart per annotare coordinate cromosomiche (chromosome, start, end) in `adata.var`."""
    import infercnvpy as cnv

    adata.var_names = (
        adata.var["gene_symbol"].astype(str).values
        if "gene_symbol" in adata.var
        else adata.var_names
    )
    adata.var_names_make_unique()
    try:
        print("[..] posizioni genomiche da biomart ...")
        cnv.io.genomic_position_from_biomart(
            adata, biomart_gene_id="hgnc_symbol", species="hsapiens"
        )
        ok = adata.var["chromosome"].notna().sum()
        print(f"[ok] posizioni per {ok:,}/{adata.n_vars:,} geni")
        return ok > 5000
    except Exception as e:  # noqa: BLE001
        print(f"[!!] biomart fallito: {e}")
        return False


def run_infercnv(adata: ad.AnnData) -> np.ndarray:
    """Esegue infercnvpy su dati log-normalizzati e restituisce lo score medio |X_cnv| per cellula."""
    import infercnvpy as cnv

    a = adata.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a = a[:, a.var["chromosome"].notna()].copy()

    ref = [c for c in IMMUNE_REF if c in a.obs["ref_group"].unique()]
    print(f"[..] infercnv: riferimento = {ref}")
    cnv.tl.infercnv(a, reference_key="ref_group", reference_cat=ref, window_size=100)
    score = np.asarray(np.abs(a.obsm["X_cnv"]).mean(axis=1)).ravel()
    return score


def run_copykat(adata: ad.AnnData) -> np.ndarray:
    """Esegue l'algoritmo CopyKAT sulle conte grezze e restituisce lo score medio |X_cnv| per cellula."""
    import infercnvpy as cnv

    a = adata.copy()
    print("[..] copykat (lento) ...")
    cnv.tl.copykat(a, n_jobs=-1)
    score = np.asarray(np.abs(a.obsm["X_cnv"]).mean(axis=1)).ravel()
    return score


def classify_lopo(X, y, groups) -> pd.DataFrame:
    """Esegue la cross-validation Leave-One-Patient-Out con RandomForestClassifier bilanciato."""
    rows, preds = [], []
    for g in np.unique(groups):
        tr, te = groups != g, groups == g
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 1:
            continue
        clf = RandomForestClassifier(
            n_estimators=300,
            max_features="sqrt",
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=SEED,
        )
        clf.fit(X[tr], y[tr])
        p = clf.predict(X[te])
        proba = (
            clf.predict_proba(X[te])[:, 1]
            if len(clf.classes_) == 2
            else np.zeros(te.sum())
        )
        preds.append(
            pd.DataFrame({"patient": g, "y_true": y[te], "y_pred": p, "proba": proba})
        )
        row = {
            "patient": g,
            "n_test": int(te.sum()),
            "balanced_acc": balanced_accuracy_score(y[te], p)
            if len(np.unique(y[te])) > 1
            else np.nan,
            "f1_malig": f1_score(y[te], p, pos_label=1, zero_division=0),
        }
        if len(np.unique(y[te])) > 1:
            try:
                row["auc"] = roc_auc_score(y[te], proba)
            except ValueError:
                row["auc"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows), pd.concat(preds) if preds else pd.DataFrame()


def main() -> None:
    """Punto di ingresso principale per l'inferenza CNV, la sogliatura robusta e il test LOPO."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--splits", type=Path, required=True)
    ap.add_argument("--method", choices=["infercnv", "copykat", "auto"], default="auto")
    ap.add_argument(
        "--cnv-cache", type=Path, default=Path("data/interim/cnv_scores.csv")
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    for c in ("tissue", "patients", "all_celltype"):
        adata.obs[c] = adata.obs[c].astype(str)

    # Selezione delle cellule duttali (Tumor + Adjacent_normal) e del riferimento immunitario diploide
    is_duct = adata.obs.all_celltype == "Ductal"
    is_imm = adata.obs.all_celltype.isin(IMMUNE_REF)
    keep = (is_duct | is_imm) & adata.obs.tissue.isin(["Tumor", "Adjacent_normal"])
    sub = adata[keep].copy()
    sub.obs["ref_group"] = np.where(
        sub.obs.all_celltype.isin(IMMUNE_REF), sub.obs.all_celltype, "Ductal"
    )
    print(
        f"[ok] sottoinsieme: {sub.n_obs:,} cellule "
        f"({(sub.obs.all_celltype == 'Ductal').sum():,} duttali, "
        f"{is_imm[keep].sum():,} immuni di riferimento)"
    )

    duct = sub.obs[sub.obs.all_celltype == "Ductal"]
    print("\nduttali per paziente x tessuto:")
    print(pd.crosstab(duct.patients, duct.tissue, margins=True).to_string())

    # Calcolo o caricamento da cache dello score CNV globale per cellula
    if args.cnv_cache.exists():
        print(f"\n[ok] CNV score da cache {args.cnv_cache}")
        cs = pd.read_csv(args.cnv_cache, index_col=0)
        sub.obs["cnv_score"] = cs.reindex(sub.obs_names).cnv_score.values
    else:
        method = args.method
        if method == "auto":
            method = "infercnv" if add_positions(sub) else "copykat"
            print(f"[ok] metodo scelto: {method}")
        elif method == "infercnv":
            if not add_positions(sub):
                sys.exit("[ERRORE] posizioni non disponibili; usa --method copykat")
        t0 = time.time()
        score = run_infercnv(sub) if method == "infercnv" else run_copykat(sub)
        sub.obs["cnv_score"] = score
        print(f"[ok] CNV in {time.time() - t0:.0f}s")
        sub.obs[["cnv_score"]].to_csv(args.cnv_cache)

    # Sogliatura statistica intra-paziente: mediana + 3 * MAD normalizzata sugli immuni del donatore
    print("\n" + "=" * 74)
    print("ETICHETTA DI MALIGNITA' (soglia per paziente sugli immuni di riferimento)")
    print("=" * 74)
    sub.obs["malignant"] = 0
    rows = []
    for pat, g in sub.obs.groupby("patients"):
        imm = g[g.all_celltype.isin(IMMUNE_REF)].cnv_score
        if len(imm) < 20:
            continue
        thr = imm.median() + 3 * (imm - imm.median()).abs().median() * 1.4826
        dmask = (sub.obs.patients == pat) & (sub.obs.all_celltype == "Ductal")
        sub.obs.loc[dmask, "malignant"] = (
            sub.obs.loc[dmask, "cnv_score"] > thr
        ).astype(int)
        dd = sub.obs[dmask]
        for tis in ["Tumor", "Adjacent_normal"]:
            t = dd[dd.tissue == tis]
            if len(t):
                rows.append(
                    {
                        "patient": pat,
                        "tissue": tis,
                        "n": len(t),
                        "malignant": int(t.malignant.sum()),
                        "frac_malig": t.malignant.mean(),
                    }
                )
    lab = pd.DataFrame(rows)
    print(lab.round(3).to_string(index=False))
    lab.to_csv(TAB / "malignancy_labels.csv", index=False)

    # Identificazione di eventuali pazienti privi di contrasto binario interno
    tum = lab[lab.tissue == "Tumor"]
    degenerate = tum[(tum.frac_malig == 0) | (tum.frac_malig == 1)]
    if len(degenerate):
        print(
            f"\n[!!] {len(degenerate)} pazienti con duttali tumorali tutte uguali "
            f"(0% o 100% maligne): {degenerate.patient.tolist()}"
        )
        print("     non forniscono contrasto binario intra-paziente.")

    # Addestramento e valutazione Leave-One-Patient-Out sulle cellule duttali del tessuto tumorale
    print("\n" + "=" * 74)
    print("CLASSIFICATORE LOPO: predire la malignita' dal trascrittoma")
    print("=" * 74)
    train = sub[(sub.obs.all_celltype == "Ductal") & (sub.obs.tissue == "Tumor")].copy()
    counts = train.obs.patients.value_counts()
    good = counts[counts >= MIN_DUCTAL].index
    train = train[train.obs.patients.isin(good)].copy()
    print(
        f"pazienti con >={MIN_DUCTAL} duttali tumorali: {len(good)} -> {sorted(good)}"
    )
    print(
        f"cellule di training: {train.n_obs:,} | maligne: {train.obs.malignant.mean():.1%}"
    )

    # Normalizzazione e selezione dei 2.000 geni ad alta variabilita' (Seurat flavor)
    sc.pp.normalize_total(train, target_sum=1e4)
    sc.pp.log1p(train)
    sc.pp.highly_variable_genes(train, n_top_genes=2000, flavor="seurat")
    Xtr = sp.csr_matrix(train[:, train.var.highly_variable].X)
    ytr = train.obs.malignant.values
    gtr = train.obs.patients.values

    res, preds = classify_lopo(Xtr, ytr, gtr)
    print("\nper paziente (leave-one-patient-out):")
    print(res.round(3).to_string(index=False))
    if not res.empty:
        print(
            f"\nbalanced accuracy media: {res.balanced_acc.mean():.3f} ± {res.balanced_acc.std():.3f}"
        )
        if "auc" in res:
            print(f"AUC media: {res.auc.mean():.3f} ± {res.auc.std():.3f}")
    res.to_csv(TAB / "malignancy_lopo.csv", index=False)

    # Validazione esterna sulle cellule duttali del tessuto adiacente non neoplastico
    print("\n" + "=" * 74)
    print("TEST ESTERNO: duttali del tessuto ADIACENTE NORMALE (mai in training)")
    print("=" * 74)
    ext = sub[
        (sub.obs.all_celltype == "Ductal") & (sub.obs.tissue == "Adjacent_normal")
    ].copy()
    if ext.n_obs > 20:
        clf = RandomForestClassifier(
            n_estimators=300,
            max_features="sqrt",
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=SEED,
        )
        clf.fit(Xtr, ytr)
        sc.pp.normalize_total(ext, target_sum=1e4)
        sc.pp.log1p(ext)
        Xext = sp.csr_matrix(ext[:, train.var.highly_variable].X)
        pred = clf.predict(Xext)
        print(f"duttali normali adiacenti: {ext.n_obs:,}")
        print(
            f"predette MALIGNE (atteso ~0 in assenza di campo di cancerizzazione): {pred.mean():.1%}"
        )
        print(f"predette non-maligne: {(1 - pred.mean()):.1%}")
        print(
            "\nIl tasso di predizione sul tessuto adiacente normale misura la separabilita'"
        )
        print("tra epitelio duttale tumorale ed epitelio duttale peri-tumorale.")
        pd.DataFrame({"pred": pred}, index=ext.obs_names).to_csv(
            TAB / "external_test_predictions.csv"
        )
    else:
        print("numero insufficiente di duttali normali per il test esterno.")

    # Generazione della figura diagnostica (distribuzione CNV score e AUC per paziente)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for tis, c in [("Tumor", "#d95f02"), ("Adjacent_normal", "#1b9e77")]:
        d = sub.obs[(sub.obs.all_celltype == "Ductal") & (sub.obs.tissue == tis)]
        axes[0].hist(
            d.cnv_score,
            bins=40,
            alpha=0.6,
            label=f"Ductal {tis}",
            color=c,
            density=True,
        )
    imm = sub.obs[sub.obs.all_celltype.isin(IMMUNE_REF)]
    axes[0].hist(
        imm.cnv_score,
        bins=40,
        alpha=0.4,
        label="Immuni (rif)",
        color="gray",
        density=True,
    )
    axes[0].set_xlabel("CNV score")
    axes[0].set_ylabel("densità")
    axes[0].legend(fontsize=8)
    axes[0].set_title("distribuzione dello score CNV")
    if not res.empty and "auc" in res:
        axes[1].bar(res.patient, res.auc, color="#7570b3", alpha=0.85)
        axes[1].axhline(0.5, ls="--", c="k", lw=0.8)
        axes[1].set_ylabel("AUC (LOPO)")
        axes[1].tick_params(axis="x", rotation=90)
        axes[1].set_title("classificazione malignità per paziente")
    fig.tight_layout()
    fig.savefig(FIG / "malignancy.png", dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: results/figures/malignancy.png")  # noqa: F541


if __name__ == "__main__":
    main()

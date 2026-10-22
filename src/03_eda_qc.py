#!/usr/bin/env python3
"""
03_eda_qc.py - Controllo di qualita' descrittivo per libreria, composizione e supporto per paziente.

Scopo
-----
Legge l'h5ad in modalita' backed (senza caricare X in memoria) e produce tabelle e
figure descrittive. Non applica soglie di QC e non filtra cellule: le cellule usate a
valle sono quelle pubblicate dagli autori. In una validazione per paziente una classe
e' valutabile solo se e' presente in abbastanza pazienti: per questo lo script conta
anche i pazienti che supportano ciascuna classe con almeno 1, 10, 50, 100 cellule.

Input
-----
  --h5ad  data/interim/{sc,sn}_raw.h5ad (da 02_build_h5ad.py)

Output
------
  results/tables/dataset_summary_{tag}.csv   cellule, geni, librerie, pazienti; cellule e pazienti per tessuto
  results/tables/qc_by_library_{tag}.csv     mediane di n_genes, n_counts, pct_mt per libreria
  results/tables/composition_{tag}.csv       tipo cellulare x compartimento
  results/tables/patient_support_{tag}.csv   pazienti con almeno 1/10/50/100 cellule per classe
  results/figures/qc_by_library_{tag}.png, qc_complexity_{tag}.png,
  composition_by_tissue_{tag}.png, patient_by_celltype_{tag}.png

Argomenti
---------
  --h5ad (obbligatorio), --tag (default sc)

Uso
---
    python src/03_eda_qc.py --h5ad data/interim/sc_raw.h5ad --tag sc
    python src/03_eda_qc.py --h5ad data/interim/sn_raw.h5ad --tag sn
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

FIG = Path("results/figures")
TAB = Path("results/tables")

EXPECTED_LIBRARIES = {"sc": 29, "sn": 8}
EXPECTED_PATIENTS = {"sc": 14, "sn": 8}


def check_counts(adata: ad.AnnData) -> None:
    """Ispeziona un blocco di 2.000 righe della matrice backed per verificare la natura intera delle conte UMI."""
    try:
        x = adata.X[:2000]
        d = x.data if hasattr(x, "data") else np.asarray(x).ravel()
        is_int = bool(np.allclose(d, np.round(d)))
        print(f"[qc] X sembra {'CONTE INTERE' if is_int else 'GIA NORMALIZZATA (!)'}")
        print(f"[qc] min={d.min():.3g} max={d.max():.3g} media={d.mean():.3g}")
        if not is_int:
            print(
                "[!!] X non sembra fatta di conteggi interi: LogCPM, HVGSelector e inferCNV "
                "a valle assumono conteggi grezzi. (03 non applica soglie di QC.)"
            )
    except Exception as e:  # noqa: BLE001
        print(f"[!] impossibile ispezionare X: {e}")


def validate(obs: pd.DataFrame, tag: str) -> None:
    """Verifica che il numero di librerie e di soggetti estratti coincida con le specifiche della coorte."""
    n_lib, n_pat = obs.library.nunique(), obs.patients.nunique()
    print(f"[chk] librerie={n_lib} pazienti={n_pat}")
    if tag in EXPECTED_LIBRARIES:
        assert n_lib == EXPECTED_LIBRARIES[tag], (
            f"{n_lib} librerie, attese {EXPECTED_LIBRARIES[tag]}. Parsing degli ID rotto: "
            "rigenera l'h5ad con 02_build_h5ad.py corretto."
        )
        assert n_pat == EXPECTED_PATIENTS[tag], (
            f"{n_pat} pazienti, attesi {EXPECTED_PATIENTS[tag]}."
        )


def dataset_summary(obs: pd.DataFrame, n_genes: int, tag: str) -> None:
    """Numeri descrittivi del dataset: cellule, geni, librerie, pazienti e cellule per tessuto."""
    row = {
        "n_cells": len(obs),
        "n_genes": n_genes,
        "n_libraries": obs.library.nunique(),
        "n_patients": obs.patients.nunique(),
    }
    if "tissue" in obs:
        for t, n in obs.tissue.astype(str).value_counts().items():
            row[f"n_cells_{t}"] = int(n)
            row[f"n_patients_{t}"] = int(
                obs.loc[obs.tissue.astype(str) == t, "patients"].nunique()
            )
    pd.DataFrame([row]).to_csv(TAB / f"dataset_summary_{tag}.csv", index=False)
    print("[ok] riepilogo: " + ", ".join(f"{k}={v:,}" for k, v in row.items()))


def qc_tables(obs: pd.DataFrame, tag: str) -> None:
    """Calcola e salva le statistiche di qualita (mediana di geni, conte UMI e percentuale mitocondriale) per libreria."""
    by_lib = obs.groupby("library").agg(
        n_cells=("n_genes", "size"),
        median_genes=("n_genes", "median"),
        median_counts=("n_counts", "median"),
        median_pct_mt=("pct_mt", "median"),
        p95_pct_mt=("pct_mt", lambda s: s.quantile(0.95)),
    )
    by_lib.sort_values("n_cells", ascending=False).to_csv(
        TAB / f"qc_by_library_{tag}.csv"
    )

    print(
        f"\n[qc] geni/cellula (mediana per libreria): "
        f"{by_lib.median_genes.min():.0f} - {by_lib.median_genes.max():.0f}"
    )
    print(
        f"[qc] pct_mt (mediana per libreria): "
        f"{by_lib.median_pct_mt.min():.2f} - {by_lib.median_pct_mt.max():.2f}"
    )
    print("\n[qc] librerie con pct_mt piu alto:")
    print(
        by_lib.sort_values("median_pct_mt", ascending=False)
        .head(3)[["n_cells", "median_genes", "median_pct_mt"]]
        .to_string()
    )


def qc_figures(obs: pd.DataFrame, tag: str) -> None:
    """Genera i grafici diagnostici di distribuzione QC per libreria e il diagramma di complessita UMI vs geni."""
    order = sorted(obs.library.unique())
    fig, axes = plt.subplots(
        3, 1, figsize=(max(10, 0.45 * len(order)), 11), sharex=True
    )
    for ax, col, log in zip(
        axes, ["n_genes", "n_counts", "pct_mt"], [True, True, False]
    ):
        sns.violinplot(
            data=obs, x="library", y=col, order=order, ax=ax, cut=0, linewidth=0.4
        )
        if log:
            ax.set_yscale("log")
        ax.set_xlabel("")
        ax.set_ylabel(col)
    axes[-1].tick_params(axis="x", rotation=90)
    fig.suptitle(f"QC per libreria: {tag}")
    fig.tight_layout()
    fig.savefig(FIG / f"qc_by_library_{tag}.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sub = obs.sample(min(40000, len(obs)), random_state=0)
    sc_ = ax.scatter(
        sub.n_counts,
        sub.n_genes,
        s=1,
        alpha=0.15,
        c=sub.pct_mt,
        cmap="viridis",
        vmin=0,
        vmax=25,
    )
    fig.colorbar(sc_, ax=ax, label="pct_mt")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("n_counts")
    ax.set_ylabel("n_genes")
    ax.set_title(f"complessita: {tag}")
    fig.tight_layout()
    fig.savefig(FIG / f"qc_complexity_{tag}.png", dpi=140)
    plt.close(fig)


def composition(obs: pd.DataFrame, tag: str) -> None:
    """Costruisce le tabelle e le heatmap di abbondanza cellulare incrociata per tessuto e per paziente."""
    ct = pd.crosstab(obs.all_celltype, obs.tissue)
    ct["total"] = ct.sum(axis=1)
    ct = ct.sort_values("total", ascending=False)
    ct.to_csv(TAB / f"composition_{tag}.csv")
    print(f"\n[eda] composizione celltype x tissue\n{ct.to_string()}")

    frac = pd.crosstab(obs.tissue, obs.all_celltype, normalize="index")
    fig, ax = plt.subplots(figsize=(11, 4))
    frac.plot(kind="barh", stacked=True, ax=ax, width=0.7, colormap="tab20")
    ax.set_xlabel("frazione di cellule")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_title(f"composizione per tessuto: {tag}")
    fig.tight_layout()
    fig.savefig(FIG / f"composition_by_tissue_{tag}.png", dpi=140)
    plt.close(fig)

    m = pd.crosstab(obs.patients, obs.all_celltype)
    fig, ax = plt.subplots(figsize=(12, 5))
    sns.heatmap(
        np.log10(m + 1),
        annot=m.values,
        fmt="d",
        annot_kws={"size": 6},
        cmap="Blues",
        ax=ax,
        cbar_kws={"label": "log10(n+1)"},
    )
    ax.set_title(f"cellule per paziente x tipo: {tag}")
    fig.tight_layout()
    fig.savefig(FIG / f"patient_by_celltype_{tag}.png", dpi=140)
    plt.close(fig)


def patient_support(obs: pd.DataFrame, tag: str, thresholds=(1, 10, 50, 100)) -> None:
    """Quantifica per ogni popolazione cellulare il numero di soggetti che superano soglie crescenti di numerosita."""
    m = pd.crosstab(obs.patients, obs.all_celltype)
    sup = pd.DataFrame({f"n_patients_ge_{k}": (m >= k).sum(axis=0) for k in thresholds})
    sup["n_cells"] = m.sum(axis=0)
    sup["n_patients_total"] = obs.patients.nunique()
    sup = sup.sort_values("n_cells", ascending=False)
    sup.to_csv(TAB / f"patient_support_{tag}.csv")

    print(f"\n[split] supporto per paziente (n pazienti = {obs.patients.nunique()})")
    print(sup.to_string())

    # Classificazione dell'idoneita di ciascuna classe rispetto agli schemi di partizionamento cross-patient.
    print("\n[split] valutazione di supporto per classe:")
    for cls, r in sup.iterrows():
        n = int(r["n_patients_ge_10"])
        if n >= 5:
            v = "OK per StratifiedGroupKFold a 5 fold"
        elif n >= 3:
            v = "fragile: solo LOPO, riporta la varianza"
        else:
            v = "NON APPRENDIBILE cross-patient -> accorpa o escludi"
        print(
            f"  {cls:14s} {int(r['n_cells']):>7,} cellule | {n:>2} pazienti con >=10 | {v}"
        )


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Mostra il default in --help solo quando esiste (non per argomenti obbligatori o None)."""

    def _get_help_string(self, action):
        """Omette il default dall'help per argomenti obbligatori o con default None."""
        if action.default is None or action.required:
            return action.help
        return super()._get_help_string(action)


def main() -> None:
    """Apre l'archivio AnnData in modalita backed='r' (sola lettura dei metadati) ed esegue l'intera suite QC/EDA."""
    ap = argparse.ArgumentParser(
        formatter_class=_HelpFormatter,
        description="QC descrittivo per libreria, composizione e supporto per paziente (nessun filtro).",
    )
    ap.add_argument(
        "--h5ad", type=Path, required=True, help="AnnData prodotto da 02_build_h5ad.py"
    )
    ap.add_argument("--tag", default="sc", help="suffisso dei file di output (sc o sn)")
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    # Apertura in modalita backed='r' per evitare il caricamento in RAM della matrice sparsa X.
    adata = ad.read_h5ad(args.h5ad, backed="r")
    print(f"[ok] {args.h5ad.name}: {adata.n_obs:,} cellule x {adata.n_vars:,} geni")
    check_counts(adata)

    obs = adata.obs.copy()
    obs["all_celltype"] = obs.all_celltype.astype(str)
    obs["library"] = obs.library.astype(str)

    validate(obs, args.tag)
    dataset_summary(obs, adata.n_vars, args.tag)
    qc_tables(obs, args.tag)
    qc_figures(obs, args.tag)
    composition(obs, args.tag)
    patient_support(obs, args.tag)

    print(f"\n[ok] tabelle in {TAB}/  figure in {FIG}/")


if __name__ == "__main__":
    main()

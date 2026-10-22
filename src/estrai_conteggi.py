#!/usr/bin/env python3
"""
estrai_conteggi.py - Cellule per paziente della coorte appaiata (Tumor + Adjacent_normal).

Scopo
-----
Utilita' di riepilogo: legge i metadati cellulari scRNA-seq di GEO (GSE278688),
esclude Platelet ed Endocrine come 04_splits.py, individua i pazienti che hanno sia
Tumor sia Adjacent_normal e stampa il numero di cellule di ciascuno nei due tessuti
solidi (sono le cellule del task task1a_paired; totale 67.662 nella versione
pubblicata). Non legge la matrice di espressione.

Input
-----
  --metadata  data/raw/GSE278688/GSE278688_sc_metadata.csv.gz
              (colonne usate: tissue, patients, all_celltype)

Output
------
  solo testo a video (nessun file scritto)

Argomenti
---------
  --metadata (default data/raw/GSE278688/GSE278688_sc_metadata.csv.gz)

Uso
---
    python src/estrai_conteggi.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Popolazioni escluse a priori, come in 04_splits.py
EXCLUDE_ALWAYS = ["Platelet", "Endocrine"]


def paired_counts(obs: pd.DataFrame) -> pd.Series:
    """Cellule per paziente in Tumor + Adjacent_normal, per i soli pazienti appaiati."""
    obs = obs[~obs["all_celltype"].isin(EXCLUDE_ALWAYS)]
    tis = pd.crosstab(obs["patients"], obs["tissue"])
    paired = sorted(
        tis.index[(tis.get("Tumor", 0) > 0) & (tis.get("Adjacent_normal", 0) > 0)]
    )
    t1 = obs[obs["tissue"].isin(["Tumor", "Adjacent_normal"])]
    return t1[t1["patients"].isin(paired)]["patients"].value_counts().sort_index()


def main() -> None:
    """Legge i metadati, controlla le colonne e stampa le cellule per paziente
    della coorte appaiata con il totale."""
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Cellule per paziente della coorte appaiata Tumor + Adjacent_normal.",
    )
    ap.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/raw/GSE278688/GSE278688_sc_metadata.csv.gz"),
        help="metadati per cellula di GEO (paziente, tessuto)",
    )
    args = ap.parse_args()
    if not args.metadata.exists():
        raise SystemExit(f"[ERRORE] file assente: {args.metadata}")

    print("Caricamento metadati...")
    obs = pd.read_csv(args.metadata)
    missing = {"tissue", "patients", "all_celltype"} - set(obs.columns)
    if missing:
        raise SystemExit(f"[ERRORE] colonne mancanti nei metadati: {sorted(missing)}")
    conteggi = paired_counts(obs)
    if conteggi.empty:
        raise SystemExit("[ERRORE] nessun paziente con Tumor e Adjacent_normal.")

    print("\n--- RIEPILOGO CONTEGGI PER PAZIENTE (COORTE APPAIATA) ---")
    print("Paziente\tCellule totali")
    print("-" * 30)
    for paziente, n_celle in conteggi.items():
        print(f"{paziente}\t\t{n_celle:,}".replace(",", "."))
    print("-" * 30)
    print(f"TOTALE\t\t{conteggi.sum():,}".replace(",", "."))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
15_inventory.py - Inventario delle colonne di obs e dei file VDJ, con confondimento rispetto al paziente.

Scopo
-----
Documenta, prima di qualsiasi analisi secondaria, tutte le variabili disponibili:
  1. per ogni colonna di `adata.obs`: tipo, cardinalita', valori piu' frequenti (o
     intervallo per le numeriche) e se e' fra le colonne considerate dalla pipeline
     (`already_used`: elenco scritto a mano `KNOWN_USED`; comprende anche colonne
     tecniche come `barcode`, `probe` e `pct_mt`, non solo variabili di analisi);
  2. per le colonne categoriali con <= 30 livelli, il confondimento con il paziente:
     frazione dei livelli osservati in un solo paziente (1.0 = la variabile coincide
     con il paziente e non e' confrontabile fra pazienti);
  3. le colonne dei file di contig VDJ, usate e non usate.
Le colonne categoriali (dtype `category`) sono gestite con
`pd.api.types.is_float_dtype` / `is_numeric_dtype` (prima `np.issubdtype` andava in
errore su di esse).

Input
-----
  --h5ad    data/interim/sc_raw.h5ad (solo obs, modalita' backed)
  --tcrdir  cartella con i file VDJ (si legge l'intestazione del primo)

Output
------
  results/tables/obs_inventory.csv   una colonna di obs per riga (colonne tutte NaN: n_unique = 0)
  results/tables/vdj_inventory.csv   un file VDJ per riga: GSM, paziente, tessuto, numero di contig

Argomenti
---------
  --h5ad (obbligatorio), --tcrdir (default data/raw/tcr)

Uso
---
    python src/15_inventory.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

TAB = Path("results/tables")
# Colonne di obs considerate dalla pipeline (elenco scritto a mano, da aggiornare con gli script)
KNOWN_USED = {
    "tissue",
    "patients",
    "all_celltype",
    "library",
    "barcode",
    "probe",
    "n_counts",
    "n_genes",
    "pct_mt",
    "cnv_score",
    "malignant",
    "ref_group",
}


def confounding_with_patient(obs: pd.DataFrame, col: str) -> float:
    """Calcola la frazione di categorie di `col` osservate in un solo paziente (1.0 = confondimento totale)."""
    ct = pd.crosstab(obs[col].astype(str), obs["patients"].astype(str))
    per_level_patients = (ct > 0).sum(axis=1)
    return float((per_level_patients == 1).mean())


def main() -> None:
    """Esegue la ricognizione delle colonne di `obs` e dei file VDJ e salva `obs_inventory.csv`."""
    ap = argparse.ArgumentParser(
        description="Inventario delle colonne di obs e dei file VDJ."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tcrdir", type=Path, default=Path("data/raw/tcr"))
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs
    if obs.shape[0] == 0 or obs.shape[1] == 0:
        raise SystemExit(f"[ERRORE] obs vuoto in {args.h5ad}")

    print("=" * 78)
    print(f"INVENTARIO obs: {obs.shape[0]:,} cellule x {obs.shape[1]} colonne")
    print("=" * 78)

    rows = []
    for col in obs.columns:
        s = obs[col]
        card = int(s.nunique(dropna=True))
        dtype = str(s.dtype)
        used = col in KNOWN_USED
        info = {"column": col, "dtype": dtype, "n_unique": card, "already_used": used}
        if card == 0:
            info["confound_patient"] = np.nan
            info["top_values"] = "colonna vuota (tutti i valori mancanti)"
        elif card <= 30 and not pd.api.types.is_float_dtype(s.dtype):
            info["confound_patient"] = (
                confounding_with_patient(obs, col) if "patients" in obs else np.nan
            )
            vals = s.value_counts(dropna=False).head(8).to_dict()
            info["top_values"] = "; ".join(f"{k}={v}" for k, v in vals.items())
        else:
            info["confound_patient"] = np.nan
            if pd.api.types.is_numeric_dtype(s.dtype):
                info["top_values"] = (
                    f"range [{s.min():.3g}, {s.max():.3g}], median {s.median():.3g}"
                )
            else:
                info["top_values"] = f"{card} valori distinti (alta cardinalita')"
        rows.append(info)

    inv = pd.DataFrame(rows)
    inv.to_csv(TAB / "obs_inventory.csv", index=False)

    print("\n--- COLONNE CONSIDERATE dalla pipeline (elenco KNOWN_USED) ---")
    print(inv[inv.already_used][["column", "dtype", "n_unique"]].to_string(index=False))

    unused = inv[~inv.already_used]
    print(f"\n--- COLONNE NON CONSIDERATE ({len(unused)}) ---")
    if unused.empty:
        print("  nessuna: tutte le colonne di obs sono nell'elenco KNOWN_USED.")
    else:
        for _, r in unused.iterrows():
            print(f"\n  [{r.column}]  ({r['dtype']}, {r.n_unique} valori)")
            print(f"    {r.top_values}")
            if pd.notna(r.confound_patient):
                c = r.confound_patient
                verdict = (
                    "CONFUSA col paziente: inutilizzabile per associazioni"
                    if c > 0.8
                    else "parzialmente confusa col paziente: cautela"
                    if c > 0.4
                    else "NON confusa col paziente: potenzialmente utilizzabile"
                )
                print(f"    confondimento col paziente = {c:.2f} -> {verdict}")

    # Ricognizione degli attributi strutturali presenti nei file VDJ (es. segmenti genici V/D/J)
    print("\n" + "=" * 78)
    print("COLONNE DEI FILE TCR (usate: barcode, chain, cdr3_nt, umis, productive...)")
    print("=" * 78)
    f = next(
        iter(sorted(args.tcrdir.glob("*_filtered_contig_annotations.csv.gz"))), None
    )
    if f:
        with gzip.open(f, "rt") as fh:
            header = fh.readline().strip().split(",")
        used_tcr = {
            "barcode",
            "chain",
            "cdr3_nt",
            "cdr3",
            "umis",
            "productive",
            "is_cell",
            "high_confidence",
            "full_length",
            "raw_clonotype_id",
        }
        print("colonne presenti:")
        for h in header:
            mark = "  (usata)" if h in used_tcr else "  [non utilizzata]"
            print(f"    {h}{mark}")
        print(
            "\nLe colonne dei segmenti genici V/D/J (v_gene, j_gene) consentono analisi"
        )
        print(
            "descrittive del repertorio recettoriale subordinatamente alla validazione cross-patient."
        )
    else:
        print("nessun file TCR trovato.")

    # Elenco dei file VDJ: paziente, tessuto e numero di contig per file
    vdj = []
    for p in sorted(args.tcrdir.glob("*_filtered_contig_annotations.csv.gz")):
        m = re.match(r"(GSM\d+)_PA_(\d+)-(\w+)-VDJ", p.name)
        with gzip.open(p, "rt") as fh:
            n_contigs = sum(1 for _ in fh) - 1
        vdj.append(
            {
                "file": p.name,
                "gsm": m.group(1) if m else "",
                "patient": f"PA{int(m.group(2)):02d}" if m else "",
                "tissue": m.group(3) if m else "",
                "n_contigs": n_contigs,
            }
        )
    pd.DataFrame(vdj, columns=["file", "gsm", "patient", "tissue", "n_contigs"]).to_csv(
        TAB / "vdj_inventory.csv", index=False
    )
    print(f"\n[ok] file VDJ: {len(vdj)} -> {TAB / 'vdj_inventory.csv'}")

    print("\n" + "=" * 78)
    print("CRITERI DI VALIDAZIONE METODOLOGICA")
    print("=" * 78)
    print("""
  Una variabile non confusa con l'identita' del paziente rappresenta un'ipotesi
  esplorativa. Per costituire un'evidenza statistica robusta richiede:
    1. Formulazione a priori del meccanismo biologico.
    2. Validazione su partizioni o coorti indipendenti.
    3. Controllo rigoroso del tasso di falsi positivi per test multipli (Holm-Bonferroni).
""")


if __name__ == "__main__":
    main()

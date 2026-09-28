#!/usr/bin/env python3
"""
15_inventory.py - Censimento strutturato dei metadati osservazionali e verifica di confondimento.

Questo modulo ispeziona sistematicamente tutte le colonne presenti nella tabella
`adata.obs` e nelle annotazioni dei contig scTCR-seq con tre obiettivi:
  1. Documentare tipo di dato, cardinalita' e valori modali di ciascun attributo.
  2. Distinguere le variabili gia' impiegate nei task supervisionati da quelle
     disponibili per analisi secondarie.
  3. Quantificare per ogni variabile categoriale il grado di collinearita' con
     l'identita' del donatore (`confounding_with_patient`), definito come la
     frazione di livelli della variabile che compaiono in un unico paziente
     (dove un valore pari a 1.0 indica confondimento perfetto con il blocco
     paziente e conseguente non identificabilita' cross-patient).

Uso:
    python src/15_inventory.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

TAB = Path("results/tables")
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tcrdir", type=Path, default=Path("data/raw/tcr"))
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs

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
        if card <= 30 and not np.issubdtype(s.dtype, np.floating):
            info["confound_patient"] = (
                confounding_with_patient(obs, col) if "patients" in obs else np.nan
            )
            vals = s.value_counts(dropna=False).head(8).to_dict()
            info["top_values"] = "; ".join(f"{k}={v}" for k, v in vals.items())
        else:
            info["confound_patient"] = np.nan
            if np.issubdtype(s.dtype, np.number):
                info["top_values"] = (
                    f"range [{s.min():.3g}, {s.max():.3g}], median {s.median():.3g}"
                )
            else:
                info["top_values"] = f"{card} valori distinti (alta cardinalita')"
        rows.append(info)

    inv = pd.DataFrame(rows)
    inv.to_csv(TAB / "obs_inventory.csv", index=False)

    print("\n--- COLONNE GIA' USATE nell'analisi ---")
    print(inv[inv.already_used][["column", "dtype", "n_unique"]].to_string(index=False))

    unused = inv[~inv.already_used]
    print(f"\n--- COLONNE NON ANCORA USATE ({len(unused)}) ---")
    if unused.empty:
        print("  nessuna: ogni colonna di obs e' gia' entrata nell'analisi.")
    else:
        for _, r in unused.iterrows():
            print(f"\n  [{r.column}]  ({r.dtype}, {r.n_unique} valori)")
            print(f"    {r.top_values}")
            if not np.isnan(r.confound_patient):
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

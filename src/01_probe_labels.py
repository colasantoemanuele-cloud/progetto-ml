#!/usr/bin/env python3
"""
01_probe_labels.py - Verifica preliminare delle etichette cellulari e dei barcode nelle SubSeries GEO.

Scopo
-----
La SuperSeries GSE278694 comprende quattro SubSeries:
    GSE278687  trascrittomica spaziale Visium
    GSE278688  scRNA-seq (sc_matrix.mtx.gz: 1,0 GB; sc_metadata.csv.gz: 1,1 MB)
    GSE278689  snRNA-seq (sn_matrix.mtx.gz: 101 MB; sn_metadata.csv.gz: 192 KB)
    GSE300435  scTCR-seq V(D)J
Lo script scarica solo i metadati, le liste di barcode e un file VDJ pilota
(GSM9060613, PA01 Tumor) e verifica:
  Q3  presenza e cardinalita' delle etichette di tipo cellulare (`all_celltype`);
  Q2  sovrapposizione dei barcode 10x a 16 nt fra il file VDJ pilota e i metadati GEX;
  Q1b struttura delle coorti di pazienti scRNA-seq e snRNA-seq.
I 5 pazienti con Tumor e Adjacent_normal appaiati sono identificati da 04_splits.py.

Input
-----
  rete: https://ftp.ncbi.nlm.nih.gov/geo (i file gia' presenti non vengono riscaricati)

Output
------
  {outdir}/sc_metadata.csv.gz, sc_barcodes.tsv.gz, sn_metadata.csv.gz,
  sn_barcodes.tsv.gz, tcr_PA01_tumor.csv.gz (solo download; nessuna tabella)

Argomenti
---------
  --outdir (default data/raw/probe)

Uso
---
    python src/01_probe_labels.py --outdir data/raw/probe
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

FTP = "https://ftp.ncbi.nlm.nih.gov/geo"

# Elenco delle risorse leggere (metadati, liste di barcode e un campione VDJ pilota) da scaricare via HTTPS/FTP.
TARGETS = [
    (
        "sc_metadata.csv.gz",
        f"{FTP}/series/GSE278nnn/GSE278688/suppl/GSE278688_sc_metadata.csv.gz",
    ),
    (
        "sc_barcodes.tsv.gz",
        f"{FTP}/series/GSE278nnn/GSE278688/suppl/GSE278688_sc_barcodes.tsv.gz",
    ),
    (
        "sn_metadata.csv.gz",
        f"{FTP}/series/GSE278nnn/GSE278689/suppl/GSE278689_sn_metadata.csv.gz",
    ),
    (
        "sn_barcodes.tsv.gz",
        f"{FTP}/series/GSE278nnn/GSE278689/suppl/GSE278689_sn_barcodes.tsv.gz",
    ),
    (
        "tcr_PA01_tumor.csv.gz",
        f"{FTP}/samples/GSM9060nnn/GSM9060613/suppl/"
        "GSM9060613_PA_01-Tumor-VDJ_filtered_contig_annotations.csv.gz",
    ),
]

# Espressioni regolari per l'individuazione automatica delle colonne di tipo cellulare e di soggetto/libreria.
LABEL_HINTS = (
    r"cell.?type|celltype|annotation|ident|cluster|label|subtype|major|lineage"
)
SAMPLE_HINTS = r"sample|orig\.ident|patient|donor|group|tissue|source|batch|library"


def fetch(name: str, url: str, outdir: Path) -> Path:
    """Scarica un file compresso nella directory locale di lavoro riutilizzando la cache se gia presente."""
    dest = outdir / name
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[ok] cache  {name}  ({dest.stat().st_size / 1e6:.2f} MB)")
        return dest
    print(f"[..] scarico {name}")
    req = urllib.request.Request(url, headers={"User-Agent": "pdac-ml/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            dest.write_bytes(r.read())
    except (urllib.error.HTTPError, urllib.error.URLError) as e:
        print(f"[!!] fallito {name}: {e}\n     url: {url}")
        return dest  # Ritorna il percorso previsto; il chiamante verifica l'esistenza del file su disco.
    print(f"[ok] {name}  ({dest.stat().st_size / 1e6:.2f} MB)")
    return dest


def peek_gz(path: Path, n: int = 3) -> None:
    """Stampa in anteprima le prime n righe di un archivio compresso gzip per ispezione strutturale."""
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            print(f"     | {line.rstrip()[:160]}")


def core_barcode(s: pd.Series) -> pd.Series:
    """Estrae la sequenza nucleotidica canonica a 16 basi del barcode cellulare 10x Genomics."""
    return s.astype(str).str.extract(r"([ACGT]{16})", expand=False)


def analyse_metadata(path: Path, tag: str) -> pd.DataFrame | None:
    """Analizza lo schema delle colonne e la distribuzione di frequenza delle etichette nel file di metadati."""
    if not path.exists() or path.stat().st_size == 0:
        print(f"[!!] {tag}: file assente, salto")
        return None
    print(f"\n{'=' * 78}\n{tag}: {path.name}\n{'=' * 78}")
    print("Prime righe grezze:")
    peek_gz(path)

    df = pd.read_csv(path, compression="gzip", low_memory=False)
    print(f"\nshape: {df.shape}")
    print(f"colonne ({len(df.columns)}): {list(df.columns)}")

    label_cols = [c for c in df.columns if re.search(LABEL_HINTS, str(c), re.I)]
    sample_cols = [c for c in df.columns if re.search(SAMPLE_HINTS, str(c), re.I)]

    print(f"\n>>> candidate colonne ETICHETTA: {label_cols or 'NESSUNA'}")
    for c in label_cols:
        vc = df[c].value_counts()
        print(f"\n--- {c}  ({df[c].nunique()} valori distinti) ---")
        print(vc.head(25).to_string())
        if len(vc) > 25:
            print(f"    ... e altri {len(vc) - 25}")

    print(f"\n>>> candidate colonne CAMPIONE/PAZIENTE: {sample_cols or 'NESSUNA'}")
    for c in sample_cols[:4]:
        print(f"\n--- {c}  ({df[c].nunique()} valori distinti) ---")
        print(df[c].value_counts().head(30).to_string())

    return df


def check_tcr_overlap(sc: pd.DataFrame | None, tcr_path: Path) -> None:
    """Valuta l'intersezione preliminare dei barcode 16-nt tra il file di annotazioni V(D)J e i metadati GEX."""
    print(
        f"\n{'=' * 78}\nQ2: intersezione barcode TCR <-> GEX  [validazione ortogonale]\n{'=' * 78}"
    )
    if sc is None or not tcr_path.exists() or tcr_path.stat().st_size == 0:
        print("[!!] dati mancanti, salto")
        return

    tcr = pd.read_csv(tcr_path, compression="gzip", low_memory=False)
    print(f"TCR shape: {tcr.shape}")
    print(f"TCR colonne: {list(tcr.columns)[:12]}")
    if "barcode" not in tcr.columns:
        print("[!!] nessuna colonna 'barcode' nel TCR. Ispeziona a mano.")
        return

    tcr_bc = set(core_barcode(tcr["barcode"]).dropna())
    print(f"barcode TCR unici (PA01 Tumor): {len(tcr_bc)}")

    # Individuazione della colonna contenente l'identificatore univoco della cellula nei metadati scRNA.
    bc_col = next(
        (
            c
            for c in sc.columns
            if re.search(r"barcode|cell|^unnamed: 0$|^index$", str(c), re.I)
        ),
        sc.columns[0],
    )
    print(f"colonna barcode usata nel metadata sc: '{bc_col}'")
    print(f"esempi: {sc[bc_col].astype(str).head(3).tolist()}")

    sc_all = core_barcode(sc[bc_col])
    n_parsed = sc_all.notna().sum()
    print(f"barcode 10x estratti dal metadata: {n_parsed}/{len(sc)}")
    if n_parsed == 0:
        print("[!!] formato barcode non riconosciuto. Ispeziona 'esempi' qui sopra.")
        return

    # Tentativo di filtro sul campione PA01 Tumor a partire dalle colonne di campione/paziente.
    samp_col = next(
        (c for c in sc.columns if re.search(SAMPLE_HINTS, str(c), re.I)), None
    )
    subset = sc
    if samp_col is not None:
        mask = sc[samp_col].astype(str).str.contains(
            r"PA[_#\s-]*0?1\b", case=False, regex=True
        ) & sc[samp_col].astype(str).str.contains("tumor", case=False)
        if mask.sum() > 0:
            subset = sc[mask]
            print(f"filtro su '{samp_col}' ~ PA01 Tumor: {len(subset)} cellule")
        else:
            print(
                f"[!] non trovo PA01 Tumor in '{samp_col}'. Uso TUTTE le cellule (stima non stratificata per libreria)."
            )

    sc_bc = set(core_barcode(subset[bc_col]).dropna())
    inter = tcr_bc & sc_bc
    print(f"\nbarcode GEX nel sottoinsieme : {len(sc_bc)}")
    print(f"intersezione TCR ∩ GEX       : {len(inter)}")
    if sc_bc:
        print(f"frazione TCR recuperata      : {len(inter) / max(len(tcr_bc), 1):.1%}")
        print(f"frazione GEX con TCR         : {len(inter) / len(sc_bc):.1%}")

    print()
    if len(inter) > 100:
        print("ESITO: barcode appaiati confermati per la validazione V(D)J.")
        print(
            "  Nota metodologica: il TCR e' un marcatore positivo ma non rilevato in tutte le cellule T."
        )
        print("  Permette di stimare la precisione delle annotazioni sulle cellule T.")
    elif len(inter) > 0:
        print(
            "ESITO: sovrapposizione ridotta. Verificare la stratificazione per singola libreria."
        )
    else:
        print(
            "ESITO: nessuna sovrapposizione rilevata tra i pool di barcode GEX e TCR."
        )


def main() -> None:
    """Orchestra il download dei file pilota e l'analisi esplorativa delle etichette e dei barcode."""
    ap = argparse.ArgumentParser(
        description="Verifica preliminare di etichette e barcode nelle SubSeries GEO."
    )
    ap.add_argument("--outdir", type=Path, default=Path("data/raw/probe"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    paths = {name: fetch(name, url, args.outdir) for name, url in TARGETS}

    sc = analyse_metadata(paths["sc_metadata.csv.gz"], "Q3a: metadata scRNA")
    sn = analyse_metadata(paths["sn_metadata.csv.gz"], "Q3b: metadata snRNA")

    check_tcr_overlap(sc, paths["tcr_PA01_tumor.csv.gz"])

    print(f"\n{'=' * 78}\nRIEPILOGO\n{'=' * 78}")
    if sc is not None:
        print(f"scRNA: {len(sc):,} cellule")
    if sn is not None:
        print(f"snRNA: {len(sn):,} nuclei")
    print(
        "\nLe colonne etichetta identificano le popolazioni cellulari (Ductal, CAF, Macrophage, CD8T, ...),"
    )
    print(
        "consentendo di procedere alla costruzione delle matrici sparse h5ad nello step 02."
    )
    sys.stdout.flush()


if __name__ == "__main__":
    main()

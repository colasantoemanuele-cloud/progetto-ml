#!/usr/bin/env python3
"""
02_build_h5ad.py - Conversione delle matrici MatrixMarket di GEO in AnnData (.h5ad) a memoria controllata.

Scopo
-----
`scipy.io.mmread` su una matrice `.mtx.gz` con ~360 milioni di elementi non nulli
richiede decine di GB. Lo script legge il file compresso a blocchi
(`CHUNK = 20_000_000` righe) con il parser C di pandas, riempie tre vettori
pre-allocati (`int32` per gli indici, `float32` per i conteggi UMI) e costruisce una
matrice CSR cellule x geni; il picco di memoria misurato e' di circa 8 GB con
`--no-counts-layer` (7,7 GB per sc, `revisione/log/_status.txt`), di piu' con il layer `counts`.

Identificativi delle cellule nelle serie GEO:
    scRNA-seq        PDAC01-Tumor-GEX_out_AAACCTGAGATAGTCA-1      (barcode 16 nt)
                     PDAC05_PBMC_GEX_out_AAACCTGAGATAGTCA-1       (separatore "_")
    snRNA-seq (4/8)  1869485_FFPE_AAACAAGCAAGCTATG-1              (barcode 16 nt)
    snRNA-seq (4/8)  2021124_FFPE_AAACAAGCAATCGTGAATGTTGAC-1      (16 nt + 8 nt di sonda 10x Flex)
Libreria, barcode e sonda sono estratti con `CELL_ID_RE`. Vengono calcolati
`n_counts`, `n_genes` e `pct_mt` per cellula. h5py non accetta NaN nelle colonne di
stringhe: le colonne testuali diventano stringhe senza mancanti e poi `category`.
Il file e' scritto su un temporaneo e poi rinominato (scrittura atomica).

Input
-----
  {outdir}/raw/GSE278688/ (sc) o {outdir}/raw/GSE278689/ (sn): matrix.mtx.gz,
  barcodes.tsv.gz, features.tsv.gz, metadata.csv.gz di GEO; scaricati se mancano.

Output
------
  {outdir}/interim/{sc,sn}_raw.h5ad   ({sc,sn}_raw_dev.h5ad con --subsample)

Argomenti
---------
  --which {sc,sn} (obbligatorio), --outdir (default data),
  --subsample N (cellule casuali per sviluppo, 0 = tutte), --no-counts-layer

Uso
---
    python src/02_build_h5ad.py --which sc --outdir data --no-counts-layer
    python src/02_build_h5ad.py --which sn --outdir data --no-counts-layer
"""

from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn"

# Specifiche strutturali e cardinalita attese (numero di librerie e di pazienti) per le due SubSeries.
SPEC = {
    "sc": {
        "gse": "GSE278688",
        "n_libraries": 29,
        "n_patients": 14,
        "files": {
            "matrix": "GSE278688_sc_matrix.mtx.gz",
            "barcodes": "GSE278688_sc_barcodes.tsv.gz",
            "features": "GSE278688_sc_features.tsv.gz",
            "metadata": "GSE278688_sc_metadata.csv.gz",
        },
    },
    "sn": {
        "gse": "GSE278689",
        "n_libraries": 8,
        "n_patients": 8,
        "files": {
            "matrix": "GSE278689_sn_matrix.mtx.gz",
            "barcodes": "GSE278689_sn_barcodes.tsv.gz",
            "features": "GSE278689_sn_features.tsv.gz",
            "metadata": "GSE278689_sn_metadata.csv.gz",
        },
    },
}

CHUNK = (
    20_000_000  # Numero di righe non nulle lette per iterazione dal file Matrix Market.
)

# Espressione regolare per separare il nome della libreria, il barcode 16-nt e l'eventuale probe 8-nt.
CELL_ID_RE = r"^(?P<library>.+)_(?P<barcode>[ACGT]{16})(?P<probe>[ACGT]{8})?(?:-\d+)?$"

# Colonne di metadati da armonizzare esplicitamente al tipo categoriale per compatibilita con HDF5/h5py.
STR_COLS = ("tissue", "patients", "all_celltype", "library", "barcode", "probe")


def fetch(url: str, dest: Path) -> Path:
    """Scarica un file a blocchi di 1 MB con riutilizzo della cache locale se gia presente."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[ok] cache {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[..] scarico {dest.name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "pdac-ml/0.1"})
    with urllib.request.urlopen(req, timeout=1800) as r, open(dest, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    print(f"[ok] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def read_mtx_header(path: Path) -> tuple[int, int, int, int]:
    """Legge l'intestazione del file Matrix Market restituendo (n_geni, n_cellule, nnz, righe_header)."""
    n_header = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            n_header += 1
            if line.startswith("%"):
                continue
            r, c, nnz = (int(x) for x in line.split())
            return r, c, nnz, n_header
    raise ValueError("header mtx non trovato")


def load_matrix(
    path: Path, n_genes: int, n_cells: int, nnz: int, skip: int
) -> sp.csr_matrix:
    """Costruisce la matrice sparsa CSR (cellule x geni, float32) leggendo il file .mtx.gz a blocchi."""
    print(f"[..] mtx: {n_genes:,} geni x {n_cells:,} cellule, {nnz:,} nonzero")
    # Pre-allocazione contigua dei tre array COO per evitare riallocazioni dinamiche durante lo streaming.
    rows = np.empty(nnz, dtype=np.int32)
    cols = np.empty(nnz, dtype=np.int32)
    vals = np.empty(nnz, dtype=np.float32)

    pos = 0
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        skiprows=skip,
        header=None,
        names=["g", "c", "v"],
        dtype={"g": np.int32, "c": np.int32, "v": np.float32},
        chunksize=CHUNK,
        compression="gzip",
        engine="c",
    )
    for i, chunk in enumerate(reader):
        n = len(chunk)
        rows[pos : pos + n] = (
            chunk["g"].values - 1
        )  # Conversione da indicizzazione 1-based (MTX) a 0-based.
        cols[pos : pos + n] = chunk["c"].values - 1
        vals[pos : pos + n] = chunk["v"].values
        pos += n
        print(f"     chunk {i + 1}: {pos:,}/{nnz:,} ({pos / nnz:.0%})", flush=True)

    if pos != nnz:
        print(
            f"[!!] letti {pos:,} nonzero, header ne dichiara {nnz:,}. Continuo con {pos:,}."
        )
        rows, cols, vals = rows[:pos], cols[:pos], vals[:pos]

    # Trasposizione implicita durante la creazione della COO: righe = cellule (cols), colonne = geni (rows).
    print("[..] costruisco CSR (cellule x geni) ...")
    X = sp.coo_matrix(
        (vals, (cols, rows)), shape=(n_cells, n_genes), dtype=np.float32
    ).tocsr()
    del rows, cols, vals
    X.sort_indices()
    print(f"[ok] X {X.shape}, {X.nnz:,} nnz, {X.data.nbytes / 1e9:.2f} GB (data)")
    return X


def parse_cell_ids(idx: pd.Index) -> pd.DataFrame:
    """Scompone ciascun identificatore cellulare in library, barcode (16 nt) e probe (8 nt se presente)."""
    parts = idx.to_series().astype(str).str.extract(CELL_ID_RE)
    parts.index = idx
    n_bad = int(parts["library"].isna().sum())
    if n_bad:
        print(
            f"[!!] {n_bad:,} ID non parsati. Esempi: {idx[parts['library'].isna()][:5].tolist()}"
        )
    return parts


def harmonize_str_cols(obs: pd.DataFrame) -> pd.DataFrame:
    """Sostituisce i valori mancanti con stringa vuota e converte le colonne testuali al tipo category."""
    for c in STR_COLS:
        if c in obs.columns:
            obs[c] = (
                obs[c]
                .astype("object")
                .where(obs[c].notna(), "")
                .astype(str)
                .astype("category")
            )
    return obs


def check_counts(X: sp.csr_matrix) -> None:
    """Verifica che gli elementi non nulli della matrice rappresentino conte intere UMI non pre-normalizzate."""
    d = X[:2000].data
    is_int = bool(np.allclose(d, np.round(d)))
    print(
        f"[qc] X sembra {'CONTE INTERE' if is_int else 'GIA NORMALIZZATA (!)'}"
        f"  min={d.min():.3g} max={d.max():.3g} media={d.mean():.3g}"
    )
    if not is_int:
        print(
            "[!!] Se X e normalizzata, i metodi basati su conte grezze richiedono adattamento."
        )


def main() -> None:
    """Esegue la pipeline completa di download, assemblaggio CSR, calcolo metriche QC e scrittura atomica .h5ad."""
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Conversione delle matrici MatrixMarket di GEO in AnnData (.h5ad).",
    )
    ap.add_argument(
        "--which",
        choices=["sc", "sn"],
        required=True,
        help="scRNA-seq (sc) o snRNA-seq (sn)",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("data"),
        help="radice dei dati: scrive in <outdir>/interim, scarica in <outdir>/raw",
    )
    ap.add_argument(
        "--subsample", type=int, default=0, help="n cellule per sviluppo (0 = tutte)"
    )
    ap.add_argument(
        "--no-counts-layer",
        action="store_true",
        help="non duplicare X in layers['counts'] (dimezza RAM e disco)",
    )
    args = ap.parse_args()

    spec = SPEC[args.which]
    raw = args.outdir / "raw" / spec["gse"]
    paths = {
        k: fetch(f"{FTP}/{spec['gse']}/suppl/{fname}", raw / fname)
        for k, fname in spec["files"].items()
    }

    n_genes, n_cells, nnz, skip = read_mtx_header(paths["matrix"])

    barcodes = pd.read_csv(paths["barcodes"], header=None, compression="gzip")[
        0
    ].astype(str)
    features = pd.read_csv(paths["features"], header=None, sep="\t", compression="gzip")
    print(f"[ok] barcodes: {len(barcodes):,} | features: {features.shape}")

    if len(barcodes) != n_cells:
        sys.exit(f"[ERRORE] barcodes ({len(barcodes)}) != colonne mtx ({n_cells})")
    if len(features) != n_genes:
        sys.exit(f"[ERRORE] features ({len(features)}) != righe mtx ({n_genes})")

    X = load_matrix(paths["matrix"], n_genes, n_cells, nnz, skip)
    check_counts(X)

    var = pd.DataFrame(index=pd.Index(features.iloc[:, 0].astype(str), name="gene"))
    if features.shape[1] > 1:
        var["gene_symbol"] = features.iloc[:, 1].astype(str).values

    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_id"))
    meta = pd.read_csv(paths["metadata"], index_col=0, low_memory=False)
    print(f"[..] celle senza metadata: {len(obs.index.difference(meta.index)):,}")
    obs = obs.join(meta, how="left").join(parse_cell_ids(obs.index))
    obs = harmonize_str_cols(obs)

    n_lib = int(obs.library.nunique())
    print(f"\n[ok] librerie: {n_lib} | pazienti: {obs.patients.nunique()}")
    print(obs.library.value_counts().to_string())

    has_probe = (obs.probe.astype(str) != "").any()
    if has_probe:
        sub = obs[obs.probe.astype(str) != ""].copy()
        sub["probe"] = sub.probe.astype(str)
        print("\n[ok] probe barcode (10x Flex), atteso nunique=1 per libreria:")
        print(
            sub.groupby("library", observed=True)
            .probe.agg(["nunique", "first"])
            .to_string()
        )
    else:
        print("\n[ok] nessun probe barcode (chimica non multiplexata)")

    # Controlli di integrita strutturale rispetto al disegno sperimentale noto.
    assert n_lib == spec["n_libraries"], (
        f"{n_lib} librerie, attese {spec['n_libraries']}: il parsing degli ID e rotto."
    )
    assert obs.patients.nunique() == spec["n_patients"], (
        f"{obs.patients.nunique()} pazienti, attesi {spec['n_patients']}."
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    if not args.no_counts_layer:
        adata.layers["counts"] = adata.X.copy()

    # Calcolo vettorizzato delle metriche di controllo di qualita per singola cellula (n_counts, n_genes, pct_mt).
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    sym = (
        adata.var["gene_symbol"]
        if "gene_symbol" in adata.var
        else adata.var_names.to_series()
    )
    mt = sym.astype(str).str.upper().str.startswith("MT-").values
    adata.obs["pct_mt"] = (
        np.asarray(adata[:, mt].X.sum(axis=1)).ravel()
        / np.maximum(adata.obs.n_counts, 1)
        * 100
    )
    print(
        f"[ok] geni mitocondriali: {mt.sum()} | mediana pct_mt: {adata.obs.pct_mt.median():.2f}"
    )

    if args.subsample and args.subsample < adata.n_obs:
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(adata.n_obs, args.subsample, replace=False))
        adata = adata[keep].copy()
        print(f"[ok] subsample -> {adata.n_obs:,} cellule")

    out = (
        args.outdir
        / "interim"
        / f"{args.which}_raw{'_dev' if args.subsample else ''}.h5ad"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    # Scrittura su file temporaneo seguita da rinomina POSIX atomica per prevenire archivi HDF5 corrotti.
    tmp = out.with_suffix(".h5ad.tmp")
    adata.write_h5ad(tmp, compression="gzip")
    tmp.replace(out)
    print(f"\n[ok] scritto {out} ({out.stat().st_size / 1e6:.0f} MB)")
    print(adata)


if __name__ == "__main__":
    main()

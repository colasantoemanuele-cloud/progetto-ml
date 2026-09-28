#!/usr/bin/env python3
"""
06_clonotype_consistency.py - Analisi della coerenza clonale scTCR-seq attraverso i compartimenti.

Questo modulo utilizza la ricombinazione somatica V(D)J del recettore delle
cellule T (sequenza nucleotidica CDR3 delle catene TRA e TRB) come invariante
biologico ortogonale al trascrittoma per valutare la stabilita' delle etichette
fenotipiche (CD4T, CD8T, NK) assegnate nei diversi compartimenti anatomici
(PBMC, Adjacent_normal, Tumor) di ciascun paziente.

L'analisi articola due metriche principali:
  1. Eccesso di discordanza inter-compartimento: differenza tra la frazione di
     cloni la cui etichetta di maggioranza varia tra due compartimenti anatomici
     distinti e il tasso basale di discordanza intra-compartimento (stimato
     mediante bipartizione casuale ripetuta delle cellule appartenenti al
     medesimo clone all'interno dello stesso tessuto).
  2. Stima preliminare del limite superiore di contaminazione T all'interno
     della popolazione annotata come NK: poiche' le cellule NK canoniche non
     esprimono un recettore TCR riarrangiato, l'eccesso di sequenze CDR3 nelle
     cellule NK rispetto al rumore di fondo non-T fornisce una misura della
     frazione di linfociti T citotossici classificati come NK.

I cloni sono definiti in modo univoco all'interno di ciascun paziente tramite
la sequenza nucleotidica esatta del segmento CDR3 (evitando gli identificativi
`raw_clonotype_id` di Cell Ranger, che sono indicizzati localmente per singolo
campione).

Uso:
    bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr
    python src/06_clonotype_consistency.py \
        --h5ad data/interim/sc_raw.h5ad \
        --filelist data/raw/metadata/filelist.txt \
        --outdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import gzip
import re
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIG = Path("results/figures")
TAB = Path("results/tables")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples"
HOST = "ftp.ncbi.nlm.nih.gov"
SEED = 0
N_EXPECTED = 30

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}

T_LABELS = {"CD4T", "CD8T"}
LYMPHOID = {"CD4T", "CD8T", "NK"}
NON_T = {
    "Ductal",
    "Endothelial",
    "CAF",
    "Stellate",
    "Myeloid",
    "B",
    "Acinar",
    "Schwann",
}

FNAME_RE = re.compile(
    r"(GSM\d+)_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz"
)


def check_dns() -> None:
    """Verifica la risoluzione DNS dell'endpoint NCBI prima di tentare connessioni HTTP."""
    try:
        socket.gethostbyname(HOST)
    except socket.gaierror:
        sys.exit(
            f"[ERRORE] DNS non risolve {HOST}.\n"
            "  Diagnostica:  getent hosts ftp.ncbi.nlm.nih.gov\n"
            "  Rimedio:      sudo systemctl restart systemd-resolved\n"
            "  Oppure usa --skip-download dopo aver scaricato con download_tcr.sh"
        )


def gsm_stub(gsm: str) -> str:
    """Costruisce il prefisso di directory NCBI GEO (es. GSM8553nnn)."""
    return gsm[:-3] + "nnn"


def targets_from_filelist(filelist: Path) -> list[tuple[str, str, str, str]]:
    """Estrae dal manifesto GEO l'elenco completo dei 30 file VDJ richiesti."""
    fl = pd.read_csv(filelist, sep="\t")
    namecol = next(c for c in fl.columns if "name" in c.lower())
    out = []
    for n in fl[namecol].astype(str):
        m = FNAME_RE.match(n)
        if not m:
            continue
        gsm, pnum, tis = m.groups()
        if tis not in TISSUE_MAP:
            print(f"[!] tessuto sconosciuto in {n}: {tis}")
            continue
        out.append((gsm, f"PA{int(pnum):02d}", TISSUE_MAP[tis], n))
    if len(out) != N_EXPECTED:
        sys.exit(f"[ERRORE] {len(out)} file VDJ nel filelist, attesi {N_EXPECTED}")
    return out


def valid_gz(p: Path) -> bool:
    """Controlla che il file esista, non sia vuoto e superi la verifica dell'header gzip."""
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        with gzip.open(p, "rb") as fh:
            fh.read(1024)
        return True
    except OSError:
        return False


def fetch(url: str, dest: Path, retries: int = 4) -> bool:
    """Scarica un file VDJ da NCBI GEO con backoff esponenziale e validazione gzip."""
    if valid_gz(dest):
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "pdac-ml/0.1"})
            with urllib.request.urlopen(req, timeout=120) as r:
                dest.write_bytes(r.read())
            if valid_gz(dest):
                return True
            print(f"[!] {dest.name}: gzip corrotto, riprovo")
            dest.unlink(missing_ok=True)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            wait = 2**k
            print(f"[!] {dest.name}: {e} - riprovo tra {wait}s ({k + 1}/{retries})")
            dest.unlink(missing_ok=True)
            time.sleep(wait)
    return False


def load_contigs(paths: list[tuple[str, str, Path]]) -> pd.DataFrame:
    """Carica i contig VDJ produttivi e seleziona per ogni (cellula, catena) quello a maggior numero di UMI."""
    frames = []
    for patient, tissue, p in paths:
        df = pd.read_csv(p, compression="gzip", low_memory=False)
        # Filtro di qualita' Cell Ranger VDJ: contig cellulari, ad alta confidenza, full-length e produttivi
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df["patient"] = patient
        df["tissue"] = tissue
        # Estrazione del barcode cellulare a 16 nt per il join con l'oggetto AnnData GEX
        df["barcode16"] = df.barcode.astype(str).str.extract(
            r"([ACGT]{16})", expand=False
        )
        frames.append(
            df[["patient", "tissue", "barcode16", "chain", "cdr3_nt", "umis"]]
        )
    c = pd.concat(frames, ignore_index=True).dropna(subset=["barcode16"])
    return c.sort_values("umis", ascending=False).drop_duplicates(
        ["patient", "tissue", "barcode16", "chain"]
    )


def build_clones(contigs: pd.DataFrame) -> pd.DataFrame:
    """Riorganizza i contig per cellula definendo i clonotipi tramite CDR3(TRB) e coppia CDR3(TRA|TRB)."""
    w = contigs.pivot_table(
        index=["patient", "tissue", "barcode16"],
        columns="chain",
        values="cdr3_nt",
        aggfunc="first",
    ).reset_index()
    for ch in ("TRA", "TRB"):
        if ch not in w.columns:
            w[ch] = np.nan
    w["clone_b"] = w.TRB
    w["clone_ab"] = np.where(w.TRA.notna() & w.TRB.notna(), w.TRA + "|" + w.TRB, np.nan)
    return w


def majority(labels: pd.Series) -> str:
    """Restituisce l'etichetta modale per un insieme di cellule."""
    return labels.value_counts().index[0]


def cross_compartment_switching(
    df: pd.DataFrame, clone_col: str, min_cells: int = 2
) -> pd.DataFrame:
    """Valuta la concordanza dell'etichetta di maggioranza per i cloni presenti in almeno due compartimenti."""
    rows = []
    for (pat, clone), g in df.groupby(["patient", clone_col], observed=True):
        comp = {
            t: s
            for t, s in g.groupby("tissue", observed=True).all_celltype
            if len(s) >= min_cells
        }
        if len(comp) < 2:
            continue
        ts = sorted(comp)
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                a, b = ts[i], ts[j]
                rows.append(
                    {
                        "patient": pat,
                        "clone": clone,
                        "comp_a": a,
                        "comp_b": b,
                        "n_a": len(comp[a]),
                        "n_b": len(comp[b]),
                        "lab_a": majority(comp[a]),
                        "lab_b": majority(comp[b]),
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["switch"] = out.lab_a != out.lab_b
    return out


def within_compartment_baseline(
    df: pd.DataFrame, clone_col: str, min_cells: int = 4, n_rep: int = 20
) -> pd.DataFrame:
    """Stima il tasso basale di discordanza intra-compartimento tramite split-half casuale dei cloni."""
    rng = np.random.default_rng(SEED)
    rows = []
    for (pat, clone, tis), g in df.groupby(
        ["patient", clone_col, "tissue"], observed=True
    ):
        if len(g) < min_cells:
            continue
        lab = g.all_celltype.values
        for _ in range(n_rep):
            perm = rng.permutation(len(lab))
            h = len(lab) // 2
            rows.append(
                {
                    "patient": pat,
                    "clone": clone,
                    "tissue": tis,
                    "lab_a": pd.Series(lab[perm[:h]]).value_counts().index[0],
                    "lab_b": pd.Series(lab[perm[h : 2 * h]]).value_counts().index[0],
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["switch"] = out.lab_a != out.lab_b
    return out


def nk_contamination(gex: pd.DataFrame) -> pd.DataFrame:
    """Calcola lo stimatore lineare c_hat = (f - q) / (pi - q) della frazione di cellule T nella classe NK."""
    rows = []
    for pat, g in gex.groupby("patients", observed=True):
        t = g[g.all_celltype.isin(T_LABELS)]
        nt = g[g.all_celltype.isin(NON_T)]
        nk = g[g.all_celltype == "NK"]
        if min(len(t), len(nt), len(nk)) < 30:
            continue
        pi, q, f = t.has_tcr.mean(), nt.has_tcr.mean(), nk.has_tcr.mean()
        c = (f - q) / (pi - q) if pi > q else np.nan
        rows.append(
            {
                "patient": pat,
                "n_NK": len(nk),
                "pi": pi,
                "q": q,
                "f": f,
                "c_hat": np.clip(c, 0, 1),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    """Esegue l'intera pipeline di verifica dei clonotipi TCR e genera tabelle e figure."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--filelist", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=Path("data/raw/tcr"))
    ap.add_argument("--clone-def", choices=["clone_b", "clone_ab"], default="clone_b")
    ap.add_argument(
        "--skip-download", action="store_true", help="usa solo i file gia' su disco"
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    tg = targets_from_filelist(args.filelist)
    print(f"[ok] {len(tg)} file VDJ nel filelist")

    if not args.skip_download and not all(valid_gz(args.outdir / f) for *_, f in tg):
        check_dns()

    # Verifica di completezza dei 30 file VDJ (arresto immediato in caso di archivi mancanti)
    paths, missing = [], []
    for gsm, pat, tis, fname in tg:
        dest = args.outdir / fname
        if valid_gz(dest):
            paths.append((pat, tis, dest))
        elif args.skip_download or not fetch(
            f"{FTP}/{gsm_stub(gsm)}/{gsm}/suppl/{fname}", dest
        ):
            missing.append(fname)
        else:
            paths.append((pat, tis, dest))

    print(f"[ok] {len(paths)}/{len(tg)} file validi in {args.outdir}")
    if missing:
        sys.exit(
            f"[ERRORE] {len(missing)} file mancanti. Esecuzione interrotta.\n"
            "  L'analisi richiede la copertura integrale dei 30 campioni VDJ.\n"
            f"  Mancanti (primi 5): {missing[:5]}\n"
            "  Rimedio: bash src/download_tcr.sh, poi rilancia con --skip-download"
        )

    contigs = load_contigs(paths)
    cells = build_clones(contigs)
    print(f"\n[ok] {len(cells):,} cellule con almeno un contig produttivo")
    print(
        f"     con TRB: {cells.clone_b.notna().sum():,} | con TRA+TRB: {cells.clone_ab.notna().sum():,}"
    )

    # Apertura in sola lettura dei metadati GEX e allineamento tramite chiave composita (paziente, tessuto, barcode16)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype", "barcode"):
        obs[c] = obs[c].astype(str)
    obs = obs.rename(columns={"barcode": "barcode16"})

    gex = obs.merge(
        cells[["patient", "tissue", "barcode16", "clone_b", "clone_ab"]],
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    print(
        f"\n[ok] cellule GEX: {len(gex):,} | con TCR appaiato: {gex.has_tcr.sum():,} "
        f"({gex.has_tcr.mean():.1%})"
    )
    print("\n[qc] tasso di TCR per etichetta (riferimento ortogonale VDJ):")
    print(
        gex.groupby("all_celltype", observed=True)
        .has_tcr.agg(["mean", "sum", "size"])
        .sort_values("mean", ascending=False)
        .round(4)
        .to_string()
    )

    nk = nk_contamination(gex)
    nk.to_csv(TAB / "nk_contamination.csv", index=False)
    print("\n" + "=" * 74)
    print("CONTAMINAZIONE NK  (limite SUPERIORE: le NKT reali hanno un TCR)")
    print("=" * 74)
    print(nk.round(4).to_string(index=False))
    if not nk.empty:
        print(
            f"\nc_hat mediano su {len(nk)} pazienti: {nk.c_hat.median():.3f} "
            f"[IQR {nk.c_hat.quantile(0.25):.3f}-{nk.c_hat.quantile(0.75):.3f}]"
        )

    cc = args.clone_def
    lym = gex[gex[cc].notna() & gex.all_celltype.isin(LYMPHOID)].copy()

    print("\n" + "=" * 74)
    print(f"CLONI CONDIVISI FRA COMPARTIMENTI  (definizione: {cc})")
    print("=" * 74)
    spanning = (
        lym.groupby(["patient", cc], observed=True)
        .agg(n_cells=("all_celltype", "size"), n_comp=("tissue", "nunique"))
        .reset_index()
    )
    print(f"cloni totali (linfoidi)   : {len(spanning):,}")
    print(f"cloni con >=2 cellule     : {(spanning.n_cells >= 2).sum():,}")
    print(f"cloni in >=2 compartimenti: {(spanning.n_comp >= 2).sum():,}")
    print(
        f"cellule in quei cloni     : {spanning.loc[spanning.n_comp >= 2, 'n_cells'].sum():,}"
    )
    print("\nper paziente:")
    print(
        spanning.assign(multi=spanning.n_comp >= 2)
        .groupby("patient", observed=True)
        .agg(cloni=("multi", "size"), multi_comp=("multi", "sum"))
        .to_string()
    )

    across = cross_compartment_switching(lym, cc)
    within = within_compartment_baseline(lym, cc)

    print("\n" + "=" * 74)
    print("SINTESI DISCORDANZA CLONALE")
    print("=" * 74)
    if across.empty or within.empty:
        print("Numero insufficiente di cloni condivisi per la stima.")
        return

    a_rate, w_rate = across.switch.mean(), within.switch.mean()
    print(f"cloni-coppie fra compartimenti : {len(across):,}")
    print(f"discordanza FRA compartimenti  : {a_rate:.3f}")
    print(
        f"discordanza ENTRO compartimento: {w_rate:.3f}   (controllo, {len(within):,} repliche)"
    )
    print(f"ECCESSO                        : {a_rate - w_rate:+.3f}")

    # Bootstrap non parametrico sulle coppie di compartimenti
    rng = np.random.default_rng(SEED)
    boot = [
        across.switch.sample(len(across), replace=True, random_state=int(s)).mean()
        - within.switch.sample(len(within), replace=True, random_state=int(s)).mean()
        for s in rng.integers(0, 10**6, 1000)
    ]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"IC 95% bootstrap               : [{lo:+.3f}, {hi:+.3f}]")

    print("\nper coppia di compartimenti:")
    print(
        across.groupby(["comp_a", "comp_b"])
        .switch.agg(["mean", "size"])
        .round(3)
        .to_string()
    )
    print("\ntransizioni di etichetta piu' frequenti:")
    print(
        across[across.switch]
        .groupby(["lab_a", "lab_b"])
        .size()
        .sort_values(ascending=False)
        .head(10)
        .to_string()
    )

    across.to_csv(TAB / "clone_switching_across.csv", index=False)
    within.to_csv(TAB / "clone_switching_within.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(
        ["entro\ncompartimento\n(rumore)", "fra\ncompartimenti"],
        [w_rate, a_rate],
        color=["#7570b3", "#d95f02"],
        alpha=0.85,
    )
    ax.axhline(w_rate, ls="--", c="k", lw=0.8)
    ax.set_ylabel("frazione di cloni con etichetta discordante")
    ax.set_title(f"eccesso = {a_rate - w_rate:+.3f}   IC95% [{lo:+.3f}, {hi:+.3f}]")
    fig.tight_layout()
    fig.savefig(FIG / "clone_label_switching.png", dpi=140)
    plt.close(fig)

    print("\n" + "=" * 74)
    if lo > 0:
        print("ESITO STATISTICO: eccesso significativamente positivo.")
        print(
            "  L'etichetta fenotipica varia sistematicamente con lo stato compartimentale."
        )
    else:
        print("ESITO STATISTICO: eccesso non significativamente diverso da zero.")
        print("  L'annotazione risulta statisticamente coerente all'interno dei cloni.")
    print("=" * 74)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
13_infercnv.py - Score CNV (inferCNV), etichetta di malignita' e classificazione LOPO delle duttali.

Scopo
-----
L'annotazione di GSE278694 riunisce tutte le cellule epiteliali duttali nella classe
`Ductal`. Lo script costruisce un'etichetta di malignita' derivata dal profilo CNV
inferito dal trascrittoma e misura quanto un classificatore la predice dal
trascrittoma stesso, in tre stadi:

1. Score CNV con `infercnvpy`: finestre genomiche di 100 geni, UNA sola esecuzione
   con `reference_key="ref_group"`. Filtro del rumore: infercnvpy e' chiamato con
   `dynamic_threshold=None` e poi le voci con |x| < 1,5 x sd di TUTTA `X_cnv` sono
   azzerate (`global_noise_filter`). Il default di infercnvpy calcola questa soglia per
   blocco di 5.000 cellule: con le cellule in ordine di libreria i blocchi sono confusi
   con tessuto e quota di duttali, e lo score dipenderebbe dall'ordine delle righe. Il riferimento
   diploide sono le cellule immunitarie (`CD8T`, `CD4T`, `NK`, `B`, `Myeloid`,
   `Plasma`) di TUTTI i pazienti insieme, non dello stesso paziente. Solo la soglia
   e' calcolata per paziente:
       soglia(paziente) = mediana(score_immuni) + 3 * MAD_normalizzata(score_immuni)
   La soglia e' stimata sulle stesse cellule immunitarie su cui inferCNV centra lo
   score: un programma trascrizionale epiteliale qualsiasi (non solo
   un'aneuploidia) sposta lo score delle duttali sopra quello degli immuni e tende a
   superarla: in `malignancy_labels.csv` risultano "maligne" anche molte duttali di
   `Adjacent_normal` (la frazione per paziente e' nella colonna `frac_malig`).
   I pazienti con meno di 20 immunitarie non hanno soglia: le loro duttali restano
   senza etichetta (NaN) e sono escluse dalla classificazione.
   RNA ambientale acinare: i geni degli enzimi acinari sono raggruppati nel genoma
   (CTRC, CELA2A/2B in 1p36.21 e CELA3A/3B in 1p36.12, dentro il locus `1p_loss` di
   14_cnv_loci.py; PRSS1/2 in 7q34). L'RNA ambientale acinare, abbondante nel tessuto
   adiacente, puo' creare pseudo-guadagni in queste regioni nelle duttali e quindi
   influenzare lo score e l'AUC tumore contro adiacente.

   Coordinate dei geni: Ensembl release 116 (GRCh38), la stessa delle cache precedenti
   dell'11/07/2026 (ora in `data/interim/cache_cnv_20260711/`: le coordinate delle loro
   finestre coincidono con quelle dei geni di Ensembl 116); le cache attuali sono
   ricalcolate con questa tabella e con la regola dei cromosomi qui sotto. La tabella e'
   scaricata da BioMart via https solo da un host che dichiara la release 116 e
   salvata in `data/interim/gene_positions_ensembl116.csv`, poi riletta da li'.
   Se non e' disponibile lo script si ferma: `copykat` non e' un ripiego automatico
   (richiede `rpy2` e i pacchetti R `copykat` e `stringr`, assenti da
   `envs/environment.yml`) e si usa solo con `--method copykat`.
   Cromosomi: entrano in inferCNV solo i geni dei cromosomi 1-22, X e Y (`add_positions`
   scarta il genoma mitocondriale `MT` e i contig non standard come `GL000194.1`);
   infercnvpy esclude poi X e Y (`exclude_chromosomes` di default), quindi le finestre
   sono solo autosomiche, e lo script si ferma se ne trova altre. Il genoma
   mitocondriale non e' una copia del DNA nucleare e la sua espressione segue la
   frazione di letture mitocondriali (metrica di QC, `pct_mt`): includerlo farebbe
   dipendere lo score CNV dalla qualita' delle cellule.
   Cache: `cnv_scores.csv` e' riusata solo se la sua impronta, salvata in
   `data/interim/cnv_cache_fingerprint.json`, coincide con quella calcolata ora
   (hash delle cellule, dei conteggi, dei geni, della tabella delle coordinate, del
   gruppo di riferimento, della versione di infercnvpy e dei parametri di inferCNV:
   window_size, step, lfc_clip, dynamic_threshold, exclude_chromosomes, regola del filtro
   del rumore e delle coordinate delle finestre); altrimenti lo
   script si ferma. Se il sottoinsieme non ha cellule immunitarie di riferimento, o ha
   solo quelle, lo script si ferma prima di inferCNV.

2. Classificazione Leave-One-Patient-Out: Random Forest bilanciato sulle duttali del
   tessuto `Tumor` dei pazienti con almeno 150 duttali. Normalizzazione (log1p CP10k)
   e selezione dei geni piu' dispersi (`HVGSelector`, default 2.000) stanno dentro
   una `Pipeline` stimata sui soli pazienti di training di ciascun fold.

3. Applicazione alle duttali di `Adjacent_normal`: un modello finale stimato su tutte
   le duttali tumorali di training (geni scelti solo su di esse) viene applicato alle
   duttali adiacenti, mai viste in training. Non e' una validazione indipendente da
   inferCNV: le duttali adiacenti sono nella stessa esecuzione di inferCNV, con lo
   stesso riferimento e la stessa regola di soglia che definisce l'etichetta.
   I file `external_test_*` hanno un nome storico: NON sono una validazione esterna.

Input
-----
  --h5ad       data/interim/sc_raw.h5ad (conteggi grezzi; obs: tissue, patients, all_celltype)
  --cnv-cache  data/interim/cnv_scores.csv (letto se esiste, contiene esattamente le cellule
               selezionate e ha l'impronta attesa; altrimenti calcolato e scritto se manca,
               mentre una cache diversa ferma lo script)
  data/interim/cnv_cache_fingerprint.json (impronte delle cache di 13 e 14)
  data/interim/gene_positions_ensembl116.csv (letta se esiste, altrimenti scaricata da BioMart)

Output
------
  data/interim/cnv_scores.csv                   score CNV per cellula (cache)
  data/interim/cnv_cache_fingerprint.json       impronta della cache (voce `cnv_scores.csv`)
  data/interim/gene_positions_ensembl116.csv    coordinate dei geni (solo se scaricate)
  results/tables/malignancy_labels.csv          duttali maligne per paziente x tessuto
  results/tables/malignancy_lopo.csv            metriche LOPO per paziente (BA, F1, AUC)
  results/tables/malignancy_lopo_summary.csv    media e deviazione standard delle metriche LOPO
  results/tables/external_test_predictions.csv  predizione per ogni duttale adiacente
  results/tables/external_test_summary.csv      frazione di duttali adiacenti predette maligne
  results/figures/malignancy.png
  (se nessun fold LOPO e' valutabile, `malignancy_lopo.csv` resta vuoto e
  `malignancy_lopo_summary.csv` non viene scritto, e quello di un'esecuzione precedente
  viene cancellato; se le duttali adiacenti sono 20 o meno, i file `external_test_*`
  non vengono scritti e quelli precedenti vengono cancellati)

Argomenti
---------
  --h5ad (obbligatorio), --method {infercnv,copykat,auto} (default auto = infercnv,
  senza ripiego automatico su copykat), --cnv-cache, --n-hvg (default 2000)

Uso
---
    python src/13_infercnv.py --h5ad data/interim/sc_raw.h5ad
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import re
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import HVGSelector, LogCPM  # noqa: E402

from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0

IMMUNE_REF = ["CD8T", "CD4T", "NK", "B", "Myeloid", "Plasma"]
# Supporto minimo di cellule duttali tumorali per l'inclusione del paziente nella validazione LOPO
MIN_DUCTAL = 150
N_HVG = 2000

# Coordinate dei geni: release fissa di Ensembl (vedi docstring, punto 1)
ENSEMBL_RELEASE = 116
BIOMART_HOSTS = ("https://jun2026.archive.ensembl.org", "https://www.ensembl.org")
GENE_POS = Path(f"data/interim/gene_positions_ensembl{ENSEMBL_RELEASE}.csv")
BIOMART_QUERY = (
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
    '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" '
    'completionStamp="1" datasetConfigVersion="0.6">'
    '<Dataset name="hsapiens_gene_ensembl" interface="default">'
    '<Attribute name="hgnc_symbol"/><Attribute name="start_position"/>'
    '<Attribute name="end_position"/><Attribute name="chromosome_name"/>'
    "</Dataset></Query>"
)
POS_COLS = ["hgnc_symbol", "start_position", "end_position", "chromosome_name"]
# cromosomi annotati: 1-22, X, Y (esclusi MT e i contig non standard; vedi docstring)
STANDARD_CHROM = re.compile(r"^(\d{1,2}|X|Y)$")
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
# parametri di inferCNV, comuni a 13 e 14 e registrati nell'impronta delle cache
WINDOW_SIZE = 100
# filtro del rumore: una sola soglia globale su tutta X_cnv (vedi `global_noise_filter`)
NOISE_K = 1.5
NOISE_RULE = (
    "dynamic_threshold=None in infercnvpy; poi |x| < 1.5 * sd(tutta X_cnv) -> 0 "
    "(soglia globale, indipendente dall'ordine delle cellule)"
)
WINDOW_RULE = (
    "finestra i di un cromosoma = geni ordinati per start [i*step, i*step+window_size); "
    "coordinate = start del primo gene, end massimo, mediana degli start"
)
FINGERPRINT = Path("data/interim/cnv_cache_fingerprint.json")
CHROM_RULE = "1-22, X, Y in add_positions; X e Y esclusi da infercnvpy"


def gene_positions(path: Path = GENE_POS) -> pd.DataFrame | None:
    """Coordinate dei geni della release ENSEMBL_RELEASE (colonne POS_COLS).

    Legge `path` se esiste; altrimenti interroga BioMart via https sul primo host di
    BIOMART_HOSTS il cui registro dichiara esattamente quella release, controlla che
    la risposta sia completa (`[success]` finale) e la salva in `path`. Restituisce
    None se nessun host risponde con la release giusta.
    """
    if path.exists():
        return pd.read_csv(path, dtype={"hgnc_symbol": str, "chromosome_name": str})
    import requests

    for host in BIOMART_HOSTS:
        url = f"{host}/biomart/martservice"
        try:
            reg = requests.get(url, params={"type": "registry"}, timeout=60)
            if f'displayName="Ensembl Genes {ENSEMBL_RELEASE}"' not in reg.text:
                print(f"[!!] {host}: release {ENSEMBL_RELEASE} non servita")
                continue
            r = requests.get(url, params={"query": BIOMART_QUERY}, timeout=900)
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"[!!] {host}: {type(e).__name__}: {e}")
            continue
        text = r.text.rstrip("\n")
        if not text.endswith("[success]"):
            print(f"[!!] {host}: risposta incompleta o errore BioMart")
            continue
        tab = pd.read_csv(io.StringIO(text[: -len("[success]")]), sep="\t", dtype=str)
        if tab.shape[1] != len(POS_COLS):
            print(f"[!!] {host}: colonne inattese {list(tab.columns)}")
            continue
        tab.columns = POS_COLS
        for c in ("start_position", "end_position"):
            tab[c] = pd.to_numeric(tab[c])
        path.parent.mkdir(parents=True, exist_ok=True)
        tab.to_csv(path, index=False)
        print(f"[ok] coordinate Ensembl {ENSEMBL_RELEASE} da {host} -> {path}")
        return tab
    return None


def require_copykat() -> None:
    """Esce con un messaggio chiaro se mancano le dipendenze di copykat (rpy2 + R)."""
    if importlib.util.find_spec("rpy2") is None:
        raise SystemExit(
            "[ERRORE] --method copykat richiede rpy2 e i pacchetti R copykat e stringr, "
            "assenti da envs/environment.yml: usa --method infercnv."
        )


def add_positions(adata: ad.AnnData) -> bool:
    """Annota `chromosome`, `start`, `end` in `adata.var` dalla release Ensembl fissa.

    Stessa regola di `infercnvpy.io.genomic_position_from_biomart`: i simboli HGNC
    duplicati nella tabella sono scartati e il cromosoma prende il prefisso `chr`.
    Restano solo i cromosomi 1-22, X e Y (STANDARD_CHROM): i geni mitocondriali e
    dei contig non standard restano senza coordinate e non entrano in inferCNV.
    """
    adata.var_names = (
        adata.var["gene_symbol"].astype(str).values
        if "gene_symbol" in adata.var
        else adata.var_names
    )
    adata.var_names_make_unique()
    print(f"[..] coordinate dei geni (Ensembl {ENSEMBL_RELEASE}) ...")
    tab = gene_positions()
    if tab is None:
        return False
    tab = tab.dropna(subset=["hgnc_symbol"])
    tab = tab[~tab.hgnc_symbol.duplicated(keep=False)].set_index("hgnc_symbol")
    standard = tab.chromosome_name.astype(str).str.match(STANDARD_CHROM)
    n_out = int((~standard & tab.index.isin(adata.var_names)).sum())
    tab = tab[standard]
    print(f"[ok] {n_out} geni su cromosomi non standard (MT, contig) esclusi")
    adata.var["chromosome"] = (
        ("chr" + tab.chromosome_name).reindex(adata.var_names).values
    )
    adata.var["start"] = tab.start_position.reindex(adata.var_names).values
    adata.var["end"] = tab.end_position.reindex(adata.var_names).values
    ok = adata.var["chromosome"].notna().sum()
    print(f"[ok] posizioni per {ok:,}/{adata.n_vars:,} geni")
    return ok > 5000


def run_infercnv(adata: ad.AnnData) -> np.ndarray:
    """Esegue infercnvpy su dati log-normalizzati e restituisce lo score medio |X_cnv| per cellula.

    Una sola esecuzione su tutte le cellule: il riferimento sono le immunitarie di
    tutti i pazienti insieme (categorie di `ref_group` in IMMUNE_REF); soglia del rumore
    globale (`run_infercnv_core`).
    """
    a = adata.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a = a[:, a.var["chromosome"].notna()].copy()

    ref = reference_categories(a.obs["ref_group"])
    print(f"[..] infercnv: riferimento = {ref} (tutti i pazienti)")
    run_infercnv_core(a, ref, WINDOW_SIZE)
    score = np.asarray(np.abs(a.obsm["X_cnv"]).mean(axis=1)).ravel()
    return score


def global_noise_filter(X, k: float = NOISE_K):
    """Azzera le voci con |x| < k * sd, con la sd calcolata su TUTTA la matrice.

    infercnvpy (`dynamic_threshold`) calcola la stessa soglia separatamente per ogni
    blocco di `chunksize` cellule: con le cellule in ordine di libreria i blocchi sono
    confusi con tessuto e tipo cellulare, e la soglia dipende dall'ordine delle righe.
    Qui la soglia e' una sola (equivale al default di infercnvpy con un solo blocco).
    La sd e' quella di popolazione (ddof = 0), come `np.std`, e include gli zeri
    impliciti di una matrice sparsa. Restituisce (matrice filtrata, soglia).
    """
    n = X.shape[0] * X.shape[1]
    if sp.issparse(X):
        X = sp.csr_matrix(X, copy=True)
        vals = X.data
    else:
        X = np.array(X, copy=True)
        vals = X.ravel()
    s1 = vals.sum(dtype=np.float64)
    s2 = np.square(vals, dtype=np.float64).sum()
    sd = float(np.sqrt(max(s2 / n - (s1 / n) ** 2, 0.0)))
    thr = k * sd
    vals[np.abs(vals) < thr] = 0
    if sp.issparse(X):
        X.eliminate_zeros()
    return X, thr


def run_infercnv_core(a: ad.AnnData, ref: list[str], window_size: int) -> float:
    """inferCNV senza soglia per blocco, poi la soglia globale (`global_noise_filter`).

    `a` deve essere log-normalizzato e avere le coordinate dei geni; modifica
    `a.obsm["X_cnv"]` sul posto e restituisce la soglia applicata. Gli altri passi di
    infercnvpy (log fold change rispetto al riferimento, clip, media mobile per
    cromosoma, centratura per cellula) sono per cellula e non dipendono dai blocchi.
    """
    import infercnvpy as cnv

    cnv.tl.infercnv(
        a,
        reference_key="ref_group",
        reference_cat=ref,
        window_size=window_size,
        dynamic_threshold=None,
    )
    check_autosomal_windows(a.uns["cnv"]["chr_pos"])
    a.obsm["X_cnv"], thr = global_noise_filter(a.obsm["X_cnv"])
    print(f"[ok] soglia del rumore globale: {thr:.4f} (1,5 x sd di tutta X_cnv)")
    return thr


def window_table(a: ad.AnnData, window_size: int = WINDOW_SIZE) -> pd.DataFrame:
    """Coordinate genomiche reali delle finestre di `X_cnv` (stessa costruzione di infercnvpy).

    In infercnvpy 0.6.1 i geni di ciascun cromosoma (esclusi quelli senza coordinate e
    i cromosomi di `exclude_chromosomes`) sono ordinati per `start`; la media mobile
    di `window_size` geni e' calcolata a ogni posizione e se ne tiene una ogni `step`:
    la finestra i copre i geni [i*step, i*step + window_size) (tutti i geni se il
    cromosoma ne ha meno di `window_size`). Per ogni finestra: cromosoma, start del
    primo gene, end massimo dei suoi geni, mediana degli start, numero di geni, primo e
    ultimo gene. Si ferma se il numero di finestre non coincide con `X_cnv`.
    """
    import inspect

    import infercnvpy as cnv

    sig = inspect.signature(cnv.tl.infercnv).parameters
    step = sig["step"].default
    excl = set(sig["exclude_chromosomes"].default or [])
    var = a.var.loc[a.var["chromosome"].notna(), ["chromosome", "start", "end"]]
    var = var[~var["chromosome"].isin(excl)]
    chr_pos = a.uns["cnv"]["chr_pos"]
    n_win = a.obsm["X_cnv"].shape[1]
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    rows = []
    for i, (chrom, start_col) in enumerate(order):
        end_col = order[i + 1][1] if i + 1 < len(order) else n_win
        # stessa espressione di infercnvpy (_running_mean_for_chromosome)
        genes = var.loc[var["chromosome"] == chrom].sort_values("start")
        n = len(genes)
        starts = np.arange(0, n - window_size + 1, step) if n > window_size else [0]
        width = window_size if n > window_size else n
        if len(starts) != end_col - start_col:
            raise SystemExit(
                f"[ERRORE] {chrom}: {end_col - start_col} finestre in X_cnv ma "
                f"{len(starts)} ricostruite dai geni"
            )
        for j, g0 in enumerate(starts):
            w = genes.iloc[g0 : g0 + width]
            rows.append(
                {
                    "col": start_col + j,
                    "chrom": chrom.replace("chr", ""),
                    "start": float(w.start.iloc[0]),
                    "end": float(w.end.max()),
                    "pos": float(np.median(w.start.values)),
                    "n_genes": len(w),
                    "first_gene": w.index[0],
                    "last_gene": w.index[-1],
                }
            )
    return pd.DataFrame(rows).sort_values("col").reset_index(drop=True)


def reference_categories(ref_group: pd.Series) -> list[str]:
    """Categorie immunitarie di riferimento presenti; si ferma se il riferimento e'
    vuoto o se tutte le cellule sono di riferimento (nessuna cellula da valutare)."""
    present = set(ref_group.astype(str).unique())
    ref = [c for c in IMMUNE_REF if c in present]
    if not ref:
        raise SystemExit(
            "[ERRORE] nessuna cellula immunitaria di riferimento nel sottoinsieme "
            f"(attese categorie fra {IMMUNE_REF}): inferCNV non puo' stimare il profilo diploide."
        )
    if present <= set(ref):
        raise SystemExit(
            "[ERRORE] tutte le cellule sono di riferimento: nessuna duttale da valutare."
        )
    return ref


def infercnv_params() -> dict:
    """Versione di infercnvpy, parametri di `infercnv` (passati e default), filtro del
    rumore e regola delle coordinate delle finestre."""
    import inspect

    import infercnvpy as cnv

    sig = inspect.signature(cnv.tl.infercnv).parameters
    params = {k: sig[k].default for k in ("lfc_clip", "step", "exclude_chromosomes")}
    params["exclude_chromosomes"] = list(params["exclude_chromosomes"] or [])
    params["window_size"] = WINDOW_SIZE
    params["dynamic_threshold"] = None  # passato esplicitamente (run_infercnv_core)
    params["noise_filter"] = NOISE_RULE
    params["window_coordinates"] = WINDOW_RULE
    return {"infercnvpy_version": cnv.__version__, "infercnv_params": params}


def check_autosomal_windows(chr_pos: dict) -> None:
    """Si ferma se `X_cnv` contiene finestre fuori dai cromosomi 1-22."""
    other = sorted(set(chr_pos) - AUTOSOMES)
    if other:
        raise SystemExit(
            f"[ERRORE] finestre inferCNV su cromosomi non autosomici: {other} "
            "(attesi solo chr1-chr22; vedi add_positions)"
        )


def _digest(*parts: bytes) -> str:
    """Hash BLAKE2b a 16 byte, in esadecimale, dei blocchi `parts` in sequenza."""
    h = hashlib.blake2b(digest_size=16)
    for b in parts:
        h.update(b)
    return h.hexdigest()


def cnv_fingerprint(sub: ad.AnnData, method: str) -> dict:
    """Impronta di una cache CNV: cellule, conteggi, geni, coordinate, riferimento, parametri.

    Va calcolata sullo stesso sottoinsieme `sub` quando la cache viene scritta e quando
    viene riletta; i geni sono presi da `gene_symbol` (se presente), quindi l'impronta non
    dipende dalla rinomina fatta da `add_positions`.
    """
    X = sp.csr_matrix(sub.X)
    genes = (
        sub.var["gene_symbol"].astype(str)
        if "gene_symbol" in sub.var
        else sub.var_names
    )
    return {
        "method": method,
        "ensembl_release": ENSEMBL_RELEASE,
        "chromosomes": CHROM_RULE,
        "window_size": WINDOW_SIZE,
        **infercnv_params(),
        "reference": IMMUNE_REF,
        "cells": _digest("\n".join(sub.obs_names).encode()),
        "ref_group": _digest("\n".join(sub.obs["ref_group"].astype(str)).encode()),
        "counts": _digest(X.data.tobytes(), X.indices.tobytes(), X.indptr.tobytes()),
        "genes": _digest("\n".join(map(str, genes)).encode()),
        "gene_positions": _digest(GENE_POS.read_bytes()) if GENE_POS.exists() else None,
    }


def check_fingerprint(key: str, fp: dict, cache: Path) -> None:
    """Si ferma se l'impronta salvata per la cache `key` manca o non coincide con `fp`."""
    saved = (
        json.loads(FINGERPRINT.read_text()).get(key) if FINGERPRINT.exists() else None
    )
    if saved is None:
        raise SystemExit(
            f"[ERRORE] la cache {cache} non ha un'impronta in {FINGERPRINT}: "
            "cancellala per ricalcolare inferCNV."
        )
    diff = sorted(k for k in fp if saved.get(k) != fp[k])
    if diff:
        raise SystemExit(
            f"[ERRORE] la cache {cache} e' stata calcolata con input o parametri diversi "
            f"({', '.join(diff)}): cancellala per ricalcolare inferCNV."
        )


def save_fingerprint(key: str, fp: dict) -> None:
    """Registra l'impronta della cache `key` in FINGERPRINT (le altre voci restano)."""
    allfp = json.loads(FINGERPRINT.read_text()) if FINGERPRINT.exists() else {}
    allfp[key] = fp
    FINGERPRINT.parent.mkdir(parents=True, exist_ok=True)
    FINGERPRINT.write_text(json.dumps(allfp, indent=2) + "\n")


def run_copykat(adata: ad.AnnData) -> np.ndarray:
    """Esegue l'algoritmo CopyKAT sulle conte grezze e restituisce lo score medio |X_cnv| per cellula."""
    import infercnvpy as cnv

    a = adata.copy()
    print("[..] copykat (lento) ...")
    cnv.tl.copykat(a, n_jobs=-1)
    score = np.asarray(np.abs(a.obsm["X_cnv"]).mean(axis=1)).ravel()
    return score


def make_pipeline(n_hvg: int = N_HVG) -> Pipeline:
    """Pipeline log1p(CP10k) -> geni piu' dispersi -> Random Forest bilanciato.

    Normalizzazione e selezione dei geni sono stimate dentro `fit()`, quindi solo
    sulle cellule di training su cui la pipeline viene addestrata.
    """
    return Pipeline(
        [
            ("logcpm", LogCPM()),
            ("hvg", HVGSelector(n_top=n_hvg)),
            (
                "clf",
                RandomForestClassifier(
                    n_estimators=300,
                    max_features="sqrt",
                    class_weight="balanced_subsample",
                    n_jobs=-1,
                    random_state=SEED,
                ),
            ),
        ]
    )


def classify_lopo(X, y, groups, n_hvg: int = N_HVG, return_models: bool = False):
    """Cross-validation Leave-One-Patient-Out sui conteggi grezzi `X`.

    Per ogni paziente escluso, la pipeline (normalizzazione, HVG, Random Forest) e'
    stimata sui soli pazienti di training. Restituisce (metriche per paziente,
    predizioni) e, con `return_models=True`, anche il dizionario paziente -> pipeline.
    I fold il cui training ha una sola classe sono saltati; l'AUC e la balanced
    accuracy sono NaN se il paziente escluso ha una sola classe.
    """
    rows, preds, models = [], [], {}
    for g in np.unique(groups):
        tr, te = groups != g, groups == g
        if len(np.unique(y[tr])) < 2 or te.sum() == 0:
            print(f"  [!] fold {g} saltato: training con una sola classe o test vuoto")
            continue
        pipe = make_pipeline(n_hvg)
        pipe.fit(X[tr], y[tr])
        models[g] = pipe
        p = pipe.predict(X[te])
        proba = pipe.predict_proba(X[te])[:, list(pipe.classes_).index(1)]
        preds.append(
            pd.DataFrame({"patient": g, "y_true": y[te], "y_pred": p, "proba": proba})
        )
        two = len(np.unique(y[te])) > 1
        rows.append(
            {
                "patient": g,
                "n_test": int(te.sum()),
                "balanced_acc": balanced_accuracy_score(y[te], p) if two else np.nan,
                "f1_malig": f1_score(y[te], p, pos_label=1, zero_division=0),
                "auc": roc_auc_score(y[te], proba) if two else np.nan,
            }
        )
    res = pd.DataFrame(rows)
    pr = pd.concat(preds) if preds else pd.DataFrame()
    return (res, pr, models) if return_models else (res, pr)


def fit_final(X, y, n_hvg: int = N_HVG) -> Pipeline:
    """Modello finale su tutte le duttali tumorali di training (geni scelti solo su di esse)."""
    return make_pipeline(n_hvg).fit(X, y)


def main() -> None:
    """Punto di ingresso principale per l'inferenza CNV, la sogliatura per paziente e il test LOPO."""
    ap = argparse.ArgumentParser(
        description="Score CNV, etichetta di malignita' e classificazione LOPO delle duttali.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--h5ad", type=Path, required=True, help="AnnData con conteggi grezzi (da 02)"
    )
    ap.add_argument(
        "--method",
        choices=["infercnv", "copykat", "auto"],
        default="auto",
        help="auto = infercnv, senza ripiego su copykat (che richiede rpy2 e R)",
    )
    ap.add_argument(
        "--cnv-cache",
        type=Path,
        default=Path("data/interim/cnv_scores.csv"),
        help="cache dello score CNV per cellula (con impronta in cnv_cache_fingerprint.json)",
    )
    ap.add_argument(
        "--n-hvg",
        type=int,
        default=N_HVG,
        help="geni piu' dispersi nella pipeline LOPO",
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    for c in ("tissue", "patients", "all_celltype"):
        adata.obs[c] = adata.obs[c].astype(str)

    # Selezione delle cellule duttali (Tumor + Adjacent_normal) e del riferimento immunitario
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
    if (sub.obs.all_celltype == "Ductal").sum() == 0:
        raise SystemExit("[ERRORE] nessuna cellula duttale in Tumor/Adjacent_normal.")

    duct = sub.obs[sub.obs.all_celltype == "Ductal"]
    print("\nduttali per paziente x tessuto:")
    print(pd.crosstab(duct.patients, duct.tissue, margins=True).to_string())

    # Calcolo o caricamento da cache dello score CNV per cellula
    if args.cnv_cache.exists():
        print(f"\n[ok] CNV score da cache {args.cnv_cache}")
        cs = pd.read_csv(args.cnv_cache, index_col=0)
        score = cs.reindex(sub.obs_names).cnv_score
        if score.isna().any():
            raise SystemExit(
                f"[ERRORE] la cache {args.cnv_cache} non copre {int(score.isna().sum()):,} "
                f"cellule su {sub.n_obs:,}: cancellala per ricalcolare inferCNV."
            )
        # cache piu' grande del sottoinsieme: calcolata su un altro riferimento immunitario
        extra = len(cs.index.difference(sub.obs_names))
        if extra:
            raise SystemExit(
                f"[ERRORE] la cache {args.cnv_cache} contiene {extra:,} cellule non selezionate "
                "(altro sottoinsieme o riferimento): cancellala per ricalcolare inferCNV."
            )
        method = "infercnv" if args.method == "auto" else args.method
        check_fingerprint(
            args.cnv_cache.name, cnv_fingerprint(sub, method), args.cnv_cache
        )
        sub.obs["cnv_score"] = score.values
    else:
        method = "infercnv" if args.method == "auto" else args.method
        if method == "infercnv":
            if not add_positions(sub):
                raise SystemExit(
                    f"[ERRORE] coordinate dei geni (Ensembl {ENSEMBL_RELEASE}) non "
                    f"disponibili: nessun host BioMart risponde e {GENE_POS} non esiste. "
                    "copykat non e' usato come ripiego (richiede rpy2 e R)."
                )
        else:
            require_copykat()
        t0 = time.time()
        score = run_infercnv(sub) if method == "infercnv" else run_copykat(sub)
        sub.obs["cnv_score"] = score
        print(f"[ok] CNV in {time.time() - t0:.0f}s")
        sub.obs[["cnv_score"]].to_csv(args.cnv_cache)
        save_fingerprint(args.cnv_cache.name, cnv_fingerprint(sub, method))

    # Soglia per paziente: mediana + 3 * MAD normalizzata sulle immunitarie del paziente.
    # Le immunitarie sono anche il riferimento su cui inferCNV centra lo score: le duttali
    # con un qualsiasi programma epiteliale tendono a superarla (vedi docstring, punto 1).
    print("\n" + "=" * 74)
    print("ETICHETTA DI MALIGNITA' (soglia per paziente sugli immuni di riferimento)")
    print("=" * 74)
    # NaN = paziente senza soglia (troppo poche immunitarie): duttali senza etichetta
    sub.obs["malignant"] = np.nan
    rows = []
    for pat, g in sub.obs.groupby("patients"):
        imm = g[g.all_celltype.isin(IMMUNE_REF)].cnv_score
        if len(imm) < 20:
            print(
                f"  [!] {pat}: {len(imm)} immuni (< 20), nessuna soglia: "
                "duttali senza etichetta, escluse dalla classificazione"
            )
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
    if lab.empty:
        raise SystemExit(
            "[ERRORE] nessun paziente con almeno 20 immunitarie: nessuna soglia."
        )
    print(lab.round(3).to_string(index=False))
    lab.to_csv(TAB / "malignancy_labels.csv", index=False)

    # Pazienti privi di contrasto binario interno
    tum = lab[lab.tissue == "Tumor"]
    degenerate = tum[(tum.frac_malig == 0) | (tum.frac_malig == 1)]
    if len(degenerate):
        print(
            f"\n[!!] {len(degenerate)} pazienti con duttali tumorali tutte uguali "
            f"(0% o 100% maligne): {degenerate.patient.tolist()}"
        )
        print("     non forniscono contrasto binario intra-paziente.")

    # Classificazione LOPO sulle duttali del tessuto tumorale (conteggi grezzi: la
    # normalizzazione e la selezione dei geni avvengono dentro la pipeline di ogni fold)
    print("\n" + "=" * 74)
    print("CLASSIFICATORE LOPO: predire la malignita' dal trascrittoma")
    print("=" * 74)
    train = sub[
        (sub.obs.all_celltype == "Ductal")
        & (sub.obs.tissue == "Tumor")
        & sub.obs.malignant.notna()
    ]
    counts = train.obs.patients.value_counts()
    good = counts[counts >= MIN_DUCTAL].index
    train = train[train.obs.patients.isin(good)]
    print(
        f"pazienti con >={MIN_DUCTAL} duttali tumorali: {len(good)} -> {sorted(good)}"
    )
    if train.n_obs == 0 or train.obs.malignant.nunique() < 2:
        raise SystemExit("[ERRORE] training LOPO vuoto o con una sola classe.")
    print(
        f"cellule di training: {train.n_obs:,} | maligne: {train.obs.malignant.mean():.1%}"
    )
    Xtr = sp.csr_matrix(train.X)
    ytr = train.obs.malignant.values.astype(int)
    gtr = train.obs.patients.values

    res, preds = classify_lopo(Xtr, ytr, gtr, n_hvg=args.n_hvg)
    print("\nper paziente (leave-one-patient-out):")
    print(res.round(3).to_string(index=False))
    res.to_csv(TAB / "malignancy_lopo.csv", index=False)
    if not res.empty:
        summ = res[["balanced_acc", "f1_malig", "auc"]].agg(["mean", "std", "count"])
        print(f"\nAUC media: {res.auc.mean():.3f} ± {res.auc.std():.3f}")
        print(
            f"balanced accuracy media: {res.balanced_acc.mean():.3f} ± {res.balanced_acc.std():.3f}"
        )
        summ.T.rename_axis("metric").to_csv(TAB / "malignancy_lopo_summary.csv")
    else:
        # niente riepilogo di un'esecuzione precedente in results/
        print(
            "[!!] nessun fold LOPO valutabile: malignancy_lopo_summary.csv non scritto"
        )
        (TAB / "malignancy_lopo_summary.csv").unlink(missing_ok=True)

    # Applicazione alle duttali del tessuto adiacente: modello finale stimato sulle sole
    # duttali tumorali di training, compresa la scelta dei geni
    print("\n" + "=" * 74)
    print("DUTTALI DEL TESSUTO ADIACENTE (mai in training; stessa esecuzione inferCNV)")
    print("=" * 74)
    ext = sub[
        (sub.obs.all_celltype == "Ductal") & (sub.obs.tissue == "Adjacent_normal")
    ]
    if ext.n_obs > 20:
        final = fit_final(Xtr, ytr, n_hvg=args.n_hvg)
        pred = final.predict(sp.csr_matrix(ext.X))
        print(f"duttali normali adiacenti: {ext.n_obs:,}")
        print(
            f"predette maligne: {pred.mean():.1%} | non maligne: {1 - pred.mean():.1%}"
        )
        print(
            "La frazione predetta maligna misura quante duttali adiacenti hanno il profilo"
        )
        print(
            "trascrizionale delle duttali tumorali etichettate maligne da inferCNV; non e' un"
        )
        print("tasso di falsi positivi")
        print("rispetto a una verita' indipendente (vedi docstring, punto 3).")
        pd.DataFrame({"pred": pred}, index=ext.obs_names).to_csv(
            TAB / "external_test_predictions.csv"
        )
        ext_rows = [
            {
                "patient": "_all",
                "n": int(ext.n_obs),
                "frac_pred_malignant": float(pred.mean()),
            }
        ]
        for p_, m_ in pd.Series(pred, index=ext.obs.patients.values).groupby(level=0):
            ext_rows.append(
                {"patient": p_, "n": len(m_), "frac_pred_malignant": float(m_.mean())}
            )
        pd.DataFrame(ext_rows).to_csv(TAB / "external_test_summary.csv", index=False)
    else:
        print("numero insufficiente di duttali normali (<= 20): passo saltato.")
        for name in ("external_test_predictions.csv", "external_test_summary.csv"):
            (TAB / name).unlink(missing_ok=True)

    # Figura diagnostica (distribuzione CNV score e AUC per paziente)
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
    if not res.empty:
        axes[1].bar(res.patient, res.auc, color="#7570b3", alpha=0.85)
        axes[1].axhline(0.5, ls="--", c="k", lw=0.8)
        axes[1].set_ylabel("AUC (LOPO)")
        axes[1].tick_params(axis="x", rotation=90)
        axes[1].set_title("classificazione malignità per paziente")
    fig.tight_layout()
    fig.savefig(FIG / "malignancy.png", dpi=140)
    plt.close(fig)
    print("\n[ok] figura: results/figures/malignancy.png")


if __name__ == "__main__":
    main()

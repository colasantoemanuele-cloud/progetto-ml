#!/usr/bin/env python3
"""
10_loco.py - Tasso d'errore delle etichette T con tre schemi di riferimento clonale e frazione di doppietti NK-TCR+.

Scopo
-----
1. Arbitrato con tre schemi di riferimento per l'identita' CD4/CD8 del clone (vedi
   09_arbitrate.py per le frazioni f8/f4 e il margine minimo 0.20):
   - `naive`: il riferimento include il compartimento giudicato (circolare). A
     differenza di 09_arbitrate.py, che tiene solo i cloni con >= 3 cellule in almeno
     2 compartimenti, `naive` tiene anche i cloni presenti in un solo compartimento
     (per i quali il riferimento e' solo il compartimento giudicato): per questo i
     tassi di `naive` (error_rate_schemes.csv) e di 09 (error_rate_by_compartment.csv)
     sono diversi;
   - `loco` (leave-one-compartment-out): il riferimento sono gli altri compartimenti
     del clone. Con cloni presenti in due soli compartimenti, giudicare il sangue
     usa il tumore come riferimento (in simulazione questo gonfia l'errore del sangue,
     vedi src/sim/sim_b_reference_schemes.py);
   - `blood` (schema primario): il riferimento e' solo il PBMC, e si giudicano solo
     Adjacent_normal e Tumor.
   Per ogni schema: errore per compartimento e differenze appaiate sugli stessi cloni
   (Tumor - Adjacent_normal, Tumor - PBMC), con IC 95% da cluster bootstrap sui pazienti
   (`scpipe.cluster_bootstrap`, "expanded percentile"). Con 4-5 pazienti l'intervallo
   non e' calibrato: si riportano anche il numero di ricampionature distinte e la
   frazione di repliche con differenza > 0.
   Etichetta di maggioranza del compartimento: moda (`scpipe.majority_label`). Se due
   o piu' etichette hanno lo stesso conteggio massimo l'etichetta non e' definita: la
   coppia clone-compartimento e' esclusa dal calcolo principale (colonna `tie` nelle
   tabelle di arbitrato, `n_ties` nei tassi) e, come sensibilita', contata come errore
   oppure come corretta (`error_rate_tie_sensitivity.csv`). La regola non dipende
   dall'ordine dei barcode.
   Presupposto: il lignaggio CD4/CD8 di un clone T maturo e' di norma stabile fra
   tessuti, con eccezioni note (CD4 citotossiche, CD4+CD8αα+ e doppi positivi nei
   tessuti, MAIT, NKT, γδ). Il controllo dei cloni a TCR semi-invariante (MAIT-like:
   TRAV1-2 con TRAJ33/12/20; iNKT: TRAV10 con TRAJ18) fra le unita' del Tumor e' in
   `clone_invariant_tcr_check.csv`.

2. Soglia di CD4 sopra l'RNA ambientale: 95-esimo percentile di CD4 nelle cellule
   Ductal ed Endothelial di ciascun compartimento. Nel PBMC queste cellule sono
   pochissime (11 nella versione pubblicata): la soglia del PBMC e' stimata su pochi
   punti ed e' instabile. Il numero di cellule e' stampato e salvato.

3. Robustezza: LOCO ripetuto escludendo le cellule con n_counts sopra il 90-esimo
   percentile della propria (libreria, classe): un proxy dei doppietti, non doppietti
   identificati.

4. Frazione di doppietti fisici nelle NK con TCR, da un unico modello di miscela:
   una NK-TCR+ e' con probabilita' d un doppietto NK + CD8T (UMI totali s + c, TCR
   della CD8T) e con probabilita' 1 - d una singola cellula con UMI totali s e TCR
   alla densita' di una T (rapporto r_CD8); s = n_counts(NK TCR-), c = n_counts(CD8T TCR+).
       d1 = (n_counts(NK TCR+) - s) / c                       (UMI totali)
       d2 = (r_CD8 - r_NK) / (r_CD8 - r_dbl),  r = UMI_TCR per 1000 UMI totali,
            r_dbl = r_CD8 * c / (s + c)                       (diluizione del TCR)
   d2 non e' indipendente da d1: usa le stesse s e c. La miscela di mediane (o di
   rapporti) e' un'approssimazione: le due stime sono indicative. Le stime
   principali usano le mediane pooled su tutti i tessuti; come sensibilita' si
   riportano le medie e le stime per compartimento. IC con bootstrap sui pazienti
   (500 repliche, "expanded percentile"). Le stesse stime per le sole NK-TCR+ dei
   cloni del Tumor contati come errore (schema blood) sono in
   `doublet_mixture_error_nk.csv`.

5. Sensibilita' (blood): scomposizione dell'errore nel Tumor (cloni CD4 nel sangue
   etichettati CD8T, cloni etichettati NK, altri), tasso escludendo le etichette NK,
   e regole d'identita' piu' severe per CD8+ (solo CD8B > 0; CD8A >= 2 UMI;
   CD8B > 0 oppure CD8A >= 2 UMI).

Input
-----
  --h5ad    data/interim/sc_raw.h5ad
  --tcrdir  cartella con i 30 file VDJ

Output
------
  results/tables/arbitration_{naive,loco,blood}.csv   una riga per (clone, compartimento giudicato);
                                                      `tie` = etichetta in pareggio (errore NaN)
  results/tables/error_rate_schemes.csv               errore per schema x compartimento (+ varianti),
                                                      pareggi esclusi (`n_ties`)
  results/tables/error_rate_tie_sensitivity.csv       per schema: errore per compartimento e
                                                      differenza Tumor - Adjacent_normal con i
                                                      pareggi esclusi, contati errore, contati corretti
  results/tables/clone_invariant_tcr_check.csv        cloni MAIT-like / iNKT fra le unita' del Tumor
                                                      (errore, corretto, pareggio) e frazione fra le
                                                      NK-TCR+
  results/tables/doublet_mixture_error_nk.csv         d1, d2 con IC per le NK-TCR+ dei cloni del
                                                      Tumor contati come errore (schema blood)
  results/tables/error_rate_paired_diff.csv           differenze appaiate per schema, con cloni e pazienti
  results/tables/error_decomposition_blood.csv        scomposizione dell'errore nel Tumor (blood)
  results/tables/error_rate_identity_rules.csv        errore (blood) con regole CD8+ alternative
  results/tables/ambient_cd4_threshold.csv            soglia CD4 per compartimento e numero di cellule
  results/tables/doublet_mixture.csv                  d1, d2 con IC: mediane/medie, pooled/per compartimento
  results/tables/doublet_mixture_boot.csv             repliche bootstrap di d1, d2 (mediane, pooled;
                                                      cancellata se la stima pooled non e' possibile)
  results/figures/error_rate_loco.png  (schema loco; barre con <= 5 pazienti grigie
                                       tratteggiate: IC non calibrato; cancellata se lo schema
                                       non ha compartimenti valutabili)

Argomenti
---------
  --h5ad (obbligatorio), --tcrdir (default data/raw/tcr)

Uso
---
    python src/10_loco.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import (  # noqa: E402
    SMALL_K,
    cluster_bootstrap,
    majority_label,
    percentile_ci,
)
from scpipe import small_k_diagnostics as _small_k  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_BOOT = 2000
N_BOOT_MIX = 500
# Con n pazienti <= SMALL_K (scpipe, = 5) l'IC del cluster bootstrap non e' calibrato

# Colonne di diagnostica aggiunte ai contrasti con <= SMALL_K pazienti
DIAG_COLS = ["n_distinct_resamples", "frac_boot_gt0"]

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}
COMPS = ["PBMC", "Adjacent_normal", "Tumor"]
MARKERS = ["CD3D", "CD4", "CD8A", "CD8B", "NKG7", "GNLY"]
# Popolazioni epiteliali/endoteliali prive di espressione endogena di CD4/CD8 (proxy di mRNA ambientale)
AMBIENT_PROXY = ["Ductal", "Endothelial"]
FNAME_RE = re.compile(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz")

MIN_CELLS = 3
MIN_MARGIN = 0.20
# Regole per chiamare una cellula CD8+ (la prima e' quella principale)
CD8_RULES = {
    "cd8a_or_cd8b": "CD8A > 0 oppure CD8B > 0 (principale)",
    "cd8b_only": "CD8B > 0",
    "cd8a_ge2": "CD8A >= 2 UMI",
    "cd8b_or_cd8a_ge2": "CD8B > 0 oppure CD8A >= 2 UMI",
}


def load_tcr(tcrdir: Path) -> pd.DataFrame:
    """Carica le annotazioni VDJ produttive e restituisce il clonotipo CDR3(TRB) e gli UMI TCR totali."""
    frames = []
    for p in sorted(tcrdir.glob("*_filtered_contig_annotations.csv.gz")):
        m = FNAME_RE.match(p.name)
        if not m:
            continue
        pnum, tis = m.groups()
        df = pd.read_csv(p, compression="gzip", low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df["patient"] = f"PA{int(pnum):02d}"
        df["tissue"] = TISSUE_MAP[tis]
        df["barcode16"] = df.barcode.astype(str).str.extract(
            r"([ACGT]{16})", expand=False
        )
        frames.append(
            df[["patient", "tissue", "barcode16", "chain", "cdr3_nt", "umis"]]
        )
    if not frames:
        raise SystemExit(f"[ERRORE] nessun file VDJ riconosciuto in {tcrdir}")
    c = pd.concat(frames, ignore_index=True).dropna(subset=["barcode16"])
    c = c.sort_values("umis", ascending=False).drop_duplicates(
        ["patient", "tissue", "barcode16", "chain"]
    )
    seq = c.pivot_table(
        index=["patient", "tissue", "barcode16"],
        columns="chain",
        values="cdr3_nt",
        aggfunc="first",
    )
    umi = c.pivot_table(
        index=["patient", "tissue", "barcode16"],
        columns="chain",
        values="umis",
        aggfunc="first",
    )
    for ch in ("TRA", "TRB"):
        if ch not in seq.columns:
            seq[ch], umi[ch] = np.nan, np.nan
    w = (
        seq.rename(columns={"TRB": "clone_b"})
        .join(umi.rename(columns={"TRA": "umi_TRA", "TRB": "umi_TRB"}))
        .reset_index()
    )
    w["umi_tcr"] = w[["umi_TRA", "umi_TRB"]].sum(axis=1, min_count=1)
    return w[["patient", "tissue", "barcode16", "clone_b", "umi_tcr"]]


def cluster_boot(by_pat: dict[str, np.ndarray], n: int = N_BOOT) -> tuple:
    """Media pooled e IC 95% con cluster bootstrap sui pazienti (pazienti vuoti esclusi)."""
    return cluster_bootstrap({p: v for p, v in by_pat.items() if len(v)}, n, SEED)


def small_k_diagnostics(by_pat: dict[str, np.ndarray], n: int = N_BOOT) -> dict:
    """Ricampionature distinte e frazione di repliche > 0 (`scpipe.small_k_diagnostics`).

    Stesso seme e stessa sequenza di estrazioni di `cluster_boot`.
    """
    return _small_k({p: v for p, v in by_pat.items() if len(v)}, n, SEED)


def ambient_cd4(gex: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """Soglia di fondo di CD4 per tessuto: 95-esimo percentile su cellule Ductal ed Endothelial."""
    print("\n" + "=" * 78)
    print("(2) AMBIENT: CD4 nelle cellule che NON possono esprimerlo")
    print("=" * 78)
    amb = gex[gex.all_celltype.isin(AMBIENT_PROXY)]
    tab = amb.groupby(["tissue", "all_celltype"], observed=True).g_CD4.agg(
        n="size",
        mean="mean",
        frac_pos=lambda s: (s > 0).mean(),
        p95=lambda s: s.quantile(0.95),
    )
    print(tab.round(4).to_string())

    # La soglia e' un quantile empirico: con poche cellule proxy (nel PBMC sono una
    # decina) e' stimata male; il numero di cellule usato viene riportato accanto.
    g = amb.groupby("tissue", observed=True).g_CD4
    thr = g.quantile(0.95)
    out = pd.DataFrame({"threshold_p95": thr, "n_cells": g.size()})
    print("\nsoglia CD4 per compartimento (95o pct dell'ambient) e cellule usate:")
    print(out.round(4).to_string())
    for tis, n in out.n_cells.items():
        if n < 30:
            print(f"  [!] {tis}: soglia stimata su {n} cellule, instabile")

    print(
        "\nCD4 medio nelle CD8T (tutte), per compartimento (verifica rispetto alla soglia):"
    )
    print(
        gex[gex.all_celltype == "CD8T"]
        .groupby("tissue", observed=True)
        .g_CD4.mean()
        .round(4)
        .to_string()
    )
    print(
        "\nLe popolazioni mieloidi esprimono CD4 endogeno e sono escluse dal proxy ambientale."
    )
    return thr, out


def cd8_positive(x: pd.DataFrame, rule: str) -> pd.Series:
    """Cellule CD8+ secondo la regola `rule` (chiavi di CD8_RULES)."""
    if rule == "cd8a_or_cd8b":
        return (x.g_CD8A > 0) | (x.g_CD8B > 0)
    if rule == "cd8b_only":
        return x.g_CD8B > 0
    if rule == "cd8a_ge2":
        return x.u_CD8A >= 2
    if rule == "cd8b_or_cd8a_ge2":
        return (x.g_CD8B > 0) | (x.u_CD8A >= 2)
    raise ValueError(rule)


def per_comp_stats(
    lym: pd.DataFrame, cd4_thr: pd.Series | None, rule: str = "cd8a_or_cd8b"
) -> pd.DataFrame:
    """Per ogni (paziente, clone, tessuto) con >= MIN_CELLS cellule: frazioni CD8+ e CD4+ ed etichetta di maggioranza."""
    x = lym.copy()
    x["cd8_pos"] = cd8_positive(x, rule)
    if cd4_thr is None:
        x["cd4_pos"] = (x.g_CD4 > 0) & ~x.cd8_pos
    else:
        x["cd4_pos"] = (x.g_CD4 > x.tissue.map(cd4_thr)) & ~x.cd8_pos
    per = (
        x.groupby(["patient", "clone_b", "tissue"], observed=True)
        .agg(
            n=("cd8_pos", "size"),
            f8=("cd8_pos", "mean"),
            f4=("cd4_pos", "mean"),
            # moda; None se piu' etichette hanno lo stesso conteggio massimo (pareggio)
            lab=("all_celltype", majority_label),
            lab_top=("all_celltype", top_labels),
        )
        .reset_index()
    )
    return per[per.n >= MIN_CELLS]


def top_labels(labels) -> str:
    """Etichette con il conteggio massimo, ordinate e unite da '|' (piu' di una = pareggio)."""
    vc = pd.Series(list(labels)).value_counts()
    return "|".join(sorted(vc.index[vc == vc.iloc[0]].astype(str)))


def arbitrate(per: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """Confronta l'etichetta di ciascun compartimento con l'identita' del clone ('naive', 'loco', 'blood')."""
    rows = []
    for (pat, clone), g in per.groupby(["patient", "clone_b"], observed=True):
        for _, r in g.iterrows():
            target = r.tissue
            if scheme == "naive":
                ref = g
            elif scheme == "loco":
                ref = g[g.tissue != target]
            elif scheme == "blood":
                if target == "PBMC":
                    continue
                ref = g[g.tissue == "PBMC"]
            else:
                raise ValueError(scheme)
            if ref.empty:
                continue
            # Media non ponderata delle frazioni sui compartimenti di riferimento
            f8, f4 = ref.f8.mean(), ref.f4.mean()
            if abs(f8 - f4) < MIN_MARGIN:
                continue
            ident = "CD8T" if f8 > f4 else "CD4T"
            tie = r.lab is None or pd.isna(r.lab)
            if tie:
                # pareggio: etichetta di maggioranza non definita, errore non definito;
                # l'unita' e' esclusa dal calcolo principale (vedi tie_sensitivity)
                if not set(r.lab_top.split("|")) & {"CD4T", "CD8T", "NK"}:
                    continue
            elif r.lab not in ("CD4T", "CD8T", "NK"):
                continue
            rows.append(
                {
                    "patient": pat,
                    "clone": clone,
                    "tissue": target,
                    "lab": None if tie else r.lab,
                    "identity": ident,
                    "error": np.nan if tie else float(r.lab != ident),
                    "tie": tie,
                    "tie_labels": r.lab_top if tie else "",
                    "n_ref_comp": ref.tissue.nunique(),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "patient",
            "clone",
            "tissue",
            "lab",
            "identity",
            "error",
            "tie",
            "tie_labels",
            "n_ref_comp",
        ],
    )


def untied(j: pd.DataFrame) -> pd.DataFrame:
    """Unita' con etichetta di maggioranza definita (calcolo principale)."""
    if "tie" not in j.columns:
        return j
    return j[~j.tie.astype(bool)]


def report(j: pd.DataFrame, title: str) -> pd.DataFrame:
    """Errore per compartimento con IC (cluster bootstrap); stampa anche le differenze appaiate."""
    print(f"\n--- {title} ---")
    if j.empty:
        print("  nessun clone risolvibile.")
        return pd.DataFrame(columns=["tissue", "n", "n_pat", "error", "lo", "hi"])
    rows = []
    for tis in COMPS:
        g_all = j[j.tissue == tis]
        g = untied(g_all)
        n_ties = len(g_all) - len(g)
        if g.empty or g.patient.nunique() < 3:
            continue
        by = {p: gg.error.values.astype(float) for p, gg in g.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        row = {
            "tissue": tis,
            "n": len(g),
            "n_ties": n_ties,
            "n_pat": g.patient.nunique(),
            "error": m,
            "lo": lo,
            "hi": hi,
        }
        if row["n_pat"] <= SMALL_K:
            # per un tasso (sempre >= 0) la frazione di repliche > 0 non e' informativa
            diag = small_k_diagnostics(by)
            row["n_distinct_resamples"] = diag["n_distinct_resamples"]
        rows.append(row)
        print(
            f"  {tis:16s} errore = {m:.3f}  IC95% [{lo:.3f}, {hi:.3f}]  "
            f"(n={len(g)}, {g.patient.nunique()} pz; {n_ties} pareggi esclusi)"
        )
    paired_diffs(j, verbose=True)
    return pd.DataFrame(rows)


def paired_diffs(j: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """Differenze appaiate dell'errore sugli stessi cloni (Tumor - Adjacent_normal, Tumor - PBMC)."""
    cols = ["contrast", "n_clones", "n_patients", "diff", "lo", "hi"]
    rows = []
    j = untied(j)  # le unita' in pareggio non hanno errore definito
    if j.empty:
        return pd.DataFrame(rows, columns=cols)
    piv = j.pivot_table(
        index=["patient", "clone"], columns="tissue", values="error", aggfunc="max"
    )
    for a, b in [("Tumor", "Adjacent_normal"), ("Tumor", "PBMC")]:
        if not {a, b}.issubset(piv.columns):
            continue
        d = piv.dropna(subset=[a, b]).reset_index()
        if len(d) < 20 or d.patient.nunique() < 3:
            if verbose:
                print(
                    f"  {a} - {b}: n insufficiente ({len(d)} cloni, {d.patient.nunique()} pz)"
                )
            continue
        d["diff"] = d[a].astype(float) - d[b].astype(float)
        by = {p: g["diff"].values for p, g in d.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        row = {
            "contrast": f"{a} - {b}",
            "n_clones": len(d),
            "n_patients": d.patient.nunique(),
            "diff": m,
            "lo": lo,
            "hi": hi,
        }
        if row["n_patients"] <= SMALL_K:
            row.update(small_k_diagnostics(by))
        rows.append(row)
        if verbose:
            v = ci_verdict(lo, hi, row["n_patients"])
            print(
                f"  {a} - {b} (stessi cloni): {m:+.3f}  IC95% [{lo:+.3f}, {hi:+.3f}]  {v}"
                f"  (n={len(d)} cloni, {d.patient.nunique()} pz)"
            )
    return pd.DataFrame(
        rows, columns=cols + [c for c in DIAG_COLS if any(c in r for r in rows)]
    )


def tie_sensitivity(j: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """Errore e differenza appaiata con i pareggi esclusi, contati come errore e come corretti.

    Versione principale: pareggi esclusi. Sensibilita': ogni unita' in pareggio contata
    come errore (1) oppure come corretta (0). Stesso bootstrap sui pazienti di `report`.
    """
    rows = []
    if j.empty:
        return pd.DataFrame(rows)
    tie = j.tie.astype(bool) if "tie" in j.columns else pd.Series(False, index=j.index)
    versions = {
        "ties_excluded": j[~tie],
        "ties_as_error": j.assign(error=j.error.where(~tie, 1.0)),
        "ties_as_correct": j.assign(error=j.error.where(~tie, 0.0)),
    }
    for version, jj in versions.items():
        for tis in COMPS:
            g = jj[jj.tissue == tis]
            if g.empty or g.patient.nunique() < 3:
                continue
            by = {p: gg.error.values.astype(float) for p, gg in g.groupby("patient")}
            m, lo, hi = cluster_boot(by)
            rows.append(
                {
                    "scheme": scheme,
                    "version": version,
                    "quantity": f"error_{tis}",
                    "n_units": len(g),
                    "n_ties": int((tie & (j.tissue == tis)).sum()),
                    "n_patients": g.patient.nunique(),
                    "value": m,
                    "lo": lo,
                    "hi": hi,
                }
            )
        # differenza appaiata Tumor - Adjacent_normal sugli stessi cloni
        piv = jj.pivot_table(
            index=["patient", "clone"], columns="tissue", values="error", aggfunc="max"
        )
        if {"Tumor", "Adjacent_normal"}.issubset(piv.columns):
            d = piv.dropna(subset=["Tumor", "Adjacent_normal"]).reset_index()
            if len(d) >= 20 and d.patient.nunique() >= 3:
                d["diff"] = d.Tumor - d.Adjacent_normal
                by = {p: g["diff"].values for p, g in d.groupby("patient")}
                m, lo, hi = cluster_boot(by)
                rows.append(
                    {
                        "scheme": scheme,
                        "version": version,
                        "quantity": "diff_Tumor_minus_Adjacent_normal",
                        "n_units": len(d),
                        "n_ties": np.nan,
                        "n_patients": d.patient.nunique(),
                        "value": m,
                        "lo": lo,
                        "hi": hi,
                    }
                )
    return pd.DataFrame(rows)


def load_tra_genes(tcrdir: Path) -> pd.DataFrame:
    """Geni V e J della catena TRA dominante (piu' UMI) per cellula, dai file VDJ."""
    frames = []
    for p in sorted(tcrdir.glob("*_filtered_contig_annotations.csv.gz")):
        m = FNAME_RE.match(p.name)
        if not m:
            continue
        pnum, tis = m.groups()
        df = pd.read_csv(p, compression="gzip", low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[df.chain == "TRA"]
        df["patient"] = f"PA{int(pnum):02d}"
        df["tissue"] = TISSUE_MAP[tis]
        df["barcode16"] = df.barcode.astype(str).str.extract(
            r"([ACGT]{16})", expand=False
        )
        frames.append(
            df[["patient", "tissue", "barcode16", "v_gene", "j_gene", "umis"]]
        )
    if not frames:
        return pd.DataFrame(
            columns=["patient", "tissue", "barcode16", "v_gene", "j_gene"]
        )
    c = pd.concat(frames, ignore_index=True).dropna(subset=["barcode16"])
    c = c.sort_values("umis", ascending=False).drop_duplicates(
        ["patient", "tissue", "barcode16"]
    )
    return c[["patient", "tissue", "barcode16", "v_gene", "j_gene"]]


def invariant_tcr_check(
    lym: pd.DataFrame, j: pd.DataFrame, tra: pd.DataFrame
) -> pd.DataFrame:
    """Cloni a TCR semi-invariante fra le unita' del Tumor (schema blood).

    MAIT-like: TRAV1-2 con TRAJ33, TRAJ12 o TRAJ20; iNKT: TRAV10 con TRAJ18 (catena TRA
    dominante del clone nelle cellule linfoidi). Controlla se le eccezioni note alla
    stabilita' del lignaggio CD4/CD8 spiegano gli errori.
    """
    x = lym.merge(tra, on=["patient", "tissue", "barcode16"], how="left")
    x["mait"] = (x.v_gene == "TRAV1-2") & x.j_gene.isin(["TRAJ33", "TRAJ12", "TRAJ20"])
    x["inkt"] = (x.v_gene == "TRAV10") & (x.j_gene == "TRAJ18")
    # un clone e' MAIT-like / iNKT se lo e' la maggioranza delle sue cellule con TRA
    cl = (
        x[x.v_gene.notna()]
        .groupby(["patient", "clone_b"], observed=True)[["mait", "inkt"]]
        .mean()
        .gt(0.5)
        .reset_index()
        .rename(columns={"clone_b": "clone"})
    )
    t = j[j.tissue == "Tumor"].merge(cl, on=["patient", "clone"], how="left")
    t[["mait", "inkt"]] = t[["mait", "inkt"]].fillna(False).astype(bool)
    tie = t.tie.astype(bool)
    rows = []
    for name, mask in (
        ("MAIT_like", t.mait),
        ("iNKT", t.inkt),
        ("altri", ~t.mait & ~t.inkt),
    ):
        rows.append(
            {
                "category": name,
                "n_units_tumor": int(mask.sum()),
                "n_error": int((mask & ~tie & (t.error == 1)).sum()),
                "n_correct": int((mask & ~tie & (t.error == 0)).sum()),
                "n_tie": int((mask & tie).sum()),
            }
        )
    nk = x[(x.all_celltype == "NK") & x.v_gene.notna()]
    for name, col in (("MAIT_like", "mait"), ("iNKT", "inkt")):
        rows.append(
            {
                "category": f"{name}_fraction_of_NK_TCR+_cells",
                "n_units_tumor": np.nan,
                "n_error": np.nan,
                "n_correct": np.nan,
                "n_tie": np.nan,
                "fraction": float(nk[col].mean()) if len(nk) else np.nan,
                "n_cells": len(nk),
            }
        )
    return pd.DataFrame(rows)


def ci_verdict(lo: float, hi: float, n_patients: int) -> str:
    """Lettura dell'IC; con <= SMALL_K pazienti l'IC non e' calibrato e non si da' verdetto."""
    if n_patients <= SMALL_K:
        return f"IC non calibrato ({n_patients} pazienti)"
    return "IC esclude 0" if lo > 0 or hi < 0 else "IC include 0"


def decompose_tumor_error(j: pd.DataFrame) -> pd.DataFrame:
    """Scompone l'errore nel Tumor (schema blood) per tipo di etichetta sbagliata."""
    t_all = j[j.tissue == "Tumor"]
    t = untied(t_all)
    t = t.assign(error=t.error.astype(bool))
    if t.empty:
        return pd.DataFrame()
    n = len(t)
    cats = {
        "CD4_in_blood_labelled_CD8T": (t.identity == "CD4T") & (t.lab == "CD8T"),
        "labelled_NK": (t.lab == "NK"),
        "CD8_in_blood_labelled_CD4T": (t.identity == "CD8T") & (t.lab == "CD4T"),
    }
    rows = [
        {
            "category": k,
            "n_errors": int(v.sum()),
            "rate_over_all_clones": v.sum() / n,
            "n_clones": n,
            "n_patients": t.patient.nunique(),
        }
        for k, v in cats.items()
    ]
    other = t.error & ~np.logical_or.reduce(list(cats.values()))
    rows.append(
        {
            "category": "other",
            "n_errors": int(other.sum()),
            "rate_over_all_clones": other.sum() / n,
            "n_clones": n,
            "n_patients": t.patient.nunique(),
        }
    )
    # Transizioni identita' (dal sangue) -> etichetta nel tumore, per tutti gli errori
    for (ident, lab), v in t[t.error].groupby(["identity", "lab"]).size().items():
        rows.append(
            {
                "category": f"transition_{ident}_to_{lab}",
                "n_errors": int(v),
                "rate_over_all_clones": v / n,
                "n_clones": n,
                "n_patients": t.patient.nunique(),
            }
        )
    # Unita' in pareggio, escluse dal denominatore
    rows.append(
        {
            "category": "ties_excluded",
            # errore non definito per le unita' in pareggio
            "n_errors": np.nan,
            "rate_over_all_clones": np.nan,
            "n_clones": len(t_all) - n,
            "n_patients": t_all[
                t_all.tie.astype(bool) if "tie" in t_all else t_all.index.isin([])
            ].patient.nunique(),
        }
    )
    # Tasso escludendo i cloni etichettati NK, con IC sui pazienti
    nn = t[t.lab != "NK"]
    by = {p: g.error.values.astype(float) for p, g in nn.groupby("patient")}
    m, lo, hi = cluster_boot(by)
    rows.append(
        {
            "category": "error_rate_excluding_NK_labels",
            "n_errors": int(nn.error.sum()),
            "rate_over_all_clones": m,
            "lo": lo,
            "hi": hi,
            "n_clones": len(nn),
            "n_patients": nn.patient.nunique(),
        }
    )
    return pd.DataFrame(rows)


def mix_est(sub: pd.DataFrame, stat: str = "median") -> tuple[float, float]:
    """d1 e d2 su un sottoinsieme (NK e CD8T con colonna `ratio`); NaN se < 10 cellule per gruppo."""
    nk_pos = sub[(sub.all_celltype == "NK") & sub.has_tcr]
    nk_neg = sub[(sub.all_celltype == "NK") & ~sub.has_tcr]
    cd8 = sub[(sub.all_celltype == "CD8T") & sub.has_tcr]
    if min(len(nk_pos), len(nk_neg), len(cd8)) < 10:
        return np.nan, np.nan
    f = (lambda s: s.median()) if stat == "median" else (lambda s: s.mean())
    s, c, o = f(nk_neg.n_counts), f(cd8.n_counts), f(nk_pos.n_counts)
    # Modello comune (docstring, punto 4): doppietto = NK + CD8T con UMI s + c e il TCR
    # della CD8T; non doppietto = una cellula con UMI s e TCR alla densita' r_CD8.
    # d1: eccesso di UMI totali delle NK-TCR+ rispetto al singoletto NK, in unita' di CD8T
    d1 = (o - s) / max(c, 1)
    r_cd8, r_nk = f(cd8.ratio), f(nk_pos.ratio)
    # rapporto atteso in un doppietto: gli UMI TCR della CD8T diluiti su s + c
    r_dbl = r_cd8 * c / max(s + c, 1)
    # d2: r_NK = (1 - d) r_CD8 + d r_dbl  ->  d = (r_CD8 - r_NK) / (r_CD8 - r_dbl)
    d2 = (r_cd8 - r_nk) / max(r_cd8 - r_dbl, 1e-9)
    return float(np.clip(d1, 0, 1)), float(np.clip(d2, 0, 1))


def mix_boot(sub: pd.DataFrame, stat: str) -> tuple[list, list, int]:
    """Repliche di d1 e d2 ricampionando i pazienti (N_BOOT_MIX, seme SEED)."""
    pats = sub.patients.unique()
    parts = {p: sub[sub.patients == p] for p in pats}
    rng = np.random.default_rng(SEED)
    b1, b2 = [], []
    for _ in range(N_BOOT_MIX):
        draw = rng.choice(pats, len(pats), replace=True)
        a, b = mix_est(pd.concat([parts[p] for p in draw]), stat)
        if not np.isnan(a):
            b1.append(a)
            b2.append(b)
    return b1, b2, len(pats)


def mixture_error_nk(gex: pd.DataFrame, jblood: pd.DataFrame) -> pd.DataFrame:
    """d1 e d2 per le cellule NK-TCR+ dei cloni del Tumor contati come errore (schema blood).

    Le NK-TCR- e le CD8T-TCR+ di riferimento sono quelle del Tumor; IC con bootstrap
    sui pazienti come in `mixture`. I cloni NK in errore vengono da pochi pazienti
    (`n_patients_with_error_nk`): l'IC non e' calibrato.
    """
    err = jblood[
        (jblood.tissue == "Tumor") & (jblood.lab == "NK") & (jblood.error == 1)
    ][["patient", "clone"]]
    g = gex[(gex.tissue == "Tumor") & gex.all_celltype.isin(["NK", "CD8T"])].copy()
    g["ratio"] = g.umi_tcr / np.maximum(g.n_counts, 1) * 1000
    key = set(zip(err.patient, err.clone))
    in_err = pd.Series(
        [(p, c) in key for p, c in zip(g.patients, g.clone_b)], index=g.index
    )
    nk_pos = (g.all_celltype == "NK") & g.has_tcr
    sub = g[~nk_pos | in_err]  # NK-TCR+ solo dei cloni in errore
    rows = []
    for stat in ("median", "mean"):
        d1, d2 = mix_est(sub, stat)
        b1, b2, k = mix_boot(sub, stat)
        lo1, hi1 = percentile_ci(np.array(b1), k) if b1 else (np.nan, np.nan)
        lo2, hi2 = percentile_ci(np.array(b2), k) if b2 else (np.nan, np.nan)
        rows.append(
            {
                "statistic": stat,
                "n_nk_tcr_pos_cells": int((nk_pos & in_err).sum()),
                "n_clones": len(key),
                # pazienti ricampionati (tutto il Tumor) e pazienti con cloni NK in errore
                "n_patients": k,
                "n_patients_with_error_nk": err.patient.nunique(),
                "ic_calibrato": False,
                "d1": d1,
                "d1_lo": lo1,
                "d1_hi": hi1,
                "d2": d2,
                "d2_lo": lo2,
                "d2_hi": hi2,
            }
        )
    out = pd.DataFrame(rows)
    print(
        "\nd1/d2 per le NK-TCR+ dei cloni del Tumor contati come errore (schema blood):"
    )
    print(out.round(3).to_string(index=False))
    return out


def mixture(gex: pd.DataFrame) -> None:
    """Frazione di doppietti nelle NK-TCR+ con le stime d1 e d2 (non indipendenti)."""
    print("\n" + "=" * 78)
    print("(4) FRAZIONE DI DOPPIETTI NELLE NK-TCR+  (due stime, non indipendenti)")
    print("=" * 78)
    # Solo NK e CD8T entrano nelle stime: si restringe la tabella per velocizzare il bootstrap
    g = gex[gex.all_celltype.isin(["NK", "CD8T"])].copy()
    g["ratio"] = g.umi_tcr / np.maximum(g.n_counts, 1) * 1000

    rows = []
    strata = [("all", g)] + [(t, g[g.tissue == t]) for t in COMPS]
    for stat in ("median", "mean"):
        for name, sub in strata:
            d1, d2 = mix_est(sub, stat)
            if np.isnan(d1):
                print(f"  {stat:6s} {name:16s} cellule insufficienti (<10 per gruppo)")
                continue
            b1, b2, k = mix_boot(sub, stat)
            lo1, hi1 = percentile_ci(np.array(b1), k)
            lo2, hi2 = percentile_ci(np.array(b2), k)
            rows.append(
                {
                    "statistic": stat,
                    "stratum": name,
                    "n_patients": k,
                    "n_boot_valid": len(b1),
                    "d1": d1,
                    "d1_lo": lo1,
                    "d1_hi": hi1,
                    "d2": d2,
                    "d2_lo": lo2,
                    "d2_hi": hi2,
                }
            )
            print(
                f"  {stat:6s} {name:16s} d1 = {d1:.3f} [{lo1:.3f}, {hi1:.3f}]  "
                f"d2 = {d2:.3f} [{lo2:.3f}, {hi2:.3f}]"
            )
            if stat == "median" and name == "all":
                pd.DataFrame({"d1": b1, "d2": b2}).to_csv(
                    TAB / "doublet_mixture_boot.csv", index=False
                )
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "doublet_mixture.csv", index=False)
    has_pooled = (
        len(out) > 0 and ((out.statistic == "median") & (out.stratum == "all")).any()
    )
    if not has_pooled:
        # nessuna stima pooled: niente repliche di un'esecuzione precedente in results/
        (TAB / "doublet_mixture_boot.csv").unlink(missing_ok=True)
    if out.empty:
        print("  nessuno strato con abbastanza cellule: d1 e d2 non stimati.")
        return
    main_row = out[(out.statistic == "median") & (out.stratum == "all")]
    if len(main_row):
        r = main_row.iloc[0]
        print(
            f"\nStima principale (mediane, tutti i tessuti): d1 = {r.d1:.3f} "
            f"[{r.d1_lo:.3f}, {r.d1_hi:.3f}], d2 = {r.d2:.3f} [{r.d2_lo:.3f}, {r.d2_hi:.3f}]."
        )
        if min(r.d1_lo, r.d2_lo) > 0.5:
            print("Entrambe le stime indicano che la maggioranza delle NK-TCR+ e'")
            print("compatibile con doppietti NK + T (sotto il modello di miscela).")
        elif max(r.d1_hi, r.d2_hi) < 0.5:
            print(
                "Entrambe le stime indicano che i doppietti sono una minoranza delle NK-TCR+."
            )
        else:
            print("Le stime non permettono di dire se i doppietti siano la maggioranza")
            print("delle NK-TCR+.")
        print("d1 e d2 usano le stesse mediane di n_counts: il loro accordo non e' una")
        print(
            "conferma indipendente. Le mediane sono pooled fra tessuti (vedi strati)."
        )
    # Strati (compartimenti, mediane e medie) in cui la maggioranza non e' sostenuta:
    # la conclusione pooled non vale per ogni compartimento
    low = out[(out.stratum != "all") | (out.statistic != "median")]
    low = low[
        (low[["d1", "d2"]].min(axis=1) < 0.5)
        | (low[["d1_lo", "d2_lo"]].min(axis=1) < 0.5)
    ]
    if len(low):
        print("\n[!!] strati con una stima o un limite inferiore dell'IC sotto 0,5:")
        for r in low.itertuples():
            print(
                f"  {r.statistic:6s} {r.stratum:16s} d1 = {r.d1:.3f} [{r.d1_lo:.3f}, {r.d1_hi:.3f}]  "
                f"d2 = {r.d2:.3f} [{r.d2_lo:.3f}, {r.d2_hi:.3f}]"
            )
        print(
            "  In questi strati la maggioranza di doppietti non e' sostenuta dai dati."
        )


def bar_with_calibration(ax, labels, est, lo, hi, n_pat) -> None:
    """Barre con IC; quelle con <= SMALL_K pazienti grigie e tratteggiate (IC non calibrato).

    Usata anche da 21_fig_blood_scheme.py.
    """
    from matplotlib.patches import Patch

    est, lo, hi = (np.asarray(v, dtype=float) for v in (est, lo, hi))
    cal = np.asarray(n_pat) > SMALL_K
    bars = ax.bar(
        list(labels),
        est,
        yerr=[est - lo, hi - est],
        capsize=5,
        color=["#d95f02" if c else "#bbbbbb" for c in cal],
        alpha=0.85,
    )
    for b, c in zip(bars, cal):
        if not c:
            b.set_hatch("//")
    if not cal.all():
        ax.legend(
            handles=[
                Patch(
                    color="#d95f02", alpha=0.85, label=f"IC 95% (> {SMALL_K} pazienti)"
                ),
                Patch(
                    facecolor="#bbbbbb",
                    hatch="//",
                    label=f"IC non calibrato (<= {SMALL_K} pazienti)",
                ),
            ],
            fontsize=8,
        )


def main() -> None:
    """Punto di ingresso principale per l'arbitrato naive/LOCO/blood e i controlli di robustezza."""
    ap = argparse.ArgumentParser(
        description="Errore delle etichette T con tre schemi di riferimento e doppietti NK-TCR+.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--h5ad",
        type=Path,
        required=True,
        help="AnnData con conteggi grezzi (da 02_build_h5ad.py)",
    )
    ap.add_argument(
        "--tcrdir",
        type=Path,
        default=Path("data/raw/tcr"),
        help="cartella con i file VDJ filtered_contig_annotations",
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype", "barcode", "library"):
        obs[c] = obs[c].astype(str)
    obs = obs.rename(columns={"barcode": "barcode16"})
    sym = adata.var.gene_symbol.astype(str).values
    for gname in MARKERS:
        i = np.where(sym == gname)[0]
        if len(i):
            x = np.asarray(adata.X[:, i[0]].todense()).ravel()
            obs[f"g_{gname}"] = np.log1p(x / np.maximum(obs.n_counts.values, 1) * 1e4)
            if gname == "CD8A":
                obs["u_CD8A"] = x  # UMI grezzi, per la regola CD8A >= 2 UMI
    missing = [g for g in ("CD4", "CD8A", "CD8B") if f"g_{g}" not in obs]
    if missing:
        raise SystemExit(f"[ERRORE] marcatori assenti dall'h5ad: {missing}")

    tcr = load_tcr(args.tcrdir)
    gex = obs.merge(
        tcr,
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    print(f"[ok] {len(gex):,} cellule | TCR: {gex.has_tcr.sum():,}")

    thr, thr_tab = ambient_cd4(gex)
    thr_tab.rename_axis("tissue").to_csv(TAB / "ambient_cd4_threshold.csv")

    lym = gex[gex.has_tcr & gex.all_celltype.isin(["CD4T", "CD8T", "NK"])].copy()
    if lym.empty:
        raise SystemExit("[ERRORE] nessuna cellula linfoide con TCR.")

    print("\n" + "=" * 78)
    print("(1) SCHEMI DI RIFERIMENTO: naive, leave-one-compartment-out, solo sangue")
    print("=" * 78)
    per = per_comp_stats(lym, None)
    res, diffs, js, ties = [], [], {}, []
    for scheme, title in [
        ("naive", "NAIVE (riferimento include il target, circolare)"),
        ("loco", "LOCO (riferimento = gli altri compartimenti)"),
        ("blood", "BLOOD (riferimento = solo PBMC)"),
    ]:
        j = arbitrate(per, scheme)
        js[scheme] = j
        res.append(report(j, title).assign(scheme=scheme, variant="main"))
        diffs.append(paired_diffs(j).assign(scheme=scheme, variant="main"))
        ties.append(tie_sensitivity(j, scheme))
        j.to_csv(TAB / f"arbitration_{scheme}.csv", index=False)

    print("\n" + "=" * 78)
    print("(2b) RIPETUTO con soglia CD4 sopra l'ambient del compartimento")
    print("=" * 78)
    j = arbitrate(per_comp_stats(lym, thr), "loco")
    res.append(
        report(j, "LOCO + soglia CD4 corretta per ambient").assign(
            scheme="loco", variant="cd4_ambient_threshold"
        )
    )

    print("\n" + "=" * 78)
    print(
        "(3) RIPETUTO senza il decile superiore di n_counts per libreria x classe"
        " (proxy dei doppietti, non doppietti identificati)"
    )
    print("=" * 78)
    cut = gex.groupby(["library", "all_celltype"], observed=True).n_counts.transform(
        lambda s: s.quantile(0.90)
    )
    keep = gex.n_counts <= cut
    print(f"  cellule escluse: {(~keep).sum():,} / {len(gex):,}")
    lym_s = gex[
        keep & gex.has_tcr & gex.all_celltype.isin(["CD4T", "CD8T", "NK"])
    ].copy()
    j = arbitrate(per_comp_stats(lym_s, None), "loco")
    res.append(
        report(
            j, "LOCO senza il decile superiore di n_counts (proxy dei doppietti)"
        ).assign(scheme="loco", variant="no_top_decile")
    )

    print("\n" + "=" * 78)
    print("(5) SENSIBILITA' dello schema BLOOD")
    print("=" * 78)
    dec = decompose_tumor_error(js["blood"])
    print("scomposizione dell'errore nel Tumor:")
    print(dec.round(4).to_string(index=False))
    dec.to_csv(TAB / "error_decomposition_blood.csv", index=False)

    ts = pd.concat(ties, ignore_index=True)
    ts.to_csv(TAB / "error_rate_tie_sensitivity.csv", index=False)
    print("\nsensibilita' ai pareggi dell'etichetta di maggioranza (schema blood):")
    print(ts[ts.scheme == "blood"].round(4).to_string(index=False))

    inv = invariant_tcr_check(lym, js["blood"], load_tra_genes(args.tcrdir))
    inv.to_csv(TAB / "clone_invariant_tcr_check.csv", index=False)
    print("\ncloni a TCR semi-invariante nel Tumor (schema blood):")
    print(inv.round(4).to_string(index=False))

    rule_rows = []
    for rule, desc in CD8_RULES.items():
        jr = arbitrate(per_comp_stats(lym, None, rule=rule), "blood")
        r = report(jr, f"BLOOD con regola CD8+: {desc}")
        for row in r.to_dict("records"):
            rule_rows.append({"rule": rule, "description": desc, **row})
    pd.DataFrame(rule_rows).to_csv(TAB / "error_rate_identity_rules.csv", index=False)

    pd.concat(res, ignore_index=True).to_csv(
        TAB / "error_rate_schemes.csv", index=False
    )
    diffs = pd.concat(diffs, ignore_index=True)
    diffs.to_csv(TAB / "error_rate_paired_diff.csv", index=False)

    mixture(gex)
    mixture_error_nk(gex, js["blood"]).to_csv(
        TAB / "doublet_mixture_error_nk.csv", index=False
    )

    r = res[1]
    if r.empty:
        # niente figura di un'esecuzione precedente in results/
        (FIG / "error_rate_loco.png").unlink(missing_ok=True)
    else:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        bar_with_calibration(ax, r.tissue, r.error, r.lo, r.hi, r.n_pat)
        ax.set_ylabel("tasso d'errore dell'etichetta")
        ax.set_title("LOCO: riferimento disgiunto dal compartimento giudicato")
        fig.tight_layout()
        fig.savefig(FIG / "error_rate_loco.png", dpi=140)
        plt.close(fig)

    # Lettura condizionata: la differenza Tumor - Adjacent_normal regge in tutti gli schemi?
    print("\n" + "=" * 78)
    ta = diffs[diffs.contrast == "Tumor - Adjacent_normal"]
    for row in ta.itertuples():
        v = ci_verdict(row.lo, row.hi, row.n_patients)
        print(
            f"  {row.scheme:6s} Tumor - Adjacent_normal = {row.diff:+.3f} "
            f"[{row.lo:+.3f}, {row.hi:+.3f}]  {v}"
        )
    if ta.empty:
        print("Nessuno schema ha abbastanza cloni appaiati Tumor / Adjacent_normal.")
    elif (ta.n_patients <= SMALL_K).all():
        print(f"Tutti i contrasti hanno <= {SMALL_K} pazienti: l'IC non e' calibrato")
        print("(copertura sotto il 95% in simulazione), quindi nessun verdetto; vedi")
        print("n_distinct_resamples e frac_boot_gt0 in error_rate_paired_diff.csv.")
        if (ta["diff"] > 0).all():
            print("La stima puntuale e' positiva in tutti gli schemi.")
        else:
            print("La stima puntuale non ha lo stesso segno in tutti gli schemi.")
    else:
        # verdetto solo sugli schemi con IC calibrato (piu' di SMALL_K pazienti)
        cal = ta[ta.n_patients > SMALL_K]
        print(f"Schemi con piu' di {SMALL_K} pazienti: {', '.join(cal.scheme)}.")
        if (cal.lo > 0).all():
            print(
                "In questi schemi la differenza Tumor - Adjacent_normal ha IC sopra 0."
            )
        else:
            print("In almeno uno di questi schemi l'IC della differenza include 0:")
            print("l'asimmetria dipende dallo schema di riferimento.")
    print("=" * 78)


if __name__ == "__main__":
    main()

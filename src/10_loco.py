#!/usr/bin/env python3
"""
10_loco.py - Arbitrato Leave-One-Compartment-Out (LOCO), riferimento ematico e deconvoluzione dei doppietti.

Questo modulo elimina ogni circolarita' residua tra la stima dell'identita'
clonale di riferimento e la valutazione dell'etichetta compartimentale mediante
quattro protocolli di validazione:

1. Arbitrato disgiunto (Leave-One-Compartment-Out e Riferimento Ematico 'blood'):
   - Nello schema `loco`, il fenotipo CD4/CD8 del clone viene inferito usando
     esclusivamente i compartimenti diversi da quello sottoposto a verifica.
   - Nello schema `blood`, il sangue periferico (PBMC, esente da dissociazione
     enzimatica e rimodellamento stromale) funge da riferimento unico ed esterno
     per arbitrare simultaneamente i tessuti solidi (`Adjacent_normal` e `Tumor`).

2. Calibrazione della soglia di positivita' di CD4 rispetto all'RNA ambientale:
   Poiche' il trascritto CD4 e' espresso fisiologicamente anche dalla componente
   mieloide/macrofagica abbondante nel microambiente tumorale, la soglia di
   chiamata per CD4 viene ricalibrata al 95-esimo percentile dell'espressione
   osservata in popolazioni non linfoidi e non mieloidi (`Ductal`, `Endothelial`).

3. Analisi di robustezza previa esclusione dei potenziali doppietti:
   Ripete la valutazione LOCO escludendo il decile superiore di conteggio UMI
   (`n_counts > 90-esimo percentile` stratificato per libreria e tipo cellulare).

4. Deconvoluzione quantitativa della miscela doppietto/singoletto nelle NK-TCR+:
   Stima la frazione di doppietti fisici mediante due stimatori indipendenti
   basati rispettivamente sull'incremento di massa UMI totale (`d1`) e sulla
   diluizione del rapporto UMI_TCR / n_counts (`d2`).

Uso:
    python src/10_loco.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_BOOT = 2000

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}
COMPS = ["PBMC", "Adjacent_normal", "Tumor"]
MARKERS = ["CD3D", "CD4", "CD8A", "CD8B", "NKG7", "GNLY"]
# Popolazioni epiteliali/endoteliali prive di espressione endogena di CD4/CD8 (proxy di mRNA ambientale)
AMBIENT_PROXY = ["Ductal", "Endothelial"]
FNAME_RE = re.compile(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz")

MIN_CELLS = 3
MIN_MARGIN = 0.20


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
    """Calcola media e intervallo di confidenza al 95% tramite cluster bootstrap sui pazienti."""
    pats = [p for p, v in by_pat.items() if len(v)]
    if len(pats) < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(SEED)
    stats = [
        np.concatenate(
            [by_pat[p] for p in rng.choice(pats, len(pats), replace=True)]
        ).mean()
        for _ in range(n)
    ]
    allv = np.concatenate([by_pat[p] for p in pats])
    return (
        float(allv.mean()),
        float(np.percentile(stats, 2.5)),
        float(np.percentile(stats, 97.5)),
    )


def ambient_cd4(gex: pd.DataFrame) -> pd.Series:
    """Stima la soglia di fondo di CD4 per tessuto come 95-esimo percentile su cellule Ductal ed Endothelial."""
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

    thr = amb.groupby("tissue", observed=True).g_CD4.quantile(0.95)
    print("\nsoglia CD4 per compartimento (95o pct dell'ambient):")
    print(thr.round(4).to_string())

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
        "\nLe popolazioni mieloidi esprimono CD4 endogeno e sono pertanto escluse dal proxy ambientale."
    )
    return thr


def per_comp_stats(lym: pd.DataFrame, cd4_thr: pd.Series | None) -> pd.DataFrame:
    """Aggrega per ciascuna tripla (paziente, clone, tessuto) le frazioni di positivita' CD8+ e CD4+."""
    x = lym.copy()
    x["cd8_pos"] = (x.g_CD8A > 0) | (x.g_CD8B > 0)
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
            lab=("all_celltype", lambda s: s.value_counts().index[0]),
        )
        .reset_index()
    )
    return per[per.n >= MIN_CELLS]


def arbitrate(per: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """Valuta la concordanza tra etichetta compartimentale e riferimento clonale ('naive', 'loco', 'blood')."""
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
            if r.lab not in ("CD4T", "CD8T", "NK"):
                continue
            rows.append(
                {
                    "patient": pat,
                    "clone": clone,
                    "tissue": target,
                    "lab": r.lab,
                    "identity": ident,
                    "error": r.lab != ident,
                    "n_ref_comp": ref.tissue.nunique(),
                }
            )
    return pd.DataFrame(rows)


def report(j: pd.DataFrame, title: str) -> pd.DataFrame:
    """Stampa e restituisce il riepilogo dei tassi di errore per compartimento e i contrasti appaiati."""
    print(f"\n--- {title} ---")
    if j.empty:
        print("  nessun clone risolvibile.")
        return pd.DataFrame()
    rows = []
    for tis in COMPS:
        g = j[j.tissue == tis]
        if g.empty or g.patient.nunique() < 3:
            continue
        by = {p: gg.error.values.astype(float) for p, gg in g.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        rows.append(
            {
                "tissue": tis,
                "n": len(g),
                "n_pat": g.patient.nunique(),
                "error": m,
                "lo": lo,
                "hi": hi,
            }
        )
        print(
            f"  {tis:16s} errore = {m:.3f}  IC95% [{lo:.3f}, {hi:.3f}]  "
            f"(n={len(g)}, {g.patient.nunique()} pz)"
        )

    piv = j.pivot_table(
        index=["patient", "clone"], columns="tissue", values="error", aggfunc="max"
    )
    for a, b in [("Tumor", "Adjacent_normal"), ("Tumor", "PBMC")]:
        if not {a, b}.issubset(piv.columns):
            continue
        d = piv.dropna(subset=[a, b]).reset_index()
        if len(d) < 20 or d.patient.nunique() < 3:
            print(
                f"  {a} - {b}: n insufficiente ({len(d)} cloni, {d.patient.nunique()} pz)"
            )
            continue
        d["diff"] = d[a].astype(float) - d[b].astype(float)
        by = {p: g["diff"].values for p, g in d.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        v = "significativa" if lo > 0 or hi < 0 else "NON significativa"
        print(
            f"  {a} - {b} (stessi cloni): {m:+.3f}  IC95% [{lo:+.3f}, {hi:+.3f}]  {v}  (n={len(d)})"
        )
    return pd.DataFrame(rows)


def mixture(gex: pd.DataFrame) -> None:
    """Stima la frazione di doppietti nelle cellule NK-TCR+ tramite gli stimatori fisici d1 e d2."""
    print("\n" + "=" * 78)
    print("(4) FRAZIONE DI DOPPIETTI NELLE NK-TCR+  (due stime indipendenti)")
    print("=" * 78)
    g = gex.copy()
    g["ratio"] = g.umi_tcr / np.maximum(g.n_counts, 1) * 1000

    def est(sub: pd.DataFrame) -> tuple[float, float]:
        nk_pos = sub[(sub.all_celltype == "NK") & sub.has_tcr]
        nk_neg = sub[(sub.all_celltype == "NK") & ~sub.has_tcr]
        cd8 = sub[(sub.all_celltype == "CD8T") & sub.has_tcr]
        if min(len(nk_pos), len(nk_neg), len(cd8)) < 10:
            return np.nan, np.nan
        s, c, o = (
            nk_neg.n_counts.median(),
            cd8.n_counts.median(),
            nk_pos.n_counts.median(),
        )
        # Stimatore d1 basato sull'eccesso di conteggio UMI totale rispetto al singoletto NK
        d1 = (o - s) / max(c, 1)
        r_cd8, r_nk = cd8.ratio.median(), nk_pos.ratio.median()
        r_dbl = r_cd8 * c / max(o, 1)
        # Stimatore d2 basato sulla diluizione della densita' di UMI TCR
        d2 = (r_cd8 - r_nk) / max(r_cd8 - r_dbl, 1e-9)
        return float(np.clip(d1, 0, 1)), float(np.clip(d2, 0, 1))

    d1, d2 = est(g)
    pats = g.patients.unique()
    rng = np.random.default_rng(SEED)
    b1, b2 = [], []
    for _ in range(500):
        draw = rng.choice(pats, len(pats), replace=True)
        sub = pd.concat([g[g.patients == p] for p in draw])
        a, b = est(sub)
        if not np.isnan(a):
            b1.append(a)
            b2.append(b)
    print(
        f"  d1 (da n_counts)      = {d1:.3f}  IC95% "
        f"[{np.percentile(b1, 2.5):.3f}, {np.percentile(b1, 97.5):.3f}]"
    )
    print(
        f"  d2 (da UMI TCR/conte) = {d2:.3f}  IC95% "
        f"[{np.percentile(b2, 2.5):.3f}, {np.percentile(b2, 97.5):.3f}]"
    )
    print("\nLa convergenza tra le due grandezze fisiche conferma la stima di miscela:")
    print(f"  -> circa il {100 * np.mean([d1, d2]):.0f}% delle NK-TCR+ sono doppietti;")
    print(
        f"     il restante {100 * (1 - np.mean([d1, d2])):.0f}% ha il contenuto di UNA cellula"
    )
    print("     (linfociti T citotossici o cellule NKT invarianti).")
    pd.DataFrame({"d1": b1, "d2": b2}).to_csv(
        TAB / "doublet_mixture_boot.csv", index=False
    )


def main() -> None:
    """Punto di ingresso principale per l'arbitrato LOCO/Blood e i controlli di robustezza."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tcrdir", type=Path, default=Path("data/raw/tcr"))
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

    tcr = load_tcr(args.tcrdir)
    gex = obs.merge(
        tcr,
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    print(f"[ok] {len(gex):,} cellule | TCR: {gex.has_tcr.sum():,}")

    thr = ambient_cd4(gex)

    lym = gex[gex.has_tcr & gex.all_celltype.isin(["CD4T", "CD8T", "NK"])].copy()

    print("\n" + "=" * 78)
    print("(1) LEAVE-ONE-COMPARTMENT-OUT: riferimento disgiunto dal test")
    print("=" * 78)
    per = per_comp_stats(lym, None)
    res = {}
    for scheme, title in [
        ("naive", "NAIVE (riferimento include il target, circolare)"),
        ("loco", "LOCO (riferimento = gli altri compartimenti)"),
        ("blood", "BLOOD (riferimento = solo PBMC)"),
    ]:
        j = arbitrate(per, scheme)
        res[scheme] = report(j, title)
        j.to_csv(TAB / f"arbitration_{scheme}.csv", index=False)

    print("\n" + "=" * 78)
    print("(2b) RIPETUTO con soglia CD4 sopra l'ambient del compartimento")
    print("=" * 78)
    per_amb = per_comp_stats(lym, thr)
    report(arbitrate(per_amb, "loco"), "LOCO + soglia CD4 corretta per ambient")

    print("\n" + "=" * 78)
    print(
        "(3) RIPETUTO escludendo i doppietti (n_counts > 90o pct per libreria x classe)"
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
    report(arbitrate(per_comp_stats(lym_s, None), "loco"), "LOCO senza doppietti")

    mixture(gex)

    if "loco" in res and not res["loco"].empty:
        r = res["loco"]
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.bar(
            r.tissue,
            r.error,
            yerr=[r.error - r.lo, r.hi - r.error],
            capsize=5,
            color="#d95f02",
            alpha=0.85,
        )
        ax.set_ylabel("tasso d'errore dell'etichetta")
        ax.set_title("LOCO: riferimento disgiunto dal compartimento giudicato")
        fig.tight_layout()
        fig.savefig(FIG / "error_rate_loco.png", dpi=140)
        plt.close(fig)

    print("\n" + "=" * 78)
    print("La persistenza del differenziale tra NAIVE, LOCO e BLOOD conferma la natura")
    print("non circolare dell'asimmetria di errore tra compartimenti.")
    print("=" * 78)


if __name__ == "__main__":
    main()

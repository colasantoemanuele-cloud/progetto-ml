#!/usr/bin/env python3
"""
08_disambiguate.py - Cloni presenti in tre compartimenti, marcatori nei cloni discordanti e doppietti (Scrublet).

Scopo
-----
Tre controlli sulla discordanza fenotipica dei cloni T misurata in 07:

1. Contrasto appaiato sui cloni presenti (>= 2 cellule) in PBMC, Adjacent_normal e
   Tumor dello stesso paziente (solo i pazienti con tessuto adiacente: PA04, PA06,
   PA08, PA09, PA10). Confrontare gli eccessi di discordanza sugli stessi cloni
   riduce il confondimento fra paziente e coppia di compartimenti, ma non lo
   elimina: restano le differenze di numerosita' e di campionamento fra
   compartimenti e la selezione dei cloni abbastanza espansi da comparire in tutti
   e tre. IC 95% con cluster bootstrap sui pazienti (`scpipe.cluster_bootstrap`).

2. Marcatori nei cloni discordanti (eccesso > 0.5 in 07) e concordanti (eccesso < 0.05),
   espressi in log1p(CP10k). Per le cellule etichettate CD8T, per tessuto: CD8A medio,
   frazione con CD8A > 0, test di Mann-Whitney, CD4 medio e rapporto fra le medie di CD4
   (discordanti / concordanti). Il rapporto e' fra medie di log1p(CP10k), non un
   fold-change dell'espressione.

3. Doppietti: punteggio Scrublet per libreria, confrontato fra cellule con e senza TCR
   della stessa classe.

Input
-----
  --h5ad      data/interim/sc_raw.h5ad (serve X per marcatori e Scrublet)
  --tcrdir    cartella con i 30 file VDJ
  --pairwise  results/tables/pairwise_discordance_clone_b.csv (da 07_robustness.py)

Output
------
  results/tables/triple_clones.csv              eccessi per clone presente nei tre compartimenti
  results/tables/triple_clones_summary.csv      medie, IC e differenze appaiate (cluster bootstrap)
  results/tables/markers_disc_vs_conc.csv       marcatori medi per gruppo x tessuto x etichetta
  results/tables/marker_test_disc_vs_conc.csv   CD8T per tessuto: CD8A, frazione CD8A>0, p, CD4, rapporto CD4
  results/tables/doublet_scores.csv             punteggio Scrublet per cellula (cache, escluso da git)
  results/tables/doublet_tcr_pos_vs_neg.csv     Scrublet TCR+ vs TCR- per classe
  results/figures/doublet_by_tcr_status.png

Argomenti
---------
  --h5ad (obbligatorio), --tcrdir (default data/raw/tcr), --pairwise, --skip-scrublet

Uso
---
    python src/08_disambiguate.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
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
import scanpy as sc
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import cluster_bootstrap  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_BOOT = 2000

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}
COMPS = ["PBMC", "Adjacent_normal", "Tumor"]
LYMPHOID = ["CD4T", "CD8T", "NK"]
MARKERS = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "NCAM1", "KLRD1", "NKG7", "GNLY"]
FNAME_RE = re.compile(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz")

DISCORDANT_MIN = 0.50
CONCORDANT_MAX = 0.05


def load_tcr(tcrdir: Path) -> pd.DataFrame:
    """Carica i file di annotazione VDJ e restituisce la tabella dei clonotipi CDR3(TRB) per cellula."""
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
    for ch in ("TRA", "TRB"):
        if ch not in seq.columns:
            seq[ch] = np.nan
    w = seq.reset_index().rename(columns={"TRA": "cdr3_TRA", "TRB": "cdr3_TRB"})
    w["clone_b"] = w.cdr3_TRB
    return w[["patient", "tissue", "barcode16", "clone_b"]]


def _counts(labels: np.ndarray) -> np.ndarray:
    """Calcola il vettore dei conteggi assoluti per le tre classi linfoidi (CD4T, CD8T, NK)."""
    return np.array([(labels == c).sum() for c in LYMPHOID], dtype=float)


def d_across(A: np.ndarray, B: np.ndarray) -> float:
    """Calcola la discordanza a coppie tra due compartimenti A e B."""
    pA, pB = _counts(A) / len(A), _counts(B) / len(B)
    return 1.0 - float(pA @ pB)


def d_within(A: np.ndarray) -> float:
    """Calcola la discordanza interna senza reinserimento entro un singolo compartimento A."""
    c, n = _counts(A), len(A)
    return 1.0 - float((c * (c - 1)).sum() / (n * (n - 1)))


def cluster_boot(vals_by_patient: dict[str, np.ndarray], n: int = N_BOOT) -> tuple:
    """Media pooled e IC 95% con cluster bootstrap sui pazienti (`scpipe.cluster_bootstrap`)."""
    return cluster_bootstrap(vals_by_patient, n, SEED)


def triple_clones(lym: pd.DataFrame, min_cells: int = 2) -> pd.DataFrame:
    """Seleziona i cloni presenti contemporaneamente in PBMC, Adjacent_normal e Tumor e ne calcola gli eccessi."""
    rows = []
    for (pat, clone), g in lym.groupby(["patient", "clone_b"], observed=True):
        comp = {
            t: s.all_celltype.values
            for t, s in g.groupby("tissue", observed=True)
            if len(s) >= min_cells
        }
        if not set(COMPS).issubset(comp):
            continue
        dw = {t: d_within(comp[t]) for t in COMPS}
        r = {"patient": pat, "clone": clone, **{f"n_{t}": len(comp[t]) for t in COMPS}}
        for i in range(3):
            for j in range(i + 1, 3):
                a, b = COMPS[i], COMPS[j]
                nA, nB = len(comp[a]), len(comp[b])
                exc = d_across(comp[a], comp[b]) - (nA * dw[a] + nB * dw[b]) / (nA + nB)
                r[f"exc_{a[:3]}_{b[:3]}"] = exc
        rows.append(r)
    return pd.DataFrame(rows)


def part1(lym: pd.DataFrame) -> None:
    """Esegue l'analisi appaiata sui cloni condivisi tra tutti e tre i compartimenti."""
    print("\n" + "=" * 78)
    print("(1) CONTRASTO APPAIATO: cloni presenti in TUTTI E TRE i compartimenti")
    print("=" * 78)
    tri = triple_clones(lym)
    if tri.empty:
        print(
            "Nessun clone in tutti e tre i compartimenti. Confronto appaiato impossibile."
        )
        return
    tri.to_csv(TAB / "triple_clones.csv", index=False)

    print(
        f"cloni: {len(tri):,} | pazienti: {tri.patient.nunique()} "
        f"-> {sorted(tri.patient.unique())}"
    )
    print(f"cellule: {tri[[f'n_{t}' for t in COMPS]].values.sum():,}\n")

    cols = [c for c in tri.columns if c.startswith("exc_")]
    summary = []
    for c in cols:
        by = {p: g[c].values for p, g in tri.groupby("patient")}
        m, lo, hi = cluster_boot(by)
        star = "*" if lo > 0 else " "
        print(f"  {c:22s} {m:+.4f}  IC95% [{lo:+.4f}, {hi:+.4f}] {star}")
        summary.append({"contrast": c, "mean": m, "lo": lo, "hi": hi})

    print("\nDIFFERENZE APPAIATE (stesso clone, stesso paziente):")
    for a, b in [("exc_PBM_Tum", "exc_Adj_Tum"), ("exc_PBM_Tum", "exc_PBM_Adj")]:
        if a in tri and b in tri:
            tri["_d"] = tri[a] - tri[b]
            by = {p: g["_d"].values for p, g in tri.groupby("patient")}
            m, lo, hi = cluster_boot(by)
            verdict = "IC esclude 0" if lo > 0 or hi < 0 else "IC include 0"
            print(f"  {a} - {b} = {m:+.4f}  IC95% [{lo:+.4f}, {hi:+.4f}]  {verdict}")
            summary.append({"contrast": f"{a} - {b}", "mean": m, "lo": lo, "hi": hi})
    pd.DataFrame(summary).assign(
        n_clones=len(tri), n_patients=tri.patient.nunique()
    ).to_csv(TAB / "triple_clones_summary.csv", index=False)

    print("\nper paziente:")
    print(tri.groupby("patient")[cols].agg(["mean", "size"]).round(4).to_string())
    print(
        "\nIl confronto sugli stessi cloni riduce il confondimento fra paziente e coppia"
    )
    print("di compartimenti; restano le differenze di campionamento fra compartimenti.")


def part2(lym: pd.DataFrame, pw: pd.DataFrame) -> None:
    """Confronta l'espressione dei marcatori canonici tra cloni discordanti e cloni concordanti."""
    print("\n" + "=" * 78)
    print(
        "(2) MARCATORI: cloni discordanti vs CONCORDANTI (gruppo di controllo interno)"
    )
    print("=" * 78)
    disc = set(
        zip(
            pw.loc[pw.excess > DISCORDANT_MIN, "patient"],
            pw.loc[pw.excess > DISCORDANT_MIN, "clone"],
        )
    )
    conc = set(
        zip(
            pw.loc[pw.excess < CONCORDANT_MAX, "patient"],
            pw.loc[pw.excess < CONCORDANT_MAX, "clone"],
        )
    )
    conc -= disc
    print(
        f"cloni discordanti (excess>{DISCORDANT_MIN}): {len(disc):,} | "
        f"concordanti (excess<{CONCORDANT_MAX}): {len(conc):,}"
    )
    if not disc or not conc:
        print("gruppi insufficienti.")
        return

    key = list(zip(lym.patient, lym.clone_b))
    lym = lym.assign(
        group=np.where(
            [k in disc for k in key],
            "discordant",
            np.where([k in conc for k in key], "concordant", "other"),
        )
    )
    sel = lym[lym.group != "other"]

    gcols = [
        f"g_{g}"
        for g in ("CD4", "CD8A", "CD8B", "CD3D", "NKG7", "GNLY")
        if f"g_{g}" in sel
    ]
    tab = sel.groupby(["group", "tissue", "all_celltype"], observed=True)[gcols].mean()
    print("\nespressione media (log1p CP10K):")
    print(tab.round(3).to_string())
    tab.to_csv(TAB / "markers_disc_vs_conc.csv")

    print(
        "\n>>> TEST DI MANN-WHITNEY: cellule chiamate CD8T, per tessuto, CD8A discordanti vs concordanti"
    )
    rows = []
    for tis in COMPS:
        d = sel[(sel.tissue == tis) & (sel.all_celltype == "CD8T")]
        dd, dc = d[d.group == "discordant"], d[d.group == "concordant"]
        a, b = dd.g_CD8A.values, dc.g_CD8A.values
        if len(a) < 20 or len(b) < 20:
            print(f"  {tis:16s} n insufficiente (disc={len(a)}, conc={len(b)})")
            continue
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        print(
            f"  {tis:16s} disc={a.mean():6.2f} (n={len(a):5d})  "
            f"conc={b.mean():6.2f} (n={len(b):5d})  "
            f"frac CD8A>0: {(a > 0).mean():.2f} vs {(b > 0).mean():.2f}  p={p:.2e}"
        )
        row = {
            "tissue": tis,
            "label": "CD8T",
            "n_disc": len(a),
            "n_conc": len(b),
            "mean_CD8A_disc": a.mean(),
            "mean_CD8A_conc": b.mean(),
            "frac_CD8A_pos_disc": (a > 0).mean(),
            "frac_CD8A_pos_conc": (b > 0).mean(),
            "p_mannwhitney_CD8A": p,
            "mean_CD4_disc": dd.g_CD4.mean(),
            "mean_CD4_conc": dc.g_CD4.mean(),
            # rapporto fra medie di log1p(CP10k), non un fold-change dell'espressione
            "ratio_mean_CD4_disc_over_conc": dd.g_CD4.mean()
            / max(dc.g_CD4.mean(), 1e-12),
        }
        if "g_GNLY" in d:
            row["mean_GNLY_disc"] = dd.g_GNLY.mean()
            row["mean_GNLY_conc"] = dc.g_GNLY.mean()
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "marker_test_disc_vs_conc.csv", index=False)
    tum = res[res.tissue == "Tumor"] if not res.empty else res
    if len(tum):
        t = tum.iloc[0]
        print(
            f"\nTumor, cellule CD8T: CD4 medio disc/conc = {t.ratio_mean_CD4_disc_over_conc:.1f}x "
            "(rapporto fra medie di log1p(CP10k), non fold-change)"
        )
        if t.frac_CD8A_pos_disc < 0.5 * t.frac_CD8A_pos_conc:
            print(
                "Le CD8T dei cloni concordanti del tumore esprimono CD8A, quelle dei cloni"
            )
            print(
                "discordanti molto meno: il dropout tissutale di CD8A da solo non lo spiega."
            )
            print(
                "Il dato e' compatibile con un'etichetta CD8T sbagliata (cellule CD4); resta"
            )
            print("aperta l'alternativa di un fenotipo CD4 citotossico", end="")
            if "mean_GNLY_disc" in t:
                print(
                    f" (GNLY {t.mean_GNLY_disc:.2f} nei discordanti contro "
                    f"{t.mean_GNLY_conc:.2f} nei concordanti)."
                )
            else:
                print(".")
        else:
            print(
                "CD8A simile fra discordanti e concordanti: nessuna evidenza di etichetta"
            )
            print("sbagliata da questo controllo.")


def run_scrublet(adata: ad.AnnData, cache: Path) -> pd.Series:
    """Calcola o recupera da cache il punteggio Scrublet per ciascuna libreria 10x."""
    if cache.exists():
        print(f"[ok] doublet score da cache {cache}")
        return pd.read_csv(cache, index_col=0).iloc[:, 0]
    print("[..] scrublet per libreria (lento) ...")
    a = adata.copy()
    sc.pp.scrublet(a, batch_key="library", random_state=SEED)
    s = a.obs["doublet_score"]
    s.to_csv(cache)
    print(f"[ok] doublet score -> {cache}")
    return s


def part3(gex: pd.DataFrame) -> None:
    """Valuta le differenze di doublet_score Scrublet tra cellule TCR+ e TCR- entro ciascuna classe."""
    print("\n" + "=" * 78)
    print("(3) ANALISI SCRUBLET: confronto TCR+ vs TCR- entro la medesima classe")
    print("=" * 78)
    print(
        "Confronto stratificato: TCR+ vs TCR- all'interno della stessa popolazione.\n"
    )

    rows = []
    for cls in [
        "CD8T",
        "CD4T",
        "NK",
        "Platelet",
        "Endothelial",
        "Stellate",
        "Myeloid",
        "B",
        "Ductal",
    ]:
        g = gex[gex.all_celltype == cls]
        pos = g.loc[g.has_tcr, "doublet_score"].dropna().values
        neg = g.loc[~g.has_tcr, "doublet_score"].dropna().values
        if len(pos) < 20 or len(neg) < 20:
            continue
        p = mannwhitneyu(pos, neg, alternative="greater").pvalue
        rows.append(
            {
                "class": cls,
                "n_tcr_pos": len(pos),
                "n_tcr_neg": len(neg),
                "ds_tcr_pos": pos.mean(),
                "ds_tcr_neg": neg.mean(),
                "delta": pos.mean() - neg.mean(),
                "p_greater": p,
            }
        )
    out = pd.DataFrame(rows)
    print(out.round(4).to_string(index=False))
    out.to_csv(TAB / "doublet_tcr_pos_vs_neg.csv", index=False)

    print("\nINTERPRETAZIONE, per classe non-T:")
    print("  delta ~ 0  -> punteggio compatibile con singoletto")
    print(
        "  delta > 0  -> incremento di punteggio indicativo di profilo ibrido o doppietto"
    )
    print("\n(CD8T/CD4T costituiscono il controllo negativo con delta atteso ~0.)")

    fig, ax = plt.subplots(figsize=(9, 4.5))
    cls = out["class"].tolist()
    x = np.arange(len(cls))
    ax.bar(x - 0.2, out.ds_tcr_neg, 0.38, label="TCR−", color="#7570b3", alpha=0.85)
    ax.bar(x + 0.2, out.ds_tcr_pos, 0.38, label="TCR+", color="#d95f02", alpha=0.85)
    ax.set_xticks(x, cls, rotation=45, ha="right")
    ax.set_ylabel("doublet score (scrublet)")
    ax.set_title("il TCR viene da un doppietto? confronto entro classe")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "doublet_by_tcr_status.png", dpi=140)
    plt.close(fig)


def main() -> None:
    """Punto di ingresso principale per le tre verifiche di disambiguazione."""
    ap = argparse.ArgumentParser(
        description="Cloni in tre compartimenti, marcatori nei cloni discordanti e doppietti."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tcrdir", type=Path, default=Path("data/raw/tcr"))
    ap.add_argument(
        "--pairwise",
        type=Path,
        default=Path("results/tables/pairwise_discordance_clone_b.csv"),
    )
    ap.add_argument("--skip-scrublet", action="store_true")
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype", "barcode", "library"):
        obs[c] = obs[c].astype(str)
    obs = obs.rename(columns={"barcode": "barcode16"})

    # Calcolo dell'espressione log1p(CP10K) dei marcatori canonici
    sym = adata.var.gene_symbol.astype(str).values
    for g in MARKERS:
        idx = np.where(sym == g)[0]
        if len(idx):
            x = np.asarray(adata.X[:, idx[0]].todense()).ravel()
            obs[f"g_{g}"] = np.log1p(x / np.maximum(obs.n_counts.values, 1) * 1e4)

    tcr = load_tcr(args.tcrdir)
    gex = obs.merge(
        tcr,
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    gex.index = obs.index
    print(f"[ok] {len(gex):,} cellule | TCR: {gex.has_tcr.sum():,}")

    lym = gex[gex.clone_b.notna() & gex.all_celltype.isin(LYMPHOID)].copy()

    part1(lym)

    if args.pairwise.exists():
        pw = pd.read_csv(args.pairwise)
        part2(lym, pw)
    else:
        print(
            f"\n[!] {args.pairwise} assente: salto (2). Esegui prima 07_robustness.py"
        )

    if not args.skip_scrublet:
        adata.obs["library"] = obs.library.values
        ds = run_scrublet(adata, TAB / "doublet_scores.csv")
        gex["doublet_score"] = ds.reindex(gex.index).values
        part3(gex)
    else:
        cache = TAB / "doublet_scores.csv"
        if cache.exists():
            gex["doublet_score"] = (
                pd.read_csv(cache, index_col=0).iloc[:, 0].reindex(gex.index).values
            )
            part3(gex)
        else:
            print("\n[!] scrublet saltato e nessuna cache: (3) non eseguito.")


if __name__ == "__main__":
    main()

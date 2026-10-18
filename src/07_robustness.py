#!/usr/bin/env python3
"""
07_robustness.py - Discordanza clonale con stimatore a coppie non distorto, cluster bootstrap e controlli.

Scopo
-----
1. Stimatore a coppie della discordanza (sostituisce lo split-half di 06). Per ogni
   clone presente con almeno 2 cellule in due compartimenti A e B:
       D_AxB    = 1 - sum_k p_A(k) * p_B(k)                     (una cellula da A, una da B)
       D_A      = 1 - sum_k n_k (n_k - 1) / [n_A (n_A - 1)]     (due cellule da A, senza reinserimento)
       D_within = (n_A * D_A + n_B * D_B) / (n_A + n_B)
       eccesso  = D_AxB - D_within
   D_A e' non distorto per 1 - sum_k p_k^2 (verificato in tests/test_stats.py), quindi
   l'eccesso ha valore atteso zero se le etichette hanno la stessa distribuzione nei
   due compartimenti. Analisi di sensibilita': eccesso a pesi uguali
   1/2 * ||p_A - p_B||^2 = D_AxB - (G_A + G_B) / 2, con G l'indice di Gini con
   reinserimento (plug-in, distorto per cloni piccoli).
2. IC 95% con cluster bootstrap sui pazienti (`scpipe.cluster_bootstrap`, 2.000
   repliche, intervallo "expanded percentile"): le coppie dello stesso paziente non
   sono indipendenti. Per confronto si riporta anche il bootstrap ingenuo sulle
   coppie (unita' = coppia di compartimenti), che sottostima l'incertezza.
3. RNA ambientale e stima c_hat della contaminazione T nelle NK: confronto di
   UMI(TRB) e frazione TRA+TRB appaiate fra classi; sensibilita' di c_hat alla classe
   non-T usata come fondo q, con tutti i TCR e con quelli ad alta confidenza
   (TRA+TRB appaiate e UMI_TRB >= --umi-min).
4. Marcatori (CD4, CD8A, CD8B, CD3D, NKG7, GNLY, in CP10k senza log) nelle cellule dei
   cloni fortemente discordanti (eccesso > 0.5).

Input
-----
  --h5ad    data/interim/sc_raw.h5ad (X necessaria per i marcatori, salvo --no-markers)
  --tcrdir  cartella con i 30 file *_filtered_contig_annotations.csv.gz

Output ({cc} = definizione di clone, clone_b o clone_ab)
------
  results/tables/pairwise_discordance_{cc}.csv      una riga per (clone, coppia di compartimenti)
  results/tables/excess_summary_{cc}.csv            D_across, D_within, eccesso e versione a pesi
                                                    uguali, con IC cluster e IC ingenuo
  results/tables/excess_by_compartment_pair_{cc}.csv eccesso per coppia di compartimenti; `significativo`
                                                    (lo > 0) solo se `ic_calibrato` (> 5 pazienti), altrimenti NaN
  results/tables/discordant_clone_markers_{cc}.csv  marcatori nei cloni con eccesso > 0.5
  results/tables/ambient_discriminant.csv           UMI TRB e appaiamento per classe, pazienti con log(UMI)
                                                    piu' alto di CD8T / del fondo Ductal+B (descrittivo,
                                                    nessun test; non dipende da {cc})
  results/tables/chat_sensitivity_raw.csv, chat_sensitivity_hc.csv (non dipendono da {cc})
  results/figures/ambient_vs_real_tcr.png, results/figures/excess_by_pair_{cc}.png
                                                    (IC con <= 5 pazienti in grigio tratteggiato:
                                                    non calibrati)
  Se una sezione non ha dati, le sue tabelle e figure non vengono scritte e quelle di
  un'esecuzione precedente vengono cancellate.

Argomenti
---------
  --h5ad (obbligatorio), --tcrdir (default data/raw/tcr),
  --clone-def {clone_b,clone_ab} (default clone_b), --umi-min (default 3), --no-markers

Uso
---
    python src/07_robustness.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
    python src/07_robustness.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr --clone-def clone_ab
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
from scpipe import SMALL_K  # noqa: E402
from scpipe import cluster_bootstrap as _cluster_bootstrap  # noqa: E402

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_BOOT = 2000
# Cellule TCR+ minime per (paziente, classe) nei confronti descrittivi per paziente
MIN_CELLS_PATIENT = 5

TISSUE_MAP = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}
LYMPHOID = ["CD4T", "CD8T", "NK"]
FNAME_RE = re.compile(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz")

# Classi non-T candidate per la stima del tasso di fondo q (valutazione di sensibilita')
Q_REFS = ["B", "Ductal", "Myeloid", "CAF", "Stellate", "Endothelial", "Plasma"]
MARKERS = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "NCAM1", "KLRD1", "NKG7", "GNLY"]


def load_tcr(tcrdir: Path) -> pd.DataFrame:
    """Estrae per ciascuna cellula le sequenze CDR3(TRA/TRB), i conteggi UMI e lo stato di appaiamento."""
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
    # Selezione del contig dominante per (cellula, catena) in base al supporto UMI
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
        seq.rename(columns={"TRA": "cdr3_TRA", "TRB": "cdr3_TRB"})
        .join(umi.rename(columns={"TRA": "umi_TRA", "TRB": "umi_TRB"}))
        .reset_index()
    )

    w["paired"] = w.cdr3_TRA.notna() & w.cdr3_TRB.notna()
    w["clone_b"] = w.cdr3_TRB
    w["clone_ab"] = np.where(w.paired, w.cdr3_TRA + "|" + w.cdr3_TRB, np.nan)
    return w


def load_gex(h5ad: Path, with_markers: bool) -> pd.DataFrame:
    """Carica i metadati GEX e, se richiesto, l'espressione dei marcatori in CP10k (senza log)."""
    adata = ad.read_h5ad(h5ad, backed=None if with_markers else "r")
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype", "barcode"):
        obs[c] = obs[c].astype(str)
    obs = obs.rename(columns={"barcode": "barcode16"})

    if with_markers:
        sym = adata.var.gene_symbol.astype(str).values
        for g in MARKERS:
            idx = np.where(sym == g)[0]
            if len(idx) == 0:
                print(f"[!] marcatore assente: {g}")
                continue
            x = np.asarray(adata.X[:, idx[0]].todense()).ravel()
            # CP10k per cellula (conteggi / n_counts * 1e4), senza log1p: non e' CPM
            obs[f"g_{g}"] = x / np.maximum(obs.n_counts.values, 1) * 1e4
    return obs


def _counts(labels: np.ndarray, classes: list[str]) -> np.ndarray:
    """Restituisce il vettore delle frequenze assolute per le classi specificate."""
    return np.array([(labels == c).sum() for c in classes], dtype=float)


def pairwise_discordance(
    df: pd.DataFrame, clone_col: str, min_cells: int = 2
) -> pd.DataFrame:
    """Calcola analiticamente D_across, D_within ed eccesso di discordanza per ogni coppia (clone, compartimenti)."""
    rows = []
    for (pat, clone), g in df.groupby(["patient", clone_col], observed=True):
        comp = {
            t: s.all_celltype.values
            for t, s in g.groupby("tissue", observed=True)
            if len(s) >= min_cells
        }
        if len(comp) < 2:
            continue
        ts = sorted(comp)
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                A, B = comp[ts[i]], comp[ts[j]]
                nA, nB = len(A), len(B)
                cA, cB = _counts(A, LYMPHOID), _counts(B, LYMPHOID)
                pA, pB = cA / nA, cB / nB

                # Probabilita' di discordanza estraendo una cellula dal compartimento A e una da B
                d_across = 1.0 - float(pA @ pB)
                # Probabilita' di discordanza interna senza reinserimento entro ciascun compartimento
                d_A = 1.0 - float((cA * (cA - 1)).sum() / (nA * (nA - 1)))
                d_B = 1.0 - float((cB * (cB - 1)).sum() / (nB * (nB - 1)))
                d_within = (nA * d_A + nB * d_B) / (nA + nB)
                # Sensibilita': pesi uguali, 1/2 ||pA - pB||^2 (Gini plug-in con reinserimento)
                excess_eq = 0.5 * float(((pA - pB) ** 2).sum())

                rows.append(
                    {
                        "patient": pat,
                        "clone": clone,
                        "comp_a": ts[i],
                        "comp_b": ts[j],
                        "n_a": nA,
                        "n_b": nB,
                        "d_across": d_across,
                        "d_within": d_within,
                        "excess": d_across - d_within,
                        "excess_eq": excess_eq,
                    }
                )
    return pd.DataFrame(rows)


def cluster_bootstrap(df: pd.DataFrame, col: str = "excess", n: int = N_BOOT) -> tuple:
    """Media di `col` e IC 95% con cluster bootstrap sui pazienti (`scpipe.cluster_bootstrap`).

    I pazienti entrano nell'ordine di prima comparsa in `df`; stesso seme per ogni
    colonna, quindi `excess` ed `excess_eq` usano le stesse estrazioni di pazienti.
    """
    by_pat = {p: df.loc[df.patient == p, col].values for p in df.patient.unique()}
    return _cluster_bootstrap(by_pat, n, SEED)


def naive_bootstrap(df: pd.DataFrame, col: str = "excess", n: int = N_BOOT) -> tuple:
    """IC 95% percentile ricampionando le singole coppie (unita' sbagliata, solo per confronto)."""
    rng = np.random.default_rng(SEED)
    v = df[col].values
    stats = np.array([rng.choice(v, len(v), replace=True).mean() for _ in range(n)])
    return (
        float(v.mean()),
        float(np.percentile(stats, 2.5)),
        float(np.percentile(stats, 97.5)),
    )


def ambient_discriminant(gex: pd.DataFrame) -> pd.DataFrame:
    """Confronta la distribuzione di UMI(TRB) e la frazione TRA+TRB appaiata rispetto a CD8T e controlli non-T."""
    t = gex[gex.has_tcr]
    out = (
        t.groupby("all_celltype", observed=True)
        .agg(
            n=("paired", "size"),
            frac_paired=("paired", "mean"),
            median_umi_TRB=("umi_TRB", "median"),
            mean_log_umi=("log_umi", "mean"),
        )
        .sort_values("n", ascending=False)
    )
    # Confronto descrittivo per paziente (la cellula non e' l'unita' indipendente): per
    # ogni classe, pazienti con >= MIN_CELLS_PATIENT cellule TCR+ nella classe e nel
    # riferimento, e quanti di loro hanno log(UMI) medio piu' alto del riferimento.
    pm = t.groupby(["patients", "all_celltype"], observed=True).log_umi.agg(
        ["mean", "size"]
    )
    pm = pm[pm["size"] >= MIN_CELLS_PATIENT]["mean"].unstack()
    refs = {"vs_CD8T": ["CD8T"], "vs_ambient": ["Ductal", "B"]}
    for name, ref_cls in refs.items():
        rm = t[t.all_celltype.isin(ref_cls)].groupby("patients", observed=True).log_umi
        rm = rm.agg(["mean", "size"])
        rm = rm[rm["size"] >= MIN_CELLS_PATIENT]["mean"]
        n_eval, n_higher = {}, {}
        for cls in out.index:
            if cls in ref_cls or cls not in pm.columns:
                n_eval[cls], n_higher[cls] = np.nan, np.nan
                continue
            both = pd.concat([pm[cls], rm], axis=1, keys=["cls", "ref"]).dropna()
            n_eval[cls] = len(both)
            n_higher[cls] = int((both.cls > both.ref).sum())
        out[f"n_patients_eval_{name}"] = pd.Series(n_eval)
        out[f"n_patients_higher_{name}"] = pd.Series(n_higher)
    return out


def chat_sensitivity(gex: pd.DataFrame, hiconf: bool) -> pd.DataFrame:
    """Valuta la sensibilita' dello stimatore c_hat = (f - q)/(pi - q) al variare della classe di riferimento q."""
    col = "has_tcr_hc" if hiconf else "has_tcr"
    rows = []
    for ref in Q_REFS:
        sub = []
        for pat, g in gex.groupby("patients", observed=True):
            t = g[g.all_celltype.isin(["CD4T", "CD8T"])]
            r = g[g.all_celltype == ref]
            nk = g[g.all_celltype == "NK"]
            if min(len(t), len(r), len(nk)) < 30:
                continue
            pi, q, f = t[col].mean(), r[col].mean(), nk[col].mean()
            if pi <= q:
                continue
            sub.append(np.clip((f - q) / (pi - q), 0, 1))
        if sub:
            rows.append(
                {
                    "q_from": ref,
                    "n_patients": len(sub),
                    "c_hat_median": float(np.median(sub)),
                    "c_hat_q25": float(np.percentile(sub, 25)),
                    "c_hat_q75": float(np.percentile(sub, 75)),
                }
            )
    return pd.DataFrame(rows)


def plot_pairs(per_pair: pd.DataFrame, cc: str) -> None:
    """Figura dell'eccesso per coppia di compartimenti con IC (salta se la tabella e' vuota)."""
    if per_pair.empty:
        print("[!] nessuna coppia di compartimenti con >=3 pazienti: figura saltata.")
        (FIG / f"excess_by_pair_{cc}.png").unlink(missing_ok=True)
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    lab = [f"{a[:4]}↔{b[:4]}" for a, b in zip(per_pair.comp_a, per_pair.comp_b)]
    # IC con <= SMALL_K pazienti: non calibrati, disegnati in grigio e tratteggiati
    cal = (per_pair.n_patients > SMALL_K).values
    for ok, color, label in (
        (True, "#d95f02", f"IC 95% (> {SMALL_K} pazienti)"),
        (False, "#999999", f"IC non calibrato (<= {SMALL_K} pazienti)"),
    ):
        m = cal == ok
        if not m.any():
            continue
        pp = per_pair[m]
        bars = ax.errorbar(
            pp.excess,
            np.arange(len(per_pair))[m],
            xerr=[pp.excess - pp.lo, pp.hi - pp.excess],
            fmt="o",
            capsize=4,
            color=color,
            label=label,
        )
        if not ok:
            bars[-1][0].set_linestyle("--")
    ax.legend(fontsize=8, loc="best")
    ax.axvline(0, ls="--", c="k", lw=0.8)
    ax.set_yticks(range(len(per_pair)), lab)
    ax.set_xlabel("eccesso di discordanza (IC95%, cluster bootstrap sui pazienti)")
    fig.tight_layout()
    fig.savefig(FIG / f"excess_by_pair_{cc}.png", dpi=140)
    plt.close(fig)


def main() -> None:
    """Punto di ingresso principale per i test di robustezza statistica e biologica."""
    ap = argparse.ArgumentParser(
        description="Discordanza clonale a coppie, cluster bootstrap e controlli di robustezza.",
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
    ap.add_argument(
        "--clone-def",
        choices=["clone_b", "clone_ab"],
        default="clone_b",
        help="definizione di clone: catena TRB (clone_b) o coppia TRA+TRB (clone_ab)",
    )
    ap.add_argument(
        "--umi-min",
        type=int,
        default=3,
        help="soglia UMI(TRB) per TCR ad alta confidenza",
    )
    ap.add_argument(
        "--no-markers",
        action="store_true",
        help="salta il controllo CD4/CD8 (non carica X)",
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    tcr = load_tcr(args.tcrdir)
    print(f"[ok] TCR: {len(tcr):,} cellule | appaiate TRA+TRB: {tcr.paired.mean():.1%}")

    gex = load_gex(args.h5ad, with_markers=not args.no_markers)
    gex = gex.merge(
        tcr[
            [
                "patient",
                "tissue",
                "barcode16",
                "clone_b",
                "clone_ab",
                "paired",
                "umi_TRA",
                "umi_TRB",
            ]
        ],
        left_on=["patients", "tissue", "barcode16"],
        right_on=["patient", "tissue", "barcode16"],
        how="left",
    )
    gex["has_tcr"] = gex.clone_b.notna()
    gex["has_tcr_hc"] = (
        gex.has_tcr & gex.paired & (gex.umi_TRB.fillna(0) >= args.umi_min)
    )
    gex["log_umi"] = np.log10(gex.umi_TRB.fillna(0) + 1)
    print(
        f"[ok] GEX: {len(gex):,} | TCR: {gex.has_tcr.sum():,} | "
        f"TCR alta confidenza: {gex.has_tcr_hc.sum():,}"
    )

    # Sezione (C): discriminazione tra cellule T reali e contaminazione da RNA ambientale
    print("\n" + "=" * 78)
    print("(C) VERA CELLULA T MAL ETICHETTATA  vs  CONTAMINAZIONE AMBIENTALE")
    print("=" * 78)
    print(
        "Attesa: T reale -> molti UMI, catene appaiate. Ambient -> pochi UMI, 1 catena.\n"
    )
    amb = ambient_discriminant(gex)
    print("(colonne n_patients_*: confronto descrittivo per paziente, nessun test)")
    print(amb.round(4).to_string())
    amb.to_csv(TAB / "ambient_discriminant.csv")

    print("\n>>> tasso di TCR per etichetta: grezzo vs alta confidenza")
    cmp = pd.DataFrame(
        {
            "raw": gex.groupby("all_celltype", observed=True).has_tcr.mean(),
            "high_conf": gex.groupby("all_celltype", observed=True).has_tcr_hc.mean(),
            "n": gex.groupby("all_celltype", observed=True).size(),
        }
    ).sort_values("raw", ascending=False)
    cmp["ratio"] = (cmp.high_conf / cmp.raw.replace(0, np.nan)).round(3)
    print(cmp.round(4).to_string())
    print(
        "\nUna classe il cui tasso si abbatte ad alta confidenza e' compatibile con segnale ambientale;"
    )
    print(
        "un tasso che persiste ad alta confidenza e' compatibile con TCR di cellule T reali."
    )

    # Analisi di sensibilita' dello stimatore c_hat al variare della classe di riferimento q
    print("\n" + "=" * 78)
    print("c_hat DELLA CONTAMINAZIONE NK - sensibilita' alla scelta di q")
    print("=" * 78)
    for hc in (False, True):
        s = chat_sensitivity(gex, hiconf=hc)
        print(f"\nTCR {'ALTA CONFIDENZA' if hc else 'grezzo'}:")
        print(s.round(4).to_string(index=False))
        s.to_csv(TAB / f"chat_sensitivity_{'hc' if hc else 'raw'}.csv", index=False)
    print(
        "\nLa dipendenza di c_hat dalla classe di riferimento q documenta la sensibilita' al rumore di fondo."
    )
    print("Si riporta pertanto l'intero intervallo inter-riferimento.")

    # Sezioni (A) e (B): stimatore non distorto a coppie e cluster bootstrap sui pazienti
    cc = args.clone_def
    lym = gex[gex[cc].notna() & gex.all_celltype.isin(LYMPHOID)].copy()
    print("\n" + "=" * 78)
    print(f"(A) ECCESSO DI DISCORDANZA - stimatore a coppie, non distorto ({cc})")
    print("=" * 78)

    pw = pairwise_discordance(lym, cc)
    if pw.empty:
        print(
            "[!] nessun clone con >=2 cellule in due compartimenti: sezioni A-B saltate."
        )
        # niente tabelle o figure di un'esecuzione precedente in results/
        for path in (
            TAB / f"pairwise_discordance_{cc}.csv",
            TAB / f"excess_summary_{cc}.csv",
            TAB / f"excess_by_compartment_pair_{cc}.csv",
            TAB / f"discordant_clone_markers_{cc}.csv",
            FIG / f"excess_by_pair_{cc}.png",
        ):
            path.unlink(missing_ok=True)
        return
    pw.to_csv(TAB / f"pairwise_discordance_{cc}.csv", index=False)

    print(
        f"coppie clone x compartimento: {len(pw):,} | pazienti: {pw.patient.nunique()}"
    )
    print(f"D_across  medio : {pw.d_across.mean():.4f}")
    print(
        f"D_within  medio : {pw.d_within.mean():.4f}   (basale appaiato, stesso clone)"
    )
    print(f"ECCESSO   medio : {pw.excess.mean():+.4f}")

    print("\n" + "=" * 78)
    print("(B) INTERVALLO DI CONFIDENZA - l'unita' indipendente e' il PAZIENTE")
    print("=" * 78)
    m, lo_n, hi_n = naive_bootstrap(pw)
    _, lo_c, hi_c = cluster_bootstrap(pw)
    print(
        f"bootstrap ingenuo (sui cloni)   : {m:+.4f}  IC95% [{lo_n:+.4f}, {hi_n:+.4f}]"
    )
    print(
        f"cluster bootstrap (sui pazienti): {m:+.4f}  IC95% [{lo_c:+.4f}, {hi_c:+.4f}]"
    )
    print(f"l'IC si allarga di {((hi_c - lo_c) / max(hi_n - lo_n, 1e-12)):.1f}x")
    m_eq, lo_eq, hi_eq = cluster_bootstrap(pw, "excess_eq")
    print(
        f"sensibilita', pesi uguali 1/2||pA-pB||^2: {m_eq:+.4f}  IC95% [{lo_eq:+.4f}, {hi_eq:+.4f}]"
    )
    pd.DataFrame(
        [
            {
                "clone_def": cc,
                "n_pairs": len(pw),
                "n_patients": pw.patient.nunique(),
                "d_across": pw.d_across.mean(),
                "d_within": pw.d_within.mean(),
                "excess": m,
                "lo_cluster": lo_c,
                "hi_cluster": hi_c,
                "lo_naive": lo_n,
                "hi_naive": hi_n,
                "excess_eq": m_eq,
                "lo_eq_cluster": lo_eq,
                "hi_eq_cluster": hi_eq,
            }
        ]
    ).to_csv(TAB / f"excess_summary_{cc}.csv", index=False)

    print("\nper coppia di compartimenti (cluster bootstrap):")
    rows = []
    for (a, b), g in pw.groupby(["comp_a", "comp_b"]):
        if g.patient.nunique() < 3:
            continue
        mm, l, h = cluster_bootstrap(g)  # noqa: E741
        rows.append(
            {
                "comp_a": a,
                "comp_b": b,
                "n_pairs": len(g),
                "n_patients": g.patient.nunique(),
                "excess": mm,
                "lo": l,
                "hi": h,
                # con <= SMALL_K pazienti l'IC non e' calibrato: niente verdetto
                "ic_calibrato": g.patient.nunique() > SMALL_K,
                "significativo": bool(l > 0)
                if g.patient.nunique() > SMALL_K
                else np.nan,
            }
        )
    per_pair = pd.DataFrame(rows)
    print(per_pair.round(4).to_string(index=False))
    per_pair.to_csv(TAB / f"excess_by_compartment_pair_{cc}.csv", index=False)

    print("\nper paziente (eccesso medio):")
    print(pw.groupby("patient").excess.agg(["mean", "size"]).round(4).to_string())

    # Sezione (D): verifica dell'espressione (RNA, CP10k) di CD4 e CD8A/CD8B nei cloni discordanti
    if not args.no_markers and "g_CD8A" in gex.columns:
        print("\n" + "=" * 78)
        print("(D) TRANSIZIONE CD4T -> CD8T: verifica sui marcatori canonici")
        print("=" * 78)
        disc = pw[pw.excess > 0.5]
        keys = set(zip(disc.patient, disc.clone))
        sel = lym[[(p, c) in keys for p, c in zip(lym.patient, lym[cc])]]
        if len(sel) > 20:
            tab = sel.groupby(["tissue", "all_celltype"], observed=True)[
                ["g_CD4", "g_CD8A", "g_CD8B", "g_CD3D", "g_NKG7", "g_GNLY"]
            ].mean()
            print(
                f"cloni fortemente discordanti: {len(keys)} | cellule: {len(sel):,}\n"
            )
            print(tab.round(3).to_string())
            # Lettura condizionata ai dati: CD8A nelle CD4T rispetto alle CD8T degli stessi cloni
            cd4 = sel.loc[sel.all_celltype == "CD4T", "g_CD8A"]
            cd8 = sel.loc[sel.all_celltype == "CD8T", "g_CD8A"]
            if len(cd4) and len(cd8):
                r = cd4.mean() / max(cd8.mean(), 1e-9)
                print(
                    f"\nCD8A medio (CP10k): CD4T {cd4.mean():.2f} vs CD8T {cd8.mean():.2f} "
                    f"(rapporto {r:.2f})"
                )
                if r > 0.5:
                    print(
                        "  Le CD4T di questi cloni esprimono CD8A quasi quanto le CD8T:"
                    )
                    print("  compatibile con un'etichetta CD4T sbagliata.")
                else:
                    print(
                        "  Le CD4T di questi cloni esprimono poco CD8A: l'etichetta CD4T e'"
                    )
                    print(
                        "  coerente con i marcatori; la discordanza va cercata nelle cellule"
                    )
                    print("  etichettate CD8T o NK (vedi 08_disambiguate.py).")
            tab.to_csv(TAB / f"discordant_clone_markers_{cc}.csv")
        else:
            print(
                "numero insufficiente di cloni fortemente discordanti per la verifica."
            )
            (TAB / f"discordant_clone_markers_{cc}.csv").unlink(missing_ok=True)

    # Generazione delle figure diagnostiche (UMI TCR e intervalli di confidenza per coppia di compartimenti)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    order = [
        c
        for c in [
            "CD8T",
            "CD4T",
            "NK",
            "Stellate",
            "Endothelial",
            "Platelet",
            "Myeloid",
            "Ductal",
            "B",
        ]
        if c in gex.all_celltype.unique()
    ]
    d = gex[gex.has_tcr & gex.all_celltype.isin(order)]
    axes[0].boxplot(
        [d.loc[d.all_celltype == c, "log_umi"] for c in order],
        tick_labels=order,
        showfliers=False,
    )
    axes[0].set_ylabel("log10(UMI TRB + 1)")
    axes[0].set_title("UMI del TCR: T reale vs ambient")
    axes[0].tick_params(axis="x", rotation=60)

    fp = d.groupby("all_celltype", observed=True).paired.mean().reindex(order)
    axes[1].bar(order, fp.values, color="#1b9e77", alpha=0.85)
    axes[1].set_ylabel("frazione con TRA+TRB appaiate")
    axes[1].set_title("catene appaiate")
    axes[1].tick_params(axis="x", rotation=60)
    fig.tight_layout()
    fig.savefig(FIG / "ambient_vs_real_tcr.png", dpi=140)
    plt.close(fig)

    plot_pairs(per_pair, cc)

    print("\n" + "=" * 78)
    if lo_c > 0:
        print(
            f"ESITO STATISTICO: eccesso {m:+.4f}, IC95% [{lo_c:+.4f}, {hi_c:+.4f}] con cluster bootstrap."
        )
        print(
            "  Con il paziente come unita' di ricampionamento l'IC al 95% esclude lo zero."
        )
    else:
        print(f"ESITO STATISTICO: IC95% [{lo_c:+.4f}, {hi_c:+.4f}] include lo zero.")
        print(
            "  Con il ricampionamento clusterizzato per paziente l'intervallo include lo zero."
        )
    print("=" * 78)


if __name__ == "__main__":
    main()

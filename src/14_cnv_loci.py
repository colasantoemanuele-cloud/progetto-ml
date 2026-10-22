#!/usr/bin/env python3
"""
14_cnv_loci.py - Score CNV sui loci ricorrenti del PDAC: separazione fra duttali tumorali e adiacenti.

Scopo
-----
Lo score CNV globale (`mean |X_cnv|`) aggrega l'intero genoma e risente dei
programmi trascrizionali tessuto-specifici. Questo modulo restringe il calcolo ai
loci con alterazioni ricorrenti nell'adenocarcinoma duttale pancreatico e misura
l'AUC tumore-vs-adiacente delle duttali, sia pooled sia entro paziente:
  - perdite attese (`direction = -1`): 9p21 (`CDKN2A`), 17p13 (`TP53`), 18q21 (`SMAD4`),
    `1p_loss` (0-60 Mb) e `6q_loss` (90-170 Mb): questi due coprono la parte distale
    dei bracci, non i bracci interi (le finestre attribuite, larghe, arrivano a 84 Mb
    su 1p e partono da 64 Mb su 6q, vedi `cnv_windows.csv`)
  - guadagni attesi (`direction = +1`): 8q24 (`MYC`), 12p12 (`KRAS`), 20q
Calcola uno score direzionale per locus e un indice composito (`PDAC_COMPOSITE`):
la somma degli z-score dei loci, con media e deviazione standard di ciascun locus
stimate su tutte le cellule di `X_cnv` (duttali e immunitarie di riferimento).

Come in `13_infercnv.py`, inferCNV gira una sola volta con il riferimento formato
dalle immunitarie di TUTTI i pazienti e con la stessa soglia del rumore globale
(`global_noise_filter` di 13, indipendente dall'ordine delle cellule). Le coordinate
delle finestre non sono restituite da infercnvpy: sono ricostruite dai geni
(`window_table` di 13, stessa costruzione di infercnvpy: la finestra i di un
cromosoma copre i geni ordinati per start [i*step, i*step + 100)), con start del
primo gene, end massimo e mediana degli start; sono salvate in `cnv_windows.csv`. Una
finestra appartiene a un locus se il suo intervallo genomico si sovrappone al locus;
le finestre sono larghe (100 geni: nelle regioni povere di geni, come 9p, fino a
decine di Mb), quindi un locus include anche geni vicini e il segnale di un locus
piccolo (CDKN2A) e' in pratica quello del braccio cromosomico circostante. Lo
script controlla che le finestre che contengono il gene marcatore di ogni locus
(CDKN2A, TP53, SMAD4, MYC, KRAS) siano attribuite al locus, e si ferma se no.
RNA ambientale acinare: CTRC, CELA2A/2B (1p36.21) e CELA3A/3B (1p36.12) cadono dentro
`1p_loss`, PRSS1/2 in 7q34. L'RNA ambientale acinare, abbondante nel tessuto adiacente,
puo' creare pseudo-guadagni in queste regioni nelle duttali adiacenti e influenzare
l'AUC tumore contro adiacente, in particolare per `1p_loss`.
L'AUC pooled mescola differenze fra pazienti e fra tessuti; l'AUC entro paziente
(pazienti con almeno 10 duttali per tessuto) elimina le differenze fra pazienti
(resta il contrasto fra campioni e librerie diversi dello stesso paziente). Il
messaggio finale tiene conto di entrambe: una separazione pooled che non si ritrova
entro paziente non e' attribuita ai loci, e il confronto fra loci PDAC e score medio
si fa sulle AUC entro paziente di entrambi (non sulle pooled).

Coordinate dei geni e metodo: le stesse di `13_infercnv.py` (Ensembl release 116,
tabella in `data/interim/gene_positions_ensembl116.csv`, funzioni importate da 13);
senza coordinate lo script si ferma. Con `--method copykat` le coordinate delle
finestre non sono ricostruibili e lo script si ferma (copykat richiede anche `rpy2`
e R, assenti dall'ambiente). Cromosomi: come in 13, entrano solo i geni dei
cromosomi 1-22, X e Y (esclusi `MT`, la cui espressione segue la frazione
mitocondriale `pct_mt`, e i contig non standard); infercnvpy esclude X e Y, e lo
script si ferma se `X_cnv` contiene finestre fuori da chr1-22. La cache e' riusata solo
se la sua impronta (voce `cnv_matrix.npz` di `data/interim/cnv_cache_fingerprint.json`,
stesso calcolo di 13) coincide.

Input
-----
  --h5ad  data/interim/sc_raw.h5ad
  cache   data/interim/cnv_matrix.npz (letta se esiste, contiene esattamente le cellule
          selezionate e ha l'impronta attesa; se manca e' calcolata e scritta, una cache
          diversa ferma lo script)
  data/interim/cnv_cache_fingerprint.json (impronte delle cache di 13 e 14)
  data/interim/gene_positions_ensembl116.csv (letta se esiste, altrimenti scaricata da BioMart)

Output
------
  data/interim/cnv_matrix.npz                    matrice X_cnv e coordinate delle finestre (cache)
  results/tables/cnv_windows.csv                 coordinate reali di ogni finestra (colonna di X_cnv,
                                                 cromosoma, start, end, mediana, n geni, primo e
                                                 ultimo gene) e loci a cui e' attribuita
  data/interim/cnv_cache_fingerprint.json        impronta della cache (voce `cnv_matrix.npz`)
  data/interim/gene_positions_ensembl116.csv     coordinate dei geni (solo se scaricate)
  results/tables/cnv_loci_auc.csv                AUC pooled per segnale (score medio, loci, composito);
                                                 n_windows = finestre usate, n_loci = loci nel segnale
  results/tables/cnv_loci_auc_by_patient.csv     AUC entro paziente di PDAC_COMPOSITE e dello score medio
                                                 (ultima riga `_mean`: media delle AUC sui pazienti
                                                 valutabili; n_tumor e n_normal sono le somme)
  results/figures/cnv_loci_separation.png, results/figures/cnv_loci_auc.png

Argomenti
---------
  --h5ad (obbligatorio), --method {infercnv,copykat,auto} (default auto = infercnv,
  senza ripiego automatico su copykat), --cache (default data/interim/cnv_matrix.npz)

Uso
---
    python src/14_cnv_loci.py --h5ad data/interim/sc_raw.h5ad
"""

from __future__ import annotations

import argparse
import importlib.util
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from sklearn.metrics import roc_auc_score

FIG = Path("results/figures")
TAB = Path("results/tables")
CACHE = Path("data/interim/cnv_matrix.npz")
SEED = 0
HERE = Path(__file__).resolve().parent

IMMUNE_REF = ["CD8T", "CD4T", "NK", "B", "Myeloid", "Plasma"]
# colonne delle coordinate delle finestre (da `window_table` di 13), salvate nella cache
WIN_COLS = ["col", "chrom", "start", "end", "pos", "n_genes", "first_gene", "last_gene"]

# Coordinate genomiche (GRCh38/hg38, in bp) e direzionalita' attesa dei loci citogenetici PDAC
PDAC_LOCI = [
    {
        "name": "CDKN2A_9p21",
        "gene": "CDKN2A",
        "chrom": "9",
        "start": 20e6,
        "end": 24e6,
        "direction": -1,
    },
    {
        "name": "TP53_17p13",
        "gene": "TP53",
        "chrom": "17",
        "start": 0,
        "end": 10e6,
        "direction": -1,
    },
    {
        "name": "SMAD4_18q21",
        "gene": "SMAD4",
        "chrom": "18",
        "start": 45e6,
        "end": 60e6,
        "direction": -1,
    },
    {"name": "1p_loss", "chrom": "1", "start": 0, "end": 60e6, "direction": -1},
    {"name": "6q_loss", "chrom": "6", "start": 90e6, "end": 170e6, "direction": -1},
    {
        "name": "MYC_8q24",
        "gene": "MYC",
        "chrom": "8",
        "start": 120e6,
        "end": 145e6,
        "direction": +1,
    },
    {
        "name": "KRAS_12p12",
        "gene": "KRAS",
        "chrom": "12",
        "start": 20e6,
        "end": 30e6,
        "direction": +1,
    },
    {"name": "20q_gain", "chrom": "20", "start": 30e6, "end": 65e6, "direction": +1},
]


def _infercnv13():
    """Importa `13_infercnv.py` (stesse coordinate dei geni e stesso controllo di copykat)."""
    spec = importlib.util.spec_from_file_location("infercnv13", HERE / "13_infercnv.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def compute_cnv(sub: ad.AnnData, method: str, m13):
    """Calcola `X_cnv` (soglia del rumore globale di 13) e le coordinate reali delle finestre.

    `sub.var` deve avere le coordinate (`chromosome`, `start`, `end`) di `add_positions`
    di 13. Le coordinate delle finestre sono ricostruite dai geni con `m13.window_table`
    (stessa costruzione di infercnvpy); con copykat non sono ricostruibili e lo script
    si ferma.
    """
    if "chromosome" not in sub.var:
        raise SystemExit(
            "[ERRORE] coordinate dei geni assenti: chiamare add_positions di 13."
        )
    if method != "infercnv":
        raise SystemExit(
            "[ERRORE] le coordinate delle finestre sono ricostruibili solo con inferCNV: "
            "usa --method infercnv."
        )
    a = sub.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a = a[:, a.var["chromosome"].notna()].copy()
    ref = m13.reference_categories(a.obs["ref_group"])
    print(f"[..] infercnv, riferimento = {ref} di tutti i pazienti (lento) ...")
    m13.run_infercnv_core(a, ref, m13.WINDOW_SIZE)
    X = a.obsm["X_cnv"]
    X = np.asarray(X.todense()) if issparse(X) else np.asarray(X)
    windows = m13.window_table(a, m13.WINDOW_SIZE)
    return X, windows, a.obs_names


def locus_windows(windows: pd.DataFrame, locus: dict) -> np.ndarray:
    """Finestre del locus: quelle il cui intervallo genomico (start del primo gene, end
    massimo dei suoi geni) si sovrappone al locus."""
    return (
        (windows.chrom == locus["chrom"])
        & (windows.end >= locus["start"])
        & (windows.start <= locus["end"])
    ).values


def check_marker_windows(windows: pd.DataFrame, genes: pd.DataFrame) -> None:
    """Controllo esplicito: le finestre che contengono il gene marcatore di un locus
    (es. CDKN2A, chr9:21,97-22,00 Mb) sono fra quelle attribuite al locus."""
    for L in PDAC_LOCI:
        g = L.get("gene")
        if g is None:
            continue
        if g not in genes.index:
            print(f"  {L['name']}: gene {g} senza coordinate, controllo non eseguito")
            continue
        gs, ge = float(genes.loc[g, "start"]), float(genes.loc[g, "end"])
        has_gene = (
            (windows.chrom == L["chrom"]) & (windows.start <= gs) & (windows.end >= ge)
        ).values
        in_locus = locus_windows(windows, L)
        ok = has_gene.any() and bool(np.all(in_locus[has_gene]))
        print(
            f"  {L['name']}: {g} a chr{L['chrom']}:{gs / 1e6:.2f}-{ge / 1e6:.2f} Mb in "
            f"{int(has_gene.sum())} finestre, tutte attribuite al locus: {ok} "
            f"(finestre del locus: {int(in_locus.sum())})"
        )
        if not ok:
            raise SystemExit(
                f"[ERRORE] il gene {g} non cade nelle finestre di {L['name']}"
            )


def locus_score(X, windows, locus) -> np.ndarray:
    """Calcola lo score CNV medio moltiplicato per la direzione citogenetica attesa (+1 guadagno, -1 perdita)."""
    m = locus_windows(windows, locus)
    if m.sum() == 0:
        return np.full(X.shape[0], np.nan)
    return X[:, m].mean(axis=1) * locus["direction"]


def main() -> None:
    """Punto di ingresso principale per la valutazione citogenetica locus-specifica."""
    ap = argparse.ArgumentParser(
        description="Score CNV sui loci PDAC: AUC tumore vs adiacente, pooled ed entro paziente.",
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
        "--cache",
        type=Path,
        default=CACHE,
        help="cache di X_cnv (con impronta in cnv_cache_fingerprint.json)",
    )
    args = ap.parse_args()
    cache = args.cache

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    for c in ("tissue", "patients", "all_celltype"):
        adata.obs[c] = adata.obs[c].astype(str)

    is_duct = adata.obs.all_celltype == "Ductal"
    is_imm = adata.obs.all_celltype.isin(IMMUNE_REF)
    keep = (is_duct | is_imm) & adata.obs.tissue.isin(["Tumor", "Adjacent_normal"])
    sub = adata[keep].copy()
    sub.obs["ref_group"] = np.where(
        sub.obs.all_celltype.isin(IMMUNE_REF), sub.obs.all_celltype, "Ductal"
    )

    # Caricamento o calcolo con salvataggio compresso (.npz) della matrice delle finestre genomiche X_cnv
    m13 = _infercnv13()
    method = "infercnv" if args.method == "auto" else args.method
    if cache.exists():
        print(f"[ok] X_cnv da cache {cache}")
        npz = np.load(cache, allow_pickle=True)
        X, obs_names = npz["X"], npz["obs_names"]
        if "win_start" not in npz:
            raise SystemExit(
                f"[ERRORE] la cache {cache} non ha le coordinate reali delle finestre "
                "(versione precedente): cancellala per ricalcolare."
            )
        windows = pd.DataFrame({c: npz[f"win_{c}"] for c in WIN_COLS})
        cached, wanted = pd.Index(obs_names), sub.obs_names
        if len(cached.difference(wanted)) or len(wanted.difference(cached)):
            raise SystemExit(
                f"[ERRORE] la cache {cache} non corrisponde alle cellule selezionate "
                f"({len(wanted.difference(cached)):,} mancanti, "
                f"{len(cached.difference(wanted)):,} in piu'): cancellala per ricalcolare."
            )
        m13.check_fingerprint(cache.name, m13.cnv_fingerprint(sub, method), cache)
        sub = sub[obs_names].copy()
    else:
        if method == "copykat":
            m13.require_copykat()
        # coordinate necessarie con entrambi i metodi (mappa finestre -> loci)
        if not m13.add_positions(sub):
            raise SystemExit(
                f"[ERRORE] coordinate dei geni (Ensembl {m13.ENSEMBL_RELEASE}) non "
                f"disponibili: nessun host BioMart risponde e {m13.GENE_POS} non esiste. "
                "copykat non e' usato come ripiego (richiede rpy2 e R)."
            )
        t0 = time.time()
        X, windows, obs_names = compute_cnv(sub, method, m13)
        print(f"[ok] X_cnv {X.shape} in {time.time() - t0:.0f}s")
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            cache,
            X=X,
            obs_names=np.array(obs_names),
            **{f"win_{c}": windows[c].values for c in WIN_COLS},
        )
        m13.save_fingerprint(cache.name, m13.cnv_fingerprint(sub, method))
        sub = sub[obs_names].copy()

    # Coordinate delle finestre: tabella e controllo dei geni marcatori dei loci
    pos_tab = m13.gene_positions()
    if pos_tab is None:
        raise SystemExit(
            f"[ERRORE] {m13.GENE_POS} assente: servono le coordinate dei geni."
        )
    pos_tab = pos_tab.dropna(subset=["hgnc_symbol"])
    pos_tab = pos_tab[~pos_tab.hgnc_symbol.duplicated(keep=False)].set_index(
        "hgnc_symbol"
    )
    gene_coord = pos_tab.rename(
        columns={"start_position": "start", "end_position": "end"}
    )
    win_out = windows.copy()
    win_out["loci"] = [
        ";".join(L["name"] for L in PDAC_LOCI if locus_windows(windows, L)[k])
        for k in range(len(windows))
    ]
    win_out.to_csv(TAB / "cnv_windows.csv", index=False)
    print(
        "\ncontrollo delle finestre dei loci (intervallo della finestra sovrapposto al locus):"
    )
    check_marker_windows(windows, gene_coord)

    obs = sub.obs
    duct = obs.all_celltype.values == "Ductal"
    tumor = duct & (obs.tissue.values == "Tumor")
    normal = duct & (obs.tissue.values == "Adjacent_normal")
    print(f"\nduttali: {tumor.sum():,} tumorali | {normal.sum():,} normali adiacenti")
    if tumor.sum() == 0 or normal.sum() == 0:
        raise SystemExit(
            "[ERRORE] servono duttali sia in Tumor sia in Adjacent_normal."
        )

    print("\n" + "=" * 74)
    print("SEPARAZIONE TUMORE vs NORMALE: score medio vs loci PDAC specifici")
    print("=" * 74)
    print("AUC = capacita' discriminativa tra duttali TUMORALI e NORMALI adiacenti.")
    print("0.5 = distribuzioni sovrapposte; l'AUC pooled mescola pazienti diversi\n")

    def auc_tn(score):
        """AUC di `score` fra duttali tumorali e normali adiacenti (NaN esclusi);
        NaN se restano <= 10 valori validi o una sola classe."""
        y = np.r_[np.ones(tumor.sum()), np.zeros(normal.sum())]
        idx = np.r_[np.where(tumor)[0], np.where(normal)[0]]
        s = score[idx]
        ok = ~np.isnan(s)
        return (
            roc_auc_score(y[ok], s[ok])
            if ok.sum() > 10 and len(np.unique(y[ok])) > 1
            else np.nan
        )

    mean_score = np.abs(X).mean(axis=1)
    rows = [
        {
            "signal": "MEAN_|CNV|",
            "auc_tumor_vs_normal": auc_tn(mean_score),
            "n_windows": X.shape[1],
            "n_loci": np.nan,
        }
    ]
    loci_scores = {}
    for L in PDAC_LOCI:
        s = locus_score(X, windows, L)
        loci_scores[L["name"]] = s
        nwin = int(locus_windows(windows, L).sum())
        rows.append(
            {
                "signal": L["name"],
                "auc_tumor_vs_normal": auc_tn(s),
                "n_windows": nwin,
                "n_loci": 1,
            }
        )

    # Somma standardizzata (z-score) dei contributi direzionali dei singoli loci PDAC:
    # media e deviazione standard di ogni locus su tutte le righe di X (duttali e immuni)
    valid_names = [k for k, s in loci_scores.items() if not np.all(np.isnan(s))]
    valid = [loci_scores[k] for k in valid_names]
    if not valid:
        raise SystemExit("[ERRORE] nessun locus PDAC coperto dalle finestre di X_cnv.")
    comp = np.nansum(
        [(s - np.nanmean(s)) / (np.nanstd(s) + 1e-9) for s in valid], axis=0
    )
    rows.append(
        {
            "signal": "PDAC_COMPOSITE",
            "auc_tumor_vs_normal": auc_tn(comp),
            "n_windows": int(
                sum(r["n_windows"] for r in rows if r["signal"] in valid_names)
            ),
            "n_loci": len(valid),
        }
    )

    res = pd.DataFrame(rows)
    print(res.round(3).to_string(index=False))
    res.to_csv(TAB / "cnv_loci_auc.csv", index=False)

    mean_auc = res.loc[res.signal == "MEAN_|CNV|", "auc_tumor_vs_normal"].iloc[0]
    comp_auc = res.loc[res.signal == "PDAC_COMPOSITE", "auc_tumor_vs_normal"].iloc[0]

    print(
        "\nAUC tumore-vs-normale PDAC_COMPOSITE, per paziente con abbastanza normali:"
    )
    pat = obs.patients.values
    by_pat = []
    for p in np.unique(pat[normal]):
        t = tumor & (pat == p)
        n = normal & (pat == p)
        if n.sum() < 10 or t.sum() < 10:
            print(
                f"  {p}: non valutabile (tumor n={t.sum()}, normal n={n.sum()}, minimo 10)"
            )
            continue
        idx = np.r_[np.where(t)[0], np.where(n)[0]]
        y = np.r_[np.ones(t.sum()), np.zeros(n.sum())]
        row = {"patient": p, "n_tumor": int(t.sum()), "n_normal": int(n.sum())}
        for name, score in (("PDAC_COMPOSITE", comp), ("MEAN_|CNV|", mean_score)):
            s = score[idx]
            ok = ~np.isnan(s)
            row[f"auc_{name}"] = (
                roc_auc_score(y[ok], s[ok]) if len(np.unique(y[ok])) > 1 else np.nan
            )
        by_pat.append(row)
        print(
            f"  {p}: AUC={row['auc_PDAC_COMPOSITE']:.3f}  "
            f"(tumor n={t.sum()}, normal n={n.sum()})"
        )
    by_pat = pd.DataFrame(by_pat)
    if not by_pat.empty:
        # riga di riepilogo: media delle AUC entro paziente, somma delle cellule
        mean_row = by_pat.drop(columns="patient").mean(numeric_only=True).to_dict()
        mean_row["n_tumor"] = int(by_pat.n_tumor.sum())
        mean_row["n_normal"] = int(by_pat.n_normal.sum())
        by_pat = pd.concat(
            [by_pat, pd.DataFrame([{"patient": "_mean", **mean_row}])],
            ignore_index=True,
        )
    by_pat.to_csv(TAB / "cnv_loci_auc_by_patient.csv", index=False)
    if not by_pat.empty:
        by_pat = by_pat[by_pat.patient != "_mean"]
    if not by_pat.empty:
        print(
            f"  AUC entro paziente PDAC_COMPOSITE: media {by_pat.auc_PDAC_COMPOSITE.mean():.3f}, "
            f"range [{by_pat.auc_PDAC_COMPOSITE.min():.3f}, {by_pat.auc_PDAC_COMPOSITE.max():.3f}] "
            f"su {len(by_pat)} pazienti"
        )
    within = by_pat.auc_PDAC_COMPOSITE.mean() if not by_pat.empty else np.nan
    # confronto con lo score medio alla stessa scala (entro paziente, non pooled)
    within_mean = by_pat["auc_MEAN_|CNV|"].mean() if not by_pat.empty else np.nan

    # Esito: la separazione pooled e' attribuita ai loci solo se si ritrova entro paziente
    print("\n" + "=" * 74)
    print(f"AUC PDAC_COMPOSITE pooled {comp_auc:.2f} (score medio {mean_auc:.2f})")
    if np.isnan(within):
        print("ESITO: nessun paziente valutabile entro paziente (minimo 10 duttali per")
        print(
            "  tessuto): l'AUC pooled mescola differenze fra pazienti, nessuna conclusione"
        )
        print("  sui loci.")
    else:
        print(
            f"AUC entro paziente (media su {len(by_pat)} pazienti): PDAC_COMPOSITE "
            f"{within:.2f}, score medio {within_mean:.2f}"
        )
        if within >= 0.65 and within > within_mean + 0.05:
            print(
                "ESITO: i loci PDAC separano le duttali tumorali dalle adiacenti anche"
            )
            print("  entro paziente, meglio dello score medio entro paziente.")
        elif within >= 0.65:
            print(
                "ESITO: entro paziente le duttali tumorali e adiacenti si separano, ma i"
            )
            print("  loci PDAC non fanno meglio dello score medio entro paziente: la")
            print("  separazione non e' attribuibile ai loci.")
        elif within < 0.6:
            print(
                "ESITO: entro paziente le duttali tumorali e adiacenti si sovrappongono"
            )
            print(
                "  anche sui loci PDAC; l'eventuale separazione pooled riflette soprattutto"
            )
            print(
                "  differenze fra pazienti (compatibile con effetto di campo peri-tumorale"
            )
            print("  o con la risoluzione limitata dell'mRNA, non distinguibili qui).")
        else:
            print(
                "ESITO: separazione entro paziente parziale; con pochi pazienti valutabili"
            )
            print("  non e' attribuibile con sicurezza ai loci.")
    print("=" * 74)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].hist(
        mean_score[tumor],
        bins=40,
        alpha=0.55,
        density=True,
        label="Tumor",
        color="#d95f02",
    )
    axes[0].hist(
        mean_score[normal],
        bins=40,
        alpha=0.55,
        density=True,
        label="Normal adj",
        color="#1b9e77",
    )
    axes[0].set_title(f"score MEDIO |CNV|  (AUC {mean_auc:.2f})")
    axes[0].set_xlabel("mean |CNV|")
    axes[0].legend(fontsize=8)

    axes[1].hist(
        comp[tumor], bins=40, alpha=0.55, density=True, label="Tumor", color="#d95f02"
    )
    axes[1].hist(
        comp[normal],
        bins=40,
        alpha=0.55,
        density=True,
        label="Normal adj",
        color="#1b9e77",
    )
    axes[1].set_title(f"score LOCI PDAC  (AUC {comp_auc:.2f})")
    axes[1].set_xlabel("score composito loci")
    axes[1].legend(fontsize=8)
    fig.suptitle("Separabilita' sui loci citogenetici PDAC vs score medio")
    fig.tight_layout()
    fig.savefig(FIG / "cnv_loci_separation.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    d = (
        res[res.signal != "MEAN_|CNV|"]
        .set_index("signal")
        .auc_tumor_vs_normal.sort_values()
    )
    colors = ["#d95f02" if v > 0.6 else "#999999" for v in d.values]
    ax.barh(d.index, d.values, color=colors, alpha=0.85)
    ax.axvline(0.5, ls="--", c="k", lw=0.8)
    ax.axvline(
        mean_auc, ls=":", c="#d95f02", lw=1.2, label=f"score medio ({mean_auc:.2f})"
    )
    ax.set_xlabel("AUC tumore vs normale adiacente")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "cnv_loci_auc.png", dpi=140)
    plt.close(fig)
    print("\n[ok] figure: results/figures/cnv_loci_separation.png, cnv_loci_auc.png")


if __name__ == "__main__":
    main()

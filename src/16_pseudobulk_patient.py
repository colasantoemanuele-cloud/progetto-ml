#!/usr/bin/env python3
"""
16_pseudobulk_patient.py: Analisi pseudobulk a livello di paziente e test di ipotesi pre-specificate.

Inquadramento Metodologico
--------------------------
Mentre gli script precedenti modellano la singola cellula come unita' osservazionale
correggendo l'inferenza a posteriori per la struttura di clustering (split per gruppo,
cluster bootstrap), questo modulo aggrega ciascun paziente in un unico profilo
pseudobulk prima dell'inferenza statistica. Riducendo le 199.184 cellule alle unita'
cliniche effettivamente indipendenti (n = 10-14 pazienti valutabili), la dipendenza
intra-paziente viene eliminata alla radice per costruzione sperimentale.

Controllo della Molteplicita' su Piccoli Campioni
-------------------------------------------------
Con una numerosita' campionaria pari a n = 10-14 pazienti, una scansione esplorativa
non vincolata sull'intero trascrittoma (~20.000 geni) produrrebbe circa 1.000 falsi
positivi nominali a soglia alfa = 0.05 sotto l'ipotesi nulla globale, rendendo
impossibile distinguere il segnale biologico dal rumore campionario. Per preservare
la validita' dell'inferenza, l'analisi e' rigorosamente circoscritta a tre ipotesi
meccanicistiche formulate a priori sulla base delle evidenze emerse nei moduli 08-10:

  - H1 (Esaurimento CD8T tumorali): il punteggio medio di espressione dei marker di
    esaurimento (PDCD1, LAG3, HAVCR2, TOX, CTLA4, TIGIT) nelle cellule CD8T infiltranti
    correla positivamente con il tasso di errore di annotazione nel tessuto tumorale,
    poiche' la convergenza trascrizionale verso fenotipi disfunzionali attenua il
    confine rispetto alle popolazioni CD4T e NK.
  - H2 (Frazione mieloide tumorale): la proporzione di cellule mieloidi nel
    microambiente tumorale correla positivamente con il tasso di errore, riflettendo
    il contributo di RNA ambientale e doppietti eterotipici mieloide-T.
  - H3 (Frazione stromale CAF): la proporzione di fibroblasti associati al cancro (CAF)
    correla positivamente con il tasso di errore, in quanto una maggiore densita'
    desmoplastica richiede protocolli di dissociazione enzimatica piu' intensi che
    degradano l'integrita' del profilo trascrizionale.

La significativita' statistica e' valutata mediante test di permutazione esatto a due
code (20.000 permutazioni del vettore dei ranghi, esente da assunzioni asintotiche di
normalita') con correzione sequenziale di Holm-Bonferroni sulla famiglia dei 3 test.

Esecuzione:
    python src/16_pseudobulk_patient.py \
        --h5ad data/interim/sc_raw.h5ad \
        --arbitration results/tables/arbitration_blood.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

# Percorsi di output e iperparametri del protocollo inferenziale
FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_PERM = 20000
MIN_CD8T_TUMOR = 30
MIN_CLONES_OUTCOME = 5

# Firma genica canonica di esaurimento funzionale dei linfociti T CD8+
EXHAUSTION_GENES = ["PDCD1", "LAG3", "HAVCR2", "TOX", "CTLA4", "TIGIT"]

# Definizione formale delle tre ipotesi a priori e della direzione attesa dell'effetto
HYPOTHESES = {
    "exhaustion_score": (
        "H1: esaurimento CD8T tumorali -> piu' errore (deriva CD8->CD4/NK)",
        "positive",
    ),
    "myeloid_fraction": (
        "H2: piu' contesto mieloide -> piu' errore (doppietti/ambient)",
        "positive",
    ),
    "caf_fraction": (
        "H3: piu' fibrosi -> qualita' peggiore -> piu' errore (speculativa)",
        "positive",
    ),
}


# --------------------------------------------------------------------------- #
def permutation_corr_test(
    x: np.ndarray, y: np.ndarray, n_perm: int = N_PERM, seed: int = SEED
):
    """Calcola il coefficiente di correlazione di Spearman e il relativo p-value
    a due code mediante test di permutazione vettorizzato sui ranghi.
    L'implementazione calcola simultaneamente la correlazione di Pearson tra i
    ranghi osservati di x e una matrice di n_perm permutazioni indipendenti dei
    ranghi di y, garantendo esattezza a campioni finiti senza approssimazioni normali."""
    rng = np.random.default_rng(seed)
    n = len(x)
    # Trasformazione nei ranghi per ottenere la statistica non parametrica di Spearman
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rho_obs = float(np.corrcoef(rx, ry)[0, 1])

    # Generazione vettorizzata di n_perm permutazioni casuali degli indici tramite argsort
    perm_idx = np.argsort(rng.random((n_perm, n)), axis=1)
    ry_perm = ry[perm_idx]  # (n_perm, n)
    # Centratura dei ranghi e calcolo vettorizzato del coseno angolare (correlazione)
    rx_c = rx - rx.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    perm_rhos = num / np.maximum(den, 1e-12)

    # Stima conservativa di Phipson-Smyth (+1 al numeratore e denominatore) del p-value a due code
    p = (np.sum(np.abs(perm_rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p), perm_rhos


def bootstrap_ci_corr(
    x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = SEED
):
    """Stima l'intervallo di confidenza al 95% per la correlazione di Spearman
    mediante ricampionamento bootstrap non parametrico delle unita' paziente."""
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        # Campionamento con reinserimento delle osservazioni a livello di paziente
        idx = rng.integers(0, n, n)
        # Scarto delle repliche degenerate prive di variabilita' su uno dei due assi
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        r, _ = spearmanr(x[idx], y[idx])
        if np.isfinite(r):
            vals.append(r)
    if len(vals) < 20:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def holm_bonferroni(pvals: dict) -> dict:
    """Applica la procedura step-down di Holm-Bonferroni per il controllo
    forte del Family-Wise Error Rate (FWER) su una famiglia finita di test."""
    # Ordinamento crescente dei p-value grezzi e applicazione dei moltiplicatori (m - i)
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, running_max = {}, 0.0
    for i, (name, p) in enumerate(items):
        # Imposizione del vincolo di monotonia cumulativa sui p-value aggiustati
        val = min((m - i) * p, 1.0)
        running_max = max(running_max, val)
        adj[name] = running_max
    return adj


# --------------------------------------------------------------------------- #
def build_outcome(arb_path: Path) -> pd.DataFrame:
    """Costruisce la variabile di esito a livello di paziente calcolando il tasso
    medio di errore di annotazione nel tessuto tumorale secondo lo schema primario
    con riferimento ematico (PBMC) generato da 10_loco.py."""
    if arb_path.exists():
        j = pd.read_csv(arb_path)
        # Filtro esclusivo sul compartimento tumorale e aggregazione per paziente
        j = j[j.tissue == "Tumor"]
        out = j.groupby("patient").agg(
            error_rate=("error", "mean"), n_clones=("error", "size")
        )
        # Esclusione dei pazienti con un numero di cloni tracciati inferiore alla soglia di stabilita'
        out = out[out.n_clones >= MIN_CLONES_OUTCOME]
        print(f"[ok] outcome da {arb_path.name} (schema 'blood', solo tessuto Tumor)")
        return out.reset_index()
    # Ripiego automatico sull'eccesso di discordanza medio calcolato da 07_robustness.py
    fb = TAB / "pairwise_discordance_clone_b.csv"
    if fb.exists():
        pw = pd.read_csv(fb)
        out = pw.groupby("patient").agg(
            error_rate=("excess", "mean"), n_clones=("excess", "size")
        )
        print(
            f"[!] {arb_path.name} assente, fallback su {fb.name} (eccesso di discordanza, non tasso d'errore)"
        )
        return out.reset_index()
    raise SystemExit(
        f"[ERRORE] ne' {arb_path} ne' {fb} trovati. Esegui prima 10_loco.py (preferito) o 07_robustness.py."
    )


def build_predictors(adata: ad.AnnData) -> pd.DataFrame:
    """Estrae i tre predittori biologici a livello di paziente (punteggio di
    esaurimento CD8T, frazione mieloide e frazione CAF) dal compartimento tumorale."""
    obs = adata.obs
    rows = []  # noqa: F841
    tumor = obs[obs.tissue == "Tumor"]

    # Calcolo delle frequenze relative dei tipi cellulari nel tessuto tumorale di ciascun paziente
    frac = (
        tumor.groupby("patients")
        .all_celltype.value_counts(normalize=True)
        .unstack(fill_value=0.0)
    )
    myeloid_fraction = frac.get("Myeloid", pd.Series(0.0, index=frac.index))
    caf_fraction = frac.get("CAF", pd.Series(0.0, index=frac.index))

    # Mappatura dei simboli genici della firma di esaurimento sugli indici di colonna di AnnData
    sym = (
        adata.var.gene_symbol.astype(str).values
        if "gene_symbol" in adata.var
        else adata.var_names.values
    )
    idx = {g: np.where(sym == g)[0] for g in EXHAUSTION_GENES}
    missing = [g for g, ix in idx.items() if len(ix) == 0]
    if missing:
        print(f"[!] geni di esaurimento assenti dal pannello: {missing}")
    found = [g for g in EXHAUSTION_GENES if g not in missing]
    col_idx = np.concatenate([idx[g] for g in found])

    # Calcolo del punteggio medio di esaurimento log1p(CPM/100) sulle cellule CD8T tumorali per paziente
    exh_scores = {}
    for pat, g in tumor[tumor.all_celltype == "CD8T"].groupby("patients"):
        if len(g) < MIN_CD8T_TUMOR:
            continue
        rows_idx = adata.obs_names.get_indexer(g.index)
        X = adata.X[rows_idx][:, col_idx]
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        # Normalizzazione per profondita' di sequenziamento individuale della cellula (n_counts)
        n_counts = g.n_counts.values.reshape(-1, 1)
        norm = np.log1p(X / np.maximum(n_counts, 1) * 1e4)
        per_cell_score = norm.mean(axis=1)
        exh_scores[pat] = float(per_cell_score.mean())

    # Assemblaggio della tabella finale dei predittori a livello di paziente
    out = pd.DataFrame(
        {
            "patient": list(exh_scores.keys()),
            "exhaustion_score": list(exh_scores.values()),
        }
    ).set_index("patient")
    out["myeloid_fraction"] = myeloid_fraction.reindex(out.index)
    out["caf_fraction"] = caf_fraction.reindex(out.index)
    return out.reset_index()


# --------------------------------------------------------------------------- #
def main() -> None:
    """Orchestra la costruzione della tabella pseudobulk per paziente, l'esecuzione
    dei test di permutazione sulle tre ipotesi a priori e la PCA descrittiva."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument(
        "--arbitration", type=Path, default=Path("results/tables/arbitration_blood.csv")
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    # Caricamento dell'oggetto AnnData e conversione delle colonne categoriche in stringhe
    print("[..] carico h5ad ...")
    adata = ad.read_h5ad(args.h5ad)
    for c in ("tissue", "patients", "all_celltype"):
        adata.obs[c] = adata.obs[c].astype(str)

    # Costruzione e inner join tra outcome di errore e predittori microambientali per paziente
    outcome = build_outcome(args.arbitration)
    predictors = build_predictors(adata)

    df = outcome.merge(predictors, on="patient", how="inner")
    print(
        f"\n[ok] pazienti con outcome + predittori: {len(df)} -> {sorted(df.patient)}"
    )
    dropped_out = set(outcome.patient) - set(df.patient)
    dropped_pred = set(predictors.patient) - set(df.patient)
    if dropped_out:
        print(f"[i] esclusi (mancano predittori): {sorted(dropped_out)}")
    if dropped_pred:
        print(
            f"[i] esclusi (mancano outcome/troppo pochi cloni): {sorted(dropped_pred)}"
        )

    # Verifica della numerosita' minima richiesta per procedere con i test di correlazione
    if len(df) < 6:
        print(
            "\n[!!] troppo pochi pazienti per un test con un minimo di senso. Fermo qui."
        )
        return

    df.to_csv(TAB / "pseudobulk_patient_table.csv", index=False)

    print("\n" + "=" * 78)
    print(f"TEST DI CORRELAZIONE A LIVELLO DI PAZIENTE (n={len(df)})")
    print("=" * 78)
    print("Tre ipotesi, dichiarate PRIMA di guardare i dati. Permutation test,")
    print(f"{N_PERM} permutazioni. Correzione di Holm-Bonferroni sui 3 test.\n")

    # Valutazione di ciascuna ipotesi tramite correlazione di Spearman, permutazione e IC bootstrap
    results = {}
    for col, (desc, direction) in HYPOTHESES.items():
        x, y = df[col].values, df.error_rate.values
        rho, p, _ = permutation_corr_test(x, y)
        lo, hi = bootstrap_ci_corr(x, y)
        results[col] = {
            "desc": desc,
            "rho": rho,
            "p_perm": p,
            "ci_lo": lo,
            "ci_hi": hi,
            "direction_attesa": direction,
            "coerente_con_ipotesi": (rho > 0) == (direction == "positive"),
        }
        print(
            f"  {col:20s} rho={rho:+.3f}  p_perm={p:.4f}  IC95% boot=[{lo:+.3f}, {hi:+.3f}]"
        )
        print(f"      {desc}")
        print(
            f"      direzione osservata coerente con l'ipotesi: {results[col]['coerente_con_ipotesi']}\n"
        )

    # Correzione per molteplicita' secondo Holm-Bonferroni ed esportazione della tabella risultati
    pvals = {k: v["p_perm"] for k, v in results.items()}
    adj = holm_bonferroni(pvals)
    print("-" * 78)
    print("p-value corretti (Holm-Bonferroni, su 3 test):")
    for k, v in adj.items():
        sig = "SIGNIFICATIVO" if v < 0.05 else "non significativo"
        print(f"  {k:20s} p_adj = {v:.4f}   {sig}")

    out_tab = pd.DataFrame(results).T
    out_tab["p_holm"] = pd.Series(adj)
    out_tab.to_csv(TAB / "pseudobulk_hypothesis_tests.csv")

    print("\n" + "=" * 78)
    print(f"n={len(df)} pazienti: qualsiasi risultato qui e' IPOTESI-GENERANTE, non")
    print("confermativo. Un IC bootstrap su 10-14 punti e' largo per costruzione.")
    print("Un p<0.05 con n così piccolo va replicato prima di essere creduto;")
    print("un p>0.05 non e' evidenza di assenza di effetto, e' assenza di potenza.")
    print("=" * 78)

    # Analisi delle componenti principali (PCA) puramente descrittiva sui profili pseudobulk T tumorali
    print("\n[..] PCA pseudobulk (descrittiva, non un test) ...")
    tumor = adata.obs[
        (adata.obs.tissue == "Tumor") & (adata.obs.all_celltype.isin(["CD8T", "CD4T"]))
    ]
    pats = [p for p in df.patient if (tumor.patients == p).sum() >= MIN_CD8T_TUMOR]
    if len(pats) >= 5:
        sym = (  # noqa: F841
            adata.var.gene_symbol.astype(str).values
            if "gene_symbol" in adata.var
            else adata.var_names.values
        )
        pb = []
        for p in pats:
            # Estrazione e normalizzazione log1p(CPM/100) del profilo medio per singolo paziente
            idxp = adata.obs_names.get_indexer(tumor.index[tumor.patients == p])
            Xp = adata.X[idxp]
            nc = adata.obs.n_counts.values[idxp].reshape(-1, 1)
            norm = Xp.multiply(1e4 / np.maximum(nc, 1))
            mean_expr = np.asarray(norm.mean(axis=0)).ravel()
            pb.append(np.log1p(mean_expr))
        PB = np.vstack(pb)
        # Selezione dei 2.000 geni a maggiore varianza inter-paziente e proiezione sulle prime 2 PC
        var = PB.var(axis=0)
        top = np.argsort(-var)[:2000]
        pca = PCA(n_components=2, random_state=SEED)
        xy = pca.fit_transform(PB[:, top])

        errs = df.set_index("patient").loc[pats, "error_rate"].values
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        sc_ = ax.scatter(
            xy[:, 0], xy[:, 1], c=errs, cmap="viridis", s=140, edgecolor="k"
        )
        for i, p in enumerate(pats):
            ax.annotate(
                p,
                (xy[i, 0], xy[i, 1]),
                fontsize=8,
                xytext=(4, 4),
                textcoords="offset points",
            )
        fig.colorbar(sc_, ax=ax, label="tasso d'errore (Tumor)")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
        ax.set_title("Pseudobulk per paziente (CD8T+CD4T, tumore): analisi descrittiva")
        fig.tight_layout()
        fig.savefig(FIG / "pseudobulk_patient_pca.png", dpi=140)
        plt.close(fig)
        print(
            f"[ok] figura: results/figures/pseudobulk_patient_pca.png ({len(pats)} pazienti)"
        )
    else:
        print(f"[!] solo {len(pats)} pazienti con abbastanza cellule: salto la PCA.")

    print("\n" + "=" * 78)
    print("NOTA METODOLOGICA SULL'ESCLUSIONE DI SCANSIONI TRANSCRITTOMICHE")
    print("=" * 78)
    print("""
  Valutare la correlazione col tasso d'errore su tutti i ~20.000 geni per n=10-14
  pazienti e selezionare i p-value minori produrrebbe un'inflazione massiva di
  falsi positivi. Con 20.000 test indipendenti sotto ipotesi nulla, circa 1.000
  geni risulterebbero nominalmente "p<0.05" per puro effetto del campionamento,
  generando una lista priva di validita' statistica. Un'esplorazione genome-wide
  richiederebbe sia il controllo del False Discovery Rate (es. Benjamini-Hochberg)
  sia una potenza statistica adeguata (tipicamente n>50 campioni indipendenti per
  screening non guidati da ipotesi a priori), condizioni non compatibili con una
  coorte di 14 pazienti.
""")


if __name__ == "__main__":
    main()

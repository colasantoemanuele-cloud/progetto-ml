#!/usr/bin/env python3
"""
17_clonal_expansion.py: Associazione tra espansione clonale T e discordanza fenotipica intratumorale.

Razionale Biologico e Formulazione dell'Ipotesi
-----------------------------------------------
I cloni T maggiormente espansi (caratterizzati da un numero superiore di cellule che
condividono il medesimo recettore TCR alpha/beta) rappresentano popolazioni effettrici
sottoposte a stimolazione antigenica prolungata e differenziamento terminale. Poiche'
la transizione verso stati effettori/esausti attenua il profilo trascrizionale canonico
del lignaggio (convergenza CD8T -> CD4T/NK documentata nei moduli 08 e 09), si formula
l'ipotesi a priori che la dimensione clonale complessiva (sangue + tumore) correli
positivamente con la probabilita' di errore dell'annotazione trascrizionale nel
compartimento tumorale rispetto all'identita' di riferimento nel sangue periferico.

Integrazione dei Dati e Disegno Gerarchico dell'Inferenza
---------------------------------------------------------
Lo script effettua il join relazionale su `(patient, clone)` tra due tabelle gia'
validate nei moduli precedenti:
  - `arbitration_blood.csv` (da `10_loco.py`): esito binario di errore nel tumore
    rispetto al riferimento imparziale del sangue periferico (PBMC);
  - `pairwise_discordance_clone_b.csv` (da `07_robustness.py`): abbondanza clonale nei
    due compartimenti (`n_a` in PBMC, `n_b` in Tumor) ed eccesso continuo di discordanza.

Poiche' i cloni sono annidati all'interno dei pazienti, una correlazione calcolata sul
pool complessivo dei cloni assumendo indipendenza osservazionale sarebbe vulnerabile a
confondimento da livello paziente (ad esempio, variazioni sistematiche nella profondita'
di sequenziamento o nella cellularita' del campione). Per isolare rigorosamente
l'effetto clone-specifico dal confondente inter-paziente, l'analisi viene condotta su
due livelli complementari:
  1. Livello paziente (aggregato pseudobulk): correlazione di Spearman tra dimensione
     clonale media e tasso d'errore medio per paziente, valutata tramite test di
     permutazione esatto (immune da confondimento di cluster ma a potenza limitata da n=10);
  2. Livello clone (pooled stratificato): correlazione sul pool dei cloni valutata
     mediante Block Permutation Test (permutazione del vettore di esito esclusivamente
     all'interno dei confini di ciascun paziente, preservando i margini individuali) e
     intervalli di confidenza tramite Cluster Bootstrap sui pazienti interi.

Esecuzione:
    python src/17_clonal_expansion.py \
        --arbitration results/tables/arbitration_blood.csv \
        --pairwise results/tables/pairwise_discordance_clone_b.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Percorsi di output e costanti per i test di permutazione e filtro per paziente
FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_PERM = 20000
MIN_CLONES_PER_PATIENT = 5


# --------------------------------------------------------------------------- #
# Stimatori a livello di paziente (aggregazione pseudobulk per unita' clinica)
# --------------------------------------------------------------------------- #
def permutation_corr_test(
    x: np.ndarray, y: np.ndarray, n_perm: int = N_PERM, seed: int = SEED
):
    """Calcola la correlazione di Spearman e il p-value esatto a due code tramite
    permutazione vettorizzata dei ranghi sulle medie per paziente."""
    rng = np.random.default_rng(seed)
    n = len(x)
    # Conversione in ranghi e calcolo della statistica osservata
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rho_obs = float(np.corrcoef(rx, ry)[0, 1])
    # Generazione della distribuzione nulla di permutazione mediante argsort matriciale
    perm_idx = np.argsort(rng.random((n_perm, n)), axis=1)
    ry_perm = ry[perm_idx]
    rx_c = rx - rx.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    rhos = num / np.maximum(den, 1e-12)
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p)


def bootstrap_ci_corr(
    x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = SEED
):
    """Costruisce l'intervallo di confidenza al 95% per la correlazione di Spearman
    a livello di paziente mediante ricampionamento bootstrap con reinserimento."""
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        r = np.corrcoef(pd.Series(x[idx]).rank(), pd.Series(y[idx]).rank())[0, 1]
        if np.isfinite(r):
            vals.append(r)
    if len(vals) < 20:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# --------------------------------------------------------------------------- #
# Stimatori a livello di clone (pooled con struttura a blocchi per paziente)
# --------------------------------------------------------------------------- #
def block_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    group: np.ndarray,
    n_perm: int = N_PERM,
    seed: int = SEED,
):
    """Esegue un test di permutazione vincolato a blocchi (Block Permutation Test)
    per la correlazione di Spearman sul pool complessivo dei cloni.
    Permutando i valori dell'esito y esclusivamente all'interno di ciascun paziente
    (blocco), la distribuzione marginale di ogni paziente rimane invariata sotto
    l'ipotesi nulla, eliminando qualsiasi confondimento inter-paziente."""
    rng = np.random.default_rng(seed)
    # Ordinamento stabile delle osservazioni per paziente per indicizzare i blocchi contigui
    order = np.argsort(group, kind="stable")
    x_s, g_s = x[order], group[order]
    ry_s = pd.Series(y[order]).rank().values
    rx_s = pd.Series(x_s).rank().values
    rho_obs = float(np.corrcoef(rx_s, ry_s)[0, 1])

    # Costruzione della matrice di indici permutati blocco per blocco (intra-paziente)
    n = len(x_s)
    perm_idx = np.tile(np.arange(n), (n_perm, 1))
    uniq, starts = np.unique(g_s, return_index=True)
    bounds = list(starts) + [n]
    for i in range(len(uniq)):
        s, e = bounds[i], bounds[i + 1]
        size = e - s
        if size <= 1:
            continue
        # Permutazione indipendente ristretta all'intervallo [s, e) del singolo paziente
        keys = rng.random((n_perm, size))
        local = np.argsort(keys, axis=1)
        perm_idx[:, s:e] = s + local

    # Valutazione vettorizzata del coefficiente di correlazione sulle n_perm repliche a blocchi
    ry_perm = ry_s[perm_idx]
    rx_c = rx_s - rx_s.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    rhos = num / np.maximum(den, 1e-12)
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p)


def cluster_bootstrap_corr(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_col: str,
    n_boot: int = 2000,
    seed: int = SEED,
):
    """Calcola l'intervallo di confidenza al 95% della correlazione di Spearman
    pooled mediante Cluster Bootstrap, ricampionando con reinserimento interi
    pazienti anziche' singoli cloni per rispettare la dipendenza intra-cluster."""
    rng = np.random.default_rng(seed)
    pats = df[group_col].unique()
    vals = []
    for _ in range(n_boot):
        # Estrazione con reinserimento dei cluster paziente e concatenazione dei rispettivi cloni
        draw = rng.choice(pats, len(pats), replace=True)
        sub = pd.concat([df[df[group_col] == p] for p in draw])
        if sub[x_col].nunique() < 2 or sub[y_col].nunique() < 2:
            continue
        r = np.corrcoef(sub[x_col].rank(), sub[y_col].rank())[0, 1]
        if np.isfinite(r):
            vals.append(r)
    if len(vals) < 20:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# --------------------------------------------------------------------------- #
def run_analysis(df: pd.DataFrame, outcome_col: str, outcome_label: str) -> None:
    """Esegue e confronta l'inferenza a livello paziente, a livello clone con
    permutazione a blocchi e mediante test ingenuo non stratificato."""
    print(
        f"\n{'=' * 78}\nOUTCOME: {outcome_label}  (n_cloni={len(df)}, n_pazienti={df.patient.nunique()})\n{'=' * 78}"
    )

    # Analisi (1): aggregazione delle medie per paziente e test di permutazione esatto
    print(
        "\n--- (1) LIVELLO PAZIENTE: immune da confondimento di cluster, minore potenza ---"
    )
    agg = df.groupby("patient").agg(
        log_size=("log_size", "mean"),
        outcome=(outcome_col, "mean"),
        n=(outcome_col, "size"),
    )
    agg = agg[agg.n >= MIN_CLONES_PER_PATIENT]
    if len(agg) < 6:
        print(
            f"  solo {len(agg)} pazienti con >={MIN_CLONES_PER_PATIENT} cloni: salto."
        )
    else:
        rho, p = permutation_corr_test(agg.log_size.values, agg.outcome.values)
        lo, hi = bootstrap_ci_corr(agg.log_size.values, agg.outcome.values)
        print(
            f"  n={len(agg)} pazienti | rho={rho:+.3f}  p_perm={p:.4f}  IC95% boot=[{lo:+.3f}, {hi:+.3f}]"
        )
        print(f"  {'COERENTE con H1 (rho>0)' if rho > 0 else 'direzione OPPOSTA a H1'}")

    # Analisi (2): livello clone pooled con permutazione vincolata ai blocchi paziente
    print(
        "\n--- (2) LIVELLO CLONE, POOLED: maggiore potenza, corretto per paziente ---"
    )
    x, y, g = df.log_size.values, df[outcome_col].values, df.patient.values
    rho_p, p_p = block_permutation_test(x, y, g)
    lo_p, hi_p = cluster_bootstrap_corr(df, "log_size", outcome_col, "patient")
    print(
        f"  n={len(df)} cloni | rho={rho_p:+.3f}  p_block_perm={p_p:.4f}  "
        f"IC95% cluster-boot={[round(lo_p, 3), round(hi_p, 3)]}"
    )
    print(f"  {'COERENTE con H1 (rho>0)' if rho_p > 0 else 'direzione OPPOSTA a H1'}")

    # Confronto diagnostico con il test di permutazione non stratificato (che ignora i blocchi paziente)
    rng = np.random.default_rng(SEED)
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rho_naive = float(np.corrcoef(rx, ry)[0, 1])
    perm_idx = np.argsort(rng.random((N_PERM, len(x))), axis=1)
    ry_perm = ry[perm_idx]
    rx_c = rx - rx.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    rhos_naive = (ry_c * rx_c).sum(axis=1) / np.maximum(
        np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1)), 1e-12
    )
    p_naive = (np.sum(np.abs(rhos_naive) >= abs(rho_naive) - 1e-12) + 1) / (N_PERM + 1)
    print(
        f"\n  [confronto] permutation NAIVE (ignora i pazienti): rho={rho_naive:+.3f} p={p_naive:.4f}"
    )
    if p_naive < 0.05 <= p_p:
        print("  -> il naive trova un effetto che il block-permutation NON conferma:")
        print("     probabile CONFONDENTE DI PAZIENTE. Non fidarsi del naive.")
    elif p_p < 0.05:
        print(
            "  -> l'effetto sopravvive anche correggendo per paziente: piu' credibile."
        )


def make_figure(
    df: pd.DataFrame, outcome_col: str, outcome_label: str, fname: str
) -> None:
    """Genera il diagramma a dispersione dell'esito di discordanza in funzione
    della dimensione clonale log-trasformata, sovrapponendo i centroidi per paziente."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    # Visualizzazione dei singoli cloni (punti semitrasparenti raggruppati per paziente)
    for pat, g in df.groupby("patient"):
        ax.scatter(g.log_size, g[outcome_col], alpha=0.35, s=22, label=None)
    # Calcolo e sovrapposizione dei valori medi aggregati per ciascun paziente
    agg = df.groupby("patient").agg(
        log_size=("log_size", "mean"), outcome=(outcome_col, "mean")
    )
    ax.scatter(
        agg.log_size,
        agg.outcome,
        s=180,
        c="crimson",
        marker="D",
        edgecolor="k",
        zorder=5,
        label="media per paziente",
    )
    for pat, r in agg.iterrows():
        ax.annotate(
            pat,
            (r.log_size, r.outcome),
            fontsize=8,
            xytext=(5, 5),
            textcoords="offset points",
        )
    ax.set_xlabel("log(dimensione del clone, sangue+tumore)")
    ax.set_ylabel(outcome_label)
    ax.set_title(
        "Espansione clonale vs errore di annotazione (punti chiari = cloni singoli)"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=140)
    plt.close(fig)
    print(f"[ok] figura: results/figures/{fname}")


# --------------------------------------------------------------------------- #
def main() -> None:
    """Punto di ingresso dello script: unisce le tabelle di arbitrato ed espansione
    clonale ed esegue i test stratificati sia sull'errore binario sia sull'eccesso continuo."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--arbitration", type=Path, default=Path("results/tables/arbitration_blood.csv")
    )
    ap.add_argument(
        "--pairwise",
        type=Path,
        default=Path("results/tables/pairwise_discordance_clone_b.csv"),
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    for p in (args.arbitration, args.pairwise):
        if not p.exists():
            raise SystemExit(
                f"[ERRORE] {p} assente. Esegui prima 07_robustness.py e 10_loco.py."
            )

    # Estrazione degli esiti di errore nel tumore e delle abbondanze clonali appaiate PBMC-Tumor
    arb = pd.read_csv(args.arbitration)
    arb = arb[arb.tissue == "Tumor"][["patient", "clone", "error"]]
    pw = pd.read_csv(args.pairwise)
    pw = pw[(pw.comp_a == "PBMC") & (pw.comp_b == "Tumor")][
        ["patient", "clone", "n_a", "n_b", "excess"]
    ]

    # Inner join su (patient, clone) e calcolo della dimensione clonale totale log1p-trasformata
    df = arb.merge(pw, on=["patient", "clone"], how="inner")
    df["total_size"] = df.n_a + df.n_b
    df["log_size"] = np.log1p(df.total_size)
    print(
        f"[ok] join arbitration x pairwise: {len(df):,} cloni, {df.patient.nunique()} pazienti"
    )
    print(
        f"     dimensione clonale (sangue+tumore): mediana={df.total_size.median():.0f}, "
        f"range [{df.total_size.min()}, {df.total_size.max()}]"
    )
    df.to_csv(TAB / "clonal_expansion_joined.csv", index=False)

    print("\n" + "#" * 78)
    print("# IPOTESI H1 (unica, dichiarata prima di guardare i risultati):")
    print("# cloni T piu' espansi (sangue+tumore) -> piu' probabile che l'etichetta")
    print(
        "# nel tumore sia sbagliata rispetto all'identita' di riferimento dal sangue."
    )
    print("#" * 78)

    # Valutazione dell'ipotesi primaria (errore binario) e di robustezza (eccesso di discordanza)
    run_analysis(
        df, "error", "errore di annotazione nel tumore (0/1, riferimento sangue)"
    )
    make_figure(df, "error", "errore (0/1)", "clonal_expansion_vs_error.png")

    run_analysis(df, "excess", "eccesso di discordanza sangue-tumore (continuo, da 07)")
    make_figure(
        df, "excess", "eccesso di discordanza", "clonal_expansion_vs_excess.png"
    )

    print("\n" + "=" * 78)
    print("CRITERIO INTERPRETATIVO: esaminare prioritariamente la coerenza direzionale")
    print("al livello paziente (immune da confondimento ma a bassa potenza) insieme")
    print("al test pooled con permutazione a blocchi intra-paziente.")
    print("=" * 78)


if __name__ == "__main__":
    main()

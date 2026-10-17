#!/usr/bin/env python3
"""
17_clonal_expansion.py - Espansione clonale T contro errore di etichetta nel tumore.

Scopo
-----
Ipotesi dell'espansione clonale (qui "ipotesi E", distinta da H1-H4 di 16 e 19),
formulata prima di calcolare l'associazione: i cloni T piu' espansi (sangue + tumore)
hanno piu' probabilita' che l'etichetta nel tumore sia sbagliata rispetto all'identita'
stimata dal sangue. Due esiti per clone:
  - `error`  : errore 0/1 nel Tumor (schema `blood`, da `arbitration_blood.csv` di 10);
  - `excess` : eccesso di discordanza PBMC-Tumor (da `pairwise_discordance_clone_b.csv` di 07).
La covariata e' log1p(n_PBMC + n_Tumor) del clone.

I cloni sono annidati nei pazienti, quindi l'analisi ha due livelli:
  1. livello paziente: Spearman fra medie per paziente (pazienti con >= 5 cloni),
     permutazione Monte Carlo (20.000 permutazioni casuali, non esatta), IC bootstrap
     sui pazienti;
  2. livello clone: Spearman pooled su tutti i cloni con block permutation test
     (l'esito e' permutato solo entro paziente) e IC con cluster bootstrap sui pazienti.
     Il p e' condizionato al paziente (la distribuzione nulla conserva le differenze
     fra pazienti), ma il rho riportato e' pooled: contiene anche la parte di
     associazione fra pazienti. Per questo si riportano anche la media della
     distribuzione nulla e il p sulla statistica centrata su quella media.
Come confronto, una permutazione che ignora i pazienti (unita' sbagliata).

Analisi di sensibilita' del block permutation test sull'errore binario: 5 semi
diversi (20.000 permutazioni) e un seme con 100.000 permutazioni.

Input
-----
  --arbitration  results/tables/arbitration_blood.csv (da 10_loco.py)
  --pairwise     results/tables/pairwise_discordance_clone_b.csv (da 07_robustness.py)

Output
------
  results/tables/clonal_expansion_joined.csv            cloni con dimensione ed esiti
  results/tables/clonal_expansion_tests.csv             rho, p, p centrato, media nulla, IC per esito x livello
  results/tables/clonal_expansion_perm_sensitivity.csv  p del block permutation per seme e numero di permutazioni
  results/figures/clonal_expansion_vs_error.png, results/figures/clonal_expansion_vs_excess.png

Argomenti
---------
  --arbitration, --pairwise (default nei percorsi sopra)

Uso
---
    python src/17_clonal_expansion.py \
        --arbitration results/tables/arbitration_blood.csv \
        --pairwise results/tables/pairwise_discordance_clone_b.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

# Backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import (  # noqa: E402
    block_permutation_test,
    percentile_ci,
    permutation_corr_test,
)

FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_PERM = 20000
MIN_CLONES_PER_PATIENT = 5
SENS_SEEDS = (0, 1, 2, 3, 4)
SENS_N_PERM_LARGE = 100_000


def bootstrap_ci_corr(
    x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = SEED
):
    """IC 95% di Spearman fra medie per paziente, bootstrap sui pazienti ("expanded percentile")."""
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
    return percentile_ci(np.array(vals), n)


def cluster_bootstrap_corr(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_col: str,
    n_boot: int = 2000,
    seed: int = SEED,
):
    """IC 95% di Spearman pooled ricampionando interi pazienti ("expanded percentile")."""
    rng = np.random.default_rng(seed)
    pats = df[group_col].unique()
    parts = {p: df[df[group_col] == p] for p in pats}
    vals = []
    for _ in range(n_boot):
        draw = rng.choice(pats, len(pats), replace=True)
        sub = pd.concat([parts[p] for p in draw])
        if sub[x_col].nunique() < 2 or sub[y_col].nunique() < 2:
            continue
        r = np.corrcoef(sub[x_col].rank(), sub[y_col].rank())[0, 1]
        if np.isfinite(r):
            vals.append(r)
    if len(vals) < 20:
        return np.nan, np.nan
    return percentile_ci(np.array(vals), len(pats))


# --------------------------------------------------------------------------- #
def run_analysis(df: pd.DataFrame, outcome_col: str, outcome_label: str) -> list:
    """Analisi a livello paziente, a livello clone (block permutation) e permutazione ingenua."""
    print(
        f"\n{'=' * 78}\nESITO: {outcome_label}  (n_cloni={len(df)}, n_pazienti={df.patient.nunique()})\n{'=' * 78}"
    )
    rows = []

    print("\n--- (1) LIVELLO PAZIENTE: medie per paziente, potenza bassa ---")
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
        rho, p = permutation_corr_test(
            agg.log_size.values, agg.outcome.values, N_PERM, SEED
        )
        lo, hi = bootstrap_ci_corr(agg.log_size.values, agg.outcome.values)
        print(
            f"  n={len(agg)} pazienti | rho={rho:+.3f}  p_perm={p:.4f}  IC95% boot=[{lo:+.3f}, {hi:+.3f}]"
        )
        print(f"  {'direzione attesa (rho>0)' if rho > 0 else 'direzione opposta'}")
        rows.append(
            {
                "outcome": outcome_col,
                "level": "patient_means",
                "n": len(agg),
                "n_patients": len(agg),
                "rho": rho,
                "p": p,
                "ci_lo": lo,
                "ci_hi": hi,
            }
        )

    print("\n--- (2) LIVELLO CLONE, POOLED: permutazione entro paziente ---")
    x, y, g = df.log_size.values, df[outcome_col].values, df.patient.values
    bp = block_permutation_test(x, y, g, n_perm=N_PERM, seed=SEED)
    lo_p, hi_p = cluster_bootstrap_corr(df, "log_size", outcome_col, "patient")
    print(
        f"  n={len(df)} cloni | rho pooled={bp['rho']:+.3f}  p_block_perm={bp['p']:.4f}  "
        f"IC95% cluster-boot=[{lo_p:+.3f}, {hi_p:+.3f}]"
    )
    print(
        f"  media della distribuzione nulla = {bp['null_mean']:+.3f} | "
        f"p sulla statistica centrata = {bp['p_centered']:.4f}"
    )
    print(f"  {'direzione attesa (rho>0)' if bp['rho'] > 0 else 'direzione opposta'}")
    rows.append(
        {
            "outcome": outcome_col,
            "level": "clone_block_permutation",
            "n": len(df),
            "n_patients": df.patient.nunique(),
            "rho": bp["rho"],
            "p": bp["p"],
            "null_mean": bp["null_mean"],
            "p_centered": bp["p_centered"],
            "ci_lo": lo_p,
            "ci_hi": hi_p,
        }
    )

    # Confronto: permutazione che ignora i pazienti (unita' sbagliata)
    rho_naive, p_naive = permutation_corr_test(x, y, N_PERM, SEED)
    print(
        f"\n  [confronto] permutazione che ignora i pazienti: rho={rho_naive:+.3f} p={p_naive:.4f}"
    )
    rows.append(
        {
            "outcome": outcome_col,
            "level": "clone_naive_permutation",
            "n": len(df),
            "n_patients": df.patient.nunique(),
            "rho": rho_naive,
            "p": p_naive,
        }
    )
    if p_naive < 0.05 <= bp["p"]:
        print("  -> l'effetto trovato ignorando i pazienti non regge alla permutazione")
        print("     entro paziente: probabile confondente di paziente.")
    elif bp["p"] < 0.05:
        print("  -> l'associazione regge anche permutando solo entro paziente.")
    return rows


def perm_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    """Block permutation sull'errore binario con semi diversi e con piu' permutazioni."""
    x, y, g = df.log_size.values, df["error"].values, df.patient.values
    rows = []
    settings = [(s, N_PERM) for s in SENS_SEEDS] + [(SEED, SENS_N_PERM_LARGE)]
    print("\n--- sensibilita' del block permutation (errore binario) ---")
    for seed, n_perm in settings:
        bp = block_permutation_test(x, y, g, n_perm=n_perm, seed=seed)
        rows.append({"outcome": "error", "seed": seed, "n_perm": n_perm, **bp})
        print(
            f"  seme {seed} | {n_perm:>7,} permutazioni | p={bp['p']:.4f} "
            f"| p centrato={bp['p_centered']:.4f} | media nulla={bp['null_mean']:+.4f}"
        )
    return pd.DataFrame(rows)


def make_figure(
    df: pd.DataFrame, outcome_col: str, outcome_label: str, fname: str
) -> None:
    """Dispersione dell'esito contro log(dimensione clonale) con le medie per paziente."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    # Singoli cloni (punti semitrasparenti, un colore per paziente)
    for pat, g in df.groupby("patient"):
        ax.scatter(g.log_size, g[outcome_col], alpha=0.35, s=22, label=None)
    # Medie per paziente
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
    """Unisce arbitrato ed espansione clonale ed esegue i test su errore ed eccesso."""
    ap = argparse.ArgumentParser(
        description="Espansione clonale T contro errore di etichetta nel tumore."
    )
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

    # Esiti nel tumore e abbondanze clonali appaiate PBMC-Tumor
    arb = pd.read_csv(args.arbitration)
    arb = arb[arb.tissue == "Tumor"][["patient", "clone", "error"]]
    pw = pd.read_csv(args.pairwise)
    pw = pw[(pw.comp_a == "PBMC") & (pw.comp_b == "Tumor")][
        ["patient", "clone", "n_a", "n_b", "excess"]
    ]

    # Inner join su (patient, clone) e dimensione clonale totale log1p-trasformata
    df = arb.merge(pw, on=["patient", "clone"], how="inner")
    if df.empty or df.patient.nunique() < 2:
        raise SystemExit(
            "[ERRORE] join vuoto o con un solo paziente: nulla da testare."
        )
    df["error"] = df.error.astype(float)
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
    print(
        "# IPOTESI E (espansione clonale), formulata prima di calcolare l'associazione:"
    )
    print("# cloni piu' espansi (sangue+tumore) -> etichetta nel tumore piu' spesso")
    print("# diversa dall'identita' stimata dal sangue.")
    print("#" * 78)

    rows = run_analysis(
        df, "error", "errore di annotazione nel tumore (0/1, riferimento sangue)"
    )
    make_figure(df, "error", "errore (0/1)", "clonal_expansion_vs_error.png")
    rows += run_analysis(
        df, "excess", "eccesso di discordanza sangue-tumore (continuo, da 07)"
    )
    make_figure(
        df, "excess", "eccesso di discordanza", "clonal_expansion_vs_excess.png"
    )
    pd.DataFrame(rows).to_csv(TAB / "clonal_expansion_tests.csv", index=False)
    perm_sensitivity(df).to_csv(
        TAB / "clonal_expansion_perm_sensitivity.csv", index=False
    )

    print("\n" + "=" * 78)
    print("Lettura: guardare insieme la direzione a livello paziente (bassa potenza) e")
    print("il test pooled con permutazione entro paziente (p condizionato al paziente,")
    print("rho pooled: confrontarlo con la media della distribuzione nulla).")
    print("=" * 78)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
18_loyo_check.py - Sensibilita' leave-one-patient-out dell'associazione fra espansione clonale ed errore.

Scopo
-----
In 17_clonal_expansion.py l'associazione fra dimensione clonale ed errore (o eccesso
di discordanza) e' stimata su 287 cloni di 10 pazienti (versione pubblicata), con un
numero di cloni per paziente molto diseguale. Questo script verifica che il
risultato non dipenda da un solo paziente: per ciascun paziente rimuove tutti i suoi
cloni e ricalcola sui restanti il rho di Spearman pooled e il p del block
permutation test (permutazione Monte Carlo entro paziente, `scpipe.block_permutation_test`,
con gli stessi parametri di 17: 20.000 permutazioni, seme 0; sul campione completo
il p coincide quindi con `clonal_expansion_tests.csv`). Per ogni esclusione riporta
segno, variazione di rho, p non centrato, media della distribuzione nulla e p sulla
statistica centrata. Sensibilita' (errore binario): p di ogni esclusione con i semi
0-4 a 20.000 permutazioni e con 100.000 permutazioni. (I nomi dei file contengono
"loyo" per ragioni storiche: lo schema e' leave-one-patient-out.)

Input
-----
  --input  results/tables/clonal_expansion_joined.csv (da 17_clonal_expansion.py)

Output
------
  results/tables/loyo_clonal_expansion_error.csv    una riga per paziente escluso (esito errore 0/1):
                                                    rho, p, media nulla, p centrato; rho e p sul
                                                    campione completo come riferimento
  results/tables/loyo_clonal_expansion_excess.csv   idem per l'eccesso di discordanza
  results/tables/loyo_clonal_expansion_perm_sensitivity.csv
                                                    errore 0/1: per paziente escluso, p con semi 0-4
                                                    (20.000 permutazioni) e con 100.000 permutazioni
  results/figures/loyo_clonal_expansion_error.png, results/figures/loyo_clonal_expansion_excess.png

Argomenti
---------
  --input (default results/tables/clonal_expansion_joined.csv)

Uso
---
    python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scpipe import block_permutation_test as _bpt  # noqa: E402

# Percorsi di output e parametri del test di permutazione a blocchi
FIG = Path("results/figures")
TAB = Path("results/tables")
# Stessi parametri di 17_clonal_expansion.py (N_PERM, SEED, semi e permutazioni della sensibilita')
SEED = 0
N_PERM = 20000
SENS_SEEDS = (0, 1, 2, 3, 4)
SENS_N_PERM_LARGE = 100_000


def block_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    group: np.ndarray,
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> dict:
    """Block permutation test (`scpipe.block_permutation_test`): rho, p, null_mean, p_centered."""
    return _bpt(x, y, group, n_perm=n_perm, seed=seed)


def loyo(df: pd.DataFrame, outcome_col: str, label: str) -> pd.DataFrame:
    """Rimuove un paziente alla volta e ricalcola rho pooled e p sui pazienti restanti."""
    print(f"\n{'=' * 78}\n{label}\n{'=' * 78}")
    pats = sorted(df.patient.unique())

    # Valore di riferimento su tutti i pazienti (stessi parametri di 17)
    full = block_permutation_test(
        df.log_size.values, df[outcome_col].values, df.patient.values
    )
    rho_full, p_full = full["rho"], full["p"]
    print(
        f"tutti i {len(pats)} pazienti:  rho={rho_full:+.3f}  p={p_full:.4f}  "
        f"media nulla={full['null_mean']:+.3f}  p centrato={full['p_centered']:.4f}\n"
    )

    rows = []
    for held_out in pats:
        # Esclusione di tutti i cloni associati al paziente held_out
        sub = df[df.patient != held_out]
        bp = block_permutation_test(
            sub.log_size.values, sub[outcome_col].values, sub.patient.values
        )
        rho, p = bp["rho"], bp["p"]
        n_clones_removed = (df.patient == held_out).sum()
        rows.append(
            {
                "patient_removed": held_out,
                "n_clones_removed": n_clones_removed,
                "rho_without": rho,
                "p_without": p,
                "null_mean_without": bp["null_mean"],
                "p_centered_without": bp["p_centered"],
                "delta_rho": rho - rho_full,
                # NaN (esito costante senza il paziente) non conta come significativo
                "still_significant": bool(p < 0.05),
                "rho_full": rho_full,
                "p_full": p_full,
                "null_mean_full": full["null_mean"],
                "p_centered_full": full["p_centered"],
            }
        )
        flag = "" if p < 0.05 else "  [non piu' significativo senza questo paziente]"
        print(
            f"  senza {held_out} ({n_clones_removed:3d} cloni): rho={rho:+.3f}  p={p:.4f}{flag}"
        )

    # Sintesi su tutte le esclusioni
    out = pd.DataFrame(rows)
    n_robust_c = int((out.p_centered_without < 0.05).sum())
    print(
        f"\ncon il p centrato: p < 0.05 in {n_robust_c}/{len(pats)} versioni leave-one-out."
    )
    n_robust = out.still_significant.sum()
    print(
        f"\nrisultato resta significativo (p<0.05) in {n_robust}/{len(pats)} versioni leave-one-out."
    )
    if out.p_without.isna().any():
        print(
            f"[!!] p non definito in {int(out.p_without.isna().sum())} esclusioni (esito costante)"
        )
    if n_robust == len(pats):
        print("VERDETTO: nessuna esclusione di un singolo paziente porta p sopra 0.05.")
    elif n_robust >= len(pats) - 1:
        print(
            "VERDETTO: p < 0.05 tranne che rimuovendo un paziente specifico (da ispezionare)."
        )
    else:
        print(
            "VERDETTO: il segnale dipende da pochi pazienti specifici. Non generalizzabile"
        )
        print("          senza cautela: da riportare come tale.")

    # Identificazione del paziente la cui esclusione determina la massima riduzione di rho
    if out.rho_without.isna().all():
        print(
            "\nrho non definito in tutte le esclusioni (esito costante): nessun confronto."
        )
        return out
    worst = out.loc[out.rho_without.idxmin()]
    print(
        f"\npaziente il cui contributo pesa di piu' (rho piu' basso se rimosso): "
        f"{worst.patient_removed} (rho scende a {worst.rho_without:+.3f})"
    )
    return out


def perm_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    """Errore binario: p di ogni esclusione con semi diversi e con piu' permutazioni."""
    settings = [(s, N_PERM) for s in SENS_SEEDS] + [(SEED, SENS_N_PERM_LARGE)]
    print(
        "\n--- sensibilita' del block permutation (errore binario, leave-one-out) ---"
    )
    rows = []
    for held_out in sorted(df.patient.unique()):
        sub = df[df.patient != held_out]
        x, y, g = sub.log_size.values, sub["error"].values, sub.patient.values
        for seed, n_perm in settings:
            bp = block_permutation_test(x, y, g, n_perm=n_perm, seed=seed)
            rows.append(
                {"patient_removed": held_out, "seed": seed, "n_perm": n_perm, **bp}
            )
    out = pd.DataFrame(rows)
    rng = out.groupby("patient_removed").p.agg(["min", "max"])
    print(
        f"  p fra semi e numero di permutazioni: da {out.p.min():.4f} a {out.p.max():.4f}; "
        f"p < 0.05 in {int((rng['max'] < 0.05).sum())}/{len(rng)} esclusioni con tutte le impostazioni"
    )
    return out


def plot_loyo(out: pd.DataFrame, rho_full: float, fname: str, title: str) -> None:
    """Rappresenta graficamente la stabilita' dei coefficienti di correlazione
    ottenuti escludendo ciascun paziente rispetto alla stima sulla coorte completa."""
    fig, ax = plt.subplots(figsize=(7.5, 4))
    order = out.sort_values("rho_without")
    colors = ["#d95f02" if not s else "#1b9e77" for s in order.still_significant]
    ax.barh(order.patient_removed, order.rho_without, color=colors, alpha=0.85)
    # Linea verticale tratteggiata (parametro grafico ls="--") per la correlazione di riferimento
    ax.axvline(
        rho_full,
        ls="--",
        c="k",
        lw=1.2,
        label=f"rho con tutti i pazienti ({rho_full:+.3f})",
    )
    ax.axvline(0, ls="-", c="gray", lw=0.6)
    ax.set_xlabel("rho pooled, rimosso quel paziente")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=140)
    plt.close(fig)
    print(f"[ok] figura: results/figures/{fname}")


def main() -> None:
    """Orchestra il test di sensibilita' Leave-One-Patient-Out per entrambi gli
    indicatori di discordanza clonale (errore binario ed eccesso continuo)."""
    ap = argparse.ArgumentParser(
        description="Sensibilita' leave-one-patient-out dell'associazione espansione-errore."
    )
    ap.add_argument(
        "--input", type=Path, default=Path("results/tables/clonal_expansion_joined.csv")
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        raise SystemExit(
            f"[ERRORE] {args.input} assente. Esegui prima 17_clonal_expansion.py."
        )
    df = pd.read_csv(args.input)
    if df.patient.nunique() < 3:
        raise SystemExit(
            "[ERRORE] servono almeno 3 pazienti per il leave-one-patient-out."
        )
    df["error"] = df.error.astype(float)
    print(
        f"[ok] {len(df):,} cloni, {df.patient.nunique()} pazienti, da {args.input.name}"
    )

    # Iterazione sui due outcome: errore binario (riferimento PBMC) ed eccesso di discordanza
    for outcome_col, label, fname_suffix in [
        ("error", "leave-one-patient-out: errore di annotazione (0/1)", "error"),
        (
            "excess",
            "leave-one-patient-out: eccesso di discordanza (continuo)",
            "excess",
        ),
    ]:
        out = loyo(df, outcome_col, label)
        rho_full = float(out.rho_full.iloc[0])
        out.to_csv(TAB / f"loyo_clonal_expansion_{fname_suffix}.csv", index=False)
        plot_loyo(
            out,
            rho_full,
            f"loyo_clonal_expansion_{fname_suffix}.png",
            label,
        )

    perm_sensitivity(df).to_csv(
        TAB / "loyo_clonal_expansion_perm_sensitivity.csv", index=False
    )

    print("\n" + "=" * 78)
    print("Se entrambi gli outcome (errore ed eccesso) restano con p < 0.05 in tutte")
    print("le versioni leave-one-out, il risultato di 17 non dipende da un singolo")
    print("paziente (resta l'associazione fra 10 pazienti, non una generalizzazione).")
    print("Se non lo sono, questo costituisce il risultato da riportare.")
    print("=" * 78)


if __name__ == "__main__":
    main()

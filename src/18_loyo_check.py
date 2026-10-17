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
permutation test (permutazione Monte Carlo entro paziente, 10.000 permutazioni,
`scpipe.block_permutation_test`). Riporta segno, variazione di rho e p < 0.05 per
ogni esclusione. (I nomi dei file contengono "loyo" per ragioni storiche: lo schema
e' leave-one-patient-out.)

Input
-----
  --input  results/tables/clonal_expansion_joined.csv (da 17_clonal_expansion.py)

Output
------
  results/tables/loyo_clonal_expansion_error.csv    una riga per paziente escluso (esito errore 0/1),
                                                    con rho e p sul campione completo come riferimento
  results/tables/loyo_clonal_expansion_excess.csv   idem per l'eccesso di discordanza
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
SEED = 0
N_PERM = 10000


def block_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    group: np.ndarray,
    n_perm: int = N_PERM,
    seed: int = SEED,
):
    """Spearman pooled e p del block permutation test (`scpipe.block_permutation_test`)."""
    r = _bpt(x, y, group, n_perm=n_perm, seed=seed)
    return r["rho"], r["p"]


def loyo(df: pd.DataFrame, outcome_col: str, label: str) -> pd.DataFrame:
    """Rimuove un paziente alla volta e ricalcola rho pooled e p sui pazienti restanti."""
    print(f"\n{'=' * 78}\n{label}\n{'=' * 78}")
    pats = sorted(df.patient.unique())

    # Valore di riferimento su tutti i pazienti
    rho_full, p_full = block_permutation_test(
        df.log_size.values, df[outcome_col].values, df.patient.values
    )
    print(
        f"tutti i {len(pats)} pazienti:  rho={rho_full:+.3f}  p={p_full:.4f}   (riferimento, da 17)\n"
    )

    rows = []
    for held_out in pats:
        # Esclusione di tutti i cloni associati al paziente held_out
        sub = df[df.patient != held_out]
        rho, p = block_permutation_test(
            sub.log_size.values, sub[outcome_col].values, sub.patient.values
        )
        n_clones_removed = (df.patient == held_out).sum()
        rows.append(
            {
                "patient_removed": held_out,
                "n_clones_removed": n_clones_removed,
                "rho_without": rho,
                "p_without": p,
                "delta_rho": rho - rho_full,
                "still_significant": p < 0.05,
                "rho_full": rho_full,
                "p_full": p_full,
            }
        )
        flag = "" if p < 0.05 else "  [non piu' significativo senza questo paziente]"
        print(
            f"  senza {held_out} ({n_clones_removed:3d} cloni): rho={rho:+.3f}  p={p:.4f}{flag}"
        )

    # Sintesi su tutte le esclusioni
    out = pd.DataFrame(rows)
    n_robust = out.still_significant.sum()
    print(
        f"\nrisultato resta significativo (p<0.05) in {n_robust}/{len(pats)} versioni leave-one-out."
    )
    if n_robust == len(pats):
        print("VERDETTO: nessuna esclusione di un singolo paziente porta p sopra 0.05.")
    elif n_robust >= len(pats) - 1:
        print(
            "VERDETTO: robusto tranne che rimuovendo un paziente specifico (da ispezionare),"
        )
        print("          ma non e' un artefatto di un solo campione.")
    else:
        print(
            "VERDETTO: il segnale dipende da pochi pazienti specifici. Non generalizzabile"
        )
        print("          senza cautela: da riportare come tale.")

    # Identificazione del paziente la cui esclusione determina la massima riduzione di rho
    worst = out.loc[out.rho_without.idxmin()]
    print(
        f"\npaziente il cui contributo pesa di piu' (rho piu' basso se rimosso): "
        f"{worst.patient_removed} (rho scende a {worst.rho_without:+.3f})"
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
        rho_full, _ = block_permutation_test(
            df.log_size.values, df[outcome_col].values, df.patient.values
        )
        out = loyo(df, outcome_col, label)
        out.to_csv(TAB / f"loyo_clonal_expansion_{fname_suffix}.csv", index=False)
        plot_loyo(
            out,
            rho_full,
            f"loyo_clonal_expansion_{fname_suffix}.png",
            label,
        )

    print("\n" + "=" * 78)
    print("Se entrambi gli outcome (errore ed eccesso) restano significativi in")
    print("tutte le versioni leave-one-out, il risultato di 17 e' da considerare")
    print("solido: non e' l'artefatto di un paziente con molti cloni o un pattern")
    print("estremo. Se non lo sono, questo costituisce il risultato da riportare.")
    print("=" * 78)


if __name__ == "__main__":
    main()

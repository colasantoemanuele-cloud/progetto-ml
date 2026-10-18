#!/usr/bin/env python3
"""
21_fig_blood_scheme.py - Tabella e figura del tasso d'errore nello schema primario (riferimento = solo PBMC).

Scopo
-----
Riassume l'errore delle etichette T in `Adjacent_normal` e `Tumor` quando l'identita'
clonale di riferimento e' stimata dal solo sangue periferico (schema `blood` di
10_loco.py). Per avere esattamente gli stessi numeri di 10_loco.py importa quel
modulo e ne usa la funzione `report()` (cluster bootstrap sui pazienti, intervallo
"expanded percentile") sulla tabella gia' calcolata `arbitration_blood.csv`.

Input
-----
  --arbitration  results/tables/arbitration_blood.csv (da 10_loco.py)

Output
------
  results/tables/error_rate_blood.csv   errore per compartimento: n cloni, n pazienti, errore, IC;
                                        con <= 5 pazienti anche il numero di ricampionature distinte
  results/figures/error_rate_blood.png  barre con IC; con <= 5 pazienti grigie tratteggiate
                                        (IC non calibrato)

Argomenti
---------
  --arbitration (default results/tables/arbitration_blood.csv)

Uso
---
    python src/21_fig_blood_scheme.py
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

# Backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
TAB, FIG = Path("results/tables"), Path("results/figures")


def load_loco():
    """Carica 10_loco.py come modulo (il nome inizia con una cifra e non e' importabile)."""
    spec = importlib.util.spec_from_file_location("loco", HERE / "10_loco.py")
    loco = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loco)
    return loco


def main() -> None:
    """Calcola la tabella dello schema blood e disegna la figura."""
    ap = argparse.ArgumentParser(
        description="Tabella e figura dell'errore con riferimento solo PBMC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--arbitration",
        type=Path,
        default=TAB / "arbitration_blood.csv",
        help="identita' clonale e errore per (clone, compartimento), schema 'blood' (da 10_loco.py)",
    )
    args = ap.parse_args()
    if not args.arbitration.exists():
        raise SystemExit(
            f"[ERRORE] {args.arbitration} assente: esegui prima 10_loco.py"
        )

    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    loco = load_loco()
    j = pd.read_csv(args.arbitration)
    r = loco.report(j, "BLOOD (riferimento = solo PBMC)")
    if r.empty:
        raise SystemExit(
            "Nessun compartimento valutabile: controlla arbitration_blood.csv"
        )

    r.to_csv(TAB / "error_rate_blood.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    # barre con <= SMALL_K pazienti grigie tratteggiate: IC non calibrato
    loco.bar_with_calibration(ax, r.tissue, r.error, r.lo, r.hi, r.n_pat)
    for i, row in enumerate(r.itertuples()):
        ax.text(
            i,
            row.hi + 0.008,
            f"{row.error:.3f}\n({row.n_pat} pz)",
            ha="center",
            fontsize=9,
        )
    ax.set_ylabel("tasso d'errore dell'etichetta")
    ax.set_title("Schema primario: riferimento = solo sangue (PBMC)")
    ax.set_ylim(0, max(r.hi) * 1.3)
    fig.tight_layout()
    fig.savefig(FIG / "error_rate_blood.png", dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: {FIG / 'error_rate_blood.png'}")
    print(f"[ok] tabella: {TAB / 'error_rate_blood.csv'}")


if __name__ == "__main__":
    main()

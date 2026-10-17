#!/usr/bin/env python3
"""
sim_d_xai_stability.py - Simulazione (d): il confronto LOPO contro controllo di 20 con e senza effetto paziente.

Scopo
-----
Applica la procedura di 20_xai_per_patient.py (stesse funzioni `lopo_explanations`,
`control_explanations`, `fit_explanation`: logistica L1, top-50 geni) a conteggi
simulati con 5 pazienti, e misura la differenza Jaccard(controllo) - Jaccard(LOPO):
  - con effetto paziente: ogni paziente ha 25 geni propri alzati di 3 volte nelle sue
    cellule tumorali (programma tumorale specifico del paziente) e 40 geni con un
    livello di base proprio (impronta del paziente, uguale nei due tessuti);
  - senza effetto paziente: nessun gene specifico e nessuna impronta.
In entrambi i casi 25 geni condivisi da tutti i pazienti sono alzati di 2 volte nel
tumore. 5 pazienti x 600 cellule (meta' Tumor, meta' Adjacent_normal), 400 geni,
conteggi Poisson con profondita' variabile per cellula. 5 repliche per condizione,
semi fissi. La differenza attesa e' piu' grande con effetto paziente; la simulazione
misura di quanto.

Input
-----
  nessuno (dati simulati); importa src/20_xai_per_patient.py

Output
------
  results/tables/sim_d_xai_stability.csv  condizione x replica: Jaccard LOPO, controllo, differenza
                                          (e una riga di media per condizione)

Argomenti
---------
  --reps (default 5)

Uso
---
    python src/sim/sim_d_xai_stability.py
"""

from __future__ import annotations

import argparse
import importlib.util
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
TAB = Path("results/tables")
N_PAT, N_CELLS, N_GENES = 5, 600, 400
N_SHARED, N_OWN, N_FP = 25, 25, 40


def load_xai():
    """Carica src/20_xai_per_patient.py come modulo."""
    spec = importlib.util.spec_from_file_location(
        "xai", HERE.parent / "20_xai_per_patient.py"
    )
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def simulate(rng, patient_effect: bool):
    base = rng.gamma(0.8, 0.5, N_GENES)
    shared = rng.choice(N_GENES, N_SHARED, replace=False)
    X, y, pats = [], [], []
    for p in range(N_PAT):
        lam_n = base.copy()
        if patient_effect:
            fp = rng.choice(N_GENES, N_FP, replace=False)
            lam_n[fp] *= rng.gamma(4.0, 0.5, N_FP)  # impronta del paziente
        lam_t = lam_n.copy()
        lam_t[shared] *= 2.0
        if patient_effect:
            own = rng.choice(
                np.setdiff1d(np.arange(N_GENES), shared), N_OWN, replace=False
            )
            lam_t[own] *= 3.0
        depth = rng.gamma(5.0, 0.2, N_CELLS)[:, None]
        half = N_CELLS // 2
        lam = np.vstack(
            [np.tile(lam_t, (half, 1)), np.tile(lam_n, (N_CELLS - half, 1))]
        )
        X.append(rng.poisson(lam * depth))
        y += ["Tumor"] * half + ["Adjacent_normal"] * (N_CELLS - half)
        pats += [f"P{p}"] * N_CELLS
    X = sp.csr_matrix(np.vstack(X).astype(np.float32))
    pats = np.array(pats)
    folds = np.array([int(p[1:]) for p in pats])
    return X, np.array(y), pats, folds


def run(reps: int) -> pd.DataFrame:
    xai = load_xai()
    genes = np.array([f"G{i}" for i in range(N_GENES)])
    rows = []
    for cond, eff in (("patient_effect", True), ("no_patient_effect", False)):
        for r in range(reps):
            rng = np.random.default_rng(500 + r + (0 if eff else 1000))
            X, y, pats, folds = simulate(rng, eff)
            xrng = np.random.default_rng(0)
            _, tops, _ = xai.lopo_explanations(X, y, folds, pats, genes, xrng)
            fs = sorted(tops)
            j_lopo = np.mean(
                [xai.jaccard(tops[a], tops[b]) for a, b in combinations(fs, 2)]
            )
            n_train = int(np.median([np.sum(folds != f) for f in fs]))
            ctrl = xai.control_explanations(X, y, n_train, genes, xrng)
            j_ctrl = np.mean([xai.jaccard(a, b) for a, b in combinations(ctrl, 2)])
            rows.append(
                {
                    "condition": cond,
                    "rep": r,
                    "jaccard_lopo": j_lopo,
                    "jaccard_control": j_ctrl,
                    "difference": j_ctrl - j_lopo,
                }
            )
            print(rows[-1])
    d = pd.DataFrame(rows)
    mean = d.groupby("condition")[
        ["jaccard_lopo", "jaccard_control", "difference"]
    ].mean()
    mean = mean.reset_index().assign(rep="mean")
    return pd.concat([d, mean], ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Simulazione (d): stabilita' XAI LOPO contro controllo, con e senza effetto paziente."
    )
    ap.add_argument("--reps", type=int, default=5)
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)
    out = run(args.reps)
    out.to_csv(TAB / "sim_d_xai_stability.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

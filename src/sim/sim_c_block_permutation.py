#!/usr/bin/env python3
"""
sim_c_block_permutation.py - Simulazione (c): falsi positivi del block permutation test con confondente di paziente.

Scopo
-----
Verifica `scpipe.block_permutation_test` (usato da 17 e 18) sotto H0 con un
confondente di paziente: un effetto di paziente z ~ N(0, 1) aumenta sia la
dimensione clonale (x = c * z + N(0, 1), log dimensione) sia la probabilita'
d'errore (y ~ Bernoulli(logit^-1(-1.5 + c * z))), ma entro paziente x e y sono
indipendenti. 10 pazienti con il numero di cloni dei dati (5-59 per paziente,
287 in totale). Per ogni replica: p non centrato (usato in 17/18), p centrato sulla
media della distribuzione nulla e, per confronto, la permutazione che ignora i
pazienti. Frazione di p < 0.05 = tasso di falsi positivi.
Forza del confondente c in {0, 0.5, 1, 2}; 1.000 repliche per livello, 2.000 permutazioni
per replica; intervallo binomiale esatto (Clopper-Pearson) al 95% per ogni tasso. Le prime
100 repliche di ogni livello sono quelle della versione precedente di questo script
(100 repliche, semi 300 + i), riportate a parte; il dato "2 su 100" citato nella
presentazione non e' riprodotto.

Input
-----
  nessuno (dati simulati, semi fissi)

Output
------
  results/tables/sim_c_block_permutation.csv  c, repliche, falsi positivi (conteggio, frazione, IC binomiale)
                                               per metodo, falsi positivi nelle prime 100 repliche

Argomenti
---------
  --n-perm (default 2000), --reps (repliche per livello, default 1000)

Uso
---
    python src/sim/sim_c_block_permutation.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scpipe import block_permutation_test, permutation_corr_test  # noqa: E402

TAB = Path("results/tables")
SIZES = np.array([17, 26, 47, 42, 22, 10, 30, 59, 29, 5])


def run(conf: float, reps: int, n_perm: int, seed: int) -> dict:
    """Falsi positivi (p < 0.05) su `reps` repliche sotto H0 con confondente di forza `conf`."""
    rng = np.random.default_rng(seed)
    gi = np.repeat(np.arange(len(SIZES)), SIZES)
    g = np.array([f"P{i:02d}" for i in gi])
    fp = {"block": 0, "block_centered": 0, "naive": 0}
    fp100 = dict(fp)
    for r in range(reps):
        z = rng.normal(0, 1, len(SIZES))
        x = conf * z[gi] + rng.normal(0, 1, len(gi))
        pr = 1 / (1 + np.exp(-(-1.5 + conf * z[gi])))
        y = (rng.random(len(gi)) < pr).astype(float)
        bp = block_permutation_test(x, y, g, n_perm=n_perm, seed=r)
        _, p_naive = permutation_corr_test(x, y, n_perm=n_perm, seed=r)
        fp["block"] += bp["p"] < 0.05
        fp["block_centered"] += bp["p_centered"] < 0.05
        fp["naive"] += p_naive < 0.05
        if r == 99:
            fp100 = dict(fp)
    out = {"confounder": conf, "reps": reps}
    for k, v in fp.items():
        ci = binomtest(int(v), reps).proportion_ci(
            confidence_level=0.95, method="exact"
        )
        out |= {
            f"fp_{k}": int(v),
            f"rate_{k}": v / reps,
            f"lo_{k}": ci.low,
            f"hi_{k}": ci.high,
        }
    return out | {f"fp_first100_{k}": int(v) for k, v in fp100.items()}


def main() -> None:
    """Esegue la simulazione e scrive la tabella."""
    ap = argparse.ArgumentParser(
        description="Simulazione (c): block permutation test sotto H0 con confondente di paziente."
    )
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--reps", type=int, default=1000)
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)
    rows = [
        run(c, args.reps, args.n_perm, seed=300 + i)
        for i, c in enumerate((0, 0.5, 1.0, 2.0))
    ]
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "sim_c_block_permutation.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

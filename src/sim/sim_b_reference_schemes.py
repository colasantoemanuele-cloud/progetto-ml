#!/usr/bin/env python3
"""
sim_b_reference_schemes.py - Simulazione (b): tre schemi di riferimento clonale a confronto.

Scopo
-----
Confronta su dati simulati i tre schemi di 10_loco.py (`naive`, `loco`, `blood`),
usando le funzioni di quel modulo (`per_comp_stats`, `arbitrate`), in due condizioni:
  - H0: nessun artefatto; l'errore di etichetta e' lo stesso in tutti i compartimenti:
    ogni (clone, compartimento) e' "mal etichettato" con probabilita' 0.10, e in quel
    caso l'80% delle sue cellule riceve l'etichetta sbagliata (l'annotazione per
    clustering sbaglia gruppi di cellule, non cellule isolate); nelle altre il 3%;
  - artefatto atteso nel tumore, in due intensita': perdita tecnica del segnale CD8
    (rilevamento di CD8A e CD8B moltiplicato per 0.4 "moderato" o 0.15 "forte") ed
    espressione ambientale di CD4 (CD4 > 0 nel 35% o nel 60% delle cellule tumorali di
    qualsiasi lignaggio). Le etichette hanno lo stesso errore vero della condizione H0:
    cambiano solo i marcatori.
Uno schema corretto dovrebbe stimare circa l'errore vero in ogni compartimento che
giudica. Lo schema `loco`, per i cloni presenti solo in PBMC e Tumor, giudica il sangue
con il tumore come riferimento: con l'artefatto attribuisce al sangue errori che non ci
sono.

Disegno: 10 pazienti x 40 cloni; meta' dei pazienti ha anche Adjacent_normal e in loro
il 30% dei cloni e' presente in tutti e tre i compartimenti (gli altri solo in PBMC e
Tumor); lignaggio CD8 o CD4 con probabilita' 0.5; 3 + Poisson(4) cellule per clone e
compartimento. Rilevamento per cellula: lignaggio CD8: CD8A 0.75, CD8B 0.55, CD4 0.03;
lignaggio CD4: CD8A 0.03, CD8B 0.02, CD4 0.45. Etichetta sbagliata: CD4 -> CD8T,
CD8 -> CD4T o NK con uguale probabilita'. 30 repliche per condizione, semi fissi.
I parametri sono scelti a priori (non tarati per riprodurre un numero atteso).

Limite: dato il lignaggio, etichette e marcatori sono simulati indipendenti (l'errore
di etichetta non dipende dai marcatori osservati). Nei dati reali l'etichetta del
clustering e i marcatori usati per l'identita' vengono dalle stesse cellule: la
circolarita' dello schema `naive` (il compartimento giudicato entra nel proprio
riferimento) qui NON e' simulata, quindi il basso bias di `naive` in questa
simulazione non dice che `naive` sia corretto sui dati reali.

Input
-----
  nessuno (dati simulati); importa src/10_loco.py

Output
------
  results/tables/sim_b_reference_schemes.csv  condizione x schema x compartimento:
                                              errore stimato medio, errore vero medio, repliche

Argomenti
---------
  --reps (default 30)

Uso
---
    python src/sim/sim_b_reference_schemes.py
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TAB = Path("results/tables")
N_PAT, N_CLONES = 10, 40
P_TRIPLE = 0.30
P_BAD_COMP, BAD_CELL_ERR, GOOD_CELL_ERR = 0.10, 0.80, 0.03
ARTIFACTS = {
    "H0": None,
    "artifact_moderate": (0.4, 0.35),
    "artifact_strong": (0.15, 0.60),
}
DET = {
    "CD8": {"CD8A": 0.75, "CD8B": 0.55, "CD4": 0.03},
    "CD4": {"CD8A": 0.03, "CD8B": 0.02, "CD4": 0.45},
}


def load_loco():
    """Carica src/10_loco.py come modulo."""
    spec = importlib.util.spec_from_file_location("loco", HERE.parent / "10_loco.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def simulate(rng, artifact) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cellule simulate (colonne come in 10_loco) e verita' per (clone, compartimento)."""
    cells, truth = [], []
    for p in range(N_PAT):
        has_adj = p < N_PAT // 2
        for c in range(N_CLONES):
            lin = "CD8" if rng.random() < 0.5 else "CD4"
            comps = ["PBMC", "Tumor"]
            if has_adj and rng.random() < P_TRIPLE:
                comps = ["PBMC", "Adjacent_normal", "Tumor"]
            for t in comps:
                n = 3 + rng.poisson(4)
                det = dict(DET[lin])
                if artifact is not None and t == "Tumor":
                    factor, amb = artifact
                    det["CD8A"] *= factor
                    det["CD8B"] *= factor
                    det["CD4"] = max(det["CD4"], amb)
                cd8a = rng.random(n) < det["CD8A"]
                cd8b = rng.random(n) < det["CD8B"]
                cd4 = rng.random(n) < det["CD4"]
                right = "CD8T" if lin == "CD8" else "CD4T"
                e = BAD_CELL_ERR if rng.random() < P_BAD_COMP else GOOD_CELL_ERR
                wrong = rng.random(n) < e
                if lin == "CD8":
                    alt = np.where(rng.random(n) < 0.5, "CD4T", "NK")
                else:
                    alt = np.full(n, "CD8T")
                lab = np.where(wrong, alt, right)
                for i in range(n):
                    cells.append(
                        {
                            "patient": f"P{p:02d}",
                            "clone_b": f"C{p}_{c}",
                            "tissue": t,
                            "all_celltype": lab[i],
                            "g_CD8A": float(cd8a[i]),
                            "g_CD8B": float(cd8b[i]),
                            "g_CD4": float(cd4[i]),
                            "u_CD8A": int(cd8a[i]),
                        }
                    )
                truth.append(
                    {
                        "patient": f"P{p:02d}",
                        "clone": f"C{p}_{c}",
                        "tissue": t,
                        "true_identity": right,
                    }
                )
    return pd.DataFrame(cells), pd.DataFrame(truth)


def run(reps: int) -> pd.DataFrame:
    """Errore stimato e vero medio per condizione x schema x compartimento."""
    loco = load_loco()
    rows = []
    for ci, (cond, art) in enumerate(ARTIFACTS.items()):
        for r in range(reps):
            rng = np.random.default_rng(7000 + r + 10_000 * ci)
            cells, truth = simulate(rng, art)
            per = loco.per_comp_stats(cells, None)
            for scheme in ("naive", "loco", "blood"):
                j = loco.arbitrate(per, scheme)
                if j.empty:
                    continue
                j = j.merge(truth, on=["patient", "clone", "tissue"], how="left")
                # errore vero: etichetta di maggioranza diversa dal lignaggio simulato
                j["true_error"] = j.lab != j.true_identity
                for t, g in j.groupby("tissue"):
                    rows.append(
                        {
                            "condition": cond,
                            "rep": r,
                            "scheme": scheme,
                            "tissue": t,
                            "est_error": g.error.mean(),
                            "true_error": g.true_error.mean(),
                            "n": len(g),
                        }
                    )
    d = pd.DataFrame(rows)
    out = (
        d.groupby(["condition", "scheme", "tissue"])
        .agg(
            est_error=("est_error", "mean"),
            est_error_max=("est_error", "max"),
            true_error=("true_error", "mean"),
            reps=("rep", "nunique"),
            mean_n=("n", "mean"),
        )
        .reset_index()
    )
    out["bias"] = out.est_error - out.true_error
    return out


def main() -> None:
    """Esegue la simulazione e scrive la tabella."""
    ap = argparse.ArgumentParser(
        description="Simulazione (b): schemi di riferimento naive, LOCO e blood."
    )
    ap.add_argument("--reps", type=int, default=30)
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)
    out = run(args.reps)
    out.to_csv(TAB / "sim_b_reference_schemes.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

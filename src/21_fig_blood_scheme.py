#!/usr/bin/env python3
"""
21_fig_blood_scheme.py: Visualizzazione del tasso di errore di annotazione secondo lo schema primario (riferimento PBMC).

Descrizione Tecnica
-------------------
Questo modulo genera la rappresentazione grafica e la tabella riassuntiva del tasso
di errore dell'annotazione trascrizionale nei compartimenti tissutali (`Tumor` e
`Adjacent_normal`) utilizzando come identita' clonale di riferimento esclusivamente
il sangue periferico (`PBMC`).

Per garantire la perfetta invarianza numerica rispetto alle stime inferenziali
prodotte da `10_loco.py`, lo script importa dinamicamente il modulo `10_loco.py`
e invoca la funzione `report()` (inclusiva di Cluster Bootstrap per paziente)
direttamente sulla tabella pre-calcolata `results/tables/arbitration_blood.csv`.

Esecuzione:
    python src/21_fig_blood_scheme.py
"""

import importlib.util
from pathlib import Path

import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Caricamento dinamico di 10_loco.py per riutilizzare le funzioni di reportistica e cluster bootstrap
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("loco", HERE / "10_loco.py")
loco = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loco)

# Lettura della tabella di arbitrato con riferimento ematico e calcolo degli intervalli di confidenza
TAB, FIG = Path("results/tables"), Path("results/figures")
j = pd.read_csv(TAB / "arbitration_blood.csv")
r = loco.report(j, "BLOOD (riferimento = solo PBMC)")

if r.empty:
    raise SystemExit("Nessun compartimento valutabile: controlla arbitration_blood.csv")

# Esportazione della tabella riassuntiva e generazione del grafico a barre con barre d'errore asimmetriche
r.to_csv(TAB / "error_rate_blood.csv", index=False)
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(
    r.tissue,
    r.error,
    yerr=[r.error - r.lo, r.hi - r.error],
    capsize=5,
    color="#d95f02",
    alpha=0.85,
)
for i, row in enumerate(r.itertuples()):
    ax.text(
        i, row.hi + 0.008, f"{row.error:.3f}\n({row.n_pat} pz)", ha="center", fontsize=9
    )
ax.set_ylabel("tasso d'errore dell'etichetta")
ax.set_title("Schema primario: riferimento = solo sangue (PBMC)")
ax.set_ylim(0, max(r.hi) * 1.3)
fig.tight_layout()
fig.savefig(FIG / "error_rate_blood.png", dpi=140)
print(f"\n[ok] figura: {FIG / 'error_rate_blood.png'}")
print(f"[ok] tabella: {TAB / 'error_rate_blood.csv'}")

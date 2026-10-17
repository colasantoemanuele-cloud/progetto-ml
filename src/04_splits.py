#!/usr/bin/env python3
"""
04_splits.py - Definizione dei task di classificazione e delle partizioni di validazione.

Scopo
-----
Costruisce i task supervisionati dalla sola tabella `obs` (X non viene letta) e, per
ciascuno, tre partizioni:
  - fold_group_5 : StratifiedGroupKFold a 5 fold con il paziente come gruppo
                   (partizione primaria, salvo i task in LOPO_PRIMARY);
  - fold_lopo    : LeaveOneGroupOut (un paziente per fold), primaria per
                   task1a_paired e rare;
  - fold_random_5: StratifiedKFold a 5 fold sulle cellule, senza gruppi: controllo
                   negativo che misura quanto la correlazione intra-paziente gonfia la stima.
Un controllo verifica che nessun paziente stia in due fold di fold_group_5.

Task: task1a_all e task1a_paired (Tumor vs Adjacent_normal; solo 5 pazienti hanno
entrambi i tessuti, PA04, PA06, PA08, PA09, PA10, e sono identificati qui),
task2_v1 (9 tipi cellulari del tumore), task2_v2 (CAF+Stellate riuniti), task2_hier1
(Immune/Epithelial/Stromal), task2_hier3 (CD4T/CD8T/NK), rare.
Platelet (anucleate, quasi solo nel sangue) ed Endocrine sono escluse a priori: le
Endocrine sono presenti in 12 pazienti, ma solo 2 ne hanno almeno 10 cellule.

Input
-----
  --h5ad  data/interim/sc_raw.h5ad (solo obs, modalita' backed)

Output
------
  data/interim/splits_{tag}.csv          una riga per (task, cellula) con i tre fold
  results/tables/task_summary_{tag}.csv  cellule, classi, pazienti e valutazione primaria per task

Argomenti
---------
  --h5ad (obbligatorio), --tag (default sc)

Uso
---
    python src/04_splits.py --h5ad data/interim/sc_raw.h5ad --tag sc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.model_selection import (
    LeaveOneGroupOut,
    StratifiedGroupKFold,
    StratifiedKFold,
)

TAB = Path("results/tables")
OUT = Path("data/interim")
SEED = 0
N_FOLDS = 5

# Esclusione sistematica di popolazioni confondenti o non valutabili cross-patient
EXCLUDE_ALWAYS = ["Platelet", "Endocrine"]
RARE = ["Acinar", "Schwann", "Plasma"]

# Definizioni tassonomiche per i task di classificazione gerarchica e granulare
CORE_V1 = [
    "CD8T",
    "CD4T",
    "NK",
    "B",
    "Myeloid",
    "Ductal",
    "CAF",
    "Stellate",
    "Endothelial",
]
FIBRO = {"CAF": "Fibroblast", "Stellate": "Fibroblast"}

LINEAGE = {
    "CD8T": "Immune",
    "CD4T": "Immune",
    "NK": "Immune",
    "B": "Immune",
    "Myeloid": "Immune",
    "Plasma": "Immune",
    "Ductal": "Epithelial",
    "Acinar": "Epithelial",
    "CAF": "Stromal",
    "Stellate": "Stromal",
    "Endothelial": "Stromal",
    "Schwann": "Stromal",
}
LYMPHOID = ["CD4T", "CD8T", "NK"]

# Task a supporto sparso per i quali la copertura in 5-fold non e' garantita (valutazione LOPO)
EXPECT_INCOMPLETE = {"rare", "task1a_all"}
LOPO_PRIMARY = {"task1a_paired", "rare"}


def build_tasks(obs: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Costruisce i DataFrame di mapping (cell_id -> paziente, etichetta) per ciascun task."""
    # Filtraggio globale delle popolazioni escluse e isolamento del microambiente tumorale
    o = obs[~obs.all_celltype.isin(EXCLUDE_ALWAYS)].copy()
    tumor = o[o.tissue == "Tumor"]

    # Identificazione dei pazienti con campionamento appaiato (Tumor e Adjacent_normal)
    # colonne fisse: un tessuto assente vale 0 per tutti i pazienti (invece di un KeyError)
    tis = pd.crosstab(o.patients, o.tissue).reindex(
        columns=sorted(set(o.tissue.astype(str)) | {"Tumor", "Adjacent_normal"}),
        fill_value=0,
    )
    paired = sorted(tis.index[(tis["Tumor"] > 0) & (tis["Adjacent_normal"] > 0)])
    print(
        f"[i] pazienti con Tumor: {int((tis['Tumor'] > 0).sum())} | "
        f"con Adjacent_normal: {int((tis['Adjacent_normal'] > 0).sum())} | "
        f"con entrambi: {len(paired)} -> {paired}"
    )

    tasks: dict[str, pd.DataFrame] = {}

    # Task 1a: classificazione del compartimento tissutale (coorte completa vs coorte appaiata)
    t1all = o[o.tissue.isin(["Tumor", "Adjacent_normal"])]
    tasks["task1a_all"] = pd.DataFrame(
        {"patient": t1all.patients, "y": t1all.tissue.astype(str)}
    )

    t1p = t1all[t1all.patients.isin(paired)]
    tasks["task1a_paired"] = pd.DataFrame(
        {"patient": t1p.patients, "y": t1p.tissue.astype(str)}
    )

    # Task 2 (granulare): 9 popolazioni principali (v1) e aggregazione mesenchimale CAF+Stellate (v2)
    core = tumor[tumor.all_celltype.isin(CORE_V1)]
    tasks["task2_v1"] = pd.DataFrame(
        {"patient": core.patients, "y": core.all_celltype.astype(str)}
    )
    tasks["task2_v2"] = pd.DataFrame(
        {"patient": core.patients, "y": core.all_celltype.astype(str).replace(FIBRO)}
    )

    # Task 2 (gerarchico livello 1): discriminazione dei tre compartimenti principali (Immune/Epithelial/Stromal)
    hier = tumor[tumor.all_celltype.isin(LINEAGE)]
    tasks["task2_hier1"] = pd.DataFrame(
        {"patient": hier.patients, "y": hier.all_celltype.astype(str).map(LINEAGE)}
    )

    # Task 2 (gerarchico livello 3): risoluzione fine del compartimento linfoide citotossico/helper
    lym = tumor[tumor.all_celltype.isin(LYMPHOID)]
    tasks["task2_hier3"] = pd.DataFrame(
        {"patient": lym.patients, "y": lym.all_celltype.astype(str)}
    )

    # Task classi rare: popolazioni a bassa prevalenza valutate congiuntamente a due classi di riferimento
    rare = tumor[tumor.all_celltype.isin(RARE + ["Ductal", "Myeloid"])]
    tasks["rare"] = pd.DataFrame(
        {"patient": rare.patients, "y": rare.all_celltype.astype(str)}
    )

    for k, df in tasks.items():
        df["patient"] = df.patient.astype(str)
        df["y"] = df.y.astype(str)
        tasks[k] = df
    return tasks, paired


def assign_folds(df: pd.DataFrame) -> pd.DataFrame:
    """Assegna gli indici di fold per i tre schemi di partizionamento (Group 5-fold, LOPO, Random 5-fold)."""
    y, g = df.y.values, df.patient.values
    df = df.copy()

    # Schema 1: StratifiedGroupKFold vincolato all'identita' del paziente (prevenzione data leakage)
    k = min(N_FOLDS, len(np.unique(g)))
    sgk = StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=SEED)
    df["fold_group_5"] = -1
    for i, (_, te) in enumerate(sgk.split(df, y, groups=g)):
        df.iloc[te, df.columns.get_loc("fold_group_5")] = i

    # Schema 2: Leave-One-Patient-Out (massimizzazione del supporto di training per classi rare)
    df["fold_lopo"] = -1
    for i, (_, te) in enumerate(LeaveOneGroupOut().split(df, y, groups=g)):
        df.iloc[te, df.columns.get_loc("fold_lopo")] = i

    # Schema 3: StratifiedKFold casuale sulle singole cellule (controllo negativo per quantificare il leakage)
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    df["fold_random_5"] = -1
    for i, (_, te) in enumerate(skf.split(df, y)):
        df.iloc[te, df.columns.get_loc("fold_random_5")] = i

    return df


def audit(name: str, df: pd.DataFrame) -> bool:
    """Esegue la verifica formale di disgiunzione dei gruppi e copertura delle classi per ogni fold."""
    primary = "LOPO" if name in LOPO_PRIMARY else "5-fold per paziente"
    print(
        f"\n{'=' * 74}\n{name}: {len(df):,} cellule | {df.y.nunique()} classi | "
        f"{df.patient.nunique()} pazienti | valutazione primaria: {primary}\n{'=' * 74}"
    )

    supp = pd.DataFrame(
        {
            "n_cells": df.y.value_counts(),
            "n_patients": df.groupby("y").patient.nunique(),
        }
    ).sort_values("n_cells", ascending=False)
    print(supp.to_string())

    # Asserzione formale di assenza di data leakage: ogni paziente deve risiedere in un unico fold di test
    pat = pd.crosstab(df.patient, df.fold_group_5)
    leaks = (pat > 0).sum(axis=1)
    assert (leaks == 1).all(), (
        f"paziente in piu' fold: {leaks[leaks > 1].index.tolist()}"
    )
    print("\n[ok] ogni paziente sta in esattamente un fold (nessun leakage di gruppo)")

    # Verifica della presenza di tutte le classi nei fold di test della partizione Stratificata per Gruppo
    ct = pd.crosstab(df.y, df.fold_group_5)
    print(f"\nclasse x fold_group_5:\n{ct.to_string()}")

    empty = ct.columns[(ct == 0).any(axis=0)].tolist()
    ok = True
    if empty:
        bad = {int(c): ct.index[ct[c] == 0].tolist() for c in empty}
        tag = "ATTESO" if name in EXPECT_INCOMPLETE else "[!!] NON ATTESO"
        print(f"\n{tag}: fold con classi assenti: {bad}")
        print("     -> 5-fold non valutabile per queste classi. Usa fold_lopo.")
        ok = name in EXPECT_INCOMPLETE
    else:
        print(
            f"[ok] nessuna classe assente da alcun fold (min = {int(ct.values.min())} cellule)"
        )

    # Diagnostica di copertura per lo schema Leave-One-Patient-Out
    lopo = pd.crosstab(df.y, df.fold_lopo)
    zeros = int((lopo == 0).sum().sum())
    if zeros:
        print(
            f"[i]  LOPO: {zeros}/{lopo.size} celle (classe, paziente) vuote -> "
            "quei pazienti non contribuiscono a quella classe. Metrica per-classe, non macro."
        )
    else:
        print(f"[ok] LOPO: ogni paziente ha tutte le {df.y.nunique()} classi")

    # Misurazione della dispersione dei pazienti nello split casuale (controllo negativo: quanti pazienti per fold)
    rnd = pd.crosstab(df.patient, df.fold_random_5)
    n_pf = (rnd > 0).sum(axis=0).mean()
    print(
        f"[i]  split casuale: {n_pf:.1f}/{df.patient.nunique()} pazienti per fold "
        "-> ogni paziente e' sia in train sia in test (controllo negativo)."
    )
    return ok


def main() -> None:
    """Punto di ingresso per la generazione e validazione delle tabelle di partizionamento."""
    ap = argparse.ArgumentParser(
        description="Task di classificazione e partizioni di validazione (per paziente e casuali)."
    )
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--tag", default="sc")
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    # Apertura in modalita' backed='r' per leggere esclusivamente la tabella obs
    adata = ad.read_h5ad(args.h5ad, backed="r")
    if adata.n_obs == 0:
        raise SystemExit(f"[ERRORE] {args.h5ad} non contiene cellule")
    obs = adata.obs.copy()
    for c in ("tissue", "patients", "all_celltype"):
        obs[c] = obs[c].astype(str)
    print(f"[ok] {args.h5ad.name}: {adata.n_obs:,} cellule")

    tasks, paired = build_tasks(obs)

    # Assegnazione dei fold e audit di integrita' per ciascun task
    all_ok, frames, rows = True, [], []
    for name, df in tasks.items():
        # servono almeno 2 pazienti e 2 classi per uno split per paziente
        if df.patient.nunique() < 2 or df.y.nunique() < 2:
            print(
                f"[!!] {name}: {df.patient.nunique()} pazienti, {df.y.nunique()} classi: "
                "task saltato (servono almeno 2 pazienti e 2 classi)"
            )
            continue
        df = assign_folds(df)
        all_ok &= audit(name, df)
        df.insert(0, "task", name)
        frames.append(df.reset_index().rename(columns={"index": "cell_id"}))
        rows.append(
            {
                "task": name,
                "n_cells": len(df),
                "n_classes": df.y.nunique(),
                "n_patients": df.patient.nunique(),
                "primary_eval": "LOPO" if name in LOPO_PRIMARY else "group_5fold",
            }
        )

    if not frames:
        raise SystemExit("[ERRORE] nessun task valutabile: nessuno split scritto.")
    out = OUT / f"splits_{args.tag}.csv"
    pd.concat(frames).to_csv(out, index=False)
    print(f"\n[ok] scritto {out}")

    summary = pd.DataFrame(rows)
    summary.to_csv(TAB / f"task_summary_{args.tag}.csv", index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\n[i] pazienti appaiati (Tumor + Adjacent_normal): {paired}")

    if not all_ok:
        print(
            "\n[!!] Un task CORE ha fold incompleti. Non addestrare prima di aver letto sopra."
        )
        sys.exit(2)
    print("\n[ok] audit superato.")


if __name__ == "__main__":
    main()

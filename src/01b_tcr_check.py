#!/usr/bin/env python3
"""
01b_tcr_check.py - Allineamento dei barcode TCR e GEX su una libreria pilota (PDAC01-Tumor).

Scopo
-----
Usa il file VDJ pilota scaricato da 01_probe_labels.py (PA01 Tumor) e i metadati GEX.
L'intersezione dei barcode 10x a 16 nt si calcola dentro la stessa libreria
(PDAC01-Tumor), non sull'insieme di tutti i barcode: la whitelist 10x 5' v1/v2 ha
737.280 sequenze, quindi un confronto globale contro ~175.000 barcode GEX darebbe
circa il 24% di coincidenze casuali. Lo script riporta l'intersezione osservata, quella
attesa per caso e il rapporto; le etichette GEX delle cellule con TCR produttivo; la
frazione di cellule etichettate T con TCR rilevato; le duttali per paziente e tessuto.
Gli identificativi GEX hanno due forme, entrambe riconosciute:
"PDAC01-Tumor-GEX_out_<barcode>-1" e "PDAC05_PBMC_GEX_out_<barcode>-1".

Input
-----
  --probe  cartella con sc_metadata.csv.gz e tcr_PA01_tumor.csv.gz (da 01_probe_labels.py)

Output
------
  solo testo a video (nessun file scritto)

Argomenti
---------
  --probe (default data/raw/probe)

Uso
---
    python src/01b_tcr_check.py --probe data/raw/probe
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Cardinalita teorica della whitelist di barcode cellulari per chimica 10x Genomics 5' v1/v2.
WHITELIST_5P = 737_280
T_LABELS = {
    "CD4T",
    "CD8T",
}  # Popolazioni linfocitarie T canoniche dotate di recettore alfa/beta.
NKT_LABELS = {
    "CD4T",
    "CD8T",
    "NK",
}  # Insieme linfoide citotossico esteso (include cellule NK/NKT).


def main() -> None:
    """Esegue la verifica di sovrapposizione dei barcode TCR/GEX entro libreria e il controllo delle etichette."""
    ap = argparse.ArgumentParser(
        description="Allineamento dei barcode TCR e GEX sulla libreria pilota PDAC01-Tumor."
    )
    ap.add_argument("--probe", type=Path, default=Path("data/raw/probe"))
    args = ap.parse_args()

    sc = pd.read_csv(args.probe / "sc_metadata.csv.gz", index_col=0, low_memory=False)
    sc.index.name = "cell_id"
    tcr = pd.read_csv(args.probe / "tcr_PA01_tumor.csv.gz", low_memory=False)

    # 1. Estrazione del prefisso di libreria e del barcode cellulare a 16 nt dall'indice composito GEX
    # Formati del dataset: PDAC01-Tumor-GEX_out_AAACCTGAGATAGTCA-1 e
    # PDAC05_PBMC_GEX_out_AAACCTGAGATAGTCA-1 (separatore "_" invece di "-").
    parts = sc.index.to_series().str.extract(
        r"^(?P<library>.+?)[-_]GEX_out_(?P<bc>[ACGT]{16})"
    )
    sc = sc.join(parts)
    n_bad = sc.library.isna().sum()
    print(f"celle totali: {len(sc):,} | indice non parsato: {n_bad:,}")
    print(f"librerie GEX distinte: {sc.library.nunique()}")
    print(sc.library.value_counts().head(6).to_string())

    # Isolamento delle sole cellule appartenenti alla libreria PDAC01-Tumor per evitare collisioni inter-libreria.
    LIB = "PDAC01-Tumor"
    gex = sc[sc.library == LIB]
    if gex.empty:
        raise SystemExit(
            f"[ERRORE] libreria '{LIB}' non trovata. Guarda i nomi qui sopra."
        )

    # 2. Filtraggio di qualita dei contig V(D)J: ritenzione dei soli riarrangiamenti produttivi TRA/TRB in cellule reali.
    print(f"\ncontig TCR grezzi: {len(tcr):,}")
    for col, val in (
        ("is_cell", True),
        ("high_confidence", True),
        ("productive", True),
    ):
        if col in tcr.columns:
            keep = tcr[col].astype(str).str.lower().isin(["true", "1"])
            tcr = tcr[keep]
            print(f"  dopo {col}=True: {len(tcr):,}")
    if "chain" in tcr.columns:
        tcr = tcr[tcr.chain.isin(["TRA", "TRB"])]
        print(f"  dopo chain in TRA/TRB: {len(tcr):,}")

    tcr_bc = set(
        tcr.barcode.astype(str).str.extract(r"([ACGT]{16})", expand=False).dropna()
    )
    gex_bc = set(gex.bc)

    print(f"\nbarcode TCR unici (PA01-Tumor, filtrati): {len(tcr_bc):,}")
    print(f"barcode GEX unici (libreria {LIB}):        {len(gex_bc):,}")

    # 3. Calcolo dell'intersezione osservata rispetto al valore atteso per collisione ipergeometrica sulla whitelist 5'.
    inter = tcr_bc & gex_bc
    exp_chance = len(tcr_bc) * len(gex_bc) / WHITELIST_5P

    print("\n" + "=" * 70)
    print(f"intersezione osservata : {len(inter):,}")
    print(f"attesa per puro caso   : {exp_chance:,.0f}   (whitelist {WHITELIST_5P:,})")
    print(f"arricchimento          : {len(inter) / max(exp_chance, 1e-9):.1f}x")
    print(f"frazione TCR ritrovata in GEX : {len(inter) / max(len(tcr_bc), 1):.1%}")
    print(f"frazione GEX con un TCR       : {len(inter) / max(len(gex_bc), 1):.1%}")
    print("=" * 70)

    # 4. Distribuzione delle annotazioni trascrittomiche sulle cellule dotate di recettore TCR produttivo.
    sub = gex[gex.bc.isin(inter)]
    print(
        f"\n>>> Validazione su PA01-Tumor: etichette delle {len(sub):,} cellule con TCR produttivo\n"
    )
    vc = sub.all_celltype.value_counts()
    pct = (vc / len(sub) * 100).round(1)
    print(pd.DataFrame({"n": vc, "%": pct}).to_string())

    prec_t = sub.all_celltype.isin(T_LABELS).mean()
    prec_tnk = sub.all_celltype.isin(NKT_LABELS).mean()
    print(f"\nprecision (etichetta in CD4T/CD8T)     : {prec_t:.1%}")
    print(f"precision (etichetta in CD4T/CD8T/NK)   : {prec_tnk:.1%}")
    print(
        "\nUn riarrangiamento TCR produttivo identifica una cellula T. Le cellule non-T"
    )
    print(
        "in questa tabella riflettono doppietti etorotipici, RNA ambientale o discordanza di annotazione."
    )
    print(
        "Questo rapporto fornisce una stima ortogonale del tasso di discordanza dell'annotazione."
    )

    # Stima della sensibilita di cattura V(D)J sulle cellule annotate trascrizionalmente come linfociti T.
    t_cells = gex[gex.all_celltype.isin(T_LABELS)]
    recall_like = t_cells.bc.isin(inter).mean()
    print(f"\ncellule etichettate T con TCR rilevato  : {recall_like:.1%}")
    print(
        "(misura l'efficienza sperimentale di cattura V(D)J 5', tipicamente compresa tra 60% e 85%)"
    )

    # 5. Conteggio delle cellule epiteliali duttali per paziente e compartimento tissutale.
    print("\n" + "=" * 70)
    print("Analisi cellule DUTTALI per tessuto e paziente")
    print("=" * 70)
    duct = sc[(sc.all_celltype == "Ductal") & (sc.tissue != "PBMC")]
    print(pd.crosstab(duct.patients, duct.tissue, margins=True).to_string())
    print(
        "\nI pazienti dotati di cellule duttali in entrambi i compartimenti solidi (Tumor e Adjacent_normal)"
    )
    print(
        "consentono il contrasto intra-paziente controllato per l'effetto individuale."
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
19_ni_score.py: Test di associazione tra punteggio clinico di invasione neurale (NI score) ed errore di annotazione.

Razionale Clinico e Ipotesi H4
------------------------------
La Tabella Supplementare S1A (`mmc2.xlsx` dello studio di riferimento Chen et al.,
Cancer Cell 2025) riporta i metadati clinici e istopatologici dei pazienti, incluso
il punteggio quantitativo di invasione neurale (`NI score`). I 14 pazienti con
`scRNA-seq = Yes` nella tabella clinica corrispondono biunivocamente ai 14 pazienti
del dataset single-cell.

Si formula l'ipotesi a priori (H4) che un `NI score` piu' elevato, indicativo di un
microambiente tumorale maggiormente rimodellato dall'invasione perineurale e
associato a disfunzione immunitaria, correli positivamente con il tasso di errore
di annotazione dei linfociti T nel compartimento tumorale. Poiche' la dicotomizzazione
clinica (`NI group`: 11 High vs 3 Low) presenta un marcato sbilanciamento di classe,
l'inferenza viene condotta direttamente sulla variabile continua `NI score`.

Controllo Congiunto della Molteplicita' (Famiglia di 4 Ipotesi)
---------------------------------------------------------------
Poiche' H4 viene valutata contro il medesimo outcome a livello di paziente (tasso
d'errore nel tessuto tumorale da `arbitration_blood.csv`) gia' impiegato per le tre
ipotesi trascrizionali H1-H3 in `16_pseudobulk_patient.py`, testare H4 isolatamente
senza aggiornare la correzione per test multipli determinerebbe un'inflazione del
Family-Wise Error Rate attraverso sessioni successive. Questo modulo carica i p-value
di permutazione di H1-H3 e ricalcola la procedura step-down di Holm-Bonferroni
sull'intera famiglia congiunta di 4 ipotesi (H1, H2, H3, H4).

Esecuzione:
    python src/19_ni_score.py \
        --table-s1 mmc2.xlsx \
        --arbitration results/tables/arbitration_blood.csv \
        --previous-tests results/tables/pseudobulk_hypothesis_tests.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd

# Percorsi di output e costanti del protocollo inferenziale
FIG = Path("results/figures")
TAB = Path("results/tables")
SEED = 0
N_PERM = 20000
MIN_CLONES_OUTCOME = 5

# Schema atteso delle colonne nella Tabella Supplementare S1A (mmc2.xlsx)
EXPECTED_COLS = [
    "PatientID",
    "Age",
    "Gender",
    "Clinical Stage",
    "pTNM",
    "scRNA-seq",
    "NI score",
    "NI group",
]


# --------------------------------------------------------------------------- #
def clean_patient_id(raw: str) -> str:
    """Normalizza l'identificativo grezzo del paziente nel formato canonico PAxx
    a due cifre per garantire l'allineamento con le chiavi di scRNA-seq."""
    m = re.search(r"(\d+)", str(raw))
    if not m:
        raise ValueError(f"ID paziente non riconosciuto: {raw!r}")
    return f"PA{int(m.group(1)):02d}"


def load_table_s1a(path: Path) -> pd.DataFrame:
    """Carica e valida il foglio 'Table S1A' dal file Excel supplementare,
    verificando l'integrita' delle intestazioni e normalizzando gli ID paziente."""
    wb = openpyxl.load_workbook(path, data_only=True)
    if "Table S1A" not in wb.sheetnames:
        raise SystemExit(
            f"[ERRORE] foglio 'Table S1A' assente. Fogli trovati: {wb.sheetnames}"
        )
    ws = wb["Table S1A"]
    rows = list(ws.iter_rows(values_only=True))
    # Estrazione delle intestazioni dalla seconda riga e verifica delle colonne obbligatorie
    header = [str(h).strip() if h else h for h in rows[1]]
    missing = [c for c in EXPECTED_COLS if c not in header]
    if missing:
        raise SystemExit(
            f"[ERRORE] colonne attese assenti: {missing}\n"
            f"colonne trovate: {header}\n"
            "Il formato del supplementare potrebbe essere cambiato: controlla a mano."
        )
    data = [dict(zip(header, r)) for r in rows[2:] if r[0] is not None]
    df = pd.DataFrame(data)
    df["patient"] = df.PatientID.apply(clean_patient_id)
    return df


def permutation_corr_test(
    x: np.ndarray, y: np.ndarray, n_perm: int = N_PERM, seed: int = SEED
):
    """Esegue il test di permutazione esatto a due code (20.000 permutazioni
    vettorizzate dei ranghi) per il coefficiente di correlazione di Spearman."""
    rng = np.random.default_rng(seed)
    n = len(x)
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rho_obs = float(np.corrcoef(rx, ry)[0, 1])
    perm_idx = np.argsort(rng.random((n_perm, n)), axis=1)
    ry_perm = ry[perm_idx]
    rx_c = rx - rx.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    rhos = num / np.maximum(den, 1e-12)
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p)


def bootstrap_ci_corr(
    x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = SEED
):
    """Calcola l'intervallo di confidenza al 95% per la correlazione di Spearman
    mediante ricampionamento bootstrap non parametrico sui pazienti."""
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
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def holm_bonferroni(pvals: dict) -> dict:
    """Applica la correzione sequenziale di Holm-Bonferroni all'insieme completo
    dei p-value della famiglia di ipotesi."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, running_max = {}, 0.0
    for i, (name, p) in enumerate(items):
        val = min((m - i) * p, 1.0)
        running_max = max(running_max, val)
        adj[name] = running_max
    return adj


def build_outcome(arb_path: Path) -> pd.DataFrame:
    """Calcola il tasso medio di errore di annotazione per paziente nel tessuto
    tumorale utilizzando esattamente la medesima soglia di 16_pseudobulk_patient.py."""
    if not arb_path.exists():
        raise SystemExit(f"[ERRORE] {arb_path} assente. Esegui prima 10_loco.py.")
    j = pd.read_csv(arb_path)
    j = j[j.tissue == "Tumor"]
    out = j.groupby("patient").agg(
        error_rate=("error", "mean"), n_clones=("error", "size")
    )
    out = out[out.n_clones >= MIN_CLONES_OUTCOME]
    return out.reset_index()


# --------------------------------------------------------------------------- #
def main() -> None:
    """Orchestra l'allineamento dei metadati clinici NI score, il test di
    permutazione per H4 e il ricalcolo congiunto di Holm-Bonferroni su H1-H4."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1", type=Path, required=True)
    ap.add_argument(
        "--arbitration", type=Path, default=Path("results/tables/arbitration_blood.csv")
    )
    ap.add_argument(
        "--previous-tests",
        type=Path,
        default=Path("results/tables/pseudobulk_hypothesis_tests.csv"),
    )
    args = ap.parse_args()

    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    # Caricamento dei metadati clinici e selezione della coorte sottoposta a scRNA-seq
    clin = load_table_s1a(args.table_s1)
    scrna = clin[clin["scRNA-seq"] == "Yes"].copy()
    print(f"[ok] pazienti con scRNA-seq=Yes in Table S1A: {len(scrna)}")
    print(
        scrna[["patient", "Age", "Gender", "Clinical Stage", "NI score", "NI group"]]
        .sort_values("patient")
        .to_string(index=False)
    )

    print(
        f"\n[i] NI group: {scrna['NI group'].value_counts().to_dict()} "
        "-> troppo sbilanciato per un test binario, uso solo NI score (continuo)."
    )

    # Join relazionale tra tasso di errore nel tumore e variabili cliniche per paziente
    outcome = build_outcome(args.arbitration)
    df = outcome.merge(
        scrna[["patient", "NI score", "Clinical Stage", "Age"]],
        on="patient",
        how="inner",
    )
    print(f"\n[ok] join outcome x NI score: {len(df)} pazienti -> {sorted(df.patient)}")
    df.to_csv(TAB / "ni_score_joined.csv", index=False)

    print("\n" + "#" * 78)
    print("# IPOTESI H4 (dichiarata prima di guardare il risultato):")
    print("# NI score piu' alto -> microambiente piu' immunosoppresso -> piu'")
    print("# T esauste -> piu' errore di annotazione nel tumore. Attesa: rho > 0.")
    print("#" * 78)

    # Test di correlazione per H4 mediante permutazione esatta e intervallo bootstrap
    x, y = df["NI score"].values, df.error_rate.values
    rho, p = permutation_corr_test(x, y)
    lo, hi = bootstrap_ci_corr(x, y)
    print(
        f"\nn={len(df)} pazienti | rho={rho:+.3f}  p_perm={p:.4f}  IC95% boot=[{lo:+.3f}, {hi:+.3f}]"
    )
    print(f"{'COERENTE con H4 (rho>0)' if rho > 0 else 'direzione OPPOSTA a H4'}")

    # Correzione congiunta di Holm-Bonferroni integrando i p-value di H1, H2 e H3
    print("\n" + "=" * 78)
    print(
        "CORREZIONE CONGIUNTA: H4 insieme a H1 (esaurimento), H2 (mieloidi), H3 (CAF)"
    )
    print("=" * 78)
    all_p = {"H4_ni_score": p}
    if args.previous_tests.exists():
        prev = pd.read_csv(args.previous_tests, index_col=0)
        for h in prev.index:
            all_p[h] = float(prev.loc[h, "p_perm"])
        print(f"[ok] caricati H1-H3 da {args.previous_tests.name}")
    else:
        print(
            f"[!] {args.previous_tests} assente: correggo solo H4 da sola (non ideale, esegui prima 16)."
        )

    adj = holm_bonferroni(all_p)
    print(f"\np-value grezzi: {all_p}")
    print("\np-value corretti (Holm-Bonferroni, famiglia completa):")
    for k, v in adj.items():
        sig = "SIGNIFICATIVO" if v < 0.05 else "non significativo"
        marker = "  [H4: questo script]" if k == "H4_ni_score" else ""
        print(f"  {k:20s} p_adj = {v:.4f}   {sig}{marker}")

    pd.Series(adj, name="p_holm_family_of_4").to_csv(TAB / "ni_score_holm_joint.csv")

    print("\n" + "=" * 78)
    print(f"n={len(df)} pazienti: come per H1-H3, questo e' un test a bassa potenza.")
    print("Un p<0.05 qui e' un indizio da prendere sul serio (variabile clinica")
    print("indipendente, non costruita da noi) ma NON prova causalita', e andrebbe")
    print("verificato su una coorte esterna prima di essere riportato come risultato.")
    print("=" * 78)

    # Generazione del grafico NI score vs tasso di errore colorato per stadio clinico TNM
    fig, ax = plt.subplots(figsize=(7, 5.5))
    stage_order = ["IA", "IB", "IIA", "IIB", "III"]
    cmap = plt.cm.viridis
    colors = [
        cmap(stage_order.index(s) / (len(stage_order) - 1))
        if s in stage_order
        else "gray"
        for s in df["Clinical Stage"]
    ]
    sc_ = ax.scatter(  # noqa: F841
        df["NI score"], df.error_rate, c=colors, s=160, edgecolor="k", zorder=5
    )
    for _, r in df.iterrows():
        ax.annotate(
            r.patient,
            (r["NI score"], r.error_rate),
            fontsize=8,
            xytext=(5, 5),
            textcoords="offset points",
        )
    ax.set_xlabel("NI score (invasione neurale, dallo studio originale)")
    ax.set_ylabel("tasso d'errore di annotazione (Tumor)")
    ax.set_title(f"NI score vs errore: rho={rho:+.3f}, p={p:.4f} (n={len(df)})")
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=cmap(i / (len(stage_order) - 1)),
            markersize=9,
            markeredgecolor="k",
            label=s,
        )
        for i, s in enumerate(stage_order)
    ]
    ax.legend(handles=handles, title="stadio clinico", fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(FIG / "ni_score_vs_error.png", dpi=140)
    plt.close(fig)
    print(f"\n[ok] figura: results/figures/ni_score_vs_error.png")  # noqa: F541


if __name__ == "__main__":
    main()

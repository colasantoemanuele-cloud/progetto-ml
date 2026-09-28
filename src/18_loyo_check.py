#!/usr/bin/env python3
"""
18_loyo_check.py: Analisi di sensibilita' Leave-One-Patient-Out (LOYO) sull'espansione clonale.

Razionale Metodologico
----------------------
Nel modulo 17 (`17_clonal_expansion.py`), l'analisi stratificata per blocchi ha
evidenziato una correlazione positiva tra la dimensione clonale T e sia il tasso di
errore di annotazione nel tumore (rho = +0.177, p = 0.0009) sia l'eccesso continuo di
discordanza (rho = +0.281, p < 0.0001) su 287 cloni distribuiti in 10 pazienti.
Tuttavia, quando il numero di cluster e' ridotto (n = 10 pazienti) e la distribuzione
del numero di cloni per paziente e' asimmetrica, anche un test di permutazione a
blocchi puo' risultare influenzato in modo sproporzionato da un singolo paziente con
elevata numerosita' clonale.

Protocollo di Validazione Leave-One-Patient-Out
-----------------------------------------------
Per accertare che l'associazione osservata costituisca una proprieta' riproducibile
inter-paziente e non un artefatto guidato da un singolo individuo influente, questo
script esegue una procedura di sensibilita' Leave-One-Patient-Out completa:
  - Per ciascuno dei 10 pazienti della coorte, vengono rimossi tutti i cloni
    appartenenti a quel paziente;
  - Sul sottoinsieme dei 9 pazienti residui viene ricalcolata la correlazione di
    Spearman pooled e il relativo p-value tramite Block Permutation Test (10.000
    permutazioni intra-paziente);
  - Viene verificata la stabilita' del segno, della magnitudine dell'effetto (`delta_rho`)
    e della significativita' statistica (`p < 0.05`) lungo tutte le 10 iterazioni.

Esecuzione:
    python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

# Configurazione del backend grafico non interattivo per esecuzione batch
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    """Calcola la correlazione di Spearman pooled e il p-value a due code tramite
    permutazione ristretta ai blocchi paziente (identico a 17_clonal_expansion.py)."""
    rng = np.random.default_rng(seed)
    # Ordinamento per gruppo paziente per identificare i limiti di ciascun blocco
    order = np.argsort(group, kind="stable")
    x_s, g_s = x[order], group[order]
    ry_s = pd.Series(y[order]).rank().values
    rx_s = pd.Series(x_s).rank().values
    rho_obs = float(np.corrcoef(rx_s, ry_s)[0, 1])

    # Generazione delle permutazioni indipendenti all'interno di ciascun paziente
    n = len(x_s)
    perm_idx = np.tile(np.arange(n), (n_perm, 1))
    uniq, starts = np.unique(g_s, return_index=True)
    bounds = list(starts) + [n]
    for i in range(len(uniq)):
        s, e = bounds[i], bounds[i + 1]
        size = e - s
        if size <= 1:
            continue
        keys = rng.random((n_perm, size))
        local = np.argsort(keys, axis=1)
        perm_idx[:, s:e] = s + local

    # Calcolo vettorizzato della statistica di correlazione sulle repliche permutate
    ry_perm = ry_s[perm_idx]
    rx_c = rx_s - rx_s.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    rhos = num / np.maximum(den, 1e-12)
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p)


def loyo(df: pd.DataFrame, outcome_col: str, label: str) -> pd.DataFrame:
    """Esegue il ciclo Leave-One-Patient-Out rimuovendo iterativamente un paziente
    alla volta e ricalcolando correlazione e significativita' sui pazienti residui."""
    print(f"\n{'=' * 78}\n{label}\n{'=' * 78}")
    pats = sorted(df.patient.unique())

    # Calcolo del valore di riferimento sull'intera coorte di 10 pazienti
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
            }
        )
        flag = "" if p < 0.05 else "  [non piu' significativo senza questo paziente]"
        print(
            f"  senza {held_out} ({n_clones_removed:3d} cloni): rho={rho:+.3f}  p={p:.4f}{flag}"
        )

    # Sintesi della stabilita' inferenziale attraverso le 10 iterazioni leave-one-out
    out = pd.DataFrame(rows)
    n_robust = out.still_significant.sum()
    print(
        f"\nrisultato resta significativo (p<0.05) in {n_robust}/{len(pats)} versioni leave-one-out."
    )
    if n_robust == len(pats):
        print(
            "VERDETTO: nessun singolo paziente e' indispensabile. Il segnale e' cross-patient."
        )
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
    ap = argparse.ArgumentParser()
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
    print(
        f"[ok] {len(df):,} cloni, {df.patient.nunique()} pazienti, da {args.input.name}"
    )

    # Iterazione sui due outcome: errore binario (riferimento PBMC) ed eccesso di discordanza
    for outcome_col, label, fname_suffix in [
        ("error", "LOYO: errore di annotazione (0/1)", "error"),
        ("excess", "LOYO: eccesso di discordanza (continuo)", "excess"),
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
            f"Leave-one-patient-out: {label}",
        )

    print("\n" + "=" * 78)
    print("Se entrambi gli outcome (errore ed eccesso) restano significativi in")
    print("tutte le versioni leave-one-out, il risultato di 17 e' da considerare")
    print("solido: non e' l'artefatto di un paziente con molti cloni o un pattern")
    print("estremo. Se non lo sono, questo costituisce il risultato da riportare.")
    print("=" * 78)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
sim_a_calibration.py - Simulazione (a): Nadeau-Bengio e cluster bootstrap sotto H0 e con effetto noto.

Scopo
-----
Verifica su dati simulati i due stimatori usati per l'inferenza nella pipeline.

1. t corretto di Nadeau-Bengio (`scpipe.nadeau_bengio_ttest`). Disegno ripreso da
   `genomic-audit-mvp/tests/test_stats_calibration.py` (stesso autore):
   - H0: due DummyClassifier(strategy="stratified") con semi diversi, confrontati su
     una StratifiedKFold a k fold (k = 3 e 5), n = 90, y ~ Bernoulli(0.5), X costante.
     I due modelli hanno la stessa prestazione attesa: la frazione di p < 0.05 e' il
     tasso di falsi positivi. 1.000 repliche per k.
   - effetto noto: LogisticRegression contro Dummy su n = 120, 8 feature N(0, 1),
     pesi N(0, 1) * 8, y logistica; k = 5; 300 repliche. Frazione di p < 0.05 = potenza.
   Per confronto si riporta anche il t appaiato non corretto.
2. Cluster bootstrap della media pooled (`scpipe.cluster_bootstrap`, intervallo
   "expanded percentile") e, per confronto, l'intervallo percentile semplice usato
   nella versione precedente. Pazienti k in {4, 5, 10, 11} con numerosita' diseguale
   (5-59 osservazioni per paziente: per k <= 10 i primi k valori di `SIZES`); effetto
   paziente casuale.
   Numerosita': per k = 10 i numeri di cloni dei 10 pazienti di
   `clonal_expansion_joined` (5-59 per paziente); per k = 11 si aggiunge un undicesimo
   valore (40) scelto a mano. Le numerosita' reali delle coppie dell'eccesso (11
   pazienti, 13-282 coppie per paziente in `pairwise_discordance_clone_b`) sono piu'
   sbilanciate: negli scenari principali la copertura e' verificata con `SIZES`, non con
   le numerosita' reali delle coppie. Per queste c'e' uno scenario a parte
   (`coverage_pairs_real_sizes_{continuous,binary}`, k = 11, `PAIR_SIZES` = coppie per
   paziente di `pairwise_discordance_clone_b`, seme 2011), che non cambia gli altri.
   - esito continuo: media vera 0.06, effetto paziente N(0, 0.05), rumore N(0, 0.15);
   - esito binario: p del paziente ~ Beta(2, 8) (media 0.2);
   - H0 sulla media: media vera 0 (continuo); frazione di IC che escludono 0;
   - effetto noto: media vera 0.06 (continuo); frazione di IC che escludono 0 (potenza).
   400 repliche per scenario, 1.000 repliche bootstrap ciascuna (le analisi ne usano
   2.000; la copertura non e' stata verificata con 2.000 repliche).
   Gli scenari continui non sono simulazioni indipendenti: per ogni k il generatore
   riparte dallo stesso seme (1000 + k), quindi `H0_mean0` usa gli stessi dati di
   `coverage_continuous` traslati (media vera 0 invece di 0.06) ed `effect_mean006` usa
   gli stessi dati, identici. Per costruzione la frazione di
   IC che escludono 0 sotto H0 e' 1 - copertura (es. 0,18 = 1 - 0,82 con 4 pazienti):
   non e' una verifica separata della copertura.

Input
-----
  nessuno (dati simulati, semi fissi)

Output
------
  results/tables/sim_a_nadeau_bengio.csv     scenario, k, repliche, frazione p < 0.05 (corretto e non)
  results/tables/sim_a_cluster_bootstrap.csv scenario, k, metodo, copertura / frazione che esclude 0
  (con --quick i file hanno il suffisso `_quick` e le tabelle pubblicate non vengono toccate)

Argomenti
---------
  --reps-nb (default 1000), --reps-cb (default 400), --quick (repliche ridotte a 1/10,
  uscita in `sim_a_*_quick.csv`);
  i default sono anche in `--help`

Uso
---
    python src/sim/sim_a_calibration.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scpipe import (  # noqa: E402
    cluster_bootstrap_replicates,
    nadeau_bengio_ttest,
    percentile_ci,
)

from sklearn.dummy import DummyClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import f1_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

TAB = Path("results/tables")
SIZES = np.array([17, 26, 47, 42, 22, 10, 30, 59, 29, 5, 40])
# coppie clone-compartimenti per paziente dell'eccesso (pairwise_discordance_clone_b, 11 pazienti)
PAIR_SIZES = np.array([62, 24, 184, 75, 282, 118, 146, 57, 97, 56, 13])
ALPHA = 0.05


def naive_paired_t(a, b) -> float:
    """t appaiato senza correzione per la sovrapposizione dei training (confronto)."""
    d = np.asarray(a) - np.asarray(b)
    if np.allclose(d, 0):
        return 1.0
    s = d.std(ddof=1)
    if s == 0:
        return 0.0
    t = d.mean() / (s / np.sqrt(len(d)))
    return float(2 * student_t.sf(abs(t), len(d) - 1))


def nb_replicate(rep: int, seed0: int, n: int, k: int, effect: bool):
    """Una replica: punteggi per fold dei due modelli, p corretto e non corretto."""
    rng = np.random.default_rng(seed0 + rep)
    if effect:
        X = rng.normal(size=(n, 8))
        w = rng.normal(size=8) * 8.0
        y = (rng.random(n) < 1 / (1 + np.exp(-(X @ w)))).astype(int)
    else:
        y = rng.integers(0, 2, n)
        X = np.zeros((n, 1))
    if np.bincount(y, minlength=2).min() < k:
        return np.nan, np.nan
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=rep)
    fa, fb, ntr, nte = [], [], [], []
    for tr, te in skf.split(X, y):
        if effect:
            a = LogisticRegression(max_iter=1000).fit(X[tr], y[tr])
        else:
            a = DummyClassifier(strategy="stratified", random_state=2 * rep).fit(
                X[tr], y[tr]
            )
        b = DummyClassifier(strategy="stratified", random_state=2 * rep + 1).fit(
            X[tr], y[tr]
        )
        fa.append(f1_score(y[te], a.predict(X[te]), average="macro", zero_division=0))
        fb.append(f1_score(y[te], b.predict(X[te]), average="macro", zero_division=0))
        ntr.append(len(tr))
        nte.append(len(te))
    p_nb = nadeau_bengio_ttest(fa, fb, int(np.mean(ntr)), int(np.mean(nte)))
    return p_nb, naive_paired_t(fa, fb)


def sim_nb(reps: int) -> pd.DataFrame:
    """Frazione di p < 0.05 di Nadeau-Bengio e del t non corretto, sotto H0 e con effetto."""
    rows = []
    scen = [("H0", 90, 3, False, 55_003), ("H0", 90, 5, False, 55_005)]
    scen.append(("effect", 120, 5, True, 51_005))
    for name, n, k, effect, seed0 in scen:
        r = reps if not effect else max(reps * 3 // 10, 30)
        ps = np.array([nb_replicate(i, seed0, n, k, effect) for i in range(r)])
        ps = ps[~np.isnan(ps[:, 0])]
        rows.append(
            {
                "scenario": name,
                "n": n,
                "k": k,
                "reps": len(ps),
                "frac_p_lt_005_nadeau_bengio": float(np.mean(ps[:, 0] < ALPHA)),
                "frac_p_lt_005_uncorrected": float(np.mean(ps[:, 1] < ALPHA)),
            }
        )
        print(rows[-1])
    return pd.DataFrame(rows)


def draw_groups(rng, k: int, kind: str, mean: float, sizes=None):
    """Osservazioni simulate di `k` pazienti (dizionario paziente -> array), continue o binarie.

    `sizes` (facoltativo) fissa le numerosita' per paziente; di default i primi k di SIZES.
    """
    if sizes is None:
        sizes = SIZES[:k] if k <= len(SIZES) else rng.integers(5, 60, k)
    if kind == "binary":
        p = rng.beta(2, 8, k)
        return {f"P{i}": (rng.random(sizes[i]) < p[i]).astype(float) for i in range(k)}
    u = rng.normal(0, 0.05, k)
    return {f"P{i}": mean + u[i] + rng.normal(0, 0.15, sizes[i]) for i in range(k)}


def sim_cb(reps: int, n_boot: int = 1000) -> pd.DataFrame:
    """Copertura e frazione di IC che escludono 0 del cluster bootstrap, per scenario e k."""
    rows = []
    scen = [
        ("coverage_continuous", "continuous", 0.06, 0.06),
        ("coverage_binary", "binary", None, 0.2),
        ("H0_mean0", "continuous", 0.0, 0.0),
        ("effect_mean006", "continuous", 0.06, 0.0),
    ]
    for name, kind, mean, target in scen:
        for k in (4, 5, 10, 11):
            rng = np.random.default_rng(1000 + k)
            hit = {"expanded_percentile": 0, "percentile": 0}
            for r in range(reps):
                by = draw_groups(rng, k, kind, mean if mean is not None else 0.0)
                stats, _ = cluster_bootstrap_replicates(by, n_boot, seed=r)
                cis = {
                    "expanded_percentile": percentile_ci(stats, k),
                    "percentile": tuple(np.percentile(stats, [2.5, 97.5])),
                }
                for meth, (lo, hi) in cis.items():
                    if name.startswith("coverage"):
                        hit[meth] += lo <= target <= hi
                    else:  # frazione di IC che escludono lo zero
                        hit[meth] += (lo > 0) or (hi < 0)
            for meth, h in hit.items():
                rows.append(
                    {
                        "scenario": name,
                        "k_patients": k,
                        "method": meth,
                        "reps": reps,
                        "rate": h / reps,
                    }
                )
                print(rows[-1])
    # copertura con le numerosita' reali delle coppie dell'eccesso (k = 11, seme proprio)
    k = len(PAIR_SIZES)
    for kind, target in (("continuous", 0.06), ("binary", 0.2)):
        rng = np.random.default_rng(2011)
        hit = {"expanded_percentile": 0, "percentile": 0}
        for r in range(reps):
            mean = 0.06 if kind == "continuous" else 0.0
            by = draw_groups(rng, k, kind, mean, sizes=PAIR_SIZES)
            stats, _ = cluster_bootstrap_replicates(by, n_boot, seed=r)
            cis = {
                "expanded_percentile": percentile_ci(stats, k),
                "percentile": tuple(np.percentile(stats, [2.5, 97.5])),
            }
            for meth, (lo, hi) in cis.items():
                hit[meth] += lo <= target <= hi
        for meth, h in hit.items():
            rows.append(
                {
                    "scenario": f"coverage_pairs_real_sizes_{kind}",
                    "k_patients": k,
                    "method": meth,
                    "reps": reps,
                    "rate": h / reps,
                }
            )
            print(rows[-1])
    return pd.DataFrame(rows)


def main() -> None:
    """Esegue le due simulazioni e scrive le tabelle."""
    ap = argparse.ArgumentParser(
        description="Simulazione (a): Nadeau-Bengio e cluster bootstrap sotto H0 e con effetto.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--reps-nb",
        type=int,
        default=1000,
        help="repliche per k nel test di Nadeau-Bengio",
    )
    ap.add_argument(
        "--reps-cb",
        type=int,
        default=400,
        help="repliche per scenario del cluster bootstrap",
    )
    ap.add_argument(
        "--quick",
        action="store_true",
        help="repliche ridotte a 1/10; scrive sim_a_*_quick.csv (non le tabelle pubblicate)",
    )
    args = ap.parse_args()
    if args.quick:
        args.reps_nb, args.reps_cb = args.reps_nb // 10, args.reps_cb // 10
    TAB.mkdir(parents=True, exist_ok=True)
    suffix = "_quick" if args.quick else ""
    nb = sim_nb(args.reps_nb)
    nb.to_csv(TAB / f"sim_a_nadeau_bengio{suffix}.csv", index=False)
    cb = sim_cb(args.reps_cb)
    cb.to_csv(TAB / f"sim_a_cluster_bootstrap{suffix}.csv", index=False)
    print("\n" + nb.to_string(index=False))
    print(
        "\n"
        + cb.pivot_table(
            index=["scenario", "k_patients"], columns="method", values="rate"
        ).to_string()
    )


if __name__ == "__main__":
    main()

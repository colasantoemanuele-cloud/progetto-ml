"""Test di calibrazione degli stimatori statistici, su dati simulati con semi fissi.

Bande dichiarate PRIMA dell'esecuzione (non si allargano: se un test esce dalla
banda si corregge lo stimatore):

1. Cluster bootstrap sui pazienti (IC 95% della media pooled), con 10 e con 11
   pazienti di numerosita' diseguale come nei dati (5-59 cloni per paziente),
   esito continuo e binario: copertura empirica dell'IC in [0.90, 0.99].
   400 repliche per scenario. (Fase 0, questo test con i semi 100 + k + 1 se l'esito e'
   binario, log in `revisione/log/fase0_pytest_slow.log`: l'intervallo percentile
   semplice dava 0.877-0.892, fuori banda; lo stimatore e' stato corretto con
   l'intervallo "expanded percentile", la banda e' rimasta invariata. L'88,8-92,3% del
   README (0,8875-0,9225 in `results/tables/sim_a_cluster_bootstrap.csv`) viene invece
   da `src/sim/sim_a_calibration.py`, semi 1000 + k: simulazione diversa, stesso
   disegno.)
2. Block permutation test (permutazione dell'esito entro paziente) sotto H0 con
   confondente di paziente (x e y dipendono entrambi da un effetto di paziente,
   nessuna associazione entro paziente): frazione di p < 0.05 in [0.02, 0.08],
   sia per il p non centrato (quello usato in 17 e 18) sia per il p centrato.
   400 repliche, 2.000 permutazioni ciascuna.
3. Holm-Bonferroni: p corretti monotoni nell'ordine dei p grezzi, >= p grezzi, <= 1
   (proprieta' esatte, nessuna banda).
4. Discordanza a coppie senza reinserimento: valore atteso esatto (enumerazione di
   tutte le composizioni multinomiali) uguale a 1 - sum_k p_k^2; eccesso atteso
   esattamente 0 quando i due compartimenti hanno la stessa distribuzione
   (tolleranza numerica 1e-12, nessuna banda).
5. t di Nadeau-Bengio (disegno ripreso da genomic-audit-mvp): sotto H0 (due
   DummyClassifier stratificati con semi diversi, StratifiedKFold k = 3 e k = 5,
   n = 90) frazione di p < 0.05 in [0.00, 0.08] per ciascun k (il test puo' essere conservativo, non
   deve superare il livello); con effetto noto (logistica contro Dummy, n = 120)
   frazione di p < 0.05 >= 0.50 (k = 5). 400 repliche per k sotto H0, 100 con effetto.
6. Regola della maggioranza, macro-F1 sulle classi presenti, livello espanso:
   proprieta' esatte (nessuna banda).
7. Test di permutazione di Spearman (`scpipe.permutation_corr_test`, unita' = righe,
   usato per H1-H4 in 16 e 19 e per le medie per paziente in 17, con n = 10
   pazienti): sotto H0 (x e y indipendenti, n = 10) frazione di p < 0.05 in
   [0.02, 0.08]. 400 repliche, 2.000 permutazioni.
"""

from __future__ import annotations

from itertools import product
from math import factorial

import numpy as np
import pandas as pd
import pytest

import scpipe
from conftest import load_script

SIZES = np.array([17, 26, 47, 42, 22, 10, 30, 59, 29, 5, 40])


def cluster_ci(by_pat: dict, n_boot: int, seed: int):
    """IC 95% del cluster bootstrap come implementato in 07_robustness.py."""
    m = load_script("07_robustness")
    df = pd.DataFrame(
        {
            "patient": np.concatenate([[p] * len(v) for p, v in by_pat.items()]),
            "excess": np.concatenate(list(by_pat.values())),
        }
    )
    old = m.SEED
    m.SEED = seed
    try:
        return m.cluster_bootstrap(df, "excess", n_boot)
    finally:
        m.SEED = old


def simulate_coverage(k: int, binary: bool, reps: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    sizes = SIZES[:k]
    hits = 0
    for r in range(reps):
        if binary:
            p = rng.beta(2, 8, k)  # media 0.2
            by = {
                f"P{i}": (rng.random(sizes[i]) < p[i]).astype(float) for i in range(k)
            }
            true = 0.2
        else:
            u = rng.normal(0, 0.05, k)
            by = {
                f"P{i}": 0.06 + u[i] + rng.normal(0, 0.15, sizes[i]) for i in range(k)
            }
            true = 0.06
        _, lo, hi = cluster_ci(by, n_boot=1000, seed=r)
        hits += lo <= true <= hi
    return hits / reps


@pytest.mark.slow
@pytest.mark.parametrize("k", [10, 11])
@pytest.mark.parametrize("binary", [False, True])
def test_cluster_bootstrap_coverage(k, binary):
    cov = simulate_coverage(k, binary, reps=400, seed=100 + k + int(binary))
    assert 0.90 <= cov <= 0.99, f"copertura {cov:.3f} fuori banda [0.90, 0.99]"


# --------------------------------------------------------------------------- #
def block_perm(x, y, g, n_perm, seed):
    """p non centrato e centrato di `scpipe.block_permutation_test` (usato da 17 e 18)."""
    r = scpipe.block_permutation_test(x, y, g, n_perm=n_perm, seed=seed)
    return r["p"], r["p_centered"]


@pytest.mark.slow
def test_block_permutation_type1_with_patient_confounder():
    rng = np.random.default_rng(2024)
    sizes = SIZES[:10]
    g = np.repeat(np.array([f"P{i:02d}" for i in range(10)]), sizes)
    pos = pos_c = 0
    reps = 400
    for r in range(reps):
        z = rng.normal(0, 1, 10)
        gi = np.repeat(np.arange(10), sizes)
        x = z[gi] + rng.normal(0, 1, len(gi))  # log dimensione clonale
        pr = 1 / (1 + np.exp(-(-1.5 + z[gi])))  # errore: dipende solo dal paziente
        y = (rng.random(len(gi)) < pr).astype(float)
        p, p_c = block_perm(x, y, g, n_perm=2000, seed=r)
        pos += p < 0.05
        pos_c += p_c < 0.05
    rate, rate_c = pos / reps, pos_c / reps
    assert 0.02 <= rate <= 0.08, f"falsi positivi {rate:.3f} fuori banda [0.02, 0.08]"
    assert 0.02 <= rate_c <= 0.08, f"falsi positivi (centrato) {rate_c:.3f} fuori banda"


def test_permutation_corr_type1_n10():
    """Banda 7: frazione di p < 0.05 sotto H0 con n = 10 in [0.02, 0.08]."""
    rng = np.random.default_rng(77)
    reps, pos = 400, 0
    for r in range(reps):
        x, y = rng.normal(size=10), rng.normal(size=10)
        _, p = scpipe.permutation_corr_test(x, y, n_perm=2000, seed=r)
        pos += p < 0.05
    rate = pos / reps
    assert 0.02 <= rate <= 0.08, f"falsi positivi {rate:.3f} fuori banda [0.02, 0.08]"


# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("script", ["16_pseudobulk_patient", "19_ni_score"])
def test_holm_monotone_and_bounded(script):
    holm = load_script(script).holm_bonferroni
    rng = np.random.default_rng(0)
    for _ in range(200):
        m = rng.integers(1, 8)
        raw = {f"h{i}": float(v) for i, v in enumerate(rng.random(m) ** 2)}
        adj = holm(raw)
        assert set(adj) == set(raw)
        order = sorted(raw, key=raw.get)
        vals = [adj[k] for k in order]
        assert all(a <= b + 1e-15 for a, b in zip(vals, vals[1:]))
        assert all(adj[k] >= raw[k] - 1e-15 for k in raw)
        assert all(0 <= adj[k] <= 1 for k in raw)


@pytest.mark.parametrize("script", ["16_pseudobulk_patient", "19_ni_score"])
def test_holm_hand_example(script):
    holm = load_script(script).holm_bonferroni
    adj = holm({"a": 0.01, "b": 0.04, "c": 0.03, "d": 0.5})
    # ordinati: a 0.01*4=0.04, c 0.03*3=0.09, b 0.04*2=0.08 -> max(0.09,0.08)=0.09, d 0.5*1 -> 0.5
    assert adj == pytest.approx({"a": 0.04, "c": 0.09, "b": 0.09, "d": 0.5})


# --------------------------------------------------------------------------- #
LYM = ["CD4T", "CD8T", "NK"]


def compositions(n: int, k: int = 3):
    for c in product(range(n + 1), repeat=k):
        if sum(c) == n:
            yield np.array(c)


def multinom_pmf(c, p) -> float:
    n = int(c.sum())
    coef = factorial(n)
    for ci in c:
        coef //= factorial(int(ci))
    return float(coef * np.prod(np.asarray(p) ** c))


def labels(c):
    return np.repeat(LYM, c)


@pytest.mark.parametrize("n_a,n_b", [(2, 2), (2, 5), (3, 4)])
@pytest.mark.parametrize("p", [(0.7, 0.2, 0.1), (0.4, 0.4, 0.2)])
def test_pairwise_discordance_unbiased(n_a, n_b, p):
    m = load_script("07_robustness")
    target = 1 - sum(pi**2 for pi in p)
    e_within = e_across = e_excess = 0.0
    for ca in compositions(n_a):
        for cb in compositions(n_b):
            w = multinom_pmf(ca, p) * multinom_pmf(cb, p)
            df = pd.DataFrame(
                {
                    "patient": "P1",
                    "clone_b": "X",
                    "tissue": ["A"] * n_a + ["B"] * n_b,
                    "all_celltype": np.r_[labels(ca), labels(cb)],
                }
            )
            r = m.pairwise_discordance(df, "clone_b").iloc[0]
            e_within += w * r.d_within
            e_across += w * r.d_across
            e_excess += w * r.excess
    assert e_within == pytest.approx(target, abs=1e-12)
    assert e_across == pytest.approx(target, abs=1e-12)
    assert e_excess == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("n", [2, 3, 6])
def test_d_within_08_unbiased(n):
    m = load_script("08_disambiguate")
    p = (0.6, 0.3, 0.1)
    e = sum(multinom_pmf(c, p) * m.d_within(labels(c)) for c in compositions(n))
    assert e == pytest.approx(1 - sum(pi**2 for pi in p), abs=1e-12)
    # il plug-in con reinserimento sarebbe distorto di un fattore (1 - 1/n)
    plug = sum(multinom_pmf(c, p) * (1 - ((c / n) ** 2).sum()) for c in compositions(n))
    assert plug == pytest.approx((1 - sum(pi**2 for pi in p)) * (1 - 1 / n), abs=1e-12)


# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_nadeau_bengio_calibration_and_power():
    sim = load_script("sim_a_calibration")
    for k in (3, 5):
        p0 = np.array(
            [sim.nb_replicate(r, 55_000 + k, 90, k, False)[0] for r in range(400)]
        )
        p0 = p0[~np.isnan(p0)]
        fpr = float(np.mean(p0 < 0.05))
        assert 0.0 <= fpr <= 0.08, (
            f"k = {k}: falsi positivi {fpr:.3f} fuori banda [0, 0.08]"
        )
    p1 = np.array([sim.nb_replicate(r, 51_005, 120, 5, True)[0] for r in range(100)])
    power = float(np.mean(p1[~np.isnan(p1)] < 0.05))
    assert power >= 0.5, f"potenza {power:.3f} < 0.50"


def test_majority_label_ties_first_appearance():
    assert scpipe.majority_label(["NK", "CD8T", "CD8T", "NK"]) == "NK"
    assert scpipe.majority_label(["CD8T", "NK", "NK", "CD8T"]) == "CD8T"
    assert scpipe.majority_label(["CD4T", "CD8T", "CD8T"]) == "CD8T"
    # pareggio fra 4 etichette distinte: vince quella che compare per prima
    lab = ["Endothelial", "CD8T", "Ductal", "Ductal", "Ductal", "NK", "NK", "NK", "B",
           "CD4T", "Endothelial", "B", "CD8T", "CD4T"]  # fmt: skip
    assert scpipe.majority_label(lab) == "Ductal"

    def reference(labels):
        counts = {}
        for x in labels:
            counts[x] = counts.get(x, 0) + 1
        top = max(counts.values())
        return next(x for x in labels if counts[x] == top)

    rng = np.random.default_rng(0)
    names = ["CD4T", "CD8T", "NK", "B", "Ductal", "Myeloid"]
    for _ in range(500):
        lab = list(rng.choice(names, rng.integers(1, 15)))
        assert scpipe.majority_label(lab) == reference(lab)


def test_macro_f1_uses_classes_present_in_test():
    y_true = np.array(["A", "A", "B", "B"])
    y_pred = np.array(["A", "B", "B", "C"])
    # classi presenti: A (P=1, R=0.5 -> F1 2/3), B (P=0.5, R=0.5 -> F1 0.5); C esclusa
    assert scpipe.macro_f1(y_true, y_pred) == pytest.approx((2 / 3 + 0.5) / 2)


def test_expanded_alpha():
    from scipy.stats import norm, t

    # formula esplicita per k = 10
    a10 = 2 * norm.cdf(np.sqrt(10 / 9) * t.ppf(0.025, 9))
    assert scpipe.expanded_alpha(10) == pytest.approx(a10)
    assert scpipe.expanded_alpha(10) < 0.05
    # per k grande l'intervallo torna al percentile semplice
    assert scpipe.expanded_alpha(10_000) == pytest.approx(0.05, abs=1e-3)
    stats = np.arange(1001, dtype=float)
    lo, hi = scpipe.percentile_ci(stats, 10_000)
    assert lo == pytest.approx(25, abs=0.5) and hi == pytest.approx(975, abs=0.5)


def test_cluster_bootstrap_edge_cases():
    m, lo, hi = scpipe.cluster_bootstrap({"P1": np.array([1.0, 0.0])})
    assert m == 0.5 and np.isnan(lo) and np.isnan(hi)
    m, lo, hi = scpipe.cluster_bootstrap({})
    assert np.isnan(m)
    # pazienti vuoti ignorati
    m, lo, hi = scpipe.cluster_bootstrap(
        {"P1": np.array([1.0]), "P2": np.array([0.0]), "P3": np.array([])}, n_boot=200
    )
    assert m == 0.5 and 0.0 <= lo <= hi <= 1.0


def test_small_k_diagnostics():
    # 4 pazienti: al massimo C(2k-1, k) = 35 multiinsiemi distinti; con 2.000 repliche
    # compaiono tutti. Differenze tutte positive -> frazione di repliche > 0 = 1.
    by = {f"P{i}": np.array([0.1 * (i + 1)]) for i in range(4)}
    d = scpipe.small_k_diagnostics(by, n_boot=2000, seed=0)
    assert d["n_distinct_resamples"] == 35
    assert d["frac_boot_gt0"] == 1.0
    # stesse estrazioni di cluster_bootstrap: stesso seme -> stesse repliche
    stats, _ = scpipe.cluster_bootstrap_replicates(by, 2000, 0)
    assert np.mean(stats > 0) == d["frac_boot_gt0"]
    assert np.isnan(
        scpipe.small_k_diagnostics({"P0": np.array([1.0])})["frac_boot_gt0"]
    )


def test_permutation_tests_constant_input_are_nan():
    # esito costante: rho non definito, nessun p minimo spurio
    with pytest.warns(UserWarning):
        rho, p = scpipe.permutation_corr_test(np.arange(10.0), np.zeros(10), n_perm=200)
    assert np.isnan(rho) and np.isnan(p)
    with pytest.warns(UserWarning):
        r = scpipe.block_permutation_test(
            np.arange(10.0), np.zeros(10), np.repeat(["A", "B"], 5), n_perm=200
        )
    assert all(np.isnan(r[k]) for k in ("rho", "p", "null_mean", "p_centered"))


def test_holm_excludes_undefined_p():
    # un p NaN non diventa 0 e non riduce il fattore delle altre ipotesi a caso
    with pytest.warns(UserWarning):
        adj = scpipe.holm_bonferroni({"a": np.nan, "b": 0.01, "c": 0.04})
    assert np.isnan(adj["a"])
    assert adj["b"] == pytest.approx(0.02) and adj["c"] == pytest.approx(0.04)
    assert list(adj) == ["a", "b", "c"]

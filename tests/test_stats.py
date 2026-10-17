"""Test di calibrazione degli stimatori statistici, su dati simulati con semi fissi.

Bande dichiarate PRIMA dell'esecuzione (non si allargano: se un test esce dalla
banda si corregge lo stimatore):

1. Cluster bootstrap sui pazienti (IC 95% della media pooled), con 10 e con 11
   pazienti di numerosita' diseguale come nei dati (5-59 cloni per paziente),
   esito continuo e binario: copertura empirica dell'IC in [0.90, 0.99].
   400 repliche per scenario.
2. Block permutation test (permutazione dell'esito entro paziente) sotto H0 con
   confondente di paziente (x e y dipendono entrambi da un effetto di paziente,
   nessuna associazione entro paziente): frazione di p < 0.05 in [0.02, 0.08].
   400 repliche, 2.000 permutazioni ciascuna.
3. Holm-Bonferroni: p corretti monotoni nell'ordine dei p grezzi, >= p grezzi, <= 1
   (proprieta' esatte, nessuna banda).
4. Discordanza a coppie senza reinserimento: valore atteso esatto (enumerazione di
   tutte le composizioni multinomiali) uguale a 1 - sum_k p_k^2; eccesso atteso
   esattamente 0 quando i due compartimenti hanno la stessa distribuzione
   (tolleranza numerica 1e-12, nessuna banda).
"""

from __future__ import annotations

from itertools import product
from math import factorial

import numpy as np
import pandas as pd
import pytest

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
    m = load_script("17_clonal_expansion")
    return m.block_permutation_test(x, y, g, n_perm=n_perm, seed=seed)[:2]


@pytest.mark.slow
def test_block_permutation_type1_with_patient_confounder():
    rng = np.random.default_rng(2024)
    sizes = SIZES[:10]
    g = np.repeat(np.array([f"P{i:02d}" for i in range(10)]), sizes)
    pos = 0
    reps = 400
    for r in range(reps):
        z = rng.normal(0, 1, 10)
        gi = np.repeat(np.arange(10), sizes)
        x = z[gi] + rng.normal(0, 1, len(gi))  # log dimensione clonale
        pr = 1 / (1 + np.exp(-(-1.5 + z[gi])))  # errore: dipende solo dal paziente
        y = (rng.random(len(gi)) < pr).astype(float)
        _, p = block_perm(x, y, g, n_perm=2000, seed=r)
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

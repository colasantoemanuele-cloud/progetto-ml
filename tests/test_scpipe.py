"""Test delle formule di `src/scpipe.py`.

Famiglia "Formule": confronto con calcoli a mano su esempi piccoli.
Nessuna banda statistica: i confronti sono esatti (tolleranza numerica 1e-9).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.stats import t as student_t

import scpipe

# --------------------------------------------------------------------------- #
# Nadeau-Bengio


def test_nadeau_bengio_hand_computed():
    # d = a - b = [0.02, 0.02, -0.01, 0.05, 0.02]; media 0.02
    # scarti dalla media: [0, 0, -0.03, 0.03, 0] -> SS = 0.0018, var (ddof=1) = 0.00045
    # k = 5, n_test / n_train = 200 / 800 = 0.25 -> fattore 1/5 + 0.25 = 0.45
    # t = 0.02 / sqrt(0.45 * 0.00045) = 1.405456..., gdl = 4
    a = [0.80, 0.82, 0.79, 0.85, 0.81]
    b = [0.78, 0.80, 0.80, 0.80, 0.79]
    t_hand = 0.02 / math.sqrt(0.45 * 0.00045)
    p_hand = 2 * student_t.sf(t_hand, df=4)
    assert t_hand == pytest.approx(1.4054567, abs=1e-6)
    assert p_hand == pytest.approx(0.2325953, abs=1e-6)
    assert scpipe.nadeau_bengio_ttest(a, b, n_train=800, n_test=200) == pytest.approx(
        p_hand, abs=1e-9
    )


def test_nadeau_bengio_symmetric_in_sign():
    a, b = [0.9, 0.8, 0.85, 0.7], [0.88, 0.81, 0.80, 0.69]
    p1 = scpipe.nadeau_bengio_ttest(a, b, 100, 25)
    p2 = scpipe.nadeau_bengio_ttest(b, a, 100, 25)
    assert p1 == pytest.approx(p2, abs=1e-12)


def test_nadeau_bengio_zero_variance_nonzero_mean():
    # differenza costante e non nulla in tutti i fold: varianza zero, effetto certo -> p = 0
    a = [0.9, 0.8, 0.7, 0.6, 0.5]
    b = [0.8, 0.7, 0.6, 0.5, 0.4]
    assert scpipe.nadeau_bengio_ttest(a, b, 800, 200) == 0.0


def test_nadeau_bengio_identical_models():
    # differenze tutte nulle: nessuna evidenza -> p = 1
    a = [0.9, 0.8, 0.7]
    assert scpipe.nadeau_bengio_ttest(a, a, 800, 200) == 1.0


def test_nadeau_bengio_single_fold():
    assert scpipe.nadeau_bengio_ttest([0.9], [0.8], 800, 200) == 1.0


# --------------------------------------------------------------------------- #
# Jaccard


def test_jaccard_basic():
    assert scpipe.jaccard({1, 2, 3}, {2, 3, 4}) == pytest.approx(0.5)
    assert scpipe.jaccard({1, 2}, {1, 2}) == 1.0
    assert scpipe.jaccard({1}, {2}) == 0.0


def test_jaccard_empty_sets():
    # due insiemi vuoti sono identici: J = 1 per convenzione
    assert scpipe.jaccard(set(), set()) == 1.0
    assert scpipe.jaccard(set(), {1}) == 0.0


def test_mean_pairwise_jaccard():
    sets = [{1, 2}, {1, 2}, {3}]
    # coppie: (0,1) = 1, (0,2) = 0, (1,2) = 0 -> media 1/3
    assert scpipe.mean_pairwise_jaccard(sets) == pytest.approx(1 / 3)
    assert np.isnan(scpipe.mean_pairwise_jaccard([{1}]))


# --------------------------------------------------------------------------- #
# LogCPM (log1p di CP10k)


def test_logcpm_row_sums_before_log():
    rng = np.random.default_rng(0)
    D = rng.poisson(1.0, size=(30, 50)).astype(np.float32)
    D[3] = 0  # cellula vuota: deve restare vuota, senza divisioni per zero
    X = sp.csr_matrix(D)
    Y = scpipe.LogCPM().fit(X).transform(X)
    raw = Y.copy()
    raw.data = np.expm1(raw.data)  # inverte log1p -> CP10k
    sums = np.asarray(raw.sum(axis=1)).ravel()
    nonempty = np.asarray(X.sum(axis=1)).ravel() > 0
    np.testing.assert_allclose(sums[nonempty], 1e4, rtol=1e-4)
    assert sums[3] == 0


def test_logcpm_hand_example():
    X = sp.csr_matrix(np.array([[1.0, 3.0, 0.0], [2.0, 0.0, 2.0]]))
    Y = scpipe.LogCPM(target_sum=1e4).transform(X).toarray()
    expected = np.log1p(np.array([[2500.0, 7500.0, 0.0], [5000.0, 0.0, 5000.0]]))
    np.testing.assert_allclose(Y, expected, rtol=1e-6)


def test_logcpm_is_stateless():
    X1 = sp.csr_matrix(np.array([[1.0, 1.0]]))
    X2 = sp.csr_matrix(np.array([[100.0, 0.0]]))
    a = scpipe.LogCPM().fit(X1).transform(X1).toarray()
    b = scpipe.LogCPM().fit(X2).transform(X1).toarray()
    np.testing.assert_array_equal(a, b)


# --------------------------------------------------------------------------- #
# HVGSelector (rapporto varianza / media)


def test_hvg_selects_highest_dispersion():
    # colonne con media nota e varianza nota; dispersione = var / media
    X = np.array(
        [
            [1, 0, 2, 5, 0],
            [1, 4, 2, 5, 0],
            [1, 0, 2, 5, 0],
            [1, 4, 2, 5, 0],
        ],
        dtype=float,
    )
    # colonna 0: media 1, var 0 -> 0
    # colonna 1: media 2, var 4 -> 2
    # colonna 2: media 2, var 0 -> 0
    # colonna 3: media 5, var 0 -> 0
    # colonna 4: media 0 -> esclusa (sotto min_mean)
    X[0, 2], X[1, 2] = 0.0, 4.0  # colonna 2: [0, 4, 2, 2] media 2, var 2 -> 1
    sel = scpipe.HVGSelector(n_top=2).fit(sp.csr_matrix(X))
    assert list(sel.idx_) == [1, 2]
    out = sel.transform(sp.csr_matrix(X))
    assert out.shape == (4, 2)


def test_hvg_excludes_low_mean_and_caps_k():
    X = sp.csr_matrix(np.array([[0.0, 1.0], [0.0, 3.0]]))
    sel = scpipe.HVGSelector(n_top=10).fit(X)
    assert list(sel.idx_) == [1]


def test_hvg_matches_numpy_dense():
    rng = np.random.default_rng(1)
    D = rng.poisson(rng.gamma(1.0, 1.0, size=40), size=(200, 40)).astype(float)
    sel = scpipe.HVGSelector(n_top=7).fit(sp.csr_matrix(D))
    mean, var = D.mean(0), D.var(0)
    disp = np.where(mean > 1e-3, var / (mean + 1e-12), -np.inf)
    assert list(sel.idx_) == sorted(np.argsort(-disp)[:7])

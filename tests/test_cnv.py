"""Test del filtro del rumore globale e delle coordinate delle finestre CNV (13, 14).

Famiglia "Formule": confronti esatti su dati sintetici, nessuna banda statistica.
- `global_noise_filter`: soglia = 1,5 x sd di tutta la matrice (come `np.std`), stesso
  risultato su matrice densa e sparsa, indipendente dall'ordine delle righe.
- `window_table`: le finestre ricostruite dai geni coincidono con quelle di infercnvpy
  (una cellula con un solo gene alterato cambia solo le finestre che contengono quel
  gene) e `locus_windows` / `check_marker_windows` di 14 attribuiscono al locus le
  finestre che contengono il gene marcatore.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import load_script


def test_global_noise_filter_matches_numpy_and_is_order_free():
    m = load_script("13_infercnv")
    rng = np.random.default_rng(0)
    D = rng.normal(0, 1, size=(300, 40))
    D[rng.random(D.shape) < 0.3] = 0.0
    out, thr = m.global_noise_filter(D)
    assert thr == pytest.approx(1.5 * np.std(D), rel=1e-12)
    expected = np.where(np.abs(D) < 1.5 * np.std(D), 0.0, D)
    np.testing.assert_array_equal(out, expected)
    # stessa soglia e stesso risultato su CSR (zeri impliciti inclusi nella sd)
    out_sp, thr_sp = m.global_noise_filter(sp.csr_matrix(D))
    assert thr_sp == pytest.approx(thr, rel=1e-12)
    np.testing.assert_array_equal(out_sp.toarray(), expected)
    # permutare le righe permuta il risultato (nessuna dipendenza dall'ordine)
    perm = rng.permutation(D.shape[0])
    out_p, thr_p = m.global_noise_filter(D[perm])
    assert thr_p == pytest.approx(thr, rel=1e-12)
    np.testing.assert_array_equal(out_p, expected[perm])
    # la matrice di ingresso non viene modificata
    assert np.count_nonzero(D) > np.count_nonzero(out)


def _synthetic_cnv_adata(rng):
    """Due cromosomi con geni a posizioni irregolari; 40 immuni di riferimento e 1 cellula
    con il solo gene `G_TARGET` molto espresso."""
    n1, n2 = 57, 33
    starts = np.r_[
        np.sort(rng.integers(0, 50_000_000, n1)), np.sort(rng.integers(0, 9e6, n2))
    ]
    chrom = ["chr1"] * n1 + ["chr2"] * n2
    names = [f"G{i}" for i in range(n1 + n2)]
    target = names[20]
    X = np.full((41, n1 + n2), 1.0)
    X[-1, names.index(target)] = 50.0
    a = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"ref_group": ["CD4T"] * 40 + ["Ductal"]},
            index=[f"c{i}" for i in range(41)],
        ),
        var=pd.DataFrame(
            {
                "chromosome": chrom,
                "start": starts.astype(float),
                "end": starts + 5_000.0,
            },
            index=names,
        ),
    )
    return a, target


def test_window_table_matches_infercnvpy_windows():
    import infercnvpy as cnv

    m = load_script("13_infercnv")
    a, target = _synthetic_cnv_adata(np.random.default_rng(1))
    cnv.tl.infercnv(
        a,
        reference_key="ref_group",
        reference_cat=["CD4T"],
        window_size=10,
        dynamic_threshold=None,
        exclude_chromosomes=None,
        n_jobs=1,
    )
    w = m.window_table(a, window_size=10)
    X = a.obsm["X_cnv"].toarray()
    assert len(w) == X.shape[1]
    # finestre che contengono il gene alterato secondo la tabella
    gs = a.var.loc[target, "start"]
    contains = ((w.chrom == "1") & (w.start <= gs) & (w.end >= gs)).values
    # finestre in cui la cellula alterata si discosta dal suo valore modale
    row = X[-1]
    vals, counts = np.unique(np.round(row, 10), return_counts=True)
    mode = vals[np.argmax(counts)]
    changed = np.abs(row - mode) > 1e-8
    np.testing.assert_array_equal(contains, changed)
    assert w.n_genes.max() == 10


def test_locus_windows_and_marker_check():
    m14 = load_script("14_cnv_loci")
    windows = pd.DataFrame(
        {
            "chrom": ["9", "9", "9", "17"],
            "start": [15e6, 21.5e6, 30e6, 1e6],
            "end": [21e6, 23e6, 36e6, 9e6],
        }
    )
    locus = {"name": "CDKN2A_9p21", "chrom": "9", "start": 20e6, "end": 24e6}
    np.testing.assert_array_equal(
        m14.locus_windows(windows, locus), [True, True, False, False]
    )
    genes = pd.DataFrame(
        {"start": [21_966_924.0], "end": [21_995_301.0]}, index=["CDKN2A"]
    )
    # CDKN2A e' nella seconda finestra, attribuita al locus: il controllo passa
    m14.check_marker_windows(windows, genes)
    # se nessuna finestra contiene il gene (coordinate sbagliate), lo script si ferma
    bad = windows.assign(start=[15e6, 25e6, 30e6, 1e6], end=[19e6, 29e6, 36e6, 9e6])
    with pytest.raises(SystemExit):
        m14.check_marker_windows(bad, genes)

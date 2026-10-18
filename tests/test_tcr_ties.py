"""Test della regola dei pareggi dell'etichetta di maggioranza in `src/10_loco.py`.

Famiglia "Formule": dati sintetici costruiti a mano, controlli esatti (nessuna banda).
Un'unita' clone-compartimento in pareggio non ha etichetta di maggioranza: e' esclusa
dal calcolo principale e, nella sensibilita', contata come errore oppure come corretta.
Il risultato non dipende dall'ordine delle cellule.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from conftest import load_script


def _cells(rng, patient, clone, tissue, labels, cd8):
    """Cellule di un clone in un compartimento: etichette e marcatori CD4/CD8 (log-scala)."""
    n = len(labels)
    return pd.DataFrame(
        {
            "patient": patient,
            "clone_b": clone,
            "tissue": tissue,
            "all_celltype": labels,
            "g_CD8A": np.full(n, 1.0 if cd8 else 0.0),
            "g_CD8B": np.zeros(n),
            "g_CD4": np.full(n, 0.0 if cd8 else 1.0),
            "u_CD8A": np.full(n, 2.0 if cd8 else 0.0),
        }
    )


def _toy(rng, order_seed=None):
    """3 pazienti x 8 cloni CD8 nel sangue; nel Tumor alcuni concordanti, alcuni sbagliati,
    uno in pareggio per paziente (2 CD8T + 2 CD4T)."""
    parts = []
    for p in ("PA01", "PA02", "PA03"):
        for c in range(8):
            clone = f"{p}_c{c}"
            parts.append(_cells(rng, p, clone, "PBMC", ["CD8T"] * 4, cd8=True))
            if c < 5:
                lab = ["CD8T"] * 4  # corretto
            elif c < 7:
                lab = ["CD4T"] * 4  # errore
            else:
                lab = ["CD8T", "CD8T", "CD4T", "CD4T"]  # pareggio
            parts.append(_cells(rng, p, clone, "Tumor", lab, cd8=True))
    df = pd.concat(parts, ignore_index=True)
    if order_seed is not None:
        df = df.sample(frac=1, random_state=order_seed).reset_index(drop=True)
    return df


def test_ties_excluded_and_order_independent():
    m = load_script("10_loco")
    rng = np.random.default_rng(0)
    res = []
    for order in (None, 1, 2):
        lym = _toy(rng, order)
        j = m.arbitrate(m.per_comp_stats(lym, None), "blood")
        t = j[j.tissue == "Tumor"]
        assert t.tie.sum() == 3  # un pareggio per paziente
        assert t.loc[t.tie, "error"].isna().all()
        main = m.untied(t)
        res.append(main.error.mean())
        sens = m.tie_sensitivity(j, "blood")
        val = sens.set_index(["version", "quantity"]).value
        # 6 errori su 21 unita' senza pareggi; con i pareggi 9/24 o 6/24
        assert np.isclose(val[("ties_excluded", "error_Tumor")], 6 / 21)
        assert np.isclose(val[("ties_as_error", "error_Tumor")], 9 / 24)
        assert np.isclose(val[("ties_as_correct", "error_Tumor")], 6 / 24)
    assert np.allclose(
        res, res[0]
    )  # stesso risultato con qualsiasi ordine delle cellule


def test_report_counts_ties():
    m = load_script("10_loco")
    lym = _toy(np.random.default_rng(0))
    j = m.arbitrate(m.per_comp_stats(lym, None), "blood")
    r = m.report(j, "test")
    row = r[r.tissue == "Tumor"].iloc[0]
    assert row.n == 21 and row.n_ties == 3
    assert np.isclose(row.error, 6 / 21)

"""Test anti-leakage sugli script che addestrano modelli (05, 11, 12, 13, 20) e sugli split (04).

Principio: ogni trasformazione stimata dai dati (normalizzazione con parametri,
selezione dei geni, SVD, scaling, classificatore) deve dipendere solo dalle cellule
di training del fold. Il test addestra due volte lo stesso fold: la seconda volta i
valori delle cellule di TEST sono sostituiti da conteggi completamente diversi.
Se il codice e' privo di leakage, i geni scelti e tutti i parametri stimati per quel
fold restano identici (confronto esatto, nessuna banda statistica).
Controllo positivo (ogni script): in un fold il cui TRAINING contiene le cellule
alterate almeno un parametro stimato cambia, quindi il confronto non e' banalmente vero.

Per 04 si verifica che gli split per paziente non mettano lo stesso paziente in due fold.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import load_script

N_PAT, N_PER, N_GENES, N_CLASSES = 6, 40, 80, 3


def make_data(seed: int = 0):
    """AnnData sintetico: 6 pazienti x 40 cellule, 3 classi, 80 geni, conteggi Poisson."""
    rng = np.random.default_rng(seed)
    patients = np.repeat([f"P{i}" for i in range(N_PAT)], N_PER)
    y = np.tile(np.array([f"C{i % N_CLASSES}" for i in range(N_PER)]), N_PAT)
    base = rng.gamma(1.0, 1.0, size=N_GENES)
    cls_shift = rng.gamma(2.0, 1.0, size=(N_CLASSES, N_GENES))
    pat_shift = rng.gamma(4.0, 0.25, size=(N_PAT, N_GENES))
    ci = np.array([int(c[1:]) for c in y])
    pi = np.array([int(p[1:]) for p in patients])
    lam = base * cls_shift[ci] * pat_shift[pi]
    X = sp.csr_matrix(rng.poisson(lam).astype(np.float32))
    names = [f"c{i}" for i in range(len(y))]
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=names))
    adata.var_names = [f"G{i}" for i in range(N_GENES)]
    adata.var["gene_symbol"] = adata.var_names.values
    folds_group = pi % 3  # 3 fold per paziente, 2 pazienti per fold
    sp_df = pd.DataFrame(
        {
            "task": "task2_v1",
            "cell_id": names,
            "patient": patients,
            "y": y,
            "fold_group_5": folds_group,
            "fold_lopo": pi,
            "fold_random_5": rng.integers(0, 3, len(y)),
        }
    )
    return adata, sp_df


def perturb_test(adata, mask, seed: int = 99):
    """Copia di adata con i conteggi delle cellule in `mask` sostituiti da valori diversi."""
    rng = np.random.default_rng(seed)
    X = adata.X.toarray()
    lam = rng.gamma(0.3, 20.0, size=X.shape[1])
    X[mask] = rng.poisson(lam, size=(mask.sum(), X.shape[1]))
    out = adata.copy()
    out.X = sp.csr_matrix(X.astype(np.float32))
    return out


def fitted_params(pipe) -> dict:
    """Estrae tutti i parametri stimati di una Pipeline (array numerici)."""
    out = {}
    for name, step in pipe.named_steps.items():
        for attr in (
            "idx_",
            "mean_",
            "scale_",
            "var_",
            "components_",
            "coef_",
            "intercept_",
            "class_prior_",
            "feature_log_prob_",
            "feature_importances_",
            "_fit_X",
        ):
            if hasattr(step, attr):
                v = getattr(step, attr)
                if v is not None:
                    out[f"{name}.{attr}"] = np.asarray(v)
    return out


def assert_same_params(p1: dict, p2: dict):
    assert p1.keys() == p2.keys()
    assert p1, "nessun parametro stimato trovato"
    for k in p1:
        np.testing.assert_array_equal(p1[k], p2[k], err_msg=k)


def params_differ(p1: dict, p2: dict) -> bool:
    """True se almeno un parametro stimato cambia (controllo positivo)."""
    return any(
        p1[k].shape != p2[k].shape or not np.array_equal(p1[k], p2[k]) for k in p1
    )


# --------------------------------------------------------------------------- #
def test_05_baseline_no_leakage(monkeypatch):
    m = load_script("05_baseline")
    built = []
    orig = m.make_pipeline

    def rec(model, n_classes):
        p = orig(model, n_classes)
        built.append((model, p))
        return p

    monkeypatch.setattr(m, "make_pipeline", rec)
    adata, sp_df = make_data()
    folds = sp_df.fold_group_5.values
    nf = len(np.unique(folds))

    m.run_task(adata, sp_df, "task2_v1", quick=False, n_hvg=20)
    first = list(built)
    built.clear()
    m.run_task(
        perturb_test(adata, folds == 0), sp_df, "task2_v1", quick=False, n_hvg=20
    )
    second = list(built)

    # ordine dei cicli in run_task: schema primario -> (dummy, logreg_l2) -> fold
    i = nf + 0  # logreg_l2, fold 0, schema primario
    assert first[i][0] == "logreg_l2"
    assert_same_params(fitted_params(first[i][1]), fitted_params(second[i][1]))
    # controllo positivo: i fold il cui training contiene le cellule alterate cambiano
    p_other1, p_other2 = (
        fitted_params(first[nf + 1][1]),
        fitted_params(second[nf + 1][1]),
    )
    assert any(not np.array_equal(p_other1[k], p_other2[k]) for k in p_other1)


@pytest.mark.parametrize("model", ["logreg_l2", "knn", "random_forest", "naive_bayes"])
def test_11_model_zoo_no_leakage(monkeypatch, model):
    m = load_script("11_model_zoo")
    monkeypatch.setattr(m, "N_HVG", 30)
    monkeypatch.setattr(m, "N_SVD", 5)
    built = []
    orig = m.build

    def rec(name, dense):
        p = orig(name, dense)
        built.append(p)
        return p

    monkeypatch.setattr(m, "build", rec)
    adata, sp_df = make_data()
    folds = sp_df.fold_lopo.values
    m.run_task(adata, sp_df, "task2_v1", [model], 0, "fold_lopo")
    first = list(built)
    built.clear()
    m.run_task(
        perturb_test(adata, folds == 0), sp_df, "task2_v1", [model], 0, "fold_lopo"
    )
    second = list(built)
    assert_same_params(fitted_params(first[0]), fitted_params(second[0]))
    # controllo positivo: il fold 1 ha le cellule alterate (paziente del fold 0) nel training
    assert params_differ(fitted_params(first[1]), fitted_params(second[1]))


def test_12_feature_selection_no_leakage():
    m = load_script("12_feature_selection")
    adata, sp_df = make_data()
    X, y = adata.X, sp_df.y.values
    folds = sp_df.fold_group_5.values
    tr, te = np.where(folds != 0)[0], np.where(folds == 0)[0]
    genes = adata.var_names.values
    classes = np.array(sorted(np.unique(y)))
    X2 = perturb_test(adata, folds == 0).X
    out1 = m.run_fold(
        X[tr],
        y[tr],
        X[te],
        y[te],
        genes,
        classes,
        list(m.SELECTORS),
        np.random.default_rng(0),
        k_values=(5, 10),
    )
    out2 = m.run_fold(
        X2[tr],
        y[tr],
        X2[te],
        y[te],
        genes,
        classes,
        list(m.SELECTORS),
        np.random.default_rng(0),
        k_values=(5, 10),
    )
    assert out1["chosen"].keys() == out2["chosen"].keys()
    for key in out1["chosen"]:
        assert out1["chosen"][key] == out2["chosen"][key], key
    # controllo positivo: nel fold 1 le cellule alterate sono nel training
    tr1, te1 = np.where(folds != 1)[0], np.where(folds == 1)[0]
    outs = [
        m.run_fold(
            Xa[tr1],
            y[tr1],
            Xa[te1],
            y[te1],
            genes,
            classes,
            list(m.SELECTORS),
            np.random.default_rng(0),
            k_values=(5, 10),
        )
        for Xa in (X, X2)
    ]
    assert any(
        outs[0]["chosen"][key] != outs[1]["chosen"][key] for key in outs[0]["chosen"]
    )


def test_13_infercnv_no_leakage():
    m = load_script("13_infercnv")
    rng = np.random.default_rng(0)
    adata, _ = make_data()
    groups = np.repeat([f"P{i}" for i in range(N_PAT)], N_PER)
    yb = rng.integers(0, 2, len(groups))
    X = adata.X
    g0 = groups == "P0"
    X2 = perturb_test(adata, g0).X
    _, _, models1 = m.classify_lopo(X, yb, groups, n_hvg=20, return_models=True)
    _, _, models2 = m.classify_lopo(X2, yb, groups, n_hvg=20, return_models=True)
    assert_same_params(fitted_params(models1["P0"]), fitted_params(models2["P0"]))
    # controllo positivo: il modello di P1 ha P0 (alterato) nel training
    assert params_differ(fitted_params(models1["P1"]), fitted_params(models2["P1"]))
    # il modello finale per Adjacent_normal sceglie i geni sul solo training:
    # coincidono con quelli di un HVGSelector stimato a parte sulle sole duttali tumorali,
    # e applicarlo a cellule esterne diverse non li modifica
    scpipe = load_script("scpipe")
    final = m.fit_final(X, yb, n_hvg=20)
    ref = scpipe.HVGSelector(n_top=20).fit(scpipe.LogCPM().transform(X))
    np.testing.assert_array_equal(final.named_steps["hvg"].idx_, ref.idx_)
    before = fitted_params(final)
    final.predict(X2)
    assert_same_params(before, fitted_params(final))


def test_20_xai_no_leakage():
    m = load_script("20_xai_per_patient")
    adata, sp_df = make_data()
    folds = sp_df.fold_lopo.values
    patients = sp_df.patient.values
    genes = adata.var_names.values
    y = np.where(sp_df.y.values == "C0", "Tumor", "Adjacent_normal")
    X2 = perturb_test(adata, folds == 0).X
    r1 = m.lopo_explanations(
        adata.X, y, folds, patients, genes, np.random.default_rng(0)
    )
    r2 = m.lopo_explanations(X2, y, folds, patients, genes, np.random.default_rng(0))
    coefs1, tops1 = r1[0], r1[1]
    coefs2, tops2 = r2[0], r2[1]
    assert tops1[0] == tops2[0]
    np.testing.assert_array_equal(coefs1[0].values, coefs2[0].values)
    # controllo positivo: il fold 1 ha le cellule alterate nel training
    assert not np.array_equal(coefs1[1].values, coefs2[1].values)


# --------------------------------------------------------------------------- #
def test_04_splits_patient_disjoint():
    m = load_script("04_splits")
    rng = np.random.default_rng(0)
    n = 900
    patients = rng.choice([f"PA{i:02d}" for i in range(1, 14)], n)
    y = rng.choice(["A", "B", "C"], n)
    df = pd.DataFrame({"patient": patients, "y": y}, index=[f"c{i}" for i in range(n)])
    out = m.assign_folds(df)
    for scheme in ("fold_group_5", "fold_lopo"):
        assert (out[scheme] >= 0).all()
        per_patient = out.groupby("patient")[scheme].nunique()
        assert (per_patient == 1).all(), scheme
    assert out.fold_group_5.nunique() == 5
    assert out.fold_lopo.nunique() == 13
    # controllo positivo: lo split casuale mette ogni paziente in piu' fold (il controllo puo' fallire)
    assert (out.groupby("patient").fold_random_5.nunique() > 1).all()

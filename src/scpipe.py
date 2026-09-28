#!/usr/bin/env python3
"""
scpipe.py: Trasformatori Scikit-learn esenti da data leakage e test statistici per la validazione incrociata.

Architettura Anti-Leakage
-------------------------
Questo modulo fornisce i componenti condivisi da `11_model_zoo.py`,
`12_feature_selection.py` e `20_xai_per_patient.py`. Ogni operazione di
pre-elaborazione che stima parametri statistici attraverso le osservazioni
(selezione dei geni altamente variabili tramite dispersione media-varianza,
standardizzazione, decomposizione SVD troncata) e' incapsulata in trasformatori
conformi all'interfaccia `BaseEstimator` / `TransformerMixin` di Scikit-learn,
garantendo che il metodo `fit()` venga invocato esclusivamente sul sottoinsieme
di training di ciascun fold.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import t as student_t
from sklearn.base import BaseEstimator, TransformerMixin


class LogCPM(BaseEstimator, TransformerMixin):
    """Trasformatore stateless per la normalizzazione di libreria per-cellula
    secondo la trasformazione `log1p(conteggi / totale_cellula * target_sum)`.
    Operando in modo indipendente su ciascuna riga sparsa CSR, non trasferisce
    informazione tra osservazioni di training e di test."""

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        # Calcolo della profondita' di sequenziamento totale per ciascuna cellula (riga)
        counts = np.asarray(X.sum(axis=1)).ravel()
        counts[counts == 0] = 1.0
        # Risigillatura in formato CSR e applicazione in-place di log1p sui soli elementi non nulli
        X = sp.csr_matrix(X.multiply(self.target_sum / counts[:, None]))
        X.data = np.log1p(X.data)
        return X


class HVGSelector(BaseEstimator, TransformerMixin):
    """Selettore non supervisionato dei `n_top` geni altamente variabili (HVG)
    basato sul rapporto di dispersione (varianza / media) calcolato in formato
    sparso esclusivamente sulle cellule del fold di training."""

    def __init__(self, n_top: int = 2000, min_mean: float = 1e-3):
        self.n_top = n_top
        self.min_mean = min_mean

    def fit(self, X, y=None):  # noqa: ARG002
        # Stima di media e momento secondo sulle colonne sparse senza densificazione
        mean = np.asarray(X.mean(axis=0)).ravel()
        sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = np.maximum(sq - mean**2, 0.0)
        # Calcolo dell'indice di Fano (var/mean) per i soli geni con media superiore a min_mean
        disp = np.where(mean > self.min_mean, var / (mean + 1e-12), -np.inf)
        k = min(self.n_top, int(np.isfinite(disp).sum()))
        self.idx_ = np.sort(np.argsort(-disp)[:k])
        return self

    def transform(self, X):
        return X[:, self.idx_]


class IndexSelector(BaseEstimator, TransformerMixin):
    """Trasformatore di proiezione colonnare che filtra una matrice sparsa CSR
    sulla base di un vettore di indici di feature pre-calcolato sul training fold."""

    def __init__(self, idx: np.ndarray | None = None):
        self.idx = idx

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        return X if self.idx is None else X[:, self.idx]


def nadeau_bengio_ttest(
    a: np.ndarray, b: np.ndarray, n_train: int, n_test: int
) -> float:
    """Esegue il test t appaiato corretto per la covarianza tra fold secondo
    Nadeau & Bengio (2003).

    Poiche' i training set dei k fold di una validazione incrociata condividono
    una frazione considerevole di osservazioni, le differenze di prestazione
    `d_i = a_i - b_i` non sono statisticamente indipendenti e il test t di Student
    standard sottostima l'errore standard producendo un eccesso di falsi positivi.
    Il fattore correttivo `(1/k + n_test / n_train)` compensa la correlazione
    indotta dalla sovrapposizione dei training set su `k - 1` gradi di liberta'.
    """
    d = np.asarray(a) - np.asarray(b)
    k = len(d)
    if k < 2 or np.allclose(d, 0):
        return 1.0
    # Calcolo della varianza campionaria delle differenze e applicazione del fattore (1/k + n_test/n_train)
    var = d.var(ddof=1)
    denom = np.sqrt((1.0 / k + n_test / max(n_train, 1)) * var)
    if denom == 0:
        return 1.0
    stat = d.mean() / denom
    return float(2 * (1 - student_t.cdf(abs(stat), df=k - 1)))


def jaccard(a: set, b: set) -> float:
    """Calcola il coefficiente di similarita' di Jaccard tra due insiemi di feature."""
    return len(a & b) / max(len(a | b), 1)


def mean_pairwise_jaccard(sets: list[set]) -> float:
    """Calcola la media dell'indice di Jaccard su tutte le coppie possibili di
    sottoinsiemi di feature selezionati nei diversi fold della cross-validation,
    quantificando la stabilita' strutturale della selezione."""
    if len(sets) < 2:
        return np.nan
    vals = [
        jaccard(sets[i], sets[j])
        for i in range(len(sets))
        for j in range(i + 1, len(sets))
    ]
    return float(np.mean(vals))

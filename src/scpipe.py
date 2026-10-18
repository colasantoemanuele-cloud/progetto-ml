#!/usr/bin/env python3
"""
scpipe.py: trasformatori Scikit-learn senza leakage e funzioni statistiche condivise.

Scopo
-----
Modulo infrastrutturale importato dagli script della pipeline (05-13, 16-20) e da
`src/sim/sim_a_calibration.py` e `src/sim/sim_c_block_permutation.py`
(`sim_b` e `sim_d` lo usano attraverso 10 e 20). Non e' uno script eseguibile e non
scrive file.

Contenuto
---------
Trasformatori (il metodo `fit()` va invocato solo sul training del fold):
  - `LogCPM`        : log1p(CP10k), cioe' log1p(conteggi / totale_cellula * 1e4).
                      Il nome storico "LogCPM" e' improprio: la scala e' 10^4
                      (CP10k), non 10^6 (CPM). Non stima parametri.
  - `HVGSelector`   : i `n_top` geni con rapporto varianza / media piu' alto,
                      stimato sulle sole cellule di training.

Metriche e test per il confronto di modelli:
  - `macro_f1`             : macro-F1 sulle sole classi presenti nel fold di test.
  - `nadeau_bengio_ttest`  : t appaiato corretto di Nadeau & Bengio (2003).
  - `jaccard`, `mean_pairwise_jaccard`.

Stimatori con il paziente come unita' (calibrazione in `tests/test_stats.py`):
  - `cluster_bootstrap`      : IC 95% della media pooled ricampionando i pazienti,
                               con intervallo "expanded percentile". La copertura e'
                               verificata solo per la media pooled con 10-11 pazienti
                               (banda [0.90, 0.99]); non e' verificata per altre
                               statistiche ricampionate con lo stesso intervallo
                               (Spearman in 16, 17, 19; d1/d2 in 10) ne' con 4-5
                               pazienti, dove in simulazione resta sotto il 95%
                               (`src/sim/sim_a_calibration.py`).
  - `small_k_diagnostics`    : con <= `SMALL_K` pazienti, ricampionature distinte e
                               frazione di repliche > 0 (descrittive, non un test).
  - `block_permutation_test` : correlazione di Spearman pooled, permutando l'esito
                               solo entro paziente; p grezzo e p sulla statistica centrata.

Test con unita' = righe (una riga per osservazione, es. un paziente in H1-H4):
  - `permutation_corr_test`  : Spearman con permutazione Monte Carlo; frazione di
                               p < 0.05 sotto H0 verificata in `tests/test_stats.py`.
  - `holm_bonferroni`        : correzione step-down di Holm.

Uso
---
    import sys; sys.path.insert(0, "src")
    from scpipe import LogCPM, HVGSelector, cluster_bootstrap
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm
from scipy.stats import t as student_t
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import f1_score


class LogCPM(BaseEstimator, TransformerMixin):
    """Normalizzazione per cellula log1p(conteggi / totale_cellula * target_sum).

    Con il default `target_sum = 1e4` la scala e' CP10k (non CPM, nonostante il
    nome). Opera riga per riga e non stima parametri: `fit()` non fa nulla, quindi
    non trasferisce informazione fra training e test. Le cellule con zero conteggi
    restano a zero.
    """

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        # Profondita' di sequenziamento di ciascuna cellula (somma della riga)
        counts = np.asarray(X.sum(axis=1)).ravel()
        counts[counts == 0] = 1.0
        # Riscalamento a target_sum e log1p sui soli elementi non nulli del CSR
        X = sp.csr_matrix(X.multiply(self.target_sum / counts[:, None]))
        X.data = np.log1p(X.data)
        return X


class HVGSelector(BaseEstimator, TransformerMixin):
    """Selettore dei `n_top` geni con rapporto varianza / media piu' alto.

    Media e varianza sono stimate in formato sparso sulle sole cellule passate a
    `fit()` (il training del fold). I geni con media <= `min_mean` sono esclusi.
    """

    def __init__(self, n_top: int = 2000, min_mean: float = 1e-3):
        self.n_top = n_top
        self.min_mean = min_mean

    def fit(self, X, y=None):  # noqa: ARG002
        # Media e momento secondo per colonna senza densificare la matrice
        mean = np.asarray(X.mean(axis=0)).ravel()
        sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = np.maximum(sq - mean**2, 0.0)
        # Indice di dispersione (var/mean) per i soli geni con media sopra min_mean
        disp = np.where(mean > self.min_mean, var / (mean + 1e-12), -np.inf)
        k = min(self.n_top, int(np.isfinite(disp).sum()))
        self.idx_ = np.sort(np.argsort(-disp)[:k])
        return self

    def transform(self, X):
        return X[:, self.idx_]


# --------------------------------------------------------------------------- #
# Metriche e confronto di modelli


def macro_f1(y_true, y_pred) -> float:
    """Macro-F1 calcolato sulle sole classi presenti nel fold di test (`y_true`).

    Una classe assente dal test non ha F1 definito: includerla con F1 = 0
    abbasserebbe artificialmente la media. Le predizioni di classi assenti dal
    test restano errori (abbassano il recall delle classi presenti).
    """
    labels = np.unique(y_true)
    return float(
        f1_score(y_true, y_pred, average="macro", labels=labels, zero_division=0)
    )


def nadeau_bengio_ttest(
    a: np.ndarray, b: np.ndarray, n_train: int, n_test: int
) -> float:
    """t appaiato corretto di Nadeau & Bengio (2003), p bilaterale.

    I training set dei k fold condividono gran parte delle osservazioni, quindi le
    differenze d_i = a_i - b_i non sono indipendenti e il t di Student standard
    sottostima l'errore standard. La varianza viene moltiplicata per
    (1/k + n_test / n_train) e il test usa k - 1 gradi di liberta'.

    Casi limite: k < 2 o differenze tutte nulle -> p = 1; differenze tutte uguali
    e non nulle (varianza zero) -> p = 0.
    """
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    k = len(d)
    if k < 2 or np.allclose(d, 0):
        return 1.0
    if np.allclose(d, d.mean()):
        return 0.0
    var = d.var(ddof=1)
    denom = np.sqrt((1.0 / k + n_test / max(n_train, 1)) * var)
    stat = d.mean() / denom
    return float(2 * student_t.sf(abs(stat), df=k - 1))


def jaccard(a: set, b: set) -> float:
    """Indice di Jaccard |A ∩ B| / |A ∪ B|; per due insiemi vuoti vale 1."""
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def mean_pairwise_jaccard(sets: list[set]) -> float:
    """Media dell'indice di Jaccard su tutte le coppie di insiemi (NaN se < 2 insiemi)."""
    if len(sets) < 2:
        return np.nan
    vals = [
        jaccard(sets[i], sets[j])
        for i in range(len(sets))
        for j in range(i + 1, len(sets))
    ]
    return float(np.mean(vals))


# --------------------------------------------------------------------------- #
# Stimatori con il paziente come unita'


def expanded_alpha(k: int, alpha: float = 0.05) -> float:
    """Livello corretto dell'intervallo "expanded percentile" (Hesterberg, 2015).

    Il percentile bootstrap con pochi cluster e' troppo stretto (in simulazione
    con 10-11 pazienti copre ~87-89% invece del 95%). L'intervallo espanso usa i
    quantili alpha'/2 e 1 - alpha'/2 della distribuzione bootstrap, con
        alpha'/2 = Phi( sqrt(k / (k - 1)) * t_{alpha/2, k-1} ),
    che compensa sia il fattore (k-1)/k della varianza plug-in sia la coda t.
    Restituisce alpha' (bilaterale).
    """
    if k < 2:
        return alpha
    return float(2 * norm.cdf(np.sqrt(k / (k - 1)) * student_t.ppf(alpha / 2, k - 1)))


def percentile_ci(stats: np.ndarray, k: int, alpha: float = 0.05) -> tuple:
    """IC "expanded percentile" dalle repliche bootstrap `stats` con k unita'."""
    stats = np.asarray(stats, dtype=float)
    stats = stats[np.isfinite(stats)]
    if len(stats) == 0:
        return np.nan, np.nan
    a = expanded_alpha(k, alpha)
    lo, hi = np.percentile(stats, [100 * a / 2, 100 * (1 - a / 2)])
    return float(lo), float(hi)


def cluster_bootstrap_replicates(
    by_pat: dict, n_boot: int = 2000, seed: int = 0
) -> tuple[np.ndarray, list]:
    """Repliche del cluster bootstrap: ogni replica estrae con reinserimento tanti
    pazienti quanti sono quelli disponibili, concatena le loro osservazioni e ne
    calcola la media. Restituisce (medie, estrazioni) dove ogni estrazione e' la
    lista dei pazienti estratti (serve a contare le ricampionature distinte)."""
    pats = list(by_pat)
    rng = np.random.default_rng(seed)
    stats, draws = [], []
    for _ in range(n_boot):
        draw = rng.choice(pats, len(pats), replace=True)
        stats.append(np.concatenate([by_pat[p] for p in draw]).mean())
        draws.append(draw)
    return np.asarray(stats, dtype=float), draws


def cluster_bootstrap(by_pat: dict, n_boot: int = 2000, seed: int = 0) -> tuple:
    """Media pooled e IC 95% con cluster bootstrap sui pazienti.

    `by_pat` mappa paziente -> array delle osservazioni (cloni, coppie, ...).
    Il punto e' la media di tutte le osservazioni concatenate; l'intervallo e'
    "expanded percentile" (vedi `expanded_alpha`) con k = numero di pazienti.
    Copertura verificata in `tests/test_stats.py` con 10 e 11 pazienti; con 4-5
    pazienti l'intervallo resta troppo stretto (vedi `src/sim/sim_a_calibration.py`).
    Restituisce (media, lo, hi); NaN se ci sono meno di 2 pazienti.
    """
    by_pat = {p: np.asarray(v, dtype=float) for p, v in by_pat.items() if len(v)}
    if len(by_pat) < 2:
        allv = np.concatenate(list(by_pat.values())) if by_pat else np.array([])
        return (float(allv.mean()) if len(allv) else np.nan, np.nan, np.nan)
    stats, _ = cluster_bootstrap_replicates(by_pat, n_boot, seed)
    lo, hi = percentile_ci(stats, len(by_pat))
    allv = np.concatenate(list(by_pat.values()))
    return float(allv.mean()), lo, hi


# Con un numero di pazienti <= SMALL_K l'IC del cluster bootstrap non e' calibrato
SMALL_K = 5


def small_k_diagnostics(by_pat: dict, n_boot: int = 2000, seed: int = 0) -> dict:
    """Diagnostica descrittiva del cluster bootstrap con pochi pazienti.

    Restituisce il numero di ricampionature distinte (come multiinsiemi di pazienti:
    con k pazienti sono al massimo C(2k-1, k), 35 con k = 4) e la frazione di
    repliche con media > 0. Usa lo stesso seme e la stessa sequenza di estrazioni di
    `cluster_bootstrap` con gli stessi `n_boot` e `seed`. La frazione > 0 ha senso
    solo per differenze con segno, non per tassi (sempre >= 0).
    """
    by_pat = {p: np.asarray(v, dtype=float) for p, v in by_pat.items() if len(v)}
    if len(by_pat) < 2:
        return {"n_distinct_resamples": np.nan, "frac_boot_gt0": np.nan}
    stats, draws = cluster_bootstrap_replicates(by_pat, n_boot, seed)
    distinct = {tuple(sorted(d)) for d in draws}
    return {
        "n_distinct_resamples": len(distinct),
        "frac_boot_gt0": float(np.mean(stats > 0)),
    }


def _corr_rows(rx: np.ndarray, ry_perm: np.ndarray) -> np.ndarray:
    """Correlazione di Pearson fra il vettore rx e ogni riga di ry_perm."""
    rx_c = rx - rx.mean()
    ry_c = ry_perm - ry_perm.mean(axis=1, keepdims=True)
    num = (ry_c * rx_c).sum(axis=1)
    den = np.sqrt((rx_c**2).sum()) * np.sqrt((ry_c**2).sum(axis=1))
    return num / np.maximum(den, 1e-12)


def permutation_corr_test(
    x: np.ndarray, y: np.ndarray, n_perm: int = 20000, seed: int = 0
) -> tuple[float, float]:
    """Spearman e p bilaterale con test di permutazione Monte Carlo.

    Le permutazioni sono `n_perm` estrazioni casuali (non l'enumerazione di tutte
    le n! permutazioni: il test non e' esatto). p con la correzione +1 di
    Phipson & Smyth (2010). Restituisce (rho, p). Se x o y e' costante rho non e'
    definito: restituisce (NaN, NaN) con un avviso, invece di un p minimo spurio.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        warnings.warn(
            "permutation_corr_test: x o y costante, rho non definito", stacklevel=2
        )
        return np.nan, np.nan
    rho_obs = float(np.corrcoef(rx, ry)[0, 1])
    perm_idx = np.argsort(rng.random((n_perm, n)), axis=1)
    rhos = _corr_rows(rx, ry[perm_idx])
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p)


def block_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    group: np.ndarray,
    n_perm: int = 20000,
    seed: int = 0,
) -> dict:
    """Spearman pooled con permutazione dell'esito solo entro ciascun gruppo (paziente).

    Sotto H0 (nessuna associazione entro paziente) le permutazioni entro blocco
    conservano l'associazione fra pazienti: la distribuzione nulla e' quindi
    centrata sulla sua media (`null_mean`), non necessariamente su zero.
    Restituisce un dizionario con:
      - rho        : Spearman pooled osservato;
      - p          : p bilaterale |rho_perm| >= |rho| (statistica non centrata);
      - null_mean  : media della distribuzione nulla;
      - p_centered : p bilaterale |rho_perm - null_mean| >= |rho - null_mean|.
    Entrambi i p sono Monte Carlo con correzione +1 (Phipson & Smyth).
    Calibrazione del p non centrato in `tests/test_stats.py`.
    Se x o y e' costante rho non e' definito: tutti i valori sono NaN, con un avviso.
    """
    rng = np.random.default_rng(seed)
    # Ordinamento stabile per paziente: ogni blocco occupa un intervallo contiguo
    order = np.argsort(group, kind="stable")
    x_s, g_s = x[order], group[order]
    ry_s = pd.Series(y[order]).rank().values
    rx_s = pd.Series(x_s).rank().values
    if np.ptp(rx_s) == 0 or np.ptp(ry_s) == 0:
        warnings.warn(
            "block_permutation_test: x o y costante, rho non definito", stacklevel=2
        )
        return {"rho": np.nan, "p": np.nan, "null_mean": np.nan, "p_centered": np.nan}
    rho_obs = float(np.corrcoef(rx_s, ry_s)[0, 1])

    n = len(x_s)
    perm_idx = np.tile(np.arange(n), (n_perm, 1))
    uniq, starts = np.unique(g_s, return_index=True)
    bounds = list(starts) + [n]
    for i in range(len(uniq)):
        s, e = bounds[i], bounds[i + 1]
        size = e - s
        if size <= 1:
            continue
        # Permutazione indipendente ristretta all'intervallo [s, e) del paziente
        keys = rng.random((n_perm, size))
        perm_idx[:, s:e] = s + np.argsort(keys, axis=1)

    rhos = _corr_rows(rx_s, ry_s[perm_idx])
    p = (np.sum(np.abs(rhos) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    c = float(rhos.mean())
    p_c = (np.sum(np.abs(rhos - c) >= abs(rho_obs - c) - 1e-12) + 1) / (n_perm + 1)
    return {"rho": rho_obs, "p": float(p), "null_mean": c, "p_centered": float(p_c)}


def holm_bonferroni(pvals: dict) -> dict:
    """Correzione step-down di Holm-Bonferroni (controllo forte del FWER).

    I p ordinati in modo crescente sono moltiplicati per (m - i), limitati a 1 e
    resi monotoni con il massimo cumulato. Un p non definito (NaN, es. covariata
    costante) e' escluso dalla famiglia, con un avviso: il suo p corretto e' NaN e
    m conta solo i test con p definito.
    """
    nan = [k for k, p in pvals.items() if p is None or np.isnan(p)]
    if nan:
        warnings.warn(
            f"holm_bonferroni: p non definito per {nan}, esclusi", stacklevel=2
        )
    items = sorted(
        ((k, p) for k, p in pvals.items() if k not in nan), key=lambda kv: kv[1]
    )
    m = len(items)
    adj, running_max = {k: np.nan for k in nan}, 0.0
    for i, (name, p) in enumerate(items):
        val = min((m - i) * p, 1.0)
        running_max = max(running_max, val)
        adj[name] = running_max
    return {k: adj[k] for k in pvals}


def majority_label(labels) -> str:
    """Etichetta di maggioranza (moda) di un insieme di cellule.

    Regola dei pareggi: a parita' di conteggio vince l'etichetta che compare per
    prima nell'ordine delle cellule passato in ingresso (per gli script TCR:
    l'ordine delle righe di `obs` dopo il merge con i contig). E' la stessa regola
    di `pd.Series.value_counts().index[0]` in pandas >= 2, usata nelle versioni
    precedenti degli script, resa esplicita e indipendente dalla versione.
    """
    labels = list(labels)
    vc = pd.Series(labels).value_counts()
    tied = set(vc.index[vc == vc.iloc[0]])
    return str(next(lab for lab in labels if lab in tied))

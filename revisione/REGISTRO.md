# Registro della revisione

Memoria di lavoro della revisione del ramo `revisione`. Ogni fase e ogni giro
del loop aggiornano questo file con stato, cosa resta da fare e comandi in corso.

Ambiente: `conda` env `pdac` (`/home/nemo/miniforge3/envs/pdac`), Python 3.11,
scikit-learn 1.9.0, scanpy 1.11.5, anndata 0.12.19, lightgbm 4.6.0, infercnvpy 0.6.1,
pytest 9.1.1, ruff 0.15.20. Macchina: 12 core, 30 GB RAM.

## Stato

| Fase | Stato | Commit |
| :--- | :--- | :--- |
| Verifica iniziale | fatta: remote `progetto-ml`, `sc_raw.h5ad` e 30 file VDJ presenti, MatrixMarket grezzi presenti | - |
| 0a Criteri | fatta | 4ef60ec |
| 0b Test | fatta | 12a356d |
| 1 Correzioni | fatta (A1-A7, B1-B7, C1-C4, D1-D5); pytest veloce 98 passati, ruff pulito | (commit Fase 1) |
| 2 Loop | da fare | - |
| 3 Esecuzione finale | da fare | - |

## Note iniziali

- La cartella locale si chiama `pdac-ml/`, il remote e la radice citata nel README
  devono essere `progetto-ml/` (D3).
- Il pre-commit esegue `ruff` e `ruff format` sui file in commit.
- Tabelle versionate prodotte da comandi diversi da quelli del README:
  `model_zoo_sc.csv` usa `fold_lopo` per task2 (README: nessun `--scheme`),
  `model_zoo_task1.csv` = `--tasks task1a_paired --tag task1`,
  `baseline_sc.csv` copre 6 task (README: 3 di default).
- `17_clonal_expansion.py` stampa rho/p senza salvarli in tabella (criterio 4).

## Fase 0b - esito dei test sul codice originale

`pytest -m "not slow"`: 79 passati, 7 falliti (log `revisione/log/fase0_pytest_fast.log`):
- test di leakage 12, 13, 20: le funzioni testabili non esistevano ancora (12 e 20 senza
  funzione per fold; 13 con normalizzazione e HVG fuori dal ciclo LOPO = leakage A1);
- `jaccard(∅, ∅)` = 0 invece di 1;
- `--help` non funziona per 21 (esegue tutto all'import), `estrai_conteggi.py` e `download_tcr.sh`.

`pytest -m slow` (log `revisione/log/fase0_pytest_slow.log`): block permutation test
dentro la banda; cluster bootstrap percentile FUORI banda: copertura 0.877, 0.885,
0.892, 0.882 (10/11 pazienti, esito continuo/binario) contro [0.90, 0.99].
Correzione dello stimatore (non della banda): intervallo "expanded percentile"
(Hesterberg 2015), alpha'/2 = Phi(sqrt(k/(k-1)) t_{0.025,k-1}); in prova copre
~0.92-0.94 con 10-11 pazienti, ma ~0.85-0.88 con 4-5 pazienti (non calibrato: B2).

Simulazioni del blocco C: un sub-agente ha cercato sulla macchina. Trovato solo (a)
in `/home/nemo/Uni/PROGETTO_ML/genomic-audit-mvp/tests/test_stats_calibration.py`
(Nadeau-Bengio sotto H0 con Dummy, potenza con logistica; cluster bootstrap con 25
gruppi). (b), (c), (d) esistono solo come testo (docstring originali, documenti
d'esame): vanno scritte da zero.

## Fase 1 - decisioni

- Stimatori statistici centralizzati in `scpipe.py` (cluster bootstrap, block
  permutation, permutazione, Holm, regola della maggioranza) e importati dagli
  script; stesse sequenze casuali di prima, cambia solo l'intervallo (espanso).
- L'intervallo espanso e' applicato a tutti i bootstrap con il paziente come unita'
  (07, 08, 09, 10, 06, 16, 17, 19); la calibrazione e' verificata solo per la media.
- A3: tutti i selettori di 12 usano lo stesso sottocampione di 20.000 cellule per
  fold (seme SEED + fold); prima L1 e RF ne usavano 8.000.
- Pareggi dell'etichetta di maggioranza: regola esistente (prima comparsa,
  identica a `value_counts().index[0]` in pandas 2.3) resa esplicita, numeri invariati.
- 07: tabelle che dipendono dalla definizione di clone ora hanno il suffisso
  `_{clone_b,clone_ab}` (prima clone_ab sovrascriveva le tabelle di clone_b).
- 13: `--splits` rimosso (non usato).

## Fase 1 - numeri cambiati finora (vecchio -> nuovo)

- 13 (A1): AUC LOPO 0.825 ± 0.079 -> 0.829 ± 0.069; duttali adiacenti predette maligne
  31.7% -> 37.1%. AUC entro paziente invariate (14 non dipende da A1).
- 06 (A5): IC dell'eccesso split-half [+0.037, +0.077] (coppie) -> [+0.017, +0.127] (pazienti, 14).
- IC del cluster bootstrap (intervallo espanso): eccesso +0.059 [+0.034, +0.100] ->
  [+0.031, +0.110]; Tumor 0.195 [0.123, 0.274] -> [0.110, 0.290]; Adjacent 0.059
  [0.000, 0.084] -> [0.000, 0.093]; Tumor - Adjacent +0.080 [+0.021, +0.170] ->
  [0.000, +0.192] (ora l'IC tocca lo zero); d1 0.723 [0.579, 0.840] -> [0.556, 0.848];
  d2 0.790 [0.682, 0.863] -> [0.660, 0.876]; LOCO Tumor - Adjacent +0.027
  [-0.002, +0.086] -> [-0.010, +0.124]. Le stime puntuali non cambiano.
- snRNA-seq: il README diceva 87.693 nuclei, i dati ne hanno 32.005 (`dataset_summary_sn`).
- `clone_switching_across.csv` versionato era prodotto con `--clone-def clone_ab`
  (908 righe); con il comando del README (clone_b) ha 1.114 righe.
- 12 (A3, stesso sottocampione di 20.000 cellule per tutti i selettori): macro-F1 medio
  task2_v1 k = 50: L1 0.936 -> 0.936, RF 0.925 -> 0.926, HVG 0.655 -> 0.659, random 0.361 -> 0.367;
  k = 200: L1 0.952 -> 0.955. task1a_paired k = 50: L1 0.729 -> 0.712, RF 0.732 -> 0.735,
  HVG 0.648 -> 0.623, random 0.587 -> 0.589. Jaccard L1 task2_v1 k = 200 0.439 -> 0.480,
  task1a k = 1000 0.275 -> 0.254; chi2 0.845 -> 0.842 / 0.674 -> 0.678; ANOVA 0.817 -> 0.794 /
  0.657 -> 0.665; HVG 0.797 -> 0.793 / 0.767 -> 0.771. Geni comuni ai 5 fold (task1a, k = 50):
  presentazione HVG 26, chi2 19, ANOVA 20 -> 18, 19, 19 (RF 16, L1 4). Nuovo dato diagnostico:
  su task2_v1 la L1 (SAGA, max_iter 400) non converge in nessun fold e ha 524-583 coefficienti
  non nulli (< 1000): a k = 1000 parte dei geni scelti ha coefficiente nullo. Riportato nel README.
  Il "0,361" e il "0,952" della presentazione sono il random k = 50 e la L1 k = 200 di task2_v1.
- 05 rieseguito: baseline invariata (0.954, 0.964, 0.875 ± 0.002, 0.729 ± 0.100).
- Simulazioni: (b) LOCO gonfia l'errore del sangue a 0.126 (moderato) / 0.385 (forte),
  non 0.24; (c) 4-11 falsi positivi su 100 (5.4% su 1.000), non 2/100; (d) differenza
  +0.143 con effetto e +0.030 senza, non 0.063 / 0.040.

## Comandi in corso

Nessuno. Log della riesecuzione di Fase 1 in `revisione/log/fase1/` (picchi di memoria in
`_status.txt`: 05 9.7 GB, 14 8.6 GB, 12 6.1 GB, 20 6.0 GB, 13 5.9 GB).
Prossimo passo: Fase 2, giro 1 (revisione completa).

## Problemi respinti

(nessuno)

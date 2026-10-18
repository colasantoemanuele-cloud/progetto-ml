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
| 1 Correzioni | fatta (A1-A7, B1-B7, C1-C4, D1-D5); pytest veloce 98 passati, ruff pulito | 52ace44 |
| 2 Loop | giro 1: 1/18/29 (8058701); giro 2: 0/4/16 (50b8095); giro 3: 0/5/8 (61fb36b); chiusura 1: 2/9/38 (59098b6); chiusura 2: 2/3/37 (04ebed5); giro 6: 2/4/22 (b696162); giro 7 (nuove regole: solo dataset reale): 2/4/15; prossimo: giro 8 (ultimo) | - |
| 3 Esecuzione finale | eseguita al commit 61fb36b (tutti rc = 0); RISULTATI_FINALI.md generato; manca il giro di chiusura | (commit Fase 3) |

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

Nessun comando di Fase 1 in corso. Log della riesecuzione di Fase 1 in `revisione/log/fase1/`
(picchi di memoria in `_status.txt`: 05 9.7 GB, 14 8.6 GB, 12 6.1 GB, 20 6.0 GB, 13 5.9 GB).

Esecuzione finale di 11 anticipata (richiesta dell'utente), codice al commit 52ace44
(`revisione/log/11_final/_commit.txt`), avviata 00:12. Un processo di 11 ha picco 5.4 GB
(VmHWM, al caricamento) e 4.7 GB a regime; RAM disponibile 17.4 GB -> 2 processi
(restano >= 8 GB liberi), 4 + 6 thread (restano 2 core liberi; OMP/OPENBLAS/MKL/LOKY limitati).
- `v1_lightgbm`: task2_v1, lightgbm, tag fin_v1_lightgbm (4 thread)
- `all_but_lgbm` -> `lgbm_hier3_task1`: i 12 modelli restanti sui 3 task (tag fin_rest), poi
  lightgbm su task2_hier3 e task1a_paired (tag fin_lgbm_rest) (6 thread)
Poi: unione (`merge11.py` nello scratchpad) in `model_zoo_{sc,task1}.csv` nell'ordine di
un'esecuzione unica, `11 --from-table` (stesso codice di 11) per riepilogo, confronti e figure,
confronto con le tabelle precedenti (copie in `revisione/log/11_final/prev_*`).
Regola della Fase 3: 11 non si riesegue se dal commit 52ace44 non cambiano il calcolo di 11,
`scpipe.py` e `04_splits.py`; lo stesso per gli altri script lunghi di Fase 1.

Ribilanciamento proposto dall'utente (elasticnet ~1.700 s, logreg_l1 ~700 s, LightGBM ~400 s):
verificato nei log locali e NON confermato per il codice attuale, quindi suddivisione invariata.
I tempi di 1.688 s e 710 s vengono da `logs/11_model_zoo.log` (10/07, 18:40): task2_hier3 con
fold_group_5 e la versione precedente dei due modelli (LogisticRegression SAGA, vedi i
FutureWarning di `_logistic.py` nel log). Dopo il passaggio a SGDClassifier (`logs/11_check.log`,
10/07 23:09: 11 s ciascuno) i tempi LOPO sono: elasticnet 39 / 26 / 12 s, logreg_l1 38 / 25 / 12 s,
LightGBM 1.450 / 400 / 72 s (task2_v1 / task2_hier3 / task1a_paired; `logs/11_model_zoo_lopo.log`,
`logs/11_model_zoo_task1.log`). LightGBM su task2_v1 resta il processo piu' lungo.
I tempi effettivi di questa esecuzione si annotano a fine corsa.

Giro 1 del loop: revisore completo in corso.

## Giro 1 (revisione completa, un revisore)

Esito del revisore: 1 ALTA, 18 MEDIA, 29 BASSA; pytest veloce 98 passati, ruff pulito.
Rapporto completo: `revisione/log/giro1_revisore.md`.
- ALTA 1: d2 dei doppietti usava la n_counts osservata delle NK-TCR+ al denominatore del
  rapporto atteso del doppietto (modello incoerente con d1): corretto con s + c (10, 09).
- MEDIA 2-19: corrette (dettaglio delle modifiche in fondo al giro).
- BASSA: corrette quelle economiche; le altre motivate in "Problemi respinti".
Correzioni divise fra: fork TCR (06-10, 17, 18, 21, scpipe, test_stats, download_tcr.sh,
01, 01b), fork CNV/sim/test (13-16, 19, sim/*, test_leakage, test_smoke) e me (04, 05, 11,
12, 20, README, envs, pre-commit).
- 11: filtro dei warning limitato al messaggio sui nomi di feature; ConvergenceWarning
  registrati per fold (colonna `converged`, riepilogo `n_folds_not_converged`) per i modelli
  iterativi; controllo di get_indexer; schema unico per task; nota sulla calibrazione di NB.
  Il calcolo delle metriche non cambia: i 5 modelli iterativi sono rieseguiti (tag fin_conv)
  per registrare la convergenza, gli altri restano quelli dell'esecuzione finale.
- 04: crosstab con colonne fisse e task degeneri saltati con messaggio. Rieseguito con tag di
  prova: `splits_sc.csv` e `task_summary_sc.csv` identici byte per byte (quindi 11 e gli
  altri script che leggono gli split non vanno rieseguiti per questa modifica).
- 12, 20: controllo di get_indexer; 20: parametro inutilizzato tolto, commento sui geni comuni
  calcolato dai dati, avvertenza sul controllo non appaiato.

Riesecuzioni del giro 1 (log in `revisione/log/giro1/`):
- TCR: 06, 07 clone_ab, 07 clone_b, 08, 09, 10, 21, 17, 18, poi 15, 16, 19: tutti rc = 0.
  Numeri cambiati: d2 0.790 [0.660; 0.876] -> 0.676 [0.559; 0.781] (ALTA 1); d2 con le medie
  0.840 -> 0.693; d2 per compartimento Tumor 0.581 -> 0.419, PBMC 0.935 -> 0.851. d1 invariato.
  18: p_full 0.0021 (10.000 permutazioni) -> 0.0009 (20.000, come 17); 10/10 invariato.
  Tolte le colonne dei Mann-Whitney per cellula (sostituite da conteggi per paziente).
- 13 senza cache (coordinate Ensembl 116 scaricate in `data/interim/gene_positions_ensembl116.csv`):
  `cnv_scores.csv` identico byte per byte alla cache dell'11/07 (spostata in
  `data/interim/cache_cnv_20260711/`); AUC 0.829 +- 0.069 e 37.1% invariati. Picco 10.8 GB.
- 14 senza cache: fermato da Claude Code per memoria critica (in parallelo con 11), nessun file
  scritto; l'utente ha chiesto di rilanciarlo da solo dopo la fine di 11.

- 14 rilanciato da solo dopo la fine di 11 (rc = 0, picco 11.4 GB): `cnv_matrix.npz` identico
  alla cache; `cnv_loci_auc.csv` cambia solo per la colonna `n_loci` (e `n_windows` del composito).
- Simulazione (c) estesa su richiesta dell'utente: 1.000 repliche per livello, IC binomiale
  esatto. FP del block test 5.0% [3.7; 6.5], 5.9% [4.5; 7.5], 4.8% [3.6; 6.3], 5.1% [3.8; 6.7]
  (confondente 0; 0.5; 1; 2): nessun livello > 8%, quindi la statistica di 17/18 non cambia.
  L'11/100 era rumore delle prime 100 repliche (stesso seme). Prima: 4/6/11/5 su 100 e 5.4% su
  1.000 (c = 1, seme 399); ora c = 1 ha seme 302 come le 100 repliche.

Esecuzione finale di 11 (codice di 11 = commit 8058701; log in `revisione/log/11_final/`):
- processi: `v1_lightgbm` (4 thread, 16:25), `all_but_lgbm` (6 thread, 46:59),
  `lgbm_hier3_task1` (5:14), `conv_iterative` (5 modelli iterativi con la convergenza, 15:06),
  `dense_rerun` (knn, svm_rbf, mlp con BLAS a 1 thread, 7:00), `mlp_check_2thr` (controllo).
  Picchi 5.2-5.8 GB per processo. Tempi per task (somma dei modelli): task2_v1 ~39 min
  (LightGBM 972 s), task2_hier3 ~18 min, task1a_paired ~8 min.
- Verifica: con lo stesso codice, 4 contro 6 thread danno fold identici per tutti i modelli
  tranne MLP (fino a 0.014 per fold): MLP dipendeva dal numero di thread BLAS. Correzione in 11:
  BLAS a un thread per i modelli densi (threadpoolctl); dopo la correzione MLP con 6 e 2 thread
  coincide (differenza 0). knn e svm_rbf cambiano di < 0.0004 per fold.
- Confronto con le tabelle precedenti: task2_v1 e task2_hier3 identiche per tutti i modelli
  salvo MLP hier3 0.9032 -> 0.9038 (README 0,903 -> 0,904); task1a_paired: differenze fino a
  0.003 (elasticnet 0.7077 -> 0.7104, LightGBM 0.7224 -> 0.7205, RF 0.7161 -> 0.7146), presenti
  anche nei modelli deterministici: la tabella precedente (3 agosto) veniva da un ambiente
  diverso. dummy 0.229 / p 0.072 e range 0,496-0,725 invariati. 0,851 del README/presentazione
  era un arrotondamento sbagliato di 0.85049 -> 0,850.
- Convergenza: logreg_l1 ed elasticnet (SGD, max_iter 50) non convergono in nessun fold; logreg_l2,
  svm_linear, mlp convergono in tutti. Scritto nel README.
- Unione: `merge11.py` (scratchpad), poi `11 --from-table` per riepilogo, confronti e figure.
- Regola della Fase 3: 11 non si riesegue se dal commit 8058701 non cambiano il calcolo di 11,
  `scpipe.py` e `04_splits.py`.
- `envs/conda-lock.yml` rigenerato (requests e threadpoolctl fissati; percorso relativo).

## Giro 2 (incrementale, `git diff 52ace44 HEAD`, un revisore)

Esito: 0 ALTA, 4 MEDIA, 16 BASSA; pytest veloce 100 passati, ruff pulito; README contro tabelle:
~230 numeri, nessuna discrepanza. Rapporto: `revisione/log/giro2_revisore.md`.
- MEDIA 1: `permutation_corr_test` e `block_permutation_test` davano un p minimo spurio con
  un input costante (rho NaN): ora restituiscono NaN con un avviso; 18 gestisce rho NaN;
  nuovo test `test_permutation_tests_constant_input_are_nan`. Per input non costanti il calcolo
  non cambia: test lenti 6/6 passati (`revisione/log/giro2_pytest_slow.log`).
- MEDIA 2: README dichiara le colonne con unita' sbagliata tenute per confronto (`lo_naive`,
  `hi_naive`, `clone_naive_permutation`) e le percentuali per cellula descrittive.
- MEDIA 3: `RISULTATI_FINALI.md` riporta B4 di 18 (p centrato, semi, 100.000 permutazioni).
- MEDIA 4: docstring del generatore corretta; pazienti, fold e cellule letti dalle tabelle.
- BASSA corrette: 6-14, 16-20 (README: fonti di 35 ricampionature e IC clone_ab, tempi
  "circa 9", test di Spearman, file non versionati; 13: cache controllata nei due sensi;
  11: docstring BLAS, `--from-table` usa tutti i task della tabella; 10: messaggio sul segno
  e verdetto solo sugli schemi calibrati; 08 e 10: niente tabelle vecchie se il risultato e'
  vuoto; sim_c: docstring; generatore: p < 0,001 in notazione esponenziale, esito
  "arrotondamento" per 0,851 -> 0,850, riga MLP, mediana 0,533 per il vecchio 0,53, nota su
  NB in LOPO). BASSA 15: commit di ogni file `fin_*` in `revisione/log/11_final/_commit.txt`.
- Nessuna riesecuzione: le modifiche non cambiano il calcolo sui dati reali (controlli su input
  degeneri, messaggi, docstring). 08, 10, 13, 18 si rieseguono comunque nella Fase 3.
- 11: dopo 8058701 cambiano solo docstring e il default di `--tasks` con `--from-table`; il
  calcolo dei modelli resta quello di 8058701. scpipe: cambiano solo le due funzioni di
  permutazione (non usate da 11, che importa LogCPM, HVGSelector, macro_f1, nadeau_bengio_ttest).

## Giro 3 (incrementale, `git diff 8058701 HEAD`, un revisore)

Esito: 0 ALTA, 5 MEDIA, 8 BASSA; pytest veloce 101 passati, ruff pulito; README contro tabelle:
~240 numeri, nessuna discrepanza; RISULTATI_FINALI.md rigenerato identico. Rapporto:
`revisione/log/giro3_revisore.md`. Tutti corretti:
- M1: `holm_bonferroni` esclude i p NaN dalla famiglia (p corretto NaN, avviso); 16 e 19
  stampano "non definito"; nuovo test `test_holm_excludes_undefined_p`.
- M2: 18 con esito costante stampa un solo messaggio; il verdetto conta solo le esclusioni con
  p definito; se il campione completo ha p >= 0,05 non da' un verdetto di robustezza.
- M3: generatore: tolto il dizionario di esiti fissi; "arrotondamento" calcolato (|esatto -
  vecchio| < 0,001); numero di modelli e di pazienti letti dalle tabelle; "IC non calibrati"
  solo con <= SMALL_K pazienti.
- M4: tolta l'attribuzione dello 0,53 alla mediana (non ricostruibile dal codice di allora):
  0,53 -> 0,540 (media) e' "cambiato".
- M5: README: "maggioranza di doppietti" invece dell'intervallo 68-72% in grassetto.
- BASSA: "1 seme"; separatori delle migliaia; README:34 (stime per schema, IC non calibrati con
  4-5 pazienti, aggiunto lo schema ingenuo +0,035); nota su Nadeau-Bengio in LOPO e sulle 1.000
  repliche del bootstrap in sim_a; 11 `--h5ad`/`--splits` non obbligatori con `--from-table`;
  docstring di sim_c, 08, 10.
- Nessuna riesecuzione necessaria (nessun calcolo cambiato sui dati reali); tutto viene
  rieseguito nella Fase 3, salvo 05, 11, 12 (vedi sotto).

## Fase 3 (esecuzione finale, prima del giro di chiusura, su richiesta dell'utente)

Comandi del README in ordine, codice 61fb36b, log `revisione/log/NN_nome.log`, stato e picchi
in `revisione/log/_status.txt`; note complete in `revisione/note_esecuzione.md` (copiate in
RISULTATI_FINALI.md). Tutti rc = 0. 02 ricrea `sc_raw.h5ad` identico byte per byte; 04 ricrea
`splits_sc.csv` identico; tabelle e figure identiche all'ultimo giro salvo l'ordine delle righe
di `ni_score_holm_joint.csv` (Holm ora restituisce l'ordine di ingresso). Non rieseguiti: 05, 11,
12 (codice di calcolo invariato, verificato con git diff e confronto AST di scpipe). Cache usate:
Scrublet, inferCNV (ricalcolata senza cache nel giro 1, identica). Nessun numero del README cambia.

## Giro 4 = chiusura 1 (completo, 4 revisori in parallelo)

Aree: (a) 04, 05, 11, 12, 20, scpipe; (b) 06-10, 15-19, 21; (c) 13, 14, sim, test; (d) README,
ambiente, docstring. Esito: (a) 1/4/10, (b) 0/3/7, (c) 0/1/9, (d) 1/1/12 -> 2 ALTA, 9 MEDIA,
38 BASSA. pytest completo (lenti compresi) 108 passati; ruff pulito. Rapporti:
`revisione/log/chiusura1_revisore_{a,b,c,d}.md`.
- ALTA (a1): `--tasks`/`--models`/`--methods` vuoti e fold con una sola classe: `nargs="+"` e
  messaggi chiari in 05, 11, 12, 20. ALTA (d1): crash del generatore se manca
  `model_zoo_summary_task1`: corretto.
- MEDIA (c1): inferCNV includeva il cromosoma mitocondriale (`chrMT` non e' escluso da
  infercnvpy) e un contig non standard; ora solo cromosomi 1-22 (X e Y esclusi da infercnvpy),
  con controllo sulle finestre. 13 e 14 rieseguiti senza cache. Numeri cambiati:
  AUC LOPO 0.829 +- 0.069 -> 0.837 +- 0.071; duttali adiacenti predette maligne 37.1% -> 30.8%;
  adiacenti "maligne" (n >= 20) 45-62% -> 41-58%; AUC pooled score medio 0.492 -> 0.531, loci PDAC
  0.682 -> 0.698; entro paziente 0.540 (0.512-0.582) -> 0.557 (0.528-0.588). Presentazione:
  0,825 -> 0,837, 31,7% -> 30,8%, 0,492 -> 0,531, 0,682 -> 0,698, 0,53 -> 0,557.
- MEDIA (a2-a5, b1-b3, d2): docstring di majority_label (l'equivalenza con value_counts vale solo
  con 2-3 etichette; nuovo test con >= 4), 05, 04; 20 senza soglie fisse (confronto con sim_d);
  README: doppietti nel Tumor (maggioranza non sostenuta), docstring di 21, tabelle vecchie
  cancellate in tutte le uscite anticipate di 08; README sezione 1 riscritta (il principio del
  paziente spiega le righe 2 e 3; riga 1: circolarita' in entrambe le stime; riga 3: pooled vs
  entro paziente, l'AUC LOPO spostata nelle note).
- BASSA: corrette quasi tutte (impronta delle cache CNV, cache Scrublet controllata, colonne
  `n_patients_across/within` in 06 e `ic_calibrato` in 09, `help` e default in tutti gli
  argparse, `\\` nelle docstring, figure nel README, ruff su `revisione/`, comando del
  generatore, test NB con k = 3, fonti delle coperture, p grezzi e selezione del migliore in 11,
  asterisco per k > geni non nulli in 12).
- Riesecuzioni: 04 (identico), 06-10, 21, 15-19, 20, 13 e 14 senza cache, simulazioni, 11 con
  `--from-table`; tutti rc = 0 (`revisione/log/_status.txt`, blocco "chiusura 1").
- Respinti: (a11) LogCPM/HVGSelector con matrici dense: gli script passano sempre CSR, e
  cambiare queste classi invaliderebbe la verifica (AST identico) che permette di non rieseguire
  05, 11 e 12; (a6) rinominare `model_comparison_{task}` e le figure con il tag: cambierebbe i
  nomi pubblicati; dichiarato in docstring che dipendono solo dal task.

## Giro 5 = chiusura 2 (completo, 4 revisori in parallelo)

Esito: (a) 0/0/11, (b) 1/0/11, (c) 0/2/8, (d) 1/1/7 -> 2 ALTA, 3 MEDIA, 37 BASSA. pytest completo
108 passati (`revisione/log/chiusura2_pytest_all.log`); ruff pulito; tutti i numeri del README e di
RISULTATI_FINALI coincidono con le tabelle. Rapporti: `revisione/log/chiusura2_revisore_{a,b,c,d}.md`.
- ALTA (b1): 19 con un foglio Table S1A vuoto: messaggio chiaro. ALTA (d1): generatore ancora
  fragile con tabelle o colonne mancanti: ora tutta la costruzione e' protetta e ogni errore e'
  registrato fra i valori mancanti (provato togliendo tabelle e colonne).
- MEDIA (c1 = d2): docstring di 13 con il vecchio 45-62%: tolto il numero. MEDIA (c2): README,
  i 4 pazienti entro paziente sono quelli con >= 10 duttali per tessuto.
- BASSA: corrette tutte salvo quelle gia' respinte (uscite anticipate che cancellano le tabelle
  vecchie in 06, 07, 09, 13, 16, 19; Scrublet con controllo esatto anche con --skip-scrublet; IC
  non calibrati segnati nelle figure di 07 e 09; rho NaN in 16, 17, 19; `choices` per --scheme e
  --methods; avvisi non di convergenza riemessi in 12; --update di 11 non cancella i modelli
  falliti; nota sul p minimo di Wilcoxon con 5 fold; sim_a --quick su file separati; README:
  sim_d prima di 20, "quasi lineare" tolto, L1 non convergente su task2_v1, deriva di H1,
  2 su 100 con il p centrato, 0,24 non riprodotto; verifica parziale della cache Scrublet
  documentata).
- Riesecuzioni: 06-10, 21, 15-19, simulazioni, 20, 13 e 14 (cache con impronta), 11 con
  --from-table (blocco "chiusura 2" di `revisione/log/_status.txt`).

## Giro 6 = chiusura 3 (completo, 4 revisori)

Esito: (a) 0/2/3, (b) 1/0/7, (c) 0/2/9, (d) 1/0/3 -> 2 ALTA, 4 MEDIA, 22 BASSA; pytest completo
108 passati (`revisione/log/chiusura3_pytest_all.log`); ruff pulito; numeri del README e di
RISULTATI_FINALI coincidenti. ALTA + MEDIA: 11 (giro 4) -> 5 (giro 5) -> 6 (giro 6): primo giro
senza diminuzione. Elenco in `revisione/STATO.md`. L'utente ha portato il limite a 10 giri, con
stop se ALTA + MEDIA non diminuisce per due giri di fila o se un problema corretto ricompare
(correzioni che si annullano). Nessuno dei due casi: i due ALTA collegati a giri precedenti sono
estensioni (NaN dopo l'input costante; un calcolo rimasto fuori dal blocco protetto), non
correzioni annullate.

Correzioni del giro 6 (fork TCR, fork CNV/sim/test e coordinatore):
- ALTA: `permutation_corr_test` e `block_permutation_test` restituiscono NaN anche con valori
  mancanti o non finiti (nuovo test); 19 esclude i pazienti senza NI score stampandone gli ID e
  nomina la colonna di Holm con il numero di p definiti. Generatore: `auc_pat` protetto.
- MEDIA: 20 accetta solo `fold_lopo` e controlla che ogni fold escluda un solo paziente; README
  riga 2: "struttura per campione" (libreria = paziente x tessuto, meccanismo non misurato) al
  posto della firma del paziente; 14: confronto con lo score medio sulle AUC entro paziente;
  13 e 14: riferimento immunitario vuoto -> messaggio chiaro.
- BASSA: corrette tutte (figure di 10 e 21 con IC non calibrati segnati; uscite anticipate di
  07, 09, 10, 16; docstring di 06, 07, 09, 10, 13, sim_a, test; impronta CNV con versione e
  parametri di infercnvpy; 04 con la classe piu' numerosa < 5; 20 con intestazioni dal numero di
  fold; default in --help solo se esistono; README: TRB produttiva, schema ingenuo, IC di 06,
  score medio entro paziente 0,572, fonti di Δ = −0,229 e delle 11 cellule).
- Riesecuzioni: 04, 06-10, 21, 15-19, simulazioni, 20, 13 e 14 senza cache (impronta nuova).

## Giro 7 = chiusura 4 (completo, 4 revisori, regole dell'utente: solo dataset reale)

Esito: (a) 0/0/6, (b) 0/2/5, (c) 2/0/2, (d) 0/2/2 -> 2 ALTA, 4 MEDIA, 15 BASSA, 3 limiti noti
nuovi; pytest completo 109 passati. ALTA + MEDIA: 5 (giro 5), 6 (giro 6), 6 (giro 7): condizione di
stop dell'utente soddisfatta (due giri senza diminuzione; criteri cambiati al giro 7). L'utente ha
deciso le correzioni e aperto il giro 8 come ultimo.
- Pareggi dell'etichetta di maggioranza (decisione dell'utente): `majority_label` restituisce None
  in caso di pareggio; le unita' in pareggio sono escluse dal calcolo principale in 06, 09, 10 e nei
  consumatori 16-19, 21; sensibilita' "pareggi = errore" / "pareggi = corretti" in
  `error_rate_tie_sensitivity` e `arbitration_tie_sensitivity`. Numeri (schema solo sangue):
  Tumor 0.195 [0.110; 0.290] -> 0.170 [0.095; 0.252] (277 cloni, 10 pareggi; 0.199 e 0.164 nelle
  sensibilita'); Adjacent_normal 0.059 [0.000; 0.093] -> 0.060 [0.000; 0.095] (117, 1 pareggio);
  differenza +0.080 [0.000; 0.192] -> +0.059 [0.000; 0.160] (85 cloni). Scomposizione 38/18 ->
  32/15; senza NK 0.141 -> 0.122. LOCO +0.027 -> +0.016, naive +0.035 -> +0.024; regole CD8
  0.206 -> 0.180/0.181; split-half di 06 +0.057 -> +0.063; espansione clonale (errore) rho 0.177
  -> 0.220, p 0.0009 -> 0.0004; H1-H4 rho -0.45/-0.50/-0.13/-0.13, p_Holm 0.58/0.58/1/1.
- Soglia di inferCNV (decisione dell'utente): filtro del rumore con una sola soglia globale
  (1,5 sd di tutta la matrice) invece che per blocco di 5.000 cellule. Finestre CNV: coordinate
  reali dei geni (`cnv_windows.csv`); CDKN2A in 10 finestre, tutte attribuite a 9p21.
  Numeri: AUC LOPO 0.837 +- 0.071 -> 0.836 +- 0.062; adiacenti predette maligne 30.8% -> 37.4%;
  adiacenti "maligne" 41-58% -> 38-65%; AUC pooled score medio 0.531 -> 0.598, loci 0.698 ->
  0.717; entro paziente loci 0.557 -> 0.566, score medio 0.572 -> 0.616.
- MEDIA (b1 = d1): eccezioni al lignaggio stabile dichiarate (README, 09, 10) e verificate
  (`clone_invariant_tcr_check`: 12 MAIT-like tutti concordanti, 0 iNKT). MEDIA (d2): README sulle
  duttali adiacenti "maligne": i dati non distinguono artefatto, cellule tumorali e campo.
- BASSA corrette; nuovi limiti noti (TISSUE_MAP, top-50 di 20) in `revisione/limiti_noti.md`.
- Riesecuzioni: 05 (solo legenda; tabelle identiche), 06-10, 21, 15-19, simulazioni, 20, 13 e 14
  senza cache; tutti rc = 0.

## Problemi respinti

- Giro 1, n. 13 (RISULTATI_FINALI.md assente): il file e' un prodotto della Fase 3; la bozza e'
  ora generata da `revisione/genera_risultati.py` (solo dalle tabelle). B6 aggiunto al README;
  la parte su 18 (B4) e' accettata e corretta.
- Giro 1, n. 43 (unione dei CSV di 11 non e' un comando del README): l'esecuzione in processi
  paralleli e' solo un'ottimizzazione di tempo; le tabelle unite si confrontano con quelle di
  un'esecuzione unica (stessi semi, stessi fold) e il comando del README le riproduce.
  Il "circa un'ora per task" e' sostituito dai tempi misurati.
- Giro 1, n. 48 (cache inferCNV dell'11 luglio): accettato; 13 rieseguito senza cache nel giro 1.
- Giro 1, n. 39 per 05 (`n_classes` non usato in `make_pipeline`): respinto. Il parametro non
  cambia il calcolo; toglierlo cambierebbe il codice di 05 (17 minuti, 9.7 GB) e
  obbligherebbe a rieseguirlo in Fase 3 senza alcun effetto sui numeri. Corretto in 20.
- Giro 2, n. 5 (valori della versione precedente nel README: d2 0,790, 0,24, 2 su 100,
  0,063/0,040): respinto. Le regole della revisione chiedono di riportare il vecchio e il nuovo
  valore di ogni numero pubblicato che cambia; nel README sono dichiarati come valori della
  versione precedente, non come risultati.

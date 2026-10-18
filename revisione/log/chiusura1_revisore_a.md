# Giro di chiusura 1 - area (a) ML

**1. `src/05_baseline.py:296-298`, `src/11_model_zoo.py:453-459`, `src/12_feature_selection.py:290-302`; anche `05:182`, `12:212`, `20:136`**
- gravità: ALTA (crash su input vuoti o degeneri; i risultati pubblicati non cambiano)
- criterio violato: 8
- motivo: con `--tasks` passato senza valori, `args.tasks` vale `[]` e `pd.concat([])` fallisce con `ValueError: No objects to concatenate`. Succede in 05 e in 11, dopo il caricamento completo dell'h5ad.
  - In 12, `--methods` senza valori lascia `rows` vuoto. Lo script scrive allora `feature_selection_{task}.csv` e `feature_selection_diag_{task}.csv` vuoti sopra le tabelle pubblicate, e poi va in crash con `KeyError` in `res.pivot_table(index="method")`.
  - Se il training di un fold ha una sola classe, 05, 12 e 20 si fermano con l'errore di sklearn (`This solver needs samples of at least 2 classes...`, riprodotto chiamando `run_task` / `run_fold` su dati sintetici). Il messaggio non dice quale task o fold è il problema. 11 invece intercetta l'errore (FALLITO).
- correzione proposta: `nargs="+"` per `--tasks` e `--methods` (argparse dà un errore chiaro). In 05, 12 e 20 un controllo `len(np.unique(y[tr])) < 2` → `SystemExit` con task e fold.

**2. `src/scpipe.py:385-389` (docstring di `majority_label`)**
- gravità: MEDIA
- criterio violato: 6
- motivo: il docstring dice «E' la stessa regola di `pd.Series.value_counts().index[0]` in pandas >= 2». È falso con pandas 2.3.3 (quello dell'ambiente).
  - Ho provato 5.000 liste casuali per ogni numero di etichette distinte: 0 discrepanze con 2-3 etichette, da 33 a 222 con 4-14 etichette.
  - Esempio: `[Endothelial, CD8T, Ductal×3, NK×3, B, CD4T, Endothelial, B, CD8T, CD4T]`. `majority_label` dà Ductal (pareggio, compare prima), `value_counts` dà NK.
  - Il test `tests/test_stats.py:238-242` usa solo 3 etichette, quindi non lo vede.
- correzione proposta: limitare l'equivalenza a ≤ 3 etichette distinte, oppure togliere la frase. Aggiungere al test un caso con ≥ 4 etichette.

**3. `src/05_baseline.py:9-11`**
- gravità: MEDIA
- criterio violato: 6
- motivo: la docstring dice «`LogCPM`, `HVGSelector` ... e `StandardScaler` senza centratura, tutti da `scpipe.py`». `StandardScaler` invece è importato da `sklearn.preprocessing` (riga 72).
- correzione proposta: «`LogCPM` e `HVGSelector` da `scpipe.py`, `StandardScaler(with_mean=False)` di scikit-learn».

**4. `src/04_splits.py:96-97`**
- gravità: MEDIA
- criterio violato: 6
- motivo: il commento dice «Task a supporto sparso ... (valutazione LOPO)» per `EXPECT_INCOMPLETE = {"rare", "task1a_all"}`. Ma task1a_all non è in `LOPO_PRIMARY` (riga 98), e `task_summary_sc.csv` riporta `task1a_all,...,group_5fold`.
- correzione proposta: «task per cui la copertura di tutte le classi in ogni fold di fold_group_5 non è garantita (rare è valutato in LOPO, task1a_all in 5 fold per paziente)».

**5. `src/20_xai_per_patient.py:312-356`**
- gravità: MEDIA (la lettura stampata per il dato reale, +0,159, resta sostenuta)
- criterio violato: 5
- motivo: le soglie 0,05 e 0,02 contraddicono la simulazione dello script stesso. In `sim_d_xai_stability.csv`, senza effetto paziente, la differenza media è +0,030 e arriva a +0,052 (rep 1).
  - Con un valore nullo come 0,052 il ramo > 0,05 stamperebbe «La variabilità osservata in LOPO eccede quella di puro campionamento».
  - Il ramo ≤ 0,02 afferma «la variabilità è rumore, non effetto-paziente»: è un'assenza di evidenza presentata come conclusione, e la simulazione dice che ≤ 0,02 è sotto la media nulla.
- correzione proposta: niente soglie fisse. Stampare la differenza accanto all'intervallo di `sim_d` per i due scenari, oppure tarare le soglie su quell'intervallo. Rendere condizionale il ramo «else».

**6. `src/11_model_zoo.py:382, 514-517`; `src/12_feature_selection.py:291-367`; `src/20_xai_per_patient.py:255-274`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: `model_comparison_{task}.csv` e `model_zoo_{task}.png` non dipendono da `--tag`. Seguendo il messaggio «usa un --tag diverso per ogni schema» (riga 471), un'esecuzione con `--scheme fold_random_5 --tag rnd` sovrascrive i file pubblicati.
  - Inoltre l'asse x della figura dice sempre «split per paziente», anche con `fold_random_5`.
  - In 12 e 20 i file sono indicizzati solo per task: un'esecuzione con altro `--scheme`, o con un sottoinsieme di `--methods`, sovrascrive le tabelle del README.
- correzione proposta: mettere il tag (o lo schema) nei nomi di questi file e un'etichetta dell'asse che dipenda dallo schema.

**7. `src/20_xai_per_patient.py:8-11` e `README.md:27`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: dicono che il controllo «toglie lo stesso numero di cellule a caso». In realtà entrambi sono sottocampioni di 20.000 cellule: quelli LOPO vengono da circa 54.000 cellule di 4 pazienti, quelli di controllo da tutte le 67.662 cellule dei 5 pazienti. È uguale la dimensione del training, non il numero di cellule tolte.
- correzione proposta: «sottoinsiemi casuali di 20.000 cellule, la stessa dimensione dei training LOPO sottocampionati».

**8. `README.md:25` (task2_hier3)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: «la logistica L2 ottiene 0,850 e i modelli non lineari di più» suggerisce che il guadagno venga dalla non linearità. Ma in `model_zoo_summary_sc` L2 è il modello più basso escluso il dummy (decision_tree 0,853). Anche modelli lineari la superano: logreg_l1 0,877 (+2,7), svm_linear 0,865 (+1,5), elasticnet 0,862.
- correzione proposta: dire che L2 è il riferimento più debole su hier3 e che anche i lineari L1 e SVM lineare lo superano.

**9. `README.md:25` (task1a_paired)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: «escluso il dummy, nessuna differenza ... è significativa (p ≥ 0,080; dummy p = 0,072)» si legge come se il dummy fosse significativo. Invece neanche il dummy (Δ = −0,229, p = 0,072 in `model_comparison_task1a_paired`) supera 0,05, ed è l'indizio più forte della bassa potenza.
- correzione proposta: «nessuna differenza dal migliore è significativa, nemmeno quella del dummy (Δ −0,229, p = 0,072); escluso il dummy p ≥ 0,080».

**10. `src/scpipe.py:184-185, 12-43, 232`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: tre imprecisioni nelle docstring.
  - `expanded_alpha` dice che il percentile semplice copre «~87-89%» con 10-11 pazienti; `sim_a_cluster_bootstrap.csv` dà 0,8875-0,9225 (README: 88,8-92,3%).
  - La sezione «Contenuto» del modulo non elenca `majority_label`, `expanded_alpha`, `percentile_ci` e `cluster_bootstrap_replicates`.
  - Il docstring di `cluster_bootstrap` dice «NaN se ci sono meno di 2 pazienti», ma con 1 paziente restituisce la media e solo l'IC è NaN.
- correzione proposta: allineare i numeri alla tabella e completare l'elenco e la descrizione del valore restituito.

**11. `src/scpipe.py:79-112`**
- gravità: BASSA
- criterio violato: nessuno (gli script passano sempre CSR)
- motivo: `LogCPM.transform` e `HVGSelector.fit` vanno in crash con un `ndarray` denso (`AttributeError: 'numpy.ndarray' object has no attribute 'multiply'`, verificato). `HVGSelector`, se nessun gene supera `min_mean`, seleziona 0 colonne senza avviso.
- correzione proposta: `sp.csr_matrix(X)` all'ingresso e un errore esplicito se k = 0.

**12. `src/12_feature_selection.py:188-201`; `RISULTATI_FINALI.md` (tabella task2_v1, riga l1_logreg)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: a k = 1000 su task2_v1 tra 417 e 476 geni scelti da L1 hanno punteggio zero (n_nonzero 524-583), e il loro ordine nei pareggi di `argsort` è arbitrario.
  - F1 0,953 e Jaccard 0,568 a k = 1000 sono riportati senza un avviso nella riga; c'è solo la nota sui coefficienti non nulli.
  - Inoltre il docstring di `feature_selection_summary_{task}` («media e sd») non elenca `balanced_accuracy_mean`, `n_folds` e `any_k_gt_nonzero`, e il parametro `classes` di `run_fold` non è usato.
- correzione proposta: segnare quei valori come non interpretabili (o NaN) quando `k_gt_nonzero` è vero, e completare la docstring.

**13. `src/11_model_zoo.py:344-382`; `RISULTATI_FINALI.md` («8 modelli su 12 con p di Nadeau-Bengio < 0,05»)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: il «migliore» è scelto sugli stessi fold su cui poi viene testato contro gli altri 12 modelli, senza correzione per confronti multipli. Il conteggio di p grezzi < 0,05 risente quindi della selezione.
- correzione proposta: dichiarare che i p sono grezzi e condizionati alla scelta del migliore, oppure riportare anche p di Holm.

**14. `src/20_xai_per_patient.py:141-147, 246`; `README.md:27`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: lo Spearman (0,234) è calcolato sui vettori estesi a tutti i 30.543 geni di `var`. I geni esclusi dal prefiltro in tutti i fold sono zeri a pari merito in ogni vettore, quindi il valore dipende dall'insieme dei geni considerati. Il README non lo specifica.
- correzione proposta: calcolare lo Spearman sull'unione dei geni passati dal prefiltro (o sui soli coefficienti non nulli), oppure dichiararne la definizione.

**15. `src/04_splits.py:173-188, 212`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: un task con meno di 5 cellule manda in crash `StratifiedKFold(n_splits=5)` (meno di 2 pazienti o di 2 classi sono già gestiti). Il controllo anti-leakage è un `assert`, che non viene eseguito con `python -O`.
- correzione proposta: saltare i task con n < N_FOLDS dando un messaggio, e sostituire l'`assert` con `raise SystemExit`.

**Verifiche senza problemi:**
- Tutti i numeri delle mie sezioni del README e di RISULTATI_FINALI (Leakage, 13 modelli, Selezione delle feature, Spiegazioni, righe delle fonti) coincidono con `results/tables`. Controllati tra gli altri: 0,954/0,964/0,875±0,002/0,729±0,100/Δ −0,146/50×; 0,943/0,962/0,850/0,923/+7,2/+5,6/+5,3/+5,1; 0,496-0,725, p ≥ 0,080, dummy 0,072; i valori della selezione delle feature, 524-583 e 177-268; 0,300/0,459/+0,159/0,234, 0,220-0,408, 9 geni; sim_d +0,143/+0,030.
- L'unione dei file `fin_*` in `revisione/log/11_final/` ricostruisce esattamente `model_zoo_sc.csv` (338 righe) e `model_zoo_task1.csv` (65 righe): differenze 0, convergenza identica. Le affermazioni sui thread nella docstring di 11 sono confermate dai log: MLP fino a 0,0140, knn e svm_rbf sotto 0,0004.
- Con `git diff`: 05 è invariato da 52ace44; 12 da 52ace44 cambia solo per il controllo dei `cell_id`; 11 da 8058701 cambia solo in docstring e `--from-table`; in scpipe cambiano solo le funzioni non usate da 05/11/12; in 04 cambiano solo i controlli.
- Correzioni A3, A6 (scpipe, 05, 20, 11:242-264), A7 (20), C(d), D2 (scpipe/05 LogCPM, 11:21, 11:207, 11:316-317, 20:279, geni di 20, 04:27; le Endocrine sono in 12 pazienti, 2 con ≥ 10 cellule, verificato su `obs`) e D3 per la mia area: tutte applicate.
- Leakage e unità statistica: nessun problema nei sei file. I pacchetti importati sono tutti in `envs/environment.yml`.

**Conteggio:** ALTA 1 / MEDIA 4 / BASSA 10

**pytest** (`-m "not slow"`): 102 passati, 0 falliti, 0 saltati (6 lenti esclusi dalla selezione).
**ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni.

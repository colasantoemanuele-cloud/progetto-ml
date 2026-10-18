# Giro di chiusura 2 - area (a) ML

**1. `src/11_model_zoo.py:370-373` (anche `26-29`); `results/tables/model_comparison_task1a_paired.csv`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: con 5 fold LOPO il p esatto bilaterale di Wilcoxon non può scendere sotto 2/2^5 = 0,0625. Su task1a_paired il test non può quindi mai rifiutare. In tabella dummy, decision_tree e naive_bayes hanno `p_wilcoxon` = 0,0625, cioè il minimo possibile, e il valore si legge come "quasi significativo". La docstring parla solo di "potenza bassa", ma qui il rifiuto è impossibile per costruzione.
- correzione proposta: scrivere nella docstring (e nella NOTA stampata a riga 390) che con k fold il p minimo di Wilcoxon è 2/2^k (0,0625 con k = 5). In alternativa lasciare `p_wilcoxon` a NaN quando k < 6.

**2. `src/12_feature_selection.py:193, 261-266, 288`; `src/11_model_zoo.py:298, 413-415`; `src/20_xai_per_patient.py:205-207, 228`**
- gravità: BASSA
- criterio violato: 8 (marginale: il traceback nomina la chiave)
- motivo: `--methods` di 12 e `--scheme` di 11, 12 e 20 non sono validati.
  - Un nome sbagliato dà `KeyError` grezzo: in 12 in `SELECTORS[meth]` e in `sub[args.scheme]`, in 11 in `sub[scheme]`, in 20 in `sub[args.scheme]`.
  - Il crash arriva dopo il caricamento completo dell'h5ad (minuti, GB di memoria).
  - 11 invece valida `--models` (righe 432-434), quindi il comportamento non è uniforme.
- correzione proposta: `choices=list(SELECTORS)` per `--methods` e `choices=["fold_group_5", "fold_lopo", "fold_random_5"]` per `--scheme`, oppure un controllo prima di `read_h5ad`.

**3. `src/12_feature_selection.py:139-143`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: durante il fit di SAGA il blocco `catch_warnings(record=True)` registra tutti gli avvisi, ma lo script legge solo i `ConvergenceWarning`. Gli altri (FutureWarning, RuntimeWarning e simili) vengono scartati in silenzio. 11 (righe 278-282) li riemette con `warnings.showwarning`. La correzione A3 chiedeva di non sopprimere gli avvisi.
- correzione proposta: riemettere gli avvisi non di convergenza come in `fit_counting_convergence` di 11, o riusare quella funzione.

**4. `src/20_xai_per_patient.py:343-348`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: il messaggio stampato (e presente in `revisione/log/20_xai_per_patient.log`) dice «togliere un paziente cambia la regola appresa piu' che togliere cellule a caso». La docstring (righe 11-13) e il README (riga 28) sono stati corretti nel giro precedente: «e' uguale la dimensione del training, non il numero di cellule tolte»; «addestrare su un sottoinsieme casuale della stessa dimensione». Il messaggio è rimasto alla formulazione vecchia.
- correzione proposta: «... piu' che addestrare su un sottoinsieme casuale di cellule della stessa dimensione».

**5. `src/20_xai_per_patient.py:41-42, 328-370`; `README.md:131` e `README.md:161`**
- gravità: BASSA
- criterio violato: nessuno (nessuna tabella cambia)
- motivo: 20 legge `sim_d_xai_stability.csv` per la "LETTURA" finale. Nel README però 20 si esegue nel blocco 3 (riga 131) e `sim_d` solo nel blocco 7 (riga 161). Eseguendo i comandi in ordine da zero, 20 stampa sempre il ramo «assente». Il log della Fase 3 mostra il confronto solo perché `sim_d` esisteva già da un'esecuzione precedente.
- correzione proposta: spostare 20 dopo `sim_d` nel README, oppure dichiarare nella nota dei comandi che la lettura di 20 richiede `sim_d` eseguito prima.

**6. `src/20_xai_per_patient.py:338-363`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: se `sim_d_xai_stability.csv` esiste ma manca la condizione `no_patient_effect`, `lo0` e `hi0` valgono NaN. Allora `diff > hi0` e `diff >= lo0` sono entrambi falsi e lo script stampa «la differenza e' sotto le repliche senza effetto paziente», che è falso.
- correzione proposta: se `no_patient_effect` manca o è NaN, stampare che il confronto non è disponibile.

**7. `src/11_model_zoo.py:45-48` (sezione Input)**
- gravità: BASSA
- criterio violato: nessuno (docstring incompleta)
- motivo: con `--from-table` e con `--update` lo script legge `results/tables/model_zoo_{tag}.csv` (righe 439-444 e 472-473). La sezione Input elenca solo `--h5ad` e `--splits`.
- correzione proposta: aggiungere `results/tables/model_zoo_{tag}.csv` (input con `--from-table` e `--update`).

**8. `src/11_model_zoo.py:334-340, 472-476`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: con `--update`, se un modello rieseguito FALLISCE, le sue righe vecchie vengono comunque tolte (`drop = old.task.isin(args.tasks) & old.model.isin(args.models)`) e non sostituite. Il modello sparisce dalla tabella; il messaggio dice solo «sostituite N righe».
- correzione proposta: calcolare `drop` solo sulle coppie (task, modello) presenti in `res`.

**9. `revisione/RISULTATI_FINALI.md:168` (generata da `revisione/genera_risultati.py:632-636`)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: la nota «\* k maggiore dei geni con punteggio > 0 ...» si stampa sotto ogni tabella. Nella tabella di task1a_paired non c'è nessun asterisco: in `feature_selection_summary_task1a_paired.csv` `any_k_gt_nonzero` è sempre False, perché L1 ha 1616-1984 coefficienti non nulli.
- correzione proposta: stampare la nota solo se almeno una cella della tabella è marcata.

**10. `README.md:27` (stabilità di L1 su task2_v1)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: «La stabilità ... è più bassa per L1 (0,480 a k = 200 su task2_v1 ...)». Su task2_v1 SAGA non converge in nessun fold (0/5, 400 iterazioni, `feature_selection_diag_task2_v1`), quindi parte dell'instabilità può dipendere dall'ottimizzazione interrotta. Il README cita la non convergenza nella frase precedente, ma non come possibile causa del Jaccard basso. La conclusione resta sostenuta da task1a_paired, dove L1 converge (177-268 iterazioni) e il Jaccard è 0,254.
- correzione proposta: aggiungere «su task2_v1 il valore viene da fit non convergenti; su task1a_paired, dove L1 converge, il Jaccard resta basso (0,254)».

**11. `README.md:26` («Su `task2_v1` la separazione e' quasi lineare»)**
- gravità: BASSA
- criterio violato: nessuno
- motivo: è un'interpretazione non misurata. In `model_comparison_task2_v1` tutti i modelli lineari sono significativamente sotto LightGBM al t di Nadeau-Bengio: logreg_l2 p = 0,025, logreg_l1 0,041, elasticnet 0,005, svm_linear 0,0007.
- correzione proposta: «i modelli lineari arrivano a 0,936-0,943 contro 0,962 di LightGBM», senza «quasi lineare».

**Verifiche senza problemi:**
- **Numeri.** Tutti i numeri del README (sezione 1, «Confronto dei modelli», righe delle fonti) e di RISULTATI_FINALI (Leakage, 13 modelli, Selezione delle feature, Spiegazioni) coincidono con `results/tables`. Ne ho controllati 60 circa, tra cui:
  - 0,954 / 0,964 / 0,875 ± 0,002 / 0,729 ± 0,100, Δ −0,009 / −0,146, 50,3×;
  - 0,943 / 0,962 / 0,850 / 0,923 e +7,2 / +5,6 / +5,3 / +5,1 / +2,7 / +1,5;
  - 0,496-0,725, p ≥ 0,080 (decision_tree 0,0802), dummy 0,072; 8/12, 12/12, 0/12;
  - selezione delle feature: tutti i valori, 524-583, 177-268, 1616-1984, geni comuni 18/19/19/16/4;
  - XAI: 0,300 / 0,459 / +0,159 / 0,234, 0,220-0,408, 9 geni; sim_d +0,143 (+0,116; +0,181) e +0,030 (−0,013; +0,052).
- **Esecuzione finale di 11.** L'unione dei file `model_zoo_fin_*` secondo `_commit.txt` ricostruisce esattamente `model_zoo_sc.csv` (338 righe) e `model_zoo_task1.csv` (65 righe): differenza massima 0 in macro_f1, BA, MCC, n_train e n_test, `converged` identica, un solo schema (fold_lopo).
- **File non rieseguiti in Fase 3.** Confronto dell'AST senza docstring:
  - in `scpipe.py` da 52ace44 cambiano solo `SMALL_K`, `small_k_diagnostics`, `permutation_corr_test`, `block_permutation_test` e `holm_bonferroni`, nessuna usata da 05, 11 e 12;
  - 05 e 12 da 52ace44 cambiano solo in `main`/`run_task`, con i controlli sugli input;
  - 11 da 8058701 cambia solo in `main`;
  - 04 cambia solo nei controlli.
  - `baseline_*_sc.csv` è identico alla versione 4d42a8f. In `baseline_sc.csv` non c'è nessun F1 NaN, quindi il passaggio di macro_f1 alle sole classi presenti non poteva cambiare i numeri. In `results/` non ci sono file orfani per 05, 11, 12 e 20.
- **Correzioni dell'elenco per la mia area: tutte applicate.** A3, A6 (scpipe, 05, 20, 11:242-264), A7 (20), C(d), D2 (LogCPM, 11:21, 11:207, 11:316-317, 20:279, geni di 20, 04:27), D3 (Jaccard su fold_group_5, sottocampionamento non stratificato, +7 solo LightGBM, intervallo di task1a, comandi 11 e 12).
- **Leakage e unità statistica.** Nessun problema: LogCPM, HVG, SVD, scaler, prefiltro dei geni e selettori sono stimati sul solo training del fold; i confronti sono fra fold per paziente. Le penalità L1 con `l1_ratio=1.0` in sklearn 1.9.0 sono davvero L1 (verificato su dati sintetici con saga e liblinear). Tutti i pacchetti importati sono in `envs/environment.yml`.

**Conteggio:** ALTA 0 / MEDIA 0 / BASSA 11

**pytest** (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 lenti esclusi dalla selezione).
**ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni.

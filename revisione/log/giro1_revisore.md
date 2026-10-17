# Giro 1 - rapporto del revisore

Ho trovato 1 problema ALTA, 18 MEDIA e 29 BASSA. Il problema ALTA è la stima d2 dei doppietti: il codice dà 0,790 invece di circa 0,676. `pytest` passa e `ruff` non segnala nulla.

Esiti:
- `pytest -m "not slow" -p no:cacheprovider`: 98 passati, 0 falliti, 0 saltati, 6 deselezionati (i test lenti, 52,7 s). Il log `revisione/log/fase1/pytest_slow.log` riporta i 6 test lenti passati.
- `ruff check src tests`: "All checks passed!".

Il controllo 2 (numeri del README contro le tabelle) ha verificato circa 150 numeri; quelli che non tornano sono i problemi 14 e 15.

## Correzioni verificate come applicate
- **A1-A4:** applicate; il test di leakage su 13 esiste.
- **A5:** il bootstrap ora ricampiona i 14 pazienti; IC +0,017; +0,127.
- **A6:** tutte le voci. `baseline_sc.csv` è identico a 12a356d. Nessun fold di test manca di una classe, quindi il nuovo `macro_f1` non cambia nessun numero, nemmeno quelli di 11.
- **A7:** tabelle presenti per 08, 10, 14, 19 e 20.
- **B1-B3, B5, B6:** applicate. **B4:** solo in 17, non in 18 (problema 13).
- **C:** sim_a, sim_b, sim_c e sim_d esistono e girano; il README riporta i numeri ottenuti al posto di quelli citati (0,24; 2 su 100; 0,063 e 0,040).
- **D1, D2:** tutte le voci elencate.
- **D3:** applicato, salvo i problemi indicati sotto.

## ALTA

**1. `src/10_loco.py:475` (stessa formula in `src/09_arbitrate.py:327`)**
- Criterio: 6 (anche 5).
- Motivo: la docstring di 09 (righe 19-20) dice che in un doppietto T+NK gli UMI attesi sono n_counts(NK TCR-) + n_counts(CD8T), cioè s + c. Il codice invece usa la mediana osservata della miscela: `r_dbl = r_cd8 * c / max(o, 1)`. Le due formule coincidono solo se tutte le NK-TCR+ sono doppietti; con d < 1 d2 risulta gonfiato.
- Evidenza (da `doublet_physical.csv`: r_cd8 = 4,0087; r_nk = 2,6101; c = 3277; s = 3500; o = 5868):
  - d2 del codice = 1,3986 / 1,770 = **0,790**;
  - d2 con s + c = 1,3986 / 2,070 = **0,676**.
- Conseguenze:
  - README:35 dice «72-79% doppietti»;
  - 06:22-24 dice che 10 «mostra» che sono soprattutto doppietti;
  - c'è anche un'incoerenza di modello: d1 tratta le non-doppietto come singoletti NK, d2 dà loro il rapporto TCR di una T.
- La formula è uguale in 12a356d: non è stata introdotta dalle correzioni.
- Correzione: usare s + c nel denominatore (in 10 e nel valore atteso stampato da 09), oppure un solo modello di miscela coerente per d1 e d2. Poi rieseguire 10, aggiornare `doublet_mixture*.csv` e il README, e attenuare «mostra» in 06.

## MEDIA

**2. `src/scpipe.py:24-29`**
- Criterio: 3.
- Motivo: la docstring mette `permutation_corr_test` fra gli stimatori «validati in tests/test_stats.py», ma nessun test lo usa (grep: solo 16, 17, 19 e sim_c) e non c'è una banda dichiarata. Inoltre sta nella sezione «paziente come unità», mentre la sua docstring dice «unità = righe». Lo usano H1-H4 e il livello paziente di 17.
- Correzione: aggiungere un test sotto H0 con n = 10 e una banda nel docstring (es. frazione di p < 0,05 in [0,02; 0,08]), oppure toglierlo dagli stimatori validati.

**3. `src/08_disambiguate.py:194-207`**
- Criterio: 5.
- Motivo: il contrasto sui cloni presenti nei tre compartimenti ha solo 4 pazienti, ma il codice stampa «IC esclude 0» per +0,066 [+0,012; +0,126].
  - L'IC coincide con il minimo e il massimo delle 4 differenze per paziente.
  - Con 4 pazienti l'IC copre l'82-83% (`sim_a_cluster_bootstrap.csv`).
  - Manca la diagnostica `small_k_diagnostics` che 10 invece riporta.
- Correzione: aggiungere la diagnostica e, con 5 pazienti o meno, togliere il verdetto o dichiararlo non calibrato.

**4. `src/07_robustness.py:487-499`**
- Criterio: 5.
- Motivo: la colonna `significativo = lo > 0` viene calcolata anche con 4-5 pazienti. In `excess_by_compartment_pair_clone_b.csv` la coppia Adjacent_normal-Tumor (5 pazienti) risulta True con [0,0034; 0,0465], cioè un IC che `scpipe` stesso dichiara non calibrato.
- Correzione: togliere la colonna, oppure segnarla come non calibrata quando i pazienti sono 5 o meno.

**5. `src/07_robustness.py:265-268`; `src/08_disambiguate.py:279` e `:379`**
- Criterio: 2.
- Motivo: test di Mann-Whitney con la cellula come unità, salvati in tabella senza che il README dichiari l'unità né il motivo:
  - `p_umi_*` in `ambient_discriminant.csv`;
  - `p_mannwhitney_CD8A` in `marker_test_disc_vs_conc.csv` (p = 0,0 su 166 contro 5.860 cellule);
  - `p_greater` in `doublet_tcr_pos_vs_neg.csv`.
- Correzione: test aggregato per paziente o per clone, oppure togliere i p, oppure dichiarare nel README unità e motivo.

**6. `src/07_robustness.py:111`, `src/08_disambiguate.py:108`, `src/09_arbitrate.py:99`**
- Criterio: 8.
- Motivo: con `--tcrdir` vuoto, `pd.concat([])` dà `ValueError: No objects to concatenate` senza un messaggio chiaro (riprodotto).
- Correzione: `if not frames: raise SystemExit(...)`, come in 10:143.

**7. `src/10_loco.py:674`**
- Criterio: 8.
- Motivo: se `paired_diffs` è vuoto per tutti gli schemi, `diffs.contrast` dà `AttributeError` (riprodotto).
- Correzione: colonne fisse nel DataFrame vuoto, oppure un controllo prima dell'uso.

**8. `src/08_disambiguate.py:403`; `src/09_arbitrate.py:245` e `:303`**
- Criterio: 8.
- Motivo: su tabelle vuote:
  - `part3` di 08 dà `KeyError: 'class'` (riprodotto);
  - `partB` di 09 dà `KeyError: 'class'` (riprodotto);
  - in `partA` di 09, `res.tissue` va in crash se nessun compartimento ha almeno 3 pazienti.
- Correzione: uscire con un messaggio chiaro.

**9. `src/08_disambiguate.py:9-11`**
- Criterio: 6.
- Motivo: la docstring elenca PA04, PA06, PA08, PA09 e PA10 per il contrasto sui tre compartimenti, ma PA04 non ha PBMC. I pazienti effettivi sono 4 (`triple_clones_summary.csv`: n_patients = 4; log di 08).
- Correzione: elencare PA06, PA08, PA09 e PA10.

**10. `src/18_loyo_check.py:83` (colonna `p_full` di `loyo_clonal_expansion_*.csv`)**
- Criterio: 6.
- Motivo: stampa «p=0.0021 (riferimento, da 17)», ma 17 e README:36 riportano p = 0,0009. 18 lo ricalcola con 10.000 permutazioni invece di 20.000.
- Correzione: leggere il p da `clonal_expansion_tests.csv`, oppure usare 20.000 permutazioni, oppure scrivere «ricalcolato con 10.000».

**11. `README.md:35` e `src/09_arbitrate.py:19`**
- Criterio: 5.
- Motivo: «misure fisiche che non usano il clustering». Le misure sono fisiche, ma i gruppi confrontati (NK TCR-/TCR+, CD8T TCR+) vengono dalle etichette del clustering.
- Correzione: «misure fisiche, confrontate fra classi dell'annotazione».

**12. `README.md:243`**
- Criteri: 4 e 5.
- Motivo: «copriva solo l'88-89% con 10-11 pazienti … copre il 92-94%». Questi valori non sono in nessuna tabella: vengono dalla fase 0 (`REGISTRO.md`). Contraddicono README:232 e `sim_a_cluster_bootstrap.csv` per k = 10 e 11:
  - percentile semplice: 0,8875-0,9225;
  - intervallo espanso: 0,935-0,9475.
- Correzione: «88,8-92,3%» e «93,5-94,8%».

**13. `RISULTATI_FINALI.md` (assente) e `src/18_loyo_check.py`**
- Criterio: 4.
- Motivo: il blocco B chiede che le sensibilità B1-B6 siano riportate in `RISULTATI_FINALI.md`, e il criterio 4 lo nomina. Il file non esiste: né nella repo, né in git, né sulla macchina.
  - B6 non è riportato da nessuna parte. In `error_rate_identity_rules.csv`, per il Tumor: CD8B > 0 dà 0,206 [0,114; 0,310]; CD8A ≥ 2 dà 0,206 [0,113; 0,307].
  - B4 manca per 18: niente media nulla, niente p centrato, nessuna prova con più semi o con 100.000 permutazioni.
- Correzione: creare il file con B1-B6 presi dalle tabelle; in 18 salvare `null_mean` e `p_centered`.

**14. `README.md:25`**
- Criterio: 4.
- Motivo: «la logistica L2 ottiene 0,851» su task2_hier3, ma `model_zoo_summary_sc.csv` (task2_hier3, fold_lopo, logreg_l2) dà 0,85049, che si arrotonda a 0,850. Anche `model_comparison_task2_hier3.csv` dà 0,8504863. I delta (+7,2 ecc.) sono corretti.
- Correzione: 0,850.

**15. `src/11_model_zoo.py:97-98`**
- Criteri: 6 (anche 5).
- Motivo: `warnings.filterwarnings("ignore", category=UserWarning)` è globale. Il commento dice che serve per i nomi di feature di LightGBM, ma in sklearn 1.9.0 `ConvergenceWarning` deriva da `UserWarning` (verificato: MRO `ConvergenceWarning -> UserWarning`). Così sparisce ogni mancata convergenza di SGD (`max_iter=50`), LinearSVC, MLP e lbfgs. È lo stesso problema che A3 ha corretto in 12, e la convergenza non è registrata da nessuna parte.
- Correzione: filtrare solo il messaggio di LightGBM; registrare la convergenza per fold.

**16. `src/14_cnv_loci.py:285-289`**
- Criterio: 5.
- Motivo: stampa «i loci PDAC separano (AUC 0.68) … Il segnale citogenetico e' locus-specifico». Ma l'AUC entro paziente dello stesso composito è 0,540 (da 0,512 a 0,582, `cnv_loci_auc_by_patient.csv`); lo 0,68 pooled mescola pazienti diversi.
- Correzione: messaggio condizionale, basato anche sull'AUC entro paziente.

**17. `src/13_infercnv.py:138-146` e `:264-266`; `src/14_cnv_loci.py:116-118` e `:196-197`**
- Criteri: 7 (anche 8).
- Motivo:
  - Con `--method auto`, se BioMart fallisce lo script passa in silenzio a copykat. `infercnvpy.tl.copykat` richiede `rpy2` e i pacchetti R `copykat` e `stringr`, assenti da `environment.yml` e dall'ambiente (verificato: `No module named 'rpy2'`). Quindi il ripiego va in ImportError.
  - BioMart interroga la release Ensembl corrente: le posizioni dei geni non sono riproducibili esattamente.
- Correzione: fallire con un messaggio chiaro invece di passare a copykat (oppure aggiungere rpy2 e R all'ambiente); fissare la release Ensembl o versionare la tabella delle posizioni.

**18. `src/04_splits.py:110-115`**
- Criterio: 8.
- Motivo: se manca il tessuto Adjacent_normal, `tis.get('Adjacent_normal', 0) > 0` è un bool e `.sum()` solleva `AttributeError` (riprodotto su h5ad sintetico). Un task vuoto arriverebbe poi a `StratifiedGroupKFold` senza un messaggio chiaro.
- Correzione: `tis.reindex(columns=[...], fill_value=0)`; saltare con un messaggio i task con meno di 2 pazienti.

**19. `src/11_model_zoo.py:251`, `src/12_feature_selection.py:262`, `src/20_xai_per_patient.py:185`**
- Criterio: 8.
- Motivo: il risultato di `get_indexer` non viene controllato. Un cell_id assente dà -1, e `X[-1]` prende in silenzio l'ultima cellula (es. split del dataset completo con un h5ad creato da `02 --subsample`). 05:150 ha già questo controllo.
- Correzione: `if (idx < 0).any(): raise SystemExit(...)`.

## BASSA (criterio violato: nessuno)

20. **`src/sim/sim_b_reference_schemes.py:100-110`**
    - Motivo: etichette e marcatori sono indipendenti dato il lignaggio, quindi la circolarità dello schema naive non è simulata. Naive risulta il meno distorto (bias 0-0,0003), e README:233 lo riporta senza dichiarare il limite.
    - Correzione: dichiarare il limite.

21. **`src/10_loco.py:100-101`, `:343`, `:383`**
    - Motivo: il commento dice «sotto questo numero», ma il codice usa `<= SMALL_K`. `frac_boot_gt0` sui tassi d'errore, che sono sempre ≥ 0, non dice nulla (0,9425 in `error_rate_blood.csv`).
    - Correzione: correggere il commento; calcolarlo solo sulle differenze.

22. **`src/19_ni_score.py:10-12` contro `:205-207`**
    - Motivo: H4 è formulata in due modi diversi («rimodellato» e «immunosoppresso → T esauste»).
    - Correzione: una sola formulazione.

23. **`src/06_clonotype_consistency.py:12-17`, `:230`, `:265`, `:531`**
    - Motivo: le soglie (cloni con almeno 2 cellule; almeno 4 nello split-half) non sono documentate. La figura non ha il suffisso `{cc}`, quindi `clone_ab` sovrascrive quella di `clone_b`.
    - Correzione: documentare le soglie e aggiungere il suffisso.

24. **`src/07_robustness.py:257-262`**
    - Motivo: condizione con precedenza fra `or` e `and` poco chiara e un test `in` ridondante.
    - Correzione: semplificarla.

25. **`src/download_tcr.sh:31-36` e `:72`**
    - Motivo: senza rete esce con 1 anche se i 30 file sono già presenti; il `find` finale conta anche file VDJ che non sono nella lista.
    - Correzione: controllare la rete solo se c'è da scaricare; contare solo i file della lista.

26. **`src/sim/sim_a_calibration.py:118`, `:140`, `:149`, `:198`; `src/sim/sim_c_block_permutation.py:51`**
    - Motivo: funzioni senza docstring; UserWarning soppressi globalmente.
    - Correzione: aggiungere le docstring e limitare il filtro.

27. **`README.md:100` e `src/scpipe.py` (`cluster_bootstrap`)**
    - Motivo: l'intervallo «expanded percentile» rimanda alla §6, ma lì la calibrazione riguarda solo la media pooled con 10-11 pazienti. È applicato anche agli IC di Spearman (16, 17, 19), a d1/d2 e ai casi con 4-5 pazienti.
    - Correzione: dichiarare per quali stimatori la copertura non è verificata.

28. **`src/13_infercnv.py:282-289`**
    - Motivo: `malignant` parte da 0, quindi le duttali dei pazienti saltati (< 20 immunitarie) diventerebbero «non maligne». Oggi nessun paziente viene saltato.
    - Correzione: inizializzare a NaN ed escludere quelle cellule.

29. **`src/13_infercnv.py:258-261`; `src/14_cnv_loci.py:188-193`**
    - Motivo: con una cache incompleta, in 13 le cellule mancanti diventano NaN e quindi «non maligne» senza avviso; in 14 una cache vecchia dà KeyError.
    - Correzione: verificare che la cache copra tutte le cellule.

30. **`src/12_feature_selection.py:123-138` e `README.md:26`**
    - Motivo: il Jaccard più basso della L1 su task2_v1 viene da soluzioni SAGA mai arrivate a convergenza, quindi il confronto è confuso con la mancata convergenza.
    - Correzione: dichiararlo, oppure aumentare `max_iter`.

31. **`src/20_xai_per_patient.py:270-272`; `src/sim/sim_d_xai_stability.py:109-110`**
    - Motivo: il controllo casuale condivide più cellule fra repliche di quante ne condividano due fold LOPO (dati: circa 5.900 contro 5.400; simulazione: 1.920 contro 1.800). Per questo la differenza è +0,030 anche senza effetto paziente.
    - Correzione: dichiararlo, oppure costruire il controllo con la stessa sovrapposizione.

32. **`src/15_inventory.py:49-62` e `:122-123`**
    - Motivo: `KNOWN_USED` è scritto a mano. La frase «ogni colonna di obs e' gia' entrata nell'analisi» è esagerata (`probe`, `pct_mt`, `barcode`). Una colonna tutta NaN viene descritta come «0 valori distinti (alta cardinalita')».
    - Correzione: scrivere «colonne considerate» e gestire a parte il caso con 0 valori distinti.

33. **`README.md:174`**
    - Motivo: la riga di 05 nella struttura dei risultati non elenca `baseline_gap_sc`, che pure è la fonte di Δ = −0,146 e del «50 volte».
    - Correzione: aggiungerla.

34. **`src/14_cnv_loci.py:270-276`**
    - Motivo: nella riga `PDAC_COMPOSITE`, `n_windows` contiene il numero di loci (8), non di finestre.
    - Correzione: una colonna separata `n_loci`.

35. **`src/13_infercnv.py:46-47` e `:374-388`**
    - Motivo: i file `external_test_*` hanno un nome che suggerisce una validazione esterna, mentre la docstring dice che non è indipendente.
    - Correzione: rinominarli (es. `adjacent_ductal_*`) o commentare il nome.

36. **`src/01b_tcr_check.py:140-150`; `src/01_probe_labels.py:214`**
    - Motivo: «stima ortogonale del tasso di discordanza» poggia su una sola libreria; «efficienza di cattura 60-85%» non ha fonte.
    - Correzione: attenuare il primo messaggio; citare la fonte o togliere l'intervallo.

37. **`src/11_model_zoo.py:426-427` e `:307`**
    - Motivo: con `--update` su uno stesso task con schemi diversi, `set_index("task")` ha duplicati e `compare()` mescola gli schemi.
    - Correzione: indicizzare per (task, schema).

38. **`src/11_model_zoo.py:322`**
    - Motivo: Nadeau-Bengio è calibrato in sim_a con StratifiedKFold k = 5 e n = 90; in 11 è usato su LOPO con test molto diseguali e training sottocampionato.
    - Correzione: una riga nella docstring o nel README.

39. **`src/05_baseline.py:91` (`n_classes`) e `rng` in `fit_explanation` di 20**
    - Motivo: parametri non usati.
    - Correzione: toglierli.

40. **`src/scpipe.py:7`**
    - Motivo: l'elenco degli script che importano scpipe è sbagliato (manca 06; 14 e 15 non lo importano; sim_b e sim_d non lo importano direttamente).
    - Correzione: allinearlo.

41. **`README.md:35`**
    - Motivo: l'elenco dei tassi di TCR nelle classi non-T (NK, Endothelial, Stellate, Platelet) sembra completo, ma omette Endocrine 11,7% (7/60) e Acinar 8,0% (24/300) di `tcr_rate_by_celltype.csv`.
    - Correzione: «fra cui», oppure aggiungerle.

42. **`README.md:33`**
    - Motivo: «simile con la sola catena TRB e con la coppia TRA+TRB». Il +0,059 è già la sola catena TRB, quindi la frase fa pensare a due varianti in più quando ce n'è una sola.
    - Correzione: «con la coppia TRA+TRB: +0,065».

43. **`results/tables/model_zoo_{sc,task1}.csv`; `README.md:165`**
    - Motivo: in fase 1 lo script 11 è stato eseguito solo con `--from-table`, quindi le tabelle per fold vengono dal codice precedente. Ho verificato che i numeri sono equivalenti: nessun fold manca di una classe. Il «circa un'ora per task» non ha un log che lo confermi.
    - L'esecuzione completa di 11 risulta in corso secondo `REGISTRO.md`, ma l'unione dei risultati (`merge11.py` nello scratchpad) non è un comando del README.
    - Correzione: completare l'esecuzione e ottenere le tabelle con un comando del README.

44. **`envs/conda-lock.yml:24` (`sources: ../../../../../home/nemo/...`)**
    - Motivo: percorso legato alla macchina locale. Inoltre `.pre-commit-config.yaml` usa ruff v0.15.21 mentre l'ambiente ha 0.15.20.
    - Correzione: rigenerare il lock dalla radice con un percorso relativo; allineare le versioni di ruff.

45. **`tests/test_leakage.py:136-237`**
    - Motivo: solo il test di 05 ha un controllo positivo (alterare il training cambia i parametri); i test di 11, 12, 13 e 20 no.
    - Correzione: aggiungere il controllo positivo.

46. **`tests/test_smoke.py:275`**
    - Motivo: `path.stem if path.parent == SRC else path.stem` ha due rami identici.
    - Correzione: semplificare.

47. **`README.md:49`**
    - Motivo: «ogni script … risponde a `--help`», ma `scpipe.py` no (è escluso anche dallo smoke test).
    - Correzione: «ogni script eseguibile».

48. **`src/13_infercnv.py` e `src/14_cnv_loci.py` (cache di inferCNV)**
    - Motivo: `cnv_scores.csv` e `cnv_matrix.npz` sono dell'11 luglio, e la riesecuzione di fase 1 le ha riusate. Il passaggio inferCNV non è stato rieseguito con il codice attuale; il README lo dice (riga 165) ma non lo segnala fra i limiti.
    - Correzione: rieseguire senza cache nella fase finale.

## Conteggio
- **ALTA:** 1
- **MEDIA:** 18
- **BASSA:** 29
- **pytest** (`-m "not slow"`): 98 passati, 0 falliti, 0 saltati (6 deselezionati, quelli lenti)
- **ruff:** 0 segnalazioni

Le prove su input vuoti e su dati sintetici sono in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/reviewer/tcr/check_empty.py` e in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/reviewer/ml/`.

Non ho modificato nessun file della repository. `revisione/REGISTRO.md` risulta modificato nel working tree, ma non da questa revisione.

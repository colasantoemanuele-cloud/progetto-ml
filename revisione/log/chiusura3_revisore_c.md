# Giro di chiusura 3 (giro 6) - area (c) CNV, simulazioni, test

**Problemi trovati, area (c): 13, 14, src/sim, tests, pyproject**

1. `src/14_cnv_loci.py:419-431`
- gravita': MEDIA
- criterio violato: 6 (il messaggio contraddice il codice; tocca anche il 5)
- motivo: il primo ramo dell'ESITO stampa "i loci PDAC separano le duttali tumorali dalle adiacenti anche entro paziente, meglio dello score medio". Il secondo stampa "entro paziente ... i loci PDAC non fanno meglio dello score medio". La condizione pero' confronta le AUC **pooled**: `comp_auc > mean_auc + 0.05`. L'AUC entro paziente dello score medio, quella che servirebbe, e' gia' calcolata (`auc_MEAN_|CNV|` in `cnv_loci_auc_by_patient`) ma non viene usata. Sui dati attuali i due criteri danno risposte opposte: pooled 0,698 contro 0,531 (i loci "vincono"), entro paziente 0,557 contro 0,572 (lo score medio fa meglio). Il ramo non scatta con i dati attuali (within = 0,557 < 0,6), quindi i log e le tabelle sono giusti; il messaggio e' sbagliato.
- correzione proposta: nella condizione usare `by_pat["auc_MEAN_|CNV|"].mean()` al posto di `mean_auc`. In alternativa, scrivere nel messaggio che il confronto con lo score medio e' pooled.

2. `src/13_infercnv.py:253-257` e `src/14_cnv_loci.py:140-144`
- gravita': MEDIA
- criterio violato: 8
- motivo: se il sottoinsieme non contiene cellule immunitarie, `ref` e' una lista vuota. Se tutte le cellule sono di riferimento, le categorie di `ref_group` sono tutte nella lista. In questi casi `cnv.tl.infercnv(..., reference_cat=[])` va in crash con un errore poco chiaro. Prova su un AnnData sintetico con soli `Ductal` (script nello scratchpad): `ValueError: need at least one array to concatenate`. Il controllo "nessuna duttale" (13:461) c'e', quello sul riferimento manca, anche in 14. Con i dati attuali non succede.
- correzione proposta: prima di inferCNV aggiungere `if not ref: raise SystemExit("[ERRORE] nessuna cellula immunitaria di riferimento ...")`, in 13 e in 14.

3. `tests/test_stats.py:6-8`, `src/sim/sim_a_calibration.py:20-21`, `README.md:237,249-250`
- gravita': BASSA
- criterio violato: nessuno (nessun numero cambia)
- motivo: secondo il testo la copertura con 11 pazienti e' verificata "con numerosita' diseguale come nei dati (5-59 cloni per paziente)". `SIZES` = numeri di cloni dei 10 pazienti di `clonal_expansion_joined` (verificato: 17, 26, 47, 42, 22, 10, 30, 59, 29, 5), piu' un undicesimo valore inventato (40). L'IC con 11 pazienti del README e' quello dell'eccesso, calcolato su `pairwise_discordance_clone_b`: 1.114 coppie, da 13 a 282 per paziente (62, 24, 184, 75, 282, 118, 146, 57, 97, 56, 13), molto piu' sbilanciate. Ho ripetuto la simulazione con queste numerosita' (400 repliche, 1.000 bootstrap, semi 1011 e 2011). La copertura e' 0,938-0,943 (continuo) e 0,928-0,930 (binario): dentro la banda, quindi la conclusione regge, ma "come nei dati" non e' vero per k = 11.
- correzione proposta: per k = 11 usare le numerosita' reali delle coppie (o dichiarare la differenza) e scrivere da dove vengono le numerosita' per k = 10 e per k = 11.

4. `README.md:237` (sez. 6, punto a) e `tests/test_stats.py:257-264`
- gravita': BASSA
- criterio violato: nessuno
- motivo: il README dice che Nadeau-Bengio e' "verificata (simulazione e test con banda)", cioe' due verifiche. Il test pero' usa gli stessi semi della simulazione: `55_000 + k` (55_003, 55_005) e `51_005`, come `sim_nb`. Siccome `nb_replicate(rep, seed0, ...)` dipende solo da `seed0 + rep` e da `rep`, le 400 repliche H0 e le 100 con effetto del test sono le prime della simulazione: stessi dati, non una seconda verifica. Il docstring del test non lo dice (per il cluster bootstrap invece la differenza di semi e' dichiarata).
- correzione proposta: dichiarare nel docstring del test e nel README che il test e' un sottoinsieme della simulazione, oppure usare semi diversi nel test.

5. `README.md:17,21` (riga 3 della tabella e nota)
- gravita': BASSA
- criterio violato: nessuno
- motivo: il confronto "pooled 0,698 → entro paziente ≈ 0,56" non e' a parita' di pazienti: il pooled usa 13 pazienti, l'entro paziente 4. Sugli stessi 4 pazienti (PA06, PA08, PA09, PA10; 4.870/333 duttali) l'AUC pooled dei loci e' 0,664 (calcolata da `cnv_matrix.npz`). La conclusione regge, ma il divario dovuto ai pazienti e' 0,664 → 0,557. Il README cita anche lo score medio pooled (0,531) e non quello entro paziente: 0,572, piu' alto dei loci (0,557; PA09: 0,766 contro 0,528). Chi legge conclude che i loci aggiungono informazione, ma entro paziente non e' cosi'.
- correzione proposta: aggiungere nella nota della riga 3 l'AUC entro paziente dello score medio (0,572, gia' in tabella). Se si cita lo 0,664 pooled sui 4 pazienti, salvarlo prima in `cnv_loci_auc` (criterio 4).

6. `src/sim/sim_a_calibration.py:28-29`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice "`H0_mean0` ed `effect_mean006` usano gli stessi dati di `coverage_continuous`, traslati". `effect_mean006` ha media 0,06, come `coverage_continuous`, quindi usa dati identici, non traslati. Solo `H0_mean0` e' traslato (media 0).
- correzione proposta: "`H0_mean0` usa gli stessi dati traslati di −0,06; `effect_mean006` usa gli stessi dati, identici".

7. `src/13_infercnv.py:79-81` (contro le righe 589 e 645-676)
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice che, se nessun fold LOPO e' valutabile, "i file corrispondenti non vengono scritti e quelli di un'esecuzione precedente vengono cancellati". In quel caso pero' `malignancy_lopo.csv` viene comunque scritto (vuoto, riga 589) e `malignancy.png` viene scritta con il pannello destro vuoto. Viene cancellato solo `malignancy_lopo_summary.csv`.
- correzione proposta: precisare nella docstring "non viene scritto `malignancy_lopo_summary.csv`; `malignancy_lopo.csv` resta vuoto", oppure cancellare anche `malignancy_lopo.csv`.

8. `tests/test_leakage.py:9-10` (contro la riga 286)
- gravita': BASSA
- criterio violato: nessuno
- motivo: il docstring dice "Controllo positivo (ogni script)". Per 04 il codice lo chiama "controllo negativo: lo split casuale mette ogni paziente in piu' fold". Nella sostanza e' un controllo positivo (mostra che il controllo puo' fallire), ma i due testi non coincidono.
- correzione proposta: uniformare il nome ("controllo positivo: lo split casuale ...").

9. `revisione/RISULTATI_FINALI.md:380` (generato da `revisione/genera_risultati.py`)
- gravita': BASSA
- criterio violato: nessuno
- motivo: "+0,143 con effetto paziente (+0,116; +0,181)" ha il formato degli IC del resto del file. In realta' e' il minimo e il massimo delle 5 repliche di `sim_d_xai_stability`. Il README scrive correttamente "da +0,116 a +0,181".
- correzione proposta: "(da +0,116 a +0,181 sulle 5 repliche)".

10. `README.md:250`
- gravita': BASSA
- criterio violato: nessuno
- motivo: "La copertura e' verificata solo per la media pooled di un esito per cella". Le simulazioni (`test_stats`, `sim_a`) hanno come osservazioni cloni o coppie clone-compartimenti, non cellule.
- correzione proposta: "un esito per clone (o coppia clone-compartimenti)".

11. `src/13_infercnv.py:280-304` (impronta usata anche da 14)
- gravita': BASSA
- criterio violato: nessuno (la versione e' fissata in `environment.yml`: infercnvpy 0.6.1)
- motivo: l'impronta registra solo `window_size`, non la versione di infercnvpy e gli altri parametri di default (`step=10`, `lfc_clip=3`, `dynamic_threshold=1.5`, `exclude_chromosomes=('chrX','chrY')`, verificati con `inspect.signature`). Con un'altra versione le cache verrebbero riusate senza avviso.
- correzione proposta: aggiungere all'impronta `infercnvpy.__version__` e i parametri passati o di default.

**Controlli senza problemi**
- **Numeri:** tutti quelli del README (riga 3, note, sez. 5 e 6) e di RISULTATI_FINALI (sezioni CNV e Simulazioni) coincidono con le tabelle.
  - CNV: 8.370/334, 13/5 pazienti, 0,837 ± 0,071 su 12 pazienti, 30,8%, 41-58%, 0,531/0,698, 0,557 (0,528-0,588), 4.870/333.
  - Simulazioni: (a) 2,7/1,8%, 4,7/6,1%, 99,3%, 93,5-94,8%, 88,8-92,3%, 91,5-92,5%, 82-83%; (b) entro 0,005, 0,126/0,093, 0,385/0,092/0,442, entro 0,008; (c) tutte le percentuali con gli IC, 4/6/11/5 e 2; (d) +0,143, +0,030 e gli intervalli.
- **Cache:** `mean|X_cnv|` di `cnv_matrix.npz` coincide con `cnv_scores.csv` (differenza massima 1,4e-16). Le finestre sono solo chr1-22.
- **Riesecuzioni:** 13, 14 e le simulazioni sono stati modificati prima dell'ultima esecuzione (date dei file contro `_status.txt`). Le modifiche successive al calcolo delle cache (chiusura 1) toccano solo docstring e cancellazione dei file. L'albero e' pulito; nessun file estraneo in `results/` per la mia area.
- **Correzioni:** A1, A2, A7 (14), D2 13:24 e blocco C applicate. Le correzioni della chiusura 2, area (c) (1-10), sono tutte applicate.
- **Test:** ogni stimatore citato come validato ha un test con banda: cluster bootstrap con k = 10/11, block permutation, Spearman con n = 10, Nadeau-Bengio con k = 3/5 e potenza, discordanza/`d_within` (verifiche esatte). I test di leakage coprono 05, 11 (percorsi sparso e SVD denso), 12, 13 (fold e modello finale), 20 e gli split di 04, con controllo positivo. La banda di Fase 0 (0,877-0,892) corrisponde a `fase0_pytest_slow.log`.
- **Pacchetti:** tutti quelli importati nella mia area sono in `envs/environment.yml` (anndata, scanpy, scipy, sklearn, matplotlib, infercnvpy, requests, pytest).

**Conteggio**
- ALTA: 0
- MEDIA: 2
- BASSA: 9
- pytest `-m "not slow" -p no:cacheprovider`: 102 passati, 0 falliti, 0 saltati (6 deselezionati, i lenti), in 109,7 s.
- ruff `check src tests revisione`: "All checks passed!"

Script di controllo (nessuna scrittura nella repo):
- `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_c/chk_cnv.py`
- `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_c/chk_cb.py`

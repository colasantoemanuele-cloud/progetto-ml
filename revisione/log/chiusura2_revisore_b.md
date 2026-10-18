# Giro di chiusura 2 - area (b) TCR

ELENCO DEI PROBLEMI (area b, giro di chiusura 2)

1. `src/19_ni_score.py:110-122` (`load_table_s1a`)
- gravita': ALTA (crash; si verifica solo con un foglio vuoto o malformato, non con il `mmc2.xlsx` pubblicato)
- criterio violato: 8
- motivo: se il foglio `Table S1A` ha l'intestazione ma nessuna riga di dati, `data` e' vuoto, `pd.DataFrame(data)` non ha colonne e la riga 122 va in crash con `AttributeError: 'DataFrame' object has no attribute 'PatientID'`. L'ho riprodotto con un xlsx sintetico (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing2_b/edge.py`). Con meno di 2 righe, `rows[1]` (riga 112) da' invece `IndexError`.
- correzione proposta: prima di leggere `rows[1]`, controllare `len(rows) < 3`; dopo la riga 120, controllare `if not data`. In entrambi i casi `SystemExit` con un messaggio chiaro ("Table S1A senza righe di dati").

2. `src/16_pseudobulk_patient.py:21-24`, `README.md:39`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la frase "la deriva ipotizzata (CD8 -> CD4/NK) va nel senso opposto agli errori osservati" e' vera solo in parte. La parte CD8 -> CD4 e' opposta (0 casi osservati, mentre CD4 -> CD8T sono 38). Ma i 12 errori CD8 -> NK (21% dei 56 errori, `error_decomposition_blood.csv`, `transition_CD8T_to_NK`) vanno proprio nel senso ipotizzato. La docstring li elenca e nella stessa frase dice "opposto".
- correzione proposta: "va nel senso opposto alla maggior parte degli errori (38 CD4 -> CD8T, 0 CD8 -> CD4); solo i 12 CD8 -> NK vanno nel senso ipotizzato".

3. Uscite anticipate che lasciano in `results/` i file di un'esecuzione precedente. In tutti questi casi lo script esce senza cancellare i propri file, mentre 08 e 16 li cancellano:
- `src/06_clonotype_consistency.py:495-497`: restano `clone_switching_{summary,across,within}_{cc}` e la figura.
- `src/07_robustness.py:458-462` (`pw.empty`): restano `pairwise_discordance_{cc}.csv`, che 08 e 17 leggerebbero senza accorgersene, e `excess_*`.
- `src/07_robustness.py:578-581`: resta `discordant_clone_markers_{cc}`.
- `src/09_arbitrate.py:197-199` e `484-487`: restano `arbitration.csv`, `error_rate_by_compartment.csv` e la figura.
- `src/16_pseudobulk_patient.py:376-377`: resta `pseudobulk_patient_pca.png` se la PCA viene saltata.
- `src/19_ni_score.py:226-227`: restano `ni_score_test.csv` e `ni_score_holm_joint.csv`.
- gravita': BASSA
- criterio violato: nessuno (con i dati attuali nessuna di queste uscite si verifica)
- correzione proposta: come in 08 (`drop_stale`), `unlink(missing_ok=True)` delle tabelle e figure della parte non eseguita.

4. `src/15_inventory.py:167-178`
- gravita': BASSA
- criterio violato: nessuno
- motivo: `used_tcr` contiene `cdr3` e `raw_clonotype_id`, che vengono stampati come "(usata)". Nessuno script li usa: `grep` li trova solo in 15. La docstring di 06 (riga 34) dice esplicitamente che i cloni non si definiscono da `raw_clonotype_id`.
- correzione proposta: togliere `cdr3` e `raw_clonotype_id` da `used_tcr`.

5. `src/08_disambiguate.py:47-49` e `README.md:167`, contro il codice di `src/08_disambiguate.py:602-611`
- gravita': BASSA
- criterio violato: nessuno
- motivo: docstring e README dicono che la cache Scrublet si riusa "solo se contiene esattamente le cellule dell'h5ad". Con `--skip-scrublet` pero' il codice controlla solo la copertura (`gex.index.difference(ds.index)`): una cache con cellule in piu' viene accettata. Il controllo esatto c'e' solo in `run_scrublet`, cioe' senza `--skip-scrublet`.
- correzione proposta: usare lo stesso controllo `set(ds.index) != set(gex.index)` anche nel ramo `--skip-scrublet`.

6. `src/08_disambiguate.py:21-24, 46`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la descrizione di `marker_test_disc_vs_conc.csv` non cita le colonne `mean_GNLY_disc` e `mean_GNLY_conc`, da cui vengono i valori "GNLY 1,78 contro 0,44" del README (riga 34). Non dice neppure che i tessuti con meno di 20 cellule CD8T per gruppo vengono saltati (righe 339-341).
- correzione proposta: aggiungere GNLY all'elenco delle colonne e la soglia di 20 cellule.

7. `src/09_arbitrate.py:7-15` (codice alla riga 187)
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring descrive l'identita' come media "sui compartimenti in cui il clone compare". Non dice che `arbitration.csv` e `error_rate_by_compartment.csv` usano solo i cloni con almeno 3 cellule in almeno 2 compartimenti (`ident.n_comp >= 2`), cioe' il filtro aggiunto con la correzione D2 "09:171".
- correzione proposta: aggiungere "solo cloni con >= 3 cellule in almeno 2 compartimenti".

8. `README.md:178` (tabella "Struttura dei risultati", riga 06)
- gravita': BASSA
- criterio violato: nessuno
- motivo: `clone_switching_within_clone_b` e' descritto come "una riga per coppia entro compartimento". In realta' ha una riga per replica split-half: 45.120 righe con colonne `patient, clone, tissue, lab_a, lab_b, switch`, cioe' (clone, tessuto) x 20 repliche, come dice la docstring di 06.
- correzione proposta: "una riga per replica split-half (clone x compartimento x 20)".

9. `revisione/RISULTATI_FINALI.md:187` (generato da `revisione/genera_risultati.py:680`)
- gravita': BASSA
- criterio violato: nessuno
- motivo: per lo split-half di 06 scrive "(14 pazienti)" senza la distinzione del README (riga 33). Da `clone_switching_summary_clone_b.csv`: `n_patients` = 14, `n_patients_across` = 11, `n_patients_within` = 14. Il problema 4 della chiusura 1 e' stato corretto nel README ma non nel generatore.
- correzione proposta: nel generatore, "14 pazienti: 11 nel termine fra compartimenti, 14 nel basale".

10. Messaggi quando rho non e' definito: `src/17_clonal_expansion.py:155, 181`, `src/19_ni_score.py:234`, `src/16_pseudobulk_patient.py:296`
- gravita': BASSA
- criterio violato: nessuno
- motivo: con esito o covariata costante, `permutation_corr_test` e `block_permutation_test` restituiscono NaN, ma 17 stampa "direzione opposta", 19 stampa "direzione OPPOSTA a H4" e 16 scrive `coerente_con_ipotesi = False`. L'ho verificato su dati sintetici con esito costante (`run_analysis` di 17: rho = nan, poi "direzione opposta").
- correzione proposta: se `np.isnan(rho)`, stampare "rho non definito" e scrivere NaN nella colonna `coerente_con_ipotesi`.

11. `revisione/note_esecuzione.md:40` (e `RISULTATI_FINALI.md:424`)
- gravita': BASSA
- criterio violato: nessuno (nessun numero del README dipende da `doublet_tcr_pos_vs_neg.csv`)
- motivo: la cache `results/tables/doublet_scores.csv` ha data 2026-07-10, prima della revisione. Il "controllo passato" verifica solo che l'insieme delle cellule coincida, non che i punteggi siano quelli del codice attuale. La verifica parziale della chiusura 1 (2 librerie ricalcolate, differenza circa 1e-16) non e' documentata.
- correzione proposta: documentare quella verifica in `note_esecuzione.md`, oppure rieseguire 08 senza cache.

12. Figure con IC non calibrati senza avviso: `src/07_robustness.py:313-333`, `src/09_arbitrate.py:287-300`
- gravita': BASSA
- criterio violato: nessuno
- motivo: `excess_by_pair_{cc}.png` mostra l'IC anche per Adjacent_normal-PBMC (4 pazienti) e Adjacent_normal-Tumor (5 pazienti), mentre la tabella mette `significativo` = NaN e `ic_calibrato` = False. `error_rate_by_compartment.png` mostra l'IC di Adjacent_normal (5 pazienti, `ic_calibrato` = False). Nessuna delle due figure segnala che quegli IC non sono calibrati.
- correzione proposta: segnare con un marcatore o una nota le barre con <= `SMALL_K` pazienti.

VERIFICHE SENZA PROBLEMI
- Numeri dell'area b confrontati con `results/tables/`, circa 100, nessuna discrepanza:
  - README, sezione "Validazione clonale con il TCR": eccesso, rumore di base, clone_ab, pesi uguali, split-half, schema solo sangue, differenze appaiate, scomposizione, CD8A, CD4, GNLY, regole d'identita', tassi di TCR, d1 e d2 in tutte le varianti, espansione clonale, leave-one-patient-out, sensibilita' di semi e permutazioni;
  - README, H1-H4 e righe corrispondenti della tabella delle fonti;
  - sezione 6 del README: valori naive +0,048/+0,070 e p = 0,0030;
  - `RISULTATI_FINALI.md`, sezioni TCR, doppietti, espansione clonale e ipotesi, compresi i valori di leave-one-patient-out con 100.000 permutazioni (0,0011-0,011, p centrato fino a 0,016 e 0,017).
- Correzioni della mia area applicate e corrette:
  - blocco A: A4, A5, A6 (`download_tcr.sh`, 16:217, 07:511, 10:319), A7 (08, 10, 19);
  - blocco B: B1-B6;
  - blocco D2: Monte Carlo al posto di "esatta", leave-one-patient-out, "ipotesi E", CP10k, limiti di H1 (circolarita', non pseudobulk), messaggi condizionali di 07, 08 e 10, "per costruzione" attenuato, 06 allineato a 10, d2 non indipendente, soglia su 11 cellule, regola dei pareggi, filtro di 09.
  - Restano a posto anche i problemi della chiusura 1 per l'area b: niente ripiego di 16, nome della colonna di Holm in 19, ordine di 07 nel README, ripresa dei download in `download_tcr.sh`.
- Codice eseguito: tutti gli script dell'area e `scpipe.py` hanno mtime 02:42-02:45. L'ultima esecuzione (chiusura 1) e' iniziata alle 02:51, il commit 59098b6 e' delle 03:05 e l'albero di lavoro e' pulito. Tabelle e figure dell'area (02:51-02:59) vengono quindi dal codice attuale, tranne la cache Scrublet (punto 11). Nessun .err dell'area contiene errori: in 19 ci sono solo avvisi di openpyxl; tutti rc = 0.
- Comandi del README per 06-10, 15-19 e 21: gli argomenti esistono. `download_tcr.sh --help` funziona. Tutti i pacchetti importati sono in `envs/environment.yml`.
- Join GEX-VDJ: nessun barcode duplicato per (paziente, tessuto) in `obs` (lettura `backed`), e c'e' una sola libreria per (paziente, tessuto).

CONTEGGIO: ALTA 1 / MEDIA 0 / BASSA 11
- pytest (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 deselezionati perche' lenti).
- ruff (`check src tests revisione`): "All checks passed!".

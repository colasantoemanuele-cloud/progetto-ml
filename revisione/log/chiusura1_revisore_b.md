# Giro di chiusura 1 - area (b) TCR

ELENCO DEI PROBLEMI (area b)

1. `README.md:35` (e `src/10_loco.py:559-561`)
- gravita': MEDIA
- criterio violato: 5
- motivo: il README dice che le misure "sono compatibili con una **maggioranza di doppietti T+NK** fra le NK-TCR+". Questo vale solo per la stima pooled su tutti i tessuti (d1 0,723 [0,556; 0,848], d2 0,676 [0,559; 0,781]). Nel Tumor, cioe' il compartimento dell'analisi d'errore, `doublet_mixture.csv` da':
  - mediane: d1 = 0,510 [0,186; 0,658] e d2 = 0,419 [0,289; 0,567];
  - medie: d1 = 0,381, con IC [0,170; 0,473] tutto sotto 0,5.
  Il README riporta l'intervallo per compartimento solo per le mediane ("da 0,51 a 0,84", "da 0,42 a 0,85"). Delle medie riporta solo la stima pooled (0,682 / 0,693), non il caso del Tumor in cui l'IC esclude la maggioranza. Anche il messaggio di 10 alle righe 559-561 ("la maggioranza delle NK-TCR+ e' compatibile con doppietti") guarda solo la riga median/all.
- correzione proposta: limitare l'affermazione al pool e aggiungere: "nel Tumor le stime vanno da 0,38 a 0,52 e con le medie l'IC di d1 (0,17-0,47) e' sotto 0,5". In 10, dare il messaggio per strato, oppure avvisare quando uno strato ha un IC sotto 0,5.

2. `src/21_fig_blood_scheme.py:19-20`
- gravita': MEDIA
- criterio violato: 6
- motivo: la docstring dice che `error_rate_blood.csv` contiene, "con <= 5 pazienti, anche ricampionature distinte e frazione > 0". Ma `10_loco.report()` (righe 343-346) aggiunge di proposito solo `n_distinct_resamples`: il commento dice che per un tasso la frazione > 0 non e' informativa. La tabella infatti ha le colonne `tissue,n,n_pat,error,lo,hi,n_distinct_resamples`, senza `frac_boot_gt0`.
- correzione proposta: "con <= 5 pazienti anche il numero di ricampionature distinte".

3. `src/08_disambiguate.py:49-50`
- gravita': MEDIA
- criterio violato: 6 (e 4: possono restare tabelle di un'esecuzione precedente)
- motivo: la docstring dice "se un'analisi non ha dati, la sua tabella non viene scritta e quella di un'esecuzione precedente viene cancellata". E' vero solo per la parte 1 (righe 196-197) e per la parte 3 con tabella vuota (riga 471). Negli altri casi:
  - `part2` con gruppi insufficienti (righe 286-288) esce senza cancellare `markers_disc_vs_conc.csv` e `marker_test_disc_vs_conc.csv`;
  - con `res` vuoto, la riga 368 scrive comunque un `marker_test_disc_vs_conc.csv` vuoto;
  - con `--pairwise` assente (riga 552) non si cancella nulla;
  - con `--skip-scrublet` e senza cache (riga 569) resta il vecchio `doublet_tcr_pos_vs_neg.csv`.
- correzione proposta: in ognuno di questi rami fare `unlink(missing_ok=True)` delle tabelle della parte, e spostare la scrittura della riga 368 dopo il controllo `res.empty`. In alternativa, limitare la frase della docstring alle parti 1 e 3.

4. `README.md:33`, `src/06_clonotype_consistency.py:301-320, 473, 483`
- gravita': BASSA
- criterio violato: nessuno
- motivo: l'eccesso split-half di 06 confronta due popolazioni diverse di pazienti e di cloni.
  - La discordanza fra compartimenti viene da 11 pazienti (1.114 coppie di cloni con almeno 2 cellule in due compartimenti).
  - Il basale viene da 14 pazienti (45.120 repliche, cloni con almeno 4 cellule in un tessuto, condivisi o no).
  - PA01, PA02 e PA15 entrano solo nel basale (verificato su `clone_switching_{across,within}_clone_b.csv`).
  Il "14 pazienti" del README e di `RISULTATI_FINALI.md`, e il k = 14 dell'intervallo espanso, contano quindi anche pazienti che non contribuiscono al termine fra compartimenti. Il dato e' dichiarato esplorativo e lo stimatore principale e' quello di 07.
- correzione proposta: scrivere "11 pazienti per la discordanza fra compartimenti, 14 per il basale", oppure restringere il basale ai pazienti e ai cloni del confronto fra compartimenti.

5. `README.md:165`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la nota "`07` va eseguito con `clone_b` dopo `clone_ab` perche' `08` e `17` leggono `pairwise_discordance_clone_b.csv`" da' una ragione che non giustifica nessun ordine:
  - le tabelle di 07 che dipendono dalla definizione di clone hanno il suffisso `{cc}`;
  - quelle in comune (`ambient_discriminant`, `chat_sensitivity_*`, `ambient_vs_real_tcr.png`) non dipendono da `--clone-def`.
  L'ordine e' indifferente; basta eseguire entrambe le chiamate prima di 08.
- correzione proposta: "08 e 17 leggono `pairwise_discordance_clone_b.csv`: eseguire 07 con `clone_b` prima di 08".

6. `src/07_robustness.py:523`
- gravita': BASSA
- criterio violato: nessuno
- motivo: il commento parla di "espressione proteica/trascrizionale" di CD4 e CD8A/CD8B, ma lo script usa solo conteggi di RNA (CP10k da `adata.X`). Non ci sono dati di proteina.
- correzione proposta: "espressione (RNA, CP10k) di CD4 e CD8A/CD8B".

7. `src/download_tcr.sh:4-5, 76` (righe 113-114 e 184 dell'output di `cat -n`)
- gravita': BASSA
- criterio violato: nessuno
- motivo: lo script promette "ripresa dei download parziali" e dice "rilancia lo script: riprende da dove ha lasciato". Pero' dopo un errore il file parziale viene cancellato (`rm -f "$dest"` sia per un file corrotto sia per un download fallito). Al nuovo avvio quel file riparte da zero: l'opzione `-C -` non trova mai un parziale da riprendere. La ripresa avviene solo per file interi, perche' i file gia' integri vengono saltati.
- correzione proposta: scrivere "i file gia' scaricati e integri vengono saltati". In alternativa, tenere il parziale quando curl fallisce, cosi' che `-C -` lo riprenda.

8. `src/08_disambiguate.py:411-413` (e `revisione/note_esecuzione.md:42`, `RISULTATI_FINALI.md`, "Scorciatoie usate")
- gravita': BASSA
- criterio violato: nessuno (vedi la verifica qui sotto)
- motivo: 08 riusa `results/tables/doublet_scores.csv` se esiste, senza controllare che copra le cellule attuali. Le cellule mancanti diventano NaN e vengono scartate in silenzio con `dropna`. La cache in uso ha data 07-10, anteriore alla revisione. Per inferCNV la cache e' stata ricalcolata e confrontata; per Scrublet no.
  Ho ricalcolato Scrublet con il codice attuale (`sc.pp.scrublet(batch_key="library", random_state=0)`) sulle due librerie piu' piccole, PDAC10-PBMC (1.801 cellule) e PDAC14-Tumor (2.127 cellule). La differenza massima dalla cache e' circa 1e-16, nessuna cellula manca. La cache quindi sembra riprodotta, ma la verifica non e' documentata.
- correzione proposta: controllare la copertura dell'indice come in 13 (fermarsi se l'indice della cache non coincide con `obs_names`). Documentare la verifica in `note_esecuzione.md`, oppure rieseguire 08 senza cache.

9. `src/09_arbitrate.py:203-228` (`error_rate_by_compartment.csv`)
- gravita': BASSA
- criterio violato: nessuno (la tabella non e' citata nel README)
- motivo: l'IC dell'errore viene salvato anche per Adjacent_normal, che ha 5 pazienti, senza segnalare che l'IC non e' calibrato. Al contrario, 10 (`report`) e 08 aggiungono la diagnostica con <= `SMALL_K` pazienti. La docstring dice solo "IC 95% con cluster bootstrap sui pazienti". La differenza appaiata Tumor - PBMC di 09 (righe 233-249) e' solo stampata.
- correzione proposta: aggiungere `n_distinct_resamples` o una colonna `ic_calibrato` con <= `SMALL_K` pazienti, e la stessa avvertenza nella docstring.

10. `src/16_pseudobulk_patient.py:143-153, 259-263`; `src/19_ni_score.py:246-270`
- gravita': BASSA
- criterio violato: nessuno
- motivo:
  - (a) Se `arbitration_blood.csv` manca, 16 passa da solo a un altro esito (l'eccesso medio di 07), con un semplice avviso. I file e le colonne mantengono i nomi `error_rate` e `pseudobulk_hypothesis_tests.csv`. 19 invece non ha questo ripiego, quindi la correzione di Holm congiunta mescolerebbe due esiti diversi.
  - (b) Se ci sono meno di 6 pazienti, 16 esce senza cancellare il vecchio `pseudobulk_hypothesis_tests.csv`, e 19 lo ricaricherebbe come H1-H3.
  - (c) Senza `--previous-tests`, 19 chiama comunque la colonna `p_holm_family_of_4`.
  - (d) La docstring di 16 non dice che le covariate richiedono anche almeno 30 CD8T tumorali (`MIN_CD8T_TUMOR`).
- correzione proposta: togliere il ripiego, o fermarsi con un errore; cancellare le tabelle vecchie nell'uscita anticipata; dare alla colonna un nome che dipenda da quanti test ci sono; citare `MIN_CD8T_TUMOR` nella docstring.

VERIFICHE SENZA PROBLEMI
- Numeri del README e di `RISULTATI_FINALI.md` dell'area b confrontati con `results/tables/`: TCR, doppietti, espansione clonale, H1-H4 e le righe corrispondenti della tabella delle fonti. Circa 90 numeri, nessuna discrepanza.
- Correzioni dell'elenco applicate: A4, A5, A6 (download_tcr.sh, 16:217, 07:511, 10:319), A7 (08, 10, 19), B1-B6, e le D2 della mia area (Monte Carlo, leave-one-patient-out, "ipotesi E", CP10k, H1, messaggi condizionali di 07, 08 e 10, pareggi, filtro di 09).
- Tutte le tabelle e figure della mia area sono state rigenerate in Fase 3 al commit 61fb36b, che e' l'ultimo commit che modifica questi script.
- Tutti i pacchetti importati sono in `envs/environment.yml`.

CONTEGGIO: ALTA 0 / MEDIA 3 / BASSA 7
- pytest (`-m "not slow"`): 102 passati, 0 falliti, 0 saltati (6 deselezionati perche' lenti).
- ruff (`check src tests revisione`): "All checks passed!".
- Script di controllo (Scrublet): `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing_b/scrub.py`.

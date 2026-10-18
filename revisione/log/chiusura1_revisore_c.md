# Giro di chiusura 1 - area (c) CNV, simulazioni, test

## Problemi trovati, area (c): 13, 14, src/sim, tests, pyproject

1. `src/13_infercnv.py:177-202` (`add_positions`), effetto in `:205-223` e in `src/14_cnv_loci.py:105-119`
- gravita': MEDIA
- criterio violato: 5. Lo script chiama "score CNV / profilo CNV" una misura che contiene anche il genoma mitocondriale; la docstring non lo dice, quindi tocca anche il 6.
- motivo: la tabella Ensembl chiama il mitocondrio `MT`, e il codice lo trasforma in `chrMT`. infercnvpy 0.6.1 esclude solo `chrM`, `chrX` e `chrY` (`_infercnv.py:327`, `exclude_chromosomes`). Per questo i 13 geni MT- e un contig con un solo gene (`GL000194.1`, MAFIP) entrano in inferCNV come cromosomi. Evidenza:
  - In `cnv_matrix.npz` i nomi delle finestre comprendono `{'MT': 1, 'GL000194.1': 1}` su 1.674.
  - La finestra MT ha media |x| = 0,344, contro una media di 0,0078 su tutte le finestre (44 volte di piu'). La sua correlazione con `pct_mt` (metrica di QC) e' r = 0,90.
  - Ho ricalcolato senza le due finestre, con la stessa soglia mediana + 3·MAD (script `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing_c/mt_check.py`, che usa `classify_lopo` di 13 e non scrive in `results/`):
    - cambiano 393 etichette su 8.704 duttali;
    - le duttali adiacenti "maligne" passano da 45-62% a 42-59% (PA06 0,615→0,538, PA08 0,452→0,424, PA10 0,578→0,590);
    - l'AUC LOPO passa da 0,829 ± 0,069 a 0,818 ± 0,074.
  - I numeri del README (righe 17 e 21) e di RISULTATI_FINALI.md dipendono quindi in parte dalla frazione mitocondriale. Le conclusioni qualitative restano le stesse.
- correzione proposta: in `add_positions` tenere solo i cromosomi `^(\d+|X|Y)$` (oppure passare `exclude_chromosomes=("chrX","chrY","chrMT")`). Poi cancellare le cache, rieseguire 13 e 14 e aggiornare README e RISULTATI_FINALI. In alternativa, dichiarare l'inclusione e il suo effetto.

2. `src/14_cnv_loci.py:120-122` e `:140`
- gravita': BASSA
- criterio violato: nessuno. Nell'ambiente dichiarato `require_copykat` esce prima con un messaggio chiaro.
- motivo: con `--method copykat` il ramo `else` di `main` (`:215-216`) non chiama `add_positions`, quindi `a.var` non ha la colonna `chromosome`. Dopo copykat, `var["chromosome"]` (riga 140) solleva KeyError. Il crash avviene solo se `rpy2` e R sono installati.
- correzione proposta: chiamare `m13.add_positions(sub)` anche per copykat, oppure ricostruire `win_pos` solo quando le coordinate esistono.

3. `src/14_cnv_loci.py:370-390`
- gravita': BASSA
- criterio violato: nessuno
- motivo: se l'AUC entro paziente e' ≥ 0,65 ma il composito pooled non supera lo score medio di 0,05, il codice finisce nel ramo "separazione entro paziente parziale". Il messaggio e' sbagliato per quel caso. Con i dati attuali non scatta (0,54).
- correzione proposta: aggiungere un ramo per "entro paziente ≥ 0,65 ma loci non migliori dello score medio".

4. `src/14_cnv_loci.py:339-346` e `revisione/RISULTATI_FINALI.md:313`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la riga `_mean` fa la media anche di `n_tumor` e `n_normal` (1.217,5 e 83,25). RISULTATI_FINALI la riporta come "Duttali Tumor 1.218 / Adjacent 83", che non sono conteggi di cellule.
- correzione proposta: mettere le somme o lasciare vuoti `n_tumor`/`n_normal` nella riga `_mean`.

5. `src/14_cnv_loci.py:285-292`
- gravita': BASSA
- criterio violato: nessuno
- motivo: gli z-score del composito sono standardizzati su tutte le righe di `X`, comprese le 92.515 immunitarie di riferimento, non sulle sole duttali. La docstring (righe 13-14) dice solo "somma degli z-score".
- correzione proposta: scrivere nel commento e nella docstring su quali cellule sono stimate media e deviazione standard.

6. `src/sim/sim_a_calibration.py:160-165`
- gravita': BASSA
- criterio violato: nessuno
- motivo: `rng = default_rng(1000 + k)` riparte uguale per ogni scenario. `H0_mean0` ed `effect_mean006` usano quindi gli stessi dati simulati di `coverage_continuous` (solo traslati, o con un altro valore vero). Per costruzione il tasso sotto H0 vale 1 − copertura: 0,18 = 1 − 0,82; 0,0625 = 1 − 0,9375; 0,065 = 1 − 0,935. Docstring e RISULTATI_FINALI le presentano come verifiche separate.
- correzione proposta: dichiararlo nella docstring, oppure usare semi diversi per scenario.

7. `src/sim/sim_b_reference_schemes.py:44-45`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dell'output elenca "errore stimato medio, errore vero medio, repliche". Il CSV contiene anche `est_error_max`, `mean_n` e `bias`, che README (0,442 "in una replica") e RISULTATI_FINALI ("Distorsione") usano.
- correzione proposta: elencare tutte le colonne.

8. `README.md:234` contro `tests/test_stats.py:222-231`
- gravita': BASSA
- criterio violato: nessuno. Lo stimatore ha comunque un test con banda per k = 5.
- motivo: il README dice che la calibrazione di Nadeau-Bengio e' "verificata con ... k = 3 e 5". Il test con banda copre solo k = 5. Il caso k = 3 (2,7%) e' solo in `sim_a`, senza banda.
- correzione proposta: aggiungere k = 3 al test, oppure scrivere "testata con k = 5, simulata anche con k = 3".

9. `tests/test_stats.py:9-11` contro `README.md:247`
- gravita': BASSA
- criterio violato: nessuno
- motivo: per la copertura dell'intervallo percentile semplice il docstring del test dice 0,877-0,892, il README 88,8-92,3%. Sono due simulazioni diverse (Fase 0, semi 100+k; sim_a, semi 1000+k), ma nessuno dei due testi dice da quale simulazione viene il proprio valore. Inoltre in sim_a un caso (0,9225, binario, k = 11) rientra nella banda [0,90; 0,99].
- correzione proposta: indicare in entrambi i punti la fonte (`fase0_pytest_slow.log` oppure `sim_a_cluster_bootstrap.csv`).

10. `src/13_infercnv.py:346-362`, `src/14_cnv_loci.py:192-204`, `README.md:165`
- gravita': BASSA
- criterio violato: nessuno
- motivo: il controllo della cache confronta solo gli identificativi delle cellule. Se cambiano i conteggi di `sc_raw.h5ad`, le coordinate dei geni o i parametri di inferCNV, una cache vecchia viene riusata senza avviso. Il README presenta il controllo come garanzia ("solo se coprono esattamente le cellule selezionate"). Oggi la cache e' stata ricalcolata nel giro 1 ed e' identica byte per byte, quindi i numeri non sono toccati.
- correzione proposta: salvare nella cache un'impronta (hash di h5ad e coordinate, `window_size`, riferimento) e fermarsi se non coincide.

## Controlli senza problemi
- **Numeri:** README (righe CNV 17, 21, 222-226; sezione 6) e RISULTATI_FINALI.md (sezioni CNV e Simulazioni) coincidono con `results/tables`.
- **Correzioni A1, A2, A7 (14) e D2 (13:24):** applicate.
  - LogCPM e HVG sono dentro la Pipeline, anche nel modello finale.
  - `--splits` e' stato tolto.
  - Il testo sul riferimento comune e sulla soglia stimata sulle stesse cellule e' corretto; il 45-62% e' verificato.
  - Il commento sull'interpolazione delle finestre c'e'.
  - Le AUC per paziente sono salvate.
- **Esecuzione:** la Fase 3 ha girato con il codice 61fb36b. Da 61fb36b a HEAD non cambia nulla in `src`, `tests` o `pyproject.toml`.
- **Test di leakage:** coprono 05, 11, 12, 13, 20 e gli split di 04, ciascuno con il suo controllo.
- **Bande di calibrazione:** ci sono per cluster bootstrap, block permutation, Spearman, Nadeau-Bengio e discordanza.
- **Pacchetti:** tutti quelli importati nella mia area sono in `envs/environment.yml`.

## Conteggio
- ALTA: 0
- MEDIA: 1
- BASSA: 9
- pytest `-m "not slow"`: 102 passati, 0 falliti, 0 saltati (6 deselezionati, i lenti).
- ruff `check src tests revisione`: "All checks passed!"

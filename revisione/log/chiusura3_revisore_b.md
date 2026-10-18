# Giro di chiusura 3 (giro 6) - area (b) TCR

ELENCO DEI PROBLEMI (area b, giro di chiusura 3)

1. `src/19_ni_score.py:245-246` (con `src/scpipe.py:301-311` e `334-346`, e la correzione di Holm a `19_ni_score.py:275-302`)
- gravita': ALTA. Non succede con il `mmc2.xlsx` attuale: i 14 NI score ci sono tutti.
- criterio violato: 8
- motivo: se manca un valore di NI score per un paziente nella Table S1A, `astype(float)` lo trasforma in NaN. Il controllo di `permutation_corr_test` (`np.ptp(rx) == 0`) non se ne accorge, perche' `ptp` di un vettore con NaN e' NaN. Quindi rho = NaN e tutte le rho permutate sono NaN. Il confronto `|rhos| >= |rho_obs|` e' sempre falso, e il p vale 1/(n_perm+1) = 0,00005, cioe' il minimo possibile.
  - L'ho riprodotto su una copia di `mmc2.xlsx` con il NI score di PA12 cancellato, usando `load_table_s1a`, `build_outcome` e Holm di 19 (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_b/nan19.py`).
  - Risultato: `rho nan, p 4.9998e-05`. Holm sulla famiglia di 4 da' H4 p_adj = 0,0002, che lo script stamperebbe come "SIGNIFICATIVO". Cambiano anche i p_Holm di H1-H3 (0,352 / 0,339 / 0,673).
  - `block_permutation_test`, che usano 17 e 18, fa lo stesso: `p = p_centered = 0,0005` con rho NaN (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_b/nan_test.py`).
  - Inoltre, se un p di H1-H3 fosse NaN, Holm lo esclude (m = 3), ma la colonna si chiamerebbe comunque `p_holm_family_of_4`.
- correzione proposta:
  - in `permutation_corr_test` e `block_permutation_test`, se x o y contiene NaN o valori non finiti: avviso e (NaN, NaN), come gia' si fa per l'input costante;
  - in 19, prima del test, togliere i pazienti con NI score mancante stampandone gli ID, o fermarsi con un messaggio chiaro;
  - chiamare la colonna con il numero di p definiti.

2. `src/21_fig_blood_scheme.py:86-102`, `src/10_loco.py:718-733`
- gravita': BASSA
- criterio violato: nessuno
- motivo: nella chiusura 2 gli IC non calibrati sono stati segnati nelle figure di 07 e 09, ma non in queste due:
  - `error_rate_blood.png`, la figura dello schema primario, mostra l'IC di Adjacent_normal con 4 pazienti (`error_rate_blood.csv`: n_pat = 4, 35 ricampionature distinte) come gli altri; scrive solo "(4 pz)";
  - `error_rate_loco.png` fa lo stesso con Adjacent_normal a 5 pazienti (`error_rate_schemes.csv`, loco/main: n_pat = 5).
- correzione proposta: barra grigia tratteggiata e legenda quando `n_pat <= SMALL_K`, come in `09_arbitrate.py:293-322`.

3. Uscite in cui restano file di un'esecuzione precedente, o viene scritto un file vuoto, nonostante le docstring:
- `src/07_robustness.py:477-490`: con `pw.empty` lo script esce prima di disegnare `ambient_vs_real_tcr.png` (righe 612-646). La sezione C ha dati e le sue tabelle vengono riscritte, ma resta la figura vecchia. La docstring (righe 49-50) dice che le figure vecchie si cancellano solo per le sezioni senza dati.
- `src/07_robustness.py:568`: con `--no-markers` resta il vecchio `discordant_clone_markers_{cc}.csv`.
- `src/07_robustness.py:560-562`: se nessuna coppia ha almeno 3 pazienti, viene scritto un `excess_by_compartment_pair_{cc}.csv` vuoto, senza colonne.
- `src/16_pseudobulk_patient.py:263-270`: con meno di 6 pazienti cancella le due tabelle ma non `pseudobulk_patient_pca.png`.
- `src/10_loco.py:718-719`: se `res[1]` e' vuoto resta il vecchio `error_rate_loco.png`.
- `src/09_arbitrate.py:259-263`: se nessun compartimento ha almeno 3 pazienti, scrive `error_rate_by_compartment.csv` con la sola intestazione. Alle righe 389-392, con `out` vuoto, resta il vecchio `doublet_physical.png`.
- gravita': BASSA
- criterio violato: nessuno (con i dati attuali nessuno di questi rami si esegue)
- correzione proposta: `unlink(missing_ok=True)` dei file della parte non eseguita; generare la figura di 07 (C) prima del `return` della sezione A.

4. `src/10_loco.py:7-9, 278-279` (contro `src/09_arbitrate.py:12-13, 191`)
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring di 10 rimanda a 09 per f8/f4 e per il margine, e chiama `naive` lo schema circolare. Le due versioni pero' non sono uguali:
  - 09 tiene solo i cloni con almeno 3 cellule in almeno 2 compartimenti;
  - `naive` di 10 tiene anche i cloni presenti in un solo compartimento, il cui riferimento e' solo il compartimento giudicato.
  Per questo i numeri divergono: Tumor 0,144 (457 cloni, 11 pazienti) in `error_rate_by_compartment.csv` contro 0,092 (1.519 cloni, 13 pazienti) in `error_rate_schemes.csv` (naive). `RISULTATI_FINALI.md` riporta lo 0,092 come "naive" senza spiegare la differenza.
- correzione proposta: scriverlo nella docstring di 10 (e in una riga del README), oppure applicare anche a `naive` il filtro `n_comp >= 2`, che pero' cambierebbe +0,035.

5. `src/07_robustness.py:335-338` (log `revisione/log/07_robustness_clone_{ab,b}.err`)
- gravita': BASSA
- criterio violato: nessuno (pandas 2.3.3 e' fissato in `envs/environment.yml`)
- motivo: `ax.errorbar` riceve `xerr` come `Series` di un solo elemento (la coppia PBMC-Tumor, l'unica calibrata). Matplotlib emette `FutureWarning: Calling float on a single element Series ... will raise a TypeError in the future` (riprodotto in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_b/warn07.py`). Con un pandas piu' recente diventerebbe un crash della figura.
- correzione proposta: passare `pp.excess.to_numpy()` e `xerr=[(pp.excess - pp.lo).to_numpy(), (pp.hi - pp.excess).to_numpy()]`.

6. `src/09_arbitrate.py:150`
- gravita': BASSA
- criterio violato: nessuno
- motivo: `clone_identity` e' annotata `-> pd.DataFrame` ma restituisce la tupla `(per, ident)` (riga 182).
- correzione proposta: `-> tuple[pd.DataFrame, pd.DataFrame]`.

7. `README.md:250`, `src/06_clonotype_consistency.py:20-23`
- gravita': BASSA
- criterio violato: nessuno (06 e' dichiarato esplorativo)
- motivo: il README elenca gli IC la cui copertura non e' verificata (Spearman, d1/d2, 4-5 pazienti) ma non l'IC split-half di 06. Questo IC e' una differenza fra due medie pooled calcolate su insiemi di pazienti diversi (11 e 14), con k = 14 nell'intervallo espanso, e non e' la "media pooled" verificata in `tests/test_stats.py`. La docstring di 06 rimanda a `scpipe.cluster_bootstrap`, ma 06 usa una funzione propria, `cluster_excess_ci`.
- correzione proposta: aggiungere l'IC di 06 all'elenco del README e scrivere nella docstring che la copertura di questa statistica non e' verificata.

8. `src/06_clonotype_consistency.py:439, 48`, `src/07_robustness.py:416`, `README.md:36`
- gravita': BASSA
- criterio violato: nessuno
- motivo: `has_tcr` e' `clone_b.notna()`, cioe' richiede una catena TRB produttiva. Le 2.690 cellule con sola TRA produttiva (su 125.970 con contig, `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing3_b/tra_only.py`) contano come TCR-. README ("con un TCR produttivo", NK 5,9% ecc.) e docstring ("frazione di cellule con TCR") non lo dicono.
- correzione proposta: "con catena TRB produttiva" nel README e nelle docstring di 06 e 07 (`tcr_rate_by_celltype`, c_hat).

VERIFICHE SENZA PROBLEMI
- Numeri del README e di `RISULTATI_FINALI.md` dell'area b confrontati con `results/tables/`, circa 110, nessuna discrepanza:
  - README "Validazione clonale con il TCR": 0,107; eccesso +0,059 [+0,031; +0,110]; clone_ab, pesi uguali, split-half con 14/11/14 pazienti; schema solo sangue; differenze appaiate (35 ricampionature, 99,7%); scomposizione 38/18/12/6/0 e 0,141; CD8A 7%/88% con 5/5 pazienti; CD4 16x; GNLY 1,78/0,44; regole CD8 0,206/0,206; tassi di TCR; d1/d2 in tutte le varianti; espansione clonale (rho, p, media nulla, p centrato, semi, 100.000 permutazioni); leave-one-patient-out (10/10, intervalli di rho, p fino a 0,012);
  - H1-H4 (rho, p, p_Holm) e le righe corrispondenti della tabella delle fonti; bootstrap ingenuo +0,048/+0,070 e p = 0,0030;
  - `RISULTATI_FINALI.md`: sezioni TCR, doppietti, espansione clonale, sensibilita' B4, leave-one-patient-out (p centrato fino a 0,016/0,017) e ipotesi.
- Correzioni della mia area applicate e corrette:
  - blocco A: A4, A5, A6 (`download_tcr.sh`, 16:217, 07:511, 10:319), A7 (08, 10, 19);
  - blocco B: B1-B6;
  - blocco D2: Monte Carlo, leave-one-patient-out, "ipotesi E", CP10k, limiti di H1 (deriva, circolarita', non pseudobulk), messaggi condizionali di 07/08/10, "per costruzione" attenuato, 06 allineato a 10, d2 non indipendente e mediane pooled, soglia su 11 cellule, pareggi e cella scartata con n dispari, filtro `n_comp >= 2` in 09.
  - Restano a posto anche le correzioni della chiusura 2 (foglio S1A vuoto, controllo esatto della cache Scrublet, rho NaN in 16/17/19, figure di 07 e 09).
- Codice e risultati:
  - tutti gli script dell'area e `scpipe.py` hanno mtime anteriore alla riesecuzione della chiusura 2 (03:25-03:35), e il commit 04ebed5 e' delle 03:36 con albero pulito;
  - tutte le tabelle e figure dell'area vengono da quella esecuzione, salvo la cache Scrublet gia' documentata;
  - nessun errore nei .err (solo avvisi di openpyxl e il FutureWarning del punto 5); tutti rc = 0.
- Comandi del README per 06-10, 15-19, 21 e `download_tcr.sh`: gli argomenti esistono e `--help` funziona. Tutti i pacchetti importati sono in `envs/environment.yml`.
- Nessun crash su input vuoti in `ambient_discriminant`, `chat_sensitivity` e `nk_contamination` (provato con un input senza TCR).

CONTEGGIO: ALTA 1 / MEDIA 0 / BASSA 7
- pytest (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 deselezionati perche' lenti).
- ruff (`check src tests revisione`): "All checks passed!".

# Giro 7 (chiusura 4) - area (b) TCR

ELENCO DEI PROBLEMI (area b, giro 7 = chiusura 4)

1. `README.md:32-35`, `src/09_arbitrate.py:7-8` (con `src/10_loco.py:7-20` e `src/15_inventory.py:182-185`)
- gravità: MEDIA, lettera e (interpretazione biologica non sostenuta), criterio 5
- criterio violato: 5
- motivo:
  - 09:7 afferma senza eccezioni che "il lignaggio CD4/CD8 di un clone T maturo non cambia fra tessuti". Il README (riga 32) usa lo stesso presupposto per dire che la CDR3 "permette di misurare l'affidabilità dell'annotazione".
  - Le eccezioni note non sono dichiarate da nessuna parte: CD4 che acquisiscono CD8αα (CD4+CD8αα+ e doppi positivi nei tessuti e nei tumori), MAIT (CD8αα, CD8αβ o doppi negativi), iNKT (CD4+ o doppi negativi). `grep -i "MAIT|NKT|CD8αα"` su README, RISULTATI_FINALI e src/06-21 trova solo il commento di 06:460 su c_hat. Le CD4 citotossiche compaiono solo come alternativa in 08 e nel README:35.
  - La regola d'identità (CD8A>0 o CD8B>0, in 09:156 e 10:245) classificherebbe come CD8 un clone CD4 con CD8αα.
  - Ho verificato che sul dataset reale le eccezioni non spostano i numeri (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_b/mait_check.py`):
    - fra i 287 cloni del Tumor (schema blood) i cloni MAIT-like (TRAV1-2 con TRAJ33/12/20) sono 12, tutti CD8T etichettati CD8T e quindi concordanti; i cloni iNKT (TRAV10-TRAJ18) sono 0;
    - le CD8T discordanti del tumore hanno CD8A>0 solo nel 6,6% (`marker_test_disc_vs_conc`), quindi non sono CD4 con CD8αα.
  - Manca quindi la dichiarazione, non la sostanza.
- correzione proposta:
  - in 09 (docstring), in 10 e nel README:32-35 scrivere "di norma stabile, con eccezioni (CD4 citotossiche, CD4+CD8αα+ e doppi positivi nei tessuti, MAIT, NKT)";
  - riportare il controllo con v_gene e j_gene (12 cloni MAIT, tutti concordanti; 0 iNKT) e il dato CD8A 7% come argomento contro CD8αα.

2. `README.md:35` ("l'etichetta è sbagliata nel 19,5% dei cloni nel Tumor"), `src/10_loco.py:271-272` + `src/scpipe.py:410-423` (stessa regola in 09:166, 16/19 come esito, 17)
- gravità: MEDIA, lettera d (affermazione nel README più forte dei dati), criterio 5
- criterio violato: 5
- motivo:
  - Nel Tumor l'etichetta di maggioranza è un pareggio, risolto dalla "prima comparsa", cioè dall'ordine dei barcode di `obs` (alfabetico, `AAACCTG...`), in 10 dei 287 cloni (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_b/tie_check.py`).
  - 9 di questi 10 cloni sono contati come errori, cioè 9 dei 56 errori: per esempio 4 cloni di PA07 con CD4T=CD8T (2/2, 2/2, 2/2, 4/4) e un clone di PA12 con CD4T=NK=19.
  - In 9 dei 10 pareggi l'etichetta pari al massimo comprende l'identità del sangue, quindi l'esito dipende solo dall'ordine dei barcode; 8 di questi 9 sono finiti come errore. Uno solo (PA11, CD8T=NK con identità CD4T) è errore comunque.
  - Ricalcolo:
    - pareggi esclusi: 0,173 [0,101; 0,255] (278 cloni);
    - pareggi contati 0,5: 0,183 [0,108; 0,266];
    - valore riportato: 0,195 [0,110; 0,290];
    - nella scomposizione, i cloni CD4 → CD8T non decisi da un pareggio sono 33 su 38.
  - La regola è documentata solo nelle docstring (correzione D2). Il README e RISULTATI_FINALI non la citano (`grep -i "pareggi|parit"`: nessun riscontro pertinente) e non riportano questa sensibilità.
- correzione proposta:
  - aggiungere a `error_rate_identity_rules.csv`, o a `error_decomposition_blood.csv`, il tasso con i pareggi esclusi (o contati 0,5) e il numero di errori decisi da un pareggio;
  - nel README:35 scrivere "0,195 (0,173 escludendo i 9 cloni con pareggio nell'etichetta di maggioranza)".

3. `src/06_clonotype_consistency.py:448`, `revisione/genera_risultati.py:837` → `revisione/RISULTATI_FINALI.md:250`
- gravità: BASSA
- criterio violato: nessuno
- motivo: dopo la correzione del giro 6 (TRB produttiva) restano due testi con "TCR produttivo": il messaggio stampato da 06 e la frase di RISULTATI_FINALI. `has_tcr` è `clone_b.notna()`, quindi le 2.690 cellule con sola TRA produttiva contano come TCR-.
- correzione proposta: "con catena TRB produttiva" in entrambi.

4. `src/10_loco.py:721, 734`
- gravità: BASSA
- criterio violato: nessuno
- motivo: la variante `no_top_decile` toglie le cellule sopra il 90° percentile di n_counts per (libreria, classe), ma l'output la chiama "escludendo i possibili doppietti" e "LOCO senza doppietti". Non sono doppietti identificati, sono le cellule con più UMI.
- correzione proposta: titolo "LOCO senza il decile superiore di n_counts (proxy dei doppietti)".

5. `src/16_pseudobulk_patient.py:93, 219`
- gravità: BASSA
- criterio violato: nessuno
- motivo: i commenti parlano di "tre ipotesi a priori", mentre la docstring (riga 15), l'output (riga 281) e il README:38 dicono "formulate dopo i moduli 08-10 ma prima di calcolare le covariate".
- correzione proposta: "tre ipotesi (formulate prima di calcolare le covariate)".

6. `src/06_clonotype_consistency.py:28-31, 108-117`
- gravità: BASSA
- criterio violato: nessuno
- motivo: la docstring definisce q come la frazione con TCR "fra le non-T", ma `NON_T` esclude Platelet (13,8% TCR+), Plasma ed Endocrine. Il valore di q ne dipende (vedi `chat_sensitivity_*`).
- correzione proposta: elencare nella docstring le classi usate per q.

7. `README.md:35-36`
- gravità: BASSA
- criterio violato: nessuno
- motivo: il README riporta i 18 cloni etichettati NK fra gli errori nel Tumor e l'errore senza di loro (0,141), ma non collega direttamente questi cloni alla stima dei doppietti, che nel Tumor non sostiene la maggioranza (d1 0,510, d2 0,419).
  - Ho fatto un controllo grezzo, senza IC, sulle 171 cellule NK dei 18 cloni (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_b/nk_err.py`): mediana n_counts 5.345 contro 3.449 delle NK TCR- e 3.150 delle CD8T TCR+ del Tumor, cioè d1 ≈ 0,60; Scrublet medio 0,062 contro 0,037.
- correzione proposta: una riga con d1 (e IC) per le cellule NK dei cloni in errore, oppure la frase "la quota NK degli errori può essere in parte dovuta a doppietti".

LIMITE NOTO (nuovo)

8. `src/07_robustness.py:119`, `src/08_disambiguate.py:118`, `src/09_arbitrate.py:110`, `src/10_loco.py:151`
- problema: `TISSUE_MAP[tis]` va in KeyError con un nome di file VDJ con un tessuto diverso da Tumor, Normal o PBMC. 06 invece stampa un avviso e salta il file (06:152-154).
- verifica sul dataset reale: i 30 file di `vdj_inventory.csv` hanno solo i tessuti Tumor, Normal e PBMC, e tutti gli script sono terminati con rc = 0 (`revisione/log/_status.txt`, giro 6).

VERIFICHE SENZA PROBLEMI
- Numeri, circa 120, tutti coincidenti con `results/tables/`:
  - README "Validazione clonale con il TCR": 0,107; +0,059 [+0,031; +0,110]; +0,065; +0,074; +0,057 con 14/11/14 pazienti; 0,195; 0,059; +0,080 con 35 ricampionature e 99,7%; 38/18/12/6/0 e 0,141; LOCO +0,027, naive +0,035; CD8B e CD8A≥2 a 0,206; 7%/88% con 5/5 pazienti; 16×; GNLY 1,78/0,44; tassi TCR; d1/d2 in tutte le varianti; espansione clonale (ρ, p, media nulla, p centrato, 5 semi, 100.000 permutazioni); leave-one-patient-out con 10/10, intervalli di ρ e p fino a 0,012;
  - H1-H4 e le righe corrispondenti della tabella delle fonti;
  - RISULTATI_FINALI, sezioni TCR, doppietti, espansione e ipotesi (p centrato fino a 0,016 e 0,017, massimo 9,0e-04).
- Correzioni A4, A5, A6 (download_tcr.sh, 16, 07, 10), A7 (08, 10, 19), B1-B6 e D2 dell'area: applicate. Sono a posto anche le correzioni del giro 6: NaN in permutation_corr_test e block_permutation_test, esclusione degli NI score mancanti in 19, `xerr` di 07 (FutureWarning sparito dai .err), uscite anticipate di 07, 09, 10 e 16, annotazione di 09, figure con IC non calibrati di 10 e 21.
- Codice eseguito uguale a HEAD: tutti gli script dell'area hanno mtime anteriore alla riesecuzione del giro 6 (04:00-04:08); tabelle e figure dell'area sono di quella esecuzione; nessun file orfano; i .err sono puliti (solo avvisi openpyxl in 19).
- `--help` funziona per tutti gli 11 script e per `download_tcr.sh`. Gli argomenti dei comandi del README esistono e l'ordine dei comandi rispetta le dipendenze. I pacchetti importati sono tutti in `envs/environment.yml`.
- La definizione di clone con la sola TRB non gonfia l'errore: la TRA dominante coincide fra PBMC e Tumor nel 96,4% dei cloni in errore e nel 96,0% dei concordanti (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_b/tra_check.py`).

CONTEGGIO: ALTA 0 / MEDIA 2 / BASSA 5 / LIMITI NOTI nuovi 1

Controlli biologici (area b):
- identità CD4/CD8 stabile nel clone, con le eccezioni: NON superato. Le eccezioni non sono dichiarate (punto 1); sul dataset reale non cambiano i numeri (MAIT 12, tutti concordanti; iNKT 0; CD8A 7% nelle discordanti).
- riferimento dell'inferenza CNV: non pertinente all'area b (13 e 14). Il README:22 dichiara il riferimento immunitario di tutti i pazienti con soglia per paziente.
- direzione dei loci PDAC: non pertinente all'area b.
- geni di stress da dissociazione fra i geni XAI: non pertinente all'area b (dichiarati nel README:28).
- confronto tumore-adiacente su 4 pazienti: superato. +0,080 con "IC non calibrato", 35 ricampionature e frazione > 0; in 08, 09 e 10 nessun verdetto con 5 pazienti o meno; nella figura `error_rate_blood.png` la barra dell'Adjacent_normal è grigia e tratteggiata.
- quota degli errori nel tumore dovuta alle cellule etichettate NK: superato nella dichiarazione (18 su 56, 0,063; errore senza NK 0,141; nel Tumor la maggioranza di doppietti è dichiarata non sostenuta). Il collegamento diretto con i doppietti non è quantificato (punto 7).

pytest (`-m "not slow" -p no:cacheprovider`): 103 passati, 0 falliti, 0 saltati (6 deselezionati perché lenti; 3.000 UserWarning di sklearn in test_leakage).
ruff (`check src tests revisione`): "All checks passed!".

# Giro 8 - area (d) README, ambiente, biologia

**1. `README.md:244` (con `README.md:213`)**
- gravità: ALTA. Verifica (b), criteri 4 e 5.
- criterio violato: 4, 5
- motivo: il §6 dice "permutazione dei cloni ignorando il paziente: p = 0,0030 contro 0,0009 del block permutation test" e conclude che le due colonne "servono a mostrare quanto cambia il risultato ignorando la dipendenza". In `clonal_expansion_tests.csv` (esito `error`) si legge invece: `clone_naive_permutation` p = 0,00045 e `clone_block_permutation` p = 0,00040. I due valori sono di fatto uguali, quindi la conclusione non regge più per l'errore. Sono valori di prima del giro 7, mai aggiornati (`git diff da4f0c2 cb2b492 -- README.md` non tocca queste righe). La riga 213 della tabella delle fonti ripete "p = 0,0030".
- correzione proposta: scrivere "p = 0,00045 contro 0,0004 del block permutation test (errore; per l'eccesso entrambi 0,00005, il minimo)". Togliere o riformulare la frase "quanto cambia" per il p, e lasciare l'esempio sull'IC (+0,048; +0,070 contro +0,031; +0,110), che è ancora vero.

**2. `README.md:212, 214, 215, 216, 217, 223` (tabella "Da dove viene ogni numero")**
- gravità: MEDIA. Verifica (b), criterio 4.
- criterio violato: 4
- motivo: le righe della tabella delle fonti riportano i numeri di prima del giro 7, che non coincidono né con il testo del README né con le tabelle:
  - 212: "+0,057 [+0,017; +0,127]". `clone_switching_summary_clone_b` dà +0,063 [+0,023; +0,134].
  - 214: "0,195 [0,110; 0,290]; 0,059 [0,000; 0,093]". `error_rate_blood` dà 0,170 [0,095; 0,252] e 0,060 [0,000; 0,095].
  - 215: "+0,080 [0,000; +0,192]; 99,7%; +0,027 [−0,010; +0,124]; +0,035 [+0,005; +0,138]". `error_rate_paired_diff` dà +0,059 [0,000; 0,160], 0,933, +0,016 [−0,020; +0,095] e +0,024 [0,000; +0,110].
  - 216: "38, 18, 12, 6, 0 cloni; 0,141 [0,058; 0,261]". `error_decomposition_blood` dà 32, 15, 11, 4, 0 e 0,122 [0,049; 0,224].
  - 217: "0,206 [0,114; 0,310]; 0,206 [0,113; 0,307]". `error_rate_identity_rules` dà 0,180 [0,099; 0,271] e 0,181 [0,098; 0,269].
  - 223: "p fino a 0,012". In `loyo_clonal_expansion_perm_sensitivity` il p massimo è 0,00155, che è il "0,0016" del testo (riga 37).

  Mancano anche le fonti dei numeri nuovi: 0,199/0,164 e +0,068/+0,057 (`error_rate_tie_sensitivity`), 10/1 pareggi, 277 cloni senza pareggi, 12 MAIT-like (`clone_invariant_tcr_check`), d1 0,617 e d2 0,458 (`doublet_mixture_error_nk`), 90 coppie e 3.943 repliche escluse.
- correzione proposta: aggiornare le righe 212-217 e 223 con i valori attuali e aggiungere le righe per le tabelle nuove.

**3. `README.md:262` (§7, punto "CNV")**
- gravità: MEDIA. Verifica (a), criteri 5 e 6.
- criterio violato: 5, 6
- motivo: la riga dice "Le coordinate delle finestre di inferCNV sono interpolate in modo uniforme per cromosoma". Dal giro 7 non è più vero:
  - le coordinate vengono dai geni reali (`window_table` di 13, `cnv_windows.csv`);
  - lo dichiarano anche `README.md:21` e la docstring di `14_cnv_loci.py:19-23`;
  - il log di 14 conferma che CDKN2A, a chr9:21,97-22,00 Mb, cade in 10 finestre, tutte attribuite a `CDKN2A_9p21`.

  La riga contraddice quindi il codice e la riga 21 dello stesso README.
- correzione proposta: "Le finestre hanno le coordinate reali dei loro 100 geni (`cnv_windows.csv`), ma sono larghe (mediana 14 Mb; 15-28 Mb su 9p): il segnale di un locus piccolo è quello della regione circostante."

**4. `src/sim/sim_b_reference_schemes.py:160-161`; `README.md:240`; `revisione/RISULTATI_FINALI.md:351-374` (colonne "Errore vero" e "Distorsione")**
- gravità: MEDIA. Verifica (a), cioè un consumatore di `arbitrate` di 10, e verifica (c). Criteri 5 e 6.
- criterio violato: 5, 6
- motivo: dopo il giro 7 `arbitrate` restituisce le unità in pareggio con `lab = None` ed `error = NaN`. In sim_b, `est_error = g.error.mean()` le esclude, mentre `true_error = j.lab != j.true_identity` vale True per ogni `None`: ogni pareggio conta come errore vero. Stima e verità usano quindi denominatori diversi, e la distorsione negativa è prodotta da questo confronto. Ho riprodotto H0 con gli stessi semi (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing5_d/simb_check.py`):
  - con il codice attuale la distorsione va da −0,0120 a −0,0063, uguale alla tabella;
  - misurando l'errore vero sulle stesse unità senza pareggio va da 0,0000 a +0,0050.

  La spiegazione del README ("le stime escludono i pareggi ..., che nella simulazione sono più spesso errori") è falsa: nel codice i pareggi sono sempre errori per definizione. In più "−0,013" non è nella tabella: il minimo sotto H0 è −0,01205, cioè −0,012.
- correzione proposta: calcolare `true_error` solo sulle unità con `~j.tie` (oppure trattare i pareggi come in `tie_sensitivity`), rieseguire sim_b e aggiornare README §6 (b) e RISULTATI. In alternativa, dichiarare nel README e nella docstring che l'errore vero conta i pareggi come errori.

**5. `README.md:181, 182, 188` (§5, struttura dei risultati)**
- gravità: BASSA
- criterio violato: 6
- motivo: la tabella non elenca le uscite aggiunte nel giro 7: `arbitration_tie_sensitivity` (09), `error_rate_tie_sensitivity`, `clone_invariant_tcr_check` e `doublet_mixture_error_nk` (10), `cnv_windows` (14). Le docstring di 09, 10 e 14 invece le elencano.
- correzione proposta: aggiungerle alle righe 09, 10 e 14.

**6. `README.md:273-279` e `revisione/limiti_noti.md` (copiato in `RISULTATI_FINALI.md:441-447`)**
- gravità: BASSA
- criterio violato: 6
- motivo: dopo le modifiche del giro 7 diversi riferimenti di riga sono superati:
  - `09_arbitrate.py:110` → 117, `10_loco.py:151` → 173 e 532 (`TISSUE_MAP[tis]`);
  - `20_xai_per_patient.py:147` → 148 (`argsort`);
  - `11_model_zoo.py:477` → 478 e `20:228` → 229 (`csr_matrix`); `11:390` → 391;
  - `sim_a_calibration.py:85` → 88 e `tests/test_stats.py:60` → 62 (`SIZES`);
  - `test_stats.py:240, :247` → 242, 249; "righe 37-38" del docstring → 39-40.
- correzione proposta: aggiornare i numeri di riga, o citare i nomi delle funzioni e delle costanti al posto delle righe.

**7. `src/20_xai_per_patient.py:312-317` (con l'insieme `IMMEDIATE_EARLY`, righe ~100-118)**
- gravità: BASSA
- criterio violato: 6
- motivo: lo script stampa "geni di risposta immediata, associati allo stress di dissociazione: BTG2, CD69, CREM, DUSP4, FOS, HSPH1" (`revisione/log/20_xai_per_patient.log:42`). Il README (riga 28) e la docstring di 20 (righe 32-35), corretti nel giro 7, separano invece la firma di dissociazione (FOS, BTG2, HSPH1) dai geni che marcano anche attivazione o residenza tissutale (CD69, DUSP4, CREM). In più HSPH1 è un gene da shock termico, non di risposta immediata.
- correzione proposta: due insiemi, firma di dissociazione e geni di risposta immediata o attivazione T, con messaggi come nella docstring.

**8. `src/18_loyo_check.py:7-9`**
- gravità: BASSA
- criterio violato: 6
- motivo: la docstring dice "stimata su 287 cloni di 10 pazienti" e non dice che i cloni con esito in pareggio sono esclusi. Il codice li esclude (`loyo` e `perm_sensitivity`, `dropna`): per l'errore i cloni sono 277.
- correzione proposta: "287 cloni (277 per l'errore: i 10 con etichetta in pareggio sono esclusi, come in 17)".

**9. `README.md:94`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: la cartella radice è indicata come `progetto-ml/`, ma la repository è `pdac-ml`.
- correzione proposta: scrivere `pdac-ml/`, oppure solo "dalla radice del repository".

**10. `src/06_clonotype_consistency.py:242`**
- gravità: BASSA
- criterio violato: 6
- motivo: `def majority(labels: pd.Series) -> str:` restituisce None in caso di pareggio, come dice la sua stessa docstring. L'annotazione contraddice il comportamento.
- correzione proposta: `-> str | None`.

**Conteggio: ALTA 1 / MEDIA 3 / BASSA 6**

**Esito delle verifiche**
- **a) Correzioni del giro 7: PARZIALMENTE SUPERATA.**
  - Pareggi: esclusi in 06, 09, 10, 16, 17, 18, 19 e 21, con `n_ties` e conteggi stampati (16 e 19: "10 coppie ... escluse"; 17: 277 cloni per l'errore). Le sensibilità "pareggi = errore" e "pareggi = corretti" ci sono solo in 09 (`arbitration_tie_sensitivity`) e in 10 (`error_rate_tie_sensitivity`); non ci sono in 06 né in 16-19, e il README non dichiara che ci siano.
  - Il consumatore sim_b confronta male stima e verità (punto 4).
  - Soglia di inferCNV: una sola soglia globale, indipendente dall'ordine (`global_noise_filter`, 1,5 sd di tutta `X_cnv`).
  - Finestre CNV: coordinate reali; CDKN2A a chr9:21,97-22,00 Mb in 10 finestre, tutte in `CDKN2A_9p21`, con controllo che ferma lo script se non è così (`check_marker_windows`). Il README §7 non è stato aggiornato (punto 3).
- **b) Numeri: NON SUPERATA** per il README (punti 1 e 2, più "−0,013" del punto 4). Il testo del §1 (righe 3-44) e le simulazioni (a), (c), (d) coincidono con le tabelle. `RISULTATI_FINALI.md` è SUPERATO: rigenerato in scratch con `genera_risultati.py --strict`, è identico byte per byte al file versionato, con 0 valori mancanti, e la tabella prima/dopo coincide con le tabelle. I comandi del README corrispondono a `--help` di tutti gli script e l'ordine è eseguibile. Tutti i pacchetti importati da `src/`, `tests/` e `revisione/` sono in `envs/environment.yml`, e le versioni coincidono con il lock e con l'ambiente installato (ruff 0.15.20, come in `.pre-commit-config.yaml`). 05, 11 e 12 non rieseguiti: confermato. L'AST di `LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard` e `mean_pairwise_jaccard` è identico a 52ace44. In 11 e 12 cambiano solo argparse, controlli, `compare` (solo stampe) e la riemissione degli avvisi; 05 è stato rieseguito nel giro 7 con tabelle identiche.
- **c) Correttezza biologica, punto per punto:**
  1. Identità CD4/CD8 stabile nel clone con eccezioni dichiarate: SUPERATO. Eccezioni in README:32 e 260, 09:7-10, 10:31-36. Verifica in `clone_invariant_tcr_check`: 12 MAIT-like, tutti corretti; 0 iNKT; TRAV1-2 con TRAJ33/12/20 e TRAV10-TRAJ18 sono corretti.
  2. Marcatori: SUPERATO. CD8 = CD8A o CD8B, CD4 = CD4 senza CD8; NK = NKG7, GNLY, KLRD1, NCAM1; esaurimento = PDCD1, LAG3, HAVCR2, TOX, CTLA4, TIGIT; coerenti fra 07, 08, 09 e 10.
  3. Riferimento CNV: SUPERATO. Immunitarie di tutti i pazienti in una sola esecuzione, soglia per paziente (mediana + 3·MAD).
  4. Loci PDAC: SUPERATO. Perdite 9p21, 17p13, 18q21, 1p e 6q con direzione −1; guadagni 8q24, 12p12 e 20q con +1. Le coordinate GRCh38 contengono CDKN2A, TP53, SMAD4, MYC e KRAS.
  5. Malignità e campo di cancerizzazione come ipotesi: SUPERATO (README:22, 13:21-25, esito stampato da 14).
  6. Quota degli errori dovuta alle etichette NK: SUPERATO. 15 cloni su 277 (0,054); d1 0,617 [0,058; 0,849]; "IC troppo largo per dire che siano in gran parte doppietti".
  7. Geni di stress da dissociazione fra i geni comuni XAI: SUPERATO nel README e nella docstring; resta la nota BASSA sull'output di 20 (punto 7).
  8. Confronto tumore-adiacente su 4 pazienti: SUPERATO. "IC non calibrato", 35 ricampionature distinte; in RISULTATI asterisco sugli IC con ≤ 5 pazienti.
  9. Ipotesi H1-H4: SUPERATO. Direzione dichiarata (tutte con ρ < 0 contro un'attesa > 0), circolarità di H1 dichiarata, "non significativo non esclude un effetto".

**pytest** (`-m "not slow" -p no:cacheprovider`): 108 passati, 0 falliti, 0 saltati, 6 deselezionati (i lenti), 112 s.
**ruff** (`check src tests revisione`): "All checks passed!".

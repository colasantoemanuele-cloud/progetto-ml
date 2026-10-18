# Giro 8 - area (a) ML

**1. `src/20_xai_per_patient.py:100-118` (insieme `IMMEDIATE_EARLY`) e `:312-319` (stampa)**
- gravità: BASSA
- criterio violato: 6
- motivo: nel giro 7 la correzione è stata applicata al docstring (`:32-35`) e al README (`:28`), ma non al codice. Docstring e README separano due gruppi: FOS, BTG2 e HSPH1 appartengono alla firma di dissociazione; CD69, DUSP4 e CREM marcano attivazione, residenza tissutale ed esaurimento delle T. Il codice invece li mette in un solo insieme («Geni di risposta immediata (stress di dissociazione)»). L'output dell'esecuzione finale (`revisione/log/20_xai_per_patient.log`) stampa: «geni di risposta immediata, associati allo stress di dissociazione: BTG2, CD69, CREM, DUSP4, FOS, HSPH1». Inoltre HSPH1 è un gene da shock termico, non di risposta immediata. L'output contraddice quindi il docstring e il README. I numeri non cambiano.
- correzione proposta: dividere l'insieme in due. DISSOCIAZIONE (firma di van den Brink: FOS, FOSB, JUN, JUNB, EGR1, BTG2, DUSP1, IER2, NR4A1, HSPA1A, HSPA1B, HSPH1) e ATTIVAZIONE_T (CD69, DUSP4, CREM), con due righe di stampa separate. Poi rieseguire 20 (circa 2,5 minuti, numeri invariati).

**2. `README.md:181, 184, 187` e `revisione/limiti_noti.md:7, 10, 13` (ricopiati in `RISULTATI_FINALI.md`, "Limiti noti")**
- gravità: BASSA
- criterio violato: nessuno
- motivo: dopo le modifiche al docstring del giro 7 quattro numeri di riga sono spostati di una:
  - `src/11_model_zoo.py:477` è `ad.read_h5ad`; la conversione CSR è alla riga 478;
  - `src/11_model_zoo.py:390` è una `)`; il `to_csv` di `model_comparison` è alla riga 391;
  - `src/20_xai_per_patient.py:228` è `read_h5ad`; la conversione CSR è alla riga 229;
  - `src/20_xai_per_patient.py:147` è una riga vuota; `np.argsort(-coef)[:TOP_K]` è alla riga 148.
  - Corretti: `05:309`, `12:289`.
- correzione proposta: scrivere 478, 391, 229, 148 e rigenerare `RISULTATI_FINALI.md`.

**3. `README.md:244` (e riga della tabella delle fonti `:213`)** (fuori dalla mia area, trovato verificando b)
- gravità: ALTA (verifica b; criteri 4 e 5)
- criterio violato: 4, 5
- motivo: il testo dice «permutazione dei cloni ignorando il paziente: p = 0,0030 contro 0,0009 del block permutation test. Servono a mostrare quanto cambia il risultato ignorando la dipendenza…». Dopo il giro 7, `clonal_expansion_tests.csv` (esito `error`) dà `clone_naive_permutation` p = 0,00045 e `clone_block_permutation` p = 0,0004. Sono i valori di prima del giro 7, e il contrasto che la frase illustra non c'è più: per il p dell'errore ignorare il paziente non cambia quasi nulla. L'IC +0,048; +0,070 contro +0,031; +0,110 è invece corretto.
- correzione proposta: riportare 0,00045 contro 0,0004 e riformulare: per il p l'effetto dell'unità sbagliata non si vede più; resta visibile nell'IC dell'eccesso. Aggiornare anche la riga `:213`.

**4. `README.md:212, 214, 215, 216, 217, 223` (tabella «Da dove viene ogni numero»)** (fuori dalla mia area)
- gravità: MEDIA (verifica b, criterio 4)
- criterio violato: 4
- motivo: queste righe riportano i valori di prima del giro 7, che non sono più in nessuna tabella e non coincidono con il testo del README:

| Riga | Nel README | Nella tabella attuale |
| :--- | :--- | :--- |
| `:212` | +0,057 [+0,017; +0,127] | +0,063 [+0,023; +0,134] (`clone_switching_summary_clone_b`) |
| `:214` | 0,195 [0,110; 0,290]; 0,059 [0,000; 0,093] | 0,170 [0,095; 0,252]; 0,060 [0,000; 0,095] (`error_rate_blood`) |
| `:215` | +0,080 [0,000; +0,192]; 99,7%; +0,027 […]; +0,035 […] | +0,059 [0,000; 0,160]; 93,3%; LOCO +0,016 [−0,020; 0,095]; naive +0,024 [0,000; 0,110] |
| `:216` | 38, 18, 12, 6; 0,141 [0,058; 0,261] | 32, 15, 11, 4; 0,122 [0,049; 0,224] |
| `:217` | 0,206 / 0,206 | 0,180 [0,099; 0,271] / 0,181 [0,098; 0,269] |
| `:223` | «p fino a 0,012» | max `p` di `loyo_clonal_expansion_perm_sensitivity` = 0,00155; il testo (riga 37) dice 0,0016 |

- correzione proposta: allineare queste righe ai valori del testo, che sono già quelli delle tabelle.

**5. `README.md:262` (Limiti, «CNV»)** (fuori dalla mia area)
- gravità: MEDIA (verifica a, finestre CNV reali; criteri 5 e 6)
- criterio violato: 5, 6
- motivo: il README dice «Le coordinate delle finestre di inferCNV sono interpolate in modo uniforme per cromosoma». Dal giro 7 `14_cnv_loci.py` usa le coordinate reali dei geni (`:23`, `:60`, `:167`, `cnv_windows.csv`), come dice anche `README.md:21`. Il README si contraddice. Inoltre `cnv_windows` è citata nella tabella delle fonti ma manca nella riga 14 della tabella degli output (§5).
- correzione proposta: «coordinate reali dei geni di ogni finestra; su 9p le finestre con CDKN2A coprono 15-28 Mb, quindi il "9p21" è in pratica il braccio». Aggiungere `cnv_windows` alla riga 14 di §5.

**Conteggio:** ALTA 1 / MEDIA 2 / BASSA 2. Gli ALTA e i MEDIA (3, 4, 5) sono fuori dalla mia area: nella mia area ALTA 0, MEDIA 0, BASSA 2.

**Esito delle verifiche**
- **a)** Nella mia area `scpipe.majority_label` è corretta: pareggio → None, nessuna dipendenza dall'ordine, coperta da `test_majority_label_ties_are_undefined`. 06, 09 e 10 la usano; 20 usa `.mode()` solo per il nome del paziente escluso, con un controllo che impone un paziente per fold. Per 13 e 14 ho trovato solo il testo vecchio del punto 5; il resto è delle altre aree.
- **b)** Nella mia area tutti i numeri coincidono con le tabelle:
  - README sezione 1, righe 1-2 (`baseline_summary_sc`, `baseline_gap_sc`: 0,954 / 0,964 / Δ −0,009; 0,875 ± 0,002; 0,729 ± 0,100; Δ −0,146; rapporto 50,3);
  - «Confronto dei modelli» (11, 12, 20) e righe 203-210 della tabella delle fonti.
  - `revisione/genera_risultati.py --strict`, che ho rieseguito scrivendo in scratch, rigenera `RISULTATI_FINALI.md` identico.
  - `11 --from-table`, che ho rieseguito con TAB/FIG in scratch, riproduce identici byte per byte `model_zoo_{sc,task1}`, i due summary e i tre `model_comparison`.
  - L'unione dei `11_final/model_zoo_fin_*` secondo `_commit.txt` ricostruisce le 403 righe con differenza 0. Confermate anche le frasi sui thread: MLP 0,01396; knn 0,00012; svm_rbf 0,00038.
  - `git diff`: 12 dopo 52ace44 e 11 dopo 8058701 cambiano solo controlli, messaggi, docstring e help. Le funzioni di scpipe usate da 05, 11 e 12 non cambiano. 05 è stato rieseguito nel giro 7 con il codice attuale; le sue tabelle sono invariate in git dal commit 52ace44.
  - Fuori area: i punti 3 e 4.
- **c) punto per punto**
  - Identità CD4/CD8 stabile nel clone, con le eccezioni: non nella mia area. Il README (32, 168) dichiara le eccezioni e cita `clone_invariant_tcr_check`. OK.
  - Marcatori: nella mia area si usano le etichette degli autori. Platelet sono tutte PBMC (1.587); Endocrine in 12 pazienti, 2 con almeno 10 cellule; `ACINAR` in 20 corretto. OK.
  - Riferimento CNV (immunitarie di tutti i pazienti, soglia per paziente): README:22 corretto. Non nella mia area.
  - Loci PDAC e coordinate: il testo sulle coordinate a README:262 è vecchio (punto 5).
  - Malignità e campo di cancerizzazione presentati come ipotesi: OK (README:22).
  - Errori nel tumore da cellule etichettate NK: OK. Il README dice che l'IC (d1 0,617 [0,058; 0,849]) è troppo largo per dire «in gran parte doppietti».
  - Geni di stress da dissociazione fra i geni comuni XAI: presenti e corretti in README e docstring; l'output di 20 è incoerente (punto 1).
  - Confronto tumore-adiacente: `task1a_paired` ha 5 pazienti, riportati con media ± sd e non come IC; i loci CNV entro paziente su 4 pazienti sono dati come «≈ 0,57» senza IC. OK.
  - H1-H4: direzione opposta, circolarità di H1 e «non esclude un effetto» sono dichiarati (README 38-44, 171). OK.
- **Altro nella mia area:**
  - nessun leakage: LogCPM, HVG, SVD, scaler, prefiltro, selettori e L1 di 20 sono stimati sul solo training;
  - Nadeau-Bengio dichiarato non verificato per fold LOPO, p grezzi, minimo di Wilcoxon dichiarato;
  - docstring complete e coerenti; tutti i pacchetti importati sono in `environment.yml`;
  - correzioni A3, A6, A7 e D2/D3 della mia area applicate.

**pytest** (`-m "not slow" -p no:cacheprovider`): 108 passati, 0 falliti, 0 saltati (6 lenti esclusi). **ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni.

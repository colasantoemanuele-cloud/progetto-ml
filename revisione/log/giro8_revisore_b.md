# Giro 8 - area (b) TCR

ELENCO DEI PROBLEMI (area b, giro 8)

1. `README.md:212-217`, `README.md:223` (tabella "Da dove viene ogni numero del README") e `README.md:244` (§6, "Colonne con un'unita' diversa dal paziente")
- gravità: ALTA. Verifica (b), criterio 4.
- criterio violato: 4
- motivo: il giro 7 ha aggiornato i numeri nel testo del §1 ma non la tabella delle fonti né il §6. Questi passaggi del README riportano ancora i valori di prima del giro 7, e nessuna tabella attuale li contiene:
  - riga 212: +0,057 [+0,017; +0,127]. In `clone_switching_summary_clone_b` c'è +0,063 [+0,023; +0,134].
  - riga 213: p = 0,0030 (permutazione dei cloni). In `clonal_expansion_tests`, riga `clone_naive_permutation`, esito error, c'è 0,00045.
  - riga 214: 0,195 [0,110; 0,290] e 0,059 [0,000; 0,093]. In `error_rate_blood` ci sono 0,170 [0,095; 0,252] e 0,060 [0,000; 0,095].
  - riga 215: +0,080 [0,000; +0,192], 99,7%, +0,027 [−0,010; +0,124], +0,035 [+0,005; +0,138]. In `error_rate_paired_diff` ci sono +0,059 [0,000; +0,160], 0,933, +0,016 [−0,020; +0,095], +0,024 [0,000; +0,110].
  - riga 216: 38, 18, 12, 6, 0 e 0,141 [0,058; 0,261]. In `error_decomposition_blood` ci sono 32, 15, 11, 4, 0 e 0,122 [0,049; 0,224].
  - riga 217: 0,206 [0,114; 0,310] e 0,206 [0,113; 0,307]. In `error_rate_identity_rules` ci sono 0,180 [0,099; 0,271] e 0,181 [0,098; 0,269].
  - riga 223: p fino a 0,012. In `loyo_clonal_expansion_perm_sensitivity` il massimo è 0,00155.
  - riga 244: "p = 0,0030 contro 0,0009 del block permutation test". Oggi i valori sono 0,00045 contro 0,0004.
  - Mancano inoltre le righe di fonte per numeri nuovi del §1: 0,199 e 0,164, +0,068 e +0,057 (`error_rate_tie_sensitivity`); 12 MAIT-like e 0 iNKT (`clone_invariant_tcr_check`); 0,617 e 0,458 con i loro IC (`doublet_mixture_error_nk`); 277 cloni senza pareggi.
- correzione proposta: aggiornare queste righe con i valori attuali e aggiungere le righe di fonte mancanti. Meglio ancora, generare la tabella delle fonti con `genera_risultati.py --strict`.

2. `README.md:262` (§7 Limiti, "CNV")
- gravità: MEDIA. Verifica (a), criterio 6, applicato al README: un testo che contraddice il codice. Fuori dall'area b, ma è nel README.
- criterio violato: 6
- motivo: il §7 dice "Le coordinate delle finestre di inferCNV sono interpolate in modo uniforme per cromosoma". Dopo il giro 7, `14_cnv_loci.py:compute_cnv` usa le coordinate reali dei geni (`m13.window_table`, `results/tables/cnv_windows.csv`), e il README stesso alla riga 20 dice "ogni finestra ha le coordinate reali dei suoi geni".
- correzione proposta: "Ogni finestra di inferCNV (100 geni) ha le coordinate reali dei suoi geni, ma copre molti Mb (su 9p 15-28 Mb): i confini dei loci restano approssimati".

3. `README.md:180-182` (tabella "Struttura dei risultati", righe 09 e 10)
- gravità: BASSA
- criterio violato: nessuno
- motivo: mancano i file nuovi del giro 7, che le docstring invece elencano: `arbitration_tie_sensitivity` (09); `error_rate_tie_sensitivity`, `clone_invariant_tcr_check`, `doublet_mixture_error_nk` (10).
- correzione proposta: aggiungerli alle righe 09 e 10.

4. `src/18_loyo_check.py:7-9`
- gravità: BASSA
- criterio violato: nessuno
- motivo: la docstring dice che l'associazione con l'errore "e' stimata su 287 cloni". Dopo l'esclusione dei pareggi l'errore usa 277 cloni (`clonal_expansion_tests`, n = 277); 287 vale solo per l'eccesso.
- correzione proposta: "277 cloni per l'errore (287 per l'eccesso; 10 pareggi esclusi)".

5. `revisione/genera_risultati.py:838` → `revisione/RISULTATI_FINALI.md:251`
- gravità: BASSA
- criterio violato: nessuno
- motivo: c'è ancora "Tassi di TCR produttivo per classe". `has_tcr` è `clone_b.notna()`, cioè una TRB produttiva. La correzione era stata chiesta al giro 7 (chiusura 4 b, punto 3): è stata fatta in 06 ma non qui.
- correzione proposta: "Tassi di cellule con catena TRB produttiva per classe".

6. `src/06_clonotype_consistency.py:504-508`, `src/16_pseudobulk_patient.py:145-146`, `src/17_clonal_expansion.py:143-144`, `src/18_loyo_check.py:86-87`, `src/19_ni_score.py:165-166`
- gravità: BASSA
- criterio violato: nessuno
- motivo: i pareggi sono esclusi correttamente, ma la sensibilità "pareggi = errore / pareggi = corretti" esiste solo in 09 e 10. 06 e i consumatori 16-19 non la riportano. L'ho ricalcolata (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing5_b/ties_consumers.py`) e le conclusioni non cambiano:
  - 17, errore: ρ 0,220 / p 0,0004 con i pareggi esclusi; 0,186 / 0,00045 contati come errore; 0,217 / 0,0001 contati come corretti;
  - H1-H4, pareggi contati come errore: ρ −0,47 / −0,54 / −0,15 / −0,14, p da 0,11 a 0,71; contati come corretti: valori identici a quelli principali;
  - 06: eccesso +0,063 con i pareggi esclusi, +0,052 contati come discordanti, +0,058 contati come concordanti.
- correzione proposta: una riga di sensibilità in `clonal_expansion_perm_sensitivity` e in `clone_switching_summary`, oppure una frase nel README che riporti questi valori.

7. `src/06_clonotype_consistency.py:221`, `src/07_robustness.py:130`, `src/08_disambiguate.py:128`, `src/09_arbitrate.py:127`, `src/10_loco.py:183` e `:544` (scelta del contig dominante)
- gravità: BASSA
- criterio violato: nessuno
- motivo: `sort_values("umis", ascending=False).drop_duplicates(...)` usa quicksort, che non è stabile. Quando una cellula ha due CDR3 diverse con lo stesso numero massimo di UMI, il clone assegnato dipende dall'ordine delle righe. È lo stesso tipo di arbitrarietà corretto al giro 7 per l'etichetta di maggioranza.
  - Sul dataset riguarda 346 cellule su 123.280 per la TRB e 2.062 su 106.316 per la TRA (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing5_b/umi_ties.py`).
  - Esempio: la cellula PA12 Tumor `GGCGACTCATGCTAGT-1` ha due TRB con 2 UMI ciascuna. Nel mio ricalcolo indipendente, con una scelta diversa del contig, il suo clone ha 3 cellule e diventa un pareggio: i pareggi del Tumor salgono da 10 a 11 e la sensibilità "pareggi = errore" passa da 0,199 a 0,201.
  - I numeri principali sono identici: Tumor 277 cloni, 0,1697; Adjacent_normal 117 cloni, 0,0598 (`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing5_b/indep_blood.py`).
- correzione proposta: usare `kind="stable"` con un secondo criterio deterministico (per esempio `cdr3_nt`), oppure escludere le cellule con due catene TRB a pari UMI (sono possibili doppietti); documentare la regola.

8. `src/10_loco.py:746-748` (`mixture_error_nk`, colonna `n_patients`) → `results/tables/doublet_mixture_error_nk.csv`
- gravità: BASSA
- criterio violato: nessuno
- motivo: la tabella riporta `n_patients = 13`, che sono tutti i pazienti del Tumor ricampionati. Ma i 15 cloni NK in errore, con le loro 149 cellule, vengono da 6 pazienti soli (PA12 7, PA05 2, PA06 2, PA11 2, PA09 1, PA13 1). Inoltre l'intervallo "expanded percentile" è calcolato con k = 13.
- correzione proposta: aggiungere la colonna `n_patients_with_error_nk` = 6 e dichiarare che l'IC non è calibrato.

9. `src/10_loco.py:656-665` (`decompose_tumor_error`, riga `ties_excluded`) → `revisione/RISULTATI_FINALI.md:225` circa
- gravità: BASSA
- criterio violato: nessuno
- motivo: nella riga dei pareggi `n_errors = 0`, e RISULTATI_FINALI la stampa come "Errori 0". In realtà l'errore dei pareggi non è definito. Inoltre `n_clones` contiene il numero di pareggi (10) e `n_patients` = 10 conta tutti i pazienti, mentre i pareggi vengono da 6.
- correzione proposta: `n_errors = NaN`, una colonna `n_ties` e il numero dei pazienti con pareggi.

10. `src/06_clonotype_consistency.py:602-611`
- gravità: BASSA
- criterio violato: nessuno
- motivo: l'esito stampato dice "l'IC 95% ... esclude lo zero. Le etichette dei cloni condivisi cambiano fra compartimenti piu' che entro". La docstring (righe 24-25) dice però che la copertura di questo IC non è verificata (differenza di due medie su insiemi di pazienti diversi).
- correzione proposta: rendere il messaggio condizionale, per esempio "IC non verificato; vedi 07 per lo stimatore principale".

11. `src/06_clonotype_consistency.py:242`
- gravità: BASSA
- criterio violato: nessuno
- motivo: `def majority(labels) -> str` restituisce None in caso di pareggio.
- correzione proposta: `-> str | None`.

CONTEGGIO: ALTA 1 / MEDIA 1 / BASSA 9

ESITO DELLE VERIFICHE
- a) Superata per l'area b:
  - `majority_label` restituisce None in caso di pareggio, e i pareggi sono esclusi in 06 (90 coppie e 3.943 repliche), 09, 10, 16, 17, 18, 19 e 21.
  - Le sensibilità ci sono in `arbitration_tie_sensitivity` e `error_rate_tie_sensitivity`; mancano in 06 e nei consumatori, senza effetto sulle conclusioni (punto 6).
  - Il ricalcolo indipendente dello schema blood con i pareggi esclusi riproduce 277 / 0,1697 e 117 / 0,0598, e la scomposizione 32 / 4 / 11.
  - Per 13 e 14 (fuori area) ho controllato solo la coerenza del README: la riga 262 contraddice le coordinate reali (punto 2). Il codice di 13 usa la soglia globale 1,5·sd; le 10 finestre attribuite a 9p21 contengono 21,97 Mb.
- b) Non superata per la tabella delle fonti e per il §6 del README (punto 1). Coincidono con le tabelle, circa 130 numeri controllati:
  - tutti i numeri del README §1 "Validazione clonale con il TCR", compresi H1-H4;
  - le sezioni TCR, doppietti, espansione clonale e ipotesi di RISULTATI_FINALI;
  - le righe TCR e H1-H4 della tabella "prima/dopo".
- c) Correttezza biologica, punto per punto:
  - identità CD4/CD8 stabile con eccezioni: superato. Le eccezioni sono dichiarate in README, 09 e 10; `clone_invariant_tcr_check` dà 12 MAIT-like tutti concordanti e 0 iNKT; CD8A 7% nelle discordanti contro CD8αα.
  - marcatori: superato. CD8+ = CD8A>0 o CD8B>0 e CD4+ = CD4>0 senza CD8 in 09 e 10; log1p(CP10k) in 08, 09, 10 e CP10k senza log in 07, dichiarati; NKG7 e GNLY usati solo in modo descrittivo; MAIT TRAV1-2 con TRAJ33/12/20 e iNKT TRAV10-TRAJ18 corretti.
  - riferimento CNV (immunitarie di tutti i pazienti, soglia per paziente): superato nel README (riga 22); il codice è fuori area.
  - loci PDAC e direzioni: superato. Coordinate hg38 in 14: CDKN2A 9:20-24 Mb, TP53 17:0-10, SMAD4 18:45-60, MYC 8:120-145, KRAS 12:20-30, 1p, 6q, 20q; direzioni corrette. Fuori area.
  - malignità e campo di cancerizzazione come ipotesi: superato (README:22).
  - quota degli errori dovuta alle etichette NK e doppietti: superato. 15 su 47 errori, errore senza NK 0,122, d1 0,617 [0,058; 0,849] e d2 0,458 [0,008; 0,653] presentati come non conclusivi; restano da correggere i dettagli dei punti 8 e 9.
  - geni da stress da dissociazione nelle spiegazioni XAI: dichiarati nel README:28 (fuori area).
  - confronto tumore-adiacente su 4 pazienti: superato. IC indicati come non calibrati, 35 ricampionature e 93,3% di repliche > 0, barre grigie nelle figure, nessun verdetto in 08, 09 e 10.
  - ipotesi H1-H4: superato. Direzione dichiarata, circolarità e deriva opposta di H1 dichiarate, risultati negativi formulati come "non esclude un effetto".
- Altre verifiche:
  - le correzioni A4, A5, A6 (download_tcr.sh, 16, 07, 10), A7 (08, 10, 19), B1-B6 e D2 dell'area sono applicate;
  - gli script dell'area e `scpipe.py` hanno mtime anteriore alla riesecuzione del giro 7 (04:44-04:57), tutte le esecuzioni hanno rc = 0 e i .err sono puliti;
  - i pacchetti importati sono tutti in `envs/environment.yml`;
  - un mio script ha scritto per errore `indep_x.csv` nella radice della repo: l'ho spostato subito nello scratchpad e `git status` è pulito.

pytest (`-m "not slow" -p no:cacheprovider`): 108 passati, 0 falliti, 0 saltati (6 deselezionati perché lenti).
ruff (`check src tests revisione`): "All checks passed!".

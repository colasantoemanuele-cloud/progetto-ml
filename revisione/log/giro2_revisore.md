# Giro 2 - rapporto del revisore

Giro 2 (incrementale, `git diff 52ace44 HEAD`): 0 ALTA, 4 MEDIA, 16 BASSA. `pytest` (veloci) passa, `ruff` è pulito. Rigenerato nello scratchpad, `RISULTATI_FINALI.md` è identico byte per byte a quello versionato. I problemi ALTA/MEDIA del giro 1 sono tutti risolti; nessuna correzione ha introdotto crash, file orfani o numeri del README non aggiornati. Non ho modificato nessun file della repo; i miei script di controllo sono in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/reviewer2/` (readme/, gen/, tcr/, ml/).

## MEDIA

**1. `src/scpipe.py:288-292` (`permutation_corr_test`) e `:322-327` (`block_permutation_test`); conseguenza in `src/18_loyo_check.py:150`**
- Gravità: MEDIA
- Criterio: 8
- Motivo: se x o y è costante, rho è NaN. Tutti i confronti `abs(rhos) >= NaN` sono falsi, quindi il p diventa il minimo, 1/(n_perm+1), e il test risulta «significativo» senza avviso.
  - Riprodotto: `permutation_corr_test(np.arange(10), np.zeros(10), 2000)` dà `(nan, 0.0004997...)`.
  - `block_permutation_test` dà rho = NaN e p = p_centered = 0,0005.
  - In 18, con un esito costante, `loyo()` stampa «p < 0.05 in 3/3» e il verdetto, poi va in crash su `out.rho_without.idxmin()` con `KeyError: nan`.
  - Sui dati reali non succede, ma questo caso limite produce un risultato sbagliato senza messaggio.
- Correzione:
  - se `np.isnan(rho_obs)`, restituire p = NaN con un avviso, oppure fermarsi con un messaggio chiaro;
  - in 18, gestire `rho_without` tutto NaN prima di `idxmin`;
  - aggiungere un test in `tests/test_stats.py`.

**2. Colonne a unità cellula/clone non dichiarate nel README**
- Dove: `src/17_clonal_expansion.py:198-212` (righe `clone_naive_permutation` in `clonal_expansion_tests.csv`) e `src/07_robustness.py:234` (`naive_bootstrap`, colonne `lo_naive`/`hi_naive` in `excess_summary_{cc}.csv`)
- Gravità: MEDIA
- Criterio: 2
- Motivo: sono un test e un IC con il clone o la coppia come unità, salvati in tabella. Il codice li etichetta come confronto («unita' sbagliata, solo per confronto»), ma il criterio 2 chiede che sia il README a dichiarare l'unità e il motivo, e il README non li nomina. Lo stesso tipo di problema era il n. 5 del giro 1 (Mann-Whitney per cellula).
  - Valori: eccesso `lo_naive`/`hi_naive` = 0,048; 0,070, contro il cluster 0,031; 0,110.
  - Errore: `clone_naive_permutation` dà p = 0,0030, contro 0,0009 del block test.
- Correzione: una riga nel README (§5 o §6) che dichiari queste colonne come confronto con l'unità sbagliata, oppure toglierle.

**3. `revisione/genera_risultati.py:823-828` (e `:136`); `revisione/RISULTATI_FINALI.md:286-288`**
- Gravità: MEDIA
- Criterio: 4 (stessa lettura del giro 1, n. 13, la cui parte su 18 è stata accettata)
- Motivo: B4 per 18 manca nel file. Per il leave-one-patient-out il generatore riporta solo «10/10» e l'intervallo di ρ.
  - Non legge `loyo_clonal_expansion_perm_sensitivity.csv`.
  - Non legge le colonne `null_mean_without` e `p_centered_without` di `loyo_clonal_expansion_{error,excess}.csv`.
  - Il README:36 cita «p centrato … 5 semi o 100.000 permutazioni (p dell'errore fino a 0,012)», che in `RISULTATI_FINALI.md` non compare.
- Correzione: aggiungere per 18, per ciascun esito, il p centrato massimo e il minimo e il massimo di p fra semi e permutazioni, letti dalla tabella delle sensibilità, con la fonte.

**4. `revisione/genera_risultati.py:7-9` contro le righe 444, 513, 551, 681, 745, 844 e 366-371**
- Gravità: MEDIA
- Criterio: 6 (la docstring contraddice il codice)
- Motivo: la docstring dice che il file «non contiene cifre scritte a mano, tranne i valori della presentazione nella colonna "prima"». Il codice invece scrive a mano:
  - «13 modelli … 20.000 cellule di training per fold»;
  - «Con 5-13 fold»;
  - «con <= 5 pazienti»;
  - «4 pazienti: IC non calibrati»;
  - «n = 10 pazienti»;
  - «≈ 0,53»;
  - il dizionario `h_old` nel testo.
  - Oggi coincidono con le tabelle (`triple_clones_summary.n_patients` = 4, `ni_score_test.n_patients` = 10, `n_train_mean` = 20000), ma con tabelle diverse resterebbero vecchi senza che nessuno se ne accorga.
- Correzione: leggerli dalle tabelle, oppure correggere la docstring dicendo che i parametri del disegno sono scritti a mano.

## BASSA (criterio: nessuno, salvo dove indicato)

5. **`README.md:35`, `:234`, `:235`, `:236`**
   - Criterio 4 solo alla lettera.
   - Motivo: i valori della versione precedente («d2 era 0,790», «0,24», «2 falsi positivi su 100», «0,063 e 0,040») non sono in nessuna tabella. Sono dichiarati come storici.
   - Correzione: spostarli in REGISTRO, oppure scrivere «valore precedente, non riprodotto» senza la cifra.

6. **`README.md:36`**
   - Motivo: «… per entrambi gli esiti …, anche con il p centrato e con 5 semi o 100.000 permutazioni». `perm_sensitivity` (`src/18_loyo_check.py:158`) usa solo `sub["error"]`: per l'eccesso c'è solo il p centrato con 20.000 permutazioni e seme 0 (massimo 0,00045).
   - Correzione: «per l'errore, anche con 5 semi o 100.000 permutazioni (p fino a 0,012)».

7. **`README.md:211`**
   - Motivo: le «35 ricampionature» della riga 34 sono quelle della differenza appaiata (`error_rate_paired_diff.n_distinct_resamples` = 35). La tabella delle fonti le attribuisce a `error_rate_blood`, che contiene 35 ma per la riga Adjacent_normal.
   - Correzione: citare `error_rate_paired_diff`.

8. **`README.md:165` (tempi di 11)**
   - Motivo: per task1a_paired la somma dei log è 476 + 52 s = 8,8 min, che si arrotonda a «circa 9», non a 8. Le somme mescolano tempi di `all_but_lgbm` e di `dense_rerun`. La tabella delle fonti non cita i log.
   - Correzione: «circa 9» e aggiungere la fonte (i log).

9. **`README.md:209`**
   - Motivo: l'IC di `clone_ab` +0,035; +0,121 (riga 33) non ha fonte nella tabella; il valore coincide con `excess_summary_clone_ab` (0,0345; 0,1214).
   - Correzione: «+0,065 [+0,035; +0,121] da `excess_summary_clone_ab`».

10. **`README.md:241`**
    - Motivo: la descrizione di `test_stats.py` non cita il nuovo `test_permutation_corr_type1_n10` (banda [0,02; 0,08], n = 10), cioè proprio il test che rende validato `permutation_corr_test` per H1-H4.
    - Correzione: aggiungerlo.

11. **`README.md:176` e `:178`**
    - Motivo: la struttura dei risultati non elenca `clone_switching_within_clone_b.csv` (scritto da 06) né `doublet_scores.csv` (scritto da 08). Sono esclusi da git; la cache Scrublet è del 10 luglio.
    - Correzione: aggiungerli con «(non versionato)».

12. **`README.md:165` contro `src/13_infercnv.py:346-352`**
    - Motivo: il README dice che le cache sono riusate «solo se coprono esattamente le cellule». 14 lo controlla nei due sensi (righe 197-203); 13 rifiuta solo le cellule mancanti (`score.isna().any()`) e accetta una cache più grande, calcolata con un altro riferimento.
    - Correzione: in 13 lo stesso controllo bidirezionale di 14, oppure distinguere i due script nel README.

13. **`src/11_model_zoo.py:35-36` contro `:311-312`**
    - Motivo: la docstring attribuisce il BLAS a un thread al solo MLP; il commento cita «TruncatedSVD e l'MLP». Secondo REGISTRO cambiano anche knn e svm_rbf (< 0,0004 per fold).
    - Correzione: allineare la docstring al commento.

14. **`src/11_model_zoo.py:422-427`**
    - Motivo: con `--from-table --tag task1` senza `--tasks`, la lista di default viene filtrata a vuota: confronti e figure vengono saltati in silenzio. Il log `from_table_task1.log` mostra che `--tasks` è stato passato.
    - Correzione: con `--from-table`, di default tutti i task della tabella, oppure un avviso.

15. **`revisione/log/11_final/_commit.txt` («8058701») e `REGISTRO.md`**
    - Motivo: `model_zoo_fin_rest.csv` e `model_zoo_fin_v1_lightgbm.csv` non hanno la colonna `converged`, quindi vengono dal codice di 52ace44 (avvio 00:12-00:14). Da lì vengono dummy, naive_bayes, decision_tree, bagging e RF sui tre task e LightGBM su task2_v1.
      - REGISTRO cita sia 52ace44 sia 8058701.
      - I numeri sono equivalenti: per quei modelli le modifiche non cambiano il calcolo.
    - Correzione: registrare il commit che ha prodotto ciascun file `fin_*`.

16. **`src/sim/sim_c_block_permutation.py:16-18`**
    - Motivo: la docstring dice che le prime 100 repliche sono «quelle del disegno … citato nella presentazione (stesso seme)». Sono invece quelle della versione precedente dello script (semi 300+i): il 2/100 della presentazione non è riprodotto (4, 6, 11, 5). README:235 lo dice correttamente.
    - Correzione: allineare la docstring al README.

17. **`src/10_loco.py:711-716`**
    - Motivo:
      - il messaggio «La stima puntuale e' di segno non costante in tutti gli schemi» si contraddice;
      - il ramo `(ta.lo > 0).all()` dà un verdetto anche mescolando schemi con 5 pazienti o meno (IC non calibrati).
    - Correzione: «non ha lo stesso segno in tutti gli schemi»; valutare `lo > 0` solo sugli schemi con `n_patients > SMALL_K`.

18. **`src/08_disambiguate.py:464-467`** (anche `triple_clones.csv` in part1 e `doublet_mixture_boot.csv` in 10)
    - Motivo: con un risultato vuoto escono prima di `to_csv`, quindi una tabella di un'esecuzione precedente resterebbe in `results/` senza avviso.
    - Correzione: scrivere la tabella vuota con colonne fisse, come fanno 09 e 10 negli altri casi.

19. **`tests/test_stats.py:30-32`**
    - Motivo: «usato per H1-H4 in 16 e 19»: manca 17, che lo usa a livello paziente (n = 10).
    - Correzione: aggiungere 17.

20. **`revisione/genera_risultati.py`: formattazione ed esiti in `RISULTATI_FINALI.md`**
    - `:84-87`: il p dell'eccesso, 4,99975e-05, è stampato «0,0000».
    - `:232-237` e `:407-428`: l'esito «cambiato» non distingue una correzione di arrotondamento da un risultato cambiato:
      - 0,851 → 0,850 era già 0,85049;
      - il vecchio 0,53 è `auc_MEAN_|CNV|` (0,533), non la media dei 4 pazienti (0,540);
      - manca MLP task2_hier3 0,903 → 0,904.
    - `:443-447`: la nota dichiara non calibrati solo gli IC con 4 pazienti, ma anche quelli con 5 lo sono (copertura 0,915-0,925).
    - `:541-552`: «8 (12) modelli su 12 con p di Nadeau-Bengio < 0,05» su 13 fold LOPO, senza dire che lì la calibrazione non è verificata.
    - `:492` e `:623`: Δ = −0,00017 stampato «-0,000»; conteggi senza separatore delle migliaia («1114»).
    - Correzione: `p < 0,001` in notazione `1e`; un esito «arrotondamento/trascrizione» separato; aggiungere la riga dell'MLP; «4-5 pazienti»; una nota su NB in LOPO; usare `n_()`.

## Problemi segnalati e respinti

- **Convergenza dell'MLP in `model_zoo_*`:** un sotto-controllo riteneva che `converged` venisse da un addestramento diverso da quello delle metriche. È falso: `revisione/log/11_final/model_zoo_fin_dense.csv` (dense_rerun, dopo la correzione dei thread) ha `converged` = True per l'MLP, e i valori per fold coincidono con `model_zoo_task1.csv`.
- **Tabelle di 07 senza suffisso `{cc}`:** `ambient_discriminant.csv`, `chat_sensitivity_{raw,hc}.csv` e `ambient_vs_real_tcr.png` non hanno il suffisso. Però dipendono solo da `has_tcr` = `clone_b.notna()`, quindi sono identiche per le due definizioni di clone: la sovrascrittura non è un problema.

## Altre verifiche
- **README contro tabelle:** circa 230 numeri verificati, nessuna discrepanza di arrotondamento, nessun file orfano in `results/`.
- **Doppietti:** la formula d2 con s + c è in `10_loco.py:489` e in `09:357`. Ricalcolato da `doublet_physical.csv`: d1 = 0,7226, d2 = 0,6756, uguali a `doublet_mixture.csv`.
- **Ambiente:** tutti gli import sono in `environment.yml` (threadpoolctl, requests e openpyxl compresi); lock con percorso relativo.
- **Comandi del README:** tutti gli argomenti esistono nei rispettivi `--help`.
- **Problemi del giro 1:** ALTA 1 e MEDIA 2-19 tutti risolti; BASSA risolte salvo i respinti già motivati.

## Conteggio
- **ALTA:** 0
- **MEDIA:** 4
- **BASSA:** 16
- **pytest** (`-m "not slow" -p no:cacheprovider`): 100 passati, 0 falliti, 0 saltati, 6 deselezionati (i test lenti), 54 s
- **ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni

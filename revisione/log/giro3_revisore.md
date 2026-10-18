# Giro 3 - rapporto del revisore

Giro 3, revisione incrementale su `git diff 8058701 HEAD`. Non ho modificato nessun file della repo; i miei script di controllo sono in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/reviewer3/`.

**Verifiche fatte**
- **README contro tabelle:** ho controllato circa 240 numeri, compresi tutti quelli delle sezioni 1, 3, 6 e 7. Nessuna discrepanza, anche tenendo conto dell'arrotondamento. Alcuni casi limite controllati con il valore esatto:
  - 0,00185 in tabella è 0,0018499 e si arrotonda correttamente a 0,0018;
  - 7% viene da 0,066;
  - 0,220 viene da 0,2195;
  - 45-62% viene da 0,452 e 0,615 (pazienti con n ≥ 20).
- **RISULTATI_FINALI.md:** l'ho rigenerato nello scratchpad (copia dello script con un altro percorso di uscita). È identico byte per byte a quello versionato. Le colonne lette sono quelle giuste: `p_holm_family_of_4`, `lo_cluster`/`hi_cluster`, `n_distinct_resamples`, `p_centered_without`, B4 di 18 da `loyo_clonal_expansion_perm_sensitivity`.
- **Cache CNV:** `data/interim/cnv_scores.csv` contiene esattamente le 101.219 cellule selezionate da 13, quindi il nuovo controllo nei due sensi non blocca il comando del README.
- **Ambiente:** tutti gli import sono in `envs/environment.yml`; l'unico nuovo è `warnings`, della libreria standard.
- **Problemi ALTA/MEDIA del giro 2:** risolti tutti e quattro (1-4). La correzione della MEDIA 1, però, introduce i problemi M1 e M2 qui sotto, e quella della MEDIA 4 è incompleta (M3).

## MEDIA

**M1. `src/scpipe.py:357-370` (`holm_bonferroni`), usata in `src/16_pseudobulk_patient.py:301-306` e `src/19_ni_score.py:256-262`**
- Gravità: MEDIA (problema nuovo, introdotto dalla correzione del giro 2)
- Criterio violato: 8
- Motivo: ora `permutation_corr_test` (`scpipe.py:291-295`) restituisce p = NaN se la covariata è costante. `holm_bonferroni` non gestisce i NaN:
  - `sorted` con chiavi NaN dà un ordine qualsiasi;
  - `min(k*NaN, 1)` dà NaN;
  - `max(0.0, NaN)` dà 0,0.
- Riprodotto:
  - `holm_bonferroni({'exhaustion_score': nan, 'myeloid_fraction': 0.113, 'caf_fraction': 0.673})` dà `{'exhaustion_score': 0.0, 'myeloid_fraction': 0.226, ...}`;
  - in 16 l'ipotesi viene quindi stampata «SIGNIFICATIVO» e salvata con `p_holm` = 0, e il fattore delle altre ipotesi scende da 3 a 2;
  - nella famiglia di 4 di 19 l'esito dipende dall'ordine del dizionario.
- È lo stesso tipo di errore della MEDIA 1 del giro 2 (p spurio senza avviso), spostato nella correzione di Holm. Sui dati reali non succede.
- Correzione: in `holm_bonferroni`, escludere i p NaN dalla famiglia (adjusted = NaN, con un avviso) oppure fermarsi con un messaggio. Aggiungere il caso a `tests/test_stats.py`, dove oggi c'è solo «Holm monotono e ≤ 1».

**M2. `src/18_loyo_check.py:123, 141-151`**
- Gravità: MEDIA (residuo della MEDIA 1 del giro 2)
- Criterio violato: 8
- Motivo: con un esito costante il crash non c'è più, ma il messaggio resta sbagliato. Ho eseguito `main()` su dati simulati (4 pazienti, `error` = 0), con FIG e TAB reindirizzati nello scratchpad. Lo script stampa, in quest'ordine:
  - «rho=+nan p=nan [non piu' significativo senza questo paziente]» per ogni esclusione;
  - «[!!] p non definito in 4 esclusioni»;
  - «VERDETTO: il segnale dipende da pochi pazienti specifici. Non generalizzabile».

  Il verdetto è falso: non c'è nessun segnale, perché il ρ è NaN già su tutti i pazienti.
- Correzione: se `full["rho"]` è NaN, stampare un solo messaggio («esito costante, nessun test») e saltare verdetto e flag. Contare `still_significant` solo sulle esclusioni con p definito.

**M3. `revisione/genera_risultati.py:7-12` (docstring) contro le righe 48-52 (`ESITO`), 439, 541 e 775**
- Gravità: MEDIA
- Criterio violato: 6
- Motivo: la docstring dice che sono scritti a mano «solo» i valori della presentazione, `h_old` e SMALL_K. Sono invece scritti a mano anche:
  - il dizionario `ESITO`, che sostituisce il confronto calcolato (`same`) con verdetti fissi, «arrotondamento» e «uguale»: se le tabelle cambiano, il verdetto resta;
  - «13 modelli» (541);
  - «dei 4 pazienti» (51, 439);
  - «IC non calibrati» (775), stampato qualunque sia `tri.n_patients`.
- È la stessa classe di problema della MEDIA 4 del giro 2, corretta solo in parte.
- Correzione: applicare `ESITO` solo quando i valori coincidono con quelli attesi, oppure dichiararlo nella docstring. Leggere il numero di modelli e di pazienti dalle tabelle. Scrivere «IC non calibrati» solo se `n_patients <= SMALL_K`.

**M4. `revisione/genera_risultati.py:51, 439-443, 471-472`; in `RISULTATI_FINALI.md` la riga «AUC loci PDAC entro paziente … 0,53 → 0,533 / 0,540, uguale (il valore precedente era la mediana)»**
- Gravità: MEDIA
- Criteri violati: 5 e 4
- Motivo: l'affermazione «il valore precedente era la mediana» non ha riscontro.
  - Il codice che ha prodotto lo 0,53 (`git show 4d42a8f:src/14_cnv_loci.py`, righe 274-286) stampava solo l'AUC di ogni paziente, senza media né mediana.
  - Il README pubblicato (`main:README.md:15`) scriveva solo «≈ 0,53».
  - Il README attuale riporta «≈ 0,54», cioè la media (0,540): il numero pubblicato è cambiato, ma la tabella dice «uguale».
  - La mediana 0,533 è calcolata dal generatore (`auc_pat.median()`) e non è in nessuna tabella di `results/tables`. Coincide per caso con `auc_MEAN_|CNV|` della riga `_mean`, che è un'altra grandezza.
- Correzione:
  - esito «cambiato (0,53 → 0,54, media dei 4 pazienti; la statistica dello 0,53 non è ricostruibile)», senza attribuirlo alla mediana;
  - oppure far scrivere la mediana a 14 in `cnv_loci_auc_by_patient`, se si vuole citarla.

**M5. `README.md:35`, «compatibili con una quota di doppietti T+NK fra il 68% e il 72%» (in grassetto)**
- Gravità: MEDIA
- Criterio violato: 5
- Motivo: l'intervallo 68-72% è solo la distanza fra due stime puntuali non indipendenti (d1 0,723 e d2 0,676, `doublet_mixture`, riga median/all). Non è un intervallo d'incertezza. Gli IC vanno da 0,556 a 0,848; per compartimento le stime vanno da 0,42 (d2, Tumor) a 0,85 (PBMC). Il grassetto suggerisce una precisione di ±2 punti che i dati non hanno; le cautele arrivano solo dopo, nello stesso paragrafo.
- Correzione: «compatibili con una maggioranza di doppietti (d1 = 0,72, IC 0,56-0,85; d2 = 0,68, IC 0,56-0,78; da 0,42 a 0,85 secondo il compartimento)», senza l'intervallo 68-72%.

## BASSA (criterio violato: nessuno)

1. **`revisione/genera_risultati.py:865`**
   - Motivo: stampa «100.000 permutazioni (1 semi)». Inoltre la riga 855 apre l'elenco con «\n- », quindi in `RISULTATI_FINALI.md` il primo punto del leave-one-patient-out è separato dagli altri da una riga vuota.
   - Correzione: «1 seme» / «5 semi», e togliere il «\n» iniziale.
2. **`revisione/genera_risultati.py:681, 815, 902, 929, 931, 1003`**
   - Motivo: restano conteggi senza separatore delle migliaia («1519», «2354», «1000»). Il punto 20 del giro 2 è corretto solo in parte.
   - Correzione: usare `n_()`.
3. **`README.md:34` contro `src/10_loco.py:717-724`**
   - Motivo: il README conclude che «l'asimmetria dipende dallo schema di riferimento». Ma tutti e tre i contrasti hanno 5 pazienti o meno (`error_rate_paired_diff`: 5, 5, 4), e per questo caso 10 stampa ora «nessun verdetto».
   - Correzione: «le stime puntuali variano con lo schema (+0,027 LOCO, +0,035 naive, +0,080 solo sangue); con ≤ 5 pazienti nessun IC è calibrato».
4. **`README.md:41` e `:233`**
   - Motivo: il README usa Nadeau-Bengio su fold LOPO di dimensione diversa (p ≥ 0,080). La nota che la calibrazione è verificata solo con StratifiedKFold a k = 3 e 5, presente in `genera_risultati.py:579-581`, nel README manca.
   - Correzione: aggiungere la stessa nota alla sezione 6 o ai Limiti.
5. **`src/11_model_zoo.py:395-396`**
   - Motivo: con `--from-table` restano obbligatori `--h5ad` e `--splits`, che non vengono letti.
   - Correzione: `required=False`, con il controllo solo senza `--from-table`.
6. **`src/sim/sim_c_block_permutation.py:16-19` contro `:99`**
   - Motivo: «semi 300 + i» è ambiguo. i è l'indice del livello di confondente (`seed=300 + i`), non della replica.
   - Correzione: «seme 300 + indice del livello».
7. **`src/08_disambiguate.py:193-195, 469` e `src/10_loco.py:542-547`**
   - Motivo: con risultato vuoto ora le tabelle vengono cancellate, ma le docstring (sezione Output) non lo dicono.
   - Correzione: una riga nella sezione Output.
8. **`src/sim/sim_a_calibration.py:150`**
   - Motivo: il cluster bootstrap è calibrato con `n_boot=1000`, mentre le analisi usano 2.000 repliche (README:93). L'effetto sulla copertura è trascurabile, ma il disegno non è identico.
   - Correzione: dichiararlo nella docstring, oppure usare 2.000.

## Conteggio
- **ALTA:** 0
- **MEDIA:** 5 (M1 e M2 introdotti dalla correzione della MEDIA 1 del giro 2; M3 correzione incompleta della MEDIA 4 del giro 2; M4 introdotto nel giro 2; M5 non in questo diff)
- **BASSA:** 8
- **pytest** (`-m "not slow" -p no:cacheprovider`): 101 passati, 0 falliti, 0 saltati, 6 deselezionati (i test lenti), 56 s
- **ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni

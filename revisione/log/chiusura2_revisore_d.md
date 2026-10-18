# Giro di chiusura 2 - area (d) README, ambiente, docstring

## Area (d), giro di chiusura 2: problemi trovati

**1. `revisione/genera_risultati.py:648`**
- gravita': ALTA (crash; succede solo se manca una tabella, quindi si puo' declassare)
- criterio violato: 8
- motivo: è lo stesso tipo di crash corretto nel giro precedente (d1), ma su un'altra tabella. La riga `f"... in tutti i {fs_t1.n_folds.max():.0f} fold LOPO ..."` sta dentro `if not core_t1.empty`, quindi non controlla `fs_t1`.
  - Ho tolto le tabelle una alla volta, in una copia nella scratchpad fatta di link simbolici. Se manca solo `feature_selection_summary_task1a_paired.csv`, lo script termina con `AttributeError: 'DataFrame' object has no attribute 'n_folds'`, anche senza `--strict`. La docstring invece promette di registrare la tabella mancante in `MISSING`. Con tutte le altre 93 tabelle tolte una alla volta lo script non va in crash.
  - Stesso comportamento se manca una colonna. Provato togliendo `n_folds` da `model_zoo_summary_task1`, `rep` da `sim_d_xai_stability`, `still_significant` da `loyo_clonal_expansion_error` e `d1_lo` da `doublet_mixture`: si ha un traceback invece di una voce in `MISSING`. Il codice accede alle colonne come attributi (righe 590, 213-215, 882, 817).
- correzione proposta: `fs_t1.n_folds.max() if not fs_t1.empty else np.nan`. In generale, prima di usare le colonne come attributi, verificare con `col in df` e registrare in `MISSING` quelle assenti.

**2. `src/13_infercnv.py:21-22`**
- gravita': MEDIA
- criterio violato: 6 (un commento contraddice i dati attuali)
- motivo: la docstring dice "In `malignancy_labels.csv` risultano 'maligne' anche il 45-62% delle duttali di `Adjacent_normal` nei pazienti con almeno 20 duttali adiacenti".
  - Dopo la riesecuzione della chiusura 1 (senza il cromosoma MT) la tabella dà PA06 0,577, PA08 0,410 e PA10 0,530, cioè 41-58%.
  - Il README (riga 22) e il REGISTRO ("45-62% -> 41-58%") sono già aggiornati. È l'unica occorrenza rimasta del vecchio intervallo in `src/`, `tests/`, README e generatore.
- correzione proposta: "41-58%" (oppure non citare il numero e rimandare alla tabella).

**3. `src/13_infercnv.py:26-28`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: "Ensembl release 116 ..., la release con cui sono state calcolate le cache dell'11/07/2026 (le coordinate delle finestre di `cnv_matrix.npz` coincidono con quelle dei geni di Ensembl 116)". Nella chiusura 1 le cache sono state ricalcolate senza cache (REGISTRO, giro 4). La motivazione si riferisce quindi a cache che non esistono più.
- correzione proposta: "release 116, la stessa usata per le cache precedenti; le cache attuali sono ricalcolate con questa tabella".

**4. `README.md:17` (riga 3 della tabella)**
- gravita': BASSA
- criterio violato: nessuno
- motivo:
  - Il testo dice "entro paziente i 4 pazienti con entrambi i tessuti". In realtà i pazienti con entrambi i tessuti sono 5: in `malignancy_labels`, PA04 ha 152 duttali tumorali e 1 adiacente. I 4 pazienti sono quelli valutabili, cioè con almeno 10 duttali per tessuto (14, riga 367), come la nota a riga 22 dice correttamente.
  - Nella colonna "Stima che rispetta la struttura di gruppo" compare anche "lo score CNV medio ha AUC pooled 0,531", che è una stima pooled.
- correzione proposta: "i 4 pazienti valutabili (almeno 10 duttali per tessuto)" e spostare lo 0,531 nella colonna "pratica corrente" o nelle note.

**5. `README.md:178` (Struttura dei risultati, riga 06)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: `clone_switching_within_clone_b` è descritto come "una riga per coppia entro compartimento". La tabella ha invece una riga per replica split-half di ogni (clone, compartimento): colonne `patient, clone, tissue, lab_a, lab_b, switch`, 45.120 righe, pari a `n_reps_within`. La docstring di 06 dice correttamente "repliche split-half".
- correzione proposta: "una riga per replica split-half di ogni (clone, compartimento)".

**6. `README.md:131` e `README.md:161` (ordine dei comandi)**
- gravita': BASSA
- criterio violato: nessuno (le tabelle non cambiano)
- motivo: 20 legge `results/tables/sim_d_xai_stability.csv`, se presente, per la "LETTURA" finale stampata (20, righe 328-366). Nel README pero' 20 viene eseguito al passo 3, mentre `sim_d` è al passo 7. Partendo da zero, 20 stampa quindi "sim_d_xai_stability.csv assente ... nessuna conclusione" invece del confronto riportato nel log. `xai_summary` e la figura non dipendono da `sim_d`.
- correzione proposta: spostare `python src/sim/sim_d_xai_stability.py` prima di 20, oppure annotare che la lettura di 20 richiede `sim_d`.

**7. `revisione/genera_risultati.py:11-12` (contro la riga 455)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice che l'esito "arrotondamento" è calcolato quando il valore esatto della tabella dista meno di 0,001 dal valore precedente. Il codice applica la regola solo alle righe che passano il valore esatto (`raw`), cioè solo alla logistica L2 di hier3.
  - La riga dell'MLP (0,903 contro 0,90376, distanza 0,0008) soddisfarebbe la regola dichiarata, ma è etichettata "cambiato".
  - L'etichetta è giusta nel merito: il valore è cambiato davvero, da 0,9032 a 0,9038, per la correzione sui thread (REGISTRO, riga 180). È la docstring a non descrivere il codice.
- correzione proposta: "l'esito 'arrotondamento' si applica solo alle righe in cui il valore precedente è la stessa stima arrotondata diversamente (indicate nel codice)".

**8. `revisione/genera_risultati.py:381` e `:391`**
- gravita': BASSA
- criterio violato: nessuno (il valore oggi coincide)
- motivo: "su 1.000" è scritto a mano, mentre il numero di repliche è nella colonna `reps` di `sim_c_block_permutation`. La docstring (righe 9-11) dice che numeri di pazienti, fold e repliche vengono letti dalle tabelle. Lo stesso vale per l'etichetta "geni comuni ai 5 fold (task1a, k = 50)" alla riga 311.
- correzione proposta: `f" su {n_(sc_.reps.iloc[0])}"`; per l'etichetta della riga 311 usare `k_core` e `fs_t1.n_folds.max()`.

**9. `README.md:238` (sim b)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: "il valore 0,24 citato nella versione precedente non e' riprodotto esattamente" fa pensare a una riproduzione approssimata. In `sim_b_reference_schemes` (loco, PBMC) i valori sono 0,126 con artefatto moderato e 0,385 con artefatto forte: nessuno dei due è vicino a 0,24.
- correzione proposta: "non e' riprodotto (0,126 e 0,385 con le due intensita' scelte a priori)".

## Verifiche senza problemi
- Ho confrontato con `results/tables/` ogni numero del README (sezioni 1-7 e tabella "Da dove viene ogni numero") e di `RISULTATI_FINALI.md`: tutti coincidono entro l'arrotondamento.
- `genera_risultati.py --strict`, rieseguito in una copia nella scratchpad, ricrea `RISULTATI_FINALI.md` identico byte per byte, con 0 valori mancanti.
- I comandi della sezione 4 del README corrispondono a `--help` di ogni script (tutti rispondono) e l'ordine è eseguibile.
- Ogni pacchetto di terze parti importato in `src/`, `tests/` e `revisione/` è in `envs/environment.yml`.
- `envs/` è invariato da 8058701. Le versioni chiave del lock coincidono con quelle dell'ambiente.
- `ruff format --check` passa (35 file), quindi l'hook di `.pre-commit-config.yaml` è soddisfatto. La versione dell'hook (ruff 0.15.20) è quella dell'ambiente.
- Docstring di tutti gli script di `src/` e `src/sim/`: scopo, input, output (confrontati con ogni `to_csv`, `savefig`, `savez` e `write_text`), argomenti e uso sono completi. Le righe "Uso" hanno `\\`. Le correzioni D1 e D2 delle docstring e D3 del README sono applicate. I problemi 2-14 del mio rapporto di chiusura 1 sono corretti.
- `git diff`:
  - 00-03 e `download_tcr.sh` da 61fb36b: cambiano solo argparse, `help` e docstring.
  - 05 da 52ace44: cambiano solo il controllo sul training con una sola classe, `nargs="+"` e `help`.
  - 11 da 8058701 e 12 da 52ace44: cambiano solo `main` (controlli, `--from-table`, etichetta dell'asse).
  - Le funzioni e classi di `scpipe.py` usate da 05, 11 e 12 hanno AST identico a 52ace44 (cambiano solo `block_permutation_test`, `permutation_corr_test`, `holm_bonferroni` e `small_k_diagnostics`, usate da script rieseguiti).

## Conteggio
- ALTA: 1 / MEDIA: 1 / BASSA: 7
- pytest (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 lenti deselezionati)
- ruff (`check src tests revisione`): All checks passed!

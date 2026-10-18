# Giro di chiusura 3 (giro 6) - area (a) ML

**1. `src/20_xai_per_patient.py:166, 205-210, 257-264`**
- gravità: MEDIA
- criterio violato: 8 (anche 6: l'help dice «un fold per paziente», ma `choices` accetta schemi che non lo sono)
- motivo: `--scheme` accetta `fold_random_5` (la scelta è stata aggiunta nella chiusura 2). Ogni fold viene però etichettato con il paziente più frequente (`pd.Series(patients[folds == f]).mode()[0]`) e le matrici J e R sono indicizzate da queste etichette.
  - Con `fold_random_5` su task1a_paired tutte e 5 le etichette valgono PA04 (verificato su `splits_sc.csv`: modalità `['PA04', 'PA04', 'PA04', 'PA04', 'PA04']`).
  - `J.loc[na, nb] = ...` con etichette duplicate scrive allora l'intera matrice, diagonale compresa. Riprodotto: tutte le celle prendono il valore dell'ultima coppia, e `j_mean` diventa il Jaccard dell'ultima coppia.
  - Lo script stampa «escluso PA04» per ogni fold, scrive `xai_jaccard/xai_spearman/xai_summary` sbagliati e non dà nessun messaggio.
  - Con `fold_group_5` le etichette nominano un solo paziente anche quando il fold ne contiene più di uno.
  - I numeri pubblicati non cambiano: sono tutti con `fold_lopo`.
- correzione proposta: `choices=["fold_lopo"]`, oppure `SystemExit` se un fold contiene più di un paziente; in alternativa usare l'indice del fold come etichetta.

**2. `README.md:16` (riga 2, colonna «Origine dell'artefatto»)**
- gravità: MEDIA
- criterio violato: 5
- motivo: il README scrive «Confondimento paziente ↔ tessuto: ... il modello sfrutta la firma del paziente e la composizione cellulare». È un meccanismo che nessuna tabella misura, e i dati non lo sostengono.
  - In task1a_paired il disegno è incrociato: ogni paziente ha entrambi i tessuti, quasi bilanciati (Adjacent_normal 0,506-0,613 per paziente: PA04 0,613, PA06 0,602, PA08 0,513, PA09 0,570, PA10 0,506). Quindi la sola «firma del paziente» quasi non predice il tessuto.
  - La composizione cellulare per tessuto è un segnale presente anche in LOPO, quindi non spiega il divario.
  - Quello che lo split casuale condivide fra training e test è la libreria: in `obs` c'è esattamente 1 libreria per paziente × tessuto e ogni libreria ha un solo paziente e un solo tessuto (verificato). Il modello può quindi sfruttare caratteristiche della libreria o del campione che non si trasferiscono a un paziente nuovo.
- correzione proposta: «Struttura per campione: ogni libreria corrisponde a un paziente × tessuto; lo split casuale mette in training e test cellule della stessa libreria, e il modello può sfruttarne le caratteristiche, che non si trasferiscono a un paziente nuovo (meccanismo non misurato)».

**3. `src/04_splits.py:174-188, 291`**
- gravità: BASSA
- criterio violato: 8 (marginale: caso non raggiungibile con i dati reali)
- motivo: il controllo aggiunto (`len(df) < N_FOLDS`, almeno 2 pazienti e 2 classi) non basta. `StratifiedKFold(n_splits=5)` fallisce quando tutte le classi hanno meno di 5 cellule. Riprodotto: 7 cellule, 3 pazienti, 2 classi (4 + 3) superano il controllo e `assign_folds` va in crash con `ValueError: n_splits=5 cannot be greater than the number of members in each class.` senza il nome del task.
- correzione proposta: saltare il task con un messaggio se `df.y.value_counts().max() < N_FOLDS`, oppure intercettare il `ValueError` e citare il task.

**4. `src/20_xai_per_patient.py:223-237, 254, 279-281, 379`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: quattro difetti di messaggi e testo.
  - Con un `--task` assente da `splits` lo script carica tutto l'h5ad e poi esce con «servono almeno 2 fold in fold_lopo» invece di «task assente». 12 invece dice «task ... assente».
  - Intestazioni fisse «DEI 5 FOLD», «unione delle 5 spiegazioni» e «TUTTE e 5»: sono sbagliate con un altro task (per esempio task2_v1, 13 fold LOPO).
  - La frase stampata «Il controllo condivide piu' cellule fra repliche di due training LOPO» è sgrammaticata. Il calcolo atteso conferma il contenuto: circa 5.912 cellule in comune fra due repliche del controllo contro circa 5.539 fra due training LOPO.
- correzione proposta: controllo `sub.empty` prima di `get_indexer` con il messaggio «task assente»; usare `len(fs)` nelle intestazioni; scrivere «le repliche del controllo condividono fra loro piu' cellule di due training LOPO».

**5. `README.md:213` (tabella delle fonti, riga «p ≥ 0,080, p = 0,072»)**
- gravità: BASSA
- criterio violato: nessuno (il numero è in tabella)
- motivo: il README a riga 26 cita «Δ = −0,229» del dummy su task1a_paired. Il valore viene da `model_comparison_task1a_paired` (`delta_vs_best` = −0,2287), ma la riga delle fonti nomina solo `p_nadeau_bengio`.
- correzione proposta: «Δ = −0,229, p ≥ 0,080, p = 0,072 | `model_comparison_task1a_paired` (`delta_vs_best`, `p_nadeau_bengio`)».

**Verifiche senza problemi:**
- **Numeri.** Tutti i numeri della sezione 1, di «Confronto dei modelli», delle righe delle fonti e delle sezioni Leakage, 13 modelli, Selezione delle feature e Spiegazioni di RISULTATI_FINALI coincidono con `results/tables`. Controllati, tra gli altri:
  - 05: 0,954 / 0,964, Δ −0,009; 0,875 ± 0,002, 0,729 ± 0,100, Δ −0,146, rapporto 50,3; tutta la tabella Leakage.
  - 11: 0,936-0,943 e 0,962; 0,850 (0,850486) e i Δ da +7,2 a +1,5; 0,496-0,725; dummy Δ −0,229, p 0,072; p ≥ 0,080 (0,0802); 8/12, 12/12, 0/12; la convergenza.
  - 12: tutti i valori di F1 e Jaccard; 524-583, 177-268, 1616-1984; geni comuni 18/19/19/16/4.
  - 20: 0,300 / 0,459 / +0,159 / 0,234, 0,220-0,408, 9 geni; sim_d +0,143 (+0,116; +0,181) e +0,030 (−0,013; +0,052).
- **Esecuzione finale di 11.** L'unione dei file `model_zoo_fin_*` secondo `_commit.txt` ricostruisce esattamente `model_zoo_sc` + `model_zoo_task1`: 403 righe, tutte presenti in entrambi, differenza 0 su ogni metrica, `converged` identica. Ho anche verificato su dati sintetici che LightGBM dà predizioni identiche con 1, 4 e 12 thread.
- **File non rieseguiti in Fase 3.** Confronto dell'AST senza docstring:
  - `scpipe` da 52ace44 cambia solo in `block_permutation_test`, `holm_bonferroni`, `permutation_corr_test`, `small_k_diagnostics` e `SMALL_K`, nessuno usato da 05, 11 e 12;
  - 05 cambia solo in `run_task` (controllo di una sola classe) e `main`;
  - 12 cambia solo in `main` e in `sel_l1`, dove gli avvisi sono riemessi e il calcolo non cambia;
  - 11 da 8058701 cambia solo in `compare` (stampa) e `main`;
  - 04 cambia solo nei controlli e `splits_sc.csv` è stato ricreato identico.
  - Riepilogo, confronti e figure di 11 e le tabelle di 20 sono stati rigenerati nella chiusura 2 dopo l'ultima modifica dei sorgenti; i messaggi nei log coincidono con il codice attuale. In `results/` non ci sono file orfani.
- **Correzioni dell'elenco per la mia area: tutte applicate.** A3, A6 (scpipe, 05, 20, 11), A7 (20), D2 (11:21, 11:207, 11:316-317, 20:279, geni di 20, 04:27, LogCPM in scpipe e 05), D3 (Jaccard su fold_group_5, sottocampionamento uniforme, +7 solo LightGBM, comandi di 05, 11 e 12). Applicate anche le correzioni BASSA della chiusura 2: `choices`, avvisi riemessi in 12, `--update`, p minimo di Wilcoxon, caso NaN di sim_d, sim_d prima di 20.
- **Leakage e unità statistica: nessun problema.** LogCPM, HVG, SVD, scaler, prefiltro e selettori sono stimati sul solo training del fold, e le spiegazioni di 20 solo sul training. I confronti sono fra fold per paziente. Nadeau-Bengio non è presentato come validato per fold LOPO.
- **Ambiente e `--help`.** Tutti i pacchetti importati sono in `envs/environment.yml`, `threadpoolctl` compreso. `--help` funziona per 04, 05, 11, 12 e 20.

**Conteggio:** ALTA 0 / MEDIA 2 / BASSA 3

**pytest** (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 lenti esclusi dalla selezione).
**ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni.

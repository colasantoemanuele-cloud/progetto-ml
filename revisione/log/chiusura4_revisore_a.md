# Giro 7 (chiusura 4) - area (a) ML

**1. `src/20_xai_per_patient.py:274-281, 323-329, 417-418` (con `:46-47` del docstring)**
- gravità: BASSA
- criterio violato: nessuno (marginale 6)
- motivo: righe e colonne di `xai_jaccard_task1a_paired.csv` e `xai_spearman_task1a_paired.csv`, e gli assi della heatmap, portano il nome del paziente **escluso** senza dirlo. Il docstring parla solo di «matrice di Jaccard fra fold». Nella sezione «geni esclusivi di un solo fold (spiegazione paziente-specifica)» il log stampa «PA04: 16/50 geni non condivisi...». Quei 16 geni vengono dal modello addestrato **senza** PA04, ma l'etichetta fa pensare a geni propri di PA04. Il legame fold-paziente escluso compare solo più sopra nel log («fold 0 (escluso PA04)»). I numeri del README non cambiano.
- correzione proposta: etichette «senza PA04» (o «−PA04») nelle tabelle e negli assi della figura; nel docstring scrivere «indice = paziente escluso dal training»; al posto di «spiegazione paziente-specifica» scrivere «geni presenti solo nella spiegazione addestrata senza quel paziente».

**2. `src/04_splits.py:21`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: il docstring dice «Platelet (anucleate, quasi solo nel sangue)». In `obs` le 1.587 Platelet sono tutte PBMC: 0 in Tumor e 0 in Adjacent_normal (crosstab `all_celltype` × `tissue`, letto in modalità backed). La frase sulle Endocrine è invece corretta: presenti in 12 pazienti, 2 con almeno 10 cellule (verificato).
- correzione proposta: «Platelet (anucleate, solo nel sangue: 1.587 cellule, tutte PBMC)».

**3. `src/05_baseline.py:239, 248`**
- gravità: BASSA
- criterio violato: nessuno (marginale 5)
- motivo: nelle figure `gap_sc_{task}` le legende dicono «split casuale (leakage)» e «split per paziente (onesto)». Per `task2_v1` il README (riga 1) dice che la stima per paziente resta circolare. «Onesto» va oltre ciò che lo split garantisce: l'assenza di leakage di paziente.
- correzione proposta: «split per paziente (nessun paziente in train e test)».

**4. `src/20_xai_per_patient.py:203-204, 292`**
- gravità: BASSA
- criterio violato: nessuno
- motivo: nel giro 6 le intestazioni stampate sono passate a `len(fs)`. Restano però fissi «sui 5 fold LOPO» nel docstring di `main` e «tutti i 5 fold» nel commento a riga 292.
- correzione proposta: «sui fold LOPO» e «tutti i fold».

**5. `src/11_model_zoo.py:69, 443, 511`**
- gravità: BASSA
- criterio violato: nessuno (marginale 6)
- motivo: docstring e `--help` dicono che `--from-table` «non addestra: rigenera riepilogo, confronti e figure». In realtà riscrive anche `model_zoo_{tag}.csv` (riga 511; il log `revisione/log/11_from_table_task1.log` stampa «[ok] scritto results/tables/model_zoo_task1.csv»). Il contenuto non cambia: ho rieseguito `main --from-table` con TAB/FIG in scratch, e summary e i 3 `model_comparison` coincidono con `results/tables` (differenza massima 0,0).
- correzione proposta: scrivere «riscrive anche model_zoo_{tag}.csv, invariato», oppure saltare `to_csv` con `--from-table`.

**6. `README.md:42` (frase finale di «Spiegazioni per paziente») e `src/20_xai_per_patient.py:32-35`**
- gravità: BASSA
- criterio violato: nessuno (testo prescritto dalla correzione D2)
- motivo: CD69, DUSP4 e CREM sono presentati come «geni di risposta immediata associati allo stress di dissociazione». Della firma di dissociazione standard (van den Brink 2017) fanno parte FOS, BTG2 e HSPH1. CD69, DUSP4 e CREM sono geni precoci, ma marcano anche attivazione, residenza tissutale (CD69) ed esaurimento (DUSP4) delle T infiltranti: cioè biologia che può differire davvero fra tumore e adiacente. La riserva «non necessariamente biologia tumore/normale» c'è già, quindi non è un'affermazione falsa. GAPDH e VIM, fra i 9 geni comuni, non sono commentati.
- correzione proposta: «FOS, BTG2 e HSPH1 appartengono alla firma di dissociazione; CD69, DUSP4 e CREM sono geni di risposta immediata che marcano anche attivazione o residenza tissutale delle T».

**7. `src/20_xai_per_patient.py:147` (LIMITE NOTO, nuovo)**
- gravità: LIMITE NOTO
- criterio violato: nessuno sul dataset reale
- motivo: i top-50 si prendono con `np.argsort(-coef)[:TOP_K]`. Se una logistica L1 avesse meno di 50 coefficienti non nulli, i geni mancanti sarebbero scelti fra i pareggi a zero in ordine arbitrario. È lo stesso caso che 12 dichiara (`n_nonzero`, `k_gt_nonzero`), ma 20 non lo registra e non avvisa.
- verifica sul dataset reale: il caso non si presenta.
  - Fold LOPO: da 3.887 a 4.380 coefficienti non nulli (`revisione/log/20_xai_per_patient.log`).
  - Repliche del controllo: da 4.261 a 4.566 (le ho ricalcolate leggendo solo le righe di task1a_paired in modalità backed).
  - Lo stesso ricalcolo riproduce esattamente Jaccard LOPO 0,30028 e controllo 0,45893.
- correzione proposta: aggiungere la riga al README §7 e a `revisione/limiti_noti.md`, oppure un avviso se `(coef > 0).sum() < TOP_K`.

**Verifiche senza problemi**
- **Numeri.** Tutti i numeri della mia area nel README coincidono con `results/tables`: sezione 1, «Confronto dei modelli», righe della tabella delle fonti, tempi di 11 (40/16/18/9 min dalla somma dei tempi per modello nei log). Coincidono anche le sezioni Dataset, Leakage, 13 modelli, Selezione delle feature e Spiegazioni di RISULTATI_FINALI (8/12, 12/12, 0/12; 524-583; 177-268; 1616-1984; geni comuni 18/19/19/16/4; 0,300/0,459/+0,159/0,234; sim_d).
- **Esecuzione finale di 11.** L'unione dei `model_zoo_fin_*` secondo `_commit.txt` ricostruisce le 403 righe di `model_zoo_sc` + `model_zoo_task1` con differenza 0 su tutte le metriche. I dati per le frasi sui thread nel docstring (MLP fino a 0,014; knn e svm_rbf meno di 0,001) sono confermati: MLP 0,01396, knn 0,00012, svm_rbf 0,00038.
- **File non rieseguiti in Fase 3.** Confronto dell'AST senza docstring:
  - 05 e 12 rispetto a 52ace44 cambiano solo in `run_task`/`main` (controlli) e in `sel_l1` (avvisi riemessi);
  - 11 rispetto a 8058701 cambia solo in `compare` (stampa) e `main`;
  - in scpipe cambiano solo funzioni non usate da 05, 11 e 12;
  - `splits_sc.csv` attuale è identico a quello ricalcolato con `04_splits.py` di 52ace44 (446.773 righe).
  - Il summary di 12 è ricalcolabile esattamente dalla tabella per fold, e i log di Fase 1 coincidono con le tabelle. Nessun file orfano in `results/` per 04, 05, 11, 12 e 20.
- **Correzioni dell'elenco per la mia area: tutte applicate**, comprese quelle del giro 6 (20 solo `fold_lopo` con controllo un paziente per fold; 04 con la classe più numerosa < 5; intestazioni di 20; fonti di Δ = −0,229).
- **Leakage e unità statistica.** Nessun leakage: LogCPM, HVG, SVD, scaler, prefiltro e selettori sono stimati sul solo training; le spiegazioni di 20 solo sul training. Unità corretta: Nadeau-Bengio è dichiarato non verificato per fold LOPO, i p sono dichiarati grezzi e Wilcoxon ha il minimo dichiarato (verificato su scipy 1.17: anche con pareggi resta 0,0625).
- **Docstring, ambiente, `--help`.** Docstring complete e coerenti con il codice; tutti i pacchetti importati sono in `envs/environment.yml`; `--help` funziona per 04, 05, 11, 12 e 20.

**Conteggio:** ALTA 0 / MEDIA 0 / BASSA 6 / LIMITI NOTI nuovi 1

**Controlli biologici**
- **Identità CD4/CD8 stabile nel clone:** non pertinente alla mia area (09/10). In `task2_hier3` si usano le etichette degli autori e non si fanno affermazioni sull'identità clonale. OK.
- **Riferimento dell'inferenza CNV:** non pertinente (13/14).
- **Direzione dei loci PDAC:** non pertinente (14).
- **Geni di stress da dissociazione fra i geni comuni XAI:** presenti e segnalati nel README, nel docstring e nell'output di 20, con la riserva «non necessariamente». OK, con la sfumatura del punto 6.
- **Confronto tumore-adiacente:** nella mia area `task1a_paired` ha 5 pazienti in LOPO. Il README dichiara la potenza bassa, nessun intervallo è presentato come affidabile e la riga 2 dichiara il meccanismo «non misurato». OK.
- **Errori nel tumore dovuti alle NK, in gran parte doppietti:** non pertinente (09/10).

**pytest** (`-m "not slow" -p no:cacheprovider`): 103 passati, 0 falliti, 0 saltati (6 lenti esclusi dalla selezione).
**ruff** (`check src tests revisione`): «All checks passed!», 0 segnalazioni.

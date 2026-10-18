# Giro di chiusura 2 - area (c) CNV, simulazioni, test

**Problemi trovati, area (c): 13, 14, src/sim, tests, pyproject**

1. `src/13_infercnv.py:21-22`
- gravita': MEDIA
- criterio violato: 6 (il commento contraddice la tabella prodotta dal codice attuale); tocca anche il 4.
- motivo: la docstring dice ancora che in `malignancy_labels.csv` risultano "maligne" il 45-62% delle duttali di `Adjacent_normal` (pazienti con almeno 20 duttali adiacenti). Sono i valori della versione con il cromosoma MT (`git show 61fb36b:results/tables/malignancy_labels.csv`: PA06 0,615, PA08 0,452, PA10 0,578). La tabella attuale dice PA06 0,577, PA08 0,410, PA10 0,530, cioe' 41-58%. Il README (righe 22 e 227) e il REGISTRO (riga 262: "45-62% -> 41-58%") sono gia' aggiornati, la docstring no.
- correzione proposta: sostituire "45-62%" con "41-58%".

2. `README.md:17` (riga 3 della tabella dei tre problemi, colonna "entro paziente")
- gravita': MEDIA
- criterio violato: 4 (la descrizione dei 4 pazienti non corrisponde alla tabella che li seleziona).
- motivo: il testo dice "entro paziente i 4 pazienti con entrambi i tessuti (4.870 e 333 duttali)". In realta' i pazienti con entrambi i tessuti sono 5: PA04 ha 152 duttali Tumor e 1 Adjacent (`malignancy_labels.csv`; log di 14: "PA04: non valutabile (tumor n=152, normal n=1, minimo 10)"). I 4 di `cnv_loci_auc_by_patient` sono quelli con almeno 10 duttali per tessuto (`14_cnv_loci.py:365`), come conferma 333 = 334 − 1. La stessa cella dice anche "334 adiacenti (5 pazienti)", quindi il testo si contraddice. I numeri sono giusti; e' sbagliata la regola di selezione.
- correzione proposta: "entro paziente i 4 pazienti con almeno 10 duttali per tessuto (4.870 e 333 duttali)".

3. `README.md:240` (sez. 6, punto c) e `src/sim/sim_c_block_permutation.py:18-19`
- gravita': BASSA
- criterio violato: nessuno
- motivo: README e docstring dicono che il dato "2 falsi positivi su 100" non e' riprodotto. Pero' in `sim_c_block_permutation.csv` con confondente 2 vale `fp_first100_block_centered = 2`, cioe' esattamente 2 su 100 con il p centrato. Le prime 100 repliche coincidono con la tabella di 52ace44 (verificato: 4/6/11/5 non centrato, 6/7/4/2 centrato). L'affermazione e' vera per il p non centrato, quello usato in 17 e 18, ma il testo non lo precisa.
- correzione proposta: scrivere "non e' riprodotto con il p non centrato usato in 17/18; con il p centrato e confondente 2 le prime 100 repliche danno 2".

4. `tests/test_stats.py:9` (contro la riga 97)
- gravita': BASSA
- criterio violato: nessuno
- motivo: il docstring dice "semi 100 + k", ma il codice usa `seed=100 + k + int(binary)`, che era anche il codice di 12a356d. In piu' (k = 10, binario) e (k = 11, continuo) usano lo stesso seme, 111. Alla riga 12 "I valori 0,8875-0,9225 del README" non corrisponde a cio' che il README scrive (88,8-92,3%).
- correzione proposta: scrivere "semi 100 + k + (1 se binario)" (oppure usare semi distinti, per esempio 100 + 2k + binario) e "88,8-92,3% del README (0,8875-0,9225 in `sim_a_cluster_bootstrap.csv`)".

5. `src/13_infercnv.py:26-28`
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice "la release con cui sono state calcolate le cache dell'11/07/2026 (le coordinate delle finestre di `cnv_matrix.npz` coincidono con quelle dei geni di Ensembl 116)". Le cache attuali (`cnv_scores.csv`, `cnv_matrix.npz`) sono state ricalcolate il 18/10/2026, senza MT. Quelle dell'11/07 sono ora in `data/interim/cache_cnv_20260711/`. Chi legge pensa che le cache in uso siano quelle di luglio.
- correzione proposta: "release 116, la stessa delle cache dell'11/07/2026 (ora in `data/interim/cache_cnv_20260711/`); le cache attuali sono ricalcolate con questa release".

6. `src/14_cnv_loci.py:22-23`
- gravita': BASSA
- criterio violato: nessuno
- motivo: "l'AUC entro paziente ... isola il contrasto tissutale" e' troppo forte. Entro paziente Tumor e Adjacent_normal restano campioni e librerie diversi, quindi il contrasto comprende anche le differenze di campionamento e lavorazione. Toglie il confondimento con il paziente, non isola il tessuto.
- correzione proposta: "elimina le differenze fra pazienti (resta il contrasto fra campioni dello stesso paziente)".

7. `src/sim/sim_a_calibration.py:26-27`
- gravita': BASSA
- criterio violato: nessuno
- motivo: "1.000 repliche bootstrap ciascuna (le analisi ne usano 2.000: l'effetto sulla copertura e' trascurabile)". Nessuna esecuzione o tabella verifica che l'effetto sia trascurabile. Il README (riga 240) giustamente lo dice senza questa affermazione.
- correzione proposta: togliere "l'effetto ... trascurabile" oppure scrivere "non verificato con 2.000".

8. `src/sim/sim_a_calibration.py:214-222` (e docstring riga 45)
- gravita': BASSA
- criterio violato: nessuno
- motivo: `--quick` riduce le repliche a 1/10 ma scrive negli stessi `results/tables/sim_a_*.csv` citati nel README e in RISULTATI_FINALI, sovrascrivendo i numeri pubblicati. La docstring non lo avverte.
- correzione proposta: con `--quick` scrivere in file con un suffisso (per esempio `_quick`) oppure non scrivere nulla, e dichiararlo nella docstring.

9. `src/13_infercnv.py:586-592` e `:602-630`
- gravita': BASSA
- criterio violato: nessuno (non succede con i dati attuali)
- motivo: se tutti i fold LOPO sono saltati (`res` vuoto), `malignancy_lopo_summary.csv` non viene riscritto. Se le duttali adiacenti sono 20 o meno, `external_test_predictions.csv` e `external_test_summary.csv` non vengono scritti. In entrambi i casi restano in `results/tables` le versioni di un'esecuzione precedente, incoerenti con il nuovo `malignancy_lopo.csv`. In 08 il problema analogo e' stato corretto cancellando le tabelle vecchie (REGISTRO, giro 4).
- correzione proposta: cancellare questi file (`unlink(missing_ok=True)`) nei rami in cui non vengono prodotti, con un messaggio.

10. `src/sim/sim_b_reference_schemes.py:14-15` (contro la riga 109)
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice "CD4 > 0 nel 35% o nel 60% delle cellule tumorali di qualsiasi lignaggio". Il codice fa `det["CD4"] = max(det["CD4"], amb)`, quindi nel lignaggio CD4 con artefatto moderato il rilevamento resta 0,45, non 0,35.
- correzione proposta: "rilevamento di CD4 portato almeno al 35% (o 60%) nelle cellule tumorali di ogni lignaggio".

**Controlli senza problemi**
- **Numeri:** README (righe 17, 21, 22, 224-229 e sezione 6) e RISULTATI_FINALI.md (sezioni CNV e Simulazioni) coincidono con `results/tables`. Valori controllati: 8.370/334, 13/5 pazienti, 0,837 ± 0,071 su 12 pazienti, 30,8%, 41-58%, 0,531/0,698, 0,557 (0,528-0,588), 4.870/333. Per le simulazioni: (a) 2,7%/1,8%, 4,7%/6,1%, 99,3%; coperture 93,5-94,8%, 88,8-92,3%, 91,5-92,5%, 82-83%; (b) 0,126/0,093, 0,385/0,092/0,442, entro 0,005 e 0,008; (c) tutte le percentuali con gli IC, 4/6/11/5; (d) +0,143 (+0,116; +0,181), +0,030 (−0,013; +0,052).
- **Riproducibilita':** con un mio script nello scratchpad, che non scrive nulla nella repo, ho verificato tre cose.
  - L'impronta delle cache ricalcolata ora coincide con `cnv_cache_fingerprint.json`, per entrambe le voci.
  - Le soglie ricalcolate da `cnv_scores.csv` riproducono `malignancy_labels.csv` cella per cella.
  - `mean|X_cnv|` di `cnv_matrix.npz` coincide con `cnv_scores.csv` (differenza massima 1e-16) e da' AUC 0,5305. Le finestre sono solo chr1-22.
- **Esecuzione:** 13, 14, le simulazioni, `scpipe.py`, `10_loco.py` e `20_xai_per_patient.py` sono stati modificati prima delle riesecuzioni della chiusura 1 (date dei file contro `_status.txt`). L'albero git e' pulito e i log corrispondono ai messaggi del codice attuale.
- **Correzioni:** applicate.
  - A1: LogCPM e HVG sono nella `Pipeline` per fold e nel modello finale; `--splits` e' stato tolto.
  - A2: riferimento comune a tutti i pazienti, soglia stimata sulle stesse cellule, interpolazione delle finestre in 14.
  - A7: le AUC per paziente sono salvate in `cnv_loci_auc_by_patient`.
  - D2 13:24: Adjacent_normal dichiarato non indipendente da inferCNV.
  - Blocco C: le simulazioni (a)-(d) sono in `src/sim/` con docstring completa.
  - Chiusura 1 (c1-c10): corretti tutti, tranne la percentuale della docstring di 13 (problema 1).
- **Test:** ogni stimatore citato come validato ha un test con banda: cluster bootstrap con k = 10 e 11, block permutation, permutazione di Spearman con n = 10, Nadeau-Bengio con k = 3 e 5 e la potenza, discordanza e `d_within` (verifiche esatte). I test di leakage coprono 05, 11 (percorso sparso, denso con SVD, RF e NB), 12, 13, 20 e gli split di 04, ciascuno con un controllo positivo. pytest non scrive file nella repo (`git status` pulito dopo l'esecuzione).
- **Pacchetti:** tutti quelli importati nella mia area (anndata, scanpy, scipy, sklearn, matplotlib, infercnvpy, requests, pytest) sono in `envs/environment.yml`.
- **Sim d:** ho controllato il limite "vale anche per 20" sui dati reali di `task1a_paired`. Cellule in comune fra due training LOPO sottocampionati: 5.539; fra due sottoinsiemi del controllo: 5.912. Il rapporto e' 0,937, lo stesso della simulazione (1.800/1.920).

**Conteggio**
- ALTA: 0
- MEDIA: 2
- BASSA: 8
- pytest `-m "not slow" -p no:cacheprovider`: 102 passati, 0 falliti, 0 saltati (6 deselezionati, i lenti), in 96 s.
- ruff `check src tests revisione`: "All checks passed!"

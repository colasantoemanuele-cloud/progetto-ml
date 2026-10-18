## Esecuzione finale (Fase 3)

Comandi del README nell'ordine; log in `revisione/log/NN_nome.log`, stato e picchi di memoria in `revisione/log/_status.txt` (misurati con `/usr/bin/time -v`; durate e memoria vengono dai log, non da `results/tables/`). Prima esecuzione completa al commit `61fb36b`: dopo 00-04 tre rami, TCR (06-10, 21, 15-19) e 20 + simulazioni in parallelo, poi CNV (13, 14) da solo per la memoria. Dopo le correzioni dei giri 4-7 sono stati rieseguiti 04, 05 (giro 7: solo la legenda delle figure; tabelle identiche byte per byte), 06-10, 13-21 e le simulazioni (colonna "Ultima esecuzione"); il log di ogni script corrisponde all'ultima esecuzione.

| Script | Esito | Picco di memoria | Durata | Ultima esecuzione |
| :--- | :--- | :--- | :--- | :--- |
| 00_fetch_metadata | rc = 0 | 0,1 GB | 0:00,44 | Fase 3 (61fb36b) |
| 00b_download_tcr | rc = 0 | 0,0 GB | 0:00,88 | Fase 3 (61fb36b) |
| 01_probe_labels | rc = 0 | 0,2 GB | 0:01,02 | Fase 3 (61fb36b) |
| 01b_tcr_check | rc = 0 | 0,2 GB | 0:00,98 | Fase 3 (61fb36b) |
| 02_build_h5ad_sc | rc = 0 | 7,7 GB | 1:49,19 | Fase 3 (61fb36b) |
| 02_build_h5ad_sn | rc = 0 | 1,0 GB | 0:11,34 | Fase 3 (61fb36b) |
| 03_eda_qc_sc | rc = 0 | 0,4 GB | 0:07,01 | Fase 3 (61fb36b) |
| 03_eda_qc_sn | rc = 0 | 0,3 GB | 0:03,75 | Fase 3 (61fb36b) |
| 04_splits | rc = 0 | 0,4 GB | 0:04,74 | giro 6 |
| 06_clonotype_consistency | rc = 0 | 0,5 GB | 0:45,53 | giro 7 |
| 07_robustness_clone_ab | rc = 0 | 3,3 GB | 0:44,03 | giro 7 |
| 20_xai_per_patient | rc = 0 | 6,0 GB | 2:26,09 | giro 7 |
| 07_robustness_clone_b | rc = 0 | 3,3 GB | 0:50,16 | giro 7 |
| 08_disambiguate | rc = 0 | 3,5 GB | 0:46,07 | giro 7 |
| 09_arbitrate | rc = 0 | 3,4 GB | 1:03,74 | giro 7 |
| sim_a_calibration | rc = 0 | 0,2 GB | 4:33,42 | giro 7 |
| sim_b_reference_schemes | rc = 0 | 0,3 GB | 2:36,21 | giro 7 |
| 10_loco | rc = 0 | 3,4 GB | 7:57,35 | giro 7 |
| 21_fig_blood_scheme | rc = 0 | 0,3 GB | 0:03,87 | giro 7 |
| 15_inventory | rc = 0 | 0,3 GB | 0:03,10 | giro 7 |
| 16_pseudobulk_patient | rc = 0 | 3,6 GB | 0:21,68 | giro 7 |
| 17_clonal_expansion | rc = 0 | 1,2 GB | 0:13,23 | giro 7 |
| 18_loyo_check | rc = 0 | 1,1 GB | 0:24,28 | giro 7 |
| 19_ni_score | rc = 0 | 0,2 GB | 0:03,43 | giro 7 |
| sim_c_block_permutation | rc = 0 | 0,2 GB | 3:06,04 | giro 7 |
| sim_d_xai_stability | rc = 0 | 0,3 GB | 0:16,57 | giro 7 |
| 13_infercnv | rc = 0 | 10,8 GB | 1:39,14 | giro 7 |
| 14_cnv_loci | rc = 0 | 12,6 GB | 1:08,03 | giro 7 |
| 05_baseline | rc = 0 | 9,7 GB | 17:09,67 | giro 7 |

**Non rieseguiti** (calcolo invariato dall'esecuzione precedente, verificato con `git diff`):
- `12_feature_selection.py` (`52ace44`, 6,1 GB) e `11_model_zoo.py` (`8058701`, processi paralleli in `revisione/log/11_final/`, 5,2-5,8 GB): dopo quei commit cambiano solo controlli sugli argomenti e sugli input (fold con una sola classe, `cell_id` assenti, `nargs`), messaggi, docstring e `help`; per 11 riepilogo, confronti e figure sono stati rigenerati con `--from-table` (stesso codice).
- Le funzioni di `scpipe.py` usate da 05, 11 e 12 (`LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard`, `mean_pairwise_jaccard`) sono identiche a quelle del commit `52ace44` (confronto dell'AST senza docstring); `04_splits.py` e' cambiato solo nei controlli e ricrea `splits_sc.csv` identico byte per byte.

**Scorciatoie usate:** download gia' presenti non ripetuti (00, `download_tcr.sh`, 01, 02); cache di Scrublet (`doublet_scores.csv`, 08), ora riusata solo se contiene esattamente le cellule dell'h5ad (199.184, controllo passato); la cache e' del 10/07 e non e' stata ricalcolata per intero: il revisore della chiusura 1 ha ricalcolato Scrublet con il codice attuale su due librerie (PDAC10-PBMC e PDAC14-Tumor) e ha trovato una differenza massima di circa 1e-16 (`revisione/log/chiusura1_revisore_b.md`, punto 8). Nessun numero del README dipende da questa cache. inferCNV ricalcolato senza cache nel giro di chiusura 1, dopo l'esclusione del cromosoma mitocondriale e dei contig non standard.

**Riproducibilita':** `02_build_h5ad.py` ricrea `sc_raw.h5ad` identico byte per byte e `04_splits.py` ricrea `splits_sc.csv` identico. Nella chiusura 1 sono cambiati i risultati di 13 e 14 (cromosomi), le colonne aggiunte a `clone_switching_summary_clone_b` e `error_rate_by_compartment` e le figure di 11 (etichetta dell'asse); nella chiusura 2 solo le figure `excess_by_pair_{cc}` e `error_rate_by_compartment`, nel giro 6 solo `error_rate_blood` ed `error_rate_loco` (IC non calibrati segnati); nel giro 6 inferCNV e' stato ricalcolato senza cache (impronta con versione e parametri) con risultati identici byte per byte.

Nel giro 7 sono cambiati i risultati del ramo TCR (pareggi dell'etichetta di maggioranza esclusi) e del ramo CNV (soglia del rumore globale, coordinate reali delle finestre), e lo scenario `coverage_pairs_real_sizes` di `sim_a`; vedi la tabella prima/dopo.

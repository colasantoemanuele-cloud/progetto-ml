## Esecuzione finale (Fase 3)

Comandi del README nell'ordine; log in `revisione/log/NN_nome.log`, stato e picchi di memoria in `revisione/log/_status.txt` (misurati con `/usr/bin/time -v`; durate e memoria vengono dai log, non da `results/tables/`). Prima esecuzione completa al commit `61fb36b`: dopo 00-04 tre rami, TCR (06-10, 21, 15-19) e 20 + simulazioni in parallelo, poi CNV (13, 14) da solo per la memoria. Dopo le correzioni dei giri di chiusura 1 e 2 sono stati rieseguiti 04 (chiusura 1), 06-10, 13-21 e le simulazioni (colonna "Ultima esecuzione"); il log di ogni script corrisponde all'ultima esecuzione.

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
| 04_splits | rc = 0 | 0,4 GB | 0:04,68 | chiusura 1 |
| 06_clonotype_consistency | rc = 0 | 0,5 GB | 0:46,28 | chiusura 2 |
| 07_robustness_clone_ab | rc = 0 | 3,3 GB | 0:36,97 | chiusura 2 |
| 20_xai_per_patient | rc = 0 | 6,0 GB | 1:41,95 | chiusura 2 |
| 07_robustness_clone_b | rc = 0 | 3,3 GB | 0:39,28 | chiusura 2 |
| 08_disambiguate | rc = 0 | 3,5 GB | 0:36,15 | chiusura 2 |
| 09_arbitrate | rc = 0 | 3,4 GB | 0:32,12 | chiusura 2 |
| sim_a_calibration | rc = 0 | 0,2 GB | 3:08,38 | chiusura 2 |
| sim_b_reference_schemes | rc = 0 | 0,3 GB | 1:26,46 | chiusura 2 |
| 10_loco | rc = 0 | 3,4 GB | 3:54,20 | chiusura 2 |
| 21_fig_blood_scheme | rc = 0 | 0,3 GB | 0:02,18 | chiusura 2 |
| 15_inventory | rc = 0 | 0,3 GB | 0:01,77 | chiusura 2 |
| 16_pseudobulk_patient | rc = 0 | 3,6 GB | 0:14,10 | chiusura 2 |
| 17_clonal_expansion | rc = 0 | 1,2 GB | 0:08,74 | chiusura 2 |
| 18_loyo_check | rc = 0 | 1,2 GB | 0:19,31 | chiusura 2 |
| 19_ni_score | rc = 0 | 0,2 GB | 0:02,54 | chiusura 2 |
| sim_c_block_permutation | rc = 0 | 0,2 GB | 2:22,74 | chiusura 2 |
| sim_d_xai_stability | rc = 0 | 0,3 GB | 0:11,68 | chiusura 2 |
| 13_infercnv | rc = 0 | 6,0 GB | 0:56,21 | chiusura 2 |
| 14_cnv_loci | rc = 0 | 8,6 GB | 0:16,79 | chiusura 2 |

**Non rieseguiti** (calcolo invariato dall'esecuzione precedente, verificato con `git diff`):
- `05_baseline.py` (eseguito al commit `52ace44`, picco 9,7 GB), `12_feature_selection.py` (`52ace44`, 6,1 GB) e `11_model_zoo.py` (`8058701`, processi paralleli in `revisione/log/11_final/`, 5,2-5,8 GB): dopo quei commit cambiano solo controlli sugli argomenti e sugli input (fold con una sola classe, `cell_id` assenti, `nargs`), messaggi, docstring e `help`; per 11 riepilogo, confronti e figure sono stati rigenerati con `--from-table` (stesso codice).
- Le funzioni di `scpipe.py` usate da 05, 11 e 12 (`LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard`, `mean_pairwise_jaccard`) sono identiche a quelle del commit `52ace44` (confronto dell'AST senza docstring); `04_splits.py` e' cambiato solo nei controlli e ricrea `splits_sc.csv` identico byte per byte.

**Scorciatoie usate:** download gia' presenti non ripetuti (00, `download_tcr.sh`, 01, 02); cache di Scrublet (`doublet_scores.csv`, 08), ora riusata solo se contiene esattamente le cellule dell'h5ad (199.184, controllo passato); la cache e' del 10/07 e non e' stata ricalcolata per intero: il revisore della chiusura 1 ha ricalcolato Scrublet con il codice attuale su due librerie (PDAC10-PBMC e PDAC14-Tumor) e ha trovato una differenza massima di circa 1e-16 (`revisione/log/chiusura1_revisore_b.md`, punto 8). Nessun numero del README dipende da questa cache. inferCNV ricalcolato senza cache nel giro di chiusura 1, dopo l'esclusione del cromosoma mitocondriale e dei contig non standard.

**Riproducibilita':** `02_build_h5ad.py` ricrea `sc_raw.h5ad` identico byte per byte e `04_splits.py` ricrea `splits_sc.csv` identico. Nella chiusura 1 sono cambiati i risultati di 13 e 14 (cromosomi), le colonne aggiunte a `clone_switching_summary_clone_b` e `error_rate_by_compartment` e le figure di 11 (etichetta dell'asse); nella chiusura 2 solo le figure `excess_by_pair_{cc}` e `error_rate_by_compartment` (IC non calibrati segnati).

## Esecuzione finale (Fase 3)

Codice: commit `61fb36b` (`revisione/log/_fase3_commit.txt`); comandi del README nell'ordine, log in `revisione/log/NN_nome.log` (stato e picchi di memoria in `revisione/log/_status.txt`, misurati con `/usr/bin/time -v`). Dopo 00-04 tre rami: TCR (06-10, 21, 15-19) e 20 + simulazioni in parallelo, poi CNV (13, 14) da solo per la memoria.

| Script | Esito | Picco di memoria | Durata |
| :--- | :--- | :--- | :--- |
| 00_fetch_metadata | rc = 0 | 0,1 GB | 0:00,44 |
| 00b_download_tcr | rc = 0 | 0,0 GB | 0:00,88 |
| 01_probe_labels | rc = 0 | 0,2 GB | 0:01,02 |
| 01b_tcr_check | rc = 0 | 0,2 GB | 0:00,98 |
| 02_build_h5ad_sc | rc = 0 | 7,7 GB | 1:49,19 |
| 02_build_h5ad_sn | rc = 0 | 1,0 GB | 0:11,34 |
| 03_eda_qc_sc | rc = 0 | 0,4 GB | 0:07,01 |
| 03_eda_qc_sn | rc = 0 | 0,3 GB | 0:03,75 |
| 04_splits | rc = 0 | 0,4 GB | 0:04,63 |
| 06_clonotype_consistency | rc = 0 | 0,5 GB | 0:52,10 |
| 07_robustness_clone_ab | rc = 0 | 3,3 GB | 0:44,19 |
| 20_xai_per_patient | rc = 0 | 6,0 GB | 1:38,96 |
| 07_robustness_clone_b | rc = 0 | 3,3 GB | 0:41,05 |
| 08_disambiguate | rc = 0 | 3,5 GB | 0:36,27 |
| 09_arbitrate | rc = 0 | 3,4 GB | 0:32,38 |
| sim_a_calibration | rc = 0 | 0,2 GB | 3:13,13 |
| sim_b_reference_schemes | rc = 0 | 0,3 GB | 1:40,64 |
| 10_loco | rc = 0 | 3,4 GB | 3:44,39 |
| 21_fig_blood_scheme | rc = 0 | 0,3 GB | 0:02,35 |
| 15_inventory | rc = 0 | 0,3 GB | 0:01,97 |
| 16_pseudobulk_patient | rc = 0 | 3,6 GB | 0:15,42 |
| 17_clonal_expansion | rc = 0 | 1,2 GB | 0:08,64 |
| 18_loyo_check | rc = 0 | 1,2 GB | 0:18,91 |
| 19_ni_score | rc = 0 | 0,2 GB | 0:02,51 |
| sim_c_block_permutation | rc = 0 | 0,2 GB | 2:13,98 |
| sim_d_xai_stability | rc = 0 | 0,3 GB | 0:08,43 |
| 13_infercnv | rc = 0 | 6,0 GB | 0:54,27 |
| 14_cnv_loci | rc = 0 | 8,6 GB | 0:14,56 |

**Non rieseguiti** (codice di calcolo invariato dall'esecuzione precedente, verificato con `git diff`):
- `05_baseline.py`: nessuna modifica dal commit `52ace44` (eseguito in Fase 1, picco 9,7 GB).
- `12_feature_selection.py`: dal commit `52ace44` cambia solo il controllo sui `cell_id` assenti (eseguito in Fase 1, picco 6,1 GB).
- `11_model_zoo.py`: dal commit `8058701` cambiano solo docstring e argomenti di `--from-table`; esecuzione finale in processi paralleli in `revisione/log/11_final/` (picchi 5,2-5,8 GB), unita e confrontata con un'esecuzione unica.
- Le funzioni di `scpipe.py` usate da 05, 11 e 12 (`LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard`, `mean_pairwise_jaccard`) sono identiche a quelle del commit `52ace44` (confronto dell'AST senza docstring); `04_splits.py` e' cambiato solo nei controlli e ricrea `splits_sc.csv` identico byte per byte.

**Scorciatoie usate:** download gia' presenti non ripetuti (00, `download_tcr.sh`, 01, 02); cache di Scrublet (`doublet_scores.csv`, 08); cache di inferCNV (`cnv_scores.csv`, `cnv_matrix.npz`), ricalcolate senza cache nel giro 1 con il codice attuale e risultate identiche byte per byte a quelle dell'11/07 (13 senza cache: picco 10,8 GB; 14: 11,4 GB).

**Riproducibilita':** `02_build_h5ad.py` ricrea `sc_raw.h5ad` identico byte per byte e `04_splits.py` ricrea `splits_sc.csv` identico; tutte le tabelle e figure rigenerate in Fase 3 coincidono byte per byte con quelle dell'ultimo giro, tranne l'ordine delle righe di `ni_score_holm_joint.csv` (stessi valori). Nessun numero del README e' cambiato.

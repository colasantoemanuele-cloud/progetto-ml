# Risultati finali della revisione

Generato da `revisione/genera_risultati.py`: ogni numero viene dalla tabella di `results/tables/` indicata accanto. "Prima" e' il valore della presentazione (o del README pubblicato, dove indicato); i metadati di esecuzione in fondo (tempi, memoria) vengono dai log in `revisione/log/`, non da `results/tables/`.

## Numeri della presentazione: prima e dopo

| Numero | Prima | Dopo | Esito | Tabella |
| :--- | :--- | :--- | :--- | :--- |
| macro-F1 task2_v1, 5 fold per paziente (logistica L2) | **0,954** | 0,954 | uguale | `baseline_summary_sc` |
| macro-F1 task1a_paired, split casuale | **0,875±0,002** | 0,875±0,002 | uguale | `baseline_summary_sc` |
| macro-F1 task1a_paired, LOPO | **0,729±0,100** | 0,729±0,100 | uguale | `baseline_summary_sc` |
| LightGBM task2_v1 (LOPO) | **0,962** | 0,962 | uguale | `model_zoo_summary_sc` |
| logistica L2 task2_hier3 | **0,851** | **0,850** | arrotondamento (valore esatto a meno di 0,001 dal precedente) | `model_zoo_summary_sc` |
| LightGBM task2_hier3 | **0,923** | 0,923 | uguale | `model_zoo_summary_sc` |
| MLP task2_hier3 (README pubblicato) | **0,903** | **0,904** | cambiato | `model_zoo_summary_sc` |
| RF vs LightGBM task2_v1: p Wilcoxon / p Nadeau-Bengio | **p 0,002/0,074** | p 0,002/0,074 | uguale | `model_comparison_task2_v1` |
| dummy task1a_paired: Δ dal migliore / p Nadeau-Bengio | **0,229/0,072** | 0,229/0,072 | uguale | `model_comparison_task1a_paired` |
| selezione feature, k casuali k = 50, task2_v1 | **0,361** | **0,367** | cambiato | `feature_selection_summary_task2_v1` |
| selezione feature, L1 k = 50, task2_v1 | **0,936** | 0,936 | uguale | `feature_selection_summary_task2_v1` |
| selezione feature, L1 k = 200, task2_v1 | **0,952** | **0,955** | cambiato | `feature_selection_summary_task2_v1` |
| Jaccard L1 k = 200, task2_v1 | **0,439** | **0,480** | cambiato | `feature_stability_task2_v1` |
| geni comuni ai 5 fold (task1a, k = 50): HVG / chi2 / ANOVA | **26/50, 19/50, 20/50** | **18/50, 19/50, 19/50** | cambiato | `feature_core_genes_task1a_paired` |
| XAI: Jaccard LOPO | **0,300** | 0,300 | uguale | `xai_summary_task1a_paired` |
| XAI: Jaccard controllo | **0,459** | 0,459 | uguale | `xai_summary_task1a_paired` |
| sim (d): differenza con / senza effetto paziente | **0,063/0,040** | **0,143/0,030** | cambiato | `sim_d_xai_stability` |
| rumore di base d_within | **0,107** | 0,107 | uguale | `excess_summary_clone_b` |
| eccesso fra compartimenti [IC 95%] | **+0,059 [0,034;0,100]** | **+0,059 [0,031; 0,110]** | cambiato | `excess_summary_clone_b` |
| errore Tumor, schema solo sangue [IC] | **0,195 [0,123;0,274]** | **0,170 [0,095; 0,252]** | cambiato | `error_rate_blood` |
| errore Adjacent_normal, schema solo sangue | **0,059** | **0,060** | cambiato | `error_rate_blood` |
| Tumor − Adjacent appaiata, solo sangue [IC] | **+0,080 [+0,021;+0,170]** | **+0,059 [+0,000; +0,160]** | cambiato | `error_rate_paired_diff` |
| CD8A > 0 nei cloni discordanti / concordanti (Tumor) | **7%/88%** | 7%/88% | uguale | `marker_test_disc_vs_conc` |
| doppietti d1 | **0,723** | 0,723 | uguale | `doublet_mixture` |
| doppietti d2 | **0,790** | **0,676** | cambiato | `doublet_mixture` |
| espansione clonale: ρ errore / ρ eccesso | **ρ 0,18/0,28** | **ρ 0,22/0,28** | cambiato | `clonal_expansion_tests` |
| sim (c): falsi positivi su 100 repliche / su 1.000 (confondente 0; 0,5; 1; 2) | **2/100** | **4/6/11/5 su 100; 5,0%/5,9%/4,8%/5,1% su 1.000** | cambiato | `sim_c_block_permutation` |
| H1: ρ, p_Holm (famiglia di 4) | **ρ -0,47, p_Holm 0,53** | **ρ -0,45, p_Holm 0,58** | cambiato | `pseudobulk_hypothesis_tests, ni_score_test, ni_score_holm_joint` |
| H2: ρ, p_Holm (famiglia di 4) | **ρ -0,54, p_Holm 0,45** | **ρ -0,50, p_Holm 0,58** | cambiato | `pseudobulk_hypothesis_tests, ni_score_test, ni_score_holm_joint` |
| H3: ρ, p_Holm (famiglia di 4) | **ρ -0,15, p_Holm 1,00** | **ρ -0,13, p_Holm 1,00** | cambiato | `pseudobulk_hypothesis_tests, ni_score_test, ni_score_holm_joint` |
| H4: ρ, p_Holm (famiglia di 4) | **ρ -0,14, p_Holm 1,00** | **ρ -0,13, p_Holm 1,00** | cambiato | `pseudobulk_hypothesis_tests, ni_score_test, ni_score_holm_joint` |
| AUC LOPO malignita' (13) | **0,825** | **0,836** | cambiato | `malignancy_lopo_summary` |
| duttali adiacenti predette maligne | **31,7%** | **37,4%** | cambiato | `external_test_summary` |
| AUC pooled score CNV medio | **0,492** | **0,598** | cambiato | `cnv_loci_auc` |
| AUC pooled loci PDAC | **0,682** | **0,717** | cambiato | `cnv_loci_auc` |
| AUC loci PDAC entro paziente (media dei 4 pazienti valutabili) | **0,53** | **0,566** | cambiato | `cnv_loci_auc_by_patient` |
| sim (b): errore PBMC stimato da LOCO con artefatto (moderato / forte) | **0,24** | **0,119 / 0,383** | cambiato | `sim_b_reference_schemes` |

Note: lo 0,53 e' l'AUC entro paziente dei loci PDAC ("≈ 0,53" nel README pubblicato; il codice di allora stampava solo le AUC per paziente, quindi la statistica riassuntiva non e' ricostruibile; ora si riporta la media). I p grezzi di H1-H4 sono nella sezione delle ipotesi. Gli IC con 5 pazienti o meno (Adjacent_normal, differenza appaiata) non sono calibrati (vedi Simulazioni).

## Dataset

scRNA-seq: 199.184 cellule, 14 pazienti, 30.543 geni, 29 librerie; PBMC 89.506 cellule (11 pazienti), Tumor 71.414 (13), Adjacent_normal 38.264 (5) (`dataset_summary_sc`). snRNA-seq: 32.005 nuclei, 8 pazienti (`dataset_summary_sn`).

| Task | Cellule | Classi | Pazienti | Valutazione |
| :--- | :--- | :--- | :--- | :--- |
| task1a_all | 109.618 | 2 | 13 | group_5fold |
| task1a_paired | 67.662 | 2 | 5 | LOPO |
| task2_v1 | 70.547 | 9 | 13 | group_5fold |
| task2_v2 | 70.547 | 8 | 13 | group_5fold |
| task2_hier1 | 71.377 | 3 | 13 | group_5fold |
| task2_hier3 | 42.168 | 3 | 13 | group_5fold |
| rare | 14.854 | 5 | 13 | LOPO |

Fonte: `task_summary_sc`.

## Leakage: split casuale contro split per paziente (05)

| Task | Per paziente | Casuale | Δ (paziente − casuale) | Rapporto sd |
| :--- | :--- | :--- | :--- | :--- |
| task1a_all | 0,682 | 0,878 | -0,196 | 86,0 |
| task1a_paired | 0,729 | 0,875 | -0,146 | 50,3 |
| task2_hier1 | 0,999 | 0,999 | -0,0002 | 2,9 |
| task2_hier3 | 0,866 | 0,885 | -0,019 | 13,9 |
| task2_v1 | 0,954 | 0,964 | -0,009 | 12,8 |
| task2_v2 | 0,951 | 0,961 | -0,010 | 19,1 |

Fonte: `baseline_gap_sc` (per paziente = LOPO per task1a_paired, 5 fold per paziente altrove).

## 13 modelli (11, LOPO, 20.000 cellule di training per fold)

**task2_v1** (`model_zoo_summary_sc`)

| Modello | macro-F1 | sd | Δ da logreg_l2 | fold non convergenti |
| :--- | :--- | :--- | :--- | :--- |
| lightgbm | 0,962 | 0,010 | +0,019 |  |
| mlp | 0,958 | 0,014 | +0,015 | 0/13 |
| svm_rbf | 0,957 | 0,016 | +0,013 |  |
| random_forest | 0,956 | 0,012 | +0,012 |  |
| naive_bayes | 0,954 | 0,014 | +0,010 |  |
| logreg_l2 | 0,943 | 0,016 | +0,000 | 0/13 |
| logreg_l1 | 0,943 | 0,019 | -0,001 | 13/13 |
| elasticnet | 0,939 | 0,017 | -0,004 | 13/13 |
| svm_linear | 0,936 | 0,014 | -0,007 | 0/13 |
| knn | 0,931 | 0,027 | -0,013 |  |
| bagging | 0,919 | 0,023 | -0,025 |  |
| decision_tree | 0,865 | 0,033 | -0,078 |  |
| dummy | 0,099 | 0,008 | -0,845 |  |

**task2_hier3** (`model_zoo_summary_sc`)

| Modello | macro-F1 | sd | Δ da logreg_l2 | fold non convergenti |
| :--- | :--- | :--- | :--- | :--- |
| lightgbm | 0,923 | 0,028 | +0,072 |  |
| svm_rbf | 0,906 | 0,035 | +0,056 |  |
| mlp | 0,904 | 0,032 | +0,053 | 0/13 |
| random_forest | 0,901 | 0,033 | +0,051 |  |
| naive_bayes | 0,882 | 0,052 | +0,031 |  |
| bagging | 0,881 | 0,043 | +0,030 |  |
| logreg_l1 | 0,877 | 0,038 | +0,027 | 13/13 |
| svm_linear | 0,865 | 0,041 | +0,015 | 0/13 |
| elasticnet | 0,862 | 0,047 | +0,012 | 13/13 |
| knn | 0,856 | 0,043 | +0,006 |  |
| decision_tree | 0,853 | 0,044 | +0,003 |  |
| logreg_l2 | 0,850 | 0,045 | +0,000 | 0/13 |
| dummy | 0,329 | 0,011 | -0,521 |  |

**task1a_paired** (`model_zoo_summary_task1`)

| Modello | macro-F1 | sd | Δ da logreg_l2 | fold non convergenti |
| :--- | :--- | :--- | :--- | :--- |
| svm_rbf | 0,725 | 0,101 | +0,003 |  |
| logreg_l2 | 0,722 | 0,082 | +0,000 | 0/5 |
| svm_linear | 0,721 | 0,085 | -0,001 | 0/5 |
| mlp | 0,721 | 0,100 | -0,001 | 0/5 |
| lightgbm | 0,720 | 0,127 | -0,002 |  |
| random_forest | 0,715 | 0,107 | -0,008 |  |
| logreg_l1 | 0,712 | 0,075 | -0,010 | 5/5 |
| elasticnet | 0,710 | 0,082 | -0,012 | 5/5 |
| knn | 0,710 | 0,086 | -0,012 |  |
| bagging | 0,681 | 0,119 | -0,041 |  |
| decision_tree | 0,666 | 0,098 | -0,056 |  |
| naive_bayes | 0,645 | 0,116 | -0,077 |  |
| dummy | 0,496 | 0,002 | -0,226 |  |

- task2_v1: migliore macro-F1 0,962; 8 modelli su 12 con p di Nadeau-Bengio < 0,05 rispetto al migliore (`model_comparison_task2_v1`).
- task2_hier3: migliore macro-F1 0,923; 12 modelli su 12 con p di Nadeau-Bengio < 0,05 rispetto al migliore (`model_comparison_task2_hier3`).
- task1a_paired: migliore macro-F1 0,725; 0 modelli su 12 con p di Nadeau-Bengio < 0,05 rispetto al migliore (`model_comparison_task1a_paired`).

Fold per task: 5-13. Con pochi fold (task1a_paired, 5 fold) la potenza e' bassa: un p non significativo non dimostra l'equivalenza. I p sono grezzi (nessuna correzione per i confronti multipli) e il migliore e' scelto sugli stessi fold. La calibrazione di Nadeau-Bengio e' verificata in `sim_a` con k = 3, 5 fold stratificati di uguale dimensione, non con fold LOPO di dimensione diversa.

## Selezione delle feature (12)

**task2_v1** (`feature_selection_summary_task2_v1`, `feature_stability_task2_v1`)

| Metodo | F1 k=50 | F1 k=200 | F1 k=1000 | Jaccard k=50 | Jaccard k=200 | Jaccard k=1000 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| l1_logreg | 0,936 | 0,955 | 0,953* | 0,514 | 0,480 | 0,568* |
| rf_import | 0,926 | 0,950 | 0,953 | 0,534 | 0,595 | 0,662 |
| chi2 | 0,908 | 0,944 | 0,950 | 0,787 | 0,842 | 0,851 |
| anova | 0,908 | 0,938 | 0,949 | 0,816 | 0,794 | 0,846 |
| hvg | 0,659 | 0,879 | 0,946 | 0,528 | 0,793 | 0,780 |
| random | 0,367 | 0,648 | 0,874 | 0,001 | 0,003 | 0,025 |

\* k maggiore dei geni con punteggio > 0 in almeno un fold (`any_k_gt_nonzero`): i geni oltre quelli con punteggio > 0 sono scelti fra pareggi in ordine arbitrario, quindi F1 e Jaccard a quel k non sono interpretabili.

L1 (SAGA): 0/5 fold convergenti, 524-583 coefficienti non nulli (`feature_selection_diag_task2_v1`).

**task1a_paired** (`feature_selection_summary_task1a_paired`, `feature_stability_task1a_paired`)

| Metodo | F1 k=50 | F1 k=200 | F1 k=1000 | Jaccard k=50 | Jaccard k=200 | Jaccard k=1000 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| l1_logreg | 0,712 | 0,751 | 0,751 | 0,273 | 0,264 | 0,254 |
| rf_import | 0,735 | 0,744 | 0,748 | 0,469 | 0,539 | 0,646 |
| chi2 | 0,696 | 0,716 | 0,742 | 0,522 | 0,594 | 0,678 |
| anova | 0,706 | 0,718 | 0,733 | 0,519 | 0,572 | 0,665 |
| hvg | 0,623 | 0,655 | 0,716 | 0,458 | 0,725 | 0,771 |
| random | 0,589 | 0,678 | 0,675 | 0,000 | 0,004 | 0,027 |

L1 (SAGA): 5/5 fold convergenti, 1616-1984 coefficienti non nulli (`feature_selection_diag_task1a_paired`).

- hvg: 18 geni in tutti i 5 fold LOPO di task1a_paired (k = 50): APOD, CCL19, CGA, CPA3, CST1, CXCL13, GNLY, IGFBP7, IGHV3-15, IGKV3-15, LUM, MGP, SFRP2, SPP1, TPSAB1, TPSB2, TPSD1, TRBV7-2
- chi2: 19 geni in tutti i 5 fold LOPO di task1a_paired (k = 50): C19orf33, CD69, CD9, CLDN4, CST3, CYSTM1, ELF3, FCMR, IFI27, IGFBP7, KRT18, KRT19, KRT7, KRT8, LMNA, LTB, S100P, SPINK1, TPM1
- anova: 19 geni in tutti i 5 fold LOPO di task1a_paired (k = 50): CD151, CD37, CD48, CD69, CD9, CORO1A, CXCR4, CYSTM1, EEF1A1, FCMR, GPSM3, HLA-E, KRT19, LGALS3, LIMD2, LMNA, LTB, S100A11, S100A6
- rf_import: 16 geni in tutti i 5 fold LOPO di task1a_paired (k = 50): CD69, CXCR4, EEF1A1, EEF1G, FAU, FCMR, HSPA1A, HSPA1B, LIMD2, LMNA, LTB, NACA, SELL, TAGAP, UBA52, VIM
- l1_logreg: 4 geni in tutti i 5 fold LOPO di task1a_paired (k = 50): CD69, DUSP4, MT-CO2, SPP1

Fonte: `feature_core_genes_task1a_paired`.

## Spiegazioni per paziente (20)

Jaccard medio fra i top-50 dei 5 fold LOPO 0,300 (0,220-0,408 fra coppie, `xai_jaccard_task1a_paired`), controllo casuale 0,459, differenza +0,159, Spearman medio 0,234 (`xai_summary_task1a_paired`). Geni comuni a tutti i fold (9): BTG2, CD69, CREM, CTRB1, DUSP4, FOS, GAPDH, HSPH1, VIM (`xai_common_genes_task1a_paired`). In simulazione la differenza e' +0,143 con effetto paziente e +0,030 senza (`sim_d_xai_stability`).

## Validazione clonale con il TCR

- Rumore di base (stesso clone, stesso compartimento): 0,107; discordanza fra compartimenti 0,166; eccesso +0,059 IC [+0,031; +0,110] (11 pazienti, 1.114 coppie; `excess_summary_clone_b`).
- Sensibilita' (B3): pesi uguali +0,074 IC [+0,050; +0,119]; clone TRA+TRB +0,065 IC [+0,035; +0,121] (`excess_summary_clone_ab`); split-half di 06 +0,063 IC [+0,023; +0,134] (14 pazienti: 11 nel termine fra compartimenti, 14 nel basale; `clone_switching_summary_clone_b`).

**Tasso d'errore per compartimento e schema di riferimento** (`error_rate_schemes`; \* IC non calibrato, <= 5 pazienti)

| Schema | Compartimento | Errore | IC 95% | Cloni | Pazienti |
| :--- | :--- | :--- | :--- | :--- | :--- |
| blood | Adjacent_normal | 0,060 | [0,000; 0,095]* | 117 | 4 |
| blood | Tumor | 0,170 | [0,095; 0,252] | 277 | 10 |
| loco | PBMC | 0,036 | [0,013; 0,063] | 309 | 10 |
| loco | Adjacent_normal | 0,062 | [0,031; 0,094]* | 292 | 5 |
| loco | Tumor | 0,128 | [0,076; 0,214] | 453 | 11 |
| naive | PBMC | 0,041 | [0,025; 0,064] | 857 | 11 |
| naive | Adjacent_normal | 0,072 | [0,046; 0,160]* | 639 | 5 |
| naive | Tumor | 0,084 | [0,058; 0,115] | 1.502 | 13 |

**Differenze appaiate** (`error_rate_paired_diff`; B2: ricampionature distinte e frazione > 0 con <= 5 pazienti)

| Schema | Variante | Contrasto | Diff. | IC 95% | Cloni | Pazienti | Ricamp. distinte | Frac > 0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| naive | main | Tumor - Adjacent_normal | +0,024 | [+0,000; +0,110] | 253 | 5 | 122 | 0,989 |
| naive | main | Tumor - PBMC | +0,129 | [+0,029; +0,212] | 272 | 10 |  |  |
| loco | main | Tumor - Adjacent_normal | +0,016 | [-0,020; +0,095] | 249 | 5 | 122 | 0,813 |
| loco | main | Tumor - PBMC | +0,121 | [+0,021; +0,205] | 265 | 10 |  |  |
| blood | main | Tumor - Adjacent_normal | +0,059 | [+0,000; +0,160] | 85 | 4 | 35 | 0,933 |

**Scomposizione dell'errore nel Tumor, schema solo sangue (B1)** (`error_decomposition_blood`)

| Categoria | Errori | Tasso | Cloni | IC 95% |
| :--- | :--- | :--- | :--- | :--- |
| CD4_in_blood_labelled_CD8T | 32 | 0,116 | 277 |  |
| labelled_NK | 15 | 0,054 | 277 |  |
| CD8_in_blood_labelled_CD4T | 0 | 0,000 | 277 |  |
| other | 0 | 0,000 | 277 |  |
| transition_CD4T_to_CD8T | 32 | 0,116 | 277 |  |
| transition_CD4T_to_NK | 4 | 0,014 | 277 |  |
| transition_CD8T_to_NK | 11 | 0,040 | 277 |  |
| ties_excluded | ⟦n.d.⟧ | ⟦n.d.⟧ | 10 |  |
| error_rate_excluding_NK_labels | 32 | 0,122 | 262 | [0,049; 0,224] |

**Regola d'identita' CD4/CD8 (B6)** (`error_rate_identity_rules`)

| Regola | Compartimento | Errore | IC 95% | Cloni |
| :--- | :--- | :--- | :--- | :--- |
| CD8A > 0 oppure CD8B > 0 (principale) | Adjacent_normal | 0,060 | [0,000; 0,095] | 117 |
| CD8A > 0 oppure CD8B > 0 (principale) | Tumor | 0,170 | [0,095; 0,252] | 277 |
| CD8B > 0 | Adjacent_normal | 0,052 | [0,000; 0,071] | 116 |
| CD8B > 0 | Tumor | 0,180 | [0,099; 0,271] | 272 |
| CD8A >= 2 UMI | Adjacent_normal | 0,052 | [0,000; 0,071] | 115 |
| CD8A >= 2 UMI | Tumor | 0,181 | [0,098; 0,269] | 277 |
| CD8B > 0 oppure CD8A >= 2 UMI | Adjacent_normal | 0,052 | [0,000; 0,071] | 116 |
| CD8B > 0 oppure CD8A >= 2 UMI | Tumor | 0,178 | [0,098; 0,265] | 276 |

Marcatori nelle cellule etichettate CD8T del Tumor (`marker_test_disc_vs_conc`): CD8A > 0 nel 6,6% dei discordanti contro 88,4% dei concordanti; rapporto delle medie di log1p(CP10k) di CD4 16,3; GNLY 1,78 contro 0,44.

**Cloni presenti nei tre compartimenti (08)** (`triple_clones_summary`; 4 pazienti: IC non calibrati)

| Contrasto | Media | IC 95% | Cloni | Pazienti |
| :--- | :--- | :--- | :--- | :--- |
| exc_PBM_Adj | +0,036 | [-0,021; +0,057] | 123 | 4 |
| exc_PBM_Tum | +0,084 | [-0,002; +0,176] | 123 | 4 |
| exc_Adj_Tum | +0,018 | [-0,014; +0,049] | 123 | 4 |
| exc_PBM_Tum - exc_Adj_Tum | +0,066 | [+0,012; +0,126] | 123 | 4 |
| exc_PBM_Tum - exc_PBM_Adj | +0,048 | [-0,006; +0,122] | 123 | 4 |

### Doppietti fra le NK con TCR (09, 10)

Tassi di cellule con catena TRB produttiva per classe: `tcr_rate_by_celltype` (14 classi).

| Statistica | Strato | Pazienti | d1 | IC d1 | d2 | IC d2 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| median | all | 14 | 0,723 | [0,556; 0,848] | 0,676 | [0,559; 0,781] |
| median | PBMC | 11 | 0,842 | [0,656; 1,000] | 0,851 | [0,623; 0,983] |
| median | Adjacent_normal | 5 | 0,690 | [0,233; 1,000] | 0,782 | [0,051; 1,000] |
| median | Tumor | 13 | 0,510 | [0,186; 0,658] | 0,419 | [0,289; 0,567] |
| mean | all | 14 | 0,682 | [0,582; 0,775] | 0,693 | [0,563; 0,790] |
| mean | PBMC | 11 | 0,800 | [0,690; 1,000] | 0,855 | [0,659; 0,997] |
| mean | Adjacent_normal | 5 | 0,559 | [0,299; 0,862] | 0,553 | [0,126; 1,000] |
| mean | Tumor | 13 | 0,381 | [0,170; 0,473] | 0,515 | [0,325; 0,675] |

Fonte: `doublet_mixture` (B5: medie e strati per compartimento). d1 e d2 non sono indipendenti; la copertura dei loro IC non e' verificata.

### Espansione clonale ed errore (17, 18)

| Esito | Livello | n | ρ | p | media nulla | p centrato |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| error | patient_means | 10 | 0,467 | 0,175 |  |  |
| error | clone_block_permutation | 277 | 0,220 | 4,0e-04 | 0,0209 | 5,5e-04 |
| error | clone_naive_permutation | 277 | 0,220 | 4,5e-04 |  |  |
| excess | patient_means | 10 | 0,321 | 0,368 |  |  |
| excess | clone_block_permutation | 287 | 0,281 | 5,0e-05 | 0,0265 | 5,0e-05 |
| excess | clone_naive_permutation | 287 | 0,281 | 5,0e-05 |  |  |

Fonte: `clonal_expansion_tests`.

**Sensibilita' del block permutation test (B4)** (`clonal_expansion_perm_sensitivity`)

| Esito | Seme | Permutazioni | p | media nulla | p centrato |
| :--- | :--- | :--- | :--- | :--- | :--- |
| error | 0 | 20.000 | 4,0e-04 | 0,0209 | 5,5e-04 |
| error | 1 | 20.000 | 1,5e-04 | 0,0200 | 2,5e-04 |
| error | 2 | 20.000 | 2,5e-04 | 0,0207 | 4,0e-04 |
| error | 3 | 20.000 | 3,0e-04 | 0,0202 | 4,0e-04 |
| error | 4 | 20.000 | 1,5e-04 | 0,0207 | 3,0e-04 |
| error | 0 | 100.000 | 1,8e-04 | 0,0204 | 3,1e-04 |

- Leave-one-patient-out, errore: p < 0,05 in 10/10 casi, ρ fra 0,195 e 0,252; p centrato < 0,05 in 10/10 (massimo 0,0024) (`loyo_clonal_expansion_error`).
- Leave-one-patient-out, eccesso: p < 0,05 in 10/10 casi, ρ fra 0,236 e 0,368; p centrato < 0,05 in 10/10 (massimo 9,0e-04) (`loyo_clonal_expansion_excess`).
- Leave-one-patient-out, errore, 20.000 permutazioni (5 semi): p fra 5,0e-05 e 0,0015, p centrato fino a 0,0024 (`loyo_clonal_expansion_perm_sensitivity`).
- Leave-one-patient-out, errore, 100.000 permutazioni (1 seme): p fra 1,1e-04 e 0,0013, p centrato fino a 0,0021 (`loyo_clonal_expansion_perm_sensitivity`).

### Quattro ipotesi per paziente (16, 19)

| Ipotesi | Covariata | ρ | p | p_Holm (4 test) |
| :--- | :--- | :--- | :--- | :--- |
| H1 | exhaustion_score | -0,455 | 0,190 | 0,580 |
| H2 | myeloid_fraction | -0,503 | 0,145 | 0,580 |
| H3 | caf_fraction | -0,127 | 0,724 | 1,000 |
| H4 | H4_ni_score | -0,127 | 0,735 | 1,000 |

Fonti: `pseudobulk_hypothesis_tests`, `ni_score_test`, `ni_score_holm_joint`. n = 10 pazienti: potenza bassa.

## CNV e malignita' (13, 14)

AUC LOPO del Random Forest sull'etichetta CNV: 0,836 ± 0,062 (12 pazienti; `malignancy_lopo_summary`); duttali adiacenti predette maligne 37,4% di 334 (`external_test_summary`, non e' una validazione esterna). AUC tumore vs adiacente pooled: score CNV medio 0,598, loci PDAC 0,717 (`cnv_loci_auc`).

| Paziente | Duttali Tumor | Duttali Adjacent | AUC loci PDAC | AUC score medio |
| :--- | :--- | :--- | :--- | :--- |
| PA06 | 2.354 | 26 | 0,608 | 0,531 |
| PA08 | 809 | 210 | 0,576 | 0,494 |
| PA09 | 981 | 14 | 0,526 | 0,798 |
| PA10 | 726 | 83 | 0,556 | 0,640 |
| _mean | 4.870 | 333 | 0,566 | 0,616 |

Fonte: `cnv_loci_auc_by_patient` (AUC entro paziente).

## Simulazioni

**(a) Nadeau-Bengio** (`sim_a_nadeau_bengio`)

| Scenario | n | k | Repliche | p < 0,05 (NB) | p < 0,05 (t non corretto) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| H0 | 90 | 3 | 1.000 | 0,027 | 0,047 |
| H0 | 90 | 5 | 1.000 | 0,018 | 0,061 |
| effect | 120 | 5 | 300 | 0,993 | 1,000 |

**(a) Cluster bootstrap: copertura / frazione di IC che escludono 0** (`sim_a_cluster_bootstrap`)

| Scenario | Metodo | 4 pazienti | 5 pazienti | 10 pazienti | 11 pazienti |
| :--- | :--- | :--- | :--- | :--- | :--- |
| H0_mean0 | expanded_percentile | 0,180 | 0,085 | 0,062 | 0,065 |
| H0_mean0 | percentile | 0,258 | 0,165 | 0,113 | 0,105 |
| coverage_binary | expanded_percentile | 0,830 | 0,925 | 0,935 | 0,948 |
| coverage_binary | percentile | 0,762 | 0,838 | 0,887 | 0,922 |
| coverage_continuous | expanded_percentile | 0,820 | 0,915 | 0,938 | 0,935 |
| coverage_continuous | percentile | 0,743 | 0,835 | 0,887 | 0,895 |
| coverage_pairs_real_sizes_binary | expanded_percentile | ⟦n.d.⟧ | ⟦n.d.⟧ | ⟦n.d.⟧ | 0,927 |
| coverage_pairs_real_sizes_binary | percentile | ⟦n.d.⟧ | ⟦n.d.⟧ | ⟦n.d.⟧ | 0,882 |
| coverage_pairs_real_sizes_continuous | expanded_percentile | ⟦n.d.⟧ | ⟦n.d.⟧ | ⟦n.d.⟧ | 0,938 |
| coverage_pairs_real_sizes_continuous | percentile | ⟦n.d.⟧ | ⟦n.d.⟧ | ⟦n.d.⟧ | 0,890 |
| effect_mean006 | expanded_percentile | 0,570 | 0,550 | 0,810 | 0,855 |
| effect_mean006 | percentile | 0,715 | 0,723 | 0,892 | 0,915 |

**(b) Schemi di riferimento** (`sim_b_reference_schemes`)

| Condizione | Schema | Compartimento | Errore stimato | Errore vero | Distorsione |
| :--- | :--- | :--- | :--- | :--- | :--- |
| H0 | blood | Adjacent_normal | 0,092 | 0,088 | +0,004 |
| H0 | blood | Tumor | 0,089 | 0,084 | +0,005 |
| H0 | loco | Adjacent_normal | 0,084 | 0,084 | +0,000 |
| H0 | loco | PBMC | 0,087 | 0,084 | +0,003 |
| H0 | loco | Tumor | 0,088 | 0,084 | +0,004 |
| H0 | naive | Adjacent_normal | 0,087 | 0,087 | +0,000 |
| H0 | naive | PBMC | 0,084 | 0,084 | +0,0002 |
| H0 | naive | Tumor | 0,084 | 0,084 | +0,0002 |
| artifact_moderate | blood | Adjacent_normal | 0,090 | 0,085 | +0,005 |
| artifact_moderate | blood | Tumor | 0,093 | 0,088 | +0,006 |
| artifact_moderate | loco | Adjacent_normal | 0,084 | 0,084 | +0,000 |
| artifact_moderate | loco | PBMC | 0,119 | 0,086 | +0,034 |
| artifact_moderate | loco | Tumor | 0,093 | 0,088 | +0,005 |
| artifact_moderate | naive | Adjacent_normal | 0,083 | 0,083 | +0,000 |
| artifact_moderate | naive | PBMC | 0,084 | 0,084 | +0,000 |
| artifact_moderate | naive | Tumor | 0,087 | 0,087 | +0,000 |
| artifact_strong | blood | Adjacent_normal | 0,085 | 0,078 | +0,007 |
| artifact_strong | blood | Tumor | 0,086 | 0,081 | +0,006 |
| artifact_strong | loco | Adjacent_normal | 0,084 | 0,083 | +0,001 |
| artifact_strong | loco | PBMC | 0,383 | 0,086 | +0,297 |
| artifact_strong | loco | Tumor | 0,085 | 0,080 | +0,005 |
| artifact_strong | naive | Adjacent_normal | 0,080 | 0,080 | +0,000 |
| artifact_strong | naive | PBMC | 0,086 | 0,085 | +0,0003 |
| artifact_strong | naive | Tumor | 0,084 | 0,084 | +0,0001 |

**(c) Block permutation sotto H0 con confondente di paziente** (`sim_c_block_permutation`)

| Confondente | Repliche | FP block (IC 95%) | FP block centrato (IC 95%) | FP ingenuo | FP nelle prime 100 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 0,0 | 1.000 | 5,0% [3,7; 6,5] | 4,9% [3,6; 6,4] | 5,1% | 4 |
| 0,5 | 1.000 | 5,9% [4,5; 7,5] | 5,1% [3,8; 6,7] | 25,7% | 6 |
| 1,0 | 1.000 | 4,8% [3,6; 6,3] | 4,1% [3,0; 5,5] | 83,1% | 11 |
| 2,0 | 1.000 | 5,1% [3,8; 6,7] | 4,5% [3,3; 6,0] | 99,4% | 5 |

**(d) Stabilita' XAI**: differenza controllo − LOPO +0,143 con effetto paziente (da +0,116 a +0,181 sulle 5 repliche) e +0,030 senza (`sim_d_xai_stability`).


## Esecuzione finale (Fase 3)

Comandi del README nell'ordine; log in `revisione/log/NN_nome.log`, stato e picchi di memoria in `revisione/log/_status.txt` (misurati con `/usr/bin/time -v`; durate e memoria vengono dai log, non da `results/tables/`). Prima esecuzione completa al commit `61fb36b`: dopo 00-04 tre rami, TCR (06-10, 21, 15-19) e 20 + simulazioni in parallelo, poi CNV (13, 14) da solo per la memoria. Dopo le correzioni dei giri 4-8 sono stati rieseguiti 00-04 (giro 8, con il codice finale), 05 (giro 7: solo la legenda delle figure; tabelle identiche byte per byte), 06-10, 13-21 e le simulazioni; 11 ha rigenerato riepilogo, confronti e figure con `--from-table` nel giro 8 (colonna "Ultima esecuzione"); il log di ogni script corrisponde all'ultima esecuzione.

| Script | Esito | Picco di memoria | Durata | Ultima esecuzione |
| :--- | :--- | :--- | :--- | :--- |
| 00_fetch_metadata | rc = 0 | 0,1 GB | 0:00,35 | giro 8c |
| 00b_download_tcr | rc = 0 | 0,0 GB | 0:00,90 | giro 8c |
| 01_probe_labels | rc = 0 | 0,2 GB | 0:01,02 | giro 8c |
| 01b_tcr_check | rc = 0 | 0,2 GB | 0:00,98 | giro 8c |
| 02_build_h5ad_sc | rc = 0 | 7,7 GB | 1:50,14 | giro 8c |
| 02_build_h5ad_sn | rc = 0 | 1,0 GB | 0:11,43 | giro 8c |
| 03_eda_qc_sc | rc = 0 | 0,4 GB | 0:07,19 | giro 8c |
| 03_eda_qc_sn | rc = 0 | 0,3 GB | 0:03,84 | giro 8c |
| 04_splits | rc = 0 | 0,4 GB | 0:04,72 | giro 8c |
| 06_clonotype_consistency | rc = 0 | 0,5 GB | 0:40,66 | giro 8 |
| 07_robustness_clone_ab | rc = 0 | 3,3 GB | 0:37,23 | giro 8 |
| 20_xai_per_patient | rc = 0 | 6,0 GB | 2:00,90 | giro 8 |
| 07_robustness_clone_b | rc = 0 | 3,3 GB | 0:39,38 | giro 8 |
| 08_disambiguate | rc = 0 | 3,5 GB | 0:36,05 | giro 8 |
| 09_arbitrate | rc = 0 | 3,4 GB | 0:45,46 | giro 8 |
| sim_a_calibration | rc = 0 | 0,2 GB | 3:38,20 | giro 8 |
| sim_b_reference_schemes | rc = 0 | 0,3 GB | 1:44,61 | giro 8 |
| 10_loco | rc = 0 | 3,4 GB | 5:43,02 | giro 8 |
| 21_fig_blood_scheme | rc = 0 | 0,3 GB | 0:02,47 | giro 8 |
| 15_inventory | rc = 0 | 0,3 GB | 0:02,11 | giro 8 |
| 16_pseudobulk_patient | rc = 0 | 3,6 GB | 0:15,72 | giro 8 |
| 17_clonal_expansion | rc = 0 | 1,2 GB | 0:09,37 | giro 8 |
| 18_loyo_check | rc = 0 | 1,1 GB | 0:17,31 | giro 8 |
| 19_ni_score | rc = 0 | 0,2 GB | 0:02,10 | giro 8 |
| sim_c_block_permutation | rc = 0 | 0,2 GB | 2:02,95 | giro 8 |
| sim_d_xai_stability | rc = 0 | 0,3 GB | 0:13,50 | giro 8 |
| 13_infercnv | rc = 0 | 6,0 GB | 0:54,77 | giro 8b |
| 14_cnv_loci | rc = 0 | 8,6 GB | 0:19,95 | giro 8b |
| 05_baseline | rc = 0 | 9,7 GB | 17:09,67 | giro 7 |

**Non rieseguiti** (calcolo invariato dall'esecuzione precedente, verificato con `git diff`):
- `12_feature_selection.py` (`52ace44`, 6,1 GB) e `11_model_zoo.py` (`8058701`, processi paralleli in `revisione/log/11_final/`, 5,2-5,8 GB): dopo quei commit cambiano solo controlli sugli argomenti e sugli input (fold con una sola classe, `cell_id` assenti, `nargs`), messaggi, docstring e `help`; per 11 riepilogo, confronti e figure sono stati rigenerati con `--from-table` (stesso codice).
- Le funzioni di `scpipe.py` usate da 05, 11 e 12 (`LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard`, `mean_pairwise_jaccard`) sono identiche a quelle del commit `52ace44` (confronto dell'AST senza docstring); `04_splits.py` e' cambiato solo nei controlli e ricrea `splits_sc.csv` identico byte per byte.

**Scorciatoie usate:** download gia' presenti non ripetuti (00, `download_tcr.sh`, 01, 02); cache di Scrublet (`doublet_scores.csv`, 08), ora riusata solo se contiene esattamente le cellule dell'h5ad (199.184, controllo passato); la cache e' del 10/07 e non e' stata ricalcolata per intero: il revisore della chiusura 1 ha ricalcolato Scrublet con il codice attuale su due librerie (PDAC10-PBMC e PDAC14-Tumor) e ha trovato una differenza massima di circa 1e-16 (`revisione/log/chiusura1_revisore_b.md`, punto 8). Nessun numero del README dipende da questa cache. inferCNV ricalcolato senza cache nel giro di chiusura 1, dopo l'esclusione del cromosoma mitocondriale e dei contig non standard.

**Riproducibilita':** `02_build_h5ad.py` ricrea `sc_raw.h5ad` identico byte per byte e `04_splits.py` ricrea `splits_sc.csv` identico. Nella chiusura 1 sono cambiati i risultati di 13 e 14 (cromosomi), le colonne aggiunte a `clone_switching_summary_clone_b` e `error_rate_by_compartment` e le figure di 11 (etichetta dell'asse); nella chiusura 2 solo le figure `excess_by_pair_{cc}` e `error_rate_by_compartment`, nel giro 6 solo `error_rate_blood` ed `error_rate_loco` (IC non calibrati segnati); nel giro 6 inferCNV e' stato ricalcolato senza cache (impronta con versione e parametri) con risultati identici byte per byte.

Nel giro 7 sono cambiati i risultati del ramo TCR (pareggi dell'etichetta di maggioranza esclusi) e del ramo CNV (soglia del rumore globale, coordinate reali delle finestre), e lo scenario `coverage_pairs_real_sizes` di `sim_a`; vedi la tabella prima/dopo.

Nel giro 8 sono cambiati solo `sim_b_reference_schemes` (stima ed errore vero sulle stesse unita'), `error_decomposition_blood` (riga dei pareggi) e `doublet_mixture_error_nk` (colonne aggiunte). Unici file di `results/` non prodotti dall'ultima versione del codice: le uscite di 12 (eseguito al commit `52ace44`; poi cambiano solo controlli sugli input, `choices`, `help` e la riemissione degli avvisi non di convergenza) e la cache di Scrublet (`doublet_scores.csv`).


## Limiti noti

Problemi del codice che riguardano dati simulati, test sintetici o input anomali (NaN, input vuoti, matrici dense, argomenti non usati dai comandi del README). Sul dataset reale GSE278694 non si presentano; per ciascuno e' indicato come e' stato verificato. Non sono corretti nel loop di revisione (decisione dell'utente).

| Dove | Problema | Verifica sul dataset reale |
| :--- | :--- | :--- |
| `src/scpipe.py:72` (`LogCPM`), `:97` (`HVGSelector`) | richiedono una matrice sparsa CSR: con un `ndarray` denso vanno in errore | tutti gli script convertono prima in CSR (`src/05_baseline.py:309`, `src/11_model_zoo.py:478`, `src/12_feature_selection.py:289`, `src/20_xai_per_patient.py:216`, `src/13_infercnv.py` prima del Random Forest e sulle duttali adiacenti); i test di leakage passano matrici CSR |
| `src/sim/sim_a_calibration.py:88`, `tests/test_stats.py:62` (`SIZES`) | con k = 11 l'undicesima numerosita' (40) e' aggiunta a mano; le coppie reali dell'eccesso (13-282 per paziente) sono piu' sbilanciate | riguarda solo la simulazione di calibrazione; lo scenario `coverage_pairs_real_sizes` di `sim_a_cluster_bootstrap.csv` usa le numerosita' reali delle 11 coppie dell'eccesso (copertura dell'IC espanso 93,8% con esito continuo e 92,8% con esito binario, dentro la banda) |
| `tests/test_stats.py:242`, `:249` | il test di Nadeau-Bengio usa gli stessi semi di `sim_a` (55.000 + k, 51.005): e' un sottoinsieme della stessa simulazione, non una seconda verifica | dichiarato nel docstring del test (righe 39-40); nessun numero del dataset reale dipende dal test |
| `src/11_model_zoo.py:391` (e figure `model_zoo_{task}`), `src/12_feature_selection.py`, `src/20_xai_per_patient.py` | i file di confronto e le figure dipendono solo dal task: un'esecuzione con un altro `--tag` o schema li sovrascrive | i comandi del README eseguono ogni task una sola volta, con uno schema solo; `model_comparison_{task2_v1,task2_hier3,task1a_paired}` coincidono con le tabelle per fold (verificato dai revisori dei giri di chiusura) |
| argparse degli script TCR e CNV | `--help` mostra "(default: None)" per gli argomenti obbligatori | solo testo d'aiuto; i comandi del README passano tutti gli argomenti obbligatori |
| `src/07_robustness.py:119`, `src/08_disambiguate.py:118`, `src/09_arbitrate.py:117`, `src/10_loco.py:173` e `:532` (`TISSUE_MAP[tis]`) | un file VDJ con un tessuto diverso da Tumor, Normal o PBMC manda lo script in `KeyError` (06 invece lo salta con un avviso) | i 30 file di `vdj_inventory.csv` hanno solo i tessuti Tumor, Normal e PBMC; tutti gli script terminano con rc = 0 (`revisione/log/_status.txt`) |
| `src/20_xai_per_patient.py:135` (`np.argsort(-coef)[:TOP_K]`) | con meno di 50 coefficienti non nulli i geni mancanti dei top-50 sarebbero scelti fra pareggi a zero in ordine arbitrario | sul dataset reale i fold LOPO hanno da 3.887 a 4.380 coefficienti non nulli (`revisione/log/20_xai_per_patient.log`) |
| scelta del contig dominante (`sort_values("umis").drop_duplicates` in 06-10) | con due catene TRB a pari UMI massimi il clone assegnato dipende dall'ordine delle righe (ordinamento non stabile) | riguarda 346 cellule su 123.280 per la TRB; un ricalcolo indipendente con una scelta diversa del contig da' gli stessi numeri principali (Tumor 277 cloni, 0,170; Adjacent_normal 117, 0,060) e cambia solo la sensibilita' "pareggi = errore" (0,199 -> 0,201): `revisione/log/giro8_revisore_b.md`, punto 7 |


## Non risolti

Problemi segnalati dai revisori dell'ultimo giro (giro 8) e non corretti, con il motivo. Nessuno cambia un numero principale sul dataset reale.

| Dove | Problema | Motivo |
| :--- | :--- | :--- |
| `src/06_clonotype_consistency.py` (riepilogo split-half), `src/16_pseudobulk_patient.py`, `src/17_clonal_expansion.py`, `src/18_loyo_check.py`, `src/19_ni_score.py` | i pareggi dell'etichetta di maggioranza sono esclusi, ma la sensibilita' "pareggi = errore / pareggi = corretti" e' tabulata solo in 09 e 10 (`arbitration_tie_sensitivity`, `error_rate_tie_sensitivity`) | il revisore l'ha ricalcolata e le conclusioni non cambiano (17, errore: ρ 0,186 e 0,217 contro 0,220; H1-H4 con segno e p simili; 06: +0,052 e +0,058 contro +0,063): `revisione/log/giro8_revisore_b.md`, punto 6. Questi valori non sono in `results/tables` e non sono citati come risultati |
| `README.md:94` | la radice indicata e' `progetto-ml/`, la cartella locale si chiama `pdac-ml` | `progetto-ml` e' il nome del repository remoto, richiesto dalla correzione D3 |


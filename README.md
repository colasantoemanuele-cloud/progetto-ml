# Machine learning su scRNA-seq e validazione clonale TCR nel carcinoma duttale pancreatico (PDAC)

Studio computazionale e metodologico sulla SuperSeries NCBI GEO **`GSE278694`** (Chen M-M et al., *Cancer Cell* 2025): **199.184 cellule**, **14 pazienti**, **30.543 geni**, **29 librerie scRNA-seq** (`GSE278688`) e **30 librerie scTCR-seq V(D)J** (`GSE300435`), in tre compartimenti anatomici (`Tumor`: 71.414 cellule, 13 pazienti; `Adjacent_normal`: 38.264 cellule, 5 pazienti; `PBMC`: 89.506 cellule, 11 pazienti). La coorte snRNA-seq (`GSE278689`, **32.005 nuclei**, 8 pazienti) e' esclusa dalle analisi supervisionate: i suoi pazienti sono diversi da quelli scRNA-seq, quindi la modalita' tecnica coinciderebbe con il paziente.

Fonti dei numeri di questo paragrafo: `results/tables/dataset_summary_sc.csv`, `dataset_summary_sn.csv`, `vdj_inventory.csv` (vedi [Struttura dei risultati](#5-struttura-dei-risultati)).

---

## 1. I tre problemi di studio e la sintesi dei risultati

Il lavoro esamina tre problemi di classificazione supervisionata sul PDAC a singola cellula. In tutti e tre la stima "di pratica corrente" e' piu' alta di quella che tiene conto della struttura dei dati, per tre ragioni diverse riconducibili allo stesso principio: **l'unita' statistica indipendente e' il paziente, non la cellula**.

| Problema | Popolazione e dimensione | Stima di pratica corrente | Stima che rispetta la struttura di gruppo | Origine dell'artefatto |
| :--- | :--- | :--- | :--- | :--- |
| **1. Tipo cellulare** (`task2_v1`, 9 classi) | solo tessuto tumorale: 70.547 cellule, 13 pazienti | split casuale: macro-F1 **0,964** (logistica L2) | split per paziente: macro-F1 **0,954** (logistica L2, 5 fold per paziente); LightGBM 0,962 (LOPO) | **Circolarita'**: le etichette di riferimento derivano da clustering e marcatori sulla stessa matrice usata per addestrare i classificatori. L'accuratezza e' reale ma misura la riproduzione dell'annotazione, non una verita' indipendente. |
| **2. Compartimento** (`task1a_paired`, Tumor vs Adjacent_normal) | 67.662 cellule, 5 pazienti con entrambi i tessuti | split casuale: macro-F1 **0,875 ± 0,002** | LOPO: macro-F1 **0,729 ± 0,100** (Δ = −0,146; deviazione standard fra fold 50 volte piu' grande) | **Confondimento paziente ↔ tessuto**: lo split casuale mette cellule dello stesso paziente in training e test, e il modello sfrutta la firma del paziente e la composizione cellulare. |
| **3. Malignita' delle duttali** (`13_infercnv.py`, `14_cnv_loci.py`) | duttali tumorali di 12 pazienti (LOPO); 334 duttali adiacenti | AUC LOPO **0,829 ± 0,069** del Random Forest sull'etichetta CNV derivata dal trascrittoma; il 37,1% delle duttali adiacenti e' predetto maligno | sui loci PDAC l'AUC tumore-vs-adiacente pooled e' 0,682 ma scende a **≈ 0,54 entro paziente** (0,51-0,58 nei 4 pazienti valutabili); lo score CNV medio ha AUC pooled 0,492 | **Etichetta derivata dal trascrittoma e confondimento con il paziente**: l'etichetta di malignita' viene da inferCNV sugli stessi dati; il contrasto tumore-adiacente pooled mescola differenze fra pazienti. |

Note alla tabella:
- Riga 1: 0,954 e 0,964 vengono da `05_baseline.py` (5 fold per paziente e 5 fold casuali); 0,962 da `11_model_zoo.py` in leave-one-patient-out.
- Riga 3: in inferCNV il riferimento diploide sono le cellule immunitarie di **tutti** i pazienti (una sola esecuzione); solo la soglia (mediana + 3·MAD degli immuni) e' calcolata per paziente, ed e' stimata sulle stesse cellule su cui inferCNV centra lo score. Per questo anche il 45-62% delle duttali adiacenti dei pazienti con almeno 20 duttali adiacenti risulta "maligno" (`malignancy_labels.csv`). La sovrapposizione fra duttali tumorali e adiacenti e' **compatibile** con un effetto di campo di cancerizzazione, ma con questi limiti (riferimento comune, soglia sulle stesse cellule, confini dei loci approssimati, 4 pazienti valutabili) e' un'interpretazione, non una dimostrazione.

### Confronto dei modelli, selezione delle feature e spiegazioni per paziente

* **13 classificatori (`11_model_zoo.py`, LOPO, 20.000 cellule di training per fold)**. Su `task2_v1` la separazione e' quasi lineare (logistica L2 0,943, LightGBM 0,962). Su `task2_hier3` (CD4T / CD8T / NK, 42.168 cellule) la logistica L2 ottiene 0,851 e i modelli non lineari di piu': LightGBM 0,923 (+7,2 punti), SVM-RBF 0,906 (+5,6), MLP 0,903 (+5,3), Random Forest 0,901 (+5,1). Su `task1a_paired` i modelli vanno da 0,496 (dummy) a 0,725 (SVM-RBF); escluso il dummy, nessuna differenza dal migliore e' significativa al t corretto di Nadeau-Bengio (p ≥ 0,080; dummy p = 0,072). Con 5 fold la potenza e' bassa: un p non significativo **non** dimostra che i modelli siano equivalenti.
* **Selezione delle feature (`12_feature_selection.py`)**. Tutti i selettori di un fold usano le stesse 20.000 cellule di training (sottocampione casuale uniforme del solo training). A k = 50 la logistica L1 ottiene 0,936 su `task2_v1` e 0,712 su `task1a_paired`, il Random Forest 0,926 e 0,735, i geni piu' dispersi (HVG) 0,659 e 0,623 (k scelti a caso: 0,367 e 0,589). Su `task2_v1` la logistica L1 non converge in 400 iterazioni in nessun fold e ha 524-583 coefficienti non nulli, quindi a k = 1000 una parte dei geni scelti ha coefficiente nullo (`feature_selection_diag_task2_v1.csv`); su `task1a_paired` converge (177-268 iterazioni). La stabilita' (Jaccard medio fra fold: i 5 fold per paziente di `fold_group_5` su `task2_v1`, i 5 fold LOPO su `task1a_paired`) e' piu' bassa per L1 (0,480 a k = 200 su `task2_v1`, 0,254 a k = 1000 su `task1a_paired`) che per i filtri univariati (chi-quadro 0,842 e 0,678; ANOVA 0,794 e 0,665; HVG 0,793 e 0,771). Su `task1a_paired`, a k = 50, i geni scelti in tutti e 5 i fold sono 18 per HVG, 19 per chi-quadro, 19 per ANOVA, 16 per il Random Forest e 4 per L1 (`feature_core_genes_task1a_paired.csv`).
* **Spiegazioni per paziente (`20_xai_per_patient.py`)**. Sui 5 fold LOPO di `task1a_paired`, il Jaccard medio fra i top-50 geni di una logistica L1 e' **0,300** (0,220-0,408 fra coppie di fold; Spearman medio fra i vettori dei coefficienti 0,234); 9 geni sono nei top-50 di tutti i fold: BTG2, CD69, CREM, CTRB1, DUSP4, FOS, GAPDH, HSPH1, VIM. Il controllo che toglie lo stesso numero di cellule a caso, senza struttura per paziente, da' **0,459** (differenza +0,159). Togliere un paziente cambia la regola appresa piu' che togliere lo stesso numero di cellule a caso, come atteso se le cellule dello stesso paziente sono ridondanti. Due cautele: in simulazione (`src/sim/sim_d_xai_stability.py`, 5 pazienti, 5 repliche per scenario) la differenza era in media +0,143 (da +0,116 a +0,181) con un effetto paziente costruito e +0,030 (da −0,013 a +0,052) senza; il valore osservato, +0,159, cade nell'intervallo dello scenario con effetto paziente, ma la simulazione ha una sola forma di effetto e pochi pazienti, quindi e' un'indicazione e non una calibrazione; dei 9 geni comuni, FOS, CD69, BTG2, HSPH1, DUSP4 e CREM sono geni di risposta immediata associati allo stress di dissociazione e CTRB1 e' un marcatore acinare (probabile RNA ambientale), non necessariamente biologia tumore/normale.

### Validazione clonale con il TCR

Il riarrangiamento V(D)J delle catene α/β precede l'espansione periferica e non dipende dallo stato trascrizionale: la sequenza CDR3 identifica un clone indipendentemente dal clustering, e permette di misurare l'affidabilita' dell'annotazione.

* **Rumore di base ed eccesso fra compartimenti (`07_robustness.py`)**. Due cellule dello stesso clone nello stesso compartimento hanno etichette diverse con probabilita' **0,107** (rumore di annotazione di base, stimatore senza reinserimento). Fra compartimenti diversi la discordanza e' piu' alta: **eccesso +0,059** (IC 95% cluster bootstrap sui pazienti **+0,031; +0,110**, 11 pazienti, 1.114 coppie clone-compartimenti), simile con la sola catena TRB e con la coppia TRA+TRB (+0,065, IC +0,035; +0,121). Con pesi uguali per i due compartimenti (½‖p_A − p_B‖², stesso bootstrap) l'eccesso e' +0,074 (IC +0,050; +0,119). La versione esplorativa split-half di `06` da' +0,057 (IC cluster bootstrap +0,017; +0,127, 14 pazienti).
* **Tasso d'errore con riferimento solo sangue (`10_loco.py`, `21_fig_blood_scheme.py`)**. Con l'identita' CD4/CD8 del clone stimata dai marcatori del solo PBMC, l'etichetta e' sbagliata nel **19,5%** dei cloni nel `Tumor` (**0,195**, IC 95% **0,110; 0,290**, 287 cloni, 10 pazienti) e nel **5,9%** nell'`Adjacent_normal` (0,059, IC 0,000; 0,093, 118 cloni, 4 pazienti). La differenza appaiata sugli stessi cloni e' **+0,080** (IC **0,000; +0,192**, 88 cloni, 4 pazienti): con 4 pazienti esistono solo 35 ricampionature distinte e l'IC non e' calibrato; il 99,7% delle repliche bootstrap e' > 0. Scomposizione dell'errore nel Tumor: 38 cloni CD4 nel sangue etichettati CD8T (0,132), 18 cloni etichettati NK (0,063: 12 CD8 e 6 CD4 nel sangue), 0 cloni CD8 etichettati CD4T; escludendo le etichette NK l'errore e' 0,141 (IC 0,058; 0,261). Nelle cellule etichettate CD8T del tumore, CD8A e' rilevato nel **7%** delle cellule dei cloni discordanti contro l'**88%** dei concordanti (`08_disambiguate.py`), e CD4 medio e' 16 volte piu' alto (rapporto fra medie di log1p(CP10k), non un fold-change). I dati sono compatibili con cellule CD4 etichettate CD8T; resta aperta l'alternativa di un fenotipo CD4 citotossico (GNLY 1,78 nei discordanti contro 0,44 nei concordanti). Con lo schema che esclude il compartimento giudicato (LOCO) la differenza Tumor − Adjacent_normal e' +0,027 (IC −0,010; +0,124): l'asimmetria dipende dallo schema di riferimento.
* **Doppietti fra le NK con TCR (`09_arbitrate.py`, `10_loco.py`)**. Fra le cellule etichettate come non-T con un TCR produttivo (NK 5,9%, Endothelial 5,6%, Stellate 5,5%, Platelet 13,8% delle cellule della classe), due misure fisiche che non usano il clustering stimano che **72-79% delle NK-TCR+ siano doppietti T+NK**: dall'eccesso di UMI totali d1 = **0,723** (IC 0,556; 0,848) e dalla diluizione degli UMI del TCR d2 = **0,790** (IC 0,660; 0,876). Le due stime non sono indipendenti (usano le stesse mediane di n_counts) e sono pooled fra tessuti: per compartimento vanno da 0,51 (Tumor) a 0,84 (PBMC) per d1; con le medie al posto delle mediane d1 = 0,682 e d2 = 0,840.
* **Espansione clonale ed errore (`17_clonal_expansion.py`, `18_loyo_check.py`)**. Su 287 cloni condivisi fra sangue e tumore (10 pazienti) la dimensione clonale e' **associata all'errore entro paziente**: errore binario ρ = +0,177 (p = 0,0009, block permutation test entro paziente con 20.000 permutazioni Monte Carlo), eccesso continuo ρ = +0,281 (p = 0,00005, il minimo raggiungibile con 20.000 permutazioni). Il p e' condizionato al paziente ma il ρ e' pooled: la distribuzione nulla ha media +0,012 (errore) e +0,027 (eccesso), e il p sulla statistica centrata e' 0,0015 e 0,00005. Con 5 semi il p dell'errore va da 0,0009 a 0,0018; con 100.000 permutazioni e' 0,0017. Fra le medie per paziente il segno e' lo stesso ma l'associazione non e' significativa (ρ = +0,467, p = 0,18, n = 10). Togliendo un paziente alla volta (leave-one-patient-out) l'associazione resta con p < 0,05 in **10/10** casi per entrambi gli esiti (ρ fra 0,146 e 0,197 per l'errore; fra 0,236 e 0,368 per l'eccesso).
* **Quattro covariate per paziente (`16_pseudobulk_patient.py`, `19_ni_score.py`)**. Per spiegare le differenze di errore fra pazienti (n = 10) senza cercare fra migliaia di geni, sono state testate quattro ipotesi, formulate dopo i moduli 08-10 ma prima di calcolare le covariate (H3 indicata fin dall'inizio come speculativa), con Spearman, permutazione Monte Carlo (20.000) e Holm-Bonferroni sulla famiglia di 4 test. Nessuna e' significativa e tutte hanno segno opposto a quello ipotizzato:
  1. `H1` esaurimento delle CD8T tumorali: ρ = −0,47, p = 0,176, p_Holm = 0,53. Limiti: la deriva CD8 → CD4/NK ipotizzata va nel senso opposto agli errori osservati (soprattutto CD4 → CD8T), e il punteggio e' calcolato sulle cellule etichettate CD8T, cioe' sulle stesse etichette il cui errore e' l'esito;
  2. `H2` frazione mieloide nel tumore: ρ = −0,54, p = 0,113, p_Holm = 0,45;
  3. `H3` frazione di CAF nel tumore: ρ = −0,15, p = 0,673, p_Holm = 1,00;
  4. `H4` NI score clinico (Table S1A di Chen et al.): ρ = −0,14, p = 0,705, p_Holm = 1,00.

  Con n = 10 questi test hanno poca potenza: "non significativo" non esclude un effetto.

---

## 2. Moduli sorgente (`src/`, `src/sim/`, `tests/`)

Gli script lavorano su matrici sparse `scipy.sparse.csr_matrix` `float32`; le trasformazioni stimate dai dati sono dentro `Pipeline` di scikit-learn o funzioni che ricevono il solo training (verificato in `tests/test_leakage.py`). Ogni script ha una docstring con scopo, input, output, argomenti e uso, e risponde a `--help`.

| Modulo | Funzione |
| :--- | :--- |
| [`src/00_fetch_metadata.py`](src/00_fetch_metadata.py) | Scarica l'archivio MINiML e `filelist.txt` di `GSE278694`, classifica i campioni per modalita' e paziente e verifica la fattibilita' del disegno. |
| [`src/download_tcr.sh`](src/download_tcr.sh) | Scarica i 30 file VDJ di `GSE300435` elencati in `filelist.txt`, con retry e controllo `gzip`; esce con errore se i file non sono 30. |
| [`src/01_probe_labels.py`](src/01_probe_labels.py) | Scarica metadati, barcode e un file VDJ pilota (PA01 Tumor) e verifica etichette e sovrapposizione dei barcode. |
| [`src/01b_tcr_check.py`](src/01b_tcr_check.py) | Sul file VDJ pilota: intersezione dei barcode TCR e GEX dentro la stessa libreria, contro l'atteso per caso. |
| [`src/02_build_h5ad.py`](src/02_build_h5ad.py) | Conversione a blocchi (20 milioni di righe) da MatrixMarket (`.mtx.gz`) ad AnnData CSR `float32`. |
| [`src/03_eda_qc.py`](src/03_eda_qc.py) | QC descrittivo per libreria, composizione e pazienti per classe (nessun filtro di QC). |
| [`src/04_splits.py`](src/04_splits.py) | Task e partizioni: 5 fold per paziente, leave-one-patient-out, 5 fold casuali (controllo); identifica i 5 pazienti appaiati. |
| [`src/05_baseline.py`](src/05_baseline.py) | Dummy e logistica L2: divario fra split casuale e split per paziente. |
| [`src/06_clonotype_consistency.py`](src/06_clonotype_consistency.py) | Discordanza delle etichette nei cloni (split-half, IC cluster bootstrap) e stima c_hat della contaminazione T nelle NK. |
| [`src/07_robustness.py`](src/07_robustness.py) | Stimatore a coppie non distorto della discordanza, cluster bootstrap sui pazienti, RNA ambientale e marcatori. |
| [`src/08_disambiguate.py`](src/08_disambiguate.py) | Cloni presenti nei tre compartimenti, marcatori nei cloni discordanti e concordanti, Scrublet. |
| [`src/09_arbitrate.py`](src/09_arbitrate.py) | Identita' CD4/CD8 del clone dai marcatori, errore per compartimento (schema circolare), misure fisiche dei doppietti. |
| [`src/10_loco.py`](src/10_loco.py) | Schemi di riferimento naive, leave-one-compartment-out e solo sangue; sensibilita'; frazione di doppietti NK-TCR+. |
| [`src/21_fig_blood_scheme.py`](src/21_fig_blood_scheme.py) | Tabella e figura dello schema solo sangue (stessa funzione `report()` di 10). |
| [`src/11_model_zoo.py`](src/11_model_zoo.py) | 13 classificatori in LOPO a parita' di cellule di training, Nadeau-Bengio e Wilcoxon. |
| [`src/12_feature_selection.py`](src/12_feature_selection.py) | 6 criteri di selezione dei geni a k = 50, 200, 1000: prestazione e Jaccard fra fold. |
| [`src/20_xai_per_patient.py`](src/20_xai_per_patient.py) | Stabilita' delle spiegazioni L1 fra fold LOPO contro un controllo senza struttura per paziente. |
| [`src/13_infercnv.py`](src/13_infercnv.py) | Score CNV (inferCNV, riferimento: immunitarie di tutti i pazienti), etichetta di malignita' con soglia per paziente, Random Forest LOPO, applicazione alle duttali adiacenti. |
| [`src/14_cnv_loci.py`](src/14_cnv_loci.py) | Score CNV sui loci ricorrenti del PDAC; AUC tumore vs adiacente pooled ed entro paziente. |
| [`src/15_inventory.py`](src/15_inventory.py) | Inventario delle colonne di `obs` e dei file VDJ, con confondimento rispetto al paziente. |
| [`src/16_pseudobulk_patient.py`](src/16_pseudobulk_patient.py) | Ipotesi H1-H3 per paziente, permutazione Monte Carlo e Holm. |
| [`src/17_clonal_expansion.py`](src/17_clonal_expansion.py) | Espansione clonale contro errore: livello paziente, block permutation entro paziente, sensibilita' a semi e permutazioni. |
| [`src/18_loyo_check.py`](src/18_loyo_check.py) | Sensibilita' leave-one-patient-out di 17 (i file si chiamano `loyo_*` per ragioni storiche). |
| [`src/19_ni_score.py`](src/19_ni_score.py) | Ipotesi H4 (NI score dalla Table S1A) e Holm sulla famiglia H1-H4. |
| [`src/scpipe.py`](src/scpipe.py) | Modulo condiviso: `LogCPM` (log1p di CP10k), `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard`, cluster bootstrap, block permutation, permutazione, Holm, regola della maggioranza. |
| [`src/estrai_conteggi.py`](src/estrai_conteggi.py) | Utilita': cellule per paziente della coorte appaiata (67.662 in totale). |
| [`src/sim/`](src/sim/) | Simulazioni (a)-(d), vedi [§6](#6-simulazioni-e-test). |
| [`tests/`](tests/) | Test pytest: formule, anti-leakage, calibrazione statistica, smoke test. |

---

## 3. Ambiente, dati e citazione

**Ambiente.** Definito in [`envs/environment.yml`](envs/environment.yml) (ambiente `pdac`, Python 3.11, versioni fissate a quelle che hanno prodotto `results/`) e bloccato con checksum in [`envs/conda-lock.yml`](envs/conda-lock.yml). Tra i pacchetti: `scanpy`, `anndata`, `scikit-learn`, `lightgbm`, `infercnvpy`, `pybiomart`, `openpyxl`, `scipy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pytest`, `ruff`.

```bash
mamba env create -f envs/environment.yml      # oppure: conda-lock install -n pdac envs/conda-lock.yml
mamba activate pdac
```

**Dati.** Tutti i comandi si eseguono dalla cartella radice del repository, `progetto-ml/`.
- Matrici scRNA-seq e snRNA-seq: `src/02_build_h5ad.py` scarica da GEO i file supplementari di `GSE278688` (scRNA) e `GSE278689` (snRNA) in `data/raw/GSE278688/` e `data/raw/GSE278689/` (matrix, barcodes, features, metadata; ~1,2 GB). Sono le SubSeries della SuperSeries `GSE278694`: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278694>.
- File VDJ: `src/00_fetch_metadata.py` scarica `filelist.txt` di `GSE278694`; `src/download_tcr.sh` scarica i 30 file `*_VDJ_filtered_contig_annotations.csv.gz` dei campioni di `GSE300435` in `data/raw/tcr/`.
- Metadati clinici: `mmc2.xlsx` e' la Table S1 del materiale supplementare di Chen et al. 2025 (pagina dell'articolo su *Cancer Cell*, doi sotto, sezione "Supplemental information"). Va scaricato a mano e messo in `data/raw/paper_supplement/mmc2.xlsx`; `19_ni_score.py` ne legge il foglio `Table S1A`.

**Citazione del dataset.** Chen M-M et al. *Cancer Cell* 2025;43(9):1656-1676.e10. doi:[10.1016/j.ccell.2025.06.020](https://doi.org/10.1016/j.ccell.2025.06.020).

**Casualita'.** Tutti gli script usano semi fissi (`SEED = 0`; nei sottocampionamenti per fold di 11 e 12, `SEED + fold`). Sottocampionamento: casuale uniforme, solo sul training, non stratificato (20.000 cellule per fold in 11, 12 e 20). Cluster bootstrap sui pazienti: 2.000 repliche (06, 07, 08, 09, 10, 21), 500 per d1 e d2 (10), 2.000 per gli IC di Spearman (16, 17, 19); intervallo "expanded percentile" (vedi [§6](#6-simulazioni-e-test)). Permutazioni Monte Carlo: 20.000 (16, 17, 19), 10.000 (18), 100.000 nella sensibilita' di 17.

---

## 4. Esecuzione

I comandi, eseguiti in quest'ordine dalla radice del repository, ricreano tutte le tabelle e le figure di `results/`.

```bash
# 1. Metadati GEO, file VDJ, verifiche preliminari
python src/00_fetch_metadata.py --outdir data/raw/metadata
bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr
python src/01_probe_labels.py --outdir data/raw/probe
python src/01b_tcr_check.py --probe data/raw/probe

# 2. AnnData, QC descrittivo, task e split
python src/02_build_h5ad.py --which sc --outdir data --no-counts-layer
python src/02_build_h5ad.py --which sn --outdir data --no-counts-layer
python src/03_eda_qc.py --h5ad data/interim/sc_raw.h5ad --tag sc
python src/03_eda_qc.py --h5ad data/interim/sn_raw.h5ad --tag sn
python src/04_splits.py --h5ad data/interim/sc_raw.h5ad --tag sc

# 3. Modelli: divario casuale/paziente, 13 modelli, selezione delle feature, XAI
H=data/interim/sc_raw.h5ad; S=data/interim/splits_sc.csv
python src/05_baseline.py --h5ad $H --splits $S \
    --tasks task2_hier1 task2_hier3 task2_v1 task2_v2 task1a_paired task1a_all
python src/11_model_zoo.py --h5ad $H --splits $S --tasks task2_v1 task2_hier3 --scheme fold_lopo
python src/11_model_zoo.py --h5ad $H --splits $S --tasks task1a_paired --tag task1 --scheme fold_lopo
python src/12_feature_selection.py --h5ad $H --splits $S --task task2_v1
python src/12_feature_selection.py --h5ad $H --splits $S --task task1a_paired --scheme fold_lopo
python src/20_xai_per_patient.py --h5ad $H --splits $S --task task1a_paired

# 4. Validazione clonale TCR
python src/06_clonotype_consistency.py --h5ad $H --filelist data/raw/metadata/filelist.txt \
    --outdir data/raw/tcr --skip-download
python src/07_robustness.py --h5ad $H --tcrdir data/raw/tcr --clone-def clone_ab
python src/07_robustness.py --h5ad $H --tcrdir data/raw/tcr --clone-def clone_b
python src/08_disambiguate.py --h5ad $H --tcrdir data/raw/tcr
python src/09_arbitrate.py --h5ad $H --tcrdir data/raw/tcr
python src/10_loco.py --h5ad $H --tcrdir data/raw/tcr
python src/21_fig_blood_scheme.py

# 5. CNV e malignita'
python src/13_infercnv.py --h5ad $H
python src/14_cnv_loci.py --h5ad $H

# 6. Inventario, ipotesi per paziente, espansione clonale
python src/15_inventory.py --h5ad $H --tcrdir data/raw/tcr
python src/16_pseudobulk_patient.py --h5ad $H --arbitration results/tables/arbitration_blood.csv
python src/17_clonal_expansion.py --arbitration results/tables/arbitration_blood.csv \
    --pairwise results/tables/pairwise_discordance_clone_b.csv
python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
python src/19_ni_score.py --table-s1 data/raw/paper_supplement/mmc2.xlsx \
    --arbitration results/tables/arbitration_blood.csv \
    --previous-tests results/tables/pseudobulk_hypothesis_tests.csv   # dopo 16

# 7. Simulazioni e test
python src/sim/sim_a_calibration.py
python src/sim/sim_b_reference_schemes.py
python src/sim/sim_c_block_permutation.py
python src/sim/sim_d_xai_stability.py
pytest                      # tutti i test (circa 3 minuti); solo i veloci: pytest -m "not slow"
ruff check src tests
```

Note: `07` va eseguito con `clone_b` dopo `clone_ab` perche' `08` e `17` leggono `pairwise_discordance_clone_b.csv`; `13` e `14` riusano le cache `data/interim/cnv_scores.csv` e `cnv_matrix.npz` se presenti (vanno cancellate per ricalcolare inferCNV); `08` riusa `results/tables/doublet_scores.csv` (Scrublet) se presente. `11` impiega circa un'ora per task.

---

## 5. Struttura dei risultati

| Script | Tabelle in `results/tables/` |
| :--- | :--- |
| 03 | `dataset_summary_{sc,sn}`, `qc_by_library_{sc,sn}`, `composition_{sc,sn}`, `patient_support_{sc,sn}` |
| 04 | `task_summary_sc` (e `data/interim/splits_sc.csv`) |
| 05 | `baseline_sc` (per fold), `baseline_summary_sc` (media e sd) |
| 06 | `nk_contamination`, `tcr_rate_by_celltype`, `clone_switching_across_clone_b`, `clone_switching_summary_clone_b` |
| 07 | `pairwise_discordance_{clone_b,clone_ab}`, `excess_summary_{cc}`, `excess_by_compartment_pair_{cc}`, `discordant_clone_markers_{cc}`, `ambient_discriminant`, `chat_sensitivity_{raw,hc}` |
| 08 | `triple_clones`, `triple_clones_summary`, `markers_disc_vs_conc`, `marker_test_disc_vs_conc`, `doublet_tcr_pos_vs_neg` |
| 09 | `arbitration`, `error_rate_by_compartment`, `doublet_physical` |
| 10 | `arbitration_{naive,loco,blood}`, `error_rate_schemes`, `error_rate_paired_diff`, `error_decomposition_blood`, `error_rate_identity_rules`, `ambient_cd4_threshold`, `doublet_mixture`, `doublet_mixture_boot` |
| 21 | `error_rate_blood` |
| 11 | `model_zoo_{sc,task1}` (per fold), `model_zoo_summary_{sc,task1}`, `model_comparison_{task2_v1,task2_hier3,task1a_paired}` |
| 12 | `feature_selection_{task}`, `feature_selection_summary_{task}`, `feature_stability_{task}`, `feature_selection_diag_{task}`, `feature_core_genes_{task}` |
| 20 | `xai_jaccard_task1a_paired`, `xai_spearman_task1a_paired`, `xai_summary_task1a_paired`, `xai_common_genes_task1a_paired` |
| 13 | `malignancy_labels`, `malignancy_lopo`, `malignancy_lopo_summary`, `external_test_predictions`, `external_test_summary` |
| 14 | `cnv_loci_auc`, `cnv_loci_auc_by_patient` |
| 15 | `obs_inventory`, `vdj_inventory` |
| 16 | `pseudobulk_patient_table`, `pseudobulk_hypothesis_tests` |
| 17 | `clonal_expansion_joined`, `clonal_expansion_tests`, `clonal_expansion_perm_sensitivity` |
| 18 | `loyo_clonal_expansion_{error,excess}` |
| 19 | `ni_score_joined`, `ni_score_test`, `ni_score_holm_joint` |
| sim | `sim_a_nadeau_bengio`, `sim_a_cluster_bootstrap`, `sim_b_reference_schemes`, `sim_c_block_permutation`, `sim_d_xai_stability` |

Da dove viene ogni numero del README:

| Numero | Tabella (colonna / riga) |
| :--- | :--- |
| 199.184 cellule, 14 pazienti, 30.543 geni, 29 librerie, cellule e pazienti per tessuto | `dataset_summary_sc` |
| 32.005 nuclei, 8 pazienti | `dataset_summary_sn` |
| 30 librerie VDJ | `vdj_inventory` (30 righe) |
| 70.547 / 67.662 / 42.168 cellule, 13 / 5 pazienti dei task | `task_summary_sc` |
| 0,954; 0,964; 0,875 ± 0,002; 0,729 ± 0,100 | `baseline_summary_sc` (`logreg_l2`, `macro_f1`) |
| Δ = −0,146 (segno opposto di `gap_random_minus_honest`), deviazione standard 50 volte piu' grande | `baseline_gap_sc` (`gap_random_minus_honest`, `sd_ratio_honest_over_random`, riga `task1a_paired`) |
| 0,943, 0,962, 0,851, 0,923, 0,906, 0,903, 0,901, 0,496-0,725; +7,2, +5,6, +5,3, +5,1 punti | `model_zoo_summary_sc`, `model_zoo_summary_task1` (`macro_f1_mean`, `delta_vs_logreg_l2`) |
| p ≥ 0,080, p = 0,072 | `model_comparison_task1a_paired` (`p_nadeau_bengio`) |
| macro-F1 e Jaccard della selezione delle feature | `feature_selection_summary_{task}` (`macro_f1_mean`), `feature_stability_{task}` |
| convergenza e coefficienti non nulli di L1; geni comuni ai 5 fold | `feature_selection_diag_{task}`, `feature_core_genes_task1a_paired` |
| 0,300; 0,459; +0,159; 0,234; 9 geni | `xai_summary_task1a_paired`, `xai_common_genes_task1a_paired`; 0,220-0,408 da `xai_jaccard_task1a_paired` |
| 0,107; +0,059 [+0,031; +0,110]; +0,074; 11 pazienti; 1.114 coppie | `excess_summary_clone_b`; +0,065 da `excess_summary_clone_ab` |
| +0,057 [+0,017; +0,127] | `clone_switching_summary_clone_b` |
| 0,195 [0,110; 0,290]; 0,059 [0,000; 0,093]; 35 ricampionature | `error_rate_blood` |
| +0,080 [0,000; +0,192]; 99,7%; +0,027 [−0,010; +0,124] | `error_rate_paired_diff` |
| 38, 18, 12, 6, 0 cloni; 0,141 [0,058; 0,261] | `error_decomposition_blood` |
| 7% / 88%; 16×; GNLY 1,78 / 0,44 | `marker_test_disc_vs_conc` (riga Tumor) |
| NK 5,9%, Endothelial 5,6%, Stellate 5,5%, Platelet 13,8% | `tcr_rate_by_celltype` |
| d1 0,723, d2 0,790 e varianti | `doublet_mixture` |
| ρ, p, media nulla, p centrato, 287 cloni | `clonal_expansion_tests`, `clonal_expansion_perm_sensitivity` |
| 10/10, intervalli di ρ | `loyo_clonal_expansion_{error,excess}` |
| H1-H3 | `pseudobulk_hypothesis_tests`; H4 `ni_score_test`; p di Holm `ni_score_holm_joint` |
| 0,829 ± 0,069; 12 pazienti | `malignancy_lopo_summary` |
| 37,1%; 334 duttali | `external_test_summary` (riga `_all`) |
| 45-62% | `malignancy_labels` (Adjacent_normal, pazienti con n ≥ 20) |
| 0,492; 0,682 | `cnv_loci_auc` |
| ≈ 0,54; 0,51-0,58; 4 pazienti | `cnv_loci_auc_by_patient` |
| numeri delle simulazioni | `sim_*` |

---

## 6. Simulazioni e test

**Simulazioni (`src/sim/`)**, eseguite con semi fissi; i risultati sono in `results/tables/sim_*.csv`.
- (a) `sim_a_calibration.py`: t di Nadeau-Bengio sotto H0 (due modelli dummy) e con effetto noto (disegno ripreso da `genomic-audit-mvp`, stesso autore); cluster bootstrap con 4, 5, 10 e 11 pazienti, con l'intervallo percentile semplice e con quello "expanded percentile". Risultati: sotto H0 Nadeau-Bengio da' p < 0,05 nel 2,7% (k = 3) e 1,8% (k = 5) delle repliche (t non corretto: 4,7% e 6,1%), con effetto noto nel 99,3%; l'IC espanso del cluster bootstrap copre il 93,5-94,8% con 10-11 pazienti (percentile semplice 88,8-92,3%), il 91,5-92,5% con 5 pazienti e solo l'82-83% con 4.
- (b) `sim_b_reference_schemes.py`: i tre schemi di riferimento di `10_loco.py` su cloni simulati, sotto H0 e con un artefatto atteso nel tumore (perdita del segnale CD8 ed espressione ambientale di CD4). Risultati: sotto H0 i tre schemi stimano l'errore vero entro 0,005; con l'artefatto lo schema che esclude il compartimento giudicato (LOCO) gonfia l'errore del sangue a 0,126 (errore vero 0,093) con artefatto moderato e a 0,385 (vero 0,092; fino a 0,442 in una replica) con artefatto forte, mentre gli schemi solo sangue e naive restano entro 0,008 dall'errore vero. I parametri sono stati scelti a priori: il valore 0,24 citato nella versione precedente non e' riprodotto esattamente.
- (c) `sim_c_block_permutation.py`: falsi positivi del block permutation test sotto H0 con confondente di paziente. Risultati: su 100 repliche il p non centrato e' < 0,05 in 4, 6, 11 e 5 casi (confondente di forza 0; 0,5; 1; 2), su 1.000 repliche (forza 1) nel 5,4%; il p centrato nel 5,3%; la permutazione che ignora i pazienti nell'85,5%. Il dato "2 falsi positivi su 100" della versione precedente non e' riprodotto: con 100 repliche il conteggio oscilla fra 4 e 11.
- (d) `sim_d_xai_stability.py`: la procedura di `20` su conteggi simulati con e senza effetto paziente. Risultati: differenza controllo - LOPO +0,143 con effetto paziente e +0,030 senza (5 repliche per condizione; nella versione precedente erano citati 0,063 e 0,040).

**Test (`tests/`, pytest).**
- `test_scpipe.py`: formule confrontate con calcoli a mano (Nadeau-Bengio, Jaccard, LogCPM, HVGSelector).
- `test_leakage.py`: per 05, 11, 12, 13 e 20, i parametri stimati in un fold non cambiano quando si alterano le cellule di test; gli split di 04 non mettono lo stesso paziente in due fold.
- `test_stats.py` (bande dichiarate nel docstring): copertura del cluster bootstrap in [0,90; 0,99] con 10 e 11 pazienti; falsi positivi del block permutation test in [0,02; 0,08]; Nadeau-Bengio sotto H0 ≤ 0,08 e potenza ≥ 0,5; Holm monotono e ≤ 1; stimatore della discordanza senza reinserimento non distorto (enumerazione esatta).
- `test_smoke.py`: ogni script si importa e risponde a `--help`.

L'intervallo percentile semplice del cluster bootstrap usato nella versione precedente copriva solo l'88-89% con 10-11 pazienti ed e' stato sostituito dall'intervallo "expanded percentile" (Hesterberg 2015), che copre il 92-94%. Con 4-5 pazienti nemmeno l'intervallo espanso raggiunge il 95% (vedi `sim_a_cluster_bootstrap.csv`).

---

## 7. Limiti

- **Pochi pazienti.** Le conclusioni fra pazienti poggiano su 5 (task1a_paired, Adjacent_normal) - 14 pazienti; con 4 pazienti (Adjacent_normal nello schema solo sangue) l'IC del cluster bootstrap non e' calibrato e le ricampionature distinte sono 35.
- **Etichette di riferimento.** Le etichette `all_celltype` sono quelle degli autori, ottenute dalla stessa matrice: per `task2` le prestazioni misurano la riproducibilita' dell'annotazione, non la correttezza biologica. L'etichetta di malignita' e' derivata da inferCNV sugli stessi dati, con riferimento comune a tutti i pazienti e soglia stimata sulle stesse cellule di riferimento.
- **Identita' clonale dai marcatori.** CD4/CD8 da pochi UMI per cellula; i risultati principali usano CD8A > 0 o CD8B > 0, con regole piu' severe in `error_rate_identity_rules.csv`. La soglia ambientale di CD4 nel PBMC poggia su 11 cellule.
- **Doppietti.** d1 e d2 non sono indipendenti e dipendono da mediane pooled fra tessuti.
- **CNV.** Le coordinate delle finestre di inferCNV sono interpolate in modo uniforme per cromosoma: i confini dei loci sono approssimati.
- **Ipotesi per paziente.** n = 10, potenza bassa; H1 ha una circolarita' (punteggio calcolato sulle cellule etichettate CD8T) e una direzione ipotizzata opposta agli errori osservati.
- **Coorte snRNA-seq** esclusa: pazienti diversi dalla coorte scRNA-seq.
- **Una sola coorte**: nessuna validazione esterna.

---

## Autore

Emanuele Colasanto, matricola 859776 — Corso di Machine Learning, A.A. 2025/2026.

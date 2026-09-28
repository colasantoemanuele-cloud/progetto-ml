# Analisi Multi-Modale di Machine Learning e Validazione Clonale TCR nel Carcinoma Duttale Pancreatico (PDAC)

Studio computazionale e metodologico basato sulla **NCBI GEO SuperSeries `GSE278694`** (Chen et al., *Cancer Cell* 2025): **199.184 cellule**, **14 pazienti**, **30.543 geni**, **29 librerie scRNA-seq** (`GSE278688`) e **30 librerie scTCR-seq V(D)J** (`GSE300435`) distribuite su tre compartimenti anatomici (`Tumor`, `Adjacent_normal`, `PBMC`). La coorte a singolo nucleo **snRNA-seq (`GSE278689`)** e' stata esclusa dalle analisi comparative supervisionate in quanto acquisita su un insieme disgiunto di soggetti, evitando il confondimento strutturale tra modalita' tecnica di sequenziamento e variabilita' biologica inter-paziente.

---

## 1. I Tre Problemi di Studio e la Sintesi dei Risultati

Il lavoro esamina criticamente tre task canonici di apprendimento automatico applicato alla trascrittomica a singola cellula del PDAC, dimostrando come le metriche predittive apparentemente elevate ottenute con protocolli standard derivino da tre distinti artefatti del disegno sperimentale, riconducibili a un unico **principio unificante**: **l'unita' statistica indipendente e' il paziente, non la singola cellula**.

| Problema di Studio | Popolazione e Dimensione | Metrica in Pratica Corrente (Split Casuale / Naive) | Stima Corretta per la Struttura di Gruppo (LOPO / Entro Paziente) | Origine Metodologica dell'Artefatto |
| :--- | :--- | :--- | :--- | :--- |
| **1. Annotazione del tipo cellulare** (`task2_v1`, 9 classi) | 199.184 cellule, 14 pazienti, 30.543 geni | **Macro-F1 = 0,954** | **Reale ma tautologica** (Macro-F1 LOPO $\approx 0,938$) | **Circolarita' definitoria**: le etichette di riferimento derivano da clustering non supervisionato e marker sulla medesima matrice di espressione usata per addestrare i classificatori. |
| **2. Classificazione del compartimento tissutale** (`task1a_paired`, Tumore vs Normale) | 67.662 cellule, 5 pazienti con tessuti appaiati (`Tumor` + `Adjacent_normal`) | **Macro-F1 = 0,875 ($\pm 0,002$)** | **Macro-F1 = 0,729 ($\pm 0,100$)** ($\Delta = -0,146$; incertezza reale **$50\times$ superiore**) | **Confondimento paziente $\leftrightarrow$ tessuto e fuga di informazione**: lo split casuale assegna cellule dello stesso individuo a train e test, memorizzando la firma paziente-specifica anziche' la biologia tumorale. |
| **3. Malignita' delle cellule duttali** (`13_infercnv.py` / `14_cnv_loci.py`) | Cellule duttali neoplastiche vs epitelio non maligno (profilo CNV) | **AUC LOPO = 0,825** (classificatore su CNV globale) | **AUC $\approx 0,53$ entro paziente** (test su duttali di `Adjacent_normal`) | **Gradiente epiteliale ed effetto di campo di cancerizzazione**: le cellule duttali del tessuto adiacente istologicamente normale condividono gia' alterazioni molecolari e trascrizionali con il tumore primario dello stesso paziente. |

### Validazione Ortogonale tramite Tracciamento Clonale TCR e Spiegabilita' (XAI)

Per superare la circolarita' delle etichette trascrizionali senza ricorrere alla matrice di espressione genica, lo studio utilizza il **recettore delle cellule T (scTCR-seq V(D)J)** come verita' biologica esterna (poiche' il riarrangiamento somatico delle catene $\alpha/\beta$ precede l'espansione periferica ed e' invariante lungo il lignaggio cellulare):

* **Discordanza compartimento-specifica (Schema Primario a Riferimento Ematico, `10_loco.py`)**: assumendo esclusivamente il sangue periferico (`PBMC`) come riferimento esterno privo di circolarita', le cellule T appartenenti allo stesso clone mostrano un tasso di errore dell'annotazione trascrizionale significativamente superiore nel compartimento tumorale (`Tumor`) rispetto al tessuto adiacente (`Adjacent_normal`).
* **Decomposizione tra plasticita' fenotipica e doppietti tecnici (`08_disambiguate.py`, `09_arbitrate.py`, `10_loco.py`)**: l'arbitrato sui marcatori continui (`CD8A`/`CD8B` vs `CD4`, moduli NK e mieloidi) e il discriminante fisico delle conte totali UMI dimostrano che la discordanza `CD8T <-> CD4T` riflette l'attenuazione trascrizionale da esaurimento funzionale intratumorale (`PDCD1`, `HAVCR2`, `TOX`, `TIGIT`), mentre tra il **72% e il 79% delle cellule annotate come `NK` ma positive per un recettore TCR produttivo (`NK-TCR+`) costituisce un artefatto da doppietti fisici T-NK**.
* **Espansione clonale e sensibilita' Leave-One-Patient-Out (`17_clonal_expansion.py`, `18_loyo_check.py`)**: la dimensione del clone T ($\text{PBMC} + \text{Tumor}$) correla positivamente sia con il tasso di errore binario nel tumore ($\rho = +0,177$, $p = 0,0009$ via *Block Permutation Test* intra-paziente) sia con l'eccesso continuo di discordanza ($\rho = +0,281$, $p < 0,0001$), risultando **robusta in 10/10 iterazioni Leave-One-Patient-Out**.
* **Spiegabilita' per paziente (`20_xai_per_patient.py`)**: la similarita' di Jaccard sui `top-50` geni selezionati dalla Regressione Logistica L1 crolla da **0,459** (controllo senza struttura per paziente) a **0,300** nei fold LOPO, confermando che l'esclusione di un singolo paziente altera sistematicamente la frontiera decisionale appresa dal classificatore.

---

## 2. Architettura della Pipeline e Moduli Sorgente (`src/`)

La cartella [`src/`](src/) contiene i moduli Python della pipeline analitica e di modellazione, progettati per operare in formato sparso (`scipy.sparse.csr_matrix` `float32`) con trasformatori Scikit-learn esenti da data leakage:

| Modulo Sorgente | Funzione Tecnica e Metodologica nel Flusso di Studio |
| :--- | :--- |
| [`src/00_fetch_metadata.py`](src/00_fetch_metadata.py) | Triage e download selettivo dei metadati XML/SOFT e dell'indice dei file supplementari dalla NCBI GEO SuperSeries `GSE278694`. |
| [`src/01_probe_labels.py`](src/01_probe_labels.py) | Verifica preliminare delle etichette cellulari su disco e individuazione dei 5 pazienti dotati di tessuti solidi appaiati (`Tumor` + `Adjacent_normal`). |
| [`src/01b_tcr_check.py`](src/01b_tcr_check.py) | Mappatura biunivoca e allineamento dei barcode nucleotidici 10x Genomics tra i campioni trascrizionali scRNA-seq e i 30 file di contig V(D)J (`GSE300435`). |
| [`src/02_build_h5ad.py`](src/02_build_h5ad.py) | Conversione in streaming a blocchi (*chunk* da 20M di righe) da formato MatrixMarket compresso (`.mtx.gz`) a matrice sparsa `CSR float32` salvata in `AnnData` (`.h5ad`). |
| [`src/03_eda_qc.py`](src/03_eda_qc.py) | Controllo di qualita' (QC) per singola libreria di sequenziamento, analisi della complessita' UMI/geni e verifica del supporto per paziente per ciascuna classe cellulare. |
| [`src/04_splits.py`](src/04_splits.py) | Generazione deterministica delle partizioni di validazione (`StratifiedKFold` casuale per cellula usato come controllo negativo, `StratifiedGroupKFold` a 5 fold e `Leave-One-Patient-Out`). |
| [`src/05_baseline.py`](src/05_baseline.py) | Quantificazione del divario di prestazione e di incertezza tra split casuale e split per paziente mediante classificatori `Dummy` e `LogisticRegression` con penalizzazione L2. |
| [`src/06_clone_feasibility.py`](src/06_clonotype_consistency.py) / [`src/06_clonotype_consistency.py`](src/06_clonotype_consistency.py) | Valutazione della fattibilita' dell'ancoraggio clonale sangue-tumore, stima preliminare della contaminazione NK e misura della discordanza intra- e inter-compartimentale nei cloni TCR condivisi. |
| [`src/07_robustness.py`](src/07_robustness.py) | Stimatore a coppie $1\times 1$ non distorto per la cardinalita' campionaria, *Cluster Bootstrap* sui pazienti, discriminante basato su UMI TCR ad alta confidenza e verifica dei marcatori `CD4`/`CD8`. |
| [`src/08_disambiguate.py`](src/08_disambiguate.py) | Contrasto appaiato sui cloni tripli presenti simultaneamente in tutti e tre i compartimenti (`PBMC`, `Adjacent_normal`, `Tumor`), profilo dei marcatori su cloni concordanti e analisi doppietti con Scrublet. |
| [`src/09_arbitrate.py`](src/09_arbitrate.py) | Arbitrato dell'identita' clonale basato sui moduli di espressione dei marcatori canonici e discriminante fisico dei doppietti (conteggi UMI totali e rapporto UMI TCR). |
| [`src/10_loco.py`](src/10_loco.py) | Implementazione dello schema *Leave-One-Compartment-Out* (LOCO) e dello schema primario *Blood-only* (riferimento fissato al solo compartimento `PBMC`), con stima della frazione di doppietti nelle cellule `NK-TCR+` (72-79%). |
| [`src/11_model_zoo.py`](src/11_model_zoo.py) | Benchmark sistematico di **13 classificatori** in validazione LOPO (modelli lineari L1/L2/ElasticNet/LinearSVC, alberi, Random Forest, ExtraTrees, LightGBM e modelli densi KNN, SVM-RBF, MLP, GaussianNB, QDA su proiezione `TruncatedSVD` a 50 componenti) con test $t$ corretto di Nadeau-Bengio e test di Wilcoxon. |
| [`src/12_feature_selection.py`](src/12_feature_selection.py) | Benchmark di **6 criteri di selezione delle feature** (`Random`, `HVG`, `ANOVA`, `chi2`, `L1 logreg`, `Random Forest`) valutati su $k \in \{50, 200, 1000\}$, con misura della stabilita' strutturale inter-fold tramite indice di Jaccard. |
| [`src/13_infercnv.py`](src/13_infercnv.py) | Inferenza del profilo di Copy Number Variation (CNV) a finestra mobile cromosomica assumendo come riferimento interno le **cellule immunitarie dello stesso paziente** (`CD8T`, `CD4T`, `NK`, `B`, `Myeloid`, `Plasma`), addestramento LOPO e test esterno sulle cellule duttali adiacenti normali. |
| [`src/14_cnv_loci.py`](src/14_cnv_loci.py) | Decomposizione dello score CNV sui loci citogenetici ricorrenti del PDAC (`9p21`, `17p`, `18q`, `1p`, `6q`, `8q`, `12p`, `20q`) e dimostrazione quantitativa dell'effetto di campo di cancerizzazione. |
| [`src/15_inventory.py`](src/15_inventory.py) | Inventario programmatico completo delle variabili disponibili nei metadati e nei file V(D)J per documentare a priori lo spazio delle variabili ed evitare selezioni esplorative post-hoc. |
| [`src/16_pseudobulk_patient.py`](src/16_pseudobulk_patient.py) | Aggregazione pseudobulk per paziente ($n = 10\text{-}14$) e test di permutazione esatto sulle tre ipotesi biologiche a priori (`H1`: esaurimento CD8T, `H2`: frazione mieloide, `H3`: frazione CAF) con correzione sequenziale di Holm-Bonferroni. |
| [`src/17_clonal_expansion.py`](src/17_clonal_expansion.py) | Test dell'associazione positiva tra dimensione di espansione clonale ed errore di annotazione nel tumore mediante *Block Permutation Test* intra-paziente e *Cluster Bootstrap*. |
| [`src/18_loyo_check.py`](src/18_loyo_check.py) | Analisi di sensibilita' *Leave-One-Patient-Out* (LOYO) sull'associazione tra espansione clonale ed errore di annotazione (significativita' preservata in 10/10 iterazioni). |
| [`src/19_ni_score.py`](src/19_ni_score.py) | Test dell'ipotesi clinica `H4` sul punteggio quantitativo di invasione neurale (`NI score`, da `Table S1A` / `mmc2.xlsx`) e ricalcolo congiunto di Holm-Bonferroni sulla famiglia completa `H1-H4`. |
| [`src/20_xai_per_patient.py`](src/20_xai_per_patient.py) | Analisi XAI per paziente: quantificazione della similarita' delle spiegazioni (Jaccard sui `top-50` geni della Regressione Logistica L1) nei fold LOPO (`0,300`) rispetto al controllo casuale privo di struttura per paziente (`0,459`). |
| [`src/21_fig_blood_scheme.py`](src/21_fig_blood_scheme.py) | Generazione della figura riassuntiva e della tabella di errore per lo schema primario di arbitrato a riferimento ematico (`PBMC`). |
| [`src/scpipe.py`](src/scpipe.py) | Modulo infrastrutturale condiviso contenente i trasformatori Scikit-learn (`LogCPM`, `HVGSelector`, `IndexSelector`) e le funzioni statistiche per il confronto di modelli (`nadeau_bengio_ttest`, `jaccard`, `mean_pairwise_jaccard`). |
| [`src/estrai_conteggi.py`](src/estrai_conteggi.py) / [`src/download_tcr.sh`](src/download_tcr.sh) | Utilità di supporto per il riepilogo dei conteggi cellulari nella coorte appaiata e script Bash idempotente con controllo di integrita' `gzip` per il download dei 30 file V(D)J. |

---

## 3. Ambiente di Esecuzione e Riproducibilita'

L'ambiente computazionale e' definito in modo dichiarativo in [`envs/environment.yml`](envs/environment.yml) e vincolato alle versioni esatte dei pacchetti binari in [`envs/conda-lock.yml`](envs/conda-lock.yml) (**Python 3.11**, `scanpy`, `anndata`, `scikit-learn`, `lightgbm`, `scipy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pybiomart`, `openpyxl`).

```bash
# Creazione dell'ambiente riproducibile tramite Mamba o Conda
mamba env create -f envs/environment.yml
mamba activate pdac-ml
```

Tutti gli algoritmi stocastici (split di cross-validation, sottocampionamento stratificato, inizializzazione dei classificatori, *Block Permutation Test* da 20.000 permutazioni e *Cluster Bootstrap* da 2.000-5.000 repliche) utilizzano un seme pseudocasuale fissato globalmente (`SEED = 0`).

---

## 4. Guida all'Esecuzione (Sintassi CLI Reale Verificata)

Tutti i comandi seguenti devono essere eseguiti dalla directory radice del repository (`pdac-ml/`) e utilizzano gli argomenti CLI effettivamente definiti tramite `argparse` nei rispettivi moduli sorgente:

```bash
# 1. Triage dei metadati GEO, verifica delle etichette e allineamento barcode scRNA-TCR
python src/00_fetch_metadata.py --outdir data/raw/metadata
bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr
python src/01_probe_labels.py --outdir data/raw/probe
python src/01b_tcr_check.py --probe data/raw/probe

# 2. Conversione MatrixMarket -> AnnData CSR float32, QC per libreria e generazione degli split
python src/02_build_h5ad.py --which sc --outdir data --no-counts-layer
python src/03_eda_qc.py --h5ad data/interim/sc_raw.h5ad --tag sc
python src/04_splits.py --h5ad data/interim/sc_raw.h5ad --tag sc

# 3. Quantificazione del divario di generalizzazione, benchmark dei 13 modelli, feature selection e XAI
python src/05_baseline.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv
python src/11_model_zoo.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv --tasks task2_v1 task2_hier3
python src/12_feature_selection.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv --task task2_v1
python src/20_xai_per_patient.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv --task task1a_paired

# 4. Tracciamento clonale TCR, controlli di robustezza, arbitrato multi-marker e schema primario PBMC
python src/06_clonotype_consistency.py --h5ad data/interim/sc_raw.h5ad --filelist data/raw/metadata/filelist.txt --outdir data/raw/tcr
python src/07_robustness.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
python src/08_disambiguate.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
python src/09_arbitrate.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
python src/10_loco.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
python src/21_fig_blood_scheme.py

# 5. Inferenza CNV (riferimento immunitario intra-paziente), loci PDAC ed effetto di campo
python src/13_infercnv.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv
python src/14_cnv_loci.py --h5ad data/interim/sc_raw.h5ad

# 6. Analisi pseudobulk per paziente (ipotesi H1-H3), espansione clonale, LOYO e NI score (H4)
python src/15_inventory.py
python src/16_pseudobulk_patient.py --h5ad data/interim/sc_raw.h5ad --arbitration results/tables/arbitration_blood.csv
python src/17_clonal_expansion.py --arbitration results/tables/arbitration_blood.csv --pairwise results/tables/pairwise_discordance_clone_b.csv
python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
python src/19_ni_score.py --table-s1 mmc2.xlsx --arbitration results/tables/arbitration_blood.csv --previous-tests results/tables/pseudobulk_hypothesis_tests.csv
```

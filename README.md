# Analisi Multi-Modale di Machine Learning e Validazione Clonale TCR nel Carcinoma Duttale Pancreatico (PDAC)

Studio computazionale e metodologico basato sulla **NCBI GEO SuperSeries `GSE278694`** (Chen et al., *Cancer Cell* 2025): **199.184 cellule**, **14 pazienti**, **30.543 geni**, **29 librerie scRNA-seq** (`GSE278688`) e **30 librerie scTCR-seq V(D)J** (`GSE300435`) distribuite su tre compartimenti anatomici (`Tumor`, `Adjacent_normal`, `PBMC`). La coorte a singolo nucleo **snRNA-seq (`GSE278689`, 87.693 nuclei)** e' stata esclusa dalle analisi comparative supervisionate in quanto acquisita su un insieme disgiunto di soggetti, evitando il confondimento strutturale tra modalita' tecnica di sequenziamento e variabilita' biologica inter-paziente.

---

## 1. I Tre Problemi di Studio e la Sintesi dei Risultati

Il lavoro esamina criticamente tre problemi canonici di apprendimento automatico supervisionato applicato alla trascrittomica a singola cellula del PDAC, dimostrando come le stime di accuratezza apparentemente elevate ottenute secondo la pratica corrente riflettano tre distinti artefatti del disegno di validazione, riconducibili a un unico **principio unificante**: **l'unita' statistica indipendente e' il paziente, non la singola cellula**.

| Problema di Studio | Popolazione e Dimensione | Metrica in Pratica Corrente (Split Casuale / Naive) | Stima Corretta per la Struttura di Gruppo (LOPO / Entro Paziente) | Origine Metodologica dell'Artefatto |
| :--- | :--- | :--- | :--- | :--- |
| **1. Annotazione del tipo cellulare** (`task2_v1`, 9 classi) | 199.184 cellule, 14 pazienti, 30.543 geni | **Macro-F1 = 0,954** (`LightGBM = 0,962`) | **Reale ma tautologica** | **Circolarita' definitoria**: le etichette di riferimento derivano da clustering non supervisionato e marcatori sulla medesima matrice trascrizionale usata per addestrare i classificatori. |
| **2. Classificazione del compartimento tissutale** (`task1a_paired`, Tumore vs Normale) | 67.662 cellule, 5 pazienti con tessuti appaiati (`Tumor` + `Adjacent_normal`) | **Macro-F1 = 0,875 ($\pm 0,002$)** | **Macro-F1 = 0,729 ($\pm 0,100$)** ($\Delta = -0,146$; incertezza reale **$50\times$ superiore**) | **Confondimento paziente $\leftrightarrow$ tessuto e fuga di informazione**: lo split casuale assegna cellule dello stesso individuo a training e test, memorizzando la firma paziente-specifica e la composizione cellulare anziche' la trasformazione neoplastica. |
| **3. Malignita' delle cellule duttali** (`13_infercnv.py` / `14_cnv_loci.py`) | Cellule duttali neoplastiche vs epitelio non maligno (profilo CNV) | **AUC LOPO = 0,825 ($\pm 0,079$)** (31,7% di falsi positivi su `Adjacent_normal`) | **AUC $\approx 0,53$ entro paziente** (range `0,51-0,58` nei 4 pazienti valutabili) | **Gradiente epiteliale ed effetto di campo di cancerizzazione**: le cellule duttali del tessuto adiacente istologicamente normale presentano alterazioni cromosomiche a larga scala sovrapponibili a quelle del tumore primario. |

### Benchmark dei Modelli, Selezione delle Feature e Spiegabilita' per Paziente (XAI)

* **Geometria del problema e confronto tra 13 algoritmi (`11_model_zoo.py`)**: nella classificazione a 9 tipi cellulari (`task2_v1`), la separazione e' pressoche' lineare (`Logistic L2 = 0,943` vs `LightGBM = 0,962`), mentre restringendo l'analisi al sottoproblema con confine trascrizionale sfumato (`CD4T` / `CD8T` / `NK`, `task2_hier3`) gli ensemble ad albero e i kernel non lineari guadagnano oltre **7 punti di Macro-F1** (`LightGBM = 0,923`, `SVM-RBF = 0,906`, `MLP = 0,903` vs `Logistic L2 = 0,851`). Sul task `Tumor` vs `Adjacent_normal` in validazione LOPO, tutti i modelli convergono nel range `0,710-0,725` senza differenze statisticamente significative al test $t$ appaiato corretto di Nadeau-Bengio ($p \ge 0,080$).
* **Prestazione predittiva contro stabilita' strutturale della selezione (`12_feature_selection.py`)**: a bassa cardinalita' ($k = 50$), i selettori supervisionati (`Logistic L1 = 0,936` su `task2_v1` e `0,729` su `task1a_paired`; `Random Forest = 0,925` e `0,732`) superano nettamente i criteri non supervisionati (`HVG = 0,655` e `0,648`). Tuttavia, valutando la stabilita' della firma genica attraverso i fold LOPO mediante l'indice di Jaccard medio a coppie, la penalizzazione `L1` risulta strutturalmente instabile (`Jaccard = 0,439` a $k = 200$ su `task2_v1` e `0,275` a $k = 1000$ su `task1a_paired`) rispetto ai filtri univariati (`Chi-quadro = 0,845` e `0,674`; `ANOVA = 0,817` e `0,657`; `HVG = 0,797` e `0,767`).
* **Divergenza delle spiegazioni XAI tra pazienti (`20_xai_per_patient.py`)**: sui 5 fold LOPO del task `task1a_paired`, l'indice di Jaccard medio sui `top-50` geni selezionati dalla Regressione Logistica L1 e' pari a **0,300** (range `0,220-0,408`; correlazione di Spearman sui vettori completi dei coefficienti = **0,234**; soli **9 geni su 50** condivisi da tutti i 5 fold: `BTG2`, `CD69`, `CREM`, `CTRB1`, `DUSP4`, `FOS`, `GAPDH`, `HSPH1`, `VIM`), mentre sale a **0,459** ($\Delta = +0,159$) nel controllo ottenuto ricampionando il medesimo numero di cellule senza struttura per paziente. La rimozione di un singolo paziente altera dunque la regola appresa dal classificatore molto piu' della rimozione casuale dello stesso numero di cellule, dimostrando per via interpretativa la ridondanza informativa intra-soggetto.

### Validazione Ortogonale tramite Tracciamento Clonale TCR e Test Meccanicistici

Poiche' il riarrangiamento somatico V(D)J delle catene $\alpha/\beta$ precede l'espansione periferica ed e' invariante rispetto allo stato trascrizionale, il recettore delle cellule T (`scTCR-seq`) fornisce un riferimento biologico indipendente dal trascrittoma per quantificare l'affidabilita' dell'annotazione:

* **Rumore di fondo e dipendenza dal compartimento (`06_clonotype_consistency.py`, `07_robustness.py`)**: la discordanza di etichetta tra cellule appartenenti al medesimo clone, nello stesso paziente e nello stesso tessuto, attesta un **rumore di annotazione di base pari a 0,107 (~11%)**. Tra compartimenti anatomici differenti si osserva un **eccesso sistematico di discordanza pari a +0,059** (IC 95% Cluster Bootstrap: `+0,034; +0,100`, $n = 11$ pazienti; robusto sia sulla sola catena `TRB`: `+0,059`, sia sulla coppia `TRA+TRB`: `+0,065`).
* **Tasso di errore nello schema primario a riferimento ematico (`10_loco.py`, `21_fig_blood_scheme.py`)**: fissando esclusivamente il sangue periferico (`PBMC`) come riferimento esterno privo di circolarita', il **tasso di errore dell'annotazione nel tessuto tumorale (`Tumor`) raggiunge 0,195** (IC 95%: `0,123; 0,274`, $n = 10$ pazienti), oltre tre volte superiore a quello osservato nel **tessuto adiacente normale (`Adjacent_normal`), pari a 0,059** (IC 95%: `0,000; 0,084`, $n = 4$ pazienti), con una differenza appaiata sugli stessi cloni pari a **+0,080** (IC 95%: `+0,021; +0,170`, 88 cloni). La verifica sui marcatori canonici (`08_disambiguate.py`, `09_arbitrate.py`) conferma che le cellule dei cloni discordanti etichettate come `CD8T` nel tumore esprimono `CD8A` solo nel **7%** dei casi (contro l'**88%** dei cloni concordanti) ed esprimono `CD4` a un livello **16 volte superiore**.
* **Distinzione fisica tra eterodoppietti e cellule T reali nelle popolazioni `NK-TCR+` (`09_arbitrate.py`, `10_loco.py`)**: tra le cellule annotate come non-T ma portatrici di un contig TCR produttivo (`NK 5,9%`, `Endothelial 5,6%`, `Stellate 5,5%`, `Platelet 13,8%`), due indicatori fisici indipendenti dal clustering trascrizionale (il numero totale di UMI per cellula `n_counts`, che nei doppietti approssima la somma delle due cellule incapsulate, e il rapporto di diluizione degli UMI TCR sul totale) convergono nell'attribuire a **doppietti fisici T-NK tra il 72% e il 79% delle cellule `NK-TCR+`** (`n_counts: 0,723`, IC `0,579-0,840`; rapporto UMI: `0,790`, IC `0,682-0,863`).
* **Espansione clonale come predittore indipendente e robustezza LOYO (`17_clonal_expansion.py`, `18_loyo_check.py`)**: su 287 cloni condivisi tra sangue periferico e tumore ($n = 10$ pazienti), la dimensione clonale correla positivamente sia con l'errore di annotazione binario ($\rho = +0,177$, $p = 0,0009$, *Block Permutation Test* intra-paziente) sia con l'eccesso continuo di discordanza ($\rho = +0,281$, $p < 0,0001$), presentando il medesimo segno anche al livello aggregato di paziente ($\rho = +0,467$, $p = 0,18$). L'analisi di sensibilita' *Leave-One-Patient-Out* (`18_loyo_check.py`) dimostra che l'associazione resta statisticamente significativa in **10/10 iterazioni** per entrambi gli esiti ($\rho \in [0,146; 0,197]$ per l'errore binario; $\rho \in [0,236; 0,368]$ per l'eccesso continuo).
* **Test ed esclusione delle quattro ipotesi biologiche a priori (`16_pseudobulk_patient.py`, `19_ni_score.py`)**: per spiegare l'eterogeneita' del tasso di errore tra pazienti ($n = 10$) senza incorrere in selezioni post-hoc su migliaia di geni, sono state testate mediante permutazione esatta quattro covariate meccanicistiche dichiarate a priori e corrette congiuntamente tramite la procedura step-down di **Holm-Bonferroni**:
  1. `H1` (Punteggio di esaurimento `CD8T` tumorali): $\rho = -0,47$, $p_{\text{Holm}} = 0,53$ (esclusa, direzione opposta);
  2. `H2` (Frazione di cellule mieloidi nel tumore): $\rho = -0,54$, $p_{\text{Holm}} = 0,45$ (esclusa, direzione opposta);
  3. `H3` (Frazione stromale di fibroblasti `CAF`): $\rho = -0,15$, $p_{\text{Holm}} = 1,00$ (esclusa, direzione opposta);
  4. `H4` (Punteggio clinico di invasione neurale `NI score`, `Table S1A`): $\rho = -0,14$, $p_{\text{Holm}} = 1,00$ (esclusa, direzione opposta).

---

## 2. Architettura della Pipeline e Moduli Sorgente (`src/`)

La cartella [`src/`](src/) contiene i 25 moduli Python della pipeline analitica e di modellazione (oltre allo script [`src/download_tcr.sh`](src/download_tcr.sh)), progettati per operare in formato sparso (`scipy.sparse.csr_matrix` `float32`) con trasformatori Scikit-learn esenti da data leakage:

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
| [`src/11_model_zoo.py`](src/11_model_zoo.py) | Benchmark sistematico di **13 classificatori** in validazione LOPO (Dummy, Naive Bayes, KNN, Logistic L2, Logistic L1, ElasticNet, Albero decisionale, Bagging, Random Forest, LightGBM, SVM lineare, SVM-RBF, MLP; con TruncatedSVD a 50 componenti per i modelli densi basati su distanze euclidee) con test t appaiato corretto di Nadeau-Bengio e test di Wilcoxon. |
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

Tutti gli algoritmi stocastici (partizioni di cross-validation, sottocampionamento stratificato a budget fisso di 20.000 cellule per fold, inizializzazione dei classificatori, *Block Permutation Test* da 20.000 permutazioni e *Cluster Bootstrap* sui pazienti da 2.000-5.000 repliche) utilizzano un seme pseudocasuale fissato globalmente (`SEED = 0`).

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

# 6. Inventario delle variabili, pseudobulk per paziente (H1-H3), espansione clonale, LOYO e NI score (H4)
python src/15_inventory.py --h5ad data/interim/sc_raw.h5ad --tcrdir data/raw/tcr
python src/16_pseudobulk_patient.py --h5ad data/interim/sc_raw.h5ad --arbitration results/tables/arbitration_blood.csv
python src/17_clonal_expansion.py --arbitration results/tables/arbitration_blood.csv --pairwise results/tables/pairwise_discordance_clone_b.csv
python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
python src/19_ni_score.py --table-s1 mmc2.xlsx --arbitration results/tables/arbitration_blood.csv --previous-tests results/tables/pseudobulk_hypothesis_tests.csv
```

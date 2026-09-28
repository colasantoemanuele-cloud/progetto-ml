# Analisi Multi-Modale di Apprendimento Automatico e Validazione Clonale TCR nel Carcinoma Duttale Pancreatico (PDAC)

## 1. Inquadramento Scientifico e Obiettivi dello Studio

Questo repository implementa una pipeline riproducibile di Machine Learning e bioinformatica statistica applicata all'atlante trascrizionale e recettoriale a singola cellula del carcinoma duttale pancreatico umano (**PDAC**, *Pancreatic Ductal Adenocarcinoma*, NCBI GEO SuperSeries **GSE278694**; Chen et al., *Cancer Cell* 2025). Il dataset integra tre modalita' sperimentali complementari acquisite su una coorte di 14 pazienti:

* **scRNA-seq (`GSE278688`)**: 199.184 profili trascrizionali a singola cellula (33.538 geni) campionati su tre compartimenti anatomici (`Tumor`, `Adjacent_normal`, `PBMC`);
* **snRNA-seq (`GSE278689`)**: 64.938 profili a singolo nucleo acquisiti su tessuto pancreatico congelato;
* **scTCR-seq (`GSE300435`)**: 30 librerie V(D)J 10x Genomics accoppiate alle medesime sospensioni cellulari scRNA-seq, con ricostruzione produttiva delle catene $\alpha$ e $\beta$ del recettore delle cellule T (CDR3 aminoacidico e nucleotidico).

Lo studio affronta due questioni metodologiche e biologiche strettamente interconnesse:

1. **Quantificazione della pseudoreplicazione e capacita' di generalizzazione inter-paziente**: nella classificazione supervisionata di stati cellulari e compartimenti tissutali da dati scRNA-seq, il campionamento di decine di migliaia di cellule da un numero limitato di individui induce una forte correlazione intra-soggetto. La pipeline misura sistematicamente il divario di prestazione tra schemi di validazione ingenui a livello di cellula (`StratifiedKFold`) e schemi gerarchici che vincolano l'indipendenza del soggetto di test (`StratifiedGroupKFold` e `Leave-One-Patient-Out`, LOPO), valutando famiglie di modelli lineari e non lineari, metodi di selezione delle feature e stabilita' delle spiegazioni XAI.
2. **Validazione ortogonale delle etichette trascrizionali tramite tracciamento clonale TCR**: poiche' il riarrangiamento somatico V(D)J delle catene $\alpha/\beta$ precede l'espansione periferica ed e' invariante rispetto allo stato di attivazione, tutte le cellule appartenenti al medesimo clone T condividono per discendenza biologica la stessa identita' di lignaggio (`CD8T` vs `CD4T`). Utilizzando il recettore TCR come verita' di riferimento indipendente dal trascrittoma, la pipeline misura il tasso di discordanza delle annotazioni basate su clustering trascrizionale durante l'infiltrazione dal sangue periferico (`PBMC`) al microambiente tumorale (`Tumor`), discriminando tra plasticita' trascrizionale indotta da esaurimento funzionale e contaminazione tecnica da doppietti eterotipici o RNA ambientale.

---

## 2. Impostazione Metodologica e Principi di Ingegneria Statistica

L'architettura del codice rispetta quattro principi rigorosi di progettazione sperimentale e computazionale:

* **Isolamento dei cluster paziente e prevenzione del data leakage**: ogni trasformazione che stima parametri sulla distribuzione dei dati (selezione dei geni altamente variabili tramite indice di dispersione Fano in `HVGSelector`, pre-filtraggio di espressione minima, standardizzazione sparsa `StandardScaler(with_mean=False)`, riduzione dimensionale `TruncatedSVD` e selettori supervisionati `F-ANOVA`, `Mutual Information`, `L1-Lasso`) e' incapsulata in pipeline Scikit-learn (`src/scpipe.py`) e addestrata esclusivamente sul fold di training.
* **Inferenza statistica corretta per la dipendenza tra osservazioni**:
  * Nel confronto tra algoritmi di classificazione e metodi di feature selection sui fold di validazione incrociata, la significativita' delle differenze di Macro-F1 e Balanced Accuracy e' valutata mediante il **test $t$ appaiato corretto di Nadeau & Bengio (2003)**, che incorpora il fattore $(1/k + n_{\text{test}}/n_{\text{train}})$ per compensare la covarianza indotta dalla sovrapposizione dei training set, associato alla procedura step-down di **Holm-Bonferroni** per il controllo forte del *Family-Wise Error Rate* (FWER).
  * Nell'analisi della discordanza clonale TCR, gli intervalli di confidenza al 95% sono stimati tramite **Cluster Bootstrap** (ricampionamento con reinserimento di interi pazienti, preservando la correlazione intra-paziente), mentre l'associazione tra espansione clonale ed errore di annotazione nel tumore e' valutata tramite **Block Permutation Test** (20.000 permutazioni ristrette ai confini di ciascun paziente).
* **Eliminazione della circolarita' nell'arbitrato compartimentale (`10_loco.py`)**: la definizione dell'identita' di riferimento di un clone T non include mai il compartimento sottoposto a valutazione. Lo schema primario utilizza esclusivamente il sangue periferico (`PBMC`) come riferimento esterno e indipendente per arbitrare le cellule osservate in `Tumor` e `Adjacent_normal`.
* **Efficienza di memoria su matrici sparse CSR**: l'intero flusso di pre-elaborazione e modellazione opera su matrici sparse `scipy.sparse.csr_matrix` a precisione singola (`float32`), evitando la densificazione dell'atlante da $199.184 \times 33.538$ e consentendo l'esecuzione completa su workstation con 16 GB di RAM.

---

## 3. Struttura del Repository e Architettura dei Moduli (`src/`)

```text
pdac-ml/
├── envs/
│   ├── environment.yml          # Specifica dichiarativa dell'ambiente Conda/Mamba
│   └── conda-lock.yml           # Lockfile deterministico multipiattaforma delle dipendenze
├── src/
│   ├── scpipe.py                # Trasformatori Scikit-learn anti-leakage e test di Nadeau-Bengio
│   ├── download_tcr.sh          # Script Bash idempotente per il download dei 30 file VDJ (GSE300435)
│   ├── 00_fetch_metadata.py     # Acquisizione e verifica dei metadati GEO (GSE278694)
│   ├── 01_probe_labels.py       # Audit strutturale delle etichette e della fattibilita' dei task
│   ├── 01b_tcr_check.py         # Allineamento dei barcode cellulari tra scRNA-seq e scTCR-seq
│   ├── 02_build_h5ad.py         # Costruzione streaming da .mtx.gz a .h5ad sparso (CSR)
│   ├── 03_eda_qc.py             # Analisi esplorativa, controllo di qualita' e verifica di legge di potenza
│   ├── 04_splits.py             # Generazione deterministica degli split (Stratified, Group, LOPO)
│   ├── 05_baseline.py           # Quantificazione del divario di generalizzazione (split casuale vs paziente)
│   ├── 06_clonotype_consistency.py # Consistenza delle annotazioni cellulari nei cloni TCR espansi
│   ├── 07_robustness.py         # Analisi di sensibilita' su definizioni clonali, doppietti (Scrublet) e soglie
│   ├── 08_disambiguate.py       # Analisi dei marker canonici e di esaurimento nelle cellule discordanti
│   ├── 09_arbitrate.py          # Arbitrato multi-marker (CD8 vs CD4, T vs NK, T vs Mieloide)
│   ├── 10_loco.py               # Arbitrato Leave-One-Compartment-Out e schema primario a riferimento PBMC
│   ├── 11_model_zoo.py          # Confronto sistematico di 5 famiglie di classificatori in LOPO
│   ├── 12_feature_selection.py  # Confronto di 5 selettori di feature (prestazione e stabilita' di Jaccard)
│   ├── 13_infercnv.py           # Inferenza del profilo di Copy Number Variation (CNV) a finestra mobile
│   ├── 14_cnv_loci.py           # Analisi di specificita' per locus cromosomico (chr7p, chr18q/SMAD4, chr6p/HLA)
│   ├── 15_inventory.py          # Inventario programmatico completo delle colonne e dei metadati grezzi
│   ├── 16_pseudobulk_patient.py # Aggregazione pseudobulk per paziente e test delle ipotesi a priori H1-H3
│   ├── 17_clonal_expansion.py   # Associazione tra espansione clonale ed errore (Block Permutation e Cluster Boot)
│   ├── 18_loyo_check.py         # Analisi di sensibilita' Leave-One-Patient-Out sull'espansione clonale
│   ├── 19_ni_score.py           # Associazione con il punteggio clinico di invasione neurale (NI score, H4)
│   ├── 20_xai_per_patient.py    # Stabilita' inter-paziente delle spiegazioni XAI (Regressione Logistica L1)
│   ├── 21_fig_blood_scheme.py   # Generazione della figura di sintesi per lo schema primario (riferimento PBMC)
│   └── estrai_conteggi.py       # Utilità di estrazione dei conteggi cellulari per la coorte appaiata
└── results/
    ├── figures/                 # Grafici ad alta risoluzione generati dalla pipeline (.png)
    └── tables/                  # Tabelle numeriche, test di ipotesi e metriche di validazione (.csv)
```

### Sintesi Funzionale dei Moduli

| Modulo | Funzione Principale e Output Scientifico |
| :--- | :--- |
| [`src/00_fetch_metadata.py`](src/00_fetch_metadata.py) | Scarica e indicizza i metadati delle serie GEO (`GSE278694`, `GSE278688`, `GSE278689`, `GSE300435`), producendo l'inventario dei campioni in `data/raw/metadata/`. |
| [`src/01_probe_labels.py`](src/01_probe_labels.py) | Verifica la distribuzione paziente $\times$ tessuto $\times$ tipo cellulare senza caricare la matrice di conteggio, identificando i 5 pazienti dotati di campioni appaiati `Tumor` + `Adjacent_normal`. |
| [`src/01b_tcr_check.py`](src/01b_tcr_check.py) | Ricostruisce la corrispondenza esatta tra i prefissi di campione scRNA-seq (`PAxx_Tissue`) e i 30 file `VDJ_filtered_contig_annotations.csv.gz` tramite intersezione dei barcode 10x. |
| [`src/02_build_h5ad.py`](src/02_build_h5ad.py) | Converte in streaming i file `.mtx.gz` in formato `AnnData` compresso (`sc_raw.h5ad`), preservando i conteggi interi grezzi in formato `csr_matrix` (`int32`) e validando le somme marginali. |
| [`src/03_eda_qc.py`](src/03_eda_qc.py) | Esegue l'analisi esplorativa della composizione cellulare, dei conteggi UMI e della frazione mitocondriale, verificando analiticamente (criterio di Clauset-Shalizi-Newman) il regime delle taglie clonali. |
| [`src/04_splits.py`](src/04_splits.py) | Costruisce la tabella deterministica delle partizioni (`splits_sc.csv`), affiancando per le medesime cellule `StratifiedKFold`, `StratifiedGroupKFold` e `LeaveOneGroupOut` (LOPO). |
| [`src/05_baseline.py`](src/05_baseline.py) | Misura la sovrastima sistematica di Macro-F1 e Balanced Accuracy indotta dallo split casuale per cellula rispetto allo split per paziente su due algoritmi di riferimento (Regressione Logistica L2 e HistGradientBoosting). |
| [`src/06_clonotype_consistency.py`](src/06_clonotype_consistency.py) | Quantifica la discordanza delle etichette trascrizionali all'interno dei cloni T espansi ($\ge 2$ cellule con identica coppia CDR3 $\alpha+\beta$ nello stesso paziente), confrontando cloni intra-compartimento e trans-compartimento. |
| [`src/07_robustness.py`](src/07_robustness.py) | Sottopone il risultato di discordanza clonale a quattro controlli di robustezza: tre definizioni di clonotipo (aa, nt, $\beta$-only), stima dei doppietti con Scrublet per campione, curva al variare della soglia di dimensione clonale e confronto a coppie a cardinalita' bilanciata $1\times 1$. |
| [`src/08_disambiguate.py`](src/08_disambiguate.py) | Analizza il profilo continuo di espressione dei marker di lignaggio (`CD8A`, `CD8B`, `CD4`, `FOXP3`) e di esaurimento (`PDCD1`, `HAVCR2`, `TOX`, `TIGIT`, `CTLA4`) nelle cellule clonalmente discordanti. |
| [`src/09_arbitrate.py`](src/09_arbitrate.py) | Risolve ciascun evento di discordanza confrontando i moduli continui di espressione delle famiglie cellulari coinvolte (`CD8T` vs `CD4T`, `T` vs `NK`, `T` vs `Myeloid`). |
| [`src/10_loco.py`](src/10_loco.py) | Formalizza gli schemi di arbitrato esenti da circolarita': schema **Leave-One-Compartment-Out (LOCO)** e schema primario **BLOOD** (riferimento fissato esclusivamente su `PBMC`), con inferenza via Cluster Bootstrap. |
| [`src/11_model_zoo.py`](src/11_model_zoo.py) | Confronta 5 architetture supervisionate (`Logistic L2`, `Logistic L1`, `LinearSVC`, `RandomForest`, `HistGradientBoosting`) in validazione LOPO tramite pipeline anti-leakage e test di Nadeau-Bengio. |
| [`src/12_feature_selection.py`](src/12_feature_selection.py) | Valuta 5 strategie di selezione delle feature (`HVG`, `F-ANOVA`, `Mutual Information`, `L1-Lasso`, `Random`) al variare di $k \in \{50, 200, 500, 2000\}$, misurando sia il Macro-F1 di test sia la stabilita' di Jaccard tra fold. |
| [`src/13_infercnv.py`](src/13_infercnv.py) | Stima il numero di copie genomiche (CNV) mediante media mobile di 100 geni ordinati lungo i cromosomi autosomici rispetto alle cellule endoteliali e fibroblastiche di riferimento. |
| [`src/14_cnv_loci.py`](src/14_cnv_loci.py) | Decompone il segnale CNV sui loci cromosomici ricorrenti del PDAC (`chr7p` EGFR, `chr18q` SMAD4) e sul complesso MHC (`chr6p`), distinguendo aneuploidia somatica da programmi trascrizionali coordinati. |
| [`src/15_inventory.py`](src/15_inventory.py) | Esegue la ricognizione sistematica di tutte le variabili presenti nei metadati scRNA-seq, snRNA-seq e nelle tabelle V(D)J per mappare il perimetro informativo del dataset. |
| [`src/16_pseudobulk_patient.py`](src/16_pseudobulk_patient.py) | Aggrega i profili a livello di paziente ($n=10\text{-}14$) e sottopone a test di permutazione esatto tre ipotesi meccanicistiche pre-specificate (H1: esaurimento CD8T; H2: frazione mieloide; H3: frazione CAF) con correzione di Holm-Bonferroni. |
| [`src/17_clonal_expansion.py`](src/17_clonal_expansion.py) | Testa l'associazione tra dimensione clonale T ($\text{PBMC} + \text{Tumor}$) e tasso di errore di annotazione nel tumore mediante Block Permutation Test intra-paziente e Cluster Bootstrap. |
| [`src/18_loyo_check.py`](src/18_loyo_check.py) | Verifica la robustezza dell'associazione tra espansione clonale ed errore rimuovendo iterativamente un paziente alla volta (*Leave-One-Patient-Out* sui 10 pazienti). |
| [`src/19_ni_score.py`](src/19_ni_score.py) | Valuta la correlazione tra il punteggio clinico di invasione neurale (`NI score` da `Table S1A`, `mmc2.xlsx`) e il tasso d'errore di annotazione (ipotesi H4), ricalcolando Holm-Bonferroni sulla famiglia congiunta H1-H4. |
| [`src/20_xai_per_patient.py`](src/20_xai_per_patient.py) | Quantifica la stabilita' delle spiegazioni XAI (coefficienti di Regressione Logistica L1) attraverso i fold LOPO rispetto a un controllo di campionamento casuale di pari numerosita'. |
| [`src/21_fig_blood_scheme.py`](src/21_fig_blood_scheme.py) | Genera la figura e la tabella riassuntiva del tasso di errore secondo lo schema primario a riferimento ematico (`PBMC`). |

---

## 4. Configurazione dell'Ambiente e Riproducibilita'

L'ambiente software completo e' specificato in [`envs/environment.yml`](envs/environment.yml) e bloccato alle versioni esatte in [`envs/conda-lock.yml`](envs/conda-lock.yml) (Python 3.11, `scanpy`, `anndata`, `scikit-learn`, `scipy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pybiomart`, `openpyxl`).

```bash
# Creazione dell'ambiente tramite Mamba / Conda
mamba env create -f envs/environment.yml
mamba activate pdac-ml
```

---

## 5. Esecuzione Sequenziale della Pipeline

Tutti gli script devono essere eseguiti dalla directory radice del repository (`pdac-ml/`). Il seme dei generatori pseudocasuali e' fissato globalmente (`SEED = 0`) in ogni modulo per garantire la riproducibilita' bit-a-bit di partizioni, permutazioni e ricampionamenti bootstrap.

```bash
# 1. Acquisizione metadati GEO, download V(D)J e verifica preliminare delle etichette
python src/00_fetch_metadata.py --out data/raw/metadata
bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr
python src/01_probe_labels.py --raw-dir data/raw
python src/01b_tcr_check.py --raw-dir data/raw

# 2. Costruzione della matrice sparsa AnnData, controllo di qualita' e generazione degli split
python src/02_build_h5ad.py --modality sc --raw-dir data/raw/GSE278688 --out data/interim/sc_raw.h5ad
python src/03_eda_qc.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/04_splits.py --h5ad data/interim/sc_raw.h5ad --out data/interim/splits_sc.csv

# 3. Quantificazione della pseudoreplicazione, confronto modelli, feature selection e XAI
python src/05_baseline.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv
python src/11_model_zoo.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv
python src/12_feature_selection.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv --task task1a_paired
python src/20_xai_per_patient.py --h5ad data/interim/sc_raw.h5ad --splits data/interim/splits_sc.csv --task task1a_paired

# 4. Validazione clonale TCR, controlli di robustezza e arbitrato multi-compartimento
python src/06_clonotype_consistency.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/07_robustness.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/08_disambiguate.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/09_arbitrate.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/10_loco.py --sc data/interim/sc_raw.h5ad --tcr-dir data/raw/tcr
python src/21_fig_blood_scheme.py

# 5. Inferenza CNV e analisi a livello di paziente / espansione clonale
python src/13_infercnv.py --h5ad data/interim/sc_raw.h5ad
python src/14_cnv_loci.py --h5ad data/interim/sc_raw.h5ad
python src/15_inventory.py --raw-dir data/raw
python src/16_pseudobulk_patient.py --h5ad data/interim/sc_raw.h5ad --arbitration results/tables/arbitration_blood.csv
python src/17_clonal_expansion.py --arbitration results/tables/arbitration_blood.csv --pairwise results/tables/pairwise_discordance_clone_b.csv
python src/18_loyo_check.py --input results/tables/clonal_expansion_joined.csv
python src/19_ni_score.py --table-s1 mmc2.xlsx --arbitration results/tables/arbitration_blood.csv --previous-tests results/tables/pseudobulk_hypothesis_tests.csv
```

# Machine learning su scRNA-seq e validazione clonale TCR nel carcinoma duttale pancreatico (PDAC)

Studio computazionale e metodologico sulla SuperSeries NCBI GEO **`GSE278694`** (Chen M-M et al., *Cancer Cell* 2025): **199.184 cellule**, **14 pazienti**, **30.543 geni**, **29 librerie scRNA-seq** (`GSE278688`) e **30 librerie scTCR-seq V(D)J** (`GSE300435`), in tre compartimenti anatomici (`Tumor`: 71.414 cellule, 13 pazienti; `Adjacent_normal`: 38.264 cellule, 5 pazienti; `PBMC`: 89.506 cellule, 11 pazienti). La coorte snRNA-seq (`GSE278689`, **32.005 nuclei**, 8 pazienti) e' esclusa dalle analisi supervisionate: i suoi pazienti sono diversi da quelli scRNA-seq, quindi la modalita' tecnica coinciderebbe con il paziente.

Fonti dei numeri di questo paragrafo: `results/tables/dataset_summary_sc.csv`, `dataset_summary_sn.csv`, `vdj_inventory.csv` (vedi [Struttura dei risultati](#5-struttura-dei-risultati)).

---

## 1. I tre problemi di studio e la sintesi dei risultati

Il lavoro esamina tre problemi di classificazione supervisionata sul PDAC a singola cellula. In ognuno la stima "di pratica corrente", che tratta le cellule come osservazioni indipendenti, e' confrontata con una stima che tiene conto della struttura dei dati, ma il problema e' diverso nei tre casi. Nel compartimento (riga 2) e nella malignita' (riga 3) il divario nasce dal fatto che **l'unita' statistica indipendente e' il paziente, non la cellula**. Nel tipo cellulare (riga 1) il divario fra split casuale e split per paziente e' piccolo: il problema principale e' la circolarita' delle etichette, che resta in entrambe le stime.

| Problema | Popolazione e dimensione | Stima di pratica corrente | Stima che rispetta la struttura di gruppo | Origine dell'artefatto |
| :--- | :--- | :--- | :--- | :--- |
| **1. Tipo cellulare** (`task2_v1`, 9 classi) | solo tessuto tumorale: 70.547 cellule, 13 pazienti | split casuale: macro-F1 **0,964** (logistica L2) | split per paziente: macro-F1 **0,954** (logistica L2, 5 fold per paziente; Δ = −0,009); LightGBM 0,962 (LOPO) | **Circolarita'**: le etichette di riferimento derivano da clustering e marcatori sulla stessa matrice usata per addestrare i classificatori. Lo split per paziente non la corregge: l'accuratezza e' reale ma misura la riproduzione dell'annotazione, non una verita' indipendente. |
| **2. Compartimento** (`task1a_paired`, Tumor vs Adjacent_normal) | 67.662 cellule, 5 pazienti con entrambi i tessuti | split casuale: macro-F1 **0,875 ± 0,002** | LOPO: macro-F1 **0,729 ± 0,100** (Δ = −0,146; deviazione standard fra fold 50 volte piu' grande) | **Struttura per campione**: ogni libreria corrisponde a un paziente x tessuto, e lo split casuale mette in training e test cellule della stessa libreria; il modello puo' sfruttarne caratteristiche che non si trasferiscono a un paziente nuovo (meccanismo non misurato: il disegno e' incrociato, ogni paziente ha entrambi i tessuti). |
| **3. Malignita' delle duttali** (`13_infercnv.py`, `14_cnv_loci.py`) | 8.370 duttali tumorali (13 pazienti) e 334 adiacenti (5 pazienti); entro paziente i 4 pazienti valutabili, con almeno 10 duttali per tessuto (4.870 e 333 duttali) | sui loci CNV ricorrenti del PDAC, AUC tumore vs adiacente pooled **0,717** | entro paziente **≈ 0,57** (0,53-0,61 nei 4 pazienti) | **Etichetta derivata dal trascrittoma e confondimento con il paziente**: il contrasto tumore-adiacente pooled mescola differenze fra pazienti; l'etichetta di malignita' viene da inferCNV sugli stessi dati (vedi le note). |

Note alla tabella:
- Riga 1: 0,954 e 0,964 vengono da `05_baseline.py` (5 fold per paziente e 5 fold casuali); 0,962 da `11_model_zoo.py` in leave-one-patient-out.
- Riga 3: un Random Forest addestrato sull'etichetta di malignita' derivata da inferCNV la riproduce con AUC LOPO **0,836 ± 0,062** (duttali tumorali di 12 pazienti) e predice maligne il 37,4% delle 334 duttali adiacenti; questa AUC misura la riproducibilita' dell'etichetta, non la malignita'. Lo score CNV medio, senza i loci, ha AUC pooled 0,598 ed entro paziente 0,616 (media dei 4 pazienti): entro paziente i loci PDAC non fanno meglio dello score medio. Il filtro del rumore di inferCNV usa una sola soglia su tutte le cellule (non dipende dall'ordine delle cellule) e ogni finestra ha le coordinate reali dei suoi geni (`cnv_windows.csv`); su 9p, povero di geni, le 10 finestre che contengono CDKN2A coprono 15-28 Mb ciascuna, quindi il segnale "9p21" e' in pratica quello del braccio.
- Riga 3: in inferCNV il riferimento diploide sono le cellule immunitarie di **tutti** i pazienti (una sola esecuzione); solo la soglia (mediana + 3·MAD degli immuni) e' calcolata per paziente, ed e' stimata sulle stesse cellule su cui inferCNV centra lo score. Questo puo' contribuire al fatto che anche il 38-65% delle duttali adiacenti dei pazienti con almeno 20 duttali adiacenti risulta "maligno" (`malignancy_labels.csv`); i dati non distinguono fra artefatto della soglia, cellule tumorali nel campione adiacente ed effetto di campo. La sovrapposizione fra duttali tumorali e adiacenti e' **compatibile** con un effetto di campo di cancerizzazione, ma con questi limiti (riferimento comune, soglia sulle stesse cellule, confini dei loci approssimati, 4 pazienti valutabili) e' un'interpretazione, non una dimostrazione.

### Confronto dei modelli, selezione delle feature e spiegazioni per paziente

* **13 classificatori (`11_model_zoo.py`, LOPO, 20.000 cellule di training per fold)**. Su `task2_v1` i modelli lineari arrivano a 0,936-0,943 contro 0,962 di LightGBM. Su `task2_hier3` (CD4T / CD8T / NK, 42.168 cellule) la logistica L2 (0,850) e' il riferimento piu' debole dopo il dummy: la superano i modelli non lineari, LightGBM 0,923 (+7,2 punti), SVM-RBF 0,906 (+5,6), MLP 0,904 (+5,3), Random Forest 0,901 (+5,1), ma anche lineari come la logistica L1 0,877 (+2,7) e la SVM lineare 0,865 (+1,5). La logistica L1 e l'elastic net (SGD, 50 epoche) non arrivano a convergenza in nessun fold; gli altri modelli iterativi convergono in tutti (`n_folds_not_converged`). Su `task1a_paired` i modelli vanno da 0,496 (dummy) a 0,725 (SVM-RBF); nessuna differenza dal migliore e' significativa al t corretto di Nadeau-Bengio, nemmeno quella del dummy (Δ = −0,229, p = 0,072; esclusi il dummy, p ≥ 0,080). I p sono grezzi e il migliore e' scelto sugli stessi fold. Con 5 fold la potenza e' bassa: un p non significativo **non** dimostra che i modelli siano equivalenti.
* **Selezione delle feature (`12_feature_selection.py`)**. Tutti i selettori di un fold usano le stesse 20.000 cellule di training (sottocampione casuale uniforme del solo training). A k = 50 la logistica L1 ottiene 0,936 su `task2_v1` e 0,712 su `task1a_paired`, il Random Forest 0,926 e 0,735, i geni piu' dispersi (HVG) 0,659 e 0,623 (k scelti a caso: 0,367 e 0,589). Su `task2_v1` la logistica L1 non converge in 400 iterazioni in nessun fold e ha 524-583 coefficienti non nulli, quindi a k = 1000 una parte dei geni scelti ha coefficiente nullo (`feature_selection_diag_task2_v1.csv`); su `task1a_paired` converge (177-268 iterazioni). La stabilita' (Jaccard medio fra fold: i 5 fold per paziente di `fold_group_5` su `task2_v1`, i 5 fold LOPO su `task1a_paired`) e' piu' bassa per L1 (0,480 a k = 200 su `task2_v1`, da fit che non convergono; 0,254 a k = 1000 su `task1a_paired`, dove L1 converge) che per i filtri univariati (chi-quadro 0,842 e 0,678; ANOVA 0,794 e 0,665; HVG 0,793 e 0,771). Su `task1a_paired`, a k = 50, i geni scelti in tutti e 5 i fold sono 18 per HVG, 19 per chi-quadro, 19 per ANOVA, 16 per il Random Forest e 4 per L1 (`feature_core_genes_task1a_paired.csv`).
* **Spiegazioni per paziente (`20_xai_per_patient.py`)**. Sui 5 fold LOPO di `task1a_paired`, il Jaccard medio fra i top-50 geni di una logistica L1 e' **0,300** (0,220-0,408 fra coppie di fold; Spearman medio fra i vettori dei coefficienti su tutti i geni di `var`, con 0 per i geni esclusi dal prefiltro, 0,234); 9 geni sono nei top-50 di tutti i fold: BTG2, CD69, CREM, CTRB1, DUSP4, FOS, GAPDH, HSPH1, VIM. Il controllo, con sottoinsiemi casuali di 20.000 cellule (la stessa dimensione dei training LOPO sottocampionati) senza struttura per paziente, da' **0,459** (differenza +0,159). Togliere un paziente cambia la regola appresa piu' che addestrare su un sottoinsieme casuale della stessa dimensione, come atteso se le cellule dello stesso paziente sono ridondanti. Due cautele: in simulazione (`src/sim/sim_d_xai_stability.py`, 5 pazienti, 5 repliche per scenario) la differenza era in media +0,143 (da +0,116 a +0,181) con un effetto paziente costruito e +0,030 (da −0,013 a +0,052) senza; il valore osservato, +0,159, cade nell'intervallo dello scenario con effetto paziente, ma la simulazione ha una sola forma di effetto e pochi pazienti, quindi e' un'indicazione e non una calibrazione; dei 9 geni comuni, FOS, BTG2 e HSPH1 appartengono alla firma dello stress da dissociazione; CD69, DUSP4 e CREM sono geni di risposta immediata che marcano anche attivazione, residenza tissutale (CD69) ed esaurimento (DUSP4) delle T infiltranti; CTRB1 e' un marcatore acinare (probabile RNA ambientale); GAPDH e VIM non hanno un'interpretazione specifica: non sono necessariamente biologia tumore/normale.

### Validazione clonale con il TCR

Il riarrangiamento V(D)J delle catene α/β precede l'espansione periferica e non dipende dallo stato trascrizionale: la sequenza CDR3 identifica un clone indipendentemente dal clustering. Il lignaggio CD4/CD8 di un clone e' di norma stabile, con eccezioni note (CD4 citotossiche, CD4+CD8αα+ e doppi positivi nei tessuti, MAIT, NKT, γδ); con questa riserva il confronto fra compartimenti permette di misurare l'affidabilita' dell'annotazione.

* **Rumore di base ed eccesso fra compartimenti (`07_robustness.py`)**. Due cellule dello stesso clone nello stesso compartimento hanno etichette diverse con probabilita' **0,107** (rumore di annotazione di base, stimatore senza reinserimento). Fra compartimenti diversi la discordanza e' piu' alta: **eccesso +0,059** (IC 95% cluster bootstrap sui pazienti **+0,031; +0,110**, 11 pazienti, 1.114 coppie clone-compartimenti), simile con la coppia TRA+TRB come definizione di clone (`clone_ab`: +0,065, IC +0,035; +0,121; il valore principale usa la sola catena TRB, `clone_b`). Con pesi uguali per i due compartimenti (½‖p_A − p_B‖², stesso bootstrap) l'eccesso e' +0,074 (IC +0,050; +0,119). La versione esplorativa split-half di `06` da' +0,063 (IC cluster bootstrap +0,023; +0,134, ricampionando 14 pazienti: 11 contribuiscono al termine fra compartimenti, tutti e 14 al basale; escluse 90 coppie e 3.943 repliche con un pareggio nell'etichetta di maggioranza).
* **Tasso d'errore con riferimento solo sangue (`10_loco.py`, `21_fig_blood_scheme.py`)**. Con l'identita' CD4/CD8 del clone stimata dai marcatori del solo PBMC, l'etichetta di maggioranza del compartimento e' sbagliata nel **17,0%** dei cloni nel `Tumor` (**0,170**, IC 95% **0,095; 0,252**, 277 cloni, 10 pazienti) e nel **6,0%** nell'`Adjacent_normal` (0,060, IC non calibrato 0,000; 0,095, 117 cloni, 4 pazienti). Un clone-compartimento in cui due etichette hanno lo stesso numero di cellule non ha un'etichetta di maggioranza ed e' escluso (10 nel Tumor, 1 nell'Adjacent_normal); contando questi pareggi come errori l'errore nel Tumor e' 0,199 (IC 0,117; 0,291), come corretti 0,164 (IC 0,093; 0,240) (`error_rate_tie_sensitivity`). La differenza appaiata sugli stessi cloni e' **+0,059** (IC **0,000; +0,160**, 85 cloni, 4 pazienti; +0,068 e +0,057 con i pareggi contati come errori o come corretti): con 4 pazienti esistono solo 35 ricampionature distinte e l'IC non e' calibrato; il 93,3% delle repliche bootstrap e' > 0. Scomposizione dell'errore nel Tumor: 32 cloni CD4 nel sangue etichettati CD8T (0,116), 15 cloni etichettati NK (0,054: 11 CD8 e 4 CD4 nel sangue), 0 cloni CD8 etichettati CD4T; escludendo le etichette NK l'errore e' 0,122 (IC 0,049; 0,224). Le cellule NK con TCR dei 15 cloni etichettati NK hanno d1 = 0,617 (IC 0,058; 0,849) e d2 = 0,458 (IC 0,008; 0,653) (`doublet_mixture_error_nk`): sono compatibili con una quota di doppietti T+NK, ma l'IC e' troppo largo per dire che siano in gran parte doppietti; questi errori possono quindi essere in parte doppietti e non errori di annotazione. Nelle cellule etichettate CD8T del tumore, CD8A e' rilevato nel **7%** delle cellule dei cloni discordanti contro l'**88%** dei concordanti (`08_disambiguate.py`; frazione piu' bassa nei discordanti in 5 pazienti valutabili su 5; percentuali per cellula, descrittive), e CD4 medio e' 16 volte piu' alto (rapporto fra medie di log1p(CP10k), non un fold-change). I dati sono compatibili con cellule CD4 etichettate CD8T; resta aperta l'alternativa di un fenotipo CD4 citotossico (GNLY 1,78 nei discordanti contro 0,44 nei concordanti); la bassa frequenza di CD8A nei discordanti rende improbabile che siano CD4+CD8αα+. Fra i cloni del Tumor, i 12 MAIT-like (TRAV1-2 con TRAJ33/12/20) sono tutti concordanti e non ci sono cloni iNKT (TRAV10-TRAJ18) (`clone_invariant_tcr_check`): le eccezioni MAIT/NKT non spiegano gli errori. Con lo schema che esclude il compartimento giudicato (LOCO) la differenza Tumor − Adjacent_normal e' +0,016 (IC −0,020; +0,095) e con lo schema ingenuo +0,024 (IC 0,000; +0,110; questo schema include anche i cloni presenti in un solo compartimento, a differenza di `09`): la stima puntuale cambia con lo schema di riferimento, e con 4-5 pazienti nessuno di questi IC e' calibrato. Con regole d'identita' CD8 piu' strette l'errore nel Tumor resta simile: 0,180 richiedendo CD8B > 0 (IC 0,099; 0,271) e 0,181 richiedendo CD8A ≥ 2 UMI (IC 0,098; 0,269).
* **Doppietti fra le NK con TCR (`09_arbitrate.py`, `10_loco.py`)**. Fra le cellule etichettate come non-T con una catena TRB produttiva (fra cui NK 5,9%, Endothelial 5,6%, Stellate 5,5%, Platelet 13,8% delle cellule della classe), due misure fisiche (UMI per goccia e UMI del TCR), confrontate fra classi dell'annotazione, sono compatibili, nella stima pooled su tutti i tessuti, con una **maggioranza di doppietti T+NK** fra le NK-TCR+: dall'eccesso di UMI totali d1 = **0,723** (IC 0,556; 0,848) e dalla diluizione degli UMI del TCR d2 = **0,676** (IC 0,559; 0,781). Il modello di miscela e' lo stesso per le due stime (non-doppietto: una cellula con TCR alla densita' di una T e gli UMI di una NK; doppietto: NK + CD8T) ed e' applicato alle mediane, quindi e' un'approssimazione. Le due stime non sono indipendenti (usano le stesse mediane di n_counts), sono pooled fra tessuti e la copertura dei loro IC non e' verificata: per compartimento vanno da 0,51 (Tumor) a 0,84 (PBMC) per d1 e da 0,42 a 0,85 per d2; con le medie al posto delle mediane d1 = 0,682 e d2 = 0,693. Nel Tumor, il compartimento dell'analisi d'errore, la maggioranza non e' sostenuta: d1 = 0,510 (IC 0,186; 0,658) e d2 = 0,419 (IC 0,289; 0,567) con le mediane, d1 = 0,381 (IC 0,170; 0,473) con le medie. Nella versione precedente d2 era 0,790, per un errore nel denominatore del rapporto atteso nei doppietti.
* **Espansione clonale ed errore (`17_clonal_expansion.py`, `18_loyo_check.py`)**. Su 287 cloni condivisi fra sangue e tumore (10 pazienti; 277 senza pareggi per l'errore) la dimensione clonale e' **associata all'errore entro paziente**: errore binario ρ = +0,220 (p = 0,0004, block permutation test entro paziente con 20.000 permutazioni Monte Carlo), eccesso continuo ρ = +0,281 (p = 0,00005, il minimo raggiungibile con 20.000 permutazioni). Il p e' condizionato al paziente ma il ρ e' pooled: la distribuzione nulla ha media +0,021 (errore) e +0,027 (eccesso), e il p sulla statistica centrata e' 0,00055 e 0,00005. Con 5 semi il p dell'errore va da 0,00015 a 0,0004; con 100.000 permutazioni e' 0,00018. Fra le medie per paziente il segno e' lo stesso ma l'associazione non e' significativa (errore: ρ = +0,467, p = 0,18; eccesso: ρ = +0,321, p = 0,37; n = 10). Togliendo un paziente alla volta (leave-one-patient-out) l'associazione resta con p < 0,05 in **10/10** casi per entrambi gli esiti (ρ fra 0,195 e 0,252 per l'errore; fra 0,236 e 0,368 per l'eccesso; 20.000 permutazioni, come in 17), anche con il p centrato; per l'errore anche con 5 semi o 100.000 permutazioni (p fino a 0,0016).
* **Quattro covariate per paziente (`16_pseudobulk_patient.py`, `19_ni_score.py`)**. Per spiegare le differenze di errore fra pazienti (n = 10) senza cercare fra migliaia di geni, sono state testate quattro ipotesi, formulate dopo i moduli 08-10 ma prima di calcolare le covariate (H3 indicata fin dall'inizio come speculativa), con Spearman, permutazione Monte Carlo (20.000) e Holm-Bonferroni sulla famiglia di 4 test. Nessuna e' significativa e tutte hanno segno opposto a quello ipotizzato:
  1. `H1` esaurimento delle CD8T tumorali: ρ = −0,45, p = 0,190, p_Holm = 0,58. Limiti: la deriva CD8 → CD4/NK ipotizzata va nel senso opposto alla maggior parte degli errori osservati (32 cloni CD4 → CD8T e 0 CD8 → CD4; solo gli 11 CD8 → NK vanno nel senso ipotizzato, `error_decomposition_blood`), e il punteggio e' calcolato sulle cellule etichettate CD8T, cioe' sulle stesse etichette il cui errore e' l'esito;
  2. `H2` frazione mieloide nel tumore: ρ = −0,50, p = 0,145, p_Holm = 0,58;
  3. `H3` frazione di CAF nel tumore: ρ = −0,13, p = 0,724, p_Holm = 1,00;
  4. `H4` NI score clinico (Table S1A di Chen et al.): ρ = −0,13, p = 0,735, p_Holm = 1,00.

  Con n = 10 questi test hanno poca potenza: "non significativo" non esclude un effetto.

---

## 2. Moduli sorgente (`src/`, `src/sim/`, `tests/`)

Gli script lavorano su matrici sparse `scipy.sparse.csr_matrix` `float32`; le trasformazioni stimate dai dati sono dentro `Pipeline` di scikit-learn o funzioni che ricevono il solo training (verificato in `tests/test_leakage.py`). Ogni script eseguibile ha una docstring con scopo, input, output, argomenti e uso, e risponde a `--help` (`scpipe.py` e' un modulo di funzioni comuni, non uno script).

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

**Ambiente.** Definito in [`envs/environment.yml`](envs/environment.yml) (ambiente `pdac`, Python 3.11, versioni fissate a quelle che hanno prodotto `results/`) e bloccato con checksum in [`envs/conda-lock.yml`](envs/conda-lock.yml). Tra i pacchetti: `scanpy`, `anndata`, `scikit-learn`, `lightgbm`, `infercnvpy`, `requests` (coordinate dei geni da BioMart), `openpyxl`, `scipy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pytest`, `ruff`.

```bash
mamba env create -f envs/environment.yml      # oppure: conda-lock install -n pdac envs/conda-lock.yml
mamba activate pdac
```

**Dati.** Tutti i comandi si eseguono dalla cartella radice del repository, `progetto-ml/`.
- Matrici scRNA-seq e snRNA-seq: `src/02_build_h5ad.py` scarica da GEO i file supplementari di `GSE278688` (scRNA) e `GSE278689` (snRNA) in `data/raw/GSE278688/` e `data/raw/GSE278689/` (matrix, barcodes, features, metadata; ~1,2 GB). Sono le SubSeries della SuperSeries `GSE278694`: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278694>.
- File VDJ: `src/00_fetch_metadata.py` scarica `filelist.txt` di `GSE278694`; `src/download_tcr.sh` scarica i 30 file `*_VDJ_filtered_contig_annotations.csv.gz` dei campioni di `GSE300435` in `data/raw/tcr/`.
- Metadati clinici: `mmc2.xlsx` e' la Table S1 del materiale supplementare di Chen et al. 2025 (pagina dell'articolo su *Cancer Cell*, doi sotto, sezione "Supplemental information"). Va scaricato a mano e messo in `data/raw/paper_supplement/mmc2.xlsx`; `19_ni_score.py` ne legge il foglio `Table S1A`.

**Citazione del dataset.** Chen M-M et al. *Cancer Cell* 2025;43(9):1656-1676.e10. doi:[10.1016/j.ccell.2025.06.020](https://doi.org/10.1016/j.ccell.2025.06.020).

**Casualita'.** Tutti gli script usano semi fissi (`SEED = 0`; nei sottocampionamenti per fold di 11 e 12, `SEED + fold`). Sottocampionamento: casuale uniforme, solo sul training, non stratificato (20.000 cellule per fold in 11, 12 e 20). Cluster bootstrap sui pazienti: 2.000 repliche (06, 07, 08, 09, 10, 21), 500 per d1 e d2 (10), 2.000 per gli IC di Spearman (16, 17, 19); intervallo "expanded percentile" (vedi [§6](#6-simulazioni-e-test)). Permutazioni Monte Carlo: 20.000 (16, 17, 18, 19), 100.000 nelle sensibilita' di 17 e 18.

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
python src/sim/sim_d_xai_stability.py    # prima di 20: la lettura finale di 20 la confronta
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
pytest                      # tutti i test (circa 3 minuti); solo i veloci: pytest -m "not slow"
ruff check src tests revisione
python revisione/genera_risultati.py --strict   # revisione/RISULTATI_FINALI.md dalle tabelle
```

Note: `07` con `--clone-def clone_b` va eseguito prima di `08` e `17`, che leggono `pairwise_discordance_clone_b.csv` (l'ordine fra le due chiamate di `07` e' indifferente: ognuna scrive file con il proprio suffisso); `13` e `14` riusano le cache `data/interim/cnv_scores.csv` e `cnv_matrix.npz` solo se coprono esattamente le cellule selezionate e la loro impronta in `data/interim/cnv_cache_fingerprint.json` (cellule, conteggi, geni, coordinate, riferimento, parametri) coincide; altrimenti si fermano e le cache vanno cancellate per ricalcolare inferCNV; inferCNV usa solo le finestre dei cromosomi 1-22 (esclusi MT, la cui espressione segue la frazione mitocondriale, i contig non standard, X e Y); le coordinate dei geni sono quelle di Ensembl 116, scaricate una volta via https da BioMart in `data/interim/gene_positions_ensembl116.csv` e poi rilette da li' (serve la rete alla prima esecuzione); `--method copykat` richiede `rpy2` e R, non inclusi nell'ambiente; `08` riusa `results/tables/doublet_scores.csv` (Scrublet) solo se contiene esattamente le cellule dell'h5ad, altrimenti si ferma. `11` addestra i 13 modelli in circa 40 minuti su `task2_v1` (16 per il solo LightGBM), 18 su `task2_hier3` e 9 su `task1a_paired` (somma dei tempi per modello stampati nei log `revisione/log/11_final/*.log`, 4-6 thread).

---

## 5. Struttura dei risultati

| Script | Tabelle in `results/tables/` | Figure in `results/figures/` (`.png`) |
| :--- | :--- | :--- |
| 03 | `dataset_summary_{sc,sn}`, `qc_by_library_{sc,sn}`, `composition_{sc,sn}`, `patient_support_{sc,sn}` | `qc_by_library_{sc,sn}`, `qc_complexity_{sc,sn}`, `composition_by_tissue_{sc,sn}`, `patient_by_celltype_{sc,sn}` |
| 04 | `task_summary_sc` (e `data/interim/splits_sc.csv`) | - |
| 05 | `baseline_sc` (per fold), `baseline_summary_sc` (media e sd), `baseline_gap_sc` (divario casuale - paziente e rapporto delle deviazioni standard) | `gap_sc_{task}` (6 task) |
| 06 | `nk_contamination`, `tcr_rate_by_celltype`, `clone_switching_across_clone_b`, `clone_switching_summary_clone_b`, `clone_switching_within_clone_b` (non versionato: una riga per replica split-half di ogni clone e compartimento) | `clone_label_switching_clone_b` |
| 07 | `pairwise_discordance_{clone_b,clone_ab}`, `excess_summary_{cc}`, `excess_by_compartment_pair_{cc}`, `discordant_clone_markers_{cc}`, `ambient_discriminant`, `chat_sensitivity_{raw,hc}` | `excess_by_pair_{cc}`, `ambient_vs_real_tcr` |
| 08 | `triple_clones`, `triple_clones_summary`, `markers_disc_vs_conc`, `marker_test_disc_vs_conc`, `doublet_tcr_pos_vs_neg`, `doublet_scores` (cache di Scrublet, non versionata) | `doublet_by_tcr_status` |
| 09 | `arbitration`, `error_rate_by_compartment`, `doublet_physical`, `arbitration_tie_sensitivity` | `doublet_physical`, `error_rate_by_compartment` |
| 10 | `arbitration_{naive,loco,blood}`, `error_rate_schemes`, `error_rate_paired_diff`, `error_decomposition_blood`, `error_rate_identity_rules`, `ambient_cd4_threshold`, `doublet_mixture`, `doublet_mixture_boot`, `error_rate_tie_sensitivity`, `clone_invariant_tcr_check`, `doublet_mixture_error_nk` | `error_rate_loco` |
| 21 | `error_rate_blood` | `error_rate_blood` |
| 11 | `model_zoo_{sc,task1}` (per fold), `model_zoo_summary_{sc,task1}`, `model_comparison_{task2_v1,task2_hier3,task1a_paired}` | `model_zoo_{task2_v1,task2_hier3,task1a_paired}` |
| 12 | `feature_selection_{task}`, `feature_selection_summary_{task}`, `feature_stability_{task}`, `feature_selection_diag_{task}`, `feature_core_genes_{task}` | `feature_selection_{task}` |
| 20 | `xai_jaccard_task1a_paired`, `xai_spearman_task1a_paired`, `xai_summary_task1a_paired`, `xai_common_genes_task1a_paired` | `xai_per_patient_task1a_paired` |
| 13 | `malignancy_labels`, `malignancy_lopo`, `malignancy_lopo_summary`, `external_test_predictions`, `external_test_summary` | `malignancy` |
| 14 | `cnv_loci_auc`, `cnv_loci_auc_by_patient`, `cnv_windows` | `cnv_loci_auc`, `cnv_loci_separation` |
| 15 | `obs_inventory`, `vdj_inventory` | - |
| 16 | `pseudobulk_patient_table`, `pseudobulk_hypothesis_tests` | `pseudobulk_patient_pca` |
| 17 | `clonal_expansion_joined`, `clonal_expansion_tests`, `clonal_expansion_perm_sensitivity` | `clonal_expansion_vs_{error,excess}` |
| 18 | `loyo_clonal_expansion_{error,excess}`, `loyo_clonal_expansion_perm_sensitivity` | `loyo_clonal_expansion_{error,excess}` |
| 19 | `ni_score_joined`, `ni_score_test`, `ni_score_holm_joint` | `ni_score_vs_error` |
| sim | `sim_a_nadeau_bengio`, `sim_a_cluster_bootstrap`, `sim_b_reference_schemes`, `sim_c_block_permutation`, `sim_d_xai_stability` | - |

Da dove viene ogni numero del README:

| Numero | Tabella (colonna / riga) |
| :--- | :--- |
| 199.184 cellule, 14 pazienti, 30.543 geni, 29 librerie, cellule e pazienti per tessuto | `dataset_summary_sc` |
| 32.005 nuclei, 8 pazienti | `dataset_summary_sn` |
| 30 librerie VDJ | `vdj_inventory` (30 righe) |
| 70.547 / 67.662 / 42.168 cellule, 13 / 5 pazienti dei task | `task_summary_sc` |
| 0,954; 0,964; 0,875 ± 0,002; 0,729 ± 0,100 | `baseline_summary_sc` (`logreg_l2`, `macro_f1`) |
| Δ = −0,009 (riga `task2_v1`) e −0,146 (segno opposto di `gap_random_minus_honest`), deviazione standard 50 volte piu' grande | `baseline_gap_sc` (`gap_random_minus_honest`, `sd_ratio_honest_over_random`, riga `task1a_paired`) |
| 0,936-0,943, 0,962, 0,850, 0,923, 0,906, 0,904, 0,901, 0,877, 0,865, 0,496-0,725; +7,2, +5,6, +5,3, +5,1, +2,7, +1,5 punti; convergenza | `model_zoo_summary_sc`, `model_zoo_summary_task1` (`macro_f1_mean`, `delta_vs_logreg_l2`, `n_folds_not_converged`) |
| Δ = −0,229, p ≥ 0,080, p = 0,072 | `model_comparison_task1a_paired` (`delta_vs_best`, `p_nadeau_bengio`) |
| macro-F1 e Jaccard della selezione delle feature | `feature_selection_summary_{task}` (`macro_f1_mean`), `feature_stability_{task}` |
| convergenza e coefficienti non nulli di L1; geni comuni ai 5 fold | `feature_selection_diag_{task}`, `feature_core_genes_task1a_paired` |
| 0,300; 0,459; +0,159; 0,234; 9 geni | `xai_summary_task1a_paired`, `xai_common_genes_task1a_paired`; 0,220-0,408 da `xai_jaccard_task1a_paired` |
| 0,107; +0,059 [+0,031; +0,110]; +0,074 [+0,050; +0,119]; 11 pazienti; 1.114 coppie | `excess_summary_clone_b`; +0,065 [+0,035; +0,121] da `excess_summary_clone_ab` |
| +0,063 [+0,023; +0,134]; 14, 11 e 14 pazienti; 90 coppie e 3.943 repliche escluse | `clone_switching_summary_clone_b` (`n_patients`, `n_patients_across`, `n_patients_within`, `n_pairs_ties_excluded`, `n_reps_ties_excluded`) |
| +0,048; +0,070 (bootstrap sulle coppie); p = 0,00045 e 0,0004 (permutazione dei cloni e block permutation) | `excess_summary_clone_b` (`lo_naive`, `hi_naive`); `clonal_expansion_tests` (righe `clone_naive_permutation` e `clone_block_permutation`, esito `error`) |
| 0,170 [0,095; 0,252]; 0,060 [0,000; 0,095]; 277 e 117 cloni; 10 e 1 pareggi | `error_rate_blood` (`n`, `n_ties`) |
| 0,199 [0,117; 0,291], 0,164 [0,093; 0,240]; +0,068 e +0,057 (pareggi contati come errori / corretti) | `error_rate_tie_sensitivity` (schema `blood`) |
| +0,059 [0,000; +0,160]; 85 cloni; 35 ricampionature; 93,3%; +0,016 [−0,020; +0,095]; +0,024 [0,000; +0,110] | `error_rate_paired_diff` |
| 32, 15, 11, 4, 0 cloni; 0,116; 0,054; 0,122 [0,049; 0,224] | `error_decomposition_blood` |
| d1 0,617 [0,058; 0,849], d2 0,458 [0,008; 0,653]; 15 cloni NK in errore | `doublet_mixture_error_nk` (riga `median`) |
| 12 cloni MAIT-like tutti concordanti, 0 iNKT | `clone_invariant_tcr_check` |
| 0,180 [0,099; 0,271] (CD8B > 0); 0,181 [0,098; 0,269] (CD8A ≥ 2) | `error_rate_identity_rules` (riga Tumor) |
| 11 cellule (soglia ambientale di CD4 nel PBMC) | `ambient_cd4_threshold` (riga PBMC, `n_cells`) |
| 7% / 88%; 5 su 5 pazienti; 16×; GNLY 1,78 / 0,44 | `marker_test_disc_vs_conc` (riga Tumor) |
| NK 5,9%, Endothelial 5,6%, Stellate 5,5%, Platelet 13,8% | `tcr_rate_by_celltype` |
| d1 0,723, d2 0,676 e varianti (medie, compartimenti) | `doublet_mixture` |
| ρ, p, media nulla, p centrato, 287 cloni; medie per paziente (ρ = +0,467 e +0,321) | `clonal_expansion_tests`, `clonal_expansion_perm_sensitivity` |
| 10/10, intervalli di ρ, p centrato | `loyo_clonal_expansion_{error,excess}` (`p_without`, `p_centered_without`); p fino a 0,0016 da `loyo_clonal_expansion_perm_sensitivity` |
| H1-H3 | `pseudobulk_hypothesis_tests`; H4 `ni_score_test`; p di Holm `ni_score_holm_joint` |
| 0,836 ± 0,062; 12 pazienti | `malignancy_lopo_summary` (riga `auc`) |
| 37,4%; 334 duttali | `external_test_summary` (riga `_all`) |
| 8.370 e 334 duttali, 13 e 5 pazienti | `malignancy_labels` (somme di `n` e pazienti per tessuto) |
| 38-65% | `malignancy_labels` (Adjacent_normal, pazienti con n ≥ 20) |
| 0,598; 0,717 | `cnv_loci_auc` |
| 10 finestre con CDKN2A, 15-28 Mb ciascuna | `cnv_windows` |
| ≈ 0,57; 0,53-0,61; 0,616; 4 pazienti; 4.870 e 333 duttali | `cnv_loci_auc_by_patient` (riga `_mean`: media delle AUC dei loci e dello score medio, colonne `auc_PDAC_COMPOSITE` e `auc_MEAN_\|CNV\|`; somme delle cellule) |
| numeri delle simulazioni | `sim_*` |

---

## 6. Simulazioni e test

**Simulazioni (`src/sim/`)**, eseguite con semi fissi; i risultati sono in `results/tables/sim_*.csv`.
- (a) `sim_a_calibration.py`: t di Nadeau-Bengio sotto H0 (due modelli dummy) e con effetto noto (disegno ripreso da `genomic-audit-mvp`, stesso autore); cluster bootstrap con 4, 5, 10 e 11 pazienti, con l'intervallo percentile semplice e con quello "expanded percentile". Risultati: sotto H0 Nadeau-Bengio da' p < 0,05 nel 2,7% (k = 3) e 1,8% (k = 5) delle repliche (t non corretto: 4,7% e 6,1%), con effetto noto nel 99,3%; l'IC espanso del cluster bootstrap copre il 93,5-94,8% con 10-11 pazienti (percentile semplice 88,8-92,3%), il 91,5-92,5% con 5 pazienti e solo l'82-83% con 4. La calibrazione di Nadeau-Bengio e' verificata (simulazione e test con banda) con fold stratificati di uguale dimensione (k = 3 e 5, n = 90); in `11` il test e' applicato a fold LOPO di dimensione diversa con training sottocampionato, caso non coperto da questa verifica. Il cluster bootstrap e' verificato con 1.000 repliche per IC (le analisi ne usano 2.000).
- (b) `sim_b_reference_schemes.py`: i tre schemi di riferimento di `10_loco.py` su cloni simulati, sotto H0 e con un artefatto atteso nel tumore (perdita del segnale CD8 ed espressione ambientale di CD4). Risultati (stima ed errore vero calcolati sulle stesse unita', esclusi i pareggi dell'etichetta di maggioranza): sotto H0 i tre schemi stimano l'errore vero entro 0,005; con l'artefatto lo schema che esclude il compartimento giudicato (LOCO) gonfia l'errore del sangue a 0,119 (errore vero 0,086) con artefatto moderato e a 0,383 (vero 0,086; fino a 0,439 in una replica) con artefatto forte, mentre gli schemi solo sangue e naive restano entro 0,008 dall'errore vero. Nella simulazione etichette e marcatori sono indipendenti dato il lignaggio: la circolarita' dello schema naive non e' simulata, quindi il suo basso bias qui non lo convalida sui dati reali. I parametri sono stati scelti a priori: il valore 0,24 citato nella versione precedente non e' riprodotto (0,119 e 0,383 con le due intensita' dell'artefatto).
- (c) `sim_c_block_permutation.py`: falsi positivi del block permutation test sotto H0 con confondente di paziente, 1.000 repliche per livello (IC binomiale esatto al 95%). Risultati: il p non centrato e' < 0,05 nel 5,0% [3,7; 6,5], 5,9% [4,5; 7,5], 4,8% [3,6; 6,3] e 5,1% [3,8; 6,7] delle repliche (confondente di forza 0; 0,5; 1; 2); il p centrato fra il 4,1% e il 5,1%; la permutazione che ignora i pazienti fra il 5,1% e il 99,4%. Nelle prime 100 repliche di ogni livello (il disegno da 100 repliche della versione precedente) i falsi positivi del p non centrato (quello usato in 17 e 18) sono 4, 6, 11 e 5: il dato "2 falsi positivi su 100" non e' riprodotto (con il p centrato e confondente 2 le prime 100 repliche danno 2), e con 100 repliche il conteggio oscilla per il solo rumore binomiale.
- (d) `sim_d_xai_stability.py`: la procedura di `20` su conteggi simulati con e senza effetto paziente. Risultati: differenza controllo - LOPO +0,143 con effetto paziente e +0,030 senza (5 repliche per condizione; nella versione precedente erano citati 0,063 e 0,040). La differenza senza effetto non e' nulla anche perche' i sottoinsiemi casuali del controllo condividono piu' cellule fra repliche di quante ne condividano due fold LOPO.

**Colonne con un'unita' diversa dal paziente.** Alcune tabelle conservano, solo come confronto, la versione con l'unita' sbagliata: `lo_naive`/`hi_naive` in `excess_summary_{cc}` (bootstrap sulle coppie clone-compartimenti: IC piu' stretto, +0,048; +0,070 contro +0,031; +0,110 sui pazienti) e la riga `clone_naive_permutation` di `clonal_expansion_tests` (permutazione dei cloni ignorando il paziente: per l'errore p = 0,00045 contro 0,0004 del block permutation test, per l'eccesso entrambi 0,00005, il minimo). Servono a mostrare quanto cambia il risultato ignorando la dipendenza fra cellule e cloni dello stesso paziente: con questi dati la differenza si vede nell'IC dell'eccesso, non nel p dell'espansione clonale; i numeri del README usano sempre il paziente come unita'. Le percentuali per cellula di `marker_test_disc_vs_conc`, `doublet_tcr_pos_vs_neg` e `ambient_discriminant` sono descrittive, senza test.

**Test (`tests/`, pytest).**
- `test_scpipe.py`: formule confrontate con calcoli a mano (Nadeau-Bengio, Jaccard, LogCPM, HVGSelector).
- `test_leakage.py`: per 05, 11, 12, 13 e 20, i parametri stimati in un fold non cambiano quando si alterano le cellule di test; gli split di 04 non mettono lo stesso paziente in due fold.
- `test_stats.py` (bande dichiarate nel docstring): copertura del cluster bootstrap in [0,90; 0,99] con 10 e 11 pazienti; falsi positivi del block permutation test in [0,02; 0,08]; falsi positivi del test di permutazione di Spearman con n = 10 (H1-H4, medie per paziente di 17) in [0,02; 0,08]; p non definito (NaN) con un esito costante; Nadeau-Bengio sotto H0 ≤ 0,08 con k = 3 e k = 5 e potenza ≥ 0,5; Holm monotono e ≤ 1; stimatore della discordanza senza reinserimento non distorto (enumerazione esatta).
- `test_cnv.py`: filtro del rumore globale di inferCNV (indipendente dall'ordine delle cellule), finestre identiche a quelle di infercnvpy, attribuzione delle finestre ai loci e controllo dei geni marcatori (CDKN2A, TP53, SMAD4, MYC, KRAS).
- `test_tcr_ties.py`: unita' in pareggio escluse dall'errore principale, sensibilita' "pareggi = errore" e "pareggi = corretti", risultato indipendente dall'ordine delle cellule.
- `test_smoke.py`: ogni script si importa e risponde a `--help`.

L'intervallo percentile semplice del cluster bootstrap usato nella versione precedente copriva l'88,8-92,3% con 10-11 pazienti ed e' stato sostituito dall'intervallo "expanded percentile" (Hesterberg 2015), che copre il 93,5-94,8% (`sim_a_cluster_bootstrap.csv`). La copertura e' verificata solo per la media pooled di un esito per clone (o per coppia clone-compartimenti); per gli IC di Spearman (16, 17, 19), di d1/d2 (10), dello split-half di 06 (differenza di due medie su insiemi di pazienti diversi) e con 4-5 pazienti non e' verificata o non e' calibrata. Con 4-5 pazienti nemmeno l'intervallo espanso raggiunge il 95% (vedi `sim_a_cluster_bootstrap.csv`).

---

## 7. Limiti

- **Pochi pazienti.** Le conclusioni fra pazienti poggiano su 5 (task1a_paired, Adjacent_normal) - 14 pazienti; con 4 pazienti (Adjacent_normal nello schema solo sangue) l'IC del cluster bootstrap non e' calibrato e le ricampionature distinte sono 35.
- **Etichette di riferimento.** Le etichette `all_celltype` sono quelle degli autori, ottenute dalla stessa matrice: per `task2` le prestazioni misurano la riproducibilita' dell'annotazione, non la correttezza biologica. L'etichetta di malignita' e' derivata da inferCNV sugli stessi dati, con riferimento comune a tutti i pazienti e soglia stimata sulle stesse cellule di riferimento.
- **Identita' clonale dai marcatori.** Il lignaggio CD4/CD8 del clone e' assunto stabile, con le eccezioni note (CD4 citotossiche, CD4+CD8αα+, MAIT, NKT, γδ); CD4/CD8 da pochi UMI per cellula; i risultati principali usano CD8A > 0 o CD8B > 0, con regole piu' severe in `error_rate_identity_rules.csv`. La soglia ambientale di CD4 nel PBMC poggia su 11 cellule.
- **Doppietti.** d1 e d2 non sono indipendenti e dipendono da mediane pooled fra tessuti; nel Tumor la maggioranza di doppietti non e' sostenuta.
- **CNV.** Ogni finestra di inferCNV (100 geni consecutivi) ha le coordinate reali dei suoi geni (`cnv_windows.csv`), ma copre molti Mb (le 10 con CDKN2A 15-28 Mb): il segnale di un locus piccolo e' quello della regione circostante, e alcune finestre di 12p e 20q contengono anche geni dell'altro braccio; i confini dei loci sono quindi approssimati. Entrano solo i cromosomi 1-22: nella versione precedente entravano anche il genoma mitocondriale (la cui finestra segue la frazione mitocondriale delle cellule, una metrica di qualita') e un contig non standard; i numeri della riga 3 sono ricalcolati senza.
- **Ipotesi per paziente.** n = 10, potenza bassa; H1 ha una circolarita' (punteggio calcolato sulle cellule etichettate CD8T) e una direzione ipotizzata opposta agli errori osservati.
- **Coorte snRNA-seq** esclusa: pazienti diversi dalla coorte scRNA-seq.
- **Una sola coorte**: nessuna validazione esterna.

### Limiti noti del codice

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

### Non risolti

Problemi segnalati dai revisori dell'ultimo giro (giro 8) e non corretti, con il motivo. Nessuno cambia un numero principale sul dataset reale.

| Dove | Problema | Motivo |
| :--- | :--- | :--- |
| `src/06_clonotype_consistency.py` (riepilogo split-half), `src/16_pseudobulk_patient.py`, `src/17_clonal_expansion.py`, `src/18_loyo_check.py`, `src/19_ni_score.py` | i pareggi dell'etichetta di maggioranza sono esclusi, ma la sensibilita' "pareggi = errore / pareggi = corretti" e' tabulata solo in 09 e 10 (`arbitration_tie_sensitivity`, `error_rate_tie_sensitivity`) | il revisore l'ha ricalcolata e le conclusioni non cambiano (17, errore: ρ 0,186 e 0,217 contro 0,220; H1-H4 con segno e p simili; 06: +0,052 e +0,058 contro +0,063): `revisione/log/giro8_revisore_b.md`, punto 6. Questi valori non sono in `results/tables` e non sono citati come risultati |
| `README.md:94` | la radice indicata e' `progetto-ml/`, la cartella locale si chiama `pdac-ml` | `progetto-ml` e' il nome del repository remoto, richiesto dalla correzione D3 |

---

## Autore

Emanuele Colasanto, matricola 859776 — Corso di Machine Learning, A.A. 2025/2026.

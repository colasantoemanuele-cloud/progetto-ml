# Giro 7 (chiusura 4) - area (d) README, ambiente, biologia

**1. `README.md:32`, `README.md:35`, `README.md:259`; `src/09_arbitrate.py:7-8`**
- gravità: MEDIA (e)
- criterio violato: 5
- motivo: tutta la validazione clonale assume che il lignaggio CD4/CD8 di un clone non cambi. Lo dicono `09_arbitrate.py:7-8` ("Il lignaggio CD4/CD8 di un clone T maturo non cambia fra tessuti") e `README.md:32` ("permette di misurare l'affidabilità dell'annotazione"). Le eccezioni previste dal controllo obbligatorio non sono dichiarate da nessuna parte. Non compaiono CD8αα (CD8A+ CD8B−, incluse le CD4 di tessuto che acquisiscono CD8αα), MAIT/NKT e γδ: `grep -i "mait|nkt|αα"` su README, src e RISULTATI non trova nulla. Le CD4 citotossiche compaiono solo come alternativa in `README.md:35`. Il §7 "Identità clonale dai marcatori" (riga 259) parla solo della scarsità di UMI. La mia verifica sui 30 file VDJ dice che sul dataset reale i numeri non cambiano (script `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_d/mait.py` e `.../closing4_d/nk_mait.py`):
  - fra i 56 cloni in errore nel Tumor (`arbitration_blood`) ci sono 0 cloni MAIT (TRAV1-2 con TRAJ33, TRAJ12 o TRAJ20) e 0 iNKT (TRAV10 con TRAJ18);
  - i 12 cloni MAIT-like stanno tutti fra i 231 corretti, con identità CD8T;
  - fra le 1.673 NK-TCR+ le MAIT-like sono il 2,6% e le iNKT lo 0,06% (nelle CD8T il 3,8%);
  - CD8A è rilevato nel 7% delle cellule dei cloni discordanti, contro l'ipotesi CD4+CD8αα.
- correzione proposta: in `09:7` e in README §1 (riga 32) e §7 (riga 259) scrivere "il lignaggio è stabile, con eccezioni note (CD4 citotossiche, CD8αα, MAIT/NKT)". Se si vuole riportare il controllo con i numeri, salvarlo prima in una tabella di `results/tables` (criterio 4), per esempio da 10. Altrimenti citarlo senza numeri.

**2. `README.md:22`**
- gravità: MEDIA (e)
- criterio violato: 5
- motivo: la frase "Per questo anche il 41-58% delle duttali adiacenti ... risulta 'maligno'" presenta come certa una spiegazione mai verificata: il riferimento comune e la soglia stimata sulle stesse cellule immunitarie causerebbero il 41-58%. Nessuna analisi la mette alla prova; per esempio 13 non calcola la frazione oltre soglia in un tipo non immunitario e non epiteliale. La riga stessa elenca subito dopo un'altra spiegazione, il campo di cancerizzazione, e ne esiste una terza: cellule tumorali nel campione adiacente. `13_infercnv.py:19-23` è più cauto ("tende a superarla"). I valori sono corretti: PA06 0,577, PA08 0,410, PA10 0,530 in `malignancy_labels`.
- correzione proposta: "Questo può contribuire al fatto che il 41-58% ... risulti 'maligno'; i dati non distinguono fra artefatto della soglia, cellule tumorali nel campione adiacente ed effetto di campo."

**3. `README.md:273` (tabella "Limiti noti del codice"); `revisione/RISULTATI_FINALI.md:434` (viene da `revisione/limiti_noti.md`)**
- gravità: BASSA
- criterio violato: 4
- motivo: la copertura "0,928-0,943" con le numerosità reali delle 1.114 coppie esiste solo nel rapporto `revisione/log/chiusura3_revisore_c.md:20`. Non è in nessuna tabella di `results/tables`, né prodotta da codice della repository.
- correzione proposta: aggiungere a `sim_a_calibration.py` lo scenario con le numerosità reali (62, 24, 184, 75, 282, 118, 146, 57, 97, 56, 13) e una riga in `sim_a_cluster_bootstrap.csv`. In alternativa scrivere "verifica esterna di un revisore, non riprodotta dal codice" senza il numero.

**4. `README.md:35`; `revisione/RISULTATI_FINALI.md:191, 194, 197` (righe Adjacent_normal della tabella degli schemi)**
- gravità: BASSA
- criterio violato: nessuno (il §7, riga 257, e RISULTATI riga 47 lo dichiarano in generale)
- motivo: gli IC dell'errore in Adjacent_normal sono riportati accanto al valore senza dire, nel punto, che non sono calibrati: 0,059 [0,000; 0,093] con 4 pazienti; 0,061 e 0,076 con 5 pazienti (`n_distinct_resamples` 35 e 122 in `error_rate_schemes`). Per la differenza appaiata lo si dice esplicitamente.
- correzione proposta: in README:35 aggiungere "(IC non calibrato)" dopo l'IC dell'Adjacent_normal. In `genera_risultati.py`, nella tabella degli schemi, marcare le righe con ≤ `SMALL_K` pazienti (asterisco o colonna "ricamp. distinte").

**Conteggio:** ALTA 0 / MEDIA 2 / BASSA 2 / LIMITI NOTI nuovi 0.

**Controlli biologici, uno per uno:**
1. **Identità CD4/CD8 stabile nel clone con le eccezioni dichiarate: NON SUPERATO.** Le eccezioni non sono dichiarate (punto 1). Sui dati reali i numeri non cambiano: 0 MAIT e 0 iNKT fra i cloni in errore.
2. **Riferimento CNV: SUPERATO.** Immunitarie di tutti i pazienti, una sola esecuzione, soglia per paziente: README:22, 71, 258; `13:12-17`; `14:16-17`; codice `13:491` e `13:546-563`.
3. **Direzione dei loci: SUPERATO.** Perdite 9p21, 17p13, 18q21, 1p, 6q (direction −1); guadagni 8q24, 12p12, 20q (+1) in `14:98-106`. Le coordinate GRCh38 contengono CDKN2A, TP53, SMAD4, MYC e KRAS. Lo score è moltiplicato per la direzione (`14:200`).
4. **Geni di stress da dissociazione fra i geni comuni XAI: SUPERATO.** FOS, CD69, BTG2, HSPH1, DUSP4 e CREM sono indicati come geni di risposta immediata, CTRB1 come marcatore acinare (README:28; docstring di 20).
5. **Confronto tumore-adiacente su 4 pazienti: SUPERATO.** +0,080 [0,000; +0,192], 35 ricampionature distinte, "IC non calibrato", 99,7% delle repliche > 0 (README:35 e 257). Resta la nota BASSA sull'IC del solo Adjacent_normal (punto 4).
6. **Quota degli errori nel tumore dovuta alle cellule etichettate NK: SUPERATO.** README:35 dà 18 cloni su 287 (0,063: 12 CD8 e 6 CD4 nel sangue) e il tasso senza etichette NK, 0,141 [0,058; 0,261]. README:36 dice che nel Tumor la maggioranza di doppietti non è sostenuta (d1 = 0,510, d2 = 0,419), coerente con `doublet_mixture`. Le NK-TCR+ non sono spiegate da MAIT o iNKT (2,6% e 0,06%).

**Altre verifiche (nessun problema):**
- **Numeri:** ogni numero di README e RISULTATI_FINALI coincide con le tabelle, tenendo conto dell'arrotondamento, e la tabella "da dove viene ogni numero" è corretta.
- **RISULTATI_FINALI.md:** l'ho rigenerato in scratch con `genera_risultati.py --strict`: identico byte per byte al file versionato, 0 valori mancanti.
- **Comandi del README:** corrispondono a `--help` di tutti gli script e l'ordine è eseguibile.
- **Ambiente:**
  - tutti i pacchetti importati da `src/`, `tests/` e `revisione/` sono in `envs/environment.yml`;
  - il `content_hash` di `conda-lock.yml` (a0b5cf54…) coincide con quello ricalcolato da `environment.yml`;
  - le versioni nel lock, nel file e nell'ambiente installato coincidono;
  - il ruff di `.pre-commit-config.yaml` (v0.15.20) è quello dell'ambiente.
- **Docstring e argparse di `src/` e `src/sim/`:** complete e coerenti con il codice (file scritti, argomenti, default).
- **Correzioni A6, A7, D1, D2 e D3 della mia area:** applicate.
- **05, 11 e 12 non rieseguiti:** `git diff` (52ace44 e 8058701 → HEAD) mostra solo controlli sugli argomenti e sugli input, `help`, docstring e la riemissione di avvisi diversi dalla convergenza. L'AST di `LogCPM`, `HVGSelector`, `macro_f1`, `nadeau_bengio_ttest`, `jaccard` e `mean_pairwise_jaccard` è identico. Anche 00-03 e `download_tcr.sh`, dopo 61fb36b, cambiano solo in `help` e messaggi.

**pytest** (`-m "not slow" -p no:cacheprovider`): 103 passati, 0 falliti, 0 saltati, 6 deselezionati (i lenti), 90,6 s.

**ruff:** `ruff check src tests revisione` → "All checks passed!". In più `ruff format --check`: 35 file già formattati.

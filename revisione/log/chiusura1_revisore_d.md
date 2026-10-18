# Giro di chiusura 1 - area (d) README, ambiente, docstring

## Area (d): elenco dei problemi (giro di chiusura)

**1. `revisione/genera_risultati.py:539`**
- gravita': ALTA (è un crash; si verifica solo se manca una tabella, quindi si può declassare)
- criterio violato: 8
- motivo: `nf = pd.concat([zs.n_folds, z1.n_folds]) if not zs.empty else ...` controlla solo `zs`. Ho tolto le tabelle una alla volta, in una copia nella scratchpad con link simbolici, e ho rieseguito lo script. Se manca solo `model_zoo_summary_task1.csv` lo script termina con `AttributeError: 'DataFrame' object has no attribute 'n_folds'`, invece di segnalare la tabella in `MISSING` come promette la docstring (`--strict`). È l'unica tabella su cui succede: senza nessuna tabella lo script funziona, come senza una qualsiasi delle altre.
- correzione proposta: `nf = pd.concat([d.n_folds for d in (zs, z1) if not d.empty])`, e se il risultato è vuoto `zs_folds = "⟦n.d.⟧"`.

**2. `README.md:11` (e tabella righe 13-17)**
- gravita': MEDIA
- criterio violato: 5
- motivo: la frase "In tutti e tre la stima di pratica corrente e' piu' alta di quella che tiene conto della struttura dei dati, per tre ragioni ... riconducibili allo stesso principio: l'unita' ... e' il paziente" non regge per le righe 1 e 3.
  - Riga 1: la tabella indica come artefatto la circolarità. Le due stime (0,964 contro 0,954, divario 0,0093 in `baseline_gap_sc`) sono entrambe circolari, e lo split per paziente non corregge la circolarità.
  - Riga 3: la stima "di pratica corrente" (AUC 0,829, `malignancy_lopo_summary`) è già una LOPO, cioè per paziente. La colonna "struttura" riporta un'altra grandezza (AUC dei loci 0,682 pooled contro 0,540 entro paziente, `cnv_loci_auc`, `cnv_loci_auc_by_patient`). Le due colonne non confrontano la stessa stima.
- correzione proposta: limitare il principio del paziente alla riga 2 e al contrasto pooled/entro paziente della riga 3. Nella riga 3 mettere 0,682 (pooled) come "pratica corrente" e ≈0,54 come stima entro paziente, e spostare 0,829 nelle note. Nella riga 1 dire che il divario casuale/paziente è piccolo e che la circolarità resta in entrambe le stime.

**3. `README.md:165` (Note)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: "07 va eseguito con clone_b dopo clone_ab perche' 08 e 17 leggono pairwise_discordance_clone_b.csv" è una motivazione sbagliata.
  - Ogni esecuzione scrive file con il proprio suffisso `{cc}`.
  - Le tabelle comuni (`ambient_discriminant`, `chat_sensitivity_*`, `ambient_vs_real_tcr.png`) non dipendono da `--clone-def` (07, righe 361-428 e 568-600: usano `gex`/`has_tcr = clone_b.notna()`).
  - L'ordine fra le due esecuzioni è quindi indifferente: conta solo che clone_b giri prima di 08 e 17.
- correzione proposta: "07 con `--clone-def clone_b` va eseguito prima di 08 e 17, che leggono `pairwise_discordance_clone_b.csv`".

**4. `README.md:165` (tempi di 11) e `revisione/RISULTATI_FINALI.md:381-414`**
- gravita': BASSA
- criterio violato: 4 (solo in senso letterale: sono metadati di esecuzione, non risultati)
- motivo: questi numeri non sono in nessuna tabella di `results/tables/`, ma solo nei log.
  - README: "circa 40 minuti ... 16 ... 18 ... 9". Ho verificato che sono corretti: le somme dei tempi in `revisione/log/11_final/*.log` danno 2.458 s, 1.085 s e 528 s.
  - RISULTATI_FINALI.md: i picchi di memoria e le durate copiati da `note_esecuzione.md`.
- correzione proposta: dichiarare esplicitamente che vengono dai log (`revisione/log/11_final/`, `_status.txt`), oppure aggiungere una riga nel testo del criterio 4.

**5. `src/15_inventory.py:23`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice "--tcrdir ... (si legge l'intestazione del primo)", ma alle righe 179-183 il codice apre e conta per intero ogni file VDJ per `n_contigs` di `vdj_inventory.csv`.
- correzione proposta: "si legge l'intestazione del primo file e si contano i contig di tutti i file".

**6. `revisione/genera_risultati.py:9-11` (e righe 193, 311, 582-585, 630, 643)**
- gravita': BASSA
- criterio violato: nessuno (i valori oggi coincidono con le tabelle)
- motivo:
  - La docstring dice che "numeri di pazienti, fold e cellule di training sono letti dalle tabelle", ma alcuni sono scritti a mano: "in tutti i 5 fold LOPO" (630), "/50" e `k == 50` (193, 311), "top-50" (643, mentre la tabella ha la colonna `top_k`), "5 fold stratificati" (583, ma `sim_a_nadeau_bengio` ha k = 3 e 5).
  - La frase "Con 5-13 fold la potenza e' bassa" (582) è smentita per 13 fold dallo stesso output (task2_hier3: 12 modelli su 12 con p < 0,05).
- correzione proposta: leggere k, `top_k` e il numero di fold dalle tabelle; scrivere "k = 3 e 5"; limitare la frase sulla potenza a task1a_paired (5 fold).

**7. `src/05_baseline.py:41`, `06_clonotype_consistency.py:59`, `11_model_zoo.py:65,67`, `16_pseudobulk_patient.py:52`, `17_clonal_expansion.py:48-49`, `19_ni_score.py:40-41`, `20_xai_per_patient.py:52`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: nella sezione "Uso" di queste docstring (stringhe non raw) c'è un `\` a fine riga. In `__doc__` diventa una continuazione di riga, quindi i comandi si uniscono con spazi multipli. Si vede con `ast.get_docstring`, per esempio "...splits_sc.csv         --tasks ...". In 12 è scritto correttamente `\\`.
- correzione proposta: usare `\\` come in 12, oppure `r"""`.

**8. `src/scpipe.py:12-44`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: la sezione "Contenuto" della docstring non elenca:
  - `majority_label`, usata da 06, 09 e 10 e citata nelle loro docstring e nel README (riga 77);
  - `percentile_ci` ed `expanded_alpha`, usate da 06;
  - `cluster_bootstrap_replicates` e `SMALL_K`.
- correzione proposta: aggiungerle a "Contenuto".

**9. `src/02_build_h5ad.py:11`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: "il picco di memoria resta intorno a 10-12 GB", ma con il comando del README (`--no-counts-layer`) il picco misurato è 7,7 GB (`revisione/log/_status.txt`: `02_build_h5ad_sc peakKB=7655848`).
- correzione proposta: "circa 8 GB con --no-counts-layer, di più con il layer counts".

**10. argparse di quasi tutti gli script (per esempio `12_feature_selection.py:246`, `13_infercnv.py:314`, `07_robustness.py:343`)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: quasi nessun argomento ha `help` o mostra il default in `--help`. Per esempio `--scheme` di 12 (default `fold_group_5`) decide su quali fold si calcola lo Jaccard, ma da `--help` non si vede.
- correzione proposta: `formatter_class=argparse.ArgumentDefaultsHelpFormatter` più un `help` breve per ogni argomento.

**11. `README.md:169-192` (Struttura dei risultati)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: la tabella elenca solo le tabelle prodotte; non dice quale script produce ciascuna figura di `results/figures/` (37 file). La correzione D3 chiedeva "quale script produce cosa".
- correzione proposta: aggiungere una colonna "Figure" (per esempio 03: `qc_by_library_*`, `qc_complexity_*`, ...; 05: `gap_sc_*`; ...).

**12. `revisione/RISULTATI_FINALI.md:3,7,15` (generato da `genera_risultati.py:266-271` e `462-463`)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: l'intestazione dice che "Prima" è il valore della presentazione, ma la riga dell'MLP (0,903) viene dal "README pubblicato".
- correzione proposta: intestazione "Prima (presentazione o README pubblicato)", oppure una nota.

**13. `README.md:156-162`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: `ruff check src tests` esclude `revisione/` (che contiene `genera_risultati.py`, comunque pulito). Il README non indica come rigenerare `RISULTATI_FINALI.md`: va lanciato `python revisione/genera_risultati.py --strict` dalla radice, perché usa percorsi relativi.
- correzione proposta: `ruff check src tests revisione`, più il comando di `genera_risultati.py` dopo le simulazioni.

**14. `README.md:36`**
- gravita': BASSA
- criterio violato: nessuno
- motivo: "Fra le medie per paziente il segno e' lo stesso ma l'associazione non e' significativa (ρ = +0,467, p = 0,18, n = 10)" vale solo per l'errore binario. Per l'eccesso il livello paziente dà ρ = +0,321, p = 0,368 (`clonal_expansion_tests`, riga excess/patient_means), e il dato non è riportato.
- correzione proposta: aggiungere "(eccesso: ρ = +0,321, p = 0,37)".

## Verifiche senza problemi

- Ogni numero del README e di `RISULTATI_FINALI.md` coincide con `results/tables/` entro l'arrotondamento, inclusa la tabella "da dove viene ogni numero".
- `genera_risultati.py --strict`, rieseguito in una copia nella scratchpad, ricrea `RISULTATI_FINALI.md` identico byte per byte, con 0 valori mancanti.
- I comandi della sezione 4 corrispondono a `--help` di ogni script e l'ordine è eseguibile.
- Ogni pacchetto importato in `src/`, `tests/` e `revisione/` è in `envs/environment.yml`.
- Il `content_hash` di `conda-lock.yml` coincide con quello ricalcolato da `environment.yml` (`a0b5cf54…`); le versioni del lock coincidono con quelle dell'ambiente.
- `ruff format --check` passa, quindi l'hook ruff-format di `.pre-commit-config.yaml` è soddisfatto.
- Le correzioni D1, D2 (docstring) e D3 della mia area sono applicate.
- I riepiloghi 05, 11 e 12 coincidono con le tabelle per fold.
- `git diff` conferma: 05 invariato da 52ace44; in 12 cambia solo il controllo sui `cell_id`; in 11, da 8058701, cambiano solo docstring e `--from-table`; le funzioni di `scpipe.py` usate da 05, 11 e 12 sono invariate.

## Conteggio

- ALTA: 1 / MEDIA: 1 / BASSA: 12
- pytest (`-m "not slow"`): 102 passati, 0 falliti, 0 saltati (6 lenti deselezionati)
- ruff (`check src tests revisione`): All checks passed!

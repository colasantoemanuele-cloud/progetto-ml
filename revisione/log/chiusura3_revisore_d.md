# Giro di chiusura 3 (giro 6) - area (d) README, ambiente, docstring

## Area (d), giro di chiusura 3: problemi trovati

**1. `revisione/genera_risultati.py:162-166`**
- gravita': ALTA (è un crash, ma succede solo se una tabella è malformata: si può declassare a MEDIA)
- criterio violato: 8
- motivo: il calcolo di `auc_pat` (`loci_p.loc[loci_p.patient != "_mean", "auc_PDAC_COMPOSITE"]`) sta prima del blocco `try` (riga 175), quindi non è protetto come il resto della costruzione.
  - Prova: ho ripetuto il generatore 409 volte in copie nella scratchpad fatte di link simbolici. Ogni volta ho tolto una delle 50 tabelle lette, oppure una delle loro colonne.
  - Solo 2 casi vanno in crash, anche senza `--strict`, invece di finire in `MISSING`:
    - `cnv_loci_auc_by_patient` senza `patient`: `AttributeError: 'DataFrame' object has no attribute 'patient'`;
    - `cnv_loci_auc_by_patient` senza `auc_PDAC_COMPOSITE`: `KeyError: 'auc_PDAC_COMPOSITE'`.
  - Le altre 407 prove, e le prove con ogni tabella ridotta alla sola intestazione, terminano senza traceback.
  - Il REGISTRO (Giro 5, ALTA d1) dice invece "ora tutta la costruzione e' protetta".
- correzione proposta: spostare il calcolo di `auc_pat` dentro il `try`, oppure proteggerlo con `if not loci_p.empty and {"patient", "auc_PDAC_COMPOSITE"} <= set(loci_p.columns)`.

**2. `src/13_infercnv.py:79-81` (docstring, Output)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: la docstring dice "se nessun fold LOPO e' valutabile [...] i file corrispondenti non vengono scritti e quelli di un'esecuzione precedente vengono cancellati". Il codice però scrive sempre `malignancy_lopo.csv` (riga 588, `res.to_csv(...)` senza condizione, anche se `res` è vuoto). Cancella soltanto `malignancy_lopo_summary.csv` (righe 589-602).
- correzione proposta: "...non vengono scritti `malignancy_lopo_summary.csv` (cancellato se presente) ne' i file `external_test_*`; `malignancy_lopo.csv` resta vuoto".

**3. `README.md:258` contro la tabella "Da dove viene ogni numero" (righe 196-230)**
- gravita': BASSA
- criterio violato: nessuno (il numero è in una tabella e coincide)
- motivo: "La soglia ambientale di CD4 nel PBMC poggia su 11 cellule" non ha una riga nella tabella delle fonti. Il valore è in `ambient_cd4_threshold.csv` (PBMC, `n_cells` = 11). La tabella dice di elencare la fonte di ogni numero del README.
- correzione proposta: aggiungere la riga "11 cellule (soglia CD4 nel PBMC) | `ambient_cd4_threshold` (`n_cells`, riga PBMC)".

**4. `src/11_model_zoo.py:414` e `:424`; `src/02_build_h5ad.py:223` (testo di `--help`)**
- gravita': BASSA
- criterio violato: nessuno
- motivo: con `ArgumentDefaultsHelpFormatter`, il `--help` di questi argomenti mostra due default diversi:
  - `--tasks`: "default task2_hier3 task2_v1; [...] (default: None)";
  - `--scheme`: "(default: quello primario del task) (default: None)";
  - `--which` di 02 è obbligatorio ma compare come "(default: None)".
- correzione proposta: per questi argomenti usare `default=argparse.SUPPRESS` nel testo di aiuto, oppure un `help` che finisca con `%(default)s` solo dove il default è reale.

## Verifiche senza problemi
- **Numeri:** ho confrontato con `results/tables/` ogni numero del README (sezioni 1-7 e tabella "Da dove viene ogni numero") e di `RISULTATI_FINALI.md`, compresi i casi limite di arrotondamento (0,1375; 0,00185; 0,9225; 0,0265). Coincidono tutti.
- **Generatore:** `genera_risultati.py --strict`, rieseguito in una copia nella scratchpad, ricrea `RISULTATI_FINALI.md` identico byte per byte, con 0 valori mancanti.
- **Comandi della sezione 4:**
  - corrispondono al `--help` di ogni script, e tutti gli script rispondono, compreso `download_tcr.sh`;
  - l'ordine è eseguibile (sim_d prima di 20, 07 clone_b prima di 08 e 17, 10 prima di 16, 17, 19 e 21, 16 prima di 19);
  - ogni tabella e figura di `results/` (94 tabelle, 38 figure) è attribuita a uno script nella sezione 5.
- **Docstring di tutti gli script di `src/` e `src/sim/`:** scopo, input, output (confrontati con ogni `to_csv`, `savefig`, `savez` e `write_text`), argomenti e uso sono completi e coerenti con il codice.
- **Valori citati nelle docstring:** controllati sulle tabelle i numeri e le soglie citati:
  - Endocrine in 12 pazienti, di cui 2 con almeno 10 cellule;
  - 11 cellule per la soglia CD4;
  - 5-59 cloni per paziente, 287 in totale;
  - 12 pazienti con almeno 150 duttali;
  - PA04 senza PBMC;
  - 6 marcatori in 07;
  - 200 alberi e C = 0,05 in 12.
- **Ambiente:**
  - ogni pacchetto di terze parti importato in `src/`, `tests/` e `revisione/` è in `envs/environment.yml`;
  - `envs/` e `.pre-commit-config.yaml` sono invariati da 8058701;
  - le versioni chiave del lock (python 3.11.15, numpy, scipy, pandas, scanpy, anndata, scikit-learn, lightgbm, infercnvpy, ruff 0.15.20, pytest) coincidono con l'ambiente installato;
  - `ruff format --check` passa (35 file).
- **Correzioni dei rapporti precedenti:** i 9 problemi della chiusura 2 (d) e le correzioni D1-D3 della mia area sono applicati.
- **`git diff` dalla chiusura 1:**
  - 11 cambia solo nella nota su Wilcoxon, nei `choices` e nella logica di `--update`;
  - 12 cambia solo nella riemissione degli avvisi non di convergenza e nei `choices`;
  - 05 e `scpipe.py` sono invariati;
  - 04 da 52ace44 cambia solo nei controlli.
  - È confermato che 05, 11 e 12 non andavano rieseguiti.

## Conteggio
- ALTA: 1 / MEDIA: 0 / BASSA: 3
- pytest (`-m "not slow" -p no:cacheprovider`): 102 passati, 0 falliti, 0 saltati (6 lenti deselezionati)
- ruff (`check src tests revisione`): All checks passed!

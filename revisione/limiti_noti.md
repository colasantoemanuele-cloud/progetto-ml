## Limiti noti

Problemi del codice che riguardano dati simulati, test sintetici o input anomali (NaN, input vuoti, matrici dense, argomenti non usati dai comandi del README). Sul dataset reale GSE278694 non si presentano; per ciascuno e' indicato come e' stato verificato. Non sono corretti nel loop di revisione (decisione dell'utente).

| Dove | Problema | Verifica sul dataset reale |
| :--- | :--- | :--- |
| `src/scpipe.py:72` (`LogCPM`), `:97` (`HVGSelector`) | richiedono una matrice sparsa CSR: con un `ndarray` denso vanno in errore | tutti gli script convertono prima in CSR (`src/05_baseline.py:309`, `src/11_model_zoo.py:477`, `src/12_feature_selection.py:289`, `src/20_xai_per_patient.py:228`, `src/13_infercnv.py:324`); i test di leakage passano matrici CSR |
| `src/sim/sim_a_calibration.py:85`, `tests/test_stats.py:60` (`SIZES`) | con k = 11 l'undicesima numerosita' (40) e' aggiunta a mano; le coppie reali dell'eccesso (13-282 per paziente) sono piu' sbilanciate | riguarda solo la simulazione di calibrazione; un revisore ha ripetuto la copertura con le numerosita' reali delle 1.114 coppie (0,928-0,943, dentro la banda): `revisione/log/chiusura3_revisore_c.md`, punto 3 |
| `tests/test_stats.py:240`, `:247` | il test di Nadeau-Bengio usa gli stessi semi di `sim_a` (55.000 + k, 51.005): e' un sottoinsieme della stessa simulazione, non una seconda verifica | dichiarato nel docstring del test (righe 37-38); nessun numero del dataset reale dipende dal test |
| `src/11_model_zoo.py:390` (e figure `model_zoo_{task}`), `src/12_feature_selection.py`, `src/20_xai_per_patient.py` | i file di confronto e le figure dipendono solo dal task: un'esecuzione con un altro `--tag` o schema li sovrascrive | i comandi del README eseguono ogni task una sola volta, con uno schema solo; `model_comparison_{task2_v1,task2_hier3,task1a_paired}` coincidono con le tabelle per fold (verificato dai revisori dei giri di chiusura) |
| argparse degli script TCR e CNV | `--help` mostra "(default: None)" per gli argomenti obbligatori | solo testo d'aiuto; i comandi del README passano tutti gli argomenti obbligatori |

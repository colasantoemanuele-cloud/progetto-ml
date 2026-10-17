# Registro della revisione

Memoria di lavoro della revisione del ramo `revisione`. Ogni fase e ogni giro
del loop aggiornano questo file con stato, cosa resta da fare e comandi in corso.

Ambiente: `conda` env `pdac` (`/home/nemo/miniforge3/envs/pdac`), Python 3.11,
scikit-learn 1.9.0, scanpy 1.11.5, anndata 0.12.19, lightgbm 4.6.0, infercnvpy 0.6.1,
pytest 9.1.1, ruff 0.15.20. Macchina: 12 core, 30 GB RAM.

## Stato

| Fase | Stato | Commit |
| :--- | :--- | :--- |
| Verifica iniziale | fatta: remote `progetto-ml`, `sc_raw.h5ad` e 30 file VDJ presenti, MatrixMarket grezzi presenti | - |
| 0a Criteri | fatta | 4ef60ec |
| 0b Test | in corso | - |
| 1 Correzioni | da fare | - |
| 2 Loop | da fare | - |
| 3 Esecuzione finale | da fare | - |

## Note iniziali

- La cartella locale si chiama `pdac-ml/`, il remote e la radice citata nel README
  devono essere `progetto-ml/` (D3).
- Il pre-commit esegue `ruff` e `ruff format` sui file in commit.
- Tabelle versionate prodotte da comandi diversi da quelli del README:
  `model_zoo_sc.csv` usa `fold_lopo` per task2 (README: nessun `--scheme`),
  `model_zoo_task1.csv` = `--tasks task1a_paired --tag task1`,
  `baseline_sc.csv` copre 6 task (README: 3 di default).
- `17_clonal_expansion.py` stampa rho/p senza salvarli in tabella (criterio 4).

## Comandi in corso

(nessuno)

## Problemi respinti

(nessuno)

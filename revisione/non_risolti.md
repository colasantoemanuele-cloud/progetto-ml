## Non risolti

Problemi segnalati dai revisori dell'ultimo giro (giro 8) e non corretti, con il motivo. Nessuno cambia un numero principale sul dataset reale.

| Dove | Problema | Motivo |
| :--- | :--- | :--- |
| `src/06_clonotype_consistency.py` (riepilogo split-half), `src/16_pseudobulk_patient.py`, `src/17_clonal_expansion.py`, `src/18_loyo_check.py`, `src/19_ni_score.py` | i pareggi dell'etichetta di maggioranza sono esclusi, ma la sensibilita' "pareggi = errore / pareggi = corretti" e' tabulata solo in 09 e 10 (`arbitration_tie_sensitivity`, `error_rate_tie_sensitivity`) | il revisore l'ha ricalcolata e le conclusioni non cambiano (17, errore: ρ 0,186 e 0,217 contro 0,220; H1-H4 con segno e p simili; 06: +0,052 e +0,058 contro +0,063): `revisione/log/giro8_revisore_b.md`, punto 6. Questi valori non sono in `results/tables` e non sono citati come risultati |
| `README.md:94` | la radice indicata e' `progetto-ml/`, la cartella locale si chiama `pdac-ml` | `progetto-ml` e' il nome del repository remoto, richiesto dalla correzione D3 |

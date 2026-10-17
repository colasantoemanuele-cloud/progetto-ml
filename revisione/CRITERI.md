## Criteri di accettazione (da copiare in `revisione/CRITERI.md`)

1. **Nessun leakage:** ogni trasformazione stimata dai dati (normalizzazione con parametri, selezione dei geni, SVD, scaling, soglie) è stimata solo sul training del fold.
2. **Unità corretta:** ogni intervallo di confidenza e ogni test usa il paziente come unità di ricampionamento o di permutazione, salvo dove il README dichiara esplicitamente un'altra unità e ne spiega il motivo.
3. **Calibrazione:** ogni stimatore citato come validato ha un test in `tests/` che lo verifica su dati simulati, con la banda dichiarata nel docstring.
4. **Numeri tracciabili:** ogni numero del README e di `RISULTATI_FINALI.md` è in una tabella di `results/tables` prodotta dal codice attuale; nessun numero esiste solo a video.
5. **Affermazioni proporzionate:** nessuna affermazione nel codice o nel README va oltre i dati («esatto», «indipendente», «esclude», «dimostra» solo se letteralmente vero).
6. **Commenti:** ogni script ha una docstring con scopo, input, output (tutti i file scritti), argomenti e uso; ogni passaggio non ovvio ha un commento; nessun commento contraddice il codice.
7. **Riproducibilità:** i comandi del README, eseguiti nell'ordine, ricreano esattamente `results/` partendo dai dati grezzi; `environment.yml` contiene tutti i pacchetti importati.
8. **Robustezza:** nessuno script va in crash su input vuoti o colonne categoriali; i casi limite producono un messaggio chiaro.
9. **Test e stile:** `pytest` passa tutto; `ruff check` non segnala nulla.

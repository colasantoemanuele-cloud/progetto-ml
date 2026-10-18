# Giro 7 (chiusura 4) - area (c) CNV, simulazioni, test

**Problemi trovati nell'area (c): 13, 14, src/sim, tests, pyproject**

1. `src/13_infercnv.py:259-261` (`cnv.tl.infercnv` chiamato senza `chunksize` e con `dynamic_threshold` di default), `src/14_cnv_loci.py:153-155`, `src/13_infercnv.py:284-297` (impronta)
- gravita': ALTA (a; tocca anche d nella docstring 13:12-17, "UNA sola esecuzione ... Solo la soglia e' calcolata per paziente")
- criterio violato: 5 e 6. I numeri dipendono dall'ordine delle cellule, che non fa parte del metodo, e la docstring descrive un trattamento uniforme che il codice non fa. Tocca anche il 7, perche' l'impronta non registra il parametro che conta.
- motivo: in infercnvpy 0.6.1 il filtro del rumore (`_infercnv_chunk`, passo 5) usa `noise_thres = 1.5 * np.std(x_res)`, calcolata **per blocco di 5.000 cellule** (`chunksize` di default). Le cellule entrano in ordine di libreria, quindi i blocchi sono confusi con tessuto e quota di duttali. Dalla cache `cnv_matrix.npz` (minimo |valore| non nullo per blocco) la soglia va da 0,0503 a 0,0590. E' piu' alta nei blocchi ricchi di duttali tumorali (blocco 5: 1.318 duttali, 0,0590; blocco 7: 1.663, 0,0574) e piu' bassa in quelli solo Adjacent_normal (0,0506-0,0525). Ho rieseguito inferCNV nello scratchpad con `dynamic_threshold=None` e la soglia applicata per blocco: riproduce la cache esattamente (differenza massima 0,0). Con una sola soglia globale (1,5·sd su tutto X_cnv = 0,0534), che non dipende dall'ordine:
  - AUC pooled dello score medio: 0,531 → **0,598**; entro paziente: 0,572 → **0,616**;
  - frazione maligna: PA06 Tumor 0,508 → 0,714, PA09 Tumor 0,613 → 0,690; duttali adiacenti maligne "41-58%" → **38-65%**;
  - maligne nel training LOPO: 55,3% → 62,4%; adiacenti predette maligne: **30,8% → 37,4%**;
  - AUC LOPO 0,837 → 0,836 (± 0,062), BA 0,736 → 0,746; composito pooled 0,698 → 0,707, entro paziente 0,557 → 0,564.

  Cambiano quindi i numeri della riga 3 e delle note del README e di RISULTATI_FINALI. La conclusione "entro paziente i loci non fanno meglio dello score medio" regge, anzi si rafforza. `infercnv_params()` non registra `chunksize`.
- correzione proposta: in 13 e 14 chiamare `infercnv(..., dynamic_threshold=None)` e applicare dopo una soglia globale, `X[np.abs(X) < 1.5*X.std()] = 0` (equivale a un solo blocco, senza la memoria densa di `chunksize=n_obs`). Registrare la scelta nell'impronta, correggere la docstring e rieseguire 13 e 14 (cache nuove), poi README e RISULTATI_FINALI.

2. `src/14_cnv_loci.py:164-187` (posizioni delle finestre interpolate in bp), docstring 14:18-21, `README.md:261`
- gravita': ALTA (a, e)
- criterio violato: 5 e 6 ("per qualche Mb" e' falso; un locus PDAC dichiarato non viene misurato)
- motivo: le finestre di infercnvpy sono equispaziate **per indice di gene**: la finestra i copre i geni ordinati [10i, 10i+99]. Il codice invece le interpola in modo lineare **in bp** fra il primo e l'ultimo gene. Ho ricostruito le posizioni vere dalla lista dei geni (stesso numero di finestre per cromosoma, verificato), con centro = gene 10i+50. Lo scarto fra posizione vera e interpolata ha mediana 5,9 Mb, 90° percentile 17,0 Mb, massimo 34,8 Mb. In particolare:
  - **CDKN2A_9p21**: le 2 finestre usate stanno in realta' a chr9:34,7-35,1 Mb (9p13); con le vere 9p21 hanno 0 finestre in comune. Il segnale "9p21" non misura 9p21.
  - Finestre vere contro usate: SMAD4_18q21 7 contro 4 (3 in comune); 1p 69 contro 47 (le 47 usate stanno davvero a 1,7-39,8 Mb); 20q 32 contro 26; TP53 14 contro 11.

  Con le posizioni vere: 1p 0,617 → 0,653, SMAD4 0,580 → 0,564; composito pooled 0,698 → 0,695, entro paziente 0,557 → 0,563. Con anche la correzione 1: 0,703 e 0,570. Cambiano `cnv_loci_auc.csv` (AUC e `n_windows` per locus); i numeri principali cambiano poco.
- correzione proposta: calcolare la posizione di ogni finestra dai geni ordinati per `start` di ciascun cromosoma (`start` del gene 10i+50, oppure media fra il gene 10i e il gene 10i+99) al posto di `np.linspace`. Aggiornare docstring e README (confini esatti a meno della larghezza della finestra), poi rieseguire 14.

3. `revisione/RISULTATI_FINALI.md:380` (da `revisione/genera_risultati.py:1092-1093`)
- gravita': BASSA
- criterio violato: nessuno
- motivo: e' il punto 9 della chiusura 3 e non e' stato corretto, anche se REGISTRO dice "BASSA: corrette tutte". "+0,143 con effetto paziente (+0,116; +0,181)" ha il formato di un IC, ma sono minimo e massimo delle 5 repliche di `sim_d_xai_stability`.
- correzione proposta: "(da +0,116 a +0,181 sulle 5 repliche)".

4. `src/sim/sim_a_calibration.py:26-27` e `tests/test_stats.py:12-13`
- gravita': BASSA
- criterio violato: nessuno
- motivo: "la copertura e' verificata con queste numerosita', non con quelle" segue subito la frase sulle coppie reali (13-282). Letto cosi', "queste" indica le coppie reali, cioe' l'opposto di quello che fa il codice (`SIZES`, 5-59 piu' il 40 aggiunto).
- correzione proposta: "la copertura e' verificata con `SIZES`, non con le numerosita' reali delle coppie".

5. `README.md:272` e `revisione/limiti_noti.md` (riga LogCPM/HVGSelector): `src/13_infercnv.py:324`
- gravita': LIMITE NOTO (descritto male)
- criterio violato: nessuno
- motivo: 13:324 e' `X = sp.csr_matrix(sub.X)` dentro `cnv_fingerprint`, non l'input della pipeline. Le conversioni in CSR prima di `LogCPM`/`HVGSelector` sono a 13:620 (`Xtr = sp.csr_matrix(train.X)`) e 13:652 (duttali adiacenti). Il limite non si presenta sul dataset reale, ma il riferimento con cui lo si verifica e' sbagliato.
- correzione proposta: sostituire `:324` con `:620, :652`.

**Altri controlli, senza problemi**
- **Numeri:** README (riga 3, note, sez. 6) e RISULTATI_FINALI (CNV, Simulazioni) coincidono con le tabelle attuali:
  - CNV: 8.370/334, 13/5 pazienti, 0,837 ± 0,071 su 12, 30,8%, 41-58%, 0,531/0,698, 0,557 (0,528-0,588), 0,572, 4.870/333;
  - simulazioni: (a) 2,7/1,8, 4,7/6,1, 99,3, 93,5-94,8, 88,8-92,3, 91,5-92,5, 82-83; (b) entro 0,005, 0,126/0,093, 0,385/0,092/0,442, entro 0,008; (c) tutte le percentuali, gli IC, 4/6/11/5 e 2; (d) +0,143/+0,030.

  I numeri CNV sono pero' da aggiornare dopo le correzioni 1 e 2.
- **Cache e riesecuzioni:** `mean|X_cnv|` di `cnv_matrix.npz` coincide con `cnv_scores.csv` (1,4e-16). 13, 14 e le simulazioni sono stati rieseguiti dopo l'ultima modifica degli script e dei moduli che importano (date dei file). Nessun file estraneo in `results/` per l'area.
- **Correzioni:** A1, A2, A7 (14), D2 13:24 e blocco C applicate. `--splits` di 13 tolto. Applicate anche le correzioni della chiusura 3 per l'area (c), punti 1, 2, 4, 6, 7, 8, 10, 11; il punto 3 e' passato fra i limiti noti; il punto 9 no (problema 3).
- **Test:** ogni stimatore citato come validato ha un test con la banda del docstring: cluster bootstrap con k = 10/11, block permutation, Spearman con n = 10, Nadeau-Bengio, discordanza e `d_within` (verifiche esatte). I test di leakage coprono 05, 11, 12, 13 (fold e modello finale), 20 e gli split di 04, ciascuno con controllo positivo.
- **Pacchetti:** tutti quelli importati nell'area sono in `envs/environment.yml`.

**Conteggio**
- ALTA: 2
- MEDIA: 0
- BASSA: 2
- LIMITI NOTI nuovi o descritti male: 1

**Controlli biologici**
- Riferimento CNV con le immunitarie di tutti i pazienti e soglia per paziente: OK nel codice (13:257-281, 555-567) e nel README (riga 22). Il filtro del rumore pero' e' stimato per blocco di 5.000 cellule (problema 1).
- Direzione dei loci PDAC: OK (perdite 9p21, 17p, 18q, 1p, 6q con `direction = -1`; guadagni 8q, 12p, 20q con `+1`; coordinate GRCh38 corrette). Le finestre assegnate pero' sono mal posizionate, e 9p21 di fatto misura 9p13 (problema 2).
- Confronto tumore-adiacente su 4 pazienti, senza intervallo presentato come affidabile: OK. README "≈ 0,56 (0,53-0,59 nei 4 pazienti)", nessun IC; "campo di cancerizzazione" presentato come "compatibile ... un'interpretazione, non una dimostrazione"; ESITO di 14 condizionale.
- Identita' CD4/CD8 nel clone con le eccezioni: fuori area (revisore d). Le CD4 citotossiche sono citate (README:35); una ricerca di CD8αα, MAIT e NKT in README e `src/` non trova nulla.
- Geni di stress da dissociazione fra i geni comuni XAI: OK nel README (riga 28: FOS, CD69, BTG2, HSPH1, DUSP4, CREM; CTRB1 acinare); area d.
- Quota degli errori nel tumore dovuta a cellule etichettate NK: il README riporta 18 cloni NK (0,063) e il tasso senza NK (0,141), e dice che nel Tumor la maggioranza di doppietti non e' sostenuta; da verificare dal revisore (d).

**pytest e ruff**
- `pytest -m "not slow" -p no:cacheprovider`: 103 passati, 0 falliti, 0 saltati (6 deselezionati), 96,0 s.
- `ruff check src tests revisione`: "All checks passed!"

Script di controllo (nessuna scrittura nella repo), in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing4_c/`:
- `chunks.py`: soglia per blocco dalla cache
- `rerun_cnv.py`: inferCNV senza soglia; replica esatta e soglia globale
- `globthr.py`, `lopo_global.py` (output in `lopo_global.out`): effetto su AUC, etichette, LOPO
- `winpos.py`, `truepos_auc.py`, `both.py`: posizioni vere delle finestre e AUC per locus

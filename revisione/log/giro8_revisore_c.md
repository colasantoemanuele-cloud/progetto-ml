# Giro 8 - area (c) CNV, simulazioni, test

**Problemi trovati, area (c): 13, 14, src/sim, tests, pyproject.toml**

1. `src/sim/sim_b_reference_schemes.py:161` (ricaduta in `README.md:240` e nella tabella (b) di `revisione/RISULTATI_FINALI.md:349-374`)
- gravità: MEDIA. Verifica a (pareggi esclusi nei consumatori di 10) e b. Criteri 4 e 5.
- criterio violato: 4, 5
- motivo: dal giro 7 `arbitrate` di 10 restituisce le unità in pareggio con `lab=None` ed `error=NaN`. In sim_b:
  - `est_error = g.error.mean()` salta i NaN, quindi esclude i pareggi;
  - `j["true_error"] = j.lab != j.true_identity` vale True per ogni pareggio (None != "CD8T"), quindi l'errore vero conta tutti i pareggi come errori.
  Stima e verità sono calcolate su insiemi di unità diversi, e la distorsione negativa pubblicata viene da questo.
  - Prova: ho rieseguito con gli stessi semi e 30 repliche (script in `/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-pdac-ml/04562894-cc7b-4d1c-962b-2f6af2123ba9/scratchpad/closing5_c/simb_check.py`). Riproduce esattamente la colonna `bias` della tabella. Con l'errore vero calcolato sulle sole unità senza pareggio:
    - H0: la distorsione passa da −0,012…−0,006 a 0,000…+0,005;
    - LOCO PBMC: errore vero 0,086 (non 0,096 / 0,095); distorsione +0,034 con artefatto moderato e +0,297 con artefatto forte;
    - schemi solo sangue e naive: entro 0,007 dall'errore vero.
  - La spiegazione del README «le stime escludono i pareggi …, che nella simulazione sono piu' spesso errori» è falsa: i pareggi sono contati come errori per costruzione, sempre.
  - Inoltre «fra −0,013 e −0,006» sotto H0 non coincide con la tabella: il minimo è −0,01205, cioè −0,012. Il −0,013 compare solo con artefatto.
  - Il docstring di sim_b non dice nulla sui pareggi.
- correzione proposta: in sim_b calcolare `true_error` solo dove `~j.tie` (oppure `est_error` con i pareggi come errore, ma lo stesso insieme per entrambi). Aggiungere `n_ties` all'output e documentare la regola nel docstring. Rieseguire sim_b, rigenerare RISULTATI_FINALI e correggere README:240 (scarti, «errore vero 0,096 / 0,095», frase sui pareggi).

2. `README.md:262` (sezione 7 «Limiti», voce CNV)
- gravità: MEDIA. Verifica c (finestre CNV, loci) e b. Criteri 5 e 6.
- criterio violato: 5, 6
- motivo: «Le coordinate delle finestre di inferCNV sono interpolate in modo uniforme per cromosoma» descrive il codice di prima del giro 7. Ora `window_table` (13:328-378) ricostruisce le coordinate reali dei geni di ogni finestra, e `README.md:21` lo dice («ogni finestra ha le coordinate reali dei suoi geni»). Il README si contraddice. È lo stesso problema segnalato in `chiusura4_revisore_c.md:17` per README:261: corretto nel codice ma non in questa riga.
- correzione proposta: sostituire con «ogni finestra copre 100 geni consecutivi (coordinate reali in `cnv_windows.csv`); le finestre sono larghe (fino a ~28 Mb su 9p), quindi un locus include geni vicini e i confini dei loci restano approssimati».

3. `src/14_cnv_loci.py:197-204` (`locus_windows`), loci `KRAS_12p12` e `20q_gain`
- gravità: BASSA
- criterio violato: nessuno (5 in senso lato)
- motivo: una finestra entra nel locus se si sovrappone anche di poco, quindi alcune finestre sono in gran parte sull'altro braccio cromosomico. Ricostruito dai geni dell'h5ad e da `cnv_windows.csv`:
  - in `KRAS_12p12` (15 finestre) fino al 72% dei geni di una finestra è su 12q (finestre che finiscono a 38-48 Mb, centromero ~35 Mb), in media il 20%;
  - in `20q_gain` (37 finestre) fino al 91% dei geni è su 20p, in media il 12%;
  - 9p21, 17p13, 8q24, 1p, 6q: 0%; 18q21: al massimo 13%.
  La direzione e le coordinate dei loci sono corrette (CDKN2A 21,97-22,00 Mb, TP53 7,66 Mb, SMAD4 51,03 Mb, MYC 127,74 Mb, KRAS 25,2 Mb; GRCh38). Il README:21 avverte solo per 9p.
- correzione proposta: attribuire le finestre per mediana degli start (`pos`) o per maggioranza dei geni nel locus, oppure dichiarare nel README la contaminazione fra bracci per 12p e 20q.

4. `src/13_infercnv.py:772-776` (messaggio stampato)
- gravità: BASSA
- criterio violato: 5
- motivo: il messaggio «La frazione predetta maligna misura quanto il classificatore separa le duttali tumorali dalle duttali peri-tumorali» è impreciso. Il classificatore è addestrato sull'etichetta CNV entro le sole duttali tumorali, non sul contrasto tumore/adiacente. Il 37,4% misura quante duttali adiacenti somigliano alle duttali tumorali con CNV sopra soglia.
- correzione proposta: «misura quante duttali adiacenti hanno il profilo trascrizionale delle duttali tumorali etichettate maligne da inferCNV».

5. `src/14_cnv_loci.py:198-199` (docstring di `locus_windows`)
- gravità: BASSA
- criterio violato: 6
- motivo: il docstring dice «(start del primo gene, end dell'ultimo)», ma `window_table` usa `end = w.end.max()` (end massimo dei geni della finestra; 13:371, come scritto anche nel docstring di 14:23).
- correzione proposta: «end massimo dei suoi geni».

6. `README.md:246-250` (elenco dei test nella sezione 6)
- gravità: BASSA
- criterio violato: 6
- motivo: mancano `tests/test_cnv.py` (filtro del rumore globale, finestre contro infercnvpy, controllo CDKN2A) e `tests/test_tcr_ties.py` (pareggi esclusi, sensibilità, indipendenza dall'ordine), aggiunti nel giro 7.
- correzione proposta: aggiungere le due voci.

7. `tests/test_scpipe.py:4`
- gravità: BASSA
- criterio violato: 6
- motivo: il docstring dice «confronti esatti (tolleranza numerica 1e-9)», ma il file usa anche `abs=1e-6` (righe 31-32), `rtol=1e-4` (riga 99) e `rtol=1e-6` (riga 107).
- correzione proposta: «tolleranze numeriche fra 1e-9 e 1e-4 (arrotondamento float32)».

**Conteggio: ALTA 0 / MEDIA 2 / BASSA 5**

**Esito delle verifiche**

- **a)** Nell'area è corretta, salvo il problema 1.
  - Soglia inferCNV: `global_noise_filter` riproduce lo step 5 di infercnvpy 0.6.1 con un solo blocco; ho controllato il sorgente: gli altri passi sono per cellula, e la sd include gli zeri impliciti. La soglia è indipendente dall'ordine, coperta dal test. 13 e 14 sono coerenti: `cnv_scores.csv` coincide con mean|X_cnv| di `cnv_matrix.npz` (scarto massimo 1,4e-16, 101.219 cellule).
  - Finestre: `window_table` coincide con `_running_mean_for_chromosome`, stessa condizione `n > window_size`. CDKN2A (chr9:21,97-22,00 Mb) cade in 10 finestre, tutte attribuite a `CDKN2A_9p21` (10 finestre del locus, larghe 14,5-28,4 Mb).
  - Pareggi: in sim_b il consumatore di 10 non è stato adeguato (problema 1).
- **b)** I numeri CNV del README (righe 17, 21, 22, 225-231) e di RISULTATI_FINALI coincidono con le tabelle: 8.370/334, 13/5 pazienti, 0,836 ± 0,062, 12 pazienti, 37,4%, 38-65%, 0,598, 0,717, ≈0,57 (0,566), 0,53-0,61, 0,616, 4.870/333, 10 finestre di 15-28 Mb. Anche le righe CNV e sim della tabella prima/dopo coincidono. Per le simulazioni (a), (c) e (d) del README §6 e di RISULTATI tutto coincide, compresi 93,8%/92,8% e le dimensioni `SIZES`/`PAIR_SIZES` contro `clonal_expansion_joined` e `pairwise_discordance_clone_b`. Per (b) i valori stampati sono coerenti con la tabella tranne −0,013 (in realtà −0,012), ma la tabella stessa è distorta dal problema 1. Le tabelle e le figure di 13, 14 e sim sono più recenti dell'ultima modifica degli script. Test di leakage: coprono 05, 11, 12, 13, 20 e gli split di 04, ciascuno con controllo positivo. Stimatori citati come validati, con banda nel docstring di test_stats: cluster bootstrap (k = 10, 11), block permutation (p e p centrato), Nadeau-Bengio, Spearman con n = 10.
- **c)** Punto per punto:
  - identità CD4/CD8 stabile nel clone: fuori area; sim_b la simula coerentemente con 10;
  - marcatori: sim_b e test_tcr_ties usano CD8A/CD8B/CD4 come 10; nessun problema;
  - riferimento CNV: immunitarie di tutti i pazienti, soglia per paziente; corretto in codice, docstring, output e README;
  - loci PDAC e direzione: corretti (perdite 9p21, 17p, 18q, 1p, 6q; guadagni 8q, 12p, 20q; coordinate GRCh38 corrette), con la contaminazione fra bracci per 12p e 20q (problema 3) e la riga del README non aggiornata (problema 2);
  - malignità e campo di cancerizzazione presentati come ipotesi: sì (13 docstring punto 3, 14 ESITO «compatibile», README:22); il messaggio di 13 è impreciso (problema 4);
  - errori nel tumore dovuti alle NK e doppietti: fuori area;
  - geni di stress fra i geni comuni XAI: fuori area (sim_d non ne parla);
  - tumore-adiacente su 4 pazienti: dichiarato in 14 e nel README, nessun intervallo presentato;
  - H1-H4: fuori area.

**pytest** (`-m "not slow"`): 108 passati, 0 falliti, 0 saltati, 6 deselezionati (lenti). **ruff check src tests revisione**: «All checks passed!».

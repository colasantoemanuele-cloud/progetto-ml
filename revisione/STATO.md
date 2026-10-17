# Stato della revisione

**Giro 1** (revisione completa), chiuso il 18/10 alle 01:45: 1 ALTA, 18 MEDIA, 29 BASSA; tutti corretti salvo 4 motivati in `REGISTRO.md` ("Problemi respinti"). Resta: giro 2 incrementale (un revisore sul `git diff`), poi il giro di chiusura con 4 revisori in parallelo e i test lenti.

**Test e stile.** `pytest -m "not slow"`: 100 passati; test lenti: 6/6 passati in Fase 1 (si ripetono al giro di chiusura). `ruff check src tests revisione`: pulito. Calibrazioni (`tests/test_stats.py`, `sim_a`): Nadeau-Bengio sotto H0 p < 0,05 nel 2,7% (k = 3) e 1,8% (k = 5), potenza 99,3%; block permutation sotto H0 5,0-5,9% (1.000 repliche per livello, IC sup. ≤ 7,5%); permutazione di Spearman con n = 10: 4,75% (banda [0,02; 0,08]).

**Numeri della presentazione cambiati** (tabella completa in `RISULTATI_FINALI.md`):

| Numero | Prima | Dopo |
| :--- | :--- | :--- |
| logistica L2 task2_hier3 | 0,851 | 0,850 (arrotondamento) |
| k casuali / L1 k = 200 / Jaccard L1 (12) | 0,361 / 0,952 / 0,439 | 0,367 / 0,955 / 0,480 |
| geni comuni HVG/chi2/ANOVA | 26/19/20 | 18/19/19 |
| IC eccesso +0,059 | [0,034; 0,100] | [0,031; 0,110] |
| IC errore Tumor 0,195 | [0,123; 0,274] | [0,110; 0,290] |
| IC Tumor − Adjacent +0,080 | [+0,021; +0,170] | [0,000; +0,192] |
| doppietti d2 | 0,790 [0,660; 0,876] | 0,676 [0,559; 0,781] (d1 0,723 [0,556; 0,848]) |
| AUC LOPO 13 / duttali adiacenti maligne | 0,825 / 31,7% | 0,829 / 37,1% |
| AUC entro paziente loci PDAC | 0,53 | 0,540 |
| sim (b) / (c) / (d) | 0,24 / 2 su 100 / 0,063-0,040 | 0,126-0,385 / 4-11 su 100, 4,8-5,9% su 1.000 / 0,143-0,030 |

**Bootstrap.** Tutti gli IC con il paziente come unita' (0,195; +0,059; +0,080; d1/d2; 06-10, 16, 17, 19) usano il cluster bootstrap con intervallo "expanded percentile" (Hesterberg 2015; `scpipe.percentile_ci`), 2.000 repliche (500 per d1/d2). Copertura al 95% nominale (`sim_a_cluster_bootstrap`, media pooled): 4 pazienti 82-83%, 5 pazienti 91,5-92,5%, 10 pazienti 93,5-93,8%, 11 pazienti 93,5-94,8%. Con 4-5 pazienti (Adjacent_normal, +0,080) l'IC non e' calibrato; per Spearman e d1/d2 la copertura non e' verificata.

**Simulazioni.** (a) trovata in `genomic-audit-mvp` e riscritta; (b), (c), (d) scritte da zero (esistevano solo come testo). Numeri sopra; dettagli in `RISULTATI_FINALI.md`.

**Esecuzione finale di 11**: completata al commit del giro 1 (unione di processi paralleli, identica a un'esecuzione unica; MLP reso deterministico).

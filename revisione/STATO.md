# Stato della revisione

**Giri:** 1 completo 1/18/29 (8058701); 2 incrementale 0/4/16 (50b8095); 3 incrementale 0/5/8 (61fb36b); Fase 3 eseguita; 4 = chiusura 1 2/9/38 (59098b6); 5 = chiusura 2 2/3/37 (04ebed5); **6 = chiusura 3 (03:53): 2 ALTA, 4 MEDIA, 22 BASSA** (ALTA + MEDIA: 11 al giro 4, 5 al giro 5, 6 al giro 6). Limite portato a 10 giri dall'utente; stop se ALTA + MEDIA non scende per due giri di fila o se un problema corretto ricompare. Resta: correzioni del giro 6 e giro 7 (chiusura 4).

**ALTA e MEDIA del giro 6** (rapporti in `revisione/log/chiusura3_revisore_*.md`):
- ALTA, `src/19_ni_score.py:245-246` + `src/scpipe.py:301-346`, criterio 8: con un NI score mancante (NaN) `np.ptp` non rileva il problema, rho = NaN e il p diventa il minimo (5e-05), "SIGNIFICATIVO" dopo Holm. Caso nuovo della funzione corretta al giro 2 (input costanti), non un conflitto.
- ALTA, `revisione/genera_risultati.py:162-166`, criterio 8: `auc_pat` calcolato prima del blocco protetto: crash se `cnv_loci_auc_by_patient` manca di `patient` o `auc_PDAC_COMPOSITE`. Terza segnalazione della stessa classe (crash del generatore), ogni volta in un punto diverso: protezione incompleta, non correzioni che si annullano.
- MEDIA, `src/20_xai_per_patient.py:166, 205-210, 257-264`, criteri 8 e 6: `--scheme` accetta schemi senza un paziente per fold (introdotto con `choices` al giro 5); con `fold_random_5` le etichette dei fold coincidono e le matrici di Jaccard vengono sovrascritte senza avviso.
- MEDIA, `README.md:16`, criterio 5: "il modello sfrutta la firma del paziente e la composizione cellulare" non e' misurato; il disegno e' incrociato (ogni paziente ha entrambi i tessuti) e lo split casuale condivide la libreria (paziente x tessuto).
- MEDIA, `src/14_cnv_loci.py:419-431`, criteri 6 e 5: il messaggio "meglio dello score medio entro paziente" confronta le AUC pooled; entro paziente lo score medio (0,572) supera i loci (0,557).
- MEDIA, `src/13_infercnv.py:253-257`, `src/14_cnv_loci.py:140-144`, criterio 8: nessuna cellula di riferimento -> errore poco chiaro di infercnvpy.

**Test e stile.** pytest completo: 108 passati (lenti compresi, giro 6). `ruff check src tests revisione`: pulito. Calibrazioni (`tests/test_stats.py`, `sim_a`): Nadeau-Bengio sotto H0 p < 0,05 nel 2,7% (k = 3) e 1,8% (k = 5), potenza 99,3%; block permutation sotto H0 5,0-5,9% (1.000 repliche per livello, IC sup. ≤ 7,5%); permutazione di Spearman con n = 10: 4,75% (banda [0,02; 0,08]).

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
| MLP task2_hier3 (README) | 0,903 | 0,904 (MLP reso deterministico) |
| AUC LOPO 13 / duttali adiacenti maligne | 0,825 / 31,7% | 0,837 / 30,8% (senza cromosoma MT) |
| AUC CNV: score medio / loci pooled / entro paziente | 0,492 / 0,682 / 0,53 | 0,531 / 0,698 / 0,557 |
| sim (b) / (c) / (d) | 0,24 / 2 su 100 / 0,063-0,040 | 0,126-0,385 / 4-11 su 100, 4,8-5,9% su 1.000 / 0,143-0,030 |

**Bootstrap.** Tutti gli IC con il paziente come unita' (0,195; +0,059; +0,080; d1/d2; 06-10, 16, 17, 19) usano il cluster bootstrap con intervallo "expanded percentile" (Hesterberg 2015; `scpipe.percentile_ci`), 2.000 repliche (500 per d1/d2). Copertura al 95% nominale (`sim_a_cluster_bootstrap`, media pooled): 4 pazienti 82-83%, 5 pazienti 91,5-92,5%, 10 pazienti 93,5-93,8%, 11 pazienti 93,5-94,8%. Con 4-5 pazienti (Adjacent_normal, +0,080) l'IC non e' calibrato; per Spearman e d1/d2 la copertura non e' verificata.

**Simulazioni.** (a) trovata in `genomic-audit-mvp` e riscritta; (b), (c), (d) scritte da zero (esistevano solo come testo). Numeri sopra; dettagli in `RISULTATI_FINALI.md`.


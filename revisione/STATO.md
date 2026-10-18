# Stato della revisione

**Giri:** 1 completo 1/18/29; 2-3 incrementali 0/4/16, 0/5/8; Fase 3; 4-6 chiusure 2/9/38, 2/3/37, 2/4/22; **7 (regole: solo dataset reale) 2/4/15**, corretti con le decisioni dell'utente (pareggi esclusi; soglia inferCNV globale; finestre CNV con coordinate reali). ALTA + MEDIA: 11, 5, 6, 6. Resta: giro 8 (ultimo), poi Fase 3 finale.

**Numeri cambiati al giro 7:** errore Tumor 0,195 -> 0,170 [0,095; 0,252]; Adjacent 0,059 -> 0,060; differenza +0,080 -> +0,059 [0,000; +0,160]; AUC LOPO 13 0,837 -> 0,836; adiacenti predette maligne 30,8% -> 37,4%; AUC loci pooled 0,698 -> 0,717, entro paziente 0,557 -> 0,566 (score medio 0,616).

**Test e stile.** pytest completo: 108 passati (lenti compresi, giro 6). `ruff check src tests revisione`: pulito. Calibrazioni (`tests/test_stats.py`, `sim_a`): Nadeau-Bengio sotto H0 p < 0,05 nel 2,7% (k = 3) e 1,8% (k = 5), potenza 99,3%; block permutation sotto H0 5,0-5,9% (1.000 repliche per livello, IC sup. ≤ 7,5%); permutazione di Spearman con n = 10: 4,75% (banda [0,02; 0,08]).

**Numeri della presentazione cambiati** (tabella completa in `RISULTATI_FINALI.md`):

| Numero | Prima | Dopo |
| :--- | :--- | :--- |
| logistica L2 / MLP task2_hier3 | 0,851 / 0,903 | 0,850 (arrotondamento) / 0,904 (MLP reso deterministico) |
| k casuali / L1 k = 200 / Jaccard L1 (12) | 0,361 / 0,952 / 0,439 | 0,367 / 0,955 / 0,480 |
| geni comuni HVG/chi2/ANOVA | 26/19/20 | 18/19/19 |
| eccesso +0,059 / errore Tumor 0,195 / Tumor − Adjacent +0,080 | [0,034; 0,100] / [0,123; 0,274] / [+0,021; +0,170] | +0,059 [0,031; 0,110] / 0,170 [0,095; 0,252] / +0,059 [0,000; +0,160] |
| doppietti d2 | 0,790 [0,660; 0,876] | 0,676 [0,559; 0,781] (d1 0,723 [0,556; 0,848]) |
| AUC LOPO 13 / duttali adiacenti maligne | 0,825 / 31,7% | 0,836 / 37,4% |
| AUC CNV: score medio / loci pooled / entro paziente | 0,492 / 0,682 / 0,53 | 0,598 / 0,717 / 0,566 |
| sim (b) / (c) / (d) | 0,24 / 2 su 100 / 0,063-0,040 | 0,126-0,385 / 4-11 su 100, 4,8-5,9% su 1.000 / 0,143-0,030 |

**Bootstrap.** Tutti gli IC con il paziente come unita' (0,195; +0,059; +0,080; d1/d2; 06-10, 16, 17, 19) usano il cluster bootstrap con intervallo "expanded percentile" (Hesterberg 2015; `scpipe.percentile_ci`), 2.000 repliche (500 per d1/d2). Copertura al 95% nominale (`sim_a_cluster_bootstrap`, media pooled): 4 pazienti 82-83%, 5 pazienti 91,5-92,5%, 10 pazienti 93,5-93,8%, 11 pazienti 93,5-94,8%. Con 4-5 pazienti (Adjacent_normal, +0,080) l'IC non e' calibrato; per Spearman e d1/d2 la copertura non e' verificata.

**Simulazioni.** (a) trovata in `genomic-audit-mvp` e riscritta; (b), (c), (d) scritte da zero (esistevano solo come testo); numeri sopra.

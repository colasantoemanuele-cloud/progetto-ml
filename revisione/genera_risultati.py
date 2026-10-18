#!/usr/bin/env python3
"""
genera_risultati.py - Scrive revisione/RISULTATI_FINALI.md leggendo solo results/tables/.

Scopo
-----
Ogni risultato di RISULTATI_FINALI.md viene da una tabella di `results/tables/`, citata
accanto al numero (criterio 4). Sono scritti a mano solo i valori della presentazione
(colonna "prima" della tabella iniziale e dizionario `h_old`) e la soglia SMALL_K = 5
pazienti sotto cui gli IC non sono calibrati; numeri di pazienti, fold e cellule di
training sono letti dalle tabelle. L'esito "arrotondamento" si applica solo alle righe
in cui il valore precedente e' la stessa stima arrotondata diversamente (indicate nel
codice passando il valore esatto) e il valore esatto dista meno di 0,001 da quello
precedente. Le note
sull'esecuzione (commit, memoria, script non rieseguiti) sono lette da
`revisione/note_esecuzione.md`, se esiste, e copiate in fondo, seguite dai limiti noti del
codice (`revisione/limiti_noti.md`).

Input
-----
  results/tables/*.csv, revisione/note_esecuzione.md e revisione/limiti_noti.md (facoltativi)

Output
------
  revisione/RISULTATI_FINALI.md

Argomenti
---------
  --strict  esce con errore se manca una tabella o una colonna

Uso
---
    python revisione/genera_risultati.py --strict
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TAB = Path("results/tables")
OUT = Path("revisione/RISULTATI_FINALI.md")
NOTE = Path("revisione/note_esecuzione.md")
LIMITS = Path("revisione/limiti_noti.md")
MISSING: list[str] = []
SMALL_K = (
    5  # pazienti: con <= SMALL_K l'IC del cluster bootstrap non e' calibrato (sim_a)
)


def t(name: str) -> pd.DataFrame:
    """Legge una tabella; se manca la registra e restituisce un DataFrame vuoto."""
    p = TAB / f"{name}.csv"
    if not p.exists() or p.stat().st_size <= 1:
        MISSING.append(name)
        return pd.DataFrame()
    return pd.read_csv(p)


def v(df: pd.DataFrame, col: str, **where) -> float:
    """Valore di `col` nell'unica riga che soddisfa le condizioni `where`."""
    if df.empty or col not in df.columns:
        MISSING.append(f"{col} {where}")
        return np.nan
    m = np.ones(len(df), dtype=bool)
    for k, x in where.items():
        m &= (df[k].astype(str) == str(x)).values
    if m.sum() != 1:
        MISSING.append(f"{col} {where} ({int(m.sum())} righe)")
        return np.nan
    return float(df.loc[m, col].iloc[0])


def f(x: float, nd: int = 3, sign: bool = False) -> str:
    """Numero in formato italiano (virgola decimale)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "⟦n.d.⟧"
    if x != 0 and abs(x) < 0.5 * 10**-nd:
        nd += 1  # evita "-0,000" per differenze piccolissime ma non nulle
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    return s.replace(".", ",")


def n_(x: float) -> str:
    """Intero con il punto come separatore delle migliaia."""
    return "⟦n.d.⟧" if np.isnan(x) else f"{x:,.0f}".replace(",", ".")


def pct(x: float, nd: int = 1) -> str:
    return f(100 * x, nd) + "%"


def ci(lo: float, hi: float, nd: int = 3, sign: bool = False) -> str:
    return f"[{f(lo, nd, sign)}; {f(hi, nd, sign)}]"


def p_fmt(p: float) -> str:
    if np.isnan(p):
        return "⟦n.d.⟧"
    if p < 0.001:
        return f"{p:.1e}".replace(".", ",")
    return f(p, 4) if p < 0.01 else f(p, 3)


def md_table(rows: list[list[str]], header: list[str]) -> str:
    out = [
        "| " + " | ".join(header) + " |",
        "|" + "|".join([" :--- "] * len(header)) + "|",
    ]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Scrive RISULTATI_FINALI.md dalle tabelle."
    )
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    # --- tabelle
    ds, dsn = t("dataset_summary_sc"), t("dataset_summary_sn")
    tasks = t("task_summary_sc")
    bs, gap = t("baseline_summary_sc"), t("baseline_gap_sc")
    zs, z1 = t("model_zoo_summary_sc"), t("model_zoo_summary_task1")
    cv1, ch3, c1 = (
        t("model_comparison_task2_v1"),
        t("model_comparison_task2_hier3"),
        t("model_comparison_task1a_paired"),
    )
    fs_v1, fs_t1 = (
        t("feature_selection_summary_task2_v1"),
        t("feature_selection_summary_task1a_paired"),
    )
    st_v1, st_t1 = t("feature_stability_task2_v1"), t("feature_stability_task1a_paired")
    core_t1 = t("feature_core_genes_task1a_paired")
    dg_v1, dg_t1 = (
        t("feature_selection_diag_task2_v1"),
        t("feature_selection_diag_task1a_paired"),
    )
    xs, xg = t("xai_summary_task1a_paired"), t("xai_common_genes_task1a_paired")
    xj = t("xai_jaccard_task1a_paired")
    exb, exab = t("excess_summary_clone_b"), t("excess_summary_clone_ab")
    sw = t("clone_switching_summary_clone_b")
    eb, pdiff = t("error_rate_blood"), t("error_rate_paired_diff")
    dec, rules = t("error_decomposition_blood"), t("error_rate_identity_rules")
    schemes = t("error_rate_schemes")
    mk, tri = t("marker_test_disc_vs_conc"), t("triple_clones_summary")
    tcr_rate, mix = t("tcr_rate_by_celltype"), t("doublet_mixture")
    ce, cps = t("clonal_expansion_tests"), t("clonal_expansion_perm_sensitivity")
    ly_e, ly_x = t("loyo_clonal_expansion_error"), t("loyo_clonal_expansion_excess")
    ly_s = t("loyo_clonal_expansion_perm_sensitivity")
    hyp, holm, ni = (
        t("pseudobulk_hypothesis_tests"),
        t("ni_score_holm_joint"),
        t("ni_score_test"),
    )
    mal, ext = t("malignancy_lopo_summary"), t("external_test_summary")
    loci, loci_p = t("cnv_loci_auc"), t("cnv_loci_auc_by_patient")
    if {"patient", "auc_PDAC_COMPOSITE"} <= set(loci_p.columns):
        auc_pat = loci_p.loc[loci_p.patient != "_mean", "auc_PDAC_COMPOSITE"]
    else:
        if not loci_p.empty:
            MISSING.append(
                "cnv_loci_auc_by_patient: colonne patient / auc_PDAC_COMPOSITE"
            )
        auc_pat = pd.Series(dtype=float)
    sa_nb, sa_cb = t("sim_a_nadeau_bengio"), t("sim_a_cluster_bootstrap")
    sb, sc_, sd = (
        t("sim_b_reference_schemes"),
        t("sim_c_block_permutation"),
        t("sim_d_xai_stability"),
    )

    L: list[str] = []
    try:

        def bsv(task, scheme, stat="mean"):
            return v(
                bs, stat, task=task, scheme=scheme, model="logreg_l2", metric="macro_f1"
            )

        def zv(df, task, model, col="macro_f1_mean"):
            return v(df, col, task=task, model=model)

        def fsv(df, method, k):
            return v(df, "macro_f1_mean", method=method, k=k)

        def stv(df, method, k):
            return v(
                df.rename(columns={df.columns[0]: "method"}) if not df.empty else df,
                str(k),
                method=method,
            )

        k_core = int(core_t1.k.min()) if not core_t1.empty else 50

        n_fold_t1 = f"{fs_t1.n_folds.max():.0f}" if not fs_t1.empty else "⟦n.d.⟧"
        n_reps_c = n_(sc_.reps.iloc[0]) if not sc_.empty else "⟦n.d.⟧"

        def core_n(method):
            if core_t1.empty:
                return np.nan
            return float(((core_t1.method == method) & (core_t1.k == k_core)).sum())

        if not hyp.empty:
            hyp = hyp.rename(columns={hyp.columns[0]: "hypothesis"})
        hcode = {
            "exhaustion_score": "H1",
            "myeloid_fraction": "H2",
            "caf_fraction": "H3",
            "H4_ni_score": "H4",
        }

        def h_rho_p(h):
            if h == "H4_ni_score":
                return v(ni, "rho", hypothesis=h), v(ni, "p_perm", hypothesis=h)
            return v(hyp, "rho", hypothesis=h), v(hyp, "p_perm", hypothesis=h)

        d1 = v(mix, "d1", statistic="median", stratum="all")
        d2 = v(mix, "d2", statistic="median", stratum="all")
        sdd = sd[sd.rep.astype(str) == "mean"] if not sd.empty else sd
        pe = (
            sd[(sd.condition == "patient_effect") & (sd.rep.astype(str) != "mean")]
            if not sd.empty
            else sd
        )
        sim_d_eff = float(pe.difference.mean()) if len(pe) else np.nan
        sim_d_null = v(sdd, "difference", condition="no_patient_effect")
        loco_pbmc = {
            c: (
                v(sb, "est_error", condition=c, scheme="loco", tissue="PBMC"),
                v(sb, "true_error", condition=c, scheme="loco", tissue="PBMC"),
            )
            for c in ("artifact_moderate", "artifact_strong")
        }

        # --- tabella prima / dopo (prima = presentazione)
        pres = [
            (
                "macro-F1 task2_v1, 5 fold per paziente (logistica L2)",
                "0,954",
                f(bsv("task2_v1", "fold_group_5")),
                "baseline_summary_sc",
            ),
            (
                "macro-F1 task1a_paired, split casuale",
                "0,875±0,002",
                f"{f(bsv('task1a_paired', 'fold_random_5'))}±{f(bsv('task1a_paired', 'fold_random_5', 'std'))}",
                "baseline_summary_sc",
            ),
            (
                "macro-F1 task1a_paired, LOPO",
                "0,729±0,100",
                f"{f(bsv('task1a_paired', 'fold_lopo'))}±{f(bsv('task1a_paired', 'fold_lopo', 'std'))}",
                "baseline_summary_sc",
            ),
            (
                "LightGBM task2_v1 (LOPO)",
                "0,962",
                f(zv(zs, "task2_v1", "lightgbm")),
                "model_zoo_summary_sc",
            ),
            (
                "logistica L2 task2_hier3",
                "0,851",
                f(zv(zs, "task2_hier3", "logreg_l2")),
                "model_zoo_summary_sc",
                zv(zs, "task2_hier3", "logreg_l2"),
            ),
            (
                "LightGBM task2_hier3",
                "0,923",
                f(zv(zs, "task2_hier3", "lightgbm")),
                "model_zoo_summary_sc",
            ),
            (
                "MLP task2_hier3 (README pubblicato)",
                "0,903",
                f(zv(zs, "task2_hier3", "mlp")),
                "model_zoo_summary_sc",
            ),
            (
                "RF vs LightGBM task2_v1: p Wilcoxon / p Nadeau-Bengio",
                "p 0,002/0,074",
                f"p {f(v(cv1, 'p_wilcoxon', model='random_forest'))}/{f(v(cv1, 'p_nadeau_bengio', model='random_forest'))}",
                "model_comparison_task2_v1",
            ),
            (
                "dummy task1a_paired: Δ dal migliore / p Nadeau-Bengio",
                "0,229/0,072",
                f"{f(-v(c1, 'delta_vs_best', model='dummy'))}/{f(v(c1, 'p_nadeau_bengio', model='dummy'))}",
                "model_comparison_task1a_paired",
            ),
            (
                "selezione feature, k casuali k = 50, task2_v1",
                "0,361",
                f(fsv(fs_v1, "random", 50)),
                "feature_selection_summary_task2_v1",
            ),
            (
                "selezione feature, L1 k = 50, task2_v1",
                "0,936",
                f(fsv(fs_v1, "l1_logreg", 50)),
                "feature_selection_summary_task2_v1",
            ),
            (
                "selezione feature, L1 k = 200, task2_v1",
                "0,952",
                f(fsv(fs_v1, "l1_logreg", 200)),
                "feature_selection_summary_task2_v1",
            ),
            (
                "Jaccard L1 k = 200, task2_v1",
                "0,439",
                f(stv(st_v1, "l1_logreg", 200)),
                "feature_stability_task2_v1",
            ),
            (
                f"geni comuni ai {n_fold_t1} fold (task1a, k = {k_core}): HVG / chi2 / ANOVA",
                "26/50, 19/50, 20/50",
                f"{core_n('hvg'):.0f}/{k_core}, {core_n('chi2'):.0f}/{k_core}, {core_n('anova'):.0f}/{k_core}",
                "feature_core_genes_task1a_paired",
            ),
            (
                "XAI: Jaccard LOPO",
                "0,300",
                f(v(xs, "jaccard_lopo")),
                "xai_summary_task1a_paired",
            ),
            (
                "XAI: Jaccard controllo",
                "0,459",
                f(v(xs, "jaccard_controllo")),
                "xai_summary_task1a_paired",
            ),
            (
                "sim (d): differenza con / senza effetto paziente",
                "0,063/0,040",
                f"{f(sim_d_eff)}/{f(sim_d_null)}",
                "sim_d_xai_stability",
            ),
            (
                "rumore di base d_within",
                "0,107",
                f(v(exb, "d_within")),
                "excess_summary_clone_b",
            ),
            (
                "eccesso fra compartimenti [IC 95%]",
                "+0,059 [0,034;0,100]",
                f"{f(v(exb, 'excess'), sign=True)} {ci(v(exb, 'lo_cluster'), v(exb, 'hi_cluster'))}",
                "excess_summary_clone_b",
            ),
            (
                "errore Tumor, schema solo sangue [IC]",
                "0,195 [0,123;0,274]",
                f"{f(v(eb, 'error', tissue='Tumor'))} {ci(v(eb, 'lo', tissue='Tumor'), v(eb, 'hi', tissue='Tumor'))}",
                "error_rate_blood",
            ),
            (
                "errore Adjacent_normal, schema solo sangue",
                "0,059",
                f(v(eb, "error", tissue="Adjacent_normal")),
                "error_rate_blood",
            ),
            (
                "Tumor − Adjacent appaiata, solo sangue [IC]",
                "+0,080 [+0,021;+0,170]",
                f"{f(v(pdiff, 'diff', scheme='blood', variant='main'), sign=True)} "
                f"{ci(v(pdiff, 'lo', scheme='blood', variant='main'), v(pdiff, 'hi', scheme='blood', variant='main'), sign=True)}",
                "error_rate_paired_diff",
            ),
            (
                "CD8A > 0 nei cloni discordanti / concordanti (Tumor)",
                "7%/88%",
                f"{pct(v(mk, 'frac_CD8A_pos_disc', tissue='Tumor'), 0)}/{pct(v(mk, 'frac_CD8A_pos_conc', tissue='Tumor'), 0)}",
                "marker_test_disc_vs_conc",
            ),
            ("doppietti d1", "0,723", f(d1), "doublet_mixture"),
            ("doppietti d2", "0,790", f(d2), "doublet_mixture"),
            (
                "espansione clonale: ρ errore / ρ eccesso",
                "ρ 0,18/0,28",
                f"ρ {f(v(ce, 'rho', outcome='error', level='clone_block_permutation'), 2)}/"
                f"{f(v(ce, 'rho', outcome='excess', level='clone_block_permutation'), 2)}",
                "clonal_expansion_tests",
            ),
            (
                f"sim (c): falsi positivi su 100 repliche / su {n_reps_c} (confondente 0; 0,5; 1; 2)",
                "2/100",
                "/".join(
                    f"{v(sc_, 'fp_first100_block', confounder=c):.0f}"
                    for c in (0.0, 0.5, 1.0, 2.0)
                )
                + " su 100; "
                + "/".join(
                    pct(v(sc_, "rate_block", confounder=c))
                    for c in (0.0, 0.5, 1.0, 2.0)
                )
                + f" su {n_reps_c}",
                "sim_c_block_permutation",
            ),
        ]
        # valori della versione pubblicata del README (stessi della presentazione)
        h_old = {
            "exhaustion_score": "ρ -0,47, p_Holm 0,53",
            "myeloid_fraction": "ρ -0,54, p_Holm 0,45",
            "caf_fraction": "ρ -0,15, p_Holm 1,00",
            "H4_ni_score": "ρ -0,14, p_Holm 1,00",
        }
        for h in (
            "exhaustion_score",
            "myeloid_fraction",
            "caf_fraction",
            "H4_ni_score",
        ):
            r, p = h_rho_p(h)
            pres.append(
                (
                    f"{hcode[h]}: ρ, p_Holm (famiglia di 4)",
                    h_old[h],
                    f"ρ {f(r, 2)}, p_Holm {f(v(holm, 'p_holm_family_of_4', hypothesis=h), 2)}",
                    "pseudobulk_hypothesis_tests, ni_score_test, ni_score_holm_joint",
                )
            )
        pres += [
            (
                "AUC LOPO malignita' (13)",
                "0,825",
                f(v(mal, "mean", metric="auc")),
                "malignancy_lopo_summary",
            ),
            (
                "duttali adiacenti predette maligne",
                "31,7%",
                pct(v(ext, "frac_pred_malignant", patient="_all")),
                "external_test_summary",
            ),
            (
                "AUC pooled score CNV medio",
                "0,492",
                f(v(loci, "auc_tumor_vs_normal", signal="MEAN_|CNV|")),
                "cnv_loci_auc",
            ),
            (
                "AUC pooled loci PDAC",
                "0,682",
                f(v(loci, "auc_tumor_vs_normal", signal="PDAC_COMPOSITE")),
                "cnv_loci_auc",
            ),
            (
                f"AUC loci PDAC entro paziente (media dei {len(auc_pat)} pazienti valutabili)",
                "0,53",
                f(v(loci_p, "auc_PDAC_COMPOSITE", patient="_mean")),
                "cnv_loci_auc_by_patient",
            ),
            (
                "sim (b): errore PBMC stimato da LOCO con artefatto (moderato / forte)",
                "0,24",
                f"{f(loco_pbmc['artifact_moderate'][0])} / {f(loco_pbmc['artifact_strong'][0])}",
                "sim_b_reference_schemes",
            ),
        ]
        rows = []
        for lab, old, new, src, *raw in pres:
            same = old.replace(" ", "") == new.replace(" ", "")
            esito = "uguale" if same else "cambiato"
            # stessa stima arrotondata diversamente: il valore esatto dista dal vecchio < 0,001
            if not same and raw and abs(raw[0] - float(old.replace(",", "."))) < 1e-3:
                esito = "arrotondamento (valore esatto a meno di 0,001 dal precedente)"
            rows.append(
                [lab, f"**{old}**", new if same else f"**{new}**", esito, f"`{src}`"]
            )

        L.append("# Risultati finali della revisione\n")
        L.append(
            "Generato da `revisione/genera_risultati.py`: ogni numero viene dalla tabella di "
            '`results/tables/` indicata accanto. "Prima" e\' il valore della presentazione (o del README '
            "pubblicato, dove indicato); i metadati di esecuzione in fondo (tempi, memoria) vengono dai log "
            "in `revisione/log/`, non da `results/tables/`.\n"
        )
        L.append("## Numeri della presentazione: prima e dopo\n")
        L.append(md_table(rows, ["Numero", "Prima", "Dopo", "Esito", "Tabella"]))
        L.append(
            "\nNote: lo 0,53 e' l'AUC entro paziente dei loci PDAC (\"≈ 0,53\" nel README pubblicato; "
            "il codice di allora stampava solo le AUC per paziente, quindi la statistica riassuntiva non "
            "e' ricostruibile; ora si riporta la media). I p grezzi di H1-H4 sono nella sezione delle ipotesi. "
            f"Gli IC con {SMALL_K} pazienti o meno (Adjacent_normal, differenza appaiata) non sono "
            "calibrati (vedi Simulazioni).\n"
        )

        # --- dataset
        L.append("## Dataset\n")

        L.append(
            f"scRNA-seq: {n_(v(ds, 'n_cells'))} cellule, {n_(v(ds, 'n_patients'))} pazienti, "
            f"{n_(v(ds, 'n_genes'))} geni, {n_(v(ds, 'n_libraries'))} librerie; PBMC {n_(v(ds, 'n_cells_PBMC'))} cellule "
            f"({n_(v(ds, 'n_patients_PBMC'))} pazienti), Tumor {n_(v(ds, 'n_cells_Tumor'))} ({n_(v(ds, 'n_patients_Tumor'))}), "
            f"Adjacent_normal {n_(v(ds, 'n_cells_Adjacent_normal'))} ({n_(v(ds, 'n_patients_Adjacent_normal'))}) "
            f"(`dataset_summary_sc`). snRNA-seq: {n_(v(dsn, 'n_cells'))} nuclei, {n_(v(dsn, 'n_patients'))} pazienti "
            "(`dataset_summary_sn`).\n"
        )
        if not tasks.empty:
            L.append(
                md_table(
                    [
                        [
                            r.task,
                            f"{r.n_cells:,}".replace(",", "."),
                            str(r.n_classes),
                            str(r.n_patients),
                            r.primary_eval,
                        ]
                        for r in tasks.itertuples()
                    ],
                    ["Task", "Cellule", "Classi", "Pazienti", "Valutazione"],
                )
            )
            L.append("\nFonte: `task_summary_sc`.\n")

        # --- leakage
        L.append("## Leakage: split casuale contro split per paziente (05)\n")
        if not gap.empty:
            rr = []
            for r in gap.itertuples():
                honest = r.fold_lopo if not np.isnan(r.fold_lopo) else r.fold_group_5
                rr.append(
                    [
                        r.task,
                        f(honest),
                        f(r.fold_random_5),
                        f(-r.gap_random_minus_honest, sign=True),
                        f(r.sd_ratio_honest_over_random, 1),
                    ]
                )
            L.append(
                md_table(
                    rr,
                    [
                        "Task",
                        "Per paziente",
                        "Casuale",
                        "Δ (paziente − casuale)",
                        "Rapporto sd",
                    ],
                )
            )
            L.append(
                "\nFonte: `baseline_gap_sc` (per paziente = LOPO per task1a_paired, 5 fold per paziente altrove).\n"
            )

        # --- modelli
        n_tr = zs.n_train_mean.max() if not zs.empty else np.nan
        parts_nf = [d.n_folds for d in (zs, z1) if not d.empty]
        nf = pd.concat(parts_nf) if parts_nf else pd.Series([np.nan])
        zs_folds = f"{nf.min():.0f}-{nf.max():.0f}"
        n_mod = zs[zs.task == "task2_v1"].model.nunique() if not zs.empty else 0
        L.append(
            f"## {n_mod} modelli (11, LOPO, {n_(n_tr)} cellule di training per fold)\n"
        )
        for df, task in ((zs, "task2_v1"), (zs, "task2_hier3"), (z1, "task1a_paired")):
            if df.empty:
                continue
            d = df[df.task == task].sort_values("macro_f1_mean", ascending=False)
            cols = ["Modello", "macro-F1", "sd", "Δ da logreg_l2"]
            has_conv = "n_folds_not_converged" in d.columns
            if has_conv:
                cols.append("fold non convergenti")
            rr = []
            for r in d.itertuples():
                row = [
                    r.model,
                    f(r.macro_f1_mean),
                    f(r.macro_f1_std),
                    f(r.delta_vs_logreg_l2, sign=True),
                ]
                if has_conv:
                    row.append(
                        ""
                        if np.isnan(r.n_folds_not_converged)
                        else f"{r.n_folds_not_converged:.0f}/{r.n_folds}"
                    )
                rr.append(row)
            L.append(
                f"**{task}** (`model_zoo_summary_{'task1' if task == 'task1a_paired' else 'sc'}`)\n"
            )
            L.append(md_table(rr, cols) + "\n")
        for task, df in (
            ("task2_v1", cv1),
            ("task2_hier3", ch3),
            ("task1a_paired", c1),
        ):
            if df.empty:
                continue
            best = df.macro_f1.iloc[0] - df.delta_vs_best.iloc[0]
            n_sig = int((df.p_nadeau_bengio < 0.05).sum())
            L.append(
                f"- {task}: migliore macro-F1 {f(best)}; {n_sig} modelli su {len(df)} con p di Nadeau-Bengio < 0,05 "
                f"rispetto al migliore (`model_comparison_{task}`)."
            )
        k_nb = (
            ", ".join(f"{k:.0f}" for k in sorted(sa_nb.k.unique()))
            if not sa_nb.empty
            else "⟦n.d.⟧"
        )
        L.append(
            f"\nFold per task: {zs_folds}. Con pochi fold (task1a_paired, {z1.n_folds.max() if not z1.empty else np.nan:.0f} fold) "
            "la potenza e' bassa: un p non significativo non dimostra l'equivalenza. I p sono grezzi "
            "(nessuna correzione per i confronti multipli) e il migliore e' scelto sugli stessi fold. "
            f"La calibrazione di Nadeau-Bengio e' verificata in `sim_a` con k = {k_nb} fold stratificati "
            "di uguale dimensione, non con fold LOPO di dimensione diversa.\n"
        )

        # --- feature selection
        L.append("## Selezione delle feature (12)\n")
        for task, fs, st, dg in (
            ("task2_v1", fs_v1, st_v1, dg_v1),
            ("task1a_paired", fs_t1, st_t1, dg_t1),
        ):
            if fs.empty:
                continue
            rr = []
            ks = sorted(fs.k.unique())

            def mark(m, k, txt):
                # F1 e Jaccard a un k maggiore dei geni con punteggio > 0: pareggi arbitrari
                gt = (
                    v(fs, "any_k_gt_nonzero", method=m, k=k)
                    if "any_k_gt_nonzero" in fs
                    else 0
                )
                return f"{txt}*" if gt == 1 else txt

            for m in ("l1_logreg", "rf_import", "chi2", "anova", "hvg", "random"):
                rr.append(
                    [m]
                    + [mark(m, k, f(fsv(fs, m, k))) for k in ks]
                    + [mark(m, k, f(stv(st, m, k))) for k in ks]
                )
            L.append(
                f"**{task}** (`feature_selection_summary_{task}`, `feature_stability_{task}`)\n"
            )
            L.append(
                md_table(
                    rr,
                    ["Metodo"]
                    + [f"F1 k={k}" for k in ks]
                    + [f"Jaccard k={k}" for k in ks],
                )
            )
            if any("*" in c for r in rr for c in r):
                L.append(
                    "\n\\* k maggiore dei geni con punteggio > 0 in almeno un fold "
                    "(`any_k_gt_nonzero`): i geni oltre quelli con punteggio > 0 sono scelti fra "
                    "pareggi in ordine arbitrario, quindi F1 e Jaccard a quel k non sono interpretabili."
                )
            if not dg.empty:
                l1 = dg[dg.method == "l1_logreg"]
                conv = l1.converged.astype(str).str.lower().eq("true")
                L.append(
                    f"\nL1 (SAGA): {int(conv.sum())}/{len(l1)} fold convergenti, {l1.n_nonzero.min():.0f}-{l1.n_nonzero.max():.0f} "
                    f"coefficienti non nulli (`feature_selection_diag_{task}`).\n"
                )
        if not core_t1.empty:
            for m in ("hvg", "chi2", "anova", "rf_import", "l1_logreg"):
                g = core_t1[(core_t1.method == m) & (core_t1.k == k_core)].gene.tolist()
                L.append(
                    f"- {m}: {len(g)} geni in tutti i {n_fold_t1} fold LOPO di task1a_paired "
                    f"(k = {k_core}): {', '.join(g)}"
                )
            L.append("\nFonte: `feature_core_genes_task1a_paired`.\n")

        # --- XAI
        L.append("## Spiegazioni per paziente (20)\n")
        if not xj.empty:
            m = xj.set_index(xj.columns[0]).values
            off = m[np.triu_indices_from(m, 1)]
            jr = f"{f(off.min())}-{f(off.max())}"
        else:
            jr = "⟦n.d.⟧"
        L.append(
            f"Jaccard medio fra i top-{v(xs, 'top_k'):.0f} dei {len(xj)} fold LOPO {f(v(xs, 'jaccard_lopo'))} ({jr} fra coppie, `xai_jaccard_task1a_paired`), "
            f"controllo casuale {f(v(xs, 'jaccard_controllo'))}, differenza {f(v(xs, 'differenza'), sign=True)}, "
            f"Spearman medio {f(v(xs, 'spearman_lopo'))} (`xai_summary_task1a_paired`). "
            f"Geni comuni a tutti i fold ({len(xg)}): {', '.join(xg.gene_symbol.astype(str)) if not xg.empty else ''} "
            "(`xai_common_genes_task1a_paired`). In simulazione la differenza e' "
            f"{f(sim_d_eff, sign=True)} con effetto paziente e {f(sim_d_null, sign=True)} senza (`sim_d_xai_stability`).\n"
        )

        # --- TCR
        L.append("## Validazione clonale con il TCR\n")
        L.append(
            f"- Rumore di base (stesso clone, stesso compartimento): {f(v(exb, 'd_within'))}; discordanza fra compartimenti "
            f"{f(v(exb, 'd_across'))}; eccesso {f(v(exb, 'excess'), sign=True)} IC {ci(v(exb, 'lo_cluster'), v(exb, 'hi_cluster'), sign=True)} "
            f"({v(exb, 'n_patients'):.0f} pazienti, {n_(v(exb, 'n_pairs'))} coppie; `excess_summary_clone_b`)."
        )
        L.append(
            f"- Sensibilita' (B3): pesi uguali {f(v(exb, 'excess_eq'), sign=True)} IC {ci(v(exb, 'lo_eq_cluster'), v(exb, 'hi_eq_cluster'), sign=True)}; "
            f"clone TRA+TRB {f(v(exab, 'excess'), sign=True)} IC {ci(v(exab, 'lo_cluster'), v(exab, 'hi_cluster'), sign=True)} "
            f"(`excess_summary_clone_ab`); split-half di 06 {f(v(sw, 'excess'), sign=True)} IC {ci(v(sw, 'lo'), v(sw, 'hi'), sign=True)} "
            f"({v(sw, 'n_patients'):.0f} pazienti: {v(sw, 'n_patients_across'):.0f} nel termine fra "
            f"compartimenti, {v(sw, 'n_patients_within'):.0f} nel basale; `clone_switching_summary_clone_b`)."
        )
        rr = []
        for sch in ("blood", "loco", "naive"):
            for tis in ("PBMC", "Adjacent_normal", "Tumor"):
                d = (
                    schemes[
                        (schemes.scheme == sch)
                        & (schemes.variant == "main")
                        & (schemes.tissue == tis)
                    ]
                    if not schemes.empty
                    else schemes
                )
                if len(d) == 1:
                    r = d.iloc[0]
                    rr.append(
                        [
                            sch,
                            tis,
                            f(r.error),
                            ci(r.lo, r.hi) + ("*" if r.n_pat <= SMALL_K else ""),
                            n_(r.n),
                            f"{r.n_pat:.0f}",
                        ]
                    )
        L.append(
            "\n**Tasso d'errore per compartimento e schema di riferimento** (`error_rate_schemes`; "
            f"\\* IC non calibrato, <= {SMALL_K} pazienti)\n"
        )
        L.append(
            md_table(
                rr, ["Schema", "Compartimento", "Errore", "IC 95%", "Cloni", "Pazienti"]
            )
        )
        rr = []
        for r in pdiff.itertuples() if not pdiff.empty else []:
            rr.append(
                [
                    r.scheme,
                    r.variant,
                    r.contrast,
                    f(r.diff, sign=True),
                    ci(r.lo, r.hi, sign=True),
                    n_(r.n_clones),
                    f"{r.n_patients:.0f}",
                    ""
                    if np.isnan(r.n_distinct_resamples)
                    else f"{r.n_distinct_resamples:.0f}",
                    "" if np.isnan(r.frac_boot_gt0) else f(r.frac_boot_gt0),
                ]
            )
        L.append(
            f"\n**Differenze appaiate** (`error_rate_paired_diff`; B2: ricampionature distinte e frazione > 0 con <= {SMALL_K} pazienti)\n"
        )
        L.append(
            md_table(
                rr,
                [
                    "Schema",
                    "Variante",
                    "Contrasto",
                    "Diff.",
                    "IC 95%",
                    "Cloni",
                    "Pazienti",
                    "Ricamp. distinte",
                    "Frac > 0",
                ],
            )
        )
        rr = (
            [
                [
                    r.category,
                    n_(r.n_errors),
                    f(r.rate_over_all_clones),
                    n_(r.n_clones),
                    "" if np.isnan(r.lo) else ci(r.lo, r.hi),
                ]
                for r in dec.itertuples()
            ]
            if not dec.empty
            else []
        )
        L.append(
            "\n**Scomposizione dell'errore nel Tumor, schema solo sangue (B1)** (`error_decomposition_blood`)\n"
        )
        L.append(md_table(rr, ["Categoria", "Errori", "Tasso", "Cloni", "IC 95%"]))
        rr = (
            [
                [r.description, r.tissue, f(r.error), ci(r.lo, r.hi), n_(r.n)]
                for r in rules.itertuples()
            ]
            if not rules.empty
            else []
        )
        L.append(
            "\n**Regola d'identita' CD4/CD8 (B6)** (`error_rate_identity_rules`)\n"
        )
        L.append(md_table(rr, ["Regola", "Compartimento", "Errore", "IC 95%", "Cloni"]))
        L.append(
            f"\nMarcatori nelle cellule etichettate CD8T del Tumor (`marker_test_disc_vs_conc`): CD8A > 0 nel "
            f"{pct(v(mk, 'frac_CD8A_pos_disc', tissue='Tumor'))} dei discordanti contro {pct(v(mk, 'frac_CD8A_pos_conc', tissue='Tumor'))} "
            f"dei concordanti; rapporto delle medie di log1p(CP10k) di CD4 {f(v(mk, 'ratio_mean_CD4_disc_over_conc', tissue='Tumor'), 1)}; "
            f"GNLY {f(v(mk, 'mean_GNLY_disc', tissue='Tumor'), 2)} contro {f(v(mk, 'mean_GNLY_conc', tissue='Tumor'), 2)}.\n"
        )
        if not tri.empty:
            rr = [
                [
                    r.contrast,
                    f(r.mean, sign=True),
                    ci(r.lo, r.hi, sign=True),
                    n_(r.n_clones),
                    f"{r.n_patients:.0f}",
                ]
                for r in tri.itertuples()
            ]
            L.append(
                f"**Cloni presenti nei tre compartimenti (08)** (`triple_clones_summary`; {tri.n_patients.max():.0f} pazienti"
                + (": IC non calibrati" if tri.n_patients.max() <= SMALL_K else "")
                + ")\n"
            )
            L.append(
                md_table(rr, ["Contrasto", "Media", "IC 95%", "Cloni", "Pazienti"])
                + "\n"
            )

        L.append("### Doppietti fra le NK con TCR (09, 10)\n")
        if not tcr_rate.empty:
            L.append(
                f"Tassi di TCR produttivo per classe: `tcr_rate_by_celltype` ({len(tcr_rate)} classi).\n"
            )
        rr = (
            [
                [
                    r.statistic,
                    r.stratum,
                    f"{r.n_patients:.0f}",
                    f(r.d1),
                    ci(r.d1_lo, r.d1_hi),
                    f(r.d2),
                    ci(r.d2_lo, r.d2_hi),
                ]
                for r in mix.itertuples()
            ]
            if not mix.empty
            else []
        )
        L.append(
            md_table(
                rr, ["Statistica", "Strato", "Pazienti", "d1", "IC d1", "d2", "IC d2"]
            )
        )
        L.append(
            "\nFonte: `doublet_mixture` (B5: medie e strati per compartimento). d1 e d2 non sono indipendenti; la copertura dei loro IC non e' verificata.\n"
        )

        L.append("### Espansione clonale ed errore (17, 18)\n")
        rr = (
            [
                [
                    r.outcome,
                    r.level,
                    n_(r.n),
                    f(r.rho),
                    p_fmt(r.p),
                    "" if np.isnan(r.null_mean) else f(r.null_mean, 4),
                    "" if np.isnan(r.p_centered) else p_fmt(r.p_centered),
                ]
                for r in ce.itertuples()
            ]
            if not ce.empty
            else []
        )
        L.append(
            md_table(
                rr, ["Esito", "Livello", "n", "ρ", "p", "media nulla", "p centrato"]
            )
        )
        L.append("\nFonte: `clonal_expansion_tests`.\n")
        rr = (
            [
                [
                    r.outcome,
                    str(r.seed),
                    f"{r.n_perm:,}".replace(",", "."),
                    p_fmt(r.p),
                    f(r.null_mean, 4),
                    p_fmt(r.p_centered),
                ]
                for r in cps.itertuples()
            ]
            if not cps.empty
            else []
        )
        L.append(
            "**Sensibilita' del block permutation test (B4)** (`clonal_expansion_perm_sensitivity`)\n"
        )
        L.append(
            md_table(
                rr, ["Esito", "Seme", "Permutazioni", "p", "media nulla", "p centrato"]
            )
        )
        L.append("")
        for nm, ly in (("errore", ly_e), ("eccesso", ly_x)):
            if not ly.empty:
                src = f"loyo_clonal_expansion_{'error' if nm == 'errore' else 'excess'}"
                L.append(
                    f"- Leave-one-patient-out, {nm}: p < 0,05 in {int(ly.still_significant.sum())}/{len(ly)} casi, "
                    f"ρ fra {f(ly.rho_without.min())} e {f(ly.rho_without.max())}; p centrato < 0,05 in "
                    f"{int((ly.p_centered_without < 0.05).sum())}/{len(ly)} (massimo {p_fmt(ly.p_centered_without.max())}) "
                    f"(`{src}`)."
                )
        if not ly_s.empty:
            for n_perm, g in ly_s.groupby("n_perm"):
                L.append(
                    f"- Leave-one-patient-out, errore, {n_(n_perm)} permutazioni ({g.seed.nunique()} {'seme' if g.seed.nunique() == 1 else 'semi'}): p fra "
                    f"{p_fmt(g.p.min())} e {p_fmt(g.p.max())}, p centrato fino a {p_fmt(g.p_centered.max())} "
                    "(`loyo_clonal_expansion_perm_sensitivity`)."
                )
        L.append("\n### Quattro ipotesi per paziente (16, 19)\n")
        rr = []
        for h in (
            "exhaustion_score",
            "myeloid_fraction",
            "caf_fraction",
            "H4_ni_score",
        ):
            r, p = h_rho_p(h)
            rr.append(
                [
                    hcode[h],
                    h,
                    f(r, 3),
                    f(p, 3),
                    f(v(holm, "p_holm_family_of_4", hypothesis=h), 3),
                ]
            )
        L.append(md_table(rr, ["Ipotesi", "Covariata", "ρ", "p", "p_Holm (4 test)"]))
        L.append(
            "\nFonti: `pseudobulk_hypothesis_tests`, `ni_score_test`, `ni_score_holm_joint`. "
            f"n = {v(ni, 'n_patients', hypothesis='H4_ni_score'):.0f} pazienti: potenza bassa.\n"
        )

        # --- CNV
        L.append("## CNV e malignita' (13, 14)\n")
        L.append(
            f"AUC LOPO del Random Forest sull'etichetta CNV: {f(v(mal, 'mean', metric='auc'))} ± {f(v(mal, 'std', metric='auc'))} "
            f"({v(mal, 'count', metric='auc'):.0f} pazienti; `malignancy_lopo_summary`); duttali adiacenti predette maligne "
            f"{pct(v(ext, 'frac_pred_malignant', patient='_all'))} di {v(ext, 'n', patient='_all'):.0f} (`external_test_summary`, "
            "non e' una validazione esterna). "
            f"AUC tumore vs adiacente pooled: score CNV medio {f(v(loci, 'auc_tumor_vs_normal', signal='MEAN_|CNV|'))}, loci PDAC "
            f"{f(v(loci, 'auc_tumor_vs_normal', signal='PDAC_COMPOSITE'))} (`cnv_loci_auc`).\n"
        )
        if not loci_p.empty:
            rr = [
                [
                    r.patient,
                    n_(r.n_tumor),
                    n_(r.n_normal),
                    f(r.auc_PDAC_COMPOSITE),
                    f(getattr(r, "_5")),
                ]
                for r in loci_p.itertuples()
            ]
            L.append(
                md_table(
                    rr,
                    [
                        "Paziente",
                        "Duttali Tumor",
                        "Duttali Adjacent",
                        "AUC loci PDAC",
                        "AUC score medio",
                    ],
                )
            )
            L.append("\nFonte: `cnv_loci_auc_by_patient` (AUC entro paziente).\n")

        # --- simulazioni
        L.append("## Simulazioni\n")
        rr = (
            [
                [
                    r.scenario,
                    n_(r.n),
                    f"{r.k:.0f}",
                    n_(r.reps),
                    f(r.frac_p_lt_005_nadeau_bengio),
                    f(r.frac_p_lt_005_uncorrected),
                ]
                for r in sa_nb.itertuples()
            ]
            if not sa_nb.empty
            else []
        )
        L.append("**(a) Nadeau-Bengio** (`sim_a_nadeau_bengio`)\n")
        L.append(
            md_table(
                rr,
                [
                    "Scenario",
                    "n",
                    "k",
                    "Repliche",
                    "p < 0,05 (NB)",
                    "p < 0,05 (t non corretto)",
                ],
            )
        )
        if not sa_cb.empty:
            piv = sa_cb.pivot_table(
                index=["scenario", "method"], columns="k_patients", values="rate"
            ).reset_index()
            rr = [
                [r["scenario"], r["method"]]
                + [f(r[k]) for k in sorted(sa_cb.k_patients.unique())]
                for _, r in piv.iterrows()
            ]
            L.append(
                "\n**(a) Cluster bootstrap: copertura / frazione di IC che escludono 0** (`sim_a_cluster_bootstrap`)\n"
            )
            L.append(
                md_table(
                    rr,
                    ["Scenario", "Metodo"]
                    + [f"{k} pazienti" for k in sorted(sa_cb.k_patients.unique())],
                )
            )
        if not sb.empty:
            rr = [
                [
                    r.condition,
                    r.scheme,
                    r.tissue,
                    f(r.est_error),
                    f(r.true_error),
                    f(r.bias, sign=True),
                ]
                for r in sb.itertuples()
            ]
            L.append("\n**(b) Schemi di riferimento** (`sim_b_reference_schemes`)\n")
            L.append(
                md_table(
                    rr,
                    [
                        "Condizione",
                        "Schema",
                        "Compartimento",
                        "Errore stimato",
                        "Errore vero",
                        "Distorsione",
                    ],
                )
            )
        if not sc_.empty:
            rr = [
                [
                    f(r.confounder, 1),
                    n_(r.reps),
                    f"{pct(r.rate_block)} {ci(100 * r.lo_block, 100 * r.hi_block, 1)}",
                    f"{pct(r.rate_block_centered)} {ci(100 * r.lo_block_centered, 100 * r.hi_block_centered, 1)}",
                    f"{pct(r.rate_naive)}",
                    f"{r.fp_first100_block:.0f}",
                ]
                for r in sc_.itertuples()
            ]
            L.append(
                "\n**(c) Block permutation sotto H0 con confondente di paziente** (`sim_c_block_permutation`)\n"
            )
            L.append(
                md_table(
                    rr,
                    [
                        "Confondente",
                        "Repliche",
                        "FP block (IC 95%)",
                        "FP block centrato (IC 95%)",
                        "FP ingenuo",
                        "FP nelle prime 100",
                    ],
                )
            )
        if not sd.empty:
            L.append(
                f"\n**(d) Stabilita' XAI**: differenza controllo − LOPO {f(sim_d_eff, sign=True)} con effetto paziente "
                f"(da {f(pe.difference.min(), sign=True)} a {f(pe.difference.max(), sign=True)} sulle {len(pe)} repliche) e {f(sim_d_null, sign=True)} senza "
                "(`sim_d_xai_stability`).\n"
            )

    except (AttributeError, KeyError, IndexError, ValueError, TypeError) as e:
        # tabella o colonna mancante/malformata: si registra e si scrive quanto gia' generato
        MISSING.append(f"interruzione della generazione: {type(e).__name__}: {e}")
    if NOTE.exists():
        L.append("\n" + NOTE.read_text())
    if LIMITS.exists():
        L.append("\n" + LIMITS.read_text())
    if MISSING:
        L.append(
            "\n## Valori mancanti\n\n" + "\n".join(f"- {m}" for m in MISSING) + "\n"
        )
    OUT.write_text("\n".join(L) + "\n")
    print(f"[ok] scritto {OUT}; valori mancanti: {len(MISSING)}")
    for m in MISSING:
        print("  -", m)
    if args.strict and MISSING:
        sys.exit(1)


if __name__ == "__main__":
    main()

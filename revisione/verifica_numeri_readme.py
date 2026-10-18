#!/usr/bin/env python3
"""
verifica_numeri_readme.py - Controllo automatico dei numeri del README contro le tabelle.

Due controlli:
  1. per riga della tabella delle fonti (README §5, "Da dove viene ogni numero"): ogni
     numero della colonna di sinistra deve comparire, all'arrotondamento con cui e'
     scritto, in una delle tabelle citate nella colonna di destra (valore, valore x 100
     per le percentuali, valore con segno cambiato; per gli interi anche numero di
     righe, somme di colonna e numero di valori distinti);
  2. su tutto il testo delle sezioni 1, 6 e 7: ogni numero decimale (con la virgola o
     con il %) deve comparire, allo stesso arrotondamento, in almeno una tabella di
     `results/tables/` (controllo piu' debole: dice solo che il numero esiste).
Stampa i numeri non trovati, da verificare a mano (parametri come "20.000
permutazioni", soglie, numeri della versione precedente).

Uso
---
    python revisione/verifica_numeri_readme.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results/tables"
README = ROOT / "README.md"

NUM_RE = re.compile(
    r"(?<![\w.,])[−-]?\+?\d{1,3}(?:\.\d{3})+(?:,\d+)?%?|(?<![\w.,])[−+-]?\d+(?:,\d+)?%?"
)


def parse(tok: str) -> tuple[float, int, bool]:
    """Valore, numero di decimali e flag percentuale di un numero scritto all'italiana."""
    pct = tok.endswith("%")
    t = tok.rstrip("%").replace("−", "-").replace("+", "")
    t = t.replace(".", "").replace(",", ".")
    dec = len(t.split(".")[1]) if "." in t else 0
    return float(t), dec, pct


def table_values(names: list[str]) -> tuple[np.ndarray, set]:
    """Valori numerici delle tabelle e conteggi derivati (righe, somme, valori distinti)."""
    vals, counts = [], set()
    for name in names:
        for path in sorted(TAB.glob(name + ".csv")):
            try:
                df = pd.read_csv(path, low_memory=False)
            except Exception:
                continue
            counts.add(len(df))
            # righe di riepilogo ('_mean', '_all') escluse dal conteggio
            first = df.iloc[:, 0].astype(str) if df.shape[1] else pd.Series(dtype=str)
            counts.add(int((~first.str.startswith("_")).sum()))
            # somme e valori distinti per gruppo (colonne con pochi valori, es. tessuto)
            for g in [c for c in df.columns if 1 < df[c].nunique() <= 10][:5]:
                for _, sub in df.groupby(g, observed=True):
                    for c in df.columns:
                        counts.add(sub[c].nunique())
                        s_ = pd.to_numeric(sub[c], errors="coerce")
                        if s_.notna().any():
                            counts.add(float(np.round(s_.sum(), 6)))
            for c in df.columns:
                s = pd.to_numeric(df[c], errors="coerce")
                if s.notna().any():
                    v = s.dropna().to_numpy(dtype=float)
                    vals.append(v)
                    if len(v) <= 5000:
                        counts.add(float(np.round(v.sum(), 6)))
                counts.add(df[c].nunique())
    allv = np.concatenate(vals) if vals else np.array([])
    return allv, counts


def found(tok: str, vals: np.ndarray, counts: set, scale100: bool = False) -> bool:
    """Il numero compare (allo stesso arrotondamento) fra i valori o i conteggi?

    x 100 per le percentuali, per i punti percentuali (un decimale) e, con
    `scale100`, per ogni numero (controllo debole sul testo).
    """
    x, dec, pct = parse(tok)
    cands = [vals, -vals]
    if pct or dec == 1 or scale100:
        cands += [vals * 100, -vals * 100]
    tol = 0.5 * 10**-dec + 1e-9
    for c in cands:
        if len(c) and np.any(np.abs(c - x) <= tol):
            return True
    if dec == 0 and (x in counts or abs(x) in counts):
        return True
    return False


def source_rows(text: str) -> list[tuple[str, list[str]]]:
    """Righe della tabella delle fonti: (numeri, tabelle citate)."""
    start = text.index("Da dove viene ogni numero del README")
    end = text.index("## 6.", start)
    rows = []
    for line in text[start:end].splitlines():
        if not line.startswith("| ") or line.startswith("| Numero") or ":---" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        left, right = cells[0], " ".join(cells[1:])
        names = re.findall(r"`([A-Za-z0-9_{},*]+)`", right)
        exp = []
        for n in names:
            m = re.match(r"(.*)\{(.*)\}(.*)", n)
            if m:
                exp += [m.group(1) + o + m.group(3) for o in m.group(2).split(",")]
            elif n == "sim_*":
                exp.append("sim_*")
            elif re.match(r"^[a-z][a-z0-9_]+$", n):
                exp.append(n)
        rows.append((left, exp))
    return rows


def split_ranges(s: str) -> str:
    """'0,53-0,61' -> '0,53 0,61'; '38-65%' -> '38% 65%' (il % vale per i due estremi)."""
    s = re.sub(r"(\d+(?:,\d+)?)-(\d+(?:,\d+)?)%", r"\1% \2%", s)
    return re.sub(r"(\d)-(\d)", r"\1 \2", s)


def main() -> None:
    text = split_ranges(README.read_text())
    print("== 1. tabella delle fonti ==")
    n_bad = 0
    for left, names in source_rows(text):
        if not names:
            continue
        vals, counts = table_values(names)
        clean = re.sub(r"`[^`]*`", " ", left)
        toks = [t for t in NUM_RE.findall(clean)]
        miss = [t for t in toks if not found(t, vals, counts)]
        if miss:
            n_bad += 1
            print(f"  [?] {left[:70]!r} -> {names}: non trovati {miss}")
    print(f"  righe con numeri non trovati: {n_bad}")

    print("\n== 2. testo delle sezioni 1, 6, 7 (solo decimali e percentuali) ==")
    all_vals, _ = table_values(["*"])
    sec = []
    for a, b in (("## 1.", "## 2."), ("## 6.", "### Limiti noti del codice")):
        sec.append(text[text.index(a) : text.index(b)])
    body = re.sub(r"`[^`]*`", " ", "\n".join(sec))
    toks = [t for t in NUM_RE.findall(body) if "," in t or t.endswith("%")]
    miss = sorted({t for t in toks if not found(t, all_vals, set(), scale100=True)})
    print(
        f"  numeri controllati: {len(toks)}, non trovati in nessuna tabella: {len(miss)}"
    )
    for t in miss:
        i = body.find(t)
        ctx = body[max(0, i - 60) : i + 20].replace("\n", " ")
        print(f"  [?] {t:>10s}  ...{ctx}...")


if __name__ == "__main__":
    main()

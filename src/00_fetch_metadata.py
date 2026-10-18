#!/usr/bin/env python3
"""
00_fetch_metadata.py - Metadati GEO della SuperSeries GSE278694 e verifica di fattibilita' del disegno.

Scopo
-----
Scarica solo gli archivi descrittivi leggeri della SuperSeries GSE278694:
  - GSE278694_family.xml.tgz  archivio MINiML (titolo, caratteristiche e file
                              supplementari di ogni campione GSM);
  - filelist.txt              elenco dei file supplementari della serie (serve a
                              download_tcr.sh e 06 per i 30 file VDJ).
Classifica ogni campione per modalita' (TCR, Visium, snRNA, scRNA) e paziente con
regole lessicali sul titolo (es. "PA#01, Tumor, scRNA-seq" -> PA01; riconosce anche
le forme PDAC01, PA_01, PA01) e risponde a tre domande:
  Q1  scRNA-seq e snRNA-seq sono sugli stessi pazienti?
  Q2  le librerie scTCR-seq V(D)J sono appaiate a librerie GEX dello stesso paziente?
  Q3  i file supplementari contengono annotazioni di tipo cellulare?

Input
-----
  rete: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278694/ (i file gia'
  presenti in --outdir non vengono riscaricati)

Output
------
  {outdir}/GSE278694_family.xml.tgz, {outdir}/filelist.txt (cache dei download)
  {outdir}/samples.tsv   un campione per riga: gsm, titolo, modalita', paziente, file

Argomenti
---------
  --outdir (default data/raw/metadata)

Uso
---
    python src/00_fetch_metadata.py --outdir data/raw/metadata
"""

from __future__ import annotations

import argparse
import io
import re
import sys
import tarfile
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd

# Identificatore NCBI GEO della SuperSeries e costruzione degli endpoint FTP ufficiali.
GSE = "GSE278694"
STUB = GSE[:-3] + "nnn"  # Prefisso di partizionamento directory GEO (GSE278nnn)
BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{STUB}/{GSE}"
MINIML_URL = f"{BASE}/miniml/{GSE}_family.xml.tgz"
FILELIST_URL = f"{BASE}/suppl/filelist.txt"

# Namespace XML standard per lo schema MINiML di NCBI GEO.
NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

# Espressioni regolari ordinate per specificita decrescente per inferire la modalita di sequenziamento.
MODALITY_RULES = [
    ("TCR", r"\btcr\b|vdj|v\(d\)j|immune[\s_-]*profil"),
    ("Visium", r"visium|spatial|cytassist|space[\s_-]*ranger|\bst\b"),
    ("snRNA", r"\bsn(rna)?\b|single[\s_-]*nucle|nuclei|\bnuc\b"),
    ("scRNA", r"\bsc(rna)?\b|single[\s_-]*cell|\bgex\b|\bcell\b"),
]

# Pattern lessicali indicativi di tabelle contenenti metadati o annotazioni cellulari.
ANNOTATION_HINTS = (
    r"meta|annot|celltype|cell[\s_-]*type|label|cluster|barcode|obs\b|seurat|umap"
)


def fetch(url: str, timeout: int = 120) -> bytes:
    """Esegue il download HTTP/FTP di una risorsa remota con gestione esplicita degli errori di rete."""
    req = urllib.request.Request(url, headers={"User-Agent": "pdac-ml/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError as e:
        sys.exit(f"[ERRORE] HTTP {e.code} su {url}\n{e.reason}")
    except urllib.error.URLError as e:
        sys.exit(f"[ERRORE] rete non raggiungibile: {url}\n{e.reason}")


def load_miniml(outdir: Path) -> ET.Element:
    """Scarica (o recupera dalla cache locale) l'archivio MINiML e restituisce il nodo radice dell'albero XML."""
    tgz = outdir / f"{GSE}_family.xml.tgz"
    if not tgz.exists():
        print(f"[..] scarico {MINIML_URL}")
        tgz.write_bytes(fetch(MINIML_URL))
    else:
        print(f"[ok] uso cache {tgz}")

    # Estrazione in streaming del documento XML dall'archivio compresso tar.gz senza scritture intermedie su disco.
    with tarfile.open(tgz, "r:gz") as tar:
        member = next((m for m in tar.getmembers() if m.name.endswith(".xml")), None)
        if member is None:
            sys.exit("[ERRORE] nessun .xml dentro l'archivio MINiML")
        fh = tar.extractfile(member)
        assert fh is not None
        return ET.parse(io.BytesIO(fh.read())).getroot()


def text(node: ET.Element | None) -> str:
    """Estrae e normalizza il contenuto testuale di un nodo XML, restituendo stringa vuota se assente."""
    return (node.text or "").strip() if node is not None else ""


def parse_samples(root: ET.Element) -> pd.DataFrame:
    """Costruisce un DataFrame strutturato con un record per ciascun campione GSM censito nel documento MINiML."""
    rows = []
    for s in root.findall("m:Sample", NS):
        gsm = s.attrib.get("iid", "")
        title = text(s.find("m:Title", NS))
        strategy = text(s.find("m:Library-Strategy", NS))
        source = " ".join(text(n) for n in s.findall(".//m:Source", NS))
        descr = " ".join(text(n) for n in s.findall("m:Description", NS))

        # Estrazione dinamica dei tag di caratteristica clinico-sperimentale associati al campione.
        chars: dict[str, str] = {}
        for c in s.findall(".//m:Characteristics", NS):
            tag = (c.attrib.get("tag") or "char").strip().lower()
            tag = re.sub(r"\s+", "_", tag)
            chars[f"ch_{tag}"] = (c.text or "").strip()

        supp = [text(n) for n in s.findall("m:Supplementary-Data", NS)]

        rows.append(
            {
                "gsm": gsm,
                "title": title,
                "library_strategy": strategy,
                "source": source,
                "description": descr,
                "n_supp_files": len(supp),
                "supp_files": " | ".join(Path(u).name for u in supp if u),
                **chars,
            }
        )
    if not rows:
        sys.exit("[ERRORE] nessun campione trovato nel MINiML.")
    return pd.DataFrame(rows)


def infer_modality(row: pd.Series) -> str:
    """Inferisce la modalita sperimentale del campione applicando in sequenza le regole lessicali su campi multipli."""
    hay = " ".join(
        str(row.get(k, ""))
        for k in ("title", "library_strategy", "source", "description", "supp_files")
    ).lower()
    # L'ordine delle regole preserva la priorita delle librerie specifiche (TCR, Visium) rispetto a scRNA/snRNA.
    for name, pattern in MODALITY_RULES:
        if re.search(pattern, hay):
            return name
    return "UNKNOWN"


def infer_patient(row: pd.Series) -> str:
    """Identifica il codice del paziente dai campi strutturati 'Characteristics' o tramite espressione regolare sul titolo."""
    for k, v in row.items():
        if str(k).startswith("ch_") and re.search(
            r"patient|subject|donor|individual|case|sample[_ ]?id", str(k)
        ):
            v = str(v).strip()
            if v:
                return v
    # Forme riconosciute: "PA#01", "PDAC01", "PA_01", "PA01", "Patient 1" -> PA01
    m = re.search(
        r"(?<![A-Za-z])(?:PDAC|PA|PT|Patient|P)[\s#_-]*0*(\d{1,3})(?!\d)",
        str(row.get("title", "")),
        re.I,
    )
    return f"PA{int(m.group(1)):02d}" if m else "UNKNOWN"


def load_filelist(outdir: Path) -> pd.DataFrame | None:
    """Scarica e legge il manifesto filelist.txt dei file supplementari della serie GEO."""
    fl = outdir / "filelist.txt"
    if not fl.exists():
        print(f"[..] scarico {FILELIST_URL}")
        try:
            fl.write_bytes(fetch(FILELIST_URL))
        except SystemExit:
            print("[!!] filelist.txt non disponibile (non tutte le serie ce l'hanno)")
            return None
    try:
        return pd.read_csv(fl, sep="\t")
    except Exception as e:  # noqa: BLE001
        print(f"[!!] filelist.txt non parsabile: {e}")
        return None


def main() -> None:
    """Esegue la pipeline di acquisizione metadati, classificazione dei campioni GSM e verifica di fattibilita."""
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Metadati GEO di GSE278694 e verifica di fattibilita' del disegno.",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("data/raw/metadata"),
        help="cartella dei metadati GEO (filelist.txt, tabelle dei campioni)",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Caricamento dell'albero XML MINiML ed estrazione tabellare dei metadati per campione.
    root = load_miniml(args.outdir)
    df = parse_samples(root)
    df["modality"] = df.apply(infer_modality, axis=1)
    df["patient"] = df.apply(infer_patient, axis=1)

    out_tsv = args.outdir / "samples.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"\n[ok] {len(df)} campioni -> {out_tsv}\n")
    print("=" * 78)
    print(f"GSE278694: {len(df)} GSM")
    print("=" * 78)

    print("\nModalita inferite (verificare in samples.tsv):")
    print(df["modality"].value_counts().to_string())

    unknown = (df["modality"] == "UNKNOWN").sum()
    if unknown:
        print(f"\n[!!] {unknown} campioni non classificati. Titoli:")
        print(df.loc[df.modality == "UNKNOWN", ["gsm", "title"]].to_string(index=False))

    # Verifica Q1: sovrapposizione per paziente tra profili scRNA-seq e snRNA-seq.
    print("\n" + "-" * 78)
    print("Q1  scRNA e snRNA sugli STESSI pazienti?   [analisi cross-modalita]")
    print("-" * 78)
    pivot = pd.crosstab(df["patient"], df["modality"])
    print(pivot.to_string())

    have = {m for m in ("scRNA", "snRNA") if m in pivot.columns}
    if have == {"scRNA", "snRNA"}:
        paired = pivot[(pivot["scRNA"] > 0) & (pivot["snRNA"] > 0)]
        n = len(paired.drop(index="UNKNOWN", errors="ignore"))
        print(f"\n=> {n} pazienti con ENTRAMBE le modalita.")
        if n >= 3:
            print(
                "   ESITO: confronto cross-modalita scRNA vs snRNA applicabile su soggetti appaiati."
            )
        elif n >= 1:
            print(
                "   ESITO: confronto cross-modalita limitato a pochi soggetti appaiati."
            )
        else:
            print(
                "   ESITO: coorti disgiunte; la modalita tecnica risulta confusa col paziente."
            )
            print(
                "          Escludere il trasferimento diretto scRNA -> snRNA non controllato."
            )
    else:
        print(
            f"\n=> Modalita RNA trovate: {have or 'nessuna'}. Controlla le euristiche."
        )

    # Verifica Q2: presenza di librerie V(D)J TCR appaiate a librerie trascrittomiche GEX.
    print("\n" + "-" * 78)
    print("Q2  Librerie TCR appaiate a librerie GEX?  [validazione clonale V(D)J]")
    print("-" * 78)
    if "TCR" in pivot.columns:
        gex_cols = [c for c in ("scRNA", "snRNA") if c in pivot.columns]
        tcr_pat = set(pivot.index[pivot["TCR"] > 0])
        gex_pat = (
            set(pivot.index[pivot[gex_cols].sum(axis=1) > 0]) if gex_cols else set()
        )
        both = (tcr_pat & gex_pat) - {"UNKNOWN"}
        print(
            f"pazienti con TCR: {len(tcr_pat)} | con GEX: {len(gex_pat)} | con entrambi: {len(both)}"
        )
        print(f"   -> {sorted(both)}")
        if both:
            print("\n   ESITO: coorte con doppia libreria TCR + GEX disponibile.")
            print(
                "   Procedere con la verifica di intersezione dei barcode cellulari 16-nt nello step 01."
            )
        else:
            print("\n   ESITO: nessun paziente con doppia libreria TCR + GEX.")
    else:
        print("Nessuna libreria TCR riconosciuta. Controlla i titoli in samples.tsv")
        print("prima di concludere: potrebbe essere un problema di euristica.")

    # Verifica Q3: censimento dei file supplementari contenenti annotazioni di tipo cellulare.
    print("\n" + "-" * 78)
    print(
        "Q3  Annotazioni di cell type nei supplementary?  [classificazione citologica]"
    )
    print("-" * 78)
    fl = load_filelist(args.outdir)
    candidates: list[str] = []
    if fl is not None:
        namecol = next((c for c in fl.columns if "name" in c.lower()), fl.columns[0])
        for n in fl[namecol].astype(str):
            if re.search(ANNOTATION_HINTS, n, re.I):
                candidates.append(n)
        print(f"[ok] filelist.txt: {len(fl)} file di serie")

    per_sample = defaultdict(list)
    for _, r in df.iterrows():
        for f in str(r["supp_files"]).split(" | "):
            if (
                f
                and re.search(ANNOTATION_HINTS, f, re.I)
                and not f.endswith((".h5", ".tar.gz"))
            ):
                per_sample[r["gsm"]].append(f)

    if candidates:
        print("\nCandidati a livello di SERIE:")
        for c in sorted(set(candidates)):
            print(f"   - {c}")
    if per_sample:
        print("\nCandidati a livello di CAMPIONE (primi 10):")
        for gsm, fs in list(per_sample.items())[:10]:
            print(f"   - {gsm}: {', '.join(fs)}")

    if not candidates and not per_sample:
        print("\n[!!] Nessun file di annotazione evidente nella serie principale GEO.")
        print(
            "     ESITO: verificare le SubSeries dedicate (GSE278688, GSE278689, GSE300435)."
        )
    else:
        print("\n   ESITO: individuati file candidati per i metadati cellulari.")
        print(
            "   Nota: i file 'barcodes.tsv' contengono gli identificatori di colonna della matrice, non le etichette."
        )

    print("\n" + "=" * 78)
    print(f"Prossimo passo: verificare {out_tsv} per la corrispondenza tra")
    print("le colonne `modality`, `patient` e `title`.")
    print("=" * 78)


if __name__ == "__main__":
    main()

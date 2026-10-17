#!/usr/bin/env bash
# download_tcr.sh - Download dei 30 file VDJ (filtered_contig_annotations) di GSE300435.
#
# Scopo:   scarica da NCBI GEO i file *_VDJ_filtered_contig_annotations.csv.gz elencati in
#          filelist.txt, con retry e ripresa dei download parziali. Idempotente: i file gia'
#          presenti e integri (gzip -t) non vengono riscaricati.
# Input:   FILELIST (default data/raw/metadata/filelist.txt, da 00_fetch_metadata.py)
# Output:  OUTDIR/<GSM>_PA_<nn>-<tessuto>-VDJ_filtered_contig_annotations.csv.gz (30 file)
# Errori:  esce con codice 1 se il filelist non elenca esattamente 30 file VDJ, se un
#          download fallisce o se alla fine in OUTDIR non ci sono 30 file validi.
# Uso:     bash src/download_tcr.sh [FILELIST] [OUTDIR]
#          bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr

set -euo pipefail

N_EXPECTED=30
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "uso: bash src/download_tcr.sh [FILELIST] [OUTDIR]"
  echo "     scarica i ${N_EXPECTED} file VDJ di GSE300435 elencati in FILELIST"
  echo "     (default data/raw/metadata/filelist.txt) dentro OUTDIR (default data/raw/tcr)"
  exit 0
fi

FILELIST="${1:-data/raw/metadata/filelist.txt}"
OUTDIR="${2:-data/raw/tcr}"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"

[[ -f "$FILELIST" ]] || { echo "[ERRORE] filelist assente: $FILELIST" >&2; exit 1; }
mkdir -p "$OUTDIR"

# controllo DNS prima di iniziare: meglio fallire subito che 30 volte
if ! getent hosts ftp.ncbi.nlm.nih.gov >/dev/null 2>&1; then
  echo "[ERRORE] DNS non risolve ftp.ncbi.nlm.nih.gov" >&2
  echo "         prova: sudo systemctl restart systemd-resolved" >&2
  exit 1
fi

mapfile -t FILES < <(awk -F'\t' '$2 ~ /VDJ_filtered_contig_annotations\.csv\.gz$/ {print $2}' "$FILELIST")
echo "[ok] ${#FILES[@]} file VDJ da scaricare"
if [[ ${#FILES[@]} -ne $N_EXPECTED ]]; then
  echo "[ERRORE] attesi $N_EXPECTED file VDJ nel filelist, trovati ${#FILES[@]}" >&2
  exit 1
fi

ok=0; fail=0
for f in "${FILES[@]}"; do
  gsm="${f%%_*}"                 # GSM9060613
  stub="${gsm:0:$((${#gsm}-3))}nnn"   # GSM9060nnn
  url="$BASE/$stub/$gsm/suppl/$f"
  dest="$OUTDIR/$f"

  if [[ -s "$dest" ]] && gzip -t "$dest" 2>/dev/null; then
    echo "[ok]   cache  $f"; ok=$((ok+1)); continue
  fi

  printf '[..]   %s ... ' "$f"
  if curl -fsSL --retry 5 --retry-delay 2 --retry-connrefused \
          --connect-timeout 20 --max-time 300 -C - -o "$dest" "$url"; then
    if gzip -t "$dest" 2>/dev/null; then
      echo "ok"; ok=$((ok+1))
    else
      echo "CORROTTO"; rm -f "$dest"; fail=$((fail+1))
    fi
  else
    echo "FALLITO"; rm -f "$dest"; fail=$((fail+1))
  fi
done

echo
echo "[ok] scaricati/validi: $ok | falliti: $fail"
[[ $fail -eq 0 ]] || { echo "[ERRORE] rilancia lo script: riprende da dove ha lasciato" >&2; exit 1; }
n_ok=$(find "$OUTDIR" -name '*VDJ_filtered_contig_annotations.csv.gz' | wc -l)
if [[ $n_ok -ne $N_EXPECTED ]]; then
  echo "[ERRORE] in $OUTDIR ci sono $n_ok file VDJ, attesi $N_EXPECTED" >&2
  exit 1
fi
du -sh "$OUTDIR"

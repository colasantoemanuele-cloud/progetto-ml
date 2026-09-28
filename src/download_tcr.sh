#!/usr/bin/env bash
# download_tcr.sh: acquisizione automatizzata dei 30 file VDJ di GSE300435 con retry e ripresa.
#
# Disaccoppiato dall'analisi statistica per garantire l'integrita' dei file di input.
# Idempotente e rilanciabile: verifica l'integrita' gzip dei file in cache e riprende i download parziali.
#
#   bash src/download_tcr.sh data/raw/metadata/filelist.txt data/raw/tcr

set -euo pipefail

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
[[ ${#FILES[@]} -eq 30 ]] || echo "[!!] attesi 30, trovati ${#FILES[@]}"

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
du -sh "$OUTDIR"

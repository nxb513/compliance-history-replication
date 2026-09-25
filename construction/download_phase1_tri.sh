#!/usr/bin/env bash
# Phase 1 TRI Basic Data Files, national (US), calendar years 2010-2024.
set -u
ROOT="."; LOG="$ROOT/logs/download_phase1_tri.log"; SUM="$ROOT/logs/checksums_phase1.csv"
B="https://data.epa.gov/api-cache/downloads/tri/mv_tri_basic_download"
[ -f "$SUM" ] || echo "file,path,url,bytes,sha256,last_modified,downloaded_utc" > "$SUM"
echo "=== TRI phase1 start $(date -u +%FT%TZ) ===" >> "$LOG"
for y in $(seq 2010 2024); do
  url="$B/${y}_us.csv"; out="$ROOT/data/raw/tri/TRI_${y}_US.csv"
  echo "--- $(date -u +%H:%M:%SZ) TRI $y" >> "$LOG"
  curl -sS --fail --retry 5 --retry-delay 5 -o "$out" "$url" >> "$LOG" 2>&1 || { echo "  ERR $y" >> "$LOG"; continue; }
  bytes=$(stat -c %s "$out"); sha=$(sha256sum "$out" | awk '{print $1}')
  echo "TRI_${y}_US.csv,$out,$url,$bytes,$sha,\"\",$(date -u +%FT%TZ)" >> "$SUM"
  echo "  OK $y bytes=$bytes" >> "$LOG"
done
echo "=== TRI phase1 done $(date -u +%FT%TZ) ===" >> "$LOG"; echo "TRI_PHASE1_DONE" >> "$LOG"

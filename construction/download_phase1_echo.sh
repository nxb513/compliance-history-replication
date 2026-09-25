#!/usr/bin/env bash
# Phase 1 ECHO bulk downloads (FRS, ICIS-FEC, DMR FY2010-2014). Raw files kept immutable.
set -u
ROOT="."; LOG="$ROOT/logs/download_phase1_echo.log"; SUM="$ROOT/logs/checksums_phase1.csv"
B="https://echo.epa.gov/files/echodownloads"
[ -f "$SUM" ] || echo "file,path,url,bytes,sha256,last_modified,downloaded_utc" > "$SUM"
echo "=== ECHO phase1 start $(date -u +%FT%TZ) ===" >> "$LOG"
dl () { # $1 url  $2 destdir
  local url="$1" dest="$2" f; f=$(basename "$url"); local out="$dest/$f"
  echo "--- $(date -u +%H:%M:%SZ) $f" >> "$LOG"
  curl -sS --fail --retry 5 --retry-delay 5 -C - -o "$out" "$url" >> "$LOG" 2>&1 || { echo "  ERR $f" >> "$LOG"; return; }
  local bytes sha lm; bytes=$(stat -c %s "$out"); sha=$(sha256sum "$out" | awk '{print $1}')
  lm=$(curl -sI "$url" | grep -i last-modified | tr -d '\r' | sed 's/^[Ll]ast-[Mm]odified: //')
  echo "$f,$out,$url,$bytes,$sha,\"$lm\",$(date -u +%FT%TZ)" >> "$SUM"
  echo "  OK $f bytes=$bytes" >> "$LOG"
}
dl "$B/frs_downloads.zip"  "$ROOT/data/raw/frs"
dl "$B/case_downloads.zip" "$ROOT/data/raw/icis_fec"
for y in 2010 2011 2012 2013 2014; do dl "$B/npdes_dmrs_fy${y}.zip" "$ROOT/data/raw/echo"; done
echo "=== ECHO phase1 done $(date -u +%FT%TZ) ===" >> "$LOG"; echo "ECHO_PHASE1_DONE" >> "$LOG"

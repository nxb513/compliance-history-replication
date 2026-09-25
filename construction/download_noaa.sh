#!/usr/bin/env bash
# Sequential download of NOAA NCEI Storm Events detail files 2009-2025 (author-approved 2026-09-25), with size and gzip checks.
set -u
cd "$(dirname "$0")"
B="https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles"
LIST=$(curl -s --max-time 60 "$B/" | grep -o 'StormEvents_details-ftp_v1.0_d20\(09\|1[0-9]\|2[0-5]\)_c[0-9]*.csv.gz' | sort -u)
for f in $LIST; do
  len=$(curl -sI --max-time 60 "$B/$f" | grep -i content-length | tr -dc '0-9')
  ok=0
  for a in 1 2 3 4; do
    curl -s --max-time 600 -o "$f" "$B/$f"
    got=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [ "$got" -eq "$len" ] && gzip -t "$f" 2>/dev/null; then ok=1; break; fi
    sleep 5
  done
  if [ "$ok" -eq 1 ]; then echo "$(date +%T) DONE $f $got"; else echo "$(date +%T) ERROR $f"; fi
done
sha256sum StormEvents_details-*.csv.gz > SHA256SUMS.txt
echo "$(date +%T) ALL_DONE $(ls StormEvents_details-*.csv.gz | wc -l) files"

#!/usr/bin/env bash
# Sequential ranged download (one connection at a time) of EPA ECHO bulk files, with integrity checks.
set -u
cd "$(dirname "$0")"
BASE="https://echo.epa.gov/files/echodownloads"
CH=10000000
for f in ICIS-AIR_downloads.zip rcra_downloads.zip echo_exporter.zip; do
  len=$(curl -sIL --max-time 60 "$BASE/$f" | grep -i '^content-length' | tail -1 | tr -dc '0-9')
  echo "$(date +%T) START $f length=$len"
  rm -f "$f" "$f".part.*
  start=0; i=0
  while [ "$start" -lt "$len" ]; do
    end=$((start + CH - 1)); [ "$end" -ge "$len" ] && end=$((len - 1))
    part=$(printf "%s.part.%05d" "$f" "$i"); want=$((end - start + 1)); ok=0
    for attempt in 1 2 3 4 5 6; do
      curl -s --max-time 120 -r "$start-$end" -o "$part" "$BASE/$f"
      got=$(stat -c %s "$part" 2>/dev/null || echo 0)
      if [ "$got" -eq "$want" ]; then ok=1; break; fi
      sleep 5
    done
    if [ "$ok" -ne 1 ]; then echo "$(date +%T) ERROR $f chunk $i failed after retries"; exit 1; fi
    start=$((end + 1)); i=$((i + 1))
    [ $((i % 10)) -eq 0 ] && echo "$(date +%T) PROGRESS $f $start/$len"
  done
  cat "$f".part.* > "$f" && rm -f "$f".part.*
  size=$(stat -c %s "$f")
  if [ "$size" -ne "$len" ]; then echo "$(date +%T) ERROR $f size $size != $len"; exit 1; fi
  if unzip -tq "$f" > /dev/null 2>&1; then echo "$(date +%T) DONE $f size=$size zip_ok"; else echo "$(date +%T) ERROR $f zip test failed"; exit 1; fi
done
sha256sum ICIS-AIR_downloads.zip rcra_downloads.zip echo_exporter.zip > SHA256SUMS.txt
echo "$(date +%T) ALL_DONE"; cat SHA256SUMS.txt

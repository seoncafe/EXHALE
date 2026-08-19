#!/bin/bash
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for p in hd189 hd209 wasp121 wasp52; do
  cd "$BASE/insitu_$p"
  echo "[$(date +%H:%M)] START $p" >> "$BASE/insitu_driver.log"
  mpirun -n 72 /home/kiseon/LaRT/combine/LaRT_v2.00/LaRT_calcJPP.x ${p}_insitu.in > lart.log 2>&1
  echo "[$(date +%H:%M)] DONE $p (h5=$(ls ${p}_insitu.h5 2>/dev/null))" >> "$BASE/insitu_driver.log"
done
echo "ALL DONE" >> "$BASE/insitu_driver.log"

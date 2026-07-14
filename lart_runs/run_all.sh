#!/bin/bash
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALCP=/nfs/lart4/kiseon/LaRT/combine/LaRT_v2.00/LaRT_calcP.x
for p in hd209 hd189 wasp121 wasp52; do
  cd "$BASE/$p"
  echo "[$(date +%H:%M)] START $p" >> "$BASE/driver.log"
  mpirun -n 72 $CALCP ${p}_lya.in > lart.log 2>&1
  echo "[$(date +%H:%M)] DONE $p (h5=$(ls ${p}_lya.h5 2>/dev/null))" >> "$BASE/driver.log"
done
echo "ALL DONE" >> "$BASE/driver.log"

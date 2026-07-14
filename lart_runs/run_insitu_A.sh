#!/bin/bash
# Driver A: in-situ Lya LaRT for wasp121 + wasp52 (winds unchanged today).
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for p in wasp121 wasp52; do
  cd $BASE/insitu_$p
  echo "[$(date +%H:%M)] START insitu $p" >> $BASE/insitu_driver.log
  mpirun -n 72 /nfs/lart4/kiseon/LaRT/combine/LaRT_v2.00/LaRT_calcP.x ${p}_insitu.in > lart.log 2>&1
  echo "[$(date +%H:%M)] DONE insitu $p (h5=$(ls ${p}_insitu.h5 2>/dev/null))" >> $BASE/insitu_driver.log
done
echo "[$(date +%H:%M)] DRIVER A DONE" >> $BASE/insitu_driver.log

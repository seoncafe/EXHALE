#!/bin/bash
# Driver C: the two hd209 LaRT runs that failed when the shared LaRT_calcP.x
# vanished (rebuilt from commit dab6f28 = same physics as all completed runs).
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC="$BASE/LaRT_build_dab6f28/LaRT_calcP.x"
run() {
  cd $BASE/$1
  echo "[$(date +%m-%d\ %H:%M)] START $3" >> $BASE/insitu_driver.log
  mpirun -n 72 $EXEC $2 > lart.log 2>&1
  echo "[$(date +%m-%d\ %H:%M)] DONE $3 (exit $?)" >> $BASE/insitu_driver.log
}
run insitu_hd209 hd209_insitu.in "insitu hd209 (C)"
run bm_hd209     hd209_lya.in    "stellar hd209 (C)"
echo "[$(date +%m-%d\ %H:%M)] DRIVER C DONE" >> $BASE/insitu_driver.log

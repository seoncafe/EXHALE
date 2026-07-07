#!/bin/bash
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/lart_runs
EXEC=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/lart_runs/LaRT_build_dab6f28/LaRT_calcP.x
run() {
  cd $BASE/$1
  echo "[$(date +%m-%d\ %H:%M)] START $3" >> $BASE/insitu_driver.log
  mpirun -n 72 $EXEC $2 > lart.log 2>&1
  echo "[$(date +%m-%d\ %H:%M)] DONE $3 (exit $?)" >> $BASE/insitu_driver.log
}
run bm_hd209     hd209_lya.in    "stellar hd209 (C2, critical)"
run insitu_hd209 hd209_insitu.in "insitu hd209 (C2, 5e4)"
echo "[$(date +%m-%d\ %H:%M)] DRIVER C2 DONE" >> $BASE/insitu_driver.log

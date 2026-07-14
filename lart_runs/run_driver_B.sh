#!/bin/bash
# Driver B: hd189/hd209 stellar reruns (new winds) + their in-situ runs.
source /data/opt/oneapi_2025.3.1/setvars.sh >/dev/null 2>&1 || true
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
run() { # dir infile tag
  cd $BASE/$1
  echo "[$(date +%m-%d\ %H:%M)] START $3" >> $BASE/insitu_driver.log
  mpirun -n 72 /nfs/lart4/kiseon/LaRT/combine/LaRT_v2.00/LaRT_calcP.x $2 > lart.log 2>&1
  echo "[$(date +%m-%d\ %H:%M)] DONE $3" >> $BASE/insitu_driver.log
}
run bm_hd189     hd189_lya.in     "stellar hd189"
run insitu_hd189 hd189_insitu.in  "insitu hd189"
run insitu_hd209 hd209_insitu.in  "insitu hd209"
run bm_hd209     hd209_lya.in     "stellar hd209"
echo "[$(date +%m-%d\ %H:%M)] DRIVER B DONE" >> $BASE/insitu_driver.log

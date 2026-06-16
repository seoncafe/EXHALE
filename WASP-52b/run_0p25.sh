#!/bin/bash
EXE=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/EXHALE.x
export OMP_NUM_THREADS=8
for d in fxuv0p25_solar fxuv0p25_he98; do
  echo "===== $(date '+%H:%M:%S') START $d ====="
  ( cd "$d" && "$EXE" > run.log 2>&1 )
  echo "===== $(date '+%H:%M:%S') END $d rc=$? ====="
  grep -E "Log10 of steady-state Mdot" "$d/run.log" | tail -1
done
echo "0p25 DONE"

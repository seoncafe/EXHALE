#!/bin/bash
# Run all 6 WASP-52b cases to convergence, sequentially.
EXE=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/EXHALE.x
export OMP_NUM_THREADS=8
for d in fxuv1p0_solar fxuv0p5_solar fxuv1p0_he98 fxuv0p5_he98 fxuv1p0_solar_met fxuv0p5_solar_met; do
  echo "===== $(date '+%H:%M:%S') START $d ====="
  ( cd "$d" && "$EXE" > run.log 2>&1 )
  echo "===== $(date '+%H:%M:%S') END   $d  rc=$? ====="
  tail -n 2 "$d/run.log"
done
echo "ALL DONE"

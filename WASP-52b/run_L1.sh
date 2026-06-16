#!/bin/bash
EXE=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/EXHALE.x
export OMP_NUM_THREADS=8
for d in fxuv0p25_solar_L1 fxuv0p25_he98_L1 fxuv0p5_solar_L1 fxuv0p5_he98_L1 fxuv1p0_solar_L1 fxuv1p0_he98_L1 fxuv0p5_solar_met_L1 fxuv1p0_solar_met_L1; do
  echo "===== $(date '+%H:%M:%S') START $d ====="
  ( cd "$d" && "$EXE" > run.log 2>&1 )
  echo "===== $(date '+%H:%M:%S') END $d rc=$? ====="
  grep -E "Log10 of steady-state Mdot" "$d/run.log" | awk 'END{print "  logMdot="$NF}'
done
echo "L1 ALL DONE"

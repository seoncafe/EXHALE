#!/bin/bash
EXE=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/EXHALE.x
export OMP_NUM_THREADS=8
for d in hescan_heh005_L1 hescan_heh0010101_L1 hescan_heh0005025_L1; do
  echo "===== $(date '+%H:%M:%S') START $d ====="
  ( cd "$d" && "$EXE" > run.log 2>&1 )
  echo "===== $(date '+%H:%M:%S') END $d rc=$? ====="
done
echo "HESCAN DONE"

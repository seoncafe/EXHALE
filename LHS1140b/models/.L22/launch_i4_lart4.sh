#!/bin/bash
# The T3 solves of increment I4 on lart4 (the idle host; user instruction of
# 2026-09-17). Each reloads the state the campaign left and runs the
# alternation with the handover of "Coupled carrier solve: On stall" armed.
L=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L22
X=$L/EXHALE_L22e.x
for d in i4_wm0083 i4_wm055 i4_kz0083; do
  cd $L/$d
  nohup env OMP_NUM_THREADS=8 OMP_WAIT_POLICY=passive EXHALE_PTC_DTAU0=1.0e8 \
       EXHALE_OUTER_PASSES=40 EXHALE_JUDGED_ROWS=1 $X > run.log 2>&1 &
  echo "$! $L/$d"
done

#!/bin/bash
# Continue the He/H = 2.13 XUV chain once the 0.05 rung has written its verdict.
M=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models
until grep -q 'the stationary solve returned info' $M/.L14/x005_HeH2.13/run.log 2>/dev/null; do sleep 60; done
sleep 10
if grep -q 'certified=T' $M/.L14/x005_HeH2.13/output/Hydro_ioniz.txt 2>/dev/null; then
   exec $M/.L14/ladder.sh 2.13 $M/.L14/x005_HeH2.13/output x003 x002 x0015 x001
fi
echo "x005_HeH2.13 did not certify; chain not continued"

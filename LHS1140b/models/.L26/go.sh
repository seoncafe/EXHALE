#!/bin/bash
M=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models
cd $M/.L26/$1
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=40 \
   $M/.L26/EXHALE_L26r.x > run.log 2>&1
echo "exit $? $1" >> $M/.L26/done.txt

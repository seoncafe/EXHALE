#!/bin/bash
# The L22 step 3 solves of increments I3 and I4. Each reloads a state the
# campaign wrote and solves it on the route its input names.
L=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L22
X=$L/EXHALE_L22e.x
go() {  # go <dir> <dtau0> [extra env]
  d=$L/$1; shift; dt=$1; shift
  cd $d
  nohup env OMP_NUM_THREADS=8 OMP_WAIT_POLICY=passive EXHALE_PTC_DTAU0=$dt EXHALE_OUTER_PASSES=40 "$@" \
        $X > run.log 2>&1 &
  echo "$! $d"
}
go i3_block  1.0
go i3_alt    1.0
go i4_wm0083 1.0e8 EXHALE_JUDGED_ROWS=1
go i4_wm055  1.0e8 EXHALE_JUDGED_ROWS=1
go i4_kz0083 1.0e8 EXHALE_JUDGED_ROWS=1

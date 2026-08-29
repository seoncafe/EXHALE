#!/bin/bash
# Wait for the fixed-reservoir ladder to deliver the He/H = 11.1 arm at one
# K_zz, then start that decade's reservoir search.  The decades are
# independent once their first arm exists, so they run side by side.
#
# The readiness test is the arm's own product -- a converged closure history
# and a transit spectrum.  The only reason to give up is that the ladder is
# no longer running and the arm still is not there; log strings are not used,
# because run_arm.sh retries a refused arm and the refusal of the first
# attempt stays in the log the retry appends to.
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
ktag=$1; kval=$2; nmax=$3
d=$SC/kzz$ktag/heh11p1
while true; do
  if [ -f $d/closure_history.txt ]; then
    last=$(ls -d $d/k[0-9][0-9] 2>/dev/null | sort | tail -n 1)
    if [ -n "$last" ] && [ -f $last/tpm_He10830.txt ] &&
       grep -q "CONVERGED" $SC/kzz${ktag}_heh11p1.log 2>/dev/null; then
      break
    fi
  fi
  if ! ps -eo args --no-headers | grep -q "[l]adder_kzz.sh 11p1"; then
    echo "the ladder is no longer running and $kval has no converged arm"
    exit 1
  fi
  sleep 20
done
exec $SC/search_heh.sh $ktag $kval $d $nmax

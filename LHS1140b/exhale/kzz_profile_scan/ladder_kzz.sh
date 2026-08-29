#!/bin/bash
# Step 1 of the scan: hold the reservoir He/H fixed and walk K_zz down, each
# arm seeded from the converged wind of the arm above it.  This gives EW(K_zz)
# and the match-surface He/H at a fixed reservoir before any reservoir search
# starts.  Half-decade rungs are available because a full decade below 1e8
# moves the wind's helium profile more than the diffusion outer loop absorbs
# in one restart (kzz1e7/heh11p1_attempt_decade_step).
#
#   ladder_kzz.sh <heh_label> <heh_value> <seed_output_dir> <phi_H> <phi_He> <K_zz> ...
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
hlab=$1; hval=$2; seed=$3; phiH=$4; phiHe=$5; shift 5
for kval in "$@"; do
  ktag=$(echo "$kval" | sed 's/^1\.0e/1e/; s/\./p/g')
  echo "=== K_zz $kval (tag $ktag), He/H $hval, seed $seed ==="
  $SC/run_arm.sh $ktag $kval $hlab $hval $seed $phiH $phiHe || { echo "STOP at $kval"; exit 1; }
  d=$SC/kzz$ktag/heh$hlab
  last=$(ls -d $d/k[0-9][0-9] | sort | tail -n 1)
  read phiH phiHe <<< $(awk '!/^#/{a=$4; b=$5} END{print a, b}' $d/closure_history.txt)
  seed=$last/output
done

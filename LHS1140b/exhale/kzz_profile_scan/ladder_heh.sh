#!/bin/bash
# A plain reservoir ladder at one K_zz: each arm seeded from the arm before
# it, no target.  Used to resolve where the wind changes branch, which a
# secant search steps straight over.
#
#   ladder_heh.sh <kzz_tag> <kzz_value> <start_case_dir> <He/H> ...
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
ktag=$1; kval=$2; cur=$3; shift 3
read h ew seed fH fHe <<< $(python3 $SC/arm_state.py $cur)
echo "start: He/H=$h EW=$ew"
# Each argument is <He/H> or <He/H>:<label>.  The label is needed when the
# same reservoir ratio is solved twice at one K_zz -- once on each branch of
# the wind -- so the two arms do not overwrite each other.
for spec in "$@"; do
  nh=${spec%%:*}
  lab=${spec#*:}
  [ "$lab" = "$spec" ] && lab=$(echo $nh | tr '.' 'p')
  $SC/run_arm.sh $ktag $kval $lab $nh $seed $fH $fHe || { echo "STOP at $nh"; exit 1; }
  read h ew seed fH fHe <<< $(python3 $SC/arm_state.py $SC/kzz$ktag/heh$lab)
  echo "arm: He/H=$h EW=$ew"
done

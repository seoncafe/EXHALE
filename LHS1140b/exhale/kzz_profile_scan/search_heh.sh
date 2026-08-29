#!/bin/bash
# Step 2 of the scan: at one K_zz, walk the reservoir He/H until the red-pair
# equivalent width brackets the observed 1.108 %A.  The step is a secant in
# log EW vs log(He/H) and is capped at a factor 1.5 per arm, because a larger
# jump makes the new wind inherit the outer state of a distant seed
# (Update_EXHALE.md section 83).  Every arm is seeded from the arm before it.
#
#   search_heh.sh <kzz_tag> <kzz_value> <start_case_dir> <max_new_arms>
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
ktag=$1; kval=$2; cur=$3; nmax=$4
EW_OBS=1.108

prev_h=""; prev_ew=""
for i in $(seq 1 $nmax); do
  read h ew seed fH fHe <<< $(python3 $SC/arm_state.py $cur)
  echo "arm: He/H=$h  EW=$ew  ($cur)"
  # bracketed?  stop as soon as the two most recent arms straddle the target
  if [ -n "$prev_ew" ]; then
    lo=$(python3 -c "print(1 if ($prev_ew-$EW_OBS)*($ew-$EW_OBS) <= 0 else 0)")
    [ "$lo" = "1" ] && { echo "BRACKETED: $prev_h ($prev_ew) and $h ($ew)"; break; }
  fi
  # secant slope in log-log; the first step uses the ladder's 0.38
  nh=$(python3 - "$h" "$ew" "${prev_h:-}" "${prev_ew:-}" <<'PY'
import math, sys
h, ew = float(sys.argv[1]), float(sys.argv[2])
ph = sys.argv[3]; pe = sys.argv[4]
s = 0.38
if ph and pe and abs(math.log(float(ph)/h)) > 1e-9:
    s = math.log(float(pe)/ew)/math.log(float(ph)/h)
    s = min(max(s, 0.15), 1.2)
f = math.exp(math.log(1.108/ew)/s)
# The response is not single-valued below K_zz ~ 2e7: a 1.42x jump at
# 1.78e7 landed on the high-Mdot branch (kzz1p78e7/heh15p383, log10
# Mdot 7.83 against 7.42).  The step is capped at 1.25 so the ladder
# stays on the branch it started from.
f = min(max(f, 1.0/1.25), 1.25)
print('%.3f' % (h*f))
PY
)
  echo "next: He/H=$nh"
  prev_h=$h; prev_ew=$ew
  # clima's surface-temperature root solve fails on scattered values of the
  # reservoir ratio and succeeds a per cent away (measured in clima_probe:
  # 15.6 and 16.0 fail, 18, 20 and 30 solve), so a refused point is retried
  # beside itself rather than abandoned.  Which value actually ran is the
  # arm's directory name.
  ok=0
  for f in 1.000 1.010 0.990 1.020; do
    nh2=$(python3 -c "print('%.3f' % ($nh*$f))")
    lab=$(echo $nh2 | tr '.' 'p')
    if $SC/run_arm.sh $ktag $kval $lab $nh2 $seed $fH $fHe; then ok=1; break; fi
    echo "  He/H=$nh2 did not solve; trying beside it"
  done
  [ $ok -eq 1 ] || { echo "STOP: no reservoir near $nh solved"; exit 1; }
  cur=$SC/kzz$ktag/heh$lab
done

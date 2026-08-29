#!/bin/bash
# Reservoir He/H ladder under elemental-flux closure (Phase F item 2).
# Each point starts from the converged flux and converged wind of the one
# before it, so the loop begins near its own fixed point.
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
HERE=$EX/LHS1140b/exhale/flux_closure
DRV=$EX/src/utils/element_flux_closure.py

seed_dir=$1; shift
phi_H=$1; shift
phi_He=$1; shift

for h in "$@"; do
  d=$HERE/heh$h
  echo "=== He/H = $h : seed $seed_dir, Phi_H=$phi_H Phi_He=$phi_He ==="
  resume=""
  [ -f $d/closure_history.txt ] && resume="--resume"
  python3 $DRV $d --phi0-H $phi_H --phi0-He $phi_He \
      --config $d/closure.json --seed $seed_dir $resume \
      --tol 0.05 --kmax 8 > $HERE/heh$h.log 2>&1
  rc=$?
  echo "closure rc=$rc"
  last=$(ls -d $d/k[0-9][0-9] 2>/dev/null | sort | tail -n 1)
  if [ $rc -ne 0 ] || [ -z "$last" ]; then echo "STOP at He/H=$h"; exit 1; fi
  # transit spectrum on the converged iterate
  ( cd $last && MPLBACKEND=Agg python3 $EX/EXHALE_transit.py > transit.log 2>&1 )
  echo "transit rc=$? in $last"
  read phi_H phi_He <<< $(awk '!/^#/{fh=$4; fhe=$5} END{print fh, fhe}' $d/closure_history.txt)
  seed_dir=$last/output
done

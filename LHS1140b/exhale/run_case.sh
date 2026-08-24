#!/bin/bash
# One He/H case of the LHS 1140b scan: converge the wind, then synthesize the
# He 10830 transit spectrum in the same directory.  A per-case lock keeps a
# relaunch from competing with a run already in flight.
set -u
c=$1
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
exec 9>"/tmp/exhale_lhs_${c}.lock"
flock -n 9 || { echo "[$c] locked; exiting"; exit 0; }

echo "[$c] WIND START $(date +%H:%M:%S)"
( cd "$c" && "$EX/EXHALE.x" > run.log 2>&1 ); rc=$?
echo "[$c] WIND DONE rc=$rc $(date +%H:%M:%S)  $(grep -i 'steady-state Mdot' "$c/run.log" | tail -n 1)"
[ $rc -eq 0 ] || exit $rc
echo "[$c] TRANSIT START $(date +%H:%M:%S)"
( cd "$c" && MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1 )
echo "[$c] TRANSIT DONE rc=$? $(date +%H:%M:%S)"
grep -i 'He 10830 metrics' "$c/transit.log" | tail -n 1
echo "[$c] CASE DONE $(date +%H:%M:%S)"

#!/bin/bash
# Converge one He/H case and synthesize its transit spectrum.
#
# One case per invocation, so the four ratios run side by side.  A per-case
# lock keeps a relaunch from competing with a run already in flight, and a
# wind already running in the case directory (started by hand or by an
# earlier driver) is waited on rather than duplicated.
#
# usage: ./run_case.sh heh_10
set -u
c=$1
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

exec 9>"/tmp/exhale_${c}.lock"
if ! flock -n 9; then echo "[$c] another instance holds the lock; exiting"; exit 0; fi

# An EXHALE.x already working in this directory?
here=$(cd "$c" && pwd -P)
running=
for p in $(pgrep -x EXHALE.x); do
  [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$here" ] && running=$p
done

if [ -n "$running" ]; then
  echo "[$c] WIND ALREADY RUNNING pid=$running -- waiting  $(date +%H:%M:%S)"
  while kill -0 "$running" 2>/dev/null; do sleep 20; done
  rc=0
else
  echo "[$c] WIND START $(date +%H:%M:%S)"
  ( cd "$c" && "$EX/EXHALE.x" > run.log 2>&1 )
  rc=$?
fi
echo "[$c] WIND DONE rc=$rc $(date +%H:%M:%S)  $(grep -i 'steady-state Mdot' "$c/run.log" | tail -n 1)"

if [ $rc -eq 0 ]; then
  echo "[$c] TRANSIT START $(date +%H:%M:%S)"
  ( cd "$c" && MPLBACKEND=Agg PYTHONPATH="$EX" \
      python3 "$EX/EXHALE_transit.py" > transit.log 2>&1 )
  echo "[$c] TRANSIT DONE rc=$? $(date +%H:%M:%S)"
fi
echo "[$c] CASE DONE $(date +%H:%M:%S)"

#!/bin/bash
# THE COMPOSITION IS AN ELIMINATED VARIABLE, SO THE RESIDUAL DOES NOT DEPEND
# ON THE COMPOSITION THE CALLER HANDS THE ELIMINATION.
#
# eval_residual is given a state and a seed composition. It runs the
# equilibrium at THAT state until the composition stops moving, so a seed
# that is not the equilibrium of the state must be equilibrated first and the
# residual returned must be the residual at the equilibrated state. This runs
# the solver's own probe (EXHALE_RESID_SC_PROBE=1), which evaluates F at the
# hand-off state twice: once from the run's own composition and once from a
# composition displaced by five cells, and reports the difference through the
# acceptance measure resid_relnorm, i.e. in the units `Resid tol` is written
# in.
#
#   GREEN: the elimination is converged, so the difference is the accuracy of
#          the composition fixed point and is far below `Resid tol`.
#   RED:   with EXHALE_RESID_SC_MAX=1 the elimination is one pass and the
#          residual still carries its seed; the difference is then of the
#          order of the residual itself.
#
# The RED is the same binary with EXHALE_RESID_SC_MAX=1 exported, which caps
# the elimination at one pass and is the lagged residual this item removed;
# the variable reaches the run through the environment.
#
# Usage: ./run_seed_independence.sh <binary> <case-directory> [tolerance]
# Default tolerance 1.0e-7, one hundredth of the 1e-5 `Resid tol` the steady
# solves are accepted on.
set -u
BIN="${1:?give the EXHALE binary}"
CASE="${2:?give a case directory whose run reaches the steady solver}"
TOL="${3:-1.0e-7}"
WORK="$(mktemp -d /tmp/exhale_seedindep.XXXX)"
trap 'rm -rf "$WORK"' EXIT

for f in "$CASE"/*.inp "$CASE"/*.dat; do [ -e "$f" ] && cp -f "$f" "$WORK/"; done
mkdir -p "$WORK/output"
for f in "$CASE/output"/*_IC.txt; do [ -e "$f" ] && cp -f "$f" "$WORK/output/"; done
cp -f "$BIN" "$WORK/EXHALE.x"

# The probe runs at the ENTRY of the steady solve and stops the run there, so
# the marching before it is only the way to a real hand-off state and the case
# is used exactly as it stands. MEASURED on wasp_full_newton: 4679 marching
# steps, about six minutes single-threaded, of which the probe itself is two
# residual evaluations.

( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_RESID_SC_PROBE=1 \
      ./EXHALE.x ) > "$WORK/probe.log" 2>&1 || true

grep 'resid_sc_probe' "$WORK/probe.log"
d=$(grep 'seed dependence of the residual' "$WORK/probe.log" | tail -n 1 \
    | sed -n 's/.*Resid tol: *\([^ ]*\).*/\1/p')
if [ -z "$d" ]; then
   echo "FAIL residual_is_seed_independent measured=probe_did_not_run reference=0 tol=$TOL"
   exit 1
fi
ok=$(awk -v d="$d" -v t="$TOL" 'BEGIN{if(d<0)d=-d; print (d<=t+0)?1:0}')
if [ "$ok" = "1" ]; then
   echo "PASS residual_is_seed_independent measured=$d reference=0 tol=$TOL"
   exit 0
fi
echo "FAIL residual_is_seed_independent measured=$d reference=0 tol=$TOL"
echo "     the residual of one state changed when the elimination was given a"
echo "     different starting composition, so the composition is not being"
echo "     eliminated at the state being evaluated (eval_residual)."
exit 1

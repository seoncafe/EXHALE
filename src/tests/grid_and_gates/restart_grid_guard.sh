#!/bin/bash
# A restart onto a different grid is refused, or else keeps the run's grid.
#
# QUANTITY UNDER TEST
#   The cell centers r(j) of a run started with "Load IC? True" from state
#   files written on a DIFFERENT grid.  A state file carries one radius per
#   row while the faces r_edg, the widths dr_j and the window indices j_min
#   and j_flux come from the run's own "Grid type:", so a state whose centers
#   are not this run's cannot be marched: the file's centers and the run's
#   faces are two different constructions of one grid.  load_IC compares the
#   centers of the physical cells 1..N against the ones define_grid built,
#   with a tolerance of 1e-10 relative (the writer emits 17 significant
#   digits, and the smallest real grid difference in the tree is 3.3e-4), and
#   stops the run when they disagree.
#
# WHAT IS RUN (three short runs, none of them in the regression directory)
#   A  backup/regression/roundtrip as shipped (Grid type: Mixed, 40 steps):
#      its Hydro_ioniz.txt and Ion_species.txt are the state to restart from,
#      exactly as roundtrip_check.sh builds them.
#   B  the same input with "Grid type: Uniform" and no restart, 1 step: its
#      r column IS the run's own grid, built by the production define_grid.
#   C  "Grid type: Uniform" and "Load IC? True", with the stage-A files handed
#      back as output/*_IC.txt, 1 step.
#
# ASSERTION
#   Either C exits nonzero (the mismatch is refused), or C's r column equals
#   B's -- the run's own grid -- to 1e-12 relative.  Adopting stage A's Mixed
#   radii is the failure.  Both distances are printed either way.
#
# MEASURED 2026-09-10: the first branch is the one that fires.  C stops in
# load_IC with "holds cell centers of a different grid from the one this run
# built" and exit 1, so the row passes as restart_grid_mismatch_refused and
# the r-column comparison below is not reached.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" restart_grid_guard
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
CASE="$ROOT/backup/regression/roundtrip"

if [ ! -x "$EXE" ]; then
   echo "FAIL restart_grid_guard measured=no_binary reference=$EXE tol=0"
   exit 1
fi

for d in rgA rgB rgC; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/$d/"
done
sed -i 's/^Grid type:.*/Grid type: Uniform/' "$WORK/rgB/input.inp" \
        "$WORK/rgC/input.inp"
sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/rgC/input.inp"

( cd "$WORK/rgA" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40 "$EXE" \
     > run.log 2>&1 ) || {
   echo "FAIL restart_grid_guard_stageA measured=nonzero_exit reference=exit_0 tol=0"
   exit 1; }
( cd "$WORK/rgB" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
     > run.log 2>&1 ) || {
   echo "FAIL restart_grid_guard_stageB measured=nonzero_exit reference=exit_0 tol=0"
   exit 1; }

cp "$WORK/rgA/output/Hydro_ioniz.txt" "$WORK/rgC/output/Hydro_ioniz_IC.txt"
cp "$WORK/rgA/output/Ion_species.txt" "$WORK/rgC/output/Ion_species_IC.txt"

( cd "$WORK/rgC" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
     > run.log 2>&1 )
rcC=$?
echo "  stage C (Uniform grid, IC files written on the Mixed grid): exit $rcC"
# Exit status 2 = the run declared a stationary state that the A2 certification
# refused (docs/a2_certification_contract_20260906.md section 8); the run wrote
# its outputs in full, which is what this gate reads. Only 1 (a Fortran error
# stop) or a signal is a failed run here.
if [ $rcC -ne 0 ] && [ $rcC -ne 2 ]; then
   echo "PASS restart_grid_mismatch_refused measured=exit_$rcC reference=nonzero tol=0"
   grep -i -m 3 -E 'error|refus|mismatch' "$WORK/rgC/run.log" || true
   exit 0
fi
if [ ! -s "$WORK/rgC/output/Hydro_ioniz.txt" ]; then
   echo "FAIL restart_grid_guard_stageC measured=no_output reference=written tol=0"
   exit 1
fi

python3 - "$WORK" <<'PY'
import sys
import numpy as np

work = sys.argv[1]
rA = np.loadtxt(work + "/rgA/output/Hydro_ioniz.txt", comments="#")[:, 0]
rB = np.loadtxt(work + "/rgB/output/Hydro_ioniz.txt", comments="#")[:, 0]
rC = np.loadtxt(work + "/rgC/output/Hydro_ioniz.txt", comments="#")[:, 0]

d_own = float(np.max(np.abs(rC - rB) / rB))
d_file = float(np.max(np.abs(rC - rA) / rA))
print("  max relative distance of the restart's r column")
print("     from the run's own Uniform grid : %.6e" % d_own)
print("     from the Mixed grid of the file : %.6e" % d_file)

tol = 1.0e-12
ok = d_own <= tol
print("%s restart_keeps_own_grid measured=%.6e reference=%.6e tol=%.2e"
      % ("PASS" if ok else "FAIL", d_own, 0.0, tol))
sys.exit(0 if ok else 1)
PY
exit $?

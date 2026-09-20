#!/bin/bash
# A state keeps loading on its own base grid, and a changed default refuses it
# visibly.
#
# QUANTITY UNDER TEST
#   Whether load_IC accepts a state pair, and what it prints when it does not.
#   load_IC compares the cell centers of the state with those the run built
#   and refuses a difference above 1e-10 relative; the base cell width of the
#   Mixed grid is what moves them.
#
# WHAT IS RUN (copies of backup/regression/roundtrip, Grid type: Mixed)
#   W  the input with the pinned line "Base grid [dr,cells]:
#      1.9999999494757503e-4 50", the width every input without the key was
#      run on before 2026-09-19, marched 40 steps: its output pair is a state
#      written on that grid (base_grid_pinned_width.f90 asserts that the
#      pinned width builds the pre-2026-09-19 default grid to the bit).
#   L  the same pinned input with "Load IC? True" and W's pair as the IC,
#      1 step.
#   D  the same input WITHOUT any "Base grid" line (the present default,
#      2.0e-4) with "Load IC? True" and W's pair as the IC, 1 step.
#
# ASSERTIONS
#   1. L loads the state: it exits 0 or 2 (2 is a stationary claim the
#      certification refused, which still writes every output) and its log
#      carries no grid refusal.
#   2. D is refused: it exits 1 with load_IC's "holds cell centers of a
#      different grid" message.
#   3. The refusal names the run's base cell width at round-trip precision
#      (2.00000000000000010E-004, the default, "the default (no "Base grid"
#      key)") and the line that restores the old grid, so the contract is
#      visible where it acts rather than inferred from a tolerance.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" base_grid_restart_contract
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}/base_grid_restart"
CASE="$ROOT/backup/regression/roundtrip"
PIN='Base grid [dr,cells]: 1.9999999494757503e-4 50'

nfail=0
say_pass() { echo "PASS $1 measured=$2 reference=$3 tol=$4"; }
say_fail() { echo "FAIL $1 measured=$2 reference=$3 tol=$4"; nfail=$((nfail+1)); }

if [ ! -x "$EXE" ]; then
   say_fail base_grid_restart_binary no_binary "$EXE" 0
   exit 1
fi

rm -rf "$WORK"
for d in W L D; do
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/$d/"
   sed -i '/^Base grid/d' "$WORK/$d/input.inp"
done
for d in W L; do
   sed -i "s/^Grid type:.*/&\n$PIN/" "$WORK/$d/input.inp"
done
for d in L D; do
   sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/$d/input.inp"
done

( cd "$WORK/W" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40 "$EXE" > run.log 2>&1 )
rc=$?
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   say_fail base_grid_restart_write "exit_$rc" exit_0 0
   exit 1
fi
for d in L D; do
   cp "$WORK/W/output/Hydro_ioniz.txt" "$WORK/$d/output/Hydro_ioniz_IC.txt"
   cp "$WORK/W/output/Ion_species.txt" "$WORK/$d/output/Ion_species_IC.txt"
done

( cd "$WORK/L" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" > run.log 2>&1 )
rcL=$?
( cd "$WORK/D" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" > run.log 2>&1 )
rcD=$?
echo "  L (pinned input, pinned state): exit $rcL"
echo "  D (default input, pinned state): exit $rcD"

if { [ $rcL -eq 0 ] || [ $rcL -eq 2 ]; } && \
   ! grep -q 'holds cell centers of a different grid' "$WORK/L/run.log"; then
   say_pass base_grid_pinned_state_loads_on_pinned_input "exit_$rcL" loaded 0
else
   say_fail base_grid_pinned_state_loads_on_pinned_input "exit_$rcL" loaded 0
   grep -m 3 -E 'ERROR|refus' "$WORK/L/run.log"
fi

if [ $rcD -eq 1 ] && \
   grep -q 'holds cell centers of a different grid' "$WORK/D/run.log"; then
   say_pass base_grid_old_state_refused_on_default "exit_$rcD" refused 0
else
   say_fail base_grid_old_state_refused_on_default "exit_$rcD" refused 0
fi

if grep -q 'base grid: width *2.00000000000000010E-004 R_p x 50 uniform cells, the default' \
        "$WORK/D/run.log" && \
   grep -qF '"Base grid [dr,cells]: 1.9999999494757503e-4 50"' "$WORK/D/run.log"; then
   say_pass base_grid_refusal_names_the_width stated stated 0
else
   say_fail base_grid_refusal_names_the_width absent stated 0
fi
grep -A 9 'holds cell centers of a different grid' "$WORK/D/run.log" \
   | sed 's/^/  /'

if [ $nfail -gt 0 ]; then exit 1; fi
exit 0

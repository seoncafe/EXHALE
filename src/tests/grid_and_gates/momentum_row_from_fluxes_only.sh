#!/bin/bash
# The momentum row of the conserved state is written by the operator split
# and by nothing else.
#
# QUANTITY UNDER TEST
#   u(2,j) = rho v, the momentum row of the conserved state, in code units.
#   The finite-volume update moves it by the Riemann face fluxes, the
#   geometric and gravitational sources, the semi-implicit energy step, the
#   transport and diffusion operators and Apply_BC. Anything else that
#   assigns to u adds momentum that no flux and no source produced, and that
#   increment is then carried into the next step's hydro, into the steady
#   residual, and into the mass flux rho v r^2 the convergence test measures.
#
#   src/EXHALE_main.f90 carried, between the du evaluation and the state
#   change norm,
#       u(2,:) = u(2,:) + 1.0e-16 ! To avoid division by zero
#   a constant increment applied to every cell and every ghost of the
#   momentum row on every step, whose stated purpose was to keep the NEXT
#   step's dtu denominator away from zero. The denominator belongs to the
#   diagnostic, so the guard belongs at the point of use and the state stays
#   the state.
#
# ASSERTION 1 (source): no statement adds a constant to a whole row of the
#   conserved state. measured = the number of such statements found in
#   src/EXHALE_main.f90 outside comment lines, reference = 0, tol = 0.
#   This is a source assertion and not a run measurement because a
#   1.0e-16 increment of the code momentum unit is absorbed into a state
#   that stays self-consistent (u, rho and v are re-derived from each other
#   every step), so no single run's output files expose it; what a pair of
#   runs shows is a changed state, which is a movement measurement and not
#   a criterion.
#
# ASSERTION 2 (run): the relative state change dtu the run reports is a
#   finite number at every step. This is the guard that replaces the
#   increment: where a conserved row was exactly zero at the start of a
#   step the ratio 1 - u/u_old is not a number, and the run must report an
#   unbounded change rather than a NaN that no comparison can act on.
#   measured = the number of steps whose dtu is NaN or Infinity in a
#   200-step run, reference = 0, tol = 0.
#
# WHAT IS RUN
#   A copy of backup/regression/mol_base_handoff (its input.inp and
#   base.inp) in build/tests/grid_and_gates/momentum_row, capped at 200
#   steps. The regression directory is never written to.
#
# EXPECTED AT HEAD 35d9dd5: assertion 1 RED, assertion 2 GREEN.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
SRC="$ROOT/src/EXHALE_main.f90"
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}/momentum_row"
CASE="$ROOT/backup/regression/mol_base_handoff"
n_fail=0

# ---- assertion 1: the source of the marching loop ----
# Comments are stripped, then a whole-row assignment of the form
#     u(<row>,:) = u(<row>,:) +/- <numeric literal>
# is matched: a CONSTANT increment, with no state, no time step and no rate
# in it. The explicit energy source u(3,:) = u(3,:) + dt_loc*(heat - cool)
# is an operator of the split and carries dt_loc, so it does not match.
PATTERN='u\([0-9]+,:\)[[:space:]]*=[[:space:]]*u\([0-9]+,:\)[[:space:]]*[-+][[:space:]]*[0-9][0-9.]*([deDE][-+]?[0-9]+)?[[:space:]]*$'
n_add=$(sed -e 's/!.*$//' -e 's/[[:space:]]*$//' "$SRC" \
        | grep -cE "$PATTERN" || true)
if [ "$n_add" -ne 0 ]; then
   echo "  the statements found:"
   sed -e 's/!.*$//' -e 's/[[:space:]]*$//' "$SRC" \
      | grep -nE "$PATTERN" | sed 's/^/     /'
fi
if [ "$n_add" -eq 0 ]; then
   echo "PASS momentum_row_has_no_added_constant measured=$n_add reference=0 tol=0"
else
   echo "FAIL momentum_row_has_no_added_constant measured=$n_add reference=0 tol=0"
   n_fail=$((n_fail+1))
fi

# ---- assertion 2: the state change norm the run reports ----
if [ ! -x "$EXE" ]; then
   echo "FAIL state_change_norm_is_finite measured=no_binary reference=$EXE tol=0"
   exit 1
fi
rm -rf "$WORK"
mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/"
( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=200 "$EXE" > run.log 2>&1 )
rc=$?
# Exit status 2 = the run declared a stationary state that the A2 certification
# refused (docs/a2_certification_contract_20260906.md section 8); the run wrote
# its outputs in full, which is what this gate reads. Only 1 (a Fortran error
# stop) or a signal is a failed run here.
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL state_change_norm_is_finite measured=exit_$rc reference=exit_0 tol=0"
   echo "     see $WORK/run.log"
   exit 1
fi

# The step lines are "<count> <du> <dtu> <lev_rel>", written by the loop.
n_steps=$(awk '$1 ~ /^[0-9]+$/ && NF == 4 { n++ } END { print n+0 }' "$WORK/run.log")
n_bad=$(awk 'BEGIN{IGNORECASE=1}
             $1 ~ /^[0-9]+$/ && NF == 4 {
                if ($3 ~ /nan/ || $3 ~ /inf/) n++
             } END { print n+0 }' "$WORK/run.log")
echo "  step lines read from the run: $n_steps"
if [ "$n_bad" -eq 0 ]; then
   echo "PASS state_change_norm_is_finite measured=$n_bad reference=0 tol=0"
else
   echo "FAIL state_change_norm_is_finite measured=$n_bad reference=0 tol=0"
   n_fail=$((n_fail+1))
fi

exit $((n_fail > 0))

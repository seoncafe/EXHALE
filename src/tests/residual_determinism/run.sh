#!/bin/bash
# THE RESIDUAL CONTRACT, IN THREE PARTS (PLAN_20260909_rev1 item N5;
# PLAN_20260909_review F9; ISSUES_20260909_review 5.1).
#
# The residual of the stationary system is R(Y) = L(Y) + S(Y, c*(Y)), and
# c*(Y) is a fixed point of the composition elimination. Three statements
# about it are separate, because the physics allows one to fail while
# another holds, and a test that ran them together would hide a real
# bistability behind a tolerance:
#
#   1  REPLAY DETERMINISM. Identical complete inputs and identical branch
#      data return the same F(Y) bitwise, however many other states were
#      evaluated and discarded in between.
#   2  CONDITIONAL CLOSURE CONVERGENCE. Seeds inside one branch give the
#      same outer residual within a measured budget of the ROW MAXIMA of
#      the certification -- the row maxima and not a window sum, because
#      the window integrals of the assembled residual are not Lipschitz in
#      the state in the odd-even band above the base (N9).
#   3  BRANCH DISCOVERY. A seed that reaches a distinct root is recorded
#      with its composition and its residual, never averaged and never
#      forced onto a preferred root.
#
# Contract 1 is measured by EXHALE_RESID_DETERMINISM=1, contracts 2 and 3 by
# EXHALE_RESID_BRANCH_REPORT=1; both stop the run where they are written.
#
# CASES. Two live in the tree and are run by default:
#   mol_base_handoff     the molecular gate rung, three unknowns per cell;
#   atomic_elem_newton   the atomic element reload, eleven unknowns per
#                        cell (three hydrodynamic and eight element rows),
#                        loaded from its own IC/ directory.
# A carrier reload is not in the tree (it is a session state, and its
# restart is 600 kB), so it is given as a run directory:
#   EXHALE_RESID_RELOADS=<dir>[:<dir>...]
# Each such directory holds input.inp, any base.inp/metals.inp and its
# restart in output/*_IC.txt.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line. Exit status is nonzero if any row fails.
#
# EXHALE_EXE selects the binary (default $ROOT/EXHALE.x); EXHALE_RESID_EXE
# is its accepted alias and a conflicting pair is refused. The policy, and
# the identity block this suite prints before it runs anything, are in
# src/tests/exhale_exe.sh.
#
# Usage: src/tests/residual_determinism/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" residual_determinism
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/residual_determinism}"
mkdir -p "$OUT"
rc=0

if [ ! -x "$EXE" ]; then
   echo "FAIL residual_determinism_binary measured=missing reference=$EXE tol=0"
   exit 1
fi

CASES="$ROOT/backup/regression/mol_base_handoff $ROOT/backup/regression/atomic_elem_newton"
if [ -n "${EXHALE_RESID_RELOADS:-}" ]; then
   CASES="$CASES $(echo "$EXHALE_RESID_RELOADS" | tr ':' ' ')"
fi

# ------------------------------------------------------------------ #
# CONTRACT 1: the replay, on every case.
#
# Two rows each. The verdict row is the contract itself. The second row is
# the one that was RED before this item: the species-row registry is empty
# until a solve registers it, so the replay used to evaluate the
# three-unknown residual in EVERY configuration -- 1500 entries on the
# carrier, element and atomic-element reloads alike -- and said nothing
# about the rows those configurations solve. It now replays the system the
# configuration carries, and the row states how many unknowns per cell that
# is.
# ------------------------------------------------------------------ #
for CASE in $CASES; do
   name="$(basename "$CASE")"
   log="$OUT/replay_$name.log"
   EXHALE_EXE="$EXE" "$HERE/run_residual_determinism.sh" "$CASE" replay \
       > "$log" 2>&1 || true
   verdict="$(grep -c 'VERDICT: the residual is a state function' "$log" || true)"
   nvar="$(sed -n 's/.*system replayed: \([0-9]*\) unknowns per cell.*/\1/p' "$log" | head -n 1)"
   [ -z "$nvar" ] && nvar=0
   if [ "$verdict" -ge 1 ]; then
      echo "PASS replay_is_a_state_function_$name measured=state_function reference=state_function tol=0"
   else
      echo "FAIL replay_is_a_state_function_$name measured=history_dependent reference=state_function tol=0"
      echo "     see $log"
      rc=1
   fi
   echo "     $name replays $nvar unknowns per cell"
   # The binary the case actually executed is the copy the per-case script
   # made in its working directory, and the row it printed there names that
   # copy's md5. This is the evidence of which build produced the numbers
   # above; a passing replay row is not.
   ranmd5="$(sed -n 's/^PASS residual_determinism_binary_identity measured=\([0-9a-f]*\) .*/\1/p' \
             "$log" | head -n 1)"
   [ -z "$ranmd5" ] && ranmd5=not_reported
   if [ "$ranmd5" = "$EXHALE_RUN_EXE_MD5" ]; then
      echo "PASS the_binary_that_ran_is_the_one_requested_$name measured=$ranmd5 reference=$EXHALE_RUN_EXE_MD5 tol=0"
   else
      echo "FAIL the_binary_that_ran_is_the_one_requested_$name measured=$ranmd5 reference=$EXHALE_RUN_EXE_MD5 tol=0"
      rc=1
   fi
done

# The atomic element reload is the case whose registry the replay used to
# miss entirely: three hydrodynamic rows plus the element rows of the
# elements metals.inp carries. A replay of three unknowns per cell there is
# the old defect back.
alog="$OUT/replay_atomic_elem_newton.log"
anvar="$(sed -n 's/.*system replayed: \([0-9]*\) unknowns per cell.*/\1/p' "$alog" 2>/dev/null | head -n 1)"
[ -z "$anvar" ] && anvar=0
if [ "$anvar" -gt 3 ]; then
   echo "PASS replay_carries_the_species_rows measured=$anvar reference=>3 tol=0"
else
   echo "FAIL replay_carries_the_species_rows measured=$anvar reference=>3 tol=0"
   echo "     the replay evaluated the three-unknown residual on a"
   echo "     configuration that registers element rows"
   rc=1
fi

# ------------------------------------------------------------------ #
# CONTRACTS 2 AND 3: the seed family, on the cases that register a species
# row. A configuration with no species row has no branch data to state, so
# the family is not run there.
# ------------------------------------------------------------------ #
for CASE in $CASES; do
   name="$(basename "$CASE")"
   nvar="$(sed -n 's/.*system replayed: \([0-9]*\) unknowns per cell.*/\1/p' \
           "$OUT/replay_$name.log" 2>/dev/null | head -n 1)"
   [ -z "$nvar" ] && nvar=0
   [ "$nvar" -le 3 ] && continue
   log="$OUT/branch_$name.log"
   EXHALE_EXE="$EXE" "$HERE/run_residual_determinism.sh" "$CASE" branch \
       > "$log" 2>&1 || true
   if ! grep -q 'resid_branch) VERDICT' "$log"; then
      echo "FAIL seed_family_is_evaluated_$name measured=no_report reference=report tol=0"
      echo "     see $log"
      rc=1
      continue
   fi
   echo "PASS seed_family_is_evaluated_$name measured=report reference=report tol=0"

   # CONTRACT 2. The budget is the certification tolerance of each row
   # class: the spread of a row maximum over the SCALED seeds of the
   # reference branch -- small admissible perturbations of one eliminated
   # partition, which is the family contract 2 is about -- divided by that
   # class's tolerance, must be below one. Anything above
   # one means the closure moves the judged state by more than the
   # tolerance the state is judged with, which is a closure to tighten or a
   # slow mode to promote into Y -- it is not a tolerance to widen.
   worst="$(grep 'resid_branch) scaled seeds of root  1,' "$log" \
            | sed -n 's/.*spread over tolerance *//p' \
            | sort -g | tail -n 1)"
   [ -z "$worst" ] && worst=0
   ok="$(awk -v w="$worst" 'BEGIN{print (w+0 < 1.0) ? 1 : 0}')"
   if [ "$ok" -eq 1 ]; then
      echo "PASS closure_spread_within_the_row_tolerance_$name measured=$worst reference=<1 tol=0"
   else
      echo "FAIL closure_spread_within_the_row_tolerance_$name measured=$worst reference=<1 tol=0"
      rc=1
   fi

   # CONTRACT 3. Every root the family reached is named on its own line
   # with the cells the branch datum is read at and the composition there.
   # A report that says it found n roots and prints fewer has averaged
   # them, which is the thing this contract forbids.
   nroot="$(sed -n 's/.*resid_branch) roots found \([0-9]*\) among.*/\1/p' "$log" | head -n 1)"
   [ -z "$nroot" ] && nroot=0
   nlines="$(grep -c 'resid_branch) root  *[0-9]* holds' "$log" || true)"
   if [ "$nroot" -ge 1 ] && [ "$nlines" -eq "$nroot" ]; then
      echo "PASS every_root_is_recorded_$name measured=$nlines reference=$nroot tol=0"
   else
      echo "FAIL every_root_is_recorded_$name measured=$nlines reference=$nroot tol=0"
      rc=1
   fi
done

exit $rc

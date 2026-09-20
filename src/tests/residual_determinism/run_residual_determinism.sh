#!/bin/bash
# src/tests/residual_determinism/run_residual_determinism.sh
#
# ONE CASE OF THE RESIDUAL CONTRACT, run and printed. run.sh beside this
# file turns the printed blocks into PASS/FAIL rows; this script exists so
# that a single configuration can be looked at on its own.
#
# CONTRACT 1, REPLAY DETERMINISM: the steady residual is a function of its
# argument. Evaluating F at a state, then at other states, then at the first
# state again must return the SAME BITS. Anything else means the Jacobian
# columns and the line-search comparisons of a solve were taken against
# moving ground.
#
# What carries history is the equilibrium sweep's starting composition, and
# since section 147 the interface forbids it: eval_residual takes the seed as
# one argument and writes the swept composition to another, and the two may not
# be the same array, so a caller cannot hand the sweep its own previous answer
# even by accident. This test is what keeps that true.
#
# The hook evaluates F(Y0) twice with nothing in between (the control), then
# F(Y0) after three trials, then F(Y0) after a discarded trust-region trial
# and a Jacobian-vector probe, reusing ONE output workspace so that it
# carries each evaluation's composition into the next call. The control
# absorbs the one-time table building of the first evaluation in a run; the
# replays that follow must agree bit for bit. The system replayed is the one
# a solve of this configuration would carry, species rows included.
#
# CONTRACT 2 (conditional closure convergence) and CONTRACT 3 (branch
# discovery) are the other half, under EXHALE_RESID_BRANCH_REPORT=1:
# a family of admissible seeds on one state, the spread of each certified row
# maximum inside one branch, and every distinct root recorded.
#
# Usage:  ./run_residual_determinism.sh [case-directory] [replay|branch]
# Default case: backup/regression/mol_base_handoff with the gate rung's
# restart, i.e. a molecular state where the caloric EOS makes the ghost
# conversion composition-dependent and the invariant is hardest.
#
# EXHALE_EXE selects the binary (default $ROOT/EXHALE.x), which is what a
# concurrent item's private build needs; EXHALE_RESID_EXE is its accepted
# alias and a conflicting pair is refused. The script states the path and
# md5 of the binary before it runs it, and asserts the same md5 on the copy
# it executes. The policy is in src/tests/exhale_exe.sh.
#
# Exit status 0 if the residual is bitwise reproducible, 1 otherwise.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
CASE="${1:-$ROOT/backup/regression/mol_base_handoff}"
WHAT="${2:-replay}"
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" residual_determinism
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="$(mktemp -d /tmp/exhale_determinism.XXXX)"
trap 'rm -rf "$WORK"' EXIT

for f in "$CASE"/*.inp "$CASE"/*.dat; do [ -e "$f" ] && cp -f "$f" "$WORK/"; done
mkdir -p "$WORK/output"
# The restart of a case that keeps one: beside its outputs, or in an IC/
# directory of its own. A case whose input.inp does not already load it is
# switched to, because the state the contract is measured on is that restart
# and not a cold start.
nic=0
for f in "$CASE/output"/*_IC.txt; do
    [ -e "$f" ] && cp -f "$f" "$WORK/output/" && nic=$((nic+1))
done
if [ "$nic" -eq 0 ]; then
    for f in "$CASE/IC"/*_IC.txt; do
        [ -e "$f" ] && cp -f "$f" "$WORK/output/" && nic=$((nic+1))
    done
    if [ "$nic" -gt 0 ]; then
        sed -i -e 's/^Load IC?.*/Load IC? True/' "$WORK/input.inp"
    fi
fi
cp -f "$EXE" "$WORK/EXHALE.x"
# The copy is the file that executes, so it is the file whose identity is
# asserted.
exhale_assert_exe_identity "$WORK/EXHALE.x" || exit 1

if [ "$WHAT" = "branch" ]; then
    ( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_RESID_BRANCH_REPORT=1 \
          ./EXHALE.x ) > "$WORK/determinism.log" 2>&1 || true
    grep '(resid_branch)' "$WORK/determinism.log" || true
    grep -q 'resid_branch) VERDICT' "$WORK/determinism.log" && exit 0
    echo "==> BRANCH REPORT MISSING: the run did not reach the stationary"
    echo "    solve, so no seed family was evaluated."
    exit 1
fi

( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_RESID_DETERMINISM=1 ./EXHALE.x ) \
    > "$WORK/determinism.log" 2>&1 || true

sed -n '/resid_determinism/,/VERDICT/p' "$WORK/determinism.log"
if grep -q 'VERDICT: the residual is a state function' "$WORK/determinism.log"; then
    echo "==> DETERMINISM PASS"
    exit 0
fi
echo "==> DETERMINISM FAIL: F(Y) moved MORE when another state was evaluated"
echo "    in between than it moves when nothing is (the control line above)."
echo "    The carrier to look at first is the equilibrium sweep's starting"
echo "    composition: since section 147 eval_residual takes it as one"
echo "    argument and writes the swept composition to another, so a caller"
echo "    cannot seed the sweep from its own previous answer.  If that is"
echo "    still intact, something else is holding state between evaluations."
exit 1

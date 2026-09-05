#!/bin/bash
# src/tests/residual_determinism/run_residual_determinism.sh
#
# THE INVARIANT: the steady residual is a function of its argument.
# Evaluating F at a state, then at two other states, then at the first state
# again must return the SAME BITS. Anything else means the Jacobian columns and
# the line-search comparisons of a solve were taken against moving ground.
#
# What carries history is the equilibrium sweep's starting composition, and
# since section 147 the interface forbids it: eval_residual takes the seed as
# one argument and writes the swept composition to another, and the two may not
# be the same array, so a caller cannot hand the sweep its own previous answer
# even by accident. This test is what keeps that true.
#
# The hook evaluates F(Y0) twice with nothing in between (the control), then
# F(Y0), a trial, F(Y0), another trial, F(Y0), reusing ONE output workspace so
# that it carries each evaluation's composition into the next call. The control
# absorbs the one-time table building of the first evaluation in a run; the
# three that follow must agree bit for bit.
#
# Usage:  ./run_residual_determinism.sh [case-directory]
# Default case: backup/regression/mol_base_handoff with the gate rung's
# restart, i.e. a molecular state where the caloric EOS makes the ghost
# conversion composition-dependent and the invariant is hardest.
#
# Exit status 0 if the residual is bitwise reproducible, 1 otherwise.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
CASE="${1:-$ROOT/backup/regression/mol_base_handoff}"
WORK="$(mktemp -d /tmp/exhale_determinism.XXXX)"
trap 'rm -rf "$WORK"' EXIT

for f in "$CASE"/*.inp "$CASE"/*.dat; do [ -e "$f" ] && cp -f "$f" "$WORK/"; done
mkdir -p "$WORK/output"
for f in "$CASE/output"/*_IC.txt; do [ -e "$f" ] && cp -f "$f" "$WORK/output/"; done
cp -f "$ROOT/EXHALE.x" "$WORK/"

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

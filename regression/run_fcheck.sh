#!/bin/bash
# Periodic runtime-checked regression (code review 2026-07-02 recommendation):
# rebuild EXHALE with gfortran runtime checks in a scratch copy of build/,
# run a short bounded HD 209458 b case (He 2^3S + metals if metals.inp is
# given), and fail loudly on any bounds/uninitialized-memory trap.
#
# Usage:  ./regression/run_fcheck.sh [steps]     (default 1500)
# Exit 0 = clean; nonzero = runtime check fired or run failed.
set -u
STEPS=${1:-1500}
EX=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d /tmp/exhale_fcheck.XXXX)
trap 'rm -rf "$WORK"' EXIT

echo "== fcheck regression: rebuilding with runtime checks =="
cd "$EX"
make distclean > /dev/null
make FFLAGS='-O1 -fopenmp -g -fcheck=bounds,do,mem' \
     > "$WORK/build.log" 2>&1 || { echo "BUILD FAILED"; exit 2; }

mkdir -p "$WORK/run/output"
cp EXHALE.x "$WORK/run/"
cp HD209458b/input.inp "$WORK/run/"
[ -f HD209458b/metals.inp ] && cp HD209458b/metals.inp "$WORK/run/"
sed -i 's/Include He23S? False/Include He23S? True/' "$WORK/run/input.inp"

echo "== running $STEPS bounded steps =="
( cd "$WORK/run" && ATES_MAXSTEPS=$STEPS ./EXHALE.x > run.log 2>&1 )
RC=$?
grep -iE "Fortran runtime error|Backtrace" "$WORK/run/run.log" \
     && { echo "RUNTIME CHECK FIRED (see above)"; RC=3; }
grep -i "steady-state Mdot" "$WORK/run/run.log"

echo "== restoring production build =="
make distclean > /dev/null && make > /dev/null 2>&1 && echo "production EXHALE.x restored"

[ $RC -eq 0 ] && echo "fcheck regression: CLEAN" || echo "fcheck regression: FAILED (rc=$RC)"
exit $RC

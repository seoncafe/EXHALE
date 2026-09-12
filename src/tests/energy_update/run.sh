#!/bin/bash
# Build and run the tests of the two implicit temperature updates of the
# marching loop: the residual-controlled solve of the energy source step
# (docs/PLAN_20260906_rev2.md step B2) and the Crank-Nicolson transport
# stage that carries the second temperature floor (step B2b).
# See README.md in this directory for what each one asserts.
#
# Every test compares a measured number or a status against a reference and
# prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion. Lines that begin with two spaces or with DIAGNOSTIC are
# context, not verdicts. Exit status is nonzero if any assertion fails.
#
# The driver links the PRODUCTION objects, so the solver it tests is the one
# the binary was built from. It supplies its own analytic cooling laws in
# place of eval_cool, which is what gives each case a closed form or an
# independently integrable root.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped, because the
# objects then need not match the default build/. EXHALE_TEST_OUT selects
# where the test executable is written.
#
# Usage: src/tests/energy_update/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/energy_update}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

# LAPACK, resolved the way the Makefile resolves it: the library of the
# compiler's own prefix when that prefix carries an OpenBLAS, else -llapack.
FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL energy_update_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   # The driver links the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL energy_update_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
# Every object but the main program: the driver brings its own.
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/energy_update_tests.x" "$OUT/conduction_floor_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/energy_update_tests.x" \
    "$HERE/cooling_models.f90" "$HERE/energy_update_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL energy_update_build measured=compile_error reference=ok tol=0"
   exit 1; }

# The transport stage that carries the second temperature floor of the
# marching loop (viscous_conduction): what it returns, what it refuses to
# return, and what it conserves.
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/conduction_floor_tests.x" \
    "$HERE/conduction_floor_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL conduction_floor_build measured=compile_error reference=ok tol=0"
   exit 1; }

rc=0
env OMP_NUM_THREADS=1 "$OUT/energy_update_tests.x"  || rc=1
echo ""
env OMP_NUM_THREADS=1 "$OUT/conduction_floor_tests.x" || rc=1
echo ""
if [ $rc -ne 0 ]; then
   echo "energy_update: a test program reported a failing assertion"
   exit 1
fi
echo "energy_update: every test program passed"

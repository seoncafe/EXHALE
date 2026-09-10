#!/bin/bash
# Build and run the element-operator boundary tests
# (docs/PLAN_20260909_rev1.md item N26b).
#
# The Fortran driver links the PRODUCTION objects and states, on a synthetic
# column, that element_transport_residual poses its own upper boundary: the
# zero-gradient condition the fixed-wind relaxation writes on its working
# arrays, so that the two spellings of the element balance are one operator
# whatever the caller left in the composition ghosts.  See the header of
# element_operator_tests.f90 for what each row asserts.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Lines beginning with two spaces or with DIAGNOSTIC are context.
# Exit status is nonzero if any row fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped.
# EXHALE_TEST_OUT selects where the test executable is written.
#
# Usage: src/tests/element_operator/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/element_operator}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

# LAPACK, resolved the way the Makefile resolves it: the library of the
# compiler's own prefix when that prefix carries an OpenBLAS, else -llapack.
FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib"
else
   LAPACK="-llapack"
fi

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL element_operator_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   # The driver links the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL element_operator_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/element_operator_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/element_operator_tests.x" \
    "$HERE/element_operator_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL element_operator_build measured=compile_error reference=ok tol=0"
   exit 1; }

rc=0
# The driver runs where it stands and reads no input file.
( cd "$OUT" && OMP_NUM_THREADS=1 "$OUT/element_operator_tests.x" ) || rc=1
exit $rc

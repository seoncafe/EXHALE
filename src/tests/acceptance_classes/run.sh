#!/bin/bash
# Build and run the acceptance-class contract test: one meaning of the
# ionization acceptance classes for every consumer (classes 1, 2, 3 and 5 are
# roots of the requested equations; classes 4 and 6 are not, and a state
# carrying one of those passes no final acceptance).
#
# The Fortran driver links the PRODUCTION objects, so the counting function
# and the acceptance predicate it tests are the ones the binary was built
# from.
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a parallel build
# uses); with it set the staleness check is skipped.
#
# Usage: src/tests/acceptance_classes/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/acceptance_classes}"
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
   echo "FAIL acceptance_classes_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     build the code first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL acceptance_classes_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; build first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/acceptance_class_contract.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/acceptance_class_contract.x" \
    "$HERE/acceptance_class_contract.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL acceptance_class_contract_build measured=compile_error reference=ok tol=0"
   exit 1; }

env OMP_NUM_THREADS=1 "$OUT/acceptance_class_contract.x"

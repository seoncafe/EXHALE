#!/bin/bash
# Build and run the ionization-stage-flux identity tests
# (docs/lhs1140b_stationary_L12b_derivation_20260916.md section 4;
#  PLAN_20260916_rev3 section 8, item L12 stage B).
#
# The driver links the PRODUCTION objects and states, on synthetic columns,
# that the stage flux of the derivation memo sums over ALL stages of an
# element -- neutral and ionized -- to the element flux at every face, and
# that the stage eddy term telescopes to zero.  It tests the formula, the
# face-fraction rule, the production module's own rows and the two end rows
# the transport operator folds the ghost onto.  See the header
# of ionization_stage_flux_tests.f90 for what each row asserts and for the
# one thing the rows cannot reach (the operator's own diffusive face
# coefficients are private).
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
# Usage: src/tests/ionization_stage_flux/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/ionization_stage_flux}"
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
   echo "FAIL ionization_stage_flux_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   # The driver links the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL ionization_stage_flux_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

# src/tests/test_columns.f90 builds the synthetic column whose species
# carry their own density; it belongs to no production object, so it is
# compiled here beside the driver.
rm -f "$OUT/ionization_stage_flux_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/ionization_stage_flux_tests.x" \
    "$ROOT/src/tests/test_columns.f90" \
    "$HERE/ionization_stage_flux_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL ionization_stage_flux_build measured=compile_error reference=ok tol=0"
   exit 1; }

rc=0
# The driver runs where it stands and reads no input file.
( cd "$OUT" && OMP_NUM_THREADS=1 "$OUT/ionization_stage_flux_tests.x" ) || rc=1
exit $rc

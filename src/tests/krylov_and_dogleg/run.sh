#!/bin/bash
# THE LINEAR SOLVE AND THE STEP GEOMETRY OF THE STATIONARY SOLVER
# (docs/PLAN_20260909_rev1.md section 4, N0 to N3;
#  docs/ISSUES_20260909_review.md R1, R2, R3).
#
# The Fortran driver links the PRODUCTION objects and calls the production
# pgmres, dogleg_step, cauchy_length_along_the_banded_gradient,
# predicted_model_decrease and dogleg_step on small linear
# systems it states itself.  What each row asserts is in the header of
# krylov_and_dogleg_tests.f90.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Exit status is nonzero if any row fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped.
# EXHALE_TEST_OUT selects where the test executable is written.
# FFLAGS_TEST overrides the driver's own flags, which is how the bounds and
# memory checked build of this suite is taken:
#   FFLAGS_TEST='-O0 -g -fbacktrace -fopenmp -fcheck=bounds,do,mem' run.sh
#
# Usage: src/tests/krylov_and_dogleg/run.sh
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/krylov_and_dogleg}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib"
else
   LAPACK="-llapack"
fi

rc=0

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL krylov_and_dogleg_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL krylov_and_dogleg_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/krylov_and_dogleg_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/krylov_and_dogleg_tests.x" \
    "$HERE/krylov_and_dogleg_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL krylov_and_dogleg_build measured=compile_error reference=ok tol=0"
   exit 1; }

"$OUT/krylov_and_dogleg_tests.x" || rc=1

exit $rc

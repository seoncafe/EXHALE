#!/bin/bash
# Build and run the stage-row balance tests (PLAN_20260918_rev2 items D3a
# and D7b).
#
# The driver links the PRODUCTION objects and asks the two reaction kernels
# the stage sources come from -- mol_heh_rows in molecular gas and
# heh_tr_rows in atomic gas -- for the gross channels of the H II, He II and
# He III rows, then checks that the signed channel sums are the assembled
# source of each stage after the caller's charge exchange: the He <-> H pair
# of Table 4 group B in every cell, and on a metal-bearing cell the metal
# charge exchange of groups A, C and E as well, which the transport operator
# adds as stoichiometric stage sources.
# It also derives and validates the floating-point bound of the
# stage-nucleus-sum identity.
#
# With
#   EXHALE_STAGE_ROW_RECORD=<run>/output/stage_channels.txt
#   EXHALE_STAGE_ROW_TERMS=<run>/output/carrier_row_terms.txt
# it re-evaluates the rows at the state a run recorded
# (EXHALE_STAGE_CHANNELS=<cell list> on that run) and prints the D3a step 2
# tables for every recorded cell.
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
# Usage: src/tests/stage_row_balance/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/stage_row_balance}"
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
   echo "FAIL stage_row_balance_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   # The driver links the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL stage_row_balance_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/stage_row_balance_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/stage_row_balance_tests.x" \
    "$HERE/stage_row_balance_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL stage_row_balance_build measured=compile_error reference=ok tol=0"
   exit 1; }

rc=0
# The driver runs where it stands and reads no input file.
( cd "$OUT" && OMP_NUM_THREADS=1 "$OUT/stage_row_balance_tests.x" ) || rc=1
exit $rc

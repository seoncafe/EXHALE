#!/bin/bash
# Build and run the coupled-block handover contract test: the transition
# from "the wind-composition alternation stopped at the movement bound" to
# "the coupled carrier block is entered", and the one entry state it
# carries.
#
# The Fortran driver links the PRODUCTION objects, so the transition, the
# snapshot and the acceptance it tests are the ones the binary was built
# from (coupled_block_handover.f90).
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a parallel build
# uses); with it set the staleness check is skipped.
#
# Usage: src/tests/coupled_block_handover/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/coupled_block_handover}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fcheck=all -fbacktrace -fopenmp}"

mkdir -p "$OUT"

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL coupled_block_handover_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     build the code first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL coupled_block_handover_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; build first"
      exit 1
   fi
fi

# The module under test uses nothing else, so the object of that module is
# the whole link: no LAPACK, no MINPACK, no grid.
MOD_OBJ="$OBJDIR/coupled_block_handover.o"
if [ ! -f "$MOD_OBJ" ]; then
   echo "FAIL coupled_block_handover_build measured=no_object reference=$MOD_OBJ tol=0"
   exit 1
fi

rm -f "$OUT/coupled_block_handover_contract.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/coupled_block_handover_contract.x" \
    "$HERE/coupled_block_handover_contract.f90" "$MOD_OBJ" || {
   echo "FAIL coupled_block_handover_contract_build measured=compile_error reference=ok tol=0"
   exit 1; }

env OMP_NUM_THREADS=1 "$OUT/coupled_block_handover_contract.x"

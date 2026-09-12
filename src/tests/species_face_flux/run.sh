#!/bin/bash
# Build and run the tests of the species advective transport on the face
# mass fluxes (docs/PLAN_20260906_rev2.md step B4-1;
# docs/b4_spatial_operator_design_20260906.md T-B4.1 to T-B4.3, rows B4-b
# and B4-j).
#
# Rows.
#   * The Fortran driver links the PRODUCTION objects and states the
#     identities on a synthetic column: a uniform composition preserved, the
#     face identity, the bounds and the normalization, the spatial order and
#     the element nucleus budget.  See the header of
#     species_face_flux_tests.f90.
#   * One row is a separate invocation, because the operator STOPS on an
#     injected violation of the face identity and a process that stops
#     cannot report its own row.
#   * The whole-binary rows state what only a run can show: the new
#     operation is refusable at its injection point and the step is retaken
#     from the checkpoint.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Exit status is nonzero if any row fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped.
# EXHALE_TEST_OUT selects where the test executable is written.
# EXHALE_SPECIES_EXE selects the binary the whole-binary rows run; without
# it those rows are skipped rather than building the production binary here.
#
# Usage: src/tests/species_face_flux/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/species_face_flux}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

rc=0

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL species_face_flux_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL species_face_flux_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/species_face_flux_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/species_face_flux_tests.x" \
    "$HERE/species_face_flux_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL species_face_flux_build measured=compile_error reference=ok tol=0"
   exit 1; }

"$OUT/species_face_flux_tests.x" || rc=1

# ------------------------------------------------------------------ #
# THE ASSERTION FIRES.  An assertion no run has ever made fire is not known
# to fire, so one face is given a deliberate violation and the operator is
# expected to stop with the face named.
# ------------------------------------------------------------------ #
if "$OUT/species_face_flux_tests.x" inject > "$OUT/inject.log" 2>&1; then
   echo "FAIL face_identity_assertion_fires measured=exit_0 reference=stop tol=0"
   rc=1
elif grep -q 'do not sum to the mass flux at face' "$OUT/inject.log"; then
   echo "PASS face_identity_assertion_fires measured=stop reference=stop tol=0"
else
   echo "FAIL face_identity_assertion_fires measured=other_stop reference=stop tol=0"
   echo "     see $OUT/inject.log"
   rc=1
fi

# ------------------------------------------------------------------ #
# THE WHOLE-BINARY ROWS.  mol_diffusion is the matrix case that runs the
# element transport, so it is the case in which the new operation exists at
# all.  A refusal is injected after the element-transport operation
# (as_op_diffusion = 4), which is the injection point the advected
# composition is written into f_sp at; the run has to report the rejection
# with that operation named and go on, which is the checkpoint restoring the
# composition the attempt started from.
# ------------------------------------------------------------------ #
EXE="${EXHALE_SPECIES_EXE:-}"
if [ -z "$EXE" ]; then
   echo "  whole-binary rows skipped: set EXHALE_SPECIES_EXE to a built binary"
   exit $rc
fi
if [ ! -x "$EXE" ]; then
   echo "FAIL species_face_flux_binary measured=not_executable reference=$EXE tol=0"
   exit 1
fi

CASE="$ROOT/backup/regression/mol_diffusion"
WORK="$OUT/mol_diffusion"
rm -rf "$WORK"; mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/"
( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40 \
     EXHALE_REJECT_AFTER_OP=4 EXHALE_REJECT_AT_STEP=20 \
     EXHALE_REJECT_ATTEMPTS=1 "$EXE" > run.log 2>&1 )
brc=$?
if [ $brc -ne 0 ] && [ $brc -ne 2 ]; then
   echo "FAIL species_transport_refusable measured=exit_$brc reference=exit_0 tol=0"
   echo "     see $WORK/run.log"
   exit 1
fi
if grep -q 'element diffusion step' "$WORK/run.log"; then
   echo "PASS species_transport_refusable measured=rejected reference=rejected tol=0"
else
   echo "FAIL species_transport_refusable measured=no_rejection reference=rejected tol=0"
   echo "     see $WORK/run.log"
   rc=1
fi

exit $rc

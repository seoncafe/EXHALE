#!/bin/bash
# Build and run the spectrum-type tests of decision 13 (one spectrum type
# builds every band), docs/development_plan_20260905_rev3.md section 10.5.
# See README.md in this directory for what each one asserts.
#
# Every test compares a measured number, an exit status or a recorded line
# against a reference and prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion.  Lines that begin with two spaces or with DIAGNOSTIC
# are context, not verdicts.  Exit status is nonzero if any assertion fails.
#
# The Fortran driver links the PRODUCTION objects, so the photon grid it
# tests is the one the binary was built from; the shell gate runs the binary
# itself on copies of backup/regression/wasp_full.
#
# EXHALE_OBJDIR and EXHALE_EXE select another build (a private OBJDIR/EXE
# pair, as the ifx build uses): with either of them set the staleness check
# is skipped, because the objects then need not match the default build/.
#
# Usage: src/tests/spectrum_type/run.sh [test ...]
#        names: planck_field balmer_field balmer_quad n2_floor sed_edges
#               fuv_quad sed_semantics wasp121_sed spectrum_gate
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="$ROOT/build/tests/spectrum_type"
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
   echo "FAIL spectrum_type_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}${EXHALE_EXE:-}" ]; then
   # The driver links the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL spectrum_type_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
# Every object but the main program: the driver brings its own.
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

WANT="${*:-planck_field balmer_field balmer_quad fuv_quad n2_floor sed_edges sed_semantics wasp121_sed spectrum_gate}"
n_fail=0

want() { case " $WANT " in *" $1 "*) return 0;; *) return 1;; esac; }

run_one() {
   local name="$1"; shift
   echo ""
   echo "---- $name ----"
   "$@"
   local rc=$?
   if [ $rc -ne 0 ]; then n_fail=$((n_fail+1)); fi
   return 0
}

if want planck_field; then
   rm -f "$OUT/planck_photon_field.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/planck_photon_field.x" \
       "$HERE/planck_photon_field.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL planck_photon_field_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver runs where it stands: set_energy_vectors reads an optional
   # opacity.inp from the working directory, and there must be none.
   rm -rf "$OUT/planck_run"
   mkdir -p "$OUT/planck_run"
   [ -x "$OUT/planck_photon_field.x" ] && \
      run_one planck_photon_field env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/planck_run' && '$OUT/planck_photon_field.x'"
fi

if want balmer_field; then
   rm -f "$OUT/balmer_continuum_field.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/balmer_continuum_field.x" \
       "$HERE/balmer_continuum_field.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL balmer_continuum_field_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/balmer_run"
   mkdir -p "$OUT/balmer_run"
   [ -x "$OUT/balmer_continuum_field.x" ] && \
      run_one balmer_continuum_field env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/balmer_run' && '$OUT/balmer_continuum_field.x'"
fi

if want balmer_quad; then
   rm -f "$OUT/balmer_band_quadrature.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/balmer_band_quadrature.x" \
       "$HERE/balmer_band_quadrature.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL balmer_band_quadrature_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver writes its own synthetic table into the directory it runs
   # in, and set_energy_vectors reads an optional opacity.inp from there.
   rm -rf "$OUT/balmer_quad_run"
   mkdir -p "$OUT/balmer_quad_run"
   [ -x "$OUT/balmer_band_quadrature.x" ] && \
      run_one balmer_band_quadrature env OMP_NUM_THREADS=1 \
              EXHALE_TEST_ROOT="$ROOT" \
              sh -c "cd '$OUT/balmer_quad_run' && '$OUT/balmer_band_quadrature.x'"
fi

if want fuv_quad; then
   rm -f "$OUT/fuv_band_quadrature.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/fuv_band_quadrature.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/fuv_band_quadrature.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL fuv_band_quadrature_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver writes its own small spectrum tables into the directory
   # it runs in.
   rm -rf "$OUT/fuv_quad_run"
   mkdir -p "$OUT/fuv_quad_run"
   [ -x "$OUT/fuv_band_quadrature.x" ] && \
      run_one fuv_band_quadrature env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/fuv_quad_run' && '$OUT/fuv_band_quadrature.x'"
fi

if want n2_floor; then
   rm -f "$OUT/photon_grid_n2_floor.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/photon_grid_n2_floor.x" \
       "$HERE/photon_grid_n2_floor.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL photon_grid_n2_floor_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/n2_floor_run"
   mkdir -p "$OUT/n2_floor_run"
   [ -x "$OUT/photon_grid_n2_floor.x" ] && \
      run_one photon_grid_n2_floor env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/n2_floor_run' && '$OUT/photon_grid_n2_floor.x'"
fi

if want sed_edges; then
   rm -f "$OUT/loaded_sed_threshold_edges.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/loaded_sed_threshold_edges.x" \
       "$HERE/loaded_sed_threshold_edges.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL loaded_sed_threshold_edges_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/sed_edges_run"
   mkdir -p "$OUT/sed_edges_run"
   [ -x "$OUT/loaded_sed_threshold_edges.x" ] && \
      run_one loaded_sed_threshold_edges env OMP_NUM_THREADS=1 \
              EXHALE_TEST_ROOT="$ROOT" \
              sh -c "cd '$OUT/sed_edges_run' && '$OUT/loaded_sed_threshold_edges.x'"
fi

if want sed_semantics; then
   rm -f "$OUT/loaded_sed_table_semantics.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/loaded_sed_table_semantics.x" \
       "$HERE/loaded_sed_table_semantics.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL loaded_sed_table_semantics_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver writes its own synthetic tables into the directory it runs
   # in, and set_energy_vectors reads an optional opacity.inp from there.
   rm -rf "$OUT/sed_semantics_run"
   mkdir -p "$OUT/sed_semantics_run"
   [ -x "$OUT/loaded_sed_table_semantics.x" ] && \
      run_one loaded_sed_table_semantics env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/sed_semantics_run' && '$OUT/loaded_sed_table_semantics.x'"
fi

if want wasp121_sed; then
   rm -f "$OUT/wasp121b_sed_sensitivity.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/wasp121b_sed_sensitivity.x" \
       "$HERE/wasp121b_sed_sensitivity.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL wasp121b_sed_sensitivity_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver runs where it stands: set_energy_vectors reads an optional
   # opacity.inp from the working directory, and there must be none.
   rm -rf "$OUT/wasp121_sed_run"
   mkdir -p "$OUT/wasp121_sed_run"
   [ -x "$OUT/wasp121b_sed_sensitivity.x" ] && \
      run_one wasp121b_sed_sensitivity env OMP_NUM_THREADS=1 \
              EXHALE_TEST_ROOT="$ROOT" \
              sh -c "cd '$OUT/wasp121_sed_run' && '$OUT/wasp121b_sed_sensitivity.x'"
fi

if want spectrum_gate; then
   run_one spectrum_type_gate env EXHALE_EXE="$EXE" \
           EXHALE_TEST_WORK="$OUT" bash "$HERE/spectrum_type_gate.sh"
fi

echo ""
if [ $n_fail -gt 0 ]; then
   echo "spectrum_type: $n_fail test program(s) reported a failing assertion"
   exit 1
fi
echo "spectrum_type: every test program passed"

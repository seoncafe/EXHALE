#!/bin/bash
# Build and run the charge-exchange row tests of item D7 (the equation basis
# of the helium row, and who owns the cell's rate coefficients).
# See README.md in this directory for what each one asserts.
#
# Every test compares a measured number or an exit status against a
# reference and prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion.  Lines that begin with two spaces or with DIAGNOSTIC
# are context, not verdicts.  Exit status is nonzero if any assertion fails.
#
# The Fortran drivers link the PRODUCTION objects, so the reaction set and
# the rate coefficients they test are the ones the binary was built from.
#
# EXHALE_OBJDIR and EXHALE_EXE select another build (a private OBJDIR/EXE
# pair): with either of them set the staleness check is skipped, because the
# objects then need not match the default build/.
#
# Usage: src/tests/charge_exchange_rows/run.sh [test ...]
#        names: reaction_sources stage_sources pair_reservoir
#               advection_pair metal_helium stale_cell
#               stale_cell_stage stale_cell_turnover turnover_rates
#               transport_refusal
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" charge_exchange_rows
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="$ROOT/build/tests/charge_exchange_rows"
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
   echo "FAIL charge_exchange_rows_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}${EXHALE_EXE:-}" ]; then
   # The drivers link the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL charge_exchange_rows_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
# Every object but the main program: the drivers bring their own.
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

WANT="${*:-reaction_sources stage_sources pair_reservoir advection_pair metal_helium stale_cell stale_cell_stage stale_cell_turnover turnover_rates transport_refusal}"
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

if want reaction_sources; then
   rm -f "$OUT/reaction_source_orientation.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/reaction_source_orientation.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/reaction_source_orientation.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL reaction_source_orientation_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/reaction_run"
   mkdir -p "$OUT/reaction_run"
   [ -x "$OUT/reaction_source_orientation.x" ] && \
      run_one reaction_source_orientation env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/reaction_run' && '$OUT/reaction_source_orientation.x'"
fi

if want stage_sources; then
   rm -f "$OUT/transported_stage_sources.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/transported_stage_sources.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/transported_stage_sources.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL transported_stage_sources_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/stage_source_run"
   mkdir -p "$OUT/stage_source_run"
   [ -x "$OUT/transported_stage_sources.x" ] && \
      run_one transported_stage_sources env OMP_NUM_THREADS=4 \
              sh -c "cd '$OUT/stage_source_run' && '$OUT/transported_stage_sources.x'"
fi

if want pair_reservoir; then
   rm -f "$OUT/he_h_pair_reservoir.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/he_h_pair_reservoir.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/he_h_pair_reservoir.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL he_h_pair_reservoir_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/pair_run"
   mkdir -p "$OUT/pair_run"
   [ -x "$OUT/he_h_pair_reservoir.x" ] && \
      run_one he_h_pair_reservoir env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/pair_run' && '$OUT/he_h_pair_reservoir.x'"
fi

if want advection_pair; then
   rm -f "$OUT/advection_pair_reservoir.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/advection_pair_reservoir.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/advection_pair_reservoir.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL advection_pair_reservoir_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/advection_run"
   mkdir -p "$OUT/advection_run"
   [ -x "$OUT/advection_pair_reservoir.x" ] && \
      run_one advection_pair_reservoir env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/advection_run' && '$OUT/advection_pair_reservoir.x'"
fi

if want metal_helium; then
   rm -f "$OUT/metal_helium_reservoir.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/metal_helium_reservoir.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/metal_helium_reservoir.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL metal_helium_reservoir_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/metal_helium_run"
   mkdir -p "$OUT/metal_helium_run"
   [ -x "$OUT/metal_helium_reservoir.x" ] && \
      run_one metal_helium_reservoir env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/metal_helium_run' && '$OUT/metal_helium_reservoir.x'"
fi

if want stale_cell; then
   rm -f "$OUT/stale_cell_rates_refused.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/stale_cell_rates_refused.x" \
       "$HERE/stale_cell_rates_refused.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL stale_cell_rates_refused_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/stale_cell_run"
   mkdir -p "$OUT/stale_cell_run"
   if [ -x "$OUT/stale_cell_rates_refused.x" ]; then
      echo ""
      echo "---- stale_cell_rates_refused ----"
      ( cd "$OUT/stale_cell_run" &&                                       \
        OMP_NUM_THREADS=1 "$OUT/stale_cell_rates_refused.x" ) \
         > "$OUT/stale_cell_run/out.log" 2>&1
      rc=$?
      sed -n '1,40p' "$OUT/stale_cell_run/out.log" | sed 's/^/  /'
      ok=no
      if [ $rc -ne 0 ] &&                                                 \
         grep -q '(charge_exchange) ERROR' "$OUT/stale_cell_run/out.log"; then
         ok=yes
      fi
      if [ $ok = yes ]; then
         echo "PASS stale_cell_rates_refused measured=rc=$rc reference=nonzero_and_named tol=0"
      else
         echo "FAIL stale_cell_rates_refused measured=rc=$rc reference=nonzero_and_named tol=0"
         n_fail=$((n_fail+1))
      fi
   fi
fi

if want stale_cell_stage; then
   rm -f "$OUT/stale_cell_stage_sources_refused.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/stale_cell_stage_sources_refused.x" \
       "$HERE/stale_cell_stage_sources_refused.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL stale_cell_stage_sources_refused_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/stale_cell_stage_run"
   mkdir -p "$OUT/stale_cell_stage_run"
   if [ -x "$OUT/stale_cell_stage_sources_refused.x" ]; then
      echo ""
      echo "---- stale_cell_stage_sources_refused ----"
      ( cd "$OUT/stale_cell_stage_run" &&                                 \
        OMP_NUM_THREADS=1 "$OUT/stale_cell_stage_sources_refused.x" ) \
         > "$OUT/stale_cell_stage_run/out.log" 2>&1
      rc=$?
      sed -n '1,40p' "$OUT/stale_cell_stage_run/out.log" | sed 's/^/  /'
      ok=no
      if [ $rc -ne 0 ] &&                                                 \
         grep -q '(charge_exchange) ERROR' "$OUT/stale_cell_stage_run/out.log"; then
         ok=yes
      fi
      if [ $ok = yes ]; then
         echo "PASS stale_cell_stage_sources_refused measured=rc=$rc reference=nonzero_and_named tol=0"
      else
         echo "FAIL stale_cell_stage_sources_refused measured=rc=$rc reference=nonzero_and_named tol=0"
         n_fail=$((n_fail+1))
      fi
   fi
fi

if want stale_cell_turnover; then
   rm -f "$OUT/stale_cell_turnover_refused.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/stale_cell_turnover_refused.x" \
       "$HERE/stale_cell_turnover_refused.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL stale_cell_turnover_refused_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/stale_cell_turnover_run"
   mkdir -p "$OUT/stale_cell_turnover_run"
   if [ -x "$OUT/stale_cell_turnover_refused.x" ]; then
      echo ""
      echo "---- stale_cell_turnover_refused ----"
      ( cd "$OUT/stale_cell_turnover_run" &&                              \
        OMP_NUM_THREADS=1 "$OUT/stale_cell_turnover_refused.x" ) \
         > "$OUT/stale_cell_turnover_run/out.log" 2>&1
      rc=$?
      sed -n '1,40p' "$OUT/stale_cell_turnover_run/out.log" | sed 's/^/  /'
      ok=no
      if [ $rc -ne 0 ] &&                                                 \
         grep -q '(charge_exchange) ERROR' "$OUT/stale_cell_turnover_run/out.log"; then
         ok=yes
      fi
      if [ $ok = yes ]; then
         echo "PASS stale_cell_turnover_refused measured=rc=$rc reference=nonzero_and_named tol=0"
      else
         echo "FAIL stale_cell_turnover_refused measured=rc=$rc reference=nonzero_and_named tol=0"
         n_fail=$((n_fail+1))
      fi
   fi
fi

if want turnover_rates; then
   rm -f "$OUT/turnover_cell_rates.x"
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/turnover_cell_rates.x" \
       "$HERE/../physics_probe/assertion_report.f90" \
       "$HERE/turnover_cell_rates.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL turnover_cell_rates_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   rm -rf "$OUT/turnover_rates_run"
   mkdir -p "$OUT/turnover_rates_run"
   [ -x "$OUT/turnover_cell_rates.x" ] && \
      run_one turnover_cell_rates env OMP_NUM_THREADS=4 \
              sh -c "cd '$OUT/turnover_rates_run' && '$OUT/turnover_cell_rates.x'"
fi

if want transport_refusal; then
   run_one ionization_transport_with_metals env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" bash "$HERE/ionization_transport_with_metals.sh"
fi

echo ""
if [ $n_fail -gt 0 ]; then
   echo "charge_exchange_rows: $n_fail test program(s) reported a failing assertion"
   exit 1
fi
echo "charge_exchange_rows: every test program passed"

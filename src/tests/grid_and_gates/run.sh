#!/bin/bash
# Build and run the grid and gate tests of section 10.4 of
# docs/development_plan_20260905_rev3.md.  See README.md in this directory
# for what each one asserts and which verdict is expected at HEAD.
#
# Every test compares a measured number against a reference and prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion.  Lines that begin with two spaces or with DIAGNOSTIC
# are context, not verdicts.  Exit status is nonzero if any assertion fails.
#
# The Fortran drivers link the PRODUCTION objects in build/, so the grid, the
# photon grid and the right-hand side they test are the ones EXHALE.x was
# built from.  The shell tests run EXHALE.x itself on copies of regression
# cases; no regression directory is written to.
#
# EXHALE_OBJDIR and EXHALE_EXE select another build (a private OBJDIR/EXE
# pair, as the ifx build uses): with either of them set the staleness check
# is skipped, because the objects then need not match the default build/.
# EXHALE_TEST_OUT selects where the test executables are written.
#
# base_continuity is a DIAGNOSTIC that needs a run directory: it calls
# input_read and init on the input.inp and the restart pair that directory
# holds, so that the boundary it probes carries that case's grid, reservoir
# and options.  EXHALE_L26_STATE names it, and the probe writes
# EXHALE_setup.out and an ./output directory into it, so name a scratch copy
# of a run and not the run itself.  Without that variable the probe is
# skipped and the suite says so.  Its two printed groups differ in kind: the
# ROW, DERIV and JUMP lines are measurements for the memo, and the PASS|FAIL
# rows are invariants: that the boundary is a function of its argument
# (repeated evaluation, an evaluation of another state in between, a complete
# flux evaluation in between, a changed call order), and that the stationary
# operator puts the wind's own mass flux through the base face and through
# every face of a stationary state.  The base-face bound is 3e-5 of the
# window mean, the offset a resolved wind carries between a face value and
# the mean of a cell-centered product; EXHALE_L26_BASE_FLUX_TOL moves it for
# a state whose wind is resolved differently.
#
# Usage: src/tests/grid_and_gates/run.sh [test ...]
#        names: grid_width grid_window photon_quadrature threshold_edges
#               base_cell_width base_grid_key base_grid_pinned base_grid_restart
#               hydrostatic_residual base_branch base_continuity
#               free_outflow_boundary
#               flux_spread
#               output_state base_level restart_grid restart_round_trip
#               restart_intent evaluate_products restart_option_change
#               direct_steady_setup
#               sed_coverage coupled_carrier_h2 carrier_transport_inert
#               iontrans_atomic momentum_row fpe_traps carrier_bound_hold
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" grid_and_gates
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
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
   echo "FAIL grid_and_gates_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}${EXHALE_EXE:-}" ]; then
   # The drivers link the production objects, so a build the sources have
   # moved past would be tested instead of the sources.
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL grid_and_gates_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
# Every object but the main program: the drivers bring their own.
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

WANT="${*:-grid_width grid_window base_cell_width base_grid_key base_grid_pinned base_grid_restart photon_quadrature threshold_edges hydrostatic_residual base_branch base_continuity free_outflow_boundary flux_spread output_state base_level restart_grid restart_round_trip restart_intent evaluate_products restart_option_change direct_steady_setup sed_coverage coupled_carrier_h2 carrier_transport_inert iontrans_atomic momentum_row fpe_traps carrier_bound_hold}"
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

if want grid_width; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/grid_width_identity.x" \
       "$HERE/grid_width_identity.f90" \
       "$OBJDIR/parameters.o" "$OBJDIR/define_grid.o" || {
      echo "FAIL grid_width_identity_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   [ -x "$OUT/grid_width_identity.x" ] && \
      run_one grid_width_identity env OMP_NUM_THREADS=1 \
              "$OUT/grid_width_identity.x"
fi

if want grid_window; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/grid_window_indices.x" \
       "$HERE/grid_window_indices.f90" \
       "$OBJDIR/parameters.o" "$OBJDIR/define_grid.o" || {
      echo "FAIL grid_window_indices_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   [ -x "$OUT/grid_window_indices.x" ] && \
      run_one grid_window_indices env OMP_NUM_THREADS=1 \
              "$OUT/grid_window_indices.x"
fi

if want base_cell_width; then
   # The serializer lives in setup_report, which pulls in the rest of the
   # production objects; parameters.o alone would not link it.
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/base_cell_width_provenance.x" \
       "$HERE/base_cell_width_provenance.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL base_cell_width_provenance_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   [ -x "$OUT/base_cell_width_provenance.x" ] && \
      run_one base_cell_width_provenance env OMP_NUM_THREADS=1 \
              "$OUT/base_cell_width_provenance.x"
fi

if want base_grid_key; then
   run_one base_grid_key_reproduces_default env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/base_grid_key_reproduces_default.sh"
fi

if want base_grid_pinned; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/base_grid_pinned_width.x" \
       "$HERE/base_grid_pinned_width.f90" \
       "$OBJDIR/parameters.o" "$OBJDIR/define_grid.o" || {
      echo "FAIL base_grid_pinned_width_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   [ -x "$OUT/base_grid_pinned_width.x" ] && \
      run_one base_grid_pinned_width env OMP_NUM_THREADS=1 \
              "$OUT/base_grid_pinned_width.x"
fi

if want base_grid_restart; then
   run_one base_grid_restart_contract env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/base_grid_restart_contract.sh"
fi

if want photon_quadrature; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/photon_grid_quadrature.x" \
       "$HERE/photon_grid_quadrature.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL photon_grid_quadrature_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver runs where it stands: set_energy_vectors reads an optional
   # opacity.inp from the working directory, and there must be none.
   rm -rf "$OUT/photon_run"
   mkdir -p "$OUT/photon_run"
   [ -x "$OUT/photon_grid_quadrature.x" ] && \
      run_one photon_grid_quadrature env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/photon_run' && '$OUT/photon_grid_quadrature.x'"
fi

if want threshold_edges; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/photon_grid_threshold_edges.x" \
       "$HERE/photon_grid_threshold_edges.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL photon_grid_threshold_edges_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver runs where it stands: set_energy_vectors reads an optional
   # opacity.inp from the working directory, and there must be none.
   rm -rf "$OUT/threshold_run"
   mkdir -p "$OUT/threshold_run"
   [ -x "$OUT/photon_grid_threshold_edges.x" ] && \
      run_one photon_grid_threshold_edges env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/threshold_run' && '$OUT/photon_grid_threshold_edges.x'"
fi

if want hydrostatic_residual; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/hydrostatic_residual.x" \
       "$HERE/hydrostatic_residual.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL hydrostatic_residual_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # One grid per invocation: define_grid's Newton convergence flag is SAVEd,
   # so a second grid built in the same process is not the grid its N asks
   # for (the driver header states it).  Each run appends its records to
   # hydrostatic_ladder.dat in this directory; the last invocation, with no
   # argument, reads them back and gives the verdicts.
   rm -rf "$OUT/hydrostatic_run"
   mkdir -p "$OUT/hydrostatic_run"
   [ -x "$OUT/hydrostatic_residual.x" ] && \
      run_one hydrostatic_residual env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/hydrostatic_run' &&
                     for n in 250 500 1000 2000; do
                        '$OUT/hydrostatic_residual.x' \$n || exit 1
                     done
                     '$OUT/hydrostatic_residual.x'"
fi

if want base_branch; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/base_branch_discriminant.x" \
       "$HERE/base_branch_discriminant.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL base_branch_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   [ -x "$OUT/base_branch_discriminant.x" ] && \
      run_one base_branch env OMP_NUM_THREADS=1 \
              "$OUT/base_branch_discriminant.x"
fi

if want base_continuity; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/base_boundary_continuity_probe.x" \
       "$HERE/base_boundary_continuity_probe.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL base_continuity_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   if [ -z "${EXHALE_L26_STATE:-}" ]; then
      echo ""
      echo "---- base_continuity ----"
      echo "DIAGNOSTIC base_continuity skipped: set EXHALE_L26_STATE to a"
      echo "  scratch copy of a run directory (input.inp plus"
      echo "  output/Hydro_ioniz_IC.txt and output/Ion_species_IC.txt)"
   elif [ ! -f "$EXHALE_L26_STATE/input.inp" ]; then
      echo "FAIL base_continuity_state measured=no_input_inp reference=$EXHALE_L26_STATE tol=0"
      n_fail=$((n_fail+1))
   elif [ -x "$OUT/base_boundary_continuity_probe.x" ]; then
      run_one base_continuity env OMP_NUM_THREADS=1 \
              sh -c "cd '$EXHALE_L26_STATE' &&
                     '$OUT/base_boundary_continuity_probe.x' \
                     '$(basename "$EXHALE_L26_STATE")'"
   fi
fi

if want free_outflow_boundary; then
   $FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
       -o "$OUT/free_outflow_boundary.x" \
       "$HERE/free_outflow_boundary.f90" $PROD_OBJ $LAPACK || {
      echo "FAIL free_outflow_boundary_build measured=compile_error reference=ok tol=0"
      n_fail=$((n_fail+1)); }
   # The driver reads no input file and writes none; it runs where it stands.
   rm -rf "$OUT/outflow_run"
   mkdir -p "$OUT/outflow_run"
   [ -x "$OUT/free_outflow_boundary.x" ] && \
      run_one free_outflow_boundary env OMP_NUM_THREADS=1 \
              sh -c "cd '$OUT/outflow_run' && '$OUT/free_outflow_boundary.x'"
fi

if want flux_spread; then
   run_one mass_flux_spread_functional env EXHALE_EXE="$EXE" \
           python3 "$HERE/mass_flux_spread_functional.py"
fi

if want output_state; then
   # EXHALE_TEST_OUT as well: this program runs the binary in work
   # directories of its own, and without it they are written under the
   # tree's build/tests/ even when the caller named a scratch directory.
   run_one output_state_consistency env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/output_state_consistency.sh"
fi

if want base_level; then
   run_one base_level_single_statement env EXHALE_EXE="$EXE" \
           bash "$HERE/base_level_single_statement.sh"
fi

if want restart_grid; then
   run_one restart_grid_guard env EXHALE_EXE="$EXE" \
           bash "$HERE/restart_grid_guard.sh"
fi

if want restart_round_trip; then
   run_one restart_round_trip env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/restart_round_trip.sh"
fi

if want restart_intent; then
   run_one restart_intent_and_metadata env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/restart_intent_and_metadata.sh"
fi

if want evaluate_products; then
   run_one stationary_evaluate_products env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/stationary_evaluate_products.sh"
fi

if want restart_option_change; then
   run_one restart_option_change env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/restart_option_change.sh"
fi

if want direct_steady_setup; then
   run_one direct_steady_setup_report env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/direct_steady_setup_report.sh"
fi

if want sed_coverage; then
   run_one sed_coverage_stop env EXHALE_EXE="$EXE" \
           bash "$HERE/sed_coverage_stop.sh"
fi

if want coupled_carrier_h2; then
   run_one coupled_carrier_h2_row env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/coupled_carrier_h2_row.sh"
fi

if want carrier_transport_inert; then
   run_one carrier_transport_inert_report env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/carrier_transport_inert_report.sh"
fi

if want iontrans_atomic; then
   run_one ionization_transport_atomic env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/ionization_transport_atomic.sh"
fi

if want momentum_row; then
   run_one momentum_row_from_fluxes_only env EXHALE_EXE="$EXE" \
           bash "$HERE/momentum_row_from_fluxes_only.sh"
fi

if want carrier_bound_hold; then
   run_one carrier_movement_bound_hold env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/carrier_movement_bound_hold.sh"
fi

if want fpe_traps; then
   run_one floating_point_traps env EXHALE_EXE="$EXE" \
           EXHALE_TEST_OUT="$OUT" \
           bash "$HERE/floating_point_traps_clean.sh"
fi

echo ""
if [ $n_fail -gt 0 ]; then
   echo "grid_and_gates: $n_fail test program(s) reported a failing assertion"
   exit 1
fi
echo "grid_and_gates: every test program passed"

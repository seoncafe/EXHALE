#!/bin/bash
# The boundary-state suite: what a residual's lower boundary is a function
# of, and whether every consumer of one evaluation reads the same one.
#
# Every test prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion.  Lines beginning with two spaces or with DIAGNOSTIC are
# context, not verdicts.  Exit status is nonzero if any assertion fails.
#
# The tests here run the binary on copies of regression cases; no regression
# directory is written to.  EXHALE_EXE selects the binary and EXHALE_TEST_OUT
# where the work directories are written.
#
# Usage: src/tests/boundary_state/run.sh [test ...]
#        names: installed_composition declared_inputs molecular_seed contact_law
#
# contact_law is a Fortran driver and links the production objects, so it
# needs an object directory: EXHALE_OBJDIR, default build/.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" boundary_state
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

# LAPACK, resolved the way the Makefile resolves it: the library of the
# compiler's own prefix when that prefix carries an OpenBLAS, else -llapack.
FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

mkdir -p "$OUT"

WANT="${*:-installed_composition declared_inputs molecular_seed contact_law}"
n_fail=0

want() { case " $WANT " in *" $1 "*) return 0;; *) return 1;; esac; }

run_one() {
   local name="$1"; shift
   echo ""
   echo "---- $name ----"
   "$@"
   if [ $? -ne 0 ]; then n_fail=$((n_fail+1)); fi
   return 0
}

if want installed_composition; then
   run_one boundary_of_the_installed_composition \
      env EXHALE_EXE="$EXE" EXHALE_TEST_OUT="$OUT" \
      bash "$HERE/boundary_of_the_installed_composition.sh"
fi

if want declared_inputs; then
   run_one boundary_from_the_declared_inputs \
      env EXHALE_EXE="$EXE" EXHALE_TEST_OUT="$OUT" \
      bash "$HERE/boundary_from_the_declared_inputs.sh"
fi

if want molecular_seed; then
   run_one molecular_seed_of_an_atomic_state \
      env EXHALE_EXE="$EXE" EXHALE_TEST_OUT="$OUT" \
      bash "$HERE/molecular_seed_of_an_atomic_state.sh"
fi

if want contact_law; then
   # Every object but the main program: the driver brings its own.
   PROD_OBJ="$(ls "$OBJDIR"/*.o 2>/dev/null | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"
   if [ -z "$PROD_OBJ" ]; then
      echo "FAIL contact_law_build measured=no_objects reference=$OBJDIR tol=0"
      n_fail=$((n_fail+1))
   else
      mkdir -p "$OUT/contact_law"
      $FC $FFLAGS_TEST -J"$OUT/contact_law" -I"$OBJDIR" \
          -o "$OUT/contact_law/contact_law_at_the_base.x" \
          "$HERE/contact_law_at_the_base.f90" $PROD_OBJ $LAPACK || {
         echo "FAIL contact_law_build measured=compile_error reference=ok tol=0"
         n_fail=$((n_fail+1)); }
      [ -x "$OUT/contact_law/contact_law_at_the_base.x" ] && \
         run_one contact_law env OMP_NUM_THREADS=1 \
                 sh -c "cd '$OUT/contact_law' &&
                        './contact_law_at_the_base.x'"
   fi
fi

echo ""
if [ $n_fail -gt 0 ]; then
   echo "boundary_state: $n_fail test program(s) reported a failing assertion"
   exit 1
fi
echo "boundary_state: every test program passed"

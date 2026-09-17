#!/bin/bash
# The discrete kinetic-energy quadratic form of the gated fourth-difference
# momentum stress of src/modules/flux/low_mach_dissipation.f90
# (PLAN_20260916_rev3 section 5, step 1).
#
# The driver is standalone: it re-forms the stress from the equations the
# module header states and needs no production object, so it can be run
# while the tree is being edited.  It links LAPACK for dsyev.
#
# Each assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Lines beginning with DIAGNOSTIC or with two spaces are context.
# The measured quantity is the smallest eigenvalue of the symmetric part of
# the operator A of E' = -v^T A v, divided by its largest: nonnegative
# dissipation is measured >= 0.
#
# The uniform-grid and varying-coefficient configurations are verdicts.  The
# configuration built on the actual grid and state of an EXHALE run is a
# DIAGNOSTIC for the existing stress and a verdict for the corrected form.
#
# Usage: src/tests/low_mach_stress_energy/run.sh [Hydro_ioniz.txt]
# with the default state the LHS 1140 b 0.02-XUV transient of item L25.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="${EXHALE_TEST_OUTDIR:-$ROOT/build/tests/low_mach_stress_energy}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O2 -g -fbacktrace}"
STATE="${1:-$ROOT/LHS1140b/models/.L14/x002_HeH2.13/output/Hydro_ioniz.txt}"

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

if ! $FC $FFLAGS_TEST "$HERE/low_mach_stress_energy_form.f90" \
        -o "$OUT/low_mach_stress_energy_form.x" $LAPACK; then
   echo "FAIL low_mach_stress_energy_build measured=compile_error reference=0 tol=0"
   exit 1
fi

if [ ! -f "$STATE" ]; then
   echo "DIAGNOSTIC state file absent, running the synthetic configurations only: $STATE"
   STATE=""
fi

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
   "$OUT/low_mach_stress_energy_form.x" $STATE

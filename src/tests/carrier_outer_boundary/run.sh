#!/usr/bin/env bash
# Build and run the carrier outer-boundary tests.
#
# EXHALE_TEST_CURVED=1 gives the carrier profile a quadratic part, which is
# what makes a ghost that is not a copy reach the outflow face and what makes
# the limited reconstruction differ from the donor cell average.  The driver
# is run once per profile.
#
# The link line follows src/tests/carrier_boundary_jacobian/run.sh: the
# driver reaches the frozen cell state bg_cell, whose module closure pulls in
# the equilibrium solvers, and those call MINPACK and LAPACK as free
# subroutines that source_closure.py cannot follow.
#
# EXHALE_TEST_OBJDIR selects another object directory, so that two builds of
# this suite can run side by side.  EXHALE_OBJDIR names the build tree the
# generated build_stamp.f90 is taken from.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/carrier_outer_boundary}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -f "$prefix/lib/libopenblas.so" ]; then
  lapack=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib" -ldl)
else
  lapack=(-llapack -ldl)
fi

extra=(
  "$probe/assertion_report.f90"
  "$root/src/modules/nonlinear_system_solver/dpmpar.f90"
  "$root/src/modules/nonlinear_system_solver/enorm.f90"
  "$root/src/modules/nonlinear_system_solver/qform.f90"
  "$root/src/modules/nonlinear_system_solver/qrfac.f90"
  "$root/src/modules/nonlinear_system_solver/r1mpyq.f90"
  "$root/src/modules/nonlinear_system_solver/r1updt.f90"
  "$root/src/modules/nonlinear_system_solver/dogleg.f90"
  "$root/src/modules/nonlinear_system_solver/fdjac1.f90"
  "$root/src/modules/nonlinear_system_solver/hybrd.f90"
  "$root/src/modules/nonlinear_system_solver/hybrd1.f90"
)

objdir="${EXHALE_OBJDIR:-$root/build}"
stamp="$objdir/build_stamp.f90"
if [ ! -f "$stamp" ]; then
  stamp="$work/build_stamp.f90"
  cat > "$stamp" <<'STAMP'
      module build_stamp
      ! Written by src/tests/carrier_outer_boundary/run.sh when the Makefile
      ! has not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

driver="$here/carrier_outer_boundary.f90"
order=$(python3 "$probe/source_closure.py" "$root" "$driver") || { echo "FAIL source_closure"; exit 1; }
order=$(printf '%s\n' "$order" | grep -v '/build_stamp\.f90$')
order="$stamp
$order"

objects=()
for src in "${extra[@]}" $order; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  if [ "$src" != "$driver" ]; then objects+=("$obj"); fi
done

exe="$work/carrier_outer_boundary.x"
if ! "$FC" "${flags[@]}" "$work/carrier_outer_boundary.o" "${objects[@]}" \
       -o "$exe" "${lapack[@]}"; then
  echo "FAIL link carrier_outer_boundary"
  exit 1
fi

status=0
run_one() {
  # $1 a label for the setting, the rest the environment of the invocation
  local label=$1; shift
  echo "--- $label"
  env OMP_NUM_THREADS=1 "$@" "$exe" || status=1
}

run_one "the constant profile" EXHALE_TEST_CURVED=0
run_one "the curved profile"   EXHALE_TEST_CURVED=1

if [ "$status" -ne 0 ]; then
  echo "carrier_outer_boundary: FAILED"
else
  echo "carrier_outer_boundary: PASSED"
fi
exit "$status"

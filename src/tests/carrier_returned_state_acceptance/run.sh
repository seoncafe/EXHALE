#!/usr/bin/env bash
# Build and run the carrier returned-state acceptance test.
#
# WHY THIS SUITE IS NOT IN physics_probe.  Its driver uses
# diffusive_photochemistry, whose source closure reaches
# ionization_equilibrium and through it the equilibrium solvers, and those
# call MINPACK (hybrd1, hybrd, dpmpar) and LAPACK (dgesvd) as free
# subroutines.  source_closure.py follows `use` statements, so it cannot
# find them, and the physics_probe link line carries no library.  This
# script adds the MINPACK sources and the LAPACK the Makefile resolves,
# and nothing else.
#
# EXHALE_TEST_OBJDIR selects another object directory, so that two builds
# of this suite can run side by side.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/carrier_returned_state_acceptance}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

# The LAPACK the Makefile resolves for the production link.  The library is
# needed because constrained_chemical_equilibrium calls dgesvd; nothing in
# this suite calls it.
prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -f "$prefix/lib/libopenblas.so" ]; then
  lapack=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib")
else
  lapack=(-llapack)
fi

# MINPACK and the free subroutines the equilibrium solvers call by name,
# plus the shared verdict module of the Phase 0 test directories (it lives
# beside the physics_probe drivers, so source_closure.py does not find it).
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

# utilities.f90 reads the revision strings of build/build_stamp.f90, which
# the Makefile generates.  Without a build tree of its own to read, write a
# stamp that says so; nothing tested here depends on its contents.
if [ ! -f "$root/build/build_stamp.f90" ]; then
  cat > "$work/build_stamp.f90" <<'STAMP'
      module build_stamp
      ! Written by src/tests/carrier_returned_state_acceptance/run.sh when the
      ! Makefile has not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

driver="$here/carrier_returned_state_acceptance.f90"
order=$(python3 "$probe/source_closure.py" "$root" "$driver") || { echo "FAIL source_closure"; exit 1; }

objects=()
for src in "${extra[@]}" $order; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  if [ "$src" != "$driver" ]; then objects+=("$obj"); fi
done

exe="$work/carrier_returned_state_acceptance.x"
if ! "$FC" "${flags[@]}" "$work/carrier_returned_state_acceptance.o" "${objects[@]}" \
       -o "$exe" "${lapack[@]}"; then
  echo "FAIL link carrier_returned_state_acceptance"
  exit 1
fi

status=0
"$exe" || status=1

if [ "$status" -ne 0 ]; then
  echo "carrier_returned_state_acceptance: FAILED"
else
  echo "carrier_returned_state_acceptance: PASSED"
fi
exit "$status"

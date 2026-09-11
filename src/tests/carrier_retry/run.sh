#!/usr/bin/env bash
# Build and run the carrier checkpoint and retry-controller test.
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
# EXHALE_OBJDIR names the build tree the generated build_stamp.f90 is
# taken from; without one the script writes its own.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/carrier_retry}"
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

# utilities.f90 reads the revision strings of the build_stamp.f90 the
# Makefile generates in the object directory.  Without a build tree to read
# one from, write a stamp that says so; nothing tested here depends on its
# contents.  Whichever of the two it is, it is compiled below as the first
# file of the closure, so that build_stamp.mod stands in the module
# directory before utilities.f90 asks for it.
objdir="${EXHALE_OBJDIR:-$root/build}"
stamp="$objdir/build_stamp.f90"
if [ ! -f "$stamp" ]; then
  stamp="$work/build_stamp.f90"
  cat > "$stamp" <<'STAMP'
      module build_stamp
      ! Written by src/tests/carrier_retry/run.sh when the
      ! Makefile has not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

driver="$here/carrier_retry.f90"
order=$(python3 "$probe/source_closure.py" "$root" "$driver") || { echo "FAIL source_closure"; exit 1; }
# source_closure.py looks for the stamp in build/; compile the one selected
# above instead, first, since it uses no other module.
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

exe="$work/carrier_retry.x"
if ! "$FC" "${flags[@]}" "$work/carrier_retry.o" "${objects[@]}" \
       -o "$exe" "${lapack[@]}"; then
  echo "FAIL link carrier_retry"
  exit 1
fi

status=0
"$exe" || status=1

# THE THERMODYNAMIC CO CEILING IS GONE FROM THE SOURCES.  It clamped the
# transported CO to co_equilibrium_density with no rate behind it; the CO row
# now carries published destruction rates and the domain record of that model
# replaces the ceiling's counters.  A name below reappearing in src/ means a
# clamp or a counter came back with it, so the scan is an assertion and not a
# tidiness check.
ceiling_names='co_ceiling|pct_co_ceiling|carrier_co_ceiling_record|carrier_co_ceiling_cells|co_ceiling_seen|co_ceiling_applications|co_ceiling_cells_hit|co_ceiling_CO_removed'
self=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")
hits=$(grep -rIEl "$ceiling_names" "$root/src" 2>/dev/null | grep -vxF "$self" | tr '\n' ' ')
if [ -n "$hits" ]; then
  echo "FAIL no_CO_ceiling_identifier_in_src measured=$hits reference=none tol=0"
  status=1
else
  echo "PASS no_CO_ceiling_identifier_in_src measured=none reference=none tol=0"
fi

if [ "$status" -ne 0 ]; then
  echo "carrier_retry: FAILED"
else
  echo "carrier_retry: PASSED"
fi
exit "$status"

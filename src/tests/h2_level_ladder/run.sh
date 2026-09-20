#!/usr/bin/env bash
# Build and run the reduced statistical-equilibrium model of the H2
# ground-state ladder (item L31 step 4 of docs/PLAN_20260917.md).
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Lines that begin with '#' or with two spaces are context, not verdicts.
# Exit status is nonzero if any assertion fails or the build fails.
#
# The driver compiles the PRODUCTION sources it uses through
# src/tests/physics_probe/source_closure.py, so the ladder it compares
# against is the one the binary carries, and it needs no build tree of its
# own.  EXHALE_TEST_OUT (or EXHALE_TEST_WORK) selects a private object
# directory, which is what a concurrent item uses.
#
# The published collision tables of Lique (2015) and Jozwiak et al. (2024)
# are third-party material distributed in the workspace references/ tree
# beside the repository, not inside it.  Their absence is reported as a
# SKIP with the path that was looked for.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
work="${EXHALE_TEST_OUT:-${EXHALE_TEST_WORK:-$root/build/tests/h2_level_ladder}}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O2 -g -fbacktrace -J"$work" -I"$work")

# LAPACK, resolved the way the Makefile resolves it: the library of the
# compiler's own prefix when that prefix carries an OpenBLAS, else
# -llapack.  The statistical equilibrium is one dense linear solve.
fc_path="$(command -v "$FC" 2>/dev/null || true)"
prefix="$(cd "$(dirname "$fc_path")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$prefix/lib/libopenblas.so" ]; then
   lapack=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib" -ldl)
else
   lapack=(-llapack -ldl)
fi

driver="$here/h2_ladder_statistical_equilibrium.f90"
report="$root/src/tests/physics_probe/assertion_report.f90"

# The production modules the driver and the model use, together with the
# model itself and the assertion reporter, in an order in which each file
# follows the modules it uses.
order=$(python3 "$root/src/tests/physics_probe/source_closure.py" \
                "$root" "$report" "$driver") \
   || { echo "FAIL source_closure"; exit 1; }

# utilities.f90 reads the revision strings of the build_stamp.f90 the
# Makefile generates.  Without a build tree, write one that says so.
stamp="$root/build/build_stamp.f90"
if [ ! -f "$stamp" ]; then
  stamp="$work/build_stamp.f90"
  cat > "$stamp" <<'STAMP'
      module build_stamp
      ! Written by src/tests/h2_level_ladder/run.sh when the Makefile has
      ! not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi
order=$(printf '%s\n' "$order" | grep -v '/build_stamp\.f90$')
order="$stamp
$order"

objects=()
for src in $order; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  case "$src" in
    "$driver") ;;
    *) objects+=("$obj") ;;
  esac
done

exe="$work/h2_ladder_statistical_equilibrium.x"
if ! "$FC" "${flags[@]}" "$work/h2_ladder_statistical_equilibrium.o" \
     "${objects[@]}" "${lapack[@]}" -o "$exe"; then
  echo "FAIL link h2_ladder_statistical_equilibrium"
  exit 1
fi

status=0
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 EXHALE_TEST_ROOT="$root" \
   "$exe" || status=1

if [ "$status" -ne 0 ]; then
  echo "h2_level_ladder: FAILED"
else
  echo "h2_level_ladder: PASSED"
fi
exit "$status"

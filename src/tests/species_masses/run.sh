#!/usr/bin/env bash
# Build and run the species-mass tests.
#
# The driver links the PRODUCTION sources under src/modules, never a copy of
# them: src/tests/physics_probe/source_closure.py follows the `use` statements
# and prints the closure in an order in which each file follows the modules it
# uses.  Two things that closure cannot give, because they are not reached
# through a `use`, are added here:
#   * the MINPACK routines under src/modules/nonlinear_system_solver, which are
#     bare external subprograms and not modules;
#   * LAPACK, which the constrained-equilibrium solve calls directly.
# The carrier transport pulls the whole ionization sweep in behind it, which is
# why this suite exists beside src/tests/physics_probe rather than inside it:
# the physics_probe drivers link no library.
#
# Objects and module files go to build/tests/species_masses/, which this
# script creates, or to the directory EXHALE_TEST_OBJDIR names, so that two
# builds of this suite can run side by side.  EXHALE_OBJDIR names the build
# tree the generated build_stamp.f90 is taken from; without one the script
# writes its own.  One line per assertion:
#     PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/species_masses}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

# The LAPACK of the compiler's own prefix, by the rule the Makefile uses: a
# LAPACK built against another libgfortran runtime resolves the symbols and
# then misbehaves at run time.
prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -f "$prefix/lib/libopenblas.so" ]; then
  ldlibs=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib")
else
  ldlibs=(-llapack)
fi

drivers=(
  carrier_background_masses.f90
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
      ! Written by src/tests/species_masses/run.sh when the Makefile has not
      ! generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

# assertion_report lives with the physics probes; it is the same verdict line
# and the same failure count.  It is compiled first, by hand, because the
# closure script only walks src/modules and the driver's own directory.
if ! "$FC" "${flags[@]}" -c "$probe/assertion_report.f90" \
     -o "$work/assertion_report.o"; then
  echo "FAIL compile assertion_report.f90"
  exit 1
fi

driver_paths=()
for d in "${drivers[@]}"; do driver_paths+=("$here/$d"); done

order=$(python3 "$probe/source_closure.py" "$root" "${driver_paths[@]}") \
  || { echo "FAIL source_closure"; exit 1; }
# source_closure.py looks for the stamp in build/; compile the one selected
# above instead, first, since it uses no other module.
order=$(printf '%s\n' "$order" | grep -v '/build_stamp\.f90$')
order="$stamp
$order"

# The bare (non-module) subprograms of the MINPACK solver directory.  A file
# that is a PROGRAM of its own is skipped: that directory also holds the
# standalone constrained-equilibrium probe, whose main would collide with the
# driver's.
extra=()
for f in "$root"/src/modules/nonlinear_system_solver/*.f90; do
  if grep -qiE '^[[:space:]]*module[[:space:]]+[a-z]' "$f"; then continue; fi
  if grep -qiE '^[[:space:]]*program[[:space:]]+[a-z]' "$f"; then continue; fi
  extra+=("$f")
done

objects=()
status=0
for src in $order "${extra[@]}"; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  case " ${driver_paths[*]} " in
    *" $src "*) ;;
    *) objects+=("$obj") ;;
  esac
done
objects+=("$work/assertion_report.o")

for d in "${drivers[@]}"; do
  exe="$work/${d%.f90}.x"
  if ! "$FC" "${flags[@]}" "$work/${d%.f90}.o" "${objects[@]}" -o "$exe" \
       "${ldlibs[@]}"; then
    echo "FAIL link ${d%.f90}"
    status=1
    continue
  fi
  "$exe" || status=1
done

if [ "$status" -ne 0 ]; then
  echo "species_masses: FAILED"
else
  echo "species_masses: PASSED"
fi
exit "$status"

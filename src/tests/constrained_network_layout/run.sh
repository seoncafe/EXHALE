#!/usr/bin/env bash
# Build and run the constrained-network layout tests.
#
# The driver links the PRODUCTION sources under src/modules, never a copy of
# them: source_closure.py (kept with the physics probe tests) follows the
# `use` statements of the driver and the closure it prints is what gets
# compiled, in an order in which each file follows the modules it uses.  The
# in-tree MINPACK routines are plain external subroutines rather than
# modules, so no `use` statement reaches them and they are named here.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="$root/build/tests/constrained_network_layout"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

# The module carries an opt-in singular-value diagnostic of its Jacobian
# (dgesvd), so the driver needs the same LAPACK the Makefile picks: the one
# of the compiler's own prefix when it has an OpenBLAS, the system one
# otherwise. Override with LAPACK_LIBS, exactly as for `make`.
prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -z "${LAPACK_LIBS:-}" ]; then
  if [ -f "$prefix/lib/libopenblas.so" ]; then
    LAPACK_LIBS="-L$prefix/lib -lopenblas -Wl,-rpath,$prefix/lib"
  else
    LAPACK_LIBS="-llapack"
  fi
fi
read -r -a lapack <<< "$LAPACK_LIBS"

drivers=(
  transported_proton_constraint_layout.f90
)

# MINPACK hybrd and everything it calls: reached by a CALL, not by a USE.
minpack=(
  dpmpar.f90 enorm.f90 qrfac.f90 qform.f90 r1mpyq.f90 r1updt.f90
  dogleg.f90 fdjac1.f90 hybrd.f90 hybrd1.f90
)

# utilities.f90 reads the revision strings of build/build_stamp.f90, which
# the Makefile generates.  Without a build tree of its own to read, write a
# stamp that says so; nothing tested here depends on its contents.
if [ ! -f "$root/build/build_stamp.f90" ]; then
  cat > "$work/build_stamp.f90" <<'STAMP'
      module build_stamp
      ! Written by src/tests/constrained_network_layout/run.sh when the
      ! Makefile has not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

status=0
objects=()

# The verdict lines are the physics probe's; compile that module first so the
# closure below finds its module file.
if ! "$FC" "${flags[@]}" -c "$probe/assertion_report.f90" \
        -o "$work/assertion_report.o"; then
  echo "FAIL compile assertion_report.f90"
  exit 1
fi
objects+=("$work/assertion_report.o")

for m in "${minpack[@]}"; do
  src="$root/src/modules/nonlinear_system_solver/$m"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$work/${m%.f90}.o"; then
    echo "FAIL compile $m"
    exit 1
  fi
  objects+=("$work/${m%.f90}.o")
done

driver_paths=()
for d in "${drivers[@]}"; do driver_paths+=("$here/$d"); done

order=$(python3 "$probe/source_closure.py" "$root" "${driver_paths[@]}") \
  || { echo "FAIL source_closure"; exit 1; }

for src in $order; do
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

for d in "${drivers[@]}"; do
  exe="$work/${d%.f90}.x"
  if ! "$FC" "${flags[@]}" "$work/${d%.f90}.o" "${objects[@]}" "${lapack[@]}" -o "$exe"; then
    echo "FAIL link ${d%.f90}"
    status=1
    continue
  fi
  "$exe" || status=1
done

if [ "$status" -ne 0 ]; then
  echo "constrained_network_layout: FAILED"
else
  echo "constrained_network_layout: PASSED"
fi
exit "$status"

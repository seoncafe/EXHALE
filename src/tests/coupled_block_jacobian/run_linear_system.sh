#!/usr/bin/env bash
# Build and run coupled_block_linear_system.f90, the L32 measurement of the
# coupled block's linear system on one frozen iterate.
#
# WHY IT IS BUILT AND NOT RUN AS EXHALE.x. The four measurements it makes
# need several samples of ONE operator, and the operator of a block
# iteration includes the radiation and chemistry caches of modules that
# have no accessor a second process could be put back into. So the samples
# are taken inside one process, between one initialization and one exit,
# by a program that links the production objects and calls the same
# eval_residual and jacobian_action_of_direction the solve calls.
#
# The program runs in a run directory (input.inp, base.inp, output/ with
# the state). EXHALE_L32_CASE names it; the default is the checkpoint the
# L32 memo states.
#
# EXHALE_OBJDIR names the build tree whose objects are linked (a private
# OBJDIR, as a concurrent item's build uses). EXHALE_TEST_OBJDIR selects
# where this program's own objects go.
#
# One line per assertion:  PASS|FAIL <name> measured= reference=
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
objdir="${EXHALE_OBJDIR:-$root/build}"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/coupled_block_linear_system}"
mkdir -p "$work"

FC=${FC:-gfortran}
# The same optimization the production objects were built with, since they
# are linked unchanged; -J points at this program's own module directory
# and -I at the one the production objects wrote.
flags=(-O2 -fopenmp -J"$work" -I"$objdir")

prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -f "$prefix/lib/libopenblas.so" ]; then
  lapack=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib" -ldl)
else
  lapack=(-llapack -lblas -ldl)
fi

if [ ! -d "$objdir" ] || [ -z "$(ls "$objdir"/*.o 2>/dev/null)" ]; then
  echo "FAIL coupled_block_linear_system measured=no_objects reference=$objdir"
  exit 1
fi

# Every production object but the main program's.
objects=()
for o in "$objdir"/*.o; do
  case "$(basename "$o")" in
    EXHALE_main.o) continue ;;
  esac
  objects+=("$o")
done

driver="$here/coupled_block_linear_system.f90"
if ! "$FC" "${flags[@]}" -c "$driver" -o "$work/coupled_block_linear_system.o"; then
  echo "FAIL compile coupled_block_linear_system"
  exit 1
fi
exe="$work/coupled_block_linear_system.x"
if ! "$FC" "${flags[@]}" "$work/coupled_block_linear_system.o" \
       "${objects[@]}" -o "$exe" "${lapack[@]}"; then
  echo "FAIL link coupled_block_linear_system"
  exit 1
fi

case_dir="${EXHALE_L32_CASE:-}"
if [ -z "$case_dir" ]; then
  echo "     set EXHALE_L32_CASE to the run directory of the checkpoint"
  echo "PASS coupled_block_linear_system_built measured=$exe reference=built"
  exit 0
fi

status=0
( cd "$case_dir" && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$exe" ) || status=1
if [ "$status" -ne 0 ]; then
  echo "coupled_block_linear_system: FAILED"
else
  echo "coupled_block_linear_system: PASSED"
fi
exit "$status"

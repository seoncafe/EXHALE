#!/bin/bash
# "Ionization transport: True" runs with metals in the mixture and without
# them.
#
# QUANTITY UNDER TEST
#   Whether a run whose transported ionization stages and whose local sweep
#   answer the same H+ and He+ balance of a cell is allowed to start.  The
#   stage sources take the metal charge exchange of Huang et al. (2023)
#   Table 4 -- group A with hydrogen, group C with helium under cx_full,
#   and the group E electron capture -- through
#   charge_exchange::charge_exchange_stage_sources, the same reaction set
#   and the same rate coefficients the sweep's rows take through
#   cx_add_to_fvec (item D7b).  While those terms were missing the
#   combination was refused, and this fixture is what that refusal was
#   asserted by.
#
# WHAT IS RUN (two short runs, on scratch copies)
#   A  examples/14_diffusion with its metals.inp and "Ionization
#      transport: True": must run.
#   B  the same case with metals.inp removed: must run.  This is the
#      configuration the atomic ionization-transport fixture of
#      src/tests/grid_and_gates/ionization_transport_atomic.sh exercises.
#
# ASSERTIONS (two)
#   ionization_transport_with_metals_runs   A is accepted
#   ionization_transport_metal_free_runs    B is accepted
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" ionization_transport_with_metals
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/charge_exchange_rows}"
SRC="$ROOT/examples/14_diffusion"

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL ionization_transport_with_metals_exe no_binary "$EXE" 0
   exit 1
fi
if [ ! -f "$SRC/input.inp" ] || [ ! -f "$SRC/metals.inp" ]; then
   verdict FAIL ionization_transport_with_metals_fixture no_input "$SRC" 0
   exit 1
fi

run_case() {  # run_case <dir name> <keep metals: yes|no>
   local d="$WORK/$1"
   rm -rf "$d"; mkdir -p "$d"
   cp "$SRC"/*.inp "$d/"
   [ "$2" = no ] && rm -f "$d/metals.inp"
   printf '%s\n' 'Ionization transport: True' >> "$d/input.inp"
   ( cd "$d" && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        EXHALE_MAXSTEPS=3 "$EXE" > run.log 2>&1 )
   RC_LAST=$?
   return 0
}

run_case iontrans_metals yes
rcA=$RC_LAST
LA="$WORK/iontrans_metals/run.log"
ok=no
if [ $rcA -eq 0 ] && ! grep -q '(input_read) ERROR' "$LA" 2>/dev/null; then
   ok=yes
fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_with_metals_runs "rc=$rcA" rc=0_and_no_refusal 0

run_case iontrans_metal_free no
rcB=$RC_LAST
LB="$WORK/iontrans_metal_free/run.log"
ok=no
if [ $rcB -eq 0 ] && ! grep -q '(input_read) ERROR' "$LB" 2>/dev/null; then
   ok=yes
fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_metal_free_runs "rc=$rcB" rc=0_and_no_refusal 0

exit $(( n_fail > 0 ))

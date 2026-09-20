#!/bin/bash
# "Ionization transport: True" in an ATOMIC gas: the rows are carried, and
# the configuration the operator does not carry is still refused.
#
# QUANTITY UNDER TEST
#   What a run does with "Ionization transport: True" when the molecular
#   network is off.  The three ionization stages x(H II) per hydrogen
#   nucleus and x(He II), x(He III) per helium nucleus are stages of an
#   element, transported on that element's own nucleus face flux
#   (ionization_stage_transport, equation 1); a hydrogen and helium mixture
#   has them whether or not it has molecules, and their source in an atomic
#   gas is the H/He ionization balance of ion_residual_core, the rows the
#   local equilibrium sweep of that same gas solves.
#
# WHAT IS RUN (two short runs, on scratch copies; neither writes into an
# examples directory)
#   A  examples/14_diffusion (atomic, helium, He 2^3S, "He_diffusion: True")
#      plus "Ionization transport: True", capped at five steps, WITHOUT the
#      example's metals.inp: the transported stage sources are the H/He rows
#      of ion_residual_core and do not carry the metal charge exchange the
#      local sweep adds to the same rows, so a metal-bearing mixture is
#      refused (item D7, src/tests/charge_exchange_rows).  The run must be
#      accepted, the setup report must name the carried set, and the
#      certification of the state it writes must carry the three stage rows
#      and the two stage-sum entries -- which together say that the operator
#      was entered, that its rows were assembled on a frozen background, and
#      that the stage fluxes of each element sum to that element's nucleus
#      flux at every face.
#   B  examples/15_molecular ("Molecular chemistry: True", carriers NOT
#      transported) plus the same line.  Must be refused: with the molecular
#      carriers on their local root the molecular network would be evaluated
#      at two different molecular compositions of the same cell, once for
#      the stage rows and once for the sweep.
#
# ASSERTIONS (five)
#   ionization_transport_atomic_runs              A is accepted and exits 0
#   ionization_transport_atomic_setup_names_rows  A's setup report names the
#                                                 three carried fractions and
#                                                 the atomic row source
#   ionization_transport_atomic_stage_sum         A's certification carries
#                                                 both stage-sum entries
#   ionization_transport_atomic_rows_measured     A's certification evaluates
#                                                 the three stage rows
#   ionization_transport_molecular_local_refused  B is refused by name
#
# EXPECTED BEFORE THIS ITEM: RED on the first four (A aborted in input_read,
# which refused the key in an atomic gas); B was refused then as now.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" ionization_transport_atomic
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
ATOMIC="$ROOT/examples/14_diffusion"
MOL="$ROOT/examples/15_molecular"

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL ionization_transport_atomic_exe no_binary "$EXE" 0
   exit 1
fi
if [ ! -f "$ATOMIC/input.inp" ] || [ ! -f "$MOL/input.inp" ]; then
   verdict FAIL ionization_transport_atomic_fixture no_input "$ATOMIC $MOL" 0
   exit 1
fi

run_case() {  # run_case <dir name> <source case> <extra line> <maxsteps> [drop]
   local d="$WORK/$1"; local src="$2"; local extra="$3"; local ms="$4"
   rm -rf "$d"; mkdir -p "$d"
   cp "$src"/*.inp "$d/" 2>/dev/null || true
   # A fifth argument names an input file the case is run without.
   [ $# -ge 5 ] && rm -f "$d/$5"
   printf '%s\n' "$extra" >> "$d/input.inp"
   ( cd "$d" && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        EXHALE_MAXSTEPS="$ms" "$EXE" > run.log 2>&1 )
   RC_LAST=$?
   return 0
}

# ---- A: the atomic gas carries the three stages -----------------------!
run_case iontrans_atomic "$ATOMIC" 'Ionization transport: True' 5 metals.inp
rcA=$RC_LAST
LA="$WORK/iontrans_atomic/run.log"
SA="$WORK/iontrans_atomic/EXHALE_setup.out"

ok=no
if [ $rcA -eq 0 ] && ! grep -q '(input_read) ERROR' "$LA" 2>/dev/null; then
   ok=yes
fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_atomic_runs "rc=$rcA" rc=0_and_no_refusal 0

ok=no
if grep -q 'IONIZATION STATE of hydrogen and helium is TRANSPORTED' \
        "$SA" 2>/dev/null &&                                              \
   grep -q 'x(He III)' "$SA" 2>/dev/null &&                               \
   grep -q 'H/He ionization balance of' "$SA" 2>/dev/null; then ok=yes; fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_atomic_setup_names_rows "$ok" yes 0

nsum=$(grep -c 'ionization stage nucleus sum' "$LA" 2>/dev/null || true)
ok=no
if [ "${nsum:-0}" -ge 2 ] &&                                              \
   ! grep -E 'ionization stage nucleus sum (H|He) +' "$LA" 2>/dev/null     \
     | grep -q 'unavailable'; then ok=yes; fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_atomic_stage_sum "entries=${nsum:-0}" \
        at_least_2_evaluated 0

nrow=$(grep -E 'carrier balance (H\+|He\+|He\+\+) ' "$LA" 2>/dev/null     \
       | grep -c 'evaluated' || true)
ok=no
[ "${nrow:-0}" -ge 3 ] && ok=yes
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_atomic_rows_measured "evaluated=${nrow:-0}" \
        3 0

# ---- B: a molecular gas whose carriers are not transported ------------!
run_case iontrans_mol_local "$MOL" 'Ionization transport: True' 1
rcB=$RC_LAST
LB="$WORK/iontrans_mol_local/run.log"
ok=no
if [ $rcB -ne 0 ] &&                                                      \
   grep -q 'Molecular carrier transport: True' "$LB" 2>/dev/null &&       \
   grep -q '(input_read) ERROR' "$LB" 2>/dev/null; then ok=yes; fi
verdict $( [ $ok = yes ] && echo PASS || echo FAIL ) \
        ionization_transport_molecular_local_refused "rc=$rcB" \
        refused_by_name 0

exit $(( n_fail > 0 ))

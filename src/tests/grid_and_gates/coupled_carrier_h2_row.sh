#!/bin/bash
# A coupled steady solve on a molecular configuration must carry the H2 row.
#
# QUANTITY UNDER TEST
#   The exit status of EXHALE.x, and the text it prints, when
#   "Coupled carrier solve: True" is asked for together with
#   "Molecular chemistry: True" while "Molecular carrier transport" is off.
#   The check is in input_read (src/modules/files_IO/input_read.f90), after
#   every key is resolved, so it also sees the default that
#   "Molecular carrier transport" takes from the oxygen chemistry.
#
# WHY IT IS REFUSED
#   With the carriers eliminated, n(H2) is not an unknown of the stationary
#   system: it is whatever the local-equilibrium sweep returns for the
#   composition it was seeded with.  In the shielded layer the sweep cannot
#   return a content at all, because the fast chemistry cycles
#   H2 -> H2+ -> H3+ -> H2 without changing the number of H2 nuclei: the
#   local balance rows fix the partition among the molecular species and
#   leave their sum where the seed put it.  MEASURED on the hot Uranus
#   element state, 0.77 of any seed perturbation of the layer's H2 content
#   survives every pass, and the base cell's energy row amplifies that by
#   ~4e3, so two evaluations of ONE state differ by 3.2e5 times "Resid tol"
#   per unit relative seed change.  A residual that is not a function of its
#   unknowns has no root.  With "Molecular carrier transport: True" the H2
#   content is a Newton unknown solved from its own transport balance and the
#   same measurement reads 9.7e-10, below "Resid tol" = 1e-8.
#
# REFERENCE
#   docs/To_be_determined_by_user_20260906.md section 10, decided (a) on
#   2026-09-08: the combination is refused at startup with the remedy key
#   named, the way the SED-coverage and retired-key refusals work.  The
#   measurement is docs/steady_solver_design.md (item B5h).
#
# WHAT IS RUN (three short runs, none of them in an examples directory;
# each capped at one step)
#   A  examples/15_molecular ("Molecular chemistry: True", no base.inp,
#      carriers not transported) plus "Coupled carrier solve: True".
#      Must exit 1 and name the remedy key.
#   B  the same plus "Molecular carrier transport: True".  Must not print
#      the refusal; only its absence is asserted, since the run may stop
#      later for a reason of its own.
#   C  examples/14_diffusion ("He_diffusion: True", atomic) plus
#      "Coupled carrier solve: True".  The refusal is on the eliminated H2
#      of a molecular layer, not on the coupled route itself: an element
#      row carried in the same Newton is untouched.  Must not print the
#      refusal.
#
# ASSERTIONS (five)
#   coupled_carrier_eliminated_h2_stops      A exits 1
#   coupled_carrier_refusal_names_remedy     A's message names what was
#                                            asked, why it has no root, and
#                                            "Molecular carrier transport:
#                                            True" on one line
#   coupled_carrier_transported_h2_starts    B does not print the refusal
#   coupled_carrier_transported_h2_ran       B took its step and wrote a
#                                            profile
#   coupled_carrier_element_row_starts       C does not print the refusal
#
# EXPECTED BEFORE THE B5i CHANGE: RED on the first two assertions (A ran to
# completion and printed nothing, exit 0).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
MOL="$ROOT/examples/15_molecular"
ELEM="$ROOT/examples/14_diffusion"

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL coupled_carrier_binary no_binary "$EXE" 0
   exit 1
fi
for f in "$MOL/input.inp" "$ELEM/input.inp"; do
   if [ ! -f "$f" ]; then
      verdict FAIL coupled_carrier_case missing "$f" 0
      exit 1
   fi
done

for d in ccA ccB ccC; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d/output"
done
cp "$MOL/input.inp"  "$WORK/ccA/input.inp"
cp "$MOL/input.inp"  "$WORK/ccB/input.inp"
cp "$ELEM/input.inp" "$WORK/ccC/input.inp"
[ -f "$ELEM/metals.inp" ] && cp "$ELEM/metals.inp" "$WORK/ccC/"
echo 'Coupled carrier solve: True'      >> "$WORK/ccA/input.inp"
echo 'Coupled carrier solve: True'      >> "$WORK/ccB/input.inp"
echo 'Molecular carrier transport: True'>> "$WORK/ccB/input.inp"
echo 'Coupled carrier solve: True'      >> "$WORK/ccC/input.inp"

run_stage() {  # run_stage <dir>; echoes the exit status
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   echo $?
}

rcA=$(run_stage ccA)
rcB=$(run_stage ccB)
rcC=$(run_stage ccC)
echo "  exit status: A(H2 eliminated)=$rcA  B(H2 carried)=$rcB" \
     " C(element row, atomic)=$rcC"

# ---- stage A: the refusal ----
if [ "$rcA" -eq 1 ]; then
   verdict PASS coupled_carrier_eliminated_h2_stops "exit_$rcA" exit_1 0
else
   verdict FAIL coupled_carrier_eliminated_h2_stops "exit_$rcA" exit_1 0
   tail -n 5 "$WORK/ccA/run.log" | sed 's/^/     /'
fi
msgA="$WORK/ccA/run.log"
missA=""
grep -q 'Coupled carrier solve: True'         "$msgA" || missA="$missA asked_key"
grep -q 'Molecular chemistry: True'           "$msgA" || missA="$missA molecular_key"
grep -q 'not a function of'                   "$msgA" || missA="$missA reason"
grep -q 'Molecular carrier transport: True'   "$msgA" || missA="$missA remedy_key"
grep -q 'steady_solver_design.md'             "$msgA" || missA="$missA reference"
if [ -z "$missA" ]; then
   verdict PASS coupled_carrier_refusal_names_remedy all_named \
           "asked,molecular,reason,remedy,reference" 0
else
   verdict FAIL coupled_carrier_refusal_names_remedy "missing:${missA# }" \
           "asked,molecular,reason,remedy,reference" 0
fi

# ---- stage B: the carried H2 passes the check ----
# Only the absence of the refusal is asserted: what the run does after
# input_read is the business of the solver tests, not of this gate.
nB=$(grep -c 'needs the molecular carriers transported' "$WORK/ccB/run.log")
if [ "$nB" -eq 0 ]; then
   verdict PASS coupled_carrier_transported_h2_starts not_refused not_refused 0
else
   verdict FAIL coupled_carrier_transported_h2_starts refused not_refused 0
fi
if [ -s "$WORK/ccB/output/Hydro_ioniz.txt" ]; then
   verdict PASS coupled_carrier_transported_h2_ran written written 0
else
   verdict FAIL coupled_carrier_transported_h2_ran no_output written 0
   tail -n 5 "$WORK/ccB/run.log" | sed 's/^/     /'
fi

# ---- stage C: an element row in the same Newton is not the refused case ----
nC=$(grep -c 'needs the molecular carriers transported' "$WORK/ccC/run.log")
if [ "$nC" -eq 0 ]; then
   verdict PASS coupled_carrier_element_row_starts not_refused not_refused 0
else
   verdict FAIL coupled_carrier_element_row_starts refused not_refused 0
fi

exit $(( n_fail > 0 ))

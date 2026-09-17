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
#   D  examples/15_molecular plus "Molecular carrier transport: True" and
#      "Coupled carrier solve: On stall".  The third value of the key: the
#      alternation runs first and the block takes the state from the pass
#      at which the alternation stops approaching a joint fixed point.  It
#      must not be refused and it must be echoed as itself, since a value
#      silently read as False leaves an input file asking for the block and
#      a run that never enters it.
#   E  the same plus "Coupled carrier solve: Sometimes", a word the key has
#      no meaning for.  Must exit 1 and name the three values.
#   F  examples/15_molecular (H2 eliminated) plus "Coupled carrier solve:
#      On stall".  "On stall" reaches the same block one pass later, so the
#      eliminated-H2 refusal holds for it too.
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
#   coupled_carrier_on_stall_starts          D is not refused
#   coupled_carrier_on_stall_echoed          D's resolved report carries
#                                            carrier_newton_on_stall T with
#                                            carrier_in_newton F
#   coupled_carrier_unknown_value_stops      E exits 1 and names the three
#                                            values
#   coupled_carrier_on_stall_eliminated_h2_stops
#                                            F exits 1 on the same refusal
#   coupled_carrier_handover_message_present the handover the outer loop
#                                            prints is in the binary's text
#                                            (that it FIRES is measured in a
#                                            stalling solve, not here)
#
# EXPECTED BEFORE THE I2 CHANGE: RED on the five rows above.  The entry-text
# parser reads word 4 alone and tests it against True, so "On stall" falls
# through to False: D is not refused but is echoed as the alternation, E is
# accepted silently, and F runs the alternation instead of refusing.
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

for d in ccA ccB ccC ccD ccE ccF; do
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
cp "$MOL/input.inp"  "$WORK/ccD/input.inp"
cp "$MOL/input.inp"  "$WORK/ccE/input.inp"
cp "$MOL/input.inp"  "$WORK/ccF/input.inp"
echo 'Molecular carrier transport: True'>> "$WORK/ccD/input.inp"
echo 'Coupled carrier solve: On stall'  >> "$WORK/ccD/input.inp"
echo 'Molecular carrier transport: True'>> "$WORK/ccE/input.inp"
echo 'Coupled carrier solve: Sometimes' >> "$WORK/ccE/input.inp"
echo 'Coupled carrier solve: On stall'  >> "$WORK/ccF/input.inp"

run_stage() {  # run_stage <dir>; echoes the exit status
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   echo $?
}

rcA=$(run_stage ccA)
rcB=$(run_stage ccB)
rcC=$(run_stage ccC)
rcD=$(run_stage ccD)
rcE=$(run_stage ccE)
rcF=$(run_stage ccF)
echo "  exit status: A(H2 eliminated)=$rcA  B(H2 carried)=$rcB" \
     " C(element row, atomic)=$rcC"
echo "  exit status: D(On stall)=$rcD  E(unknown value)=$rcE" \
     " F(On stall, H2 eliminated)=$rcF"

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

# ---- stage D: the third value of the key ----
nD=$(grep -c 'is not a value this key takes' "$WORK/ccD/run.log")
if [ "$nD" -eq 0 ] && [ "$rcD" -ne 1 ]; then
   verdict PASS coupled_carrier_on_stall_starts not_refused not_refused 0
else
   verdict FAIL coupled_carrier_on_stall_starts refused not_refused 0
   tail -n 5 "$WORK/ccD/run.log" | sed 's/^/     /'
fi
# The value has to be echoed as ITSELF: the block from the first pass and
# the handover on a stall are different runs, and a report that shows only
# carrier_in_newton cannot tell them apart.
onstall=$(awk '/^carrier_newton_on_stall/ { print $2 }' \
          "$WORK/ccD/EXHALE_resolved.out" 2>/dev/null)
innewton=$(awk '/^carrier_in_newton/ { print $2 }' \
          "$WORK/ccD/EXHALE_resolved.out" 2>/dev/null)
if [ "${onstall:-}" = T ] && [ "${innewton:-}" = F ]; then
   verdict PASS coupled_carrier_on_stall_echoed "on_stall=T,in_newton=F" \
           "on_stall=T,in_newton=F" 0
else
   verdict FAIL coupled_carrier_on_stall_echoed \
           "on_stall=${onstall:-unread},in_newton=${innewton:-unread}" \
           "on_stall=T,in_newton=F" 0
fi

# ---- stage E: a word the key has no meaning for ----
missE=""
grep -q 'is not a value this key takes' "$WORK/ccE/run.log" || missE="$missE reason"
grep -q 'False'    "$WORK/ccE/run.log" || missE="$missE false"
grep -q 'True'     "$WORK/ccE/run.log" || missE="$missE true"
grep -q 'On stall' "$WORK/ccE/run.log" || missE="$missE on_stall"
if [ "$rcE" -eq 1 ] && [ -z "$missE" ]; then
   verdict PASS coupled_carrier_unknown_value_stops "exit_${rcE}_all_named" \
           "exit_1,false,true,on_stall" 0
else
   verdict FAIL coupled_carrier_unknown_value_stops \
           "exit_${rcE}_missing:${missE# }" "exit_1,false,true,on_stall" 0
   tail -n 5 "$WORK/ccE/run.log" | sed 's/^/     /'
fi

# ---- stage F: "On stall" reaches the same block, so the same refusal ----
if [ "$rcF" -eq 1 ] && \
   grep -q 'needs the molecular carriers transported' "$WORK/ccF/run.log"; then
   verdict PASS coupled_carrier_on_stall_eliminated_h2_stops "exit_$rcF" \
           exit_1 0
else
   verdict FAIL coupled_carrier_on_stall_eliminated_h2_stops "exit_$rcF" \
           exit_1 0
   tail -n 5 "$WORK/ccF/run.log" | sed 's/^/     /'
fi

# ---- the handover the outer loop prints ----
# WHAT THIS ROW IS AND IS NOT.  It asserts that the sentence the outer loop
# prints at the handover exists in the binary; whether it FIRES is a
# property of a solve that stalls, which takes tens of outer passes and is
# measured in the item's own runs, not in a startup gate.
if strings "$EXE" | grep -q 'HANDOVER to the coupled block'; then
   verdict PASS coupled_carrier_handover_message_present present present 0
else
   verdict FAIL coupled_carrier_handover_message_present absent present 0
fi

exit $(( n_fail > 0 ))

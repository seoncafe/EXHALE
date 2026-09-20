#!/bin/bash
# "Molecular carrier transport: True" with the molecular network off says so.
#
# QUANTITY UNDER TEST
#   The exit status of EXHALE.x, and the text it prints, when
#   "Molecular carrier transport: True" is asked for in a run whose
#   "Molecular chemistry" is off.  The report is in input_read
#   (src/modules/files_IO/input_read.f90), after every key is resolved, so
#   it also sees the default the key takes from the oxygen chemistry.
#
# WHY IT IS REPORTED AND NOT REFUSED
#   The carriers this key transports are species of the molecular network
#   (H2, and OH/H2O/CO under the oxygen chemistry), and every consumer of
#   carrier_transport in the code is guarded by thereis_mol.  With the
#   network off the key therefore changes nothing: the state the run
#   produces is the atomic one it would have produced without the line.
#   The three ionization stages are NOT among them: they are stages of an
#   element and are carried on their own key, "Ionization transport", in
#   any gas.  What is wrong is the input file, not the
#   result, so the run continues and the reader is told -- the same
#   treatment "Molecular IR bands" and "Stellar LW flux" get when the
#   network they act on is absent.  A refusal would reject a run whose
#   physics is correct.
#
#   docs/input_schema.md K15d and the user manual both used to state this as
#   "Needs Molecular chemistry: True" while nothing checked it, which is the
#   defect this test closes (item HYG-CHECKED).
#
# WHAT IS RUN (two short runs, neither of them in an examples directory;
# each capped at one step)
#   A  examples/14_diffusion (atomic, "He_diffusion: True") plus
#      "Molecular carrier transport: True".  Must print the report, name
#      the remedy, and NOT stop.
#   B  examples/15_molecular ("Molecular chemistry: True") plus the same
#      line.  Must not print the report: there the carriers exist.
#
# ASSERTIONS (four)
#   carrier_transport_inert_reports        A prints the report
#   carrier_transport_inert_names_remedy   A's text names the key asked,
#                                          that it has no effect, and
#                                          "Molecular chemistry: True"
#   carrier_transport_inert_does_not_stop  A exits 0 and writes a profile
#   carrier_transport_with_network_silent  B does not print the report
#
# EXPECTED BEFORE THE HYG-CHECKED CHANGE: RED on the first two assertions
# (A ran to completion printing nothing about the key).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" carrier_transport_inert_report
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
   verdict FAIL carrier_transport_inert_binary no_binary "$EXE" 0
   exit 1
fi
for f in "$ATOMIC/input.inp" "$MOL/input.inp"; do
   if [ ! -f "$f" ]; then
      verdict FAIL carrier_transport_inert_case missing "$f" 0
      exit 1
   fi
done

for d in ctiA ctiB; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d/output"
done
cp "$ATOMIC/input.inp" "$WORK/ctiA/input.inp"
[ -f "$ATOMIC/metals.inp" ] && cp "$ATOMIC/metals.inp" "$WORK/ctiA/"
cp "$MOL/input.inp"    "$WORK/ctiB/input.inp"
echo 'Molecular carrier transport: True' >> "$WORK/ctiA/input.inp"
echo 'Molecular carrier transport: True' >> "$WORK/ctiB/input.inp"

run_stage() {  # run_stage <dir>; echoes the exit status
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   echo $?
}

rcA=$(run_stage ctiA)
rcB=$(run_stage ctiB)
echo "  exit status: A(atomic, key set)=$rcA  B(molecular, key set)=$rcB"

logA="$WORK/ctiA/run.log"
nA=$(grep -c 'Molecular carrier transport: True" is set but' "$logA")
if [ "$nA" -ge 1 ]; then
   verdict PASS carrier_transport_inert_reports reported reported 0
else
   verdict FAIL carrier_transport_inert_reports silent reported 0
fi

missA=""
grep -q 'Molecular carrier transport: True' "$logA" || missA="$missA asked_key"
grep -q 'Molecular chemistry" is off'       "$logA" || missA="$missA state"
grep -q 'no effect'                         "$logA" || missA="$missA consequence"
grep -q 'Molecular chemistry: True'         "$logA" || missA="$missA remedy_key"
if [ -z "$missA" ]; then
   verdict PASS carrier_transport_inert_names_remedy all_named \
           "asked,state,consequence,remedy" 0
else
   verdict FAIL carrier_transport_inert_names_remedy "missing:${missA# }" \
           "asked,state,consequence,remedy" 0
fi

# The atomic run is the correct one, so the report must not end it.
if [ "$rcA" -eq 0 ] && [ -s "$WORK/ctiA/output/Hydro_ioniz.txt" ]; then
   verdict PASS carrier_transport_inert_does_not_stop "exit_${rcA}_written" \
           exit_0_written 0
else
   verdict FAIL carrier_transport_inert_does_not_stop "exit_${rcA}" \
           exit_0_written 0
   tail -n 5 "$logA" | sed 's/^/     /'
fi

nB=$(grep -c 'Molecular carrier transport: True" is set but' \
          "$WORK/ctiB/run.log")
if [ "$nB" -eq 0 ]; then
   verdict PASS carrier_transport_with_network_silent silent silent 0
else
   verdict FAIL carrier_transport_with_network_silent reported silent 0
fi

exit $(( n_fail > 0 ))

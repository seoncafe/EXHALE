#!/bin/bash
# EXHALE_CARRIER_TRUST_HOLD: unset it changes nothing, set it holds the
# carrier movement bound at its value for every outer pass.
#
# QUANTITY UNDER TEST
#   The composition movement bound of the stationary outer iteration, the
#   largest fraction by which one pass may move the transported carriers.
#   The progress control halves it (to a floor of 1e-3) on a pass in which
#   neither the joint distance from certification nor the composition
#   residual fell, so the bound falls monotonically along a solve and the
#   movement a pass makes falls with it.  EXHALE_CARRIER_TRUST_HOLD=<value>
#   sets the bound to that value and takes the halving out of the pass.
#
# WHY IT EXISTS
#   A carrier row that decays while the bound is being halved is either a
#   relaxation the bound throttles, whose rate is proportional to the
#   bound, or a mode the alternation cannot damp, whose rate is the same
#   whatever the bound.  Holding the bound separates them
#   (docs/lhs1140b_stationary_L33_20260917.md section 6.3;
#   docs/lhs1140b_stationary_D8bound_20260918.md is the measurement).
#   It is a diagnostic and not an input key: nothing reads it unless the
#   environment names it.
#
# WHAT IS RUN (scratch copies of backup/regression/carrier_model_a_newton,
# the hot Uranus with the three ionization stages transported, entered from
# its pinned IC pair)
#   A  one outer pass, key unset, against one with the key set to an EMPTY
#      value.  The read path is not reached in either, so the two states
#      must be bitwise equal.
#   B  fourteen outer passes, key unset, against fourteen with the key held
#      at 1.0e-2, the value the run starts from.  MEASURED: the control's
#      bound is halved to 5.0e-3 at pass 14 and the held run's is not, and
#      the thirteen passes before that differ in nothing.
#
# ASSERTIONS (five)
#   carrier_bound_absent_is_inert          A's two states are bitwise equal
#   carrier_bound_absent_is_silent         neither A run announces a held bound
#   carrier_bound_hold_is_announced        B's held run states the held
#                                          bound (whether the control halves
#                                          is printed as a DIAGNOSTIC: it
#                                          turns on rounding of the mass row)
#   carrier_bound_key_holds                B's held run never leaves it and
#                                          the halving line never appears
#   carrier_bound_passes_before_agree      the passes before the control's
#                                          first halving agree in every
#                                          printed measure
#
# EXPECTED BEFORE THE KEY WAS ADDED: RED on carrier_bound_key_holds (the
# held run halved with the control), GREEN on the other four.
#
# COST: about four minutes at one thread, the two fourteen-pass runs in
# parallel.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" carrier_movement_bound_hold
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
CASE="$ROOT/backup/regression/carrier_model_a_newton"
PASSES="${EXHALE_BOUND_HOLD_PASSES:-14}"

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL carrier_bound_binary no_binary "$EXE" 0
   exit 1
fi
if [ ! -f "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   verdict FAIL carrier_bound_case missing "$CASE/IC" 0
   exit 1
fi

stage() {  # stage <dir>: a scratch copy entered from the pinned IC pair
   rm -rf "$WORK/$1"
   mkdir -p "$WORK/$1/output"
   cp "$CASE/input.inp" "$WORK/$1/input.inp"
   [ -f "$CASE/base.inp" ] && cp "$CASE/base.inp" "$WORK/$1/"
   cp "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
      "$WORK/$1/output/"
}

# ---- A: the key absent and the key empty ----------------------------
for d in cbhA1 cbhA2; do stage "$d"; done
( cd "$WORK/cbhA1" && OMP_NUM_THREADS=1 EXHALE_OUTER_PASSES=1 \
     "$EXE" > run.log 2>&1 )
( cd "$WORK/cbhA2" && OMP_NUM_THREADS=1 EXHALE_OUTER_PASSES=1 \
     EXHALE_CARRIER_TRUST_HOLD= "$EXE" > run.log 2>&1 )

# The provenance header stamps the RUN time, which two runs cannot share.
same=yes
for f in Hydro_ioniz.txt Ion_species.txt; do
   cmp -s <(grep -v '^# provenance:' "$WORK/cbhA1/output/$f") \
          <(grep -v '^# provenance:' "$WORK/cbhA2/output/$f") || \
      same="differs:$f"
done
if [ "$same" = yes ]; then
   verdict PASS carrier_bound_absent_is_inert bitwise_equal bitwise_equal 0
else
   verdict FAIL carrier_bound_absent_is_inert "$same" bitwise_equal 0
fi

nsil=$(cat "$WORK/cbhA1/run.log" "$WORK/cbhA2/run.log" | \
       grep -c 'movement bound is HELD at' || true)
if [ "$nsil" -eq 0 ]; then
   verdict PASS carrier_bound_absent_is_silent 0 0 0
else
   verdict FAIL carrier_bound_absent_is_silent "$nsil" 0 0
fi

# ---- B: the bound halved against the bound held ---------------------
for d in cbhB1 cbhB2; do stage "$d"; done
( cd "$WORK/cbhB1" && OMP_NUM_THREADS=1 EXHALE_OUTER_PASSES="$PASSES" \
     "$EXE" > run.log 2>&1 ) &
b1=$!
( cd "$WORK/cbhB2" && OMP_NUM_THREADS=1 EXHALE_OUTER_PASSES="$PASSES" \
     EXHALE_CARRIER_TRUST_HOLD=1.0e-2 "$EXE" > run.log 2>&1 ) &
b2=$!
wait $b1; wait $b2

# The bound each pass ran under, in order, as the pass summary prints it.
bounds() { grep -oE 'trust [0-9]\.[0-9]E-[0-9]{2},' "$WORK/$1/run.log" | \
           sed 's/trust //; s/,//'; }
bc=$(bounds cbhB1 | tr '\n' ' ')
bh=$(bounds cbhB2 | tr '\n' ' ')
echo "  bound by pass, control: $bc"
echo "  bound by pass, held:    $bh"

# Whether the CONTROL halves within the pass budget is a DIAGNOSTIC and not
# an assertion: the progress control halves the bound on a pass whose joint
# distance did not fall, and on this fixture that decision turns on the mass
# row moving at the level of its own rounding (MEASURED 2026-09-19: 3.13e-11
# to 7.16e-11 on one build, 4.06e-11 to 4.03e-11 on the next, so one halved
# at pass 14 and the other never did). A verdict on it would be a verdict on
# rounding. What the key does is asserted instead by the held run itself.
nlow=$(bounds cbhB1 | grep -c -v '^1.0E-02$' || true)
echo "  DIAGNOSTIC the control left 1.0E-02 on $nlow of $PASSES passes"
nann=$(grep -c 'bound is HELD at' "$WORK/cbhB2/run.log" || true)
if [ "$nann" -ge 1 ]; then
   verdict PASS carrier_bound_hold_is_announced "$nann" ge_1 0
else
   verdict FAIL carrier_bound_hold_is_announced 0 ge_1 0
fi

nheld=$(bounds cbhB2 | grep -c -v '^1.0E-02$' || true)
nhalv=$(grep -c 'carrier movement bound =' "$WORK/cbhB2/run.log" || true)
if [ "$nheld" -eq 0 ] && [ "$nhalv" -eq 0 ]; then
   verdict PASS carrier_bound_key_holds held_every_pass held_every_pass 0
else
   verdict FAIL carrier_bound_key_holds \
           "passes_off_the_value=$nheld,halvings=$nhalv" \
           held_every_pass 0
fi

# Every printed measure of the passes that ran under the same bound.  The
# trust field and the wall clock are the two that may differ.
first_low=$(bounds cbhB1 | grep -n -v '^1.0E-02$' | head -n 1 | cut -d: -f1)
[ -z "$first_low" ] && first_low=$((PASSES+1))
same_head=$( { for d in cbhB1 cbhB2; do
      grep -E 'outer pass .*: hydro info' "$WORK/$d/run.log" | \
        head -n $((first_low-1)) | \
        sed -E 's/, trust [0-9.E+-]+,//; s/[0-9]+\.[0-9]+ s//' | md5sum
   done; } | sort -u | wc -l)
if [ "$same_head" -eq 1 ]; then
   verdict PASS carrier_bound_passes_before_agree \
           "$((first_low-1))_passes_equal" equal 0
else
   verdict FAIL carrier_bound_passes_before_agree differ equal 0
   diff <(grep -E 'outer pass .*: hydro info' "$WORK/cbhB1/run.log") \
        <(grep -E 'outer pass .*: hydro info' "$WORK/cbhB2/run.log") | \
        head -n 6 | sed 's/^/     /'
fi

exit $(( n_fail > 0 ))

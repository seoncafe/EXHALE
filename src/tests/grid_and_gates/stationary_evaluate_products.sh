#!/bin/bash
# THE PRODUCTS OF "Restart intent: stationary evaluate", AND WHICH STATE EACH
# OF THEM DESCRIBES (item L9 of docs/PLAN_20260913_lhs_stationary.md,
# docs/PLAN_20260916_rev3.md section 9).
#
# WHAT IS UNDER TEST
#   1. The evaluation produces the complete set of files a run of this case
#      is read for: the equilibrium state pair, the channel breakdowns, the
#      advection-corrected pair and the mass-loss line.  The output directory
#      is created EMPTY except for a stale sentinel written into the file
#      each producer is supposed to write, so a producer that is skipped
#      leaves its sentinel behind and is caught; an assertion on file
#      existence alone would pass on a leftover from an earlier run.
#   2. The three states of the route are kept apart.  The LOADED state is
#      untouched by the evaluation; the WORK state, which the products are
#      made on, is not moved by the post-process that reads it; the
#      ADVECTION-DERIVED product carries the work state's certificate as
#      provenance and makes no certification claim of its own.
#   3. The two answers are reported apart: whether the file's own stationary
#      claim reproduces, and what verdict the work state gets.  A state whose
#      file claims a certification, evaluated under a different stellar EUV
#      luminosity from the one it was solved at, must say that the claim did
#      NOT reproduce and must still state the work state's own verdict.
#
# WHY THE COMPARISON IS IN MEMORY
#   The loader never writes its input files back, so an unchanged
#   Hydro_ioniz_IC.txt says nothing about whether the arrays the chemistry
#   and the post-process read were moved.  EXHALE_EVAL_STATE_DUMP makes the
#   evaluation write the conserved state, the composition and the live
#   composition-derived module state at four named points of the route, and
#   the rows below compare those points.
#
# THE TOLERANCES
#   There are none.  Every comparison here is an identity: an array that a
#   routine must not write to comes back with the same bits, or it was
#   written to.
#
# THE FIXTURE
#   backup/regression/wasp_full_newton/IC/, the certified state pair the
#   restart rows of this suite already read, on a copy.  The regression
#   directory is never written to.  That state does not re-certify on this
#   code generation (its cool column belongs to an earlier one,
#   docs/lhs1140b_stationary_L18_20260915.md section 7), so the run exits 2;
#   the products are written whatever the verdict is, which is the point.
#
# Usage: stationary_evaluate_products.sh   (EXHALE_EXE selects the binary)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
WORK="$OUT/stationary_evaluate_products"
CASE="$ROOT/backup/regression/wasp_full_newton"
SENTINEL='# STALE SENTINEL of a producer that did not run'

n_fail=0
check_eq() {   # check_eq <name> <measured> <reference>
   if [ "$2" = "$3" ]; then
      echo "PASS $1 measured=$2 reference=$3 tol=0"
   else
      echo "FAIL $1 measured=$2 reference=$3 tol=0"
      n_fail=$((n_fail + 1))
   fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL stationary_evaluate_products measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   echo "FAIL stationary_evaluate_products measured=no_fixture reference=$CASE/IC tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"

# stage <dir> [extra input lines...] : a run directory whose output/ holds
# the restart pair and a stale sentinel in every file a producer must write.
stage() {
   local d="$1"; shift
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/metals.inp" "$WORK/$d/"
   sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' \
      "$WORK/$d/input.inp"
   printf 'Restart intent: stationary evaluate\n' >> "$WORK/$d/input.inp"
   local line
   for line in "$@"; do printf '%s\n' "$line" >> "$WORK/$d/input.inp"; done
   cp "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
      "$WORK/$d/output/"
   local f
   for f in Hydro_ioniz.txt Ion_species.txt Hydro_ioniz_adv.txt \
            Ion_species_adv.txt Cooling_breakdown.txt Heating_breakdown.txt; do
      printf '%s\n' "$SENTINEL" > "$WORK/$d/output/$f"
   done
}

run_it() {   # run_it <dir> [env assignments...]
   local d="$1"; shift
   ( cd "$WORK/$d" && env OMP_NUM_THREADS=1 "$@" "$EXE" > run.log 2>&1 )
   RC_LAST=$?
   return 0
}

# ---- A: the products of one evaluation ----------------------------------
stage A
run_it A EXHALE_EVAL_STATE_DUMP=1
rcA=$RC_LAST
if [ "$rcA" != 0 ] && [ "$rcA" != 2 ]; then
   echo "FAIL stationary_evaluate_products_stage_A measured=exit_$rcA reference=0_or_2 tol=0"
   tail -n 20 "$WORK/A/run.log"
   exit 1
fi

for f in Hydro_ioniz.txt Ion_species.txt Hydro_ioniz_adv.txt \
         Ion_species_adv.txt Cooling_breakdown.txt Heating_breakdown.txt; do
   name="produced_$(echo "$f" | sed 's/\.txt$//')"
   if [ ! -s "$WORK/A/output/$f" ]; then
      check_eq "$name" missing written
   elif grep -q 'STALE SENTINEL' "$WORK/A/output/$f"; then
      check_eq "$name" stale_sentinel written
   else
      check_eq "$name" written written
   fi
done

# The mass-loss line the runner reads out of the log of this pass.
mdot=$(sed -n 's/.*steady-state Mdot *= *\([0-9.ED+-]*\).*/\1/p' \
       "$WORK/A/run.log" | tail -n 1)
check_eq mass_loss_line_written "$(test -n "$mdot" && echo present || echo absent)" present

# The certification metadata of the state the evaluation wrote.
check_eq certification_metadata_written \
   "$(grep -c '^# coupling: .*certified=' "$WORK/A/output/Hydro_ioniz.txt")" 1

# ---- the derived product's provenance, and its silence on certification --
check_eq adv_header_states_derived_from \
   "$(grep -c '^# derived_from: certified=' "$WORK/A/output/Hydro_ioniz_adv.txt")" 1
# The negative assertion: 'certified=' appears in the derived file ONLY on
# the provenance line, never as a coupling header of its own rows.
check_eq adv_header_makes_no_coupling_claim \
   "$(grep -c '^# coupling:' "$WORK/A/output/Hydro_ioniz_adv.txt")" 0
check_eq adv_certified_only_as_provenance \
   "$(grep -c 'certified=' "$WORK/A/output/Hydro_ioniz_adv.txt")" \
   "$(grep -c '^# derived_from: certified=' "$WORK/A/output/Hydro_ioniz_adv.txt")"

# ---- the three states, compared in memory -------------------------------
python3 - "$WORK/A/output/eval_state_dump.txt" <<'PY'
import sys
path = sys.argv[1]
blocks, order, cur, ns = {}, [], None, None
try:
    fh = open(path)
except OSError:
    print('FAIL in_memory_dump_written measured=missing reference=written tol=0')
    raise SystemExit(1)
for ln in fh:
    if ln.startswith('# block'):
        w = ln.split()
        cur = w[2]; ns = int(w[6]); blocks[cur] = []; order.append(cur)
        continue
    if ln.startswith('#') or not ln.strip():
        continue
    if cur is not None:
        blocks[cur].append(ln.split())

n_fail = 0
def check(name, measured, reference):
    global n_fail
    ok = measured == reference
    print('%s %s measured=%s reference=%s tol=0'
          % ('PASS' if ok else 'FAIL', name, measured, reference))
    if not ok:
        n_fail += 1

check('in_memory_dump_has_four_points', ','.join(order),
      'loaded,work_before_products,work_after_products,loaded_kept')
if n_fail:
    raise SystemExit(1)

# Columns 1..(5 + n_species) are the state the block names; the rest are the
# live composition-derived module state at the moment of the call.
st = 1 + 5 + ns
def state(b):  return [r[1:st] for r in blocks[b]]
def module(b): return [r[st:] for r in blocks[b]]

# The loaded state is kept: the arrays the route copied at entry are the
# arrays it still holds at exit, bit for bit.
check('loaded_state_untouched_in_memory',
      state('loaded') == state('loaded_kept'), True)
# ... and the comparison is not vacuous: the live module state at exit is the
# WORK state's and not the loaded state's, so the two blocks are two states.
check('loaded_block_and_work_block_are_two_states',
      module('loaded') != module('loaded_kept'), True)
# The post-process is a function of the state it is handed: it moves neither
# the conserved state and composition it describes nor the module state the
# next reader of that state would read.
check('products_do_not_move_the_work_state',
      state('work_before_products') == state('work_after_products'), True)
check('products_do_not_move_the_module_state',
      module('work_before_products') == module('work_after_products'), True)
raise SystemExit(1 if n_fail else 0)
PY
if [ $? -ne 0 ]; then n_fail=$((n_fail + 1)); fi

# ---- B: a certified file evaluated under a changed option ---------------
# The state claims certified=T; the model it is measured under is not the one
# it was solved under (the stellar EUV luminosity is doubled), so the claim
# cannot reproduce. The two answers must be reported as two.
stage B
sed -i -E 's/^Log10 of EUV luminosity \[erg\/s\]:.*/Log10 of EUV luminosity [erg\/s]: 30.72/' \
   "$WORK/B/input.inp"
run_it B
check_eq changed_option_original_claim_not_reproduced \
   "$(grep -c 'original claim: NOT REPRODUCED' "$WORK/B/run.log")" 1
check_eq changed_option_work_state_verdict_reported \
   "$(grep -c 'work state verdict:' "$WORK/B/run.log")" 1
check_eq changed_option_written_pair_reported \
   "$(grep -c 'the pair written into the state: certified=' "$WORK/B/run.log")" 1
check_eq changed_option_products_written \
   "$(test -s "$WORK/B/output/Hydro_ioniz_adv.txt" && \
      ! grep -q 'STALE SENTINEL' "$WORK/B/output/Hydro_ioniz_adv.txt" && \
      echo written || echo missing)" written

exit $((n_fail > 0))

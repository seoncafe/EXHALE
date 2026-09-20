#!/bin/bash
# Build and run the acceptance-class contract test: one meaning of the
# ionization acceptance classes for every consumer (classes 1, 2, 3 and 5 are
# roots of the requested equations; class 4 is not, and a state carrying one
# passes no final acceptance).
#
# The Fortran driver links the PRODUCTION objects, so the counting function
# and the acceptance predicate it tests are the ones the binary was built
# from.
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a parallel build
# uses); with it set the staleness check is skipped.
#
# Usage: src/tests/certification/run.sh
#   EXHALE_CERT_RUNLOG        one molecular run log (carrier round trip)
#   EXHALE_CERT_STATE_LOGS    run logs, colon-separated (closure entries)
#   EXHALE_CERT_CLAIM_LOGS    <run log>=<measured shell status>, colon-
#                             separated (the exit-status rule below)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/certification}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

# LAPACK, resolved the way the Makefile resolves it: the library of the
# compiler's own prefix when that prefix carries an OpenBLAS, else -llapack.
FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL certification_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     build the code first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL certification_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; build first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/certification_contexts.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/certification_contexts.x" \
    "$HERE/certification_contexts.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL certification_contexts_build measured=compile_error reference=ok tol=0"
   exit 1; }

n_fail=0
env OMP_NUM_THREADS=1 "$OUT/certification_contexts.x" || n_fail=$((n_fail+1))

# THE CAPACITY CONTRACT OF A REPORT (PLAN_20260919_rev1, item P5a).
# One equation past a full report stops the run: the refusal names the
# equation that did not fit and the capacity, and it happens before any
# field is written, so the entry that stands last is untouched. A driver
# cannot assert a stop from inside the process it stops, so the assertions
# are made here on the driver's exit status and its text.
rm -f "$OUT/certification_capacity.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/certification_capacity.x" \
    "$HERE/certification_capacity.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL certification_capacity_build measured=compile_error reference=ok tol=0"
   exit 1; }
cap_out="$(env OMP_NUM_THREADS=1 "$OUT/certification_capacity.x" 2>&1)"
cap_st=$?
cap_cap=$(echo "$cap_out" | sed -n 's/^filled entry_count=\([0-9]*\)$/\1/p')
if [ "$cap_st" -ne 0 ]; then
   echo "PASS one_equation_past_a_full_report_stops_the_run measured=$cap_st reference=nonzero tol=0"
else
   echo "FAIL one_equation_past_a_full_report_stops_the_run measured=$cap_st reference=nonzero tol=0"
   echo "     $(echo "$cap_out" | grep -E 'after |silent_overwrite=' | tr '\n' ' ')"
   n_fail=$((n_fail+1))
fi
if echo "$cap_out" | grep -q 'REPORT FULL: the equation "ionization stage nucleus sum He" does not fit' \
   && echo "$cap_out" | grep -q "capacity cert_max_entries = ${cap_cap:-0} entries"; then
   echo "PASS the_refusal_names_the_equation_and_the_capacity measured=named reference=named tol=0"
else
   echo "FAIL the_refusal_names_the_equation_and_the_capacity measured=not_named reference=named tol=0"
   n_fail=$((n_fail+1))
fi
# Nothing the report already holds moved: the entry that stands last is the
# sentinel the driver put there, with its own name and its own measure.
if echo "$cap_out" | grep -q 'last entry unchanged: ionization stage nucleus sum H  max= 4.200000E+01'; then
   echo "PASS a_refused_entry_overwrites_nothing measured=sentinel_intact reference=sentinel_intact tol=0"
else
   echo "FAIL a_refused_entry_overwrites_nothing measured=$(echo "$cap_out" | grep -c 'last entry unchanged') reference=sentinel_intact tol=0"
   n_fail=$((n_fail+1))
fi

# The isolated evaluation of the carrier balance asserts its own round trip
# at every call (certification.f90, carrier_rows_of_state): a run in which
# the module state was not reinstated prints a WARNING line.  A molecular run
# log given in EXHALE_CERT_RUNLOG is scanned for it; without one the check
# says so rather than passing silently.
LOG="${EXHALE_CERT_RUNLOG:-}"
if [ -n "$LOG" ] && [ -f "$LOG" ]; then
   nwarn=$(grep -c 'carrier module .*state was NOT reinstated' "$LOG" || true)
   if [ "$nwarn" -eq 0 ]; then
      echo "PASS carrier_module_state_reinstated_in_a_real_run measured=0 reference=0 tol=0"
   else
      echo "FAIL carrier_module_state_reinstated_in_a_real_run measured=$nwarn reference=0 tol=0"
      n_fail=$((n_fail+1))
   fi
else
   echo "  DIAGNOSTIC carrier_module_state_reinstated_in_a_real_run: no run log given (EXHALE_CERT_RUNLOG)"
fi

# THE STEP-4 EVALUATORS ON REAL RUNS.  EXHALE_CERT_STATE_LOGS is a
# colon-separated list of run logs.  For each one the final-state
# certification block is read and three things are asserted:
#   - the eliminated-species closure is `evaluated`, not `unavailable`;
#   - its row measure is a finite number;
#   - it reads `within` its tolerance exactly when the same block reports
#     zero cells without a chemical root, and `ABOVE` exactly when it
#     reports some.  That is the contract between the acceptance classes
#     and the closure measure: classes 1, 2, 3 and 5 are at or below
#     ieq_res_tol and class 4 is above it.
LOGS="${EXHALE_CERT_STATE_LOGS:-}"
if [ -n "$LOGS" ]; then
   IFS=':' read -r -a _logs <<< "$LOGS"
   for L in "${_logs[@]}"; do
      nm=$(basename "$(dirname "$L")")
      if [ ! -f "$L" ]; then
         echo "FAIL closure_entry_of_$nm measured=no_log reference=$L tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      blk=$(awk '/final state, as written/{f=1} f' "$L")
      line=$(echo "$blk" | grep 'eliminated-species closure' | head -n 1)
      if ! echo "$line" | grep -q 'evaluated'; then
         echo "FAIL closure_entry_of_$nm measured=not_evaluated reference=evaluated tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      echo "PASS closure_entry_of_$nm measured=evaluated reference=evaluated tol=0"
      val=$(echo "$line" | sed -n 's/.*max=\( *[^ ]*\).*/\1/p' | tr -d ' ')
      if echo "$val" | grep -qi 'nan\|infinity'; then
         echo "FAIL closure_measure_of_${nm}_is_finite measured=$val reference=finite tol=0"
         n_fail=$((n_fail+1))
      else
         echo "PASS closure_measure_of_${nm}_is_finite measured=$val reference=finite tol=0"
      fi
      nroot=$(echo "$blk" | sed -n 's/.*cells without a chemical root: \([0-9]*\).*/\1/p' | head -n 1)
      if echo "$line" | grep -q 'within'; then within=yes; else within=no; fi
      if [ "${nroot:-0}" -eq 0 ]; then want=yes; else want=no; fi
      if [ "$within" = "$want" ]; then
         echo "PASS closure_agrees_with_the_acceptance_classes_of_$nm measured=$within reference=$want tol=0"
      else
         echo "FAIL closure_agrees_with_the_acceptance_classes_of_$nm measured=$within reference=$want tol=0"
         n_fail=$((n_fail+1))
      fi
   done
else
   echo "  DIAGNOSTIC closure_entry_on_a_real_run: no run logs given (EXHALE_CERT_STATE_LOGS)"
fi

# The closure evaluation and the elemental transport residual assert their
# own round trip on every call: a run in which either was not reinstated
# prints a WARNING line.
if [ -n "$LOGS" ]; then
   IFS=':' read -r -a _logs <<< "$LOGS"
   for L in "${_logs[@]}"; do
      [ -f "$L" ] || continue
      nm=$(basename "$(dirname "$L")")
      w=$(grep -c 'was NOT reinstated' "$L" || true)
      if [ "$w" -eq 0 ]; then
         echo "PASS module_state_reinstated_in_$nm measured=0 reference=0 tol=0"
      else
         echo "FAIL module_state_reinstated_in_$nm measured=$w reference=0 tol=0"
         n_fail=$((n_fail+1))
      fi
   done
fi

# THE EXIT STATUS OF A RUN THAT DECLARED A STATIONARY STATE.
# A run that DECLARED one -- the du stop, the marching residual gate, or a
# JFNK finish -- and then wrote a state the certification refused exits 2 and
# names the refusing entries; a run that declared nothing exits 0. The status
# follows the CLAIM and not the run mode: the `init` exemption of the
# run-mode contract covers the physical-step claims (the clock, the budgets,
# the histories) and not a stationarity claim
# (docs/a0_run_mode_contract_20260906.md section 6,
# docs/a2_certification_contract_20260906.md section 8).
#
# A JFNK finish returning info = 2 is the case this gate was written for: the
# solve converged, the state it hands back does not meet the gate its last
# iterate met, that state IS the state the run writes, and the run declares
# it. `wasp_full_newton` ends there.
#
# EXHALE_CERT_CLAIM_LOGS is a colon-separated list of <run log>=<status>,
# where <status> is the shell status the caller measured when it produced
# that log. Without it the gate says so rather than passing silently.
CLAIMS="${EXHALE_CERT_CLAIM_LOGS:-}"
if [ -n "$CLAIMS" ]; then
   IFS=':' read -r -a _claims <<< "$CLAIMS"
   for E in "${_claims[@]}"; do
      L="${E%=*}"; ST="${E##*=}"
      nm=$(basename "$(dirname "$L")")
      if [ ! -f "$L" ]; then
         echo "FAIL exit_status_of_$nm measured=no_log reference=$L tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      # What the run reported about its own claim and its own certification.
      if grep -q 'exits with status 2 because it is not certified' "$L"
      then claimed_and_refused=yes; else claimed_and_refused=no; fi
      if grep -q 'made no stationary claim' "$L"
      then no_claim=yes; else no_claim=no; fi
      if awk '/final state, as written/{f=1} f' "$L" | grep -q 'NOT CERTIFIED'
      then final_refused=yes; else final_refused=no; fi
      # An uncertified final state leaves exactly one of the two lines above:
      # the claim was refused (status 2) or there was no claim (status 0).
      # Neither of them on an uncertified state means the status decision was
      # never taken on this run.
      if [ "$claimed_and_refused" = yes ]; then got="$ST"; want=2
      elif [ "$no_claim" = yes ]; then got="$ST"; want=0
      elif [ "$final_refused" = yes ]; then got=no_decision; want=decided_by_the_claim
      else got="$ST"; want=0; fi        # certified: nothing to refuse
      if [ "$got" = "$want" ]; then
         echo "PASS exit_status_of_$nm measured=$got reference=$want tol=0"
      else
         echo "FAIL exit_status_of_$nm measured=$got reference=$want tol=0"
         n_fail=$((n_fail+1))
      fi
      # The JFNK info = 2 route, where the claim is the solver's own: the
      # state is written, the certification refuses it, and the status is 2.
      if grep -q 'JFNK returned info = 2' "$L"; then
         if grep -q 'NOT CERTIFIED' "$L"; then r=refused; else r=certified; fi
         if [ "$r" = refused ] && [ "$ST" = 2 ]; then v=refused_and_2
         else v="${r}_and_${ST}"; fi
         echo "$([ "$v" = refused_and_2 ] && echo PASS || echo FAIL)" \
              "jfnk_info_2_is_a_refused_claim_in_$nm measured=$v" \
              "reference=refused_and_2 tol=0"
         [ "$v" = refused_and_2 ] || n_fail=$((n_fail+1))
      fi
      # The status is a verdict about a complete set of files, never a
      # reason for an incomplete one.
      d=$(dirname "$L")
      if [ -s "$d/output/Hydro_ioniz.txt" ] && [ -s "$d/output/Hydro_ioniz_adv.txt" ]
      then v=written; else v=missing; fi
      if [ "$v" = written ]; then
         echo "PASS outputs_complete_in_$nm measured=$v reference=written tol=0"
      else
         echo "FAIL outputs_complete_in_$nm measured=$v reference=written tol=0"
         n_fail=$((n_fail+1))
      fi
   done
else
   echo "  DIAGNOSTIC exit_status_of_a_declared_stationary_state: no run logs given (EXHALE_CERT_CLAIM_LOGS)"
fi

if [ $n_fail -gt 0 ]; then exit 1; fi

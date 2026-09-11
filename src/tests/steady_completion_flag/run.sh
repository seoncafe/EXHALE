#!/bin/bash
# The completion flag of a steady solve is a statement about the state it
# hands back, and the test the iteration stops on is the test that state is
# accepted by.
#
# What is asserted, on the log of a real run (the logs are given in
# EXHALE_STEADY_RUNLOGS, colon separated; each must contain one steady
# solve):
#
#   1. info = 0 exactly when every hydrodynamic row of the certification of
#      the state handed back is `within` its tolerance; a row `ABOVE` it or
#      `UNAVAILABLE` on it makes the flag nonzero. Which nonzero flag it is
#      says how the solve ended: 1 on the iteration budget, 2 on a state
#      the return refused or on an abort.
#   2. The row the solver names when it refuses carries the SAME number
#      and the SAME cell as that row in the certification block: one
#      record, read twice, never two evaluations.
#   3. THE ITERATION DOES NOT STOP FOR SUCCESS WHILE A CERTIFIED ROW
#      REFUSES. The stop test of the iteration is the acceptance test of
#      the state: the acceptance gate AND every row the certification
#      judges the carried system by, each against its own tolerance (mass
#      3e-12, momentum 1e-8, energy 1e-6). The gate alone is one number
#      against the run's `Resid tol`, so a solve stopping on it stops at
#      states the return then refuses. A solve that ends info = 0 must
#      therefore say that its stop test read those rows, and a refusal at
#      return -- which the composition the final evaluation adopts can
#      still cause -- must stand beside that same statement.
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# Usage: EXHALE_STEADY_RUNLOGS=<log>[:<log>...] src/tests/steady_completion_flag/run.sh
set -u

LOGS="${EXHALE_STEADY_RUNLOGS:-}"
if [ -z "$LOGS" ]; then
   echo "FAIL steady_completion_flag measured=no_log reference=EXHALE_STEADY_RUNLOGS tol=0"
   echo "     give it the log of a run whose steady solver ran"
   exit 1
fi

n_fail=0
IFS=':' read -r -a _logs <<< "$LOGS"
for L in "${_logs[@]}"; do
   nm=$(basename "$(dirname "$L")")
   if [ ! -f "$L" ]; then
      echo "FAIL steady_flag_of_$nm measured=no_log reference=$L tol=0"
      n_fail=$((n_fail+1)); continue
   fi
   # The last solve of the run: its certification block and its done line.
   # The LAST certification block of a steady solve, header to verdict line.
   blk=$(awk '/\(certification\) \((JFNK|PTC)\) state /{n=NR} {l[NR]=$0} \
              END{if(n) for(i=n;i<=NR;i++){print l[i]; if(l[i] ~ /CERTIFIED/) exit}}' "$L")
   done_line=$(grep -E '^ \((JFNK|PTC)\) done info=' "$L" | tail -n 1)
   if [ -z "$blk" ] || [ -z "$done_line" ]; then
      echo "FAIL steady_flag_of_$nm measured=no_steady_solve reference=one_solve tol=0"
      n_fail=$((n_fail+1)); continue
   fi
   info=$(echo "$done_line" | sed -n 's/.*done info=\([0-9-]*\).*/\1/p')
   rows=$(echo "$blk" | grep 'hydrodynamic .* row ')
   nrow=$(echo "$rows" | grep -c 'evaluated')
   nabove=$(echo "$rows" | grep -c 'ABOVE')
   nunavail=$(echo "$rows" | grep -c 'UNAVAILABLE')
   if [ "$nrow" -ne 3 ]; then
      echo "FAIL hydro_rows_of_$nm measured=$nrow reference=3 tol=0"
      n_fail=$((n_fail+1))
   else
      echo "PASS hydro_rows_of_$nm measured=3 reference=3 tol=0"
   fi
   # 1. the flag follows the rows of the state handed back.
   if [ "$nabove" -eq 0 ] && [ "$nunavail" -eq 0 ]; then
      if [ "$info" = "0" ]; then
         echo "PASS flag_follows_the_returned_state_of_$nm measured=$info reference=0 tol=0"
      else
         echo "FAIL flag_follows_the_returned_state_of_$nm measured=$info reference=0 tol=0"
         echo "     every hydrodynamic row of the state handed back is within its tolerance"
         n_fail=$((n_fail+1))
      fi
   else
      if [ "$info" != "0" ]; then
         echo "PASS flag_follows_the_returned_state_of_$nm measured=$info reference=nonzero tol=0"
      else
         echo "FAIL flag_follows_the_returned_state_of_$nm measured=$info reference=nonzero tol=0"
         echo "     rows above tolerance: $nabove, unavailable: $nunavail"
         n_fail=$((n_fail+1))
      fi
   fi
   # 2. The refusal the solver prints and the row of the certification.
   ref=$(grep -E 'refusing row: hydrodynamic' "$L" | tail -n 1)
   if [ -n "$ref" ]; then
      rname=$(echo "$ref" | sed -n 's/.*refusing row: \(hydrodynamic [a-z]* row\).*/\1/p')
      rval=$(echo "$ref" | sed -n 's/.* measure *\([^ ]*\) above.*/\1/p')
      rcell=$(echo "$ref" | sed -n 's/.*at cell \([0-9]*\).*/\1/p')
      cline=$(echo "$blk" | grep "$rname " | head -n 1)
      cval=$(echo "$cline" | sed -n 's/.*max=\( *[^ ]*\).*/\1/p' | tr -d ' ')
      ccell=$(echo "$cline" | sed -n 's/.*cell=\([0-9]*\).*/\1/p')
      if [ "$rval" = "$cval" ] && [ "$rcell" = "$ccell" ]; then
         echo "PASS refusal_is_the_certification_row_of_$nm measured=$rval@$rcell reference=$cval@$ccell tol=0"
      else
         echo "FAIL refusal_is_the_certification_row_of_$nm measured=$rval@$rcell reference=$cval@$ccell tol=0"
         n_fail=$((n_fail+1))
      fi
   elif [ "$info" = "2" ]; then
      # info = 2 with no refusing row named: allowed only when the solve
      # never claimed its stop test (it aborted instead).
      if grep -qE '\((JFNK|PTC)\) (STAGNATED|no descent|line search found no descent)' "$L"; then
         echo "PASS refusal_is_the_certification_row_of_$nm measured=aborted reference=aborted tol=0"
      else
         echo "FAIL refusal_is_the_certification_row_of_$nm measured=no_row_named reference=a_named_row tol=0"
         n_fail=$((n_fail+1))
      fi
   fi
   # 3. the stop test of the iteration read the certified rows.
   nstop=$(grep -cE '^ \((JFNK|PTC)\) loop-top stop: the acceptance gate is met and every one of the' "$L")
   nread=$(grep -E '^ \((JFNK|PTC)\) loop-top stop: ' "$L" | tail -n 1 \
           | sed -n 's/.*every one of the \([0-9]*\) certified row.*/\1/p')
   if [ "$info" = "0" ]; then
      if [ "$nstop" -ge 1 ] && [ -n "$nread" ] && [ "$nread" -ge 3 ]; then
         echo "PASS stop_test_reads_the_certified_rows_of_$nm measured=${nread}_rows reference=at_least_3 tol=0"
      else
         echo "FAIL stop_test_reads_the_certified_rows_of_$nm measured=${nread:-no_statement} reference=at_least_3 tol=0"
         echo "     an info=0 solve must state that its stop test read the"
         echo "     certified rows of the iterate it stopped on, over at least"
         echo "     the three hydrodynamic rows"
         n_fail=$((n_fail+1))
      fi
   elif [ -n "$ref" ]; then
      if [ "$nstop" -ge 1 ]; then
         echo "PASS stop_test_reads_the_certified_rows_of_$nm measured=refusal_after_a_stop_statement reference=a_stop_statement tol=0"
      else
         echo "FAIL stop_test_reads_the_certified_rows_of_$nm measured=refusal_with_no_stop_statement reference=a_stop_statement tol=0"
         echo "     the state handed back refuses a hydrodynamic row and the"
         echo "     solve made no statement that its stop test read that row:"
         echo "     the iteration stopped for success on the gate alone"
         n_fail=$((n_fail+1))
      fi
   fi
done

if [ $n_fail -gt 0 ]; then exit 1; fi

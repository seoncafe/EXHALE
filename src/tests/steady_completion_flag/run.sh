#!/bin/bash
# The completion flag of a steady solve is a statement about the state it
# hands back.
#
# What is asserted, on the log of a real run (the logs are given in
# EXHALE_STEADY_RUNLOGS, colon separated; each must contain one steady
# solve):
#
#   1. info = 0 only when every hydrodynamic row of the certification of
#      the state handed back is `within` its tolerance, and info = 2 when
#      one of them is `ABOVE` it. The loop-top gate reads the residual of
#      the iterate under the composition of the previous evaluation, so
#      the two can disagree; the flag must follow the state, not the gate.
#   2. The row the solver names when it refuses carries the SAME number
#      and the SAME cell as that row in the certification block: one
#      record, read twice, never two evaluations.
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
   if [ "$nabove" -eq 0 ] && [ "$nunavail" -eq 0 ]; then want=0; else want=2; fi
   if [ "$info" = "$want" ]; then
      echo "PASS flag_follows_the_returned_state_of_$nm measured=$info reference=$want tol=0"
   else
      echo "FAIL flag_follows_the_returned_state_of_$nm measured=$info reference=$want tol=0"
      echo "     rows above tolerance: $nabove, unavailable: $nunavail"
      n_fail=$((n_fail+1))
   fi
   # 2. The refusal the solver prints and the row of the certification.
   ref=$(grep -E 'refusing row: hydrodynamic' "$L" | tail -n 1)
   if [ "$want" = "2" ] && [ -n "$ref" ]; then
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
   elif [ "$want" = "2" ]; then
      # info = 2 with no refusing row named: allowed only when the solve
      # never claimed the gate at the loop top (it aborted instead).
      if grep -qE '\((JFNK|PTC)\) (STAGNATED|no descent|line search found no descent)' "$L"; then
         echo "PASS refusal_is_the_certification_row_of_$nm measured=aborted reference=aborted tol=0"
      else
         echo "FAIL refusal_is_the_certification_row_of_$nm measured=no_row_named reference=a_named_row tol=0"
         n_fail=$((n_fail+1))
      fi
   fi
done

if [ $n_fail -gt 0 ]; then exit 1; fi

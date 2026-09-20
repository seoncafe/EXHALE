#!/bin/bash
# THE RESIDUAL A STEADY SOLVE ACCEPTS ON IS THE RESIDUAL OF THE STATE IT
# HANDS BACK, AT THAT STATE'S OWN COMPOSITION.
#
# The composition is an eliminated variable: R(U) = L(U) + S(U, c*(U)) with
# c*(U) the equilibrium composition OF U. When the elimination is lagged --
# S formed with the composition of the previous iterate -- the Newton drives
# a different system to zero, and the state it returns is not a root of the
# coupled one. The three assertions below are what says the elimination is
# not lagged, taken from the log of a real run:
#
#   1. THE GATE IS THE CERTIFICATION. The ||R|| the solve reports on its
#      `done info=` line and the largest hydrodynamic row of the
#      certification of the state handed back are one measurement of one
#      state, so they agree to 1e-10 relative. (The solver asserts the same
#      equality at 1e-12 internally and refuses the solve if it fails; this
#      reads the two printed numbers, so that the assertion cannot be
#      satisfied by the assertion.)
#   2. THE FLAG FOLLOWS THAT MEASUREMENT: a hydrodynamic row above its own
#      tolerance, or one the certification could not measure, forbids
#      info = 0, so info = 0 states that every certified hydrodynamic row of
#      the state handed back is within its tolerance. The converse is not
#      asserted, and must not be: a solve whose hydrodynamic rows are all
#      within tolerance can still hand back a nonzero flag, info = 1 when it
#      ends on its pass budget with the joint gates unmet and info = 2 when
#      the chemistry at return refuses the state, so a nonzero flag is not by
#      itself a statement about these rows.
#   3. THE HAND-BACK MEASUREMENT AGREES WITH THE ACCEPTED ITERATE'S OWN
#      NUMBER, both read from the log. With the elimination inside the
#      residual evaluation the iterate and the state it becomes are one
#      state, so the ratio is 1 to the accuracy of the composition fixed
#      point; with the elimination lagged the Newton stops where the LAGGED
#      residual is small and the ratio is the size of the lag. The bound is
#      a factor 2. MEASURED on wasp_full_newton: 7.95 with the lagged
#      residual (accepted iterate 6.056e-06, state handed back 4.814e-05).
#
# Assertions 1 and 2 hold on a lagged build too and are not the item's RED:
# they are what must not break. Assertion 3 is the one the lag fails.
#
# ALL THREE READ ONE SOLVE. A run writes the records of every stationary
# solve it makes into one log, so a reader that selects the iteration line,
# the best-iterate restore, the certification block and the completion line
# by four independent searches over the whole file can take them from four
# different solves; the ratio it then forms measures nothing. The scope
# below is the LAST COMPLETED solve, and every record is taken from that
# solve alone:
#   * the solve is named by the last `done info=` line;
#   * a record carrying a `solve=<n>` token belongs to the solve with that
#     n, and nothing else does;
#   * a log written before that token existed is delimited instead: the
#     solve spans from the previous completion line to its own, and the
#     records inside that span are its records. Every archived log is in
#     that form, so this path is not a compatibility detail but the one
#     that reads the record;
#   * the certification block carries no token of its own, and is matched
#     positionally: the solver prints the block and then the completion
#     line, so the block of the chosen solve is the last block whose header
#     stands before that line with no other completion line between them;
#   * a record printed AFTER the last completion line belongs to a solve
#     that did not complete, and a chosen solve with no iterate record or no
#     certification block has no evidence to read. Both print
#     `incomplete_evidence` and fail. A ratio is never built across solves.
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# Usage: EXHALE_STEADY_RUNLOGS=<log>[:<log>...] \
#            src/tests/steady_selfconsistent_residual/run.sh
set -u

LOGS="${EXHALE_STEADY_RUNLOGS:-}"
if [ -z "$LOGS" ] && [ -z "${EXHALE_SEED_PROBE_LOGS:-}" ]; then
   echo "FAIL steady_selfconsistent_residual measured=no_log reference=EXHALE_STEADY_RUNLOGS tol=0"
   echo "     give it the log of a run whose steady solver ran, or"
   echo "     EXHALE_SEED_PROBE_LOGS for assertion 4 alone"
   exit 1
fi

n_fail=0
_logs=()
if [ -n "$LOGS" ]; then IFS=':' read -r -a _logs <<< "$LOGS"; fi
for L in ${_logs[@]+"${_logs[@]}"}; do
   nm=$(basename "$(dirname "$L")")
   if [ ! -f "$L" ]; then
      echo "FAIL selfconsistent_residual_of_$nm measured=no_log reference=$L tol=0"
      n_fail=$((n_fail+1)); continue
   fi
   # THE SCOPE: the records of the last completed stationary solve, and of
   # no other. The scalars come out prefixed `scope `, the certification
   # block after the marker line, so that a block line can never be read as
   # a scalar.
   scoped=$(awk '
      function solve_id(s,   t) {
         t = s
         if (match(t, /(^|[ \t])solve=[0-9]+([ \t]|$)/)) {
            t = substr(t, RSTART, RLENGTH); gsub(/[^0-9]/, "", t); return t
         }
         return ""
      }
      function norm_of(s,   t) {
         t = s
         sub(/^.*\|\|R\|\|=[ \t]*/, "", t)
         sub(/[ \t].*$/, "", t)
         return t
      }
      # A record of the chosen solve stands between the previous completion
      # line and this one, and, where the records carry the token, names
      # this solve. Neither test alone is enough: the span alone would take
      # in a solve that printed records and no completion line of its own,
      # and the token alone would take in a same-numbered solve of another
      # run whose log was appended to this one.
      function in_solve(n) {
         if (n <= lo || n >= d) return 0
         if (fmt == "identified") return (solve_id(l[n]) == id)
         return 1
      }
      { l[NR] = $0 }
      /^ \((JFNK|PTC)\) done info=/             { nd++; dl[nd] = NR }
      /^ \((JFNK|PTC)\) returning best iterate/  { nb++; bl[nb] = NR }
      /^ \((JFNK|PTC)\) it +[0-9]+ +\|\|R\|\|=/  { ni++; il[ni] = NR }
      /\(certification\) \((JFNK|PTC)\) state /  { nc++; cl[nc] = NR }
      END {
         if (nd == 0) { print "scope status=no_steady_solve"; exit }
         d = dl[nd]
         id = solve_id(l[d])
         fmt = (id == "") ? "legacy" : "identified"
         lo = 0
         if (nd > 1) lo = dl[nd-1]
         # A solve whose records stand after the last completion line never
         # completed: its evidence is a fragment.
         if ((nb > 0 && bl[nb] > d) || (ni > 0 && il[ni] > d)) {
            print "scope status=incomplete_evidence"
            print "scope reason=records_after_the_last_completion_line"
            exit
         }
         ib = 0; for (k = 1; k <= nb; k++) if (in_solve(bl[k])) ib = bl[k]
         ii = 0; for (k = 1; k <= ni; k++) if (in_solve(il[k])) ii = il[k]
         if (ib > 0)      { kind = "best_iterate"; iline = ib }
         else if (ii > 0) { kind = "iteration";    iline = ii }
         else {
            print "scope status=incomplete_evidence"
            print "scope reason=no_iterate_record_in_that_solve"
            exit
         }
         cb = 0; for (k = 1; k <= nc; k++) if (cl[k] < d) cb = cl[k]
         if (cb == 0) {
            print "scope status=incomplete_evidence"
            print "scope reason=no_certification_block_before_that_completion_line"
            exit
         }
         for (k = 1; k <= nd; k++) if (dl[k] > cb && dl[k] < d) {
            print "scope status=incomplete_evidence"
            print "scope reason=a_completion_line_stands_between_the_block_and_the_solve"
            exit
         }
         print "scope status=ok"
         print "scope format=" fmt
         print "scope solve=" (fmt == "identified" ? id : "-")
         print "scope done_line=" d
         print "scope iterate_kind=" kind
         print "scope iterate_line=" iline
         print "scope iterate_norm=" norm_of(l[iline])
         print "scope done_text=" l[d]
         print "--- certification block"
         for (i = cb; i <= d; i++) { print l[i]; if (l[i] ~ /CERTIFIED/) break }
      }' "$L")
   status=$(printf '%s\n' "$scoped" | sed -n 's/^scope status=//p')
   if [ "$status" = "no_steady_solve" ] || [ -z "$status" ]; then
      echo "FAIL selfconsistent_residual_of_$nm measured=no_steady_solve reference=one_solve tol=0"
      n_fail=$((n_fail+1)); continue
   fi
   if [ "$status" = "incomplete_evidence" ]; then
      reason=$(printf '%s\n' "$scoped" | sed -n 's/^scope reason=//p')
      echo "FAIL selfconsistent_residual_of_$nm measured=incomplete_evidence reference=one_completed_solve tol=0"
      echo "     $reason"
      echo "     the records of one solve are not all present, and the rows of"
      echo "     this suite are not built from records of different solves"
      n_fail=$((n_fail+1)); continue
   fi
   fmt=$(printf '%s\n' "$scoped" | sed -n 's/^scope format=//p')
   sid=$(printf '%s\n' "$scoped" | sed -n 's/^scope solve=//p')
   it_kind=$(printf '%s\n' "$scoped" | sed -n 's/^scope iterate_kind=//p')
   it_r=$(printf '%s\n' "$scoped" | sed -n 's/^scope iterate_norm=//p')
   done_line=$(printf '%s\n' "$scoped" | sed -n 's/^scope done_text=//p')
   blk=$(printf '%s\n' "$scoped" | sed -n '/^--- certification block$/,$p' | tail -n +2)
   info=$(echo "$done_line" | sed -n 's/.*done info=\([0-9-]*\).*/\1/p')
   gate=$(echo "$done_line" | sed -n 's/.*||R||= *\([^ ]*\).*/\1/p')

   # 1. the gate number and the largest certified hydrodynamic row.
   cert=$(echo "$blk" | grep 'hydrodynamic .* row ' | sed -n 's/.*max=\( *[^ ]*\).*/\1/p' \
          | tr -d ' ' | sort -g | tail -n 1)
   rel=$(awk -v a="$gate" -v b="$cert" 'BEGIN{ if(a+0==0){print (b+0==0)?0:1; exit}
                                               d=(a-b); if(d<0)d=-d; print d/(a<0?-a:a)}')
   ok=$(awk -v r="$rel" 'BEGIN{print (r<=1e-10)?1:0}')
   if [ -n "$cert" ] && [ "$ok" = "1" ]; then
      echo "PASS gate_is_the_certification_of_$nm measured=$rel reference=0 tol=1e-10"
   else
      echo "FAIL gate_is_the_certification_of_$nm measured=${rel:-no_row} reference=0 tol=1e-10"
      echo "     gate ||R||=$gate, largest certified hydrodynamic row=${cert:-none}"
      n_fail=$((n_fail+1))
   fi

   # 2. the flag follows that certification: a row above its tolerance
   #    forbids info = 0. Which nonzero flag it is depends on where the
   #    solve stopped (the pass budget, or the chemistry at return) and is
   #    not this assertion's business.
   rows=$(echo "$blk" | grep 'hydrodynamic .* row ')
   nabove=$(echo "$rows" | grep -c 'ABOVE')
   nunavail=$(echo "$rows" | grep -c 'UNAVAILABLE')
   if [ -z "$info" ]; then
      echo "FAIL flag_of_$nm measured=no_flag reference=a_completion_flag tol=0"
      n_fail=$((n_fail+1))
   elif [ "$nabove" -eq 0 ] && [ "$nunavail" -eq 0 ]; then
      echo "PASS flag_of_$nm measured=info=$info reference=no_hydrodynamic_row_above tol=0"
      if [ "$info" != "0" ]; then
         echo "     every certified hydrodynamic row is within its tolerance;"
         echo "     this flag was not set by one of them"
      fi
   else
      wrow=$(echo "$rows" | grep 'ABOVE' | head -n 1 \
             | sed -n 's/.*\(hydrodynamic [a-z]* row\).*max=\( *[^ ]*\).*/\1 \2/p' | tr -s ' ')
      if [ "$info" != "0" ]; then
         echo "PASS flag_of_$nm measured=info=$info reference=info/=0 tol=0 (${wrow:-unavailable})"
      else
         echo "FAIL flag_of_$nm measured=info=0 reference=info/=0 tol=0"
         echo "     rows above tolerance: $nabove, unavailable: $nunavail (${wrow:-unavailable})"
         n_fail=$((n_fail+1))
      fi
   fi

   # 3. the accepted iterate's number and the state it became.
   #    The accepted iterate's own ||R|| is the last one THAT SOLVE printed
   #    for an iterate (its best-iterate restore line when it has one, else
   #    its last iteration line); the state's own is the one on that solve's
   #    completion line. Read from the log in both cases, so the assertion
   #    does not depend on a line only the self-consistent build prints.
   if [ -z "$it_r" ]; then
      echo "FAIL handback_matches_the_accepted_iterate_of_$nm measured=incomplete_evidence reference=1 tol=2"
      n_fail=$((n_fail+1))
   else
      ratio=$(awk -v a="$gate" -v b="$it_r" 'BEGIN{ if(b+0==0){print 0; exit} print (a+0)/(b+0)}')
      ok=$(awk -v r="$ratio" 'BEGIN{if(r<0)r=-r; print (r<=2.0 && r>=0.5)?1:0}')
      if [ "$ok" = "1" ]; then
         echo "PASS handback_matches_the_accepted_iterate_of_$nm measured=$ratio reference=1 tol=2"
         echo "     read from the $fmt scope of solve $sid, iterate record $it_kind"
      else
         echo "FAIL handback_matches_the_accepted_iterate_of_$nm measured=$ratio reference=1 tol=2"
         echo "     read from the $fmt scope of solve $sid, iterate record $it_kind"
         echo "     the accepted iterate scored $it_r and the state it became scores $gate:"
         echo "     the composition elimination is lagged, so the Newton drove a"
         echo "     different system to zero than the one the state is judged on"
         n_fail=$((n_fail+1))
      fi
   fi
done

# ----------------------------------------------------------------------
# 4. THE RESIDUAL OF ONE STATE IS ONE NUMBER, whatever composition the
#    elimination was seeded with.
#
# R(U) = L(U) + S(U, c*(U)) is a function of U only if the elimination
# reaches one composition c*(U). Where it does not, the state is judged on a
# number that depends on which evaluation produced it, and a Newton that
# accepts at `Resid tol` is accepting noise five decades larger.
#
# The seed perturbation this asserts on is SMALL and in one direction: the
# molecular hydrogen partition scaled by (1 + 1e-6), which is far below the
# difference between the seeds a real solve hands over (the iterate's own
# composition against the previous trial's). A residual that moves by more
# than `Resid tol` under it is not a function of its arguments at the
# accuracy the solve is asked for.
#
# MEASURED on the hot Uranus element state (mol_diffusion reloaded,
# `Coupled carrier solve: True`, `Resid tol 1e-8`, one thread):
#   H2 eliminated  (Molecular carrier transport off): 3.155e-03   FAIL
#   H2 carried     (Molecular carrier transport on):  9.736e-10   PASS
# The elimination cannot determine the molecular hydrogen content of the
# shielded layer -- the fast chemistry conserves H2 nuclei there -- so with
# H2 eliminated the sweep returns the seed's content (0.77 of any seed
# perturbation survives, MEASURED identically at 1e-6 and 1e-2) and the base
# cell's energy row turns that into 3.2e3 times the relative seed change.
#
# The DISPLACED-seed number the run also prints, from a seed that hands every
# cell the composition of a neighbour five cells out, is reported and not
# asserted: at that size the seed crosses a root boundary of the ionization
# front (two roots 19 per cent apart in x(H I) around r = 1.09-1.11), which is
# a question about which root the network is allowed to land on and not about
# the accuracy of one fixed point.
#
# Usage: EXHALE_SEED_PROBE_LOGS=<log>[:<log>...] on a run made with
#        EXHALE_RESID_SC_BASIN=<passes>.
PROBE_LOGS="${EXHALE_SEED_PROBE_LOGS:-}"
if [ -n "$PROBE_LOGS" ]; then
   IFS=':' read -r -a _plogs <<< "$PROBE_LOGS"
   for L in "${_plogs[@]}"; do
      nm=$(basename "$(dirname "$L")")
      if [ ! -f "$L" ]; then
         echo "FAIL residual_is_one_number_for_$nm measured=no_log reference=$L tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      # `Resid tol` of that very run, so no tolerance is written here.
      rtol=$(grep -E 'Residual-based convergence, tol =' "$L" | tail -n 1 \
             | sed -n 's/.*tol = *\([^ ]*\).*/\1/p')
      seeddep=$(grep -E '^ \(sc_basin\) seed x\(H2\)\*\(1\+ *1\.00E-06\)' "$L" \
                | tail -n 1 | sed -n 's/.*||R_d-R_0|| *\([^ ,]*\).*/\1/p')
      if [ -z "$rtol" ] || [ -z "$seeddep" ]; then
         echo "FAIL residual_is_one_number_for_$nm measured=${seeddep:-no_basin_line} reference=${rtol:-no_resid_tol} tol=0"
         echo "     run the case with EXHALE_RESID_SC_BASIN=<passes>"
         n_fail=$((n_fail+1)); continue
      fi
      ok=$(awk -v d="$seeddep" -v t="$rtol" 'BEGIN{print ((d+0)<=(t+0))?1:0}')
      if [ "$ok" = "1" ]; then
         echo "PASS residual_is_one_number_for_$nm measured=$seeddep reference=0 tol=$rtol"
      else
         echo "FAIL residual_is_one_number_for_$nm measured=$seeddep reference=0 tol=$rtol"
         echo "     a 1e-6 relative change of the seed's molecular hydrogen"
         echo "     partition moves the residual of the SAME state by $seeddep,"
         echo "     so the composition elimination does not determine it"
         n_fail=$((n_fail+1))
      fi
      shifted=$(grep -E '^ \(resid_sc_probe\) seed dependence' "$L" | tail -n 1 \
                | sed -n 's/.*Resid tol: *\([^ ]*\).*/\1/p')
      if [ -n "$shifted" ]; then
         echo "     reported, not asserted: displaced-seed dependence $shifted"
      fi
   done
fi

# ------------------------------------------------------------------ #
# 5. THE ELIMINATION REACHES THE TOLERANCE IT IS ASKED FOR (B5j part 2).
#
# The stopping test is on the quantities the residual reads the composition
# for: the particle count of the thermal equation of state, the H2 share of
# it the caloric map weights its ladder by, and the two radiative terms of
# the energy row (increment_of_what_the_residual_reads, steady_newton.f90).
#
# The measure this replaced was the largest relative change of any SPECIES
# FRACTION of any cell with an absolute floor of 1e-12 under its
# denominator, and in the fully ionized wind the He I and H2+ fractions are
# ~1e-20: their denominator is the floor, the quotient is their round-off,
# and it plateaus at 6e-11. MEASURED on that measure, `eq_sweep_reltol` of
# 1e-12 was never met from either seed at pass caps of 200 and 500 (reports
# B5f, B5g section 5.2, B5h section 2): the elimination ran the whole cap
# every time and the passes bought round-off.
#
# So the assertion is the one that measure cannot satisfy: with a tolerance
# below its floor, every evaluation of the run STOPS on the test rather than
# on the cap. `EXHALE_EQ_MEASURE=species` on the same binary is the RED.
#
# Usage: EXHALE_EQ_MEASURE_LOGS=<log>[:<log>...] on runs made with
#        EXHALE_EQ_MEASURE_TRACE=1, EXHALE_RESID_EQ_TOL=<t> and
#        EXHALE_RESID_SC_MAX=<cap>; the tolerance and the cap are read out of
#        the log's own trace so that no number is written here.
EQ_LOGS="${EXHALE_EQ_MEASURE_LOGS:-}"
if [ -n "$EQ_LOGS" ]; then
   IFS=':' read -r -a _elogs <<< "$EQ_LOGS"
   for L in "${_elogs[@]}"; do
      nm=$(basename "$(dirname "$L")")
      if [ ! -f "$L" ]; then
         echo "FAIL eq_measure_stops_on_its_test_for_$nm measured=no_log reference=$L tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      etol="${EXHALE_EQ_MEASURE_TOL:-}"
      ecap="${EXHALE_EQ_MEASURE_CAP:-}"
      if [ -z "$etol" ] || [ -z "$ecap" ]; then
         echo "FAIL eq_measure_stops_on_its_test_for_$nm measured=no_setting reference=EXHALE_EQ_MEASURE_TOL,EXHALE_EQ_MEASURE_CAP tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      # Every traced evaluation is a run of passes numbered from 1; the last
      # pass of each is the one the loop left on. An evaluation that left on
      # the CAP did not reach the tolerance.
      read -r ngrp natcap worst <<< "$(awk -v t="$etol" -v cap="$ecap" '
         /\(eq_measure\) pass/ {
            p = $3 + 0;  v = $5 + 0
            if (have && p <= pp) { grp++; if (pp >= cap) atcap++;
                                   if (pv > worst) worst = pv }
            pp = p;  pv = v;  have = 1
         }
         END{ if (have) { grp++; if (pp >= cap) atcap++;
                          if (pv > worst) worst = pv }
              printf "%d %d %.6e\n", grp, atcap, worst }' "$L")"
      if [ "${ngrp:-0}" -eq 0 ]; then
         echo "FAIL eq_measure_stops_on_its_test_for_$nm measured=no_trace reference=EXHALE_EQ_MEASURE_TRACE=1 tol=0"
         n_fail=$((n_fail+1)); continue
      fi
      if [ "$natcap" -eq 0 ]; then
         echo "PASS eq_measure_stops_on_its_test_for_$nm measured=0_of_$ngrp reference=0 tol=$ecap"
      else
         echo "FAIL eq_measure_stops_on_its_test_for_$nm measured=${natcap}_of_$ngrp reference=0 tol=$ecap"
         echo "     ${natcap} of ${ngrp} evaluations ran the whole pass cap"
         echo "     instead of reaching eq_sweep_reltol = $etol"
         n_fail=$((n_fail+1))
      fi
      ok=$(awk -v w="$worst" -v t="$etol" 'BEGIN{print ((w+0)<=(t+0))?1:0}')
      if [ "$ok" = "1" ]; then
         echo "PASS eq_measure_last_pass_is_within_the_tolerance_for_$nm measured=$worst reference=0 tol=$etol"
      else
         echo "FAIL eq_measure_last_pass_is_within_the_tolerance_for_$nm measured=$worst reference=0 tol=$etol"
         n_fail=$((n_fail+1))
      fi
   done
fi

# ------------------------------------------------------------------ #
# 6. THE READER TAKES THE RECORDS OF ONE SOLVE, on logs whose records are
#    known.
#
# Each fixture under fixtures/ is a written-out log; its numbers are chosen
# to be read and are not physical. The solvable three all carry the same two
# numbers, a final completion at ||R|| = 1.208E-08 and an accepted iterate at
# 1.433E-08, so the row must read 1.208E-08 / 1.433E-08 = 0.842987 on each of
# them; what differs is WHICH record carries the second number and how the
# reader has to find it:
#   multiple_solves_identified        three complete solves, the last with
#                                     its own best-iterate restore;
#   restore_only_in_an_earlier_solve  the only restore of the run is an
#                                     earlier solve's, and the final solve's
#                                     last iteration line is the record;
#   no_restore_in_the_final_solve     the same, written without the solve
#                                     token, so the span between completion
#                                     lines is the only delimiter;
#   truncated_final_solve             a solve's records with no completion
#                                     line of their own: no ratio exists and
#                                     the evidence is incomplete.
# An unscoped reader takes the earlier restore, 4.267E-01, on the middle two
# and reports 2.83E-08.
FIXTURES="$(cd "$(dirname "$0")" && pwd)/fixtures"
if [ -z "${EXHALE_STEADY_FIXTURE_CHILD:-}" ] && [ -d "$FIXTURES" ]; then
   for spec in \
      "multiple_solves_identified|0.842987|identified scope of solve 3, iterate record best_iterate" \
      "restore_only_in_an_earlier_solve|0.842987|identified scope of solve 2, iterate record iteration" \
      "no_restore_in_the_final_solve|0.842987|legacy scope of solve -, iterate record iteration" \
      "truncated_final_solve|incomplete_evidence|records_after_the_last_completion_line"
   do
      fx="${spec%%|*}";  rest="${spec#*|}"
      want="${rest%%|*}";  scope_txt="${rest#*|}"
      out=$(EXHALE_STEADY_FIXTURE_CHILD=1 \
            EXHALE_STEADY_RUNLOGS="$FIXTURES/$fx/run.log" \
            EXHALE_SEED_PROBE_LOGS= EXHALE_EQ_MEASURE_LOGS= \
            bash "$0" 2>&1)
      if [ "$want" = "incomplete_evidence" ]; then
         if printf '%s\n' "$out" | grep -q "measured=incomplete_evidence" \
            && printf '%s\n' "$out" | grep -q -- "$scope_txt"; then
            echo "PASS reader_reports_incomplete_evidence_on_$fx measured=incomplete_evidence reference=incomplete_evidence tol=0"
         else
            echo "FAIL reader_reports_incomplete_evidence_on_$fx measured=a_verdict reference=incomplete_evidence tol=0"
            printf '%s\n' "$out" | sed 's/^/     /'
            n_fail=$((n_fail+1))
         fi
         continue
      fi
      got=$(printf '%s\n' "$out" \
            | sed -n "s/^PASS handback_matches_the_accepted_iterate_of_$fx measured=\([^ ]*\).*/\1/p")
      scope_seen=$(printf '%s\n' "$out" | grep -c -- "$scope_txt")
      if [ -z "$got" ]; then
         echo "FAIL reader_reads_one_solve_of_$fx measured=no_row reference=$want tol=1e-5"
         printf '%s\n' "$out" | sed 's/^/     /'
         n_fail=$((n_fail+1)); continue
      fi
      ok=$(awk -v a="$got" -v b="$want" 'BEGIN{ d=(a-b); if(d<0)d=-d;
                                                print (d/(b+0)<=1e-5)?1:0}')
      if [ "$ok" = "1" ] && [ "$scope_seen" -ge 1 ]; then
         echo "PASS reader_reads_one_solve_of_$fx measured=$got reference=$want tol=1e-5"
      else
         echo "FAIL reader_reads_one_solve_of_$fx measured=$got reference=$want tol=1e-5"
         [ "$scope_seen" -ge 1 ] || echo "     the records came from another scope than: $scope_txt"
         n_fail=$((n_fail+1))
      fi
   done
fi

if [ $n_fail -gt 0 ]; then exit 1; fi

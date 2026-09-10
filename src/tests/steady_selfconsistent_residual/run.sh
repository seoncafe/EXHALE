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
#   2. THE FLAG FOLLOWS THAT MEASUREMENT: info = 0 only when every
#      hydrodynamic row of that certification is within its tolerance, and
#      info = 2 with the refusing row and its number otherwise.
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
   blk=$(awk '/\(certification\) \((JFNK|PTC)\) state /{n=NR} {l[NR]=$0} \
              END{if(n) for(i=n;i<=NR;i++){print l[i]; if(l[i] ~ /CERTIFIED/) exit}}' "$L")
   done_line=$(grep -E '^ \((JFNK|PTC)\) done info=' "$L" | tail -n 1)
   if [ -z "$blk" ] || [ -z "$done_line" ]; then
      echo "FAIL selfconsistent_residual_of_$nm measured=no_steady_solve reference=one_solve tol=0"
      n_fail=$((n_fail+1)); continue
   fi
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

   # 2. the flag follows that certification.
   rows=$(echo "$blk" | grep 'hydrodynamic .* row ')
   nabove=$(echo "$rows" | grep -c 'ABOVE')
   nunavail=$(echo "$rows" | grep -c 'UNAVAILABLE')
   if [ "$nabove" -eq 0 ] && [ "$nunavail" -eq 0 ]; then want=0; else want=2; fi
   if [ "$info" = "$want" ]; then
      if [ "$want" = "0" ]; then
         echo "PASS flag_of_$nm measured=info=0 reference=every_row_within tol=0"
      else
         wrow=$(echo "$rows" | grep 'ABOVE' | head -n 1 \
                | sed -n 's/.*\(hydrodynamic [a-z]* row\).*max=\( *[^ ]*\).*/\1 \2/p' | tr -s ' ')
         echo "PASS flag_of_$nm measured=info=2 reference=a_row_above tol=0 (${wrow:-unavailable})"
      fi
   else
      echo "FAIL flag_of_$nm measured=info=$info reference=info=$want tol=0"
      echo "     rows above tolerance: $nabove, unavailable: $nunavail"
      n_fail=$((n_fail+1))
   fi

   # 3. the accepted iterate's number and the state it became.
   #    The accepted iterate's own ||R|| is the last one the solve printed
   #    for an iterate (the best-iterate restore line when there is one, else
   #    the last iteration line); the state's own is the one on the done
   #    line. Read from the log in both cases, so the assertion does not
   #    depend on a line only the self-consistent build prints.
   it_r=$(grep -E '^ \((JFNK|PTC)\) returning best iterate' "$L" | tail -n 1 \
          | sed -n 's/.*||R||= *\([^ ]*\).*/\1/p')
   if [ -z "$it_r" ]; then
      it_r=$(grep -E '^ \((JFNK|PTC)\) it +[0-9]+ +\|\|R\|\|=' "$L" | tail -n 1 \
             | sed -n 's/.*||R||= *\([^ ]*\).*/\1/p')
   fi
   if [ -z "$it_r" ]; then
      echo "FAIL handback_matches_the_accepted_iterate_of_$nm measured=no_iterate_norm reference=1 tol=2"
      n_fail=$((n_fail+1))
   else
      ratio=$(awk -v a="$gate" -v b="$it_r" 'BEGIN{ if(b+0==0){print 0; exit} print (a+0)/(b+0)}')
      ok=$(awk -v r="$ratio" 'BEGIN{if(r<0)r=-r; print (r<=2.0 && r>=0.5)?1:0}')
      if [ "$ok" = "1" ]; then
         echo "PASS handback_matches_the_accepted_iterate_of_$nm measured=$ratio reference=1 tol=2"
      else
         echo "FAIL handback_matches_the_accepted_iterate_of_$nm measured=$ratio reference=1 tol=2"
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

if [ $n_fail -gt 0 ]; then exit 1; fi

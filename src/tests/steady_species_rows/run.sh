#!/bin/bash
# THE STATIONARY SYSTEM CARRIES ONE ROW FOR EVERY TRANSPORTED BALANCE
# (docs/PLAN_20260906_rev2.md B5; docs/b1_target_system_20260906.md T1.7).
#
# Two kinds of row.
#   * The Fortran driver links the PRODUCTION objects and asserts the row
#     registry of steady_newton against the transported set carrier_set_init
#     fixes, configuration by configuration: the row count, the unknown
#     count per cell, the band half-width, the coloring stride, which
#     carrier each row names, the f_sp column each carrier names, and the
#     element headroom of each carrier.  See the header of
#     steady_species_rows_tests.f90.
#   * The log rows state what only a run can show: a solve that entered the
#     JFNK with species rows reports the geometry the formula gives, its
#     Jacobian action matches a directional difference of the same map, and
#     a species row that is out of tolerance is named by the completion
#     flag rather than passed over.  Give them the log of such a run
#     through EXHALE_SPECIES_ROW_LOGS.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Exit status is nonzero if any row fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped.
# EXHALE_TEST_OUT selects where the test executable is written.
# EXHALE_SPECIES_ROW_LOGS is a colon-separated list of run logs; a run whose
# log carries a "(species_jac_test)" line was made with
# EXHALE_SPECIES_JAC_TEST=1.
# EXHALE_CARRIER_UNKNOWN_PAIR gives two run directories of one case that
# differ only in EXHALE_CARRIER_LOG_UNKNOWN, and
# EXHALE_CARRIER_FLOOR_LOGS a list of coupled-solve logs whose adopted
# carriers must stay above the floor of their own unknown space, and
# EXHALE_CARRIER_BUDGET_LOGS a list of them whose adopted carriers must stay
# inside the element budget of their own cell.
# EXHALE_INVENTORY_GOLDEN_DIR is the directory of golden case outputs the
# elemental map classifies (default backup/regression/golden);
# EXHALE_INVENTORY_COUPLED_STATES is a colon-separated list of run
# directories whose output/Ion_species.txt is a state a JFNK with carrier
# rows returned (no default: those rows are skipped without it); and
# EXHALE_INVENTORY_SHARED=0 is the control that judges a state by
# the separate box sides instead of the shared element constraint, and
# EXHALE_ELEMENT_CONSTRAINT_ROWS=0 the one that carries the shared element
# budget as a coordinate face of the unknown box instead of as a row of the
# step.
#
# Usage: src/tests/steady_species_rows/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.  The reload rows at the end
# of this suite run only when a binary was asked for.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" steady_species_rows EXHALE_SPECIES_EXE
if [ "$EXHALE_RUN_EXE_REQUESTED" = "1" ]; then exhale_announce_exe; fi
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/steady_species_rows}"
# The reference states the elemental map classifies (N4a): the golden
# Ion_species.txt of the molecular regression cases, read through their own
# schema headers.
export EXHALE_INVENTORY_GOLDEN_DIR="${EXHALE_INVENTORY_GOLDEN_DIR:-$ROOT/backup/regression/golden}"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"

mkdir -p "$OUT"

FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

rc=0

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL steady_species_rows_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL steady_species_rows_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

# src/tests/test_columns.f90 builds the synthetic column whose species
# carry their own density, shared with the other driver that measures an
# operator against the mass it was handed; it is compiled here because it
# belongs to no production object.
rm -f "$OUT/steady_species_rows_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/steady_species_rows_tests.x" \
    "$ROOT/src/tests/test_columns.f90" \
    "$HERE/steady_species_rows_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL steady_species_rows_build measured=compile_error reference=ok tol=0"
   exit 1; }

"$OUT/steady_species_rows_tests.x" || rc=1

# ------------------------------------------------------------------ #
# THE SHARED BUDGET AS A ROW OF THE STEP IS RED WITHOUT IT AND GREEN WITH
# IT (N4b).
#
# A carrier's own ceiling is one corner of a half-space two carriers of one
# element share, so the budget is a row of the step and not a face of the
# unknown box: the step is projected onto the rows that are active, a
# candidate that violates the nonlinear constraint is restored, and neither
# exists when the budget is a coordinate face.  Running the same driver with
# EXHALE_ELEMENT_CONSTRAINT_ROWS=0 puts it back under a face and forms no
# row, which is the box the B5j to B5l measurements were made in.  The rows
# named below are the ones whose reference distinguishes the two, so they
# must fail in that pass; a suite in which they pass either way would be
# asserting nothing.
# ------------------------------------------------------------------ #
RED_OUT="$OUT/element_rows_off.log"
EXHALE_ELEMENT_CONSTRAINT_ROWS=0 "$OUT/steady_species_rows_tests.x" \
   > "$RED_OUT" 2>&1 || true
nred=0
for row in element_budget_is_a_row_of_the_step \
           species_box_upper_face_is_the_cell_density \
           element_row_budget_is_the_operator_budget \
           the_shared_row_of_the_cell_is_active \
           restoration_measures_the_violation_before \
           restoration_moves_every_infeasible_cell ; do
   if grep -q "^FAIL $row " "$RED_OUT"; then
      nred=$((nred+1))
   else
      echo "     $row does not fail with the budget as a coordinate face"
   fi
done
if [ "$nred" -eq 6 ]; then
   echo "PASS element_constraint_rows_are_red_without_them measured=$nred reference=6 tol=0"
else
   echo "FAIL element_constraint_rows_are_red_without_them measured=$nred reference=6 tol=0"
   rc=1
fi

# ------------------------------------------------------------------ #
# THE SHARED ELEMENT CONSTRAINT IS RED WITHOUT IT AND GREEN WITH IT (N4a).
#
# The rows of shared_element_constraint_is_the_feasible_set are about a
# constraint the separate box sides of the carriers do not state: two
# carriers of one element spend one budget, so a state on both their
# ceilings holds two hydrogen nuclei where the cell has one.  Running the
# same driver with EXHALE_INVENTORY_SHARED=0 judges a state by those
# separate sides, which is what the shared constraint is measured
# against, and the three named below must then call every one of the
# infeasible states feasible; a suite in which they pass either way would
# be asserting nothing.
# ------------------------------------------------------------------ #
SHARED_OUT="$OUT/shared_constraint_off.log"
EXHALE_INVENTORY_SHARED=0 "$OUT/steady_species_rows_tests.x" \
   > "$SHARED_OUT" 2>&1 || true
nred=0
for row in inventory_separate_ceilings_are_infeasible \
           inventory_negative_remainder_is_infeasible \
           inventory_shared_oxygen_is_infeasible ; do
   if grep -q "^FAIL $row " "$SHARED_OUT"; then
      nred=$((nred+1))
   else
      echo "     $row does not fail without the shared constraint"
   fi
done
if [ "$nred" -eq 3 ]; then
   echo "PASS inventory_shared_rows_are_red_without_the_constraint measured=$nred reference=3 tol=0"
else
   echo "FAIL inventory_shared_rows_are_red_without_the_constraint measured=$nred reference=3 tol=0"
   rc=1
fi

# ------------------------------------------------------------------ #
# A SELF-TEST DOES NOT CHANGE THE RUN IT MEASURES (B5g).
#
# EXHALE_SPECIES_JAC_TEST=1 evaluates the residual at Y and at Y + eps r.
# Those evaluations write module state -- the elemental row scales, the
# carrier row and column scales, the carrier measures, the frozen
# background -- that the outer iteration reads back, so before B5g the
# self-test moved the solve it was measuring: MEASURED on the hot Uranus
# reload with one element row, the tested and the untested run agreed to
# every printed digit for five outer iterations, differed at iteration 6
# (||R|| 1.957 against 1.956) and were on different paths by iteration 13,
# one converging and one stalling.  A probe evaluation now puts back what
# it overwrote, so the two runs are the same run.
#
# EXHALE_SPECIES_SELFTEST_PAIR=<tested dir>:<untested dir> gives the two run
# directories of ONE case that differ in nothing but the environment
# variable.  Both rows are exact: same iteration lines, same output files.
# ------------------------------------------------------------------ #
PAIR="${EXHALE_SPECIES_SELFTEST_PAIR:-}"
if [ -n "$PAIR" ]; then
   IFS=':' read -r DT DU <<< "$PAIR"
   if [ ! -f "$DT/run.log" ] || [ ! -f "$DU/run.log" ]; then
      echo "FAIL species_rows_self_test_pair measured=missing reference=$PAIR tol=0"
      rc=1
   elif ! grep -q '(species_jac_test)' "$DT/run.log"; then
      echo "FAIL species_rows_self_test_pair measured=no_self_test reference=$DT tol=0"
      rc=1
   elif grep -q '(species_jac_test)' "$DU/run.log"; then
      echo "FAIL species_rows_self_test_pair measured=self_test_in_both reference=$DU tol=0"
      rc=1
   else
      grep '(JFNK) it' "$DT/run.log" > "$OUT/selftest_on.it"  || true
      grep '(JFNK) it' "$DU/run.log" > "$OUT/selftest_off.it" || true
      nit=$(wc -l < "$OUT/selftest_off.it")
      ndiff=$(diff "$OUT/selftest_on.it" "$OUT/selftest_off.it" | grep -c '^[<>]' || true)
      if [ "$nit" -lt 2 ]; then
         echo "FAIL species_rows_self_test_changes_no_iteration measured=$nit reference=>=2 tol=0"
         rc=1
      elif [ "$ndiff" -eq 0 ]; then
         echo "PASS species_rows_self_test_changes_no_iteration measured=0 reference=0 tol=0"
      else
         echo "FAIL species_rows_self_test_changes_no_iteration measured=$ndiff reference=0 tol=0"
         diff "$OUT/selftest_on.it" "$OUT/selftest_off.it" | head -n 4 | sed 's/^/     /'
         rc=1
      fi
      # The provenance header carries the wall-clock time of the run, which
      # differs between any two runs of anything; everything else in the
      # file, header and columns, has to be the same bytes.
      nbad=0
      for f in Hydro_ioniz.txt Ion_species.txt; do
         if ! diff -q <(grep -v '^# provenance: git=' "$DT/output/$f") \
                      <(grep -v '^# provenance: git=' "$DU/output/$f") \
                 >/dev/null 2>&1; then
            nbad=$((nbad+1))
            echo "     $f differs"
         fi
      done
      if [ "$nbad" -eq 0 ]; then
         echo "PASS species_rows_self_test_changes_no_output measured=0 reference=0 tol=0"
      else
         echo "FAIL species_rows_self_test_changes_no_output measured=$nbad reference=0 tol=0"
         rc=1
      fi
   fi
fi

# ------------------------------------------------------------------ #
# THE LOG ROWS.
# ------------------------------------------------------------------ #
# The geometry rows below want a self-test log, and the rows AFTER them do
# not: a solve that carried a species row reports its bounds, its unknown
# space and its floor whether or not EXHALE_SPECIES_JAC_TEST was set.  So an
# absent EXHALE_SPECIES_ROW_LOGS skips this block and nothing else -- it used
# to end the script, which left the bound rows unreachable without a variable
# they ask nothing of.
LOGS="${EXHALE_SPECIES_ROW_LOGS:-}"
if [ -z "$LOGS" ]; then
   echo "  geometry log rows skipped: set EXHALE_SPECIES_ROW_LOGS to the"
   echo "  log(s) of a run that entered the JFNK with EXHALE_SPECIES_JAC_TEST=1"
   LOGARR=()
else
   IFS=':' read -r -a LOGARR <<< "$LOGS"
fi
for LOG in "${LOGARR[@]+"${LOGARR[@]}"}"; do
   TAG="$(basename "$(dirname "$LOG")")"
   if [ ! -f "$LOG" ]; then
      echo "FAIL species_rows_log_$TAG measured=missing reference=$LOG tol=0"
      rc=1; continue
   fi

   LINE="$(grep -m1 '(species_jac_test) nvar' "$LOG" || true)"
   if [ -z "$LINE" ]; then
      echo "FAIL species_rows_geometry_$TAG measured=no_line reference=species_jac_test tol=0"
      rc=1; continue
   fi
   NVAR=$(echo "$LINE"  | sed -n 's/.*nvar = \([0-9]*\).*/\1/p')
   NROW=$(echo "$LINE"  | sed -n 's/.*species rows = \([0-9]*\).*/\1/p')
   KL=$(echo "$LINE"    | sed -n 's/.*kl = ku = \([0-9]*\).*/\1/p')
   NCOL=$(echo "$LINE"  | sed -n 's/.*colors = \([0-9]*\).*/\1/p')

   # nvar = 3 + one unknown per transported balance.
   if [ "$NVAR" -eq $((3 + NROW)) ]; then
      echo "PASS species_rows_unknown_count_$TAG measured=$NVAR reference=$((3 + NROW)) tol=0"
   else
      echo "FAIL species_rows_unknown_count_$TAG measured=$NVAR reference=$((3 + NROW)) tol=0"
      rc=1
   fi
   # |row - col| <= 2*nvar (the WENO3 reach) + nvar - 1 (the cell itself).
   if [ "$KL" -eq $((3 * NVAR - 1)) ]; then
      echo "PASS species_rows_bandwidth_$TAG measured=$KL reference=$((3 * NVAR - 1)) tol=0"
   else
      echo "FAIL species_rows_bandwidth_$TAG measured=$KL reference=$((3 * NVAR - 1)) tol=0"
      rc=1
   fi
   # the coloring stride is kl + ku + 1, so same-color columns have disjoint
   # row supports.
   if [ "$NCOL" -eq $((2 * KL + 1)) ]; then
      echo "PASS species_rows_coloring_$TAG measured=$NCOL reference=$((2 * KL + 1)) tol=0"
   else
      echo "FAIL species_rows_coloring_$TAG measured=$NCOL reference=$((2 * KL + 1)) tol=0"
      rc=1
   fi
   # the geometry the solve reports at its end is the one it was built with
   COSTV="$(grep -m1 '(JFNK) cost:' "$LOG" | sed -n 's/.*cost: \([0-9]*\) unknowns.*/\1/p' || true)"
   COSTC="$(grep -m1 '(JFNK) cost:' "$LOG" | sed -n 's/.*, \([0-9]*\) colors.*/\1/p' || true)"
   if [ -n "$COSTV" ]; then
      if [ "$COSTV" -eq "$NVAR" ] && [ "$COSTC" -eq "$NCOL" ]; then
         echo "PASS species_rows_cost_geometry_$TAG measured=${COSTV}x${COSTC} reference=${NVAR}x${NCOL} tol=0"
      else
         echo "FAIL species_rows_cost_geometry_$TAG measured=${COSTV}x${COSTC} reference=${NVAR}x${NCOL} tol=0"
         rc=1
      fi
   fi

   # THE JACOBIAN ACTION.  J r against [F(Y + eps r) - F(Y)]/eps of the same
   # map, species rows measured on their own so a mismatch confined to them
   # cannot be averaged away by the hydrodynamic rows.  What a WRONG band
   # layout or coloring gives is an O(1) mismatch (MEASURED, 1.0), so the
   # bound is what separates the layout from the physics and not a
   # convergence criterion.
   #
   # THE BOUND IS 1e-2 AND NOT 1e-6, and an element row is why.  A carrier
   # unknown is local: MEASURED with the H2 and H+ rows, all rows 7.19e-07
   # and species rows 4.60e-05.  An ELEMENT unknown is not: perturbing the
   # helium mass fraction of one cell changes that cell's opacity and with
   # it the radiation field of every cell above it, which is a coupling the
   # banded Jacobian does not carry by construction and the matrix-free
   # Krylov action does.  MEASURED with the He/H row, all rows 4.92e-03 and
   # species rows 1.29e-03 -- two decades from an O(1) mismatch and three
   # decades above the local case.
   AERR="$(grep -m1 '(species_jac_test) ||J r' "$LOG" | sed -n 's/.*all rows *\([0-9.E+-]*\).*/\1/p' || true)"
   SERR="$(grep -m1 '(species_jac_test) ||J r' "$LOG" | sed -n 's/.*species rows *\([0-9.E+-]*\).*/\1/p' || true)"
   NUNR="$(grep -m1 '(species_jac_test) ||J r' "$LOG" | sed -n 's/.*unresolved colors *\([0-9]*\).*/\1/p' || true)"
   if [ -z "$AERR" ]; then
      echo "FAIL species_rows_jacobian_action_$TAG measured=no_line reference=species_jac_test tol=0"
      rc=1
   else
      awk -v a="$AERR" -v s="$SERR" -v t="$TAG" 'BEGIN{
         if (a+0 < 1.0e-2) printf "PASS species_rows_jacobian_all_%s measured=%s reference=0 tol=1e-2\n", t, a;
         else { printf "FAIL species_rows_jacobian_all_%s measured=%s reference=0 tol=1e-2\n", t, a; exit 1 }
      }' || rc=1
      awk -v s="$SERR" -v t="$TAG" 'BEGIN{
         if (s+0 < 1.0e-2) printf "PASS species_rows_jacobian_species_%s measured=%s reference=0 tol=1e-2\n", t, s;
         else { printf "FAIL species_rows_jacobian_species_%s measured=%s reference=0 tol=1e-2\n", t, s; exit 1 }
      }' || rc=1
   fi
   # A COLOR THAT COULD NOT BE PROBED LEAVES ITS COLUMNS AT ZERO, so a
   # Jacobian with every color unresolved is no Newton model at all and the
   # solve returns the state it started from.
   if [ -n "$NUNR" ]; then
      if [ "$NUNR" -lt "$NCOL" ]; then
         echo "PASS species_rows_colors_resolved_$TAG measured=$NUNR reference=<$NCOL tol=0"
      else
         echo "FAIL species_rows_colors_resolved_$TAG measured=$NUNR reference=<$NCOL tol=0"
         rc=1
      fi
   fi

   # THE COMPLETION FLAG READS THE SPECIES ROWS TOO.  A solve whose species
   # row is above its tolerance may not end info = 0, and the refusal names
   # the row.
   if grep -q '(JFNK) done info=0' "$LOG"; then
      if grep -q 'gate NOT met: carrier row' "$LOG"; then
         echo "FAIL species_rows_flag_reads_species_$TAG measured=info0_with_refused_row reference=info2 tol=0"
         rc=1
      else
         echo "PASS species_rows_flag_reads_species_$TAG measured=info0_no_refused_row reference=info0 tol=0"
      fi
   else
      echo "PASS species_rows_flag_reads_species_$TAG measured=info2 reference=info2_or_clean_info0 tol=0"
   fi

   # THE BANDED MODEL STATES SOMETHING ABOUT EVERY UNKNOWN (B5d).
   #
   # The band is the preconditioner of the Krylov cycle, the Cauchy
   # direction of the trust region and the matrix the damped Gauss-Newton
   # escape factors.  A row with no entry inside it is an equation that
   # responds to no unknown, and a column with none an unknown no equation
   # responds to; either makes the banded system singular, and a singular
   # model cannot produce a descent direction however the step is
   # globalized.  The solve states this itself, once per registry, and it is
   # a structural property so the row is exact and carries no tolerance.
   #
   # RED, MEASURED (hot Uranus reload, one element row): 1 empty row, the
   # equation of the base helium unknown, which element_transport_residual
   # posed as a zero because cell 1 is the Dirichlet reservoir of the
   # element operator and carries no transport balance.  The banded LU then
   # grew from a largest entry of 5.229e+03 to a pivot of 8.413e+19, and the
   # damped Gauss-Newton escape reported ||grad merit|| = 1.10e+44 and
   # aborted at outer iteration 8.
   #
   # Both spellings are counted: the solve's own sentence, and the count line
   # of the EXHALE_BAND_REPORT=1 report, which carries the same fact.
   n=$(grep -cE 'the banded model has a structurally empty|\[band\] structurally empty rows [1-9]|\[band\] structurally empty rows [0-9]+, empty columns [1-9]' "$LOG" || true)
   if [ "$n" -eq 0 ]; then
      echo "PASS species_rows_band_states_every_unknown_$TAG measured=0 reference=0 tol=0"
   else
      echo "FAIL species_rows_band_states_every_unknown_$TAG measured=$n reference=0 tol=0"
      grep -m2 'the banded model has a structurally empty' "$LOG" | sed 's/^/     /'
      rc=1
   fi

   # AND THE MERIT GRADIENT IS A NUMBER THE MODEL'S OWN ENTRIES CAN PRODUCE.
   #
   # The damped Gauss-Newton escape prints ||grad merit|| = ||A^T s|| of the
   # SCALED banded model whenever it finds nothing downhill.  A row whose
   # scale has collapsed below the size on which its own row responds turns
   # an ordinary derivative into an enormous entry of A, and the gradient
   # reports it: RED, MEASURED on the hot Uranus reload with the element row
   # scaled by the operator's own floor, 1.104E+44, from a band whose largest
   # entry was 4.890E+22 for a dF/dY of 46.7.  GREEN, MEASURED with the row
   # scale floored at the flow-time rate of the transported quantity,
   # 4.813E+03 from a band whose largest entry stayed at 5.2E+03.  The bound
   # sits fourteen decades above the second and fourteen below the first; it
   # is a discriminator between a scaled model and a manufactured one, not a
   # convergence criterion, and the solve may still report no descent under
   # it -- that is a statement about the state.
   BIG=0
   while read -r g; do
      awk -v g="$g" 'BEGIN{ if (g+0 > 1.0e30) exit 1 }' || BIG=$((BIG+1))
   done < <(grep 'no descent along the Newton' "$LOG" |
            sed -n 's/.*||grad merit||= *\([0-9.E+-]*\).*/\1/p')
   if [ "$BIG" -eq 0 ]; then
      echo "PASS species_rows_merit_gradient_is_scaled_$TAG measured=$BIG reference=0 tol=1e30"
   else
      echo "FAIL species_rows_merit_gradient_is_scaled_$TAG measured=$BIG reference=0 tol=1e30"
      grep -m2 'no descent along the Newton' "$LOG" | sed 's/^/     /'
      rc=1
   fi

   # THE STATE HANDED BACK IS THE BEST ONE ON THE MEASURE IT IS JUDGED BY
   # (B5e).
   #
   # The merit the step control descends on is a 2-norm on the Newton's own
   # column scales; the state is judged by a maximum over cells of each row
   # against its own largest term, and by every SPECIES row the registry
   # carried as well as the three hydrodynamic ones.  A best-iterate ledger
   # ranking on ||R|| covers three of those rows: it can keep the iterate
   # whose element or carrier balance is the worse of two and hand it back
   # to be refused on exactly that row.  The ledger now ranks on
   # distance_from_certification, each row over the tolerance it is judged
   # against.
   #
   # The row is exact and needs no tolerance: the number the ledger reports
   # at the end IS one of the numbers of the iterations, so it is the minimum of
   # them, and it would not be under a ranking on a different measure.
   BEST="$(grep -m1 '(JFNK) best judged iterate: distance' "$LOG" |
           sed -n 's/.*distance *\([0-9.E+-]*\).*/\1/p' || true)"
   n=$(grep -c '(JFNK) judged rows:' "$LOG" || true)
   if [ -z "$BEST" ] || [ "$n" -eq 0 ]; then
      echo "FAIL species_rows_best_iterate_is_the_best_judged_$TAG measured=no_line reference=judged_rows tol=0"
      rc=1
   else
      MINSEEN=$(grep '(JFNK) judged rows:' "$LOG" |
                sed -n 's/.*distance *\([0-9.E+-]*\).*/\1/p' |
                awk 'NR==1||$1<m{m=$1} END{print m}')
      awk -v b="$BEST" -v m="$MINSEEN" -v t="$TAG" 'BEGIN{
         if (b+0 <= m+0) printf "PASS species_rows_best_iterate_is_the_best_judged_%s measured=%s reference=<=%s tol=0\n", t, b, m;
         else { printf "FAIL species_rows_best_iterate_is_the_best_judged_%s measured=%s reference=<=%s tol=0\n", t, b, m; exit 1 }
      }' || rc=1
   fi

   # AND A TRUST-REGION ITERATION THAT TAKES NO STEP SAYS WHY (B5e).
   #
   # trust_region_step has eleven named exits and the summary line used to
   # call every one of them "model refused by the ray test".  An
   # inadmissible state on the ray, a model that over-predicts and a finite
   # difference below its own cancellation floor are three different
   # failures asking for three different repairs, so every refused iteration
   # names its own.
   NREF=$(grep -c 'accepted=F' "$LOG" || true)
   if [ "$NREF" -gt 0 ]; then
      NWHY=$(grep -c '\[TR\]      no step:' "$LOG" || true)
      if [ "$NWHY" -eq "$NREF" ]; then
         echo "PASS species_rows_tr_refusal_names_its_reason_$TAG measured=$NWHY reference=$NREF tol=0"
      else
         echo "FAIL species_rows_tr_refusal_names_its_reason_$TAG measured=$NWHY reference=$NREF tol=0"
         rc=1
      fi
   fi
done

# ------------------------------------------------------------------ #
# A SPECIES UNKNOWN NEVER LEAVES THE ADMISSIBLE SET (B5j).
#
# A species density and a species mass fraction are non-negative, and a
# carrier density is at most the cell's own.  The trust region writes its
# trial onto the faces of that box before it evaluates it
# (species_unknowns_outside_their_bounds), the probe step is cut to the same
# box, and the screen of eval_residual refuses a state outside it, so no
# state the solve ADOPTS may carry a species unknown below zero.  The solve
# reports the smallest one any adopted iterate carried, and the row is exact:
# zero is a bound, not a tolerance.
#
# RED, MEASURED before the projection existed: the negative carrier of a
# Krylov sample reached f_sp unclamped, and 527 residual samples of one
# `mol_carrier` solve were refused for a negative species unknown while the
# dogleg cut its trials back by 2^-25 to 2^-49 against the same one unknown
# (reports B5g sections 6.3 and 8, B5j section 3).
#
# EXHALE_SPECIES_BOUND_LOGS is its own variable, and not the geometry rows'
# one, because this row asks nothing of the self-test: any solve that carried
# a species row reports the number.
# ------------------------------------------------------------------ #
BOUND_LOGS="${EXHALE_SPECIES_BOUND_LOGS:-}"
if [ -n "$BOUND_LOGS" ]; then
   IFS=':' read -r -a BLOGARR <<< "$BOUND_LOGS"
   for BLOG in "${BLOGARR[@]}"; do
      BTAG="$(basename "$(dirname "$BLOG")")"
      if [ ! -f "$BLOG" ]; then
         echo "FAIL species_rows_adopted_unknown_is_admissible_$BTAG measured=missing reference=$BLOG tol=0"
         rc=1; continue
      fi
      # The end-of-solve statement where the solve completed, and the
      # one of each iteration where it was stopped: the two are the same running
      # minimum, so the last of either will do.
      MINSP="$(grep 'smallest species unknown any adopted iterate carried' "$BLOG" |
               tail -n 1 | sed -n 's/.*carried: *\([0-9.E+-]*\).*/\1/p' || true)"
      if [ -z "$MINSP" ]; then
         MINSP="$(grep '(JFNK) species unknowns: smallest adopted' "$BLOG" |
                  tail -n 1 | sed -n 's/.*smallest adopted *\([0-9.E+-]*\).*/\1/p' || true)"
      fi
      if [ -z "$MINSP" ]; then
         echo "FAIL species_rows_adopted_unknown_is_admissible_$BTAG measured=no_line reference=a_solve_with_a_species_row tol=0"
         rc=1; continue
      fi
      awk -v v="$MINSP" -v t="$BTAG" 'BEGIN{
         if (v+0 >= 0.0) printf "PASS species_rows_adopted_unknown_is_admissible_%s measured=%s reference=>=0 tol=0\n", t, v;
         else { printf "FAIL species_rows_adopted_unknown_is_admissible_%s measured=%s reference=>=0 tol=0\n", t, v; exit 1 }
      }' || rc=1
   done
fi

# ------------------------------------------------------------------ #
# THE CARRIER UNKNOWN OF THE COUPLED SOLVE IS ln n (B5k).
#
# The molecular carrier densities of the coupled steady solve are carried as
# their logarithms; EXHALE_CARRIER_LOG_UNKNOWN=0 restores the density unknown,
# which exists to be measured against (user decision, 2026-09-08, on the
# measurement of item B5j).
#
# Which space a solve was in is a line of its own report, so the pair of runs
# that differ in nothing but that variable states both halves: the default is
# the logarithm, and the switch reaches the density.
#
# EXHALE_CARRIER_UNKNOWN_PAIR=<default run dir>:<EXHALE_CARRIER_LOG_UNKNOWN=0
# run dir>.
# ------------------------------------------------------------------ #
CPAIR="${EXHALE_CARRIER_UNKNOWN_PAIR:-}"
if [ -n "$CPAIR" ]; then
   IFS=':' read -r CD CS <<< "$CPAIR"
   for f in "$CD/run.log" "$CS/run.log"; do
      if [ ! -f "$f" ]; then
         echo "FAIL carrier_unknown_space_pair measured=missing reference=$f tol=0"
         rc=1
      fi
   done
   if [ -f "$CD/run.log" ] && [ -f "$CS/run.log" ]; then
      LD="$(grep -m1 '(JFNK) species unknown space:' "$CD/run.log" || true)"
      LS="$(grep -m1 '(JFNK) species unknown space:' "$CS/run.log" || true)"
      if [ -z "$LD" ] || [ -z "$LS" ]; then
         echo "FAIL carrier_unknown_space_named measured=no_line reference=both_runs tol=0"
         rc=1
      else
         case "$LD" in
            *"a carrier is ln n"*)
               echo "PASS carrier_unknown_default_is_the_logarithm measured=ln_n reference=ln_n tol=0" ;;
            *) echo "FAIL carrier_unknown_default_is_the_logarithm measured=not_ln_n reference=ln_n tol=0"
               echo "     $LD" | sed 's/^/     /'
               rc=1 ;;
         esac
         case "$LS" in
            *"a carrier is its density"*)
               echo "PASS carrier_unknown_switch_restores_the_density measured=density reference=density tol=0" ;;
            *) echo "FAIL carrier_unknown_switch_restores_the_density measured=not_density reference=density tol=0"
               echo "     $LS" | sed 's/^/     /'
               rc=1 ;;
         esac
      fi
      # AND THE FLOOR IS THE SAME NUMBER IN BOTH, because it is formed from
      # the state and not from the unknown space: the density unknown and
      # the logarithm are then measured against one yardstick.
      FD="$(grep -m1 '(JFNK) carrier floor' "$CD/run.log" |
            sed -n 's/.*smallest *\([0-9.E+-]*\), largest *\([0-9.E+-]*\).*/\1 \2/p' || true)"
      FS="$(grep -m1 '(JFNK) carrier floor' "$CS/run.log" |
            sed -n 's/.*smallest *\([0-9.E+-]*\), largest *\([0-9.E+-]*\).*/\1 \2/p' || true)"
      if [ -n "$FD" ] && [ "$FD" = "$FS" ]; then
         echo "PASS carrier_floor_is_the_same_in_both_spaces measured=$FD reference=$FS tol=0"
      else
         echo "FAIL carrier_floor_is_the_same_in_both_spaces measured=$FD reference=$FS tol=0"
         rc=1
      fi
   fi
fi

# ------------------------------------------------------------------ #
# AND NO ADOPTED CARRIER SITS AT THE FLOOR OF ITS OWN UNKNOWN SPACE (B5k).
#
# The floor under ln n is 1e-20 of the cell's own element budget -- the
# fraction of its element below which carrier_residual's absolute floor
# already calls a row converged -- so it is a bound the space HAS and not one
# the solve is meant to work against.  The solve reports the smallest ratio of
# an adopted carrier density to that floor and the count of adopted carriers at
# or below it; above one and zero is the statement that the space never had to
# bound anything.  The row is exact: the floor is a bound, not a tolerance.
#
# RED and GREEN, MEASURED on the coupled `mol_carrier` reload, one binary and
# one environment variable apart: with the density unknown the Newton drives
# the carrier of cell 215 to exactly zero, so the ratio is 0 and 500 adopted
# carriers sit at or below the floor; with ln n the ratio is above 1e6 and the
# count is zero.
#
# EXHALE_CARRIER_FLOOR_LOGS is a colon-separated list of run logs.
# ------------------------------------------------------------------ #
FLOOR_LOGS="${EXHALE_CARRIER_FLOOR_LOGS:-}"
if [ -n "$FLOOR_LOGS" ]; then
   IFS=':' read -r -a FLARR <<< "$FLOOR_LOGS"
   for FLOG in "${FLARR[@]}"; do
      FTAG="$(basename "$(dirname "$FLOG")")"
      if [ ! -f "$FLOG" ]; then
         echo "FAIL carrier_is_above_its_own_floor_$FTAG measured=missing reference=$FLOG tol=0"
         rc=1; continue
      fi
      RAT="$(grep 'smallest adopted carrier density over the floor' "$FLOG" |
             tail -n 1 |
             sed -n 's/.*of its own cell: *\([0-9.E+-]*\).*/\1/p' || true)"
      NAT="$(grep 'smallest adopted carrier density over the floor' "$FLOG" |
             tail -n 1 |
             sed -n 's/.*at or below that floor: *\([0-9]*\).*/\1/p' || true)"
      if [ -z "$RAT" ]; then
         # A run stopped before its end prints the same two numbers on every
         # outer iteration; the last of those is the same running minimum.
         RAT="$(grep 'smallest carrier over its own floor' "$FLOG" | tail -n 1 |
                sed -n 's/.*over its own floor *\([0-9.E+-]*\).*/\1/p' || true)"
         NAT="$(grep 'smallest carrier over its own floor' "$FLOG" | tail -n 1 |
                sed -n 's/.*carriers at that floor *\([0-9]*\).*/\1/p' || true)"
      fi
      if [ -z "$RAT" ] || [ -z "$NAT" ]; then
         echo "FAIL carrier_is_above_its_own_floor_$FTAG measured=no_line reference=a_coupled_carrier_solve tol=0"
         rc=1; continue
      fi
      awk -v v="$RAT" -v t="$FTAG" 'BEGIN{
         if (v+0 > 1.0) printf "PASS carrier_is_above_its_own_floor_%s measured=%s reference=>1 tol=0\n", t, v;
         else { printf "FAIL carrier_is_above_its_own_floor_%s measured=%s reference=>1 tol=0\n", t, v; exit 1 }
      }' || rc=1
      if [ "$NAT" -eq 0 ]; then
         echo "PASS carrier_never_reached_its_floor_$FTAG measured=$NAT reference=0 tol=0"
      else
         echo "FAIL carrier_never_reached_its_floor_$FTAG measured=$NAT reference=0 tol=0"
         rc=1
      fi
   done
fi

# ------------------------------------------------------------------ #
# AND NO ADOPTED CARRIER EXCEEDS THE ELEMENT BUDGET OF ITS OWN CELL (B5l).
#
# A carrier holds nuclei of one element, so it cannot hold more of that
# element than the cell has free (carrier_element_headroom).  Before this
# item the budget was not a face of the unknown box: it was enforced only by
# refusing an already evaluated state, so the step walked into it and was
# halved.  With the budget a face the solve reports two numbers over every
# state it adopted -- the largest carrier density over its own budget, and
# how many adopted carriers stood above it -- and the statement is that the
# first is at or below one and the second is zero.  Both rows are exact: a
# budget is a bound, not a tolerance.
#
# The refusal census of the same log is the other half: with the budget a
# face, no residual sample can be refused for leaving it, so that entry of
# the census must be zero.
#
# EXHALE_CARRIER_BUDGET_LOGS is a colon-separated list of run logs.
# ------------------------------------------------------------------ #
BUDGET_LOGS="${EXHALE_CARRIER_BUDGET_LOGS:-}"
if [ -n "$BUDGET_LOGS" ]; then
   IFS=':' read -r -a BUARR <<< "$BUDGET_LOGS"
   for BULOG in "${BUARR[@]}"; do
      BUTAG="$(basename "$(dirname "$BULOG")")"
      if [ ! -f "$BULOG" ]; then
         echo "FAIL carrier_is_inside_its_element_budget_$BUTAG measured=missing reference=$BULOG tol=0"
         rc=1; continue
      fi
      RAT="$(grep 'largest adopted carrier density over the element budget' "$BULOG" |
             tail -n 1 |
             sed -n 's/.*of its own cell: *\([0-9.E+-]*\).*/\1/p' || true)"
      NAT="$(grep 'largest adopted carrier density over the element budget' "$BULOG" |
             tail -n 1 |
             sed -n 's/.*above that budget: *\([0-9]*\).*/\1/p' || true)"
      if [ -z "$RAT" ]; then
         # A run stopped before its end prints the same two numbers on every
         # outer iteration; the last of those is the same running maximum.
         RAT="$(grep 'largest carrier over its own element budget' "$BULOG" | tail -n 1 |
                sed -n 's/.*over its own element budget *\([0-9.E+-]*\).*/\1/p' || true)"
         NAT="$(grep 'largest carrier over its own element budget' "$BULOG" | tail -n 1 |
                sed -n 's/.*above that budget *\([0-9]*\).*/\1/p' || true)"
      fi
      if [ -z "$RAT" ] || [ -z "$NAT" ]; then
         echo "FAIL carrier_is_inside_its_element_budget_$BUTAG measured=no_line reference=a_coupled_carrier_solve tol=0"
         rc=1; continue
      fi
      awk -v v="$RAT" -v t="$BUTAG" 'BEGIN{
         if (v+0 <= 1.0) printf "PASS carrier_is_inside_its_element_budget_%s measured=%s reference=<=1 tol=0\n", t, v;
         else { printf "FAIL carrier_is_inside_its_element_budget_%s measured=%s reference=<=1 tol=0\n", t, v; exit 1 }
      }' || rc=1
      if [ "$NAT" -eq 0 ]; then
         echo "PASS carrier_never_left_its_element_budget_$BUTAG measured=$NAT reference=0 tol=0"
      else
         echo "FAIL carrier_never_left_its_element_budget_$BUTAG measured=$NAT reference=0 tol=0"
         rc=1
      fi
      NREF="$(grep 'residual samples: admitted' "$BULOG" | tail -n 1 |
              sed -n 's/.*cells outside the element budget \([0-9]*\).*/\1/p' || true)"
      if [ -z "$NREF" ]; then
         echo "FAIL carrier_budget_refuses_no_sample_$BUTAG measured=no_line reference=a_completed_solve tol=0"
         rc=1
      elif [ "$NREF" -eq 0 ]; then
         echo "PASS carrier_budget_refuses_no_sample_$BUTAG measured=$NREF reference=0 tol=0"
      else
         echo "FAIL carrier_budget_refuses_no_sample_$BUTAG measured=$NREF reference=0 tol=0"
         rc=1
      fi
   done
fi

# ------------------------------------------------------------------ #
# THE RELOAD ROWS (B5b item 4).
#
# A state written by one run and read back by another with "Load IC? True"
# is certified again at the t = 0 handoff, and the species rows are among
# the equations that certification measures.  What the rows below state is
# that a reloaded state which does not satisfy a species balance is REFUSED,
# that the refusal names the balance, and that the run exits 2 -- because a
# du stop is a stationarity claim about the state it writes and the exit
# status is what tells a caller whether the claim was accepted.
#
# The case is the hot Uranus molecular gate with the carrier transport on,
# which is the configuration that carries a transported balance at all.  Run
# a: a step-capped snapshot, which declares nothing and exits 0.  Run b:
# that snapshot reloaded with a du threshold above its own du, so the stop
# fires at the first evaluation and the claim is made on a relaxation
# snapshot the carrier balance is nowhere near.
#
# EXHALE_EXE gives the binary, EXHALE_SPECIES_EXE being its alias here;
# without either of them these rows are skipped rather than building one
# here.  The selection is made at the head of this file.
# ------------------------------------------------------------------ #
if [ "$EXHALE_RUN_EXE_REQUESTED" != "1" ]; then
   echo "  reload rows skipped: set EXHALE_EXE (or its alias EXHALE_SPECIES_EXE)"
   echo "  to the binary to test"
   exit $rc
fi
EXE_RUN="$EXHALE_RUN_EXE"
if [ ! -x "$EXE_RUN" ]; then
   echo "FAIL species_rows_reload_binary measured=no_binary reference=$EXE_RUN tol=0"
   exit 1
fi

REG="$ROOT/backup/regression"
WORK="$OUT/reload"
rm -rf "$WORK"; mkdir -p "$WORK/a/output" "$WORK/b/output"
for f in input.inp base.inp metals.inp opacity.inp; do
   [ -f "$REG/mol_carrier/$f" ] && cp "$REG/mol_carrier/$f" "$WORK/a/"
   [ -f "$REG/mol_carrier/$f" ] && cp "$REG/mol_carrier/$f" "$WORK/b/"
done
# The staged secondary-ionization coupling is switched off in both runs:
# staged, the first du stop flips the coupling on and holds every stop for
# N_stall steps instead, and this pair is about the certification of a
# reloaded state and not about that staging.
for d in a b; do
   echo 'Secondary_ionization: False' >> "$WORK/$d/input.inp"
done
( cd "$WORK/a" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40 "$EXE_RUN" \
     > run.log 2>&1 )
ra=$?
if [ "$ra" -eq 0 ]; then
   echo "PASS species_rows_reload_snapshot measured=0 reference=0 tol=0"
else
   echo "FAIL species_rows_reload_snapshot measured=$ra reference=0 tol=0"
   rc=1
fi

sed -i -e 's/^Load IC?.*/Load IC? True/' \
       -e 's/^du_th .*/du_th [PLM,WENO3]: 1.0e3 1.0e3/' \
       -e 's/^Solver:.*/Solver: None/' "$WORK/b/input.inp"
cp "$WORK/a/output/Hydro_ioniz.txt" "$WORK/b/output/Hydro_ioniz_IC.txt"
cp "$WORK/a/output/Ion_species.txt" "$WORK/b/output/Ion_species_IC.txt"
( cd "$WORK/b" && OMP_NUM_THREADS=1 "$EXE_RUN" > run.log 2>&1 )
rb=$?

n=$(grep -c 'stopped: mass flux' "$WORK/b/run.log" || true)
if [ "$n" -ge 1 ]; then
   echo "PASS species_rows_reload_declares_a_stationary_state measured=$n reference=>=1 tol=0"
else
   echo "FAIL species_rows_reload_declares_a_stationary_state measured=$n reference=>=1 tol=0"
   rc=1
fi
# The refusal names the row as "gated row measure" since the carrier rows
# gained the absent-carrier floor (2026-09-14, L7b); both spellings count.
n=$(grep -cE 'carrier balance H2: (gated )?row measure' "$WORK/b/run.log" || true)
if [ "$n" -ge 1 ]; then
   echo "PASS species_rows_reload_names_the_refused_row measured=$n reference=>=1 tol=0"
else
   echo "FAIL species_rows_reload_names_the_refused_row measured=$n reference=>=1 tol=0"
   rc=1
fi
if [ "$rb" -eq 2 ]; then
   echo "PASS species_rows_reload_refusal_exits_2 measured=$rb reference=2 tol=0"
else
   echo "FAIL species_rows_reload_refusal_exits_2 measured=$rb reference=2 tol=0"
   rc=1
fi

exit $rc

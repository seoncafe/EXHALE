#!/bin/bash
# THE GHOST COMPOSITION IS A FIXED POINT OF THE SWEEP THAT RETURNS IT, AND
# THE BOUNDARY STATES WHICH COMPOSITION THAT SWEEP STARTS FROM (item P1
# step 2 of docs/PLAN_20260919_rev1.md; the measurements are
# docs/lhs1140b_p1_step2_20260920.md).
#
# WHAT IS UNDER TEST, row by row.
#
#   1  THE RETURNED GHOST DOES NOT DEPEND ON THE SEED. The same state is
#      evaluated three times: with the seed the boundary model states (the
#      reservoir row, which load_IC installs in place of the ghost rows the
#      restart pair carries), with the pair's own ghost rows, and with a
#      seed whose ionized hydrogen is thirty times the ghost's own, taken
#      out of its neutral hydrogen at fixed hydrogen nuclei. The last is the
#      decisive one: the ionized seed is the direction the returned ghost
#      responds to, with a gain of 4.9e-3 measured over one application
#      (READ, docs/lhs1140b_p1_step1b_20260919.md section 6.1), and it is
#      the size of the difference between the reservoir row and a solved
#      ghost on the state that opened item P1. The two lower ghost rows of
#      the written state must agree within the accuracy the contract states.
#
#   2  AND NEITHER DO THE PHYSICAL CELLS. The same comparison over the 500
#      physical rows, which the seed has no business moving at all.
#
#      MEASURED: on THIS fixture the entry text of item P1 step 2 already
#      passes both, its returned ghost standing 1.3e-7 from the one another
#      seed gives, so the rows guard the property rather than reproduce the
#      failure. The state that fails them is the LHS 1140 b hot-Uranus
#      molecular generation of the catalog, whose two seeds give ghosts
#      5.5e-2 apart and verdicts 8.335e-08 against 2.923e-09; a catalog
#      generation is not a fixture and the measurement is in
#      docs/lhs1140b_p1_step2_20260920.md.
#
#   3  THE BOUNDARY REPORTS WHAT THE CONTRACT REACHED. The boundary model
#      block carries the move of the ghost's species densities over the last
#      application of the map, in the composition's own units, and it is
#      within the accuracy.
#
#   4  THE THERMODYNAMIC LEG CLOSES WITH IT. The same block carries the move
#      of the ghost's heavy-particle plus electron count, the quantity that
#      turns the ghost's temperature into its pressure and its internal
#      energy, and it is within its own accuracy.
#
#   5  THE SEED IS PART OF THE BOUNDARY MODEL. The model identity the run
#      states names the ghost fixed point and the seed, and the boundary
#      report states the seed in words.
#
# THE FIXTURE
#   backup/regression/carrier_model_a_newton/IC/, the hot-Uranus molecular
#   carrier state, evaluated with "Restart intent: stationary evaluate" on a
#   copy. The regression directory is never written to.
#
# THE TOLERANCES
#   1e-6 relative on rows 1 and 2, which compare two written state files
#   column by column.
#   1e-11 on the reported composition move and 1e-6 on the reported count
#   move, the two accuracies the boundary states
#   (ionization_equilibrium, ghost_composition_fixed_point_move and
#   ghost_count_fixed_point_move). The first is in the composition's own
#   units, a quarter of the rounding floor of the base continuity row that
#   reads the ghost (MEASURED: 107.4 floors per 4.1e-9 of composition,
#   docs/lhs1140b_p1_step1_20260919.md sections 5.1 and 5.2) and above the
#   plateau the composition solve pins the ghost to (MEASURED at 6.5e-12,
#   docs/lhs1140b_p1_step2b_20260920.md). The second is against the count's
#   own value, which is an O(1) quantity.
#
# Usage: ghost_composition_fixed_point.sh   (EXHALE_EXE selects the binary,
#        EXHALE_TEST_OUT the work directory)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" ghost_composition_fixed_point
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/ghost_composition_fixed_point"
CASE="$ROOT/backup/regression/carrier_model_a_newton"
TOL=1.0e-6
TOL_FP=1.0e-11
TOL_NP=1.0e-6

n_fail=0
check_le() {   # check_le <name> <measured> <tolerance>
   local v
   v=$(awk -v m="${2:-}" -v t="$3" 'BEGIN {
          if (m == "") { print "missing" }
          else if (m + 0 <= t + 0) { print "within" } else { print "above" } }')
   if [ "$v" = "within" ]; then
      echo "PASS $1 measured=${2:-missing} reference=at_most tol=$3"
   else
      echo "FAIL $1 measured=${2:-missing} reference=at_most tol=$3"
      n_fail=$((n_fail + 1))
   fi
}
check_eq() {   # check_eq <name> <measured> <reference>
   if [ "$2" = "$3" ]; then
      echo "PASS $1 measured=$2 reference=$3 tol=0"
   else
      echo "FAIL $1 measured=$2 reference=$3 tol=0"
      n_fail=$((n_fail + 1))
   fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL ghost_composition_fixed_point measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   echo "FAIL ghost_composition_fixed_point measured=no_fixture reference=$CASE/IC tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"

stage() {   # stage <dir>
   local d="$WORK/$1"
   mkdir -p "$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$d/"
   sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' \
      "$d/input.inp"
   printf 'Restart intent: stationary evaluate\n' >> "$d/input.inp"
   cp "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
      "$d/output/"
}

run_it() {   # run_it <dir> [env assignments...]
   local d="$1"; shift
   ( cd "$WORK/$d" && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        "$@" "$EXE" > run.log 2>&1 )
   return 0
}

# The largest relative difference between two state files over a row range.
# Values whose larger magnitude is below 1e-30 are skipped: a run without
# metals writes its metal columns at a floor of order 1e-301 cm^-3 whose last
# digits are not a composition.
worst_over() {   # worst_over <file a> <file b> <first data row> <last data row>
   awk -v lo="$3" -v hi="$4" '
      FNR == NR { if ($0 !~ /^#/ && NF) { n++; a[n] = $0 } ; next }
      $0 !~ /^#/ && NF {
         m++
         if (m < lo || m > hi) next
         split(a[m], x); split($0, y)
         for (k = 1; k <= NF; k++) {
            u = x[k] + 0; v = y[k] + 0
            s = (u < 0 ? -u : u); t = (v < 0 ? -v : v)
            if (t > s) s = t
            if (s < 1e-30) continue
            d = (u > v ? u - v : v - u)/s
            if (d > w) w = d
         }
      }
      END { printf "%.6e\n", w + 0 }' "$1" "$2"
}

stage seed_of_record
run_it seed_of_record EXHALE_GHOST_SEED_WRITE="$WORK/seed_of_record.txt"
stage seed_state_file
run_it seed_state_file EXHALE_GHOST_COMPOSITION_SEED=the_state_file_ghost_rows

# The ionized seed: thirty times the ghost's own H II, taken out of its H I
# at fixed hydrogen nuclei. f(1) is H I and f(2) is H II of the f_sp layout,
# which the seed file writes after the cell index and the species count.
if [ -s "$WORK/seed_of_record.txt" ]; then
   awk 'BEGIN { OFMT = "%.17g"; CONVFMT = "%.17g" }
        /^#/ { print; next }
        $1 == "ghost" { d = 29*$5; $5 = $5 + d; $4 = $4 - d; print; next }
        { print }' "$WORK/seed_of_record.txt" > "$WORK/seed_ionized.txt"
   stage seed_ionized
   run_it seed_ionized \
      EXHALE_GHOST_COMPOSITION_SEED="$WORK/seed_ionized.txt"
else
   echo "FAIL ghost_seed_rows_written measured=missing reference=written tol=0"
   n_fail=$((n_fail + 1))
fi

A="$WORK/seed_of_record/output"
B="$WORK/seed_state_file/output"
if [ ! -s "$A/Ion_species.txt" ] || [ ! -s "$B/Ion_species.txt" ]; then
   echo "FAIL ghost_composition_fixed_point measured=no_state reference=written tol=0"
   exit 1
fi

# The two lower ghost cells are data rows 1 and 2; the physical cells follow.
nrow=$(grep -vc '^#' "$A/Ion_species.txt")
check_le ghost_rows_are_independent_of_the_state_file_seed_Ion_species \
   "$(worst_over "$A/Ion_species.txt" "$B/Ion_species.txt" 1 2)" "$TOL"
check_le ghost_rows_are_independent_of_the_state_file_seed_Hydro_ioniz \
   "$(worst_over "$A/Hydro_ioniz.txt" "$B/Hydro_ioniz.txt" 1 2)" "$TOL"
check_le physical_cells_are_independent_of_the_state_file_seed \
   "$(worst_over "$A/Ion_species.txt" "$B/Ion_species.txt" 3 "$nrow")" "$TOL"
C="$WORK/seed_ionized/output"
if [ -s "$C/Ion_species.txt" ]; then
   check_le ghost_rows_are_independent_of_the_ionized_seed_Ion_species \
      "$(worst_over "$A/Ion_species.txt" "$C/Ion_species.txt" 1 2)" "$TOL"
   check_le ghost_rows_are_independent_of_the_ionized_seed_Hydro_ioniz \
      "$(worst_over "$A/Hydro_ioniz.txt" "$C/Hydro_ioniz.txt" 1 2)" "$TOL"
   check_le physical_cells_are_independent_of_the_ionized_seed \
      "$(worst_over "$A/Ion_species.txt" "$C/Ion_species.txt" 3 "$nrow")" "$TOL"
else
   echo "FAIL ghost_rows_are_independent_of_the_ionized_seed_Ion_species measured=no_state reference=at_most tol=$TOL"
   n_fail=$((n_fail + 1))
fi

# The boundary's own report of what the contract reached.
fp=$(awk '/ghost composition is a fixed point of its own solve to/ \
          { print $12; exit }' "$WORK/seed_of_record/run.log")
np=$(awk '/its particle plus electron count to/ { v = $7; sub(",", "", v);
          print v; exit }' "$WORK/seed_of_record/run.log")
check_le ghost_composition_fixed_point_is_reached "$fp" "$TOL_FP"
check_le ghost_particle_and_electron_count_closes "$np" "$TOL_NP"

# The seed and the fixed point travel in the model identity and the report.
if grep -q 'model  characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3' \
     "$WORK/seed_of_record/run.log"; then id=stated; else id=absent; fi
check_eq boundary_model_identity_carries_the_ghost_seed "$id" stated
if grep -q 'ghost composition seed: the reservoir row at the first sweep' \
     "$WORK/seed_of_record/run.log"; then sd=stated; else sd=absent; fi
check_eq boundary_report_states_the_ghost_seed "$sd" stated

echo ""
echo "  applications of the composition map, seed of record:" \
     "$(awk '/its particle plus electron count to/ { print $(NF-1); exit }' \
        "$WORK/seed_of_record/run.log")"
echo "  applications, seeded with the state file ghost rows:" \
     "$(awk '/its particle plus electron count to/ { print $(NF-1); exit }' \
        "$WORK/seed_state_file/run.log")"

echo ""
if [ $n_fail -gt 0 ]; then
   echo "ghost_composition_fixed_point: $n_fail assertion(s) failed"
   exit 1
fi
exit 0

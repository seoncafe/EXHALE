#!/bin/bash
# THE LOWER BOUNDARY IS A FUNCTION OF THE DECLARED INPUTS, AND OF NOTHING
# ELSE (item D5b-2 of docs/PLAN_20260918_rev2.md; the measurements behind it
# are docs/lhs1140b_stationary_D5a_20260918.md and
# docs/lhs1140b_stationary_D5b2_20260918.md).
#
# WHAT IS UNDER TEST, row by row.
#
#   1  THE GHOST ROWS OF A RESTART FILE ARE NOT AN INPUT. The boundary-state
#      operation takes the physical conserved state and composition of the
#      column, the prescribed reservoir, the radiation context and the model
#      options; the two rows a file carries below the base are the
#      operation's own output. Measured with them in (D5a table 3): holding
#      every physical cell and every declared input fixed and exchanging
#      only those two species rows for another admissible ghost moved the
#      cell-1 continuity row of one hot-Uranus state over eight decades and
#      reversed the sign of the base face mass flux. The two runs below
#      differ in exactly those rows and must agree in every data line.
#
#   2  THE GHOST'S MOLECULAR PARTITION IS SOLVED, NOT PRESCRIBED FROM THE
#      COMPOSITION THE STATE WAS ENTERED WITH. x_H2 = x2 (1 - x_ion) is a
#      condition on the ghost coupled to the ghost's own ionization balance;
#      the sweep closes the pair and reports the residual it reached.
#
#   3  THE PRESCRIBED RESERVOIR AND THE SOLVED GHOST COUNT ARE TWO
#      QUANTITIES. The reservoir the boundary stands on is the one the input
#      resolved, before and after any sweep; the ghost's own particle count
#      is a measurement of the state and is never written back. On this
#      fixture the two differ by per cent, so a run that confused them would
#      not pass by coincidence.
#
#   4  EVERY READ OF THE CACHED FACE STATE BELONGS TO THE COMPOSITION
#      INSTALLED AT THE READ. The cache carries the input set it was derived
#      from, and a read that finds the composition moved derives it again
#      before the Riemann solve sees it. Zero recomputes on this route is
#      the call order being right; the count is the statement that nothing
#      stale was used either way.
#
#   5  THE BOUNDARY MODEL AND THE PRESCRIBED RESERVOIR TRAVEL WITH THE
#      STATE. A state file carries the identity of the boundary it was
#      produced under and the reservoir that boundary was prescribed with,
#      which is the only boundary input the physical column cannot
#      reconstruct; a restart of that file recognizes them.
#
#   6  THE TWO EVALUATION SITES OF A SOLVE DERIVE THE BOUNDARY AFTER THE
#      SWEEP. The joint test inside steady_wind_with_element_diffusion and
#      the final certification of a stationary restart assemble their
#      residual on the boundary of the composition they are taken at, as the
#      three evaluation routes of D5b-1 already do.
#
# THE FIXTURE
#   backup/regression/carrier_model_a_newton, the hot-Uranus molecular
#   carrier state, on copies; the regression directory is never written to.
#   It is a relaxation snapshot and not a stationary solution, so nothing
#   here asserts anything about its rows or its certification. A molecular
#   base is required: rows 2 and 3 have no subject without one.
#
# Usage: boundary_from_the_declared_inputs.sh   (EXHALE_EXE selects the
#        binary, EXHALE_TEST_OUT the work directory)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" boundary_from_the_declared_inputs
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/boundary_from_the_declared_inputs"
CASE="$ROOT/backup/regression/carrier_model_a_newton"

n_fail=0
check_eq() {   # check_eq <name> <measured> <reference>
   if [ "$2" = "$3" ]; then
      echo "PASS $1 measured=$2 reference=$3 tol=0"
   else
      echo "FAIL $1 measured=$2 reference=$3 tol=0"
      n_fail=$((n_fail + 1))
   fi
}
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

if [ ! -x "$EXE" ]; then
   echo "FAIL boundary_from_the_declared_inputs measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   echo "FAIL boundary_from_the_declared_inputs measured=no_fixture reference=$CASE/IC tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"

stage() {   # stage <dir> <intent>
   local d="$WORK/$1"
   mkdir -p "$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$d/"
   sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' \
      "$d/input.inp"
   printf 'Restart intent: %s\n' "$2" >> "$d/input.inp"
   cp "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
      "$d/output/"
}

run_it() {   # run_it <dir> [env assignments...]
   local d="$1"; shift
   ( cd "$WORK/$d" && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        "$@" "$EXE" > run.log 2>&1 )
   return 0
}

data_md5() {   # data_md5 <dir> <file>
   grep -v '^ *#' "$WORK/$1/output/$2" 2>/dev/null | md5sum | cut -d' ' -f1
}

# ---- the two runs that differ only in the file's lower ghost rows --------
stage own   'stationary evaluate'
stage other 'stationary evaluate'
# An admissible ghost that is not the file's: the composition of the first
# physical cell, rescaled to the ghost's own heavy-particle count. It is what
# an "extrapolate the interior" rule would supply, and D5a measured it
# reversing the sign of the base face mass flux when the file's ghost was an
# input.
python3 - "$WORK/other/output/Ion_species_IC.txt" <<'PY'
import sys
path = sys.argv[1]
head, rows = [], []
for line in open(path):
    (head if line.lstrip().startswith('#') else rows).append(line.rstrip('\n'))
rows = [r for r in rows if r.strip()]
c1 = [float(x) for x in rows[2].split()]
for g in (0, 1):
    v = [float(x) for x in rows[g].split()]
    s = sum(v[1:])/sum(c1[1:])
    v = [v[0]] + [x*s for x in c1[1:]]
    rows[g] = ''.join('  %22.16E' % x for x in v)
open(path, 'w').write(''.join(l + '\n' for l in head) + '\n'.join(rows) + '\n')
PY
run_it own   EXHALE_BOUNDARY_TRACE=1
run_it other EXHALE_BOUNDARY_TRACE=1

if [ ! -s "$WORK/own/output/Hydro_ioniz.txt" ]; then
   echo "FAIL evaluation_wrote_its_state measured=missing reference=written tol=0"
   tail -n 20 "$WORK/own/run.log"
   exit 1
fi

for f in Hydro_ioniz.txt Ion_species.txt; do
   a="$(data_md5 own "$f")"; b="$(data_md5 other "$f")"
   check_eq "ghost_rows_of_the_file_are_not_an_input_$f" \
      "$( [ "$a" = "$b" ] && echo same || echo different )" same
done
# The signed continuity rows of the base cells and the two face mass fluxes
# each of them differences, as the trace writes them: the first residual of
# the two runs at sixteen digits, where the printed norms carry five.
rows_own="$(awk '/^row /' "$WORK/own/output/boundary_trace.txt" 2>/dev/null \
            | md5sum | cut -d' ' -f1)"
rows_oth="$(awk '/^row /' "$WORK/other/output/boundary_trace.txt" 2>/dev/null \
            | md5sum | cut -d' ' -f1)"
have_rows="$(awk '/^row /' "$WORK/own/output/boundary_trace.txt" 2>/dev/null \
             | wc -l)"
check_eq ghost_rows_of_the_file_are_not_an_input_first_residual \
   "$( [ "$have_rows" -gt 0 ] && [ "$rows_own" = "$rows_oth" ] \
       && echo same || echo different )" same

# ---- the ghost's molecular partition, closed against its own ionization --
clos="$(awk '/ghost H2 partition closed to/ { print $6 }' \
        "$WORK/own/run.log" | tail -n 1)"
check_le ghost_h2_partition_closure_residual "$clos" 1.0e-10
echo "  passes the closure took: $(awk '/ghost H2 partition closed to/ \
   { print $(NF-1) }' "$WORK/own/run.log" | tail -n 1)"

# ---- the prescribed reservoir against the solved ghost count -------------
# The reservoir the written state carries, after every sweep of the run, must
# be the one the input resolved, and it must not have become the ghost's own
# particle count.
res_start="$(awk '/^ *p \[p0\]/ { print $3; exit }' "$WORK/own/run.log")"
res_file="$(awk '$2 == "boundary_reservoir" { print $9 }' \
            "$WORK/own/output/Hydro_ioniz.txt" | head -n 1)"
# Compared as numbers: the report and the file render the same double with
# different edit descriptors.
check_eq prescribed_reservoir_is_serialized_unchanged \
   "$(awk -v a="${res_file:-}" -v b="${res_start:-}" 'BEGIN {
        if (a == "" || b == "") print "missing";
        else if (a + 0 == b + 0) print "same"; else print "different" }')" same
ghost_n="$(awk '/solved ghost counts/ { print $6 }' "$WORK/own/run.log" \
           | tail -n 1)"
# Both have to be there: two quantities neither of which was reported is not
# the statement this row makes.
check_eq prescribed_reservoir_and_solved_ghost_count_are_two_quantities \
   "$(awk -v a="${ghost_n:-}" -v b="${res_start:-}" 'BEGIN {
        if (a == "" || b == "") print "missing";
        else if (a + 0 == b + 0) print "one"; else print "two" }')" two
echo "  prescribed reservoir pressure ${res_start:-missing}, solved ghost"
echo "  heavy-particle count ${ghost_n:-missing}"

# ---- every read of the cached face state is of the installed composition -
rec="$(awk '/face states recomputed at a read/ { print $NF }' \
       "$WORK/own/run.log" | tail -n 1)"
check_eq face_state_cache_reads_are_all_current "${rec:-missing}" 0

# ---- the boundary model travels with the state --------------------------
mod_h="$(awk '$2 == "boundary_model" { print $3 }' \
         "$WORK/own/output/Hydro_ioniz.txt" | head -n 1)"
mod_i="$(awk '$2 == "boundary_model" { print $3 }' \
         "$WORK/own/output/Ion_species.txt" | head -n 1)"
check_eq boundary_model_written_on_both_halves \
   "$( [ -n "$mod_h" ] && [ "$mod_h" = "$mod_i" ] && echo both || echo no )" both
stage reload 'stationary evaluate'
cp "$WORK/own/output/Hydro_ioniz.txt" "$WORK/reload/output/Hydro_ioniz_IC.txt"
cp "$WORK/own/output/Ion_species.txt" "$WORK/reload/output/Ion_species_IC.txt"
run_it reload
check_eq boundary_model_recognized_on_reload \
   "$(grep -c "boundary model of the restart: ${mod_h:-none} (this run's)" \
      "$WORK/reload/run.log")" 1

# ---- the two evaluation sites of a solve --------------------------------
# One outer pass and a two-iteration hydrodynamic solve: enough that the
# joint test and the final certification are reached, and nothing is asserted
# about what the solve achieved.
stage solve 'stationary'
run_it solve EXHALE_BOUNDARY_TRACE=1 EXHALE_OUTER_PASSES=1 \
   EXHALE_JFNK_MAXIT=2
gap_at() {   # gap_at <dir> <point label>
   awk -v pat="$2" '$1 == "face_state_gap" && $2 == pat { print $3 }' \
      "$WORK/$1/output/boundary_trace.txt" 2>/dev/null | head -n 1
}
for point in steady_wind_joint_test_at_assemble_residual \
             stationary_restart_final_certification; do
   g="$(gap_at solve "$point")"
   check_eq "boundary_at_$point" \
      "$(awk -v g="${g:-}" 'BEGIN { if (g == "") print "missing";
          else if (g + 0 == 0) print "zero"; else print "nonzero" }')" zero
done

echo ""
if [ $n_fail -gt 0 ]; then
   echo "boundary_from_the_declared_inputs: $n_fail assertion(s) failed"
   exit 1
fi
exit 0

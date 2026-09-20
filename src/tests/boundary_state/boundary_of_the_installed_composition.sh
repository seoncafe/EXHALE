#!/bin/bash
# THE BOUNDARY A RESIDUAL STANDS ON BELONGS TO THE COMPOSITION THE RESIDUAL
# IS ASSEMBLED WITH (item D5b-1 of docs/PLAN_20260918_rev2.md; the
# measurements are docs/lhs1140b_stationary_D5a_20260918.md and
# docs/lhs1140b_stationary_D5b1_20260918.md).
#
# WHAT IS UNDER TEST
#   The evaluation routes apply the lower boundary BEFORE the composition
#   sweep, because U_to_W needs a ghost density to divide by and ioniz_eq a
#   ghost pressure.  The ghost conserved state and the base face state that
#   Apply_BC then caches in BC_Apply therefore carry the composition the
#   state was loaded with, while the residual's sources, pressure map and
#   sound speeds carry the composition the sweep returned; Rec_BC puts the
#   cached face state straight into the left slot of the base face.  Under
#   the caloric equation of state those are two different gases.
#
#   The route derives the boundary once more, from the composition the
#   residual is assembled with, so that the two coincide.  The measure below
#   is the largest relative difference between the CACHED boundary (the base
#   face state, the ghost cell averages and the state at the next face down)
#   and the boundary base_boundary_states gives for the state installed at
#   that instant.  It is zero exactly where the invariant holds.
#
# WHY THERE ARE TWO RUNS
#   A row that measures zero because nothing was measured passes silently.
#   The second run sets EXHALE_TRACE_EXPERIMENT to leave the boundary of the
#   entry composition standing, which is the order this item corrects, and
#   the same measure must then be strictly positive on a molecular base.
#
# THE FIXTURE
#   backup/regression/carrier_model_a_newton/IC/, the hot-Uranus molecular
#   carrier state, evaluated with "Restart intent: stationary evaluate" on a
#   copy.  The regression directory is never written to.  It is a relaxation
#   snapshot and not a stationary solution, so its rows are far from zero and
#   its certification refuses; neither is asserted here.  The measure under
#   test is a property of the evaluation and not of the state's stationarity.
#   A molecular base is required: on an atomic mixture the ghost continuation
#   carries no composition (adiabatic_index_from_state returns gamma_ad), so
#   the measure is zero in both orders and the second run would assert
#   nothing.
#
# THE TOLERANCES
#   None on the first row: one composition and one boundary is an identity,
#   and the derivation is the one Apply_BC performs.
#
# Usage: boundary_of_the_installed_composition.sh   (EXHALE_EXE selects the
#        binary)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" boundary_of_the_installed_composition
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/boundary_of_the_installed_composition"
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

if [ ! -x "$EXE" ]; then
   echo "FAIL boundary_of_the_installed_composition measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   echo "FAIL boundary_of_the_installed_composition measured=no_fixture reference=$CASE/IC tol=0"
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
        EXHALE_BOUNDARY_TRACE=1 "$@" "$EXE" > run.log 2>&1 )
   return 0
}

# The gap the trace wrote at the named point, as printed.
gap_at() {   # gap_at <dir> <point label>
   awk -v pat="$2" '$1 == "face_state_gap" && $2 == pat { print $3 }' \
      "$WORK/$1/output/boundary_trace.txt" | head -n 1
}

# "zero" / "nonzero" of a printed value, without a locale-dependent printf.
sign_of() {   # sign_of <value>
   awk -v g="${1:-}" 'BEGIN {
      if (g == "") { print "missing" }
      else if (g + 0 == 0) { print "zero" }
      else { print "nonzero" } }'
}

# ---- the boundary of the composition the residual is assembled with ------
stage installed
run_it installed
if [ ! -s "$WORK/installed/output/boundary_trace.txt" ]; then
   echo "FAIL boundary_trace_written measured=missing reference=written tol=0"
   tail -n 20 "$WORK/installed/run.log"
   exit 1
fi
g_installed="$(gap_at installed evaluate_route_at_assemble_residual)"
check_eq face_state_gap_at_assemble_residual "$(sign_of "$g_installed")" zero
echo "  the cached base face state against the boundary of the installed"
echo "  composition, at assemble_residual: ${g_installed:-missing}"

# ---- the same measure with the boundary of the entry composition ---------
stage entry
run_it entry EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition
g_entry="$(gap_at entry evaluate_route_at_assemble_residual)"
check_eq face_state_gap_of_the_entry_composition_boundary \
   "$(sign_of "$g_entry")" nonzero
echo "  the same measure with the boundary of the entry composition left"
echo "  standing: ${g_entry:-missing}"

# ---- reported, not gated: what the face cache is worth to the next row ---
# The report that follows the residual installs a boundary of its own on a
# copy of the conserved array, and the cached face state carries no validity
# (D5b-2).  The continuity row of cell 1 re-assembled afterwards is written
# beside the row the residual gave, in floors of that row's own rounding
# floor, so that the movement is on the record and gates nothing here.
for d in installed entry; do
   python3 - "$WORK/$d/output/boundary_trace.txt" "$d" <<'PY'
import sys
path, tag = sys.argv[1], sys.argv[2]
rows, cur = {}, None
for line in open(path):
    t = line.strip()
    if t.startswith('# massrow'):
        cur = t[10:].strip(); rows[cur] = {}
    elif t.startswith('row ') and cur:
        w = t.split()
        rows[cur][int(w[1])] = (float(w[2]), float(w[5]), float(w[6]))
out = []
for k in ('C_residual_rows', 'F_residual_rows_after_report'):
    if k in rows and 1 in rows[k]:
        R1, sc, fl = rows[k][1]
        out.append(abs(R1)/sc/fl)
if len(out) == 2:
    print('  DIAGNOSTIC %s: cell-1 continuity row %.4E floors at the'
          ' residual, %.4E after the report' % (tag, out[0], out[1]))
PY
done

echo ""
if [ $n_fail -gt 0 ]; then
   echo "boundary_of_the_installed_composition: $n_fail assertion(s) failed"
   exit 1
fi
exit 0

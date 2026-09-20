#!/bin/bash
# THE CELLS ABOVE THE LOWER GHOSTS, HANDED BACK AS THEY CAME: A DECLARED
# DIAGNOSTIC POLICY THAT MOVES NOTHING UNLESS IT IS ASKED FOR (item P1 step 1
# of docs/PLAN_20260919_rev1.md; the measurements are
# docs/lhs1140b_p1_step1b_20260919.md).
#
# WHAT IS UNDER TEST, row by row.
#
#   1  THE KEY IS INERT. With EXHALE_INTERIOR_COMPOSITION_HELD unset the
#      sweep returns the composition it solved on every cell, which is what
#      it has always returned. The row reads the data rows of both state
#      files of a run without the key against a second run without it and
#      asks for byte identity, so that the row fails if the branch is ever
#      entered by default.
#
#   2  WHAT THE KEY HOLDS. With the key set the cells above the two lower
#      ghosts are returned at the composition the sweep was GIVEN, so the
#      state written reproduces the loaded pair on those cells to the round
#      trip of its own printed digits.
#
#   3  WHAT THE KEY DOES NOT HOLD. The two lower ghost cells keep the
#      composition the sweep SOLVED, which is the whole point of the policy:
#      a base row measured after the call carries the ghost's refresh and not
#      the interior's. The row asks that the ghost rows of the written state
#      stand away from the loaded pair's by more than the round trip of row
#      2, so a policy that held the ghosts too would fail it.
#
# THE FIXTURE
#   backup/regression/carrier_model_a_newton/IC/, the hot-Uranus molecular
#   carrier state, evaluated with "Restart intent: stationary evaluate" on a
#   copy. The regression directory is never written to.
#
# THE TOLERANCES
#   Row 1 is an identity and carries none. Row 2 holds the physical cells to
#   1e-12 relative, the round trip of the state file's own printed digits
#   (the writer emits 16 significant digits and the reader rebuilds the
#   fractions from them). Row 3 asks only that the ghost rows exceed that
#   same 1e-12, which is a statement that they were solved and not copied.
#   Both rows read densities above 1e-30 cm^-3 alone: the metal columns of
#   this fixture, which carries no metals, hold the writer's floor of order
#   1e-301 cm^-3, and the last digits of that floor are not a composition.
#
# Usage: interior_composition_held.sh   (EXHALE_EXE selects the binary,
#        EXHALE_TEST_OUT the work directory)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" interior_composition_held
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/interior_composition_held"
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
check_gt() {   # check_gt <name> <measured> <tolerance>
   local v
   v=$(awk -v m="${2:-}" -v t="$3" 'BEGIN {
          if (m == "") { print "missing" }
          else if (m + 0 > t + 0) { print "above" } else { print "within" } }')
   if [ "$v" = "above" ]; then
      echo "PASS $1 measured=${2:-missing} reference=more_than tol=$3"
   else
      echo "FAIL $1 measured=${2:-missing} reference=more_than tol=$3"
      n_fail=$((n_fail + 1))
   fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL interior_composition_held measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Ion_species_IC.txt" ]; then
   echo "FAIL interior_composition_held measured=no_fixture reference=$CASE/IC tol=0"
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

# The data rows of both state files, the '#' lines dropped.
same_state() {   # same_state <dir a> <dir b>
   local f a b
   for f in Hydro_ioniz.txt Ion_species.txt; do
      a="$WORK/$1/output/$f"
      b="$WORK/$2/output/$f"
      if [ ! -s "$a" ] || [ ! -s "$b" ]; then echo missing; return; fi
      if ! diff -q <(grep -v '^#' "$a") <(grep -v '^#' "$b") > /dev/null
      then echo different; return; fi
   done
   echo identical
}

# The largest relative difference between the composition a run wrote and the
# one it loaded, over a band of rows: 'interior' is every row but the two at
# each end, 'ghost' is the two lower ghost rows.
loaded_gap() {   # loaded_gap <dir> <interior|ghost>
   awk -v band="$2" '
      FNR == NR { if ($0 ~ /^#/) next
                  n++; nc[n] = NF
                  for (k = 1; k <= NF; k++) a[n,k] = $k
                  next }
      $0 !~ /^#/ { m++; nb[m] = NF
                   for (k = 1; k <= NF; k++) b[m,k] = $k }
      END {
         lo = 3; hi = n - 2
         if (band == "ghost") { lo = 1; hi = 2 }
         for (i = lo; i <= hi; i++) {
            for (k = 1; k <= nc[i]; k++) {
               x = a[i,k]; y = b[i,k]
               dx = (x < 0 ? -x : x); dy = (y < 0 ? -y : y)
               s = (dx > dy ? dx : dy)
               # A density below this is not a composition: in a run that
               # carries no metals the metal columns hold a floor of order
               # 1e-301 cm^-3, whose last digits say nothing about what the
               # sweep returned.
               if (s <= 1.0e-30) continue
               d = (x > y ? x - y : y - x)/s
               if (d > w) w = d } }
         printf "%.17g\n", w }' \
      "$WORK/$1/output/Ion_species_IC.txt" "$WORK/$1/output/Ion_species.txt"
}

# ---- row 1: the key unset ------------------------------------------------
stage plain
run_it plain
stage plain2
run_it plain2
check_eq interior_hold_key_is_inert "$(same_state plain plain2)" identical

# ---- rows 2 and 3: the key set -------------------------------------------
stage held
run_it held EXHALE_INTERIOR_COMPOSITION_HELD=1
d_int=$(loaded_gap held interior)
d_gho=$(loaded_gap held ghost)
d_free=$(loaded_gap plain interior)
check_le interior_hold_returns_the_given_composition "$d_int" 1.0e-12
check_gt interior_hold_leaves_the_ghost_swept "$d_gho" 1.0e-12
echo "  the cells above the lower ghosts, against the loaded pair:"
echo "  held ${d_int:-missing}, freely refreshed ${d_free:-missing};"
echo "  the two lower ghost rows of the held run ${d_gho:-missing}"

echo ""
if [ $n_fail -gt 0 ]; then
   echo "interior_composition_held: $n_fail assertion(s) failed"
   exit 1
fi
exit 0

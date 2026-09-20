#!/bin/bash
# THE GHOST SEED KEY AND THE GHOST RECORD: THEY MEASURE AND THEY MOVE
# NOTHING (item P1 step 1 of docs/PLAN_20260919_rev1.md; the measurements
# are docs/lhs1140b_p1_step1_20260919.md).
#
# WHAT IS UNDER TEST, row by row.
#
#   1  THE RECORD KEY IS INERT. EXHALE_GHOST_RECORD adds one record per
#      evaluation to output/ghost_record.txt and reads state only. The row
#      reads the data rows of both state files of a run with the key set
#      against a run without it and asks for byte identity.
#
#   2  THE RECORD IS WRITTEN ONLY WHERE IT IS ASKED FOR. With no key set
#      neither output/ghost_record.txt nor a seed file exists.
#
#   3  THE SEED KEY IS INERT WHERE IT STATES THE COMPOSITION THE SWEEP WAS
#      ALREADY ENTERED AT. The first run writes the lower ghost rows it
#      installed (EXHALE_GHOST_SEED_WRITE); the second is seeded from that
#      file. The sweep then solves the same ghost system from the same
#      starting point, so the two states must agree bit for bit. The row
#      also exercises the file reader, which nothing else reaches.
#
#   4  THE SEED KEY IS INERT WHERE IT HAS NOTHING TO STATE. Seeded with
#      the_previous_sweep_ghost on the evaluate route, which runs one sweep
#      and therefore has no previous one, the run falls back to the entry
#      state and writes the same state.
#
#   5  THE RECORD'S NUMBERS ARE THE STATE'S. The electron count of lower
#      ghost cell 0 in the record is compared with the count made of the
#      Ion_species.txt row the same run wrote: H II, He II, twice He III and
#      the three molecular ions H2+, H3+ and HeH+, each with its charge.
#      That sum is calc_ne's policy, and the row is what item R2 of the plan
#      asks of the solved ghost electron count. The fixture carries no
#      metals (no metals.inp), so no metal stage enters the sum.
#
# THE FIXTURE
#   backup/regression/carrier_model_a_newton/IC/, the hot-Uranus molecular
#   carrier state, evaluated with "Restart intent: stationary evaluate" on a
#   copy. The regression directory is never written to.
#
# THE TOLERANCES
#   Rows 1, 3 and 4 are identities and carry none. Row 5 holds the record
#   against the written state to 1e-12 relative, the round trip of the
#   state's own printed digits (the file carries 16 significant digits, the
#   sum is over seven terms of widely different size, and the record is
#   written before the state file is formatted).
#
# Usage: ghost_seed_and_ghost_record.sh   (EXHALE_EXE selects the binary,
#        EXHALE_TEST_OUT the work directory)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" ghost_seed_and_ghost_record
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/ghost_seed_and_ghost_record"
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
   echo "FAIL ghost_seed_and_ghost_record measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -s "$CASE/IC/Hydro_ioniz_IC.txt" ]; then
   echo "FAIL ghost_seed_and_ghost_record measured=no_fixture reference=$CASE/IC tol=0"
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

# ---- row 1 and row 2: the record key --------------------------------------
stage plain
run_it plain
stage recorded
run_it recorded EXHALE_GHOST_RECORD=1

check_eq ghost_record_key_is_inert "$(same_state plain recorded)" identical

if [ -e "$WORK/plain/output/ghost_record.txt" ]; then
   present=written
else
   present=absent
fi
check_eq ghost_record_absent_without_the_key "$present" absent

if [ ! -s "$WORK/recorded/output/ghost_record.txt" ]; then
   echo "FAIL ghost_record_written_with_the_key measured=missing reference=written tol=0"
   n_fail=$((n_fail + 1))
else
   echo "PASS ghost_record_written_with_the_key measured=written reference=written tol=0"
fi

# ---- row 3: the seed key on the composition already installed -------------
stage seedwrite
run_it seedwrite EXHALE_GHOST_SEED_WRITE=ghost_seed.txt
if [ ! -s "$WORK/seedwrite/ghost_seed.txt" ]; then
   echo "FAIL ghost_seed_rows_written measured=missing reference=written tol=0"
   n_fail=$((n_fail + 1))
else
   echo "PASS ghost_seed_rows_written measured=written reference=written tol=0"
   stage seedback
   run_it seedback \
      EXHALE_GHOST_COMPOSITION_SEED="$WORK/seedwrite/ghost_seed.txt"
   if grep -q 'was not read' "$WORK/seedback/run.log"; then
      echo "FAIL ghost_seed_file_is_read measured=not_read reference=read tol=0"
      n_fail=$((n_fail + 1))
   else
      echo "PASS ghost_seed_file_is_read measured=read reference=read tol=0"
   fi
   check_eq ghost_seed_of_the_installed_rows_is_inert \
      "$(same_state plain seedback)" identical
fi

# ---- row 4: a named source with nothing stored ----------------------------
stage seedempty
run_it seedempty EXHALE_GHOST_COMPOSITION_SEED=the_previous_sweep_ghost
check_eq ghost_seed_without_a_previous_sweep_is_inert \
   "$(same_state plain seedempty)" identical

# ---- row 5: the record against the state written --------------------------
# Lower ghost cell 0 is the second data row of Ion_species.txt (the rows run
# 1-Ng ... N+Ng).  The species columns of the '# columns' header are, in
# order, r HI HII HeI HeII HeIII HeITR then the metal stages then H2 H2p H3p
# HeHp, so the electron count of the row is HII + HeII + 2 HeIII + H2p +
# H3p + HeHp with no metals in this fixture.
ne_state=$(awk 'BEGIN { OFMT = "%.17g" }
   !/^#/ { n++; if (n == 2) {
      print $3 + $5 + 2*$6 + $(NF-2) + $(NF-1) + $NF; exit } }' \
   "$WORK/recorded/output/Ion_species.txt")
ne_record=$(awk '$1 == "ghost" && $3 == 0 { print $8; exit }' \
   "$WORK/recorded/output/ghost_record.txt")
d_ne=$(awk -v a="${ne_state:-}" -v b="${ne_record:-}" 'BEGIN { OFMT = "%.17g" 
      if (a == "" || b == "") { print "" }
      else { s = (a > b ? a : b); if (s <= 0) print 0;
             else print (a > b ? a - b : b - a)/s } }')
check_le ghost_record_electron_count_is_the_states "$d_ne" 1.0e-12
echo "  ghost cell 0 electron density [cm^-3]: record ${ne_record:-missing},"
echo "  counted from the written composition ${ne_state:-missing}"

echo ""
if [ $n_fail -gt 0 ]; then
   echo "ghost_seed_and_ghost_record: $n_fail assertion(s) failed"
   exit 1
fi
exit 0

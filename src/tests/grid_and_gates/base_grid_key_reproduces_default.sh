#!/bin/bash
# The resolved base cell width reproduces the grid, and so does the default
# spelled as a reader writes it.
#
# QUANTITY UNDER TEST
#   The cell centers r(j) of the Mixed grid, column 1 of output/Hydro_ioniz.txt,
#   written list-directed and therefore carrying the full double.
#
# WHAT IS RUN (three one-step runs of backup/regression/roundtrip, on copies)
#   A  as shipped: no "Base grid [dr,cells]:" key, so the run builds its grid
#      on the code's default width.
#   B  the same input with "Base grid [dr,cells]: <w> <n>", where <w> and <n>
#      are read back out of A's EXHALE_resolved.out (base_cell_width_Rp and
#      base_uniform_cells).  Nothing is rounded and nothing is retyped.
#   C  the same input with "Base grid [dr,cells]: 2.0e-4 <n>" -- the width as a
#      reader would spell the default. The default is the double 2.0d-4, so
#      this is the default.
#   P  the same input with the line src/utils/pin_base_grid.py wrote into the
#      inputs run before 2026-09-19, "Base grid [dr,cells]:
#      1.9999999494757503e-4 50": the width those runs were made on.
#
# ASSERTIONS
#   1. A's record says the width came from the default, B's says it came from
#      the key.  Provenance is a flag, not a comparison of values: B states the
#      default's own digits and must still be reported as a key.
#   2. B's r column is identical to A's BYTE FOR BYTE.  This is what makes the
#      record a reproduction instruction: a run can be put back on its grid by
#      stating the width the record carries.
#   3. C's r column is identical to A's BYTE FOR BYTE: an absent key and the
#      key stating 2.0e-4 build one grid.
#   4. P's r column is NOT identical to A's: the pinned width keeps a pinned
#      input on the grid its stored results were written on, which is not
#      the present default grid.  The measured maximum relative difference
#      of the cell centers is REPORTED as an observation of this grid; it
#      depends on the grid the width difference is rescaled onto.  That P
#      builds the pre-2026-09-19 default grid to the bit is asserted on the
#      coordinates themselves by base_grid_pinned_width.f90.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" base_grid_key_reproduces_default
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
CASE="$ROOT/backup/regression/roundtrip"

nfail=0
say_pass() { echo "PASS $1 measured=$2 reference=$3 tol=$4"; }
say_fail() { echo "FAIL $1 measured=$2 reference=$3 tol=$4"; nfail=$((nfail+1)); }

if [ ! -x "$EXE" ]; then
   say_fail base_grid_key_binary no_binary "$EXE" 0
   exit 1
fi

for d in bgA bgB bgC bgP; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d"
   cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/$d/"
done
# A must be the shipped configuration: no key of this name anywhere in it.
sed -i '/^Base grid/d' "$WORK/bgA/input.inp" "$WORK/bgB/input.inp" \
        "$WORK/bgC/input.inp" "$WORK/bgP/input.inp"

run_one() {   # $1 = directory
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   local rc=$?
   # Exit 2 is a refused stationary claim, which still writes every output;
   # only 1 or a signal is a failed run here.
   if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
      say_fail "base_grid_key_run_$1" "exit_$rc" exit_0 0
      return 1
   fi
   return 0
}

field() {     # $1 = directory, $2 = key of EXHALE_resolved.out
   awk -v k="$2" '$1 == k { print $2 }' "$WORK/$1/EXHALE_resolved.out"
}

radii() {     # $1 = directory: the r column, one value per line, verbatim
   grep -v '^#' "$WORK/$1/output/Hydro_ioniz.txt" | awk 'NF > 1 { print $1 }'
}

run_one bgA || exit 1
W="$(field bgA base_cell_width_Rp)"
NCELL="$(field bgA base_uniform_cells)"
SRC_A="$(field bgA base_cell_width_source)"
if [ -z "$W" ] || [ -z "$NCELL" ] || [ -z "$SRC_A" ]; then
   say_fail base_grid_resolved_record_present missing_rows \
            base_cell_width_Rp+base_uniform_cells+base_cell_width_source 0
   exit 1
fi
echo "  A resolved: width=$W cells=$NCELL source=$SRC_A"

if [ "$SRC_A" = "default" ]; then
   say_pass base_grid_width_source_without_key default default 0
else
   say_fail base_grid_width_source_without_key "$SRC_A" default 0
fi

echo "Base grid [dr,cells]: $W $NCELL" >> "$WORK/bgB/input.inp"
echo "Base grid [dr,cells]: 2.0e-4 $NCELL" >> "$WORK/bgC/input.inp"
echo "Base grid [dr,cells]: 1.9999999494757503e-4 50" >> "$WORK/bgP/input.inp"
run_one bgB || exit 1
run_one bgC || exit 1
run_one bgP || exit 1

SRC_B="$(field bgB base_cell_width_source)"
if [ "$SRC_B" = "key" ]; then
   say_pass base_grid_width_source_with_key key key 0
else
   say_fail base_grid_width_source_with_key "$SRC_B" key 0
fi

radii bgA > "$WORK/bgA/r.txt"
radii bgB > "$WORK/bgB/r.txt"
radii bgC > "$WORK/bgC/r.txt"
radii bgP > "$WORK/bgP/r.txt"
if [ ! -s "$WORK/bgA/r.txt" ]; then
   say_fail base_grid_radius_column empty non_empty 0
   exit 1
fi

if cmp -s "$WORK/bgA/r.txt" "$WORK/bgB/r.txt"; then
   say_pass base_grid_resolved_width_reproduces_default identical identical 0
else
   say_fail base_grid_resolved_width_reproduces_default differs identical 0
   diff "$WORK/bgA/r.txt" "$WORK/bgB/r.txt" | head -n 4
fi

if cmp -s "$WORK/bgA/r.txt" "$WORK/bgC/r.txt"; then
   say_pass base_grid_spelled_default_is_the_default identical identical 0
else
   say_fail base_grid_spelled_default_is_the_default differs identical 0
   diff "$WORK/bgA/r.txt" "$WORK/bgC/r.txt" | head -n 4
fi

# The measured displacement of the pinned grid, reported and not gated.
paste "$WORK/bgA/r.txt" "$WORK/bgP/r.txt" | awk '
   { if ($1 != 0) { d = ($2 - $1)/$1; if (d < 0) d = -d;
                    if (d > m) { m = d; j = NR } } }
   END { printf "  A vs P: max relative difference of the cell centers %.6e at row %d of %d\n", m, j, NR }'
if cmp -s "$WORK/bgA/r.txt" "$WORK/bgP/r.txt"; then
   say_fail base_grid_pinned_width_keeps_its_own_grid identical differs 0
else
   say_pass base_grid_pinned_width_keeps_its_own_grid differs differs 0
fi

if [ $nfail -gt 0 ]; then exit 1; fi
exit 0

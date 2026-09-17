#!/usr/bin/env bash
# T1 of docs/lhs1140b_stationary_L22_step3_design_20260916.md section 7:
# the coupling of the coupled block, measured on a frozen state against a
# central difference of the same map.
#
# WHAT IS RUN.  EXHALE.x on a scratch copy of the frozen state
# LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55 (its input.inp,
# its base.inp and its output/, with the spectrum path made absolute),
# twice: once with a direction concentrated on the CARRIER unknown of cells
# 300 to 312 and once with a direction on the ENERGY unknown of the same
# cells.  EXHALE_COUPLED_JAC_ACTION registers the block's rows for the
# measurement, writes the file and stops the run without a step, so neither
# invocation takes a solve.
#
# WHAT THE FILE HOLDS, row by row for every unknown of the block (mass,
# momentum, energy and each species row): the action the solve uses
# (jacobian_action_of_direction, a forward difference of the full residual
# at the probe length the route sets), a CENTRAL difference of the same
# residual at that length and at a tenth of it, and the row of the banded
# preconditioner on the same direction.
#
# THE ROWS ASSERTED
#   the_action_of_a_carrier_direction_moves_the_energy_row
#   the_action_of_a_hydrodynamic_direction_moves_the_carrier_row
#       the two cross terms of section 3.2: the composition reaches the
#       energy row through the caloric equation of state and through the
#       heating and the cooling of the sweep, and the hydrodynamic state
#       reaches the carrier row through the advective face flux, the
#       temperature of every rate, the diffusivities and the photolysis.
#   the_assembled_action_reproduces_the_central_difference
#       over the rows the direction MOVES, at or below 1.0e-2.  The
#       tolerance is stated before the measurement and is the order of the
#       1.542e-03 that L22 step 2b section 5 measured for the carrier
#       operator's own assembled action and accepted.  A row whose two
#       sides both stand at the rounding floor carries no derivative to
#       compare (the mass row of a cell does not move with the composition
#       at all), so the mean over EVERY row of the support is reported
#       beside it and is not the gate.
#   the_band_holds_the_two_cell_reconstruction_entry
#       d res_c(j+2)/d f_c(j), the entry the carrier relaxation's
#       block-tridiagonal matrix has no place for, stands at a flat
#       distance of 2*nvar_jac inside a band of kl_jac = 3*nvar_jac - 1.
#       Every one of those entries must be nonzero.
#   the_species_rows_scale_is_the_certifications
#       the carrier row of eval_residual, divided by the code time scale
#       tscale_code it was converted by, is the row carrier_steady_residual
#       returns, bitwise, at every cell of the direction.
#   the_outer_ghost_rule_is_the_same_in_both_evaluations
#       the same comparison at cells N-1 and N.  EXPECTED TO FAIL while two
#       outer ghost rules are in force (L22 step 2b section 1); it is
#       written to be read, not widened.
#
# RED BEFORE THE INCREMENT: the entry text has no
# EXHALE_COUPLED_JAC_ACTION, so no file is written and no row of this suite
# can be stated against it.
#
# EXHALE_EXE selects the binary, EXHALE_TEST_OUT the working directory.
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/coupled_block_jacobian}"
CASE="${EXHALE_JAC_CASE:-$ROOT/LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55}"
SED="${EXHALE_JAC_SED:-$ROOT/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt}"
LO=300
HI=312
TOL=1.0e-2

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL coupled_block_jacobian_binary no_binary "$EXE" 0
   exit 1
fi
if [ ! -f "$CASE/input.inp" ] || [ ! -f "$CASE/output/Hydro_ioniz.txt" ]; then
   verdict FAIL coupled_block_jacobian_case missing "$CASE" 0
   exit 1
fi

stage() {  # stage <dir>
   rm -rf "$OUT/$1"
   mkdir -p "$OUT/$1/output"
   cp -f "$CASE/input.inp" "$CASE/base.inp" "$OUT/$1/"
   cp -f "$CASE"/output/*.txt "$OUT/$1/output/"
   sed -i "s|^Spectrum file: .*|Spectrum file: $SED|" "$OUT/$1/input.inp"
   sed -i "s|^Restart intent: .*|Restart intent: stationary|" "$OUT/$1/input.inp"
}

run_direction() {  # run_direction <dir> <carrier|energy>
   stage "$1"
   # THE MEASUREMENT STOPS THE RUN. A binary that does not carry the key
   # ignores it and goes on solving, which is the state of the entry text,
   # so the invocation is bounded and a binary that does not stop is a
   # failed row rather than a suite that never returns.
   ( cd "$OUT/$1" && OMP_NUM_THREADS=1 \
        EXHALE_COUPLED_JAC_ACTION="$OUT/$1/jac.txt" \
        EXHALE_COUPLED_JAC_CELLS="$LO,$HI" \
        EXHALE_COUPLED_JAC_DIR="$2" \
        timeout "${EXHALE_JAC_TIMEOUT:-1800}" "$EXE" > run.log 2>&1 )
}

run_direction dcarrier carrier
run_direction denergy  energy

for d in dcarrier denergy; do
   if [ ! -s "$OUT/$d/jac.txt" ]; then
      verdict FAIL coupled_block_jacobian_written no_file "$OUT/$d/jac.txt" 0
      tail -n 5 "$OUT/$d/run.log" | sed 's/^/     /'
      exit 1
   fi
done

CF="$OUT/dcarrier/jac.txt"
EF="$OUT/denergy/jac.txt"

# ---- the two cross terms ----
# The smallest magnitude the named row kind takes over the cells of the
# direction; a zero anywhere in that window is no coupling.
smallest_row() {  # smallest_row <file> <kind text>
   awk -v lo="$LO" -v hi="$HI" -v kind="$2" '
      $1 ~ /^[0-9]+$/ && $1 >= lo && $1 <= hi {
         line = $0
         if (index(line, kind) == 0) next
         v = $(NF-3); if (v < 0) v = -v
         if (n == 0 || v < m) m = v
         n++
      }
      END { if (n == 0) print "none"; else printf "%.6e", m }' "$1"
}
e_in_carrier=$(smallest_row "$CF" "energy of cell")
c_in_energy=$(smallest_row "$EF" "carrier")
nz() { awk -v v="$1" 'BEGIN { exit !(v+0 > 0) }'; }
if [ "$e_in_carrier" != none ] && nz "$e_in_carrier"; then
   verdict PASS the_action_of_a_carrier_direction_moves_the_energy_row \
           "$e_in_carrier" nonzero 0
else
   verdict FAIL the_action_of_a_carrier_direction_moves_the_energy_row \
           "$e_in_carrier" nonzero 0
fi
if [ "$c_in_energy" != none ] && nz "$c_in_energy"; then
   verdict PASS the_action_of_a_hydrodynamic_direction_moves_the_carrier_row \
           "$c_in_energy" nonzero 0
else
   verdict FAIL the_action_of_a_hydrodynamic_direction_moves_the_carrier_row \
           "$c_in_energy" nonzero 0
fi

# ---- the action against the central difference ----
for pair in "carrier:$CF" "energy:$EF"; do
   name="${pair%%:*}"; f="${pair#*:}"
   moved=$(awk '/^# mean relative error over the rows the direction moves/ \
                { for (i = 1; i <= NF; i++) if ($i == "moves") print $(i+1) }' "$f")
   allrows=$(awk '/^# mean relative error of the action against the central difference over every row of the support/ \
                { for (i = 1; i <= NF; i++) if ($i == "support") print $(i+1) }' "$f")
   echo "  $name direction: mean over every row of the support $allrows," \
        "over the rows the direction moves $moved"
   if [ -n "$moved" ] && awk -v v="$moved" -v t="$TOL" 'BEGIN { exit !(v+0 <= t+0) }'; then
      verdict PASS "the_assembled_action_reproduces_the_central_difference_${name}" \
              "$moved" 0 "$TOL"
   else
      verdict FAIL "the_assembled_action_reproduces_the_central_difference_${name}" \
              "${moved:-unread}" 0 "$TOL"
   fi
done

# ---- the band holds the two-cell reconstruction entry ----
zeros=$(awk '/^# two-cell reconstruction entries that are exactly zero/ { print $9 }' "$CF")
total=$(awk '/^# two-cell reconstruction entries that are exactly zero/ { print $11 }' "$CF")
if [ "${zeros:-1}" = 0 ] && [ "${total:-0}" -gt 0 ]; then
   verdict PASS the_band_holds_the_two_cell_reconstruction_entry \
           "0_of_$total" "0_of_$total" 0
else
   verdict FAIL the_band_holds_the_two_cell_reconstruction_entry \
           "${zeros:-unread}_of_${total:-unread}" "0_of_$total" 0
fi

# ---- the species rows' scale, and the outer ghost rule ----
bitwise() {  # bitwise <file> <header text>; the number of rows that are not bitwise equal
   awk -v k="$1" '$0 ~ k { if ($NF != "T") n++ } END { print n+0 }' "$2"
}
nrows_scale=$(awk '/^# species row scale, cell/ { n++ } END { print n+0 }' "$CF")
bad_scale=$(bitwise "^# species row scale, cell" "$CF")
if [ "$bad_scale" = 0 ] && [ "$nrows_scale" -gt 0 ]; then
   verdict PASS the_species_rows_scale_is_the_certifications \
           "0_of_$nrows_scale" "0_of_$nrows_scale" 0
else
   verdict FAIL the_species_rows_scale_is_the_certifications \
           "${bad_scale}_of_$nrows_scale" "0_of_$nrows_scale" 0
fi
nrows_ghost=$(awk '/^# outer ghost rule, cell/ { n++ } END { print n+0 }' "$CF")
bad_ghost=$(bitwise "^# outer ghost rule, cell" "$CF")
if [ "$bad_ghost" = 0 ] && [ "$nrows_ghost" -gt 0 ]; then
   verdict PASS the_outer_ghost_rule_is_the_same_in_both_evaluations \
           "0_of_$nrows_ghost" "0_of_$nrows_ghost" 0
else
   verdict FAIL the_outer_ghost_rule_is_the_same_in_both_evaluations \
           "${bad_ghost}_of_$nrows_ghost" "0_of_$nrows_ghost" 0
fi

# ---- the three numbers the coupled route was set aside for ----
grep -h '^# fraction to the boundary of the species box' "$CF" "$EF" \
   | sed 's/^# /  /'
grep -h '^# the action was sampled' "$CF" "$EF" | sed 's/^# /  /'

if [ "$n_fail" -ne 0 ]; then
   echo "coupled_block_jacobian: FAILED"
else
   echo "coupled_block_jacobian: PASSED"
fi
exit $(( n_fail > 0 ))

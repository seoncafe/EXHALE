#!/usr/bin/env bash
# T1 of docs/lhs1140b_stationary_L22_step3_design_20260916.md section 7:
# the coupling of the coupled block, measured on a frozen state against a
# central difference of the same map.
#
# WHAT IS RUN.  EXHALE.x on a scratch copy of the frozen state
# LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55 (its input.inp,
# its base.inp and its output/, with the spectrum path made absolute), four
# times: a direction concentrated on the CARRIER unknown of cells 300 to
# 312 and a direction on the ENERGY unknown of the same cells, each of them
# at the probe arc the route sets and again at ARC_LONG times that arc
# (EXHALE_JV_PROBE_ARC multiplies probe_length_of_the_jacobian_action, so
# the second pair needs no code of its own).  EXHALE_COUPLED_JAC_ACTION
# registers the block's rows for the measurement, writes the file and stops
# the run without a step, so no invocation takes a solve.
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
#   the_distance_to_the_central_difference_scales_as_the_first_power_of_the_arc
#       The mean, over the rows the direction MOVES, of the distance
#       between the action the solve uses and a central difference of the
#       same residual AT THE SAME ARC is, to leading order,
#       (h/2)|F''|/|F'| plus a rounding term of order eps|F|/h: the
#       truncation of the FORWARD difference, and not a property of the
#       assembly.  Its value at one arc therefore states the arc and
#       nothing else (MEASURED over four decades, energy direction:
#       1.2390e-06, 1.2388e-05, 1.2396e-04, 1.2381e-03, 1.2379e-02 at arc
#       multipliers 1e-4 to 1; docs/lhs1140b_coupled_action_20260921.md
#       section 10.2).  What is arc-free is the POWER: the distance moves
#       with the first power of the arc, upward on the truncation branch
#       and downward on the rounding branch, and a constant error of the
#       assembly is exactly what cannot do that.  So the statistic is taken
#       at the production arc and at ARC_LONG times it and
#       |d ln e / d ln h| is gated at 1 within EXPONENT_TOL, per direction.
#       A row whose two sides both stand at the rounding floor carries no
#       derivative to compare (the mass row of a cell does not move with
#       the composition at all), so the mean over EVERY row of the support
#       is reported beside it and is not the gate.
#   the_band_holds_the_two_cell_reconstruction_entry
#       d res_c(j+2)/d f_c(j), the entry the carrier relaxation's
#       block-tridiagonal matrix has no place for, stands at a flat
#       distance of 2*nvar_jac inside a band of kl_jac = 3*nvar_jac - 1.
#       Every one of those entries must be nonzero.
#   the_species_rows_scale_is_the_certifications
#       the carrier row of eval_residual is the row carrier_steady_residual
#       returns converted by the code time scale tscale_code, bitwise, at
#       every cell of the direction.  Both sides are read in the units of
#       the residual vector, the certification's row multiplied by
#       tscale_code and never the residual divided back: for a double x and
#       a scale t, (x*t)/t is not in general x again, so a bitwise
#       statement made across that round trip would fail on the rounding of
#       the conversion alone.
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
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" coupled_block_jacobian
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/coupled_block_jacobian}"
CASE="${EXHALE_JAC_CASE:-$ROOT/LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55}"
SED="${EXHALE_JAC_SED:-$ROOT/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt}"
LO=300
HI=312
# THE TWO ARCS, and the tolerance on the power between them.  ARC_LONG is a
# decade, which is far enough that the two statistics are separated by many
# digits and near enough that neither direction crosses the turn of its own
# V between them (the energy direction is on the truncation branch from
# 1e-4 to 10 and the carrier one on the rounding branch from 1e-3 to 10,
# MEASURED, docs/lhs1140b_coupled_action_20260921.md section 10.2).
# EXPONENT_TOL 0.30 is stated from that ladder, before this suite was run
# with it: it admits every decade of it, the widest being 0.289 at the
# carrier decade 1e-1 to 1 where the V turns, while a constant error of the
# assembly reads an exponent of 0 and is refused by 0.70.  The curvature it
# admits is a factor 2.0 in the ratio over one decade.
ARC_LONG=10
EXPONENT_TOL=0.30

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

run_direction() {  # run_direction <dir> <carrier|energy> <arc multiplier>
   stage "$1"
   # THE MEASUREMENT STOPS THE RUN. A binary that does not carry the key
   # ignores it and goes on solving, which is the state of the entry text,
   # so the invocation is bounded and a binary that does not stop is a
   # failed row rather than a suite that never returns.
   ( cd "$OUT/$1" && OMP_NUM_THREADS=1 \
        EXHALE_COUPLED_JAC_ACTION="$OUT/$1/jac.txt" \
        EXHALE_COUPLED_JAC_CELLS="$LO,$HI" \
        EXHALE_COUPLED_JAC_DIR="$2" \
        EXHALE_JV_PROBE_ARC="$3" \
        timeout "${EXHALE_JAC_TIMEOUT:-1800}" "$EXE" > run.log 2>&1 )
}

run_direction dcarrier      carrier 1
run_direction dcarrier_long carrier "$ARC_LONG"
run_direction denergy       energy  1
run_direction denergy_long  energy  "$ARC_LONG"

for d in dcarrier dcarrier_long denergy denergy_long; do
   if [ ! -s "$OUT/$d/jac.txt" ]; then
      verdict FAIL coupled_block_jacobian_written no_file "$OUT/$d/jac.txt" 0
      tail -n 5 "$OUT/$d/run.log" | sed 's/^/     /'
      exit 1
   fi
done

# Every row but the power below is read at the PRODUCTION arc, so that the
# suite states the cross terms, the band, the species-row scale and the
# outer ghost rule of the arc the solve actually takes.
CF="$OUT/dcarrier/jac.txt"
EF="$OUT/denergy/jac.txt"
CFL="$OUT/dcarrier_long/jac.txt"
EFL="$OUT/denergy_long/jac.txt"

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

# ---- the action against the central difference, against the arc ----
# The number the file carries at one arc is that arc's truncation, so what
# is gated is the power of the arc it moves with and never the number
# itself.  The two statistics and the two probe scalars are printed, so
# that a later drift can be attributed to the arc or to the curvature.
moved_mean() {  # moved_mean <file>; the mean over the rows the direction moves
   awk '/^# mean relative error over the rows the direction moves/ \
        { for (i = 1; i <= NF; i++) if ($i == "moves") print $(i+1) }' "$1"
}
support_mean() {  # support_mean <file>; the same over every row of the support
   awk '/^# mean relative error of the action against the central difference over every row of the support/ \
        { for (i = 1; i <= NF; i++) if ($i == "support") print $(i+1) }' "$1"
}
probe_arc() {  # probe_arc <file>; the probe scalar the file was written at
   # The writer calls it `# probe scalar`, because the number is the scalar
   # the direction is multiplied by and the displacement is that times ||v||.
   # `# probe length` is the wording of files written before 2026-09-22 and is
   # read too, so an archived jac.txt still parses.
   awk '/^# probe scalar/ { print $4 }
        /^# probe length/ { print $4 }' "$1" | head -n 1
}
for pair in "carrier:$CF:$CFL" "energy:$EF:$EFL"; do
   name="${pair%%:*}"; rest="${pair#*:}"
   f="${rest%%:*}"; g="${rest#*:}"
   moved=$(moved_mean "$f");       moved_long=$(moved_mean "$g")
   arc=$(probe_arc "$f");          arc_long=$(probe_arc "$g")
   allrows=$(support_mean "$f")
   echo "  $name direction at probe scalar $arc: mean over every row of" \
        "the support $allrows, over the rows the direction moves $moved"
   echo "  $name direction at probe scalar $arc_long (a factor $ARC_LONG" \
        "above): over the rows the direction moves $moved_long"
   # The power, over an interval of ARC_LONG in the arc.  Both statistics
   # have to be positive numbers for a power to exist at all: a zero is a
   # forward difference that is not a forward difference, and it refuses
   # here rather than reaching the logarithm.
   power=$(awk -v a="$moved" -v b="$moved_long" -v s="$ARC_LONG" \
      'BEGIN { if (a+0 <= 0 || b+0 <= 0) { print "unread"; exit }
               p = log((b+0)/(a+0))/log(s+0); if (p < 0) p = -p
               printf "%.4f", p }')
   if [ "$power" != unread ] && \
      awk -v p="$power" -v t="$EXPONENT_TOL" \
          'BEGIN { d = p - 1; if (d < 0) d = -d; exit !(d <= t+0) }'; then
      verdict PASS \
         "the_distance_to_the_central_difference_scales_as_the_first_power_of_the_arc_${name}" \
         "$power" 1 "$EXPONENT_TOL"
   else
      verdict FAIL \
         "the_distance_to_the_central_difference_scales_as_the_first_power_of_the_arc_${name}" \
         "$power" 1 "$EXPONENT_TOL"
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

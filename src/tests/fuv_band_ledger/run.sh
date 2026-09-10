#!/bin/bash
# The FUV band ledger of output/FUV_bands.txt says what it measures.
#
# QUANTITY UNDER TEST
#   the three statements the gate-G4 trailer of output/FUV_bands.txt makes
#   about each of the four FUV photolysis bands, written by
#   write_output.f90 (write_oxygen_chemistry):
#
#   (a) THE RATES AND THE COLUMNS DESCRIBE THE SAME ABSORBERS.
#       cont_absorbed_ph is the H2O and OH rates applied to the absorbers
#       the optical-depth columns record, cont_beam_loss the same cells
#       written as the beam's loss between their two faces,
#       (j_H2O/sigma_H2O) dtau.  The two agree as long as the band depth is
#       exactly sigma_H2O dN_H2O + sigma_OH dN_OH, so rel_diff is round-off
#       and a residual says a column has gone negative or the depths and the
#       columns are not one state.  It holds whatever the rate form is.
#
#   (b) BEAM BUDGET. unrated_frac is the share of the beam's loss that no
#       rate accounts for.  For B3 and B4 there is no line absorber and
#       every absorber is rated, so it is a CLOSURE residual, and it is the
#       test of the cell-mean rate form of water_photolysis_rate: the cell
#       sums telescope to N_b (1 - exp(-tau)) exactly, at any grid spacing,
#       only because the rate is the mean over the cell.  A face-value rate
#       breaks it one-signed.  For B2 the band IS the H I Ly-alpha line and
#       the H I that scatters the stellar line out of the beam dissociates
#       nothing, so it is the SHARE H I takes and is a real number near 1.
#       For LW it is a CLOSURE residual as well, since 2026-09-06: the
#       band is 912-1201 A on both sides, and both the rated absorption and
#       the beam's T_line come from the column integral of sigma_pump =
#       sigma_diss/p_eff of one level-resolved table.  It is not round-off,
#       and the reason is stated rather than absorbed into a tolerance:
#       this case's star-ward H2 column reaches 6.4e21 cm^-2, 1.28 times
#       the top of the table's column axis, where the table is clamped at
#       its edge cross section; A then grows past 1 and the beam is capped
#       at 1 while the rate is not, so the rated side stands above the beam
#       loss by the part of the column the table cannot describe.  The
#       tolerance below is set at 5e-2 for that reason.  MEASURED on this
#       case: -3.62e-1 before the band merge, -1.18e-2 after.
#       lyman_werner.f90 sec. 2g and 3a and
#       src/tests/physics_probe/lyman_werner_cell_mean.f90 sec. 6 are the
#       measurement of the closure itself, on a column inside the axis,
#       where it is 1e-6.
#
#   (c) the trailer names all of that in its own header lines, so a reader
#       of the file is not left to guess which column is a residual.
#
# WHAT IS RUN
#   one short run on a copy of backup/regression/oxygen_chemistry, the only
#   named case that turns the oxygen chemistry on and therefore the only one
#   that writes output/FUV_bands.txt.  EXHALE_MAXSTEPS keeps it to a few
#   seconds; the ledger is a property of whatever state is written, so a
#   relaxation snapshot tests it exactly as a converged state would.
#
# TOLERANCES
#   1e-6 on every quantity the trailer calls round-off, and 5e-2 on the LW
#   beam budget, for the reason paragraph (b) gives.  Bands whose flux is
#   zero contribute zeros and are skipped, not passed by a loose tolerance.
#
# EXPECTED BEFORE THE CHANGE OF ITEM LEDGER-FUV: RED.  The trailer compared
# an absorbed total that carries the B2 Ly-alpha line factor and the LW H2
# share against a closed form that carries neither, and its one rel_diff
# column read 4.5e-1 (LW), 2.9e-8 (B1), 1.0 (B2), 4.2e-7 (B3) and 6.7e-6
# (B4) on this case at 45 steps.  After: GREEN.  Band B1 was merged into
# LW on 2026-09-06 (item LW-NORM-B) and is no longer a row of the trailer.
#
# Usage: src/tests/fuv_band_ledger/run.sh
#        EXHALE_EXE       binary to run (default $ROOT/EXHALE.x)
#        EXHALE_TEST_WORK directory to run in (default build/tests/fuv_band_ledger)
#        EXHALE_TEST_STEPS steps to march (default 5)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="${EXHALE_TEST_WORK:-$ROOT/build/tests/fuv_band_ledger}"
STEPS="${EXHALE_TEST_STEPS:-5}"
CASE="$ROOT/backup/regression/oxygen_chemistry"

nfail=0
say() {  # say PASS|FAIL name measured reference tol
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && nfail=$((nfail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   say FAIL fuv_band_ledger_binary no_binary "$EXE" 0
   exit 1
fi

rm -rf "$WORK/run"
mkdir -p "$WORK/run/output"
cp "$CASE/input.inp" "$WORK/run/"
[ -f "$CASE/metals.inp" ] && cp "$CASE/metals.inp" "$WORK/run/"
( cd "$WORK/run" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS="$STEPS" "$EXE" \
     > run.log 2>&1 )
rc=$?
LED="$WORK/run/output/FUV_bands.txt"
if [ "$rc" -ne 0 ] || [ ! -f "$LED" ]; then
   say FAIL fuv_band_ledger_run "exit=$rc" "exit=0 with $LED" 0
   sed -n '$p' "$WORK/run/run.log" 2>/dev/null
   exit 1
fi
say PASS fuv_band_ledger_run "exit=0" "exit=0" 0

# The trailer is three labeled blocks; read each band row of each block.
band_row() {  # band_row <header regexp> <band>
   awk -v hdr="$1" -v b="$2" '
      $0 ~ hdr        { on = 1; next }
      on && $2 == b   { $1 = ""; $2 = ""; print; exit }
      on && $0 ~ /^# \(/ { exit }
   ' "$LED"
}

abs_lt() {  # abs_lt <value> <tol> -> 0 when |value| < tol
   awk -v v="$1" -v t="$2" 'BEGIN { exit !((v < 0 ? -v : v) < t) }'
}
gt() {  # gt <value> <bound>
   awk -v v="$1" -v t="$2" 'BEGIN { exit !(v > t) }'
}

# ---- (a) the closure of the rate form, every band -----------------------
for b in LW B2 B3 B4; do
   row=$(band_row '^# band cont_absorbed_ph' "$b")
   set -- $row
   if [ $# -ne 3 ]; then
      say FAIL "closure_row_$b" "$# fields" "3 fields" 0
      continue
   fi
   loss="$2"; rel="$3"
   if awk -v x="$loss" 'BEGIN { exit !(x == 0) }'; then
      say PASS "closure_${b}_no_flux" "beam_loss=0" "skipped" 0
   elif abs_lt "$rel" 1e-6; then
      say PASS "closure_$b" "$rel" 0 1e-6
   else
      say FAIL "closure_$b" "$rel" 0 1e-6
   fi
done

# ---- (b) the beam budget: round-off where the header says closure -------
for b in B3 B4; do
   row=$(band_row '^# band rated_ph beam_loss_ph' "$b")
   set -- $row
   if [ $# -ne 5 ]; then
      say FAIL "budget_row_$b" "$# fields" "5 fields" 0
      continue
   fi
   loss="$2"; frac="$4"
   if awk -v x="$loss" 'BEGIN { exit !(x == 0) }'; then
      say PASS "unrated_${b}_no_flux" "beam_loss=0" "skipped" 0
   elif abs_lt "$frac" 1e-6; then
      say PASS "unrated_$b" "$frac" 0 1e-6
   else
      say FAIL "unrated_$b" "$frac" 0 1e-6
   fi
done

# ---- (b) the beam budget: the LW closure ---------------------------------
# The band the H2 lines pump and the band the flux is normalized over are
# one interval, so the photons the rate spends and the photons the beam
# loses are one number.  Before the merge this read -0.362 and changed sign
# along the column.
row=$(band_row '^# band rated_ph beam_loss_ph' LW)
set -- $row
if [ $# -ne 5 ]; then
   say FAIL budget_row_LW "$# fields" "5 fields" 0
else
   loss="$2"; frac="$4"
   if awk -v x="$loss" 'BEGIN { exit !(x == 0) }'; then
      say PASS unrated_LW_no_flux "beam_loss=0" "skipped" 0
   elif abs_lt "$frac" 5e-2; then
      say PASS unrated_LW "$frac" 0 5e-2
   else
      say FAIL unrated_LW "$frac" 0 5e-2
   fi
fi

# ---- (b) the beam budget: a SHARE where the header says share -----------
# Band B2 is the H I Ly-alpha line.  The case carries a stellar Ly-alpha
# flux and a molecular base, so the line is thick and H I takes nearly the
# whole band; what is asserted is that the number is a share of the beam,
# strictly inside (0, 1], and not a residual that has been reported as one.
row=$(band_row '^# band rated_ph beam_loss_ph' B2)
set -- $row
if [ $# -ne 5 ]; then
   say FAIL budget_row_B2 "$# fields" "5 fields" 0
else
   frac="$4"
   if gt "$frac" 0.5 && awk -v v="$frac" 'BEGIN { exit !(v <= 1.0) }'; then
      say PASS unrated_share_B2 "$frac" "in (0.5, 1]" 0
   else
      say FAIL unrated_share_B2 "$frac" "in (0.5, 1]" 0
   fi
fi

# ---- (c) the header states what each column is --------------------------
need() {  # need <name> <fixed string>
   if grep -qF "$2" "$LED"; then
      say PASS "$1" present present 0
   else
      say FAIL "$1" absent present 0
   fi
}
need header_names_closure   'it is NOT the test of the cell mean'
need header_names_cellmean  'only because the rate is the mean over the cell'
need header_names_share     'this is the SHARE H I takes'
need header_names_lw_closure 'CLOSURE residual like B3 and'
need header_names_state     'state_drift'

# ---- (d) CO, the fourth absorber of the shared band ---------------------
# CO predissociates in 37 lines between 912.7 and 1076.1 A, all inside the
# Lyman-Werner interval, and it shields itself in them.  It contributes no
# term to tau_cont -- its equivalent width is inside the Visser shielding
# function -- so block (a) stays a closure of the H2O and OH continua and
# it appears only in the rated total of block (b) and in three per-cell
# columns of the file.  The two checks below are the CO column list and
# the domain record of the destruction model the rate belongs to.
need header_names_co_column  'N_CO[cm^-2]'
need header_names_co_rate    'k_CO[1/s] Theta_CO'

OXY="$WORK/run/output/Oxygen_chemistry.txt"
for f in dom_f_dom dom_cells_out dom_worst_ratio dom_form_ratio \
         dom_cells_hot dom_cells_HeP; do
   if grep -q "^# $f " "$OXY" 2>/dev/null; then
      say PASS "co_domain_record_$f" present present 0
   else
      say FAIL "co_domain_record_$f" absent present 0
   fi
done

# Every cell must carry a shielding function in (0, 1]: it is a suppression
# factor of the unshielded rate, and a value above 1 would mean the table
# had been extrapolated rather than clamped.
bad=$(awk '/^#/ { next } NF > 3 { t = $(NF-1); if (t <= 0 || t > 1) n++ }
           END { print n+0 }' "$LED")
row=$(grep -m1 '^# co_absorbed_ph ' "$LED")
if [ -n "$row" ]; then
   say PASS co_ledger_row present present 0
else
   say FAIL co_ledger_row absent present 0
fi
if [ "$bad" = "0" ]; then
   say PASS co_shielding_in_unit_interval 0 0 0
else
   say FAIL co_shielding_in_unit_interval "$bad" 0 0
fi

echo ""
if [ $nfail -gt 0 ]; then
   echo "fuv_band_ledger: $nfail assertion(s) failed"
   exit 1
fi
echo "fuv_band_ledger: every assertion passed"

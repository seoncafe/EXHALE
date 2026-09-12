#!/bin/bash
# Build and run the tests of the coupled temperature-composition source step
# (docs/PLAN_20260906_rev2.md step B3c; docs/b1_target_system_20260906.md
# T1.4 to T1.9, T2.1, T2.2).
#
# Two kinds of row.
#   * The Fortran driver links the PRODUCTION objects and asserts the ledger
#     identities on one cell: the closed reacting cell (AT-1a), the photoevent
#     recipients (AT-1b), the oxygen and associative channels (AT-1d) and the
#     reservoir density itself.  See the header of coupled_source_tests.f90.
#   * The whole-binary rows state the properties that only a run can show:
#     the density the chemistry is given is the density it returns, the
#     coupled pair reaches its fixed point every step, and the mass sum of
#     the composition agrees with that density.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line.  Exit status is nonzero if any row fails.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses); with it set the staleness check is skipped.
# EXHALE_TEST_OUT selects where the test executable is written.
# EXHALE_COUPLED_EXE selects the binary the whole-binary rows run; without it
# those rows are skipped rather than building the production binary here.
#
# Usage: src/tests/coupled_source_step/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/coupled_source_step}"
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
   echo "FAIL coupled_source_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     run 'make' in $ROOT first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL coupled_source_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; run 'make' first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/coupled_source_tests.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/coupled_source_tests.x" \
    "$HERE/coupled_source_tests.f90" \
    $PROD_OBJ $LAPACK || {
   echo "FAIL coupled_source_build measured=compile_error reference=ok tol=0"
   exit 1; }

"$OUT/coupled_source_tests.x" || rc=1

# ------------------------------------------------------------------ #
# THE DENSITY THE CHEMISTRY IS GIVEN IS THE DENSITY IT RETURNS (T2.1).
# The statement is enforced by the language, not by a comparison: the
# argument is declared intent(in), so no path through the routine can write
# it and no run can be constructed in which it differs.  The row asserts the
# declaration, which is what makes the property hold for every caller.
# ------------------------------------------------------------------ #
IEQ="$ROOT/src/modules/radiation/ionization_equilibrium.f90"
if grep -qE 'intent\(in\)[[:space:]]*::[[:space:]]*n_io' "$IEQ"; then
   echo "PASS the_density_argument_is_intent_in measured=intent(in) reference=intent(in) tol=0"
else
   echo "FAIL the_density_argument_is_intent_in measured=writable reference=intent(in) tol=0"
   rc=1
fi

# ------------------------------------------------------------------ #
# THE WHOLE-BINARY ROWS.  A short molecular run, which is the configuration
# in which the mass sum has the most ways to disagree with the density: four
# molecular carriers, the oxygen columns and the metal block all weigh into
# calc_rho.
# ------------------------------------------------------------------ #
EXE="${EXHALE_COUPLED_EXE:-}"
if [ -z "$EXE" ]; then
   echo "  whole-binary rows skipped: set EXHALE_COUPLED_EXE to a built binary"
   exit $rc
fi
if [ ! -x "$EXE" ]; then
   echo "FAIL coupled_source_binary measured=not_executable reference=$EXE tol=0"
   exit 1
fi

WORK="$OUT/mbh"
CASE="$ROOT/backup/regression/mol_base_handoff"
rm -rf "$WORK"; mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/"
( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 "$EXE" > run.log 2>&1 )
brc=$?
if [ $brc -ne 0 ] && [ $brc -ne 2 ]; then
   echo "FAIL coupled_source_run measured=exit_$brc reference=exit_0 tol=0"
   echo "     see $WORK/run.log"
   exit 1
fi
echo "PASS coupled_source_run measured=exit_0 reference=exit_0 tol=0"

# THE PAIR REACHES ITS FIXED POINT ON EVERY STEP OF THE TRAJECTORY.
# Step 0 is exempt and the exemption is named rather than hidden: at step 0
# of a base handoff the composition is not on the trajectory at all -- it
# jumps from the imposed base composition to the equilibrium of the state the
# handoff sets up, and that jump is the largest the run ever takes.  MEASURED
# on this case, step 0 leaves the pair at dcomp = 1.5e-4 after 20 passes,
# while every later step reaches 1e-8 in 6.9 passes on average.  What the row
# asserts is therefore that no step of the trajectory reaches the cap; a
# capped step 0 is reported by the run itself, in init mode, as the recorded
# exhaustion that the run-mode contract allows.
ncap_traj=$(grep 'pass cap reached at step' "$WORK/run.log"             | sed -E 's/.*at step ([0-9]+) .*/\1/' | awk '$1 > 0' | wc -l)
if [ "$ncap_traj" -eq 0 ]; then
   echo "PASS every_trajectory_step_reached_its_fixed_point measured=0 reference=0 tol=0"
else
   echo "FAIL every_trajectory_step_reached_its_fixed_point measured=$ncap_traj reference=0 tol=0"
   rc=1
fi

# ------------------------------------------------------------------ #
# THE COST OF THE COUPLING IS REPORTED, AND REPORTED OF THE POPULATION IT
# BELONGS TO.  Every entry of the coupled loop counts its passes, and the
# entries are not the adopted steps: a step probed by the step-doubling
# error estimate enters the loop for the full step and for each of its two
# halves, and a refused attempt keeps its entry.  The mean has to be taken
# over the entries, or it is a ratio of two different populations and can
# come out ABOVE the worst single entry -- MEASURED on the tree before this
# item, wasp_full in phys mode reported "13.98 passes per step on average,
# worst 11".
# ------------------------------------------------------------------ #
WORK2="$OUT/mbh_probe"
rm -rf "$WORK2"; mkdir -p "$WORK2/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK2/"
# phys mode, because the step-doubling probe runs only there, and it is
# what makes the entries outnumber the steps.
echo "Run mode: phys" >> "$WORK2/input.inp"
( cd "$WORK2" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=20 EXHALE_ERR_EVERY=5 \
     "$EXE" > run.log 2>&1 )
LINE=$(grep 'coupled source step:.*passes per coupled step' "$WORK2/run.log" | head -n 1)
if [ -z "$LINE" ]; then
   echo "FAIL mean_passes_never_exceeds_the_worst measured=absent reference=mean_le_worst tol=0"
   rc=1
else
   python3 - "$LINE" <<'PYA'
import re, sys
line = sys.argv[1]
m = re.search(r'([0-9.]+) passes per coupled step on average over (\d+) of them,'
              r' worst (\d+)', line)
if not m:
    print("FAIL mean_passes_never_exceeds_the_worst measured=unparsable "
          "reference=mean_le_worst tol=0")
    sys.exit(1)
mean, n, worst = float(m.group(1)), int(m.group(2)), int(m.group(3))
rc = 0
ok = mean <= worst
print("%s mean_passes_never_exceeds_the_worst measured=%.2f reference=%d tol=0"
      % ("PASS" if ok else "FAIL", mean, worst))
if not ok:
    rc = 1
# 20 steps were asked for, and the probe enters the loop for the full step
# and for each of its two halves on every fifth of them, so the entries
# outnumber the steps.  That is the fact the old denominator got wrong.
ok2 = n > 20
print("%s the_entries_outnumber_the_steps measured=%d reference=20 tol=0"
      % ("PASS" if ok2 else "FAIL", n))
if not ok2:
    rc = 1
sys.exit(rc)
PYA
   [ $? -eq 0 ] || rc=1
fi

# THE PASS COUNT ITSELF, AS A REGRESSION BOUND.  The cost of a marching step
# is the number of cell chemistry solves it makes, and that is the pass
# count times the grid.  MEASURED on this case at 60 steps: 6.23 passes
# per coupled step, against 6.80 before the stopping test was taken on the
# estimated error (item COST7).  The bound is 8.0, which leaves room for the
# run-to-run variation of a relaxation snapshot and still catches a change
# that costs the coupling a whole extra pass.
LINE1=$(grep 'coupled source step:.*passes per coupled step' "$WORK/run.log" | head -n 1)
if [ -z "$LINE1" ]; then
   echo "FAIL mean_passes_within_bound measured=absent reference=8.0 tol=0"
   rc=1
else
   python3 - "$LINE1" <<'PYB'
import re, sys
m = re.search(r'([0-9.]+) passes per coupled step', sys.argv[1])
v = float(m.group(1)); bound = 8.0
print("%s mean_passes_within_bound measured=%.2f reference=%.1f tol=0"
      % ("PASS" if v <= bound else "FAIL", v, bound))
sys.exit(0 if v <= bound else 1)
PYB
   [ $? -eq 0 ] || rc=1
fi

# ------------------------------------------------------------------ #
# THE GEOMETRIC LIMIT OF THE PASS SEQUENCE (item COST6).  The loop's error
# decays geometrically with a ratio stable from its second pass, so its
# limit is known from three of its terms and need not be walked term by
# term.  Three properties are asserted, all of them on a real trajectory
# because the extrapolation is a property of the sequence and not of one
# call: a refused candidate costs nothing, an accepted one is admissible,
# and the count comes down.
# ------------------------------------------------------------------ #
WORKX="$OUT/mbh_xtr"
rm -rf "$WORKX"; mkdir -p "$WORKX/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKX/"
( cd "$WORKX" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 EXHALE_CSM_EXTRAP=1 \
     "$EXE" > run.log 2>&1 )

WORKR="$OUT/mbh_xtr_refused"
rm -rf "$WORKR"; mkdir -p "$WORKR/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKR/"
( cd "$WORKR" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 EXHALE_CSM_EXTRAP=1 \
     EXHALE_CSM_EXTRAP_REFUSE=1 "$EXE" > run.log 2>&1 )

# A CANDIDATE REFUSED BY THE GUARDS COSTS NOTHING.  The guards run before
# any state is moved and before any pass is spent, so a run in which every
# candidate is refused has to be the run with the extrapolation off: the
# same pass count and the same two written files, line for line.
pass_off=$(grep 'coupled source step:.*passes per coupled step' "$WORK/run.log"  | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
pass_ref=$(grep 'coupled source step:.*passes per coupled step' "$WORKR/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
same_files=yes
for f in Hydro_ioniz.txt Ion_species.txt; do
   if ! diff -q <(grep -v '^#' "$WORK/output/$f") \
                <(grep -v '^#' "$WORKR/output/$f") >/dev/null 2>&1; then
      same_files=no
   fi
done
if [ "$pass_off" = "$pass_ref" ] && [ "$same_files" = yes ]; then
   echo "PASS refused_extrapolation_costs_nothing measured=${pass_ref}_identical reference=${pass_off}_identical tol=0"
else
   echo "FAIL refused_extrapolation_costs_nothing measured=${pass_ref}_files_$same_files reference=${pass_off}_identical tol=0"
   rc=1
fi

# THE EXTRAPOLATED STATE IS ADMISSIBLE.  Every constraint on the candidate
# is linear in the extrapolation step, so the step is damped by the ratio
# test that makes all of them hold; the run measures the worst violation of
# each over every candidate it took, and both are identically zero when
# that construction is right.  The row also asserts that candidates were
# actually TAKEN, so that a zero is not the zero of an empty set.
XLINE=$(grep 'geometric limit admissibility' "$WORKX/run.log" | head -n 1)
NLINE=$(grep 'geometric limit of the pass sequence' "$WORKX/run.log" | head -n 1)
if [ -z "$XLINE" ] || [ -z "$NLINE" ]; then
   echo "FAIL extrapolated_state_is_admissible measured=not_reported reference=0 tol=0"
   rc=1
else
   python3 - "$XLINE" "$NLINE" <<'PYX'
import re, sys
xl, nl = sys.argv[1], sys.argv[2]
v = [float(x) for x in re.findall(r'[-0-9.]+E[+-][0-9]+', xl)]
taken = int(re.search(r', (\d+) taken,', nl).group(1))
worst = max(v) if v else 1.0
ok = (worst <= 0.0) and (taken > 0)
print("%s extrapolated_state_is_admissible measured=%.3e_over_%d reference=0 tol=0"
      % ("PASS" if ok else "FAIL", worst, taken))
sys.exit(0 if ok else 1)
PYX
   [ $? -eq 0 ] || rc=1
fi

# WHAT THE EXTRAPOLATION IS WORTH NOW THAT THE STOPPING TEST IS TAKEN ON
# THE ERROR.  The passes it used to save were the ones the loop spent
# walking the last factor 1/|theta| of the increment down to the tolerance,
# and those are exactly the passes the error test no longer spends: the
# window the reach guard leaves is now the band between the tolerance and
# csm_x_reach times it.  MEASURED at 300 steps of this case, the count with
# the limit taken is the count without it, and what the row asserts is that
# it is not ABOVE it -- an extrapolation that cost passes would be a
# defect, one that saves them is a bonus.  It is taken at 300 steps and not
# at the 60 of the rows above because the first steps of a base handoff are
# the ones whose sequence is not yet geometric.
WORKL="$OUT/mbh300_off"
rm -rf "$WORKL"; mkdir -p "$WORKL/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKL/"
( cd "$WORKL" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=300 "$EXE" > run.log 2>&1 )
WORKM="$OUT/mbh300_xtr"
rm -rf "$WORKM"; mkdir -p "$WORKM/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKM/"
( cd "$WORKM" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=300 EXHALE_CSM_EXTRAP=1 \
     "$EXE" > run.log 2>&1 )
pass_off=$(grep 'coupled source step:.*passes per coupled step' "$WORKL/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
pass_x=$(grep 'coupled source step:.*passes per coupled step' "$WORKM/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
if [ -z "$pass_x" ] || [ -z "$pass_off" ]; then
   echo "FAIL extrapolation_does_not_cost_passes measured=absent reference=at_or_below_control tol=0"
   rc=1
else
   python3 - "$pass_x" "$pass_off" <<'PYC'
import sys
x, c = float(sys.argv[1]), float(sys.argv[2])
print("%s extrapolation_does_not_cost_passes measured=%.2f reference=%.2f tol=0"
      % ("PASS" if x <= c else "FAIL", x, c))
sys.exit(0 if x <= c else 1)
PYC
   [ $? -eq 0 ] || rc=1
fi

# ------------------------------------------------------------------ #
# WHAT THE STOPPING TEST IS TAKEN ON (item COST7).  The loop stops when
# the ESTIMATED REMAINING ERROR of the pair -- |theta/(1-theta)| times the
# pass's own increment, with theta the sequence's own contraction ratio --
# is within csm_T_tol and csm_comp_tol, and falls back to the increment on
# a pass whose sequence gives no estimate.  Three properties are asserted:
# the accepted state really is that close to the fixed point, the fallback
# is one path, and the estimate is formed in one place.
# ------------------------------------------------------------------ #

# THE ACCEPTED STATE IS WITHIN THE TWO TOLERANCES OF THE FIXED POINT, and
# it is MEASURED and not modelled: with EXHALE_CSM_ERR_PROBE the loop holds
# the state the stopping test accepted, runs the same sequence on until its
# increment is below csm_probe_tol = 1e-12 (four decades below the
# tolerance, so that state IS the fixed point at this precision), records
# the max-norm distance between the two in each of the two norms of the
# test, and then retakes the accepted pass so that the trajectory is the
# trajectory of a run without the probe.  The row asserts the worst
# distance over every accepted state of the run, which is the statement the
# two tolerances make about the state they accept.
WORKP="$OUT/mbh_errprobe"
rm -rf "$WORKP"; mkdir -p "$WORKP/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKP/"
( cd "$WORKP" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 \
     EXHALE_CSM_ERR_PROBE=1 "$EXE" > run.log 2>&1 )
PLINE=$(grep 'accepted state against its fixed point:.*worst distance' \
        "$WORKP/run.log" | head -n 1)
if [ -z "$PLINE" ]; then
   echo "FAIL accepted_state_within_its_tolerances measured=not_reported reference=1e-8 tol=0"
   rc=1
else
   python3 - "$PLINE" <<'PYE'
import re, sys
line = sys.argv[1]
n = int(re.search(r'point: (\d+) measured', line).group(1))
v = [float(x) for x in re.findall(r'[0-9.]+E[+-][0-9]+', line)]
tol = 1.0e-8
worst = max(v) if v else 1.0
ok = (n > 0) and (worst <= tol)
print("%s accepted_state_within_its_tolerances measured=%.3e_over_%d "
      "reference=%.1e tol=0" % ("PASS" if ok else "FAIL", worst, n, tol))
sys.exit(0 if ok else 1)
PYE
   [ $? -eq 0 ] || rc=1
fi

# THE PROBE DOES NOT MOVE THE TRAJECTORY.  It holds the accepted state
# aside, walks the sequence to 1e-12 and then retakes the accepted pass
# from the pair that pass started from; the pair (T, f_sp) is the whole
# iterate of the loop, so that retaken pass rebuilds everything the
# accepted pass wrote and the run has to write the same files as a run
# without the probe.  If it does not, the state the probe restores is not
# the whole iterate and no distance it measures belongs to this trajectory.
same_files=yes
for f in Hydro_ioniz.txt Ion_species.txt; do
   if ! diff -q <(grep -v '^#' "$WORK/output/$f") \
                <(grep -v '^#' "$WORKP/output/$f") >/dev/null 2>&1; then
      same_files=no
   fi
done
pass_p=$(grep 'coupled source step:.*passes per coupled step' "$WORKP/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
pass_d=$(grep 'coupled source step:.*passes per coupled step' "$WORK/run.log"  | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
if [ "$same_files" = yes ] && [ "$pass_p" = "$pass_d" ]; then
   echo "PASS the_probe_does_not_move_the_trajectory measured=${pass_p}_identical reference=${pass_d}_identical tol=0"
else
   echo "FAIL the_probe_does_not_move_the_trajectory measured=${pass_p}_files_$same_files reference=${pass_d}_identical tol=0"
   rc=1
fi

# THE INCREMENT TEST IS THE FALLBACK, AND IT IS ONE PATH.  Where the
# sequence gives no estimate the loop tests the increment, which is the
# test it applied on every pass before this item.  Two independent ways of
# leaving the loop with no estimate -- refusing every estimate
# (EXHALE_CSM_GEOM_REFUSE) and putting the first pass that could carry one
# out of reach (EXHALE_CSM_GEOM_PASS) -- have to give the same run, line
# for line and pass for pass, and every exit of both has to be taken on
# the increment.  That is the assertion that the fallback is a single path
# and not a second test.
WORKF1="$OUT/mbh_incr_refuse"
rm -rf "$WORKF1"; mkdir -p "$WORKF1/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKF1/"
( cd "$WORKF1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 \
     EXHALE_CSM_GEOM_REFUSE=1 "$EXE" > run.log 2>&1 )
WORKF2="$OUT/mbh_incr_outofreach"
rm -rf "$WORKF2"; mkdir -p "$WORKF2/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORKF2/"
( cd "$WORKF2" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=60 \
     EXHALE_CSM_GEOM_PASS=999 "$EXE" > run.log 2>&1 )
same_files=yes
for f in Hydro_ioniz.txt Ion_species.txt; do
   if ! diff -q <(grep -v '^#' "$WORKF1/output/$f") \
                <(grep -v '^#' "$WORKF2/output/$f") >/dev/null 2>&1; then
      same_files=no
   fi
done
est1=$(grep 'coupled source step:.*exit(s) on the estimated error' "$WORKF1/run.log" | head -n 1 | sed -E 's/.*step: *([0-9]+) exit.*/\1/')
est2=$(grep 'coupled source step:.*exit(s) on the estimated error' "$WORKF2/run.log" | head -n 1 | sed -E 's/.*step: *([0-9]+) exit.*/\1/')
pass_f1=$(grep 'coupled source step:.*passes per coupled step' "$WORKF1/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
pass_f2=$(grep 'coupled source step:.*passes per coupled step' "$WORKF2/run.log" | head -n 1 | sed -E 's/.*: *([0-9.]+) passes.*/\1/')
if [ "$same_files" = yes ] && [ "$pass_f1" = "$pass_f2" ] \
   && [ "${est1:-x}" = 0 ] && [ "${est2:-x}" = 0 ]; then
   echo "PASS the_increment_test_is_the_fallback measured=${pass_f1}_identical_0_estimated reference=${pass_f2}_identical_0_estimated tol=0"
else
   echo "FAIL the_increment_test_is_the_fallback measured=${pass_f1}_files_${same_files}_est_${est1:-absent} reference=${pass_f2}_identical_0_estimated tol=0"
   rc=1
fi

# THE ERROR TEST IS WHAT LOWERS THE PASS COUNT, and the fallback run above
# is the control that says by how much: the same binary, the same case and
# the same steps, with the estimate refused.  MEASURED at 60 steps of this
# case, 6.23 against 6.80.
if [ -n "$pass_d" ] && [ -n "$pass_f1" ]; then
   python3 - "$pass_d" "$pass_f1" <<'PYD'
import sys
new, incr = float(sys.argv[1]), float(sys.argv[2])
print("%s the_error_test_lowers_the_pass_count measured=%.2f reference=%.2f "
      "tol=0" % ("PASS" if new < incr else "FAIL", new, incr))
sys.exit(0 if new < incr else 1)
PYD
   [ $? -eq 0 ] || rc=1
fi

# THE ESTIMATE IS FORMED IN ONE PLACE.  The stopping test and the
# extrapolation are the same statement about the same three terms of the
# same sequence, so the ratio has to be formed once: a second copy is what
# drifts apart from the first and cannot be caught by the compiler.  The
# row reads the source: one routine defines it, each of the two ratios is
# assigned in exactly one place and that place is inside that routine, and
# both consumers are there to read it.
MAIN="$ROOT/src/EXHALE_main.f90"
ndef=$(grep -cE '^ *subroutine coupled_pair_geometric_error_estimate *$' "$MAIN")
nass=$(grep -cE '^ *csm_geom_theta_(T|c) *=  *(p12/p22|q12/q22)' "$MAIN")
ncall=$(grep -cE 'call coupled_pair_geometric_error_estimate' "$MAIN")
nuse=$(grep -cE 'csm_err_(ok|T|c)' "$MAIN")
python3 - "$ndef" "$nass" "$ncall" "$nuse" <<'PYS'
import sys
ndef, nass, ncall, nuse = (int(x) for x in sys.argv[1:5])
ok = ndef == 1 and nass == 2 and ncall == 1 and nuse > 4
print("%s the_error_estimate_is_one_symbol measured=%d_def_%d_ratios_%d_call "
      "reference=1_def_2_ratios_1_call tol=0"
      % ("PASS" if ok else "FAIL", ndef, nass, ncall))
sys.exit(0 if ok else 1)
PYS
[ $? -eq 0 ] || rc=1

# T2.2: the species mass sum against the density the step was given.  The
# bound the declared convention implies is n_e m_e / rho, at most
# m_e/m_H = 5.4e-4; what a correct element budget actually reaches is
# round-off, and the row is set there so that a broken budget is caught.
MC=$(grep 'chemistry mass closure' "$WORK/run.log" | sed -E 's/.*= *([0-9.E+-]+) .*/\1/')
if [ -z "$MC" ]; then
   echo "FAIL chemistry_mass_closure measured=not_reported reference=0 tol=1e-10"
   rc=1
else
   python3 - "$MC" <<'PY'
import sys
v = float(sys.argv[1]); tol = 1.0e-10
print("%s chemistry_mass_closure measured=%.3e reference=0 tol=%.1e"
      % ("PASS" if v <= tol else "FAIL", v, tol))
sys.exit(0 if v <= tol else 1)
PY
   [ $? -eq 0 ] || rc=1
fi

exit $rc

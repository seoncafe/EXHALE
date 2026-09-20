#!/bin/bash
# ONE STATE PER OUTPUT FILE.  Two statements, each with its own run:
#   0. the cool column of Hydro_ioniz.txt is the cooling of the state written
#      beside it, recomputed channel by channel from (T, rho, f_sp) into
#      Cooling_breakdown.txt (item L20);
#   1. the heat column of Hydro_ioniz.txt is the heat the channel breakdown
#      of the same run adds up to (the '# coupling:' header names the physics
#      that produced it);
#   2. the state an ENDING of the stationary outer iteration hands back is one
#      state: the composition of the file, the particle densities and the
#      temperature beside it belong together, so re-evaluating the state as
#      written reproduces it to round-off.
#
# QUANTITY UNDER TEST 1
#   Hydro_ioniz.txt column 6, `heat` [erg cm^-3 s^-1], against
#   Heating_breakdown.txt column 4, `heat_total`, cell by cell.  The
#   breakdown file's own header states the relation: "Channel sum reproduces
#   the heat_total column (and the Hydro_ioniz.txt heat column up to
#   convergence)", so on one state the two columns are the same number.
#
#   The two columns come from ONE assembly of the heating,
#   heating_of_composition (src/modules/radiation/util_ion_eq.f90 1503),
#   which fills a channel array and returns the running sum of its columns.
#   The ionization sweep calls it on the composition it returns
#   (src/modules/radiation/ionization_equilibrium.f90 2192 and following),
#   the `heat` column is that total, and write_heat_breakdown_eq
#   (src/modules/radiation/util_ion_eq.f90 2613) writes the very array the
#   sweep filled.  Nothing is recomputed for the dump, so what separates the
#   two columns is the decimal round trip of the two files and the order in
#   which seventeen channels are added, not the size of one sweep's rate
#   lag.
#
#   Rebuilding the rates on the written state instead, which the dump used to
#   do, put a rate lag one sweep wide between the two columns (measured at
#   1.24e-5 in the median cell of this case) and left out any deposit the
#   dump's own copy of the sum did not carry (measured at 6.89e-3 in the
#   worst cell, the associative He(2^3S) branch and the collisional oxygen
#   channels).
#
#   The secondary-ionization coupling is not the cause any more.  The final
#   write leaves sec_ion_active exactly as the marching loop left it
#   (src/EXHALE_main.f90 1826-1848); the only activations are the staged one
#   and the pre-Newton one, both inside the loop (:1649, :1698).  A run whose
#   staged activation never fired now writes sec_ion=F beside a heat column
#   produced without the coupling, so write_coupling_state_header
#   (src/modules/functions/utilities.f90 141) no longer labels one state with
#   another state's physics.
#
# WHAT IS RUN
#   A copy of backup/regression/mol_base_handoff (its input.inp and base.inp)
#   in build/tests/grid_and_gates/mbh, capped at 200 steps.  The shipped case
#   is a 12000-step snapshot; 200 steps is enough to write both files, and
#   the difference is present at any step count because the rate lag is a
#   property of one sweep and not of how far the run has relaxed.  The
#   regression directory is never written to.
#
# REFERENCE AND TOLERANCE
#   reference = 0: the largest relative difference over the physical cells,
#   |heat - heat_total| / max(|heat|,|heat_total|), must be at most 1e-6.
#   Two names for one number on one state agree to round-off; 1e-6 leaves
#   room for the decimal round trip of the two files.
#
# MEASURED on the one assembly (2026-09-06): 1.68e-16, with the header
# reading sec_ion=F.  MEASURED on the two hand-maintained copies it
# replaced: 6.89e-3.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" output_state_consistency
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
WORK="$OUT/mbh"
CASE="$ROOT/backup/regression/mol_base_handoff"
STAG="$OUT/outer_iteration_ending"
CASE2="$ROOT/backup/regression/carrier_elem_newton"

if [ ! -x "$EXE" ]; then
   echo "FAIL output_state_consistency measured=no_binary reference=$EXE tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/"

( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=200 "$EXE" > run.log 2>&1 )
rc=$?
# Exit status 2 = the run declared a stationary state that the A2 certification
# refused (docs/a2_certification_contract_20260906.md section 8); the run wrote
# its outputs in full, which is what this gate reads. Only 1 (a Fortran error
# stop) or a signal is a failed run here.
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL output_state_consistency_run measured=exit_$rc reference=exit_0 tol=0"
   echo "     see $WORK/run.log"
   exit 1
fi
for f in output/Hydro_ioniz.txt output/Heating_breakdown.txt \
         output/Cooling_breakdown.txt; do
   if [ ! -s "$WORK/$f" ]; then
      echo "FAIL output_state_consistency_files measured=missing_$f reference=written tol=0"
      exit 1
   fi
done
echo "PASS output_state_consistency_run measured=exit_0 reference=exit_0 tol=0"

python3 - "$WORK" <<'PY'
import re, sys
import numpy as np

work = sys.argv[1]
hydro = work + "/output/Hydro_ioniz.txt"
brk = work + "/output/Heating_breakdown.txt"

coupling = ""
n_cells = None
n_ghost = None
with open(hydro) as fh:
    for line in fh:
        if not line.startswith("#"):
            break
        if line.startswith("# coupling:"):
            coupling = line.strip()
        m = re.search(r"rows\s+\d+:\s+(\d+)\s+ghost cells", line)
        if m:
            n_ghost = int(m.group(1))
        m = re.search(r"\bN=(\d+)", line)
        if m:
            n_cells = int(m.group(1))
print("  header written by the run: %s" % coupling)

h = np.loadtxt(hydro, comments="#")
b = np.loadtxt(brk, comments="#")
if h.shape[0] != b.shape[0]:
    print("FAIL output_state_consistency_rows measured=%d reference=%d tol=0"
          % (h.shape[0], b.shape[0]))
    sys.exit(1)

lo, hi = n_ghost, n_ghost + n_cells      # physical cells only
heat = h[lo:hi, 5]
heat_total = b[lo:hi, 3]
scale = np.maximum(np.abs(heat), np.abs(heat_total))
good = scale > 0.0
rel = np.zeros_like(scale)
rel[good] = np.abs(heat[good] - heat_total[good]) / scale[good]

ratio = np.full_like(scale, np.nan)
nz = heat != 0.0
ratio[nz] = heat_total[nz] / heat[nz]
print("  heat_total/heat over the %d physical cells: median %.6f  min %.6f"
      "  max %.6f" % (hi - lo, np.nanmedian(ratio), np.nanmin(ratio),
                      np.nanmax(ratio)))

worst = float(np.max(rel))
tol = 1.0e-6
ok = worst <= tol
print("%s heat_column_matches_breakdown_total measured=%.6e reference=%.6e "
      "tol=%.2e" % ("PASS" if ok else "FAIL", worst, 0.0, tol))

# THE SAME STATEMENT FOR THE COOLING, which nothing asserted until item L20.
#
# The two columns are not guarded the same way by construction, and that is
# why this row exists.  Heating_breakdown.txt is WRITTEN FROM the channel
# array the sweep filled (heat_channel_state), so the heat column and that
# file cannot disagree by more than the copy; Cooling_breakdown.txt RECOMPUTES
# every channel from the written (T, rho, f_sp) through eval_cool, so the
# comparison below is a real statement about the state: it says that the cool
# column of the file is the cooling of the state written beside it.
#
# It caught a stale fixture the first time it was applied: the pinned pair
# backup/regression/wasp_full_newton/IC/, written on 2026-09-08 by an earlier
# code generation, carries a cool column standing at 1.3 to 2.9 times the
# cooling of its own state while its heat column reproduces exactly
# (docs/lhs1140b_stationary_L19_L20_20260915.md).  Every state the current
# code writes agrees to between 3.4e-15 and 3.0e-10, the documented one-sweep
# lag of the temperature; the tolerance is the same 1e-6 the heating row uses.
cbrk = work + "/output/Cooling_breakdown.txt"
c = np.loadtxt(cbrk, comments="#")
if c.shape[0] != h.shape[0]:
    print("FAIL cool_column_rows measured=%d reference=%d tol=0"
          % (c.shape[0], h.shape[0]))
    sys.exit(1)
cool = h[lo:hi, 6]
cool_total = c[lo:hi, 3]
cscale = np.maximum(np.abs(cool), np.abs(cool_total))
cgood = cscale > 0.0
crel = np.zeros_like(cscale)
crel[cgood] = np.abs(cool[cgood] - cool_total[cgood]) / cscale[cgood]
cratio = np.full_like(cscale, np.nan)
cnz = cool != 0.0
cratio[cnz] = cool_total[cnz] / cool[cnz]
print("  cool_total/cool over the %d physical cells: median %.6f  min %.6f"
      "  max %.6f" % (hi - lo, np.nanmedian(cratio), np.nanmin(cratio),
                      np.nanmax(cratio)))
cworst = float(np.max(crel))
cok = cworst <= tol
print("%s cool_column_matches_breakdown_total measured=%.6e reference=%.6e "
      "tol=%.2e" % ("PASS" if cok else "FAIL", cworst, 0.0, tol))
sys.exit(0 if (ok and cok) else 1)
PY
rc1=$?

# ---------------------------------------------------------------------------
# QUANTITY UNDER TEST 2
#   The temperature column of the state the stationary outer iteration hands
#   back, against the temperature the same conserved state and the same
#   composition give when they are read back and evaluated, over the physical
#   cells.  T is p/(n_tot + n_e) of the composition beside it, so the two
#   agree on ONE state and part wherever the file's composition has moved past
#   the particle count and the temperature written with it.
#
#   Every ending of steady_wind_with_element_diffusion (src/EXHALE_main.f90)
#   has to hand back one state, because the caller certifies f_sp and writes
#   n_HI ... n_m and T: a composition updated after the certification would be
#   certified as one state and written as another.  The ending under test is
#   the stagnation ending, the one whose composition update stood between the
#   certification and the return.
#
# WHAT IS RUN (two runs; nothing in backup/regression is written to)
#   S  a copy of backup/regression/carrier_elem_newton, the hot-Uranus carrier
#      reload, partitioned (`Coupled carrier solve: False`) with
#      `Restart intent: stationary`, run until it stagnates: at
#      EXHALE_JFNK_MAXIT=40 the hydrodynamic solve reaches its root
#      (info = 0, the three hydrodynamic rows inside their tolerances) and
#      EXHALE_CARRIER_TRUST=1e-6 bounds the carrier pass so tightly that it
#      keeps no transport step at all, so the H2 row, the entry the joint
#      distance is then set by, stands at 7.397E-02 of 1.0E-05 from pass 2
#      to the end.
#
#      The ending asks for outer_no_fall_max = 3 CONSECUTIVE passes in which
#      neither the joint distance nor the composition residual falls, and the
#      budget has to outlast the one thing that is still moving: the
#      composition elimination sweep the handed-back state carries, which
#      moves the composition a little every pass and therefore moves the H2
#      row in a digit below the four the log prints.  The budget is the
#      setting that has moved with the tree, not the assertion: a budget of 8
#      reached the ending at pass 4 when this row was written (READ,
#      2026-09-12), 50 reached it at pass 44 (MEASURED 2026-09-18), and since
#      the joint test of the outer iteration derives the lower boundary from
#      the composition it is taken at (item D5b-2, boundary model v1) the
#      joint distance carries one more digit of movement per pass and the
#      ending comes later still: MEASURED 2026-09-18, single-threaded, pass 89
#      announces it, so the budget here is 100.  If a future change pushes it
#      past that, RAISE THE BUDGET or pick a state that stagnates sooner; do
#      not let the row accept `pass_budget`, which is the ending it exists to
#      tell apart.
#
#      Loosening the movement bound does not substitute for the budget:
#      at EXHALE_CARRIER_TRUST=1e-4 the carrier does take steps and the joint
#      distance creeps down monotonically, 7.40E+03 to 7.36E+03 over 12
#      passes, without one pass that fails to fall (MEASURED 2026-09-18,
#      with and without the JFNK cap).  The still earlier setting,
#      EXHALE_JFNK_MAXIT=5 with EXHALE_CARRIER_TRUST=1e-4, reached the ending
#      only while the progress measure was the worst species row alone: on
#      the joint measure the starved hydrodynamic rows (energy 4.9E-01
#      falling to 3.9E-02 over 8 passes) keep falling and that configuration
#      runs out its budget (MEASURED 2026-09-12).
#   R  the state S wrote, handed back as the _IC pair with
#      `Restart intent: stationary evaluate`, the route that measures a loaded
#      state and writes it back unchanged.  No step and no solve are taken, so
#      what separates R's columns from S's is what S's own state carries.
#
#   If a future change makes this configuration converge or run out its
#   budget instead, the first assertion says so by name: pick a setting that
#   reaches the ending again rather than loosening the second assertion.
#
# REFERENCE AND TOLERANCE
#   reference = 0: the largest relative difference of the T column over the
#   physical cells must be at most 1e-12.  MEASURED on the state this fixture
#   hands back at the stagnation ending: 6.4e-16 (2026-09-18), the decimal
#   round trip of the file's own digits, with the worst species column of the
#   same two files moving by 1.5e-14 (HI).  MEASURED where the ending stood
#   after the update instead: 1.3e-8, four decades above the allowance.
if [ ! -f "$CASE2/input.inp" ]; then
   echo "FAIL outer_iteration_ending_case measured=no_case reference=$CASE2 tol=0"
   exit 1
fi
rm -rf "$STAG"
mkdir -p "$STAG/solve/output" "$STAG/reeval/output"
cp "$CASE2/base.inp" "$STAG/solve/"
sed 's/^Coupled carrier solve:.*/Coupled carrier solve: False/' \
    "$CASE2/input.inp" > "$STAG/solve/input.inp"
printf 'Restart intent: stationary\n' >> "$STAG/solve/input.inp"
cp "$CASE2"/IC/*.txt "$STAG/solve/output/"
( cd "$STAG/solve" && env OMP_NUM_THREADS=1 EXHALE_CARRIER_TRUST=1e-6 \
  EXHALE_OUTER_PASSES=100 EXHALE_JFNK_MAXIT=40 "$EXE" > run.log 2>&1 )
rc=$?
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL outer_iteration_ending_run measured=exit_$rc reference=exit_0 tol=0"
   echo "     see $STAG/solve/run.log"
   exit 1
fi

# The ending this row is about, read from the line the routine prints for it.
if grep -q 'has not fallen in 3 consecutive passes' "$STAG/solve/run.log"; then
   ending=no_progress
elif grep -q 'spent its budget of' "$STAG/solve/run.log"; then
   ending=pass_budget
elif grep -q 'ACCEPTED -- every active equation' "$STAG/solve/run.log"; then
   ending=certified
else
   ending=unrecognized
fi
# WHAT THIS ROW IS FOR. It is not a statement about the solve: it names the
# ENDING the row below is about. The stagnation ending is the one whose
# composition update stood between the certification and the return, so only
# on that ending does the next row measure anything; reaching the pass budget
# or certifying instead means the fixture no longer exercises it, and the
# answer is a fixture that stagnates again, not a row that accepts another
# ending.
if [ "$ending" = no_progress ]; then
   echo "PASS outer_iteration_ending_is_the_stagnation_one measured=$ending"\
        "reference=no_progress tol=0"
else
   echo "FAIL outer_iteration_ending_is_the_stagnation_one measured=$ending"\
        "reference=no_progress tol=0"
   echo "     this configuration no longer reaches the ending under test;"
   echo "     see $STAG/solve/run.log"
   exit 1
fi

cp "$CASE2/base.inp" "$STAG/reeval/"
sed 's/^Restart intent:.*/Restart intent: stationary evaluate/' \
    "$STAG/solve/input.inp" > "$STAG/reeval/input.inp"
cp "$STAG/solve/output/Hydro_ioniz.txt" "$STAG/reeval/output/Hydro_ioniz_IC.txt"
cp "$STAG/solve/output/Ion_species.txt" "$STAG/reeval/output/Ion_species_IC.txt"
( cd "$STAG/reeval" && env OMP_NUM_THREADS=1 "$EXE" > run.log 2>&1 )
rc=$?
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL outer_iteration_ending_reeval measured=exit_$rc reference=exit_0 tol=0"
   echo "     see $STAG/reeval/run.log"
   exit 1
fi
grep 'the loaded composition against' "$STAG/reeval/run.log" | sed 's/^ */  /'

python3 - "$STAG" <<'PY'
import re, sys
import numpy as np

stag = sys.argv[1]


def read(path):
    labels, n_cells, n_ghost = None, None, None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.lstrip("#").split()[:1] == ["columns"]:
                labels = [c.split("[")[0]
                          for c in line.replace("#", "", 1).split()[1:]]
            m = re.search(r"rows\s+\d+:\s+(\d+)\s+ghost cells", line)
            if m:
                n_ghost = int(m.group(1))
            m = re.search(r"\bN=(\d+)", line)
            if m:
                n_cells = int(m.group(1))
    return labels, np.loadtxt(path, comments="#"), n_ghost, n_cells


la, A, ng, nc = read(stag + "/solve/output/Hydro_ioniz.txt")
lb, B, ngb, ncb = read(stag + "/reeval/output/Hydro_ioniz.txt")
if A.shape != B.shape or ng != ngb or nc != ncb:
    print("FAIL outer_iteration_ending_rows measured=%s reference=%s tol=0"
          % (str(A.shape), str(B.shape)))
    sys.exit(1)
lo, hi = ng, ng + nc                      # physical cells only
iT = la.index("T")
a, b = A[lo:hi, iT], B[lo:hi, iT]
rel = np.abs(a - b) / np.maximum(np.abs(a), 1.0e-99)
worst = float(np.max(rel))
jw = int(np.argmax(rel))
print("  worst T cell %d: written %.12e  re-evaluated %.12e" % (jw + 1, a[jw],
                                                                b[jw]))

# Context, not a verdict: the species columns of the same two files.  A trace
# stage answers to the temperature with a steep exponential, so its column is
# a sensitive indicator and not a reference of its own.
ls, SA, _, _ = read(stag + "/solve/output/Ion_species.txt")
_, SB, _, _ = read(stag + "/reeval/output/Ion_species.txt")
wsp, wname = 0.0, "none"
for k, nm in enumerate(ls):
    if nm == "r":
        continue
    x, y = SA[lo:hi, k], SB[lo:hi, k]
    sc = np.maximum(np.abs(x), np.abs(y))
    g = sc > 0.0
    if not g.any():
        continue
    d = float(np.max(np.abs(x[g] - y[g]) / sc[g]))
    if d > wsp:
        wsp, wname = d, nm
print("  worst species column over the same cells: %.3e (%s)" % (wsp, wname))

tol = 1.0e-12
ok = worst <= tol
print("%s outer_iteration_ending_hands_back_one_state measured=%.6e "
      "reference=%.6e tol=%.2e" % ("PASS" if ok else "FAIL", worst, 0.0, tol))
sys.exit(0 if ok else 1)
PY
rc2=$?

[ $rc1 -eq 0 ] && [ $rc2 -eq 0 ]
exit $?

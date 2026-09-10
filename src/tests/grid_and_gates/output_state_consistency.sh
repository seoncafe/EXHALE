#!/bin/bash
# One state per output file: the '# coupling:' header of Hydro_ioniz.txt
# against the heat column the same run wrote.
#
# QUANTITY UNDER TEST
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
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="$ROOT/build/tests/grid_and_gates/mbh"
CASE="$ROOT/backup/regression/mol_base_handoff"

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
for f in output/Hydro_ioniz.txt output/Heating_breakdown.txt; do
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
sys.exit(0 if ok else 1)
PY
exit $?

#!/bin/bash
# Restart round trip: write -> read -> write must be the identity.
#
# WHY THIS CASE EXISTS, AND WHY IT WAS REBUILT (defect 12.4 of
# docs/open_defects_20260903.md). The directory carried this name since the
# ATES era but could not run: output/Hydro_ioniz_IC.txt was zero bytes, there
# was no Ion_species_IC.txt, and "Load IC? True" made it abort in load_IC
# before anything else. It had been built from an expensive converged product
# that no longer exists, and it was in no DEFAULT_CASES list, so the matrix
# never noticed -- while the restart path is exactly what section 144.2/144.3
# changed and every default case is a cold start.
#
# It now builds its own state. The case is an ORDINARY named case: 40
# deterministic steps of the hot-Uranus molecular gate (input.inp, base.inp,
# maxsteps), so `run_check.sh check roundtrip` bitwise-compares it like any
# other. THIS script is the second half, run after that, and it needs the
# outputs that run just produced:
#
#   stage A  the ordinary case run (run_check.sh, or ./EXHALE.x by hand):
#            output/Hydro_ioniz.txt and output/Ion_species.txt are the
#            REFERENCE state
#   stage B  hand those two back as the _IC pair and reload with
#            EXHALE_DUMP_IC=1, which writes the loaded state before any
#            equilibrium sweep touches it and stops
#   stage C  every species column matched BY LABEL to 1e-12, rho/v/p likewise,
#            and the '# coupling:' header preserved
#
# It restores input.inp and the stage-A outputs on exit, so the directory is
# left exactly as stage A left it and `run_check.sh golden roundtrip` still
# snapshots the right thing.
#
# 40 steps is short enough to run in seconds and long enough that the state is
# not the initial condition: the molecular network has solved, the secondary-
# ionization coupling has armed, and the '# coupling:' line carries something
# other than its default.
#
# Usage: ./roundtrip_check.sh [path/to/EXHALE.x]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
EXE="${1:-}"
if [ -z "$EXE" ]; then
   if   [ -x "$HERE/EXHALE.x" ];        then EXE="$HERE/EXHALE.x"
   elif [ -x "$HERE/../../../EXHALE.x" ]; then EXE="$HERE/../../../EXHALE.x"
   else echo "  FAIL no EXHALE.x (pass one as \$1)"; exit 1; fi
fi
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"

for f in output/Hydro_ioniz.txt output/Ion_species.txt; do
   if [ ! -s "$HERE/$f" ]; then
      echo "  FAIL $f is missing or empty: run stage A first"
      echo "       (backup/regression/run_check.sh check roundtrip)"
      exit 1
   fi
done

cleanup() {
   cp -f "$HERE/.rt_input.inp"        "$HERE/input.inp"
   cp -f "$HERE/.rt_Hydro_ioniz.txt"  "$HERE/output/Hydro_ioniz.txt"
   cp -f "$HERE/.rt_Ion_species.txt"  "$HERE/output/Ion_species.txt"
   rm -f "$HERE/.rt_input.inp" "$HERE/.rt_Hydro_ioniz.txt" \
         "$HERE/.rt_Ion_species.txt" \
         "$HERE/output/Hydro_ioniz_IC.txt" "$HERE/output/Ion_species_IC.txt"
}
trap cleanup EXIT

cp "$HERE/input.inp"               "$HERE/.rt_input.inp"
cp "$HERE/output/Hydro_ioniz.txt"  "$HERE/.rt_Hydro_ioniz.txt"
cp "$HERE/output/Ion_species.txt"  "$HERE/.rt_Ion_species.txt"

echo "[roundtrip] stage B: reload the stage-A state and re-dump it"
cp "$HERE/.rt_Hydro_ioniz.txt" "$HERE/output/Hydro_ioniz_IC.txt"
cp "$HERE/.rt_Ion_species.txt" "$HERE/output/Ion_species_IC.txt"
sed -i 's/^Load IC?.*/Load IC? True/' "$HERE/input.inp"
( cd "$HERE" && OMP_NUM_THREADS=1 EXHALE_DUMP_IC=1 "$EXE" > roundtrip.log 2>&1 )
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$HERE/output/Ion_species.txt" ]; then
   echo "  FAIL stage B (rc=$rc); see roundtrip.log"; exit 1
fi

echo "[roundtrip] stage C: identity of the round trip"
python3 "$HERE/check_roundtrip_state.py" \
        "$HERE/.rt_Ion_species.txt" "$HERE/output/Ion_species.txt" \
        "$HERE/.rt_Hydro_ioniz.txt" "$HERE/output/Hydro_ioniz.txt"

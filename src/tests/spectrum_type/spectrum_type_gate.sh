#!/bin/bash
# The "Spectrum type" key states the spectrum of the WHOLE grid, and the run
# says which one it used.
#
# QUANTITY UNDER TEST
#   the acceptance of "Spectrum type: Planck" by input_read, the two stellar
#   quantities that type needs (Stellar Teff, Stellar radius), and the lines
#   write_setup_report writes about the spectrum for every run.  Decision 13
#   of docs/development_plan_20260905_rev3.md section 10.5: one spectrum type
#   builds every band, so the type and the source of the band below 13.6 eV
#   are part of the record of a run.
#
# WHAT IS RUN (three one-step runs on copies of backup/regression/wasp_full;
# the regression directory is not written to)
#   A  the case as shipped (Spectrum type: Power-law).  ASSERT exit 0, and
#      EXHALE_setup.out states the spectrum type, names the source of the
#      band below 13.6 eV, and prints the integrated grid flux.
#   B  Spectrum type: Planck with the Stellar Teff / Stellar radius lines of
#      the case.  ASSERT exit 0 and the setup report names the Planck field
#      as the source of every band.
#   C  Spectrum type: Planck with both stellar lines removed.  ASSERT exit
#      nonzero and a message naming both keys: the field cannot be built.
#
# TOLERANCE: exit status and the presence of the named strings, tol=0.
#
# EXPECTED BEFORE THE CHANGE: RED (B and C both stop with "unknown spectrum
# type", and A's setup report has no sub-13.6 eV line).  After: GREEN.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="${EXHALE_TEST_WORK:-$ROOT/build/tests/spectrum_type}"
CASE="$ROOT/backup/regression/wasp_full"

nfail=0
say() {  # say PASS|FAIL name measured reference
   echo "$1 $2 measured=$3 reference=$4 tol=0"
   [ "$1" = FAIL ] && nfail=$((nfail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   say FAIL spectrum_type_gate no_binary "$EXE"
   exit 1
fi

for d in stA stB stC; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$WORK/$d/"
   [ -f "$CASE/metals.inp" ] && cp "$CASE/metals.inp" "$WORK/$d/"
done
sed -i 's/^Spectrum type:.*/Spectrum type: Planck/' "$WORK/stB/input.inp" \
        "$WORK/stC/input.inp"
sed -i '/^Stellar Teff/d;/^Stellar radius/d' "$WORK/stC/input.inp"

run_case() {  # run_case <dir>; echoes the exit status
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   echo $?
}

rcA=$(run_case stA)
rcB=$(run_case stB)
rcC=$(run_case stC)

# --- A: the power-law run states its spectrum and its sub-Lyman source ----
if [ "$rcA" -eq 0 ]; then
   say PASS power_law_run_exit 0 0
else
   say FAIL power_law_run_exit "$rcA" 0
   sed -n '$p' "$WORK/stA/run.log"
fi
rep="$WORK/stA/EXHALE_setup.out"
if [ -f "$rep" ] && grep -qi 'Spectrum type' "$rep"; then
   say PASS power_law_report_states_type present present
else
   say FAIL power_law_report_states_type absent present
fi
if [ -f "$rep" ] && grep -qi 'below 13.6 eV' "$rep"; then
   say PASS power_law_report_states_sub_lyman_source present present
else
   say FAIL power_law_report_states_sub_lyman_source absent present
fi
if [ -f "$rep" ] && grep -qi 'Integrated grid flux' "$rep"; then
   say PASS power_law_report_states_grid_flux present present
else
   say FAIL power_law_report_states_grid_flux absent present
fi

# --- B: the Planck type runs and is recorded ------------------------------
if [ "$rcB" -eq 0 ]; then
   say PASS planck_run_exit 0 0
else
   say FAIL planck_run_exit "$rcB" 0
   tail -n 3 "$WORK/stB/run.log"
fi
repB="$WORK/stB/EXHALE_setup.out"
if [ -f "$repB" ] && grep -qi 'Planck' "$repB"; then
   say PASS planck_report_states_type present present
else
   say FAIL planck_report_states_type absent present
fi

# --- C: the Planck type without the stellar quantities is refused ---------
if [ "$rcC" -ne 0 ]; then
   say PASS planck_without_stellar_input_refused "$rcC" nonzero
else
   say FAIL planck_without_stellar_input_refused 0 nonzero
fi
if grep -qi 'Stellar Teff' "$WORK/stC/run.log" && \
   grep -qi 'Stellar radius' "$WORK/stC/run.log"; then
   say PASS planck_refusal_names_both_keys present present
else
   say FAIL planck_refusal_names_both_keys absent present
fi

if [ "$nfail" -gt 0 ]; then
   echo "spectrum_type_gate: $nfail assertion(s) FAILED"
   exit 1
fi
echo "spectrum_type_gate: all assertions PASSED"

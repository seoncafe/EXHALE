#!/bin/bash
# The base level of a run with a lower-atmosphere profile is stated once.
#
# QUANTITY UNDER TEST
#   The base pressure the run reports in two places:
#     EXHALE_setup.out, the "Base level:" line written by
#       src/modules/files_IO/write_setup_report.f90, which prints n0 from the
#       "Log10 lower boundary number density" key and the pressure it implies;
#     the run log, the "profile: p_base ->" line written when
#       apply_lower_atmosphere_profile (src/modules/files_IO/input_read.f90
#       2326) adopts the profile's matching level p_match_bar as p_base_bar.
#   They are the same physical quantity, the pressure at which the wind
#   starts, so they must be one number.  input_read.f90 1652 skips the
#   base-level consistency check whenever a profile is in use, so nothing in
#   the code compares them.
#
# WHAT IS RUN
#   A copy of backup/regression/lower_profile (input.inp, base.inp and
#   lower_atmosphere_profile.dat) in build/tests/grid_and_gates/lp, capped at
#   2 steps: both lines are written during startup.  The regression directory
#   is never written to.
#
# REFERENCE AND TOLERANCE
#   reference = the profile's p_match_bar; tolerance 1% relative, which is
#   looser than the round trip of the two printed formats (5 and 3 significant
#   figures) and far tighter than any disagreement that matters physically.
#
# EXPECTED AT HEAD 35d9dd5: RED, by a factor of about 20.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" base_level_single_statement
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}/lp"
CASE="$ROOT/backup/regression/lower_profile"

if [ ! -x "$EXE" ]; then
   echo "FAIL base_level_single_statement measured=no_binary reference=$EXE tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$CASE/lower_atmosphere_profile.dat" \
   "$WORK/"

( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=2 "$EXE" > run.log 2>&1 )
rc=$?
# Exit status 2 = the run declared a stationary state that the A2 certification
# refused (docs/a2_certification_contract_20260906.md section 8); the run wrote
# its outputs in full, which is what this gate reads. Only 1 (a Fortran error
# stop) or a signal is a failed run here.
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL base_level_single_statement_run measured=exit_$rc reference=exit_0 tol=0"
   echo "     see $WORK/run.log"
   exit 1
fi
echo "PASS base_level_single_statement_run measured=exit_0 reference=exit_0 tol=0"

python3 - "$WORK" <<'PY'
import re, sys

work = sys.argv[1]

p_setup = None
with open(work + "/EXHALE_setup.out") as fh:
    for line in fh:
        if "Base level" in line:
            print("  EXHALE_setup.out:%s" % line.rstrip())
            m = re.search(r"->\s*p\s*=\s*([0-9.eE+-]+)\s*bar", line)
            if m:
                p_setup = float(m.group(1))

p_profile = None
with open(work + "/run.log") as fh:
    for line in fh:
        if "profile: p_base" in line:
            print("  run log:%s" % line.rstrip())
            m = re.search(r"p_base\s*->\s*([0-9.eE+-]+)\s*bar", line)
            if m:
                p_profile = float(m.group(1))

if p_setup is None or p_profile is None:
    print("FAIL base_level_stated_once measured=unparsed reference=one_value "
          "tol=1.00e-02")
    sys.exit(1)

rel = abs(p_setup - p_profile) / p_profile
tol = 1.0e-2
ok = rel <= tol
print("%s base_level_stated_once measured=%.6e reference=%.6e tol=%.2e "
      "(setup/profile = %.4f)"
      % ("PASS" if ok else "FAIL", p_setup, p_profile, tol,
         p_setup / p_profile))
sys.exit(0 if ok else 1)
PY
exit $?

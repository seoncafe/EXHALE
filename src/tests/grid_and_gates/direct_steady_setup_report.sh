#!/bin/bash
# The direct steady route states its resolved configuration.
#
# QUANTITY UNDER TEST
#   EXHALE_setup.out, the record of the configuration a run resolved from its
#   input file. Every consumer of a run directory reads it to answer "which
#   physics was on, which scheme, which spectrum", and the run that writes no
#   report leaves a directory whose outputs cannot be attributed to a
#   configuration. The marching route writes it before it starts integrating;
#   the direct steady route (EXHALE_PTC=1) solves and stops inside its own
#   branch and never reaches that call, so the file it opens at startup stays
#   empty.
#
# WHAT IS RUN
#   backup/regression/roundtrip, copied out, with EXHALE_PTC=1 and the JFNK
#   arm of that route bounded to one outer iteration (EXHALE_PTC_JFNK=1,
#   EXHALE_JFNK_MAXIT=1): the report is written before the solve, so one
#   iteration is enough to reach the point under test and the row costs a
#   few seconds.
#
# ASSERTION
#   The run's EXHALE_setup.out carries the report: its header line and the
#   reconstruction-method line the marching route's report carries.
#
# MEASURED 2026-09-10: RED on the entry binary (the file is written and
# closed empty, 0 lines), GREEN with the report called on this route
# (66 lines).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}/direct_steady_setup"
CASE="$ROOT/backup/regression/roundtrip"

if [ ! -x "$EXE" ]; then
   echo "FAIL direct_steady_setup_report measured=no_binary reference=$EXE tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/"

( cd "$WORK" && env OMP_NUM_THREADS=1 EXHALE_PTC=1 EXHALE_PTC_JFNK=1 \
     EXHALE_JFNK_MAXIT=1 "$EXE" > run.log 2>&1 )
rc=$?
echo "  the direct steady route (EXHALE_PTC=1): exit $rc"
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL direct_steady_setup_report measured=exit_$rc reference=exit_0_or_2 tol=0"
   tail -n 5 "$WORK/run.log"
   exit 1
fi

n=$(wc -l < "$WORK/EXHALE_setup.out" 2>/dev/null || echo 0)
echo "  EXHALE_setup.out: $n lines"
fail=0
if [ "$n" -eq 0 ]; then
   echo "FAIL direct_steady_setup_report_written measured=0_lines reference=nonempty tol=0"
   fail=1
else
   echo "PASS direct_steady_setup_report_written measured=${n}_lines reference=nonempty tol=0"
fi
for pat in 'Simulation for' 'Reconstruction method'; do
   if grep -q "$pat" "$WORK/EXHALE_setup.out" 2>/dev/null; then
      echo "PASS direct_steady_setup_report_states_$(echo "$pat" | tr ' ' '_') measured=present reference=present tol=0"
   else
      echo "FAIL direct_steady_setup_report_states_$(echo "$pat" | tr ' ' '_') measured=absent reference=present tol=0"
      fail=1
   fi
done
exit $fail

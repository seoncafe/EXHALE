#!/bin/bash
# The run raises no arithmetic exception.
#
# QUANTITY UNDER TEST
#   Whether the binary at EXHALE_EXE completes a bounded hydrostatic_column
#   run without SIGFPE.  A production build masks every floating-point
#   exception, so an Inf or a NaN made in an intermediate is silent there and
#   is only visible in a binary built with -ffpe-trap; an integer division by
#   zero raises SIGFPE in every build, trap flags or not, because the x86
#   divide instruction faults on its own.
#
#   The gate is therefore worth what the binary handed to it is worth: with
#   the production EXHALE.x it asserts only that the case runs, and with
#
#     make OBJDIR=build_trap EXE=EXHALE_trap.x \
#          FFLAGS='-O0 -g -gz=none -fopenmp \
#                  -ffpe-trap=invalid,zero,overflow -fbacktrace'
#     EXHALE_EXE=$PWD/EXHALE_trap.x src/tests/grid_and_gates/run.sh fpe_traps
#
#   it asserts that no statement of the marching loop, the ionization sweep or
#   the post processing forms an undefined or out-of-range value.
#
# WHAT IS RUN
#   backup/regression/hydrostatic_column as shipped, 50 steps, in a scratch
#   copy under build/tests/grid_and_gates.  Nothing in backup/ is written to.
#   The case is an isothermal hydrostatic column whose post-processed profile
#   reaches T = 17 K at the top of the domain, which is what takes the
#   2^3S critical density out of the double-precision range.
#
# ASSERTION
#   The run exits 0 or 2 (2 = the stationary state was refused by the A2
#   certification, which still writes the outputs in full), and its log
#   carries no "Backtrace for this error".  A signal, hence a shell status
#   of 128 + n, is the failure.
#
# EXPECTED AT HEAD 35d9dd5, MEASURED 2026-09-06
#   production EXHALE.x                     PASS
#   -ffpe-trap build at -O3                 PASS
#   -ffpe-trap build at -O0                 FAIL, and not at a trap.  The
#      controller's duty-cycle test in src/EXHALE_main.f90 reads
#         n_err_every .gt. 0 .and. ... .and. mod(count, n_err_every) .eq. 0
#      and Fortran does not promise that the first operand short-circuits the
#      third.  At -O3 gfortran short-circuits it; at -O0 it evaluates the
#      integer mod with n_err_every = 0, which the divide instruction faults
#      on in any build, trap flags or not.  n_err_every is 0 whenever the
#      step-doubling duty cycle is off, its default outside phys mode.  The
#      fix is to nest the test rather than to rely on an evaluation order the
#      standard does not give.
#
#      A wasp_full run under the same build then reaches a second site of the
#      same kind, the four copies of
#         abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol
#      in the stall detector: du_prev is huge(1.0d0) on the first step of each
#      stage, so the quotient overflows to +Inf.  The comparison against +Inf
#      is false, which is the reset the code wants, so the value is discarded
#      -- and the statement still makes one.  Guarding the sentinel
#      (du_prev .lt. huge(1.0d0)) keeps the arithmetic of every other step.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" floating_point_traps_clean
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}/fpe_traps"
CASE="$ROOT/backup/regression/hydrostatic_column"

if [ ! -x "$EXE" ]; then
   echo "FAIL floating_point_traps measured=no_binary reference=$EXE tol=0"
   exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK/output"
cp "$CASE/input.inp" "$WORK/"

( cd "$WORK" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=50 "$EXE" \
     > run.log 2>&1 )
rc=$?

if [ $rc -gt 128 ]; then
   echo "  the run was killed by signal $((rc-128)); its last frames:"
   sed -n '/Backtrace for this error/,$p' "$WORK/run.log" | head -n 8
   echo "FAIL floating_point_traps measured=signal_$((rc-128)) reference=no_signal tol=0"
   exit 1
fi
if grep -q "Backtrace for this error" "$WORK/run.log"; then
   sed -n '/Backtrace for this error/,$p' "$WORK/run.log" | head -n 8
   echo "FAIL floating_point_traps measured=backtrace reference=no_backtrace tol=0"
   exit 1
fi
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   echo "FAIL floating_point_traps measured=exit_$rc reference=exit_0_or_2 tol=0"
   tail -n 5 "$WORK/run.log"
   exit 1
fi
echo "  50 steps of hydrostatic_column, exit $rc, no arithmetic exception"
echo "PASS floating_point_traps measured=exit_$rc reference=exit_0_or_2 tol=0"
exit 0

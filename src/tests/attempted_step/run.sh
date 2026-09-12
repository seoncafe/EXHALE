#!/usr/bin/env bash
# Build and run the attempted-step controller test (PLAN_20260906_rev2 B3a):
# the checkpoint of the whole step, the three counter classes, the step-size
# policy of a rejection, the integration-error estimate and the
# formation-energy table, plus whole-binary rows that inject a refusal after
# each of the twelve operations of the trial.
#
# WHY THIS SUITE IS NOT IN physics_probe.  Its driver reaches
# diffusive_photochemistry and ionization_equilibrium and through them the
# equilibrium solvers, and those call MINPACK (hybrd1, hybrd, dpmpar) and
# LAPACK (dgesvd) as free subroutines.  source_closure.py follows `use`
# statements, so it cannot find them, and the physics_probe link line carries
# no library.  This script adds the MINPACK sources and the LAPACK the
# Makefile resolves, and nothing else.
#
# EXHALE_TEST_OBJDIR selects another object directory, so that two builds
# of this suite can run side by side.
# EXHALE_OBJDIR names the build tree the generated build_stamp.f90 is
# taken from; without one the script writes its own.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
probe="$root/src/tests/physics_probe"
work="${EXHALE_TEST_OBJDIR:-$root/build/tests/attempted_step}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

# The LAPACK the Makefile resolves for the production link.  The library is
# needed because constrained_chemical_equilibrium calls dgesvd; nothing in
# this suite calls it.
prefix=$(dirname "$(dirname "$(command -v "$FC")")")
if [ -f "$prefix/lib/libopenblas.so" ]; then
  lapack=(-L"$prefix/lib" -lopenblas -Wl,-rpath,"$prefix/lib" -ldl)
else
  lapack=(-llapack -ldl)
fi

# MINPACK and the free subroutines the equilibrium solvers call by name,
# plus the shared verdict module of the Phase 0 test directories (it lives
# beside the physics_probe drivers, so source_closure.py does not find it).
extra=(
  "$probe/assertion_report.f90"
  "$root/src/modules/nonlinear_system_solver/dpmpar.f90"
  "$root/src/modules/nonlinear_system_solver/enorm.f90"
  "$root/src/modules/nonlinear_system_solver/qform.f90"
  "$root/src/modules/nonlinear_system_solver/qrfac.f90"
  "$root/src/modules/nonlinear_system_solver/r1mpyq.f90"
  "$root/src/modules/nonlinear_system_solver/r1updt.f90"
  "$root/src/modules/nonlinear_system_solver/dogleg.f90"
  "$root/src/modules/nonlinear_system_solver/fdjac1.f90"
  "$root/src/modules/nonlinear_system_solver/hybrd.f90"
  "$root/src/modules/nonlinear_system_solver/hybrd1.f90"
)

# utilities.f90 reads the revision strings of the build_stamp.f90 the
# Makefile generates in the object directory.  Without a build tree to read
# one from, write a stamp that says so; nothing tested here depends on its
# contents.  Whichever of the two it is, it is compiled below as the first
# file of the closure, so that build_stamp.mod stands in the module
# directory before utilities.f90 asks for it.
objdir="${EXHALE_OBJDIR:-$root/build}"
stamp="$objdir/build_stamp.f90"
if [ ! -f "$stamp" ]; then
  stamp="$work/build_stamp.f90"
  cat > "$stamp" <<'STAMP'
      module build_stamp
      ! Written by src/tests/attempted_step/run.sh when the
      ! Makefile has not generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

driver="$here/attempted_step_tests.f90"
order=$(python3 "$probe/source_closure.py" "$root" "$driver") || { echo "FAIL source_closure"; exit 1; }
# source_closure.py looks for the stamp in build/; compile the one selected
# above instead, first, since it uses no other module.
order=$(printf '%s\n' "$order" | grep -v '/build_stamp\.f90$')
order="$stamp
$order"

objects=()
for src in "${extra[@]}" $order; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  if [ "$src" != "$driver" ]; then objects+=("$obj"); fi
done

exe="$work/attempted_step_tests.x"
if ! "$FC" "${flags[@]}" "$work/attempted_step_tests.o" "${objects[@]}" \
       -o "$exe" "${lapack[@]}"; then
  echo "FAIL link attempted_step"
  exit 1
fi

status=0
"$exe" || status=1

# ---------------------------------------------------------------------- #
# THE WHOLE-BINARY ROWS.
#
# The unit driver above measures the checkpoint and the step-size policy on
# a constructed column. What it cannot measure is the controller running the
# fourteen operations of a real step, so the rows below inject a refusal
# after each of operations 1 to 12 of a real run and assert, for each:
#
#   * the state comes back (the run continues and produces the same output
#     it produces with no injection, because the retry succeeds);
#   * the clock does not advance for the rejected attempt: t_phys equals the
#     sum of the ACCEPTED global dt, over a run in which attempts were
#     refused;
#   * the physical ledgers do not carry the refused attempt, and the attempt
#     statistics do;
#   * exhaustion stops in phys mode and continues, flagged, in init.
#
# EXHALE_ATTEMPTED_STEP_EXE names the binary. Without it the rows are
# skipped and only the unit driver runs, so this suite never builds the
# production binary itself (other workers hold that tree).
# ---------------------------------------------------------------------- #

exe_run=${EXHALE_ATTEMPTED_STEP_EXE:-}
work_run=${EXHALE_ATTEMPTED_STEP_WORK:-$work/runs}
case_dir=${EXHALE_ATTEMPTED_STEP_CASE:-$root/backup/regression/hydrostatic_column}

say() { echo "$1 $2 measured=$3 reference=$4 tol=$5"; }
fail() { status=1; }

run_case() {
  # $1 = tag, rest = environment assignments. case_dir may be overridden for
  # one call by setting it in the call's own environment.
  local tag=$1; shift
  local d="$work_run/$tag"
  rm -rf "$d"; mkdir -p "$d/output"
  cp "$case_dir/input.inp" "$d/"
  for f in metals.inp base.inp opacity.inp; do
    [ -f "$case_dir/$f" ] && cp "$case_dir/$f" "$d/"
  done
  [ -n "${RUN_MODE:-}" ] && { grep -v "^Run mode" "$d/input.inp" > "$d/i2"; \
      mv "$d/i2" "$d/input.inp"; echo "Run mode: $RUN_MODE" >> "$d/input.inp"; }
  ( cd "$d" && env OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=${STEPS:-20} \
        EXHALE_STEP_CLOCK=1 "$@" "$exe_run" > run.log 2>&1; \
    echo "exit=$?" >> run.log )
  echo "$d"
}

if [ -n "$exe_run" ] && [ -x "$exe_run" ]; then
  mkdir -p "$work_run"

  # -- the reference run: nothing injected ----------------------------- #
  ref=$(STEPS=20 run_case reference)
  ref_acc=$(grep -c "step-clock" "$ref/run.log")
  say PASS the_reference_run_takes_its_steps "$ref_acc" 20 0
  [ "$ref_acc" = "20" ] || { say FAIL the_reference_run_takes_its_steps \
      "$ref_acc" 20 0; fail; }

  # -- one injected refusal after each operation of the trial ---------- #
  for op in 1 2 3 4 5 6 7 8 9 10 11 12; do
    d=$(STEPS=20 run_case "op$op" EXHALE_REJECT_AFTER_OP=$op \
                                  EXHALE_REJECT_AT_STEP=5 \
                                  EXHALE_REJECT_ATTEMPTS=1)
    ex=$(grep -o "exit=[0-9]*" "$d/run.log" | tail -n 1)
    ref_txt="$ref/output/Hydro_ioniz.txt"
    got_txt="$d/output/Hydro_ioniz.txt"
    # The run completed at all.
    if [ ! -f "$got_txt" ]; then
      say FAIL "op${op}_the_run_completes" missing present 0; fail; continue
    fi
    # It refused exactly once, at the operation asked for.
    nref=$(grep -c "REFUSED at" "$d/run.log" || true)
    if [ "$nref" = "1" ]; then
      say PASS "op${op}_refused_exactly_once" "$nref" 1 0
    else
      say FAIL "op${op}_refused_exactly_once" "$nref" 1 0; fail
    fi
    # The clock did not advance for it: the accepted steps are still 20.
    nacc=$(grep -c "step-clock" "$d/run.log" || true)
    if [ "$nacc" = "20" ]; then
      say PASS "op${op}_the_rejected_attempt_is_not_a_step" "$nacc" 20 0
    else
      say FAIL "op${op}_the_rejected_attempt_is_not_a_step" "$nacc" 20 0; fail
    fi
    # The restore put every item back: the controller asserts its own
    # round trip on every rejection and stops loudly if it did not.
    if grep -q "RESTORE INCOMPLETE" "$d/run.log"; then
      say FAIL "op${op}_the_restore_is_complete" incomplete complete 0; fail
    else
      say PASS "op${op}_the_restore_is_complete" complete complete 0
    fi
    # The outer attempt count carries it and the accepted count does not.
    natt=$(grep "outer attempts:" "$d/run.log" | \
           sed "s/.*outer attempts: \\([0-9]*\\),.*/\\1/")
    nrej=$(grep "outer attempts:" "$d/run.log" | \
           sed "s/.*of which \\([0-9]*\\) were refused.*/\\1/")
    if [ "$nrej" = "1" ] && [ "$natt" = "21" ]; then
      say PASS "op${op}_the_attempt_statistics_carry_it" "${natt}_${nrej}" 21_1 0
    else
      say FAIL "op${op}_the_attempt_statistics_carry_it" "${natt}_${nrej}" 21_1 0
      fail
    fi
  done

  # -- t_phys is the sum of the accepted dt, with refusals in the run -- #
  d=$(RUN_MODE=phys STEPS=20 run_case clock EXHALE_REJECT_AFTER_OP=7 \
                                            EXHALE_REJECT_AT_STEP=5 \
                                            EXHALE_REJECT_ATTEMPTS=1 \
                                            EXHALE_ERR_EVERY=0)
  sum=$(grep "step-clock" "$d/run.log" | \
        sed "s/.* dt= *\\([-0-9.E+]*\\) .*/\\1/" | \
        python3 -c "import sys; print(repr(sum(float(x) for x in sys.stdin)))")
  last=$(grep "step-clock" "$d/run.log" | tail -n 1 | \
         sed "s/.* t_phys= *\\([-0-9.E+]*\\).*/\\1/")
  same=$(python3 -c "print('equal' if abs($sum-$last) <= 1e-12*abs($last) else 'differs')")
  if [ "$same" = "equal" ]; then
    say PASS t_phys_is_the_sum_of_the_accepted_dt equal equal 1e-12
  else
    say FAIL t_phys_is_the_sum_of_the_accepted_dt differs equal 1e-12; fail
  fi

  # -- exhaustion: exits 2 in BOTH modes, with the retry history ------- #
  # A step refused at every dt down to the floor is the same statement in a
  # physical run and in a relaxation: the code could not produce an
  # acceptable state. The status is the 2 certification.f90 uses for a
  # refused stationary claim, so one status means one thing to a caller.
  for m in phys init; do
    d=$(RUN_MODE=$m STEPS=20 run_case "exhaust_$m" \
          EXHALE_REJECT_AFTER_OP=7 EXHALE_REJECT_AT_STEP=5 \
          EXHALE_REJECT_ATTEMPTS=99 EXHALE_ERR_EVERY=0)
    ex=$(grep -o "exit=[0-9]*" "$d/run.log" | tail -n 1)
    if grep -q "ATTEMPTED STEP EXHAUSTED" "$d/run.log" && \
       [ "$ex" = "exit=2" ]; then
      say PASS "exhaustion_exits_2_in_$m" "$ex" exit=2 0
    else
      say FAIL "exhaustion_exits_2_in_$m" "$ex" exit=2 0; fail
    fi
    # The retry history: one line per attempt, all nine of them.
    nh=$(grep -c "^     attempt " "$d/run.log" || true)
    if [ "$nh" = "9" ]; then
      say PASS "exhaustion_prints_the_retry_history_in_$m" "$nh" 9 0
    else
      say FAIL "exhaustion_prints_the_retry_history_in_$m" "$nh" 9 0; fail
    fi
    if [ "$m" = "phys" ]; then
      nacc=$(grep -c "step-clock" "$d/run.log" || true)
      if [ "$nacc" = "5" ]; then
        say PASS the_exhausted_step_is_not_accepted "$nacc" 5 0
      else
        say FAIL the_exhausted_step_is_not_accepted "$nacc" 5 0; fail
      fi
    fi
  done

  # -- an energy failure refuses the step, in BOTH run modes ---------- #
  # The injected refusals above say where a refusal is turned into a
  # rejection; this row says that the energy update raises one with its own
  # reason, and that it does so in init mode as well: a failed update
  # produced no temperature, so there is no iterate to march on and the
  # leniency of the initialization mode does not apply. One leading call is
  # forced to report a failure, so the retry at half dt succeeds and the run
  # still takes its twenty steps.
  for m in phys init; do
    d=$(RUN_MODE=$m STEPS=20 run_case "energy_$m" \
          EXHALE_ENERGY_FAIL_AT_STEP=5 EXHALE_ENERGY_FAIL_CALLS=1 \
          EXHALE_ERR_EVERY=0)
    if grep -q "REFUSED at the energy update: the energy update failed" \
         "$d/run.log"; then
      say PASS "an_energy_failure_refuses_the_step_in_$m" refused refused 0
    else
      say FAIL "an_energy_failure_refuses_the_step_in_$m" not_refused \
          refused 0; fail
    fi
    # It reported and did not stop: the run went on and took its steps.
    nacc=$(grep -c "step-clock" "$d/run.log" || true)
    nrep=$(grep -c "energy update FAILURE" "$d/run.log" || true)
    if [ "$nrep" -ge 1 ] && [ "$nacc" = "20" ]; then
      say PASS "the_energy_failure_is_reported_and_retried_in_$m" \
          "reported${nrep}_accepted${nacc}" reported1_accepted20 0
    else
      say FAIL "the_energy_failure_is_reported_and_retried_in_$m" \
          "reported${nrep}_accepted${nacc}" reported1_accepted20 0; fail
    fi
    # The end-of-run breakdown names it.
    if grep -q "the energy update failed: 1" "$d/run.log"; then
      say PASS "the_reason_breakdown_counts_it_in_$m" 1 1 0
    else
      say FAIL "the_reason_breakdown_counts_it_in_$m" 0 1 0; fail
    fi
  done

  # -- the step-doubling estimate is taken and reported ---------------- #
  d=$(RUN_MODE=phys STEPS=25 run_case err_estimate EXHALE_ERR_EVERY=10)
  n_est=$(grep -c "integration-error estimate:" "$d/run.log" || true)
  if [ "$n_est" -ge 1 ]; then
    say PASS the_step_doubling_estimate_is_taken "$n_est" ">=1" 0
  else
    say FAIL the_step_doubling_estimate_is_taken "$n_est" ">=1" 0; fail
  fi
  # ... and it costs what the duty cycle says it costs.
  if grep -q "seconds spent in them" "$d/run.log"; then
    say PASS the_estimate_reports_its_cost reported reported 0
  else
    say FAIL the_estimate_reports_its_cost "not_reported" reported 0; fail
  fi

  # -- the thermal identity of the source step ------------------------- #
  # The instrument closes the SOURCE STEP: the thermal energy the source
  # operators produced against the thermal source they were given, with the
  # transport contribution MEASURED from the two marks the marching loop
  # takes on either side of the coupled source step. Without those marks
  # there is no measured transport term, the printed residual is the whole
  # step's imbalance instead of a closure, and the gate is inert; the rows
  # below say so rather than reporting a number that is not a closure.
  main="$root/src/EXHALE_main.f90"
  if grep -q "attempted_step_thermal_energy_before_sources" "$main" && \
     grep -q "attempted_step_thermal_energy_after_sources"  "$main"; then
    d=$(RUN_MODE=phys STEPS=20 run_case source_energy EXHALE_ERR_EVERY=0)
    worst=$(grep "energy-identity" "$d/run.log" | \
            sed "s/.*residual\/scale= *\([-0-9.E+]*\).*/\1/" | \
            python3 -c "import sys; v=[abs(float(x)) for x in sys.stdin]; print(repr(max(v)) if v else 'none')")
    if [ "$worst" = "none" ]; then
      say FAIL the_identity_is_printed_every_step none printed 0; fail
    else
      ok=$(python3 -c "print('yes' if $worst <= 1e-9 else 'no')")
      if [ "$ok" = "yes" ]; then
        say PASS the_source_step_closes_its_thermal_energy "$worst" 0 1e-9
      else
        say FAIL the_source_step_closes_its_thermal_energy "$worst" 0 1e-9; fail
      fi
    fi
    # ... and a thermal energy the sources did not produce is refused.
    d=$(RUN_MODE=phys STEPS=20 run_case source_energy_injected \
          EXHALE_ERR_EVERY=0 EXHALE_INJECT_SOURCE_ENERGY=1.0d-6)
    if grep -q "did not close its thermal energy" "$d/run.log"; then
      say PASS an_injected_inconsistency_is_refused refused refused 0
    else
      say FAIL an_injected_inconsistency_is_refused not_refused refused 0; fail
    fi
  else
    echo "SKIP the thermal-identity rows: EXHALE_main.f90 does not call" \
         "attempted_step_thermal_energy_before_sources /" \
         "..._after_sources, so the transport contribution is not measured" \
         "and the gate is inert"
  fi

  # -- an operation outside the trial is refused as a request ---------- #
  d=$(STEPS=3 run_case outside_trial EXHALE_REJECT_AFTER_OP=13)
  if grep -q "OUTSIDE THE TRIAL" "$d/run.log"; then
    say PASS an_operation_outside_the_trial_is_refused refused refused 0
  else
    say FAIL an_operation_outside_the_trial_is_refused accepted refused 0; fail
  fi

  # ------------------------------------------------------------------- #
  # THE ERROR-CONTROL MACROSTEP ON THE PRODUCTION ROUTE
  # (PLAN_20260909_rev1 item N16a; review findings F1, F2, F3).
  #
  # The rows above measure a REFUSAL inside one trial. The rows below
  # measure the transaction the sampled step-doubling comparison is: the
  # full pass and the two half passes covering ONE macro-interval H, the
  # estimate of their difference as a rejection criterion, and the clock
  # that may move only for a macrostep that was accepted.
  #
  # Every run below prescribes its step (EXHALE_AS_FIXED_DT), so that the
  # interval, and with it the size of the integration error, is a property
  # of the row and not of the state the CFL condition happens to give. At
  # the prescribed 0.1 s the case's own estimate is 0.2 to 0.4, well inside
  # the tolerance, so a rejection in these rows is the injected one and
  # nothing else. EXHALE_ERR_EVERY=5 puts one sampled macrostep at step 5
  # and another at step 10 of a twelve-step run.
  # ------------------------------------------------------------------- #

  errenv=(EXHALE_ERR_EVERY=5 EXHALE_AS_FIXED_DT=0.1)

  # The interval every pass of a sampled macrostep covers: the full pass
  # carries H and the two halves 0.5*H of the SAME H (F2's endpoint).
  d=$(RUN_MODE=phys STEPS=12 run_case macro_endpoint "${errenv[@]}")
  bad=$(python3 - "$d/run.log" <<'PYEND'
import re, sys
bad = 0
seen = 0
for line in open(sys.argv[1]):
    m = re.search(r"err-pass: step=(\d+) pass=(.+?) dt= *(\S+) H= *(\S+)", line)
    if not m:
        continue
    seen += 1
    pas, dt, H = m.group(2), float(m.group(3)), float(m.group(4))
    want = H if pas == "the full step" else 0.5 * H
    if abs(dt - want) > 1e-12 * H:
        bad += 1
print(bad if seen else "no_pass_line")
PYEND
)
  if [ "$bad" = "0" ]; then
    say PASS the_three_passes_cover_one_macro_interval 0 0 1e-12
  else
    say FAIL the_three_passes_cover_one_macro_interval "$bad" 0 1e-12; fail
  fi

  # F1: the estimate REJECTS. The hook multiplies the first estimate of the
  # run, so the macrostep is refused for its integration error with every
  # local solve succeeding.
  d=$(RUN_MODE=phys STEPS=12 run_case forced_error "${errenv[@]}" \
        EXHALE_AS_ERR_INJECT=1.0d6)
  nrej=$(grep -c "macrostep-reject:.*the integration error was too big" \
         "$d/run.log" || true)
  if [ "$nrej" = "1" ]; then
    say PASS forced_error_rejects_the_macrostep "$nrej" 1 0
  else
    say FAIL forced_error_rejects_the_macrostep "$nrej" 1 0; fail
  fi
  # ... and the run took every step it was asked for, so the rejection was
  # a retry and not a stop.
  nacc=$(grep -c "step-clock" "$d/run.log" || true)
  if [ "$nacc" = "12" ]; then
    say PASS the_rejected_macrostep_is_retried_not_fatal "$nacc" 12 0
  else
    say FAIL the_rejected_macrostep_is_retried_not_fatal "$nacc" 12 0; fail
  fi

  # The clock does not carry a rejected macrostep: t_phys is the sum of the
  # accepted increments, and no accepted increment is an interval that was
  # refused.
  clock_is_the_sum() {
    # $1 = run directory, $2 = row name
    local dd=$1 nm=$2
    local verdict
    verdict=$(python3 - "$dd/run.log" <<'PYEND'
import re, sys
acc, rej = [], []
last = None
for line in open(sys.argv[1]):
    m = re.search(r"step-clock: count=(\d+) .* dt= *(\S+) t_phys= *(\S+)", line)
    if m:
        acc.append((int(m.group(1)), float(m.group(2))))
        last = float(m.group(3))
        continue
    m = re.search(r"macrostep-reject: step=(\d+) .* H= *(\S+) H_new=", line)
    if m:
        rej.append((int(m.group(1)), float(m.group(2))))
if not acc or last is None:
    print("no_clock")
    sys.exit()
tot = sum(d for _, d in acc)
if abs(tot - last) > 1e-12 * abs(last):
    print("sum_differs")
    sys.exit()
if not rej:
    print("no_rejection")
    sys.exit()
for c, H in rej:
    for ca, da in acc:
        if ca == c and abs(da - H) <= 1e-12 * H:
            print("refused_interval_accepted")
            sys.exit()
print("equal")
PYEND
)
    if [ "$verdict" = "equal" ]; then
      say PASS "$nm" equal equal 1e-12
    else
      say FAIL "$nm" "$verdict" equal 1e-12; fail
    fi
  }
  clock_is_the_sum "$d" rejected_macrostep_leaves_the_clock

  # The retry is SHORTER than the interval that was refused, and the clock
  # accepted the interval the retried full pass actually covered.
  verdict=$(python3 - "$d/run.log" <<'PYEND'
import re, sys
rej = {}
acc = {}
for line in open(sys.argv[1]):
    m = re.search(r"macrostep-reject: step=(\d+) .* H= *(\S+) H_new= *(\S+)", line)
    if m:
        rej.setdefault(int(m.group(1)), []).append(
            (float(m.group(2)), float(m.group(3))))
    m = re.search(r"step-clock: count=(\d+) .* dt= *(\S+) t_phys=", line)
    if m:
        acc[int(m.group(1))] = float(m.group(2))
if not rej:
    print("no_rejection")
    sys.exit()
for c, pairs in rej.items():
    for H, Hnew in pairs:
        if not Hnew < H:
            print("retry_not_shorter")
            sys.exit()
    if c not in acc:
        print("no_accepted_step_there")
        sys.exit()
    if abs(acc[c] - pairs[-1][1]) > 1e-12 * pairs[-1][1]:
        print("clock_is_not_the_retried_interval")
        sys.exit()
print("shorter_and_accepted")
PYEND
)
  if [ "$verdict" = "shorter_and_accepted" ]; then
    say PASS retry_is_shorter_and_then_accepted shorter_and_accepted \
        shorter_and_accepted 1e-12
  else
    say FAIL retry_is_shorter_and_then_accepted "$verdict" \
        shorter_and_accepted 1e-12; fail
  fi

  # F2: a refusal INSIDE a pass rejects the whole macrostep, so the three
  # passes of the retry share the endpoint the clock then records. One row
  # per pass, the refusal aimed at that pass.
  pass_name() { case $1 in 1) echo "the full step";; 2) echo "the first half";;
                           3) echo "the second half";; esac; }
  for ip in 1 2 3; do
    case $ip in
      1) row=retry_inside_the_full_pass_keeps_a_common_endpoint;;
      2) row=rejection_in_the_first_half;;
      3) row=rejection_in_the_second_half;;
    esac
    d=$(RUN_MODE=phys STEPS=12 run_case "in_pass$ip" "${errenv[@]}" \
          EXHALE_REJECT_AFTER_OP=7 EXHALE_REJECT_AT_STEP=5 \
          EXHALE_REJECT_ATTEMPTS=1 EXHALE_AS_REJECT_IN_PASS=$ip)
    want=$(pass_name $ip)
    if grep -q "macrostep-reject: step=5 pass=$want reason=" "$d/run.log"; then
      say PASS "${row}_names_the_pass" "$want" "$want" 0
    else
      say FAIL "${row}_names_the_pass" absent "$want" 0; fail
    fi
    # The three passes of every attempt still share one interval, the run
    # reached its steps, and the clock is the sum of what it accepted.
    bad=$(python3 - "$d/run.log" <<'PYEND'
import re, sys
bad = 0
for line in open(sys.argv[1]):
    m = re.search(r"err-pass: step=(\d+) pass=(.+?) dt= *(\S+) H= *(\S+)", line)
    if not m:
        continue
    pas, dt, H = m.group(2), float(m.group(3)), float(m.group(4))
    want = H if pas == "the full step" else 0.5 * H
    if abs(dt - want) > 1e-12 * H:
        bad += 1
print(bad)
PYEND
)
    nacc=$(grep -c "step-clock" "$d/run.log" || true)
    if [ "$bad" = "0" ] && [ "$nacc" = "12" ]; then
      say PASS "${row}_keeps_the_endpoint_and_the_steps" "0_${nacc}" 0_12 1e-12
    else
      say FAIL "${row}_keeps_the_endpoint_and_the_steps" "${bad}_${nacc}" \
          0_12 1e-12; fail
    fi
    clock_is_the_sum "$d" "${row}_leaves_the_clock"
  done

  # EXHAUSTION IN A PASS: refused in the same pass at every interval the
  # policy offers. The named failure stops the run with status 2, says
  # which pass refused it, and the state it restored is the last one the
  # clock accepted: the same hash, not merely a plausible profile.
  for ip in 2 3; do
    d=$(RUN_MODE=phys STEPS=12 run_case "exhaust_pass$ip" "${errenv[@]}" \
          EXHALE_REJECT_AFTER_OP=7 EXHALE_REJECT_AT_STEP=5 \
          EXHALE_REJECT_ATTEMPTS=99 EXHALE_AS_REJECT_IN_PASS=$ip)
    want=$(pass_name $ip)
    ex=$(grep -o "exit=[0-9]*" "$d/run.log" | tail -n 1)
    if grep -q "the pass that refused it: $want of the error-control" \
         "$d/run.log" && [ "$ex" = "exit=2" ]; then
      say PASS "exhaustion_in_${want// /_}_is_named" "$ex" exit=2 0
    else
      say FAIL "exhaustion_in_${want// /_}_is_named" "$ex" exit=2 0; fail
    fi
    got=$(grep "restored-state: state=" "$d/run.log" | tail -n 1 | \
          sed "s/.*state= *//")
    ref=$(grep "step-clock" "$d/run.log" | tail -n 1 | sed "s/.*state= *//")
    if [ -n "$got" ] && [ "$got" = "$ref" ]; then
      say PASS "exhaustion_in_${want// /_}_restores_the_accepted_state" \
          "$got" "$ref" 0
    else
      say FAIL "exhaustion_in_${want// /_}_restores_the_accepted_state" \
          "${got:-absent}" "$ref" 0; fail
    fi
  done

  # EXHAUSTION ON THE INTEGRATION ERROR ITSELF: the estimate of every
  # interval the policy offers is above one, so the macrostep is refused
  # for its integration error until the budget is spent.
  d=$(RUN_MODE=phys STEPS=12 run_case exhaust_error "${errenv[@]}" \
        EXHALE_AS_ERR_INJECT=1.0d12 EXHALE_AS_ERR_INJECT_COUNT=99)
  ex=$(grep -o "exit=[0-9]*" "$d/run.log" | tail -n 1)
  if grep -q "the pass that refused it: the error estimate" "$d/run.log" \
     && [ "$ex" = "exit=2" ]; then
    say PASS exhaustion_on_the_integration_error_is_named "$ex" exit=2 0
  else
    say FAIL exhaustion_on_the_integration_error_is_named "$ex" exit=2 0; fail
  fi
  got=$(grep "restored-state: state=" "$d/run.log" | tail -n 1 | \
        sed "s/.*state= *//")
  ref=$(grep "step-clock" "$d/run.log" | tail -n 1 | sed "s/.*state= *//")
  if [ -n "$got" ] && [ "$got" = "$ref" ]; then
    say PASS exhaustion_on_the_integration_error_restores_the_state \
        "$got" "$ref" 0
  else
    say FAIL exhaustion_on_the_integration_error_restores_the_state \
        "${got:-absent}" "$ref" 0; fail
  fi

  # AN UNRESOLVED ESTIMATE IS NAMED AND JUDGES NOTHING (item N16b). The
  # hook makes every class that carries a difference report an inner
  # nonlinear error that swallows it, while the error injection puts the
  # reported estimate six decades above the tolerance. A difference that is
  # the solvers' own error may not refuse a step, so the macrostep is
  # adopted, the outcome is named, and the run goes on.
  d=$(RUN_MODE=phys STEPS=12 run_case unresolved "${errenv[@]}" \
        EXHALE_AS_ERR_INJECT=1.0d6 EXHALE_AS_ERR_INJECT_COUNT=99 \
        EXHALE_AS_FORCE_UNRESOLVED=99)
  nun=$(grep -c "integration-error UNRESOLVED:" "$d/run.log" || true)
  nrej=$(grep -c "macrostep-reject:.*the integration error was too big" \
         "$d/run.log" || true)
  nacc=$(grep -c "step-clock" "$d/run.log" || true)
  if [ "$nun" -ge 1 ] && [ "$nrej" = "0" ] && [ "$nacc" = "12" ]; then
    say PASS unresolved_estimate_is_named_not_judged \
        "named${nun}_rejected${nrej}_accepted${nacc}" \
        "named_ge1_rejected0_accepted12" 0
  else
    say FAIL unresolved_estimate_is_named_not_judged \
        "named${nun}_rejected${nrej}_accepted${nacc}" \
        "named_ge1_rejected0_accepted12" 0; fail
  fi
  # ... and the line says which class could not be trusted and by how much.
  if grep -q "integration-error UNRESOLVED:.*first unresolved class=.*inner/difference=" \
       "$d/run.log"; then
    say PASS the_unresolved_outcome_names_its_class named named 0
  else
    say FAIL the_unresolved_outcome_names_its_class unnamed named 0; fail
  fi
  # The end-of-run block says how much of the trajectory was not
  # error-controlled, so that a run cannot be read as one that was.
  if grep -q "were UNRESOLVED" "$d/run.log"; then
    say PASS the_run_reports_its_unresolved_estimates reported reported 0
  else
    say FAIL the_run_reports_its_unresolved_estimates absent reported 0; fail
  fi

  # THE ESTIMATE NAMES THE CLASS THAT SETS THE REDUCTION AND ITS ORDER.
  # The retained-step factor and the exponent of the reduction are the
  # order of the integrator that carries the deciding class, so a run that
  # refused a step has to say which class refused it and at what order,
  # or the interval it asks for next cannot be read back.
  nd=$(grep -c "integration-error estimate:.* deciding=.* p=[0-9]" \
       "$d/run.log" || true)
  ne_=$(grep -c "integration-error estimate:" "$d/run.log" || true)
  if [ "$nd" = "$ne_" ] && [ "$nd" -ge 1 ]; then
    say PASS the_estimate_names_the_deciding_class_and_its_order \
        "${nd}_of_${ne_}" all 0
  else
    say FAIL the_estimate_names_the_deciding_class_and_its_order \
        "${nd}_of_${ne_}" all 0; fail
  fi

  # THE TEMPORAL ORDER OF THE ROW CLASSES (item N16b deliverable 3). The
  # macrostep is refused at every interval the policy offers, so the same
  # entry state is retaken at H, then 0.2 H, then 0.04 H, and every class
  # is measured at each. For a complete split of order p the difference
  # scales as H^(p+1); the rows below read p + 1 off consecutive levels of
  # the SPECIES class and of the class the electron fraction and the
  # temperature share, over three levels each.
  d=$(RUN_MODE=phys STEPS=8 run_case order_by_class "${errenv[@]}" \
        EXHALE_AS_ERR_INJECT=1.0d12 EXHALE_AS_ERR_INJECT_COUNT=99)
  order_of_class() {
    # $1 = run directory, $2 = class name as printed, $3 = how many
    # consecutive ratios to read, $4 = how many leading levels to drop.
    # Prints the ratios space separated, or a word.
    python3 - "$1/run.log" "$2" "$3" "$4" <<'PYEND'
import re, sys, math
log, want, nwant = sys.argv[1], sys.argv[2], int(sys.argv[3])
nskip = int(sys.argv[4])
H = None
levels = []
for line in open(log):
    m = re.search(r"err-class: (.+?) e= *(\S+) diff= *(\S+) ", line)
    if m and m.group(1).strip() == want and H is not None:
        levels.append((H, float(m.group(3))))
        continue
    m = re.search(r"err-pass: step=\d+ pass=the full step dt= *(\S+) H=", line)
    if m:
        H = float(m.group(1))
levels = levels[nskip:]
if len(levels) < nwant + 1:
    print("too_few_levels")
    sys.exit()
out = []
for i in range(nwant):
    H0, d0 = levels[i]
    H1, d1 = levels[i + 1]
    if d0 <= 0.0 or d1 <= 0.0 or H0 == H1:
        print("no_ratio")
        sys.exit()
    out.append(math.log(d0 / d1) / math.log(H0 / H1))
print(" ".join("%.4f" % v for v in out))
PYEND
  }
  # TWO CLASSES, TWO ORDERS, AND THE REFERENCE OF EACH IS THE ORDER OF THE
  # INTEGRATOR THAT CARRIES IT. These two ladders are the measurement
  # err_order_p(class) is declared from, so a change that breaks either
  # order breaks the retained-step factor and the reduction exponent of
  # that class with it.
  #
  #  * the electron fraction and the temperature ride on the conservative
  #    update, which is the three-stage Runge-Kutta of the hydrodynamics.
  #    A local error of order H^4 makes the difference between one step of
  #    H and two of H/2 fall as H^4, so p + 1 = 4. MEASURED 3.9867 on the
  #    factor-5 level below the coarsest, and 3.9962 / 3.9926 on the
  #    controller-formula levels of the ladder of the same case with
  #    EXHALE_AS_ERR_INJECT=1.0d6 (item N16d).
  #  * the transported species fractions ride on their own advection
  #    stage, which is split from the source at first order, so their
  #    difference falls as H^2 and p + 1 = 2. MEASURED 1.9987 and 1.9799.
  #
  # THE TOLERANCE IS THE MEASUREMENT'S OWN: 5 percent admits what was
  # measured and refuses anything that is not the named order.
  #
  # THE COARSEST INTERVAL OF THE LADDER IS NOT IN THE ASYMPTOTIC REGIME
  # for the rows the base cell holds: at H = 0.1 s the electron and
  # temperature rows give 4.85 and the conservative rows 5.35, while the
  # finer levels of the same ladder give 4. The species rows are in the
  # regime from the first level, so nothing is dropped for them. Below the
  # named level the electron and temperature difference reaches the
  # round-off of its own normalization (1e-15 at H = 8e-4 s on this case)
  # and stops falling, which is why one ratio is read for that class and
  # two for the species. Both numbers are stated here, per class.
  for spec in "the transported species fractions:0:2:2.0" \
              "the electron fraction and T:1:1:4.0"; do
    IFS=: read -r cls nskip nwant pref <<< "$spec"
    ratios=$(order_of_class "$d" "$cls" "$nwant" "$nskip")
    ok=$(RATIOS="$ratios" PREF="$pref" python3 -c "
import os
r = os.environ['RATIOS'].split()
ref = float(os.environ['PREF'])
try:
    v = [float(x) for x in r]
except ValueError:
    print('no_ratio')
else:
    print('yes' if v and all(abs(x - ref) <= 0.05 * ref for x in v) else 'no')")
    nm=$(echo "$cls" | tr ' ' '_')
    if [ "$ok" = "yes" ]; then
      say PASS "temporal_order_of_${nm}" "$ratios" "$pref" 0.05
    else
      say FAIL "temporal_order_of_${nm}" "${ratios:-none}" "$pref" 0.05; fail
    fi
  done

  # A REJECTED PASS MUST LEAVE NOTHING BEHIND (item N16b deliverable 4).
  # Three runs that differ ONLY in which pass of the sampled macrostep was
  # refused all end by adopting the SAME computation from the SAME restored
  # state, so their outputs have to agree to round-off.
  #
  # WHAT THE ROW HOLDS TO ACCOUNT are the two module states that a restore
  # has to put back and that are not in the checkpoint because they are
  # functions of (u, f_sp): the caloric composition of caloric_eos
  # (caloric_mixture_active, nk_per_mass, x_h2, molecular_cell), which the
  # energy-to-pressure map reads, and the base-face states of BC_Apply
  # (base_face_W, base_face_lower_W), which Rec_BC puts into the left slot
  # of the base face in every Reconstruct. rebuild_state_from_checkpoint
  # rebuilds both, composition first and the boundary second, and this row
  # is what says it does. Both output files are compared, so a difference
  # in the composition columns cannot hide behind the hydrodynamic ones.
  leak=""
  for ip in 1 2 3; do
    d=$(RUN_MODE=phys STEPS=6 run_case "leak_p$ip" EXHALE_ERR_EVERY=5 \
          EXHALE_AS_FIXED_DT=0.1 EXHALE_REJECT_AFTER_OP=7 \
          EXHALE_REJECT_AT_STEP=5 EXHALE_REJECT_ATTEMPTS=1 \
          EXHALE_AS_REJECT_IN_PASS=$ip)
    leak="$leak $d/output"
  done
  worst=$(python3 - $leak <<'PYEND'
import sys
def cols(path):
    out = []
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.split()
        if f:
            out.append([float(x) for x in f])
    return out
w = 0.0
for name in ('Hydro_ioniz.txt', 'Ion_species.txt'):
    a, b, c = (cols(d + '/' + name) for d in sys.argv[1:4])
    for x, y in ((a, b), (a, c), (b, c)):
        for ra, rb in zip(x, y):
            for va, vb in zip(ra, rb):
                s = max(abs(va), abs(vb))
                if s > 0.0:
                    w = max(w, abs(va - vb) / s)
print(repr(w))
PYEND
)
  ok=$(WORST="$worst" python3 -c "
import os
print('yes' if float(os.environ['WORST']) <= 1e-12 else 'no')")
  if [ "$ok" = "yes" ]; then
    say PASS rejected_pass_leaves_no_state_behind "$worst" 0 1e-12
  else
    say FAIL rejected_pass_leaves_no_state_behind "$worst" 0 1e-12
    echo "  Which pass a macrostep was refused in changed the state it" \
         "finally adopts, so something an attempt wrote survived the" \
         "restore. The two known channels are the caloric composition of" \
         "caloric_eos, read by the energy-to-pressure map, and the" \
         "base-face states of BC_Apply, read by Rec_BC; both are rebuilt" \
         "by rebuild_state_from_checkpoint, so a failure here is either" \
         "that order broken or a third quantity of the same kind."
    fail
  fi

  # INITIALIZATION MODE IS NOT ERROR-CONTROLLED. A relaxation makes no
  # statement about a trajectory, so it has no integration error to bound
  # and takes no comparison pass.
  d=$(RUN_MODE=init STEPS=12 run_case init_no_error_control)
  nest=$(grep -c "integration-error estimate:" "$d/run.log" || true)
  npas=$(grep -c "err-pass:" "$d/run.log" || true)
  if [ "$nest" = "0" ] && [ "$npas" = "0" ]; then
    say PASS no_error_control_in_initialization_mode "${nest}_${npas}" 0_0 0
  else
    say FAIL no_error_control_in_initialization_mode "${nest}_${npas}" 0_0 0
    fail
  fi
else
  echo "SKIP the whole-binary rows: set EXHALE_ATTEMPTED_STEP_EXE"
fi

if [ "$status" -ne 0 ]; then
  echo "attempted_step: FAILED"
else
  echo "attempted_step: PASSED"
fi
exit "$status"

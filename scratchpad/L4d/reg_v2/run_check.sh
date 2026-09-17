#!/bin/bash
# Regression check for no-physics-change refactor steps.
#
# Usage:
#   ./run_check.sh golden [case...]   # snapshot current <case>/output as reference
#   ./run_check.sh check  [case...]   # rebuild, re-run each case single-thread,
#                                     #   bitwise-compare vs its golden
#
# Environment:
#   REGRESSION_REL_TOL=<x>      relative tolerance of the "within tolerance"
#                               verdict (default 1e-3; 0 = strict bitwise).
#   REGRESSION_GOLDEN_DIR=<dir> the directory "check" reads its reference
#                               files from (default $HERE/golden). A scratch
#                               baseline lives in baseline_post170_<date>/ and
#                               is read through this variable while golden/ is
#                               deliberately left stale (development plan
#                               rev 3 section 3.2). "golden" mode is NOT
#                               affected: it always writes to golden/, so a
#                               snapshot can never land in a scratch baseline
#                               by an environment variable left set.
#   REGRESSION_EXE=<path>       run the matrix with that binary and skip the
#                               rebuild (another compiler, an older build).
#
# Single-threaded (OMP_NUM_THREADS=1) => fully deterministic => bitwise compare.
# A no-physics-change refactor must produce byte-identical Hydro_ioniz.txt and
# Ion_species.txt for EVERY case in the feature matrix.
#
# The cases are independent -- each runs in its own directory on its own copy
# of the binary -- so "check" starts them all at once and waits; the bitwise
# guarantee is a property of the single-threaded run, not of the order the
# cases start in. Each case writes its report to <case>/check.log, and the
# reports are printed in list order once every case has finished.
#
# A case directory holds its own input.inp (plus metals.inp / base.inp / ...),
# an output/ dir, and optionally a one-line "maxsteps" file for cases that are
# relaxation snapshots rather than converged solutions.
#
# Directories here that are NOT cases:
#   golden/, golden_*/   the reference snapshots.  golden/ holds exactly the
#                        DEFAULT_CASES below (eleven since 2026-09-04); the
#                        dated golden_block*/ directories are earlier snapshots
#                        and hold the subset of the matrix that existed when
#                        they were taken.
#   parse_golden*/       the parser corpus, refreshed by run_parse_corpus.sh.
#   windae_oracle/       a campaign directory: it holds its own sub-runs and
#                        cannot be passed as an argument to this script.
#                        valve_sens/ was one too, and moved to _quarantined/ on
#                        2026-09-04: it scanned `Valve eps`, which section 152
#                        retires, so there is no epsilon left to scan.
#   _quarantined/        cases withdrawn from the matrix, kept for the record
#                        with a NOTE.txt saying why.  Do not pass these as
#                        arguments.  _quarantined/armD_D1 (moved 2026-09-03)
#                        cannot start at all: it is refused by the P35 check
#                        in input_read.f90 ("Molecular base: True" with
#                        "Molecular chemistry: False").
#
# Deleted 2026-09-03 (user decision), each byte-identical to the case named
# beside it -- same input.inp and metals.inp, only the stored outputs differed:
#   crit_warm          -> use wasp_hybrid_finish
#   ptc_wasp           -> use resid_on
#   solver_newton_cold -> use wasp_full_newton
# None of the three was ever in DEFAULT_CASES or in any golden/ directory.
#
# Matrix cases (add more dirs as the matrix grows):
#   wasp_full         He23S on,  metals on   (exercises HeITR + metal paths)
#   wasp_he23off      He23S off, metals on   (exercises the HeITR-off branch)
#   mol_base_handoff  molecular chemistry + base.inp handoff on the Tier-2 gate
#                     hot Uranus AT 1 MICROBAR (the Koskinen 2022 handoff level;
#                     base.inp p_base fixes the level and input.inp carries no
#                     density key), q_H2_base = 0.84 from that same paper:
#                     q_H2_base > 0 drives the photochemical branch of the
#                     molecular base particle count (comp_ntot_bc);
#                     12000-step snapshot, the Tier-2 gate convention
#   mol_metals        the same hot-Uranus Tier-2 gate + solar 7-element trace
#                     metals: guards the molecular+metals coupled system
#                     (metal electrons dominate the shielded molecular base);
#                     12000-step snapshot
#   mol_lyman_werner  mol_base_handoff plus a Lyman-Werner band flux
#                     (Stellar LW flux 343 erg/cm2/s): guards the H2
#                     photodissociation path -- star-ward H2 column,
#                     Draine & Bertoldi self-shielding, the H2-budget loss
#                     row and the 0.4 eV heating; 12000-step snapshot
#   mol_diffusion     mol_base_handoff plus binary H/He element diffusion
#                     (He_diffusion: True, He_Kzz 1e9): the widest single
#                     guard on binary_element_diffusion -- molecular-carrier
#                     closure, stage-resolved friction pairs, projection back
#                     onto f_sp, and the Coulomb friction of the ionized
#                     region, all in one run; 12000-step snapshot
#   mol_ir_bands      mol_base_handoff plus "Base IR field: True" and
#                     "Molecular IR bands: True": guards the H2 / H2O / CO
#                     infrared cooling channels of the molecular layer;
#                     12000-step snapshot
#   mol_sec_ion       mol_base_handoff plus "Secondary_ionization: Immediate":
#                     the ONLY case in which a photoelectron is partitioned in
#                     molecular gas at all.  Secondary_ionization defaults to
#                     STAGED -- the coupling switches on only after the wind
#                     first converges without it -- and every molecular gate is
#                     a fixed-step relaxation snapshot that never converges, so
#                     in all of them sec_ion_active never flips and the whole
#                     partition path is bypassed.  This case guards the H2
#                     secondary-ionization target, the H2 dissociation heat
#                     return, and the Dalgarno heating fraction in the gas
#                     they were written for; 12000-step snapshot
#   mol_carrier       mol_base_handoff plus "Molecular carrier transport: True":
#                     the ONLY case that runs the H2 carrier transport
#                     operator (diffusive_photochemistry with no oxygen
#                     chemistry).  It guards the second-order limited
#                     advection and its deferred correction, the base inflow
#                     condition taken from the wind's mass flux rather than
#                     from the base cell's velocity, and the hydrogen headroom
#                     the chemistry rows and the element limiter now share;
#                     12000-step snapshot
#   lower_profile     the lower atmosphere handed over as a PROFILE instead of
#                     the scalars of base.inp ("Lower atmosphere profile:",
#                     docs/phase_e_flux_closure_design.md): pins the profile
#                     reader, the matching-level base state (T0, R0, p_base,
#                     q_H2, He/H), the elemental reservoirs the profile
#                     carries (solar C/N/O, so the case is metals-on with no
#                     metals.inp), K_zz(p) interpolated onto the grid instead
#                     of the scalar He_Kzz, the accepting branch of the
#                     profile / base.inp pair enforcement (a provenance-only
#                     base.inp whose solution_id matches), and the elemental
#                     flux window statistics; 12000-step snapshot
#   wasp_full_newton  the ONLY default case that reaches the steady solver.
#                     Added 2026-09-04 because a named-case audit found that
#                     none of the ten cases above ever runs one: the two
#                     atomic cases have no "Solver:" key and stop on the flux
#                     gate, and the eight molecular/profile cases stop at
#                     EXHALE_MAXSTEPS=12000 before the JFNK hand-off can fire,
#                     so the whole steady path (solve_steady_jfnk /
#                     solve_steady_ptc, and the acceptance measures of
#                     Update_EXHALE.md sections 126 to 155) was rebuilt
#                     repeatedly with no golden guarding it.  This case
#                     converges (info=0), so it guards the solver AND is a
#                     byte-identity test of it.

# Named cases OUTSIDE the matrix (pass them as arguments; they are not in
# DEFAULT_CASES):
#   roundtrip         the restart round trip, rebuilt 2026-09-03 from a SHORT
#                     deterministic run (defect 12.4 of
#                     docs/open_defects_20260903.md).  The version that carried
#                     this name could not run at all -- its
#                     output/Hydro_ioniz_IC.txt was zero bytes, there was no
#                     Ion_species_IC.txt, and "Load IC? True" made it abort in
#                     load_IC -- because it had been built from a converged
#                     product that no longer exists.  It now builds its own
#                     state: 40 deterministic steps of the hot-Uranus molecular
#                     gate, so `check roundtrip` bitwise-compares it like any
#                     other case.  The round trip itself is the second half and
#                     is run separately, after that:
#                         ( cd roundtrip && ./roundtrip_check.sh )
#                     which hands the stage-A outputs back as the _IC pair,
#                     reloads with EXHALE_DUMP_IC=1, and requires every species
#                     column (matched BY LABEL), rho, v, p and the
#                     '# coupling:' header to round-trip to 1e-12.  It restores
#                     input.inp and the stage-A outputs on exit, so `golden
#                     roundtrip` still snapshots the right state.  It is out of
#                     DEFAULT_CASES because its second half is not a bitwise
#                     comparison and needs its own invocation; the retired
#                     contents, check_roundtrip.py among them, are kept in
#                     roundtrip/.superseded_20260903/.

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
# Repo root is two levels up (backup/regression -> backup -> EXHALE), where
# the Makefile and EXHALE.x live.
ROOT="$(cd "$HERE/../.." && pwd)"
# The files compared. The two *_adv.txt profiles were added 2026-09-05
# (item 10.3 of docs/development_plan_20260905_rev3.md): post_process_adv.f90
# writes the profiles the transit tool and the analysis notebooks read, and
# nothing in the matrix compared them, so the whole post-processing pass was
# outside the regression. They carry the same '# columns' schema header as
# the two files beside them, so the same header-skipping comparison applies.
#
# A reference directory written before that date has no *_adv.txt in it.
# Such a file is reported as MISS and is NOT a failure: the reference of
# record (golden/) is refreshed once at the end of Phase 1 and picks the two
# files up then, because "golden" mode snapshots this same list. The count of
# missing references is printed with the verdict, so it is never silent.
FILES="output/Hydro_ioniz.txt output/Ion_species.txt \
       output/Hydro_ioniz_adv.txt output/Ion_species_adv.txt"
# Relative tolerance below which a difference from the golden counts as
# identical (see the comparison block in check_case). Override for one run
# with REGRESSION_REL_TOL=... ; 1e-3 is the 0.1 per cent the standing rule
# names. Set it to 0 to demand strict bitwise agreement.
REL_TOL="${REGRESSION_REL_TOL:-1e-3}"
# Where "check" reads its reference files. Default golden/, the reference of
# record; a scratch baseline is read by setting REGRESSION_GOLDEN_DIR.
GOLDEN_DIR="${REGRESSION_GOLDEN_DIR:-$HERE/golden}"

# ONE RUN AT A TIME, and the reason is a false verdict that actually happened.
# Every case runs IN its own case directory and writes its outputs and its
# check.log there, so two invocations of this script share those files: the
# second run's executable overwrites the first's, the two write the same
# output/*.txt while the other is comparing them, and their check.log lines
# interleave inside a single line. On 2026-09-05 an overlapping strict-mode
# run (REGRESSION_REL_TOL=0) left a "FAIL ... EXCEEDS 0.0e+00" line in a log
# that the tolerance run then counted, and the matrix reported FAIL on a case
# that had passed. Refuse instead: a stale lock names the process still
# holding it, and if that process is gone the lock can be removed by hand.
LOCK="$HERE/.run_check.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  holder=$(cat "$LOCK/pid" 2>/dev/null || echo unknown)
  if [ "$holder" != unknown ] && kill -0 "$holder" 2>/dev/null; then
    echo "refusing to start: $0 is already running as PID $holder." >&2
    echo "  Two runs share the case directories and corrupt each other's" >&2
    echo "  outputs and verdicts. Wait for it, or stop it by PID." >&2
  else
    echo "refusing to start: a stale lock is present at $LOCK" >&2
    echo "  (it names PID $holder, which is not running)." >&2
    echo "  Check that no run is in flight, then: rmdir $LOCK/pid $LOCK 2>/dev/null; rm -rf $LOCK" >&2
  fi
  exit 3
fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK"' EXIT INT TERM

DEFAULT_CASES="wasp_full wasp_he23off wasp_full_newton mol_base_handoff mol_metals mol_lyman_werner mol_diffusion mol_ir_bands mol_sec_ion mol_carrier lower_profile hydrostatic_column oxygen_chemistry hp_zero_seed hp_trace_seed hp_front"

# Cases built in Phase 0 (2026-09-05). They joined DEFAULT_CASES with the
# golden refresh at the end of Phase 1 batch 2c (docs/Update_EXHALE.md
# section 6); `phase0-cases` still lists them on their own:
#
#   ./run_check.sh check $(./run_check.sh phase0-cases)
#
#   hydrostatic_column  a mechanical column: no irradiation, no chemistry,
#                       no conduction, no viscosity. It exists so that the
#                       discrete hydrostatic residual and the spurious
#                       velocity of item D2 (Phase 4) have a case to be
#                       measured on; the measurement itself is Phase 4 work
#                       and this directory only pins the configuration and
#                       its output.
#   oxygen_chemistry    the configuration of examples/18_oxygen_chemistry
#                       (HD 209458 b, molecular + oxygen chemistry, the five
#                       FUV band fluxes), capped to a step count that runs in
#                       minutes. Item 10.3: the whole oxygen path -- the
#                       oxygen rate set, water photolysis, the OH/H2O
#                       carriers -- was outside the regression.
#   hp_zero_seed        transported ionization state ("Ionization transport:
#   hp_trace_seed       True"), restarted from the mol_base_handoff state
#   hp_front            with the H II column set to zero, to 1e-12 of the H
#                       nuclei, and left as written. The seed reaches the
#                       solver only through a restart: the carrier constraint
#                       is imposed on the sweep only when the background is
#                       ready or an IC is loaded (ionization_equilibrium.f90
#                       1203-1211), so a cold start would initialize the
#                       proton from the local root and the seed would be
#                       nothing. 100-step snapshots.
PHASE0_CASES="hydrostatic_column oxygen_chemistry hp_zero_seed hp_trace_seed hp_front"

mode="$1"; shift || true
CASES="${*:-$DEFAULT_CASES}"

# Run one case and compare it against its golden. Writes the report to stdout
# and exits 1 on any FAIL, so the caller can read the verdict off the exit
# status as well as the text.
check_case() {
  local c="$1" cap="" crc=0 f b
  # A case that is a relaxation snapshot rather than a converged solution
  # pins its step count in <case>/maxsteps (read by EXHALE_MAXSTEPS).
  # No such file => no cap, i.e. the run the case had before.
  [ -f "$HERE/$c/maxsteps" ] && cap="EXHALE_MAXSTEPS=$(cat "$HERE/$c/maxsteps")"
  echo "[$c] running (OMP_NUM_THREADS=1 $cap) ..."
  # A case that cannot be given the binary must not go on to compare the
  # outputs of the PREVIOUS run sitting in its directory against the
  # reference and report PASS. The copy fails on a missing case directory, a
  # read-only tree, or a full disk, and none of those is visible afterwards.
  if ! cp "${REGRESSION_EXE:-$ROOT/EXHALE.x}" "$HERE/$c/EXHALE.x"; then
    echo "       FAIL cannot copy ${REGRESSION_EXE:-$ROOT/EXHALE.x} into $HERE/$c/"
    echo "            (missing case directory, read-only tree, or no space)"
    return 1
  fi
  # A RESTART CASE CARRIES ITS OWN IC PAIR, AND output/ IS CLEARED FIRST.
  # load_IC reads output/Hydro_ioniz_IC.txt and output/Ion_species_IC.txt
  # (load_IC.f90 119, 161, 181), and the rm below clears output/*.txt before
  # every run, so a case with "Load IC? True" that kept its pair in output/
  # would run exactly once and fail for the rest of the directory's life.
  # Such a case keeps the pair in <case>/IC/ instead, and it is restored here.
  # A case with no IC/ directory is unaffected.
  ( cd "$HERE/$c" && rm -f output/*.txt \
    && if [ -d IC ]; then cp IC/*_IC.txt output/; fi \
    && env $cap OMP_NUM_THREADS=1 ./EXHALE.x > run.log 2>&1 )
  echo "       $(grep -E 'final:' "$HERE/$c/run.log" | tail -n1)"
  for f in $FILES; do
    b=$(basename "$f")
    # A reference directory written before a file entered $FILES has no copy
    # of it (the *_adv.txt pair, added 2026-09-05, is the case that made this
    # reachable). Say so and move on: comparing against nothing is not a
    # verdict, and failing on it would make every `check` against the
    # unrefreshed golden/ red for a reason that is not a code change.
    if [ ! -f "$GOLDEN_DIR/$c/$b" ]; then
      echo "       MISS $b (no reference in $GOLDEN_DIR/$c/)"
      continue
    fi
    # Compare numeric content only: '#' schema-header lines are excluded,
    # so adding/changing headers never breaks the bitwise data guarantee.
    if cmp -s <(grep -v '^ *#' "$HERE/$c/$f") <(grep -v '^ *#' "$GOLDEN_DIR/$c/$b"); then
      echo "       PASS $b (data identical)"
    else
      # NOT IDENTICAL: ask HOW FAR, and let the answer decide.
      #
      # Bitwise is still the first question, because a refactor with no
      # intended physics change must not move a digit. But some changes move
      # results by an amount no physics can notice -- a compiler or LAPACK
      # swap, a constant that two files had defined 0.05 per cent apart, a
      # thread reordering the summation of a reduction -- and spending a
      # golden refresh on those destroys the reference the whole matrix is
      # measured against. The standing rule (2026-09-05) is that a
      # numerical or technical method difference under REL_TOL counts as
      # identical and the golden stays. So the verdict is measured, never
      # assumed: compare_within_tolerance.py reports the worst RELATIVE
      # difference and the row, column and pair of values that carry it, and
      # every within-tolerance PASS below is printed with that number beside
      # it. A real physics change lands far above REL_TOL and still fails.
      tolmsg=$("$HERE/compare_within_tolerance.py" \
                 "$HERE/$c/$f" "$GOLDEN_DIR/$c/$b" "$REL_TOL" 2>&1)
      tolrc=$?
      if [ $tolrc -eq 0 ]; then
        echo "       PASS $b ($tolmsg)"
      else
        echo "       FAIL $b ($tolmsg)"
        crc=1
      fi
    fi
  done
  return $crc
}

case "$mode" in
  golden)
    # NB: golden mode does NOT re-run the cases -- it snapshots whatever
    # output/ currently sits in each case directory. After a code change,
    # run "check" first (regenerates the outputs with the current binary),
    # then "golden", then "check" again to verify. Snapshotting stale
    # outputs silently baselines the PREVIOUS binary (bitten 2026-08-13).
    for c in $CASES; do
      mkdir -p "$HERE/golden/$c"
      for f in $FILES; do cp "$HERE/$c/$f" "$HERE/golden/$c/$(basename $f)"; done
      cp "$HERE/$c/run.log" "$HERE/golden/$c/run.log" 2>/dev/null || true
      echo "golden[$c]: $(grep -E 'final:' "$HERE/$c/run.log" | tail -n1)"
    done
    ;;
  check)
    # REGRESSION_EXE=<path>: run the matrix with THAT binary instead of the
    # one the sources make (another compiler, an older build), and skip
    # the rebuild and the up-to-date test, which are about $ROOT/EXHALE.x.
    # The verdict then certifies that binary, and the log says which.
    if [ -n "${REGRESSION_EXE:-}" ]; then
      echo "[build] skipped: REGRESSION_EXE=$REGRESSION_EXE ($(md5sum "$REGRESSION_EXE" | cut -c1-12))"
    else
      echo "[build] rebuilding EXHALE.x ..."
      # THE REBUILD IS FATAL, AND ITS OUTPUT IS KEPT.
      #
      # This line used to be `( cd "$ROOT" && make >/dev/null 2>&1 ) && echo
      # "  build OK"`, which is two defects in one: `set -e` does not act on a
      # command that is the left operand of `&&`, so a failed compile did not
      # stop the script, and the compiler's diagnostics went to /dev/null, so
      # nothing said WHY. The matrix then ran on whatever EXHALE.x was already
      # there. (The `make -q` guard below caught that particular case, because
      # an unbuilt object leaves work to do -- but only after the reason had
      # been thrown away, and it says nothing about a build that fails while
      # leaving every target up to date.) `if ! ...; then` puts the command in
      # a condition, where its status is the thing being tested and cannot be
      # swallowed.
      BUILDLOG="$HERE/.build.log"
      if ! ( cd "$ROOT" && make ) > "$BUILDLOG" 2>&1; then
        echo "  BUILD FAILED: 'make' in $ROOT returned nonzero."
        echo "  Compiler output ($BUILDLOG), last 40 lines:"
        tail -n 40 "$BUILDLOG" | sed 's/^/    /'
        echo "  Refusing to run the matrix: a verdict from an old binary"
        echo "  certifies code that is not in the tree."
        exit 2
      fi
      echo "  build OK (compiler output kept in $BUILDLOG)"
      # THE BINARY THE MATRIX RUNS MUST BE THE ONE THE SOURCES MAKE.
      #
      # `make` above can succeed and still leave EXHALE.x older than a source
      # whose rebuild the dependency graph missed, and a bitwise "identical to
      # golden" verdict taken from such a binary certifies the wrong code. On
      # 2026-09-03 a check was run on an incremental build whose provenance
      # could no longer be reconstructed, precisely because nothing here ever
      # asked. `make -q` is exactly that question: it exits 0 when every target
      # is up to date and 1 when something still needs making, and it builds
      # nothing. This is a hard stop, not a warning: a stale binary makes every
      # line below meaningless.
      if ! ( cd "$ROOT" && make -q >/dev/null 2>&1 ); then
        echo "  BUILD NOT UP TO DATE: 'make -q' still reports work to do after"
        echo "  a successful 'make'. The dependency graph missed something;"
        echo "  run 'make distclean && make' and investigate before trusting"
        echo "  any result. Refusing to run the matrix."
        exit 2
      fi
      echo "  build up to date (make -q)"
    fi
    # Every case at once; each report goes to its own file so the outputs do
    # not interleave, and "set -e" is lifted inside the subshell so a FAIL is
    # reported rather than aborting the case mid-comparison.
    for c in $CASES; do
      ( set +e; check_case "$c" ) > "$HERE/$c/check.log" 2>&1 &
    done
    wait
    # THE SUMMARY IS BUILT FROM THE CASES OF THIS RUN ONLY.
    #
    # The tolerance question used to be asked of "$HERE"/*/check.log, i.e. of
    # every case directory that has ever been checked, including the ones this
    # invocation did not run. A `check wasp_full` that was byte-identical then
    # printed "within $REL_TOL relative -- see the max rel values above" on the
    # strength of a months-old log of another case, and the max rel values it
    # points at are not above. The three counters are accumulated in the loop
    # that reads the logs, so they can describe nothing else.
    rc=0; ntol=0; nmiss=0; ncmp=0
    for c in $CASES; do
      cat "$HERE/$c/check.log"
      grep -q '^       FAIL' "$HERE/$c/check.log" && rc=1
      if grep -q '^       PASS .*max rel' "$HERE/$c/check.log"; then ntol=1; fi
      m=$(grep -c '^       MISS' "$HERE/$c/check.log" || true)
      k=$(grep -c -E '^       (PASS|FAIL)' "$HERE/$c/check.log" || true)
      nmiss=$(( nmiss + m ))
      ncmp=$(( ncmp + k ))
    done
    if [ $nmiss -gt 0 ]; then
      echo "==> $nmiss file(s) were NOT compared: no reference in $GOLDEN_DIR"
      echo "    (a reference directory older than an entry of FILES; the"
      echo "     verdict below is about the files that do have one)"
    fi
    # A run in which nothing was compared has no verdict to give, and "all
    # cases byte-identical" is exactly the wrong thing to print: it is what a
    # new case reports on its first invocation, before the reference exists.
    if [ $ncmp -eq 0 ]; then
      echo "==> NO VERDICT: not one file had a reference in $GOLDEN_DIR."
      echo "    The cases ran and their outputs are in <case>/output/;"
      echo "    nothing was compared against anything."
      exit 0
    fi
    if [ $rc -eq 0 ]; then
      if [ $ntol -eq 1 ]; then
        echo "==> REGRESSION PASS (identical, or within $REL_TOL relative -- see the max rel values above)"
      else
        echo "==> REGRESSION PASS (all cases byte-identical)"
      fi
    else
      echo "==> REGRESSION FAIL (a case exceeds $REL_TOL relative)"
    fi
    exit $rc
    ;;
  phase0-cases)
    # The Phase 0 cases as one word list, so that a `check` over them does not
    # have to repeat the names (see the comment beside PHASE0_CASES).
    echo "$PHASE0_CASES" ;;
  *)
    echo "usage: $0 {golden|check|phase0-cases} [case...]"; exit 2 ;;
esac

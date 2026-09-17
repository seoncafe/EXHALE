#!/bin/bash
# Solve one prescribed-composition LHS 1140 b case and synthesize its transit
# spectrum, by the recipe MODELS.md section 6 adopted on 2026-09-13.
#
#     models/run_case.sh <group>/<case>
#
# Four passes, in the order the result depends on:
#
#   0. the seed.  An atomic case starts from a solved state of the same
#      physics, chosen by models/pick_seed.py: a CERTIFIED case of models/
#      where one exists, an archived state of archive_20260830/ otherwise.
#      It is interpolated onto the current cell centers by
#      src/utils/map_state_to_grid.py (the 2026-09-10 grid change and the
#      +2.26 percent Jupiter radius put every archived state on centers this
#      binary refuses).  The grid is the same in R_p for every case of the
#      tree, so one file, models/current_grid_Hydro_ioniz.txt, names it.
#      THE RADIAL SCALE R0 IS NOT SHARED. A case with a lower-atmosphere
#      profile builds its grid on the profile's radius at the matching level
#      times R_J, not on "Planet radius [R_J]", and load_IC compares the
#      "grid" metadata field as text, so such a case gets its own target
#      header: target_grid_Hydro_ioniz.txt, the shared file with that one
#      field replaced, written by src/utils/profile_match_level.py, which is
#      the same routine the elemental-flux closure seeds an iteration with.
#      A seed that states a reservoir He/H -- every state this code wrote --
#      and whose reservoir is not this case's is carried onto it by
#      --reservoir He/H, which multiplies every helium column by the ratio
#      of the two: load_IC refuses a state whose reservoir field is not the
#      run's, and the base rows of a diffused profile then arrive at this
#      case's composition with the profile's shape kept.  An archived state
#      carries no reservoir field and load_IC rescales it by its own rule.
#      A MOLECULAR case has no archived state of its physics and cannot be
#      given one: "mol" and "carrier" decide how many unknowns a state has,
#      so load_IC refuses an atomic file for it by name.  Its seed is built
#      instead, by the binary, out of the CERTIFIED atomic case of the same
#      name and He/H: EXHALE_MOLECULAR_SEED reads that pair, forms the
#      molecular state, writes output/*_IC.txt and stops
#      (docs/input_schema.md appendix D.3,
#      docs/lhs1140b_stationary_L7_20260913.md).
#
#   1. the wind.  "Restart intent: stationary" measures the state as loaded,
#      takes no time step, and enters the partitioned stationary solve: the
#      hydrodynamic rows by JFNK and the element partition by its relaxation,
#      alternating until every active equation is within its own tolerance.
#      What makes it reach that state on this planet is "Well balanced:
#      True": at the base Mach number of 3e-7 the plain HLLC contact speed is
#      set by the reconstruction's pressure truncation error rather than by
#      the flow, and the energy row then reads its own scale
#      (docs/lhs1140b_stationary_L5c_20260913.md).  A cold molecular case
#      marches instead, two-stage PLM then WENO3, with the JFNK hand-off.
#      The verdict is the solve's own info and the certification block.
#
#   2. the advection-corrected profiles.  EXHALE reads back its own solved
#      state ("Load IC? True", "Do only PP: True") and writes the *_adv.txt
#      files and the steady-state Mdot the transit synthesis consumes.  That
#      pass takes one time step before it stops, which would move a certified
#      state, so it runs at CFL 1e-12: measured, the state it writes then
#      differs from the state it read by 4e-14 in rho and 2e-12 in v, against
#      1.2e-6 and 0.59 at the default CFL.  input.inp is put back the way it
#      was afterwards, so the file in the case directory is always the file
#      the solution was started from.
#
#   3. the transit spectrum, through the WINERED HIRES-Y kernel of the
#      measurement (R = 68000, LHS1140b/winered_hires_y.sh).
#
#   4. the record.  models/write_reproduce.py writes REPRODUCE.md in the case
#      directory: the physics, the binary and its md5, the seed and how it was
#      mapped, every command with its environment, the outer passes and the
#      certification as run.log states them, the line measures, and the clock.
#      It is written by the run it describes, and also when a pass fails after
#      the wind, so a case that did not finish still says how far it got.
#
# Environment:
#   OMP_NUM_THREADS  threads for the wind pass (default 8)
#   EXHALE_BIN       the binary (default <repo>/EXHALE.x)
#   SEED             an output/ directory to map in place of pick_seed's choice.
#                    For a MOLECULAR case it names a molecular solution to
#                    continue from and REPLACES the atomic->molecular
#                    conversion; the mapper matches the "# grid" header and
#                    carries the reservoir with --reservoir He/H.
#   NOSEED=1         skip pass 0 and start from whatever output/ holds
#   FORCE=1          re-run a case that already carries tpm_He10830_metrics.txt
#   DTAU0_CONTINUATION  the pseudo-time start of the second solve when the
#                    first ends refused (default 1.0e8; L4e)
#   SEED_ATTEMPTS    how many of pick_seed's ranked seeds are tried before the
#                    case is given up (default 3); each failed attempt is kept
#                    in attempt_<n>/
#   ATOMIC_SEED      for a molecular case, the atomic case directory the seed
#                    is built from (default: this case's name with the leading
#                    "molecular" replaced by "atomic")
#   SEED_X2          which extension of the base partition the molecular seed
#                    carries: "local" (default) gives every cell the H2
#                    content that HOLDS there -- in EVERY cell, the base layer
#                    included, the smaller of the thermochemical fit q_H2(p,T)
#                    and the root of that cell's own H2 carrier row, with the
#                    photodissociation, the photoionization and the ion
#                    channels in it (item L7e; the fit alone asserts "every H
#                    nucleus is in H2" over both ends of the column and in the
#                    irradiated outer wind that is 80 to 360 times the root).
#                    The root is SOLVED for and not divided out: the row is
#                    quadratic in n(H2), because the three-body association
#                    carries it and atomic hydrogen is closed out of the
#                    element budget (item L7f).
#                    "handoff" carries the base x2 to the outer boundary, a
#                    number in [0,1] is that fraction everywhere.  MEASURED on HeH2.13
#                    (docs/lhs1140b_stationary_L7c_20260914.md section 3.2):
#                    at the handoff's own x2 = 0.99998 the column has no
#                    neutral hydrogen left to photoionize above the base, the
#                    energy row of the seed is its radiative imbalance, and
#                    the solve makes no outer pass; "local" keeps the wind's
#                    own H2 front and is the default for that reason.
#   SEED_INVARIANT   p (default) or T, the primitive the conversion keeps
#
# Exit status 0 when every pass finished, 1 otherwise.

set -u

MODELS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LHS=$(dirname "$MODELS")
EX=$(dirname "$LHS")

if [ $# -ne 1 ]; then
   echo "usage: $(basename "$0") <group>/<case>" >&2
   exit 1
fi
CASE=${1%/}
DIR="$MODELS/$CASE"
: "${OMP_NUM_THREADS:=8}"
: "${EXHALE_BIN:=$EX/EXHALE.x}"
: "${FORCE:=0}"
: "${NOSEED:=0}"
: "${SEED:=}"
export OMP_NUM_THREADS
# The pseudo-time the stationary solve starts its first outer pass at; 1.0 is
# the value the HD 209458 b element-diffusion control and every LHS 1140 b
# solve of 2026-09-13 were made at, and is the default here.
#
# WHEN TO RAISE IT (L14). On the He/H = 2.13 column of LHS 1140 b below 0.10
# of the GJ 1132 spectrum the dtau0 = 1 solve falls into the ramp collapse of
# plan item L4h -- MEASURED at XUV steps of a factor 2 and a factor 1.43
# alike, lam 9.5e-7 and dtau 6e-4 with the residual unmoved over 130 of the
# 3000 iterations it is allowed -- while the SAME mapped seed solved at
# dtau0 = 1e8 reaches its root in about 40 iterations. Where that is the
# case, state EXHALE_PTC_DTAU0=1.0e8 in the environment and the first solve
# is taken there; the continuation below is then not needed and does not run.
: "${EXHALE_PTC_DTAU0:=1.0}"
export EXHALE_PTC_DTAU0
# HOW MANY OUTER PASSES THE ALTERNATION MAY SPEND (L14). The element
# relaxation removes its distance to its own fixed point by a factor
# 1 - omega = 0.5 per pass, so a seed one XUV step away needs as many passes
# as that distance has decades: MEASURED on the certified
# atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13, 4.89e-1 at the first pass to
# 2.90e-4 at the twelfth and certified at the thirteenth, and the lower-XUV
# cases enter further out still (6.62e-1 at 0.10). The binary's own default
# is 20, which leaves no room for the first passes of a far seed, on which
# the gated row rises while the composition is still travelling: MEASURED on
# atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13, whose row peaks at pass 5 and is
# ACCEPTED at pass 30 exactly.
: "${EXHALE_OUTER_PASSES:=40}"
export EXHALE_OUTER_PASSES
# One paragraph for this case's REPRODUCE.md, stating what this run had to be
# told and why. It is used where a run replaces an earlier attempt whose logs
# are no longer in the case directory (the logs go to models/.stopped/).
: "${RUN_NOTE:=}"

GRID="$MODELS/current_grid_Hydro_ioniz.txt"

STARTED=$(date '+%Y-%m-%d %H:%M:%S')
T_START=$(date +%s)
COMMANDS=""
record () {
   # The record of this run, written by the run itself.  It is not allowed to
   # fail the case: a missing record is reported and nothing else.
   [ -n "${WROTE_RECORD:-}" ] && return 0
   # Nothing to record before the solver has written anything.
   [ -f "$DIR/run.log" ] || return 0
   WROTE_RECORD=1
   local wall=$(( $(date +%s) - T_START ))
   python3 "$MODELS/write_reproduce.py" "$DIR" --kind case \
      --case-name "$CASE" --binary "$EXHALE_BIN" \
      --threads "$OMP_NUM_THREADS" \
      --seed-line "$SEED_NOTE" --seed-command "$SEED_COMMAND" \
      --seed-rescale "$SEED_RESCALE" --note "$RUN_NOTE" \
      --grid-target "${GRID_TARGET:-}" --grid-r0 "${GRID_R0:-}" \
      --commands "$COMMANDS" \
      --started "$STARTED" --finished "$(date '+%Y-%m-%d %H:%M:%S')" \
      --wall "$(printf '%dm%02ds' $((wall/60)) $((wall%60)))" \
      > /dev/null || echo "[$CASE] note: REPRODUCE.md was not written"
}

fail () { record; echo "[$CASE] FAILED $*"; exit 1; }

[ -d "$DIR" ]            || fail "no such case directory $DIR"
[ -f "$DIR/input.inp" ]  || fail "$DIR has no input.inp (a closure rung uses run_closure.sh)"
[ -x "$EXHALE_BIN" ]     || fail "no executable $EXHALE_BIN"

if [ -f "$DIR/tpm_He10830_metrics.txt" ] && [ "$FORCE" != 1 ]; then
   echo "[$CASE] SKIP (already synthesized; FORCE=1 to redo)"
   exit 0
fi

cd "$DIR" || fail "cannot enter $DIR"
mkdir -p output

# ---------------------------------------------------------------- the seed --
# A case that states "Load IC? True" is solved from a seed; one that states
# False marches from the code's own initial condition and wants none.
LOADS=$(sed -n 's/^Load IC? *//p' input.inp | tail -n 1)
SEED_NOTE="cold start"
SEED_COMMAND=""
SEED_RESCALE=""
: "${SEED_GIVEN:=$([ -n "$SEED" ] && [ -z "${ATTEMPT:-}" ] && echo 1 || echo 0)}"
MOLECULAR=0
grep -qi '^Molecular chemistry: *True' input.inp && MOLECULAR=1
# A MOLECULAR CASE TAKES A STATED SEED LIKE ANY OTHER CASE.  With SEED given
# the state named is a MOLECULAR solution to be interpolated onto this grid
# (the branch below), not an atomic one to be converted: a continuation from
# the nearest solved molecular state is the campaign's own recipe, and for a
# base layer whose composition is set by eddy transport rather than by local
# chemistry it is the only seed that starts near the fixed point (item L7f
# section 5.2 measures what the local-chemistry seed costs there).  Without
# SEED the case is built from the certified atomic solution of the same name
# and He/H, which is what a molecular case with no molecular state to
# continue from has.
if [ "$MOLECULAR" = 1 ] && [ -z "$SEED" ] && [ "$NOSEED" != 1 ]; then

   # ------------------------------------------------ the molecular seed --
   # The certified atomic solution of the same name and He/H, converted by
   # the binary itself.  The conversion calls the production census,
   # equation of state, boundary and metadata routines and writes an
   # UNCERTIFIED initialization pair (certified=F cert_reason=molecular_seed
   # mode=init) with this case's own option, species and reservoir metadata
   # and a "# molecular_partition:" line in both halves.
   #
   # MEASURED 2026-09-13/14 (docs/lhs1140b_stationary_L7_20260913.md,
   # docs/lhs1140b_stationary_L7c_20260914.md): the conversion is exact (the
   # H-nucleus, helium and mass budgets to 3.5e-16, the equation-of-state
   # round trip to 3.4e-16) and the pair loads with no option permitted to
   # differ.  Which extension of the base partition it carries is SEED_X2
   # above.
   : "${ATOMIC_SEED:=$MODELS/${CASE/#molecular/atomic}}"
   : "${SEED_X2:=local}"
   [ -d "$ATOMIC_SEED/output" ] || \
      fail "no atomic case $ATOMIC_SEED to build the molecular seed from"
   grep -q 'certified=T' "$ATOMIC_SEED/output/Hydro_ioniz_IC.txt" 2>/dev/null || \
      fail "the atomic state in $ATOMIC_SEED/output is not certified; a seed is only as good as the wind it is made of"
   # "handoff" is the conversion's own default and is asked for by leaving
   # the variable out of the environment, so it is a name here and not a value.
   X2_ENV=""
   [ "$SEED_X2" = handoff ] || X2_ENV="EXHALE_MOLECULAR_SEED_X2=$SEED_X2"
   INV_ENV=""
   [ -z "${SEED_INVARIANT:-}" ] || INV_ENV="EXHALE_MOLECULAR_SEED_INVARIANT=$SEED_INVARIANT"
   SEED_NOTE="molecular seed from $ATOMIC_SEED/output, x2 $SEED_X2"
   SEED_COMMAND="# with \"Load IC? True\" in input.inp, which this case states
OMP_NUM_THREADS=$OMP_NUM_THREADS EXHALE_MOLECULAR_SEED=$ATOMIC_SEED/output \\
    ${X2_ENV:+$X2_ENV }${INV_ENV:+$INV_ENV }\\
    $EXHALE_BIN > seed.log 2>&1"
   COMMANDS="${COMMANDS}# the molecular seed, built by the binary from the certified atomic state
$SEED_COMMAND
"
   echo "[$CASE] seed $SEED_NOTE"
   env EXHALE_MOLECULAR_SEED="$ATOMIC_SEED/output" $X2_ENV $INV_ENV \
      "$EXHALE_BIN" > seed.log 2>&1
   rc=$?
   if [ $rc -ne 0 ]; then
      tail -n 20 seed.log
      fail "the molecular seed could not be built; see $DIR/seed.log"
   fi
   [ -f output/Ion_species_IC.txt ] || fail "the seed pass wrote no output/Ion_species_IC.txt"

elif [ "$LOADS" = True ] && [ "$NOSEED" != 1 ]; then
   if [ -z "$SEED" ]; then
      [ -f "$GRID" ] || fail "no $GRID: the seed cannot be mapped without the current grid"
      PICK=$(python3 "$MODELS/pick_seed.py" "$CASE" 2>seed.log)
      rc=$?
      if [ $rc -eq 3 ]; then
         cat seed.log
         fail "this case states \"Load IC? True\" and the archive carries no state of its physics"
      fi
      [ $rc -eq 0 ] || { cat seed.log; fail "pick_seed.py exited $rc"; }
      SEED=${PICK%% *}
      SEED_NOTE=$PICK
   else
      SEED_NOTE="$SEED (given)"
      [ -n "${ATTEMPT:-}" ] && [ "$ATTEMPT" -gt 1 ] && SEED_NOTE="$SEED (pick_seed's next candidate, attempt $ATTEMPT)"
   fi
   [ -n "${PRIOR_SEED_NOTE:-}" ] && SEED_NOTE="${PRIOR_SEED_NOTE}then $SEED_NOTE"
   echo "[$CASE] seed $SEED_NOTE"
   # THE COMPOSITION the seed was solved at, from its own metadata, against
   # the one this case is solved at.  A seed with no reservoir field is an
   # archived state and is mapped as it stands.
   #
   # `load_IC` compares the `# reservoir` line ELEMENT BY ELEMENT and refuses
   # a state whose ratio differs from the run's by more than 1e-6 in ANY of
   # them.  Where the case states its composition in `input.inp` that is one
   # number, He/H, and the rest of the elements are the code's own defaults,
   # which both states share.  Where the case takes a LOWER-ATMOSPHERE
   # PROFILE, every ratio is the profile's `X_<El>` column at the matching
   # level, so a seed solved under a different column has moved in C/H, N/H
   # and O/H as well: MEASURED 2026-09-16, a seed of this same case solved
   # under the previous flux-closed column was refused on C/H at 3.7e-6 while
   # He/H had been carried.  `src/utils/profile_match_level.py` states which
   # elements moved and by how much -- the same routine the elemental-flux
   # closure seeds an iteration with -- and every one of them is carried.
   SEED_STATE="$SEED/Hydro_ioniz.txt"
   [ -f "$SEED_STATE" ] || SEED_STATE="$SEED/Hydro_ioniz_IC.txt"
   PROFILE=$(sed -n 's/^Lower atmosphere profile: *//p' input.inp | tail -n 1)
   RES_OPT=""
   SEED_RESCALE=""
   if [ -n "$PROFILE" ] && [ -f "$PROFILE" ]; then
      RES_OPT=$(python3 "$EX/src/utils/profile_match_level.py" \
                   --reservoir-options "$SEED_STATE" "$PROFILE" | tr '\n' ' ')
      [ -n "$RES_OPT" ] && echo "[$CASE] the seed is carried onto the profile's composition: $RES_OPT"
   else
      SEED_HEH=$(awk '$1=="#" && $2=="reservoir" {for(i=3;i<NF;i++) if($i=="He/H"){print $(i+1); exit}}' "$SEED_STATE" 2>/dev/null)
      CASE_HEH=$(sed -n 's/^He\/H number ratio: *//p' input.inp | tail -n 1)
      if [ -n "$SEED_HEH" ] && [ -n "$CASE_HEH" ]; then
         SEED_RESCALE=$(python3 -c '
import sys
old, new = float(sys.argv[1]), float(sys.argv[2])
if old > 0.0 and new > 0.0 and abs(new/old - 1.0) > 1.0e-6:
    print("%s %s %.12f" % (sys.argv[1], sys.argv[2], new/old))
' "$SEED_HEH" "$CASE_HEH")
         [ -n "$SEED_RESCALE" ] && RES_OPT="--reservoir He/H $CASE_HEH"
      fi
      [ -n "$RES_OPT" ] && echo "[$CASE] the seed reservoir He/H $SEED_HEH is carried onto $CASE_HEH"
   fi
   # THE TARGET HEADER.  The mapper copies the "# grid" line of its target into
   # the state it writes, and load_IC compares that field against the run's as
   # text.  With a lower-atmosphere profile in use the run's R0 is the
   # profile's radius at the matching level times R_J, not the planet radius
   # the shared grid file carries, so the case gets its own copy of that file
   # with the one field replaced.  Without a profile the shared file already
   # states the run's R0 and is used as it stands.
   GRID_TARGET="$GRID"
   GRID_R0=""
   if [ -n "$PROFILE" ] && [ -f "$PROFILE" ]; then
      GRID_R0=$(python3 "$EX/src/utils/profile_match_level.py" "$PROFILE" "$GRID" target_grid_Hydro_ioniz.txt) \
         || fail "the target grid header could not be built from $PROFILE"
      GRID_TARGET="target_grid_Hydro_ioniz.txt"
      echo "[$CASE] the target grid header states R0[cm] $GRID_R0, the radius $PROFILE carries at its matching level times R_J"
   fi
   SEED_COMMAND=""
   if [ -n "$GRID_R0" ]; then
      SEED_COMMAND="# the target grid header, this case's R0 in place of the shared file's
python3 $EX/src/utils/profile_match_level.py \\
    $PROFILE $GRID target_grid_Hydro_ioniz.txt
"
   fi
   SEED_COMMAND="${SEED_COMMAND}python3 $EX/src/utils/map_state_to_grid.py \\
    $SEED \\
    $GRID_TARGET output --ic $RES_OPT"
   COMMANDS="${COMMANDS}# the seed, interpolated onto the current cell centers
$SEED_COMMAND
"
   python3 "$EX/src/utils/map_state_to_grid.py" "$SEED" "$GRID_TARGET" output --ic $RES_OPT \
      >> seed.log 2>&1 || { tail -n 20 seed.log; fail "the seed could not be mapped onto the current grid; see $DIR/seed.log"; }
fi

# ---------------------------------------------------------------- the wind --
COMMANDS="${COMMANDS}
# the wind, from the input.inp of this directory
OMP_NUM_THREADS=$OMP_NUM_THREADS EXHALE_PTC_DTAU0=$EXHALE_PTC_DTAU0 \\
    EXHALE_OUTER_PASSES=$EXHALE_OUTER_PASSES \\
    $EXHALE_BIN > run.log 2>&1
"
"$EXHALE_BIN" > run.log 2>&1
rc=$?
# Status 2 is the binary's way of saying the state was written in full but did
# not pass the stationary certification (certification.f90); the files are
# complete and the run is not a failure of this script.
if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
   tail -n 25 run.log
   fail "the wind pass exited $rc; see $DIR/run.log"
fi

# The partitioned route states its outcome once, at the end; the marching
# route states it on the solver's own line.  Each outer pass of the first
# writes a "done info=" line of its own, so that line is not the solve.
INFO=$(sed -n 's/.*stationary solve returned info *= *\(-\?[0-9]*\).*/\1/p' run.log | tail -n 1)
if [ -z "$INFO" ]; then
   DONE=$(grep -E '^ \((JFNK|PTC)\) done info=' run.log | tail -n 1)
   INFO=$(echo "$DONE" | sed -n 's/.*done info=\([0-9-]*\).*/\1/p')
   RNORM=$(echo "$DONE" | sed -n 's/.*||R||= *\([0-9.ED+-]*\).*/\1/p')
else
   RNORM=$(grep -E '^ \((JFNK|PTC)\) done info=' run.log | tail -n 1 \
           | sed -n 's/.*||R||= *\([0-9.ED+-]*\).*/\1/p')
fi
if [ -z "$INFO" ]; then
   tail -n 25 run.log
   fail "the wind pass printed no solver verdict; see $DIR/run.log"
fi
# THE CONTINUATION IN THE PSEUDO-TIME START (L4e, docs/lhs1140b_stationary_L4e_20260914.md):
# the ramp of the pseudo-transient hydro solve is tied to the line-search
# merit, which is nearly flat on a state that is already close, so from such
# a state dtau never leaves dtau0 = 1.0 and the outer loop ends refused with
# the hydrodynamic rows stuck; the same state restarted at dtau0 = 1e8
# certifies in a few Newton iterations, while a raw seed at 1e8 fails.  So:
# the first solve starts at 1.0; if it does not return info = 0, the state it
# wrote is reloaded and solved once more at DTAU0_CONTINUATION (default 1e8).
# Both runs and both logs are kept and named in REPRODUCE.md.
: "${DTAU0_CONTINUATION:=1.0e8}"
CONTINUED=0
if [ "$INFO" != 0 ] && [ -f output/Hydro_ioniz.txt ] && [ -f output/Ion_species.txt ]; then
   echo "[$CASE] the stationary solve returned info=$INFO at dtau0=$EXHALE_PTC_DTAU0; continuing from the written state at dtau0=$DTAU0_CONTINUATION"
   \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
   \cp -f output/Ion_species.txt output/Ion_species_IC.txt
   \mv -f run.log run_dtau0_first.log
   COMMANDS="${COMMANDS}
# the first solve returned info=$INFO: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=$OMP_NUM_THREADS EXHALE_PTC_DTAU0=$DTAU0_CONTINUATION \\
    $EXHALE_BIN > run.log 2>&1
"
   EXHALE_PTC_DTAU0=$DTAU0_CONTINUATION "$EXHALE_BIN" > run.log 2>&1
   rc=$?
   if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
      tail -n 25 run.log
      fail "the continuation pass exited $rc; see $DIR/run.log"
   fi
   INFO=$(sed -n 's/.*stationary solve returned info *= *\(-\?[0-9]*\).*/\1/p' run.log | tail -n 1)
   RNORM=$(grep -E '^ \((JFNK|PTC)\) done info=' run.log | tail -n 1 \
           | sed -n 's/.*||R||= *\([0-9.ED+-]*\).*/\1/p')
   CONTINUED=1
   SEED_NOTE="$SEED_NOTE; then a continuation from the first solve's written state at dtau0=$DTAU0_CONTINUATION (its log run_dtau0_first.log)"
fi
# WHAT A REFUSED SOLVE REFUSED (L14).  A state whose three hydrodynamic rows
# are all inside their own tolerances is a stationary wind; what refuses it is
# then the composition, and no seed can be nearer to the wind than the wind
# the solve already found.  MEASURED 2026-09-14 on
# atomic_scalar_gj1132x0.{10,15,20}_kzz1e9/HeH2.13: the solve from the
# certified 0.25 state ended with the three rows inside (mass 2.3e-8,
# momentum 3.4e-14, energy 8.5e-7 of 1e-6) and only the gated elemental He/H
# row refusing, and the attempt from the next seed, one further XUV step
# away, spent its passes at lam 1e-6 and dtau 1e-4 without moving the
# residual (plan item L4h).  So the seed walk below is skipped for such a
# state.
HYDRO_ROWS=$(grep -E '^   hydrodynamic (mass|momentum|energy) row' run.log \
             | tail -n 3)
HYDRO_WITHIN=0
if [ "$(echo "$HYDRO_ROWS" | grep -c 'within')" = 3 ]; then HYDRO_WITHIN=1; fi

# THE NEXT SEED. Whether a solve from a mapped seed reaches the root is a
# basin question (2026-09-14: from the certified He/H 1.60 state, 2.13 and
# 1.55 certify and 1.50 does not, while 1.50 certifies from the archived
# state of its own composition).  So when a solve ends refused and the seed
# was pick_seed's choice, the next candidate of its ranked list is tried, up
# to SEED_ATTEMPTS seeds in all; the failed attempt keeps its logs and state
# in attempt_<n>/ and REPRODUCE.md names every seed tried.
: "${SEED_ATTEMPTS:=3}"
: "${ATTEMPT:=1}"
: "${TRIED_SEEDS:=}"
if [ "$INFO" != 0 ] && [ "$HYDRO_WITHIN" = 1 ]; then
   echo "[$CASE] the solve returned info=$INFO with every hydrodynamic row inside its tolerance; what refuses is the composition, so no further seed is tried"
fi
if [ "$INFO" != 0 ] && [ "$HYDRO_WITHIN" != 1 ] && [ "$SEED_GIVEN" != 1 ] && [ "$LOADS" = True ] \
   && [ "$MOLECULAR" != 1 ] && [ "$ATTEMPT" -lt "$SEED_ATTEMPTS" ]; then
   TRIED="$TRIED_SEEDS $SEED"
   NEXT=""
   while read -r cand rest; do
      skip=0
      for t in $TRIED; do [ "$t" = "$cand" ] && skip=1; done
      [ $skip -eq 0 ] && { NEXT=$cand; break; }
   done < <(python3 "$MODELS/pick_seed.py" "$CASE" --list $((SEED_ATTEMPTS + 2)) 2>/dev/null)
   if [ -n "$NEXT" ]; then
      echo "[$CASE] the solve from seed $SEED ended info=${INFO:-?}; attempt $((ATTEMPT + 1)) from $NEXT"
      mkdir -p "attempt_$ATTEMPT"
      \mv -f output run.log run_dtau0_first.log seed.log EXHALE_setup.out EXHALE_resolved.out "attempt_$ATTEMPT/" 2>/dev/null
      echo "seed $SEED ($SEED_NOTE): stationary solve returned info=${INFO:-?}" > "attempt_$ATTEMPT/OUTCOME"
      cd "$MODELS" || exit 1
      SEED="$NEXT" ATTEMPT=$((ATTEMPT + 1)) TRIED_SEEDS="$TRIED" PRIOR_SEED_NOTE="${PRIOR_SEED_NOTE:-}attempt $ATTEMPT from $SEED ended info=${INFO:-?} (kept in attempt_$ATTEMPT/); " \
         exec "$0" "$CASE"
   fi
fi
if [ "$INFO" != 0 ]; then
   fail "the stationary solve returned info=${INFO:-?} (||R||=${RNORM:-?}); see $DIR/run.log"
fi
# The certification is a separate verdict on the state the solve returned.
# Every outer pass writes one, so it is the LAST verdict in the log that
# belongs to the state that was written.
CERT=$(awk '/NOT CERTIFIED:/ {v="uncertified"} /^ *CERTIFIED:/ {v="certified"} \
            END {print (v == "" ? "no-verdict" : v)}' run.log)

# ---------------------------------- the advection-corrected profiles --------
# Hydro_ioniz_IC.txt and Ion_species_IC.txt are the two files load_IC.f90
# reads; the other outputs of the solved state are products, not state.
#
# THE SOLVED STATE IS HANDED BACK AND MEASURED, AND NO STEP IS TAKEN.  The
# pass used to take one time step at CFL 1e-12, small enough to leave the
# state approximately where it was, because that was the only route on which
# the advection post-process ran at all.  A step is another state: it is
# taken with the reconstruction the input names rather than the one the
# solution was reached with, it re-stages the switches the solve had armed,
# and the state it writes back over the solved one is a relaxation snapshot,
# written certified=F with no stationary claim whatever the solve had
# certified.  "Restart intent: stationary evaluate" measures the state as it
# stands, writes the same conserved state back with the certification it
# just made on it, and produces the advection-corrected profiles and the
# mass-loss line from that measurement.
for f in Hydro_ioniz Ion_species; do
   [ -f "output/$f.txt" ] || fail "the wind pass wrote no output/$f.txt"
   \cp -f "output/$f.txt" "output/${f}_IC.txt"
done
\cp -f input.inp input.inp.solved
COMMANDS="${COMMANDS}
# the advection-corrected profiles: the solved state is handed back, measured
# as it stands with no step taken, and the products are written from it
for f in Hydro_ioniz Ion_species; do cp -f output/\$f.txt output/\${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' input.inp
echo 'Restart intent: stationary evaluate' >> input.inp
OMP_NUM_THREADS=$OMP_NUM_THREADS $EXHALE_BIN > pp.log 2>&1
mv -f input.inp.solved input.inp
"
sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' input.inp
echo 'Restart intent: stationary evaluate' >> input.inp
"$EXHALE_BIN" > pp.log 2>&1
pprc=$?
\mv -f input.inp.solved input.inp
# Exit 2 is the evaluation refusing the state it was handed: the solve
# certified it and the re-measurement does not, which is a statement about
# the state and not a failure of the pass (the row measures of a subsonic
# base do not round-trip bitwise,
# docs/lhs1140b_stationary_L18_20260915.md section 8).  The products are
# written either way, and the verdict of the pass is recorded below.
if [ $pprc -ne 0 ] && [ $pprc -ne 2 ]; then
   tail -n 25 pp.log
   fail "the post-processing pass exited $pprc; see $DIR/pp.log"
fi
# What the re-measurement answered, kept apart from the solve's own verdict.
PP_CLAIM=$(sed -n 's/^ *original claim: \([A-Z ]*[A-Z]\).*/\1/p' pp.log | tail -n 1)
PP_CERT=$(sed -n 's/^ *work state verdict: \([A-Z ]*[A-Z]\).*/\1/p' pp.log | tail -n 1)
MDOT=$(sed -n 's/.*steady-state Mdot *= *\([0-9.ED+-]*\).*/\1/p' pp.log | tail -n 1)

# ------------------------------------------------------ the transit spectrum -
# The WINERED HIRES-Y kernel of the Cherubim et al. (2026) observation.
COMMANDS="${COMMANDS}
# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. $LHS/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$EX python3 $EX/EXHALE_transit.py > transit.log 2>&1
"
. "$LHS/winered_hires_y.sh"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
trc=$?
if [ $trc -ne 0 ] || [ ! -f tpm_He10830_metrics.txt ]; then
   tail -n 25 transit.log
   fail "the transit synthesis exited $trc; see $DIR/transit.log"
fi

# The red-pair equivalent width over the measurement's own vacuum window,
# read the same way `LHS1140b/make_memo_figures.py` reads it.
EW=$(python3 - <<'PY'
import numpy as np
AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20
s = np.loadtxt('tpm_He10830.txt')
lam = s[:, 0]*AIR
exc = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
m = (lam >= EW_LO) & (lam <= EW_HI)
print('%.4f' % np.trapz(exc[m], lam[m]))
PY
)
DEPTH=$(awk '$1=="red_depth"{printf "%.4f", $2}' tpm_He10830_metrics.txt)

record

echo "[$CASE] DONE info=$INFO ||R||=${RNORM:-?} $CERT Mdot=${MDOT:-?} EW=${EW:-?} %A depth=${DEPTH:-?} %"
# The re-measurement of the solved state, as two answers.
echo "[$CASE] re-measured: original claim ${PP_CLAIM:-?}; the state the post-processing pass evaluated ${PP_CERT:-?}"

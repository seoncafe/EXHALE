#!/bin/bash
# Run one elemental-flux closure rung of the LHS 1140 b tree and synthesize
# the transit spectrum of the iterate it converged on.
#
#     models/run_closure.sh <group>/<case>
#
# The rung is the fixed point Phi_El = F_El(Phi_El) of the elemental flux at
# the microbar match: `src/utils/element_flux_closure.py` alternates the
# Photochem column (told a trial flux through its top) with the EXHALE wind
# (which carries a flux of its own), under-relaxed, until the two agree to
# `--tol`.  The composition is therefore an OUTPUT of the rung, and the
# `HeH<value>` of the case name is only where it started.
#
# The driver needs a solved `output/` to seed iteration 0, and does not
# support a cold start (element_flux_closure.py: "--seed is required for
# iteration 0"; every iteration is written with "Load IC? True" and the seed
# state copied in as *_IC.txt).  That seed is built here the way run_case.sh
# builds its own: models/pick_seed.py names the archived state of this rung's
# physics nearest in composition, and src/utils/map_state_to_grid.py
# interpolates it onto the current cell centers, into `seed/` beside the rung.
#
# Each iteration's wind is solved by the partitioned stationary route, which
# the rung's closure.json states under `input_keys` (MODELS.md section 6).
#
# Environment:
#   PHI0_H, PHI0_HE  starting trial fluxes [g/s].  The defaults are the
#                    converged pair of the He/H = 9.7 rung of the 2026-08-30
#                    ladder, which is the nearest solved point to this tree's
#                    starting reservoirs.
#   SEED             an output/ directory to map in place of pick_seed's choice
#   OMP_NUM_THREADS  threads for each wind solve (default 8; overrides the
#                    value in closure.json)
#   TOL, KMAX        closure tolerance and iteration limit (0.05, 8)
#   FORCE=1          re-synthesize a rung that already carries a spectrum
#
# Exit status 0 when the closure converged and the spectrum was written.

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
: "${PHI0_H:=4.7682069882E+06}"
: "${PHI0_HE:=2.0433648550E+07}"
: "${SEED:=}"
: "${OMP_NUM_THREADS:=8}"
: "${TOL:=0.05}"
: "${KMAX:=8}"
: "${FORCE:=0}"
export OMP_NUM_THREADS

STARTED=$(date '+%Y-%m-%d %H:%M:%S')
T_START=$(date +%s)
COMMANDS=""
SEED_NOTE="cold start"
SEED_COMMAND=""
record () {
   # The record of this rung, written by the run itself; a missing record is
   # reported and never fails the rung.
   [ -n "${WROTE_RECORD:-}" ] && return 0
   # Nothing to record before the solver has written anything.
   [ -f "$DIR/closure.log" ] || return 0
   WROTE_RECORD=1
   local wall=$(( $(date +%s) - T_START ))
   local iters=$(ls -d "$DIR"/k[0-9][0-9] 2>/dev/null | xargs -n1 basename 2>/dev/null | tr '\n' ' ')
   python3 "$MODELS/write_reproduce.py" "$DIR" --kind closure \
      --case-name "$CASE" --binary "$EX/EXHALE.x" \
      --threads "$OMP_NUM_THREADS" \
      --seed-line "$SEED_NOTE" --seed-command "$SEED_COMMAND" \
      --commands "$COMMANDS" --iterations "$iters" \
      --started "$STARTED" --finished "$(date '+%Y-%m-%d %H:%M:%S')" \
      --wall "$(printf '%dm%02ds' $((wall/60)) $((wall%60)))" \
      > /dev/null || echo "[$CASE] note: REPRODUCE.md was not written"
}

fail () { record; echo "[$CASE] FAILED $*"; exit 1; }

[ -d "$DIR" ]                   || fail "no such case directory $DIR"
[ -f "$DIR/closure.json" ]      || fail "$DIR has no closure.json (a prescribed case uses run_case.sh)"
[ -f "$DIR/input_template.inp" ] || fail "$DIR has no input_template.inp"

last_iterate () { ls -d "$DIR"/k[0-9][0-9] 2>/dev/null | sort | tail -n 1; }

LAST=$(last_iterate)
if [ -n "$LAST" ] && [ -f "$LAST/tpm_He10830_metrics.txt" ] && [ "$FORCE" != 1 ]; then
   echo "[$CASE] SKIP (already synthesized in $(basename "$LAST"); FORCE=1 to redo)"
   exit 0
fi

# --resume picks up the last completed iteration and continues from ITS
# solution, so the seed belongs to iteration 0 alone and is left off there.
GRID="$MODELS/current_grid_Hydro_ioniz.txt"
START=()
if [ -f "$DIR/closure_history.txt" ]; then
   START=(--resume)
   echo "[$CASE] resuming from closure_history.txt"
else
   if [ -z "$SEED" ]; then
      [ -f "$GRID" ] || fail "no $GRID: the seed cannot be mapped without the current grid"
      PICK=$(python3 "$MODELS/pick_seed.py" "$CASE" 2>"$DIR/seed.log") \
         || { cat "$DIR/seed.log"; fail "pick_seed.py found no archived state of this rung's physics"; }
      SEED=${PICK%% *}
      SEED_NOTE=$PICK
      echo "[$CASE] seed $PICK"
   else
      SEED_NOTE="$SEED (given)"
   fi
   rm -rf "$DIR/seed" && mkdir -p "$DIR/seed"
   SEED_COMMAND="python3 $EX/src/utils/map_state_to_grid.py \\
    $SEED \\
    $GRID seed"
   COMMANDS="${COMMANDS}# the seed of iteration 0, on the current cell centers
mkdir -p seed
$SEED_COMMAND
"
   python3 "$EX/src/utils/map_state_to_grid.py" "$SEED" "$GRID" "$DIR/seed" \
      >> "$DIR/seed.log" 2>&1 \
      || { tail -n 20 "$DIR/seed.log"; fail "the seed could not be mapped onto the current grid; see $DIR/seed.log"; }
   START=(--seed "$DIR/seed")
fi

cd "$DIR" || fail "cannot enter $DIR"

# The driver takes the thread count from the configuration and overwrites
# whatever the environment says, so OMP_NUM_THREADS is substituted into the
# copy that actually runs.  closure.json stays the record of the rung.
sed -e "s/\"omp_num_threads\": *[0-9]*/\"omp_num_threads\": $OMP_NUM_THREADS/" \
   closure.json > closure_run.json

COMMANDS="${COMMANDS}
# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/\"omp_num_threads\": *[0-9]*/\"omp_num_threads\": $OMP_NUM_THREADS/' \\
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 $EX/src/utils/element_flux_closure.py . \\
    --config closure_run.json \\
    --phi0-H $PHI0_H --phi0-He $PHI0_HE \\
    --tol $TOL --kmax $KMAX ${START[*]} > closure_stdout.log 2>&1
"
/opt/miniconda3/bin/python3 "$EX/src/utils/element_flux_closure.py" . \
   --config closure_run.json \
   --phi0-H "$PHI0_H" --phi0-He "$PHI0_HE" \
   --tol "$TOL" --kmax "$KMAX" \
   "${START[@]}" > closure_stdout.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then
   tail -n 25 closure.log
   fail "the closure driver exited $rc; see $DIR/closure.log"
fi

LAST=$(last_iterate)
[ -n "$LAST" ] || fail "the driver wrote no iteration directory"

# ------------------------------------------------------ the transit spectrum -
COMMANDS="${COMMANDS}
# the transit spectrum of the iterate the rung converged on
. $LHS/winered_hires_y.sh
cd \$(ls -d k[0-9][0-9] | sort | tail -n 1)
MPLBACKEND=Agg PYTHONPATH=$EX python3 $EX/EXHALE_transit.py > transit.log 2>&1
"
. "$LHS/winered_hires_y.sh"
cd "$LAST" || fail "cannot enter $LAST"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
trc=$?
if [ $trc -ne 0 ] || [ ! -f tpm_He10830_metrics.txt ]; then
   tail -n 25 transit.log
   fail "the transit synthesis exited $trc; see $LAST/transit.log"
fi

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

# The converged reservoir and the mass-loss rate are the last row of the
# history: k, trial fluxes, measured fluxes, residuals, omega, window, spreads,
# HeH_match, log10 Mdot, solution_id, info.
read -r K HEH MDOT INFO <<EOF
$(awk '!/^#/ && NF>=15 {k=$1; heh=$12; mdot=$13; info=$15} END {print k, heh, mdot, info}' "$DIR/closure_history.txt")
EOF

record

echo "[$CASE] DONE k=$K info=$INFO He/H=$HEH Mdot=$MDOT EW=${EW:-?} %A depth=${DEPTH:-?} %"

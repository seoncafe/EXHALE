#!/bin/bash
# Solve the XUV-scaled cases at the RAISED pseudo-time start.
#
#     models/run_lowxuv_dtau0.sh <list-file> [JOBS] [THREADS]
#
# WHY THIS EXISTS (item L14, and `run_case.sh` says the same at its
# EXHALE_PTC_DTAU0 default). Below about a third of the fiducial GJ 1132
# spectrum the composition enters far from its fixed point, and from there
# the dtau0 = 1 solve falls into the ramp collapse of item L4h: MEASURED on
# the corrected base boundary 2026-09-15/16, forty outer passes with the
# gated species row already at 1e-12 and the hydrodynamic ENERGY row of the
# OUTERMOST physical cell above its tolerance (item L4f), thirty to ninety
# minutes a pass, no root. The SAME mapped seed at dtau0 = 1e8 reaches its
# root in about forty Newton iterations.
#
# So these cases are taken there directly. `SEED_ATTEMPTS=1`, because what
# refuses is the pseudo-time ramp and not the seed, and `FORCE=1`, because
# the case may already carry the products of an attempt that was stopped.
# Each case's record states why the run was made this way; the logs of the
# stopped attempt belong in `models/.stopped/<case>_dtau0-1_<timestamp>/`.
#
# The list file holds one `<group>/HeH<value>` per line.
set -u
M=$(cd "$(dirname "$0")" && pwd)
LIST=${1:?usage: run_lowxuv_dtau0.sh <list-file> [JOBS] [THREADS]}
JOBS=${2:-6}; THREADS=${3:-8}
cd "$M" || exit 1
LOG=$M/lowxuv_dtau0_$(date +%Y%m%d).log
echo "low-XUV solve start $(date '+%F %T') on $(hostname -s): $(wc -l < "$LIST") cases at dtau0=1.0e8, $JOBS jobs x $THREADS threads" | tee -a "$LOG"
export OMP_NUM_THREADS=$THREADS
export EXHALE_PTC_DTAU0=1.0e8
export FORCE=1
export SEED_ATTEMPTS=1
xargs -a "$LIST" -P "$JOBS" -I{} bash -c '
   c={}
   s=$(ls -d "'"$M"'/.stopped/"$(echo "$c" | tr / _)"_dtau0-1_"* 2>/dev/null | sort | tail -n 1)
   RUN_NOTE="The first solve of this re-run was taken at the default pseudo-time start (\`EXHALE_PTC_DTAU0=1.0\`) and refused: forty outer passes with the gated species row already at 1e-12 and the hydrodynamic ENERGY row of the outermost physical cell above its tolerance (item L4f), which is the ramp collapse of item L4h. It was stopped${s:+, and its logs are in \`models/.stopped/$(basename "$s")/\`}. This run is what \`run_case.sh\` prescribes for that situation (item L14): the same mapped seed at \`EXHALE_PTC_DTAU0=1.0e8\`."
   export RUN_NOTE
   r=$("'"$M"'/run_case.sh" "$c" 2>&1 | tail -n 1)
   echo "$(date +%T) $c: $r"' >> "$LOG" 2>&1
echo "low-XUV solve end $(date '+%F %T') on $(hostname -s)" | tee -a "$LOG"

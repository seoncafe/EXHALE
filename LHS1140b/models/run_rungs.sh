#!/bin/bash
# Run every flux-closure rung of the tree that has no line yet, JOBS at a
# time, each rung at THREADS OpenMP threads (a rung alternates Photochem and
# EXHALE, so its wall time is several EXHALE solves).
#
#     models/run_rungs.sh [JOBS] [THREADS] [pattern]
#
# pattern (default 'atomic_photochem*') selects the groups; run_closure.sh
# skips a rung whose last iterate already carries tpm_He10830_metrics.txt.
# One line per rung goes to rungs_<date>.log beside this script.
set -u
JOBS=${1:-4}; THREADS=${2:-8}; PAT=${3:-atomic_photochem*}
HERE=$(cd "$(dirname "$0")" && pwd)
LOG=$HERE/rungs_$(date +%Y%m%d).log
cd "$HERE" || exit 1
# Named for the host and the process: `xargs` reads the list as the run
# proceeds and the tree is reachable from every machine of this NFS mount.
LIST=$HERE/.rung_list_$(hostname -s)_$$
# CAMPAIGN_LIST names a file of `<group>/HeH<value>` lines to run instead of
# selecting them with a pattern; the pattern `atomic_photochem*` also picks up
# the XUV-scaled PRESCRIBED cases above the fiducial column, which are not
# rungs and which `run_closure.sh` refuses by name.
: "${CAMPAIGN_LIST:=}"
if [ -n "$CAMPAIGN_LIST" ]; then
   sed '/^[[:space:]]*$/d' "$CAMPAIGN_LIST" > "$LIST"
else
   ls -d $PAT/HeH* | sort > "$LIST"
fi
echo "rungs start $(date '+%F %T') on $(hostname -s): $(wc -l < "$LIST") rungs, $JOBS jobs x $THREADS threads" | tee -a "$LOG"
export OMP_NUM_THREADS=$THREADS
xargs -a "$LIST" -P "$JOBS" -I{} bash -c 'r=$("'"$HERE"'/run_closure.sh" {} 2>&1 | tail -n 1); echo "$(date +%T) {}: $r"' >> "$LOG" 2>&1
echo "rungs end $(date '+%F %T') on $(hostname -s)" | tee -a "$LOG"

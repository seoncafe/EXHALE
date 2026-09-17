#!/bin/bash
# Run every prescribed-composition atomic case of the tree that has no line
# yet, JOBS at a time, each at THREADS OpenMP threads (default 8 x 8 = 64 of
# the 72 cores, leaving room for the closure rungs and diagnostics).
#
#     models/run_campaign.sh [JOBS] [THREADS] [pattern]
#
# pattern (default 'atomic_scalar*') selects the groups; run_case.sh skips a
# case that already carries tpm_He10830_metrics.txt. One line per case goes
# to campaign_<date>.log beside this script; each case keeps its own logs and
# REPRODUCE.md.
#
# CAMPAIGN_LIST names a file of `<group>/HeH<value>` lines to run instead of
# selecting them with a pattern, for a run over cases that no pattern picks
# out -- the final pass of 2026-09-16 re-solved the cases that had certified
# on the previous binary, which is a list and not a glob.
#
# CAMPAIGN_EXCLUDE, an extended regular expression, drops cases the pattern
# selected: a group the pattern cannot separate is named here instead of
# splitting the run. The 2026-09-15 re-run leaves out the three cases at 0.01
# of the fiducial XUV with CAMPAIGN_EXCLUDE='x0\.01' (they never solved, item
# L17, and carry no state to continue from; models/README.md says so).
set -u
JOBS=${1:-8}; THREADS=${2:-8}; PAT=${3:-atomic_scalar*}
: "${CAMPAIGN_EXCLUDE:=}"
HERE=$(cd "$(dirname "$0")" && pwd)
LOG=$HERE/campaign_$(date +%Y%m%d).log
cd "$HERE" || exit 1
# The list is named for the host and the process, because `xargs` reads it
# as the run proceeds and two campaigns can be in flight at once -- the same
# tree is reachable from every machine of this NFS mount.
LIST=$HERE/.campaign_list_$(hostname -s)_$$
: "${CAMPAIGN_LIST:=}"
if [ -n "$CAMPAIGN_LIST" ]; then
   sed '/^[[:space:]]*$/d' "$CAMPAIGN_LIST" > "$LIST.selected"
else
   ls -d $PAT/HeH* | sort > "$LIST.selected"
fi
if [ -n "$CAMPAIGN_EXCLUDE" ]; then
   grep -v -E "$CAMPAIGN_EXCLUDE" "$LIST.selected" > "$LIST"
   echo "CAMPAIGN_EXCLUDE='$CAMPAIGN_EXCLUDE' drops $(( $(wc -l < "$LIST.selected") - $(wc -l < "$LIST") )) of $(wc -l < "$LIST.selected") selected cases" | tee -a "$LOG"
else
   \cp -f "$LIST.selected" "$LIST"
fi
echo "campaign start $(date '+%F %T') on $(hostname -s): $(wc -l < "$LIST") cases, $JOBS jobs x $THREADS threads" | tee -a "$LOG"
export OMP_NUM_THREADS=$THREADS
xargs -a "$LIST" -P "$JOBS" -I{} bash -c 'r=$("'"$HERE"'/run_case.sh" {} 2>&1 | tail -n 1); echo "$(date +%T) {}: $r"' >> "$LOG" 2>&1
echo "campaign end $(date '+%F %T') on $(hostname -s)" | tee -a "$LOG"

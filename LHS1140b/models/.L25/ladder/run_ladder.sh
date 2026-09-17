#!/bin/bash
# One continuation in XUV: each rung solved from the certified state of the
# rung above it (L25 step 4).  A rung that does not certify stops the ladder.
#
#     run_ladder.sh <ladder> <seed-state-dir> <tag> [tag ...]
#
# Every rung is the catalog input of the ladder's TARGET case with only
# `Spectrum file:` changed to the rung's scaled spectrum, so the flux (ROE),
# the base boundary, the tolerance and the lower-atmosphere profile are the
# target's throughout; the XUV is the only thing that moves.  The recipe is
# the catalog's: `EXHALE_PTC_DTAU0=1.0e8` (item L14, the low-XUV cases),
# forty outer passes, one seed, 8 threads.
set -u
LAD=$(cd "$(dirname "$0")" && pwd)
M=$(dirname "$(dirname "$LAD")")
NAME=$1; SEED=$2; shift 2
LOG=$LAD/$NAME.ladder.log
HOST=$(hostname -s)
echo "=== ladder $NAME start $(date '+%F %T') on $HOST, seed $SEED" >> "$LOG"
export OMP_NUM_THREADS=8
export EXHALE_PTC_DTAU0=1.0e8
export EXHALE_OUTER_PASSES=40
export SEED_ATTEMPTS=1
export FORCE=1
for tag in "$@"; do
   D=$LAD/$NAME/x$tag
   export SEED
   export RUN_NOTE="L25 step 4, continuation in XUV: rung x$tag of ladder $NAME, the catalog input of the target case with only the spectrum changed, solved from the certified state of the rung above it ($SEED). Host $HOST, 8 threads."
   t0=$(date +%s)
   "$M/run_case.sh" ".L25/ladder/$NAME/x$tag" > "$D.case.log" 2>&1
   rc=$?
   r=$(tail -n 3 "$D.case.log" | tr '\n' ' ')
   w=$(( $(date +%s) - t0 ))
   cert=$(awk '/NOT CERTIFIED:/ {v="uncertified"} /^ *CERTIFIED:/ {v="certified"} END {print (v == "" ? "no-verdict" : v)}' "$D/run.log" 2>/dev/null)
   echo "$(date '+%F %T') $HOST $NAME x$tag rc=$rc cert=$cert wall=${w}s :: $r" >> "$LOG"
   info=$(sed -n 's/.*stationary solve returned info *= *\(-\?[0-9]*\).*/\1/p' "$D/run.log" 2>/dev/null | tail -n 1)
   echo "$(date '+%F %T') $HOST $NAME x$tag info=${info:-?}" >> "$LOG"
   if [ "$cert" != certified ] || [ "${info:-1}" != 0 ]; then
      echo "$(date '+%F %T') $HOST ladder $NAME STOPS at x$tag ($cert)" >> "$LOG"
      exit 1
   fi
   SEED=$D/output
done
echo "=== ladder $NAME complete $(date '+%F %T') on $HOST" >> "$LOG"

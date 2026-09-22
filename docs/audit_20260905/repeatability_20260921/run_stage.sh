#!/bin/bash
# One diagnostic invocation in its own directory, with the sidecar the gate of
# PLAN_20260920_rev9 section 4.2 requires.
# usage: run_stage.sh <atomic|molecular> <name> <OMP> <OPENBLAS> <stage> [QUAD] [coupled]
#   stage: determinism | residual | dumpic
set -u
SC="$(cd "$(dirname "$0")" && pwd)"
CASE="$1"; NAME="$2"; OMP="$3"; OBLAS="$4"; STAGE="$5"; QUAD="${6:-0}"; COUP="${7:-no}"
case "$CASE" in
  atomic)    SRC="$SC/tree/models/atomic/HeH9.7";;
  molecular) SRC="$SC/tree/models/molecular/HeH0.083";;
  *) echo "unknown case $CASE"; exit 2;;
esac
case "$STAGE" in
  determinism) STAGEVAR="EXHALE_RESID_DETERMINISM=1";;
  residual)    STAGEVAR="EXHALE_RESIDUAL=1";;
  dumpic)      STAGEVAR="EXHALE_DUMP_IC=1";;
  newton)      STAGEVAR="EXHALE_NEWTON_TEST=1";;
  *) echo "unknown stage $STAGE"; exit 2;;
esac
RUN="$SC/runs/$NAME"
rm -rf "$RUN"
mkdir -p "$RUN/models/$CASE/work/output"
ln -sfn /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed "$RUN/sed"
W="$RUN/models/$CASE/work"
for f in "$SRC"/*.inp; do \cp -f "$f" "$W/"; done
for f in "$SRC"/output/*_IC.txt; do \cp -f "$f" "$W/output/"; done
if [ "$COUP" = "coupled" ]; then
   echo 'Coupled carrier solve: True' >> "$W/input.inp"
fi
\cp -f "$SC/EXHALE.x" "$W/EXHALE.x"
START=$(date -u +%Y-%m-%dT%H:%M:%SZ); T0=$(date +%s.%N)
( cd "$W" && env -i PATH=/usr/bin:/bin HOME="$HOME" \
      OMP_NUM_THREADS="$OMP" OPENBLAS_NUM_THREADS="$OBLAS" \
      OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 \
      $STAGEVAR EXHALE_RESID_QUAD="$QUAD" \
      ./EXHALE.x ) > "$RUN/stdout.log" 2> "$RUN/stderr.log"
RC=$?
T1=$(date +%s.%N); END=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{
  echo "invocation        $NAME"
  echo "mode              R (current operator on a stored generation; no physical step taken, nothing published)"
  echo "case              $CASE"
  echo "stage             $STAGEVAR"
  echo "route change      $COUP"
  echo "run directory     $W"
  echo "executable        $W/EXHALE.x  md5 $(md5sum "$W/EXHALE.x" | cut -d' ' -f1)"
  echo "environment in force (env -i, so this is the complete set):"
  echo "  PATH=/usr/bin:/bin  HOME=$HOME"
  echo "  OMP_NUM_THREADS=$OMP  OPENBLAS_NUM_THREADS=$OBLAS  OMP_DYNAMIC=FALSE  OMP_MAX_ACTIVE_LEVELS=1"
  echo "  $STAGEVAR  EXHALE_RESID_QUAD=$QUAD"
  echo "input state (md5 MEASURED):"
  for f in "$W"/*.inp "$W"/output/*_IC.txt; do
     echo "  $(md5sum "$f" | sed "s#$W/##")"
  done
  echo "started           $START"
  echo "ended             $END"
  echo "wall seconds      $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.2f", b-a}')"
  echo "exit status       $RC"
  echo "files this invocation left (size, path, md5 MEASURED):"
  for f in "$RUN/stdout.log" "$RUN/stderr.log" "$W"/output/*; do
     [ -e "$f" ] || continue
     echo "  $(stat -c '%s  %n' "$f")  $(md5sum "$f" | cut -d' ' -f1)  rows=$(wc -l < "$f")"
  done
  echo "verdict line printed:"
  grep 'VERDICT' "$RUN/stdout.log" | sed 's/^/  /' || echo "  (none printed by this stage)"
} > "$RUN/SIDECAR.txt"
if [ "$STAGE" = determinism ]; then
   sed -n '/resid_determinism/,$p' "$RUN/stdout.log" > "$RUN/determinism_block.txt"
fi
echo "== $NAME rc=$RC $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.1fs", b-a}')  $(grep -m1 'VERDICT' "$RUN/stdout.log" || echo '(no verdict line)')"

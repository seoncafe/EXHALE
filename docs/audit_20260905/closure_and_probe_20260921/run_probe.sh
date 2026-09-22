#!/bin/bash
# One closure-and-arc invocation in its own directory, with the sidecar the
# gate of PLAN_20260920_rev9 section 4.2 requires.
# usage: run_probe.sh <atomic|molecular> <name> <OMP> <OPENBLAS> <stage> [EXE]
set -u
SC="$(cd "$(dirname "$0")" && pwd)"
TREE=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
CASE="$1"; NAME="$2"; OMP="$3"; OBLAS="$4"; STAGE="$5"; EXE="${6:-$SC/EXHALE_test.x}"
case "$CASE" in
  atomic)    SRC="$TREE/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7"
             GEN="g0002_20260919T004843Z_f7485b14";;
  molecular) SRC="$TREE/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083"
             GEN="g0004_20260920T053734Z_79b42a03";;
  *) echo "unknown case $CASE"; exit 2;;
esac
case "$STAGE" in
  closure)     STAGEVAR="EXHALE_CLOSURE_PROBE=1";;
  determinism) STAGEVAR="EXHALE_RESID_DETERMINISM=1";;
  residual)    STAGEVAR="EXHALE_RESIDUAL=1";;
  *) echo "unknown stage $STAGE"; exit 2;;
esac
RUN="$SC/runs/$NAME"
rm -rf "$RUN"; mkdir -p "$RUN/models/$CASE/work/output"
ln -sfn "$TREE/LHS1140b/sed" "$RUN/sed"
W="$RUN/models/$CASE/work"
for f in "$SRC"/*.inp; do \cp -f "$f" "$W/"; done
\cp -f "$SRC/states/$GEN/Hydro_ioniz.txt" "$W/output/Hydro_ioniz_IC.txt"
\cp -f "$SRC/states/$GEN/Ion_species.txt" "$W/output/Ion_species_IC.txt"
\cp -f "$EXE" "$W/EXHALE.x"
START=$(date -u +%Y-%m-%dT%H:%M:%SZ); T0=$(date +%s.%N)
( cd "$W" && env -i PATH=/usr/bin:/bin HOME="$HOME" \
      OMP_NUM_THREADS="$OMP" OPENBLAS_NUM_THREADS="$OBLAS" \
      OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 \
      $STAGEVAR EXHALE_RESID_QUAD=0 \
      ./EXHALE.x ) > "$RUN/stdout.log" 2> "$RUN/stderr.log"
RC=$?
T1=$(date +%s.%N); END=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{
  echo "invocation        $NAME"
  echo "mode              R (current operator on a stored generation; no physical step taken, nothing published)"
  echo "case              $CASE"
  echo "generation        $GEN"
  echo "stage             $STAGEVAR"
  echo "run directory     $W"
  echo "executable        $W/EXHALE.x  md5 $(md5sum "$W/EXHALE.x" | cut -d' ' -f1)"
  echo "environment in force (env -i, so this is the complete set):"
  echo "  PATH=/usr/bin:/bin  HOME=$HOME"
  echo "  OMP_NUM_THREADS=$OMP  OPENBLAS_NUM_THREADS=$OBLAS  OMP_DYNAMIC=FALSE  OMP_MAX_ACTIVE_LEVELS=1"
  echo "  $STAGEVAR  EXHALE_RESID_QUAD=0"
  echo "input state (md5 MEASURED):"
  for f in "$W"/*.inp "$W"/output/*_IC.txt; do
     echo "  $(md5sum "$f" | sed "s#$W/##")"
  done
  echo "started           $START"
  echo "ended             $END"
  echo "wall seconds      $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.2f", b-a}')"
  echo "exit status       $RC"
  echo "evaluation stage reached (the block the stage prints):"
  if grep -q 'closure_probe) ---- end of the grid ----' "$RUN/stdout.log"; then
     echo "  the grid completed"
  elif grep -q 'closure_probe' "$RUN/stdout.log"; then
     echo "  the grid STARTED and did not reach its end"
  else
     echo "  no diagnostic block printed"
  fi
  echo "files this invocation left (size, path, md5 MEASURED):"
  for f in "$RUN/stdout.log" "$RUN/stderr.log" "$W"/output/*; do
     [ -e "$f" ] || continue
     echo "  $(stat -c '%s  %n' "$f")  $(md5sum "$f" | cut -d' ' -f1)  rows=$(wc -l < "$f")"
  done
} > "$RUN/SIDECAR.txt"
grep 'closure_probe' "$RUN/stdout.log" > "$RUN/closure_block.txt" 2>/dev/null
echo "== $NAME rc=$RC $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.1fs", b-a}') lines=$(wc -l < "$RUN/closure_block.txt")"

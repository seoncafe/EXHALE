#!/bin/bash
# One coupled closure-and-arc invocation in its own directory, with the
# sidecar the gate of PLAN_20260920_rev9 section 4.2 requires.
# usage: run_coupled_probe.sh <case> <name> <OMP> <OPENBLAS> <coupled T|F> [named direction]
set -u
SC="$(cd "$(dirname "$0")" && pwd)"
TREE=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
CASE="$1"; NAME="$2"; OMP="$3"; OBLAS="$4"; COUP="$5"; NAMEDDIR="${6:-}"; EXTRA="${7:-EXHALE_UNUSED_MARKER2=0}"
EXE="$SC/EXHALE_test.x"
case "$CASE" in
  kzz1e9)   SRC="$TREE/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083"
            GEN="g0004_20260920T053734Z_79b42a03"
            HY="$SRC/states/$GEN/Hydro_ioniz.txt"; IO="$SRC/states/$GEN/Ion_species.txt";;
  wellmixed) SRC="$TREE/LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55"
            GEN="the suite fixture, output/*_IC.txt = state g0002_20260915T221824Z_d7347336"
            HY="$SRC/output/Hydro_ioniz_IC.txt"; IO="$SRC/output/Ion_species_IC.txt";;
  *) echo "unknown case $CASE"; exit 2;;
esac
RUN="$SC/runs/$NAME"
rm -rf "$RUN"; mkdir -p "$RUN/work/output"
W="$RUN/work"
for f in "$SRC"/*.inp; do \cp -f "$f" "$W/"; done
\cp -f "$HY" "$W/output/Hydro_ioniz_IC.txt"
\cp -f "$IO" "$W/output/Ion_species_IC.txt"
sed -i "s|^Spectrum file: .*|Spectrum file: $TREE/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" "$W/input.inp"
sed -i "/^Coupled carrier solve:/d" "$W/input.inp"
if [ "$COUP" = T ]; then echo "Coupled carrier solve: True" >> "$W/input.inp"; fi
\cp -f "$EXE" "$W/EXHALE.x"
if [ -n "$NAMEDDIR" ]; then DIRVAR="EXHALE_CLOSURE_PROBE_DIRECTION=$NAMEDDIR"; else DIRVAR="EXHALE_UNUSED_MARKER=0"; fi
START=$(date -u +%Y-%m-%dT%H:%M:%SZ); T0=$(date +%s.%N)
( cd "$W" && env -i PATH=/usr/bin:/bin HOME="$HOME" \
      OMP_NUM_THREADS="$OMP" OPENBLAS_NUM_THREADS="$OBLAS" \
      OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 \
      EXHALE_CLOSURE_PROBE=1 EXHALE_RESID_QUAD=0 $DIRVAR $EXTRA \
      ./EXHALE.x ) > "$RUN/stdout.log" 2> "$RUN/stderr.log"
RC=$?
T1=$(date +%s.%N); END=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{
  echo "invocation        $NAME"
  echo "mode              R (current operator on an immutable stored state; no physical step taken, nothing published)"
  echo "case              $CASE   $SRC"
  echo "state             $GEN"
  echo "coupled           Coupled carrier solve: $COUP"
  echo "named direction   ${NAMEDDIR:-none}"
  echo "run directory     $W"
  echo "executable        $W/EXHALE.x  md5 $(md5sum "$W/EXHALE.x" | cut -d' ' -f1)"
  echo "environment in force (env -i, so this is the complete set):"
  echo "  PATH=/usr/bin:/bin  HOME=$HOME"
  echo "  OMP_NUM_THREADS=$OMP  OPENBLAS_NUM_THREADS=$OBLAS  OMP_DYNAMIC=FALSE  OMP_MAX_ACTIVE_LEVELS=1"
  echo "  EXHALE_CLOSURE_PROBE=1  EXHALE_RESID_QUAD=0  $DIRVAR  $EXTRA"
  echo "input state (md5 MEASURED):"
  for f in "$W"/*.inp "$W"/output/*_IC.txt; do echo "  $(md5sum "$f" | sed "s#$W/##")"; done
  echo "started           $START"
  echo "ended             $END"
  echo "wall seconds      $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.2f", b-a}')"
  echo "exit status       $RC"
  echo "evaluation stage reached:"
  if grep -q 'closure_probe) ---- end of the grid ----' "$RUN/stdout.log"; then echo "  the grid completed"
  elif grep -q 'closure_probe' "$RUN/stdout.log"; then echo "  the grid STARTED and did not reach its end"
  else echo "  no diagnostic block printed"; fi
  echo "files this invocation left (size, path, md5 MEASURED):"
  for f in "$RUN/stdout.log" "$RUN/stderr.log" "$W"/output/*; do
     [ -e "$f" ] || continue
     echo "  $(stat -c '%s  %n' "$f")  $(md5sum "$f" | cut -d' ' -f1)"
  done
} > "$RUN/SIDECAR.txt"
grep 'closure_probe' "$RUN/stdout.log" > "$RUN/closure_block.txt" 2>/dev/null
echo "== $NAME rc=$RC $(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.1fs", b-a}') lines=$(wc -l < "$RUN/closure_block.txt")"

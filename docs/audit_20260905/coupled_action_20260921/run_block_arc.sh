#!/bin/bash
# The coupled-block Jacobian action of src/tests/coupled_block_jacobian at a
# named probe arc, on the suite's own fixture and with the suite's own
# statistic, so that the failing row can be read as a function of the arc.
# usage: run_block_arc.sh <carrier|energy> <arc>
set -u
SC="$(cd "$(dirname "$0")" && pwd)"
TREE=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
DIR="$1"; ARC="$2"
CASE="$TREE/LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55"
SED="$TREE/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt"
OUT="$SC/blockarc/${DIR}_arc${ARC}"
rm -rf "$OUT"; mkdir -p "$OUT/output"
\cp -f "$CASE/input.inp" "$CASE/base.inp" "$OUT/"
\cp -f "$CASE"/output/*.txt "$OUT/output/"
sed -i "s|^Spectrum file: .*|Spectrum file: $SED|" "$OUT/input.inp"
sed -i "s|^Restart intent: .*|Restart intent: stationary|" "$OUT/input.inp"
\cp -f "$SC/EXHALE_test.x" "$OUT/EXHALE.x"
T0=$(date +%s.%N)
( cd "$OUT" && env -i PATH=/usr/bin:/bin HOME="$HOME" OMP_NUM_THREADS=1 \
     EXHALE_COUPLED_JAC_ACTION="$OUT/jac.txt" \
     EXHALE_COUPLED_JAC_CELLS=300,312 EXHALE_COUPLED_JAC_DIR="$DIR" \
     EXHALE_JV_PROBE_ARC="$ARC" ./EXHALE.x ) > "$OUT/run.log" 2>&1
RC=$?; T1=$(date +%s.%N)
PL=$(awk '/^# probe length/ {print $4}' "$OUT/jac.txt" 2>/dev/null)
MV=$(awk '/^# mean relative error over the rows the direction moves/ {print $(NF-4), $(NF-2), $NF}' "$OUT/jac.txt" 2>/dev/null)
SP=$(awk '/over every row of the support/ {print $(NF-4), $(NF-2), $NF}' "$OUT/jac.txt" 2>/dev/null)
RM=$(awk '/^# largest response in the support/ {print $7}' "$OUT/jac.txt" 2>/dev/null)
echo "dir=$DIR arc=$ARC rc=$RC wall=$(awk -v a=$T0 -v b=$T1 'BEGIN{printf "%.0f", b-a}')s probe_length=$PL moved(mean worst n)=$MV support(mean worst n)=$SP largest_response=$RM"

#!/bin/bash
# One arm of the elemental-flux-closure reservoir bracket, run on the
# environment that carries the epsfcn repair (env/photochem), so that the
# arm can be differenced against the stored ../crossings_pc090 arm, which
# was run on photochem_090_fix -- the same photochem 0.9.0 source without
# that repair.
#
# Identical to ../crossings_pc090/run_closure.sh in every respect but two:
# the "python" entry of closure.json points at env/photochem, and
# omp_num_threads is 3 rather than 4 because three other jobs share this
# machine.  The ionization sweep is byte-identical across thread counts.
#
# usage: ./run_arm_installed.sh <tag> <reservoir He/H>
#   TOL=<value>  closure tolerance (default 0.05, the stored ladder's value)
#   KMAX=<n>     iteration ceiling (default 8)
set -u
tag=$1; resv=$2
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
FC=$EX/LHS1140b/exhale/flux_closure
seed=${3:-$FC/heh9p7/k03/output}
phiH=${4:-4.7682069882E+06}
phiHe=${5:-2.0433648550E+07}
PY=$EX/env/photochem/bin/python
DRV=$EX/src/utils/element_flux_closure.py
. "$EX/LHS1140b/winered_hires_y.sh"          # He I 10830 at R = 68,000
cd "$(dirname "$0")" || exit 1
HERE=$PWD
d=$HERE/$tag
mkdir -p "$d"

sed -e "s/He=9.7/He=$resv/" \
    -e "s|$FC/heh9p7/input_template.inp|$d/input_template.inp|" \
    -e "s/He\/H = 9.7 as the STARTING/He\/H = $resv as the STARTING/" \
    -e "s|/home/kiseon/.conda/envs/photochem_cmp/bin/python|$PY|" \
    -e "s/\"omp_num_threads\": 4/\"omp_num_threads\": 3/" \
    "$FC/heh9p7/closure.json" > "$d/closure.json"
sed -e "s|^He/H number ratio: .*|He/H number ratio: $resv|" \
    -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    "$FC/heh9p7/input_template.inp" > "$d/input_template.inp"

echo "[$tag] CLOSURE START reservoir=$resv seed=$seed  $(date +%H:%M:%S)"
resume=""; [ -f "$d/closure_history.txt" ] && resume="--resume"
python3 "$DRV" "$d" --phi0-H "$phiH" --phi0-He "$phiHe" \
    --config "$d/closure.json" --seed "$seed" $resume \
    --tol ${TOL:-0.05} --kmax ${KMAX:-8} > "$HERE/$tag.closure.log" 2>&1
rc=$?
last=$(\ls -d --color=never $d/k[0-9][0-9] 2>/dev/null | sort | tail -n 1)
echo "[$tag] CLOSURE rc=$rc last=$last"
[ $rc -eq 0 ] && [ -n "$last" ] || { echo "[$tag] CLOSURE FAILED"; exit 1; }
( cd "$last" && MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" \
    > transit.log 2>&1 )
echo "[$tag] TRANSIT rc=$? in $last  $(date +%H:%M:%S)"

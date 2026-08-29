#!/bin/bash
# usage: ./run_tight.sh <tag> <guess> <longdy> <longdydt> <eqtime>
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
PY=$EX/env/photochem/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=1
tag=$1; guess=$2; longdy=$3; longdydt=$4; eqt=$5
d=$HERE/runs/$tag
[ -f "$d/lower_atmosphere_profile.dat" ] && { echo "[$tag] have it"; exit 0; }
mkdir -p "$d"
( cd "$d" && $PY $HERE/run_tightened.py "$longdy" "$longdydt" "$eqt" -- . \
    --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 \
    --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 \
    --climate-t-deep-guess "$guess" \
    --stellar-flux $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
    --flux-at-planet --wavelength-unit A --toa 1.0e-2 \
    --atoms H,He,N,O,C --abundances He=9.0 --kzz-const 1.0e9 \
    --trial-flux-H 4768206.9882 --trial-flux-He 20433648.55 \
    --blocks 400 --iteration 0 ) > "$d/adapter.log" 2>&1
echo "[$tag] guess=$guess longdy=$longdy rc=$? $(date +%H:%M:%S)"

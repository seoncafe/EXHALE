#!/bin/bash
# One handoff column per initial guess for the climate solve.
#
# The configuration is that of ../nh_refusal_diagnosis/scan_reservoir_pc090.sh
# at He/H = 9.0 -- the closure crossing's own arm -- in every respect but
# --climate-t-deep-guess.  The guess does not enter the physics: it is the
# starting point MINPACK is given for the deep boundary temperature, so every
# run below is the SAME radiative-convective solution, reached from a
# different place.  What separates the columns is therefore the solve's own
# convergence floor and nothing else.
#
# usage: ./run_guess.sh <tag> <guess> [<tag> <guess> ...]
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
PY=$EX/env/photochem/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=1
while [ $# -ge 2 ]; do
    tag=$1; guess=$2; shift 2
    d=$HERE/runs/$tag
    [ -f "$d/lower_atmosphere_profile.dat" ] && { echo "[$tag] have it"; continue; }
    mkdir -p "$d"
    ( cd "$d" && $PY $EX/src/utils/photochem_to_lower_profile.py . \
        --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 \
        --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 \
        --climate-t-deep-guess "$guess" \
        --stellar-flux $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
        --flux-at-planet --wavelength-unit A --toa 1.0e-2 \
        --atoms H,He,N,O,C --abundances He=9.0 --kzz-const 1.0e9 \
        --trial-flux-H 4768206.9882 --trial-flux-He 20433648.55 \
        --iteration 0 ) > "$d/adapter.log" 2>&1
    echo "[$tag] guess=$guess rc=$? $(date +%H:%M:%S)"
done

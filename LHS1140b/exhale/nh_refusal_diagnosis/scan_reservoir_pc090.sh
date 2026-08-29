#!/bin/bash
# The scan of scan_reservoir.sh, repeated on the corrected Photochem build.
#
# Same configuration in every respect but the interpreter: photochem 0.9.0
# with the element-relative Equilibrate convergence test and the bracketed,
# scaled Clima solves, instead of the 0.8.4 conda package.  The conservation
# check is again disarmed (--abundance-tol 1.0) so the profile is written and
# the departure can be read from its deepest row; the DEFAULT tolerance in
# the code is not touched here.
#
# usage: ./scan_reservoir_pc090.sh <heh> [<heh> ...]
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
PY=/home/kiseon/.conda/envs/photochem_090_fix/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
for resv in "$@"; do
    tag=$(echo "$resv" | tr '.' 'p')
    d=$HERE/scan_pc090/heh$tag
    [ -f "$d/lower_atmosphere_profile.dat" ] && { echo "[$resv] have it"; continue; }
    mkdir -p "$d"
    ( cd "$d" && $PY $EX/src/utils/photochem_to_lower_profile.py . \
        --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 \
        --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 \
        --stellar-flux $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
        --flux-at-planet --wavelength-unit A --toa 1.0e-2 \
        --atoms H,He,N,O,C --abundances He=$resv --kzz-const 1.0e9 \
        --trial-flux-H 4768206.9882 --trial-flux-He 20433648.55 \
        --iteration 0 --abundance-tol 1.0 ) > "$d/adapter.log" 2>&1
    echo "[$resv] rc=$? $(grep -c STEADY "$d/adapter.log") $(date +%H:%M:%S)"
done

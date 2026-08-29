#!/bin/bash
# Does the chemistry+climate handoff exist at a given reservoir He/H?  The
# 15.60 arm at K_zz = 1e8 stopped inside clima's surface-temperature root
# solve, and the crossings at lower K_zz are expected to sit higher still, so
# the reachable range of the reservoir ratio is measured here on its own.
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
PY=/home/kiseon/.conda/envs/photochem_cmp/bin/python
hval=$1; kval=$2
d=$SC/clima_probe/heh${hval}_k${kval}
mkdir -p $d
( cd $d && $PY $EX/src/utils/photochem_to_lower_profile.py . \
    --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 \
    --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 \
    --stellar-flux $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
    --flux-at-planet --wavelength-unit A --toa 1.0e-2 \
    --atoms H,He,N,O,C --abundances He=$hval --kzz-const $kval \
    --trial-flux-H 4.0e6 --trial-flux-He 2.2e7 --iteration 0 \
    > probe.log 2>&1 )
rc=$?
echo "He/H=$hval K_zz=$kval rc=$rc  $(tail -n 1 $d/probe.log | cut -c1-140)"

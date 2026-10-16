#!/bin/bash
# usage: run_variant.sh <case> ; CD_* from env
D=/home/kiseon/.claude/jobs/f24ea156/tmp/condensate
case=$1
mkdir -p $D/$case && cd $D/$case || exit 1
exec env -i HOME=$HOME PATH=/opt/miniconda3/bin:/usr/local/bin:/usr/bin:/bin OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} \
  $(env | grep '^CD_' | tr '\n' ' ') \
  timeout 4h /opt/miniconda3/bin/python3 $D/condensate_driver.py . --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 --stellar-flux /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt --flux-at-planet --wavelength-unit A --atoms H,He,N,O,C --abundances He=2.09 --kzz-const 1.0e9 --column-energy-treatment nlte --blocks 1600 --trial-flux-H 13114328.6265 --trial-flux-He 9748630.636500001 --iteration 5 --conducted-heat-flux 2.3595432 --trial-flux-H2 31463504.996 > $D/$case/run.log 2>&1

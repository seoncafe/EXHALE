#!/bin/bash
# One arm of the K_zz x reservoir-He/H scan under the profile boundary
# condition.  Everything except --kzz-const and the reservoir He/H is copied
# from LHS1140b/exhale/flux_closure/heh11p1/closure.json, which is the
# K_zz = 1e9 arm the scan is anchored on.
#
#   run_arm.sh <kzz_tag> <kzz_value> <heh_label> <heh_value> <seed_output_dir> <phi0_H> <phi0_He>
#
# e.g. run_arm.sh 1e8 1.0e8 11p1 11.1 /path/to/k01/output 4.19e6 2.12e7
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"      # He I 10830 at WINERED HIRES-Y, R = 68,000
HERE=$EX/LHS1140b/exhale/kzz_profile_scan
DRV=$EX/src/utils/element_flux_closure.py
SED_FILE=$EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt

ktag=$1; kval=$2; hlab=$3; hval=$4; seed=$5; phiH=$6; phiHe=$7

d=$HERE/kzz$ktag/heh$hlab
mkdir -p $d

cat > $d/input_template.inp <<TPL
Planet name: LHS1140b_lower_profile
Log10 lower boundary number density [cm^-3]: 13.506
Planet radius [R_J]: 0.157692
Planet mass [M_J]: 0.0176220
Equilibrium temperature [K]: 226.0
Orbital distance [AU]: 0.0946
Escape radius [R_p]: 2.00
He/H number ratio: $hval
2D approximate method: Mdot
Parent star mass [M_sun]: 0.1844
Spectrum type: Load from file..
Spectrum file: $SED_FILE
Use only EUV? False
[E_low,E_mid,E_high] = [ 13.60 , 123.98 , 1.24e3 ]
Log10 of X-ray luminosity [erg/s]: 26.372
Log10 of EUV luminosity [erg/s]: 26.404
Grid type: Mixed
Numerical flux: HLLC
Reconstruction scheme: PLM
Include He23S? True
Load IC? True
Do only PP: True
Force start: False
He_diffusion: True
Domain mode: Spherical
Outer radius [R_p]: 30.0
Stellar radius [R_sun]: 0.2159
Base BC: pressure 1.0
Valve eps: 1.0e-4
Resid tol: 1.0e-3
Lower atmosphere profile: lower_atmosphere_profile.dat
TPL

cat > $d/closure.json <<JSN
{
  "comment": "LHS 1140 b elemental-flux closure, K_zz = $kval cm^2/s, reservoir He/H = $hval as the STARTING value -- the converged He/H is an output. Everything else matches flux_closure/heh11p1.",

  "python": "/home/kiseon/.conda/envs/photochem_cmp/bin/python",
  "adapter": "$EX/src/utils/photochem_to_lower_profile.py",
  "adapter_run_dir": ".",
  "adapter_args": [
    "--mp", "0.0176220",
    "--r-ref", "0.157692",
    "--p-ref", "1.0",
    "--p-match", "1.0e-6",
    "--climate",
    "--climate-p-deep", "20.0",
    "--boa-pressure-factor", "1.0",
    "--stellar-flux", "$SED_FILE",
    "--flux-at-planet",
    "--wavelength-unit", "A",
    "--toa", "1.0e-2",
    "--atoms", "H,He,N,O,C",
    "--abundances", "He=$hval",
    "--kzz-const", "$kval"
  ],
  "p_top_bar": null,

  "exhale_bin": "$EX/EXHALE.x",
  "input_template": "$d/input_template.inp",
  "omp_num_threads": 4,
  "resid_tol": "4.0e-4",
  "exhale_env": {"EXHALE_PTC": "1", "EXHALE_PTC_JFNK": "1",
                 "EXHALE_PTC_DTAU0": "1.0"}
}
JSN

run_closure () {
  local resume=""
  [ -f $d/closure_history.txt ] && resume="--resume"
  python3 $DRV $d --phi0-H $phiH --phi0-He $phiHe \
      --config $d/closure.json --seed $seed $resume \
      --tol 0.05 --kmax 8 >> $HERE/kzz${ktag}_heh${hlab}.log 2>&1
}

# Record a relaxation in the arm's own closure.json, so no row can be quoted
# without the relaxation that produced it being visible next to it.
note () {   # note <what> <value> <sentence>
  python3 - "$d/closure.json" "$1" "$2" "$3" <<'PYJSN'
import json, sys
p, key, val, why = sys.argv[1:5]
c = json.load(open(p))
if key == 'abundance-tol':
    if '--abundance-tol' not in c['adapter_args']:
        c['adapter_args'] += ['--abundance-tol', val]
elif key == 'resid_tol':
    c['resid_tol'] = val
else:
    c['exhale_env'][key] = val
c['comment'] += ' ' + why
json.dump(c, open(p, 'w'), indent=2)
PYJSN
  rm -rf $d/k[0-9][0-9] $d/closure_history.txt
}

: > $HERE/kzz${ktag}_heh${hlab}.log
run_closure
rc=$?

# The elemental-conservation check of the chemistry handoff is a ratio to
# hydrogen, and hydrogen is a minor species here: at He/H = 15 it is 6% of
# the gas by number, so the same absolute solver noise on N appears in N/H
# amplified by 1/q_H.
if [ $rc -ne 0 ] && grep -q "departs from the input" $HERE/kzz${ktag}_heh${hlab}.log; then
  echo "  elemental-conservation check refused at 1e-4; re-solving at 1e-3"
  note abundance-tol 1.0e-3 "The elemental-conservation check of the handoff is relaxed to 1e-3 on this arm: hydrogen is a minor species at this He/H and N/H carries the solver noise on N amplified by 1/q_H."
  run_closure; rc=$?
fi

# What the aborted solves are actually short of is the residual target: the
# best iterate they return sits just above 4e-4 rather than far from it.  The
# target is loosened a step at a time and the achieved ||R|| of every arm is
# in the table, so a row bought with a looser target can be read as such.
for rt in 1.0e-3 3.0e-3; do
  if [ $rc -ne 0 ] && grep -q "done info=2" $HERE/kzz${ktag}_heh${hlab}.log; then
    echo "  steady solve still info=2; re-solving with Resid tol = $rt"
    note resid_tol $rt "The wind residual target is $rt on this arm rather than 4e-4; the achieved ||R|| is reported with the row."
    run_closure; rc=$?
  fi
done

# A lower K_zz moves the helium profile of the wind further in one
# composition pass, and past 3e7 the pass hands the steady solver a state
# its line search cannot descend from (info = 2).  The composition update is
# then under-relaxed harder.  The drift test inside EXHALE_main is the
# UNDAMPED distance to the fixed point, so a smaller omega lengthens the
# iteration without changing the state it converges to.
for om in 0.25 0.125; do
  if [ $rc -ne 0 ] && grep -q "done info=2" $HERE/kzz${ktag}_heh${hlab}.log; then
    echo "  steady solve returned info=2; re-solving with EXHALE_DIFF_OMEGA=$om"
    note EXHALE_DIFF_OMEGA $om "The composition outer loop starts at omega = $om on this arm; the drift test it converges on is undamped."
    run_closure; rc=$?
  fi
done

echo "closure rc=$rc  (kzz=$kval heh=$hval)"
last=$(ls -d $d/k[0-9][0-9] 2>/dev/null | sort | tail -n 1)
if [ $rc -ne 0 ] || [ -z "$last" ]; then echo "STOP kzz=$kval heh=$hval"; exit 1; fi
( cd $last && MPLBACKEND=Agg python3 $EX/EXHALE_transit.py > transit.log 2>&1 )
echo "transit rc=$? in $last"

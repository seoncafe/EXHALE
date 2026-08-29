#!/bin/bash
# One point of the LHS 1140 b eddy-diffusion scan: build the case directory
# from the Kzz = 0 control, restart from the control's converged state, take
# the direct-steady (JFNK) route, then the post-processing pass and the
# transit synthesis.  Everything except "He_Kzz" is the control's setup.
# usage: ./run_kzz_case.sh <tag> <Kzz value>
set -u
tag=$1; kzz=$2
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(dirname "$0"); cd "$here" || exit 1
ctrl=heh0p55_diff_ctrl
d=$tag

mkdir -p "$d/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    "$ctrl/input.inp" > "$d/input.inp"
echo "He_Kzz: $kzz" >> "$d/input.inp"

\cp -f "$ctrl/output/Hydro_ioniz.txt" "$d/output/Hydro_ioniz_IC.txt"
\cp -f "$ctrl/output/Ion_species.txt" "$d/output/Ion_species_IC.txt"

cd "$d" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
echo "[$tag] JFNK START $(date +%H:%M:%S)  He_Kzz=$kzz"
EXHALE_DIFFUSION_CHECK=1 EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 \
   "$BIN" > run.log 2> diffcheck.log
rc=$?
echo "[$tag] JFNK DONE rc=$rc $(grep -o 'done info=[0-9]*' run.log | tail -n 1)" \
     "$(grep -o '||R||= *[0-9.E+-]*' run.log | tail -n 1)"
grep -q "done info=0" run.log || { echo "[$tag] JFNK FAILED"; exit 1; }

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log | tail -n 1)"

MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"

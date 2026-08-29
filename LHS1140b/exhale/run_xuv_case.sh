#!/bin/bash
# One point of the LHS 1140 b XUV grid: the scalar-base diffusion
# configuration of the composition scan, with the stellar spectrum scaled by
# a fixed factor, to bracket the 2025 non-detection (Phase F item 3).
#
# Identical to run_heh_diff_case.sh in every respect except the spectrum
# file, which carries the scaling in its flux column -- read_sed recomputes
# LX and LEUV from the file, so the luminosity lines of input.inp are
# provenance only and are rewritten here to match.
#
# usage: ./run_xuv_case.sh <seed_tag> <tag> <He/H> <sed_file> [n_repeats]
#   KZZ=<value>    eddy coefficient      (default 1.0e9)
#   RESID=<value>  wind residual target  (default 5.0e-3)
set -u
seed=$1; tag=$2; heh=$3; sed_file=$4; nrep=${5:-4}
kzz=${KZZ:-1.0e9}
resid=${RESID:-5.0e-3}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(dirname "$0"); cd "$here" || exit 1
ctrl=heh0p55_diff_ctrl

mkdir -p "$tag/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s|^He/H number ratio: .*|He/H number ratio: $heh|" \
    -e "s|^Spectrum file: .*|Spectrum file: $sed_file|" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    -e "s/^Resid tol: .*/Resid tol: $resid/" \
    "$ctrl/input.inp" > "$tag/input.inp"
echo "He_Kzz: $kzz" >> "$tag/input.inp"
\cp -f "$seed/output/Hydro_ioniz.txt" "$tag/output/Hydro_ioniz_IC.txt"
\cp -f "$seed/output/Ion_species.txt" "$tag/output/Ion_species_IC.txt"

cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
echo "[$tag] START from $seed, He/H=$heh, SED=$sed_file  $(date +%H:%M:%S)"
for i in $(seq 1 "$nrep"); do
   [ $i -gt 1 ] && { \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
                     \cp -f output/Ion_species.txt output/Ion_species_IC.txt; }
   sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
   EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > "run$i.log" 2>&1
   echo "[$tag] pass-set $i: $(grep -o 'done info=[0-9]*' run$i.log | tail -n 1)" \
        "$(grep -o '||R||= *[0-9.E+-]*' run$i.log | tail -n 1)" \
        "$(grep -o 'composition drift = *[0-9.E+-]*' run$i.log | tail -n 1)"
done

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"

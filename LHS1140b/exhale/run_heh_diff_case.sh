#!/bin/bash
# One point of the LHS 1140 b composition scan *with* the Phase-D binary
# diffusion operator active at the adopted eddy coefficient He_Kzz = 1e9.
#
# The configuration is the one of heh0p55_diff_ctrl in every respect except
# the He/H reservoir ratio, the eddy coefficient and the wind residual
# target.  The case restarts from the converged state of the seed case named
# on the command line -- the same-composition diffusion-off solution where
# one exists, otherwise the neighbouring composition already solved with
# diffusion on -- takes the direct-steady (JFNK) route repeatedly until the
# composition drift stops moving, then one "Do only PP" pass and the transit
# synthesis.
#
# usage: ./run_heh_diff_case.sh <seed_tag> <tag> <He/H> [n_repeats]
#   KZZ=<value>    eddy coefficient           (default 1.0e9)
#   RESID=<value>  wind residual target       (default 5.0e-3)
set -u
seed=$1; tag=$2; heh=$3; nrep=${4:-6}
kzz=${KZZ:-1.0e9}
resid=${RESID:-5.0e-3}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(dirname "$0"); cd "$here" || exit 1
ctrl=heh0p55_diff_ctrl

mkdir -p "$tag/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s|^He/H number ratio: .*|He/H number ratio: $heh|" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    -e "s/^Resid tol: .*/Resid tol: $resid/" \
    "$ctrl/input.inp" > "$tag/input.inp"
echo "He_Kzz: $kzz" >> "$tag/input.inp"
\cp -f "$seed/output/Hydro_ioniz.txt" "$tag/output/Hydro_ioniz_IC.txt"
\cp -f "$seed/output/Ion_species.txt" "$tag/output/Ion_species_IC.txt"

cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
echo "[$tag] START from $seed, He/H=$heh, He_Kzz=$kzz, Resid tol=$resid  $(date +%H:%M:%S)"
for i in $(seq 1 "$nrep"); do
   [ $i -gt 1 ] && { \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
                     \cp -f output/Ion_species.txt output/Ion_species_IC.txt; }
   sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
   EXHALE_DIFFUSION_CHECK=1 EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 \
      "$BIN" > "run$i.log" 2> diffcheck.log
   echo "[$tag] pass-set $i: $(grep -o 'done info=[0-9]*' run$i.log | tail -n 1)" \
        "$(grep -o '||R||= *[0-9.E+-]*' run$i.log | tail -n 1)" \
        "passes=$(grep -c 'outer pass' run$i.log)" \
        "$(grep -o 'composition drift = *[0-9.E+-]*' run$i.log | tail -n 1)"
done

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"

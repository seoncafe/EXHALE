#!/bin/bash
# Walk the eddy-diffusion scan upward, each case restarted from the converged
# solution of the case one decade below it, and with a wind residual target
# the solver can actually meet on this configuration.
#
# Two things stop the large-K_zz cases otherwise.  Restarting them from the
# K_zz = 0 control asks the wind solver to absorb the whole helium profile in
# one step; one decade at a time is a step it takes.  And on this wind the
# JFNK line search settles at ||R|| ~ 1.3e-3 rather than the 1.0e-3 the
# control reaches, which the composition outer loop reads as a solver failure
# and stops on after a single pass -- so the residual target is stated per
# case through RESID (default 3.0e-3) and the achieved ||R|| is reported with
# every row rather than assumed.
#
# usage: ./ladder_kzz.sh <seed_tag> <tag> <Kzz> [n_repeats]
#   RESID=<value>  wind residual target written into input.inp
set -u
seed=$1; tag=$2; kzz=$3; nrep=${4:-6}
resid=${RESID:-3.0e-3}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(dirname "$0"); cd "$here" || exit 1
ctrl=heh0p55_diff_ctrl

# keep the attempt that started from the control, if there is one
if [ -d "$tag" ] && [ ! -d "$tag/attempt_seed_ctrl" ]; then
   mkdir -p "$tag/attempt_seed_ctrl"
   for f in "$tag"/*.log "$tag"/*.txt "$tag"/input.inp "$tag"/output; do
      [ -e "$f" ] && \cp -rf "$f" "$tag/attempt_seed_ctrl/"
   done
fi

mkdir -p "$tag/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    -e "s/^Resid tol: .*/Resid tol: $resid/" \
    "$ctrl/input.inp" > "$tag/input.inp"
echo "He_Kzz: $kzz" >> "$tag/input.inp"
\cp -f "$seed/output/Hydro_ioniz.txt" "$tag/output/Hydro_ioniz_IC.txt"
\cp -f "$seed/output/Ion_species.txt" "$tag/output/Ion_species_IC.txt"

cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
echo "[$tag] LADDER START from $seed, He_Kzz=$kzz, Resid tol=$resid  $(date +%H:%M:%S)"
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

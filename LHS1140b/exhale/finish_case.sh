#!/bin/bash
# Finish one LHS 1140b case from its marched snapshot:
#   JFNK (EXHALE_PTC) -> post-processing pass -> EXHALE_transit.py.
# usage: ./finish_case.sh <case> [seed_dir]
# With seed_dir, the IC is that directory's converged output (load_IC carries
# it onto this case's He/H); otherwise the case's own marched snapshot.
set -u
c=$1; seed=${2:-}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd "$(dirname "$0")/$c" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

[ -f run.log ] && mv run.log run_march.log
src=${seed:+../$seed/output}; src=${src:-output}
\cp -f $src/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f $src/Ion_species.txt output/Ion_species_IC.txt

# solver configuration: Load IC + PLM + Valve/Resid, no marching keys
sed -i 's/^Load IC? False$/Load IC? True/;
        s/^Reconstruction scheme: .*/Reconstruction scheme: PLM/;
        s/^du_th \[PLM,WENO3\]: .*/Valve eps: 1.0e-4/;
        s/^Solver: Newton$/Resid tol: 1.0e-3/;
        /^IC mode: auto$/d; s/^Do only PP: True$/Do only PP: False/' input.inp

echo "[$c] JFNK START $(date +%H:%M:%S) ${seed:+(seeded from $seed)}"
EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$EX/EXHALE.x" > run.log 2>&1
rc=$?
info=$(grep -o "done info=[0-9]*" run.log | tail -n 1)
echo "[$c] JFNK DONE rc=$rc $info  $(grep -o '||R||= *[0-9.E+-]*' run.log | tail -n 1)"
grep -q "done info=0" run.log || { echo "[$c] JFNK FAILED"; exit 1; }

# post-processing pass on the solved state
\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$EX/EXHALE.x" > pp.log 2>&1
echo "[$c] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log | tail -n 1)"

MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$c] TRANSIT DONE rc=$?"
grep -i 'He 10830 metrics' transit.log | tail -n 1
echo "[$c] CASE DONE $(date +%H:%M:%S)"

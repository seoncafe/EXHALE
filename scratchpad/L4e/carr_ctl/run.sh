#!/bin/bash
# Reload of the carrier reload on the MATCHED Model A configuration (the fixture
# of docs/koskinen2022_model_a_comparison.tex section 5, 2026-09-12). Runs in place: EXHALE reads
# ./input.inp and writes ./output/. The IC/ pair is copied into output/ first,
# so the run enters the Newton from the pinned state without marching.
# Usage: cd backup/regression/carrier_elem_newton && ./run.sh [binary]
cd "$(dirname "$0")"
exe=${1:-../../../EXHALE.x}
mkdir -p output
\cp -f IC/Hydro_ioniz_IC.txt IC/Ion_species_IC.txt output/
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$exe" > run.log 2>&1
echo "exit=$?"

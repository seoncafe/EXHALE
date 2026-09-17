#!/bin/bash
# Reproducer of the atomic element-row coupled solve that produces no Krylov
# direction (ISSUES_20260909 section 3.1 d; PLAN_20260909_rev1 item N7).
# Runs in place: EXHALE reads ./input.inp and writes ./output/.
# Usage: cd backup/regression/atomic_elem_newton && ./run.sh [binary]
cd "$(dirname "$0")"
exe=${1:-../../../EXHALE.x}
OMP_NUM_THREADS=${OMP_NUM_THREADS:-16} "$exe" > run.log 2>&1
echo "exit=$?"

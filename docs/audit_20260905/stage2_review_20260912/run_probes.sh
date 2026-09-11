#!/usr/bin/env bash
# Reproduce the diagnostic measurements without changing production products.
set -euo pipefail
audit_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "$audit_dir/../../.." && pwd)
audit_build=${1:?Supply the private build root containing build/*.o}
audit_objects=()
for obj in "$audit_build"/build/*.o; do
  if [[ "${obj##*/}" != EXHALE_main.o ]]; then audit_objects+=("$obj"); fi
done
for probe in carrier_state_contract_probe molecular_wind_state_probe stationary_retry_history_probe; do
  extra=()
  if [[ "$probe" != molecular_wind_state_probe ]]; then
    extra+=("$project_dir/src/tests/physics_probe/assertion_report.f90")
  fi
  gfortran -O0 -g -fcheck=all -fopenmp -ffree-line-length-none \
    -I"$audit_build/build" -J"$audit_build/build" \
    "${extra[@]}" "$audit_dir/$probe.f90" "${audit_objects[@]}" \
    -L/opt/miniconda3/lib -lopenblas -Wl,-rpath,/opt/miniconda3/lib \
    -o "$audit_build/$probe.x"
done
run_dir=$(mktemp -d "$audit_dir/current_probe_run.XXXXXX")
mkdir "$run_dir/output"
cp "$project_dir/backup/regression/carrier_elem_newton/input.inp" "$run_dir/"
for name in metals.inp base.inp; do
  if [[ -f "$project_dir/backup/regression/carrier_elem_newton/$name" ]]; then
    cp "$project_dir/backup/regression/carrier_elem_newton/$name" "$run_dir/"
  fi
done
cp "$project_dir/backup/regression/carrier_elem_newton/IC/Hydro_ioniz_IC.txt" "$run_dir/output/"
cp "$project_dir/backup/regression/carrier_elem_newton/IC/Ion_species_IC.txt" "$run_dir/output/"
(
  git -C "$project_dir" rev-parse HEAD
  sha256sum "$audit_build"/*probe.x "$audit_dir"/*probe.f90
  git -C "$project_dir" diff --stat -- src
) > "$run_dir/provenance.log"
(
  cd "$project_dir"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$audit_build/carrier_state_contract_probe.x" \
    > "$run_dir/carrier_state_contract_probe.log" 2>&1
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$audit_build/stationary_retry_history_probe.x" \
    > "$run_dir/stationary_retry_history_probe.log" 2>&1
)
(
  cd "$run_dir"
  for env_name in ${!EXHALE_@}; do unset "$env_name"; done
  export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE
  set +e
  timeout 180s "$audit_build/molecular_wind_state_probe.x" partition 1 40 0.5 0.01 \
    > molecular_wind_state_probe.log 2>&1
  result=$?
  printf 'exit_status=%s\n' "$result" > exit_status.log
  exit "$result"
)
printf 'Measurements: %s\n' "$run_dir"

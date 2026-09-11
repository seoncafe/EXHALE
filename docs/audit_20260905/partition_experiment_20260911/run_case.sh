#!/usr/bin/env bash
# Run a diagnostic in a fresh directory, preserving the original fixture.
set -euo pipefail
experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "$experiment_dir/../../.." && pwd)
experiment_exe=${1:?Supply the absolute path to transport_wind_experiment.x}
fixture=${2:?Supply atomic_elem_newton or carrier_elem_newton}
mode=${3:-partition}
passes=${4:-3}
maxit=${5:-40}
omega=${6:-0.5}
trust=${7:-0.01}
experiment_threads=${8:-1}
case "$fixture" in
  atomic_elem_newton|carrier_elem_newton) ;;
  *) echo 'Unsupported fixture' >&2; exit 2 ;;
esac
mkdir -p "$experiment_dir/runs"
run_dir=$(mktemp -d "$experiment_dir/runs/${fixture}_${mode}.XXXXXX")
mkdir "$run_dir/output"
cp "$project_dir/backup/regression/$fixture/input.inp" "$run_dir/input.inp"
# This mechanical input change selects the saved state; the standalone
# driver performs no physical time steps before its first measurement.
sed -i 's/Load IC? False/Load IC? True/' "$run_dir/input.inp"
for input_file in metals.inp base.inp; do
  if [[ -f "$project_dir/backup/regression/$fixture/$input_file" ]]; then
    cp "$project_dir/backup/regression/$fixture/$input_file" "$run_dir/"
  fi
done
cp "$project_dir/backup/regression/$fixture/IC/Hydro_ioniz_IC.txt" "$run_dir/output/"
cp "$project_dir/backup/regression/$fixture/IC/Ion_species_IC.txt" "$run_dir/output/"
cp "$experiment_dir/transport_wind_experiment.f90" "$run_dir/driver_snapshot.f90"
(
  git -C "$project_dir" rev-parse HEAD
  sha256sum "$experiment_exe" "$run_dir/driver_snapshot.f90" "$run_dir/input.inp" \
    "$run_dir/output/Hydro_ioniz_IC.txt" "$run_dir/output/Ion_species_IC.txt"
  printf 'mode=%s passes=%s maxit=%s omega=%s trust=%s threads=%s\n' \
    "$mode" "$passes" "$maxit" "$omega" "$trust" "$experiment_threads"
) > "$run_dir/provenance.log"
printf 'Run directory: %s\n' "$run_dir"
(
  cd "$run_dir"
  # Disable inherited diagnostic overrides. Only this experiment's arguments
  # set its budgets and OpenMP thread count; BLAS uses one thread.
  for env_name in ${!EXHALE_@}; do unset "$env_name"; done
  export OMP_NUM_THREADS="$experiment_threads" OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE
  set +e
  timeout 180s "$experiment_exe" "$mode" "$passes" "$maxit" "$omega" "$trust" > run.log 2>&1
  result=$?
  printf 'exit_status=%s\n' "$result" > exit_status.log
  exit "$result"
)
printf 'Completed: %s\n' "$run_dir"

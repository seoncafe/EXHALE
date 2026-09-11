#!/usr/bin/env bash
# Build current production modules and a separate diagnostic driver.
set -euo pipefail
experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "$experiment_dir/../../.." && pwd)
build_root=${1:-$(mktemp -d /tmp/exhale_partition_20260911.XXXXXX)}
cd "$project_dir"
make -j4 OBJDIR="$build_root/build" EXE="$build_root/EXHALE_control.x"
experiment_objects=()
for obj in "$build_root"/build/*.o; do
  if [[ "${obj##*/}" != EXHALE_main.o ]]; then experiment_objects+=("$obj"); fi
done
gfortran -O3 -fopenmp -ffree-line-length-none -fcheck=bounds \
  -I"$build_root/build" -J"$build_root/build" \
  "$experiment_dir/transport_wind_experiment.f90" "${experiment_objects[@]}" \
  -L/opt/miniconda3/lib -lopenblas -Wl,-rpath,/opt/miniconda3/lib \
  -o "$build_root/transport_wind_experiment.x"
printf 'Diagnostic executable: %s/transport_wind_experiment.x\n' "$build_root"

#!/usr/bin/env bash
# Compile only the production modules exercised by this review.
# All compiler products stay in a newly created temporary directory.
set -euo pipefail
review_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "$review_dir/../.." && pwd)
review_build=$(mktemp -d /tmp/exhale_rev1_review_XXXXXX)
cd -- "$project_dir"
gfortran -O0 -fopenmp -fcheck=all -ffree-line-length-none \
    -J"$review_build" -I"$review_build" \
    src/modules/init/parameters.f90 \
    src/modules/init/species_table.f90 \
    src/modules/radiation/charge_exchange.f90 \
    "$review_dir/plan_20260918_rev1_source_probe.f90" \
    -o "$review_build/source_probe"
{
    gfortran --version | head -n 1
    "$review_build/source_probe"
    python3 "$review_dir/plan_20260918_rev1_workflow_checks.py"
} | tee "$review_dir/plan_20260918_rev1_review_checks.log"

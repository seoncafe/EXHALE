#!/usr/bin/env bash
# Tests of src/utils/map_state_to_grid.py (item S3 of PLAN_20260913):
# identity mapping, a shifted grid round trip, the refusals (mismatched
# species radii, out-of-support physical targets, non-monotone radii), the
# ghost extrapolation, the He/H invariance and the seed metadata. One line
# per assertion: PASS|FAIL <name> measured= reference= tol=
set -u
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
work="${EXHALE_TEST_OUT:-$root/build/tests/state_mapper}"
mkdir -p "$work"
python3 "$here/state_mapper_tests.py" "$root/src/utils/map_state_to_grid.py" "$work"
status=$?
if [ "$status" -ne 0 ]; then echo "state_mapper: FAILED"; else echo "state_mapper: PASSED"; fi
exit "$status"

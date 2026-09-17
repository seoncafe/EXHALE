#!/usr/bin/env bash
# Tests of src/utils/map_state_to_grid.py (item S3 of PLAN_20260913):
# identity mapping, a shifted grid round trip, the refusals (mismatched
# species radii, out-of-support physical targets, non-monotone radii), the
# ghost extrapolation, the He/H invariance and the seed metadata; and the
# reservoir rescaling --reservoir <El>/H <value> (the element's columns
# multiplied by the ratio of the two reservoirs, the hydrogen, pressure and
# velocity columns left as the option-free mapping writes them, the density
# and temperature following the new particle count, the header lines, and
# the refusals when the source states no reservoir, when the value is not
# positive and when the state carries HeH+), for helium and for a metal
# element -- one element, several in one call, and the refusals of an
# element the '# reservoir' line does not name and of one the species file
# carries no column for. One line per assertion:
# PASS|FAIL <name> measured= reference= tol=
set -u
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
work="${EXHALE_TEST_OUT:-$root/build/tests/state_mapper}"
mkdir -p "$work"
python3 "$here/state_mapper_tests.py" "$root/src/utils/map_state_to_grid.py" "$work"
status=$?
if [ "$status" -ne 0 ]; then echo "state_mapper: FAILED"; else echo "state_mapper: PASSED"; fi
exit "$status"

#!/usr/bin/env bash
# FROZEN RECORD (2026-09-05): these probes call the five-argument speed_estimate_ROE of
# revision 35d9dd5. Since Phase 1 item A2 (docs/a2_roe_interface.md) the production routine
# has a different signature, so this script no longer compiles against the tree; it is kept
# as the record of what was measured at 35d9dd5 and is not rebuilt.
# Reproduce the second review's local diagnostics using production routines.
set -euo pipefail
audit_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
audit_root=$(cd -- "$audit_dir/../.." && pwd)
audit_src="$audit_root/src/modules"
audit_work=$(mktemp -d /tmp/exhale-review2-checks-XXXXXX)
cd -- "$audit_work"
audit_fc=${FC:-gfortran}
"$audit_fc" --version | head -n 1
audit_flags=(-O0 -g -fcheck=all -fbacktrace -J. -I.)
"$audit_fc" "${audit_flags[@]}" \
  "$audit_src/init/parameters.f90" \
  "$audit_src/flux/speed_estimate_ROE.f90" \
  "$audit_src/lower_atmosphere/h3p_cooling.f90" \
  "$audit_dir/roe_equal_pressure.f90" -o roe_equal_pressure.x
./roe_equal_pressure.x
"$audit_fc" "${audit_flags[@]}" \
  "$audit_src/init/parameters.f90" \
  "$audit_src/functions/cross_sec.f90" \
  "$audit_src/functions/h2_photo_channels.f90" \
  "$audit_dir/h2_detector_ratio.f90" -o h2_detector_ratio.x
./h2_detector_ratio.x
python3 "$audit_dir/photoevent_energy_check.py"
printf 'Temporary build directory: %s\n' "$audit_work"
printf '%s\n' 'Diagnostics completed; inspect values. Exit status is not a physical acceptance verdict.'

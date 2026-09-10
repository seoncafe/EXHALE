#!/usr/bin/env bash
# FROZEN RECORD (2026-09-05): these probes call the five-argument speed_estimate_ROE of
# revision 35d9dd5. Since Phase 1 item A2 (docs/a2_roe_interface.md) the production routine
# has a different signature, so this script no longer compiles against the tree; it is kept
# as the record of what was measured at 35d9dd5 and is not rebuilt.
# Trace the production ROE flux downstream of its auxiliary speed estimate.
set -euo pipefail
audit_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
audit_root=$(cd -- "$audit_dir/../.." && pwd)
audit_src="$audit_root/src/modules"
audit_work=$(mktemp -d /tmp/exhale-review3-checks-XXXXXX)
cd -- "$audit_work"
audit_fc=${FC:-gfortran}
"$audit_fc" --version | head -n 1
audit_flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J. -I.)
"$audit_fc" "${audit_flags[@]}" \
  "$audit_src/init/parameters.f90" \
  "$audit_src/init/species_table.f90" \
  "$audit_src/lower_atmosphere/molecular_infrared_data.f90" \
  "$audit_src/states/caloric_eos.f90" \
  "$audit_src/functions/UW_conversions.f90" \
  "$audit_src/flux/speed_estimate_ROE.f90" \
  "$audit_src/flux/Num_Fluxes.f90" \
  "$audit_dir/roe_flux_vacuum_probe.f90" -o roe_flux_vacuum_probe.x
./roe_flux_vacuum_probe.x
python3 "$audit_dir/review3_ledger_checks.py"
printf 'Temporary build directory: %s\n' "$audit_work"
printf '%s\n' 'Diagnostics completed; this is not a full RK stage or a run of the proposed revision.'

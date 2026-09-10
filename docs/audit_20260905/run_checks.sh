#!/usr/bin/env bash
# FROZEN RECORD (2026-09-05): these probes call the five-argument speed_estimate_ROE of
# revision 35d9dd5. Since Phase 1 item A2 (docs/a2_roe_interface.md) the production routine
# has a different signature, so this script no longer compiles against the tree; it is kept
# as the record of what was measured at 35d9dd5 and is not rebuilt.
# Compile the inspected production routines without changing the project build.
set -euo pipefail
audit_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
audit_root=$(cd -- "$audit_dir/../.." && pwd)
audit_src="$audit_root/src/modules"
audit_work=$(mktemp -d /tmp/exhale-physics-checks-XXXXXX)
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
  "$audit_src/lower_atmosphere/h3p_cooling.f90" \
  "$audit_src/lower_atmosphere/mol_rates.f90" \
  "$audit_src/lower_atmosphere/molecular_infrared_cooling.f90" \
  "$audit_src/functions/cross_sec.f90" \
  "$audit_src/functions/h2_photo_channels.f90" \
  "$audit_dir/audit_probe.f90" -o audit_probe.x
./audit_probe.x > probe.out
"$audit_fc" "${audit_flags[@]}" \
  "$audit_src/init/parameters.f90" \
  "$audit_src/functions/cross_sec.f90" \
  "$audit_src/functions/h2_photo_channels.f90" \
  "$audit_root/src/tests/e1_h2/e1_h2_channel_check.f90" -o e1.x
./e1.x > e1.out
"$audit_fc" "${audit_flags[@]}" \
  "$audit_src/lower_atmosphere/oxygen_rates.f90" \
  "$audit_src/lower_atmosphere/h2_self_shielding_table.f90" \
  "$audit_src/lower_atmosphere/lyman_werner.f90" \
  "$audit_src/lower_atmosphere/water_photolysis.f90" \
  "$audit_root/src/tests/a2_m2/a2_m2_kinetics_check.f90" -o oxygen.x
./oxygen.x > oxygen.out
printf 'Diagnostic outputs: %s\n' "$audit_work"
# These are diagnostic reports, not pass/fail acceptance tests. Inspect values.

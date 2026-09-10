#!/usr/bin/env bash
# Build and run the physics probe tests.
#
# Every driver links the PRODUCTION sources under src/modules, never a copy
# of them: source_closure.py follows the `use` statements of each driver,
# and the closure it prints is what gets compiled, in an order in which each
# file follows the modules it uses.  Objects and module files go to
# build/tests/physics_probe/, which this script creates.
#
# One line per assertion:  PASS|FAIL <name> measured= reference= tol=
# Exit status is nonzero if any assertion failed or any build failed.
#
# EXHALE_OBJDIR selects another build (a private OBJDIR, as a concurrent
# item's build uses).  The drivers compile the production sources themselves,
# so the only file this suite takes from a build tree is the generated
# build_stamp.f90, and that is the one it then reads.  EXHALE_TEST_OUT
# selects where the objects, module files and executables are written;
# EXHALE_TEST_WORK, the name the drivers read for their scratch directory,
# selects the same thing.  No test here runs EXHALE.x, so EXHALE_EXE has no
# effect on this suite.
set -uo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../../.." && pwd)
objdir="${EXHALE_OBJDIR:-$root/build}"
# A private object directory, so that two builds of this suite in one tree
# do not overwrite each other.
work="${EXHALE_TEST_OUT:-${EXHALE_TEST_WORK:-$root/build/tests/physics_probe}}"
mkdir -p "$work"

FC=${FC:-gfortran}
flags=(-O0 -g -fcheck=all -fbacktrace -fopenmp -J"$work" -I"$work")

drivers=(
  riemann_wave_speeds.f90
  weno3_reconstruction_order.f90
  positivity_limiter_scaling.f90
  h2_channel_detector_ratio.f90
  h3p_cooling_limits.f90
  h2_rovibrational_identity.f90
  h2_infrared_line_populations.f90
  photoevent_energy_ledger.f90
  wind_ae_bridge_composition.f90
  opacity_model_p_ledger.f90
  lya_band_transmission.f90
  lya_beam_cell_mean.f90
  lya_escape_probability.f90
  fine_structure_escape_probability.f90
  lyman_werner_cell_mean.f90
  helium_level_cooling_ledger.f90
  voronov_carbon_ionization.f90
  atomic_mass_and_radius_constants.f90
  ionization_threshold_turn_on.f90
  xuv_cell_mean_attenuation.f90
  photoionization_field_substitution.f90
  photoionization_field_self_consistency.f90
  absorbed_fraction_switch.f90
  species_formation_energy_table.f90
  heating_channel_closure.f90
  co_destruction.f90
  co_helium_ion_sink.f90
  h3p_recombination_branching.f90
  water_photolysis_lyman_alpha_yields.f90
  metal_photoionization_fits.f90
  recombination_coefficient_fits.f90
)

# utilities.f90 reads the revision strings of the build_stamp.f90 that the
# Makefile generates in the object directory.  Without a build tree to read
# one from, write a stamp that says so; nothing tested here depends on its
# contents.
stamp="$objdir/build_stamp.f90"
if [ ! -f "$stamp" ]; then
  stamp="$work/build_stamp.f90"
  cat > "$stamp" <<'STAMP'
      module build_stamp
      ! Written by src/tests/physics_probe/run.sh when the Makefile has not
      ! generated one. The tests do not read these strings.
      implicit none
      character(len=*), parameter :: build_git   = 'tests'
      character(len=*), parameter :: build_dirty = 'unknown'
      end module build_stamp
STAMP
fi

driver_paths=()
for d in "${drivers[@]}"; do driver_paths+=("$here/$d"); done

order=$(python3 "$here/source_closure.py" "$root" "${driver_paths[@]}") || { echo "FAIL source_closure"; exit 1; }
# source_closure.py looks for the stamp in build/; compile the one selected
# above instead, first, since it uses no other module.
order=$(printf '%s\n' "$order" | grep -v '/build_stamp\.f90$')
order="$stamp
$order"

objects=()
status=0
for src in $order; do
  obj="$work/$(basename "${src%.f90}").o"
  if ! "$FC" "${flags[@]}" -c "$src" -o "$obj"; then
    echo "FAIL compile $(basename "$src")"
    exit 1
  fi
  case " ${driver_paths[*]} " in
    *" $src "*) ;;
    *) objects+=("$obj") ;;
  esac
done

for d in "${drivers[@]}"; do
  exe="$work/${d%.f90}.x"
  if ! "$FC" "${flags[@]}" "$work/${d%.f90}.o" "${objects[@]}" -o "$exe"; then
    echo "FAIL link ${d%.f90}"
    status=1
    continue
  fi
  EXHALE_TEST_WORK="$work" "$exe" || status=1
done

python3 "$here/radiation_exchange_ledger.py" || status=1
python3 "$here/jupiter_radius_python_tools.py" || status=1
python3 "$here/threshold_literal_uniqueness.py" || status=1
python3 "$here/metal_photoion_table_transcription.py" || status=1
python3 "$here/h2_partition_sum_uniqueness.py" || status=1
python3 "$here/heating_sum_uniqueness.py" || status=1
EXHALE_TEST_WORK="$work" python3 "$here/species_formation_energy_python_table.py" || status=1

if [ "$status" -ne 0 ]; then
  echo "physics_probe: FAILED"
else
  echo "physics_probe: PASSED"
fi
exit "$status"

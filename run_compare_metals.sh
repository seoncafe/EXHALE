#!/bin/bash
# ATES-metal metal on/off comparison. Metals are solved INSIDE the
# MINPACK system, with Badnell RR+DR recombination and Kingdon&Ferland
# charge exchange. Abundances are now set at RUNTIME via metals.inp, so
# the executable is built ONCE and the two cases differ only by whether
# metals.inp is present (no source edits, no rebuild between cases).
#
# Runs the same planet (./input.inp) with C+N+O on, then off, saving
#   output_metals_on/   and   output_metals_off/
#
# Usage: ./run_compare_metals.sh

set -e
PWD_DIR=$(pwd)
S=src/modules
MOD=src/mod
IR=$S/files_IO/input_read.f90

mkdir -p "$MOD"

build() {
   rm -f "$MOD"/*.mod ATES.x
   gfortran -O3 -J"$MOD" -I"$MOD" -fopenmp -ffree-line-length-none \
      $S/init/parameters.f90 $S/files_IO/metals_input_read.f90 \
      $S/files_IO/opacity_input_read.f90 $S/files_IO/input_read.f90 \
      $S/files_IO/load_IC.f90 \
      $S/files_IO/write_output.f90 $S/files_IO/write_setup_report.f90 \
      $S/functions/grav_field.f90 $S/functions/cross_sec.f90 \
      $S/radiation/opacity_models.f90 \
      $S/functions/UW_conversions.f90 $S/functions/utilities.f90 \
      $S/nonlinear_system_solver/dogleg.f90 $S/nonlinear_system_solver/enorm.f90 \
      $S/nonlinear_system_solver/hybrd1.f90 $S/nonlinear_system_solver/qform.f90 \
      $S/nonlinear_system_solver/r1mpyq.f90 $S/nonlinear_system_solver/System_HeH.f90 \
      $S/nonlinear_system_solver/System_HeHCO.f90 $S/nonlinear_system_solver/System_HeH_TR.f90 \
      $S/nonlinear_system_solver/System_implicit_adv_HeH.f90 \
      $S/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90 \
      $S/nonlinear_system_solver/dpmpar.f90 $S/nonlinear_system_solver/fdjac1.f90 \
      $S/nonlinear_system_solver/hybrd.f90 $S/nonlinear_system_solver/qrfac.f90 \
      $S/nonlinear_system_solver/r1updt.f90 $S/nonlinear_system_solver/System_H.f90 \
      $S/nonlinear_system_solver/System_implicit_adv_H.f90 \
      $S/radiation/sed_read.f90 $S/radiation/J_inc.f90 $S/radiation/Cool_coeff.f90 \
      $S/radiation/util_ion_eq.f90 $S/radiation/ionization_equilibrium.f90 \
      $S/nonlinear_system_solver/T_equation.f90 $S/post_process/post_process_adv.f90 \
      $S/states/Apply_BC.f90 $S/states/PLM_rec.f90 $S/states/Source.f90 \
      $S/states/Reconstruction.f90 $S/flux/speed_estimate_HLLC.f90 \
      $S/flux/speed_estimate_ROE.f90 $S/flux/Num_Fluxes.f90 \
      $S/time_step/RK_rhs.f90 $S/time_step/eval_dt.f90 $S/time_step/energy_semi_implicit.f90 \
      $S/init/define_grid.f90 $S/init/set_energy_vectors.f90 \
      $S/init/set_gravity_grid.f90 $S/init/set_IC.f90 $S/init/init.f90 \
      ATES_main.f90 -o ATES.x
   [ -f ATES.x ] || { echo "BUILD FAILED"; exit 1; }
}

echo "=== ATES-metal metal comparison (HD209458b) ==="

# Build once (abundances come from metals.inp at runtime).
build

# 1. metal ON (solar C, N, O via metals.inp)
echo "--- [1/2] metals ON  (C=2.69e-4, N=6.76e-5, O=4.90e-4) ---"
printf 'CI  2.69e-4\nNI  6.76e-5\nOI  4.90e-4\n' > metals.inp
mkdir -p output; ./ATES.x
rm -rf output_metals_on; mv output output_metals_on
echo "  -> output_metals_on/"

# 2. metal OFF (no metals.inp)
echo "--- [2/2] metals OFF (no metals.inp) ---"
rm -f metals.inp
mkdir -p output; ./ATES.x
rm -rf output_metals_off; mv output output_metals_off
echo "  -> output_metals_off/"

echo "Done. output_metals_on/ and output_metals_off/ written (single build)."

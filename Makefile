# =====================================================================
#  Makefile for EXHALE  (replaces the inline gfortran call that
#  run_EXHALE.sh used to perform).
#
#  Quick start:
#     make              # build ./EXHALE.x with gfortran (default)
#     make FC=ifort     # Intel classic        (or:  make ifort)
#     make FC=ifx       # Intel LLVM compiler   (or:  make ifx)
#     make -j           # parallel build
#     make clean        # remove objects + .mod files (KEEP EXHALE.x)
#     make distclean    # remove the whole build/ dir and EXHALE.x
#
#  Incremental: only the sources you edit -- and the modules that depend
#  on them -- are recompiled. The module (.o -> .o) dependency graph is
#  generated automatically by src/utils/fortdep.py into $(DEPFILE), so it
#  stays correct when sources are added or their `use` statements change.
#
#  All build artifacts (objects, .mod files, dep/compiler stamps) live in
#  $(OBJDIR)/; source basenames are unique, so a flat object dir is safe.
# =====================================================================

OBJDIR := build
MODDIR := $(OBJDIR)
EXE    := EXHALE.x

# ---- choose compiler -------------------------------------------------
# Override make's built-in default (f77), but keep any value passed on
# the command line, e.g.  make FC=ifort
ifeq ($(origin FC),default)
  FC := gfortran
endif

# ---- compiler-specific optimization / OpenMP / module-path flags -----
# (override with e.g.  make FFLAGS='-O2 -g')
ifneq (,$(filter ifort ifx,$(FC)))
  FFLAGS  := -O3 -xHost -qopenmp -no-wrap-margin
  MODFLAG := -module $(MODDIR)
else
  FFLAGS  := -O3 -fopenmp
  MODFLAG := -J$(MODDIR) -I$(MODDIR)
endif

# LAPACK (banded LU dgbtrf/dgbtrs) for the steady-state Newton/PTC solver.
# Override with e.g.  make LDLIBS='-L/path -llapack -lblas'
LDLIBS ?= -llapack

# ---- source list (canonical build order, mirrors run_EXHALE.sh) -------
SRC := \
  src/modules/init/parameters.f90 \
  src/modules/init/species_table.f90 \
  src/modules/radiation/charge_exchange.f90 \
  src/modules/files_IO/metals_input_read.f90 \
  src/modules/files_IO/opacity_input_read.f90 \
  src/modules/files_IO/input_read.f90 \
  src/modules/files_IO/load_IC.f90 \
  src/modules/files_IO/write_output.f90 \
  src/modules/files_IO/write_setup_report.f90 \
  src/modules/functions/grav_field.f90 \
  src/modules/functions/cross_sec.f90 \
  src/modules/radiation/opacity_models.f90 \
  src/modules/functions/UW_conversions.f90 \
  src/modules/functions/utilities.f90 \
  src/modules/functions/composition.f90 \
  src/modules/functions/species_diffusion.f90 \
  src/modules/lower_atmosphere/lower_column.f90 \
  src/modules/lower_atmosphere/h3p_cooling.f90 \
  src/modules/lower_atmosphere/mol_rates.f90 \
  src/modules/time_step/steady_residual.f90 \
  src/modules/time_step/steady_newton.f90 \
  src/modules/nonlinear_system_solver/dogleg.f90 \
  src/modules/nonlinear_system_solver/enorm.f90 \
  src/modules/nonlinear_system_solver/hybrd1.f90 \
  src/modules/nonlinear_system_solver/newton_solver.f90 \
  src/modules/nonlinear_system_solver/qform.f90 \
  src/modules/nonlinear_system_solver/r1mpyq.f90 \
  src/modules/nonlinear_system_solver/System_HeH.f90 \
  src/modules/nonlinear_system_solver/System_HeH_metals.f90 \
  src/modules/nonlinear_system_solver/System_HeH_mol.f90 \
  src/modules/nonlinear_system_solver/System_HeH_TR.f90 \
  src/modules/nonlinear_system_solver/System_HeH_TR_metals.f90 \
  src/modules/nonlinear_system_solver/System_implicit_adv_HeH.f90 \
  src/modules/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90 \
  src/modules/nonlinear_system_solver/dpmpar.f90 \
  src/modules/nonlinear_system_solver/fdjac1.f90 \
  src/modules/nonlinear_system_solver/hybrd.f90 \
  src/modules/nonlinear_system_solver/qrfac.f90 \
  src/modules/nonlinear_system_solver/r1updt.f90 \
  src/modules/nonlinear_system_solver/System_H.f90 \
  src/modules/nonlinear_system_solver/System_implicit_adv_H.f90 \
  src/modules/radiation/sed_read.f90 \
  src/modules/radiation/J_inc.f90 \
  src/modules/radiation/Cool_coeff.f90 \
  src/modules/radiation/util_ion_eq.f90 \
  src/modules/radiation/ionization_equilibrium.f90 \
  src/modules/radiation/lya_rt.f90 \
  src/modules/radiation/excited_hydrogen.f90 \
  src/modules/nonlinear_system_solver/T_equation.f90 \
  src/modules/post_process/post_process_adv.f90 \
  src/modules/states/Apply_BC.f90 \
  src/modules/states/PLM_rec.f90 \
  src/modules/states/Source.f90 \
  src/modules/states/Reconstruction.f90 \
  src/modules/flux/speed_estimate_HLLC.f90 \
  src/modules/flux/speed_estimate_ROE.f90 \
  src/modules/flux/Num_Fluxes.f90 \
  src/modules/time_step/RK_rhs.f90 \
  src/modules/time_step/eval_dt.f90 \
  src/modules/time_step/energy_semi_implicit.f90 \
  src/modules/init/define_grid.f90 \
  src/modules/init/set_energy_vectors.f90 \
  src/modules/init/set_gravity_grid.f90 \
  src/modules/init/set_IC.f90 \
  src/modules/init/init.f90 \
  $(wildcard src/modules/wind_ae/wae_*.f90) \
  src/EXHALE_main.f90

# objects (flat in $(OBJDIR)); let make find the sources in their subdirs
OBJ     := $(addprefix $(OBJDIR)/,$(notdir $(SRC:.f90=.o)))
DEPFILE := $(OBJDIR)/.deps.mk
vpath %.f90 $(sort $(dir $(SRC)))

# ---- Wind-AE IC generator (standalone tree; not built by `all`) ------
# wind_ae_ic.x reads an EXHALE input.inp, solves the Murray-Clay/Broome
# Wind-AE wind (ported to Fortran under src/modules/wind_ae/), and writes
# EXHALE Load-IC files. The wae_* modules are self-contained (no ATES
# module deps); inter-module order is auto-resolved by fortdep.py. Build
# with:  make wind_ae_ic.x
WAE_DIR := src/modules/wind_ae
# standalone build EXCLUDES wae_exhale_bridge.f90 (it alone uses EXHALE's
# global_parameters; it is linked only into EXHALE.x for "IC mode: windae").
WAE_SRC := $(filter-out $(WAE_DIR)/wae_exhale_bridge.f90, \
             $(sort $(wildcard $(WAE_DIR)/wae_*.f90))) $(WAE_DIR)/wind_ae_ic.f90
WAE_OBJ := $(addprefix $(OBJDIR)/,$(notdir $(WAE_SRC:.f90=.o)))
WAE_EXE := wind_ae_ic.x
WAE_DEPFILE := $(OBJDIR)/.deps_wae.mk
vpath %.f90 $(WAE_DIR)

# Rebuild everything when the compiler changes: the stamp file name
# encodes $(FC), so a different compiler makes the previous objects stale.
CSTAMP  := $(OBJDIR)/.compiler-$(FC)

# ---------------------------------------------------------------------
.PHONY: all clean distclean ifort ifx wind_ae_ic
all: $(EXE)
wind_ae_ic: $(WAE_EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(MODFLAG) $(OBJ) -o $@ $(LDLIBS)
	@echo "built $@"

# Wind-AE IC generator (separate executable; no LAPACK needed)
$(WAE_EXE): $(WAE_OBJ)
	$(FC) $(FFLAGS) $(MODFLAG) $(WAE_OBJ) -o $@
	@echo "built $@"
$(WAE_OBJ): $(CSTAMP)

# compile each source to $(OBJDIR)/<base>.o (also writes its .mod there)
$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# compiler-change stamp: a normal prerequisite of every object so that
# switching compilers (the stamp file name changes) forces a full rebuild
$(OBJ): $(CSTAMP)
$(CSTAMP): | $(OBJDIR)
	@rm -f $(OBJDIR)/.compiler-* && touch $@

$(OBJDIR):
	@mkdir -p $@

# convenience aliases matching the run_EXHALE.sh flags
ifort: ; @$(MAKE) --no-print-directory FC=ifort
ifx:   ; @$(MAKE) --no-print-directory FC=ifx

clean:
	@rm -f $(OBJ) $(WAE_OBJ) $(MODDIR)/*.mod
	@echo "cleaned objects and .mod files in $(OBJDIR)/ (kept $(EXE), $(WAE_EXE))"

distclean:
	@rm -rf $(OBJDIR) $(EXE) $(WAE_EXE)
	@echo "removed $(OBJDIR)/, $(EXE) and $(WAE_EXE)"

# ---- auto-generated module dependencies (skip when only cleaning) ---
$(DEPFILE): $(SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(SRC) > $@

$(WAE_DEPFILE): $(WAE_SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(WAE_SRC) > $@

ifeq ($(filter clean distclean,$(MAKECMDGOALS)),)
-include $(DEPFILE)
ifneq ($(filter wind_ae_ic wind_ae_ic.x,$(MAKECMDGOALS)),)
-include $(WAE_DEPFILE)
endif
endif

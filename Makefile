# =====================================================================
#  Makefile for EXHALE  (replaces the inline gfortran call that
#  run_EXHALE.sh used to perform).
#
#  Quick start:
#     make              # build ./EXHALE.x with gfortran (default)
#     make test         # run every test suite; nonzero exit if any failed
#     make check        # byte-identical regression over the whole matrix
#     make check CASES='wasp_full mol_base_handoff'   # ... over those cases
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
#  All build artifacts (objects, .mod files, dep/flag stamps) live in
#  $(OBJDIR)/; source basenames are unique, so a flat object dir is safe.
#
#  Changing FC, FFLAGS, or MODFLAG forces a full rebuild: the objects
#  depend on a stamp file whose NAME encodes a hash of the effective flag
#  string, so a different flag set names a different (missing) stamp.
# =====================================================================

OBJDIR := build
MODDIR := $(OBJDIR)
comma  := ,
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
  # The Intel compilers bring their own LAPACK: MKL, built for their own
  # runtime, so the same rule as below (a LAPACK of the compiler's own
  # prefix) is the one-line -qmkl.  Sequential: the code threads itself.
  # The compiler's runtime (libiomp5, libimf, ...) and MKL are recorded in
  # the rpath, as the OpenBLAS prefix is below, so the binary runs without
  # the oneAPI environment (the run scripts start it under env -i).
  LAPACK_LIBS ?= -qmkl=sequential -Wl,-rpath,$(TOOLCHAIN_PREFIX)/lib \
                 $(if $(MKLROOT),-Wl$(comma)-rpath$(comma)$(MKLROOT)/lib)
else
  FFLAGS  := -O3 -fopenmp
  MODFLAG := -J$(MODDIR) -I$(MODDIR)
endif

# LAPACK (banded LU dgbtrf/dgbtrs) for the steady-state Newton/PTC solver.
#
# The LAPACK must be built against the same libgfortran runtime as this
# build: the system LAPACK of this machine was built for libgfortran.so.4,
# and a binary linked against it maps BOTH .so.4 and .so.5 (measured
# 2026-09-05), which is a runtime ABI hazard. So the library is taken from
# the prefix the compiler itself lives in when that prefix carries an
# OpenBLAS (conda-forge: /opt/miniconda3, gfortran 16.2, OpenBLAS with the
# LAPACK entry points), and that prefix is recorded in the rpath so the
# binary finds its runtime without an environment. Otherwise the plain
# -llapack of the system. Override without editing this file:
#     make LAPACK_LIBS='-L/path -llapack -lblas'
FC_PATH    := $(shell command -v $(FC) 2>/dev/null)
# The first line of --version, not -dumpfullversion: it carries the package
# build number (conda-forge "gcc 16.2.0-5"), which -dumpfullversion drops, and
# two builds of one release are different compilers for the flag stamp below.
FC_VERSION := $(shell $(FC) --version 2>/dev/null | head -n 1)
TOOLCHAIN_PREFIX ?= $(abspath $(dir $(FC_PATH))..)
ifneq (,$(wildcard $(TOOLCHAIN_PREFIX)/lib/libopenblas.so))
  LAPACK_LIBS ?= -L$(TOOLCHAIN_PREFIX)/lib -lopenblas -Wl,-rpath,$(TOOLCHAIN_PREFIX)/lib -ldl
else
  LAPACK_LIBS ?= -llapack -ldl
endif
LDLIBS ?= $(LAPACK_LIBS)

# ---- source list (canonical build order, mirrors run_EXHALE.sh) -------
# The build stamp: the git revision and dirty flag the executable was built
# from, generated at Makefile PARSE time so that the file it writes is an
# ordinary source with an ordinary timestamp. It carries NO build time --
# a timestamp in a compiled constant makes the executable differ between two
# builds of the same source, which is exactly what made a stale-object
# suspicion impossible to test on 2026-09-03; the provenance header reports
# the RUN time instead, from date_and_time. The file is rewritten only when
# its content changes, so `make` reaches a fixed point and `make -q` is a
# usable up-to-dateness test.
BUILDSTAMP := $(OBJDIR)/build_stamp.f90
$(shell mkdir -p $(OBJDIR); \
  rev=$$(git rev-parse --short=12 HEAD 2>/dev/null || echo unknown); \
  if git rev-parse --git-dir >/dev/null 2>&1 && \
     [ -n "$$(git status --porcelain 2>/dev/null)" ]; then \
    dirty=dirty; else dirty=clean; fi; \
  { echo "      module build_stamp"; \
    echo "      ! GENERATED by the Makefile. Do not edit."; \
    echo "      ! The run never calls git: these are the values as of the"; \
    echo "      ! build that produced this executable. No build time is"; \
    echo "      ! stamped in -- see the Makefile comment above."; \
    echo "      implicit none"; \
    echo "      character(len=*), parameter :: build_git   = '$$rev'"; \
    echo "      character(len=*), parameter :: build_dirty = '$$dirty'"; \
    echo "      end module build_stamp"; } > $(OBJDIR)/.build_stamp.tmp; \
  cmp -s $(OBJDIR)/.build_stamp.tmp $(BUILDSTAMP) 2>/dev/null || \
    mv -f $(OBJDIR)/.build_stamp.tmp $(BUILDSTAMP); \
  rm -f $(OBJDIR)/.build_stamp.tmp)


SRC := \
  $(BUILDSTAMP) \
  src/modules/init/parameters.f90 \
  src/modules/init/species_table.f90 \
  src/modules/radiation/charge_exchange.f90 \
  src/modules/files_IO/metals_input_read.f90 \
  src/modules/files_IO/opacity_input_read.f90 \
  src/modules/files_IO/lower_atmosphere_profile.f90 \
  src/modules/files_IO/input_read.f90 \
  src/modules/files_IO/load_IC.f90 \
  src/modules/files_IO/conserved_state_restart.f90 \
  src/modules/files_IO/write_output.f90 \
  src/modules/files_IO/write_setup_report.f90 \
  src/modules/functions/grav_field.f90 \
  src/modules/functions/cross_sec.f90 \
  src/modules/functions/h2_photo_channels.f90 \
  src/modules/radiation/opacity_models.f90 \
  src/modules/functions/UW_conversions.f90 \
  src/modules/functions/utilities.f90 \
  src/modules/functions/composition.f90 \
  src/modules/functions/element_census.f90 \
  src/modules/functions/binary_element_diffusion.f90 \
  src/modules/functions/ionization_stage_transport.f90 \
  src/modules/lower_atmosphere/lower_column.f90 \
  src/modules/lower_atmosphere/h3p_cooling.f90 \
  src/modules/lower_atmosphere/molecular_infrared_data.f90 \
  src/modules/lower_atmosphere/molecular_infrared_cooling.f90 \
  src/modules/lower_atmosphere/mol_rates.f90 \
  src/modules/lower_atmosphere/h2_vibrational_relaxation.f90 \
  src/modules/lower_atmosphere/h2_self_shielding_table.f90 \
  src/modules/lower_atmosphere/molecular_reaction_heat.f90 \
  src/modules/lower_atmosphere/lyman_werner.f90 \
  src/modules/lower_atmosphere/oxygen_rates.f90 \
  src/modules/lower_atmosphere/water_photolysis.f90 \
  src/modules/lower_atmosphere/co_self_shielding_table.f90 \
  src/modules/lower_atmosphere/co_photodissociation.f90 \
  src/modules/lower_atmosphere/element_inventory.f90 \
  src/modules/lower_atmosphere/diffusive_photochemistry.f90 \
  src/modules/time_step/viscous_conduction.f90 \
  src/modules/time_step/hydrodynamic_rows.f90 \
  src/modules/time_step/conservation_budget.f90 \
  src/modules/time_step/steady_residual.f90 \
  src/modules/time_step/certification.f90 \
  src/modules/time_step/attempted_step.f90 \
  src/modules/time_step/steady_newton.f90 \
  src/modules/time_step/coupled_block_handover.f90 \
  src/modules/nonlinear_system_solver/ion_cell_state.f90 \
  src/modules/nonlinear_system_solver/ion_residual_core.f90 \
  src/modules/nonlinear_system_solver/dogleg.f90 \
  src/modules/nonlinear_system_solver/enorm.f90 \
  src/modules/nonlinear_system_solver/hybrd1.f90 \
  src/modules/nonlinear_system_solver/newton_solver.f90 \
  src/modules/nonlinear_system_solver/qform.f90 \
  src/modules/nonlinear_system_solver/r1mpyq.f90 \
  src/modules/nonlinear_system_solver/System_HeH.f90 \
  src/modules/nonlinear_system_solver/System_HeH_metals.f90 \
  src/modules/nonlinear_system_solver/System_HeH_mol.f90 \
  src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90 \
  src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90 \
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
  src/modules/radiation/electron_energy_degradation.f90 \
  src/modules/radiation/util_ion_eq.f90 \
  src/modules/radiation/ionization_equilibrium.f90 \
  src/modules/radiation/hydrogen_n2_rates.f90 \
  src/modules/radiation/lya_rt.f90 \
  src/modules/radiation/excited_hydrogen.f90 \
  src/modules/nonlinear_system_solver/T_equation.f90 \
  src/modules/post_process/post_process_adv.f90 \
  src/modules/states/caloric_eos.f90 \
  src/modules/states/base_boundary.f90 \
  src/modules/states/Apply_BC.f90 \
  src/modules/states/boundary_state_trace.f90 \
  src/modules/states/PLM_rec.f90 \
  src/modules/states/Source.f90 \
  src/modules/states/Reconstruction.f90 \
  src/modules/states/stationary_operator.f90 \
  src/modules/flux/speed_estimate_ROE.f90 \
  src/modules/flux/Num_Fluxes.f90 \
  src/modules/flux/low_mach_dissipation.f90 \
  src/modules/flux/species_face_flux.f90 \
  src/modules/time_step/RK_rhs.f90 \
  src/modules/time_step/eval_dt.f90 \
  src/modules/time_step/energy_semi_implicit.f90 \
  src/modules/init/define_grid.f90 \
  src/modules/init/set_energy_vectors.f90 \
  src/modules/init/set_gravity_grid.f90 \
  src/modules/init/set_IC.f90 \
  src/modules/init/molecular_seed_from_atomic_state.f90 \
  src/modules/init/blas_thread_policy.f90 \
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

# ---- diffusion unit tests (standalone; not built by `all`) -----------
# diffusion_tests.x exercises binary_element_diffusion on synthetic columns
# (acceptance tests T1a/T1b/T3-T7/T9/T10/T11/T12).
# It links only the module and what it uses, so it needs no LAPACK.
DIFT_SRC := \
  $(BUILDSTAMP) \
  src/modules/init/parameters.f90 \
  src/modules/init/species_table.f90 \
  src/modules/functions/utilities.f90 \
  src/modules/lower_atmosphere/lower_column.f90 \
  src/modules/lower_atmosphere/molecular_infrared_data.f90 \
  src/modules/states/caloric_eos.f90 \
  src/modules/functions/composition.f90 \
  src/modules/functions/grav_field.f90 \
  src/modules/files_IO/lower_atmosphere_profile.f90 \
  src/modules/functions/UW_conversions.f90 \
  src/modules/states/base_boundary.f90 \
  src/modules/states/Apply_BC.f90 \
  src/modules/states/PLM_rec.f90 \
  src/modules/states/Reconstruction.f90 \
  src/modules/flux/species_face_flux.f90 \
  src/modules/init/define_grid.f90 \
  src/modules/functions/binary_element_diffusion.f90 \
  src/tests/diffusion_tests.f90
DIFT_OBJ := $(addprefix $(OBJDIR)/,$(notdir $(DIFT_SRC:.f90=.o)))
DIFT_EXE := diffusion_tests.x
DIFT_DEPFILE := $(OBJDIR)/.deps_dift.mk
vpath %.f90 src/tests

# ---- element census / conservation tests (standalone) ----------------
# element_census_tests.x exercises the elemental and charge invariants of
# item A3: the stoichiometry table, the closure of one carrier write-back, the H2
# photoevent ledger, and the chemical-equilibrium base H2 branch. It links
# every module the production binary does (bar the main program), because
# the carrier write-back reads the ionization sweep's cell state, so it
# needs LAPACK like cce_probe.x.
ECT_SRC := $(filter-out src/EXHALE_main.f90,$(SRC)) \
  src/tests/element_census_tests.f90
ECT_OBJ := $(addprefix $(OBJDIR)/,$(notdir $(ECT_SRC:.f90=.o)))
ECT_EXE := element_census_tests.x
ECT_DEPFILE := $(OBJDIR)/.deps_ect.mk

# ---------------------------------------------------------------------
# cce_probe.x re-evaluates ONE saved constrained-equilibrium cell state
# without the hydrodynamics, for the derivative and conditioning diagnostics
# of the charge-exchange cancellation limit. The state is written by a
# normal run with EXHALE_CCE_DUMP set to a file name. Not part of the
# default build; build it with:  make cce_probe
# It links the same modules the equilibrium sweep uses (everything but the
# main program), so the coefficients it evaluates are the production ones,
# and it needs LAPACK for the singular values.
CCE_SRC := $(filter-out src/EXHALE_main.f90,$(SRC)) \
  src/modules/nonlinear_system_solver/constrained_equilibrium_probe.f90
CCE_OBJ := $(addprefix $(OBJDIR)/,$(notdir $(CCE_SRC:.f90=.o)))
CCE_EXE := cce_probe.x
CCE_DEPFILE := $(OBJDIR)/.deps_cce.mk

# Rebuild everything when the effective build flags change. The stamp file
# NAME encodes a hash of the full flag string -- the compiler's resolved
# path and version line, package build number included (two gfortrans of the same name on different PATHs are
# different compilers; objects of one must not be linked by the other) and
# $(FFLAGS) $(MODFLAG) --
# so a different flag set names a different stamp: the previous one becomes
# a missing prerequisite of every object and forces a rebuild, while an
# unchanged flag set names the same (already-present) stamp and rebuilds
# nothing. Deriving the name from the flags -- rather than rewriting a
# fixed-name file -- means a dry run (make -n) with other flags leaves the
# real build state untouched.
BUILDFLAGS := $(FC_PATH) $(FC_VERSION) $(FFLAGS) $(MODFLAG)
FLAGHASH   := $(firstword $(shell printf '%s' '$(BUILDFLAGS)' | cksum))
FLAGSTAMP  := $(OBJDIR)/.buildflags-$(FLAGHASH)

# Relink the executables when the LINK command changes. The compile stamp
# above does not see LDLIBS (the LAPACK the executable is linked against,
# LAPACK_LIBS), so a changed library left the old executable standing as
# "up to date" (code audit of 2026-09-29, md/CODE_AUDIT_20260929.md F7).
# The same device, on the link command: its hash names a stamp that the
# executables depend on, so a changed library relinks them and recompiles
# nothing, and an unchanged command names the present stamp.
LINKFLAGS  := $(FC_PATH) $(FC_VERSION) $(FFLAGS) $(MODFLAG) $(LDLIBS)
LINKHASH   := $(firstword $(shell printf '%s' '$(LINKFLAGS)' | cksum))
LINKSTAMP  := $(OBJDIR)/.linkflags-$(LINKHASH)

# ---------------------------------------------------------------------
.PHONY: all clean distclean ifort ifx wind_ae_ic check test diffusion_tests \
        cce_probe residual_determinism element_census_tests
all: $(EXE)
wind_ae_ic: $(WAE_EXE)
diffusion_tests: $(DIFT_EXE)
cce_probe: $(CCE_EXE)
element_census_tests: $(ECT_EXE)

# Byte-identical regression harness: rebuilds and re-runs each case of the
# matrix single-thread, comparing against the reference snapshots.
# Non-fatal when the harness is absent (e.g. a checkout without backup/).
#
# CASES limits the run to the cases named, e.g.
#     make check CASES='wasp_full mol_base_handoff'
# Empty (the default) leaves the case list to the script, which then runs its
# DEFAULT_CASES. The reference directory and the tolerance are the script's
# own environment variables (REGRESSION_GOLDEN_DIR, REGRESSION_REL_TOL), so
#     make check CASES=wasp_full REGRESSION_GOLDEN_DIR=$(PWD)/backup/regression/baseline_post170_20260905
# compares that one case against the scratch baseline instead of golden/.
CASES ?=
check: $(EXE)
	@if [ -x backup/regression/run_check.sh ]; then \
	   backup/regression/run_check.sh check $(CASES) ; \
	 else \
	   echo "regression harness not found (backup/regression/run_check.sh)"; \
	 fi

# ---- every test suite in one command --------------------------------
# `make test` runs the whole test set and prints one verdict line per suite.
#
# EVERY suite runs even after one has failed, and the target exits nonzero if
# any of them did. A test set whose job is to say which defects are still open
# -- which is exactly what it is asked to do in Phase 0 -- is useless if the
# first red suite hides the state of the rest.
#
# The two Fortran suites are built through a recursive $(MAKE) naming the
# suite as the goal, because their generated module dependency file is
# included only when the goal names it (the -include block at the end of this
# file). Their executables are <suite>.x by the naming above.
#
# The directory suites are discovered by the shell, in the recipe, not by a
# $(wildcard) at parse time: any src/tests/<name>/run.sh that is executable is
# run. A new suite is therefore added by adding its directory, with no edit
# here. The convention (Phase 0) is that run.sh builds into build/tests/<name>/,
# prints one `PASS|FAIL <name> measured=... reference=... tol=...` line per
# assertion, and exits nonzero on any failure.
test:
	@rc=0; summary=""; \
	 for suite in element_census_tests diffusion_tests; do \
	   echo ""; echo "=== $$suite ==="; \
	   if $(MAKE) --no-print-directory $$suite && ./$$suite.x; then \
	     summary="$$summary\n  PASS $$suite"; \
	   else \
	     summary="$$summary\n  FAIL $$suite"; rc=1; \
	   fi; \
	 done; \
	 echo ""; echo "=== residual_determinism ==="; \
	 if $(MAKE) --no-print-directory residual_determinism; then \
	   summary="$$summary\n  PASS residual_determinism"; \
	 else \
	   summary="$$summary\n  FAIL residual_determinism"; rc=1; \
	 fi; \
	 for sh in src/tests/*/run.sh; do \
	   [ -x "$$sh" ] || continue; \
	   name=`basename \`dirname "$$sh"\``; \
	   echo ""; echo "=== $$name ==="; \
	   if "$$sh"; then \
	     summary="$$summary\n  PASS $$name"; \
	   else \
	     summary="$$summary\n  FAIL $$name"; rc=1; \
	   fi; \
	 done; \
	 echo ""; echo "=== make test summary ==="; \
	 printf '%b\n' "$$summary"; \
	 if [ $$rc -ne 0 ]; then echo "==> TEST FAIL"; else echo "==> TEST PASS"; fi; \
	 exit $$rc

# The steady residual is a function of its argument: evaluating F at a state,
# at two others, then at the first again must return the same bits.  What used
# to carry history was the equilibrium sweep's starting composition, which the
# interface of section 147 now forbids; this test is what keeps that true.
residual_determinism: $(EXE)
	@if [ -x src/tests/residual_determinism/run_residual_determinism.sh ]; then \
	   src/tests/residual_determinism/run_residual_determinism.sh ; \
	 else \
	   echo "determinism test not found (src/tests/residual_determinism/)"; \
	 fi

$(EXE): $(OBJ) $(LINKSTAMP)
	$(FC) $(FFLAGS) $(MODFLAG) $(OBJ) -o $@ $(LDLIBS)
	@echo "built $@"

# Wind-AE IC generator (separate executable; no LAPACK needed)
$(WAE_EXE): $(WAE_OBJ)
	$(FC) $(FFLAGS) $(MODFLAG) $(WAE_OBJ) -o $@
	@echo "built $@"
$(WAE_OBJ): $(FLAGSTAMP)

# Diffusion unit tests (separate executable; no LAPACK needed)
$(DIFT_EXE): $(DIFT_OBJ)
	$(FC) $(FFLAGS) $(MODFLAG) $(DIFT_OBJ) -o $@
	@echo "built $@"
$(DIFT_OBJ): $(FLAGSTAMP)

# Constrained-equilibrium probe (separate executable; needs LAPACK for the
# singular values).
$(CCE_EXE): $(CCE_OBJ) $(LINKSTAMP)
	$(FC) $(FFLAGS) $(MODFLAG) $(CCE_OBJ) -o $@ $(LDLIBS)
	@echo "built $@"
$(CCE_OBJ): $(FLAGSTAMP)

# Element census / conservation tests (separate executable; needs LAPACK
# for the same reason cce_probe.x does).
$(ECT_EXE): $(ECT_OBJ) $(LINKSTAMP)
	$(FC) $(FFLAGS) $(MODFLAG) $(ECT_OBJ) -o $@ $(LDLIBS)
	@echo "built $@"
$(ECT_OBJ): $(FLAGSTAMP)

# compile each source to $(OBJDIR)/<base>.o (also writes its .mod there)
$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# flag-change stamp: a normal prerequisite of every object so that a change
# in FC / FFLAGS / MODFLAG (which renames the stamp) forces a full rebuild.
# The recipe drops any stale sibling stamps, so only the current flag set is
# ever present in $(OBJDIR)/.
$(OBJ): $(FLAGSTAMP)
$(FLAGSTAMP): | $(OBJDIR)
	@rm -f $(OBJDIR)/.buildflags-* && touch $@
$(LINKSTAMP): | $(OBJDIR)
	@rm -f $(OBJDIR)/.linkflags-* && touch $@

$(OBJDIR):
	@mkdir -p $@

# convenience aliases matching the run_EXHALE.sh flags
ifort: ; @$(MAKE) --no-print-directory FC=ifort
ifx:   ; @$(MAKE) --no-print-directory FC=ifx

clean:
	@rm -f $(OBJ) $(WAE_OBJ) $(DIFT_OBJ) $(CCE_OBJ) $(ECT_OBJ) $(MODDIR)/*.mod
	@echo "cleaned objects and .mod files in $(OBJDIR)/ (kept $(EXE), $(WAE_EXE))"

distclean:
	@rm -rf $(OBJDIR) $(EXE) $(WAE_EXE) $(DIFT_EXE) $(CCE_EXE) $(ECT_EXE)
	@echo "removed $(OBJDIR)/, $(EXE) and $(WAE_EXE)"

# ---- auto-generated module dependencies (skip when only cleaning) ---

$(DEPFILE): $(SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(SRC) > $@

$(WAE_DEPFILE): $(WAE_SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(WAE_SRC) > $@

$(DIFT_DEPFILE): $(DIFT_SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(DIFT_SRC) > $@

$(CCE_DEPFILE): $(CCE_SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(CCE_SRC) > $@

$(ECT_DEPFILE): $(ECT_SRC) src/utils/fortdep.py | $(OBJDIR)
	@python3 src/utils/fortdep.py --objdir $(OBJDIR) $(ECT_SRC) > $@

ifeq ($(filter clean distclean,$(MAKECMDGOALS)),)
-include $(DEPFILE)
ifneq ($(filter wind_ae_ic wind_ae_ic.x,$(MAKECMDGOALS)),)
-include $(WAE_DEPFILE)
endif
ifneq ($(filter diffusion_tests diffusion_tests.x,$(MAKECMDGOALS)),)
-include $(DIFT_DEPFILE)
endif
ifneq ($(filter cce_probe cce_probe.x,$(MAKECMDGOALS)),)
-include $(CCE_DEPFILE)
endif
ifneq ($(filter element_census_tests element_census_tests.x,$(MAKECMDGOALS)),)
-include $(ECT_DEPFILE)
endif
endif

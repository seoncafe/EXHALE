      module global_parameters
      ! Definition of global parameters and vectors
      
      implicit none
      
      integer, parameter :: outfile = 99  ! Unit number of report file
      integer, parameter :: N = 500       ! Number of computational cells
      integer, parameter :: Ng = 2        ! Number of ghost cells
      integer, parameter :: Nl_fix = 200  ! Number of default energy bins
      ! Number of ion-fraction species carried in f_sp:
      !  1-6   HI,HII,HeI,HeII,HeIII,HeITR
      !  7-9   CI,CII,CIII    10-12 OI,OII,OIII    13-15 NI,NII,NIII
      !  16-18 MgI,MgII,MgIII 19-21 SiI,SiII,SiIII  22-24 CaI,CaII,CaIII
      !  25-26 NaI,NaII       27-28 KI,KII          29-30 SI,SII
      !  31-33 FeI,FeII,FeIII
      !  34-37 H2,H2+,H3+,HeH+   (molecular extension; zero unless
      !                           thereis_mol)
      ! (Si/Ca/Fe carry three stages like Mg; Na/K/S carry two stages.)
      integer, parameter :: n_species = 37
      integer :: Nl, NlTR                 ! Number of points for energy integrations
      integer :: N_eq                     ! Numbers of equations in NL solver
      integer :: lwa                      ! Working array length for NL solver
      integer :: info                     ! Output info variable of NL solver
      integer :: j_min
      integer :: count
      
      character(len = 9), parameter   :: inp_file = 'input.inp'
      character(len = :), allocatable :: p_name
      character(len = :), allocatable :: grid_type

      ! ----- Base grid resolution ('Grid type: Mixed' only) -----
      ! The Mixed grid stacks N_low_cells uniform cells of size dr_base [R_p]
      ! on the base and fills the rest of the domain with N - N_low_cells
      ! geometrically stretched cells out to r_max. Their product
      ! dr_base*N_low_cells is the radial extent of the uniform region
      ! (0.01 R_p with the defaults below), so refining the base at fixed
      ! extent means dividing dr_base and multiplying N_low_cells by the same
      ! factor. The Uniform and Stretched grid types ignore both values.
      !
      ! Physical requirement: the uniform cells must resolve the base density
      ! scale height H = kT/(mu g). Where H/dr_base is of order a few cells,
      ! the discretization supports a stationary 2*dr entropy (contact) mode
      ! that nothing in the scheme damps -- HLLC resolves a contact of zero
      ! speed exactly, the gravity source is cell-local, and the WENO3
      ! pressure gradient only sees interface pressures. The measured
      ! dependence is a 1e-2 alternating amplitude in ln(rho) for H/dr < 5
      ! and 1e-4 for H/dr > 100 (docs/hd189_base_checkerboard.md). High-
      ! gravity planets have the least margin: on the default grid the value
      ! write_setup_report echoes is 26.9 cells per H at T_eq for
      ! HD 189733 b against 102.3 for WASP-121 b.
      ! The default is written with a DEFAULT-REAL literal (2.0e-4, not 2.0d-4)
      ! because that is what the hardcoded local in define_grid.f90 was: the
      ! value stored is the single-precision neighbour of 2e-4, 2.5e-8 relative
      ! below it. Writing 2.0d-4 here moves every base cell by that amount and
      ! changes the last few digits of a converged solution. The literal is kept
      ! as-is so that an input.inp without the key reproduces earlier runs
      ! bit-for-bit; note that spelling the default out in input.inp
      ! ("Base grid [dr,cells]: 2.0e-4 50") does NOT reproduce it, because the
      ! list-directed read into a real*8 gives the exact double 2e-4.
      real*8  :: dr_base     = 2.0e-4   ! uniform base cell size [R_p]
      integer :: N_low_cells = 50       ! number of uniform base cells

      character(len = :), allocatable :: flux
      character(len = :), allocatable :: rec_method 
      character(len = :), allocatable :: appx_mth
      character(len = :), allocatable :: sp_type
      character(len = :), allocatable :: sed_file
      ! Wind-AE warm-start IC seed (IC mode: windae). "Wind-AE seed:" sets the
      ! starting windsoln CSV (default inputdata/windae_seed.csv); the special
      ! value 'grid' auto-picks the nearest solution in inputdata/windae_grid/.
      ! "Wind-AE seed out:" optionally dumps the converged soln as a new seed.
      character(len = :), allocatable :: windae_seed_file
      character(len = :), allocatable :: windae_seed_out

      logical :: is_mom_const  = .false.  ! Is momentum constant within tolerance
      logical :: is_zero_dt    = .false.  ! Is time derivative really zero 
      logical :: force_start   = .false.  ! Force to do first 1000 iterations
      logical :: do_only_pp    = .false.  ! Do only the post processing
      logical :: do_load_IC    = .false.  ! Load existing IC
      logical :: thereis_He    = .false.  ! Is He included in computations
      logical :: do_read_sed   = .false.  ! Is numerical SED read from file
      logical :: is_monochr    = .false.  ! Is monochromatic radiation selected
      logical :: is_PL_sed     = .false.  ! Is the SED a power law
      logical :: thereis_Xray  = .false.  ! Include only EUV band	
      logical :: thereis_HeITR = .false.  ! Include calculations for He triplet
      ! He I (1^1S) photoionization source (default .false. = Verner+1996;
      !  .true. = legacy ATES two-term fit).  (No recombination switch: the
      !  Benjamin+1999 He recombination already matches modern data.)
      logical :: ates_photoion_rate = .false.
      ! H/He recombination + collisional-ionization rate model:
      !  .false. (default) = Badnell RR (+ He II DR) minus Mao & Kaastra 2016
      !   alpha_1 for case-B recombination, and Voronov 1997 collisional
      !   ionization; .true. = legacy ATES fits (Hui & Gnedin 1997 recombination,
      !   Abel+1997/HG97 collisional ionization).  Free-free always uses the
      !   van Hoof et al. 2014 Gaunt table regardless of this flag.
      logical :: legacy_hhe_rates = .false.
      ! Secondary ionization by fast photoelectrons (Shull & van Steenberg 1985).
      !  .true. (default) = high-energy photoelectrons (E0 > 40 eV) partition their
      !   excess energy into heating f_heat(x), H I secondary ionization, and He I
      !   secondary ionization following the SvS85 asymptotic fits; x is the ionized
      !   fraction of the H+He nuclei. .false. = legacy full-thermalization
      !   (bit-identical to the pre-2026 behavior).
      logical :: use_sec_ion = .true.
      ! Runtime state of the SvS85 coupling: .true. only once it is actually
      ! applied. From a cold IC the secondary-ionization base feedback amplifies
      ! the startup transient into a runaway, so the coupling is switched on only
      ! after the wind has first converged without it (set in serial code only --
      ! no threadprivate).
      logical :: sec_ion_active = .false.
      ! Input override "Secondary_ionization: Immediate": apply the coupling from
      ! step 0 (pre-staging behavior), for A/B tests only.
      logical :: sec_ion_immediate = .false.
      ! He recombination radiation ionizing H (Draine 2011 y/z parametrization,
      ! on-the-spot; see docs/QUESTIONS_2026-07-17.md).
      !  .true. (default) = the >= 24.6 eV ground-capture continuum ionizes H
      !   with the local fraction y (Draine Eq. 14.16) and the < 24.6 eV cascade
      !   photons ionize H with the density-dependent fraction z (Draine Sec.
      !   15.5); couples an extra H I photoionization rate and its photoelectron
      !   heating, and corrects the He II recombination to alpha_B + y alpha_1.
      !   In TR mode this also restores the singlet-excited capture channel
      !   (0.25 alpha_B) that the alpha_1-only network omits -- with the
      !   coupling off, the TR singlet recombination is neither case A nor
      !   case B. The photons are real; default on (2026-07-23, Update §39).
      !  .false. = He II -> He I recombination photons are all lost locally
      !   (pure case B, y=0), the legacy path.
      logical :: use_he_rec_coupling = .true.
      ! He/H diffusive separation (docs/design_hehe_diffusion.md):
      !  .false. (default) = He/H frozen at the input HeH everywhere (legacy);
      !  .true. = evolve the He element ratio with advection + molecular
      !  diffusion (He settles, He/H falls with altitude).  Metals stay frozen
      !  to H.  Off = byte-identical to legacy.
      logical :: he_diffusion = .false.
      ! Eddy (turbulent) diffusion coefficient [cm^2/s] used with he_diffusion.
      ! Mixes the He/H ratio toward uniform below the homopause (n where the
      ! molecular D equals he_kzz), preventing runaway molecular settling from
      ! the dense base; only the higher, tenuous layers separate.  Default
      ! 1e9 cm^2/s ~ Taylor et al. (2025) K_zz = 1e5 m^2/s.  Set 0 for pure
      ! molecular diffusion.  Runtime key "He_Kzz: <value>".
      real*8  :: he_kzz = 1.0d9
      ! P2b: ambipolar-corrected effective settling mass (ionized wind lifts
      !  He ions, reducing settling).  Default .true.; .false. = neutral Dm=3.
      logical :: he_ambipolar = .true.
      ! P2c: thermal-diffusion factor alpha_T for He (settling ~ (1+alpha_T)).
      !  Default 0 (off).  Runtime key "He_alphaT: <value>".
      real*8  :: he_alphaT = 0.0d0
      ! P2d: also diffuse the trace metals (each element with its own mass and
      !  binary diffusion coefficient); default .false. = metals frozen to H.
      logical :: he_metal_diffusion = .false.
      ! analytic lower column (Koskinen+2022; docs/lower_atmosphere_
      ! coupling.*): radius of the 1-bar level [R_J].  If > 0, on startup the
      ! isothermal-Teq hypsometric column with chemical-equilibrium H2/H/He
      ! is integrated from 1 bar to 1 ubar and the derived base radius
      ! (equilibrium + fully-atomic bracket), base H2/H/He fractions and mu
      ! are REPORTED next to the input "Planet radius" (consistency check;
      ! nothing is overridden).  <= 0 (default) = off.
      real*8  :: lower_col_r1bar = -1.0d0
      ! Lower-atmosphere pre-step ("Lower atmosphere: vulcan|analytic <R_1bar>"):
      !  0 = off (default; classic base, or a hand-made base.inp),
      !  1 = analytic chemical-equilibrium column (src/utils/run_lower.py),
      !  2 = VULCAN photochemistry (src/utils/vulcan_driver.py; the
      !      first run takes HOURS, later runs reuse the cached .vul).
      ! When set and no base.inp exists, EXHALE invokes the generator itself
      ! (EXECUTE_COMMAND_LINE) and then reads the produced base.inp -- i.e.
      ! VULCAN acts as a subroutine-style pre-step of the wind solve.
      integer :: lower_atm_mode  = 0
      real*8  :: lower_atm_r1bar = -1.0d0
      ! passive molecular base: reduce the base particle count
      ! (ntot_bc) by the H nuclei bound into H2 at (p_base_bar, T0), taken
      ! from the photochemical handoff when one was supplied and from the
      ! chemical-equilibrium fit otherwise.  EOS-only correction (lower base
      ! pressure, heavier base mu); the chemistry stays atomic -- crude,
      ! documented in docs/lower_atmosphere_coupling.*.
      ! Key "Molecular base: True".
      logical :: molecular_base = .false.
      ! H2 volume mixing ratio q_H2 = n_H2/(n_H2+n_H+n_He) at the base, as
      ! determined by the lower-atmosphere photochemistry and carried by the
      ! "q_H2_base" key of base.inp.  This is the same quantity the
      ! chemical-equilibrium fit q_h2_equilibrium returns, so no conversion
      ! is involved: when supplied it REPLACES the fit in the molecular-base
      ! particle count (docs/base_composition_handoff_plan.md).  Negative
      ! (default) = no photochemical value, the fit is used, so a base.inp
      ! without the key behaves exactly as before.
      real*8  :: q_h2_base  = -1.0d0
      ! Pressure level of the lower-atmosphere handoff [bar] ("p_base" in
      ! base.inp; --pbase of src/utils/vulcan_to_base.py).  The
      ! chemical-equilibrium fit is evaluated at this level, so the fit and
      ! the handoff always describe the same level.  Default 1 microbar.
      real*8  :: p_base_bar = 1.0d-6
      ! molecular chemistry (H2/H2+/H3+/HeH+ in the coupled ionization
      ! equilibrium + H2 photoionization opacity/heating + H3+ IR cooling).
      ! Key "Molecular chemistry: True".  v1 constraints: requires He; not
      ! combined with trace metals (input_read errors out).  Default off =
      ! byte-identical legacy.  docs/lower_atmosphere_coupling.*.
      logical :: thereis_mol = .false.
      logical :: thereis_metals = .false. ! Include trace-metal species
                                          !  (C/N/O/Mg/Si/Ca/Na/K/S/Fe)
      ! EOS mass/electron/particle policy. .true. (default) = metals enter
      ! the bulk gas budget consistently: their mass in rho_bc/calc_rho
      ! (hence mu, v0, p0, b0), their electrons in calc_ne/dp_bc, and
      ! their nuclei in calc_ntot/the ghost pressure (ntot_bc) -- matching
      ! the charge balance already used inside the MINPACK ionization
      ! systems. .false. = legacy trace approximation (H/He-only budget;
      ! metals enter ionization/cooling but not the bulk gas). Runtime
      ! key: 'eos_metals 0|1' in metals.inp. Metals-off runs are
      ! identical either way (all metal sums vanish).
      logical :: eos_include_metals = .true.
      ! Composition constants derived in input_read (legacy values when
      ! eos_include_metals is off or metals are absent):
      real*8 :: mass_per_H = 1.0d0  ! gas mass per H nucleus [m_H]:
                                    !  1 + 4*HeH (+ sum melem_ab*melem_A)
      real*8 :: ntot_bc    = 1.0d0  ! total nuclei density at base [n0]:
                                    !  (1 + HeH + sum melem_ab)/(1 + HeH)
      logical :: thereis_lowIP_metal = .false. ! An active metal whose neutral
                                          !  ionization potential lies below the
                                          !  13.6 eV HI edge (e.g. Mg I, 7.646
                                          !  eV) is present; triggers extension
                                          !  of the energy grid/SED below 13.6 eV
      logical :: use_2lev_cool  = .false. ! Two-level fine-structure metal
                                          !  cooling ([O I] 63um, [C II] 158um)
      logical :: cno_chianti    = .true.  ! C/N/O line cooling source:
                                          !  .true. (default) = CHIANTI v11
                                          !   closed-form fits incl. N I/N II,
                                          !  .false. = legacy AIOLOS analytic
                                          !   fits (C I/C II/O I/O II only;
                                          !   no N cooling; O off by 40-70%
                                          !   vs CHIANTI in the wind region)
                                          !   (metals.inp key 'cno_cool 0|1')
      integer :: pp_metal_mode  = 1       ! Metal treatment in the advection
                                          !  post-process (post_process_adv):
                                          !  0 = metal-free (legacy: metals
                                          !      zeroed in _adv heating/cooling);
                                          !  1 = frozen (use converged eq metal
                                          !      densities, held fixed) (default);
                                          !  2 = re-solve metal ionization at the
                                          !      post-process T/n_e (not yet
                                          !      implemented; warns and falls
                                          !      back to mode 1).
                                          !  Default is 1 (frozen): the frozen
                                          !  _adv profile already tracks the eq
                                          !  solution closely, so the re-solve
                                          !  would change little.
                                          !  Set at runtime via 'pp_metals' in
                                          !  metals.inp. No effect when metals
                                          !  are off (nm = 0 either way).
      logical :: use_weno3     = .false.  ! Use WENO3 reconstruction
      logical :: use_plm       = .false.  ! Use PLM reconstruction
      logical :: recon_two_stage = .false.  ! "Reconstruction scheme: PLM+WENO3": run PLM (stage 1) then WENO3 (stage 2), using BOTH du_th values. PLM/WENO3 alone are single-stage and use only the first du_th value.
      logical :: is_stalled    = .false.  ! Convergence stalled at a du plateau
      logical :: spherical_domain = .false. ! Domain mode (set via input.inp):
                                          !  .false. = full Roche potential
                                          !   truncated at the Hill/L1 radius
                                          !   (default ATES behavior);
                                          !  .true.  = pure planetary potential
                                          !   (-b0/r), tidal + centrifugal terms
                                          !   dropped, domain extended to
                                          !   r_out_user [R_p] (Huang Case A-like)
      integer :: ic_mode       = 0        ! IC selection mode (input.inp
                                          !  "IC mode: <word>"): 0 = cold
                                          !  hydrostatic (default), 1 =
                                          !  transonic, 2 = hot_parker,
                                          !  3 = auto (select_IC_auto picks
                                          !  the family from the cold
                                          !  sonic-point topology; see
                                          !  docs/auto_ic_design.md);
                                          !  4 = windae (init.f90 builds an
                                          !  in-process Wind-AE warm-start
                                          !  IC via wae_exhale_bridge, then
                                          !  load_IC reads it). The explicit
                                          !  legacy keys below take
                                          !  precedence over auto.
      logical :: transonic_ic  = .false.  ! IC type (set via input.inp,
                                          !  "Transonic IC: True"; default off):
                                          !  .false. = isothermal hydrostatic
                                          !   atmosphere + v~0 seed (default);
                                          !  .true.  = transonic isothermal-wind
                                          !   IC (Parker-type, solved from the
                                          !   ATES potential via the Bernoulli
                                          !   integral). Needed for deep-RLOF
                                          !   bases (sonic point near L1) that
                                          !   launch no wind from the hydrostatic
                                          !   IC. Auto-enables force_start so the
                                          !   cold wind heats to its hot steady
                                          !   state before the convergence test
                                          !   (the IC already has rho*v*r^2=const,
                                          !   so du~0 would exit prematurely).
                                          !  Falls back to the hydrostatic IC if
                                          !   no interior sonic point exists.

      ! Hot-Parker IC (set via input.inp "Hot Parker IC: <Twind[K]>"; default
      ! off): a transonic isothermal-wind IC built at a WARM effective wind
      ! temperature T_wind_ic (instead of the cold base T0), with T and H/He
      ! ionization ramped from the cold neutral base to the hot ionized wind, so
      ! v, rho, T, and ionization all start close to the final escaping-wind state
      ! (a p-winds-style Parker initial guess). Auto-enables force_start.
      logical :: hot_parker_ic = .false.
      real*8  :: T_wind_ic     = 1.0d4    ! effective wind temperature [K] for the IC
      ! Cold-base transition radius [R_p] for the hot-Parker IC: below it the IC
      ! is forced to the cold, static, cold-pressure base state (matching the
      ! lower BC), above it the warm Parker wind is ramped in. This removes the
      ! inverted base pressure gradient that otherwise drives a base inflow and
      ! a thermal runaway (see docs/initial_condition_benchmark, root-cause).
      real*8  :: hp_base_rtr   = 1.3d0    ! transition radius [R_p]; >1

      ! Nonlinear-solver selection (set via input.inp; defaults = the upgraded
      ! solvers). use_newton_ieq: ionization equilibrium via the analytic-
      ! Jacobian Newton (+hybrd1 fallback) vs. legacy MINPACK hybrd1
      ! ("Newton solver: False"). use_brent_tsolve: the post-process energy
      ! equation via the bracketing Brent root-finder vs. legacy hybrd1 + the
      ! 2x-band reject ("Brent solver: False"). See docs/Update_EXHALE_solver.
      logical :: use_newton_ieq  = .true.
      logical :: use_brent_tsolve = .true.

      !------- Global constants -------!
      
      ! Physical constants
      real*8,parameter ::  pi      = 3.1415926536     ! pi   
      real*8,parameter ::  kb_erg  = 1.38e-16         ! Boltzmann constant in CGS units
      real*8,parameter ::  kb_eV   = 8.6167e-05       ! Boltzmann constant (eV/K)
      real*8,parameter ::  mu      = 1.673e-24        ! Hydrogen mass (g)
      real*8,parameter ::  g       = 1.666666666667   ! Polytropic index
      real*8,parameter ::  Gc      = 6.67259e-8       ! Gravitational constant (CGS)   
      real*8,parameter ::  erg2eV  = 6.241509075e11   ! 1 erg measured in eV
      real*8,parameter ::  hp_erg  = 6.62620e-27      ! Planck constant in CGS units
      real*8,parameter ::  hp_eV   = 4.1357e-15       ! Planck constant (eV*s)
      real*8,parameter ::  c_light = 2.99792458e10    ! Speed of light in cm/s
      real*8,parameter ::  parsec  = 3.08567758147e18 ! 1 pc in cm
      real*8,parameter ::  AU      = 1.495978707e13   ! Astronomical unit
      real*8,parameter ::  RJ      = 6.9911e9         ! Jupiter radius (cm)
      real*8,parameter ::  MJ      = 1.898e30         ! Jupiter mass (g)
      real*8,parameter ::  Msun    = 1.989e33         ! Sun mass (g)
      real*8,parameter ::  Rsun    = 6.957e10         ! Sun radius (cm)
      real*8,parameter ::  R_earth = 6.3725e8         ! Earth radius (cm)
      real*8,parameter ::  M_earth = 5.9726e27        ! Earth mass (g)
      real*8,parameter ::  ih      = 1.0              ! Atomic number of Hydrogen
      real*8,parameter ::  ihe     = 2.0              ! Atomic number of Helium
	
	   ! -- - Energy constants

      ! Energy intervals
      real*8 ::  e_top
      real*8 ::  e_mid
      real*8 ::  e_low
      
      ! Threshold energies
      real*8,parameter ::  e_th_HI   = 13.6      ! Threshold for HI ionization
      real*8, parameter ::  e_th_H2  = 15.4d0  ! H2 photoionization threshold [eV]
      real*8,parameter ::  e_th_HeI  = 24.6      ! Threshold for HeI ionization
      real*8,parameter ::  e_th_HeII = 54.4      ! Threshold for HeII ionization
      real*8,parameter ::  e_th_HeTR = 4.80      ! Threshold for HeI triplet ionization
      ! Photoelectron energy threshold above which the SvS85 secondary-ionization
      ! partition is applied (40 eV, matching the wind_ae X-ray cutoff). Below it a
      ! photoelectron thermalizes fully. Note SvS85 is strictly an E0 >~ 100 eV
      ! asymptotic fit; using it down to 40 eV is a deliberate approximation.
      real*8, parameter :: E_sec_ion = 40.0d0
      real*8,parameter ::  e_th_MgI  = 7.646     ! Threshold for MgI ionization
      real*8,parameter ::  e_th_MgII = 15.035    ! Threshold for MgII ionization

	   ! Numerical constants
      real*8 ::  CFL    = 0.6         ! CFL number; settable in input.inp via "CFL:" (lower = smaller dt, may damp a numerical limit cycle)
      real*8 ::  du_th     = 1.0d-3   ! final (stage-2 / WENO3) escape-momentum threshold; settable in input.inp via "du_th [PLM,WENO3]:". 1e-3 is the original ATES-Code-main value (EXHALE had loosened it to 2e-2, accepting ~2% mass-flux spread).
      real*8 ::  du_th_plm = -1.0d0   ! stage-1 (PLM) threshold; if > du_th the run is two-stage: PLM until du<du_th_plm, then switch reconstruction to WENO3 and converge at du<du_th. <=0 => single-stage at du_th.
      real*8,parameter ::  dtu_th = 1.0d-8      ! Threshold variation of time deriv.
      real*8           ::  du                   ! Initial momentum variation
      real*8           ::  dtu                  ! Norm of time derivative

      ! Convergence stall detection (ported from ATES_extended): stop runs
      ! whose du settles on a plateau above du_th (steady state reached but
      ! the strict thresholds physically unreachable for the given setup).
      integer,parameter :: count_max = 1000000  ! hard cap on iterations
      integer :: N_stall   = 2000     ! window of stalled steps; settable via
                                      !   input line "Stall [tol,N]: <tol> <N>"
      real*8  :: stall_tol = 1.0d-6   ! rel. du change defining a stall (same line)

      ! Mass-flux LEVEL stability tolerance: in stage 2 a run only counts as
      ! converged/stalled when the mean |rho*v*r^2| over [j_min:N] changed by
      ! less than lev_th (relative) across the last N_stall steps. du is the
      ! SPREAD of the mass flux and is blind to a uniform drift of the whole
      ! profile; without this gate, runs can stop on spatially flat but still
      ! level-drifting states that differ by tens of percent in Mdot
      ! (path-dependent quasi-steady snapshots; see refactor_stage1_log.md).
      ! "Level tol: <val>" in input.inp; <= 0 disables (DEFAULT, legacy
      ! behavior preserved byte-identically). Kept opt-in until a steady-
      ! state Newton/PTC reference calibrates the production tolerance: the
      ! WASP-121b study showed the true converged state needs >35k marching
      ! steps (du dips below 1e-3 long before, at a premature ~0.08-dex-low
      ! snapshot), so a default-on level gate would change every result
      ! without a validated "truth" to anchor the tolerance to.
      real*8  :: lev_th = -1.0d0

      ! Production steady-solver wiring ("Solver: Newton [du_switch]").
      ! When .true., the normal marching loop (including the automatic
      ! two-stage PLM->WENO3) runs as a WARM-UP; the JFNK steady solver then
      ! finishes the run to ||R|| < resid_th (default 1e-3 if "Resid tol:" was
      ! not given) and the standard final outputs / post-processing follow.
      ! Requires a smooth base valve: if "Valve eps:" was not set, 1e-4 is
      ! adopted.
      logical :: use_newton_solver = .false.
      ! Hand-off to JFNK is keyed on the FLUX metric du (the radial spread of
      ! rho*v*r^2), consistent with the flux-based convergence decision: once the
      ! cheap marching has flattened the wind to du < newton_du_switch the Newton
      ! solver polishes it the rest of the way (Newton's quadratic convergence
      ! tightens du 1e-2 -> 1e-3 far faster than continued marching). Settable as
      ! "Solver: Newton [du_switch]". (Earlier revisions keyed the hand-off on the
      ! residual itself, ||R|| < 5e-2; the du-based key is equivalent in practice
      ! and matches the flux-based convergence decision.)
      real*8  :: newton_du_switch  = 1.0d-2

      ! WENO3 weight-freezing mode for the steady Newton solver (lagged /
      ! frozen nonlinear weights, the standard FV steady-solve remedy):
      !   0 = compute weights fresh every call (DEFAULT, byte-identical)
      !   1 = compute AND store the smoothness factors S0/S1
      !   2 = reuse the stored factors (linearizes the reconstruction so the
      !       inner Newton problem is mildly nonlinear; refreshed per outer
      !       iteration by a mode-1 evaluation)
      integer :: weno_mode = 0

      ! Base-valve smoothing ("Valve eps: <v_eps>", code velocity units).
      ! The lower-BC ghost velocity is v_g = max(v_1, 0) (one-way valve); the
      ! kink at v_1 = 0 sits exactly where the breathing base lives (steady
      ! states hug v_1 -> 0+), making the steady residual non-differentiable
      ! there and blocking Newton line searches. With valve_eps > 0 the valve
      ! becomes the softplus 0.5*(v + sqrt(v^2 + eps^2)) - 0.5*eps, which is
      ! smooth, -> max(v,0) as eps -> 0, and = v - O(eps^2/v) for v >> eps.
      ! Default <= 0 keeps the exact legacy valve (byte-identical).
      real*8  :: valve_eps = -1.0d0

      ! Momentum-consistent base pressure ("Hydrostatic base: True"). The
      ! legacy lower BC pins the ghost pressure to ntot_bc+dp_bc (T=T0), so
      ! dp/dr is flattened to ~0 at the base -- but steady momentum balance
      ! needs dp/dr = -rho g. For strongly bound planets (large b0, e.g.
      ! HD189733b) the resulting base-face momentum residual is O(10-100) and
      ! drives the breathing limit cycle. With this flag the ghost pressure is
      ! instead a linear extrapolation of the interior pressure gradient
      ! (cells 1,2), so dp/dr stays CONTINUOUS through the base and can balance
      ! gravity; rho stays anchored at rho_bc and T_base floats slightly off
      ! T0. Default off = byte-identical legacy behavior.
      logical :: hydrostatic_base = .false.

      ! Base boundary-condition anchor ("Base BC: density" / "pressure [<p_ubar>]").
      !   density  (0, default, legacy): the base number density is fixed at n0
      !            ( = 10^"Log10 lower boundary number density" ), so the ghost
      !            pins rho = rho_bc and p = ntot_bc + dp_bc (=> T = T0).
      !   pressure (1, CETIMB-style): the base is anchored by PRESSURE instead.
      !            n0 is derived so the base pressure n0*kb*T0*ntot_bc equals
      !            base_p_ubar [microbar] (1 microbar = 1 erg/cm^3). The same
      !            isothermal ghost then pins p (= base_p_ubar) + T0 and derives
      !            rho. Because the base is over-determined isothermal, this is a
      !            re-parameterization of the SAME ghost: its physical effect is a
      !            much less dense base (1 microbar ~ 20x below n0=1e14), shrinking
      !            the dense-base gravity source rho*g that drives the breathing.
      integer :: base_bc_mode = 0
      real*8  :: base_p_ubar   = 1.0d0   ! target base pressure [microbar] (mode 1)

      ! Periodic Shapiro (1970) low-pass filter ("Shapiro filter: <eps> [<every>]").
      ! Damps the gravity-unbalanced sound-wave (base-breathing) instability the
      ! same way CETIMB (Koskinen et al. 2013a) does: every shapiro_every steps a
      ! weak 1-2-1 filter u_j += (eps/4)(u_{j-1}-2u_j+u_{j+1}) is applied to the
      ! conservative variables. shapiro_eps in [0,1]; eps <= 0 disables it.
      ! OFF BY DEFAULT (opt-in). It damps the HD189733b breathing TRANSIENT, but
      ! the HD209458b cold-IC S x V sweep (docs/base_breathing_progress.md,
      ! HD209458b_test/) showed that with the filter ON the clean cold-IC solution
      ! is driven OFF the transonic-wind saddle into the INFALL attractor (v < 0
      ! everywhere) regardless of the velocity BC, while filter-OFF relaxes to a
      ! clean outflow. So it must NOT be a global default; enable for each run with
      ! "Shapiro filter: <eps> [<every>]" only for cases that actually breathe.
      real*8  :: shapiro_eps   = -1.0d0
      integer :: shapiro_every = 4

      ! CETIMB-style base velocity ("Base velocity: massflux" / "valve"). The
      ! legacy lower BC valves v (max(v1,0)); CETIMB (Koskinen 2013a) instead sets
      ! the base velocity from the steady mass-flux continuity rho0*v0*r0^2 = F_c,
      ! with F_c the wind's flux constant. Here F_c = mean(rho*v*r^2) over
      ! [j_min:N] (the escape / constant-momentum region, NOT the base where
      ! rho*v*r^2 is not yet flat), updated each step, and v_ghost = F_c/(rho_bc
      ! r^2). OFF BY DEFAULT (opt-in). It is physically benign (never causes
      ! infall on its own) and gives the cleanest base (suppresses the +-m/s base
      ! oscillation), but it slows convergence ~5x (554k vs 107k steps on the
      ! HD209458b cold IC), so the fast legacy valve stays the default. Enable
      ! with "Base velocity: massflux".
      logical :: base_v_massflux = .false.
      real*8  :: base_flux_const = -1.0d0   ! F_c [code units], updated each step

      ! Molecular transport (Navier-Stokes level), the physical damping of
      ! the near-base momentum imbalance that CETIMB carries and EXHALE's
      ! inviscid HLLC scheme lacks. Full derivation, discretization and
      ! boundary treatment: src/modules/time_step/viscous_conduction.f90 and
      ! docs/viscosity_conduction.md. Both switches default OFF, so a run
      ! without the keys is byte-identical to the inviscid code.
      !
      ! "Viscosity: True" -- radial viscous force (div tau)_r plus its
      !  dissipation q_mu = (4/3) mu (dv/dr - v/r)^2 in the energy equation,
      !  with mu(T) = (4/15)(m_H/k_B) kappa(T) (monatomic Chapman-Enskog,
      !  Prandtl 2/3) tied to the Watson et al. (1981) atomic-hydrogen
      !  conductivity.
      logical :: visc_on = .false.
      ! "Conduction: True" -- heat conduction (1/r^2) d/dr(r^2 kappa dT/dr)
      !  with kappa(T) = 4.45e4 (T/1000 K)^0.7 erg cm^-1 s^-1 K^-1.
      logical :: cond_on = .false.
      ! "Viscosity: <mu0> [<s>]" -- diagnostic power-law override in CODE
      !  units, mu = visc_mu0*T^visc_s (visc_mu0 > 0 takes precedence over
      !  the calibrated form). Kept for sensitivity scans.
      real*8 :: visc_mu0 = 0.0d0
      real*8 :: visc_s   = 0.7d0

      ! Residual-based convergence ("Resid tol: <val>"). When > 0, the run
      ! converges when the finite-volume steady residual
      !   R = dF - S        (mass, momentum)
      !   R = dF_E - S_E - (heat - cool)   (energy)
      ! has max_j |R(j,k)|/max|u(:,k)| over [j_min:N] below resid_th, REPLACING
      ! the du<du_th test (du measures only the mass-flux spread and is blind
      ! to an operator-split energy imbalance: at the premature WASP golden,
      ! du/dtu were tiny while the energy residual was ~30). Evaluated every
      ! N_resid steps (cheap: one Reconstruct + RK_rhs, no extra ioniz_eq).
      ! Default <= 0 disables (legacy du-based stop, byte-identical).
      real*8  :: resid_th = -1.0d0
      integer :: N_resid  = 500
      ! Residual NORM used for the convergence / Newton-trigger test.
      ! .true. (DEFAULT) = volume-weighted  sum_j|R(j,k)|V_j / sum_j|u(j,k)|V_j
      !   (V_j = r_j^2 dr_j): the fractional drift rate of the volume-integrated
      !   conserved quantity. This is the physically meaningful steady measure
      !   and, unlike the L-inf max, is NOT inflated by isolated small near-base
      !   cells on the non-uniform grid (where R ~ 1/dr).
      ! .false. = legacy L-inf  max_j|R(j,k)| / max_j|u(:,k)|.
      ! Set "Resid norm: Linf" in input.inp to revert.
      logical :: resid_vol = .true.

      ! Energy source-term integrator: .true. = semi-implicit backward-Euler
      ! cell solve (EXHALE default); .false. = original explicit forward
      ! Euler ("Energy solver: Explicit" in input.inp). Runtime switch kept for
      ! solver-component isolation tests.
      logical :: use_semi_implicit_energy = .true.

      ! Local (cell-by-cell) pseudo-time stepping ("Time stepping: Local" in
      ! input.inp). Steady-state acceleration only: each cell marches with its
      ! own CFL step dt_j = CFL*dr_j/(|v|+cs), removing the global-dt
      ! bottleneck set by the smallest base cell. Time accuracy is lost but
      ! the steady fixed point (dF = S, heating = cooling) is unchanged.
      ! Default .false. = original global-dt marching, byte-identical.
      logical :: use_local_dt = .false.
      
      !------- Planetary parameters -------!

      real*8  ::  n0          ! Density at lower boundary (cm^-3)
      real*8  ::  R0          ! Planetary radius (cm)
      real*8  ::  Mp          ! Planet mass (g)
      real*8  ::  T0          ! Temperature at lower boundary (K)
      real*8  ::  LX          ! Log10 X-Ray luminosity (erg/s)
      real*8  ::  LEUV        ! Log10 EUV luminosity (erg/s)
      real*8  ::  J_XUV       ! Bolometric star XUV flux (erg/cm^2*s)
      real*8  ::  Mstar       ! Mass of companion star
      real*8  ::  Mrapp       ! Ratio M_star/M_p
      real*8  ::  atilde      ! Orbital radius in unit of R0 (a/R0)
      real*8  ::  r_esc       ! Escape radius for constant momentum
      real*8  ::  a_orb       ! Orbital distance
      real*8  ::  r_max       ! Maximum radius = Roche Lobe dimension
      real*8  ::  r_out_user  ! User-set outer radius [R_p], spherical mode only
      real*8  ::  HeH         ! He/H ratio
      real*8  ::  X_C         ! C/H number ratio (solar ~ 2.69e-4)
      real*8  ::  X_N         ! N/H number ratio (solar ~ 6.76e-5)
      real*8  ::  X_O         ! O/H number ratio (solar ~ 4.90e-4)
      real*8  ::  X_Mg        ! Mg/H number ratio (solar ~ 3.98e-5)
      real*8  ::  X_Si        ! Si/H number ratio (solar ~ 3.24e-5)
      real*8  ::  X_Ca        ! Ca/H number ratio (solar ~ 2.19e-6)
      real*8  ::  X_Na        ! Na/H number ratio (solar ~ 1.74e-6)
      real*8  ::  X_K         ! K/H  number ratio (solar ~ 1.07e-7)
      real*8  ::  X_S         ! S/H  number ratio (solar ~ 1.32e-5)
      real*8  ::  X_Fe        ! Fe/H number ratio (solar ~ 3.16e-5)
      real*8, allocatable :: melem_ab(:) ! abundance for each element, canonical
                                          !  element order (iel_*); set in
                                          !  input_read so metal code can index
                                          !  abundance by element, not by name
      real*8  ::  rho_bc      ! Adimensional number density at origin
      real*8  ::  a_tau       ! Rate correction coefficient
      real*8  ::  PLind       ! Index of spectral power law
      real*8  ::  dp_bc       ! Boundary condition for pressure
	
      !------- Normalizations -------!

      real*8 :: v0      ! Velocity normalization
      real*8 :: t_s     ! Time normalization
      real*8 :: p0      ! Pressure normalization
      real*8 :: q0      ! Scale normalization 
      real*8 :: b0      ! Jeans parameter at planet surface

      !------- Global vectors -------!
      
      real*8, dimension(:), allocatable :: e_v, de_v
      real*8, dimension(:), allocatable :: s_hi,s_hei,s_heii,s_heiTR
      ! H2 photoionization cross section on the energy grid (molecular;
      ! Yan+1998 fit, filled in set_energy_vectors)
      real*8, dimension(:), allocatable :: s_h2
      ! Metal photoionization cross sections, one column per photo-ionizable
      ! metal ion in species_table iphot order (1=CI,2=CII,3=OI,4=OII,
      ! 5=NI,6=NII,7=MgI,8=MgII,9=SiI,10=SiII,11=CaI,12=CaII,13=NaI,14=KI,
      ! 15=SI). Replaces the former loose s_ci..s_mgii.
      real*8, dimension(:,:), allocatable :: sigma_tab

      !------- Opacity-model dispatcher (set via opacity.inp) -------!
      ! See src/modules/radiation/opacity_models.f90 for the dispatcher.
      character(len = 1) :: opacity_model = 'A'  ! 'A','C','P','T'
      ! Constant-model multiplicative factors (default 1.0):
      real*8 :: opa_const_HI    = 1.0d0
      real*8 :: opa_const_HeI   = 1.0d0
      real*8 :: opa_const_HeII  = 1.0d0
      real*8 :: opa_const_HeITR = 1.0d0
      ! Robinson & Catling pressure broadening f(p)=1+a*(p/p_pivot)^n:
      real*8 :: opa_pb_factor   = 0.0d0    ! a
      real*8 :: opa_pb_exponent = 1.0d0    ! n
      real*8 :: opa_pb_pivot    = 1.0d5    ! dyne/cm^2 = 0.1 bar
      ! Tabulated-model (.opa) file paths for each species:
      character(len = :), allocatable :: opa_file_HI
      character(len = :), allocatable :: opa_file_HeI
      character(len = :), allocatable :: opa_file_HeII
      character(len = :), allocatable :: opa_file_HeITR
      real*8, dimension(:), allocatable :: F_XUV
      ! Cell-by-cell pressure-broadening multiplier for the opacity ('P'
      ! model). 1.0 everywhere unless opacity_model='P'; set before each
      ! photoheating call and used to weight the opacity column density.
      real*8, dimension(1-Ng:N+Ng) :: opa_pf = 1.0d0
      real*8, dimension(1-Ng:N+Ng) :: r,r_edg,dr_j
      real*8, dimension(1-Ng:N+Ng) :: Gphi_c,Gphi_i

      !------- excited hydrogen H(n=2) coupling -------!
      ! Christie+2013 / Huang+2017 n=2 (2s/2p) model feeding back into
      ! the H ionization balance and the energy equation. All terms are
      ! ZERO unless use_excited_H = .true. (enabled at runtime when a
      ! stellar T_eff is supplied in input.inp), so the default build is
      ! byte-identical to the result with excited-H off. See excited_hydrogen.f90.
      logical :: use_excited_H   = .false. ! master switch (set if T_star_eff>0)
      logical :: incl_deexc_heat = .false. ! add collisional de-excitation
                                           !  heating (overlaps the existing HI
                                           !  coex cooling; off by default)
      real*8  :: T_star_eff = 0.0d0        ! stellar effective temperature [K]
                                           !  (<=0 disables the Balmer continuum)
      real*8  :: R_star     = 0.0d0        ! stellar radius [cm] (Balmer dilution)
      real*8  :: gamma2_bal = 0.0d0        ! n=2 Balmer photoion. rate [s^-1]
      real*8  :: hpe2_bal   = 0.0d0        ! photoelec. heat per n=2 atom [erg/s]
      ! Ly-alpha mean-intensity source for the n=2 pumping term:
      !  jlya_mode = 0 -> parameterized J_Lya = 0.1*F_LyC/Dnu_D with a
      !                    top-down 1/(1+tau_Lya) line-center attenuation;
      !  jlya_mode = 1 -> read a J_Lya(r) profile from jlya_rt_file (an
      !                    external RT result, e.g. a real Monte Carlo) directly;
      !  jlya_mode = 2 -> in-line escape-probability RT (Neufeld/Harrington wing
      !                    escape), computed every timestep from the current state
      !                    in excited_H_update (see lya_rt.f90).
      integer :: jlya_mode  = 0
      character(len=200) :: jlya_rt_file = 'jlya_rt.txt'
      ! ----- Ly-alpha radiative transfer ----- !
      !  F_Lya_star -> incident stellar Ly-alpha flux at the planet
      !  [erg cm^-2 s^-1] (the stellar source of J_lya in jlya_mode 2; Case D
      !  scales it x0.35). The "cooling trapping" is realized through the
      !  H(n=2) heating channels (photoelectric + collisional de-excitation), as in
      !  Huang et al. 2023, NOT a separate escape-probability factor on the Ly-alpha
      !  cooling. OFF by default => byte-identical to the no-Lya result.
      real*8  :: F_Lya_star = 0.0d0
      ! Stellar Ly-alpha line half-width [km/s] (broad plateau ~+-70 km/s,
      ! Huang+2017); sets how deep the (trapped) stellar beam penetrates via the
      ! line wings. Editable in input.inp ("Lya stellar halfwidth [km/s]:").
      real*8  :: dv_star_lya = 70.0d0
      ! Bounded stellar trapping buildup: the penetrating stellar Ly-alpha photons
      ! scatter and build the mean intensity above free-streaming by E=min(boost,
      ! 1/beta) (-> 1 in the thin outer wind, capped at lya_star_boost in the thick
      ! region; the full 1/beta over-counts). Tuned to Huang+2023 Fig. 11.
      ! Editable: "Lya stellar boost [-]:".
      real*8  :: lya_star_boost = 5.0d0
      ! Cell-by-cell feedback arrays injected into ioniz_eq (zero unless enabled):
      real*8, dimension(1-Ng:N+Ng) :: gph_balmer_HI = 0.0d0 ! extra HI photoion [s^-1]
      real*8, dimension(1-Ng:N+Ng) :: heat_balmer   = 0.0d0 ! extra heat [erg cm^-3 s^-1]
      ! Cell-by-cell diagnostics for output/Excited_H.txt:
      real*8, dimension(1-Ng:N+Ng) :: Jlya_arr  = 0.0d0 ! Lya mean intensity [cgs]
      real*8, dimension(1-Ng:N+Ng) :: n2s_arr   = 0.0d0 ! H(2s) density [cm^-3]
      real*8, dimension(1-Ng:N+Ng) :: n2p_arr   = 0.0d0 ! H(2p) density [cm^-3]
      real*8, dimension(1-Ng:N+Ng) :: Sproton_arr = 0.0d0 ! Balmer proton src [cm^-3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: Hpe_arr   = 0.0d0 ! photoelec. heat [erg cm^-3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: Hdx_arr   = 0.0d0 ! de-excit. heat [erg cm^-3 s^-1]
      ! Ground-state H proton-budget rate coefficients, captured by ioniz_eq on
      ! the converged pass so write_excited_H can compare the n=2 photoionization
      ! proton source against the ground-state channels (Huang Figs. 11/27):
      real*8, dimension(1-Ng:N+Ng) :: gph_ground_HI = 0.0d0 ! ground-state HI photoion. [s^-1]
      real*8, dimension(1-Ng:N+Ng) :: cion_HI       = 0.0d0 ! HI collisional ioniz. [cm^3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: arec_HII      = 0.0d0 ! HII recombination [cm^3 s^-1]

      ! NL solver vectors. These are cell-by-cell scratch for the ionization
      ! equilibrium solve, which now runs OpenMP-parallel over cells, so each
      ! thread needs its own copy (info is the status flag for each solve). In serial
      ! regions (input_read setup, post_process_adv) they resolve to the master
      ! thread's copy = the original behavior. The allocatables are allocated
      ! per thread inside ioniz_eq (the master's also in input_read).
      real*8, dimension(:), allocatable :: sys_sol, sys_x
      real*8, dimension(:), allocatable :: wa
      !$omp threadprivate(sys_x, sys_sol, wa, info)

      contains

      ! End of module      
      end module global_parameters

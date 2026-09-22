      module global_parameters
      ! Definition of global parameters and vectors
      
      implicit none
      
      integer, parameter :: outfile = 99  ! Unit number of report file
      ! Number of computational cells. Runtime value, set from the optional
      ! input.inp key "Grid cells: <N>" and fixed for the rest of the run;
      ! without the key it keeps the 500 that used to be compiled in, so an
      ! existing input.inp reproduces its earlier result. Every grid-sized
      ! array below is allocated by allocate_grid_arrays once N is known;
      ! grid-sized arrays inside procedures are automatic (sized on entry),
      ! exactly as the energy-grid arrays dimension(Nl) already were.
      integer :: N = 500
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
      !  38-40 OH,H2O,CO         (oxygen chemistry; zero unless
      !                           thereis_oxychem)
      ! (Si/Ca/Fe carry three stages like Mg; Na/K/S carry two stages.)
      integer, parameter :: n_species = 40
      integer :: Nl, NlTR                 ! Number of points for energy integrations
      integer :: N_eq                     ! Numbers of equations in NL solver
      integer :: lwa                      ! Working array length for NL solver
      integer :: info                     ! Output info variable of NL solver
      integer :: j_min
      ! First cell of the FLUX window: the innermost cell with r >= r_flux,
      ! over which the flux gate measures the spread of the face mass flux
      ! (flux_spread_of_state, steady_residual.f90). Set in define_grid.
      integer :: j_flux
      ! Index of the marching step the run is on: set to 0 or 1 in init and
      ! raised by one per pass of the marching loop in EXHALE_main. It counts
      ! passes, not accepted steps; n_steps_attempted and n_steps_accepted
      ! below are the two step ledgers. The run prints it as "count=" and
      ! write_setup_report keys the cap below as "count_max".
      integer :: marching_step

      ! ----- WHAT A RUN IS DOING, AND THE CLOCK THAT GOES WITH IT -----
      !
      ! Three different things share one marching loop: reaching an
      ! admissible state from a guess, relaxing to a stationary one, and
      ! advancing a state in time. Only the third has an elapsed time. Local
      ! pseudo-time (use_local_dt) advances neighboring cells by different
      ! intervals, so the material sum across an internal face is not
      ! conserved and the update is a relaxation iterate rather than one
      ! physical step; the same holds for the PTC route and for a stationary
      ! Newton finish, whose trials are numerical iterates. A run therefore
      ! states which of the two it is, and the state it writes carries that
      ! statement.
      !
      !   run_mode_init  initialization / continuation: no claim about
      !                  elapsed time; local pseudo-time, PTC and non-root
      !                  chemistry are all permitted as numerical devices
      !   run_mode_phys  physical integration: one global dt per step, and
      !                  the clock below advances only after a complete
      !                  accepted step
      integer, parameter :: run_mode_init = 1
      integer, parameter :: run_mode_phys = 2
      integer :: run_mode = run_mode_init
      ! Was the mode stated in input.inp, or derived from the other keys?
      ! A defaulted mode is a reading of the input and is reported as such.
      logical :: run_mode_given = .false.

      ! PHYSICAL ELAPSED TIME OF THE STATE THE RUN HOLDS, in SECONDS.
      ! Advanced by the global dt of a step, and only after that step has
      ! been accepted in full; a rejected trial advances nothing. Meaningless
      ! in run_mode_init, where it stays at its initial value and no output
      ! reports it. The code's time unit is R0/v0, so the increment is
      ! dt*R0/v0.
      real*8 :: t_phys = 0.0d0

      ! THE THREE COUNTERS OF CONTRACT SECTION 5, never conflated.
      ! `marching_step` above is the marching-loop index. These two separate the steps the
      ! run tried from the steps it kept: a step re-taken at half dt after a
      ! positivity violation is one more attempt and not one more accepted
      ! step, and the difference is what tells a plateau of accepted steps
      ! from a run spending its time in bisections.
      integer :: n_steps_attempted = 0
      integer :: n_steps_accepted  = 0


      ! WHICH FAMILY OF LEDGERS A DIAGNOSTIC BELONGS TO (contract section 5).
      ! The same fields are kept twice: once for initialization and
      ! continuation, where they are diagnostics of a relaxation, and once
      ! for physical integration, where they are the history of accepted
      ! steps. A quantity counted while the run was relaxing is not part of
      ! the history of a trajectory, and adding the two would produce a
      ! budget belonging to no single run state. Producers index their
      ! counters by this value; it is set by the marching loop and by the
      ! entry and exit of the stationary solver, which is continuation
      ! whatever mode the run is in.
      ! B6 CATEGORY 4, THE UNBUDGETED ACCEPTED CORRECTIONS THAT ARE STILL
      ! INSIDE THE ATTEMPTED STEP (b1 section 7.4).
      ! Row 12's Shapiro filter alters the adopted state with no source term
      ! behind the change. It is carried until B3b reaches it, and is
      ! counted here so that the certification can refuse a state whose
      ! history contains one. It lives in this module rather than in
      ! attempted_step because certification.f90 must read it and
      ! attempted_step already uses certification.
      integer :: n_shapiro_applied    = 0

      integer, parameter :: ledger_family_init = 1
      integer, parameter :: ledger_family_phys = 2
      integer :: ledger_family = ledger_family_init

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
      ! and 1e-4 for H/dr > 100. High-
      ! gravity planets have the least margin: on the default grid the value
      ! write_setup_report echoes is 26.9 cells per H at T_eq for
      ! HD 189733 b against 102.3 for WASP-121 b.
      !
      ! The default width is the double 2.0d-4 written with a d exponent, so
      ! the value is the one a reader of the key would spell: an input.inp
      ! WITHOUT the key and one stating "Base grid [dr,cells]: 2.0e-4 50"
      ! build the same grid to the bit (the list-directed read of 2.0e-4 into
      ! a real*8 is this double). The grid is part of a stored
      ! state: load_IC refuses a state whose cell centers differ from the
      ! run's by more than 1e-10 relative. An input written before 2026-09-19
      ! was run on the width 1.9999999494757503d-4 (the default-real literal
      ! 2.0e-4), whose centers differ from this grid's by 3.8e-9 to 6.6e-9
      ! relative; every such input beside stored results carries that width
      ! as "Base grid [dr,cells]: 1.9999999494757503e-4 50"
      ! (src/utils/pin_base_grid.py), and each preserved tree whose inputs
      ! were not edited has a GRID_DEFAULT_NOTE.md at its root stating the
      ! same line. The width a run used is recorded at round-trip precision
      ! by write_resolved_config.
      real*8, parameter :: dr_base_default = 2.0d-4
      real*8  :: dr_base     = dr_base_default  ! uniform base cell size [R_p]
      integer, parameter :: N_low_cells_default = 50
      integer :: N_low_cells = N_low_cells_default ! number of uniform base cells
      ! Provenance of dr_base: .true. when "Base grid [dr,cells]:" supplied
      ! the width, .false. when it still holds dr_base_default. Value
      ! equality cannot answer this, because a key may state the default's
      ! own digits; the resolved-configuration record reports the flag.
      logical :: dr_base_from_key = .false.

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

      ! ---- What stopped the marching loop -----------------------------!
      ! One flag per stop CONDITION, each true only when the condition it
      ! names actually held, so the loop's exit and the line the run prints
      ! are the same statement. The residual/flux pair is one flag because
      ! steady_gates_met is one test (steady_residual.f90).
      logical :: mass_flux_converged    = .false. ! du  < du_th: the radial
                                          !  spread of the face mass flux
      logical :: step_change_converged  = .false. ! dtu < dtu_th: the largest
                                          !  relative change of u in one step
      logical :: du_plateaued           = .false. ! du settled on a plateau
                                          !  (stall net), target NOT reached
      logical :: steady_gates_converged = .false. ! "Resid tol:": BOTH gates
                                          !  met on the marched state
      logical :: hit_max_steps          = .false. ! env EXHALE_MAXSTEPS cap
      logical :: force_start   = .false.  ! Force to do first 1000 iterations
      logical :: do_only_pp    = .false.  ! Do only the post processing
      logical :: do_load_IC    = .false.  ! Load existing IC
      logical :: thereis_He    = .false.  ! Is He included in computations
      logical :: do_read_sed   = .false.  ! Is numerical SED read from file
      logical :: is_monochr    = .false.  ! Is monochromatic radiation selected
      logical :: is_PL_sed     = .false.  ! Is the SED a power law
      logical :: thereis_Xray  = .false.  ! Include only EUV band	
      ! He metastable triplet: ON by default.  The 2^3S state carries the
      ! 10830 A observable and its channels feed the energy and electron
      ! budgets, so a run that does not say otherwise gets it.  The opt-out
      ! "Include He23S? False" exists for the deliberate HeITR-off branch
      ! check (backup/regression/wasp_he23off).  Forced back to .false. in
      ! input_read when the gas carries no helium (thereis_He false).
      logical :: thereis_HeITR = .true.   ! Include calculations for He triplet
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
      ! Atomic H/He rate set for the four reactions Koskinen et al. (2022,
      ! ApJ 929, 52) list in their Table 1 as R1-R4: radiative recombination
      ! of H+ and He+, and electron-impact ionization of H and He.
      !  .false. (default) = whichever set legacy_hhe_rates selects above.
      !  .true. ("Atomic rate set: Koskinen2022") = their Table 1 entries,
      !   so that a run can be compared like for like with their Model A.
      !   Their recombination is the Storey & Hummer (1995) power law
      !   4.0e-12 (300/T)^0.64 (H+) and 4.6e-12 (300/T)^0.64 (He+), in place
      !   of the default Badnell RR minus Mao & Kaastra alpha_1 case B; their
      !   collisional ionization is the same Voronov (1997) fit EXHALE
      !   already uses by default, so R3/R4 move nothing unless
      !   legacy_hhe_rates is set as well. The default set is the physically
      !   preferred one here (case B, the Lyman continuum being optically
      !   thick), so this key exists to reproduce their choice, not to
      !   replace ours. It swaps the RATE coefficients only: the
      !   recombination COOLING rates are untouched.
      logical :: atomic_rate_set_k22 = .false.
      ! "Caloric EOS: monatomic" -- every particle, H2 included, stores
      ! (3/2) k T (gamma = 5/3 everywhere), the state of the code before the
      ! H2 rovibrational ladder of section P53. A COMPARISON option: a model
      ! whose energy equation is u = c_v T with a monatomic c_v cannot be
      ! matched with the ladder on. Default .false. = the ladder.
      logical :: caloric_eos_monatomic = .false.
      ! "Photoelectron heating: full" -- every photoionization deposits the
      ! WHOLE photon energy h nu as heat instead of the photoelectron's
      ! h nu - I (the ionization energy I is then counted twice: once here
      ! and once when it is radiated away by recombination). Unphysical, and
      ! a COMPARISON option only: Koskinen et al. (2022) Figure 9 shows a
      ! stellar heating rate 2.5-3 times what h nu - I gives for their own
      ! ionization rate and spectrum, and this option tests whether that
      ! accounting reproduces their profile. Default .false. = h nu - I.
      logical :: photoheat_full_photon_energy = .false.
      ! "Photoelectron heating: <f>" -- a stated FRACTION f of the photon
      ! energy per ionization is deposited as heat (f = 1 is "full"; the
      ! physical accounting h nu - I is f < 0, the default). Comparison only:
      ! the Koskinen et al. (2022) Model A heating profile sits between the
      ! two accountings, at about f = 0.65, and this reproduces it without
      ! claiming a physics.
      real*8  :: photoheat_photon_fraction = -1.0d0
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
      ! Step at which sec_ion_active was last switched on (-1 = it was never
      ! staged in, i.e. it is either off or was on from step 0). Written into
      ! the restart file's coupling header and read back from it.
      integer :: sec_ion_armed_step = -1
      ! Input override "Secondary_ionization: Immediate": apply the coupling from
      ! step 0 (pre-staging behavior), for A/B tests only.
      logical :: sec_ion_immediate = .false.
      ! He recombination radiation ionizing H (Draine 2011 y/z parametrization,
      ! on-the-spot).
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
      ! He/H diffusive separation:
      !  .false. (default) = He/H frozen at the input HeH everywhere (legacy);
      !  .true. = evolve the He element ratio with advection + molecular
      !  diffusion (He settles, He/H falls with altitude).  Metals stay frozen
      !  to H.  Off = byte-identical to legacy.
      logical :: he_diffusion = .false.
      ! Eddy (turbulent) diffusion coefficient [cm^2/s] used with he_diffusion.
      ! Mixes the composition toward a uniform MASS fraction below the
      ! homopause (where the molecular D_12 equals he_kzz) and carries no
      ! settling term of its own; above it the elements separate.  Default 0,
      ! i.e. pure molecular diffusion: an eddy term is a property of the
      ! atmosphere being modelled, so it is stated, not inherited.  Taylor et
      ! al. (2025) use K_zz = 1e5 m^2/s = 1e9 cm^2/s.  Runtime key
      ! "He_Kzz: <value>"; base.inp may override it (Kzz_base).  A run that
      ! carries a lower-atmosphere profile takes K_zz from the profile
      ! instead and this scalar is inert.
      real*8  :: he_kzz = 0.0d0
      ! Eddy diffusion coefficient of every cell [cm^2/s].  This is the array
      ! the element-diffusion operator reads; he_kzz is only the constant a
      ! run states when it has no profile.  Filled once the radial grid
      ! exists (eddy_diffusion_on_grid, lower_atmosphere_profile.f90): from
      ! the interpolated "Kzz" column of a lower-atmosphere profile when one
      ! is given, and from he_kzz in every cell otherwise.  The uniform case
      ! is exactly the old scalar: 0.5*(a+a) = a in IEEE double.
      real*8, dimension(:), allocatable :: kzz_cell
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
      ! "q_H2_base" key of base.inp (or by the profile at its matching
      ! level).  This is the same quantity the chemical-equilibrium fit
      ! q_h2_equilibrium returns, so no conversion is involved.
      !
      ! It is the COMPOSITION OF THE INFLOWING GAS, not only an equation-of-
      ! state anchor (section 117 of docs/Update_EXHALE.*): when supplied it
      ! replaces the fit in the molecular-base particle count AND is imposed
      ! on the H2 partition of the lower ghost species, with the base
      ! particle count then taken from that same species state.  The two
      ! therefore describe one gas.  The value is upstream information: the
      ! lower atmosphere is where photodissociation and mixing set the
      ! partition, and the shielded base cell cannot derive it.
      !
      ! Must not exceed 0.5/(0.5 + HeH), the fully molecular limit; startup
      ! refuses a larger value instead of capping it.  Negative (default) =
      ! no handoff value, the fit is used and the species partition is left
      ! free, so a base.inp without the key behaves as before.
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
      ! Band-integrated stellar flux in the H2 Lyman-Werner bands
      ! (912-1201 A, i.e. 10.3-13.6 eV) AT THE PLANET'S ORBIT
      ! [erg cm^-2 s^-1].  Key "Stellar LW flux [erg cm^-2 s^-1]: <F>".
      ! It is a separate input because the code's own energy grid does not
      ! carry this band: for a numerical SED the grid stops at the 13.6 eV
      ! HI edge (only the He 2^3S / low-IP-metal thresholds push it lower),
      ! and for a power-law SED anything below 13.6 eV is an extrapolation
      ! of the XUV law, which has nothing to do with a star's FUV.
      ! 0 (default) = no Lyman-Werner photodissociation, i.e. the molecular
      ! network as it was before this key existed, bit for bit.
      !
      ! THE KEY IS THE FLUX AT THE PLANET, NOT THE FLUX THE MOLECULES SEE.
      ! The run-wide dayside convention of "2D approximate method" is applied
      ! to it where the beam is used, by dayside_dilution() in fuv_band_flux,
      ! exactly as it is to the XUV grid and to the stellar Ly-alpha beam.
      ! State the band flux at the orbit and let the run dilute it; do not
      ! pre-divide it here (Update_EXHALE_stage1 section 150).
      ! See src/modules/lower_atmosphere/lyman_werner.f90.
      real*8  :: F_LW_star = 0.0d0
      ! Whether input.inp STATED "Stellar LW flux" (any value, zero
      ! included). A stated zero switches the band off even when a spectrum
      ! file could supply it -- the comparison with a published model that
      ! excludes H2 photodissociation needs exactly that.
      logical :: lw_flux_stated = .false.
      ! Was F_LW_star computed from the numerical spectrum rather than
      ! stated by the run?  Reported, so the source of the number is visible.
      logical :: lw_from_spectrum = .false.
      ! The same two statements for the FUV continuum bands B3 (1231-1450 A)
      ! and B4 (1451-2304 A) of the oxygen photolysis: whether input.inp
      ! STATED the key (a stated zero switches the band off), and whether
      ! the value was integrated from the spectrum file instead. Band B2 is
      ! the Ly-alpha line and is always stated ("Stellar Lya flux").
      logical :: fuv_b3_flux_stated = .false., fuv_b4_flux_stated = .false.
      logical :: fuv_b3_from_spectrum = .false., fuv_b4_from_spectrum = .false.
      ! Oxygen chemistry (the A2 option):
      ! OH / H2O / CO added to the coupled molecular ionization equilibrium,
      ! with the H2O and OH photolysis of the FUV bands.  It is what lets
      ! EXHALE compute its own base H2/H partition instead of importing it
      ! through q_H2_base or taking it from the chemical-equilibrium fit.
      ! Key "Oxygen chemistry: True". Requires the molecular network, helium
      ! and a non-zero oxygen abundance (input_read refuses otherwise).
      ! Default off = byte-identical to a run without the key.
      logical :: thereis_oxychem = .false.
      ! Vertical transport of the molecular carriers -- H2 always, plus OH,
      ! H2O and CO when the oxygen cycle is on -- solved implicitly with
      ! their chemistry by diffusive_photochemistry.
      !
      ! Key "Molecular carrier transport: True|False". Its DEFAULT is not a
      ! constant: it is on whenever the oxygen chemistry is on, because a
      ! local steady state is the wrong physics at the cool base that option
      ! exists for -- tau_chem(H2)/tau_adv is of order unity there. Resolved
      ! in input_read
      ! once every key is parsed; carrier_transport_stated records whether
      ! the run said so itself, so the default can change without silently
      ! overriding a stated value.
      !
      ! False restores the local-kinetics limit of milestone M2, which is a
      ! test of the chemistry alone and not a model of a base.
      logical :: carrier_transport = .false.
      ! Solve the carrier continuity equation TOGETHER with the wind, as a
      ! fourth Newton unknown per cell, instead of alternating the two
      ! (section 139). Meaningless without carrier_transport, and default
      ! off: the coupled route changes the size and the band geometry of the
      ! steady system, so a run that does not ask for it must not pay for it.
      logical :: carrier_in_newton = .false.
      ! THE CARRIER GATE. When the carrier row counts as steady: 0.1 percent
      ! of the row's own largest terms, volume-weighted over the layer and
      ! the wind separately (carrier_steady_residual). It is a THIRD gate
      ! beside the residual and the flux gate of section 133, and not part of
      ! ||R||, for the reason section 133 gives for keeping those two apart:
      ! the numbers are not commensurable. `Resid tol` was renormalized in
      ! section 133 against a row scale that BOUNDS each hydrodynamic row's
      ! largest term rather than being it, and converged states sit at 1e-6
      ! on it; the carrier row is measured against the terms themselves, so
      ! a converged carrier row sits near 1e-3. One threshold cannot serve
      ! both. Read by the Picard loop and by steady_gates_met.
      real*8  ::  carrier_resid_th = 1.0d-3
      logical :: carrier_transport_stated = .false.
      ! TRANSPORT THE HYDROGEN IONIZATION STATE: H+ carried with the flow as
      ! a fifth transported species, instead of being re-solved every step
      ! as a local photoionization/recombination equilibrium.
      !
      ! Key "Ionization transport: True|False", DEFAULT FALSE.
      !
      ! WHY IT EXISTS, AND WHERE THE LOCAL EQUILIBRIUM IS WRONG. The local
      ! partition is the right answer only where a parcel is ionized faster
      ! than it leaves the shell it sits in, P r/|v| >> 1. Measured on the
      ! Koskinen 2022 Model A comparison (a 0.0457 M_J planet at 0.048 au
      ! with the flux quartered, docs/k22_electron_density_excess.md sec. 7)
      ! that number is 0.15-0.35 above 1.5 r_base: the gas leaves each shell
      ! three to seven times faster than it can be photoionized, so the
      ! ionization fraction is not the local root but whatever the parcel
      ! accumulated on the way up. The local equilibrium then gives
      ! x(H+) = 0.129/0.327/0.585 at 2.0/2.4/3.0 r_base where integrating
      ! dx/dt along the same flow gives 0.057/0.078/0.106 and Model A, which
      ! advects every species, has 0.020/0.055/0.080.
      !
      ! It is OFF by default because for the hot Jupiters the local closure
      ! was built for P r/|v| is large and the equilibrium holds, and because
      ! turning it on changes every ionization-dependent number of a run.
      ! Requires "Molecular carrier transport: True": the operator that
      ! carries it is the molecular carrier solve.
      logical :: ionization_transport = .false.
      ! Band-integrated stellar flux AT THE PLANET'S ORBIT in the two FUV
      ! continuum bands of the oxygen chemistry [erg cm^-2 s^-1]:
      !   B3 1231-1450 A   "Stellar FUV B3 flux [erg/cm2/s]: <F>"
      !   B4 1451-2304 A   "Stellar FUV B4 flux [erg/cm2/s]: <F>"
      ! B2 is the Ly-alpha line, supplied by the existing F_Lya_star, and
      ! the 912-1201 A band is supplied by F_LW_star, which the H2
      ! Lyman-Werner absorber shares with H2O and OH -- one wavelength
      ! interval, one incident flux, one beam. There was a fourth key,
      ! "Stellar FUV B1 flux", for a 1110-1201 A band between it and B2; on
      ! 2026-09-06 that band was merged back into the Lyman-Werner interval,
      ! because the H2 lines pump across 1110 A at the temperature of a
      ! planetary base and a band edge there normalized the pumping per
      ! photon of a band narrower than the one the lines drink from. The key
      ! is retired and input_read stops a run that states it. Separate keys
      ! rather than one flux plus an assumed shape: decision D1 of the
      ! design, taken because the single-key form is weakest exactly on the
      ! line-dominated FUV of an M dwarf. The band edges are fixed by the H2O
      ! branching ratios and by the Ly-alpha line, not chosen; see
      ! src/modules/lower_atmosphere/water_photolysis.f90.
      ! 0 (default) = no photolysis in that band.
      ! Like F_LW_star, these are the fluxes AT THE ORBIT; the dayside
      ! convention is applied in fuv_band_flux (section 150).
      real*8  :: F_FUV_B3 = 0.0d0
      real*8  :: F_FUV_B4 = 0.0d0
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
      ! Roll-off width of the coronal-excitation guard, in units of the
      ! fractional temperature deficit below the 1e3 K CHIANTI fit floor
      ! ("Coronal cutoff width: <w>"; see coronal_excitation_cutoff in
      ! Cool_coeff.f90). The guard is exp(-x^2) with x = (T_floor/T - 1)/w, so
      ! w is a modeling choice, not a measured quantity: the temperature the
      ! base settles at depends on it at the ~100 K level. The measured
      ! justified window is w = 0.08-0.13 (docs/coronal_cutoff_width.md);
      ! 0.1 is the default.
      real*8  :: coronal_cutoff_width = 0.1d0
      ! Thermal infrared field of the lower atmosphere, seen by the
      ! ground-term fine-structure lines ("Base IR field: T|F"). The gas
      ! below the base is optically thick in [C I] 609/370um and
      ! [O I] 63/145um (tau = 1 within 1-3 pressure scale heights below the
      ! base), so it radiates a diluted blackbody at T0 into the lower
      ! hemisphere of every cell. .false. (default) leaves those lines with
      ! no incident field, i.e. the lower atmosphere is treated as cold and
      ! the lines can only cool; .true. solves the ground term with the
      ! incident field in it, which gives a radiative-equilibrium floor near
      ! the temperature where B_nu(T) = W B_nu(T0). See
      ! fine_structure_line_transfer in Cool_coeff.f90.
      logical :: base_ir_field  = .false.
      ! Molecular infrared bands (`Molecular IR bands`). .false. is the state
      ! before 2026-08-30: the only infrared coolants below the H2 -> H front
      ! were H3+ and the ground-term fine-structure lines, so a converged
      ! molecular layer had nothing holding it and collapsed to 190-400 K.
      ! .true. adds the H2 quadrupole and magnetic
      ! dipole line spectrum and the H2O and CO bands, each exchanging with the
      ! same diluted B_nu(T0) `Base IR field` supplies, so each stops at its own
      ! radiative equilibrium temperature. See molecular_infrared_cooling.f90.
      logical :: mol_ir_bands   = .false.
      ! H2 photoabsorption channels beyond the single dissociative one.
      ! Both default to the pre-E1 arithmetic, which is why an old run
      ! reproduces bit for bit.
      !
      ! h2_double_ionization selects the model for the share of the
      ! photo-released protons that comes from H2 + hv -> H+ + H+ + 2e-,
      ! a channel with a vertical threshold of 51.4 eV (Yan, Sadeghpour &
      ! Dalgarno 1998, sec. 4). Accepted values:
      !   'chung80'  DEFAULT. Chung, Lee, Masuoka & Samson (1993),
      !              discussion of their Fig. 4: "by 80 eV about 20% of
      !              sigma(H+) comes from double ionization. This
      !              percentage remains fairly constant towards higher
      !              energies", ramped from the 51.4 eV threshold.
      !   'yan_rho'  built from the double-to-single ionization ratio of
      !              Yan et al. (1998) sec. 4 (0.038 at 110 eV, asymptotic
      !              0.0225 from Sadeghpour & Dalgarno 1993).
      !   'off'      the channel is folded back into the single
      !              dissociative one, d(E) = 0, which is what the code did
      !              before this channel existed. Use it to reproduce a
      !              pre-E1 result, not as a physical statement: the
      !              reaction happens.
      !
      ! WHY 'chung80' IS THE DEFAULT. The reaction is real and measured;
      ! what is not measured as a function of energy is d(E) itself. The
      ! 20 percent is the number the source states outright, so it is the
      ! one that carries. RANGE: the ramp between 51.4 and 80 eV is NOT in
      ! the source, which says only that the contribution is "small" near
      ! threshold. UNCERTAINTY: 'yan_rho' rests on different published
      ! numbers and gives a double-to-single ratio 1.6 times larger at
      ! 110 eV (0.038 against 0.023); running both is how that is reported,
      ! and it moves n(H2+) by 25 percent of the channel's own effect.
      character(len=16) :: h2_double_ionization = 'chung80'
      ! h2_neutral_dissociation resolves H2 + hv -> H + H, the absorptions
      ! that leave no ion. It is nonzero ONLY over the 33-41 eV window of
      ! Chung et al. (1993) Table 1, where the measured photoionization
      ! yield falls below unity; outside that window the source ASSUMES a
      ! unit yield rather than measuring one, so the channel is set to
      ! zero there.
      !
      ! DEFAULT True: sigma_n is a MEASUREMENT (Chung et al. Table 1,
      ! twelve rows), and without it up to 7.4 percent of the absorptions
      ! at 37.5 eV are given an H2+ and an electron that the event does not
      ! make. Setting it False folds that share back into the ionizing
      ! channels in their own proportion, which is what the code did before
      ! this channel existed; use it to reproduce a pre-E1 result, not as a
      ! physical statement.
      logical :: h2_neutral_dissociation = .true.
      ! `Molecular reaction heat: True|False`.  The energy the COLLISIONAL
      ! reactions of the H2/He network release into the gas -- above all the
      ! dissociative recombination of H3+ and H2+, which returns the H2
      ! ionization energy to the gas instead of to a photon as radiative
      ! recombination does in an atomic gas.
      !
      ! DEFAULT ON, because leaving it out is physically wrong rather than
      ! approximate.  Follow the closed cycle a photon drives:
      !     H2 + hv -> H2+ + e        (photon pays I(H2) = 15.43 eV; the
      !                                photoelectron keeps hv - I(H2), which
      !                                PH_heat_HHe already deposits)
      !     H2+ + H2 -> H3+ + H       (+1.70 eV)
      !     H3+ + e  -> H2 + H        (+9.25 eV)
      ! The cycle returns to H2 having converted one H2 into H + H and left
      ! I(H2) - D0(H2) = 15.43 - 4.48 = 10.95 eV in the gas as kinetic energy
      ! of the fragments -- NOT as a photon, which is what makes this
      ! different from the atomic case, where radiative recombination carries
      ! the ionization energy out of the gas.  Measured on the converged
      ! He/H = 0.0793 hot Uranus, that sum is 81 percent of the total heating
      ! rate at 1.02 r_base, so omitting it is not a small error in the
      ! molecular layer's energy budget; it is most of it.
      !
      ! An atomic run is untouched: the term is gated on thereis_mol as well,
      ! so every non-molecular golden is byte-identical either way.
      ! See molecular_reaction_heat.f90.
      logical :: mol_reaction_heat = .true.
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
      ! ---- well-balanced pressure/gravity pair ("Well balanced:") ----
      ! .true. carries the DEPARTURE from a local hydrostatic equilibrium
      ! through the reconstruction, the Riemann jumps and the pressure force,
      ! instead of the state itself: within cell j the equilibrium is the one
      ! of constant density,
      !     p_eq,j(r) = p_j - rho_j (phi(r) - phi(r_j)),
      ! which is the mechanical balance dp/dr = -rho dphi/dr to second order
      ! in the cell width and assumes NO thermal stratification, so an
      ! arbitrary entropy and composition profile is preserved (Kaeppeli and
      ! Mishra 2016, A&A 587, A94, their sections 2.1.1 and 2.1.3; the
      ! face-pressure-difference form of the momentum source is their 2014
      ! paper's eq. 2.26, J. Comput. Phys. 259, 199).  The equilibrium's flux
      ! difference and its source then cancel ANALYTICALLY rather than in
      ! floating point, and the Riemann dissipation acts on the departure.
      ! Exact preservation needs a numerical flux that resolves a stationary
      ! contact discontinuity (ROE and HLLC do; LLF does not).
      ! Default .false.: with it the operator is the one every golden was
      ! taken with, to the bit.
      logical :: well_balanced = .false.
      logical :: recon_two_stage = .false.  ! "Reconstruction scheme: PLM+WENO3": run PLM (stage 1) then WENO3 (stage 2), using BOTH du_th values. PLM/WENO3 alone are single-stage and use only the first du_th value.
      ! ---- PLM -> WENO3 continuation ("Reconstruction continuation:") ----
      ! The two-stage recipe changes the discrete operator in ONE step: the
      ! face reconstruction, the form of the momentum equation (pressure
      ! inside the flux plus a geometric source, versus a face-pressure
      ! difference) and the outer free-outflow ghost all change together, and
      ! the state marched to that point under PLM is handed to WENO3 as if it
      ! were a solution of it. Measured at the hand-off of the hot-Uranus
      ! molecular case, the two operators' steady residuals differ by 6.7
      ! times the largest PLM residual in the mass row at the base and by 62
      ! percent in the momentum row at the outer boundary, while agreeing to
      ! a part in 1e2 to 1e4 through the wind.
      !
      ! The continuation replaces the jump by a homotopy between the two
      ! discretizations,
      !
      !     R_lambda(u) = (1 - lambda) R_PLM(u) + lambda R_WENO3(u),
      !
      ! assembled by reconstruction_continuation_rhs (steady_residual.f90) and
      ! walked from lambda = 0 to lambda = 1 with step control. The
      ! conservative ghost fill (Apply_BC_W) is blended by the same lambda, so
      ! lambda = 0 reproduces the pure PLM operator and lambda = 1 the pure
      ! WENO3 operator to the bit, and the endpoints cost one right-hand side
      ! rather than two. Once lambda = 1 is reached the continuation disarms
      ! itself and the flags are set to WENO3, so the JFNK finish and every
      ! acceptance gate see the production operator with no blending in it.
      !
      ! recon_lambda_step0 <= 0 (the default, and what an absent key leaves)
      ! is the shipped one-step switch. See docs/input_schema.md.
      logical :: recon_lambda_on   = .false. ! continuation armed for this run
      real*8  :: recon_lambda      = 0.0d0   ! current lambda in [0,1]
      real*8  :: recon_lambda_step = 0.0d0   ! current step in lambda
      real*8  :: recon_lambda_step0 = 0.0d0  ! "Reconstruction continuation:" dlambda; <= 0 = off
      real*8  :: recon_lambda_dtu_tol = 1.2d0 ! lambda advances while dtu <= tol * dtu at the ramp start
      logical :: recon_lambda_adaptive = .true. ! .false. = ramp with no step control
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
                                          !  sonic-point topology);
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
      ! 2x-band reject ("Brent solver: False").
      logical :: use_newton_ieq  = .true.
      logical :: use_brent_tsolve = .true.

      !------- Global constants -------!
      
      ! Physical constants
      real*8,parameter ::  pi      = 3.1415926536d0   ! pi
      ! Boltzmann constant in CGS units (CODATA/SI exact value 1.380649e-16;
      ! updated 2026-08-15 from the truncated ATES literal 1.38e-16, a 4.7e-4
      ! relative change that moves every thermal quantity -- goldens were
      ! re-snapshotted with this value).
      real*8,parameter ::  kb_erg  = 1.380649d-16
      real*8,parameter ::  kb_eV   = 8.617333262d-05  ! Boltzmann constant (eV/K), CODATA exact
      !
      ! ----- Particle masses -----
      !
      ! ONE DEFINITION EACH, and every mass the code uses is either one of
      ! these measured quantities or a ratio formed from them below.  The
      ! MASS UNIT of the code is the hydrogen ATOM mu: the adimensional
      ! density rho is a number of hydrogen atoms per unit volume times n0,
      ! so a species mass entering rho is m_species/mu and NOT its atomic
      ! weight in u.  The two differ by 0.78 percent (amu_over_m_H), which
      ! is why the conversion is written here once instead of being carried
      ! implicitly by whichever unit a table happened to be transcribed in.
      !
      ! Mass of the hydrogen ATOM (m_p + m_e - 13.6 eV/c^2), CODATA 2018.
      ! This is the mass unit of the density normalization, so it also sets
      ! m_H2 in lyman_werner.f90 -- keep the two in step.
      real*8,parameter ::  mu      = 1.67353284d-24   ! Hydrogen atom mass (g)
      ! Mass of the helium-4 ATOM: 4.002603254 u (AME2020 atomic mass of
      ! 4He) times the atomic mass unit below, equivalently the alpha
      ! particle plus two electrons less the 79 eV electronic binding
      ! (2.6e-8 of the mass, below the digits kept).  It is the mass the
      ! species table gives He I, He II, He III and the He 2^3S metastable:
      ! removing one electron changes the atom's mass by 1.4e-4, below the
      ! precision at which any helium mass enters this code.
      real*8,parameter ::  m_He_atom = 6.6464790722d-24 ! Helium-4 atom mass (g)
      ! Unified atomic mass unit, 1/12 of the mass of a neutral 12C atom,
      ! CODATA 2018.  Standard atomic weights (species_table melem_A_u) are
      ! quoted in this unit and must be converted before they are added to a
      ! mass sum carried in hydrogen atoms.
      real*8,parameter ::  amu     = 1.66053906660d-24 ! atomic mass unit (g), CODATA 2018
      real*8,parameter ::  m_e     = 9.1093837015d-28  ! electron mass (g), CODATA 2018
      real*8,parameter ::  m_p     = 1.67262192369d-24 ! proton mass (g), CODATA 2018
      ! Derived mass ratios, formed from the constants above so that a mass
      ! ratio cannot drift from the masses it is a ratio of.
      !   m_He_over_m_H = 3.9715259 : the weight of one helium nucleus in
      !     the species table (bsp_mass), in calc_rho, in the mean molecular
      !     weight, and in the analytic lower column.
      !   amu_over_m_H  = 0.99223573 : the factor that turns an atomic
      !     weight in u into the code's mass unit.
      real*8,parameter ::  m_He_over_m_H = m_He_atom/mu
      real*8,parameter ::  amu_over_m_H  = amu/mu
      ! Adiabatic index of a monatomic gas, written as the exact rational so
      ! that the value is 5/3 to full double precision rather than to the
      ! digits a literal happens to carry.  This is the reference value: where
      ! molecules are present the caloric EOS uses the mixture value
      ! gamma_eff(T, composition) of module caloric_eos instead, and this
      ! constant is what that module returns, verbatim, for a molecule-free
      ! cell.
      real*8,parameter ::  gamma_ad = 5.0d0/3.0d0      ! monatomic reference
      real*8,parameter ::  Gc      = 6.67430d-8       ! Gravitational constant (CGS), CODATA 2018
      real*8,parameter ::  erg2eV  = 6.241509075d11   ! 1 erg measured in eV
      real*8,parameter ::  hp_erg  = 6.62607015d-27   ! Planck constant (CGS), CODATA exact
      real*8,parameter ::  hp_eV   = 4.135667696d-15  ! Planck constant (eV*s), CODATA exact
      real*8,parameter ::  c_light = 2.99792458d10    ! Speed of light in cm/s
      real*8,parameter ::  parsec  = 3.08567758147d18 ! 1 pc in cm
      real*8,parameter ::  AU      = 1.495978707d13   ! Astronomical unit
      ! Jupiter and solar mass/radius: the IAU 2015 nominal values (Prsa et
      ! al. 2016, AJ 152, 41, Table 1), which are the units transiting-planet
      ! radii and masses are quoted in. R_J is the nominal EQUATORIAL radius
      ! R_J^N(eq), the one a transit depth measures; the masses follow from
      ! the nominal GM with the CODATA G above. These constants scale the
      ! `Planet radius [R_J]`, `Planet mass [M_J]` and `Parent star mass
      ! [M_sun]` keys, and are the DEFINITION the Wind-AE front end
      ! (wae_exhale_input.f90), which cannot use this module, mirrors.
      real*8,parameter ::  RJ      = 7.1492d9         ! Jupiter equatorial radius (cm), IAU 2015 nominal
      real*8,parameter ::  MJ      = 1.8982d30        ! Jupiter mass (g), IAU 2015 nominal
      real*8,parameter ::  Msun    = 1.98842d33       ! Sun mass (g), IAU 2015 nominal
      real*8,parameter ::  Rsun    = 6.957d10         ! Sun radius (cm)
      real*8,parameter ::  R_earth = 6.3725d8         ! Earth radius (cm)
      real*8,parameter ::  M_earth = 5.9726d27        ! Earth mass (g)
      real*8,parameter ::  ih      = 1.0d0            ! Atomic number of Hydrogen
      real*8,parameter ::  ihe     = 2.0d0            ! Atomic number of Helium
	
	   ! -- - Energy constants

      ! Energy intervals
      real*8 ::  e_top
      real*8 ::  e_mid
      real*8 ::  e_low
      
      ! ----- Ionization thresholds -----
      !
      ! ONE DEFINITION EACH.  Every consumer of an ionization threshold reads
      ! the constant below: the photon-grid band edges (set_energy_vectors),
      ! the photoelectron energy h nu - e_th, the secondary-ionization targets
      ! (electron_energy_degradation), the collisional-ionization cooling
      ! (util_ion_eq, T_equation) and the turn-on of the photoionization cross
      ! sections themselves (cross_sec).  A cross section that turned on at
      ! its own copy of the threshold left the band between the two copies
      ! integrated as zero: with e_th_HeI = 24.6 against the He I fit's own
      ! 24.59, the band [24.59, 24.60] eV, 0.115 percent of the He I
      ! photoionization rate of the default power law, was lost.
      !
      ! The values are the measured ionization energies, NIST Atomic Spectra
      ! Database (Kramida, Ralchenko, Reader & NIST ASD Team), levels and
      ! ionization-energy data, and for H2 the adiabatic ionization energy of
      ! the NIST Chemistry WebBook (Herzberg & Jungen 1972 series limit).
      real*8,parameter ::  e_th_HI   = 13.598434599d0  ! H I  1s -> H+ , NIST ASD [eV]
      ! H2 -> H2+ + e-, adiabatic (v'=0 <- v''=0) ionization energy, NIST.
      real*8, parameter ::  e_th_H2  = 15.425927d0     ! H2 photoionization threshold [eV]
      ! Threshold of the DISSOCIATIVE H2 photoionization channel,
      ! H2 + hv -> H + H+ + e- (Chung, Lee, Masuoka & Samson 1993,
      ! J. Chem. Phys. 99, 885, Table II, whose first row is the
      ! threshold 18.076 eV). The non-dissociative channel above
      ! it keeps e_th_H2; the two share one cross section, split by
      ! frac_H2_dissociative_ionization.
      real*8, parameter ::  e_th_H2_di = 18.076d0
      ! Threshold of the DOUBLE ionization channel,
      ! H2 + hv -> H+ + H+ + 2e- . Yan, Sadeghpour & Dalgarno (1998)
      ! sec. 4: the channel "has a vertical threshold of 51.4 eV".
      ! Charged against that threshold the same way every other channel
      ! is charged its own ionization potential.
      real*8, parameter ::  e_th_H2_dd = 51.400d0
      real*8,parameter ::  e_th_HeI  = 24.587389d0    ! He I 1^1S -> He+ , NIST ASD [eV]
      real*8,parameter ::  e_th_HeII = 54.417765d0    ! He II 1s -> He++ , NIST ASD [eV]
      ! He 2^3S metastable: the ionization energy of the ground singlet minus
      ! the excitation energy of the 1s2s 3S1 level, 159855.9743 cm^-1 =
      ! 19.819614 eV (NIST ASD), i.e. 24.587389 - 19.819614 eV.  This is the
      ! lowest threshold in the code and it sets the floor of the photon grid
      ! whenever the metastable is carried.
      real*8,parameter ::  e_th_HeTR = 4.767775d0     ! He 2^3S -> He+ [eV]
      ! The same thresholds IN ERG -- the energy the electron gas loses per
      ! collisional ionization. They are the eV constants above divided by
      ! erg2eV, so a threshold is written once and the collisional cooling
      ! assembly (util_ion_eq eval_cool) and the cell-by-cell temperature root
      ! of the advection post-process (T_equation) charge the same energy per
      ! event as the photon grid charges per photoionization.
      ! Values: 2.178709e-11, 3.939334e-11, 8.718687e-11 and 7.638818e-12 erg.
      ! (The Hui & Gnedin 1997, MNRAS 292, 27 collisional-ionization fit for
      ! He II carries the same threshold in its exponent as 2*631515/T; that
      ! fit parameter stays as published in Cool_coeff, and 631515 K is
      ! 0.0037 percent above the measured potential.)
      real*8,parameter ::  e_th_HI_erg   = e_th_HI/erg2eV
      real*8,parameter ::  e_th_HeI_erg  = e_th_HeI/erg2eV
      real*8,parameter ::  e_th_HeII_erg = e_th_HeII/erg2eV
      real*8,parameter ::  e_th_HeTR_erg = e_th_HeTR/erg2eV
      ! Photoelectron energy above which the energy partition of
      ! electron_energy_degradation is applied; below it the photoelectron is
      ! taken to deposit all of its energy as heat.
      !
      ! 30 eV is the LOWEST PRIMARY ENERGY AT WHICH ANY COEFFICIENT OF THAT
      ! PARTITION IS PUBLISHED: Dalgarno, Yan & Liu (1999) tabulate 30, 50,
      ! 100, 200, 500 and 1000 eV, so below 30 eV every coefficient would be
      ! an extrapolation off the end of their tables. (It was 40 eV, carried
      ! over from the wind_ae X-ray cutoff -- a numerical convention of
      ! another code, not a physical threshold, and it sat between the two
      ! sources rather than on either.)
      !
      ! WHAT IS STILL DISCARDED, AND WHY. A photoelectron of 13.6-30 eV can
      ! ionize hydrogen, and in the hot-Uranus band those photons carry 18%
      ! of the incident XUV energy; here they are still treated as depositing
      ! all of their energy as heat. That is a limit of the tables, not a
      ! physical threshold: no source consulted resolves the partition below
      ! 30 eV, and extrapolating into it would be unsupported.
      real*8, parameter :: E_sec_ion = 30.0d0
      real*8,parameter ::  e_th_MgI  = 7.646d0   ! Threshold for MgI ionization
      real*8,parameter ::  e_th_MgII = 15.035d0  ! Threshold for MgII ionization

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
      ! Hard cap on marching iterations, settable as "Max steps: <N>".
      ! (The env variable EXHALE_MAXSTEPS is a separate, lower deterministic
      ! cap used by the regression harness; it exits the loop without the
      ! "reached count_max" message and never raises this one.)
      integer :: marching_step_max = 1000000   ! printed as "count_max"
      integer :: N_stall   = 2000     ! window of stalled steps; settable via
                                      !   input line "Stall [tol,N]: <tol> <N>"
      real*8  :: stall_tol = 1.0d-6   ! rel. du change defining a stall (same line)

      ! Mass-flux LEVEL stability tolerance: in stage 2 a run only counts as
      ! converged/stalled when the mean |rho*v*r^2| over [j_min:N] changed by
      ! less than lev_th (relative) across the last N_stall steps. du is the
      ! SPREAD of the mass flux and is blind to a uniform drift of the whole
      ! profile; without this gate, runs can stop on spatially flat but still
      ! level-drifting states that differ by tens of percent in Mdot
      ! (path-dependent quasi-steady snapshots).
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
      ! finishes the run to ||R|| < resid_th (default 1e-5 if "Resid tol:" was
      ! not given, EXHALE_main) and the standard final outputs /
      ! post-processing follow.  The base valve of the early versions is
      ! retired ("Valve eps:" is refused by input_read); the base boundary is
      ! the characteristic face state of base_boundary.f90.
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
      ! ACCEPTANCE IS MEASURED AT THE STATE'S OWN COMPOSITION (section 155).
      ! At hand-back the residual is re-evaluated from the composition the
      ! RETURNED state carries -- iterating the sweep at fixed Y until it stops
      ! moving -- and the gate is applied to that.
      !
      ! WHY THIS IS THE DEFAULT. `eval_residual` is a function of two things,
      ! the unknowns and the composition the equilibrium sweep starts from
      ! (section 147 made that explicit in its interface). The line search
      ! measures a trial with the PREVIOUS iterate's composition, which is the
      ! right thing to do while comparing trials; but the number then reported
      ! as the accepted root's residual is a residual the returned state does
      ! not have. Measured: at its own composition it is 4 to 56 times larger
      ! (hot Uranus 2.2e-6 -> 1.3e-5; WASP-121b 6.3e-6 -> 3.9e-5; and 8.9e-6 ->
      ! 1.9e-4 on `wasp_full_newton` before section 152). A state that is
      ! accepted must satisfy its own residual, so this is not an option to be
      ! opted into -- turning it OFF is what needs a reason, and the only one
      ! is reproducing a pre-section-155 number.
      !
      ! WHAT IT COSTS. One extra sweep sequence at hand-back, 3 iterations on
      ! the hot Uranus and 9 on WASP-121b, once per solve. `Resid tol` values
      ! calibrated against the old measure are correspondingly loose under this
      ! one; the flux gate is unaffected, reading conserved face fluxes that do
      ! not depend on the composition at all.
      logical :: resid_at_own_composition = .true.
      ! Cap on the fixed-Y sweep iteration, and A WORKING LIMIT rather than a
      ! guard. MEASURED by B5a on wasp_full_newton (2026-09-06): the sweep
      ! takes 13 to 15 passes from the marching hand-off and 14.09 on
      ! average per residual evaluation over the whole solve, so the margin
      ! to this cap is under a factor of two, not the factor of three to
      ! eight an earlier note claimed. An evaluation that reaches the cap is
      ! reported rather than silently accepted, and the number to watch when
      ! a solve slows down is this one.
      integer :: n_selfconsistent_max = 25

      ! Particle density (n_tot + n_e) of the first interior cell in units of
      ! n0, refreshed by the composition solve (get_species_densities, the
      ! single policy point). The lower boundary needs it to turn cell 1's
      ! pressure into a temperature and a particle count per unit mass,
      ! T(1) = p(1)/n_part(1) and nhat(1) = n_part(1)/rho(1), without
      ! duplicating the electron/nuclei bookkeeping (base_boundary.f90).
      ! Initialized in input_read to the base value ntot_bc + dp_bc.
      real*8  :: n_part_cell1 = 1.0d0

      ! Base boundary-condition anchor ("Base BC: density" / "pressure [<p_ubar>]").
      !   density  (0, default, legacy): the base number density is fixed at n0
      !            ( = 10^"Log10 lower boundary number density" ), and the
      !            reservoir of base_boundary is the isentrope through
      !            (p = ntot_bc + dp_bc, T = T0) at that density.
      !   pressure (1, CETIMB-style): the base is anchored by PRESSURE instead.
      !            n0 is derived so the base pressure n0*kb*T0*ntot_bc equals
      !            base_p_ubar [microbar] (1 microbar = 1 erg/cm^3). The
      !            reservoir pressure is then base_p_ubar and its temperature
      !            T0, at the base LEVEL r = 1 -- not at the first face, which
      !            base_boundary reaches by continuing the reservoir isentrope,
      !            so the level a user states does not move with the grid.
      !            Its physical effect is a much less dense base (1 microbar
      !            ~ 20x below n0 = 1e14), shrinking the dense-base gravity
      !            source rho*g that drives the breathing.
      integer :: base_bc_mode = 0
      real*8  :: base_p_ubar   = 1.0d0   ! target base pressure [microbar] (mode 1)

      ! Periodic Shapiro (1970) low-pass filter ("Shapiro filter: <eps> [<every>]").
      ! Damps the gravity-unbalanced sound-wave (base-breathing) instability the
      ! same way CETIMB (Koskinen et al. 2013a) does: every shapiro_every steps a
      ! weak 1-2-1 filter u_j += (eps/4)(u_{j-1}-2u_j+u_{j+1}) is applied to the
      ! conservative variables. shapiro_eps in [0,1]; eps <= 0 disables it.
      ! OFF BY DEFAULT (opt-in). It damps the HD189733b breathing TRANSIENT, but
      ! the HD209458b cold-IC S x V sweep showed that with the filter ON the
      ! clean cold-IC solution is driven OFF the transonic-wind saddle into
      ! the INFALL attractor (v < 0
      ! everywhere) regardless of the velocity BC, while filter-OFF relaxes to a
      ! clean outflow. So it must NOT be a global default; enable for each run with
      ! "Shapiro filter: <eps> [<every>]" only for cases that actually breathe.
      real*8  :: shapiro_eps   = -1.0d0
      integer :: shapiro_every = 4

      ! Low-Mach contact-mode dissipation ("Low-Mach damping: <eps4> [<M_th>]").
      ! The HLLC contact wave is carried at S* ~ v, so the dissipation the flux
      ! applies to the entropy/contact family vanishes as the flow stagnates and
      ! a 2 dr mode in v and T becomes marginally damped. A gated fourth-
      ! difference (Jameson-Schmidt-Turkel) flux restores that damping where
      ! M < lowmach_damp_mach_th and vanishes identically above it. Added to the
      ! numerical flux inside RK_rhs, so the marching loop and the steady
      ! (Newton) residual see the SAME equation -- unlike the Shapiro filter,
      ! which touches only the marching state. Full statement, stability bound
      ! eps4 < 1/(16 CFL), and measured magnitudes:
      ! src/modules/flux/low_mach_dissipation.f90.
      ! OFF BY DEFAULT: eps4 <= 0 disables, and the code path is then skipped
      ! entirely, so a run without the key is byte-identical to the code
      ! without the term.
      real*8  :: lowmach_damp_eps     = -1.0d0
      real*8  :: lowmach_damp_mach_th =  1.0d-3

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
      ! has its cellwise maximum over [j_min:N] of each row divided by that
      ! row's own largest terms (Update log sections 143, 145: mass by
      ! max(|F| r^2 at the two faces)/dV, momentum by max(|dF_2|,|S_2|), energy
      ! by max(|dF_3|,|S_3|,heat,cool)) below resid_th, REPLACING
      ! the du<du_th test (du measures only the mass-flux spread and is blind
      ! to an operator-split energy imbalance: at the premature WASP golden,
      ! du/dtu were tiny while the energy residual was ~30). Evaluated every
      ! N_resid steps (cheap: one Reconstruct + RK_rhs, no extra ioniz_eq).
      ! Default <= 0 disables (legacy du-based stop, byte-identical).
      real*8  :: resid_th = -1.0d0
      integer :: N_resid  = 500

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
      ! Flux gate (docs/Update_EXHALE_stage1.pdf section 133). The steady solve is
      ! accepted only when BOTH the residual ||R|| < resid_tol AND the radial
      ! spread of the Riemann FACE mass flux over r >= r_flux is below
      ! flux_spread_th (section 145). The two
      ! say different things: the residual that nothing is changing, the flux
      ! that the wind carries one mass flux at every altitude.
      !   flux_spread_th  "Flux spread tol: <tol> [<r_flux>]"; <=0 disables.
      !   r_flux          inner edge of the window, in R_p.
      ! Defaults are measured (section 133): every converged state in hand is
      ! flat to 2.4e-3 - 4.8e-4 over r >= 1.2 R_p, while the state section
      ! 126.6 showed the old measure could not reject reads 7.6e-2 there.
      ! RE-DERIVED for the face-flux definition of section 145. The old
      ! 5.0e-3 was calibrated against the CELL-CENTRED product; on the
      ! conserved face flux it is not a gate at all. Measured face spread at
      ! r >= r_flux over every state available on 2026-09-03:
      !   accepted by the JFNK   1.0e-10 (the section 143 hot-Uranus root),
      !                          3.9e-08, 4.0e-06, 3.8e-06
      !   stopped on du only     1.1e-04, 1.0e-03, 3.2e-02, 3.2e-02
      !   marching, 600 steps    4.0 to 1.4e+01
      ! The two populations are separated by the interval 4e-6 to 1.1e-4;
      ! 2e-5 is its geometric middle to one digit -- five times above the
      ! loosest accepted state and six times below the tightest du stop --
      ! and is not fitted to any one of them.
      real*8, parameter :: flux_spread_th_default = 2.0d-5
      real*8  ::  flux_spread_th = flux_spread_th_default
      real*8  ::  r_flux         = 1.2d0
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
      ! True for an element whose reservoir the lower-atmosphere handoff
      ! states itself -- a "<El>_H_base" key of base.inp, or an elemental
      ! ratio of the file named by "Lower atmosphere profile:". Set in
      ! set_element_abundance, the single door both handoffs go through, and
      ! false for an abundance that comes from metals.inp alone. load_IC uses
      ! it to decide whether a restart column may be renormalized onto the
      ! reservoir the handoff states (see load_IC).
      logical, allocatable :: melem_from_handoff(:)
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
      ! The loaded SED as the FILE states it, ascending in photon energy:
      ! e_sed_node the tabulated energies [eV] and F_sed_node the flux per
      ! unit photon energy [erg cm^-2 s^-1 eV^-1] at the planet, before the
      ! dayside dilution. This is the FIELD; e_v/de_v above is the
      ! quadrature partition the field is integrated on, whose points are
      ! bin centres and not table rows. Filled by read_sed, unallocated for
      ! every other spectrum type.
      real*8, dimension(:), allocatable :: e_sed_node, F_sed_node
      real*8, dimension(:), allocatable :: s_hi,s_hei,s_heii,s_heiTR
      ! H2 photoionization cross section on the energy grid (molecular;
      ! Yan+1998 fit, filled in set_energy_vectors)
      real*8, dimension(:), allocatable :: s_h2
      ! The DISSOCIATIVE part of that cross section, s_h2 times the
      ! branching of frac_H2_dissociative_ionization: the sub-channel that
      ! leaves H + H+ instead of H2+. It is a share of s_h2, not an
      ! addition to it, so the H2 opacity is unchanged.
      real*8, dimension(:), allocatable :: s_h2_di
      ! The DOUBLE-ionization part, H2 + hv -> H+ + H+ + 2e- (threshold
      ! 51.4 eV): the sub-channel that leaves TWO protons and no bound
      ! fragment. Like s_h2_di it is a share of s_h2, not an addition to
      ! it, so the H2 opacity is unchanged; it is zero unless
      ! h2_double_ionization selects a model for it.
      real*8, dimension(:), allocatable :: s_h2_dd
      ! The NEUTRAL-dissociation part, H2 + hv -> H + H: the share of the
      ! absorptions that make no ion at all, i.e. the sub-unity
      ! photoionization yield measured over 33-41 eV. Again a share of
      ! s_h2 and not an addition to it, so the H2 opacity is unchanged;
      ! it is zero unless h2_neutral_dissociation is set.
      real*8, dimension(:), allocatable :: s_h2_nd
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
      real*8, dimension(:), allocatable :: opa_pf
      real*8, dimension(:), allocatable :: r,r_edg,dr_j
      real*8, dimension(:), allocatable :: Gphi_c,Gphi_i

      !------- excited hydrogen H(n=2) coupling -------!
      ! Christie+2013 / Huang+2017 n=2 (2s/2p) model feeding back into
      ! the H ionization balance and the energy equation. All terms are
      ! ZERO unless use_excited_H = .true. (enabled at runtime when a
      ! stellar T_eff is supplied in input.inp), so the default build is
      ! byte-identical to the result with excited-H off. See excited_hydrogen.f90.
      logical :: use_excited_H   = .false. ! master switch (set if T_star_eff>0)
      ! Collisional de-excitation of H(n=2) returns 10.2 eV to the electron
      ! gas (Hdx_arr). It is NOT a double count of the H I collisional-
      ! excitation cooling: the Cen (1992) coefficient in Cool_coeff.f90
      ! (coex_rate_HI) is a one-way, Boltzmann-suppressed excitation rate --
      ! the coronal limit, in which every excitation is assumed to escape --
      ! and carries no density-dependent de-excitation term. Subtracting the
      ! de-excitation is the correction to that limit. Most of the n=2
      ! population is maintained by Ly-alpha pumping rather than by a
      ! collision, so the term is best read as absorbed Ly-alpha thermalized
      ! by a collision. It is the same size as the photoelectric heating the
      ! code already applies unconditionally, so it is on by default; set
      ! "Deexc heat: False" in input.inp for the one-way coronal ledger.
      logical :: incl_deexc_heat = .true.  ! collisional de-excitation heating
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
      !  jlya_mode = 2 -> in-line escape-probability RT (the static
      !                    plane-parallel damping-wing slab solution of
      !                    Neufeld 1990 eq. 3.27 / Harrington 1973 eq. 40),
      !                    computed every timestep from the current state
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
      ! scatter and build the mean intensity above free-streaming by
      ! E = 1 + (boost - 1)(1 - beta) T_star (-> 1 in the thin outer wind,
      ! -> lya_star_boost where the beam is fully trapped; lya_rt.f90, the
      ! stellar-beam block). Tuned to Huang+2023 Fig. 11.
      ! Editable: "Lya stellar boost [-]:".
      real*8  :: lya_star_boost = 5.0d0
      ! Absorbing (pure-sink) lower boundary for the Ly-alpha field. The local
      ! closure Jbar = S(1-beta) has no bottom boundary: a photon travelling
      ! downward is assumed to come back. Huang et al. (2017) instead terminate
      ! their Monte Carlo domain with a purely absorbing bottom, because at
      ! N(H2) ~ 1e14 cm^-2 the accidental resonances between Ly-alpha and the
      ! H2 Lyman/Werner bands give true (non-scattering) absorption, so the
      ! molecular layer under the wind is a photon sink rather than a mirror.
      ! With this flag jlya_escape_prob puts the planet-ward face of the
      ! trapping slab at the bottom of the domain instead of mirroring the
      ! star-ward one, which shortens the random walk and lowers Jbar (and
      ! the pumped n=2 density) in the few scale heights above the base.
      ! OFF by default => the reflecting-bottom closure.
      ! Editable: "Lya absorbing bottom:".
      logical :: lya_bottom_absorber = .false.
      ! Cell-by-cell feedback arrays injected into ioniz_eq (zero unless enabled):
      real*8, dimension(:), allocatable :: gph_balmer_HI ! extra HI photoion [s^-1]
      real*8, dimension(:), allocatable :: heat_balmer   ! extra heat [erg cm^-3 s^-1]
      ! Cell-by-cell diagnostics for output/Excited_H.txt:
      real*8, dimension(:), allocatable :: Jlya_arr    ! Lya mean intensity [cgs]
      real*8, dimension(:), allocatable :: n2s_arr     ! H(2s) density [cm^-3]
      real*8, dimension(:), allocatable :: n2p_arr     ! H(2p) density [cm^-3]
      real*8, dimension(:), allocatable :: Sproton_arr ! Balmer proton src [cm^-3 s^-1]
      real*8, dimension(:), allocatable :: Hpe_arr     ! photoelec. heat [erg cm^-3 s^-1]
      real*8, dimension(:), allocatable :: Hdx_arr     ! de-excit. heat [erg cm^-3 s^-1]
      ! Ground-state H proton-budget rate coefficients, captured by ioniz_eq on
      ! the converged pass so write_excited_H can compare the n=2 photoionization
      ! proton source against the ground-state channels (Huang Figs. 11/27):
      real*8, dimension(:), allocatable :: gph_ground_HI ! ground-state HI photoion. [s^-1]
      real*8, dimension(:), allocatable :: cion_HI       ! HI collisional ioniz. [cm^3 s^-1]
      real*8, dimension(:), allocatable :: arec_HII      ! HII recombination [cm^3 s^-1]

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

      ! Dayside dilution of the stellar beams -- the ATES "2D approximate
      ! method", which is the run's one statement about how a 1-D radial
      ! model stands in for an irradiated sphere.  Rate/4 -> 1/4,
      ! Rate/2 -> 1/2, anything else -> 1.
      !
      ! ONE DEFINITION.  Every band the star supplies is diluted through this
      ! function and nowhere else: the XUV grid (set_energy_vectors), the
      ! stellar Ly-alpha beam (lya_rt, excited_hydrogen) and the four FUV
      ! bands (fuv_band_flux).  Before section 150 the XUV used an exact
      ! string match while the Ly-alpha beam used a substring test and the
      ! FUV bands used neither, so one run could carry three conventions at
      ! once; the substring form is kept because it is the one that survives
      ! a keyword the input normalization does not rewrite.
      !
      ! THE STATED FLUXES KEEP THEIR MEANING.  F_LW_star, F_FUV_B3/B4 and
      ! F_Lya_star are the band flux AT THE PLANET'S ORBIT, whether stated by
      ! a key or integrated from the spectrum file; the dilution is a
      ! run-wide convention applied where the beam is USED, so the setup
      ! report and the resolved dump still state the flux the run was given.
      !
      ! RANGE.  1/2 is the illuminated hemisphere seen by an optically thin
      ! absorber.  Where the band is shielded the true shell average is
      ! lower, because the slant columns away from the substellar point are
      ! longer: 0.20-0.33 for the Lyman-Werner band over a hot-Uranus
      ! molecular layer.  THAT REFINEMENT
      ! IS NOT IMPLEMENTED: there is no key that asks for it, and this
      ! function is the only dilution in the code, so a shielded band is
      ! diluted by the optically thin factor.
      double precision function dayside_dilution() result(xi)
      if      (index(appx_mth,'Rate/4') .gt. 0) then
         xi = 0.25d0
      else if (index(appx_mth,'Rate/2') .gt. 0) then
         xi = 0.5d0
      else
         xi = 1.0d0
      endif
      end function dayside_dilution

      function reconstruction_operator_label() result(lbl)
      ! THE NAME OF THE DISCRETE OPERATOR THE RUN IS CURRENTLY SOLVING WITH.
      !
      ! rec_method alone is not that name. While the PLM -> WENO3 continuation
      ! is armed the run marches R_lambda, a combination of both operators,
      ! and rec_method is left on 'PLM' -- so a file header built from
      ! rec_method would say PLM of a state produced by neither scheme. Every
      ! writer that records the discretization calls this instead, and it is
      ! defined here once because there is one such name.
      character(len=32) :: lbl
      if (recon_lambda_on) then
         write(lbl,'(A,F6.4,A)') 'PLM+WENO3(lambda=', recon_lambda, ')'
      else
         lbl = rec_method
      endif
      end function reconstruction_operator_label

      logical function assembled_reconstruction_is_plm() result(is_plm)
      ! WHICH RECONSTRUCTION THE ROW OF ONE EVALUATION IS BUILT WITH.
      !
      ! The flag pair is not that answer while the continuation is armed: its
      ! endpoints are PLM at lambda <= 0 and WENO3 at lambda >= 1 whatever
      ! use_plm holds, because every assembly selects on lambda there
      ! (reconstruction_continuation_rhs and the kind-generic text of
      ! hydrodynamic_rows).  An intermediate lambda is neither scheme, and the
      ! callers that cannot carry the blend refuse it before they ask.
      !
      ! It is defined here, beside reconstruction_operator_label, because the
      ! assembly of a row and the terms that row is read against must state
      ! this once: the kind-generic assembly selected the endpoint while
      ! equilibrium_pressure_force_of_state and momentum_row_terms_of_state
      ! read the flags, and the two disagreed silently when the flags named
      ! the other endpoint.
      if (recon_lambda_on) then
         is_plm = (recon_lambda .le. 0.0d0)
      else
         is_plm = use_plm
      endif
      end function assembled_reconstruction_is_plm

      subroutine allocate_grid_arrays
      ! Allocate the grid-sized module arrays once N is known (called from
      ! input_read, right after the "Grid cells" key has been resolved and
      ! before anything builds the grid). The lower bound 1-Ng and the upper
      ! bound N+Ng are the ones the declarations used to carry. The values
      ! assigned here are the initializers those declarations carried; the
      ! arrays that carried none (r, r_edg, dr_j, Gphi_c, Gphi_i) were in
      ! static storage and therefore started at zero, which is reproduced.

      allocate(opa_pf(1-Ng:N+Ng))
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      allocate(Gphi_c(1-Ng:N+Ng), Gphi_i(1-Ng:N+Ng))
      allocate(gph_balmer_HI(1-Ng:N+Ng), heat_balmer(1-Ng:N+Ng))
      allocate(Jlya_arr(1-Ng:N+Ng), n2s_arr(1-Ng:N+Ng), n2p_arr(1-Ng:N+Ng))
      allocate(Sproton_arr(1-Ng:N+Ng))
      allocate(Hpe_arr(1-Ng:N+Ng), Hdx_arr(1-Ng:N+Ng))
      allocate(gph_ground_HI(1-Ng:N+Ng), cion_HI(1-Ng:N+Ng))
      allocate(arec_HII(1-Ng:N+Ng))
      allocate(kzz_cell(1-Ng:N+Ng))

      opa_pf        = 1.0d0
      ! Placeholder until eddy_diffusion_on_grid runs (after the grid exists
      ! and after any base.inp / profile override of he_kzz).
      kzz_cell      = 0.0d0
      r             = 0.0d0
      r_edg         = 0.0d0
      dr_j          = 0.0d0
      Gphi_c        = 0.0d0
      Gphi_i        = 0.0d0
      gph_balmer_HI = 0.0d0
      heat_balmer   = 0.0d0
      Jlya_arr      = 0.0d0
      n2s_arr       = 0.0d0
      n2p_arr       = 0.0d0
      Sproton_arr   = 0.0d0
      Hpe_arr       = 0.0d0
      Hdx_arr       = 0.0d0
      gph_ground_HI = 0.0d0
      cion_HI       = 0.0d0
      arec_HII      = 0.0d0

      end subroutine allocate_grid_arrays

      ! End of module
      end module global_parameters

   module  Read_input
   ! Read input planetary parameters adn define

   use global_parameters
   use metals_input        ! optional metals.inp abundance reader
   use charge_exchange, only: cx_init,       &  ! build active charge-exchange set
                              he_h_charge_exchange  ! He <-> H pair (group B) switch
   use species_table, only: n_melem, iel_C, iel_O, iel_N, iel_Mg,  &
                            iel_Si, iel_Ca, iel_Na, iel_K, iel_S,  &
                            iel_Fe, mion_ethr, melem_i0, melem_A,  &
                            melem_name
   use composition, only: comp_mass_per_H, comp_ntot_bc, comp_rho_bc,   &
                          h2_mixing_ratio_base, h2_mixing_ratio_ceiling, &
                          base_h2_nuclei_fraction, h2_bound_fraction
   use lower_atmosphere_profile, only: lap_file, lap_in_use,             &
                          lap_solution_id, lap_p_match_bar,              &
                          lap_source_code, lap_iteration,                &
                          read_lower_atmosphere_profile,                 &
                          lap_value_at_match, lap_element_ratio_at_match
   ! FUV photolysis thresholds of the oxygen chemistry, filled once here
   ! (serially) because the OpenMP cell sweep only reads them.
   use water_photolysis, only: water_photolysis_init, ib_LW, ib_B3, ib_B4, n_fuv_band
   use sed_reader, only: sed_band_fluxes
   use constrained_chemical_equilibrium, only: cce_dump_paths_read
   use newton_solver, only: newton_solver_read_environment
   use oxygen_rates, only: fuv_band_lo_A, fuv_band_hi_A
   use diffusive_photochemistry, only: carrier_set_init
   use base_boundary, only: set_base_reservoir
   use Numerical_Fluxes, only: low_mach_velocity_jump
   ! The vocabulary of the restart metadata block's 'options' field, defined
   ! once where that field is written, so a token this file accepts on a
   ! "Restart option change:" line is a token the comparison knows.
   use IC_load, only: n_opt, opt_name, opt_changes_layout,               &
                      restart_option_change_named,                      &
                      restart_option_change_given

   implicit none

   ! ---- who fixes the base level ------------------------------------- !
   ! The lower boundary is one level, so exactly one input states it. A
   ! lower-atmosphere handoff -- a base.inp carrying "p_base", or a
   ! lower-atmosphere profile, whose matching pressure p_match is the level
   ! its T, r, q_H2 and elemental ratios are read at -- is written AT that
   ! pressure, so when one is present that pressure IS the base level and n0
   ! follows from it; otherwise the legacy density key states it. These two
   ! record which of them spoke, so the pair can be checked against each
   ! other instead of one silently winning.
   ! WHERE THE COUPLED BLOCK IS ENTERED FROM. "Coupled carrier solve" is
   ! three-valued: False, the alternation of a hydrodynamic solve with the
   ! fixed-wind transport relaxations; True, the block from the first pass;
   ! "On stall", the alternation first and the block from the pass at which
   ! the alternation stops approaching a joint fixed point. The third value
   ! sets this and leaves carrier_in_newton false, because the route of the
   ! first pass is the alternation.
   logical, save, public :: carrier_newton_on_stall = .false.
   logical, save, public :: base_level_from_handoff  = .false.
   logical, save, public :: base_density_key_given   = .false.
   ! Which input stated the level, named once here so the setup report and
   ! the refusal messages below cannot disagree about it.
   character(len=64), save, public :: base_level_source = 'the density key'

   ! ---- what a restart continues ------------------------------------- !
   ! WHAT A RESTART IS FOR, stated once and checked against the other keys.
   ! A loaded state can
   ! be continued in three different senses, and they are not variants of one
   ! path: a physical trajectory carries the clock the file records and every
   ! step obeys the physical-mode contract; a relaxation carries no clock and
   ! marches toward stationarity; a stationary evaluation takes no step at all
   ! and measures the state as it stands. One selector says which, because a
   ! run that does not state its purpose cannot be checked against it.
   integer, parameter, public :: restart_intent_relaxation = 0
   integer, parameter, public :: restart_intent_trajectory = 1
   integer, parameter, public :: restart_intent_stationary = 2
   integer, save, public :: restart_intent = restart_intent_relaxation
   logical, save, public :: restart_intent_given = .false.
   ! The second word of a stationary intent. "evaluate" returns the state
   ! unchanged with its certification and takes no solver step; "equilibrate"
   ! asks for the loaded composition to be put on its own fixed point BEFORE
   ! the evaluation, which is a change of the state and is reported when it
   ! is taken. Neither is the default: the default stationary intent measures
   ! the state as loaded and then enters the stationary solve.
   logical, save, public :: stationary_evaluate_only      = .false.
   logical, save, public :: stationary_equilibrate_loaded = .false.
   ! ---- what a restart is allowed to change -------------------------- !
   ! WHICH PHYSICS OPTIONS A RESTART MAY CHANGE. A state is a
   ! solution of an equation set, so the restart contract refuses a state
   ! whose option set is not the run's; the way this project reaches its
   ! solutions is to converge without an option and restart with it on, so
   ! "Restart option change: <token>[, <token> ...]" names the tokens that
   ! are allowed to differ. Everything else still refuses by name, and the
   ! change is written into the state the run produces. The tokens are the
   ! ones of the '# options' line of a state file (IC_load); the flags they
   ! set live there too, because that is where the comparison is made.

   contains
      
   subroutine input_read
   real*8,  dimension(n_fuv_band) :: band_flux_sed
   logical, dimension(n_fuv_band) :: band_covered_sed
   ! Subroutine to read the input file and assign names and values
   !    to global constants

      character(len = :), allocatable :: str
      ! The second word of a two-word value ("On stall").
      character(len = :), allocatable :: str2
      character(len = 250)            :: line
      character(len = 250), allocatable :: filelines(:)
      integer                         :: ios
      integer                         :: im
      integer                         :: i, ios_pe, nlines
      integer                         :: kk
      logical                         :: is_known
      ! Resolved base H2 mixing ratio and the ceiling it is checked against
      real*8                          :: q_h2_resolved, q_h2_max
      real*8                          :: p_base_cgs, n0_from_p_base
      real*8                          :: p_from_density_key
      ! The forced base reservoir particle count and the value the base
      ! composition gives, kept apart so the log carries both.
      character(len = 64)             :: ntot_bc_forced_env
      real*8                          :: ntot_bc_of_composition
      real*8                          :: ntot_bc_forced
      integer                         :: ios_ntot_bc
      integer                         :: iw
      ! The token list of "Restart option change" and its two counters
      character(len = 250)            :: optline
      integer                         :: iopt, ktok

   ! Every label input_read recognizes: the core block followed by the
   ! keyword-extension block, in the same order as the reads below. It has two
   ! readers, and one of them stops the run:
   !   - refuse_duplicate_keys, called on the loaded file before anything is
   !     consumed, which REFUSES a file stating any of these keys twice. A key
   !     absent from this list is therefore outside that refusal as well as
   !     outside the warning below, so a key parsed by the keyword loop and
   !     missing here would silently accept two answers to one question.
   !   - the trailing unknown-line scan, which WARNS (never stops) on a
   !     non-blank, non-'#' line matching no key -- e.g. a "Newton Solver:"
   !     (capital S) typo that anchored matching would otherwise ignore.
   ! The energy-band line ("[E_low...") is a bracket prefix, checked
   ! separately there. base.inp has its own list, built in read_base_inp.
   character(len=32), parameter :: known_keys(*) = [ character(len=32) ::    &
      'Planet name', 'Log10 lower boundary', 'Planet radius', 'Planet mass', &
      'Equilibrium temperature', 'Orbital distance', 'Escape radius',        &
      'He/H number ratio', '2D approximate method', 'Parent star mass',      &
      'Spectrum type', 'Spectrum file', 'Power-law index', 'Photon energy',  &
      'Use only EUV', 'Log10 of X-ray luminosity', 'Log10 of EUV luminosity',&
      'Grid type', 'Base grid', 'Grid cells',                                &
      'Reconstruction continuation',                                         &
      'Numerical flux', 'Reconstruction scheme', 'Include He23S',            &
      'Load IC', 'Do only PP', 'Force start',                                &
      'Domain mode', 'Outer radius', 'Stellar Teff', 'Stellar radius',       &
      'Caloric EOS', 'Photoelectron heating',                               &
      'Deexc heat', 'Wind-AE seed out', 'Wind-AE seed', 'Jlya RT file',      &
      'Jlya escape-prob', 'Stellar Lya flux', 'Lya stellar halfwidth',       &
      'Lya stellar boost', 'Lya absorbing bottom',                           &
      'du_th', 'ATES_photoionization_rate',                                  &
      'Legacy_HHe_rates', 'Secondary_ionization', 'He_rec_coupling',         &
      'He_H_charge_exchange', 'Atomic rate set',                             &
      'Molecular chemistry', 'Molecular base', 'Stellar LW flux',           &
      'Oxygen chemistry', 'Molecular carrier transport',                   &
      'Ionization transport',                                                  &
      'Coupled carrier solve',                                             &
      'Oxygen transport',                                                  &
      'Stellar FUV B1 flux', 'Stellar FUV B3 flux',                        &
      'Stellar FUV B4 flux',                                               &
      'Lower atmosphere', 'Lower atmosphere profile',                        &
      'Lower column', 'He_Kzz', 'He_alphaT', 'He_ambipolar',                 &
      'He_metal_diffusion', 'He_diffusion', 'Stall', 'Energy solver',        &
      'Time stepping', 'Level tol', 'Solver', 'Valve eps', 'Hydrostatic base',&
      'Shapiro filter', 'Low-Mach damping', 'Well balanced',                 &
      'Low Mach velocity jump',                                              &
      'Base BC', 'Base velocity', 'Viscosity',                               &
      'Base ghost temperature', 'Max steps', 'Coronal cutoff width',         &
      'Base IR field', 'Molecular IR bands', 'Molecular reaction heat',       &
      'H2 double ionization', 'H2 neutral dissociation',                     &
      'Conduction', 'Resid tol',                                             &
      'Resid norm', 'Flux spread tol', 'CFL', 'Transonic IC',            &
      'Hot Parker IC', 'IC mode',                                            &
      'Newton solver', 'Brent solver', 'Run mode', 'Restart intent',       &
      'Restart option change' ]

   ! ----- Read planetary parameters from input file ----- !

   ! Metal abundances (default: no metals). Set at runtime by an optional
   ! metals.inp file (no recompile) and, for the elements a lower-atmosphere
   ! handoff carries, overridden afterwards by the "<El>_H_base" keys of
   ! base.inp. thereis_metals and melem_ab are derived from the final values
   ! in the composition block near the end of this subroutine.
   X_C  = 0.0d0
   X_N  = 0.0d0
   X_O  = 0.0d0
   X_Mg = 0.0d0
   X_Si = 0.0d0
   X_Ca = 0.0d0
   X_Na = 0.0d0
   X_K  = 0.0d0
   X_S  = 0.0d0
   X_Fe = 0.0d0
   call read_metals_input

   ! ----- Read every line of input.inp into memory -----
   ! The core block below is matched by LABEL: each key is anchored at the
   ! start of the left-trimmed line and terminated by its value separator
   ! (':', '?', whitespace, or '='), so the physical line order no longer
   ! matters and blank / '#'-comment lines are ignored. Legacy positional
   ! files parse unchanged because their lines are self-labeling (e.g.
   ! "Planet radius [R_J]: 1.401"); within a line the value is still taken by
   ! word position with get_word, exactly as before. The optional keyword-
   ! extension block that follows uses the same anchored matching. A KEY MAY
   ! APPEAR ONLY ONCE: two lines stating the same quantity are two answers to
   ! one question, and silently keeping the last of them is how a run acquires
   ! a tolerance nobody chose (measured: valve_sens/eps5 carried "Resid tol"
   ! twice and ran on the second). refuse_duplicate_keys below stops the run
   ! and names both lines.
   write(*,*) '(input_read.f90) Reading the input.inp file..'
   open(unit = 11, file = inp_file)
   nlines = 0
   do
      read(11,'(A)',iostat = ios) line
      if (ios .ne. 0) exit
      nlines = nlines + 1
   enddo
   allocate(character(len=250) :: filelines(nlines))
   rewind(11)
   do i = 1, nlines
      read(11,'(A)') filelines(i)
   enddo
   close(11)

   ! ----- One line per key -----
   ! Checked BEFORE anything is consumed, so a file that states a quantity
   ! twice never gets as far as producing a run with one of the two values.
   call refuse_duplicate_keys(inp_file, filelines, nlines,                &
                              known_keys, size(known_keys))

   ! ----- Core block (label-matched; order-independent) -----

   ! Planet name (single token at word 3)
   p_name = get_word(req('Planet name'), 3)

   ! Log10 of n0. Optional: a base.inp that states p_base fixes the base
   ! level instead, and then this key is redundant (and is refused if it
   ! disagrees -- see "the base level has one source" below). n0 is left
   ! non-positive here so that a run stating neither is caught there.
   line = find_lbl('Log10 lower boundary', base_density_key_given)
   if (base_density_key_given) then
      str = get_word(line, 7);  read(str,*) n0
   else
      n0 = -1.0d0
   endif

   ! Planet radius
   str = get_word(req('Planet radius'), 4);  read(str,*) R0

   ! Planet mass
   str = get_word(req('Planet mass'), 4);  read(str,*) Mp

   ! Equilibrium temperature
   str = get_word(req('Equilibrium temperature'), 4);  read(str,*) T0

   ! Orbital distance
   str = get_word(req('Orbital distance'), 4);  read(str,*) a_orb

   ! Escape radius
   str = get_word(req('Escape radius'), 4);  read(str,*) r_esc

   ! He/H number ratio
   str = get_word(req('He/H number ratio'), 4);  read(str,*) HeH
   if (HeH .gt. 0.0e0) thereis_He = .true.

   ! 2D approximate method
   line = req('2D approximate method')
   appx_mth = get_word(line, 4)
   ! Read alpha if selected
   if (appx_mth .eq. 'alpha') then
      str = get_word(line, 6)
      read(str,*) a_tau
   else
      a_tau = 0.0
   endif
   ! Correct appx_meth keywords
   if (appx_mth .eq. 'Rate/4') appx_mth = 'Rate/4 + Mdot'
   if (appx_mth .eq. 'Rate/2') appx_mth = 'Rate/2 + Mdot/2'

   ! Parent star mass
   str = get_word(req('Parent star mass'), 5);  read(str,*) Mstar

   ! Spectrum type selects WHICH property line is consumed (Spectrum file /
   ! Power-law index / Photon energy); the value word position within that
   ! line is unchanged from the legacy layout.
   sp_type = get_word(req('Spectrum type'), 3)
   select case (sp_type)

      case ('Load')          ! Load from file
         sed_file    = get_word(req('Spectrum file'), 3)
         do_read_sed = .true.

      case ('Power-law')
         str = get_word(req('Power-law index'), 3)
         read(str,*) PLind
         is_PL_sed = .true.

      case ('Planck')
         ! The whole photon grid is the photospheric blackbody
         ! pi B_nu(T_eff) (R_star/a)^2, built from "Stellar Teff [K]:" and
         ! "Stellar radius [R_sun]:" (development plan rev 3, section 10.5
         ! decision 13: one spectrum type builds every band). Those two
         ! lines belong to the optional keyword block below, which has not
         ! been scanned yet; that they are present is checked once it has
         ! been ("Spectrum type: Planck needs both stellar quantities").
         ! No property line of its own, so nothing is consumed here.

      case ('Monochromatic')
         is_monochr = .true.
         ! Read photon energy
         str = get_word(req('Photon energy'), 4)
         read(str,*) e_low
         ! Remove helium if monochromatic and photon energy lower than the
         ! helium ionization threshold. Also zero the He/H ratio so the whole
         ! downstream (mass_per_H, rho_bc, set_IC He fractions, EOS) is a
         ! self-consistent H-only gas rather than carrying He mass with
         ! thereis_He = .false.
         if (e_low .lt. e_th_HeI) then
            thereis_He = .false.
            HeH        = 0.0d0
            write(*,*) '(input_read) Monochromatic photon energy below'//&
               ' the He I ionization threshold: gas treated as pure'//&
               ' hydrogen (He/H set to 0).'
         endif

      case default
         write(*,*) '(input_read) ERROR: unknown spectrum type "'//&
            trim(sp_type)//'".'
         write(*,*) '   Allowed values: Load, Power-law, Planck,'//&
            ' Monochromatic.'
         error stop 1

   end select

   ! Only EUV status: word 4 == 'False' includes the X-rays (see schema).
   str = get_word(req('Use only EUV'), 4)
   if (str .eq. 'False') thereis_Xray = .true.

   ! Energy bands, present unless the spectrum is monochromatic. Word
   ! positions 4/6/8 land on the numbers regardless of the value separator
   ! used inside the brackets.
   if (.not. is_monochr) then
      line = req_eband()
      str = get_word(line, 4);  read(str,*) e_low
      str = get_word(line, 6);  read(str,*) e_mid
      if (thereis_Xray) then
         str = get_word(line, 8);  read(str,*) e_top
      else
         e_top = 1.24e3
      endif
   endif

   ! X-ray luminosity, only when X-rays are included
   if (thereis_Xray) then
      str = get_word(req('Log10 of X-ray luminosity'), 6)
      read(str,*) LX
   else
      LX = 0.0
   endif

   ! EUV luminosity
   str = get_word(req('Log10 of EUV luminosity'), 6);  read(str,*) LEUV

   ! Grid type
   grid_type = get_word(req('Grid type'), 3)

   ! Numerical flux
   flux = get_word(req('Numerical flux'), 3)

   ! Reconstruction scheme. "PLM" and "WENO3" are single-stage and use only
   ! the FIRST du_th value; "PLM+WENO3" is two-stage (PLM then WENO3) and
   ! uses BOTH du_th values (see the du_th parsing below).
   !
   ! rec_method is the scheme Reconstruct() dispatches on, and it is ALWAYS
   ! left as 'PLM' or 'WENO3' -- never the input word 'PLM+WENO3', which
   ! names a two-stage RECIPE and not a reconstruction. Which stage the run
   ! starts in is decided once both du_th values are known, at the end of
   ! this routine; the two-stage intent travels in recon_two_stage.
   str = get_word(req('Reconstruction scheme'), 3)
   if (str .eq. 'WENO3') then
      rec_method = 'WENO3';  use_weno3 = .true.
   else if (str .eq. 'PLM') then
      rec_method = 'PLM';    use_plm   = .true.
   else if (str .eq. 'PLM+WENO3') then
      rec_method      = 'PLM'     ! provisional; resolved at the end
      use_plm         = .true.
      recon_two_stage = .true.
   else
      write(*,*) '(input_read) ERROR: unknown "Reconstruction scheme": ',   &
                 trim(str), '  (expected PLM, WENO3 or PLM+WENO3)'
      error stop 1
   endif

   ! Include He23S -- OPTIONAL.  The triplet defaults to ON (parameters.f90);
   ! an input.inp that omits the line gets it.  "Include He23S? False" is the
   ! deliberate opt-out and is what the HeITR-off branch check uses.
   line = find_lbl('Include He23S', is_known)
   if (is_known) then
      str = get_word(line, 3)
      if (str .eq. 'True'  .or. str .eq. 'true' ) thereis_HeITR = .true.
      if (str .eq. 'False' .or. str .eq. 'false') thereis_HeITR = .false.
   endif

   ! IC status
   str = get_word(req('Load IC'), 3)
   if (str .eq. 'True')  do_load_IC = .true.

   ! Do only post-processing. Kept BEFORE "Force start" so that, when both
   ! are True, force_start wins (matching the legacy sequential order).
   str = get_word(req('Do only PP'), 4)
   if (str .eq. 'True')  then
      do_only_pp  = .true.
      force_start = .false. ! Set to false to avoid overlap
   endif

   ! Force start of sim.
   str = get_word(req('Force start'), 3)
   if (str .eq. 'True')  then
      force_start = .true.
      do_only_pp  = .false. ! Set to false to avoid overlap
   endif

		! ---- Optional domain-extent option ----
		! Appended at the end of input.inp (after "Force start:") so older
		! files lacking these lines keep the default Roche/Hill behavior.
		! "Domain mode: Spherical" + "Outer radius [R_p]: <value>" switches
		! to a pure planetary potential (-b0/r) extended to <value> R_p
		! (Huang Case A-like). Scanned by keyword (index) so blank trailing
		! lines and line ordering do not matter.
		! excited-H option (also appended after "Force start:" so
		! older files are unaffected): "Stellar Teff [K]: <T>" + "Stellar
		! radius [R_sun]: <R>" supply the diluted-blackbody Balmer continuum
		! that photoionizes/heats H(n=2). The coupling is enabled iff both are
		! given (T_star_eff>0, R_star>0). The collisional de-excitation heating
		! of the pumped n=2 population is on with them; "Deexc heat: False"
		! drops back to the one-way coronal ledger (see parameters.f90).
		spherical_domain = .false.
		r_out_user       = 0.0d0
		T_star_eff       = 0.0d0
		R_star           = 0.0d0
		incl_deexc_heat  = .true.
		jlya_mode        = 0
		lya_bottom_absorber = .false.  ! reflecting bottom (legacy closure)
		transonic_ic     = .false.
		hot_parker_ic    = .false.
		T_wind_ic        = 1.0d4
		use_newton_ieq   = .true.    ! upgraded solvers are the default
		use_brent_tsolve = .true.
		windae_seed_file = 'inputdata/windae_seed.csv'
		windae_seed_out  = ''
		dr_base          = dr_base_default  ! uniform base cell size [R_p]
		                                   ! (Mixed grid), see parameters.f90
		dr_base_from_key = .false.    ! set by the "Base grid" key below
		N_low_cells      = N_low_cells_default ! number of uniform base cells
		base_bc_mode     = 0          ! density-anchored base (legacy) by default
		ates_photoion_rate = .false.  ! default: Verner+1996 He I (1^1S) photoion.
		legacy_hhe_rates   = .false.  ! default: Badnell/Mao + Voronov H/He rates
		atomic_rate_set_k22 = .false. ! default: EXHALE's own atomic H/He rates
		use_sec_ion        = .true.   ! default: SvS85 secondary ionization ON
		sec_ion_immediate  = .false.  ! default: staged (applied after 1st converge)
		use_he_rec_coupling = .true.  ! default: He rec. photons ionize/heat H
		he_h_charge_exchange = .true. ! default: He <-> H charge exchange (group B) ON
		do i = 1, nlines
			line = filelines(i)
			if (len_trim(line) .eq. 0) cycle
			if (index(adjustl(line),'#') .eq. 1) cycle
			if (lbl_match(line, 'Domain mode')) then
				str = get_word(line, 3)
				if (str .eq. 'Spherical') spherical_domain = .true.
			else if (lbl_match(line, 'Outer radius')) then
				str = get_word(line, 4)
				read(str,*) r_out_user
			else if (lbl_match(line, 'Stellar Teff')) then
				str = get_word(line, 4)
				read(str,*) T_star_eff
			else if (lbl_match(line, 'Stellar radius')) then
				str = get_word(line, 4)
				read(str,*) R_star
			else if (lbl_match(line, 'Deexc heat')) then
				str = get_word(line, 3)
				if (str .eq. 'True'  .or. str .eq. 'true' ) incl_deexc_heat = .true.
				if (str .eq. 'False' .or. str .eq. 'false') incl_deexc_heat = .false.
			else if (lbl_match(line, 'Wind-AE seed out')) then
				windae_seed_out = trim(get_word(line, 4))
			else if (lbl_match(line, 'Wind-AE seed')) then
				windae_seed_file = trim(get_word(line, 3))
			else if (lbl_match(line, 'Jlya RT file')) then
				jlya_rt_file = get_word(line, 4)
				jlya_mode    = 1
			else if (lbl_match(line, 'Jlya escape-prob')) then
				str = get_word(line, 3)
				if (str .eq. 'True') jlya_mode = 2
			else if (lbl_match(line, 'Stellar Lya flux')) then
				str = get_word(line, 5)
				read(str,*) F_Lya_star
			else if (lbl_match(line, 'Lya stellar halfwidth')) then
				str = get_word(line, 5)
				read(str,*) dv_star_lya
			else if (lbl_match(line, 'Lya stellar boost')) then
				str = get_word(line, 5)
				read(str,*) lya_star_boost
			else if (lbl_match(line, 'Lya absorbing bottom')) then
				! "Lya absorbing bottom: True" terminates the Ly-alpha domain
				! with a pure absorber (the H2 layer below the base; Huang
				! et al. 2017), adding the downward escape as a loss channel.
				str = get_word(line, 4)
				if (str .eq. 'True'  .or. str .eq. 'true' ) lya_bottom_absorber = .true.
				if (str .eq. 'False' .or. str .eq. 'false') lya_bottom_absorber = .false.
			else if (lbl_match(line, 'du_th')) then
				! "du_th [PLM,WENO3]: <du1> [<du2>]". How the numbers are used is
				! decided by "Reconstruction scheme:" (parsed above):
				!   PLM+WENO3 -> two-stage: PLM until du<du1, then WENO3 until du<du2.
				!   PLM or WENO3 -> single-stage at du1; any second number is ignored.
				str = get_word(line, 3);  read(str,*) du_th      ! first threshold
				str = get_word(line, 4)                          ! optional second
				if (recon_two_stage) then
					if (len_trim(str) .gt. 0) then
						du_th_plm = du_th   ! first number = stage-1 (PLM) threshold
						read(str,*) du_th   ! second number = stage-2 (WENO3) threshold
					else
						write(*,*) '(input_read) WARNING: "Reconstruction scheme:'//&
						   ' PLM+WENO3" needs two du_th values; only one given'//&
						   ' -> the PLM stage is empty, running single-stage'//  &
						   ' WENO3 at that threshold.'
						du_th_plm = -1.0d0
					endif
				else
					! PLM or WENO3: single-stage; ignore any second du_th value.
					du_th_plm = -1.0d0
				endif
			else if (lbl_match(line, 'Reconstruction continuation')) then
				! "Reconstruction continuation: <dlambda> [<dtu_tol>] [fixed]"
				! Walk the PLM -> WENO3 hand-off along the homotopy
				! R_lambda = (1-lambda) R_PLM + lambda R_WENO3 instead of changing
				! the discrete operator in one step (see recon_lambda_step0 in
				! parameters.f90 and docs/input_schema.md). Only "Reconstruction
				! scheme: PLM+WENO3" has a hand-off to walk.
				!   dlambda   step in lambda per marching step; <= 0 disables,
				!             1.0 is the one-step switch this replaces
				!   dtu_tol   lambda advances while dtu <= tol * dtu at the start
				!             of the ramp (default 1.2)
				!   fixed     ramp unconditionally, with no step control
				! The two optional fields are told apart by content, not position:
				! 'fixed'/'adaptive' is the mode, anything else is the tolerance.
				str = get_word(line, 3);  read(str,*) recon_lambda_step0
				do iw = 4, 5
					str = get_word(line, iw)
					if (len_trim(str) .eq. 0) cycle
					if (str .eq. 'fixed') then
						recon_lambda_adaptive = .false.
					else if (str .eq. 'adaptive') then
						recon_lambda_adaptive = .true.
					else
						read(str,*) recon_lambda_dtu_tol
					endif
				enddo
				if (recon_lambda_step0 .gt. 0.0d0)                             &
					write(*,'(A,ES9.2,A,ES9.2,A,A)') ' (input_read) PLM -> '//  &
						'WENO3 continuation: dlambda =', recon_lambda_step0,     &
						', dtu tolerance =', recon_lambda_dtu_tol, ', ',         &
						trim(merge('adaptive', 'fixed   ', recon_lambda_adaptive))
			else if (lbl_match(line, 'ATES_photoionization_rate')) then
				! Revert He I (1^1S) photoionization to the legacy ATES 2-term fit
				! (default is Verner+1996). "ATES_photoionization_rate: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') ates_photoion_rate = .true.
			else if (lbl_match(line, 'Legacy_HHe_rates')) then
				! Revert H/He recombination + collisional ionization to the legacy
				! ATES fits (default is Badnell/Mao case B + Voronov 1997).
				! "Legacy_HHe_rates: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') legacy_hhe_rates = .true.
			else if (lbl_match(line, 'Atomic rate set')) then
				! Swap the four atomic H/He rate coefficients Koskinen et al.
				! (2022, ApJ 929, 52) list in their Table 1 as R1-R4 --
				! radiative recombination of H+ and He+, electron-impact
				! ionization of H and He -- for the published fits, so that a
				! run can be compared like for like with their Model A. Only
				! those four move: the recombination COOLING rates, the metal
				! rates and every other coefficient stay as they are. The
				! default set is the physically preferred one here (case B,
				! the Lyman continuum being optically thick), so the key is
				! for reproducing their choice, not for replacing ours.
				! "Atomic rate set: Koskinen2022" (value at word 4).
				str = get_word(line, 4)
				if (str .eq. 'Koskinen2022' .or. str .eq. 'koskinen2022') then
					atomic_rate_set_k22 = .true.
					write(*,'(A)') ' (input_read) Atomic rate set:'//         &
					   ' Koskinen et al. (2022) Table 1 R1-R4 (H+/He+'//      &
					   ' recombination 4.0e-12/4.6e-12 (300/T)^0.64,'//       &
					   ' Voronov 1997 collisional ionization)'
				else if (str .eq. 'default' .or. str .eq. 'Default') then
					atomic_rate_set_k22 = .false.
				else
					write(*,*) '(input_read) ERROR: "Atomic rate set" must'
					write(*,*) '  be Koskinen2022 or default. Got: '//trim(str)
					error stop 1
				endif
			else if (lbl_match(line, 'Secondary_ionization')) then
				! SvS85 secondary ionization (default ON, staged: applied only
				! after the wind first converges without the coupling).
				! "Secondary_ionization: False"     -> off entirely
				! "Secondary_ionization: Immediate"  -> on from step 0
				!    (restores the pre-staging behavior, for A/B tests only)
				str = get_word(line, 2)
				if (str .eq. 'False' .or. str .eq. 'false') use_sec_ion = .false.
				if (str .eq. 'True'  .or. str .eq. 'true' ) use_sec_ion = .true.
				if (str .eq. 'Immediate' .or. str .eq. 'immediate') then
					use_sec_ion       = .true.
					sec_ion_immediate = .true.
				endif
			else if (lbl_match(line, 'He_rec_coupling')) then
				! Couple He II -> He I recombination radiation to H ionization
				! (Draine 2011 y/z, on-the-spot). Default ON (the photons are
				! real; "He_rec_coupling: False" restores the legacy lost-photon
				! path, which in TR mode leaves the singlet recombination
				! neither case A nor case B).
				str = get_word(line, 2)
				if (str .eq. 'True'  .or. str .eq. 'true' ) use_he_rec_coupling = .true.
				if (str .eq. 'False' .or. str .eq. 'false') use_he_rec_coupling = .false.
			else if (lbl_match(line, 'He_H_charge_exchange')) then
				! He <-> H charge exchange (Huang 2023 Table 4 group B, rates
				! from Koskinen 2013): He0+H+ <-> He++H0. Default ON in every
				! ionization system that contains He. "He_H_charge_exchange:
				! False" restores the legacy no-He-CX path (Group B not
				! assembled anywhere).
				str = get_word(line, 2)
				if (str .eq. 'True'  .or. str .eq. 'true' ) he_h_charge_exchange = .true.
				if (str .eq. 'False' .or. str .eq. 'false') he_h_charge_exchange = .false.
			else if (lbl_match(line, 'Molecular chemistry')) then
				! molecular network (docs/lower_atmosphere_*). Solving molecular
				! chemistry implies the molecular-base particle count for the
				! base pressure BC: with an atomic ntot_bc the base temperature
				! is inflated by the particle-count ratio (~1.8x for a fully
				! molecular base), thermally dissociating the H2 layer the
				! chemistry just built. molecular_base is therefore coupled on
				! below (after all keys are parsed).
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true') thereis_mol = .true.
			else if (lbl_match(line, 'Oxygen chemistry')) then
				! In-code oxygen chemistry (the A2 option): OH / H2O / CO
				! solved in the coupled molecular ionization equilibrium,
				! with the H2O and OH photolysis of the FUV bands. It is
				! what lets the code compute its own base H2/H partition
				! instead of importing it.
				! Requires the molecular network, helium and oxygen; the
				! checks are below, after every key is parsed.
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true')                  &
					thereis_oxychem = .true.
			else if (lbl_match(line, 'Molecular carrier transport')) then
				! "Molecular carrier transport: True|False" -- vertical
				! transport of the molecular carriers, solved implicitly
				! with their chemistry (diffusive_photochemistry). The set
				! this key fixes is H2 alone, or H2 with OH, H2O and CO when
				! the oxygen cycle is on; the three ionization stages are
				! carried on "Ionization transport" instead, which is a
				! separate set of rows of the same operator and exists in an
				! atomic gas too (carrier_set_init). False restores the
				! local-kinetics limit of milestone M2, which isolates the
				! chemistry for testing and is not a model of a base. Left
				! unstated, the default is resolved below from the chemistry
				! the run carries. This key means nothing with the molecular
				! network off, because there is no molecular carrier to
				! transport; that is reported below, not refused.
				str = get_word(line, 4)
				carrier_transport_stated = .true.
				carrier_transport = (str .eq. 'True' .or. str .eq. 'true')
			else if (lbl_match(line, 'Ionization transport')) then
				! "Ionization transport: True|False" -- carry the
				! ionization state of hydrogen and helium with the flow
				! instead of re-solving each cell's partition as a local
				! equilibrium every step. The carried set is x(H II) per
				! hydrogen nucleus and x(He II), x(He III) per helium
				! nucleus, each transported on its own element's nucleus
				! face flux by the transport-chemistry operator; the neutral
				! stage of each element closes its simplex and is not a row.
				! It applies to an atomic gas as well as a molecular one,
				! and each stage's source is the ionization balance the
				! local sweep of that gas solves. Default False. The
				! requirements are checked below, after every key is parsed.
				str = get_word(line, 3)
				ionization_transport = (str .eq. 'True' .or. str .eq. 'true')
			else if (lbl_match(line, 'Coupled carrier solve')) then
				! "Coupled carrier solve: True|False", DEFAULT False -- solve
				! n(H2) (and the diffused element fractions) as Newton
				! unknowns beside the three hydrodynamic ones, instead of
				! the default alternation of a hydrodynamic solve with the
				! fixed-wind transport relaxations, judged jointly by the
				! certification of the refreshed state
				! (steady_wind_with_element_diffusion). MEASURED on the two
				! reloads of backup/regression (2026-09-11): the alternation
				! brings every elemental wind
				! row of the HD 209458 b reload inside 1e-5 in 12 passes
				! and the hot-Uranus H2 wind row from 7.4e-2 to 2.2e-2 in
				! 40 passes, while the coupled solve certifies neither
				! (stage 2, items N0 to N38). True keeps the coupled solve
				! as a control to measure against. On a molecular configuration it
				! requires "Molecular carrier transport: True", because H2
				! must be an unknown and not an eliminated variable;
				! checked below, after every key is parsed.
				! THREE VALUES, AND A WORD THAT IS NONE OF THEM STOPS THE
				! RUN. "On stall" is two words, so the value is read as the
				! pair (word 4, word 5); a silent fall-through to False
				! would leave an input file asking for the block and a run
				! that never enters it.
				str = get_word(line, 4)
				str2 = get_word(line, 5)
				if (str .eq. 'True' .or. str .eq. 'true') then
					carrier_in_newton = .true.
				else if (str .eq. 'False' .or. str .eq. 'false') then
					carrier_in_newton = .false.
				else if ((str .eq. 'On' .or. str .eq. 'on') .and.        &
				         (str2 .eq. 'stall' .or. str2 .eq. 'Stall')) then
					carrier_newton_on_stall = .true.
				else
					write(*,*) '(input_read) ERROR: "Coupled carrier'
					write(*,*) '  solve: '//trim(str)//' '//trim(str2)//'"'
					write(*,*) '  is not a value this key takes. It takes'
					write(*,*) '  False (the alternation of a hydrodynamic'
					write(*,*) '  solve with the fixed-wind transport'
					write(*,*) '  relaxations, the default), True (the'
					write(*,*) '  coupled block from the first pass), or'
					write(*,*) '  "On stall" (the alternation first, the'
					write(*,*) '  block from the pass at which the'
					write(*,*) '  alternation stops approaching a joint'
					write(*,*) '  fixed point). Aborting.'
					error stop 1
				endif
			else if (lbl_match(line, 'Oxygen transport')) then
				! RETIRED 2026-09-02. The operator this key switched was
				! never about oxygen: H2 is its first carrier and is
				! transported with or without the oxygen cycle, so the name
				! stated something the code does not do. Refused rather than
				! aliased, because a silent alias would leave stored input
				! files describing a switch that no longer exists.
				write(*,*) '(input_read) ERROR: "Oxygen transport" is no'
				write(*,*) '  longer a key. It was renamed to "Molecular'
				write(*,*) '  carrier transport: True|False": the operator'
				write(*,*) '  transports H2 whether or not the oxygen cycle'
				write(*,*) '  is on, so the old name named the wrong thing.'
				write(*,*) '  Replace the line in input.inp. Aborting.'
				error stop 1
			else if (lbl_match(line, 'Stellar FUV B1 flux')) then
				! RETIRED 2026-09-06. Band B1 was 1110-1201 A, between the
				! Lyman-Werner interval and Ly-alpha. The H2 Lyman and
				! Werner lines pump on both sides of 1110 A at the
				! temperature of a planetary base, so an edge there
				! normalized the H2 pumping per photon of a band narrower
				! than the one the lines absorb from, and the self-shielding
				! table rated 45 percent more absorptions than the beam
				! lost. B1 is now part of the Lyman-Werner band, 912-1201 A.
				! Refused rather than aliased: the two numbers have to be
				! added, which the code cannot do for a file that states
				! only one of them.
				write(*,*) '(input_read) ERROR: "Stellar FUV B1 flux" is'
				write(*,*) '  no longer a key. On 2026-09-06 band B1'
				write(*,*) '  (1110-1201 A) was merged into the'
				write(*,*) '  Lyman-Werner band, which is now 912-1201 A,'
				write(*,*) '  because the H2 pumping lines cross 1110 A at'
				write(*,*) '  the temperature of a planetary base.'
				write(*,*) '  To restate: delete this line and set'
				write(*,*) '  "Stellar LW flux" to the sum of the two, i.e.'
				write(*,*) '  the stellar flux at the planet integrated'
				write(*,*) '  over 912-1201 A. Aborting.'
				error stop 1
			else if (lbl_match(line, 'Stellar FUV B3 flux')) then
				! "Stellar FUV B3 flux [erg/cm2/s]: <F>" -- 1231-1450 A.
				str = get_word(line, 6)
				read(str,*) F_FUV_B3
				fuv_b3_flux_stated = .true.
			else if (lbl_match(line, 'Stellar FUV B4 flux')) then
				! "Stellar FUV B4 flux [erg/cm2/s]: <F>" -- 1451-2304 A.
				! B2 is the Ly-alpha line and is supplied by "Stellar Lya
				! flux"; the 912-1201 A band is supplied by "Stellar LW
				! flux". The edges are fixed by the H2O branching ratios and
				! by the Ly-alpha line (water_photolysis.f90), not chosen.
				str = get_word(line, 6)
				read(str,*) F_FUV_B4
				fuv_b4_flux_stated = .true.
			else if (lbl_match(line, 'Stellar LW flux')) then
				! "Stellar LW flux [erg/cm2/s]: <F>" -- band-integrated
				! stellar flux over 912-1201 A at the planet's orbit. Drives
				! H2 photodissociation in the molecular network
				! (lyman_werner.f90) AND, with the oxygen chemistry on, the
				! H2O and OH photolysis of the same interval: the three
				! absorbers share one beam, so the interval has one flux.
				! 0 = off (default).
				str = get_word(line, 5)
				read(str,*) F_LW_star
				lw_flux_stated = .true.
			else if (lbl_match(line, 'Molecular base')) then
				! EOS-only molecular base (docs/lower_atmosphere_*).
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true') molecular_base = .true.
			else if (lbl_match(line, 'Lower atmosphere profile')) then
				! Lower-atmosphere solution handed over as a table over an
				! interval of pressure instead of the single-level scalars of
				! base.inp.
				! "Lower atmosphere profile: lower_atmosphere_profile.dat".
				! Tested BEFORE 'Lower atmosphere', whose label is a prefix of
				! this one (longest / most-specific first, as for the two
				! "Wind-AE seed" keys).
				lap_file = get_word(line, 4)
			else if (lbl_match(line, 'Lower atmosphere')) then
				! Lower-atmosphere pre-step (docs/lower_atmosphere_*).
				str = get_word(line, 3)
				if (str .eq. 'vulcan')   lower_atm_mode = 2
				if (str .eq. 'analytic') lower_atm_mode = 1
				if (str .eq. 'none')     lower_atm_mode = 0
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) lower_atm_r1bar
			else if (lbl_match(line, 'Lower column')) then
				! analytic lower column: "Lower column: <R_1bar in R_J>"
				str = get_word(line, 3);  read(str,*) lower_col_r1bar
			else if (lbl_match(line, 'He_Kzz')) then
				! Eddy diffusion coefficient [cm^2/s]. It fills every entry of
				! kzz_cell, which the element transport (He_diffusion) and the
				! molecular carrier transport both read, so it is inert only
				! when neither of those runs. "He_Kzz: 1.0e9"
				str = get_word(line, 2);  read(str,*) he_kzz
			else if (lbl_match(line, 'He_alphaT')) then
				! Thermal-diffusion factor alpha_T for He (P2c). "He_alphaT: 0.0"
				str = get_word(line, 2);  read(str,*) he_alphaT
			else if (lbl_match(line, 'He_ambipolar')) then
				! Ambipolar-corrected settling mass (P2b). "He_ambipolar: False"
				str = get_word(line, 2)
				if (str .eq. 'False' .or. str .eq. 'false') he_ambipolar = .false.
			else if (lbl_match(line, 'He_metal_diffusion')) then
				! Diffuse trace metals too (P2d). "He_metal_diffusion: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') he_metal_diffusion = .true.
			else if (lbl_match(line, 'He_diffusion')) then
				! He/H diffusive separation (default off).
				! "He_diffusion: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') he_diffusion = .true.
			else if (lbl_match(line, 'Stall')) then
				! Stall-detector override: "Stall [tol,N]: <rel_tol> <N_steps>"
				! (smaller tol and/or larger N = harder to declare a plateau)
				str = get_word(line, 3);  read(str,*) stall_tol
				str = get_word(line, 4);  read(str,*) N_stall
				write(*,'(A,ES9.2,A,I0)') ' (input_read) Stall override: tol =', &
				                          stall_tol, ', N =', N_stall
			else if (lbl_match(line, 'Energy solver')) then
				! "Energy solver: Explicit" reverts to the original forward-
				! Euler source update (solver-component isolation tests).
				str = get_word(line, 3)
				if (str .eq. 'Explicit') then
					use_semi_implicit_energy = .false.
					write(*,*) '(input_read) Energy solver: explicit forward Euler'
				endif
			else if (lbl_match(line, 'Run mode')) then
				! WHAT THIS RUN IS DOING, stated rather than inferred.
				! "Run mode: init" = initialization / continuation, no claim
				! about elapsed time; "Run mode: phys" = physical
				! integration, one global dt per step and a clock that
				! advances only on an accepted step. init is the default in
				! every configuration; phys is asked for. The refusal rule
				! that goes with it is applied after every key is read.
				str = get_word(line, 3)
				if (str .eq. 'init' .or. str .eq. 'Init') then
					run_mode = run_mode_init;  run_mode_given = .true.
				else if (str .eq. 'phys' .or. str .eq. 'Phys') then
					run_mode = run_mode_phys;  run_mode_given = .true.
				else
					write(*,*) '(input_read) ERROR: "Run mode" takes '//   &
					           'init or phys, not "'//trim(str)//'"'
					error stop 1
				endif
			else if (lbl_match(line, 'Restart intent')) then
				! "Restart intent: trajectory | relaxation | stationary
				!  [evaluate|equilibrate]"
				! Only meaningful for a run that loads a state; the
				! consistency block below refuses it otherwise, and
				! refuses the two combinations that contradict themselves.
				str = get_word(line, 3)
				if (str .eq. 'trajectory') then
					restart_intent = restart_intent_trajectory
				else if (str .eq. 'relaxation') then
					restart_intent = restart_intent_relaxation
				else if (str .eq. 'stationary') then
					restart_intent = restart_intent_stationary
				else
					write(*,*) '(input_read) ERROR: "Restart intent" takes'
					write(*,*) '  trajectory, relaxation or stationary, not'
					write(*,*) '  "'//trim(str)//'". Aborting.'
					error stop 1
				endif
				restart_intent_given = .true.
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) then
					if (restart_intent .ne. restart_intent_stationary) then
						write(*,*) '(input_read) ERROR: "Restart intent"'
						write(*,*) '  takes a second word only with'
						write(*,*) '  stationary; "'//trim(str)//'" stands'
						write(*,*) '  beside another intent. Aborting.'
						error stop 1
					endif
					if (str .eq. 'evaluate') then
						stationary_evaluate_only = .true.
					else if (str .eq. 'equilibrate') then
						stationary_equilibrate_loaded = .true.
					else
						write(*,*) '(input_read) ERROR: the second word of'
						write(*,*) '  "Restart intent: stationary" is'
						write(*,*) '  evaluate or equilibrate, not "'//     &
						           trim(str)//'". Aborting.'
						error stop 1
					endif
				endif
			else if (lbl_match(line, 'Restart option change')) then
				! "Restart option change: <token>[, <token> ...]"
				! (decision 21). The tokens are separated by commas, blanks
				! or tabs; each one is checked against the vocabulary of
				! the state file's '# options' line, and a token that
				! decides how many unknowns the state has is refused there
				! and not here, because naming it cannot make the file's
				! rows this run's rows.
				optline = adjustl(line)
				optline = optline(len('Restart option change')+1:)
				if (len_trim(optline) .gt. 0) then
					if (is_sep(optline(1:1))) optline(1:1) = ' '
				endif
				do iopt = 1, len_trim(optline)
					if (optline(iopt:iopt) .eq. ','  .or.                  &
					    optline(iopt:iopt) .eq. char(9))                   &
						optline(iopt:iopt) = ' '
				enddo
				if (len_trim(optline) .eq. 0) then
					write(*,*) '(input_read) ERROR: "Restart option'
					write(*,*) '  change" names no token. Name the'
					write(*,*) '  option tokens that are allowed to'
					write(*,*) '  differ between the state file and this'
					write(*,*) '  run, or remove the line. Aborting.'
					error stop 1
				endif
				ktok = 1
				do
					str = get_word(optline, ktok)
					if (len_trim(str) .eq. 0) exit
					call name_restart_option_change_token(trim(str))
					ktok = ktok + 1
				enddo
			else if (lbl_match(line, 'Time stepping')) then
				! "Time stepping: Local" = cell-by-cell pseudo-time steps
				! (steady-state convergence acceleration; not time-accurate).
				str = get_word(line, 3)
				if (str .eq. 'Local') then
					use_local_dt = .true.
					write(*,*) '(input_read) Time stepping: local (cell-by-cell) pseudo-dt'
				endif
			else if (lbl_match(line, 'Level tol')) then
				! "Level tol: <val>" overrides the mass-flux level-stability
				! tolerance (<= 0 disables the level gate; legacy stops).
				str = get_word(line, 3);  read(str,*) lev_th
				write(*,'(A,ES9.2)') ' (input_read) Level-stability tol =', lev_th
			else if (lbl_match(line, 'Solver')) then
				! "Solver: Newton [du_switch]" = marching warm-up until the
				! flux metric du (radial spread of rho*v*r^2) < du_switch
				! (default 1e-2), then the JFNK steady solve polishes to du<du_th.
				str = get_word(line, 2)
				if (str .eq. 'Newton') then
					use_newton_solver = .true.
					str = get_word(line, 3)
					if (len_trim(str) .gt. 0) read(str,*) newton_du_switch
					write(*,'(A,ES9.2)') ' (input_read) Solver: Newton, '// &
						'JFNK hand-off at du <', newton_du_switch
				endif
			else if (lbl_match(line, 'Valve eps')) then
				call retired_base_key('Valve eps')
			else if (lbl_match(line, 'Hydrostatic base')) then
				call retired_base_key('Hydrostatic base')
			else if (lbl_match(line, 'Well balanced')) then
				! "Well balanced: True|False" (default False) -- carry the
				! DEPARTURE from the cell's own hydrostatic equilibrium
				! through the reconstruction, the Riemann jumps and the
				! pressure force, so that the equilibrium's flux difference
				! and its source cancel in the algebra instead of in
				! floating point (Kaeppeli and Mishra 2016, A&A 587, A94).
				! Exact preservation needs a flux that resolves a stationary
				! contact, which ROE and HLLC do and LLF does not.
				str = get_word(line, 3)
				well_balanced = (str .eq. 'True' .or. str .eq. 'true')
			else if (lbl_match(line, 'Low Mach velocity jump')) then
				! "Low Mach velocity jump: True|False" (default False) --
				! on the ROE branch, scale the normal velocity jump of the
				! Roe dissipation by min(|U_Roe|/a_Roe, 1), so that the
				! artificial viscosity of the momentum is of the order of
				! the momentum update instead of one order in the Mach
				! number larger (Rieper 2011, J. Comput. Phys. 230, 5263,
				! his eq. 3.15-3.16).  An accuracy correction for the
				! low-Mach regime, not a stiffness one, and derived for the
				! Roe flux alone; it has no effect under HLLC or LLF.
				! The value is word 5: the key itself is four words.
				str = get_word(line, 5)
				low_mach_velocity_jump = (str .eq. 'True' .or. str .eq. 'true')
			else if (lbl_match(line, 'Shapiro filter')) then
				str = get_word(line, 3);  read(str,*) shapiro_eps
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) shapiro_every
				if (shapiro_eps .gt. 0.0d0) write(*,'(A,ES9.2,A,I0,A)')  &
				   ' (input_read) Shapiro filter eps =', shapiro_eps,    &
				   ', every ', shapiro_every, ' steps'
			else if (lbl_match(line, 'Low-Mach damping')) then
				! "Low-Mach damping: <eps4> [<M_th>]" -- gated fourth-
				! difference dissipation of the 2 dr contact/entropy mode
				! the HLLC flux stops damping as v -> 0. eps4 <= 0 disables.
				! See src/modules/flux/low_mach_dissipation.f90.
				str = get_word(line, 3);  read(str,*) lowmach_damp_eps
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) lowmach_damp_mach_th
				if (lowmach_damp_eps .gt. 0.0d0) then
					write(*,'(A,ES9.2,A,ES9.2)') ' (input_read) Low-Mach'// &
					   ' contact-mode damping eps4 =', lowmach_damp_eps,  &
					   ', M_th =', lowmach_damp_mach_th
					if (lowmach_damp_mach_th .le. 0.0d0)                  &
					   write(*,*) '(input_read.f90) WARNING: "Low-Mach '//&
					   'damping" M_th <= 0 switches the term off everywhere.'
				endif
			else if (lbl_match(line, 'Base BC')) then
				str = get_word(line, 3)
				if (str .eq. 'density')  base_bc_mode = 0
				if (str .eq. 'pressure') base_bc_mode = 1
				if (base_bc_mode .eq. 1) then
					str = get_word(line, 4)
					if (len_trim(str) .gt. 0) read(str,*) base_p_ubar
					write(*,'(A,ES9.2,A)')                               &
					   ' (input_read) Base BC: pressure-anchored, '//    &
					   'p_base =', base_p_ubar, ' microbar (n0 derived)'
				else
					write(*,'(A)') ' (input_read) Base BC: density '//   &
					   '(legacy, n0 from input)'
				endif
			else if (lbl_match(line, 'Base ghost temperature')) then
				call retired_base_key('Base ghost temperature')
			else if (lbl_match(line, 'Max steps')) then
				! "Max steps: <N>" overrides the hard cap on marching
				! iterations (default 1000000).
				str = get_word(line, 3);  read(str,*) marching_step_max
				write(*,'(A,I0)') ' (input_read) Max marching steps =', marching_step_max
			else if (lbl_match(line, 'Coronal cutoff width')) then
				! "Coronal cutoff width: <w>" sets the roll-off width of the
				! coronal-excitation guard below the 1e3 K CHIANTI fit floor
				! (Cool_coeff.f90). Default 0.1. The earlier 0.08-0.13 window
				! is superseded: with the C/N/O ground-term fine-structure
				! floors solved in statistical equilibrium no coronal cooling
				! survives below the fit floor for the guard to remove, and the
				! base-cell balance temperature of all four paper planets is
				! identical over w = 0.02-1.2
				! (docs/coronal_cutoff_width.md section 7.2).
				str = get_word(line, 4);  read(str,*) coronal_cutoff_width
				write(*,'(A,F6.3)') ' (input_read) Coronal excitation cutoff'// &
				   ' width w =', coronal_cutoff_width
			else if (lbl_match(line, 'Base IR field')) then
				! "Base IR field: True|False" lets the infrared coolants of a
				! molecular layer -- the eight ground-term fine-structure lines
				! of C I, C II, N II, O I and the H3+ bands -- see the thermal
				! radiation of the lower atmosphere (a blackbody at T0,
				! geometrically diluted) instead of only emitting into vacuum,
				! so each stops cooling at its own radiative equilibrium
				! temperature. Default False. See fine_structure_line_transfer
				! (Cool_coeff.f90) and h3p_net_cooling_rate (h3p_cooling.f90).
				str = get_word(line, 4)
				if (str .eq. 'True' .or. str .eq. 'true') base_ir_field = .true.
				if (base_ir_field) write(*,'(A)') ' (input_read) Base IR '//   &
				   'field on: infrared coolants see B_nu(T0) from below'
			else if (lbl_match(line, 'Molecular IR bands')) then
				! "Molecular IR bands: True|False" adds the infrared coolants a
				! real H2 atmosphere carries below the H2 -> H front and the
				! code did not: the H2 quadrupole and magnetic dipole line
				! spectrum, and the H2O and CO vibration-rotation bands. Each
				! emits in LTE and absorbs the same diluted B_nu(T0) the
				! `Base IR field` closure supplies, so each stops cooling at
				! its own radiative equilibrium temperature instead of running
				! the layer down. Default False.
				! H2O and CO need `Oxygen chemistry: True` to exist at all; H2
				! needs `Molecular chemistry: True`. See
				! molecular_infrared_cooling.f90.
				str = get_word(line, 4)
				if (str .eq. 'True' .or. str .eq. 'true') mol_ir_bands = .true.
				if (mol_ir_bands) write(*,'(A)') ' (input_read) Molecular IR'// &
				   ' bands on: H2 lines + H2O/CO bands exchange with B_nu(T0)'
			else if (lbl_match(line, 'Molecular reaction heat')) then
				! "Molecular reaction heat: True|False" deposits the energy
				! the collisional reactions of the H2/He network release,
				! built from one species-enthalpy table so that a closed
				! chemical cycle releases exactly zero. The photon-driven
				! reactions, the radiative recombinations and the collisional
				! ionizations are excluded: their energy is already in the
				! ledger. The term lives inside the molecular network
				! (`with_molecules` in util_ion_eq): an atomic gas has no
				! collisional molecular reactions, so there the key is
				! inert rather than refused. DEFAULT True --
				! the reason is at the declaration in parameters.f90 -- so
				! this key exists to turn the term OFF, for A/B work against
				! the state the code had before it.
				str = get_word(line, 4)
				mol_reaction_heat = (str .eq. 'True' .or. str .eq. 'true')
				if (.not. mol_reaction_heat) write(*,'(A)') ' (input_read)'//&
				   ' Molecular reaction heat OFF: the collisional network'// &
				   ' chemical heating is NOT deposited'
			else if (lbl_match(line, 'H2 double ionization')) then
				! "H2 double ionization: off|chung80|yan_rho" resolves
				! H2 + hv -> H+ + H+ + 2e- as a channel of its own. It has a
				! vertical threshold of 51.4 eV (Yan, Sadeghpour & Dalgarno
				! 1998 sec. 4) and releases TWO protons per event, where the
				! single dissociative channel that would otherwise
				! absorb it releases one. DEFAULT 'chung80', the 20 percent the
				! source states outright; the key exists to SELECT THE
				! MODEL or to TURN THE CHANNEL OFF, not to turn it on.
				! 'off' reproduces a pre-E1 result and is not a physical
				! statement -- the reaction happens. Neither model is a
				! measurement of the branching itself; 'chung80' and
				! 'yan_rho' rest on different published numbers and
				! disagree by a factor 1.6 at 110 eV, so running both is
				! how the uncertainty is reported. See the model list in
				! h2_photo_channels.f90.
				str = get_word(line, 4)
				if (str .eq. 'off' .or. str .eq. 'chung80' .or.            &
				    str .eq. 'yan_rho') then
					h2_double_ionization = str
				else
					write(*,*) '(input_read) ERROR: "H2 double'
					write(*,*) '  ionization" must be one of off,'
					write(*,*) '  chung80, yan_rho. Got: '//trim(str)
					error stop 1
				endif
				if (h2_double_ionization .ne. 'off') then
					write(*,'(A)') ' (input_read) H2 double ionization'//   &
					   ' model '//trim(h2_double_ionization)//': H+ + H+'// &
					   ' + 2e- resolved above 51.4 eV'
				else
					write(*,'(A)') ' (input_read) WARNING H2 double'//      &
					   ' ionization OFF: the channel is folded back into'
					write(*,'(A)') '   the single dissociative one, which'//&
					   ' gives each event one proton where it makes two'
					write(*,'(A)') '   and an H atom it does not make.'//   &
					   ' Pre-E1 behavior, kept for reproduction only.'
				endif
			else if (lbl_match(line, 'H2 neutral dissociation')) then
				! "H2 neutral dissociation: True|False" resolves
				! H2 + hv -> H + H, the absorptions that make no ion. It is
				! nonzero ONLY over 33-41 eV, the window in which Chung,
				! Lee, Masuoka & Samson (1993) Table 1 measures a
				! photoionization yield below unity; outside it their source
				! ASSUMES unit yield rather than measuring it, so the
				! channel is set to zero there. DEFAULT True: sigma_n is a
				! measurement. The key exists to TURN THE CHANNEL OFF, and
				! False folds that share back into the ionizing channels in
				! their own proportion -- the unit yield the code assumed
				! before, kept for reproducing a pre-E1 result.
				str = get_word(line, 4)
				h2_neutral_dissociation =                                  &
					(str .eq. 'True' .or. str .eq. 'true')
				if (.not. h2_neutral_dissociation) then
					write(*,'(A)') ' (input_read) WARNING H2 neutral'//     &
					   ' dissociation OFF: over 33-41 eV up to 7.4 per'
					write(*,'(A)') '   cent of the absorptions are given'// &
					   ' an H2+ and an electron the event does not make.'
				endif
			else if (lbl_match(line, 'Base velocity')) then
				call retired_base_key('Base velocity')
			else if (lbl_match(line, 'Base grid')) then
				! "Base grid [dr,cells]: <dr_base> [<N_low_cells>]" sets the
				! resolution of the uniform region of the Mixed grid: N_low_cells
				! cells of size dr_base [R_p] stacked on the base. The two numbers
				! are NOT independent -- their product is the radial extent of the
				! uniform region (0.01 R_p by default) -- so they share one line,
				! as "du_th [PLM,WENO3]" does. Refining the base at fixed extent
				! means dividing dr_base and multiplying N_low_cells by the same
				! factor: "Base grid [dr,cells]: 5.0e-5 200" is the 4x refinement.
				! Requirement: dr_base must resolve the base scale height
				! H = kT/(mu g). It is checked at startup by
				! write_setup_report, which prints H(T_eq)/dr in cells and
				! warns below 10; the undamped stationary 2 dr entropy mode
				! sets in near 5. Ignored
				! by the Uniform and Stretched grid types.
				! The key is the only writer of dr_base_from_key: the
				! resolved-configuration record says "key" or "default"
				! from it, and the two cannot be told apart by value (a
				! key may state the default's own digits).
				str = get_word(line, 4);  read(str,*) dr_base
				dr_base_from_key = .true.
				str = get_word(line, 5)
				if (len_trim(str) .gt. 0) read(str,*) N_low_cells
				write(*,'(A,ES9.2,A,I0,A,ES9.2,A)')                       &
				   ' (input_read) Base grid: dr =', dr_base,              &
				   ' R_p x ', N_low_cells, ' cells (uniform region ',     &
				   dr_base*N_low_cells, ' R_p)'
			else if (lbl_match(line, 'Grid cells')) then
				! "Grid cells: <N>" sets the number of computational cells of
				! the radial domain (ghost cells are added on top and are not
				! counted). Omitting the key keeps the 500 cells that used to
				! be a compile-time constant, so an existing input.inp is
				! unaffected. For the Mixed grid the split between the uniform
				! base region and the stretched region is set separately by
				! "Base grid [dr,cells]:", and define_grid checks that the two
				! are compatible.
				str = get_word(line, 3);  read(str,*) N
				if (N .lt. 10) then
					write(*,'(A,I0,A)') ' (input_read.f90) ERROR: "Grid '//  &
					   'cells: ', N, '" is not a usable number of '//        &
					   'computational cells (at least 10 are needed).'
					error stop 1
				endif
				write(*,'(A,I0,A)') ' (input_read) Grid cells: ', N,         &
				   ' computational cells'
			else if (lbl_match(line, 'Viscosity')) then
				! "Viscosity: True" selects the calibrated mu(T) (Watson+1981
				! conductivity through the monatomic Chapman-Enskog relation)
				! and enables the viscous dissipation q_mu; the numeric form
				! "Viscosity: <mu0> [<s>]" keeps its meaning, a diagnostic
				! power law mu = mu0*T^s in code units. Both are one-word
				! keys, so the value sits at word 2 (as for "Solver"/"CFL").
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') then
					visc_on = .true.
					write(*,'(A)') ' (input_read) Viscosity: calibrated '// &
					   'mu(T) (Chapman-Enskog / Watson+1981), q_mu included'
				else if (str .eq. 'False' .or. str .eq. 'false' .or.      &
				         len_trim(str) .eq. 0) then
					visc_on = .false.;  visc_mu0 = 0.0d0
				else
					read(str,*) visc_mu0
					str = get_word(line, 3)
					if (len_trim(str) .gt. 0) read(str,*) visc_s
					if (visc_mu0 .gt. 0.0d0) write(*,'(A,ES9.2,A,F5.2)')  &
					   ' (input_read) Viscosity power-law override mu0 =',  &
					   visc_mu0, ', s =', visc_s
				endif
			else if (lbl_match(line, 'Photoelectron heating')) then
				str = get_word(line, 3)
				if (str .eq. 'full') then
					photoheat_full_photon_energy = .true.
					photoheat_photon_fraction    = 1.0d0
				else if (str .eq. 'excess') then
					photoheat_full_photon_energy = .false.
					photoheat_photon_fraction    = -1.0d0
				else
					read(str,*,iostat=ios_pe) photoheat_photon_fraction
					if (ios_pe .ne. 0 .or. photoheat_photon_fraction .le. 0.0d0  &
					    .or. photoheat_photon_fraction .gt. 1.0d0) then
						write(*,*) '(input_read.f90) ERROR: "Photoelectron heating" takes excess, full or a fraction in (0,1], not "'//trim(str)//'"'
						error stop 1
					endif
					photoheat_full_photon_energy = (photoheat_photon_fraction .ge. 1.0d0)
				endif
			else if (lbl_match(line, 'Caloric EOS')) then
				str = get_word(line, 3)
				if (str .eq. 'monatomic') then
					caloric_eos_monatomic = .true.
				else if (str .eq. 'ladder') then
					caloric_eos_monatomic = .false.
				else
					write(*,*) '(input_read.f90) ERROR: "Caloric EOS" takes ladder or monatomic, not "'//trim(str)//'"'
					error stop 1
				endif
			else if (lbl_match(line, 'Conduction')) then
				! "Conduction: True" enables heat conduction with the
				! Watson et al. (1981) atomic-hydrogen kappa(T).
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') cond_on = .true.
				if (str .eq. 'False' .or. str .eq. 'false') cond_on = .false.
				if (cond_on) write(*,'(A)') ' (input_read) Heat conduction: '// &
				   'kappa(T) = 4.45e4 (T/1000 K)^0.7 erg/cm/s/K'
			else if (lbl_match(line, 'Resid tol')) then
				! "Resid tol: <val>" = converge on the steady residual ||R||
				! instead of du (<= 0 disables; legacy du-based stop).
				str = get_word(line, 3);  read(str,*) resid_th
				write(*,'(A,ES9.2)') ' (input_read) Residual-based convergence, tol =', resid_th
			else if (lbl_match(line, 'Resid norm')) then
				! RETIRED (section 145). The residual norm is the maximum over
				! cells of a cell's own scaled residual and has no second
				! form, so this key selects nothing. It is still matched, and
				! warned about, so that an input file carrying it says so
				! instead of failing the unknown-key check.
				write(*,'(A)') ' (input_read) WARNING: "Resid norm" is'//    &
				   ' retired and ignored -- the residual norm is the'//      &
				   ' maximum over cells of a'
				write(*,'(A)') '     cell''s own scaled residual, and there'//&
				   ' is no other form (Update_EXHALE_stage1.pdf section 145).'
			else if (lbl_match(line, 'Flux spread tol')) then
				! "Flux spread tol: <tol> [<r_flux [R_p]>]" -- the FLUX gate
				! (section 133). The steady solve is accepted only when the
				! radial spread of the Riemann FACE mass flux over
				! r >= r_flux is below <tol> (section 145)
				! AND the residual is below "Resid tol". <tol> <= 0 disables
				! the flux gate and leaves the residual alone in charge.
				str = get_word(line, 4);  read(str,*) flux_spread_th
				! A pre-section-145 value. The gate used to measure the
				! CELL-CENTRED rho*v*r^2 and its tolerance was 5e-3; it now
				! measures the Riemann face flux, on which every state that
				! stops on du alone is already below 5e-3. Warned, not
				! refused: the value may be deliberate.
				if (flux_spread_th .ge. 1.0d-3) then
					write(*,'(A,ES9.2,A)') ' (input_read) WARNING: "Flux'//  &
					   ' spread tol" =', flux_spread_th, ' is a pre-145'//  &
					   ' value. The gate now measures the'
					write(*,'(A)') '     Riemann FACE mass flux, on which'//&
					   ' converged states sit at 1e-10 to 4e-6 and states'//&
					   ' that stopped on du'
					write(*,'(A,ES9.2,A)') '     alone sit at 1e-4 to'//    &
					   ' 3e-2, so this admits states the gate is meant'//   &
					   ' to refuse (default', flux_spread_th_default, ').'
				endif
				str = get_word(line, 5)
				if (len_trim(str) .gt. 0) read(str,*) r_flux
				if (flux_spread_th .gt. 0.0d0) then
					write(*,'(A,ES9.2,A,F6.3,A)') ' (input_read) Flux gate: '// &
					   'spread of the face mass flux <', flux_spread_th,        &
					   ' over r >=', r_flux, ' R_p'
				else
					write(*,'(A)') ' (input_read) Flux gate disabled '//        &
					   '(Flux spread tol <= 0)'
				endif
			else if (lbl_match(line, 'CFL')) then
				! Override the CFL number ("CFL: <value>"); lower = smaller dt.
				str = get_word(line, 2);  read(str,*) CFL
			else if (lbl_match(line, 'Transonic IC')) then
				str = get_word(line, 3)
				if (str .eq. 'True') transonic_ic = .true.
			else if (lbl_match(line, 'Hot Parker IC')) then
				str = get_word(line, 4)
				read(str,*) T_wind_ic
				hot_parker_ic = .true.
			else if (lbl_match(line, 'IC mode')) then
				! "IC mode: <cold|transonic|hot_parker|auto>". The named
				! modes are synonyms for the legacy keys; 'auto' defers the
				! choice to select_IC_auto (set_IC.f90), which probes the
				! cold sonic-point topology of the actual potential. The
				! explicit legacy keys take precedence over 'auto'.
				str = get_word(line, 3)
				if (str .eq. 'cold') then
					ic_mode = 0
				else if (str .eq. 'transonic') then
					ic_mode = 1
					transonic_ic = .true.
				else if (str .eq. 'hot_parker') then
					ic_mode = 2
					hot_parker_ic = .true.
				else if (str .eq. 'auto') then
					ic_mode = 3
				else if (str .eq. 'windae') then
					! In-process Wind-AE warm-start IC: init.f90 calls the
					! ported Wind-AE solver (src/modules/wind_ae/) to build
					! an IC on the EXHALE grid, then loads it. No separate
					! wind_ae_ic.x run needed. Standalone tool is unaffected.
					ic_mode = 4
				else
					write(*,*) '(input_read.f90) WARNING: unknown "IC ' // &
					           'mode: ', str, '"; using cold hydrostatic.'
					ic_mode = 0
				endif
			else if (lbl_match(line, 'Newton solver')) then
				str = get_word(line, 3)
				if (str .eq. 'False') use_newton_ieq = .false.
			else if (lbl_match(line, 'Brent solver')) then
				str = get_word(line, 3)
				if (str .eq. 'False') use_brent_tsolve = .false.
			endif
		enddo

		! ----- Unknown-line warning -----
		! Now that every core and keyword line has been consumed, flag any
		! non-blank, non-'#' line that matches NO known label (and is not the
		! energy-band "[E_low" line). Catches typos like "Newton Solver:"
		! (capital S) that anchored matching silently ignores. WARN only --
		! never stop, so a stray line never aborts the run.
		do i = 1, nlines
			if (len_trim(filelines(i)) .eq. 0) cycle
			line = adjustl(filelines(i))
			if (line(1:1) .eq. '#') cycle
			if (line(1:6) .eq. '[E_low') cycle
			is_known = .false.
			do kk = 1, size(known_keys)
				if (lbl_match(filelines(i), trim(known_keys(kk)))) then
					is_known = .true.
					exit
				endif
			enddo
			if (.not. is_known)                                          &
				write(*,*) '(input_read.f90) WARNING: unrecognized input '// &
				   'line (matches no known key): '//trim(line)
		enddo

		! ----- Low-Mach damping vs the explicit stability bound -----
		! The 2 dr mode decays at 16 eps4 lambda/dr while the marching step is
		! dt = CFL dr/lambda, so the damping factor per step is 16 eps4 CFL and
		! the explicit update is unstable beyond 1. Checked here, after the
		! keyword loop, because "CFL:" may appear on either side of this key.
		if (lowmach_damp_eps .gt. 0.0d0 .and.                            &
		    16.0d0*lowmach_damp_eps*CFL .gt. 1.0d0) then
			write(*,'(A,ES9.2,A,F7.4,A)') ' (input_read.f90) WARNING: '// &
			   '"Low-Mach damping" eps4 =', lowmach_damp_eps,            &
			   ' exceeds the explicit stability bound 1/(16 CFL) =',      &
			   1.0d0/(16.0d0*CFL), '; the marching loop may diverge.'
		endif
		if (coronal_cutoff_width .le. 0.0d0) then
			write(*,*) '(input_read.f90) ERROR: "Coronal cutoff width" must '// &
			   'be > 0 (it divides the fractional temperature deficit).'
			error stop 1
		endif

		! The transonic-wind IC already satisfies steady mass conservation
		! (rho*v*r^2 = const), so du ~ 0 at step 0 would trip the "momentum
		! constant" exit before the cold wind heats to its hot steady state.
		! Force the first iterations so the heating develops first (the normal
		! convergence test then resumes; see EXHALE_main.f90). Skipped in
		! post-processing-only runs (force_start/do_only_pp are exclusive).
		if ((transonic_ic .or. hot_parker_ic) .and. .not. do_only_pp) &
			force_start = .true.

		! Convert the stellar radius to cm (Balmer dilution R_star/a_orb) and
		! enable the excited-H coupling only when both stellar inputs are set.
		R_star        = R_star*Rsun
		use_excited_H = (T_star_eff .gt. 0.0d0) .and. (R_star .gt. 0.0d0)

		! "Spectrum type: Planck" builds the flux on EVERY point of the
		! photon grid out of pi B_nu(T_eff) (R_star/a)^2, so without both
		! stellar quantities there is no spectrum at all -- not a default
		! one. Refused rather than run on a zero field.
		if (sp_type .eq. 'Planck' .and. .not. use_excited_H) then
			write(*,*) '(input_read.f90) ERROR: "Spectrum type: Planck"'
			write(*,*) '  needs "Stellar Teff [K]:" > 0 and'
			write(*,*) '  "Stellar radius [R_sun]:" > 0. They are the only'
			write(*,*) '  source of the field pi B_nu(T_eff) (R_star/a)^2'
			write(*,*) '  that this type puts on the whole photon grid.'
			write(*,*) '  State both lines, or select another spectrum'
			write(*,*) '  type. Aborting.'
			error stop 1
		endif

		! Guard: the in-line Ly-alpha escape-probability RT ("Jlya escape-prob:
		! True", jlya_mode = 2) builds J_lya = J_int + J_star, where the stellar
		! beam J_star is proportional to F_Lya_star (set by "Stellar Lya flux
		! [erg/cm2/s]:"). If F_Lya_star is left at its default 0, the stellar
		! beam vanishes and the "Lya stellar halfwidth"/"Lya stellar boost"
		! settings become silent no-ops. Refuse to run that misconfiguration
		! rather than produce incomplete physics without warning.
		if (jlya_mode .eq. 2 .and. F_Lya_star .le. 0.0d0) then
			write(*,*) '(input_read.f90) ERROR: "Jlya escape-prob: True" was'
			write(*,*) '  set but "Stellar Lya flux [erg/cm2/s]:" is missing or'
			write(*,*) '  <= 0. The stellar Ly-alpha beam (J_star) would be zero'
			write(*,*) '  and the halfwidth/boost settings would do nothing.'
			write(*,*) '  Add e.g. "Stellar Lya flux [erg/cm2/s]: 1.0e5" to'
			write(*,*) '  input.inp, or disable escape-prob mode. Aborting.'
			error stop 1
		endif

		! The absorbing lower boundary is a term of the escape-probability
		! closure only, so it does nothing for jlya_mode 0 (parameterized) or 1
		! (imported field). Say so rather than let the key look effective.
		if (lya_bottom_absorber .and. jlya_mode .ne. 2) then
			write(*,*) '(input_read.f90) WARNING: "Lya absorbing bottom: True"'
			write(*,*) '  only acts on the in-line escape-probability RT'
			write(*,*) '  ("Jlya escape-prob: True"); ignored for jlya_mode ='
			write(*,*) '  0 (parameterized) and 1 (imported J_Lya profile).'
		endif

	write(*,*) '(input_read.f90) Done'

   !------ Definition of physical parameters ------!
      
   if (base_density_key_given) n0 = 10.0**(n0)
   ! ---- Lower-atmosphere pre-step ("Lower atmosphere: vulcan|analytic").
   ! If requested and no base.inp exists yet, generate it now by invoking
   ! the generator (VULCAN photochemistry or the analytic column).
   ! The wind solve then proceeds on the produced base -- VULCAN as a
   ! subroutine.  Opt-out: omit the key (default off).
   call run_lower_atm_prestep

   ! ---- optional external base.inp (written by src/utils/run_lower.py or a
   ! lower-atmosphere model): overrides the base temperature, base radius
   ! [R_J], He/H ratio and eddy K_zz BEFORE the derived constants below.
   ! Absent file = no-op (byte-identical legacy).
   ! ---- optional lower-atmosphere PROFILE ("Lower atmosphere profile:").
   ! Read before base.inp so read_base_inp can refuse the scalar keys the
   ! profile now owns, and applied after it so nothing scalar can overwrite
   ! the profile.  Absent key = no-op (byte-identical legacy).
   ! Which elements the handoff states itself is recorded as the two readers
   ! run: set_element_abundance is the only door either of them uses, so the
   ! flag cannot drift from the abundance it belongs to. An element left
   ! false keeps whatever metals.inp gave it, and a restart of that element
   ! is not renormalized (load_IC).
   if (.not. allocated(melem_from_handoff)) allocate(melem_from_handoff(n_melem))
   melem_from_handoff = .false.
   call read_lower_atmosphere_profile
   call read_base_inp
   call apply_lower_atmosphere_profile

   ! ---- Composition reconciliation and validation. Placed here so every
   ! flag (thereis_He, HeH, thereis_HeITR, thereis_mol, he_diffusion,
   ! thereis_metals) has its final value: base.inp above can still flip
   ! thereis_He/HeH and the elemental abundances X_*, and the trailing
   ! keyword scan sets thereis_mol and he_diffusion. This is the single
   ! authoritative check before N_eq sizing.

   ! ---- Elemental reservoirs. The abundances reach this point from
   ! metals.inp (read at the top of input_read) and, for the elements the
   ! lower-atmosphere handoff carries, from the "<El>_H_base" keys of
   ! base.inp, which override metals.inp because they describe the
   ! composition at the base of THIS wind. Everything derived from the
   ! elemental reservoirs -- which elements are active, the canonical
   ! melem_ab array, and the below-threshold energy grid -- is therefore
   ! derived here, after the handoff, and nowhere else.

   ! Activate metal species if any metal abundance is set
   thereis_metals = (X_C .gt. 0.0d0) .or. (X_N .gt. 0.0d0)        &
                                      .or. (X_O .gt. 0.0d0)        &
                                      .or. (X_Mg .gt. 0.0d0)       &
                                      .or. (X_Si .gt. 0.0d0)       &
                                      .or. (X_Ca .gt. 0.0d0)       &
                                      .or. (X_Na .gt. 0.0d0)       &
                                      .or. (X_K  .gt. 0.0d0)       &
                                      .or. (X_S  .gt. 0.0d0)       &
                                      .or. (X_Fe .gt. 0.0d0)

   ! Abundances for each element in canonical element order (iel_*), so the
   ! grid/solver code can index metals by element rather than by named
   ! scalar. Extend this block (and the metals.inp reader) when adding
   ! elements.
   allocate(melem_ab(n_melem))
   melem_ab(iel_C)  = X_C
   melem_ab(iel_O)  = X_O
   melem_ab(iel_N)  = X_N
   melem_ab(iel_Mg) = X_Mg
   melem_ab(iel_Si) = X_Si
   melem_ab(iel_Ca) = X_Ca
   melem_ab(iel_Na) = X_Na
   melem_ab(iel_K)  = X_K
   melem_ab(iel_S)  = X_S
   melem_ab(iel_Fe) = X_Fe

   ! An active metal whose neutral ionization threshold lies below the
   ! 13.6 eV HI edge (e.g. Mg I at 7.646 eV) needs the photon grid carried
   ! below that edge. The floor is the LOWEST threshold over every active
   ! sub-Lyman absorber -- the He 2^3S metastable at 4.768 eV, a low-IP
   ! metal, and H(n=2) at 3.400 eV when the excited-hydrogen coupling is
   ! armed -- so the three are carried together and none of them excludes
   ! another (photon_grid_floor_eV in sed_read.f90; the band layout in
   ! set_energy_vectors.f90). Triggered by any such active element.
   do im = 1, n_melem
      if (melem_ab(im) .gt. 0.0d0 .and.                       &
          mion_ethr(melem_i0(im)) .lt. e_th_HI)               &
         thereis_lowIP_metal = .true.
   enddo


   ! Resolve which stage a two-stage run starts in, now that both du_th
   ! values are known. "PLM+WENO3" with a stage-1 threshold that does not
   ! exceed the stage-2 one has an EMPTY PLM stage -- the marching loop's
   ! in_plm_stage = (du_th_plm > du_th) is false from step 0 -- so the run is
   ! single-stage WENO3 and must be dispatched as WENO3. Leaving rec_method
   ! at the input word made the first Reconstruct() fall through to
   ! `case default` and abort with "unknown reconstruction scheme:
   ! PLM+WENO3" (Reconstruction.f90); "du_th [PLM,WENO3]: -1.0 1.0e-3" is
   ! the form that hit it (wasp_hybrid_finish, wasp_localdt_cont).
   if (recon_two_stage .and. du_th_plm .le. du_th) then
      rec_method = 'WENO3';  use_plm = .false.;  use_weno3 = .true.
      write(*,'(A,ES9.2,A,ES9.2,A)') ' (input_read) "PLM+WENO3" with a '//  &
           'stage-1 threshold ', du_th_plm, ' not above the stage-2 one ',  &
           du_th, ': the PLM stage is empty -> single-stage WENO3.'
   endif

   ! The continuation walks the PLM -> WENO3 hand-off; a run that has no
   ! hand-off has nothing for it to walk. Said rather than silently ignored,
   ! because a key that is present and inert is the kind of thing a reader of
   ! the input file has no way to notice.
   if (recon_lambda_step0 .gt. 0.0d0 .and.                                 &
       .not. (recon_two_stage .and. du_th_plm .gt. du_th)) then
      recon_lambda_step0 = 0.0d0
      write(*,'(A)') ' (input_read) WARNING: "Reconstruction continuation"'//&
           ' needs a two-stage run with a non-empty PLM stage'//            &
           ' ("Reconstruction scheme: PLM+WENO3" and two du_th values):'//  &
           ' ignored.'
   endif

   ! Remove HeITR chemistry if He is not included
   if (.not. thereis_He) thereis_HeITR = .false.

   ! molecular chemistry constraints: requires He. Trace metals are solved
   ! together with the molecular network (System_HeH_mol_metals), so the two
   ! are not exclusive. Neither is He/H diffusion any more: the element
   ! transport of binary_element_diffusion closes over the molecular carriers
   ! (Blanc friction, mean carrier mass and charge, mole-fraction driver).
   if (thereis_mol .and. .not. thereis_He) then
      write(*,*) '(input_read) ERROR: Molecular chemistry needs He/H>0.'
      error stop 1
   endif

   ! The Roe flux is derived for ONE ideal gas with a constant adiabatic
   ! index: its average state (a_avg from the averaged enthalpy) and the
   ! star-state estimates of speed_estimate_ROE both assume it, and a
   ! two-index Roe average needs the Vinokur-Montagne/Glaister extension the
   ! code does not carry. With molecular chemistry and the ladder caloric
   ! EOS the index varies across the H2 front, so that derivation does not
   ! hold there.
   if (flux .eq. 'ROE' .and. thereis_mol .and.                            &
       .not. caloric_eos_monatomic) then
      write(*,*) '(input_read) ERROR: "Numerical flux: ROE" conflicts'//  &
                 ' with "Molecular chemistry: True" while "Caloric EOS"'
      write(*,*) '  is not "monatomic": the Roe average and the star-state'
      write(*,*) '  estimate are derived for one constant adiabatic index,'
      write(*,*) '  which the molecular ladder EOS does not provide.'
      write(*,*) '  Use "Numerical flux: HLLC" or "LLF", or set'
      write(*,*) '  "Caloric EOS: monatomic".'
      error stop 1
   endif

   ! ---- oxygen chemistry (the A2 option) ----
   ! Name the key, name the other owner of the quantity, name the fix, stop.
   ! The option is a third
   ! producer of the base H2/H partition, so it joins the single-source rule
   ! that the lower-atmosphere profile and base.inp already follow rather
   ! than inventing one of its own.
   ! THE FUV BAND FLUXES, when the run did not state them and a numerical
   ! spectrum is available. Done here, BEFORE the oxygen-chemistry report
   ! and its warnings below, so that those judge the fluxes the run will
   ! actually use.
   !
   ! H2 photodissociation and the H2O/OH photolysis are not options of the
   ! physics: the bands exist whenever the star does. What was optional was
   ! our knowing the numbers, and a run with a spectrum file knows them --
   ! so each band is integrated by the documented prescription (the band's
   ! interval of the SED, at the planet) instead of being left at zero and
   ! silently switching the channel off. A stated key always wins: it is
   ! the more specific statement, and a band-integrated measurement can be
   ! better than our trapezoid over whatever grid the file happens to have
   ! -- and a stated ZERO wins too: it is the run saying the band is not to
   ! be used (a comparison with a model that has no H2 photodissociation,
   ! or no continuum photolysis), not the run not knowing. Band B2 is the
   ! Ly-alpha line and stays a key ("Stellar Lya flux"): a line flux
   ! reconstructed from observations is a better number than a trapezoid
   ! over the file's rows across 1202-1230 A. The Lyman-Werner band is
   ! integrated for any molecular run (H2 absorbs it); B3 and B4 only with
   ! the oxygen chemistry, the only consumer of those two.
   ! Serial, once, before any parallel sweep: the opt-in cell-dump paths of
   ! the constrained chemical equilibrium and the Newton validation switch
   ! (review P1/P2, 2026-09-12: both used to be read lazily inside the sweep).
   call cce_dump_paths_read()
   call newton_solver_read_environment()

   if (do_read_sed .and. thereis_mol) then
      ! One read of the file for the three continuum bands; a band the run
      ! needs that the file does not reach stops the run and asks for the
      ! key (a stated value, zero included, is the run's own decision; a
      ! silent zero is not).
      call sed_band_fluxes(n_fuv_band, fuv_band_lo_A, fuv_band_hi_A,      &
                           band_flux_sed, band_covered_sed)
      if (.not. lw_flux_stated) then
         call band_from_spectrum('Stellar LW flux', ib_LW, F_LW_star,     &
                                 lw_from_spectrum)
      endif
      if (thereis_oxychem) then
         if (.not. fuv_b3_flux_stated)                                    &
            call band_from_spectrum('Stellar FUV B3 flux', ib_B3,         &
                                    F_FUV_B3, fuv_b3_from_spectrum)
         if (.not. fuv_b4_flux_stated)                                    &
            call band_from_spectrum('Stellar FUV B4 flux', ib_B4,         &
                                    F_FUV_B4, fuv_b4_from_spectrum)
      endif
   endif

   if (thereis_oxychem) then
      if (.not. thereis_mol) then
         write(*,*) '(input_read) ERROR: "Oxygen chemistry: True" needs'// &
                    ' "Molecular chemistry: True".'
         write(*,*) '  The oxygen cycle acts on H2: OH + H2 -> H2O + H is'
         write(*,*) '  95% of its net rate, and without the molecular'
         write(*,*) '  network there is no H2 to act on. Turn the'
         write(*,*) '  molecular chemistry on, or the oxygen chemistry off.'
         error stop 1
      endif
      if (.not. thereis_He) then
         write(*,*) '(input_read) ERROR: "Oxygen chemistry: True" needs'// &
                    ' He/H > 0 (it extends the molecular network, which'
         write(*,*) '  itself needs helium).'
         error stop 1
      endif
      if (melem_ab(iel_O) .le. 0.0d0) then
         write(*,*) '(input_read) ERROR: "Oxygen chemistry: True" but'//   &
                    ' the oxygen abundance is zero.'
         write(*,*) '  The option solves the partition of the oxygen'
         write(*,*) '  element among O I/II/III, OH, H2O and CO, so there'
         write(*,*) '  has to be an oxygen element. Set it in metals.inp'
         write(*,*) '  (X_O) or through the base.inp key "O_H_base", or'
         write(*,*) '  turn the oxygen chemistry off.'
         error stop 1
      endif
      if (q_h2_base .gt. 0.0d0 .and. lap_in_use) then
         ! A lower-atmosphere PROFILE is a different case from the scalar
         ! key, and section 4.6 of the design accepts it: the profile owns
         ! the region below the matching level, where its q_H2 sets the base
         ! particle count, and the chemistry owns the partition above it.
         ! They are not two owners of one quantity, but the reader has to be
         ! told that both are in play.
         write(*,'(A,F8.5,A)') ' (input_read) the lower-atmosphere'//      &
            ' profile states q_H2 =', q_h2_base, ' at the matching level;'
         write(*,'(A)') '   it sets the base particle count there, and'//  &
            ' the oxygen chemistry computes the'
         write(*,'(A)') '   partition above it.'
      else if (q_h2_base .gt. 0.0d0) then
         write(*,*) '(input_read) ERROR: base.inp key "q_H2_base" is'//    &
                    ' refused while "Oxygen chemistry: True".'
         write(*,*) '  q_H2_base is the imported base H2 fraction, and the'
         write(*,*) '  oxygen chemistry COMPUTES that partition; accepting'
         write(*,*) '  both would build the base particle count and the'
         write(*,*) '  base composition from two different H2 fractions,'
         write(*,*) '  and would make the comparison against the imported'
         write(*,*) '  value circular. Delete q_H2_base from base.inp, or'
         write(*,*) '  turn the oxygen chemistry off.'
         error stop 1
      endif
      ! The H/He element operator is now compatible with the oxygen
      ! chemistry and no longer refused: project_elements scales every
      ! H-bearing species by the same rH, and the metal ion columns with it,
      ! so an element that is partly in OH, H2O or CO keeps its ratio to
      ! hydrogen and the mass closure m_1 n_H + m_He n_He = rho still holds.
      !
      ! What is still refused is TRACE-METAL DIFFUSION. It transports each
      ! metal element's own mixing ratio, so it has to move the element's
      ! molecular carriers with the ion stages -- and CO carries one oxygen
      ! AND one carbon, so a step that moves O and C by different factors
      ! has no single factor to give it. That is the position HeH+ occupies
      ! for the H/He pair, where the code uses rBoth = min(rH, rHe); writing
      ! the same rule for a two-element molecule is its own piece of work
      ! and is not done here.
      if (he_metal_diffusion) then
         write(*,*) '(input_read) ERROR: "Oxygen chemistry: True" and'//   &
                    ' "He_metal_diffusion: True" are not solved together.'
         write(*,*) '  The trace-metal diffusion operator counts each metal'
         write(*,*) '  element over its ion stages alone, so with the'
         write(*,*) '  oxygen chemistry on it would transport an oxygen'
         write(*,*) '  reservoir missing everything bound into OH, H2O'
         write(*,*) '  and CO -- about half the element. Moving the'
         write(*,*) '  carriers with it needs a rule for CO, which'
         write(*,*) '  carries an oxygen AND a carbon nucleus and so'
         write(*,*) '  cannot follow two different element factors.'
         write(*,*) '  "He_diffusion: True" alone is accepted: the'
         write(*,*) '  molecular carriers already move with hydrogen'
         write(*,*) '  there. Run one or the other.'
         error stop 1
      endif
      ! A lower-atmosphere PROFILE is accepted beside the option: the
      ! profile owns the region below the matching level and the option owns
      ! the region above it, so they are not two owners of one quantity.
      ! Its solution_id pairing is checked as usual, elsewhere.
      write(*,'(A)') ' (input_read) Oxygen chemistry on: OH, H2O and CO'// &
         ' solved with the molecular network.'
      write(*,'(A,4ES10.3)') '   FUV band fluxes at the planet'//          &
         ' [erg cm^-2 s^-1] LW/B2(Lya)/B3/B4: ',                          &
         F_LW_star, F_Lya_star, F_FUV_B3, F_FUV_B4
      if (F_LW_star + F_Lya_star + F_FUV_B3 + F_FUV_B4                    &
          .le. 0.0d0) then
         write(*,'(A)') ' (input_read) WARNING: the oxygen chemistry is'// &
            ' on but every FUV band flux is zero, so there is no'
         write(*,'(A)') '   H2O or OH photolysis. The cycle then has no'// &
            ' way back from H2O and the partition it'
         write(*,'(A)') '   returns is the chemical equilibrium of the'//  &
            ' O/OH/H2O family, not a photochemical one.'
         write(*,'(A)') '   Set "Stellar LW flux", "Stellar FUV'//         &
            ' B3/B4 flux" and "Stellar Lya flux" for the'//                &
            ' photochemical result.'
      else if (F_LW_star .le. 0.0d0 .and.                                 &
               F_FUV_B3 .le. 0.0d0 .and. F_FUV_B4 .le. 0.0d0) then
         write(*,'(A)') ' (input_read) WARNING: no FUV continuum band'//   &
            ' flux is set (only Ly-alpha), so the run gets no oxygen'
         write(*,'(A)') '   photolysis outside the Ly-alpha line -- and'// &
            ' the continuum bands carry most of the'
         write(*,'(A)') '   H2O loss. Set "Stellar LW flux" and'//         &
            ' "Stellar FUV B3 flux" at least.'
      else if (F_LW_star .le. 0.0d0) then
         write(*,'(A)') ' (input_read) WARNING: "Stellar LW flux" is'//    &
            ' zero, so the run gets no H2O or OH photolysis over'
         write(*,'(A)') '   912-1201 A, where their cross sections'//      &
            ' peak. That band is entered through the'
         write(*,'(A)') '   Lyman-Werner key because H2 shares it.'
      endif
   else
      if (F_FUV_B3 + F_FUV_B4 .gt. 0.0d0)                                &
         write(*,'(A)') ' (input_read) WARNING: a "Stellar FUV B*'//       &
            ' flux" is set but "Oxygen chemistry" is off; there is'//      &
            ' no H2O or OH to photolyse, so the key has no effect.'
   endif

   ! The molecular infrared bands act on H2, H2O and CO. H2 is a solved
   ! species only with the molecular network on; H2O and CO only with the
   ! oxygen chemistry on. Report rather than stop -- the key is then inert on
   ! the species that are missing.
   if (mol_ir_bands) then
      if (.not. thereis_mol) then
         write(*,'(A)') ' (input_read) WARNING: "Molecular IR bands" is'//  &
            ' set but "Molecular chemistry" is off; there is no H2, H2O'//  &
            ' or CO, so the key has no effect.'
      else
         if (.not. thereis_oxychem)                                        &
            write(*,'(A)') ' (input_read) "Molecular IR bands": H2 lines'// &
               ' only -- H2O and CO exist only with "Oxygen chemistry".'
         if (.not. base_ir_field) then
            write(*,'(A)') ' (input_read) WARNING: "Molecular IR bands"'//  &
               ' is on but "Base IR field" is off, so the new bands emit'
            write(*,'(A)') '   into vacuum with no incident field. That'//  &
               ' is the configuration item (G) is ABOUT: it deepens the'
            write(*,'(A)') '   collapse of the molecular layer instead'//   &
               ' of holding it. Set "Base IR field: True" as well unless'
            write(*,'(A)') '   the emission-only limit is what you want.'
         endif
      endif
   endif

   ! Lyman-Werner photodissociation acts on H2, which only exists as a
   ! solved species with the molecular network on. Report rather than stop:
   ! the key is then simply inert.
   if (F_LW_star .gt. 0.0d0) then
      if (thereis_mol) then
         write(*,'(A,ES10.3,A)') ' (input_read) Lyman-Werner band flux at'//&
            ' the planet (912-1201 A; also the first oxygen photolysis'//   &
            ' band): ', F_LW_star, ' erg cm^-2 s^-1'
      else
         write(*,'(A)') ' (input_read) WARNING: "Stellar LW flux" is set'// &
            ' but "Molecular chemistry" is off; there is no H2 to'//        &
            ' photodissociate, so the key has no effect.'
      endif
   endif

   R0     = R0*RJ
   Mp     = Mp*MJ
   a_orb  = a_orb*AU
   Mstar  = Mstar*Msun
   Mrapp  = Mstar/Mp
   atilde = a_orb/R0
   if (spherical_domain) then
      ! Spherical mode: outer boundary set explicitly by the user [R_p];
      ! grid is normalized to R0 = R_p so r_max = r_out_user directly.
      if (r_out_user .le. 1.0d0) then
         write(*,*) '(input_read.f90) ERROR: Domain mode = Spherical '   // &
                    'requires "Outer radius [R_p]:" > 1.0 in input.inp.'
         error stop 1
      endif
      r_max = r_out_user
   else
      ! Roche mode: full tidal+centrifugal potential retained (grav_field
      ! else-branch). Default outer boundary is the Hill/L1 radius; an
      ! explicit "Outer radius [R_p]:" > 1 extends the domain PAST L1 while
      ! KEEPING the tidal potential (Yan 2022 / Huang 2023-style tidal +
      ! extended Rmax). Beyond L1 the Roche potential is past its maximum
      ! (net outward force), so the 1D radial flow there is a spherical
      ! approximation to the L1 funnel -- matching the CETIMB extended-domain
      ! treatment.
      r_max = (3.0*Mrapp)**(-1.0/3.0)*atilde
      if (r_out_user .gt. 1.0d0) r_max = r_out_user
   endif
            
	!------ Normalization constants ------!

   ! Solving molecular chemistry implies the molecular-base particle count
   ! for the base pressure BC (see the 'Molecular chemistry' key note above);
   ! an atomic ntot_bc with a molecular base is physically inconsistent.
   if (thereis_mol .and. .not. molecular_base) then
      molecular_base = .true.
      write(*,*) '(input_read.f90) Molecular chemistry on: enabling the'//   &
         ' molecular base particle count (Molecular base: True).'
   endif

   ! Composition factors. mass_per_H = gas mass per H nucleus [m_H];
   ! ntot_bc = total nuclei density at the base in units of n0 (n0 keeps
   ! its H/He-nuclei meaning, so n_H = n0/(1+HeH) is unchanged). With
   ! eos_metals 1 (default) the trace metals contribute their mass and
   ! their nuclei; with eos_metals 0 (or no metals) both reduce to the
   ! H/He-only values (mass_per_H = 1 + (m_He/m_H)*HeH with m_He/m_H =
   ! bsp_mass(He I) of species_table.f90, ntot_bc = 1).  With
   ! "Molecular base: True" comp_ntot_bc also removes the H nuclei bound
   ! into H2 at the base (passive molecular base, EOS-only: the species
   ! arrays stay atomic; docs/lower_atmosphere_coupling.*).
   ! Routed through the composition module (single source of the base
   ! composition policy); the eos_metals, metals-present and molecular_base
   ! branching lives inside it.
   ! ---- the base H2 fraction must be attainable at this He/H ---------- !
   ! q_H2 = n_H2/(n_H2+n_H+n_He) cannot exceed 0.5/(0.5+He/H), the value
   ! reached when every H nucleus is bound into H2. A larger number does not
   ! describe a more molecular gas; it describes a mixture that the element
   ! ratio cannot supply, so it is refused here rather than quietly reduced
   ! to the ceiling further down (which is what h2_bound_fraction used to do,
   ! leaving the run to march on a base composition nobody requested).
   !
   ! Checked on the resolved value, whichever source produced it, so the two
   ! sources obey one rule. Both can violate it: a handoff q_H2_base is an
   ! external number, and the Visscher/Koskinen fit is calibrated for
   ! solar-like composition, so its asymptote (0.8384) lies above the ceiling
   ! for any He/H >= 0.167.
   !
   ! The comparison carries a round-off tolerance only. A handoff written
   ! from a fully molecular lower-atmosphere column sits just under its own
   ! ceiling by construction -- LHS1140b/examples/scalar_base_cno passes at
   ! 99.88% of it -- so anything looser would refuse a legitimate handoff.
   if (molecular_base .or. q_h2_base .gt. 0.0d0) then
      q_h2_resolved = h2_mixing_ratio_base()
      q_h2_max      = h2_mixing_ratio_ceiling()
      if (q_h2_resolved .gt. q_h2_max*(1.0d0 + 1.0d-12)) then
         write(*,*) '(input_read) ERROR: the base H2 fraction is not'//     &
                    ' attainable at this He/H.'
         if (q_h2_base .gt. 0.0d0) then
            write(*,*) '  source: lower-atmosphere handoff'//               &
                       ' (q_H2_base)'
         else
            write(*,*) '  source: chemical-equilibrium fit at'//            &
                       ' (p_base, T_base)'
         endif
         write(*,'(A,F10.6)') '   requested q_H2 = ', q_h2_resolved
         write(*,'(A,F10.6)') '   ceiling  q_H2  = ', q_h2_max
         write(*,'(A,F10.6)') '   He/H           = ', HeH
         write(*,*) '  The ceiling is 0.5/(0.5 + He/H): every H nucleus'//  &
                    ' bound into H2.'
         if (q_h2_base .gt. 0.0d0) then
            write(*,*) '  Fix the handoff: q_H2_base and HeH_base must'//   &
                       ' come from the same'
            write(*,*) '  lower-atmosphere solution, and this pair'//       &
                       ' cannot.'
         else
            write(*,*) '  The fit is calibrated for solar-like'//           &
                       ' composition and is out of range'
            write(*,*) '  here. Supply q_H2_base from a lower-atmosphere'// &
                       ' solution instead.'
         endif
         error stop 1
      endif
   endif

   mass_per_H = comp_mass_per_H()
   ntot_bc    = comp_ntot_bc()
   rho_bc     = comp_rho_bc()

   ! Echo which H2 source set the base particle count: the photochemical
   ! handoff (base.inp key q_H2_base) or the chemical-equilibrium fit.
   if (molecular_base) then
      if (q_h2_base .gt. 0.0d0) then
         write(*,'(A,F6.3,A,F6.3)') ' (input_read) Molecular base: '//     &
            'q_H2(base, photochemical) =', h2_mixing_ratio_base(),         &
            ' -> ntot_bc =', ntot_bc
      else
         write(*,'(A,F6.3,A,F6.3)') ' (input_read) Molecular base: '//     &
            'q_H2(base, chem.eq. fit)  =', h2_mixing_ratio_base(),         &
            ' -> ntot_bc =', ntot_bc
      endif
   endif

   ! ---- an EOS-only molecular base describes gas the species state does not -!
   ! "Molecular base: True" removes the H nuclei bound into H2 from the base
   ! particle count (comp_ntot_bc above). Without the molecular network the
   ! species arrays stay atomic, and the ntot_bc recomputation that reconciles
   ! the two (end of ioniz_eq) is gated on the network being solved, so it
   ! never fires. The base ghost pressure is then (ntot_bc + dp_bc)*T0 at the
   ! pinned density while the gas in it counts one particle per nucleus, and
   ! the isothermal base boundary silently sits at ntot_bc*T0 instead of T0.
   !
   ! This refusal was asked for in exactly these terms ("If the intended
   ! quantity remains EOS-only, the code must instead refuse a mismatch
   ! between that EOS state and the species state used at the same ghost");
   ! it is item P35, decided 2026-09-02 and recorded in section 120 of
   ! docs/Update_EXHALE_stage1.pdf. Measured before the refusal existed: the
   ! HD 209458 b VULCAN handoff runs ran at ntot_bc = 0.555-0.994, i.e. base
   ! ghosts from 0.555 to 0.994 of the temperature their input.inp asked for.
   if (molecular_base .and. .not. thereis_mol) then
      write(*,*) '(input_read) ERROR: a molecular base particle count with'// &
                 ' an atomic species state.'
      write(*,*) '  Molecular base: True removes the H nuclei bound into'//   &
                 ' H2 from the base'
      write(*,*) '  particle count, but Molecular chemistry is off, so'//     &
                 ' there are no H2 species'
      write(*,*) '  to hold them and the base ghost carries one particle'//   &
                 ' per nucleus. The'
      write(*,*) '  isothermal base boundary is then not isothermal:'
      write(*,'(A,F10.6)') '   base H2 mixing ratio q_H2 = ',                 &
                           h2_mixing_ratio_base()
      write(*,'(A,F10.6)') '   H nuclei bound into H2  x2 = ',                &
                           base_h2_nuclei_fraction()
      write(*,'(A,F10.6)') '   base particle count ntot_bc = ', ntot_bc
      write(*,'(A,F10.6,A,F10.3,A)') '   ghost temperature = ', ntot_bc,      &
                           ' x T0 = ', ntot_bc*T0, ' K'
      write(*,'(A,F10.3,A)') '   requested base temperature T0 = ', T0, ' K'
      write(*,*) '  Two ways to fix it, and which one is right depends on'//  &
                 ' what the base is:'
      write(*,*) '   (1) the base IS molecular -- set'//                      &
                 ' "Molecular chemistry: True" so the'
      write(*,*) '       network carries the H2 the particle count'//         &
                 ' already assumes;'
      write(*,*) '   (2) the base is ATOMIC -- remove'//                      &
                 ' "Molecular base: True" from input.inp'
      write(*,*) '       and q_H2_base from base.inp, so neither the'//       &
                 ' particle count nor the'
      write(*,*) '       ghost composition claims molecules.'
      error stop 1
   endif

   ! ---- forcing the base reservoir particle count, a diagnostic ------ !
   ! EXHALE_BASE_PARTICLE_COUNT = <positive real> | atomic replaces
   ! ntot_bc, the free particles the base reservoir holds per (H+He)
   ! nucleus, by the stated count. It exists to separate the two statements
   ! a molecular base makes at once, which comp_ntot_bc above combines:
   !
   !   (a) the H nuclei bound into H2 are not free particles, so at the
   !       prescribed base pressure and temperature the reservoir holds
   !       1/ntot_bc times as many nuclei, and n0 = p_base/(k_B T0
   !       ntot_bc) below carries that into the base mass density;
   !   (b) the lower ghost carries the same H2 in its composition.
   !
   ! The value "atomic" gives the count of the SAME composition with no H2
   ! binding, ntot_bc + h2_bound_fraction(), which is what comp_ntot_bc
   ! returns with molecular_base off. The base then has the mass density an
   ! atomic reservoir would have at that (p_base, T0) while (b) is
   ! untouched, so the face condition can be read against one changed
   ! quantity. Anything else is read as the count itself.
   !
   ! Off by default and in every production run: unset, empty, unreadable
   ! or not a positive number leaves ntot_bc exactly as comp_ntot_bc
   ! returned it, writes no line and touches nothing. When it is set the
   ! base state is not the one this input.inp describes, so the two counts
   ! are announced together and the run says what it is.
   ntot_bc_of_composition = ntot_bc
   call get_environment_variable('EXHALE_BASE_PARTICLE_COUNT',            &
                                 ntot_bc_forced_env)
   if (len_trim(ntot_bc_forced_env) .gt. 0) then
      if (trim(adjustl(ntot_bc_forced_env)) .eq. 'atomic') then
         ntot_bc_forced = ntot_bc
         if (molecular_base)                                             &
            ntot_bc_forced = ntot_bc + h2_bound_fraction()
      else
         read(ntot_bc_forced_env, *, iostat=ios_ntot_bc) ntot_bc_forced
         if (ios_ntot_bc .ne. 0) ntot_bc_forced = -1.0d0
      endif
      if (.not. (ntot_bc_forced .gt. 0.0d0)) then
         write(*,*) '(input_read) ERROR: EXHALE_BASE_PARTICLE_COUNT ='//  &
                    ' "'//trim(ntot_bc_forced_env)//'" is neither a'//    &
                    ' positive number nor "atomic".'
         write(*,*) '  It states the free particles per (H+He) nucleus'// &
                    ' the base reservoir holds.'
         error stop 1
      endif
      ntot_bc = ntot_bc_forced
      write(*,'(A)') ' (input_read) EXHALE_BASE_PARTICLE_COUNT: the'//    &
           ' base reservoir particle count is FORCED; this run is not'//  &
           ' the configuration its input.inp describes.'
      write(*,'(A,ES23.15E3)') '   ntot_bc from the base composition = ', &
           ntot_bc_of_composition
      write(*,'(A,ES23.15E3)') '   ntot_bc forced                    = ', &
           ntot_bc
   endif

   ! ---- the base level has one source -------------------------------- !
   ! A lower-atmosphere handoff is written AT one pressure: its composition,
   ! its q_H2 and its temperature all refer to that level, so the level is
   ! the handoff's to state and n0 follows,
   !
   !     n0 = p_base/(k_B T0 ntot_bc),
   !
   ! the same relation "Base BC: pressure" uses (1 bar = 1e6 erg/cm^3, and
   ! the small electron term dp_bc is added later by the ghost BC).  The
   ! legacy density key then says the same thing twice, and the two can
   ! disagree -- which is what put the hot-Uranus gate at 9 microbar while
   ! its handoff and its base radius both said 1.  So: if both are given they
   ! must agree, and if
   ! they do not the run stops here rather than marching on a base level
   ! nobody chose.  1% is the tolerance; a handoff and a density key that
   ! describe one level agree far better than that.
   !
   ! Both handoff forms are one rule: a base.inp carrying "p_base", and a
   ! lower-atmosphere PROFILE, whose matching pressure p_match is the level
   ! apply_lower_atmosphere_profile read T0, R0, q_H2 and every elemental
   ! ratio at (rev 3 section 10.2 item 4). A base placed anywhere else is a
   ! base whose state was taken from a different level.
   if (base_level_from_handoff) then
      p_base_cgs     = p_base_bar*1.0d6
      n0_from_p_base = p_base_cgs/(kb_erg*T0*ntot_bc)
      if (base_bc_mode .eq. 1) then
         if (abs(base_p_ubar/(p_base_bar*1.0d6) - 1.0d0) .gt. 1.0d-6) then
            write(*,*) '(input_read) ERROR: two different base levels were'//&
                       ' given.'
            write(*,'(A,A,A,ES12.4,A)') '   ', trim(base_level_source),    &
               '  = ', p_base_bar, ' bar'
            write(*,'(A,ES12.4,A)') '   "Base BC: pressure"   = ',          &
               base_p_ubar*1.0d-6, ' bar'
            write(*,*) '  The lower boundary is one level, so one input'//  &
                       ' states it. Fix by deleting'
            write(*,*) '  the "Base BC: pressure" line from '//             &
               trim(inp_file)//' (the handoff already'
            write(*,*) '  fixes the level), or by making the two numbers'// &
                       ' the same.'
            error stop 1
         endif
      endif
      if (base_density_key_given) then
         p_from_density_key = n0*kb_erg*T0*ntot_bc
         if (abs(p_from_density_key/p_base_cgs - 1.0d0) .gt. 1.0d-2) then
            write(*,*) '(input_read) ERROR: the base level stated by the'// &
                       ' lower-atmosphere handoff and'
            write(*,*) '  the one implied by "Log10 lower boundary number'//&
                       ' density" disagree by more than 1%.'
            write(*,'(A,A,A,ES12.4,A)') '   level from ',                   &
               trim(base_level_source), ' = ', p_base_bar, ' bar'
            write(*,'(A,ES12.4,A)') '   density key implies p          = ', &
               p_from_density_key*1.0d-6, ' bar'
            write(*,'(A,F10.4)')    '   ratio (density key / handoff)  = ', &
               p_from_density_key/p_base_cgs
            write(*,'(A,ES12.4,A)') '   n0 given                       = ', &
               n0, ' cm^-3'
            write(*,'(A,ES12.4,A)') '   n0 implied by the handoff      = ', &
               n0_from_p_base, ' cm^-3'
            write(*,'(A,F10.6)')    '   base particle count ntot_bc    = ', &
               ntot_bc
            write(*,*) '  The lower boundary is one level and these two'//  &
                       ' put it in different places.'
            write(*,*) '  Fix by deleting ONE of them:'
            write(*,*) '   - delete "Log10 lower boundary number'//         &
                       ' density" from '//trim(inp_file)//','
            write(*,*) '     and the handoff level fixes the base (the'//   &
                       ' usual choice: the handoff'
            write(*,*) '     composition refers to that level);'
            write(*,*) '   - or detach the handoff ("p_base" in base.inp,'//&
                       ' or the "Lower atmosphere'
            write(*,*) '     profile:" key in '//trim(inp_file)//'), and'// &
                       ' the density key fixes the level,'
            write(*,*) '     with the handoff composition then referring'//&
                       ' to whatever level that is.'
            error stop 1
         endif
      endif
      n0 = n0_from_p_base
      write(*,'(A,A,A,ES12.4,A,ES12.4,A)') ' (input_read) base level from ',&
         trim(base_level_source), ' =', p_base_bar, ' bar -> n0 =', n0,     &
         ' cm^-3'
   endif

   ! Pressure-anchored base (Base BC: pressure): override n0 so that the base
   ! pressure n0*kb*T0*ntot_bc matches the target base_p_ubar [microbar].
   ! 1 microbar = 1 erg/cm^3. (dp_bc, the small electron term, is negligible
   ! here and is added by the ghost BC later.) This is the CETIMB 1-microbar
   ! lower boundary: a much less dense base, hence a weaker rho*g source.
   if (base_bc_mode .eq. 1) then
      n0 = base_p_ubar/(kb_erg*T0*ntot_bc)
      if (.not. base_level_from_handoff)                                  &
         base_level_source = '"Base BC: pressure"'
      write(*,'(A,ES12.4,A)') ' (input_read) Base BC pressure mode: '//   &
         'derived n0 =', n0, ' cm^-3'
   endif

   ! Neither source stated the base level.
   if (.not. (n0 .gt. 0.0d0)) then
      write(*,*) '(input_read) ERROR: no base level was given.'
      write(*,*) '  Supply exactly one of:'
      write(*,*) '   - "Log10 lower boundary number density [cm^-3]:'//    &
                 ' <log10 n0>" in '//trim(inp_file)//','
      write(*,*) '   - "p_base <bar>" in base.inp (a lower-atmosphere'//   &
                 ' handoff states its own level),'
      write(*,*) '   - "Base BC: pressure <microbar>" in '//               &
                 trim(inp_file)//'.'
      error stop 1
   endif

   v0     = sqrt(kb_erg*T0/mu)
   t_s    = R0/v0
   p0     = n0*mu*v0*v0
   q0     = n0*mu*v0*v0*v0/R0
   b0     = (Gc*Mp*mu)/(kb_erg*T0*R0)
   dp_bc  = 1.0e-10
   ! Cell-1 particle count seen by the lower boundary before the first
   ! composition solve refreshes it (Apply_BC can run first, e.g. in the
   ! IC/residual paths). The base value makes cell 1 start at the base
   ! composition, which is the only statement available before a solve.
   n_part_cell1 = ntot_bc + dp_bc

   ! ----- The lower-boundary reservoir -----
   ! (p, s) at the base LEVEL r = 1, carried as the isentrope through
   ! (p = ntot_bc + dp_bc, T = T0) at the base composition. base_boundary
   ! continues it to the first face itself, so the level the user states does
   ! not move when the grid does. "Base BC: pressure" has already set n0 so
   ! that this pressure is the requested one; "Base BC: density" states n0 and
   ! the pressure follows.
   call set_base_reservoir(ntot_bc + dp_bc, 1.0d0,                        &
                           (ntot_bc + dp_bc)/rho_bc, 1.0d0)

   !------ Allocations ------!

   ! Grid-sized module arrays of global_parameters. N is final here (the
   ! optional "Grid cells:" key was resolved in the keyword block above) and
   ! nothing before this point touches the radial grid.
   call allocate_grid_arrays

   ! Allocate variables according to composition
   if (.not.thereis_He) then
      	N_eq = 1
   else
		if (thereis_HeITR) then
			N_eq = 4
		else
			N_eq = 3
		endif
		! If metals are present, the system grows by 2*n_melem variables:
		! HII, HeII, HeIII (+ HeITR if the triplet is on), then two fractions
		! (X+, X++) per metal element. Absent elements are force-zeroed.
		! HeITR + metals is solved by the merged System_HeH_TR_metals: the
		! triplet keeps x(4) and the metals shift to x(5+2*(e-1)).
		if (thereis_metals .and. .not.thereis_HeITR) N_eq = 3 + 2*n_melem
		if (thereis_metals .and.      thereis_HeITR) N_eq = 4 + 2*n_melem

	! molecular system: H+/He+/He++ + H2/H2+/H3+/HeH+ (+ He 2^3S), and, when
	! metals are also present, two more unknowns per metal element appended
	! above them (System_HeH_mol_metals; metal_row_base = 8 or 9).
	if (thereis_mol) then
		N_eq = 7
		if (thereis_HeITR) N_eq = 8
		! Two more unknowns for the oxygen carriers OH and H2O, between the
		! molecular block and the metals (System_HeH_mol_metals:
		! oxygen_row_base / metal_row_base).
		if (thereis_oxychem) N_eq = N_eq + 2
		if (thereis_metals) N_eq = N_eq + 2*n_melem
	endif
	endif
	
   lwa  = (N_eq*(3*N_eq+13))/2
   allocate (sys_sol(N_eq))
   allocate (sys_x(N_eq))
   allocate (wa(lwa))

   ! Build the active charge-exchange reaction set (Huang Table 4). cx_full
   ! was set by read_metals_input; the default is the metal-H group only.
   if (thereis_metals) call cx_init

   ! Fill the FUV photolysis threshold energies from the thermodynamic
   ! table. Done here, serially, because the cell sweep that reads them runs
   ! OpenMP-parallel and must find them already written.
   if (thereis_oxychem) call water_photolysis_init

   ! Whether the carriers are transported, when the run did not say.
   !
   ! The oxygen cycle brings its own answer: it exists to compute a base
   ! partition that a local steady state cannot produce, so it turns the
   ! transport on. A molecular run without it keeps the local-equilibrium
   ! closure it has always had; that closure is violated on its own solution
   ! across the H2 front, and
   ! changing this default is a deliberate step with a golden refresh, not a
   ! side effect of adding the operator.
   if (.not. carrier_transport_stated) carrier_transport = thereis_oxychem

   ! WHAT "Ionization transport: True" REQUIRES, refused rather than repaired.
   ! Each of these is a configuration in which the option would silently
   ! mean something other than what it says.
   if (ionization_transport) then
      ! WHAT THE KEY CARRIES.  The ionized stages of an element as
      ! fractions of that element's NUCLEI, transported on that element's
      ! own nucleus face flux (ionization_stage_transport, equation 1).
      ! The carried set is x(H II) per hydrogen nucleus and x(He II),
      ! x(He III) per helium nucleus; the neutral stage of each element
      ! closes its simplex and is not a row.
      if (.not. thereis_He) then
         write(*,*) '(input_read) ERROR: "Ionization transport: True" needs'
         write(*,*) '  helium in the mixture. The flux a stage rides on is'
         write(*,*) '  the element nucleus flux of the binary H/He'
         write(*,*) '  transport operator, which returns nothing for a'
         write(*,*) '  mixture with one element, and the carried set'
         write(*,*) '  includes the two ionized stages of helium.'
         write(*,*) '  Aborting.'
         error stop 1
      endif
      ! AN ATOMIC GAS CARRIES THEM TOO.  The stages are stages of an
      ! element, not molecular carriers: their unknowns are fractions per
      ! element nucleus, their flux is the element's own nucleus face flux,
      ! and their source in an atomic gas is the H/He ionization balance the
      ! local sweep of that gas solves (ion_residual_core, through
      ! diffusive_photochemistry's carrier_source).  The entry points of the
      ! transport-chemistry operator, the frozen background its rows are
      ! evaluated on and the alternation of the stationary outer iteration
      ! are all conditioned on "a transported row exists"
      ! (ionization_equilibrium, transported_rows_exist), so the molecular
      ! network is no longer what lets them run.
      if (thereis_mol .and. .not. carrier_transport) then
         write(*,*) '(input_read) ERROR: "Ionization transport: True" in a'
         write(*,*) '  molecular gas ("Molecular chemistry: True") needs'
         write(*,*) '  "Molecular carrier transport: True". The stage'
         write(*,*) '  sources are then rows (1), (2) and (3) of the'
         write(*,*) '  molecular H/He network, whose molecular sinks are'
         write(*,*) '  written at the molecular densities of the same'
         write(*,*) '  transported composition; with the molecular'
         write(*,*) '  carriers left on their local root the two halves of'
         write(*,*) '  that network would describe different states of the'
         write(*,*) '  same cell. Set both keys, or neither. Aborting.'
         error stop 1
      endif
      ! METALS IN THE MIXTURE ARE SUPPORTED. The transported stage sources
      ! carry the metal charge exchange of Huang et al. (2023) Table 4 --
      ! group A, metal + H and H+, active whenever metals are present;
      ! group C, metal + He and He+, under "cx_full 1"; and the group E
      ! electron capture O2+ + H0 -> O+ + H+ under its own scale -- through
      ! charge_exchange::charge_exchange_stage_sources, the same reaction
      ! set and the same rate coefficients the local sweep's rows take
      ! through cx_add_to_fvec. The stage the flow carries and the
      ! composition the sweep returns therefore answer one H+ and He+
      ! balance of a cell, which is what the refusal that stood here
      ! protected while the terms were missing.
      if (carrier_in_newton) then
         write(*,*) '(input_read) ERROR: "Ionization transport: True" with'
         write(*,*) '  "Coupled carrier solve: True" is refused. The'
         write(*,*) '  coupled row registry carries every carrier unknown'
         write(*,*) '  as the species mass fraction of its f_sp column,'
         write(*,*) '  and a stage unknown is a fraction per element'
         write(*,*) '  nucleus; the two are different variables and the'
         write(*,*) '  bound the coupled solve puts on the unknown would'
         write(*,*) '  be the bound of the other one. Use the partitioned'
         write(*,*) '  stationary route. Aborting.'
         error stop 1
      endif
      ! THE STATIONARY SOLVE AND THE TRANSPORTED PROTON.
      !
      ! The stationary route is the partitioned
      ! alternation (steady_wind_with_element_diffusion, item P6 of
      ! 2026-09-11): the hydrodynamic solve holds the composition, and
      ! every equilibrium sweep it makes is handed the transported stage
      ! fractions of that composition (ionization_equilibrium, the imposed
      ! block: written whenever the key is on and the state is a restart or
      ! carries a background), so no stage is put back on its local root;
      ! the carrier relaxation between the hydrodynamic solves then
      ! transports them. Marching is unaffected.
      if (use_newton_solver)                                              &
         write(*,'(A)') ' (input_read) Ionization transport with the'//    &
            ' partitioned stationary route: the sweeps of the'//          &
            ' hydrodynamic solve are handed the transported x(H II),'//   &
            ' x(He II) and x(He III); the carrier relaxation'//           &
            ' transports them.'
   endif

   ! A COUPLED STEADY SOLVE MAY NOT ELIMINATE H2.
   !
   ! With the carriers eliminated, n(H2) is not an unknown of the stationary
   ! system: it is whatever the local-equilibrium sweep returns for the
   ! composition it was seeded with. In the shielded layer the sweep cannot
   ! return a content at all. The fast chemistry there cycles
   ! H2 -> H2+ -> H3+ -> H2 without changing the total number of H2 nuclei,
   ! so the local balance rows fix only the PARTITION among the molecular
   ! species and leave their sum where the seed put it; the content is set by
   ! the slow formation and dissociation and by transport, neither of which a
   ! local equilibrium sees. MEASURED on the hot Uranus element state: 0.77
   ! of any seed perturbation
   ! of the layer's H2 content survives every pass of the sweep, the same 0.77
   ! at perturbations of 1e-6 and of 1e-2, and the base cell's energy row,
   ! being a near-cancellation of the fluxes its continuous-temperature ghost
   ! produces, amplifies that by ~4e3. Two evaluations of ONE state then
   ! differ by 3.2e5 times "Resid tol" per unit relative seed change. A
   ! residual that is not a function of its unknowns has no root, and a state
   ! accepted at "Resid tol" on it is accepting its seed.
   !
   ! With "Molecular carrier transport: True" the H2 row of the network
   ! becomes x - x_fix with x_fix a Newton unknown, so the content is solved
   ! from its transport balance instead of eliminated, and the same
   ! measurement reads 9.7e-10, below "Resid tol" = 1e-8. That is the
   ! configuration this refusal names. The marching path is not affected:
   ! it never eliminates a quantity it also has to determine.
   ! "On stall" reaches the same block, one pass later, so the same
   ! statement holds for it.
   if ((carrier_in_newton .or. carrier_newton_on_stall) .and. thereis_mol &
       .and. .not. carrier_transport) then
      write(*,*) '(input_read) ERROR: "Coupled carrier solve: True" on a'
      write(*,*) '  (and "Coupled carrier solve: On stall", which'
      write(*,*) '  reaches the same block one pass later) on a'
      write(*,*) '  molecular configuration ("Molecular chemistry: True")'
      write(*,*) '  needs the molecular carriers transported. With them'
      write(*,*) '  eliminated, n(H2) is not an unknown of the stationary'
      write(*,*) '  system, and the local-equilibrium elimination does not'
      write(*,*) '  determine the molecular hydrogen content of the'
      write(*,*) '  shielded layer: the fast chemistry conserves H2 nuclei'
      write(*,*) '  there, so the sweep returns the content it was seeded'
      write(*,*) '  with. The stationary residual is then not a function of'
      write(*,*) '  its unknowns -- MEASURED, its seed dependence is 3.2e5'
      write(*,*) '  times "Resid tol" per unit relative seed change -- and'
      write(*,*) '  there is no root for the Newton to converge to.'
      write(*,*) '  Set "Molecular carrier transport: True", which makes'
      write(*,*) '  n(H2) a Newton unknown solved from its own transport'
      write(*,*) '  balance, or drop "Coupled carrier solve" and march.'
      write(*,*) '  Aborting.'
      error stop 1
   endif

   ! A CARRIER TRANSPORT WITH NOTHING TO TRANSPORT.
   !
   ! The transported carriers are species of the molecular network, and
   ! every consumer of carrier_transport is guarded by thereis_mol, so with
   ! the network off the key changes nothing: the run is the atomic one it
   ! would have been without the line. That is why this is reported and not
   ! refused -- what is wrong is the input file, not the state it produces.
   ! Said rather than left silent, because a key that is present and inert
   ! is the kind of thing a reader of the input file has no way to notice
   ! (the same reason "Molecular IR bands" and "Stellar LW flux" say it).
   !
   ! Only an explicit "Molecular carrier transport: True" reaches here with
   ! the network off: the unstated default is the oxygen chemistry, and that
   ! is refused above without the network.
   if (carrier_transport .and. .not. thereis_mol) then
      write(*,'(A)') ' (input_read) WARNING: "Molecular carrier'//        &
         ' transport: True" is set but "Molecular chemistry" is off;'
      write(*,'(A)') '   the transported carriers are species of the'//   &
         ' molecular network, so the run carries none of them and'
      write(*,'(A)') '   the key has no effect. Set "Molecular'//         &
         ' chemistry: True" to transport H2, or remove the line.'
   endif

   ! THE RUN MODE, AND WHY ITS DEFAULT IS init.
   !
   ! The ordinary use of this code is a stationary solution: march until the
   ! flux functional is flat, finish with a stationary solve, certify the
   ! state. Every step of that is a relaxation iterate and none of it is a
   ! trajectory, so in the contract's own terms the ordinary run is
   ! initialization / continuation. A physical time integration is the
   ! exception and is asked for; nothing is inferred, and a run that says
   ! nothing claims no elapsed time.
   !
   ! A physical integration additionally needs one global dt per step,
   ! because cells advanced by different intervals do not form one
   ! trajectory: in a closed two-cell system with one internal flux, unequal
   ! steps leave the material sum changed, so a local-dt update is a
   ! relaxation iterate however long it is run. That is refused rather than
   ! reinterpreted.
   if (run_mode .eq. run_mode_phys .and. use_local_dt) then
      write(*,*) '(input_read) ERROR: "Run mode: phys" and'
      write(*,*) '  "Time stepping: Local" cannot both be set. Local'
      write(*,*) '  pseudo-time advances each cell by its own interval, so'
      write(*,*) '  the material sum across an internal face is not'
      write(*,*) '  conserved and the update is a relaxation iterate rather'
      write(*,*) '  than one physical step; there is no elapsed time to'
      write(*,*) '  report for it. Use "Run mode: init" with local time'
      write(*,*) '  stepping, or a global dt with "Run mode: phys".'
      write(*,*) '  Aborting.'
      error stop 1
   endif

   ! WHAT THIS RESTART CONTINUES, AND THE THREE STATEMENTS THAT CONTRADICT
   ! THEMSELVES.
   !
   ! The intent is a statement about a state read from a file, so a run that
   ! reads none has nothing to state; a trajectory is a physical integration,
   ! so it cannot be continued by a run that claims no elapsed time; and a
   ! stationary intent enters the stationary solve, so it needs one to enter.
   ! Each is refused rather than reinterpreted: a key that is silently
   ! ignored is a run doing something other than what its input file says.
   if (restart_intent_given .and. .not. do_load_IC) then
      write(*,*) '(input_read) ERROR: "Restart intent" is set but'
      write(*,*) '  "Load IC?" is False, so there is no loaded state for'
      write(*,*) '  the intent to be about. Set "Load IC? True", or'
      write(*,*) '  remove the line. Aborting.'
      error stop 1
   endif
   ! The same statement about the same absent state: naming the options a
   ! restart may change says nothing when there is no state being restarted.
   if (restart_option_change_given .and. .not. do_load_IC) then
      write(*,*) '(input_read) ERROR: "Restart option change" is set but'
      write(*,*) '  "Load IC?" is False, so no state is being restarted'
      write(*,*) '  and there is no option set to compare with. Set'
      write(*,*) '  "Load IC? True", or remove the line. Aborting.'
      error stop 1
   endif
   if (restart_intent_given .and.                                          &
       restart_intent .eq. restart_intent_trajectory .and.                 &
       run_mode .ne. run_mode_phys) then
      write(*,*) '(input_read) ERROR: "Restart intent: trajectory" with'
      write(*,*) '  "Run mode: init". A trajectory is a physical'
      write(*,*) '  integration: its clock, its accepted-step budgets and'
      write(*,*) '  its histories exist only in physical mode, and an'
      write(*,*) '  initialization run claims no elapsed time at all.'
      write(*,*) '  Set "Run mode: phys" to continue the trajectory, or'
      write(*,*) '  "Restart intent: relaxation". Aborting.'
      error stop 1
   endif
   if (restart_intent .eq. restart_intent_stationary .and.                 &
       .not. use_newton_solver) then
      write(*,*) '(input_read) ERROR: "Restart intent: stationary" needs'
      write(*,*) '  "Solver: Newton". The intent evaluates the stationary'
      write(*,*) '  residual of the loaded state and then enters the'
      write(*,*) '  stationary solve; with no stationary solver selected'
      write(*,*) '  there is nothing to enter. Set "Solver: Newton", or'
      write(*,*) '  "Restart intent: relaxation" to march. Aborting.'
      error stop 1
   endif
   ! THE DEFAULT IS THE RUN MODE'S OWN MEANING. A physical run continues the
   ! trajectory the file records; an initialization or continuation run
   ! continues relaxing, which is what every restart did before this key
   ! existed. No configuration changes behavior by the key being absent.
   if (.not. restart_intent_given) then
      if (run_mode .eq. run_mode_phys) then
         restart_intent = restart_intent_trajectory
      else
         restart_intent = restart_intent_relaxation
      endif
   endif
   if (do_load_IC) then
      select case (restart_intent)
      case (restart_intent_trajectory)
         write(*,'(A)') ' (input_read) Restart intent: trajectory'//       &
              ' (the clock continues from the state file)'
      case (restart_intent_stationary)
         if (stationary_evaluate_only) then
            write(*,'(A)') ' (input_read) Restart intent: stationary'//    &
                 ' evaluate (measure the state as loaded, take no step)'
         else if (stationary_equilibrate_loaded) then
            write(*,'(A)') ' (input_read) Restart intent: stationary'//    &
                 ' equilibrate (put the loaded composition on its own'//   &
                 ' fixed point first)'
         else
            write(*,'(A)') ' (input_read) Restart intent: stationary'//    &
                 ' (measure the state as loaded, then solve; no CFL step)'
         endif
      case default
         write(*,'(A)') ' (input_read) Restart intent: relaxation'//       &
              ' (march toward stationarity, then the solver hand-off)'
      end select
      if (restart_option_change_given) then
         optline = ' '
         do iopt = 1, n_opt
            if (restart_option_change_named(iopt))                         &
               optline = trim(optline)//' '//trim(opt_name(iopt))
         enddo
         write(*,'(A)') ' (input_read) Restart option change:'//           &
              trim(optline)//' (these options may differ between the'//    &
              ' state and this run; every other'
         write(*,'(A)') '              difference still refuses the load)'
      endif
   endif

   ! Fix the transported molecular carrier set (H2 alone, or H2 with the
   ! three oxygen carriers, plus H+ when the ionization state is carried).
   ! One owner, evaluated after every key is parsed.
   call carrier_set_init

   ! End of subroutine

   contains

   function find_lbl(key, found) result(res)
   ! Last line matching `key` (last occurrence wins, as in the keyword loop).
   character(len=*), intent(in)  :: key
   logical,          intent(out) :: found
   character(len=250) :: res
   integer :: j
   found = .false.
   res   = ''
   do j = 1, nlines
      if (lbl_match(filelines(j), key)) then
         res   = filelines(j)
         found = .true.
      endif
   enddo
   end function find_lbl

   function req(key) result(res)
   ! Mandatory key: return its (last) line or abort naming the key.
   character(len=*), intent(in) :: key
   character(len=250) :: res
   logical :: found
   res = find_lbl(key, found)
   if (.not. found) then
      write(*,*) '(input_read.f90) ERROR: mandatory input line "'//   &
         trim(key)//'" not found in '//trim(inp_file)//'.'
      error stop 1
   endif
   end function req

   function req_eband() result(res)
   ! Mandatory energy-band line "[E_low,E_mid(,E_high)] = [ ... ]", anchored
   ! by its bracket-delimited "[E_low" prefix (the char after the prefix is
   ! ',', so this uses a plain prefix test rather than lbl_match).
   character(len=250) :: res, t
   integer :: j
   logical :: found
   found = .false.
   res   = ''
   do j = 1, nlines
      t = adjustl(filelines(j))
      if (len_trim(t) .ge. 6) then
         if (t(1:6) .eq. '[E_low') then
            res   = filelines(j)
            found = .true.
         endif
      endif
   enddo
   if (.not. found) then
      write(*,*) '(input_read.f90) ERROR: mandatory energy-band line'//   &
         ' "[E_low,E_mid,E_high] = [ ... ]" not found in '//              &
         trim(inp_file)//'.'
      error stop 1
   endif
   end function req_eband

   subroutine band_from_spectrum(keyname, ib, F, from_spectrum)
   ! One band of the loaded spectrum handed to its key: the integral when
   ! the file covers the band, a stop naming the key when it does not.
   character(len=*), intent(in)  :: keyname
   integer,          intent(in)  :: ib
   real*8,           intent(out) :: F
   logical,          intent(out) :: from_spectrum
   if (.not. band_covered_sed(ib)) then
      write(*,'(A,A,A)') ' (input_read) ERROR: "', keyname, '" is not'//   &
           ' stated and the spectrum file does not reach the band'
      write(*,'(A,F7.1,A,F7.1,A)') '   ', fuv_band_lo_A(ib), ' - ',       &
           fuv_band_hi_A(ib), ' A it would be integrated over. State'//   &
           ' the key (a stated 0 switches the band off), or supply a'//   &
           ' file that covers the band. Aborting.'
      error stop 1
   endif
   F = band_flux_sed(ib)
   from_spectrum = .true.
   write(*,'(A,A,A,ES10.3,A,F7.1,A,F7.1,A)') ' (input_read) ', keyname,   &
        ' from the spectrum file: ', F, ' erg cm^-2 s^-1 (',              &
        fuv_band_lo_A(ib), '-', fuv_band_hi_A(ib), ' A)'
   end subroutine band_from_spectrum

   end subroutine input_read

   ! ------------------------------------------------------------------- !

   ! --------------------------------------------------------------------- !
   ! Label matching, shared by the input.inp reader and read_base_inp.
   ! A key matches a line ONLY as a label: anchored at the start of the
   ! left-trimmed line and terminated by a value separator (':', '?',
   ! whitespace, '=', or end-of-line). This makes line order irrelevant and
   ! removes the old substring-collision hazards (e.g.
   ! "Newton solver:" no longer false-matches the "Solver" key). The one
   ! remaining prefix collision, "Wind-AE seed out" vs "Wind-AE seed" (they
   ! share a whitespace-separated prefix), is resolved by testing the longer
   ! key first in the keyword loop above (longest / most-specific first).
   ! --------------------------------------------------------------------- !

   subroutine name_restart_option_change_token(tok)
   ! ONE TOKEN OF "Restart option change".
   !
   ! A token is either one of the physics switches the state file's
   ! '# options' line carries, in which case naming it allows that switch
   ! to differ between the file and this run; or one of the switches that
   ! decide HOW MANY UNKNOWNS the state has, in which case naming it
   ! cannot make the file's rows this run's rows and the run stops; or a
   ! field of the state's configuration that is not an option at all (the
   ! grid, the reservoir, the constants), which no restart carries across;
   ! or a word this code has no meaning for, which stops the run rather
   ! than being ignored.
   character(len=*), intent(in) :: tok
   ! The words that name a part of the configuration OTHER than the
   ! options field: naming them would ask for a state built on another
   ! discretization, another composition or another constant set to be
   ! continued, which needs a conservative remap and not a key.
   character(len=16), parameter :: not_an_option(11) = [ character(len=16)::&
        'N', 'grid', 'R0', 'r_min', 'r_max', 'mode', 'reservoir',           &
        'constants', 'species_columns', 't_phys', 'source' ]
   integer :: i
   do i = 1, n_opt
      if (trim(tok) .ne. trim(opt_name(i))) cycle
      if (opt_changes_layout(i)) then
         write(*,*) '(input_read) ERROR: "Restart option change" names'
         write(*,*) '  "'//trim(tok)//'", which decides how many unknowns'
         write(*,*) '  the state has. The rows of the state file are then'
         write(*,*) '  not the rows of this run, so this is not a restart'
         write(*,*) '  of that state whatever is named: start this'
         write(*,*) '  configuration cold. The tokens that may not be'
         write(*,*) '  named are:'
         call write_option_tokens(.true.)
         error stop 1
      endif
      restart_option_change_named(i) = .true.
      restart_option_change_given    = .true.
      return
   enddo
   do i = 1, size(not_an_option)
      if (trim(tok) .ne. trim(not_an_option(i))) cycle
      write(*,*) '(input_read) ERROR: "Restart option change" names'
      write(*,*) '  "'//trim(tok)//'", which is not a physics option but'
      write(*,*) '  part of the discretization, the composition or the'
      write(*,*) '  constant set the state was built on. A state is not'
      write(*,*) '  carried across such a change by naming it; that'
      write(*,*) '  needs a conservative remap or a cold start.'
      error stop 1
   enddo
   write(*,*) '(input_read) ERROR: "Restart option change" names the'
   write(*,*) '  unknown token "'//trim(tok)//'". The tokens are the'
   write(*,*) '  ones of the "# options" line of a state file:'
   call write_option_tokens(.false.)
   error stop 1
   end subroutine name_restart_option_change_token

   subroutine write_option_tokens(layout_only)
   ! The vocabulary, so that a refusal states what could have been said.
   logical, intent(in) :: layout_only
   character(len=250) :: s
   integer :: i
   s = ' '
   do i = 1, n_opt
      if (layout_only .and. .not. opt_changes_layout(i)) cycle
      s = trim(s)//' '//trim(opt_name(i))
   enddo
   write(*,*) '  '//trim(s)
   end subroutine write_option_tokens

   logical function lbl_match(line, key)
   ! .true. iff adjustl(line) begins with `key` followed by a value separator.
   character(len=*), intent(in) :: line, key
   character(len=250) :: t
   integer :: lk, lt
   lbl_match = .false.
   t  = adjustl(line)
   lk = len(key)
   lt = len_trim(t)
   if (lt .lt. lk) return
   if (t(1:lk) .ne. key) return
   if (lk .eq. lt) then
      lbl_match = .true.                 ! key is the whole trimmed line
   else
      lbl_match = is_sep(t(lk+1:lk+1))    ! next char must be a separator
   endif
   end function lbl_match

   logical function is_sep(c)
   character(len=1), intent(in) :: c
   is_sep = (c .eq. ':' .or. c .eq. '?' .or. c .eq. ' ' .or.   &
             c .eq. char(9) .or. c .eq. '=')
   end function is_sep

   ! ------------------------------------------------------------------- !

   subroutine refuse_duplicate_keys(fname, lines, nl, keys, nk)
   ! Refuse a settings file that states the same key twice.
   !
   ! One key is one quantity, so two lines carrying it are two answers to one
   ! question and the file does not say which is meant. Every reader here
   ! resolves such a pair by letting the LAST line win, which is a rule about
   ! file order rather than about intent -- the run then proceeds on a value
   ! nobody chose and nothing in its log says a choice was made. This stops
   ! instead, and prints both line numbers and both values so the file can be
   ! repaired by whoever knows which was meant.
   !
   ! Each non-blank, non-comment line is attributed to the LONGEST key that
   ! matches it, because the key set is not prefix-free ("Lower atmosphere" is
   ! a prefix of "Lower atmosphere profile", "Wind-AE seed" of "Wind-AE seed
   ! out"): the longest match is the key the parser itself consumes, since the
   ! keyword loop tests the longer label first. A line matching no key is left
   ! to the unknown-line warning.
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: lines(:)
   integer,          intent(in) :: nl, nk
   character(len=*), intent(in) :: keys(:)
   character(len=250) :: t
   integer :: i, kk, lbest, kbest, first_i, n_dup, n_bad
   n_bad = 0
   do kk = 1, nk
      first_i = 0
      n_dup   = 0
      do i = 1, nl
         if (len_trim(lines(i)) .eq. 0) cycle
         t = adjustl(lines(i))
         if (t(1:1) .eq. '#') cycle
         ! longest matching key for this line
         call longest_key_of_line(lines(i), keys, nk, kbest, lbest)
         if (kbest .ne. kk) cycle
         n_dup = n_dup + 1
         if (n_dup .eq. 1) then
            first_i = i
         else
            if (n_dup .eq. 2) then
               write(*,*) '(input_read.f90) ERROR: "'//trim(keys(kk))//     &
                  '" appears more than once in '//trim(fname)//':'
               write(*,'(A,I0,A)') '     line ', first_i, ': '//            &
                  trim(adjustl(lines(first_i)))
            endif
            write(*,'(A,I0,A)') '     line ', i, ': '//                     &
               trim(adjustl(lines(i)))
            n_bad = n_bad + 1
         endif
      enddo
   enddo
   if (n_bad .gt. 0) then
      write(*,*) '   A key states one quantity, so it may appear only'//    &
                 ' once. Delete the line that is not meant'
      write(*,*) '   (a superseded value belongs in a "#" comment, which'// &
                 ' is not parsed).'
      error stop 1
   endif
   end subroutine refuse_duplicate_keys

   ! ------------------------------------------------------------------ !

   subroutine longest_key_of_line(line, keys, nk, kbest, lbest)
   ! Index and length of the LONGEST key in `keys` that matches `line` as a
   ! label; kbest = 0 if none does.
   character(len=*), intent(in)  :: line
   character(len=*), intent(in)  :: keys(:)
   integer,          intent(in)  :: nk
   integer,          intent(out) :: kbest, lbest
   integer :: kk
   kbest = 0;  lbest = 0
   do kk = 1, nk
      if (lbl_match(line, trim(keys(kk)))) then
         if (len_trim(keys(kk)) .gt. lbest) then
            lbest = len_trim(keys(kk));  kbest = kk
         endif
      endif
   enddo
   end subroutine longest_key_of_line

   ! ------------------------------------------------------------------ !

   subroutine read_base_inp
   ! Lower-atmosphere handoff file (optional; a missing file is a no-op).
   ! One "key value" line per entry, '#' comments and unknown keys ignored,
   ! keys matched as labels (lbl_match), exactly as in input.inp. Written by
   ! src/utils/run_lower.py (analytic column) or src/utils/vulcan_to_base.py
   ! (photochemistry); see docs/lower_atmosphere_coupling.* and
   ! docs/input_schema.md section 2c, which carries the same table.
   !
   ! Every key belongs to one of five categories, and the category says what
   ! the value is allowed to do to the wind:
   !
   !  provenance         which code / network / profile produced the file.
   !                     Comments only today; no key is parsed.
   !  EOS boundary       state of the gas at the handoff level:
   !                     T_base [K]     -> T0
   !                     r_base [R_J]   -> R0
   !                     p_base [bar]   -> p_base_bar, the level all of the
   !                                       above refer to
   !                     q_H2_base      -> q_h2_base, the H2 volume mixing
   !                                       ratio of the base gas. It sets the
   !                                       base particle count through
   !                                       comp_ntot_bc AND is imposed on the
   !                                       H2 partition of the inflowing ghost
   !                                       species, the particle count then
   !                                       coming from that same species
   !                                       state, so the equation of state and
   !                                       the chemistry describe one gas
   !                                       (section 117 of Update_EXHALE_stage1).
   !                                       Refused above 0.5/(0.5+He/H).
   !  elemental          the reservoirs the wind transports and redistributes:
   !   reservoir         HeH_base       -> HeH (He/H nuclei)
   !                     <El>_H_base    -> X_<El> (El/H nuclei) for the ten
   !                                       elements of species_table, e.g.
   !                                       "O_H_base 4.90e-4". Overrides
   !                                       metals.inp for that element, and
   !                                       activates the metal system if it
   !                                       is the only nonzero abundance.
   !  initial guess      none today. A key in this category would seed a
   !                     profile the solver is free to move away from.
   !  boundary           Kzz_base [cm2/s] -> he_kzz, the eddy diffusion
   !   constraint        coefficient at the base. It fills kzz_cell, which
   !                     the element transport and the molecular carrier
   !                     transport both read, so it is inert only when
   !                     neither of those runs.
   !
   ! Species mixing ratios other than q_H2_base (q_H2O, q_CO, ...) stay
   ! comments: no part of the code consumes them, so they are diagnostic
   ! metadata and must not be presented as physics.
   !
   ! With a lower-atmosphere PROFILE in use ("Lower atmosphere profile:"),
   ! the EOS-boundary, elemental-reservoir and boundary-constraint keys are
   ! REFUSED: the profile states all of them at the matching level, and a
   ! scalar accepted beside it would be a second source of the same
   ! quantities. What remains to check is then the pair as a whole, and the
   ! one provenance comment that is parsed, "# solution_id <hash>", is what
   ! makes "these two files are the same lower-atmosphere solution" a
   ! statement rather than an assumption.
   character(len=250) :: line
   character(len=:), allocatable :: str
   character(len=32)  :: key
   character(len=128) :: id_base
   logical :: ex
   integer :: ios, ub, ie, ip
   integer :: nb, ib, nkb
   real*8  :: ab
   character(len=250), allocatable :: baselines(:)
   character(len=32),  allocatable :: base_keys(:)

   inquire(file='base.inp', exist=ex)
   if (.not. ex) return
   write(*,*) '(input_read) Reading base.inp (lower-atmosphere handoff)..'

   ! ----- One line per key, here too -----
   ! Same rule and same routine as input.inp: a handoff that states the base
   ! temperature (or pressure, or an elemental ratio) twice does not say
   ! which level it was written at, and last-line-wins would pick one.
   nkb = 6 + n_melem
   allocate(base_keys(nkb))
   base_keys(1) = 'T_base';     base_keys(2) = 'r_base'
   base_keys(3) = 'HeH_base';   base_keys(4) = 'Kzz_base'
   base_keys(5) = 'q_H2_base';  base_keys(6) = 'p_base'
   do ie = 1, n_melem
      base_keys(6+ie) = trim(melem_name(ie))//'_H_base'
   enddo
   open(newunit=ub, file='base.inp', status='old')
   nb = 0
   do
      read(ub,'(A)',iostat=ios) line
      if (ios .ne. 0) exit
      nb = nb + 1
   enddo
   allocate(character(len=250) :: baselines(nb))
   rewind(ub)
   do ib = 1, nb
      read(ub,'(A)') baselines(ib)
   enddo
   close(ub)
   call refuse_duplicate_keys('base.inp', baselines, nb, base_keys, nkb)
   deallocate(baselines, base_keys)
   id_base = ''
   open(newunit=ub, file='base.inp', status='old')
   do
      read(ub,'(A)',iostat=ios) line
      if (ios .ne. 0) exit
      if (len_trim(line) .eq. 0) cycle
      if (index(adjustl(line),'#') .eq. 1) then
         ! Provenance comments. The one that is parsed is "# solution_id
         ! <hash>": it is what makes "the profile and the scalar file are the
         ! same lower-atmosphere solution" a checkable statement rather than
         ! an assumption.
         ip = index(line, 'solution_id')
         if (ip .gt. 0) id_base = adjustl(line(ip+len('solution_id'):))
         cycle
      endif
      if (lbl_match(line,'T_base')) then
         call refuse_scalar_key('T_base', 'EOS boundary')
         str = get_word(line,2);  read(str,*) T0
         write(*,'(A,F9.1,A)') '   base.inp: T0 -> ', T0, ' K'
      else if (lbl_match(line,'r_base')) then
         call refuse_scalar_key('r_base', 'EOS boundary')
         str = get_word(line,2);  read(str,*) R0
         write(*,'(A,F8.4,A)') '   base.inp: R0 -> ', R0, ' R_J'
      else if (lbl_match(line,'HeH_base')) then
         call refuse_scalar_key('HeH_base', 'elemental reservoir')
         str = get_word(line,2);  read(str,*) HeH
         if (HeH .gt. 0.0d0) thereis_He = .true.
         write(*,'(A,F8.5)') '   base.inp: He/H -> ', HeH
      else if (lbl_match(line,'Kzz_base')) then
         call refuse_scalar_key('Kzz_base', 'boundary constraint')
         str = get_word(line,2);  read(str,*) he_kzz
         write(*,'(A,ES9.2,A)') '   base.inp: He_Kzz -> ', he_kzz, ' cm2/s'
      else if (lbl_match(line,'q_H2_base')) then
         call refuse_scalar_key('q_H2_base', 'EOS boundary')
         str = get_word(line,2);  read(str,*) q_h2_base
         write(*,'(A,F8.5)') '   base.inp: q_H2(base) -> ', q_h2_base
      else if (lbl_match(line,'p_base')) then
         call refuse_scalar_key('p_base', 'EOS boundary')
         str = get_word(line,2);  read(str,*) p_base_bar
         base_level_from_handoff = .true.
         base_level_source       = 'base.inp p_base'
         write(*,'(A,ES9.2,A)') '   base.inp: p_base -> ', p_base_bar, ' bar'
      else
         ! Elemental reservoirs "<El>_H_base": El/H nuclei ratio at the
         ! handoff level, for any element species_table knows. Element
         ! symbols are unique, so the ten labels cannot collide with each
         ! other or with the keys above.
         do ie = 1, n_melem
            key = trim(melem_name(ie))//'_H_base'
            if (lbl_match(line, trim(key))) then
               call refuse_scalar_key(trim(key), 'elemental reservoir')
               str = get_word(line,2);  read(str,*) ab
               call set_element_abundance(ie, ab)
               write(*,'(A,A,A,ES10.3)') '   base.inp: ',                  &
                  trim(melem_name(ie)), '/H -> ', ab
               exit
            endif
         enddo
      endif
   enddo
   close(ub)

   ! Consistency echo, printed once with the final values (the keys may come
   ! in any order).  q_H2_base and HeH_base are NOT independent -- both are
   ! read off the same lower-atmosphere solution -- but only HeH enters the
   ! elemental budget, so a mismatched pair would otherwise pass unnoticed.
   ! Print only: an inconsistent handoff is the user's to judge, never a
   ! reason to refuse to run.
   if (q_h2_base .gt. 0.0d0) then
      write(*,'(A,F8.5,A,F8.5)') '   base.inp: photochemical q_H2 =',      &
         q_h2_base, ' is used with He/H =', HeH
      write(*,*) '     (both must come from the same lower-atmosphere'//   &
                 ' solution)'
   endif
   if (abs(p_base_bar/1.0d-6 - 1.0d0) .gt. 1.0d-6) then
      write(*,'(A,ES9.2,A)') '   base.inp: NOTE the handoff level p_base =',&
         p_base_bar, ' bar is not the'
      write(*,*) '     standard 1 microbar; the base composition and the'
      write(*,*) '     chemical-equilibrium H2 fit both refer to that level.'
   endif

   ! ---- the pair must be one solution -------------------------------- !
   ! A profile and a base.inp side by side describe the same gas at the same
   ! level. If they came from different lower-atmosphere solutions the run is
   ! not a model of anything, so this is a refusal and not a message.
   if (lap_in_use) then
      if (len_trim(id_base) .eq. 0) then
         write(*,*) '(input_read) ERROR: base.inp sits beside the lower-'// &
                    'atmosphere profile'
         write(*,*) '  '//trim(lap_file)//' but carries no'//              &
                    ' "# solution_id <hash>" line, so there is no way to'
         write(*,*) '  tell whether the two describe the same lower-'//    &
                    'atmosphere solution. Add the'
         write(*,*) '  line (the adapter writes it), or delete base.inp.'
         error stop 1
      endif
      if (trim(id_base) .ne. trim(lap_solution_id)) then
         write(*,*) '(input_read) ERROR: the lower-atmosphere profile and'//&
                    ' base.inp are different solutions.'
         write(*,*) '  profile  solution_id: '//trim(lap_solution_id)
         write(*,*) '  base.inp solution_id: '//trim(id_base)
         error stop 1
      endif
      write(*,*) '   base.inp: solution_id matches the profile.'
   endif

   contains

   subroutine refuse_scalar_key(kname, cat)
   ! With a profile in use, the three physics categories of base.inp have a
   ! single source and it is the profile. Accepting a scalar beside it would
   ! reintroduce exactly the same-solution problem the profile removes by
   ! construction.
   character(len=*), intent(in) :: kname, cat
   if (.not. lap_in_use) return
   write(*,*) '(input_read) ERROR: base.inp key "'//trim(kname)//          &
              '" ('//trim(cat)//')'
   write(*,*) '  is refused because the lower-atmosphere profile'
   write(*,*) '  '//trim(lap_file)//' is in use and already states that'
   write(*,*) '  quantity at the matching level. Delete the key (or the'
   write(*,*) '  whole base.inp) and let the profile be the single source;'
   write(*,*) '  provenance comments and diagnostic keys are still allowed.'
   error stop 1
   end subroutine refuse_scalar_key

   end subroutine read_base_inp

   ! ------------------------------------------------------------------- !

   subroutine apply_lower_atmosphere_profile
   ! Set the base state and the elemental reservoirs from the lower-
   ! atmosphere solution at the matching level.
   !
   ! The values go in through exactly the doors the scalar base.inp keys use
   ! -- direct assignment for the EOS anchors, set_element_abundance for the
   ! elements -- so comp_mass_per_H, comp_ntot_bc and comp_rho_bc build the
   ! base EOS from the numbers the lower model reported and nothing in those
   ! functions changes. n_tot and rho are carried by the file as the lower
   ! model's own values and are NOT imposed here: the base density is n0 of
   ! input.inp, and the base particle count and mass follow from T0, HeH,
   ! q_H2 and the elemental ratios exactly as they do on the scalar path.
   ! Placed after read_base_inp so nothing scalar can overwrite a profile
   ! value, and before the composition block so melem_ab, thereis_metals and
   ! thereis_lowIP_metal are derived from the handoff abundances.
   real*8  :: val
   logical :: got
   integer :: ie

   if (.not. lap_in_use) return

   call lap_value_at_match('T', val, got)
   if (got) then
      T0 = val
      write(*,'(A,F9.1,A)') '   profile: T0 -> ', T0, ' K'
   endif
   call lap_value_at_match('r', val, got)
   if (got) then
      R0 = val
      write(*,'(A,F8.4,A)') '   profile: R0 -> ', R0, ' R_J'
   endif
   ! The matching level is the base level: T, r, q_H2 and every elemental
   ! ratio above were read AT p_match, so the wind starts there and n0
   ! follows from it in the base-level block of input_read, exactly as it
   ! follows from base.inp's p_base (rev 3 section 10.2 item 4).
   p_base_bar              = lap_p_match_bar
   base_level_from_handoff = .true.
   base_level_source       = 'the profile matching level p_match'
   write(*,'(A,ES9.2,A)') '   profile: p_base -> ', p_base_bar, ' bar'
   call lap_value_at_match('q_H2', val, got)
   if (got) then
      q_h2_base = val
      write(*,'(A,F8.5)') '   profile: q_H2(base) -> ', q_h2_base
   endif
   call lap_element_ratio_at_match('He', val, got)
   if (got) then
      HeH = val
      if (HeH .gt. 0.0d0) thereis_He = .true.
      write(*,'(A,F8.5)') '   profile: He/H -> ', HeH
   endif
   do ie = 1, n_melem
      call lap_element_ratio_at_match(trim(melem_name(ie)), val, got)
      if (got) then
         call set_element_abundance(ie, val)
         write(*,'(A,A,A,ES10.3)') '   profile: ',                        &
            trim(melem_name(ie)), '/H -> ', val
      endif
   enddo

   write(*,'(A,A,A,I0,A)') '   profile: source ', trim(lap_source_code),  &
      ', closure iteration ', lap_iteration,                              &
      ' (solution_id in EXHALE_resolved.out)'

   ! "He_Kzz:" in input.inp is the constant a run states when it has no
   ! profile. It is not refused -- it is a different file with a different
   ! contract from base.inp -- but with a profile in use it is inert, and a
   ! key that looks effective and is not must say so.
   if (he_kzz .gt. 0.0d0) then
      write(*,'(A)') ' (input_read) WARNING: "He_Kzz:" is set but the'//   &
         ' lower-atmosphere profile carries'
      write(*,'(A)') '   the eddy diffusion coefficient as a profile;'//   &
         ' the constant is ignored.'
   endif

   end subroutine apply_lower_atmosphere_profile

   ! ------------------------------------------------------------------- !

   subroutine set_element_abundance(ie, ab)
   ! Store an El/H nuclei ratio in the named scalar of element index ie
   ! (canonical iel_* order of species_table). The named scalars are what
   ! metals.inp writes and what the melem_ab array is built from, so a
   ! handoff abundance enters through exactly the same door as metals.inp.
   ! Passing through this door is also what marks the element as stated by the
   ! handoff (melem_from_handoff), which is what lets a restart column be
   ! renormalized onto the new reservoir in load_IC.
   integer, intent(in) :: ie
   real*8,  intent(in) :: ab
   if (allocated(melem_from_handoff)) melem_from_handoff(ie) = .true.
   select case (ie)
      case (iel_C);  X_C  = ab
      case (iel_O);  X_O  = ab
      case (iel_N);  X_N  = ab
      case (iel_Mg); X_Mg = ab
      case (iel_Si); X_Si = ab
      case (iel_Ca); X_Ca = ab
      case (iel_Na); X_Na = ab
      case (iel_K);  X_K  = ab
      case (iel_S);  X_S  = ab
      case (iel_Fe); X_Fe = ab
   end select
   end subroutine set_element_abundance

   ! ------------------------------------------------------------------- !

   subroutine run_lower_atm_prestep
   ! Invoke the lower-atmosphere generator (VULCAN or the analytic column)
   ! when "Lower atmosphere: vulcan|analytic <R_1bar[R_J]>" is set and no
   ! base.inp is present.  The EXHALE code root is taken from the
   ! EXHALE_ROOT environment variable; it must be set to the EXHALE install
   ! directory so the generator scripts under src/utils/ can be located.
   character(len=1024) :: root, cmd
   character(len=32)   :: r1str
   logical :: ex
   integer :: rc, cst

   if (lower_atm_mode .le. 0) return
   inquire(file='base.inp', exist=ex)
   if (ex) then
      write(*,*) '(input_read) Lower atmosphere: existing base.inp found'//&
                 ' -- using it (delete it to regenerate).'
      return
   endif
   if (lower_atm_r1bar .le. 0.0d0) then
      write(*,*) '(input_read) ERROR: "Lower atmosphere:" needs the 1-bar'//&
                 ' (transit) radius, e.g. "Lower atmosphere: vulcan 1.36".'
      error stop 1
   endif

   call get_environment_variable('EXHALE_ROOT', root)
   if (len_trim(root) .eq. 0) then
      write(*,*) '(input_read) ERROR: "Lower atmosphere:" requires the'//  &
                 ' EXHALE_ROOT environment variable to locate the'
      write(*,*) '  generator scripts under src/utils/. Set it to the'//   &
                 ' EXHALE install directory, e.g.'
      write(*,*) '  export EXHALE_ROOT=/path/to/EXHALE'
      error stop 1
   endif
   write(r1str,'(F0.5)') lower_atm_r1bar

   if (lower_atm_mode .eq. 2) then
      write(*,*) '(input_read) Lower atmosphere: running VULCAN'//        &
                 ' photochemistry (first run takes hours; cached after).'
      cmd = 'python3 '//trim(root)//'/src/utils/vulcan_driver.py . '//   &
            '--r1bar '//trim(r1str)
   else
      write(*,*) '(input_read) Lower atmosphere: analytic column.'
      cmd = 'python3 '//trim(root)//'/src/utils/run_lower.py . '//       &
            '--r1bar '//trim(r1str)
   endif
   call execute_command_line(trim(cmd), exitstat=rc, cmdstat=cst)
   inquire(file='base.inp', exist=ex)
   if (rc .ne. 0 .or. cst .ne. 0 .or. .not. ex) then
      write(*,*) '(input_read) ERROR: lower-atmosphere generator failed'//&
                 ' (see messages above); no base.inp produced.'
      error stop 1
   endif
   end subroutine run_lower_atm_prestep
      
   ! ------------------------------------------------------- !
      
   function get_word(string_in,n_word)
	! Function to read the nth_word in the current string
	! 	"Words" are separated by spaces
	
	character(len = *), intent(in) :: string_in
	integer, intent(in) :: n_word
	
	character(len = :), allocatable  :: string
	character(len = 300) :: c_string	
	character :: p_char,c_char
	integer :: counter
	integer :: c_word_counter
	integer :: str_len
	
	character(len = :), allocatable :: get_word
	
	! Initialize counters and strings. The result is allocated here, so a
	! string that holds fewer than n_word words returns the EMPTY string
	! rather than an unallocated deferred-length result, which the standard
	! leaves undefined; the callers test the result with len_trim.
	get_word       = ''
	counter        = 1
	c_word_counter = 0
	string   = trim(string_in)
	str_len  = len(string)
	p_char = ''
	c_char = ''
	c_string = ''
	
	! Loop inside the string	
	do while (counter .ge. 0 .and. counter .le. str_len)

		! Characters
		if (counter .ge. 2) then	! Skip if its the first iteration
			p_char = string(counter - 1:counter - 1)
		endif
		c_char = string(counter:counter)
		
		! If a character is found
		if (c_char .ne. '') then	
			
			! Attach character to current string
			c_string = trim(c_string) // c_char
			
			! Update counter and continue
			counter = counter + 1 
		
			continue			
		
		else
		
			! If it's a first space after a character 		
			if (p_char .ne. '') then		
				
				! Update word counter
				c_word_counter = c_word_counter + 1
				
				! Exit from loop if word counter 
				! is equal to n_word in input				
				if (c_word_counter .eq. n_word) then 
					get_word = trim(c_string)
					return
				endif
			
				! Reset current string	
				c_string = ''
			
				! Update counter
				counter = counter + 1
				
			else	! If multiple spaces
				
				counter = counter + 1
				continue
			endif
				
		endif
		
		! If last character, return
		if (counter .eq. str_len) then 
			
			! Update string 
			c_char = string(counter:counter)
			c_string = trim(c_string) // c_char
			
			! Update word counter
			c_word_counter = c_word_counter + 1
				
			! Exit from loop if word counter 
			! is equal to n_word in input				
			if (c_word_counter .eq. n_word) then 
				get_word = trim(c_string)
				return
			endif
		endif
		
	enddo	! End of while loop
	
	! End of get_word function
	end function
      
	!------------------------------------------!

	subroutine retired_base_key(key)
	! A lower-boundary key that no longer exists is an ERROR, not a warning.
	!
	! Every one of these named a component of the ghost-cell closure that the
	! characteristic face condition replaced: a ghost velocity copied from the
	! interior and smoothed (Valve eps), a ghost pressure extrapolated from the
	! interior gradient (Hydrostatic base), a ghost temperature copied from
	! cell 1 (Base ghost temperature: continuous), a ghost velocity taken from
	! the wind's flux constant (Base velocity: massflux). None of them has a
	! meaning in a boundary condition that states a reservoir at the face and
	! takes one relation from the interior, and a run that silently ignored one
	! would answer a question the input did not ask. So it stops.
	!
	! The base LEVEL is unaffected: "Base BC: density|pressure" still states
	! where the boundary is, and it is the key that carries over.
	character(len=*), intent(in) :: key
	write(*,*)
	write(*,'(A)') ' (input_read.f90) ERROR: "'//trim(key)//':" is no'//   &
	   ' longer a key of this code.'
	write(*,'(A)') '   The lower boundary is now a characteristic'//       &
	   ' condition imposed at the face r_edg(0):'
	write(*,'(A)') '   the reservoir states the pressure and the entropy'//&
	   ' of the lower atmosphere at the'
	write(*,'(A)') '   base level, and the outgoing acoustic invariant of'//&
	   ' the first interior cell states'
	write(*,'(A)') '   the velocity. The ghost cells are the volume'//     &
	   ' averages of that face state continued'
	write(*,'(A)') '   below it, so there is no separate ghost velocity,'//&
	   ' pressure or temperature closure'
	write(*,'(A)') '   left to select.'
	write(*,'(A)') '   Remove the line. Use "Base BC: pressure'//          &
	   ' [<p_ubar>]" to state the base LEVEL, and'
	write(*,'(A)') '   base.inp / "Lower atmosphere profile:" to state'//  &
	   ' the base temperature and composition.'
	write(*,*)
	error stop 1
	end subroutine retired_base_key

    ! End of module
	end module Read_input

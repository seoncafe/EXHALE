   module  Read_input
   ! Read input planetary parameters adn define

   use global_parameters
   use metals_input        ! optional metals.inp abundance reader
   use charge_exchange, only: cx_init,       &  ! build active charge-exchange set
                              he_h_charge_exchange  ! He <-> H pair (group B) switch
   use species_table, only: n_melem, iel_C, iel_O, iel_N, iel_Mg,  &
                            iel_Si, iel_Ca, iel_Na, iel_K, iel_S,  &
                            iel_Fe, mion_ethr, melem_i0, melem_A
   use composition, only: comp_mass_per_H, comp_ntot_bc, comp_rho_bc,   &
                          h2_mixing_ratio_base

   implicit none
      
   contains
      
   subroutine input_read
   ! Subroutine to read the input file and assign names and values
   !    to global constants

      character(len = :), allocatable :: str
      character(len = 250)            :: line
      character(len = 250), allocatable :: filelines(:)
      integer                         :: ios
      integer                         :: im
      integer                         :: i, nlines
      integer                         :: kk
      logical                         :: is_known

   ! Every label input_read recognizes: the core block followed by the
   ! keyword-extension block, in the same order as the reads below. Used only
   ! by the trailing unknown-line scan to WARN (never stop) on a non-blank,
   ! non-'#' line that matches no key -- e.g. a "Newton Solver:" (capital S)
   ! typo that anchored matching would otherwise silently ignore. The energy-
   ! band line ("[E_low...") is a bracket prefix, checked separately there.
   character(len=32), parameter :: known_keys(*) = [ character(len=32) ::    &
      'Planet name', 'Log10 lower boundary', 'Planet radius', 'Planet mass', &
      'Equilibrium temperature', 'Orbital distance', 'Escape radius',        &
      'He/H number ratio', '2D approximate method', 'Parent star mass',      &
      'Spectrum type', 'Spectrum file', 'Power-law index', 'Photon energy',  &
      'Use only EUV', 'Log10 of X-ray luminosity', 'Log10 of EUV luminosity',&
      'Grid type', 'Base grid',                                              &
      'Numerical flux', 'Reconstruction scheme', 'Include He23S',            &
      'Load IC', 'Do only PP', 'Force start',                                &
      'Domain mode', 'Outer radius', 'Stellar Teff', 'Stellar radius',       &
      'Deexc heat', 'Wind-AE seed out', 'Wind-AE seed', 'Jlya RT file',      &
      'Jlya escape-prob', 'Stellar Lya flux', 'Lya stellar halfwidth',       &
      'Lya stellar boost', 'du_th', 'ATES_photoionization_rate',             &
      'Legacy_HHe_rates', 'Secondary_ionization', 'He_rec_coupling',         &
      'He_H_charge_exchange',                                                &
      'Molecular chemistry', 'Molecular base', 'Stellar LW flux',           &
      'Lower atmosphere',                                                   &
      'Lower column', 'He_Kzz', 'He_alphaT', 'He_ambipolar',                 &
      'He_metal_diffusion', 'He_diffusion', 'Stall', 'Energy solver',        &
      'Time stepping', 'Level tol', 'Solver', 'Valve eps', 'Hydrostatic base',&
      'Shapiro filter', 'Base BC', 'Base velocity', 'Viscosity',             &
      'Base ghost temperature', 'Max steps', 'Coronal cutoff width',         &
      'Base IR field',                                                       &
      'Conduction', 'Resid tol',                                             &
      'Resid norm', 'CFL', 'Transonic IC', 'Hot Parker IC', 'IC mode',       &
      'Newton solver', 'Brent solver' ]

   ! ----- Read planetary parameters from input file ----- !

   ! Metal abundances (default: no metals). Overridden at runtime by an
   ! optional metals.inp file (no recompile); thereis_metals is set from
   ! the resulting values below.
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
   ! extension block that follows uses the same anchored matching. A
   ! duplicated key resolves to its LAST occurrence (find_lbl below), matching
   ! the keyword loop, which overwrites on each match.
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

   ! ----- Core block (label-matched; order-independent) -----

   ! Planet name (single token at word 3)
   p_name = get_word(req('Planet name'), 3)

   ! Log10 of n0
   str = get_word(req('Log10 lower boundary'), 7);  read(str,*) n0

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
   ! 13.6 eV HI edge (e.g. Mg I at 7.646 eV) needs the below-threshold
   ! sub-grid extension in set_energy_vectors (mutually exclusive with
   ! the HeI triplet). Triggered by any such active element.
   do im = 1, n_melem
      if (melem_ab(im) .gt. 0.0d0 .and.                       &
          mion_ethr(melem_i0(im)) .lt. e_th_HI)               &
         thereis_lowIP_metal = .true.
   enddo

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
         write(*,*) '   Allowed values: Load, Power-law, Monochromatic.'
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
   rec_method = get_word(req('Reconstruction scheme'), 3)
   if (rec_method.eq.'WENO3') use_weno3 = .true.
   if (rec_method.eq.'PLM')   use_plm = .true.
   if (rec_method.eq.'PLM+WENO3') then
      use_plm         = .true.    ! start in PLM; switch to WENO3 mid-run
      recon_two_stage = .true.
   endif

   ! Include He23S
   str = get_word(req('Include He23S'), 3)
   if (str .eq. 'True')  thereis_HeITR = .true.

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
		transonic_ic     = .false.
		hot_parker_ic    = .false.
		T_wind_ic        = 1.0d4
		use_newton_ieq   = .true.    ! upgraded solvers are the default
		use_brent_tsolve = .true.
		windae_seed_file = 'inputdata/windae_seed.csv'
		windae_seed_out  = ''
		hydrostatic_base = .false.
		dr_base          = 2.0e-4     ! uniform base cell size [R_p] (Mixed grid);
		                              ! default-real literal on purpose, see parameters.f90
		N_low_cells      = 50         ! number of uniform base cells
		base_bc_mode     = 0          ! density-anchored base (legacy) by default
		resid_vol        = .true.     ! volume-weighted residual norm by default
		ates_photoion_rate = .false.  ! default: Verner+1996 He I (1^1S) photoion.
		legacy_hhe_rates   = .false.  ! default: Badnell/Mao + Voronov H/He rates
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
						   ' -> single-stage PLM.'
						du_th_plm = -1.0d0
					endif
				else
					! PLM or WENO3: single-stage; ignore any second du_th value.
					du_th_plm = -1.0d0
				endif
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
			else if (lbl_match(line, 'Stellar LW flux')) then
				! "Stellar LW flux [erg/cm2/s]: <F>" -- band-integrated
				! stellar flux in the H2 Lyman-Werner bands (912-1110 A) at
				! the planet's orbit. Drives H2 photodissociation in the
				! molecular network (lyman_werner.f90). 0 = off (default).
				str = get_word(line, 5)
				read(str,*) F_LW_star
			else if (lbl_match(line, 'Molecular base')) then
				! EOS-only molecular base (docs/lower_atmosphere_*).
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true') molecular_base = .true.
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
				! Eddy diffusion coefficient [cm^2/s] for He/H separation.
				! "He_Kzz: 1.0e9"
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
				! "Valve eps: <v_eps>" smooths the base one-way valve
				! (softplus; <= 0 keeps the exact legacy max(v,0)).
				str = get_word(line, 3);  read(str,*) valve_eps
				write(*,'(A,ES9.2)') ' (input_read) Smooth base valve, eps =', valve_eps
			else if (lbl_match(line, 'Hydrostatic base')) then
				str = get_word(line, 3)
				if (str .eq. 'True') hydrostatic_base = .true.
				if (hydrostatic_base) write(*,'(A)') ' (input_read) '//   &
				   'Hydrostatic base ghost cells enabled'
			else if (lbl_match(line, 'Shapiro filter')) then
				str = get_word(line, 3);  read(str,*) shapiro_eps
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) shapiro_every
				if (shapiro_eps .gt. 0.0d0) write(*,'(A,ES9.2,A,I0,A)')  &
				   ' (input_read) Shapiro filter eps =', shapiro_eps,    &
				   ', every ', shapiro_every, ' steps'
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
				! "Base ghost temperature: isothermal|continuous" selects the
				! temperature closure of the lower ghost cells (the density
				! anchor rho_bc is unaffected). isothermal (default) pins
				! T = T0; continuous imposes dT/dr = 0 at the base face, so the
				! ghost carries T(cell 1). See parameters.f90 and
				! docs/hd189_base_checkerboard.md sec. 4.2 and 13.
				str = get_word(line, 4)
				if (str .eq. 'continuous') base_ghost_T_continuous = .true.
				if (str .eq. 'isothermal') base_ghost_T_continuous = .false.
				if (str .ne. 'continuous' .and. str .ne. 'isothermal')       &
					write(*,*) '(input_read.f90) WARNING: unknown "Base '//  &
					   'ghost temperature: ', trim(str), '"; keeping isothermal.'
				if (base_ghost_T_continuous) write(*,'(A)') ' (input_read) '// &
				   'Base ghost temperature: continuous (dT/dr = 0, T_ghost = T_1)'
			else if (lbl_match(line, 'Max steps')) then
				! "Max steps: <N>" overrides the hard cap on marching
				! iterations (default 1000000).
				str = get_word(line, 3);  read(str,*) count_max
				write(*,'(A,I0)') ' (input_read) Max marching steps =', count_max
			else if (lbl_match(line, 'Coronal cutoff width')) then
				! "Coronal cutoff width: <w>" sets the roll-off width of the
				! coronal-excitation guard below the 1e3 K CHIANTI fit floor
				! (Cool_coeff.f90). Default 0.1; justified window 0.08-0.13
				! (docs/coronal_cutoff_width.md).
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
			else if (lbl_match(line, 'Base velocity')) then
				str = get_word(line, 3)
				if (str .eq. 'valve')    base_v_massflux = .false.
				if (str .eq. 'massflux') base_v_massflux = .true.
				if (base_v_massflux) then
					write(*,'(A)') ' (input_read) Base velocity from '//  &
					   'mass-flux F_c (CETIMB-style)'
				else
					write(*,'(A)') ' (input_read) Base velocity: legacy valve'
				endif
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
				! H = kT/(mu g); see docs/hd189_base_checkerboard.md. Ignored by
				! the Uniform and Stretched grid types.
				str = get_word(line, 4);  read(str,*) dr_base
				str = get_word(line, 5)
				if (len_trim(str) .gt. 0) read(str,*) N_low_cells
				write(*,'(A,ES9.2,A,I0,A,ES9.2,A)')                       &
				   ' (input_read) Base grid: dr =', dr_base,              &
				   ' R_p x ', N_low_cells, ' cells (uniform region ',     &
				   dr_base*N_low_cells, ' R_p)'
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
				! "Resid norm: vol|Linf" -- residual norm for convergence.
				! vol (default) = volume-weighted; Linf = legacy max-over-cells.
				str = get_word(line, 3)
				if (str .eq. 'Linf' .or. str .eq. 'linf' .or. str .eq. 'LINF') &
					resid_vol = .false.
				if (str .eq. 'vol' .or. str .eq. 'volume') resid_vol = .true.
				write(*,'(A,L1)') ' (input_read) Volume-weighted residual norm: ', resid_vol
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

		! ----- Base ghost-pressure closures are mutually exclusive -----
		! Both keys set the SAME quantity, the ghost pressure: hydrostatic_base
		! extrapolates the interior gradient (fixing neither T nor p), while
		! base_ghost_T_continuous fixes T_ghost = T_1. hydrostatic_base is
		! tested first in BC_component_constrho, so say so rather than let the
		! continuous-T key look effective.
		if (hydrostatic_base .and. base_ghost_T_continuous) then
			write(*,*) '(input_read.f90) WARNING: "Hydrostatic base: True" '// &
			   'and "Base ghost temperature: continuous" both set; the'
			write(*,*) '  hydrostatic ghost pressure wins and the temperature'//&
			   ' closure is ignored.'
		endif
		if (base_bc_mode .eq. 1 .and. base_ghost_T_continuous) then
			write(*,*) '(input_read.f90) NOTE: "Base BC: pressure" derives n0'//&
			   ' from the target base pressure AT T0; with a continuous-T'
			write(*,*) '  ghost the base pressure then floats with T(cell 1),'//&
			   ' so only the base DENSITY stays anchored.'
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

   close(unit = 1)
	write(*,*) '(input_read.f90) Done'

   !------ Definition of physical parameters ------!
      
   n0     = 10.0**(n0)
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
   call read_base_inp

   ! ---- Composition reconciliation and validation. Placed here so every
   ! flag (thereis_He, HeH, thereis_HeITR, thereis_mol, he_diffusion,
   ! thereis_metals) has its final value: base.inp above can still flip
   ! thereis_He/HeH, and the trailing keyword scan sets thereis_mol and
   ! he_diffusion. This is the single authoritative check before N_eq sizing.

   ! Remove HeITR chemistry if He is not included
   if (.not. thereis_He) thereis_HeITR = .false.

   ! molecular chemistry constraints: requires He. Trace metals are solved
   ! together with the molecular network (System_HeH_mol_metals), so the two
   ! are no longer exclusive; He/H diffusion still is.
   if (thereis_mol .and. .not. thereis_He) then
      write(*,*) '(input_read) ERROR: Molecular chemistry needs He/H>0.'
      error stop 1
   endif
   if (thereis_mol .and. he_diffusion) then
      write(*,*) '(input_read) ERROR: Molecular chemistry + He_diffusion'//&
                 ' not supported yet.'
      error stop 1
   endif

   ! Lyman-Werner photodissociation acts on H2, which only exists as a
   ! solved species with the molecular network on. Report rather than stop:
   ! the key is then simply inert.
   if (F_LW_star .gt. 0.0d0) then
      if (thereis_mol) then
         write(*,'(A,ES10.3,A)') ' (input_read) Lyman-Werner band flux at'//&
            ' the planet: ', F_LW_star, ' erg cm^-2 s^-1'
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
   ! legacy H/He-only values (mass_per_H = 1+4*HeH, ntot_bc = 1).  With
   ! "Molecular base: True" comp_ntot_bc also removes the H nuclei bound
   ! into H2 at the base (passive molecular base, EOS-only: the species
   ! arrays stay atomic; docs/lower_atmosphere_coupling.*).
   ! Routed through the composition module (single source of the base
   ! composition policy). comp_* reproduce the legacy expressions bitwise;
   ! eos_metals / metals-present / molecular_base branching lives inside.
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

   ! Pressure-anchored base (Base BC: pressure): override n0 so that the base
   ! pressure n0*kb*T0*ntot_bc matches the target base_p_ubar [microbar].
   ! 1 microbar = 1 erg/cm^3. (dp_bc, the small electron term, is negligible
   ! here and is added by the ghost BC later.) This is the CETIMB 1-microbar
   ! lower boundary: a much less dense base, hence a weaker rho*g source.
   if (base_bc_mode .eq. 1) then
      n0 = base_p_ubar/(kb_erg*T0*ntot_bc)
      write(*,'(A,ES12.4,A)') ' (input_read) Base BC pressure mode: '//   &
         'derived n0 =', n0, ' cm^-3'
   endif

   v0     = sqrt(kb_erg*T0/mu)
   t_s    = R0/v0
   p0     = n0*mu*v0*v0
   q0     = n0*mu*v0*v0*v0/R0
   b0     = (Gc*Mp*mu)/(kb_erg*T0*R0)
   dp_bc  = 1.0e-10
   ! Base particle count seen by the continuous-temperature ghost before the
   ! first composition solve refreshes it (Apply_BC can run first, e.g. in the
   ! IC/residual paths). With the base value the ghost then simply copies the
   ! cell-1 pressure, which is the T_ghost = T_1 statement for base composition.
   n_part_cell1 = ntot_bc + dp_bc

   !------ Allocations ------!
      
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

   ! End of subroutine

   contains

   ! --------------------------------------------------------------------- !
   ! input.inp label matching. A key matches a line ONLY as a label:
   ! anchored at the start of the left-trimmed line and terminated by a value
   ! separator (':', '?', whitespace, '=', or end-of-line). This makes line
   ! order irrelevant and removes the old substring-collision hazards (e.g.
   ! "Newton solver:" no longer false-matches the "Solver" key). The one
   ! remaining prefix collision, "Wind-AE seed out" vs "Wind-AE seed" (they
   ! share a whitespace-separated prefix), is resolved by testing the longer
   ! key first in the keyword loop above (longest / most-specific first).
   ! --------------------------------------------------------------------- !

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

   end subroutine input_read

   ! ------------------------------------------------------------------- !

   subroutine read_base_inp
   ! lower-atmosphere handoff file (optional).  Keyword lines:
   !   T_base    <K>      -> overrides T0 (base temperature)
   !   r_base    <R_J>    -> overrides the "Planet radius" (1-ubar radius)
   !   HeH_base  <ratio>  -> overrides the He/H number ratio
   !   Kzz_base  <cm2/s>  -> sets he_kzz (used by He_diffusion)
   !   q_H2_base <ratio>  -> photochemical H2 volume mixing ratio at the base;
   !                         replaces the chemical-equilibrium fit in the
   !                         molecular-base particle count
   !   p_base    <bar>    -> pressure level the handoff describes (default
   !                         1e-6 bar); the equilibrium fit is evaluated there
   ! '#' comments and unknown keys are ignored.  Written by
   ! src/utils/run_lower.py (analytic column) or by an external
   ! photochemical/RC model; see docs/lower_atmosphere_coupling.*.
   character(len=250) :: line
   character(len=:), allocatable :: str
   logical :: ex
   integer :: ios, ub

   inquire(file='base.inp', exist=ex)
   if (.not. ex) return
   write(*,*) '(input_read) Reading base.inp (lower-atmosphere handoff)..'
   open(newunit=ub, file='base.inp', status='old')
   do
      read(ub,'(A)',iostat=ios) line
      if (ios .ne. 0) exit
      if (len_trim(line) .eq. 0) cycle
      if (index(adjustl(line),'#') .eq. 1) cycle
      if (index(line,'T_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) T0
         write(*,'(A,F9.1,A)') '   base.inp: T0 -> ', T0, ' K'
      else if (index(line,'r_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) R0
         write(*,'(A,F8.4,A)') '   base.inp: R0 -> ', R0, ' R_J'
      else if (index(line,'HeH_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) HeH
         if (HeH .gt. 0.0d0) thereis_He = .true.
         write(*,'(A,F8.5)') '   base.inp: He/H -> ', HeH
      else if (index(line,'Kzz_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) he_kzz
         write(*,'(A,ES9.2,A)') '   base.inp: He_Kzz -> ', he_kzz, ' cm2/s'
      else if (index(line,'q_H2_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) q_h2_base
         write(*,'(A,F8.5)') '   base.inp: q_H2(base) -> ', q_h2_base
      else if (index(line,'p_base') .gt. 0) then
         str = get_word(line,2);  read(str,*) p_base_bar
         write(*,'(A,ES9.2,A)') '   base.inp: p_base -> ', p_base_bar, ' bar'
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
   end subroutine read_base_inp

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
	
	! Initialize counters and strings
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
      
    ! End of module
	end module Read_input

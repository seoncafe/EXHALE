   module  Read_input
   ! Read input planetary parameters adn define

   use global_parameters
   use lower_column, only: q_h2_equilibrium   ! Tier-2a molecular base
   use metals_input        ! optional metals.inp abundance reader
   use charge_exchange, only: cx_init       ! build active charge-exchange set
   use species_table, only: n_melem, iel_C, iel_O, iel_N, iel_Mg,  &
                            iel_Si, iel_Ca, iel_Na, iel_K, iel_S,  &
                            iel_Fe, mion_ethr, melem_i0, melem_A

   implicit none
      
   contains
      
   subroutine input_read
   ! Subroutine to read the input file and assign names and values 
   !	to global constants
	
	character(len = :), allocatable :: str
	character(len = 250) 		    :: line
	integer                         :: ios
	integer                         :: im
      
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

   ! Open file for reading
	write(*,*) '(input_read.f90) Reading the input.inp file..'
   open(unit = 11, file = inp_file)

	! --- Go line by line and read
	
		! Planet name
		read(11,'(A)') line
		p_name = get_word(line, 3)

     		! Log10 of n0
      	read(11,'(A)') line
      	str = get_word(line, 7)
     		read(str,*) n0
    
	     	! Planet radius
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) R0 
	     	
	     	! Planet mass
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) Mp 
	     	
	     	! Equilibrium temperature
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) T0 
	     	
	     	! Orbital distance
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) a_orb 
	     	
	     	! Escape radius
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) r_esc 
	     	
	     	! He/H number ratio
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) HeH 
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

	! Per-element abundances in canonical element order (iel_*), so the
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
		read(11,'(A)') line
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
		read(11,'(A)') line
		str = get_word(line, 5)
	     	read(str,*) Mstar 
	     	
		! Spectrum type 
		read(11,'(A)') line
		sp_type = get_word(line, 3)

		! Next read properties of spectrum
		select case (sp_type)
		
			case ('Load')			! Load from file
				read(11,'(A)') line
				sed_file = get_word(line, 3)
				do_read_sed = .true.
			
			case ('Power-law')
				read(11,'(A)') line
				str = get_word(line, 3)
				read(str,*) PLind 
				is_PL_sed = .true.
	
			case ('Monochromatic')
			
				! Set corresponding logical to true
				is_monochr = .true.
				
				! Read photon enerrgy
				read(11,'(A)') line
				str = get_word(line, 4)
				read(str,*) e_low  
				
				! Remove helium if monochromatic and
				!	photon energy lower than helium ionization threshold
				if (e_low .lt. e_th_HeI) thereis_He = .false.
			     	
		end select
		
		! Only EUV status
		read(11,'(A)') line
		str = get_word(line, 4)
		if (str .eq. 'False') thereis_Xray = .true.
		
		! If not monochromatic, read energy bands
		if (.not. is_monochr ) then 
			
			if (.not.thereis_Xray) then 
			
				! Read e_low
				read(11,'(A)') line
				str = get_word(line, 4)
			     	read(str,*) e_low  
			     	
				! Read e_mid
				str = get_word(line, 6)
				read(str,*) e_mid 
				
				! Set e_top to default
				e_top = 1.24e3 
			     	
			else
				! Read e_low
				read(11,'(A)') line
				str = get_word(line, 4)
			     	read(str,*) e_low  
			     	
				! Read e_mid
				str = get_word(line, 6)
			     	read(str,*) e_mid 
				
				! Read e_top
				str = get_word(line, 8)
			     	read(str,*) e_top 
			     	
			endif
			
		endif
			
		! Read X-ray luminosity if included
		if (thereis_Xray) then
			read(11,'(A)') line
			str = get_word(line, 6)
		     	read(str,*) LX  
		else
			LX = 0.0
		endif
		
		! Read LEUV luminosity
		read(11,'(A)') line
		str = get_word(line, 6)
		read(str,*) LEUV  
		
		! Read grid type
		read(11,'(A)') line
		grid_type = get_word(line, 3)
		
		! Read numerical flux
		read(11,'(A)') line
		flux = get_word(line, 3)
		
		! Read reconstruction scheme
		read(11,'(A)') line
		rec_method = get_word(line, 3)
		if (rec_method.eq.'WENO3') use_weno3 = .true.
		if (rec_method.eq.'PLM')   use_plm = .true.
		
		! Include He23S
		read(11,'(A)') line
		str = get_word(line, 3)
		if (str .eq. 'True')  thereis_HeITR = .true.

		! Remove HeITR chemistry if He is not included
		if (.not. thereis_He) thereis_HeITR = .false.

		! Tier-2 molecular chemistry constraints (v1): requires He;
		! exclusive with trace metals (merged mol+metals = later work item).
		if (thereis_mol .and. .not. thereis_He) then
			write(*,*) '(input_read) ERROR: Molecular chemistry needs He/H>0.'
			stop
		endif
		if (thereis_mol .and. he_diffusion) then
			write(*,*) '(input_read) ERROR: Molecular chemistry + He_diffusion'//&
			           ' not supported yet.'
			stop
		endif
		if (thereis_mol .and. thereis_metals) then
			write(*,*) '(input_read) ERROR: Molecular chemistry + metals '//&
			           'not supported yet (remove metals.inp).'
			stop
		endif

		! IC status
		read(11,'(A)') line
		str = get_word(line, 3)
		if (str .eq. 'True')  do_load_IC = .true.
		
		! Do only post-processing
		read(11,'(A)') line
		str = get_word(line, 4)
		if (str .eq. 'True')  then
			do_only_pp  = .true.
			force_start = .false. ! Set to false to avoid overlap
		endif

		! Force start of sim.
		read(11,'(A)') line
		str = get_word(line, 3)
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
		! Phase 3a excited-H option (also appended after "Force start:" so
		! older files are unaffected): "Stellar Teff [K]: <T>" + "Stellar
		! radius [R_sun]: <R>" supply the diluted-blackbody Balmer continuum
		! that photoionizes/heats H(n=2). The coupling is enabled iff both are
		! given (T_star_eff>0, R_star>0). "Deexc heat: True" additionally turns
		! on the (overlapping) collisional de-excitation heating term.
		spherical_domain = .false.
		r_out_user       = 0.0d0
		T_star_eff       = 0.0d0
		R_star           = 0.0d0
		incl_deexc_heat  = .false.
		jlya_mode        = 0
		transonic_ic     = .false.
		hot_parker_ic    = .false.
		T_wind_ic        = 1.0d4
		use_newton_ieq   = .true.    ! upgraded solvers are the default
		use_brent_tsolve = .true.
		windae_seed_file = 'inputdata/windae_seed.csv'
		windae_seed_out  = ''
		hydrostatic_base = .false.
		base_bc_mode     = 0          ! density-anchored base (legacy) by default
		resid_vol        = .true.     ! volume-weighted residual norm by default
		ates_photoion_rate = .false.  ! default: Verner+1996 He I (1^1S) photoion.
		do
			read(11,'(A)',iostat = ios) line
			if (ios .ne. 0) exit
			if (len_trim(line) .eq. 0) cycle
			if (index(line,'Domain mode') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'Spherical') spherical_domain = .true.
			else if (index(line,'Outer radius') .gt. 0) then
				str = get_word(line, 4)
				read(str,*) r_out_user
			else if (index(line,'Stellar Teff') .gt. 0) then
				str = get_word(line, 4)
				read(str,*) T_star_eff
			else if (index(line,'Stellar radius') .gt. 0) then
				str = get_word(line, 4)
				read(str,*) R_star
			else if (index(line,'Deexc heat') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'True') incl_deexc_heat = .true.
			else if (index(line,'Wind-AE seed out') .gt. 0) then
				windae_seed_out = trim(get_word(line, 4))
			else if (index(line,'Wind-AE seed') .gt. 0) then
				windae_seed_file = trim(get_word(line, 3))
			else if (index(line,'Jlya RT file') .gt. 0) then
				jlya_rt_file = get_word(line, 4)
				jlya_mode    = 1
			else if (index(line,'Jlya escape-prob') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'True') jlya_mode = 2
			else if (index(line,'Stellar Lya flux') .gt. 0) then
				str = get_word(line, 5)
				read(str,*) F_Lya_star
			else if (index(line,'Lya stellar halfwidth') .gt. 0) then
				str = get_word(line, 5)
				read(str,*) dv_star_lya
			else if (index(line,'Lya stellar boost') .gt. 0) then
				str = get_word(line, 5)
				read(str,*) lya_star_boost
			else if (index(line,'du_th') .gt. 0) then
				! Two convergence thresholds: "du_th [PLM,WENO3]: <du_plm> <du_final>"
				! Run PLM until du < du_plm, then switch to WENO3 and converge at
				! du < du_final. If du_plm <= du_final, single-stage at du_final.
				str = get_word(line, 3);  read(str,*) du_th_plm
				str = get_word(line, 4);  read(str,*) du_th
			else if (index(line,'ATES_photoionization_rate') .gt. 0) then
				! Revert He I (1^1S) photoionization to the legacy ATES 2-term fit
				! (default is Verner+1996). "ATES_photoionization_rate: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') ates_photoion_rate = .true.
			else if (index(line,'Molecular chemistry') .gt. 0) then
				! Tier-2 molecular network (docs/lower_atmosphere_*).
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true') thereis_mol = .true.
			else if (index(line,'Molecular base') .gt. 0) then
				! Tier-2a: EOS-only molecular base (docs/lower_atmosphere_*).
				str = get_word(line, 3)
				if (str .eq. 'True' .or. str .eq. 'true') molecular_base = .true.
			else if (index(line,'Lower column') .gt. 0) then
				! Tier-1 analytic lower column: "Lower column: <R_1bar in R_J>"
				str = get_word(line, 3);  read(str,*) lower_col_r1bar
			else if (index(line,'He_Kzz') .gt. 0) then
				! Eddy diffusion coefficient [cm^2/s] for He/H separation.
				! "He_Kzz: 1.0e9"
				str = get_word(line, 2);  read(str,*) he_kzz
			else if (index(line,'He_alphaT') .gt. 0) then
				! Thermal-diffusion factor alpha_T for He (P2c). "He_alphaT: 0.0"
				str = get_word(line, 2);  read(str,*) he_alphaT
			else if (index(line,'He_ambipolar') .gt. 0) then
				! Ambipolar-corrected settling mass (P2b). "He_ambipolar: False"
				str = get_word(line, 2)
				if (str .eq. 'False' .or. str .eq. 'false') he_ambipolar = .false.
			else if (index(line,'He_metal_diffusion') .gt. 0) then
				! Diffuse trace metals too (P2d). "He_metal_diffusion: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') he_metal_diffusion = .true.
			else if (index(line,'He_diffusion') .gt. 0) then
				! Phase-1 He/H diffusive separation (default off).
				! "He_diffusion: True"
				str = get_word(line, 2)
				if (str .eq. 'True' .or. str .eq. 'true') he_diffusion = .true.
			else if (index(line,'Stall') .gt. 0) then
				! Stall-detector override: "Stall [tol,N]: <rel_tol> <N_steps>"
				! (smaller tol and/or larger N = harder to declare a plateau)
				str = get_word(line, 3);  read(str,*) stall_tol
				str = get_word(line, 4);  read(str,*) N_stall
				write(*,'(A,ES9.2,A,I0)') ' (input_read) Stall override: tol =', &
				                          stall_tol, ', N =', N_stall
			else if (index(line,'Energy solver') .gt. 0) then
				! "Energy solver: Explicit" reverts to the original forward-
				! Euler source update (solver-component isolation tests).
				str = get_word(line, 3)
				if (str .eq. 'Explicit') then
					use_semi_implicit_energy = .false.
					write(*,*) '(input_read) Energy solver: explicit forward Euler'
				endif
			else if (index(line,'Time stepping') .gt. 0) then
				! "Time stepping: Local" = per-cell pseudo-time steps
				! (steady-state convergence acceleration; not time-accurate).
				str = get_word(line, 3)
				if (str .eq. 'Local') then
					use_local_dt = .true.
					write(*,*) '(input_read) Time stepping: local (per-cell) pseudo-dt'
				endif
			else if (index(line,'Level tol') .gt. 0) then
				! "Level tol: <val>" overrides the mass-flux level-stability
				! tolerance (<= 0 disables the level gate; legacy stops).
				str = get_word(line, 3);  read(str,*) lev_th
				write(*,'(A,ES9.2)') ' (input_read) Level-stability tol =', lev_th
			else if (index(line,'Solver') .gt. 0) then
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
			else if (index(line,'Valve eps') .gt. 0) then
				! "Valve eps: <v_eps>" smooths the base one-way valve
				! (softplus; <= 0 keeps the exact legacy max(v,0)).
				str = get_word(line, 3);  read(str,*) valve_eps
				write(*,'(A,ES9.2)') ' (input_read) Smooth base valve, eps =', valve_eps
			else if (index(line,'Hydrostatic base') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'True') hydrostatic_base = .true.
				if (hydrostatic_base) write(*,'(A)') ' (input_read) '//   &
				   'Hydrostatic base ghost cells enabled'
			else if (index(line,'Shapiro filter') .gt. 0) then
				str = get_word(line, 3);  read(str,*) shapiro_eps
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) shapiro_every
				if (shapiro_eps .gt. 0.0d0) write(*,'(A,ES9.2,A,I0,A)')  &
				   ' (input_read) Shapiro filter eps =', shapiro_eps,    &
				   ', every ', shapiro_every, ' steps'
			else if (index(line,'Base BC') .gt. 0) then
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
			else if (index(line,'Base velocity') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'valve')    base_v_massflux = .false.
				if (str .eq. 'massflux') base_v_massflux = .true.
				if (base_v_massflux) then
					write(*,'(A)') ' (input_read) Base velocity from '//  &
					   'mass-flux F_c (CETIMB-style)'
				else
					write(*,'(A)') ' (input_read) Base velocity: legacy valve'
				endif
			else if (index(line,'Viscosity') .gt. 0) then
				str = get_word(line, 3);  read(str,*) visc_mu0
				str = get_word(line, 4)
				if (len_trim(str) .gt. 0) read(str,*) visc_s
				if (visc_mu0 .gt. 0.0d0) write(*,'(A,ES9.2,A,F5.2,A)')  &
				   ' (input_read) Viscosity (Phase-1) mu0 =', visc_mu0,  &
				   ', s =', visc_s, '  [un-validated]'
			else if (index(line,'Resid tol') .gt. 0) then
				! "Resid tol: <val>" = converge on the steady residual ||R||
				! instead of du (<= 0 disables; legacy du-based stop).
				str = get_word(line, 3);  read(str,*) resid_th
				write(*,'(A,ES9.2)') ' (input_read) Residual-based convergence, tol =', resid_th
			else if (index(line,'Resid norm') .gt. 0) then
				! "Resid norm: vol|Linf" -- residual norm for convergence.
				! vol (default) = volume-weighted; Linf = legacy max-over-cells.
				str = get_word(line, 3)
				if (str .eq. 'Linf' .or. str .eq. 'linf' .or. str .eq. 'LINF') &
					resid_vol = .false.
				if (str .eq. 'vol' .or. str .eq. 'volume') resid_vol = .true.
				write(*,'(A,L1)') ' (input_read) Volume-weighted residual norm: ', resid_vol
			else if (index(line,'CFL') .gt. 0) then
				! Override the CFL number ("CFL: <value>"); lower = smaller dt.
				str = get_word(line, 2);  read(str,*) CFL
			else if (index(line,'Transonic IC') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'True') transonic_ic = .true.
			else if (index(line,'Hot Parker IC') .gt. 0) then
				str = get_word(line, 4)
				read(str,*) T_wind_ic
				hot_parker_ic = .true.
			else if (index(line,'IC mode') .gt. 0) then
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
			else if (index(line,'Newton solver') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'False') use_newton_ieq = .false.
			else if (index(line,'Brent solver') .gt. 0) then
				str = get_word(line, 3)
				if (str .eq. 'False') use_brent_tsolve = .false.
			endif
		enddo

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
			stop 1
		endif

   close(unit = 1)
	write(*,*) '(input_read.f90) Done'

   !------ Definition of physical parameters ------!
      
   n0     = 10.0**(n0)
   ! ---- Tier-3 optional base.inp (written by src/utils/run_lower.py or a
   ! lower-atmosphere model): overrides the base temperature, base radius
   ! [R_J], He/H ratio and eddy K_zz BEFORE the derived constants below.
   ! Absent file = no-op (byte-identical legacy).
   call read_base_inp

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
         stop
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

   ! Composition factors. mass_per_H = gas mass per H nucleus [m_H];
   ! ntot_bc = total nuclei density at the base in units of n0 (n0 keeps
   ! its H/He-nuclei meaning, so n_H = n0/(1+HeH) is unchanged). With
   ! eos_metals 1 (default) the trace metals contribute their mass and
   ! their nuclei; with eos_metals 0 (or no metals) both reduce to the
   ! legacy H/He-only values (mass_per_H = 1+4*HeH, ntot_bc = 1).
   mass_per_H = 1.0 + 4.0*HeH
   ntot_bc    = 1.0
   if (eos_include_metals .and. thereis_metals) then
      mass_per_H = mass_per_H + sum(melem_ab*melem_A)
      ntot_bc    = (1.0 + HeH + sum(melem_ab))/(1.0 + HeH)
   endif

   rho_bc = mass_per_H/(1.0 + HeH)

   ! Tier-2a passive molecular base (docs/lower_atmosphere_coupling.*):
   ! remove from the base particle budget the H nuclei bound into H2 at
   ! (1 ubar, T0) according to the chemical-equilibrium fit; per n0
   ! (H+He nuclei) that is (x2/2)/(1+HeH) particles.  Lowers the base
   ! pressure / raises the base mean molecular weight.  EOS-only: the
   ! species arrays stay atomic (H2 chemistry is Tier-2 proper).
   if (molecular_base) then
      block
         real*8 :: qmb, x2mb
         qmb  = q_h2_equilibrium(1.0d-6, T0)
         x2mb = 2.0d0*qmb*(1.0d0 + HeH)/(1.0d0 + qmb)
         if (x2mb .gt. 1.0d0) x2mb = 1.0d0
         ntot_bc = ntot_bc - 0.5d0*x2mb/(1.0d0 + HeH)
         write(*,'(A,F6.3,A,F6.3)') ' (input_read) Molecular base: '//   &
            'q_H2(1ubar,T0) =', qmb, ' -> ntot_bc =', ntot_bc
      end block
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

	! Tier-2 molecular system: H+/He+/He++ + H2/H2+/H3+/HeH+ (+ He 2^3S)
	if (thereis_mol) then
		N_eq = 7
		if (thereis_HeITR) N_eq = 8
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
   end subroutine input_read

   ! ------------------------------------------------------------------- !

   subroutine read_base_inp
   ! Tier-3 lower-atmosphere handoff file (optional).  Keyword lines:
   !   T_base    <K>      -> overrides T0 (base temperature)
   !   r_base    <R_J>    -> overrides the "Planet radius" (1-ubar radius)
   !   HeH_base  <ratio>  -> overrides the He/H number ratio
   !   Kzz_base  <cm2/s>  -> sets he_kzz (used by He_diffusion)
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
      endif
   enddo
   close(ub)
   end subroutine read_base_inp
      
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

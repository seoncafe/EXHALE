	module setup_report
	! Contains subroutine to write the setup of the current 
	!	simulation to file
	
	use global_parameters
	use charge_exchange, only: he_h_charge_exchange, cx_o2p_h_scale
	use species_table,   only: melem_name
	use IC_load,         only: melem_from_abundance
	use composition,     only: h2_mixing_ratio_base, h2_mixing_ratio_ceiling
	use lower_atmosphere_profile, only: lap_in_use, lap_report_provenance, &
	                     lap_flux_measured, lap_flux_window_empty,        &
	                     lap_flux_nface, lap_FH_median, lap_FH_spread,    &
	                     lap_FHe_median, lap_FHe_spread,                  &
	                     lap_Mdot_median, lap_Mdot_spread,                &
	                     lap_flux_r_lo_Rp, lap_flux_r_hi_Rp,              &
	                     lap_flux_r_lo_measured,                          &
	                     lap_steady_r_lo_Rp, lap_steady_nface,            &
	                     lap_steady_FH_median,  lap_steady_FH_spread,     &
	                     lap_steady_FHe_median, lap_steady_FHe_spread,    &
	                     lap_steady_Mdot_median, lap_steady_Mdot_spread

	implicit none
	contains
	
	! ------------------------------ !
	
	subroutine write_setup_report
	! Write summary of the current simulation setup to file

	! Number of metal elements with a non-zero abundance. NOTE: the intrinsic
	! COUNT cannot be used here -- `count` is the global step counter in
	! global_parameters, so the name is shadowed module-wide.
	integer :: imet, n_met_active

	write(*,*) '(write_setup_report.f90) Writing the setup report on EXHALE_setup.out..'

	n_met_active = 0
	if (thereis_metals) then
		do imet = 1,size(melem_ab)
			if (melem_ab(imet) .gt. 0.0d0) n_met_active = n_met_active + 1
		enddo
	endif
	
	write(outfile,*) '######## Simulation for ', p_name, ' ########'
	write(outfile,*) ' '
	write(outfile,*) ' ----- Planetary parameters ----- '
	write(outfile,*) ' '
	write(outfile,1)  &
      ' - Planet mass: ', Mp/MJ, ' [M_J], ', Mp/M_earth, ' [M_earth]'
	write(outfile,2)  &
      ' - Planet radius: ', R0/RJ, ' [R_J], ', R0/R_earth, ' [R_earth]'
	write(outfile,3) ' - Orbital distance: ', a_orb/AU, ' [AU]'
	write(outfile,4) ' - Equilibrium temperature: ', T0, ' [K]'
	write(outfile,5) ' - Jeans parameter (beta_0): ', b0
	write(outfile,15) & 
      ' - Log surface gravitational potential: ', log10(Gc*Mp/R0), ' erg/g'	 
	write(outfile,*) ' '
	write(outfile,*) ' ----- Stellar parameters ----- '
	write(outfile,*) ' '
	write(outfile,6) ' - Star mass: ', Mstar/Msun, ' [M_sun]'
	write(outfile,7) ' - Log EUV luminosity: ', LEUV, ' [erg/s]'
	write(outfile,8) ' - Log X-ray luminosity: ', LX, ' [erg/s]'
	write(outfile,9)  &
      ' - Log XUV flux at the planet distance: ', log10(J_XUV), ' [erg/(cm^2 s)]'
	write(outfile,*)	
	write(outfile,*) ' ----- Simulation setup parameters -----'
	write(outfile,*)
	write(outfile,10) &
      ' - Upper boundary of the domain: ', r_max, ' [R_p]'
	if (spherical_domain) then
		write(outfile,*) &
         '- Domain mode: Spherical (pure planetary potential, '          // &
         'tidal/centrifugal dropped, extended to user outer radius)'
	else
		write(outfile,*) &
         '- Domain mode: Roche (full tidal potential, truncated at '     // &
         'the Hill/L1 radius)'
	endif
	if (thereis_He) then
		write(outfile,11) & 
         ' - Simulating a H/He atmosphere with',' He/H ratio of ', HeH
		if (thereis_HeITR) 	&
			write(outfile,*) '- Including helium triplet chemistry'
	else
		write(outfile,*) '- Simulating a pure H atmosphere'
	endif
	! Which species blocks the coupled equilibrium system carries, and how
	! large that makes it. Molecules and metals are solved together when both
	! are present (System_HeH_mol_metals), so say so here rather than leaving
	! the reader to infer it from N_eq.
	if (thereis_mol) write(outfile,*) &
      '- Including molecular chemistry (H2, H2+, H3+, HeH+)'
	if (thereis_mol) then
		if (F_LW_star .gt. 0.0d0) then
			write(outfile,'(A,ES10.3,A)') &
      ' - H2 Lyman-Werner photodissociation on: band flux ', F_LW_star,   &
      ' erg cm^-2 s^-1 (912-1110 A, at the planet)'
		else
			write(outfile,*) &
      '- H2 Lyman-Werner photodissociation off ("Stellar LW flux" unset)'
		endif
	endif
	if (thereis_oxychem) then
		write(outfile,*) &
      '- Including oxygen chemistry (OH, H2O, CO) in the same system'
		if (oxygen_transport) then
			write(outfile,*) &
      '- Molecular carriers (H2, OH, H2O, CO) are TRANSPORTED: implicit'// &
      ' diffusion-advection solved with their chemistry'
			if (maxval(kzz_cell) .le. 0.0d0) write(outfile,*) &
      '  WARNING K_zz = 0 everywhere, so the transport is pure molecular'//&
      ' diffusion -- the wrong limit for a lower atmosphere,'
			if (maxval(kzz_cell) .le. 0.0d0) write(outfile,*) &
      '  where eddy mixing is what holds the composition well mixed'//     &
      ' below the homopause. Set "He_Kzz:" or a profile.'
		else
			write(outfile,*) &
      '- Molecular carriers are a LOCAL steady state ("Oxygen transport:'//&
      ' False"): the chemistry alone, not a model of a base'
		endif
		write(outfile,'(A,5ES10.3)') &
      ' - FUV band fluxes at the planet [erg cm^-2 s^-1]'//               &
      ' LW/B1/B2(Lya)/B3/B4: ', F_LW_star, F_FUV_B1, F_Lya_star,          &
      F_FUV_B3, F_FUV_B4
		write(outfile,*) &
      '- The 912-1110 A band is the Lyman-Werner interval and carries'//  &
      ' "Stellar LW flux": H2, H2O and OH share that beam'
		write(outfile,*) &
      '- B2 is the INCIDENT stellar Ly-alpha flux: only the H2O/OH'//     &
      ' continuum attenuates it, not H I resonance scattering,'
		write(outfile,*) &
      '  so the Ly-alpha photolysis rate is an upper bound'//             &
      ' (see output/FUV_bands.txt)'
		write(outfile,*) &
      '- Base H2/H partition: COMPUTED by the oxygen chemistry'//         &
      ' (q_H2_base is refused)'
	endif
	if (thereis_metals) write(outfile,'(A,I0,A)') &
      ' - Including ', n_met_active, ' trace metal element(s)'
	if (thereis_mol .and. thereis_metals) write(outfile,*) &
      '- Molecules and metals solved in one system (shared electron density)'
	write(outfile,'(A,I0)') &
      ' - Coupled equilibrium system size: N_eq = ', N_eq
	if (do_read_sed) &
		write(outfile,*) & 
         ' - Spectrum read from external file: ', sed_file
	if (is_PL_sed) &
      	write(outfile,12)  & 
            ' - Using power-law spectrum ', 'with index ', PLind
	if (is_monochr) & 
		write(outfile,13)  &
         ' - Using monochromatic radiation with energy ', e_low
	if (appx_mth.eq.'alpha') then 
		write(outfile,14) ' - 2D approximation used: alpha =, with alpha = ',a_tau
	else
		write(outfile,*) '- 2D approximation used: ', appx_mth
	endif
	write(outfile,*) 
	write(outfile,*) '----- Numerical parameters -----'	
	write(outfile,*)
	write(outfile,*) '- Grid type: ', grid_type
	write(outfile,'(A,I0,A,I0,A)') ' - Grid cells: ', N,                    &
		' computational cells (+ ', Ng, ' ghost cells on each side)'
	if (grid_type .eq. 'Mixed') then
		write(outfile,16) '- Base grid: ', N_low_cells,                    &
			' uniform cells of ', dr_base, ' R_p (uniform region ',        &
			dr_base*N_low_cells, ' R_p)'
	else
		write(outfile,*) '- Base grid: "Base grid [dr,cells]" is ignored'//&
			' by grid type '//trim(grid_type)//' (Mixed only)'
	endif
	! Resolution of the base density scale height H = kT_eq/(mu g) in cells:
	! b0 = R_p/H(T_eq) is the Jeans parameter, dr_j(1) the first cell size
	! after the Mixed-grid smoothing. Below ~5 cells per H the scheme carries
	! an undamped stationary 2*dr entropy mode (docs/hd189_base_checkerboard.md).
	write(outfile,17) '- Base scale-height resolution: H(T_eq)/dr = ',      &
		1.0d0/(b0*dr_j(1)), ' cells'
	if (1.0d0/(b0*dr_j(1)) .lt. 10.0d0)                                     &
		write(outfile,*) '  WARNING: base scale height spans < 10 cells;'//&
			' expect a stationary cell-to-cell entropy mode at the base.'
	write(outfile,*) '- Numerical flux: ', flux
	write(outfile,*) '- Reconstruction method: ', rec_method
	! Artificial dissipation of the 2*dr contact/entropy mode that the
	! contact-resolving upwind flux stops damping as v -> 0. It enters the
	! numerical flux, so the marching loop and the steady residual see the
	! same equation; the Shapiro filter above does not.
	if (lowmach_damp_eps .gt. 0.0d0) then
		write(outfile,20) '- Low-Mach contact-mode damping: eps4 = ',        &
			lowmach_damp_eps, ', gate closes at M = ', lowmach_damp_mach_th
		if (16.0d0*lowmach_damp_eps*CFL .gt. 1.0d0)                          &
			write(outfile,*) '  WARNING: eps4 exceeds the explicit'//        &
				' stability bound 1/(16 CFL).'
	else
		write(outfile,*) '- Low-Mach contact-mode damping: off'
	endif
	! Lower boundary: which closure sets the ghost pressure, and the hard cap
	! on marching steps.
	if (hydrostatic_base) then
		write(outfile,*) '- Base ghost pressure: interior gradient '//       &
			'extrapolation (Hydrostatic base)'
	else if (base_ghost_T_continuous) then
		write(outfile,*) '- Base ghost temperature: continuous '//           &
			'(dT/dr = 0, T_ghost = T_1)'
	else
		write(outfile,*) '- Base ghost temperature: isothermal (T_ghost = T0)'
	endif
	write(outfile,18) '- Coronal cooling cutoff width: w = ',                &
		coronal_cutoff_width
	if (base_ir_field) then
		write(outfile,*) '- Base IR field: on (fine-structure lines see '//   &
			'a diluted B_nu(T0) from the lower atmosphere)'
	else
		write(outfile,*) '- Base IR field: off (fine-structure lines '//      &
			'emit into vacuum)'
	endif
	if (mol_ir_bands) then
		write(outfile,*) '- Molecular IR bands: on (H2 quadrupole and '//     &
			'magnetic dipole lines, H2O and CO bands, in LTE, exchanging'
		write(outfile,*) '    with the same diluted B_nu(T0); each stops '//  &
			'at its own radiative equilibrium temperature)'
	else
		write(outfile,*) '- Molecular IR bands: off (no H2, H2O or CO '//     &
			'infrared channel below the H2 -> H front)'
	endif
	write(outfile,19) '- Max marching steps: ', count_max
	write(outfile,*) 
	if (.not.do_load_IC) &
		write(outfile,*) '----- Starting a new simulation ----- '
	if (do_load_IC) 		&
		write(outfile,*) '----- Continuing existing simulation ----- '
	if (do_load_IC .and. allocated(melem_from_abundance)) then
		do imet = 1,size(melem_from_abundance)
			if (melem_from_abundance(imet))                                 &
				write(outfile,*) '- WARNING: element ',                      &
					trim(melem_name(imet)), ' was absent (or identically '// &
					'zero) in the restart file and was initialized as '//    &
					'neutral at the input abundance'
		enddo
	endif
	if (.not.do_load_IC) then
		! IC family actually in effect (this report is written after
		! set_IC, so an "IC mode: auto" selection has already run).
		if (hot_parker_ic) then
			write(outfile,*) '- IC: hot-Parker warm seed'
		else if (transonic_ic) then
			write(outfile,*) '- IC: transonic isothermal wind'
		else
			write(outfile,*) '- IC: cold hydrostatic'
		endif
		if (ic_mode .eq. 3) &
			write(outfile,*) '  (chosen automatically: IC mode = auto)'
	endif
	if (do_only_pp)		&
		write(outfile,*) '- Evaluating post processing only'
	if (force_start) &
		write(outfile,*) '- Forcing the simulation for the first 1000 steps'
	   
1	format(A16,F7.3,A8,F9.2,A10)
2	format(A18,F7.3,A8,F9.2,A10)	
3	format(A21,F6.4,A5)	
4	format(A28,F6.1,A4)	
5	format(A29,F7.2)
6	format(A14,F5.3,A8)	
7	format(A23,F6.3,A8)	
8	format(A25,F6.3,A8)	
9	format(A40,F5.3,A15)	
10	format(A33,F5.2,A6)	
11	format(A36,A15,ES11.4)		
12	format(A28,A11,F5.2)	
13 format(A45,F7.2)	
14	format(A48,F6.3)		
15 format(A40,F5.2,A6)
16 format(A14,I4,A19,ES9.2,A22,ES9.2,A6)
17 format(A46,F8.1,A7)
18 format(A,F6.3)
19 format(A,I0)
20 format(A,ES9.2,A,ES9.2)

	write(*,*) '(write_setup_report.f90) Done.'

	end subroutine write_setup_report

	! ------------------------------ !

	subroutine write_parse_dump
	! Dump every variable input_read derives from input.inp (plus any base.inp
	! override) to parse_dump.txt, one "name = value" line per variable, in the
	! docs/input_schema.md key order. Gated by EXHALE_PARSE_DUMP=1 in
	! EXHALE_main and used by the parser-refactor regression corpus
	! (backup/regression/run_parse_corpus.sh). metals.inp / opacity.inp
	! variables are out of scope (separate parsers). Values are the final
	! post-input_read state: n0, R0, Mp, a_orb, Mstar in cgs, and T0/R0/HeH/
	! he_kzz after any base.inp override, so the dump captures the whole
	! derived-parameter chain deterministically.
	! Formats: reals ES23.15E3, integers plain, logicals T/F, strings trimmed.
	integer :: u

	open(newunit=u, file='parse_dump.txt', status='replace', action='write')

	write(u,'(A)') '# EXHALE parse dump (EXHALE_PARSE_DUMP=1): variables set by input_read'
	write(u,'(A)') '# order follows docs/input_schema.md; cgs where input_read converts'

	! ----- core block (fixed order) -----
	call put_s('p_name', p_name)
	call put_r('n0', n0)
	call put_r('R0', R0)
	call put_r('Mp', Mp)
	call put_r('T0', T0)
	call put_r('a_orb', a_orb)
	call put_r('r_esc', r_esc)
	call put_r('HeH', HeH)
	call put_l('thereis_He', thereis_He)
	call put_s('appx_mth', appx_mth)
	call put_r('a_tau', a_tau)
	call put_r('Mstar', Mstar)
	call put_s('sp_type', sp_type)
	if (do_read_sed) then
	   call put_s('sed_file', sed_file)
	else
	   call put_s('sed_file', '(unset)')
	endif
	call put_l('do_read_sed', do_read_sed)
	if (is_PL_sed) then
	   call put_r('PLind', PLind)
	else
	   call put_s('PLind', '(unset)')
	endif
	call put_l('is_PL_sed', is_PL_sed)
	call put_l('is_monochr', is_monochr)
	call put_r('e_low', e_low)
	call put_l('thereis_Xray', thereis_Xray)
	if (.not. is_monochr) then
	   call put_r('e_mid', e_mid)
	   call put_r('e_top', e_top)
	else
	   call put_s('e_mid', '(unset)')
	   call put_s('e_top', '(unset)')
	endif
	call put_r('LX', LX)
	call put_r('LEUV', LEUV)
	call put_s('grid_type', grid_type)
	call put_i('N', N)
	call put_r('dr_base', dr_base)
	call put_i('N_low_cells', N_low_cells)
	call put_s('flux', flux)
	call put_s('rec_method', rec_method)
	call put_l('use_weno3', use_weno3)
	call put_l('use_plm', use_plm)
	call put_l('recon_two_stage', recon_two_stage)
	call put_l('thereis_HeITR', thereis_HeITR)
	call put_l('do_load_IC', do_load_IC)
	call put_l('do_only_pp', do_only_pp)
	call put_l('force_start', force_start)

	! ----- keyword-extension block -----
	call put_l('spherical_domain', spherical_domain)
	call put_r('r_out_user', r_out_user)
	call put_r('T_star_eff', T_star_eff)
	call put_r('R_star', R_star)
	call put_l('incl_deexc_heat', incl_deexc_heat)
	call put_s('windae_seed_out', windae_seed_out)
	call put_s('windae_seed_file', windae_seed_file)
	call put_s('jlya_rt_file', jlya_rt_file)
	call put_i('jlya_mode', jlya_mode)
	call put_r('F_Lya_star', F_Lya_star)
	call put_r('dv_star_lya', dv_star_lya)
	call put_r('lya_star_boost', lya_star_boost)
	call put_l('lya_bottom_absorber', lya_bottom_absorber)
	call put_r('du_th', du_th)
	call put_r('du_th_plm', du_th_plm)
	call put_l('ates_photoion_rate', ates_photoion_rate)
	call put_l('legacy_hhe_rates', legacy_hhe_rates)
	call put_l('use_sec_ion', use_sec_ion)
	call put_l('use_he_rec_coupling', use_he_rec_coupling)
	call put_l('he_h_charge_exchange', he_h_charge_exchange)
	call put_r('cx_o2p_h_scale', cx_o2p_h_scale)
	call put_l('thereis_mol', thereis_mol)
	call put_r('F_LW_star', F_LW_star)
	call put_l('thereis_oxychem', thereis_oxychem)
	call put_l('oxygen_transport', oxygen_transport)
	call put_r('F_FUV_B1', F_FUV_B1)
	call put_r('F_FUV_B3', F_FUV_B3)
	call put_r('F_FUV_B4', F_FUV_B4)
	call put_l('molecular_base', molecular_base)
	call put_r('q_h2_base', q_h2_base)
	call put_r('p_base_bar', p_base_bar)
	call put_i('lower_atm_mode', lower_atm_mode)
	call put_r('lower_atm_r1bar', lower_atm_r1bar)
	call put_r('lower_col_r1bar', lower_col_r1bar)
	call put_r('he_kzz', he_kzz)
	call put_r('he_alphaT', he_alphaT)
	call put_l('he_ambipolar', he_ambipolar)
	call put_l('he_metal_diffusion', he_metal_diffusion)
	call put_l('he_diffusion', he_diffusion)
	call put_r('stall_tol', stall_tol)
	call put_i('N_stall', N_stall)
	call put_l('use_semi_implicit_energy', use_semi_implicit_energy)
	call put_l('use_local_dt', use_local_dt)
	call put_r('lev_th', lev_th)
	call put_l('use_newton_solver', use_newton_solver)
	call put_r('newton_du_switch', newton_du_switch)
	call put_r('valve_eps', valve_eps)
	call put_l('hydrostatic_base', hydrostatic_base)
	call put_r('shapiro_eps', shapiro_eps)
	call put_i('shapiro_every', shapiro_every)
	call put_r('lowmach_damp_eps', lowmach_damp_eps)
	call put_r('lowmach_damp_mach_th', lowmach_damp_mach_th)
	call put_l('base_ghost_T_continuous', base_ghost_T_continuous)
	call put_i('count_max', count_max)
	call put_r('coronal_cutoff_width', coronal_cutoff_width)
	call put_l('base_ir_field', base_ir_field)
	call put_l('mol_ir_bands', mol_ir_bands)
	call put_i('base_bc_mode', base_bc_mode)
	call put_r('base_p_ubar', base_p_ubar)
	call put_l('base_v_massflux', base_v_massflux)
	call put_l('visc_on', visc_on)
	call put_l('cond_on', cond_on)
	call put_r('visc_mu0', visc_mu0)
	call put_r('visc_s', visc_s)
	call put_r('resid_th', resid_th)
	call put_l('resid_vol', resid_vol)
	call put_r('CFL', CFL)
	call put_l('transonic_ic', transonic_ic)
	call put_r('T_wind_ic', T_wind_ic)
	call put_l('hot_parker_ic', hot_parker_ic)
	call put_i('ic_mode', ic_mode)
	call put_l('use_newton_ieq', use_newton_ieq)
	call put_l('use_brent_tsolve', use_brent_tsolve)

	! ----- derived flags / normalization constants -----
	call put_l('use_excited_H', use_excited_H)
	call put_r('Mrapp', Mrapp)
	call put_r('atilde', atilde)
	call put_r('r_max', r_max)
	call put_r('mass_per_H', mass_per_H)
	call put_r('ntot_bc', ntot_bc)
	call put_r('rho_bc', rho_bc)
	call put_r('v0', v0)
	call put_r('t_s', t_s)
	call put_r('p0', p0)
	call put_r('q0', q0)
	call put_r('b0', b0)
	call put_r('dp_bc', dp_bc)
	call put_i('N_eq', N_eq)
	call put_i('lwa', lwa)

	close(u)
	write(*,*) '(write_parse_dump) wrote parse_dump.txt'

	contains

	subroutine put_r(name, val)   ! real*8
	character(len=*), intent(in) :: name
	real*8,           intent(in) :: val
	write(u,'(A,ES23.15E3)') name//' = ', val
	end subroutine put_r

	subroutine put_i(name, val)   ! integer
	character(len=*), intent(in) :: name
	integer,          intent(in) :: val
	write(u,'(A,I0)') name//' = ', val
	end subroutine put_i

	subroutine put_l(name, val)   ! logical
	character(len=*), intent(in) :: name
	logical,          intent(in) :: val
	write(u,'(A,L1)') name//' = ', val
	end subroutine put_l

	subroutine put_s(name, val)   ! string
	character(len=*), intent(in) :: name
	character(len=*), intent(in) :: val
	write(u,'(A)') name//' = '//trim(val)
	end subroutine put_s

	end subroutine write_parse_dump

	! ------------------------------ !

	subroutine write_resolved_config
	! Machine-readable record of the configuration the wind actually uses,
	! written once the input.inp values and the optional base.inp handoff
	! overrides (r_base/T_base/HeH_base) are resolved.  Consumers
	! (EXHALE_transit.py) read these values instead of re-parsing input.inp,
	! so a base.inp override of the radius/temperature/He ratio reaches the
	! transit geometry as well (docs/lhs1140b_lower_atmosphere_plan_new.md,
	! Phase B).  Format: '# ' comments, then one 'key  value' pair per line.
	integer :: u, ie
	logical :: base_present

	inquire(file='base.inp', exist=base_present)
	open(newunit=u, file='EXHALE_resolved.out', status='replace',        &
	     action='write')
	write(u,'(A)') '# EXHALE resolved configuration (machine-readable).'
	write(u,'(A)') '# Values in effect after input.inp + base.inp resolution;'
	write(u,'(A)') '# these are what the wind solver uses, and what'
	write(u,'(A)') '# EXHALE_transit.py should use instead of input.inp.'
	write(u,'(A)') '# Reals are written to full double precision so a budget'
	write(u,'(A)') '# check can use them without a round-off floor.'
	write(u,'(A,ES23.15E3)') 'planet_radius_RJ          ', R0/RJ
	write(u,'(A,ES23.15E3)') 'planet_mass_MJ            ', Mp/MJ
	write(u,'(A,ES23.15E3)') 'equilibrium_temperature_K ', T0
	write(u,'(A,ES23.15E3)') 'HeH_number_ratio          ', HeH
	write(u,'(A,ES23.15E3)') 'orbital_distance_AU       ', a_orb/AU
	write(u,'(A,ES23.15E3)') 'star_mass_Msun            ', Mstar/Msun
	write(u,'(A,L1)')     'base_inp_present          ', base_present
	! With element diffusion on, He/H is a solved profile and the number
	! above is the reservoir the base is held at, not a column invariant.
	! The budget check has to know which of the two it is testing.
	write(u,'(A,L1)')     'he_diffusion              ', he_diffusion
	! Oxygen chemistry (the A2 option). oxygen_chemistry says whether the
	! O I column means FREE ATOMIC oxygen (it does when this is T) and
	! whether the OH / H2O / CO columns of Ion_species.txt exist;
	! oxygen_base_partition is the provenance of the base H2/H partition,
	! which is the whole point of the option. The five band fluxes are the
	! photon input the oxygen photochemistry actually ran on; the first of
	! them, fuv_band_LW_flux, is the 912-1110 A interval shared with the H2
	! Lyman-Werner absorber.
	write(u,'(A,L1)')     'mol_ir_bands              ', mol_ir_bands
	write(u,'(A,L1)')     'oxygen_chemistry          ', thereis_oxychem
	if (thereis_oxychem) then
		write(u,'(A)')     'oxygen_reaction_set       a2_v1'
		write(u,'(A,L1)')  'oxygen_transport          ', oxygen_transport
		write(u,'(A,ES23.15E3)') 'fuv_band_B1_flux          ', F_FUV_B1
		write(u,'(A,ES23.15E3)') 'fuv_band_Lya_flux         ', F_Lya_star
		write(u,'(A,ES23.15E3)') 'fuv_band_B3_flux          ', F_FUV_B3
		write(u,'(A,ES23.15E3)') 'fuv_band_B4_flux          ', F_FUV_B4
		write(u,'(A,ES23.15E3)') 'fuv_band_LW_flux          ', F_LW_star
		write(u,'(A)')     'oxygen_base_partition     computed'
	else if (q_h2_base .gt. 0.0d0) then
		write(u,'(A)')     'oxygen_base_partition     handoff'
	else if (thereis_mol) then
		write(u,'(A)')     'oxygen_base_partition     equilibrium_fit'
	endif
	! Resolved elemental reservoirs and the EOS factors built from them, so
	! the element-budget check (src/utils/element_budget.py) can compare the
	! solved profiles against the abundances the run actually used, whether
	! they came from metals.inp or from the "<El>_H_base" handoff keys.
	write(u,'(A,ES23.15E3)') 'mass_per_H_amu            ', mass_per_H
	write(u,'(A,ES23.15E3)') 'ntot_bc_per_H             ', ntot_bc
	! The base H2 fraction the run resolved, the ceiling the element ratio
	! allows, and what the pair implies for the hydrogen nuclei. The mixing
	! ratio alone does not say how molecular the base is -- the same q_H2
	! means different things at different He/H -- so the implied fraction of
	! H nuclei bound into H2 is written beside it. Startup refuses a
	! requested value above the ceiling, so these two always satisfy
	! q_H2 <= ceiling; they are recorded because the margin between them is
	! what a lower-atmosphere handoff has to be read against.
	if (molecular_base .or. q_h2_base .gt. 0.0d0) then
		write(u,'(A,ES23.15E3)') 'q_H2_base_resolved        ',              &
			h2_mixing_ratio_base()
		write(u,'(A,ES23.15E3)') 'q_H2_base_ceiling         ',              &
			h2_mixing_ratio_ceiling()
		write(u,'(A,ES23.15E3)') 'H_nuclei_in_H2_fraction   ',              &
			2.0d0*h2_mixing_ratio_base()*(1.0d0 + HeH)                     &
			/(1.0d0 + h2_mixing_ratio_base())
	endif
	do ie = 1, size(melem_ab)
		write(u,'(A,A,A,ES23.15E3)') 'abundance_',                          &
			trim(melem_name(ie)), repeat(' ', 15 - len_trim(melem_name(ie))), &
			melem_ab(ie)
	enddo
	! Provenance of the lower-atmosphere profile, if one is in use, so the
	! closure driver and element_budget.py read one authority for which
	! solution the wind was built on.
	call lap_report_provenance(u)
	! The elemental fluxes measured over the overlap window (section 3.4 of
	! docs/phase_e_flux_closure_design.md). Reported only when a diffusion
	! step has actually produced them: this routine also runs before the wind,
	! and an unmeasured flux must say so rather than print a zero.
	if (lap_in_use) then
		if (.not. lap_flux_measured) then
			write(u,'(A)') 'lower_profile_flux_state  unmeasured'
		else
			! (a) the overlap window: both edges, whether the lower one was
			! measured or fell back, and the reduction over it.
			write(u,'(A,ES23.15E3)')  'lower_profile_flux_r_lo_Rp ',        &
				lap_flux_r_lo_Rp
			if (lap_flux_r_lo_measured) then
				write(u,'(A)') 'lower_profile_flux_r_lo_source spread_rule'
			else
				write(u,'(A)') 'lower_profile_flux_r_lo_source default_1.02'
			endif
			write(u,'(A,ES23.15E3)')  'lower_profile_flux_r_hi_Rp ',        &
				lap_flux_r_hi_Rp
			if (lap_flux_window_empty) then
				write(u,'(A)') 'lower_profile_flux_state  window_empty'
				write(u,'(A,I0)') 'lower_profile_flux_nface  ', 0
				write(u,'(A)') '# The overlap window is empty: the profile'//&
					' stops at or below the radius'
				write(u,'(A)') '# where the base sound wave leaves the'//    &
					' mass flux flat, so no face of the'
				write(u,'(A)') '# interval both models describe carries a'// &
					' usable elemental flux.  The'
				write(u,'(A)') '# steady_* window below is then the only'//  &
					' measurement of the handoff flux.'
			else
				write(u,'(A)') 'lower_profile_flux_state  measured'
				write(u,'(A,I0)')         'lower_profile_flux_nface  ',     &
					lap_flux_nface
				write(u,'(A,ES23.15E3)')  'lower_profile_F_H_median  ',     &
					lap_FH_median
				write(u,'(A,ES23.15E3)')  'lower_profile_F_H_spread  ',     &
					lap_FH_spread
				write(u,'(A,ES23.15E3)')  'lower_profile_F_He_median ',     &
					lap_FHe_median
				write(u,'(A,ES23.15E3)')  'lower_profile_F_He_spread ',     &
					lap_FHe_spread
				write(u,'(A,ES23.15E3)')  'lower_profile_Mdot_median ',     &
					lap_Mdot_median
				write(u,'(A,ES23.15E3)')  'lower_profile_Mdot_spread ',     &
					lap_Mdot_spread
			endif
			! (b) the steady-flux window r >= r_esc.  At a steady state the
			! elemental flux does not depend on radius, so this IS the flux
			! through the matching level, measured where the solution is flat.
			write(u,'(A,ES23.15E3)')  'steady_flux_window_r_lo_Rp ',        &
				lap_steady_r_lo_Rp
			write(u,'(A,I0)')         'steady_flux_window_nface  ',         &
				lap_steady_nface
			write(u,'(A,ES23.15E3)')  'steady_F_H_median         ',         &
				lap_steady_FH_median
			write(u,'(A,ES23.15E3)')  'steady_F_H_spread         ',         &
				lap_steady_FH_spread
			write(u,'(A,ES23.15E3)')  'steady_F_He_median        ',         &
				lap_steady_FHe_median
			write(u,'(A,ES23.15E3)')  'steady_F_He_spread        ',         &
				lap_steady_FHe_spread
			write(u,'(A,ES23.15E3)')  'steady_Mdot_median        ',         &
				lap_steady_Mdot_median
			write(u,'(A,ES23.15E3)')  'steady_Mdot_spread        ',         &
				lap_steady_Mdot_spread
		endif
	endif
	close(u)
	end subroutine write_resolved_config

	! End of module
	end module setup_report

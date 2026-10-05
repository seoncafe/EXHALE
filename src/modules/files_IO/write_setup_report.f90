	module setup_report
	! Contains subroutine to write the setup of the current 
	!	simulation to file
	
	use global_parameters
	use base_boundary, only: base_reservoir_prescription_version,       &
	                         base_ghost_counts_measured,                &
	                         base_ghost_particle_count,                 &
	                         base_ghost_electron_count,                 &
	                         r_base_level, base_reservoir_p,            &
                            base_reservoir_T, base_face_mach_blend
	use charge_exchange, only: he_h_charge_exchange, cx_o2p_h_scale,      &
	                           cx_n2p_h_scale
	! Which spectrum built the grid; the single record of the choice.
	use J_incident,      only: spectrum_is_planck, loaded_table_floor_eV
	use sed_reader,      only: photon_grid_floor_eV
	use caloric_eos,     only: caloric_eos_state_line
	use Numerical_Fluxes, only: low_mach_velocity_jump
	use binary_element_diffusion, only: interdiffusion_enthalpy_scale
	use viscous_conduction, only: conduction_scale, conduction_active,  &
		lower_atmosphere_heat_measured, lower_atmosphere_heat_flux_cgs,   &
		lower_atmosphere_bath_T_K, lower_atmosphere_cell1_T_K
	use species_table,   only: melem_name
	use IC_load,         only: melem_from_abundance,                    &
	                     ic_coupling_present, ic_sec_ion_active,        &
	                     ic_provenance_unknown,                         &
	                     ic_restart_schema_present, n_opt, opt_name,    &
	                     restart_option_change_named,                   &
	                     restart_option_change_given,                   &
	                     ic_option_change_applied, ic_option_change_inert, &
	                     ic_legacy_factor_migrated
	use Read_input,      only: base_level_source, carrier_newton_on_stall,  &
	                           base_density_key_overridden_n0
	use grid_construction, only: domain_outer_radius, mixed_stretch_ratio, &
	                     outer_shells_width_ratio
	use composition,     only: h2_mixing_ratio_base, h2_mixing_ratio_ceiling
	use diffusive_photochemistry, only: carrier_co_domain_record,        &
	                     carrier_co_domain_f_dom
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

	! Number of metal elements with a non-zero abundance.
	integer :: imet, n_met_active
	! The flux the photon grid actually carries, and the nominal
	! (10^LX + 10^LEUV)/(4 pi a^2) at the same dayside dilution.
	real*8  :: f_grid, f_nominal, e_floor
	! The option tokens this run was allowed to change at the restart
	character(len=250) :: optlist
	integer :: iopt

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
	! WHAT THIS RUN IS DOING, and whether the input said so or the code took
	! its default. The default is init in every configuration: the ordinary
	! run of this code relaxes to a stationary solution, which is
	! initialization / continuation and claims no elapsed time. A physical
	! integration is asked for, carries
	! a clock that advances only on an accepted global step, and says so in
	! the header of the state it writes.
	if (run_mode .eq. run_mode_phys) then
		write(outfile,*) '- Run mode: phys (given) -- physical integration:'// &
		                 ' one global dt per step, a handoff test on the'//    &
		                 ' initial state, and an elapsed time'
	else if (run_mode_given) then
		write(outfile,*) '- Run mode: init (given) -- initialization /'//      &
		                 ' continuation: no physical elapsed time is claimed'//&
		                 ' for the state this run writes'
	else
		write(outfile,*) '- Run mode: init (default) -- initialization /'//    &
		                 ' continuation: no physical elapsed time is claimed'//&
		                 ' for the state this run writes'
	endif
	write(outfile,10) &
      ' - Upper boundary of the domain: ', domain_outer_radius(), ' [R_p]'
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
      ' erg cm^-2 s^-1 (912-1201 A, at the planet)'
			if (lw_from_spectrum) write(outfile,*) &
      '  (integrated from the numerical spectrum, not stated by the run)'
		else if (lw_flux_stated) then
			write(outfile,*) &
      '- H2 Lyman-Werner photodissociation OFF: "Stellar LW flux" stated'//&
      ' as zero (a comparison with a model that has no H2 photodissociation).'
		else
			write(outfile,*) &
      '- WARNING H2 Lyman-Werner photodissociation OFF: no "Stellar LW'//  &
      ' flux" key and no spectrum to integrate one from.'
			write(outfile,*) &
      '  The band exists whenever the star does, so this is a MISSING'//   &
      ' CHANNEL, not a modelling choice: H2 keeps a sink it should not.'
		endif
	endif
	! The final-state channels of the H2 photoabsorption. Which of them are
	! resolved changes the fragments, never the opacity: all four are shares
	! of one cross section (h2_photo_channels.f90).
	if (thereis_mol) then
		if (h2_double_ionization .ne. 'off') then
			write(outfile,*) &
      '- H2 double ionization RESOLVED (model '//                          &
      trim(h2_double_ionization)//'): H2 + hv -> H+ + H+ + 2e- above'//    &
      ' 51.4 eV, two protons per event'
		else
			write(outfile,*) &
      '- WARNING H2 double ionization NOT resolved: the protons'//        &
      ' it would release stay inside the single dissociative channel,'
			write(outfile,*) &
      '  which counts one proton per event instead of two above 51.4 eV'//&
      ' and adds an H atom the event does not make. This is the'
			write(outfile,*) &
      '  pre-E1 behavior and is kept only for reproducing an older run.'
		endif
		if (h2_neutral_dissociation) then
			write(outfile,*) &
      '- H2 neutral dissociation RESOLVED: H2 + hv -> H + H over the'//    &
      ' 33-41 eV window, no ion and no photoelectron'
		else
			write(outfile,*) &
      '- WARNING H2 neutral dissociation NOT resolved: the'//              &
      ' 33-41 eV photoionization yield is taken to be unity, as before,'
			write(outfile,*) &
      '  so up to 7.4 percent of the absorptions at 37.5 eV are counted'//&
      ' as ionizations rather than as two neutral H atoms. This is'
			write(outfile,*) &
      '  the pre-E1 behavior and is kept only for reproducing an'//        &
      ' older run.'
		endif
	endif
	! Whether the carriers are transported is a property of the molecular
	! network, not of the oxygen cycle: H2 is carrier 1 in both.
	if (thereis_mol) then
		if (carrier_transport) then
			if (thereis_oxychem) then
				write(outfile,*) &
      '- Molecular carriers (H2, OH, H2O, CO) are TRANSPORTED: implicit'// &
      ' diffusion-advection solved with their chemistry'
			else
				write(outfile,*) &
      '- The molecular carrier H2 is TRANSPORTED: implicit diffusion'//    &
      '-advection solved with the H2 balance of the same network'
			endif
			if (carrier_in_newton) then
				write(outfile,*) &
      '  and n(H2) is a FOURTH NEWTON UNKNOWN per cell: the wind and'//    &
      ' the carriers are solved together, not alternated (sec. 139)'
			endif
			if (carrier_newton_on_stall) then
				write(outfile,*) &
      '  and "Coupled carrier solve: On stall": the alternation runs'//    &
      ' first, and n(H2) becomes a FOURTH NEWTON UNKNOWN per cell from'
				write(outfile,*) &
      '  the pass at which the joint distance of the state has not'//      &
      ' fallen in three consecutive passes whose carrier relaxation'
				write(outfile,*) &
      '  ended on the composition movement bound. The handover is'//       &
      ' printed when it fires; the acceptance is the certification'
				write(outfile,*) &
      '  of the refreshed state either way.'
			endif
			! Stated only when it is not the default, so the report of a
			! run without the key is the one it always was.
			if (composition_update_holds_pressure) then
				write(outfile,*) &
      '  and "Composition update holds: pressure": the composition'//      &
      ' update of the stationary alternation keeps rho, v and p of'
				write(outfile,*) &
      '  every cell, T follows from p at the new particle count and E'//    &
      ' is rebuilt from the caloric EOS of the new composition'
				write(outfile,*) &
      '  (the default holds the conserved E and moves p). The energy'//     &
      ' the rebuild adds is printed at every pass.'
			endif
			if (maxval(kzz_cell) .le. 0.0d0) write(outfile,*) &
      '  WARNING K_zz = 0 everywhere, so the transport is pure molecular'//&
      ' diffusion -- the wrong limit for a lower atmosphere,'
			if (maxval(kzz_cell) .le. 0.0d0) write(outfile,*) &
      '  where eddy mixing is what holds the composition well mixed'//     &
      ' below the homopause. Set "He_Kzz:" or a profile.'
		else
			write(outfile,*) &
      '- Molecular carriers are a LOCAL steady state ("Molecular carrier'//&
      ' transport: False"): the chemistry alone, not a model of a base'
		endif
	endif
	if (ionization_transport) then
		write(outfile,*) &
      '- The IONIZATION STATE of hydrogen and helium is TRANSPORTED'//     &
      ' too: the carried set is x(H II) = n(H II)/n(H nuclei),'//         &
      ' x(He II) = n(He II)/n(He nuclei) and x(He III) ='//               &
      ' n(He III)/n(He nuclei), each a fraction per element NUCLEUS,'//   &
      ' so the ionization partition is what the flow accumulated and'//   &
      ' not the local equilibrium'
		write(outfile,*) &
      '  a stage flux is x N_el - n_el K_zz dx/dr with N_el that'//        &
      ' element''s nucleus flux from the element operator, so the'//      &
      ' stage fluxes of an element sum to its element flux face by face'
		write(outfile,*) &
      '  the neutral stage of each element closes its simplex and is'//    &
      ' not a row: H I with the hydrogen the molecules hold, He I with'// &
      ' the helium of HeH+'
		write(outfile,*) &
      '  (for a wind in which P r/|v| < 1 the local root over-ionizes;'//  &
      ' md/k22_electron_density_excess.md sec. 7)'
		if (thereis_mol) then
			write(outfile,*) &
      '  the source of each stage row is rows (1), (2) and (3) of the'//   &
      ' molecular H/He network, the same rows the equilibrium sweep of'// &
      ' this gas solves'
		else
			write(outfile,*) &
      '  the source of each stage row is the H/He ionization balance of'//&
      ' an atomic gas -- photoionization with its secondaries,'//         &
      ' collisional ionization, radiative and dielectronic'//             &
      ' recombination, He <-> H charge exchange and the He 2^3S'//        &
      ' channels where the level is tracked -- the same rows the'//       &
      ' equilibrium sweep of this gas solves'
		endif
		if (he23s_transport) then
			write(outfile,*) &
      '- The He 2^3S POPULATION is TRANSPORTED as well: x3 ='//            &
      ' n(2^3S)/n(He nuclei) rides on the helium nucleus flux with its'// &
      ' own eddy term, its source is the level balance the sweep solves'//&
      ' (tr_triplet_row), and the ground singlet closes the helium'//     &
      ' simplex x(He II) + x(He III) + x3'
		endif
	endif
	if (thereis_oxychem) then
		write(outfile,*) &
      '- Including oxygen chemistry (OH, H2O, CO) in the same system'
		write(outfile,'(A,4ES10.3)') &
      ' - FUV band fluxes at the planet [erg cm^-2 s^-1]'//               &
      ' LW/B2(Lya)/B3/B4: ', F_LW_star, F_Lya_star,                       &
      F_FUV_B3, F_FUV_B4
		if (fuv_b3_from_spectrum .or. fuv_b4_from_spectrum) then
			write(outfile,'(A,L1,A,L1,A)') &
      '  (B3 integrated from the numerical spectrum: ', fuv_b3_from_spectrum, &
      '; B4: ', fuv_b4_from_spectrum, '; a stated key always wins)'
		endif
		if (fuv_b3_flux_stated .and. F_FUV_B3 .le. 0.0d0) write(outfile,*) &
      '  - band B3 OFF: "Stellar FUV B3 flux" stated as zero'
		if (fuv_b4_flux_stated .and. F_FUV_B4 .le. 0.0d0) write(outfile,*) &
      '  - band B4 OFF: "Stellar FUV B4 flux" stated as zero'
		write(outfile,*) &
      '- The 912-1201 A band is the Lyman-Werner interval and carries'//  &
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
	! Photoelectron secondary ionization. The STAGED default is reported
	! explicitly because it is not a synonym for "on": the coupling is
	! applied only after the wind first converges without it, so a run that
	! is stopped at a fixed step count before converging never applies it at
	! all. That is what silently excluded every molecular regression case
	! from this physics until 2026-09-02.
	if (.not. use_sec_ion) then
		write(outfile,*) &
         ' - Photoelectron secondary ionization: off'
	else if (sec_ion_immediate) then
		write(outfile,*) &
         ' - Photoelectron secondary ionization: on from step 0 (Immediate)'
	else
		write(outfile,*) &
         ' - Photoelectron secondary ionization: STAGED -- applied only'// &
         ' after the first du convergence, so a run stopped before that'// &
         ' never applies it'
	endif
	! The photoelectrons below the 30 eV threshold of that partition. Only
	! the non-default choice is reported, so a run without the key writes
	! exactly what it wrote before.
	if (use_low_energy_partition) then
		if (use_sec_ion) then
			write(outfile,*) &
         ' - Photoelectrons of 10.2-30 eV: Furlanetto & Stoever (2010)'// &
         ' partition in cells without H2, applied with the secondary'// &
         ' ionization above'
		else
			write(outfile,*) &
         ' - Photoelectrons of 10.2-30 eV: Furlanetto & Stoever (2010)'// &
         ' partition selected, but secondary ionization is off, so it'// &
         ' is never applied'
		endif
	endif
	! Whether the staging was overruled by the restart file. A file written
	! by a run that had already armed the coupling restarts with it armed;
	! the line says so, and says when the file could not be asked.
	if (do_load_IC .and. use_sec_ion .and. .not. sec_ion_immediate) then
		if (ic_coupling_present) then
			if (ic_sec_ion_active) then
				write(outfile,*) &
            ' - Restart coupling header: the file states secondary'// &
            ' ionization ON; the staging is overruled and it is applied'// &
            ' from step 0'
			else
				write(outfile,*) &
            ' - Restart coupling header: the file states secondary'// &
            ' ionization OFF; the run starts staged'
			endif
		else
			write(outfile,*) &
         ' - Restart coupling header: absent (the IC predates it);'// &
         ' the run starts staged, as it did before'
		endif
	endif
	! WHAT IS KNOWN ABOUT THE CONFIGURATION THE LOADED STATE WAS PRODUCED
	! UNDER, and what this run was allowed to change about it.
	! A state whose configuration cannot be compared with this run's is not
	! a state this run can be held to, and a state reloaded under changed
	! physics options is a starting point and not a solution: both are
	! statements about the run, so both belong in the report of the resolved
	! configuration and not only in the log.
	if (do_load_IC) then
		if (ic_provenance_unknown .and. .not. ic_restart_schema_present) then
			write(outfile,*) &
         ' - Restart provenance: UNKNOWN. The state files carry no'// &
         ' configuration block, so the reservoir, the grid, the'// &
         ' constants and the physics options they were produced'// &
         ' under were not compared with this run''s; every file this'// &
         ' run writes says so'
		else if (ic_provenance_unknown) then
			write(outfile,*) &
         ' - Restart provenance: UNKNOWN by descent. The state files'// &
         ' state their own configuration and it agrees with this'// &
         ' run, but the state they carry came from a file that'// &
         ' carried no configuration block, so what it was reached'// &
         ' under is not known; the mark stays with every file this'// &
         ' run writes'
		else
			write(outfile,*) &
         ' - Restart provenance: the state files state the'// &
         ' configuration they were produced under, and it agrees'// &
         ' with this run'
		endif
		if (restart_option_change_given) then
			optlist = ' '
			do iopt = 1, n_opt
				if (restart_option_change_named(iopt)) &
					optlist = trim(optlist)//' '//trim(opt_name(iopt))
			enddo
			write(outfile,*) &
         ' - Restart option change: allowed for'//trim(optlist)// &
         '; every other option difference refuses the load'
			if (ic_option_change_applied) then
				write(outfile,*) &
         '   the load DID change a named option (the lines'// &
         ' "# option_change" of the output state say which): the'// &
         ' loaded state is a solution of the former option set'
			endif
			if (ic_option_change_inert) then
				write(outfile,*) &
         '   at least one named option does not differ between the'// &
         ' state and this run'
			endif
		endif
		! A legacy F7.5 factor token read as a change of representation
		! (compare_options_field): the state is a seed, whatever the key.
		if (ic_legacy_factor_migrated) then
			write(outfile,*) &
         ' - Restart legacy factor token: the state files carry a'// &
         ' continuation factor in the F7.5 form (the "#'// &
         ' legacy_factor_token" line of the output state names it); the'// &
         ' requested factor lies inside its rounding interval, so the'// &
         ' state is loaded as a seed and its inherited claim is dropped'
		endif
	endif
	! Atomic H/He rate set. Reported only when the published set is in
	! force: it is the statement that this run is not on EXHALE's own rates.
	if (photoheat_photon_fraction .gt. 0.0d0) then
		write(outfile,'(A,F6.3,A)') '  - Photoelectron heating: a fraction ',      &
			photoheat_photon_fraction, ' of the photon energy h nu per '//      &
			'ionization (comparison option; the physical default is h nu - I)'
	endif
	if (atomic_rate_set_k22) then
		write(outfile,*) &
         ' - Atomic H/He rate set: Koskinen et al. (2022) Table 1 R1-R4.'// &
         ' H+ and He+ radiative recombination are the Storey & Hummer'
		write(outfile,*) &
         '   (1995) power laws 4.0e-12 and 4.6e-12 (300/T)^0.64 in place'// &
         ' of the Badnell/Milne case B; H and He collisional ionization'
		write(outfile,*) &
         '   are the Voronov (1997) fit, which is the default set as'//     &
         ' well. The recombination cooling follows the coefficients.'
	endif
	if (do_read_sed) &
		write(outfile,*) & 
         ' - Spectrum read from external file: ', sed_file
	if (is_PL_sed) &
      	write(outfile,12)  & 
            ' - Using power-law spectrum ', 'with index ', PLind
	if (is_monochr) &
		write(outfile,13)  &
         ' - Using monochromatic radiation with energy ', e_low
	! ONE SPECTRUM TYPE BUILDS EVERY BAND. Stated for every run: which type
	! filled the photon grid, what the band BELOW 13.6 eV was built from --
	! the band where the He 2^3S metastable (4.80 eV) and the low-IP metals
	! absorb, and the one band a reader is most likely to assume came from
	! somewhere else -- and what the grid integrates to against the nominal
	! XUV flux.
	write(outfile,*) '- Spectrum type: '//trim(sp_type)//                  &
      ', and it builds every band of the photon grid'
	! The floor is the lower EDGE of the first bin (e_v(1) is its center):
	! the lowest active threshold for the analytic types, the table's
	! lowest read energy for a loaded SED.
	e_floor = photon_grid_floor_eV()
	if (do_read_sed) e_floor = loaded_table_floor_eV()
	if (.not. is_monochr .and. allocated(e_v)) then
		if (NlTR .le. 0) then
			write(outfile,'(A,F7.2,A)') '   Below 13.6 eV: no grid'//     &
            ' point. The grid starts at ', e_floor, ' eV, no absorber'//   &
            ' below the H I threshold being active.'
		else if (is_PL_sed) then
			write(outfile,'(A,F6.2,A)') '   Below 13.6 eV (down to ',     &
            e_floor, ' eV): the SAME power law, extrapolated below the'
			write(outfile,*) '     [e_low, e_mid] band it is normalized'//&
            ' on. A photospheric near-ultraviolet field is a different'
			write(outfile,*) '     spectrum type ("Planck", or a "Load"'//&
            'ed SED reaching that wavelength), not a correction to this one.'
		else if (spectrum_is_planck()) then
			write(outfile,'(A,F6.2,A)') '   Below 13.6 eV (down to ',     &
            e_floor, ' eV): the same photospheric blackbody as every'//    &
            ' other band.'
		else if (do_read_sed) then
			write(outfile,'(A,F6.2,A)') '   Below 13.6 eV (down to ',     &
            e_floor, ' eV): the same SED file as every other band.'
		endif
	endif
	if (allocated(F_XUV) .and. allocated(de_v)) then
		f_grid    = sum(F_XUV*de_v)
		f_nominal = dayside_dilution()*J_XUV
		write(outfile,'(A,ES11.4,A)') '   Integrated grid flux'//         &
         ' sum(F dE) = ', f_grid, ' erg cm^-2 s^-1 (after the dilution)'
		write(outfile,'(A,ES11.4,A,ES10.3)') '     against the nominal'//   &
         ' (10^LX + 10^LEUV)/(4 pi a^2) = ', f_nominal,                   &
         ' at the same dilution: ratio ', f_grid/max(f_nominal,1.0d-300)
	endif
	if (appx_mth.eq.'alpha') then 
		write(outfile,14) ' - 2D approximation used: alpha =, with alpha = ',a_tau
	else
		write(outfile,*) '- 2D approximation used: ', appx_mth
	endif
	! The one number that convention turns into, stated once: every stellar
	! beam of the run -- XUV grid, Ly-alpha, and the four FUV bands -- is
	! multiplied by it (global_parameters dayside_dilution).
	write(outfile,'(A,F6.3)') &
      ' - Dayside dilution applied to every stellar beam: ',              &
      dayside_dilution()
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
		! The width to all its digits, and where it came from. Two inputs
		! that read alike build two grids whose states do not load into
		! each other (load_IC compares cell centers to 1e-10), so the
		! report states the value that reproduces this grid, not a
		! rounded spelling of it.
		if (dr_base_from_key) then
			write(outfile,'(A,A,A)') '   dr_base = ',                      &
				trim(round_trip_decimal(dr_base)), ' (from the key)'
		else
			write(outfile,'(A,A,A)') '   dr_base = ',                      &
				trim(round_trip_decimal(dr_base)), ' (the default)'
		endif
	else
		write(outfile,*) '- Base grid: "Base grid [dr,cells]" is ignored'//&
			' by grid type '//trim(grid_type)//' (Mixed only)'
	endif
	! The shells beyond the constructed grid, where the two faces are and
	! the three width ratios that say whether the stretch continues across
	! the old outer face: the ratio of the shells (set by their cell
	! count), the ratio of the constructed grid's stretched region, and the
	! jump dr(nc+1)/dr(nc) at the old outer face itself.
	if (n_outer_shells .gt. 0) then
		write(outfile,'(A,I0,A,I0,A)') ' - Outer shells: ',               &
			n_outer_shells, ' cells beyond the constructed grid of ',      &
			N - n_outer_shells, ' cells'
		write(outfile,'(A,ES23.16,A,ES23.16,A)')                           &
			'   outer face of the constructed grid ',                      &
			r_edg(N - n_outer_shells), ' R_p, of the shells ', r_edg(N),   &
			' R_p'
		write(outfile,'(A,F12.8,A,F12.8,A,F12.8)')                         &
			'   width ratio of the shells ', outer_shells_width_ratio,      &
			'; of the Mixed stretched region ', mixed_stretch_ratio,       &
			'; across the old outer face ',                                &
			dr_j(N - n_outer_shells + 1)/dr_j(N - n_outer_shells)
	endif
	! Resolution of the base density scale height H = kT_eq/(mu g) in cells:
	! b0 = R_p/H(T_eq) is the Jeans parameter, dr_j(1) the first cell size
	! after the Mixed-grid smoothing. Below ~5 cells per H the scheme carries
	! an undamped stationary 2*dr entropy mode.
	write(outfile,17) '- Base scale-height resolution: H(T_eq)/dr = ',      &
		1.0d0/(b0*dr_j(1)), ' cells'
	if (1.0d0/(b0*dr_j(1)) .lt. 10.0d0)                                     &
		write(outfile,*) '  WARNING: base scale height spans < 10 cells;'//&
			' expect a stationary cell-to-cell entropy mode at the base.'
	! Which input fixed the lower boundary, and where it put it. n0 and the
	! base pressure are one statement -- p = n0 k T0 ntot_bc at the base
	! composition -- so the report gives the density the run marches with,
	! the pressure that density is, and the input those came from, in one
	! line. The pressure is recomputed from n0 rather than echoed from the
	! input, so a level stated twice cannot be reported as two levels.
	write(outfile,71) '- Base level: n0 = ', n0, ' cm^-3 -> p = ',          &
		n0*kb_erg*T0*ntot_bc*1.0d-6, ' bar (level from '//                  &
		trim(base_level_source)//')'
71	format(1X,A,ES11.4,A,ES11.4,A)
	if (base_density_key_overridden_n0 .gt. 0.0d0)                          &
		write(outfile,71) '  WARNING: the density key asked for n0 = ',   &
			base_density_key_overridden_n0, ' cm^-3 -> p = ',               &
			base_density_key_overridden_n0*kb_erg*T0*ntot_bc*1.0d-6,        &
			' bar; "Base BC: pressure" overrides it (no handoff)'
	write(outfile,*) '- Caloric EOS (energy <-> pressure): ',              &
		trim(caloric_eos_state_line())
	if (thereis_mol) then
		if (mol_reaction_heat) then
			write(outfile,*) &
      '- Molecular reaction heat ON: the collisional H2/He network'//     &
      ' deposits its chemical energy (I(H2) - D0(H2) = 10.95 eV per'//    &
      ' photon-driven cycle, via H3+/H2+ dissociative recombination)'
		else
			write(outfile,*) &
      '- WARNING Molecular reaction heat OFF: the ionization energy the'//&
      ' H2 photoabsorption spends is returned to the gas by dissociative'
			write(outfile,*) &
      '  recombination, not to a photon, and with this off the code'//    &
      ' keeps none of it. Measured at 81 percent of the total heating'
			write(outfile,*) &
      '  rate at 1.02 r_base on the He/H = 0.0793 rung, so this is a'//   &
      ' MISSING CHANNEL, not a modelling choice.'
		endif
	endif
	write(outfile,*) '- Numerical flux: ', flux
	write(outfile,*) '- Reconstruction method: ', rec_method
	! The pressure/gravity pair: whether the reconstruction, the Riemann
	! jumps and the pressure force carry the state or its departure from the
	! cell's own hydrostatic equilibrium.
	if (well_balanced) then
		write(outfile,*) '- Well balanced: on (the reconstruction and the'// &
			' Riemann jumps carry the departure from the local'
		write(outfile,*) '    hydrostatic equilibrium of constant density,'//&
			' and the momentum source is that equilibrium'
		write(outfile,*) '    face-pressure difference, so a discrete'//     &
			' equilibrium is preserved to rounding.'
		if (trim(flux) .eq. 'LLF') write(outfile,*) '    WARNING: the LLF'// &
			' flux resolves no stationary contact, so the'//                 &
			' well-balanced property does not hold with it.'
	else
		write(outfile,*) '- Well balanced: off (opt-in key "Well balanced")'
	endif
	! Heat conduction and the conductivity it uses.
	if (cond_on) then
		write(outfile,*) '- Heat conduction: on, kappa = (n_e 1.2e-6 T^2.5'// &
			' + n_HI 379 T^0.69 + n_HeI 299 T^0.69 + n_H2 k_H2(T))/n'//      &
			' erg/cm/s/K'
		write(outfile,*) '    (Banks & Kockarts 1973 as Salz et al. 2015'//  &
			' Eqs. 24-26 and Sutton et al. 2015 Eq. 5 state them; k_H2 the'
		write(outfile,*) '    Eucken-form fit to Incropera et al. 2007 Table A.4,'// &
			' valid 200-2000 K)'
		if (conduction_scale() .ne. 1.0d0)                                 &
			write(outfile,'(A,ES13.6,A)') '  - CONTINUATION FACTOR on the'// &
			' conductivity: ', conduction_scale(),                         &
			' (EXHALE_CONDUCTION_SCALE; a step, not the model)'
	endif
	! The energy carried by the element fluxes. Stated only where the
	! elements move, so the report of a run without He_diffusion is the one
	! it always was.
	if (he_diffusion) then
		if (interdiffusion_enthalpy_flux) then
			write(outfile,*) '- Interdiffusion enthalpy flux: on (the'//      &
				' energy equation carries q_d = sum_s h_s J_s of the'
			write(outfile,*) '    element fluxes, Cook 2009 eqs. 11-13,'//   &
				' in the marching update and the stationary energy row)'
			if (interdiffusion_enthalpy_scale() .ne. 1.0d0)                &
				write(outfile,'(A,ES13.6,A)') '  - CONTINUATION FACTOR on'//  &
				' that term: ', interdiffusion_enthalpy_scale(),               &
				' (EXHALE_INTERDIFF_ENTH_SCALE; a step, not the model)'
		else
			write(outfile,*) '- Interdiffusion enthalpy flux: OFF by'//       &
				' request (the elements move but the energy equation'
			write(outfile,*) '    omits the enthalpy they carry, as'//       &
				' Koskinen et al. 2013, 2022 and Yelle 2004 do)'
		endif
	endif
	! Rieper (2011) low-Mach correction of the Roe dissipation: the normal
	! velocity jump entering the two acoustic expansion coefficients is
	! scaled by min(|U_Roe|/a_Roe, 1).  It is defined on the ROE branch only.
	if (low_mach_velocity_jump) then
		write(outfile,*) '- Low Mach velocity jump: on (the normal'//     &
			' velocity jump of the Roe dissipation is scaled by'
		write(outfile,*) '    min(|U_Roe|/a_Roe, 1), Rieper 2011 J.'//    &
			' Comput. Phys. 230, 5263, eq. 3.15-3.16; an accuracy'
		write(outfile,*) '    correction of the low-Mach regime, not a'//&
			' stiffness one, and exactly 1 at a supersonic face.'
		if (trim(flux) .ne. 'ROE') write(outfile,*) '    WARNING: the'//  &
			' correction is derived for the Roe flux and this run'//        &
			' selects '//trim(flux)//', so the key has NO effect.'
	else
		write(outfile,*) '- Low Mach velocity jump: off (opt-in key'//    &
			' "Low Mach velocity jump", ROE branch only)'
	endif
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
	! Lower boundary: one characteristic condition at the face, so the report
	! states WHERE it is and WHAT the reservoir says, not which ghost closure
	! is selected -- there is no longer a choice of closure to report.
	write(outfile,*) '- Base boundary: characteristic condition at the'//   &
		' face r_edg(0), which is the base level'
	write(outfile,23) '    reservoir (p, s) at r = ', r_base_level,         &
		': p = ', base_reservoir_p, ' p0, T = ', base_reservoir_T, ' T0'
	write(outfile,21) '    supersonic-outflow regularization width, face Mach = ',      &
		base_face_mach_blend
	if (shapiro_eps .gt. 0.0d0) then
		write(outfile,22) '- Shapiro filter: on, eps = ', shapiro_eps,       &
			', applied every ', shapiro_every, ' steps'
	else
		write(outfile,*) '- Shapiro filter: off (opt-in key "Shapiro filter")'
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
	write(outfile,19) '- Max marching steps: ', marching_step_max
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
21 format(A,ES10.3)
22 format(A,ES10.3,A,I0,A)
23 format(A,F8.4,A,ES11.4,A,F8.4,A)

	write(*,*) '(write_setup_report.f90) Done.'

	end subroutine write_setup_report

	! ------------------------------ !

	subroutine write_parse_dump
	! Dump every variable input_read derives from input.inp (plus any base.inp
	! override) to parse_dump.txt, one "name = value" line per variable, in the
	! md/input_schema.md key order. Gated by EXHALE_PARSE_DUMP=1 in
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
	write(u,'(A)') '# order follows md/input_schema.md; cgs where input_read converts'

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
	! The photospheric-blackbody type has no logical of its own: sp_type,
	! dumped above, is the record of the choice, and T_star_eff and R_star,
	! dumped below, are the field it builds.
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
	! Only when "Outer shells" is given, so that every parse-corpus case
	! without the key dumps exactly what it dumped before.
	if (n_outer_shells .gt. 0) then
		call put_i('n_outer_shells', n_outer_shells)
		call put_r('r_outer_shells_face', r_outer_shells_face)
	endif
	call put_s('flux', flux)
	call put_s('rec_method', rec_method)
	call put_l('use_weno3', use_weno3)
	call put_l('use_plm', use_plm)
	call put_l('recon_two_stage', recon_two_stage)
	call put_l('well_balanced', well_balanced)
	! Only when the PLM -> WENO3 continuation is actually asked for.
	! An absent key leaves the shipped one-step hand-off, which has no
	! continuation settings to resolve, and every parse-corpus case
	! that does not use it dumps exactly what it dumped before.
	if (recon_lambda_step0 .gt. 0.0d0) then
		call put_r('recon_lambda_step0', recon_lambda_step0)
		call put_r('recon_lambda_dtu_tol', recon_lambda_dtu_tol)
		call put_l('recon_lambda_adaptive', recon_lambda_adaptive)
	endif
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
	call put_l('atomic_rate_set_k22', atomic_rate_set_k22)
	call put_l('use_sec_ion', use_sec_ion)
	! Only when selected, so that every parse-corpus case without the key
	! dumps exactly what it dumped before.
	if (use_low_energy_partition)                                        &
		call put_l('use_low_energy_partition', use_low_energy_partition)
	call put_l('use_he_rec_coupling', use_he_rec_coupling)
	call put_l('use_h_rec_escape', use_h_rec_escape)
	call put_l('he_h_charge_exchange', he_h_charge_exchange)
	call put_r('cx_o2p_h_scale', cx_o2p_h_scale)
	call put_r('cx_n2p_h_scale', cx_n2p_h_scale)
	call put_l('thereis_mol', thereis_mol)
	call put_s('h2_double_ionization', h2_double_ionization)
	call put_l('h2_neutral_dissociation', h2_neutral_dissociation)
	call put_r('F_LW_star', F_LW_star)
	call put_l('thereis_oxychem', thereis_oxychem)
	call put_l('carrier_transport', carrier_transport)
	call put_l('carrier_in_newton', carrier_in_newton)
	call put_l('carrier_newton_on_stall', carrier_newton_on_stall)
	call put_l('ionization_transport', ionization_transport)
	! Written only for the non-default value, so the dump of every run
	! without the key is the one it was.
	if (he23s_transport) call put_l('he23s_transport', he23s_transport)
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
	call put_l('interdiffusion_enthalpy_flux', interdiffusion_enthalpy_flux)
	call put_r('stall_tol', stall_tol)
	call put_i('N_stall', N_stall)
	call put_l('use_semi_implicit_energy', use_semi_implicit_energy)
	call put_l('use_local_dt', use_local_dt)
	call put_r('lev_th', lev_th)
	call put_l('use_newton_solver', use_newton_solver)
	call put_r('newton_du_switch', newton_du_switch)
	call put_r('shapiro_eps', shapiro_eps)
	call put_i('shapiro_every', shapiro_every)
	call put_r('lowmach_damp_eps', lowmach_damp_eps)
	call put_r('lowmach_damp_mach_th', lowmach_damp_mach_th)
	call put_i('count_max', marching_step_max)
	call put_r('coronal_cutoff_width', coronal_cutoff_width)
	call put_l('base_ir_field', base_ir_field)
	call put_l('mol_ir_bands', mol_ir_bands)
	call put_i('base_bc_mode', base_bc_mode)
	call put_r('base_p_ubar', base_p_ubar)
	call put_r('base_reservoir_p', base_reservoir_p)
	call put_r('base_reservoir_T', base_reservoir_T)
	call put_l('visc_on', visc_on)
	call put_l('cond_on', cond_on)
	call put_r('visc_mu0', visc_mu0)
	call put_r('visc_s', visc_s)
	call put_r('resid_th', resid_th)
	call put_r('flux_spread_th', flux_spread_th)
	call put_r('r_flux', r_flux)
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
	! WHICH COUNT THIS IS. ntot_bc and dp_bc below are the PRESCRIBED base
	! reservoir's counts, version 1 of that prescription: the base nuclei
	! count per unit n0 corrected for H2 binding, and the electron term of
	! the same level. They are resolved once at startup and no sweep writes
	! them back. What a composition sweep MEASURES in the lower ghost is a
	! different quantity under a similar name, and it is written beside them
	! in EXHALE_resolved.out (base_ghost_count_solved); this dump is a
	! parse record compared line by line, so the statement is made here and
	! not in it.
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
	! transit geometry as well.  Format: '# ' comments, then one
	! 'key  value' pair per line.
	integer :: u, ie
	logical :: base_present
	! Domain record of the one-sided CO destruction model over the run so
	! far (zero at the startup call, the whole run at the end-of-run call).
	integer :: co_dom_out, co_dom_hot, co_dom_hep
	real*8  :: co_dom_ratio, co_dom_r, co_dom_form

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
	! The discretization of the pressure/gravity pair, which decides what a
	! profile written by this run is a steady state OF.
	write(u,'(A,L1)')     'well_balanced             ', well_balanced
	! THE GRID THIS RUN IS SOLVED ON, stated so that it can be rebuilt from
	! this record alone. The width is written at round-trip precision
	! because a grid is not reproduced by a rounded width: two widths that
	! read alike to a human build two grids whose states do not load into
	! each other (load_IC compares cell centers to 1e-10 relative), and a
	! reproduction input must carry the digits, not the spelling of the key
	! that was typed. base_cell_width_source says whether the width came
	! from "Base grid [dr,cells]:" or from the code's default; value
	! equality cannot answer that, because a key may state the default's own
	! digits. base_grid_in_effect is F for the Uniform and Stretched grid
	! types, which ignore the width and the cell count.
	write(u,'(A,A)')      'grid_type                 ', trim(grid_type)
	! The "Grid cells" of the constructed grid; the shells, when there
	! are any, are the two lines after base_grid_in_effect, and the run
	! solves on grid_cells + outer_shells_cells cells.
	write(u,'(A,I0)')     'grid_cells                ', N - n_outer_shells
	write(u,'(A,A)')      'outer_radius_Rp           ',                  &
		trim(round_trip_decimal(r_max))
	write(u,'(A,A)')      'base_cell_width_Rp        ',                  &
		trim(round_trip_decimal(dr_base))
	write(u,'(A,I0)')     'base_uniform_cells        ', N_low_cells
	if (dr_base_from_key) then
		write(u,'(A)')     'base_cell_width_source    key'
	else
		write(u,'(A)')     'base_cell_width_source    default'
	endif
	write(u,'(A,L1)')     'base_grid_in_effect       ',                  &
		grid_type .eq. 'Mixed'
	if (n_outer_shells .gt. 0) then
		write(u,'(A,I0)')  'outer_shells_cells        ', n_outer_shells
		write(u,'(A,A)')   'outer_shells_face_Rp      ',                  &
			trim(round_trip_decimal(r_outer_shells_face))
	endif
	! Which set of atomic H/He rate coefficients the run used: 'K22' = the
	! four Koskinen et al. (2022) Table 1 entries R1-R4 selected by "Atomic
	! rate set:", 'default' = EXHALE's own (Badnell/Milne case B, or the
	! legacy ATES fits under Legacy_HHe_rates). The recombination cooling
	! follows the coefficients either way.
	if (atomic_rate_set_k22) then
		write(u,'(A)')     'atomic_rate_set           K22'
	else
		write(u,'(A)')     'atomic_rate_set           default'
	endif
		if (photoheat_photon_fraction .gt. 0.0d0) then
			write(u,'(A,F8.4)') 'photoelectron_heating     ', photoheat_photon_fraction
		else
			write(u,'(A)')     'photoelectron_heating     excess'
		endif
	! The partition of the 10.2-30 eV photoelectrons. Written only when it is
	! selected, so that a run without the key writes the record it wrote
	! before; absent means none (deposited whole as heat).
	if (use_low_energy_partition)                                        &
		write(u,'(A)')     'low_energy_electron_partition Furlanetto2010'
	! Oxygen chemistry (the A2 option). oxygen_chemistry says whether the
	! O I column means FREE ATOMIC oxygen (it does when this is T) and
	! whether the OH / H2O / CO columns of Ion_species.txt exist;
	! oxygen_base_partition is the provenance of the base H2/H partition,
	! which is the whole point of the option. The four band fluxes are the
	! photon input the oxygen photochemistry actually ran on; the first of
	! them, fuv_band_LW_flux, is the 912-1201 A interval shared with the H2
	! Lyman-Werner absorber.  There was a fifth, fuv_band_B1_flux, for the
	! 1110-1201 A band merged into it on 2026-09-06.
	write(u,'(A,L1)')     'mol_ir_bands              ', mol_ir_bands
	write(u,'(A,L1)')     'carrier_transport         ', carrier_transport
	write(u,'(A,L1)')     'carrier_in_newton         ', carrier_in_newton
	write(u,'(A,L1)')     'carrier_newton_on_stall   ', carrier_newton_on_stall
	! Written only for the non-default value; a file without the line is
	! a run whose composition update held the conserved energy.
	if (composition_update_holds_pressure)                              &
		write(u,'(A)')     'composition_update_holds  pressure'
	write(u,'(A,L1)')     'ionization_transport      ', ionization_transport
	! Written only for the non-default value; a file without the line is a
	! run whose He 2^3S population is the local root of each cell.
	if (he23s_transport)                                                 &
		write(u,'(A,L1)') 'he23s_transport           ', he23s_transport
	write(u,'(A,L1)')     'oxygen_chemistry          ', thereis_oxychem
	if (thereis_oxychem) then
		write(u,'(A)')     'oxygen_reaction_set       a2_v1'
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
	! THE PRESCRIBED RESERVOIR COUNT, and not what a sweep measured in the
	! ghost. ntot_bc is the base nuclei count per unit n0 corrected for H2
	! binding, resolved once at startup from the base mixing ratio; version 1
	! of the reservoir prescription (base_boundary). The ghost's own counts,
	! where a sweep has measured them, are the two lines below it: they are a
	! statement about the state and are never fed back into the reservoir.
	write(u,'(A)') '# ntot_bc_per_H is the PRESCRIBED base reservoir'//   &
	     ' count (reservoir prescription version 1);'
	write(u,'(A)') '#   base_ghost_count_solved is what the composition'//&
	     ' sweep measured in the lower ghost.'
	write(u,'(A,ES23.15E3)') 'ntot_bc_per_H             ', ntot_bc
	write(u,'(A,I0)')        'base_reservoir_version    ',               &
	     base_reservoir_prescription_version
	if (base_ghost_counts_measured) then
		write(u,'(A,ES23.15E3)') 'base_ghost_count_solved   ',            &
		     base_ghost_particle_count
		write(u,'(A,ES23.15E3)') 'base_ghost_electrons_solved',           &
		     base_ghost_electron_count
	endif
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
	! The elemental fluxes measured over the overlap window. Reported only
	! when the certification has measured a state and produced them: this
	! routine also runs before the wind, and an unmeasured flux must say so
	! rather than print a zero.
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
				write(u,'(A)') '# from which the face mass flux of the'//   &
					' state is constant, so no face of the'
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
	! THE HEAT THE WIND CONDUCTS INTO THE LOWER ATMOSPHERE through the base
	! face, on the state the certification measured
	! (heat_conducted_into_lower_atmosphere, viscous_conduction.f90).
	! Written whenever conduction is on; before the wind exists it says
	! unmeasured rather than print a zero.  The energy coupling of the
	! elemental-flux closure (element_flux_closure.py) reads the flux.
	if (conduction_active()) then
		if (lower_atmosphere_heat_measured) then
			write(u,'(A)') 'base_conduction_state     measured'
			write(u,'(A)') '# base_conductive_flux_cgs [erg cm^-2 s^-1]:'//  &
				' q = -kappa_face (T_1 - T_bath)/(r_1 - r_0),'
			write(u,'(A)') '# positive = heat leaving the wind into the'//  &
				' lower atmosphere; T_bath at the ghost centre r(0)'
			write(u,'(A,ES23.15E3)') 'base_conductive_flux_cgs  ',         &
				lower_atmosphere_heat_flux_cgs
			write(u,'(A,ES23.15E3)') 'base_conduction_T_bath_K  ',         &
				lower_atmosphere_bath_T_K
			write(u,'(A,ES23.15E3)') 'base_conduction_T_cell1_K ',         &
				lower_atmosphere_cell1_T_K
		else
			write(u,'(A)') 'base_conduction_state     unmeasured'
		endif
	endif
	! THE DOMAIN OF THE ONE-SIDED CO DESTRUCTION MODEL, CUMULATIVE OVER
	! THE RUN.  The CO row destroys CO and never forms it, which is
	! legitimate in a cell where tau_dest << tau_res << tau_form.  Both
	! inequalities are evaluated on the state the accepted carrier steps
	! hand on, and the counts are cell visits summed over the run's carrier
	! intervals.  The record is informational: the rates are on everywhere,
	! and where the ordering fails it is the omitted formation that fails
	! with it, so the transported value stands.  Written even when nothing
	! was out of domain, because a zero in a cumulative record IS the
	! statement that the model stayed inside its domain.
	if (carrier_transport .and. thereis_oxychem) then
		call carrier_co_domain_record(co_dom_out, co_dom_hot,           &
			co_dom_hep, co_dom_ratio, co_dom_r, co_dom_form)
		write(u,'(A,ES23.15E3)') 'co_domain_f_dom           ',          &
			carrier_co_domain_f_dom()
		write(u,'(A,I0)')        'co_domain_cells_out       ',          &
			co_dom_out
		write(u,'(A,ES23.15E3)') 'co_domain_worst_ratio     ',          &
			co_dom_ratio
		write(u,'(A,ES23.15E3)') 'co_domain_worst_r         ',          &
			co_dom_r
		write(u,'(A,ES23.15E3)') 'co_domain_form_ratio      ',          &
			co_dom_form
		write(u,'(A,I0)')        'co_domain_cells_hot       ',          &
			co_dom_hot
		write(u,'(A,I0)')        'co_domain_cells_HeII      ',          &
			co_dom_hep
	endif
	close(u)
	end subroutine write_resolved_config

	! ------------------------------ !

	function round_trip_decimal(val) result(s)
	! The decimal spelling of a binary64 value that a list-directed read
	! returns to the same bits. 17 significant decimal digits suffice for
	! binary64; the 18 printed here are within the field, and the field is
	! wide enough for the sign of the value and a three-digit exponent of
	! either sign (ES24.17E3 is not: it fills with asterisks for a negative
	! value). Both signs, zero, the exponent range and the parse back to the
	! same bits have been verified.
	real*8, intent(in) :: val
	character(len=26)  :: s
	write(s,'(ES26.17E3)') val
	s = adjustl(s)
	end function round_trip_decimal

	! End of module
	end module setup_report

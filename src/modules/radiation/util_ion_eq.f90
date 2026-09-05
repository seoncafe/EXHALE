   module utils_ion_eq

   use global_parameters
   use species_table, only: n_mion, n_mphot, n_melem,                 &
                            mion_isphot, mion_iphot, mion_ethr,       &
                            mion_z2, mion_elem, mion_iscool,          &
                            mion_stage, mion_name, mion_fsp,          &
                            melem_Z, melem_top, im_FeII,              &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,          &
                            isp_OH, isp_H2O, isp_CO
   use utils
   use Cooling_Coefficients      ! Various functions for cooling coefficients
   ! Miller+2013 H3+ infrared cooling, emitted and net of the lower-atmosphere field
   use h3p_cooling, only: h3p_cooling_rate, h3p_net_cooling_rate
   ! H2 quadrupole/magnetic-dipole lines and the H2O and CO bands, net of the
   ! infrared field of the lower atmosphere (`Molecular IR bands`)
   use molecular_reaction_heat, only: molecular_chemical_heating
   use molecular_infrared_cooling, only: h2_line_net_cooling_rate,     &
                            h2o_band_net_cooling_rate,                 &
                            co_band_net_cooling_rate,                  &
                            molecular_planck_cross_section
   use Cross_sections, only: sigma, sigma_HeI, sigma_H2,              &
                             metal_photoion_sigma  ! sigma_H(E,Z),
                                        ! sigma_HeI(E), sigma_H2(E),
                                        ! metal_photoion_sigma(k,E)
   use composition, only: he_ground_singlet_density
   ! D0(H2) in eV. The bond energy has one definition in this code
   ! (mol_rates.f90); the neutral-dissociation heat below reads it from
   ! there rather than writing 4.478 of its own.
   use mol_rates, only: h2_dissociation_energy_eV
   ! Degradation of a fast photoelectron in a partly ionized H/H2/He gas:
   ! Shull & van Steenberg (1985) for the ionization energy budget, Dalgarno,
   ! Yan & Liu (1999) for everything that depends on molecular hydrogen.
   use h2_vibrational_relaxation, only: h2_vibrational_heat_fraction,     &
                                       h2_energy_per_bound_fluorescence_erg, &
                                       h2_energy_per_bound_fluorescence_eV
   use electron_energy_degradation, only: photoelectron_partition_t,      &
                            photoelectron_energy_partition,                 &
                            photoelectron_shares, n_abs_fixed,              &
                            iabs_HI, iabs_HeI, iabs_HeII, iabs_HeTR,     &
                            iabs_H2, iabs_H2_di, dissoc_ion_per_H2p
   ! H2 Lyman-Werner photodissociation, for the heating breakdown diagnostic
   use water_photolysis, only: n_fuv_band, ib_LW, ib_B1, ib_B2,      &
                       ib_B3, ib_B4,                                    &
                       water_photolysis_rate, hydroxyl_photolysis_rate, &
                       fuv_band_optical_depth,                          &
                       heat_per_water_dissociation,                     &
                       heat_per_hydroxyl_dissociation
   use lyman_werner_photodissociation, only:                          &
                            lyman_werner_dissociation_rate,             &
                            e_lw_fragment_erg,                          &
                            h2_lw_dissociation_per_pump,                &
                            h2_lw_dissociation_per_absorbed_photon,     &
                            h2_self_shielding_level_resolved,           &
                            h2_shield_max_column,                       &
                            h2_band_equivalent_width,                   &
                            h2_doppler_parameter
   use omp_lib                   ! OMP libraries
	
	! Move here photoionization and photoheating
	! Two subroutines for H and He+H
	
   implicit none
	! eval_cool sub-block timers (EXHALE_PROFILE=1 via ec_prof_on); diagnostic only
	logical, save :: ec_prof_on = .false.
	real*8,  save :: ec_t(5) = 0.0d0, ec_t0 = 0.0d0
	character(len=40), parameter :: ec_name(5) = (/ 'recomb/ion/coex coefficient fits        ', &
	   'bremsstrahlung Gaunt                    ', 'metal rec/ion/cool tables               ', &
	   'fine-structure transfer + CNO           ', 'H3+ / molecular IR                      ' /)

	! ----- He I 2^1S -> 1^1S two-photon continuum ----- !
	! The 2^1S term decays by emitting a photon PAIR summing to 20.62 eV, with
	! the Drake, Victor & Dalgarno (1969) spectral shape peaking at half that.
	! The part above the H I edge is therefore neither a round fraction of the
	! pair nor flat within the window; integrating the tabulated shape over
	! 13.598-20.62 eV gives 0.5564 H-ionizing photons per decay carrying
	! 8.9646 eV, i.e. a mean in-band photon energy of 16.110 eV and hence
	! 2.512 eV of photoelectron energy per photon. A flat distribution over the
	! same window would instead put the mean at 17.109 eV and the photoelectron
	! at 3.511 eV, overestimating the deposited energy by 40%.
	real*8, parameter :: f_2q_HeI  = 0.5564d0   ! ionizing photons per 2^1S decay
	real*8, parameter :: Ee_2q_HeI = 2.512d0    ! photoelectron energy [eV]
	! Photoelectron energies of the He cascade exits that are single lines:
	! 584 A resonance (21.2-13.6) and the 2^3S 19.8 eV line (19.8-13.6).
	real*8, parameter :: Ee_584_HeI  = 7.6d0
	! Singlet-excited captures: 2/3 go to 2^1P -> 584 A (always ionizing),
	! 1/3 to 2^1S -> two-photon.
	real*8, parameter :: f_sing_HeI  = (2.0d0 + f_2q_HeI)/3.0d0
	real*8, parameter :: Ee_sing_HeI = (2.0d0*Ee_584_HeI                    &
	                                    + f_2q_HeI*Ee_2q_HeI)/3.0d0
	! Low-density channel-weighted average over the case-B cascade exits
	! (19.8 eV line / 584 A / two-photon), used in atomic mode.
	real*8, parameter :: Ee_casc_HeI = 0.75d0*6.2d0 + 0.17d0*Ee_584_HeI     &
	                                   + 0.08d0*Ee_2q_HeI

	! ----- Representative photon energies of the He recombination channels -----
	! Each channel of he_rec_coupling emits at one energy, at which the
	! absorbers (H I, H2 and, above 24.6 eV, He I) compete for the photon.
	real*8, parameter :: E_gnd_HeI  = 24.6d0   ! ground-capture continuum edge
	real*8, parameter :: E_584_HeI  = 21.2d0   ! 2^1P -> 1^1S resonance line
	real*8, parameter :: E_19_HeI   = 19.8d0   ! 2^3S -> 1^1S line
	! Mean energy of the 2^1S two-photon photons that lie above the H I edge:
	! the continuum is not a line, and this single energy stands for it. The
	! shape (Drake, Victor & Dalgarno 1969) is not carried in the code, only
	! the two integrals f_2q_HeI and Ee_2q_HeI over the H I window, so the
	! energy below is 13.598 + Ee_2q_HeI and cannot be re-integrated over the
	! narrower H2 window; see EeH2_2q_HeI.
	real*8, parameter :: E_2q_HeI   = 13.598d0 + Ee_2q_HeI      ! 16.110 eV
	! Cascade-averaged exit energy used by the atomic (case-B) branch, the same
	! 0.75/0.17/0.08 weighting that defines Ee_casc_HeI.
	real*8, parameter :: E_casc_HeI = 13.598d0 + Ee_casc_HeI    ! 19.741 eV
	! The photoelectron each of those photons leaves in H2 (threshold e_th_H2 =
	! 15.4 eV), the H2 counterparts of 11.0 / 7.6 / 6.2 / Ee_2q_HeI /
	! Ee_casc_HeI.
	real*8, parameter :: EeH2_gnd_HeI  = E_gnd_HeI  - e_th_H2   ! 9.200 eV
	real*8, parameter :: EeH2_584_HeI  = E_584_HeI  - e_th_H2   ! 5.800 eV
	real*8, parameter :: EeH2_19_HeI   = E_19_HeI   - e_th_H2   ! 4.400 eV
	! APPROXIMATE: the H2-ionizing photon COUNT of the two-photon continuum is
	! the integral of the same shape over 15.4-20.62 eV, which is smaller than
	! f_2q_HeI (the integral over 13.598-20.62 eV); the shape is not in the
	! code, so f_2q_HeI is used for H2 as well and this sub-channel's H2 share
	! is overestimated. It is 1/3 of one of five channels, and R2_2q = 1.2 is
	! the smallest H2/H I cross-section ratio of the set.
	real*8, parameter :: EeH2_2q_HeI   = E_2q_HeI   - e_th_H2   ! 0.710 eV
	real*8, parameter :: EeH2_casc_HeI = E_casc_HeI - e_th_H2   ! 4.341 eV

	! ----- He recombination coupling diagnostic -----
	! Cells in which the He recombination photons ionize H I faster than the
	! stellar field does by more than this factor. The coupling is a CORRECTION
	! to the direct photoionization integral, so it dominating by three decades
	! is a statement that the on-the-spot budget has run away -- which is what
	! the unweighted 1/n_HI of the earlier form did in a fully ionized cell
	! (docs/supersonic_molecular_base.md section 12). Counted over cells and
	! steps, like the positivity counters, and reported at the end of the run
	! only when it is not zero.
	real*8, parameter :: he_rec_dominant_ratio = 1.0d3
	integer :: n_cells_he_rec_photoionization_dominant = 0
	real*8  :: he_rec_photoionization_ratio_max = 0.0d0

	contains

	! Incident stellar flux of FUV band ib at the planet [erg cm^-2 s^-1],
	! in the band order of water_photolysis (LW, B1, B2 = Ly-alpha, B3, B4).
	! Single definition, read by the equilibrium solve and by write_output.
	!
	! THE FIRST BAND IS THE LYMAN-WERNER INTERVAL and carries the flux of the
	! existing "Stellar LW flux" key, because 912-1110 A is one interval with
	! one incident flux whether the absorber is H2 in lines or H2O and OH in
	! a continuum. Supplying it twice -- once for H2, once inside a band
	! reaching across it -- would count the energy of the interval twice, and
	! that is what B1 starting at 1110 A removes.
	!
	! B2 IS THE INCIDENT STELLAR LY-ALPHA FLUX, NOT A SOLVED FIELD, and that
	! is the weakest part of the band treatment. The H2O and OH attenuation
	! applied to it is their own continuum, as for every other band; what is
	! NOT applied is the H I resonance scattering that actually decides how
	! much Ly-alpha reaches a molecular base under an ionized wind. lya_rt.f90
	! solves that field for the n = 2 pumping, and section 2.6 of
	! docs/a2_oxygen_option_design.md says the solved J_Lya(r) is the field
	! B2 should use; wiring it is not done here, because the conversion of a
	! mean intensity into a photodissociation rate is a different quantity
	! from the one lya_rt returns and getting it wrong is worse than leaving
	! it out. So B2's photolysis rate is an UPPER BOUND, and a run that finds
	! B2 dominating its oxygen chemistry should say so. The band is reported
	! separately in output/FUV_bands.txt precisely so that it can be read off.
	!
	! Note that setting "Stellar Lya flux" for the n = 2 pumping therefore
	! also drives B2, and setting "Stellar LW flux" for the H2 network also
	! drives the oxygen photolysis of the LW band; the setup report states
	! all five band fluxes at startup so this is visible rather than implicit.
	! The dayside dilution is applied HERE and not at the keys, so that
	! F_LW_star and the rest keep their stated meaning -- the band flux at
	! the planet's orbit -- in the setup report and the resolved dump, and
	! so that a flux integrated from the spectrum file is diluted by the
	! same convention as a stated one.  Before section 150 these five bands
	! carried no dilution at all while the same run halved its XUV and its
	! stellar Ly-alpha beam, so the Lyman-Werner rate stood a factor
	! 1/dayside_dilution() above the run's own convention.
	double precision function fuv_band_flux(ib) result(F)
	integer, intent(in) :: ib
	select case (ib)
	case (ib_LW)
		F = F_LW_star
	case (ib_B1)
		F = F_FUV_B1
	case (ib_B2)
		F = F_Lya_star
	case (ib_B3)
		F = F_FUV_B3
	case (ib_B4)
		F = F_FUV_B4
	case default
		F = 0.0d0
	end select
	F = dayside_dilution()*F
	end function fuv_band_flux

	! ------------------------------------------------------------- !

	! The 912-2304 A photon field of the molecular layer, solved once for
	! ALL of its absorbers: H2 in the Lyman and Werner lines, H2O and OH in
	! their continua.  Single definition, called by the equilibrium sweep
	! and by the heating-breakdown dump, so the two cannot drift apart.
	!
	! WHY ONE ROUTINE.  Over 912-1110 A the three absorbers share one beam.
	! An H2O molecule cannot absorb a photon an H2 line has already taken,
	! and the H2 pumping rate is reduced by whatever continuum sits above
	! it -- that continuum term is the exp(-tau) of Draine & Bertoldi (1996)
	! eq. (40), identically 1 before this option existed and not 1 now.
	! Solving the two sides in two places is how the energy of the interval
	! comes to be counted twice, so they are solved here together.
	!
	! The beam's transmission down to a face is
	!
	!     T = (1 - A) exp(-tau_c) ,
	!     A(j) = int_j^top sigma_pump f_shield(N_H2) n_H2 dr ,
	!
	! and across one cell the identity
	!
	!   T_out - T_in = e^{-tau_out} [ (1-A_out)(1 - e^{-dtau}) + dA e^{-dtau} ]
	!
	! is exact, so the photons the cell removes split into a continuum share
	! and a line share with nothing left over (water_photolysis.f90 sec. 3).
	! The continuum share is (1-A_out) times what it would be alone, and the
	! line share is exp(-tau_c at the cell's own depth) times what IT would
	! be alone -- which is why the two factors below are evaluated at
	! different faces and that is not an inconsistency.
	!
	! Every output is zero, or the neutral value 1, for a run that does not
	! carry the absorber in question, so a molecular run without the oxygen
	! chemistry gets exactly the rate it got before.
	subroutine fuv_lw_photon_field(nH2, nH2O, nOH, T_K, nH_nuc,           &
	                               NH2col, NH2Ocol, NOHcol,               &
	                               f_shield, tr_lines, tau_b,             &
	                               k_lw, p_lw_single, p_lw_absorbed,      &
	                               j_h2o, j_oh, a_lines_max,              &
	                               col_over_overlap)
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nH2, nH2O, nOH, T_K
	! Total hydrogen NUCLEUS density, the density axis of the self-shielding
	! table (CLOUDY's hden).  It enters nothing else here.
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nH_nuc
	real*8, dimension(1-Ng:N+Ng), intent(out) :: NH2col, NH2Ocol, NOHcol
	real*8, dimension(1-Ng:N+Ng), intent(out) :: f_shield, tr_lines, k_lw
	! The two branching ratios of one Lyman-Werner absorption, both from the
	! same level-resolved table as k_lw and both varying with depth.
	!   p_lw_single    dissociations per PUMP with every fluorescent photon
	!                  free to escape.  The number of fluorescent decays per
	!                  dissociation is (1 - p)/p with THIS p: a trapped decay
	!                  and the re-absorption that follows it move one molecule
	!                  down and another up and deposit nothing, so trapping
	!                  cancels out of that count.  The vibrational heat return
	!                  uses it.
	!   p_lw_absorbed  dissociations per band photon REMOVED FROM THE BEAM.
	!                  Only fresh pumps take photons out of the stellar beam,
	!                  so the photon ledger of the band uses this one.
	! Both are 0 for a run with no band flux, like k_lw.
	real*8, dimension(1-Ng:N+Ng), intent(out) :: p_lw_single, p_lw_absorbed
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(out) :: tau_b
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(out) :: j_h2o, j_oh
	real*8, intent(out) :: a_lines_max
	! Largest ratio of the star-ward H2 column to the column at which the
	! self-shielding table stops being a value and becomes an upper bound
	! (lyman_werner.f90 sec. 2d).  Above 1 the deepest cells are in the
	! line-overlap regime the table does not carry.
	real*8, intent(out) :: col_over_overlap
	real*8, dimension(1-Ng:N+Ng) :: a_lines
	real*8  :: tau_out, dtau, tr_out
	integer :: j, ib

	NH2col   = 0.0d0
	NH2Ocol  = 0.0d0
	NOHcol   = 0.0d0
	f_shield = 0.0d0
	tr_lines = 1.0d0
	k_lw     = 0.0d0
	p_lw_single   = 0.0d0
	p_lw_absorbed = 0.0d0
	tau_b    = 0.0d0
	j_h2o    = 0.0d0
	j_oh     = 0.0d0
	a_lines_max = 0.0d0
	col_over_overlap = 0.0d0

	! ---- H2 lines: self-shielding and the fraction of the LW band they
	! take out of the shared beam.  Both are functions of the H2 column
	! alone, so they exist whether or not the oxygen chemistry is on.
	!
	! A is the summed dimensionless equivalent width of the pumping lines,
	! int_0^N sigma_pump f_shield dN', which Draine & Bertoldi (1996)
	! give in closed form as their eq. (39) -- h2_band_equivalent_width.
	! Taking it from there rather than re-integrating on the radial grid
	! makes the line transmission independent of the grid and exactly
	! consistent with the fit whose integral it is.  A is an equivalent
	! width and not an optical depth, so the transmission of the beam is
	! 1 - A, not exp(-A).  The clamp is a floor on that transmission:
	! eq. (39) has an asymptote slightly above 1 (lyman_werner.f90), which
	! is fit slack, not physics.
	!
	! THE TWO LINES BELOW DELIBERATELY COME FROM DIFFERENT PLACES, and that
	! is the point of lyman_werner.f90 sec. 2: f_shield is the suppression of
	! the photodissociation RATE, taken from the level-resolved CLOUDY table
	! because neither published fit follows that calculation at the
	! temperature of this layer; A is the SHARE of the band the lines remove,
	! which exists in closed form only as the integral of DB96's own eq. (37),
	! and whose normalization the same level-resolved calculation confirms to
	! 2 per cent.  f_shield is what the run reports in
	! output/Lyman_Werner.txt, i.e. the factor the rate actually carries.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		call calc_column_dens_one(nH2, NH2col)
		do j = 1-Ng,N+Ng
			f_shield(j) = h2_self_shielding_level_resolved(NH2col(j),     &
			                            T_K(j), nH_nuc(j))
			a_lines(j) = min(h2_band_equivalent_width(NH2col(j),          &
			                 h2_doppler_parameter(T_K(j))), 1.0d0)
			tr_lines(j) = 1.0d0 - a_lines(j)
			col_over_overlap = max(col_over_overlap,                       &
			                       NH2col(j)/h2_shield_max_column())
		enddo
		a_lines_max = maxval(a_lines)
	endif

	! ---- H2O and OH continua, on every band.
	if (thereis_oxychem) then
		call calc_column_dens_one(nH2O, NH2Ocol)
		call calc_column_dens_one(nOH,  NOHcol)
		do ib = 1,n_fuv_band
			do j = 1-Ng,N+Ng
				tau_b(j,ib) = fuv_band_optical_depth(ib, NH2Ocol(j),      &
				                                     NOHcol(j))
			enddo
		enddo
	endif

	! ---- H2 photodissociation, with the continuum of the SAME interval
	! attenuating it: DB96 eq. (40) in structure, the level-resolved
	! self-shielding table in it.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		do j = 1-Ng,N+Ng
			k_lw(j) = lyman_werner_dissociation_rate(                     &
			              fuv_band_flux(ib_LW),                           &
			              NH2col(j), T_K(j), nH_nuc(j), tau_b(j,ib_LW))
			p_lw_single(j)   = h2_lw_dissociation_per_pump(NH2col(j),     &
			                       T_K(j), nH_nuc(j))
			p_lw_absorbed(j) = h2_lw_dissociation_per_absorbed_photon(    &
			                       NH2col(j), T_K(j), nH_nuc(j))
		enddo
	endif

	! ---- H2O and OH photodissociation, band by band.  tau_b(j,ib) is the
	! depth at the INNER face of cell j (the column at j already contains
	! cell j), so the star-ward face of cell j carries tau_b(j+1,ib) and the
	! line transmission tr_lines(j+1).  The rate is the MEAN over the cell,
	! which makes the photons the model absorbs exactly the photons the beam
	! loses at any grid spacing (water_photolysis.f90).
	if (thereis_oxychem) then
		do ib = 1,n_fuv_band
			do j = 1-Ng,N+Ng
				if (j .lt. N+Ng) then
					tau_out = tau_b(j+1,ib)
					tr_out  = tr_lines(j+1)
				else
					tau_out = 0.0d0   ! nothing above the top cell
					tr_out  = 1.0d0
				endif
				if (ib .ne. ib_LW) tr_out = 1.0d0
				dtau = max(tau_b(j,ib) - tau_out, 0.0d0)
				j_h2o(j,ib) = water_photolysis_rate(fuv_band_flux(ib),    &
				                            ib, tau_out, dtau, tr_out)
				j_oh(j,ib)  = hydroxyl_photolysis_rate(fuv_band_flux(ib), &
				                            ib, tau_out, dtau, tr_out)
			enddo
		enddo
	endif

	end subroutine fuv_lw_photon_field

	
	! ------------------------------------------------------------- !

	! The share of a photon's energy that a photoionization event turns into
	! heat: the photoelectron's (E - I)/E, or the whole photon (1) under the
	! comparison option "Photoelectron heating: full" (see parameters.f90).
	elemental double precision function photoelectron_share(e_th, e) result(f)
	real*8, intent(in) :: e_th, e
	if (photoheat_photon_fraction .gt. 0.0d0) then
		f = photoheat_photon_fraction
	else
		f = 1.0d0 - e_th/e
	endif
	end function photoelectron_share

	subroutine PH_heat_H(nhi, xion, P_HI,heat,q, heat_of_one_HI)
	! Computes photoionization rates and heating rates for a PURE HYDROGEN
	! atmosphere (thereis_He = .false.); the He/H+metals path is PH_heat_HHe.
	! The helium column densities below are held at zero so that the shared
	! calc_column_dens can be called unchanged.
	!
	! What is returned, and at which state. P_HI and heat_of_one_HI are
	! properties of the RADIATION FIELD and of the photoelectron partition:
	! the frequency integrals with the absorber density divided out, i.e.
	! the photoionization rate [s^-1] and the photoheating rate [erg s^-1] of
	! ONE H I atom sitting in the attenuated field of this cell. `heat` is
	! their contraction with the composition passed in, heat = n_HI *
	! heat_of_one_HI, and so describes THAT composition and no other; a
	! caller that changes the composition afterwards must re-form it from
	! heat_of_one_HI (photoheating_of_composition does exactly this).
	! `q` is the heating efficiency of the same composition.

	integer :: j

	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng),intent(in) :: xion

	! Secondary-ionization scratch (H-only: no He I and no H2 channel).
	! fiHI is the branching per H I atom [cm^3] and Psec_HI the rate it gives
	! [1/s]; fh_v is the heat fraction of the photoelectron. Both are resolved
	! on the spectral grid through the absorber's excess energy e_v - e_th_HI.
	real*8, dimension(Nl) :: acc_secHI, fhv, fh_v, fiHI
	real*8, dimension(Nl) :: fiHeI_dum, fiH2_dum
	type(photoelectron_partition_t) :: pep
	real*8 :: Psec_HI
	! SvS85 coupling applied only when enabled AND staged on (see EXHALE_main).
	logical :: sec_on

	! Dummy zero (set below; the number of cells is a runtime value, so these
	! cannot be named constants)
	real*8, dimension(1-Ng:N+Ng) :: nhei, nheii, nheiii, nheiTR
	real*8, dimension(1-Ng:N+Ng) :: N15, N2, NTR

   real*8 :: PIR_1		            ! Photoionization rates
   real*8 :: Hea_1 		        ! Heating rates  
	real*8 :: q_abs           	! Absorbed energy
	
	! Integral variables
	real*8, dimension(Nl) :: tauE			
	real*8, dimension(Nl) :: int_f,int_1,int_q,int_H

	! Column densities                                 
   real*8, dimension(1-Ng:N+Ng) ::  N1 
	
   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  P_HI 
      
   ! Heating efficiency 
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q
      
   ! Heating rate at the composition passed in [erg cm^-3 s^-1]
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat

   ! Photoheating rate of ONE H I atom [erg s^-1] (see the header).
   real*8, dimension(1-Ng:N+Ng),intent(out),optional ::  heat_of_one_HI
   real*8, dimension(1-Ng:N+Ng) ::  h1_HI

	!----------------------------------!

	sec_on = use_sec_ion .and. sec_ion_active

	nhei   = 0.0
	nheii  = 0.0
	nheiii = 0.0
	nheiTR = 0.0

	! Evaluate the column density
	call calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)

   ! Evaluate photoionization rates and photoheating rates
	do j = 1-Ng,N+Ng

		! Initialization of integrands
		Hea_1  = 0.0
		PIR_1  = 0.0
		q_abs  = 0.0

		tauE = s_hi*N1(j)*1.0e-18	

		! Initial integrands
		int_f = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		! Secondary-ionization energy partition for this cell. Hydrogen is the
		! only target on this path, so this is the n_HeI -> 0, n_H2 -> 0 limit
		! of photoelectron_energy_partition: the whole ionization energy goes to
		! hydrogen, and fiHI comes back per H I atom [cm^3]. With no molecular
		! hydrogen the H2 terms are identically absent and the heat fraction is
		! Dalgarno's H-He heating efficiency, closed at x = 1.
		if (sec_on) then
			call photoelectron_energy_partition(xion(j), nhi(j), 0.0d0,     &
			                                  0.0d0, 0.0d0, 0.0d0, pep)
			call photoelectron_shares(pep, iabs_HI, fh_v, fiHI,            &
			                          fiHeI_dum, fiH2_dum)
		else
			fh_v = 1.0d0; fiHI = 0.0d0
		endif
		! Heating fraction: fh_v above the E_sec_ion photoelectron threshold, 1
		! (full thermalization) below it. fhv = 1 when the coupling is off, so
		! the heating integrand is bit-identical to the legacy path.
		fhv = 1.0d0
		if (sec_on) fhv = merge(fh_v, 1.0d0, e_v > e_th_HI + E_sec_ion)
		! The heating integrand carries NO density: Hea_1 is the heating rate
		! of one H I atom, and the contraction with the composition is done
		! below. The absorbed-energy integrand keeps its density, because the
		! heating efficiency is a property of the state passed in.
		int_H = int_f*photoelectron_share(e_th_HI,e_v)*fhv*s_hi
		int_1 = int_f*s_hi/e_v
		int_q = int_f*s_hi*nhi(j)

		! Value of integrals
		Hea_1 = sum(int_H*de_v)
		PIR_1 = sum(int_1*de_v)
		q_abs = sum(int_q*de_v)

		! Multiply for the dimensional coefficient
		h1_HI(j)  = Hea_1*1.0e-18
		heat(j)   = h1_HI(j)*nhi(j)
		P_HI(j)   = PIR_1*1.0e-18*erg2eV
		! Add the H I secondary-ionization rate from fast photoelectrons.
		! fiHI is per H I atom, so the integral is already a rate [1/s].
		if (sec_on) then
			acc_secHI = int_f*s_hi*nhi(j)/e_v * &
			     merge(fiHI*(e_v-e_th_HI)/e_th_HI, 0.0d0, e_v > e_th_HI + E_sec_ion)
			Psec_HI = sum(acc_secHI*de_v)*1.0e-18*erg2eV
			P_HI(j) = P_HI(j) + Psec_HI
		endif
		! q_abs = int F sigma n_HI dE is the energy absorbed per unit volume
		! and time: it is non-negative for any physical state, and it vanishes
		! together with Hea_1 in a transparent or unilluminated cell. The
		! heating efficiency is undefined there and zero is its physical value
		! (nothing absorbed, nothing deposited). A non-positive q_abs with a
		! non-zero Hea_1 can only come from a negative neutral density, i.e.
		! from an unphysical ionization root; the ionization solve rejects
		! those (ionization_fractions_physical), and the test here keeps the
		! efficiency finite and signed correctly if one ever survives.
		if (q_abs .gt. 0.0d0) then
			q(j)  = Hea_1*nhi(j)/q_abs
		else
			q(j)  = 0.0d0
		endif

	enddo

	if (present(heat_of_one_HI)) heat_of_one_HI = h1_HI

	end subroutine PH_heat_H

	! ------------------------------------------------------------- !
	
	subroutine PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,       &
				     P_HI,P_HeI,P_HeII,P_HeITR, P_m,       &
				     heat,q, nh2,P_H2, P_H2_di, P_H2_dd, P_H2_nd,  &
				     f_vibq, e_vibq,                                &
				     heat_chan,                                     &
				     heat_of_one_HI, heat_of_one_HeI,               &
				     heat_of_one_HeII, heat_of_one_HeTR,            &
				     heat_of_one_H2, heat_of_one_mion)
	! Computes photoionization rates and heating rates for an
	!	atmosphere composed of H, He, and (optionally) metals.
	! Metal ion densities arrive as nm(:,1:n_mion) in canonical
	! species_table order; the photoionization rates for each ion leave as
	! P_m(:,1:n_mion). Only photo-ionizable metal ions (mion_isphot)
	! contribute to opacity, photoheating, absorbed energy, and have a
	! nonzero P_m; top-stage ions (C III, N III, O III, Mg III) are inert.
	!
	! What is returned, and at which state. The attenuated field, the
	! photoionization rates P_* and the heat_of_one_* rates are properties of
	! the RADIATION FIELD (through the columns and the photoelectron
	! partition of the composition passed in), not of the abundance of the
	! absorber that receives them: heat_of_one_X is the frequency integral
	! with the density of X divided out, i.e. the photoheating rate
	! [erg s^-1] of ONE particle of X in this cell's field. `heat` and
	! `heat_chan` are their contraction with the composition passed in
	! (photoheating_of_composition, the single definition of that sum), and
	! so describe THAT composition and no other. A caller that changes the
	! composition afterwards -- an equilibrium sweep, an advection
	! correction -- must re-form the heating from heat_of_one_* and the new
	! densities rather than keep this `heat`. `q` is the heating efficiency
	! of the composition passed in.

	integer :: i,j,k

	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhei,nheii
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng),intent(in) :: xion
	! Optional H2 (molecular extension): adds the H2 opacity,
	! photoionization rate P_H2 [1/s] and photoelectric heating using the
	! Yan+1998 cross section s_h2 (threshold e_th_H2 = 15.4 eV).
	real*8, dimension(1-Ng:N+Ng), intent(in),  optional :: nh2
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: P_H2
	! The DISSOCIATIVE part of P_H2, H2 + hv -> H + H+ + e- (18.08 eV): a
	! subset of P_H2, not an addition to it. The H2 destruction rate is
	! P_H2; the H2+ production rate is P_H2 - P_H2_di, and P_H2_di is a
	! proton and an H atom source instead. Present only when P_H2 is.
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: P_H2_di
	! The DOUBLE-ionization part of P_H2, H2 + hv -> H+ + H+ + 2e-
	! (51.4 eV): again a subset of P_H2, and also a subset of nothing else
	! -- it is disjoint from P_H2_di. One event releases TWO protons, so
	! the proton source it feeds carries a factor 2 while the H2 it
	! consumes is one. Zero unless a double-ionization model is selected.
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: P_H2_dd
	! The NEUTRAL-dissociation part of P_H2, H2 + hv -> H + H: the
	! absorptions that make no ion at all. A subset of P_H2 like the two
	! above, disjoint from both, and the only one of the four channels
	! that leaves no photoelectron. Zero unless h2_neutral_dissociation.
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: P_H2_nd
	! Fraction of an H2 vibrational excitation that is collisionally
	! de-excited into heat rather than radiated (h2_vibrational_relaxation).
	! It needs the local temperature, which this routine does not carry, so
	! the caller evaluates it. Absent => the two vibrational heat channels
	! are off, which is also what n_H2 = 0 gives.
	real*8, dimension(1-Ng:N+Ng), intent(in),  optional :: f_vibq
	! Mean internal energy [eV] left in X by one B or C fluorescence, the
	! same h2_energy_per_bound_fluorescence_eV the Lyman-Werner heat uses.
	! Supplied by the molecular path together with f_vibq; where it is
	! absent f_vibq is absent too and the channel it feeds is zero.
	real*8, dimension(1-Ng:N+Ng), intent(in),  optional :: e_vibq
	real*8, dimension(1-Ng:N+Ng) :: nheiS
   real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2,NTR
   real*8, dimension(1-Ng:N+Ng,n_mphot) :: Nm_col

   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR     ! Photoionization rates (H/He)
   real*8 :: PIR_H2                        ! H2 photoionization rate
   real*8 :: PIR_H2_di                     ! ... its dissociative part
   real*8 :: PIR_H2_dd                     ! ... its double-ionization part
   real*8 :: PIR_H2_nd                     ! ... its neutral-dissociation part
   real*8, dimension(1-Ng:N+Ng) :: NH2col  ! H2 column density
   real*8, dimension(Nl) :: int_h2, int_h2_di, int_h2_dd, int_h2_nd
   real*8 :: q_abs                         ! Absorbed energy
   real*8 :: Pm_loc(n_mion)                ! Metal photoion. rates in each cell

	! Integral variables
	real*8, dimension(Nl) :: tauE,tau_m
	real*8, dimension(Nl) :: int_f,int_1,int_15,int_2,int_TR
	real*8, dimension(Nl) :: int_m
	real*8, dimension(Nl) :: int_q,acc_q
	! Secondary-ionization scratch. fiHI/fiHeI/fiH2 are the branchings per
	! target particle [cm^3] (photoelectron_energy_partition) and
	! Psec_HI/Psec_HeI/Psec_H2 the rates they give [1/s]. fh_v is the heat
	! fraction of the photoelectron of the absorber currently being summed,
	! resolved on the spectral grid through that absorber's own excess
	! energy e_v - e_th; it is recomputed inside each absorber block.
	real*8, dimension(Nl) :: acc_secHI,acc_secHeI,acc_secH2,fhv,fh_v
	real*8, dimension(Nl) :: fiHI,fiHeI,fiH2
	type(photoelectron_partition_t) :: pep
	real*8 :: Psec_HI,Psec_HeI,Psec_H2,Psec_H2_di
	! Secondary-ionization coupling applied only when enabled AND staged on
	! (see EXHALE_main). mol_sec additionally requires the run to carry H2,
	! which is what makes the molecular target channel cost nothing, and
	! change nothing, in a run without molecular hydrogen.
	logical :: sec_on, mol_sec
	! Local copy of f_vibq, zero when the caller supplied none.
	real*8, dimension(1-Ng:N+Ng) :: fvq, evq
	! D0(H2) in eV, hoisted out of the cell loop: it is a constant of the
	! molecule, read from the one place this code defines it.
	real*8 :: D0_H2_eV

   ! Photo ionization rates
	real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HI
   real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) ::  P_m

   ! Heating efficiency
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q

   ! Heating rate at the composition passed in [erg cm^-3 s^-1]
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat

	! Photoheating rate of ONE particle of each absorber [erg s^-1]: H I,
	! He I in its ground singlet, He II, He 2^3S, H2 (its four final-state
	! channels summed, since one absorption destroys one molecule whichever
	! channel it takes), and each metal ion in canonical order. See the
	! header for why these, and not `heat`, are what a caller with a
	! different composition needs.
	real*8, dimension(1-Ng:N+Ng),intent(out),optional :: heat_of_one_HI,   &
	         heat_of_one_HeI, heat_of_one_HeII, heat_of_one_HeTR,          &
	         heat_of_one_H2
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(out),optional ::            &
	         heat_of_one_mion
	! The same quantities as this routine's own working arrays, and the
	! absorbed energy of the composition passed in, which is what `q` is
	! measured against.
	real*8, dimension(1-Ng:N+Ng) :: h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2
	real*8, dimension(1-Ng:N+Ng,n_mion) :: h1_m
	real*8, dimension(1-Ng:N+Ng) :: q_abs_cell
	! Zero H2 density for the contraction below when the caller tracks no
	! molecular hydrogen (h1_H2 is then zero as well).
	real*8, dimension(1-Ng:N+Ng) :: nh2_loc
	real*8, dimension(Nl) :: acc_mion
	real*8 :: h1m_loc(n_mion)

	! Optional photoheating breakdown by absorber (cgs erg cm^-3 s^-1).
	! Columns: 1 H I, 2 He I (singlet ground), 3 He II, 4 He 2^3S,
	! 5 H2 (molecular, 0 when absent), 6 metals (sum over photo-ionizable
	! metal ions). `heat` IS the sum of these six columns: both come from the
	! same contraction of heat_of_one_* with the composition, so the split is
	! exact rather than a separately accumulated approximation of it.
	real*8, dimension(1-Ng:N+Ng,6),intent(out),optional :: heat_chan
	! Frequency integrands of each absorber, with the absorber density
	! divided out (that density enters only in the contraction below).
	real*8, dimension(Nl) :: acc_HI,acc_HeI,acc_HeII,acc_HeTR,acc_H2

	!----------------------------------!

	! Ground singlet n(1^1S): the summed neutral helium when the metastable
	! is not tracked, and the clamped difference n(He I) - n(2^3S) when it
	! is. composition.f90 is the only place that difference is taken.
	nheiS = nhei
	if (thereis_HeITR) nheiS = he_ground_singlet_density(nhei, nheiTR)

	! Evaluate the column density
	call calc_column_dens(nhi,nheiS,nheii,nheiTR,N1,N15,N2,NTR)
	call calc_column_dens_metals(nm, Nm_col)

	! H2 column (molecular)
	if (present(nh2)) call calc_column_dens_one(nh2, NH2col)

	! Metal-free P_m entries (top-stage ions) stay zero
	P_m = 0.0

	D0_H2_eV = h2_dissociation_energy_eV()
	sec_on  = use_sec_ion .and. sec_ion_active
	mol_sec = sec_on .and. present(nh2) .and. present(P_H2)
	fvq     = 0.0d0
	evq     = 0.0d0
	if (present(f_vibq)) fvq = f_vibq
	if (present(e_vibq)) evq = e_vibq

    !----------------------------------!
	!$OMP PARALLEL DO &
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,P_HeITR,P_m, sec_on, mol_sec, fvq, evq, &
	!$OMP          D0_H2_eV ) &
	!$OMP PRIVATE ( PIR_1,PIR_15,PIR_2,PIR_TR,PIR_H2,int_h2,                     &
	!$OMP           PIR_H2_di,int_h2_di,Psec_H2_di,                             &
	!$OMP           PIR_H2_dd,int_h2_dd,PIR_H2_nd,int_h2_nd,                    &
	!$OMP           Pm_loc,h1m_loc,acc_mion,                                     &
	!$OMP           int_1,int_15,int_2,int_TR,int_m,                            &
	!$OMP           acc_secHI,acc_secHeI,acc_secH2,fhv,fh_v,pep,               &
	!$OMP           fiHI,fiHeI,fiH2,Psec_HI,Psec_HeI,Psec_H2,                  &
	!$OMP           acc_HI,acc_HeI,acc_HeII,acc_HeTR,acc_H2,                    &
	!$OMP           int_f,int_q,acc_q,tauE,tau_m,q_abs,i,k,j)

	do j = 1-Ng,N+Ng

      PIR_1   = 0.0
      PIR_15  = 0.0
      PIR_2   = 0.0
      PIR_TR  = 0.0
      q_abs   = 0.0

		! The photoheating integrands default to zero so the He 2^3S and H2
		! rates stay 0 in cells/runs where those absorbers are absent.
		acc_HeTR = 0.0d0
		acc_H2   = 0.0d0
		h1m_loc  = 0.0d0

		! Calculate optical depth. Accumulate the metal block separately
		! in iphot order before applying the 1e-18 factor, matching the
		! original (sum)*1e-18 association exactly.
		tauE = (s_hi*N1(j) + s_hei*N15(j) + s_heii*N2(j))*1.0e-18
		if (thereis_HeITR) tauE = tauE + s_heiTR*NTR(j)*1.0e-18
		if (present(nh2)) tauE = tauE + s_h2*NH2col(j)*1.0e-18
		tau_m = 0.0
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			tau_m = tau_m + sigma_tab(:,k)*Nm_col(j,k)
		enddo
		tauE = tauE + tau_m*1.0e-18

		! Calculate photoionization integrals
		int_f  = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		int_1  = int_f*s_hi/e_v
		int_15 = int_f*s_hei/e_v
		int_2  = int_f*s_heii/e_v
		if (thereis_HeITR) int_TR =  int_f*s_heiTR/e_v
		if (present(nh2)) then
			int_h2    = int_f*s_h2/e_v
			int_h2_di = int_f*s_h2_di/e_v
			int_h2_dd = int_f*s_h2_dd/e_v
			int_h2_nd = int_f*s_h2_nd/e_v
		endif

		! Secondary-ionization energy partition for this cell. The split of
		! the ionization energy between the H I, He I and H2 channels is
		! renormalized onto the cell's own neutral densities -- the fits carry
		! Shull & van Steenberg's n(He I)/n(H I) = 0.1 and no H2 at all -- and
		! comes back per target particle, so the absorber loops below build
		! rates in [1/s] directly and nothing is divided by a vanishing
		! neutral density. The H2 share uses the collision weight of Dalgarno,
		! Yan & Liu (1999) eqs. (9) and (10); it and the molecular part of the
		! heat fraction vanish identically when the run carries no H2.
		! Derivation, limits and the parts of the partition that are still
		! composition-blind: electron_energy_degradation.f90.
		if (sec_on) then
			if (present(nh2)) then
				call photoelectron_energy_partition(xion(j), nhi(j),     &
				                       nheiS(j), nh2(j), fvq(j), evq(j),  &
				                       pep)
			else
				call photoelectron_energy_partition(xion(j), nhi(j),     &
				                       nheiS(j), 0.0d0, 0.0d0, 0.0d0, pep)
			endif
		else
			fiHI = 0.0d0; fiHeI = 0.0d0; fiH2 = 0.0d0
		endif
		acc_secHI  = 0.0d0
		acc_secHeI = 0.0d0
		acc_secH2  = 0.0d0

		! Photoheating integrand of each absorber, WITHOUT that absorber's
		! density: what is built here is the heating rate of one particle of
		! it, and the composition enters only in the contraction after the
		! loop. Where a photoelectron energy E0 = e_v - E_th exceeds
		! E_sec_ion, only f_heat(x) of its excess is deposited as heat (fhv)
		! and the balance drives H I / He I secondary ionizations; below the
		! threshold it thermalizes fully. fhv = 1 when the coupling is off.
		! He I triplet photoionization (threshold e_th_HeTR = 4.8 eV)
		! deposits its photoelectron energy here as well, consistently with
		! its opacity and its P_HeITR rate. The secondary-ionization
		! integrands below DO carry the absorber densities: they are rates of
		! the field the entry composition makes, not one-particle quantities.
		fhv = 1.0d0
		if (sec_on) then
			call photoelectron_shares(pep, iabs_HI, fh_v, fiHI, fiHeI, fiH2)
			fhv = merge(fh_v, 1.0d0, e_v > e_th_HI + E_sec_ion)
		endif
		acc_HI = photoelectron_share(e_th_HI,e_v)*fhv*s_hi
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_hi*nhi(j)/e_v *                       &
			     merge(fiHI *(e_v-e_th_HI)/e_th_HI , 0.0d0, e_v > e_th_HI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hi*nhi(j)/e_v *                       &
			     merge(fiHeI*(e_v-e_th_HI)/e_th_HeI, 0.0d0, e_v > e_th_HI + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_hi*nhi(j)/e_v *                &
			     merge(fiH2 *(e_v-e_th_HI)/e_th_H2 , 0.0d0, e_v > e_th_HI + E_sec_ion)

		fhv = 1.0d0
		if (sec_on) then
			call photoelectron_shares(pep, iabs_HeI, fh_v, fiHI, fiHeI, fiH2)
			fhv = merge(fh_v, 1.0d0, e_v > e_th_HeI + E_sec_ion)
		endif
		acc_HeI = photoelectron_share(e_th_HeI,e_v)*fhv*s_hei
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_hei*nheiS(j)/e_v *                    &
			     merge(fiHI *(e_v-e_th_HeI)/e_th_HI , 0.0d0, e_v > e_th_HeI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hei*nheiS(j)/e_v *                    &
			     merge(fiHeI*(e_v-e_th_HeI)/e_th_HeI, 0.0d0, e_v > e_th_HeI + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_hei*nheiS(j)/e_v *             &
			     merge(fiH2 *(e_v-e_th_HeI)/e_th_H2 , 0.0d0, e_v > e_th_HeI + E_sec_ion)

		fhv = 1.0d0
		if (sec_on) then
			call photoelectron_shares(pep, iabs_HeII, fh_v, fiHI, fiHeI, fiH2)
			fhv = merge(fh_v, 1.0d0, e_v > e_th_HeII + E_sec_ion)
		endif
		acc_HeII = photoelectron_share(e_th_HeII,e_v)*fhv*s_heii
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_heii*nheii(j)/e_v *                   &
			     merge(fiHI *(e_v-e_th_HeII)/e_th_HI , 0.0d0, e_v > e_th_HeII + E_sec_ion)
			acc_secHeI = acc_secHeI + s_heii*nheii(j)/e_v *                   &
			     merge(fiHeI*(e_v-e_th_HeII)/e_th_HeI, 0.0d0, e_v > e_th_HeII + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_heii*nheii(j)/e_v *            &
			     merge(fiH2 *(e_v-e_th_HeII)/e_th_H2 , 0.0d0, e_v > e_th_HeII + E_sec_ion)

		! He I 2^3S (triplet): photoelectron energy hv - 4.8 eV, same
		! secondary partition as the other absorbers.
		if (thereis_HeITR) then
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, iabs_HeTR, fh_v, fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > e_th_HeTR + E_sec_ion)
			endif
			acc_HeTR = photoelectron_share(e_th_HeTR,e_v)*fhv*s_heiTR
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_heiTR*nheiTR(j)/e_v *             &
				     merge(fiHI *(e_v-e_th_HeTR)/e_th_HI , 0.0d0, e_v > e_th_HeTR + E_sec_ion)
				acc_secHeI = acc_secHeI + s_heiTR*nheiTR(j)/e_v *             &
				     merge(fiHeI*(e_v-e_th_HeTR)/e_th_HeI, 0.0d0, e_v > e_th_HeTR + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_heiTR*nheiTR(j)/e_v *      &
				     merge(fiH2 *(e_v-e_th_HeTR)/e_th_H2 , 0.0d0, e_v > e_th_HeTR + E_sec_ion)
		endif

		! The molecular absorber, in its four final-state channels. Each
		! consumes one H2 and all four are inside the SAME cross section
		! (s_h2_di, s_h2_dd and s_h2_nd are shares of s_h2, not additions to
		! it), so the opacity, the absorbed energy and the H2 destruction
		! rate are what they were; what differs is the threshold each
		! channel charges -- 15.4 eV to leave H2+, 18.08 eV to leave H + H+,
		! 51.4 eV to leave two protons, and the bond energy alone where no
		! ion is made -- and therefore the energy left as heat.
		! s_h2 - s_h2_di - s_h2_dd - s_h2_nd is the non-dissociative channel
		! that leaves H2+; the last two vanish identically unless their
		! options are on, so the arithmetic below is unchanged by default.
		if (present(nh2)) then
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, iabs_H2, fh_v, fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > e_th_H2 + E_sec_ion)
			endif
			acc_H2 = photoelectron_share(e_th_H2,e_v)*fhv*(s_h2 - s_h2_di - s_h2_dd - s_h2_nd)
			if (sec_on) then
				acc_secHI  = acc_secHI  + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2(j)/e_v *       &
				     merge(fiHI *(e_v-e_th_H2)/e_th_HI , 0.0d0, e_v > e_th_H2 + E_sec_ion)
				acc_secHeI = acc_secHeI + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2(j)/e_v *       &
				     merge(fiHeI*(e_v-e_th_H2)/e_th_HeI, 0.0d0, e_v > e_th_H2 + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2(j)/e_v *&
				     merge(fiH2 *(e_v-e_th_H2)/e_th_H2 , 0.0d0, e_v > e_th_H2 + E_sec_ion)

			! Dissociative ionization H2 + hv -> H + H+ + e-. The 2.68 eV
			! between the two thresholds goes into breaking the bond and is
			! not available as heat, exactly as every other channel here is
			! charged its own ionization potential.
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, iabs_H2_di, fh_v, fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > e_th_H2_di + E_sec_ion)
			endif
			acc_H2 = acc_H2 + photoelectron_share(e_th_H2_di,e_v)*fhv*s_h2_di
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_h2_di*nh2(j)/e_v *                &
				     merge(fiHI *(e_v-e_th_H2_di)/e_th_HI , 0.0d0, e_v > e_th_H2_di + E_sec_ion)
				acc_secHeI = acc_secHeI + s_h2_di*nh2(j)/e_v *                &
				     merge(fiHeI*(e_v-e_th_H2_di)/e_th_HeI, 0.0d0, e_v > e_th_H2_di + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_h2_di*nh2(j)/e_v *         &
				     merge(fiH2 *(e_v-e_th_H2_di)/e_th_H2 , 0.0d0, e_v > e_th_H2_di + E_sec_ion)

			! Double ionization H2 + hv -> H+ + H+ + 2e-, charged its own
			! 51.4 eV threshold. The excess e_v - 51.4 eV is the TOTAL
			! kinetic energy of the two photoelectrons, and it is partitioned
			! between local heat and secondary ionization exactly as the two
			! channels above partition theirs. It must not be claimed whole
			! as heat: this channel lives only above 51.4 eV, which is the
			! band where the secondary-ionization share is largest, so
			! charging fhv = 1 here would overstate the local heating and
			! lose the ionizations the fast electrons make.
			!
			! APPROXIMATION, with its range. photoelectron_shares is indexed
			! by absorber and returns the partition of ONE electron carrying
			! the whole excess; this channel makes TWO electrons that share
			! it. The shares of the dissociative absorber are used at this
			! photon energy, which is the closest of the tabulated ones (same
			! molecular target, adjacent threshold). Because the partition
			! varies slowly with electron energy above E_sec_ion = 30 eV,
			! the error is second order in the difference between the mean
			! electron energy and the full excess; it is NOT zero, and it
			! biases toward too little heating and too much ionization,
			! since a slower electron heats more. A partition evaluated at
			! (e_v - 51.4)/2 would remove it and needs its own absorber row
			! in electron_energy_degradation.
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, iabs_H2_di, fh_v, fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > e_th_H2_dd + E_sec_ion)
			endif
			acc_H2 = acc_H2 + photoelectron_share(e_th_H2_dd,e_v)*fhv*s_h2_dd
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_h2_dd*nh2(j)/e_v *                &
				     merge(fiHI *(e_v-e_th_H2_dd)/e_th_HI , 0.0d0, e_v > e_th_H2_dd + E_sec_ion)
				acc_secHeI = acc_secHeI + s_h2_dd*nh2(j)/e_v *                &
				     merge(fiHeI*(e_v-e_th_H2_dd)/e_th_HeI, 0.0d0, e_v > e_th_H2_dd + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_h2_dd*nh2(j)/e_v *         &
				     merge(fiH2 *(e_v-e_th_H2_dd)/e_th_H2 , 0.0d0, e_v > e_th_H2_dd + E_sec_ion)

			! Neutral dissociation H2 + hv -> H + H. NO photoelectron: the
			! channel is not an ionization, so no fhv factor and no
			! secondary-ionization terms. The bond energy D0(H2) is spent
			! breaking the molecule and the whole remainder of the photon
			! goes into the kinetic energy of the two H atoms, i.e. into
			! heat. The channel is nonzero only over 33-41 eV, where
			! e_v >> D0(H2), so the bracket is positive throughout.
			acc_H2 = acc_H2 + (1.0-D0_H2_eV/e_v)*s_h2_nd
		endif

		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, n_abs_fixed+i, fh_v,       &
				                          fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > mion_ethr(i) + E_sec_ion)
			endif
			acc_mion   = photoelectron_share(mion_ethr(i),e_v)*fhv*sigma_tab(:,k)
			h1m_loc(i) = sum(int_f*acc_mion*de_v)*1.0e-18
			if (sec_on) then
				acc_secHI  = acc_secHI  + sigma_tab(:,k)*nm(j,i)/e_v *        &
				     merge(fiHI *(e_v-mion_ethr(i))/e_th_HI , 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
				acc_secHeI = acc_secHeI + sigma_tab(:,k)*nm(j,i)/e_v *        &
				     merge(fiHeI*(e_v-mion_ethr(i))/e_th_HeI, 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + sigma_tab(:,k)*nm(j,i)/e_v * &
				     merge(fiH2 *(e_v-mion_ethr(i))/e_th_H2 , 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
		enddo
		! Absorbed energy integral (this one is a property of the entry
		! composition: it is what the heating efficiency q is measured
		! against, so it keeps the absorber densities)
		acc_q = s_hi *nhi  (j)
		acc_q = acc_q + s_hei *nheiS(j)
		acc_q = acc_q + s_heii*nheii(j)
		if (thereis_HeITR) acc_q = acc_q + s_heiTR*nheiTR(j)
		if (present(nh2)) acc_q = acc_q + s_h2*nh2(j)
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			acc_q = acc_q + sigma_tab(:,k)*nm(j,i)
		enddo
		int_q = int_f*acc_q

		! Midpoint rule for integrals
		PIR_1   = sum(int_1  *de_v)
		PIR_15  = sum(int_15 *de_v)
		PIR_2   = sum(int_2  *de_v)
		if(thereis_HeITR) PIR_TR = sum(int_TR*de_v)
		PIR_H2 = 0.0
		PIR_H2_di = 0.0
		PIR_H2_dd = 0.0
		PIR_H2_nd = 0.0
		if (present(nh2)) then
			PIR_H2    = sum(int_h2   *de_v)
			PIR_H2_di = sum(int_h2_di*de_v)
			PIR_H2_dd = sum(int_h2_dd*de_v)
			PIR_H2_nd = sum(int_h2_nd*de_v)
		endif
		do i = 1,n_mion
			if (.not. mion_isphot(i)) then
				Pm_loc(i) = 0.0
				cycle
			endif
			k = mion_iphot(i)
			int_m = int_f*sigma_tab(:,k)/e_v
			Pm_loc(i) = sum(int_m*de_v)*1.0e-18*erg2eV
		enddo
		q_abs   = sum(int_q  *de_v)

		! Photoheating rate of one particle of each absorber [erg s^-1] in
		! this cell's attenuated field (the metal ions were integrated in the
		! loop above). The contraction with the composition follows the loop.
		h1_HI(j)   = sum(int_f*acc_HI  *de_v)*1.0e-18
		h1_HeI(j)  = sum(int_f*acc_HeI *de_v)*1.0e-18
		h1_HeII(j) = sum(int_f*acc_HeII*de_v)*1.0e-18
		h1_HeTR(j) = sum(int_f*acc_HeTR*de_v)*1.0e-18
		h1_H2(j)   = sum(int_f*acc_H2  *de_v)*1.0e-18
		h1_m(j,:)  = h1m_loc(:)
		q_abs_cell(j) = q_abs*1.0e-18

		! Save into vectors. No synchronization needed: each thread writes
		! only its own j elements of the shared arrays (the former OMP
		! CRITICAL serialized the loop for no correctness benefit; removed
		! per the 2026-07-02 review's performance note -- values unchanged).
    	P_HI(j)    = PIR_1  *1.0e-18*erg2eV
		if (present(P_H2)) P_H2(j) = PIR_H2*1.0e-18*erg2eV
		if (present(P_H2_di)) P_H2_di(j) = PIR_H2_di*1.0e-18*erg2eV
		if (present(P_H2_dd)) P_H2_dd(j) = PIR_H2_dd*1.0e-18*erg2eV
		if (present(P_H2_nd)) P_H2_nd(j) = PIR_H2_nd*1.0e-18*erg2eV
    	P_HeI(j)   = PIR_15 *1.0e-18*erg2eV
		! Add the H I / He I / H2 secondary-ionization rates from fast
		! photoelectrons. fiHI, fiHeI and fiH2 are per target particle, so
		! the integrals are already rates [1/s]: no division by a neutral
		! density, and each stays finite as its own target vanishes
		! (electron_energy_degradation.f90).
		if (sec_on) then
			Psec_HI  = sum(int_f*acc_secHI *de_v)*1.0e-18*erg2eV
			Psec_HeI = sum(int_f*acc_secHeI*de_v)*1.0e-18*erg2eV
			P_HI(j)  = P_HI(j)  + Psec_HI
			P_HeI(j) = P_HeI(j) + Psec_HeI
		endif
		! The molecular target. P_H2 drives BOTH the H2 destruction row and
		! the H2+ production row of the molecular system, exactly as the
		! primary photoionization rate it is added to does, so the secondary
		! ionizations enter the chemistry through the channel that already
		! exists rather than through a parallel one.
		if (mol_sec) then
			Psec_H2 = sum(int_f*acc_secH2*de_v)*1.0e-18*erg2eV
			! The secondaries dissociatively ionize H2 as well: one proton
			! (and one H atom) per 22 H2+ ions, Dalgarno, Yan & Liu (1999)
			! after their eq. (10). It is an ADDITIONAL yield, not a share of
			! Psec_H2 -- their harmonic-mean statement makes the two H+
			! sources add -- so the H2 destroyed by secondaries is the sum of
			! the two and only the second of them makes protons.
			Psec_H2_di = dissoc_ion_per_H2p*Psec_H2
			P_H2(j) = P_H2(j) + Psec_H2 + Psec_H2_di
			if (present(P_H2_di)) P_H2_di(j) = P_H2_di(j) + Psec_H2_di
		endif
		P_HeII(j)  = PIR_2  *1.0e-18*erg2eV
		P_HeITR(j) = PIR_TR *1.0e-18*erg2eV
		P_m(j,:)   = Pm_loc(:)

	enddo
	!$OMP END PARALLEL DO

	! Heating of the composition this routine was given: the one-particle
	! rates above contracted with the densities, through the single
	! definition of that sum. `heat` and the six-column split are therefore
	! the same object, and the split is exact.
	nh2_loc = 0.0d0
	if (present(nh2)) nh2_loc = nh2
	call photoheating_of_composition(h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2,  &
	                                 h1_m, nhi,nhei,nheii,nheiTR,nh2_loc, &
	                                 nm, heat, heat_chan)

	! Heating efficiency: heat deposited over energy absorbed, both of the
	! composition passed in. Absorbed energy is non-negative; see PH_heat_H
	! for why a non-positive q_abs gives a zero heating efficiency.
	where (q_abs_cell .gt. 0.0d0)
		q = heat/q_abs_cell
	elsewhere
		q = 0.0d0
	endwhere

	if (present(heat_of_one_HI))   heat_of_one_HI   = h1_HI
	if (present(heat_of_one_HeI))  heat_of_one_HeI  = h1_HeI
	if (present(heat_of_one_HeII)) heat_of_one_HeII = h1_HeII
	if (present(heat_of_one_HeTR)) heat_of_one_HeTR = h1_HeTR
	if (present(heat_of_one_H2))   heat_of_one_H2   = h1_H2
	if (present(heat_of_one_mion)) heat_of_one_mion = h1_m

	! End of subroutine
	end subroutine PH_heat_HHe

	! ------------------------------------------------------------- !

	subroutine photoheating_of_composition(h1_HI,h1_HeI,h1_HeII,h1_HeTR,  &
	                                       h1_H2, h1_m,                   &
	                                       nhi,nhei,nheii,nheiTR,nh2,nm,  &
	                                       heat, heat_chan)
	! Photoheating rate of a composition [erg cm^-3 s^-1]: the heating rate
	! of ONE particle of each absorber (h1_*, from PH_heat_HHe / PH_heat_H,
	! a property of the radiation field alone) contracted with the density
	! of that absorber.
	!
	! This is the ONE place the contraction is written. The rates are lagged
	! quantities -- they are built from the columns and the photoelectron
	! partition of the state that entered the ionization sweep -- while the
	! densities are whichever composition the caller wants the heating of;
	! ionization_equilibrium calls this after its sweep so that the heating
	! it returns belongs to the composition it returns.
	!
	! He I absorbs in its GROUND SINGLET, so the He I density is reduced by
	! the 2^3S metastable through the single definition of that difference
	! (he_ground_singlet_density); the metastable has its own column.

	real*8, dimension(1-Ng:N+Ng),intent(in) :: h1_HI,h1_HeI,h1_HeII,     &
	                                            h1_HeTR,h1_H2
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: h1_m
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhei,nheii,nheiTR,nh2
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm

	real*8, dimension(1-Ng:N+Ng),intent(out) :: heat
	! Optional split by absorber: 1 H I, 2 He I (ground singlet), 3 He II,
	! 4 He 2^3S, 5 H2, 6 metals (summed over the ions). The six columns sum
	! to `heat` exactly.
	real*8, dimension(1-Ng:N+Ng,6),intent(out),optional :: heat_chan

	real*8, dimension(1-Ng:N+Ng) :: nheiS, h_HI,h_HeI,h_HeII,h_HeTR,     &
	                                 h_H2,h_mtl
	integer :: i

	nheiS = nhei
	if (thereis_HeITR) nheiS = he_ground_singlet_density(nhei, nheiTR)

	h_HI   = h1_HI  *nhi
	h_HeI  = h1_HeI *nheiS
	h_HeII = h1_HeII*nheii
	h_HeTR = h1_HeTR*nheiTR
	h_H2   = h1_H2  *nh2
	h_mtl  = 0.0d0
	do i = 1,n_mion
		h_mtl = h_mtl + h1_m(:,i)*nm(:,i)
	enddo

	heat = h_HI + h_HeI + h_HeII + h_HeTR + h_H2 + h_mtl

	if (present(heat_chan)) then
		heat_chan(:,1) = h_HI
		heat_chan(:,2) = h_HeI
		heat_chan(:,3) = h_HeII
		heat_chan(:,4) = h_HeTR
		heat_chan(:,5) = h_H2
		heat_chan(:,6) = h_mtl
	endif

	end subroutine photoheating_of_composition
	
	!----------------------------------!
	
	subroutine eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,        &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,             &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,      &
				   cool, cool_chan, nheiTR, a_ion_HeITR, nmol, nox)

	! Evaluate the cooling rate contributions to energy and
	!	rate equations. Includes H, He, and metal channels. Metals are
	!	driven from species_table metadata (canonical ion order); nm,
	!	rec_m and aion_m carry the densities and rates for each ion, so the
	!	argument list no longer grows when a metal is added.

	! Structure. Everything below is a function of the cell's own temperature
	! and densities, with ONE exception: the fine-structure line transfer,
	! whose escape probabilities are integrals of the columns above and below
	! the cell. That exception is solved here for the whole grid; the
	! cell-local remainder runs in eval_cool_cells over contiguous blocks of
	! cells inside a single OpenMP parallel region. Each cell therefore gets
	! the same arithmetic as a serial sweep, and the result is bitwise
	! reproducible at any thread count.

	! Block decomposition of the cell range and the sub-block timers a block
	! returns (see ec_t / ec_name).
	integer :: ib, nblk, j_lo, j_hi, ncell, nthr
	real*8  :: ect(5)

	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,           &
	                                            nhei,nheii,nheiii
	! Metal ion densities (canonical species_table order)
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm

	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

   ! Line transfer of the ground-term fine-structure lines solved explicitly
   ! (Cool_coeff: fine_structure_line_transfer): escape probabilities and the
   ! photon occupation number of the field incident from the lower atmosphere
   real*8, dimension(1-Ng:N+Ng,n_fsline) :: beta_fs, nbar_fs
	real*8, dimension(1-Ng:N+Ng) :: ne		  			 ! Electron number density

   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB
   ! Metal recombination rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: rec_m

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII
   ! Metal collisional ionization rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: aion_m

	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  cool

	! Optional cooling breakdown in each channel (cgs erg cm^-3 s^-1, same
	! units as `cool`). Columns 1-6 = H/He recombination, collisional
	! ionization, then the collisional-excitation channel split into its
	! three absorbers -- H I (the Lyman-alpha-dominated H-line cooling),
	! He I, He II -- and finally bremsstrahlung (incl. metal-ion charges);
	! column 7 = H3+ infrared cooling (0 unless the caller supplies nmol);
	! columns 8-10 = the molecular infrared bands under `Molecular IR bands`
	! -- H2 lines (needs nmol), H2O and CO bands (need nox) -- each the NET
	! rate, emission minus absorption of the field from below, so a column is
	! negative wherever that channel heats;
	! columns 10+i = metal ion i line cooling (0 for non-coolant ions).
	! The He I column also carries the He 2^3S metastable collisional cooling
	! (10830 A + singlet-conversion terms), so columns 3-5 sum exactly to
	! ne*coex. This is an exact decomposition of `cool` in the default
	! (.not.use_2lev_cool) branch; in the two-level branch the metal terms for
	! each ion are the resonance-line approximation and need not sum to cool_M.
	real*8, dimension(1-Ng:N+Ng,10+n_mion),intent(out),optional :: cool_chan

	! He 2^3S metastable density [cm^-3], present only for the triplet-tracking
	! callers. When supplied it adds the collisional-ionization cooling of the
	! 2^3S state (4.8 eV per event, ci_HeI23S) to the CI channel. a_ion_HeITR
	! returns the He(2^3S) collisional-ionization rate coefficient [cm^3 s^-1]
	! for the ionization equations, mirroring a_ion_HI/HeI/HeII.
	real*8, dimension(1-Ng:N+Ng),intent(in),optional  :: nheiTR
	real*8, dimension(1-Ng:N+Ng),intent(out),optional :: a_ion_HeITR

	! Molecular densities [cm^-3], canonical order H2, H2+, H3+, HeH+ (the
	! same layout calc_ne takes). Supplied by every caller that tracks the
	! molecular network, so the electron density used by the cooling is the
	! one the equilibrium solver itself uses, and so the H3+ infrared cooling
	! (which needs n_H3+ and the n_H2 collider density) is part of the same
	! `cool` every caller gets. Omitted only by callers that model a
	! molecule-free gas (see the note at calc_ne below).
	real*8, dimension(1-Ng:N+Ng,4),intent(in),optional :: nmol

	! Oxygen-carrier densities [cm^-3], canonical order OH, H2O, CO (the
	! layout nox_eq uses). Present only when the oxygen chemistry is on, which
	! is the only configuration in which H2O and CO exist at all. Used for the
	! H2O and CO infrared bands; OH has no cross-section table and is left out
	! (it carries a few percent of the oxygen where the water does the rest).
	real*8, dimension(1-Ng:N+Ng,3),intent(in),optional :: nox

   ! Free electron density: ONE definition, the calc_ne charge sum over every
   ! tracked ion -- H+, He+, He++, the metal stages under eos_metals, and the
   ! molecular ions H2+, H3+, HeH+ when the caller tracks them. Inside a deep
   ! molecular base H3+ is the dominant ion, so dropping the molecular donors
   ! would under-count every ne-scaling cooling channel (recombination,
   ! collisional excitation/ionization, bremsstrahlung, metal line cooling)
   ! there, and would disagree with the ne the equilibrium solver balances
   ! ionization against. Callers that model a molecule-free gas (the _adv
   ! advection post-process) omit nmol and get the atomic charge sum.
   if (present(nmol)) then
      call calc_ne(nhii,nheii,nheiii,ne,nm,nmol)
   else
      call calc_ne(nhii,nheii,nheiii,ne,nm)
   endif

	! Line trapping of the ground-term fine-structure lines this module
	! solves explicitly. Earlier versions built a GRAY escape
	! probability from the lowest XUV band opacity over one cell width
	! (AIOLOS chemistry.cpp:1006 scaled that by an arbitrary 1e8, driving
	! beta -> 0 and switching metal-line cooling off; EXHALE replaced it
	! with beta = 1 everywhere, the optically thin limit). Neither is a
	! line optical depth, and beta = 1 overestimates the cooling of a
	! dense base where [O I] 63um reaches tau ~ 3. beta is now the
	! line-center escape probability of each line, from the columns above
	! and below the cell, and nbar_fs carries the thermal infrared field of
	! the lower atmosphere the same lines absorb ("Base IR field", off by
	! default) (Cool_coeff.f90: fine_structure_line_transfer). Every other
	! metal ion keeps the optically thin limit; see the scope note there.
	beta_fs = 1.0d0
	nbar_fs = 0.0d0
	!$ if (ec_prof_on) ec_t0 = omp_get_wtime()
	if (thereis_metals) call fine_structure_line_transfer(T_K, nm,      &
	                                                  beta_fs, nbar_fs)
	!$ if (ec_prof_on) ec_t(4) = ec_t(4) + (omp_get_wtime() - ec_t0)

	! One parallel region per call, over contiguous blocks of cells. The
	! blocks are disjoint and every quantity is cell-local, so the thread
	! count changes neither the arithmetic nor the result. The sub-block
	! timers of the blocks are summed, so slots 1-3, 5 and the CNO part of
	! slot 4 report thread time rather than wall time once threads are on.
	ncell = N + 2*Ng
	nthr  = 1
	!$ nthr = omp_get_max_threads()
	nblk  = max(1, min(2*nthr, ncell))
	!$omp parallel do schedule(dynamic) default(shared)                 &
	!$omp    private(ib,j_lo,j_hi,ect) reduction(+:ec_t)
	do ib = 1,nblk
	   j_lo = (1-Ng) + ((ib-1)*ncell)/nblk
	   j_hi = (1-Ng) + (ib*ncell)/nblk - 1
	   if (j_hi .lt. j_lo) cycle
	   ect = 0.0d0
	   call eval_cool_cells(j_lo,j_hi, ect,                             &
	          T_K,nhi,nhii,nhei,nheii,nheiii, nm, ne, beta_fs, nbar_fs, &
	          rchiiB,rcheiiB,rcheiiiB, rec_m,                           &
	          a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,                    &
	          cool, cool_chan, nheiTR, a_ion_HeITR, nmol, nox)
	   ec_t = ec_t + ect
	enddo
	!$omp end parallel do

	! End of subroutine
	end subroutine eval_cool

	!---------------------------------------------------!

	subroutine eval_cool_cells(j_lo,j_hi, ect,                       &
	                   T_K,nhi,nhii,nhei,nheii,nheiii, nm, ne,       &
	                   beta_fs, nbar_fs,                             &
	                   rchiiB,rcheiiB,rcheiiiB, rec_m,               &
	                   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,        &
	                   cool, cool_chan, nheiTR, a_ion_HeITR, nmol, nox)

	! The cell-local part of eval_cool, over the cells j_lo:j_hi. Called
	! from inside the parallel region of eval_cool, once per block of
	! cells. Every array keeps the explicit shape 1-Ng:N+Ng of the whole
	! grid and only the elements j_lo:j_hi are read or written, so the
	! compiled arithmetic of a cell is what the serial whole-array form
	! produced (assumed-shape dummies moved it; Update_EXHALE.md 145).
	! ne, beta_fs and nbar_fs come in from eval_cool: they are not
	! cell-local. ect returns this block's sub-block times.

	integer, intent(in) :: j_lo, j_hi
	real*8, intent(inout) :: ect(5)
	real*8 :: ec_tstart

	integer :: i,j

	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,           &
	                                            nhei,nheii,nheiii
	! Metal ion densities (canonical species_table order)
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm

	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

   real*8, dimension(1-Ng:N+Ng) :: brem,coex,coio,reco  ! Cooling rates
   real*8, dimension(1-Ng:N+Ng) :: cool_M               ! Metal cooling
   real*8, dimension(1-Ng:N+Ng) :: cool_H3p             ! H3+ infrared cooling
   real*8, dimension(1-Ng:N+Ng) :: cool_H2, cool_H2O, cool_CO ! molecular bands
   real*8 :: w_ir                                       ! incident-field dilution
   real*8, dimension(1-Ng:N+Ng) :: brem_acc,coolm_acc   ! sum accumulators
   real*8, dimension(1-Ng:N+Ng) :: metal_col            ! dispatcher scratch
   real*8, dimension(1-Ng:N+Ng,n_mion)  :: c_metal      ! metal line-cool coeffs
   ! Line transfer of the ground-term fine-structure lines solved explicitly
   ! (Cool_coeff: fine_structure_line_transfer): escape probabilities and the
   ! photon occupation number of the field incident from the lower atmosphere
	real*8, dimension(1-Ng:N+Ng,n_fsline),intent(in) :: beta_fs, nbar_fs
	real*8, dimension(1-Ng:N+Ng),intent(in) :: ne        ! Electron density
	real*8, dimension(1-Ng:N+Ng) :: GF_z1,GF_z2			 ! free-free Gaunt at Z_ion=1,2
	real*8 :: Cdex_OI,Cdex_CII                           ! 2-level collis. de-exc.
	! Named bridges for the (verbatim) two-level cooling branch
	real*8, dimension(1-Ng:N+Ng) :: nci,ncii,noi,noii,nmgi,nmgii
	real*8, dimension(1-Ng:N+Ng) :: c_CI,c_OII,c_MgI,c_MgII

   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(inout) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB
   ! Metal recombination rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(inout) :: rec_m

	! Recombination cooling coefficients
	real*8, dimension(1-Ng:N+Ng) :: coeff_rec_cool_HII,  &
									coeff_rec_cool_HeII, &
									coeff_rec_cool_HeIII

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(inout) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII
   ! Metal collisional ionization rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(inout) :: aion_m

	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_rate_HI,    &
	 								coeff_coex_rate_HeI,   &
									coeff_coex_rate_HeII
	! He 2^3S metastable cooling coefficients (triplet-tracking callers only):
	! 10830 A collisional-excitation cooling, and the q31a/q31b conversion
	! rate coefficients reused for their thermal-energy ledger.
	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_HeI23S_10830, q31a_l, q31b_l

	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(inout) ::  cool

	! Optional cooling breakdown in each channel (cgs erg cm^-3 s^-1, same
	! units as `cool`). Columns 1-6 = H/He recombination, collisional
	! ionization, then the collisional-excitation channel split into its
	! three absorbers -- H I (the Lyman-alpha-dominated H-line cooling),
	! He I, He II -- and finally bremsstrahlung (incl. metal-ion charges);
	! column 7 = H3+ infrared cooling (0 unless the caller supplies nmol);
	! columns 8-10 = the molecular infrared bands under `Molecular IR bands`
	! -- H2 lines (needs nmol), H2O and CO bands (need nox) -- each the NET
	! rate, emission minus absorption of the field from below, so a column is
	! negative wherever that channel heats;
	! columns 10+i = metal ion i line cooling (0 for non-coolant ions).
	! The He I column also carries the He 2^3S metastable collisional cooling
	! (10830 A + singlet-conversion terms), so columns 3-5 sum exactly to
	! ne*coex. This is an exact decomposition of `cool` in the default
	! (.not.use_2lev_cool) branch; in the two-level branch the metal terms for
	! each ion are the resonance-line approximation and need not sum to cool_M.
	real*8, dimension(1-Ng:N+Ng,10+n_mion),intent(inout),optional :: cool_chan

	! He 2^3S metastable density [cm^-3], present only for the triplet-tracking
	! callers. When supplied it adds the collisional-ionization cooling of the
	! 2^3S state (4.8 eV per event, ci_HeI23S) to the CI channel. a_ion_HeITR
	! returns the He(2^3S) collisional-ionization rate coefficient [cm^3 s^-1]
	! for the ionization equations, mirroring a_ion_HI/HeI/HeII.
	real*8, dimension(1-Ng:N+Ng),intent(in),optional  :: nheiTR
	real*8, dimension(1-Ng:N+Ng),intent(inout),optional :: a_ion_HeITR

	! Molecular densities [cm^-3], canonical order H2, H2+, H3+, HeH+ (the
	! same layout calc_ne takes). Supplied by every caller that tracks the
	! molecular network, so the electron density used by the cooling is the
	! one the equilibrium solver itself uses, and so the H3+ infrared cooling
	! (which needs n_H3+ and the n_H2 collider density) is part of the same
	! `cool` every caller gets. Omitted only by callers that model a
	! molecule-free gas (see the note at calc_ne below).
	real*8, dimension(1-Ng:N+Ng,4),intent(in),optional :: nmol

	! Oxygen-carrier densities [cm^-3], canonical order OH, H2O, CO (the
	! layout nox_eq uses). Present only when the oxygen chemistry is on, which
	! is the only configuration in which H2O and CO exist at all. Used for the
	! H2O and CO infrared bands; OH has no cross-section table and is left out
	! (it carries a few percent of the oxygen where the water does the rest).
	real*8, dimension(1-Ng:N+Ng,3),intent(in),optional :: nox

	! He 2^3S collisional-ionization rate coefficient (always computed; only
	! exported / applied through the optional arguments above).
	real*8, dimension(1-Ng:N+Ng) :: aion_HeITR

	!-- Recombination --!

	! Rate coefficients
	!$ if (ec_prof_on) ec_tstart = omp_get_wtime()
	call rec_HII_B_range(T_K,rchiiB,j_lo,j_hi)      ! HII
	call rec_HeII_B_range(T_K,rcheiiB,j_lo,j_hi)    ! HeII
	call rec_HeIII_B_range(T_K,rcheiiiB,j_lo,j_hi)  ! HeIII
	
	! Cooling rate coefficients
	call rec_cool_HII_range(T_K,coeff_rec_cool_HII,j_lo,j_hi)
	call rec_cool_HeII_range(T_K,coeff_rec_cool_HeII,j_lo,j_hi)
	call rec_cool_HeIII_range(T_K,coeff_rec_cool_HeIII,j_lo,j_hi)

	! Cooling rate
	reco(j_lo:j_hi)  = coeff_rec_cool_HII(j_lo:j_hi)*nhii(j_lo:j_hi)     & ! HII
		   + coeff_rec_cool_HeII(j_lo:j_hi)*nheii(j_lo:j_hi)   & ! HeII
		   + coeff_rec_cool_HeIII(j_lo:j_hi)*nheiii(j_lo:j_hi)   ! HeIII

	!-- Collisional ionization --!
	
	! Rate coefficients
	call ion_coeff_HI_range(T_K,a_ion_HI,j_lo,j_hi)      ! HI
	call ion_coeff_HeI_range(T_K,a_ion_HeI,j_lo,j_hi)    ! HeI
	call ion_coeff_HeII_range(T_K,a_ion_HeII,j_lo,j_hi)  ! HeII
	call ci_HeI23S_range(T_K,aion_HeITR,j_lo,j_hi)       ! He 2^3S metastable (4.8 eV threshold)
	if (present(a_ion_HeITR)) a_ion_HeITR(j_lo:j_hi) = aion_HeITR(j_lo:j_hi)

	! Cooling rate. The prefactors are the ionization potentials in erg
	! (13.6/24.6/54.4 eV for HI/HeI/HeII; e_th_HeTR = 4.8 eV for the 2^3S
	! metastable, e_th_HeTR/erg2eV = 7.69e-12 erg). He II reads its potential
	! from the named e_th_HeII_erg, so this assembly and the cell-by-cell
	! temperature root of the post-process charge the same energy per event;
	! the H I and He I prefactors are still literals here, and they agree
	! between the two sites. The 2^3S term (added only
	! when nheiTR is supplied) reproduces the Black (1981) form
	! 6.41e-21 sqrt(T) exp(-55338/T) n_e n_23S once multiplied by n_e below.
	coio(j_lo:j_hi) =  2.179e-11*a_ion_HI(j_lo:j_hi)*nhi(j_lo:j_hi)  	 & ! HI
		  + 3.940e-11*a_ion_HeI(j_lo:j_hi)*nhei(j_lo:j_hi) 	 & ! HeI
		  + e_th_HeII_erg*a_ion_HeII(j_lo:j_hi)*nheii(j_lo:j_hi)   ! HeII
	if (present(nheiTR)) coio(j_lo:j_hi) = coio(j_lo:j_hi)                                     &
		  + (e_th_HeTR/erg2eV)*aion_HeITR(j_lo:j_hi)*nheiTR(j_lo:j_hi)   ! He 2^3S
	
	!-- Bremsstrahlung --!

	! Free-free scales with the ion NET charge Z_ion (not the nuclear number):
	! H II, He II and singly-ionized metals are Z_ion = 1; He III and doubly-
	! ionized metals are Z_ion = 2. The Gaunt factor is evaluated at Z_ion, so
	! only the charge-1 and charge-2 values are needed.
	!$ if (ec_prof_on) then
	!$    ect(1) = ect(1) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif
	call GF_range(T_K, 1.0d0, GF_z1,j_lo,j_hi)
	call GF_range(T_K, 2.0d0, GF_z2,j_lo,j_hi)

	! Cooling rate: sum n_ion * Z_ion^2 * gbar(Z_ion) over all charged ions.
	! mion_z2 = mion_stage^2 already holds the metal charge^2 (neutral -> 0,
	! an exact +0 term); the Gaunt table is selected by mion_stage.
	! The molecular ions H2+, H3+ and HeH+ are Z_ion = 1 and DO donate to ne
	! above, but are deliberately left out of this charge sum: they exist only
	! in the cold (T ~ 1e3 K) molecular base, where this hot-plasma free-free
	! expression is an extrapolation whose emission comes out at radio/IR
	! frequencies the atmosphere is not thin to. Free-free is negligible against
	! the H3+ infrared and metal line cooling there, so the ion charge sum is
	! smaller than ne by the molecular-ion density inside the molecular layer.
	brem_acc(j_lo:j_hi) = GF_z1(j_lo:j_hi)*nhii(j_lo:j_hi)                       ! HII   (Z_ion = 1)
	brem_acc(j_lo:j_hi) = brem_acc(j_lo:j_hi) + GF_z1(j_lo:j_hi)*nheii(j_lo:j_hi)           ! HeII  (Z_ion = 1)
	brem_acc(j_lo:j_hi) = brem_acc(j_lo:j_hi) + 4.0*GF_z2(j_lo:j_hi)*nheiii(j_lo:j_hi)      ! HeIII (Z_ion = 2)
	! Metal ions only when the run carries metals: with none, every term
	! below is an exact +0 and the loop is skipped (see the note before the
	! metal rate block).
	if (thereis_metals) then
	do i = 1,n_mion
		if (mion_stage(i) == 2) then
			brem_acc(j_lo:j_hi) = brem_acc(j_lo:j_hi) + mion_z2(i)*GF_z2(j_lo:j_hi)*nm(j_lo:j_hi,i)
		else
			brem_acc(j_lo:j_hi) = brem_acc(j_lo:j_hi) + mion_z2(i)*GF_z1(j_lo:j_hi)*nm(j_lo:j_hi,i)
		endif
	enddo
	endif
	brem(j_lo:j_hi) = 1.426e-27*sqrt(T_K(j_lo:j_hi))*brem_acc(j_lo:j_hi)

	!-- Collisional excitation --!

	! Rate coefficients
	call coex_rate_HI_range(T_K,coeff_coex_rate_HI,j_lo,j_hi) 		! HI
	call coex_rate_HeI_range(T_K,coeff_coex_rate_HeI,j_lo,j_hi)   	! HeI
	call coex_rate_HeII_range(T_K,coeff_coex_rate_HeII,j_lo,j_hi)  	! HeII

	! Cooling rate
	coex(j_lo:j_hi) = coeff_coex_rate_HI(j_lo:j_hi)*nhi(j_lo:j_hi)       &    ! HI
		  + coeff_coex_rate_HeI(j_lo:j_hi)*nhei(j_lo:j_hi)     &    ! HeI
		  + coeff_coex_rate_HeII(j_lo:j_hi)*nheii(j_lo:j_hi)        ! HeII

	! He 2^3S metastable collisional cooling (triplet-tracking callers only).
	! Both terms scale with the EXPLICITLY computed n_23S (nheiTR); the common
	! n_e factor is applied together with the other coex terms in `cool` below.
	!  (1) 2^3S -> 2^3P collisional excitation, then 10830 A radiative decay
	!      and photon escape (Black 1981 / Allan 2024 / Falorca & Vidotto 2026
	!      Table A2). The 10830 photon leaves the gas radially (the wind is thin
	!      to it); its line optical depth is what the transit sees. Black's
	!      implicit steady-state triplet form (~T^-0.6687 n_e n_He+) is NOT used.
	!  (2) 2^3S -> 2^1S (0.80 eV) and 2^3S -> 2^1P (1.40 eV) collisional
	!      conversions each remove their threshold energy from the electron gas;
	!      the excited singlet then decays radiatively (that photon is not
	!      thermal). Rate coefficients reused from coex_HeI_23S_21S / _21P.
	! The ground -> triplet excitation (q13, 19.82 eV) is deliberately EXCLUDED:
	! that channel is already carried by the Cen-1992 He I coex term above
	! (coeff_coex_rate_HeI*nhei); adding q13 here would double count it.
	if (present(nheiTR)) then
		call coex_rate_HeI23S_10830_range(T_K,coeff_coex_HeI23S_10830,j_lo,j_hi)
		call coex_HeI_23S_21S_range(T_K,q31a_l,j_lo,j_hi)
		call coex_HeI_23S_21P_range(T_K,q31b_l,j_lo,j_hi)
		coex(j_lo:j_hi) = coex(j_lo:j_hi) + ( coeff_coex_HeI23S_10830(j_lo:j_hi)                          &
		              + (0.80d0*q31a_l(j_lo:j_hi) + 1.40d0*q31b_l(j_lo:j_hi))/erg2eV )*nheiTR(j_lo:j_hi)
	endif

	!-- Metal recombination + collisional ionization rates --!
	! (rates for ionization equilibrium; not part of cool here.)
	! Filled per canonical ion via the metadata dispatchers, which call
	! the same routines for each ion as before; ions with no entry (inert top
	! stage, neutral non-recombiner) return 0.
	!$ if (ec_prof_on) then
	!$    ect(2) = ect(2) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif
	! THE METAL BLOCKS RUN ONLY WHEN THE RUN CARRIES METALS.  Measured on the
	! H/He hot Uranus (docs/marching_step_performance_plan.md): with every
	! metal density identically zero, the metal table interpolations, the
	! fine-structure transfer and the CNO cooling below were 79 per cent of
	! this routine's time -- and this routine, called four times a step, was
	! half of a 16-thread marching step.  Skipping them is bitwise identical
	! for a metals-off run: each skipped term is a product with a zero
	! density, and the outputs a caller could read are set here to what they
	! would have multiplied into (zero rates, beta = 1, no infrared field).
	! The ionization sweep already gates its metal solve on the same flag.
	rec_m(j_lo:j_hi,:)   = 0.0d0
	aion_m(j_lo:j_hi,:)  = 0.0d0
	c_metal(j_lo:j_hi,:) = 0.0d0
	cool_M(j_lo:j_hi)  = 0.0d0
	if (thereis_metals) then
	do i = 1,n_mion
		call rec_coeff_by_ion_range(i,T_K,metal_col,j_lo,j_hi)
		rec_m(j_lo:j_hi,i)  = metal_col(j_lo:j_hi)
		call ion_coeff_by_ion_range(i,T_K,metal_col,j_lo,j_hi)
		aion_m(j_lo:j_hi,i) = metal_col(j_lo:j_hi)
	enddo

	!-- Metal radiative cooling (forbidden/fine-structure lines) --!
	do i = 1,n_mion
		call cool_coeff_by_ion_range(i,T_K,metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,i) = metal_col(j_lo:j_hi)
	enddo

	! Density-dependent override for Fe II line cooling. The 1-D coronal
	! cool_coeff_metal('FeII',...) above overestimates cooling at the dense
	! base by ~1e4x because the forbidden a6D fine-structure / metastable
	! lines (n_crit ~ 1e4-1e7 cm^-3) are collisionally saturated there
	! (n_e >> n_crit). Replace c_metal(:,FeII) with the multilevel
	! statistical-equilibrium coefficient Lambda_eff(T,ne) = (sum_u n_u A_ul
	! dE_ul)/ne. In the assembly cool_M = ne*sum_i nm(:,i)*c_metal,
	! the ne cancels the 1/ne in Lambda_eff, leaving the correct LTE-saturated
	! cooling for each ion (collider-independent, so the electron-only SE solve is
	! exact in this limit). At low ne it reduces to the coronal rate.
	call cool_FeII_ne_range(T_K, ne, metal_col,j_lo,j_hi)
	c_metal(j_lo:j_hi,im_FeII) = metal_col(j_lo:j_hi)

	!$ if (ec_prof_on) then
	!$    ect(3) = ect(3) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif

	! Density-dependent override for the ground-term fine-structure floors
	! of C I, C II, N II and O I (CHIANTI mode only; the legacy AIOLOS fits
	! keep their own constant floors). Same Lambda_eff = W_FS/ne +
	! remainder convention as Fe II above; the statistical-equilibrium
	! solution saturates the floor (n_crit,e([C II] 158um) ~ 20 cm^-3!) and
	! adds the H-collision excitation channel the electron-only coronal
	! curve misses. beta enters as A_ul -> beta*A_ul inside that solution.
	! See cool_CI_ne_func / cooling_data/fit_fs_saturation.py.
	if (cno_chianti) then
		call cool_CI_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_CI)  = metal_col(j_lo:j_hi)
		call cool_CII_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_CII) = metal_col(j_lo:j_hi)
		call cool_NII_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_NII) = metal_col(j_lo:j_hi)
		call cool_OI_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_OI)  = metal_col(j_lo:j_hi)
	endif

	if (use_2lev_cool) then
		! Bridge the metadata arrays to the named scalars used below.
		nci(j_lo:j_hi)   = nm(j_lo:j_hi,1)
		ncii(j_lo:j_hi)  = nm(j_lo:j_hi,2)
		noi(j_lo:j_hi)   = nm(j_lo:j_hi,4)
		noii(j_lo:j_hi)  = nm(j_lo:j_hi,5)
		nmgi(j_lo:j_hi)  = nm(j_lo:j_hi,10)
		nmgii(j_lo:j_hi) = nm(j_lo:j_hi,11)
		c_CI(j_lo:j_hi)   = c_metal(j_lo:j_hi,1)
		c_OII(j_lo:j_hi)  = c_metal(j_lo:j_hi,5)
		c_MgI(j_lo:j_hi)  = c_metal(j_lo:j_hi,10)
		c_MgII(j_lo:j_hi) = c_metal(j_lo:j_hi,11)
		! Two-level fine-structure cooling for the dominant coolants
		! [O I] 63um and [C II] 158um (critical-density saturation +
		! H-atom collisions). C I and O II keep the Black-1981 fit; the
		! exponential ("forbidden") part of the C II / O I fit is retained
		! on top of the two-level ground term (matching ATES_extended).
		! N has no line cooling.
		! Trapping enters as A_ul -> beta*A_ul in the two fine-structure
		! lambda_2level calls (the physically correct place; see
		! cool_OI_ne_func). The exponential "forbidden" add-ons and every
		! other ion stay optically thin, and this branch takes no incident
		! field even under "Base IR field": it is the legacy AIOLOS two-level
		! form kept for comparison, and the ground-term statistical
		! equilibrium that replaced it is where the field belongs.
		do j = j_lo,j_hi
			! Collisional de-excitation rates [s^-1]
			Cdex_OI  = nhi(j)*4.2d-11*(T_K(j)/100.0d0)**0.67          ! H
			Cdex_CII = ne(j) *8.7d-8 *(T_K(j)/2000.0d0)**(-0.37)      & ! e
			         + nhi(j)*4.0d-11                                   ! H
			cool_M(j) =                                                    &
			    ne(j)*nci(j)*c_CI(j)                                       &
			  + noi(j) *( lambda_2level(beta_fs(j,ifs_OI63)*8.91d-5,227.7d0,0.6d0,Cdex_OI ,T_K(j)) &
			              + ne(j)*1.1d-20*exp(-30162.0d0/T_K(j))           &
			                     *(1.0d0+(T_K(j)/0.75d4)**0.5) )           &
			  + ncii(j)*( lambda_2level(beta_fs(j,ifs_CII158)*2.29d-6,91.21d0,2.0d0,Cdex_CII,T_K(j)) &
			              + ne(j)*3.1d-20*exp(-45162.0d0/T_K(j))           &
			                     *(1.0d0+(T_K(j)/0.75d4)**1.5) )           &
			  + ne(j)*noii(j)*c_OII(j)                                  &
			  + ne(j)*nmgi(j)*c_MgI(j) + ne(j)*nmgii(j)*c_MgII(j)                          &
                  + ne(j)*nm(j,17)*c_metal(j,17)               &
                  + ne(j)*nm(j,19)*c_metal(j,19)               &
                  + ne(j)*nm(j,26)*c_metal(j,26)
		enddo
	else
		! Sum the line-cooling metal ions (mion_iscool) in canonical
		! order, then apply the ne prefactor once. This is
		! bit-identical to the explicit eight-term expression: the
		! iscool ions, in canonical order, are exactly
		! CI,CII,OI,OII,NI,NII,MgI,MgII. Line trapping is already inside
		! c_metal for [O I] 63um / [C II] 158um; the rest are thin.
		coolm_acc(j_lo:j_hi) = 0.0d0
		do i = 1,n_mion
			if (.not. mion_iscool(i)) cycle
			coolm_acc(j_lo:j_hi) = coolm_acc(j_lo:j_hi) + nm(j_lo:j_hi,i)*c_metal(j_lo:j_hi,i)
		enddo
		cool_M(j_lo:j_hi) = ne(j_lo:j_hi) * coolm_acc(j_lo:j_hi)
	endif
	endif   ! thereis_metals

	!-- H3+ infrared cooling (molecular layer) --!

	! Optically thin rotational-vibrational emission of H3+, Miller et al.
	! (2013) LTE emission per molecule with their Table-6 non-LTE departure
	! factor s(T, n_H2) (h3p_cooling module -- ONE definition, shared with
	! every caller of eval_cool). Inside a molecular base at T ~ 1e3 K this is
	! the dominant coolant: the atomic channels above are all exponentially
	! suppressed there, so leaving it out of `cool` leaves that gas with no
	! radiative loss at all. It lives here, not on top of eval_cool's return
	! value, so that the temperature update in the marching loop
	! (energy_semi_implicit) and the steady residual (ioniz_eq ->
	! assemble_residual) balance the SAME cooling function.
	! Zero for callers that model a molecule-free gas (nmol absent) and for
	! cells with no H3+; the rate is exactly proportional to n_H3+, so
	! skipping those cells is not an approximation.
	! With "Base IR field" on the same 3-4 um bands also ABSORB the thermal
	! infrared of the lower atmosphere, a blackbody at T0 covering half the
	! sky at the base; W_dil = 0 (the default) leaves the emission-only rate
	! untouched. Approximations and their range: h3p_net_cooling_rate.
	!$ if (ec_prof_on) then
	!$    ect(4) = ect(4) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif
	cool_H3p(j_lo:j_hi) = 0.0d0
	if (present(nmol)) then
		do j = j_lo,j_hi
			if (nmol(j,3) .ne. 0.0d0) then
				if (base_ir_field) then
					cool_H3p(j) = h3p_net_cooling_rate(T_K(j),         &
					                 nmol(j,3), nmol(j,1), T0,         &
					                 0.5d0*base_sky_fraction(r(j),1.0d0))
				else
					cool_H3p(j) = h3p_cooling_rate(T_K(j), nmol(j,3),  &
					                               nmol(j,1))
				endif
			endif
		enddo
	endif

	!-- Molecular infrared bands (H2 lines, H2O and CO bands) --!

	! The infrared coolants a real H2 atmosphere carries below the H2 -> H
	! front and this code did not: the H2 quadrupole plus magnetic dipole line
	! spectrum (Roueff et al. 2019) and the H2O and CO vibration-rotation bands
	! (HITEMP through the Photochem k-coefficients). Each is the NET rate --
	! LTE emission minus absorption of the diluted B_nu(T0) the lower
	! atmosphere presents -- so each vanishes at its own radiative equilibrium
	! temperature instead of running the layer down to nothing. That fixed
	! point is what TO_BE_DONE.md item (G) asks for; the emission magnitudes
	! alone would only deepen the collapse.
	! Off by default (`Molecular IR bands`). The dilution is the same
	! 0.5*base_sky_fraction the H3+ closure uses, and it is zero when
	! `Base IR field` is off, which reduces the channels to pure emitters --
	! the configuration input_read warns about.
	! Physics, validity range and sources: molecular_infrared_cooling.f90.
	cool_H2(j_lo:j_hi)  = 0.0d0
	cool_H2O(j_lo:j_hi) = 0.0d0
	cool_CO(j_lo:j_hi)  = 0.0d0
	if (mol_ir_bands) then
		do j = j_lo,j_hi
			if (base_ir_field) then
				w_ir = 0.5d0*base_sky_fraction(r(j),1.0d0)
			else
				w_ir = 0.0d0
			endif
			if (present(nmol)) cool_H2(j) =                            &
				h2_line_net_cooling_rate(T_K(j), nmol(j,1), w_ir)
			if (present(nox)) then
				cool_H2O(j) = h2o_band_net_cooling_rate(T_K(j),        &
				                                        nox(j,2), w_ir)
				cool_CO(j)  = co_band_net_cooling_rate(T_K(j),         &
				                                       nox(j,3), w_ir)
			endif
		enddo
	endif

	! Total cooling rate
	cool(j_lo:j_hi) = ne(j_lo:j_hi)*(brem(j_lo:j_hi) + coex(j_lo:j_hi)  &
	                  + reco(j_lo:j_hi) + coio(j_lo:j_hi))              &
	       + cool_M(j_lo:j_hi) + cool_H3p(j_lo:j_hi)                    &
	       + cool_H2(j_lo:j_hi) + cool_H2O(j_lo:j_hi) + cool_CO(j_lo:j_hi)

	! Breakdown by channel for the diagnostic (Huang Fig. 10).
	! Read straight from the arrays already computed above, so the sum of
	! all channels reproduces `cool` exactly in the default branch.
	if (present(cool_chan)) then
		cool_chan(j_lo:j_hi,1) = ne(j_lo:j_hi)*reco(j_lo:j_hi)
		cool_chan(j_lo:j_hi,2) = ne(j_lo:j_hi)*coio(j_lo:j_hi)
		! Collisional excitation split by absorber. H I is the
		! Lyman-alpha-dominated H-line cooling; He II is its own term. The
		! He I column is taken as the remainder ne*coex - HI - HeII so that it
		! also absorbs the He 2^3S metastable terms folded into coex above,
		! keeping columns 3-5 an exact split of ne*coex.
		cool_chan(j_lo:j_hi,3) = ne(j_lo:j_hi)*(coeff_coex_rate_HI(j_lo:j_hi)*nhi(j_lo:j_hi))      ! coex_HI [Lya]
		cool_chan(j_lo:j_hi,5) = ne(j_lo:j_hi)*(coeff_coex_rate_HeII(j_lo:j_hi)*nheii(j_lo:j_hi))  ! coex_HeII
		cool_chan(j_lo:j_hi,4) = ne(j_lo:j_hi)*coex(j_lo:j_hi) - cool_chan(j_lo:j_hi,3) - cool_chan(j_lo:j_hi,5)  ! coex_HeI
		cool_chan(j_lo:j_hi,6) = ne(j_lo:j_hi)*brem(j_lo:j_hi)
		cool_chan(j_lo:j_hi,7)  = cool_H3p(j_lo:j_hi)
		cool_chan(j_lo:j_hi,8)  = cool_H2(j_lo:j_hi)
		cool_chan(j_lo:j_hi,9)  = cool_H2O(j_lo:j_hi)
		cool_chan(j_lo:j_hi,10) = cool_CO(j_lo:j_hi)
		do i = 1,n_mion
			if (mion_iscool(i)) then
				cool_chan(j_lo:j_hi,10+i) = ne(j_lo:j_hi)*nm(j_lo:j_hi,i)*c_metal(j_lo:j_hi,i)
			else
				cool_chan(j_lo:j_hi,10+i) = 0.0d0
			endif
		enddo
	endif

	! End of subroutine
	!$ if (ec_prof_on) then
	!$    ect(5) = ect(5) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif
	end subroutine eval_cool_cells

	! ------------------------------------------------------------- !

	subroutine write_cool_breakdown_eq(T_in,n_in,f_sp_in)
	! Diagnostic. Dump the radiative cooling rate in each channel vs
	! radius for the converged equilibrium state, reusing eval_cool's exact
	! coefficients (no offline re-derivation). Columns: H/He recombination,
	! collisional ionization, collisional excitation, bremsstrahlung, H3+
	! infrared cooling, then one column per metal ion line-cooling channel
	! (canonical species_table order). All in cgs erg cm^-3 s^-1; the channel
	! sum reproduces the Hydro_ioniz.txt `cool` column. The printed max
	! relative residual is the internal consistency check. Lets the user
	! identify the dominant coolant in 1.15 <~ r/Rp <~ 1.4 against Huang et
	! al. (2023) Fig. 10, and the H3+ column the coolant of the molecular base.

	integer :: j,i,im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) :: n_dim,T_K,ne
	real*8, dimension(1-Ng:N+Ng) :: nhi,nhii,nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	! Molecular densities [cm^-3] (H2, H2+, H3+, HeH+); zero for an atomic run
	real*8, dimension(1-Ng:N+Ng,4) :: nmol
	! Oxygen carriers [cm^-3] (OH, H2O, CO); zero without the oxygen chemistry
	real*8, dimension(1-Ng:N+Ng,3) :: nox
	! Planck-mean optical depth of the H2O and CO columns above each cell
	real*8, dimension(1-Ng:N+Ng) :: tau_h2o,tau_co
	real*8, dimension(1-Ng:N+Ng) :: cool,csum,rel
	real*8, dimension(1-Ng:N+Ng,10+n_mion) :: chan
	! Throwaway eval_cool rate outputs (not needed for the dump)
	real*8, dimension(1-Ng:N+Ng) :: rchiiB,rcheiiB,rcheiiiB
	real*8, dimension(1-Ng:N+Ng) :: a_ion_HI,a_ion_HeI,a_ion_HeII
	real*8, dimension(1-Ng:N+Ng,n_mion) :: rec_m,aion_m
	real*8 :: maxrel

	! Dimensionalize exactly as ioniz_eq does
	n_dim = n_in*n0
	T_K   = T_in*T0
	nhi   = f_sp_in(:,1)*n_dim
	nhii  = f_sp_in(:,2)*n_dim
	if (thereis_He) then
		nhei   = f_sp_in(:,3)*n_dim
		nheii  = f_sp_in(:,4)*n_dim
		nheiii = f_sp_in(:,5)*n_dim
	else
		nhei = 0.0d0; nheii = 0.0d0; nheiii = 0.0d0
	endif
	! He 2^3S metastable, so the collisional-ionization channel matches the
	! main cool column (0 when the triplet is not tracked).
	if (thereis_HeITR) then
		nheiTR = f_sp_in(:,6)*n_dim
	else
		nheiTR = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	! Molecular ions, so the dumped ne and the cooling channels are the ones
	! ioniz_eq itself uses (all zero for an atomic run).
	nmol = 0.0d0
	if (thereis_mol) then
		nmol(:,1) = f_sp_in(:,isp_H2)  *n_dim
		nmol(:,2) = f_sp_in(:,isp_H2p) *n_dim
		nmol(:,3) = f_sp_in(:,isp_H3p) *n_dim
		nmol(:,4) = f_sp_in(:,isp_HeHp)*n_dim
	endif
	! Oxygen carriers, so the H2O and CO infrared columns are the ones
	! ioniz_eq itself used (all zero without the oxygen chemistry).
	nox = 0.0d0
	if (thereis_oxychem) then
		nox(:,1) = f_sp_in(:,isp_OH) *n_dim
		nox(:,2) = f_sp_in(:,isp_H2O)*n_dim
		nox(:,3) = f_sp_in(:,isp_CO) *n_dim
	endif
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol)

	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,            &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,                &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
				   cool, cool_chan = chan, nheiTR = nheiTR,       &
				   nmol = nmol, nox = nox)

	! Internal consistency: channel sum vs total cool (default branch -> ~eps)
	csum = 0.0d0
	do i = 1,10+n_mion
		csum = csum + chan(:,i)
	enddo
	rel    = abs(csum - cool)/max(abs(cool),1.0d-99)
	maxrel = maxval(rel(1:N))
	write(*,'(a,es9.2)')                                                  &
		' (write_cool_breakdown_eq) max |sum(channels)/cool - 1| = ', maxrel

	open(unit = 71, file = './output/Cooling_breakdown.txt')
	write(71,'(a)') '# Radiative cooling rate in each channel [cgs erg cm^-3 s^-1] vs radius.'
	write(71,'(a)') '# Channel sum reproduces the Hydro_ioniz.txt cool column.'
	write(71,'(a)') '# col1 r/Rp  col2 T[K]  col3 ne  col4 cool_total  col5 reco'  &
	             // '  col6 coio  col7 coex_HI[Lya]  col8 coex_HeI  col9 coex_HeII'  &
	             // '  col10 brem  col11 H3p_IR  col12 H2_IR  col13 H2O_IR'          &
	             // '  col14 CO_IR  then one col per metal ion:'
	write(71,'(a)') '#   the three molecular-band columns are NET rates'          &
	             // ' (emission minus absorption of the field from below), so'
	write(71,'(a)') '#   a negative value is that band heating the gas; they'     &
	             // ' are identically zero unless "Molecular IR bands" is on.'
	write(71,'(a)',advance='no') '#   metal-ion columns (canonical order):'
	do i = 1,n_mion
		write(71,'(1x,a)',advance='no') trim(mion_name(i))
	enddo
	write(71,*)
	call write_row_layout_header(71)
	do j = 1-Ng,N+Ng
		write(71,*) r(j), T_K(j), ne(j), cool(j),                        &
		            chan(j,1), chan(j,2), chan(j,3), chan(j,4),           &
		            chan(j,5), chan(j,6), chan(j,7), chan(j,8),           &
		            chan(j,9), chan(j,10),                                &
		            (chan(j,10+i), i = 1,n_mion)
	enddo
	! Whether the optically thin limit the molecular bands are computed in is
	! actually where the run sits: the Planck-mean optical depth of the H2O
	! and CO columns, integrated from each cell to the top of the domain, and
	! its maximum over the wind. This is a measurement of the approximation
	! stated at molecular_infrared_cooling.f90, not an input to it.
	if (mol_ir_bands .and. thereis_oxychem) then
		tau_h2o = 0.0d0
		tau_co  = 0.0d0
		do j = N+Ng-1,1-Ng,-1
			tau_h2o(j) = tau_h2o(j+1)                                    &
			  + molecular_planck_cross_section(1, T_K(j))*nox(j,2)        &
			    *(r(j+1) - r(j))*R0
			tau_co(j)  = tau_co(j+1)                                     &
			  + molecular_planck_cross_section(2, T_K(j))*nox(j,3)        &
			    *(r(j+1) - r(j))*R0
		enddo
		write(71,'(a,es10.3,a,f9.5)') '# max Planck-mean tau(H2O) to the'     &
		   // ' top of the domain = ', maxval(tau_h2o(1:N)),                  &
		   ' at r/Rp = ', r(maxloc(tau_h2o(1:N),1))
		write(71,'(a,es10.3,a,f9.5)') '# max Planck-mean tau(CO)  to the'     &
		   // ' top of the domain = ', maxval(tau_co(1:N)),                   &
		   ' at r/Rp = ', r(maxloc(tau_co(1:N),1))
		write(*,'(a,2es10.3)') ' (write_cool_breakdown_eq) max Planck-mean'   &
		   // ' tau(H2O), tau(CO) upward = ', maxval(tau_h2o(1:N)),           &
		   maxval(tau_co(1:N))
	endif
	close(71)

	end subroutine write_cool_breakdown_eq

	! ------------------------------------------------------------- !

	subroutine write_heat_breakdown_eq(T_in,n_in,f_sp_in)
	! Diagnostic. Dump the volumetric heating rate in each channel vs
	! radius for the converged equilibrium state, recomputing the same
	! photoheating (PH_heat_HHe) the solver uses plus the excited-H Balmer,
	! He-recombination, Penning and Lyman-Werner heating terms added in
	! ionization_equilibrium. All in cgs erg cm^-3 s^-1; the channel sum
	! reproduces the total heating (heat_total column) and, up to convergence,
	! the Hydro_ioniz.txt heat column. Photoheating columns: H I, He I,
	! He II, He 2^3S, H2, metals (sum over photo-ionizable metal ions). Then
	! the excited-H photoelectric (Hpe) and Lyman-alpha de-excitation (Hdx)
	! heating, He-recombination-driven H heating, He(2^3S)+H and He(2^3S)+H2
	! Penning heating, and H2 Lyman-Werner photodissociation heating. The last
	! three molecular-run channels (H2 photoheating, He(2^3S)+H2 Penning,
	! Lyman-Werner) are identically zero for an atomic run. The printed max
	! relative residual is the internal consistency check on the photoheating
	! split.

	integer :: j,im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) :: n_dim,T_K,ne
	real*8, dimension(1-Ng:N+Ng) :: nhi,nhii,nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng) :: nh,nhe,xion
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	! Molecular densities [cm^-3] (H2, H2+, H3+, HeH+); zero for an atomic run
	real*8, dimension(1-Ng:N+Ng,4) :: nmol
	real*8, dimension(1-Ng:N+Ng,6) :: hchan
	real*8, dimension(1-Ng:N+Ng) :: heat_ph,heat_tot,csum,rel
	real*8, dimension(1-Ng:N+Ng) :: h_hrc,h_penning,h_penning_h2,h_lw
	! Chemical heat of the collisional molecular reactions, and the total
	! gas-particle density (electrons excluded) the three-body rates in it
	! read as the third body M.
	real*8, dimension(1-Ng:N+Ng) :: h_chem, ntot_dump
	! Oxygen-chemistry carriers [cm^-3] (OH, H2O, CO; zero without the
	! option), their star-ward FUV columns and the photolysis heating.
	real*8, dimension(1-Ng:N+Ng,3) :: nox
	real*8, dimension(1-Ng:N+Ng) :: NH2Oc,NOHc,h_fuv
	! Scratch for the shared-beam field routine: this dump needs only the
	! two heating terms out of it, but the routine solves the whole field.
	real*8, dimension(1-Ng:N+Ng) :: fsh_dump, trl_dump
	real*8, dimension(1-Ng:N+Ng,n_fuv_band) :: tau_dump, jh2o_dump, joh_dump
	real*8 :: a_dump, ovl_dump
	integer :: ib
	! Throwaway PH_heat_HHe rate outputs (not needed for the dump)
	real*8, dimension(1-Ng:N+Ng) :: P_HI,P_HeI,P_HeII,P_HeITR,q,P_H2
	real*8, dimension(1-Ng:N+Ng) :: P_H2_di
	! Collisional quench fraction of an H2 vibrational excitation and the mean
	! internal energy one B/C fluorescence leaves in X: PH_heat_HHe needs both
	! from the caller because it does not carry the temperature, and without
	! them its two H2 vibrational heat channels are silently off.
	real*8, dimension(1-Ng:N+Ng) :: f_vib_quench,e_vib_bound
	real*8, dimension(1-Ng:N+Ng,n_mion) :: P_m
	! Star-ward H2 column [cm^-2] and its Lyman-Werner dissociation rate [s^-1]
	real*8, dimension(1-Ng:N+Ng) :: NH2col,k_lw,p_lw_s,p_lw_a
	! He-recombination coupling / triplet scratch
	real*8, dimension(1-Ng:N+Ng) :: rcheiTR,rcheii,q13,q31a,q31b,Q31
	real*8, dimension(1-Ng:N+Ng) :: rcheiiB_hrc,dP_HI_hrc,dP_H2_hrc
	real*8, dimension(1-Ng:N+Ng,n_mion) :: dP_m_hrc
	real*8 :: A31,maxrel

	! Dimensionalize exactly as ioniz_eq does
	n_dim = n_in*n0
	T_K   = T_in*T0
	nhi   = f_sp_in(:,1)*n_dim
	nhii  = f_sp_in(:,2)*n_dim
	if (thereis_He) then
		nhei   = f_sp_in(:,3)*n_dim
		nheii  = f_sp_in(:,4)*n_dim
		nheiii = f_sp_in(:,5)*n_dim
	else
		nhei = 0.0d0; nheii = 0.0d0; nheiii = 0.0d0
	endif
	if (thereis_HeITR) then
		nheiTR = f_sp_in(:,6)*n_dim
	else
		nheiTR = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	! Oxygen carriers, so the FUV photolysis heating channel is the one
	! ioniz_eq itself applies (all zero without the oxygen chemistry).
	nox = 0.0d0
	if (thereis_oxychem) then
		nox(:,1) = f_sp_in(:,isp_OH) *n_dim
		nox(:,2) = f_sp_in(:,isp_H2O)*n_dim
		nox(:,3) = f_sp_in(:,isp_CO) *n_dim
	endif
	! Molecular ions, so the dumped ne and the molecular heating channels are
	! the ones ioniz_eq itself uses (all zero for an atomic run).
	nmol = 0.0d0
	if (thereis_mol) then
		nmol(:,1) = f_sp_in(:,isp_H2)  *n_dim
		nmol(:,2) = f_sp_in(:,isp_H2p) *n_dim
		nmol(:,3) = f_sp_in(:,isp_H3p) *n_dim
		nmol(:,4) = f_sp_in(:,isp_HeHp)*n_dim
	endif
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol)

	! H and He NUCLEI totals and the ionized fraction handed to the
	! photoelectron partition: the TOTAL free electron density over those
	! nuclei, the quantity Dalgarno, Yan & Liu (1999) section 7 define.
	! The one shared definition ionization_equilibrium itself calls, so the
	! oxygen carriers OH and H2O contribute their H nuclei here too -- both
	! to xion and to the H nucleus density the FUV/Lyman-Werner beam below
	! is scaled by.
	call hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,       &
	                                    nh,nhe,nmol,nox)
	xion = min(max(ne/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Photoheating split (same call the solver makes; H2 adds its opacity and
	! its photoelectric heating, hchan(:,5), when the run is molecular). The
	! absorber columns are weighted by the opa_pf the last equilibrium pass
	! left in place, so they are the solver's own up to convergence.
	if (thereis_He) then
		if (thereis_mol) then
			! Same arguments the solver passes in ionization_equilibrium.
			! f_vibq/e_vibq are what switch on the two H2 vibrational heat
			! channels; dropping them made this dump under-report the
			! molecular photoheating by up to a factor 2 in the H2 layer.
			f_vib_quench = h2_vibrational_heat_fraction(T_K, nhi,     &
			                                           nmol(:,1))
			e_vib_bound  = h2_energy_per_bound_fluorescence_eV(T_K)
			call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,         &
			         P_HI,P_HeI,P_HeII,P_HeITR, P_m,                  &
			         heat_ph,q, nmol(:,1),P_H2, P_H2_di,              &
			         f_vibq = f_vib_quench, e_vibq = e_vib_bound,     &
			         heat_chan = hchan)
		else
			call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,         &
			         P_HI,P_HeI,P_HeII,P_HeITR, P_m,                  &
			         heat_ph,q, heat_chan = hchan)
		endif
	else
		call PH_heat_H(nhi, xion, P_HI, heat_ph, q)
		hchan = 0.0d0
		hchan(:,1) = heat_ph
	endif

	! He-recombination-driven H heating and He(2^3S)+H Penning heating,
	! reconstructed exactly as ionization_equilibrium adds them to `heat`.
	A31 = 0.0d0; q31a = 0.0d0; q31b = 0.0d0; Q31 = 0.0d0
	if (thereis_HeITR) &
		call HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31a,q31b,Q31)
	h_hrc = 0.0d0
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nmol(:,1), nhei, nheii, nheiTR, &
		                     ne, nm, A31, q31a, q31b,                  &
		                     rcheiiB_hrc, dP_HI_hrc, dP_H2_hrc,        &
		                     dP_m_hrc, h_hrc)
	endif
	h_penning = 0.0d0
	if (thereis_HeITR) h_penning =                                    &
		f_penning_HeI23S*nheiTR*nhi*Q31                                &
		*(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

	! He(2^3S)+H2 -> He(1^1S)+H2+ + e- Penning ionization heating: the electron
	! carries away (e_th_HeI - e_th_HeTR) - e_th_H2 = 4.4 eV. ioniz_HeI23S_H2
	! is the total, so the Penning branch alone carries this exothermicity.
	! Zero unless a molecular run also tracks the triplet.
	h_penning_h2 = 0.0d0
	if (thereis_mol .and. thereis_HeITR) h_penning_h2 =               &
		f_penning_HeI23S*nheiTR*nmol(:,1)*ioniz_HeI23S_H2(T_K)         &
		*((e_th_HeI - e_th_HeTR) - e_th_H2)/erg2eV

	! Lyman-Werner heating, both halves of it: the 13.5% of pumps that
	! dissociate leave the fragment pair about 0.4 eV of kinetic energy (the
	! 4.48 eV bond energy is paid by the photon, not by the gas), and the
	! other 86.5% fluoresce back into vibrationally excited bound levels
	! whose energy is collisionally de-excited into heat at these densities.
	! Same two terms, in the same order, as the energy equation of
	! ionization_equilibrium. And the FUV photolysis heating of
	! H2O and OH (the oxygen chemistry), whose ledger is the same: the
	! excess of the absorbed photon over the bond energy. Both come from the
	! ONE shared-beam field routine the equilibrium solve uses, rebuilt here
	! from f_sp_in because this dump must be a function of the state it is
	! given rather than of the last state solved.
	h_lw  = 0.0d0
	h_fuv = 0.0d0
	if (thereis_mol .or. thereis_oxychem) then
		call fuv_lw_photon_field(nmol(:,1), nox(:,2), nox(:,1), T_K,   &
		                         nh,                                   &
		                         NH2col, NH2Oc, NOHc,                  &
		                         fsh_dump, trl_dump, tau_dump,         &
		                         k_lw, p_lw_s, p_lw_a,                 &
		                         jh2o_dump, joh_dump, a_dump,          &
		                         ovl_dump)
		h_lw = k_lw*nmol(:,1)*e_lw_fragment_erg                       &
		     + k_lw*(1.0d0 - p_lw_s)/max(p_lw_s, 1.0d-30)*nmol(:,1)   &
		       *h2_energy_per_bound_fluorescence_erg(T_K)             &
		       *h2_vibrational_heat_fraction(T_K, nhi, nmol(:,1))
		do ib = 1,n_fuv_band
			h_fuv = h_fuv                                              &
			  + jh2o_dump(:,ib)*nox(:,2)                               &
			    *heat_per_water_dissociation(ib)                       &
			  + joh_dump(:,ib)*nox(:,1)                                &
			    *heat_per_hydroxyl_dissociation(ib)
		enddo
	endif

	! Total heating (independent of the channel columns; the residual
	! below checks the photoheating decomposition against heat_ph).
	! Chemical heat of the collisional molecular reactions, rebuilt here
	! from the state this dump was given, exactly as ionization_equilibrium
	! adds it to `heat`. Zero unless the run switched it on.
	h_chem = 0.0d0
	call calc_ntot(nhi,nhii,nhei,nheii,nheiii,ntot_dump,nm,nmol,nox)
	if (thereis_mol .and. mol_reaction_heat)                           &
		call molecular_chemical_heating(T_K, nhi, nhii, nheii, nmol,   &
		                                ne, ntot_dump, h_chem)

	heat_tot = heat_ph + Hpe_arr + Hdx_arr + h_hrc + h_penning        &
	         + h_penning_h2 + h_lw + h_fuv + h_chem

	! Internal consistency of the photoheating split.
	csum   = hchan(:,1) + hchan(:,2) + hchan(:,3) + hchan(:,4)        &
	       + hchan(:,5) + hchan(:,6)
	rel    = abs(csum - heat_ph)/max(abs(heat_ph),1.0d-99)
	maxrel = maxval(rel(1:N))
	write(*,'(a,es9.2)')                                                  &
		' (write_heat_breakdown_eq) max |sum(photo channels)/heat_photo - 1| = ', maxrel

	open(unit = 72, file = './output/Heating_breakdown.txt')
	write(72,'(a)') '# Volumetric heating rate in each channel [cgs erg cm^-3 s^-1] vs radius.'
	write(72,'(a)') '# Channel sum reproduces the heat_total column (and the Hydro_ioniz.txt'  &
	             // ' heat column up to convergence).'
	write(72,'(a)') '# col1 r/Rp  col2 T[K]  col3 ne  col4 heat_total  col5 heat_HI'  &
	             // '  col6 heat_HeI  col7 heat_HeII  col8 heat_He23S  col9 heat_H2'  &
	             // '  col10 heat_metals  col11 heat_Hpe[excitedH]'                    &
	             // '  col12 heat_Hdx[Lya-deexc]  col13 heat_He_recomb  col14 heat_He23S_Penning'  &
	             // '  col15 heat_He23S_H2_Penning  col16 heat_H2_LW'   &
	             // '  col17 heat_FUV_photolysis  col18 heat_mol_chem'
	call write_row_layout_header(72)
	do j = 1-Ng,N+Ng
		write(72,*) r(j), T_K(j), ne(j), heat_tot(j),                    &
		            hchan(j,1), hchan(j,2), hchan(j,3), hchan(j,4),       &
		            hchan(j,5), hchan(j,6),                               &
		            Hpe_arr(j), Hdx_arr(j), h_hrc(j), h_penning(j),       &
		            h_penning_h2(j), h_lw(j), h_fuv(j), h_chem(j)
	enddo
	close(72)

	end subroutine write_heat_breakdown_eq

	!----------------------------------!

	! Various coefficients for HeI triplet chemistry
	subroutine HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31a,q31b,Q31)
	
	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

	real*8, intent(out) :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiTR,rcheii,   &
								   q13,q31a,q31b,Q31

	call rec_HeII_23S(T_K,rcheiTR)
	call rec_HeII_11S(T_K,rcheii)
	call coex_HeI_1S_23S(T_K,q13)
	call coex_HeI_23S_21S(T_K,q31a)
	call coex_HeI_23S_21P(T_K,q31b)
	A31 = 1.272e-4

	! He(2^3S)+H total ionization (Penning + associative), Garcia Munoz (2025)
	! from the Movre & Meyer (1997) cross sections, used unconditionally (the
	! legacy temperature-independent 5e-10 constant is no longer selectable).
	! Q31 is the TOTAL: it is what removes the metastable. The consumers that
	! create a lasting proton or deposit the Penning exothermicity scale it by
	! f_penning_HeI23S.
	Q31 = ioniz_HeI23S_H(T_K)

   ! End of subroutine
	end subroutine HeITR_coeffs

	! ------------------------------------------------------------- !

	! Fraction of the photons emitted in a cell that the cell itself absorbs,
	! 1 - exp(-tau), for a cell optical depth tau >= 0. Evaluated by its series
	! below tau = 1e-4, where 1 - exp(-tau) loses all but a few digits to
	! cancellation; the series is truncated where its next term is below the
	! double-precision epsilon of the result. Returns 0 for a non-positive tau,
	! which the ionization solve can hand it as a small negative trace density.
	pure function absorbed_photon_fraction(tau) result(f_abs)
	real*8, intent(in) :: tau
	real*8 :: f_abs
	if (tau .le. 0.0d0) then
		f_abs = 0.0d0
	else if (tau .lt. 1.0d-4) then
		f_abs = tau*(1.0d0 - 0.5d0*tau*(1.0d0 - tau/3.0d0))
	else
		f_abs = 1.0d0 - exp(-tau)
	endif
	end function absorbed_photon_fraction

	! ------------------------------------------------------------- !

	! He recombination radiation ionizing H I and H2 (Draine 2011 on-the-spot
	! emission, absorbed locally; docs/QUESTIONS_2026-07-17.md,
	! docs/supersonic_molecular_base.md section
	! 12). Given the pre-solve (lagged) densities and the rate coefficients,
	! returns the He II recombination coefficient the ionization balance should
	! use (rcheiiB_new), the extra H I and H2 photoionization rates [s^-1]
	! (dP_HI, dP_H2), and the extra photoelectron heating [erg cm^-3 s^-1]
	! (dheat). All four are zero when the flag is off; the caller applies them.
	!
	! WHICH SPECIES ABSORBS THE PHOTON. Every channel below emits a photon of a
	! known energy E_c, and the on-the-spot assumption is that the photon is
	! absorbed inside the same cell. WHO absorbs it is then a competition among
	! the species whose ionization threshold lies below E_c, in the ratio of
	! their absorption coefficients,
	!
	!     w_s(E_c) = n_s sigma_s(E_c) / sum_s' n_s' sigma_s'(E_c),
	!
	! with s running over H I (13.6 eV), H2 (15.4 eV), every photo-ionizable
	! metal ion whose threshold lies below E_c, and, for the >= 24.6 eV
	! ground-capture continuum only, He I. The H I share ionizes H I, the H2
	! share ionizes H2 (it is returned in dP_H2 and added to the H2
	! photoionization rate by the caller), each metal share ionizes that metal
	! ion (returned in dP_m, added to its photoionization rate P_m), and the
	! He I share re-ionizes He and is therefore excluded from the effective
	! He II recombination coefficient.
	! In a molecular gas H2 is not a small competitor: sigma_H2/sigma_HI is 1.2
	! at 16 eV and 3.5-3.7 from 20 to 25 eV, so where H2 outnumbers H I the He
	! recombination photons go to H2, not to H I.
	!
	! THE METALS ARE NOT NEGLIGIBLE EITHER (item P34). Their abundance is
	! 1e-4 to 1e-3, but their cross sections in this band are large: at
	! 24.6 eV sigma is 4.9 (C II), 5.1 (Fe II) and 11.9 (O I) against 1.24
	! for H I, and at 16.11 eV 13.4 (C I) and 6.6 (O I) against 4.0. Measured
	! on the converged WASP-121b cases the metal share of the absorption
	! coefficient is 3e-3 to 4e-3 at the base (O I) and reaches 2.3e-2 at
	! 1.61 R_p, where hydrogen is ionized and C II carries 89% of it; on the
	! hot-Uranus molecular cases it is 3e-3 in the molecular layer and 1.1e-2
	! in the wind. That is a share of the photons, and it was going to
	! hydrogen.
	!
	! The competition is written in the ratio form w_HI = 1/(1 + R_He n_HeI/n_HI
	! + R_H2 n_H2/n_HI + sum_m R_m n_m/n_HI), with R_s = sigma_s(E_c)/
	! sigma_HI(E_c) evaluated once at entry. A run with no metals adds
	! EXACTLY 0.0 to every one of those sums, which is what keeps it bitwise
	! unchanged.
	!
	! HOW MANY OF THE PHOTONS STAY. The cell does not necessarily absorb the
	! photons it emits. Each channel has its own cell optical depth and its own
	! absorbed fraction,
	!
	!     tau_c  = (sum_s n_s sigma_s(E_c)) dr_j,      dr_j = cell width [cm],
	!     f_abs  = 1 - exp(-tau_c),
	!
	! the remaining 1 - f_abs leaving the cell, and the rate per particle of
	! absorber s is the absorbed photons shared out by absorption coefficient,
	!
	!     dP_s = P_c f_abs w_s(E_c)/n_s
	!          = P_c sigma_s(E_c) dr_j [1 - exp(-tau_c)]/tau_c.
	!
	! The second form is the one to reason about: the bracket is 1 at tau -> 0
	! and 1/tau at tau -> infinity, so the rate has NO division singularity
	! anywhere -- at tau -> 0 it is the optically thin limit P_c sigma_s dr_j,
	! the rate a single absorber in a transparent cell would see. It is
	! evaluated in the first form, share divided by density and multiplied by
	! f_abs, because f_abs is exactly 1.0 in double precision for tau > 37 and
	! the expression then reduces bit for bit to the optically thick limit
	! (which is what a molecular base and a dense atomic base are). The
	! vanishing-absorber divergence this replaces is item (P40): with n_HI = 0
	! and n_H2 a 1e-151 numerical residue the thick-limit expression alone
	! returns 1e100 s^-1, while tau_c = 1e-160 and the photon has in fact left.
	!
	! An escaped photon does not re-ionize He either, so the He I share is
	! removed from the effective He II recombination coefficient only in
	! proportion to the photons that stay: alpha_eff = alpha_B + alpha_1
	! (1 - f_abs w_HeI), written below as (1 - f_abs) + f_abs (w_HI + w_H2) so
	! that the optically thick limit is again bit for bit the earlier form.
	!
	! VALIDITY. (i) LOCAL absorption: a photon that leaves the cell is dropped
	! rather than followed, so its absorption in some outer cell is not counted.
	! The wind thins outward, so most of what leaves a thin cell escapes the
	! domain; in a thick layer, where the neighbours would absorb it, f_abs is
	! already 1 and nothing leaves. (ii) A photoelectron released by one of
	! these photons is deposited entirely as heat: it is not passed through
	! the Shull & van Steenberg secondary-ionization partition that the
	! stellar photoelectrons of PH_heat_HHe go through. The energies involved
	! are 0.7-11 eV, at or below the 11 eV threshold below which that
	! partition returns pure heat anyway, except for the 11.0 eV H I share of
	! the ground-capture channel. (iii) Each channel is collapsed onto one representative photon energy; the
	! 2^1S two-photon continuum, which is not a line, is discussed at E_2q_HeI
	! below. (iv) tau_c uses the cell width dr_j, i.e. a radially escaping
	! photon; the true mean chord of an isotropically emitted photon in a plane
	! slab is longer by a factor of order 2.
	!
	! atomic mode (thereis_HeITR = .false.):
	!   alpha_1 = Mao & Kaastra ground (1s^2) capture; alpha_B = active He II
	!   case B (rec_HeII_B). rcheiiB_new = alpha_B + (w_HI + w_H2) alpha_1
	!   (Draine Eq. 14.17 with H2 added to the competition); the case-B cascade
	!   fraction z alpha_B is split between H I and H2 at the cascade-averaged
	!   photon energy E_casc_HeI; the heating uses the ground photoelectron
	!   energies E_gnd - 13.6 = 11.0 eV (H I) and E_gnd - 15.4 = 9.2 eV (H2) and
	!   the cascade-averaged Ee_casc_HeI ~ 6.14 eV (H I) / EeH2_casc_HeI (H2).
	! TR mode (thereis_HeITR = .true.): channels are explicit -- alpha_1 is the
	!   1^1S channel (rec_HeII_11S, which HeITR_coeffs writes into rcheiiB) and
	!   the 2^3S / singlet cascade rates (A31, q31a, q31b, n_2^3S = nheiTR) drive
	!   the H-ionizing photon production directly, so z is not used (a 2^3S
	!   destroyed by photoionization/Penning emits no 19.8 eV photon).
	subroutine he_rec_coupling(T_K, nhi, nh2, nhei, nheii, nheiTR, ne,   &
	                           nm, A31, q31a, q31b,                      &
	                           rcheiiB_new, dP_HI, dP_H2, dP_m, dheat)
	! He recombination radiation absorbed locally by H I, H2 and the metal
	! ions (Draine 2011 emission; see the call site in
	! ionization_equilibrium for the physics).
	!
	! What is returned, and at which state. rcheiiB_new, dP_HI, dP_H2 and
	! dP_m are RATES and coefficients: an effective He II recombination
	! coefficient [cm^3 s^-1] and the photoionization rate [s^-1] each
	! absorber takes, i.e. the quantities the ionization system is solved
	! with. dheat is the photoelectron HEATING [erg cm^-3 s^-1] of the
	! composition passed in, since every channel of it carries the density
	! of the absorber that receives the photon. A caller whose composition
	! changes afterwards must call this again at the new densities and take
	! dheat from that call; ionization_equilibrium does exactly that, so the
	! rates its sweep uses stay the lagged ones while the heating it returns
	! belongs to the composition it returns.

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi, nh2, nhei,   &
	                                              nheii, nheiTR, ne,     &
	                                              q31a, q31b
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in)  :: nm
	real*8,                       intent(in)  :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiiB_new, dP_HI,    &
	                                              dP_H2, dheat
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) :: dP_m

	real*8, dimension(1-Ng:N+Ng) :: alpha1, alphaB
	real*8 :: Rratio, R2_gnd, R2_584, R2_19, R2_2q, R2_casc
	real*8 :: s_gnd, s_584, s_19, s_2q, s_casc
	real*8 :: s2_gnd, s2_584, s2_19, s2_2q, s2_casc, sHe_gnd
	real*8 :: dl, f_gnd, f_584, f_19, f_2q, f_casc, y_net
	real*8 :: y, y2, u2, z, T4, ncrit, P_add, P2_add, H_add, H2_add
	real*8 :: wc, wc2, w584, w584_2, w19, w19_2, w2q, w2q_2
	real*8 :: fs_HI, fs_R2, Ees_HI, Ees_H2, nhig, n2
	logical :: absorb_low, absorb_gnd
	integer :: j
	! --- metal absorbers (item P34) ---
	! Channel index: 1 ground capture (24.6 eV), 2 the 584 A resonance
	! (21.2), 3 the 2^3S line (19.8), 4 the 2^1S two-photon continuum
	! (16.110), 5 the case-B cascade average (19.741) used by the atomic
	! branch. For every photo-ionizable metal ion: its cross section at those
	! five energies (sm), the ratio to sigma_HI there (Rm), and the
	! photoelectron energy E_c - E_th it leaves (Em). mabs lists the ions
	! that absorb at least one channel, so a run with metals off has
	! n_mabs = 0 and every metal sum below is untouched from 0.0.
	integer, parameter :: n_chan_hrc = 5
	real*8  :: E_chan(n_chan_hrc), sHI_chan(n_chan_hrc)
	real*8  :: sm(n_mion,n_chan_hrc), Rm(n_mion,n_chan_hrc)
	real*8  :: Em(n_mion,n_chan_hrc), gmv(n_mion)
	real*8  :: nmtot, gm, Pm_add, Hm_add, fs_Rm, Ees_m
	real*8  :: tm_gnd, tm_584, tm_19, tm_2q, tm_casc
	real*8  :: um_gnd, um_584, um_19, um_2q, um_casc
	integer :: mabs(n_mion), n_mabs, im, ic, kph, ii

	! Photoionization cross-section ratios sigma_s/sigma_HI at the
	! representative energy of each channel, used to split that channel's
	! photons among the absorbers (Draine Eq. 14.16 generalized to H2). Their
	! (small) kT dependence is neglected. Evaluated once here, not per cell.
	s_gnd   = sigma(E_gnd_HeI ,1.0d0)
	s_584   = sigma(E_584_HeI ,1.0d0)
	s_19    = sigma(E_19_HeI  ,1.0d0)
	s_2q    = sigma(E_2q_HeI  ,1.0d0)
	s_casc  = sigma(E_casc_HeI,1.0d0)
	s2_gnd  = sigma_H2(E_gnd_HeI)
	s2_584  = sigma_H2(E_584_HeI)
	s2_19   = sigma_H2(E_19_HeI)
	s2_2q   = sigma_H2(E_2q_HeI)
	s2_casc = sigma_H2(E_casc_HeI)
	sHe_gnd = sigma_HeI(E_gnd_HeI)
	Rratio  = sHe_gnd/max(s_gnd , 1.0d-99)
	R2_gnd  = s2_gnd /max(s_gnd , 1.0d-99)
	R2_584  = s2_584 /max(s_584 , 1.0d-99)
	R2_19   = s2_19  /max(s_19  , 1.0d-99)
	R2_2q   = s2_2q  /max(s_2q  , 1.0d-99)
	R2_casc = s2_casc/max(s_casc, 1.0d-99)

	! Metal absorbers of the same photons (item P34), evaluated once here.
	! metal_photoion_sigma returns zero below the ion's own threshold, so an
	! ion that cannot absorb a channel carries Rm = 0 for it and drops out of
	! that channel's competition without a test. Em is the photoelectron
	! energy left behind, meaningful only where Rm > 0.
	E_chan   = [E_gnd_HeI, E_584_HeI, E_19_HeI, E_2q_HeI, E_casc_HeI]
	sHI_chan = [s_gnd, s_584, s_19, s_2q, s_casc]
	sm     = 0.0d0
	Rm     = 0.0d0
	Em     = 0.0d0
	gmv    = 0.0d0
	dP_m   = 0.0d0
	n_mabs = 0
	if (thereis_metals) then
		do im = 1,n_mion
			if (.not. mion_isphot(im)) cycle
			kph = mion_iphot(im)
			do ic = 1,n_chan_hrc
				sm(im,ic) = metal_photoion_sigma(kph, E_chan(ic))
				Rm(im,ic) = sm(im,ic)/max(sHI_chan(ic), 1.0d-99)
				Em(im,ic) = E_chan(ic) - mion_ethr(im)
			enddo
			if (maxval(sm(im,:)) .gt. 0.0d0) then
				n_mabs       = n_mabs + 1
				mabs(n_mabs) = im
			endif
		enddo
	endif

	! alpha_B is the active He II case B (whichever fit eval_cool uses).
	call rec_HeII_B(T_K, alphaB)

	if (.not. thereis_HeITR) then
		! atomic (case-B) mode: alpha_1 = Mao & Kaastra ground capture.
		alpha1 = alpha1_HeII_mao(T_K)
		do j = 1-Ng,N+Ng
			n2   = nh2(j)
			nhig = max(nhi(j),1.0d-99)
			! Metal opacity (tm) and metal/H I absorption ratio (um) of the
			! two channels this branch uses, item P34. Both are untouched
			! zeros without metals, so every expression below reduces to the
			! H I / H2 / He I form bit for bit.
			nmtot   = 0.0d0
			tm_gnd  = 0.0d0
			tm_casc = 0.0d0
			um_gnd  = 0.0d0
			um_casc = 0.0d0
			do ii = 1,n_mabs
				im      = mabs(ii)
				gm      = nm(j,im)
				gmv(im) = gm/nhig
				nmtot   = nmtot   + gm
				tm_gnd  = tm_gnd  + gm*sm(im,1)
				tm_casc = tm_casc + gm*sm(im,5)
				um_gnd  = um_gnd  + gmv(im)*Rm(im,1)
				um_casc = um_casc + gmv(im)*Rm(im,5)
			enddo
			absorb_gnd = (nhi(j) + n2 + nhei(j) + nmtot) .gt. 0.0d0
			absorb_low = (nhi(j) + n2 + nmtot)           .gt. 0.0d0
			! Fraction of each channel's photons absorbed within the cell,
			! 1 - exp(-tau_c) over the cell width; the rest leaves. The
			! sigma_* functions return the cross section in 1e-18 cm^2, the
			! unit the stellar optical depth of PH_heat_HHe also carries.
			dl     = dr_j(j)*R0*1.0d-18
			f_gnd  = absorbed_photon_fraction((nhi(j)*s_gnd               &
			         + n2*s2_gnd + nhei(j)*sHe_gnd + tm_gnd)*dl)
			f_casc = absorbed_photon_fraction((nhi(j)*s_casc              &
			         + n2*s2_casc + tm_casc)*dl)
			! y, y2: H I and H2 shares of the >= 24.6 eV ground-capture
			! continuum; the rest is re-absorbed by He I.
			if (absorb_gnd) then
				u2 = R2_gnd*n2/nhig
				y  = 1.0d0/(1.0d0 + Rratio*nhei(j)/nhig + u2 + um_gnd)
				y2 = y*u2
			else
				y  = 0.0d0
				y2 = 0.0d0
			endif
			! wc, wc2: H I and H2 shares of the case-B cascade exits, all of
			! which lie below the He I edge, at the cascade-averaged energy.
			if (absorb_low) then
				u2   = R2_casc*n2/nhig
				wc   = 1.0d0/(1.0d0 + u2 + um_casc)
				wc2  = wc*u2
			else
				wc   = 0.0d0
				wc2  = 0.0d0
			endif
			! z: density-dependent fraction of case-B cascade photons above the
			! H I edge; 0.96 (low density, 19.8 eV line ionizes H) -> 0.67 (2^3S
			! collisionally converted to singlets) via the 2^3S critical density
			! (documented interpolation between Draine's two limits). The
			! fraction above the H2 edge differs slightly (the part of the
			! two-photon continuum between 13.6 and 15.4 eV cannot ionize H2);
			! z is used for both, as the two-photon exit is 8% of the cascade.
			T4    = T_K(j)/1.0d4
			ncrit = 1100.0d0*exp(1.2d0/T4)*sqrt(T4)
			z     = 0.67d0 + 0.29d0/(1.0d0 + ne(j)/ncrit)
			! He II recombination: alpha_eff = alpha_B + alpha_1
			! (1 - f_abs w_HeI); a ground-capture photon that leaves the cell
			! does not re-ionize He, so it is a net recombination too.
			y_net = (1.0d0 - f_gnd) + f_gnd*(y + y2 + y*um_gnd)
			rcheiiB_new(j) = alphaB(j) + y_net*alpha1(j)
			! Extra H I and H2 photoionization rates [s^-1], each channel's
			! share of the photons it keeps. The H2 rate is the H I one channel
			! by channel times sigma_H2/sigma_HI, because both shares carry the
			! same 1/(n_HI sigma_HI + n_H2 sigma_H2) factor.
			dP_HI(j) = nheii(j)*ne(j)*(z*alphaB(j)*wc*f_casc                &
			           + y*alpha1(j)*f_gnd)/nhig
			dP_H2(j) = nheii(j)*ne(j)*(z*alphaB(j)*wc*R2_casc*f_casc        &
			           + y*alpha1(j)*R2_gnd*f_gnd)/nhig
			! The same two channels absorbed by each metal ion: a rate per
			! ion of that stage, added to its photoionization rate by the
			! caller, and its photoelectron heating at E_c - E_th of that ion
			! (item P34). Hm_add stays 0.0 without metals.
			Hm_add = 0.0d0
			do ii = 1,n_mabs
				im = mabs(ii)
				dP_m(j,im) = nheii(j)*ne(j)                                &
				             *(z*alphaB(j)*wc*Rm(im,5)*f_casc              &
				               + y*alpha1(j)*Rm(im,1)*f_gnd)/nhig
				Hm_add = Hm_add + gmv(im)*nheii(j)*ne(j)                   &
				         *(z*alphaB(j)*wc*Rm(im,5)*Em(im,5)*f_casc         &
				           + y*alpha1(j)*Rm(im,1)*Em(im,1)*f_gnd)
			enddo
			! Photoelectron heating [erg cm^-3 s^-1]. Ee_casc_HeI = 6.14 eV is a
			! low-density channel-weighted average, approximate;
			! E_gnd = 11.0 eV (24.6-13.6) for H I, 9.2 eV (24.6-15.4) for H2.
			dheat(j) = (nheii(j)*ne(j)*(z*alphaB(j)*wc*Ee_casc_HeI*f_casc  &
			           + y*alpha1(j)*11.0d0*f_gnd)                          &
			           + nheii(j)*ne(j)*(z*alphaB(j)*wc2*EeH2_casc_HeI     &
			             *f_casc                                            &
			           + y2*alpha1(j)*EeH2_gnd_HeI*f_gnd)                  &
			           + Hm_add)/erg2eV
		enddo
	else
		! TR mode: alpha_1 = 1^1S channel (rec_HeII_11S).
		call rec_HeII_11S(T_K, alpha1)
		do j = 1-Ng,N+Ng
			n2   = nh2(j)
			nhig = max(nhi(j),1.0d-99)
			! Metal opacity and metal/H I absorption ratio of the four
			! channels this branch resolves (item P34); untouched zeros
			! without metals.
			nmtot  = 0.0d0
			tm_gnd = 0.0d0
			tm_584 = 0.0d0
			tm_19  = 0.0d0
			tm_2q  = 0.0d0
			um_gnd = 0.0d0
			um_584 = 0.0d0
			um_19  = 0.0d0
			um_2q  = 0.0d0
			do ii = 1,n_mabs
				im      = mabs(ii)
				gm      = nm(j,im)
				gmv(im) = gm/nhig
				nmtot   = nmtot  + gm
				tm_gnd  = tm_gnd + gm*sm(im,1)
				tm_584  = tm_584 + gm*sm(im,2)
				tm_19   = tm_19  + gm*sm(im,3)
				tm_2q   = tm_2q  + gm*sm(im,4)
				um_gnd  = um_gnd + gmv(im)*Rm(im,1)
				um_584  = um_584 + gmv(im)*Rm(im,2)
				um_19   = um_19  + gmv(im)*Rm(im,3)
				um_2q   = um_2q  + gmv(im)*Rm(im,4)
			enddo
			absorb_gnd = (nhi(j) + n2 + nhei(j) + nmtot) .gt. 0.0d0
			absorb_low = (nhi(j) + n2 + nmtot)           .gt. 0.0d0
			! Fraction of each channel's photons absorbed within the cell
			! (dl carries the 1e-18 cm^2 unit of the sigma_* functions).
			dl    = dr_j(j)*R0*1.0d-18
			f_gnd = absorbed_photon_fraction((nhi(j)*s_gnd                &
			        + n2*s2_gnd + nhei(j)*sHe_gnd + tm_gnd)*dl)
			f_584 = absorbed_photon_fraction((nhi(j)*s_584 + n2*s2_584    &
			        + tm_584)*dl)
			f_19  = absorbed_photon_fraction((nhi(j)*s_19  + n2*s2_19     &
			        + tm_19 )*dl)
			f_2q  = absorbed_photon_fraction((nhi(j)*s_2q  + n2*s2_2q     &
			        + tm_2q )*dl)
			if (absorb_gnd) then
				u2 = R2_gnd*n2/nhig
				y  = 1.0d0/(1.0d0 + Rratio*nhei(j)/nhig + u2 + um_gnd)
				y2 = y*u2
			else
				y  = 0.0d0
				y2 = 0.0d0
			endif
			! H I / H2 shares of the three sub-24.6 eV exits: the 584 A
			! resonance, the 19.8 eV 2^3S line, and the 2^1S two-photon
			! continuum at its mean in-band energy.
			if (absorb_low) then
				u2     = R2_584*n2/nhig
				w584   = 1.0d0/(1.0d0 + u2 + um_584)
				w584_2 = w584*u2
				u2     = R2_19*n2/nhig
				w19    = 1.0d0/(1.0d0 + u2 + um_19)
				w19_2  = w19*u2
				u2     = R2_2q*n2/nhig
				w2q    = 1.0d0/(1.0d0 + u2 + um_2q)
				w2q_2  = w2q*u2
			else
				w584   = 0.0d0
				w584_2 = 0.0d0
				w19    = 0.0d0
				w19_2  = 0.0d0
				w2q    = 0.0d0
				w2q_2  = 0.0d0
			endif
			! Singlet-excited captures: 2/3 go to 2^1P -> 584 A (always
			! ionizing), 1/3 to 2^1S -> two-photon (f_2q_HeI of which lies above
			! the H I edge). Absorber-weighted photon counts and deposited
			! energies per capture; with n_H2 = 0 these reduce exactly to the
			! constants f_sing_HeI and Ee_sing_HeI.
			fs_HI  = (2.0d0*w584*f_584                                      &
			          + f_2q_HeI*w2q*f_2q  )/3.0d0
			fs_R2  = (2.0d0*w584*R2_584*f_584                               &
			          + f_2q_HeI*w2q*R2_2q*f_2q)/3.0d0
			Ees_HI = (2.0d0*Ee_584_HeI*w584*f_584                           &
			          + f_2q_HeI*Ee_2q_HeI*w2q*f_2q)/3.0d0
			Ees_H2 = (2.0d0*EeH2_584_HeI*w584_2*f_584                       &
			          + f_2q_HeI*EeH2_2q_HeI*w2q_2*f_2q)/3.0d0
			! (1) 1^1S channel coefficient: net ground capture, alpha_1
			! (1 - f_abs w_HeI), plus the singlet-excited capture channel
			! (0.25 alpha_B) missing from the current network. Those exits are
			! all below the He I edge, so they are a net recombination whether
			! they are absorbed or not.
			y_net = (1.0d0 - f_gnd) + f_gnd*(y + y2 + y*um_gnd)
			rcheiiB_new(j) = y_net*alpha1(j) + 0.25d0*alphaB(j)
			! (2) H-ionizing photon production [cm^-3 s^-1], each channel
			! reduced to the photons the cell keeps:
			P_add = y*alpha1(j)*nheii(j)*ne(j)*f_gnd            ! ground (>=24.6)
			P_add = P_add + fs_HI*0.25d0*alphaB(j)*nheii(j)*ne(j)
			! 2^3S radiative decay (19.8 eV line).
			P_add = P_add + A31*nheiTR(j)*w19*f_19
			! 2^3S collisionally converted to singlets, then decaying: 2^1S
			! two-photon (f_2q_HeI ionizing) + 2^1P -> 584 A (1.0 ionizing).
			P_add = P_add + ne(j)*nheiTR(j)                                &
			        *(q31a(j)*f_2q_HeI*w2q*f_2q + q31b(j)*w584*f_584)
			dP_HI(j) = P_add/nhig
			! (2b) the same photons absorbed by H2 instead, as a rate per H2
			! molecule (the n_H2 of the share cancels, as for H I above).
			P2_add = y*alpha1(j)*nheii(j)*ne(j)*R2_gnd*f_gnd
			P2_add = P2_add + fs_R2*0.25d0*alphaB(j)*nheii(j)*ne(j)
			P2_add = P2_add + A31*nheiTR(j)*w19*R2_19*f_19
			P2_add = P2_add + ne(j)*nheiTR(j)                              &
			         *(q31a(j)*f_2q_HeI*w2q*R2_2q*f_2q                     &
			           + q31b(j)*w584*R2_584*f_584)
			dP_H2(j) = P2_add/nhig
			! (2c) the same four channels absorbed by each metal ion, as a
			! rate per ion of that stage (item P34). Same structure as (2b):
			! the metal density of the share cancels, leaving Rm/n_HI.
			do ii = 1,n_mabs
				im = mabs(ii)
				fs_Rm = (2.0d0*w584*Rm(im,2)*f_584                        &
				         + f_2q_HeI*w2q*Rm(im,4)*f_2q)/3.0d0
				Pm_add = y*alpha1(j)*nheii(j)*ne(j)*Rm(im,1)*f_gnd
				Pm_add = Pm_add + fs_Rm*0.25d0*alphaB(j)*nheii(j)*ne(j)
				Pm_add = Pm_add + A31*nheiTR(j)*w19*Rm(im,3)*f_19
				Pm_add = Pm_add + ne(j)*nheiTR(j)                          &
				         *(q31a(j)*f_2q_HeI*w2q*Rm(im,4)*f_2q              &
				           + q31b(j)*w584*Rm(im,2)*f_584)
				dP_m(j,im) = Pm_add/nhig
			enddo
			! (3) photoelectron heating [erg cm^-3 s^-1], channel E_dep [eV]:
			! ground 11.0; singlet-excited Ee_sing_HeI = 5.53; 19.8 eV line
			! -> 6.2; 2^3S coll. -> singlet: two-photon Ee_2q_HeI + 584 A 7.6.
			H_add = y*alpha1(j)*nheii(j)*ne(j)*11.0d0*f_gnd
			H_add = H_add + 0.25d0*alphaB(j)*nheii(j)*ne(j)*Ees_HI
			H_add = H_add + A31*nheiTR(j)*w19*6.2d0*f_19
			H_add = H_add + ne(j)*nheiTR(j)                                &
			        *(q31a(j)*f_2q_HeI*Ee_2q_HeI*w2q*f_2q                  &
			          + q31b(j)*Ee_584_HeI*w584*f_584)
			! and the same channels deposited in H2, 15.4 eV lower each.
			H2_add = y2*alpha1(j)*nheii(j)*ne(j)*EeH2_gnd_HeI*f_gnd
			H2_add = H2_add + 0.25d0*alphaB(j)*nheii(j)*ne(j)*Ees_H2
			H2_add = H2_add + A31*nheiTR(j)*w19_2*EeH2_19_HeI*f_19
			H2_add = H2_add + ne(j)*nheiTR(j)                              &
			         *(q31a(j)*f_2q_HeI*EeH2_2q_HeI*w2q_2*f_2q             &
			           + q31b(j)*EeH2_584_HeI*w584_2*f_584)
			! and in each metal ion, E_c - E_th of that ion per photon; the
			! shares carry the metal density, as the H2 ones carry n_H2.
			! Hm_add stays 0.0 without metals (item P34).
			Hm_add = 0.0d0
			do ii = 1,n_mabs
				im    = mabs(ii)
				Ees_m = (2.0d0*Em(im,2)*w584*Rm(im,2)*gmv(im)*f_584        &
				         + f_2q_HeI*Em(im,4)*w2q*Rm(im,4)*gmv(im)*f_2q)    &
				        /3.0d0
				Hm_add = Hm_add                                            &
				         + y*Rm(im,1)*gmv(im)*alpha1(j)*nheii(j)*ne(j)     &
				           *Em(im,1)*f_gnd                                 &
				         + 0.25d0*alphaB(j)*nheii(j)*ne(j)*Ees_m           &
				         + A31*nheiTR(j)*w19*Rm(im,3)*gmv(im)*Em(im,3)     &
				           *f_19                                           &
				         + ne(j)*nheiTR(j)                                 &
				           *(q31a(j)*f_2q_HeI*Em(im,4)*w2q*Rm(im,4)        &
				             *gmv(im)*f_2q                                 &
				             + q31b(j)*Em(im,2)*w584*Rm(im,2)*gmv(im)      &
				               *f_584)
			enddo
			dheat(j) = (H_add + H2_add + Hm_add)/erg2eV
		enddo
	endif

	! End of subroutine
	end subroutine he_rec_coupling


	! End of module
	end module utils_ion_eq

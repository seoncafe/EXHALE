   module utils_ion_eq

   use global_parameters
   use species_table, only: n_mion, n_mphot, n_melem,                 &
                            mion_isphot, mion_iphot, mion_ethr,       &
                            mion_z2, mion_elem, mion_iscool,          &
                            mion_stage, mion_name, mion_fsp,          &
                            melem_Z, melem_top, im_FeII, im_OI,       &
                            isp_HeTR,                                    &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,          &
                            isp_OH, isp_H2O, isp_CO
   use utils
   use Cooling_Coefficients      ! Various functions for cooling coefficients
   ! Miller+2013 H3+ infrared cooling, emitted and net of the lower-atmosphere field
   use h3p_cooling, only: h3p_cooling_rate, h3p_net_cooling_rate
   ! H2 quadrupole/magnetic-dipole lines and the H2O and CO bands, net of the
   ! infrared field of the lower atmosphere (`Molecular IR bands`)
   use molecular_reaction_heat, only: molecular_chemical_heating,      &
                            oxygen_chemical_heating,                    &
                            oxygen_reaction_energy_eV,                  &
                            species_formation_energy, ir_assoc_HeTR, &
                            ir_D1
   ! Quantum yield of the H2O + hv -> H2 + O(1D) branch, band by band: the
   ! O(1D) production that its local steady state carries into the O6 sink.
   use oxygen_rates, only: qy_H2O_H2_O1D, rk_D1_Hep_CO
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
   ! Where the energy of ONE H2 photoevent goes, channel by channel:
   ! h2_channel_energy_recipients (h2_photo_channels.f90) is the single
   ! definition, and the two scalars below are the entries of it that the
   ! heating integrand needs.  h2_double_fragment_kinetic_energy is the
   ! kinetic energy release of the two protons of the double channel, fixed
   ! with no free parameter by the VERTICAL 51.4 eV threshold of Yan,
   ! Sadeghpour & Dalgarno (1998): at that photon energy the two electrons
   ! come off at rest and the asymptotic product energy is D0 + 2 I(H) =
   ! 31.675 eV, so 19.725 eV is on the receding nuclei.  e_rad_H2_neutral is
   ! the internal excitation of the H(2p) + H(2s) pair the neutral window
   ! leaves (Chung, Lee, Masuoka & Samson 1993), which departs as prompt
   ! Ly-alpha and two-photon continuum and is not heat.
   use h2_photo_channels, only: h2_double_fragment_kinetic_energy,        &
                                e_rad_H2_neutral
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
   use water_photolysis, only: n_fuv_band, ib_LW, ib_B2,             &
                       ib_B3, ib_B4,                                    &
                       water_photolysis_rate, hydroxyl_photolysis_rate, &
                       fuv_band_optical_depth,                          &
                       heat_per_water_dissociation,                     &
                       heat_per_hydroxyl_dissociation,                  &
                       sigma_H2O_band, e_photon_flat_band
   ! H I Ly-alpha: the line-centre optical depth of the atomic hydrogen
   ! column and the MEAN over one cell of the stellar beam's transmission
   ! through it.  Band B2 of the FUV photolysis set IS that resonance line,
   ! so it reads the same two quantities the n = 2 pumping field is built
   ! from.
   use lya_rt, only: lya_line_center_optical_depth,                    &
                     lya_stellar_beam_transmission_cell_mean
   use lyman_werner_photodissociation, only:                          &
                            lyman_werner_dissociation_rate_cell_mean,   &
                            lyman_werner_band_absorption_rate_cell_mean, &
                            e_lw_fragment_erg, e_lw_photon_erg,          &
                            h2_lw_dissociation_per_pump,                &
                            h2_lw_dissociation_per_absorbed_photon,     &
                            h2_self_shielding_level_resolved,           &
                            h2_shield_max_column,                       &
                            h2_lw_band_photon_fraction_absorbed
   ! CO on the same beam: the shielded photodissociation rate and its cell
   ! mean, the Visser shielding function itself for the reported column, the
   ! kinetic energy of one dissociation event, and the He+ charge-transfer
   ! rate coefficient that is the other CO destruction channel.
   use co_photodissociation, only:                                     &
                            co_photodissociation_rate_cell_mean,        &
                            heat_per_co_dissociation, e_co_photon_erg
   use co_self_shielding_table, only: co_self_shielding
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
	! The ground-capture continuum starts AT the He I ionization threshold,
	! so this energy is that threshold and is read from it.
	real*8, parameter :: E_gnd_HeI  = e_th_HeI ! ground-capture continuum edge
	real*8, parameter :: E_584_HeI  = 21.2d0   ! 2^1P -> 1^1S resonance line
	real*8, parameter :: E_19_HeI   = 19.8d0   ! 2^3S -> 1^1S line
	! Mean energy of the 2^1S two-photon photons that lie above the H I edge:
	! the continuum is not a line, and this single energy stands for it. The
	! shape (Drake, Victor & Dalgarno 1969) is not carried in the code, only
	! the two integrals f_2q_HeI and Ee_2q_HeI over the H I window, so the
	! energy below is e_th_HI + Ee_2q_HeI and cannot be re-integrated over the
	! narrower H2 window; see EeH2_2q_HeI.
	real*8, parameter :: E_2q_HeI   = e_th_HI + Ee_2q_HeI       ! 16.110 eV
	! Cascade-averaged exit energy used by the atomic (case-B) branch, the same
	! 0.75/0.17/0.08 weighting that defines Ee_casc_HeI.
	real*8, parameter :: E_casc_HeI = e_th_HI + Ee_casc_HeI     ! 19.741 eV
	! The photoelectron each of those photons leaves in H2 (threshold
	! e_th_H2), the H2 counterparts of 11.0 / 7.6 / 6.2 / Ee_2q_HeI /
	! Ee_casc_HeI.
	real*8, parameter :: EeH2_gnd_HeI  = E_gnd_HeI  - e_th_H2   ! 9.161 eV
	real*8, parameter :: EeH2_584_HeI  = E_584_HeI  - e_th_H2   ! 5.774 eV
	real*8, parameter :: EeH2_19_HeI   = E_19_HeI   - e_th_H2   ! 4.374 eV
	! APPROXIMATE: the H2-ionizing photon COUNT of the two-photon continuum is
	! the integral of the same shape over e_th_H2 to 20.62 eV, which is smaller
	! than f_2q_HeI (the integral over e_th_HI to 20.62 eV); the shape is not in the
	! code, so f_2q_HeI is used for H2 as well and this sub-channel's H2 share
	! is overestimated. It is 1/3 of one of five channels, and R2_2q = 1.2 is
	! the smallest H2/H I cross-section ratio of the set.
	real*8, parameter :: EeH2_2q_HeI   = E_2q_HeI   - e_th_H2   ! 0.685 eV
	real*8, parameter :: EeH2_casc_HeI = E_casc_HeI - e_th_H2   ! 4.315 eV

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

	! ----- Channels of the volumetric heating rate -----
	! The heating of a cell is a sum of physically distinct deposits. There
	! is ONE list of them: heating_of_composition fills every column of
	! heat_chan, the total it returns is the sum of those columns in this
	! order, and the header of output/Heating_breakdown.txt is built from the
	! names below. A deposit therefore cannot exist in the energy equation
	! and be missing from the breakdown, which is what two hand-maintained
	! copies of the sum allowed.
	!
	!   1-6  photoionization of each absorber (photoheating_of_composition)
	!   7    photoelectric heating of H(n=2)
	!   8    collisional de-excitation of H(n=2) by the Lyman-alpha field
	!   9    photoelectrons of the H I / H2 / metal ionizations driven by He
	!        recombination radiation
	!  10    He(2^3S) + H  -> He + H+ + e   (Penning branch)
	!  11    He(2^3S) + H  -> HeH+ + e      (associative branch)
	!  12    He(2^3S) + H2 -> He + H2+ + e  (Penning branch)
	!  13    kinetic energy of the H2 Lyman-Werner dissociation fragments
	!  14    collisional de-excitation of the H2 Lyman-Werner fluorescence
	!  15    collisional reactions of the H2/He network
	!  16    excess energy of the H2O and OH FUV photolysis events
	!  17    collisional reactions of the oxygen network, O(1D) sink included
	!  18    He+ + CO -> C+ + O + He, the charge-transfer destruction of CO
	!  19    kinetic energy of the CO photodissociation fragments
	integer, parameter :: n_heat_channel = 19
	! Named index of the one channel a consumer outside this module reads
	! on its own: output/Lyman_Werner.txt writes the Lyman-Werner fragment
	! deposit beside the rate it belongs to, and takes it from the array
	! the sweep filled rather than re-forming it there.
	integer, parameter :: ih_H2_LW_dissoc = 13
	character(len=24), parameter ::                                        &
	   heat_channel_name(n_heat_channel) = (/                              &
	      'heat_HI                 ', 'heat_HeI                ',          &
	      'heat_HeII               ', 'heat_He23S              ',          &
	      'heat_H2                 ', 'heat_metals             ',          &
	      'heat_Hpe[excitedH]      ', 'heat_Hdx[Lya-deexc]     ',          &
	      'heat_He_recomb          ', 'heat_He23S_Penning      ',          &
	      'heat_He23S_assoc        ', 'heat_He23S_H2_Penning   ',          &
	      'heat_H2_LW_dissoc       ', 'heat_H2_LW_fluor        ',          &
	      'heat_mol_chem           ', 'heat_FUV_photolysis     ',          &
	      'heat_oxygen_collisional ', 'heat_CO_Hep_transfer    ',          &
	      'heat_CO_photodissoc     ' /)

	! The channel array of the state the ionization sweep RETURNED. The
	! sweep passes this array to heating_of_composition, so it holds the
	! deposits the heat column of Hydro_ioniz.txt was built from, at the
	! composition and with the lagged rates of that same sweep.
	! Heating_breakdown.txt writes this array and recomputes nothing: one
	! state per output file, so the channel sum and the heat column are the
	! same number to round-off rather than to the size of one sweep's rate
	! lag. Allocated with the other grid-sized arrays of the sweep.
	real*8, allocatable :: heat_channel_state(:,:)

	! Length of the cell blocks the two parallel sweeps of this module
	! (chemical_rate_coefficients, eval_cool) are cut into. EVEN, and not a
	! function of the thread count: the vectorized rate loops call libmvec's
	! two-lane exp/log/pow on cell pairs counted from each block's first
	! cell, and the pair variant differs from the scalar remainder in the
	! last bit, so a block start at an odd offset, or one that moved with
	! the thread count, moved the last bit of the result. 32 keeps 16 blocks
	! on the 504-cell grids of the regression matrix.
	integer, parameter :: xuv_rate_block = 32

	contains

	! Incident stellar flux of FUV band ib at the planet [erg cm^-2 s^-1],
	! in the band order of water_photolysis (LW, B2 = Ly-alpha, B3, B4).
	! Single definition, read by the equilibrium solve and by write_output.
	!
	! THE FIRST BAND IS THE LYMAN-WERNER INTERVAL and carries the flux of the
	! "Stellar LW flux" key, because 912-1201 A is one interval with one
	! incident flux whether the absorber is H2 in lines or H2O and OH in a
	! continuum. Supplying it twice -- once for H2, once inside a band
	! reaching across it -- would count the energy of the interval twice, and
	! that is what B2 starting at 1202 A removes.
	!
	! B2 IS THE INCIDENT STELLAR LY-ALPHA FLUX AT THE PLANET'S ORBIT, and
	! this function returns it undepleted: the attenuation of that beam on
	! its way down belongs to the field routine, not to the flux. Two
	! absorbers deplete it there. The H2O and OH continua, as in every band,
	! and -- because B2 is the H I resonance line itself -- the atomic
	! hydrogen column, through the fraction of the broad stellar line whose
	! wings penetrate to a given line-centre depth
	! (lya_stellar_beam_transmission in lya_rt.f90, the same expression that
	! module's own stellar term uses for the n = 2 pumping, and averaged
	! over a cell by the same routine that module averages it with). That
	! is a transmission of the beam and not the trapped mean intensity, so
	! the photons the band's ledger removes are still the photons its
	! absorbers take; the trapping buildup lya_rt applies on top of it
	! raises the LOCAL mean intensity of the scattered field and would
	! double-count photons already removed if it were applied here.
	!
	! Section 2.6 of docs/a2_oxygen_option_design.md asks instead for the
	! solved J_Lya(r) of lya_rt as B2's field. That remains the fuller
	! treatment: it would also carry the internally generated Ly-alpha of
	! the recombination cascade, which a stellar transmission cannot. The
	! band is reported separately in output/FUV_bands.txt so that a run
	! whose oxygen chemistry turns on B2 can be read off.
	!
	! Note that setting "Stellar Lya flux" for the n = 2 pumping therefore
	! also drives B2, and setting "Stellar LW flux" for the H2 network also
	! drives the oxygen photolysis of the LW band; the setup report states
	! all four band fluxes at startup so this is visible rather than implicit.
	! The dayside dilution is applied HERE and not at the keys, so that
	! F_LW_star and the rest keep their stated meaning -- the band flux at
	! the planet's orbit -- in the setup report and the resolved dump, and
	! so that a flux integrated from the spectrum file is diluted by the
	! same convention as a stated one.  Before section 150 these bands
	! carried no dilution at all while the same run halved its XUV and its
	! stellar Ly-alpha beam, so the Lyman-Werner rate stood a factor
	! 1/dayside_dilution() above the run's own convention.
	double precision function fuv_band_flux(ib) result(F)
	integer, intent(in) :: ib
	select case (ib)
	case (ib_LW)
		F = F_LW_star
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
	! WHY ONE ROUTINE.  Over 912-1201 A the three absorbers share one beam.
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
	subroutine fuv_lw_photon_field(nH2, nH2O, nOH, nCO, nHI, T_K, nH_nuc, &
	                               NH2col, NH2Ocol, NOHcol, NCOcol,       &
	                               f_shield, tr_lines, tau_b,             &
	                               k_lw, p_lw_single, p_lw_absorbed,      &
	                               k_co, theta_co,                        &
	                               j_h2o, j_oh, a_lines_max,              &
	                               col_over_overlap)
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nH2, nH2O, nOH, T_K
	! Carbon monoxide [cm^-3].  It is the fourth absorber of this beam: CO
	! predissociates in 37 lines between 912.7 and 1076.1 A, all inside the
	! Lyman-Werner interval, and it shields itself in them
	! (co_photodissociation.f90).  Zero for a run without the oxygen
	! chemistry.
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nCO
	! Neutral atomic hydrogen [cm^-3].  It is an absorber of this beam and
	! not a spectator: band B2 is the H I Ly-alpha resonance line, and the
	! column above a molecular base decides how much of the stellar line
	! reaches it.
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nHI
	! Total hydrogen NUCLEUS density, the density axis of the self-shielding
	! table (CLOUDY's hden).  It enters nothing else here.
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nH_nuc
	real*8, dimension(1-Ng:N+Ng), intent(out) :: NH2col, NH2Ocol, NOHcol
	! Star-ward CO column [cm^-2], on the same radial points and by the same
	! rectangle rule as the other three, so NCOcol(j) is the column at the
	! inner face of cell j and NCOcol(j+1) that at its star-ward face.  That
	! face convention is what the cell mean of the CO rate depends on.
	real*8, dimension(1-Ng:N+Ng), intent(out) :: NCOcol
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
	! CO photodissociation rate [s^-1], the mean over each cell, and the
	! Visser shielding function at the cell's own two columns, which is what
	! the run reports.  Both 0 (rate) and 1 (shielding, the neutral value)
	! for a run with no CO or no band flux.
	real*8, dimension(1-Ng:N+Ng), intent(out) :: k_co, theta_co
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(out) :: tau_b
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(out) :: j_h2o, j_oh
	real*8, intent(out) :: a_lines_max
	! Largest ratio of the star-ward H2 column to the column at which the
	! self-shielding table stops being a value and becomes an upper bound
	! (lyman_werner.f90 sec. 2d).  Above 1 the deepest cells are in the
	! line-overlap regime the table does not carry.
	real*8, intent(out) :: col_over_overlap
	real*8, dimension(1-Ng:N+Ng) :: a_lines
	! Line-centre H I Ly-alpha optical depth for band B2: the depth at each
	! cell's INNER face, the depth at its star-ward face, and the cell's own
	! depth, so that the beam can be averaged across the cell it crosses.
	real*8, dimension(1-Ng:N+Ng) :: tau_lya, tau_lya_out, dtau_lya
	! tr_out2 receives the cell mean of the SQUARE of the beam transmission,
	! which the same routine returns for the quadratic trapping term of the
	! n = 2 pumping field.  The band is linear in the beam, so it is not
	! read here.
	real*8  :: tau_out, dtau, tr_out, tr_out2
	! Star-ward face values of the H2 column and of the Lyman-Werner
	! continuum depth, the outer end of the cell the rate is averaged over.
	real*8  :: NH2col_out, tau_lw_out
	! Star-ward face value of the CO column, the outer end of the cell the
	! CO rate is averaged over.
	real*8  :: NCOcol_out
	integer :: j, ib

	NH2col   = 0.0d0
	NH2Ocol  = 0.0d0
	NOHcol   = 0.0d0
	NCOcol   = 0.0d0
	k_co     = 0.0d0
	theta_co = 1.0d0
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
	tau_lya  = 0.0d0
	tau_lya_out = 0.0d0
	dtau_lya = 0.0d0

	! ---- H2 lines: self-shielding and the fraction of the LW band they
	! take out of the shared beam.  Both are read from the table at the H2
	! column above the cell and the cell's own T and n_H (the table is built
	! for a homogeneous column, so the temperature and density of the gas
	! above are approximated by the local ones in both), and both exist
	! whether or not the oxygen chemistry is on.
	!
	! A is the fraction of the band photons the pumping lines have taken
	! out of the beam by that column, int_0^N sigma_pump dN' with
	! sigma_pump = sigma_diss/p_eff, served in closed form from the same
	! level-resolved table the dissociation rate comes from
	! (h2_lw_band_photon_fraction_absorbed).  Taking it from the table
	! rather than re-integrating on the radial grid makes the line
	! transmission independent of the grid and makes the photons the beam
	! loses exactly the photons the rate spends: one absorption, one
	! normalization.  A is a share of the band and not an optical depth, so
	! the transmission of the beam is 1 - A, not exp(-A).  The clamp is a
	! floor on that transmission: the table is clamped at its edge cross
	! section above the top of its column axis, so A goes on growing
	! linearly there, and a capped A means the lines have taken the band.
	!
	! THE TWO LINES BELOW COME FROM ONE TABLE and are the two faces of one
	! absorption: f_shield is the suppression of the photodissociation RATE,
	! which the run reports in output/Lyman_Werner.txt at each cell's own
	! column, and A is the SHARE of the band the same lines remove.  Until
	! 2026-09-06 A was the Draine & Bertoldi (1996) eq. (39) equivalent
	! width of a narrower band, and the two disagreed by 45 per cent
	! (lyman_werner.f90 sec. 3a).  The rate below is the mean of the same
	! table over the cell.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		call calc_column_dens_one(nH2, NH2col)
		do j = 1-Ng,N+Ng
			f_shield(j) = h2_self_shielding_level_resolved(NH2col(j),     &
			                            T_K(j), nH_nuc(j))
			a_lines(j) = min(h2_lw_band_photon_fraction_absorbed(         &
			                 NH2col(j), T_K(j), nH_nuc(j)), 1.0d0)
			tr_lines(j) = 1.0d0 - a_lines(j)
			col_over_overlap = max(col_over_overlap,                       &
			                       NH2col(j)/h2_shield_max_column())
		enddo
		a_lines_max = maxval(a_lines)
	endif

	! ---- H I Ly-alpha, the line that band B2 IS.  The photolysis bands are
	! all longward of the 912 A Lyman edge, so the H I photoionization
	! CONTINUUM does not touch them; the resonance LINE at 1215.67 A does,
	! and it is the strongest absorber in the atmosphere at that wavelength
	! (a star-ward column of 1e19-1e21 cm^-2 gives a line-centre optical
	! depth of 1e6-1e8).  The stellar line is scattered rather than
	! destroyed, so its transmission is not exp(-tau) but the fraction of
	! the stellar profile whose wings penetrate to that depth --
	! lya_stellar_beam_transmission, the same expression and the same
	! line-centre depth that lya_rt.f90 uses for the n = 2 pumping beam, so
	! the two treatments of the one stellar line agree by construction.
	!
	! The depth is the rectangle rule of lya_line_center_optical_depth, the
	! same rule and the same widths as the continuum columns above, so
	! tau_lya_out(j) and tau_b(j+1,ib) are the depths at one face and
	! dtau_lya(j) and tau_b(j,ib) - tau_b(j+1,ib) are two depths of one
	! cell.
	if (thereis_oxychem .and. F_Lya_star .gt. 0.0d0)                       &
		call lya_line_center_optical_depth(T_K, nHI, tau_lya,              &
		                                   dtau_lya, tau_lya_out)

	! ---- H2O and OH continua, on every band.
	if (thereis_oxychem) then
		call calc_column_dens_one(nH2O, NH2Ocol)
		call calc_column_dens_one(nOH,  NOHcol)
		! The CO column is built by the same routine on the same points, but
		! it is NOT added to tau_b: CO is a LINE absorber and its equivalent
		! width is already inside the Visser shielding function, so a
		! continuum term for it would shield the H2O and OH continua with
		! photons the shielding function has already removed.  Its lines
		! shield CO, and the H2 lines shield CO through N_H2 inside Theta.
		call calc_column_dens_one(nCO,  NCOcol)
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
	!
	! The rate is the MEAN over the cell, as the H2O and OH rates of the
	! same beam are.  NH2col(j) already contains the whole of cell j, so it
	! is the column at that cell's INNER face and NH2col(j+1) is the column
	! at its star-ward face (nothing sits above the outermost cell); the
	! continuum depths tau_b(j,ib_LW) and tau_b(j+1,ib_LW) are the two
	! faces of the same cell.  The cross section falls by four decades
	! across the self-shielding transition and fastest at the H2 front,
	! where one cell spans a large fraction of a decade of column, so the
	! inner-face value applied to the whole cell is low there, in one
	! direction, by whatever that fraction is worth
	! (lyman_werner_dissociation_rate_cell_mean states how the mean is
	! taken and to what accuracy).
	!
	! The two branching ratios stay at the cell's own column: they are
	! ratios of one absorption's outcomes and not rates, so they enter the
	! heat return and the band photon ledger as multipliers of a rate that
	! is now the cell mean, and their variation with column is a factor 12
	! (p_single) and 25 (p_absorbed) over the WHOLE column axis against the
	! four decades of the cross section.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		do j = 1-Ng,N+Ng
			if (j .lt. N+Ng) then
				NH2col_out  = NH2col(j+1)
				tau_lw_out  = tau_b(j+1,ib_LW)
			else
				NH2col_out  = 0.0d0
				tau_lw_out  = 0.0d0
			endif
			k_lw(j) = lyman_werner_dissociation_rate_cell_mean(           &
			              fuv_band_flux(ib_LW),                           &
			              NH2col_out, NH2col(j), T_K(j), nH_nuc(j),       &
			              tau_lw_out, tau_b(j,ib_LW))
			p_lw_single(j)   = h2_lw_dissociation_per_pump(NH2col(j),     &
			                       T_K(j), nH_nuc(j))
			p_lw_absorbed(j) = h2_lw_dissociation_per_absorbed_photon(    &
			                       NH2col(j), T_K(j), nH_nuc(j))
		enddo
	endif

	! ---- CO photodissociation on the same beam.  Same faces, same
	! continuum depth and the same reason for a cell mean as the H2 rate
	! above: the Visser shielding function falls by more than three decades
	! along the CO column axis and by more than six along the H2 axis, and
	! one cell at the CO front carries a large fraction of a decade of both.
	!
	! Theta is a function of TWO columns, so both are handed over at both
	! faces.  The reported theta_co is the value at the cell's own inner
	! face -- the same convention f_shield follows for H2 -- and it is a
	! diagnostic; the rate carries the mean.
	if (thereis_oxychem .and. F_LW_star .gt. 0.0d0) then
		do j = 1-Ng,N+Ng
			if (j .lt. N+Ng) then
				NH2col_out  = NH2col(j+1)
				NCOcol_out  = NCOcol(j+1)
				tau_lw_out  = tau_b(j+1,ib_LW)
			else
				NH2col_out  = 0.0d0
				NCOcol_out  = 0.0d0
				tau_lw_out  = 0.0d0
			endif
			k_co(j) = co_photodissociation_rate_cell_mean(                &
			              fuv_band_flux(ib_LW),                           &
			              NCOcol_out, NCOcol(j), NH2col_out, NH2col(j),   &
			              tau_lw_out, tau_b(j,ib_LW))
			theta_co(j) = co_self_shielding(NCOcol(j), NH2col(j))
		enddo
	endif

	! ---- H2O and OH photodissociation, band by band.  tau_b(j,ib) is the
	! depth at the INNER face of cell j (the column at j already contains
	! cell j), so the star-ward face of cell j carries tau_b(j+1,ib) and the
	! H2 line transmission tr_lines(j+1).  The rate is the MEAN over the
	! cell, which makes the photons the model absorbs exactly the photons
	! the beam loses at any grid spacing (water_photolysis.f90).  The B2
	! line factor below is a cell mean as well, so in that band both factors
	! of the rate are means and neither is a face value.
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
				! Line transmission of the beam, band by band.  tr_lines
				! carries the H2 Lyman-Werner lines, which exist only in
				! the LW band; band B2 is the H I Ly-alpha line and carries
				! the transmission of the stellar line through the atomic
				! hydrogen column; the remaining three bands have no line
				! absorber and keep the neutral value 1.
				!
				! B2 TAKES THE MEAN OF THAT TRANSMISSION OVER THE CELL, not
				! its value at the cell's star-ward face, for the reason the
				! continuum factor is already a cell mean: the rate is what
				! the cell's molecules undergo averaged over the cell, and
				! the beam falls across the cell by that cell's own depth.
				! The stellar Ly-alpha line is resonantly SCATTERED, so its
				! transmission is erfc(c sqrt(tau)) and the fall across a
				! cell is set by the change of the erfc ARGUMENT, i.e. by
				! c[sqrt(tau_in) - sqrt(tau_out)] and not by dtau: it is
				! negligible for a thin cell high in the column and a factor
				! at the base of a molecular layer, where one cell holds
				! most of the atomic hydrogen column
				! (lya_stellar_beam_transmission_cell_mean states the
				! quadrature; the driver lya_beam_cell_mean measures both).
				!
				! THE PRODUCT OF THE TWO MEANS IS TAKEN, not the mean of the
				! product: water_photolysis_rate multiplies this factor by
				! the cell mean of exp(-tau_cont).  Both factors fall
				! monotonically across the cell, so the product of the means
				! is a LOWER bound on the mean of the product, short of it
				! by their covariance, which is at most a quarter of the
				! product of their two relative variations across the cell
				! -- dtau_cont for the exponential factor, and
				! (T_out - T_in)/<T> for this one.  The continuum depth of
				! ONE cell is a fraction of unity wherever the band still
				! carries photons, so that bound is small; the driver
				! lya_band_transmission measures the product of means
				! against an exact joint integration over the cell.  Taking
				! the mean of the product would need one joint quadrature of
				! both absorbers inside water_photolysis_rate.
				!
				! Above the top cell there is no column, tau_lya_out is zero
				! there, and erfc(0) = 1 is the neutral value.
				if (ib .eq. ib_B2) then
					call lya_stellar_beam_transmission_cell_mean(        &
					         T_K(j), tau_lya_out(j), dtau_lya(j),        &
					         tr_out, tr_out2)
				else if (ib .ne. ib_LW) then
					tr_out = 1.0d0
				endif
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

	subroutine advance_starward_columns(j_hi, j_lo,                       &
	              nhi, nheiS, nheii, nheiTR, nh2, nm, has_h2,             &
	              N1_c, N15_c, N2_c, NTR_c, NH2_c, Nm_c,                  &
	              N1_face, N15_face, N2_face, NTR_face, NH2_face, Nm_face)
	! Carries the column of each absorber inward across the cells j_hi down
	! to j_lo, by the same rectangle rule and the same opa_pf opacity weight
	! calc_column_dens, calc_column_dens_one and calc_column_dens_metals use,
	! and in the same order: a column entering at j_hi and leaving at j_lo
	! reproduces those routines cell for cell, bit for bit.
	!
	! N*_c enter as the column of everything OUTSIDE j_hi and leave as the
	! column of everything outside j_lo, the cells of the range included.
	! The optional N*_face arrays record, for each cell of the range, the
	! column outside THAT cell, which is the depth its star-ward face sees
	! (photoionization_field_at_cell_HHe).
	!
	! With the metals off no metal ion is present, so every metal slot of
	! Nm_c would take a zero increment: the sweep over the 27 slots is
	! skipped on thereis_metals and the column each caller reads is the
	! same number it was.
	!
	! Why this exists as a sweep rather than a whole-grid integral: it lets
	! a caller advance the columns with the composition it has already
	! solved, one block of cells at a time, instead of building them once
	! from a composition the whole traversal then lags behind.

	integer, intent(in) :: j_hi, j_lo
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi, nheiS, nheii, nheiTR
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nh2
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm
	logical, intent(in) :: has_h2
	real*8, intent(inout) :: N1_c, N15_c, N2_c, NTR_c, NH2_c
	real*8, dimension(n_mphot), intent(inout) :: Nm_c
	real*8, dimension(1-Ng:N+Ng), intent(out), optional ::                &
	              N1_face, N15_face, N2_face, NTR_face, NH2_face
	real*8, dimension(1-Ng:N+Ng,n_mphot), intent(out), optional :: Nm_face

	integer :: j, i, k
	real*8  :: dr

	do j = j_hi, j_lo, -1

		if (present(N1_face))  N1_face(j)  = N1_c
		if (present(N15_face)) N15_face(j) = N15_c
		if (present(N2_face))  N2_face(j)  = N2_c
		if (present(NTR_face)) NTR_face(j) = NTR_c
		if (present(NH2_face)) NH2_face(j) = NH2_c
		if (present(Nm_face))  Nm_face(j,:) = Nm_c

		dr = dr_j(j)*R0*opa_pf(j)

		! The outermost cell of the grid is the one calc_column_dens starts
		! its integral at, and it groups the width and the density the other
		! way round; the column entering there is zero, so adding its term in
		! that grouping reproduces that first value exactly.
		if (j .eq. N+Ng) then
			N1_c = N1_c + dr_j(j)*R0*nhi(j)*opa_pf(j)
			if (thereis_He) then
				N15_c = N15_c + dr_j(j)*R0*nheiS(j)*opa_pf(j)
				N2_c  = N2_c  + dr_j(j)*R0*nheii(j)*opa_pf(j)
				if (thereis_HeITR) NTR_c = NTR_c                          &
				                 + dr_j(j)*R0*nheiTR(j)*opa_pf(j)
			endif
			if (has_h2) NH2_c = NH2_c + dr_j(j)*R0*nh2(j)*opa_pf(j)
			if (thereis_metals) then
				do i = 1,n_mion
					if (.not. mion_isphot(i)) cycle
					k = mion_iphot(i)
					Nm_c(k) = Nm_c(k) + dr_j(j)*R0*nm(j,i)*opa_pf(j)
				enddo
			endif
		else
			N1_c = N1_c + nhi(j)*dr
			if (thereis_He) then
				N15_c = N15_c + nheiS(j)*dr
				N2_c  = N2_c  + nheii(j)*dr
				if (thereis_HeITR) NTR_c = NTR_c + nheiTR(j)*dr
			endif
			! calc_column_dens_one multiplies the width in cell by cell
			! rather than hoisting it; the same grouping is kept here.
			if (has_h2) NH2_c = NH2_c                                     &
			                  + nh2(j)*dr_j(j)*R0*opa_pf(j)
			if (thereis_metals) then
				do i = 1,n_mion
					if (.not. mion_isphot(i)) cycle
					k = mion_iphot(i)
					Nm_c(k) = Nm_c(k) + nm(j,i)*dr
				enddo
			endif
		endif

	enddo

	end subroutine advance_starward_columns

	! ------------------------------------------------------------- !

	subroutine photoionization_field_at_cell_H(j, N1_out, nhi_j, xion_j,  &
	                          sec_on, P_HI_j, h1_HI_j, heat_j, q_j)
	! The attenuated stellar XUV field of ONE cell of a PURE HYDROGEN
	! atmosphere (thereis_He = .false.), and the rates it drives there.  The
	! H+He+metals form is photoionization_field_at_cell_HHe; PH_heat_H is
	! this routine's whole-grid form.
	!
	! The field at this cell is set by the H I column OUTSIDE it, which the
	! caller passes in, and by the cell's own column, formed here from the
	! density passed in.  Nothing in it depends on any cell inside j, so a
	! caller sweeping outside in can hand it the column of the composition it
	! has already solved and get a field self-consistent with that
	! composition in one traversal; and a caller iterating one cell can hand
	! it the density of the ITERATE, which makes the cell's own attenuation
	! the attenuation of the composition it is being solved for
	! (ionization_equilibrium, xuv_self_field_passes).  The cell mean is
	! formed ONCE, from the two depths together, so an iterating caller
	! re-evaluates that one expression and never composes a second
	! attenuation on top of it.  See photoionization_field_at_cell_HHe for
	! why the fixed point is the same either way.

	integer, intent(in) :: j
	! H I column outside the cell [cm^-2], opa_pf-weighted.
	real*8,  intent(in) :: N1_out
	real*8,  intent(in) :: nhi_j, xion_j
	! Is the secondary-ionization coupling enabled AND staged on.
	logical, intent(in) :: sec_on
	! Photoionization rate [s^-1], photoheating rate of ONE H I atom
	! [erg s^-1], the heating of the composition passed in [erg cm^-3 s^-1]
	! and its heating efficiency.
	real*8,  intent(out) :: P_HI_j, h1_HI_j, heat_j, q_j

	real*8, dimension(Nl) :: acc_secHI, fhv, fh_v, fiHI
	real*8, dimension(Nl) :: fiHeI_dum, fiH2_dum
	type(photoelectron_partition_t) :: pep
	real*8 :: Psec_HI
	real*8 :: PIR_1, Hea_1, q_abs
	real*8, dimension(Nl) :: tauE_out, dtauE
	real*8 :: dr_cell
	real*8, dimension(Nl) :: int_f,int_1,int_q,int_H

	!----------------------------------!

		! Initialization of integrands
		Hea_1  = 0.0
		PIR_1  = 0.0
		q_abs  = 0.0

		! The depth at the star-ward face comes from the column of
		! everything OUTSIDE the cell, which the caller passes in; nothing
		! sits above the outermost cell, so its caller passes zero and the
		! depth there is exactly zero.  The cell's own depth is formed from
		! its own density and width rather than as a difference of two
		! columns, which is the same number in exact arithmetic and loses
		! its digits to cancellation for a thin cell.
		tauE_out = s_hi*N1_out*1.0d-18
		dr_cell = dr_j(j)*R0*opa_pf(j)
		dtauE = s_hi*nhi_j*dr_cell*1.0d-18

		! Initial integrands. int_f is the attenuated stellar flux of this
		! cell weighted by the cell's cross-section factor opa_pf(j) = f(p):
		! the pressure broadening of opacity model 'P' (opacity_pT_factor;
		! f = 1 for every other model, so this weight is then exactly 1).
		! INVARIANT: the photons a cell removes from the beam are the photons
		! it absorbs. The depths above already carry f (calc_column_dens),
		! so the beam loses f*sigma*n*dr across the cell, and the local
		! ionization, heating and absorbed energy must carry the same f, i.e.
		! the cross section that acts locally is f*sigma as well.  Every
		! integrand below is linear in exactly one cross section, so folding f
		! once into the flux weight they share is identical to multiplying each
		! cross section by it (docs/development_plan_20260905_rev3.md section
		! 10.2 item 7).
		!
		! The attenuation is the MEAN over the cell, exp(-tau_out)
		! (1 - exp(-dtau))/dtau, not the inner-face value: the rate of the
		! cell is the rate averaged over the cell, and only the mean makes
		! the photons the cell absorbs equal the photons the beam loses
		! across it (utils, cell_mean_attenuation).
		int_f = F_XUV*cell_mean_attenuation(tauE_out, dtauE)*opa_pf(j)
		! Secondary-ionization energy partition for this cell. Hydrogen is the
		! only target on this path, so this is the n_HeI -> 0, n_H2 -> 0 limit
		! of photoelectron_energy_partition: the whole ionization energy goes to
		! hydrogen, and fiHI comes back per H I atom [cm^3]. With no molecular
		! hydrogen the H2 terms are identically absent and the heat fraction is
		! Dalgarno's H-He heating efficiency, closed at x = 1.
		if (sec_on) then
			call photoelectron_energy_partition(xion_j, nhi_j, 0.0d0,     &
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
		int_q = int_f*s_hi*nhi_j

		! Value of integrals
		Hea_1 = sum(int_H*de_v)
		PIR_1 = sum(int_1*de_v)
		q_abs = sum(int_q*de_v)

		! Multiply for the dimensional coefficient
		h1_HI_j  = Hea_1*1.0d-18
		heat_j   = h1_HI_j*nhi_j
		P_HI_j   = PIR_1*1.0d-18*erg2eV
		! Add the H I secondary-ionization rate from fast photoelectrons.
		! fiHI is per H I atom, so the integral is already a rate [1/s].
		if (sec_on) then
			acc_secHI = int_f*s_hi*nhi_j/e_v * &
			     merge(fiHI*(e_v-e_th_HI)/e_th_HI, 0.0d0, e_v > e_th_HI + E_sec_ion)
			Psec_HI = sum(acc_secHI*de_v)*1.0d-18*erg2eV
			P_HI_j = P_HI_j + Psec_HI
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
			q_j  = Hea_1*nhi_j/q_abs
		else
			q_j  = 0.0d0
		endif


	end subroutine photoionization_field_at_cell_H

	! ------------------------------------------------------------- !


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

	! SvS85 coupling applied only when enabled AND staged on (see EXHALE_main).
	logical :: sec_on

	! Dummy zero (set below; the number of cells is a runtime value, so these
	! cannot be named constants)
	real*8, dimension(1-Ng:N+Ng) :: nhei, nheii, nheiii, nheiTR
	real*8, dimension(1-Ng:N+Ng) :: N15, N2, NTR

	! The H I column outside the cell being evaluated.
	real*8 :: N1o

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

   ! Evaluate photoionization rates and photoheating rates.  Every column is
	! already known here, so the cells are independent and the loop is the
	! lagged form of the sweep described at photoionization_field_at_cell_H.
	do j = 1-Ng,N+Ng
		N1o = 0.0d0
		if (j .lt. N+Ng) N1o = N1(j+1)
		call photoionization_field_at_cell_H(j, N1o, nhi(j), xion(j),     &
		         sec_on, P_HI(j), h1_HI(j), heat(j), q(j))
	enddo

	if (present(heat_of_one_HI)) heat_of_one_HI = h1_HI

	end subroutine PH_heat_H

	! ------------------------------------------------------------- !
	
	subroutine photoionization_field_at_cell_HHe(j,                       &
	             N1_out,N15_out,N2_out,NTR_out,NH2_out,Nm_out,            &
	             nhi_j,nheiS_j,nheii_j,nheiTR_j,nh2_j,nm_j,xion_j,        &
	             fvq_j,evq_j, has_h2, sec_on, mol_sec,                    &
	             D0_H2_eV, E_ker_H2_dd_eV,                                &
	             P_HI_j,P_HeI_j,P_HeII_j,P_HeITR_j,Pm_j,                  &
	             P_H2_j,P_H2_di_j,P_H2_dd_j,P_H2_nd_j,                    &
	             h1_HI_j,h1_HeI_j,h1_HeII_j,h1_HeTR_j,h1_H2_j,h1m_j,      &
	             heat_j, chan_j, q_j, q_abs_j)
	! The attenuated stellar XUV field of ONE cell, and every rate it drives
	! there: the photoionization rate of each absorber, the photoheating rate
	! of one particle of each of them, the absorbed energy and the heating
	! efficiency.  This is the single definition of that field; PH_heat_HHe
	! is its whole-grid form and the equilibrium sweep calls it cell by cell.
	!
	! WHAT THE FIELD DEPENDS ON, AND WHY THAT MAKES IT SEPARABLE.  The
	! attenuation at this cell is set by the column of each absorber OUTSIDE
	! it, which the caller passes in (N1_out ... Nm_out, the columns
	! evaluated at j+1), and by the cell's OWN column, which it forms here
	! from the densities passed in.  Nothing in it depends on any cell inside
	! j.  A caller sweeping the grid outside in can therefore hand this
	! routine the columns of the composition it has ALREADY solved, and the
	! field each cell is solved at is then self-consistent with that
	! composition in one traversal; a caller that hands it the columns of one
	! whole-grid composition gets the lagged (Jacobi) field instead.  Both
	! have the same fixed point, because at the fixed point the entry and the
	! returned composition are the same state.
	!
	! THE CELL'S OWN DEPTH IS THE SAME STATEMENT ONE CELL DOWN.  dtau is
	! built here from the densities passed in, so a caller that iterates one
	! cell -- solve, re-form the field from what it returned, solve again
	! (ionization_equilibrium, xuv_self_field_passes) -- gets a field whose
	! own attenuation is that of the composition the cell is being solved
	! for.  The cell mean over the two depths is formed ONCE, in one
	! expression, so that iteration re-evaluates that expression and never
	! stacks a second attenuation on the first.
	!
	! Zero columns are what the outermost cell is handed: nothing sits above
	! it, and its star-ward depth is then exactly zero.
	!
	! P_* and h1_* are properties of the FIELD (the density of the absorber
	! that receives them is divided out of h1_*); heat_j and q_j are their
	! contraction with the composition passed in, and describe that
	! composition and no other.

	integer, intent(in) :: j
	! Columns of each absorber outside the cell [cm^-2], opa_pf-weighted.
	real*8, intent(in) :: N1_out,N15_out,N2_out,NTR_out,NH2_out
	real*8, dimension(n_mphot), intent(in) :: Nm_out
	! This cell's own densities.  nheiS_j is the He I GROUND SINGLET, the
	! caller having taken the metastable out of the summed neutral helium
	! (he_ground_singlet_density is the one place that difference is made).
	real*8, intent(in) :: nhi_j,nheiS_j,nheii_j,nheiTR_j,nh2_j,xion_j
	real*8, dimension(n_mion), intent(in) :: nm_j
	! Share of an H2 vibrational excitation collisionally de-excited into
	! heat, and the mean internal energy [eV] one B or C fluorescence leaves.
	real*8, intent(in) :: fvq_j,evq_j
	! Does this run carry molecular hydrogen; is the secondary-ionization
	! coupling on; does it have a molecular target.  Resolved once by the
	! caller, not per cell.
	logical, intent(in) :: has_h2, sec_on, mol_sec
	! D0(H2) and the kinetic energy release of the double-ionization
	! fragments [eV], both constants of the molecule the caller hoisted.
	real*8, intent(in) :: D0_H2_eV, E_ker_H2_dd_eV

	real*8, intent(out) :: P_HI_j,P_HeI_j,P_HeII_j,P_HeITR_j
	real*8, dimension(n_mion), intent(out) :: Pm_j
	real*8, intent(out) :: P_H2_j,P_H2_di_j,P_H2_dd_j,P_H2_nd_j
	real*8, intent(out) :: h1_HI_j,h1_HeI_j,h1_HeII_j,h1_HeTR_j,h1_H2_j
	real*8, dimension(n_mion), intent(out) :: h1m_j
	! Photoheating rate [erg cm^-3 s^-1] and heating efficiency of the
	! composition passed in, and the energy it absorbs [erg cm^-3 s^-1].
	real*8, intent(out) :: heat_j, q_j, q_abs_j
	! The same heating split by absorber: 1 H I, 2 He I (ground singlet),
	! 3 He II, 4 He 2^3S, 5 H2, 6 metals.  The six sum to heat_j exactly,
	! both coming from the one contraction below.
	real*8, dimension(6), intent(out) :: chan_j

	integer :: i,k
	real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR
	real*8 :: PIR_H2, PIR_H2_di, PIR_H2_dd, PIR_H2_nd
	real*8, dimension(Nl) :: int_h2, int_h2_di, int_h2_dd, int_h2_nd
	real*8 :: q_abs
	real*8 :: Pm_loc(n_mion), h1m_loc(n_mion)
	real*8, dimension(Nl) :: tauE_out,dtauE,tau_m
	real*8 :: dr_cell
	real*8, dimension(Nl) :: int_f,int_1,int_15,int_2,int_TR,int_m
	real*8, dimension(Nl) :: int_q,acc_q
	real*8, dimension(Nl) :: acc_secHI,acc_secHeI,acc_secH2,fhv,fh_v
	real*8, dimension(Nl) :: fiHI,fiHeI,fiH2
	real*8, dimension(Nl) :: acc_HI,acc_HeI,acc_HeII,acc_HeTR,acc_H2
	real*8, dimension(Nl) :: acc_mion
	type(photoelectron_partition_t) :: pep
	real*8 :: Psec_HI,Psec_HeI,Psec_H2,Psec_H2_di

	!----------------------------------!

	P_H2_j    = 0.0d0
	P_H2_di_j = 0.0d0
	P_H2_dd_j = 0.0d0
	P_H2_nd_j = 0.0d0
	Pm_j      = 0.0d0

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

		! Optical depth at the cell's STAR-WARD face, from the columns of
		! everything OUTSIDE the cell.  The metal block is accumulated
		! separately in iphot order before the 1e-18 factor.  Nothing sits
		! above the outermost cell, so its caller passes zero columns and the
		! depth is exactly zero there.
		tauE_out = (s_hi*N1_out + s_hei*N15_out                           &
		            + s_heii*N2_out)*1.0d-18
		if (thereis_He .and. thereis_HeITR) tauE_out = tauE_out           &
		                            + s_heiTR*NTR_out*1.0d-18
		if (has_h2) tauE_out = tauE_out                                   &
		                           + s_h2*NH2_out*1.0d-18
		tau_m = 0.0
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			tau_m = tau_m + sigma_tab(:,k)*Nm_out(k)
		enddo
		tauE_out = tauE_out + tau_m*1.0d-18

		! The cell's OWN optical depth, from its own densities and width --
		! the same absorbers, cross sections and opa_pf weight the columns
		! sum.  Formed directly rather than as the difference of the two
		! columns, which is the same number in exact arithmetic and loses
		! its digits to cancellation for a thin cell.
		dr_cell = dr_j(j)*R0*opa_pf(j)
		if (thereis_He) then
			dtauE = (s_hi*nhi_j + s_hei*nheiS_j                         &
			         + s_heii*nheii_j)*dr_cell*1.0d-18
			if (thereis_HeITR) dtauE = dtauE                              &
			                         + s_heiTR*nheiTR_j*dr_cell*1.0d-18
		else
			dtauE = s_hi*nhi_j*dr_cell*1.0d-18
		endif
		if (has_h2) dtauE = dtauE + s_h2*nh2_j*dr_cell*1.0d-18
		tau_m = 0.0
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			tau_m = tau_m + sigma_tab(:,k)*nm_j(i)
		enddo
		dtauE = dtauE + tau_m*dr_cell*1.0d-18

		! Calculate photoionization integrals. int_f is the attenuated
		! stellar flux of this cell weighted by the cell's cross-section
		! factor opa_pf(j) = f(p), the pressure broadening of opacity model
		! 'P' (opacity_pT_factor; f = 1 for every other model, so the weight
		! is then exactly 1).  INVARIANT: the photons a cell removes from the
		! beam are the photons it absorbs.  The depths above already carry f
		! through the columns (calc_column_dens, calc_column_dens_one,
		! calc_column_dens_metals), so the beam loses f*sigma*n*dr across the
		! cell; the local photoionization rates, the photoheating of one
		! particle of each absorber, the secondary-ionization rates and the
		! absorbed energy must carry the same f, i.e. the cross section that
		! acts locally is f*sigma too.  Every integrand below is linear in
		! exactly one cross section, so folding f once into the flux weight
		! they all share is identical to multiplying each cross section by it
		! (docs/development_plan_20260905_rev3.md section 10.2 item 7).
		!
		! The attenuation is the MEAN over the cell, exp(-tau_out)
		! (1 - exp(-dtau))/dtau, and not the inner-face value: the rate of
		! the cell is the rate averaged over the cell, and only the mean
		! makes the photons the cell absorbs equal the photons the beam
		! loses across it, cell by cell and summed over the column (utils,
		! cell_mean_attenuation).
		int_f  = F_XUV*cell_mean_attenuation(tauE_out, dtauE)*opa_pf(j)
		int_1  = int_f*s_hi/e_v
		int_15 = int_f*s_hei/e_v
		int_2  = int_f*s_heii/e_v
		if (thereis_HeITR) int_TR =  int_f*s_heiTR/e_v
		if (has_h2) then
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
			if (has_h2) then
				call photoelectron_energy_partition(xion_j, nhi_j,     &
				                       nheiS_j, nh2_j, fvq_j, evq_j,  &
				                       pep)
			else
				call photoelectron_energy_partition(xion_j, nhi_j,     &
				                       nheiS_j, 0.0d0, 0.0d0, 0.0d0, pep)
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
			acc_secHI  = acc_secHI  + s_hi*nhi_j/e_v *                       &
			     merge(fiHI *(e_v-e_th_HI)/e_th_HI , 0.0d0, e_v > e_th_HI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hi*nhi_j/e_v *                       &
			     merge(fiHeI*(e_v-e_th_HI)/e_th_HeI, 0.0d0, e_v > e_th_HI + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_hi*nhi_j/e_v *                &
			     merge(fiH2 *(e_v-e_th_HI)/e_th_H2 , 0.0d0, e_v > e_th_HI + E_sec_ion)

		fhv = 1.0d0
		if (sec_on) then
			call photoelectron_shares(pep, iabs_HeI, fh_v, fiHI, fiHeI, fiH2)
			fhv = merge(fh_v, 1.0d0, e_v > e_th_HeI + E_sec_ion)
		endif
		acc_HeI = photoelectron_share(e_th_HeI,e_v)*fhv*s_hei
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_hei*nheiS_j/e_v *                    &
			     merge(fiHI *(e_v-e_th_HeI)/e_th_HI , 0.0d0, e_v > e_th_HeI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hei*nheiS_j/e_v *                    &
			     merge(fiHeI*(e_v-e_th_HeI)/e_th_HeI, 0.0d0, e_v > e_th_HeI + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_hei*nheiS_j/e_v *             &
			     merge(fiH2 *(e_v-e_th_HeI)/e_th_H2 , 0.0d0, e_v > e_th_HeI + E_sec_ion)

		fhv = 1.0d0
		if (sec_on) then
			call photoelectron_shares(pep, iabs_HeII, fh_v, fiHI, fiHeI, fiH2)
			fhv = merge(fh_v, 1.0d0, e_v > e_th_HeII + E_sec_ion)
		endif
		acc_HeII = photoelectron_share(e_th_HeII,e_v)*fhv*s_heii
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_heii*nheii_j/e_v *                   &
			     merge(fiHI *(e_v-e_th_HeII)/e_th_HI , 0.0d0, e_v > e_th_HeII + E_sec_ion)
			acc_secHeI = acc_secHeI + s_heii*nheii_j/e_v *                   &
			     merge(fiHeI*(e_v-e_th_HeII)/e_th_HeI, 0.0d0, e_v > e_th_HeII + E_sec_ion)
		endif
		if (mol_sec) acc_secH2 = acc_secH2 + s_heii*nheii_j/e_v *            &
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
				acc_secHI  = acc_secHI  + s_heiTR*nheiTR_j/e_v *             &
				     merge(fiHI *(e_v-e_th_HeTR)/e_th_HI , 0.0d0, e_v > e_th_HeTR + E_sec_ion)
				acc_secHeI = acc_secHeI + s_heiTR*nheiTR_j/e_v *             &
				     merge(fiHeI*(e_v-e_th_HeTR)/e_th_HeI, 0.0d0, e_v > e_th_HeTR + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_heiTR*nheiTR_j/e_v *      &
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
		if (has_h2) then
			fhv = 1.0d0
			if (sec_on) then
				call photoelectron_shares(pep, iabs_H2, fh_v, fiHI, fiHeI, fiH2)
				fhv = merge(fh_v, 1.0d0, e_v > e_th_H2 + E_sec_ion)
			endif
			acc_H2 = photoelectron_share(e_th_H2,e_v)*fhv*(s_h2 - s_h2_di - s_h2_dd - s_h2_nd)
			if (sec_on) then
				acc_secHI  = acc_secHI  + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2_j/e_v *       &
				     merge(fiHI *(e_v-e_th_H2)/e_th_HI , 0.0d0, e_v > e_th_H2 + E_sec_ion)
				acc_secHeI = acc_secHeI + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2_j/e_v *       &
				     merge(fiHeI*(e_v-e_th_H2)/e_th_HeI, 0.0d0, e_v > e_th_H2 + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + (s_h2 - s_h2_di - s_h2_dd - s_h2_nd)*nh2_j/e_v *&
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
				acc_secHI  = acc_secHI  + s_h2_di*nh2_j/e_v *                &
				     merge(fiHI *(e_v-e_th_H2_di)/e_th_HI , 0.0d0, e_v > e_th_H2_di + E_sec_ion)
				acc_secHeI = acc_secHeI + s_h2_di*nh2_j/e_v *                &
				     merge(fiHeI*(e_v-e_th_H2_di)/e_th_HeI, 0.0d0, e_v > e_th_H2_di + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_h2_di*nh2_j/e_v *         &
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
			! The photon buys three things here: the 31.675 eV of chemical
			! and ionization energy the products carry (D0 + 2 I(H)), the
			! kinetic energy of the two electrons, and the kinetic energy
			! release of the two protons.  The 51.4 eV charged above is the
			! VERTICAL threshold, which exceeds the asymptotic product
			! energy by exactly that proton kinetic energy, so the term
			! E_ker_H2_dd_eV/e_v below is the share of the photon the
			! nuclei take and it is thermal from the instant it is
			! released: it carries no fhv, because it never passes through
			! the photoelectron degradation partition.  Without it the
			! channel loses 19.725 eV per event to nowhere
			! (h2_channel_energy_recipients).
			acc_H2 = acc_H2 + (photoelectron_share(e_th_H2_dd,e_v)*fhv    &
			                   + E_ker_H2_dd_eV/e_v)*s_h2_dd
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_h2_dd*nh2_j/e_v *                &
				     merge(fiHI *(e_v-e_th_H2_dd)/e_th_HI , 0.0d0, e_v > e_th_H2_dd + E_sec_ion)
				acc_secHeI = acc_secHeI + s_h2_dd*nh2_j/e_v *                &
				     merge(fiHeI*(e_v-e_th_H2_dd)/e_th_HeI, 0.0d0, e_v > e_th_H2_dd + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + s_h2_dd*nh2_j/e_v *         &
				     merge(fiH2 *(e_v-e_th_H2_dd)/e_th_H2 , 0.0d0, e_v > e_th_H2_dd + E_sec_ion)

			! Neutral dissociation H2 + hv -> H(2p) + H(2s). NO
			! photoelectron: the channel is not an ionization, so no fhv
			! factor and no secondary-ionization terms. The bond energy
			! D0(H2) is spent breaking the molecule, and the fragments of
			! this window are BOTH in n = 2 -- the Q2 1Pi_u(1) state whose
			! fluorescence curve overlaps this cross section -- so a further
			! 2 E(n=1 to 2) = e_rad_H2_neutral leaves the cell as prompt
			! Ly-alpha and two-photon continuum. What is left, the remainder
			! of the photon, is the kinetic energy of the two atoms, i.e.
			! heat (h2_channel_energy_recipients). The channel is nonzero
			! only over 33-41 eV, well above D0 + 2 E(n=1 to 2) = 24.876 eV,
			! so the bracket is positive throughout.
			acc_H2 = acc_H2                                               &
			       + (1.0d0-(D0_H2_eV+e_rad_H2_neutral)/e_v)*s_h2_nd
		endif

		! sigma_tab is the TOTAL photoabsorption of the ion, inner shells
		! included, and every absorption in it is charged here to one
		! ionization of stage i and to a photoelectron of e_v - mion_ethr(i).
		! Above a K or L edge that is the model's approximation of an Auger
		! event, not the event itself; what it under-counts, and why it is
		! what the published photochemistry models do, is at the table in
		! cross_sec.f90.
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
			h1m_loc(i) = sum(int_f*acc_mion*de_v)*1.0d-18
			if (sec_on) then
				acc_secHI  = acc_secHI  + sigma_tab(:,k)*nm_j(i)/e_v *        &
				     merge(fiHI *(e_v-mion_ethr(i))/e_th_HI , 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
				acc_secHeI = acc_secHeI + sigma_tab(:,k)*nm_j(i)/e_v *        &
				     merge(fiHeI*(e_v-mion_ethr(i))/e_th_HeI, 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
			endif
			if (mol_sec) acc_secH2 = acc_secH2 + sigma_tab(:,k)*nm_j(i)/e_v * &
				     merge(fiH2 *(e_v-mion_ethr(i))/e_th_H2 , 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
		enddo
		! Absorbed energy integral (this one is a property of the entry
		! composition: it is what the heating efficiency q is measured
		! against, so it keeps the absorber densities).  The cross sections
		! here are the unbroadened ones; the model 'P' factor f(p) of this
		! cell reaches this integral through int_f in int_q below, the same
		! way it reaches every other absorption integrand, so the energy
		! taken out of the beam and the energy absorbed here are one number.
		acc_q = s_hi *nhi_j
		acc_q = acc_q + s_hei *nheiS_j
		acc_q = acc_q + s_heii*nheii_j
		if (thereis_HeITR) acc_q = acc_q + s_heiTR*nheiTR_j
		if (has_h2) acc_q = acc_q + s_h2*nh2_j
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			acc_q = acc_q + sigma_tab(:,k)*nm_j(i)
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
		if (has_h2) then
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
			Pm_loc(i) = sum(int_m*de_v)*1.0d-18*erg2eV
		enddo
		q_abs   = sum(int_q  *de_v)

		! Photoheating rate of one particle of each absorber [erg s^-1] in
		! this cell's attenuated field (the metal ions were integrated in the
		! loop above). The contraction with the composition follows the loop.
		h1_HI_j   = sum(int_f*acc_HI  *de_v)*1.0d-18
		h1_HeI_j  = sum(int_f*acc_HeI *de_v)*1.0d-18
		h1_HeII_j = sum(int_f*acc_HeII*de_v)*1.0d-18
		h1_HeTR_j = sum(int_f*acc_HeTR*de_v)*1.0d-18
		h1_H2_j   = sum(int_f*acc_H2  *de_v)*1.0d-18
		h1m_j(:)  = h1m_loc(:)
		q_abs_j = q_abs*1.0d-18

		! Save into vectors. No synchronization needed: each thread writes
		! only its own j elements of the shared arrays (the former OMP
		! CRITICAL serialized the loop for no correctness benefit; removed
		! per the 2026-07-02 review's performance note -- values unchanged).
    	P_HI_j    = PIR_1  *1.0d-18*erg2eV
		if (has_h2) P_H2_j = PIR_H2*1.0d-18*erg2eV
		if (has_h2) P_H2_di_j = PIR_H2_di*1.0d-18*erg2eV
		if (has_h2) P_H2_dd_j = PIR_H2_dd*1.0d-18*erg2eV
		if (has_h2) P_H2_nd_j = PIR_H2_nd*1.0d-18*erg2eV
    	P_HeI_j   = PIR_15 *1.0d-18*erg2eV
		! Add the H I / He I / H2 secondary-ionization rates from fast
		! photoelectrons. fiHI, fiHeI and fiH2 are per target particle, so
		! the integrals are already rates [1/s]: no division by a neutral
		! density, and each stays finite as its own target vanishes
		! (electron_energy_degradation.f90).
		if (sec_on) then
			Psec_HI  = sum(int_f*acc_secHI *de_v)*1.0d-18*erg2eV
			Psec_HeI = sum(int_f*acc_secHeI*de_v)*1.0d-18*erg2eV
			P_HI_j  = P_HI_j  + Psec_HI
			P_HeI_j = P_HeI_j + Psec_HeI
		endif
		! The molecular target. P_H2 drives BOTH the H2 destruction row and
		! the H2+ production row of the molecular system, exactly as the
		! primary photoionization rate it is added to does, so the secondary
		! ionizations enter the chemistry through the channel that already
		! exists rather than through a parallel one.
		if (mol_sec) then
			Psec_H2 = sum(int_f*acc_secH2*de_v)*1.0d-18*erg2eV
			! The secondaries dissociatively ionize H2 as well: one proton
			! (and one H atom) per 22 H2+ ions, Dalgarno, Yan & Liu (1999)
			! after their eq. (10). It is an ADDITIONAL yield, not a share of
			! Psec_H2 -- their harmonic-mean statement makes the two H+
			! sources add -- so the H2 destroyed by secondaries is the sum of
			! the two and only the second of them makes protons.
			Psec_H2_di = dissoc_ion_per_H2p*Psec_H2
			P_H2_j = P_H2_j + Psec_H2 + Psec_H2_di
			if (has_h2) P_H2_di_j = P_H2_di_j + Psec_H2_di
		endif
		P_HeII_j  = PIR_2  *1.0d-18*erg2eV
		P_HeITR_j = PIR_TR *1.0d-18*erg2eV
		Pm_j(:)   = Pm_loc(:)

	! The heating of the composition passed in, through the single
	! definition of that contraction, and the efficiency it gives.
	! Absorbed energy is non-negative; see PH_heat_H for why a non-positive
	! q_abs gives a zero heating efficiency.
	call photoheating_of_cell(h1_HI_j,h1_HeI_j,h1_HeII_j,h1_HeTR_j,       &
	                          h1_H2_j,h1m_j,                              &
	                          nhi_j,nheiS_j,nheii_j,nheiTR_j,nh2_j,nm_j,  &
	                          heat_j, chan_j)
	if (q_abs_j .gt. 0.0d0) then
		q_j = heat_j/q_abs_j
	else
		q_j = 0.0d0
	endif

	end subroutine photoionization_field_at_cell_HHe

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
	! (photoheating_of_cell, the single definition of that sum), and
	! so describe THAT composition and no other. A caller that changes the
	! composition afterwards -- an equilibrium sweep, an advection
	! correction -- must re-form the heating from heat_of_one_* and the new
	! densities rather than keep this `heat`. `q` is the heating efficiency
	! of the composition passed in.

	integer :: j

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

   real*8, dimension(1-Ng:N+Ng) :: NH2col  ! H2 column density

	! The columns of everything outside the cell being evaluated, and the
	! cell quantities the field routine returns one at a time.
	real*8 :: N1o,N15o,N2o,NTRo,NH2o
	real*8 :: Nmo(n_mphot)
	real*8 :: nh2_cell
	real*8 :: Pm_row(n_mion), h1m_row(n_mion)
	real*8 :: P_H2_c,P_H2_di_c,P_H2_dd_c,P_H2_nd_c
	real*8 :: chan_c(6)
	! Secondary-ionization coupling applied only when enabled AND staged on
	! (see EXHALE_main). mol_sec additionally requires the run to carry H2,
	! which is what makes the molecular target channel cost nothing, and
	! change nothing, in a run without molecular hydrogen.
	logical :: sec_on, mol_sec
	! Local copy of f_vibq, zero when the caller supplied none.
	real*8, dimension(1-Ng:N+Ng) :: fvq, evq
	! D0(H2) in eV, hoisted out of the cell loop: it is a constant of the
	! molecule, read from the one place this code defines it.
	real*8 :: D0_H2_eV, E_ker_H2_dd_eV

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

	! Optional photoheating breakdown by absorber (cgs erg cm^-3 s^-1).
	! Columns: 1 H I, 2 He I (singlet ground), 3 He II, 4 He 2^3S,
	! 5 H2 (molecular, 0 when absent), 6 metals (sum over photo-ionizable
	! metal ions). `heat` IS the sum of these six columns: both come from the
	! same contraction of heat_of_one_* with the composition, so the split is
	! exact rather than a separately accumulated approximation of it.
	real*8, dimension(1-Ng:N+Ng,6),intent(out),optional :: heat_chan

	!----------------------------------!

	! Ground singlet n(1^1S): the summed neutral helium when the metastable
	! is not tracked, and the clamped difference n(He I) - n(2^3S) when it
	! is. composition.f90 is the only place that difference is taken.
	nheiS = nhei
	if (thereis_HeITR) nheiS = he_ground_singlet_density(nhei, nheiTR)

	! Evaluate the column density
	call calc_column_dens(nhi,nheiS,nheii,nheiTR,N1,N15,N2,NTR)
	call calc_column_dens_metals(nm, Nm_col)

	! H2 column (molecular).  Zero, and read as zero, in a run without it.
	NH2col = 0.0d0
	if (present(nh2)) call calc_column_dens_one(nh2, NH2col)

	! Metal-free P_m entries (top-stage ions) stay zero
	P_m = 0.0

	D0_H2_eV = h2_dissociation_energy_eV()
	E_ker_H2_dd_eV = h2_double_fragment_kinetic_energy()
	sec_on  = use_sec_ion .and. sec_ion_active
	mol_sec = sec_on .and. present(nh2) .and. present(P_H2)
	fvq     = 0.0d0
	evq     = 0.0d0
	if (present(f_vibq)) fvq = f_vibq
	if (present(e_vibq)) evq = e_vibq

    !----------------------------------!

	! The field of every cell, from the columns of the ONE composition this
	! routine was given: the lagged form of the sweep described at
	! photoionization_field_at_cell_HHe.  Every column is already known here,
	! so the cells are independent and the loop is parallel.  No
	! synchronization is needed: each thread writes only its own j elements
	! of the shared arrays.
	!$OMP PARALLEL DO &
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,P_HeITR,P_m, sec_on, mol_sec, fvq, evq, &
	!$OMP          D0_H2_eV, E_ker_H2_dd_eV ) &
	!$OMP PRIVATE ( j, N1o,N15o,N2o,NTRo,NH2o,Nmo, nh2_cell,                 &
	!$OMP           Pm_row,h1m_row, P_H2_c,P_H2_di_c,P_H2_dd_c,P_H2_nd_c,    &
	!$OMP           chan_c )
	do j = 1-Ng,N+Ng

		! The columns of everything outside the cell.  Nothing sits above the
		! outermost cell, so it is handed zeros.
		if (j .lt. N+Ng) then
			N1o  = N1(j+1)
			N15o = N15(j+1)
			N2o  = N2(j+1)
			NTRo = NTR(j+1)
			NH2o = NH2col(j+1)
			Nmo  = Nm_col(j+1,:)
		else
			N1o  = 0.0d0
			N15o = 0.0d0
			N2o  = 0.0d0
			NTRo = 0.0d0
			NH2o = 0.0d0
			Nmo  = 0.0d0
		endif

		nh2_cell = 0.0d0
		if (present(nh2)) nh2_cell = nh2(j)

		call photoionization_field_at_cell_HHe(j,                          &
		         N1o,N15o,N2o,NTRo,NH2o,Nmo,                               &
		         nhi(j),nheiS(j),nheii(j),nheiTR(j),nh2_cell,nm(j,:),      &
		         xion(j), fvq(j),evq(j), present(nh2), sec_on, mol_sec,    &
		         D0_H2_eV, E_ker_H2_dd_eV,                                 &
		         P_HI(j),P_HeI(j),P_HeII(j),P_HeITR(j),Pm_row,             &
		         P_H2_c,P_H2_di_c,P_H2_dd_c,P_H2_nd_c,                     &
		         h1_HI(j),h1_HeI(j),h1_HeII(j),h1_HeTR(j),h1_H2(j),h1m_row,&
		         heat(j), chan_c, q(j), q_abs_cell(j))

		P_m(j,:)  = Pm_row
		h1_m(j,:) = h1m_row
		if (present(P_H2))    P_H2(j)    = P_H2_c
		if (present(P_H2_di)) P_H2_di(j) = P_H2_di_c
		if (present(P_H2_dd)) P_H2_dd(j) = P_H2_dd_c
		if (present(P_H2_nd)) P_H2_nd(j) = P_H2_nd_c
		if (present(heat_chan)) heat_chan(j,:) = chan_c

	enddo
	!$OMP END PARALLEL DO

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

	real*8, dimension(1-Ng:N+Ng) :: nheiS
	real*8  :: chan_j(6)
	integer :: j

	nheiS = nhei
	if (thereis_HeITR) nheiS = he_ground_singlet_density(nhei, nheiTR)

	do j = 1-Ng,N+Ng
		call photoheating_of_cell(h1_HI(j),h1_HeI(j),h1_HeII(j),          &
		         h1_HeTR(j),h1_H2(j),h1_m(j,:),                           &
		         nhi(j),nheiS(j),nheii(j),nheiTR(j),nh2(j),nm(j,:),       &
		         heat(j), chan_j)
		if (present(heat_chan)) heat_chan(j,:) = chan_j
	enddo

	end subroutine photoheating_of_composition

	!----------------------------------!

	pure subroutine photoheating_of_cell(h1_HI,h1_HeI,h1_HeII,h1_HeTR,    &
	                                     h1_H2, h1_m,                     &
	                                     nhi,nheiS,nheii,nheiTR,nh2,nm,   &
	                                     heat, heat_chan)
	! Photoheating rate of one cell's composition [erg cm^-3 s^-1]: the
	! heating rate of ONE particle of each absorber (a property of the
	! radiation field alone) times the density of that absorber.
	!
	! THIS IS THE ONE PLACE THE CONTRACTION IS WRITTEN.  The whole-grid form
	! above and the field routine that reports the efficiency of the state it
	! was handed both come here, so `heat` and its six-column split are one
	! object and the split is exact rather than a separate approximation.
	!
	! He I absorbs in its GROUND SINGLET: the density asked for here is
	! n(1^1S), the caller having taken the metastable out of the summed
	! neutral helium (he_ground_singlet_density is the one place that
	! difference is made).  The metastable has its own column.

	real*8, intent(in) :: h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2
	real*8, dimension(n_mion), intent(in) :: h1_m
	real*8, intent(in) :: nhi,nheiS,nheii,nheiTR,nh2
	real*8, dimension(n_mion), intent(in) :: nm
	real*8, intent(out) :: heat
	! Split by absorber: 1 H I, 2 He I (ground singlet), 3 He II,
	! 4 He 2^3S, 5 H2, 6 metals (summed over the ions).
	real*8, dimension(6), intent(out) :: heat_chan

	real*8  :: h_HI,h_HeI,h_HeII,h_HeTR,h_H2,h_mtl
	integer :: i

	h_HI   = h1_HI  *nhi
	h_HeI  = h1_HeI *nheiS
	h_HeII = h1_HeII*nheii
	h_HeTR = h1_HeTR*nheiTR
	h_H2   = h1_H2  *nh2
	h_mtl  = 0.0d0
	do i = 1,n_mion
		h_mtl = h_mtl + h1_m(i)*nm(i)
	enddo

	heat = h_HI + h_HeI + h_HeII + h_HeTR + h_H2 + h_mtl

	heat_chan(1) = h_HI
	heat_chan(2) = h_HeI
	heat_chan(3) = h_HeII
	heat_chan(4) = h_HeTR
	heat_chan(5) = h_H2
	heat_chan(6) = h_mtl

	end subroutine photoheating_of_cell

	!----------------------------------!

	subroutine heating_of_composition(T_K,                                &
	         nhi,nhii,nhei,nheii,nheiii,nheiTR, nm, nmol, nox, ne, n_tot, &
	         h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2,h1_m,                     &
	         A31,q31a,q31b,Q31, k_lw, p_lw, k_co, j_h2o, j_oh,            &
	         with_molecules, with_oxygen, heat, heat_chan)
	! Volumetric heating rate of a composition [erg cm^-3 s^-1], channel by
	! channel, and its total.
	!
	! THIS IS THE ONE PLACE THE HEATING IS ASSEMBLED. Every consumer calls
	! it: the ionization sweep for the heat it returns, the breakdown dump
	! for what it writes, the advection-corrected post-process for the
	! heating its energy solve balances. Each channel is formed exactly once
	! and the total is the running sum of the channels, so a deposit added
	! to the energy equation appears in the breakdown by construction.
	!
	! WHAT IS A RATE AND WHAT IS A DENSITY. The arguments split in two. The
	! photoheating of ONE particle of each absorber (h1_*), the He(2^3S)
	! coefficients (A31, q31a, q31b, Q31), the Lyman-Werner dissociation
	! rate and its single-pump branching (k_lw, p_lw) and the FUV band
	! photolysis rates (j_h2o, j_oh) are properties of the radiation field
	! and of T. The densities are whichever composition the caller wants the
	! heating of. Each channel is a contraction of the two, so the caller
	! that has just solved a sweep passes its post-sweep densities against
	! the rates the sweep was solved with, and the lag of one sweep is
	! documented there and not repeated here.
	!
	! The H(n=2) channels are the exception and are read from the module
	! state: Hpe_arr and Hdx_arr are already-contracted volumetric rates
	! filled by excited_H_update from the previous converged outer pass, so
	! there is nothing for a caller to contract.
	!
	! WHAT THE COMPOSITION CARRIES. with_molecules and with_oxygen say
	! whether the densities passed include the molecular carriers (H2, H2+,
	! H3+, HeH+) and the oxygen carriers (OH, H2O, CO). The deposits of a
	! network whose carriers are not represented are left at zero rather
	! than evaluated against zero densities, because two of them do not
	! vanish there: the three-body H2 formation of the collisional network
	! runs on n(H I) squared, and the oxygen network's O2 channel on the
	! free atomic oxygen. The ionization sweep passes the configuration
	! switches; the advection-corrected post-process reconstructs an
	! H/He + metals composition only (see its header) and passes .false.,
	! the same approximation its densities, its n_e and its n_tot carry.

	real*8, dimension(1-Ng:N+Ng),intent(in) :: T_K
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhii,nhei,nheii,       &
	                                            nheiii,nheiTR, ne, n_tot
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm
	! Molecular carriers (H2, H2+, H3+, HeH+) and oxygen carriers (OH, H2O,
	! CO) [cm^-3]; zero for a run without them.
	real*8, dimension(1-Ng:N+Ng,4),intent(in) :: nmol
	real*8, dimension(1-Ng:N+Ng,3),intent(in) :: nox
	! Photoheating of one particle of each absorber [erg s^-1]
	real*8, dimension(1-Ng:N+Ng),intent(in) :: h1_HI,h1_HeI,h1_HeII,      &
	                                            h1_HeTR,h1_H2
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: h1_m
	! He(2^3S) radiative and collisional coefficients (HeITR_coeffs): A31 is
	! the 2^3S -> 1^1S decay rate, q31a/q31b its collisional de-excitation
	! and Q31 the TOTAL He(2^3S) + H ionization rate coefficient.
	real*8, intent(in) :: A31
	real*8, dimension(1-Ng:N+Ng),intent(in) :: q31a,q31b,Q31
	! Lyman-Werner dissociation rate [s^-1] and dissociations per pump
	real*8, dimension(1-Ng:N+Ng),intent(in) :: k_lw, p_lw
	! CO photodissociation rate [s^-1], the cell mean fuv_lw_photon_field
	! returns; zero for a composition with no CO or no band flux.
	real*8, dimension(1-Ng:N+Ng),intent(in) :: k_co
	! Band-resolved H2O and OH photodissociation rates [s^-1]
	real*8, dimension(1-Ng:N+Ng,n_fuv_band),intent(in) :: j_h2o, j_oh
	logical, intent(in) :: with_molecules, with_oxygen

	real*8, dimension(1-Ng:N+Ng),intent(out) :: heat
	! Deposits in the order of heat_channel_name. Not optional: the total is
	! their sum, so a caller that wants the total has the columns too.
	real*8, dimension(1-Ng:N+Ng,n_heat_channel),intent(out) :: heat_chan

	! Scratch of the He recombination coupling: only its heating is kept
	! here, the rate corrections belong to the sweep that solved with them.
	! Ground-singlet neutral helium of this composition, formed here for
	! the one chemical-heat channel that has a neutral-helium reactant.
	real*8, dimension(1-Ng:N+Ng) :: nheiS_chem
	real*8, dimension(1-Ng:N+Ng) :: rcheiiB_hrc,dP_HI_hrc,dP_H2_hrc
	real*8, dimension(1-Ng:N+Ng,n_mion) :: dP_m_hrc
	! Production of O(1D) [cm^-3 s^-1] through the H2O + hv -> H2 + O(1D)
	! branch, the flux its local steady state carries into the O6 sink.
	real*8, dimension(1-Ng:N+Ng) :: flux_o1d
	integer :: ib

	heat_chan = 0.0d0

	! Photoheating: the heating rate of one particle of each absorber times
	! the density of that absorber, split by absorber into the first six
	! channels.
	call photoheating_of_composition(h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2,  &
	                                 h1_m, nhi,nhei,nheii,nheiTR,         &
	                                 nmol(:,1), nm, heat,                 &
	                                 heat_chan = heat_chan(:,1:6))

	! H(n=2) Balmer heating: photoelectric heating of the excited level plus
	! its Lyman-alpha de-excitation. Unlike every other term here it is an
	! already-contracted volumetric rate, filled by excited_H_update from the
	! previous converged outer pass, so it carries that documented lag and is
	! added as it stands. Zero unless use_excited_H.
	if (use_excited_H) then
		heat_chan(:,7) = Hpe_arr
		heat_chan(:,8) = Hdx_arr
		heat = heat + (heat_chan(:,7) + heat_chan(:,8))
	endif

	! He recombination radiation absorbed by H I, H2 and the metal ions. The
	! RATE corrections of this coupling belong to the state that entered the
	! sweep and are not rebuilt here; its photoelectron heating carries the
	! densities of the absorbers, so it is evaluated at the composition this
	! routine was given and only the heating is kept.
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nmol(:,1), nhei, nheii, nheiTR,     &
		                     ne, nm, A31, q31a, q31b,                      &
		                     rcheiiB_hrc, dP_HI_hrc, dP_H2_hrc,            &
		                     dP_m_hrc, heat_chan(:,9))
		heat = heat + heat_chan(:,9)
	endif

	! Penning ionization heating: He(2^3S)+H0 -> He(1^1S)+H+ + e- releases the
	! electron kinetic energy e_th_HeI - e_th_HeTR - e_th_HI (= 6.2 eV) into the
	! gas. Q31 is the total ionization rate, so only its Penning branch
	! (f_penning_HeI23S) carries this exothermicity; the associative branch
	! ends in HeH+ and has a different one.
	if (thereis_HeITR) then
		heat_chan(:,10) = f_penning_HeI23S*nheiTR*nhi*Q31                  &
		                  *(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV
		heat = heat + heat_chan(:,10)
	endif

	! ASSOCIATIVE branch of the same collision, He(2^3S)+H0 -> HeH+ + e-,
	! the remaining (1 - f_penning_HeI23S) of Q31 (D0 C25, b1 T1.9 item 3).
	! Its heat is a difference of the ONE species formation-energy table and
	! is not written down here; what the two branches below choose is which
	! reaction the configuration actually runs.
	!
	!   MOLECULAR run: HeH+ is a species of the network (System_HeH_mol row
	!   7), so the event that is deposited here is exactly
	!   He(2^3S) + H -> HeH+ + e, 8.07 eV.  What becomes of that HeH+ is the
	!   business of R16, R18 and R19, whose heats molecular_chemical_heating
	!   deposits, so nothing is counted twice.
	!
	!   ATOMIC run: this system carries no HeH+, and its own closure
	!   (ion_residual_core.f90, the Penning comment) is that in an atomic gas
	!   the HeH+ dissociatively recombines straight back to He + H, faster
	!   than any competing reaction, so the branch is a metastable sink and
	!   not a lasting ion source.  The reaction the code solves is therefore
	!   He(2^3S) + H -> He + H, whose heat is eps(He 2^3S) alone, 19.82 eV.
	!   Depositing only the 8.07 eV there would strand the 11.75 eV that the
	!   recombination step returns and that no other term of an atomic run
	!   carries.
	!
	! The branch below is thereis_mol and not with_molecules: what it asks
	! is whether the SYSTEM carries HeH+ as a species, which is a property
	! of the network the run solves and not of whether the composition
	! handed to this routine lists the molecular carriers.
	if (thereis_HeITR) then
		if (thereis_mol) then
			heat_chan(:,11) = (1.0d0 - f_penning_HeI23S)*nheiTR*nhi*Q31    &
			                  *oxygen_reaction_energy_eV(ir_assoc_HeTR)    &
			                  /erg2eV
		else
			heat_chan(:,11) = (1.0d0 - f_penning_HeI23S)*nheiTR*nhi*Q31    &
			                  *species_formation_energy(isp_HeTR)/erg2eV
		endif
		heat = heat + heat_chan(:,11)
	endif

	! Molecular Penning ionization heating: He(2^3S)+H2 -> He(1^1S)+H2+ + e-
	! releases the electron kinetic energy (e_th_HeI - e_th_HeTR) - e_th_H2
	! (= 24.5874 - 4.7678 - 15.4259 = 4.394 eV) into the gas. nmol(:,1) is the
	! neutral-H2 number density. ioniz_HeI23S_H2 (Cool_coeff.f90) is the
	! total, scaled here to the Penning branch. Zero unless a molecular run
	! also tracks the triplet.
	if (with_molecules .and. thereis_HeITR) then
		heat_chan(:,12) = f_penning_HeI23S*nheiTR*nmol(:,1)                &
		                  *ioniz_HeI23S_H2(T_K)                            &
		                  *((e_th_HeI - e_th_HeTR) - e_th_H2)/erg2eV
		heat = heat + heat_chan(:,12)
	endif

	! Lyman-Werner photodissociation heating: H2 + hv -> H + H leaves the
	! fragment pair with about 0.4 eV of kinetic energy (Black & Dalgarno
	! 1977, ApJS 34, 405, p. 418). The 4.48 eV bond energy is paid by the
	! absorbed photon, not by the gas, so it is NOT a thermal sink of this
	! channel.
	if (with_molecules .and. F_LW_star .gt. 0.0d0) then
		heat_chan(:,13) = k_lw*nmol(:,1)*e_lw_fragment_erg
		heat = heat + heat_chan(:,13)
	endif

	! Lyman-Werner FLUORESCENCE heating, the other and much larger half of
	! the same absorption. k_lw is a DISSOCIATION rate, and every pump
	! that does not dissociate fluoresces back into a bound, vibrationally
	! excited level of the ground state carrying the mean landing-level energy
	! h2_energy_per_bound_fluorescence_erg(T) -- 2.06 eV at 700 K to 2.13 eV
	! at 3200 K, computed from the Abgrall, Roueff & Drira (2000) transition
	! probabilities rather than adopted, and replacing the 2.0 eV Burton,
	! Hollenbach & Tielens (1990) Appendix A adopt without a derivation
	! (docs/p39_lw_cross_section_sources.md sec. 5.2). The count of such
	! decays per dissociation is (1 - p)/p. At the density of a molecular
	! base that energy is collisionally de-excited and becomes heat; at low
	! density it is radiated away in the infrared quadrupole lines, and
	! h2_vibrational_heat_fraction is the ratio between the two.
	!
	! p IS THE SINGLE-PUMP BRANCHING, not the effective one. Once the layer
	! is thick a fluorescent photon can be re-absorbed, but that pair of
	! events -- one molecule down, another up -- deposits nothing, so the
	! trapping cancels out of the count and what survives is the branching a
	! lone pump would have. It is not the constant 0.135 either: shielding
	! removes the strongest pumping lines first, and the lines that survive
	! to depth have a different branching, measured as p_lw
	! (docs/p39_lw_cross_section_sources.md).
	! THE GROUND-SINGLET NEUTRAL HELIUM, formed once for every molecular
	! channel of this assembly that needs it.  It is the third collider of
	! the H2 vibrational cascade (Jozwiak et al. 2024) and the neutral
	! reactant of H2+ + He -> HeH+ + H; a helium ION quenches an H2
	! vibration through a different interaction and carries no coefficient
	! here, and the metastable's own channels are the He(2^3S) ones.
	nheiS_chem = nhei
	if (thereis_HeITR) nheiS_chem = he_ground_singlet_density(nhei, nheiTR)

	if (with_molecules .and. F_LW_star .gt. 0.0d0) then
		heat_chan(:,14) = k_lw*(1.0d0 - p_lw)                              &
		                  /max(p_lw, 1.0d-30)                              &
		                  *nmol(:,1)                                       &
		                  *h2_energy_per_bound_fluorescence_erg(T_K)       &
		                  *h2_vibrational_heat_fraction(T_K, nhi,          &
		                                       nmol(:,1), nheiS_chem)
		heat = heat + heat_chan(:,14)
	endif

	! Chemical heat of the COLLISIONAL reactions of the H2/He network. The
	! photon-driven reactions, the radiative recombinations and the
	! collisional ionizations are excluded there, so nothing above is counted
	! twice; see molecular_reaction_heat.f90 for the exclusion list and for
	! why an atomic gas has no such term. Default off.
	if (with_molecules .and. mol_reaction_heat) then
		! nheiS_chem, the ground-singlet neutral helium, is formed above:
		! it is the collision partner of H2+ + He -> HeH+ + H and the
		! monatomic third body of the R12/R15 pair.
		call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS_chem, &
		                                nmol, ne, n_tot, heat_chan(:,15))
		heat = heat + heat_chan(:,15)
	endif

	if (with_oxygen) then
		! FUV photolysis heating: each H2O or OH dissociation leaves the
		! fragments with the excess of the absorbed photon over the bond
		! energy, the same ledger the Lyman-Werner channel above uses (the
		! bond energy is paid by the photon, not by the gas, so it is not a
		! thermal sink). The H2 + O(1D) branch keeps its 1.96 eV of
		! electronic excitation out of this sum: that energy leaves as O(1D)
		! and is released later, in the O(1D) + H2 -> OH + H reaction. THAT
		! exothermicity is NOT deposited by this network -- an omission of
		! the same kind as the missing thermal dissociation sink of R12/R14,
		! and of the same size (a few percent of the photolysis heat),
		! recorded here rather than hidden.
		do ib = 1,n_fuv_band
			heat_chan(:,16) = heat_chan(:,16)                              &
			  + j_h2o(:,ib)*nox(:,2)                                       &
			    *heat_per_water_dissociation(ib)                           &
			  + j_oh(:,ib)*nox(:,1)                                        &
			    *heat_per_hydroxyl_dissociation(ib)
		enddo
		heat = heat + heat_chan(:,16)

		! COLLISIONAL oxygen channels O1, O1r, O2, O2r and the O(1D) sink
		! O6 (D0 C25, b1 T1.9).  Each is a difference of the one species
		! formation-energy table, so nothing is transcribed and the forward
		! and reverse channels are exact negatives.  The photolysis channels
		! are NOT in this sum: their enthalpy was paid by the absorbed photon
		! and channel 16 above deposits the excess.
		!
		! O(1D) is eliminated by its local steady state (b1 T1.7), and an
		! elimination transfers the reservoir with the nuclei (T1.9 item 1):
		! the flux through its single sink O6 equals its production
		! oj4*n(H2O), and it carries eps(O(1D)) to the products.
		! System_HeH_mol builds the same oj4 from the same quantum yields for
		! its rows.
		!
		! nm(:,im_OI) is the FREE atomic oxygen, the O I column with the
		! oxygen bound in OH, H2O and CO already removed, which is the
		! density the O2 channel runs on in oxygen_carrier_rows.
		flux_o1d = 0.0d0
		do ib = 1,n_fuv_band
			flux_o1d = flux_o1d + qy_H2O_H2_O1D(ib)*j_h2o(:,ib)
		enddo
		flux_o1d = flux_o1d*nox(:,2)
		call oxygen_chemical_heating(T_K, nhi, nmol(:,1),                  &
		                             nm(:,im_OI), nox, flux_o1d,           &
		                             heat_chan(:,17))
		heat = heat + heat_chan(:,17)

		! The two CO destruction channels of the one-sided CO model
		! (docs/b3b_co_destruction_design_20260906.md).  Both are formed
		! HERE and nowhere else, and both take their energy from the one
		! formation-energy table, so they cannot state two different C=O
		! bond energies.
		!
		! 18: He+ + CO -> C+ + O + He, UMIST RATE22 entry 4068.  The heat is
		! eps(He+) + eps(CO) - eps(C+) - eps(O) - eps(He) = +2.2117 eV, the
		! difference of the helium and carbon ionization potentials less the
		! C=O bond energy; the products carry it as kinetic energy.  nox(:,3)
		! is CO and nheii the He+ of the state.
		!
		! THE EVENT RATE IS THE ONE THE BALANCE ROWS USE.  k_D1 n(He+) n(CO)
		! with the same rk_D1_Hep_CO is the He+ sink of row (2) of
		! mol_heh_rows and the CO sink of the carrier row
		! (diffusive_photochemistry::carrier_source), so the number of
		! events this energy is deposited for is the number of events the
		! species ledger performs.
		heat_chan(:,18) = rk_D1_Hep_CO()*nheii*nox(:,3)                    &
		                  *oxygen_reaction_energy_eV(ir_D1)/erg2eV
		heat = heat + heat_chan(:,18)

		! 19: CO + hv -> C + O.  The photon pays the 11.1157 eV bond and the
		! fragments keep the rest, which is the same ledger the
		! Lyman-Werner and the H2O/OH photolysis channels above use.  The
		! bond energy is NOT a thermal sink of this channel.
		heat_chan(:,19) = k_co*nox(:,3)*heat_per_co_dissociation()
		heat = heat + heat_chan(:,19)
	endif

	end subroutine heating_of_composition

	!----------------------------------!

	subroutine fuv_band_absorption_ledger(T_K, nhi, nh2, nh2o, nheiS,    &
	         nH_nuc,                                                     &
	         NH2col, NH2Ocol, NOHcol, NCOcol, tau_b, j_h2o, j_oh,        &
	         k_lw, p_lw_single, k_co,                                    &
	         absph, absen, heat_col, bond_col, cont_ph, cont_beam,       &
	         co_ph, co_en, co_heat, e_lw_abs, drift_worst, drift_r)
	! The COLUMN-INTEGRATED photon and energy ledger of the four FUV bands,
	! per unit area of the star-ward column: what each band's absorbers take
	! out of the beam, what that carries, and how much of it reaches the gas.
	! output/FUV_bands.txt writes what this returns and forms none of it.
	!
	! WHY IT LIVES BESIDE THE HEATING ASSEMBLY. The energy of one absorption
	! event -- the excess of a photolysis photon over the bond, the kinetic
	! energy of a Lyman-Werner fragment pair, the fluorescence that follows
	! the pumps that do not dissociate, the CO fragment energy -- is the same
	! energy the energy equation is charged. Written a second time in the
	! output module it can drift from the equation exactly as the three
	! copies of the heating sum did, and that is what happened: until this
	! routine existed the scanner
	! src/tests/physics_probe/heating_sum_uniqueness.py had to record
	! write_output.f90 as an exception. It is one file, one set of imports
	! and one place those energies are read.
	!
	! WHAT IT IS NOT. It is not a second heating assembly. This is a sum over
	! the GRID at fixed band; heating_of_composition is a sum over the BANDS
	! at fixed cell. The two answer different questions and neither can be
	! got from the other: the ledger needs each band's own absorbed photons
	! to compare with that band's beam loss, and the energy equation needs
	! each cell's total. What they share, and what this move makes single, is
	! the energy of one event, which each of them multiplies a rate by.
	!
	! THE STATE IS THE CALLER'S. Every density and every column is passed in,
	! because the ledger describes the state the output file writes and that
	! is the caller's to choose. drift_worst reports how far the H2O density
	! written has moved from the one the photon field was built on.
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_K, nhi, nh2, nh2o
	! Ground-singlet neutral helium, the third collider of the H2
	! vibrational cascade in the fluorescence heat of the Lyman-Werner band
	! (Jozwiak et al. 2024).  The same collider the heating assembly and the
	! chemistry pass, so one collider sum decides the thermalized share.
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiS
	! Total hydrogen NUCLEUS density, the density axis of the H2
	! self-shielding table the pump cross section is read from.
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nH_nuc
	real*8, dimension(1-Ng:N+Ng), intent(in) :: NH2col, NH2Ocol, NOHcol
	real*8, dimension(1-Ng:N+Ng), intent(in) :: NCOcol
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(in) :: tau_b
	real*8, dimension(1-Ng:N+Ng,n_fuv_band), intent(in) :: j_h2o, j_oh
	real*8, dimension(1-Ng:N+Ng), intent(in) :: k_lw, p_lw_single, k_co
	real*8, dimension(n_fuv_band), intent(out) :: absph, absen
	real*8, dimension(n_fuv_band), intent(out) :: heat_col, bond_col
	real*8, dimension(n_fuv_band), intent(out) :: cont_ph, cont_beam
	real*8, intent(out) :: co_ph, co_en, co_heat, e_lw_abs
	real*8, intent(out) :: drift_worst, drift_r
	real*8  :: dr_cm, a_cell, dtau_col, dn_h2o, dn_oh, dn_h2, dn_co
	real*8  :: dstate, ph_lw, NH2_out_lw, tau_out_lw
	integer :: j, ib

      absph    = 0.0d0
      absen    = 0.0d0
      heat_col = 0.0d0
      bond_col = 0.0d0
      cont_ph   = 0.0d0
      cont_beam = 0.0d0
      e_lw_abs = 0.0d0
      co_ph    = 0.0d0
      co_en    = 0.0d0
      co_heat  = 0.0d0
      drift_worst = 0.0d0
      drift_r     = 0.0d0
      do j = 1-Ng,N+Ng
         dr_cm = dr_j(j)*R0*opa_pf(j)
         ! THE ABSORBERS OF THE BEAM ARE THE ONES THE COLUMNS RECORD, so a
         ! cell's absorber count here is its own column increment and not
         ! its written density times the cell width. The column at cell j
         ! already contains cell j, so that increment is N(j) - N(j+1), and
         ! N(N+Ng) is the outermost cell's own content because nothing sits
         ! above it. The two forms agree wherever the state written is the
         ! state the photon field was built on; where they do not, the beam
         ! and its own record are the pair that has to balance, and the
         ! drift between the two states is reported on its own line below.
         if (j .lt. N+Ng) then
            dn_h2o = NH2Ocol(j)   - NH2Ocol(j+1)
            dn_oh  = NOHcol(j)    - NOHcol(j+1)
            dn_h2  = NH2col(j) - NH2col(j+1)
            dn_co  = NCOcol(j)    - NCOcol(j+1)
         else
            dn_h2o = NH2Ocol(j)
            dn_oh  = NOHcol(j)
            dn_h2  = NH2col(j)
            dn_co  = NCOcol(j)
         endif
         dstate = abs(nh2o(j)*dr_cm - dn_h2o)/max(abs(dn_h2o), 1.0d-99)
         if (dn_h2o .gt. 0.0d0 .and. dstate .gt. drift_worst) then
            drift_worst = dstate
            drift_r     = r(j)
         endif
         do ib = 1,n_fuv_band
            a_cell = j_h2o(j,ib)*dn_h2o + j_oh(j,ib)*dn_oh
            absph(ib)   = absph(ib)   + a_cell
            cont_ph(ib) = cont_ph(ib) + a_cell
            ! The same cell written as the beam's loss between its two
            ! faces (paragraph (a) above). j_H2O/sigma_H2O is
            ! N_b tr exp(-tau_out) (1 - exp(-dtau))/dtau with the band's own
            ! cross section divided out, so multiplying it by the cell's
            ! continuum depth gives that loss with no reference to a density.
            if (j .lt. N+Ng) then
               dtau_col = max(tau_b(j,ib) - tau_b(j+1,ib), 0.0d0)
            else
               dtau_col = max(tau_b(j,ib), 0.0d0)
            endif
            cont_beam(ib) = cont_beam(ib)                                  &
                          + j_h2o(j,ib)/sigma_H2O_band(ib)*dtau_col
            ! Each absorbed photon carries the band's own mean energy
            ! <hv>_b, so the total can never exceed the incident flux; see
            ! the conservation argument in water_photolysis.f90 sec. 2.
            absen(ib) = absen(ib) + a_cell*e_photon_flat_band(ib)
            heat_col(ib) = heat_col(ib)                                    &
                      + j_h2o(j,ib)*dn_h2o                             &
                        *heat_per_water_dissociation(ib)                   &
                      + j_oh(j,ib)*dn_oh                               &
                        *heat_per_hydroxyl_dissociation(ib)
         enddo
         if (F_LW_star .gt. 0.0d0) then
            ! The H2 share of the shared LW beam, on the same ledger as the
            ! continuum absorbers. Two DIFFERENT branchings are needed and
            ! they are not the same number (lyman_werner.f90, and
            ! docs/p39_lw_cross_section_sources.md):
            !   p_lw_single    how many fluorescent decays accompany each
            !                  dissociation, (1 - p)/p, which is the term the
            !                  energy equation carries.  Both terms are in
            !                  the energy equation, so both belong here;
            !                  before the second was carried, bond_col
            !                  counted it as energy that never comes back.
            !   sigma_pump     sigma_diss/p_eff, the photons the lines take
            !                  out of the beam.  The dissociation rate is NOT
            !                  divided by p_eff here: a cell of a molecular
            !                  base spans decades of column over which p_eff
            !                  varies, so the cell mean of sigma_diss divided
            !                  by p_eff at one column is not the cell mean of
            !                  the ratio.
            !                  The pump rate is taken as its own cell mean,
            !                  on the same quadrature and the same faces as
            !                  the dissociation rate
            !                  (lyman_werner.f90), which is what makes the
            !                  rated photons telescope to the beam's loss.
            !                  The column and the depth are re-read from the
            !                  state this file writes, exactly as the
            !                  continuum absorbers above are.
            if (j .lt. N+Ng) then
               NH2_out_lw = NH2col(j+1)
               tau_out_lw = tau_b(j+1,ib_LW)
            else
               NH2_out_lw = 0.0d0
               tau_out_lw = 0.0d0
            endif
            ph_lw    = lyman_werner_band_absorption_rate_cell_mean(        &
                          fuv_band_flux(ib_LW), NH2_out_lw, NH2col(j), &
                          T_K(j), nH_nuc(j),             &
                          tau_out_lw, tau_b(j,ib_LW))                    &
                       *dn_h2
            e_lw_abs = e_lw_abs + ph_lw*e_lw_photon_erg
            absph(ib_LW)    = absph(ib_LW)    + ph_lw
            absen(ib_LW)    = absen(ib_LW)    + ph_lw*e_lw_photon_erg
            heat_col(ib_LW) = heat_col(ib_LW)                              &
                            + k_lw(j)*dn_h2                           &
                              *e_lw_fragment_erg                           &
                            + k_lw(j)*dn_h2                           &
                              *(1.0d0 - p_lw_single(j))                    &
                              /max(p_lw_single(j), 1.0d-30)                &
                              *h2_energy_per_bound_fluorescence_erg(T_K(j))&
                              *h2_vibrational_heat_fraction(T_K(j),        &
                                       nhi(j), nh2(j), nheiS(j))
            ! CO IS THE FOURTH ABSORBER OF THIS BEAM AND IT IS RATED,
            ! but it is accumulated on its OWN row and not into
            ! rated_ph.  The reason is the beam it would be compared
            ! against.  beam_loss_ph is N_b (1 - T_line T_cont) with
            ! T_line the H2 lines and T_cont the H2O and OH continua; CO
            ! contributes to neither, because its equivalent width is
            ! inside its own shielding function and it adds no term to
            ! tau_cont (co_photodissociation.f90, and the design decision
            ! that put it there).  Counting CO's absorptions against a
            ! beam loss that carries no CO term would compare two
            ! different beams, so they are printed separately, together
            ! with their share of the beam's loss -- which IS the size of
            ! the approximation that the other three absorbers see a beam
            ! undepleted by CO.
            !
            ! The energy is charged at the CO events' OWN mean photon
            ! energy and not at the flat-band mean, because CO absorbs at
            ! the blue end of the beam (37 lines between 912.7 and
            ! 1076.1 A) and its dissociating photons are 12.87 eV against
            ! the band's 11.74.  That is what makes absorbed = heat + bond
            ! exact for this absorber.
            co_ph   = co_ph   + k_co(j)*dn_co
            co_en   = co_en   + k_co(j)*dn_co*e_co_photon_erg
            co_heat = co_heat + k_co(j)*dn_co                         &
                                *heat_per_co_dissociation()
         endif
      enddo
      bond_col = absen - heat_col
	end subroutine fuv_band_absorption_ledger
	
	!----------------------------------!
	
	subroutine chemical_rate_coefficients(T_K,                            &
	              rchiiB, rcheiiB, rcheiiiB, rec_m,                        &
	              a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, a_ion_HeITR)
	! The recombination and collisional-ionization rate coefficients the
	! ionization equilibrium is solved with. Every one of them is a function
	! of the TEMPERATURE ALONE, which is what lets them be formed before the
	! composition sweep and used unchanged by it.
	!
	! They are the same coefficients eval_cool builds on its way to the
	! cooling, from the same range routines in the same order, so a cell
	! gets the same number from either; what this routine does not do is
	! assemble the cooling, which is a contraction with the composition and
	! therefore belongs to the state the sweep RETURNS, not to the one it
	! was handed (the second eval_cool call of ionization_equilibrium).
	!
	! Blocked over cells in one parallel region, like eval_cool: every
	! quantity is cell-local, so the thread count changes neither the
	! arithmetic nor the result.

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rchiiB, rcheiiB, rcheiiiB
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) :: rec_m, aion_m
	real*8, dimension(1-Ng:N+Ng), intent(out) :: a_ion_HI, a_ion_HeI,     &
	                                             a_ion_HeII
	! The He 2^3S collisional-ionization coefficient, for the callers that
	! track the metastable.
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: a_ion_HeITR

	integer :: ib, nblk, j_lo, j_hi, ncell, i
	real*8, dimension(1-Ng:N+Ng) :: metal_col, aion_HeITR

	! FIXED BLOCKS, NOT ONE PER THREAD. The rate-coefficient loops below are
	! vectorized by gfortran -O3 into glibc's libmvec exp/log/pow (two lanes
	! per call, MEASURED: Cool_coeff.o imports _ZGVbN2v_exp, _ZGVbN2v_log,
	! _ZGVbN2vv_pow), and the vector and the scalar libm variants agree to
	! one ulp, not bitwise. Inside a block the cells are taken in lane pairs
	! from j_lo and any odd remainder by the scalar call, so WHERE a block
	! starts decides which cell meets which variant. A block count tied to
	! omp_get_max_threads() therefore made the thread count an input to the
	! last bit of every rate (6.6e-14 between 1 and 16 threads on wasp_full,
	! COST3/HYG-MAIN). The blocks are now a fixed, even length, so the lane
	! pairing is that of the serial whole-array loop whatever the thread
	! count; the threads only decide who takes which block.
	ncell = N + 2*Ng
	nblk  = max(1, (ncell + xuv_rate_block - 1)/xuv_rate_block)
	!$omp parallel do schedule(dynamic) default(shared)                   &
	!$omp    private(ib,j_lo,j_hi,i,metal_col,aion_HeITR)
	do ib = 1,nblk
	   j_lo = (1-Ng) + (ib-1)*xuv_rate_block
	   j_hi = min((1-Ng) + ib*xuv_rate_block - 1, N+Ng)
	   if (j_hi .lt. j_lo) cycle

	   call rec_HII_B_range(T_K,rchiiB,j_lo,j_hi)      ! HII
	   call rec_HeII_B_range(T_K,rcheiiB,j_lo,j_hi)    ! HeII
	   call rec_HeIII_B_range(T_K,rcheiiiB,j_lo,j_hi)  ! HeIII

	   call ion_coeff_HI_range(T_K,a_ion_HI,j_lo,j_hi)      ! HI
	   call ion_coeff_HeI_range(T_K,a_ion_HeI,j_lo,j_hi)    ! HeI
	   call ion_coeff_HeII_range(T_K,a_ion_HeII,j_lo,j_hi)  ! HeII
	   ! He 2^3S metastable, 4.8 eV threshold
	   call ci_HeI23S_range(T_K,aion_HeITR,j_lo,j_hi)
	   if (present(a_ion_HeITR))                                          &
	      a_ion_HeITR(j_lo:j_hi) = aion_HeITR(j_lo:j_hi)

	   ! The metal stages. Skipped, and left at zero, for a metals-off run:
	   ! every rate they carry multiplies a density that is identically
	   ! zero there, which is the same gate eval_cool applies.
	   rec_m(j_lo:j_hi,:)  = 0.0d0
	   aion_m(j_lo:j_hi,:) = 0.0d0
	   if (thereis_metals) then
	      do i = 1,n_mion
	         call rec_coeff_by_ion_range(i,T_K,metal_col,j_lo,j_hi)
	         rec_m(j_lo:j_hi,i)  = metal_col(j_lo:j_hi)
	         call ion_coeff_by_ion_range(i,T_K,metal_col,j_lo,j_hi)
	         aion_m(j_lo:j_hi,i) = metal_col(j_lo:j_hi)
	      enddo
	   endif

	enddo
	!$omp end parallel do

	end subroutine chemical_rate_coefficients

	!---------------------------------------------------!

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
	integer :: ib, nblk, j_lo, j_hi, ncell
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
	! Neutral helium in its ground singlet, n(1^1S) (see below)
	real*8, dimension(1-Ng:N+Ng) :: nheiS

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

	! Neutral helium in its ground singlet, n(1^1S) = n(He I) - n(2^3S).
	! The state vector's He I column CONTAINS the metastable
	! (composition.f90), and the ground-state coefficients of the cooling --
	! the 24.6 eV collisional ionization and the Cen (1992) collisional
	! excitation, both of the 1^1S term -- act on the singlet alone: the
	! metastable sits 19.8 eV up and carries its own 4.8 eV ionization and
	! its own 10830 A and singlet-conversion channels, which eval_cool_cells
	! adds from nheiTR. Charging it the ground-state coefficients as well
	! would count it twice, once in each level's channels. Same subtraction
	! as PH_heat_HHe and photoheating_of_composition make of the same column.
	! Formed here, over the whole grid, before the parallel region below.
	! Without a metastable column there is no metastable inside n(He I) and
	! the singlet is that column unchanged.
	if (present(nheiTR)) then
		nheiS = he_ground_singlet_density(nhei, nheiTR)
	else
		nheiS = nhei
	endif

	! One parallel region per call, over contiguous blocks of cells of the
	! fixed even length xuv_rate_block (see chemical_rate_coefficients for
	! why the length is fixed and not one block per thread: the libmvec lane
	! pairing of the vectorized exp/log/pow inside eval_cool_cells). The
	! blocks are disjoint and every quantity is cell-local, so the thread
	! count decides only which thread takes which block. The sub-block
	! timers of the blocks are summed, so slots 1-3, 5 and the CNO part of
	! slot 4 report thread time rather than wall time once threads are on.
	ncell = N + 2*Ng
	nblk  = max(1, (ncell + xuv_rate_block - 1)/xuv_rate_block)
	!$omp parallel do schedule(dynamic) default(shared)                 &
	!$omp    private(ib,j_lo,j_hi,ect) reduction(+:ec_t)
	do ib = 1,nblk
	   j_lo = (1-Ng) + (ib-1)*xuv_rate_block
	   j_hi = min((1-Ng) + ib*xuv_rate_block - 1, N+Ng)
	   if (j_hi .lt. j_lo) cycle
	   ect = 0.0d0
	   call eval_cool_cells(j_lo,j_hi, ect,                             &
	          T_K,nhi,nhii,nheiS,nheii,nheiii, nm, ne, beta_fs, nbar_fs,&
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
	                   T_K,nhi,nhii,nheiS,nheii,nheiii, nm, ne,      &
	                   beta_fs, nbar_fs,                             &
	                   rchiiB,rcheiiB,rcheiiiB, rec_m,               &
	                   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,        &
	                   cool, cool_chan, nheiTR, a_ion_HeITR, nmol, nox)

	! The cell-local part of eval_cool, over the cells j_lo:j_hi. Called
	! from inside the parallel region of eval_cool, once per block of
	! cells. Every array keeps the explicit shape 1-Ng:N+Ng of the whole
	! grid and only the elements j_lo:j_hi are read or written, so the
	! compiled arithmetic of a cell is what the serial whole-array form
	! produced (assumed-shape dummies moved it; Update_EXHALE_stage1.md 145).
	! ne, beta_fs and nbar_fs come in from eval_cool: they are not
	! cell-local. ect returns this block's sub-block times.

	integer, intent(in) :: j_lo, j_hi
	real*8, intent(inout) :: ect(5)
	real*8 :: ec_tstart

	integer :: i,j

	! nheiS is the He I GROUND SINGLET density n(1^1S), formed by eval_cool
	! from the summed neutral-helium column and the metastable it contains.
	! The ground-state collisional coefficients below act on it alone; the
	! metastable's own channels are added from nheiTR.
	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,           &
	                                            nheiS,nheii,nheiii
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

	! Cooling rate. The prefactor of each term is the ionization potential of
	! that stage in erg, read from the named global constants e_th_*_erg
	! (global_parameters), which are the eV thresholds of the photon grid
	! divided by erg2eV. One definition each: this assembly, the cell-by-cell
	! temperature root of the advection post-process (T_equation) and the
	! photoelectron energy h nu - e_th all charge the same energy per event.
	! e_th_HeI is the potential of the He I GROUND SINGLET and is charged to
	! nheiS, the metastable removed; the metastable's own e_th_HeTR is the
	! term below it. The 2^3S term (added only when nheiTR is supplied)
	! reproduces the Black (1981) form 6.41e-21 sqrt(T) exp(-55338/T)
	! n_e n_23S once multiplied by n_e below.
	coio(j_lo:j_hi) =  e_th_HI_erg*a_ion_HI(j_lo:j_hi)*nhi(j_lo:j_hi)  	 & ! HI
		  + e_th_HeI_erg*a_ion_HeI(j_lo:j_hi)*nheiS(j_lo:j_hi) 	 & ! HeI (1^1S)
		  + e_th_HeII_erg*a_ion_HeII(j_lo:j_hi)*nheii(j_lo:j_hi)   ! HeII
	if (present(nheiTR)) coio(j_lo:j_hi) = coio(j_lo:j_hi)                                     &
		  + e_th_HeTR_erg*aion_HeITR(j_lo:j_hi)*nheiTR(j_lo:j_hi)   ! He 2^3S
	
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
		  + coeff_coex_rate_HeI(j_lo:j_hi)*nheiS(j_lo:j_hi)    &    ! HeI (1^1S)
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
	! (coeff_coex_rate_HeI*nheiS, an excitation OUT of the ground singlet);
	! adding q13 here would double count it.
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
	! of C I, C II, N II and O I, and for the METASTABLE terms of all six
	! C/N/O fits (CHIANTI mode only; the legacy AIOLOS fits keep their own
	! constant floors). Same Lambda_eff = W_FS/ne + remainder convention as
	! Fe II above; the statistical-equilibrium solution saturates the floor
	! (n_crit,e([C II] 158um) ~ 20 cm^-3!) and adds the H-collision
	! excitation channel the electron-only coronal curve misses. beta
	! enters as A_ul -> beta*A_ul inside that solution. The metastable
	! terms of the remainder are saturated against their LTE ceiling in the
	! same coefficients (n_crit,e is 1e4-1e9 cm^-3 for those levels). N I
	! and O II have a single-level 4S* ground term, so their coefficients
	! carry the metastable saturation alone and take no line-trapping
	! argument. See cool_CI_ne_func / cooling_data/fit_fs_saturation.py.
	if (cno_chianti) then
		call cool_CI_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_CI)  = metal_col(j_lo:j_hi)
		call cool_CII_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_CII) = metal_col(j_lo:j_hi)
		call cool_NI_ne_range(T_K, ne, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_NI)  = metal_col(j_lo:j_hi)
		call cool_NII_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_NII) = metal_col(j_lo:j_hi)
		call cool_OI_ne_range(T_K, ne, nhi, beta_fs, nbar_fs, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_OI)  = metal_col(j_lo:j_hi)
		call cool_OII_ne_range(T_K, ne, metal_col,j_lo,j_hi)
		c_metal(j_lo:j_hi,im_OII) = metal_col(j_lo:j_hi)
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
	! point is what docs/TO_BE_DONE.md item (G) asks for; the emission magnitudes
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
	! Diagnostic. Dump output/Heating_breakdown.txt: the volumetric heating
	! rate of each channel of heat_channel_name against radius, in cgs
	! erg cm^-3 s^-1.
	!
	! ONE STATE PER OUTPUT FILE. The channel columns are heat_channel_state,
	! the array the ionization sweep filled through heating_of_composition
	! for the state it returned. Nothing is recomputed here: no radiation
	! field, no rate, no contraction. The heat_total column is the sum of
	! the channels, and it is therefore the heat column of Hydro_ioniz.txt
	! to round-off, not to the size of one sweep's rate lag. Rebuilding the
	! rates on the written state instead would describe a different state
	! from the one whose heating the run integrated: the rates a sweep
	! contracts are those of the composition that ENTERED it (see the
	! heating assembly in ionization_equilibrium), so they cannot be
	! recovered from the state the sweep returned.
	!
	! T and n_e are columns of the state this routine is given, so that the
	! file carries the state its channels belong to. A run that writes this
	! file before any sweep gets the zeros the array was allocated with.

	integer :: j,im,ic
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) :: n_dim,T_K,ne,heat_tot
	real*8, dimension(1-Ng:N+Ng) :: nhii,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	! Molecular densities [cm^-3] (H2, H2+, H3+, HeH+); zero for an atomic run
	real*8, dimension(1-Ng:N+Ng,4) :: nmol
	character(len=800) :: head
	character(len=8)   :: colnum

	n_dim = n_in*n0
	T_K   = T_in*T0
	nhii  = f_sp_in(:,2)*n_dim
	if (thereis_He) then
		nheii  = f_sp_in(:,4)*n_dim
		nheiii = f_sp_in(:,5)*n_dim
	else
		nheii = 0.0d0; nheiii = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	nmol = 0.0d0
	if (thereis_mol) then
		nmol(:,1) = f_sp_in(:,isp_H2)  *n_dim
		nmol(:,2) = f_sp_in(:,isp_H2p) *n_dim
		nmol(:,3) = f_sp_in(:,isp_H3p) *n_dim
		nmol(:,4) = f_sp_in(:,isp_HeHp)*n_dim
	endif
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol)

	heat_tot = 0.0d0
	if (allocated(heat_channel_state)) then
		do ic = 1,n_heat_channel
			heat_tot = heat_tot + heat_channel_state(:,ic)
		enddo
	endif

	open(unit = 72, file = './output/Heating_breakdown.txt')
	write(72,'(a)') '# Volumetric heating rate in each channel [cgs erg cm^-3 s^-1] vs radius.'
	write(72,'(a)') '# The channels are those the ionization sweep deposited for the state'  &
	             // ' written beside it, so their sum is the heat column of Hydro_ioniz.txt.'
	! Header built from the channel list, so a channel cannot be added to
	! the heating and be missing from the columns.
	head = '# col1 r/Rp  col2 T[K]  col3 ne  col4 heat_total'
	do ic = 1,n_heat_channel
		write(colnum,'(i0)') ic + 4
		head = trim(head)//'  col'//trim(colnum)//' '                     &
		       //trim(heat_channel_name(ic))
	enddo
	write(72,'(a)') trim(head)
	call write_row_layout_header(72)
	if (allocated(heat_channel_state)) then
		do j = 1-Ng,N+Ng
			write(72,*) r(j), T_K(j), ne(j), heat_tot(j),                 &
			            (heat_channel_state(j,ic), ic = 1,n_heat_channel)
		enddo
	else
		do j = 1-Ng,N+Ng
			write(72,*) r(j), T_K(j), ne(j), heat_tot(j),                 &
			            (0.0d0, ic = 1,n_heat_channel)
		enddo
	endif
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
			! n_crit itself leaves the double-precision range at low
			! temperature: 1100 exp(1.2/T4) sqrt(T4) is +Inf below
			! T = 17.1 K (T4 = 1.71e-3 gives 1.2/T4 + ln(1100 sqrt(T4))
			! = 705.6, and the exponent range ends at 709.78). There
			! n_e/n_crit is zero and z is at its low-density limit, so
			! that limit is the value taken and n_crit is never formed.
			T4    = T_K(j)/1.0d4
			z     = 0.67d0 + 0.29d0
			if (T4 .gt. 1.71d-3) then
				ncrit = 1100.0d0*exp(1.2d0/T4)*sqrt(T4)
				z     = 0.67d0 + 0.29d0/(1.0d0 + ne(j)/ncrit)
			endif
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

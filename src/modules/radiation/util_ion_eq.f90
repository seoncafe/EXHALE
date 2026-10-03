   module utils_ion_eq

   use global_parameters
   use species_table, only: n_mion, n_mphot, n_melem,                 &
                            mion_isphot, mion_iphot, mion_ethr,       &
                            mion_z2, mion_elem, mion_iscool,          &
                            mion_stage, mion_name, mion_fsp,          &
                            melem_top, im_FeII, im_OI,                &
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
                             metal_photoion_sigma,                    &
                             metal_shell_photoion_sigma,              &
                             metal_shell_relaxation, mph_n_sub
   use composition, only: he_ground_singlet_density
   use charge_exchange, only: charge_exchange_heating
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
   ! Ly-alpha and two-photon continuum and is not heat.  The recombination
   ! photons absorbed on the spot by H2 take the same channel cross sections
   ! (h2_photoabsorption_cross_sections) and the same recipients.
   use h2_photo_channels, only: h2_double_fragment_kinetic_energy,        &
                                e_rad_H2_neutral, n_h2_channels,          &
                                ICH_S, ICH_D, ICH_N,                      &
                                h2_photoabsorption_cross_sections,        &
                                h2_channel_energy_recipients,             &
                                h2_double_ionization_model
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
                            iabs_H2, iabs_H2_di, dissoc_ion_per_H2p,     &
                            n_dal_E, dalgarno_energy_node_weights
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
	character(len=40), parameter :: ec_name(5) = (/ 'recombination/ionization coefficients   ', &
	   'radiative cooling of the cells          ', 'He ground-capture shares                ', &
	   'fine-structure line transfer            ', 'H3+ / molecular IR                      ' /)

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
	! Draine (2011, sect. 14.3.2, p. 149): "56% of these two-photon decays
	! produce a photon with hnu > 13.60 eV"; Osterbrock & Ferland (2006,
	! sect. 2.4, p. 30): "0.56 per radiative decay from He0 2 1S".
	!
	! The radiative rate of the 2^3S metastable to 1^1S [s^-1], Drake
	! (1971) as listed by Lampon et al. (2020, A&A 636, A13, Table 2, "A31
	! ... 1.272 x 10^-4"); the one value HeITR_coeffs and the atomic-mode
	! cascade of recombination_radiation_absorbed use.
	real*8, parameter :: A_HeI_23S_11S = 1.272d-4

	! ----- Representative photon energies of the He recombination channels -----
	! Each channel of recombination_radiation_absorbed emits at one energy, at which the
	! absorbers (H I, H2 and, above 24.6 eV, He I) compete for the photon.
	! The ground-capture continuum starts AT the He I ionization threshold,
	! so this energy is that threshold and is read from it.
	real*8, parameter :: E_gnd_HeI  = e_th_HeI ! ground-capture continuum edge
	! The two He I lines, at the NIST ASD level energies of Cool_coeff
	! (E_HeI_23S_eV = 19.8196 eV; 2^1P at 21.2180 eV).
	real*8, parameter :: E_584_HeI  = E_HeI_23S_eV + E_HeI_21P_23S_eV
	real*8, parameter :: E_19_HeI   = E_HeI_23S_eV
	! Mean energy of the 2^1S two-photon photons that lie above the H I edge:
	! the continuum is not a line, and this single energy stands for it. The
	! shape (Drake, Victor & Dalgarno 1969) is not carried in the code, only
	! the two integrals f_2q_HeI and Ee_2q_HeI over the H I window, so the
	! energy below is e_th_HI + Ee_2q_HeI and cannot be re-integrated over the
	! narrower H2 window. APPROXIMATE: the H2-ionizing photon COUNT of the
	! two-photon continuum is the integral of the same shape over e_th_H2 to
	! 20.62 eV, which is smaller than f_2q_HeI (the integral over e_th_HI to
	! 20.62 eV); f_2q_HeI is used for H2 as well, so the H2 share of this
	! channel is overestimated. R2_2q = 1.2 is the smallest H2/H I
	! cross-section ratio of the set.
	real*8, parameter :: E_2q_HeI   = e_th_HI + Ee_2q_HeI       ! 16.110 eV

	! ----- He I 584 A: a resonance line, destroyed on the spot ----- !
	! Radiative rates of 1s2p 1P1 [s^-1] to the ground (584.334 A) and to
	! 1s2s 1S0 (2.0587 um), READ from the NIST Atomic Spectra Database
	! (He I lines; accuracy class AAA, transition-probability reference
	! T8636c73): 1.7989e9 and 1.9746e6.  Draine (2011, eq. 14.18) prints
	! 1.799e9 and 1.976e6, a branching ratio of 1.098e-3.
	real*8, parameter :: A_HeI_21P_11S = 1.7989d9
	real*8, parameter :: A_HeI_21P_21S = 1.9746d6
	! The probability that one excitation of 2^1P ends in 2^1S (then in
	! the two-photon continuum) rather than in a new 584 A photon:
	! 1/912.0 (DERIVED).
	real*8, parameter :: eps_HeI_21P_21S = A_HeI_21P_21S                 &
	                                     /(A_HeI_21P_11S + A_HeI_21P_21S)
	! lambda^3 g_u A/(8 pi^(3/2) g_l) of the line [cm^4 s^-1], g_u/g_l = 3,
	! lambda = h c/E_584: divided by the thermal speed of the helium atom
	! it is the line-centre cross section of the Doppler profile
	! (He_I_584_mean_cross_section).
	real*8, parameter :: lam_584_cm   = hp_eV*c_light/E_584_HeI
	real*8, parameter :: sig0v_584    = 3.0d0*lam_584_cm**3*A_HeI_21P_11S &
	                                  /(8.0d0*pi**1.5d0)

	! ----- The ground-capture edges of the three recombinations ----- !
	real*8, parameter :: E_gnd_HI   = e_th_HI    ! H II -> H I
	real*8, parameter :: E_gnd_HeII = e_th_HeII  ! He III -> He II

	! ----- He III -> He II recombination radiation ----- !
	! He II level energies above the ground [eV], NIST ASD as distributed
	! in CHIANTI v11 he_2.elvlc (2s1/2 329179.767 cm^-1, 2p1/2 329179.299,
	! 2p3/2 329185.156; hc = 1.239841984e-4 eV cm): Ly-alpha at the
	! g-weighted 2p, and the ionization potential of n = 2, the edge of the
	! direct-capture continuum into n = 2 (13.6047 eV, above the H I edge
	! by 6.2 meV).
	real*8, parameter :: E_lya_HeII = 329183.204d0*1.239841984d-4  ! 40.8135 eV
	real*8, parameter :: E_2s_HeII  = 329179.767d0*1.239841984d-4  ! 40.8131 eV
	real*8, parameter :: E_n2_HeII  = e_th_HeII - E_2s_HeII        ! 13.6047 eV
	! Two-photon decay of He II 2s: Drake (1986, Phys. Rev. A 34, 2871,
	! eq. 27; hydrogen_n2_rates: A_2s1s) with Z = 2, the effective
	! radiative charge Z_r = 1 + m_e/(M + m_e) and the alpha-particle mass
	! M: 8.22938 x 64 x Z_r^4 (1 - m_e/M) x 0.99986007 (relativistic
	! factor) = 526.823 s^-1 (DERIVED); the nonrelativistic Nussbaumer &
	! Schmutz (1984, A&A 138, 495, eq. 6) scaling 8.2249 Z^6 R_Z/R_H gives
	! 526.61 s^-1, 0.04% lower. The spectrum is
	! their fit A(y) (eq. 2: alpha = 0.88, beta = 1.53, gamma = 0.8,
	! C = 202.0 s^-1, y = E/E_2s), integrated here in three bands, the
	! photons of each per decay (the fit's own integral normalizing) and
	! their mean energy: 13.598-15.426 eV (H I only), 15.426-24.587 eV
	! (H I, H2) and above 24.587 eV (H I, H2, He I): 0.11163 at 14.5160 eV,
	! 0.57725 at 20.0127 eV, 0.73620 at 31.1293 eV (DERIVED by quadrature
	! of the published fit, accurate to 0.1% over 0.07 < y < 0.93).
	real*8, parameter :: A_2q_HeII     = 526.823d0
	real*8, parameter :: f_2q_HeII_lo  = 0.11163d0, E_2q_HeII_lo  = 14.5160d0
	real*8, parameter :: f_2q_HeII_mid = 0.57725d0, E_2q_HeII_mid = 20.0127d0
	real*8, parameter :: f_2q_HeII_hi  = 0.73620d0, E_2q_HeII_hi  = 31.1293d0
	! 2s -> 2p ion collisions of He II (Pengelly & Seaton 1964, eq. 48:
	! H+ and He2+ impact): the separations of 2s1/2 from 2p1/2 and 2p3/2
	! [erg] (0.468 and 5.389 cm^-1, he_2.elvlc) and the reduced masses of
	! He+ with H+ and with He2+ [g].
	real*8, parameter :: dE_2s2p12_HeII = 0.468d0*hp_erg*c_light
	real*8, parameter :: dE_2s2p32_HeII = 5.389d0*hp_erg*c_light
	real*8, parameter :: mu_HeII_p   = m_p*(m_He_atom - m_e)              &
	                                   /(m_p + m_He_atom - m_e)
	real*8, parameter :: mu_HeII_He2p = (m_He_atom - 2.0d0*m_e)           &
	                                    *(m_He_atom - m_e)                 &
	                                    /(2.0d0*m_He_atom - 3.0d0*m_e)

	! ----- The channels of recombination_radiation_absorbed ----- !
	! 1 H ground capture, 2 He I ground capture, 3 He I 584 A, 4 He I
	! 19.8 eV, 5 He I two-photon (above 13.6 eV), 6 He II ground capture,
	! 7 He II Ly-alpha, 8-10 He II two-photon bands, 11 He II n = 2
	! continuum: their photon energies
	! [eV], and the cross sections of the absorbers at them [1e-18 cm^2]
	! (H I, He I ground, He II, H2; each metal ion), which depend on
	! nothing but the switches that select the cross sections and are
	! formed once (on_the_spot_cross_sections).
	integer, parameter :: n_otsp_ch = 11
	integer, parameter :: ic_gnd_HI = 1, ic_gnd_HeI = 2, ic_584_HeI = 3,   &
	                      ic_19_HeI = 4, ic_2q_HeI = 5, ic_gnd_HeII = 6,  &
	                      ic_lya_HeII = 7, ic_2q_HeII_lo = 8,             &
	                      ic_2q_HeII_mid = 9, ic_2q_HeII_hi = 10,         &
	                      ic_n2_HeII = 11
	real*8, parameter :: otsp_E(n_otsp_ch) = [ E_gnd_HI, E_gnd_HeI,        &
	        E_584_HeI, E_19_HeI, E_2q_HeI, E_gnd_HeII,                    &
	        E_lya_HeII, E_2q_HeII_lo, E_2q_HeII_mid, E_2q_HeII_hi,        &
	        E_n2_HeII ]
	real*8,  save :: otsp_sab(4,n_otsp_ch) = 0.0d0
	real*8,  save :: otsp_smet(n_mion,n_otsp_ch) = 0.0d0
	! What a metal absorption of each channel's photon does beyond the
	! outer-shell event (Cross_sections: metal_shell_relaxation), averaged
	! over the shells the photon can open with their cross sections as
	! weights: otsp_met_dbind is the energy [eV] the photoelectron and the
	! Auger electrons receive LESS than E - mion_ethr (zero where only the
	! outer shell is open), otsp_met_fmulti the fraction of the
	! absorptions that eject two or more electrons, kept only for an ion
	! whose next-but-one stage the balance carries (the neutrals of the
	! three-stage elements; every other ion goes one carried stage up).
	real*8,  save :: otsp_met_dbind(n_mion,n_otsp_ch)  = 0.0d0
	real*8,  save :: otsp_met_fmulti(n_mion,n_otsp_ch) = 0.0d0
	logical, save :: otsp_ready  = .false.
	! The H2 share of each channel split into the four final states of
	! h2_photo_channels (M, S, D, N), otsp_h2_sig(:,c) [1e-18 cm^2], which
	! sums to otsp_sab(4,c), and the heat [eV] one absorption in each
	! final state deposits at the channel energy, the photoelectron plus
	! the fragment kinetic energy of h2_channel_energy_recipients. The
	! channel selection depends on two run switches (the neutral window,
	! the double-ionization model), recorded with the table.
	real*8,  save :: otsp_h2_sig(n_h2_channels,n_otsp_ch)  = 0.0d0
	real*8,  save :: otsp_h2_heat(n_h2_channels,n_otsp_ch) = 0.0d0
	logical, save :: otsp_h2_neutral = .false.
	character(len=16), save :: otsp_h2_double = ''

	! ----- The metal photoabsorption of the photon grid, shell by shell ----- !
	! Built by metal_photoabsorption_spectral_tables on the grid e_v it was
	! handed, and rebuilt if the grid is not that one.  For photo-table
	! column k and bin i:
	!   mpa_sig_multi(i,k)  the cross section [Mb] of the absorptions that
	!                       eject two or more electrons (Auger decay);
	!   mpa_w(i,0,k)        the energy [eV] of the electrons at or below
	!                       E_sec_ion, which thermalize whole, and
	!   mpa_w(i,m,k)        m = 1..n_dal_E, the energy of those above it
	!                       on the node m of the Dalgarno energy grid
	!                       (dalgarno_energy_node_weights),
	! both summed over the photoelectron (e_v - E_th,s) and the Auger
	! electrons of every open shell, weighted by the shell's cross section
	! [Mb] and divided by e_v.  Contracted with a cell's partition
	! coefficients they give its metal photoheating and secondary
	! ionizations for any number of electrons (photoionization_field_at_
	! cell_HHe).
	real*8,  allocatable, save :: mpa_e_v(:)
	real*8,  allocatable, save :: mpa_sig_multi(:,:)
	real*8,  allocatable, save :: mpa_w(:,:,:)
	logical, save :: mpa_ready = .false.
	logical, save :: otsp_ates   = .false.
	logical, save :: otsp_metals = .false.

	! ----- He recombination coupling diagnostic -----
	! Cells in which the He recombination photons ionize H I faster than the
	! stellar field does by more than this factor. The coupling is a CORRECTION
	! to the direct photoionization integral, so it dominating by three decades
	! is a statement that the on-the-spot budget has run away -- which is what
	! the unweighted 1/n_HI of the earlier form did in a fully ionized cell.
	! Counted over cells and steps, like the positivity counters, and
	! reported at the end of the run only when it is not zero.
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
	!   9    photoelectrons of the H I / He I / H2 / metal ionizations driven
	!        by recombination radiation absorbed on the spot (H II, He II,
	!        He III; recombination_radiation_absorbed)
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
	!  20    energy defects of the charge-exchange reactions
	integer, parameter :: n_heat_channel = 20
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
	      'heat_CO_photodissoc     ', 'heat_charge_exchange    ' /)

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
	! The fuller treatment is instead the solved J_Lya(r) of lya_rt as B2's
	! field: it would also carry the internally generated Ly-alpha of
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
	! Band share the Lyman-Werner lines have taken out of the beam, summed
	! cell by cell from the top of the column (see where tr_lines is set).
	real*8  :: a_path
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
	! take out of the shared beam.  The table is built for a homogeneous
	! column.  The self-shielding factor is read at the H2 column above the
	! cell and the cell's own T and n_H, the gas above approximated by the
	! local one; the band fraction A is the sum over the cells above of the
	! same local integral that cell's own rate takes (below), so that the
	! beam loses exactly what the rates of the cells above have spent.  Both
	! exist whether or not the oxygen chemistry is on.
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
	!
	! A ALONG THE PATH, NOT A OF A UNIFORM SLAB.  The table is tabulated for
	! a slab of one (T, n_H), and the column above a cell of the molecular
	! layer is not one: MEASURED on the warm-started oxygen_chemistry state,
	! T runs from 462 to 3370 K over the column.  The H2 rates of every cell
	! (lyman_werner_band_absorption_rate_cell_mean) take sigma_pump at that
	! cell's own (T, n_H) across the cell, so the photons the lines have
	! taken out of the beam above a face are the sum over the cells above of
	! the same integral, A(N_in; T_k, n_k) - A(N_out; T_k, n_k) cell by
	! cell, and not A(N; T_j, n_j) of the whole column at the temperature of
	! the cell the beam has reached.  The continuum rates of the shared band
	! below multiply by this transmission, so the two readings are the
	! difference between a shared beam whose photons are counted once and
	! one whose continuum sees photons the lines have already spent (or
	! misses photons they have not): MEASURED on that state, the continuum
	! share of the band 0.361 with the slab reading against 0.544 with the
	! path one, and the band's photon ledger 1.8e-1 of its beam loss short
	! with the slab reading against -3.4e-3 with the path one, the residual
	! of the product-of-means discretization of the shared beam.  On an
	! isothermal column at uniform n_H the two readings are the same number.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		call calc_column_dens_one(nH2, NH2col)
		a_path = 0.0d0
		do j = N+Ng,1-Ng,-1
			f_shield(j) = h2_self_shielding_level_resolved(NH2col(j),     &
			                            T_K(j), nH_nuc(j))
			if (j .lt. N+Ng) then
				NH2col_out = NH2col(j+1)
			else
				NH2col_out = 0.0d0
			endif
			a_path = a_path                                               &
			       + h2_lw_band_photon_fraction_absorbed(NH2col(j),       &
			                 T_K(j), nH_nuc(j))                           &
			       - h2_lw_band_photon_fraction_absorbed(NH2col_out,      &
			                 T_K(j), nH_nuc(j))
			a_lines(j) = min(a_path, 1.0d0)
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

	! THE METAL PHOTOABSORPTION OF THE PHOTON GRID, SHELL BY SHELL: the
	! tables mpa_sig_multi and mpa_w (declared with their meaning at the
	! head of this module) on the current e_v.  Each shell s of each
	! photo-ionizable ion takes its own partial cross section
	! (metal_shell_photoion_sigma) and hands on a photoelectron of
	! e_v - E_th,s and, where its vacancy autoionizes, Auger electrons of
	! e_auger,s (metal_shell_relaxation); the fluorescence and the
	! ionization energy of the further electrons are not electron energy
	! and are not in mpa_w.  Formed once for a grid, inside a critical
	! region; a grid that has changed since (a new set_energy_vectors) is
	! detected by comparison and the tables are rebuilt.
	!
	! A FIRST CALL FROM A PARALLEL REGION IS SAFE because the flag is
	! published: mpa_ready is written with a sequentially consistent atomic
	! write after the tables are filled and read with a sequentially
	! consistent atomic read before them, so a thread that sees it set also
	! sees the tables the setting thread wrote.  A REBUILD is safe only
	! outside a parallel region, since it deallocates tables another thread
	! may be reading; e_v is formed only at initialization
	! (set_energy_vectors, sed_read), serially, so every rebuild is serial.
	subroutine metal_photoabsorption_spectral_tables
	integer :: k, is, i, ie
	real*8  :: E, sg, eth, eaug, pmul, efl, eim, Ee(2), wt(n_dal_E)
	logical :: current

	!$omp atomic read seq_cst
	current = mpa_ready
	if (current) current = (size(mpa_e_v) .eq. Nl)
	if (current) current = all(mpa_e_v .eq. e_v(1:Nl))
	if (current) return
	!$omp critical (metal_photoabsorption_table)
	!$omp atomic read seq_cst
	current = mpa_ready
	if (current) current = (size(mpa_e_v) .eq. Nl)
	if (current) current = all(mpa_e_v .eq. e_v(1:Nl))
	if (.not. current) then
		!$omp atomic write seq_cst
		mpa_ready = .false.
		if (allocated(mpa_e_v))       deallocate(mpa_e_v)
		if (allocated(mpa_sig_multi)) deallocate(mpa_sig_multi)
		if (allocated(mpa_w))         deallocate(mpa_w)
		allocate(mpa_e_v(Nl), mpa_sig_multi(Nl,n_mphot),                 &
		         mpa_w(Nl,0:n_dal_E,n_mphot))
		mpa_sig_multi = 0.0d0
		mpa_w         = 0.0d0
		do k = 1,n_mphot
			do i = 1,Nl
				E = e_v(i)
				do is = 0,mph_n_sub(k)
					sg = metal_shell_photoion_sigma(k, is, E)
					if (.not. (sg > 0.0d0)) cycle
					call metal_shell_relaxation(k, is, eth, eaug, pmul,  &
					                            efl, eim)
					mpa_sig_multi(i,k) = mpa_sig_multi(i,k) + sg*pmul
					Ee(1) = E - eth
					Ee(2) = eaug
					do ie = 1,2
						if (.not. (Ee(ie) > 0.0d0)) cycle
						if (Ee(ie) .gt. E_sec_ion) then
							call dalgarno_energy_node_weights(Ee(ie), wt)
							mpa_w(i,1:n_dal_E,k) = mpa_w(i,1:n_dal_E,k)    &
							                     + sg*Ee(ie)*wt/E
						else
							mpa_w(i,0,k) = mpa_w(i,0,k) + sg*Ee(ie)/E
						endif
					enddo
				enddo
			enddo
		enddo
		mpa_e_v   = e_v(1:Nl)
		!$omp atomic write seq_cst
		mpa_ready = .true.
	endif
	!$omp end critical (metal_photoabsorption_table)
	end subroutine metal_photoabsorption_spectral_tables

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
		Hea_1  = 0.0d0
		PIR_1  = 0.0d0
		q_abs  = 0.0d0

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
		! cross section by it.
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

	nhei   = 0.0d0
	nheii  = 0.0d0
	nheiii = 0.0d0
	nheiTR = 0.0d0

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
	             heat_j, chan_j, q_j, q_abs_j, Pm2_j)
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
	! The part of Pm_j whose absorptions eject two or more electrons (an
	! inner-shell vacancy that autoionizes) and that the balance carries as
	! a jump of two stages: nonzero only for the neutrals of the elements
	! with three carried stages (species_table melem_top = 2).  Every other
	! absorber goes one carried stage up whatever it ejects (Cross_sections
	! header, STAGES ABOVE THE CARRIED ONES).
	real*8, dimension(n_mion), intent(out), optional :: Pm2_j
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
	! The metal electron energy of the cell on each Dalgarno node,
	! sum_i n_i mpa_w(:,m,k(i)): what the metal photoelectrons and Auger
	! electrons hand the secondary-ionization partition.
	real*8, dimension(Nl,n_dal_E) :: wsec_m
	integer :: m
	type(photoelectron_partition_t) :: pep
	real*8 :: Psec_HI,Psec_HeI,Psec_H2,Psec_H2_di

	!----------------------------------!

	P_H2_j    = 0.0d0
	P_H2_di_j = 0.0d0
	P_H2_dd_j = 0.0d0
	P_H2_nd_j = 0.0d0
	Pm_j      = 0.0d0
	if (present(Pm2_j)) Pm2_j = 0.0d0
	call metal_photoabsorption_spectral_tables

      PIR_1   = 0.0d0
      PIR_15  = 0.0d0
      PIR_2   = 0.0d0
      PIR_TR  = 0.0d0
      q_abs   = 0.0d0

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
		tau_m = 0.0d0
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
		tau_m = 0.0d0
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
		! they all share is identical to multiplying each cross section by it.
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
		! included.  Each shell's absorption hands the electron gas its
		! photoelectron, e_v - E_th,s, and, where the vacancy autoionizes,
		! its Auger electrons (Cross_sections: metal_shell_relaxation); the
		! fluorescence leaves and the further ionization energy is spent.
		! mpa_w carries that electron energy node by node, so the heat of
		! one ion and the secondary ionizations of all of them are one
		! contraction with the cell's partition, for any number of
		! electrons.  Under "Photoelectron heating: full"
		! (photoheat_photon_fraction > 0, a comparison option that heats
		! with a fixed share of the photon) the outer-shell form is kept.
		if (photoheat_photon_fraction .gt. 0.0d0) then
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
		else
		wsec_m = 0.0d0
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			acc_mion = mpa_w(:,0,k)
			if (sec_on) then
				do m = 1,n_dal_E
					acc_mion    = acc_mion + pep%f_heat(m)*mpa_w(:,m,k)
					wsec_m(:,m) = wsec_m(:,m) + nm_j(i)*mpa_w(:,m,k)
				enddo
			else
				do m = 1,n_dal_E
					acc_mion = acc_mion + mpa_w(:,m,k)
				enddo
			endif
			h1m_loc(i) = sum(int_f*acc_mion*de_v)*1.0d-18
		enddo
		if (sec_on) then
			do m = 1,n_dal_E
				acc_secHI  = acc_secHI  + pep%c_ion_HI(m) /e_th_HI *wsec_m(:,m)
				acc_secHeI = acc_secHeI + pep%c_ion_HeI(m)/e_th_HeI*wsec_m(:,m)
				if (mol_sec) acc_secH2 = acc_secH2                            &
				                       + pep%c_ion_H2(m)/e_th_H2*wsec_m(:,m)
			enddo
		endif
		endif
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
		PIR_H2 = 0.0d0
		PIR_H2_di = 0.0d0
		PIR_H2_dd = 0.0d0
		PIR_H2_nd = 0.0d0
		if (has_h2) then
			PIR_H2    = sum(int_h2   *de_v)
			PIR_H2_di = sum(int_h2_di*de_v)
			PIR_H2_dd = sum(int_h2_dd*de_v)
			PIR_H2_nd = sum(int_h2_nd*de_v)
		endif
		do i = 1,n_mion
			if (.not. mion_isphot(i)) then
				Pm_loc(i) = 0.0d0
				cycle
			endif
			k = mion_iphot(i)
			int_m = int_f*sigma_tab(:,k)/e_v
			Pm_loc(i) = sum(int_m*de_v)*1.0d-18*erg2eV
			if (present(Pm2_j)) then
				if (mion_stage(i) .eq. 0 .and.                               &
				    melem_top(mion_elem(i)) .ge. 2)                           &
					Pm2_j(i) = sum(int_f*mpa_sig_multi(:,k)/e_v*de_v)        &
					           *1.0d-18*erg2eV
			endif
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
				     heat_of_one_H2, heat_of_one_mion, P_m2)
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
	real*8 :: Pm_row(n_mion), h1m_row(n_mion), Pm2_row(n_mion)
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
   ! The part of P_m that the balance carries as a two-stage jump
   ! (photoionization_field_at_cell_HHe, Pm2_j).
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(out), optional ::  P_m2

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
	P_m = 0.0d0
	if (present(P_m2)) P_m2 = 0.0d0

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
	!$OMP           chan_c, Pm2_row )
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
		         heat(j), chan_c, q(j), q_abs_cell(j), Pm2_row)

		P_m(j,:)  = Pm_row
		if (present(P_m2)) P_m2(j,:) = Pm2_row
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
	real*8, dimension(1-Ng:N+Ng) :: rchiiB_hrc,rcheiiB_hrc,rcheiiiB_hrc
	real*8, dimension(1-Ng:N+Ng) :: dP_HI_hrc,dP_HeI_hrc,dP_H2_hrc
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

	! Recombination radiation absorbed on the spot: the H I, He I, H2 and
	! metal ionizations by the ground-capture photons of H II, He II and
	! He III and by the He II and He III cascade photons
	! (recombination_radiation_absorbed). The RATE corrections of this
	! coupling belong to the state that entered the sweep and are not
	! rebuilt here; its photoelectron heating carries the densities of the
	! absorbers, so it is evaluated at the composition this routine was
	! given and only the heating is kept.
	if (use_h_rec_escape .or. (use_he_rec_coupling .and. thereis_He)) then
		call recombination_radiation_absorbed(T_K, nhi, nhii, nmol(:,1),   &
		                     nhei, nheii, nheiii, nheiTR, ne, nm,          &
		                     A31, q31a, q31b,                              &
		                     rchiiB_hrc, rcheiiB_hrc, rcheiiiB_hrc,        &
		                     dP_HI_hrc, dP_HeI_hrc, dP_H2_hrc, dP_m_hrc,   &
		                     heat_chan(:,9))
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
	! (md/p39_lw_cross_section_sources.md sec. 5.2). The count of such
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
	! (md/p39_lw_cross_section_sources.md).
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
	! why an atomic gas has no such term. Default ON (mol_reaction_heat in
	! parameters.f90; the key "Molecular reaction heat: False" turns it off).
	! On a column that carries no H2 where the network would form it, this
	! term is the three-body formation heat of that out-of-equilibrium
	! composition and it is large by construction: MEASURED 2026-09-21 on
	! molecular_photochem_gj1132_kzzprofile/HeH9 evaluated with an H2-free
	! column, the energy row of cells 3 to 5 stands at 0.98 of its own scale
	! with this term and at 7.3e-03 without it.
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

		! The two CO destruction channels of the one-sided CO model are formed
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

	! 20: the energy defects of the charge-exchange reactions the balance
	! rows apply (metal + H/He/metal and He <-> H), less the excitation
	! of product states that leaves as radiation; charge_exchange_heating
	! states the product states and the ledger. The helium reactant is
	! the ground singlet, nheiS_chem; ne sets the metastable populations.
	call charge_exchange_heating(T_K, ne, nhi, nhii, nheiS_chem, nheii,  &
	                             nheiii, nm, heat_chan(:,20))
	heat = heat + heat_chan(:,20)

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
	! routine existed the scanner had to record write_output.f90 as an
	! exception. It is one file, one set of imports and one place those
	! energies are read.
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
            ! md/p39_lw_cross_section_sources.md):
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
	! and densities, with TWO exceptions formed here for the whole grid
	! before the cell loop: the fine-structure line transfer, whose escape
	! probabilities are integrals of the columns above and below the cell,
	! and the ground-capture weight of the He II recombination, which needs
	! the cell's width and its absorber densities (ground_capture_escape_weights).
	! The cell-local remainder runs in eval_cool_cells over contiguous blocks
	! of cells inside a single OpenMP parallel region. Each cell therefore
	! gets the same arithmetic as a serial sweep, and the result is bitwise
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
	! Neutral helium in its ground singlet, n(1^1S), and the 2^3S metastable
	! (zero where the caller carries none)
	real*8, dimension(1-Ng:N+Ng) :: nheiS, nheiTR_c
	! Ground-capture escape weights of the H II, He II and He III
	! recombinations (ground_capture_escape_weights)
	real*8, dimension(1-Ng:N+Ng) :: y_HI, y_gnd, y_HeII
	real*8, dimension(1-Ng:N+Ng) :: nh2_c

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
	! units as `cool`), in the n_cool_chan columns of
	! radiative_cooling_of_cell (Cool_coeff): 1 recombination, 2 collisional
	! ionization, 3-5 collisional excitation of H I (the Lyman-alpha-dominated
	! H-line cooling), He I (with the He 2^3S metastable's own channels) and
	! He II, 6 free-free (incl. metal-ion charges); 7 H3+ infrared cooling (0
	! unless the caller supplies nmol); 8-10 the molecular infrared bands under
	! `Molecular IR bands` -- H2 lines (needs nmol), H2O and CO bands (need
	! nox) -- each the NET rate, emission minus absorption of the field from
	! below, so a column is negative wherever that channel heats; 10+i the
	! line cooling of metal ion i (0 for non-coolant ions). The columns are an
	! exact decomposition of `cool`: their sum is `cool` to round-off.
	real*8, dimension(1-Ng:N+Ng,n_cool_chan),intent(out),optional :: cool_chan

	! He 2^3S metastable density [cm^-3], present only for the triplet-tracking
	! callers. When supplied it adds the collisional channels of the 2^3S
	! state (radiative_cooling_of_cell). a_ion_HeITR returns the He(2^3S)
	! collisional-ionization rate coefficient [cm^3 s^-1] for the ionization
	! equations, mirroring a_ion_HI/HeI/HeII.
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
	! solves explicitly: beta is the line-center escape probability of each
	! line, from the columns above and below the cell, and nbar_fs carries
	! the thermal infrared field of the lower atmosphere the same lines
	! absorb ("Base IR field", off by default) (Cool_coeff.f90:
	! fine_structure_line_transfer). Every other metal ion keeps the
	! optically thin limit; see the scope note there.
	beta_fs = 1.0d0
	nbar_fs = 0.0d0
	!$ if (ec_prof_on) ec_t0 = omp_get_wtime()
	if (thereis_metals) call fine_structure_line_transfer(T_K, nm,      &
	                                                  beta_fs, nbar_fs)
	!$ if (ec_prof_on) ec_t(4) = ec_t(4) + (omp_get_wtime() - ec_t0)

	! Neutral helium in its ground singlet, n(1^1S) = n(He I) - n(2^3S).
	! The state vector's He I column CONTAINS the metastable
	! (composition.f90), and the ground-state coefficients of the cooling --
	! the 24.6 eV collisional ionization and the collisional excitation out
	! of 1^1S -- act on the singlet alone: the metastable sits 19.8 eV up and
	! carries its own channels (radiative_cooling_of_cell). Charging it the
	! ground-state coefficients as well would count it twice. Same
	! subtraction as PH_heat_HHe and photoheating_of_composition make of the
	! same column. Without a metastable column there is no metastable inside
	! n(He I) and the singlet is that column unchanged.
	if (present(nheiTR)) then
		nheiS    = he_ground_singlet_density(nhei, nheiTR)
		nheiTR_c = nheiTR
	else
		nheiS    = nhei
		nheiTR_c = 0.0d0
	endif

	! The ground-capture escape weights of the recombinations the balance
	! runs on, from the same densities recombination_radiation_absorbed is
	! handed (the summed He I column and the H2 density), so the cooling
	! charges the captures the balance performs. They do not depend on T,
	! so a caller that varies T at fixed densities (the semi-implicit
	! energy update) sees them frozen.
	y_HI   = 0.0d0
	y_gnd  = 0.0d0
	y_HeII = 0.0d0
	!$ if (ec_prof_on) ec_t0 = omp_get_wtime()
	if (use_h_rec_escape .or. (use_he_rec_coupling .and. thereis_He)) then
		if (present(nmol)) then
			nh2_c = nmol(:,1)
		else
			nh2_c = 0.0d0
		endif
		call ground_capture_escape_weights(nhi, nh2_c, nhei, nheii, nm,  &
		                                   y_HI, y_gnd, y_HeII)
	endif
	!$ if (ec_prof_on) ec_t(3) = ec_t(3) + (omp_get_wtime() - ec_t0)

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
	          T_K,nhi,nhii,nheiS,nheiTR_c,nheii,nheiii, nm, ne,         &
	          y_HI, y_gnd, y_HeII, beta_fs, nbar_fs,                    &
	          rchiiB,rcheiiB,rcheiiiB, rec_m,                           &
	          a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,                    &
	          cool, cool_chan, a_ion_HeITR, nmol, nox)
	   ec_t = ec_t + ect
	enddo
	!$omp end parallel do

	! End of subroutine
	end subroutine eval_cool

	!---------------------------------------------------!

	subroutine eval_cool_cells(j_lo,j_hi, ect,                       &
	                   T_K,nhi,nhii,nheiS,nheiTR,nheii,nheiii, nm,   &
	                   ne, y_HI, y_gnd, y_HeII, beta_fs, nbar_fs,    &
	                   rchiiB,rcheiiB,rcheiiiB, rec_m,               &
	                   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,        &
	                   cool, cool_chan, a_ion_HeITR, nmol, nox)

	! The cell-local part of eval_cool, over the cells j_lo:j_hi. Called
	! from inside the parallel region of eval_cool, once per block of
	! cells. Every array keeps the explicit shape 1-Ng:N+Ng of the whole
	! grid and only the elements j_lo:j_hi are read or written.
	! ne, the escape weights, beta_fs and nbar_fs come in from eval_cool:
	! they are not
	! cell-local. ect returns this block's sub-block times.

	integer, intent(in) :: j_lo, j_hi
	real*8, intent(inout) :: ect(5)
	real*8 :: ec_tstart

	integer :: i,j

	! nheiS is the He I GROUND SINGLET density n(1^1S), nheiTR the 2^3S
	! metastable (zero where the caller carries none).
	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,           &
	                                            nheiS,nheiTR,       &
	                                            nheii,nheiii
	! Metal ion densities (canonical species_table order)
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm

	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

	real*8, dimension(1-Ng:N+Ng,n_fsline),intent(in) :: beta_fs, nbar_fs
	real*8, dimension(1-Ng:N+Ng),intent(in) :: ne        ! Electron density
	! Ground-capture escape weights of H II, He II and He III
	real*8, dimension(1-Ng:N+Ng),intent(in) :: y_HI, y_gnd, y_HeII

   real*8, dimension(1-Ng:N+Ng) :: metal_col            ! dispatcher scratch
   real*8, dimension(1-Ng:N+Ng) :: cool_H3p             ! H3+ infrared cooling
   real*8, dimension(1-Ng:N+Ng) :: cool_H2, cool_H2O, cool_CO ! molecular bands
   real*8 :: w_ir                                       ! incident-field dilution
   ! One cell's channels, and the cell's metal densities and fine-structure
   ! transfer as contiguous vectors
   real*8 :: chan(n_cool_chan), nm_cell(n_mion)
   real*8 :: beta_cell(n_fsline), nbar_cell(n_fsline)

   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(inout) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB
   ! Metal recombination rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(inout) :: rec_m

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(inout) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII
   ! Metal collisional ionization rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(inout) :: aion_m

	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(inout) ::  cool

	! Optional cooling breakdown (see eval_cool)
	real*8, dimension(1-Ng:N+Ng,n_cool_chan),intent(inout),optional :: cool_chan

	! He(2^3S) collisional-ionization rate coefficient returned to the
	! triplet-tracking callers (see eval_cool)
	real*8, dimension(1-Ng:N+Ng),intent(inout),optional :: a_ion_HeITR

	! Molecular and oxygen-carrier densities (see eval_cool)
	real*8, dimension(1-Ng:N+Ng,4),intent(in),optional :: nmol
	real*8, dimension(1-Ng:N+Ng,3),intent(in),optional :: nox

	!-- Rate coefficients of the ionization balance --!
	! The recombination and collisional-ionization coefficients, the same
	! range routines chemical_rate_coefficients calls.
	!$ if (ec_prof_on) ec_tstart = omp_get_wtime()
	call rec_HII_B_range(T_K,rchiiB,j_lo,j_hi)      ! HII
	call rec_HeII_B_range(T_K,rcheiiB,j_lo,j_hi)    ! HeII
	call rec_HeIII_B_range(T_K,rcheiiiB,j_lo,j_hi)  ! HeIII
	call ion_coeff_HI_range(T_K,a_ion_HI,j_lo,j_hi)      ! HI
	call ion_coeff_HeI_range(T_K,a_ion_HeI,j_lo,j_hi)    ! HeI
	call ion_coeff_HeII_range(T_K,a_ion_HeII,j_lo,j_hi)  ! HeII
	if (present(a_ion_HeITR)) call ci_HeI23S_range(T_K,a_ion_HeITR,j_lo,j_hi)

	! THE METAL BLOCKS RUN ONLY WHEN THE RUN CARRIES METALS: with every
	! metal density identically zero each skipped term is a product with a
	! zero density, and the rates a caller could read are set to zero.
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
	!$ if (ec_prof_on) then
	!$    ect(1) = ect(1) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif

	!-- H3+ infrared cooling (molecular layer) --!

	! Optically thin rotational-vibrational emission of H3+, Miller et al.
	! (2013) LTE emission per molecule with their Table-6 non-LTE departure
	! factor s(T, n_H2) (h3p_cooling module -- ONE definition, shared with
	! every caller of eval_cool). Inside a molecular base at T ~ 1e3 K this is
	! the dominant coolant: the atomic channels are all exponentially
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
	! front: the H2 quadrupole plus magnetic dipole line spectrum (Roueff et
	! al. 2019) and the H2O and CO vibration-rotation bands (HITEMP through
	! the Photochem k-coefficients). Each is the NET rate -- LTE emission
	! minus absorption of the diluted B_nu(T0) the lower atmosphere presents
	! -- so each vanishes at its own radiative equilibrium temperature
	! instead of running the layer down to nothing.
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
	!$ if (ec_prof_on) then
	!$    ect(5) = ect(5) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
	!$ endif

	!-- The atomic, ionic and metal channels, cell by cell --!
	! radiative_cooling_of_cell (Cool_coeff) is the ONE assembly of them,
	! the same one the post-process temperature root balances; the total is
	! the sum of the channels, so the breakdown reproduces it to round-off.
	do j = j_lo,j_hi
		nm_cell   = nm(j,:)
		beta_cell = beta_fs(j,:)
		nbar_cell = nbar_fs(j,:)
		call radiative_cooling_of_cell(T_K(j), ne(j), nhi(j), nhii(j),     &
		        nheiS(j), nheiTR(j), nheii(j), nheiii(j), y_HI(j),         &
		        y_gnd(j), y_HeII(j), nm_cell, beta_cell, nbar_cell, chan)
		chan(7)  = cool_H3p(j)
		chan(8)  = cool_H2(j)
		chan(9)  = cool_H2O(j)
		chan(10) = cool_CO(j)
		cool(j)  = sum(chan)
		if (present(cool_chan)) cool_chan(j,:) = chan
	enddo

	! End of subroutine
	!$ if (ec_prof_on) then
	!$    ect(2) = ect(2) + (omp_get_wtime() - ec_tstart); ec_tstart = omp_get_wtime()
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
	write(71,'(a)') '#   with the He 2^3S tracked, coex_HeI carries the net'     &
	             // ' 1^1S <-> 2^3S exchange, negative where the superelastic'  &
	             // ' collisions of the metastable heat the gas.'
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

	! Various coefficients for HeI triplet chemistry. The electron-impact
	! rates are those of Cool_coeff.f90 (their block header). q13 is the
	! rate at which a 1^1S atom is put into 2^3S by electron impact: the
	! direct excitation 1^1S -> 2^3S PLUS the excitation of every higher
	! triplet level, which cascades into 2^3S (Cool_coeff:
	! excitation_rate_HeI_11S_triplets; 0.10 of the direct rate at 1e4 K,
	! 0.33 at 2e4 K). q31g is the detailed-balance reverse of the DIRECT
	! excitation, 2^3S -> 1^1S; q31a excites 2^3S -> 2^1S, and q31b
	! 2^3S -> 2^1P plus every singlet level above it (Cool_coeff:
	! excitation_rate_HeI_23S_singlets_n3; 0.023 of q31a + q31b at 1e4 K,
	! 0.09 at 2e4 K), whose cascades end in 1^1S as 2^1P's does.
	! A31 = 1.272e-4 s^-1 is the 2^3S -> 1^1S magnetic-dipole decay rate
	! (Drake 1971, as used by Oklopcic & Hirata 2018). rcheiTR is the
	! capture into the triplets (all of which end in 2^3S) and rcheii the
	! case-B capture into the singlets, i.e. into the excited singlets
	! (Cool_coeff: alpha_rec_HeII_23S, alpha_rec_HeII_excited_singlets);
	! the ground capture that escapes the cell is added to rcheii by
	! recombination_radiation_absorbed. The energy of each excitation is
	! charged by radiative_cooling_of_cell from the direct rate and the
	! He I excitation sum, not from q13 here.
	subroutine HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31g,q31a,q31b,Q31)
	
	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

	real*8, intent(out) :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiTR,rcheii,   &
								   q13,q31g,q31a,q31b,Q31

	call rec_HeII_23S(T_K,rcheiTR)
	rcheii = alpha_rec_HeII_into_singlets(T_K, 0.0d0)
	call coex_HeI_1S_23S(T_K,q13)
	q13 = q13 + excitation_rate_HeI_11S_triplets(T_K)
	call deexc_HeI_23S_1S(T_K,q31g)
	call coex_HeI_23S_21S(T_K,q31a)
	call coex_HeI_23S_21P(T_K,q31b)
	q31b = q31b + excitation_rate_HeI_23S_singlets_n3(T_K)
	A31 = A_HeI_23S_11S

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
	! 1 - exp(-tau), for a cell optical depth tau >= 0. Formed as
	! tau (1 - exp(-tau))/tau with the quotient from one_minus_exp_over_x
	! (utilities.f90): 1 - exp(-tau) itself cancels as tau -> 0, with a
	! relative error of about eps/tau and a staircase in tau of the same
	! size, while the quotient is evaluated without cancellation at every
	! tau and the product with tau costs one rounding. Returns 0 for a
	! non-positive tau, which the ionization solve can hand it as a small
	! negative trace density.
	pure function absorbed_photon_fraction(tau) result(f_abs)
	real*8, intent(in) :: tau
	real*8 :: f_abs
	if (tau .le. 0.0d0) then
		f_abs = 0.0d0
	else
		f_abs = tau*one_minus_exp_over_x(tau)
	endif
	end function absorbed_photon_fraction

	! ------------------------------------------------------------- !

	! THE GROUND-CAPTURE ESCAPE WEIGHTS of the three recombinations, cell by
	! cell. A capture into the ground level of the recombined species r
	! (H I, He I, He II) emits a photon at the ionization edge of r
	! (E_gnd_HI, E_gnd_HeI, E_gnd_HeII). The cell keeps the fraction
	! f = 1 - exp(-tau) of those photons, tau = dl sum_s n_s sigma_s(E) the
	! cell optical depth at that energy of every absorber (H I, He I, He II,
	! H2 and the photo-ionizable metal ions; the cross sections are zero
	! below each threshold), and shares what it keeps in the ratio of the
	! absorption coefficients. The fraction NOT re-absorbed by r,
	!     y_r = 1 - n_r sigma_r(E) dl (1 - exp(-tau))/tau ,
	! is the weight of the ground capture in the net recombination
	! coefficient of the balance (Cool_coeff: alpha_rec_HII_net,
	! alpha_rec_HeII_net, alpha_rec_HeIII_net) and in its cooling: y_r -> 0
	! is case B on the spot, y_r -> 1 is case A. ONE definition:
	! recombination_radiation_absorbed, eval_cool and the advection
	! post-process all call this routine. A weight whose switch is off is
	! zero (case B): y_HI under "H_rec_escape", y_gnd and y_HeII under
	! "He_rec_coupling". Cross sections in 1e-18 cm^2 (dl carries the unit).
	! VALIDITY: the local closure of recombination_radiation_absorbed.
	subroutine ground_capture_escape_weights(nhi, nh2, nhei, nheii, nm,   &
	                                         y_HI, y_gnd, y_HeII)
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi, nh2, nhei, nheii
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: y_HI, y_gnd, y_HeII
	integer, parameter :: ich(3) = [ic_gnd_HI, ic_gnd_HeI, ic_gnd_HeII]
	real*8  :: kap, dl, g_tau
	logical :: on(3)
	integer :: j, ic, im

	on  = [use_h_rec_escape,                                              &
	       use_he_rec_coupling .and. thereis_He,                          &
	       use_he_rec_coupling .and. thereis_He]
	call on_the_spot_cross_sections
	y_HI   = 0.0d0
	y_gnd  = 0.0d0
	y_HeII = 0.0d0
	!$omp parallel do default(shared) schedule(static)                    &
	!$omp    private(j, ic, im, kap, dl, g_tau)
	do j = 1-Ng,N+Ng
		dl = dr_j(j)*R0*1.0d-18
		do ic = 1,3
			if (.not. on(ic)) cycle
			kap = nhi(j)*otsp_sab(1,ich(ic)) + nhei(j)*otsp_sab(2,ich(ic))  &
			    + nheii(j)*otsp_sab(3,ich(ic)) + nh2(j)*otsp_sab(4,ich(ic))
			if (thereis_metals) then
				do im = 1,n_mion
					kap = kap + nm(j,im)*otsp_smet(im,ich(ic))
				enddo
			endif
			g_tau = absorbed_fraction_over_tau(kap*dl)*dl
			select case (ic)
				case (1); y_HI(j)   = 1.0d0 - nhi(j)*otsp_sab(1,ic_gnd_HI)*g_tau
				case (2); y_gnd(j)  = 1.0d0 - nhei(j)*otsp_sab(2,ic_gnd_HeI)*g_tau
				case (3); y_HeII(j) = 1.0d0 - nheii(j)                     &
				                      *otsp_sab(3,ic_gnd_HeII)*g_tau
			end select
		enddo
	enddo
	!$omp end parallel do

	end subroutine ground_capture_escape_weights

	! ------------------------------------------------------------- !

	! (1 - exp(-tau))/tau for tau >= 0, the fraction of a cell's photons the
	! cell absorbs divided by its optical depth: 1 at tau -> 0, 1/tau at
	! tau -> infinity, with no division singularity. Evaluated by
	! one_minus_exp_over_x (utilities.f90), which avoids the cancellation
	! of the closed form below tau = 1e-2 (the former series here switched
	! at 1e-4, where the closed form still carried a relative error of
	! 2.2e-12 and the three-term series one of 4e-14). The rate at which
	! one absorber of cross section sigma takes photons emitted at P per
	! unit volume is P sigma dl times this.
	pure function absorbed_fraction_over_tau(tau) result(g)
	real*8, intent(in) :: tau
	real*8 :: g
	g = one_minus_exp_over_x(tau)
	end function absorbed_fraction_over_tau

	! ------------------------------------------------------------- !

	! Cross sections [1e-18 cm^2] of the on-the-spot absorbers at the
	! channel energies otsp_E: otsp_sab(1..4,c) for H I, He I (ground),
	! He II and H2, otsp_smet(i,c) for metal ion i (zero for an ion that is
	! not photo-ionizable). Each is zero below its own threshold, so an
	! absorber that cannot take a photon drops out of that channel without
	! a test. Formed on the first call and again only if a switch that
	! selects a cross section ("ATES photoionization rate", the metals, the
	! H2 neutral window and double-ionization model) has changed. The fill
	! is a critical region that re-tests the flag, and the flag is published
	! with sequentially consistent atomic accesses (as in
	! metal_photoabsorption_spectral_tables), so a first call from inside a
	! parallel region is safe; both callers run it serially before their
	! own parallel loops.
	subroutine on_the_spot_cross_sections
	integer :: ic, im, k, is, ich
	real*8  :: sg, eth, eaug, pmul, efl, eim
	real*8  :: e_res, e_ele, e_frg, e_rad
	logical :: current
	!$omp atomic read seq_cst
	current = otsp_ready
	if (current .and. (otsp_ates .eqv. ates_photoion_rate)              &
	    .and. (otsp_metals .eqv. thereis_metals)                        &
	    .and. (otsp_h2_neutral .eqv. h2_neutral_dissociation)           &
	    .and. otsp_h2_double .eq. h2_double_ionization_model) return
	!$omp critical (on_the_spot_cross_section_table)
	!$omp atomic read seq_cst
	current = otsp_ready
	if (.not. (current .and. (otsp_ates .eqv. ates_photoion_rate)       &
	           .and. (otsp_metals .eqv. thereis_metals)                 &
	           .and. (otsp_h2_neutral .eqv. h2_neutral_dissociation)    &
	           .and. otsp_h2_double .eq. h2_double_ionization_model)) then
	do ic = 1,n_otsp_ch
		otsp_sab(1,ic) = sigma(otsp_E(ic), 1.0d0, e_th_HI)
		otsp_sab(2,ic) = sigma_HeI(otsp_E(ic))
		otsp_sab(3,ic) = sigma(otsp_E(ic), 2.0d0, e_th_HeII)
		otsp_sab(4,ic) = sigma_H2(otsp_E(ic))
		call h2_photoabsorption_cross_sections(otsp_E(ic), otsp_sab(4,ic), &
		                         h2_neutral_dissociation, otsp_h2_sig(:,ic))
		do ich = 1,n_h2_channels
			call h2_channel_energy_recipients(ich, otsp_E(ic), e_res,     &
			                                  e_ele, e_frg, e_rad)
			otsp_h2_heat(ich,ic) = e_ele + e_frg
		enddo
		otsp_smet(:,ic)       = 0.0d0
		otsp_met_dbind(:,ic)  = 0.0d0
		otsp_met_fmulti(:,ic) = 0.0d0
		if (thereis_metals) then
			do im = 1,n_mion
				if (.not. mion_isphot(im)) cycle
				k = mion_iphot(im)
				otsp_smet(im,ic) = metal_photoion_sigma(k, otsp_E(ic))
				if (.not. (otsp_smet(im,ic) > 0.0d0)) cycle
				! The subshells this photon opens (the outer shell
				! contributes neither term).
				do is = 1,mph_n_sub(k)
					sg = metal_shell_photoion_sigma(k, is, otsp_E(ic))
					if (.not. (sg > 0.0d0)) cycle
					call metal_shell_relaxation(k, is, eth, eaug, pmul,  &
					                            efl, eim)
					otsp_met_dbind(im,ic) = otsp_met_dbind(im,ic)        &
					     + sg*(eth - eaug - mion_ethr(im))/otsp_smet(im,ic)
					if (mion_stage(im) .eq. 0 .and.                      &
					    melem_top(mion_elem(im)) .ge. 2)                  &
						otsp_met_fmulti(im,ic) = otsp_met_fmulti(im,ic)  &
						     + sg*pmul/otsp_smet(im,ic)
				enddo
			enddo
		endif
	enddo
	otsp_ates   = ates_photoion_rate
	otsp_metals = thereis_metals
	otsp_h2_neutral = h2_neutral_dissociation
	otsp_h2_double  = h2_double_ionization_model
	!$omp atomic write seq_cst
	otsp_ready  = .true.
	endif
	!$omp end critical (on_the_spot_cross_section_table)
	end subroutine on_the_spot_cross_sections

	! ------------------------------------------------------------- !

	! THE HE I 584 A LINE AS AN ABSORBER: the mean cross section [1e-18
	! cm^2] of one He I atom for a 584 A photon, the quantity Wood, Mathis
	! & Ercolano (2004, MNRAS 348, 1337, sect. 4.1) set the mean free path
	! between two scatterings with, l_scat = [n(He0) a_nu(Ly-alpha)]^-1:
	! the Doppler profile of their eq. (15),
	!    a_nu = pi^(1/2) e^2 f exp[-(dnu/dnu_D)^2]/(m_e c dnu_D) ,
	! with exp[-(dnu/dnu_D)^2] replaced by its mean over the Maxwellian of
	! the atom that emitted the photon, 2^(-3/2) (their eq. 16).  Written
	! through the radiative rate, pi^(1/2) e^2 f/(m_e c dnu_D) =
	! 3 lambda^3 A/(8 pi^(3/2) v_th), v_th = (2 k T/m_He)^(1/2).  At 1e4 K
	! 1.33e-14 cm^2 (DERIVED; with their f = 0.29, 1.39e-14).
	elemental double precision function He_I_584_mean_cross_section(T)   &
	                                   result(s)
	real*8, intent(in) :: T
	s = 2.0d0**(-1.5d0)*sig0v_584                                         &
	    /sqrt(2.0d0*kb_erg*max(T, 1.0d0)/m_He_atom)*1.0d18
	end function He_I_584_mean_cross_section

	! ------------------------------------------------------------- !

	! ONE CHANNEL OF ONE CELL: P photons per unit volume and time of
	! channel ic (energy otsp_E(ic) plus the capture kinetic energy Ek [eV]
	! for a continuum), shared among the absorbers of the cell (densities
	! nhi, nhei, nheii, nh2, nm; cell width dl in cm times 1e-18). i_self
	! is the recombined species of a ground-capture channel (1 H I, 2 He I,
	! 3 He II), whose share is the on-the-spot cancellation and is skipped,
	! 0 for every other channel. Adds the rate per absorber particle
	! [s^-1] of H I, He I, H2 and each metal ion, the H2 rate split into
	! its four final states (dP_H2_chan, M S D N of h2_photo_channels,
	! summing to the dP_H2 increment), the part of each metal rate that
	! ejects two or more electrons (dP_m2, nonzero only where the balance
	! carries that jump) and the photoelectron heat [erg cm^-3 s^-1]: a
	! metal absorption heats with the photoelectron and the Auger
	! electrons of the shells the photon opens, E - mion_ethr -
	! otsp_met_dbind (on_the_spot_cross_sections); an H2 absorption with
	! the photoelectron(s) and the fragment kinetic energy of its final
	! state (otsp_h2_heat), the formation energy the products store and
	! the n = 2 excitation the neutral window radiates not being heat.
	!
	! A RESONANCE LINE (the He I 584 A photons: kap_line present, the line
	! opacity n(He I) sigma_bar in the units of kap) is scattered by He I
	! as it crosses the cell, and each scattering either re-emits it (the
	! fraction 1 - eps_line) or converts it to the 2^1S two-photon
	! continuum (eps_line, the 2.06 um branch).  One flight through the
	! cell, of total depth tau = (kap + kap_line) dl, ends in a continuum
	! absorption with probability c = q a, in a line scattering with
	! l = (1 - q) a, or leaves the cell with e = 1 - a, where
	! a = 1 - exp(-tau) and q = kap/(kap + kap_line).  Repeating the flight
	! after every scattering that re-emits, the photon is absorbed in the
	! continuum with probability c/D, converted with l eps_line/D and
	! escapes with e/D, D = e + c + l eps_line (the three sum to one; D
	! is formed as that sum, which has no cancellation).  With no He I
	! (kap_line = 0) this is the continuum channel exactly; for a cell
	! thick in both the line and the continuum (e -> 0) the absorbed share
	! is q/(q + eps_line (1 - q)), Wood, Mathis & Ercolano's eq. (17),
	!   P(H_OTS) = n(H0) a(H0)/[n(H0) a(H0) + n(He0) a_nu(He0)/914] ,
	! with every continuum absorber of the cell in place of H0 alone.  The
	! photons converted are returned in P_conv [cm^-3 s^-1], for the
	! caller to put through the 2^1S two-photon channel.
	pure subroutine absorb_recombination_channel(ic, P, Ek, i_self,       &
	                   nhi, nhei, nheii, nh2, nm, dl,                     &
	                   dP_HI, dP_HeI, dP_H2, dP_H2_chan, dP_m, dP_m2,     &
	                   heat, kap_line, eps_line, P_conv)
	integer, intent(in)    :: ic, i_self
	real*8,  intent(in)    :: P, Ek, nhi, nhei, nheii, nh2, nm(n_mion), dl
	real*8,  intent(inout) :: dP_HI, dP_HeI, dP_H2,                       &
	                          dP_H2_chan(n_h2_channels), dP_m(n_mion),    &
	                          dP_m2(n_mion), heat
	real*8,  intent(in),    optional :: kap_line, eps_line
	real*8,  intent(inout), optional :: P_conv
	real*8  :: kap, g_dl, Eph, rate_s, heat_c, rate_ch(n_h2_channels)
	real*8  :: tau_t, a_t, e_t, q_c, q_l, D_t
	integer :: im
	if (.not. (P > 0.0d0)) return
	kap = nhi*otsp_sab(1,ic) + nhei*otsp_sab(2,ic) + nheii*otsp_sab(3,ic) &
	    + nh2*otsp_sab(4,ic)
	if (thereis_metals) then
		do im = 1,n_mion
			kap = kap + nm(im)*otsp_smet(im,ic)
		enddo
	endif
	g_dl   = absorbed_fraction_over_tau(kap*dl)*dl
	if (present(kap_line)) then
		if (kap_line > 0.0d0) then
			tau_t = (kap + kap_line)*dl
			a_t   = absorbed_photon_fraction(tau_t)
			if (tau_t .lt. 1.0d-4) then
				e_t = 1.0d0 - a_t
			else
				e_t = exp(-tau_t)
			endif
			q_c  = kap/(kap + kap_line)
			q_l  = kap_line/(kap + kap_line)
			D_t  = e_t + q_c*a_t + q_l*a_t*eps_line
			g_dl = absorbed_fraction_over_tau(tau_t)*dl/D_t
			P_conv = P_conv + P*q_l*a_t*eps_line/D_t
		endif
	endif
	Eph    = otsp_E(ic) + Ek
	heat_c = 0.0d0
	if (i_self .ne. 1 .and. otsp_sab(1,ic) > 0.0d0) then
		rate_s = P*otsp_sab(1,ic)*g_dl
		dP_HI  = dP_HI + rate_s
		heat_c = heat_c + rate_s*nhi*(Eph - e_th_HI)
	endif
	if (i_self .ne. 2 .and. otsp_sab(2,ic) > 0.0d0) then
		rate_s = P*otsp_sab(2,ic)*g_dl
		dP_HeI = dP_HeI + rate_s
		heat_c = heat_c + rate_s*nhei*(Eph - e_th_HeI)
	endif
	if (otsp_sab(4,ic) > 0.0d0) then
		! Each final state is charged the heat of h2_channel_energy_
		! recipients at the channel energy E_c, and the capture kinetic
		! energy Ek on top of it in every final state: in M, S and D it is
		! more photoelectron energy (the Coulomb-explosion release of D is
		! fixed by the vertical threshold), in N more recoil of the two
		! atoms. Every channel with Ek > 0 lies outside the neutral window,
		! so the last case does not arise with the present channels.
		rate_s     = P*otsp_sab(4,ic)*g_dl
		rate_ch    = P*otsp_h2_sig(:,ic)*g_dl
		dP_H2      = dP_H2 + rate_s
		dP_H2_chan = dP_H2_chan + rate_ch
		heat_c     = heat_c + nh2*(sum(rate_ch*otsp_h2_heat(:,ic))       &
		                           + rate_s*Ek)
	endif
	if (thereis_metals) then
		do im = 1,n_mion
			if (.not. (otsp_smet(im,ic) > 0.0d0)) cycle
			rate_s     = P*otsp_smet(im,ic)*g_dl
			dP_m(im)   = dP_m(im) + rate_s
			dP_m2(im)  = dP_m2(im) + rate_s*otsp_met_fmulti(im,ic)
			heat_c     = heat_c + rate_s*nm(im)                           &
			           *(Eph - mion_ethr(im) - otsp_met_dbind(im,ic))
		enddo
	endif
	heat = heat + heat_c/erg2eV
	end subroutine absorb_recombination_channel

	! ------------------------------------------------------------- !

	! RECOMBINATION RADIATION ABSORBED ON THE SPOT: the ionizing photons the
	! recombinations of H II, He II and He III emit, absorbed in the cell
	! they are emitted in. Given the pre-solve (lagged) densities, returns
	! the recombination coefficients the ionization balance should use
	! (rchiiB_new, rcheiiB_new, rcheiiiB_new [cm^3 s^-1]), the extra
	! photoionization rates [s^-1] of H I, He I (ground), H2 and every
	! metal ion (dP_HI, dP_HeI, dP_H2, dP_m), the final-state subsets of
	! the H2 rate (dP_H2_di, dP_H2_dd, dP_H2_nd), and the photoelectron
	! heating [erg cm^-3 s^-1] of the composition passed in (dheat). With
	! both switches off it returns the case-B coefficients and zeros.
	!
	! WHICH SPECIES ABSORBS A PHOTON. Every channel emits photons of one
	! representative energy E_c; the cell keeps the fraction 1 - exp(-tau_c)
	! of them, tau_c = dl sum_s n_s sigma_s(E_c), and shares it among the
	! absorbers s in the ratio of n_s sigma_s(E_c). The rate at which ONE
	! absorber s takes them is then P_c sigma_s dl (1 - exp(-tau_c))/tau_c
	! (P_c the emission per unit volume), finite at tau_c -> 0 (the
	! optically thin rate of a single absorber) and -> P_c sigma_s/kappa_c
	! in a thick cell. Absorbers: H I, He I (ground), He II, H2 and the
	! photo-ionizable metal ions; each takes only what lies above its own
	! threshold. In a molecular gas H2 is not a small competitor
	! (sigma_H2/sigma_HI is 1.2 at 16 eV and 3.5-3.7 from 20 to 25 eV), and
	! the metals take a measured share of 3e-3 to 2e-2 on the WASP-121 b
	! and hot-Uranus cases despite abundances of 1e-4 to 1e-3.
	!
	! THE GROUND CAPTURES are the one channel whose photon can be taken by
	! the species that emitted it: the share re-absorbed by the recombined
	! species is a capture and a re-ionization that cancel, and is left out
	! of both the net coefficient (ground_capture_escape_weights) and the
	! rates below; the part taken by another absorber ionizes it; the part
	! that leaves the cell is a net recombination too.
	!
	! CHANNELS.
	!  H II -> H I ("H_rec_escape"): the ground capture, alpha_1(H I)
	!   (Cool_coeff: Milne relation), at 13.598 eV. The case-B captures
	!   into n >= 2 emit below 13.6 eV and ionize nothing here.
	!  He II -> He I ("He_rec_coupling"): the ground capture at 24.587 eV,
	!   and the case-B captures through their exits: the excited-singlet
	!   captures (Cool_coeff: alpha_rec_HeII_excited_singlets), 2/3 to 2^1P
	!   -> 584 A and 1/3 to 2^1S -> two photons (f_2q_HeI of them above
	!   the H I edge); and the captures into the triplets, all ending in
	!   2^3S, which leave it by the 19.8 eV 2^3S -> 1^1S line (A31), by
	!   conversion to 2^1S (q31a) or to 2^1P and the higher singlets (q31b),
	!   or by de-excitation to 1^1S (q31g, no photon). With the metastable
	!   tracked the last three act on the n(2^3S) of the network; without
	!   it (atomic mode) each triplet capture is shared among them in the
	!   ratio A31 : n_e q31a : n_e q31b : n_e q31g (the 2^3S photoionization
	!   and Penning ionization, which only the network carries, are not
	!   exits there). This is the construction of Draine (2011, sect.
	!   14.3.2: "Approximately 25% of these will be to states with total
	!   spin S = 0", "approximately 1/3 end up in 1s2s 1S0, and approximately
	!   2/3 in 1s2p 1P1o"; sect. 15.5: "z ~ 0.96 at low densities ... to
	!   z ~ 0.67 at high densities", z the H-ionizing photons of a case-B
	!   capture, with his critical density eqs. 14.19-14.20 the ratio
	!   A31/(q31a + q31b)) and of Osterbrock & Ferland (2006, sect. 2.4,
	!   p. 30, "p ~ 3/4 + 1/4 [2/3 + 1/3 (0.56)] = 0.96" and 0.66 at high
	!   density), with the triplet share of Hummer & Storey (1998;
	!   Cool_coeff: case_b_triplet_share_HeI) in place of 3/4 and the
	!   network's own Bray et al. (2000) rates in place of the fixed
	!   high-density split: z = 0.966 (n_e -> 0) and 0.660 (n_e = 1e8
	!   cm^-3) at 1e4 K (MEASURED, physics probe
	!   recombination_radiation_on_the_spot), in a cell without He I; He I
	!   scatters the 584 A photons and turns a share of them into 2^1S
	!   two-photon decays (VALIDITY (iii)), which lowers z.
	!  He III -> He II ("He_rec_coupling"): the ground capture,
	!   alpha_1(He II), at 54.418 eV (He II re-absorbs it on the spot; H I,
	!   He I, H2 and metals take the rest); the direct capture into n = 2
	!   (Cool_coeff: alpha_n2_hydrogenic_seaton), whose continuum starts at
	!   I_2 = 13.6047 eV and ionizes H I and metals only; and the case-B
	!   exits, every case-B capture ending in 2s (the fraction
	!   case_b_2s_fraction_hydrogenic) or 2p. 2p decays by He II
	!   Ly-alpha, 40.8135 eV; 2s decays by two photons (A_2q_HeII) unless an
	!   ion collision moves it to 2p first (Pengelly & Seaton 1964, the
	!   H+ and He2+ impact of their eq. 48, Cool_coeff:
	!   l_mixing_2s2p_pengelly_seaton), and the two-photon spectrum of
	!   Nussbaumer & Schmutz (1984) is carried in three bands (constants at
	!   E_2q_HeII_lo). Osterbrock & Ferland (2006, sect. 2.5, p. 34 and
	!   Table 2.6) name the same three H-ionizing products of He III -> He II
	!   recombination: He II Ly-alpha from 2p, the 2s two-photon continuum
	!   ("on the average, 1.42 ionizing photons are emitted per decay"; the
	!   three bands here sum to 1.425), and the n = 2 continuum of the
	!   direct captures. Their Table 2.6 implies a case-B 2s share of 0.28
	!   (1e4 K) and 0.31 (2e4 K) for He II (DERIVED), against 0.264 and
	!   0.293 from the Pengelly (1964) table used here. They take Ly-alpha
	!   and the n = 2 continuum as absorbed on the spot by H0 and let most
	!   two-photon photons leave the He++ zone; here every channel competes
	!   by cross section within the cell (VALIDITY (i)).
	!
	! ENERGY. A photoelectron receives E_c - I_s, plus for a continuum
	! channel the mean kinetic energy of the captured electron (E1 = beta/
	! alpha from the capture relation of Cool_coeff, capture_energy_loss_
	! rate), which the recombination cooling charges to every capture it
	! counts (lambda_rec_*: the net coefficients). All of it is deposited as
	! heat: it is not passed through the secondary-ionization partition of
	! the stellar photoelectrons, which that partition would leave as heat
	! below its 30 eV threshold anyway, except for the He III ground capture
	! (40.8 eV on H I) and the He II Ly-alpha on H I (27.2 eV).
	!
	! H2. An absorption by H2 ends in one of the four final states of
	! h2_photo_channels, H2+ + e- (M), H + H+ + e- (S), H+ + H+ + 2e- (D)
	! or H + H (N), in the ratio of the channel cross sections the stellar
	! field uses at the same energy (h2_photoabsorption_cross_sections,
	! with the run's neutral-window and double-ionization switches); the
	! total, dP_H2, is the absorption itself and is not changed by the
	! split. At the channel energies, M is open from 15.4 eV, S from
	! 18.08 eV (He I 2^3S line and 584 A, He I ground capture, the two
	! upper He II two-photon bands, He II Ly-alpha and ground capture), N
	! over the 32-41.5 eV window (He II Ly-alpha at 40.8 eV), D above
	! 51.4 eV (the He II ground capture at 54.4 eV). The heat of one
	! absorption is that of h2_channel_energy_recipients at E_c: the
	! photoelectron(s) and the fragment kinetic energy, not the formation
	! and ionization energy the products store (accounted by the reactions
	! that later consume them, as for the stellar field) nor the n = 2
	! excitation the N fragments radiate; the capture kinetic energy of a
	! continuum channel is added in every final state
	! (absorb_recombination_channel).
	!
	! VALIDITY. (i) LOCAL absorption: a photon that leaves the cell is
	! dropped rather than followed, so its absorption in some outer cell is
	! not counted; the escape weight is therefore a property of the cell
	! width as well as of the gas, and a finer grid moves it toward case A
	! wherever a cell is not thick. The wind thins outward, so most of what
	! leaves a thin cell escapes the domain; in a thick layer f is 1 and
	! nothing leaves. (ii) tau_c uses the cell width dr_j, i.e. a radially
	! escaping photon; the mean chord of an isotropically emitted photon in
	! a plane slab is longer by a factor of order 2. (iii) The He I 584 A
	! line is scattered by the cell's He I, and every scattering can end
	! it in the 2^1S two-photon continuum (absorb_recombination_channel,
	! after Wood, Mathis & Ercolano 2004, sect. 4.1). Their mean free path
	! averages the Doppler cross section over the emitting atom's
	! Maxwellian and does not follow the frequency redistribution; the
	! complete-redistribution mean of the same competition in one flight,
	! the integral over x of exp(-x^2)/pi^(1/2) beta/(beta + exp(-x^2))
	! (beta the continuum-to-line-centre opacity ratio), is 0.86 to 1.48
	! times theirs for beta = 1e-2 to 1e-6 (DERIVED here by quadrature, not
	! in the paper), and the natural damping wings (Voigt a = 1.3e-3 at
	! 1e4 K) are in neither. The escape from the cell is that of one
	! flight through the depth of the mean cross section, exp(-tau): for a
	! cell thick in the line (line-centre depth tau_0 >~ 10) a resonance
	! photon escapes by frequency diffusion into the wings with a
	! probability of order 1/(tau_0 (pi ln tau_0)^(1/2)) per scattering in
	! a static slab (the usual order-of-magnitude form, not from a source
	! read for this code; larger in a velocity gradient), which exp(-tau)
	! underestimates; that matters where the escape competes with eps_line
	! and the continuum, tau_0 of 1e1 to 1e3. Collisional transfer
	! 2^1P -> 2^1S by electrons, which competes with the 2.06 um branch once n_e
	! q(2^1P,2^1S) approaches A(2^1P,2^1S) = 1.97e6 s^-1, is not included;
	! its rate coefficient was not read. He II Ly-alpha is still a
	! continuum photon of its energy: its scattering by He II, which
	! lengthens its path and so raises the share the cell keeps, is not
	! followed. (iv) The He 2^3S metastable is not an absorber (its density
	! is 1e-6 to 1e-3 of He I), and the He I absorber density is the summed
	! He I column the caller passes. (v) The H2 final-state shares, and the
	! heat of each, are those at the representative energy E_c of the
	! channel, not averaged over the capture continuum above an edge; the
	! H2+ vibrational excitation the M channel leaves is inside the
	! measured cross section and is not followed, as in the stellar field.
	! The photoelectrons are not passed through the secondary-ionization
	! partition (ENERGY above); of the H2 events only those of the He II
	! ground capture (54.4 eV) carry an electron above its 30 eV threshold
	! (39.0 eV in M, 36.3 eV in S, plus the capture kinetic energy).
	subroutine recombination_radiation_absorbed(T_K, nhi, nhii, nh2,      &
	                           nhei, nheii, nheiii, nheiTR, ne, nm,       &
	                           A31, q31a, q31b,                           &
	                           rchiiB_new, rcheiiB_new, rcheiiiB_new,     &
	                           dP_HI, dP_HeI, dP_H2, dP_m, dheat, dP_m2, &
	                           dP_H2_di, dP_H2_dd, dP_H2_nd)
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi, nhii, nh2,    &
	                                              nhei, nheii, nheiii,   &
	                                              nheiTR, ne, q31a, q31b
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in)  :: nm
	real*8,                       intent(in)  :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rchiiB_new, rcheiiB_new, &
	                                              rcheiiiB_new, dP_HI,    &
	                                              dP_HeI, dP_H2, dheat
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) :: dP_m
	! The part of dP_m that ejects two or more electrons and that the
	! balance carries as a two-stage jump (absorb_recombination_channel).
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(out), optional :: dP_m2
	! The parts of dP_H2 that end in H + H+ + e- (S), H+ + H+ + 2e- (D)
	! and H + H (N): SUBSETS of the total H2 absorption rate dP_H2, the
	! same convention as P_H2_di, P_H2_dd and P_H2_nd of the stellar field,
	! so that dP_H2 minus the three is the H2+ (M) share.
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: dP_H2_di,     &
	                                                       dP_H2_dd, dP_H2_nd

	real*8, dimension(1-Ng:N+Ng) :: y_HI, y_gnd, y_HeII
	real*8, dimension(1-Ng:N+Ng) :: aB2, aB3, a1H, a1He, a1He2, a2He2
	real*8, dimension(1-Ng:N+Ng) :: E1_H, E1_He, E1_He2, E2_He2
	real*8  :: dl, t3, qa_c, qb_c, D_exit, f2s, mix, P2q, P_c
	real*8  :: dPm_cell(n_mion), dPm2_cell(n_mion), nm_cell(n_mion)
	real*8  :: kap_584, P_584_conv, h2_chan(n_h2_channels)
	logical :: h_on, he_on
	integer :: j

	h_on  = use_h_rec_escape
	he_on = use_he_rec_coupling .and. thereis_He

	call on_the_spot_cross_sections

	call ground_capture_escape_weights(nhi, nh2, nhei, nheii, nm,        &
	                                   y_HI, y_gnd, y_HeII)

	! The coefficients of the balance: case B plus the escaping ground
	! captures.
	rchiiB_new   = alpha_rec_HII_net(T_K, y_HI)
	rcheiiB_new  = alpha_rec_HeII_into_singlets(T_K, y_gnd)
	rcheiiiB_new = alpha_rec_HeIII_net(T_K, y_HeII)

	dP_HI  = 0.0d0
	dP_HeI = 0.0d0
	dP_H2  = 0.0d0
	dP_m   = 0.0d0
	if (present(dP_m2)) dP_m2 = 0.0d0
	if (present(dP_H2_di)) dP_H2_di = 0.0d0
	if (present(dP_H2_dd)) dP_H2_dd = 0.0d0
	if (present(dP_H2_nd)) dP_H2_nd = 0.0d0
	dheat  = 0.0d0
	if (.not. (h_on .or. he_on)) return

	! Capture coefficients and the mean kinetic energy [eV] of each
	! continuum channel's captured electrons, beta/alpha = k T (3/2 +
	! dln alpha/dln T) (Cool_coeff: capture_energy_loss_rate), which its
	! photon carries above the edge.
	a1H   = alpha_1_HI(T_K)
	a1He  = alpha_1_HeI(T_K)
	a1He2 = alpha_1_HeII(T_K)
	a2He2 = alpha_n2_hydrogenic_seaton(T_K, 2.0d0)
	aB2   = alpha_rec_HeII_B(T_K)
	aB3   = alpha_rec_HeIII_B(T_K)
	E1_H   = capture_kinetic_energy_eV(T_K, 1)
	E1_He  = capture_kinetic_energy_eV(T_K, 2)
	E1_He2 = capture_kinetic_energy_eV(T_K, 3)
	E2_He2 = capture_kinetic_energy_eV(T_K, 4)

	! Every quantity is cell-local, so the cells are shared among threads
	! and the result does not depend on their number.
	!$omp parallel do default(shared) schedule(static)                    &
	!$omp    private(j, dl, P_c, t3, qa_c, qb_c, D_exit, f2s, mix, P2q,   &
	!$omp            dPm_cell, dPm2_cell, nm_cell, kap_584, P_584_conv,   &
	!$omp            h2_chan)
	do j = 1-Ng,N+Ng
		dl = dr_j(j)*R0*1.0d-18
		nm_cell   = nm(j,:)
		dPm_cell  = 0.0d0
		dPm2_cell = 0.0d0
		h2_chan   = 0.0d0

		! ---- H II -> H I ----
		if (h_on) then
			P_c = a1H(j)*nhii(j)*ne(j)
			call absorb_recombination_channel(ic_gnd_HI, P_c, E1_H(j), 1, &
			        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,       &
			        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
		endif

		if (he_on) then
			! ---- He II -> He I ----
			! The 584 A line opacity of the cell's He I, in the units of
			! the continuum opacity of absorb_recombination_channel.
			kap_584    = nhei(j)*He_I_584_mean_cross_section(T_K(j))
			P_584_conv = 0.0d0
			P_c = a1He(j)*nheii(j)*ne(j)
			call absorb_recombination_channel(ic_gnd_HeI, P_c, E1_He(j), 2, &
			        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,       &
			        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
			if (.not. thereis_HeITR) then
				! atomic mode: the case-B captures through their exits (see
				! CHANNELS above); t3 of them end in 2^3S, which radiates
				! the 19.8 eV line (A), is converted to 2^1S (n_e q31a) or
				! 2^1P and the higher singlets (n_e q31b), or is
				! de-excited to 1^1S (n_e q31g), in the ratio of those
				! rates.
				P_c   = aB2(j)*nheii(j)*ne(j)
				t3    = case_b_triplet_share_HeI(T_K(j))
				qa_c  = excitation_rate_HeI_23S_21S(T_K(j))
				qb_c  = excitation_rate_HeI_23S_21P(T_K(j))                 &
				      + excitation_rate_HeI_23S_singlets_n3(T_K(j))
				D_exit = A_HeI_23S_11S + ne(j)*(qa_c + qb_c                 &
				       + deexcitation_rate_HeI_23S_11S(T_K(j)))
				call absorb_recombination_channel(ic_19_HeI,  &
				        t3*A_HeI_23S_11S/D_exit*P_c, 0.0d0, 0,            &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_584_HeI, &
				        ((1.0d0 - t3)*2.0d0/3.0d0                         &
				         + t3*ne(j)*qb_c/D_exit)*P_c, 0.0d0, 0,           &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j), &
				        kap_584, eps_HeI_21P_21S, P_584_conv)
				call absorb_recombination_channel(ic_2q_HeI,  &
				        f_2q_HeI*((1.0d0 - t3)/3.0d0                      &
				         + t3*ne(j)*qa_c/D_exit)*P_c, 0.0d0, 0,           &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
			else
				! metastable mode: the explicit exits
				P_c = alpha_rec_HeII_excited_singlets(T_K(j))*nheii(j)*ne(j)
				call absorb_recombination_channel(ic_584_HeI, 2.0d0/3.0d0*P_c, &
				        0.0d0, 0, nhi(j), nhei(j), nheii(j), nh2(j),      &
				        nm_cell, dl, dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan,       &
				        dPm_cell, dPm2_cell, dheat(j),                    &
				        kap_584, eps_HeI_21P_21S, P_584_conv)
				call absorb_recombination_channel(ic_2q_HeI, f_2q_HeI/3.0d0*P_c, &
				        0.0d0, 0, nhi(j), nhei(j), nheii(j), nh2(j),      &
				        nm_cell, dl, dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan,       &
				        dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_19_HeI, A31*nheiTR(j), &
				        0.0d0, 0, nhi(j), nhei(j), nheii(j), nh2(j),      &
				        nm_cell, dl, dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan,       &
				        dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_2q_HeI,  &
				        f_2q_HeI*q31a(j)*ne(j)*nheiTR(j), 0.0d0, 0,       &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_584_HeI, &
				        q31b(j)*ne(j)*nheiTR(j), 0.0d0, 0,                &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j), &
				        kap_584, eps_HeI_21P_21S, P_584_conv)
			endif
			! The 584 A photons the He I scattering converted to 2^1S
			! (each leaving a 2.06 um photon of 0.60 eV to escape) decay
			! through the 2^1S two-photon continuum, f_2q_HeI H-ionizing
			! photons each.
			call absorb_recombination_channel(ic_2q_HeI,                  &
			        f_2q_HeI*P_584_conv, 0.0d0, 0,                        &
			        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,       &
			        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))

			! ---- He III -> He II ----
			if (nheiii(j) > 0.0d0) then
				P_c = a1He2(j)*nheiii(j)*ne(j)
				call absorb_recombination_channel(ic_gnd_HeII, P_c, E1_He2(j), 3, &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_n2_HeII, &
				        a2He2(j)*nheiii(j)*ne(j), E2_He2(j), 0,           &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				f2s = case_b_2s_fraction_hydrogenic(T_K(j), 2.0d0)
				mix = nhii(j)*l_mixing_2s2p_pengelly_seaton(T_K(j),       &
				          ne(j), 2.0d0, 1.0d0, mu_HeII_p, dE_2s2p12_HeII, &
				          dE_2s2p32_HeII, A_2q_HeII)                      &
				    + nheiii(j)*l_mixing_2s2p_pengelly_seaton(T_K(j),     &
				          ne(j), 2.0d0, 2.0d0, mu_HeII_He2p,              &
				          dE_2s2p12_HeII, dE_2s2p32_HeII, A_2q_HeII)
				P2q = A_2q_HeII/(A_2q_HeII + mix)
				P_c = aB3(j)*nheiii(j)*ne(j)
				call absorb_recombination_channel(ic_lya_HeII, &
				        P_c*(1.0d0 - f2s*P2q), 0.0d0, 0,                  &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_2q_HeII_lo, &
				        f_2q_HeII_lo*f2s*P2q*P_c, 0.0d0, 0,               &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_2q_HeII_mid, &
				        f_2q_HeII_mid*f2s*P2q*P_c, 0.0d0, 0,              &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
				call absorb_recombination_channel(ic_2q_HeII_hi, &
				        f_2q_HeII_hi*f2s*P2q*P_c, 0.0d0, 0,               &
				        nhi(j), nhei(j), nheii(j), nh2(j), nm_cell, dl,   &
				        dP_HI(j), dP_HeI(j), dP_H2(j), h2_chan, dPm_cell, dPm2_cell, dheat(j))
			endif
		endif
		dP_m(j,:) = dPm_cell
		if (present(dP_m2)) dP_m2(j,:) = dPm2_cell
		if (present(dP_H2_di)) dP_H2_di(j) = h2_chan(ICH_S)
		if (present(dP_H2_dd)) dP_H2_dd(j) = h2_chan(ICH_D)
		if (present(dP_H2_nd)) dP_H2_nd(j) = h2_chan(ICH_N)
	enddo
	!$omp end parallel do

	end subroutine recombination_radiation_absorbed

	! Mean kinetic energy [eV] of the electrons captured into the ground
	! level of H I (1), He I (2), He II (3), and into n = 2 of He II (4):
	! beta/alpha = k T (3/2 + dln alpha/dln T) (Cool_coeff:
	! capture_energy_loss_rate, dlnT_capture).
	elemental double precision function capture_kinetic_energy_eV(T, ich)
	real*8,  intent(in) :: T
	integer, intent(in) :: ich
	real*8 :: a_lo, a_0, a_hi, T_lo, T_hi
	T_lo = T*exp(-dlnT_capture)
	T_hi = T*exp( dlnT_capture)
	select case (ich)
		case (1)
			a_lo = alpha_1_HI(T_lo);  a_0 = alpha_1_HI(T);  a_hi = alpha_1_HI(T_hi)
		case (2)
			a_lo = alpha_1_HeI(T_lo); a_0 = alpha_1_HeI(T); a_hi = alpha_1_HeI(T_hi)
		case (3)
			a_lo = alpha_1_HeII(T_lo); a_0 = alpha_1_HeII(T)
			a_hi = alpha_1_HeII(T_hi)
		case default
			a_lo = alpha_n2_hydrogenic_seaton(T_lo, 2.0d0)
			a_0  = alpha_n2_hydrogenic_seaton(T,    2.0d0)
			a_hi = alpha_n2_hydrogenic_seaton(T_hi, 2.0d0)
	end select
	capture_kinetic_energy_eV = capture_energy_loss_rate(T, a_lo, a_0, a_hi) &
	                            /a_0*erg2eV
	end function capture_kinetic_energy_eV


	! End of module
	end module utils_ion_eq

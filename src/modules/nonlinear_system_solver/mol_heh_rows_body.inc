!  THE EIGHT BALANCE ROWS OF THE MOLECULAR NETWORK, ONE SOURCE TEXT FOR TWO
!  ARITHMETICS: included by mol_heh_rows_real (real*8 densities: the local
!  sweep, System_HeH_mol and System_HeH_mol_metals, and the carrier
!  residual) and by mol_heh_rows_dual (dual densities: the rows and their
!  exact first derivatives, for the carrier Jacobian) of System_HeH_mol.
!  The named production and loss halves of rows (1) and (4) are defined
!  here; the optional diagnostics that read them (production/loss split,
!  H2 channels, gross rates, channel record) stay in mol_heh_rows_real.

	! Rate coefficients hoisted once per cell into module state by
	! set_mol_coeffs (they depend only on T and n_tot); read here.
	k5  = mk5
	k6  = mk6
	k7  = mk7
	k8  = mk8
	k9  = mk9
	k10 = mk10
	k11 = mk11
	k12 = mk12
	k13 = mk13
	k14 = mk14
	! THE THIRD BODY OF R15 AND OF ITS REVERSE R12, at this composition.
	! mk15 is k1(H2) [cm^6 s^-1] and n_third the H2-equivalent collider
	! density, so k15 is the two-body-equivalent rate [cm^3 s^-1] the row
	! multiplies n(H)^2 by, and k12 n_third the one it multiplies n(H2) by.
	! Both directions carry the SAME density, so the pair is an exact
	! detailed balance collider by collider (mol_rates).  The density
	! depends on n(H2) and on n(H), which are unknowns of this system, so
	! the numerical Jacobian of the solve picks up the derivative of the
	! collider sum as well as that of the densities the rate multiplies.
	n_third = h2_third_body_density(ntot, n_h2, n_hi, n_heiSI)
	k15 = mk15*n_third
	k16 = mk16
	k17 = mk17
	k18 = mk18
	k19 = mk19
	k_h2p_he = mk_h2p_he
	k23 = mk23
	! Total (Penning + associative) He(2^3S)+H2 ionization; 0 without triplet.
	k_ion_H2 = mk_ion_H2

	! He+ + CO -> C+ + O + He, as a rate per He+ ion. The CO density is the
	! cell's own reservoir ieq_cell%n_co and is a BACKGROUND density for
	! this system, not one of its unknowns: CO is a transported carrier and
	! the one-sided model gives it no formation channel, so it has no local
	! equilibrium and cannot be a row here. It is frozen for the duration of
	! the solve exactly as the ion stages are frozen for the duration of the
	! carrier step that transports CO.
	k_co_hep = 0.0d0
	if (thereis_oxychem) k_co_hep = rk_D1_Hep_CO()*ieq_cell%n_co

	! (1) H+ balance
	! g_h2_di*n_h2: the dissociative branch of the H2 photoionization,
	! H2 + hv -> H + H+ + e-. Its H atom is supplied by the H-nucleus
	! closure, and the H2 it consumes is already in row (4) through g_h2,
	! of which g_h2_di is a part.
	! 2*g_h2_dd*n_h2: the double branch, H2 + hv -> H+ + H+ + 2e-. The
	! factor 2 is the stoichiometry, TWO protons out of one event -- while
	! the H2 destroyed is still one per event, which is why row (4) is not
	! doubled with it. It leaves no H atom, so the closure again needs no
	! term. Zero unless a double-ionization model is selected.
	! f_penning_HeI23S*Q31*n_heiTR*n_hi: Penning ionization
	! He(2^3S)+H0 -> He(1^1S)+H+ + e-. Q31 is the TOTAL He(2^3S)+H ionization
	! rate; the associative 10% makes HeH+ instead of a proton and is carried
	! by row (7). (n_heiTR is 0 when thereis_HeITR is false, so the term is
	! unconditional.)
	prod_hp     = (g_hi + b_hi*n_e)*n_hi                              &
	            + (g_h2_di + 2.0d0*g_h2_dd)*n_h2                      &
	            + k9*n_h2p*n_hi + k17*n_heii*n_h2                     &
	            + f_penning_HeI23S*Q31*n_heiTR*n_hi
	loss_hp_rec = a_hii*n_e*n_hii
	loss_hp_h2  = (k10 + k13)*n_hii*n_h2
	fvec(1) = prod_hp - loss_hp_rec - loss_hp_h2

	! (2) He+ balance (atomic part consistent with System_HeH_TR rows
	!     2+3 combined; + molecular sinks R17/R23). b_heiTR*n_e*n_heiTR is
	!     electron-impact ionization of the He 2^3S metastable into He+.
	! k_co_hep*n_heii: He+ + CO -> C+ + O + He, UMIST RATE22 entry 4068,
	! measured, temperature independent at the Langevin value 1.6e-9
	! cm^3 s^-1 (oxygen_rates::rk_D1_Hep_CO). It neutralizes the helium ion,
	! and the free neutral He that closes the helium budget receives it, so
	! no other helium row moves. The carbon it ionizes is not a row of this
	! system: the metal block re-solves C I / C II / C III from its own
	! balance, so this reaction moves the helium charge here and the carbon
	! charge is set there. The same rate coefficient and the same two
	! densities carry the reaction's +2.2117 eV into heating channel 18
	! (util_ion_eq::heating_of_composition), so the energy released and the
	! ions consumed are one event counted once.
	fvec(2) = (g_hei + b_hei*n_e)*n_heiSI + g_heiTR*n_heiTR           &
	        + b_heiTR*n_e*n_heiTR                                     &
	        + a_heiii*n_e*n_heiii                                     &
	        - (a_heii + a_heiTR)*n_e*n_heii                           &
	        - (g_heii + b_heii*n_e)*n_heii                            &
	        - (k17 + k23)*n_heii*n_h2                                 &
	        - k_co_hep*n_heii

	! (3) He++ balance (verbatim atomic form + collisional ionization)
	fvec(3) = (g_heii + b_heii*n_e)*n_heii - a_heiii*n_e*n_heiii

	! (4) H2 balance
	! k_ion_H2*n_heiTR: loss of H2 to He(2^3S)+H2 ionization. The total is
	! taken, Penning and associative alike: both channels consume the H2.
	! (n_heiTR is 0 when thereis_HeITR is false, so the term is unconditional).
	! g_lw: Lyman-Werner photodissociation H2 + hv -> H + H, already
	! self-shielded (lyman_werner.f90); 0 when the run supplies no band flux.
	! g_h2 is the TOTAL H2 destruction rate and already contains all four
	! product channels, the neutral dissociation among them, so this row is
	! the same whichever of them are resolved.
	prod_h2 = k6*n_e*n_h3p + k9*n_h2p*n_hi + k11*n_h3p*n_hi           &
	        + k15*n_hi*n_hi
	loss_h2 = ( g_h2 + g_lw + (k10 + k13)*n_hii + k12*n_third + k14*n_e  &
	          + k8*n_h2p + (k17 + k23)*n_heii + k18*n_hehp            &
	          + k_ion_H2*n_heiTR )*n_h2
	fvec(4) = prod_h2 - loss_h2

	! (5) H2+ balance
	! (g_h2 - g_h2_di - g_h2_dd - g_h2_nd): only the branch that leaves the
	! molecule bound makes an H2+ ion. The two ionizing branches removed
	! here are proton sources in row (1); the neutral one makes no ion at
	! all. Both of the latter two are zero unless their options are on.
	! f_penning_HeI23S*k_ion_H2*n_heiTR*n_h2: H2+ produced by the Penning
	! branch of He(2^3S)+H2 ionization. The associative branch of the same
	! collision gives H + HeH+ + e- and appears in row (7) instead.
	! k_h2p_he*n_heiSI: H2+ + He -> HeH+ + H, the HeH+ source the network
	! carries since Koskinen R20 was retired (mol_rates: the paper R20
	! cites measures no HeH+ channel and bounds it 42 times below the value
	! Table 1 gives). The collision partner is GROUND-STATE helium; the
	! metastable's own two channels are the associative branch below.
	fvec(5) = (g_h2 - g_h2_di - g_h2_dd - g_h2_nd)*n_h2               &
	        + k10*n_hii*n_h2 + k11*n_h3p*n_hi                         &
	        + k19*n_hehp*n_hi + k23*n_heii*n_h2                       &
	        + f_penning_HeI23S*k_ion_H2*n_heiTR*n_h2                  &
	        - (k5*n_e + k8*n_h2 + k9*n_hi                             &
	           + k_h2p_he*n_heiSI)*n_h2p

	! (6) H3+ balance
	fvec(6) = k8*n_h2p*n_h2 + k13*n_hii*n_h2 + k18*n_hehp*n_h2        &
	        - ((k6 + k7)*n_e + k11*n_hi)*n_h3p

	! (7) HeH+ balance
	! k_h2p_he*n_heiSI*n_h2p: H2+ + He -> HeH+ + H (Black 1978, from the
	! Chupka et al. 1969 cross sections; see mol_rates). It is the network's
	! only HeH+ source that does not need the metastable, and it replaces
	! Koskinen R20, which was retired because the paper Table 1 cites
	! for it measures no HeH+ channel and bounds it at <= 1.0e-14 against
	! the 4.2e-13 it was given. A helium NUCLEUS moves from free neutral He
	! into HeH+ here and the helium closure of this system, free neutral
	! He = n_He(1 - x2 - x3) - n_HeH+, accounts for it with no second term.
	! The associative branch of He(2^3S) ionization, (1 - f_penning_HeI23S) of
	! the total in both collision partners, ends in HeH+ rather than in H+ or
	! H2+ (Garcia Munoz 2025, network rows 199 and 203):
	!   He(2^3S) + H  -> HeH+ + e-
	!   He(2^3S) + H2 -> H + HeH+ + e-.
	! It is resolved here because this system carries HeH+ explicitly and the
	! onward chemistry that breaks the HeH+ -> He + H recombination cycle
	! (k18: HeH+ + H2 -> H3+ + He) is present. The atomic systems, which have
	! no HeH+ row, drop the branch instead; see ion_residual_core.
	fvec(7) = k_h2p_he*n_heiSI*n_h2p                                  &
	        + (1.0d0 - f_penning_HeI23S)*n_heiTR                      &
	          *(Q31*n_hi + k_ion_H2*n_h2)                             &
	        - (k16*n_e + k18*n_h2 + k19*n_hi)*n_hehp

	! (8) He 2^3S balance (VERBATIM System_HeH_TR row 4)
	if (thereis_HeITR) then
		call tr_triplet_row(fvec(8), n_hi, n_heiSI, n_heiTR, n_heii, n_e,   &
		                    g_heiTR, a_heiTR, q13, q31g, q31a, q31b, Q31,   &
		                    A31, b_heiTR)
		! He(2^3S)+H2 ionization triplet loss, total of the Penning and the
		! associative branch: both quench the metastable, and neither makes
		! He+, so there is no He+ row term. Garcia Munoz (2025) Table A.5.
		fvec(8) = fvec(8) - k_ion_H2*n_heiTR*n_h2
	endif

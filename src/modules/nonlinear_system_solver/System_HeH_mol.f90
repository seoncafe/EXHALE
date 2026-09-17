	module System_HeH_mol
	! molecular ionization-equilibrium system: H/He (+ optional
	! He 2^3S triplet) extended with H2, H2+, H3+, HeH+ (docs/
	! lower_atmosphere_coupling.*).
	!
	! Unknowns (fractions):
	!   x(1) = n_HII / n_H(nuclei)
	!   x(2) = n_HeII / n_He          x(3) = n_HeIII / n_He
	!   x(4) = 2 n_H2  / n_H          (H nuclei bound in H2)
	!   x(5) = 2 n_H2+ / n_H          x(6) = 3 n_H3+ / n_H
	!   x(7) = n_HeH+ / n_H
	!   x(8) = n_HeITR / n_He         (only when thereis_HeITR)
	! Neutral atomic H and neutral He close the element budgets (the He
	! budget includes the He nucleus carried by HeH+).
	!
	! With the oxygen chemistry on (thereis_oxychem) two more unknowns sit
	! above those, at oxygen_row_base() of System_HeH_mol_metals:
	!   x(iox)   = n_OH  / n_O(free family)
	!   x(iox+1) = n_H2O / n_O(free family)
	! and free atomic oxygen closes the oxygen budget the way atomic H
	! closes the hydrogen one. The carriers also hold H nuclei (OH one, H2O
	! two), so the atomic-H closure gains them; that is done at the call
	! site, where n_H is known. The rows themselves are oxygen_carrier_rows
	! below. Oxygen is a metal element here, so the option is only ever
	! solved by System_HeH_mol_metals; the metals-free system below cannot
	! reach it (input_read refuses the key without oxygen).
	!
	! Rows are steady-state production-loss balances.  The atomic rows
	! use EXHALE's own rate arrays (P_HI/rchiiB/Voronov...) so the
	! molecular-free limit reproduces the atomic systems' solution; the
	! He 2^3S row follows System_HeH_TR plus the He(2^3S)+H2 ionization loss
	! (Garcia Munoz 2025 Table A.5), the dominant metastable sink toward an
	! H2-rich base.  Molecular channels are
	! the Koskinen et al. (2022) Table-1 network via mol_rates, plus the
	! Lyman-Werner photodissociation the Table-1 network omits
	! (lyman_werner.f90; opt-in, see below).  H <-> He charge exchange IS
	! carried, but not from mol_rates: it enters through he_h_cx_fvec, the
	! Huang et al. (2023) Table-4 B1/B2 pair every system with He shares, so
	! the molecular-free limit reproduces the atomic systems exactly.  The
	! Table-1 R21/R22 are therefore transcribed in mol_rates but unused.
	!
	! Where the cell state comes from: ieq_cell (ion_cell_state.f90), one
	! named field per quantity, P_H2 and P_H2_di among them.  Every
	! photoionization rate among them is a CONSTANT of this solve although
	! the field that sets it is a function of these unknowns; the reason,
	! and where that dependence is closed instead, are stated once in
	! System_HeH.  The
	! params(40) argument ion_system_HeH_mol still carries is the MINPACK
	! callback signature and nothing else -- no element of it is read here
	! -- so there is no numbered layout to keep in step.
	!
	! Lyman-Werner photodissociation H2 + hv -> H + H enters row 4 next to
	! the H2 photoionization, as the self-shielded rate k_LW carried by the
	! cell state (lyman_werner.f90).  It is zero unless the run supplies a
	! Lyman-Werner band flux, so the molecular network without it is
	! unchanged.  Both products are neutral H, which the H-nucleus closure
	! (1 - x1 - x4 - x5 - x6 - x7) supplies automatically; no other row moves.

	! THE FOUR PRODUCT CHANNELS of the H2 photoabsorption
	! (h2_photo_channels.f90) are carried as three SUBSETS of the total rate
	! P_H2: P_H2_di for H + H+ + e- (18.08 eV), P_H2_dd for H+ + H+ + 2e-
	! (51.4 eV) and P_H2_nd for the neutral H + H (33-41 eV window), the
	! remainder being the channel that leaves H2+.  Each destroys one H2, so
	! the H2 destruction row (4) sees P_H2 unchanged; the H2+ production row
	! (5) sees P_H2 minus the three of them; and row (1) gains
	! P_H2_di + 2 P_H2_dd as a proton source, the factor 2 because the
	! double channel releases two protons in one event.  The H atoms of
	! every channel are supplied by the H-nucleus closure, as the
	! Lyman-Werner fragments are: each channel conserves H nuclei, so the
	! closure needs no term of its own.  Both options are ON by default;
	! switching both off zeroes P_H2_dd and P_H2_nd and reproduces the
	! two-channel arithmetic exactly.  The dissociative
	! channel therefore removes H2 without feeding the H3+ chain (R8), which
	! is what makes it worth resolving in a layer whose proton budget and
	! H3+ cooling are the quantities of interest.  Yan, Sadeghpour &
	! Dalgarno (1998) sec. 4 for the photon branching, Dalgarno, Yan & Liu
	! (1999) after their eq. (10) for the secondary-electron one.
	!
	! The eight balance rows live in mol_heh_rows below, which takes the free
	! electron density as an INPUT. System_HeH_mol_metals calls the same
	! routine with a metal-inclusive n_e, so the molecular network is written
	! once (the rate coefficients mk5..mk23 likewise stay in this module).

	use global_parameters, only: thereis_HeITR, thereis_oxychem
	use mol_rates
	! Oxygen-chemistry coefficients (the A2 option): the three rate
	! coefficients of the audited set, their thermodynamic reverses and the
	! photolysis channels. Read only when thereis_oxychem.
	use oxygen_rates, only: rk_O1_OH_H2_water, rk_O2_O_H2_hydroxyl,       &
	                        rk_O6_O1D_H2_hydroxyl,                        &
	                        rate_from_detailed_balance,                   &
	                        ith_H, ith_H2, ith_O, ith_OH, ith_H2O,        &
	                        n_fuv_band, qy_H2O_OH_H, qy_H2O_H2_O1D,       &
	                        qy_H2O_O_H_H, rk_D1_Hep_CO
	use ion_residual_core, only: tr_triplet_row,                        &
	                            impose_transported_ionization_fractions
	use ion_cell_state, only: ieq_cell
	use charge_exchange, only: he_h_cx_fvec
	use Cooling_Coefficients, only: ioniz_HeI23S_H2,      &
	                               f_penning_HeI23S
	! The control key of the collider-resolved third body, read from the one
	! accessor that owns it, so that EXHALE_REACTION_HEAT_RECIPIENTS=0 means
	! the same thing in the chemistry and in the energy ledger.
	use molecular_reaction_heat, only: reaction_heat_recipients_corrected

	implicit none

	! Molecular rate coefficients, invariant across a cell's Newton solve
	! (depend only on T and n_tot). Computed once per cell by set_mol_coeffs
	! and read by the residual, mirroring set_metal_coeffs / cx_set_cell.
	!
	! mk15 IS THE THREE-BODY COEFFICIENT k1(H2) ITSELF [cm^6 s^-1], not a
	! two-body-equivalent rate: the third body of R15 and of its reverse R12
	! is a COLLIDER SUM, k1(H2) n(H2) + k1(H) n(H) + k1(Ar) n(He) (Cohen &
	! Westberg 1983, p. 559), and n(H2) and n(H) are unknowns of this system.
	! The H2-equivalent third-body density is therefore formed inside the
	! rows, at the composition the row is evaluated at, and mk15 and mk12
	! carry only the temperature dependence they had.  The ledger that
	! charges the same two channels their enthalpy forms the same density
	! from the same routine (molecular_reaction_heat), so one R15 rate sets
	! the H2 abundance and its heat.
	real*8, save :: mk5,mk6,mk7,mk8,mk9,mk10,mk11,mk12,mk13,mk14,mk15
	real*8, save :: mk16,mk17,mk18,mk19,mk_h2p_he,mk23
	! He(2^3S)+H2 TOTAL ionization rate coefficient, Penning plus associative
	! (Garcia Munoz 2025 Table A.5); depends only on T. Zero-effect unless the
	! triplet is present.
	real*8, save :: mk_ion_H2
	! Whether the third body of R12 and R15 is the collider sum (the
	! published physics) or the total heavy-particle density (the control
	! key).  Set with the coefficients, read by the rows.
	logical, save :: mk_third_body_resolved

	! ---------------------------------------------------------------
	! Oxygen chemistry (the A2 option, docs/a2_oxygen_option_design.md).
	! Cell-invariant coefficients, hoisted once per cell by
	! set_oxygen_coeffs exactly as the molecular ones are, and read by
	! oxygen_carrier_rows. All zero, and never read, without the option.
	!
	!   ok1  O1   OH + H2  -> H2O + H        Baulch et al. (2005) p. 1029
	!   ok1r O1r  H2O + H  -> OH + H2        detailed balance of ok1
	!   ok2  O2   O + H2   -> OH + H         Baulch et al. (2005) p. 804
	!   ok2r O2r  OH + H   -> O + H2         detailed balance of ok2
	!   ok6  O6   O(1D)+H2 -> OH + H         Atkinson et al. (2004) I.A2.18
	!   oj3  O3   H2O + hv -> OH + H         summed over the FUV bands
	!   oj4  O4   H2O + hv -> H2 + O(1D)     summed over the FUV bands
	!   oj5  O5   H2O + hv -> O + H + H      summed over the FUV bands
	!   oj7  O7   OH  + hv -> O + H          summed over the FUV bands
	!
	! ok1r and ok2r are NOT transcribed: O1 and O2 run close to cancellation
	! on a hot base, so an independently transcribed reverse produces an
	! arbitrary net rather than a small error. Detailed balance against the
	! module's Shomate table makes the pair exact by construction and makes
	! the hot limit reduce to chemical equilibrium (decision D2; the measured
	! validation against a published reverse rate is
	! docs/a2_reaction_audit.md sec. 5).
	!
	! ok6 is used ONLY for the O(1D) number-density diagnostic. O(1D)'s only
	! sink in the audited set is O6 itself, so in steady state the O6 flux
	! equals the O(1D) production rate oj4*n_H2O whatever the H2 density is,
	! and the network never divides by n_H2. That is why O4 and O6 enter the
	! rows as one combined channel; see oxygen_carrier_rows.
	real*8, save :: ok1, ok1r, ok2, ok2r, ok6
	real*8, save :: oj3, oj4, oj5, oj7

	! Reciprocal turnover scale of each balance row.  A row of this system is
	! a production-loss balance in cm^-3 s^-1; its turnover scale is the rate
	! at which the species it balances can be produced or destroyed in this
	! cell.  Every row is multiplied by the reciprocal below before hybrd1
	! sees it, so the residual the solver works on is the RELATIVE imbalance
	! of each species balance rather than an absolute reaction rate.
	!
	! Why the system needs it.  The magnitude of a row is set by the element
	! that carries it: the helium rows (2, 3, 8) run as n_He times a rate and
	! the hydrogen and molecular rows (1, 4-7) as n_H times a rate, with a
	! second density inside every bilinear term.  At He/H = 10^3 the two
	! blocks of one residual vector therefore differ by about 10^6.  MINPACK's
	! hybrd1 scales the variables (its internal diag) but never the rows: the
	! dogleg step minimizes ||J dx + f||_2, in which a block 10^6 below the
	! other carries no weight, so the unknowns that block determines -- the
	! molecular fractions -- are left wherever the helium block puts them.
	! That is the route by which a converged root leaves the physical simplex
	! and has to be clamped onto the element budget.  A positive diagonal
	! scaling of the residual has exactly the same zeros, so this changes the
	! path to the root and not the root.
	!
	! What the scale is.  The sum of the magnitudes the terms of the row reach
	! when every species is set to the whole of its element: each rate
	! coefficient of the row times the reference densities n_H, n_He, n_e and
	! n_tot of the cell.  It is an upper bound on the row, strictly positive
	! whenever the row carries any reaction, and constant across the cell's
	! solve -- so the scaled residual is the same smooth function of x that
	! the finite-difference Jacobian of hybrd1 assumes.  Rows 9.. (the metal
	! block of System_HeH_mol_metals) are filled by that module; rows left
	! unset keep the 1 that set_mol_turnover_rates writes over the whole
	! array, which leaves them untouched.  The array carries no declaration
	! initializer on purpose: it is threadprivate, where an initializer
	! reaches the master thread only, and every entry is written per cell by
	! set_mol_turnover_rates before any residual reads it.
	!
	! set_mol_turnover_rates below MIRRORS mol_heh_rows term by term.  A
	! reaction added to a row must be added to its scale as well; a missing
	! term only makes the scale a weaker bound, it cannot make the root wrong.
	! >= largest N_eq: 7 molecular + 1 triplet + 2 per metal element (10).
	integer, parameter :: n_mol_rows_max = 40
	real*8, save :: mol_inv_turnover(n_mol_rows_max)
	!$omp threadprivate(mk5,mk6,mk7,mk8,mk9,mk10,mk11,mk12,mk13,mk14,mk15, &
	!$omp                mk16,mk17,mk18,mk19,mk_h2p_he,mk23,mk_ion_H2,     &
	!$omp                mk_third_body_resolved,                           &
	!$omp                ok1,ok1r,ok2,ok2r,ok6,oj3,oj4,oj5,oj7,            &
	!$omp                mol_inv_turnover)

	contains

	! Compute the molecular rate coefficients that are invariant across a
	! cell's Newton solve (they depend only on T and n_tot). Called once per
	! cell from ioniz_eq before hybrd1, like set_metal_coeffs / cx_set_cell.
	subroutine set_mol_coeffs(T, ntot)
	real*8, intent(in) :: T, ntot
	mk5  = rk_R5_H2p_dr(T)
	mk6  = rk_R6_H3p_dr_H2(T)
	mk7  = rk_R7_H3p_dr_3H(T)
	mk8  = rk_R8_H2p_H2()
	mk9  = rk_R9_H2p_H()
	mk10 = rk_R10_Hp_H2v4(T)
	mk11 = rk_R11_H3p_H(T)
	mk12 = rk_R12_H2_thdis(T)
	mk13 = rk_R13_Hp_H2_M(ntot)
	mk14 = rk_R14_H2_edis(T)
	! k1(H2) of Cohen & Westberg (1983); the collider density it multiplies
	! is a function of the unknowns and is formed in mol_heh_rows.
	mk15 = k3b_H_H_to_H2(T)
	mk16 = rk_R16_HeHp_dr(T)
	mk17 = rk_R17_Hep_H2_diss(T)
	mk18 = rk_R18_HeHp_H2()
	mk19 = rk_R19_HeHp_H()
	mk_h2p_he = rk_H2p_He_HeHp(T)
	mk23 = rk_R23_H2_Hep_cx()
	mk_ion_H2 = ioniz_HeI23S_H2(T)     ! He(2^3S)+H2 total ionization
	! WHETHER THE THIRD BODY OF THE R12/R15 PAIR IS COLLIDER-RESOLVED.  The
	! key EXHALE_REACTION_HEAT_RECIPIENTS=0 restores the total
	! heavy-particle density, and it has to mean the same thing in the
	! chemistry as in the energy ledger or the H2 abundance is set by one
	! rate and its heat charged at another.  The accessor that owns the key
	! caches it in state shared by the threads, and this routine runs inside
	! the sweep's parallel region, so the read is serialized; it is one
	! environment-variable read per cell against a Newton solve, and the
	! cost is not measurable against it.
	!$omp critical (exhale_reaction_heat_recipients_key)
	mk_third_body_resolved = reaction_heat_recipients_corrected()
	!$omp end critical (exhale_reaction_heat_recipients_key)
	end subroutine set_mol_coeffs

	! The H2-EQUIVALENT third-body density of the R12/R15 pair at the
	! composition the row is evaluated at [cm^-3], from the one routine that
	! owns the collider sum (mol_rates::h2_association_collider_density).
	! Both directions of the pair are multiplied by it, so the pair stays an
	! exact detailed balance collider by collider whatever the mixture is.
	! The helium collider is the GROUND-SINGLET neutral, which is what the
	! energy ledger passes as well; an ion is not the third body this
	! coefficient was measured for.
	double precision function h2_third_body_density(ntot, n_h2, n_hi,     &
	                                                n_hei) result(n_eff)
	real*8, intent(in) :: ntot, n_h2, n_hi, n_hei
	if (mk_third_body_resolved) then
		n_eff = h2_association_collider_density(ieq_cell%T_K, n_h2,   &
		                                        n_hi, n_hei)
	else
		n_eff = ntot
	endif
	end function h2_third_body_density

	! Turnover scale of every balance row of the molecular system, from the
	! cell's element densities and the rate coefficients just hoisted by
	! set_mol_coeffs.  Called once per cell from ioniz_eq, after
	! set_mol_coeffs and before the hybrd1 attempts; see the declaration of
	! mol_inv_turnover above for what the scale is and why the system needs
	! one.  n_e_ref is the cell's incoming electron density -- the electron
	! density the row terms actually see, not the neutrality bound
	! n_H + 2 n_He, which in the shielded molecular base overstates them by
	! many orders and would scale those rows into insignificance.
	!
	! photo_scale multiplies every RADIATION-DRIVEN rate of the bound
	! (P_HI, P_HeI, P_HeII, P_HeITR, P_H2 and the Lyman-Werner k_LW) and
	! nothing else, so that a solve run at a scaled radiation field is
	! judged against the turnover its own rows carry. The cell solve passes
	! 1, which is exact in IEEE arithmetic and leaves the term order
	! untouched; the radiation-field continuation of
	! constrained_chemical_equilibrium passes its lambda.
	subroutine set_mol_turnover_rates(n_e_ref, photo_scale)
	real*8, intent(in) :: n_e_ref, photo_scale
	real*8 :: nH, nHe, ne, ntot, s(8), n3ref
	integer :: i, nrow

	nH   = ieq_cell%nh
	nHe  = ieq_cell%nhe
	ne   = n_e_ref
	ntot = ieq_cell%ntot

	! (1) H+ : photo- and collisional ionization of H0, the dissociative H2
	!     photoionization, radiative recombination, the H2 channels
	!     R9/R10/R13, R17, Penning, and the H <-> He charge exchange of both
	!     directions.
	s(1) = (photo_scale*ieq_cell%P_HI + photo_scale*ieq_cell%P_H2_di       &
	        + photo_scale*2.0d0*ieq_cell%P_H2_dd                          &
	        + (ieq_cell%a_ion_HI + ieq_cell%rchiiB)*ne)*nH               &
	     + (mk9 + mk10 + mk13)*nH*nH                                      &
	     + (mk17 + ieq_cell%Q31                                           &
	        + ieq_cell%kcx_He0_Hp + ieq_cell%kcx_Hep_H0)*nHe*nH

	! (2) He+ : the He I and He II photo/collisional/recombination channels
	!     of the summed helium balance, the molecular sinks R17/R23, the
	!     He+ + CO channel and the H <-> He charge exchange. The CO term is
	!     the rate coefficient times the cell's CO reservoir and the whole
	!     helium element, the same bound the other rows are built on; it is
	!     the dominant He+ sink of the shielded molecular base, where CO is
	!     the most abundant heavy species, so leaving it out would scale
	!     that row far below the reaction it is mostly made of.
	s(2) = (photo_scale*ieq_cell%P_HeI + photo_scale*ieq_cell%P_HeII      &
	        + photo_scale*ieq_cell%P_HeITR                                &
	        + (ieq_cell%a_ion_HeI + ieq_cell%a_ion_HeII                   &
	           + ieq_cell%a_ion_HeITR + ieq_cell%rcheiiB                  &
	           + ieq_cell%rcheiiiB + ieq_cell%rcheiTR)*ne)*nHe            &
	     + (mk17 + mk23                                                   &
	        + ieq_cell%kcx_He0_Hp + ieq_cell%kcx_Hep_H0)*nHe*nH
	if (thereis_oxychem)                                                  &
		s(2) = s(2) + rk_D1_Hep_CO()*ieq_cell%n_co*nHe

	! (3) He++ : He II ionization and He III recombination.
	s(3) = (photo_scale*ieq_cell%P_HeII + (ieq_cell%a_ion_HeII            &
	                           + ieq_cell%rcheiiiB)*ne)*nHe

	! (4) H2 : formation by R6/R9/R11/R15 and every destruction channel
	!     (photoionization, Lyman-Werner, R8/R10/R12/R13/R14, He+ and HeH+
	!     reactions, He(2^3S) ionization).
	!     The R12/R15 third body is the collider sum, which is a function of
	!     the composition; the bound puts each collider at the whole of its
	!     element, n(H2) = n_H/2, n(H) = n_H and n(He) = n_He.  No mixture
	!     reaches all three at once, so this is an upper bound on the
	!     collider sum, which is what a turnover scale has to be.
	n3ref = h2_third_body_density(ntot, 0.5d0*nH, nH, nHe)
	s(4) = ((mk6 + mk14)*ne + (mk9 + mk11)*nH + mk15*n3ref*nH)*nH         &
	     + (photo_scale*ieq_cell%P_H2 + photo_scale*ieq_cell%k_LW         &
	        + (mk8 + mk10 + mk13 + mk18)*nH + mk12*n3ref                  &
	        + (mk17 + mk23 + mk_ion_H2)*nHe)*nH

	! (5) H2+ : R10/R11/R19/R23, H2 photoionization and the Penning branch of
	!     He(2^3S)+H2, against R5/R8/R9 and the H2+ + He -> HeH+ + H channel
	!     that replaced R20.  Only the channel that leaves the molecule
	!     intact makes H2+, i.e. the total less the dissociative,
	!     double-ionization and neutral-dissociation parts.
	s(5) = (photo_scale*(ieq_cell%P_H2 - ieq_cell%P_H2_di                 &
	                     - ieq_cell%P_H2_dd - ieq_cell%P_H2_nd)           &
	        + (mk10 + mk11 + mk19)*nH                                     &
	        + (mk23 + f_penning_HeI23S*mk_ion_H2)*nHe)*nH                 &
	     + (mk5*ne + (mk8 + mk9)*nH + mk_h2p_he*nHe)*nH

	! (6) H3+ : R8/R13/R18 against R6/R7 dissociative recombination and R11.
	s(6) = (mk8 + mk13 + mk18)*nH*nH                                      &
	     + ((mk6 + mk7)*ne + mk11*nH)*nH

	! (7) HeH+ : H2+ + He and the associative branch of both He(2^3S)
	!     channels, against R16 dissociative recombination, R18 and R19.
	s(7) = mk_h2p_he*nHe*nH                                               &
	     + (1.0d0 - f_penning_HeI23S)*(ieq_cell%Q31 + mk_ion_H2)*nHe*nH   &
	     + (mk16*ne + (mk18 + mk19)*nH)*nH

	! (8) He 2^3S : populated from He+ recombination and 1^1S excitation,
	!     depopulated by photoionization, A31, de-excitation, electron-impact
	!     ionization and the two Penning collisions.
	if (thereis_HeITR) then
		s(8) = (photo_scale*ieq_cell%P_HeITR + ieq_cell%A31           &
		        + (ieq_cell%rcheiTR + ieq_cell%q13 + ieq_cell%q31a    &
		           + ieq_cell%q31b + ieq_cell%a_ion_HeITR)*ne)*nHe    &
		     + (ieq_cell%Q31 + mk_ion_H2)*nHe*nH
		nrow = 8
	else
		nrow = 7
	endif

	mol_inv_turnover(:) = 1.0d0
	do i = 1,nrow
		if (s(i) .gt. 0.0d0) mol_inv_turnover(i) = 1.0d0/s(i)
	enddo

	end subroutine set_mol_turnover_rates

	! ===============================================================
	!  Oxygen chemistry (the A2 option)
	! ===============================================================

	! Rate coefficients and photolysis rates of the oxygen carriers that are
	! invariant across a cell's Newton solve. Called once per cell from
	! ioniz_eq, next to set_mol_coeffs. j_h2o_band and j_oh_band are the
	! band-resolved photodissociation rates [s^-1] already carrying the
	! star-ward attenuation (water_photolysis.f90); they are summed onto the
	! four channels with the quantum yields of oxygen_rates here, so the
	! residual sees one number per channel.
	subroutine set_oxygen_coeffs(T, j_h2o_band, j_oh_band)
	real*8, intent(in) :: T
	real*8, intent(in) :: j_h2o_band(n_fuv_band), j_oh_band(n_fuv_band)
	integer :: ib

	ok1  = rk_O1_OH_H2_water(T)
	ok1r = rate_from_detailed_balance(ok1, (/ ith_OH, ith_H2 /),          &
	                                       (/ ith_H2O, ith_H /), T)
	ok2  = rk_O2_O_H2_hydroxyl(T)
	ok2r = rate_from_detailed_balance(ok2, (/ ith_O, ith_H2 /),           &
	                                       (/ ith_OH, ith_H /), T)
	ok6  = rk_O6_O1D_H2_hydroxyl()

	oj3 = 0.0d0
	oj4 = 0.0d0
	oj5 = 0.0d0
	oj7 = 0.0d0
	do ib = 1, n_fuv_band
		oj3 = oj3 + qy_H2O_OH_H(ib)  *j_h2o_band(ib)
		oj4 = oj4 + qy_H2O_H2_O1D(ib)*j_h2o_band(ib)
		oj5 = oj5 + qy_H2O_O_H_H(ib) *j_h2o_band(ib)
		oj7 = oj7 + j_oh_band(ib)
	enddo

	end subroutine set_oxygen_coeffs

	! Turnover scale of the two oxygen-carrier rows, appended to the
	! molecular ones (see mol_inv_turnover for what the scale is and why the
	! system needs one). An oxygen row runs as n_O times a rate and the
	! oxygen family is ~1e-4 of the hydrogen, so without this the block sits
	! far below the helium one for the same reason the molecular block does.
	! Called once per cell from ioniz_eq AFTER set_mol_turnover_rates, which
	! resets the whole array. Mirrors oxygen_carrier_rows term by term, with
	! every carrier set to the whole oxygen family and every H partner to the
	! whole hydrogen.
	! photo_scale multiplies the four photolysis channels oj3/oj4/oj5/oj7
	! and nothing else, with the same meaning and the same exactness at 1 as
	! in set_mol_turnover_rates above.
	subroutine set_oxygen_turnover_rates(iox, photo_scale)
	integer, intent(in) :: iox
	real*8,  intent(in) :: photo_scale
	real*8 :: nH, nO, sc

	nH = ieq_cell%nh
	nO = ieq_cell%n_ofam

	! (iox) OH: made by O1 reverse, O2 and the H2O photolysis channels that
	!       end in OH; destroyed by O1, O2 reverse and O7.
	sc = nO*((ok1r + ok2 + ok1 + ok2r)*nH + photo_scale*oj3               &
	         + photo_scale*oj4 + photo_scale*oj7)
	if (sc .gt. 0.0d0) mol_inv_turnover(iox) = 1.0d0/sc

	! (iox+1) H2O: made by O1, destroyed by O1 reverse and by all three
	!       photolysis channels.
	sc = nO*((ok1 + ok1r)*nH + photo_scale*oj3 + photo_scale*oj4          &
	         + photo_scale*oj5)
	if (sc .gt. 0.0d0) mol_inv_turnover(iox+1) = 1.0d0/sc

	! Row 4 (H2) gains the oxygen exchange terms, so its bound gains them
	! too. They are far below the hydrogen terms already in it, so this can
	! only tighten the bound slightly; a missing term would make the scale a
	! weaker bound, never the root wrong.
	sc = (ok1 + ok2 + ok1r + ok2r)*nO*nH
	if (sc .gt. 0.0d0 .and. mol_inv_turnover(4) .gt. 0.0d0)               &
		mol_inv_turnover(4) = 1.0d0/(1.0d0/mol_inv_turnover(4) + sc)

	end subroutine set_oxygen_turnover_rates

	! The two oxygen-carrier balance rows, plus the oxygen cycle's exchange
	! with the H2 row. Shares the sign convention of mol_heh_rows (production
	! positive) and is called AFTER it, so fvec(4) already holds the
	! molecular H2 balance and is added to here.
	!
	! Reactions (docs/a2_reaction_audit.md sec. 2; ids as in that table):
	!   O1   OH  + H2  -> H2O + H      O1r  H2O + H  -> OH + H2
	!   O2   O   + H2  -> OH  + H      O2r  OH  + H  -> O  + H2
	!   O3   H2O + hv  -> OH  + H
	!   O4   H2O + hv  -> H2  + O(1D)  followed by
	!   O6   O(1D)+ H2 -> OH  + H
	!   O5   H2O + hv  -> O   + H + H
	!   O7   OH  + hv  -> O   + H
	!
	! O4 AND O6 ENTER AS ONE CHANNEL, and that is the accurate treatment, not
	! a shortcut. O(1D) has exactly one sink in the audited set, O6, and a
	! chemical lifetime of 1.6e-3 s at the HD 189733 b base -- shorter than
	! every other time scale in the problem by many decades -- so its local
	! steady state is exact and the flux through O6 EQUALS the O(1D)
	! production rate oj4*n_H2O. Writing the pair as one channel therefore
	! adds no approximation and removes the division by n_H2 that the
	! explicit steady-state density would need. Their net effect on the H2
	! budget cancels exactly (O4 makes an H2, O6 consumes it), which is why
	! only O1 and O2 appear in the H2 row below; the measured budget of the
	! design's section 2.2 shows the same cancellation, -7.8% against +7.8%.
	!
	! WHAT THIS FORCES, AND WHERE IT STOPS BEING TRUE. Because O6 is the only
	! sink, every O(1D) is sent to OH + H however little H2 there is. In a
	! gas with no H2 the true fate would be radiative decay to O(3P) (the
	! [O I] 6300/6364 A pair) or collisional quenching, and the branch would
	! end in atomic O instead. Neither rate is in the audited set, so neither
	! is used. The regime is self-limiting rather than dangerous: O(1D) is
	! made only out of H2O, and H2O is made only by OH + H2, so where H2 is
	! gone the source of the channel is gone with it.
	!
	! n_o0 is the FREE ATOMIC oxygen of the cell -- the O I column with the
	! oxygen bound in OH, H2O and CO removed. That redefinition of O I is
	! section 4.2 of the design and is made by the caller.
	!
	! The four photolysis channels enter as DUMMY ARGUMENTS pj3/pj4/pj5/pj7
	! rather than being read from the module state oj3/oj4/oj5/oj7 that
	! set_oxygen_coeffs fills. The cell solves pass exactly that state, so
	! their arithmetic is unchanged; the radiation-field continuation of
	! constrained_chemical_equilibrium passes the same rates scaled by its
	! lambda, which is what lets one definition of these rows serve a solve
	! at any field strength.
	subroutine oxygen_carrier_rows(fvec, iox, n_hi, n_h2, n_oh, n_h2o,    &
	                               n_o0, pj3, pj4, pj5, pj7,               &
	                               p_OH, l_OH, p_H2O, l_H2O,               &
	                               p_H2_oxy, l_H2_oxy)
	real*8 :: fvec(*)
	integer, intent(in) :: iox
	real*8, intent(in)  :: n_hi, n_h2, n_oh, n_h2o, n_o0
	real*8, intent(in)  :: pj3, pj4, pj5, pj7
	! The production and the loss of each row separately, both positive, for
	! the same reason and the same reader mol_heh_rows states.  p_H2_oxy and
	! l_H2_oxy are what the oxygen cycle adds to the H2 row.
	real*8, optional, intent(out) :: p_OH, l_OH, p_H2O, l_H2O
	real*8, optional, intent(out) :: p_H2_oxy, l_H2_oxy
	real*8 :: prod_oh, loss_oh_o1, loss_oh_o2r, loss_oh_j7
	real*8 :: prod_h2o, loss_h2o_o1r, loss_h2o_j
	real*8 :: prod_h2_oxy, loss_h2_oxy

	! (iox) OH balance
	prod_oh     = ok1r*n_h2o*n_hi + ok2*n_o0*n_h2                        &
	            + (pj3 + pj4)*n_h2o
	loss_oh_o1  = ok1*n_oh*n_h2
	loss_oh_o2r = ok2r*n_oh*n_hi
	loss_oh_j7  = pj7*n_oh
	fvec(iox)   = prod_oh - loss_oh_o1 - loss_oh_o2r - loss_oh_j7
	if (present(p_OH)) p_OH = prod_oh
	if (present(l_OH)) l_OH = loss_oh_o1 + loss_oh_o2r + loss_oh_j7

	! (iox+1) H2O balance
	prod_h2o     = ok1*n_oh*n_h2
	loss_h2o_o1r = ok1r*n_h2o*n_hi
	loss_h2o_j   = (pj3 + pj4 + pj5)*n_h2o
	fvec(iox+1)  = prod_h2o - loss_h2o_o1r - loss_h2o_j
	if (present(p_H2O)) p_H2O = prod_h2o
	if (present(l_H2O)) l_H2O = loss_h2o_o1r + loss_h2o_j

	! H2 exchange with the oxygen cycle, added to the molecular H2 balance.
	! O1 and O2 consume an H2 each; their reverses give one back. O3, O5 and
	! O7 do not touch H2, and the O4/O6 pair cancels (see above).
	loss_h2_oxy = (ok1*n_oh + ok2*n_o0)*n_h2
	prod_h2_oxy = (ok1r*n_h2o + ok2r*n_oh)*n_hi
	fvec(4) = fvec(4)                                                    &
	        - loss_h2_oxy                                                &
	        + prod_h2_oxy
	if (present(p_H2_oxy)) p_H2_oxy = prod_h2_oxy
	if (present(l_H2_oxy)) l_H2_oxy = loss_h2_oxy

	end subroutine oxygen_carrier_rows

	! O(1D) number density [cm^-3] from its local steady state, for the
	! diagnostic output only: production oj4*n_H2O against the single sink
	! O6. Nothing in the network reads it (see oxygen_carrier_rows), so the
	! floor on n_H2 below cannot feed back into the solution.
	double precision function excited_oxygen_density(n_h2, n_h2o)         &
	                          result(n_o1d)
	real*8, intent(in) :: n_h2, n_h2o
	real*8 :: sink
	sink = ok6*max(n_h2, 0.0d0)
	if (sink .le. 0.0d0) then
		n_o1d = 0.0d0
	else
		n_o1d = oj4*max(n_h2o, 0.0d0)/sink
	endif
	end function excited_oxygen_density

	subroutine ion_system_HeH_mol(Neq,x,fvec,iflag,params)

	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: params(40)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR,g_h2,g_h2_di,g_h2_dd,g_h2_nd,g_lw
	real*8  :: b_hi,b_hei,b_heii,b_heiTR
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR
	real*8  :: A31,q13,q31a,q31b,Q31
	real*8  :: n_h,n_he,n_e,T,ntot
	real*8  :: n_hi,n_hii,n_h2,n_h2p,n_h3p,n_hehp
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI

	g_hi    = ieq_cell%P_HI
	g_hei   = ieq_cell%P_HeI
	g_heii  = ieq_cell%P_HeII
	a_hii   = ieq_cell%rchiiB
	a_heii  = ieq_cell%rcheiiB
	a_heiii = ieq_cell%rcheiiiB
	n_h     = ieq_cell%nh
	n_he    = ieq_cell%nhe
	b_hi    = ieq_cell%a_ion_HI
	b_hei   = ieq_cell%a_ion_HeI
	b_heii  = ieq_cell%a_ion_HeII
	b_heiTR = ieq_cell%a_ion_HeITR   ! He 2^3S collisional ioniz. (0 if no triplet)
	a_heiTR = ieq_cell%rcheiTR
	A31     = ieq_cell%A31
	g_heiTR = ieq_cell%P_HeITR
	q13     = ieq_cell%q13
	q31a    = ieq_cell%q31a
	q31b    = ieq_cell%q31b
	Q31     = ieq_cell%Q31
	g_h2    = ieq_cell%P_H2
	g_h2_di = ieq_cell%P_H2_di
	g_h2_dd = ieq_cell%P_H2_dd   ! double ionization (0 unless a model is on)
	g_h2_nd = ieq_cell%P_H2_nd   ! neutral dissociation (0 unless on)
	g_lw    = ieq_cell%k_LW      ! Lyman-Werner photodissociation (0 if off)
	T       = ieq_cell%T_K
	ntot    = ieq_cell%ntot

	! Species densities
	n_hi   = (1.0d0 - x(1) - x(4) - x(5) - x(6) - x(7))*n_h
	n_hii  = x(1)*n_h
	n_h2   = 0.5d0*x(4)*n_h
	n_h2p  = 0.5d0*x(5)*n_h
	n_h3p  = x(6)*n_h/3.0d0
	n_hehp = x(7)*n_h
	n_hei   = (1.0d0 - x(2) - x(3))*n_he - n_hehp   ! free neutral He
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	if (thereis_HeITR) then
		n_heiTR = x(8)*n_he
	else
		n_heiTR = 0.0d0
	endif
	n_heiSI = n_hei - n_heiTR

	! Electron density (each molecular ion carries +1)
	n_e = n_hii + n_h2p + n_h3p + n_hehp + n_heii + 2.0d0*n_heiii

	call mol_heh_rows(fvec, n_hi, n_hii, n_h2, n_h2p, n_h3p, n_hehp,   &
	                  n_heiSI, n_heiTR, n_heii, n_heiii, n_e, ntot,     &
	                  g_hi, g_hei, g_heii, g_heiTR, g_h2, g_h2_di,      &
	                  g_h2_dd, g_h2_nd, g_lw,                           &
	                  a_hii, a_heii, a_heiii, a_heiTR,                  &
	                  b_hi, b_hei, b_heii, b_heiTR,                     &
	                  q13, q31a, q31b, Q31, A31)

	! He <-> H charge exchange (Huang Table 4 group B). Row 1 (H+ balance)
	! and row 2 (He+ balance) are both written production positive, so
	! he_row_sign = +1. The ground singlet n_heiSI is the CX He I reservoir.
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,   &
	                  n_hi, n_hii, n_heiSI, n_heii, 1.0d0)

	! Each row divided by its own turnover rate (set_mol_turnover_rates), so
	! the helium and molecular blocks reach hybrd1 with the same weight.
	fvec(1:Neq) = fvec(1:Neq)*mol_inv_turnover(1:Neq)

	! Where the H2 partition is imposed rather than solved -- transported
	! carriers, or the lower-boundary reservoir composition -- this row is
	! handed the answer (see the same block in System_HeH_mol_metals, which
	! is the system the oxygen chemistry actually reaches -- oxygen is a
	! metal element, so a run without metals is refused). Applied after the
	! turnover scaling so the row is exactly x - x_fix.
	if (ieq_cell%x_h2_fixed) fvec(4) = x(4) - ieq_cell%x_h2_fix
	! The transported ionization fractions, where the flow carries them
	! and not this cell's local balance (ion_residual_core). The H2
	! partition just above is imposed for its own reasons and keeps its
	! own flag.
	call impose_transported_ionization_fractions(ieq_cell, x, fvec)

	return
	end subroutine ion_system_HeH_mol

	!----------------------------------!

	! The eight balance rows of the molecular network, shared verbatim by
	! System_HeH_mol and System_HeH_mol_metals. n_e is an INPUT so that the
	! metal-coupled system can pass its metal-inclusive electron density; the
	! molecular rate coefficients are read from this module's own cell state
	! (set_mol_coeffs). Rows 1-7 are always written; row 8 only when the He
	! 2^3S metastable is tracked. Only explicit-shape / assumed-size dummies
	! are used, as in ion_residual_core.
	subroutine mol_heh_rows(fvec, n_hi, n_hii, n_h2, n_h2p, n_h3p, n_hehp,  &
	                        n_heiSI, n_heiTR, n_heii, n_heiii, n_e, ntot,    &
	                        g_hi, g_hei, g_heii, g_heiTR, g_h2, g_h2_di,     &
	                        g_h2_dd, g_h2_nd, g_lw,                          &
	                        a_hii, a_heii, a_heiii, a_heiTR,                 &
	                        b_hi, b_hei, b_heii, b_heiTR,                    &
	                        q13, q31a, q31b, Q31, A31,                       &
	                        p_Hp, l_Hp, p_H2, l_H2, l_H2_phot, h2_chan)

	real*8 :: fvec(*)
	real*8, intent(in) :: n_hi,n_hii,n_h2,n_h2p,n_h3p,n_hehp
	real*8, intent(in) :: n_heiSI,n_heiTR,n_heii,n_heiii,n_e,ntot
	real*8, intent(in) :: g_hi,g_hei,g_heii,g_heiTR,g_h2,g_lw
	! The three resolved product channels of g_h2, each a SUBSET of it and
	! disjoint from the other two: g_h2_di leaves H + H+ + e-, g_h2_dd
	! leaves H+ + H+ + 2e- and g_h2_nd leaves H + H with no ion. g_h2 stays
	! the whole H2 destruction rate, and what is left of it after the three
	! are removed is the channel that makes H2+. g_h2_dd and g_h2_nd are
	! zero unless their options are on.
	real*8, intent(in) :: g_h2_di, g_h2_dd, g_h2_nd
	real*8, intent(in) :: a_hii,a_heii,a_heiii,a_heiTR
	real*8, intent(in) :: b_hi,b_hei,b_heii,b_heiTR
	real*8, intent(in) :: q13,q31a,q31b,Q31,A31
	! THE ROWS SPLIT INTO THEIR PRODUCTION AND THEIR LOSS, both positive,
	! for the two rows a transported carrier owns: the proton (1) and H2 (4).
	! A row is production - loss, and where the chemistry is fast the two
	! can stand orders above their difference, so a reader of the row's
	! imbalance cannot reconstruct them from it.  They are handed out here,
	! at the one place the terms are written, for the cell-by-cell record of
	! a carrier row (EXHALE_CARRIER_ROW_TERMS, diffusive_photochemistry.f90)
	! -- a diagnostic, and nothing in the solution reads them.  The rows
	! below are assembled FROM these names, in the order and grouping they
	! were written in, so the residual is unchanged to the last bit and the
	! row and its terms cannot drift apart.  l_H2_phot is the part of the H2
	! loss carried by the radiation field alone, photodissociation and
	! photoionization together.
	real*8, optional, intent(out) :: p_Hp, l_Hp, p_H2, l_H2, l_H2_phot
	! THE H2 ROW CHANNEL BY CHANNEL, each a volumetric rate [cm^-3 s^-1],
	! in the order the named parameters n_h2chan_* of
	! diffusive_photochemistry give: four formation channels, then the nine
	! destruction channels, every one of them a single term of prod_h2 or
	! loss_h2 above and written from the same factors, so the channels sum
	! to the production and the loss the row is built from.  Diagnostic
	! only (EXHALE_CARRIER_ROW_TERMS=1); nothing in the solution reads it.
	real*8, optional, intent(out) :: h2_chan(*)
	real*8  :: k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15
	real*8  :: k16,k17,k18,k19,k_h2p_he,k23,k_ion_H2
	! Rate at which one He+ ion is destroyed by the cell's CO, k_D1 n_CO
	! [s^-1]; see row (2). Zero without the oxygen chemistry, which is the
	! only option that carries CO.
	real*8  :: k_co_hep
	! Production and loss of the two carrier rows, named so that the row
	! and its scale are built from one expression (declarations above).
	real*8  :: prod_hp, loss_hp_rec, loss_hp_h2, prod_h2, loss_h2
	! H2-equivalent third-body density of the R12/R15 pair [cm^-3].
	real*8  :: n_third

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
	if (present(p_Hp)) p_Hp = prod_hp
	if (present(l_Hp)) l_Hp = loss_hp_rec + loss_hp_h2

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
	if (present(p_H2)) p_H2 = prod_h2
	if (present(l_H2)) l_H2 = loss_h2
	! The radiative part of the same loss: g_h2 is the whole photo-
	! destruction of H2 and g_lw the Lyman-Werner dissociation.
	if (present(l_H2_phot)) l_H2_phot = (g_h2 + g_lw)*n_h2
	! The same two expressions, term by term.  Grouped exactly as they are
	! grouped above, so no channel is a second statement of a rate.
	if (present(h2_chan)) then
		h2_chan(1)  = k15*n_hi*n_hi
		h2_chan(2)  = k9*n_h2p*n_hi
		h2_chan(3)  = k6*n_e*n_h3p
		h2_chan(4)  = k11*n_h3p*n_hi
		h2_chan(5)  = g_h2*n_h2
		h2_chan(6)  = g_lw*n_h2
		h2_chan(7)  = (k10 + k13)*n_hii*n_h2
		h2_chan(8)  = k12*n_third*n_h2
		h2_chan(9)  = k14*n_e*n_h2
		h2_chan(10) = k8*n_h2p*n_h2
		h2_chan(11) = (k17 + k23)*n_heii*n_h2
		h2_chan(12) = k18*n_hehp*n_h2
		h2_chan(13) = k_ion_H2*n_heiTR*n_h2
	endif

	! (5) H2+ balance
	! (g_h2 - g_h2_di - g_h2_dd - g_h2_nd): only the branch that leaves the
	! molecule bound makes an H2+ ion. The two ionizing branches removed
	! here are proton sources in row (1); the neutral one makes no ion at
	! all. Both of the latter two are zero unless their options are on.
	! f_penning_HeI23S*k_ion_H2*n_heiTR*n_h2: H2+ produced by the Penning
	! branch of He(2^3S)+H2 ionization. The associative branch of the same
	! collision gives H + HeH+ + e- and appears in row (7) instead.
	! k_h2p_he*n_heiSI: H2+ + He -> HeH+ + H, the HeH+ source the network
	! carries since item L7f retired Koskinen R20 (mol_rates: the paper R20
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
	! Koskinen R20, which item L7f retired because the paper Table 1 cites
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
		                    g_heiTR, a_heiTR, q13, q31a, q31b, Q31, A31,    &
		                    b_heiTR)
		! He(2^3S)+H2 ionization triplet loss, total of the Penning and the
		! associative branch: both quench the metastable, and neither makes
		! He+, so there is no He+ row term. Garcia Munoz (2025) Table A.5.
		fvec(8) = fvec(8) - k_ion_H2*n_heiTR*n_h2
	endif

	end subroutine mol_heh_rows

	! End of module
	end module System_HeH_mol

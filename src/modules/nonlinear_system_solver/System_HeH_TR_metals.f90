	module System_HeH_TR_metals
	! Coupled ionization-equilibrium system for H, He WITH the He 2^3S
	! metastable triplet (HeITR) AND an arbitrary number of trace metals.
	!
	! This is the merge of System_HeH_TR (H/He + triplet) and
	! System_HeH_metals (H/He + metals); the two were historically mutually
	! exclusive because both placed their extra unknown at x(4). Here the
	! He-triplet fraction keeps x(4) and the metals are shifted up by one:
	!
	!   x(1) = n_HII /n_H
	!   x(2) = n_HeII/n_He
	!   x(3) = n_HeIII/n_He
	!   x(4) = n_HeITR/n_He                         (He 2^3S, neutral)
	! and, for each metal element e = 1..met_nelem (canonical order),
	!   x(5 + 2*(e-1)) = n_(Xe)II /n_Xe
	!   x(6 + 2*(e-1)) = n_(Xe)III/n_Xe             (pinned 0 if met_top(e)<2)
	!
	! so N_eq = 4 + 2*met_nelem.
	!
	! Design choices (documented for validation):
	!  * The H/He/triplet rows (1-4) are taken VERBATIM from System_HeH_TR,
	!    so with met_nelem = 0 this system reduces exactly to the validated
	!    triplet solver. Those rows carry electron-impact ionization of
	!    H0/He(1^1S)/He+ and of the He 2^3S metastable, matching the metal
	!    rows; the H/He collisional ionization is negligible at the ~1e4 K
	!    wind temperatures but is kept for internal consistency.
	!  * The electron density couples the two blocks: it is the H/He sum PLUS
	!    the metal charges (X+ once, X++ twice). The triplet is a neutral
	!    excited state and does NOT contribute to n_e.
	!  * Metal coefficients for each cell come through the System_HeH_metals
	!    module arrays (set by set_metal_coeffs); charge exchange is added by
	!    cx_add_to_fvec with cx_metal_base = 5 (metals shifted to row 5).
	!
	! Solved with hybrd1 (numerical Jacobian), matching the triplet system;
	! no analytic Jacobian is provided for the merged residual.

	use global_parameters
	use ion_cell_state, only: ieq_cell
	use charge_exchange,    only: cx_add_to_fvec, he_h_cx_fvec
	use ion_residual_core,  only: metal_fractions, metal_electron_sum,    &
	                              metal_rows, heh_tr_rows,             &
	                              impose_transported_ionization_fractions
	use System_HeH_metals,  only: met_nelem, met_ntot, met_g0, met_g1,    &
	                              met_b0, met_b1, met_a1, met_a2, met_top, met_g02

	implicit none

	contains

	subroutine ion_system_HeH_TR_metals(N_eq,x,fvec,iflag,params)

	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: params(60)

	! H/He/triplet coefficients (params 1-18, same layout as System_HeH_TR)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR        ! photoionization
	real*8  :: b_hi,b_hei,b_heii,b_heiTR        ! collisional ionization
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR     ! recombination
	real*8  :: A31,q13,q31g,q31a,q31b,Q31            ! triplet kinetics
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_heii,n_heiii,n_heiTR,n_heiSI
	! Each element's metal densities (neutral/+/++); charge exchange added later.
	real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)

	! Unpack H/He coefficients
	g_hi    = ieq_cell%P_HI       ! = P_HI
	g_hei   = ieq_cell%P_HeI      ! = P_HeI
	g_heii  = ieq_cell%P_HeII     ! = P_HeII
	a_hii   = ieq_cell%rchiiB     ! = rchiiB
	a_heii  = ieq_cell%rcheiiB    ! = rcheiiB
	a_heiii = ieq_cell%rcheiiiB   ! = rcheiiiB
	n_h     = ieq_cell%nh         ! = nh
	n_he    = ieq_cell%nhe        ! = nhe
	b_hi    = ieq_cell%a_ion_HI      ! = a_ion_HI
	b_hei   = ieq_cell%a_ion_HeI     ! = a_ion_HeI
	b_heii  = ieq_cell%a_ion_HeII    ! = a_ion_HeII
	b_heiTR = ieq_cell%a_ion_HeITR   ! = a_ion_HeITR (He 2^3S collisional ioniz.)

	! Triplet parameters
	a_heiTR = ieq_cell%rcheiTR   ! = rcheiTR
	A31     = ieq_cell%A31       ! = A31
	g_heiTR = ieq_cell%P_HeITR   ! = P_HeITR
	q13     = ieq_cell%q13       ! = q13
	q31g    = ieq_cell%q31g      ! = q31g (reverse of q13)
	q31a    = ieq_cell%q31a      ! = q31a
	q31b    = ieq_cell%q31b      ! = q31b
	Q31     = ieq_cell%Q31       ! = Q31

	! H/He densities from fractions
	n_hi    = (1.0 - x(1))*n_h
	n_hii   = x(1)*n_h
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	n_heiSI = (1.0 - x(2) - x(3) - x(4))*n_he   ! singlet ground only
	n_heiTR = x(4)*n_he                         ! 2^3S triplet

	! Metal ion densities (neutral/+/++) from fractions, canonical order.
	! Metals occupy rows 5.. (shifted up by one from System_HeH_metals).
	call metal_fractions(x, 5, met_nelem, met_ntot, nm0, nm1, nm2)

	! Electron density: H/He first (reproducing the original sum), then the
	! metal charges. The triplet is neutral and contributes nothing.
	n_e = n_hii + n_heii + 2.0*n_heiii
	call metal_electron_sum(n_e, met_nelem, nm1, nm2)

	! --- H/He/triplet rows (verbatim System_HeH_TR, n_e now metal-inclusive) ---
	call heh_tr_rows(fvec, n_hi, n_hii, n_heiSI, n_heiTR, n_heii, n_heiii,  &
	                 n_e, g_hi, g_hei, g_heii, g_heiTR,                      &
	                 a_hii, a_heii, a_heiii, a_heiTR,                        &
	                 b_hi, b_hei, b_heii, b_heiTR,                           &
	                 q13, q31g, q31a, q31b, Q31, A31)

	! --- Metal rows (verbatim System_HeH_metals, shifted to rows 5..) ---
	call metal_rows(fvec, x, 5, met_nelem, met_ntot, met_g0, met_g1,     &
	                met_g02,                                           &
	                met_b0, met_b1, met_a1, met_a2, met_top,             &
	                nm0, nm1, nm2, n_e)

	! Charge exchange (Huang Table 4) on the H, He and metal rows. The metal
	! rows are at base 5 here; the driver sets cx_metal_base = 5 before this
	! solve. Absent reactants contribute zero, preserving the identity rows.
	! fvec(2) is the SUMMED He I balance of heh_tr_rows, written He I-gain
	! positive, so he_row_sign = -1: the group C reactions metal + He / He+
	! destroy He II and make He I, which is a gain of this row. The H row,
	! the He II <-> He III row and the metal rows are ionization positive and
	! take the generic sign. The rate coefficients are the ones cx_set_cell
	! stored for this cell, whose temperature ieq_cell%T_K is.
	! The helium reactant of the group C metal + He reactions is the GROUND
	! SINGLET n_heiSI: those rates are ground-state rates (see the group C
	! paragraph of charge_exchange), so the metastable must not be charged
	! to them, and the singlet loss they report is a loss of the summed
	! He I this row balances.
	call cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                       &
	                    n_hi, n_hii, n_heiSI, n_heii, n_heiii,           &
	                    -1.0d0, ieq_cell%T_K)

	! He <-> H charge exchange (Huang Table 4 group B). The He reactant of
	! He + H+ -> He+ + H is the GROUND SINGLET He(1^1S), n_heiSI: the rate
	! charge_exchange::he_h_cx_rates forms for it, the detailed-balance
	! reverse of He+ + H -> He(1^1S) + H+, carries the barrier
	! exp(-12.75/T4), and 12.75e4 K = 10.99 eV is the ionization-potential
	! difference 24.587 - 13.598 eV of ground-state helium against hydrogen.
	! He(2^3S) lies 19.82 eV above the singlet, so its own charge exchange
	! with H+ is exothermic and has no such barrier; it is a different
	! reaction with a different rate, and neither this system nor the
	! metastable balance tr_triplet_row carries it.
	! The summed He I row (fvec 2) is written HeI-gain positive here, so
	! he_row_sign = -1, and the singlet loss it reports is a loss of the sum.
	! Group B is excluded from cx_act, so it is applied only here (no double
	! counting with cx_add_to_fvec).
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,    &
	                  ieq_cell%kcx_Hepp_H0,                              &
	                  n_hi, n_hii, n_heiSI, n_heii, n_heiii, -1.0d0,      &
	                  .false.)

	! The transported ionization fractions, where the flow carries them
	! and not this cell's local balance (ion_residual_core).
	call impose_transported_ionization_fractions(ieq_cell, x, fvec)

	return

	! End of subroutine
	end subroutine ion_system_HeH_TR_metals

	! End of module
	end module System_HeH_TR_metals

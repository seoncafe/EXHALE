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
	! Rows are steady-state production-loss balances.  The atomic rows
	! use EXHALE's own rate arrays (P_HI/rchiiB/Voronov...) so the
	! molecular-free limit reproduces the atomic systems' solution; the
	! He 2^3S row follows System_HeH_TR plus the He(2^3S)+H2 Penning loss
	! (Garcia Munoz 2025 Table A.5), the dominant metastable sink toward an
	! H2-rich base.  Molecular channels are
	! the Koskinen et al. (2022) Table-1 network via mol_rates (R21/R22
	! H-He charge exchange excluded to preserve the atomic limit; the
	! Lyman-Werner photodissociation caveat is inherited -- see mol_rates).
	!
	! params layout (1-18 identical to System_HeH_TR):
	!   1 P_HI  2 P_HeI  3 P_HeII  4 rchiiB  5 rcheiiB  6 rcheiiiB
	!   7 n_h(nuclei)  8 n_he  9 a_ion_HI  10 a_ion_HeI  11 a_ion_HeII
	!   12 rcheiTR  13 A31  14 P_HeITR  15 q13  16 q31a  17 q31b  18 Q31
	!   19 P_H2 (photoionization rate coefficient of H2, s^-1)
	!   20 T [K]   21 n_tot (total particle density, for 3-body M)

	use global_parameters, only: thereis_HeITR
	use mol_rates
	use ion_residual_core, only: tr_triplet_row
	use ion_cell_state, only: ieq_cell
	use charge_exchange, only: he_h_cx_fvec
	use Cooling_Coefficients, only: penning_HeI23S_H2   ! He(2^3S)+H2 Penning rate

	implicit none

	! Molecular rate coefficients, invariant across a cell's Newton solve
	! (depend only on T and n_tot). Computed once per cell by set_mol_coeffs
	! and read by the residual, mirroring set_metal_coeffs / cx_set_cell.
	real*8, save :: mk5,mk6,mk7,mk8,mk9,mk10,mk11,mk12,mk13,mk14,mk15
	real*8, save :: mk16,mk17,mk18,mk19,mk20,mk23
	! He(2^3S)+H2 Penning ionization rate coefficient (Garcia Munoz 2025
	! Table A.5); depends only on T. Zero-effect unless the triplet is present.
	real*8, save :: mk_pen_H2
	!$omp threadprivate(mk5,mk6,mk7,mk8,mk9,mk10,mk11,mk12,mk13,mk14,mk15, &
	!$omp                mk16,mk17,mk18,mk19,mk20,mk23,mk_pen_H2)

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
	mk15 = rk_R15_3body_H2(T, ntot)
	mk16 = rk_R16_HeHp_dr(T)
	mk17 = rk_R17_Hep_H2_diss(T)
	mk18 = rk_R18_HeHp_H2()
	mk19 = rk_R19_HeHp_H()
	mk20 = rk_R20_Hep_H2_HeHp()
	mk23 = rk_R23_H2_Hep_cx()
	mk_pen_H2 = penning_HeI23S_H2(T)   ! He(2^3S)+H2 Penning ionization
	end subroutine set_mol_coeffs

	subroutine ion_system_HeH_mol(Neq,x,fvec,iflag,params)

	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: params(40)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR,g_h2
	real*8  :: b_hi,b_hei,b_heii,b_heiTR
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR
	real*8  :: A31,q13,q31a,q31b,Q31
	real*8  :: n_h,n_he,n_e,T,ntot
	real*8  :: n_hi,n_hii,n_h2,n_h2p,n_h3p,n_hehp
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	real*8  :: k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15
	real*8  :: k16,k17,k18,k19,k20,k23,k_pen_H2

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
	k15 = mk15
	k16 = mk16
	k17 = mk17
	k18 = mk18
	k19 = mk19
	k20 = mk20
	k23 = mk23
	k_pen_H2 = mk_pen_H2   ! He(2^3S)+H2 -> He(1^1S)+H2+ + e- (0 without triplet)

	! (1) H+ balance
	! Q31*n_heiTR*n_hi: Penning ionization He(2^3S)+H0 -> He(1^1S)+H+ + e-
	! (n_heiTR is 0 when thereis_HeITR is false, so the term is unconditional).
	fvec(1) = (g_hi + b_hi*n_e)*n_hi                                  &
	        + k9*n_h2p*n_hi + k17*n_heii*n_h2                         &
	        + Q31*n_heiTR*n_hi                                        &
	        - a_hii*n_e*n_hii - (k10 + k13)*n_hii*n_h2

	! (2) He+ balance (atomic part consistent with System_HeH_TR rows
	!     2+3 combined; + molecular sinks R17/R20/R23). b_heiTR*n_e*n_heiTR is
	!     electron-impact ionization of the He 2^3S metastable into He+.
	fvec(2) = (g_hei + b_hei*n_e)*n_heiSI + g_heiTR*n_heiTR           &
	        + b_heiTR*n_e*n_heiTR                                     &
	        + a_heiii*n_e*n_heiii                                     &
	        - (a_heii + a_heiTR)*n_e*n_heii                           &
	        - (g_heii + b_heii*n_e)*n_heii                            &
	        - (k17 + k20 + k23)*n_heii*n_h2

	! (3) He++ balance (verbatim atomic form + collisional ionization)
	fvec(3) = (g_heii + b_heii*n_e)*n_heii - a_heiii*n_e*n_heiii

	! (4) H2 balance
	! k_pen_H2*n_heiTR: He(2^3S)+H2 -> He(1^1S)+H2+ + e- Penning loss of H2
	! (n_heiTR is 0 when thereis_HeITR is false, so the term is unconditional).
	fvec(4) = k6*n_e*n_h3p + k9*n_h2p*n_hi + k11*n_h3p*n_hi           &
	        + k15*n_hi*n_hi                                           &
	        - ( g_h2 + (k10 + k13)*n_hii + k12*ntot + k14*n_e         &
	          + k8*n_h2p + (k17 + k20 + k23)*n_heii + k18*n_hehp      &
	          + k_pen_H2*n_heiTR )*n_h2

	! (5) H2+ balance
	! k_pen_H2*n_heiTR*n_h2: H2+ produced by He(2^3S)+H2 Penning ionization.
	fvec(5) = g_h2*n_h2 + k10*n_hii*n_h2 + k11*n_h3p*n_hi             &
	        + k19*n_hehp*n_hi + k23*n_heii*n_h2                       &
	        + k_pen_H2*n_heiTR*n_h2                                   &
	        - (k5*n_e + k8*n_h2 + k9*n_hi)*n_h2p

	! (6) H3+ balance
	fvec(6) = k8*n_h2p*n_h2 + k13*n_hii*n_h2 + k18*n_hehp*n_h2        &
	        - ((k6 + k7)*n_e + k11*n_hi)*n_h3p

	! (7) HeH+ balance
	fvec(7) = k20*n_heii*n_h2                                         &
	        - (k16*n_e + k18*n_h2 + k19*n_hi)*n_hehp

	! (8) He 2^3S balance (VERBATIM System_HeH_TR row 4)
	if (thereis_HeITR) then
		call tr_triplet_row(fvec(8), n_hi, n_heiSI, n_heiTR, n_heii, n_e,   &
		                    g_heiTR, a_heiTR, q13, q31a, q31b, Q31, A31,    &
		                    b_heiTR)
		! He(2^3S)+H2 Penning ionization triplet loss (product is ground
		! He I, so no He+ row term). Garcia Munoz (2025) Table A.5.
		fvec(8) = fvec(8) - k_pen_H2*n_heiTR*n_h2
	endif

	! He <-> H charge exchange (Huang Table 4 group B). Row 1 (H+ balance)
	! and row 2 (He+ balance) are both written production positive, so
	! he_row_sign = +1. The ground singlet n_heiSI is the CX He I reservoir.
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,   &
	                  n_hi, n_hii, n_heiSI, n_heii, 1.0d0)

	return
	end subroutine ion_system_HeH_mol

	! End of module
	end module System_HeH_mol

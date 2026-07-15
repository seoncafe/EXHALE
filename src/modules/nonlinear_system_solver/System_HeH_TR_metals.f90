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
	!    triplet solver. Consequently the H/He rows omit collisional
	!    ionization (the Oklopcic triplet formulation does), whereas the
	!    metal rows keep it (System_HeH_metals form). Collisional ionization
	!    of H/He is negligible at the ~1e4 K wind temperatures, so this
	!    asymmetry is immaterial; it is kept to preserve each validated block.
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
	use params_idx, only: IPE_PHI, IPE_PHEI, IPE_PHEII, IPE_AHII, IPE_AHEII,   &
	                      IPE_AHEIII, IPE_NH, IPE_NHE, IPE_ATR, IPE_A31,        &
	                      IPE_PTR, IPE_Q13, IPE_Q31A, IPE_Q31B, IPE_Q31
	use charge_exchange,    only: cx_add_to_fvec
	use ion_residual_core,  only: metal_fractions, metal_electron_sum,    &
	                              metal_rows
	use System_HeH_metals,  only: met_nelem, met_ntot, met_g0, met_g1,    &
	                              met_b0, met_b1, met_a1, met_a2, met_top

	implicit none

	contains

	subroutine ion_system_HeH_TR_metals(N_eq,x,fvec,iflag,params)

	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: params(60)

	! H/He/triplet coefficients (params 1-18, same layout as System_HeH_TR)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR        ! photoionization
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR     ! recombination
	real*8  :: A31,q13,q31a,q31b,Q31            ! triplet kinetics
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	! Each element's metal densities (neutral/+/++); charge exchange added later.
	real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)
	real*8  :: n_X
	integer :: e,ix

	! Unpack H/He coefficients
	g_hi    = params(IPE_PHI)     ! = P_HI
	g_hei   = params(IPE_PHEI)     ! = P_HeI
	g_heii  = params(IPE_PHEII)     ! = P_HeII
	a_hii   = params(IPE_AHII)     ! = rchiiB
	a_heii  = params(IPE_AHEII)     ! = rcheiiB
	a_heiii = params(IPE_AHEIII)     ! = rcheiiiB
	n_h     = params(IPE_NH)     ! = nh
	n_he    = params(IPE_NHE)     ! = nhe
	! params(9-11) = a_ion_HI/HeI/HeII are intentionally unused here: the TR
	! H/He balance omits collisional ionization (see header note).

	! Triplet parameters
	a_heiTR = params(IPE_ATR)    ! = rcheiTR
	A31     = params(IPE_A31)    ! = A31
	g_heiTR = params(IPE_PTR)    ! = P_HeITR
	q13     = params(IPE_Q13)    ! = q13
	q31a    = params(IPE_Q31A)    ! = q31a
	q31b    = params(IPE_Q31B)    ! = q31b
	Q31     = params(IPE_Q31)    ! = Q31

	! H/He densities from fractions
	n_hi    = (1.0 - x(1))*n_h
	n_hii   = x(1)*n_h
	n_hei   = (1.0 - x(2) - x(3))*n_he          ! total HeI (singlet + triplet)
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

	! HI <-> HII
	fvec(1) = n_hi*g_hi - a_hii*n_hii*n_e

	! HeI singlet ground (recombination into singlet + triplet, minus
	! photoionization of the singlet ground and the triplet)
	fvec(2) =   n_heii*(a_heiTR + a_heii)*n_e                            &
	          - n_heiSI*g_hei                                           &
	          - n_heiTR*g_heiTR

	! HeII <-> HeIII
	fvec(3) = n_heii*g_heii - a_heiii*n_heiii*n_e

	! He 2^3S triplet balance
	fvec(4) = - n_heiTR*g_heiTR                                         &
	          + n_e*( n_heii*a_heiTR                                    &
	                + n_heiSI*q13                                       &
	                - n_heiTR*(q31a + q31b))                            &
	          - n_heiTR*(A31 + n_hi*Q31)

	! --- Metal rows (verbatim System_HeH_metals, shifted to rows 5..) ---
	call metal_rows(fvec, x, 5, met_nelem, met_ntot, met_g0, met_g1,     &
	                met_b0, met_b1, met_a1, met_a2, met_top,             &
	                nm0, nm1, nm2, n_e)

	! Charge exchange (Huang Table 4) on the H, He and metal rows. The metal
	! rows are at base 5 here; the driver sets cx_metal_base = 5 before this
	! solve. Absent reactants contribute zero, preserving the identity rows.
	call cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                       &
	                    n_hi, n_hii, n_hei, n_heii, n_heiii)

	return

	! End of subroutine
	end subroutine ion_system_HeH_TR_metals

	! End of module
	end module System_HeH_TR_metals

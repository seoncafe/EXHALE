	module System_HeH_mol_metals
	! Coupled ionization-equilibrium system for H and He with the molecular
	! network (H2, H2+, H3+, HeH+), optionally the He 2^3S metastable, AND an
	! arbitrary number of trace metals.
	!
	! This is the merge of System_HeH_mol (H/He + molecules + optional
	! triplet) and System_HeH_metals (H/He + metals). It follows the merge
	! already made for the triplet and the metals (System_HeH_TR_metals): the
	! molecular unknowns keep their rows and the metals are appended above
	! them, so both blocks are written once and shared verbatim.
	!
	!   x(1) = n_HII  /n_H(nuclei)
	!   x(2) = n_HeII /n_He        x(3) = n_HeIII/n_He
	!   x(4) = 2 n_H2 /n_H         x(5) = 2 n_H2+/n_H
	!   x(6) = 3 n_H3+/n_H         x(7) = n_HeH+ /n_H
	!   x(8) = n_HeITR/n_He        (only when thereis_HeITR)
	! and, for each metal element e = 1..met_nelem (canonical order),
	!   x(mbase   + 2*(e-1)) = n_(Xe)II /n_Xe
	!   x(mbase+1 + 2*(e-1)) = n_(Xe)III/n_Xe    (pinned 0 if met_top(e) < 2)
	! with mbase = metal_row_base() = 8, or 9 when the triplet occupies x(8).
	! N_eq = 7 (+1 with the triplet) + 2*met_nelem.
	!
	! Physics of the coupling. The two blocks meet in the free electron
	! density: it is the H/He/molecular sum PLUS the metal charges (X+ once,
	! X++ twice). In the shielded, molecular base the metals are the dominant
	! electron donors, so this is not a small correction there -- it is the
	! reason the two networks cannot be solved apart from each other. Every
	! recombination and electron-impact term of BOTH blocks then sees the same
	! n_e. Metal <-> H/He charge exchange is added by cx_add_to_fvec with the
	! metal rows pointed at mbase; the atomic neutral H it reacts with is the
	! free H0 (the H nuclei bound into H2/H2+/H3+/HeH+ are excluded, exactly as
	! in the molecular rows), and the He I reservoir is the free neutral He.
	!
	! With met_nelem = 0 this system reduces exactly to System_HeH_mol, and
	! with the molecular fractions at zero its H/He/metal rows reduce to
	! System_HeH_metals; both limits are shared code, not re-derived here.
	!
	! Solved with hybrd1 (numerical Jacobian), matching the molecular system;
	! no analytic Jacobian is written for the merged residual.

	use global_parameters, only: thereis_HeITR
	use ion_cell_state,    only: ieq_cell
	use ion_residual_core, only: metal_fractions, metal_electron_sum,     &
	                             metal_rows
	use System_HeH_mol,    only: mol_heh_rows
	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,     &
	                             met_b0, met_b1, met_a1, met_a2, met_top
	use charge_exchange,   only: cx_add_to_fvec, he_h_cx_fvec

	implicit none

	contains

	! Row of the first metal unknown in the molecular layout: the molecular
	! network owns x(1..7), the He 2^3S metastable x(8) when tracked, and the
	! metals follow. Single definition, read by the residual here and by the
	! driver (unknown seeding, root validation, cx_metal_base).
	integer function metal_row_base()
	metal_row_base = 8
	if (thereis_HeITR) metal_row_base = 9
	end function metal_row_base

	subroutine ion_system_HeH_mol_metals(N_eq,x,fvec,iflag,params)

	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: params(60)

	real*8  :: g_hi,g_hei,g_heii,g_heiTR,g_h2   ! photoionization
	real*8  :: b_hi,b_hei,b_heii,b_heiTR        ! collisional ionization
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR     ! recombination
	real*8  :: A31,q13,q31a,q31b,Q31            ! triplet kinetics
	real*8  :: n_h,n_he,n_e,ntot
	real*8  :: n_hi,n_hii,n_h2,n_h2p,n_h3p,n_hehp
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	! Each element's metal densities (neutral/+/++); charge exchange later.
	real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)
	integer :: mbase

	! Unpack the cell state (same named fields as System_HeH_mol)
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
	ntot    = ieq_cell%ntot

	! Species densities (verbatim System_HeH_mol)
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

	! Metal ion densities (neutral/+/++) from fractions, canonical order.
	mbase = metal_row_base()
	call metal_fractions(x, mbase, met_nelem, met_ntot, nm0, nm1, nm2)

	! Electron density: H/He/molecular first (reproducing the molecular
	! system's sum), then the metal charges.
	n_e = n_hii + n_h2p + n_h3p + n_hehp + n_heii + 2.0d0*n_heiii
	call metal_electron_sum(n_e, met_nelem, nm1, nm2)

	! --- Molecular network rows (System_HeH_mol, metal-inclusive n_e) ---
	call mol_heh_rows(fvec, n_hi, n_hii, n_h2, n_h2p, n_h3p, n_hehp,   &
	                  n_heiSI, n_heiTR, n_heii, n_heiii, n_e, ntot,     &
	                  g_hi, g_hei, g_heii, g_heiTR, g_h2,               &
	                  a_hii, a_heii, a_heiii, a_heiTR,                  &
	                  b_hi, b_hei, b_heii, b_heiTR,                     &
	                  q13, q31a, q31b, Q31, A31)

	! --- Metal rows (System_HeH_metals, shifted to mbase..) ---
	call metal_rows(fvec, x, mbase, met_nelem, met_ntot, met_g0, met_g1, &
	                met_b0, met_b1, met_a1, met_a2, met_top,             &
	                nm0, nm1, nm2, n_e)

	! Charge exchange (Huang Table 4) on the H, He and metal rows. The
	! driver sets cx_metal_base = mbase before this solve. Rows 1 and 2 are
	! the H+ and He+ balances written production positive, the convention
	! cx_add_to_fvec assumes. Absent reactants contribute zero, preserving
	! the identity rows of absent elements.
	call cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                       &
	                    n_hi, n_hii, n_hei, n_heii, n_heiii)

	! He <-> H charge exchange (Huang Table 4 group B); rows 1 and 2 are
	! production positive here, so he_row_sign = +1. Group B is excluded
	! from cx_act, so it is applied only here (no double counting).
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,    &
	                  n_hi, n_hii, n_heiSI, n_heii, 1.0d0)

	return

	end subroutine ion_system_HeH_mol_metals

	! End of module
	end module System_HeH_mol_metals

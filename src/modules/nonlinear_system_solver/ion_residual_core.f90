	module ion_residual_core
	! Pure, stateless helpers for the standard H/He ionization residual rows
	! (with collisional ionization) and their analytic-Jacobian pieces. These
	! blocks are verbatim-shared by System_HeH and System_HeH_metals
	! (docs/refactor_plan_system_composition_parser.md, §5.2 Inc 1). Only
	! explicit-shape / assumed-size dummies are used here: assumed-shape (:)
	! dummies can change gfortran -O3 code generation and break the
	! byte-identical regression (review item 4.1). The arithmetic order is
	! kept exactly as in the original inline statements.

	implicit none

	contains

	! Standard three H/He residual rows. Writes fvec(1), fvec(2), fvec(3) in
	! this exact order. n_e is an INPUT, so the metals system can pass its own
	! metal-inclusive electron density.
	subroutine heh_rows(fvec, n_hi, n_hii, n_hei, n_heii, n_heiii, n_e,  &
	                    g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii,      &
	                    b_hi, b_hei, b_heii)
	real*8 :: fvec(*)
	real*8, intent(in) :: n_hi, n_hii, n_hei, n_heii, n_heiii, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii
	real*8, intent(in) :: b_hi, b_hei, b_heii

	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e
	fvec(2) = n_hei*g_hei + (n_hei*b_hei - a_heii*n_heii)*n_e
	fvec(3) = n_heii*g_heii + (n_heii*b_heii - a_heiii*n_heiii)*n_e
	end subroutine heh_rows

	! n_e-coefficient C_i of each standard H/He row (the factor multiplying
	! n_e in fvec_i). Sets C1, C2, C3 in this exact order.
	subroutine heh_crow(C1, C2, C3, n_hi, n_hii, n_hei, n_heii, n_heiii,  &
	                    a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
	real*8, intent(out) :: C1, C2, C3
	real*8, intent(in)  :: n_hi, n_hii, n_hei, n_heii, n_heiii
	real*8, intent(in)  :: a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii

	C1 = n_hi*b_hi - a_hii*n_hii
	C2 = n_hei*b_hei - a_heii*n_heii
	C3 = n_heii*b_heii - a_heiii*n_heiii
	end subroutine heh_crow

	! Local "direct" terms (photoionization + dC_i/dx_local * n_e) added to the
	! standard H/He Jacobian rows. Applies the five updates in this exact
	! order; call AFTER the rank-1 n_e coupling has populated fjac.
	subroutine heh_jac_local(nsz, fjac, n_h, n_he, n_e, g_hi, g_hei, g_heii, &
	                         a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
	integer, intent(in) :: nsz
	real*8 :: fjac(nsz,nsz)
	real*8, intent(in) :: n_h, n_he, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii
	real*8, intent(in) :: b_hi, b_hei, b_heii

	fjac(1,1) = fjac(1,1) - n_h*g_hi + (-n_h*b_hi - a_hii*n_h)*n_e
	fjac(2,2) = fjac(2,2) - n_he*g_hei + (-n_he*b_hei - a_heii*n_he)*n_e
	fjac(2,3) = fjac(2,3) - n_he*g_hei + (-n_he*b_hei)*n_e
	fjac(3,2) = fjac(3,2) + n_he*g_heii + (n_he*b_heii)*n_e
	fjac(3,3) = fjac(3,3) + (-a_heiii*n_he)*n_e
	end subroutine heh_jac_local

	! End of module
	end module ion_residual_core

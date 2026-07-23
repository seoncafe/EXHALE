	module ion_residual_core
	! Pure, stateless helpers for the standard H/He ionization residual rows
	! (with collisional ionization) and their analytic-Jacobian pieces, plus
	! the TR-form H/He/triplet rows (Oklopcic form; electron-impact ionization
	! of H0/He(1^1S)/He+ and of the He 2^3S metastable is included).
	! The standard blocks are verbatim-shared by System_HeH and
	! System_HeH_metals (§5.2 Inc 1); the TR-form rows and the He 2^3S triplet
	! row are shared by System_HeH_TR, System_HeH_TR_metals and System_HeH_mol
	! (§5.2 Inc 3, docs/refactor_plan_system_composition_parser.md). Only
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

	! He 2^3S triplet balance row (Oklopcic form + electron-impact ionization
	! of the metastable). Sets ftr to the single triplet residual expression;
	! verbatim-shared as row 4 of the TR / TR_metals systems and row 8 of the
	! molecular system. b_heiTR is the He(2^3S) collisional-ionization
	! coefficient (ci_HeI23S, threshold 4.8 eV): He(2^3S)+e- -> He+ + 2e-
	! removes the triplet, so it enters as a destruction term -n_e*n_heiTR*b_heiTR.
	subroutine tr_triplet_row(ftr, n_hi, n_heiSI, n_heiTR, n_heii, n_e,   &
	                          g_heiTR, a_heiTR, q13, q31a, q31b, Q31, A31, &
	                          b_heiTR)
	real*8, intent(out) :: ftr
	real*8, intent(in)  :: n_hi, n_heiSI, n_heiTR, n_heii, n_e
	real*8, intent(in)  :: g_heiTR, a_heiTR, q13, q31a, q31b, Q31, A31
	real*8, intent(in)  :: b_heiTR

	ftr = - n_heiTR*g_heiTR                                         &
	      + n_e*( n_heii*a_heiTR                                    &
	            + n_heiSI*q13                                       &
	            - n_heiTR*(q31a + q31b)                             &
	            - n_heiTR*b_heiTR)                                  &
	      - n_heiTR*(A31 + n_hi*Q31)
	end subroutine tr_triplet_row

	! TR-form four H/He/triplet residual rows (Oklopcic form, no collisional
	! ionization). Writes fvec(1), fvec(2), fvec(3) in this exact order, then
	! obtains fvec(4) from tr_triplet_row so the triplet expression lives in
	! one place. Verbatim-shared by System_HeH_TR and System_HeH_TR_metals;
	! n_e is an INPUT, so the metals system can pass its metal-inclusive
	! electron density.
	subroutine heh_tr_rows(fvec, n_hi, n_hii, n_heiSI, n_heiTR, n_heii,   &
	                       n_heiii, n_e, g_hi, g_hei, g_heii, g_heiTR,     &
	                       a_hii, a_heii, a_heiii, a_heiTR,                &
	                       b_hi, b_hei, b_heii, b_heiTR,                   &
	                       q13, q31a, q31b, Q31, A31)
	real*8 :: fvec(*)
	real*8, intent(in) :: n_hi, n_hii, n_heiSI, n_heiTR, n_heii, n_heiii, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, g_heiTR
	real*8, intent(in) :: a_hii, a_heii, a_heiii, a_heiTR
	real*8, intent(in) :: b_hi, b_hei, b_heii, b_heiTR
	real*8, intent(in) :: q13, q31a, q31b, Q31, A31

	! Penning ionization source He(2^3S)+H0 -> He(1^1S) + H+ + e- (rate Q31,
	! Taylor 2025): the same event that removes the triplet in tr_triplet_row
	! ionizes H0, so it enters here as an H+ production term. All products are
	! counted as H+; the associative branch (-> HeH+ + e-, ~10%, GM25) is not
	! resolved here. n_hi*b_hi*n_e is the electron-impact ionization of H0
	! (Voronov b_hi), restored to match the standard heh_rows.
	fvec(1) = n_hi*g_hi + n_heiTR*n_hi*Q31 + n_hi*b_hi*n_e - a_hii*n_hii*n_e

	! New equation for hei - sum of the two equations of Oklopcic. The last
	! two terms are electron-impact ionization of ground-state (b_hei) and
	! metastable (b_heiTR) He I; both remove He I and produce He+, so they sit
	! on the loss side of this summed He I balance.
	fvec(2) =   n_heii*(a_heiTR + a_heii)*n_e                            &
	          - n_heiSI*g_hei                                           &
	          - n_heiTR*g_heiTR                                         &
	          - n_heiSI*b_hei*n_e                                       &
	          - n_heiTR*b_heiTR*n_e

	! n_heii*b_heii*n_e: electron-impact ionization of He+ into He++.
	fvec(3) = n_heii*g_heii + n_heii*b_heii*n_e - a_heiii*n_heiii*n_e

	call tr_triplet_row(fvec(4), n_hi, n_heiSI, n_heiTR, n_heii, n_e,     &
	                    g_heiTR, a_heiTR, q13, q31a, q31b, Q31, A31,      &
	                    b_heiTR)
	end subroutine heh_tr_rows

	! Metal ion densities (neutral/+/++) from fractions, in canonical element
	! order. `base` is the row of the first metal unknown (4 in the metals
	! system, 5 with the He 2^3S triplet), so ix = base + 2*(e-1) locates the
	! X+/X++ fractions of element e. Sets nm0/nm1/nm2 in this exact order.
	subroutine metal_fractions(x, base, nelem, mtot, nm0, nm1, nm2)
	integer, intent(in) :: base, nelem
	real*8, intent(in)  :: x(*)
	real*8, intent(in)  :: mtot(nelem)
	real*8, intent(out) :: nm0(nelem), nm1(nelem), nm2(nelem)
	integer :: e, ix
	real*8  :: n_X

	do e = 1,nelem
		ix     = base + 2*(e-1)
		n_X    = mtot(e)
		nm1(e) = x(ix)*n_X
		nm2(e) = x(ix+1)*n_X
		nm0(e) = (1.0 - x(ix) - x(ix+1))*n_X
	enddo
	end subroutine metal_fractions

	! Add the metal electron contribution to n_e (X+ counts once, X++ twice),
	! in canonical element order. n_e accumulates in place, reproducing the
	! original loop's add order exactly.
	subroutine metal_electron_sum(n_e, nelem, nm1, nm2)
	integer, intent(in)   :: nelem
	real*8, intent(inout) :: n_e
	real*8, intent(in)    :: nm1(nelem), nm2(nelem)
	integer :: e

	do e = 1,nelem
		n_e = n_e + nm1(e) + 2.0*nm2(e)
	enddo
	end subroutine metal_electron_sum

	! Metal ionization balance rows, one element at a time (force-zero if the
	! element is absent; otherwise the normal balance). `base` locates the
	! first metal row (4 or 5); ix = base + 2*(e-1). Three-stage elements
	! (mtop >= 2) solve both X0<->X+ and X+<->X++; two-stage elements solve
	! only X0<->X+ and pin the unused upper unknown. Expression order is
	! verbatim from the System_HeH_metals residual.
	subroutine metal_rows(fvec, x, base, nelem, mtot, mg0, mg1,        &
	                      mb0, mb1, ma1, ma2, mtop, nm0, nm1, nm2, n_e)
	integer, intent(in) :: base, nelem
	real*8 :: fvec(*)
	real*8, intent(in)  :: x(*)
	real*8, intent(in)  :: mtot(nelem), mg0(nelem), mg1(nelem)
	real*8, intent(in)  :: mb0(nelem), mb1(nelem), ma1(nelem), ma2(nelem)
	integer, intent(in) :: mtop(nelem)
	real*8, intent(in)  :: nm0(nelem), nm1(nelem), nm2(nelem)
	real*8, intent(in)  :: n_e
	integer :: e, ix

	do e = 1,nelem
		ix = base + 2*(e-1)
		if (mtot(e) .le. 1.0d-30) then
			fvec(ix)   = x(ix)
			fvec(ix+1) = x(ix+1)
		else
			! X0 <-> X+
			fvec(ix)   = nm0(e)*mg0(e)                            &
			           + (nm0(e)*mb0(e) - ma1(e)*nm1(e))*n_e
			if (mtop(e) .ge. 2) then
				! X+ <-> X++
				fvec(ix+1) = nm1(e)*mg1(e)                            &
				           + (nm1(e)*mb1(e) - ma2(e)*nm2(e))*n_e
			else
				! Two-stage element: no X++, pin the unused unknown.
				fvec(ix+1) = x(ix+1)
			endif
		endif
	enddo
	end subroutine metal_rows

	! End of module
	end module ion_residual_core

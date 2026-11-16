!  THE ATOMIC H / He / He 2^3S ROWS, ONE SOURCE TEXT FOR TWO ARITHMETICS:
!  included by heh_tr_rows_real (real*8 densities) and heh_tr_rows_dual
!  (dual densities, value and first derivatives) of ion_residual_core.

	! Penning ionization source He(2^3S)+H0 -> He(1^1S) + H+ + e-: the same
	! collision that removes the triplet in tr_triplet_row ionizes H0, so it
	! enters here as an H+ production term. Q31 is the TOTAL ionization rate,
	! of which only the Penning branch f_penning_HeI23S leaves a proton behind.
	! The remaining 10% is associative ionization He(2^3S)+H0 -> HeH+ + e-;
	! this system carries no HeH+, and in an atomic gas HeH+ dissociatively
	! recombines back to He + H (1.3e-8 cm^3 s^-1 at 2000 K, faster than any
	! competing HeH+ reaction there), returning the H atom and consuming the
	! electron, so that branch is a metastable sink but not a lasting proton
	! or electron source. Where H2 is abundant that cycle is broken by
	! HeH+ + H2 -> H3+ + He, and System_HeH_mol therefore does carry the
	! branch explicitly into its HeH+ row.
	! n_hi*b_hi*n_e is the electron-impact ionization of H0 (Voronov b_hi),
	! restored to match the standard heh_rows.
	fvec(1) = n_hi*g_hi + f_penning_HeI23S*n_heiTR*n_hi*Q31            &
	        + n_hi*b_hi*n_e - a_hii*n_hii*n_e

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
	                    g_heiTR, a_heiTR, q13, q31g, q31a, q31b, Q31,     &
	                    A31, b_heiTR)

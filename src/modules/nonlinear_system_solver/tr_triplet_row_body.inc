!  THE He 2^3S LEVEL BALANCE, ONE SOURCE TEXT FOR TWO ARITHMETICS: included
!  by tr_triplet_row_real (real*8 densities) and tr_triplet_row_dual (dual
!  densities, value and first derivatives) of ion_residual_core.  The
!  physics of the row is stated above tr_triplet_row_real.

	ftr = - n_heiTR*g_heiTR                                         &
	      + n_e*( n_heii*a_heiTR                                    &
	            + n_heiSI*q13                                       &
	            - n_heiTR*(q31g + q31a + q31b)                      &
	            - n_heiTR*b_heiTR)                                  &
	      - n_heiTR*(A31 + n_hi*Q31)

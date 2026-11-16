!  THE TWO OXYGEN-CARRIER ROWS AND THE OXYGEN CYCLE'S EXCHANGE WITH THE H2 ROW,
!  ONE SOURCE TEXT FOR TWO ARITHMETICS: included by oxygen_carrier_rows_real
!  (real*8 densities) and oxygen_carrier_rows_dual (dual densities, value
!  and first derivatives) of System_HeH_mol (reactions stated above
!  oxygen_carrier_rows_real).

	! (iox) OH balance
	prod_oh     = ok1r*n_h2o*n_hi + ok2*n_o0*n_h2                        &
	            + (pj3 + pj4)*n_h2o
	loss_oh_o1  = ok1*n_oh*n_h2
	loss_oh_o2r = ok2r*n_oh*n_hi
	loss_oh_j7  = pj7*n_oh
	fvec(iox)   = prod_oh - loss_oh_o1 - loss_oh_o2r - loss_oh_j7

	! (iox+1) H2O balance
	prod_h2o     = ok1*n_oh*n_h2
	loss_h2o_o1r = ok1r*n_h2o*n_hi
	loss_h2o_j   = (pj3 + pj4 + pj5)*n_h2o
	fvec(iox+1)  = prod_h2o - loss_h2o_o1r - loss_h2o_j

	! H2 exchange with the oxygen cycle, added to the molecular H2 balance.
	! O1 and O2 consume an H2 each; their reverses give one back. O3, O5 and
	! O7 do not touch H2, and the O4/O6 pair cancels (see above).
	loss_h2_oxy = (ok1*n_oh + ok2*n_o0)*n_h2
	prod_h2_oxy = (ok1r*n_h2o + ok2r*n_oh)*n_hi
	fvec(4) = fvec(4)                                                    &
	        - loss_h2_oxy                                                &
	        + prod_h2_oxy

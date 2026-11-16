!  THE R12/R15 THIRD-BODY DENSITY, ONE SOURCE TEXT FOR TWO ARITHMETICS:
!  included by h2_third_body_density_real and h2_third_body_density_dual of
!  System_HeH_mol (physics stated above h2_third_body_density_real).
	if (mk_third_body_resolved) then
		n_eff = h2_association_collider_density(ieq_cell%T_K, n_h2,   &
		                                        n_hi, n_hei)
	else
		n_eff = ntot
	endif

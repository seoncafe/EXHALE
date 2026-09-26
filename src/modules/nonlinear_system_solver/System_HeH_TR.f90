	module System_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use ion_cell_state, only: ieq_cell
	use ion_residual_core, only: heh_tr_rows,                            &
	                            impose_transported_ionization_fractions
	use charge_exchange, only: he_h_cx_fvec

	implicit none
	
	contains
	
	subroutine ion_system_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR		! Photoionization rates
	real*8  :: b_hi,b_hei,b_heii,b_heiTR		! Collisional ionization rates
	real*8  :: A31,q13,q31g,q31a,q31b,Q31
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR	! Recombination rates
   real*8  :: params(40)
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_heii,n_heiii,n_heiTR,n_heiSI
	
	! Coefficients of the system

 	g_hi    = ieq_cell%P_HI        ! = P_HI
 	g_hei   = ieq_cell%P_HeI       ! = P_HeI
 	g_heii  = ieq_cell%P_HeII      ! = P_HeII
 	a_hii   = ieq_cell%rchiiB      ! = rchiiB
 	a_heii  = ieq_cell%rcheiiB     ! = rcheiiB
 	a_heiii = ieq_cell%rcheiiiB    ! = rcheiiiB
 	n_h     = ieq_cell%nh          ! = nh
 	n_he    = ieq_cell%nhe         ! = nhe
   b_hi    = ieq_cell%a_ion_HI    ! = a_ion_HI
 	b_hei   = ieq_cell%a_ion_HeI   ! = a_ion_HeI
 	b_heii  = ieq_cell%a_ion_HeII  ! = a_ion_HeII

 	! Triplet parameters
 	a_heiTR = ieq_cell%rcheiTR  ! = rcheiTR
 	A31     = ieq_cell%A31      ! = A31
 	g_heiTR = ieq_cell%P_HeITR  ! = P_HeITR
 	q13     = ieq_cell%q13      ! = q13
 	q31g    = ieq_cell%q31g     ! = q31g (reverse of q13)
 	q31a    = ieq_cell%q31a     ! = q31a
 	q31b    = ieq_cell%q31b     ! = q31b
 	Q31     = ieq_cell%Q31      ! = Q31
 	b_heiTR = ieq_cell%a_ion_HeITR  ! = a_ion_HeITR (He 2^3S collisional ioniz.)


 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_heii  = x(2)*n_he
 	n_heiii = x(3)*n_he
 	n_heiSI = (1.0 - x(2) - x(3) - x(4))*n_he
 	n_heiTR = x(4)*n_he
 	
 	
 	! Electron density
   n_e = n_hii + n_heii + 2.0*n_heiii

    ! System of equations (verbatim TR-form rows now live in ion_residual_core)
	call heh_tr_rows(fvec, n_hi, n_hii, n_heiSI, n_heiTR, n_heii, n_heiii,  &
	                 n_e, g_hi, g_hei, g_heii, g_heiTR,                      &
	                 a_hii, a_heii, a_heiii, a_heiTR,                        &
	                 b_hi, b_hei, b_heii, b_heiTR,                           &
	                 q13, q31g, q31a, q31b, Q31, A31)

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
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,  &
	                  ieq_cell%kcx_Hepp_H0,                            &
	                  n_hi, n_hii, n_heiSI, n_heii, n_heiii, -1.0d0,    &
	                  .false.)

	! The transported ionization fractions, where the flow carries them
	! and not this cell's local balance (ion_residual_core).
	call impose_transported_ionization_fractions(ieq_cell, x, fvec)

	return

	! End of subroutine
	end subroutine ion_system_HeH_TR
	
	! End of module
	end module System_HeH_TR

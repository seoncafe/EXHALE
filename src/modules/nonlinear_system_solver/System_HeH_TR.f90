	module System_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use ion_cell_state, only: ieq_cell
	use ion_residual_core, only: heh_tr_rows

	implicit none
	
	contains
	
	subroutine ion_system_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR		! Photoionization rates
	real*8  :: b_hi,b_hei,b_heii			! Collisional ionization rates
	real*8  :: A31,q13,q31a,q31b,Q31
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR	! Recombination rates
   real*8  :: params(40)
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	
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
 	q31a    = ieq_cell%q31a     ! = q31a
 	q31b    = ieq_cell%q31b     ! = q31b
 	Q31     = ieq_cell%Q31      ! = Q31
 	
 	
 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_hei   = (1.0 - x(2) - x(3))*n_he 
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
	                 q13, q31a, q31b, Q31, A31)

	return
	
	! End of subroutine
	end subroutine ion_system_HeH_TR
	
	! End of module
	end module System_HeH_TR

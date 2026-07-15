	module System_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use params_idx, only: IPE_PHI, IPE_PHEI, IPE_PHEII, IPE_AHII, IPE_AHEII,   &
	                      IPE_AHEIII, IPE_NH, IPE_NHE, IPE_BHI, IPE_BHEI,       &
	                      IPE_BHEII, IPE_ATR, IPE_A31, IPE_PTR, IPE_Q13,        &
	                      IPE_Q31A, IPE_Q31B, IPE_Q31

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

 	g_hi    = params(IPE_PHI)    ! = P_HI 
 	g_hei   = params(IPE_PHEI)    ! = P_HeI
 	g_heii  = params(IPE_PHEII)    ! = P_HeII 
 	a_hii   = params(IPE_AHII)    ! = rchiiB 
 	a_heii  = params(IPE_AHEII)    ! = rcheiiB 
 	a_heiii = params(IPE_AHEIII)    ! = rcheiiiB 
 	n_h     = params(IPE_NH)    ! = nh 
 	n_he    = params(IPE_NHE)    ! = nhe 
   b_hi    = params(IPE_BHI)    ! = a_ion_HI 
 	b_hei   = params(IPE_BHEI)   ! = a_ion_HeI 
 	b_heii  = params(IPE_BHEII)   ! = a_ion_HeII 
 	
 	! Triplet parameters
 	a_heiTR = params(IPE_ATR)   ! = rcheiTR
 	A31     = params(IPE_A31)   ! = A31
 	g_heiTR = params(IPE_PTR)   ! = P_HeITR
 	q13     = params(IPE_Q13)   ! = q13
 	q31a    = params(IPE_Q31A)   ! = q31a
 	q31b    = params(IPE_Q31B)   ! = q31b
 	Q31     = params(IPE_Q31)   ! = Q31
 	
 	
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

    ! System of equations      
  	fvec(1) = n_hi*g_hi - a_hii*n_hii*n_e  
  	
  	! New equation for hei - sum of the two equations of Oklopcic
  	fvec(2) =  n_heii*(a_heiTR + a_heii)*n_e	&
  		     - n_heiSI*g_hei					      &
  		     - n_heiTR*g_heiTR
  	
  	fvec(3) = n_heii*g_heii - a_heiii*n_heiii*n_e
  	
  	fvec(4) = - n_heiTR*g_heiTR 		  	      &
  		      + n_e*( n_heii*a_heiTR   	      &
			        + n_heiSI*q13   		      &
		    	    - n_heiTR*(q31a + q31b))	   &
		       - n_heiTR*(A31 + n_hi*Q31)
	
	return 
	
	! End of subroutine
	end subroutine ion_system_HeH_TR
	
	! End of module
	end module System_HeH_TR

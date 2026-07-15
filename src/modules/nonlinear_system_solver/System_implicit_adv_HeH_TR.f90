	module System_implicit_adv_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use params_idx, only: IPA_C1, IPA_XHI, IPA_XHEI, IPA_XHEIII, IPA_NH,        &
	                      IPA_PHI, IPA_PHEI, IPA_PHEII, IPA_AHII, IPA_AHEII,     &
	                      IPA_AHEIII, IPA_BHI, IPA_BHEI, IPA_BHEII, IPAT_ATR,    &
	                      IPAT_A31, IPAT_PTR, IPAT_Q13, IPAT_Q31A, IPAT_Q31B,    &
	                      IPAT_Q31, IPAT_XTR, IPAT_HEH

	implicit none
	
	contains
	
	subroutine adv_implicit_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: xhi_old,xhei_old,xheiii_old,xheiTR_old
	real*8  :: ghi,ghei,gheii,gheiTR
	real*8  :: xhi,xhii
	real*8  :: xhei,xheii,xheiii
	real*8  :: xheiTR, xheiS
	real*8  :: xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii,aheiTR
	real*8  :: ionhi,ionhei,ionheii
	real*8  :: heh_loc
   real*8  :: params(25)
   real*8  :: A31,q13,q31a,q31b,Q31
	
	! Coefficients of the system
 	c1         = params(IPA_C1)    ! = dr/v
 	xhi_old    = params(IPA_XHI)    ! = nhi/nhe 
 	xhei_old   = params(IPA_XHEI)    ! = nheii/nhe 
 	xheiii_old = params(IPA_XHEIII)    ! = nheiii/nhe
 	n_h        = params(IPA_NH)    ! = nh
 	ghi        = params(IPA_PHI)    ! = P_HI 
 	ghei       = params(IPA_PHEI)    ! = P_HeI 
 	gheii      = params(IPA_PHEII)    ! = P_HeII 
   ahii       = params(IPA_AHII)    ! = rchiiB  
 	aheii      = params(IPA_AHEII)   ! = rcheiiB 
 	aheiii     = params(IPA_AHEIII)   ! = rcheiiiB  
	ionhi	     = params(IPA_BHI)   ! = a_ion_HI
	ionhei     = params(IPA_BHEI)   ! = a_ion_HeI
	ionheii    = params(IPA_BHEII)   ! = a_ion_HEII
	aheiTR     = params(IPAT_ATR)	  ! = rcheiTR
	A31	     = params(IPAT_A31)   ! = A31
	gheiTR     = params(IPAT_PTR)   ! = P_HeITR
	q13        = params(IPAT_Q13)   ! = q13
	q31a	     = params(IPAT_Q31A)   ! = q31a
	q31b	     = params(IPAT_Q31B)   ! = q31b
	Q31 	     = params(IPAT_Q31)   ! = Q31
	xheiTR_old = params(IPAT_XTR)   ! = nheiTR/nh
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (legacy, byte-identical); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead.
	heh_loc    = params(IPAT_HEH)   ! = He/H (local when he_diffusion)

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xhei   = x(2)
   xheii  = 1.0 - x(2) - x(3)
	xheiii = x(3)
	xheiS  = x(2) - x(4)
	xheiTR = x(4)

 	! Electron density
   xe = xhii + heh_loc*(xheii + 2.0*xheiii)
      
    ! System of equations      
  	fvec(1) =  xhi_old - xhi + c1*(  		&
  	           - (ghi + ionhi*xe*n_h)*xhi 	&
  	           +  ahii*xhii*xe*n_h) 	   
  	
  	fvec(2) =  xhei_old - xhei + c1*(  	       &
  		       xheii*(aheiTR + aheii)*xe*n_h	 &
  		     - xheiS*ghei - xheiTR*gheiTR) 
  	 	      	 	    
  	fvec(3) =  xheiii_old - xheiii + c1*(	    &
  	           (gheii + ionheii*xe*n_h)*xheii  & 
  	          - aheiii*xheiii*xe*n_h) 
  	           
	fvec(4) =    xheiTR_old - xheiTR + c1*(   &
		     - gheiTR*xheiTR				         &
		     + (xheii*aheiTR + xheiS*q13    	&
		     -  xheiTR*(q31a + q31b))*xe*n_h	&
		     - xheiTR*(A31 + xhi*Q31*n_h))
	
	! End of subroutine
	end subroutine adv_implicit_HeH_TR
	
	! End of module
	end module System_implicit_adv_HeH_TR

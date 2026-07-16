	module System_implicit_adv_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use ion_cell_state, only: adv_cell

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
 	c1         = adv_cell%c1    ! = dr/v
 	xhi_old    = adv_cell%xhi_old    ! = nhi/nhe 
 	xhei_old   = adv_cell%xhei_old    ! = nheii/nhe 
 	xheiii_old = adv_cell%xheiii_old    ! = nheiii/nhe
 	n_h        = adv_cell%nh    ! = nh
 	ghi        = adv_cell%P_HI    ! = P_HI 
 	ghei       = adv_cell%P_HeI    ! = P_HeI 
 	gheii      = adv_cell%P_HeII    ! = P_HeII 
   ahii       = adv_cell%rchiiB    ! = rchiiB  
 	aheii      = adv_cell%rcheiiB   ! = rcheiiB 
 	aheiii     = adv_cell%rcheiiiB   ! = rcheiiiB  
	ionhi	     = adv_cell%a_ion_HI   ! = a_ion_HI
	ionhei     = adv_cell%a_ion_HeI   ! = a_ion_HeI
	ionheii    = adv_cell%a_ion_HeII   ! = a_ion_HEII
	aheiTR     = adv_cell%rcheiTR	  ! = rcheiTR
	A31	     = adv_cell%A31   ! = A31
	gheiTR     = adv_cell%P_HeITR   ! = P_HeITR
	q13        = adv_cell%q13   ! = q13
	q31a	     = adv_cell%q31a   ! = q31a
	q31b	     = adv_cell%q31b   ! = q31b
	Q31 	     = adv_cell%Q31   ! = Q31
	xheiTR_old = adv_cell%xheiTR_old   ! = nheiTR/nh
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (legacy, byte-identical); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead.
	heh_loc    = adv_cell%heh_loc   ! = He/H (local when he_diffusion)

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

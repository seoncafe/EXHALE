	module System_implicit_adv_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use params_idx, only: IPA_C1, IPA_XHI, IPA_XHEI, IPA_XHEIII, IPA_NH,        &
	                      IPA_PHI, IPA_PHEI, IPA_PHEII, IPA_AHII, IPA_AHEII,     &
	                      IPA_AHEIII, IPA_BHI, IPA_BHEI, IPA_BHEII, IPA_HEH

	implicit none
	
	contains
	
	subroutine adv_implicit_HeH(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: xhi_old,xhei_old,xheiii_old
	real*8  :: ghi,ghei,gheii
	real*8  :: xhi,xhii
	real*8  :: xhei,xheii,xheiii
	real*8  :: xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii
	real*8  :: ionhi,ionhei,ionheii
	real*8  :: heh_loc
   real*8  :: params(25)

	! Coefficients of the system

 	c1         = params(IPA_C1)    ! = dr/v
 	xhi_old    = params(IPA_XHI)    ! = nhi/nh
 	xhei_old   = params(IPA_XHEI)    ! = nhei/nhe
 	xheiii_old = params(IPA_XHEIII)    ! = nheiii/nhe
 	n_h        = params(IPA_NH)    ! = nh
 	ghi        = params(IPA_PHI)    ! = P_HI
 	ghei       = params(IPA_PHEI)    ! = P_HeI = nh
 	gheii      = params(IPA_PHEII)    ! = P_HeII = nhe
   ahii       = params(IPA_AHII)    ! = rchiiB
 	aheii      = params(IPA_AHEII)   ! = rcheiiB
 	aheiii     = params(IPA_AHEIII)   ! = rcheiiiB
	ionhi	   = params(IPA_BHI)   ! = a_ion_HI
	ionhei     = params(IPA_BHEI)   ! = a_ion_HeI
	ionheii    = params(IPA_BHEII)   ! = a_ion_HEII
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (legacy, byte-identical); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead (the global HeH would misstate
	! n_e by the local separation factor).
	heh_loc    = params(IPA_HEH)   ! = He/H (local when he_diffusion)

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xhei   = x(2)
   xheii  = 1.0 - x(2) - x(3)
	xheiii = x(3)

 	! Electron density
   xe = xhii + heh_loc*(xheii + 2.0*xheiii)
      
   ! System of equations      
  	fvec(1) =  xhi_old - xhi   					        		&
  	        + c1*(-(ghi+ionhi)*xhi    + ahii*xhii*xe*n_h) 	   
  		     		    
  	fvec(2) =  xhei_old - xhei 						  		    &
  		    + c1*(-(ghei+ionhei)*xhei  + aheii*xheii*xe*n_h) 
  	 	      	 	    
  	fvec(3) =  xheiii_old - xheiii					  			&
  	        + c1*((gheii+ionheii)*xheii - aheiii*xheiii*xe*n_h) 
  		    
	! End of subroutine
	end subroutine adv_implicit_HeH
	
	! End of module
	end module System_implicit_adv_HeH

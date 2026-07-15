	module System_implicit_adv_H
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use params_idx, only: IPAH_C1, IPAH_XHI, IPAH_NH, IPAH_PHI, IPAH_AHII, IPAH_BHI

	implicit none
	
	contains
	
	subroutine adv_implicit_H(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: xhi_old
	real*8  :: ghi
	real*8  :: xhi,xhii,xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii
	real*8  :: ionhi
    real*8  :: params(25)
	
	
	! Coefficients of the system
 	c1      = params(IPAH_C1)    ! = dr/v
 	xhi_old = params(IPAH_XHI)    ! = nhi/nh
 	n_h     = params(IPAH_NH)    ! = nh
 	ghi     = params(IPAH_PHI)    ! = P_HI  
   ahii    = params(IPAH_AHII)    ! = rchiiB  
	ionhi   = params(IPAH_BHI)    ! = a_ion_HI
	
	! Substitutions
	xhi = x(1)
	xhii = 1.0 - x(1)

 	! Electron density
   xe = xhii
      
      ! System of equations      
  	fvec(1) =  xhi_old - x(1)					&
  		    + c1*(-(ghi+ionhi*xe*n_h)*xhi + ahii*xhii*xe*n_h)
  		   
	! End of subroutine
	end subroutine adv_implicit_H
	
	! End of module
	end module System_implicit_adv_H

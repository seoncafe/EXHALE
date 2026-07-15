	module System_H
	! Ionization equilibrium system with only H
	
	use global_parameters
	use params_idx, only: IPH_PHI, IPH_AHII, IPH_NH, IPH_BHI

	implicit none
	
	contains
	
	subroutine ion_system_H(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: g_hi                   ! Photoionization rates
	real*8  :: b_hi                   ! Collisional ionization rates
	real*8  :: a_hii                  ! Recombination rates
   real*8  :: params(40)
	real*8  :: n_h,n_e
	real*8  :: n_hi,n_hii
	
	! Coefficients of the system

 	g_hi    = params(IPH_PHI)    ! = P_HI 
 	a_hii   = params(IPH_AHII)    ! = rchiiB 
 	n_h     = params(IPH_NH)    ! = nh 
   b_hi    = params(IPH_BHI)    ! = a_ion_HI 
 	
 	! Species densities
 	n_hi  = (1.0-x(1))*n_h
 	n_hii = x(1)*n_h
 	
 	! Electron density
   n_e = n_hii
      
    ! System of equations      
  	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e

	! End of subroutine
	end subroutine ion_system_H

	! Analytic Jacobian of ion_system_H (Task 2). f1 = n_hi*g_hi + C1*n_e
	! with n_hi=(1-x1)nh, n_hii=x1*nh, n_e=n_hii, C1=n_hi*b_hi - a_hii*n_hii.
	subroutine jac_system_H(Neq,x,fjac,params)
	integer :: Neq
	real*8  :: x(Neq), fjac(Neq,Neq), params(40)
	real*8  :: g_hi, a_hii, n_h, b_hi
	real*8  :: n_hi, n_hii, n_e, C1
	g_hi  = params(IPH_PHI)
	a_hii = params(IPH_AHII)
	n_h   = params(IPH_NH)
	b_hi  = params(IPH_BHI)
	n_hi  = (1.0 - x(1))*n_h
	n_hii = x(1)*n_h
	n_e   = n_hii
	C1    = n_hi*b_hi - a_hii*n_hii
	! direct (-nh*g_hi + dC1/dx1 * n_e) + n_e coupling (C1*dne, dne=nh)
	fjac(1,1) = -n_h*g_hi + (-n_h*b_hi - a_hii*n_h)*n_e + C1*n_h
	end subroutine jac_system_H

	! End of module
	end module System_H

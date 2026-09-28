	module System_implicit_adv_H
	! Advection-corrected ionization system for hydrogen alone: one step of
	! the backward differentiation formula for v dx/dr = R(x) in the neutral
	! fraction x = n_HI/n_H, taken by post_process_adv (the step is stated
	! at the adv_rates type of ion_cell_state).
	
	use global_parameters
	use ion_cell_state, only: adv_cell

	implicit none
	
	contains
	
	subroutine adv_implicit_H(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: xhi_hist
	real*8  :: ghi
	real*8  :: xhi,xhii,xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii
	real*8  :: ionhi
    real*8  :: params(25)
	
	
	! Coefficients of the system
 	c1      = adv_cell%c1    ! = g*h_j/v_j, the rate weight of the step
 	xhi_hist = adv_cell%xhi_hist    ! history of n_HI/n_H
 	n_h     = adv_cell%nh    ! = nh
 	ghi     = adv_cell%P_HI    ! = P_HI  
   ahii    = adv_cell%rchiiB    ! = rchiiB  
	ionhi   = adv_cell%a_ion_HI    ! = a_ion_HI
	
	! Substitutions
	xhi = x(1)
	xhii = 1.0 - x(1)

 	! Electron density, per H nucleus. The metal electrons (adv_cell%xe_metal,
 	! the same X+/X++ sum the equilibrium residual counts) are included: they
 	! dominate the electron budget of the shielded base, where the H ionized
 	! fraction is vanishingly small.
   xe = xhii + adv_cell%xe_metal
      
      ! System of equations      
  	fvec(1) =  xhi_hist - x(1)					&
  		    + c1*(-(ghi+ionhi*xe*n_h)*xhi + ahii*xhii*xe*n_h)
  		   
	! End of subroutine
	end subroutine adv_implicit_H
	
	! End of module
	end module System_implicit_adv_H

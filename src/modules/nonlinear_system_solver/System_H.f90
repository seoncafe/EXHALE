	module System_H
	! Ionization equilibrium system with only H
	
	use global_parameters
	use ion_cell_state, only: ieq_cell
	use ion_residual_core, only: impose_transported_ionization_fractions, &
	                             imposed_ionization_jacobian, &
	                             ionization_row_gross

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

 	g_hi    = ieq_cell%P_HI      ! = P_HI
 	a_hii   = ieq_cell%rchiiB    ! = rchiiB
 	n_h     = ieq_cell%nh        ! = nh
   b_hi    = ieq_cell%a_ion_HI  ! = a_ion_HI
 	
 	! Species densities
 	n_hi  = (1.0-x(1))*n_h
 	n_hii = x(1)*n_h
 	
 	! Electron density
   n_e = n_hii
      
    ! System of equations      
  	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e
	! Its gross rate, term by term (ion_residual_core, ionization_row_gross).
	ionization_row_gross(1) = abs(n_hi*g_hi) + abs(n_hi*b_hi*n_e)        &
	                        + abs(a_hii*n_hii*n_e)

	! The transported ionization fractions, where the flow carries them
	! and not this cell's local balance (ion_residual_core).
	call impose_transported_ionization_fractions(ieq_cell, x, fvec)

	end subroutine ion_system_H

	! Analytic Jacobian of ion_system_H (Task 2). f1 = n_hi*g_hi + C1*n_e
	! with n_hi=(1-x1)nh, n_hii=x1*nh, n_e=n_hii, C1=n_hi*b_hi - a_hii*n_hii.
	subroutine jac_system_H(Neq,x,fjac,params)
	integer :: Neq
	real*8  :: x(Neq), fjac(Neq,Neq), params(40)
	real*8  :: g_hi, a_hii, n_h, b_hi
	real*8  :: n_hi, n_hii, n_e, C1
	g_hi  = ieq_cell%P_HI
	a_hii = ieq_cell%rchiiB
	n_h   = ieq_cell%nh
	b_hi  = ieq_cell%a_ion_HI
	n_hi  = (1.0 - x(1))*n_h
	n_hii = x(1)*n_h
	n_e   = n_hii
	C1    = n_hi*b_hi - a_hii*n_hii
	! direct (-nh*g_hi + dC1/dx1 * n_e) + n_e coupling (C1*dne, dne=nh)
	fjac(1,1) = -n_h*g_hi + (-n_h*b_hi - a_hii*n_h)*n_e + C1*n_h
	! Apply the same transported-fraction constraints as the residual.
	call imposed_ionization_jacobian(ieq_cell, Neq, fjac)
	end subroutine jac_system_H

	! THE POSITIVE ROOT OF THE H I <-> H II BALANCE OF ion_system_H, in
	! closed form. With n_e = n(H II) = x n_H the row divided by n_H is
	!     g (1 - x) + n_H x (b (1 - x) - a x) = 0 ,
	! i.e. A x^2 + B x - g = 0 with A = n_H (a + b), B = g - n_H b. Its
	! roots have the product -g/A <= 0, so exactly one lies in [0, 1]:
	!     x = 2 g / (B + sqrt(B^2 + 4 A g))   (B >= 0),
	!     x = (-B + sqrt(B^2 + 4 A g)) / (2 A)   (B < 0),
	! each form free of cancellation on its side. Without photoionization
	! (g = 0) it is the collisional balance b/(a + b) where collisions
	! ionize (B < 0), and 0 otherwise. The starting point of the second
	! attempt of the H-only sweep (ionization_equilibrium).
	real*8 function hydrogen_ionization_balance_root() result(x)
	real*8 :: g_hi, a_hii, b_hi, n_h, qa, qb, disc
	g_hi  = ieq_cell%P_HI
	a_hii = ieq_cell%rchiiB
	n_h   = ieq_cell%nh
	b_hi  = ieq_cell%a_ion_HI
	qa    = n_h*(a_hii + b_hi)
	qb    = g_hi - n_h*b_hi
	disc  = sqrt(qb*qb + 4.0d0*qa*g_hi)
	x = 0.0d0
	if (qb .ge. 0.0d0) then
		if (qb + disc .gt. 0.0d0) x = 2.0d0*g_hi/(qb + disc)
	else
		if (qa .gt. 0.0d0) x = (disc - qb)/(2.0d0*qa)
	endif
	x = min(max(x, 0.0d0), 1.0d0)
	end function hydrogen_ionization_balance_root

	end module System_H

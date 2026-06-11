	module System_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine ion_system_HeH(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: g_hi,g_hei,g_heii		! Photoionization rates
	real*8  :: b_hi,b_hei,b_heii		! Collisional ionization rates
	real*8  :: a_hii,a_heii,a_heiii	    ! Recombination rates
   real*8  :: params(40)
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_hei,n_heii,n_heiii
	
	! Coefficients of the system

 	g_hi    = params(1)    ! = P_HI 
 	g_hei   = params(2)    ! = P_HeI
 	g_heii  = params(3)    ! = P_HeII 
 	a_hii   = params(4)    ! = rchiiB 
 	a_heii  = params(5)    ! = rcheiiB 
 	a_heiii = params(6)    ! = rcheiiiB 
 	n_h     = params(7)    ! = nh 
 	n_he    = params(8)    ! = nhe 
   b_hi    = params(9)    ! = a_ion_HI 
 	b_hei   = params(10)   ! = a_ion_HeI 
 	b_heii  = params(11)   ! = a_ion_HeII 
 	
 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_hei   = (1.0 - x(2) - x(3))*n_he 
 	n_heii  = x(2)*n_he
 	n_heiii = x(3)*n_he 
 	
 	! Electron density
   n_e = n_hii + n_heii + 2.0*n_heiii
      
      
      ! System of equations      
  	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e      
  	fvec(2) = n_hei*g_hei + (n_hei*b_hei - a_heii*n_heii)*n_e
  	fvec(3) = n_heii*g_heii + (n_heii*b_heii - a_heiii*n_heiii)*n_e
	
	return

	! End of subroutine
	end subroutine ion_system_HeH

	! Analytic Jacobian of ion_system_HeH (Task 2). Rows:
	!   f1 = n_hi*g_hi   + C1*n_e,  C1 = n_hi*b_hi   - a_hii*n_hii
	!   f2 = n_hei*g_hei + C2*n_e,  C2 = n_hei*b_hei - a_heii*n_heii
	!   f3 = n_heii*g_heii+ C3*n_e, C3 = n_heii*b_heii- a_heiii*n_heiii
	! n_e = n_hii + n_heii + 2 n_heiii ; dne = (nh, nhe, 2 nhe).
	subroutine jac_system_HeH(N_eq,x,fjac,params)
	integer :: N_eq
	real*8  :: x(N_eq), fjac(N_eq,N_eq), params(40)
	real*8  :: g_hi,g_hei,g_heii, a_hii,a_heii,a_heiii, b_hi,b_hei,b_heii
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii,n_hei,n_heii,n_heiii
	real*8  :: C1,C2,C3, dne(3)
	integer :: k

	g_hi    = params(1);  g_hei   = params(2);  g_heii  = params(3)
	a_hii   = params(4);  a_heii  = params(5);  a_heiii = params(6)
	n_h     = params(7);  n_he    = params(8)
	b_hi    = params(9);  b_hei   = params(10); b_heii  = params(11)

	n_hi    = (1.0-x(1))*n_h
	n_hii   = x(1)*n_h
	n_hei   = (1.0 - x(2) - x(3))*n_he
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	n_e     = n_hii + n_heii + 2.0*n_heiii

	C1 = n_hi*b_hi   - a_hii*n_hii
	C2 = n_hei*b_hei - a_heii*n_heii
	C3 = n_heii*b_heii - a_heiii*n_heiii
	dne(1) = n_h;  dne(2) = n_he;  dne(3) = 2.0*n_he

	! rank-1 n_e coupling: fjac(i,k) = C_i * dne(k)
	do k = 1, 3
		fjac(1,k) = C1*dne(k)
		fjac(2,k) = C2*dne(k)
		fjac(3,k) = C3*dne(k)
	enddo
	! add the local "direct" terms (photoionization + dC_i/dx_local * n_e)
	fjac(1,1) = fjac(1,1) - n_h*g_hi + (-n_h*b_hi - a_hii*n_h)*n_e
	fjac(2,2) = fjac(2,2) - n_he*g_hei + (-n_he*b_hei - a_heii*n_he)*n_e
	fjac(2,3) = fjac(2,3) - n_he*g_hei + (-n_he*b_hei)*n_e
	fjac(3,2) = fjac(3,2) + n_he*g_heii + (n_he*b_heii)*n_e
	fjac(3,3) = fjac(3,3) + (-a_heiii*n_he)*n_e
	end subroutine jac_system_HeH

	! End of module
	end module System_HeH

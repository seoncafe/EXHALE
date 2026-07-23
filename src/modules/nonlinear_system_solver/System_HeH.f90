	module System_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use ion_cell_state, only: ieq_cell
	use ion_residual_core, only: heh_rows, heh_crow, heh_jac_local
	use charge_exchange, only: he_h_cx_fvec, he_h_cx_jac

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
 	
 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_hei   = (1.0 - x(2) - x(3))*n_he 
 	n_heii  = x(2)*n_he
 	n_heiii = x(3)*n_he 
 	
 	! Electron density
   n_e = n_hii + n_heii + 2.0*n_heiii
      
      
      ! System of equations (standard H/He rows, shared helper)
	call heh_rows(fvec, n_hi, n_hii, n_hei, n_heii, n_heiii, n_e,  &
	              g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii,      &
	              b_hi, b_hei, b_heii)

	! He <-> H charge exchange (Huang Table 4 group B). He row is written
	! HeI->HeII (ionization) positive here, so he_row_sign = +1.
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,  &
	                  n_hi, n_hii, n_hei, n_heii, 1.0d0)

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

	g_hi    = ieq_cell%P_HI;  g_hei   = ieq_cell%P_HeI;  g_heii  = ieq_cell%P_HeII
	a_hii   = ieq_cell%rchiiB;  a_heii  = ieq_cell%rcheiiB;  a_heiii = ieq_cell%rcheiiiB
	n_h     = ieq_cell%nh;  n_he    = ieq_cell%nhe
	b_hi    = ieq_cell%a_ion_HI;  b_hei   = ieq_cell%a_ion_HeI; b_heii  = ieq_cell%a_ion_HeII

	n_hi    = (1.0-x(1))*n_h
	n_hii   = x(1)*n_h
	n_hei   = (1.0 - x(2) - x(3))*n_he
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	n_e     = n_hii + n_heii + 2.0*n_heiii

	call heh_crow(C1, C2, C3, n_hi, n_hii, n_hei, n_heii, n_heiii,  &
	              a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
	dne(1) = n_h;  dne(2) = n_he;  dne(3) = 2.0*n_he

	! rank-1 n_e coupling: fjac(i,k) = C_i * dne(k)
	do k = 1, 3
		fjac(1,k) = C1*dne(k)
		fjac(2,k) = C2*dne(k)
		fjac(3,k) = C3*dne(k)
	enddo
	! add the local "direct" terms (photoionization + dC_i/dx_local * n_e)
	call heh_jac_local(N_eq, fjac, n_h, n_he, n_e, g_hi, g_hei, g_heii,  &
	                   a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)

	! He <-> H charge-exchange Jacobian (rows 1,2; cols 1,2,3).
	call he_h_cx_jac(N_eq, fjac, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,  &
	                 n_h, n_he, n_hi, n_hii, n_hei, n_heii)
	end subroutine jac_system_HeH

	! End of module
	end module System_HeH

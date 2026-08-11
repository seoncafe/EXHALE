	module System_implicit_adv_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use ion_cell_state, only: adv_cell
	use charge_exchange, only: he_h_cx_fvec_adv

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

 	c1         = adv_cell%c1    ! = dr/v
 	xhi_old    = adv_cell%xhi_old    ! = nhi/nh
 	xhei_old   = adv_cell%xhei_old    ! = nhei/nhe
 	xheiii_old = adv_cell%xheiii_old    ! = nheiii/nhe
 	n_h        = adv_cell%nh    ! = nh
 	ghi        = adv_cell%P_HI    ! = P_HI
 	ghei       = adv_cell%P_HeI    ! = P_HeI = nh
 	gheii      = adv_cell%P_HeII    ! = P_HeII = nhe
   ahii       = adv_cell%rchiiB    ! = rchiiB
 	aheii      = adv_cell%rcheiiB   ! = rcheiiB
 	aheiii     = adv_cell%rcheiiiB   ! = rcheiiiB
	ionhi	   = adv_cell%a_ion_HI   ! = a_ion_HI
	ionhei     = adv_cell%a_ion_HeI   ! = a_ion_HeI
	ionheii    = adv_cell%a_ion_HeII   ! = a_ion_HEII
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (legacy, byte-identical); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead (the global HeH would misstate
	! n_e by the local separation factor).
	heh_loc    = adv_cell%heh_loc   ! = He/H (local when he_diffusion)

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xhei   = x(2)
   xheii  = 1.0 - x(2) - x(3)
	xheiii = x(3)

 	! Electron density, per H nucleus. The metal electrons (adv_cell%xe_metal,
 	! the same X+/X++ sum the equilibrium residual counts) are included: they
 	! dominate the electron budget of the shielded base, where the H/He
 	! ionized fractions are vanishingly small.
   xe = xhii + heh_loc*(xheii + 2.0*xheiii) + adv_cell%xe_metal
      
   ! System of equations
   ! The electron-impact ionization coefficients (ionhi, ionhei, ionheii,
   ! [cm^3 s^-1]) multiply the electron density xe*n_h, as in adv_implicit_H,
   ! adv_implicit_HeH_TR and the equilibrium rows (heh_rows: n_hi*b_hi*n_e).
  	fvec(1) =  xhi_old - xhi   					        		&
  	        + c1*(-(ghi + ionhi*xe*n_h)*xhi    + ahii*xhii*xe*n_h)

  	fvec(2) =  xhei_old - xhei 						  		    &
  		    + c1*(-(ghei + ionhei*xe*n_h)*xhei  + aheii*xheii*xe*n_h)

  	fvec(3) =  xheiii_old - xheiii					  			&
  	        + c1*((gheii + ionheii*xe*n_h)*xheii - aheiii*xheiii*xe*n_h)

	! He <-> H charge exchange (Huang Table 4 group B) on the H (row 1) and
	! He I (row 2) rows, both written neutral-gain positive here. xhei is the
	! He I fraction (n_HeI/n_he).
	call he_h_cx_fvec_adv(fvec, c1, xhi, xhii, xhei, xheii, heh_loc, n_h, &
	                      adv_cell%kcx_He0_Hp, adv_cell%kcx_Hep_H0)

	! End of subroutine
	end subroutine adv_implicit_HeH
	
	! End of module
	end module System_implicit_adv_HeH

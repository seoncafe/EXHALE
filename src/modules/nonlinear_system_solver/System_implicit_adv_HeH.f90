	module System_implicit_adv_HeH
	! Advection-corrected ionization system for H and He without the He 2^3S
	! metastable: one step of the backward differentiation formula for
	! v dx/dr = R(x) in the fractions, taken by post_process_adv (the step is
	! stated at the adv_rates type of ion_cell_state).
	
	use global_parameters
	use ion_cell_state, only: adv_cell
	use charge_exchange, only: he_h_cx_fvec_adv

	implicit none
	
	contains
	
	subroutine adv_implicit_HeH(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: xhi_hist,xheiS_hist,xheiii_hist
	real*8  :: ghi,ghei,gheii
	real*8  :: xhi,xhii
	real*8  :: xheiS,xheii,xheiii
	real*8  :: xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii
	real*8  :: ionhi,ionhei,ionheii
	real*8  :: heh_loc
   real*8  :: params(25)

	! Coefficients of the system

 	c1         = adv_cell%c1    ! = g*h_j/v_j, the rate weight of the step
 	xhi_hist    = adv_cell%xhi_hist    ! history of n_HI/n_H
 	xheiS_hist  = adv_cell%xheiS_hist   ! history of n_HeI/n_He (all He I is the singlet here)
 	xheiii_hist = adv_cell%xheiii_hist    ! history of n_HeIII/n_He
 	n_h        = adv_cell%nh    ! = nh
 	ghi        = adv_cell%P_HI    ! = P_HI
 	ghei       = adv_cell%P_HeI    ! = P_HeI
 	gheii      = adv_cell%P_HeII    ! = P_HeII
   ahii       = adv_cell%rchiiB    ! = rchiiB
 	aheii      = adv_cell%rcheiiB   ! = rcheiiB
 	aheiii     = adv_cell%rcheiiiB   ! = rcheiiiB
	ionhi	   = adv_cell%a_ion_HI   ! = a_ion_HI
	ionhei     = adv_cell%a_ion_HeI   ! = a_ion_HeI
	ionheii    = adv_cell%a_ion_HeII   ! = a_ion_HEII
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (the global HeH normally); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead (the global HeH would misstate
	! n_e by the local separation factor).
	heh_loc    = adv_cell%heh_loc   ! = He/H (local when he_diffusion)

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xheiS  = x(2)
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
  	fvec(1) =  xhi_hist - xhi   					        		&
  	        + c1*(-(ghi + ionhi*xe*n_h)*xhi    + ahii*xhii*xe*n_h)

  	fvec(2) =  xheiS_hist - xheiS 					  		    &
  		    + c1*(-(ghei + ionhei*xe*n_h)*xheiS + aheii*xheii*xe*n_h)

  	fvec(3) =  xheiii_hist - xheiii					  			&
  	        + c1*((gheii + ionheii*xe*n_h)*xheii - aheiii*xheiii*xe*n_h)

	! He <-> H charge exchange (Huang Table 4 group B) on the H (row 1) and
	! He I (row 2) rows, both written neutral-gain positive here, and He2+ +
	! H0 on rows 1 and 3 (He III gain positive). Without the
	! metastable the whole He I population is the ground singlet, so xheiS is
	! the He I fraction (n_HeI/n_he) the reaction sees.
	call he_h_cx_fvec_adv(fvec, c1, xhi, xhii, xheiS, xheii, xheiii,       &
	                      heh_loc, n_h, adv_cell%kcx_He0_Hp,              &
	                      adv_cell%kcx_Hep_H0, adv_cell%kcx_Hepp_H0)

	! End of subroutine
	end subroutine adv_implicit_HeH
	
	! End of module
	end module System_implicit_adv_HeH

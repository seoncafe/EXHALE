	module System_implicit_adv_HeH_TR
	! Advection-corrected ionization system for H, He and the He 2^3S
	! metastable: one step of the backward differentiation formula for
	! v_el dx/dr = R(x) in the fractions, each element's rows on the
	! velocity of its own nuclei, taken by post_process_adv (the step is
	! stated at the adv_rates type of ion_cell_state).
	!
	! Unknowns (all fractions; He fractions per He nucleus):
	!   x(1) = n_HI    /n_H
	!   x(2) = n_He(1^1S)/n_He      ground singlet He I
	!   x(3) = n_HeIII /n_He
	!   x(4) = n_He(2^3S)/n_He      metastable triplet He I
	! so the summed neutral helium is the SUM x(2)+x(4) and the once-ionized
	! fraction is 1 - x(2) - x(3) - x(4).
	!
	! The singlet is an unknown of its own rather than the difference
	! (summed He I) - (metastable): where helium is heavily ionized and the
	! little neutral helium that remains sits mostly in the metastable, that
	! difference loses every significant digit and finally evaluates to
	! exactly zero, which pins n(1^1S) at zero, degenerates the summed He I
	! row and leaves hybrd1 with a residual it cannot reduce. A sum has no
	! such failure mode: n(1^1S)/n(2^3S) may reach 1e-15 and both populations
	! keep their full relative accuracy.
	!
	! Row 2 is therefore the ground-singlet balance, i.e. the summed He I
	! balance minus the metastable balance of row 4. Its terms:
	!   gain  He+ recombination into the singlet ladder (aheii), the
	!         metastable returning to the ground state by electron
	!         collisions -- directly (q31g, the detailed-balance reverse of
	!         q13) or through 2^1S and 2^1P (q31a+q31b) --, by the
	!         2^3S -> 1^1S radiative decay (A31) and by the He(2^3S)+H0
	!         ionizing collisions (Q31, both the Penning and the associative
	!         branch leaving He in 1^1S), and the He+ + H0 charge exchange;
	!   loss  photoionization (ghei), electron-impact ionization (ionhei),
	!         collisional excitation into the metastable (q13), and the
	!         He0 + H+ charge exchange.
	
	use global_parameters
	use ion_cell_state, only: adv_cell
	use charge_exchange, only: he_h_cx_fvec_adv
	use Cooling_Coefficients, only: f_penning_HeI23S

	implicit none

	contains

	subroutine adv_implicit_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: xhi_hist,xheiS_hist,xheiii_hist,xheiTR_hist
	real*8  :: ghi,ghei,gheii,gheiTR
	real*8  :: xhi,xhii
	real*8  :: xheii,xheiii
	real*8  :: xheiTR, xheiS
	real*8  :: xe
	real*8  :: c1, c1_he
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii,aheiTR
	real*8  :: ionhi,ionhei,ionheii,ionheiTR
	real*8  :: heh_loc
   real*8  :: params(25)
   real*8  :: A31,q13,q31g,q31a,q31b,Q31
	
	! Coefficients of the system
 	c1         = adv_cell%c1    ! = g*h_j/v_H, the rate weight of the H row
 	c1_he      = adv_cell%c1_he ! = g*h_j/v_He, the rate weight of the He rows
 	xhi_hist    = adv_cell%xhi_hist    ! history of n_HI/n_H
 	xheiS_hist  = adv_cell%xheiS_hist   ! history of n_He(1^1S)/n_He (ground singlet)
 	xheiii_hist = adv_cell%xheiii_hist    ! history of n_HeIII/n_He
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
	ionheiTR   = adv_cell%a_ion_HeITR  ! = a_ion_HeITR (He 2^3S collisional ioniz.)
	aheiTR     = adv_cell%rcheiTR	  ! = rcheiTR
	A31	     = adv_cell%A31   ! = A31
	gheiTR     = adv_cell%P_HeITR   ! = P_HeITR
	q13        = adv_cell%q13   ! = q13
	q31g       = adv_cell%q31g  ! = q31g (reverse of q13)
	q31a	     = adv_cell%q31a   ! = q31a
	q31b	     = adv_cell%q31b   ! = q31b
	Q31 	     = adv_cell%Q31   ! = Q31
	xheiTR_hist = adv_cell%xheiTR_hist   ! history of n_He(2^3S)/n_He
	! Effective He/H for the electron density: packed as the global HeH by
	! post_process_adv (the global HeH normally); with He_diffusion the local,
	! radius-dependent nhe/nh is passed instead.
	heh_loc    = adv_cell%heh_loc   ! = He/H (local when he_diffusion)

	! Substitutions. The two neutral-He populations are solved for directly
	! (see the module header), so the once-ionized fraction is what is left
	! of the He nucleus after both of them.
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xheiS  = x(2)
	xheiii = x(3)
	xheiTR = x(4)
   xheii  = 1.0 - xheiS - xheiTR - xheiii

 	! Electron density, per H nucleus. The metal electrons (adv_cell%xe_metal,
 	! the same X+/X++ sum the equilibrium residual counts) are included: they
 	! dominate the electron budget of the shielded base, where the H/He
 	! ionized fractions are vanishingly small.
   xe = xhii + heh_loc*(xheii + 2.0*xheiii) + adv_cell%xe_metal
      
    ! System of equations      
  	! - f_penning_HeI23S*Q31*xheiTR*heh_loc*n_h*xhi: Penning loss of neutral
  	! H, He(2^3S)+H0 -> He(1^1S)+H+ + e-. xheiTR = n_23S/n_he and
  	! heh_loc = n_he/n_h, so xheiTR*heh_loc*n_h = n_23S; times Q31*xhi gives
  	! the volume rate divided by n_h (the normalization of this row). Q31 is
  	! the TOTAL ionization rate; the associative 10% returns its H atom when
  	! the HeH+ it makes recombines, so only the Penning branch is a net loss
  	! of neutral H here (see ion_residual_core).
  	fvec(1) =  xhi_hist - xhi + c1*(  		&
  	           - (ghi + ionhi*xe*n_h)*xhi 	&
  	           - f_penning_HeI23S*Q31*xheiTR*heh_loc*n_h*xhi &
  	           +  ahii*xhii*xe*n_h)
  	
  	! Ground-singlet He(1^1S) balance (module header). Losses: photoionization
  	! ghei, electron-impact ionization ionhei, and collisional excitation q13
  	! into the metastable. Gains: recombination of He+ into the singlet ladder
  	! (aheii; the aheiTR channel feeds the metastable and is charged to row 4),
  	! and every route by which the metastable returns to the ground state --
  	! electron collisions (q31g + q31a + q31b), the 2^3S -> 1^1S decay A31, and
  	! the He(2^3S)+H0 ionizing collisions Q31, whose Penning branch leaves
  	! He(1^1S)+H+ + e- and whose associative branch makes HeH+ that
  	! dissociatively recombines back to ground-state He (see ion_residual_core).
  	fvec(2) =  xheiS_hist - xheiS + c1_he*(                  &
  		       xheii*aheii*xe*n_h                            &
  		     - xheiS*(ghei + (ionhei + q13)*xe*n_h)          &
  		     + xheiTR*((q31g + q31a + q31b)*xe*n_h           &
  		               + A31 + xhi*Q31*n_h))
  	 	      	 	    
  	fvec(3) =  xheiii_hist - xheiii + c1_he*(  &
  	           (gheii + ionheii*xe*n_h)*xheii  & 
  	          - aheiii*xheiii*xe*n_h) 
  	           
	! - xheiTR*ionheiTR*xe*n_h: electron-impact ionization of He(2^3S) removes
	! the triplet (He(2^3S)+e- -> He+ + 2e-). The Q31 sink below takes the
	! TOTAL He(2^3S)+H ionization rate: Penning and associative both quench.
	fvec(4) =    xheiTR_hist - xheiTR + c1_he*(   &
		     - gheiTR*xheiTR				         &
		     + (xheii*aheiTR + xheiS*q13    	&
		     -  xheiTR*(q31g + q31a + q31b))*xe*n_h	&
		     - xheiTR*ionheiTR*xe*n_h           &
		     - xheiTR*(A31 + xhi*Q31*n_h))

	! He <-> H charge exchange (Huang Table 4 group B) on the H (row 1) and
	! He (row 2) rows, both written neutral-gain positive here, and He2+ +
	! H0 on rows 1 and 3 (He III gain positive). The helium
	! reactant of He + H+ -> He+ + H is the GROUND SINGLET, xheiS: the rate
	! charge_exchange::he_h_cx_rates forms for it, the detailed-balance
	! reverse of He+ + H -> He(1^1S) + H+, carries the barrier
	! exp(-12.75/T4), and 12.75e4 K = 10.99 eV is the ionization-potential
	! difference 24.587 - 13.598 eV of ground-state helium against hydrogen.
	! He(2^3S) lies 19.82 eV above the singlet, so its own charge exchange
	! with H+ is exothermic and has no such barrier; it is a different
	! reaction with a different rate, and it is carried neither here nor in
	! row 4, which gives the metastable no charge-exchange term at all.
	! Row 2 is the ground-singlet balance, so the singlet loss the pair
	! reports is charged to the row that owns the singlet, exactly as the
	! equilibrium systems charge it.
	call he_h_cx_fvec_adv(fvec, c1, c1_he, xhi, xhii, xheiS, xheii, xheiii,       &
	                      heh_loc, n_h, adv_cell%kcx_He0_Hp,              &
	                      adv_cell%kcx_Hep_H0, adv_cell%kcx_Hepp_H0)

	! End of subroutine
	end subroutine adv_implicit_HeH_TR
	
	! End of module
	end module System_implicit_adv_HeH_TR

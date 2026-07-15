	module post_processing
	! Subroutine to correct the output ionization profiles
	!	taking into account the ionization term
	
	use global_parameters
	use params_idx, only: IPAH_C1, IPAH_XHI, IPAH_NH, IPAH_PHI, IPAH_AHII,     &
	                      IPAH_BHI, IPA_C1, IPA_XHI, IPA_XHEI, IPA_XHEIII,      &
	                      IPA_NH, IPA_PHI, IPA_PHEI, IPA_PHEII, IPA_AHII,       &
	                      IPA_AHEII, IPA_AHEIII, IPA_BHI, IPA_BHEI, IPA_BHEII,  &
	                      IPA_HEH, IPAT_ATR, IPAT_A31, IPAT_PTR, IPAT_Q13,      &
	                      IPAT_Q31A, IPAT_Q31B, IPAT_Q31, IPAT_XTR, IPAT_HEH,   &
	                      IPE_NH, IPE_NHE, IPT_NHI, IPT_NHII, IPT_NHEI,         &
	                      IPT_NHEII, IPT_NHEIII, IPT_MUP, IPT_MUM, IPT_RHOV,    &
	                      IPT_COEFF, IPT_DR, IPT_TOLD, IPT_HEAOLD
	use species_table, only: n_mion, n_melem, melem_i0, melem_top
	use utils
	use System_implicit_adv_H
	use System_implicit_adv_HeH
	use System_implicit_adv_HeH_TR
	use System_HeH_metals, only: ion_system_HeH_metals, set_metal_coeffs
	use charge_exchange,   only: cx_set_cell
	use Cooling_Coefficients
	use utils_ion_eq
	use output_write
	use equation_T
	use opacity_models           ! opacity_pT_factor for the 'P' model

	implicit none

	! Advection-corrected H/He ionized fractions for the current cell, pinned
	! while the metal re-solve (pp_metals=2) adjusts only the metal stages.
	! Set per cell before each ion_system_metals_pp / hybrd1 call.
	real*8, save :: pp_xHII_fix   = 0.0d0
	real*8, save :: pp_xHeII_fix  = 0.0d0
	real*8, save :: pp_xHeIII_fix = 0.0d0

	contains
	
	subroutine post_process_adv(rho,v,p,T_in,heat,cool,eta,   &
                                  nhi_in,nhii_in,		    &
                                  nhei_in,nheii_in,nheiii_in,   &
                                  nheiTR_in, nm_in)


	real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
   real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi_in,nhii_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei_in,nheii_in,   &
      							  				        nheiii_in,nheiTR_in
   ! Converged equilibrium metal densities (dimensionless, n0 units), used
   ! by the metal-aware post-process modes (pp_metal_mode = 1 frozen, 2 re-solve).
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm_in
	
	integer i,j,k
	integer :: n_pp_reject       ! cell-by-cell T solves rejected as non-physical
	integer :: Neq_adv,lwa_adv   ! advection system size (metal-independent)
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K,p_out,T_out     ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhe,ne,n_tot 
	
	! Discard scratch: distinct locals for the discarded intent(out) slots so no
	! two out-arguments in the same call alias one another.
	real*8, dimension(1-Ng:N+Ng) ::  dum_v1,dum_v2,dum_v3,dum_v4,dum_v5,dum_v6

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   ! Metal photoionization rates for each ion (filled by PH_heat_HHe; used for the
   ! metal re-solve in mode 2, otherwise discarded).
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
   ! Working metal densities driving the post-process heating/cooling [cgs].
   ! Built from nm_in per pp_metal_mode: 0 -> zero (metal-free, legacy _adv),
   ! 1 -> frozen eq metals, 2 -> re-solved. nm_out is its dimensionless (n0)
   ! copy written to the _adv ion-species file.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm_w, nm_out

   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR
   ! Metal recombination/ionization rates for each ion returned by eval_cool.
   ! In the re-solve mode (pp_metals=2) they feed the cell-by-cell metal
   ! ionization-balance solve; in the frozen mode they are discarded.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m_pp,aion_m_pp
   ! Each element's total metal density (sum of stages, held fixed across the
   ! re-solve since ionization only redistributes among an element's stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot_pp
   ! Each element's metal coefficients handed to set_metal_coeffs for the
   ! re-solve (canonical element order); built per cell from the 2-D rates.
   real*8, dimension(n_melem) ::  meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1,  &
                                  meg_a1,meg_a2
   integer, dimension(n_melem) ::  meg_top
   integer :: i0,top,im
      
   real*8, dimension(1-Ng:N+Ng) ::  q13,q31a,q31b,Q31
	real*8 :: A31
 	
 	! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII 
      
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  theat,tcool
      
 	! Updated species densities
   real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
   real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii,nheiTR,nheiS
   real*8, dimension(1-Ng:N+Ng) :: mmw
   real*8, dimension(1-Ng:N+Ng) :: nhi_w,nhii_w
   real*8, dimension(1-Ng:N+Ng) :: nhei_w,nheii_w,nheiii_w,nheiTR_w
      							  
	      
	      
   real*8 :: TT                              ! Temperature component
   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR       ! Photoionization rates
   real*8 :: deltal                          ! Optical depth
   real*8 :: dr	                           ! Grid spacing
   real*8 :: iup_1,ilo_1,        &           ! Photoheating integral variables
             iup_15,ilo_15,      &
             iup_2,ilo_2,        &
             iup_TR,ilo_TR,	   &
             iup_f,ilo_f 
   real*8 :: elo,eup                         ! Energy parameters
   real*8 :: tol,dpmpar                      ! Equilibrium system setup
   real*8 :: Hea_1 		          	         ! Heating rates
   real*8 :: brem,coex,coio,reco             ! Cooling rates
   real*8 :: iup_H,ilo_H                     ! Heating rate integral variables         
      
      
	! Substitution in the ODE solution
	real*8 :: As


   ! 60 entries to match the ion_system_HeH_metals / ion_system_metals_pp
   ! dummy length used by the pp_metals=2 re-solve (only 1-11 are consumed;
   ! the advection systems use 1-22).
   real*8, dimension(60) :: params
   real*8, dimension(40) :: paramsT
      
	real*8 :: rhop,rhom,vp,mum,mup,vm
	real*8 :: sys_sol_T(1), sys_x_T(1)
   real*8 :: wa_T(8)
   logical :: brent_ok                       ! Task 1: Brent T-solve status
      
      
   !----------------------------------------------------------!      
      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

   ! Advection system size. The advection-correction systems
   ! (adv_implicit_H/HeH/HeH_TR) solve only H/He fractions and never
   ! the metal stages, so the count is independent of the (possibly
   ! larger) global N_eq used by the equilibrium solver. Using N_eq
   ! here would feed hybrd1 uninitialized fvec/x entries for the metal
   ! slots. The matching MINPACK workspace size is lwa_adv.
   if (.not.thereis_He) then
      Neq_adv = 1
   else if (thereis_HeITR) then
      Neq_adv = 4
   else
      Neq_adv = 3
   endif
   lwa_adv = (Neq_adv*(3*Neq_adv + 13))/2

   !----------------------------------!
	
	! Preliminary profiles extraction
	
	! Dimensional total number density profile and temperature
	T_K = T_in*T0
	         
   ! Initialize vectors
   nhi    = nhi_in*n0
	nhii   = nhii_in*n0
   if (thereis_He) then
		nhei   = nhei_in*n0
		nheii  = nheii_in*n0
		nheiii = nheiii_in*n0
		! Zero-init the triplet so the innermost ghost cell (1-Ng), which
		! the advection loop below never assigns when HeITR is off, does
		! not write uninitialized memory to the output column.
		nheiTR = 0.0
		if (thereis_HeITR) nheiTR = nheiTR_in*n0
	endif

	!----------------------------------!

	! Select the metal treatment for the post-process (set via 'pp_metals'
	! in metals.inp -> pp_metal_mode). Metals are off => zero either way.
	!   0 metal-free (legacy), 1 frozen eq metals, 2 re-solve.
	! pp_metal_on (module flag in equation_T) tells the cell-by-cell temperature
	! solve to include the metal cooling/brem/n_e terms.
	if (thereis_metals .and. pp_metal_mode >= 1) then
		nm_w        = nm_in*n0          ! eq metal densities [cgs] (frozen or seed)
		pp_metal_on = .true.
	else
		nm_w        = 0.0d0
		pp_metal_on = .false.
	endif

	! Re-solve mode (pp_metals=2): the metal ionization balance is recomputed
	! per cell inside the post-process loop at the advection-corrected H/He and
	! the post-process temperature (see the block after the H/He advection
	! solve). The element totals are conserved by ionization (only the stage
	! split changes), so freeze each element's total from the eq densities and
	! reuse it as a constant throughout. nm_w starts at the eq split and is
	! overwritten with the re-solved split each iteration.
	if (pp_metal_mode == 2 .and. thereis_metals .and. thereis_He) then
		do im = 1,n_melem
			i0 = melem_i0(im)
			if (melem_top(im) >= 2) then
				nm_tot_pp(:,im) = nm_w(:,i0) + nm_w(:,i0+1) + nm_w(:,i0+2)
			else
				nm_tot_pp(:,im) = nm_w(:,i0) + nm_w(:,i0+1)
			endif
		enddo
	endif

	!----------------------------------!

	! Iterate the post processing
	do k = 1,10	! Usually 10 gives a good convergence
	
	! Use singlet if included
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii 
	nhe = nheiS + nheii + nheiii 
	if (thereis_HeITR) nhe = nhe + nheiTR

	! Free electron density (assuming overall neutrality; nm_w adds the
	! metal electrons under the eos_metals policy)
	call calc_ne(nhii,nheii,nheiii,ne,nm_w)

	! Cell-by-cell opacity pressure factor ('P' model; =1 otherwise)
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

   ! Calculate the photoionization rates

	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w,           &
					 P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
					 dum_v1,dum_v2)
  	else
	  	call PH_heat_H(nhi,P_HI,dum_v1,dum_v2)
  	endif

   !---- Recombination rates ----!

	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm_w,            &
	  			   rchiiB,rcheiiB,rcheiiiB, rec_m_pp,             &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m_pp,          &
				   dum_v1)

 	
 	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten

	endif
	
   !----------------------------------!
      
   ! Evolve species including the advection term in the 
   ! 	ODE form
   ! Note: we are using point values here instead of 
   ! 	volume averages; they agree up to O(dr^2)
      
   ! The ionization fraction at the inner boundary are taken 
	!	from the input vectors (completely neutral atmosphere)
     

   ! Loop to solve the differential equation
   ! It is implicitly assumed that the velocity fields does 
   !	not change by including the advection term
      
      
	      
   if (.not.thereis_He) then
	      
		do j = 2-Ng,N+Ng 
			! Option-(c) guard: in metal-cooled mode keep the converged eq ionization
			! wherever the flow is not a clean outflow (v<=0: the breathing/inflow
			! base). The advection correction assumes outflow upwinding and is invalid
			! there; skipping it also breaks the upwind cascade that produces the
			! spurious base temperature spike. Metals-off mode is left byte-identical.
			if (pp_metal_on .and. (v(j) <= 0.0d0 .or. v(j-1) <= 0.0d0)) then
				nhi(j)  = nhi_in(j)*n0
				nhii(j) = nhii_in(j)*n0
				cycle
			endif
	
			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			params(IPAH_C1) = As
			params(IPAH_XHI) = nhi(j-1)/nh(j-1)
			params(IPAH_NH) = nh(j)
			params(IPAH_PHI) = P_HI(j)
			params(IPAH_AHII) = rchiiB(j)
			params(IPAH_BHI) = a_ion_HI(j)
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j)     
			
			! Call hybrd1 routine (from minpack)
			call hybrd1(adv_implicit_H,Neq_adv,sys_x,sys_sol,   &
						tol,info,wa,lwa_adv,params)
			
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)

      	enddo
      	
		   ! Force condition of zero helium
		   nhei   = 0.0
		   nheii  = 0.0
		   nheiii = 0.0			
		   nheiTR = 0.0

	else
		
		do j = 2-Ng,N+Ng 
			! Option-(c) guard (see the no-He branch): in metal-cooled mode the
			! non-outflow base (v<=0) keeps the converged eq H/He ionization.
			if (pp_metal_on .and. (v(j) <= 0.0d0 .or. v(j-1) <= 0.0d0)) then
				nhi(j)    = nhi_in(j)*n0
				nhii(j)   = nhii_in(j)*n0
				nhei(j)   = nhei_in(j)*n0
				nheii(j)  = nheii_in(j)*n0
				nheiii(j) = nheiii_in(j)*n0
				nheiTR(j) = 0.0
				if (thereis_HeITR) nheiTR(j) = nheiTR_in(j)*n0
				cycle
			endif
	
			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			params(IPA_C1)  = As
			params(IPA_XHI)  = nhi(j-1)/nh(j-1)
			params(IPA_XHEI)  = nhei(j-1)/nhe(j-1)
			params(IPA_XHEIII)  = nheiii(j-1)/nhe(j-1)
			params(IPA_NH)  = nh(j)
			params(IPA_PHI)  = P_HI(j)
			params(IPA_PHEI)  = P_HeI(j)
			params(IPA_PHEII)  = P_HeII(j)
			params(IPA_AHII)  = rchiiB(j)
			params(IPA_AHEII) = rcheiiB(j)
			params(IPA_AHEIII) = rcheiiiB(j)
			params(IPA_BHI) = a_ion_HI(j)
			params(IPA_BHEI) = a_ion_HeI(j)
			params(IPA_BHEII) = a_ion_HeII(j)
			! Effective He/H for the electron density inside the adv system:
			! the global HeH normally (byte-identical legacy), the local
			! (diffused) nhe/nh when He_diffusion is on.
			if (he_diffusion) then
				params(IPA_HEH) = nhe(j)/max(nh(j),1.0d-30)
			else
				params(IPA_HEH) = HeH
			endif

			! Add more if HeITR is present
			if (thereis_HeITR) then 
				params(IPAT_ATR) = rcheiTR(j)
				params(IPAT_A31) = A31
				params(IPAT_PTR) = P_HeITR(j)
				params(IPAT_Q13) = q13(j)
				params(IPAT_Q31A) = q31a(j)
				params(IPAT_Q31B) = q31b(j)
				params(IPAT_Q31) = Q31(j)
				params(IPAT_XTR) = nheiTR(j-1)/nhe(j-1)
				! Effective He/H for the electron density (see non-TR block).
				if (he_diffusion) then
					params(IPAT_HEH) = nhe(j)/max(nh(j),1.0d-30)
				else
					params(IPAT_HEH) = HeH
				endif
			endif
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j) 
			sys_x(2) = nhei(j)/nhe(j)
			sys_x(3) = nheiii(j)/nhe(j)
			if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j) 
			
			! Call hybrd1 routine (from minpack)
			if (thereis_HeITR) then 
				call hybrd1(adv_implicit_HeH_TR,Neq_adv,sys_x,sys_sol,   &
							tol,info,wa,lwa_adv,params)
			else
				call hybrd1(adv_implicit_HeH,Neq_adv,sys_x,sys_sol,   &
							tol,info,wa,lwa_adv,params)
			endif
				
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)
			nhei(j)   = sys_x(2)*nhe(j)
			nheii(j)  = (1.0 - sys_x(2) - sys_x(3))*nhe(j)
			nheiii(j) = sys_x(3)*nhe(j)
			if (thereis_HeITR) then
				nheiTR(j) = sys_x(4)*nhe(j) 
			else
				nheiTR(j) = 0.0
			endif
			
		enddo
		
	endif ! End if thereis_He
	      
      
   !---------------------------------------------------------
      
   !------- Fix stationarity of new pressure profile -------!
      
   !---- Update densities and temperature ----!
	
	! Number densities      
   nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii 
	nhe = nheiS + nheii + nheiii
	if (thereis_HeITR) nhe = nhe + nheiTR

   ! Total number density (incl. metal nuclei under eos_metals)
   call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot,nm_w)

   ! Free electron density (assuming overall neutrality; incl. metal
   ! electrons under eos_metals)
   call calc_ne(nhii,nheii,nheiii,ne,nm_w)

	!----------------------------------!

	!---- Re-solve metal ionization (pp_metals=2) ----!
	! Recompute the metal ionization balance per cell at the advection-
	! corrected H/He and the post-process temperature T_K. The full coupled
	! H/He+metal residual (ion_system_HeH_metals) is reused, but rows 1-3 are
	! pinned by the ion_system_metals_pp wrapper so only the metal stages move;
	! the electron density inside the residual then combines the fixed H/He
	! with the re-solved metals (matching the equilibrium solver's n_e). The
	! photo/collisional/recombination rates (P_m, aion_m_pp, rec_m_pp) and the
	! charge-exchange couplings to H+/He+ come from the current pass, so the
	! k-loop drives metals, H/He and T to a joint fixed point. Element totals
	! are conserved (nm_tot_pp), so only the stage split is updated.
	if (pp_metal_mode == 2 .and. thereis_metals .and. thereis_He) then
		do j = 1-Ng,N+Ng
			if (nh(j) <= 0.0d0) cycle

			! Charge-exchange rate coefficients for this cell temperature.
			call cx_set_cell(T_K(j))

			! Each element's metal coefficients (canonical order) from the 2-D
			! photo/collisional/recombination rate arrays.
			do im = 1,n_melem
				i0  = melem_i0(im)
				top = melem_top(im)
				meg_top(im)  = top
				meg_ntot(im) = nm_tot_pp(j,im)
				meg_g0(im)   = P_m(j,i0)
				meg_b0(im)   = aion_m_pp(j,i0)
				meg_a1(im)   = rec_m_pp(j,i0+1)
				if (top >= 2) then
					meg_g1(im) = P_m(j,i0+1)
					meg_b1(im) = aion_m_pp(j,i0+1)
					meg_a2(im) = rec_m_pp(j,i0+2)
				else
					meg_g1(im) = 0.0d0
					meg_b1(im) = 0.0d0
					meg_a2(im) = 0.0d0
				endif
			enddo
			call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1,   &
			                      meg_b0, meg_b1, meg_a1, meg_a2, meg_top)

			! Pin the advection-corrected H/He fractions for the wrapper.
			pp_xHII_fix   = nhii(j)/nh(j)
			pp_xHeII_fix  = nheii(j)/max(nhe(j),1.0d-300)
			pp_xHeIII_fix = nheiii(j)/max(nhe(j),1.0d-300)

			! H/He params for the residual: only n_h, n_he are consumed once
			! rows 1-3 are pinned (the H/He rates drop out).
			params    = 0.0d0
			params(IPE_NH) = nh(j)
			params(IPE_NHE) = nhe(j)

			! Initial guess: pinned H/He fractions + current metal split.
			sys_x(1) = pp_xHII_fix
			sys_x(2) = pp_xHeII_fix
			sys_x(3) = pp_xHeIII_fix
			do im = 1,n_melem
				i0 = melem_i0(im)
				sys_x(4+2*(im-1)) = nm_w(j,i0+1)/max(nm_tot_pp(j,im),1.0d-30)
				if (melem_top(im) >= 2) then
					sys_x(5+2*(im-1)) = nm_w(j,i0+2)/max(nm_tot_pp(j,im),1.0d-30)
				else
					sys_x(5+2*(im-1)) = 0.0d0
				endif
			enddo

			call hybrd1(ion_system_metals_pp,N_eq,sys_x,sys_sol,   &
			            tol,info,wa,lwa,params)

			! Extract the re-solved stage split (element totals conserved).
			do im = 1,n_melem
				i0 = melem_i0(im)
				nm_w(j,i0+1) = nm_tot_pp(j,im)*sys_x(4+2*(im-1))
				if (melem_top(im) >= 2) then
					nm_w(j,i0+2) = nm_tot_pp(j,im)*sys_x(5+2*(im-1))
					nm_w(j,i0)   = nm_tot_pp(j,im)                          &
					             *(1.0d0 - sys_x(4+2*(im-1)) - sys_x(5+2*(im-1)))
				else
					nm_w(j,i0)   = nm_tot_pp(j,im)*(1.0d0 - sys_x(4+2*(im-1)))
				endif
			enddo
		enddo
	endif

	!----------------------------------!

	!---- Update photoheating rate ----!
	
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w,           &
		                 dum_v1,dum_v2,dum_v3,dum_v4, P_m,          &
		                 theat,dum_v5)
  	else
	  	call PH_heat_H(nhi,dum_v1,theat,dum_v2)
  	endif

	! Adimensionalize
	theat = theat/q0		

	!----------------------------------!
	
	!---- Solve stationary energy equation ----!
	! This procedure uses the same velocity profile and 
	!	the ionization profile after the advection correction
	
	! Initialize temperature at ghost cells
	T_out = T_K/T0
	
	! Calculate mean molecular weight. nm_w adds the metal mass/nuclei under
	! the eos_metals policy, so the _adv temperature solve uses the same
	! composition as the main loop (ne above already carries the metal
	! electrons via calc_ne).
	call calc_mmw(nh,nhe,ne,mmw,nm_w)

	! Count cell-by-cell temperature solves rejected as non-physical (metal modes).
	n_pp_reject = 0

	do j = 3-Ng,N+Ng ! Start from first computational cell
		
		! Substitutions
		rhop = rho(j)
		rhom = rho(j-1)
		vm = v(j-1) 
		vp = v(j)
		! Option-(c) guard: in metal-cooled mode the non-outflow base (v<=0) keeps
		! the converged eq temperature; the advection-corrected energy solve assumes
		! outflow and would otherwise land on the spurious hot root that cascades up.
		if (pp_metal_on .and. (vp <= 0.0d0 .or. vm <= 0.0d0)) then
			T_out(j) = T_in(j)
			cycle
		endif
		dr = r(j) - r(j-1)
		mum = mmw(j-1)
		mup = mmw(j)

	 	!--- Solve equation for temperature implicitly ---!
		
		! Parameters
		paramsT(IPT_NHI)  = nhi(j)
	 	paramsT(IPT_NHII)  = nhii(j)
	 	paramsT(IPT_NHEI)  = nhei(j)
	 	paramsT(IPT_NHEII)  = nheii(j)
	 	paramsT(IPT_NHEIII)  = nheiii(j)
	 	paramsT(IPT_MUP)  = mmw(j)
	 	paramsT(IPT_MUM)  = mmw(j-1)
	 	paramsT(IPT_RHOV)  = rhop*vp
	 	paramsT(IPT_COEFF)  = mum*vp*(rhop-rhom)
	 	paramsT(IPT_DR) = dr
	 	paramsT(IPT_TOLD) = T_out(j-1)
	 	paramsT(IPT_HEAOLD) = theat(j)
	 	! Metal densities for this cell [cgs] go through the equation_T module
	 	! array (the 27-ion vector does not fit params). pp_metal_on gates
	 	! whether T_equation adds the metal cooling/brem/n_e terms.
	 	pp_nm_cell(:) = nm_w(j,:)

	 	! Initial guess of solution
		sys_x_T(1) = T_out(j)

	 	! Task 1: solve the scalar energy equation by bracketing the physical
	 	! (lowest) root + Brent when metal cooling is on (the default). The
	 	! metal-cooled residual is non-monotone and has a second, spurious *hot*
	 	! root that a Newton/Powell solve (hybrd1) could land on; bracketing from
	 	! below selects the physical root structurally. Fall back to the converged
	 	! eq T if no bracket is found. With "Brent solver: False" (use_brent_tsolve
	 	! = .false.) the legacy MINPACK solve + 2x-band reject is used instead.
	 	! Metals-off always keeps the original MINPACK solve (monotone residual,
	 	! byte-identical).
	 	if (pp_metal_on .and. use_brent_tsolve) then
	 		call solve_T_brent(paramsT, T_in(j), sys_x_T(1), brent_ok)
	 		if (.not. brent_ok) then
	 			sys_x_T(1)  = T_in(j)
	 			n_pp_reject = n_pp_reject + 1
	 		endif
	 	else
	 		! Legacy MINPACK solve.
	 		call hybrd1(T_equation,1,sys_x_T,sys_sol_T,   &
	 		            tol,info,wa_T,8,paramsT)
	 		! With metal cooling, reject a non-physical / out-of-band root (the
	 		! spurious hot root) and fall back to eq T (the legacy guard).
	 		if (pp_metal_on) then
	 			if (.not. (sys_x_T(1) > 0.0d0)  .or.   &
	 			    sys_x_T(1) > 2.0d0*T_in(j)  .or.   &
	 			    sys_x_T(1) < 0.5d0*T_in(j)) then
	 				sys_x_T(1)  = T_in(j)
	 				n_pp_reject = n_pp_reject + 1
	 			endif
	 		endif
	 	endif

		! Extract solution profiles
		T_out(j) = sys_x_T(1)
	 	
	enddo
	
	! Update pressure and temperature
	p_out = (n_tot + ne)/n0*T_out
	T_K = T_out*T0
      
	enddo ! End loop on post processing

	! Report how many base cells fell back to the eq temperature (metal modes).
	if (pp_metal_on .and. n_pp_reject > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) metal-mode T solve: ', &
		   n_pp_reject, ' of ', N+2*Ng-2,                                     &
		   ' cells fell back to eq T (stiff base band).'
	endif

   ! ---------------------------- !
      
	!---- Update cooling rates ----!	
	
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm_w,            &
	  			   dum_v1,dum_v2,dum_v3, rec_m_pp,                   &
				   dum_v4,dum_v5,dum_v6, aion_m_pp,                      &
				   tcool)

	! Adimensionalize
	tcool = tcool/q0

	!----------------------------------!
 
   ! Adimensionalize ion densities before writing
   nhi_w    = nhi/n0
   nhii_w   = nhii/n0
   nhei_w   = nhei/n0
   nheii_w  = nheii/n0
   nheiii_w = nheiii/n0
   nheiTR_w = nheiTR/n0

   !----------------------------------!
      
   ! Write updated thermodynamic and ionization profiles. Metals are written
   ! in dimensionless (n0) units, matching the other species; zero in the
   ! metal-free mode, frozen/re-solved eq densities otherwise.
   nm_out = nm_w/n0
   call write_output(rho,v,p_out,T_out,theat,tcool,eta,                &
                     nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,              &
                     nheiTR_w,nm_out,'ad')

	! End of subroutine
	end subroutine post_process_adv

	! ----------------------------------------------------------------- !

	! Residual for the metal-only re-solve (pp_metals=2). Reuses the full
	! coupled H/He+metal balance (ion_system_HeH_metals, including the Huang
	! Table-4 charge exchange), then overwrites the three H/He rows with
	! identity equations that pin x(1..3) to the advection-corrected fractions
	! stored in pp_xHII_fix / pp_xHeII_fix / pp_xHeIII_fix. hybrd1 therefore
	! leaves H/He fixed and moves only the metal stages, while the electron
	! density and charge exchange inside the coupled residual still see the
	! correct (fixed) H/He densities. The cell-by-cell metal coefficients must be
	! loaded via set_metal_coeffs and the fractions pinned before each call.
	subroutine ion_system_metals_pp(N_in,x,fvec,iflag,params)
	integer :: N_in,iflag
	real*8  :: x(N_in),fvec(N_in)
	real*8  :: params(60)

	call ion_system_HeH_metals(N_in,x,fvec,iflag,params)

	! Pin the H/He fractions (rows 1-3) to the advection-corrected values.
	fvec(1) = x(1) - pp_xHII_fix
	fvec(2) = x(2) - pp_xHeII_fix
	fvec(3) = x(3) - pp_xHeIII_fix

	end subroutine ion_system_metals_pp

	! End of module
	end module post_processing

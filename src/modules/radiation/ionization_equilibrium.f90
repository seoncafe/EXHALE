   module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature

	use global_parameters
   use species_table, only: n_mion, mion_fsp, n_melem, melem_i0,        &
                            melem_top, mion_stage
   use utils
   use utils_ion_eq
   use Cooling_Coefficients      ! eval_cool, recombination/ionization rates
   use System_HeH                ! Equilibrium equations
	use System_HeH_TR
   use System_HeH_metals
   use System_HeH_TR_metals      ! merged He-triplet + metals system
   use charge_exchange, only: cx_set_cell, cx_metal_base   ! Huang Table 4 charge exchange
   use System_H
   use newton_solver, only: solve_ieq   ! Task 2: analytic-Jacobian Newton (+ hybrd1 fallback)
   use opacity_models            ! opacity_pT_factor for the 'P' model

   implicit none

	contains 
	
	subroutine ioniz_eq(T_in,n_in,f_sp_in,n_out, &
      	 		  f_sp_out,heat_out,cool_out,q)
      	 		  
	integer :: j,im
	logical :: usednt                     ! Task 2: Newton-vs-fallback flag

	real*8, dimension(1-Ng:N+Ng),   intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,nheiTR,  &
	                                 ne,n_in_dim
   ! Metal ion densities in canonical species_table order (col im maps
   ! to f_sp column mion_fsp(im)); used throughout in place of named
   ! per-ion scalars so the driver scales with the number of metals.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm
   ! Total density of each metal element (canonical element order),
   ! held constant across the cell sweep (sum of its three stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot

   ! Per-ion metal recombination and collisional ionization rates from
   ! eval_cool (canonical order); bridged to the named rc*/a_ion_*
   ! scalars below for the params packing.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m,aion_m

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   ! Per-ion metal photoionization rates (canonical order) from PH_heat.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
                       	
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  heat,cool    
                 
   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR

	real*8, dimension(1-Ng:N+Ng) :: q13,q31a,q31b
	real*8 :: A31,Q31

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII

	! Equilibrium system setup
   real*8 :: tol,dpmpar
   real*8, dimension(60) :: params
   ! Per-element metal coefficients for ion_system_HeH_metals (canonical order),
   ! built per cell from the 2D rate arrays and handed to set_metal_coeffs.
   real*8, dimension(n_melem) :: meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1, &
                                 meg_a1,meg_a2
   ! Highest stage per element handed to set_metal_coeffs (2 = three-stage,
   ! 1 = two-stage); see species_table::melem_top.
   integer, dimension(n_melem) :: meg_top
   integer :: i0,top
   ! Base index of the first metal element's X+ unknown in sys_x: 4 normally,
   ! 5 when the He triplet occupies x(4) (merged HeITR+metals system).
   integer :: mbase
      
   ! Output density
   real*8, dimension(1-Ng:N+Ng),intent(out) :: n_out 
      
   ! Output heating,cooling and absorbed energy 
   real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out,q 
      
   ! Output species fractions
   real*8, dimension(1-Ng:N+Ng,n_species),intent(out) :: f_sp_out
	
   !----------------------------------------------------------!      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

	!----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	n_in_dim = n_in*n0
	T_K      = T_in*T0
		
	! Extract species profiles
	nhi    = f_sp_in(:,1)*n_in_dim    ! HI
	nhii   = f_sp_in(:,2)*n_in_dim    ! HII
	if (thereis_He) then
		nhei   = f_sp_in(:,3)*n_in_dim    ! HeI
		nheii  = f_sp_in(:,4)*n_in_dim    ! HeII
		nheiii = f_sp_in(:,5)*n_in_dim    ! HeIII
		nheiTR = f_sp_in(:,6)*n_in_dim    ! HeITR
	else
		! Enforce condition of zero helium
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0
	endif

    ! Metal ion densities in canonical species_table order (col im maps to
    ! f_sp column mion_fsp(im)).
    do im = 1,n_mion
       nm(:,im) = f_sp_in(:,mion_fsp(im))*n_in_dim
    enddo

    ! Total density of each metal element (sum of its three stages, in the
    ! same neutral+singly+doubly order as the original nC=nci+ncii+nciii).
    do im = 1,n_melem
       i0 = melem_i0(im)
       if (melem_top(im) .ge. 2) then
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1) + nm(:,i0+2)
       else
          ! Two-stage element: neutral + singly ionized only.
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1)
       endif
    enddo
	
	! Total number densities      
	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii
	
	! Free electron density (assuming overall neutrality; nm adds the
	! metal electrons under the eos_metals policy)
	call calc_ne(nhii,nheii,nheiii,ne,nm)

	! Per-cell pressure-broadening factor for the opacity ('P' model).
	! opacity_pT_factor returns 1.0 for all other models, so opa_pf=1
	! and the column densities are unchanged (bit-identical).
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

	!----------------------------------!

    !---- Photoionization and photoheating ----!
      
	if (thereis_He) then
      	call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm,             &
      			     P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
      			     heat,q)
	else
		call PH_heat_H(nhi,P_HI,heat,q)
		P_m = 0.0
	endif

	! Phase 3a excited-H feedback (zero unless use_excited_H). The Balmer
	! photoionization of H(n=2) adds an effective HI photoionization rate
	! [s^-1] (proton source), and the Balmer photoelectric (+ optional
	! collisional de-excitation) heating adds to the photoheating rate
	! [erg cm^-3 s^-1]. Both global arrays are filled from the previous
	! converged outer pass by excited_hydrogen::excited_H_update, so the
	! coupling is decoupled from the hydro sub-step. Default-off => no change.
	if (use_excited_H) then
		gph_ground_HI = P_HI         ! capture pure ground-state rate (pre-n=2)
		P_HI = P_HI + gph_balmer_HI
		heat = heat + heat_balmer
	endif


   !----------------------------------!
      
   ! Evaluate cooling rates and recombination/collisional 
   ! 	ionization rates
      
    call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,           &
			   	   	rchiiB,rcheiiB,rcheiiiB, rec_m,             &
			    	a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
			    	cool)

	! Phase-3a: capture the ground-state H proton-budget coefficients on
	! every pass (the converged pass is the one read out by write_excited_H).
	if (use_excited_H) then
		cion_HI  = a_ion_HI          ! collisional ionization [cm^3 s^-1]
		arec_HII = rchiiB            ! case-B recombination  [cm^3 s^-1]
	endif

	! Charge-exchange rate coefficients are evaluated per cell below
	! (cx_set_cell) before each metal ionization solve.

	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten
	endif
	
   !----------------------------------!
      
   ! Ionization equilibrium system solution	

	if (.not.thereis_He) then ! If no helium

		do j = N+Ng,1-Ng,-1
		
			! Ionization equilibrium system setup
			params(1) = P_HI(j)
			params(2) = rchiiB(j)
			params(3) = nh(j)
			params(4) = a_ion_HI(j)

			 ! Initial guess
			if (count.le.0) then
				if(r(j).le.(1.5))then
					sys_x(1) = r(j)-0.5
				else
					sys_x(1) = 1.0
				endif     
			else

				sys_x(1) = nhii(j)/nh(j)

			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			call solve_ieq(ion_system_H,jac_system_H,N_eq,sys_x,    &
				      params,tol,wa,lwa,usednt)

			! Extract solution profiles
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)

		enddo
		
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0

	else  ! If there is helium

		! Metal unknowns start at sys_x(4) normally, but shift to sys_x(5)
		! when the He triplet occupies sys_x(4) (merged HeITR+metals system).
		mbase = 4
		if (thereis_HeITR) mbase = 5

		do j = N+Ng,1-Ng,-1

			! System coefficients
			params(1)  = P_HI(j)
			params(2)  = P_HeI(j)
			params(3)  = P_HeII(j)
			params(4)  = rchiiB(j)
			params(5)  = rcheiiB(j)
			params(6)  = rcheiiiB(j)
			params(7)  = nh(j)
			params(8)  = nhe(j)
			params(9)  = a_ion_HI(j) 
			params(10) = a_ion_HeI(j) 
			params(11) = a_ion_HeII(j)  
			
			! Add more if HeITR is present
			if (thereis_HeITR) then
				params(12) = rcheiTR(j)
				params(13) = A31
				params(14) = P_HeITR(j)
				params(15) = q13(j)
				params(16) = q31a(j)
				params(17) = q31b(j)
				params(18) = Q31
			endif

			! Per-element metal coefficients are handed to
			! ion_system_HeH_metals via set_metal_coeffs; the charge-
			! exchange rate coefficients are stored for this cell by
			! cx_set_cell (used inside the residual by cx_add_to_fvec).
			if (thereis_metals) then
				call cx_set_cell(T_K(j))

				! Build per-element metal coefficients in canonical order
				! from the 2D rate arrays (col i0 = neutral, i0+1 = +,
				! i0+2 = ++) and store them for the residual.
				do im = 1,n_melem
					i0  = melem_i0(im)
					top = melem_top(im)
					meg_top(im)  = top
					meg_ntot(im) = nm_tot(j,im)
					meg_g0(im)   = P_m(j,i0)
					meg_b0(im)   = aion_m(j,i0)
					meg_a1(im)   = rec_m(j,i0+1)
					if (top .ge. 2) then
						meg_g1(im) = P_m(j,i0+1)
						meg_b1(im) = aion_m(j,i0+1)
						meg_a2(im) = rec_m(j,i0+2)
					else
						! No second ionization stage for this element.
						meg_g1(im) = 0.0d0
						meg_b1(im) = 0.0d0
						meg_a2(im) = 0.0d0
					endif
				enddo
				call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1, &
				                    meg_b0, meg_b1, meg_a1, meg_a2,     &
				                    meg_top)
			endif

			! Initial guess
			if (count .eq. 0) then
				if (j .eq. N+Ng) then
					sys_x(1) = 1.0
					sys_x(2) = 1.0
					sys_x(3) = 1.0
					if (thereis_HeITR) sys_x(4) = 0.01
					if (thereis_metals) then
						! Each metal starts fully singly ionized.
						do im = 1,n_melem
							sys_x(mbase + 2*(im-1))   = 1.0  ! X+  frac
							sys_x(mbase+1 + 2*(im-1)) = 0.0  ! X++ frac
						enddo
					endif
				else
					sys_x(1) = nhii(j+1)/nh(j+1)
					sys_x(2) = nheii(j+1)/nhe(j+1)
					sys_x(3) = nheiii(j+1)/nhe(j+1)
					if (thereis_HeITR) &
						sys_x(4) = nheiTR(j+1)/nhe(j+1)
					if (thereis_metals) then
						do im = 1,n_melem
							i0 = melem_i0(im)
							sys_x(mbase+2*(im-1)) = nm(j+1,i0+1)        &
							                  /max(nm_tot(j+1,im),1.0d-30)
							if (melem_top(im) .ge. 2) then
								sys_x(mbase+1+2*(im-1)) = nm(j+1,i0+2)        &
								                  /max(nm_tot(j+1,im),1.0d-30)
							else
								sys_x(mbase+1+2*(im-1)) = 0.0d0
							endif
						enddo
					endif
				endif
			else
				sys_x(1) = nhii(j)/nh(j)
				sys_x(2) = nheii(j)/nhe(j)
				sys_x(3) = nheiii(j)/nhe(j)
				if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j)
				if (thereis_metals) then
					do im = 1,n_melem
						i0 = melem_i0(im)
						sys_x(mbase+2*(im-1)) = nm(j,i0+1)              &
						                  /max(nm_tot(j,im),1.0d-30)
						if (melem_top(im) .ge. 2) then
							sys_x(mbase+1+2*(im-1)) = nm(j,i0+2)              &
							                  /max(nm_tot(j,im),1.0d-30)
						else
							sys_x(mbase+1+2*(im-1)) = 0.0d0
						endif
					enddo
				endif
			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			! The He metastable-triplet systems keep the MINPACK solve (no
			! analytic Jacobian written for the triplet kinetics).
			if (thereis_HeITR .and. thereis_metals) then
				! Merged He-triplet + metals: triplet at sys_x(4), metals at
				! sys_x(5..). The per-cell metal coefficients (set_metal_coeffs)
				! and charge-exchange rates (cx_set_cell) were already loaded
				! above in the thereis_metals block; here we only point charge
				! exchange at the shifted metal rows for the merged solve.
				cx_metal_base = 5
				call hybrd1(ion_system_HeH_TR_metals,N_eq,sys_x,sys_sol,  &
						    tol,info,wa,lwa,params)
				cx_metal_base = 4
			else if (thereis_HeITR) then
				call hybrd1(ion_system_HeH_TR,N_eq,sys_x,sys_sol,   &
						    tol,info,wa,lwa,params)
			else if (thereis_metals) then
				call solve_ieq(ion_system_HeH_metals,jac_system_HeH_metals,  &
				            N_eq,sys_x,params,tol,wa,lwa,usednt)
			else
				call solve_ieq(ion_system_HeH,jac_system_HeH,       &
                            N_eq,sys_x,params,tol,wa,lwa,usednt)
			endif

			! Extract solution profiles
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)
			nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
			nheii(j)  = nhe(j)*sys_x(2)
			nheiii(j) = nhe(j)*sys_x(3)
			if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(4)
			if (thereis_metals) then
				do im = 1,n_melem
					i0 = melem_i0(im)
					nm(j,i0+1) = nm_tot(j,im)*sys_x(mbase+2*(im-1))
					if (melem_top(im) .ge. 2) then
						nm(j,i0+2) = nm_tot(j,im)*sys_x(mbase+1+2*(im-1))
						nm(j,i0)   = nm_tot(j,im)                            &
						           *(1.0 - sys_x(mbase+2*(im-1)) - sys_x(mbase+1+2*(im-1)))
					else
						! Two-stage element: neutral = total - singly ionized.
						nm(j,i0)   = nm_tot(j,im)*(1.0 - sys_x(mbase+2*(im-1)))
					endif
				enddo
			endif

		enddo
	
	endif

	
	! Density with atomic numbers (nm adds the metal mass under the
	! eos_metals policy)
   call calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_out,nm)

   ! Abundancies profiles
   f_sp_out(:,1) = nhi/n_out
   f_sp_out(:,2) = nhii/n_out
   f_sp_out(:,3) = nhei/n_out
   f_sp_out(:,4) = nheii/n_out
   f_sp_out(:,5) = nheiii/n_out
   f_sp_out(:,6) = nheiTR/n_out
   ! Metal abundances in canonical order (col mion_fsp(im) of f_sp_out).
   do im = 1,n_mion
      f_sp_out(:,mion_fsp(im)) = nm(:,im)/n_out
   enddo

   ! Adimensional number density profile
	n_out = n_out/n0

   ! Adimensional heating and cooling rates
   heat_out = heat/q0
   cool_out = cool/q0
      
   ! Adjust value of pressure boundary condition (the base electron
   ! density in units of n0; with eos_metals the metal electrons are
   ! included, consistently with calc_ne)
   dp_bc = (nhii(1-Ng) + nheii(1-Ng) + 2.0*nheiii(1-Ng))/n0
   if (eos_include_metals .and. thereis_metals) then
      do im = 1,n_mion
         if (mion_stage(im) .gt. 0)                                     &
            dp_bc = dp_bc + dble(mion_stage(im))*nm(1-Ng,im)/n0
      enddo
   endif

	! End of subroutine 
	end subroutine ioniz_eq
	
	! End of module
	end module ionization_equilibrium         

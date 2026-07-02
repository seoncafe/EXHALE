   module utils_ion_eq

   use global_parameters
   use species_table, only: n_mion, n_mphot, n_melem,                 &
                            mion_isphot, mion_iphot, mion_ethr,       &
                            mion_z2, mion_elem, mion_iscool,          &
                            mion_stage, mion_name, mion_fsp,          &
                            melem_Z, melem_top, im_FeII
   use utils
   use Cooling_Coefficients      ! Various functions for cooling coefficients
   use omp_lib                   ! OMP libraries
	
	! Move here photoionization and photoheating
	! Two subroutines for H and He+H
	
   implicit none
	
	contains
	
	! ------------------------------------------------------------- !

	subroutine PH_heat_H(nhi,P_HI,heat,q)
	! Computes photoionization rates and heating rates for
	!	an atmosphere composed of H and He
	
	integer :: i,j
	
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi    

	! Dummy zero 
	real*8, dimension(1-Ng:N+Ng), parameter :: nhei   = 0.0, nheii  = 0.0
	real*8, dimension(1-Ng:N+Ng), parameter :: nheiii = 0.0, nheiTR = 0.0
	real*8, dimension(1-Ng:N+Ng) :: N15, N2, NTR

   real*8 :: dr	                ! Grid spacing
   real*8 :: PIR_1		            ! Photoionization rates
   real*8 :: Hea_1 		        ! Heating rates  
	real*8 :: q_abs           	! Absorbed energy
	
	! Integral variables
	real*8, dimension(Nl) :: tauE			
	real*8, dimension(Nl) :: int_f,int_1,int_q,int_H

	! Column densities                                 
   real*8, dimension(1-Ng:N+Ng) ::  N1 
	
   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  P_HI 
      
   ! Heating efficiency 
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q
      
   ! Heating rate
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat
      
	!----------------------------------!
	
	! Evaluate the column density
	call calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)

   ! Evaluate photoionization rates and photoheating rates
	do j = 1-Ng,N+Ng   

		! Initialization of integrands
		Hea_1  = 0.0
		PIR_1  = 0.0
		q_abs  = 0.0

		tauE = s_hi*N1(j)*1.0e-18	

		! Initial integrands
		int_f = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		int_H = int_f*(1.0-e_th_HI/e_v)*s_hi*nhi(j)
		int_1 = int_f*s_hi/e_v
		int_q = int_f*s_hi*nhi(j)

		! Value of integrals
		Hea_1 = sum(int_H*de_v)	
		PIR_1 = sum(int_1*de_v)	
		q_abs = sum(int_q*de_v)
		
		! Multiply for the dimensional coefficient 
		heat(j)   = Hea_1*1.0e-18	
		P_HI(j)   = PIR_1*1.0e-18*erg2eV	
		! Guard: q_abs (absorbed-energy normalization) can be 0 in a fully
		! transparent/unilluminated cell; avoid 0/0 -> NaN in the efficiency.
		q(j)      = Hea_1/max(q_abs, 1.0d-99)
	
	enddo
	
	end subroutine PH_heat_H

	! ------------------------------------------------------------- !
	
	subroutine PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm,             &
				     P_HI,P_HeI,P_HeII,P_HeITR, P_m,       &
				     heat,q)
	! Computes photoionization rates and heating rates for an
	!	atmosphere composed of H, He, and (optionally) metals.
	! Metal ion densities arrive as nm(:,1:n_mion) in canonical
	! species_table order; per-ion photoionization rates leave as
	! P_m(:,1:n_mion). Only photo-ionizable metal ions (mion_isphot)
	! contribute to opacity, photoheating, absorbed energy, and have a
	! nonzero P_m; top-stage ions (C III, N III, O III, Mg III) are inert.

	integer :: i,j,k

	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhei,nheii
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm
	real*8, dimension(1-Ng:N+Ng) :: nheiS
   real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2,NTR
   real*8, dimension(1-Ng:N+Ng,n_mphot) :: Nm_col

   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR     ! Photoionization rates (H/He)
   real*8 :: Hea_1                         ! Heating rate
   real*8 :: q_abs                         ! Absorbed energy
   real*8 :: Pm_loc(n_mion)                ! Per-cell metal photoion. rates

	! Integral variables
	real*8, dimension(Nl) :: tauE,tau_m
	real*8, dimension(Nl) :: int_f,int_1,int_15,int_2,int_TR
	real*8, dimension(Nl) :: int_m
	real*8, dimension(Nl) :: int_q,int_H,acc_H,acc_q

   ! Photo ionization rates
	real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HI
   real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) ::  P_m

   ! Heating efficiency
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q

   ! Heating rate
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat

	!----------------------------------!

	! Use nheiS as variable
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR

	! Evaluate the column density
	call calc_column_dens(nhi,nheiS,nheii,nheiTR,N1,N15,N2,NTR)
	call calc_column_dens_metals(nm, Nm_col)

	! Metal-free P_m entries (top-stage ions) stay zero
	P_m = 0.0

    !----------------------------------!
	!$OMP PARALLEL DO &
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,P_HeITR,P_m,heat,q )                        &
	!$OMP PRIVATE ( Hea_1,PIR_1,PIR_15,PIR_2,PIR_TR,Pm_loc,                      &
	!$OMP           int_1,int_15,int_2,int_TR,int_m,                            &
	!$OMP           int_f,int_H,int_q,acc_H,acc_q,tauE,tau_m,q_abs,i,k,j)

	do j = 1-Ng,N+Ng

		Hea_1   = 0.0
      PIR_1   = 0.0
      PIR_15  = 0.0
      PIR_2   = 0.0
      PIR_TR  = 0.0
      q_abs   = 0.0

		! Calculate optical depth. Accumulate the metal block separately
		! in iphot order before applying the 1e-18 factor, matching the
		! original (sum)*1e-18 association exactly.
		tauE = (s_hi*N1(j) + s_hei*N15(j) + s_heii*N2(j))*1.0e-18
		if (thereis_HeITR) tauE = tauE + s_heiTR*NTR(j)*1.0e-18
		tau_m = 0.0
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			tau_m = tau_m + sigma_tab(:,k)*Nm_col(j,k)
		enddo
		tauE = tauE + tau_m*1.0e-18

		! Calculate photoionization integrals
		int_f  = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		int_1  = int_f*s_hi/e_v
		int_15 = int_f*s_hei/e_v
		int_2  = int_f*s_heii/e_v
		if (thereis_HeITR) int_TR =  int_f*s_heiTR/e_v

		! Photoheating integral. Accumulate the inner bracket in the
		! original order (H, He, then metals in iphot order) and apply
		! the int_f factor once, preserving the original association.
		acc_H = (1.0-e_th_HI  /e_v)*s_hi  *nhi  (j)
		acc_H = acc_H + (1.0-e_th_HeI /e_v)*s_hei *nheiS(j)
		acc_H = acc_H + (1.0-e_th_HeII/e_v)*s_heii*nheii(j)
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			acc_H = acc_H + (1.0-mion_ethr(i)/e_v)*sigma_tab(:,k)*nm(j,i)
		enddo
		int_H = int_f*acc_H

		! Absorbed energy integral (same ordering as int_H)
		acc_q = s_hi *nhi  (j)
		acc_q = acc_q + s_hei *nheiS(j)
		acc_q = acc_q + s_heii*nheii(j)
		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			acc_q = acc_q + sigma_tab(:,k)*nm(j,i)
		enddo
		int_q = int_f*acc_q

		! Midpoint rule for integrals
		PIR_1   = sum(int_1  *de_v)
		PIR_15  = sum(int_15 *de_v)
		PIR_2   = sum(int_2  *de_v)
		if(thereis_HeITR) PIR_TR = sum(int_TR*de_v)
		do i = 1,n_mion
			if (.not. mion_isphot(i)) then
				Pm_loc(i) = 0.0
				cycle
			endif
			k = mion_iphot(i)
			int_m = int_f*sigma_tab(:,k)/e_v
			Pm_loc(i) = sum(int_m*de_v)*1.0e-18*erg2eV
		enddo
		Hea_1   = sum(int_H  *de_v)
		q_abs   = sum(int_q  *de_v)

		!$OMP CRITICAL
		! Save into vector
    	P_HI(j)    = PIR_1  *1.0e-18*erg2eV
    	P_HeI(j)   = PIR_15 *1.0e-18*erg2eV
		P_HeII(j)  = PIR_2  *1.0e-18*erg2eV
		P_HeITR(j) = PIR_TR *1.0e-18*erg2eV
		P_m(j,:)   = Pm_loc(:)
		heat(j)    = Hea_1*1.0e-18
		! Guard against q_abs = 0 (see PH_heat_H).
		q(j)       = Hea_1/max(q_abs, 1.0d-99)
		!$OMP END CRITICAL

	enddo
	!$OMP END PARALLEL DO

	! End of subroutine
	end subroutine PH_heat_HHe
	
	!----------------------------------!
	
	subroutine eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,        &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,             &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,      &
				   cool, cool_chan)

	! Evaluate the cooling rate contributions to energy and
	!	rate equations. Includes H, He, and metal channels. Metals are
	!	driven from species_table metadata (canonical ion order); nm,
	!	rec_m and aion_m carry the per-ion densities and rates, so the
	!	argument list no longer grows when a metal is added.

	integer :: i,j,e

	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,           &
	                                            nhei,nheii,nheiii
	! Metal ion densities (canonical species_table order)
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm

	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

   real*8, dimension(1-Ng:N+Ng) :: brem,coex,coio,reco  ! Cooling rates
   real*8, dimension(1-Ng:N+Ng) :: cool_M               ! Metal cooling
   real*8, dimension(1-Ng:N+Ng) :: brem_acc,coolm_acc   ! sum accumulators
   real*8, dimension(1-Ng:N+Ng) :: metal_col            ! dispatcher scratch
   real*8, dimension(1-Ng:N+Ng,n_mion)  :: c_metal      ! metal line-cool coeffs
   real*8, dimension(1-Ng:N+Ng,n_melem) :: GF_elem      ! per-element Gaunt fac.
   real*8, dimension(1-Ng:N+Ng) :: tau_eff,beta_esc     ! beta escape prob.
   real*8, dimension(1-Ng:N+Ng) :: kappa_loc,dr_cm
	real*8, dimension(1-Ng:N+Ng) :: ne		  			 ! Electron number density
	real*8, dimension(1-Ng:N+Ng) :: GF_H,GF_He			 ! Gaunt factors
	real*8 :: Cdex_OI,Cdex_CII                           ! 2-level collis. de-exc.
	! Named bridges for the (verbatim) two-level cooling branch
	real*8, dimension(1-Ng:N+Ng) :: nci,ncii,noi,noii,nmgi,nmgii
	real*8, dimension(1-Ng:N+Ng) :: c_CI,c_OII,c_MgI,c_MgII

   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB
   ! Per-ion metal recombination rates (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: rec_m

	! Recombination cooling coefficients
	real*8, dimension(1-Ng:N+Ng) :: coeff_rec_cool_HII,  &
									coeff_rec_cool_HeII, &
									coeff_rec_cool_HeIII

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII
   ! Per-ion metal collisional ionization rates (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: aion_m

	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_rate_HI,    &
	 								coeff_coex_rate_HeI,   &
									coeff_coex_rate_HeII
									 
	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  cool

	! Optional per-channel cooling breakdown (cgs erg cm^-3 s^-1, same
	! units as `cool`). Columns 1-4 = H/He recombination, collisional
	! ionization, collisional excitation, bremsstrahlung (the last incl.
	! metal-ion charges); columns 4+i = metal ion i line cooling (0 for
	! non-coolant ions). This is an exact decomposition of `cool` in the
	! default (.not.use_2lev_cool) branch; in the two-level branch the
	! per-ion metal terms are the resonance-line approximation and need
	! not sum to cool_M.
	real*8, dimension(1-Ng:N+Ng,4+n_mion),intent(out),optional :: cool_chan

   ! Free electron density (incl. metal electrons under eos_metals)
   call calc_ne(nhii,nheii,nheiii,ne,nm)
	

	!-- Recombination --!

	! Rate coefficients
	call rec_HII_B(T_K,rchiiB)      ! HII
	call rec_HeII_B(T_K,rcheiiB)    ! HeII
	call rec_HeIII_B(T_K,rcheiiiB)  ! HeIII
	
	! Cooling rate coefficients
	call rec_cool_HII(T_K,coeff_rec_cool_HII)
	call rec_cool_HeII(T_K,coeff_rec_cool_HeII)
	call rec_cool_HeIII(T_K,coeff_rec_cool_HeIII)

	! Cooling rate
	reco  = coeff_rec_cool_HII*nhii     & ! HII
		   + coeff_rec_cool_HeII*nheii   & ! HeII
		   + coeff_rec_cool_HeIII*nheiii   ! HeIII

	!-- Collisional ionization --!
	
	! Rate coefficients
	call ion_coeff_HI(T_K,a_ion_HI)      ! HI
	call ion_coeff_HeI(T_K,a_ion_HeI)    ! HeI
	call ion_coeff_HeII(T_K,a_ion_HeII)  ! HeII
	
	! Cooling rate
	coio =  2.179e-11*a_ion_HI*nhi  	 & ! HI
		  + 3.940e-11*a_ion_HeI*nhei 	 & ! HeI
		  + 8.715e-11*a_ion_HeII*nheii   ! HeII
	
	!-- Bremsstrahlung --!

	! Gaunt factors (H, He explicit; one per metal element from its Z)
	call GF(T_K,ih,GF_H)
	call GF(T_K,ihe,GF_He)
	do e = 1,n_melem
		call GF(T_K, dble(melem_Z(e)), GF_elem(:,e))
	enddo

	! Cooling rate (extended with metal ion charges). Accumulate the
	! charge^2-weighted sum in the original H, He, then canonical-metal
	! order and apply the 1.426e-27*sqrt(T_K) prefactor once, so the
	! result is bit-identical to the explicit expression. Neutral metals
	! carry z2=0 and contribute an exact +0 term.
	brem_acc = ih**2.0*GF_H*nhii                          ! HII
	brem_acc = brem_acc + ihe**2.0*GF_He*(nheii + nheiii) ! He
	do i = 1,n_mion
		brem_acc = brem_acc                                          &
		         + mion_z2(i)*GF_elem(:,mion_elem(i))*nm(:,i)
	enddo
	brem = 1.426e-27*sqrt(T_K)*brem_acc

	!-- Collisional excitation --!

	! Rate coefficients
	call coex_rate_HI(T_K,coeff_coex_rate_HI) 		! HI
	call coex_rate_HeI(T_K,coeff_coex_rate_HeI)   	! HeI
	call coex_rate_HeII(T_K,coeff_coex_rate_HeII)  	! HeII

	! Cooling rate
	coex = coeff_coex_rate_HI*nhi       &    ! HI
		  + coeff_coex_rate_HeI*nhei     &    ! HeI
		  + coeff_coex_rate_HeII*nheii        ! HeII

	!-- Metal recombination + collisional ionization rates --!
	! (rates for ionization equilibrium; not part of cool here.)
	! Filled per canonical ion via the metadata dispatchers, which call
	! the same per-ion routines as before; ions with no entry (inert top
	! stage, neutral non-recombiner) return 0.
	do i = 1,n_mion
		call rec_coeff_by_ion(i,T_K,metal_col)
		rec_m(:,i)  = metal_col
		call ion_coeff_by_ion(i,T_K,metal_col)
		aion_m(:,i) = metal_col
	enddo

	!-- Metal radiative cooling (forbidden/fine-structure lines) --!
	do i = 1,n_mion
		call cool_coeff_by_ion(i,T_K,metal_col)
		c_metal(:,i) = metal_col
	enddo

	! Density-dependent override for Fe II line cooling. The 1-D coronal
	! cool_coeff_metal('FeII',...) above overestimates cooling at the dense
	! base by ~1e4x because the forbidden a6D fine-structure / metastable
	! lines (n_crit ~ 1e4-1e7 cm^-3) are collisionally saturated there
	! (n_e >> n_crit). Replace c_metal(:,FeII) with the multilevel
	! statistical-equilibrium coefficient Lambda_eff(T,ne) = (sum_u n_u A_ul
	! dE_ul)/ne. In the assembly cool_M = beta_esc*ne*sum_i nm(:,i)*c_metal,
	! the ne cancels the 1/ne in Lambda_eff, leaving the correct LTE-saturated
	! per-ion cooling (collider-independent, so the electron-only SE solve is
	! exact in this limit). At low ne it reduces to the coronal rate.
	call cool_FeII_ne(T_K, ne, metal_col)
	c_metal(:,im_FeII) = metal_col

	! Density-dependent override for the [C II] 158um / [O I] 63um
	! ground-term fine-structure floors (CHIANTI mode only; the legacy
	! AIOLOS fits keep their own constant floors). Same Lambda_eff =
	! W_FS/ne + remainder convention as Fe II above; the two-level
	! solution saturates the floor (n_crit,e([C II]) ~ 20 cm^-3!) and
	! adds the H-collision excitation channel the electron-only coronal
	! curve misses. See cool_CII_ne_func / cooling_data/
	! fit_fs_saturation.py.
	if (cno_chianti) then
		call cool_CII_ne(T_K, ne, nhi, metal_col)
		c_metal(:,im_CII) = metal_col
		call cool_OI_ne(T_K, ne, nhi, metal_col)
		c_metal(:,im_OI) = metal_col
	endif

	! AIOLOS-style β escape probability for resonance line trapping.
	! tau_eff = local opacity (lowest XUV band) * cell width [cm].
	! s_*(1) are in 1e-18 cm^2, so multiply by 1e-18.
	! The trailing 1.0e8 factor mirrors AIOLOS chemistry.cpp:1006:
	! a boost converting gray opacity to metal line-center opacity.
	! See feedback memory: "To Be Checked/AIOLOS tuning?".
	do j = 1-Ng,N+Ng
		dr_cm(j)     = dr_j(j)*R0
		kappa_loc(j) = (s_hi(1)*nhi(j) + s_hei(1)*nhei(j) +           &
		                s_heii(1)*nheii(j))*1.0e-18
		!To Be Checked/AIOLOS tuning?
		!tau_eff(j) = kappa_loc(j) * dr_cm(j) * 1.0e8
		tau_eff(j) = kappa_loc(j) * dr_cm(j)
		if (tau_eff(j) .lt. 7.0) then
			beta_esc(j) = (1.0 - exp(-2.34*tau_eff(j)))/             &
			              (4.68*max(tau_eff(j),1.0d-30))
		else
			beta_esc(j) = 1.0/(4.0*tau_eff(j)*                       &
			              sqrt(log(tau_eff(j)/sqrt(pi))))
		endif
	enddo
	beta_esc(1-Ng) = 0.0

	! AIOLOS multiplies the escape probability by a 1e8 factor
	! (chemistry.cpp:1006), which drives tau_eff up to ~1e6 so that
	! beta_esc -> ~0 and the metal lines are treated as essentially
	! fully trapped (metal-line cooling switched off).
	! Here we instead assume 100% escape (optically-thin limit):
	! beta_esc = 1, so the full metal-line cooling is applied,
	! matching ATES_extended's coronal treatment.
	beta_esc = 1.0d0

	if (use_2lev_cool) then
		! Bridge the metadata arrays to the named scalars used below.
		nci   = nm(:,1)
		ncii  = nm(:,2)
		noi   = nm(:,4)
		noii  = nm(:,5)
		nmgi  = nm(:,10)
		nmgii = nm(:,11)
		c_CI   = c_metal(:,1)
		c_OII  = c_metal(:,5)
		c_MgI  = c_metal(:,10)
		c_MgII = c_metal(:,11)
		! Two-level fine-structure cooling for the dominant coolants
		! [O I] 63um and [C II] 158um (critical-density saturation +
		! H-atom collisions). C I and O II keep the Black-1981 fit; the
		! exponential ("forbidden") part of the C II / O I fit is retained
		! on top of the two-level ground term (matching ATES_extended).
		! N has no line cooling.
		do j = 1-Ng,N+Ng
			! Collisional de-excitation rates [s^-1]
			Cdex_OI  = nhi(j)*4.2d-11*(T_K(j)/100.0d0)**0.67          ! H
			Cdex_CII = ne(j) *8.7d-8 *(T_K(j)/2000.0d0)**(-0.37)      & ! e
			         + nhi(j)*4.0d-11                                   ! H
			cool_M(j) = beta_esc(j)*(                                       &
			    ne(j)*nci(j)*c_CI(j)                                       &
			  + noi(j) *( lambda_2level(8.91d-5,227.7d0,0.6d0,Cdex_OI ,T_K(j)) &
			              + ne(j)*1.1d-20*exp(-30162.0d0/T_K(j))           &
			                     *(1.0d0+(T_K(j)/0.75d4)**0.5) )           &
			  + ncii(j)*( lambda_2level(2.29d-6,91.21d0,2.0d0,Cdex_CII,T_K(j)) &
			              + ne(j)*3.1d-20*exp(-45162.0d0/T_K(j))           &
			                     *(1.0d0+(T_K(j)/0.75d4)**1.5) )           &
			  + ne(j)*noii(j)*c_OII(j)                                  &
			  + ne(j)*nmgi(j)*c_MgI(j) + ne(j)*nmgii(j)*c_MgII(j)                          &
                  + ne(j)*nm(j,17)*c_metal(j,17)               &
                  + ne(j)*nm(j,19)*c_metal(j,19)               &
                  + ne(j)*nm(j,26)*c_metal(j,26) )
		enddo
	else
		! Sum the line-cooling metal ions (mion_iscool) in canonical
		! order, then apply the (beta_esc*ne) prefactor once. This is
		! bit-identical to the explicit eight-term expression: the
		! iscool ions, in canonical order, are exactly
		! CI,CII,OI,OII,NI,NII,MgI,MgII.
		coolm_acc = 0.0d0
		do i = 1,n_mion
			if (.not. mion_iscool(i)) cycle
			coolm_acc = coolm_acc + nm(:,i)*c_metal(:,i)
		enddo
		cool_M = beta_esc * ne * coolm_acc
	endif

	! Total cooling rate
	cool = ne*(brem + coex + reco + coio) + cool_M

	! Per-channel breakdown for the Phase-2 diagnostic (Huang Fig. 10).
	! Read straight from the arrays already computed above, so the sum of
	! all channels reproduces `cool` exactly in the default branch.
	if (present(cool_chan)) then
		cool_chan(:,1) = ne*reco
		cool_chan(:,2) = ne*coio
		cool_chan(:,3) = ne*coex
		cool_chan(:,4) = ne*brem
		do i = 1,n_mion
			if (mion_iscool(i)) then
				cool_chan(:,4+i) = beta_esc*ne*nm(:,i)*c_metal(:,i)
			else
				cool_chan(:,4+i) = 0.0d0
			endif
		enddo
	endif

	! End of subroutine
	end subroutine eval_cool

	! ------------------------------------------------------------- !

	subroutine write_cool_breakdown_eq(T_in,n_in,f_sp_in)
	! Phase-2 diagnostic. Dump the per-channel radiative cooling rate vs
	! radius for the converged equilibrium state, reusing eval_cool's exact
	! coefficients (no offline re-derivation). Columns: H/He recombination,
	! collisional ionization, collisional excitation, bremsstrahlung, then
	! one column per metal ion line-cooling channel (canonical species_table
	! order). All in cgs erg cm^-3 s^-1; the channel sum reproduces the
	! Hydro_ioniz.txt `cool` column. The printed max relative residual is the
	! internal consistency check. Lets the user identify the dominant coolant
	! in 1.15 <~ r/Rp <~ 1.4 against Huang et al. (2023) Fig. 10.

	integer :: j,i,im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) :: n_dim,T_K,ne
	real*8, dimension(1-Ng:N+Ng) :: nhi,nhii,nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	real*8, dimension(1-Ng:N+Ng) :: cool,csum,rel
	real*8, dimension(1-Ng:N+Ng,4+n_mion) :: chan
	! Throwaway eval_cool rate outputs (not needed for the dump)
	real*8, dimension(1-Ng:N+Ng) :: rchiiB,rcheiiB,rcheiiiB
	real*8, dimension(1-Ng:N+Ng) :: a_ion_HI,a_ion_HeI,a_ion_HeII
	real*8, dimension(1-Ng:N+Ng,n_mion) :: rec_m,aion_m
	real*8 :: maxrel

	! Dimensionalize exactly as ioniz_eq does
	n_dim = n_in*n0
	T_K   = T_in*T0
	nhi   = f_sp_in(:,1)*n_dim
	nhii  = f_sp_in(:,2)*n_dim
	if (thereis_He) then
		nhei   = f_sp_in(:,3)*n_dim
		nheii  = f_sp_in(:,4)*n_dim
		nheiii = f_sp_in(:,5)*n_dim
	else
		nhei = 0.0d0; nheii = 0.0d0; nheiii = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	call calc_ne(nhii,nheii,nheiii,ne,nm)

	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,            &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,                &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
				   cool, cool_chan = chan)

	! Internal consistency: channel sum vs total cool (default branch -> ~eps)
	csum = 0.0d0
	do i = 1,4+n_mion
		csum = csum + chan(:,i)
	enddo
	rel    = abs(csum - cool)/max(abs(cool),1.0d-99)
	maxrel = maxval(rel(1:N))
	write(*,'(a,es9.2)')                                                  &
		' (write_cool_breakdown_eq) max |sum(channels)/cool - 1| = ', maxrel

	open(unit = 71, file = './output/Cooling_breakdown.txt')
	write(71,'(a)') '# Per-channel radiative cooling rate [cgs erg cm^-3 s^-1] vs radius.'
	write(71,'(a)') '# Channel sum reproduces the Hydro_ioniz.txt cool column.'
	write(71,'(a)') '# col1 r/Rp  col2 T[K]  col3 ne  col4 cool_total  col5 reco'  &
	             // '  col6 coio  col7 coex  col8 brem  then one col per metal ion:'
	write(71,'(a)',advance='no') '#   metal-ion columns (canonical order):'
	do i = 1,n_mion
		write(71,'(1x,a)',advance='no') trim(mion_name(i))
	enddo
	write(71,*)
	do j = 1-Ng,N+Ng
		write(71,*) r(j), T_K(j), ne(j), cool(j),                        &
		            chan(j,1), chan(j,2), chan(j,3), chan(j,4),           &
		            (chan(j,4+i), i = 1,n_mion)
	enddo
	close(71)

	end subroutine write_cool_breakdown_eq

	!----------------------------------!

	! Various coefficients for HeI triplet chemistry
	subroutine HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31a,q31b,Q31)
	
	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K

	real*8, intent(out) :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiTR,rcheii,   &
								   q13,q31a,q31b,Q31

	call rec_HeII_23S(T_K,rcheiTR)
	call rec_HeII_11S(T_K,rcheii)
	call coex_HeI_1S_23S(T_K,q13)
	call coex_HeI_23S_21S(T_K,q31a)
	call coex_HeI_23S_21P(T_K,q31b)
	A31 = 1.272e-4

	! Penning ionization He(2^3S)+H: temperature-dependent Taylor et al. (2025)
	! rate, used unconditionally (the legacy temperature-independent 5e-10
	! constant is no longer selectable).
	call penning_HeI_23S(T_K,Q31)

   ! End of subroutine
	end subroutine HeITR_coeffs
	
	! End of module
	end module utils_ion_eq

   module utils_ion_eq

   use global_parameters
   use species_table, only: n_mion, n_mphot, n_melem,                 &
                            mion_isphot, mion_iphot, mion_ethr,       &
                            mion_z2, mion_elem, mion_iscool,          &
                            mion_stage, mion_name, mion_fsp,          &
                            melem_Z, melem_top, im_FeII
   use utils
   use Cooling_Coefficients      ! Various functions for cooling coefficients
   use Cross_sections, only: sigma, sigma_HeI  ! sigma_H(E,Z), sigma_HeI(E)
   use omp_lib                   ! OMP libraries
	
	! Move here photoionization and photoheating
	! Two subroutines for H and He+H
	
   implicit none
	
	contains
	
	! ------------------------------------------------------------- !

	subroutine PH_heat_H(nhi, xion, P_HI,heat,q)
	! Computes photoionization rates and heating rates for
	!	an atmosphere composed of H and He

	integer :: i,j

	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng),intent(in) :: xion

	! SvS85 secondary-ionization scratch (H-only: no He I secondary channel).
	real*8, dimension(Nl) :: acc_secHI, fhv
	real*8 :: xj, fh, fiHI, R_secHI
	! SvS85 coupling applied only when enabled AND staged on (see EXHALE_main).
	logical :: sec_on

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
	
	sec_on = use_sec_ion .and. sec_ion_active

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
		! SvS85 secondary-ionization energy partition (scalars for this cell).
		if (sec_on) then
			xj   = min(max(xion(j), 0.0d0), 1.0d0)
			fh   = svs85_fheat(xj)
			fiHI = svs85_fion_HI(xj)
		else
			fh = 1.0d0; fiHI = 0.0d0
		endif
		! Heating fraction: fh above the E_sec_ion photoelectron threshold, 1
		! (full thermalization) below it. fhv = 1 when the coupling is off, so
		! the heating integrand is bit-identical to the legacy path.
		fhv = 1.0d0
		if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_HI + E_sec_ion)
		int_H = int_f*(1.0-e_th_HI/e_v)*fhv*s_hi*nhi(j)
		int_1 = int_f*s_hi/e_v
		int_q = int_f*s_hi*nhi(j)

		! Value of integrals
		Hea_1 = sum(int_H*de_v)	
		PIR_1 = sum(int_1*de_v)	
		q_abs = sum(int_q*de_v)
		
		! Multiply for the dimensional coefficient 
		heat(j)   = Hea_1*1.0e-18
		P_HI(j)   = PIR_1*1.0e-18*erg2eV
		! Add the H I secondary-ionization rate from fast photoelectrons.
		if (sec_on) then
			acc_secHI = int_f*s_hi*nhi(j)/e_v * &
			     merge(fiHI*(e_v-e_th_HI)/e_th_HI, 0.0d0, e_v > e_th_HI + E_sec_ion)
			R_secHI = sum(acc_secHI*de_v)*1.0e-18*erg2eV
			P_HI(j) = P_HI(j) + R_secHI/max(nhi(j), 1.0d-99)
		endif
		! Guard: q_abs (absorbed-energy normalization) can be 0 in a fully
		! transparent/unilluminated cell; avoid 0/0 -> NaN in the efficiency.
		q(j)      = Hea_1/max(q_abs, 1.0d-99)
	
	enddo
	
	end subroutine PH_heat_H

	! ------------------------------------------------------------- !
	
	subroutine PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,       &
				     P_HI,P_HeI,P_HeII,P_HeITR, P_m,       &
				     heat,q, nh2,P_H2, heat_chan)
	! Computes photoionization rates and heating rates for an
	!	atmosphere composed of H, He, and (optionally) metals.
	! Metal ion densities arrive as nm(:,1:n_mion) in canonical
	! species_table order; the photoionization rates for each ion leave as
	! P_m(:,1:n_mion). Only photo-ionizable metal ions (mion_isphot)
	! contribute to opacity, photoheating, absorbed energy, and have a
	! nonzero P_m; top-stage ions (C III, N III, O III, Mg III) are inert.

	integer :: i,j,k

	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhei,nheii
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion),intent(in) :: nm
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng),intent(in) :: xion
	! Optional H2 (molecular extension): adds the H2 opacity,
	! photoionization rate P_H2 [1/s] and photoelectric heating using the
	! Yan+1998 cross section s_h2 (threshold e_th_H2 = 15.4 eV).
	real*8, dimension(1-Ng:N+Ng), intent(in),  optional :: nh2
	real*8, dimension(1-Ng:N+Ng), intent(out), optional :: P_H2
	real*8, dimension(1-Ng:N+Ng) :: nheiS
   real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2,NTR
   real*8, dimension(1-Ng:N+Ng,n_mphot) :: Nm_col

   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR     ! Photoionization rates (H/He)
   real*8 :: PIR_H2                        ! H2 photoionization rate
   real*8, dimension(1-Ng:N+Ng) :: NH2col  ! H2 column density
   real*8, dimension(Nl) :: int_h2
   real*8 :: Hea_1                         ! Heating rate
   real*8 :: q_abs                         ! Absorbed energy
   real*8 :: Pm_loc(n_mion)                ! Metal photoion. rates in each cell

	! Integral variables
	real*8, dimension(Nl) :: tauE,tau_m
	real*8, dimension(Nl) :: int_f,int_1,int_15,int_2,int_TR
	real*8, dimension(Nl) :: int_m
	real*8, dimension(Nl) :: int_q,int_H,acc_H,acc_q
	! SvS85 secondary-ionization scratch.
	real*8, dimension(Nl) :: acc_secHI,acc_secHeI,fhv
	real*8 :: xj,fh,fiHI,fiHeI,R_secHI,R_secHeI
	! SvS85 coupling applied only when enabled AND staged on (see EXHALE_main).
	logical :: sec_on

   ! Photo ionization rates
	real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HI
   real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(out) ::  P_m

   ! Heating efficiency
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q

   ! Heating rate
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat

	! Optional per-absorber photoheating breakdown (cgs erg cm^-3 s^-1).
	! Columns: 1 H I, 2 He I (singlet ground), 3 He II, 4 He 2^3S,
	! 5 H2 (molecular, 0 when absent), 6 metals (sum over photo-ionizable
	! metal ions). Columns 1-6 sum to `heat` up to rounding; `heat` itself is
	! computed unchanged, so this is a diagnostic-only add-on.
	real*8, dimension(1-Ng:N+Ng,6),intent(out),optional :: heat_chan
	! Frequency-integrand accumulators for the per-absorber split (only used
	! when heat_chan is present).
	real*8, dimension(Nl) :: acc_HI,acc_HeI,acc_HeII,acc_HeTR,acc_H2,acc_mtl

	!----------------------------------!

	! Use nheiS as variable
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR

	! Evaluate the column density
	call calc_column_dens(nhi,nheiS,nheii,nheiTR,N1,N15,N2,NTR)
	call calc_column_dens_metals(nm, Nm_col)

	! H2 column (molecular)
	if (present(nh2)) call calc_column_dens_one(nh2, NH2col)

	! Metal-free P_m entries (top-stage ions) stay zero
	P_m = 0.0

	sec_on = use_sec_ion .and. sec_ion_active

    !----------------------------------!
	!$OMP PARALLEL DO &
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,P_HeITR,P_m,heat,q, sec_on )                &
	!$OMP PRIVATE ( Hea_1,PIR_1,PIR_15,PIR_2,PIR_TR,PIR_H2,int_h2,               &
	!$OMP           Pm_loc,                                                      &
	!$OMP           int_1,int_15,int_2,int_TR,int_m,                            &
	!$OMP           acc_secHI,acc_secHeI,fhv,xj,fh,fiHI,fiHeI,R_secHI,R_secHeI, &
	!$OMP           acc_HI,acc_HeI,acc_HeII,acc_HeTR,acc_H2,acc_mtl,            &
	!$OMP           int_f,int_H,int_q,acc_H,acc_q,tauE,tau_m,q_abs,i,k,j)

	do j = 1-Ng,N+Ng

		Hea_1   = 0.0
      PIR_1   = 0.0
      PIR_15  = 0.0
      PIR_2   = 0.0
      PIR_TR  = 0.0
      q_abs   = 0.0

		! Per-absorber photoheating accumulators default to zero so the He 2^3S
		! and H2 columns stay 0 in cells/runs where those absorbers are absent.
		if (present(heat_chan)) then
			acc_HeTR = 0.0d0
			acc_H2   = 0.0d0
			acc_mtl  = 0.0d0
		endif

		! Calculate optical depth. Accumulate the metal block separately
		! in iphot order before applying the 1e-18 factor, matching the
		! original (sum)*1e-18 association exactly.
		tauE = (s_hi*N1(j) + s_hei*N15(j) + s_heii*N2(j))*1.0e-18
		if (thereis_HeITR) tauE = tauE + s_heiTR*NTR(j)*1.0e-18
		if (present(nh2)) tauE = tauE + s_h2*NH2col(j)*1.0e-18
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
		if (present(nh2)) int_h2 = int_f*s_h2/e_v

		! SvS85 secondary-ionization energy partition (scalars for this cell).
		if (sec_on) then
			xj    = min(max(xion(j), 0.0d0), 1.0d0)
			fh    = svs85_fheat(xj)
			fiHI  = svs85_fion_HI(xj)
			fiHeI = svs85_fion_HeI(xj)
		else
			fh = 1.0d0; fiHI = 0.0d0; fiHeI = 0.0d0
		endif
		acc_secHI  = 0.0d0
		acc_secHeI = 0.0d0

		! Photoheating integral. Accumulate the inner bracket in the
		! original order (H, He, then metals in iphot order) and apply
		! the int_f factor once, preserving the original association. Where a
		! photoelectron energy E0 = e_v - E_th exceeds E_sec_ion, only f_heat(x)
		! of its excess is deposited as heat (fhv) and the balance drives H I /
		! He I secondary ionizations; below the threshold it thermalizes fully.
		! fhv = 1 when the coupling is off, so the heating integrand is
		! bit-identical to the legacy full-thermalization path. He I triplet
		! photoionization (threshold e_th_HeTR = 4.8 eV) deposits its
		! photoelectron energy here as well, consistently with its opacity
		! and its P_HeITR rate.
		fhv = 1.0d0
		if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_HI + E_sec_ion)
		acc_H = (1.0-e_th_HI/e_v)*fhv*s_hi*nhi(j)
		if (present(heat_chan)) acc_HI = (1.0-e_th_HI/e_v)*fhv*s_hi*nhi(j)
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_hi*nhi(j)/e_v *                       &
			     merge(fiHI *(e_v-e_th_HI)/e_th_HI , 0.0d0, e_v > e_th_HI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hi*nhi(j)/e_v *                       &
			     merge(fiHeI*(e_v-e_th_HI)/e_th_HeI, 0.0d0, e_v > e_th_HI + E_sec_ion)
		endif

		fhv = 1.0d0
		if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_HeI + E_sec_ion)
		acc_H = acc_H + (1.0-e_th_HeI/e_v)*fhv*s_hei*nheiS(j)
		if (present(heat_chan)) acc_HeI = (1.0-e_th_HeI/e_v)*fhv*s_hei*nheiS(j)
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_hei*nheiS(j)/e_v *                    &
			     merge(fiHI *(e_v-e_th_HeI)/e_th_HI , 0.0d0, e_v > e_th_HeI + E_sec_ion)
			acc_secHeI = acc_secHeI + s_hei*nheiS(j)/e_v *                    &
			     merge(fiHeI*(e_v-e_th_HeI)/e_th_HeI, 0.0d0, e_v > e_th_HeI + E_sec_ion)
		endif

		fhv = 1.0d0
		if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_HeII + E_sec_ion)
		acc_H = acc_H + (1.0-e_th_HeII/e_v)*fhv*s_heii*nheii(j)
		if (present(heat_chan)) acc_HeII = (1.0-e_th_HeII/e_v)*fhv*s_heii*nheii(j)
		if (sec_on) then
			acc_secHI  = acc_secHI  + s_heii*nheii(j)/e_v *                   &
			     merge(fiHI *(e_v-e_th_HeII)/e_th_HI , 0.0d0, e_v > e_th_HeII + E_sec_ion)
			acc_secHeI = acc_secHeI + s_heii*nheii(j)/e_v *                   &
			     merge(fiHeI*(e_v-e_th_HeII)/e_th_HeI, 0.0d0, e_v > e_th_HeII + E_sec_ion)
		endif

		! He I 2^3S (triplet): photoelectron energy hv - 4.8 eV, same
		! secondary partition as the other absorbers.
		if (thereis_HeITR) then
			fhv = 1.0d0
			if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_HeTR + E_sec_ion)
			acc_H = acc_H + (1.0-e_th_HeTR/e_v)*fhv*s_heiTR*nheiTR(j)
			if (present(heat_chan)) acc_HeTR = (1.0-e_th_HeTR/e_v)*fhv*s_heiTR*nheiTR(j)
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_heiTR*nheiTR(j)/e_v *             &
				     merge(fiHI *(e_v-e_th_HeTR)/e_th_HI , 0.0d0, e_v > e_th_HeTR + E_sec_ion)
				acc_secHeI = acc_secHeI + s_heiTR*nheiTR(j)/e_v *             &
				     merge(fiHeI*(e_v-e_th_HeTR)/e_th_HeI, 0.0d0, e_v > e_th_HeTR + E_sec_ion)
			endif
		endif

		if (present(nh2)) then
			fhv = 1.0d0
			if (sec_on) fhv = merge(fh, 1.0d0, e_v > e_th_H2 + E_sec_ion)
			acc_H = acc_H + (1.0-e_th_H2/e_v)*fhv*s_h2*nh2(j)
			if (present(heat_chan)) acc_H2 = (1.0-e_th_H2/e_v)*fhv*s_h2*nh2(j)
			if (sec_on) then
				acc_secHI  = acc_secHI  + s_h2*nh2(j)/e_v *                   &
				     merge(fiHI *(e_v-e_th_H2)/e_th_HI , 0.0d0, e_v > e_th_H2 + E_sec_ion)
				acc_secHeI = acc_secHeI + s_h2*nh2(j)/e_v *                   &
				     merge(fiHeI*(e_v-e_th_H2)/e_th_HeI, 0.0d0, e_v > e_th_H2 + E_sec_ion)
			endif
		endif

		do i = 1,n_mion
			if (.not. mion_isphot(i)) cycle
			k = mion_iphot(i)
			fhv = 1.0d0
			if (sec_on) fhv = merge(fh, 1.0d0, e_v > mion_ethr(i) + E_sec_ion)
			acc_H = acc_H + (1.0-mion_ethr(i)/e_v)*fhv*sigma_tab(:,k)*nm(j,i)
			if (present(heat_chan)) acc_mtl = acc_mtl + (1.0-mion_ethr(i)/e_v)*fhv*sigma_tab(:,k)*nm(j,i)
			if (sec_on) then
				acc_secHI  = acc_secHI  + sigma_tab(:,k)*nm(j,i)/e_v *        &
				     merge(fiHI *(e_v-mion_ethr(i))/e_th_HI , 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
				acc_secHeI = acc_secHeI + sigma_tab(:,k)*nm(j,i)/e_v *        &
				     merge(fiHeI*(e_v-mion_ethr(i))/e_th_HeI, 0.0d0, e_v > mion_ethr(i) + E_sec_ion)
			endif
		enddo
		int_H = int_f*acc_H

		! Absorbed energy integral (same ordering as int_H)
		acc_q = s_hi *nhi  (j)
		acc_q = acc_q + s_hei *nheiS(j)
		acc_q = acc_q + s_heii*nheii(j)
		if (thereis_HeITR) acc_q = acc_q + s_heiTR*nheiTR(j)
		if (present(nh2)) acc_q = acc_q + s_h2*nh2(j)
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
		PIR_H2 = 0.0
		if (present(nh2)) PIR_H2 = sum(int_h2*de_v)
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

		! Save into vectors. No synchronization needed: each thread writes
		! only its own j elements of the shared arrays (the former OMP
		! CRITICAL serialized the loop for no correctness benefit; removed
		! per the 2026-07-02 review's performance note -- values unchanged).
    	P_HI(j)    = PIR_1  *1.0e-18*erg2eV
		if (present(P_H2)) P_H2(j) = PIR_H2*1.0e-18*erg2eV
    	P_HeI(j)   = PIR_15 *1.0e-18*erg2eV
		! Add the H I / He I secondary-ionization rates from fast photoelectrons.
		if (sec_on) then
			R_secHI  = sum(int_f*acc_secHI *de_v)*1.0e-18*erg2eV
			R_secHeI = sum(int_f*acc_secHeI*de_v)*1.0e-18*erg2eV
			! x -> 1 gives f_ion -> 0, so R_sec -> 0 there; the max() is only a
			! divide-by-zero guard for an (unphysical) fully depleted cell.
			P_HI(j)  = P_HI(j)  + R_secHI /max(nhi(j)  , 1.0d-99)
			P_HeI(j) = P_HeI(j) + R_secHeI/max(nheiS(j), 1.0d-99)
		endif
		P_HeII(j)  = PIR_2  *1.0e-18*erg2eV
		P_HeITR(j) = PIR_TR *1.0e-18*erg2eV
		P_m(j,:)   = Pm_loc(:)
		heat(j)    = Hea_1*1.0e-18
		! Per-absorber photoheating split (same 1e-18 factor and int_f weight
		! as heat above). Columns 1-6 sum to heat(j) up to rounding.
		if (present(heat_chan)) then
			heat_chan(j,1) = sum(int_f*acc_HI  *de_v)*1.0e-18
			heat_chan(j,2) = sum(int_f*acc_HeI *de_v)*1.0e-18
			heat_chan(j,3) = sum(int_f*acc_HeII*de_v)*1.0e-18
			heat_chan(j,4) = sum(int_f*acc_HeTR*de_v)*1.0e-18
			heat_chan(j,5) = sum(int_f*acc_H2  *de_v)*1.0e-18
			heat_chan(j,6) = sum(int_f*acc_mtl *de_v)*1.0e-18
		endif
		! Guard against q_abs = 0 (see PH_heat_H).
		q(j)       = Hea_1/max(q_abs, 1.0d-99)

	enddo
	!$OMP END PARALLEL DO

	! End of subroutine
	end subroutine PH_heat_HHe
	
	!----------------------------------!
	
	subroutine eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,        &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,             &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,      &
				   cool, cool_chan, nheiTR, a_ion_HeITR)

	! Evaluate the cooling rate contributions to energy and
	!	rate equations. Includes H, He, and metal channels. Metals are
	!	driven from species_table metadata (canonical ion order); nm,
	!	rec_m and aion_m carry the densities and rates for each ion, so the
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
   real*8, dimension(1-Ng:N+Ng) :: tau_eff,beta_esc     ! beta escape prob.
   real*8, dimension(1-Ng:N+Ng) :: kappa_loc,dr_cm
	real*8, dimension(1-Ng:N+Ng) :: ne		  			 ! Electron number density
	real*8, dimension(1-Ng:N+Ng) :: GF_z1,GF_z2			 ! free-free Gaunt at Z_ion=1,2
	real*8 :: Cdex_OI,Cdex_CII                           ! 2-level collis. de-exc.
	! Named bridges for the (verbatim) two-level cooling branch
	real*8, dimension(1-Ng:N+Ng) :: nci,ncii,noi,noii,nmgi,nmgii
	real*8, dimension(1-Ng:N+Ng) :: c_CI,c_OII,c_MgI,c_MgII

   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB
   ! Metal recombination rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: rec_m

	! Recombination cooling coefficients
	real*8, dimension(1-Ng:N+Ng) :: coeff_rec_cool_HII,  &
									coeff_rec_cool_HeII, &
									coeff_rec_cool_HeIII

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII
   ! Metal collisional ionization rates for each ion (canonical order)
   real*8, dimension(1-Ng:N+Ng,n_mion),intent(out) :: aion_m

	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_rate_HI,    &
	 								coeff_coex_rate_HeI,   &
									coeff_coex_rate_HeII
	! He 2^3S metastable cooling coefficients (triplet-tracking callers only):
	! 10830 A collisional-excitation cooling, and the q31a/q31b conversion
	! rate coefficients reused for their thermal-energy ledger.
	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_HeI23S_10830, q31a_l, q31b_l
									 
	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  cool

	! Optional cooling breakdown in each channel (cgs erg cm^-3 s^-1, same
	! units as `cool`). Columns 1-6 = H/He recombination, collisional
	! ionization, then the collisional-excitation channel split into its
	! three absorbers -- H I (the Lyman-alpha-dominated H-line cooling),
	! He I, He II -- and finally bremsstrahlung (incl. metal-ion charges);
	! columns 6+i = metal ion i line cooling (0 for non-coolant ions).
	! The He I column also carries the He 2^3S metastable collisional cooling
	! (10830 A + singlet-conversion terms), so columns 3-5 sum exactly to
	! ne*coex. This is an exact decomposition of `cool` in the default
	! (.not.use_2lev_cool) branch; in the two-level branch the metal terms for
	! each ion are the resonance-line approximation and need not sum to cool_M.
	real*8, dimension(1-Ng:N+Ng,6+n_mion),intent(out),optional :: cool_chan

	! He 2^3S metastable density [cm^-3], present only for the triplet-tracking
	! callers. When supplied it adds the collisional-ionization cooling of the
	! 2^3S state (4.8 eV per event, ci_HeI23S) to the CI channel. a_ion_HeITR
	! returns the He(2^3S) collisional-ionization rate coefficient [cm^3 s^-1]
	! for the ionization equations, mirroring a_ion_HI/HeI/HeII.
	real*8, dimension(1-Ng:N+Ng),intent(in),optional  :: nheiTR
	real*8, dimension(1-Ng:N+Ng),intent(out),optional :: a_ion_HeITR

	! He 2^3S collisional-ionization rate coefficient (always computed; only
	! exported / applied through the optional arguments above).
	real*8, dimension(1-Ng:N+Ng) :: aion_HeITR

   ! Free electron density (incl. metal electrons under eos_metals). Molecular
   ! ions are deliberately omitted as trace electron donors (negligible in the
   ! hot, atomic gas where eval_cool operates).
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
	call ci_HeI23S(T_K,aion_HeITR)       ! He 2^3S metastable (4.8 eV threshold)
	if (present(a_ion_HeITR)) a_ion_HeITR = aion_HeITR

	! Cooling rate. The prefactors are the ionization potentials in erg
	! (13.6/24.6/54.4 eV for HI/HeI/HeII; e_th_HeTR = 4.8 eV for the 2^3S
	! metastable, e_th_HeTR/erg2eV = 7.69e-12 erg). The 2^3S term (added only
	! when nheiTR is supplied) reproduces the Black (1981) form
	! 6.41e-21 sqrt(T) exp(-55338/T) n_e n_23S once multiplied by n_e below.
	coio =  2.179e-11*a_ion_HI*nhi  	 & ! HI
		  + 3.940e-11*a_ion_HeI*nhei 	 & ! HeI
		  + 8.715e-11*a_ion_HeII*nheii   ! HeII
	if (present(nheiTR)) coio = coio                                     &
		  + (e_th_HeTR/erg2eV)*aion_HeITR*nheiTR   ! He 2^3S
	
	!-- Bremsstrahlung --!

	! Free-free scales with the ion NET charge Z_ion (not the nuclear number):
	! H II, He II and singly-ionized metals are Z_ion = 1; He III and doubly-
	! ionized metals are Z_ion = 2. The Gaunt factor is evaluated at Z_ion, so
	! only the charge-1 and charge-2 values are needed.
	call GF(T_K, 1.0d0, GF_z1)
	call GF(T_K, 2.0d0, GF_z2)

	! Cooling rate: sum n_ion * Z_ion^2 * gbar(Z_ion) over all charged ions.
	! mion_z2 = mion_stage^2 already holds the metal charge^2 (neutral -> 0,
	! an exact +0 term); the Gaunt table is selected by mion_stage.
	brem_acc = GF_z1*nhii                       ! HII   (Z_ion = 1)
	brem_acc = brem_acc + GF_z1*nheii           ! HeII  (Z_ion = 1)
	brem_acc = brem_acc + 4.0*GF_z2*nheiii      ! HeIII (Z_ion = 2)
	do i = 1,n_mion
		if (mion_stage(i) == 2) then
			brem_acc = brem_acc + mion_z2(i)*GF_z2*nm(:,i)
		else
			brem_acc = brem_acc + mion_z2(i)*GF_z1*nm(:,i)
		endif
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

	! He 2^3S metastable collisional cooling (triplet-tracking callers only).
	! Both terms scale with the EXPLICITLY computed n_23S (nheiTR); the common
	! n_e factor is applied together with the other coex terms in `cool` below.
	!  (1) 2^3S -> 2^3P collisional excitation, then 10830 A radiative decay
	!      and photon escape (Black 1981 / Allan 2024 / Falorca & Vidotto 2026
	!      Table A2). The 10830 photon leaves the gas radially (the wind is thin
	!      to it); its line optical depth is what the transit sees. Black's
	!      implicit steady-state triplet form (~T^-0.6687 n_e n_He+) is NOT used.
	!  (2) 2^3S -> 2^1S (0.80 eV) and 2^3S -> 2^1P (1.40 eV) collisional
	!      conversions each remove their threshold energy from the electron gas;
	!      the excited singlet then decays radiatively (that photon is not
	!      thermal). Rate coefficients reused from coex_HeI_23S_21S / _21P.
	! The ground -> triplet excitation (q13, 19.82 eV) is deliberately EXCLUDED:
	! that channel is already carried by the Cen-1992 He I coex term above
	! (coeff_coex_rate_HeI*nhei); adding q13 here would double count it.
	if (present(nheiTR)) then
		call coex_rate_HeI23S_10830(T_K,coeff_coex_HeI23S_10830)
		call coex_HeI_23S_21S(T_K,q31a_l)
		call coex_HeI_23S_21P(T_K,q31b_l)
		coex = coex + ( coeff_coex_HeI23S_10830                          &
		              + (0.80d0*q31a_l + 1.40d0*q31b_l)/erg2eV )*nheiTR
	endif

	!-- Metal recombination + collisional ionization rates --!
	! (rates for ionization equilibrium; not part of cool here.)
	! Filled per canonical ion via the metadata dispatchers, which call
	! the same routines for each ion as before; ions with no entry (inert top
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
	! cooling for each ion (collider-independent, so the electron-only SE solve is
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

	! Breakdown by channel for the diagnostic (Huang Fig. 10).
	! Read straight from the arrays already computed above, so the sum of
	! all channels reproduces `cool` exactly in the default branch.
	if (present(cool_chan)) then
		cool_chan(:,1) = ne*reco
		cool_chan(:,2) = ne*coio
		! Collisional excitation split by absorber. H I is the
		! Lyman-alpha-dominated H-line cooling; He II is its own term. The
		! He I column is taken as the remainder ne*coex - HI - HeII so that it
		! also absorbs the He 2^3S metastable terms folded into coex above,
		! keeping columns 3-5 an exact split of ne*coex.
		cool_chan(:,3) = ne*(coeff_coex_rate_HI*nhi)      ! coex_HI [Lya]
		cool_chan(:,5) = ne*(coeff_coex_rate_HeII*nheii)  ! coex_HeII
		cool_chan(:,4) = ne*coex - cool_chan(:,3) - cool_chan(:,5)  ! coex_HeI
		cool_chan(:,6) = ne*brem
		do i = 1,n_mion
			if (mion_iscool(i)) then
				cool_chan(:,6+i) = beta_esc*ne*nm(:,i)*c_metal(:,i)
			else
				cool_chan(:,6+i) = 0.0d0
			endif
		enddo
	endif

	! End of subroutine
	end subroutine eval_cool

	! ------------------------------------------------------------- !

	subroutine write_cool_breakdown_eq(T_in,n_in,f_sp_in)
	! Diagnostic. Dump the radiative cooling rate in each channel vs
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
	real*8, dimension(1-Ng:N+Ng) :: nhi,nhii,nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	real*8, dimension(1-Ng:N+Ng) :: cool,csum,rel
	real*8, dimension(1-Ng:N+Ng,6+n_mion) :: chan
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
	! He 2^3S metastable, so the collisional-ionization channel matches the
	! main cool column (0 when the triplet is not tracked).
	if (thereis_HeITR) then
		nheiTR = f_sp_in(:,6)*n_dim
	else
		nheiTR = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	! Molecular ions are deliberately omitted as trace electron donors
	! (negligible in this diagnostic's hot, atomic gas).
	call calc_ne(nhii,nheii,nheiii,ne,nm)

	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,            &
				   rchiiB,rcheiiB,rcheiiiB, rec_m,                &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
				   cool, cool_chan = chan, nheiTR = nheiTR)

	! Internal consistency: channel sum vs total cool (default branch -> ~eps)
	csum = 0.0d0
	do i = 1,6+n_mion
		csum = csum + chan(:,i)
	enddo
	rel    = abs(csum - cool)/max(abs(cool),1.0d-99)
	maxrel = maxval(rel(1:N))
	write(*,'(a,es9.2)')                                                  &
		' (write_cool_breakdown_eq) max |sum(channels)/cool - 1| = ', maxrel

	open(unit = 71, file = './output/Cooling_breakdown.txt')
	write(71,'(a)') '# Radiative cooling rate in each channel [cgs erg cm^-3 s^-1] vs radius.'
	write(71,'(a)') '# Channel sum reproduces the Hydro_ioniz.txt cool column.'
	write(71,'(a)') '# col1 r/Rp  col2 T[K]  col3 ne  col4 cool_total  col5 reco'  &
	             // '  col6 coio  col7 coex_HI[Lya]  col8 coex_HeI  col9 coex_HeII'  &
	             // '  col10 brem  then one col per metal ion:'
	write(71,'(a)',advance='no') '#   metal-ion columns (canonical order):'
	do i = 1,n_mion
		write(71,'(1x,a)',advance='no') trim(mion_name(i))
	enddo
	write(71,*)
	do j = 1-Ng,N+Ng
		write(71,*) r(j), T_K(j), ne(j), cool(j),                        &
		            chan(j,1), chan(j,2), chan(j,3), chan(j,4),           &
		            chan(j,5), chan(j,6),                                 &
		            (chan(j,6+i), i = 1,n_mion)
	enddo
	close(71)

	end subroutine write_cool_breakdown_eq

	! ------------------------------------------------------------- !

	subroutine write_heat_breakdown_eq(T_in,n_in,f_sp_in)
	! Diagnostic. Dump the volumetric heating rate in each channel vs
	! radius for the converged equilibrium state, recomputing the same
	! photoheating (PH_heat_HHe) the solver uses plus the excited-H Balmer,
	! He-recombination, and He(2^3S) Penning heating terms added in
	! ionization_equilibrium. All in cgs erg cm^-3 s^-1; the channel sum
	! reproduces the total heating (heat_total column) and, up to convergence,
	! the Hydro_ioniz.txt heat column. Photoheating columns: H I, He I,
	! He II, He 2^3S, H2, metals (sum over photo-ionizable metal ions). Then
	! the excited-H photoelectric (Hpe) and Lyman-alpha de-excitation (Hdx)
	! heating, He-recombination-driven H heating, and He(2^3S)+H Penning
	! heating. The printed max relative residual is the internal consistency
	! check on the photoheating split.

	integer :: j,i,im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in

	real*8, dimension(1-Ng:N+Ng) :: n_dim,T_K,ne
	real*8, dimension(1-Ng:N+Ng) :: nhi,nhii,nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng) :: nh,nhe,xion
	real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
	real*8, dimension(1-Ng:N+Ng,6) :: hchan
	real*8, dimension(1-Ng:N+Ng) :: heat_ph,heat_tot,csum,rel
	real*8, dimension(1-Ng:N+Ng) :: h_hrc,h_penning
	! Throwaway PH_heat_HHe rate outputs (not needed for the dump)
	real*8, dimension(1-Ng:N+Ng) :: P_HI,P_HeI,P_HeII,P_HeITR,q
	real*8, dimension(1-Ng:N+Ng,n_mion) :: P_m
	! He-recombination coupling / triplet scratch
	real*8, dimension(1-Ng:N+Ng) :: rcheiTR,rcheii,q13,q31a,q31b,Q31
	real*8, dimension(1-Ng:N+Ng) :: rcheiiB_hrc,dP_HI_hrc
	real*8 :: A31,maxrel

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
	if (thereis_HeITR) then
		nheiTR = f_sp_in(:,6)*n_dim
	else
		nheiTR = 0.0d0
	endif
	do im = 1,n_mion
		nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
	enddo
	call calc_ne(nhii,nheii,nheiii,ne,nm)

	! Ionized fraction for the SvS85 secondary-ionization partition (atomic
	! form; molecular donors are omitted -- this diagnostic targets atomic
	! runs). See ionization_equilibrium for the exact expression.
	if (thereis_mol) write(*,'(a)') ' (write_heat_breakdown_eq) NOTE: molecular '  &
		// 'run -- H2 and molecular Penning heating channels are omitted.'
	nh   = nhi + nhii
	nhe  = nhei + nheii + nheiii
	xion = min(max((nhii + nheii + nheiii)/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Photoheating split (same call the solver makes, atomic path).
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,             &
		         P_HI,P_HeI,P_HeII,P_HeITR, P_m,                      &
		         heat_ph,q, heat_chan = hchan)
	else
		call PH_heat_H(nhi, xion, P_HI, heat_ph, q)
		hchan = 0.0d0
		hchan(:,1) = heat_ph
	endif

	! He-recombination-driven H heating and He(2^3S)+H Penning heating,
	! reconstructed exactly as ionization_equilibrium adds them to `heat`.
	A31 = 0.0d0; q31a = 0.0d0; q31b = 0.0d0; Q31 = 0.0d0
	if (thereis_HeITR) &
		call HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31a,q31b,Q31)
	h_hrc = 0.0d0
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,       &
		                     A31, q31a, q31b,                         &
		                     rcheiiB_hrc, dP_HI_hrc, h_hrc)
	endif
	h_penning = 0.0d0
	if (thereis_HeITR) h_penning =                                    &
		nheiTR*nhi*Q31*(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

	! Total heating (independent of the per-channel columns; the residual
	! below checks the photoheating decomposition against heat_ph).
	heat_tot = heat_ph + Hpe_arr + Hdx_arr + h_hrc + h_penning

	! Internal consistency of the photoheating split.
	csum   = hchan(:,1) + hchan(:,2) + hchan(:,3) + hchan(:,4)        &
	       + hchan(:,5) + hchan(:,6)
	rel    = abs(csum - heat_ph)/max(abs(heat_ph),1.0d-99)
	maxrel = maxval(rel(1:N))
	write(*,'(a,es9.2)')                                                  &
		' (write_heat_breakdown_eq) max |sum(photo channels)/heat_photo - 1| = ', maxrel

	open(unit = 72, file = './output/Heating_breakdown.txt')
	write(72,'(a)') '# Volumetric heating rate in each channel [cgs erg cm^-3 s^-1] vs radius.'
	write(72,'(a)') '# Channel sum reproduces the heat_total column (and the Hydro_ioniz.txt'  &
	             // ' heat column up to convergence).'
	write(72,'(a)') '# col1 r/Rp  col2 T[K]  col3 ne  col4 heat_total  col5 heat_HI'  &
	             // '  col6 heat_HeI  col7 heat_HeII  col8 heat_He23S  col9 heat_H2'  &
	             // '  col10 heat_metals  col11 heat_Hpe[excitedH]'                    &
	             // '  col12 heat_Hdx[Lya-deexc]  col13 heat_He_recomb  col14 heat_He23S_Penning'
	do j = 1-Ng,N+Ng
		write(72,*) r(j), T_K(j), ne(j), heat_tot(j),                    &
		            hchan(j,1), hchan(j,2), hchan(j,3), hchan(j,4),       &
		            hchan(j,5), hchan(j,6),                               &
		            Hpe_arr(j), Hdx_arr(j), h_hrc(j), h_penning(j)
	enddo
	close(72)

	end subroutine write_heat_breakdown_eq

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

	! ------------------------------------------------------------- !

	! He recombination radiation ionizing H I (Draine 2011 on-the-spot y/z;
	! docs/QUESTIONS_2026-07-17.md). Given the pre-solve (lagged) densities and
	! the rate coefficients, returns the He II recombination coefficient the
	! ionization balance should use (rcheiiB_new), the extra H I photoionization
	! rate [s^-1] (dP_HI), and the extra photoelectron heating [erg cm^-3 s^-1]
	! (dheat). All three are zero when the flag is off; the caller applies them.
	!
	! The corrections depend only on T and the densities: alpha_1 and alpha_B are
	! recomputed here from T (matching the active case-B function eval_cool uses),
	! so the passed rcheiiB is not read -- the caller may call this both before
	! (rate/rate-limited) and after (heating) it overwrites rcheiiB.
	!
	! atomic mode (thereis_HeITR = .false.):
	!   alpha_1 = Mao & Kaastra ground (1s^2) capture; alpha_B = active He II
	!   case B (rec_HeII_B). rcheiiB_new = alpha_B + y alpha_1 (Draine Eq. 14.17);
	!   dP_HI = n_HeII n_e [z alpha_B + y alpha_1]/n_HI; the heating uses the
	!   ground photoelectron energy E_gnd = 24.6-13.6 = 11.0 eV and a
	!   cascade-averaged E_casc ~ 6.3 eV.
	! TR mode (thereis_HeITR = .true.): channels are explicit -- alpha_1 is the
	!   1^1S channel (rec_HeII_11S, which HeITR_coeffs writes into rcheiiB) and
	!   the 2^3S / singlet cascade rates (A31, q31a, q31b, n_2^3S = nheiTR) drive
	!   the H-ionizing photon production directly, so z is not used (a 2^3S
	!   destroyed by photoionization/Penning emits no 19.8 eV photon).
	subroutine he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,        &
	                           A31, q31a, q31b,                          &
	                           rcheiiB_new, dP_HI, dheat)

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi, nhei, nheii, &
	                                              nheiTR, ne, q31a, q31b
	real*8,                       intent(in)  :: A31
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiiB_new, dP_HI, dheat

	real*8, dimension(1-Ng:N+Ng) :: alpha1, alphaB
	real*8 :: Rratio, y, z, T4, ncrit, P_add, H_add
	integer :: j

	! Photoionization cross-section ratio sigma_He/sigma_H at the 24.6 eV He I
	! ground edge, used to split the >= 24.6 eV ground-capture continuum between
	! H and He (Draine Eq. 14.16). Its (small) kT dependence is neglected; R ~ 6.
	Rratio = sigma_HeI(24.6d0)/max(sigma(24.6d0,1.0d0), 1.0d-99)

	! alpha_B is the active He II case B (whichever fit eval_cool uses).
	call rec_HeII_B(T_K, alphaB)

	if (.not. thereis_HeITR) then
		! atomic (case-B) mode: alpha_1 = Mao & Kaastra ground capture.
		alpha1 = alpha1_HeII_mao(T_K)
		do j = 1-Ng,N+Ng
			! y: fraction of >= 24.6 eV ground-capture photons ionizing H.
			y = 1.0d0/(1.0d0 + Rratio*nhei(j)/max(nhi(j),1.0d-99))
			! z: density-dependent fraction of case-B cascade photons ionizing H;
			! 0.96 (low density, 19.8 eV line ionizes H) -> 0.67 (2^3S
			! collisionally converted to singlets) via the 2^3S critical density
			! (documented interpolation between Draine's two limits).
			T4    = T_K(j)/1.0d4
			ncrit = 1100.0d0*exp(1.2d0/T4)*sqrt(T4)
			z     = 0.67d0 + 0.29d0/(1.0d0 + ne(j)/ncrit)
			! He II recombination: alpha_eff = alpha_B + y alpha_1.
			rcheiiB_new(j) = alphaB(j) + y*alpha1(j)
			! Extra H I photoionization rate [s^-1].
			dP_HI(j) = nheii(j)*ne(j)*(z*alphaB(j) + y*alpha1(j))          &
			           /max(nhi(j),1.0d-99)
			! Photoelectron heating [erg cm^-3 s^-1]. E_casc = 6.3 eV is a
			! low-density channel-weighted average (~0.75*6.2 + 0.17*7.6 +
			! 0.08*3.0), approximate; E_gnd = 11.0 eV (24.6-13.6).
			dheat(j) = nheii(j)*ne(j)*(z*alphaB(j)*6.3d0                   &
			           + y*alpha1(j)*11.0d0)/erg2eV
		enddo
	else
		! TR mode: alpha_1 = 1^1S channel (rec_HeII_11S).
		call rec_HeII_11S(T_K, alpha1)
		do j = 1-Ng,N+Ng
			y = 1.0d0/(1.0d0 + Rratio*nhei(j)/max(nhi(j),1.0d-99))
			! (1) 1^1S channel coefficient: net ground capture (y alpha_1) plus
			! the singlet-excited capture channel (0.25 alpha_B) missing from the
			! current network.
			rcheiiB_new(j) = y*alpha1(j) + 0.25d0*alphaB(j)
			! (2) H-ionizing photon production [cm^-3 s^-1]:
			P_add = y*alpha1(j)*nheii(j)*ne(j)                    ! ground (>=24.6)
			! singlet-excited captures: 0.85 = 2/3*1.0 (584 A resonance) +
			! 1/3*0.56 (2^1S two-photon fraction above 13.6 eV).
			P_add = P_add + 0.85d0*0.25d0*alphaB(j)*nheii(j)*ne(j)
			! 2^3S radiative decay (19.8 eV line, always ionizes H).
			P_add = P_add + A31*nheiTR(j)
			! 2^3S collisionally converted to singlets, then decaying: 2^1S
			! two-photon (0.56 ionizing) + 2^1P -> 584 A (1.0 ionizing).
			P_add = P_add + ne(j)*nheiTR(j)*(q31a(j)*0.56d0 + q31b(j)*1.0d0)
			dP_HI(j) = P_add/max(nhi(j),1.0d-99)
			! (3) photoelectron heating [erg cm^-3 s^-1], channel E_dep [eV]:
			! ground 11.0; singlet-excited 5.6 (= 2/3*7.6 + 1/3*0.56*3.0, where
			! 584 A -> 7.6, two-photon -> 3.0); 19.8 eV line -> 6.2; 2^3S coll.
			! -> singlet: two-photon 3.0 + 584 A 7.6.
			H_add = y*alpha1(j)*nheii(j)*ne(j)*11.0d0
			H_add = H_add + 0.25d0*alphaB(j)*nheii(j)*ne(j)*5.6d0
			H_add = H_add + A31*nheiTR(j)*6.2d0
			H_add = H_add + ne(j)*nheiTR(j)                                &
			        *(q31a(j)*0.56d0*3.0d0 + q31b(j)*1.0d0*7.6d0)
			dheat(j) = H_add/erg2eV
		enddo
	endif

	! End of subroutine
	end subroutine he_rec_coupling

	! ------------------------------------------------------------- !

	! Shull & van Steenberg (1985, ApJ 298, 268) high-energy asymptotic
	! partition of a fast photoelectron's excess energy. x = ionized fraction.
	! A possible future refinement is to couple the SvS85 Ly-alpha excitation
	! channel f_exc,Lya = 0.4766*(1-x^0.2735)^1.5221 to the Ly-alpha field of
	! the excited-H model; here that energy is assumed to escape as line
	! radiation and is not put into any rate.
	elemental function svs85_fheat(x) result(f)
		real*8, intent(in) :: x
		real*8 :: f
		f = 0.9971d0*(1.0d0 - (1.0d0 - x**0.2663d0)**1.3163d0)
	end function svs85_fheat

	elemental function svs85_fion_HI(x) result(f)
		real*8, intent(in) :: x
		real*8 :: f
		f = 0.3908d0*(1.0d0 - x**0.4092d0)**1.7592d0
	end function svs85_fion_HI

	elemental function svs85_fion_HeI(x) result(f)
		real*8, intent(in) :: x
		real*8 :: f
		f = 0.0554d0*(1.0d0 - x**0.4614d0)**1.6660d0
	end function svs85_fion_HeI

	! End of module
	end module utils_ion_eq

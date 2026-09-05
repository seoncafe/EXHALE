   module energy_vectors_construct
   ! Construct energy grid for the computation of radiative
	!      contributions and radiative equilibrium
      
   use global_parameters
   use species_table, only: n_mphot, n_melem, n_mion, mion_ethr, melem_i0
   use electron_energy_degradation, only: photoelectron_energy_grid,      &
                                          n_abs_fixed
   use sed_reader
   use J_incident
   use Cross_sections
   use opacity_input            ! opacity.inp reader
   use opacity_models           ! A/C/P/T dispatcher
   ! The four mutually exclusive final-state channels of the H2 absorption
   ! cross section, and the switch that selects the double-ionization model.
   use h2_photo_channels, only: h2_channel_cross_sections, n_h2_channels,  &
                                h2_double_ionization_model,                &
                                ICH_M, ICH_S, ICH_D, ICH_N

   implicit none
      
   contains
      
   subroutine set_energy_vectors
	! Subroutine to construct the energy grid

   integer :: i,j,k
   real*8  :: xi_day
	! Ionization threshold of every absorber that produces a photoelectron,
	! in the order electron_energy_degradation indexes them.
	real*8, dimension(n_abs_fixed + n_mion) :: eth_abs
	integer, parameter :: num_HI = 50
	integer, parameter :: num_HeI = 50
	integer, parameter :: num_HeII = 50
	integer, parameter :: num_X = 50
	
	real*8 :: e_min, e_max
	real*8, dimension(:), allocatable :: dum_E,dum_F,dum_dE
	real*8 :: temp1,temp2
	real*8 :: e_sub_low    ! lower edge of the below-13.6 eV sub-grid
	! The four H2 channel cross sections at one energy, and their ionizing
	! part (channels M + S + D).
	real*8 :: sig4(n_h2_channels), sig_ion

	! Set the energy band of the spectrum
	if (is_PL_sed) then

		! Points below 13.6 eV. Two sources need them, and they can be
		! active together (the merged triplet+metals systems):
		!  - He triplet (HeI 2^3S), threshold 4.8 eV;
		!  - low-IP metals (e.g. Mg I, threshold 7.646 eV; K I, 4.341 eV).
		! The grid floor is the LOWEST threshold over whichever are active,
		! as in the loaded-SED path (sed_read): stopping at the triplet edge
		! would leave K I, whose threshold is below it, without ionizing
		! photons and hold it spuriously neutral.
		NlTR      = 0
		e_sub_low = e_th_HI
		if (thereis_HeITR) then
			NlTR      = 20
			e_sub_low = e_th_HeTR
		endif
		if (thereis_lowIP_metal) then
			NlTR      = 20
			do j = 1,n_melem
				if (melem_ab(j) .gt. 0.0d0 .and.                  &
				    mion_ethr(melem_i0(j)) .lt. e_th_HI)          &
					e_sub_low = min(e_sub_low, mion_ethr(melem_i0(j)))
			enddo
		endif
		Nl = Nl_fix + NlTR

		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(F_XUV(Nl))

		! --- Construct energy grid --- !

		! Points between e_sub_low and 13.6 eV (HeI triplet or low-IP metal)
		if (NlTR .gt. 0) then
			do j = 1,NlTR
				e_v(j) = e_sub_low*	&
					  (e_th_HI/e_sub_low)**((j-1.0)/(NlTR*1.0))
		  	enddo
		endif

		! Pure-HI grid [13.6 eV, 24.6 eV]
		e_min = e_th_HI
		e_max = e_th_HeI
		do j = 1,num_HI
			e_v(NlTR + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HI*1.0)) 
		enddo 
		 
		! HI+HeI grid [24.6 eV, 54.4 eV]
		e_min = e_th_HeI
		e_max = e_th_HeII
		do j = 1,num_HeI
			e_v(NlTR + num_HI + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HeI*1.0)) 
		enddo 
		 
		! HI+HeI+HeII grid	[54.4 eV,124 eV] 
		e_min = e_th_HeII
		e_max = e_mid
		do j = 1,num_HeII
			e_v(NlTR + num_HI + num_HeI + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HeII*1.0)) 
		enddo 
		 
		! X-ray grid [124 eV, e_XUV_top]	 
		e_min = e_mid
		e_max = e_top
		do j = 1,num_X
			e_v(NlTR + num_HI + num_HeI + num_HeII + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_X*1.0))
		enddo 

		! Construct bin width
		de_v(1)      = 0.5*(e_v(2) - e_v(1))
		de_v(2:Nl-1) = 0.5*(e_v(3:Nl) - e_v(1:Nl-2))
		de_v(Nl)     = 0.5*(e_v(Nl) - e_v(Nl-1))
	
	else if (is_monochr) then	! If monochromatic radiation
		
		! Set only one wavelength
		Nl = 1
		
		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(F_XUV(Nl))

		! Define the only energy value according to the input
		e_v(1)  = e_low
		de_v(1) = 1.0
		
		! Define the total flux at this energy
		F_XUV(1) = 10.0**LEUV/(4.0*pi*a_orb**2.0)
	
	else if (do_read_sed) then    ! If read sed
	 
		! Read the SED file
	 	call read_sed
 	
		! Flip energy vectors
		allocate(dum_E(Nl),dum_F(Nl),dum_dE(Nl))
		
		dum_E  = e_v
		dum_F  = F_XUV
		dum_dE = de_v
		
		do i = 1,Nl
			e_v(i)   = dum_E(Nl-i+1)
			F_XUV(i) = dum_F(Nl-i+1)
			de_v(i)  = dum_dE(Nl-i+1)
		enddo
	
	 	deallocate(dum_E,dum_F,dum_dE)
	endif
	

	! Calculate the integrated value of the flux
	J_XUV  = (10.0**LX + 10.0**LEUV)/(4.0*pi*a_orb**2.0)

	! Read optional opacity.inp and load any .opa tables, then
	! fill the H/He cross sections through the model dispatcher
	! (model 'A' reproduces the analytic curves exactly).
	call read_opacity_input
	call load_opacity_tables

	! HI photoioiniz. cross section
	s_hi = (/ (photoion_sigma('HI', e_v(i)), i = 1,Nl) /)

	! HeI photoioiniz. cross section
	s_hei = (/ (photoion_sigma('HeI', e_v(i)), i = 1,Nl) /)

	! HeII photoioiniz. cross section
	s_heii = (/ (photoion_sigma('HeII', e_v(i)), i = 1,Nl) /)

	! HeI triplet photoioiniz. cross section
	s_heiTR = (/ (photoion_sigma('HeITR', e_v(i)), i = 1,Nl) /)

	! H2 TOTAL photoabsorption cross section (molecular extension; the
	! Backx et al. 1976 and Samson & Haddad 1994 tables below 300 eV and
	! the Yan et al. 1998 Eq. 19 sum-rule tail above it, zero below
	! 15.4 eV -- see sigma_H2 in cross_sec.f90)
	s_h2 = (/ (sigma_H2(e_v(i)), i = 1,Nl) /)
	! The cross-section vectors are allocated HERE, once, for every branch
	! above (power law, monochromatic, loaded SED): Nl is known only now.
	! Until 2026-09-05 the two spectral branches allocated them and the SED
	! branch did not; s_hi..s_h2 survived by assignment reallocation, the two
	! channel vectors by a guard here, and s_h2_di by nothing -- every
	! molecular run with a loaded spectrum stopped at its first element
	! assignment below.
	if (.not. allocated(s_h2_di)) allocate(s_h2_di(Nl))
	if (.not. allocated(s_h2_dd)) allocate(s_h2_dd(Nl))
	if (.not. allocated(s_h2_nd)) allocate(s_h2_nd(Nl))

	! The final-state channels of that same absorption. Every one of them
	! destroys one H2, and each is a SHARE of s_h2 rather than an addition
	! to it, so s_h2 -- the opacity and the total H2 destruction rate -- is
	! the same in every branch below. What differs is which fragments the
	! source terms of the molecular system are told to make:
	!   s_h2_di  H2 + hv -> H  + H+ + e-    (threshold 18.08 eV)
	!   s_h2_dd  H2 + hv -> H+ + H+ + 2e-   (threshold 51.4  eV)
	!   s_h2_nd  H2 + hv -> H  + H          (33-41 eV window, no ion)
	! and s_h2 minus the three of them drives the H2+ row.
	h2_double_ionization_model = h2_double_ionization
	if (.not. h2_neutral_dissociation .and.                              &
	    trim(h2_double_ionization) .eq. 'off') then
		! Neither of the two channels beyond the single dissociative one is
		! resolved. The split is then the single branching
		! frac_H2_dissociative_ionization, evaluated exactly as it was
		! before the channels existed, so this default reproduces the
		! previous arithmetic bit for bit.
		s_h2_di = (/ (sigma_H2(e_v(i))                                   &
		              *frac_H2_dissociative_ionization(e_v(i)), i = 1,Nl) /)
		s_h2_dd = 0.0d0
		s_h2_nd = 0.0d0
	else
		do i = 1,Nl
			call h2_channel_cross_sections(e_v(i), s_h2(i), sig4)
			if (.not. h2_neutral_dissociation) then
				! Fold the neutral share back into the three ionizing
				! channels in their own proportion, i.e. restore the unit
				! photoionization yield the pre-E1 code assumed. Each of
				! the three is linear in the ionizing part of the cross
				! section, so rescaling them is the same as evaluating
				! them with a neutral fraction of zero.
				sig_ion = sig4(ICH_M) + sig4(ICH_S) + sig4(ICH_D)
				if (sig_ion .gt. 0.0d0) then
					sig4(ICH_M) = sig4(ICH_M)*s_h2(i)/sig_ion
					sig4(ICH_S) = sig4(ICH_S)*s_h2(i)/sig_ion
					sig4(ICH_D) = sig4(ICH_D)*s_h2(i)/sig_ion
				endif
				sig4(ICH_N) = 0.0d0
			endif
			s_h2_di(i) = sig4(ICH_S)
			s_h2_dd(i) = sig4(ICH_D)
			s_h2_nd(i) = sig4(ICH_N)
		enddo
	endif

	! Metal photoionization cross sections (Verner+1996), stored in the
	! photo cross-section table in species_table iphot order.
	! The column order is metal_photoion_sigma's own (Cross_sections), which
	! is also what he_rec_coupling evaluates off this grid, so the ordering
	! of the photo cross-section table is written down in exactly one place.
	allocate(sigma_tab(Nl, n_mphot))
	do k = 1,n_mphot
		sigma_tab(:,k) = (/ (metal_photoion_sigma(k, e_v(i)), i = 1,Nl) /)
	enddo

	! --------------------------------------------------

	! Incident flux.  The dayside dilution is global_parameters'
	! dayside_dilution(), the single definition every stellar beam uses
	! (its header, and Update_EXHALE section 150).
	xi_day = dayside_dilution()
	if (is_PL_sed) then 

		F_XUV  = (/ (xi_day*J_inc(e_v(i)), i = 1,Nl) /)
		
	else 
	
		F_XUV  = xi_day*F_XUV
	
	endif      
	
	! --------------------------------------------------

	! Where each absorber's photoelectron energy E0 = hv - E_th falls in the
	! Dalgarno energy grid of the secondary-ionization partition. Both e_v
	! and the thresholds are fixed for the run, so the bracket index and the
	! interpolation weight are built once here instead of once per cell and
	! per absorber.
	eth_abs(1) = e_th_HI
	eth_abs(2) = e_th_HeI
	eth_abs(3) = e_th_HeII
	eth_abs(4) = e_th_HeTR
	eth_abs(5) = e_th_H2
	eth_abs(6) = e_th_H2_di
	do i = 1,n_mion
		eth_abs(n_abs_fixed+i) = mion_ethr(i)
	enddo
	call photoelectron_energy_grid(e_v, eth_abs)

	! End of subroutine 
    end subroutine set_energy_vectors
      
    ! End of module
    end module energy_vectors_construct


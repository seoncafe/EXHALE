   module energy_vectors_construct
   ! Construct energy grid for the computation of radiative
	!      contributions and radiative equilibrium
      
   use global_parameters
   use species_table, only: n_mphot, n_melem, mion_ethr, melem_i0
   use sed_reader
   use J_incident
   use Cross_sections
   use opacity_input            ! opacity.inp reader
   use opacity_models           ! A/C/P/T dispatcher

   implicit none
      
   contains
      
   subroutine set_energy_vectors
	! Subroutine to construct the energy grid

   integer :: i,j
	integer, parameter :: num_HI = 50
	integer, parameter :: num_HeI = 50
	integer, parameter :: num_HeII = 50
	integer, parameter :: num_X = 50
	
	real*8 :: e_min, e_max
	real*8, dimension(:), allocatable :: dum_E,dum_F,dum_dE
	real*8 :: temp1,temp2
	real*8 :: e_sub_low    ! lower edge of the below-13.6 eV sub-grid

	! Set the energy band of the spectrum
	if (is_PL_sed) then

		! Points below 13.6 eV. Two mutually exclusive sources need them:
		!  - He triplet (HeI 2^3S), threshold 4.8 eV; or
		!  - low-IP metals (e.g. Mg I, threshold 7.646 eV).
		! HeITR + metals is unsupported, so at most one applies.
		NlTR      = 0
		e_sub_low = e_th_HI
		if (thereis_HeITR) then
			NlTR      = 20
			e_sub_low = e_th_HeTR
		else if (thereis_lowIP_metal) then
			NlTR      = 20
			! Grid floor = lowest neutral-metal ionization threshold
			! below the HI edge, over active elements. Reduces to
			! e_th_MgI when Mg is the only sub-13.6 eV species present.
			do j = 1,n_melem
				if (melem_ab(j) .gt. 0.0d0 .and.                  &
				    mion_ethr(melem_i0(j)) .lt. e_th_HI)          &
					e_sub_low = min(e_sub_low, mion_ethr(melem_i0(j)))
			enddo
		endif
		Nl = Nl_fix + NlTR

		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(s_hi(Nl), s_hei(Nl),s_heii(Nl), s_h2(Nl))
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
		allocate(s_hi(Nl), s_hei(Nl),s_heii(Nl), s_h2(Nl))
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

	! H2 photoionization cross section (Tier-2 molecular extension;
	! Yan, Sadeghpour & Dalgarno 1998 fit, zero below 15.4 eV)
	s_h2 = (/ (sigma_H2(e_v(i)), i = 1,Nl) /)

	! Metal photoionization cross sections (Verner+1996), stored in the
	! photo cross-section table in species_table iphot order.
	allocate(sigma_tab(Nl, n_mphot))
	sigma_tab(:,1)  = (/ (sigma_CI  (e_v(i)), i = 1,Nl) /) ! CI
	sigma_tab(:,2)  = (/ (sigma_CII (e_v(i)), i = 1,Nl) /) ! CII
	sigma_tab(:,3)  = (/ (sigma_OI  (e_v(i)), i = 1,Nl) /) ! OI
	sigma_tab(:,4)  = (/ (sigma_OII (e_v(i)), i = 1,Nl) /) ! OII
	sigma_tab(:,5)  = (/ (sigma_NI  (e_v(i)), i = 1,Nl) /) ! NI
	sigma_tab(:,6)  = (/ (sigma_NII (e_v(i)), i = 1,Nl) /) ! NII
	sigma_tab(:,7)  = (/ (sigma_MgI (e_v(i)), i = 1,Nl) /) ! MgI
	sigma_tab(:,8)  = (/ (sigma_MgII(e_v(i)), i = 1,Nl) /) ! MgII
	sigma_tab(:,9)  = (/ (sigma_SiI (e_v(i)), i = 1,Nl) /) ! SiI
	sigma_tab(:,10) = (/ (sigma_SiII(e_v(i)), i = 1,Nl) /) ! SiII
	sigma_tab(:,11) = (/ (sigma_CaI (e_v(i)), i = 1,Nl) /) ! CaI
	sigma_tab(:,12) = (/ (sigma_CaII(e_v(i)), i = 1,Nl) /) ! CaII
	sigma_tab(:,13) = (/ (sigma_NaI (e_v(i)), i = 1,Nl) /) ! NaI
	sigma_tab(:,14) = (/ (sigma_KI  (e_v(i)), i = 1,Nl) /) ! KI
	sigma_tab(:,15) = (/ (sigma_SI  (e_v(i)), i = 1,Nl) /) ! SI
	sigma_tab(:,16) = (/ (sigma_FeI (e_v(i)), i = 1,Nl) /) ! FeI
	sigma_tab(:,17) = (/ (sigma_FeII(e_v(i)), i = 1,Nl) /) ! FeII

	! --------------------------------------------------

	! Incident flux    
	if (is_PL_sed) then 

		if (appx_mth.eq.'Rate/2 + Mdot/2') then

			F_XUV  = (/ (0.5e0*J_inc(e_v(i)), i = 1,Nl) /)

		elseif (appx_mth.eq.'Rate/4 + Mdot') then

			F_XUV  = (/ (0.25e0*J_inc(e_v(i)), i = 1,Nl) /)

		else

			F_XUV  = (/ (J_inc(e_v(i)), i = 1,Nl) /)

		endif      
		
	else 
	
		if (appx_mth.eq.'Rate/2 + Mdot/2') then

			F_XUV  = 0.5e0*F_XUV

		elseif (appx_mth.eq.'Rate/4 + Mdot') then

			F_XUV  = 0.25e0*F_XUV

		endif 
	
	endif      
	
	! End of subroutine 
    end subroutine set_energy_vectors
      
    ! End of module
    end module energy_vectors_construct


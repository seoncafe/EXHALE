   module utils
	! Collection of auxiliary subroutines

   use global_parameters
   use species_table, only: n_mion, n_mphot, mion_isphot, mion_iphot,  &
                            mion_stage, mion_elem, melem_A

   implicit none

	contains

	! ------------------------------------------------------!

	subroutine calc_ne(nhii,nheii,nheiii,ne,nm,nmol)
	! Calculate the free electron density.
	! The optional nm (per-ion metal densities, same units as nhii) adds
	! the metal electrons when the eos_metals policy is on; omitting it
	! (or eos_metals 0) reproduces the legacy H/He-only electron count.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nhii
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nheii,nheiii
	! Optional molecular ions: cols 1 H2 (neutral), 2 H2+, 3 H3+,
	! 4 HeH+ -- each molecular ion carries one electron.
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: ne
	
	if (thereis_He) then
		ne = nhii + nheii + 2.0*nheiii
	else
		ne = nhii
	endif	

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			if (mion_stage(im) .gt. 0)                                  &
				ne = ne + dble(mion_stage(im))*nm(:,im)
		enddo
	endif
	
	if (present(nmol)) ne = ne + nmol(:,2) + nmol(:,3) + nmol(:,4)

	end subroutine calc_ne
	
	! ------------------------------------------------------!

	subroutine calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot,nm,nmol)
	! Calculate the total atomic number density.
	! The optional nm adds the metal nuclei (all stages) when the
	! eos_metals policy is on.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_tot
	
	if (thereis_He) then
		n_tot = nhi + nhii + nhei + nheii + nheiii
		if (thereis_HeITR) n_tot = n_tot + nheiTR
	else
		n_tot = nhi + nhii 
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_tot = n_tot + nm(:,im)
		enddo
	endif

	! Each molecule is ONE gas particle (pressure/EOS particle count).
	if (present(nmol)) n_tot = n_tot + nmol(:,1) + nmol(:,2)              &
	                                 + nmol(:,3) + nmol(:,4)

	end subroutine calc_ntot

	! ------------------------------------------------------!
	
	subroutine calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_out,nm,nmol)
	! Calculate the total mass density (adimensional).
	! The optional nm adds the metal mass (melem_A per nucleus, all
	! stages) when the eos_metals policy is on.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_out
	
	if (thereis_He) then
		n_out = nhi + nhii + 4.0*(nhei + nheii + nheiii)
		if (thereis_HeITR) n_out = n_out + 4.0*nheiTR
	else
		n_out = nhi + nhii 
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_out = n_out + melem_A(mion_elem(im))*nm(:,im)
		enddo
	endif

	! Molecular mass: H2/H2+ = 2 m_H, H3+ = 3 m_H, HeH+ = 5 m_H (the He
	! nucleus in HeH+ is NOT in the nhei..nheiii free-He arrays).
	if (present(nmol)) n_out = n_out + 2.0d0*(nmol(:,1) + nmol(:,2))      &
	                         + 3.0d0*nmol(:,3) + 5.0d0*nmol(:,4)

	end subroutine calc_rho

	! ------------------------------------------------------!

	subroutine calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)
	! Calculates the column densities for given ionization profiles
	! 	by method of rectangles
    
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiTR
   real*8 :: dr 
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N1
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N15,N2,NTR
	
	! Initialize outputs
	N1  = 0.0
	N15 = 0.0
	N2  = 0.0
	NTR = 0.0
	
	! Outer point (opa_pf weights the opacity for the 'P' model; =1 otherwise)
	N1(N+Ng)  = dr_j(N+Ng)*R0*nhi(N+Ng)*opa_pf(N+Ng)
	if (thereis_He) then

		N15(N+Ng) = dr_j(N+Ng)*R0*nhei(N+Ng)*opa_pf(N+Ng)
		N2(N+Ng)  = dr_j(N+Ng)*R0*nheii(N+Ng)*opa_pf(N+Ng)
		if(thereis_HeITR) NTR(N+Ng) = dr_j(N+Ng)*R0*nheiTR(N+Ng)*opa_pf(N+Ng)

	endif

	do j = N+Ng-1,1-Ng,-1

	    ! Spacing
	    dr = dr_j(j)*R0*opa_pf(j)

	      ! Evaluate new column densities by integration
	    N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI

		if (thereis_He) then
			N15(j) = N15(j+1) + nhei(j)*dr	  					  ! HeI
			N2(j)  = N2(j+1)  + nheii(j)*dr       				  ! HeII
			if(thereis_HeITR) NTR(j) = NTR(j+1) + nheiTR(j)*dr	  ! HeI triplet
		endif

	enddo
	
	! End of subroutine
	end subroutine calc_column_dens

	! ------------------------------------------------------!

	subroutine calc_column_dens_one(nsp, Ncol)
	! Column density of a single species (same rectangle rule and
	! opa_pf weighting as calc_column_dens).  Used for N_H2 (molecular).
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nsp
	real*8, dimension(1-Ng:N+Ng), intent(out) :: Ncol
	Ncol = 0.0
	Ncol(N+Ng) = dr_j(N+Ng)*R0*nsp(N+Ng)*opa_pf(N+Ng)
	do j = N+Ng-1,1-Ng,-1
		Ncol(j) = Ncol(j+1) + nsp(j)*dr_j(j)*R0*opa_pf(j)
	enddo
	end subroutine calc_column_dens_one

	! ------------------------------------------------------!

	subroutine calc_column_dens_metals(nm, Nm_col)
	! Column densities for metals (same scheme as calc_column_dens).
	! nm(:,1:n_mion) holds the metal ion densities in canonical
	! species_table order; Nm_col(:,1:n_mphot) returns the column
	! density of each photo-ionizable metal ion in iphot order.

	integer :: j,i,k
	real*8, dimension(1-Ng:N+Ng,n_mion),  intent(in)  :: nm
	real*8 :: dr
	real*8, dimension(1-Ng:N+Ng,n_mphot), intent(out) :: Nm_col

	Nm_col = 0.0

	! Outer point (opa_pf weights the opacity for the 'P' model)
	do i = 1,n_mion
	   if (.not. mion_isphot(i)) cycle
	   k = mion_iphot(i)
	   Nm_col(N+Ng,k) = dr_j(N+Ng)*R0*nm(N+Ng,i)*opa_pf(N+Ng)
	enddo

	do j = N+Ng-1,1-Ng,-1
	   dr = dr_j(j)*R0*opa_pf(j)
	   do i = 1,n_mion
	      if (.not. mion_isphot(i)) cycle
	      k = mion_iphot(i)
	      Nm_col(j,k) = Nm_col(j+1,k) + nm(j,i)*dr
	   enddo
	enddo

	end subroutine calc_column_dens_metals

	! ------------------------------------------------------!

	subroutine calc_mmw(nh,nhe,ne,mmw)
	! Calculate the mean molecular weight for a certain ionization profile

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nh,nhe,ne
	real*8, dimension(1-Ng:N+Ng), intent(out) :: mmw
	
	if (thereis_He) then
		mmw = (nh + 4.0*nhe)/(nh + nhe + ne)
	else
		mmw = nh/(nh + ne)
	endif

	! End of subroutine
	end subroutine calc_mmw

	! End of module
	end module utils 

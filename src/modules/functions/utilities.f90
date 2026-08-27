   module utils
	! Collection of auxiliary subroutines

   use global_parameters
   use species_table, only: n_mion, n_mphot, mion_isphot, mion_iphot,  &
	                         mion_stage, mion_elem, melem_A,  &
	                         bsp_charge, bsp_mass

   implicit none

	contains

	! ------------------------------------------------------!

	subroutine calc_ne(nhii,nheii,nheiii,ne,nm,nmol)
	! Calculate the free electron density.
	! The optional nm (metal densities for each ion, same units as nhii) adds
	! the metal electrons when the eos_metals policy is on; omitting it
	! (or eos_metals 0) reproduces the legacy H/He-only electron count.
	! The base H/He/molecular electrons are accumulated in the canonical bsp
	! order of species_table, each species weighted by its net charge
	! (bsp_charge = free electrons released); neutral species (charge 0) are
	! skipped as no-ops.  Accumulating in bsp order deliberately fixes the FP
	! add order (a golden re-snapshot decision, section 5.3 Inc 1).

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nhii
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nheii,nheiii
	! Optional molecular ions: cols 1 H2 (neutral), 2 H2+, 3 H3+,
	! 4 HeH+ -- each molecular ion carries one electron.
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: ne

	ne = 0.0d0

	! Base H/He electrons (bsp 2,4,5; bsp 1 HI and 3 HeI are neutral, and
	! bsp 6 HeTR is both neutral and an excited level of HeI, already inside
	! the HeI column -- see bsp_is_excited_level in species_table).
	call accum(nhii,   dble(bsp_charge(2)))        ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nheii,  dble(bsp_charge(4)))     ! HeII  (bsp 4)
		call accum(nheiii, dble(bsp_charge(5)))     ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			if (mion_stage(im) .gt. 0)                                  &
				ne = ne + dble(mion_stage(im))*nm(:,im)
		enddo
	endif

	! Molecular ions (bsp 8,9,10; bsp 7 H2 is neutral).
	if (present(nmol)) then
		call accum(nmol(:,2), dble(bsp_charge(8)))     ! H2+  (bsp 8)
		call accum(nmol(:,3), dble(bsp_charge(9)))     ! H3+  (bsp 9)
		call accum(nmol(:,4), dble(bsp_charge(10)))    ! HeH+ (bsp 10)
	endif

	contains
		subroutine accum(vec, w)
		! Accumulate w*vec into ne in canonical bsp order.
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		real*8, intent(in) :: w
		ne = ne + w*vec
		end subroutine accum

	end subroutine calc_ne
	
	! ------------------------------------------------------!

	subroutine calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol)
	! Calculate the total atomic number density.
	! Every species counts as ONE gas particle, so the base H/He/molecular
	! contribution is accumulated with unit weight in the canonical bsp order
	! of species_table; the optional nm adds the metal nuclei (all stages)
	! when the eos_metals policy is on.  bsp-order accumulation deliberately
	! fixes the FP add order (a golden re-snapshot decision, section 5.3 Inc 1).
	! The He 2^3S column is NOT an argument: it is an excited level of He I
	! (bsp_is_excited_level), so its gas particle is already the He I particle
	! counted through nhei, and adding it counted the triplet twice.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_tot

	n_tot = 0.0d0

	! Base H/He particles (bsp 1..5), one particle each; bsp 6 (HeTR) is an
	! excited level of bsp 3 (HeI) and is already inside it.
	call accum(nhi)                     ! HI    (bsp 1)
	call accum(nhii)                    ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nhei)                ! HeI   (bsp 3, He 2^3S included)
		call accum(nheii)               ! HeII  (bsp 4)
		call accum(nheiii)              ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_tot = n_tot + nm(:,im)
		enddo
	endif

	! Each molecule is ONE gas particle (bsp 7..10).
	if (present(nmol)) then
		call accum(nmol(:,1))           ! H2   (bsp 7)
		call accum(nmol(:,2))           ! H2+  (bsp 8)
		call accum(nmol(:,3))           ! H3+  (bsp 9)
		call accum(nmol(:,4))           ! HeH+ (bsp 10)
	endif

	contains
		subroutine accum(vec)
		! Accumulate vec into n_tot (unit weight; one particle per species).
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		n_tot = n_tot + vec
		end subroutine accum

	end subroutine calc_ntot

	! ------------------------------------------------------!
	
	subroutine calc_rho(nhi,nhii,nhei,nheii,nheiii,n_out,nm,nmol)
	! Calculate the total mass density (adimensional).
	! The base H/He/molecular mass is accumulated in the canonical bsp order
	! of species_table, each species weighted by bsp_mass [m_H units]; this
	! bsp-order weighted accumulation deliberately reorders the FP adds versus
	! the old factored 4.0*(nhei+...) form (a golden re-snapshot decision,
	! section 5.3 Inc 1).  The optional nm adds the metal mass (melem_A per
	! nucleus, all stages) when the eos_metals policy is on.  (HeH+ carries
	! 5 m_H: its He nucleus is NOT in the nhei..nheiii free-He arrays.)
	! The He 2^3S column is NOT an argument: it is an excited level of He I
	! (bsp_is_excited_level), so its 4 m_H are already the He I atom's mass
	! counted through nhei, and adding it put the triplet mass in rho twice.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_out

	n_out = 0.0d0

	! Base H/He mass (bsp 1..5), weighted by bsp_mass; bsp 6 (HeTR) is an
	! excited level of bsp 3 (HeI) and is already inside it.
	call accum(nhi,    bsp_mass(1))     ! HI    (bsp 1)
	call accum(nhii,   bsp_mass(2))     ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nhei,   bsp_mass(3))     ! HeI   (bsp 3, He 2^3S included)
		call accum(nheii,  bsp_mass(4))     ! HeII  (bsp 4)
		call accum(nheiii, bsp_mass(5))     ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_out = n_out + melem_A(mion_elem(im))*nm(:,im)
		enddo
	endif

	! Molecular mass (bsp 7..10): H2/H2+ = 2, H3+ = 3, HeH+ = 5 m_H.
	if (present(nmol)) then
		call accum(nmol(:,1), bsp_mass(7))     ! H2   (bsp 7)
		call accum(nmol(:,2), bsp_mass(8))     ! H2+  (bsp 8)
		call accum(nmol(:,3), bsp_mass(9))     ! H3+  (bsp 9)
		call accum(nmol(:,4), bsp_mass(10))    ! HeH+ (bsp 10)
	endif

	contains
		subroutine accum(vec, w)
		! Accumulate w*vec into n_out in canonical bsp order.
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		real*8, intent(in) :: w
		n_out = n_out + w*vec
		end subroutine accum

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

	subroutine calc_mmw(nh,nhe,ne,mmw,nm)
	! Calculate the mean molecular weight for a certain ionization profile.
	! The H/He nucleus masses come from the species_table metadata
	! (bsp_mass(1) = HI = 1, bsp_mass(3) = HeI = 4 in m_H units — for the
	! atomic species the bsp position equals the f_sp column), reproducing
	! the old literals bitwise. The optional nm adds the metal mass and
	! metal nuclei under the same eos_metals policy as calc_rho/calc_ntot,
	! so the post-process temperature solve uses the same composition as
	! the main loop (the metal mass raises mmw by ~1%; omitting nm keeps
	! the legacy H/He-only diagnostic).

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nh,nhe,ne
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: mmw
	real*8, dimension(1-Ng:N+Ng) :: mass_l, npart_l

	if (thereis_He) then
		mass_l  = bsp_mass(1)*nh + bsp_mass(3)*nhe
		npart_l = nh + nhe + ne
	else
		mass_l  = nh
		npart_l = nh + ne
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			mass_l  = mass_l + melem_A(mion_elem(im))*nm(:,im)
			npart_l = npart_l + nm(:,im)
		enddo
	endif

	mmw = mass_l/npart_l

	! End of subroutine
	end subroutine calc_mmw

	! End of module
	end module utils 

   module energy_vectors_construct
   ! Construct energy grid for the computation of radiative
	!      contributions and radiative equilibrium
      
   use global_parameters
   use species_table, only: n_mphot, n_mion, mion_ethr,                 &
                            mion_isphot, mion_elem, mion_iphot
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

   ! Bound of the active-threshold list ionization_thresholds_active
   ! returns: the fixed absorbers, + 1 for e_th_H2_dd (not in eth_abs), and
   ! for every metal ion its ionization threshold plus one turn-on per
   ! subshell whose K or L edge the grid can reach (metal_subshell_turn_on).
   ! Callers size their e_thr by this.
   integer, parameter, public :: n_thr_max = n_abs_fixed + 1               &
                                           + n_mion*(1 + n_metal_subshell)

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
	! Bins of the two bands below the H I edge: [e_abs_low, 13.6 eV], where
	! the metastable helium and the low-IP metals absorb, and
	! [e_th_HI_n2, e_abs_low], where the only absorber is hydrogen in n = 2.
	! The 89 bins of the first give it the bin ratio of the pure-H I band
	! (ln(13.598/4.768)/89 = ln(24.587/13.598)/50 = 0.0118): the neutral
	! metal cross sections there vary on a small fraction of their own
	! threshold energy (Verner et al. 1996 fit of Fe I: E_0 = 0.0546 eV, a
	! factor 2.5 over 0.4 eV), and with this budget every metal
	! photoionization rate the grid forms is within 5e-4 of the exact
	! integral (photon_grid_threshold_edges; docs/Update_EXHALE.md section 6).
	integer, parameter :: num_TR = 89
	! The n = 2 band carries no cross section of an active absorber (that is
	! what e_abs_low means, see the split below), and the H(n=2)
	! photoionization rate is an integral of stellar_flux_eV over its own
	! band (excited_hydrogen), not a sum over these bins. The count therefore
	! only has to STATE the band; 20 bins keep its widest bin ratio at 1.072
	! (e_abs_low = 13.6 eV, the widest case).
	integer, parameter :: num_n2 = 20

	real*8, dimension(:), allocatable :: dum_E,dum_F,dum_dE
	real*8 :: e_sub_low    ! lower edge of the whole below-13.6 eV sub-grid
	real*8 :: e_abs_low    ! lowest threshold of a bulk-gas absorber [eV]
	! The bands the grid is built from: their boundaries, their bin budgets,
	! and the number of them this configuration has.
	integer, parameter :: n_band_max = 6
	real*8  :: e_band_lo(n_band_max), e_band_hi(n_band_max)
	integer :: n_band_bin(n_band_max)
	integer :: n_band, ib, k_band
	! The ionization thresholds of the absorbers this run carries, the
	! ascending boundaries of the sub-intervals one band is cut into at those
	! of them that fall inside it, and the bins each sub-interval gets.
	! n_thr_max, the bound of that list, is a module parameter (above
	! contains) so that every caller of ionization_thresholds_active sizes
	! its array by the same number.
	real*8  :: e_thr(n_thr_max)
	real*8  :: e_bnd(n_thr_max + 2)
	integer :: n_sub_bin(n_thr_max + 1)
	integer :: n_thr, n_sub
	! The Nl+1 bin edges of the photon grid, and the index of the last edge
	! written so far.
	real*8, dimension(:), allocatable :: e_edge
	integer :: k_edge
	! The four H2 channel cross sections at one energy, and their ionizing
	! part (channels M + S + D).
	real*8 :: sig4(n_h2_channels), sig_ion

	! Set the energy band of the spectrum.
	!
	! The power law and the photospheric blackbody are both ANALYTIC over
	! the whole band, so they are evaluated on the same grid: the grid is a
	! property of the absorbers (their thresholds and the resolution the
	! quadrature needs), not of the spectrum. Which of the two fills F_XUV
	! is decided below, where the flux is put on the grid. One spectrum type
	! builds every band (development plan rev 3, section 10.5 decision 13).
	if (is_PL_sed .or. spectrum_is_planck()) then

		! Points below 13.6 eV. Three absorbers need them, and they can be
		! active together (the merged triplet+metals systems):
		!  - He triplet (HeI 2^3S), threshold 4.768 eV;
		!  - low-IP metals (e.g. Mg I, threshold 7.646 eV; K I, 4.341 eV);
		!  - hydrogen in n = 2, threshold e_th_HI_n2 = 3.400 eV, whenever
		!    the excited-hydrogen coupling is armed.
		! The grid floor is the LOWEST threshold over whichever are active,
		! by the one rule the loaded-SED path uses as well (sed_read's
		! photon_grid_floor_eV): stopping at the triplet edge would leave
		! K I, whose threshold is below it, without ionizing photons and
		! hold it spuriously neutral.
		!
		! The band is split at e_abs_low, the lowest threshold of an
		! absorber of the bulk gas, because that threshold is a bin edge
		! like every other (see the partition below): below it no active
		! absorber has a cross section (sub_lyman_absorber_threshold_eV
		! in sed_read states this exactly).
		e_abs_low = sub_lyman_absorber_threshold_eV()
		e_sub_low = photon_grid_floor_eV()

		! --- Construct energy grid --- !
		!
		! EVERY IONIZATION THRESHOLD OF AN ACTIVE ABSORBER IS A BIN EDGE.
		! A photoionization rate is the rectangle rule sum(F sigma/E de_v)
		! over the bins, so a bin that CONTAINS a threshold charges that
		! absorber's cross section, read at the bin centre, over the part of
		! the bin that lies BELOW the threshold as well, where the absorber
		! cannot absorb at all. The partition is therefore cut at every
		! threshold the run carries -- H I, He I, He II, the He 2^3S
		! metastable when it is tracked, H2 when the molecular chemistry is
		! on, and every photo-ionizable stage of an active metal -- and e_v
		! holds the bin CENTRES, which then lie strictly inside one interval
		! of continuity of every cross section the sum contracts with
		! (development plan rev 3, section 10.1 item 2; the audit of
		! 2026-09-05 measured the old grid, whose points sat ON the H and He
		! thresholds with central-difference widths, at P_HI 1.095x, P_HeI
		! 1.016x and P_HeII 1.030x the exact integral).
		!
		! HOW THE BINS ARE SHARED OUT. The grid keeps the six BANDS it has
		! always had -- [floor, e_abs_low], [e_abs_low, 13.6], [13.6, 24.6],
		! [24.6, 54.4], [54.4, e_mid], [e_mid, e_top] -- and each keeps its
		! own budget of bins (num_n2, num_TR, num_HI, num_HeI, num_HeII,
		! num_X), because that is what puts the resolution where the
		! absorbers are rather than spreading it over the X-ray decades. A
		! band is then cut at every threshold strictly inside it, and its
		! budget is shared out over the resulting sub-intervals IN PROPORTION
		! TO THEIR LOGARITHMIC WIDTHS, at least one bin each (largest
		! remainder, so the budget is met exactly). A band with no threshold
		! inside it is left with its budget in one piece, i.e. exactly the
		! partition it had before, and the total stays Nl_fix + NlTR. The one
		! case that could raise the total is a band cut into more
		! sub-intervals than its budget, which no configuration of this
		! species table reaches (the widest, [e_abs_low, 13.6 eV], holds at
		! most seven of them against a budget of twenty).
		!
		! The bins are contiguous and de_v is the exact bin width, so
		! sum(de_v) = e_top - e_sub_low exactly: the grid covers the whole
		! band [e_low, e_top] that J_inc is normalized over, the X-ray top
		! included.
		!
		! The centre of a bin is the GEOMETRIC mean of its two edges. The
		! edges are geometrically spaced and the integrand F(E) sigma(E)/E
		! is close to a power law in E, for which the geometric abscissa of
		! a rectangle rule is the more accurate of the two choices: on the
		! widest bins of the grid, the sub-Lyman ones (edge ratio up to
		! 1.078, the one-bin interval [4.768, 5.139] eV between the
		! metastable and Na I in a solar-metal run), an E^-4 integrand is
		! reproduced to 0.09% by the geometric centre against 0.23% by the
		! arithmetic one.
		call ionization_thresholds_active(e_thr, n_thr)

		n_band = 0
		! Band between the n=2 threshold and the lowest bulk-gas absorber
		! (H(n=2) is the only absorber of this band)
		if (e_sub_low .lt. e_abs_low)                                   &
			call append_band(e_sub_low, e_abs_low, num_n2,              &
			                 e_band_lo, e_band_hi, n_band_bin, n_band)
		! Band between e_abs_low and 13.6 eV (HeI triplet, low-IP metals)
		if (e_abs_low .lt. e_th_HI)                                     &
			call append_band(e_abs_low, e_th_HI,  num_TR,               &
			                 e_band_lo, e_band_hi, n_band_bin, n_band)
		! Pure-HI band [13.6 eV, 24.6 eV]
		call append_band(e_th_HI,   e_th_HeI,  num_HI,                  &
		                 e_band_lo, e_band_hi, n_band_bin, n_band)
		! HI+HeI band [24.6 eV, 54.4 eV]
		call append_band(e_th_HeI,  e_th_HeII, num_HeI,                 &
		                 e_band_lo, e_band_hi, n_band_bin, n_band)
		! HI+HeI+HeII band [54.4 eV, 124 eV]
		call append_band(e_th_HeII, e_mid,     num_HeII,                &
		                 e_band_lo, e_band_hi, n_band_bin, n_band)
		! X-ray band [124 eV, e_XUV_top]
		call append_band(e_mid,     e_top,     num_X,                   &
		                 e_band_lo, e_band_hi, n_band_bin, n_band)

		! The edges, band by band. NlTR counts every bin below the H I edge,
		! which is what write_setup_report reads it for.
		allocate(e_edge(Nl_fix + num_TR + num_n2 + n_thr_max + 1))
		e_edge(1) = e_band_lo(1)
		k_edge    = 1
		NlTR      = 0
		do ib = 1,n_band
			call band_sub_boundaries(e_band_lo(ib), e_band_hi(ib),      &
			                         e_thr, n_thr, e_bnd, n_sub)
			call bins_by_log_width(e_bnd, n_sub, n_band_bin(ib),        &
			                       n_sub_bin)
			k_band = k_edge
			do j = 1,n_sub
				call geometric_band_edges(e_bnd(j), e_bnd(j+1),         &
				                          n_sub_bin(j), e_edge, k_edge)
			enddo
			if (e_band_hi(ib) .le. e_th_HI) NlTR = NlTR + k_edge - k_band
		enddo
		! k_edge is the index of the LAST edge written, so the partition
		! carries k_edge - 1 bins.
		Nl = k_edge - 1

		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(F_XUV(Nl))

		! Bin centres and exact bin widths
		do j = 1,Nl
			e_v(j)  = sqrt(e_edge(j)*e_edge(j+1))
			de_v(j) = e_edge(j+1) - e_edge(j)
		enddo
		deallocate(e_edge)

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
	 
		! Read the SED file. The thresholds of the active absorbers go with
		! it: the loaded grid is cut at them exactly as the analytic one is,
		! and this is the one list (below).
		call ionization_thresholds_active(e_thr, n_thr)
	 	call read_sed(e_thr, n_thr)
 	
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
	! (its header, and Update_EXHALE_stage1 section 150).
	xi_day = dayside_dilution()
	if (is_PL_sed) then

		F_XUV  = (/ (xi_day*J_inc(e_v(i)), i = 1,Nl) /)

	else if (spectrum_is_planck()) then

		! The photospheric field pi B_nu(T_eff) (R_star/a)^2 per unit
		! photon energy, on every point of the grid -- the XUV and the
		! sub-13.6 eV band alike (decision 13). Its absolute scale comes
		! from T_eff and R_star, not from LX/LEUV, so the integrated grid
		! flux differs from the nominal J_XUV; write_setup_report states
		! both. a_orb is in cm here (input_read), as R_star is.
		F_XUV  = (/ (xi_day*planck_stellar_flux_eV(T_star_eff,          &
		             R_star/max(a_orb,1.0d-30), e_v(i)), i = 1,Nl) /)

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

    ! --------------------------------------------------

    subroutine ionization_thresholds_active(e_thr, n_thr)
    ! The ionization thresholds [eV] of every absorber the run carries: the
    ! hydrogen and helium stages, the He 2^3S metastable when it is tracked,
    ! H2 when the molecular chemistry is on, and every photo-ionizable stage
    ! of an element the run has a nonzero abundance for.  These are the
    ! energies at which a cross section of the photoionization sums turns on,
    ! and therefore the energies the photon grid has to be cut at.
    !
    ! A metal ion contributes more than one energy.  Its photoabsorption is
    ! the sum over the shells (metal_photoion_sigma), so it steps again at
    ! every subshell turn-on above the ionization threshold: the K edges of
    ! C, N and O at 291, 404.8 and 538 eV, the L edges of Na, Mg, Si, S, K,
    ! Ca and Fe, and the E_max at which a valence subshell that the 1996
    ! outer-shell fit represented below it is allowed to appear on its own.
    ! E_max is the same energy at which the outer shell itself changes fit,
    ! so cutting there also puts that step on an edge.
    !
    ! An inactive absorber is left out on purpose: the cross-section table
    ! sigma_tab is filled for all seventeen photo-ionizable metal ions
    ! whatever metals.inp carries, but a column of an element the run does
    ! not have is never contracted with anything, so cutting the grid at its
    ! threshold would only spend bins.
    !
    ! The list is unordered and may repeat an energy; band_sub_boundaries
    ! sorts and deduplicates whatever falls inside a band.
    real*8,  intent(out) :: e_thr(:)
    integer, intent(out) :: n_thr

    integer :: i, m, n_on
    real*8  :: e_on(n_metal_subshell)

    n_thr    = 3
    e_thr(1) = e_th_HI
    e_thr(2) = e_th_HeI
    e_thr(3) = e_th_HeII
    if (thereis_HeITR) then
       n_thr        = n_thr + 1
       e_thr(n_thr) = e_th_HeTR
    endif
    if (thereis_mol) then
       n_thr        = n_thr + 1
       e_thr(n_thr) = e_th_H2
       ! The H2 channel thresholds are turn-ons as well: above e_th_H2_di
       ! the dissociative-ionization share and above e_th_H2_dd the
       ! double-ionization share of one absorption cross section switch on
       ! (h2_photo_channels), so a bin straddling either charges part of
       ! the bin to a channel that is closed there.  The neutral
       ! dissociation window 32-41.5 eV is not cut: its share is anchored
       ! to zero at both ends and is continuous, so a straddling bin has a
       ! kink, not a step, in its integrand.
       n_thr        = n_thr + 1
       e_thr(n_thr) = e_th_H2_di
       n_thr        = n_thr + 1
       e_thr(n_thr) = e_th_H2_dd
    endif
    if (allocated(melem_ab)) then
       do i = 1,n_mion
          if (.not. mion_isphot(i))              cycle
          if (mion_ethr(i) .le. 0.0d0)           cycle
          if (melem_ab(mion_elem(i)) .le. 0.0d0) cycle
          n_thr        = n_thr + 1
          e_thr(n_thr) = mion_ethr(i)
          call metal_subshell_turn_on(mion_iphot(i), e_on, n_on)
          do m = 1,n_on
             n_thr        = n_thr + 1
             e_thr(n_thr) = e_on(m)
          enddo
       enddo
    endif

    end subroutine ionization_thresholds_active

    ! --------------------------------------------------

    subroutine append_band(e_lo, e_hi, nbin, e_band_lo, e_band_hi,          &
                           n_band_bin, n_band)
    ! Append one band of the photon grid, with the budget of bins it is to
    ! be partitioned into, to the ascending band list.
    real*8,  intent(in)    :: e_lo, e_hi
    integer, intent(in)    :: nbin
    real*8,  intent(inout) :: e_band_lo(:), e_band_hi(:)
    integer, intent(inout) :: n_band_bin(:)
    integer, intent(inout) :: n_band

    n_band             = n_band + 1
    e_band_lo(n_band)  = e_lo
    e_band_hi(n_band)  = e_hi
    n_band_bin(n_band) = nbin

    end subroutine append_band

    ! --------------------------------------------------

    subroutine band_sub_boundaries(e_lo, e_hi, e_thr, n_thr, e_bnd, n_sub)
    ! Cut the band [e_lo, e_hi] at every ionization threshold that lies
    ! strictly inside it, and return the n_sub+1 ascending boundaries of the
    ! sub-intervals.
    !
    ! rel_same is the relative separation below which two energies are one
    ! partition point.  Two thresholds closer than that -- or a threshold
    ! that is already a band boundary -- would otherwise leave a bin of zero
    ! or of round-off width, whose geometric centre is the threshold itself
    ! and whose weight is meaningless.  1e-10 is far below the spacing of
    ! any two thresholds of the species table (the closest pair is Fe I
    ! 7.902 eV and Mg I 7.646 eV, 3 per cent apart) and far above the
    ! round-off of the constants.
    real*8,  intent(in)  :: e_lo, e_hi
    real*8,  intent(in)  :: e_thr(:)
    integer, intent(in)  :: n_thr
    real*8,  intent(out) :: e_bnd(:)
    integer, intent(out) :: n_sub

    real*8, parameter :: rel_same = 1.0d-10
    integer :: i, j, m
    real*8  :: e
    logical :: is_same

    m        = 1
    e_bnd(1) = e_lo
    do i = 1,n_thr
       e = e_thr(i)
       if (e .le. e_lo*(1.0d0 + rel_same)) cycle
       if (e .ge. e_hi*(1.0d0 - rel_same)) cycle
       is_same = .false.
       do j = 1,m
          if (abs(e - e_bnd(j)) .le. rel_same*e) is_same = .true.
       enddo
       if (is_same) cycle
       ! insert, keeping the list ascending
       j = m
       do while (j .ge. 2)
          if (e_bnd(j) .le. e) exit
          e_bnd(j+1) = e_bnd(j)
          j = j - 1
       enddo
       e_bnd(j+1) = e
       m = m + 1
    enddo
    e_bnd(m+1) = e_hi
    n_sub      = m

    end subroutine band_sub_boundaries

    ! --------------------------------------------------

    subroutine bins_by_log_width(e_bnd, n_sub, n_budget, n_bin)
    ! Share a band's budget of bins out over its n_sub sub-intervals in
    ! proportion to their LOGARITHMIC widths, at least one bin each.
    !
    ! Logarithmic width is the measure because the bins are geometric: equal
    ! shares of ln E give every bin of the band the same edge ratio, so the
    ! resolution of the quadrature is uniform across the band whatever the
    ! thresholds cut out of it.  The exact budget is met by the largest
    ! remainder: floor the proportional shares, then hand the leftover bins
    ! to the sub-intervals with the largest fractional deficit (and take
    ! them back from the smallest when the one-bin minimum overshoots).
    !
    ! A band cut into more sub-intervals than its budget gets one bin each
    ! and so more bins than the budget; the caller's comment states that no
    ! configuration of the present species table reaches that.
    real*8,  intent(in)  :: e_bnd(:)
    integer, intent(in)  :: n_sub, n_budget
    integer, intent(out) :: n_bin(:)

    integer :: i, i_pick, n_tot, n_use
    real*8  :: w(n_sub), q(n_sub), w_tot, dev, dev_pick

    n_use = max(n_budget, n_sub)
    w_tot = 0.0d0
    do i = 1,n_sub
       w(i)  = log(e_bnd(i+1)/e_bnd(i))
       w_tot = w_tot + w(i)
    enddo

    n_tot = 0
    do i = 1,n_sub
       q(i)     = dble(n_use)*w(i)/w_tot
       n_bin(i) = max(1, int(q(i)))
       n_tot    = n_tot + n_bin(i)
    enddo

    do while (n_tot .lt. n_use)
       i_pick   = 1
       dev_pick = -huge(1.0d0)
       do i = 1,n_sub
          dev = q(i) - dble(n_bin(i))
          if (dev .gt. dev_pick) then
             dev_pick = dev
             i_pick   = i
          endif
       enddo
       n_bin(i_pick) = n_bin(i_pick) + 1
       n_tot         = n_tot + 1
    enddo

    do while (n_tot .gt. n_use)
       i_pick   = 0
       dev_pick = huge(1.0d0)
       do i = 1,n_sub
          if (n_bin(i) .le. 1) cycle
          dev = q(i) - dble(n_bin(i))
          if (dev .lt. dev_pick) then
             dev_pick = dev
             i_pick   = i
          endif
       enddo
       if (i_pick .eq. 0) exit
       n_bin(i_pick) = n_bin(i_pick) - 1
       n_tot         = n_tot - 1
    enddo

    end subroutine bins_by_log_width

    ! --------------------------------------------------

    subroutine geometric_band_edges(e_lo, e_hi, nbin, e_edge, k_edge)
    ! Append the UPPER edges of nbin geometrically spaced bins that partition
    ! the band [e_lo, e_hi] to the edge list e_edge, whose entry k_edge is
    ! already e_lo, and advance k_edge to the last one written.
    !
    ! The spacing is geometric because the photon field and the
    ! photoionization cross sections are power laws over decades of energy,
    ! so a constant ratio between neighbouring edges spends the bins where
    ! the integrand varies. The top edge is set to e_hi itself rather than
    ! to e_lo (e_hi/e_lo)^(nbin/nbin): the band boundary is an ionization
    ! threshold and must be exact, not the round-off of a power.
    real*8,  intent(in)    :: e_lo, e_hi
    integer, intent(in)    :: nbin
    real*8,  intent(inout) :: e_edge(:)
    integer, intent(inout) :: k_edge

    integer :: j

    do j = 1,nbin-1
       e_edge(k_edge+j) = e_lo*(e_hi/e_lo)**(dble(j)/dble(nbin))
    enddo
    e_edge(k_edge+nbin) = e_hi
    k_edge = k_edge + nbin

    end subroutine geometric_band_edges


    ! End of module
    end module energy_vectors_construct


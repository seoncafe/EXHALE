	module sed_reader
	! Module containing the subroutine for reading the numerical SED
		
	use global_parameters
	use species_table, only: n_melem, mion_ethr, melem_i0, melem_name
	! The n = 2 ionization threshold of hydrogen, one of the absorbers the
	! photon-grid floor below is the lowest threshold over.
	use J_incident, only: e_th_HI_n2

	implicit none

	contains

	!---------------------------------------------------!

	double precision function sub_lyman_absorber_threshold_eV()
	! Lowest ionization threshold [eV] of an active absorber whose threshold
	! lies below the H I edge, OTHER than hydrogen in n = 2:
	!   the He 2^3S metastable, e_th_HeTR = 4.768 eV = 2600 A, and
	!   the neutral stage of any active low-IP metal (K I 4.341 eV, Na I
	!   5.139 eV, Ca I 6.113 eV, Mg I 7.646 eV, Si I 8.152 eV, Fe I 7.902
	!   eV, S I 10.36 eV).
	! e_th_HI when none of them is active, i.e. when the ionizing band
	! starts at the H I edge.
	!
	! This is the energy where the ABSORBERS OF THE BULK GAS turn on, and it
	! is a bin edge of the photon grid: below it no active absorber has a
	! cross section. (set_energy_vectors fills the metal cross-section table
	! for every photo-ionizable ion whatever metals.inp carries, so entries
	! below this energy do exist -- K I at 4.341 eV, say -- but they belong
	! to an element the run does not carry; were it carried, its threshold
	! would BE this energy.)
	integer :: j

	sub_lyman_absorber_threshold_eV = e_th_HI
	if (thereis_HeITR)                                                    &
		sub_lyman_absorber_threshold_eV = e_th_HeTR
	if (thereis_lowIP_metal .and. allocated(melem_ab)) then
		do j = 1,n_melem
			if (melem_ab(j) .gt. 0.0d0 .and.                            &
			    mion_ethr(melem_i0(j)) .lt. e_th_HI)                    &
				sub_lyman_absorber_threshold_eV =                       &
				   min(sub_lyman_absorber_threshold_eV,                 &
				       mion_ethr(melem_i0(j)))
		enddo
	endif

	end function sub_lyman_absorber_threshold_eV

	!---------------------------------------------------!

	double precision function photon_grid_floor_eV()
	! Lowest photon energy [eV] the run's photon grid has to reach: the
	! lowest ionization threshold over every active absorber, hydrogen in
	! n = 2 included.
	!
	! H(n=2) IS AN ABSORBER LIKE THE OTHERS. Its threshold e_th_HI_n2 =
	! 3.400 eV = 3647 A is the head of the Balmer continuum, and the
	! excited-hydrogen coupling (armed by "Stellar Teff [K]:" and "Stellar
	! radius [R_sun]:", input_read use_excited_H) photoionizes the n = 2
	! population out of the stellar field there. So its threshold floors the
	! grid exactly as the metastable's and the low-IP metals' do: a loaded
	! table is read down to it, and the analytic types are evaluated down to
	! it (development plan rev 3, section 10.5 decisions 13 and 17).
	!
	! e_th_HI when no absorber below the H I edge is active.

	photon_grid_floor_eV = sub_lyman_absorber_threshold_eV()
	if (use_excited_H)                                                    &
		photon_grid_floor_eV = min(photon_grid_floor_eV, e_th_HI_n2)

	end function photon_grid_floor_eV

	!---------------------------------------------------!

	subroutine read_sed(e_thr, n_thr)
	! Subroutine to read the numerical SED from file, and to lay the photon
	! grid out on it.
	!	The file has to be formatted in two columns:
	!	   1) bin central wavelength [Angstrom] 
	!	   2) flux at the planet surface [erg cm^{-2} s^{-2} A^{-1}]
	!	Data have to be order with increasing wavelength
	!
	! THE GRID IS BUILT FROM BIN EDGES, AND EVERY IONIZATION THRESHOLD OF AN
	! ACTIVE ABSORBER IS ONE OF THEM. A photoionization rate is the rectangle
	! rule sum(F sigma/E de_v) over the bins, and a cross section turns ON at
	! its threshold: a bin that CONTAINS a threshold charges the absorber's
	! cross section, read at the bin centre, over the part of the bin lying
	! BELOW the threshold as well, where the absorber cannot absorb at all.
	! This is the same reasoning the analytic types' grid is built on
	! (set_energy_vectors; development plan rev 3 section 10.1 item 2), where
	! the audit of 2026-09-05 measured the straddling error at P_HI 1.095x,
	! P_HeI 1.016x and P_HeII 1.030x the exact integral.
	!
	! The layout, from the selected rows E_1 < E_2 < ... < E_N (ascending in
	! energy; the file is ascending in wavelength):
	!   - interior edges are the GEOMETRIC mean sqrt(E_k E_k+1) of
	!     neighbouring rows, the log-grid writing of the midpoint between two
	!     tabulated points;
	!   - the two end edges are E_1 and E_N themselves, so the grid spans
	!     EXACTLY the band the table states a field on and states nothing
	!     outside it. sum(de_v) = E_N - E_1;
	!   - each row's flux is the flux of its bin: a tabulated spectrum is a
	!     flux histogram, one value per row;
	!   - every threshold of e_thr strictly inside the span is inserted as a
	!     further edge, splitting the bin it falls in into two halves that
	!     both carry that row's flux. The split is therefore EXACT for the
	!     integrated flux, sum(F_XUV de_v) being unchanged by it, and each
	!     half gets its own geometric centre and its own exact width. Edges
	!     closer than rel_same in relative separation are one edge.
	! e_v is the bin CENTRE and de_v the exact bin width; the table's own
	! rows are kept separately, in e_sed_node/F_sed_node, for the consumers
	! that read the FIELD off the file rather than integrate over the grid
	! (stellar_flux_eV, loaded_table_floor_eV).
	!
	! Ordering: e_v, de_v and F_XUV are left DECREASING in energy, the order
	! of the file, and set_energy_vectors flips them to ascending as it has
	! always done. de_v > 0 throughout, and every consumer sums over Nl, so
	! nothing downstream depends on the order.

	! The ionization thresholds [eV] of the absorbers this run carries. The
	! caller supplies them (energy_vectors_construct's
	! ionization_thresholds_active), which is where the one list lives; the
	! list may be unordered and may repeat an energy.
	real*8,  intent(in) :: e_thr(:)
	integer, intent(in) :: n_thr

	! Two energies closer than this in relative separation are one edge: a
	! bin of round-off width has a meaningless weight and a geometric centre
	! equal to its own edges. The same value band_sub_boundaries uses.
	real*8, parameter :: rel_same = 1.0d-10

	integer :: io
	integer :: i,j,k,m,jj,skip
	integer :: n_node        ! rows of the table selected
	integer :: n_bin         ! bins of the photon grid after the splits
	logical :: is_dup
	real*8 :: e_lo_k,e_hi_k,e_prev,e_swap
	real*8, dimension(:), allocatable :: f_lam   ! F_lambda of the rows
	real*8, dimension(:), allocatable :: en,fn   ! rows, ascending in energy
	real*8, dimension(:), allocatable :: b_lo,b_hi,f_bin ! bins, ascending
	real*8, dimension(:), allocatable :: t_in    ! thresholds inside one bin
	integer :: nrow          ! running count of rows read from the SED file
	real*8 :: e_top_read
	real*8 :: dum_w,dum_f,dum_e
	real*8 :: w_prev         ! previous selected wavelength (monotonicity check)
	real*8 :: e_min_seen     ! lowest selected photon energy encountered
	real*8 :: LEUV_int,LX_int,Df_int
	real*8 :: e_thr_j        ! neutral-metal ionization threshold [eV]
	real*8 :: e_grid_low     ! lowest active ionization threshold [eV]
	integer :: n_absorber    ! absorbers left without photons above threshold
	integer :: n_metal_absorber ! of those, active neutral metals
	character(len=5) :: neutral_label ! element symbol + ionization stage
	logical :: hit_eof       ! true if the selection block ran to end of file
	real*8, dimension(:), allocatable :: wave_c

   write(*,*) '(sed_read.f90) Reading the numerical spectrum..'
   
	! Open file to read number of lines to be skipped
	!	according to the selected energy interval
	open(unit = 1, file = sed_file, status = 'old')
	
		! Lower the floor of the band that is read to the lowest
		! ionization threshold of an active absorber, so that every
		! absorber is given the photons it can absorb. The rule is the
		! same one the analytic types build their grid on
		! (set_energy_vectors), written once in photon_grid_floor_eV.
		! Only lowered: the [E_low,E_mid,E_high] line of input.inp states
		! the band the run asks for, and no absorber below 13.6 eV being
		! active leaves it as it stands.
		e_grid_low = photon_grid_floor_eV()
		if (e_grid_low .lt. e_th_HI) e_low = min(e_low, e_grid_low)

		! Set highest energy in SED
		e_top_read = e_top
		if (.not.thereis_Xray) e_top_read = e_mid 

		! Initialize counters
		Nl        = 0
		skip      = 0
		nrow      = 0
		hit_eof   = .false.
		e_min_seen = e_top_read

		! Read through file
		do
			! Read wavelength and flux
			if (.not. sed_next_row(1, dum_w, dum_f, io) .and. io .eq. 0)  &
				io = -1

			! Handle the read status before using the values: a clean
			! end of file stops the selection; any other error is fatal
			! and reports the SED file and the offending row.
			if (io .lt. 0) then
				hit_eof = .true.
				exit
			else if (io .gt. 0) then
				write(*,*) '(sed_read.f90) Read error in SED file "',   &
				           trim(sed_file), '" at row ', nrow + 1
				error stop 'read_sed: malformed SED file'
			endif
			nrow = nrow + 1

			! Convert wavelength to eV
			dum_e = hp_eV*c_light/(dum_w*1.0d-8)

			! Add one line to skip if current energy > e_top_read
			!	and keep reading
			if (dum_e .gt. e_top_read) then
				skip = skip + 1
				cycle
			endif

			! if current energy < e_low quit loop
			if (dum_e .lt. e_low) exit

			! Update number of selected points and record the lowest
			! selected energy reached so far
			Nl = Nl + 1
			e_min_seen = dum_e

		enddo

		! Require at least two selected points; fewer leaves the energy-bin
		! width and downstream integrals undefined.
		if (Nl .lt. 2) then
			write(*,*) '(sed_read.f90) SED file "', trim(sed_file),      &
			           '" has fewer than 2 points in [e_low, e_top_read] = [', &
			           e_low, ',', e_top_read, '] eV'
			error stop 'read_sed: SED does not cover the requested band'
		endif

		! The loaded SED must cover the whole photon grid [e_low, e_top_read]:
		! the grid floor e_low is the lowest ionization threshold any active
		! absorber has (H(n=2) at e_th_HI_n2 = 3.400 eV = 3647 A when the
		! excited-hydrogen coupling is armed, He 2^3S at e_th_HeTR = 4.768
		! eV = 2600 A, or the neutral stage of an active low-IP metal, in
		! photon_grid_floor_eV above), and a table that stops
		! above it leaves that absorber with no field over part of its own
		! band, so its photoionization rate is not the rate of the stated
		! spectrum but of a truncation of it. The run therefore stops; the
		! composition of the run, not the reader, decides the floor
		! (development plan 2026-09-05 rev 3, section 10.5 decision 17).
		if (hit_eof .and. e_min_seen .gt. e_low) then
			write(*,*) '(sed_read.f90) ERROR: the loaded SED stops above '// &
			           'the photon-grid floor the active absorbers need.'
			write(*,'(A,A)')  '    SED file             : ', trim(sed_file)
			write(*,'(A,ES12.5,A,ES12.5,A)')                               &
			   '    lowest photon in it  : ', e_min_seen, ' eV = ',        &
			   hp_eV*c_light*1.0d8/e_min_seen, ' A'
			write(*,'(A,ES12.5,A,ES12.5,A)')                               &
			   '    band not covered     : ', e_low, ' to ', e_min_seen,   &
			   ' eV'
			write(*,'(A,ES12.5,A,ES12.5,A)')                               &
			   '                         = ',                              &
			   hp_eV*c_light*1.0d8/e_min_seen, ' to ',                     &
			   hp_eV*c_light*1.0d8/e_low, ' A'

			! Name the absorbers that lose part of their band: the metastable
			! helium when it is included, and every active element whose
			! NEUTRAL stage ionizes below the lowest photon the file carries.
			write(*,*) '   absorbers whose ionization threshold '//        &
			           'is below that lowest photon:'
			n_absorber = 0
			n_metal_absorber = 0
			if (use_excited_H .and. e_th_HI_n2 .lt. e_min_seen) then
				n_absorber = n_absorber + 1
				write(*,'(A,ES12.5,A,ES12.5,A)')                           &
				   '       H(n=2), Balmer continuum          : ',          &
				   e_th_HI_n2, ' eV = ',                                   &
				   hp_eV*c_light*1.0d8/e_th_HI_n2, ' A'
			endif
			if (thereis_HeITR .and. e_th_HeTR .lt. e_min_seen) then
				n_absorber = n_absorber + 1
				write(*,'(A,ES12.5,A,ES12.5,A)')                           &
				   '       He 2^3S (Include He23S? True)     : ', e_th_HeTR, &
				   ' eV = ', hp_eV*c_light*1.0d8/e_th_HeTR, ' A'
			endif
			if (allocated(melem_ab)) then
				do j = 1,n_melem
					e_thr_j = mion_ethr(melem_i0(j))
					if (melem_ab(j) .gt. 0.0d0 .and.                    &
					    e_thr_j .lt. e_min_seen) then
						n_absorber = n_absorber + 1
						n_metal_absorber = n_metal_absorber + 1
						neutral_label = trim(melem_name(j))//' I'
						write(*,'(A,A5,A,ES12.5,A,ES12.5,A)')           &
						   '       ', neutral_label,                    &
						   ' (metals.inp)                : ',           &
						   e_thr_j, ' eV = ',                           &
						   hp_eV*c_light*1.0d8/e_thr_j, ' A'
					endif
				enddo
			endif
			if (n_absorber .eq. 0) write(*,*) '      none; the floor is '//  &
			   'the E_low of the [E_low,E_mid,E_high] line of input.inp'

			write(*,*) '   remedies (there is no key to continue):'
			if (use_excited_H .and. e_th_HI_n2 .lt. e_min_seen) then
				write(*,*) '      remove "Stellar Teff [K]:" and '//    &
				           '"Stellar radius [R_sun]:" from input.inp,'
				write(*,*) '      which is what arms the '//            &
				           'excited-hydrogen coupling '//              &
				           '(input_read.f90 use_excited_H); or'
			endif
			if (thereis_HeITR .and. e_th_HeTR .lt. e_min_seen)             &
				write(*,*) '      state "Include He23S? False" '//     &
				           'in input.inp; or'
			if (n_metal_absorber .gt. 0)                                   &
				write(*,*) '      remove the element(s) named above '//   &
				           'from metals.inp; or'
			write(*,'(A,ES12.5,A)')                                        &
			   '       supply an SED that reaches ',                       &
			   hp_eV*c_light*1.0d8/e_low, ' A'
			error stop 1
		endif

		! Return to the beginning of the file
		rewind (unit = 1)

		! Allocate the table as it is read. The photon grid built from it
		! is allocated below, once the number of bins is known.
		allocate(wave_c(Nl))
		allocate(f_lam(Nl))

		! Read lines to be skipped
		do j = 1,skip
			if (.not. sed_next_row(1, dum_w, dum_f, io) .and. io .eq. 0)  &
				io = -1
			if (io .ne. 0) then
				write(*,*) '(sed_read.f90) SED file "', trim(sed_file),  &
				           '" changed or truncated during skip pass at row ', j
				error stop 'read_sed: SED file read failed'
			endif
		enddo

		! Read lines of desired energy interval, validating that wavelengths
		! are positive and strictly increasing (energy strictly decreasing),
		! the SED file convention, which also keeps the energy-bin widths
		! well defined.
		w_prev = 0.0d0
		do j = 1,Nl
			if (.not. sed_next_row(1, wave_c(j), f_lam(j), io) .and.      &
			    io .eq. 0) io = -1
			if (io .ne. 0) then
				write(*,*) '(sed_read.f90) SED file "', trim(sed_file),  &
				           '" changed or truncated during data pass at row ', j
				error stop 'read_sed: SED file read failed'
			endif
			if (wave_c(j) .le. 0.0d0) then
				write(*,*) '(sed_read.f90) SED file "', trim(sed_file),  &
				           '" has a non-positive wavelength at selected row ', j
				error stop 'read_sed: invalid wavelength in SED file'
			endif
			if (wave_c(j) .le. w_prev) then
				write(*,*) '(sed_read.f90) SED file "', trim(sed_file),  &
				           '" wavelengths are not strictly increasing at selected row ', j
				error stop 'read_sed: non-monotonic wavelengths in SED file'
			endif
			w_prev = wave_c(j)
		enddo

	! Close SED file
	close(1)
	write(*,*) ' (sed_read.f90) Done.'

	! ---- The rows of the file, ascending in photon energy --------------
	! The file is ascending in wavelength, so the row order is reversed
	! here; F_lambda is rescaled to the flux per unit photon energy F_E by
	! F_E = F_lambda hc/E^2.
	n_node = Nl
	allocate(en(n_node),fn(n_node))
	do j = 1,n_node
		k     = n_node - j + 1
		en(j) = hp_eV*c_light/(wave_c(k)*1.0d-8)
		fn(j) = f_lam(k)*hp_eV*c_light*1.0d8/en(j)**2
	enddo

	! The table AS THE FILE STATES IT, kept for the consumers that read the
	! field off the file instead of integrating over the grid: the photon
	! grid below holds bin centres, which are not table rows.
	if (allocated(e_sed_node)) deallocate(e_sed_node)
	if (allocated(F_sed_node)) deallocate(F_sed_node)
	allocate(e_sed_node(n_node),F_sed_node(n_node))
	e_sed_node = en
	F_sed_node = fn

	! ---- The bins: one per row, split at every active threshold ---------
	! Row k owns [sqrt(E_k-1 E_k), sqrt(E_k E_k+1)], the two end edges being
	! E_1 and E_N, and carries its own flux over the whole of it. A threshold
	! strictly inside a bin cuts it in two, so at most n_thr bins are added.
	allocate(b_lo(n_node + n_thr),b_hi(n_node + n_thr),                     &
	         f_bin(n_node + n_thr))
	allocate(t_in(max(n_thr,1)))
	n_bin = 0
	do k = 1,n_node
		if (k .eq. 1) then
			e_lo_k = en(1)
		else
			e_lo_k = sqrt(en(k-1)*en(k))
		endif
		if (k .eq. n_node) then
			e_hi_k = en(n_node)
		else
			e_hi_k = sqrt(en(k)*en(k+1))
		endif

		! Thresholds strictly inside this bin, deduplicated against the two
		! edges and against each other, then sorted ascending (m is at most
		! the handful of thresholds one table row can span).
		m = 0
		do i = 1,n_thr
			if (e_thr(i) .le. e_lo_k .or. e_thr(i) .ge. e_hi_k) cycle
			if (e_thr(i) - e_lo_k .le. rel_same*e_hi_k) cycle
			if (e_hi_k - e_thr(i) .le. rel_same*e_hi_k) cycle
			is_dup = .false.
			do jj = 1,m
				if (abs(e_thr(i) - t_in(jj)) .le. rel_same*e_hi_k)          &
					is_dup = .true.
			enddo
			if (is_dup) cycle
			m       = m + 1
			t_in(m) = e_thr(i)
		enddo
		do i = 1,m-1
			do jj = i+1,m
				if (t_in(jj) .lt. t_in(i)) then
					e_swap   = t_in(i)
					t_in(i)  = t_in(jj)
					t_in(jj) = e_swap
				endif
			enddo
		enddo

		! Both halves of a split carry the row's flux: the row is a histogram
		! value, so sum(F de_v) over the halves is the row's own contribution
		! whatever the split.
		e_prev = e_lo_k
		do i = 1,m
			n_bin        = n_bin + 1
			b_lo(n_bin)  = e_prev
			b_hi(n_bin)  = t_in(i)
			f_bin(n_bin) = fn(k)
			e_prev       = t_in(i)
		enddo
		n_bin        = n_bin + 1
		b_lo(n_bin)  = e_prev
		b_hi(n_bin)  = e_hi_k
		f_bin(n_bin) = fn(k)
	enddo

	! ---- The photon grid, decreasing in energy as the file is -----------
	Nl = n_bin
	allocate(e_v(Nl),de_v(Nl),F_XUV(Nl))
	do j = 1,Nl
		k        = Nl - j + 1
		e_v(j)   = sqrt(b_lo(k)*b_hi(k))
		de_v(j)  = b_hi(k) - b_lo(k)
		F_XUV(j) = f_bin(k)
	enddo

	! Number of grid bins below the H I edge. This is what NlTR means for
	! every spectrum type -- set_energy_vectors counts it while it lays out
	! the bands of the analytic types -- and write_setup_report reads it to
	! state what the run's field below 13.6 eV is built from. e_th_HI is a
	! bin edge here, so no bin straddles it and the centre test decides.
	! (Counted in a loop: global_parameters carries an integer named
	! `count`, the marching step counter, which shadows the intrinsic of
	! that name in every scope that uses the module.)
	NlTR = 0
	do j = 1,Nl
		if (e_v(j) .lt. e_th_HI) NlTR = NlTR + 1
	enddo

	! Calculate LEUV and LX luminosities
	LEUV_int  = 0.0
	LX_int    =	0.0	
	do j = 1,Nl
		
		! Skip if energy is below ionization threshold
		if (e_v(j) .lt. e_th_HI) exit
		
		! Integrate F(E) w.r.t. energy
		Df_int = de_v(j)*F_XUV(j)*4.0e0*pi*a_orb*a_orb
					
		! Update separately EUV and X-ray luminosities
		if (e_v(j) .lt. e_mid) then
			LEUV_int = LEUV_int + Df_int
			else
			LX_int = LX_int + Df_int
		endif
		
	enddo
	
	! Convert to log10 (guard: an SED with no bins in a band leaves the
	! integral at 0; log10(0) = -Inf would poison downstream fluxes, so
	! floor to a negligible luminosity instead).
	LX   = log10(max(LX_int,   1.0d-99))
	LEUV = log10(max(LEUV_int, 1.0d-99))
	
	deallocate(wave_c,f_lam,en,fn,b_lo,b_hi,f_bin,t_in)

	! End of subroutine
	end subroutine read_sed
	
	! End of module
	
      double precision function lyman_werner_band_flux_from_sed()         &
                                result(F_band)
      ! Band-integrated flux over 912-1201 A of the numerical SED, at the
      ! planet [erg cm^-2 s^-1] -- the quantity "Stellar LW flux" states.
      !
      ! THE BAND IS THE LYMAN-WERNER INTERVAL, 912-1201 A: the H Lyman edge
      ! to the start of band B2 of the FUV band list, which is also the red
      ! end of the line list the H2 self-shielding table is built from
      ! (lyman_werner.f90 sec. 1). It ran to 1110 A until 2026-09-06, when
      ! the 1110-1201 A band B1 was merged into it.
      !
      ! The integral is carried out by the code instead of by hand: the
      ! stellar spectrum over 91.2-120.1 nm at the planet's orbit. No
      ! dilution is applied here because EXHALE's own SED file is already AT
      ! THE PLANET, which is what read_sed's header states; the (R_star/a)^2
      ! step belongs to the stellar-surface files the value used to be
      ! produced from. Checked against that route on the NARROWER band, when
      ! it was the band: Gueymard's solar spectrum integrated over
      ! 912-1110 A at the stellar surface and diluted to 0.048 AU behind a
      ! 1.155 R_sun star gave 329 erg cm^-2 s^-1 against the 343 the
      ! molecular cases carried, i.e. the two agreed to 4%.
      !
      ! The interval is OUTSIDE the ionizing range read_sed retains (912 A is
      ! 13.6 eV, the H I edge, and everything longward is below it), so the
      ! file is re-read here rather than taken from the selected arrays.
      ! Trapezoid on the bin centres; rows outside the band are skipped, and
      ! the two rows bracketing each edge are kept so a coarse grid does not
      ! lose the ends.
      real*8, parameter :: w_lo = 912.0d0, w_hi = 1201.0d0
      real*8  :: w, f, w_prev, f_prev, wa, wb
      integer :: io, nin
      logical :: have_prev
      F_band    = 0.0d0
      nin       = 0
      have_prev = .false.
      w_prev    = 0.0d0
      f_prev    = 0.0d0
      open(unit = 71, file = sed_file, status = 'old', iostat = io)
      if (io .ne. 0) return
      do
         if (.not. sed_next_row(71, w, f, io)) exit
         if (io .ne. 0) exit
         if (have_prev .and. w .gt. w_prev) then
            wa = max(w_prev, w_lo)
            wb = min(w,      w_hi)
            if (wb .gt. wa) then
               ! trapezoid of the segment, clipped to the band
               F_band = F_band + 0.5d0*(f_prev + f)*(wb - wa)
               nin    = nin + 1
            endif
         endif
         w_prev    = w
         f_prev    = f
         have_prev = .true.
         if (w .gt. w_hi) exit
      enddo
      close(71)
      if (nin .lt. 2) F_band = 0.0d0
      end function lyman_werner_band_flux_from_sed

      ! ------------------------------------------------------------- !

      logical function sed_next_row(unit, w, f, io)
      ! Next data row of an SED file: blank lines and lines whose first
      ! non-blank character is '#' are skipped, so a file may carry a
      ! provenance header. Returns .false. at end of file; io > 0 is a
      ! malformed data row and is left for the caller to report.
      integer, intent(in)  :: unit
      real*8,  intent(out) :: w, f
      integer, intent(out) :: io
      character(len=512)   :: ln
      sed_next_row = .false.
      io = 0
      do
         read(unit,'(A)',iostat = io) ln
         if (io .lt. 0) return
         if (io .gt. 0) return
         ln = adjustl(ln)
         if (len_trim(ln) .eq. 0) cycle
         if (ln(1:1) .eq. '#') cycle
         read(ln,*,iostat = io) w, f
         sed_next_row = (io .eq. 0)
         return
      enddo
      end function sed_next_row

	end module sed_reader	

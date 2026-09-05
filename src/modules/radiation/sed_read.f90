	module sed_reader
	! Module containing the subroutine for reading the numerical SED
		
	use global_parameters
	use species_table, only: n_melem, mion_ethr, melem_i0

	implicit none

	contains

	!---------------------------------------------------!
	
	subroutine read_sed
	! Subroutine to read the numerical SED from file.
	!	The file has to be formatted in two columns:
	!	   1) bin central wavelength [Angstrom] 
	!	   2) flux at the planet surface [erg cm^{-2} s^{-2} A^{-1}]
	!	Data have to be order with increasing wavelength
	
	integer :: io
	integer :: j,skip
	integer :: nrow          ! running count of rows read from the SED file
	real*8 :: e_top_read
	real*8 :: dum_w,dum_f,dum_e
	real*8 :: w_prev         ! previous selected wavelength (monotonicity check)
	real*8 :: e_min_seen     ! lowest selected photon energy encountered
	real*8 :: LEUV_int,LX_int,Df_int
	logical :: hit_eof       ! true if the selection block ran to end of file
	real*8, dimension(:), allocatable :: wave_c

   write(*,*) '(sed_read.f90) Reading the numerical spectrum..'
   
	! Open file to read number of lines to be skipped
	!	according to the selected energy interval
	open(unit = 1, file = sed_file, status = 'old')
	
		! Extend to 4.8 eV if He triplet is included
		if (thereis_HeITR) e_low = e_th_HeTR

		! Lower e_low to the smallest active neutral-metal ionization
		! threshold below the HI edge, so photons under 13.6 eV that the
		! low-IP metals (Mg, Si, Ca, Na, K, Fe) photoionize are retained.
		if (thereis_lowIP_metal) then
			do j = 1,n_melem
				if (melem_ab(j) .gt. 0.0d0 .and.                  &
				    mion_ethr(melem_i0(j)) .lt. e_th_HI)          &
					e_low = min(e_low, mion_ethr(melem_i0(j)))
			enddo
		endif

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
			dum_e = hp_eV*c_light/(dum_w*1e-8)

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

		! Warn if the file ran out before reaching e_low: the lowest
		! available photon energy is above e_low, so photoionization by
		! low-IP metals in the missing band may be underestimated.
		if (hit_eof .and. e_min_seen .gt. e_low) then
			write(*,*) '(sed_read.f90) WARNING: SED file "', trim(sed_file), &
			           '" ends above e_low; missing band [', e_low, ',',     &
			           e_min_seen, '] eV.'
		endif

		! Return to the beginning of the file
		rewind (unit = 1)

		! Allocate energy vectors
		allocate(wave_c(Nl))
		allocate(F_XUV(Nl))
		allocate(e_v(Nl),de_v(Nl))

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
			if (.not. sed_next_row(1, wave_c(j), F_XUV(j), io) .and.      &
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

	! Convert wavelength to energy
	e_v   = hp_eV*c_light/(wave_c*1e-8)
	
	! Calculate energy bins width 
	de_v(1)      = -0.5*(e_v(2) - e_v(1))
	de_v(2:Nl-1) = -0.5*(e_v(3:Nl) - e_v(1:Nl-2))
	de_v(Nl)     = -0.5*(e_v(Nl) - e_v(Nl-1))

	! Rescale from F_lambda to F_E
	F_XUV = F_XUV*hp_eV*c_light*1.0e8/e_v**2.0e0
	
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
	
	deallocate(wave_c)

	! End of subroutine
	end subroutine read_sed
	
	! End of module
	
      double precision function lyman_werner_band_flux_from_sed()         &
                                result(F_band)
      ! Band-integrated flux over 912-1110 A of the numerical SED, at the
      ! planet [erg cm^-2 s^-1] -- the quantity "Stellar LW flux" states.
      !
      ! THIS IS THE PRESCRIPTION THE BAND FLUX HAS ALWAYS BEEN BUILT BY, now
      ! carried out by the code instead of by hand: integrate the stellar
      ! spectrum over 91.2-111.0 nm and take it at the planet's orbit
      ! (docs/Update_EXHALE.md section 109). No dilution is applied here
      ! because EXHALE's own SED file is already AT THE PLANET, which is what
      ! read_sed's header states; the (R_star/a)^2 step belongs to the
      ! stellar-surface files the value used to be produced from. Checked
      ! against that route: Gueymard's solar spectrum integrated over the band
      ! at the stellar surface and diluted to 0.048 AU behind a 1.155 R_sun
      ! star gives 329 erg cm^-2 s^-1 against the 343 the molecular cases
      ! carry, i.e. the two agree to 4%.
      !
      ! The interval is OUTSIDE the ionizing range read_sed retains (912 A is
      ! 13.6 eV, the H I edge, and everything longward is below it), so the
      ! file is re-read here rather than taken from the selected arrays.
      ! Trapezoid on the bin centres; rows outside the band are skipped, and
      ! the two rows bracketing each edge are kept so a coarse grid does not
      ! lose the ends.
      real*8, parameter :: w_lo = 912.0d0, w_hi = 1110.0d0
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

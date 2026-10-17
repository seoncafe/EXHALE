      module wae_spectrum
      ! Spectrum loading for the Wind-AE port.
      ! Reads inputs/spectrum.inp (GLQ nodes, cross sections, ion. pots),
      ! exactly the fields wind-ae's init_glq() consumes for the rate
      ! calculation, with the C parsing. The mass fractions, atomic masses
      ! and Ftot below are copied from the parameter list by
      ! wae_continuation.
      use wae_config, only: wae_nspecies
      implicit none

      integer :: wae_npts = 0
      real*8, allocatable :: wae_hc_over_wl(:)     ! hc/lambda_i  [erg]
      real*8, allocatable :: wae_wPhi_wl(:)        ! w_i Phi_i / Ftot
      real*8, allocatable :: wae_sigma_wl(:,:)     ! (npts, nspecies) [cm^2]
      real*8 :: wae_ion_pot(wae_nspecies) = 0.0d0  ! [erg]

      ! parameters (PARAMLIST subset needed by the rate calc)
      real*8 :: wae_HX(wae_nspecies)          = 0.0d0   ! mass fractions
      real*8 :: wae_atomic_mass(wae_nspecies) = 0.0d0   ! [g]
      real*8 :: wae_Ftot = 0.0d0                         ! [erg/cm^2/s]

      contains

      subroutine wae_load_spectrum(fname)
      ! Mirror of init_glq(): parse the '# KEY: ...' header for NPTS and
      ! IONPOTS, then read npts data rows "hc, wPhi, sigma_HI, sigma_HeI".
      character(len=*), intent(in) :: fname
      character(len=4096) :: line
      integer :: u, ios, i, j, p, c
      character(len=64) :: tok

      open(newunit=u, file=fname, status='old', action='read')
      ! ---- header ----
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) stop '(wae_load_spectrum) unexpected EOF in header'
         if (line(1:1) .ne. '#') then
            backspace(u); exit
         end if
         if (index(line,'NPTS') .gt. 0) then
            call after_colon(line, tok); read(tok,*) wae_npts
         else if (index(line,'IONPOTS') .gt. 0) then
            ! comma-separated list after the colon
            p = index(line,':') + 1
            do j = 1, wae_nspecies
               call next_csv(line, p, tok)
               read(tok,*) wae_ion_pot(j)
            end do
         end if
      end do
      if (wae_npts .le. 0) stop '(wae_load_spectrum) NPTS not found'

      allocate(wae_hc_over_wl(wae_npts), wae_wPhi_wl(wae_npts), &
               wae_sigma_wl(wae_npts, wae_nspecies))

      ! ---- data rows ----
      i = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '#') cycle
         i = i + 1
         if (i .gt. wae_npts) exit
         p = 1
         call next_csv(line, p, tok); read(tok,*) wae_hc_over_wl(i)
         call next_csv(line, p, tok); read(tok,*) wae_wPhi_wl(i)
         do c = 1, wae_nspecies
            call next_csv(line, p, tok); read(tok,*) wae_sigma_wl(i,c)
         end do
      end do
      close(u)
      if (i .lt. wae_npts) stop '(wae_load_spectrum) too few data rows'
      end subroutine wae_load_spectrum

      ! ---- small CSV/token helpers ----

      subroutine after_colon(line, tok)
      ! token after ':' up to the first ',' (or end)
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: tok
      integer :: a, b
      a = index(line,':') + 1
      b = index(line(a:), ',')
      if (b .eq. 0) then
         tok = adjustl(line(a:))
      else
         tok = adjustl(line(a:a+b-2))
      end if
      end subroutine after_colon

      subroutine next_csv(line, p, tok)
      ! Read the next comma-separated token of `line` starting at position
      ! p; advance p past the consumed comma. Trailing token (no comma)
      ! is returned and p set past the end.
      character(len=*), intent(in)    :: line
      integer,          intent(inout) :: p
      character(len=*), intent(out)   :: tok
      integer :: b, n
      n = len_trim(line)
      if (p .gt. n) then
         tok = ''
         return
      end if
      b = index(line(p:n), ',')
      if (b .eq. 0) then
         tok = adjustl(line(p:n)); p = n + 1
      else
         tok = adjustl(line(p:p+b-2)); p = p + b
      end if
      end subroutine next_csv

      end module wae_spectrum

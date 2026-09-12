module isolated
implicit none
integer, parameter :: dp=kind(1d0), N=1, Ng=1, n_species=1
integer :: eos_calls=0, scenario=1
character(len=256) :: sed_file
type ioniz_eq_ledger
integer :: n_nonfinite=0, n_offsimplex=0
end type
contains
double precision function sed_band_integrated_flux(w_lo, w_hi)      &
                                result(F_band)
      ! Band-integrated flux of the numerical SED over [w_lo, w_hi] A, at
      ! the planet [erg cm^-2 s^-1]: the quantity the band-flux keys state
      ! ("Stellar LW flux" over 912-1201 A, "Stellar FUV B3 flux" over
      ! 1231-1450 A, "Stellar FUV B4 flux" over 1451-2304 A; the edges are
      ! fuv_band_lo_A / fuv_band_hi_A of oxygen_rates). Band B2, the
      ! Ly-alpha line, is NOT integrated here: a line flux reconstructed
      ! from observations is a better number than a trapezoid over
      ! whatever rows the file has across 1202-1230 A, so it stays a key
      ! ("Stellar Lya flux").
      !
      ! The Lyman-Werner interval is 912-1201 A: the H Lyman edge to the
      ! start of band B2, which is also the red end of the line list the
      ! H2 self-shielding table is built from (lyman_werner.f90 sec. 1).
      ! It ran to 1110 A until 2026-09-06, when the 1110-1201 A band B1 was
      ! merged into it.
      !
      ! The integral is carried out by the code instead of by hand. No
      ! dilution is applied here because EXHALE's own SED file is already AT
      ! THE PLANET, which is what read_sed's header states; the (R_star/a)^2
      ! step belongs to the stellar-surface files the value used to be
      ! produced from. Checked against that route on the NARROWER
      ! Lyman-Werner band, when it was the band: Gueymard's solar spectrum
      ! integrated over 912-1110 A at the stellar surface and diluted to
      ! 0.048 AU behind a 1.155 R_sun star gave 329 erg cm^-2 s^-1 against
      ! the 343 the molecular cases carried, i.e. the two agreed to 4%.
      !
      ! Every one of these intervals is OUTSIDE the ionizing range read_sed
      ! retains (912 A is 13.6 eV, the H I edge, and everything longward is
      ! below it), so the file is re-read here rather than taken from the
      ! selected arrays. Trapezoid on the bin centres; rows outside the band
      ! are skipped, and the two rows bracketing each edge are kept so a
      ! coarse grid does not lose the ends. Zero if fewer than two segments
      ! fall in the band (a file that does not reach it).
      real*8, intent(in) :: w_lo, w_hi
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
      end function sed_band_integrated_flux
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
subroutine equilibrate_chemistry_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, heat, cool, eta, ok,   &
                                             n_cycles)
      real(dp), dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: p, T
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: heat, cool
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: eta
      logical,                                  intent(out)   :: ok
      integer,                                  intent(out)   :: n_cycles

      integer,  parameter :: cycles_max = 5
      real(dp), parameter :: t_cycle_tol = 1.0d-6
      real(dp), dimension(1-Ng:N+Ng) :: ntot, ne, T_prev
      type(ioniz_eq_ledger) :: ledger
      integer :: k

      ok = .true.
      n_cycles = 0
      call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,     &
                                                             p, T, ntot, ne)
      do k = 1, cycles_max
         T_prev = T
         call ioniz_eq(T, u(1,:), f_sp, heat, cool, eta, ledger)
         n_cycles = k
         if (ledger%n_nonfinite .gt. 0 .or.                               &
             .not. all(f_sp(1:N,:) .eq. f_sp(1:N,:)) .or.                 &
             .not. all(abs(f_sp(1:N,:)) .le. huge(1.0d0))) then
            ok = .false.
            return
         endif
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                                             p, T, ntot, ne)
         if (maxval(abs(T(1:N) - T_prev(1:N))/max(T(1:N), 1.0d-300))      &
             .lt. t_cycle_tol) return
      enddo
      end subroutine equilibrate_chemistry_at_fixed_conserved_state
subroutine pressure_and_temperature_at_fixed_conserved_state(u,f,p,t,ntot,ne)
real(dp), intent(in) :: u(3,0:2),f(0:2,1)
real(dp), intent(out) :: p(0:2),t(0:2),ntot(0:2),ne(0:2)
eos_calls=eos_calls+1
p=1d0; t=100d0; ntot=1d0; ne=0d0
if (scenario==1) t=100d0+dble(eos_calls)
end subroutine
subroutine ioniz_eq(t,rho,f,heat,cool,eta,ledger)
real(dp), intent(in) :: t(0:2),rho(0:2)
real(dp), intent(inout) :: f(0:2,1),heat(0:2),cool(0:2),eta(0:2)
type(ioniz_eq_ledger), intent(out) :: ledger
ledger=ioniz_eq_ledger()
if (scenario==2) ledger%n_offsimplex=1
end subroutine
end module
program probe
use isolated
implicit none
real(dp) :: u(3,0:2), f(0:2,1), p(0:2),t(0:2),h(0:2),c(0:2),eta(0:2)
integer :: cycles
logical :: ok
sed_file='linear_sed.txt'
write(*,'(A,ES24.16)') 'clipped_linear_flux=',sed_band_integrated_flux(950d0,1100d0)
sed_file='single_segment_sed.txt'
write(*,'(A,ES24.16)') 'single_segment_flux=',sed_band_integrated_flux(950d0,1100d0)
u=1d0; f=1d0; h=0d0; c=0d0; eta=0d0
call equilibrate_chemistry_at_fixed_conserved_state(u,f,p,t,h,c,eta,ok,cycles)
write(*,'(A,L1,A,I0,A,F10.6)') 'cycle_exhaustion_ok=',ok,', cycles=',cycles, &
 ', last_relative_increment=',1d0/t(1)
scenario=2; eos_calls=0
call equilibrate_chemistry_at_fixed_conserved_state(u,f,p,t,h,c,eta,ok,cycles)
write(*,'(A,L1,A,I0)') 'off_simplex_ok=',ok,', cycles=',cycles
end program

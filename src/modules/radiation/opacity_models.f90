      module opacity_models
      ! AIOLOS-style opacity-model dispatcher for photoionization
      ! cross sections. Selects between
      !
      !   'A' analytic (cross_sec.f90; current ATES default)
      !   'C' constant   (sigma_threshold * user factor, above threshold)
      !   'P' physical   (constant + Robinson&Catling pressure broadening)
      !   'T' tabulated  (.opa two-column file: E[eV], sigma[1e-18 cm^2])
      !
      ! Cross sections are returned in ATES's internal unit of 1e-18 cm^2,
      ! same as cross_sec.f90, so downstream code is unchanged.
      ! A separate runtime multiplier opacity_pT_factor(p) is applied
      ! per cell in 'P' mode (see PH_heat_HHe / calc_column_dens).
      !
      ! Ported from ATES_extended; in ATES-Code-main the per-cell
      ! pressure factor is actually applied (it was a deferred hook in
      ! ATES_extended).

      use global_parameters
      use Cross_sections

      implicit none

      private
      public :: photoion_sigma
      public :: opacity_pT_factor
      public :: load_opacity_tables
      public :: free_opacity_tables

      type :: opa_table
         integer                            :: n = 0
         real*8, dimension(:), allocatable :: e
         real*8, dimension(:), allocatable :: sig
      end type opa_table

      type(opa_table), save, target :: tab_HI
      type(opa_table), save, target :: tab_HeI
      type(opa_table), save, target :: tab_HeII
      type(opa_table), save, target :: tab_HeITR

      contains

      !----------------------------------------------------------!

      subroutine load_opacity_tables
      ! Load every per-species .opa file whose path was set in
      ! opacity.inp. Empty paths are no-ops (the dispatcher will fall
      ! back to the analytic curve for that species).

      if (allocated(opa_file_HI)) then
         if (len_trim(opa_file_HI) > 0) call load_one(opa_file_HI, tab_HI)
      endif
      if (allocated(opa_file_HeI)) then
         if (len_trim(opa_file_HeI) > 0) call load_one(opa_file_HeI, tab_HeI)
      endif
      if (allocated(opa_file_HeII)) then
         if (len_trim(opa_file_HeII) > 0) call load_one(opa_file_HeII, tab_HeII)
      endif
      if (allocated(opa_file_HeITR)) then
         if (len_trim(opa_file_HeITR) > 0) call load_one(opa_file_HeITR, tab_HeITR)
      endif

      end subroutine load_opacity_tables

      !----------------------------------------------------------!

      subroutine load_one(path, tab)
      character(len=*), intent(in)    :: path
      type(opa_table),  intent(inout) :: tab
      integer :: io, n, k
      real*8  :: e_dum, s_dum

      open(unit = 37, file = path, status = 'old', action = 'read')

      n = 0
      do
         read(37, *, iostat = io) e_dum, s_dum
         if (io /= 0) exit
         n = n + 1
      enddo

      if (n < 2) then
         write(*,*) '(opacity_models) ERROR: ', trim(path),                 &
                    ' has fewer than 2 rows.'
         close(37)
         stop 1
      endif

      rewind(37)
      allocate(tab%e(n), tab%sig(n))
      do k = 1, n
         read(37, *) tab%e(k), tab%sig(k)
      enddo
      close(37)
      tab%n = n

      ! Sanity: enforce monotonically increasing energy
      do k = 2, n
         if (tab%e(k) <= tab%e(k-1)) then
            write(*,*) '(opacity_models) ERROR: ', trim(path),              &
                       ' energy column not strictly increasing at row', k
            stop 1
         endif
      enddo

      write(*,'(a,i0,a,a)') ' (opacity_models) Loaded ', n,                 &
                            ' rows from ', trim(path)

      end subroutine load_one

      !----------------------------------------------------------!

      subroutine free_opacity_tables
      if (allocated(tab_HI%e))    deallocate(tab_HI%e,    tab_HI%sig)
      if (allocated(tab_HeI%e))   deallocate(tab_HeI%e,   tab_HeI%sig)
      if (allocated(tab_HeII%e))  deallocate(tab_HeII%e,  tab_HeII%sig)
      if (allocated(tab_HeITR%e)) deallocate(tab_HeITR%e, tab_HeITR%sig)
      tab_HI%n    = 0
      tab_HeI%n   = 0
      tab_HeII%n  = 0
      tab_HeITR%n = 0
      end subroutine free_opacity_tables

      !----------------------------------------------------------!

      double precision function photoion_sigma(sp, E)
      ! Returns sigma(E) for species sp ('HI','HeI','HeII','HeITR')
      ! in 1e-18 cm^2 units, dispatched by global opacity_model.
      character(len=*), intent(in) :: sp
      real*8,           intent(in) :: E

      select case (opacity_model)
         case ('C', 'P')
            photoion_sigma = constant_sigma(sp, E)
         case ('T')
            photoion_sigma = tabulated_sigma(sp, E)
         case default
            photoion_sigma = analytic_sigma(sp, E)
      end select
      end function photoion_sigma

      !----------------------------------------------------------!

      double precision function analytic_sigma(sp, E)
      character(len=*), intent(in) :: sp
      real*8,           intent(in) :: E
      select case (sp)
         case ('HI');    analytic_sigma = sigma(E, ih)
         case ('HeI');   analytic_sigma = sigma_HeI(E)
         case ('HeII');  analytic_sigma = sigma(E, ihe)
         case ('HeITR'); analytic_sigma = sigma_HeI23S(E)
         case default;   analytic_sigma = 0.0d0
      end select
      end function analytic_sigma

      !----------------------------------------------------------!

      double precision function constant_sigma(sp, E)
      ! Constant cross section above threshold:
      !     sigma = factor * analytic_sigma(threshold + eps)
      ! Pegging to the threshold value keeps the model physically
      ! anchored regardless of which species/units convention is in use.
      character(len=*), intent(in) :: sp
      real*8,           intent(in) :: E
      real*8 :: thr, fac

      select case (sp)
         case ('HI');    thr = e_th_HI;   fac = opa_const_HI
         case ('HeI');   thr = e_th_HeI;  fac = opa_const_HeI
         case ('HeII');  thr = e_th_HeII; fac = opa_const_HeII
         case ('HeITR'); thr = e_th_HeTR; fac = opa_const_HeITR
         case default
            constant_sigma = 0.0d0
            return
      end select

      if (E < 0.99999d0*thr) then
         constant_sigma = 0.0d0
      else
         constant_sigma = fac * analytic_sigma(sp, thr*1.00001d0)
      endif
      end function constant_sigma

      !----------------------------------------------------------!

      double precision function tabulated_sigma(sp, E)
      character(len=*), intent(in) :: sp
      real*8,           intent(in) :: E
      type(opa_table), pointer :: tab

      select case (sp)
         case ('HI');    tab => tab_HI
         case ('HeI');   tab => tab_HeI
         case ('HeII');  tab => tab_HeII
         case ('HeITR'); tab => tab_HeITR
         case default
            tabulated_sigma = 0.0d0
            return
      end select

      if (tab%n == 0) then
         ! No table provided for this species: fall back to analytic.
         tabulated_sigma = analytic_sigma(sp, E)
      else
         tabulated_sigma = loglog_interp(tab%e, tab%sig, tab%n, E)
      endif
      end function tabulated_sigma

      !----------------------------------------------------------!

      double precision function loglog_interp(x, y, n, x0)
      ! Bisection + log-log linear interpolation. Falls back to linear
      ! interpolation in y when either endpoint is non-positive.
      ! Out-of-range queries return the nearest endpoint (clamp).
      integer,              intent(in) :: n
      real*8, dimension(n), intent(in) :: x, y
      real*8,               intent(in) :: x0
      integer :: lo, hi, mid
      real*8  :: t

      if (x0 <= x(1)) then
         loglog_interp = y(1)
         return
      endif
      if (x0 >= x(n)) then
         loglog_interp = y(n)
         return
      endif

      lo = 1; hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (x(mid) > x0) then
            hi = mid
         else
            lo = mid
         endif
      enddo

      if (y(lo) <= 0.0d0 .or. y(hi) <= 0.0d0) then
         t = (x0 - x(lo)) / (x(hi) - x(lo))
         loglog_interp = y(lo) + t*(y(hi) - y(lo))
      else
         t = (log(x0) - log(x(lo))) / (log(x(hi)) - log(x(lo)))
         loglog_interp = exp(log(y(lo)) + t*(log(y(hi)) - log(y(lo))))
      endif
      end function loglog_interp

      !----------------------------------------------------------!

      double precision function opacity_pT_factor(p_dyne)
      ! Pressure-broadening multiplier applied to the base cross
      ! section in 'P' mode only (1.0 otherwise). Robinson & Catling
      ! (2012) form:
      !    f(p) = 1 + a * (p / p_pivot)^n
      ! with p in dyne/cm^2; default pivot 1e5 dyne/cm^2 = 0.1 bar.
      real*8, intent(in) :: p_dyne
      if (opacity_model == 'P') then
         opacity_pT_factor = 1.0d0                                          &
              + opa_pb_factor*(max(p_dyne,0.0d0)/opa_pb_pivot)**opa_pb_exponent
      else
         opacity_pT_factor = 1.0d0
      endif
      end function opacity_pT_factor

      ! End of module
      end module opacity_models

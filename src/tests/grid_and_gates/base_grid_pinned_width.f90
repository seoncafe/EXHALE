      program base_grid_pinned_width
      ! QUANTITY UNDER TEST
      !   The cell centers r(j), j = 1-Ng..N+Ng, of the Mixed grid built by
      !   the PRODUCTION define_grid (module grid_construction), for four
      !   widths of the uniform base cells, with N_low_cells = 50 and the
      !   500-cell grid of the LHS 1140 b catalog (Domain mode: Spherical,
      !   Outer radius 30 R_p) and of a Roche-bounded 10 R_p domain.
      !
      !     default  dr_base_default, the width of an input without the key
      !     key_2e-4 the list-directed read of "2.0e-4", the width of
      !              "Base grid [dr,cells]: 2.0e-4 50"
      !     old      the default-real literal 2.0e-4 widened to double, the
      !              width every input without the key was run on before
      !              2026-09-19
      !     pinned   the list-directed read of "1.9999999494757503e-4", the
      !              width src/utils/pin_base_grid.py wrote into those inputs
      !
      ! ASSERTIONS (every comparison is on the bits of each coordinate)
      !   1. default and key_2e-4 build the same grid: the key spelled the
      !      way a reader writes the default states the default.
      !   2. pinned and old build the same grid: the pinned line reproduces
      !      the grid the stored results of a pinned input were written on.
      !   3. default and old build different grids. The largest relative
      !      displacement of a center is REPORTED for each domain, as an
      !      observation of that grid and not against a reference: it is
      !      what makes load_IC (tolerance 1e-10) refuse a state of one grid
      !      on the other.
      use global_parameters
      use grid_construction, only: define_grid

      implicit none

      integer, parameter :: nd = 2
      real*8  :: rmax_d(nd) = (/ 30.0d0, 10.0d0 /)
      real*8, allocatable :: r_default(:), r_key(:), r_old(:), r_pinned(:)
      real*8  :: w_key, w_pinned, w_old, dmax
      real    :: w_single
      integer :: k, nfail, j, jw
      character(len=64) :: tag
      character(len=32) :: spelled

      nfail = 0
      N = 500
      call allocate_grid_arrays
      N_low_cells = 50
      r_esc  = 2.0d0
      r_flux = 1.20d0
      grid_type = 'Mixed'
      allocate(r_default(1-Ng:N+Ng), r_key(1-Ng:N+Ng),                       &
               r_old(1-Ng:N+Ng), r_pinned(1-Ng:N+Ng))

      spelled = '2.0e-4'
      read(spelled, *) w_key
      spelled = '1.9999999494757503e-4'
      read(spelled, *) w_pinned
      w_single = 2.0e-4
      w_old    = w_single

      do k = 1, nd
         r_max = rmax_d(k)
         dr_base = dr_base_default;  call define_grid;  r_default = r
         dr_base = w_key;            call define_grid;  r_key     = r
         dr_base = w_old;            call define_grid;  r_old     = r
         dr_base = w_pinned;         call define_grid;  r_pinned  = r

         write(tag,'(A,F0.1,A)') '[r_max=', rmax_d(k), ']'
         call same_bits('default_grid_equals_key_2.0e-4_grid'//trim(tag),   &
                        r_default, r_key, nfail)
         call same_bits('pinned_grid_equals_pre_20260919_default_grid'//    &
                        trim(tag), r_pinned, r_old, nfail)

         dmax = 0.0d0
         jw   = 1
         do j = 1, N
            if (abs(r_default(j) - r_old(j))/r_old(j) .gt. dmax) then
               dmax = abs(r_default(j) - r_old(j))/r_old(j)
               jw   = j
            endif
         enddo
         if (dmax .gt. 0.0d0) then
            write(*,'(A,A,A)') 'PASS default_grid_differs_from_pre_20260919', &
                 trim(tag), ' measured=differs reference=differs tol=0'
         else
            write(*,'(A,A,A)') 'FAIL default_grid_differs_from_pre_20260919', &
                 trim(tag), ' measured=identical reference=differs tol=0'
            nfail = nfail + 1
         endif
         write(*,'(A,ES12.5,A,I0,A,I0)') '     max relative displacement ',   &
              dmax, ' at physical cell ', jw, ' of ', N
      enddo

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'base_grid_pinned_width: ', nfail,              &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'base_grid_pinned_width: all assertions PASSED'

      contains

      subroutine same_bits(name, a, b, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: a(1-Ng:), b(1-Ng:)
      integer,          intent(inout) :: nfail
      integer   :: j, ndiff
      integer*8 :: ia, ib
      ndiff = 0
      do j = 1-Ng, N+Ng
         ia = transfer(a(j), ia)
         ib = transfer(b(j), ib)
         if (ia .ne. ib) ndiff = ndiff + 1
      enddo
      if (ndiff .eq. 0) then
         write(*,'(A,A,A)') 'PASS ', name,                                   &
              ' measured=0_cells_differ reference=0 tol=0'
      else
         write(*,'(A,A,A,I0,A)') 'FAIL ', name, ' measured=', ndiff,         &
              '_cells_differ reference=0 tol=0'
         nfail = nfail + 1
      endif
      end subroutine same_bits

      end program base_grid_pinned_width

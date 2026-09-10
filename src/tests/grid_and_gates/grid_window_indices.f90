      program grid_window_indices
      ! WHAT IS UNDER TEST: the two index statements of the radial grid that
      ! decide WHICH CELLS a quantity is read from.
      !
      !   (1) cell_nearest_radius(x) is a SUBSCRIPT of r, i.e.
      !         |r(cell_nearest_radius(x)) - x| = min_j |r(j) - x| ,
      !       the minimum taken over the whole declared range 1-Ng .. N+Ng.
      !       minloc counts positions from 1 whatever the declared lower
      !       bound is, so on the grid arrays (declared 1-Ng:N+Ng) a bare
      !       minloc result is Ng cells below the cell it names. The initial
      !       condition reads the mid-domain density through this function
      !       (set_IC.f90), so the cell it names has to be the cell meant.
      !
      !   (2) the constant-momentum window [j_min:N] and the flux window
      !       [j_flux:N] are ranges of PHYSICAL cells: j_min >= 1 and
      !       j_flux >= 1 whatever the escape radius is. The ghosts
      !       (1-Ng .. 0) carry the boundary closure and not a solution, so a
      !       spread taken over them is not a property of the wind. An
      !       "Escape radius [R_p]:" at or below the base radius puts the
      !       ghosts inside the window unless the search result is clamped.
      !
      ! Both are exact statements about integers, so both are tested exactly
      ! (tolerance 0). The grid is built by the PRODUCTION routine
      ! define_grid on the PRODUCTION globals; nothing here rebuilds a grid.
      !
      ! Assertion (3) states that the clamp of (2) does not move an ordinary
      ! configuration: at the escape radius every regression case uses, the
      ! window starts where the unclamped search puts it.

      use global_parameters
      use grid_construction, only: define_grid, cell_nearest_radius

      implicit none

      integer :: nfail
      nfail = 0

      ! Same resolution and base grid as grid_width_identity, i.e. the
      ! input_read defaults.
      N           = 500
      call allocate_grid_arrays
      dr_base     = 2.0e-4
      N_low_cells = 50
      grid_type   = 'Mixed'
      r_max       = 10.0d0
      r_esc       = 1.50d0
      r_flux      = 1.20d0

      call define_grid

      ! (1) the nearest-cell index, at three radii: mid-domain (the one
      !     set_IC asks for), just above the base and just below the top.
      call check_nearest(0.5d0*(r_max + 1.0d0), 'mid_domain', nfail)
      call check_nearest(1.0d0 + 0.01d0*(r_max - 1.0d0), 'near_base', nfail)
      call check_nearest(r_max - 0.01d0*(r_max - 1.0d0), 'near_top',  nfail)

      ! (3) an ordinary escape radius leaves both windows where the radius
      !     search puts them: the clamp may not move a window that already
      !     starts inside the physical cells. The reference is an
      !     independent scan for the first cell with r >= the radius.
      call verdict_int('window_start_matches_search[r_esc=1.50]', j_min,  &
                       first_cell_above(r_esc), nfail)
      call verdict_int('window_start_matches_search[r_flux=1.20]',        &
                       j_flux, first_cell_above(r_flux), nfail)
      write(*,'(A,I5,A,I5)') '     j_min=', j_min, '  j_flux=', j_flux

      ! (2) the escape radius at the base radius itself. r(2-Ng) = 1 on the
      !     Mixed grid, so the unclamped search returns 2-Ng = 0 and the
      !     window would open in the ghosts.
      r_esc  = 1.0d0
      r_flux = 1.0d0
      call define_grid
      call verdict_ge('window_start_is_physical[r_esc=1.00]',  j_min, 1,  &
                      nfail)
      call verdict_ge('window_start_is_physical[r_flux=1.00]', j_flux, 1, &
                      nfail)
      write(*,'(A,I5,A,I5,A,ES12.5)') '     j_min=', j_min,               &
           '  j_flux=', j_flux, '  r(1-Ng)=', r(1-Ng)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'grid_window_indices: ', nfail,              &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'grid_window_indices: all assertions PASSED'

      contains

      subroutine check_nearest(x, tag, nfail)
      real*8,           intent(in)    :: x
      character(len=*), intent(in)    :: tag
      integer,          intent(inout) :: nfail
      integer :: j, jnear, jbrute
      real*8  :: dbest, d
      jnear  = cell_nearest_radius(x)
      ! Independent brute-force scan over the declared range.
      jbrute = 1 - Ng
      dbest  = abs(r(1-Ng) - x)
      do j = 1-Ng, N+Ng
         d = abs(r(j) - x)
         if (d .lt. dbest) then
            dbest  = d
            jbrute = j
         endif
      enddo
      call verdict_int('cell_nearest_radius['//tag//']', jnear, jbrute,   &
                       nfail)
      write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                            &
           '     target=', x, '  r(returned)=', r(jnear),                &
           '  r(nearest)=', r(jbrute)
      end subroutine check_nearest

      ! First cell whose center lies at or above a radius, scanned over the
      ! whole declared range: the unclamped answer the window search gives.
      integer function first_cell_above(x) result(jf)
      real*8, intent(in) :: x
      integer :: j
      jf = N + Ng + 1
      do j = 1-Ng, N+Ng
         if (r(j) .ge. x) then
            jf = j
            return
         endif
      enddo
      end function first_cell_above

      subroutine verdict_ge(name, measured, reference, nfail)
      ! One-sided verdict: the measured index must not lie below the
      ! reference (the first physical cell).
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .ge. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name,                     &
              ' measured=', measured, ' reference=>=', reference,        &
              ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name,                     &
              ' measured=', measured, ' reference=>=', reference,        &
              ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_ge

      subroutine verdict_int(name, measured, reference, nfail)
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name,                     &
              ' measured=', measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name,                     &
              ' measured=', measured, ' reference=', reference, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_int

      end program grid_window_indices

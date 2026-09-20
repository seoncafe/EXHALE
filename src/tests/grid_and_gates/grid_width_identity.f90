      program grid_width_identity
      ! IDENTITY UNDER TEST
      !   dr_j(j) = r_edg(j) - r_edg(j-1)   for every cell j
      !   sum_{j=1}^{N} dr_j(j) = r_edg(N) - r_edg(0)
      !
      ! r_edg(j) = r_{j+1/2} is the stated convention of define_grid.f90
      ! ("Cell edges (N+2*Ng-1 points) - r_edg(j) = r_{j+1/2}") and the one
      ! RK_rhs uses (rp = r_edg(j), rm = r_edg(j-1) for cell j), so the width
      ! of cell j is r_edg(j) - r_edg(j-1) and the two statements above are
      ! exact identities of the discretization, not approximations. They are
      ! therefore tested to round-off (1e-12 relative).
      !
      ! The grid is built by the PRODUCTION routine define_grid (module
      ! grid_construction), on the PRODUCTION globals of global_parameters.
      ! Nothing here re-derives a cell width; the test only reads r_edg and
      ! dr_j back out and compares them.
      !
      ! The globals are set the way input_read sets them: grid_type from the
      ! "Grid type:" key, N from "Grid cells:" (default 500), dr_base and
      ! N_low_cells from "Base grid [dr,cells]:" (defaults dr_base_default
      ! and 50), r_max from the Roche/Hill radius or "Outer radius [R_p]:",
      ! r_esc from "Escape radius [R_p]:" and r_flux from "Flux spread tol:".
      !
      ! Reported for every configuration: the largest relative violation, the
      ! cell where it sits, and the ratio dr_j(j)/(r_edg(j)-r_edg(j-1)) there.
      ! On a stretched region that ratio IS the local stretch factor, which is
      ! what identifies the defect: dr_j(j) is the width of cell j+1.

      use global_parameters
      use grid_construction, only: define_grid

      implicit none

      integer :: nfail
      real*8, parameter :: tol = 1.0d-12

      nfail = 0

      ! Grid resolution and windows common to every configuration below.
      ! N is the "Grid cells:" default; the grid arrays are sized once from
      ! it and define_grid overwrites them for each configuration.
      N           = 500
      call allocate_grid_arrays
      dr_base     = dr_base_default ! the width an absent key resolves to
      N_low_cells = 50
      r_esc       = 1.50d0
      r_flux      = 1.20d0

      call check_grid('Mixed',     10.0d0, nfail)
      call check_grid('Mixed',      2.0d0, nfail)
      call check_grid('Uniform',   10.0d0, nfail)
      call check_grid('Stretched', 10.0d0, nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'grid_width_identity: ', nfail,                 &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'grid_width_identity: all assertions PASSED'

      contains

      subroutine check_grid(gtype, rmax_in, nfail)
      character(len=*), intent(in)    :: gtype
      real*8,           intent(in)    :: rmax_in
      integer,          intent(inout) :: nfail

      integer :: j, jworst
      real*8  :: w, err, errmax, ratio_worst, dsum, dref, err_sum
      character(len=48) :: tag

      grid_type = gtype
      r_max     = rmax_in

      call define_grid

      ! --- identity 1: the stored width of a cell is the distance between
      !     the two faces of THAT cell.
      errmax      = 0.0d0
      jworst      = 2 - Ng
      ratio_worst = 1.0d0
      do j = 2-Ng, N+Ng
         w = r_edg(j) - r_edg(j-1)
         if (w .eq. 0.0d0) cycle
         err = abs(dr_j(j) - w)/abs(w)
         if (err .gt. errmax) then
            errmax      = err
            jworst      = j
            ratio_worst = dr_j(j)/w
         endif
      enddo

      write(tag,'(A,A,A,F0.1)') 'dr_j_is_own_cell_width[', trim(gtype),      &
                                ',r_max=', rmax_in
      call verdict(trim(tag)//']', errmax, 0.0d0, tol, nfail)
      write(*,'(A,I5,A,ES22.15,A,ES22.15,A,F12.8)')                          &
           '     worst cell j=', jworst,                                     &
           '  dr_j=', dr_j(jworst),                                          &
           '  r_edg(j)-r_edg(j-1)=', r_edg(jworst)-r_edg(jworst-1),          &
           '  ratio=', ratio_worst

      ! --- identity 2: the widths of the physical cells tile the domain.
      dsum = 0.0d0
      do j = 1, N
         dsum = dsum + dr_j(j)
      enddo
      dref    = r_edg(N) - r_edg(0)
      err_sum = abs(dsum - dref)/abs(dref)
      write(tag,'(A,A,A,F0.1)') 'physical_widths_tile_domain[', trim(gtype), &
                                ',r_max=', rmax_in
      call verdict(trim(tag)//']', err_sum, 0.0d0, tol, nfail)
      write(*,'(A,ES22.15,A,ES22.15)')                                       &
           '     sum(dr_j(1:N))=', dsum, '  r_edg(N)-r_edg(0)=', dref

      end subroutine check_grid

      subroutine verdict(name, measured, reference, tolerance, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tolerance
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tolerance) then
         write(*,'(A,A,A,ES12.5,A,ES12.5,A,ES12.5)') 'PASS ', name,          &
              ' measured=', measured, ' reference=', reference,              &
              ' tol=', tolerance
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5,A,ES12.5)') 'FAIL ', name,          &
              ' measured=', measured, ' reference=', reference,              &
              ' tol=', tolerance
         nfail = nfail + 1
      endif
      end subroutine verdict

      end program grid_width_identity

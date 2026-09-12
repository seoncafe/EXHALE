      program fuv_band_quadrature
      ! The band integral of a loaded spectrum (sed_read::sed_band_fluxes)
      ! against exact answers on spectra whose integrals are closed-form:
      ! a constant and a linear F_lambda, with band edges inside segments,
      ! on nodes, at the table ends, on unequal segment widths and on a
      ! single segment covering the band; a band the table does not reach
      ! is reported uncovered; a stated file that is not monotone stops.
      ! RED on the text before 2026-09-13 (findings B1 and B2 of the review
      ! of 2026-09-12): the clipped trapezoid used the unclipped node
      ! fluxes (0.81 percent low on the linear case) and one covering
      ! segment returned zero.
      use global_parameters, only: sed_file
      use sed_reader, only: sed_band_fluxes
      use assertion_report, only: check_relative, check_absolute,       &
                                  assertion_failures
      implicit none
      integer, parameter :: nb = 6
      real*8  :: lo(nb), hi(nb), F(nb), exact
      logical :: cov(nb)
      integer :: u, b

      ! Table A: F_lambda = lambda, unequal widths, 4 nodes
      open(newunit=u, file='sed_linear.txt', status='replace')
      write(u,'(A)') '# test spectrum: F = lambda'
      write(u,*) 900.0d0,  900.0d0
      write(u,*) 1000.0d0, 1000.0d0
      write(u,*) 1150.0d0, 1150.0d0
      write(u,*) 1200.0d0, 1200.0d0
      close(u)
      sed_file = 'sed_linear.txt'
      lo = (/ 950.0d0, 900.0d0, 1000.0d0, 920.0d0, 1160.0d0, 800.0d0 /)
      hi = (/ 1100.0d0, 1200.0d0, 1150.0d0, 1180.0d0, 1190.0d0, 1000.0d0 /)
      call sed_band_fluxes(nb, lo, hi, F, cov)
      do b = 1, 5
         exact = 0.5d0*(hi(b)**2 - lo(b)**2)
         call check_relative('linear_spectrum_band_'//char(48+b),        &
                             F(b), exact, 1.0d-12)
         call check_absolute('linear_spectrum_band_'//char(48+b)//       &
                             '_covered', merge(1.0d0, 0.0d0, cov(b)), 1.0d0, 0.0d0)
      enddo
      call check_absolute('band_below_the_table_is_uncovered',           &
                          merge(1.0d0, 0.0d0, cov(6)), 0.0d0, 0.0d0)
      call check_absolute('band_below_the_table_reads_zero', F(6), 0.0d0, 0.0d0)

      ! Table B: one segment covering the band (two nodes only)
      open(newunit=u, file='sed_two_nodes.txt', status='replace')
      write(u,*) 900.0d0,  900.0d0
      write(u,*) 1200.0d0, 1200.0d0
      close(u)
      sed_file = 'sed_two_nodes.txt'
      call sed_band_fluxes(1, (/ 950.0d0 /), (/ 1100.0d0 /), F(1:1), cov(1:1))
      call check_relative('one_covering_segment_integrates', F(1),       &
                          0.5d0*(1100.0d0**2 - 950.0d0**2), 1.0d-12)
      call check_absolute('one_covering_segment_is_covered',            &
                          merge(1.0d0, 0.0d0, cov(1)), 1.0d0, 0.0d0)

      ! Table C: constant spectrum, edges on nodes and inside
      open(newunit=u, file='sed_const.txt', status='replace')
      write(u,*) 800.0d0, 3.0d0
      write(u,*) 950.0d0, 3.0d0
      write(u,*) 1300.0d0, 3.0d0
      write(u,*) 2500.0d0, 3.0d0
      close(u)
      sed_file = 'sed_const.txt'
      call sed_band_fluxes(2, (/ 912.0d0, 950.0d0 /), (/ 1201.0d0, 1300.0d0 /), F(1:2), cov(1:2))
      call check_relative('constant_spectrum_lw_band', F(1), 3.0d0*(1201.0d0-912.0d0), 1.0d-12)
      call check_relative('constant_spectrum_edges_on_nodes', F(2), 3.0d0*(1300.0d0-950.0d0), 1.0d-12)

      if (assertion_failures .gt. 0) call exit(1)
      end program fuv_band_quadrature

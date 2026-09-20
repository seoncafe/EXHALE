      program turnover_cell_rates
      ! WHICH CELL THE CHARGE-EXCHANGE TURNOVER BOUND BELONGS TO.
      !
      ! The bound cx_add_to_turnover adds to a residual row is
      !     kc(T) N(donor element) N(acceptor element)
      ! with kc the rate coefficient of the cell whose temperature
      ! cx_set_cell filled cx_kc at. Two statements are asserted here.
      !
      !   (1) It is the loaded cell's bound and carries no history. The
      !       rate coefficients are temperature dependent, so alternating
      !       two cells on one thread must alternate the bound, and coming
      !       back to the first temperature must reproduce its bound to the
      !       last bit: the routine reads cx_kc and nothing it accumulated
      !       on the way.
      !
      !   (2) Each thread answers for its own cell. cx_kc, cx_cell_T and
      !       cx_metal_base are thread-local, so a parallel region in which
      !       every thread loads a different cell must give every thread the
      !       bound of the cell it loaded, bitwise equal to the serial value
      !       for that temperature. The refusal of a cell that was never
      !       loaded is the subject of stale_cell_turnover_refused.
      !
      ! The two temperatures span the range the metal charge-exchange fits
      ! are used over in the wind, 6000 K and 10000 K.
      use species_table,    only: iel_Si, iel_O
      use charge_exchange,  only: cx_full, cx_init, cx_set_cell,          &
                                  cx_metal_base, cx_add_to_turnover
      use assertion_report, only: check_positive, check_absolute,         &
                                  assertion_failures
!$    use omp_lib,          only: omp_get_max_threads
      implicit none
      integer, parameter :: neq = 32, ncell = 8
      real*8, parameter  :: T_cool = 6.0d3, T_warm = 1.0d4
      real*8 :: el_tot(12)
      real*8 :: s_cool(neq), s_warm(neq), s_again(neq)
      real*8 :: s_thread(neq,ncell), Tc(ncell)
      real*8 :: dev, spread_rel
      integer :: i, nthreads

      cx_full = .true.
      call cx_init()
      el_tot         = 0.0d0
      el_tot(iel_Si) = 1.0d6
      el_tot(iel_O)  = 1.0d7
      el_tot(11)     = 1.0d10          ! cx_H
      el_tot(12)     = 1.0d9           ! cx_He
      cx_metal_base  = 5

      ! (1) one thread, alternating cells.
      call cx_set_cell(T_cool)
      s_cool = 0.0d0
      call cx_add_to_turnover(s_cool, el_tot, T_cool)

      call cx_set_cell(T_warm)
      s_warm = 0.0d0
      call cx_add_to_turnover(s_warm, el_tot, T_warm)

      call cx_set_cell(T_cool)
      s_again = 0.0d0
      call cx_add_to_turnover(s_again, el_tot, T_cool)

      ! The two cells must be told apart at all: a bound that did not move
      ! with the temperature would make the assertion below vacuous.
      spread_rel = maxval(abs(s_warm - s_cool))/maxval(abs(s_cool))
      call check_positive('turnover_moves_with_the_cell_temperature',     &
                          spread_rel)
      dev = maxval(abs(s_again - s_cool))
      call check_absolute('turnover_of_a_reloaded_cell_is_its_own',       &
                          dev, 0.0d0, 0.0d0)

      ! (2) several threads, a different cell on each. Every thread loads
      ! its own cell inside the region, because cx_kc, cx_cell_T and
      ! cx_metal_base are thread-local and a worker entering the region
      ! holds no cell.
      do i = 1, ncell
         Tc(i) = T_cool
         if (mod(i,2) .eq. 0) Tc(i) = T_warm
      enddo
      s_thread = 0.0d0
      nthreads = 1
!$    nthreads = omp_get_max_threads()

!$omp parallel do default(shared) private(i) schedule(static,1)
      do i = 1, ncell
         cx_metal_base = 5
         call cx_set_cell(Tc(i))
         call cx_add_to_turnover(s_thread(:,i), el_tot, Tc(i))
      enddo
!$omp end parallel do

      dev = 0.0d0
      do i = 1, ncell
         if (Tc(i) .eq. T_cool) then
            dev = max(dev, maxval(abs(s_thread(:,i) - s_cool)))
         else
            dev = max(dev, maxval(abs(s_thread(:,i) - s_warm)))
         endif
      enddo
      write(*,'(a,i0,a)') '  the parallel region ran on ', nthreads,      &
                          ' thread(s)'
      call check_absolute('turnover_on_several_threads_is_each_cells',    &
                          dev, 0.0d0, 0.0d0)

      if (assertion_failures .gt. 0) error stop 1
      end program turnover_cell_rates

      program stale_cell_stage_sources_refused
      ! The same cell ownership for the stage-source interface the
      ! transported ionization stages take their charge exchange through.
      ! The rate coefficients are thread-local and hold one cell, so a
      ! caller that evaluates a cell whose temperature is not the one
      ! cx_set_cell filled them at would be using whichever cell the thread
      ! handled last, and the source of a cell would depend on the order the
      ! cells were handed out. charge_exchange_stage_sources refuses it.
      !
      ! HOW TO READ THIS PROGRAM: it is expected to STOP with a nonzero
      ! status. The driver script turns a clean exit into a failing verdict.
      ! The first two evaluations must succeed, so that the abort is shown
      ! to come from the temperature mismatch alone.
      use species_table,   only: n_melem, iel_Si
      use charge_exchange, only: cx_full, cx_init, cx_set_cell,           &
                                 charge_exchange_stage_sources
      implicit none
      real*8 :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8 :: sH(0:1), sHe(0:2)

      cx_full = .true.
      call cx_init()
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      nm0(iel_Si) = 1.0d0

      call cx_set_cell(1.0d4)
      call charge_exchange_stage_sources(nm0, nm1, nm2,                   &
                       0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d4, sH, sHe)
      write(*,'(a)') '  the cell that was loaded is accepted'

      call cx_set_cell(6.0d3)
      call charge_exchange_stage_sources(nm0, nm1, nm2,                   &
                       0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 6.0d3, sH, sHe)
      write(*,'(a)') '  the next cell that was loaded is accepted'

      ! The cell of the previous call, evaluated again without reloading
      ! its rates: this must stop the run.
      call charge_exchange_stage_sources(nm0, nm1, nm2,                   &
                       0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d4, sH, sHe)
      write(*,'(a)') '  REACHED: the stale-cell evaluation was not refused'
      end program stale_cell_stage_sources_refused

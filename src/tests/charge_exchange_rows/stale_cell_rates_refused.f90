      program stale_cell_rates_refused
      ! The charge-exchange rate coefficients belong to one cell at a time.
      ! A caller that assembles a cell whose temperature is not the one
      ! cx_set_cell filled cx_kc at would be using whichever cell the thread
      ! handled last, and its source would depend on the order the cells
      ! were handed out. cx_add_to_fvec refuses that call.
      !
      ! HOW TO READ THIS PROGRAM: it is expected to STOP with a nonzero
      ! status. The driver script turns a clean exit into a failing verdict.
      ! The first two assemblies must succeed, so that the abort is shown to
      ! come from the temperature mismatch alone and not from the module
      ! being unusable.
      use species_table,   only: n_melem, iel_Si
      use charge_exchange, only: cx_full, cx_init, cx_set_cell,           &
                                 cx_metal_base, cx_add_to_fvec
      implicit none
      integer, parameter :: neq = 32
      real*8 :: nm0(n_melem), nm1(n_melem), nm2(n_melem), fvec(neq)

      cx_full = .true.
      call cx_init()
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      nm0(iel_Si) = 1.0d0
      cx_metal_base = 5

      call cx_set_cell(1.0d4)
      fvec = 0.0d0
      call cx_add_to_fvec(neq, fvec, nm0, nm1, nm2,                       &
                          0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,              &
                          -1.0d0, 1.0d4)
      write(*,'(a)') '  the cell that was loaded is accepted'

      call cx_set_cell(6.0d3)
      fvec = 0.0d0
      call cx_add_to_fvec(neq, fvec, nm0, nm1, nm2,                       &
                          0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,              &
                          -1.0d0, 6.0d3)
      write(*,'(a)') '  the next cell that was loaded is accepted'

      ! The cell of the previous call, assembled again without reloading
      ! its rates: this must stop the run.
      fvec = 0.0d0
      call cx_add_to_fvec(neq, fvec, nm0, nm1, nm2,                       &
                          0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,              &
                          -1.0d0, 1.0d4)
      write(*,'(a)') '  REACHED: the stale-cell assembly was not refused'
      end program stale_cell_rates_refused

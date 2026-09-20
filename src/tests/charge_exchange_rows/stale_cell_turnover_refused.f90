      program stale_cell_turnover_refused
      ! The turnover bound of a residual row is kc(T) times two element
      ! counts, so it belongs to the cell whose temperature filled cx_kc.
      ! A caller that states another cell would normalize that cell's
      ! reaction residual by this one's turnover, and the acceptance verdict
      ! read from it would depend on the order the cells were handed out.
      ! cx_add_to_turnover refuses that call, as cx_add_to_fvec and
      ! cx_add_to_jac refuse theirs.
      !
      ! HOW TO READ THIS PROGRAM: it is expected to STOP with a nonzero
      ! status. The driver script turns a clean exit into a failing verdict.
      ! The first two bounds must be taken without complaint, so that the
      ! abort is shown to come from the temperature mismatch alone and not
      ! from the module being unusable.
      use species_table,   only: iel_Si
      use charge_exchange, only: cx_full, cx_init, cx_set_cell,           &
                                 cx_metal_base, cx_add_to_turnover
      implicit none
      integer, parameter :: neq = 32
      real*8 :: srow(neq), el_tot(12)

      cx_full = .true.
      call cx_init()
      el_tot        = 0.0d0
      el_tot(iel_Si) = 1.0d6
      el_tot(11)    = 1.0d10          ! cx_H
      el_tot(12)    = 1.0d9           ! cx_He
      cx_metal_base = 5

      call cx_set_cell(1.0d4)
      srow = 0.0d0
      call cx_add_to_turnover(srow, el_tot, 1.0d4)
      write(*,'(a)') '  the cell that was loaded is accepted'

      call cx_set_cell(6.0d3)
      srow = 0.0d0
      call cx_add_to_turnover(srow, el_tot, 6.0d3)
      write(*,'(a)') '  the next cell that was loaded is accepted'

      ! The cell of the previous call, judged again without reloading its
      ! rates: this must stop the run.
      srow = 0.0d0
      call cx_add_to_turnover(srow, el_tot, 1.0d4)
      write(*,'(a)') '  REACHED: the stale-cell turnover was not refused'
      end program stale_cell_turnover_refused

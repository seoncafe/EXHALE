      program carrier_returned_state_acceptance
      ! WHAT MAKES A RETURNED CARRIER STATE A SOLUTION, AND WHAT DOES NOT.
      !
      ! The carrier transport-chemistry solve hands a state back and that
      ! state is written into the run: it becomes the free atomic hydrogen
      ! the next ionization sweep closes on, the chemical heat of the
      ! molecular layer and the H2 opacity.  So the question the code has to
      ! answer before the write-back is whether THAT state satisfies the
      ! carrier balances, and this suite states the answer on constructed
      ! rows, where the residual of every row is known exactly.
      !
      ! THREE THINGS THAT ARE NOT THE ANSWER, each measured below against a
      ! row set built to separate it from the answer:
      !
      !  (1) how the iteration ended.  The element limiter runs AFTER the
      !      iteration, and it rescales whole carrier families in the cells
      !      it binds, so a state that leaves a converged iteration is not
      !      the state the iteration measured.
      !
      !  (2) whether the cell holding the LARGEST returned residual was
      !      constrained.  A constraint acting in one cell says something
      !      about that cell and nothing about any other, so that test lets
      !      a large constrained residual excuse a smaller unacceptable one
      !      anywhere else in the column.  This is the masking condition,
      !      and the reference predicate accepts_by_status_and_worst_cell
      !      below reproduces it so that the two verdicts can be compared in
      !      one run rather than asserted apart.
      !
      !  (3) how far the residual FELL.  A drop of six decades from a start
      !      of 1e8 leaves 1e2, and a norm of 1e2 is not an accuracy
      !      statement.  carrier_stall_class separates that from reaching
      !      the absolute floor, and only the floor is an accuracy
      !      statement.
      !
      ! THE MEASURE, and it is the one carrier_residual and
      ! carrier_steady_residual already use: the row residual divided by
      ! that row's own scale, the sum of the magnitudes of the terms the row
      ! balances plus an absolute floor tied to the free density of the
      ! carrier's element.  Every row scale here is 1, so the constructed
      ! residuals ARE the floored relative imbalances and no assertion
      ! depends on a scale this suite invents.

      use global_parameters, only: thereis_oxychem, ionization_transport, &
                                   thereis_mol, carrier_transport
      use diffusive_photochemistry, only:                                 &
              carrier_returned_state_verdict, carrier_stall_class,        &
              carrier_verdict, carrier_set_init,                          &
              carrier_accept_ok, carrier_reject_row_residual,             &
              carrier_reject_nonfinite,                                   &
              carrier_stall_none, carrier_stall_at_floor,                 &
              carrier_stall_at_relative_drop, carrier_stall_open,         &
              carrier_solve_converged, carrier_solve_iteration_cap,       &
              carrier_solve_stagnated, carrier_solve_line_search_failed,  &
              n_carrier_max, ic_H2, ic_OH, ic_H2O, ic_CO
      use assertion_report
      implicit none

      integer, parameter :: ncell = 6
      ! The state is accepted below this floored relative row imbalance; it
      ! is newton_floor, the absolute residual the solve declares to be what
      ! the discretization resolves.  Stated here as a number so that a
      ! change to it shows up as a failure and not as a silent retune.
      real*8,  parameter :: accept_tol = 1.0d-8

      real*8  :: res(ncell,n_carrier_max), terms(ncell,n_carrier_max)
      logical :: constrained(ncell)
      type(carrier_verdict) :: v
      real*8  :: nan, inf

      nan = transfer(int(z'7FF8000000000000', kind=8), 1.0d0)
      inf = transfer(int(z'7FF0000000000000', kind=8), 1.0d0)

      ! The molecular network with its carriers transported and the oxygen
      ! cycle on, so H2, OH, H2O and CO all carry an equation.  A carrier
      ! is a row of the transport-chemistry operator where the gas has it
      ! AND the run transports it, so the configuration is stated in full
      ! and not through the oxygen cycle alone.
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .true.
      ionization_transport = .false.
      call carrier_set_init()

      ! ---- (1) a converged iteration, then a limiter that leaves a bad
      !          unconstrained row ---------------------------------------!
      ! The limiter acted in cell 2.  Cell 4 it did not touch, and cell 4
      ! does not satisfy its H2 balance.
      call blank_rows()
      constrained(2) = .true.
      res(4,ic_H2)   = 1.0d-3
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('converged_then_limiter_leaves_bad_row_refused',&
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)
      call check_absolute('converged_then_limiter_reason_is_row',         &
                          dble(v%reason),                                 &
                          dble(carrier_reject_row_residual), 0.0d0)
      call check_absolute('converged_then_limiter_names_the_bad_cell',    &
                          dble(v%jworst), 4.0d0, 0.0d0)
      call check_absolute('converged_then_limiter_names_the_carrier',     &
                          dble(v%icworst), dble(ic_H2), 0.0d0)
      call check_relative('converged_then_limiter_reports_the_imbalance', &
                          v%rworst, 1.0d-3, 1.0d-14)
      ! The raw status is kept, and it is not the decision.
      call check_absolute('converged_then_limiter_keeps_the_raw_status',  &
                          dble(v%raw_status),                             &
                          dble(carrier_solve_converged), 0.0d0)
      call check_absolute('converged_then_limiter_old_rule_accepted',     &
                          logical_as_double(                              &
                             accepts_by_status_and_worst_cell(            &
                                carrier_solve_converged)), 1.0d0, 0.0d0)

      ! ---- (2) the masking condition: a failed iteration whose worst cell
      !          is constrained, and a smaller bad row elsewhere ---------!
      ! Cell 2 was clamped and carries an imbalance of 10; cell 5, which no
      ! constraint touched, carries 1e-4.  Asking only about the worst cell
      ! finds cell 2, finds it constrained, and certifies the column.
      call blank_rows()
      constrained(2) = .true.
      res(2,ic_CO)   = 1.0d1
      res(5,ic_OH)   = 1.0d-4
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_iteration_cap,    &
                                          carrier_stall_open, v)
      call check_absolute('masked_unconstrained_row_refused',             &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)
      call check_absolute('masked_row_named_by_cell',                     &
                          dble(v%jworst), 5.0d0, 0.0d0)
      call check_absolute('masked_row_named_by_carrier',                  &
                          dble(v%icworst), dble(ic_OH), 0.0d0)
      call check_relative('masked_row_imbalance_is_the_unconstrained_one',&
                          v%rworst, 1.0d-4, 1.0d-14)
      call check_absolute('masked_row_is_counted_once',                   &
                          dble(v%n_rows_above), 1.0d0, 0.0d0)
      call check_absolute('masked_row_is_not_beside_a_clamped_cell',      &
                          dble(v%n_rows_above_beside_constrained),        &
                          0.0d0, 0.0d0)
      call check_absolute('masked_row_old_rule_accepted',                 &
                          logical_as_double(                              &
                             accepts_by_status_and_worst_cell(            &
                                carrier_solve_iteration_cap)),            &
                          1.0d0, 0.0d0)

      ! A bad row in the cell NEXT to the clamped one is counted as such:
      ! a clamp moves its own cell's carrier densities and those enter the
      ! face fluxes of both its neighbours, so the two cases are separated
      ! by counting and not by argument.
      call blank_rows()
      constrained(2) = .true.
      res(3,ic_H2O)  = 1.0d-4
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_bad_row_beside_a_clamped_cell_is_counted',   &
                          dble(v%n_rows_above_beside_constrained),        &
                          1.0d0, 0.0d0)
      call check_absolute('a_bad_row_beside_a_clamp_is_still_refused',    &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)

      ! ---- (3) a large relative drop is not an accuracy statement -----!
      ! The counterexample of review R6: an iteration that started at 1e8,
      ! stood at 1e2 one step ago and stands at 60 now.  It did not halve,
      ! and it is six decades below its start.
      call check_absolute('relative_drop_is_not_the_floor',               &
                          dble(carrier_stall_class(3, 6.0d1, 1.0d2,       &
                                                   1.0d8)),               &
                          dble(carrier_stall_at_relative_drop), 0.0d0)
      call check_absolute('absolute_floor_is_the_floor',                  &
                          dble(carrier_stall_class(3, 1.0d-9, 1.5d-9,     &
                                                   1.0d8)),               &
                          dble(carrier_stall_at_floor), 0.0d0)
      call check_absolute('a_stuck_iteration_is_neither',                 &
                          dble(carrier_stall_class(9, 1.0d-2, 1.2d-2,     &
                                                   1.0d0)),               &
                          dble(carrier_stall_open), 0.0d0)
      call check_absolute('a_halving_iteration_is_not_stalling',          &
                          dble(carrier_stall_class(9, 1.0d-3, 1.0d0,      &
                                                   1.0d8)),               &
                          dble(carrier_stall_none), 0.0d0)
      call check_absolute('too_early_to_call_it_a_stall',                 &
                          dble(carrier_stall_class(2, 6.0d1, 1.0d2,       &
                                                   1.0d8)),               &
                          dble(carrier_stall_none), 0.0d0)
      ! And the state such an iteration returns is refused on its rows.
      call blank_rows()
      res(3,ic_H2) = 6.0d1
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_stagnated,        &
                                          carrier_stall_at_relative_drop, &
                                          v)
      call check_absolute('state_after_a_relative_drop_refused',          &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)
      call check_relative('state_after_a_relative_drop_imbalance',        &
                          v%rworst, 6.0d1, 1.0d-14)

      ! ---- (4) a non-finite row, constrained or not -------------------!
      call blank_rows()
      res(3,ic_H2O) = nan
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_NaN_row_is_refused',                         &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)
      call check_absolute('a_NaN_row_says_why',                           &
                          dble(v%reason), dble(carrier_reject_nonfinite), &
                          0.0d0)
      call blank_rows()
      constrained(3) = .true.
      res(3,ic_H2O)  = nan
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_constraint_does_not_certify_a_NaN',          &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)
      call blank_rows()
      terms(2,ic_CO) = inf
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('an_infinite_row_scale_is_refused',             &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)

      ! ---- (5) a constrained state that is otherwise a solution -------!
      ! Cell 2 was clamped and its four rows stand at 5; every other row in
      ! the column is at 1e-12.  The state is accepted, and the verdict
      ! says how much of it the constraint holds rather than the solver.
      call blank_rows()
      constrained(2) = .true.
      res            = 1.0d-12
      res(2,:)       = 5.0d0
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_constrained_cell_does_not_refuse_the_state', &
                          logical_as_double(v%accepted), 1.0d0, 0.0d0)
      call check_absolute('the_constrained_cells_are_reported',           &
                          dble(v%n_cells_constrained), 1.0d0, 0.0d0)
      call check_absolute('the_constrained_rows_are_reported',            &
                          dble(v%n_rows_constrained), 4.0d0, 0.0d0)
      call check_relative('the_worst_certified_row_is_the_solver'//''''// &
                          's', v%rworst, 1.0d-12, 1.0d-14)

      ! ---- (6) a state that solves its equations ----------------------!
      call blank_rows()
      res = 1.0d-13
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_solved_state_is_accepted',                   &
                          logical_as_double(v%accepted), 1.0d0, 0.0d0)
      call check_absolute('a_solved_state_says_so',                       &
                          dble(v%reason), dble(carrier_accept_ok), 0.0d0)
      ! The boundary of the condition, from both sides.
      call blank_rows()
      res(4,ic_OH) = accept_tol
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_row_at_the_floor_is_accepted',               &
                          logical_as_double(v%accepted), 1.0d0, 0.0d0)
      call blank_rows()
      res(4,ic_OH) = accept_tol*1.0d1
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('a_row_a_decade_above_the_floor_is_refused',    &
                          logical_as_double(v%accepted), 0.0d0, 0.0d0)

      ! ---- (7) a carrier this run does not solve carries no equation --!
      thereis_oxychem = .false.
      call carrier_set_init()
      call blank_rows()
      res(3,ic_OH) = 1.0d3
      call carrier_returned_state_verdict(res, terms, constrained,        &
                                          carrier_solve_converged,        &
                                          carrier_stall_at_floor, v)
      call check_absolute('an_unsolved_carrier_row_is_not_an_equation',   &
                          logical_as_double(v%accepted), 1.0d0, 0.0d0)
      thereis_oxychem = .true.
      call carrier_set_init()

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      subroutine blank_rows()
      ! A column that solves every one of its rows exactly, with a unit
      ! scale on each, so that a residual written into it is read back as
      ! its own floored relative imbalance.
      res         = 0.0d0
      terms       = 1.0d0
      constrained = .false.
      end subroutine blank_rows

      !--------------!

      logical function accepts_by_status_and_worst_cell(status)
      ! THE CONDITION THIS SUITE REPLACES, kept here so that the two
      ! verdicts are compared on one row set in one run: the state passes
      ! if the iteration reported convergence, or if the cell holding the
      ! largest residual of the whole column was constrained.  It is stated
      ! on the same arrays the verdict reads.
      integer, intent(in) :: status
      integer :: j, ic, jw
      real*8  :: rel, rmax
      accepts_by_status_and_worst_cell = .true.
      if (status .eq. carrier_solve_converged) return
      rmax = 0.0d0
      jw   = 0
      do j = 1, ncell
         do ic = 1, n_carrier_max
            rel = abs(res(j,ic))/max(terms(j,ic), 1.0d-300)
            if (rel .gt. rmax) then
               rmax = rel
               jw   = j
            endif
         enddo
      enddo
      accepts_by_status_and_worst_cell = .false.
      if (jw .ge. 1) accepts_by_status_and_worst_cell = constrained(jw)
      end function accepts_by_status_and_worst_cell

      !--------------!

      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double

      end program carrier_returned_state_acceptance

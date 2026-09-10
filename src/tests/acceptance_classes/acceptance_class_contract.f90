      program acceptance_class_contract
      ! ONE MEANING OF THE IONIZATION ACCEPTANCE CLASSES FOR EVERY CONSUMER.
      !
      ! ioniz_eq classifies each accepted cell state (ionization_equilibrium.
      ! f90, the acceptance block of each branch):
      !   1 solver-converged root, 2 root without solver convergence,
      !   3 projected/handback state whose rechecked residual still marks a
      !   root, 5 root of the constrained element-conserving continuation --
      !   all four certified by the same two tests, a state inside the
      !   element simplex and a normalized reaction residual at or below
      !   ieq_res_tol -- and
      !   4 a state that is NOT a root, kept only by the relaxation amnesty,
      !   6 a candidate that was NOT a root and whose residual was above
      !     ieq_nonroot_res_cap, so the cell kept the composition it entered
      !     the sweep with -- a state that solved an EARLIER cell state and
      !     not this one, and therefore not a root of this one either.
      !
      ! What is asserted here:
      !   (a) the count of cells whose composition is not a root of the
      !       network is classes 4 and 6 and nothing else, so equivalent
      !       roots reported through class 1 and through class 5 give the
      !       same verdict, one class-4 cell gives one, and one class-6 cell
      !       gives one as well;
      !   (b) the FINAL acceptance predicate (steady_gates_met) refuses a
      !       state carrying a class-4 or a class-6 cell and accepts one
      !       whose cells are all roots, whatever class reported them.
      !
      ! The predicate is exercised with the flux gate disabled
      ! (flux_spread_th <= 0) and the flux window placed outside the grid, so
      ! that the residual, carrier and chemistry conditions are what is being
      ! measured and no Riemann solve is needed.
      !
      ! Each assertion prints
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! and the exit status is nonzero if any of them fails.
      use global_parameters
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use steady_residual_mod, only: n_cells_without_chemical_root,      &
                                     steady_gates_met
      use ionization_equilibrium, only: nonroot_acceptance_class,        &
                                        ieq_nonroot_res_cap,             &
                                        ieq_res_tol
      implicit none

      integer :: n_fail
      integer :: acc_root1(6), acc_root5(6), acc_mixed(6), acc_one4(6)
      integer :: acc_many4(6), acc_one6(6), acc_mix46(6)

      n_fail = 0

      ! (a) the counting function. Populations of one sweep, by class.
      acc_root1 = (/ 500, 0, 0, 0,   0, 0 /) ! every cell a converged root
      acc_root5 = (/   0, 0, 0, 0, 500, 0 /) ! the same roots, found by the
                                             ! constrained continuation
      acc_mixed = (/ 300, 90, 60, 0, 50, 0 /)! roots through all four routes
      acc_one4  = (/ 499, 0, 0, 1,   0, 0 /) ! one cell without a root
      acc_many4 = (/ 300, 90, 57, 3, 50, 0 /)! three of them
      acc_one6  = (/ 499, 0, 0, 0,   0, 1 /) ! one cell that kept its entry
                                             ! composition above the cap
      acc_mix46 = (/ 297, 90, 57, 3, 50, 3 /)! three of each

      call check_int('count_class1_roots',                                &
                     n_cells_without_chemical_root(acc_root1), 0)
      call check_int('count_class5_roots',                                &
                     n_cells_without_chemical_root(acc_root5), 0)
      call check_int('count_roots_all_classes',                           &
                     n_cells_without_chemical_root(acc_mixed), 0)
      call check_int('count_one_class4_cell',                             &
                     n_cells_without_chemical_root(acc_one4), 1)
      call check_int('count_three_class4_cells',                          &
                     n_cells_without_chemical_root(acc_many4), 3)
      call check_int('count_one_class6_cell',                             &
                     n_cells_without_chemical_root(acc_one6), 1)
      call check_int('count_class4_and_class6_together',                  &
                     n_cells_without_chemical_root(acc_mix46), 6)
      ! Class 6 is not a milder class 4: one cell of either is one cell
      ! the chemistry does not describe.
      call check_int('class4_and_class6_same_verdict',                    &
                     n_cells_without_chemical_root(acc_one6),             &
                     n_cells_without_chemical_root(acc_one4))
      ! The equivalence the contract states, asserted as such: the same
      ! number of roots is the same verdict whether class 1 or class 5
      ! reported them.
      call check_int('class1_and_class5_same_verdict',                    &
                     n_cells_without_chemical_root(acc_root5),            &
                     n_cells_without_chemical_root(acc_root1))

      ! (a2) the rule that splits the two non-root classes, at its single
      ! definition. Below the cap the candidate is adopted under the
      ! amnesty (class 4); above it the cell keeps the composition it
      ! entered the sweep with (class 6); a residual that is not a finite
      ! number falls on the class-6 side, because a state that has driven a
      ! balance row out of the reals is not one to march on.
      call check_int('nonroot_just_above_the_root_tolerance_is_class4',   &
                     nonroot_acceptance_class(10.0d0*ieq_res_tol), 4)
      call check_int('nonroot_just_below_the_cap_is_class4',              &
                     nonroot_acceptance_class(0.99d0*ieq_nonroot_res_cap),&
                     4)
      call check_int('nonroot_at_the_cap_is_class4',                      &
                     nonroot_acceptance_class(ieq_nonroot_res_cap), 4)
      call check_int('nonroot_just_above_the_cap_is_class6',              &
                     nonroot_acceptance_class(1.01d0*ieq_nonroot_res_cap),&
                     6)
      call check_int('the_measured_poisoning_event_is_class6',            &
                     nonroot_acceptance_class(5.2d4), 6)
      call check_int('a_nonfinite_residual_is_class6',                    &
                     nonroot_acceptance_class(ieee_value(1.0d0,           &
                                              ieee_quiet_nan)), 6)

      ! (b) the final acceptance predicate.
      call check_gate()

      if (n_fail .gt. 0) then
         write(*,'(A,I0,A)') 'acceptance_class_contract: ', n_fail,       &
              ' assertion(s) failed'
         stop 1
      endif
      write(*,'(A)') 'acceptance_class_contract: every assertion passed'

      contains

      subroutine check_int(name, measured, reference)
      character(*), intent(in) :: name
      integer,      intent(in) :: measured, reference
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_int

      subroutine check_log(name, measured, reference)
      character(*), intent(in) :: name
      logical,      intent(in) :: measured, reference
      character(5) :: sm, sr
      sm = 'false';  sr = 'false'
      if (measured)  sm = 'true '
      if (reference) sr = 'true '
      if (measured .eqv. reference) then
         write(*,'(A,A,A,A,A,A,A)') 'PASS ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
      else
         write(*,'(A,A,A,A,A,A,A)') 'FAIL ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_log

      subroutine check_gate()
      ! A state that meets the residual gate by a wide margin, so that the
      ! chemistry condition is the only thing that can refuse it.
      real*8, dimension(:,:), allocatable :: u
      real*8  :: fspread, resid_tol, rnorm
      N              = 8
      j_flux         = N + 4        ! outside the grid: the flux functional
                                    ! returns zero without a Riemann solve
      flux_spread_th = -1.0d0       ! flux gate disabled
      allocate(u(3,1-Ng:N+Ng))
      u         = 1.0d0
      rnorm     = 1.0d-9
      resid_tol = 1.0d-5

      ! Roots only, reported through class 1, then through class 5, then
      ! through all four root classes: the state is certifiable in each case.
      call check_log('gate_accepts_class1_roots',                         &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_root1)), &
           .true.)
      call check_log('gate_accepts_class5_roots',                         &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_root5)), &
           .true.)
      call check_log('gate_accepts_roots_all_classes',                    &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_mixed)), &
           .true.)
      ! One cell without a chemical root, everything else identical: not a
      ! certified solution.
      call check_log('gate_refuses_one_class4_cell',                      &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_one4)), &
           .false.)
      call check_log('gate_refuses_three_class4_cells',                   &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_many4)), &
           .false.)
      call check_log('gate_refuses_one_class6_cell',                      &
           steady_gates_met(rnorm, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_one6)), &
           .false.)
      ! The residual gate still decides on its own: a rootless-free state
      ! that is far from balance is refused for the residual.
      call check_log('gate_refuses_large_residual',                       &
           steady_gates_met(1.0d0, u, resid_tol, fspread, .false.,        &
                            0.0d0, n_cells_without_chemical_root(acc_root1)), &
           .false.)
      deallocate(u)
      end subroutine check_gate

      end program acceptance_class_contract

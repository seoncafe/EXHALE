      program krylov_and_dogleg_tests
      ! WHAT THE LINEAR SOLVE AND THE STEP GEOMETRY OF THE STATIONARY SOLVER
      ! ARE ALLOWED TO RETURN (docs/PLAN_20260909_rev1.md section 4, items
      ! N0 to N3; docs/ISSUES_20260909_review.md R1, R2, R3).
      !
      ! Every row calls a PRODUCTION routine of steady_newton. The Krylov
      ! rows drive the production cycle on a dense operator with an identity
      ! preconditioner, installed through linear_operator_set_for_test, so
      ! that what is measured is the routine the solver runs and not a copy
      ! of its text. N, nvar_jac and the band half-widths are set here to
      ! the size of the small linear system each row states: pgmres reads
      ! them only as the length of its vectors and the shape of the banded
      ! factorization it is handed, and no species row is registered
      ! (set_transported_species_rows(.false.)), so the cycle takes its
      ! three-unknown branch. Every local carries a kd_ prefix, because a
      ! plain one-letter name here would shadow a global of the same letter.
      !
      ! WHAT EACH GROUP IS FOR.
      !   * a singular operator: a vanishing Arnoldi remainder is not by
      !     itself evidence that the system is solved, and the triangular
      !     solution may not divide by a zero diagonal. The reviewed source
      !     returned a NONFINITE step and reported relative residual zero
      !     for the operator zero with right-hand side one.
      !   * a rank-deficient reduced problem, compatible and incompatible:
      !     a finite least-squares step with its true residual reported, or
      !     a named refusal; never a nonfinite step.
      !   * an Arnoldi remainder that underflows: the norm of the remainder
      !     is not the evidence of convergence, the true residual is.
      !   * an unusable Jacobian action and a failing preconditioner: a
      !     named outcome, not arithmetic on a value that does not exist.
      !   * the step length along the banded gradient: the model must fall
      !     at the point returned. With residual one, operator one and a
      !     band ten the reviewed length raised the merit 0.5 to 40.5.
      !   * the dogleg with a radius of zero: a zero step, never a nonzero
      !     step inside a ball of radius zero. The geometry is not the
      !     defect; the caller's radius is (R1, F4).
      !   * the ray test's floor and what may overrule its verdict
      !     (ray_slope_verdict): a difference of two merit evaluations
      !     carries a verdict only above the reproducibility of the merit
      !     itself, and a slope quotient taken over a thousandth of the step
      !     may not refuse a step whose own predicted and actual reductions
      !     agree. The numbers of these rows are the refusing iteration of
      !     the atomic element solve (N7: merit 2.661e-5, model slope
      !     -9.132e-9, ray slope -1.483e-9, difference 2.97e-12 against an
      !     assumed floor of 2.661e-13, reduction ratio 1.0003).
      !   * one functional (decision 20 a): the merit the step control
      !     descends and the measure the state is judged by read every row
      !     through one scale (merit_row_scale_from_certification), so the
      !     row the merit is largest on and the row the gate reports are the
      !     same row, and an element row's share of the squared merit is its
      !     own certification measure squared.
      !   * the approximate-gradient leg: a leg shorter than
      !     cauchy_leg_min_fraction of the radius is not in the step, and
      !     then it is not in the model slope of the step either.
      !   * the trust radius: a ceiling that doubles on every good step has
      !     one absolute bound (trust_region_radius_absolute_cap), and it
      !     does not bind on an ordinary radius.
      !   * the outer iteration cap: a bounded diagnostic run has to stop
      !     where it is asked to, so EXHALE_JFNK_MAXIT replaces the caller's
      !     number and the resolved cap is what the solve prints and runs.
      !   * the restart of the trust-region state (item N21): the trigger
      !     reads the stagnation detector's own statistic, the number of
      !     consecutive outer iterations in which the best iterate did not
      !     improve, so it fires on a sequence of iterations that leave the
      !     iterate where it is and never on a sequence of improving steps;
      !     a restart resets exactly the parts of that state it names and
      !     touches no other; and a Krylov solve of one cycle is the single
      !     cycle the solver has always taken, while a second cycle of the
      !     same subspace lowers the residual the first left.
      !
      ! Every row prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line. Exit status is nonzero if any row fails.

      use global_parameters
      use steady_newton, only: pgmres, dogleg_step,                       &
                               cauchy_length_along_the_banded_gradient,   &
                               predicted_model_decrease,                  &
                               linear_operator_set_for_test,              &
                               linear_operator_clear_for_test,            &
                               lin_test_dense,                            &
                               lin_test_dense_preconditioner_fails,       &
                               set_transported_species_rows,              &
                               species_unknowns_outside_their_bounds,     &
                               hold_the_step_where_it_leaves_the_species_box, &
                               species_box_set_for_test,                  &
                               species_box_clear_for_test,                &
                               probe_length_of_the_jacobian_action,       &
                               nvar_jac, kl_jac, ku_jac,                  &
                               gm_tolerance_reached,                      &
                               gm_subspace_exhausted,                     &
                               gm_jacobian_action_unusable,               &
                               gm_preconditioner_failed,                  &
                               gm_reduced_system_singular,                &
                               gm_no_direction_sampled,                   &
                               initial_trust_region_radius,               &
                               tr_radius_from_the_krylov_step,            &
                               tr_radius_from_the_gradient,               &
                               tr_radius_no_feasible_descent,             &
                               ray_slope_verdict, ray_floor_factor,       &
                               merit_row_scale_from_certification,        &
                               trust_region_radius_absolute_cap,          &
                               grow_the_trust_region_ceiling,             &
                               cauchy_leg_min_fraction,                   &
                               cauchy_leg_image_share_max,                &
                               cauchy_leg_step_share_max,                 &
                               cauchy_leg_shares_of_step_and_image,       &
                               cauchy_leg_is_image_without_step,          &
                               cauchy_leg_is_dropped_on_its_length,       &
                               cauchy_leg_admitted_by_image_alone,        &
                               tr_ratio_below_acceptance,                 &
                               read_species_unknown_space_controls,       &
                               jfnk_outer_iteration_cap,                  &
                               tr_step_taken,                             &
                               tr_ray_probe_inadmissible,                 &
                               tr_ray_slopes_disagree_size,               &
                               tr_ray_slopes_disagree_sign,               &
                               tr_ray_below_the_noise_floor,              &
                               equilibrate_the_scaled_columns,            &
                               element_unknown_column_scale,              &
                               element_column_scale_floor,                &
                               equilibrate_the_preconditioner,            &
                               preconditioner_equilibrated,               &
                               precon_row_scale, precon_col_scale,        &
                               read_trust_region_restart_controls,        &
                               trust_region_restart_is_due,               &
                               restart_the_trust_region_state,            &
                               tr_reset_wanted, n_tr_reset_kind,          &
                               tr_reset_radius_and_ceiling,               &
                               tr_reset_pseudo_transient,                 &
                               tr_reset_closure_map,                      &
                               tr_reset_acceptance_memory,                &
                               tr_reset_back_to_best, n_tr_restart,       &
                               gm_restart_cycles,                         &
                               pgmres_with_restarts,                      &
                               gm_stopped_on_the_trust_ball,              &
                               krylov_truncated_on_the_trust_ball,        &
                               gm_reorthogonalize,                        &
                               gm_measure_orthogonality,                  &
                               gm_orthogonality_loss_cycle,               &
                               gm_reorthogonalization_loss_level,         &
                               the_point_where_the_step_leaves_the_ball,  &
                               name_of_unknown
      use certification,   only: certification_row_measure
      use ieee_arithmetic, only: ieee_is_finite, ieee_value,              &
                                 ieee_quiet_nan
      use iso_c_binding,   only: c_char, c_int, c_null_char

      implicit none
      ! The environment is what the solver's controls are read from, so a row
      ! that measures a control has to write one. setenv is the libc call
      ! get_environment_variable reads back.
      interface
         integer(c_int) function c_setenv(nm, vl, overwrite)              &
                                 bind(C, name='setenv')
         import :: c_char, c_int
         character(kind=c_char), dimension(*), intent(in) :: nm, vl
         integer(c_int), value :: overwrite
         end function c_setenv
      end interface
      ! LAPACK's banded factorization, used by the preconditioner rows to
      ! factor the manufactured band before and after its equilibration.
      interface
         subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
         integer, intent(in) :: m, n, kl, ku, ldab
         real*8, intent(inout) :: ab(ldab,*)
         integer, intent(out) :: ipiv(*), info
         end subroutine dgbtrf
      end interface
      integer :: nfail, kd_it, kd_outcome, kd_status
      real*8  :: kd_resid, kd_merit, kd_tau, kd_pred, kd_pref, kd_nan
      real*8  :: kd_snorm, kd_delta, kd_dmax
      ! The ray test's inputs and what it returned.
      real*8  :: kd_f2p, kd_f2m, kd_mslope, kd_fslope, kd_repro
      real*8  :: kd_ray_floor
      integer :: kd_verdict, kd_refuse, kd_env_status
      logical :: kd_model_ok
      logical :: kd_ascends
      real*8, allocatable :: kd_x(:), kd_A(:,:), kd_b(:)
      real*8, allocatable :: kd_r(:), kd_As(:)

      nfail = 0
      N     = 1
      call set_transported_species_rows(.false.)

      ! ---- the zero operator with a nonzero right-hand side ----
      ! One direction is sampled and its image is the zero vector, so the
      ! Krylov space of the operator is {0}: the reduced problem has rank
      ! zero and the system has no solution. Nothing about that is
      ! convergence.
      allocate(kd_A(1,1), kd_b(1))
      kd_A(1,1) = 0.0d0
      kd_b(1)   = 1.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call log_row('zero_operator_step_is_finite',                        &
                   ieee_is_finite(kd_x(1)), .true., nfail)
      call int_row('zero_operator_outcome_is_named_singular', kd_outcome, &
                   gm_reduced_system_singular, nfail)
      call rel_row('zero_operator_true_residual_is_nonzero', kd_resid,    &
                   1.0d0, 1.0d-12, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- the identity operator, the control ----
      allocate(kd_A(1,1), kd_b(1))
      kd_A(1,1) = 1.0d0
      kd_b(1)   = 1.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call rel_row('identity_operator_solution_is_one', kd_x(1), 1.0d0,   &
                   1.0d-14, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- a zero right-hand side is solved by the zero step ----
      ! No factorization and no normalization: zero is not a scale.
      allocate(kd_A(1,1), kd_b(1))
      kd_A(1,1) = 1.0d0
      kd_b(1)   = 0.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call int_row('zero_rhs_costs_no_product', kd_it, 0, nfail)
      call int_row('zero_rhs_outcome_is_tolerance_reached', kd_outcome,   &
                   gm_tolerance_reached, nfail)
      call rel_row('zero_rhs_step_is_zero', kd_x(1), 0.0d0, 0.0d0, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- a rank-one operator, right-hand side in and out of range ----
      allocate(kd_A(2,2), kd_b(2))
      kd_A    = 1.0d0
      kd_b(1) = 1.0d0
      kd_b(2) = 1.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call log_row('compatible_singular_system_has_finite_step',          &
                   ieee_is_finite(kd_x(1)) .and. ieee_is_finite(kd_x(2)), &
                   .true., nfail)
      call rel_row('compatible_singular_system_residual_is_zero',         &
                   kd_resid, 0.0d0, 1.0d-12, nfail)
      deallocate(kd_x)
      kd_b(1) = 1.0d0
      kd_b(2) = 0.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call log_row('incompatible_singular_system_step_is_finite',         &
                   ieee_is_finite(kd_x(1)) .and. ieee_is_finite(kd_x(2)), &
                   .true., nfail)
      ! The least-squares residual of this operator and this right-hand
      ! side: the image is t(1,1), (1-t)^2 + t^2 is minimal at t = 1/2, and
      ! the residual norm is then sqrt(1/2) with a right-hand side of norm
      ! one.
      call int_row('incompatible_singular_system_is_named', kd_outcome,   &
                   gm_reduced_system_singular, nfail)
      call rel_row('incompatible_singular_system_true_residual',          &
                   kd_resid, sqrt(0.5d0), 1.0d-12, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- an Arnoldi remainder that underflows ----
      ! The remainder is 1e-200, so the sum of its squares underflows to
      ! zero and its norm reads as an exact breakdown while the residual of
      ! the step it certifies is 1e-200 and not zero. What the cycle reports
      ! has to be the residual it has.
      allocate(kd_A(2,2), kd_b(2))
      kd_A      = 0.0d0
      kd_A(1,1) = 1.0d0
      kd_A(2,2) = 1.0d0
      kd_A(2,1) = 1.0d-200
      kd_b(1)   = 1.0d0
      kd_b(2)   = 0.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call int_row('near_breakdown_outcome_is_tolerance_reached',         &
                   kd_outcome, gm_tolerance_reached, nfail)
      call rel_row('near_breakdown_true_residual_is_verified', kd_resid,  &
                   1.0d-200, 1.0d-6, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- an action that is not a number ----
      ! One product only, so that the cycle cannot walk into a basis vector
      ! the Arnoldi never wrote.
      kd_nan = ieee_value(1.0d0, ieee_quiet_nan)
      allocate(kd_A(1,1), kd_b(1))
      kd_A(1,1) = kd_nan
      kd_b(1)   = 1.0d0
      call run_krylov(kd_A, kd_b, 1, 1.0d-1, lin_test_dense, kd_x,        &
                      kd_it, kd_outcome, kd_resid,      &
                      kd_snorm)
      call log_row('unusable_jacobian_action_step_is_finite',             &
                   ieee_is_finite(kd_x(1)), .true., nfail)
      call int_row('unusable_jacobian_action_is_named', kd_outcome,       &
                   gm_jacobian_action_unusable, nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- a preconditioner that reports a failure ----
      ! The vector comes back untouched with a nonzero status: the direction
      ! the cycle would sample is not the preconditioned one, so there is
      ! nothing to build a subspace on.
      allocate(kd_A(1,1), kd_b(1))
      kd_A(1,1) = 1.0d0
      kd_b(1)   = 1.0d0
      call run_krylov(kd_A, kd_b, 4, 1.0d-1,                              &
                      lin_test_dense_preconditioner_fails, kd_x,          &
                      kd_it, kd_outcome, kd_resid, kd_snorm)
      call int_row('preconditioner_failure_is_named', kd_outcome,         &
                   gm_preconditioner_failed, nfail)
      call log_row('preconditioner_failure_step_is_finite',               &
                   ieee_is_finite(kd_x(1)), .true., nfail)
      deallocate(kd_A, kd_b, kd_x)

      ! ---- the radius the trust region starts from ----
      ! A cycle that sampled no direction returns the zero step, and a
      ! hundredth of zero is not a small radius: it is the absence of one.
      ! The projected gradient is available in exactly that case, and the
      ! radius must be measured on it.
      N        = 1
      nvar_jac = 3
      allocate(kd_r(3), kd_As(3))
      kd_r  = 0.0d0
      kd_As = [ 0.0d0, -1.0d2, 0.0d0 ]
      call initial_trust_region_radius(gm_no_direction_sampled, kd_r,     &
                                       kd_As, kd_delta, kd_dmax,          &
                                       kd_status)
      call int_row('zero_first_direction_then_a_direction_is_taken',      &
                   kd_status, tr_radius_from_the_gradient, nfail)
      call rel_row('zero_first_direction_radius_is_positive', kd_delta,   &
                   1.0d0, 1.0d-14, nfail)
      call log_row('zero_first_direction_ceiling_is_above_the_radius',    &
                   kd_dmax .gt. kd_delta, .true., nfail)
      ! With a Krylov step the radius is measured on it.
      kd_r  = [ 0.0d0, 0.0d0, 2.0d2 ]
      kd_As = [ 0.0d0, -1.0d2, 0.0d0 ]
      call initial_trust_region_radius(gm_tolerance_reached, kd_r, kd_As, &
                                       kd_delta, kd_dmax, kd_status)
      call int_row('krylov_step_sets_the_first_radius', kd_status,        &
                   tr_radius_from_the_krylov_step, nfail)
      call rel_row('krylov_step_radius_is_a_hundredth_of_it', kd_delta,   &
                   2.0d0, 1.0d-14, nfail)
      ! Neither direction: there is nothing to size, and saying so is the
      ! answer.
      kd_r  = 0.0d0
      kd_As = 0.0d0
      call initial_trust_region_radius(gm_no_direction_sampled, kd_r,     &
                                       kd_As, kd_delta, kd_dmax,          &
                                       kd_status)
      call int_row('no_feasible_descent_is_named', kd_status,             &
                   tr_radius_no_feasible_descent, nfail)
      call rel_row('no_feasible_descent_radius_is_zero', kd_delta, 0.0d0, &
                   0.0d0, nfail)
      deallocate(kd_r, kd_As)

      ! ---- the dogleg inside a ball of radius zero ----
      call dogleg_row('radius_zero_dogleg_step_is_zero', nfail)

      ! ---- the step length along the banded gradient ----
      ! Residual one, operator one, band ten: the band gradient is ten, its
      ! image is ten, so the squared norms are one hundred each and the
      ! residual times the image is ten.
      call cauchy_length_along_the_banded_gradient(10.0d0, 100.0d0,       &
                                                   kd_tau, kd_ascends)
      kd_merit = 0.5d0*(1.0d0 - kd_tau*10.0d0)**2
      call rel_row('approximate_gradient_length_decreases_the_model',     &
                   kd_merit, 0.0d0, 1.0d-12, nfail)
      ! The exact gradient of the same model: band one, so the gradient is
      ! one and the length is the squared gradient norm over the squared
      ! norm of its image, unchanged.
      call cauchy_length_along_the_banded_gradient(1.0d0, 1.0d0, kd_tau,  &
                                                   kd_ascends)
      call rel_row('exact_gradient_length_is_unchanged', kd_tau, 1.0d0,   &
                   0.0d0, nfail)
      ! A negative slope: the direction ascends the model and no leg is
      ! taken along it.
      call cauchy_length_along_the_banded_gradient(-1.0d0, 1.0d0, kd_tau, &
                                                   kd_ascends)
      call log_row('ascending_direction_is_not_taken', kd_ascends,        &
                   .true., nfail)
      call rel_row('ascending_direction_length_is_zero', kd_tau, 0.0d0,   &
                   0.0d0, nfail)

      ! ---- the predicted drop is the drop of the explicit model ----
      N        = 1
      nvar_jac = 3
      allocate(kd_r(3), kd_As(3))
      kd_r  = [ 1.0d0, 2.0d0, 3.0d0 ]
      kd_As = [ 1.0d-1, -2.0d-1, 5.0d-2 ]
      kd_pred = predicted_model_decrease(kd_r, kd_As)
      kd_pref = 0.5d0*sum(kd_r*kd_r)                                      &
              - 0.5d0*sum((kd_r+kd_As)*(kd_r+kd_As))
      call rel_row('predicted_decrease_matches_the_explicit_model',       &
                   kd_pred, kd_pref, 1.0d-12, nfail)
      deallocate(kd_r, kd_As)

      ! ---- the ray test's floor is the measured reproducibility ----
      ! The refusing iteration of the atomic element solve, with the merit's
      ! reproducibility measured rather than assumed. The difference the
      ! quotient is built from, 2.97e-12, stands a factor 11 above the
      ! assumed floor eq_sweep_reltol*merit = 2.661e-13 and BELOW the
      ! reproducibility of the merit, so the quotient is a quotient of
      ! noise: -1.483e-9 against the model's -9.132e-9, a factor 6.2 apart
      ! on a step whose own reductions agree. The reduction ratio is set
      ! well outside the sound band here, so what this row measures is the
      ! floor alone.
      kd_f2p    =  2.661d-5
      kd_f2m    =  2.661d-5 + 2.97d-12
      kd_mslope = -9.132d-9
      kd_fslope = -2.97d-12/(2.0d0*1.0d-3)
      kd_repro  =  1.0d-12
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .true., 0.5d0, kd_repro, kd_ray_floor,       &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('ray_floor_is_measured_not_assumed', kd_refuse,        &
                   tr_ray_below_the_noise_floor, nfail)
      call rel_row('ray_floor_is_the_reproducibility_times_its_factor',   &
                   kd_ray_floor, ray_floor_factor*kd_repro, 1.0d-12,      &
                   nfail)
      ! With no reproducibility to read, the floor is the assumed one and
      ! the same inputs are refused as a slope disagreement: this is the
      ! verdict the measured floor replaces.
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .true., 0.5d0, 0.0d0, kd_ray_floor,          &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('assumed_floor_alone_calls_it_a_slope_disagreement',   &
                   kd_refuse, tr_ray_slopes_disagree_size, nfail)

      ! ---- a sound ratio is not overruled by the slope test ----
      ! The same disagreement of a factor 6 between the model's slope and
      ! the quotient, this time over a difference of 2.0e-11 that stands
      ! above the measured floor as well as the assumed one, and the
      ! reduction ratio the refused trials of that solve carried. The
      ! verdict is still formed and still says the slopes differ; it is not
      ! acted on, because the ratio is evidence about the same step taken
      ! without a finite difference.
      kd_f2p    =  2.661d-5
      kd_f2m    =  2.661d-5 + 2.0d-11
      kd_mslope = -6.0d-8
      kd_fslope = -2.0d-11/(2.0d0*1.0d-3)
      kd_repro  =  1.0d-12
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .true., 1.0003d0, kd_repro, kd_ray_floor,    &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('sound_ratio_is_not_overruled', kd_refuse,             &
                   tr_step_taken, nfail)
      call log_row('sound_ratio_leaves_the_model_usable', kd_model_ok,    &
                   .true., nfail)
      call int_row('sound_ratio_keeps_the_ray_verdict_reportable',        &
                   kd_verdict, tr_ray_slopes_disagree_size, nfail)

      ! ---- an over-predicting step is still refused ----
      ! The same disagreement with a ratio far from unity: the band is
      ! narrow enough that a model which promised three times what it
      ! delivered is still refused on the ray test, before the acceptance
      ! threshold is ever reached.
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .true., 0.3d0, kd_repro, kd_ray_floor,       &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('over_prediction_still_refuses', kd_refuse,            &
                   tr_ray_slopes_disagree_size, nfail)
      ! And a probe with no merit value is refused whatever the ratio: the
      ! test could not be taken, which is not the same as a noisy quotient.
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .false., 1.0d0, kd_repro, kd_ray_floor,      &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('inadmissible_probe_is_refused_at_any_ratio',          &
                   kd_refuse, tr_ray_probe_inadmissible, nfail)
      ! Slopes that agree are not refused, and then the ratio is the only
      ! thing left to decide the step.
      call ray_slope_verdict(kd_mslope, kd_mslope, kd_f2p, kd_f2m,        &
                             .true., 0.3d0, kd_repro, kd_ray_floor,       &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('agreeing_slopes_leave_the_ratio_to_decide',           &
                   kd_refuse, tr_step_taken, nfail)

      ! ---- a catastrophic over-prediction is named by the ratio ----
      ! The refusing iterate of the atomic element reload (N7b): the ratio
      ! is -533 because the actual reduction is NEGATIVE, and the two slopes
      ! disagree in sign as well. Both statements refuse the trial; the one
      ! the ledger has to record is the ratio's, which is about the whole
      ! step, and not the quotient's, which is about a thousandth of it.
      kd_f2p    = 2.661d-5
      kd_f2m    = 2.661d-5 + 3.391d-10
      kd_mslope = -5.541d-7
      kd_fslope =  1.696d-7
      kd_repro  = 0.0d0
      call ray_slope_verdict(kd_mslope, kd_fslope, kd_f2p, kd_f2m,        &
                             .true., -5.33d2, kd_repro, kd_ray_floor,     &
                             kd_verdict, kd_model_ok, kd_refuse)
      call int_row('over_prediction_is_named_by_the_ratio', kd_refuse,    &
                   tr_ratio_below_acceptance, nfail)
      ! The verdict about the quotient is still formed and still available;
      ! only the name the step control records changes.
      call int_row('over_prediction_keeps_the_ray_verdict_reportable',    &
                   kd_verdict, tr_ray_slopes_disagree_sign, nfail)
      ! And the model is still not believed: the region is cut on the same
      ! rule it was cut on before.
      call log_row('over_prediction_leaves_the_model_unusable',           &
                   kd_model_ok, .false., nfail)

      ! ---- one functional: the merit and the gate read one scale ----
      call merit_and_gate_rows(nfail)

      ! ---- the approximate-gradient leg and the slope of the step ----
      call cauchy_leg_rows(nfail)
      call cauchy_leg_image_share_rows(nfail)

      ! ---- the absolute cap on the trust radius ----
      call radius_cap_rows(nfail)

      ! ---- a manufactured operator with prescribed scales (N8a) ----
      call column_equilibration_rows(nfail)

      ! ---- the predicted decrease against the step it is taken for ----
      call predicted_decrease_scaling_rows(nfail)
      call predicted_decrease_ladder_rows(nfail)
      call step_out_of_the_box_rows(nfail)
      call jacobian_probe_length_rows(nfail)

      ! ---- the restart of the trust-region state (N21) ----
      call restart_trigger_rows(nfail)
      call restart_reset_rows(nfail)
      call krylov_restart_cycle_rows(nfail)

      ! ---- which rule admits the approximate-gradient leg (N24) ----
      call cauchy_leg_admission_rule_rows(nfail)

      ! ---- the basis, the ball and the name of a row (N25) ----
      call reorthogonalization_rows(nfail)
      call krylov_on_the_ball_rows(nfail)
      call row_name_rows(nfail)

      ! ---- the outer iteration cap ----
      kd_env_status = c_setenv('EXHALE_JFNK_MAXIT'//c_null_char,          &
                               '77'//c_null_char, 1)
      call int_row('jfnk_maxit_environment_is_writable', kd_env_status, 0, &
                   nfail)
      call read_species_unknown_space_controls
      call int_row('jfnk_maxit_hook', jfnk_outer_iteration_cap(500), 77,  &
                   nfail)
      kd_env_status = c_setenv('EXHALE_JFNK_MAXIT'//c_null_char,          &
                               ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call int_row('jfnk_maxit_default_is_the_caller',                    &
                   jfnk_outer_iteration_cap(500), 500, nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'krylov_and_dogleg: ', nfail,                &
              ' row(s) failed'
         stop 1
      endif
      write(*,'(A)') 'krylov_and_dogleg: all rows passed'

      contains

      ! ------------------------------------------------------!

      subroutine column_equilibration_rows(nf)
      ! A MANUFACTURED OPERATOR WHOSE SCALED CONDITION IS KNOWN, and what
      ! the column scaling does to it (item N8a, deliverables 2 and 5).
      !
      ! The operator is diagonal, so its scaled form is diagonal too and
      ! every quantity below is arithmetic rather than an estimate: with
      ! column scales D and row scales Drow fixed, A = Drow^-1 J D has
      ! diagonal entries J(j,j) D(j)/Drow(j) and its condition number in any
      ! norm is the ratio of the largest to the smallest of them. Four
      ! unknowns are given entries spanning eight decades, which is the kind
      ! of spread the certification row scaling produces: MEASURED on the
      ! atomic element reload, the scaled columns of the real operator span
      ! 2.6e-1 to 5.4e+6 at the first iterate.
      !
      ! WHAT IS RED BEFORE. With the state scales as the column scaling the
      ! spread survives into the operator, and it is exactly the spread of
      ! the entries: 1e8 here. That is the state the entry text hands the
      ! banded factorization and the Krylov space.
      !
      ! WHAT IS GREEN AFTER. equilibrate_the_scaled_columns brings every
      ! column into [1/sqrt(2), sqrt(2)], so the spread is at most 2, and it
      ! does so WITHOUT CHANGING THE MAP: the factor it multiplies a column
      ! by is the factor it multiplies that column's scale by, and it is a
      ! power of two, so A_after (D_after^-1 D_before) = A_before exactly
      ! and the same unknown-space step is reachable. An empty column keeps
      ! its scale, there being nothing in it to equilibrate.
      integer, intent(inout) :: nf
      integer, parameter :: nq = 5
      real*8  :: kq_ab(1,nq), kq_ab0(1,nq), kq_D(nq), kq_D0(nq)
      real*8  :: kq_s(nq), kq_before_min, kq_before_max
      real*8  :: kq_after_min, kq_after_max, kq_worst, kq_e
      integer :: kq_j, kq_nsave, kq_nvarsave, kq_klsave, kq_kusave
      kq_nsave    = N
      kq_nvarsave = nvar_jac
      kq_klsave   = kl_jac
      kq_kusave   = ku_jac
      N        = 1
      nvar_jac = nq
      kl_jac   = 0
      ku_jac   = 0
      ! The diagonal of the SCALED operator, eight decades of spread, plus
      ! one column that is empty.
      kq_ab(1,1) = 1.0d-4
      kq_ab(1,2) = 1.0d0
      kq_ab(1,3) = 1.0d4
      kq_ab(1,4) = 3.0d0
      kq_ab(1,5) = 0.0d0
      kq_ab0     = kq_ab
      kq_D       = (/ 2.0d0, 5.0d0, 7.0d0, 1.0d0, 9.0d0 /)
      kq_D0      = kq_D
      kq_before_min = huge(1.0d0);  kq_before_max = 0.0d0
      do kq_j = 1, 4
         kq_before_min = min(kq_before_min, abs(kq_ab0(1,kq_j)))
         kq_before_max = max(kq_before_max, abs(kq_ab0(1,kq_j)))
      enddo
      call rel_row('manufactured_columns_span_eight_decades',             &
                   kq_before_max/kq_before_min, 1.0d8, 1.0d-12, nf)
      call equilibrate_the_scaled_columns(kq_ab, kq_D)
      kq_after_min = huge(1.0d0);  kq_after_max = 0.0d0
      do kq_j = 1, 4
         kq_after_min = min(kq_after_min, abs(kq_ab(1,kq_j)))
         kq_after_max = max(kq_after_max, abs(kq_ab(1,kq_j)))
      enddo
      call log_row('equilibrated_columns_are_of_unit_size',               &
                   (kq_after_max .le. sqrt(2.0d0) .and.                   &
                    kq_after_min .ge. 1.0d0/sqrt(2.0d0)), .true., nf)
      call log_row('equilibrated_spread_is_at_most_two',                  &
                   (kq_after_max/kq_after_min .le. 2.0d0), .true., nf)
      ! The map is unchanged: undoing the column scale returns the entry.
      kq_worst = 0.0d0
      do kq_j = 1, 4
         kq_worst = max(kq_worst,                                         &
              abs(kq_ab(1,kq_j)*kq_D0(kq_j)/kq_D(kq_j) - kq_ab0(1,kq_j)))
      enddo
      call rel_row('equilibration_leaves_the_map_exactly_unchanged',      &
                   kq_worst, 0.0d0, 0.0d0, nf)
      ! And so is the step it can take: the same unknown-space displacement
      ! is reached by s scaled the other way, to the last bit.
      kq_worst = 0.0d0
      do kq_j = 1, 4
         kq_s(kq_j) = 0.125d0*dble(kq_j)
         kq_worst = max(kq_worst,                                         &
              abs(kq_D(kq_j)*(kq_s(kq_j)*kq_D0(kq_j)/kq_D(kq_j))          &
                  - kq_D0(kq_j)*kq_s(kq_j)))
      enddo
      call rel_row('equilibration_reaches_the_same_unknown_space_step',   &
                   kq_worst, 0.0d0, 0.0d0, nf)
      ! Every factor is a power of two, which is why the two rows above are
      ! exact and not merely close.
      kq_worst = 0.0d0
      do kq_j = 1, 4
         kq_e = kq_D(kq_j)/kq_D0(kq_j)
         kq_worst = max(kq_worst,                                         &
                        abs(kq_e - 2.0d0**nint(log(kq_e)/log(2.0d0))))
      enddo
      call rel_row('equilibration_factor_is_a_power_of_two', kq_worst,    &
                   0.0d0, 0.0d0, nf)
      ! An empty column has nothing to equilibrate and keeps its scale.
      call rel_row('empty_column_keeps_its_scale', kq_D(5), kq_D0(5),     &
                   0.0d0, nf)

      ! ---- the element column scale is tied to its reservoir ----
      ! A cell holding none of the element takes the reservoir floor, so its
      ! probe step is sqrt(eps)*1e-4*X_reservoir and the column it produces
      ! is resolvable; a cell holding more than the floor keeps its own
      ! value. The third row is the arithmetic the floor was chosen by: the
      ! probe step has to stand above the round-off of a row of size
      ! X_reservoir, which is eps*X_reservoir.
      call rel_row('empty_element_cell_takes_the_reservoir_floor',        &
                   element_unknown_column_scale(0.0d0, 1.0d-3),          &
                   1.0d-3*element_column_scale_floor, 1.0d-12, nf)
      call rel_row('occupied_element_cell_keeps_its_own_value',           &
                   element_unknown_column_scale(1.0d-2, 1.0d-3),         &
                   1.0d-2, 1.0d-12, nf)
      call log_row('element_probe_step_stands_above_the_row_round_off',   &
                   (sqrt(epsilon(1.0d0))                                  &
                      *element_unknown_column_scale(0.0d0, 1.0d0)         &
                    .gt. 1.0d3*epsilon(1.0d0)), .true., nf)
      ! ---- the band that is factorized, before and after (N8a) ----
      call preconditioner_equilibration_rows(nf)
      N        = kq_nsave
      nvar_jac = kq_nvarsave
      kl_jac   = kq_klsave
      ku_jac   = kq_kusave
      end subroutine column_equilibration_rows

      ! ------------------------------------------------------!

      subroutine preconditioner_equilibration_rows(nf)
      ! WHAT THE EQUILIBRATION OF THE FACTORIZED BAND BUYS, on a band whose
      ! scale spread is stated here (item N8a, deliverable 2).
      !
      ! A tridiagonal system of four unknowns whose rows are scaled by
      ! successive powers of a million, which is the shape of the real one:
      ! MEASURED on the atomic element reload, the rows of the scaled
      ! operator span 8.5e-1 to 5.4e+6 at the first iterate and its smallest
      ! pivot stood at 3.7e-19 against a largest of 1.05e7 at iteration 11.
      !
      ! RED: factored as the scalings leave it, the elimination's pivot
      ! ratio is the scale spread itself, and a triangular solve with it is
      ! a division by round-off.
      ! GREEN: equilibrate_the_preconditioner brings every row and column of
      ! the SAME matrix to unit size, and the pivot ratio is of order one.
      ! The band the model is built on is not touched by any of this.
      integer, intent(inout) :: nf
      integer, parameter :: np = 4
      real*8  :: kp_ab(4,np), kp_work(4,np), kp_sc
      integer :: kp_ipiv(np), kp_info, kp_j, kp_k
      real*8  :: kp_umin, kp_umax, kp_ratio_before, kp_ratio_after
      N        = 1
      nvar_jac = np
      kl_jac   = 1
      ku_jac   = 1
      ! A tridiagonal [-1, 2, -1] whose ROW i is multiplied by 1e6^(i-1),
      ! which is how a row scaling fixed by the gate reaches the band.
      kp_ab = 0.0d0
      do kp_j = 1, np
         do kp_k = max(1, kp_j-ku_jac), min(np, kp_j+kl_jac)
            kp_sc = 1.0d6**(kp_k-1)
            if (kp_k .eq. kp_j) then
               kp_ab(kl_jac+ku_jac+1 + kp_k - kp_j, kp_j) = 2.0d0*kp_sc
            else
               kp_ab(kl_jac+ku_jac+1 + kp_k - kp_j, kp_j) = -1.0d0*kp_sc
            endif
         enddo
      enddo
      kp_work = kp_ab
      call dgbtrf(np, np, kl_jac, ku_jac, kp_work, 4, kp_ipiv, kp_info)
      call int_row('manufactured_band_factors', kp_info, 0, nf)
      kp_umin = huge(1.0d0);  kp_umax = 0.0d0
      do kp_k = 1, np
         kp_umin = min(kp_umin, abs(kp_work(kl_jac+ku_jac+1, kp_k)))
         kp_umax = max(kp_umax, abs(kp_work(kl_jac+ku_jac+1, kp_k)))
      enddo
      kp_ratio_before = kp_umin/kp_umax
      call log_row('unequilibrated_pivot_ratio_is_the_scale_spread',      &
                   (kp_ratio_before .lt. 1.0d-15), .true., nf)
      kp_work = kp_ab
      call equilibrate_the_preconditioner(kp_work)
      call log_row('equilibration_is_recorded', preconditioner_equilibrated,&
                   .true., nf)
      call dgbtrf(np, np, kl_jac, ku_jac, kp_work, 4, kp_ipiv, kp_info)
      call int_row('equilibrated_band_factors', kp_info, 0, nf)
      kp_umin = huge(1.0d0);  kp_umax = 0.0d0
      do kp_k = 1, np
         kp_umin = min(kp_umin, abs(kp_work(kl_jac+ku_jac+1, kp_k)))
         kp_umax = max(kp_umax, abs(kp_work(kl_jac+ku_jac+1, kp_k)))
      enddo
      kp_ratio_after = kp_umin/kp_umax
      call log_row('equilibrated_pivot_ratio_is_of_order_one',            &
                   (kp_ratio_after .gt. 1.0d-3), .true., nf)
      ! MEASURED here: 5.0e-18 before, 5.4e-1 after, seventeen decades.
      call log_row('equilibration_improves_the_pivot_ratio',              &
                   (kp_ratio_after .gt. 1.0d12*kp_ratio_before), .true.,  &
                   nf)
      ! Both rescalings are powers of two, so nothing was rounded into the
      ! band on the way.
      kp_sc = 0.0d0
      do kp_k = 1, np
         kp_sc = max(kp_sc, abs(precon_row_scale(kp_k)                    &
                     - 2.0d0**nint(log(precon_row_scale(kp_k))            &
                                   /log(2.0d0))))
         kp_sc = max(kp_sc, abs(precon_col_scale(kp_k)                    &
                     - 2.0d0**nint(log(precon_col_scale(kp_k))            &
                                   /log(2.0d0))))
      enddo
      call rel_row('preconditioner_rescalings_are_powers_of_two', kp_sc,  &
                   0.0d0, 0.0d0, nf)
      preconditioner_equilibrated = .false.
      end subroutine preconditioner_equilibration_rows

      ! ------------------------------------------------------!

      subroutine run_krylov(kd_Amat, kd_rhs, kd_m, kd_rtol, kd_kind,      &
                            kd_xout, kd_gm_it, kd_gm_out, kd_gm_rel,      &
                            kd_gm_snorm)
      ! ONE CYCLE OF THE PRODUCTION KRYLOV SOLVE ON A STATED OPERATOR.
      ! The column and row scales are one and the pseudo-transient shift is
      ! zero, so the scaled operator the cycle builds is the stated matrix
      ! itself, and the preconditioner the test installs is the identity:
      ! what the rows read is the Arnoldi recursion, the rotations and the
      ! triangular solution with nothing in front of them.
      real*8, dimension(:,:), intent(in)  :: kd_Amat
      real*8, dimension(:),   intent(in)  :: kd_rhs
      integer,                intent(in)  :: kd_m, kd_kind
      real*8,                 intent(in)  :: kd_rtol
      real*8, allocatable,    intent(out) :: kd_xout(:)
      integer,                intent(out) :: kd_gm_it, kd_gm_out
      real*8,                 intent(out) :: kd_gm_rel, kd_gm_snorm
      real*8,  allocatable :: kd_Y(:), kd_F0(:), kd_D(:), kd_Drow(:)
      real*8,  allocatable :: kd_abf(:,:), kd_bb(:)
      integer, allocatable :: kd_ipiv(:)
      real*8  :: kd_fsp(1-Ng:1+Ng,n_species)
      integer :: kd_neq
      kd_neq   = size(kd_rhs)
      N        = 1
      nvar_jac = kd_neq
      kl_jac   = 0
      ku_jac   = 0
      allocate(kd_Y(kd_neq), kd_F0(kd_neq), kd_D(kd_neq))
      allocate(kd_Drow(kd_neq), kd_xout(kd_neq), kd_ipiv(kd_neq))
      allocate(kd_abf(1,kd_neq), kd_bb(kd_neq))
      kd_Y    = 0.0d0
      kd_F0   = 0.0d0
      kd_D    = 1.0d0
      kd_Drow = 1.0d0
      kd_ipiv = 0
      kd_abf  = 0.0d0
      kd_fsp  = 0.0d0
      kd_bb   = kd_rhs
      call linear_operator_set_for_test(kd_Amat, kd_kind)
      call pgmres(kd_Y, kd_F0, kd_fsp, kd_D, kd_Drow, kd_abf, kd_ipiv,    &
                  0.0d0, kd_bb, kd_xout, kd_m, kd_rtol, kd_gm_it,         &
                  kd_gm_out, kd_gm_rel, kd_gm_snorm)
      call linear_operator_clear_for_test
      deallocate(kd_Y, kd_F0, kd_D, kd_Drow, kd_ipiv, kd_abf, kd_bb)
      end subroutine run_krylov

      ! ------------------------------------------------------!

      subroutine dogleg_row(name, nf)
      ! A ball of radius zero contains one point. The geometry is handed two
      ! nonzero legs and must return that point.
      character(len=*), intent(in)    :: name
      integer,          intent(inout) :: nf
      real*8  :: kd_sU(1), kd_sN(1), kd_aU, kd_aN, kd_snorm
      logical :: kd_on_boundary
      N        = 1
      nvar_jac = 1
      kd_sU(1) = -1.0d0
      kd_sN(1) = -2.0d0
      call dogleg_step(kd_sU, kd_sN, 0.0d0, kd_aU, kd_aN, kd_on_boundary)
      kd_snorm = abs(kd_aU*kd_sU(1) + kd_aN*kd_sN(1))
      call rel_row(name, kd_snorm, 0.0d0, 0.0d0, nf)
      end subroutine dogleg_row

      ! ------------------------------------------------------!

      subroutine merit_and_gate_rows(nf)
      ! THE MERIT AND THE GATE ON ONE SCALING, on a three-row system with an
      ! element row (decision 20 a).
      !
      ! The gate is certification_row_measure, the maximum over cells of
      ! |res|/scale of one row. The merit is the 2-norm of every row divided
      ! by merit_row_scale_from_certification of the same scale, on the
      ! residual as the solver holds it -- in code units, so the element row
      ! carries a unit conversion the hydrodynamic rows do not. If the two
      ! read one scale, the row the merit is largest on IS the row the gate
      ! reports, and an element row that is nonzero in one cell holds
      ! exactly its own certification measure squared of the squared merit.
      !
      ! The numbers are the shape the atomic element reload has: two
      ! hydrodynamic rows a few thousandths out of balance and an element
      ! row an order of magnitude worse, which the merit could not see at
      ! all while it divided the element row by a floored transport scale
      ! (N7 report, item 7).
      integer, intent(inout) :: nf
      integer, parameter :: nc = 4
      real*8  :: kd_res(3,nc), kd_sc(3,nc), kd_unit(3), kd_q(3)
      real*8  :: kd_dr, kd_scaled, kd_worst, kd_total, kd_el
      integer :: kd_k, kd_j, kd_jw, kd_kgate, kd_kmerit
      logical :: kd_fin
      kd_res(1,:) = (/ 1.0d-6, 2.0d-6, 3.0d-6, 4.0d-6 /)
      kd_sc(1,:)  = 1.0d-3
      kd_res(2,:) = (/ 1.0d-5, 1.0d-5, 5.0d-5, 1.0d-5 /)
      kd_sc(2,:)  = 1.0d-2
      kd_res(3,:) = (/ 0.0d0,  0.0d0,  0.0d0,  9.0d-8 /)
      kd_sc(3,:)  = 1.0d-6
      ! The factor that writes each row per code time; the element row's is
      ! not one, which is the case the equality has to survive.
      kd_unit     = (/ 1.0d0, 1.0d0, 2.5d-1 /)
      kd_kgate = 0;  kd_worst = -1.0d0
      do kd_k = 1, 3
         call certification_row_measure(nc, nc+1, kd_res(kd_k,:),         &
                                        kd_sc(kd_k,:), kd_q(kd_k),        &
                                        kd_jw, kd_fin)
         if (kd_q(kd_k) .gt. kd_worst) then
            kd_worst = kd_q(kd_k);  kd_kgate = kd_k
         endif
      enddo
      kd_kmerit = 0;  kd_worst = -1.0d0;  kd_total = 0.0d0;  kd_el = 0.0d0
      do kd_k = 1, 3
         do kd_j = 1, nc
            kd_dr = merit_row_scale_from_certification(kd_sc(kd_k,kd_j),  &
                                                       kd_unit(kd_k))
            kd_scaled = abs(kd_res(kd_k,kd_j)*kd_unit(kd_k))/kd_dr
            kd_total  = kd_total + kd_scaled*kd_scaled
            if (kd_k .eq. 3) kd_el = kd_el + kd_scaled*kd_scaled
            if (kd_scaled .gt. kd_worst) then
               kd_worst = kd_scaled;  kd_kmerit = kd_k
            endif
         enddo
      enddo
      call int_row('merit_and_gate_name_the_same_row', kd_kmerit,         &
                   kd_kgate, nf)
      ! The element row's share of the squared merit against its own
      ! certification measure squared: one nonzero cell, so the sum over
      ! cells is the maximum over them.
      call rel_row('element_row_enters_the_merit_on_its_certification'//  &
                   '_scale', kd_el/kd_total,                              &
                   kd_q(3)*kd_q(3)/kd_total, 1.0d-13, nf)
      ! And the largest scaled row of the merit is the gate's own number,
      ! not a different one on a different scaling.
      call rel_row('largest_scaled_row_is_the_gate_measure', kd_worst,    &
                   kd_q(kd_kgate), 1.0d-13, nf)
      end subroutine merit_and_gate_rows

      ! ------------------------------------------------------!

      subroutine cauchy_leg_rows(nf)
      ! THE LEG THAT IS NOT IN THE STEP IS NOT IN THE SLOPE OF IT.
      !
      ! The lengths are the refusing iterate of the atomic element reload
      ! (N7): ||sU|| = 6.62e-9 against a radius of 3.911e-3, so the
      ! approximate-gradient leg is 1.7e-6 of the region, while its image
      ! A sU is not small at all and carries most of r0 . A s. Below
      ! cauchy_leg_min_fraction the leg is dropped from the step, and the
      ! model slope handed to the ray test is then r0 . A s of the step the
      ! dogleg actually returns.
      integer, intent(inout) :: nf
      real*8  :: kd_sU(2), kd_sN(2), kd_AsU(2), kd_AsN(2), kd_r0(2)
      real*8  :: kd_s(2), kd_As(2), kd_zero(2)
      real*8  :: kd_delta, kd_aU, kd_aN, kd_nU
      logical :: kd_on_boundary
      N        = 1
      nvar_jac = 2
      kd_delta  = 3.911d-3
      kd_sU     = (/ 6.62d-9, 0.0d0 /)
      kd_AsU    = (/ -5.0d-5, 0.0d0 /)
      kd_sN     = (/ 0.0d0, 1.0d-2 /)
      kd_AsN    = (/ 0.0d0, -1.0d-3 /)
      kd_r0     = (/ 7.29d-3, 7.29d-3 /)
      kd_zero   = 0.0d0
      kd_nU     = sqrt(sum(kd_sU*kd_sU))
      call log_row('cauchy_leg_is_short_beside_the_radius',               &
                   kd_nU .lt. cauchy_leg_min_fraction*kd_delta, .true., nf)
      ! The step with the leg dropped, from the production geometry.
      call dogleg_step(kd_zero, kd_sN, kd_delta, kd_aU, kd_aN,            &
                       kd_on_boundary)
      kd_s  = kd_aU*kd_zero + kd_aN*kd_sN
      kd_As = kd_aU*kd_zero + kd_aN*kd_AsN
      call rel_row('cauchy_leg_slope_is_the_step_slope', sum(kd_r0*kd_As),&
                   kd_aN*sum(kd_r0*kd_AsN), 0.0d0, nf)
      ! The step is the Krylov leg cut to the ball, so its length is the
      ! radius and nothing of the dropped leg is left in it.
      call rel_row('dropped_leg_leaves_the_step_on_the_ball',             &
                   sqrt(sum(kd_s*kd_s)), kd_delta, 1.0d-13, nf)
      ! WHY IT MATTERS THAT IT IS DROPPED FROM BOTH. Kept, the leg carries
      ! 7.8 percent of the model slope while it is 1.7e-6 of the step: its
      ! share of the slope stands four decades above its share of the step,
      ! which is what made the slope a statement about a direction the step
      ! did not take.
      call dogleg_step(kd_sU, kd_sN, kd_delta, kd_aU, kd_aN,              &
                       kd_on_boundary)
      call log_row('kept_leg_owns_more_of_the_slope_than_of_the_step',    &
                   abs(kd_aU*sum(kd_r0*kd_AsU)) .gt.                      &
                   1.0d4*(kd_nU/kd_delta)                                 &
                   *abs(kd_aN*sum(kd_r0*kd_AsN)), .true., nf)
      end subroutine cauchy_leg_rows

      ! ------------------------------------------------------!

      subroutine cauchy_leg_image_share_rows(nf)
      ! THE LEG IS ADMITTED BY ITS SHARE OF THE IMAGE, NOT BY ITS LENGTH.
      !
      ! The numbers are outer iteration 145 of the atomic element reload
      ! (N22 report, noticed item 1): r0 . A g = 8.262e13 against
      ! ||A g|| = 7.923e14, so tau = (r0.Ag)/||Ag||^2 = 1.316e-16 and the
      ! leg's image tau ||A g|| = 1.043e-1 is the whole of ||A s||, while
      ! its length is 6.9e-3 of the step and stands ABOVE
      ! cauchy_leg_min_fraction of the radius, so the length test keeps it.
      ! The model slope is then r0 . A s of a direction the step barely
      ! contains, and the predicted decrease 0.5 (r0.Ag)^2/||Ag||^2 that
      ! goes with it has no radius in it.
      integer, intent(inout) :: nf
      real*8  :: ci_sU(2), ci_sN(2), ci_AsU(2), ci_AsN(2), ci_r0(2)
      real*8  :: ci_s(2), ci_As(2), ci_zero(2)
      real*8  :: ci_delta, ci_aU, ci_aN, ci_tau, ci_step, ci_image
      real*8  :: ci_slope_kept, ci_slope_step
      logical :: ci_bnd, ci_ascends
      N        = 1
      nvar_jac = 2
      ci_delta = 1.7448d-7
      ci_zero  = 0.0d0
      ! The length along the banded gradient, from the production routine.
      call cauchy_length_along_the_banded_gradient(8.262d13,              &
                             (7.923d14)**2, ci_tau, ci_ascends)
      call rel_row('gradient_leg_length_is_the_production_one', ci_tau,   &
                   8.262d13/(7.923d14)**2, 1.0d-13, nf)
      ! sU = -tau g and A sU = -tau A g, with ||g|| chosen so that the leg
      ! is 6.8e-3 of the step and ||A g|| the measured 7.923e14.
      ci_sU  = (/ ci_tau*9.0900d6, 0.0d0 /)
      ci_AsU = (/ -ci_tau*7.923d14, 0.0d0 /)
      ci_sN  = (/ 0.0d0, 1.8d1 /)
      ci_AsN = (/ 0.0d0, -1.0d0 /)
      ci_r0  = (/ 8.262d13/7.923d14, 7.29d-3 /)
      ! THE DEFECT, STATED AS A ROW: the length test keeps this leg.
      call log_row('image_dominant_leg_is_kept_by_its_length',            &
                   sqrt(sum(ci_sU*ci_sU)) .ge.                            &
                   cauchy_leg_min_fraction*ci_delta, .true., nf)
      call dogleg_step(ci_sU, ci_sN, ci_delta, ci_aU, ci_aN, ci_bnd)
      call cauchy_leg_shares_of_step_and_image(ci_sU, ci_sN, ci_AsU,      &
                             ci_AsN, ci_aU, ci_aN, ci_step, ci_image)
      call rel_row('image_dominant_leg_step_share', ci_step, 6.856d-3,    &
                   2.0d-2, nf)
      call rel_row('image_dominant_leg_image_share', ci_image, 1.0d0,     &
                   1.0d-6, nf)
      call log_row('image_dominant_leg_is_dropped_on_its_image',          &
                   cauchy_leg_is_image_without_step(ci_step, ci_image),   &
                   .true., nf)
      ! The slope the model hands over with the leg kept, against the slope
      ! of the step that is taken once it is dropped.
      ci_slope_kept = sum(ci_r0*(ci_aU*ci_AsU + ci_aN*ci_AsN))
      call dogleg_step(ci_zero, ci_sN, ci_delta, ci_aU, ci_aN, ci_bnd)
      ci_s  = ci_aN*ci_sN
      ci_As = ci_aN*ci_AsN
      ci_slope_step = sum(ci_r0*ci_As)
      call rel_row('dropped_leg_step_is_the_krylov_leg_on_the_ball',      &
                   sqrt(sum(ci_s*ci_s)), ci_delta, 1.0d-13, nf)
      call rel_row('model_slope_is_the_slope_of_the_step_taken',          &
                   ci_slope_step, ci_aN*sum(ci_r0*ci_AsN), 0.0d0, nf)
      ! The slope of the step per unit of its length is the Krylov leg's
      ! own, the step being that leg cut to the ball.
      call rel_row('model_slope_per_unit_length_of_the_step',             &
                   ci_slope_step/sqrt(sum(ci_s*ci_s)),                     &
                   sum(ci_r0*ci_AsN)/sqrt(sum(ci_sN*ci_sN)), 1.0d-13, nf)
      ! Kept, the leg's slope is 1.087e-2 against the step's 7.066e-11, a
      ! factor 1.54e8; a million is the floor the row asserts.
      call log_row('kept_leg_slope_is_not_the_slope_of_the_step',         &
                   abs(ci_slope_kept) .gt. 1.0d6*abs(ci_slope_step),      &
                   .true., nf)
      ! THE CONTROL: the same geometry with the leg's image made small.
      ! Its image share falls below cauchy_leg_image_share_max and the leg
      ! is kept, however short it is, so nothing but an image-dominant leg
      ! is removed by this rule.
      ci_AsU = (/ -1.0d-9, 0.0d0 /)
      call dogleg_step(ci_sU, ci_sN, ci_delta, ci_aU, ci_aN, ci_bnd)
      call cauchy_leg_shares_of_step_and_image(ci_sU, ci_sN, ci_AsU,      &
                             ci_AsN, ci_aU, ci_aN, ci_step, ci_image)
      call log_row('control_leg_image_share_is_below_the_threshold',      &
                   ci_image .lt. cauchy_leg_image_share_max, .true., nf)
      call log_row('control_leg_is_kept', cauchy_leg_is_image_without_step&
                   (ci_step, ci_image), .false., nf)
      ci_s  = ci_aU*ci_sU  + ci_aN*ci_sN
      ci_As = ci_aU*ci_AsU + ci_aN*ci_AsN
      call rel_row('control_model_slope_is_the_slope_of_the_step_taken',  &
                   sum(ci_r0*ci_As),                                      &
                   ci_aU*sum(ci_r0*ci_AsU) + ci_aN*sum(ci_r0*ci_AsN),     &
                   0.0d0, nf)
      call rel_row('control_step_is_on_the_ball', sqrt(sum(ci_s*ci_s)),   &
                   ci_delta, 1.0d-6, nf)
      ! And here the kept leg's slope IS commensurate with the step's, which
      ! is why the rule leaves it alone: 1.7e-10 against 7.1e-11, a factor
      ! 2.5, where the image-dominant leg above stood eight decades away.
      call log_row('control_kept_leg_slope_is_commensurate_with_the_step',&
                   abs(sum(ci_r0*ci_As)) .lt.                             &
                   1.0d1*abs(ci_aN*sum(ci_r0*ci_AsN)), .true., nf)
      ! AND A LEG THAT IS THE STEP MAY BE ITS IMAGE. With the ball smaller
      ! than the gradient leg the dogleg is that leg cut to the ball and
      ! the Newton leg is out of the step entirely, so both shares are one
      ! and the rule keeps it: a leg that IS the step is allowed to be all
      ! of its image, which is the classical Cauchy step every convergence
      ! guarantee of a trust region rests on.
      ci_AsU = (/ -ci_tau*7.923d14, 0.0d0 /)
      ci_delta = 1.0d-10
      call dogleg_step(ci_sU, ci_sN, ci_delta, ci_aU, ci_aN, ci_bnd)
      call cauchy_leg_shares_of_step_and_image(ci_sU, ci_sN, ci_AsU,      &
                             ci_AsN, ci_aU, ci_aN, ci_step, ci_image)
      call log_row('a_leg_that_is_the_step_keeps_its_image',              &
                   ci_step .ge. cauchy_leg_step_share_max, .true., nf)
      call log_row('a_leg_that_is_the_step_is_kept',                      &
                   cauchy_leg_is_image_without_step(ci_step, ci_image),   &
                   .false., nf)
      end subroutine cauchy_leg_image_share_rows

      ! ------------------------------------------------------!

      subroutine cauchy_leg_admission_rule_rows(nf)
      ! WHICH RULE DECIDES THAT THE LEG IS NOT IN THE STEP, and what the
      ! two rules do differently. The length test asks whether the leg is
      ! small beside the RADIUS, the image test whether it is small beside
      ! the STEP while owning the model image; they agree only where the
      ! step stands on the ball, and where they part it is the length test
      ! that removes a leg the model has no complaint about. The switch
      ! cauchy_leg_admitted_by_image_alone leaves the image test alone in
      ! charge, and these rows state its default, its arming, and the two
      ! legs the two rules disagree about.
      integer, intent(inout) :: nf
      real*8  :: ar_sU(2), ar_sN(2), ar_AsU(2), ar_AsN(2), ar_r0(2)
      real*8  :: ar_s(2), ar_As(2)
      real*8  :: ar_delta, ar_aU, ar_aN, ar_tau, ar_step, ar_image
      logical :: ar_bnd, ar_ascends
      integer :: ar_status
      ! THE DEFAULT, read through the production reader with the variable
      ! absent from the environment.
      ar_status = c_setenv('EXHALE_CAUCHY_LEG_BY_IMAGE'//c_null_char,     &
                           ''//c_null_char, 1)
      call int_row('cauchy_leg_by_image_environment_is_writable',         &
                   ar_status, 0, nf)
      call read_species_unknown_space_controls
      call log_row('cauchy_leg_by_image_default_is_off',                  &
                   cauchy_leg_admitted_by_image_alone, .false., nf)
      ! AND ITS ARMING.
      ar_status = c_setenv('EXHALE_CAUCHY_LEG_BY_IMAGE'//c_null_char,     &
                           '1'//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('cauchy_leg_by_image_is_armed_by_the_environment',     &
                   cauchy_leg_admitted_by_image_alone, .true., nf)
      ar_status = c_setenv('EXHALE_CAUCHY_LEG_BY_IMAGE'//c_null_char,     &
                           ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('cauchy_leg_by_image_returns_to_its_default',          &
                   cauchy_leg_admitted_by_image_alone, .false., nf)
      ! THE LEG THE TWO RULES DISAGREE ABOUT: a hundred-thousandth of the
      ! radius in length, and little of the image. The length test drops
      ! it; the image test keeps it, because the model slope of the step
      ! is then still the slope of the step.
      N        = 1
      nvar_jac = 2
      ar_delta = 1.0d-7
      call log_row('a_short_leg_is_dropped_on_its_length',                &
                   cauchy_leg_is_dropped_on_its_length(1.0d-5*ar_delta,   &
                                       ar_delta, .false.), .true., nf)
      call log_row('the_same_leg_is_kept_when_the_image_alone_decides',   &
                   cauchy_leg_is_dropped_on_its_length(1.0d-5*ar_delta,   &
                                       ar_delta, .true.), .false., nf)
      ! A leg of no length at all is kept by the length rule under the
      ! switch, which is what "the length is not asked about" means.
      call log_row('a_leg_of_zero_length_is_not_dropped_on_its_length',   &
                   cauchy_leg_is_dropped_on_its_length(0.0d0, ar_delta,   &
                                       .true.), .false., nf)
      ! And a leg above the fraction is kept by the length rule either way.
      call log_row('a_long_leg_is_kept_on_its_length',                    &
                   cauchy_leg_is_dropped_on_its_length(1.0d-2*ar_delta,   &
                                       ar_delta, .false.), .false., nf)
      ! WHAT THE IMAGE TEST THEN SAYS ABOUT THAT LEG, and what the step
      ! and its model slope are with the leg in. The leg is 1e-5 of the
      ! radius with an image of 1e-9 against the Newton leg's own 1, so
      ! its image share is negligible and it stays in the step.
      ar_sU  = (/ 1.0d-5*ar_delta, 0.0d0 /)
      ar_AsU = (/ -1.0d-9, 0.0d0 /)
      ar_sN  = (/ 0.0d0, 1.8d1 /)
      ar_AsN = (/ 0.0d0, -1.0d0 /)
      ar_r0  = (/ 1.043d-1, 7.29d-3 /)
      call dogleg_step(ar_sU, ar_sN, ar_delta, ar_aU, ar_aN, ar_bnd)
      call cauchy_leg_shares_of_step_and_image(ar_sU, ar_sN, ar_AsU,      &
                             ar_AsN, ar_aU, ar_aN, ar_step, ar_image)
      call log_row('the_kept_leg_is_little_of_the_image',                 &
                   ar_image .lt. cauchy_leg_image_share_max, .true., nf)
      call log_row('the_kept_leg_is_not_dropped_on_its_image',            &
                   cauchy_leg_is_image_without_step(ar_step, ar_image),   &
                   .false., nf)
      ar_s  = ar_aU*ar_sU  + ar_aN*ar_sN
      ar_As = ar_aU*ar_AsU + ar_aN*ar_AsN
      call rel_row('the_model_slope_is_the_slope_of_the_step_with_the_leg',&
                   sum(ar_r0*ar_As),                                      &
                   ar_aU*sum(ar_r0*ar_AsU) + ar_aN*sum(ar_r0*ar_AsN),     &
                   0.0d0, nf)
      ! AND THE LEG THE TWO RULES AGREE ABOUT: the image-dominant one of
      ! outer iteration 145, which the image test removes at any length.
      call cauchy_length_along_the_banded_gradient(8.262d13,              &
                             (7.923d14)**2, ar_tau, ar_ascends)
      ar_delta = 1.7448d-7
      ar_sU  = (/ ar_tau*9.0900d6, 0.0d0 /)
      ar_AsU = (/ -ar_tau*7.923d14, 0.0d0 /)
      call dogleg_step(ar_sU, ar_sN, ar_delta, ar_aU, ar_aN, ar_bnd)
      call cauchy_leg_shares_of_step_and_image(ar_sU, ar_sN, ar_AsU,      &
                             ar_AsN, ar_aU, ar_aN, ar_step, ar_image)
      call log_row('the_image_dominant_leg_survives_the_length_rule',     &
                   cauchy_leg_is_dropped_on_its_length(                   &
                        sqrt(sum(ar_sU*ar_sU)), ar_delta, .false.),       &
                   .false., nf)
      call log_row('the_image_dominant_leg_is_dropped_by_the_image_rule', &
                   cauchy_leg_is_image_without_step(ar_step, ar_image),   &
                   .true., nf)
      end subroutine cauchy_leg_admission_rule_rows

      ! ------------------------------------------------------!

      subroutine radius_cap_rows(nf)
      ! THE CEILING HAS ONE ABSOLUTE BOUND, and it does not bind on an
      ! ordinary radius. In the scaled coordinates a step of length
      ! sqrt(neq) moves every unknown by one unit of its own scale, and the
      ! cap is a thousand of those; the atomic element reload ran the
      ! ceiling to 2.684e8 against a step of 1.018e-5 before it existed
      ! (N7b, noticed item 4).
      integer, intent(inout) :: nf
      real*8  :: kd_dZ(3), kd_g(3), kd_delta, kd_dmax
      integer :: kd_status, kd_i
      N        = 1
      nvar_jac = 3
      kd_g     = 0.0d0
      ! A Krylov step so long that a hundredfold of a hundredth of it stands
      ! above the cap.
      kd_dZ    = (/ 1.0d8, 0.0d0, 0.0d0 /)
      call initial_trust_region_radius(gm_tolerance_reached, kd_dZ, kd_g, &
                                       kd_delta, kd_dmax, kd_status)
      call rel_row('radius_ceiling_respects_the_absolute_cap', kd_dmax,   &
                   trust_region_radius_absolute_cap(), 1.0d-13, nf)
      ! And doubling it on every good step never passes the cap.
      do kd_i = 1, 200
         call grow_the_trust_region_ceiling(kd_dmax)
      enddo
      call rel_row('radius_cap_holds', kd_dmax,                           &
                   trust_region_radius_absolute_cap(), 1.0d-13, nf)
      ! An ordinary first radius is nowhere near it: the ceiling is a
      ! hundredfold of the radius, as it was.
      kd_dZ = (/ 1.0d0, 0.0d0, 0.0d0 /)
      call initial_trust_region_radius(gm_tolerance_reached, kd_dZ, kd_g, &
                                       kd_delta, kd_dmax, kd_status)
      call rel_row('radius_cap_does_not_bind_on_an_ordinary_radius',      &
                   kd_dmax, 1.0d0, 1.0d-13, nf)
      end subroutine radius_cap_rows

      ! ------------------------------------------------------!

      subroutine predicted_decrease_scaling_rows(nf)
      ! THE PREDICTED DECREASE FALLS WITH THE STEP IT IS TAKEN FOR.
      !
      ! The dogleg cut to a smaller ball returns a shorter step, and the
      ! model's predicted decrease of a step short enough for the quadratic
      ! term to be negligible is linear in the step:
      ! -(r0 . As) - 0.5 ||As||^2 with As proportional to the step. So
      ! quartering the radius quarters the predicted decrease to the
      ! accuracy of that quadratic term, and a ratio test that compares a
      ! predicted with an actual reduction is a comparison of two quantities
      ! of the step only while that holds.
      !
      ! MEASURED on the atomic element reload (item N21): it holds over the
      ! first 136 outer iterations, where a refusal and the accepted step of
      ! the next iteration stand at a step ratio of 0.250 and a predicted
      ! ratio of 0.16 to 0.28, and it stops holding from about iteration 140,
      ! where the step has fallen to 1e-8 and the predicted decrease stays at
      ! 7.099e-4 across a step ratio of 0.250. These rows fix the invariant
      ! on the production geometry so that the two places it could fail --
      ! the geometry, or the caller's bookkeeping of the step it took -- are
      ! separated.
      integer, intent(inout) :: nf
      real*8  :: kp_sU(2), kp_sN(2), kp_AsU(2), kp_AsN(2), kp_r0(2)
      real*8  :: kp_s(2), kp_As(2)
      real*8  :: kp_aU, kp_aN, kp_p1, kp_p4, kp_n1, kp_n4
      logical :: kp_bnd
      N        = 1
      nvar_jac = 2
      ! A Newton leg well outside the ball, so the dogleg is cut by the
      ! radius at both radii below and the step is proportional to it.
      kp_sU  = (/ -1.0d0,  0.0d0 /)
      kp_sN  = (/ -1.0d1, -1.0d1 /)
      kp_AsU = (/  1.0d0,  0.0d0 /)
      kp_AsN = (/  1.0d1,  1.0d1 /)
      kp_r0  = (/ -1.0d0, -1.0d0 /)
      call dogleg_step(kp_sU, kp_sN, 1.0d-4, kp_aU, kp_aN, kp_bnd)
      kp_s  = kp_aU*kp_sU  + kp_aN*kp_sN
      kp_As = kp_aU*kp_AsU + kp_aN*kp_AsN
      kp_n1 = sqrt(sum(kp_s*kp_s))
      kp_p1 = predicted_model_decrease(kp_r0, kp_As)
      call dogleg_step(kp_sU, kp_sN, 0.25d-4, kp_aU, kp_aN, kp_bnd)
      kp_s  = kp_aU*kp_sU  + kp_aN*kp_sN
      kp_As = kp_aU*kp_AsU + kp_aN*kp_AsN
      kp_n4 = sqrt(sum(kp_s*kp_s))
      kp_p4 = predicted_model_decrease(kp_r0, kp_As)
      call rel_row('quartered_radius_quarters_the_step', kp_n4/kp_n1,      &
                   0.25d0, 1.0d-12, nf)
      call rel_row('quartered_step_quarters_the_predicted_decrease',       &
                   kp_p4/kp_p1, 0.25d0, 1.0d-3, nf)
      call log_row('the_predicted_decrease_is_positive_at_both_radii',     &
                   (kp_p1 .gt. 0.0d0) .and. (kp_p4 .gt. 0.0d0), .true., nf)
      end subroutine predicted_decrease_scaling_rows

      ! ------------------------------------------------------!

      subroutine predicted_decrease_ladder_rows(nf)
      ! THE SAME INVARIANT DOWN TWELVE QUARTERINGS, from a radius of 1e-2 to
      ! 6e-10, which is the range of step lengths a stalling solve walks
      ! through (item N22).
      !
      ! The dogleg is cut by the radius at every rung, so the step is
      ! exactly proportional to the radius, its image is the same
      ! combination of the two leg images, and the predicted decrease is
      ! -(r0 . As) - 0.5 ||As||^2 with the quadratic term falling as the
      ! square. So the ratio of consecutive predicted decreases is 0.25 to
      ! the accuracy of that term, which at the top rung of this ladder is
      ! 0.375 of the radius, at EVERY rung and not only at the first:
      ! the geometry has no length scale of its own, and if the ratio walks
      ! away as the rungs get short the loss is in the arithmetic of the two
      ! terms.
      !
      ! What this fixes is the GEOMETRY. The step control's own bookkeeping
      ! of the step a trial took cannot be reached from here: it needs a
      ! state to evaluate the residual at. The two are separated on the
      ! atomic element reload by the model of each trial printed beside a
      ! fresh image of the same step (elem_diag_every_trial).
      integer, intent(inout) :: nf
      real*8  :: kl_sU(2), kl_sN(2), kl_AsU(2), kl_AsN(2), kl_r0(2)
      real*8  :: kl_s(2), kl_As(2)
      real*8  :: kl_aU, kl_aN, kl_delta, kl_pred, kl_prev, kl_norm
      real*8  :: kl_norm_prev, kl_worst_pred, kl_worst_step, kl_last
      integer :: kl_rung
      logical :: kl_bnd
      N        = 1
      nvar_jac = 2
      kl_sU  = (/ -1.0d0,  0.0d0 /)
      kl_sN  = (/ -1.0d1, -1.0d1 /)
      kl_AsU = (/  1.0d0,  0.0d0 /)
      kl_AsN = (/  1.0d1,  1.0d1 /)
      kl_r0  = (/ -1.0d0, -1.0d0 /)
      kl_delta      = 1.0d-2
      kl_prev       = 0.0d0
      kl_norm_prev  = 0.0d0
      kl_worst_pred = 0.0d0
      kl_worst_step = 0.0d0
      kl_last       = 0.0d0
      do kl_rung = 0, 12
         call dogleg_step(kl_sU, kl_sN, kl_delta, kl_aU, kl_aN, kl_bnd)
         kl_s    = kl_aU*kl_sU  + kl_aN*kl_sN
         kl_As   = kl_aU*kl_AsU + kl_aN*kl_AsN
         kl_norm = sqrt(sum(kl_s*kl_s))
         kl_pred = predicted_model_decrease(kl_r0, kl_As)
         if (kl_rung .gt. 0) then
            kl_worst_pred = max(kl_worst_pred,                            &
                            abs(kl_pred/kl_prev - 0.25d0)/0.25d0)
            kl_worst_step = max(kl_worst_step,                            &
                            abs(kl_norm/kl_norm_prev - 0.25d0)/0.25d0)
         endif
         kl_prev      = kl_pred
         kl_norm_prev = kl_norm
         kl_last      = kl_norm
         kl_delta     = 0.25d0*kl_delta
      enddo
      call rel_row('twelve_quarterings_reach_below_one_nanometer_of_step', &
                   kl_last, 5.96d-10, 1.0d-2, nf)
      call rel_row('every_quartering_quarters_the_step',                   &
                   kl_worst_step, 0.0d0, 1.0d-12, nf)
      call rel_row('every_quartering_quarters_the_predicted_decrease',     &
                   kl_worst_pred, 0.0d0, 1.0d-2, nf)
      call log_row('the_predicted_decrease_stays_positive_down_the_ladder',&
                   kl_prev .gt. 0.0d0, .true., nf)
      end subroutine predicted_decrease_ladder_rows

      ! ------------------------------------------------------!

      subroutine step_out_of_the_box_rows(nf)
      ! WHAT HAPPENS TO A STEP THAT CARRIES A SPECIES UNKNOWN OUT OF ITS
      ! BOX (item N22).
      !
      ! Two rules on one trial, so that they are distinguishable: the
      ! entry rule writes the unknown onto the nearest face of the box; the
      ! rule the step control uses holds the STEP there, which leaves the
      ! unknown at the value the iterate has. The second is the one a model
      ! of the step can describe: a component on a face has no admissible
      ! sample of the operator on either side, so every product the model
      ! is built from drops it, and a trial written onto the face carries a
      ! move no model of the step contains.
      !
      ! MEASURED on the atomic element reload, outer iterations 140 to 160
      ! of the entry text: exactly the eleven trials with one unknown on a
      ! face are refused and exactly the ten with none are accepted at a
      ! reduction ratio of 1.000, while the step in the unknowns themselves
      ! is the same to 0.4 percent across the pair.
      integer, intent(inout) :: nf
      real*8  :: kb_lo(8), kb_hi(8), kb_Y(8), kb_try(8), kb_face(8)
      integer :: kb_n, kb_nface
      N        = 2
      nvar_jac = 4
      ! Two cells, four unknowns each: the three hydrodynamic ones and one
      ! species unknown, whose box is [0, 1]. The species unknown of cell j
      ! is the fourth slot of that cell, which is where the box rules look.
      kb_lo = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
      kb_hi = (/ 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /)
      call species_box_set_for_test(1, kb_lo, kb_hi)
      ! The iterate holds both species unknowns inside their box; the step
      ! carries the first below its lower bound and leaves the second
      ! inside.
      kb_Y   = (/ 1.0d0, 2.0d0, 3.0d0,  2.0d-1,                          &
                  1.0d0, 2.0d0, 3.0d0,  5.0d-1 /)
      kb_try = (/ 1.5d0, 2.5d0, 3.5d0, -3.0d-1,                          &
                  1.5d0, 2.5d0, 3.5d0,  6.0d-1 /)
      kb_face = kb_try
      call species_unknowns_outside_their_bounds(kb_face, kb_nface, .true.)
      kb_n = 0
      call hold_the_step_where_it_leaves_the_species_box(kb_Y, kb_try,    &
                                                         kb_n)
      call int_row('one_unknown_leaves_the_box', kb_nface, 1, nf)
      call int_row('the_same_unknown_is_counted_by_the_hold', kb_n, 1, nf)
      call rel_row('the_entry_rule_writes_the_face', kb_face(4), 0.0d0,   &
                   1.0d-300, nf)
      call rel_row('the_step_is_held_at_the_iterate', kb_try(4),          &
                   2.0d-1, 1.0d-14, nf)
      call rel_row('an_unknown_inside_the_box_keeps_its_trial_value',     &
                   kb_try(8), 6.0d-1, 1.0d-14, nf)
      call rel_row('the_hydrodynamic_unknowns_keep_their_trial_values',   &
                   kb_try(1), 1.5d0, 1.0d-14, nf)
      ! With no species row registered neither rule touches anything: the
      ! three-unknown route is unreachable by construction.
      call species_box_clear_for_test
      kb_try = (/ 1.5d0, 2.5d0, 3.5d0, -3.0d-1,                          &
                  1.5d0, 2.5d0, 3.5d0,  6.0d-1 /)
      call hold_the_step_where_it_leaves_the_species_box(kb_Y, kb_try,    &
                                                         kb_n)
      call int_row('with_no_species_row_nothing_is_held', kb_n, 0, nf)
      call rel_row('with_no_species_row_the_trial_is_untouched',          &
                   kb_try(4), -3.0d-1, 1.0d-14, nf)
      end subroutine step_out_of_the_box_rows

      ! ------------------------------------------------------!

      subroutine jacobian_probe_length_rows(nf)
      ! THE ARC THE MATRIX-FREE ACTION SAMPLES THE RESIDUAL OVER (item
      ! N22).
      !
      ! The forward difference is taken at Y + eps0 v with
      ! eps0 = probe_length/||v||, so the DISPLACEMENT has the probe length
      ! whatever the direction and whatever its length. Two readings follow
      ! and both are asserted here: the quotient is homogeneous of degree
      ! one in v, so there is no direction short enough to drive the action
      ! into the rounding of the residual; and the arc is a property of the
      ! iterate, so a step shorter than it is a step the model of it does
      ! not reach.
      !
      ! The numbers are the atomic element reload's at outer iteration 151
      ! (MEASURED, item N22): ||Y|| = 11.3817 gives an arc of 1.845e-7, nine
      ! decades above the residual's reproducibility of 1e-15 relative, and
      ! the step there is ||D s|| = 5.0e-10, which is 2.7e-3 of the arc.
      integer, intent(inout) :: nf
      real*8  :: kj_Y(4), kj_arc, kj_step
      N        = 2
      nvar_jac = 2
      ! A state of norm 11.4, the reload's own at that iterate.
      kj_Y = (/ 1.13817d1, 0.0d0, 0.0d0, 0.0d0 /)
      kj_arc = probe_length_of_the_jacobian_action(kj_Y)
      call rel_row('the_probe_arc_is_the_scaled_interval', kj_arc,        &
                   sqrt(epsilon(1.0d0))*1.23817d1, 1.0d-12, nf)
      call rel_row('the_probe_arc_at_the_stalling_iterate', kj_arc,       &
                   1.845012d-7, 1.0d-5, nf)
      ! Nine decades above the reproducibility of the residual, so the
      ! action is not the rounding of it.
      call log_row('the_probe_arc_stands_above_the_residual_noise',       &
                   kj_arc/(1.0d-15*sqrt(sum(kj_Y*kj_Y))) .gt. 1.0d6,     &
                   .true., nf)
      ! And the step of that iterate is three decades inside the arc.
      kj_step = 5.0d-10
      call rel_row('the_stalling_step_is_a_thousandth_of_the_arc',        &
                   kj_step/kj_arc, 2.71d-3, 1.0d-2, nf)
      end subroutine jacobian_probe_length_rows

      ! ------------------------------------------------------!

      subroutine restart_trigger_rows(nf)
      ! WHEN A RESTART OF THE TRUST-REGION STATE IS DUE. The statistic is
      ! the stagnation detector's own: consecutive outer iterations in which
      ! the best iterate did not improve. The two sequences below are the
      ! solve's own bookkeeping of that count -- an iteration that leaves
      ! the iterate where it is raises it, an iteration that improves the
      ! best iterate zeroes it -- so what is asserted is the behavior on a
      ! run of refused steps and on a run of accepted improving ones.
      integer, intent(inout) :: nf
      integer :: kr_it, kr_since, kr_fired, kr_first
      kr_fired = c_setenv('EXHALE_TR_RESTART_STALL'//c_null_char,          &
                          '3'//c_null_char, 1)
      call int_row('restart_stall_environment_is_writable', kr_fired, 0, nf)
      call read_trust_region_restart_controls
      ! Twelve iterations, none of which improves the best iterate: the
      ! trigger has to fire, and first at the third.
      kr_since = 0;  kr_fired = 0;  kr_first = 0
      do kr_it = 1, 12
         kr_since = kr_since + 1
         if (trust_region_restart_is_due(kr_it, kr_since)) then
            kr_fired = kr_fired + 1
            if (kr_first .eq. 0) kr_first = kr_it
            kr_since = 0
         endif
      enddo
      call int_row('restart_fires_on_iterations_without_improvement',      &
                   kr_first, 3, nf)
      call int_row('restart_fires_once_every_three_such_iterations',       &
                   kr_fired, 4, nf)
      ! Twelve iterations, every one of which improves the best iterate:
      ! the count never reaches the threshold and the trigger never fires.
      kr_since = 0;  kr_fired = 0
      do kr_it = 1, 12
         kr_since = 0
         if (trust_region_restart_is_due(kr_it, kr_since))                 &
            kr_fired = kr_fired + 1
      enddo
      call int_row('restart_does_not_fire_on_improving_steps', kr_fired,   &
                   0, nf)
      ! Disarmed is the default: with the trigger off no count fires it.
      kr_fired = c_setenv('EXHALE_TR_RESTART_STALL'//c_null_char,          &
                          ''//c_null_char, 1)
      call read_trust_region_restart_controls
      call log_row('restart_trigger_is_disarmed_by_default',               &
                   trust_region_restart_is_due(50, 500), .false., nf)
      ! The forced iteration is one iteration and not a period.
      kr_fired = c_setenv('EXHALE_TR_RESTART_AT'//c_null_char,             &
                          '7'//c_null_char, 1)
      call read_trust_region_restart_controls
      call log_row('restart_forced_at_the_named_iteration',                &
                   trust_region_restart_is_due(7, 0), .true., nf)
      call log_row('restart_not_forced_at_the_iteration_before',           &
                   trust_region_restart_is_due(6, 0), .false., nf)
      call log_row('restart_not_forced_at_the_iteration_after',            &
                   trust_region_restart_is_due(8, 0), .false., nf)
      kr_fired = c_setenv('EXHALE_TR_RESTART_AT'//c_null_char,             &
                          ''//c_null_char, 1)
      call read_trust_region_restart_controls
      end subroutine restart_trigger_rows

      ! ------------------------------------------------------!

      subroutine restart_reset_rows(nf)
      ! WHAT A RESTART RESETS, AND WHAT IT LEAVES ALONE. The iterate is not
      ! an argument of the routine; the array below stands for it and is
      ! asserted bitwise unchanged across the call, which is the invariant
      ! the signature states.
      integer, intent(inout) :: nf
      real*8  :: kr_delta, kr_dmax, kr_dtau, kr_iterate(4), kr_before(4)
      integer :: kr_since, kr_nodesc, kr_restarts
      logical :: kr_closure, kr_window
      call read_trust_region_restart_controls
      tr_reset_wanted = .true.
      kr_delta = 3.25d0;  kr_dmax = 3.25d2;  kr_dtau = 1.0d7
      kr_since = 11;  kr_nodesc = 4
      kr_closure = .false.;  kr_window = .false.
      kr_iterate = (/ 1.0d0, -2.0d0, 3.5d0, 7.25d0 /)
      kr_before  = kr_iterate
      kr_restarts = n_tr_restart
      call restart_the_trust_region_state(kr_delta, kr_dmax, kr_dtau,      &
                     1.0d0, kr_since, kr_nodesc, kr_closure, kr_window)
      call log_row('restart_removes_the_radius', kr_delta .lt. 0.0d0,      &
                   .true., nf)
      call log_row('restart_removes_the_ceiling', kr_dmax .lt. 0.0d0,      &
                   .true., nf)
      call rel_row('restart_returns_the_shift_to_its_opening_value',       &
                   kr_dtau, 1.0d0, 0.0d0, nf)
      call log_row('restart_arms_the_closure_map', kr_closure, .true., nf)
      call log_row('restart_arms_the_acceptance_window', kr_window,        &
                   .true., nf)
      call int_row('restart_clears_the_stagnation_counts',                 &
                   kr_since + kr_nodesc, 0, nf)
      call int_row('restart_is_counted', n_tr_restart - kr_restarts, 1, nf)
      call rel_row('restart_leaves_the_iterate_unchanged',                 &
                   maxval(abs(kr_iterate - kr_before)), 0.0d0, 0.0d0, nf)
      ! One reset alone touches one part of the state. The radius kind is
      ! taken and the other four are not.
      tr_reset_wanted = .false.
      tr_reset_wanted(tr_reset_radius_and_ceiling) = .true.
      kr_delta = 3.25d0;  kr_dmax = 3.25d2;  kr_dtau = 1.0d7
      kr_since = 11;  kr_nodesc = 4
      kr_closure = .false.;  kr_window = .false.
      call restart_the_trust_region_state(kr_delta, kr_dmax, kr_dtau,      &
                     1.0d0, kr_since, kr_nodesc, kr_closure, kr_window)
      call log_row('radius_reset_alone_removes_the_radius',                &
                   kr_delta .lt. 0.0d0, .true., nf)
      call rel_row('radius_reset_alone_leaves_the_shift',                  &
                   kr_dtau, 1.0d7, 0.0d0, nf)
      call log_row('radius_reset_alone_leaves_the_closure_map',            &
                   kr_closure, .false., nf)
      call log_row('radius_reset_alone_leaves_the_acceptance_window',      &
                   kr_window, .false., nf)
      call int_row('radius_reset_alone_leaves_the_stagnation_counts',      &
                   kr_since + kr_nodesc, 15, nf)
      ! And the shift kind alone leaves the radius where it was.
      tr_reset_wanted = .false.
      tr_reset_wanted(tr_reset_pseudo_transient) = .true.
      kr_delta = 3.25d0;  kr_dmax = 3.25d2;  kr_dtau = 1.0d7
      call restart_the_trust_region_state(kr_delta, kr_dmax, kr_dtau,      &
                     1.0d0, kr_since, kr_nodesc, kr_closure, kr_window)
      call rel_row('shift_reset_alone_leaves_the_radius', kr_delta,        &
                   3.25d0, 0.0d0, nf)
      call rel_row('shift_reset_alone_returns_the_shift', kr_dtau, 1.0d0,  &
                   0.0d0, nf)
      tr_reset_wanted = .true.
      end subroutine restart_reset_rows

      ! ------------------------------------------------------!

      subroutine krylov_restart_cycle_rows(nf)
      ! A KRYLOV SOLVE OF ONE CYCLE IS THE SINGLE CYCLE THE SOLVER HAS
      ! ALWAYS TAKEN, and a second cycle of the SAME subspace size lowers
      ! the residual the first one left.
      !
      ! The operator is a 3x3 whose Krylov space from this right-hand side
      ! is the whole of R^3, and the subspace is cut to one: a single cycle
      ! then leaves a residual above the tolerance asked for, which is the
      ! shape of the plateau this apparatus was written for, and further
      ! cycles of the same one-dimensional space are what reduce it.
      integer, intent(inout) :: nf
      real*8, allocatable :: kr_A(:,:), kr_b(:), kr_x1(:), kr_xn(:)
      real*8, allocatable :: kr_Y(:), kr_F0(:), kr_D(:), kr_Drow(:)
      real*8, allocatable :: kr_abf(:,:), kr_Ax(:)
      integer, allocatable :: kr_ipiv(:)
      real*8  :: kr_fsp(1-Ng:1+Ng,n_species)
      integer :: kr_it1, kr_itn, kr_out1, kr_outn
      real*8  :: kr_rel1, kr_reln, kr_sn1, kr_snn
      N        = 1
      nvar_jac = 3
      kl_jac   = 0
      ku_jac   = 0
      allocate(kr_A(3,3), kr_b(3), kr_x1(3), kr_xn(3))
      allocate(kr_Y(3), kr_F0(3), kr_D(3), kr_Drow(3), kr_ipiv(3))
      allocate(kr_abf(1,3), kr_Ax(3))
      kr_A = 0.0d0
      kr_A(1,1) = 4.0d0;  kr_A(1,2) = 1.0d0
      kr_A(2,1) = 1.0d0;  kr_A(2,2) = 3.0d0;  kr_A(2,3) = 1.0d0
      kr_A(3,2) = 1.0d0;  kr_A(3,3) = 2.0d0
      kr_b  = (/ 1.0d0, 1.0d0, 1.0d0 /)
      kr_Y  = 0.0d0;  kr_F0 = 0.0d0
      kr_D  = 1.0d0;  kr_Drow = 1.0d0
      kr_ipiv = 0;  kr_abf = 0.0d0;  kr_fsp = 0.0d0
      call read_trust_region_restart_controls
      gm_restart_cycles = 1
      call linear_operator_set_for_test(kr_A, lin_test_dense)
      call pgmres_with_restarts(kr_Y, kr_F0, kr_fsp, kr_D, kr_Drow,        &
                  kr_abf, kr_ipiv, 0.0d0, kr_b, kr_x1, 1, 1.0d-8,          &
                  kr_it1, kr_out1, kr_rel1, kr_sn1, kr_Ax)
      call int_row('one_cycle_spends_the_whole_subspace', kr_it1, 1, nf)
      call int_row('one_cycle_ends_with_the_subspace_exhausted', kr_out1,  &
                   gm_subspace_exhausted, nf)
      ! Four cycles of the same one-dimensional subspace, and the residual
      ! of the sum is what the rows read.
      gm_restart_cycles = 4
      call pgmres_with_restarts(kr_Y, kr_F0, kr_fsp, kr_D, kr_Drow,        &
                  kr_abf, kr_ipiv, 0.0d0, kr_b, kr_xn, 1, 1.0d-8,          &
                  kr_itn, kr_outn, kr_reln, kr_snn, kr_Ax)
      call linear_operator_clear_for_test
      call int_row('four_cycles_spend_four_products', kr_itn, 4, nf)
      call log_row('restarted_cycles_lower_the_residual',                  &
                   kr_reln .lt. kr_rel1, .true., nf)
      call log_row('the_step_moves_off_the_first_cycle',                   &
                   maxval(abs(kr_xn - kr_x1)) .gt. 0.0d0, .true., nf)
      gm_restart_cycles = 1
      deallocate(kr_A, kr_b, kr_x1, kr_xn, kr_Y, kr_F0, kr_D, kr_Drow)
      deallocate(kr_ipiv, kr_abf, kr_Ax)
      end subroutine krylov_restart_cycle_rows

      ! ------------------------------------------------------!

      subroutine reorthogonalization_rows(nf)
      ! WHETHER THE ARNOLDI BASIS IS ORTHOGONALIZED TWICE, and what the
      ! measured loss of orthogonality of one cycle is on an operator this
      ! driver states (item N25). The condition for the second pass is a
      ! MEASURED loss above gm_reorthogonalization_loss_level; these rows
      ! fix the default, the arming, the level, and that a second pass over
      ! a basis that had lost nothing returns the same step.
      integer, intent(inout) :: nf
      real*8, allocatable :: ro_A(:,:), ro_b(:), ro_x1(:), ro_x2(:)
      integer :: ro_it1, ro_out1, ro_it2, ro_out2, ro_env
      real*8  :: ro_rel1, ro_sn1, ro_rel2, ro_sn2, ro_loss1, ro_loss2
      ro_env = c_setenv('EXHALE_GM_REORTHO'//c_null_char,                 &
                        ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('reorthogonalization_default_is_off',                  &
                   gm_reorthogonalize, .false., nf)
      ro_env = c_setenv('EXHALE_GM_REORTHO'//c_null_char,                 &
                        '1'//c_null_char, 1)
      call int_row('reorthogonalization_environment_is_writable', ro_env, &
                   0, nf)
      call read_species_unknown_space_controls
      call log_row('reorthogonalization_is_armed_by_the_environment',     &
                   gm_reorthogonalize, .true., nf)
      call rel_row('reorthogonalization_loss_level',                      &
                   gm_reorthogonalization_loss_level, 1.0d-8, 0.0d0, nf)
      ! A three-row operator whose Krylov space from this right-hand side
      ! is the whole space, measured with one pass and with two.
      allocate(ro_A(3,3), ro_b(3))
      ro_A = 0.0d0
      ro_A(1,1) = 4.0d0;  ro_A(1,2) = 1.0d0
      ro_A(2,1) = 1.0d0;  ro_A(2,2) = 3.0d0;  ro_A(2,3) = 1.0d0
      ro_A(3,2) = 1.0d0;  ro_A(3,3) = 2.0d0
      ro_b = (/ 1.0d0, 1.0d0, 1.0d0 /)
      gm_measure_orthogonality = .true.
      gm_reorthogonalize       = .false.
      call run_krylov(ro_A, ro_b, 3, 1.0d-12, lin_test_dense, ro_x1,      &
                      ro_it1, ro_out1, ro_rel1, ro_sn1)
      ro_loss1 = gm_orthogonality_loss_cycle
      gm_reorthogonalize = .true.
      call run_krylov(ro_A, ro_b, 3, 1.0d-12, lin_test_dense, ro_x2,      &
                      ro_it2, ro_out2, ro_rel2, ro_sn2)
      ro_loss2 = gm_orthogonality_loss_cycle
      gm_reorthogonalize       = .false.
      gm_measure_orthogonality = .false.
      call log_row('the_loss_of_orthogonality_of_one_cycle_is_measured',  &
                   ro_loss1 .ge. 0.0d0, .true., nf)
      ! On this operator the loss is below the level, so this is a basis
      ! the second pass is not needed for: the measurement decides.
      call log_row('this_basis_loses_less_than_the_level',                &
                   ro_loss1 .lt. gm_reorthogonalization_loss_level,       &
                   .true., nf)
      call log_row('the_second_pass_does_not_lose_more',                  &
                   ro_loss2 .le. ro_loss1, .true., nf)
      call log_row('the_second_pass_keeps_the_tolerance',                 &
                   ro_rel2 .le. 1.0d-12, .true., nf)
      call rel_row('the_second_pass_returns_the_same_step',               &
                   maxval(abs(ro_x2 - ro_x1)), 0.0d0, 1.0d-12, nf)
      deallocate(ro_A, ro_b, ro_x1, ro_x2)
      ro_env = c_setenv('EXHALE_GM_REORTHO'//c_null_char,                 &
                        ''//c_null_char, 1)
      call read_species_unknown_space_controls
      end subroutine reorthogonalization_rows

      ! ------------------------------------------------------!

      subroutine krylov_on_the_ball_rows(nf)
      ! THE CYCLE THAT STOPS ON THE TRUST BALL (Steihaug-Toint truncation,
      ! item N25). What the rows fix: the default, the arming, where the
      ! segment from an iterate inside the ball to one outside it crosses
      ! the boundary, that the point returned lies ON the ball, that the
      ! image handed back is the operator's image OF THAT POINT -- so the
      ! trust region's model is the model of the step it takes -- that the
      ! outcome is named, and that a ball wide enough to hold the
      ! untruncated step changes nothing.
      integer, intent(inout) :: nf
      real*8, allocatable :: bl_A(:,:), bl_b(:), bl_x(:), bl_xf(:)
      real*8, allocatable :: bl_Ax(:), bl_Y(:), bl_F0(:), bl_D(:)
      real*8, allocatable :: bl_Drow(:), bl_abf(:,:)
      integer, allocatable :: bl_ipiv(:)
      real*8  :: bl_fsp(1-Ng:1+Ng,n_species)
      real*8  :: bl_delta, bl_rel, bl_snorm, bl_relf, bl_snormf
      real*8  :: bl_in(3), bl_out(3), bl_t, bl_pt(3)
      integer :: bl_it, bl_out_code, bl_itf, bl_outf, bl_env
      bl_env = c_setenv('EXHALE_KRYLOV_ON_THE_BALL'//c_null_char,         &
                        ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('krylov_on_the_ball_default_is_off',                   &
                   krylov_truncated_on_the_trust_ball, .false., nf)
      bl_env = c_setenv('EXHALE_KRYLOV_ON_THE_BALL'//c_null_char,         &
                        '1'//c_null_char, 1)
      call int_row('krylov_on_the_ball_environment_is_writable', bl_env,  &
                   0, nf)
      call read_species_unknown_space_controls
      call log_row('krylov_on_the_ball_is_armed_by_the_environment',      &
                   krylov_truncated_on_the_trust_ball, .true., nf)
      ! From the origin to a point of length two, the unit ball is crossed
      ! half way.
      bl_in  = 0.0d0
      bl_out = (/ 2.0d0, 0.0d0, 0.0d0 /)
      bl_t   = the_point_where_the_step_leaves_the_ball(bl_in, bl_out,    &
                                                        1.0d0)
      call rel_row('the_boundary_of_the_ball_from_the_origin', bl_t,      &
                   0.5d0, 0.0d0, nf)
      ! And from an iterate already inside, the point returned has the
      ! radius for its length.
      bl_in  = (/ 0.6d0, 0.3d0, 0.0d0 /)
      bl_out = (/ 3.0d0, -2.0d0, 1.0d0 /)
      bl_t   = the_point_where_the_step_leaves_the_ball(bl_in, bl_out,    &
                                                        1.0d0)
      bl_pt  = bl_in + bl_t*(bl_out - bl_in)
      call rel_row('the_boundary_point_has_the_radius_for_its_length',    &
                   sqrt(sum(bl_pt*bl_pt)), 1.0d0, 1.0d-14, nf)
      ! The production cycle on a stated operator, with the image asked
      ! for. The first unknown of this operator is a thousand times softer
      ! than the others, so the untruncated step is a thousand long.
      N        = 1
      nvar_jac = 3
      kl_jac   = 0
      ku_jac   = 0
      allocate(bl_A(3,3), bl_b(3), bl_x(3), bl_xf(3), bl_Ax(3))
      allocate(bl_Y(3), bl_F0(3), bl_D(3), bl_Drow(3), bl_ipiv(3))
      allocate(bl_abf(1,3))
      bl_A = 0.0d0
      bl_A(1,1) = 1.0d-3
      bl_A(2,2) = 1.0d0
      bl_A(3,3) = 2.0d0
      bl_b = (/ 1.0d0, 1.0d0, 1.0d0 /)
      bl_Y = 0.0d0;  bl_F0 = 0.0d0
      bl_D = 1.0d0;  bl_Drow = 1.0d0
      bl_ipiv = 0;  bl_abf = 0.0d0;  bl_fsp = 0.0d0
      krylov_truncated_on_the_trust_ball = .false.
      call linear_operator_set_for_test(bl_A, lin_test_dense)
      call pgmres(bl_Y, bl_F0, bl_fsp, bl_D, bl_Drow, bl_abf, bl_ipiv,    &
                  0.0d0, bl_b, bl_xf, 3, 1.0d-10, bl_itf, bl_outf,        &
                  bl_relf, bl_snormf, bl_Ax, 1.0d0)
      call log_row('an_unarmed_cycle_ignores_the_ball',                   &
                   bl_snormf .gt. 1.0d0, .true., nf)
      krylov_truncated_on_the_trust_ball = .true.
      bl_delta = 1.0d0
      call pgmres(bl_Y, bl_F0, bl_fsp, bl_D, bl_Drow, bl_abf, bl_ipiv,    &
                  0.0d0, bl_b, bl_x, 3, 1.0d-10, bl_it, bl_out_code,      &
                  bl_rel, bl_snorm, bl_Ax, bl_delta)
      call int_row('the_truncated_cycle_names_its_outcome', bl_out_code,  &
                   gm_stopped_on_the_trust_ball, nf)
      call log_row('the_point_returned_is_inside_the_ball',               &
                   sqrt(sum(bl_x*bl_x)) .le. bl_delta*(1.0d0 + 1.0d-12),  &
                   .true., nf)
      call rel_row('the_point_returned_is_on_the_boundary',               &
                   sqrt(sum(bl_x*bl_x)), bl_delta, 1.0d-12, nf)
      ! THE MODEL OF THAT POINT: the image handed back is what the trust
      ! region forms its predicted decrease from, so it is the operator
      ! applied to the point returned and not to the untruncated step.
      call rel_row('the_image_is_the_image_of_the_point_returned',        &
                   maxval(abs(bl_Ax - matmul(bl_A, bl_x))), 0.0d0,        &
                   1.0d-12, nf)
      call log_row('the_products_of_the_truncated_cycle_are_counted',     &
                   bl_it .ge. 1 .and. bl_it .le. 3, .true., nf)
      ! A ball wide enough to hold the untruncated step leaves the cycle
      ! exactly as it was.
      call pgmres(bl_Y, bl_F0, bl_fsp, bl_D, bl_Drow, bl_abf, bl_ipiv,    &
                  0.0d0, bl_b, bl_x, 3, 1.0d-10, bl_it, bl_out_code,      &
                  bl_rel, bl_snorm, bl_Ax, 1.0d6)
      call int_row('a_wide_ball_leaves_the_outcome_alone', bl_out_code,   &
                   bl_outf, nf)
      call rel_row('a_wide_ball_leaves_the_step_alone',                   &
                   maxval(abs(bl_x - bl_xf)), 0.0d0, 0.0d0, nf)
      call linear_operator_clear_for_test
      krylov_truncated_on_the_trust_ball = .false.
      bl_env = c_setenv('EXHALE_KRYLOV_ON_THE_BALL'//c_null_char,         &
                        ''//c_null_char, 1)
      call read_species_unknown_space_controls
      deallocate(bl_A, bl_b, bl_x, bl_xf, bl_Ax, bl_Y, bl_F0, bl_D)
      deallocate(bl_Drow, bl_ipiv, bl_abf)
      end subroutine krylov_on_the_ball_rows

      ! ------------------------------------------------------!

      subroutine row_name_rows(nf)
      ! THE ROW OF A FLAT INDEX, BY NAME (item N25). The iteration line of
      ! the solve reports the largest scaled residual, and a slot number
      ! says neither which element a species row belongs to nor, in a field
      ! of one digit, anything at all. name_of_unknown is the one place an
      ! index is turned into a cell and a quantity, and these rows read it.
      integer, intent(inout) :: nf
      character(len=48) :: nm_txt
      ! The registry is emptied FIRST: set_transported_species_rows sets the
      ! width of the unknown space from the rows it registers, so a width
      ! written before it is overwritten. Four unknowns per cell with no
      ! registered row is the state a stale label would be read in.
      call set_transported_species_rows(.false.)
      N        = 500
      nvar_jac = 4
      call name_of_unknown(nvar_jac*(500-1)+3, nm_txt)
      call text_row('the_energy_row_of_a_cell_is_named', trim(nm_txt),     &
                    'energy of cell 500', nf)
      call name_of_unknown(nvar_jac*(172-1)+1, nm_txt)
      call text_row('the_mass_row_of_a_cell_is_named', trim(nm_txt),       &
                    'mass of cell 172', nf)
      ! A species slot with no registry behind it says so rather than
      ! printing a number.
      call name_of_unknown(nvar_jac*(500-1)+4, nm_txt)
      call text_row('a_species_slot_without_a_registry_is_named',          &
                    trim(nm_txt), 'unregistered species of cell 500', nf)
      N        = 1
      nvar_jac = 3
      end subroutine row_name_rows

      ! ------------------------------------------------------!

      subroutine text_row(name, got, want, nf)
      ! One row whose measured value is a sentence.
      character(len=*), intent(in)    :: name, got, want
      integer,          intent(inout) :: nf
      if (trim(got) .eq. trim(want)) then
         write(*,'(A,A,A,A,A,A,A)') 'PASS ', name, ' measured=',          &
              trim(got), ' reference=', trim(want), ' tol=0'
      else
         write(*,'(A,A,A,A,A,A,A)') 'FAIL ', name, ' measured=',          &
              trim(got), ' reference=', trim(want), ' tol=0'
         nf = nf + 1
      endif
      end subroutine text_row

      ! ------------------------------------------------------!

      subroutine int_row(name, got, want, nf)
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: got, want
      integer,          intent(inout) :: nf
      if (got .eq. want) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=', got,   &
              ' reference=', want, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=', got,   &
              ' reference=', want, ' tol=0'
         nf = nf + 1
      endif
      end subroutine int_row

      ! ------------------------------------------------------!

      subroutine log_row(name, got, want, nf)
      character(len=*), intent(in)    :: name
      logical,          intent(in)    :: got, want
      integer,          intent(inout) :: nf
      if (got .eqv. want) then
         write(*,'(A,A,A,L1,A,L1,A)') 'PASS ', name, ' measured=', got,   &
              ' reference=', want, ' tol=0'
      else
         write(*,'(A,A,A,L1,A,L1,A)') 'FAIL ', name, ' measured=', got,   &
              ' reference=', want, ' tol=0'
         nf = nf + 1
      endif
      end subroutine log_row

      ! ------------------------------------------------------!

      subroutine rel_row(name, got, want, tol, nf)
      ! Relative where the reference is nonzero, absolute where it is zero.
      ! A measured value that is not a number fails whatever the tolerance.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: got, want, tol
      integer,          intent(inout) :: nf
      real*8  :: dev
      logical :: ok
      if (.not. ieee_is_finite(got)) then
         ok = .false.
      else if (want .eq. 0.0d0) then
         ok = (abs(got) .le. tol)
      else
         dev = abs(got - want)/abs(want)
         ok  = (dev .le. tol)
      endif
      if (ok) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'PASS ', name,        &
              ' measured=', got, ' reference=', want, ' tol=', tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'FAIL ', name,        &
              ' measured=', got, ' reference=', want, ' tol=', tol
         nf = nf + 1
      endif
      end subroutine rel_row

      end program krylov_and_dogleg_tests

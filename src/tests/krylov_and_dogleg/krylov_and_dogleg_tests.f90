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
                               lin_test_dense_action_drifts,              &
                               lin_test_drift_after,                      &
                               model_row_equilibration_on,                &
                               model_row_equilibration,                   &
                               model_rows_equilibrated,                   &
                               unit_infinity_norm_row_scaling_of_the_band,&
                               gm_step_by_its_true_residual,              &
                               gm_true_residual_turned,                   &
                               gm_true_residual_first,                    &
                               set_transported_species_rows,              &
                               species_unknowns_outside_their_bounds,     &
                               hold_the_step_where_it_leaves_the_species_box, &
                               species_box_set_for_test,                  &
                               species_box_clear_for_test,                &
                               probe_length_of_the_jacobian_action,       &
                               probe_step_on_the_column_scales,           &
                               jv_probe_on_the_column_scales,             &
                               jv_additivity_on, jv_additivity_here,      &
                               resid_jump_scan_on, resid_jump_scan_here,  &
                               median_of_the_magnitudes,                  &
                               steps_in_a_sampled_row,                    &
                               jv_probe_length_factor,                    &
                               state_column_scale,                        &
                               ritz_values_of_a_matrix_free_operator,     &
                               column_split_against_the_band,             &
                               unknown_mask_of_the_row_class,             &
                               row_class_shares_of_a_vector,              &
                               precond_spectrum_on, krylov_size_scan_on,  &
                               band_difference_on, front_row_on,          &
                               species_row_term_set,                      &
                               species_row_from_its_terms,                &
                               attributed_jacobian_entry,                 &
                               gm_residual_history_on,                    &
                               jacobian_action_of_direction,              &
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
                               stationary_rows_of_the_returned_state,     &
                               hydrodynamic_distance_from_certification_by_cell, &
                               distance_from_certification,               &
                               name_of_unknown
      use certification,   only: certification_row_measure,               &
                                 cert_report, cert_evaluated,             &
                                 cert_unavailable, cert_not_applicable,   &
                                 cert_tol_mass, cert_tol_momentum,        &
                                 cert_tol_energy,                         &
                                 cert_tol_element_at, cert_tol_carrier_at,&
                                 cert_regime_wind_r
      use hydrodynamic_rows, only: generic_precision_rows_selected,       &
                 ROWS_PRODUCTION, ROWS_QUADRUPLE, ROWS_GENERIC_DOUBLE,    &
                 hydrodynamic_rows_in_double_precision
      use steady_residual_mod, only: reconstruction_continuation_rhs
      use hydrodynamic_rows_quadruple, only:                              &
                 quadruple_rows => hydrodynamic_flux_difference_and_source
      use BC_Apply, only: base_face_W, base_face_lower_W
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
      ! ... and its inverse, so a row can state what an ABSENT variable
      ! means and not only what a value means.
      interface
         integer(c_int) function c_unsetenv(nm) bind(C, name='unsetenv')
         import :: c_char, c_int
         character(kind=c_char), dimension(*), intent(in) :: nm
         end function c_unsetenv
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
      call jacobian_probe_rule_rows(nfail)
      call inner_solve_tolerance_rows(nfail)
      call residual_step_detector_rows(nfail)
      call extended_precision_rows(nfail)
      call well_balanced_operator_rows(nfail)
      call generic_arm_face_departure_rows(nfail)

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
      call spectrum_and_band_split_rows(nfail)

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
      ! ---- THE THREE MEASUREMENT HOOKS OF ITEM N35 ARE SILENT UNLESS
      !      THEIR ENVIRONMENT NAMES THEM.  Read after the controls were
      !      re-read with the environment carrying none of them, so the
      !      rows state the resolved default and not the declared one.
      call log_row('precond_spectrum_is_off_by_default',                 &
                   precond_spectrum_on, .false., nfail)
      call log_row('krylov_size_scan_is_off_by_default',                 &
                   krylov_size_scan_on, .false., nfail)
      call log_row('band_difference_is_off_by_default',                  &
                   band_difference_on, .false., nfail)
      call log_row('krylov_residual_history_is_off_by_default',          &
                   gm_residual_history_on, .false., nfail)
      ! ---- THE TWO OPTIONS OF ITEM N36 ARE OFF UNLESS NAMED, read in the
      !      same resolved state.
      call log_row('the_model_row_equilibration_is_off_by_default',      &
                   model_row_equilibration_on, .false., nfail)
      call log_row('the_true_residual_cycle_is_off_by_default',          &
                   gm_step_by_its_true_residual, .false., nfail)

      ! ---- the row scaling of the linear model, and the step the true
      !      residual chooses (N36) ----
      call model_row_scaling_rows(nfail)
      call true_residual_cycle_rows(nfail)
      call front_row_attribution_rows(nfail)

      ! ---- the test the iteration stops on (P3, D3) ----
      call stop_test_of_the_iteration_rows(nfail)

      ! ---- the distance the best-iterate ledger ranks on (P10) ----
      call judged_distance_reads_each_row_against_its_own_tolerance(nfail)
      call judged_distance_reads_the_mass_row_at_the_cells_tolerance(nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'krylov_and_dogleg: ', nfail,                &
              ' row(s) failed'
         stop 1
      endif
      write(*,'(A)') 'krylov_and_dogleg: all rows passed'

      contains

      ! ------------------------------------------------------!

      subroutine stop_test_of_the_iteration_rows(nf)
      ! WHAT THE ITERATION MAY STOP FOR SUCCESS ON.
      !
      ! The stop test of the stationary solve is the acceptance test of the
      ! state it would hand back: the acceptance gate AND every row the
      ! certification judges the system the solve carries by, each against
      ! its own tolerance. stationary_rows_of_the_returned_state is the
      ! second half of that test and the only reading of those rows the
      ! solver makes; the rows below state it on reports written here, so
      ! that the decision is asserted and not the report that produced it.
      !
      ! WHY IT DECIDES A STOP AT ALL. The gate is one number against the
      ! run's "Resid tol", 1e-8 on these fixtures, while the mass row is
      ! certified against 3e-12: a state can meet the gate with a mass row
      ! two orders above its own tolerance, which is what the molecular
      ! partitioned run MEASURED (mass 5.954e-12, gate met, refused at
      ! return; docs/solver_partition_experiment_20260911.md section 7.3).
      !
      ! No species row is registered here (set_transported_species_rows
      ! (.false.) at the top of this program), so the rows the test reads
      ! are the three hydrodynamic ones and a species balance of the
      ! inventory is an equation the solve was never given.
      integer, intent(inout) :: nf
      type(cert_report) :: kr_rep
      logical :: kr_ok, kr_none
      integer :: kr_worst, kr_read

      ! Three hydrodynamic rows, each within its own tolerance: the state
      ! is stationary in every equation this solve carries.
      call three_hydrodynamic_rows(kr_rep, 2.6d-13, 1.5d-9, 4.2d-9)
      call stationary_rows_of_the_returned_state(kr_rep, kr_ok, kr_worst, &
                                                 kr_none, kr_read)
      call log_row('rows_within_their_tolerances_stop_the_iteration',     &
                   kr_ok, .true., nf)
      call int_row('the_stop_test_reads_the_three_hydrodynamic_rows',     &
                   kr_read, 3, nf)
      call log_row('a_stopping_state_names_no_refusing_row',              &
                   kr_worst .eq. 0, .true., nf)

      ! The same state with the mass row at the measure the partitioned
      ! run returned: 5.954e-12 is below the 1e-8 the acceptance gate
      ! reads and above the 3e-12 the mass row is certified against, so
      ! the iteration may not stop here.
      call three_hydrodynamic_rows(kr_rep, 5.954d-12, 1.019d-9, 2.541d-9)
      call stationary_rows_of_the_returned_state(kr_rep, kr_ok, kr_worst, &
                                                 kr_none, kr_read)
      call log_row('a_mass_row_above_its_own_tolerance_does_not_stop',    &
                   kr_ok, .false., nf)
      call int_row('the_refusing_row_named_is_the_mass_row',              &
                   kr_worst, 1, nf)
      call log_row('a_refused_row_is_not_an_unmeasured_one',              &
                   kr_none, .false., nf)

      ! A row that could not be measured on the state refuses on its own
      ! and carries no number: an equation the run solves and cannot read
      ! is never a satisfied one.
      call three_hydrodynamic_rows(kr_rep, 2.6d-13, 1.5d-9, 4.2d-9)
      kr_rep%e(3)%status = cert_unavailable
      call stationary_rows_of_the_returned_state(kr_rep, kr_ok, kr_worst, &
                                                 kr_none, kr_read)
      call log_row('an_unmeasured_row_does_not_stop_the_iteration',       &
                   kr_ok, .false., nf)
      call log_row('an_unmeasured_row_is_reported_as_unmeasured',         &
                   kr_none, .true., nf)

      ! A non-finite row refuses whatever its number says.
      call three_hydrodynamic_rows(kr_rep, 2.6d-13, 1.5d-9, 4.2d-9)
      kr_rep%e(2)%finite = .false.
      call stationary_rows_of_the_returned_state(kr_rep, kr_ok, kr_worst, &
                                                 kr_none, kr_read)
      call log_row('a_nonfinite_row_does_not_stop_the_iteration',         &
                   kr_ok, .false., nf)

      ! AND WHAT THE TEST DOES NOT READ: a balance the registry does not
      ! carry is an equation this solve was never given, so the
      ! certification measures and reports it and it decides no stop.
      call three_hydrodynamic_rows(kr_rep, 2.6d-13, 1.5d-9, 4.2d-9)
      kr_rep%n                  = 4
      kr_rep%e(4)%name          = 'carrier balance H2'
      kr_rep%e(4)%status        = cert_evaluated
      kr_rep%e(4)%row_max       = 7.0d-2
      kr_rep%e(4)%tol           = 1.0d-5
      kr_rep%e(4)%finite        = .true.
      kr_rep%e(4)%within_tol    = .false.
      kr_rep%e(4)%jworst        = 292
      call stationary_rows_of_the_returned_state(kr_rep, kr_ok, kr_worst, &
                                                 kr_none, kr_read)
      call log_row('a_balance_the_solve_does_not_carry_does_not_refuse',  &
                   kr_ok, .true., nf)
      call int_row('and_it_is_not_one_of_the_rows_the_test_reads',        &
                   kr_read, 3, nf)
      end subroutine stop_test_of_the_iteration_rows

      ! ------------------------------------------------------!

      subroutine judged_distance_reads_each_row_against_its_own_tolerance(nf)
      ! THE LEDGER READS THE HYDRODYNAMIC ROWS AGAINST THEIR OWN
      ! CERTIFICATION TOLERANCES (item P10).
      !
      ! The best iterate of a stationary solve is the one closest to being
      ! certified, and a state is certified row by row: the mass row against
      ! 3e-12, the momentum row against 1e-8, the energy row against 1e-6
      ! (certification.f90), and every species row the registry carries
      ! against 1e-5 in the wind. distance_from_certification is what makes
      ! those one comparable quantity, and d < 1 is exactly the condition
      ! stationary_rows_of_the_returned_state states.
      !
      ! WHAT A SINGLE TOLERANCE IN THE HYDRODYNAMIC SLOT LOSES. |R| is the
      ! LARGEST of the three row measures with no tolerance in it, so
      ! dividing it by the run's "Resid tol" ranks iterates by whichever row
      ! is largest on residual_row_scale and not by the row that refuses.
      ! MEASURED on the HD 209458 b element reload of the partitioned route:
      ! every hydrodynamic solve ends with |R| at 1e-8, held by the energy
      ! row, while the mass row stands at about 1e-9 against its own 3e-12,
      ! so |R|/tol was about 1 for every iterate the ledger compared and the
      ! mass row, 350 times outside its tolerance, ranked nothing
      ! (docs/PLAN_20260911_partitioned_solver.md, what remains, item 3).
      !
      ! EVERY COLUMN BELOW IS ONE CELL carrying the three row measures rc of
      ! that cell, read by the one entry point the ledger and the acceptance
      ! both read, with a rounding floor below 3e-13 so that the continuity
      ! tolerance there is the fixed 3e-12
      ! (distance_of_one_cell_under_the_fixed_tolerance). Which cell's
      ! tolerance binds when the floors of a column differ is the subject of
      ! judged_distance_reads_the_mass_row_at_the_cells_tolerance. The
      ! ratios are exact arithmetic and the references are written
      ! independently of the expression under test, so the tolerances here
      ! are at the rounding of one division.
      integer, intent(inout) :: nf
      real*8 :: rc(3), rows(3), tel, tca, d_lo, d_hi
      tel = cert_tol_element_at(cert_regime_wind_r)
      tca = cert_tol_carrier_at(cert_regime_wind_r)
      ! The three tolerances this reads, stated so that a row below can be
      ! read without opening another file.
      call rel_row('the_mass_tolerance_the_distance_divides_by',          &
                   cert_tol_mass, 3.0d-12, 0.0d0, nf)
      call rel_row('the_momentum_tolerance_the_distance_divides_by',      &
                   cert_tol_momentum, 1.0d-8, 0.0d0, nf)
      call rel_row('the_energy_tolerance_the_distance_divides_by',        &
                   cert_tol_energy, 1.0d-6, 0.0d0, nf)

      ! ---- THE MEASURED STATE OF THE ELEMENT RELOAD ----
      ! mass 1e-9 against 3e-12, momentum 1e-13 against 1e-8, energy 1e-8
      ! against 1e-6: the mass row is 333 times outside its tolerance, the
      ! other two are far inside theirs, and |R| is 1e-8.
      rc = (/ 1.0d-9, 1.0d-13, 1.0d-8 /)
      call rel_row('the_mass_row_holds_the_hydrodynamic_distance',        &
                   distance_of_one_cell_under_the_fixed_tolerance(rc),    &
                   3.333333333333333d+02, 1.0d-14, nf)
      ! And it is NOT |R| over the run's residual target, which is where
      ! the ledger used to read this slot: |R| is the largest of the three
      ! measures, 1e-8, and on these fixtures "Resid tol" is 1e-8 too, so
      ! that reading returns 1 and calls the state as good as certified.
      call rel_row('the_superseded_reading_was_one',                      &
                   maxval(rc)/1.0d-8, 1.0d0, 0.0d0, nf)
      call log_row('and_the_state_is_not_within_a_factor_two_of_'//       &
                   'certification',                                       &
                   distance_of_one_cell_under_the_fixed_tolerance(rc)     &
                   .gt. 3.0d2, .true., nf)

      ! ---- d < 1 IS THE CERTIFICATION CONDITION ITSELF ----
      ! Every row inside its own tolerance.
      rc = (/ 1.0d-12, 1.0d-9, 1.0d-7 /)
      call log_row('every_hydrodynamic_row_inside_its_tolerance_is_'//    &
                   'below_one',                                           &
                   distance_of_one_cell_under_the_fixed_tolerance(rc)     &
                   .lt. 1.0d0, .true., nf)
      ! One row outside its own tolerance, the other two inside, and the
      ! plain maximum of the three unchanged: the mass row alone refuses.
      rc = (/ 1.0d-11, 1.0d-9, 1.0d-7 /)
      call log_row('one_row_outside_its_tolerance_is_above_one',          &
                   distance_of_one_cell_under_the_fixed_tolerance(rc)     &
                   .lt. 1.0d0, .false., nf)
      call rel_row('and_that_row_is_the_one_the_distance_reports',        &
                   distance_of_one_cell_under_the_fixed_tolerance(rc),    &
                   3.333333333333333d+00, 1.0d-14, nf)

      ! ---- WHICH OF TWO ITERATES OF EQUAL |R| THE LEDGER PREFERS ----
      !
      ! The ledger of solve_steady_jfnk keeps the iterate of smaller
      ! distance (keep_this_iterate) and hands that one back; the choice
      ! itself needs a whole solve and is not stated here, so what is
      ! asserted is the ORDERING the ledger reads. Two states with the same
      ! energy row, hence the same |R|, whose mass rows differ by a decade:
      ! a ranking on |R| cannot separate them at all, and the distance
      ! prefers the smaller mass row by the same decade.
      rc   = (/ 1.0d-9,  1.0d-13, 1.0d-8 /)
      d_hi = distance_of_one_cell_under_the_fixed_tolerance(rc)
      rc   = (/ 1.0d-10, 1.0d-13, 1.0d-8 /)
      d_lo = distance_of_one_cell_under_the_fixed_tolerance(rc)
      call rel_row('two_iterates_of_equal_plain_residual',                &
                   maxval((/ 1.0d-9, 1.0d-13, 1.0d-8 /))                  &
                   - maxval((/ 1.0d-10, 1.0d-13, 1.0d-8 /)),              &
                   0.0d0, 0.0d0, nf)
      call log_row('the_smaller_mass_row_is_the_smaller_distance',        &
                   d_lo .lt. d_hi, .true., nf)
      call rel_row('and_by_the_decade_that_separates_the_two_mass_rows',  &
                   d_hi/d_lo, 1.0d1, 1.0d-14, nf)

      ! ---- THE SPECIES SLOTS ARE UNCHANGED, AND THE HYDRODYNAMIC ONE
      !      ARRIVES ALREADY DIVIDED ----
      call rel_row('the_element_tolerance_the_distance_divides_by',       &
                   tel, 1.0d-5, 0.0d0, nf)
      call rel_row('the_carrier_tolerance_the_distance_divides_by',       &
                   tca, 1.0d-5, 0.0d0, nf)
      rows = (/ 0.0d0, 4.0d0*tel, 0.0d0 /)
      call rel_row('judged_distance_reads_the_element_row',               &
                   distance_from_certification(rows), 4.0d0, 0.0d0, nf)
      rows = (/ 0.0d0, 0.0d0, 8.0d0*tca /)
      call rel_row('judged_distance_reads_the_carrier_row',               &
                   distance_from_certification(rows), 8.0d0, 0.0d0, nf)
      ! The hydrodynamic slot is the distance formed above and is taken as
      ! it stands, with no second division.
      rc   = (/ 1.0d-9, 1.0d-13, 1.0d-8 /)
      rows = (/ distance_of_one_cell_under_the_fixed_tolerance(rc),       &
                0.0d0, 0.0d0 /)
      call rel_row('judged_distance_takes_the_hydrodynamic_slot_as_'//    &
                   'a_distance',                                          &
                   distance_from_certification(rows),                     &
                   3.333333333333333d+02, 1.0d-14, nf)
      ! And together it is the WORST of the three, not a sum and not an
      ! average: a state is certified when every row is within its own
      ! tolerance.
      rows = (/ 1.0d0/3.0d0, 4.0d0*tel, 2.0d0*tca /)
      call rel_row('judged_distance_is_the_worst_row',                    &
                   distance_from_certification(rows), 4.0d0, 0.0d0, nf)
      rows = (/ 1.0d0/3.0d0, 0.5d0*tel, 0.5d0*tca /)
      call log_row('judged_distance_below_one_is_certified',              &
                   distance_from_certification(rows) .lt. 1.0d0,          &
                   .true., nf)
      end subroutine judged_distance_reads_each_row_against_its_own_tolerance

      ! ------------------------------------------------------!

      real*8 function distance_of_one_cell_under_the_fixed_tolerance(rc)  &
                      result(d)
      ! THE HYDRODYNAMIC DISTANCE OF A ONE-CELL COLUMN whose continuity
      ! tolerance is the fixed cert_tol_mass: the three row measures rc of
      ! one cell, handed to the one entry point the ledger and the
      ! acceptance read, with a rounding floor of 1e-30. Ten times that
      ! floor is far below 3e-12, so the cell's tolerance is the fixed value
      ! (mass_row_cell_verdict) and the reading is
      ! max(rc(1)/3e-12, rc(2)/1e-8, rc(3)/1e-6).
      real*8, dimension(3), intent(in) :: rc
      real*8 :: q(3,1), fl(1)
      q(1:3,1) = rc
      fl(1)    = 1.0d-30
      d = hydrodynamic_distance_from_certification_by_cell(1, q, fl)
      end function distance_of_one_cell_under_the_fixed_tolerance

      ! ------------------------------------------------------!

      subroutine judged_distance_reads_the_mass_row_at_the_cells_tolerance(nf)
      ! THE LEDGER READS THE CONTINUITY ROW AGAINST THE TOLERANCE OF THE
      ! CELL THAT ROW SITS IN (item P17).
      !
      ! That tolerance is a function of the cell: 3e-12 where the cell's
      ! arithmetic resolves the row, and ten times the cell's own rounding
      ! floor where it does not, the rounding of a flux difference in a
      ! quasi-hydrostatic layer being about eps/Mach of the flux and not eps
      ! of it (mass_row_cell_verdict and cert_mass_round_margin,
      ! certification.f90; the floors are MEASURED in
      ! docs/certification_tolerance_anchoring_20260910.md, anchor 6). The
      ! verdict is therefore taken CELL BY CELL, and the cell that binds is
      ! not in general the cell of the largest measure.
      !
      ! THE COLUMN BELOW IS THE ONE THE HD 209458 b ELEMENT RELOAD PRESENTS.
      ! A BASE cell carrying a mass row of 1e-9 with a rounding floor of
      ! 1e-9, hence a tolerance of 1e-8 which it stands a decade inside; and
      ! a WIND cell carrying 1e-11 with a floor decades below 3e-13, hence
      ! the fixed tolerance, which it stands 3.33 times outside. The
      ! distance the acceptance applies is the wind cell's 3.33. A maximum
      ! taken over cells FIRST and divided afterwards reports the base
      ! cell's 333, which ranks iterates by a row the acceptance admits:
      ! the maximum of a ratio is not the ratio of the maxima once the
      ! divisor moves with the cell.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 3
      real*8 :: q(3,nc), fl(nc), flfix(nc), d, dfixed

      ! ---- THE TWO CELLS, AND WHICH OF THEM THE DISTANCE REPORTS ----
      q = 0.0d0;  fl = 1.0d-30
      q(1,1) = 1.0d-9;   fl(1) = 1.0d-9     ! base: tolerance 1e-8
      q(1,2) = 1.0d-11;  fl(2) = 1.0d-14    ! wind: tolerance 3e-12
      d = hydrodynamic_distance_from_certification_by_cell(nc, q, fl)
      call rel_row('the_wind_cell_holds_the_mass_rows_distance',          &
                   d, 3.333333333333333d+00, 1.0d-14, nf)
      call log_row('and_the_hydrodynamic_rows_are_not_certified',         &
                   d .lt. 1.0d0, .false., nf)
      ! THE SAME COLUMN WITH THE CELLS' FLOORS DISCARDED, so that every
      ! cell is read against cert_tol_mass, the FLOOR of the continuity
      ! tolerance and the tolerance of no particular cell: the base cell's
      ! measure over 3e-12, a hundred times the condition the acceptance
      ! applies here. Read that way a column stands at or above its
      ! distance, because a cell's tolerance is never below cert_tol_mass.
      flfix  = 1.0d-30
      dfixed = hydrodynamic_distance_from_certification_by_cell           &
                  (nc, q, flfix)
      call rel_row('the_fixed_tolerance_alone_reports_the_base_cell',     &
                   dfixed, 3.333333333333333d+02, 1.0d-14, nf)
      call rel_row('a_hundredfold_above_the_distance_of_that_column',     &
                   dfixed/d, 1.0d2, 1.0d-13, nf)
      call log_row('and_the_fixed_reading_is_never_below_the_distance',   &
                   dfixed .ge. d, .true., nf)

      ! ---- THE BASE CELL ON ITS OWN ----
      q = 0.0d0;  fl = 1.0d-30
      q(1,1) = 1.0d-9;  fl(1) = 1.0d-9
      call rel_row('the_base_cell_is_a_decade_inside_its_own_tolerance',  &
                   hydrodynamic_distance_from_certification_by_cell       &
                      (nc, q, fl), 1.0d-1, 1.0d-14, nf)
      ! And the anchor is ten times the floor and no more: a row at exactly
      ! that value stands AT its tolerance and is refused.
      q(1,1) = 1.0d-8
      d = hydrodynamic_distance_from_certification_by_cell(nc, q, fl)
      call rel_row('a_row_at_ten_times_the_floor_stands_at_its_tolerance',&
                   d, 1.0d0, 1.0d-14, nf)
      call log_row('and_a_row_at_its_tolerance_is_refused',               &
                   d .lt. 1.0d0, .false., nf)

      ! ---- EVERY ROW INSIDE ITS OWN TOLERANCE IS BELOW ONE ----
      ! The base mass row at a tenth of its anchored tolerance, a wind mass
      ! row at a third of the fixed one, the momentum row at a tenth of
      ! 1e-8 and the energy row at a tenth of 1e-6.
      q = 0.0d0;  fl = 1.0d-30
      q(1,1) = 1.0d-9;   fl(1) = 1.0d-9
      q(1,2) = 1.0d-12
      q(2,3) = 1.0d-9
      q(3,3) = 1.0d-7
      d = hydrodynamic_distance_from_certification_by_cell(nc, q, fl)
      call rel_row('every_row_inside_its_own_tolerance',                  &
                   d, 3.333333333333333d-01, 1.0d-14, nf)
      call log_row('is_a_certified_hydrodynamic_state',                   &
                   d .lt. 1.0d0, .true., nf)

      ! ---- THE MOMENTUM AND ENERGY SLOTS ARE UNCHANGED ----
      ! One tolerance each for the whole column, and the distance is the
      ! worst of the three rows over any cell.
      q = 0.0d0;  fl = 1.0d-30
      q(2,2) = 4.0d-8
      call rel_row('the_momentum_row_over_its_own_tolerance',             &
                   hydrodynamic_distance_from_certification_by_cell       &
                      (nc, q, fl), 4.0d0, 0.0d0, nf)
      q = 0.0d0
      q(3,2) = 8.0d-6
      call rel_row('the_energy_row_over_its_own_tolerance',               &
                   hydrodynamic_distance_from_certification_by_cell       &
                      (nc, q, fl), 8.0d0, 0.0d0, nf)
      q(1,1) = 1.0d-9;  fl(1) = 1.0d-9
      q(1,2) = 1.0d-11; fl(2) = 1.0d-14
      q(2,2) = 4.0d-8
      call rel_row('and_the_distance_is_the_worst_of_the_three_rows',     &
                   hydrodynamic_distance_from_certification_by_cell       &
                      (nc, q, fl), 8.0d0, 0.0d0, nf)
      end subroutine judged_distance_reads_the_mass_row_at_the_cells_tolerance

      ! ------------------------------------------------------!

      subroutine three_hydrodynamic_rows(rep, mass, mom, ener)
      ! A certification report of the three hydrodynamic rows alone, with
      ! the measures given and the tolerances the certification module
      ! defines for them. Every row is evaluated and finite; within_tol is
      ! the measure against that row's own tolerance, which is the verdict
      ! the evaluator itself takes.
      type(cert_report), intent(out) :: rep
      real*8,            intent(in)  :: mass, mom, ener
      rep%n = 3
      rep%e(1)%name = 'hydrodynamic mass row'
      rep%e(2)%name = 'hydrodynamic momentum row'
      rep%e(3)%name = 'hydrodynamic energy row'
      rep%e(1)%row_max = mass;  rep%e(1)%tol = cert_tol_mass
      rep%e(2)%row_max = mom;   rep%e(2)%tol = cert_tol_momentum
      rep%e(3)%row_max = ener;  rep%e(3)%tol = cert_tol_energy
      rep%e(1)%jworst = 14;  rep%e(2)%jworst = 500;  rep%e(3)%jworst = 16
      rep%e(1:3)%status     = cert_evaluated
      rep%e(1:3)%finite     = .true.
      rep%e(1:3)%within_tol = (rep%e(1:3)%row_max .lt. rep%e(1:3)%tol)
      end subroutine three_hydrodynamic_rows

      ! ------------------------------------------------------!

      subroutine spectrum_and_band_split_rows(nf)
      ! THE RITZ VALUES OF A STATED OPERATOR, AND THE SPLIT OF A STATED
      ! COLUMN AGAINST A STATED BAND (item N35, deliverables 2 and 3).
      !
      ! WHAT THE RITZ ROW ASSERTS. An Arnoldi recursion carried to the
      ! dimension of the space is a similarity transformation: the
      ! Hessenberg is orthogonally similar to the operator and its
      ! eigenvalues ARE the operator's, not an approximation of them. The
      ! operator here is upper triangular with distinct diagonal entries
      ! spanning six decades, so its eigenvalues are its diagonal and the
      ! reference is arithmetic and not another computation. The
      ! preconditioner is the identity the test operator installs, the
      ! column and row scales are one and the pseudo-transient shift is
      ! zero, so A_z M^-1 is the stated matrix itself.
      !
      ! WHAT THE SPLIT ROWS ASSERT. The band holds |row - col| <= kl_jac
      ! and nothing beyond it, so the difference between a column and the
      ! band's padded column decomposes exactly into the part inside the
      ! band, where the two disagree, and the part outside it, which the
      ! band does not have. The row states a column whose parts are known
      ! by hand and checks the identity to rounding.
      integer, intent(inout) :: nf
      integer, parameter :: nr = 6
      real*8  :: kr_A(nr,nr), kr_Y(nr), kr_F0(nr), kr_D(nr), kr_Drow(nr)
      real*8  :: kr_abf(1,nr), kr_V(nr,nr+1), kr_wr(nr), kr_wi(nr)
      real*8  :: kr_VR(nr,nr), kr_fsp(1-Ng:1+Ng,n_species)
      real*8  :: kr_eig(nr), kr_got(nr), kr_swap, kr_worst
      real*8  :: kr_col(nr), kr_ab(4,nr), kr_mask(nr), kr_share(5)
      real*8  :: kr_in, kr_bandgap, kr_out, kr_whole
      integer :: kr_ipiv(nr), kr_kdone, kr_info, kr_i, kr_j
      integer :: kr_cells(10), kr_ncell
      integer :: kr_nsave, kr_nvarsave, kr_klsave, kr_kusave
      kr_nsave    = N
      kr_nvarsave = nvar_jac
      kr_klsave   = kl_jac
      kr_kusave   = ku_jac
      N        = 1
      nvar_jac = nr
      kl_jac   = 0
      ku_jac   = 0
      ! ---- the Ritz values of a stated operator ----
      kr_A = 0.0d0
      kr_eig = (/ 1.0d-3, 2.5d-2, 1.0d0, 4.0d0, 3.7d1, 8.0d2 /)
      do kr_i = 1, nr
         kr_A(kr_i,kr_i) = kr_eig(kr_i)
         do kr_j = kr_i+1, nr
            kr_A(kr_i,kr_j) = 0.5d0*dble(kr_i) - 0.25d0*dble(kr_j)
         enddo
      enddo
      kr_Y = 0.0d0;  kr_F0 = 0.0d0;  kr_D = 1.0d0;  kr_Drow = 1.0d0
      kr_abf = 0.0d0;  kr_ipiv = 0;  kr_fsp = 0.0d0;  kr_mask = 1.0d0
      call linear_operator_set_for_test(kr_A, lin_test_dense)
      call ritz_values_of_a_matrix_free_operator(kr_Y, kr_F0, kr_fsp,    &
                kr_D, kr_Drow, kr_abf, kr_ipiv, 0.0d0, kr_mask, nr,      &
                kr_V, kr_wr, kr_wi, kr_VR, nr, kr_kdone, kr_info)
      call linear_operator_clear_for_test
      call int_row('ritz_recursion_spans_the_whole_space', kr_kdone, nr, &
                   nf)
      call int_row('ritz_dense_eigensolver_succeeds', kr_info, 0, nf)
      ! Sorted by magnitude, so that the comparison is of two sets and not
      ! of two orderings.
      do kr_i = 1, nr
         kr_got(kr_i) = hypot(kr_wr(kr_i), kr_wi(kr_i))
      enddo
      do kr_i = 1, nr-1
         do kr_j = kr_i+1, nr
            if (kr_got(kr_j) .lt. kr_got(kr_i)) then
               kr_swap = kr_got(kr_i)
               kr_got(kr_i) = kr_got(kr_j)
               kr_got(kr_j) = kr_swap
            endif
         enddo
      enddo
      kr_worst = 0.0d0
      do kr_i = 1, nr
         kr_worst = max(kr_worst,                                        &
                        abs(kr_got(kr_i) - kr_eig(kr_i))/kr_eig(kr_i))
      enddo
      call rel_row('ritz_values_are_the_stated_eigenvalues', kr_worst,   &
                   0.0d0, 1.0d-10, nf)
      ! ---- the split of a stated column against a stated band ----
      ! Half-bandwidth one, six unknowns, the column of unknown three: the
      ! band holds rows 2, 3 and 4 of it and nothing else.
      kl_jac = 1;  ku_jac = 1
      kr_ab  = 0.0d0
      kr_col = (/ 5.0d0, 3.0d0, 4.0d0, 12.0d0, 6.0d0, 8.0d0 /)
      ! ldab = 2*kl + ku + 1 = 4, and the column of unknown three is
      ! stored at ab(kl+ku+1 + irow - 3, 3) for irow in 2, 3, 4. Left
      ! empty first, so that the disagreement inside the band is the
      ! column itself and the reference is arithmetic.
      call column_split_against_the_band(kr_col, 3, kr_ab, kr_in,        &
                                         kr_bandgap, kr_out, kr_whole)
      ! Inside the band, rows 2, 3, 4: 3^2 + 4^2 + 12^2 = 169. Outside,
      ! rows 1, 5, 6: 25 + 36 + 64 = 125. The band is empty, so the
      ! disagreement inside it is the column itself.
      call rel_row('column_split_inside_the_band', kr_in, 169.0d0,       &
                   1.0d-14, nf)
      call rel_row('column_split_outside_the_band', kr_out, 125.0d0,     &
                   1.0d-14, nf)
      call rel_row('column_split_sums_to_the_whole_difference',          &
                   kr_bandgap + kr_out, kr_whole, 1.0d-14, nf)
      call rel_row('an_empty_band_leaves_the_column_as_the_difference',  &
                   kr_bandgap, kr_in, 1.0d-14, nf)
      ! And with the band carrying the column exactly, nothing is left
      ! inside it and the whole difference is what lies outside.
      call column_split_with_a_full_band(kr_col, nf)
      ! ---- the mask and the class shares ----
      N        = 2
      nvar_jac = 3
      kl_jac   = kr_klsave
      ku_jac   = kr_kusave
      call unknown_mask_of_the_row_class(2, kr_mask(1:6))
      call rel_row('the_hydrodynamic_mask_keeps_the_three_rows',         &
                   sum(kr_mask(1:6)), 6.0d0, 1.0d-14, nf)
      call unknown_mask_of_the_row_class(1, kr_mask(1:6))
      call rel_row('the_species_mask_keeps_nothing_without_species_rows',&
                   sum(kr_mask(1:6)), 0.0d0, 1.0d-14, nf)
      ! A vector on the energy row of cell 2 alone: one class, one cell.
      kr_col(1:6) = 0.0d0
      kr_col(6)   = 2.0d0
      call row_class_shares_of_a_vector(kr_col(1:6), kr_share, kr_cells, &
                                        kr_ncell)
      call rel_row('a_vector_on_one_energy_row_is_all_energy',           &
                   kr_share(3), 1.0d0, 1.0d-14, nf)
      call int_row('a_vector_on_one_row_sits_in_one_cell', kr_ncell, 1,  &
                   nf)
      call int_row('and_that_cell_is_the_one_it_sits_in', kr_cells(1), 2,&
                   nf)
      N        = kr_nsave
      nvar_jac = kr_nvarsave
      kl_jac   = kr_klsave
      ku_jac   = kr_kusave
      end subroutine spectrum_and_band_split_rows

      ! ------------------------------------------------------!

      subroutine column_split_with_a_full_band(kb_col, nf)
      ! THE SAME SPLIT WITH THE BAND CARRYING THE COLUMN IT HOLDS. Nothing
      ! is then left inside the band and the whole difference is exactly
      ! what lies beyond it, which is the statement the measurement rests
      ! on: an outside part that is small says the band is the operator.
      real*8, dimension(6), intent(in)    :: kb_col
      integer,              intent(inout) :: nf
      real*8  :: kb_ab(2*1+1+1, 6)
      real*8  :: kb_in, kb_gap, kb_out, kb_whole
      integer :: kb_i
      kb_ab = 0.0d0
      ! icol = 3, kl = ku = 1: the stored rows are 2, 3 and 4.
      do kb_i = 2, 4
         kb_ab(kl_jac+ku_jac+1 + kb_i - 3, 3) = kb_col(kb_i)
      enddo
      call column_split_against_the_band(kb_col, 3, kb_ab, kb_in,        &
                                         kb_gap, kb_out, kb_whole)
      call rel_row('a_band_that_carries_the_column_leaves_no_gap',       &
                   kb_gap, 0.0d0, 1.0d-20, nf)
      call rel_row('and_the_whole_difference_is_then_what_lies_outside', &
                   kb_whole, kb_out, 1.0d-14, nf)
      call rel_row('which_is_the_stated_outside_part', kb_out, 125.0d0,  &
                   1.0d-14, nf)
      end subroutine column_split_with_a_full_band

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
      ! Nine decades above the reproducibility of the residual AS A
      ! LENGTH OF THE WHOLE VECTOR. That is not the same statement as
      ! "the action is not the rounding of it", which item N31 MEASURED
      ! to be false along the preconditioned Krylov directions: the arc
      ! is a length of the whole displacement, and the fraction of ITS
      ! OWN scale by which a single unknown moves is one to two decades
      ! below sqrt(epsilon) for most of them (jacobian_probe_rule_rows).
      call log_row('the_probe_arc_stands_above_the_residual_noise',       &
                   kj_arc/(1.0d-15*sqrt(sum(kj_Y*kj_Y))) .gt. 1.0d6,     &
                   .true., nf)
      ! And the step of that iterate is three decades inside the arc.
      kj_step = 5.0d-10
      call rel_row('the_stalling_step_is_a_thousandth_of_the_arc',        &
                   kj_step/kj_arc, 2.71d-3, 1.0d-2, nf)
      end subroutine jacobian_probe_length_rows

      ! ------------------------------------------------------!

      subroutine jacobian_probe_rule_rows(nf)
      ! THE TWO PROBE RULES, AND WHAT THE HOOKS OF ITEM N31 DEFAULT TO.
      !
      ! probe_length_of_the_jacobian_action fixes the LENGTH of the
      ! displacement in the unknowns themselves;
      ! probe_step_on_the_column_scales fixes it in the column-scaled
      ! coordinates, so that every unknown is displaced by about
      ! sqrt(epsilon) of its own magnitude. The two agree, to the ratio
      ! asserted below, on a direction that runs along the scales, and
      ! part on one that does not, which is the case the Krylov basis
      ! produces.
      integer, intent(inout) :: nf
      real*8  :: jp_Y(4), jp_v(4), jp_eps_norm, jp_eps_col, jp_ratio
      integer :: jp_status
      N        = 4
      nvar_jac = 1
      ! Two unknowns of order one and two four decades below them, each
      ! sitting at its own scale.
      jp_Y = (/ 1.0d0, 1.0d0, 1.0d-4, 1.0d-4 /)
      if (allocated(state_column_scale)) deallocate(state_column_scale)
      ! WITHOUT A COLUMN SCALE the column-scaled rule is the unscaled one.
      jp_v        = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
      jp_eps_norm = probe_length_of_the_jacobian_action(jp_Y)             &
                    /sqrt(sum(jp_v*jp_v))
      jp_eps_col  = probe_step_on_the_column_scales(jp_Y, jp_v)
      call rel_row('the_column_scaled_probe_is_the_length_rule_'//        &
                   'without_scales', jp_eps_col, jp_eps_norm, 1.0d-14, nf)
      ! WITH ONE, a direction along the scales moves every unknown by the
      ! same fraction of itself under either rule.
      allocate(state_column_scale(4))
      state_column_scale = jp_Y
      jp_v = jp_Y
      jp_eps_col = probe_step_on_the_column_scales(jp_Y, jp_v)
      call rel_row('the_column_scaled_probe_moves_every_unknown_by_'//    &
                   'root_epsilon',                                        &
                   jp_eps_col, sqrt(epsilon(1.0d0))*3.0d0/2.0d0,          &
                   1.0d-12, nf)
      ! The unscaled rule is within an eighth of it on THAT direction:
      ! the ratio is (1 + ||Dc^-1 Y||) ||v|| / ((1 + ||Y||) ||Dc^-1 v||)
      ! = (3/2)(1.41421/2.41421), the two rules differing only in which
      ! norm counts the unknowns that sit far below the largest.
      jp_eps_norm = probe_length_of_the_jacobian_action(jp_Y)             &
                    /sqrt(sum(jp_v*jp_v))
      call rel_row('the_two_probe_rules_agree_along_the_scales',          &
                   jp_eps_col/jp_eps_norm, 8.7868d-1, 1.0d-4, nf)
      ! AND THEY PART ON A DIRECTION THAT DOES NOT RUN ALONG THE SCALES:
      ! all the weight of this one sits on the two large unknowns, so the
      ! unscaled rule sets the arc from them and the small ones are
      ! displaced by four decades less of themselves than sqrt(epsilon).
      jp_v = (/ 1.0d0, 1.0d0, 1.0d-8, 1.0d-8 /)
      jp_eps_norm = probe_length_of_the_jacobian_action(jp_Y)             &
                    /sqrt(sum(jp_v*jp_v))
      call rel_row('the_length_rule_underdisplaces_the_small_unknown',    &
                   jp_eps_norm*jp_v(3)/state_column_scale(3)              &
                   /sqrt(epsilon(1.0d0)), 1.7071d-4, 1.0d-3, nf)
      jp_eps_col = probe_step_on_the_column_scales(jp_Y, jp_v)
      jp_ratio   = jp_eps_col/jp_eps_norm
      call rel_row('the_column_scaled_probe_lengthens_that_arc',          &
                   jp_ratio, 1.2426d0, 1.0d-3, nf)
      deallocate(state_column_scale)
      ! THE LENGTH MULTIPLIER OF THE ADDITIVITY SCAN IS EXACTLY ONE
      ! outside the scan, so the standard arc is the unscaled expression.
      call rel_row('the_probe_length_multiplier_is_one',                  &
                   jv_probe_length_factor, 1.0d0, 0.0d0, nf)
      ! BOTH HOOKS ARE OFF UNLESS THE ENVIRONMENT ARMS THEM.
      jp_status = c_setenv('EXHALE_JV_ADDITIVITY'//c_null_char,           &
                           ''//c_null_char, 1)
      jp_status = c_setenv('EXHALE_JV_COLUMN_SCALE'//c_null_char,         &
                           ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_additivity_hook_is_off_by_default',               &
                   jv_additivity_on, .false., nf)
      call log_row('the_additivity_hook_speaks_at_no_iteration_by_'//     &
                   'default', jv_additivity_here, .false., nf)
      call log_row('the_column_scaled_probe_is_off_by_default',           &
                   jv_probe_on_the_column_scales, .false., nf)
      jp_status = c_setenv('EXHALE_JV_ADDITIVITY'//c_null_char,           &
                           '1'//c_null_char, 1)
      jp_status = c_setenv('EXHALE_JV_COLUMN_SCALE'//c_null_char,         &
                           '1'//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_additivity_hook_arms_from_the_environment',       &
                   jv_additivity_on, .true., nf)
      call log_row('the_column_scaled_probe_arms_from_the_environment',   &
                   jv_probe_on_the_column_scales, .true., nf)
      ! WITH BOTH HOOKS ON, A LINEAR TEST OPERATOR IS STILL RETURNED
      ! EXACTLY: neither probe rule enters where the action is stated as a
      ! matrix, so a comparison on a linear system measures nothing.
      call probe_rule_on_a_linear_operator(nf)
      jp_status = c_setenv('EXHALE_JV_ADDITIVITY'//c_null_char,           &
                           ''//c_null_char, 1)
      jp_status = c_setenv('EXHALE_JV_COLUMN_SCALE'//c_null_char,         &
                           ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_hooks_disarm_again',                              &
                   jv_additivity_on .or. jv_probe_on_the_column_scales,   &
                   .false., nf)
      end subroutine jacobian_probe_rule_rows

      ! ------------------------------------------------------!

      subroutine inner_solve_tolerance_rows(nf)
      ! THE TWO STOPPING TOLERANCES INSIDE A RESIDUAL EVALUATION, and the
      ! hooks that replace them (item N32). The composition's cell solve
      ! stops on sqrt(dpmpar(1)) = 1.49e-8 and the post-process
      ! temperature root on a bracket of 1e-10; each is a size below
      ! which the map the residual is assembled from is a PIECEWISE
      ! function of the state, so each has to be movable to be measured.
      ! The rows assert the default is returned bit for bit when nothing
      ! is set, that the hook is honored when it names a positive number,
      ! and that a non-positive or unreadable value leaves the default.
      use ionization_equilibrium, only: composition_solve_tolerance,      &
                                        ieq_inner_tol, ieq_inner_tol_read
      use equation_T,             only: temperature_bracket_tolerance,    &
                                        teq_bracket_tol
      integer, intent(inout) :: nf
      real*8, parameter :: it_minpack = 1.4901161193847656d-8
      integer :: it_status
      ! NOTHING IS ARMED UNTIL A SWEEP READS THE HOOK: the cached value of
      ! the sweep is untouched by anything this driver does.
      call log_row('the_composition_tolerance_is_unread_in_this_driver',  &
                   ieq_inner_tol_read, .false., nf)
      call rel_row('the_composition_tolerance_cache_is_empty',            &
                   ieq_inner_tol, 0.0d0, 0.0d0, nf)
      call rel_row('the_temperature_bracket_default_is_1e_10',            &
                   teq_bracket_tol, 1.0d-10, 0.0d0, nf)
      ! UNSET: both return the default they are handed, bit for bit.
      it_status = c_setenv('EXHALE_IEQ_TOL'//c_null_char,                 &
                           ''//c_null_char, 1)
      it_status = c_setenv('EXHALE_TEQ_TOL'//c_null_char,                 &
                           ''//c_null_char, 1)
      call rel_row('the_composition_tolerance_is_minpacks_when_unset',    &
                   composition_solve_tolerance(it_minpack), it_minpack,   &
                   0.0d0, nf)
      call rel_row('the_temperature_bracket_is_the_default_when_unset',   &
                   temperature_bracket_tolerance(1.0d-10), 1.0d-10,       &
                   0.0d0, nf)
      ! ARMED: the value the environment names, exactly.
      it_status = c_setenv('EXHALE_IEQ_TOL'//c_null_char,                 &
                           '1.0e-12'//c_null_char, 1)
      it_status = c_setenv('EXHALE_TEQ_TOL'//c_null_char,                 &
                           '1.0e-13'//c_null_char, 1)
      call rel_row('the_composition_tolerance_follows_the_hook',          &
                   composition_solve_tolerance(it_minpack), 1.0d-12,      &
                   1.0d-15, nf)
      call rel_row('the_temperature_bracket_follows_the_hook',            &
                   temperature_bracket_tolerance(1.0d-10), 1.0d-13,       &
                   1.0d-15, nf)
      ! A NON-POSITIVE VALUE IS NOT A TOLERANCE: the default stands.
      it_status = c_setenv('EXHALE_IEQ_TOL'//c_null_char,                 &
                           '-1.0e-12'//c_null_char, 1)
      it_status = c_setenv('EXHALE_TEQ_TOL'//c_null_char,                 &
                           '0.0'//c_null_char, 1)
      call rel_row('a_negative_composition_tolerance_is_refused',         &
                   composition_solve_tolerance(it_minpack), it_minpack,   &
                   0.0d0, nf)
      call rel_row('a_zero_temperature_bracket_is_refused',               &
                   temperature_bracket_tolerance(1.0d-10), 1.0d-10,       &
                   0.0d0, nf)
      ! AN UNREADABLE VALUE LIKEWISE.
      it_status = c_setenv('EXHALE_IEQ_TOL'//c_null_char,                 &
                           'tight'//c_null_char, 1)
      it_status = c_setenv('EXHALE_TEQ_TOL'//c_null_char,                 &
                           'tight'//c_null_char, 1)
      call rel_row('an_unreadable_composition_tolerance_is_refused',      &
                   composition_solve_tolerance(it_minpack), it_minpack,   &
                   0.0d0, nf)
      call rel_row('an_unreadable_temperature_bracket_is_refused',        &
                   temperature_bracket_tolerance(1.0d-10), 1.0d-10,       &
                   0.0d0, nf)
      ! AND THEY DISARM AGAIN.
      it_status = c_setenv('EXHALE_IEQ_TOL'//c_null_char,                 &
                           ''//c_null_char, 1)
      it_status = c_setenv('EXHALE_TEQ_TOL'//c_null_char,                 &
                           ''//c_null_char, 1)
      call rel_row('the_composition_tolerance_disarms_again',             &
                   composition_solve_tolerance(it_minpack), it_minpack,   &
                   0.0d0, nf)
      call rel_row('the_temperature_bracket_disarms_again',               &
                   temperature_bracket_tolerance(1.0d-10), 1.0d-10,       &
                   0.0d0, nf)
      end subroutine inner_solve_tolerance_rows

      ! ------------------------------------------------------!

      subroutine residual_step_detector_rows(nf)
      ! WHAT MAKES AN INTERVAL OF A SAMPLED ROW A STEP OF THE MAP (item
      ! N33). The jump scan reads the residual at equally spaced points
      ! along one direction and asks which interval holds a discontinuity.
      ! The test it uses carries no absolute scale: an interval is a step
      ! when its difference stands a stated factor above the MEDIAN of the
      ! row's own differences, so a row of any size and any units is read
      ! by the same rule. These rows state what that means on sequences
      ! whose answer is known: a straight line and a parabola are
      ! continuous and hold none, a line with one value displaced holds
      ! exactly one and at the interval where it was put, a constant row
      ! has no median to compare against and holds none, and the factor is
      ! honored rather than assumed.
      integer, intent(inout) :: nf
      integer, parameter :: nrs = 65
      real*8, dimension(nrs)   :: rs_f
      logical, dimension(nrs-1):: rs_mark
      real*8  :: rs_med
      integer :: rs_k, rs_n, rs_status
      ! THE MEDIAN OF THE MAGNITUDES, even and odd length.
      call rel_row('the_median_of_four_magnitudes_is_the_middle_pair',    &
                   median_of_the_magnitudes((/ -3.0d0, 1.0d0, -2.0d0,     &
                                               4.0d0 /)), 2.5d0,          &
                   0.0d0, nf)
      call rel_row('the_median_of_three_magnitudes_is_the_middle',        &
                   median_of_the_magnitudes((/ 1.0d0, -5.0d0, 2.0d0 /)),  &
                   2.0d0, 0.0d0, nf)
      ! A STRAIGHT LINE HOLDS NO STEP: every difference is the median.
      do rs_k = 1, nrs
         rs_f(rs_k) = 7.0d0 + 3.0d0*real(rs_k, 8)
      enddo
      call steps_in_a_sampled_row(rs_f, 1.0d2, rs_n, rs_mark)
      call int_row('a_straight_line_holds_no_step', rs_n, 0, nf)
      ! NOR DOES A PARABOLA: its differences vary by a factor 65 across
      ! the window, which is far below the factor asked for.
      do rs_k = 1, nrs
         rs_f(rs_k) = real(rs_k, 8)**2
      enddo
      call steps_in_a_sampled_row(rs_f, 1.0d2, rs_n, rs_mark)
      call int_row('a_parabola_holds_no_step', rs_n, 0, nf)
      ! ONE DISPLACED VALUE IS TWO STEPS, one into it and one out of it,
      ! and they are the intervals it was put between.
      do rs_k = 1, nrs
         rs_f(rs_k) = 7.0d0 + 3.0d0*real(rs_k, 8)
      enddo
      rs_f(21) = rs_f(21) + 3.0d4
      call steps_in_a_sampled_row(rs_f, 1.0d2, rs_n, rs_mark)
      call int_row('one_displaced_value_is_two_steps', rs_n, 2, nf)
      call log_row('the_step_is_the_interval_into_the_displacement',      &
                   rs_mark(20), .true., nf)
      call log_row('and_the_interval_out_of_it', rs_mark(21), .true., nf)
      call log_row('the_neighboring_interval_is_not_a_step',             &
                   rs_mark(19), .false., nf)
      ! A ROW THAT DOES NOT MOVE HAS NO MEDIAN TO STAND ABOVE.
      rs_f = 4.0d0
      call steps_in_a_sampled_row(rs_f, 1.0d2, rs_n, rs_mark)
      call int_row('a_constant_row_holds_no_step', rs_n, 0, nf)
      rs_med = median_of_the_magnitudes(rs_f(1:nrs-1) - rs_f(2:nrs))
      call rel_row('and_its_median_difference_is_zero', rs_med, 0.0d0,    &
                   0.0d0, nf)
      ! THE FACTOR IS HONORED: a displacement ten times the increment is
      ! not a step at a hundred and is one at five.
      do rs_k = 1, nrs
         rs_f(rs_k) = 7.0d0 + 3.0d0*real(rs_k, 8)
      enddo
      rs_f(21) = rs_f(21) + 3.0d1
      call steps_in_a_sampled_row(rs_f, 1.0d2, rs_n, rs_mark)
      call int_row('a_tenfold_displacement_is_no_step_at_a_hundred',      &
                   rs_n, 0, nf)
      call steps_in_a_sampled_row(rs_f, 5.0d0, rs_n, rs_mark)
      call int_row('and_is_two_steps_at_five', rs_n, 2, nf)
      ! THE SCAN IS OFF UNLESS THE ENVIRONMENT ARMS IT, and it speaks at
      ! no iteration until an outer iteration selects one.
      rs_status = c_setenv('EXHALE_RESID_JUMP_SCAN'//c_null_char,         &
                           ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_jump_scan_is_off_by_default',                     &
                   resid_jump_scan_on, .false., nf)
      call log_row('the_jump_scan_speaks_at_no_iteration_by_default',     &
                   resid_jump_scan_here, .false., nf)
      rs_status = c_setenv('EXHALE_RESID_JUMP_SCAN'//c_null_char,         &
                           '1'//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_jump_scan_arms_from_the_environment',             &
                   resid_jump_scan_on, .true., nf)
      call log_row('arming_it_selects_no_iteration_by_itself',            &
                   resid_jump_scan_here, .false., nf)
      rs_status = c_setenv('EXHALE_RESID_JUMP_SCAN'//c_null_char,         &
                           ''//c_null_char, 1)
      call read_species_unknown_space_controls
      call log_row('the_jump_scan_disarms_again',                         &
                   resid_jump_scan_on, .false., nf)
      end subroutine residual_step_detector_rows

      ! ------------------------------------------------------!

      subroutine probe_rule_on_a_linear_operator(nf)
      ! THE ACTION OF A STATED MATRIX IS THAT MATRIX, whichever probe rule
      ! is armed: jacobian_action_of_direction returns matmul(A, v) where
      ! a test has set the operator, and no finite difference is taken.
      integer, intent(inout) :: nf
      real*8  :: pl_A(3,3), pl_Y(3), pl_F0(3), pl_v(3), pl_Jv(3), pl_want(3)
      real*8  :: pl_fsp(1-Ng:3+Ng,n_species)
      logical :: pl_ok
      integer :: pl_i
      N        = 3
      nvar_jac = 1
      pl_A = reshape((/ 2.0d0, 1.0d0, 0.0d0,                              &
                        1.0d0, 3.0d0, 1.0d0,                              &
                        0.0d0, 1.0d0, 4.0d0 /), (/ 3, 3 /))
      call linear_operator_set_for_test(pl_A, lin_test_dense)
      pl_Y  = (/ 1.0d0, 2.0d0, 3.0d0 /)
      pl_F0 = 0.0d0
      pl_v  = (/ 1.0d0, -1.0d0, 0.5d0 /)
      pl_fsp = 0.0d0
      call jacobian_action_of_direction(pl_Y, pl_F0, pl_fsp, pl_v,        &
                                        pl_Jv, pl_ok)
      pl_want = matmul(pl_A, pl_v)
      call log_row('the_linear_test_action_is_admissible', pl_ok,         &
                   .true., nf)
      do pl_i = 1, 3
         call rel_row('the_armed_probe_reproduces_the_linear_operator',   &
                      pl_Jv(pl_i), pl_want(pl_i), 1.0d-14, nf)
      enddo
      call linear_operator_clear_for_test
      end subroutine probe_rule_on_a_linear_operator

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

      subroutine model_row_scaling_rows(nf)
      ! THE ROW SCALING OF THE LINEAR MODEL IS A FREEDOM OF THE MODEL AND
      ! OF NOTHING ELSE (item N36, decision 20 a).
      !
      ! Solving (E A) s = E r0 with E diagonal and positive is the same
      ! linear system: where the cycle spans the whole space, the step it
      ! returns is the same step. What E changes is the norm a TRUNCATED
      ! cycle minimizes. These rows state the first half, and that the
      ! relative residual handed back is the certification one whatever
      ! rows the cycle worked in.
      integer, intent(inout) :: nf
      real*8, allocatable :: mr_A(:,:), mr_b(:), mr_x(:), mr_xe(:)
      real*8, allocatable :: mr_ab(:,:)
      integer :: mr_i, mr_j, mr_n, mr_it, mr_out, mr_oute
      real*8  :: mr_rel, mr_rele, mr_snorm, mr_snorme, mr_rowmax, mr_worst
      ! --- E on a stated band: unit infinity norm, and a power of two ---
      ! Five unknowns of one cell, one sub- and one superdiagonal, the
      ! rows spanning six decades.
      N        = 1
      nvar_jac = 5
      kl_jac   = 1
      ku_jac   = 1
      allocate(mr_ab(2*kl_jac+ku_jac+1, 5))
      mr_ab = 0.0d0
      do mr_j = 1, 5
         do mr_i = max(1, mr_j-ku_jac), min(5, mr_j+kl_jac)
            mr_ab(kl_jac+ku_jac+1 + mr_i - mr_j, mr_j) =                 &
                 1.0d0*10.0d0**(mr_i - 3)*dble(mr_j)
         enddo
      enddo
      call unit_infinity_norm_row_scaling_of_the_band(mr_ab)
      mr_worst = 0.0d0
      do mr_i = 1, 5
         mr_rowmax = 0.0d0
         do mr_j = max(1, mr_i-kl_jac), min(5, mr_i+ku_jac)
            mr_rowmax = max(mr_rowmax, abs(                              &
                 mr_ab(kl_jac+ku_jac+1 + mr_i - mr_j, mr_j)))
         enddo
         mr_worst = max(mr_worst, abs(log(mr_rowmax                      &
                    *model_row_equilibration(mr_i))/log(2.0d0)))
      enddo
      ! Every row within a factor sqrt(2) of one, which is what a power of
      ! two can reach.
      call log_row('the_row_scaling_brings_every_row_to_unit_size',      &
                   mr_worst .le. 0.5d0, .true., nf)
      mr_worst = 0.0d0
      do mr_i = 1, 5
         mr_worst = max(mr_worst, abs(model_row_equilibration(mr_i)      &
              - 2.0d0**nint(log(model_row_equilibration(mr_i))           &
                            /log(2.0d0))))
      enddo
      call rel_row('the_row_scaling_is_a_power_of_two', mr_worst, 0.0d0,  &
                   0.0d0, nf)
      deallocate(mr_ab)
      model_rows_equilibrated = .false.
      ! --- the same step, whatever the rows are scaled by ---
      ! A dense six by six whose rows span six decades, solved over the
      ! whole space: the equilibrated cycle and the unequilibrated one
      ! return one step.
      mr_n = 6
      allocate(mr_A(mr_n,mr_n), mr_b(mr_n))
      do mr_i = 1, mr_n
         do mr_j = 1, mr_n
            mr_A(mr_i,mr_j) = 10.0d0**(mr_i-3)                           &
                 *(1.0d0/dble(mr_i+mr_j) + merge(2.0d0, 0.0d0,           &
                                                 mr_i .eq. mr_j))
         enddo
         mr_b(mr_i) = 10.0d0**(mr_i-3)
      enddo
      model_row_equilibration_on = .false.
      model_rows_equilibrated    = .false.
      call run_krylov(mr_A, mr_b, mr_n, 1.0d-14, lin_test_dense, mr_x,    &
                      mr_it, mr_out, mr_rel, mr_snorm)
      ! E of the dense operator, stated here because the band the cycle is
      ! preconditioned by in these rows is the identity the test installs
      ! and carries none of the operator.
      if (allocated(model_row_equilibration))                            &
         deallocate(model_row_equilibration)
      allocate(model_row_equilibration(mr_n))
      do mr_i = 1, mr_n
         mr_rowmax = maxval(abs(mr_A(mr_i,:)))
         model_row_equilibration(mr_i) =                                 &
              2.0d0**(-nint(log(mr_rowmax)/log(2.0d0)))
      enddo
      model_row_equilibration_on = .true.
      model_rows_equilibrated    = .true.
      call run_krylov(mr_A, mr_b, mr_n, 1.0d-14, lin_test_dense, mr_xe,   &
                      mr_it, mr_oute, mr_rele, mr_snorme)
      call rel_row('the_equilibrated_solve_returns_the_same_step',        &
                   maxval(abs(mr_xe - mr_x))                              &
                   /max(maxval(abs(mr_x)), 1.0d-300), 0.0d0, 1.0d-12, nf)
      ! AND THE RESIDUAL IT REPORTS IS THE CERTIFICATION ONE, not the one
      ! it minimized: read against the operator and the right-hand side as
      ! the caller states them. Taken on a TRUNCATED cycle, where the two
      ! norms differ; on the whole space both are the rounding.
      deallocate(mr_xe)
      call run_krylov(mr_A, mr_b, 3, 1.0d-14, lin_test_dense, mr_xe,      &
                      mr_it, mr_oute, mr_rele, mr_snorme)
      call rel_row('the_equilibrated_cycle_reports_the_certification'//   &
                   '_residual', mr_rele,                                  &
                   sqrt(sum((mr_b - matmul(mr_A, mr_xe))**2))             &
                   /sqrt(sum(mr_b*mr_b)), 1.0d-8, nf)
      call log_row('the_equilibrated_cycle_minimizes_another_norm',       &
                   mr_rele .gt. 1.0d-6, .true., nf)
      model_row_equilibration_on = .false.
      model_rows_equilibrated    = .false.
      deallocate(model_row_equilibration)
      deallocate(mr_A, mr_b, mr_x, mr_xe)
      N        = 1
      nvar_jac = 3
      kl_jac   = 0
      ku_jac   = 0
      end subroutine model_row_scaling_rows

      ! ------------------------------------------------------!

      subroutine true_residual_cycle_rows(nf)
      ! THE CYCLE THAT RETURNS THE STEP ITS TRUE RESIDUAL CHOOSES (N36).
      !
      ! On a LINEAR operator the residual of the reduced problem is the
      ! residual of the step, so that stop may not move anything: the step it
      ! returns is the step the same subspace gives, and the tolerance it
      ! announces is announced only where the residual reached it. On an
      ! operator whose action drifts with the product count -- what the
      ! matrix-free secant of a nonlinear residual does over a long cycle
      ! (section N35) -- the two part company, and the cycle has to stop
      ! where the residual of its step turns upward instead of running to
      ! the end of its subspace.
      integer, intent(inout) :: nf
      real*8, allocatable :: tr_A(:,:), tr_b(:), tr_x20(:), tr_xa(:)
      real*8, allocatable :: tr_xd(:)
      integer :: tr_i, tr_j, tr_n, tr_it20, tr_out20, tr_ita, tr_outa
      integer :: tr_itd, tr_outd
      real*8  :: tr_rel20, tr_rela, tr_reld, tr_sn20, tr_sna, tr_snd
      ! A diagonal spanning the size of the system with one superdiagonal,
      ! so that the Krylov space of a right-hand side of ones is full and
      ! no breakdown cuts the cycle short of its checks.
      tr_n = 30
      allocate(tr_A(tr_n,tr_n), tr_b(tr_n))
      tr_A = 0.0d0
      do tr_i = 1, tr_n
         tr_A(tr_i,tr_i) = dble(tr_i)
         if (tr_i .lt. tr_n) tr_A(tr_i,tr_i+1) = 0.3d0
         tr_b(tr_i) = 1.0d0
      enddo
      gm_step_by_its_true_residual = .false.
      call run_krylov(tr_A, tr_b, gm_true_residual_first, 1.0d-14,        &
                      lin_test_dense, tr_x20, tr_it20, tr_out20,          &
                      tr_rel20, tr_sn20)
      ! The true-residual stop, asked for a tolerance the cycle reaches at
      ! first check: it has to stop there, name the tolerance and hand
      ! back the step of that subspace.
      gm_step_by_its_true_residual = .true.
      call run_krylov(tr_A, tr_b, tr_n - 5, 2.0d0*tr_rel20,               &
                      lin_test_dense, tr_xa, tr_ita, tr_outa, tr_rela,    &
                      tr_sna)
      call int_row('the_true_residual_cycle_names_the_tolerance',         &
                   tr_outa, gm_tolerance_reached, nf)
      call rel_row('the_true_residual_cycle_returns_the_same_step',       &
                   maxval(abs(tr_xa - tr_x20))                            &
                   /max(maxval(abs(tr_x20)), 1.0d-300), 0.0d0, 1.0d-12,   &
                   nf)
      ! On a linear operator the true residual IS the reduced one.
      call rel_row('the_true_residual_of_a_linear_map_is_the_reduced'//   &
                   '_one', tr_rela, tr_rel20, 1.0d-10, nf)
      ! --- an action that stops being one operator ---
      ! The same cycle on the drifting action: the reduced problem keeps
      ! falling, the residual of the step against the action does not, and
      ! the cycle has to say so rather than spend its whole subspace.
      tr_n = 70
      deallocate(tr_A, tr_b)
      allocate(tr_A(tr_n,tr_n), tr_b(tr_n))
      tr_A = 0.0d0
      do tr_i = 1, tr_n
         tr_A(tr_i,tr_i) = dble(tr_i)
         if (tr_i .lt. tr_n) tr_A(tr_i,tr_i+1) = 0.3d0
         tr_b(tr_i) = 1.0d0
      enddo
      call run_krylov(tr_A, tr_b, tr_n, 1.0d-14,                          &
                      lin_test_dense_action_drifts, tr_xd, tr_itd,        &
                      tr_outd, tr_reld, tr_snd)
      call int_row('a_drifting_action_turns_the_true_residual',           &
                   tr_outd, gm_true_residual_turned, nf)
      ! It stopped before the subspace was spent, and after the drift
      ! began.
      call log_row('the_cycle_stopped_where_the_residual_turned',         &
                   tr_itd .gt. lin_test_drift_after .and.                 &
                   tr_itd .lt. tr_n, .true., nf)
      call log_row('the_step_returned_is_finite',                         &
                   maxval(abs(tr_xd)) .lt. huge(1.0d0), .true., nf)
      gm_step_by_its_true_residual = .false.
      deallocate(tr_A, tr_b, tr_x20, tr_xa, tr_xd)
      N        = 1
      nvar_jac = 3
      kl_jac   = 0
      ku_jac   = 0
      end subroutine true_residual_cycle_rows

      ! ------------------------------------------------------!

      subroutine front_row_attribution_rows(nf)
      ! THE ROW THAT BINDS, TERM BY TERM (N38).
      !
      ! Two statements of the arithmetic the measurement stands on, and one
      ! of its default.
      !
      !   * A term set composes the row it was split from. The split is a
      !     statement about one equation,
      !         div(J_diffusive) + div(F_rho Y^face) - S_chemical = row ,
      !     so a set whose terms do not sum to its row describes no
      !     equation; species_row_from_its_terms is the one expression of
      !     that sum and both the measurement and this row read it.
      !   * The entry a term contributes to the scaled band is the term's
      !     own difference quotient in the band's two-sided scaling. On a
      !     term that is LINEAR in the unknown with a stated coefficient,
      !     the attributed entry has to be that coefficient times the
      !     scaling and nothing else, and a step of zero has to give zero
      !     rather than a division by it.
      !   * The hook is off unless it is asked for.
      integer, intent(inout) :: nf
      type(species_row_term_set) :: kd_t
      real*8  :: kd_c, kd_h, kd_x, kd_tc, kd_dc, kd_dr, kd_a

      call log_row('the_front_row_hook_is_off_by_default',                &
                   front_row_on, .false., nf)

      ! A stated term set: the three terms and the row they compose.
      kd_t%diffusive = -3.5d0
      kd_t%adv_total =  8.25d0
      kd_t%chemical  =  1.75d0
      kd_t%row       = kd_t%diffusive + kd_t%adv_total - kd_t%chemical
      call rel_row('a_term_set_composes_its_own_row',                     &
                   species_row_from_its_terms(kd_t), kd_t%row,            &
                   1.0d-15, nf)
      ! And the two masked faces compose the advective term, which is what
      ! the masking device gives: the divergence of a face flux vanishes
      ! wherever the mass flux it multiplies does.
      kd_t%adv_left  = -1.25d0
      kd_t%adv_right =  9.5d0
      call rel_row('the_two_faces_compose_the_advective_term',            &
                   kd_t%adv_left + kd_t%adv_right, kd_t%adv_total,        &
                   1.0d-15, nf)

      ! A term linear in the unknown with a stated coefficient.
      kd_c  = -7.25d0
      kd_x  =  0.5d0
      kd_h  =  1.0d-6
      kd_tc =  3.0d0
      kd_dc =  2.0d0
      kd_dr =  8.0d0
      kd_a  = attributed_jacobian_entry(kd_c*(kd_x + kd_h), kd_c*kd_x,    &
                                        kd_h, kd_tc, kd_dc, kd_dr)
      call rel_row('the_attributed_entry_of_a_linear_term_is_its'//       &
                   '_coefficient', kd_a, kd_c*kd_tc*kd_dc/kd_dr,          &
                   1.0d-9, nf)
      ! A step of zero is not a derivative and may not be a division.
      call rel_row('a_zero_probe_step_attributes_nothing',                &
                   attributed_jacobian_entry(kd_c*kd_x, kd_c*kd_x,        &
                            0.0d0, kd_tc, kd_dc, kd_dr), 0.0d0,           &
                   1.0d-300, nf)
      ! And a row read on no scale attributes nothing either.
      call rel_row('a_row_with_no_scale_attributes_nothing',              &
                   attributed_jacobian_entry(kd_c*(kd_x + kd_h),          &
                            kd_c*kd_x, kd_h, kd_tc, kd_dc, 0.0d0),        &
                   0.0d0, 1.0d-300, nf)

      end subroutine front_row_attribution_rows

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

      ! ------------------------------------------------------!

      subroutine well_balanced_operator_rows(nf)
      ! THE WELL-BALANCED OPTION ("Well balanced:", default off; item N37,
      ! docs/well_balanced_flux_difference_design_20260910.md).  Four
      ! statements, on the same twelve-cell stretched grid the extended
      ! precision rows use:
      !
      !   * WITH NO GRAVITY it is the base scheme.  The local
      !     equilibrium is then the constant p_j, the departure data is the
      !     ordinary pressure stencil and the face pressure measured against
      !     the equilibrium is the ordinary one, so every row must come back
      !     the same to rounding.  It is NOT asserted bitwise: the pressure
      !     jump of the Riemann problem is formed additively under the key
      !     (the difference of the two cell pressures plus the departures)
      !     and as pR - pL without it, and the two roundings differ in the
      !     last bit of the jump.  Measured on a uniform state and on a Sod
      !     shock tube, for both reconstructions.
      !   * WITH GRAVITY, on the discrete equilibrium the key defines (the
      !     state whose two neighboring equilibrium extrapolations agree at
      !     every face), the momentum row falls to the rounding level while
      !     the base scheme leaves its truncation error.  This is the
      !     operator identity the key exists for, asserted here on a second
      !     grid and a second geometry from the one
      !     src/tests/grid_and_gates/hydrostatic_residual.f90 uses.
      !   * The generic double instantiation of the kind-generic operator
      !     reproduces the production operator bit for bit WITH THE KEY ON
      !     as well, which is what makes EXHALE_RESID_QUAD=2 the same
      !     experiment under both key values.
      integer, intent(inout) :: nf
      integer, parameter :: xn = 12, xg = 2
      real*8, dimension(3,1-xg:xn+xg) :: xu, xWLa, xWRa, xdFa, xSa
      real*8, dimension(3,1-xg:xn+xg) :: xWLb, xWRb, xdFb, xSb
      real*8, dimension(1-xg:xn+xg)   :: x_rho_eq, x_p_eq
      real*8  :: x_rho, x_vel, x_pre, x_dev, x_ref, x_b0, x_cs2
      real*8  :: x_wb, x_off
      integer :: xj, xk, xs
      integer :: sv_N, sv_wmode
      logical :: sv_plm, sv_weno, sv_lam_on, sv_wb
      character(len=:), allocatable :: sv_rec, sv_flux
      real*8, dimension(3) :: sv_bfW, sv_bflW
      real*8, dimension(:), allocatable :: sv_r, sv_redg, sv_dr,          &
                                           sv_gpi, sv_gpc

      sv_N = N;  sv_wmode = weno_mode;  sv_wb = well_balanced
      sv_plm = use_plm;  sv_weno = use_weno3;  sv_lam_on = recon_lambda_on
      sv_rec = rec_method;  sv_flux = flux
      sv_bfW = base_face_W;  sv_bflW = base_face_lower_W
      if (allocated(r)) then
         allocate(sv_r(lbound(r,1):ubound(r,1)));  sv_r = r
         allocate(sv_redg(lbound(r_edg,1):ubound(r_edg,1)));  sv_redg = r_edg
         allocate(sv_dr(lbound(dr_j,1):ubound(dr_j,1)));      sv_dr  = dr_j
         allocate(sv_gpi(lbound(Gphi_i,1):ubound(Gphi_i,1))); sv_gpi = Gphi_i
         allocate(sv_gpc(lbound(Gphi_c,1):ubound(Gphi_c,1))); sv_gpc = Gphi_c
         deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      endif

      N = xn;  weno_mode = 0;  recon_lambda_on = .false.;  flux = 'ROE'
      allocate(r(1-xg:xn+xg), r_edg(1-xg:xn+xg), dr_j(1-xg:xn+xg))
      allocate(Gphi_i(1-xg:xn+xg), Gphi_c(1-xg:xn+xg))
      do xj = 1-xg, xn+xg
         r_edg(xj) = 1.0d0 + 0.05d0*dble(xj)*(1.0d0 + 0.02d0*dble(xj))
      enddo
      do xj = 2-xg, xn+xg
         r(xj)    = 0.5d0*(r_edg(xj) + r_edg(xj-1))
         dr_j(xj) = r_edg(xj) - r_edg(xj-1)
      enddo
      r(1-xg)    = r_edg(1-xg) - 0.5d0*dr_j(2-xg)
      dr_j(1-xg) = dr_j(2-xg)

      ! ---- (1) no gravity, a uniform state ----
      Gphi_i = 0.0d0
      Gphi_c = 0.0d0
      x_rho = 1.25d0;  x_vel = 0.5d0;  x_pre = 0.75d0
      do xj = 1-xg, xn+xg
         xu(1,xj) = x_rho
         xu(2,xj) = x_rho*x_vel
         xu(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = x_rho
      base_face_W(2) = x_vel
      base_face_W(3) = x_pre
      base_face_lower_W = base_face_W
      call rel_row('well_balanced_uniform_flow_is_the_base_scheme[PLM]',  &
                   key_departure_of_the_rows(xu, xn, xg, .true.),         &
                   0.0d0, 1.0d-13, nf)
      call rel_row('well_balanced_uniform_flow_is_the_base_scheme[WENO3]',&
                   key_departure_of_the_rows(xu, xn, xg, .false.),        &
                   0.0d0, 1.0d-13, nf)

      ! ---- (2) no gravity, a Sod shock tube ----
      ! Sod (1978, J. Comput. Phys. 27, 1): (rho, v, p) = (1, 0, 1) on the
      ! left and (0.125, 0, 0.1) on the right of the middle of the domain,
      ! as initial data.  The row is evaluated once on it, so what is
      ! measured is the operator at a discontinuity and not the evolution.
      do xj = 1-xg, xn+xg
         if (xj .le. xn/2) then
            x_rho = 1.0d0;    x_pre = 1.0d0
         else
            x_rho = 0.125d0;  x_pre = 0.1d0
         endif
         xu(1,xj) = x_rho
         xu(2,xj) = 0.0d0
         xu(3,xj) = x_pre/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = 1.0d0
      base_face_W(2) = 0.0d0
      base_face_W(3) = 1.0d0
      base_face_lower_W = base_face_W
      call rel_row('well_balanced_shock_tube_is_the_base_scheme[PLM]',    &
                   key_departure_of_the_rows(xu, xn, xg, .true.),         &
                   0.0d0, 1.0d-13, nf)
      call rel_row('well_balanced_shock_tube_is_the_base_scheme[WENO3]',  &
                   key_departure_of_the_rows(xu, xn, xg, .false.),        &
                   0.0d0, 1.0d-13, nf)

      ! ---- (3) with gravity, on the key's own discrete equilibrium ----
      ! phi = -b0/r, and an isothermal column with p = c^2 rho whose
      ! densities solve the face-matching condition
      !   p_j+1 + rho_j+1 (phi_c(j+1) - phi_i(j))
      !         = p_j - rho_j (phi_i(j) - phi_c(j))
      ! exactly: rho_j+1 = rho_j (c^2 - a_j)/(c^2 + b_j).
      x_b0  = 3.0d0
      x_cs2 = 1.0d0
      do xj = 1-xg, xn+xg
         Gphi_i(xj) = -x_b0/r_edg(xj)
         Gphi_c(xj) = -x_b0/r(xj)
      enddo
      x_rho_eq(1) = 1.0d0
      do xj = 1, xn+xg-1
         x_rho_eq(xj+1) = x_rho_eq(xj)                                    &
            *(x_cs2 - (Gphi_i(xj)   - Gphi_c(xj)))                        &
            /(x_cs2 + (Gphi_c(xj+1) - Gphi_i(xj)))
      enddo
      do xj = 1, 2-xg, -1
         x_rho_eq(xj-1) = x_rho_eq(xj)                                    &
            *(x_cs2 + (Gphi_c(xj)   - Gphi_i(xj-1)))                      &
            /(x_cs2 - (Gphi_i(xj-1) - Gphi_c(xj-1)))
      enddo
      x_rho_eq(1-xg) = x_rho_eq(2-xg)
      x_p_eq = x_cs2*x_rho_eq
      do xj = 1-xg, xn+xg
         xu(1,xj) = x_rho_eq(xj)
         xu(2,xj) = 0.0d0
         xu(3,xj) = x_p_eq(xj)/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = x_rho_eq(1)
      base_face_W(2) = 0.0d0
      base_face_W(3) = x_p_eq(1)
      base_face_lower_W = base_face_W

      do xs = 1, 2
         if (xs .eq. 1) then
            rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
         else
            rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
         endif
         x_ref = 0.0d0
         do xj = 3, xn-2
            x_ref = max(x_ref, x_rho_eq(xj)*x_b0/(r(xj)*r(xj)))
         enddo

         well_balanced = .true.
         call reconstruction_continuation_rhs(xu, xWLb, xWRb, xdFb, xSb)
         x_wb = 0.0d0
         do xj = 3, xn-2
            x_wb = max(x_wb, abs(xdFb(2,xj) - xSb(2,xj))/x_ref)
         enddo

         well_balanced = .false.
         call reconstruction_continuation_rhs(xu, xWLa, xWRa, xdFa, xSa)
         x_off = 0.0d0
         do xj = 3, xn-2
            x_off = max(x_off, abs(xdFa(2,xj) - xSa(2,xj))/x_ref)
         enddo

         if (xs .eq. 1) then
            call rel_row('well_balanced_momentum_row_on_its_own_'//       &
                 'equilibrium[PLM]', x_wb, 0.0d0, 1.0d-13, nf)
            call log_row('the_same_column_is_not_balanced_without_'//    &
                 'the_key[PLM]', x_off .gt. 1.0d-9, .true., nf)
         else
            call rel_row('well_balanced_momentum_row_on_its_own_'//       &
                 'equilibrium[WENO3]', x_wb, 0.0d0, 1.0d-13, nf)
            call log_row('the_same_column_is_not_balanced_without_'//    &
                 'the_key[WENO3]', x_off .gt. 1.0d-9, .true., nf)
         endif
      enddo

      ! ---- (4) the generic double instantiation, with the key ON ----
      ! Same state, same gravity, a state with structure in it so that the
      ! reconstruction limits and every branch of the assembly is reached.
      rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
      do xj = 1-xg, xn+xg
         x_rho = 1.0d0/r(xj)**3
         x_vel = 0.3d0*r(xj)
         x_pre = 0.8d0/r(xj)**4
         xu(1,xj) = x_rho
         xu(2,xj) = x_rho*x_vel
         xu(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = 1.0d0/r(1)**3
      base_face_W(2) = 0.3d0*r(1)
      base_face_W(3) = 0.8d0/r(1)**4
      base_face_lower_W = base_face_W

      well_balanced = .true.
      call reconstruction_continuation_rhs(xu, xWLa, xWRa, xdFa, xSa)
      call hydrodynamic_rows_in_double_precision(xu, xWLb, xWRb, xdFb, xSb)
      x_dev = 0.0d0
      do xj = 2-xg, xn+xg
         do xk = 1, 3
            if (xdFb(xk,xj) .ne. xdFa(xk,xj))                             &
               x_dev = max(x_dev, abs(xdFb(xk,xj)-xdFa(xk,xj)))
            if (xSb(xk,xj) .ne. xSa(xk,xj))                               &
               x_dev = max(x_dev, abs(xSb(xk,xj)-xSa(xk,xj)))
         enddo
      enddo
      call rel_row('generic_double_is_the_operator_bitwise_with_the_key', &
                   x_dev, 0.0d0, 0.0d0, nf)

      ! ---- put the globals back ----
      well_balanced = sv_wb
      deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      if (allocated(sv_r)) then
         allocate(r(lbound(sv_r,1):ubound(sv_r,1)));  r = sv_r
         allocate(r_edg(lbound(sv_redg,1):ubound(sv_redg,1)))
         r_edg = sv_redg
         allocate(dr_j(lbound(sv_dr,1):ubound(sv_dr,1)));  dr_j = sv_dr
         allocate(Gphi_i(lbound(sv_gpi,1):ubound(sv_gpi,1)))
         Gphi_i = sv_gpi
         allocate(Gphi_c(lbound(sv_gpc,1):ubound(sv_gpc,1)))
         Gphi_c = sv_gpc
         deallocate(sv_r, sv_redg, sv_dr, sv_gpi, sv_gpc)
      endif
      N = sv_N;  weno_mode = sv_wmode
      use_plm = sv_plm;  use_weno3 = sv_weno;  recon_lambda_on = sv_lam_on
      rec_method = sv_rec;  flux = sv_flux
      base_face_W = sv_bfW;  base_face_lower_W = sv_bflW

      end subroutine well_balanced_operator_rows

      ! ------------------------------------------------------!

      real*8 function key_departure_of_the_rows(xu, xn, xg, plm)        &
                                                              result(dev)
      ! The rows R = dF - S of one state evaluated twice, with the
      ! well-balanced key off and on, and the largest departure between
      ! them relative to the largest row of the base scheme, over the
      ! physical cells.  With no gravity the two are the same operator, so
      ! this is the number that says so.
      integer, intent(in) :: xn, xg
      real*8, dimension(3,1-xg:xn+xg), intent(in) :: xu
      logical, intent(in) :: plm
      real*8, dimension(3,1-xg:xn+xg) :: yWL, yWR, ydF, yS
      real*8, dimension(3,1-xg:xn+xg) :: R_off, R_on
      integer :: xj, xk
      real*8  :: sc

      if (plm) then
         rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
      else
         rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
      endif

      well_balanced = .false.
      call reconstruction_continuation_rhs(xu, yWL, yWR, ydF, yS)
      R_off = ydF - yS
      well_balanced = .true.
      call reconstruction_continuation_rhs(xu, yWL, yWR, ydF, yS)
      R_on = ydF - yS
      well_balanced = .false.

      dev = 0.0d0
      sc  = 0.0d0
      do xj = 1, xn
         do xk = 1, 3
            sc = max(sc, abs(R_off(xk,xj)))
         enddo
      enddo
      if (sc .le. 0.0d0) return
      do xj = 1, xn
         do xk = 1, 3
            dev = max(dev, abs(R_on(xk,xj) - R_off(xk,xj))/sc)
         enddo
      enddo

      end function key_departure_of_the_rows


      ! ------------------------------------------------------!

      subroutine generic_arm_face_departure_rows(nf)
      ! WHAT THE KIND-GENERIC ROWS HAND BACK BESIDE THE ROW (item P15).
      !
      ! Under the well-balanced key the momentum row's pressure-gradient
      ! term is the gradient of the FACE DEPARTURE q, (A+ q_up - A- q_dn)/dV
      ! under PLM and (q_up - q_dn)/dr under WENO3
      ! (momentum_row_terms_of_cell, RK_rhs.f90), and the departures are
      ! read from the module arrays face_q_up / face_q_dn of
      ! RK_integration.  The kind-generic rows are assembled here, so
      ! unless it stores its own departures in those arrays the term, the
      ! momentum row's reference scale built on it and the certification
      ! measure read against that scale describe the last state RK_rhs was
      ! evaluated on, or the zeros of the first allocation, and not the
      ! state they were assembled on.
      !
      ! Both groups below are evaluated on a state those rows assemble AFTER
      ! a right-hand side was evaluated on a DIFFERENT state, which is the
      ! configuration in which the two can be told apart:
      !
      !   * momentum_pressure_gradient and residual_row_scale(2, j, u) of
      !     those rows are those the production routines give for the SAME
      !     state.  Bitwise for the generic-double instantiation, which is
      !     the production operator bit for bit (tolerance 0, as in N37's
      !     generic_double_is_the_operator_bitwise_with_the_key), and to
      !     1e-13 for the quadruple one, which rounds to double once on the
      !     way out (the tolerance of N37's
      !     well_balanced_uniform_flow_is_the_base_scheme rows).
      !   * the stored departures are THOSE ROWS' OWN.  A first evaluation
      !     with no prior right-hand side, where the arrays hold zeros,
      !     cannot be arranged inside one driver process, so what is
      !     asserted instead is the identity that implies it: the momentum
      !     row they returned is rebuilt here from the stored face flux
      !     and the stored departures alone, and must be that row.  With
      !     departures belonging to another state the rebuild misses it.
      use RK_integration,      only: momentum_pressure_gradient,         &
                                     face_flux, face_q_up, face_q_dn
      use steady_residual_mod, only: assemble_residual,                  &
                                     residual_row_scale
      use ionization_equilibrium, only:                                  &
                                  set_ioniz_eq_sweep_state_kind,         &
                                  ieq_state_marching,                    &
                                  ieq_state_steady_iterate
      integer, intent(inout) :: nf
      integer, parameter :: xn = 12, xg = 2
      real*8, dimension(3,1-xg:xn+xg) :: uA, uB, Rp, Rq
      real*8, dimension(1-xg:xn+xg)   :: np1, zz
      real*8, dimension(1-xg:xn+xg)   :: pg_ref, sc_ref
      real*8  :: x_rho, x_vel, x_pre, x_b0
      real*8  :: d_pg, d_sc, d_row, ref_pg, ref_sc, ref_row, row_q
      real*8  :: x_dr, x_rp, x_rm, x_dAp, x_dAm, x_dV
      integer :: xj, xs, xa, x_status
      character(len=8)  :: sch
      character(len=10) :: kind_name
      integer :: sv_N, sv_wmode
      logical :: sv_plm, sv_weno, sv_lam_on, sv_wb, sv_visc, sv_cond
      real*8  :: sv_mu0
      character(len=:), allocatable :: sv_rec, sv_flux
      real*8, dimension(3) :: sv_bfW, sv_bflW
      real*8, dimension(:), allocatable :: sv_r, sv_redg, sv_dr,          &
                                           sv_gpi, sv_gpc

      sv_N = N;  sv_wmode = weno_mode;  sv_wb = well_balanced
      sv_plm = use_plm;  sv_weno = use_weno3;  sv_lam_on = recon_lambda_on
      sv_rec = rec_method;  sv_flux = flux
      sv_bfW = base_face_W;  sv_bflW = base_face_lower_W
      sv_visc = visc_on;  sv_cond = cond_on;  sv_mu0 = visc_mu0
      if (allocated(r)) then
         allocate(sv_r(lbound(r,1):ubound(r,1)));  sv_r = r
         allocate(sv_redg(lbound(r_edg,1):ubound(r_edg,1)));  sv_redg = r_edg
         allocate(sv_dr(lbound(dr_j,1):ubound(dr_j,1)));      sv_dr  = dr_j
         allocate(sv_gpi(lbound(Gphi_i,1):ubound(Gphi_i,1))); sv_gpi = Gphi_i
         allocate(sv_gpc(lbound(Gphi_c,1):ubound(Gphi_c,1))); sv_gpc = Gphi_c
         deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      endif

      ! The momentum row has to be the flux difference alone, so that the
      ! rebuild below states the row and not the row plus a split source:
      ! no operator-split viscosity, no conduction.
      visc_on = .false.;  cond_on = .false.;  visc_mu0 = 0.0d0
      N = xn;  weno_mode = 0;  recon_lambda_on = .false.;  flux = 'ROE'
      well_balanced = .true.
      allocate(r(1-xg:xn+xg), r_edg(1-xg:xn+xg), dr_j(1-xg:xn+xg))
      allocate(Gphi_i(1-xg:xn+xg), Gphi_c(1-xg:xn+xg))
      do xj = 1-xg, xn+xg
         r_edg(xj) = 1.0d0 + 0.05d0*dble(xj)*(1.0d0 + 0.02d0*dble(xj))
      enddo
      do xj = 2-xg, xn+xg
         r(xj)    = 0.5d0*(r_edg(xj) + r_edg(xj-1))
         dr_j(xj) = r_edg(xj) - r_edg(xj-1)
      enddo
      r(1-xg)    = r_edg(1-xg) - 0.5d0*dr_j(2-xg)
      dr_j(1-xg) = dr_j(2-xg)

      ! A gravitating column, so the departures are not all zero and the
      ! well-balanced substitution has something to cancel.
      x_b0 = 3.0d0
      do xj = 1-xg, xn+xg
         Gphi_i(xj) = -x_b0/r_edg(xj)
         Gphi_c(xj) = -x_b0/r(xj)
      enddo

      ! TWO STATES THAT SHARE NOTHING BUT THE GRID.  uB is the state the
      ! key is asked about; uA, with a different density slope, the
      ! opposite sign of velocity and a different pressure slope, is the
      ! state the previous right-hand side was evaluated on, so that its
      ! face departures are nowhere near uB's.
      do xj = 1-xg, xn+xg
         x_rho = 1.0d0/r(xj)**3
         x_vel = 0.3d0*r(xj)
         x_pre = 0.8d0/r(xj)**4
         uB(1,xj) = x_rho
         uB(2,xj) = x_rho*x_vel
         uB(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
         x_rho = 2.5d0/r(xj)**2
         x_vel = -0.4d0/r(xj)
         x_pre = 1.7d0/r(xj)**3
         uA(1,xj) = x_rho
         uA(2,xj) = x_rho*x_vel
         uA(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
      enddo
      np1 = 1.0d0
      zz  = 0.0d0

      do xs = 1, 2
         if (xs .eq. 1) then
            rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
            sch = '[PLM]'
         else
            rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
            sch = '[WENO3]'
         endif
         base_face_W(1) = uB(1,1)
         base_face_W(2) = uB(2,1)/uB(1,1)
         base_face_W(3) = (gamma_ad - 1.0d0)                              &
                          *(uB(3,1) - 0.5d0*uB(2,1)**2/uB(1,1))
         base_face_lower_W = base_face_W

         ! ---- the production routines on uB: the reference ----
         x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
         call set_ioniz_eq_sweep_state_kind(ieq_state_marching)
         call assemble_residual(uB, np1, zz, zz, Rp)
         ref_pg = 0.0d0;  ref_sc = 0.0d0
         do xj = 1, xn
            pg_ref(xj) = momentum_pressure_gradient(xj)
            sc_ref(xj) = residual_row_scale(2, xj, uB)
            ref_pg = max(ref_pg, abs(pg_ref(xj)))
            ref_sc = max(ref_sc, abs(sc_ref(xj)))
         enddo

         do xa = 1, 2
            ! ---- a right-hand side on the OTHER state ----
            x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
            call set_ioniz_eq_sweep_state_kind(ieq_state_marching)
            call assemble_residual(uA, np1, zz, zz, Rp)

            ! ---- the kind-generic rows on uB ----
            if (xa .eq. 1) then
               x_status = c_setenv('EXHALE_RESID_QUAD'//c_null_char,      &
                                   '1'//c_null_char, 1)
               kind_name = '[quad]'
            else
               x_status = c_setenv('EXHALE_RESID_QUAD'//c_null_char,      &
                                   '2'//c_null_char, 1)
               kind_name = '[double]'
            endif
            call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
            call assemble_residual(uB, np1, zz, zz, Rq)
            x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
            call set_ioniz_eq_sweep_state_kind(ieq_state_marching)

            d_pg = 0.0d0;  d_sc = 0.0d0
            do xj = 1, xn
               d_pg = max(d_pg,                                          &
                          abs(momentum_pressure_gradient(xj) - pg_ref(xj)))
               d_sc = max(d_sc,                                          &
                          abs(residual_row_scale(2, xj, uB) - sc_ref(xj)))
            enddo
            if (ref_pg .gt. 0.0d0) d_pg = d_pg/ref_pg
            if (ref_sc .gt. 0.0d0) d_sc = d_sc/ref_sc

            ! ---- the row rebuilt from the stored face data ----
            ! Under the key the source is zero and no transport is active,
            ! so the row Rq(2,j) IS the flux difference dF(2,j), which the
            ! assembly builds from the face flux and the departures alone.
            d_row = 0.0d0;  ref_row = 0.0d0
            do xj = 1, xn
               x_dr  = dr_j(xj)
               x_rp  = r_edg(xj)
               x_rm  = r_edg(xj-1)
               x_dAp = x_rp*x_rp
               x_dAm = x_rm*x_rm
               x_dV  = (x_dAp*x_rp - x_dAm*x_rm)/3.0d0
               if (use_plm) then
                  row_q = (x_dAp*(face_flux(2,xj)   + face_q_up(xj))      &
                         - x_dAm*(face_flux(2,xj-1)                       &
                                + face_q_dn(xj-1)))/x_dV
               else
                  row_q = (x_dAp*face_flux(2,xj)                          &
                         - x_dAm*face_flux(2,xj-1))/x_dV                  &
                        + (face_q_up(xj) - face_q_dn(xj-1))/x_dr
               endif
               d_row   = max(d_row, abs(row_q - Rq(2,xj)))
               ref_row = max(ref_row, abs(Rq(2,xj)))
            enddo
            if (ref_row .gt. 0.0d0) d_row = d_row/ref_row

            if (xa .eq. 1) then
               call rel_row('generic_rows_momentum_pressure_gradient_'//  &
                    'is_the_production_term'//trim(kind_name)//trim(sch), &
                    d_pg, 0.0d0, 1.0d-13, nf)
               call rel_row('generic_rows_momentum_row_scale_is_the_'//   &
                    'production_scale'//trim(kind_name)//trim(sch),        &
                    d_sc, 0.0d0, 1.0d-13, nf)
            else
               call rel_row('generic_rows_momentum_pressure_gradient_'//  &
                    'is_the_production_term'//trim(kind_name)//trim(sch), &
                    d_pg, 0.0d0, 0.0d0, nf)
               call rel_row('generic_rows_momentum_row_scale_is_the_'//   &
                    'production_scale'//trim(kind_name)//trim(sch),        &
                    d_sc, 0.0d0, 0.0d0, nf)
            endif
            call rel_row('generic_rows_stored_departures_rebuild_the_'//  &
                 'momentum_row'//trim(kind_name)//trim(sch),               &
                 d_row, 0.0d0, 1.0d-13, nf)
         enddo
      enddo

      ! ---- put the globals back ----
      x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
      call set_ioniz_eq_sweep_state_kind(ieq_state_marching)
      well_balanced = sv_wb
      visc_on = sv_visc;  cond_on = sv_cond;  visc_mu0 = sv_mu0
      deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      if (allocated(sv_r)) then
         allocate(r(lbound(sv_r,1):ubound(sv_r,1)));  r = sv_r
         allocate(r_edg(lbound(sv_redg,1):ubound(sv_redg,1)))
         r_edg = sv_redg
         allocate(dr_j(lbound(sv_dr,1):ubound(sv_dr,1)));  dr_j = sv_dr
         allocate(Gphi_i(lbound(sv_gpi,1):ubound(sv_gpi,1)))
         Gphi_i = sv_gpi
         allocate(Gphi_c(lbound(sv_gpc,1):ubound(sv_gpc,1)))
         Gphi_c = sv_gpc
         deallocate(sv_r, sv_redg, sv_dr, sv_gpi, sv_gpc)
      endif
      N = sv_N;  weno_mode = sv_wmode
      use_plm = sv_plm;  use_weno3 = sv_weno;  recon_lambda_on = sv_lam_on
      rec_method = sv_rec;  flux = sv_flux
      base_face_W = sv_bfW;  base_face_lower_W = sv_bflW

      end subroutine generic_arm_face_departure_rows

      ! ------------------------------------------------------!

      subroutine extended_precision_rows(nf)
      ! THE HYDRODYNAMIC ROWS EVALUATED IN ANOTHER REAL KIND (item N34).
      ! The stationary residual's non-smoothness floor is the rounding of
      ! the flux assembly (N33), and the control experiment that separates
      ! the rounding from every other candidate is the SAME operator in a
      ! wider mantissa. That is only an experiment about the arithmetic if
      ! the two instantiations of the generic text are one operator, so
      ! these rows state:
      !
      !   * the selection is off unless the environment names it, and it
      !     names the quadruple and the generic-double instantiation
      !     separately;
      !   * the DOUBLE instantiation of the generic text reproduces
      !     Reconstruct + RK_rhs + Num_flux + source BIT FOR BIT on a state
      !     with structure in it (a bitwise reference, not a tolerance);
      !   * the QUADRUPLE instantiation gets the analytically known answer
      !     of a case where one exists. On a state uniform in (rho, v, p)
      !     with no gravity every interface flux is the physical flux of
      !     that state (the Roe jumps vanish identically), so the mass row
      !     is exactly F_rho (r_+^2 - r_-^2)/dV. That reference is formed
      !     here in quadruple from the same grid, and the row asserts the
      !     quadruple pipeline reaches it to 1e-30 while the double one
      !     stands at its own rounding.
      integer, intent(inout) :: nf
      integer, parameter :: qp = selected_real_kind(30, 300)
      integer, parameter :: xn = 12, xg = 2
      real*8, dimension(3,1-xg:xn+xg) :: xu, xWLa, xWRa, xdFa, xSa
      real*8, dimension(3,1-xg:xn+xg) :: xWLb, xWRb, xdFb, xSb
      real*8, dimension(3,1-xg:xn+xg) :: xff
      real*8, dimension(1-xg:xn+xg)   :: xfpr
      real(qp) :: q_rp, q_rm, q_dV, q_rho, q_vel, q_ref
      real(qp), dimension(3,1-xg:xn+xg) :: q_dF
      real*8  :: x_rho, x_vel, x_pre, x_worst_dF, x_worst_S, x_dev
      real*8  :: x_dev_q, x_dev_d
      integer :: xj, xk, x_scaled, x_hlle, x_llf, x_status
      ! Saved global state; every row here rewrites the grid and the
      ! scheme selection and puts them back.
      integer :: sv_N, sv_wmode
      logical :: sv_plm, sv_weno, sv_lam_on
      character(len=:), allocatable :: sv_rec, sv_flux
      real*8, dimension(3) :: sv_bfW, sv_bflW
      real*8, dimension(:), allocatable :: sv_r, sv_redg, sv_dr,          &
                                           sv_gpi, sv_gpc

      ! ---- the selection, and only from the environment ----
      x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
      call int_row('the_extended_precision_rows_are_off_by_default',      &
                   generic_precision_rows_selected(), ROWS_PRODUCTION, nf)
      x_status = c_setenv('EXHALE_RESID_QUAD'//c_null_char,               &
                          '1'//c_null_char, 1)
      call int_row('the_selector_names_the_quadruple_instantiation',      &
                   generic_precision_rows_selected(), ROWS_QUADRUPLE, nf)
      x_status = c_setenv('EXHALE_RESID_QUAD'//c_null_char,               &
                          '2'//c_null_char, 1)
      call int_row('the_selector_names_the_generic_double_'//             &
                   'instantiation',                                      &
                   generic_precision_rows_selected(), ROWS_GENERIC_DOUBLE, nf)
      x_status = c_setenv('EXHALE_RESID_QUAD'//c_null_char,               &
                          'yes'//c_null_char, 1)
      call int_row('an_unnamed_value_leaves_the_production_rows',         &
                   generic_precision_rows_selected(), ROWS_PRODUCTION, nf)
      x_status = c_unsetenv('EXHALE_RESID_QUAD'//c_null_char)
      call int_row('the_selection_is_dropped_again',                      &
                   generic_precision_rows_selected(), ROWS_PRODUCTION, nf)

      ! ---- a grid and a scheme to evaluate the operator on ----
      sv_N = N;  sv_wmode = weno_mode
      sv_plm = use_plm;  sv_weno = use_weno3;  sv_lam_on = recon_lambda_on
      sv_rec = rec_method;  sv_flux = flux
      sv_bfW = base_face_W;  sv_bflW = base_face_lower_W
      if (allocated(r)) then
         allocate(sv_r(lbound(r,1):ubound(r,1)));  sv_r = r
         allocate(sv_redg(lbound(r_edg,1):ubound(r_edg,1)));  sv_redg = r_edg
         allocate(sv_dr(lbound(dr_j,1):ubound(dr_j,1)));      sv_dr  = dr_j
         allocate(sv_gpi(lbound(Gphi_i,1):ubound(Gphi_i,1))); sv_gpi = Gphi_i
         allocate(sv_gpc(lbound(Gphi_c,1):ubound(Gphi_c,1))); sv_gpc = Gphi_c
         deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      endif

      N = xn;  weno_mode = 0
      use_plm = .false.;  use_weno3 = .true.;  recon_lambda_on = .false.
      rec_method = 'WENO3';  flux = 'ROE'
      allocate(r(1-xg:xn+xg), r_edg(1-xg:xn+xg), dr_j(1-xg:xn+xg))
      allocate(Gphi_i(1-xg:xn+xg), Gphi_c(1-xg:xn+xg))
      ! A stretched radial grid, as the production grid is; the faces and
      ! the centres are consistent with each other so the volume weights of
      ! the reconstruction are those of a real grid.
      do xj = 1-xg, xn+xg
         r_edg(xj) = 1.0d0 + 0.05d0*dble(xj)*(1.0d0 + 0.02d0*dble(xj))
      enddo
      do xj = 2-xg, xn+xg
         r(xj)    = 0.5d0*(r_edg(xj) + r_edg(xj-1))
         dr_j(xj) = r_edg(xj) - r_edg(xj-1)
      enddo
      r(1-xg)    = r_edg(1-xg) - 0.5d0*dr_j(2-xg)
      dr_j(1-xg) = dr_j(2-xg)
      Gphi_i = 0.0d0
      Gphi_c = 0.0d0

      ! ---- the double instantiation IS the production operator ----
      ! A state with structure in it: a density and a pressure that fall
      ! with radius and a velocity that grows, so the reconstruction
      ! limits, the Roe jumps are nonzero and every branch of the assembly
      ! is reached with something to do.
      do xj = 1-xg, xn+xg
         x_rho = 1.0d0/r(xj)**3
         x_vel = 0.3d0*r(xj)
         x_pre = 0.8d0/r(xj)**4
         xu(1,xj) = x_rho
         xu(2,xj) = x_rho*x_vel
         xu(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = 1.0d0/r(1)**3
      base_face_W(2) = 0.3d0*r(1)
      base_face_W(3) = 0.8d0/r(1)**4
      base_face_lower_W = base_face_W

      call reconstruction_continuation_rhs(xu, xWLa, xWRa, xdFa, xSa)
      call hydrodynamic_rows_in_double_precision(xu, xWLb, xWRb, xdFb, xSb)

      ! Bitwise: the reference is zero and the tolerance is zero, so any
      ! difference at all in any component of any cell fails the row.
      x_worst_dF = 0.0d0;  x_worst_S = 0.0d0
      do xj = 2-xg, xn+xg
         do xk = 1, 3
            if (xdFb(xk,xj) .ne. xdFa(xk,xj))                             &
               x_worst_dF = max(x_worst_dF, abs(xdFb(xk,xj)-xdFa(xk,xj)))
            if (xSb(xk,xj) .ne. xSa(xk,xj))                               &
               x_worst_S  = max(x_worst_S,  abs(xSb(xk,xj)-xSa(xk,xj)))
         enddo
      enddo
      call rel_row('generic_double_flux_difference_is_the_operator_bitwise',&
                   x_worst_dF, 0.0d0, 0.0d0, nf)
      call rel_row('generic_double_source_is_the_operator_bitwise',       &
                   x_worst_S, 0.0d0, 0.0d0, nf)
      x_dev = 0.0d0
      do xj = 1-xg, xn+xg
         do xk = 1, 3
            if (xWLb(xk,xj) .ne. xWLa(xk,xj))                             &
               x_dev = max(x_dev, abs(xWLb(xk,xj)-xWLa(xk,xj)))
            if (xWRb(xk,xj) .ne. xWRa(xk,xj))                             &
               x_dev = max(x_dev, abs(xWRb(xk,xj)-xWRa(xk,xj)))
         enddo
      enddo
      call rel_row('generic_double_face_states_are_the_operator_bitwise', &
                   x_dev, 0.0d0, 0.0d0, nf)

      ! ---- the quadruple instantiation on a case with an exact answer ----
      ! Uniform (rho, v, p) and no gravity: every reconstructed face state
      ! is that state, the Roe jumps vanish, and the interface flux is the
      ! physical flux. The mass row is then exactly F_rho (r_+^2 - r_-^2)/dV,
      ! with dV = (r_+^3 - r_-^3)/3.
      ! rho and v are exact binary fractions, so the mass flux
      ! rho (rho v)/rho of every face is the exact 0.625 in any kind and
      ! the reference below is the analytic row and not a re-rounding of
      ! the pipeline's own arithmetic.
      x_rho = 1.25d0;  x_vel = 0.5d0;  x_pre = 0.75d0
      do xj = 1-xg, xn+xg
         xu(1,xj) = x_rho
         xu(2,xj) = x_rho*x_vel
         xu(3,xj) = 0.5d0*x_rho*x_vel*x_vel + x_pre/(gamma_ad - 1.0d0)
      enddo
      base_face_W(1) = x_rho
      base_face_W(2) = x_vel
      base_face_W(3) = x_pre
      base_face_lower_W = base_face_W

      call quadruple_rows(xu, xWLb, xWRb, xdFb, xSb, xff, xfpr,           &
                          x_scaled, x_hlle, x_llf, q_dF)
      call hydrodynamic_rows_in_double_precision(xu, xWLa, xWRa, xdFa, xSa)

      q_rho = real(x_rho, qp)
      q_vel = real(x_vel, qp)
      ! Cell 1 is left out: its lower face takes the imposed base state,
      ! whose pressure is the double round of the same map and not the
      ! quadruple one, so that one face carries a jump of the size of a
      ! double ulp by construction. Every other cell of the domain sees
      ! the uniform state on both sides.
      x_dev_q = 0.0d0;  x_dev_d = 0.0d0
      do xj = 2, xn
         q_rp = real(r_edg(xj),   qp)
         q_rm = real(r_edg(xj-1), qp)
         q_dV = (q_rp*q_rp*q_rp - q_rm*q_rm*q_rm)/3.0_qp
         q_ref = q_rho*q_vel*(q_rp*q_rp - q_rm*q_rm)/q_dV
         ! The quadruple row BEFORE the rounding to double: the single
         ! rounding on the way out is what these rows hand the solver, but
         ! it is also a double ulp, which would hide the arithmetic this
         ! row is about.
         x_dev_q = max(x_dev_q,                                           &
                   real(abs(q_dF(1,xj) - q_ref)/abs(q_ref), 8))
         x_dev_d = max(x_dev_d,                                           &
                   real(abs(real(xdFa(1,xj), qp) - q_ref)/abs(q_ref), 8))
      enddo
      call rel_row('quadruple_mass_row_of_a_uniform_state_is_exact',      &
                   x_dev_q, 0.0d0, 1.0d-30, nf)
      ! And that the reference is a real test of the arithmetic: the same
      ! row in double stands at its own rounding, far above 1e-30 and far
      ! below one.
      call log_row('the_same_row_in_double_stands_at_its_own_rounding',   &
                   (x_dev_d .gt. 1.0d-30) .and. (x_dev_d .lt. 1.0d-12),   &
                   .true., nf)
      call int_row('a_uniform_state_needs_no_hlle_fallback',              &
                   x_hlle, 0, nf)
      call int_row('a_uniform_state_needs_no_positivity_scaling',         &
                   x_scaled, 0, nf)

      ! ---- put the globals back ----
      deallocate(r, r_edg, dr_j, Gphi_i, Gphi_c)
      if (allocated(sv_r)) then
         allocate(r(lbound(sv_r,1):ubound(sv_r,1)));  r = sv_r
         allocate(r_edg(lbound(sv_redg,1):ubound(sv_redg,1)))
         r_edg = sv_redg
         allocate(dr_j(lbound(sv_dr,1):ubound(sv_dr,1)));  dr_j = sv_dr
         allocate(Gphi_i(lbound(sv_gpi,1):ubound(sv_gpi,1)))
         Gphi_i = sv_gpi
         allocate(Gphi_c(lbound(sv_gpc,1):ubound(sv_gpc,1)))
         Gphi_c = sv_gpc
         deallocate(sv_r, sv_redg, sv_dr, sv_gpi, sv_gpc)
      endif
      N = sv_N;  weno_mode = sv_wmode
      use_plm = sv_plm;  use_weno3 = sv_weno;  recon_lambda_on = sv_lam_on
      rec_method = sv_rec;  flux = sv_flux
      base_face_W = sv_bfW;  base_face_lower_W = sv_bflW

      end subroutine extended_precision_rows

      end program krylov_and_dogleg_tests

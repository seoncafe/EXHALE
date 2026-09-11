! Diagnostic derived from src/tests/carrier_retry/carrier_retry.f90 on 2026-09-12.
! Production routines are unchanged. Existing fixture statements are retained;
! original full-line comments are omitted. Added AUDIT rows test returned-state contracts.
      program carrier_state_contract_probe

      use global_parameters
      use species_table
      use ionization_equilibrium, only: bg_cell, bg_ready,                &
                                        ioniz_eq_allocate_arrays, ioniz_eq
      use diffusive_photochemistry, only:                                 &
           carrier_set_init, carrier_transport_interval,                  &
           carrier_checkpoint, carrier_checkpoint_take,                   &
           carrier_checkpoint_restore, carrier_checkpoint_matches,        &
           carrier_perturb_checkpointed_state_for_test,                   &
           carrier_reject_leading_attempts_for_test,                      &
           carrier_initial_substep_for_test,                              &
           carrier_attempt_record, carrier_co_domain_record,              &
           carrier_co_domain_perturb_for_test, carrier_retry_max,         &
           carrier_verdict, n_carrier_max, ic_H2, ic_Hp,                  &
           carrier_row_scales, carrier_returned_state_verdict,            &
           carrier_solve_converged, carrier_stall_at_floor,             &
           carrier_transport_diagnostics,                               &
           carrier_last_interval_halvings, carrier_row_roundoff,        &
           carrier_reject_row_residual, carrier_refused_rows,           &
           carrier_refused_row, carrier_roundoff_limited_record,         &
           carrier_record_refused_rows, carrier_name,                     &
           photochemical_transport_step, carrier_interval_covered,        &
           carrier_interval_exhausted, carrier_last_interval_status,      &
           carrier_exhausted_record, carrier_history_certifiable,         &
           carrier_transport_stop_on_failure,                             &
           carrier_transport_stops_suppressed,                            &
           carrier_step_verdict, carrier_last_verdict_of_run,             &
           carrier_no_interval,                                           &
           relax_photochemical_composition,                               &
           carrier_relax_movement_bound, carrier_relax_fixed_point,       &
           carrier_relax_step_budget, carrier_relax_interval_refused,     &
           carrier_relax_nothing_to_advance, carrier_relax_outcome_text
      use element_census, only: element_nuclei_and_charge, n_element
      use gravity_grid_construction, only: set_gravity_grid
      use energy_vectors_construct,  only: set_energy_vectors
      use charge_exchange,           only: cx_init
      use composition,               only: get_species_densities,         &
                                           comp_p_from_T, comp_T_from_p
      use base_boundary,             only: set_base_reservoir
      use BC_Apply,                  only: Apply_BC
      use steady_residual_mod,       only: assemble_residual
      use caloric_eos, only: energy_density_from_pressure, pressure_from_energy_density
      use certification, only: mass_row_cell_verdict
      use assertion_report
      implicit none

      integer, parameter :: ncell = 12
      real*8,  allocatable :: rho(:), v(:), dt_code(:)
      real*8,  allocatable :: f_sp(:,:), f_sp0(:,:)
      real*8,  allocatable :: f_one(:,:), f_half(:,:), f_quarter(:,:)
      real*8,  allocatable :: f_retried(:,:)
      real*8,  allocatable :: fc_probe(:,:)
      type(carrier_checkpoint) :: chk0, chk_fresh, chk_after
      type(carrier_verdict)    :: verdict
      logical :: completed, ok
      real*8  :: frac_done, e1, e2
      integer :: nsub
      integer :: nattempt0, nreject0, naccept0, nretried0
      integer :: nattempt,  nreject,  naccept,  nretried
      integer :: nout_co0, nhot_co0, nhep_co0
      integer :: nout_co,  nhot_co,  nhep_co
      real*8,  allocatable :: tfull_dt(:,:), tphys_dt(:,:)
      real*8,  allocatable :: tfull_short(:,:), tphys_short(:,:)
      real*8,  allocatable :: res_dt(:,:), res_short(:,:)
      real*8,  allocatable :: dt_short(:)
      logical, allocatable :: unconstrained(:)
      type(carrier_verdict) :: vd_full_old, vd_full_new
      type(carrier_verdict) :: vd_short_old, vd_short_new
      real*8  :: ratio_co0, r_co0, form_co0
      real*8  :: ratio_co,  r_co,  form_co
      real*8,  allocatable :: dt_1024(:)
      real*8  :: rw_dt, rw_16, rw_1024
      real*8  :: rs_dt, rs_16, rs_1024, worst_lim
      integer :: nit_dt, nit_16, nit_1024, nlim
      integer :: nsub_dt, nsub_16, nsub_1024
      integer :: nhalve_1024, nrl_1024, nrl_last
      integer :: nrl_total, nrl_substeps
      integer :: nrefused, jref, icref
      real*8  :: res_ref, phys_ref, full_ref
      real*8,  allocatable :: tf_round(:,:)
      type(carrier_verdict) :: vd_round, vd_mixed
      logical :: cov_dt, cov_16, cov_1024
      integer :: st_phys, st_init, nstop0, nstop_phys, nstop_init
      integer :: nex0, nex_phys, nex_init
      integer :: ex_first, ex_last, ex_j, ex_ic
      real*8  :: ex_ratio, ex_phys, ex_full
      integer :: dout_tot, dhot_tot, dhep_tot
      integer :: dout_f1,  dhot_f1,  dhep_f1
      integer :: dout_f2,  dhot_f2,  dhep_f2
      real*8  :: drat_tot, dr_tot, dform_tot
      real*8  :: drat_f1,  dr_f1,  dform_f1
      real*8  :: drat_f2,  dr_f2,  dform_f2
      integer :: st_nothing, reason_before
      type(carrier_verdict) :: vd_step, vd_run
      real*8  :: dr_wide, dr_tight, dr_none, dr_zero, dr_ref
      integer :: ns_wide, ns_tight, ns_none, ns_zero, ns_ref
      integer :: oc_wide, oc_tight, oc_none, oc_zero, oc_ref
      real*8, allocatable :: nnuc0(:,:), nnuc1(:,:)
      real*8, allocatable :: nchg0(:), nchg1(:), rcomp0(:), rcomp1(:)
      real*8, allocatable :: v_relax(:), v_fast(:)
      real*8, allocatable :: p_col(:), T_col(:)
      real*8, allocatable :: heat_col(:), cool_col(:), eta_col(:)
      real*8, allocatable :: f_swept(:,:)
      real*8, allocatable :: nhi_r(:), nhii_r(:), nhei_r(:), nheii_r(:)
      real*8, allocatable :: nheiii_r(:), nheiTR_r(:), ne_r(:), ntot_r(:)
      real*8, allocatable :: nm_r(:,:), T_ret(:)
      real*8, allocatable :: ntot0_c(:), T_ent(:)
      real*8  :: move_pass, move_sweep, dev_ntot, dev_TK
      integer :: jj
      real*8, allocatable :: audit_energy(:)
      real*8 :: audit_tol, audit_dist, audit_p_error, audit_energy_error
      logical :: audit_within, audit_anchored

      call setup_globals()
      call build_molecular_hydrogen_column()

      allocate(fc_probe(1-Ng:N+Ng,n_carrier_max))
      fc_probe = 0.25d0
      call carrier_checkpoint_take(chk_fresh, fc_probe)
      call carrier_perturb_checkpointed_state_for_test(fc_probe)
      call check_absolute('a_perturbed_module_is_not_the_checkpoint',     &
           logical_as_double(carrier_checkpoint_matches(chk_fresh,        &
                                                        fc_probe)),       &
           0.0d0, 0.0d0)
      call carrier_checkpoint_restore(chk_fresh, fc_probe)
      call check_absolute('the_unallocated_module_is_restored_exactly',   &
           logical_as_double(carrier_checkpoint_matches(chk_fresh,        &
                                                        fc_probe)),       &
           1.0d0, 0.0d0)
      call check_absolute('the_carrier_fractions_are_restored_exactly',   &
           maxval(abs(fc_probe - 0.25d0)), 0.0d0, 0.0d0)

      f_sp = f_sp0
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      call check_absolute('the_reference_interval_is_covered',            &
           logical_as_double(completed), 1.0d0, 0.0d0)
      call check_absolute('the_reference_interval_takes_one_substep',     &
           dble(nsub), 1.0d0, 0.0d0)
      call check_absolute('the_carriers_moved_over_the_interval',         &
           logical_as_double(maxval(abs(f_sp(1:N,isp_H2)                  &
                                        - f_sp0(1:N,isp_H2)))            &
                             .gt. 1.0d-12*f_sp0(1,isp_H2)), 1.0d0, 0.0d0)
      f_one = f_sp

      call carrier_checkpoint_take(chk_after, fc_probe)
      call carrier_perturb_checkpointed_state_for_test(fc_probe)
      call check_absolute('a_perturbed_filled_module_differs',            &
           logical_as_double(carrier_checkpoint_matches(chk_after,        &
                                                        fc_probe)),       &
           0.0d0, 0.0d0)
      call carrier_checkpoint_restore(chk_after, fc_probe)
      call check_absolute('the_filled_module_is_restored_exactly',        &
           logical_as_double(carrier_checkpoint_matches(chk_after,        &
                                                        fc_probe)),       &
           1.0d0, 0.0d0)

      call carrier_checkpoint_take(chk0)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      carrier_reject_leading_attempts_for_test = 1
      call carrier_attempt_record(nattempt0, nreject0, naccept0, nretried0)
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      carrier_reject_leading_attempts_for_test = 0
      call carrier_attempt_record(nattempt, nreject, naccept, nretried)
      call check_absolute('a_refused_first_attempt_still_covers_the'//    &
           '_interval', logical_as_double(completed), 1.0d0, 0.0d0)
      call check_absolute('the_elapsed_interval_is_dt_exactly',           &
           frac_done, 1.0d0, 0.0d0)
      call check_absolute('one_refusal_costs_one_refused_and_two'//       &
           '_accepted_substeps', dble(nsub), 3.0d0, 0.0d0)
      call check_absolute('the_attempted_record_counts_the_refusal',      &
           dble(nreject - nreject0), 1.0d0, 0.0d0)
      call check_absolute('the_accepted_record_counts_the_two_halves',    &
           dble(naccept - naccept0), 2.0d0, 0.0d0)
      call check_absolute('the_interval_is_recorded_as_retried_once',     &
           dble(nretried - nretried0), 1.0d0, 0.0d0)
      f_retried = f_sp

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      carrier_initial_substep_for_test = 0.5d0
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      carrier_initial_substep_for_test = 1.0d0
      call check_absolute('the_subdivided_interval_is_covered',           &
           logical_as_double(completed), 1.0d0, 0.0d0)
      call check_absolute('the_subdivided_interval_takes_two_substeps',   &
           dble(nsub), 2.0d0, 0.0d0)
      f_half = f_sp
      call check_relative('the_retried_state_is_the_subdivided'//         &
           '_integration_H2', maxval(abs(f_retried(1:N,isp_H2))),         &
           maxval(abs(f_half(1:N,isp_H2))), 1.0d-10)
      call check_absolute('the_retried_state_is_the_subdivided'//         &
           '_integration_every_row',                                      &
           carrier_row_difference(f_retried, f_half), 0.0d0, 1.0d-10)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_attempt_record(nattempt0, nreject0, naccept0, nretried0)
      call carrier_co_domain_perturb_for_test(3, 5)
      call carrier_co_domain_record(nout_co0, nhot_co0, nhep_co0,         &
                                    ratio_co0, r_co0, form_co0)
      carrier_reject_leading_attempts_for_test = 10000
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      carrier_reject_leading_attempts_for_test = 0
      call carrier_attempt_record(nattempt, nreject, naccept, nretried)
      call carrier_co_domain_record(nout_co, nhot_co, nhep_co,            &
                                    ratio_co, r_co, form_co)
      call check_absolute('an_exhausted_interval_is_not_covered',         &
           logical_as_double(completed), 0.0d0, 0.0d0)
      call check_absolute('an_exhausted_interval_accepted_no_substep'//  &
           '_before_the_failure', frac_done, 0.0d0, 0.0d0)
      call check_absolute('exhaustion_costs_one_attempt_per_halving'//    &
           '_plus_the_first', dble(nsub), dble(carrier_retry_max + 1),    &
           0.0d0)
      call check_absolute('the_composition_is_not_written',               &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      ok = carrier_checkpoint_matches(chk0)
      call check_absolute('the_entry_module_state_is_restored_exactly',   &
           logical_as_double(ok), 1.0d0, 0.0d0)
      call check_absolute('no_accepted_substep_is_recorded',              &
           dble(naccept - naccept0), 0.0d0, 0.0d0)
      call check_absolute('every_attempt_is_recorded_as_refused',         &
           dble(nreject - nreject0), dble(carrier_retry_max + 1), 0.0d0)
      call check_absolute('the_CO_domain_record_survives_the_restore',    &
           dble(nout_co - nout_co0), 0.0d0, 0.0d0)
      call check_absolute('and_so_does_its_hot_cell_count',               &
           dble(nhot_co - nhot_co0), 0.0d0, 0.0d0)
      call check_absolute('and_its_worst_ratio',                          &
           ratio_co - ratio_co0, 0.0d0, 0.0d0)
      call check_positive('the_domain_record_read_here_is_not_zero',      &
           dble(nout_co0))

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      carrier_initial_substep_for_test = 0.25d0
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      carrier_initial_substep_for_test = 1.0d0
      call check_absolute('the_quarter_interval_takes_four_substeps',     &
           dble(nsub), 4.0d0, 0.0d0)
      f_quarter = f_sp

      e1 = carrier_row_difference(f_one,  f_half)
      e2 = carrier_row_difference(f_half, f_quarter)
      call check_positive('step_doubling_difference_dt_against_two'//     &
                          '_halves', e1)
      call check_positive('step_doubling_difference_two_halves'//         &
                          '_against_four_quarters', e2)
      call check_absolute('halving_again_reduces_the_step_doubling'//     &
           '_difference', logical_as_double(e2 .lt. e1), 1.0d0, 0.0d0)
      write(*,'(a,es12.4,a,es12.4,a,f8.3)')                               &
           ' (carrier_retry) step-doubling difference: dt vs 2 x dt/2 = ', &
           e1, ', 2 x dt/2 vs 4 x dt/4 = ', e2, ', ratio = ', e1/max(e2,  &
           1.0d-300)

      allocate(tfull_dt(1:N,n_carrier_max), tphys_dt(1:N,n_carrier_max))
      allocate(tfull_short(1:N,n_carrier_max),                            &
               tphys_short(1:N,n_carrier_max))
      allocate(res_dt(1:N,n_carrier_max), res_short(1:N,n_carrier_max))
      allocate(unconstrained(1:N), dt_short(1-Ng:N+Ng))
      unconstrained = .false.

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      call carrier_row_scales(tfull_dt, tphys_dt)

      f_sp     = f_sp0
      dt_short = dt_code/16.0d0
      call carrier_checkpoint_restore(chk0)
      call carrier_transport_interval(rho, v, f_sp, dt_short, completed,  &
                                      frac_done, nsub, verdict)
      call carrier_row_scales(tfull_short, tphys_short)

      call check_absolute('the_short_step_interval_is_covered',           &
           logical_as_double(completed), 1.0d0, 0.0d0)
      call check_absolute('the_full_term_scale_grows_when_the_step_is'//  &
           '_shortened', logical_as_double(maxval(tfull_short(1:N,ic_H2)) &
           .gt. 4.0d0*maxval(tfull_dt(1:N,ic_H2))), 1.0d0, 0.0d0)

      res_dt    = 1.0d-7*tphys_dt
      res_short = 1.0d-7*tphys_short

      call carrier_returned_state_verdict(res_dt, tfull_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_full_old)
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_full_new)
      call carrier_returned_state_verdict(res_short, tfull_short,         &
           unconstrained, carrier_solve_converged, carrier_stall_at_floor,&
           vd_short_old)
      call carrier_returned_state_verdict(res_short, tphys_short,         &
           unconstrained, carrier_solve_converged, carrier_stall_at_floor,&
           vd_short_new)

      call check_absolute('the_full_term_scale_refuses_the_imbalance'//   &
           '_at_the_full_step',                                           &
           logical_as_double(vd_full_old%accepted), 0.0d0, 0.0d0)
      call check_absolute('the_full_term_scale_accepts_the_imbalance'//   &
           '_at_a_short_step',                                            &
           logical_as_double(vd_short_old%accepted), 1.0d0, 0.0d0)
      call check_absolute('the_full_term_measure_of_the_same_defect'//    &
           '_falls_with_the_step',                                        &
           logical_as_double(vd_short_old%rworst .lt.                     &
                             0.25d0*vd_full_old%rworst), 1.0d0, 0.0d0)
      call check_absolute('the_physical_scale_refuses_it_at_the_full'//   &
           '_step', logical_as_double(vd_full_new%accepted), 0.0d0, 0.0d0)
      call check_absolute('the_physical_scale_refuses_it_at_the_short'//  &
           '_step', logical_as_double(vd_short_new%accepted), 0.0d0,      &
           0.0d0)
      call check_relative('the_physical_measure_is_the_same_at_both'//    &
           '_steps', vd_short_new%rworst, vd_full_new%rworst, 1.0d-12)
      call check_relative('the_physical_measure_is_the_imbalance'//       &
           '_itself', vd_full_new%rworst, 1.0d-7, 1.0d-12)
      write(*,'(a,es12.4,a,es12.4)')                                      &
           ' (carrier_retry) a 1e-7 physical imbalance measures, on the'// &
           ' full-term scale, ', vd_full_old%rworst, ' at dt and ',       &
           vd_short_old%rworst, ' at dt/16'

      allocate(dt_1024(1-Ng:N+Ng))
      dt_1024 = dt_code/1024.0d0

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub_dt, verdict)
      cov_dt = completed
      rw_dt  = verdict%rworst
      call carrier_transport_diagnostics(nit_dt, rs_dt, nlim, worst_lim)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_transport_interval(rho, v, f_sp, dt_short, completed,  &
                                      frac_done, nsub_16, verdict)
      cov_16 = completed
      rw_16  = verdict%rworst
      call carrier_transport_diagnostics(nit_16, rs_16, nlim, worst_lim)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_transport_interval(rho, v, f_sp, dt_1024, completed,   &
                                      frac_done, nsub_1024, verdict)
      cov_1024    = completed
      rw_1024     = verdict%rworst
      nrl_1024    = verdict%n_rows_roundoff_limited
      nhalve_1024 = carrier_last_interval_halvings()
      call carrier_roundoff_limited_record(nrl_last, nrl_total,           &
                                           nrl_substeps)
      call carrier_transport_diagnostics(nit_1024, rs_1024, nlim,         &
                                         worst_lim)

      write(*,'(a)') ' (carrier_retry) substep, covered, worst physical'//&
           ' row of the returned state, full-scale Newton residual,'//    &
           ' iterations, attempts:'
      write(*,'(a,l2,es12.4,es12.4,i5,i4)') '   dt      ', cov_dt,        &
           rw_dt, rs_dt, nit_dt, nsub_dt
      write(*,'(a,l2,es12.4,es12.4,i5,i4)') '   dt/16   ', cov_16,        &
           rw_16, rs_16, nit_16, nsub_16
      write(*,'(a,l2,es12.4,es12.4,i5,i4)') '   dt/1024 ', cov_1024,      &
           rw_1024, rs_1024, nit_1024, nsub_1024
      write(*,'(a,i0,a,i0,a,i0)') ' (carrier_retry) rows accepted as'//   &
           ' round-off limited at dt/1024: ', nrl_1024,                   &
           '; over the run so far ', nrl_total, ' in substeps ',          &
           nrl_substeps

      call check_absolute('the_unit_column_is_covered_at_dt',             &
           logical_as_double(cov_dt), 1.0d0, 0.0d0)
      call check_absolute('the_unit_column_is_covered_at_dt_over_16',     &
           logical_as_double(cov_16), 1.0d0, 0.0d0)
      call check_absolute('the_short_substep_is_covered',                 &
           logical_as_double(cov_1024), 1.0d0, 0.0d0)
      call check_absolute('the_short_substep_costs_no_halving',           &
           logical_as_double(nhalve_1024 .eq. 0 .and. nsub_1024 .eq. 1),  &
           1.0d0, 0.0d0)
      call check_absolute('the_short_substep_reports_round_off_limited'// &
           '_rows', logical_as_double(nrl_1024 .ge. 1), 1.0d0, 0.0d0)
      call check_absolute('the_run_record_carries_those_rows',            &
           logical_as_double(nrl_last .eq. nrl_1024 .and.                 &
                             nrl_total .ge. nrl_1024 .and.                &
                             nrl_substeps .ge. 1), 1.0d0, 0.0d0)
      call check_absolute('the_short_substep_state_is_at_the_round_off'// &
           '_of_the_full_terms',                                          &
           logical_as_double(rs_1024 .le. carrier_row_roundoff), 1.0d0,   &
           0.0d0)
      call check_absolute('a_covered_interval_reports_no_refused_row',    &
           dble(carrier_refused_rows()), 0.0d0, 0.0d0)
      call check_absolute('the_covered_steps_return_an_acceptable'//      &
           '_physical_measure',                                           &
           logical_as_double(rw_dt .le. 1.0d-8 .and. rw_16 .le. 1.0d-8),  &
           1.0d0, 0.0d0)

      allocate(tf_round(1:N,n_carrier_max))
      tf_round = 2.0d0*abs(res_dt)/carrier_row_roundoff
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_round,     &
           terms_full=tf_round)
      call check_absolute('a_round_off_row_is_accepted',                  &
           logical_as_double(vd_round%accepted), 1.0d0, 0.0d0)
      call check_absolute('a_round_off_row_is_counted_as_round_off'//     &
           '_limited',                                                    &
           logical_as_double(vd_round%n_rows_roundoff_limited .ge. 1      &
                             .and. vd_round%n_rows_above .eq. 0), 1.0d0,  &
           0.0d0)
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_mixed)
      call check_absolute('without_the_full_terms_the_same_rows_are'//    &
           '_refused',                                                    &
           logical_as_double(vd_mixed%accepted .or.                       &
                             vd_mixed%n_rows_roundoff_limited .ne. 0),    &
           0.0d0, 0.0d0)
      tf_round(1,ic_H2) = abs(res_dt(1,ic_H2))
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_mixed,     &
           terms_full=tf_round)
      call check_absolute('one_row_above_its_round_off_is_still'//        &
           '_refused', logical_as_double(vd_mixed%accepted), 0.0d0, 0.0d0)
      call check_absolute('and_it_is_the_ordinary_refusal',               &
           dble(vd_mixed%reason), dble(carrier_reject_row_residual),      &
           0.0d0)
      call check_absolute('the_refusal_names_the_row_above_its_round'//   &
           '_off',                                                        &
           logical_as_double(vd_mixed%jworst .eq. 1 .and.                 &
                             vd_mixed%icworst .eq. ic_H2 .and.            &
                             vd_mixed%n_rows_above .eq. 1), 1.0d0, 0.0d0)
      call check_relative('the_named_row_keeps_its_physical_measure',     &
           vd_mixed%rworst, 1.0d-7, 1.0d-12)

      call carrier_record_refused_rows(res_dt, tphys_dt, tf_round,        &
                                       roundoff_of(res_dt, tf_round),     &
                                       unconstrained)
      nrefused = carrier_refused_rows()
      call check_absolute('a_refusal_records_exactly_the_rows_it'//       &
           '_rests_on', dble(nrefused), 1.0d0, 0.0d0)
      call carrier_refused_row(1, jref, icref, res_ref, phys_ref,         &
                               full_ref)
      write(*,'(a,i0,a,a,a,es12.4,a,es12.4,a,es12.4)')                    &
           ' (carrier_retry) the refused row: cell ', jref, ', carrier ', &
           trim(carrier_name(icref)), ', |res| ', res_ref,                &
           ', physical terms ', phys_ref, ', full terms ', full_ref
      call check_absolute('the_recorded_row_is_the_one_above_its_own'//   &
           '_round_off',                                                  &
           logical_as_double(jref .eq. 1 .and. icref .eq. ic_H2), 1.0d0,  &
           0.0d0)
      call check_relative('the_recorded_row_carries_its_physical_scale',  &
           res_ref/phys_ref, 1.0d-7, 1.0d-12)
      call check_absolute('the_recorded_row_carries_its_full_scale',      &
           logical_as_double(res_ref .gt.                                 &
                             carrier_row_roundoff*full_ref), 1.0d0, 0.0d0)

      count = 77
      carrier_transport_stop_on_failure = .false.

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_co_domain_record(nout_co0, nhot_co0, nhep_co0,         &
                                    ratio_co0, r_co0, form_co0)
      call carrier_exhausted_record(nex0, ex_first, ex_last, ex_j, ex_ic, &
                                    ex_ratio, ex_phys, ex_full)
      nstop0 = carrier_transport_stops_suppressed
      run_mode = run_mode_phys
      carrier_reject_leading_attempts_for_test = 10000
      call photochemical_transport_step(rho, v, f_sp, dt_code, st_phys)
      carrier_reject_leading_attempts_for_test = 0
      nstop_phys = carrier_transport_stops_suppressed
      call carrier_exhausted_record(nex_phys, ex_first, ex_last, ex_j,    &
                                    ex_ic, ex_ratio, ex_phys, ex_full)
      call carrier_co_domain_record(nout_co, nhot_co, nhep_co,            &
                                    ratio_co, r_co, form_co)
      call check_absolute('an_exhausted_interval_in_phys_mode_is'//       &
           '_reported_exhausted', dble(st_phys),                          &
           dble(carrier_interval_exhausted), 0.0d0)
      call check_absolute('and_the_status_is_kept_for_a_caller_that'//    &
           '_did_not_take_it', dble(carrier_last_interval_status()),      &
           dble(carrier_interval_exhausted), 0.0d0)
      call check_absolute('and_the_run_reaches_the_stop',                 &
           dble(nstop_phys - nstop0), 1.0d0, 0.0d0)
      call check_absolute('the_composition_is_not_written_in_phys_mode',  &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('the_CO_domain_record_survives_the_restore'//   &
           '_in_phys_mode', dble(nout_co - nout_co0), 0.0d0, 0.0d0)
      call check_absolute('and_its_worst_ratio_in_phys_mode',             &
           ratio_co - ratio_co0, 0.0d0, 0.0d0)
      call check_absolute('the_carrier_history_is_not_marked_in'//        &
           '_phys_mode',                                                  &
           logical_as_double(carrier_history_certifiable()), 1.0d0, 0.0d0)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_co_domain_record(nout_co0, nhot_co0, nhep_co0,         &
                                    ratio_co0, r_co0, form_co0)
      run_mode = run_mode_init
      carrier_reject_leading_attempts_for_test = 10000
      call photochemical_transport_step(rho, v, f_sp, dt_code, st_init)
      carrier_reject_leading_attempts_for_test = 0
      nstop_init = carrier_transport_stops_suppressed
      call carrier_exhausted_record(nex_init, ex_first, ex_last, ex_j,    &
                                    ex_ic, ex_ratio, ex_phys, ex_full)
      call carrier_co_domain_record(nout_co, nhot_co, nhep_co,            &
                                    ratio_co, r_co, form_co)
      call check_absolute('an_exhausted_interval_in_init_mode_is'//       &
           '_reported_exhausted', dble(st_init),                          &
           dble(carrier_interval_exhausted), 0.0d0)
      call check_absolute('and_the_run_does_not_stop',                    &
           dble(nstop_init - nstop_phys), 0.0d0, 0.0d0)
      call check_absolute('the_interval_is_not_written_in_init_mode',     &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('the_attempts_ledger_counts_the_exhausted'//    &
           '_interval', dble(nex_init - nex_phys), 1.0d0, 0.0d0)
      call check_absolute('the_ledger_names_the_step',                    &
           dble(ex_last), dble(count), 0.0d0)
      call check_absolute('the_ledger_keeps_the_first_step_it_saw',       &
           dble(ex_first), dble(count), 0.0d0)
      call check_absolute('the_ledger_names_a_row',                       &
           logical_as_double(ex_j .ge. 1 .and. ex_ic .ge. 1), 1.0d0,      &
           0.0d0)
      call check_absolute('a_forced_refusal_rests_on_no_row',             &
           dble(carrier_refused_rows()), 0.0d0, 0.0d0)
      call check_absolute('so_the_ledger_names_the_verdicts_worst_row'//  &
           '_with_no_scales',                                             &
           logical_as_double(ex_phys .eq. 0.0d0 .and.                     &
                             ex_full .eq. 0.0d0), 1.0d0, 0.0d0)
      call check_absolute('and_the_ratio_it_kept_is_that_rows_measure',   &
           logical_as_double(ex_ratio .ge. 0.0d0), 1.0d0, 0.0d0)
      call check_absolute('the_carrier_history_is_marked_not'//           &
           '_certifiable',                                                &
           logical_as_double(carrier_history_certifiable()), 0.0d0, 0.0d0)
      call check_absolute('the_CO_domain_record_survives_the_restore'//   &
           '_in_init_mode', dble(nout_co - nout_co0), 0.0d0, 0.0d0)
      call check_absolute('and_its_worst_ratio_in_init_mode',             &
           ratio_co - ratio_co0, 0.0d0, 0.0d0)

      call photochemical_transport_step(rho, v, f_sp, dt_code, st_init)
      call check_absolute('the_next_interval_is_covered',                 &
           dble(st_init), dble(carrier_interval_covered), 0.0d0)
      call check_absolute('and_it_starts_from_the_entry_carriers',        &
           carrier_row_difference(f_sp, f_one), 0.0d0, 0.0d0)

      carrier_transport_stop_on_failure = .true.

      vd_run = carrier_last_verdict_of_run()
      reason_before = vd_run%reason
      call check_absolute('the_run_level_record_holds_a_real_verdict'//   &
           '_from_the_interval_before',                                   &
           logical_as_double(reason_before .ne. carrier_no_interval),     &
           1.0d0, 0.0d0)

      thereis_mol = .false.
      st_nothing  = -999
      call photochemical_transport_step(rho, v, f_sp, dt_code, st_nothing)
      vd_step = carrier_step_verdict()
      vd_run  = carrier_last_verdict_of_run()
      call check_absolute('a_step_with_no_carriers_reports_a_covered'//   &
           '_interval', dble(st_nothing), dble(carrier_interval_covered), &
           0.0d0)
      call check_absolute('and_the_kept_status_says_the_same',            &
           dble(carrier_last_interval_status()),                          &
           dble(carrier_interval_covered), 0.0d0)
      call check_absolute('and_the_step_verdict_is_no_interval',          &
           dble(vd_step%reason), dble(carrier_no_interval), 0.0d0)
      call check_absolute('and_it_is_not_an_acceptance',                  &
           logical_as_double(vd_step%accepted), 0.0d0, 0.0d0)
      call check_absolute('while_the_run_level_record_still_holds_the'//  &
           '_last_interval_of_the_run',                                   &
           dble(vd_run%reason), dble(reason_before), 0.0d0)
      thereis_mol = .true.

      carrier_transport = .false.
      st_nothing        = -999
      call photochemical_transport_step(rho, v, f_sp, dt_code, st_nothing)
      vd_step = carrier_step_verdict()
      call check_absolute('a_step_with_the_transport_option_off'//        &
           '_reports_a_covered_interval', dble(st_nothing),               &
           dble(carrier_interval_covered), 0.0d0)
      call check_absolute('and_its_step_verdict_is_no_interval_too',      &
           dble(vd_step%reason), dble(carrier_no_interval), 0.0d0)
      carrier_transport = .true.

      call carrier_co_domain_record(nout_co0, nhot_co0, nhep_co0,         &
                                    ratio_co0, r_co0, form_co0, 1)
      call carrier_co_domain_record(nout_co, nhot_co, nhep_co,            &
                                    ratio_co, r_co, form_co, 2)
      call carrier_co_domain_perturb_for_test(7, 19)
      call carrier_co_domain_record(dout_f1, dhot_f1, dhep_f1,            &
                                    drat_f1, dr_f1, dform_f1, 1)
      call carrier_co_domain_record(dout_f2, dhot_f2, dhep_f2,            &
                                    drat_f2, dr_f2, dform_f2, 2)
      call carrier_co_domain_record(dout_tot, dhot_tot, dhep_tot,         &
                                    drat_tot, dr_tot, dform_tot)
      call check_absolute('the_domain_record_keeps_the_initialization'//  &
           '_family_apart', dble(dout_f1 - nout_co0), 7.0d0, 0.0d0)
      call check_absolute('and_the_physical_family_apart',                &
           dble(dout_f2 - nout_co), 19.0d0, 0.0d0)
      call check_absolute('the_run_total_counts_the_two_families'//       &
           '_together', dble(dout_tot), dble(dout_f1 + dout_f2), 0.0d0)
      call check_absolute('and_so_do_the_hot_and_He_II_led_counts',       &
           dble(dhot_tot + dhep_tot),                                     &
           dble(dhot_f1 + dhot_f2 + dhep_f1 + dhep_f2), 0.0d0)
      call check_absolute('the_worst_ratio_of_the_run_is_the_larger'//    &
           '_of_the_two_and_not_their_sum', drat_tot,                     &
           max(drat_f1, drat_f2), 0.0d0)
      call check_absolute('and_the_radius_reported_with_it_is_that'//     &
           '_familys', dr_tot, merge(dr_f1, dr_f2,                        &
           drat_f1 .ge. drat_f2), 0.0d0)
      call carrier_checkpoint_take(chk_after)
      call carrier_checkpoint_restore(chk_after)
      call carrier_co_domain_record(dout_f1, dhot_f1, dhep_f1,            &
                                    drat_f1, dr_f1, dform_f1)
      call check_absolute('the_domain_record_is_not_a_checkpointed'//     &
           '_quantity', dble(dout_f1), dble(dout_tot), 0.0d0)


      allocate(nnuc0(1-Ng:N+Ng,n_element), nnuc1(1-Ng:N+Ng,n_element))
      allocate(nchg0(1-Ng:N+Ng), nchg1(1-Ng:N+Ng))
      allocate(rcomp0(1-Ng:N+Ng), rcomp1(1-Ng:N+Ng))
      allocate(v_relax(1-Ng:N+Ng))
      v_relax = 20.0d0*r
      run_mode = run_mode_phys
      carrier_transport_stop_on_failure = .false.
      call element_nuclei_and_charge(rho, f_sp0, nnuc0, nchg0, rcomp0)
      write(*,'(A,ES24.16)') 'AUDIT synthetic_entry_mass_closure=', &
           maxval(abs(rcomp0(1:N)/(rho(1:N)*n0)-1d0))
      call seed_mass_row_of_the_column(v_relax)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d-2, dr_wide,      &
                                           ns_wide, oc_wide)
      write(*,'(a,es10.2,a,es12.4,a,i0,a,a)') ' (carrier_retry) trust ',  &
           1.0d-2, ': drift ', dr_wide, ' in ', ns_wide, ' kept steps,'// &
           ' ending on '//trim(carrier_relax_outcome_text(oc_wide))
      call check_absolute('a_pass_returns_a_state_inside_its_bound',      &
           logical_as_double(dr_wide .le. 1.0d-2*(1.0d0 + 1.0d-12)),      &
           1.0d0, 0.0d0)
      call check_positive('and_a_shortened_trial_advanced_the_carriers',  &
           dr_wide)
      call check_absolute('and_the_bound_is_what_ended_it',               &
           dble(oc_wide), dble(carrier_relax_movement_bound), 0.0d0)
      call check_absolute('the_bounded_state_is_finite_and_nonnegative',  &
           logical_as_double(all(f_sp(1:N,:) .ge. 0.0d0) .and.            &
                             .not. any(f_sp(1:N,:) .ne. f_sp(1:N,:))),    &
           1.0d0, 0.0d0)
      call element_nuclei_and_charge(rho, f_sp, nnuc1, nchg1, rcomp1)
      call check_absolute('the_bounded_pass_conserves_every_element',     &
           element_deviation(nnuc0, nnuc1), 0.0d0, 1.0d-12)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d-3, dr_tight,     &
                                           ns_tight, oc_tight)
      write(*,'(a,es10.2,a,es12.4,a,i0,a,a)') ' (carrier_retry) trust ',  &
           1.0d-3, ': drift ', dr_tight, ' in ', ns_tight,               &
           ' kept steps, ending on '//                                    &
           trim(carrier_relax_outcome_text(oc_tight))
      call check_absolute('a_tenfold_tighter_bound_also_holds',           &
           logical_as_double(dr_tight .le. 1.0d-3*(1.0d0 + 1.0d-12)),      &
           1.0d0, 0.0d0)
      call check_positive('and_it_too_advanced_the_carriers', dr_tight)
      call check_absolute('and_the_returned_movement_falls_with_the'//    &
           '_bound', logical_as_double(dr_tight .lt. dr_wide), 1.0d0,     &
           0.0d0)
      call element_nuclei_and_charge(rho, f_sp, nnuc1, nchg1, rcomp1)
      call check_absolute('the_tightly_bounded_pass_conserves_every'//    &
           '_element', element_deviation(nnuc0, nnuc1), 0.0d0, 1.0d-12)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d-4, dr_none,      &
                                           ns_none, oc_none)
      call check_absolute('a_bound_no_trial_can_meet_keeps_no_step',      &
           dble(ns_none), 0.0d0, 0.0d0)
      call check_absolute('and_hands_back_the_entry_composition',         &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('and_reports_the_bound_as_what_ended_it',       &
           dble(oc_none), dble(carrier_relax_movement_bound), 0.0d0)

      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 0.0d0, dr_zero,       &
                                           ns_zero, oc_zero)
      call check_absolute('a_pass_that_may_not_move_keeps_no_step',       &
           dble(ns_zero), 0.0d0, 0.0d0)
      call check_absolute('and_returns_the_entry_composition_bit_for'//   &
           '_bit', maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('and_reports_no_drift', dr_zero, 0.0d0, 0.0d0)
      call check_absolute('and_names_the_bound_as_what_ended_it',         &
           dble(oc_zero), dble(carrier_relax_movement_bound), 0.0d0)

      allocate(v_fast(1-Ng:N+Ng))
      v_fast = 2000.0d0*r
      call seed_mass_row_of_the_column(v_fast)
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      run_mode = run_mode_init
      carrier_reject_leading_attempts_for_test = 1000000
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_fast, p_col, T_col,     &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d0, dr_ref,        &
                                           ns_ref, oc_ref)
      carrier_reject_leading_attempts_for_test = 0
      run_mode = run_mode_phys
      write(*,'(a,i0,a,a)') ' (carrier_retry) a pass whose every'//       &
           ' interval is refused keeps ', ns_ref, ' steps, ending on '//  &
           trim(carrier_relax_outcome_text(oc_ref))
      call check_absolute('a_refusal_at_every_length_ends_the_pass_by'//  &
           '_name', dble(oc_ref), dble(carrier_relax_interval_refused),   &
           0.0d0)
      call check_absolute('a_refused_pass_keeps_no_step', dble(ns_ref),   &
           0.0d0, 0.0d0)
      call check_absolute('a_refused_pass_returns_the_entry'//            &
           '_composition', maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('a_refused_pass_reports_no_drift', dr_ref,      &
           0.0d0, 0.0d0)


      carrier_transport = .false.
      f_sp = f_sp0
      call restore_frozen_background()
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d0, dr_zero,       &
                                           ns_zero, oc_zero)
      carrier_transport = .true.
      carrier_transport_stop_on_failure = .true.
      call check_absolute('a_pass_with_no_transported_carrier_says_so',   &
           dble(oc_zero), dble(carrier_relax_nothing_to_advance), 0.0d0)
      call check_absolute('and_moves_nothing', maxval(abs(f_sp - f_sp0)), &
           0.0d0, 0.0d0)


      allocate(f_swept(1-Ng:N+Ng,n_species))
      allocate(nhi_r(1-Ng:N+Ng), nhii_r(1-Ng:N+Ng))
      allocate(nhei_r(1-Ng:N+Ng), nheii_r(1-Ng:N+Ng))
      allocate(nheiii_r(1-Ng:N+Ng), nheiTR_r(1-Ng:N+Ng))
      allocate(ne_r(1-Ng:N+Ng), ntot_r(1-Ng:N+Ng), T_ret(1-Ng:N+Ng))
      allocate(nm_r(1-Ng:N+Ng,n_mion))

      call seed_mass_row_of_the_column(v_relax)
      call seed_background_at_the_entry_composition()
      allocate(audit_energy(1:N))
      do jj=1,N
         audit_energy(jj)=energy_density_from_pressure(jj,rho(jj),p_col(jj))
      enddo
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call relax_photochemical_composition(rho, v_relax, p_col, T_col,    &
                                           f_sp, heat_col, cool_col,      &
                                           eta_col, 1.0d-2, dr_wide,      &
                                           ns_wide, oc_wide)
      call check_positive('the_consistency_rows_read_a_pass_that_kept'//  &
           '_a_step', dble(ns_wide))

      nhei_r   = 0.0d0
      nheii_r  = 0.0d0
      nheiii_r = 0.0d0
      nheiTR_r = 0.0d0
      call get_species_densities(rho, f_sp, nhi_r, nhii_r, nhei_r,        &
                                 nheii_r, nheiii_r, nheiTR_r, nm_r,       &
                                 ne_r, ntot_r)
      call comp_T_from_p(p_col, ntot_r, ne_r, T_ret)
      write(*,'(A,ES24.16)') 'AUDIT returned_temperature_EOS_error=', &
           maxval(abs(T_col(1:N)-T_ret(1:N))/abs(T_ret(1:N)))
      audit_p_error=0d0; audit_energy_error=0d0
      do jj=1,N
         audit_p_error=max(audit_p_error, &
           abs(pressure_from_energy_density(jj,rho(jj),audit_energy(jj))-p_col(jj))/abs(p_col(jj)))
         audit_energy_error=max(audit_energy_error, &
           abs(energy_density_from_pressure(jj,rho(jj),p_col(jj))-audit_energy(jj))/abs(audit_energy(jj)))
      enddo
      write(*,'(A,ES24.16)') 'AUDIT fixed_energy_pressure_error=',audit_p_error
      write(*,'(A,ES24.16)') 'AUDIT fixed_pressure_energy_change=',audit_energy_error
      call mass_row_cell_verdict(0.5d0,1d0,audit_tol,audit_dist,audit_within,audit_anchored)
      write(*,'(A,2ES24.16,2L2)') 'AUDIT ceiling_q=0.5_floor=1 tol_distance_within_anchored=', &
           audit_tol,audit_dist,audit_within,audit_anchored


      dev_ntot  = 0.0d0
      move_pass = 0.0d0
      do jj = 1, N
         dev_ntot  = max(dev_ntot, abs(bg_cell(jj)%ntot/n0 - ntot_r(jj))  &
                                   /ntot_r(jj))
         move_pass = max(move_pass, abs(bg_cell(jj)%ntot/n0 - ntot0_c(jj))&
                                    /ntot0_c(jj))
      enddo
      write(*,'(a,2es12.4)') ' (carrier_retry) background third-body'//   &
           ' density against the returned and the entry composition: ',   &
           dev_ntot, move_pass
      call check_absolute('the_background_third_body_density_follows'//   &
           '_the_composition',                                            &
           logical_as_double(dev_ntot .lt. move_pass), 1.0d0, 0.0d0)

      dev_TK     = 0.0d0
      move_sweep = 0.0d0
      do jj = 1, N
         dev_TK     = max(dev_TK, abs(bg_cell(jj)%T_K/T0 - T_ret(jj))     &
                                  /T_ret(jj))
         move_sweep = max(move_sweep, abs(bg_cell(jj)%T_K/T0 - T_ent(jj)) &
                                      /T_ent(jj))
      enddo
      write(*,'(a,2es12.4)') ' (carrier_retry) background temperature'//  &
           ' against the returned and the entry composition: ',           &
           dev_TK, move_sweep
      call check_absolute('the_background_temperature_follows_the'//      &
           '_composition',                                                &
           logical_as_double(dev_TK .lt. move_sweep), 1.0d0, 0.0d0)

      move_pass = maxval(abs(f_sp(1:N,:) - f_sp0(1:N,:)))
      f_swept   = f_sp
      call ioniz_eq(T_ret, rho, f_swept, heat_col, cool_col, eta_col)
      move_sweep = maxval(abs(f_swept(1:N,:) - f_sp(1:N,:)))
      write(*,'(a,es12.4,a,es12.4)') ' (carrier_retry) the pass moved'//   &
           ' the composition by ', move_pass, ' and one further sweep'//  &
           ' moves it by ', move_sweep
      call check_absolute('one_sweep_after_the_pass_moves_less_than'//    &
           '_the_pass_did',                                               &
           logical_as_double(move_sweep .lt. move_pass), 1.0d0, 0.0d0)
      call check_absolute('and_moves_no_more_than_the_cell_solve_s_own'// &
           '_step_criterion',                                             &
           logical_as_double(move_sweep .lt. sqrt(epsilon(1.0d0))),       &
           1.0d0, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a)') 'carrier_retry: FAILED'
         stop 1
      endif
      write(*,'(a)') 'carrier_retry: every assertion passed'

      contains


      function roundoff_of(res, tfull) result(flag)
      real*8, intent(in) :: res(:,:), tfull(:,:)
      logical :: flag(size(res,1),size(res,2))
      flag = (abs(res) .le. carrier_row_roundoff*tfull)
      end function roundoff_of


      double precision function carrier_row_difference(a, b) result(d)
      real*8, intent(in) :: a(1-Ng:,:), b(1-Ng:,:)
      integer :: k, isp(6)
      real*8  :: sc
      isp = (/ isp_HI, isp_HII, isp_H2, isp_OH, isp_H2O, isp_CO /)
      d = 0.0d0
      do k = 1, 6
         sc = maxval(abs(f_sp0(1:N,isp(k))))
         if (sc .le. 0.0d0) cycle
         d = max(d, maxval(abs(a(1:N,isp(k)) - b(1:N,isp(k))))/sc)
      enddo
      end function carrier_row_difference


      double precision function element_deviation(a, b) result(d)
      real*8, intent(in) :: a(1-Ng:,:), b(1-Ng:,:)
      integer :: j, ie
      d = 0.0d0
      do ie = 1, n_element
         do j = 1, N
            if (a(j,ie) .le. 0.0d0) cycle
            d = max(d, abs(b(j,ie) - a(j,ie))/a(j,ie))
         enddo
      enddo
      end function element_deviation


      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double


      subroutine setup_globals()
      integer :: j
      N   = ncell
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      b0  = 10.0d0
      p_base_bar         = 1.0d-6
      spherical_domain   = .true.
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_metals     = .true.
      thereis_mol        = .true.
      thereis_oxychem    = .false.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport  = .true.
      ionization_transport = .true.
      weno_mode          = 0
      j_min              = 1
      call allocate_grid_arrays
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.02d0*dble(j)
         r_edg(j) = 1.0d0 + 0.02d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.02d0
      enddo
      kzz_cell = 0.0d0
      call set_gravity_grid
      grid_type  = 'Uniform'
      rec_method = 'WENO3'
      use_plm    = .false.
      use_weno3  = .true.
      flux       = 'ROE'
      CFL        = 0.6d0
      r_max      = r_edg(N)
      r_esc      = r(N)
      r_flux     = r(1)
      Mp         = 1.0d30
      q0         = 1.0d0
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      is_PL_sed    = .true.
      thereis_Xray = .true.
      e_low  = 13.60d0
      e_mid  = 123.98d0
      e_top  = 1.24d3
      PLind  = -1.0d0
      LX     = 27.20d0
      LEUV   = 27.93d0
      a_orb  = 0.0480d0*AU
      call set_energy_vectors
      N_eq = 7 + 2*n_melem
      lwa  = (N_eq*(3*N_eq + 13))/2
      allocate(sys_sol(N_eq), sys_x(N_eq), wa(lwa))
      call cx_init
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      bg_ready = .true.
      end subroutine setup_globals


      subroutine build_molecular_hydrogen_column()
      integer :: j
      allocate(rho(1-Ng:N+Ng), v(1-Ng:N+Ng), dt_code(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species), f_sp0(1-Ng:N+Ng,n_species))
      allocate(f_one(1-Ng:N+Ng,n_species), f_half(1-Ng:N+Ng,n_species))
      allocate(f_quarter(1-Ng:N+Ng,n_species))
      allocate(f_retried(1-Ng:N+Ng,n_species))
      do j = 1-Ng, N+Ng
         rho(j) = exp(-4.0d0*(r(j) - 1.0d0))
         v(j)   = 0.02d0*r(j)
      enddo
      dt_code = 1.0d-3
      f_sp = 0.0d0
      f_sp(:,isp_HI)   = 0.20d0
      f_sp(:,isp_HII)  = 1.0d-6
      f_sp(:,isp_H2)   = 0.40d0
      f_sp(:,isp_H2p)  = 1.0d-10
      f_sp(:,isp_H3p)  = 1.0d-10
      f_sp(:,isp_HeI)  = 0.0793d0
      f_sp(:,isp_HeII) = 1.0d-8
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      f_sp0 = f_sp
      allocate(p_col(1-Ng:N+Ng), T_col(1-Ng:N+Ng))
      allocate(heat_col(1-Ng:N+Ng), cool_col(1-Ng:N+Ng))
      allocate(eta_col(1-Ng:N+Ng))
      call seed_pressure_of_the_column()
      end subroutine build_molecular_hydrogen_column


      subroutine seed_pressure_of_the_column()
      real*8, allocatable :: nhi_s(:), nhii_s(:), nhei_s(:), nheii_s(:)
      real*8, allocatable :: nheiii_s(:), nheiTR_s(:), ne_s(:), ntot_s(:)
      real*8, allocatable :: nm_s(:,:)
      integer :: j
      allocate(nhi_s(1-Ng:N+Ng), nhii_s(1-Ng:N+Ng))
      allocate(nhei_s(1-Ng:N+Ng), nheii_s(1-Ng:N+Ng))
      allocate(nheiii_s(1-Ng:N+Ng), nheiTR_s(1-Ng:N+Ng))
      allocate(ne_s(1-Ng:N+Ng), ntot_s(1-Ng:N+Ng))
      allocate(nm_s(1-Ng:N+Ng,n_mion))
      nhei_s   = 0.0d0
      nheii_s  = 0.0d0
      nheiii_s = 0.0d0
      nheiTR_s = 0.0d0
      call get_species_densities(rho, f_sp0, nhi_s, nhii_s, nhei_s,       &
                                 nheii_s, nheiii_s, nheiTR_s, nm_s,       &
                                 ne_s, ntot_s)
      do j = 1-Ng, N+Ng
         T_col(j) = bg_cell(j)%T_K/T0
      enddo
      call comp_p_from_T(T_col, ntot_s, ne_s, p_col)
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine seed_pressure_of_the_column


      subroutine seed_background_at_the_entry_composition()
      real*8, allocatable :: nhi_e(:), nhii_e(:), nhei_e(:), nheii_e(:)
      real*8, allocatable :: nheiii_e(:), nheiTR_e(:), ne_e(:)
      real*8, allocatable :: nm_e(:,:)
      integer :: j
      if (.not. allocated(ntot0_c)) then
         allocate(ntot0_c(1-Ng:N+Ng), T_ent(1-Ng:N+Ng))
      endif
      allocate(nhi_e(1-Ng:N+Ng), nhii_e(1-Ng:N+Ng))
      allocate(nhei_e(1-Ng:N+Ng), nheii_e(1-Ng:N+Ng))
      allocate(nheiii_e(1-Ng:N+Ng), nheiTR_e(1-Ng:N+Ng))
      allocate(ne_e(1-Ng:N+Ng), nm_e(1-Ng:N+Ng,n_mion))
      nhei_e   = 0.0d0
      nheii_e  = 0.0d0
      nheiii_e = 0.0d0
      nheiTR_e = 0.0d0
      call get_species_densities(rho, f_sp0, nhi_e, nhii_e, nhei_e,       &
                                 nheii_e, nheiii_e, nheiTR_e, nm_e,       &
                                 ne_e, ntot0_c)
      call comp_T_from_p(p_col, ntot0_c, ne_e, T_ent)
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
         bg_cell(j)%ntot = ntot0_c(j)*n0
         bg_cell(j)%T_K  = T_ent(j)*T0
         T_col(j)        = T_ent(j)
      enddo
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine seed_background_at_the_entry_composition


      subroutine restore_frozen_background()
      integer :: j
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
         T_col(j) = bg_cell(j)%T_K/T0
      enddo
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine restore_frozen_background


      subroutine seed_mass_row_of_the_column(v_wind)
      real*8, intent(in)  :: v_wind(1-Ng:N+Ng)
      real*8, allocatable :: u(:,:), Res(:,:)
      real*8, allocatable :: n_part(:), heat(:), cool(:), p_s(:), T_s(:)
      integer :: j
      allocate(u(3,1-Ng:N+Ng), Res(3,1-Ng:N+Ng))
      allocate(n_part(1-Ng:N+Ng), heat(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      allocate(p_s(1-Ng:N+Ng), T_s(1-Ng:N+Ng))
      heat = 0.0d0
      cool = 0.0d0
      do j = 1-Ng, N+Ng
         T_s(j)    = bg_cell(j)%T_K/T0
         n_part(j) = rho(j)
         p_s(j)    = n_part(j)*T_s(j)
         u(1,j)    = rho(j)
         u(2,j)    = rho(j)*v_wind(j)
         u(3,j)    = 0.5d0*rho(j)*v_wind(j)**2                            &
                     + p_s(j)/(gamma_ad - 1.0d0)
      enddo
      n_part_cell1 = n_part(1)
      call set_base_reservoir(p_s(1), T_s(1), 1.0d0, 1.0d0)
      call Apply_BC(u)
      call assemble_residual(u, n_part, heat, cool, Res)
      end subroutine seed_mass_row_of_the_column


      subroutine set_frozen_cell_rates(j)
      integer, intent(in) :: j
      bg_cell(j)%P_HI        = 1.0d-6      ! [1/s]
      bg_cell(j)%P_HeI       = 1.0d-7
      bg_cell(j)%P_HeII      = 0.0d0
      bg_cell(j)%P_HeITR     = 0.0d0
      bg_cell(j)%P_H2        = 1.0d-8
      bg_cell(j)%P_H2_di     = 0.0d0
      bg_cell(j)%P_H2_dd     = 0.0d0
      bg_cell(j)%P_H2_nd     = 0.0d0
      bg_cell(j)%k_LW        = 0.0d0
      bg_cell(j)%rchiiB      = 2.6d-13     ! [cm^3/s] at ~1e4 K
      bg_cell(j)%rcheiiB     = 4.3d-13
      bg_cell(j)%rcheiiiB    = 2.2d-12
      bg_cell(j)%rcheiTR     = 0.0d0
      bg_cell(j)%a_ion_HI    = 0.0d0
      bg_cell(j)%a_ion_HeI   = 0.0d0
      bg_cell(j)%a_ion_HeII  = 0.0d0
      bg_cell(j)%a_ion_HeITR = 0.0d0
      bg_cell(j)%q13         = 0.0d0
      bg_cell(j)%q31a        = 0.0d0
      bg_cell(j)%q31b        = 0.0d0
      bg_cell(j)%Q31         = 0.0d0
      bg_cell(j)%A31         = 0.0d0
      bg_cell(j)%kcx_He0_Hp  = 0.0d0
      bg_cell(j)%kcx_Hep_H0  = 0.0d0
      bg_cell(j)%nh          = 0.0d0
      bg_cell(j)%nhe         = 0.0d0
      bg_cell(j)%n_ofam      = 0.0d0
      bg_cell(j)%n_co        = 0.0d0
      bg_cell(j)%x_h2_fixed  = .false.
      bg_cell(j)%x_ox_fixed  = .false.
      bg_cell(j)%x_hp_fixed  = .false.
      bg_cell(j)%x_h2_fix    = 0.0d0
      bg_cell(j)%x_oh_fix    = 0.0d0
      bg_cell(j)%x_h2o_fix   = 0.0d0
      bg_cell(j)%x_hp_fix    = 0.0d0
      bg_cell(j)%T_K         = 1500.0d0
      bg_cell(j)%ntot        = 0.5d0*rho(j)*n0
      end subroutine set_frozen_cell_rates

      end program carrier_state_contract_probe

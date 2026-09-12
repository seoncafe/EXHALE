      program carrier_retry
      ! THE CARRIER-LOCAL CHECKPOINT AND THE REJECTION-AND-RETRY CONTROLLER
      ! (PLAN_20260906_rev2 step A3), measured on a real molecular column.
      !
      ! WHAT IS BEING TESTED.
      !  (1) The checkpoint is complete.  Every item it covers is given a
      !      distinct change by carrier_perturb_checkpointed_state_for_test,
      !      including the allocation of arrays that were not allocated, and
      !      the restore has to put all of them back bit for bit.  An item
      !      missing from either half of the pair fails here.
      !  (2) A refused substep leaves nothing behind.  The controller is made
      !      to refuse the first attempt of an interval through the test hook,
      !      and the interval must still be covered: the entry state comes
      !      back, the substep halves, and two accepted halves cover it, with
      !      the elapsed fraction exactly 1.
      !  (3) The retried result is the subdivided integration.  The state the
      !      controller reaches after a refused full step and two accepted
      !      halves is compared with two accepted halves taken from the start
      !      over THE SAME frozen background.  They agree to 1e-10, which is
      !      restoration measured on the physics rather than on the memory.
      !  (4) Exhaustion completes nothing.  With every attempt refused, the
      !      controller halves down to dt/2**carrier_retry_max, restores the
      !      state the interval was entered with, reports the interval as not
      !      covered, and writes no composition.  The attempt statistics
      !      carry the attempts, and the domain record of the CO
      !      destruction model, which is not part of the checkpoint, is not
      !      rolled back with the state: whether a state was in domain has
      !      an answer whether or not the step that read it was accepted.
      !  (5) Time accuracy (review R10).  The carrier substeps integrate the
      !      transport-chemistry system at a background frozen at the entry
      !      state, so subdividing refines the time integration of that
      !      operator and does not change the splitting with the operators
      !      around it.  The step-doubling difference between one step of dt
      !      and two of dt/2 is measured and printed, and a second halving
      !      has to reduce it, which is first order or better.  No tolerance
      !      is asserted on its size.
      !
      ! THE COLUMN.  Twelve cells of molecular hydrogen and helium with the
      ! transported proton and no metals, the frozen cell state set here
      ! field by field, an outflowing velocity and a step long enough that
      ! the H2 and H+ rows move by more than round-off over the interval.
      ! Every assertion below is a statement about the controller, so the
      ! column only has to be one the solve can integrate.

      use global_parameters
      use species_table
      use ionization_equilibrium, only: bg_cell, bg_ready,                &
                                        ioniz_eq_allocate_arrays, ioniz_eq,&
                                        ieq_res_tol
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
           carrier_history_reset_for_test,                                &
           carrier_transport_stop_on_failure,                             &
           carrier_transport_stops_suppressed,                            &
           carrier_step_verdict, carrier_last_verdict_of_run,             &
           carrier_no_interval,                                           &
           relax_photochemical_composition,                               &
           carrier_relax_movement_bound, carrier_relax_fixed_point,       &
           carrier_relax_step_budget, carrier_relax_interval_refused,     &
           carrier_relax_nothing_to_advance, carrier_relax_outcome_text,  &
           carrier_relax_chemistry_refused,                               &
           chem_cycles_cap_for_test, n_chem_last_reason,                  &
           chem_closure_exhausted, chem_closure_reason_text,               &
           chem_last_increment, chem_last_offsimplex, chem_last_mol_clamped, &
           chem_last_viol_worst
      use test_columns,  only: column_carrying_its_own_density,           &
                               column_mass_closure
      use Conversion,     only: W_to_U
      use caloric_eos,    only: pressure_from_energy_density
      use element_census, only: element_nuclei_and_charge, n_element
      use gravity_grid_construction, only: set_gravity_grid
      use energy_vectors_construct,  only: set_energy_vectors
      use charge_exchange,           only: cx_init
      use composition,               only: get_species_densities,         &
                                           comp_p_from_T, comp_T_from_p
      use base_boundary,             only: set_base_reservoir
      use BC_Apply,                  only: Apply_BC
      use steady_residual_mod,       only: assemble_residual
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
      ! The domain record of the one-sided CO destruction model, read
      ! before and after an interval that is refused in every attempt.
      integer :: nout_co0, nhot_co0, nhep_co0
      integer :: nout_co,  nhot_co,  nhep_co
      ! The two row scales of one assembly, and the imbalance a state is
      ! judged on when it is a fixed fraction of the row's physical terms.
      real*8,  allocatable :: tfull_dt(:,:), tphys_dt(:,:)
      real*8,  allocatable :: tfull_short(:,:), tphys_short(:,:)
      real*8,  allocatable :: res_dt(:,:), res_short(:,:)
      real*8,  allocatable :: dt_short(:)
      logical, allocatable :: unconstrained(:)
      type(carrier_verdict) :: vd_full_old, vd_full_new
      type(carrier_verdict) :: vd_short_old, vd_short_new
      real*8  :: ratio_co0, r_co0, form_co0
      real*8  :: ratio_co,  r_co,  form_co
      ! (7) the shortest substep the arithmetic can certify: the same
      ! interval integrated at dt, dt/16 and dt/1024, with the physical
      ! measure of the state each one returns and the iterations it cost.
      real*8,  allocatable :: dt_1024(:)
      real*8  :: rw_dt, rw_16, rw_1024
      real*8  :: rs_dt, rs_16, rs_1024, worst_lim
      integer :: nit_dt, nit_16, nit_1024, nlim
      integer :: nsub_dt, nsub_16, nsub_1024
      integer :: nhalve_1024, nrl_1024, nrl_last
      integer :: nrl_total, nrl_substeps
      integer :: nrefused, jref, icref
      real*8  :: res_ref, phys_ref, full_ref
      ! Constructed full-term scales, for the classification of a failing
      ! row as at the round-off of its own terms or not.
      real*8,  allocatable :: tf_round(:,:)
      type(carrier_verdict) :: vd_round, vd_mixed
      logical :: cov_dt, cov_16, cov_1024
      ! (8) what an exhausted interval does in each run mode.
      integer :: st_phys, st_init, nstop0, nstop_phys, nstop_init
      integer :: nex0, nex_phys, nex_init
      integer :: ex_first, ex_last, ex_j, ex_ic
      real*8  :: ex_ratio, ex_phys, ex_full
      ! The domain record read as the run total and family by family.
      integer :: dout_tot, dhot_tot, dhep_tot
      integer :: dout_f1,  dhot_f1,  dhep_f1
      integer :: dout_f2,  dhot_f2,  dhep_f2
      real*8  :: drat_tot, dr_tot, dform_tot
      real*8  :: drat_f1,  dr_f1,  dform_f1
      real*8  :: drat_f2,  dr_f2,  dform_f2
      ! (9) a step the operator has nothing to do in.  st_nothing is set to
      ! a value that is neither status, so that the assertion reads what
      ! the operator wrote and not what this program left there.
      integer :: st_nothing, reason_before
      type(carrier_verdict) :: vd_step, vd_run
      ! (10) the movement bound of a relaxation pass: the drift each pass
      ! reports, what ended it, and the elemental content of the state it
      ! returns.
      real*8  :: dr_wide, dr_tight, dr_none, dr_zero, dr_ref
      integer :: ns_wide, ns_tight, ns_none, ns_zero, ns_ref
      integer :: oc_wide, oc_tight, oc_none, oc_zero, oc_ref
      real*8, allocatable :: nnuc0(:,:), nnuc1(:,:)
      real*8, allocatable :: nchg0(:), nchg1(:), rcomp0(:), rcomp1(:)
      real*8, allocatable :: v_relax(:), v_fast(:)
      ! The hydrodynamic state the relaxation holds fixed while it advances
      ! the carriers, and the sweep outputs it writes beside the
      ! composition: the pressure of the entry state, the temperature of
      ! whatever composition stands beside it, and the heating, cooling and
      ! heating-efficiency columns of the last equilibrium sweep.
      real*8, allocatable :: p_col(:), T_col(:)
      real*8, allocatable :: heat_col(:), cool_col(:), eta_col(:)
      ! (11) the chemistry the pass carries with it: the densities of the
      ! composition a pass hands back, the temperature they imply at the
      ! fixed pressure, and one further sweep taken on that state.
      real*8, allocatable :: f_swept(:,:)
      real*8, allocatable :: nhi_r(:), nhii_r(:), nhei_r(:), nheii_r(:)
      real*8, allocatable :: nheiii_r(:), nheiTR_r(:), ne_r(:), ntot_r(:)
      real*8, allocatable :: nm_r(:,:), T_ret(:)
      real*8, allocatable :: ntot0_c(:), T_ent(:)
      real*8  :: move_pass, move_sweep, dev_ntot, dev_TK
      integer :: jj
      ! The conserved state the relaxation holds (rho, rho v, E) formed
      ! from the column's density, the wind of the section and the entry
      ! pressure by the code's own equation of state.
      real*8, allocatable :: u_col(:,:), W_col(:,:)
      real*8  :: closure_entry, dev_p, dev_T, move_rel, fref
      ! (12) a rejected stationary trial leaves the carrier history alone.
      real*8  :: dr_hist
      integer :: ns_hist, oc_hist, st_trial, run_mode_held, isp
      logical :: hist_before, hist_after

      call setup_globals()
      call build_molecular_hydrogen_column()
      ! THE COLUMN CARRIES ITS OWN DENSITY: every quantitative statement
      ! below rests on a composition that reconstructs rho, so a fraction
      ! added on top of a closed mass budget would make them statements
      ! about a gas that does not exist.
      closure_entry = column_mass_closure(rho, f_sp0)
      write(*,'(a,es12.4)') ' (carrier_retry) mass closure of the'//      &
           ' entry column: ', closure_entry
      call check_absolute('the_column_reconstructs_its_own_density',      &
           closure_entry, 0.0d0, 1.0d-13)

      ! ---- (1) the checkpoint round trip ----------------------------- !
      !
      ! First on a module in which nothing is allocated yet, so that the
      ! deallocation half of the restore is exercised: the perturbation
      ! allocates every array, and putting the checkpoint back has to leave
      ! them unallocated again.
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

      ! Then on a module a real interval has filled, so that the same list
      ! is measured on states the code actually produces.
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

      ! The state every run below starts from: the composition of the
      ! column and the module as it stood before the first interval.
      call carrier_checkpoint_take(chk0)

      ! ---- (2) an injected refusal, and the interval still covered ---- !
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

      ! ---- (3) the same as an explicitly subdivided integration ------- !
      !
      ! Two halves over the SAME frozen background, which is what the
      ! controller took after the refusal.  Calling the transport step
      ! twice with half the step would not be this reference: it would
      ! rebuild the background from the half-advanced composition.
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

      ! ---- (4) exhaustion covers nothing ------------------------------ !
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call carrier_attempt_record(nattempt0, nreject0, naccept0, nretried0)
      ! A nonzero domain record to enter the refused interval with, so that
      ! "not rolled back" is measured and not read off two zeros.
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
      ! The domain record is NOT part of the checkpoint and the restore
      ! must not undo it: whether the state a refused attempt read was in
      ! domain has an answer whether or not that attempt was accepted.
      call check_absolute('the_CO_domain_record_survives_the_restore',    &
           dble(nout_co - nout_co0), 0.0d0, 0.0d0)
      call check_absolute('and_so_does_its_hot_cell_count',               &
           dble(nhot_co - nhot_co0), 0.0d0, 0.0d0)
      call check_absolute('and_its_worst_ratio',                          &
           ratio_co - ratio_co0, 0.0d0, 0.0d0)
      call check_positive('the_domain_record_read_here_is_not_zero',      &
           dble(nout_co0))

      ! ---- (5) time accuracy by step doubling ------------------------- !
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
      ! Printed as numbers, not bounded by a tolerance: this is the size of
      ! the time-integration error of the carrier operator over the
      ! interval, measured by step doubling on this column.
      call check_positive('step_doubling_difference_dt_against_two'//     &
                          '_halves', e1)
      call check_positive('step_doubling_difference_two_halves'//         &
                          '_against_four_quarters', e2)
      ! Halving again reduces it, which is first order or better.  Nothing
      ! stronger is claimed.
      call check_absolute('halving_again_reduces_the_step_doubling'//     &
           '_difference', logical_as_double(e2 .lt. e1), 1.0d0, 0.0d0)
      write(*,'(a,es12.4,a,es12.4,a,f8.3)')                               &
           ' (carrier_retry) step-doubling difference: dt vs 2 x dt/2 = ', &
           e1, ', 2 x dt/2 vs 4 x dt/4 = ', e2, ', ratio = ', e1/max(e2,  &
           1.0d-300)

      ! ---- (6) the scale the returned state is judged on -------------- !
      !
      ! THE MEASURE MUST NOT DEPEND ON THE STEP.  A carrier row is
      ! (n - n_old)/dt = transport + reaction.  The imbalance a returned
      ! state carries -- the flux a neighbouring cell's element clamp left
      ! unbalanced, which is the case oxygen_chemistry refuses -- is a
      ! transport term and is the same number whatever step produced the
      ! state; the time term is 1/dt times a density and is not.  So a
      ! relative measure taken against the row's FULL terms falls like dt
      ! for one fixed physical defect, and shortening the substep enough
      ! times brings any bounded defect under any fixed floor.  Below, the
      ! same column is integrated at dt and at dt/16, and in each case a
      ! residual of exactly 1e-7 of that row's PHYSICAL terms is put to the
      ! acceptance on each of the two scales in turn.  Sixteen is enough on
      ! this column because two thirds of its full row scale is already the
      ! time term at dt (MEASURED: a 1e-7 physical imbalance reads 3.34e-8
      ! on the full scale at dt and 5.58e-9 at dt/16, so the superseded
      ! measure has crossed the 1e-8 floor while the state's defect has not
      ! moved).
      !
      ! THE SHORT STEP ALSO HAS A FLOOR, and it is the point of the change.
      ! The absolute residual an exact solve can leave is round-off times
      ! the largest term of the row, which at a short step is the time term,
      ! so a substep far below the physical time scale of the row cannot be
      ! certified by any solver: MEASURED on this column at dt/1024, where
      ! the Newton reaches 1e-12 of the full scale and the returned state
      ! still carries 7.7e-5 of its physical terms, and the interval is
      ! refused down to the substep floor.  A retry that shortens the step
      ! therefore cannot buy an acceptance, which is what it could do while
      ! the measure carried the time term.
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
      ! The time term grows with the shortened step and carries the full
      ! row scale with it, which is what makes the superseded measure of an
      ! unchanged defect collapse.
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

      ! The superseded scale: the same physical imbalance is refused at dt
      ! and accepted once the step is short enough.  Executed here rather
      ! than asserted in words, so the defect is a measurement.
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
      ! The physical scale: refused at every step, and by the same number.
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

      ! ---- (7) the shortest substep that can be certified ------------- !
      !
      ! CONVERGENCE AND ACCEPTANCE HAVE TO MEASURE THE SAME THING.  The
      ! acceptance of section (6) is taken against the row's physical
      ! terms.  A Newton that stops when its residual is a fixed fraction
      ! of the FULL terms therefore stops, at a substep far below the row's
      ! physical time scale, at a state the acceptance must refuse: the
      ! full terms are then the time term, so what the iteration left
      ! behind reads 1/dt larger on the scale the state is judged on.
      !
      ! And the arithmetic sets a floor under the measure itself.  The
      ! absolute residual any solve can leave is round-off times the
      ! largest term of the row, which at a short step is the time term, so
      ! |res|/(physical terms) has a floor that grows like 1/dt.  A row
      ! there carries the residual an exact solve of it leaves, which says
      ! nothing about the state: it is accepted and counted as round-off
      ! limited, and the controller does not halve for it.
      !
      ! The three integrations below are the measurement: dt and dt/16 are
      ! covered on the physical measure alone, dt/1024 is the regime in
      ! which the round-off floor bites on this column and the acceptance
      ! rests on the classification.
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
      ! THE INTERVAL AT dt/1024 IS COVERED, and covered in one attempt.
      ! Its rows are at the round-off of their own terms, which is the
      ! residual an exact solve leaves there, so they are no evidence of a
      ! defect in the returned state and refusing them would only send the
      ! controller halving into a floor that grows with every halving.
      call check_absolute('the_short_substep_is_covered',                 &
           logical_as_double(cov_1024), 1.0d0, 0.0d0)
      call check_absolute('the_short_substep_costs_no_halving',           &
           logical_as_double(nhalve_1024 .eq. 0 .and. nsub_1024 .eq. 1),  &
           1.0d0, 0.0d0)
      ! And it is covered BECAUSE of the classification and not because the
      ! physical measure happened to pass: at least one row was accepted as
      ! round-off limited, and the run record carries it.
      call check_absolute('the_short_substep_reports_round_off_limited'// &
           '_rows', logical_as_double(nrl_1024 .ge. 1), 1.0d0, 0.0d0)
      call check_absolute('the_run_record_carries_those_rows',            &
           logical_as_double(nrl_last .eq. nrl_1024 .and.                 &
                             nrl_total .ge. nrl_1024 .and.                &
                             nrl_substeps .ge. 1), 1.0d0, 0.0d0)
      ! The number the classification rests on is measured here rather than
      ! taken on trust: the state IS at the round-off of the full terms.
      call check_absolute('the_short_substep_state_is_at_the_round_off'// &
           '_of_the_full_terms',                                          &
           logical_as_double(rs_1024 .le. carrier_row_roundoff), 1.0d0,   &
           0.0d0)
      ! Nothing is left in the refused-row record of a covered interval.
      call check_absolute('a_covered_interval_reports_no_refused_row',    &
           dble(carrier_refused_rows()), 0.0d0, 0.0d0)
      ! And the state the Newton now returns at the covered steps is one
      ! the acceptance takes on ITS measure, which is what makes the
      ! stopping test and the acceptance the same statement.
      call check_absolute('the_covered_steps_return_an_acceptable'//      &
           '_physical_measure',                                           &
           logical_as_double(rw_dt .le. 1.0d-8 .and. rw_16 .le. 1.0d-8),  &
           1.0d0, 0.0d0)

      ! THE CLASSIFICATION ON CONSTRUCTED ROWS.  The same 1e-7 physical
      ! imbalance of section (6) is given a full-term scale on which its
      ! residual sits at half the round-off bound, which is what an exact
      ! solve of that row leaves.  Every row is then round-off limited, so
      ! the state is accepted and the count says on how many rows.
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
      ! WITHOUT THE FULL TERMS NOTHING IS CLASSIFIED, so the same rows are
      ! refused: the acceptance a caller that builds rows by hand gets is
      ! the one it always got.
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_mixed)
      call check_absolute('without_the_full_terms_the_same_rows_are'//    &
           '_refused',                                                    &
           logical_as_double(vd_mixed%accepted .or.                       &
                             vd_mixed%n_rows_roundoff_limited .ne. 0),    &
           0.0d0, 0.0d0)
      ! ONE row above its own round-off, with a physical measure above the
      ! floor, still refuses: it carries more imbalance than the arithmetic
      ! that formed it explains, and a shorter substep may still reach it.
      tf_round(1,ic_H2) = abs(res_dt(1,ic_H2))
      call carrier_returned_state_verdict(res_dt, tphys_dt, unconstrained,&
           carrier_solve_converged, carrier_stall_at_floor, vd_mixed,     &
           terms_full=tf_round)
      call check_absolute('one_row_above_its_round_off_is_still'//        &
           '_refused', logical_as_double(vd_mixed%accepted), 0.0d0, 0.0d0)
      call check_absolute('and_it_is_the_ordinary_refusal',               &
           dble(vd_mixed%reason), dble(carrier_reject_row_residual),      &
           0.0d0)
      ! The refusal names THAT row, and not one of the round-off limited
      ! ones it stands among.
      call check_absolute('the_refusal_names_the_row_above_its_round'//   &
           '_off',                                                        &
           logical_as_double(vd_mixed%jworst .eq. 1 .and.                 &
                             vd_mixed%icworst .eq. ic_H2 .and.            &
                             vd_mixed%n_rows_above .eq. 1), 1.0d0, 0.0d0)
      call check_relative('the_named_row_keeps_its_physical_measure',     &
           vd_mixed%rworst, 1.0d-7, 1.0d-12)

      ! THE TWO SCALES OF A REFUSED ROW ARE KEPT FOR THE STOP.  The stop
      ! prints, for every row a refusal rests on, the residual and both
      ! scales it was measured against, so a state that is out of balance
      ! can be told from a row the arithmetic cannot resolve without
      ! either being argued from the other.  The record is taken on the
      ! rows of the mixed state above: one row above its own round-off,
      ! the rest of them at it.
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

      ! ---- (8) an exhausted interval, in each run mode ---------------- !
      !
      ! The same forced exhaustion is put to photochemical_transport_step,
      ! which is where the run mode decides what happens next
      ! (a0_run_mode_contract_20260906 sec. 2).  The stop of the physical
      ! branch is suppressed so that both branches can be read in one run;
      ! carrier_transport_stops_suppressed counts the times it was reached,
      ! which is the only statement a test can make about a stop.
      count = 77
      carrier_transport_stop_on_failure = .false.

      ! PHYSICAL INTEGRATION: the status says exhausted and the run ends.
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

      ! INITIALIZATION: nothing is written, the record is kept, the mark is
      ! set, and the march goes on from the entry carriers.
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
      ! The refusal this test forces is the hook's and rests on no row of
      ! the returned state, which the record says by keeping no refused
      ! row; the ledger then names the verdict's own worst row and carries
      ! no scales for it, because the assembly they came from was restored.
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

      ! THE MARCH CONTINUES FROM THE ENTRY CARRIERS.  The next interval is
      ! taken with nothing refused, and the state it reaches is the state
      ! one interval from the entry composition reaches, bit for bit: the
      ! exhausted interval left the carriers exactly where it found them.
      call photochemical_transport_step(rho, v, f_sp, dt_code, st_init)
      call check_absolute('the_next_interval_is_covered',                 &
           dble(st_init), dble(carrier_interval_covered), 0.0d0)
      call check_absolute('and_it_starts_from_the_entry_carriers',        &
           carrier_row_difference(f_sp, f_one), 0.0d0, 0.0d0)

      carrier_transport_stop_on_failure = .true.

      ! ---- (9) a step in which there is nothing to integrate ---------- !
      !
      ! A configuration that carries no carriers has no interval to cover,
      ! which is not a failure: the operator returns at once.  Two things
      ! are asserted about that return, and the interval just covered above
      ! is what makes them worth asserting.
      !
      !  * THE STATUS IS DEFINED.  st_nothing is set to a value that is
      !    neither carrier_interval_covered nor carrier_interval_exhausted,
      !    so an operator that left its optional argument untouched would
      !    be read here and not covered up by a caller's pre-set.
      !  * THE STEP VERDICT IS "NO INTERVAL", and not the verdict of the
      !    last interval of the run, which is still standing and is read
      !    from the run-level record beside it.
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

      ! The same, through the other gate: the molecular layer is there and
      ! the transport option is off.
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

      ! THE CO DOMAIN RECORD IS INDEXED BY THE LEDGER FAMILY, and the two
      ! readings compose by two different rules.  A cell first out of
      ! domain while the run was relaxing is not a cell the physical
      ! history touched, so family 1 (initialization and continuation) and
      ! family 2 (accepted physical steps) are kept apart.  The COUNTS of
      ! the run are their sum, because they are cell visits and a visit
      ! belongs to one family; the WORST RATIO of the run is their maximum,
      ! because it is an extreme value and not a tally.
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
      ! The record is a run-level statement and belongs to no checkpoint:
      ! taking one and restoring it leaves the record where it stands.
      call carrier_checkpoint_take(chk_after)
      call carrier_checkpoint_restore(chk_after)
      call carrier_co_domain_record(dout_f1, dhot_f1, dhep_f1,            &
                                    drat_f1, dr_f1, dform_f1)
      call check_absolute('the_domain_record_is_not_a_checkpointed'//     &
           '_quantity', dble(dout_f1), dble(dout_tot), 0.0d0)


      ! ---- (10) the movement bound of a relaxation pass -------------- !
      !
      ! A RELAXATION PASS ADVERTISES A BOUND ON THE STATE IT HANDS BACK:
      ! the largest change of any solved carrier column anywhere on the
      ! grid, divided by the largest H2 mixing ratio of the entry state, is
      ! at most `trust`.  A displacement tested only after the step has
      ! been applied bounds nothing, which is what was MEASURED in
      ! docs/solver_partition_experiment_20260911.md sec. 7.2 (trust 0.01
      ! returning 0.0281).  The rows below are that statement on this
      ! column, on the same measure the routine enforces and reports.
      !
      ! THE WIND OF THIS SECTION IS ITS OWN.  The pass builds its step from
      ! the cell crossing time of the wind it is handed, and it rides the
      ! face mass flux of the mass row of that state, so the column is
      ! given an outflow whose crossing time is the interval the earlier
      ! sections measured this operator on, and its mass row is assembled
      ! before any pass runs.
      allocate(nnuc0(1-Ng:N+Ng,n_element), nnuc1(1-Ng:N+Ng,n_element))
      allocate(nchg0(1-Ng:N+Ng), nchg1(1-Ng:N+Ng))
      allocate(rcomp0(1-Ng:N+Ng), rcomp1(1-Ng:N+Ng))
      allocate(v_relax(1-Ng:N+Ng))
      v_relax = 20.0d0*r
      run_mode = run_mode_phys
      carrier_transport_stop_on_failure = .false.
      call element_nuclei_and_charge(rho, f_sp0, nnuc0, nchg0, rcomp0)
      call seed_mass_row_of_the_column(v_relax)

      ! THE LADDER, MEASURED ON THIS COLUMN.  A trial at the full cell
      ! crossing time moves the composition by 0.441 of the largest entry
      ! H2 mixing ratio, and halving the trial halves the movement: 0.252,
      ! 0.144, 0.0779, 0.0404, 0.0205, 0.0103, 5.16e-3, 2.58e-3, 1.29e-3
      ! and 6.46e-4 at the tenth halving, which is relax_grow_min (all
      ! MEASURED).  So a bound anywhere in that range is met by a shortened
      ! trial, and the rows below read the pass at three of them.
      !
      ! A BOUND BELOW WHAT A FULL TRIAL WOULD MOVE IS STILL A BOUND ON THE
      ! RETURNED STATE, and it is the bound that ends the pass.
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
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
      ! The write-back of every covered step restores the entry element
      ! totals cell by cell, so a shortened trial cannot move a nucleus
      ! count either.
      call element_nuclei_and_charge(rho, f_sp, nnuc1, nchg1, rcomp1)
      call check_absolute('the_bounded_pass_conserves_every_element',     &
           element_deviation(nnuc0, nnuc1), 0.0d0, 1.0d-12)

      ! THE BOUND SCALES.  A tenfold tighter setting holds too, and holds
      ! at a smaller displacement, which is what an outer controller needs
      ! of it: reducing the setting reduces the update rather than leaving
      ! it where it was.
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
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

      ! A BOUND BELOW WHAT THE SHORTEST ADMISSIBLE TRIAL MOVES IS STILL
      ! HONORED, by handing back the state that is inside it: the entry
      ! state.  On this column the shortest trial moves 6.46e-4, so at
      ! 1e-4 no trial fits and nothing is kept.  That is the bound doing
      ! its job and not a failure of the pass; what it costs is the
      ! advance, which the ending and the zero drift both say.
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
                                           eta_col, 1.0d-4, dr_none,      &
                                           ns_none, oc_none)
      call check_absolute('a_bound_no_trial_can_meet_keeps_no_step',      &
           dble(ns_none), 0.0d0, 0.0d0)
      call check_absolute('and_hands_back_the_entry_composition',         &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('and_reports_the_bound_as_what_ended_it',       &
           dble(oc_none), dble(carrier_relax_movement_bound), 0.0d0)

      ! A FIRST-STEP OVERSHOOT IS CAUGHT AND NOT KEPT.  At a bound of zero
      ! every step of this operator leaves it, the shortest admissible
      ! trial included, so nothing is kept: the entry composition comes
      ! back bit for bit and the ending names the bound.  This is the row
      ! the defect fails: the first step of the operator moves the
      ! composition, and a test taken after that step is applied cannot
      ! take it back.
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
                                           eta_col, 0.0d0, dr_zero,       &
                                           ns_zero, oc_zero)
      call check_absolute('a_pass_that_may_not_move_keeps_no_step',       &
           dble(ns_zero), 0.0d0, 0.0d0)
      call check_absolute('and_returns_the_entry_composition_bit_for'//   &
           '_bit', maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('and_reports_no_drift', dr_zero, 0.0d0, 0.0d0)
      call check_absolute('and_names_the_bound_as_what_ended_it',         &
           dble(oc_zero), dble(carrier_relax_movement_bound), 0.0d0)

      ! A TRANSPORT INTERVAL THAT CANNOT BE COVERED AT ANY ADMISSIBLE
      ! LENGTH IS ITS OWN ENDING, and it leaves the entry state.  The
      ! refusal is injected through the test hook, and the run mode is the
      ! one that does not stop inside the operator so that the pass can be
      ! read here.
      allocate(v_fast(1-Ng:N+Ng))
      v_fast = 2000.0d0*r
      call seed_mass_row_of_the_column(v_fast)
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      run_mode = run_mode_init
      carrier_reject_leading_attempts_for_test = 1000000
      call restore_frozen_background()
      call conserved_of(v_fast)
      call relax_photochemical_composition(u_col, v_fast, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
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

      ! THE FIXED-POINT ENDING IS NOT REACHABLE ON THIS COLUMN, and the
      ! ladder above says why: the movement of a trial is proportional to
      ! its length, and the shortest admissible trial still moves 6.46e-4
      ! of the largest entry H2 mixing ratio, which is six decades above
      ! relax_tol = 1e-10.  No trial of this operator on this column leaves
      ! the state where it found it, so no pass of it can report a fixed
      ! point.  That is the propagating H2 front the routine's own header
      ! records, not a defect of the ending; the rows above pin every
      ! ending this fixture can reach.

      ! A CONFIGURATION THAT TRANSPORTS NO CARRIER HAS NOTHING TO ADVANCE,
      ! and that is its own ending and neither a bound nor a fixed point.
      carrier_transport = .false.
      f_sp = f_sp0
      call restore_frozen_background()
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
                                           eta_col, 1.0d0, dr_zero,       &
                                           ns_zero, oc_zero)
      carrier_transport = .true.
      carrier_transport_stop_on_failure = .true.
      call check_absolute('a_pass_with_no_transported_carrier_says_so',   &
           dble(oc_zero), dble(carrier_relax_nothing_to_advance), 0.0d0)
      call check_absolute('and_moves_nothing', maxval(abs(f_sp - f_sp0)), &
           0.0d0, 0.0d0)


      ! ---- (11) the chemistry the pass carries with it --------------- !
      !
      ! A RELAXATION PASS ADVANCES TRANSPORT AND CHEMISTRY TOGETHER.  The
      ! transport operator integrates the carrier rows on the background of
      ! ONE equilibrium sweep (bg_cell: the rates, the temperature, the
      ! third-body density and the eliminated species), so a pass that
      ! keeps many steps on one background converges to the fixed point of
      ! a different problem.  MEASURED on the hot-Uranus carrier reload
      ! (READ, docs/solver_partition_experiment_20260911.md sec. 4): a full
      ! frozen-background relaxation drove the gated H2 wind row to
      ! 2.449021e-13 and the chemical refresh that followed put the same
      ! row at 7.989781e-2, above the 7.380744e-2 the pass had started
      ! from.  The rows below state the repair on this column: after a
      ! pass, the background the NEXT transport step would read describes
      ! the composition the pass hands back, not the composition it began
      ! with.
      !
      ! THE JOINT FIXED POINT CANNOT BE ASSERTED ON THIS COLUMN, and the
      ! ladder above says why: the shortest admissible trial still moves
      ! 6.46e-4 of the largest entry H2 mixing ratio, six decades above
      ! relax_tol = 1e-10, so no pass of this operator on this column can
      ! report carrier_relax_fixed_point.  What is asserted instead is the
      ! consistency the joint fixed point rests on, on the two quantities
      ! of the background that the carrier rows read directly, and by one
      ! further sweep taken on the returned state.
      allocate(f_swept(1-Ng:N+Ng,n_species))
      allocate(nhi_r(1-Ng:N+Ng), nhii_r(1-Ng:N+Ng))
      allocate(nhei_r(1-Ng:N+Ng), nheii_r(1-Ng:N+Ng))
      allocate(nheiii_r(1-Ng:N+Ng), nheiTR_r(1-Ng:N+Ng))
      allocate(ne_r(1-Ng:N+Ng), ntot_r(1-Ng:N+Ng), T_ret(1-Ng:N+Ng))
      allocate(nm_r(1-Ng:N+Ng,n_mion))

      call seed_mass_row_of_the_column(v_relax)
      call seed_background_at_the_entry_composition()
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,     &
                                           T_col, heat_col, cool_col,    &
                                           eta_col, 1.0d-2, dr_wide,      &
                                           ns_wide, oc_wide)
      call check_positive('the_consistency_rows_read_a_pass_that_kept'//  &
           '_a_step', dble(ns_wide))
      write(*,'(A,A,A,I0,A,I0)') '  DIAGNOSTIC closure of the section-11'// &
           ' pass: ', trim(chem_closure_reason_text(n_chem_last_reason)),   &
           ', off-simplex cells ', chem_last_offsimplex, ', steps kept ', ns_wide

      ! The densities and the temperature of the composition the pass
      ! handed back, formed here from that composition and the fixed
      ! pressure, independently of anything the routine wrote.
      nhei_r   = 0.0d0
      nheii_r  = 0.0d0
      nheiii_r = 0.0d0
      nheiTR_r = 0.0d0
      call get_species_densities(rho, f_sp, nhi_r, nhii_r, nhei_r,        &
                                 nheii_r, nheiii_r, nheiTR_r, nm_r,       &
                                 ne_r, ntot_r)
      call comp_T_from_p(p_col, ntot_r, ne_r, T_ret)
      ! THE BACKGROUND DESCRIBES THE COMPOSITION HANDED BACK, absolutely:
      ! the third-body density and the temperature the next transport step
      ! reads are those of the returned composition to the tolerance the
      ! sweep itself converges to (ieq_res_tol), not merely closer to it
      ! than to the entry.
      dev_ntot = 0.0d0
      dev_TK   = 0.0d0
      do jj = 1, N
         dev_ntot = max(dev_ntot, abs(bg_cell(jj)%ntot/n0 - ntot_r(jj))   &
                                  /ntot_r(jj))
         dev_TK   = max(dev_TK, abs(bg_cell(jj)%T_K/T0 - T_col(jj))       &
                                /T_col(jj))
      enddo
      write(*,'(a,2es12.4)') ' (carrier_retry) background third-body'//   &
           ' density and temperature against the returned state: ',      &
           dev_ntot, dev_TK
      call check_absolute('the_background_third_body_density_is_that'//   &
           '_of_the_returned_composition', dev_ntot, 0.0d0,               &
           1.0d1*ieq_res_tol)
      call check_absolute('the_background_temperature_is_that_of_the'//   &
           '_returned_composition', dev_TK, 0.0d0, 1.0d1*ieq_res_tol)
      ! THE EQUATION OF STATE ROUND TRIP AT THE FIXED CONSERVED STATE: the
      ! pressure handed back is the one the caloric equation of state of
      ! the returned composition assigns to the unchanged thermal energy,
      ! and the temperature is p/(n_tot + n_e) of that composition.
      dev_p = 0.0d0
      dev_T = 0.0d0
      do jj = 1, N
         dev_p = max(dev_p, abs(p_col(jj)                                 &
                    - pressure_from_energy_density(jj, u_col(1,jj),       &
                        u_col(3,jj) - 0.5d0*u_col(2,jj)**2/u_col(1,jj)))  &
                    /p_col(jj))
         dev_T = max(dev_T, abs(T_col(jj) - T_ret(jj))/T_ret(jj))
      enddo
      write(*,'(a,2es12.4)') ' (carrier_retry) equation of state round'// &
           ' trip at the fixed conserved state, p and T: ', dev_p, dev_T
      call check_absolute('the_returned_pressure_is_that_of_the'//        &
           '_unchanged_thermal_energy', dev_p, 0.0d0, 1.0d-13)
      call check_absolute('and_the_returned_temperature_is_that_of'//     &
           '_the_returned_composition', dev_T, 0.0d0, 1.0d-13)
      ! THE CHEMICAL RESIDUAL OF THE RETURNED STATE: a further sweep at
      ! the returned temperature moves the composition by no more than the
      ! sweep's own convergence tolerance, relative to each species that
      ! is present, which is what a closed chemical state means.
      move_pass = maxval(abs(f_sp(1:N,:) - f_sp0(1:N,:)))
      f_swept   = f_sp
      call ioniz_eq(T_col, rho, f_swept, heat_col, cool_col, eta_col)
      move_rel = 0.0d0
      do jj = 1, N
         do isp = 1, n_species
            fref = f_sp(jj,isp)
            if (fref .gt. 1.0d-20)                                        &
               move_rel = max(move_rel, abs(f_swept(jj,isp) - fref)/fref)
         enddo
      enddo
      write(*,'(a,es12.4,a,es12.4)') ' (carrier_retry) the pass moved'//   &
           ' the composition by ', move_pass, ' and one further sweep'//  &
           ' moves it by (relative, present species) ', move_rel
      call check_absolute('one_further_sweep_leaves_the_returned'//       &
           '_composition_within_the_sweep_tolerance', move_rel, 0.0d0,   &
           ieq_res_tol)

      ! ---- (12) a rejected stationary trial leaves the history alone --- !
      !
      ! A trial of the relaxation is a step the caller keeps only if it
      ! likes the result.  One that the operator cannot cover is undone,
      ! composition and background, and is not an interval the run adopted
      ! without integrating it: the carrier history the certification
      ! reads stays as it was, whatever the run mode.  An adopted physical
      ! interval that fails still marks it (section (2) above).
      call carrier_history_reset_for_test()
      run_mode_held = run_mode
      run_mode = run_mode_init
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      carrier_reject_leading_attempts_for_test = 1000000
      hist_before = carrier_history_certifiable()
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,   &
                                           T_col, heat_col, cool_col,     &
                                           eta_col, 1.0d-2, dr_hist,      &
                                           ns_hist, oc_hist)
      hist_after = carrier_history_certifiable()
      carrier_reject_leading_attempts_for_test = 0
      write(*,'(a,l1,a,l1,a,i0,a,a)') ' (carrier_retry) history before'//  &
           ' the refused trials ', hist_before, ', after ', hist_after,    &
           ', kept steps ', ns_hist, ', ending on '//                     &
           trim(carrier_relax_outcome_text(oc_hist))
      call check_absolute('a_rejected_stationary_trial_keeps_the'//       &
           '_history_certifiable',                                        &
           logical_as_double(hist_before .and. hist_after), 1.0d0, 0.0d0)
      call check_absolute('and_keeps_no_step', dble(ns_hist), 0.0d0, 0.0d0)
      call check_absolute('and_hands_back_the_entry_composition_bit'//    &
           '_for_bit', maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      call check_absolute('and_names_the_refused_interval',               &
           dble(oc_hist), dble(carrier_relax_interval_refused), 0.0d0)
      ! The pass that must keep a step runs on the background seeded at the
      ! entry composition (as the consistency rows of section 11 do), not
      ! on the frozen background of the checkpoint rows: on that stale
      ! background the sweep leaves a cell off the element simplex and the
      ! closure contract of 2026-09-13 (S2) rightly refuses the step, which
      ! is the chemistry's verdict and not the history's.
      call seed_mass_row_of_the_column(v_relax)
      call seed_background_at_the_entry_composition()
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,   &
                                           T_col, heat_col, cool_col,     &
                                           eta_col, 1.0d-2, dr_hist,      &
                                           ns_hist, oc_hist)
      write(*,'(A,A,A,ES10.3,A,I0)') '  DIAGNOSTIC closure of the pass after'// &
           ' the refusals: ', trim(chem_closure_reason_text(n_chem_last_reason)), &
           ', last increment ', chem_last_increment, ', steps kept ', ns_hist
      write(*,'(A,I0,A,I0,A,ES10.3)') '  DIAGNOSTIC off-simplex cells ',    &
           chem_last_offsimplex, ' of which molecular clamps ',             &
           chem_last_mol_clamped, ', worst budget violation ', chem_last_viol_worst
      ! What the refused trials must not have done is refuse THIS pass at
      ! the interval or mark the history: the pass is judged by the
      ! chemistry on its own (on this synthetic column the closure of the
      ! first trial leaves one cell off the element simplex, a verdict of
      ! the chemistry contract of 2026-09-13, S2, and not of the history;
      ! the pass of section 11, prepared on the background of the
      ! composition it relaxes from, keeps its steps).
      call check_absolute('a_pass_after_the_refusals_is_not_refused_at'// &
           '_the_interval', logical_as_double(oc_hist .ne.                &
           carrier_relax_interval_refused), 1.0d0, 0.0d0)
      call check_absolute('and_stays_eligible',                           &
           logical_as_double(carrier_history_certifiable()), 1.0d0, 0.0d0)
      run_mode = run_mode_phys
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      carrier_reject_leading_attempts_for_test = 1000000
      call photochemical_transport_step(rho, v_relax, f_sp, dt_code,      &
                                        st_trial, trial = .true.)
      carrier_reject_leading_attempts_for_test = 0
      call check_absolute('a_refused_trial_in_phys_mode_reports'//        &
           '_exhausted', dble(st_trial), dble(carrier_interval_exhausted), &
           0.0d0)
      call check_absolute('and_leaves_the_history_certifiable',           &
           logical_as_double(carrier_history_certifiable()), 1.0d0, 0.0d0)
      call check_absolute('and_writes_no_composition',                    &
           maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      run_mode = run_mode_held

      ! ---- 13. The thermochemical closure contract (S2 of PLAN_20260913,
      ! findings B3/B4): a closure that spends its cycle budget is NOT
      ! reported as reached, the trial is refused with the chemistry named,
      ! nothing is kept, and the entry composition comes back bit for bit.
      ! The budget is forced to zero through the test knob; RED on the text
      ! before 2026-09-13, where ok started true and exhaustion fell
      ! through with it (the relaxation then kept the step).
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call restore_frozen_background()
      call conserved_of(v_relax)
      chem_cycles_cap_for_test = 0
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,   &
                                           T_col, heat_col, cool_col,     &
                                           eta_col, 1.0d-2, dr_hist,      &
                                           ns_hist, oc_hist)
      chem_cycles_cap_for_test = -1
      call check_absolute('an_exhausted_closure_refuses_the_trial',       &
           dble(oc_hist), dble(carrier_relax_chemistry_refused), 0.0d0)
      call check_absolute('and_names_the_spent_budget',                   &
           dble(n_chem_last_reason), dble(chem_closure_exhausted), 0.0d0)
      call check_absolute('and_keeps_no_step', dble(ns_hist), 0.0d0, 0.0d0)
      call check_absolute('and_hands_back_the_entry_composition_bit'//    &
           '_for_bit', maxval(abs(f_sp - f_sp0)), 0.0d0, 0.0d0)
      write(*,'(A,A)') '  DIAGNOSTIC closure reason text: ',              &
           trim(chem_closure_reason_text(n_chem_last_reason))

      if (assertion_failures .gt. 0) then
         write(*,'(a)') 'carrier_retry: FAILED'
         stop 1
      endif
      write(*,'(a)') 'carrier_retry: every assertion passed'

      contains

      !--------------!

      function roundoff_of(res, tfull) result(flag)
      ! Which of these rows carry a residual at or below the round-off of
      ! their own full terms, which is what carrier_returned_state_verdict
      ! classifies and accepts.
      real*8, intent(in) :: res(:,:), tfull(:,:)
      logical :: flag(size(res,1),size(res,2))
      flag = (abs(res) .le. carrier_row_roundoff*tfull)
      end function roundoff_of

      !--------------!

      double precision function carrier_row_difference(a, b) result(d)
      ! The largest difference between two composition states over the
      ! transported columns, relative to the entry state of each column.
      ! One number for the five carriers and the atomic hydrogen they
      ! exchange with.
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

      !--------------!

      double precision function element_deviation(a, b) result(d)
      ! The largest relative difference between two elemental nucleus
      ! censuses over the physical cells.  Repartitioning the carriers at a
      ! fixed density moves no nucleus, so this is zero to round-off for
      ! every state a relaxation pass may hand back.
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

      !--------------!

      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double

      !--------------!

      subroutine setup_globals()
      ! A twelve-cell hydrogen and helium column with the molecular
      ! chemistry and the transported proton on, and no metals in the gas.
      ! The gravity is spherical so that Dphi is b0/r^2 and the settling
      ! term of the carrier faces is defined.
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
      ! Every grid-sized array of global_parameters at once, by the one
      ! routine the run itself uses, because the equilibrium sweep this
      ! column now calls reads several of them (opa_pf among the rest).
      call allocate_grid_arrays
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.02d0*dble(j)
         r_edg(j) = 1.0d0 + 0.02d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.02d0
      enddo
      kzz_cell = 0.0d0
      call set_gravity_grid
      ! The hydrodynamic keys the mass row of this column is assembled with
      ! (seed_mass_row_of_the_column): one reconstruction, one Riemann
      ! solve, and the domain the ghost fill closes on.
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
      ! THE PHOTON GRID AND THE COUPLED SYSTEM THE SWEEP SOLVES.  The
      ! relaxation pass equilibrates the chemistry on every step it keeps,
      ! so this column needs the spectrum the rates are formed from and the
      ! unknown count of the molecular system.  The spectrum is the
      ! power law of the hot-Uranus fixture, which fixes nothing tested
      ! here: what the rows below read is the carrier movement and the
      ! consistency of the background with the composition, not a rate.
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
      ! The molecular system: H+/He+/He++ and the four molecular carriers of
      ! the hydrogen nuclei, with two stages of each metal element above
      ! them (input_read sets the same count for this configuration).
      N_eq = 7 + 2*n_melem
      lwa  = (N_eq*(3*N_eq + 13))/2
      allocate(sys_sol(N_eq), sys_x(N_eq), wa(lwa))
      call cx_init
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      bg_ready = .true.
      end subroutine setup_globals

      !--------------!

      subroutine build_molecular_hydrogen_column()
      ! Molecular base gas: 80 per cent of the hydrogen nuclei in H2, the
      ! rest atomic, helium at the solar-like ratio, a trace of the
      ! molecular ions, no oxygen or carbon anywhere, and an outflowing
      ! velocity that grows outward so the advective term is not uniform.
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
      ! A column whose species carry the density they are measured
      ! against (test_columns): 80 percent of the hydrogen nuclei in H2,
      ! helium at the reservoir mass fraction everywhere.  The trace ions
      ! are then moved out of their own element's neutral so that nothing
      ! is added on top of a closed mass budget.
      call column_carrying_its_own_density(f_sp, 0.8d0, .false., 0.0d0)
      f_sp(:,isp_HII)  = 1.0d-6
      f_sp(:,isp_HI)   = f_sp(:,isp_HI) - 1.0d-6
      f_sp(:,isp_H2p)  = 1.0d-10
      f_sp(:,isp_H3p)  = 1.0d-10
      f_sp(:,isp_H2)   = f_sp(:,isp_H2) - 2.5d-10
      f_sp(:,isp_HeII) = 1.0d-8
      f_sp(:,isp_HeI)  = f_sp(:,isp_HeI) - 1.0d-8
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      f_sp0 = f_sp
      allocate(p_col(1-Ng:N+Ng), T_col(1-Ng:N+Ng))
      allocate(heat_col(1-Ng:N+Ng), cool_col(1-Ng:N+Ng))
      allocate(eta_col(1-Ng:N+Ng))
      allocate(u_col(3,1-Ng:N+Ng), W_col(3,1-Ng:N+Ng))
      call seed_pressure_of_the_column()
      end subroutine build_molecular_hydrogen_column

      !--------------!

      subroutine seed_pressure_of_the_column()
      ! THE PRESSURE THE RELAXATION HOLDS FIXED.  A pass is taken at a wind
      ! that does not move, and the temperature of every composition it
      ! visits is p/(n_tot + n_e) at that pressure.  The pressure is
      ! therefore set once, from the entry composition and the temperature
      ! of the frozen cell state, so that the entry state is its own:
      ! recomputing T at the entry composition returns exactly the
      ! temperature the frozen rates were formed at, and what moves T later
      ! is the particle count the carriers change.
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

      ! The conserved state (rho, rho v, E) of the entry composition at
      ! the wind vv and the entry pressure, by the code's own equation of
      ! state evaluated on the ENTRY composition (the caloric state is the
      ! one get_species_densities last refreshed, so it is refreshed here
      ! on f_sp0 first).
      subroutine conserved_of(vv)
      real*8, dimension(1-Ng:N+Ng), intent(in) :: vv
      real*8, allocatable :: a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)
      real*8, allocatable :: a7(:), a8(:), am(:,:)
      allocate(a1(1-Ng:N+Ng), a2(1-Ng:N+Ng), a3(1-Ng:N+Ng))
      allocate(a4(1-Ng:N+Ng), a5(1-Ng:N+Ng), a6(1-Ng:N+Ng))
      allocate(a7(1-Ng:N+Ng), a8(1-Ng:N+Ng), am(1-Ng:N+Ng,n_mion))
      a3 = 0.0d0;  a4 = 0.0d0;  a5 = 0.0d0;  a6 = 0.0d0
      call get_species_densities(rho, f_sp0, a1, a2, a3, a4, a5, a6, am,  &
                                 a7, a8)
      W_col(1,:) = rho
      W_col(2,:) = vv
      W_col(3,:) = p_col
      call W_to_U(W_col, u_col)
      end subroutine conserved_of

      !--------------!

      subroutine seed_background_at_the_entry_composition()
      ! THE FROZEN CELL STATE OF SECTION (11), SEEDED AT THE COLUMN'S OWN
      ! COMPOSITION.  The rows of that section read how far the background
      ! stands from the entry composition and from the returned one, so the
      ! two quantities of it that the carrier rows read directly, the
      ! third-body density and the temperature, start at the values the
      ! entry composition actually has.  The rate coefficients are the
      ! section's usual frozen ones.
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

      !--------------!

      subroutine restore_frozen_background()
      ! The cell state every pass of section (10) is entered with.  A pass
      ! equilibrates the chemistry on each step it keeps, so it leaves
      ! bg_cell describing the composition it handed back; the rows below
      ! compare passes with one another, and each has to start from the
      ! same background and the same temperature.
      integer :: j
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
         T_col(j) = bg_cell(j)%T_K/T0
      enddo
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine restore_frozen_background

      !--------------!

      subroutine seed_mass_row_of_the_column(v_wind)
      ! THE FACE MASS FLUX THE CARRIER ROWS RIDE ON.  A relaxation pass at a
      ! fixed wind carries the advective term in its rows, and that term is
      ! built from the face mass flux of the MASS ROW of the state it was
      ! handed (steady_residual_mod, face_mass_flux_of_state).  So the
      ! column needs its mass row assembled before any pass can run on it,
      ! and the state assembled here is the column's own density and
      ! velocity with an ideal pressure at its frozen temperature.
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

      !--------------!

      subroutine set_frozen_cell_rates(j)
      ! Every field of the frozen cell state, set here so that none is read
      ! undefined.  The photoionization rate and the case-B recombination
      ! coefficient are the two that make the proton row depend on n(H+):
      ! production from H I and the recombination sink alpha n_e n_H+.
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

      end program carrier_retry

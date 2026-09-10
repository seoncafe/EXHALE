      program attempted_step_tests
      ! THE CHECKPOINT AND THE STEP-SIZE POLICY OF THE ATTEMPTED-STEP
      ! CONTROLLER (PLAN_20260906_rev2 step B3a).
      !
      ! WHAT IS TESTED HERE, and what is tested on a whole binary by
      ! run.sh beside this driver.
      !
      !  (0) THE UNSWEPT BACKGROUND IS THE ZERO RATE STATE. bg_cell and
      !      ieq_rate_cell are hashed, saved and restored by the checkpoint
      !      from the first step on, so what allocation leaves in them is
      !      part of the state; a cell the sweep has not reached carries no
      !      rates and no imposed carrier partition.
      !  (1) THE CHECKPOINT IS COMPLETE. Every item of the list is given a
      !      distinct value by
      !      attempted_step_perturb_checkpointed_state_for_test, the
      !      checkpoint is restored, and the result has to come back bit
      !      for bit: item by item through
      !      attempted_step_checkpoint_matches, and as one number through
      !      attempted_step_checkpoint_checksum, which also hashes the
      !      quantities operation 3 REBUILDS and the checkpoint therefore
      !      does not save. An item missing from either half of the pair
      !      fails here.
      !  (2) THE THREE COUNTER CLASSES. A restore puts the physical
      !      accumulations back (the clock, the accepted-step count) and
      !      leaves the attempt statistics and the diagnostic extrema where
      !      the refused attempt put them (b1 T7.2, decision 10).
      !  (3) THE STEP-SIZE POLICY. The halving of an admissibility refusal,
      !      the controller factor of an integration-error refusal, the
      !      floor at 0.2, and the bound a reduced step puts on the next
      !      one (advisor decision 3).
      !  (4) THE INTEGRATION-ERROR ESTIMATE. Its measure on two states with
      !      a known difference, the RETAINED-STEP factor 1/(1 - 2^-p) it
      !      carries because the trajectory the run keeps is the full step,
      !      and the ORDER it reports on a smooth pair: halving the
      !      difference halves the estimate.
      !  (5) THE FORMATION-ENERGY TABLE of the energy reservoir: each entry
      !      against the constant it is built from, so that a change to one
      !      of those constants cannot silently change the reservoir. The
      !      table has one home, molecular_reaction_heat, and this is where
      !      its entries are checked against their sources.
      !  (6) THE THERMAL IDENTITY OF THE SOURCE STEP. A synthetic step is
      !      given a known thermal source, a known composition change and a
      !      known transport change, in that order and each separately
      !      measurable. The instrument has to return the closure of the
      !      SOURCE STEP alone at round-off, the transport contribution as
      !      the number that was put in, and the reservoir change reported
      !      and NOT added. The same block records the RED reference: the
      !      quantity the instrument used to report, d_u_th + d_u_form -
      !      Q_ext_dt, is far from zero on exactly this step. The gate is
      !      then shown to refuse an injected inconsistency.
      !
      ! Every assertion prints PASS|FAIL name measured= reference= tol= and
      ! the program exits nonzero if any fails.

      use global_parameters
      use species_table
      use ionization_equilibrium, only: ioniz_eq_allocate_arrays,        &
                                        bg_ready, ieq_marching_ledger,   &
                                        ieq_nonroot_streak, bg_cell,     &
                                        ieq_rate_cell
      use ion_cell_state, only: ion_rates
      use attempted_step
      use molecular_reaction_heat, only: species_formation_energy
      use utils_ion_eq, only: heat_channel_state
      use assertion_report, only: check_relative, check_absolute,        &
                                  assertion_failures
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none

      integer, parameter :: ncell = 12
      real*8, allocatable :: u(:,:), f_sp(:,:), heat(:), cool(:), eta(:)
      real*8, allocatable :: u_a(:,:), u_b(:,:)
      ! The species, temperature and electron rows of the estimate.
      real*8, allocatable :: f_a(:,:), f_b(:,:)
      real*8, allocatable :: T_a(:), T_b(:), ne_a(:), ne_b(:)
      ! The coupled pair's remaining error at each cell, as the marching
      ! loop hands it over.
      real*8, allocatable :: inner_c(:), inner_T(:)
      type(integration_error) :: est
      ! The synthetic attempted step of block (6).
      real*8, allocatable :: W(:,:), Tk(:), dt_cell(:)
      type(attempted_step_verdict) :: verd
      real*8  :: du_transport = 0.0d0, q_expected = 0.0d0
      ! The heating channel decomposition as the checkpoint was taken with
      ! it, for the row that moves that array alone.
      real*8, allocatable :: hcs0(:,:)
      type(attempted_step_checkpoint) :: chk
      ! The state of a cell that carries no gas and no rates, against which
      ! the freshly allocated background is measured in block (0). It is
      ! SAVEd, so it is the zero state whether or not ion_rates carries
      ! default initializers.
      type(ion_rates), save :: unswept
      real*8  :: h0, h1, h2, dt_a, dt_b, e1, e2
      real*8  :: t_phys0
      integer :: n_acc0, jw, kw
      logical :: ok
      character(len=48) :: bad

      call setup_globals()
      call build_state()

      ! ---- (0) the unswept background is the zero rate state ---------- !
      ! bg_cell and ieq_rate_cell are hashed, saved and restored by the
      ! checkpoint from the first step on, so what allocation leaves in them
      ! is part of the state the controller carries. The invariant is that a
      ! cell the sweep has not reached carries no rates: no photoionization,
      ! no recombination, no third body, and no imposed carrier partition.
      ! Overwriting them with that state must not move the checksum.
      !
      ! setup_globals dirties the free lists before it allocates, so the
      ! allocation is served from a recycled block; without that, the row
      ! would measure the zero pages a fresh heap hands out rather than the
      ! initialization of ion_rates, and would pass either way.
      h0 = attempted_step_checkpoint_checksum(u, f_sp, heat, cool, eta)
      do jw = 1-Ng, N+Ng
         bg_cell(jw)       = unswept
         ieq_rate_cell(jw) = unswept
      enddo
      h1 = attempted_step_checkpoint_checksum(u, f_sp, heat, cool, eta)
      call check_absolute('the_unswept_background_carries_no_rates',      &
           h1 - h0, 0.0d0, 0.0d0)

      ! ---- (1) the checkpoint round trip ------------------------------ !
      t_phys           = 12.5d0
      n_steps_accepted = 7
      n_steps_attempted = 9
      call attempted_step_checkpoint_take(chk, u, f_sp, heat, cool, eta,  &
                                          0, 3, 2)
      h0 = attempted_step_checkpoint_checksum(u, f_sp, heat, cool, eta)

      call attempted_step_perturb_checkpointed_state_for_test(u, f_sp,    &
                                                       heat, cool, eta)
      h1 = attempted_step_checkpoint_checksum(u, f_sp, heat, cool, eta)
      call check_absolute('a_perturbed_state_is_not_the_checkpoint',      &
           logical_as_double(h1 .ne. h0), 1.0d0, 0.0d0)
      ok = attempted_step_checkpoint_matches(chk, u, f_sp, heat, cool,    &
                                             eta, bad)
      call check_absolute('the_perturbation_is_seen_item_by_item',        &
           logical_as_double(ok), 0.0d0, 0.0d0)

      ! The attempt statistics and the diagnostic extrema are moved on
      ! before the restore, so the assertion below is that they did NOT
      ! come back.
      n_steps_attempted = n_steps_attempted + 5

      call attempted_step_checkpoint_restore(chk, u, f_sp, heat, cool,    &
                                             eta)
      h2 = attempted_step_checkpoint_checksum(u, f_sp, heat, cool, eta)
      call check_absolute('the_state_is_restored_bit_for_bit',            &
           h2 - h0, 0.0d0, 0.0d0)
      ok = attempted_step_checkpoint_matches(chk, u, f_sp, heat, cool,    &
                                             eta, bad)
      call check_absolute('every_item_is_restored',                       &
           logical_as_double(ok), 1.0d0, 0.0d0)

      ! ---- (1b) the heating channel decomposition, on its own --------- !
      ! heat_channel_state is the definition the heat column and
      ! Heating_breakdown.txt are both built from, so a trial refused after
      ! the sweep must not leave it describing the refused state. It is
      ! asserted alone and not only through the whole-list perturbation
      ! above, because there an item that is NOT compared is masked by
      ! every item that is.
      call attempted_step_checkpoint_take(chk, u, f_sp, heat, cool, eta,  &
                                          0, 3, 2)
      hcs0 = heat_channel_state
      heat_channel_state = heat_channel_state + 3.0d0
      bad = ''
      ok = attempted_step_checkpoint_matches(chk, u, f_sp, heat, cool,    &
                                             eta, bad)
      call check_absolute('a_moved_heating_channel_array_is_seen',        &
           logical_as_double((.not. ok) .and.                             &
                             trim(bad) .eq. 'heat_channel_state'),        &
           1.0d0, 0.0d0)
      call attempted_step_checkpoint_restore(chk, u, f_sp, heat, cool,    &
                                             eta)
      call check_absolute('the_heating_channel_array_is_restored',        &
           maxval(abs(heat_channel_state - hcs0)), 0.0d0, 0.0d0)

      ! ---- (2) the three counter classes ------------------------------ !
      call check_absolute('the_clock_is_restored',                        &
           t_phys, 12.5d0, 0.0d0)
      call check_absolute('the_accepted_step_count_is_restored',          &
           dble(n_steps_accepted), 7.0d0, 0.0d0)
      call check_absolute('the_attempt_statistics_are_not_restored',      &
           dble(n_steps_attempted), 14.0d0, 0.0d0)

      ! ---- (3) the step-size policy ----------------------------------- !
      dt_a = 1.0d0
      call check_absolute('an_admissibility_refusal_halves_the_step',     &
           attempted_step_reduced_dt(dt_a, as_reject_positivity),         &
           0.5d0, 0.0d0)
      call check_absolute('a_carrier_refusal_halves_the_step',            &
           attempted_step_reduced_dt(dt_a, as_reject_carrier),            &
           0.5d0, 0.0d0)
      ! THE DECIDING CLASS SETS THE REDUCTION. The controller factor is
      ! 0.9 e**(-1/(p+1)) with p the order of the class whose lower bound
      ! the gate read, so the same estimate asks for two different
      ! intervals depending on which class refused the step: p = 3 for the
      ! conservative rows and for the electron fraction and the
      ! temperature, p = 1 for the transported fractions. e = 4 gives
      ! 0.9*4**(-1/4) = 0.6364 on the first and 0.9*4**(-1/2) = 0.45 on
      ! the second.
      last_error_estimate = 4.0d0
      last_error_order_p  = err_order_p(err_class_hydro)
      call check_relative('the_deciding_class_sets_the_reduction',        &
           attempted_step_reduced_dt(dt_a, as_reject_int_error),          &
           0.9d0*4.0d0**(-0.25d0), 1.0d-12)
      last_error_order_p  = err_order_p(err_class_species)
      call check_relative('a_species_refusal_reduces_by_its_own_order',   &
           attempted_step_reduced_dt(dt_a, as_reject_int_error),          &
           0.45d0, 1.0d-12)
      last_error_order_p  = err_order_p(err_class_charge_T)
      call check_relative('the_electron_class_reduces_as_the_rk_update',  &
           attempted_step_reduced_dt(dt_a, as_reject_int_error),          &
           0.9d0*4.0d0**(-0.25d0), 1.0d-12)
      last_error_order_p  = err_order_p(err_class_hydro)
      last_error_estimate = 1.0d6
      call check_absolute('the_controller_factor_is_floored_at_0p2',      &
           attempted_step_reduced_dt(dt_a, as_reject_int_error),          &
           0.2d0, 1.0d-14)
      ! A NONFINITE ESTIMATE STEERS NOTHING. It carries no statement about
      ! how much shorter the interval should be, and a negative power of it
      ! would put the nonfinite value into the step, so the policy takes
      ! the halving instead.
      last_error_estimate = ieee_value(1.0d0, ieee_quiet_nan)
      call check_absolute('a_nonfinite_estimate_halves_the_step',         &
           attempted_step_reduced_dt(dt_a, as_reject_int_error),          &
           0.5d0, 0.0d0)
      last_error_estimate = -1.0d0

      call check_absolute('a_reduced_step_bounds_the_next_one',           &
           attempted_step_bound_next_dt(10.0d0, 0.25d0), 0.5d0, 0.0d0)
      call check_absolute('an_unbounded_step_is_the_cfl_step',            &
           attempted_step_bound_next_dt(0.3d0, 0.25d0), 0.3d0, 0.0d0)
      call check_absolute('the_first_step_is_not_bounded',                &
           attempted_step_bound_next_dt(7.0d0, -1.0d0), 7.0d0, 0.0d0)
      call check_absolute('the_outer_retry_cap_matches_the_carrier_one',  &
           dble(n_step_retry_max), 8.0d0, 0.0d0)

      ! ---- (4) the integration-error estimate ------------------------- !
      ! Two states differing by a known amount: the measure is that amount
      ! over atol + rtol|u|, so it is a statement about the arithmetic and
      ! not about a case.
      allocate(u_a(3,1-Ng:N+Ng), u_b(3,1-Ng:N+Ng))
      allocate(f_a(1-Ng:N+Ng,n_species), f_b(1-Ng:N+Ng,n_species))
      allocate(T_a(1-Ng:N+Ng), T_b(1-Ng:N+Ng))
      allocate(ne_a(1-Ng:N+Ng), ne_b(1-Ng:N+Ng))
      allocate(inner_c(1-Ng:N+Ng), inner_T(1-Ng:N+Ng))
      u_a = 1.0d0;  u_b = 1.0d0
      f_a = 0.0d0;  f_b = 0.0d0
      f_a(:,isp_HI)  = 0.8d0;  f_b(:,isp_HI)  = 0.8d0
      f_a(:,isp_HII) = 0.2d0;  f_b(:,isp_HII) = 0.2d0
      T_a = 1.0d0;  T_b = 1.0d0
      ne_a = 0.2d0;  ne_b = 0.2d0
      err_rtol = 1.0d-3;  err_atol = 0.0d0
      call attempted_step_inner_error_reset
      u_a(2,5) = 1.0d0 + 1.0d-3
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      e1 = est%e
      jw = est%jworst;  kw = est%kworst
      ! The difference is one row scale, and the estimate is that times the
      ! factor the RETAINED full step carries.
      call check_relative('the_error_estimate_is_the_scaled_difference',  &
           e1, err_full_step_factor(err_class_hydro), 1.0d-12)
      call check_absolute('the_worst_row_is_named_by_its_class',          &
           dble(est%worst_class), dble(err_class_hydro), 0.0d0)
      ! THE CONVENTION ITSELF: the retained trajectory is the full step, so
      ! the difference of a class is divided by (1 - 2^-p) with the order p
      ! of the integrator that carries THAT class.
      call check_relative('the_retained_step_factor_follows_the_order',   &
           err_full_step_factor(err_class_hydro),                         &
           1.0d0/(1.0d0 - 2.0d0**(-err_order_p(err_class_hydro))),        &
           1.0d-14)
      ! ONE ORDER PER ROW CLASS, AND IT IS ITS INTEGRATOR'S. The three
      ! conservative rows ARE the three-stage Runge-Kutta and the electron
      ! fraction and the temperature are algebraic functions of that
      ! update, so both classes carry p = 3 and a retained-step factor of
      ! 8/7; the transported species and element fractions are composed
      ! with the source at first order, so both carry p = 1 and a factor
      ! of 2. The ladders that measured the two orders are at the
      ! declaration and the suite re-measures one of each below.
      ! The measure is the total departure of the four orders from the
      ! four the ladders measured, so one wrong class fails the row.
      call check_absolute('order_of_each_class_is_its_integrators',       &
           abs(dble(err_order_p(err_class_hydro))    - 3.0d0)             &
           + abs(dble(err_order_p(err_class_charge_T)) - 3.0d0)           &
           + abs(dble(err_order_p(err_class_species))  - 1.0d0)           &
           + abs(dble(err_order_p(err_class_element))  - 1.0d0),          &
           0.0d0, 0.0d0)
      call check_relative('the_conservative_class_factor_is_eight_sevenths',&
           err_full_step_factor(err_class_hydro), 8.0d0/7.0d0, 1.0d-14)
      call check_relative('the_electron_class_factor_is_eight_sevenths',  &
           err_full_step_factor(err_class_charge_T), 8.0d0/7.0d0, 1.0d-14)
      call check_relative('the_species_class_factor_is_two',              &
           err_full_step_factor(err_class_species), 2.0d0, 1.0d-14)
      call check_relative('the_element_class_factor_is_two',              &
           err_full_step_factor(err_class_element), 2.0d0, 1.0d-14)
      call check_absolute('the_estimate_names_the_cell',                  &
           dble(jw), 5.0d0, 0.0d0)
      call check_absolute('the_estimate_names_the_row',                   &
           dble(kw), 2.0d0, 0.0d0)
      ! ORDER: on a difference that halves, the estimate halves. This is
      ! the property of the estimator, measured rather than asserted at a
      ! size, and it is what makes the duty-cycle number readable as an
      ! integration error at all.
      u_a(2,5) = 1.0d0 + 0.5d-3
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      e2 = est%e
      call check_relative('halving_the_difference_halves_the_estimate',   &
           e1/e2, 2.0d0, 1.0d-12)

      ! ---- (4b) THE SPECIES ROWS ARE IN THE ESTIMATE ------------------ !
      ! A difference in a transported species alone, with the conservative
      ! rows, the temperature and the electron fraction identical. Before
      ! item N16b the estimate saw only the three conservative rows and
      ! returned exactly zero for this pair (MEASURED on the entry text;
      ! the report of PLAN_20260909_rev1 item N16b records the number).
      u_a = 1.0d0
      f_a(:,isp_HII) = 0.2d0*(1.0d0 + 1.0d-2)
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      call check_absolute('species_rows_enter_the_estimate',              &
           logical_as_double(est%e .gt. 0.0d0), 1.0d0, 0.0d0)
      call check_absolute('the_species_class_holds_the_worst_row',        &
           dble(est%worst_class), dble(err_class_species), 0.0d0)
      call check_absolute('the_species_row_names_its_f_sp_column',        &
           dble(est%kworst), dble(isp_HII), 0.0d0)
      ! A 1 percent change of a species judged relative to its own value.
      ! The reference is the scale written out: the relative term of the
      ! half-step value plus the floor, which here is the whole hydrogen of
      ! the cell because H I and H II hold all of it.
      call check_relative('the_species_row_is_relative_to_its_own_value', &
           est%e, err_full_step_factor(err_class_species)*0.2d0*1.0d-2    &
                  /(err_rtol*0.2d0 + err_species_floor_of_budget), 1.0d-12)

      ! ---- (4c) THE SPECIES SCALE IS ABUNDANCE-AWARE ------------------ !
      ! A trace carrier at 1e-20 that holds half of the hydrogen its cell
      ! has: a 1 percent change of it is judged by its own value, because
      ! the absolute floor is a fraction of what the species can be and not
      ! a fixed number in code units. With the flat 1e-12 of the
      ! conservative rows the same change would measure 2e-10 and be
      ! invisible.
      f_a = 0.0d0;  f_b = 0.0d0
      f_a(:,isp_HI)  = 1.0d-20;  f_b(:,isp_HI)  = 1.0d-20
      f_b(:,isp_HII) = 1.0d-20
      f_a(:,isp_HII) = 1.0d-20*(1.0d0 + 1.0d-2)
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      call check_relative('species_scale_is_abundance_aware',             &
           est%e, err_full_step_factor(err_class_species)*1.0d-20*1.0d-2  &
                  /(err_rtol*1.0d-20                                      &
                    + err_species_floor_of_budget*2.0d-20), 1.0d-12)
      ! THE SAME PAIR ON THE FLAT ABSOLUTE FLOOR THE CONSERVATIVE ROWS
      ! CARRY, as the reference this row exists to reject: 1e-12 in code
      ! units puts the same 1 percent change at 2e-10, four decades of
      ! tolerance below anything the gate could see.
      call check_relative('the_flat_floor_would_hide_the_same_change',    &
           err_full_step_factor(err_class_species)*1.0d-20*1.0d-2         &
           /(1.0d-12 + err_rtol*1.0d-20), 2.0d-10, 1.0d-6)
      ! The budget the floor is a fraction of: two hydrogen nuclei per unit
      ! of rho*n0 divided by one nucleus per proton, at nd = u(1)*n0.
      call check_relative('the_species_budget_is_its_element_budget',     &
           species_budget_fraction(isp_HII, 1.0d0*n0, 2.0d-20*n0,         &
                                   0.0d0, 0.0d0), 2.0d-20, 1.0d-12)

      ! ---- (4d) AN UNRESOLVED ESTIMATE IS NAMED, NOT JUDGED ----------- !
      ! The inner nonlinear error of the passes above the difference the
      ! estimate measured: the difference is the solvers' own error, so no
      ! class that carries it may decide the step. The estimate still
      ! reports its number; what it does not do is offer it to the gate.
      inner_c = 1.0d0;  inner_T = 0.0d0
      call attempted_step_note_inner_error(inner_c, inner_T, 0.0d0)
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      call check_absolute('unresolved_estimate_is_named_not_judged',      &
           logical_as_double(est%resolved), 0.0d0, 0.0d0)
      call check_absolute('the_unresolved_class_is_named',                &
           dble(est%unresolved_class), dble(err_class_species), 0.0d0)
      call check_absolute('an_unresolved_class_offers_nothing_to_the_gate',&
           est%e_resolved, 0.0d0, 0.0d0)
      call check_absolute('the_unresolved_estimate_still_reports_itself', &
           logical_as_double(est%e .gt. 0.0d0), 1.0d0, 0.0d0)
      call check_relative('the_resolution_test_is_a_tenth',               &
           err_inner_fraction, 0.1d0, 1.0d-14)
      ! ... and with the inner error back below a tenth of the difference
      ! the same pair is resolved and the gate reads it.
      call attempted_step_inner_error_reset
      inner_c = 1.0d-24;  inner_T = 0.0d0
      call attempted_step_note_inner_error(inner_c, inner_T, 0.0d0)
      call attempted_step_error_estimate(u_a, f_a, T_a, ne_a,             &
                                         u_b, f_b, T_b, ne_b, est)
      call check_absolute('a_small_inner_error_leaves_the_estimate_read', &
           logical_as_double(est%resolved), 1.0d0, 0.0d0)
      ! The gate reads the LOWER end of the interval the three inner errors
      ! leave, so only an error their tolerances cannot explain refuses a
      ! step. Here that end is still far above one, and the interval it
      ! sits in is the one the three inner errors of 1e-24 define.
      call check_relative('a_resolved_estimate_reaches_the_gate',         &
           est%e_resolved, est%e_lower, 1.0d-14)
      call check_absolute('the_gate_reads_inside_the_interval',           &
           logical_as_double(est%e_lower .le. est%e .and.                 &
                             est%e .le. est%e_upper), 1.0d0, 0.0d0)
      call check_relative('the_interval_is_the_three_inner_errors',       &
           est%e_upper - est%e_lower,                                     &
           err_full_step_factor(err_class_species)*2.0d0*3.0d0*1.0d-24    &
           /(err_rtol*1.0d-20                                             &
             + err_species_floor_of_budget*2.0d-20), 1.0d-10)
      call attempted_step_inner_error_reset

      ! ---- (5) the formation-energy table ----------------------------- !
      ! Against the constants each entry is built from, so that a change to
      ! one of those constants cannot silently change the reservoir.
      call check_absolute('the_neutral_reference_carries_no_energy',      &
           species_formation_energy(isp_HI), 0.0d0, 0.0d0)
      call check_absolute('the_helium_neutral_carries_no_energy',         &
           species_formation_energy(isp_HeI), 0.0d0, 0.0d0)
      call check_relative('the_proton_carries_the_HI_threshold',          &
           species_formation_energy(isp_HII), e_th_HI, 1.0d-14)
      call check_relative('HeII_carries_the_HeI_threshold',               &
           species_formation_energy(isp_HeII), e_th_HeI, 1.0d-14)
      call check_relative('HeIII_carries_both_thresholds',                &
           species_formation_energy(isp_HeIII), e_th_HeI + e_th_HeII, 1.0d-14)
      call check_relative('the_metastable_carries_its_excitation',        &
           species_formation_energy(isp_HeTR), e_th_HeI - e_th_HeTR, 1.0d-14)
      call check_absolute('H2_carries_a_NEGATIVE_formation_energy',       &
           logical_as_double(species_formation_energy(isp_H2) .lt. 0.0d0),     &
           1.0d0, 0.0d0)
      call check_relative('H2_carries_minus_the_bond_energy',             &
           species_formation_energy(isp_H2), -4.478075d0, 1.0d-4)
      call check_relative('H2p_is_the_ionization_minus_the_bond',         &
           species_formation_energy(isp_H2p), e_th_H2                          &
           + species_formation_energy(isp_H2), 1.0d-12)
      ! mion index 1 is neutral carbon, 2 is C II and 3 is C III.
      call check_absolute('neutral_carbon_carries_no_energy',             &
           species_formation_energy(mion_fsp(1)), 0.0d0, 0.0d0)
      call check_relative('CII_carries_the_first_carbon_threshold',       &
           species_formation_energy(mion_fsp(2)), 11.26d0, 1.0d-12)
      call check_relative('CIII_carries_both_carbon_thresholds',          &
           species_formation_energy(mion_fsp(3)), 11.26d0 + 24.38d0, 1.0d-12)

      ! ---- (6) the thermal identity of the source step ---------------- !
      ! One synthetic step, built so that every term of the identity is
      ! known in advance: a thermal source dt (heat - cool), a composition
      ! change that moves the formation reservoir, and a transport change
      ! that moves the thermal energy with no source behind it.
      run_mode = run_mode_phys
      call build_identity_step(1.0d-6, .false.)
      call check_absolute('the_transport_contribution_is_measured',        &
           logical_as_double(verd%transport_measured), 1.0d0, 0.0d0)
      call check_absolute('the_source_step_closes_its_thermal_energy',     &
           abs(verd%identity_residual)/verd%identity_scale, 0.0d0, 1.0d-14)
      call check_relative('the_transport_contribution_is_what_went_in',    &
           verd%d_u_th_transport, du_transport, 1.0d-12)
      call check_relative('the_source_term_is_dt_times_heat_minus_cool',   &
           verd%q_ext_dt, q_expected, 1.0d-12)
      call check_absolute('the_reservoir_change_is_reported_and_nonzero',  &
           logical_as_double(verd%d_u_form .ne. 0.0d0), 1.0d0, 0.0d0)
      call check_absolute('a_closing_source_step_is_accepted',             &
           logical_as_double(verd%accepted), 1.0d0, 0.0d0)
      ! THE RED REFERENCE. The quantity the instrument used to report is
      ! the whole step's d_u_th plus the reservoir change against the same
      ! source. On this step it is nowhere near zero, so the restated
      ! instrument is not a relabelling of something that already closed.
      call check_absolute('the_old_statement_does_NOT_close',              &
           logical_as_double(abs(verd%d_u_th + verd%d_u_form               &
                                 - verd%q_ext_dt)                          &
                             .gt. 1.0d-3*verd%identity_scale),             &
           1.0d0, 0.0d0)

      ! The same step with a thermal energy the sources did not produce.
      call build_identity_step(1.0d-6, .true.)
      call check_absolute('an_injected_inconsistency_is_refused',          &
           logical_as_double(.not. verd%accepted), 1.0d0, 0.0d0)
      call check_absolute('the_refusal_names_the_source_energy_row',       &
           dble(verd%reason), dble(as_reject_source_energy), 0.0d0)
      call check_absolute('the_refusal_names_the_energy_operation',        &
           dble(verd%operation), dble(as_op_energy), 0.0d0)
      inject_source_energy_offset = 0.0d0
      run_mode = run_mode_init

      if (assertion_failures .gt. 0) then
         write(*,'(a)') 'attempted_step_tests: FAILED'
         stop 1
      endif
      write(*,'(a)') 'attempted_step_tests: every assertion passed'

      contains

      !--------------!

      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double

      !--------------!

      subroutine build_identity_step(dt_step, inject)
      ! One synthetic attempted step, in the order the marching loop takes
      ! it: the checkpoint, the mark before the coupled source step, the
      ! source step itself (a thermal source and a composition change), the
      ! mark after it, and finally a transport change standing for the
      ! boundary condition, the conduction stage and the filter, all of
      ! which write the energy row after the sources.
      real*8,  intent(in) :: dt_step
      logical, intent(in) :: inject
      integer :: j

      inject_source_energy_offset = 0.0d0
      if (inject) inject_source_energy_offset = 1.0d-6

      if (allocated(u))    deallocate(u)
      if (allocated(f_sp)) deallocate(f_sp)
      if (allocated(heat)) deallocate(heat)
      if (allocated(cool)) deallocate(cool)
      if (allocated(eta))  deallocate(eta)
      call build_state()
      allocate(W(3,1-Ng:N+Ng), Tk(1-Ng:N+Ng), dt_cell(1-Ng:N+Ng))
      dt_cell = dt_step
      Tk      = T0
      call attempted_step_checkpoint_take(chk, u, f_sp, heat, cool, eta,   &
                                          0, 0, 0)

      call attempted_step_thermal_energy_before_sources(u)

      ! The source step: the thermal source the energy row is given, and a
      ! composition change beside it. The reservoir moves; the row does not
      ! see it, because heat and cool are already the net thermal source.
      q_expected = 0.0d0
      do j = 1, N
         u(3,j) = u(3,j) + dt_step*(heat(j) - cool(j))
         q_expected = q_expected                                           &
              + dt_step*(heat(j) - cool(j))*r(j)*r(j)*dr_j(j)
      enddo
      f_sp(:,isp_HI)  = f_sp(:,isp_HI)  - 0.05d0
      f_sp(:,isp_HII) = f_sp(:,isp_HII) + 0.05d0

      call attempted_step_thermal_energy_after_sources(u, dt_cell, heat,   &
                                                       cool)

      ! Transport: an energy change with no source behind it, of a size the
      ! test knows.
      du_transport = 0.0d0
      do j = 1, N
         u(3,j) = u(3,j) + 1.0d-4*dble(j)
         du_transport = du_transport + 1.0d-4*dble(j)*r(j)*r(j)*dr_j(j)
      enddo

      do j = 1-Ng, N+Ng
         W(1,j) = u(1,j)
         W(2,j) = u(2,j)/u(1,j)
         W(3,j) = u(3,j)
      enddo

      verd%accepted = .true.
      verd%reason   = as_ok
      call attempted_step_energy_identity(chk, u, W, Tk, f_sp, heat, cool, &
                                          dt_step, verd)

      deallocate(u, f_sp, heat, cool, eta, W, Tk, dt_cell)
      end subroutine build_identity_step

      !--------------!

      subroutine setup_globals()
      ! A twelve-cell hydrogen and helium column with metals, no molecules
      ! and no carrier transport, so that the checkpoint exercised here is
      ! the one an atomic run takes. The molecular half of the list, which
      ! reaches A3's carrier checkpoint, is exercised on a whole binary by
      ! run.sh beside this driver.
      integer :: j
      real*8, allocatable :: recycled(:)
      N   = ncell
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      b0  = 10.0d0
      spherical_domain   = .true.
      thereis_He         = .true.
      thereis_HeITR      = .true.
      thereis_metals     = .true.
      thereis_mol        = .false.
      thereis_oxychem    = .false.
      he_diffusion       = .false.
      carrier_transport  = .false.
      use_excited_H      = .true.
      run_mode           = run_mode_init
      j_min              = 1
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.02d0*dble(j)
         r_edg(j) = 1.0d0 + 0.02d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.02d0
      enddo
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      ! A filled block returned to the allocator, so that the arrays below
      ! are served from a recycled block instead of from a fresh zero page
      ! of the OS. Block (0) measures what allocation leaves in them, and on
      ! an untouched heap that measurement is of the OS and not of the code.
      allocate(recycled(2000))
      recycled = 8.7d244
      deallocate(recycled)
      call ioniz_eq_allocate_arrays()
      call allocate_excited_H_globals()
      bg_ready = .true.
      end subroutine setup_globals

      !--------------!

      subroutine allocate_excited_H_globals()
      ! The eight grid arrays of parameters.f90 that operation 6 writes.
      allocate(gph_balmer_HI(1-Ng:N+Ng), heat_balmer(1-Ng:N+Ng))
      allocate(Jlya_arr(1-Ng:N+Ng), n2s_arr(1-Ng:N+Ng), n2p_arr(1-Ng:N+Ng))
      allocate(Sproton_arr(1-Ng:N+Ng))
      allocate(Hpe_arr(1-Ng:N+Ng), Hdx_arr(1-Ng:N+Ng))
      gph_balmer_HI = 1.0d-3;  heat_balmer = 2.0d-3
      Jlya_arr      = 3.0d-3;  n2s_arr     = 4.0d-3
      n2p_arr       = 5.0d-3;  Sproton_arr = 6.0d-3
      Hpe_arr       = 7.0d-3;  Hdx_arr     = 8.0d-3
      end subroutine allocate_excited_H_globals

      !--------------!

      subroutine build_state()
      integer :: j
      allocate(u(3,1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(heat(1-Ng:N+Ng), cool(1-Ng:N+Ng), eta(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         u(1,j) = exp(-4.0d0*(r(j) - 1.0d0))
         u(2,j) = 0.02d0*u(1,j)
         u(3,j) = 1.5d0*u(1,j)
         heat(j) = 1.0d-5*dble(j + 3)
         cool(j) = 2.0d-5*dble(j + 3)
         eta(j)  = 0.3d0
      enddo
      f_sp = 0.0d0
      f_sp(:,isp_HI)   = 0.8d0
      f_sp(:,isp_HII)  = 0.2d0
      f_sp(:,isp_HeI)  = HeH
      f_sp(:,isp_HeII) = 1.0d-8
      ! The sweep's non-root streak, one of the class-1 items of the list.
      if (allocated(ieq_nonroot_streak)) ieq_nonroot_streak = 0
      ieq_marching_ledger(1)%n_sweep = 3
      ieq_marching_ledger(2)%n_sweep = 5
      end subroutine build_state

      end program attempted_step_tests

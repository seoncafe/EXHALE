      module attempted_step
      ! THE ATTEMPTED-STEP CONTROLLER: one checkpoint, one trial, one
      ! adoption boundary.
      !
      ! The marching loop of EXHALE_main takes one step as fourteen
      ! operations in a fixed order. Rows 1 through 12
      ! change the physical state; the adoption boundary sits immediately
      ! before update_map_end_step, which is the first point at which the
      ! state of the step is complete. This module owns
      !
      !   (a) the checkpoint of everything rows 1 to 12 write, taken before
      !       row 1 and restored on a rejected attempt;
      !   (b) the acceptance predicate evaluated at the boundary;
      !   (c) the step-size policy of a rejection and the retry cap;
      !   (d) the counters of the outer loop, kept apart from the counters
      !       of the nested hydrodynamic retry.
      !
      ! WHAT THIS MODULE DOES NOT CLAIM. A demonstration that the controller
      ! restores state is a demonstration about memory, not about physics.
      ! The energy instrument below closes the THERMAL budget of the source
      ! step; the Shapiro filter of row 12 still alters the adopted state
      ! with no source behind it, so a step the controller adopts while the
      ! filter is active is a step that was RESTORABLE AND ADMISSIBLE, not
      ! a step whose sources balance. That is the remaining part of the
      ! interim contract of the design document, and it
      ! ends when the energy the filter removes is budgeted. No tolerance
      ! anywhere is loosened to make it pass.
      !
      ! THE THREE COUNTER CLASSES.
      !   1 physical accumulations   restored on a rejected attempt, so a
      !                              rejected trial leaves no contribution:
      !                              the CO ceiling record, t_phys,
      !                              n_steps_accepted, the physical ledger
      !                              family.
      !   2 attempt statistics       never restored: they count attempts.
      !   3 diagnostic extrema       kept, not restored, and never reported
      !                              as a property of the adopted state.
      ! Classes 2 and 3 are SAVED here anyway, so that the tests can assert
      ! that they did not move back; the restore skips them.

      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use ion_cell_state,   only: ion_rates
      use species_table,    only: n_mion, mion_fsp, mion_elem,              &
                                  mion_stage, mion_ethr, melem_i0,          &
                                  isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                                  isp_HeIII, isp_HeTR, isp_H2, isp_H2p,     &
                                  isp_H3p, isp_HeHp, isp_OH, isp_H2O,       &
                                  isp_CO
      use ionization_equilibrium
      ! The channel decomposition of the heating, named explicitly because
      ! the checkpoint saves it: it is defined by the assembly that builds
      ! the heat column, not by the sweep that calls it.
      use utils_ion_eq,     only: heat_channel_state
      use excited_hydrogen, only: jlya_rt_loaded, jlya_rt_grid, Tdiag,      &
                                  nhidiag, nediag, nhiidiag, taulya,        &
                                  excited_H_field_ready
      use lya_rt,           only: jint_arr, jstar_arr
      use diffusive_photochemistry, only: carrier_checkpoint,               &
                                  carrier_checkpoint_take,                  &
                                  carrier_checkpoint_restore,               &
                                  carrier_checkpoint_matches
      use energy_semi_implicit, only: n_energy_floor_hits,                  &
                                  n_energy_floor_hits_family,               &
                                  energy_floor_first_step,                  &
                                  energy_floor_last_step,                   &
                                  energy_floor_cell_hits,                   &
                                  energy_res_tol, energy_scale_floor
      use viscous_conduction, only: n_conduction_floor_hits,                &
                                  n_conduction_floor_hits_family,           &
                                  conduction_floor_first_step,              &
                                  conduction_floor_last_step,               &
                                  n_conduction_floor_cells
      use binary_element_diffusion, only: he_fraction_over_one,             &
                                  he_fraction_under_zero,                   &
                                  he_fraction_newton_steps,                 &
                                  he_fraction_newton_resid,                 &
                                  trace_ratio_under_zero,                   &
                                  element_mass_fractions
      ! ONE STOICHIOMETRIC MAP FOR THE ROW SCALES. The absolute floor of a
      ! species row is a fraction of the largest density that species can
      ! reach in the cell, which is the element budget divided by the
      ! nuclei one particle of it holds. Both come from element_inventory,
      ! so the scale of a row and the constraint the same species is
      ! judged against cannot state different stoichiometry.
      use element_inventory, only: nuclei_per_particle,                    &
                                  carrier_element_totals,                  &
                                  ien_H, ien_O, ien_C
      ! The fraction of that largest density at which two states stop being
      ! distinguishable, from the module that owns the element bookkeeping
      ! it is the round-off of.
      use element_census,    only: element_ratio_gate
      use RK_integration,   only: n_faces_flux_positivity_limited,          &
                                  n_faces_flux_positivity_limited_accepted
      ! ONE FORMATION-ENERGY TABLE FOR THE WHOLE CODE. eps_s is written
      ! down in exactly one place, molecular_reaction_heat, and every
      ! reaction heat, every photoevent partition and every reservoir sum
      ! reads it there, so none of them can disagree with the others about
      ! how much energy a particle holds. This module used to carry a
      ! second copy of that table with no entry for OH, H2O and CO, so an
      ! oxygen-network run reported a reservoir with three species missing.
      use molecular_reaction_heat, only: species_formation_energy
      ! ONE EVALUATOR, ONE SET OF TOLERANCES. The
      ! equations a physical step is judged by live in certification.f90 as
      ! its physical-step context; this module owns the checkpoint, the
      ! step-size policy, the energy identity and the integration-error
      ! estimate, and never re-implements a row measure or a tolerance.
      use certification,    only: cert_step_verdict,                       &
                                  certification_evaluate_physical_step,    &
                                  cert_step_reason_text,                   &
                                  cert_step_ok, cert_step_hydro_row

      implicit none
      private

      ! ---------------------------------------------------------------- !
      ! THE FOURTEEN OPERATIONS, in a fixed numbering. A
      ! rejection names the operation it came from with these.
      integer, parameter, public :: as_op_dt          =  0
      integer, parameter, public :: as_op_checkpoint  =  1
      integer, parameter, public :: as_op_hydro       =  2
      integer, parameter, public :: as_op_primitives  =  3
      integer, parameter, public :: as_op_diffusion   =  4
      integer, parameter, public :: as_op_carriers    =  5
      integer, parameter, public :: as_op_excited_H   =  6
      integer, parameter, public :: as_op_ioniz_eq    =  7
      ! Row 8 held the composition projection until it was removed; the
      ! number is kept because the numbering indexes by it, and it now
      ! names the point at which the composition and the energy of the
      ! coupled source step have become one state.
      integer, parameter, public :: as_op_composition =  8
      integer, parameter, public :: as_op_energy      =  9
      integer, parameter, public :: as_op_bc          = 10
      integer, parameter, public :: as_op_conduction  = 11
      integer, parameter, public :: as_op_shapiro     = 12
      integer, parameter, public :: as_op_boundary    = 13
      integer, parameter, public :: as_op_last        = 13

      ! ---------------------------------------------------------------- !
      ! THE RETRY CAP AND THE STEP-SIZE POLICY.
      !
      ! n_step_retry_max = 8 gives the outer controller a floor of dt/2^8,
      ! the same depth at which a carrier interval gives up
      ! (carrier_retry_max), so the failure of one operator and the failure
      ! of the step are reached at the same depth. The nested hydrodynamic
      ! retry_step keeps its own cap of 20, deliberately: its criterion is
      ! the admissibility of one Runge-Kutta stage and its cost is one
      ! stage, not a whole step. The two counts are reported separately and
      ! are never added.
      integer, parameter, public :: n_step_retry_max = 8
      ! A dt reduced by a rejection BOUNDS the next step: dt_next =
      ! min(dt_CFL, 2 dt_accepted). A run that rejects nothing never enters
      ! this branch and is unchanged.
      real*8, parameter, public :: as_dt_growth = 2.0d0

      ! THE TIME-DISCRETE HYDRODYNAMIC RESIDUAL of row 2, the one new
      ! measure this item introduces, lives in certification.f90 with the
      ! other equations and their tolerances (cert_tol_hydro_step): one
      ! evaluator, one set of numbers. Its declaration there states what it
      ! bounds, why no defensible value for it exists yet, and the
      ! measurement that says so.

      ! Absolute floor under the row scale, so that a row of a cell in
      ! which nothing happens cannot divide by zero.
      real*8, parameter, public :: as_scale_floor = 1.0d-300

      ! ---------------------------------------------------------------- !
      ! WHY AN ATTEMPT WAS REFUSED.
      ! Codes 0 to 9 ARE certification's cert_step_* codes, so that the
      ! controller and the evaluator name a refusal with one number. 10 to
      ! 13 are the controller's own: the integration error, the injected
      ! refusal of the tests, the exhaustion of the nested hydrodynamic
      ! retry, and the coupled source step that did not reach its fixed
      ! point, none of which is a statement about an equation. The last one
      ! is separate from as_reject_energy, which IS a statement about the
      ! energy equation: an exhausted fixed point means the temperature and
      ! the composition do not belong to each other, while a failed energy
      ! update means no temperature was found at all.
      integer, parameter, public :: as_ok                 = 0
      integer, parameter, public :: as_reject_nonfinite   = 1
      integer, parameter, public :: as_reject_positivity  = 2
      integer, parameter, public :: as_reject_element     = 3
      integer, parameter, public :: as_reject_carrier     = 5
      integer, parameter, public :: as_reject_energy      = 6
      integer, parameter, public :: as_reject_conduction  = 7
      integer, parameter, public :: as_reject_chem_root   = 8
      integer, parameter, public :: as_reject_hydro_row   = 9
      integer, parameter, public :: as_reject_int_error   = 10
      integer, parameter, public :: as_reject_injected    = 11
      integer, parameter, public :: as_reject_hydro_stage = 12
      integer, parameter, public :: as_reject_source_fixed_point = 13
      ! 14 is the thermal closure of the source step: the change of thermal
      ! energy the source operators produced does not equal the thermal
      ! source they were given. It is separate from as_reject_energy (no
      ! temperature was found at all) and from as_reject_source_fixed_point
      ! (a temperature and a composition were found but have not stopped
      ! moving); this one says that the pair that WAS found does not satisfy
      ! the energy row summed over the grid.
      integer, parameter, public :: as_reject_source_energy = 14

      ! ---------------------------------------------------------------- !
      ! THE CHECKPOINT. Membership rule: an item is in it if some operation
      ! of rows 1 to 12 writes it and some later read of it, in this step
      ! or a later one, precedes a write. Quantities row 3 rebuilds from u
      ! and f_sp before anything reads them (W, rho, v, p, T, the species
      ! densities, n_e, n_tot) are NOT saved; the tests hash them so that
      ! this is a checked property and not an assumption.
      !
      ! FIVE MORE QUANTITIES ARE REBUILT AND NOT SAVED, and they are named
      ! here because each of them is module state that an attempt writes and
      ! the next one reads, so the rule alone does not say which side of it
      ! they fall on. All five are FUNCTIONS OF (u, f_sp), which is what
      ! makes rebuilding them the right answer: saving them would store a
      ! second copy of information the checkpoint already holds, and two
      ! copies of one thing can disagree.
      !
      !  * the caloric composition of caloric_eos (caloric_mixture_active,
      !    nk_per_mass, x_h2, molecular_cell), a function of f_sp alone
      !    (both arrays are ratios of number densities), refreshed by
      !    get_species_densities;
      !  * the lower boundary's face states in BC_Apply (base_face_W,
      !    base_face_lower_W), functions of the first interior cell average
      !    and that composition, refreshed by Apply_BC.
      !
      ! rebuild_state_from_checkpoint rebuilds both, in that order and
      ! BEFORE the primitive state, because the pressure map reads the
      ! first and the boundary reads the second. A restore that left either
      ! standing would hand the next attempt a base face and a pressure
      ! belonging to the refused one: MEASURED at 9.1162e-05 relative on
      ! the base cell velocity of backup/regression/hydrostatic_column
      ! before the order was stated, and at zero after
      ! (attempted_step suite row rejected_pass_leaves_no_state_behind).
      type, public :: attempted_step_checkpoint
         logical :: taken = .false.

         ! ---- class 1, restored: the physical state of the main program
         real*8, allocatable :: u(:,:)
         real*8, allocatable :: f_sp(:,:)
         real*8, allocatable :: heat(:), cool(:), eta(:)

         ! ---- class 1, restored: the carriers, through A3's own
         ! enumeration (it asserts its own completeness, so it is reused
         ! whole rather than listed a second time). It carries the CO
         ! ceiling record, which is a physical accumulation.
         logical :: carr_present = .false.
         type(carrier_checkpoint) :: carr

         ! ---- class 1, restored: the excited-hydrogen populations, their
         ! caches, and the Lyman-alpha field (row 6)
         real*8, allocatable :: gph_balmer_HI(:), heat_balmer(:)
         real*8, allocatable :: Jlya_arr(:), n2s_arr(:), n2p_arr(:)
         real*8, allocatable :: Sproton_arr(:), Hpe_arr(:), Hdx_arr(:)
         real*8, allocatable :: jlya_rt_grid(:)
         real*8, allocatable :: Tdiag(:), nhidiag(:), nediag(:)
         real*8, allocatable :: nhiidiag(:), taulya(:)
         real*8, allocatable :: jint_arr(:), jstar_arr(:)
         logical :: jlya_rt_loaded = .false.
         logical :: excited_H_field_ready = .false.

         ! ---- class 1, restored: the ionization sweep (row 7)
         real*8, allocatable :: nmol_eq(:,:)
         real*8, allocatable :: NH2_col_lw(:), f_shield_lw(:), k_lw_diss(:)
         real*8, allocatable :: p_lw_single(:), p_lw_absorbed(:)
         real*8, allocatable :: tr_lines_lw(:), P_H2_eq(:)
         real*8, allocatable :: nox_eq(:,:), n_o1d_eq(:)
         real*8, allocatable :: NH2O_col(:), NOH_col(:)
         real*8, allocatable :: NCO_col(:), k_co_diss(:), theta_co_shield(:)
         real*8, allocatable :: j_h2o_fuv(:,:), j_oh_fuv(:,:), tau_fuv(:,:)
         real*8, allocatable :: heat_fuv(:), heat_chem(:)
         ! The channel decomposition of the heating of the state the sweep
         ! returned. It is the definition the heat column and
         ! Heating_breakdown.txt are both built from, so it is class-1
         ! state: a trial that is refused after the sweep has to leave it
         ! describing the state that is restored, not the one that was
         ! refused. heat_fuv and heat_chem above are copies of two of its
         ! columns and are saved beside it.
         real*8, allocatable :: heat_channel_state(:,:)
         type(ion_rates), allocatable :: bg_cell(:)
         type(ion_rates), allocatable :: bg_cell_adopted(:)
         type(ion_rates), allocatable :: bg_cell_best(:)
         integer, allocatable :: ieq_nonroot_streak(:)
         integer :: ieq_sweep_state_kind = 0
         integer :: ieq_acc_nprint       = 0
         type(ioniz_eq_ledger) :: ieq_marching_ledger(2)
         type(ioniz_eq_ledger) :: ieq_steady_iterate_ledger
         type(ioniz_eq_ledger) :: ieq_steady_candidate_ledger
         ! The rate record the sweep leaves for the closure evaluator
         ! (A2 step 4): written by row 7, read by the certification of a
         ! later state.
         type(ion_rates), allocatable :: ieq_rate_cell(:)
         real*8, allocatable :: ieq_ne_cell(:), ieq_TK_cell(:)
         real*8, allocatable :: ieq_ntot_cell(:)
         real*8, allocatable :: ieq_met_coef(:,:,:)
         logical :: ieq_rates_ready  = .false.
         integer :: ieq_neq_stored   = 0
         integer :: ieq_mbase_stored = 0
         integer :: ieq_iox_stored   = 0

         ! ---- class 1, restored: the clock and the accepted-step count.
         ! They are advanced only at the adoption boundary, so a rejected
         ! attempt cannot have moved them; they are saved and restored so
         ! that this is enforced rather than relied on.
         real*8  :: t_phys = 0.0d0
         integer :: n_steps_accepted = 0
         integer :: ledger_family = 1

         ! ---- class 2, SAVED AND NOT RESTORED: attempt statistics
         integer :: n_steps_attempted = 0
         integer :: n_dt_halve = 0, n_steps_dt_halved = 0, n_dt_halvings = 0
         integer :: n_faces_flux_positivity_limited = 0
         integer :: n_faces_flux_positivity_limited_accepted = 0

         ! ---- class 3, SAVED AND NOT RESTORED: diagnostic extrema
         real*8  :: he_fraction_over_one = 0.0d0
         real*8  :: he_fraction_under_zero = 0.0d0
         real*8  :: trace_ratio_under_zero = 0.0d0
         integer :: he_fraction_newton_steps = 0
         real*8  :: he_fraction_newton_resid = 0.0d0
         integer :: n_energy_floor_hits = 0
         integer :: n_energy_floor_hits_family(2) = 0
         integer :: energy_floor_first_step = 0, energy_floor_last_step = 0
         integer, allocatable :: energy_floor_cell_hits(:)
         integer :: n_conduction_floor_hits = 0
         integer :: n_conduction_floor_hits_family(2) = 0
         integer :: conduction_floor_first_step = 0
         integer :: conduction_floor_last_step  = 0
         ! viscous_conduction exports the run totals and the cell COUNT but
         ! not its cell-by-cell array, so that array is not saved. It is a
         ! class-3 diagnostic extremum, which is not restored in any case;
         ! what is saved here is what the module makes readable.
         integer :: n_conduction_floor_cells = 0
      end type attempted_step_checkpoint

      ! ---------------------------------------------------------------- !
      ! WHAT THE ACCEPTANCE PREDICATE FOUND, at one adoption boundary.
      type, public :: attempted_step_verdict
         logical :: accepted = .true.
         integer :: reason   = as_ok
         integer :: operation = as_op_boundary
         integer :: jworst   = 0
         integer :: kworst   = 0
         real*8  :: measure  = 0.0d0
         real*8  :: tol      = 0.0d0
         character(len=88) :: what = ''
         ! The worst time-discrete hydrodynamic row of the step, reported
         ! whether or not it refused.
         real*8  :: hydro_row_max = 0.0d0
         integer :: hydro_jworst  = 0
         integer :: hydro_kworst  = 0
         logical :: hydro_row_above = .false.
         ! THE THERMAL-ENERGY IDENTITY OF THE STEP (T1.5 as this code's
         ! sources define it; see energy_identity below for the statement
         ! and its derivation from the B3c convention).
         !
         !   d_u_th          the thermal energy the WHOLE attempted step
         !                   changed, summed over the physical cells
         !   d_u_th_transport  how much of that the transport operators did
         !                   (the hydrodynamic stages, the diffusion and
         !                   carrier operators, the boundary condition, the
         !                   conduction stage and the filter), MEASURED by
         !                   differencing the thermal energy across the
         !                   source step and subtracting
         !   q_ext_dt        the thermal source the step was given,
         !                   sum_j dV_j dt (heat - cool)_j
         !   identity_residual  d_u_th - d_u_th_transport - q_ext_dt, the
         !                   closure of the SOURCE STEP alone
         !   identity_scale  sum_j dV_j [ |u_th^n| + dt(|heat|+|cool|) ],
         !                   the row scale of energy_semi_implicit summed
         !   d_u_form        the formation reservoir's change over the step.
         !                   REPORTED as the reservoir's change and NEVER
         !                   added to the row: with this code's net thermal
         !                   heat and cool it is already accounted for,
         !                   and adding it would count every
         !                   collisional transfer twice.
         logical :: identity_evaluated = .false.
         ! .true. when the two thermal-energy marks of the source step were
         ! both taken during this attempt. Without them the transport
         ! contribution is not measured, identity_residual falls back to the
         ! whole-step imbalance and the gate below is inert.
         logical :: transport_measured = .false.
         real*8  :: d_u_th    = 0.0d0
         real*8  :: d_u_th_transport = 0.0d0
         real*8  :: d_u_form  = 0.0d0
         real*8  :: q_ext_dt  = 0.0d0
         real*8  :: identity_residual = 0.0d0
         real*8  :: identity_scale    = 0.0d0
      end type attempted_step_verdict

      ! ---------------------------------------------------------------- !
      ! OUTER-LOOP COUNTERS, kept apart from the nested hydrodynamic
      ! retry's (class 2: attempt statistics, never restored).
      ! EACH LOOP COUNTS ITS OWN ATTEMPTS AND NOTHING ELSE: an outer
      ! attempt that contains three inner halvings contributes 1 here and 3
      ! to n_steps_attempted. The two are printed with their definitions
      ! and are never added.
      integer, save, public :: n_outer_attempts   = 0
      integer, save, public :: n_step_rejections  = 0
      integer, save, public :: n_steps_with_rejection = 0
      integer, save, public :: n_step_exhaustions = 0
      ! Rejections by reason, indexed by the as_reject_* codes.
      integer, save, public :: n_reject_by_reason(0:as_reject_source_energy) = 0
      ! Rejections by the operation that produced the refused state.
      integer, save, public :: n_reject_by_operation(0:as_op_last) = 0

      ! The two B6 category-4 counters of rows 8 and 12 live in
      ! global_parameters (their declaration there says why): the
      ! certification reads them to refuse a state whose history contains
      ! one, and it cannot use this module.

      ! ---------------------------------------------------------------- !
      ! TEST HOOKS. Each forces one named thing to happen; all default off,
      ! so a run without them is the run an unguarded build gives. They are
      ! named for what they do, not for the test that uses them.
      !
      ! Refuse the attempt immediately after this operation of the fourteen
      ! (0 = never). Set from EXHALE_REJECT_AFTER_OP.
      integer, save, public :: reject_after_operation = 0
      ! ... and only at this step count (-1 = at every step).
      integer, save, public :: reject_at_step = -1
      ! ... this many leading attempts of it (1 = one attempt, so the retry
      ! succeeds; a number above n_step_retry_max drives the exhaustion).
      integer, save, public :: reject_leading_attempts = 1
      ! ... and only in this pass of an error-control macrostep, 1 for the
      ! full step and 2 and 3 for the two halves (0 = in whichever pass
      ! comes first). It exists because a refusal in a half pass and a
      ! refusal in the full pass are different statements about the
      ! transaction and each has to be reachable on demand. Set from
      ! EXHALE_AS_REJECT_IN_PASS.
      integer, save, public :: reject_in_pass = 0
      ! WHICH PASS OF THE ERROR-CONTROL MACROSTEP IS RUNNING: 1 the full
      ! step, 2 and 3 the two halves of the same interval, 0 outside a
      ! sampled macrostep. Written by the marching loop of EXHALE_main,
      ! which owns the passes; read here only by the injected refusal.
      integer, save, public :: err_pass_index = 0
      ! How many injected refusals have been served at the current step.
      integer, save :: injected_served = 0
      integer, save :: injected_at_count = -2
      ! Add this relative offset to the thermal energy recorded at the exit
      ! of the source step, so that the closure below is asked a question
      ! whose answer is known: the step then carries a thermal-energy change
      ! that its sources did not produce, and the gate has to refuse it.
      ! Set from EXHALE_INJECT_SOURCE_ENERGY; 0 is the unguarded run.
      real*8, save, public :: inject_source_energy_offset = 0.0d0
      ! Multiply the FIRST integration-error estimate of the run by this,
      ! so that the acceptance gate below is asked a question whose answer
      ! is known: a large value puts the estimate above one with every
      ! local solve succeeding, and the macrostep has to be refused for its
      ! integration error and nothing else. Set from EXHALE_AS_ERR_INJECT;
      ! 0 is the unguarded run.
      real*8, save, public :: err_inject_factor = 0.0d0
      ! ... on this many leading estimates of the run (1 = the first, so
      ! the shortened retry is judged on its own error and the run goes on;
      ! a number above n_step_retry_max drives the macrostep to the named
      ! exhaustion on its integration error). Set from
      ! EXHALE_AS_ERR_INJECT_COUNT.
      integer, save, public :: err_inject_estimates = 1
      ! FORCE THE RESOLUTION TEST TO FAIL on this many leading estimates of
      ! the run, so that the named UNRESOLVED outcome is reachable on
      ! demand: every class that carries a difference is reported
      ! unresolved, whatever its inner error actually was. It exists
      ! because the rule that an estimate whose inner nonlinear error is
      ! not below the difference it measures judges NOTHING cannot be
      ! tested on a run whose solves happen to be tight enough. Set from
      ! EXHALE_AS_FORCE_UNRESOLVED; 0 is the unguarded run.
      integer, save, public :: force_unresolved_estimates = 0
      ! A PRESCRIBED STEP IN SECONDS in place of the CFL step, read by the
      ! marching loop of EXHALE_main. It exists so that one interval can be
      ! covered by one step, by two, by four and by eight from the same
      ! state and the order of the complete split update measured from the
      ! successive differences; 0 is the unguarded run.
      real*8, save, public :: as_fixed_dt_seconds = 0.0d0

      ! ---------------------------------------------------------------- !
      ! THE TWO MARKS OF THE SOURCE STEP. The marching loop takes the
      ! thermal energy of the grid immediately before the coupled source
      ! step and immediately after it; the difference is what the source
      ! operators did, and everything else the step did to the thermal
      ! energy is transport. Both are volume-weighted sums over the
      ! physical cells, in the code's energy density unit times volume.
      ! They are reset by attempted_step_checkpoint_take, which runs once
      ! per attempt, so a mark can never be read across attempts.
      real*8, save :: uth_before_sources = 0.0d0
      real*8, save :: uth_after_sources  = 0.0d0
      ! sum_j dV_j dt_j (heat - cool)_j, formed at the exit mark with the
      ! cell's own dt, because a local-dt run advances each cell by its own
      ! interval. In phys mode the two are the same number: that mode
      ! refuses local time stepping (input_read.f90).
      real*8, save :: q_ext_dt_sources   = 0.0d0
      logical, save :: marked_before_sources = .false.
      logical, save :: marked_after_sources  = .false.

      ! ---------------------------------------------------------------- !
      ! THE STEP-DOUBLING INTEGRATION-ERROR ESTIMATE. Every
      ! n_err_every accepted steps the same step
      ! is retaken from the same checkpoint as two of dt/2, and
      !
      !   e = max over cells and rows of
      !         |u(dt) - u(dt/2 twice)| / (atol + rtol |u|)
      !
      ! is formed. The two half steps run the WHOLE operator split, so the
      ! estimate is of the split and not of one operator: the carrier
      ! operator was MEASURED to be first order in its substep, so an estimator
      ! that ignored the source operators would report the accuracy of the
      ! one stage that is third order.
      !
      ! In phys mode the default is 20, which costs the +2/n_err_every =
      ! +10 percent the design states and is the price of a physical-mode
      ! claim. In init mode it is off unless asked: an initialization run
      ! makes no statement about a trajectory, so it has no integration
      ! error to bound.
      integer, save, public :: n_err_every = 0
      real*8,  save, public :: err_rtol    = 1.0d-3
      real*8,  save, public :: err_atol    = 1.0d-12
      ! The last estimate taken, the step it was taken at, and where.
      real*8,  save, public :: last_error_estimate = -1.0d0
      integer, save, public :: last_error_step     = -1
      integer, save, public :: last_error_cell     = 0
      integer, save, public :: last_error_row      = 0
      integer, save, public :: n_error_estimates   = 0
      ! Seconds spent inside the estimate, so its cost is measured and not
      ! argued.
      real*8,  save, public :: error_estimate_seconds = 0.0d0

      ! ---------------------------------------------------------------- !
      ! THE ROW CLASSES OF THE ESTIMATE. A conservative density, a species
      ! fraction, an element mass fraction and a temperature are different
      ! physical quantities held in different units, so one absolute floor
      ! cannot serve all four. Each class carries its own scale, stated at
      ! the scale, and the estimate reports which class held the worst row.
      !
      !   hydro        the three conservative rows of the mass, momentum
      !                and total-energy equations, in code units
      !   species      the carrier fractions of the transported chemistry,
      !                each a number of particles per unit of rho*n0
      !   element      the transported element mass fractions (helium, and
      !                the trace metals when trace-metal diffusion is on), which
      !                exist only where the element-diffusion operator runs
      !   charge_T     the electron fraction and the temperature, the two
      !                quantities through which a trace species reaches the
      !                opacity, the chemistry and the caloric equation of
      !                state. They are derived from the rows above, and
      !                they are a class of their own because the map that
      !                derives them is nonlinear: a composition difference
      !                too small to show in any one species row can still
      !                move the electron density and the temperature.
      integer, parameter, public :: err_class_hydro    = 1
      integer, parameter, public :: err_class_species  = 2
      integer, parameter, public :: err_class_element  = 3
      integer, parameter, public :: err_class_charge_T = 4
      integer, parameter, public :: n_err_class        = 4
      character(len=34), parameter, public ::                              &
         err_class_name(n_err_class) = [                                   &
            'the conservative hydrodynamic rows',                          &
            'the transported species fractions ',                          &
            'the transported element fractions ',                          &
            'the electron fraction and T       ' ]

      ! THE ORDER OF EACH ROW CLASS IS THE ORDER OF THE INTEGRATOR THAT
      ! CARRIES IT, AND THE ESTIMATOR CONVENTION.
      !
      ! The trajectory the run keeps is the FULL step, the one the
      ! acceptance predicate was evaluated on. For a class whose update has
      ! order p the leading local error of one step of length H is
      ! C H^(p+1), so
      !
      !   Delta = q(H) - q(H/2 twice) = C H^(p+1) (1 - 2^-p),
      !
      ! and the local error of the RETAINED full step is Delta/(1 - 2^-p),
      ! which is err_full_step_factor below. (Had the two-half result been
      ! retained the factor would be 1/(2^p - 1) instead: the two
      ! conventions describe different states.)
      !
      ! p IS NOT ONE NUMBER FOR THE WHOLE UPDATE, because the classes are
      ! not advanced by one integrator:
      !
      !  * the three conservative rows ARE the three-stage Runge-Kutta of
      !    the hydrodynamics, whose local error is H^4, and the electron
      !    fraction and the temperature are algebraic functions of that
      !    update and inherit its order. p = 3 for both classes.
      !  * the transported species fractions and the transported element
      !    mass fractions are advected in a stage composed with the source
      !    at FIRST order, and the splitting error, not the stage, is what
      !    their difference measures. p = 1 for both classes.
      !
      ! MEASURED on backup/regression/hydrostatic_column in physical mode,
      ! one macrostep refused at every interval the policy offers so that
      ! the same entry state is retaken at H, then 0.2 H, then 0.04 H, with
      ! every class measured at each. The ratios
      ! log(Delta_i/Delta_i+1)/log(H_i/H_i+1), which are p + 1:
      !
      !   the conservative rows        3.9966  3.9967  4.0104  3.9978
      !   the electron fraction and T  3.9962  3.9926  3.9867
      !   the transported species      1.9987  1.9799
      !
      ! The suite rows temporal_order_of_* re-measure one ladder of each
      ! order on every run, so a change that breaks either order is visible
      ! without taking the ladder by hand. The element class carries no
      ! ladder of its own: it exists only where the element-diffusion
      ! operator runs, and it is advected by the same first-order
      ! composition, so it is given the species class's order and that is
      ! an assignment by construction, not a measurement.
      !
      ! WHY IT WAS ONE NUMBER BEFORE, AND WHY THAT NUMBER WAS 1. A state
      ! that a refused pass left behind reached the base face of the next
      ! pass, an O(H) error in one face of one stage of a step of H/2,
      ! which is an O(H^2) contribution to every class's difference and was
      ! 40 times the real one. Every class then read as H^2 and p = 1 was
      ! calibrated on it. With that state rebuilt on every restore the
      ! ladders above are the integrators'.
      !
      ! p ALSO FIXES THE STEP-SIZE POLICY: attempted_step_reduced_dt
      ! reduces by 0.9 e^(-1/(p+1)), and the p it uses is the DECIDING
      ! class's, the class whose lower bound is the estimate the gate read.
      ! Sizing the next interval by the exponent of a class that did not
      ! refuse the step would use an order nothing about this rejection
      ! measured.
      integer, parameter, public :: err_order_p(n_err_class) =             &
                                    [ 3, 1, 1, 3 ]
      real*8,  parameter, public :: err_full_step_factor(n_err_class) =    &
                                    1.0d0/(1.0d0 - 2.0d0**(-err_order_p))
      ! The order of the class whose lower bound the last estimate handed
      ! the gate, which is the order the reduction policy reads. It is set
      ! by every estimate; the initial value is the conservative rows',
      ! which are the class that decides every estimate measured on the
      ! cases of this suite.
      integer, save, public :: last_error_order_p =                        &
                               err_order_p(err_class_hydro)

      ! THE ABSOLUTE FLOOR OF A SPECIES ROW, as a fraction of the largest
      ! density that species can reach in that cell (its element budget
      ! divided by the nuclei one of its particles holds). It is
      ! element_ratio_gate, the ratio at which the element bookkeeping of
      ! this code stops being able to tell two states apart: a difference
      ! below that fraction of the element is inside the round-off of the
      ! repartitions that produced it and is not a temporal error. Read
      ! from element_census so the two cannot drift apart.
      !
      ! WHY NOT A FLAT ABSOLUTE FLOOR. err_atol is 1e-12 on a quantity held
      ! in code units, so what it protects depends on the normalization and
      ! not on any physical smallness: with it, a species whose WHOLE
      ! element budget is 1e-12 of the gas is declared unresolvable even
      ! when the step destroys all of it, while an abundant species is
      ! judged by its own value in any case. The budget-relative floor
      ! judges every species by what it can be.
      real*8, parameter, public :: err_species_floor_of_budget =           &
                                   element_ratio_gate

      ! THE INNER NONLINEAR ERROR HAS TO BE BELOW THE TEMPORAL ONE.
      ! The three passes each stop
      ! their nonlinear solves at a finite tolerance, and the difference
      ! the estimate measures carries the errors of all three, so a
      ! difference of the size of those tolerances says nothing about the
      ! discretization. Below this fraction the interval the three inner
      ! errors leave around the difference is narrower than 30 percent of
      ! it, and the point value describes the class; above it the class
      ! reports itself unresolved and the verdict is taken from the
      ! interval instead, which can then be UNRESOLVED: neither a
      ! rejection nor an acceptance, and the step is not judged by it.
      real*8, parameter, public :: err_inner_fraction = 0.1d0
      ! How many inner errors the measured difference carries: the full
      ! pass leaves one and each half pass leaves one.
      integer, parameter, public :: n_err_passes_of_a_macrostep = 3

      ! WHAT ONE ESTIMATE FOUND.
      type, public :: integration_error
         ! The worst SCALED row over the classes that are resolved, which
         ! is the number the acceptance gate reads, and the worst over all
         ! classes, which is what the estimate would say if every class
         ! were trusted.
         real*8  :: e          = 0.0d0
         real*8  :: e_resolved = 0.0d0
         integer :: worst_class = 0
         integer :: jworst = 0
         integer :: kworst = 0
         ! Per class: the scaled measure, where it sat, and the raw
         ! difference and the scale it was divided by, so that the
         ! tolerance a state would have needed can be read off a run
         ! without taking the estimate again.
         real*8  :: e_class(n_err_class) = 0.0d0
         integer :: j_class(n_err_class) = 0
         integer :: k_class(n_err_class) = 0
         real*8  :: d_class(n_err_class) = 0.0d0
         real*8  :: s_class(n_err_class) = 0.0d0
         ! The magnitude of the quantity itself at the worst row, which is
         ! what a RELATIVE inner tolerance has to be multiplied by to
         ! become an error of the same kind as the difference.
         real*8  :: q_class(n_err_class) = 0.0d0
         logical :: class_present(n_err_class)  = .false.
         logical :: class_resolved(n_err_class) = .true.
         ! THE INTERVAL THE MEASURED DIFFERENCE LEAVES THE TRUE TEMPORAL
         ! ERROR IN. The three passes each carry their own inner nonlinear
         ! error, so the difference bounds the temporal error only to
         ! within their sum: e_lower is what the step is at least, and
         ! e_upper what it can be at most.
         real*8  :: e_lower = 0.0d0
         real*8  :: e_upper = 0.0d0
         ! .false. when the interval above straddles the acceptance limit,
         ! so that neither a rejection nor an acceptance follows from it.
         logical :: resolved = .true.
         integer :: unresolved_class = 0
         ! The class whose lower bound the gate reads. Its order is the
         ! one attempted_step_reduced_dt raises the estimate to, so that
         ! the interval a rejection asks for next is sized by the accuracy
         ! of the quantity that refused this one.
         integer :: deciding_class = 0
         ! The largest inner error over the measured difference, over the
         ! classes present, and which class held it.
         real*8  :: inner_over_difference = 0.0d0
      end type integration_error

      ! THE INNER NONLINEAR ERROR OF THE PASSES OF ONE MACROSTEP, collected
      ! by the marching loop as each pass leaves its solves and reset when
      ! the transaction (re)starts. Each is in the norm its own solver's
      ! tolerance is written in: an absolute species fraction for the
      ! coupled pair's composition, a relative temperature for its
      ! temperature, and the relative row imbalance of the carrier system
      ! for the carrier substep.
      ! CELL BY CELL, because the difference the estimate measures is at ONE
      ! cell and one row: the coupled pair states its tolerances as maxima
      ! over the grid, and comparing a grid maximum with a local difference
      ! would call a well resolved row unresolved because some other cell
      ! is still moving. The arrays carry the pair's remaining error at
      ! each cell, and the scalars beside them are their grid maxima, kept
      ! for the report.
      real*8, save, public, allocatable :: inner_error_comp_cell(:)
      real*8, save, public, allocatable :: inner_error_T_cell(:)
      real*8, save, public :: inner_error_composition = 0.0d0
      real*8, save, public :: inner_error_temperature = 0.0d0
      real*8, save, public :: inner_error_carrier_row = 0.0d0
      ! Estimates that could not be trusted, counted so that a run says how
      ! much of its trajectory was error-controlled and how much was not.
      integer, save, public :: n_error_estimates_unresolved = 0

      public :: attempted_step_checkpoint_take
      public :: attempted_step_checkpoint_restore
      public :: attempted_step_checkpoint_matches
      public :: attempted_step_checkpoint_checksum
      public :: attempted_step_perturb_checkpointed_state_for_test
      public :: attempted_step_evaluate
      public :: attempted_step_injected_refusal
      public :: attempted_step_note_step
      public :: attempted_step_reduced_dt
      public :: attempted_step_bound_next_dt
      public :: attempted_step_error_estimate
      public :: attempted_step_inner_error_reset
      public :: attempted_step_note_inner_error
      public :: attempted_step_species_rows
      public :: species_budget_fraction
      public :: attempted_step_report
      public :: attempted_step_read_environment
      public :: attempted_step_reason_text, attempted_step_operation_text
      public :: attempted_step_thermal_energy_before_sources
      public :: attempted_step_thermal_energy_after_sources
      public :: attempted_step_energy_identity

      contains

      ! ================================================================ !
      ! 1. THE CHECKPOINT
      ! ================================================================ !

      subroutine attempted_step_checkpoint_take(chk, u, f_sp, heat, cool,  &
                                                eta, n_halve, n_halvings,  &
                                                n_steps_halved)
      ! Taken before operation 1 of the step.
      type(attempted_step_checkpoint), intent(inout) :: chk
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng),            intent(in) :: heat, cool, eta
      ! The nested positivity retry's three counters, which live in the
      ! marching loop. They are class-2 attempt statistics: saved so the
      ! tests can assert the restore does NOT move them back.
      integer, intent(in) :: n_halve, n_halvings, n_steps_halved

      ! The two thermal-energy marks of the source step belong to ONE
      ! attempt: a trial that is refused before the source step must not
      ! leave a mark the next attempt could read as its own.
      marked_before_sources = .false.
      marked_after_sources  = .false.

      ! ---- the main-program state
      call save_grid_state(chk%u,    u)
      call save_grid_fractions(chk%f_sp, f_sp)
      call save_r1f(chk%heat, heat)
      call save_r1f(chk%cool, cool)
      call save_r1f(chk%eta,  eta)

      ! ---- the carriers. A3's checkpoint is taken whole; it carries the
      ! CO ceiling record and asserts its own completeness. The carrier
      ! fractions themselves live in f_sp at this point in the step and are
      ! saved above, so the optional fc argument is not passed.
      if (thereis_mol) then
         call carrier_checkpoint_take(chk%carr)
         chk%carr_present = .true.
      else
         chk%carr_present = .false.
      endif

      ! ---- the excited-hydrogen populations, caches and the Lya field
      call save_r1(chk%gph_balmer_HI, gph_balmer_HI)
      call save_r1(chk%heat_balmer,   heat_balmer)
      call save_r1(chk%Jlya_arr,      Jlya_arr)
      call save_r1(chk%n2s_arr,       n2s_arr)
      call save_r1(chk%n2p_arr,       n2p_arr)
      call save_r1(chk%Sproton_arr,   Sproton_arr)
      call save_r1(chk%Hpe_arr,       Hpe_arr)
      call save_r1(chk%Hdx_arr,       Hdx_arr)
      call save_r1(chk%jlya_rt_grid,  jlya_rt_grid)
      call save_r1(chk%Tdiag,         Tdiag)
      call save_r1(chk%nhidiag,       nhidiag)
      call save_r1(chk%nediag,        nediag)
      call save_r1(chk%nhiidiag,      nhiidiag)
      call save_r1(chk%taulya,        taulya)
      call save_r1(chk%jint_arr,      jint_arr)
      call save_r1(chk%jstar_arr,     jstar_arr)
      chk%jlya_rt_loaded       = jlya_rt_loaded
      chk%excited_H_field_ready = excited_H_field_ready

      ! ---- the ionization sweep
      call save_r2(chk%nmol_eq,       nmol_eq)
      call save_r1(chk%NH2_col_lw,    NH2_col_lw)
      call save_r1(chk%f_shield_lw,   f_shield_lw)
      call save_r1(chk%k_lw_diss,     k_lw_diss)
      call save_r1(chk%p_lw_single,   p_lw_single)
      call save_r1(chk%p_lw_absorbed, p_lw_absorbed)
      call save_r1(chk%tr_lines_lw,   tr_lines_lw)
      call save_r1(chk%P_H2_eq,       P_H2_eq)
      call save_r2(chk%nox_eq,        nox_eq)
      call save_r1(chk%n_o1d_eq,      n_o1d_eq)
      call save_r1(chk%NH2O_col,      NH2O_col)
      call save_r1(chk%NOH_col,       NOH_col)
      call save_r1(chk%NCO_col,       NCO_col)
      call save_r1(chk%k_co_diss,     k_co_diss)
      call save_r1(chk%theta_co_shield, theta_co_shield)
      call save_r2(chk%j_h2o_fuv,     j_h2o_fuv)
      call save_r2(chk%j_oh_fuv,      j_oh_fuv)
      call save_r2(chk%tau_fuv,       tau_fuv)
      call save_r1(chk%heat_fuv,      heat_fuv)
      call save_r1(chk%heat_chem,     heat_chem)
      call save_r2(chk%heat_channel_state, heat_channel_state)
      call save_ir(chk%bg_cell,         bg_cell)
      call save_ir(chk%bg_cell_adopted, bg_cell_adopted)
      call save_ir(chk%bg_cell_best,    bg_cell_best)
      call save_i1(chk%ieq_nonroot_streak, ieq_nonroot_streak)
      chk%ieq_sweep_state_kind        = ieq_sweep_state_kind
      chk%ieq_acc_nprint              = ieq_acc_nprint
      chk%ieq_marching_ledger         = ieq_marching_ledger
      chk%ieq_steady_iterate_ledger   = ieq_steady_iterate_ledger
      chk%ieq_steady_candidate_ledger = ieq_steady_candidate_ledger
      call save_ir(chk%ieq_rate_cell, ieq_rate_cell)
      call save_r1(chk%ieq_ne_cell,   ieq_ne_cell)
      call save_r1(chk%ieq_TK_cell,   ieq_TK_cell)
      call save_r1(chk%ieq_ntot_cell, ieq_ntot_cell)
      call save_r3(chk%ieq_met_coef,  ieq_met_coef)
      chk%ieq_rates_ready  = ieq_rates_ready
      chk%ieq_neq_stored   = ieq_neq_stored
      chk%ieq_mbase_stored = ieq_mbase_stored
      chk%ieq_iox_stored   = ieq_iox_stored

      ! ---- the clock and the accepted-step count
      chk%t_phys           = t_phys
      chk%n_steps_accepted = n_steps_accepted
      chk%ledger_family    = ledger_family

      ! ---- class 2 and class 3, saved for the assertions only
      chk%n_steps_attempted = n_steps_attempted
      chk%n_dt_halve        = n_halve
      chk%n_steps_dt_halved = n_steps_halved
      chk%n_dt_halvings     = n_halvings
      chk%n_faces_flux_positivity_limited = n_faces_flux_positivity_limited
      chk%n_faces_flux_positivity_limited_accepted =                       &
           n_faces_flux_positivity_limited_accepted
      chk%he_fraction_over_one      = he_fraction_over_one
      chk%he_fraction_under_zero    = he_fraction_under_zero
      chk%trace_ratio_under_zero    = trace_ratio_under_zero
      chk%he_fraction_newton_steps  = he_fraction_newton_steps
      chk%he_fraction_newton_resid  = he_fraction_newton_resid
      chk%n_energy_floor_hits        = n_energy_floor_hits
      chk%n_energy_floor_hits_family = n_energy_floor_hits_family
      chk%energy_floor_first_step    = energy_floor_first_step
      chk%energy_floor_last_step     = energy_floor_last_step
      call save_i1(chk%energy_floor_cell_hits, energy_floor_cell_hits)
      chk%n_conduction_floor_hits        = n_conduction_floor_hits
      chk%n_conduction_floor_hits_family = n_conduction_floor_hits_family
      chk%conduction_floor_first_step    = conduction_floor_first_step
      chk%conduction_floor_last_step     = conduction_floor_last_step
      chk%n_conduction_floor_cells       = n_conduction_floor_cells()

      chk%taken = .true.
      end subroutine attempted_step_checkpoint_take

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_checkpoint_restore(chk, u, f_sp, heat,     &
                                                   cool, eta)
      ! Restores class 1 and NOTHING ELSE. Class 2 (attempt statistics) and
      ! class 3 (diagnostic extrema) are deliberately left where the
      ! rejected attempt put them: they are the only
      ! trace of what the refused direction did, and they can never
      ! invalidate the adopted state because a rejected trial contributes
      ! nothing to the physical history.
      type(attempted_step_checkpoint), intent(in) :: chk
      real*8, dimension(3,1-Ng:N+Ng),          intent(out) :: u
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng),            intent(out) :: heat,cool,eta

      if (.not. chk%taken) then
         write(*,'(a)') ' attempted_step: RESTORE FROM A CHECKPOINT THAT'//&
              ' WAS NEVER TAKEN. This is a defect of the caller.'
         error stop 1
      endif

      u    = chk%u
      f_sp = chk%f_sp
      heat = chk%heat
      cool = chk%cool
      eta  = chk%eta

      if (chk%carr_present) call carrier_checkpoint_restore(chk%carr)

      call put_r1(chk%gph_balmer_HI, gph_balmer_HI)
      call put_r1(chk%heat_balmer,   heat_balmer)
      call put_r1(chk%Jlya_arr,      Jlya_arr)
      call put_r1(chk%n2s_arr,       n2s_arr)
      call put_r1(chk%n2p_arr,       n2p_arr)
      call put_r1(chk%Sproton_arr,   Sproton_arr)
      call put_r1(chk%Hpe_arr,       Hpe_arr)
      call put_r1(chk%Hdx_arr,       Hdx_arr)
      call put_r1(chk%jlya_rt_grid,  jlya_rt_grid)
      call put_r1(chk%Tdiag,         Tdiag)
      call put_r1(chk%nhidiag,       nhidiag)
      call put_r1(chk%nediag,        nediag)
      call put_r1(chk%nhiidiag,      nhiidiag)
      call put_r1(chk%taulya,        taulya)
      call put_r1(chk%jint_arr,      jint_arr)
      call put_r1(chk%jstar_arr,     jstar_arr)
      jlya_rt_loaded        = chk%jlya_rt_loaded
      excited_H_field_ready = chk%excited_H_field_ready

      call put_r2(chk%nmol_eq,       nmol_eq)
      call put_r1(chk%NH2_col_lw,    NH2_col_lw)
      call put_r1(chk%f_shield_lw,   f_shield_lw)
      call put_r1(chk%k_lw_diss,     k_lw_diss)
      call put_r1(chk%p_lw_single,   p_lw_single)
      call put_r1(chk%p_lw_absorbed, p_lw_absorbed)
      call put_r1(chk%tr_lines_lw,   tr_lines_lw)
      call put_r1(chk%P_H2_eq,       P_H2_eq)
      call put_r2(chk%nox_eq,        nox_eq)
      call put_r1(chk%n_o1d_eq,      n_o1d_eq)
      call put_r1(chk%NH2O_col,      NH2O_col)
      call put_r1(chk%NOH_col,       NOH_col)
      call put_r1(chk%NCO_col,       NCO_col)
      call put_r1(chk%k_co_diss,     k_co_diss)
      call put_r1(chk%theta_co_shield, theta_co_shield)
      call put_r2(chk%j_h2o_fuv,     j_h2o_fuv)
      call put_r2(chk%j_oh_fuv,      j_oh_fuv)
      call put_r2(chk%tau_fuv,       tau_fuv)
      call put_r1(chk%heat_fuv,      heat_fuv)
      call put_r1(chk%heat_chem,     heat_chem)
      call put_r2(chk%heat_channel_state, heat_channel_state)
      call put_ir(chk%bg_cell,         bg_cell)
      call put_ir(chk%bg_cell_adopted, bg_cell_adopted)
      call put_ir(chk%bg_cell_best,    bg_cell_best)
      call put_i1(chk%ieq_nonroot_streak, ieq_nonroot_streak)
      ieq_sweep_state_kind        = chk%ieq_sweep_state_kind
      ieq_acc_nprint              = chk%ieq_acc_nprint
      ieq_marching_ledger         = chk%ieq_marching_ledger
      ieq_steady_iterate_ledger   = chk%ieq_steady_iterate_ledger
      ieq_steady_candidate_ledger = chk%ieq_steady_candidate_ledger
      call put_ir(chk%ieq_rate_cell, ieq_rate_cell)
      call put_r1(chk%ieq_ne_cell,   ieq_ne_cell)
      call put_r1(chk%ieq_TK_cell,   ieq_TK_cell)
      call put_r1(chk%ieq_ntot_cell, ieq_ntot_cell)
      call put_r3(chk%ieq_met_coef,  ieq_met_coef)
      ieq_rates_ready  = chk%ieq_rates_ready
      ieq_neq_stored   = chk%ieq_neq_stored
      ieq_mbase_stored = chk%ieq_mbase_stored
      ieq_iox_stored   = chk%ieq_iox_stored

      t_phys           = chk%t_phys
      n_steps_accepted = chk%n_steps_accepted
      ledger_family    = chk%ledger_family
      end subroutine attempted_step_checkpoint_restore

      ! ---------------------------------------------------------------- !

      logical function attempted_step_checkpoint_matches(chk, u, f_sp,     &
                             heat, cool, eta, first_mismatch) result(ok)
      ! The round trip is not assumed. Every class-1 item is compared with
      ! the live state, allocation status included, and the first item that
      ! disagrees is named.
      type(attempted_step_checkpoint), intent(in)  :: chk
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng),            intent(in) :: heat,cool,eta
      character(len=*), intent(out), optional :: first_mismatch
      character(len=48) :: bad

      ok = .true.;  bad = ''
      call same_r2(ok, bad, 'u',    chk%u,    u)
      call same_r2(ok, bad, 'f_sp', chk%f_sp, f_sp)
      call same_r1f(ok, bad, 'heat', chk%heat, heat)
      call same_r1f(ok, bad, 'cool', chk%cool, cool)
      call same_r1f(ok, bad, 'eta',  chk%eta,  eta)

      if (chk%carr_present) then
         if (.not. carrier_checkpoint_matches(chk%carr)) then
            ok = .false.
            if (len_trim(bad) .eq. 0) bad = 'carrier_checkpoint'
         endif
      endif

      call same_r1(ok, bad, 'gph_balmer_HI', chk%gph_balmer_HI,            &
                   gph_balmer_HI)
      call same_r1(ok, bad, 'heat_balmer', chk%heat_balmer, heat_balmer)
      call same_r1(ok, bad, 'Jlya_arr',    chk%Jlya_arr,    Jlya_arr)
      call same_r1(ok, bad, 'n2s_arr',     chk%n2s_arr,     n2s_arr)
      call same_r1(ok, bad, 'n2p_arr',     chk%n2p_arr,     n2p_arr)
      call same_r1(ok, bad, 'Sproton_arr', chk%Sproton_arr, Sproton_arr)
      call same_r1(ok, bad, 'Hpe_arr',     chk%Hpe_arr,     Hpe_arr)
      call same_r1(ok, bad, 'Hdx_arr',     chk%Hdx_arr,     Hdx_arr)
      call same_r1(ok, bad, 'jlya_rt_grid',chk%jlya_rt_grid,jlya_rt_grid)
      call same_r1(ok, bad, 'Tdiag',       chk%Tdiag,       Tdiag)
      call same_r1(ok, bad, 'nhidiag',     chk%nhidiag,     nhidiag)
      call same_r1(ok, bad, 'nediag',      chk%nediag,      nediag)
      call same_r1(ok, bad, 'nhiidiag',    chk%nhiidiag,    nhiidiag)
      call same_r1(ok, bad, 'taulya',      chk%taulya,      taulya)
      call same_r1(ok, bad, 'jint_arr',    chk%jint_arr,    jint_arr)
      call same_r1(ok, bad, 'jstar_arr',   chk%jstar_arr,   jstar_arr)
      call same_l (ok, bad, 'jlya_rt_loaded', chk%jlya_rt_loaded,          &
                   jlya_rt_loaded)
      call same_l (ok, bad, 'excited_H_field_ready',                       &
                   chk%excited_H_field_ready, excited_H_field_ready)

      call same_r2(ok, bad, 'nmol_eq',      chk%nmol_eq,      nmol_eq)
      call same_r1(ok, bad, 'NH2_col_lw',   chk%NH2_col_lw,   NH2_col_lw)
      call same_r1(ok, bad, 'f_shield_lw',  chk%f_shield_lw,  f_shield_lw)
      call same_r1(ok, bad, 'k_lw_diss',    chk%k_lw_diss,    k_lw_diss)
      call same_r1(ok, bad, 'p_lw_single',  chk%p_lw_single,  p_lw_single)
      call same_r1(ok, bad, 'p_lw_absorbed',chk%p_lw_absorbed,p_lw_absorbed)
      call same_r1(ok, bad, 'tr_lines_lw',  chk%tr_lines_lw,  tr_lines_lw)
      call same_r1(ok, bad, 'P_H2_eq',      chk%P_H2_eq,      P_H2_eq)
      call same_r2(ok, bad, 'nox_eq',       chk%nox_eq,       nox_eq)
      call same_r1(ok, bad, 'n_o1d_eq',     chk%n_o1d_eq,     n_o1d_eq)
      call same_r1(ok, bad, 'NH2O_col',     chk%NH2O_col,     NH2O_col)
      call same_r1(ok, bad, 'NOH_col',      chk%NOH_col,      NOH_col)
      call same_r1(ok, bad, 'NCO_col',      chk%NCO_col,      NCO_col)
      call same_r1(ok, bad, 'k_co_diss',    chk%k_co_diss,    k_co_diss)
      call same_r1(ok, bad, 'theta_co_shield', chk%theta_co_shield,        &
                   theta_co_shield)
      call same_r2(ok, bad, 'j_h2o_fuv',    chk%j_h2o_fuv,    j_h2o_fuv)
      call same_r2(ok, bad, 'j_oh_fuv',     chk%j_oh_fuv,     j_oh_fuv)
      call same_r2(ok, bad, 'tau_fuv',      chk%tau_fuv,      tau_fuv)
      call same_r1(ok, bad, 'heat_fuv',     chk%heat_fuv,     heat_fuv)
      call same_r1(ok, bad, 'heat_chem',    chk%heat_chem,    heat_chem)
      call same_r2(ok, bad, 'heat_channel_state', chk%heat_channel_state,  &
                   heat_channel_state)
      call same_ir(ok, bad, 'bg_cell',         chk%bg_cell,   bg_cell)
      call same_ir(ok, bad, 'bg_cell_adopted', chk%bg_cell_adopted,        &
                   bg_cell_adopted)
      call same_ir(ok, bad, 'bg_cell_best',    chk%bg_cell_best,           &
                   bg_cell_best)
      call same_i1(ok, bad, 'ieq_nonroot_streak', chk%ieq_nonroot_streak,  &
                   ieq_nonroot_streak)
      call same_int(ok, bad, 'ieq_sweep_state_kind',                       &
                    chk%ieq_sweep_state_kind, ieq_sweep_state_kind)
      call same_int(ok, bad, 'ieq_acc_nprint', chk%ieq_acc_nprint,         &
                    ieq_acc_nprint)
      if (ledger_hash(chk%ieq_marching_ledger(1)) .ne.                     &
          ledger_hash(ieq_marching_ledger(1)) .or.                         &
          ledger_hash(chk%ieq_marching_ledger(2)) .ne.                     &
          ledger_hash(ieq_marching_ledger(2))) then
         ok = .false.
         if (len_trim(bad) .eq. 0) bad = 'ieq_marching_ledger'
      endif
      if (ledger_hash(chk%ieq_steady_iterate_ledger) .ne.                  &
          ledger_hash(ieq_steady_iterate_ledger)) then
         ok = .false.
         if (len_trim(bad) .eq. 0) bad = 'ieq_steady_iterate_ledger'
      endif
      if (ledger_hash(chk%ieq_steady_candidate_ledger) .ne.                &
          ledger_hash(ieq_steady_candidate_ledger)) then
         ok = .false.
         if (len_trim(bad) .eq. 0) bad = 'ieq_steady_candidate_ledger'
      endif
      call same_ir(ok, bad, 'ieq_rate_cell', chk%ieq_rate_cell,            &
                   ieq_rate_cell)
      call same_r1(ok, bad, 'ieq_ne_cell',   chk%ieq_ne_cell,   ieq_ne_cell)
      call same_r1(ok, bad, 'ieq_TK_cell',   chk%ieq_TK_cell,   ieq_TK_cell)
      call same_r1(ok, bad, 'ieq_ntot_cell', chk%ieq_ntot_cell,            &
                   ieq_ntot_cell)
      call same_r3(ok, bad, 'ieq_met_coef',  chk%ieq_met_coef, ieq_met_coef)
      call same_l (ok, bad, 'ieq_rates_ready', chk%ieq_rates_ready,        &
                   ieq_rates_ready)
      call same_int(ok, bad, 'ieq_neq_stored', chk%ieq_neq_stored,         &
                    ieq_neq_stored)
      call same_int(ok, bad, 'ieq_mbase_stored', chk%ieq_mbase_stored,     &
                    ieq_mbase_stored)
      call same_int(ok, bad, 'ieq_iox_stored', chk%ieq_iox_stored,         &
                    ieq_iox_stored)

      if (chk%t_phys .ne. t_phys) then
         ok = .false.
         if (len_trim(bad) .eq. 0) bad = 't_phys'
      endif
      call same_int(ok, bad, 'n_steps_accepted', chk%n_steps_accepted,     &
                    n_steps_accepted)
      call same_int(ok, bad, 'ledger_family', chk%ledger_family,           &
                    ledger_family)

      if (present(first_mismatch)) first_mismatch = bad
      end function attempted_step_checkpoint_matches

      ! ---------------------------------------------------------------- !

      real*8 function attempted_step_checkpoint_checksum(u, f_sp, heat,    &
                             cool, eta) result(h)
      ! A hash over the LIVE state of every class-1 item PLUS the
      ! quantities row 3 rebuilds. The rebuilt ones are in it so that the
      ! membership rule of the checkpoint is a checked property: if a
      ! derived array that is not saved were read before it is rewritten,
      ! the restore-and-rerun hash would differ.
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng),            intent(in) :: heat,cool,eta
      h = 0.0d0
      call hash_r2(h, u);          call hash_r2(h, f_sp)
      call hash_r1f(h, heat);      call hash_r1f(h, cool)
      call hash_r1f(h, eta)
      call hash_r1(h, gph_balmer_HI);  call hash_r1(h, heat_balmer)
      call hash_r1(h, Jlya_arr);       call hash_r1(h, n2s_arr)
      call hash_r1(h, n2p_arr);        call hash_r1(h, Sproton_arr)
      call hash_r1(h, Hpe_arr);        call hash_r1(h, Hdx_arr)
      call hash_r1(h, jlya_rt_grid);   call hash_r1(h, Tdiag)
      call hash_r1(h, nhidiag);        call hash_r1(h, nediag)
      call hash_r1(h, nhiidiag);       call hash_r1(h, taulya)
      call hash_r1(h, jint_arr);       call hash_r1(h, jstar_arr)
      call hash_r2(h, nmol_eq)
      call hash_r1(h, NH2_col_lw);     call hash_r1(h, f_shield_lw)
      call hash_r1(h, k_lw_diss);      call hash_r1(h, p_lw_single)
      call hash_r1(h, p_lw_absorbed);  call hash_r1(h, tr_lines_lw)
      call hash_r1(h, P_H2_eq);        call hash_r2(h, nox_eq)
      call hash_r1(h, n_o1d_eq);       call hash_r1(h, NH2O_col)
      call hash_r1(h, NOH_col);        call hash_r2(h, j_h2o_fuv)
      call hash_r1(h, NCO_col);        call hash_r1(h, k_co_diss)
      call hash_r1(h, theta_co_shield)
      call hash_r2(h, j_oh_fuv);       call hash_r2(h, tau_fuv)
      call hash_r1(h, heat_fuv);       call hash_r1(h, heat_chem)
      call hash_ir(h, bg_cell);        call hash_ir(h, bg_cell_adopted)
      call hash_ir(h, bg_cell_best);   call hash_ir(h, ieq_rate_cell)
      call hash_r1(h, ieq_ne_cell);    call hash_r1(h, ieq_TK_cell)
      call hash_r1(h, ieq_ntot_cell);  call hash_r3(h, ieq_met_coef)
      call hash_i1(h, ieq_nonroot_streak)
      h = h + dble(ieq_sweep_state_kind) + dble(ieq_acc_nprint)
      h = h + ledger_hash(ieq_marching_ledger(1))                          &
            + ledger_hash(ieq_marching_ledger(2))                          &
            + ledger_hash(ieq_steady_iterate_ledger)                       &
            + ledger_hash(ieq_steady_candidate_ledger)
      h = h + t_phys + dble(n_steps_accepted)
      end function attempted_step_checkpoint_checksum

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_perturb_checkpointed_state_for_test(u,     &
                             f_sp, heat, cool, eta)
      ! Moves EVERY item of the checkpoint list to a distinct value, so
      ! that a restore that misses one is caught. Allocates nothing: an
      ! array the configuration leaves unallocated stays unallocated, which
      ! is itself part of the state the restore has to reproduce.
      real*8, dimension(3,1-Ng:N+Ng),          intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),            intent(inout) :: heat,cool,eta
      integer :: j
      real*8  :: s
      s = 1.0d0
      call bump_r2(u, s);      call bump_r2(f_sp, s)
      call bump_r1f(heat, s);  call bump_r1f(cool, s); call bump_r1f(eta, s)
      call bump_r1(gph_balmer_HI, s);  call bump_r1(heat_balmer, s)
      call bump_r1(Jlya_arr, s);       call bump_r1(n2s_arr, s)
      call bump_r1(n2p_arr, s);        call bump_r1(Sproton_arr, s)
      call bump_r1(Hpe_arr, s);        call bump_r1(Hdx_arr, s)
      call bump_r1(jlya_rt_grid, s);   call bump_r1(Tdiag, s)
      call bump_r1(nhidiag, s);        call bump_r1(nediag, s)
      call bump_r1(nhiidiag, s);       call bump_r1(taulya, s)
      call bump_r1(jint_arr, s);       call bump_r1(jstar_arr, s)
      jlya_rt_loaded        = .not. jlya_rt_loaded
      excited_H_field_ready = .not. excited_H_field_ready
      call bump_r2(nmol_eq, s)
      call bump_r1(NH2_col_lw, s);     call bump_r1(f_shield_lw, s)
      call bump_r1(k_lw_diss, s);      call bump_r1(p_lw_single, s)
      call bump_r1(p_lw_absorbed, s);  call bump_r1(tr_lines_lw, s)
      call bump_r1(P_H2_eq, s);        call bump_r2(nox_eq, s)
      call bump_r1(n_o1d_eq, s);       call bump_r1(NH2O_col, s)
      call bump_r1(NOH_col, s);        call bump_r2(j_h2o_fuv, s)
      call bump_r1(NCO_col, s);        call bump_r1(k_co_diss, s)
      call bump_r1(theta_co_shield, s)
      call bump_r2(j_oh_fuv, s);       call bump_r2(tau_fuv, s)
      call bump_r1(heat_fuv, s);       call bump_r1(heat_chem, s)
      call bump_r2(heat_channel_state, s)
      call bump_ir(bg_cell, s);        call bump_ir(bg_cell_adopted, s)
      call bump_ir(bg_cell_best, s);   call bump_ir(ieq_rate_cell, s)
      call bump_r1(ieq_ne_cell, s);    call bump_r1(ieq_TK_cell, s)
      call bump_r1(ieq_ntot_cell, s);  call bump_r3(ieq_met_coef, s)
      if (allocated(ieq_nonroot_streak)) then
         do j = lbound(ieq_nonroot_streak,1), ubound(ieq_nonroot_streak,1)
            ieq_nonroot_streak(j) = ieq_nonroot_streak(j) + 7
         enddo
      endif
      ieq_sweep_state_kind = ieq_sweep_state_kind + 1
      ieq_acc_nprint       = ieq_acc_nprint + 3
      ieq_marching_ledger(1)%n_sweep = ieq_marching_ledger(1)%n_sweep + 11
      ieq_marching_ledger(2)%n_sweep = ieq_marching_ledger(2)%n_sweep + 13
      ieq_steady_iterate_ledger%n_sweep   =                                &
           ieq_steady_iterate_ledger%n_sweep + 17
      ieq_steady_candidate_ledger%n_sweep =                                &
           ieq_steady_candidate_ledger%n_sweep + 19
      ieq_rates_ready  = .not. ieq_rates_ready
      ieq_neq_stored   = ieq_neq_stored + 1
      ieq_mbase_stored = ieq_mbase_stored + 1
      ieq_iox_stored   = ieq_iox_stored + 1
      t_phys           = t_phys + 1.0d0
      n_steps_accepted = n_steps_accepted + 1
      ledger_family    = 3 - ledger_family
      end subroutine attempted_step_perturb_checkpointed_state_for_test

      ! ================================================================ !
      ! 2. THE ACCEPTANCE PREDICATE AT THE ADOPTION BOUNDARY
      ! ================================================================ !

      subroutine attempted_step_evaluate(chk, u, W, T, f_sp, heat, cool,   &
                             dt_g, u_hydro, L_hydro, n_no_chem_root,       &
                             carrier_completed, verd)
      ! THE ADOPTION BOUNDARY. The equations are asked of certification's
      ! physical-step context (one evaluator, one set of tolerances); what
      ! is added here is the thermal-energy identity of the step, which
      ! needs the two marks the marching loop takes across the source step
      ! and a formation-energy reservoir that module does not define.
      type(attempted_step_checkpoint), intent(in) :: chk
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u, W, L_hydro
      ! The state row 2 returned, on which the time-discrete hydrodynamic
      ! row is measured.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u_hydro
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: T, heat, cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,  intent(in) :: dt_g
      integer, intent(in) :: n_no_chem_root
      logical, intent(in) :: carrier_completed
      type(attempted_step_verdict), intent(out) :: verd
      type(cert_step_verdict) :: cv

      call certification_evaluate_physical_step(chk%u, u, u_hydro, W, T,   &
               f_sp, chk%f_sp, dt_g, L_hydro, n_no_chem_root,              &
               carrier_completed, cv)

      verd%accepted  = cv%accepted
      verd%reason    = cv%reason
      verd%operation = cv%operation
      verd%jworst    = cv%jworst
      verd%kworst    = cv%kworst
      verd%measure   = cv%measure
      verd%tol       = cv%tol
      verd%what      = cv%what
      verd%hydro_row_max = cv%hydro_row_max
      verd%hydro_jworst  = cv%hydro_jworst
      verd%hydro_kworst  = cv%hydro_kworst
      verd%hydro_row_above = cv%hydro_row_above
      if (.not. verd%accepted) return

      ! THE THERMAL-ENERGY IDENTITY of the step, evaluated and gated
      ! there.
      call attempted_step_energy_identity(chk, u, W, T, f_sp, heat, cool,  &
                                          dt_g, verd)
      end subroutine attempted_step_evaluate

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_thermal_energy_before_sources(u)
      ! The thermal energy of the grid immediately BEFORE the coupled
      ! source step of rows 7 to 9, volume-weighted and summed over the
      ! physical cells. Called by the marching loop; see
      ! energy_identity for what the pair of marks is for.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      uth_before_sources    = thermal_energy_sum(u)
      marked_before_sources = .true.
      marked_after_sources  = .false.
      end subroutine attempted_step_thermal_energy_before_sources

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_thermal_energy_after_sources(u, dt_cell,   &
                                                             heat, cool)
      ! The thermal energy of the grid immediately AFTER the coupled source
      ! step, and the thermal source that step was given, each summed the
      ! same way. The source is formed with the CELL's interval, because a
      ! local-dt run advances each cell by its own; in phys mode, which
      ! refuses local time stepping, dt_cell is one number.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: dt_cell, heat, cool
      integer :: j
      uth_after_sources = thermal_energy_sum(u)
      ! The injected inconsistency of the tests: a thermal energy the
      ! sources did not produce. Off unless asked.
      if (inject_source_energy_offset .ne. 0.0d0)                          &
         uth_after_sources = uth_after_sources                             &
              *(1.0d0 + inject_source_energy_offset)
      q_ext_dt_sources = 0.0d0
      do j = 1, N
         q_ext_dt_sources = q_ext_dt_sources                               &
              + dt_cell(j)*(heat(j) - cool(j))*dvol(j)
      enddo
      marked_after_sources = .true.
      end subroutine attempted_step_thermal_energy_after_sources

      ! ---------------------------------------------------------------- !

      real*8 function thermal_energy_sum(u) result(s)
      ! sum_j dV_j [ u(3,j) - rho v^2/2 ], the conserved energy row minus
      ! its kinetic part, over the physical cells. The kinetic part is
      ! subtracted rather than carried so that a transfer between thermal
      ! and kinetic energy is a transport contribution and not a source.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      integer :: j
      s = 0.0d0
      do j = 1, N
         if (u(1,j) .gt. 0.0d0) then
            s = s + (u(3,j) - 0.5d0*u(2,j)**2/u(1,j))*dvol(j)
         else
            s = s + u(3,j)*dvol(j)
         endif
      enddo
      end function thermal_energy_sum

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_energy_identity(chk, u, W, T, f_sp, heat,  &
                                 cool, dt_g, verd)
      ! THE THERMAL IDENTITY OF THE ATTEMPTED STEP, summed over the
      ! physical cells:
      !
      !   sum_j dV_j [ u_th^{n+1} - u_th^n ]_j
      !       =  sum_j dV_j dt (heat - cool)_j
      !        +  the transport contribution to the thermal energy
      !
      ! and the residual reported and gated is the difference of the two
      ! sides, which is the closure of the SOURCE STEP alone.
      !
      ! WHY THE RESERVOIR IS NOT ON THE LEFT.
      ! The target thermal-energy closure is written as
      ! Delta u_th + Delta u_form = integral Q_ext, with Q_ext carrying
      ! only exchanges with the radiation field. THIS CODE'S heat AND cool
      ! ARE NOT THAT Q_ext: they are the NET THERMAL source. Term by term,
      ! photoionization deposits h nu - E_th and never h nu, so the
      ! potential the photon paid never passes through the thermal pool;
      ! collisional ionization takes E_th out of the electrons as the coio
      ! term of eval_cool; radiative recombination lets the potential leave
      ! as a photon and keeps only the electron's kinetic energy in reco;
      ! and every collisional reaction heat of the molecular, Penning,
      ! associative, Lyman-Werner and oxygen channels is a difference of
      ! the one species formation-energy table. So heat - cool already
      ! carries every reservoir transfer that exchanges with u_th, and
      ! adding sum_s eps_s Delta n_s on top of it would count each
      ! collisional transfer twice AND charge the gas for the ionization
      ! energy the photons paid. Delta u_form is therefore REPORTED beside
      ! the row as the reservoir's own change, never added to it.
      !
      ! WHY THE TRANSPORT CONTRIBUTION IS ON THE RIGHT. T1.5 constrains a
      ! LOCAL source step at fixed volume and fixed mass; an attempted step
      ! is that source step plus the hydrodynamic stages, the diffusion and
      ! carrier operators, the boundary condition, the conduction stage and
      ! the filter, none of whose thermal contribution is in heat - cool.
      ! That contribution is MEASURED, not modeled: the marching loop
      ! marks the thermal energy on both sides of the source step, the
      ! difference of the marks is what the sources did, and the rest of
      ! the step's thermal change is by definition what transport did.
      !
      ! WITHOUT THE MARKS the transport contribution is not measured, and
      ! the residual falls back to the whole-step imbalance
      ! d_u_th - q_ext_dt, which is not a closed budget and is reported as
      ! such (transport_measured is .false. and the gate is inert).
      type(attempted_step_checkpoint), intent(in) :: chk
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u, W
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: T, heat, cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,  intent(in) :: dt_g
      type(attempted_step_verdict), intent(inout) :: verd
      integer :: i, j
      real*8  :: dv, uth_new, uth_old, ufm_new, ufm_old
      real*8  :: sth, sfm, sq, scal, rho_new, rho_old
      real*8  :: eps(n_species)

      do j = 1, n_species
         eps(j) = species_formation_energy(j)
      enddo

      sth = 0.0d0;  sfm = 0.0d0;  sq = 0.0d0;  scal = 0.0d0
      do j = 1, N
         dv = dvol(j)
         rho_new = W(1,j)
         rho_old = chk%u(1,j)
         ! The conserved energy row minus the kinetic part, in code units.
         uth_new = u(3,j) - 0.5d0*rho_new*W(2,j)**2
         if (rho_old .gt. 0.0d0) then
            uth_old = chk%u(3,j) - 0.5d0*chk%u(2,j)**2/rho_old
         else
            uth_old = chk%u(3,j)
         endif
         ufm_new = 0.0d0
         ufm_old = 0.0d0
         do i = 1, n_species
            ufm_new = ufm_new + eps(i)*f_sp(j,i)
            ufm_old = ufm_old + eps(i)*chk%f_sp(j,i)
         enddo
         ! rho is in units of n0*m_H and the fractions are per hydrogen
         ! nucleus, so rho*n0 is the nucleus density the fractions refer to;
         ! the division puts sum_s eps_s n_s into the code's energy density
         ! unit n0 k_B T0.
         ufm_new = ufm_new*rho_new*n0/erg2eV/(n0*kb_erg*T0)
         ufm_old = ufm_old*rho_old*n0/erg2eV/(n0*kb_erg*T0)
         sth = sth + (uth_new - uth_old)*dv
         sfm = sfm + (ufm_new - ufm_old)*dv
         ! heat and cool are volumetric rates in the code's energy density
         ! unit per code time; the step is dt_g of code time.
         sq  = sq + dt_g*(heat(j) - cool(j))*dv
         ! The row scale of energy_semi_implicit, |E_old| + dt(|heat|+|cool|),
         ! summed over the cells: the same scale the cell-local row was
         ! judged by, so the summed test cannot be stricter than the row
         ! test it aggregates.
         scal = scal + (abs(uth_old)                                       &
                        + dt_g*(abs(heat(j)) + abs(cool(j))))*dv
      enddo

      verd%identity_evaluated = .true.
      verd%d_u_th   = sth
      verd%d_u_form = sfm
      verd%identity_scale = max(scal, energy_scale_floor, as_scale_floor)

      if (marked_before_sources .and. marked_after_sources) then
         verd%transport_measured = .true.
         verd%q_ext_dt = q_ext_dt_sources
         verd%d_u_th_transport = sth                                       &
              - (uth_after_sources - uth_before_sources)
         verd%identity_residual = (uth_after_sources - uth_before_sources) &
              - q_ext_dt_sources
      else
         verd%transport_measured = .false.
         verd%q_ext_dt = sq
         verd%d_u_th_transport = 0.0d0
         verd%identity_residual = sth - sq
      endif
      ! THE GATE. In phys mode the closure of the source step refuses the
      ! attempt, under the energy row's own tolerance and its own scale: the
      ! summed residual is bounded by the sum of the cell residuals, each of
      ! which the row already holds to energy_res_tol of its own scale, so
      ! this test cannot be stricter than the row test it aggregates and no
      ! tolerance was chosen to make a run pass. It is inert while the marks
      ! are absent, because without them there is no measured transport
      ! contribution to subtract and the number is not a closure.
      if (verd%transport_measured .and. run_mode .eq. run_mode_phys) then
         if (abs(verd%identity_residual) .gt.                              &
             energy_res_tol*verd%identity_scale) then
            verd%accepted  = .false.
            verd%reason    = as_reject_source_energy
            verd%operation = as_op_energy
            verd%measure   = abs(verd%identity_residual)                   &
                             /verd%identity_scale
            verd%tol       = energy_res_tol
            verd%jworst    = 0
            verd%kworst    = 3
            verd%what      = 'the thermal energy the source step'//        &
                             ' produced is not the source it was given'
         endif
      endif
      end subroutine attempted_step_energy_identity


      ! ================================================================ !
      ! 3. THE REJECTION POLICY, THE COUNTERS AND THE PROBES
      ! ================================================================ !

      real*8 function attempted_step_reduced_dt(dt_in, reason) result(dt_o)
      ! The factor a rejection applies. 0.5 on an admissibility or a
      ! returned-state refusal, matching the nested positivity path so the
      ! two cannot disagree about what a halving is. On an
      ! integration-error refusal the standard controller factor
      ! 0.9 e**(-1/(p+1)), with p the order of the class whose lower bound
      ! the gate read (last_error_order_p, set by the estimate), floored at
      ! 0.2 so one rejection cannot collapse the step. The order is the
      ! DECIDING class's because e is that class's number: raising it to
      ! the exponent of a class of another order would ask for an interval
      ! that no measurement of this rejection supports.
      ! The estimate steers the reduction only when it is FINITE AND
      ! POSITIVE: a nonfinite estimate carries no information about how
      ! much shorter the interval should be, and raising it to a negative
      ! power would propagate the nonfinite value into the step itself, so
      ! that case takes the halving.
      real*8,  intent(in) :: dt_in
      integer, intent(in) :: reason
      real*8 :: e
      if (reason .eq. as_reject_int_error .and.                            &
          last_error_estimate .gt. 0.0d0 .and.                             &
          ieee_is_finite(last_error_estimate)) then
         e = last_error_estimate
         dt_o = dt_in*max(0.2d0,                                           &
                     0.9d0*e**(-1.0d0/dble(last_error_order_p + 1)))
      else
         dt_o = 0.5d0*dt_in
      endif
      end function attempted_step_reduced_dt

      ! ---------------------------------------------------------------- !

      real*8 function attempted_step_bound_next_dt(dt_cfl, dt_accepted)    &
                      result(dt_o)
      ! ADVISOR DECISION 3: a dt reduced by a rejection bounds the next
      ! step, dt_next = min(dt_CFL, 2 dt_accepted). A run that rejects
      ! nothing never sees this bound bite, because dt_accepted is then the
      ! CFL step of the previous state and the CFL step of two neighboring
      ! states of a marching run does not double in one step.
      real*8, intent(in) :: dt_cfl, dt_accepted
      if (dt_accepted .gt. 0.0d0) then
         dt_o = min(dt_cfl, as_dt_growth*dt_accepted)
      else
         dt_o = dt_cfl
      endif
      end function attempted_step_bound_next_dt

      ! ---------------------------------------------------------------- !

      logical function attempted_step_injected_refusal(op) result(fire)
      ! THE INJECTED FAILURE, one operation at a time. It exists because
      ! the rule that a rejected trial restores the state and advances no
      ! physical accumulation cannot be tested on a run that never rejects
      ! a step, and the run's own rejections depend on its state. Rows 0
      ! and 13 are outside the trial: asking for a refusal there is a
      ! refusal of the request, reported by the caller.
      integer, intent(in) :: op
      fire = .false.
      if (reject_after_operation .le. 0) return
      if (op .ne. reject_after_operation) return
      if (reject_at_step .ge. 0 .and. marching_step .ne. reject_at_step) return
      if (reject_in_pass .gt. 0 .and. err_pass_index .ne. reject_in_pass) &
         return
      if (marching_step .ne. injected_at_count) then
         injected_at_count = marching_step
         injected_served   = 0
      endif
      if (injected_served .ge. reject_leading_attempts) return
      injected_served = injected_served + 1
      fire = .true.
      end function attempted_step_injected_refusal

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_note_step(verd)
      ! One rejection, recorded. The physical accumulations are restored by
      ! the checkpoint; only the attempt statistics move here.
      type(attempted_step_verdict), intent(in) :: verd
      n_step_rejections = n_step_rejections + 1
      if (verd%reason .ge. 0 .and.                                         &
          verd%reason .le. as_reject_source_energy)                        &
         n_reject_by_reason(verd%reason) =                                 &
            n_reject_by_reason(verd%reason) + 1
      if (verd%operation .ge. 0 .and. verd%operation .le. as_op_last)      &
         n_reject_by_operation(verd%operation) =                           &
            n_reject_by_operation(verd%operation) + 1
      end subroutine attempted_step_note_step

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_error_estimate(u_full, f_sp_full,          &
                                               T_full, ne_full,            &
                                               u_half, f_sp_half,          &
                                               T_half, ne_half, est)
      ! THE LOCAL ERROR OF THE STEP THE RUN KEEPS, over its tolerance, ON
      ! EVERY QUANTITY THE STEP ADVANCED and not on the three conservative
      ! rows alone:
      !
      !   e = max over the physical cells, the row classes and the rows of
      !         |q(H) - q(H/2, twice)| / [ (1 - 2^-p(class)) scale(q) ].
      !
      ! The retained trajectory is the full step, whose local error is
      ! larger than the difference by err_full_step_factor = 1/(1 - 2^-p)
      ! OF THE CLASS THE ROW BELONGS TO, because the classes are advanced
      ! by integrators of different order; the derivation and the ladder
      ! that measured each p are at their declaration. The
      ! two half steps ran the WHOLE operator split from the same
      ! checkpoint and over the same interval as the full step, so this is
      ! an estimate of the split at a common endpoint and not of one
      ! operator.
      !
      ! WHY THE SPECIES ROWS ARE HERE. A backward-Euler chemical step can
      ! solve its nonlinear equations to any tolerance and still have a
      ! large temporal error, and the chemistry is where the shortest time
      ! scales of this problem are, so an estimate taken on the
      ! hydrodynamic rows alone reports the accuracy of the smoothest part
      ! of the update. The rows below are the carrier fractions the
      ! chemistry transports and evolves; every other stage reaches the
      ! estimate through the electron fraction and the temperature, which
      ! are the two quantities the opacity, the rate coefficients and the
      ! caloric equation of state are functions of.
      !
      ! THE SCALES ARE THE QUANTITY'S OWN. Each class states its scale
      ! where it is formed; none of them is the flat absolute floor of a
      ! code-unit density.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u_full, u_half
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_full
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_half
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: T_full, T_half
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: ne_full
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: ne_half
      type(integration_error), intent(out) :: est

      integer :: j, k, isp, nrow
      integer :: rows(n_species)
      real*8  :: d, sc, budget, inner, cand
      real*8, dimension(1-Ng:N+Ng) :: nd_h, nH_h, nO_h, nC_h
      real*8, allocatable :: Y_full(:,:), Y_half(:,:)

      est%e = 0.0d0;  est%e_resolved = 0.0d0
      est%worst_class = 0;  est%jworst = 0;  est%kworst = 0
      est%e_class = 0.0d0;  est%j_class = 0;  est%k_class = 0
      est%d_class = 0.0d0;  est%s_class = 0.0d0
      est%class_present = .false.
      est%class_resolved = .true.
      est%resolved = .true.
      est%unresolved_class = 0
      est%deciding_class   = 0
      est%inner_over_difference = 0.0d0

      ! ---- class 1: the conservative hydrodynamic rows ---------------- !
      ! The scale is the one this module has always used on them,
      ! err_atol + err_rtol|u|, in the code units u is held in.
      do j = 1, N
         do k = 1, 3
            d  = abs(u_full(k,j) - u_half(k,j))
            sc = err_atol + err_rtol*abs(u_half(k,j))
            call take_row(est, err_class_hydro, d, sc,             &
                          abs(u_half(k,j)), j, k)
         enddo
      enddo

      ! ---- the element budgets the absolute floors are fractions of ---- !
      ! From the one stoichiometric map, on the HALF-step state, which is
      ! also the state the relative part of every scale is taken on.
      nd_h = u_half(1,:)*n0
      call carrier_element_totals(nd_h, f_sp_half, nH_h, nO_h, nC_h)

      ! ---- class 2: the transported species fractions ----------------- !
      call attempted_step_species_rows(rows, nrow)
      if (nrow .gt. 0) then
         do k = 1, nrow
            isp = rows(k)
            do j = 1, N
               d      = abs(f_sp_full(j,isp) - f_sp_half(j,isp))
               budget = species_budget_fraction(isp, nd_h(j), nH_h(j),     &
                                                nO_h(j), nC_h(j))
               sc = err_rtol*abs(f_sp_half(j,isp))                         &
                    + err_species_floor_of_budget*budget
               call take_row(est, err_class_species, d, sc,         &
                             abs(f_sp_half(j,isp)), j, isp)
            enddo
         enddo
      endif

      ! ---- class 3: the transported element mass fractions ------------ !
      ! Only where the element-diffusion operator runs: without it the
      ! element fractions are carried by the mass row and are not an
      ! unknown of their own. Their budget is the mixture mass, so the
      ! absolute floor is err_species_floor_of_budget of unity.
      if (he_diffusion) then
         allocate(Y_full(1-Ng:N+Ng,1+n_melem), Y_half(1-Ng:N+Ng,1+n_melem))
         call element_mass_fractions(f_sp_full, Y_full)
         call element_mass_fractions(f_sp_half, Y_half)
         do k = 1, 1 + n_melem
            do j = 1, N
               d  = abs(Y_full(j,k) - Y_half(j,k))
               sc = err_rtol*abs(Y_half(j,k)) + err_species_floor_of_budget
               call take_row(est, err_class_element, d, sc,         &
                             abs(Y_half(j,k)), j, k)
            enddo
         enddo
         deallocate(Y_full, Y_half)
      endif

      ! ---- class 4: the electron fraction and the temperature --------- !
      ! The electron fraction is held in the units of f_sp, so that its
      ! difference and a species difference are the same kind of number.
      ! Its budget is the hydrogen the cell holds: a cell whose hydrogen is
      ! fully ionized and whose other elements contribute nothing has that
      ! electron fraction, and no cell of this composition can have more
      ! from hydrogen. The temperature carries no absolute floor of its
      ! own: the energy update keeps it strictly positive, so a relative
      ! scale is defined everywhere and as_scale_floor only keeps the
      ! division safe.
      do j = 1, N
         if (nd_h(j) .gt. 0.0d0) then
            d  = abs(ne_full(j) - ne_half(j))/nd_h(j)
            sc = err_rtol*abs(ne_half(j))/nd_h(j)                          &
                 + err_species_floor_of_budget*nH_h(j)/nd_h(j)
            call take_row(est, err_class_charge_T, d, sc,                  &
                          abs(ne_half(j))/nd_h(j), j, 1)
         endif
         d  = abs(T_full(j) - T_half(j))
         sc = err_rtol*abs(T_half(j)) + as_scale_floor
         call take_row(est, err_class_charge_T, d, sc,              &
                       abs(T_half(j)), j, 2)
      enddo

      ! ---- what the difference can and cannot say --------------------- !
      ! Each of the three passes stops its nonlinear solves at a finite
      ! tolerance, so the difference the estimate measured bounds the
      ! temporal error only to within the sum of the three inner errors.
      ! e_lower is what the macrostep's error is AT LEAST and e_upper what
      ! it can be AT MOST; the verdict is taken from those two and not from
      ! the point value between them, so that neither a rejection nor an
      ! acceptance can rest on a difference that is the solvers' own error.
      do k = 1, n_err_class
         if (.not. est%class_present(k)) cycle
         if (est%s_class(k) .le. 0.0d0) cycle
         inner = inner_error_of_class(k, est%j_class(k),                 &
                                      est%k_class(k), est%q_class(k))
         inner = dble(n_err_passes_of_a_macrostep)*inner
         ! THE DECIDING CLASS is the one whose lower bound is the number
         ! the gate reads, and its order is the one the reduction uses.
         ! A class that carries no difference at all still competes, so
         ! that a run whose every bound is zero names a class rather than
         ! none.
         cand = err_full_step_factor(k)                                    &
                *max(est%d_class(k) - inner, 0.0d0)/est%s_class(k)
         if (est%deciding_class .eq. 0 .or. cand .gt. est%e_lower) then
            est%deciding_class = k
            est%e_lower        = cand
         endif
         est%e_upper = max(est%e_upper, err_full_step_factor(k)            &
              *(est%d_class(k) + inner)/est%s_class(k))
         if (est%d_class(k) .gt. 0.0d0) then
            est%inner_over_difference = max(est%inner_over_difference,      &
                                            inner/est%d_class(k))
            ! THE STATED FRACTION: a class whose three inner
            ! errors are below a tenth of the difference it carries has an
            ! interval narrow enough that the point value describes it.
            if (inner .gt. err_inner_fraction*est%d_class(k)) then
               est%class_resolved(k) = .false.
               if (est%unresolved_class .eq. 0) est%unresolved_class = k
            endif
         endif
      enddo

      ! The worst row overall.
      do k = 1, n_err_class
         if (.not. est%class_present(k)) cycle
         if (est%e_class(k) .gt. est%e) then
            est%e = est%e_class(k)
            est%worst_class = k
            est%jworst = est%j_class(k)
            est%kworst = est%k_class(k)
         endif
      enddo

      ! The test hook, on the leading estimates of the run. It multiplies
      ! the bounds as well as the point value, so a forced rejection is a
      ! rejection of the gate the run reads and not only of the number it
      ! reports.
      if (err_inject_factor .gt. 0.0d0 .and.                               &
          n_error_estimates .lt. err_inject_estimates) then
         est%e       = est%e*err_inject_factor
         est%e_lower = est%e_lower*err_inject_factor
         est%e_upper = est%e_upper*err_inject_factor
         do k = 1, n_err_class
            est%e_class(k) = est%e_class(k)*err_inject_factor
         enddo
      endif

      ! THE VERDICT THE GATE READS. A macrostep whose error is above one
      ! even at the bottom of the interval is refused; one whose error is
      ! below one even at the top is accepted; between the two the
      ! difference does not decide, and the outcome is UNRESOLVED. The
      ! number the gate compares with one is e_resolved: the lower bound,
      ! so that only an error the passes' own tolerances cannot explain can
      ! refuse a step.
      est%e_resolved = est%e_lower
      est%resolved   = (est%e_lower .gt. 1.0d0) .or. (est%e_upper .le. 1.0d0)
      ! The test hook: an inner error large enough to swallow the whole
      ! difference, on this many leading estimates, so that the named
      ! UNRESOLVED outcome is reachable on a run whose solves are tight.
      if (n_error_estimates .lt. force_unresolved_estimates) then
         est%e_lower    = 0.0d0
         est%e_resolved = 0.0d0
         est%e_upper    = max(est%e_upper, 2.0d0)
         est%resolved   = .false.
         do k = 1, n_err_class
            if (est%class_present(k) .and. est%d_class(k) .gt. 0.0d0) then
               est%class_resolved(k) = .false.
               if (est%unresolved_class .eq. 0) est%unresolved_class = k
            endif
         enddo
      endif
      last_error_estimate = est%e_resolved
      if (est%deciding_class .gt. 0)                                       &
         last_error_order_p = err_order_p(est%deciding_class)
      last_error_step     = marching_step
      last_error_cell     = est%jworst
      last_error_row      = est%kworst
      n_error_estimates   = n_error_estimates + 1
      if (.not. est%resolved)                                              &
         n_error_estimates_unresolved = n_error_estimates_unresolved + 1
      end subroutine attempted_step_error_estimate

      ! ---------------------------------------------------------------- !

      subroutine take_row(est, icl, d, sc, q, j, k)
      ! One row of one class, scaled by the retained-step factor, kept if
      ! it is the worst of its class. q is the magnitude of the quantity
      ! itself, kept so that a relative inner tolerance can be turned into
      ! an error of the same kind as the difference. A row whose scale is not positive is
      ! a row of a cell in which the quantity does not exist (a species of
      ! an element the run does not carry); it is skipped rather than
      ! floored, so that an absent species cannot decide a step.
      type(integration_error), intent(inout) :: est
      integer, intent(in) :: icl, j, k
      real*8,  intent(in) :: d, sc, q
      real*8 :: m
      if (sc .le. 0.0d0) return
      est%class_present(icl) = .true.
      m = err_full_step_factor(icl)*d/sc
      if (m .gt. est%e_class(icl)) then
         est%e_class(icl) = m
         est%j_class(icl) = j
         est%k_class(icl) = k
         est%d_class(icl) = d
         est%s_class(icl) = sc
         est%q_class(icl) = q
      endif
      end subroutine take_row

      ! ---------------------------------------------------------------- !

      double precision function species_budget_fraction(isp, nd, nH, nO,   &
                               nC) result(f_max)
      ! THE LARGEST VALUE THE FRACTION OF SPECIES isp CAN TAKE IN A CELL:
      ! the smallest of the budgets of the elements it is made of, each
      ! divided by the nuclei of that element one of its particles holds,
      ! and then by the density the fraction is a fraction of. It is
      ! stoichiometry alone, read from the one map (nuclei_per_particle),
      ! and it is the quantity the absolute floor of a species row is a
      ! fraction of: a species can be no larger than this, so a difference
      ! that is a negligible fraction OF THIS is negligible for that
      ! species whatever the code units make of it.
      integer, intent(in) :: isp
      real*8,  intent(in) :: nd, nH, nO, nC
      integer :: nu
      f_max = huge(1.0d0)
      nu = nuclei_per_particle(isp, ien_H)
      if (nu .gt. 0) f_max = min(f_max, nH/dble(nu))
      nu = nuclei_per_particle(isp, ien_O)
      if (nu .gt. 0) f_max = min(f_max, nO/dble(nu))
      nu = nuclei_per_particle(isp, ien_C)
      if (nu .gt. 0) f_max = min(f_max, nC/dble(nu))
      if (f_max .ge. huge(1.0d0) .or. nd .le. 0.0d0) then
         f_max = 0.0d0
      else
         f_max = f_max/nd
      endif
      end function species_budget_fraction

      ! ---------------------------------------------------------------- !

      double precision function inner_error_of_class(icl, j, krow, q)      &
                               result(e_in)
      ! ONE PASS'S INNER NONLINEAR ERROR IN THE UNITS OF THE ROW IT WOULD
      ! CONTAMINATE. The solvers state their tolerances in their own norms
      ! (attempted_step_note_inner_error), so a relative one is multiplied
      ! here by the magnitude q of the quantity at the row.
      !
      !   hydro       only the total-energy row carries the coupled pair's
      !               temperature error; the mass and momentum rows are
      !               written by the transport stages, which solve nothing
      !               nonlinearly, so nothing of the source step's inner
      !               error reaches them
      !   species     the coupled pair's composition error, absolute in
      !               f_sp already, plus the carrier system's relative row
      !               imbalance on the species' own value
      !   element     the element mass fractions are normalized sums of the
      !               same species, so a composition error of a given size
      !               enters them with a coefficient of order one
      !   charge_T    the electron fraction is in the units of f_sp, so it
      !               takes the composition error directly; the
      !               temperature row takes the relative one on T
      integer, intent(in) :: icl, j, krow
      real*8,  intent(in) :: q
      real*8 :: d_comp, rel_T
      e_in   = 0.0d0
      d_comp = 0.0d0
      rel_T  = 0.0d0
      if (allocated(inner_error_comp_cell)) then
         if (j .ge. lbound(inner_error_comp_cell,1) .and.                  &
             j .le. ubound(inner_error_comp_cell,1)) then
            d_comp = inner_error_comp_cell(j)
            rel_T  = inner_error_T_cell(j)
         endif
      endif
      select case (icl)
      case (err_class_hydro)
         if (krow .eq. 3) e_in = rel_T*q
      case (err_class_species, err_class_element)
         e_in = d_comp + inner_error_carrier_row*q
      case (err_class_charge_T)
         if (krow .eq. 1) then
            e_in = d_comp
         else
            e_in = rel_T*q
         endif
      end select
      end function inner_error_of_class

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_species_rows(rows, nrow)
      ! THE SPECIES ROWS OF THE ESTIMATE: the carrier fractions of the
      ! transported chemistry, which are the f_sp columns the carrier
      ! operator reads and writes (carrier_state builds fc from exactly
      ! these columns and carrier_write_back returns them). The proton is
      ! always among them, transported or not: it is the ionization state,
      ! and every rate coefficient, the electron density and the mean mass
      ! per particle are functions of it. The molecular and oxygen columns
      ! are there when the run carries that chemistry and are identically
      ! zero otherwise.
      integer, intent(out) :: rows(n_species), nrow
      nrow = 0
      nrow = nrow + 1;  rows(nrow) = isp_HII
      if (thereis_mol) then
         nrow = nrow + 1;  rows(nrow) = isp_H2
      endif
      if (thereis_oxychem) then
         nrow = nrow + 1;  rows(nrow) = isp_OH
         nrow = nrow + 1;  rows(nrow) = isp_H2O
         nrow = nrow + 1;  rows(nrow) = isp_CO
      endif
      end subroutine attempted_step_species_rows

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_inner_error_reset()
      ! Called when an error-control transaction starts or restarts, so
      ! that the inner errors compared with a difference are the inner
      ! errors of the passes that produced that difference.
      if (.not. allocated(inner_error_comp_cell))                          &
         allocate(inner_error_comp_cell(1-Ng:N+Ng))
      if (.not. allocated(inner_error_T_cell))                             &
         allocate(inner_error_T_cell(1-Ng:N+Ng))
      inner_error_comp_cell   = 0.0d0
      inner_error_T_cell      = 0.0d0
      inner_error_composition = 0.0d0
      inner_error_temperature = 0.0d0
      inner_error_carrier_row = 0.0d0
      end subroutine attempted_step_inner_error_reset

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_note_inner_error(d_comp, rel_T, carr_row)
      ! What one pass's nonlinear solves left behind, CELL BY CELL, in the
      ! norms their own tolerances are written in: the coupled pair's
      ! remaining composition error as an absolute species fraction and
      ! its temperature error as a relative one (EXHALE_main, the cell
      ! arrays behind csm_err_c and csm_err_T), and the carrier system's
      ! relative row imbalance
      ! (carrier_transport_diagnostics). The worst over the passes of the
      ! transaction is what the estimate is judged against, because the
      ! difference it measures carries the error of every pass.
      !
      ! THE CARRIER NUMBER IS A ROW RESIDUAL and not a solution error: it
      ! bounds the solution error only through the conditioning of that
      ! system, which is not measured here. It enters as a RELATIVE bound
      ! and is reported with the estimate, so a class refused as unresolved
      ! on it can be told from one refused on the coupled pair.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: d_comp, rel_T
      real*8, intent(in) :: carr_row
      integer :: j
      if (.not. allocated(inner_error_comp_cell))                          &
         call attempted_step_inner_error_reset
      do j = 1-Ng, N+Ng
         if (d_comp(j) .gt. inner_error_comp_cell(j))                      &
            inner_error_comp_cell(j) = d_comp(j)
         if (rel_T(j) .gt. inner_error_T_cell(j))                          &
            inner_error_T_cell(j) = rel_T(j)
      enddo
      inner_error_composition = max(inner_error_composition,               &
                                    maxval(d_comp(1:N)))
      inner_error_temperature = max(inner_error_temperature,               &
                                    maxval(rel_T(1:N)))
      if (carr_row .gt. inner_error_carrier_row)                           &
         inner_error_carrier_row = carr_row
      end subroutine attempted_step_note_inner_error

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_read_environment()
      ! The probes and the duty cycle, all off unless asked. The duty cycle
      ! defaults to 20 in phys mode and to off in
      ! init, where a run makes no statement about a trajectory and
      ! therefore has no integration error to bound.
      character(len=64) :: env
      if (run_mode .eq. run_mode_phys) n_err_every = 20
      call get_environment_variable('EXHALE_ERR_EVERY', env)
      if (len_trim(env) .gt. 0) read(env,*) n_err_every
      call get_environment_variable('EXHALE_REJECT_AFTER_OP', env)
      if (len_trim(env) .gt. 0) read(env,*) reject_after_operation
      call get_environment_variable('EXHALE_REJECT_AT_STEP', env)
      if (len_trim(env) .gt. 0) read(env,*) reject_at_step
      call get_environment_variable('EXHALE_REJECT_ATTEMPTS', env)
      if (len_trim(env) .gt. 0) read(env,*) reject_leading_attempts
      call get_environment_variable('EXHALE_INJECT_SOURCE_ENERGY', env)
      if (len_trim(env) .gt. 0) read(env,*) inject_source_energy_offset
      call get_environment_variable('EXHALE_AS_REJECT_IN_PASS', env)
      if (len_trim(env) .gt. 0) read(env,*) reject_in_pass
      call get_environment_variable('EXHALE_AS_ERR_INJECT', env)
      if (len_trim(env) .gt. 0) read(env,*) err_inject_factor
      call get_environment_variable('EXHALE_AS_ERR_INJECT_COUNT', env)
      if (len_trim(env) .gt. 0) read(env,*) err_inject_estimates
      call get_environment_variable('EXHALE_AS_FIXED_DT', env)
      if (len_trim(env) .gt. 0) read(env,*) as_fixed_dt_seconds
      call get_environment_variable('EXHALE_AS_FORCE_UNRESOLVED', env)
      if (len_trim(env) .gt. 0) read(env,*) force_unresolved_estimates
      if (reject_after_operation .gt. 0) then
         if (reject_after_operation .lt. as_op_checkpoint .or.             &
             reject_after_operation .gt. as_op_shapiro) then
            write(*,'(a,i0,a)') ' EXHALE_REJECT_AFTER_OP = ',              &
                 reject_after_operation, ' is OUTSIDE THE TRIAL.'
            write(*,'(a)') '   Operation 0 sets the step and operation'//  &
                 ' 13 reads the adopted state; only'
            write(*,'(a)') '   operations 1 to 12 are part of the'//       &
                 ' attempted step and can be refused.'
            error stop 1
         endif
      endif
      end subroutine attempted_step_read_environment

      ! ---------------------------------------------------------------- !

      subroutine attempted_step_report(n_halvings)
      ! The end-of-run block. The two attempt counts are printed with their
      ! definitions so that they are never added.
      integer, intent(in) :: n_halvings
      integer :: i
      write(*,*)
      write(*,'(a)') '   ATTEMPTED-STEP CONTROLLER'
      write(*,'(a,i0,a,i0,a)') '     outer attempts: ', n_outer_attempts,  &
           ', of which ', n_step_rejections, ' were refused at the'//      &
           ' adoption boundary'
      write(*,'(a,i0,a)') '     steps that needed at least one retry: ',   &
           n_steps_with_rejection, ''
      write(*,'(a,i0,a,i0,a)') '     hydro stage attempts (the nested'//   &
           ' positivity retry, a DIFFERENT count): ', n_steps_attempted,   &
           ' in ', n_halvings, ' halvings'
      if (n_step_exhaustions .gt. 0)                                       &
         write(*,'(a,i0)') '     retry budgets exhausted: ',               &
              n_step_exhaustions
      if (n_step_rejections .gt. 0) then
         ! BY REASON as well as by operation. The two answer different
         ! questions -- where the step was refused, and what refused it --
         ! and an energy failure and an exhausted fixed point are both
         ! refusals of the coupled source step, so the operation alone does
         ! not say which one a run met.
         write(*,'(a)') '     refusals by reason:'
         do i = 0, as_reject_source_energy
            if (n_reject_by_reason(i) .gt. 0)                              &
               write(*,'(a,a,a,i0)') '       ',                            &
                    trim(attempted_step_reason_text(i)), ': ',             &
                    n_reject_by_reason(i)
         enddo
         write(*,'(a)') '     refusals by operation:'
         do i = 0, as_op_last
            if (n_reject_by_operation(i) .gt. 0)                           &
               write(*,'(a,i2,a,a,a,i0)') '       row ', i, '  ',          &
                    trim(attempted_step_operation_text(i)), ': ',          &
                    n_reject_by_operation(i)
         enddo
      endif
      ! B6 category 4: the one unbudgeted accepted correction still inside
      ! the trial. It is reported whether or not it fired, because its being
      ! ACTIVE is what refuses certification, not its count.
      if (n_shapiro_applied .gt. 0) then
         write(*,'(a)') '     UNBUDGETED ACCEPTED CORRECTIONS inside'//    &
              ' the trial (B6 category 4):'
         write(*,'(a,i0,a)') '       row 12, the Shapiro filter: ',        &
              n_shapiro_applied, ' step(s). It alters the adopted'//       &
              ' state with no source term'
         write(*,'(a)') '       behind it; B3b budgets the energy'//       &
              ' it removes.'
         write(*,'(a)') '       NO STATE CERTIFIES WHILE IT IS'//          &
              ' ACTIVE (design section 6).'
      endif
      if (n_error_estimates .gt. 0) then
         write(*,'(a,i0,a,es12.5,a,i0)') '     step-doubling integration'//&
              '-error estimates: ', n_error_estimates, ', last e = ',      &
              last_error_estimate, ' at step ', last_error_step
         write(*,'(a,es12.5,a)') '     seconds spent in them: ',           &
              error_estimate_seconds, ' (the duty-cycle cost)'
         ! HOW MUCH OF THE TRAJECTORY WAS ERROR-CONTROLLED. An estimate
         ! whose inner nonlinear error was not below the difference it
         ! measured judged nothing, and the steps it belongs to were
         ! adopted without a bound on their temporal error.
         write(*,'(a,i0,a,i0,a)') '     of these, ',                       &
              n_error_estimates_unresolved, ' of ', n_error_estimates,     &
              ' were UNRESOLVED: the inner nonlinear error of the passes'
         write(*,'(a)') '       was not below a tenth of the difference'// &
              ' measured, so the macrostep was'
         write(*,'(a)') '       adopted WITHOUT a bound on its'//          &
              ' temporal error.'
      endif
      end subroutine attempted_step_report

      ! ---------------------------------------------------------------- !

      function attempted_step_reason_text(reason) result(s)
      ! Codes 0 to 9 are certification's and are named there, so a refusal
      ! reads the same whichever side of the boundary prints it.
      integer, intent(in) :: reason
      character(len=64) :: s
      select case (reason)
         case (as_reject_int_error);  s = 'the integration error was too big'
         case (as_reject_injected);   s = 'a refusal was injected on demand'
         case (as_reject_source_fixed_point); s =                          &
              'the coupled source step did not reach its fixed point'
         case (as_reject_source_energy); s =                               &
              'the source step did not close its thermal energy'
         case (as_reject_hydro_stage);s =                                  &
              'a Runge-Kutta stage left the admissible set'
         case default;  s = cert_step_reason_text(reason)
      end select
      end function attempted_step_reason_text

      function attempted_step_operation_text(op) result(s)
      integer, intent(in) :: op
      character(len=56) :: s
      select case (op)
         case (as_op_dt);         s = 'the time-step evaluation'
         case (as_op_checkpoint); s = 'the checkpoint'
         case (as_op_hydro);      s = 'the Runge-Kutta stages and the fluxes'
         case (as_op_primitives); s = 'the change of variables'
         case (as_op_diffusion);  s = 'the element diffusion step'
         case (as_op_carriers);   s = 'the carrier transport step'
         case (as_op_excited_H);  s = 'the excited-hydrogen update'
         case (as_op_ioniz_eq);   s = 'the ionization sweep'
         case (as_op_composition); s =                                     &
              'the composition of the coupled source step'
         case (as_op_energy);     s = 'the energy update'
         case (as_op_bc);         s = 'the boundary conditions'
         case (as_op_conduction); s = 'the viscous and conduction step'
         case (as_op_shapiro);    s = 'the Shapiro filter'
         case (as_op_boundary);   s = 'the adoption boundary'
         case default;            s = 'unnamed'
      end select
      end function attempted_step_operation_text

      ! ================================================================ !
      ! 4. HELPERS
      ! ================================================================ !

      real*8 function dvol(j) result(dv)
      ! The cell volume in code units, r^2 dr of the spherical grid, the
      ! same weight certification uses for its volume-averaged companion.
      ! Only ratios of it enter the identity, so the constant 4 pi is left
      ! out.
      integer, intent(in) :: j
      dv = r(j)*r(j)*dr_j(j)
      end function dvol

      ! ---- allocate-and-copy pairs, one per rank and type -------------

      subroutine save_r1f(dst, src)
      ! The same as save_r1 for an array that reaches this module as a
      ! dummy argument of the caller and is therefore not allocatable here.
      real*8, allocatable, intent(inout) :: dst(:)
      real*8, intent(in) :: src(:)
      if (allocated(dst)) deallocate(dst)
      allocate(dst(size(src)))
      dst = src
      end subroutine save_r1f

      subroutine same_r1f(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      real*8, allocatable, intent(in) :: a(:)
      real*8, intent(in) :: b(:)
      if (.not. allocated(a)) then
         call flag(ok, bad, nm);  return
      endif
      if (size(a) .ne. size(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (any(a .ne. b)) call flag(ok, bad, nm)
      end subroutine same_r1f

      subroutine hash_r1f(h, a)
      real*8, intent(inout) :: h
      real*8, intent(in) :: a(:)
      integer :: j
      do j = 1, size(a)
         h = h + dble(j)*a(j)
      enddo
      end subroutine hash_r1f

      subroutine bump_r1f(a, s)
      real*8, intent(inout) :: a(:)
      real*8, intent(inout) :: s
      integer :: j
      do j = 1, size(a)
         s = s + 1.0d0
         a(j) = a(j) + s
      enddo
      end subroutine bump_r1f

      subroutine save_r1(dst, src)
      real*8, allocatable, intent(inout) :: dst(:)
      real*8, allocatable, intent(in)    :: src(:)
      if (allocated(dst)) deallocate(dst)
      if (.not. allocated(src)) return
      allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine save_r1

      subroutine save_grid_state(dst, src)
      ! The conserved-variable array of the grid, saved WITH THE GRID'S
      ! INDEX RANGE. A saver that reallocates from 1 gives back an array in
      ! which cell j sits at index j + Ng, so a reader that indexes the
      ! saved array by the cell number compares cell j with cell j - Ng.
      ! Whole-array assignment and argument association to an explicit-shape
      ! dummy are conformant either way; direct indexing is not, and the
      ! energy identity below indexes.
      real*8, allocatable, intent(inout) :: dst(:,:)
      real*8, intent(in) :: src(:,1-Ng:)
      if (allocated(dst)) deallocate(dst)
      allocate(dst(size(src,1), 1-Ng:N+Ng))
      dst = src(:,1-Ng:N+Ng)
      end subroutine save_grid_state

      subroutine save_grid_fractions(dst, src)
      ! The species-fraction array, grid index first, saved with the grid's
      ! index range for the same reason as save_grid_state.
      real*8, allocatable, intent(inout) :: dst(:,:)
      real*8, intent(in) :: src(1-Ng:,:)
      if (allocated(dst)) deallocate(dst)
      allocate(dst(1-Ng:N+Ng, size(src,2)))
      dst = src(1-Ng:N+Ng,:)
      end subroutine save_grid_fractions

      subroutine save_r2(dst, src)
      real*8, allocatable, intent(inout) :: dst(:,:)
      real*8, intent(in) :: src(:,:)
      if (allocated(dst)) deallocate(dst)
      allocate(dst(size(src,1), size(src,2)))
      dst = src
      end subroutine save_r2

      subroutine save_r3(dst, src)
      real*8, allocatable, intent(inout) :: dst(:,:,:)
      real*8, allocatable, intent(in)    :: src(:,:,:)
      if (allocated(dst)) deallocate(dst)
      if (.not. allocated(src)) return
      allocate(dst(size(src,1), size(src,2), size(src,3)))
      dst = src
      end subroutine save_r3

      subroutine save_i1(dst, src)
      integer, allocatable, intent(inout) :: dst(:)
      integer, allocatable, intent(in)    :: src(:)
      if (allocated(dst)) deallocate(dst)
      if (.not. allocated(src)) return
      allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine save_i1

      subroutine save_ir(dst, src)
      type(ion_rates), allocatable, intent(inout) :: dst(:)
      type(ion_rates), allocatable, intent(in)    :: src(:)
      if (allocated(dst)) deallocate(dst)
      if (.not. allocated(src)) return
      allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine save_ir

      subroutine put_r1(src, dst)
      real*8, allocatable, intent(in)    :: src(:)
      real*8, allocatable, intent(inout) :: dst(:)
      if (.not. allocated(src)) then
         if (allocated(dst)) deallocate(dst)
         return
      endif
      if (.not. allocated(dst))                                            &
         allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine put_r1

      subroutine put_r2(src, dst)
      real*8, allocatable, intent(in)    :: src(:,:)
      real*8, allocatable, intent(inout) :: dst(:,:)
      if (.not. allocated(src)) then
         if (allocated(dst)) deallocate(dst)
         return
      endif
      if (.not. allocated(dst))                                            &
         allocate(dst(size(src,1), size(src,2)))
      dst = src
      end subroutine put_r2

      subroutine put_r3(src, dst)
      real*8, allocatable, intent(in)    :: src(:,:,:)
      real*8, allocatable, intent(inout) :: dst(:,:,:)
      if (.not. allocated(src)) then
         if (allocated(dst)) deallocate(dst)
         return
      endif
      if (.not. allocated(dst))                                            &
         allocate(dst(size(src,1), size(src,2), size(src,3)))
      dst = src
      end subroutine put_r3

      subroutine put_i1(src, dst)
      integer, allocatable, intent(in)    :: src(:)
      integer, allocatable, intent(inout) :: dst(:)
      if (.not. allocated(src)) then
         if (allocated(dst)) deallocate(dst)
         return
      endif
      if (.not. allocated(dst))                                            &
         allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine put_i1

      subroutine put_ir(src, dst)
      type(ion_rates), allocatable, intent(in)    :: src(:)
      type(ion_rates), allocatable, intent(inout) :: dst(:)
      if (.not. allocated(src)) then
         if (allocated(dst)) deallocate(dst)
         return
      endif
      if (.not. allocated(dst))                                            &
         allocate(dst(lbound(src,1):ubound(src,1)))
      dst = src
      end subroutine put_ir

      ! ---- comparisons, allocation status included --------------------

      subroutine same_r1(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      real*8, allocatable, intent(in) :: a(:), b(:)
      if (allocated(a) .neqv. allocated(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (.not. allocated(a)) return
      if (size(a) .ne. size(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (any(a .ne. b)) call flag(ok, bad, nm)
      end subroutine same_r1

      subroutine same_r2(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      real*8, allocatable, intent(in) :: a(:,:)
      real*8, intent(in) :: b(:,:)
      if (.not. allocated(a)) then
         call flag(ok, bad, nm);  return
      endif
      if (size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2)) then
         call flag(ok, bad, nm);  return
      endif
      if (any(a .ne. b)) call flag(ok, bad, nm)
      end subroutine same_r2

      subroutine same_r3(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      real*8, allocatable, intent(in) :: a(:,:,:), b(:,:,:)
      if (allocated(a) .neqv. allocated(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (.not. allocated(a)) return
      if (size(a) .ne. size(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (any(a .ne. b)) call flag(ok, bad, nm)
      end subroutine same_r3

      subroutine same_i1(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      integer, allocatable, intent(in) :: a(:), b(:)
      if (allocated(a) .neqv. allocated(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (.not. allocated(a)) return
      if (size(a) .ne. size(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (any(a .ne. b)) call flag(ok, bad, nm)
      end subroutine same_i1

      subroutine same_ir(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      type(ion_rates), allocatable, intent(in) :: a(:), b(:)
      integer :: j
      if (allocated(a) .neqv. allocated(b)) then
         call flag(ok, bad, nm);  return
      endif
      if (.not. allocated(a)) return
      if (size(a) .ne. size(b)) then
         call flag(ok, bad, nm);  return
      endif
      do j = lbound(a,1), ubound(a,1)
         if (ion_rates_hash(a(j)) .ne. ion_rates_hash(b(j))) then
            call flag(ok, bad, nm);  return
         endif
      enddo
      end subroutine same_ir

      subroutine same_l(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      logical, intent(in) :: a, b
      if (a .neqv. b) call flag(ok, bad, nm)
      end subroutine same_l

      subroutine same_int(ok, bad, nm, a, b)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      integer, intent(in) :: a, b
      if (a .ne. b) call flag(ok, bad, nm)
      end subroutine same_int

      subroutine flag(ok, bad, nm)
      logical, intent(inout) :: ok
      character(len=*), intent(inout) :: bad
      character(len=*), intent(in) :: nm
      ok = .false.
      if (len_trim(bad) .eq. 0) bad = nm
      end subroutine flag

      ! ---- hashes -----------------------------------------------------

      subroutine hash_r1(h, a)
      real*8, intent(inout) :: h
      real*8, allocatable, intent(in) :: a(:)
      integer :: j
      if (.not. allocated(a)) return
      do j = lbound(a,1), ubound(a,1)
         h = h + dble(j)*a(j)
      enddo
      end subroutine hash_r1

      subroutine hash_r2(h, a)
      real*8, intent(inout) :: h
      real*8, intent(in) :: a(:,:)
      integer :: i, j
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            h = h + dble(i + 977*j)*a(i,j)
         enddo
      enddo
      end subroutine hash_r2

      subroutine hash_r3(h, a)
      real*8, intent(inout) :: h
      real*8, allocatable, intent(in) :: a(:,:,:)
      integer :: i, j, k
      if (.not. allocated(a)) return
      do k = lbound(a,3), ubound(a,3)
         do j = lbound(a,2), ubound(a,2)
            do i = lbound(a,1), ubound(a,1)
               h = h + dble(i + 31*j + 977*k)*a(i,j,k)
            enddo
         enddo
      enddo
      end subroutine hash_r3

      subroutine hash_i1(h, a)
      real*8, intent(inout) :: h
      integer, allocatable, intent(in) :: a(:)
      integer :: j
      if (.not. allocated(a)) return
      do j = lbound(a,1), ubound(a,1)
         h = h + dble(j)*dble(a(j))
      enddo
      end subroutine hash_i1

      subroutine hash_ir(h, a)
      real*8, intent(inout) :: h
      type(ion_rates), allocatable, intent(in) :: a(:)
      integer :: j
      if (.not. allocated(a)) return
      do j = lbound(a,1), ubound(a,1)
         h = h + dble(j)*ion_rates_hash(a(j))
      enddo
      end subroutine hash_ir

      real*8 function ion_rates_hash(x) result(h)
      ! EVERY FIELD OF ion_cell_state's ion_rates, enumerated. The
      ! enumeration follows that type: a field added there and not added
      ! here would leave a piece of the sweep's cell state unchecked, so
      ! the two are kept side by side deliberately and this comment says
      ! so. 40 real*8, 7 logical and the integer x_hetr_row, in declaration
      ! order but for the fields added after the enumeration was first
      ! written, which follow it with the next primes; the integer jcell,
      ! the label a diagnostic written from inside the sweep uses, carries
      ! no state and is left out.
      type(ion_rates), intent(in) :: x
      h =        x%P_HI       + 3.0d0*x%P_HeI      + 5.0d0*x%P_HeII       &
        + 7.0d0*x%rchiiB      + 11.0d0*x%rcheiiB   + 13.0d0*x%rcheiiiB    &
        + 17.0d0*x%nh         + 19.0d0*x%nhe       + 23.0d0*x%a_ion_HI    &
        + 29.0d0*x%a_ion_HeI  + 31.0d0*x%a_ion_HeII                       &
        + 37.0d0*x%a_ion_HeITR+ 41.0d0*x%rcheiTR   + 43.0d0*x%A31         &
        + 47.0d0*x%P_HeITR    + 53.0d0*x%q13       + 59.0d0*x%q31a        &
        + 61.0d0*x%q31b       + 67.0d0*x%Q31       + 71.0d0*x%P_H2        &
        + 73.0d0*x%P_H2_di    + 79.0d0*x%P_H2_dd   + 83.0d0*x%P_H2_nd     &
        + 89.0d0*x%k_LW       + 97.0d0*x%T_K       + 101.0d0*x%ntot       &
        + 103.0d0*x%kcx_He0_Hp+ 107.0d0*x%kcx_Hep_H0                      &
        + 109.0d0*x%n_ofam    + 113.0d0*x%n_co                            &
        + 127.0d0*x%x_h2_fix  + 131.0d0*x%x_oh_fix                        &
        + 137.0d0*x%x_h2o_fix + 139.0d0*x%x_hp_fix                        &
        + 163.0d0*x%q31g      + 167.0d0*x%x_heii_fix                      &
        + 173.0d0*x%x_heiii_fix + 193.0d0*x%kcx_Hepp_H0
      if (x%x_h2_fixed) h = h + 149.0d0
      if (x%x_ox_fixed) h = h + 151.0d0
      if (x%x_hp_fixed) h = h + 157.0d0
      if (x%x_heii_fixed)  h = h + 181.0d0
      if (x%x_heiii_fixed) h = h + 191.0d0
      ! The carried He 2^3S level (He 2^3S transport); zero terms where the
      ! option is off, so the hash of every other run is the one it had.
      h = h + 197.0d0*x%x_hetr_fix + 199.0d0*dble(x%x_hetr_row)
      if (x%x_hetr_fixed)  h = h + 211.0d0
      ! The reservoir's H2 partition of the non-ionized hydrogen (a lower
      ! ghost with a base handoff); zero terms everywhere else.
      h = h + 223.0d0*x%x_h2_neutral_partition
      if (x%x_h2_neutral_partition_fixed) h = h + 227.0d0
      end function ion_rates_hash

      real*8 function ledger_hash(g) result(h)
      ! Every field of ionization_equilibrium's ioniz_eq_ledger.
      type(ioniz_eq_ledger), intent(in) :: g
      integer :: i
      h = dble(g%n_sweep) + 3.0d0*dble(g%n_reseed)                        &
        + 5.0d0*dble(g%n_retry) + 7.0d0*dble(g%n_unphys)                  &
        + 11.0d0*dble(g%n_noroot) + 13.0d0*dble(g%n_mol_clamped)          &
        + 17.0d0*dble(g%n_cce_attempt) + 19.0d0*dble(g%n_cce_root)        &
        + 23.0d0*dble(g%n_cce_solve) + 29.0d0*g%cce_seconds               &
        + 31.0d0*dble(g%n_nonfinite) + 37.0d0*dble(g%n_offsimplex)        &
        + 41.0d0*g%viol_worst + 43.0d0*dble(g%streak_peak)                &
        + 47.0d0*dble(g%n_ghost_open)
      do i = 0, 5
         h = h + dble(53 + i)*dble(g%n_mol_info(i))
      enddo
      do i = 1, 6
         h = h + dble(59 + i)*dble(g%acc_n(i)) + dble(71 + i)*g%acc_resmax(i)
      enddo
      do i = 0, 15
         h = h + dble(83 + i)*dble(g%hist_conv(i))                        &
               + dble(107 + i)*dble(g%hist_uncv(i))
      enddo
      end function ledger_hash

      ! ---- perturbation helpers, for the rollback test ----------------

      subroutine bump_r1(a, s)
      real*8, allocatable, intent(inout) :: a(:)
      real*8, intent(inout) :: s
      integer :: j
      if (.not. allocated(a)) return
      do j = lbound(a,1), ubound(a,1)
         s = s + 1.0d0
         a(j) = a(j) + s
      enddo
      end subroutine bump_r1

      subroutine bump_r2(a, s)
      real*8, intent(inout) :: a(:,:)
      real*8, intent(inout) :: s
      integer :: i, j
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            s = s + 1.0d0
            a(i,j) = a(i,j) + s
         enddo
      enddo
      end subroutine bump_r2

      subroutine bump_r3(a, s)
      real*8, allocatable, intent(inout) :: a(:,:,:)
      real*8, intent(inout) :: s
      integer :: i, j, k
      if (.not. allocated(a)) return
      do k = lbound(a,3), ubound(a,3)
         do j = lbound(a,2), ubound(a,2)
            do i = lbound(a,1), ubound(a,1)
               s = s + 1.0d0
               a(i,j,k) = a(i,j,k) + s
            enddo
         enddo
      enddo
      end subroutine bump_r3

      subroutine bump_ir(a, s)
      type(ion_rates), allocatable, intent(inout) :: a(:)
      real*8, intent(inout) :: s
      integer :: j
      if (.not. allocated(a)) return
      do j = lbound(a,1), ubound(a,1)
         s = s + 1.0d0
         a(j)%P_HI  = a(j)%P_HI + s
         a(j)%T_K   = a(j)%T_K + s
         a(j)%ntot  = a(j)%ntot + s
         a(j)%n_co  = a(j)%n_co + s
         a(j)%x_h2_fixed = .not. a(j)%x_h2_fixed
      enddo
      end subroutine bump_ir

      end module attempted_step

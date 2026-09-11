      program Hydro_ioniz

      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use species_table, only: n_mion, mion_fsp, isp_HeI, isp_HeTR
      use Read_input
      use Initialization
      ! What the restart file said it was produced under (parsed by load_IC).
      use IC_load, only: ic_coupling_present, ic_sec_ion_active,          &
                        ic_rec_method, ic_run_mode, ic_run_mode_present,  &
                        ic_t_phys, ic_t_phys_present, ic_certified
      use setup_report
      use eval_time_step
      use energy_semi_implicit
      use caloric_eos, only: caloric_mixture_active, molecular_cell
      use molecular_reaction_heat, only: formation_energy_density
      use utils
      use composition, only: get_species_densities, comp_T_from_p,         &
                             comp_p_from_T, element_ratio_HeH,             &
                             n_cells_he_singlet_clamped
      use binary_element_diffusion, only: element_diffusion_step,          &
                                          relax_element_composition,      &
                                          element_relaxation_converged,   &
                                          element_relaxation_step_budget, &
                                          element_relaxation_failed,      &
                                          species_advection_active,       &
                                          species_advection_begin_step,   &
                                          species_advection_stage,        &
                                          species_advection_project
      use diffusive_photochemistry, only: photochemical_transport_step,   &
                                          carrier_verdict,              &
                                          carrier_interval_covered,     &
                                          carrier_interval_exhausted,   &
                                          carrier_exhausted_record,     &
                                          carrier_history_certifiable,  &
                                          carrier_transport_stop_on_failure,&
                                          carrier_roundoff_limited_record,&
                                          relax_photochemical_composition,&
                                          carrier_relax_outcome_text,    &
                                          carrier_transport_diagnostics, &
                                          carrier_steady_residual,      &
                                          carrier_name,                 &
                                          carrier_drift_location,       &
                                          carrier_co_domain_record,     &
                                          carrier_co_domain_f_dom
      use element_census, only: element_census_state, element_census_take, &
                               element_census_verify,                     &
                               element_census_reservoir
      use lower_column, only: lower_column_solve
      use molecular_infrared_cooling, only: molecular_infrared_init
      use mol_rates, only: h2_thermochemistry_init
      use steady_residual_mod, only: assemble_residual, residual_norms,  &
                                     write_residual_breakdown,          &
                                     flux_spread_of_state, steady_gates_met, &
                                     flux_spread_above_radius,          &
                                     reconstruction_continuation_rhs,   &
                                     residual_row_scale,                &
                                     face_mass_flux_of_state,           &
                                     n_cells_without_chemical_root
      ! The stationary certification of the state the marching route would
      ! declare solved (docs/a2_certification_contract_20260906.md).
      ! THE ATTEMPTED-STEP CONTROLLER (B3a): the checkpoint of the whole
      ! step, the acceptance predicate at the adoption boundary, the
      ! step-size policy of a rejection and the integration-error estimate.
      use attempted_step, only: attempted_step_checkpoint,               &
                               attempted_step_verdict,                   &
                               attempted_step_checkpoint_take,           &
                               attempted_step_checkpoint_restore,        &
                               attempted_step_checkpoint_matches,        &
                               attempted_step_evaluate,                  &
                               attempted_step_thermal_energy_before_sources,&
                               attempted_step_thermal_energy_after_sources, &
                               attempted_step_injected_refusal,          &
                               attempted_step_note_step,                 &
                               attempted_step_reduced_dt,                &
                               attempted_step_bound_next_dt,             &
                               attempted_step_error_estimate,            &
                               attempted_step_report,                    &
                               attempted_step_read_environment,          &
                               attempted_step_reason_text,               &
                               attempted_step_operation_text,            &
                               n_step_retry_max, n_outer_attempts,       &
                               n_step_rejections, n_step_exhaustions,    &
                               n_steps_with_rejection,                   &
                               n_err_every, error_estimate_seconds,      &
                               last_error_estimate,                      &
                               as_op_checkpoint, as_op_hydro,            &
                               as_op_primitives, as_op_diffusion,        &
                               as_op_carriers, as_op_excited_H,          &
                               as_op_ioniz_eq, as_op_composition,         &
                               as_op_energy, as_op_bc, as_op_conduction, &
                               as_op_shapiro, as_op_boundary,           &
                               as_fixed_dt_seconds, err_pass_index,     &
                               attempted_step_checkpoint_checksum,      &
                               as_reject_injected,                      &
                               as_reject_carrier,                        &
                               as_reject_hydro_stage, as_reject_int_error, &
                               as_reject_energy,                         &
                               as_reject_source_fixed_point,             &
                               integration_error, err_class_name,        &
                               n_err_class, err_inner_fraction,          &
                               err_order_p,                              &
                               attempted_step_inner_error_reset,         &
                               attempted_step_note_inner_error,          &
                               n_error_estimates_unresolved
      use certification, only: cert_report, cert_context_stationary,     &
                               cert_evaluated, cert_unavailable,         &
                               certification_evaluate,                   &
                               certification_report_write,               &
                               certification_entry_index,                &
                               certification_note_stationarity_claim,    &
                               certification_stop_uncertified
      use viscous_conduction, only: transport_active, viscous_conduction_step,&
                                    n_conduction_floor_hits,                &
                                    n_conduction_floor_cells,               &
                                    conduction_floor_first_step,            &
                                    conduction_floor_last_step
      use steady_newton, only: neq_newton, pack_U, unpack_U, newton_residual, &
                               eval_residual, frozen_residual,              &
                               build_banded_jac, band_matvec,               &
                               kl_jac, ku_jac, solve_steady_ptc,           &
                               solve_steady_jfnk, set_base_fix,          &
                               set_transported_species_rows,             &
                               pack_species_rows, cell_state_scales,      &
                               cell_row_scales, nvar_jac, n_species_rows, &
                               read_species_unknown_space_controls,       &
                               freeze_species_unknown_box, jv_product,    &
                               replay_distance,                           &
                               gate_rnorm_accepted, gate_fspread_accepted
      use Conversion
      use ionization_equilibrium
      use utils_ion_eq, only: ec_prof_on, ec_t, ec_name
      use utils_ion_eq, only: write_cool_breakdown_eq, write_heat_breakdown_eq, &
                              n_cells_he_rec_photoionization_dominant,      &
                              he_rec_photoionization_ratio_max,             &
                              he_rec_dominant_ratio
      use lya_rt, only: lya_rt_allocate_arrays
      use excited_hydrogen, only: excited_H_update, write_excited_H,      &
                                  excited_H_allocate_arrays
      use Reconstruction_step
      use RK_integration
      use BC_Apply
      use output_write
      use post_processing
      use ionization_equilibrium
      use newton_solver, only: nt_calls, nt_fallback   ! Task 2 usage counters
      use low_mach_dissipation, only: low_mach_damping_active,          &
                                      contact_mode_dissipation_magnitude
      
      implicit none

      ! A1 elemental census taken around the restart equilibration
      ! (element_census; inert unless EXHALE_ELEMENT_ASSERT is set).
      type(element_census_state) :: cen_main

      ! Logical variables
      logical :: l_isnan = .false.
      logical :: in_plm_stage = .false.   ! true while in the stage-1 (PLM) phase of a two-stage PLM->WENO3 run
      ! PLM -> WENO3 continuation ("Reconstruction continuation:", off by
      ! default) and the fixed-state operator-difference probe that measures
      ! what the hand-off does to the equations (EXHALE_OPDIFF).
      logical :: in_lambda_ramp = .false. ! true while lambda is being walked from 0 to 1
      integer :: lam_step_count = 0       ! marching steps taken since the ramp started
      integer :: lam_good       = 0       ! consecutive steps accepted by the dtu test
      integer :: n_lambda_backoff = 0     ! steps at which the dtu test held lambda back
      real*8  :: dtu_ref        = -1.0d0  ! dtu at the start of the ramp, the reference the test uses
      real*8, dimension(:,:), allocatable :: R_plm, R_weno, u_probe
      integer :: opdiff_mode = 0

      ! Mass-flux level-stability tracking (ring buffer over N_stall steps)
      real*8, allocatable :: lev_hist(:)
      real*8  :: lev_now, lev_rel
      integer :: lev_count
      logical :: is_level_stable

      ! Step at which the staged secondary-ionization coupling was switched on
      ! (-1 = not yet); stops are held for N_stall steps after the flip so the
      ! coupling has time to feed back into the hydro before a stop is accepted.
      integer :: sec_flip_step = -1

      ! Descending-crossing guard on the du-threshold triggers. du is the radial
      ! spread of rho*v*r^2, so a SMALL du is evidence of relaxation only for a
      ! state the marching loop has actually relaxed. On a freshly generated IC
      ! (cold hydrostatic, transonic, Wind-AE, Parker) the profile is analytic
      ! and smooth, and its du measures the smoothness of the formula rather than
      ! the wind: measured du(1) = 1.85e-5 for the WASP-121b transonic IC, which
      ! tripped both the Newton hand-off and the secondary-ionization flip on
      ! step 1 and ran away to NaN. Each of these triggers is therefore armed
      ! only once du has been seen at or above its own threshold, and fires only
      ! on the way back down. A state read back with "Load IC? True" was relaxed
      ! by the marching loop of the run that wrote it, so its du is meaningful
      ! and both triggers start armed. Cold hydrostatic starts arm on step 1
      ! anyway (measured du(1) = 0.37 for the WASP regression cases, 7.9 for
      ! mol_base_handoff), so their behavior is unchanged.
      logical :: du_stop_armed   = .false.   ! du < du_th convergence stop
      ! Set when the JFNK steady solver returned info=0 and its solution is
      ! what ends the run. It is a stop condition of its own -- the marching
      ! loop tests it directly -- because a Newton finish does not have to
      ! satisfy the du criterion: it hands over at du < newton_du_switch, two
      ! decades above du_th, and converges on the residual and flux gates
      ! instead.
      logical :: newton_finished = .false.
      ! The JFNK finish converged and its RETURNED state does not meet the
      ! gate its last iterate met (info = 2, B5pre). The state is written as
      ! the final state with certified=F and the run exits 2: a stationarity
      ! claim that failed certification.
      logical :: newton_returned_uncertified = .false.
      logical :: du_newton_armed = .false.   ! JFNK hand-off, and the secondary-
                                             ! ionization flip sharing its test
      ! Step at which the JFNK hand-off last returned info /= 0 (-1 = never),
      ! and how many times it has been offered a state. The hand-off is held
      ! for N_stall steps after a failure so the marching has moved the wind
      ! before the solver sees it again.
      integer :: newton_fail_step = -1
      integer :: newton_attempts  = 0
      integer, parameter :: n_newton_attempt_max = 3

      ! Residual-based convergence monitor (Resid tol option). resid_max and
      ! fspread_gate are the two gate quantities as last measured, held
      ! between monitor evaluations so the stop report quotes the numbers the
      ! decision was taken on rather than re-measuring them afterwards.
      real*8, dimension(:,:), allocatable :: Rres
      real*8  :: resid_c(3), resid_max, flux_spread
      real*8  :: fspread_gate
      ! Reporting-only companions of the gate window: the same spread taken
      ! from further in. The gate itself is unchanged -- see
      ! flux_spread_above_radius in steady_residual.f90 for why they are
      ! printed.
      real*8  :: fspread_103, fspread_110
      logical :: gates_met_now = .false.
      ! THE ACCEPTANCE LEDGER OF THE LAST EQUILIBRIUM SWEEP THE RUN HELD,
      ! and the certification the run makes from it. Every ioniz_eq call that
      ! solves the composition of a state the run then keeps -- the initial
      ! one, the restart equilibration, the marching loop's, and the steady
      ! route's refresh after each solve -- writes it, so the count of cells
      ! without a chemical root that reaches a certification belongs to the
      ! sweep that produced the state being judged and not to some earlier
      ! one. A Jacobian probe or a line-search trial does not write it: those
      ! sweeps are tagged as candidates and the solver keeps their ledger
      ! itself.
      type(ioniz_eq_ledger) :: last_sweep
      type(cert_report)     :: cert_now
      integer               :: n_no_root_now = 0
      logical               :: cert_first = .true.
      logical               :: cert_was_certified = .false.
      character(len=48)     :: cert_label

      ! P54/(AD) mass-flux time series (env EXHALE_P54_TS = step interval,
      ! 0 = off).  One line per sampled step: the Riemann FACE mass flux
      ! F r^2 at four fixed radii, the cell-centred rho v r^2 at the same
      ! cells, and the two window spreads item (AD) quotes.  Diagnostic
      ! only -- it changes no state and is written after the step.
      integer :: p54_ts_every = 0
      integer :: p54_ts_unit  = -1
      logical :: p54_ts_open  = .false.
      integer, dimension(4) :: p54_ts_face = -1
      real*8,  dimension(4), parameter :: p54_ts_r = (/1.005d0, 1.03d0,   &
                                                       1.10d0,  1.20d0/)

      ! Composition equilibration of a state read from file (see
      ! equilibrate_loaded_composition below). On by default; set
      ! EXHALE_RELOAD_EQ=0 to reproduce a run of the previous code.
      logical :: reload_equilibrate = .true.
      integer :: reload_eq_sweeps   = 0
      real*8  :: reload_eq_dnpart   = 0.0d0

      ! HOW MANY OUTER ITERATIONS A STATIONARY SOLVE OF THIS RUN MAY TAKE.
      ! One number for every route that enters
      ! steady_wind_with_element_diffusion -- the direct steady route, the
      ! marching hand-off and the stationary restart -- because the budget
      ! is a property of the run and not of the door it came in by: the
      ! three routes hand the same solver the same system, and the state
      ! they hand it is not more or less solvable for having been reached by
      ! marching. It is a stop of last resort, far above the iteration count
      ! any converging solve takes (MEASURED: 13 on wasp_full_newton).
      integer, parameter :: jfnk_outer_iterations_default = 3000
      ! WHETHER A CAP WAS NAMED FOR THE RUN, and how many stationary solves
      ! the run has entered. EXHALE_JFNK_MAXIT names an outer-iteration cap
      ! that steady_newton applies to EACH solve, so a run that enters the
      ! solve more than once -- the marching hand-off returns to marching on
      ! info = 1 and offers the hand-off again N_stall steps later -- would
      ! spend the cap once per entry and the number measured under it would
      ! be the sum of several solves from different states. Under a named
      ! cap the run therefore takes its FIRST stationary solve and no other,
      ! so that one measurement is one solve.
      logical :: jfnk_run_cap_named       = .false.
      integer :: n_stationary_solves_run  = 0
      character(len=32) :: jfnk_cap_env

      ! Complete update-map diagnostic (env EXHALE_UPDATE_MAP, off by
      ! default). See run_update_map_report below.
      integer :: upmap_n = 0
      real*8  :: upmap_f(8) = 0.0d0
      real*8  :: upmap_scale = 1.0d0
      real*8, dimension(:,:), allocatable :: u_um0, u_umA, u_umB, u_umC,   &
                                             u_umD, u_umSave, R_um
      real*8, dimension(:,:), allocatable :: f_sp_umSave
      real*8  :: upmap_dfsp = 0.0d0
      integer :: upmap_unit = -1
      logical :: upmap_file_open = .false.
      integer :: upmap_done = 0
      real*8  :: upmap_bc_interior = 0.0d0
      ! P158: the same comparison for the CARRIER row -- the one row the
      ! conserved-variable map above does not carry.
      real*8, dimension(:), allocatable :: nH2_um0, R_carr_um, T_carr_um
      real*8, dimension(:,:), allocatable :: res_um_all, terms_um_all
      logical :: upmap_carrier = .false.
      real*8, dimension(:), allocatable :: y_um0, nrho_um0
      real*8  :: upmap_crmax = 0.0d0
      integer :: upmap_cj = 0, upmap_cic = 1

      ! Newton-residual self-test scratch (EXHALE_NEWTON_TEST hook)
      real*8, allocatable :: Yvec(:), Fvec(:)
      real*8, dimension(:,:), allocatable :: f_sp_test
      ! Residual-determinism self-test scratch (EXHALE_RESID_DETERMINISM hook):
      ! F(Y) must not depend on which states were evaluated before it.
      real*8, allocatable :: Fa(:), Fb(:), Fc(:), Ytry1(:), Ytry2(:)
      real*8, allocatable :: Ytry3(:), Jvd(:), vdir(:), Ddet(:), Drdet(:)
      real*8, dimension(:,:), allocatable :: f_sp_w
      integer :: nbitdiff_b, nbitdiff_c, nvj_det, nbitdiff_ctrl
      real*8  :: dmax_b, dmax_c, det_rs, det_sc_ctrl
      ! The replay's distance in each row class, and the cell and the slot
      ! that carry it: three numbers say where a residual moved, one does
      ! not.
      real*8  :: det_cl_ctrl(2), det_cl_b(2), det_cl_c(2)
      integer :: det_jc_ctrl(2), det_jc_b(2), det_jc_c(2)
      integer :: det_kc_ctrl(2), det_kc_b(2), det_kc_c(2)
      integer :: icl
      logical :: jv_ok_det
      character(len=16) :: det_clname
      ! Banded-Jacobian self-test scratch (EXHALE_JAC_TEST hook)
      real*8, allocatable :: abjac(:,:), rdir(:), Jr(:), dFD(:), F0f(:)
      real*8, allocatable :: Dscale(:)
      real*8, dimension(:), allocatable :: heat0, cool0, npart0
      real*8 :: jac_eps, jac_err
      ! --- lightweight phase profiler (gated by env EXHALE_PROFILE=1) ---
      logical :: do_profile = .false.
      integer :: jj_ec
      real*8  :: tp_step0, tp_a, tp_ion = 0.0d0, tp_hyd = 0.0d0, tp_tot = 0.0d0
      ! finer phases of the marching step (EXHALE_PROFILE=1): the three RK
      ! stages together, the carrier transport operator, the semi-implicit
      ! energy update, and everything else by difference.
      real*8 :: tp_hydro = 0.0d0, tp_trans = 0.0d0, tp_energy = 0.0d0, tp_b = 0.0d0
      character(len=8) :: prof_env
      ! Deterministic step cap (env EXHALE_MAXSTEPS=N): stop after N steps and write
      ! output. Used to compare serial vs parallel runs at an identical step.
      integer :: max_steps = 0
      
      ! Integers variables
      integer :: j,k,im
      
      ! Timing variables
      integer :: n_hrs 
      integer :: n_min 
      real*8  :: start, finish 
      real*8  :: exec_time
      real*8  :: n_sec
      real*8  :: dum
      
      ! Momentum variables
      real*8, dimension(:), allocatable :: mom

      ! Row index of the conserved state in the state-change infinity norm
      ! (dtu), which is taken row by row so that a row whose entry state is
      ! exactly zero can be handled where the ratio is formed.
      integer :: krow

      ! Convergence stall detection (ported from ATES_extended)
      real*8  :: du_prev
      integer :: stall_count

      ! Optional base-cell startup diagnostic (env EXHALE_DIAG_BASE=1): dumps the
      ! first few cells (r, n, v, T, heat, cool, dC/dT sign) for the first steps,
      ! to trace IC-startup transients (e.g. warm-Parker base breakdown). Off by
      ! default => zero effect on normal runs.
      logical            :: diag_base = .false.
      character(len=64)  :: diag_env

      ! Decoupled outer iteration over the excited-H (H n=2) feedback
      real*8  :: exc_rel

      ! The damped Picard alternation of the stationary route
      ! (steady_wind_with_element_diffusion): the composition movement of
      ! one pass, reported and never an acceptance, and the under-relaxation
      ! factor of the element update.
      real*8  :: comp_drift, comp_omega, elem_drift
      ! The wind the element composition relaxes in: the face mass flux of
      ! the accepted state, read from the mass row that state assembled and
      ! handed to the element operator, which has no way of its own to reach
      ! the steady residual.
      real*8, dimension(:), allocatable :: Frho_elem
      ! Drift and step count of the carrier relaxation, the third
      ! participant of the same damped Picard iteration.
      real*8  :: carrier_drift
      real*8  :: crc_max, crc_vol, crc_leg
      ! HOW FAR ONE FIXED-WIND CARRIER PASS MAY MOVE THE COMPOSITION, as a
      ! fraction of the largest H2 mixing ratio of the entry state, enforced
      ! on the accepted step by relax_photochemical_composition. 1e-2 is the
      ! measured boundary of the steady solve's tolerance on the hot Uranus
      ! hand-off: the wind solve still accepts the handed-over state at 0.2
      ! and stops accepting at 0.4, so this stands a factor 20 inside it
      ! (docs/p50_carrier_wind_alternation.md). The outer iteration shortens
      ! it on a pass that failed to move the joint measure, and leaves this
      ! value where a later entry reads the run's own setting.
      !
      ! EXHALE_CARRIER_TRUST overrides it and EXHALE_OUTER_PASSES the pass
      ! cap, so the alternation can be run as a CONTINUATION that walks the
      ! H2 front to where the carrier equation puts it (section 163); both
      ! are diagnostics and neither is an input key.
      real*8  :: carrier_trust = 1.0d-2
      integer :: outer_pass_cap = 20

      integer :: crc_j, crc_ic, cdl_j, cdl_ic
      integer :: it_diff, kd, kc

      ! Admissibility record of the low-Mach contact-mode dissipation: how big
      ! the artificial stress got against the physical momentum flux it was
      ! added to, and how far out its Mach gate stayed open. Only written when
      ! "Low-Mach damping" is on.
      real*8  :: lowmach_ratio, lowmach_r_peak, lowmach_r_gate

      ! Temporal step (global) and cell-by-cell pseudo-time steps
      real*8 :: dt
      real*8, dimension(:), allocatable :: dt_loc

      ! Positivity step control (see the retry_step loop). The three
      ! counters themselves live in global_parameters, so the attempted-step
      ! controller can checkpoint them as attempt statistics and assert that
      ! a restore does not move them back. n_dt_halve_max bounds the
      ! bisection: 20 halvings is a factor 1e-6 on dt, past which no step
      ! size repairs the state. It stays 20 while the OUTER controller's cap
      ! is 8 (n_step_retry_max), deliberately: the inner criterion is the
      ! admissibility of one Runge-Kutta stage and its cost is one stage,
      ! not a whole step, and the two counts are reported separately and
      ! never added (advisor decision 2).
      integer :: n_dt_halve, n_dt_halvings, n_steps_dt_halved
      integer, parameter :: n_dt_halve_max = 20

      ! Did the local first-order flux correction restore the admissible set
      ! for the stage it was applied to? False hands the stage to the dt
      ! bisection above.
      logical :: flux_corr_ok
      ! Did this attempt of the step survive every stage? Read at the exit of
      ! retry_step, which reaches the same statement from an accepted step
      ! and from a step abandoned after n_dt_halve_max bisections.
      logical :: step_accepted

      ! ---- THE ATTEMPTED-STEP CONTROLLER (B3a) ----
      ! as_chk is the checkpoint one outer attempt is taken from and
      ! restored to; as_entry and as_full are the two the step-doubling
      ! estimate needs (the state the step began at, and the result of the
      ! full step held aside while the two halves are taken).
      type(attempted_step_checkpoint) :: as_chk, as_entry, as_full
      type(attempted_step_verdict)    :: as_verd
      type(carrier_verdict)           :: as_cv
      real*8, dimension(:,:), allocatable :: L_hyd, u_err_full, u_hyd
      ! THE FULL-PASS STATE THE COMPARISON IS TAKEN AGAINST, on every row
      ! class of the estimate: the conservative rows, the species
      ! fractions, and the temperature and electron density operation 3
      ! derives from the two. The derived pair is held here rather than
      ! rebuilt from the checkpoint because it is the state the full pass
      ! actually left, which is what the estimate is a statement about.
      real*8, dimension(:,:), allocatable :: f_sp_err_full
      real*8, dimension(:),   allocatable :: T_err_full, ne_err_full
      type(integration_error) :: as_est
      ! THE MACRO-INTERVAL H OF THE ERROR-CONTROL TRANSACTION. Every pass
      ! of the step-doubling comparison covers the SAME interval: the full
      ! pass runs one step of dt_macro and the two half passes one step of
      ! 0.5*dt_macro each, so all three end at t + H. Outside a sampled
      ! macrostep it is the single step's own interval, unchanged by a
      ! retry, and only the adoption boundary reads it.
      real*8, dimension(:), allocatable   :: dtloc_macro
      real*8  :: dt_macro, dt_cfl_now, dt_reduced, dt_last_accepted
      real*8  :: t_err0, e_err
      ! The carrier substep's own diagnostics, read for the inner
      ! nonlinear error the estimate is judged against.
      integer :: ct_steps, ct_nlim
      real*8  :: ct_resid, ct_worst
      ! The coupled pair's remaining error at each cell, in the two norms
      ! its tolerances are written in.
      real*8, dimension(:), allocatable :: csm_inner_c, csm_inner_T
      real*8  :: csm_inner_fac
      integer :: n_step_retry, i_err, n_err_passes, as_inject_op, i_ecl
      ! Retries of the WHOLE macrostep. A rejection inside a sampled
      ! macrostep cannot be repaired by shortening one pass, because the
      ! three passes would then no longer share an endpoint, so it restarts
      ! the transaction from as_entry at a shorter H and this counter, not
      ! n_step_retry, carries the retry budget there.
      integer :: n_macro_retry, n_retry_now
      character(len=20) :: err_pass_name
      ! A refusal raised INSIDE the trial by an operation that knows its own
      ! reason (the coupled source step of rows 7 to 9), as opposed to the
      ! injected refusals, which only know where they fired. Zero means the
      ! trial ran to the adoption boundary.
      integer :: as_trial_reason, as_trial_op
      ! THE RETRY HISTORY OF ONE STEP, printed when the budget is exhausted:
      ! the dt of each attempt and the operation and reason that refused it,
      ! so the exhaustion says how the step was narrowed and whether the same
      ! thing refused it every time.
      real*8  :: as_hist_dt(0:n_step_retry_max)
      integer :: as_hist_reason(0:n_step_retry_max)
      integer :: as_hist_op(0:n_step_retry_max)
      integer :: je_err, ke_err, k_as
      integer :: nro_last, nro_total, nro_sub
      integer :: nex_n, nex_first, nex_last, nex_j, nex_ic
      real*8  :: nex_ratio, nex_phys, nex_full
      character(len=48) :: as_bad_item
      logical :: as_carrier_ok, dt_bounded_by_rejection
      real*8  :: as_frac_done
      integer :: as_nsub, as_carrier_status
      ! The ledger family in force around a stationary solve, which is
      ! continuation and writes to the initialization family.
      integer :: ledger_family_held
      ! EXHALE_STEP_CLOCK=1: one line per step naming the accepted global dt
      ! and the clock after it, so that the elapsed time a run reports can be
      ! checked against the sum of the steps it took.
      logical :: trace_step_clock = .false.
      ! EXHALE_REJECT_STEP=<n>: refuse the first attempt of the step taken at
      ! count = n through the positivity retry path, so that the rejection
      ! branch can be exercised on demand. -1 (the default) never fires.
      integer :: reject_step_probe = -1

      ! Value of n_faces_flux_positivity_limited when the current attempt at
      ! the step began. The repairs of an attempt that the dt bisection
      ! discards are not part of the trajectory, so only the difference taken
      ! at the acceptance point of the step is added to
      ! n_faces_flux_positivity_limited_accepted.
      integer :: n_faces_positivity_at_attempt

      ! Mdot value, and the dimensional mass flux 4 pi rho v r^2 [g/s] it is
      ! the logarithm of; the logical records whether that flux is an outflow
      ! at the cell it is read from.
      real*8 :: Mdot, Mdot_cgs
      logical :: outflow_at_mdot_cell
      
      ! Vectors of thermodynamical variables
      real*8, dimension(:), allocatable :: rho,v,E,p,T,cs
      real*8, dimension(:), allocatable :: heat,cool
      real*8, dimension(:), allocatable :: eta
      real*8, dimension(:), allocatable :: nhi,nhii
      real*8, dimension(:), allocatable :: nhei,nheii,nheiii
      real*8, dimension(:), allocatable :: nheiTR
      ! Metal ion densities, 2D: column i = ion i of the species_table
      real*8, dimension(:,:), allocatable :: nm
      real*8, dimension(:), allocatable :: ne,n_tot
      real*8, dimension(:,:), allocatable :: f_sp  ! 1-6: H/He(+HeITR), 7-9: CI/II/III, 10-12: OI/II/III, 13-15: NI/II/III, 16-18: MgI/II/III
       
      ! Conservative and primitive vectors
      real*8, dimension(:,:), allocatable :: u,u1,u2,u_old
      real*8, dimension(:,:), allocatable :: W,WL,WR
          
      ! Flux and source vectors
      real*8, dimension(:,:), allocatable :: dF,S

      !------------------------------------------------!
      ! THE COUPLED SOURCE STEP (rows 7 to 9; b1 T1.4 to T1.9, T2.1).
      !
      ! csm_T_tol / csm_comp_tol are the fixed-point tolerances of the
      ! composition-temperature pair, NOT accuracy tolerances of either
      ! solve: each of the two solves keeps its own residual test (the energy
      ! row at 1e-9 of its scale, the composition at the sweep's own
      ! acceptance classes) and these say when the two have stopped moving
      ! each other.
      !
      ! WHAT THEY BOUND IS THE ERROR OF THE ACCEPTED STATE, its distance
      ! to the fixed point, and not the last step towards it.  The pass
      ! sequence decays geometrically with a ratio theta the loop measures
      ! from the sequence itself, so the error a pass leaves is
      ! |theta/(1 - theta)| times that pass's own increment, a factor of
      ! 0.07 to 0.28 at the measured |theta| of 0.07 to 0.22.  The loop
      ! therefore stops when the ESTIMATED ERROR, carried into the test's
      ! own norm by csm_err_safety, is within these two tolerances
      ! (coupled_pair_geometric_error_estimate), each in the norm its own
      ! tolerance is written in.  On a pass whose sequence gives no
      ! estimate, the first two passes of a step or a pair the guards
      ! refuse, the increment is the only bound the loop can state and the
      ! test is taken on it, which is the test this loop applied on every
      ! pass before item COST7.
      !
      ! WHAT THE ACCEPTED STATE IS THEN WORTH IS MEASURED and not
      ! modeled: the probe at csm_err_probe walks the same sequence on to
      ! an increment of 1e-12 and records the distance from there to the
      ! state the test accepted.  The worst distance and what it costs are
      ! at csm_err_safety and in docs/coupled_source_loop_cost.md,
      ! section 8.
      !
      ! THE ANCHOR IS THE TIGHTEST TEST THE ADOPTED STATE HAS TO PASS, and
      ! that is the certification's composition rows: cert_tol_carrier_at(r)
      ! and cert_tol_element_at(r), 1.0e-5 in the wind and reported below it
      ! (certification.f90, decision 22), measured on the composition this
      ! step returns.  A pair that has stopped moving only
      ! to 1e-7 carries a composition error of that size into rows whose
      ! tolerance is 1e-8, so no state could certify however long it
      ! marched.  1e-8 is therefore the loosest value the acceptance
      ! interface admits, and it is what is set here.
      !
      ! WHAT IT COSTS AND WHAT LOOSENING IT WOULD BUY was MEASURED
      ! 2026-09-06 (docs/coupled_source_loop_cost.md, section 1): a decade
      ! of tolerance is 1.3 passes and about 15 percent of the wall time,
      ! so the tolerance is NOT the cost lever, and the step's own truncation
      ! error (5.3e-3 relative) is five decades above the loosest value
      ! measured.  What the cost is made of is stated at csm_max_pass.
      real*8,  parameter :: csm_T_tol      = 1.0d-8
      real*8,  parameter :: csm_comp_tol   = 1.0d-8
      ! Pass budget.  A pass costs one composition sweep and one temperature
      ! solve, and the sweep is the expensive one, so this is the cost knob
      ! of the whole step.  Exhaustion is a status, never a clamp.
      !
      ! WHAT THE COUNT IS.  MEASURED across items COST2 to COST5
      ! (docs/coupled_source_loop_cost.md, sections 2 to 4), in order of
      ! what was learned: the attenuated field is the largest single term
      ! of the sweep and is built from the composition the pass was HANDED,
      ! so the front walks inward one cell a pass; carrying the solved
      ! column forward removes the walk (10.10 to 7.13 passes) but costs
      ! the field's thread width and does not pay at sixteen threads
      ! (xuv_field_block_cells); the count is a LADDER of nested modes and
      ! equals the SLOWEST of them, not their sum (10.09, then 7.16 with the
      ! field's composition dependence removed, 4.94 with the temperature
      ! removed too, 2.00 with the He recombination rates as well), so
      ! removing a subordinate coupling changes nothing; and it is
      ! arithmetic,
      !   passes = 2 + log10(the step's own excursion
      !                      / the increment the stopping test admits)
      !                / (decades per pass) ,
      ! every term measured (0.77 decades a pass on mol_base_handoff, 7.00
      ! predicted against 6.93; 9.2 against 10.09 on wasp_full), tested
      ! against a tolerance scan (1.34 passes a decade measured, 1.30
      ! predicted) and shown not to be the cell solve (six decades of hybrd1
      ! tolerance, no change).  The increment the test admits is the
      ! tolerance divided by csm_err_safety*|theta/(1 - theta)|, since the
      ! test is taken on the error and not on the increment, so on a pass
      ! that carries an estimate the loop walks that many FEWER decades
      ! than the tolerance alone would ask: MEASURED on mol_base_handoff
      ! the factor is 0.18 on the temperature block and 0.21 on the
      ! composition block, i.e. 0.74 and 0.69 decades.  What it is worth
      ! in passes is far LESS than those decades over that rate, because
      ! the estimate has to exist at the pass that exits and a step's pass
      ! count is an integer: MEASURED 0.14 passes of 6.96 on
      ! mol_base_handoff, where 48 of 300 exits carry an estimate, and
      ! 0.08 of 10.09 on wasp_full, where 271 of 300 do (item COST7,
      ! docs/coupled_source_loop_cost.md section 8).  So the error test is
      ! not a cost lever either; what it changes is what the two
      ! tolerances bound.  The molecular path
      ! carries two modes of
      ! nearly equal cost, the composition and the temperature; the one lag
      ! that limits the composition mode's RATE there is n_tot, the third
      ! body of the three-body reactions, built from the entry composition
      ! (worth 0.8 passes at most).
      !
      ! The tolerances above are unchanged and so is the fixed point: at the
      ! fixed point a cell's entry and returned composition are one state,
      ! so every ordering builds the same columns and the same field.
      !
      ! WHAT FOLLOWS FROM THAT LAW is that a sequence whose ratio is
      ! already measured need not be walked term by term: three of its
      ! terms give both the error it still carries and its limit.  The
      ! error is what the stopping test above is taken on; the limit as
      ! the next iterate is csm_extrap_on below, default off, with its own
      ! guards.  Both read the one estimate,
      ! coupled_pair_geometric_error_estimate.
      integer, parameter :: csm_max_pass   = 20
      real*8, dimension(:), allocatable :: u_th_old_csm, u_form_csm,       &
                                           u_form_old_csm, du_form_csm,    &
                                           T_csm_prev, rho_rec_csm
      real*8, dimension(:,:), allocatable :: f_sp_csm
      ! The last pass's movement of each cell separately: the relative move
      ! of its temperature and the largest absolute move of any of its
      ! species fractions. The exit test is the statement that no cell moved
      ! (maxval of these against the two tolerances), and the same two
      ! numbers say WHICH cells are still moving, which is what the pass
      ! count is made of.
      real*8, dimension(:), allocatable :: csm_dT_cell, csm_dc_cell
      ! Which cells were at rest after the PREVIOUS pass. A cell-local skip
      ! of the next pass would be exactly the statement that these stay at
      ! rest, so the number of them that move again -- the wake-ups below --
      ! is what says whether such a skip returns the same fixed point.
      logical, dimension(:), allocatable :: csm_at_rest
      integer :: i_csm, csm_passes, csm_energy_status
      integer :: csm_n_moving, csm_n_woke
      logical :: csm_ok
      real*8  :: csm_dT, csm_dcomp, csm_mass_resid, csm_uform_frac
      real*8  :: csm_uform_frac_run = 0.0d0
      logical :: csm_debug = .false.
      character(len=32) :: csm_dbg_env
      ! Run-level records of the coupled step. All three are attempt
      ! statistics (b1 T7.2 class 2): they count what the solver did, not
      ! what the adopted state is, and a rejected trial keeps its entry.
      integer :: n_csm_passes_total = 0, n_csm_passes_worst = 0
      integer :: n_csm_worst_step = 0, n_csm_iter_cap = 0
      ! How many times the coupled loop was ENTERED. It is not the number of
      ! adopted steps: a step that is probed by the step-doubling estimate
      ! enters it for the full step and for each half, and a refused attempt
      ! keeps its entry. The mean below is per entry, so it is a statistic
      ! of the same population as the worst, and cannot exceed it.
      integer :: n_csm_calls = 0
      real*8  :: csm_mass_resid_run = 0.0d0
      integer :: csm_mass_resid_step = 0

      !------------------------------------------------!
      ! THE GEOMETRIC DECAY OF THE PASS SEQUENCE, AND WHAT IT IS READ FOR.
      !
      ! The pass count stated at csm_max_pass is arithmetic because the
      ! loop's error decays geometrically with a ratio stable from the
      ! second pass to the last.  A sequence of ratio theta whose last pass
      ! moved by d_k stands
      !     x_inf - x_k = sum_{m>=1} theta^m d_k = theta/(1 - theta) d_k
      ! from its own limit, so THREE OF ITS TERMS GIVE BOTH the error the
      ! iterate still carries and the limit itself.  Two things read that,
      ! and they read it from ONE routine,
      ! coupled_pair_geometric_error_estimate:
      !   * the stopping test above, which is a statement about the ERROR
      !     of the state it accepts and compares the estimate against
      !     csm_T_tol and csm_comp_tol;
      !   * the extrapolation below, which moves the pair to the LIMIT and
      !     is default off.
      !
      ! THE OBJECT IS THE WHOLE STATE'S DOMINANT MODE, NOT A CELL'S RATIO.
      ! A ratio taken component by component was MEASURED WORSE (12.00
      ! against 11.05 passes, items B3c and COST2) because the contraction
      ! ratio spreads 200 percent from cell to cell while the front walks
      ! through the grid and theta/(1 - theta) is unbounded as theta goes to
      ! one.  What is stable is the Rayleigh quotient of the whole increment
      ! vector, theta = <d_k, d_{k-1}> / <d_{k-1}, d_{k-1}>, the projection
      ! of the error onto the dominant decay mode (0.07 to 0.22 in size,
      ! stable from the second pass; docs/coupled_source_loop_cost.md,
      ! section 5).  One such scalar is taken for the temperature block and
      ! one for the composition block, in the norms the stopping test itself
      ! uses.
      !
      ! THE RATIO'S NORM IS NOT THE STOPPING TEST'S NORM, and that is the
      ! one approximation in the estimate: theta is an L2 Rayleigh quotient
      ! of the increment vectors, while the two tolerances are max-norm
      ! tests that one cell and one species gate.  The estimated error is
      ! therefore a MODEL of the remaining error and is good to a factor of
      ! order one only where the modeled remainder is already small --
      ! MEASURED at csm_x_reach, a candidate placed at the modeled limit
      ! leaves twice the error the model predicts once the model is asked
      ! to reach further than a factor of three, and MEASURED with the
      ! probe the max-norm distance of a state accepted near the
      ! tolerance is up to 2.4 times the projected estimate.  What the
      ! accepted state is worth is therefore not argued from the model at
      ! all: it is MEASURED against the fixed point itself by the probe at
      ! csm_err_probe, and csm_err_safety is the factor that measurement
      ! sets.
      !
      ! csm_geom_first_pass is the pass whose increment first belongs to
      ! the geometric tail: the increment of pass 1 carries the step's
      ! operator-split excursion, so the first usable pair is passes 2 and
      ! 3 and no estimate exists before the third pass.
      integer :: csm_geom_first_pass = 3
      ! An estimate is taken only where the two increments are a decaying
      ! single mode, and these say what that means.  |theta| must be
      ! nonzero and below the bound, which caps the amplification at
      ! |theta/(1-theta)| <= 4 (and, for the negative theta this loop
      ! actually has, at 0.44); the increment must have shrunk; and the
      ! geometric model must account for the increment to within
      ! csm_geom_model_tol in relative terms,
      !     || d_k - theta d_{k-1} || <= csm_geom_model_tol || d_k || ,
      ! which for a Rayleigh-quotient theta is exactly |sin| of the angle
      ! between the two increments, so the tolerance IS an alignment
      ! tolerance: 0.5 admits an alignment of |cos| >= 0.87.  A scalar
      ! ratio taken component by component cannot be asked this question at
      ! all, which is why the earlier attempt had no way to see that its
      ! ratio was not the sequence's.
      real*8 :: csm_geom_theta_max = 0.8d0
      real*8 :: csm_geom_model_tol = 0.5d0
      ! THE ESTIMATE IS TAKEN IN A DIFFERENT NORM FROM THE TEST, AND THIS
      ! IS THE FACTOR THAT COSTS.  theta is an L2 Rayleigh quotient of the
      ! increment vectors, so |theta/(1 - theta)| times the increment is
      ! the remainder of the L2-dominant mode; the two tolerances are
      ! max-norm tests that one cell and one species gate, and that cell
      ! is not in general the one the dominant mode lives in.  MEASURED
      ! with the probe at csm_err_probe on mol_base_handoff at 300 steps,
      ! the worst distance of the accepted state to its fixed point
      ! against the factor the estimate is multiplied by, and the passes a
      ! coupled step:
      !     factor          1        2        3     increment test
      !     worst dcomp  1.8e-08  1.0e-08  7.4e-09       --
      !     over 1e-08   100/299    1/299    0/299       --
      !     passes         6.14     6.76     6.82      6.96
      ! The unfactored estimate therefore accepts states OUTSIDE
      ! csm_comp_tol on a third of the steps of that case, which is the
      ! one statement the tolerance is there to make, and the factor is
      ! the smallest of the three that leaves none: 3.  The step from 2 to
      ! 3 costs 0.06 passes, because past a factor of 2 most exits are
      ! taken on the increment anyway, the guards refusing a pair whose
      ! increments have reached round-off.  A bound exceeded on one step
      ! in 300 is not a bound, which is why 2 is not the value.  The full
      ! table is in docs/coupled_source_loop_cost.md, section 8.
      real*8 :: csm_err_safety = 3.0d0
      ! TEST HOOK, default off: report no estimate on any pass.  The
      ! stopping test then falls back to the increment on every pass, which
      ! IS the test this loop applied before item COST7, and the suite
      ! asserts that this run and a run whose first usable pass is put out
      ! of reach are the same run.
      logical :: csm_geom_refuse_all = .false.
      ! DIAGNOSTIC HOOK, default off: one line for each pass sequence the
      ! estimate is measured from and one for each candidate taken.
      logical :: csm_geom_debug = .false.
      ! The pair's increments in the stopping test's own norms, this pass
      ! and the previous one, over the cells the test is taken over (1:N);
      ! the ghost cells are iterates too, but they gate nothing and are
      ! left on the plain sequence.  csm_geom_seq_broken says that the
      ! present increment does not continue the sequence the previous one
      ! belongs to, which happens only where the extrapolation has moved
      ! the pair or undone such a move.
      real*8, dimension(:),   allocatable :: csm_geom_eT, csm_geom_eT_prev
      real*8, dimension(:,:), allocatable :: csm_geom_ec, csm_geom_ec_prev
      logical :: csm_geom_have_prev, csm_geom_seq_broken
      real*8  :: csm_geom_theta_T, csm_geom_theta_c
      ! WHAT THE STOPPING TEST IS TAKEN ON: the error the pair still
      ! carries in each of the two norms -- the estimate above where the
      ! sequence gives one, and the last increment where it does not.
      ! csm_err_ok says which of the two, and the run reports how many
      ! exits were taken on each.
      logical :: csm_err_ok
      real*8  :: csm_err_T, csm_err_c
      integer :: n_csm_exit_est = 0, n_csm_exit_incr = 0
      logical :: csm_stop
      ! The pass budget of the loop as it runs.  It is csm_max_pass except
      ! under the probe below, whose passes are diagnostic and must not
      ! spend the budget the trajectory's own passes are counted against.
      integer :: csm_pass_cap
      !
      ! WHAT THE ACCEPTED STATE IS ACTUALLY WORTH, MEASURED AND NOT
      ! MODELED.  With the probe on, the pass that met the stopping test
      ! is not the last: the state it accepted is held aside, the same
      ! sequence is run on until its increment is below csm_probe_tol --
      ! four decades below the tolerance, so the state it reaches carries
      ! an error of order 1e-13 and IS the fixed point at this precision --
      ! and the max-norm distance between the two is recorded in each of
      ! the two norms of the test.  The loop then restores the state the
      ! accepted pass STARTED from and takes that pass again, which
      ! rebuilds everything the accepted pass wrote (the pair (T, f_sp) is
      ! the whole iterate of this loop), so a probe run walks the same
      ! trajectory as a run without it and the identity is asserted by the
      ! suite.  The probe is a diagnostic and default off: it costs the
      ! probe passes and one replay pass on every step.
      logical :: csm_err_probe = .false.
      real*8,  parameter :: csm_probe_tol      = 1.0d-12
      integer, parameter :: csm_probe_pass_cap = 60
      logical :: csm_probe_active, csm_probe_replay
      integer :: csm_probe_pass_accepted
      real*8, dimension(:),   allocatable :: csm_probe_T_acc,             &
                                             csm_probe_T_entry
      real*8, dimension(:,:), allocatable :: csm_probe_f_acc,             &
                                             csm_probe_f_entry
      integer :: n_csm_probe = 0, n_csm_probe_stall = 0
      integer :: n_csm_probe_over_T = 0, n_csm_probe_over_c = 0
      real*8  :: csm_probe_worst_T = 0.0d0, csm_probe_worst_c = 0.0d0

      !------------------------------------------------!
      ! THE GEOMETRIC LIMIT OF THE PASS SEQUENCE AS THE NEXT ITERATE.
      !
      ! Taking the limit of the sequence above as the next iterate is what
      ! these constants control.  It is an acceleration of the path and NOT
      ! a change of what the loop returns: the limit is a CANDIDATE, the
      ! loop continues from it, and the state that leaves still has to
      ! satisfy the same stopping test, measured the same way, on a pass
      ! taken from it.
      !
      ! WHAT IT IS WORTH NOW THAT THE STOPPING TEST IS TAKEN ON THE ERROR.
      ! The passes it used to save were the passes the loop spent walking
      ! the last factor 1/|theta| of the increment down to the tolerance,
      ! and those are exactly the passes the error test no longer spends.
      ! MEASURED on mol_base_handoff at 300 steps with the error test in
      ! place, the count with the limit taken is the count without it,
      ! where against the increment test it was worth 6.15 against 6.96.
      ! The window the reach guard leaves is now the band between the
      ! tolerance and csm_x_reach times it, a factor of two wide, so
      ! almost no pass falls in it.  The default stays off, and the code
      ! stays because the sequence it measures is what the stopping test
      ! reads.
      logical :: csm_extrap_on    = .false.
      ! How many candidates one step may take. One, because a second is
      ! unreachable in practice and the measurement says so: an accepted
      ! candidate is met by the stopping test on the very next pass in 283
      ! of 300 steps, and a candidate that falls short of the pass it
      ! replaced blocks the step, so raising the budget to three leaves
      ! mol_base_handoff at 300 steps on exactly the same ledger and the
      ! same 6.09 passes. The budget is kept as a number rather than
      ! removed because it is what makes that statement checkable.
      integer :: csm_extrap_max   = 1
      ! TEST HOOK, default off: refuse every candidate at the guard stage,
      ! so that a run can show what a refused extrapolation costs (nothing:
      ! no pass is spent and the iterate is untouched).
      logical :: csm_extrap_refuse_all = .false.
      ! MEASUREMENT HOOK, default off: undo a candidate whenever its own
      ! next pass did not beat the pass it replaced, instead of undoing it
      ! only where that pass made no progress at all.  The two rules are
      ! measured against each other at csm_x_reach.
      logical :: csm_extrap_strict = .false.
      ! HOW FAR THE ONE-MODE MODEL REACHES, and this guard is what decides
      ! whether the extrapolation pays at all.  Removing the modeled part
      ! of the error leaves the part the model does not explain, and
      ! MEASURED on mol_base_handoff that remainder is not small early in a
      ! step: taken at the third pass the extrapolation leaves a
      ! composition increment of 0.12 of the pass's own, against the 0.055
      ! the plain pass would have left, because early in the sequence
      ! several modes are alive (the ionization front is walking through
      ! cells) and the mode the Rayleigh quotient measures is the one that
      ! dominates the L2 norm, not the one that gates the max-norm stopping
      ! test.  Late in the sequence one mode is left and the model is exact
      ! to the same measurement.
      !
      ! What separates the two ends without tuning a pass number is the
      ! size of the modeled remainder ITSELF: the candidate is taken only
      ! where the estimated error is already within csm_x_reach times the
      ! stopping tolerance, which is the stopping test itself at a
      ! tolerance loosened by that factor.  The extrapolated state is then
      ! inside the ball the loop is aiming for even if the model is wrong
      ! by a factor of order one.
      !
      ! WHAT THE FACTOR IS WORTH was MEASURED on mol_base_handoff at 300
      ! steps against a control of 6.96 passes, with the stopping test
      ! taken on the increment as it then was: 6.61, 6.15, 6.36 at reach 1,
      ! 2, 3 and 7.7 to 7.8 at 10 and above (the table is in
      ! docs/coupled_source_loop_cost.md, section 5).  Above 3 the model is
      ! asked for more than it has and the loop pays a pass for each
      ! candidate its own next pass then undoes.
      real*8 :: csm_x_reach = 2.0d0
      ! The candidate is made ADMISSIBLE by construction and not by test:
      ! every constraint on it (each species fraction in [0,1], the He 2^3S
      ! sub-population no larger than the He I population that contains it,
      ! the temperature inside the bracket the energy solve uses) is linear
      ! in the step, so the largest admissible multiple of the step is a
      ! ratio test, as a limited advection step takes one.  csm_x_backoff
      ! keeps the candidate strictly inside; a candidate that has to be
      ! damped below csm_x_damp_min is not the mode's step and is dropped.
      real*8, parameter :: csm_x_backoff  = 0.95d0
      real*8, parameter :: csm_x_damp_min = 0.1d0
      ! The unextrapolated iterate, held so that a candidate whose pass
      ! turns out worse than the pass it replaced can be undone exactly.
      ! The pair (T, f_sp) IS the whole iterate of this loop: the energy
      ! solve rewrites u(3,:) and W(3,:) in full from its own temperature
      ! and from the loop-entry anchor u_th_old_csm, so nothing else is
      ! carried from one pass to the next.
      real*8, dimension(:),   allocatable :: csm_x_T_save
      real*8, dimension(:,:), allocatable :: csm_x_f_save
      logical :: csm_x_pending, csm_x_blocked
      integer :: csm_x_used
      real*8  :: csm_x_pred_dT, csm_x_pred_dc
      real*8  :: csm_x_at_dT, csm_x_at_dc
      real*8  :: csm_x_damp
      ! Run-level ledger of the extrapolation.  Attempt statistics of the
      ! same class as the pass counters above: they count what the solver
      ! did, and a refused attempt keeps its entry.  An attempt is a pass
      ! that had an estimate and an unspent budget; the guard count is what
      ! the reach guard and the ratio test then refused.  The two "worst"
      ! entries are the admissibility assertion -- both are identically
      ! zero if the ratio test is correct, and a nonzero value is the
      ! statement that an inadmissible state was handed to the next pass.
      integer :: n_csm_x_attempt = 0, n_csm_x_guard   = 0
      integer :: n_csm_x_applied = 0, n_csm_x_refused = 0
      integer :: n_csm_x_kept    = 0, n_csm_x_exit    = 0
      integer :: n_csm_x_slack   = 0
      real*8  :: csm_x_worst_fneg = 0.0d0, csm_x_worst_T_out = 0.0d0
      real*8  :: csm_x_damp_min_seen = 1.0d0
      character(len=32) :: csm_x_env
      
      !------------------------------------------------! 
      
      ! Read planetary parameters from file
      call input_read

      ! Open the setup report only once the input has been accepted:
      ! input_read refuses a configuration with error stop, and an open
      ! placed before it truncated whatever EXHALE_setup.out a previous run
      ! in the same directory had left, so a refused run erased the last
      ! accepted run's report (found by item B5i). input_read writes nothing
      ! to this unit, so nothing is lost by opening it here.
      open(unit = outfile, file = 'EXHALE_setup.out')

      ! The number of computational cells N is now final ("Grid cells:" in
      ! input.inp, 500 without the key). input_read has already allocated the
      ! grid-sized arrays of global_parameters; allocate the ones owned by the
      ! radiation modules and the state vectors of this program, all with the
      ! ghost-padded bounds 1-Ng:N+Ng. They start at zero, as they did when
      ! they were sized at compile time and lived in static storage.
      call lya_rt_allocate_arrays
      call excited_H_allocate_arrays
      call ioniz_eq_allocate_arrays
      ! Read before the state vectors are sized: the diagnostic allocates
      ! its own snapshots inside allocate_state_vectors.
      call get_environment_variable('EXHALE_UPDATE_MAP', diag_env)
      if (len_trim(diag_env) .gt. 0) call parse_update_map_request(diag_env)

      call allocate_state_vectors

      ! Build the H2, H2O and CO infrared emission/absorption tables for this
      ! run's radiating temperature. Once, here: the cell sweep that reads them
      ! runs under OpenMP, so they must be finished before it starts.
      if (mol_ir_bands) call molecular_infrared_init(T0)

      ! Build the H + H <-> H2 equilibrium-constant table that R12 (thermal
      ! dissociation of H2, obtained by detailed balance of the R15
      ! recombination) reads. Here, and unconditionally: it depends on
      ! nothing this run supplies, and the cell sweep that reads it runs
      ! under OpenMP, so it must be finished before the first parallel
      ! region opens.
      call h2_thermochemistry_init

      ! Optional parse-dump mode (env EXHALE_PARSE_DUMP=1): write every variable
      ! input_read derived from input.inp (and any base.inp override) to
      ! parse_dump.txt and stop cleanly, BEFORE any grid/IC/init or output work.
      ! Feeds the parser-refactor regression corpus
      ! (backup/regression/run_parse_corpus.sh). Env unset => run unchanged.
      call get_environment_variable('EXHALE_PARSE_DUMP', diag_env)
      if (trim(diag_env) .eq. '1') then
         call write_parse_dump
         write(*,*) '(EXHALE_main) EXHALE_PARSE_DUMP=1: parse dump written, stopping.'
         stop
      endif

      ! analytic lower column (opt-in "Lower column: <R_1bar in R_J>"):
      ! integrate the isothermal-Teq hypsometric column (Koskinen+2022) from
      ! the 1-bar radius to the 1-ubar base and report the derived base radius
      ! (chem-equilibrium mu + fully-atomic bracket), the base H2/H/He mix and
      ! mu, next to the input "Planet radius".  Report only; no override.
      if (lower_col_r1bar .gt. 0.0d0) then
         block
            real*8 :: lc_rb, lc_q2, lc_qh, lc_qhe, lc_mu, lc_rba
            call lower_column_solve(Mp, lower_col_r1bar*RJ, T0, HeH,      &
                                    1.0d0, 1.0d-6,                        &
                                    lc_rb, lc_q2, lc_qh, lc_qhe, lc_mu,   &
                                    lc_rba)
            write(*,'(a)') ' (lower_column) analytic lower colum'//&
                           'n (Koskinen+2022, isothermal Teq):'
            write(*,'(a,f8.4,a,f8.4,a)') '   r(1 ubar) = ', lc_rb/RJ,     &
               ' R_J (chem. eq.)  /  ', lc_rba/RJ, ' R_J (atomic bracket)'
            write(*,'(a,f8.4,a)') '   input "Planet radius"     = ',      &
               R0/RJ, ' R_J  (should lie in the bracket above)'
            write(*,'(a,f6.3,a,f6.3,a,f6.3,a,f6.3)')                      &
               '   base q_H2 =', lc_q2, '   q_H =', lc_qh,                &
               '   q_He =', lc_qhe, '   mu =', lc_mu
            if (lc_q2 .gt. 0.3d0) write(*,'(a)') '   NOTE: chem.-equil'// &
               'ibrium base is strongly molecular; if the planet is '//   &
               'genuinely cool, molecular physics is required '//  &
               '(for Teq~1000-2000 K hot Jupiters photochemistry '//      &
               'dissociates H2 -- see docs/lower_atmosphere_coupling).'
         end block
      endif

      !------------------------------------------------!
      
      ! Initialize simulations
      call init(W,u,f_sp)

      ! Optional IC-dump hook (env EXHALE_DUMP_IC=1): write the state exactly
      ! as initialized/loaded and stop. Lets the restart round-trip test
      ! inspect what load_IC restored BEFORE the first ionization-
      ! equilibrium solve re-equilibrates the species fractions.
      call get_environment_variable('EXHALE_DUMP_IC', diag_env)
      if (trim(diag_env) .eq. '1') then
         ! The coupling state the file carries is part of what the restart
         ! restored, and the rule that adopts it (the sec_ion_active
         ! assignment further down) has not run at this point -- so without
         ! these two lines the dump writes sec_ion=F sec_ion_step=-1 over a
         ! file that stated T, and a round-trip test reads that as the
         ! restart having lost the coupling. Same rule, stated once each.
         if (do_load_IC .and. ic_coupling_present) then
            sec_ion_active = (use_sec_ion .and. (sec_ion_immediate         &
                              .or. do_only_pp .or. ic_sec_ion_active))
            if (sec_ion_active) sec_ion_armed_step = ic_sec_ion_armed_step
         endif
         rho = W(1,:)
         v   = W(2,:)
         p   = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         ! The molecular and carrier output columns come from the arrays the
         ! equilibrium sweep fills, and no sweep has run yet: without this
         ! the dump reports a molecular restart as atomic.
         call molecular_carrier_densities_from_state(rho,f_sp)
         heat = 0.0d0; cool = 0.0d0; eta = 0.0d0
         call write_output(rho,v,p,T,heat,cool,eta,                    &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
         write(*,*) '(EXHALE_main) EXHALE_DUMP_IC=1: IC state written, stopping.'
         stop
      endif

      ! Optional steady-residual diagnostic (env EXHALE_RESIDUAL=1): evaluate
      ! the finite-volume steady residual R = du/dt of the loaded state and
      ! stop. R = dF - S for mass/momentum and dF_E - S_E - (heat - cool) for
      ! energy; at a true steady state R = 0, INDEPENDENT of how the run was
      ! stopped (du/dtu/stall). Reported as the max relative residual rate
      ! max_j |R(j,k)| / max_j scale(k,j) over the whole physical column,
      ! wind and below-escape layer normed separately and combined by the
      ! larger of the two (residual_norms). Used to rank candidate
      ! "converged" states (premature-dip vs true steady).
      ! Evaluated in WENO3 (the production scheme the states converged under).
      call get_environment_variable('EXHALE_RELOAD_EQ', diag_env)
      if (trim(diag_env) .eq. '0') reload_equilibrate = .false.
      call get_environment_variable('EXHALE_RESIDUAL', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         if (sec_ion_active) sec_ion_armed_step = 0
         ! Same treatment the marching path gets: a state read from a file is
         ! put on its own equilibrium composition before it is measured, so
         ! that the residual reported for a reloaded state is the residual of
         ! that state and not of the composition the file happened to carry.
         if ((do_load_IC .or. ic_mode .eq. 4) .and. reload_equilibrate)     &
            call equilibrate_loaded_composition
         ! Apply_BC first so the residual depends only on the interior state
         ! (ghosts set by BC), matching the Newton residual F(Y) definition.
         ! THE COMPOSITION THE BOUNDARY AND THE PRESSURE MAP READ IS THIS
         ! STATE'S: f_sp has not moved since the last get_species_densities
         ! on either path that reaches here. init ends with that call, and
         ! equilibrate_loaded_composition, the only thing between it and
         ! this block that touches f_sp, ends with it too. Where a state is
         ! INSTALLED without such a call, the composition has to be
         ! refreshed before U_to_W, whose energy-to-pressure map reads the
         ! caloric arrays.
         call Apply_BC(u)
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call assemble_residual(u, n_tot + ne, heat, cool, Rres)
         call residual_norms(Rres, u, resid_c)
         ! cell-by-cell residual profile (localize the momentum imbalance)
         block
           integer :: jj, uu
           open(newunit=uu, file='output/residual_profile.txt',           &
                status='replace', action='write')
           write(uu,'(A)') '# r[Rp]  n[cm-3]  v[cm/s]  T[K]  '//           &
                'R_mass  R_mom  R_energy'
           do jj = 1, N
              write(uu,'(1X,7(ES16.8,1X))') r(jj), W(1,jj)*n0, W(2,jj)*v0, &
                   T(jj)*T0, Rres(1,jj), Rres(2,jj), Rres(3,jj)
           end do
           close(uu)
           write(*,'(A)') ' (EXHALE_main) wrote output/residual_profile.txt'
         end block
         write(*,'(A)') ' (EXHALE_main) EXHALE_RESIDUAL=1 steady residual ||R||:'
         write(*,'(A)') '   component   max over cells of |R|/scale  [1/t_s]'
         do k = 1,3
            write(*,'(A,I2,4X,ES16.6)') '   k=', k, resid_c(k)
         enddo
         write(*,'(A,ES12.4)') '   ||R|| = max_k : ', maxval(resid_c)
         call write_residual_breakdown(Rres, u, heat, cool,                &
                                       'EXHALE_RESIDUAL diagnostic')
         write(*,*) '(EXHALE_main) EXHALE_RESIDUAL=1: residual reported, stopping.'
         stop
      endif

      ! Optional Newton-residual self-test (env EXHALE_NEWTON_TEST=1): verify
      ! the vector residual F(Y) used by the steady solver reproduces the
      ! diagnostic residual. (a) pack/unpack are exact inverses; (b) the
      ! the same relative residual per component, computed from the vector
      ! F(Y) by the solver's own convergence measure, equals the
      ! EXHALE_RESIDUAL values. Validates increment (ii)-2 before the
      ! Jacobian/PTC driver.
      ! ---- Is the steady residual a function of its argument alone? --------
      ! (env EXHALE_RESID_DETERMINISM=1.)  CONTRACT 1 OF THE RESIDUAL: with
      ! identical complete inputs and identical branch data, F(Y) comes back
      ! the same bits however many other states were evaluated and discarded
      ! in between.
      !
      ! F is evaluated at the loaded state, then at other states, then at the
      ! loaded state again, REUSING ONE OUTPUT WORKSPACE throughout so that
      ! the workspace carries each evaluation's composition into the next
      ! call.  The seed is named in every call and is the same every time, so
      ! a residual that is a state function returns the same bits; anything
      ! else means the workspace's leftover, or a module the evaluation
      ! wrote, was read.
      !
      ! The interface is what makes this pass: eval_residual takes the seed as
      ! one argument and writes the composition to another, and the two may not
      ! be the same array, so a caller CANNOT hand the sweep its own previous
      ! answer even by accident.  Before that separation the same test read
      ! 1406 of 1500 entries differing by 3.0e-3 of the row scale
      ! (Update_EXHALE_stage1.md section 146.3).
      !
      ! THE SYSTEM REPLAYED IS THE ONE A SOLVE WOULD CARRY.  The species-row
      ! registry is empty until a solve registers it, so without the two
      ! lines below this replay measures the three-unknown system in EVERY
      ! configuration and says nothing about the carrier and element rows a
      ! reload actually solves; the same statement, and the same two lines,
      ! are at EXHALE_JAC_TEST.  The species slots of Y are then filled from
      ! the state's own composition, exactly as the solve fills them.
      !
      ! THE STATES EVALUATED IN BETWEEN ARE THE PROBES OF A SOLVE, taken
      ! through the production routines so that the replay is about the code
      ! the solve runs:
      !   trial 1, 2  hydrodynamic trial points, the size a line search takes;
      !   trial 3     a species trial, moving every species slot, which no
      !               three-unknown trial reaches;
      !   the trust-region trial, an evaluation the caller DISCARDS, which is
      !               what a point on a dogleg ray is (state_is_discarded);
      !   the Jacobian-vector probe, jv_product, which holds the components
      !               of its direction that sit on an active bound
      !               (hold_the_active_bounds_of) and cuts its step to the
      !               species box frozen just above.
      ! A REJECTED MACROSTEP is not among them and cannot be: this driver
      ! runs before the marching loop and stops the run, while a rejected
      ! attempted step is a checkpoint restore inside that loop; replaying
      ! after one needs a hook in the loop, which is not this item's.
      call get_environment_variable('EXHALE_RESID_DETERMINISM', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion
         if (sec_ion_active) sec_ion_armed_step = 0
         call set_transported_species_rows(carrier_in_newton)
         call read_species_unknown_space_controls
         allocate(Yvec(neq_newton()), Fa(neq_newton()), Fb(neq_newton()),   &
                  Fc(neq_newton()), Ytry1(neq_newton()), Ytry2(neq_newton()))
         allocate(Ytry3(neq_newton()), Jvd(neq_newton()),                   &
                  vdir(neq_newton()), Ddet(neq_newton()),                   &
                  Drdet(neq_newton()))
         allocate(f_sp_w(1-Ng:N+Ng,n_species))
         call pack_U(u, Yvec)
         call pack_species_rows(u, f_sp, Yvec)
         call freeze_species_unknown_box(Yvec)
         call cell_state_scales(Yvec, Ddet)
         call cell_row_scales(Yvec, Ddet, Drdet)
         nvj_det = nvar_jac
         ! Trial states of the size a line search takes: two hydrodynamic,
         ! and one on every species slot, which the hydrodynamic pair does
         ! not reach.
         Ytry1 = Yvec;  Ytry2 = Yvec;  Ytry3 = Yvec
         do j = 1, N
            Ytry1(nvj_det*(j-1)+3) = Ytry1(nvj_det*(j-1)+3)*(1.0d0+1.0d-6)
            Ytry2(nvj_det*(j-1)+2) = Ytry2(nvj_det*(j-1)+2)*(1.0d0+1.0d-4)
            do k = 1, n_species_rows()
               Ytry3(nvj_det*(j-1)+3+k) =                                  &
                    Ytry3(nvj_det*(j-1)+3+k)*(1.0d0+1.0d-5)
            enddo
         enddo
         ! CONTROL: the same state twice in a row, nothing in between. Any
         ! difference here is state the evaluation leaves behind for itself,
         ! not something a trial injected.
         call newton_residual(Yvec, f_sp, f_sp_w, Fa)
         call newton_residual(Yvec, f_sp, f_sp_w, Fb)
         call replay_distance(Fa, Fb, Drdet, nbitdiff_ctrl, det_cl_ctrl,    &
                              det_jc_ctrl, det_kc_ctrl)
         write(*,'(A,I0,A,I0,A,I0,A,I0)')                                   &
              ' (resid_determinism) system replayed: ', nvj_det,            &
              ' unknowns per cell over ', N, ' cells = ', neq_newton(),     &
              ' entries, of which species rows ', n_species_rows()
         write(*,'(A,I0,A,I0,A,2(ES11.3,I5,I3))')                           &
              ' (resid_determinism) control, F(Y0) twice with nothing'//    &
              ' between: ', nbitdiff_ctrl, ' of ', neq_newton(),            &
              ' entries differ; row-scale max (hydrodynamic, species)'//    &
              ' with cell and slot ',                                       &
              det_cl_ctrl(1), det_jc_ctrl(1), det_kc_ctrl(1),               &
              det_cl_ctrl(2), det_jc_ctrl(2), det_kc_ctrl(2)
         ! One workspace, reused; one seed, named every time.
         call newton_residual(Yvec,  f_sp, f_sp_w, Fa)
         call newton_residual(Ytry1, f_sp, f_sp_w, Fc)
         call newton_residual(Ytry2, f_sp, f_sp_w, Fc)
         call newton_residual(Ytry3, f_sp, f_sp_w, Fc)
         call newton_residual(Yvec,  f_sp, f_sp_w, Fb)
         ! A TRUST-REGION TRIAL: an evaluation whose state the caller throws
         ! away, which is what a point on a dogleg ray is. It puts back the
         ! products it overwrote; the replay below is what says the residual
         ! does not depend on it having happened.
         call eval_residual(Ytry3, f_sp, f_sp_w, Fc, heat0, cool0,          &
                            state_is_discarded=.true.)
         ! A JACOBIAN-VECTOR PROBE with the active bounds held: the direction
         ! is deterministic and scaled by the unknowns' own magnitudes, as
         ! the solve's own probe direction is, so the step is a step in the
         ! state and not in the units of the unknowns.
         do j = 1, neq_newton()
            vdir(j) = sin(0.1d0*dble(j))
         enddo
         vdir = vdir/maxval(abs(vdir))*Ddet
         f_sp_test = f_sp
         call jv_product(Yvec, Fb, f_sp_test, vdir, Jvd, jv_ok_det)
         call newton_residual(Yvec,  f_sp, f_sp_w, Fc)
         call replay_distance(Fa, Fb, Drdet, nbitdiff_b, det_cl_b,          &
                              det_jc_b, det_kc_b)
         call replay_distance(Fa, Fc, Drdet, nbitdiff_c, det_cl_c,          &
                              det_jc_c, det_kc_c)
         write(*,'(A)') ' (resid_determinism) F(Y0) re-evaluated, one'//    &
              ' workspace reused throughout:'
         write(*,'(A,I0,A,I0,A)') '   after three trials: ', nbitdiff_b,    &
              ' of ', neq_newton(), ' entries differ'
         write(*,'(A,I0,A,I0,A,L1)')                                        &
              '   after a discarded trust-region trial and a'//             &
              ' Jacobian-vector probe: ', nbitdiff_c, ' of ',               &
              neq_newton(), ' entries differ; the probe sampled: ',         &
              jv_ok_det
         do icl = 1, 2
            det_clname = 'hydrodynamic'
            if (icl .eq. 2) det_clname = 'species'
            write(*,'(A,A,A,3(ES11.3,I5,I3))') '   ', trim(det_clname),     &
                 ' rows, row-scale distance (control, after trials,'//      &
                 ' after probes) with cell and slot ',                      &
                 det_cl_ctrl(icl), det_jc_ctrl(icl), det_kc_ctrl(icl),      &
                 det_cl_b(icl), det_jc_b(icl), det_kc_b(icl),               &
                 det_cl_c(icl), det_jc_c(icl), det_kc_c(icl)
         enddo
         if (nbitdiff_b .gt. 0 .or. nbitdiff_c .gt. 0) then
            write(*,'(A)') '   cells and rows that moved (first 12):'
            im = 0
            do j = 1, N
               do k = 1, nvj_det
                  if (Fa(nvj_det*(j-1)+k) .ne. Fc(nvj_det*(j-1)+k)) then
                     im = im + 1
                     if (im .le. 12) write(*,'(A,I5,A,F9.5,A,I2,A,ES11.3)') &
                          '     j=', j, ' r=', r(j), ' slot=', k,           &
                          '  |dF|/s = ',                                    &
                          abs(Fa(nvj_det*(j-1)+k)-Fc(nvj_det*(j-1)+k))      &
                          /max(Drdet(nvj_det*(j-1)+k), 1.0d-300)
                  endif
               enddo
            enddo
         endif
         ! THE INVARIANT IS AGAINST THE CONTROL, NOT AGAINST ZERO.  A repeat
         ! evaluation of one state is not bit-exact in every run: the first
         ! evaluation through a run builds tables, and the control above
         ! measures exactly what that costs (cells 1-2, row-scale 1e-12).
         ! What this test exists to catch is a residual that moves because
         ! ANOTHER STATE was evaluated in between, so the trials are judged
         ! against the control and against the acceptance tolerance the
         ! solver actually uses -- not against a bitwise ideal the machine
         ! does not offer.  A floor keeps the comparison meaningful when the
         ! control happens to come out clean.
         det_rs = max(det_cl_ctrl(1), det_cl_ctrl(2), 1.0d-9)
         if (max(det_cl_b(1), det_cl_b(2)) .le. det_rs .and.                &
             max(det_cl_c(1), det_cl_c(2)) .le. det_rs) then
            if (nbitdiff_b .eq. 0 .and. nbitdiff_c .eq. 0) then
               write(*,'(A)') '   VERDICT: the residual is a state'//        &
                    ' function (bitwise identical).'
            else
               write(*,'(A,ES11.3,A)') '   VERDICT: the residual is a'//     &
                    ' state function (within the control, ', det_rs,         &
                    ' in row-scale units).'
            endif
         else
            write(*,'(A)') '   VERDICT: the residual DEPENDS ON HISTORY.'
         endif
         call set_transported_species_rows(.false.)
         write(*,*) '(EXHALE_main) EXHALE_RESID_DETERMINISM=1: stopping.'
         stop
      endif

      call get_environment_variable('EXHALE_NEWTON_TEST', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         if (sec_ion_active) sec_ion_armed_step = 0
         allocate(Yvec(neq_newton()), Fvec(neq_newton()))
         call pack_U(u, Yvec)
         W = 0.0d0                                  ! reuse W as unpack target
         call unpack_U(Yvec, W)
         write(*,'(A,ES12.3)') ' (newton_test) max|unpack(pack(u))-u| (cells 1..N) = ', &
              maxval(abs(W(:,1:N) - u(:,1:N)))
         ! The seed is named and the result goes somewhere else, so no copy
         ! is needed to protect the run's own composition.
         call newton_residual(Yvec, f_sp, f_sp_test, Fvec)
         write(*,'(A)') ' (newton_test) |F|/scale per row, whole column'//  &
              ' (cf. EXHALE_RESIDUAL L-inf):'
         Rres = 0.0d0
         do j = 1,N
            do k = 1,3
               Rres(k,j) = Fvec(3*(j-1)+k)
            enddo
         enddo
         call residual_norms(Rres, u, resid_c)
         do k = 1,3
            write(*,'(A,I2,4X,ES16.6)') '   k=', k, resid_c(k)
         enddo
         write(*,*) '(EXHALE_main) EXHALE_NEWTON_TEST=1: done, stopping.'
         stop
      endif

      ! Optional banded-Jacobian self-test (env EXHALE_JAC_TEST=1): verify the
      ! colored-FD banded Jacobian of the frozen residual reproduces a
      ! directional finite difference, J*r ~= (F(Y+eps r)-F(Y))/eps. A wrong
      ! band layout gives an O(1) mismatch; a correct one matches to ~1e-6.
      call get_environment_variable('EXHALE_JAC_TEST', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         if (sec_ion_active) sec_ion_armed_step = 0
         ! THE SYSTEM THIS TEST MEASURES IS THE ONE A SOLVE WOULD CARRY.
         ! The registry is empty until a solve registers it, so without
         ! these two lines the self-test measures the three-unknown system
         ! in every configuration and says nothing about the band the
         ! species rows widen it to.
         call set_transported_species_rows(carrier_in_newton)
         k = neq_newton()
         allocate(Yvec(k), Fvec(k), F0f(k), rdir(k), Jr(k), dFD(k))
         allocate(abjac(2*kl_jac+ku_jac+1, k), Dscale(k))
         ! pack_U fills the three hydrodynamic slots of each cell; the
         ! species slots are the registry's and are filled from the state's
         ! own composition, exactly as the solve fills them.  Without this
         ! line they would be read uninitialized.
         call pack_U(u, Yvec)
         call pack_species_rows(u, f_sp, Yvec)
         f_sp_test = f_sp
         ! F0 + frozen heat/cool at Y, then the banded Jacobian
         call eval_residual(Yvec, f_sp, f_sp_test, Fvec, heat0, cool0, npart0)
         call frozen_residual(Yvec, npart0, heat0, cool0, F0f)
         call build_banded_jac(Yvec, npart0, heat0, cool0, abjac)
         ! Deterministic probe direction, SCALED BY THE UNKNOWNS' OWN
         ! MAGNITUDES.  Unscaled, max|r| = 1 in code units is a
         ! perturbation many times the state itself wherever the density is
         ! small -- MEASURED on the hot Uranus, rho = 9.5e-09 at cell 449
         ! against a step of 1e-07 -- so the difference quotient there is
         ! not a derivative and the number this prints is set by the
         ! scaling of the unknowns rather than by the band layout. The
         ! solve's own check (EXHALE_SPECIES_JAC_TEST) scales its direction
         ! the same way and for the same reason.
         do j = 1, k
            rdir(j) = sin(0.1d0*dble(j))
         enddo
         rdir = rdir/maxval(abs(rdir))
         call cell_state_scales(Yvec, Dscale)
         rdir = rdir*Dscale
         call band_matvec(abjac, rdir, Jr)
         jac_eps = 1.0d-7
         Yvec = Yvec + jac_eps*rdir
         call frozen_residual(Yvec, npart0, heat0, cool0, dFD)
         dFD = (dFD - F0f)/jac_eps
         jac_err = maxval(abs(dFD - Jr))/max(maxval(abs(Jr)), 1.0d-30)
         write(*,'(A,ES12.4)') ' (jac_test) max|F_frozen(Y)| = ', maxval(abs(F0f))
         write(*,'(A,ES12.4)') ' (jac_test) ||J*r - dFD||_inf / ||J*r||_inf = ', jac_err
         write(*,'(A)')        '   (correct band layout => ~1e-6; wrong => O(1))'
         write(*,*) '(EXHALE_main) EXHALE_JAC_TEST=1: done, stopping.'
         call set_transported_species_rows(.false.)
         stop
      endif

      ! Optional steady-state PTC-Newton solve (env EXHALE_PTC=1): solve
      ! F(Y)=0 directly from the current IC, write the converged profiles,
      ! and stop. Validation of increment (ii)-4 against the marching
      ! reference (WASP-121b ~13.71).
      call get_environment_variable('EXHALE_PTC', diag_env)
      if (trim(diag_env) .eq. '1' .and. ionization_transport) then
         ! THE SECOND ENTRANCE TO THE STEADY SOLVE, refused for the reason
         ! input_read refuses "Solver: Newton" with the ionization state
         ! carried: this route solves F(Y) = 0 with the composition on its
         ! own local root, so its last iteration would put back exactly the
         ! local equilibrium that transporting H+ exists to leave. Refusing
         ! is the honest outcome; silently solving it would return a state
         ! labelled "ionization transport" that is the local answer.
         write(*,*) '(EXHALE_main) ERROR: EXHALE_PTC=1 with "Ionization'
         write(*,*) '  transport: True". The steady PTC route holds'
         write(*,*) '  the composition at its own local root and would'
         write(*,*) '  undo the transported ionization state. Use the'
         write(*,*) '  marching route for this option. Aborting.'
         error stop 1
      endif
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         if (sec_ion_active) sec_ion_armed_step = 0
         ! THE RESOLVED CONFIGURATION OF THIS RUN, written where the
         ! marching route writes it: before the run solves anything, so a
         ! run that ends inside the solve still states what it was asked to
         ! solve. The direct steady route stops at the end of this block and
         ! never reaches the marching route's call. It is written after the
         ! three lines above, which are this route's configuration and not
         ! the input file's: the report states the scheme the solve runs
         ! under and the coupling it starts from.
         call write_setup_report
         call U_to_W(u,W)
         call eval_dt(W, dt, dt_loc)            ! CFL dt = default PTC dtau0
         resid_max = resid_th
         ! Default residual target when no "Resid tol:" was given: 1e-5 on
         ! the section-133 scale (it was 1e-3 while the rows were divided by
         ! |u|; the renormalization is measured, not a formula).
         if (resid_max .le. 0.0d0) resid_max = 1.0d-5
         ! Optional dtau0 override for experimentation (EXHALE_PTC_DTAU0=<val>):
         ! a larger start probes the Newton regime directly.
         call get_environment_variable('EXHALE_PTC_DTAU0', diag_env)
         if (len_trim(diag_env) .gt. 0) read(diag_env,*) dt
         write(*,'(A,ES10.2)') ' (EXHALE_main) PTC dtau0 = ', dt
         ! Optional frozen-base experiment (EXHALE_PTC_NFIX=<n>): anchor the
         ! first n physical cells (non-smooth base-BC region) and solve for
         ! the wind on top of them.
         call get_environment_variable('EXHALE_PTC_NFIX', diag_env)
         if (len_trim(diag_env) .gt. 0) then
            read(diag_env,*) k
            allocate(Yvec(neq_newton()))
            call pack_U(u, Yvec)
            call set_base_fix(Yvec, k)
            deallocate(Yvec)
            write(*,'(A,I0,A)') ' (EXHALE_main) frozen-base: first ', k, ' cells anchored'
         endif
         call get_environment_variable('EXHALE_PTC_JFNK', diag_env)
         ! Same solve-and-diffuse routine as the marching hand-off below, so
         ! a diffused wind can be brought to a steady state on this route too
         ! (dt_loc for the relaxation step was built by the eval_dt above).
         call steady_wind_with_element_diffusion(                       &
                                    jfnk_outer_iterations_default, dt,    &
                                    trim(diag_env) .eq. '1', j)
         ! Same run-counter ledger the marching route prints, on the route
         ! that ends here. It is called BEFORE the stop below, which is what
         ! it was missing: this route used to leave without reporting a
         ! single counter.
         call write_run_counter_report
         call element_census_reservoir('output (direct steady route)',  &
                                       rho, f_sp)
         ! The molecular columns are read from the sweep arrays: refresh
         ! them from THIS state first, so one row is one state (the lower
         ! ghosts' rho moves after the sweep; see the routine).
         call molecular_carrier_densities_from_state(rho,f_sp)
         call write_output(rho,v,p,T,heat,cool,eta,                    &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
         call assert_written_state_is_the_accepted_one
         ! The direct steady route stops here, so the resolved-configuration
         ! record is written on this route too -- element_budget.py and the
         ! Phase-E flux closure both read it, and the elemental fluxes it
         ! carries only exist once a diffusion step has run.
         call write_resolved_config
         write(*,*) '(EXHALE_main) EXHALE_PTC=1: solver done, output written, stopping.'
         stop
      endif

      !------------------------------------------------!

      ! Generate report of the current setup
      call write_setup_report
      call write_resolved_config

	!---------------------------------------------------!

	! Close outfile
	close(unit = outfile)

      !------------------------------------------------!

      !---- Start computation ----!
      
      ! Get starting time
      write(*,*) '(EXHALE_main.f90) Starting time integration..'
      start = omp_get_wtime()
      call get_environment_variable('EXHALE_PROFILE', prof_env)
      if (trim(prof_env) .eq. '1') do_profile = .true.
      ec_prof_on = do_profile
      call get_environment_variable('EXHALE_MAXSTEPS', prof_env)
      if (len_trim(prof_env) .gt. 0) read(prof_env,*) max_steps

      ! excited-H H(n=2) feedback is a decoupled (lagged-explicit) source:
      ! at the top of every timestep the H(n=2) Balmer proton source + photo-
      ! electric heating are recomputed from the current state and frozen into
      ! the global arrays that ioniz_eq reads (never inside its Newton solve).
      ! The single hydro+ionization relaxation then converges hydro, ionization
      ! and the Balmer feedback together. Off => arrays stay zero, so the run is
      ! byte-identical to the result with excited-H off.
      exc_rel = 1.0d0

      mass_flux_converged   = .false.
      step_change_converged = .false.
      du_plateaued          = .false.
      du_prev      = huge(1.0d0)
      stall_count  = 0
      n_dt_halvings     = 0
      n_steps_dt_halved = 0
      ! The controller's step-size memory (advisor decision 3): a dt reduced
      ! by a rejection bounds the next step. Nothing bounds the first one,
      ! and a run that rejects nothing never enters the branch.
      dt_last_accepted        = -1.0d0
      dt_bounded_by_rejection = .false.
      call attempted_step_read_environment()
      call energy_update_read_environment()
      ! THE CONTROLLER OWNS THE DECISION on an uncovered carrier interval
      ! (see the call site in the step), so the operator's own stop is off
      ! for the whole run.
      carrier_transport_stop_on_failure = .false.

      ! Staged secondary ionization: the SvS85 coupling starts applied only if
      ! "Secondary_ionization: Immediate" was set; otherwise it is switched on
      ! after the wind first converges without it (see the stage flip below).
      !
      ! The staging is a startup device, not a statement about the physics:
      ! from a cold IC the base feedback of the coupling amplifies the startup
      ! transient into a runaway, so the wind is first relaxed without it. A
      ! "Do only PP" run relaxes nothing -- it loads a wind that was itself
      ! converged WITH the coupling, leaves the loop after one pass, and then
      ! post-processes it, and post_process_adv runs every PH_heat_HHe call with
      ! the coupling applied (it is armed unconditionally after the loop). Left
      ! staged, the one ionization_equilibrium call of such a run would use a
      ! different photoionization rate than the post-process beside it, and its
      ! Hydro_ioniz/Ion_species files would not be the equilibrium state the
      ! _adv files correct. So arm it here, exactly as the env-driven bypass
      ! modes above do for the same reason (they also skip the time loop).
      ! THE FILE SAYS WHAT PHYSICS ITS STATE WAS PRODUCED UNDER, AND THE
      ! RESTART FOLLOWS IT.
      !
      ! The SvS85 staging is a startup device for a COLD IC, where the base
      ! feedback of the coupling amplifies the startup transient into a
      ! runaway. A state read from a previous run has no startup transient to
      ! amplify, and re-staging it re-solves its ionization at a
      ! photoionization rate it was not converged under: measured on the
      ! WASP-121b steady solution of docs/p55_base_mode.md, whose writing run
      ! armed the coupling at step 2242, re-staging moved the first
      ! ionization sweep's particle count by 8.1e-3 instead of 2.6e-4 and
      ! took the layer 250 times further from the root over 3000 steps.
      !
      ! So the coupling state travels WITH the state, in the restart file's
      ! '# coupling:' header (write_coupling_state_header, utilities.f90). A
      ! file that carries the line is restarted under the coupling it was
      ! written under; a file that does not -- anything written before the
      ! line existed -- starts staged exactly as it did before, because
      ! nothing else can be inferred about it. That is the whole rule: no
      ! environment variable and no guess.
      sec_ion_active = (use_sec_ion .and. (sec_ion_immediate .or. do_only_pp   &
                        .or. (do_load_IC .and. ic_coupling_present            &
                              .and. ic_sec_ion_active)))
      if (sec_ion_active) sec_ion_armed_step = 0
      ! THE STEP THE COUPLING WAS ARMED AT IS A PROPERTY OF THE STATE, not of
      ! the run that reads it. A run that arms the coupling itself arms it
      ! now, which is step 0 of its own marching; a stationary EVALUATION
      ! takes no step at all and writes the loaded state back unchanged, so
      ! the state it writes was armed where the file says it was, and
      ! stamping 0 on it would give the state a provenance no run produced.
      if (sec_ion_active .and. do_load_IC .and. ic_coupling_present .and.     &
          restart_intent .eq. restart_intent_stationary .and.                 &
          stationary_evaluate_only) sec_ion_armed_step = ic_sec_ion_armed_step
      if (do_load_IC .and. use_sec_ion .and.                                  &
          .not. (sec_ion_immediate .or. do_only_pp)) then
         if (ic_coupling_present) then
            write(*,'(A,L1,A)') ' (EXHALE_main) restart: the file states'//   &
                 ' secondary ionization = ', ic_sec_ion_active,               &
                 '; the run starts there'
         else
            write(*,'(A)') ' (EXHALE_main) restart: the file carries no'//    &
                 ' coupling header (written before it existed); secondary'//  &
                 ' ionization starts STAGED'
         endif
         ! The other switches the header carries are properties of the RUN,
         ! not of the state's composition, and input.inp states them, so they
         ! are reported when they disagree and never overridden.
         if (len_trim(ic_rec_method) .gt. 0 .and.                             &
             trim(ic_rec_method) .ne. trim(rec_method))                       &
            write(*,'(A)') '     note: the file was written in '//            &
                 trim(ic_rec_method)//', this run starts in '//               &
                 trim(rec_method)
      endif
      sec_flip_step  = -1

      ! WHAT THE RESTART FILE SAYS ITS STATE IS, AND WHAT THIS RUN ASKS OF IT
      ! (docs/a0_run_mode_contract_20260906.md section 3, "Restart").
      !
      ! There are two different things a physical run can do with a loaded
      ! file, and they are told apart by what the file says about itself.
      !
      ! A file written by an initialization or continuation run -- one stating
      ! mode=init, and equally one written before the field existed, which
      ! states nothing -- passed through no trajectory: its steps were
      ! relaxation iterates and no elapsed time was recorded, because none
      ! exists. That is exactly the state the HANDOFF of section 3 is defined
      ! for: the run tests it (admissible, source-consistent, caches rebuilt
      ! from it) and declares the time origin t = 0. Nothing is invented; a
      ! new trajectory starts at zero from a state the equations describe.
      !
      ! A file written by a physical integration carries the time it reached,
      ! and this run continues from it. What cannot be done is to CONTINUE a
      ! clock that is not there: a header claiming mode=phys whose t_phys is
      ! missing, unreadable, negative or not finite is a trajectory with no
      ! time on it, and no number can be supplied for it. That is refused.
      if (do_load_IC) then
         if (run_mode .eq. run_mode_phys) then
            if (ic_run_mode .eq. run_mode_phys) then
               if (.not. ic_t_phys_present) then
                  write(*,*) '(EXHALE_main) ERROR: output/Hydro_ioniz_IC.txt'
                  write(*,*) '  states mode=phys, so it claims to stand on a'
                  write(*,*) '  trajectory, but its t_phys field is missing'
                  write(*,*) '  or is not a finite non-negative number. The'
                  write(*,*) '  elapsed time of that trajectory cannot be'
                  write(*,*) '  recovered and cannot be invented. Repair the'
                  write(*,*) '  header, or restart the state as a new'
                  write(*,*) '  trajectory by setting its mode field to init.'
                  write(*,*) '  Aborting.'
                  error stop 1
               endif
               t_phys = ic_t_phys
               write(*,'(A,ES14.7,A)') ' (EXHALE_main) restart: physical'//  &
                    ' integration continues from t_phys =', t_phys, ' s'
            else
               t_phys = 0.0d0
               if (ic_run_mode_present) then
                  write(*,'(A)') ' (EXHALE_main) restart: the file holds an'//&
                       ' initialization snapshot, which carries no elapsed'
               else
                  write(*,'(A)') ' (EXHALE_main) restart: the file carries'//&
                       ' no run-mode field (written before it existed), so'
                  write(*,'(A)') '   it holds an initialization snapshot,'// &
                       ' which carries no elapsed'
               endif
               write(*,'(A)') '   time. This run takes it as the handoff'//  &
                    ' into physical integration: the state is tested below'
               write(*,'(A)') '   and the physical time origin is t = 0.'
            endif
         else if (.not. ic_run_mode_present) then
            write(*,'(A)') ' (EXHALE_main) restart: the file carries no'//   &
                 ' run-mode field (written before it existed), so its'
            write(*,'(A)') '   state is read as an initialization snapshot'//&
                 ' with no elapsed time.'
         endif
      endif

      ! du triggers start armed only for a state loaded from a previous EXHALE
      ! run (see the declarations); a generated IC has to demonstrate a du above
      ! the threshold first.
      du_stop_armed   = do_load_IC
      du_newton_armed = do_load_IC

      ! Two-stage reconstruction setup. input_read sets du_th_plm > 0 only when
      ! "Reconstruction scheme: PLM+WENO3" was given with two du_th values (it
      ! sets du_th_plm = -1 for single-stage PLM/WENO3). So du_th_plm > du_th
      ! means run PLM until du < du_th_plm, then switch reconstruction to WENO3
      ! and converge at du < du_th; otherwise the run is single-stage with the
      ! reconstruction set by "Reconstruction scheme:".
      in_plm_stage = (du_th_plm .gt. du_th)
      if (in_plm_stage) then
         rec_method = 'PLM'
         ! Keep the discretization flags consistent with rec_method: Source,
         ! Num_Fluxes, RK_rhs and Apply_BC branch on use_plm/use_weno3, not on
         ! the rec_method string.
         use_plm    = .true.
         use_weno3  = .false.
         write(*,'(A,ES9.2,A,ES9.2)') ' (EXHALE_main) Two-stage reconstruction: PLM until du <', &
                                      du_th_plm, ', then WENO3 until du <', du_th
      endif

      ! Fixed-state operator difference at the hand-off (diagnostic, env-gated
      ! like the other EXHALE_* probes). 1 = write output/opdiff_profile.txt
      ! and carry on, 2 = write it and stop. See write_operator_difference.
      call get_environment_variable('EXHALE_OPDIFF', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) opdiff_mode

      ! Mass-flux level-stability gate setup (see lev_th in parameters):
      ! converged/stalled additionally requires the mean |rho v r^2| over the
      ! wind to be unchanged (rel. < lev_th) across the last N_stall steps.
      allocate(lev_hist(0:N_stall-1))
      lev_hist        = 0.0d0
      lev_count       = 0
      lev_rel         = -1.0d0                  ! sentinel: not yet measured
      is_level_stable = (lev_th .le. 0.0d0)     ! gate disabled => always pass

      ! Residual monitor init (large until first evaluation so it cannot
      ! trip a convergence stop before the residual has been measured)
      resid_c      = huge(1.0d0)
      resid_max    = huge(1.0d0)
      fspread_gate = huge(1.0d0)
      steady_gates_converged = .false.

      ! WHAT THIS RESTART CONTINUES (Restart intent, input_read;
      ! docs/restart_contract_design_20260909.md section 2). Two of the three
      ! intents continue below: a trajectory carries the clock read from the
      ! file's header, set above, and a relaxation is the marching path this
      ! loop is. The third takes no step at all and does not return.
      if (do_load_IC .and. restart_intent .eq. restart_intent_stationary)  &
         call stationary_state_of_the_loaded_restart

      ! A state read from a file starts with the composition the file could
      ! hold, not the composition that state has. Put it on the equilibrium
      ! the first hydro step would otherwise jump to (see the routine).
      if ((do_load_IC .or. ic_mode .eq. 4) .and. reload_equilibrate)       &
         call equilibrate_loaded_composition

      ! THE HANDOFF INTO PHYSICAL INTEGRATION (contract section 3). A
      ! trajectory can only start from a state the equations describe, so the
      ! state the first step is taken from is tested before any step is taken.
      if (run_mode .eq. run_mode_phys) call physical_handoff_check

      ! One production step per requested time step, and nothing else.
      if (upmap_n .gt. 0) then
         ! The number of steps is enforced by update_map_end_step, which
         ! caps the loop once it has written its last block; count does not
         ! start from the same value for a cold and a loaded start.
         du_th     = -1.0d0;  resid_th = -1.0d0
         use_newton_solver = .false.
         ! ONE DISCRETIZATION for the whole diagnostic. A PLM+WENO3 run takes
         ! its first step in PLM and switches at the end of it, so without
         ! this the three time steps would be compared across two different
         ! spatial operators and the map would look dt-independent for that
         ! reason alone (measured: |G_dt| 1.08 -> 11.05 from the switch, not
         ! from dt). WENO3 is the production scheme every steady state in the
         ! code is converged under.
         rec_method   = 'WENO3'
         use_weno3    = .true.
         use_plm      = .false.
         in_plm_stage = .false.
         du_th_plm    = -1.0d0
         write(*,'(A,I0,A)') ' (EXHALE_main) EXHALE_UPDATE_MAP: ',         &
              upmap_n, ' production step(s), state restored between them'
      endif

      ! The step clock and the rejection probe (see their declarations).
      call get_environment_variable('EXHALE_STEP_CLOCK', diag_env)
      if (trim(diag_env) .eq. '1') trace_step_clock = .true.
      call get_environment_variable('EXHALE_REJECT_STEP', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) reject_step_probe

      ! Activate base-cell startup diagnostic if requested
      call get_environment_variable('EXHALE_P54_TS', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) p54_ts_every
      call get_environment_variable('EXHALE_DIAG_BASE', diag_env)
      if (trim(diag_env) .eq. '1') then
         diag_base = .true.
         open(unit=778, file='output/base_diag.txt', status='replace')
         write(778,'(A)') '# count  j  r[Rp]  n[cm-3]  v[cm/s]  T[K]  heat  cool'
      endif

      !------ Main temporal loop ------!
      ! Stop on any of: the mass flux constant (du < du_th), the state no
      ! longer changing (dtu < dtu_th), the du plateau (stall net), BOTH
      ! steady gates met on the marched state ("Resid tol:"), a JFNK finish
      ! that returned info = 0, or the hard iteration cap. EXHALE_MAXSTEPS
      ! and the NaN scan leave from inside the body. Each condition has its
      ! own flag and its own line in the stop report below, so the reason
      ! printed is the reason that fired.
      do while( ( .not.mass_flux_converged    .and.                          &
                  .not.step_change_converged  .and.                          &
                  .not.du_plateaued           .and.                          &
                  .not.steady_gates_converged .and.                          &
                  .not.newton_finished        .and.                          &
                  count < count_max ) .or. force_start )

            !---- Time step evaluation ----!
            ! dt_loc = cell-by-cell pseudo-dt ("Time stepping: Local"), or
            ! uniformly the global dt (default; bit-identical updates).
            call eval_dt(W,dt,dt_loc)
            ! TEST HOOK, default off (EXHALE_AS_FIXED_DT, seconds): the
            ! step is prescribed instead of taken from the CFL condition,
            ! so that one interval can be covered by one step, by two, by
            ! four and by eight and the order of the complete split update
            ! read off the successive differences. The code's time unit is
            ! R0/v0, and a prescribed step is uniform over the grid: local
            ! time stepping is refused in physical mode, the only mode in
            ! which an order this measures means anything.
            if (as_fixed_dt_seconds .gt. 0.0d0) then
               dt     = as_fixed_dt_seconds*v0/R0
               dt_loc = dt
            endif
            if (upmap_n .gt. 0) then
               ! WHICH BLOCK THIS STEP IS is counted by the diagnostic
               ! itself, not read off the loop counter: count starts at 0 on
               ! a cold start and at 1 on a loaded one, so the loop counter
               ! indexed upmap_f(0) -- out of bounds -- on every cold start,
               ! and the step then ran at dt = 0.
               upmap_scale = upmap_f(min(upmap_done + 1, upmap_n))
               dt     = dt*upmap_scale
               dt_loc = dt_loc*upmap_scale
            endif
            ! ADVISOR DECISION 3: a dt REDUCED BY A REJECTION bounds the
            ! next step, dt_next = min(dt_CFL, 2 dt_accepted). The bound is
            ! applied only after a step that was actually refused: the CFL
            ! step of a marching run can rise by more than a factor two
            ! between two neighboring states in a startup transient, and
            ! bounding that would change a run that rejected nothing, which
            ! this controller must not do.
            if (dt_bounded_by_rejection) then
               dt_cfl_now = dt
               dt = attempted_step_bound_next_dt(dt_cfl_now,              &
                                                 dt_last_accepted)
               if (dt .ne. dt_cfl_now .and. dt_cfl_now .gt. 0.0d0)        &
                  dt_loc = dt_loc*(dt/dt_cfl_now)
               dt_bounded_by_rejection = .false.
            endif
            if (do_profile) tp_step0 = omp_get_wtime()

            !-------------------------------------------------!

            ! ============ THE ATTEMPTED STEP (B3a) ============
            ! docs/b3a_attempted_step_controller_design_20260906.md. The
            ! fourteen operations of docs/b1_target_system_20260906.md
            ! section 7.1 are one trial: rows 1 to 12 change the physical
            ! state, and the adoption boundary is immediately before
            ! update_map_end_step, the first point at which that state is
            ! complete. A refused trial is restored from the checkpoint,
            ! retaken at half dt, and leaves no contribution to any
            ! physical accumulation.
            !
            ! THE STEP-DOUBLING ESTIMATE (advisor decision 4) wraps the
            ! whole thing: every n_err_every accepted steps of a physical
            ! run the step is taken once at dt and then, from the same
            ! checkpoint, as two of dt/2, and the difference is the
            ! integration error of the WHOLE operator split. Off in
            ! initialization mode, where a run makes no statement about a
            ! trajectory and so has no integration error to bound.
            !
            ! THE UNIT OF ACCEPTANCE OF A SAMPLED STEP IS THE WHOLE
            ! MACROSTEP, not one pass of it. The controller owns the
            ! macro-interval H = dt_macro: the full pass covers H and the
            ! two half passes cover 0.5*H each, so all three end at the
            ! same t + H and their difference is a temporal discretization
            ! error at a common time. Nothing inside a pass may change the
            ! interval that pass covers. A refusal there rejects the
            ! TRANSACTION, restores the state the macrostep began at,
            ! shortens H and recomputes both comparison trajectories. The
            ! estimate itself is a rejection criterion: a scaled error
            ! above one refuses the macrostep with as_reject_int_error, so
            ! the accepted state is the one whose estimate passed.
            ! NESTED, because Fortran does not short-circuit .and.: with
            ! the duty cycle off (n_err_every = 0, the default outside phys
            ! mode) a single chained test evaluates mod(count, 0), an integer
            ! division by zero that faults on every build. gfortran elides it
            ! at -O1 and above and does not at -O0, which is why the -O0
            ! build died before the first step. The nested form computes the
            ! same value the short-circuited chain did.
            n_err_passes = 1
            if (n_err_every .gt. 0 .and. run_mode .eq. run_mode_phys .and.&
                count .gt. 0) then
               if (mod(count, n_err_every) .eq. 0) n_err_passes = 3
            endif
            dt_macro    = dt
            dtloc_macro = dt_loc
            if (n_err_passes .eq. 3) then
               t_err0 = omp_get_wtime()
               call attempted_step_checkpoint_take(as_entry, u, f_sp,     &
                                                   heat, cool, eta,       &
                                                   n_dt_halve,            &
                                                   n_dt_halvings,         &
                                                   n_steps_dt_halved)
            endif
            n_macro_retry = 0

            macrostep: do

            ! THE INNER NONLINEAR ERROR OF THIS TRANSACTION'S PASSES, from
            ! zero: the estimate below is trusted only where the errors the
            ! passes left behind are well below the difference it measures,
            ! and those are the errors of THESE passes.
            call attempted_step_inner_error_reset

            err_pass: do i_err = 1, n_err_passes

            if (n_err_passes .eq. 3) then
               err_pass_index = i_err
            else
               err_pass_index = 0
            endif

            if (i_err .eq. 1) then
               dt = dt_macro;  dt_loc = dtloc_macro
               err_pass_name = 'the full step'
            else if (i_err .eq. 2) then
               ! Hold the full-step result aside and retake the step from
               ! the state it began at, as two of half the length.
               call attempted_step_checkpoint_take(as_full, u, f_sp,      &
                                                   heat, cool, eta,       &
                                                   n_dt_halve,            &
                                                   n_dt_halvings,         &
                                                   n_steps_dt_halved)
               u_err_full    = u
               f_sp_err_full = f_sp
               ! The pair operation 3 derives from the two arrays above.
               ! They are rows of the estimate in their own right: the
               ! caloric equation of state and the electron density are
               ! nonlinear in the composition, so a composition difference
               ! too small to show in any species row can still move them.
               T_err_full    = T
               ne_err_full   = ne
               call attempted_step_checkpoint_restore(as_entry, u, f_sp,  &
                                                      heat, cool, eta)
               call rebuild_state_from_checkpoint
               dt = 0.5d0*dt_macro;  dt_loc = 0.5d0*dtloc_macro
               err_pass_name = 'the first half'
            else
               ! THE SECOND HALF OF THE SAME H, set from dt_macro and not
               ! inherited from whatever the previous pass left in dt: two
               ! halves of dt_macro end at t + H by construction.
               dt = 0.5d0*dt_macro;  dt_loc = 0.5d0*dtloc_macro
               err_pass_name = 'the second half'
            endif
            ! THE INTERVAL EACH PASS COVERS, against the macro-interval it
            ! must add up to: the full pass carries H and the two halves
            ! 0.5*H of the SAME H, so the three end at one time and their
            ! difference is a temporal error there.
            if (n_err_passes .eq. 3 .and. trace_step_clock)               &
               write(*,'(A,I0,A,A,A,ES23.16,A,ES23.16)')                 &
                    '   err-pass: step=', count, ' pass=',               &
                    trim(err_pass_name), ' dt=', dt*R0/v0,               &
                    ' H=', dt_macro*R0/v0

            n_step_retry = 0
            outer_attempt: do

            n_outer_attempts = n_outer_attempts + 1
            call attempted_step_checkpoint_take(as_chk, u, f_sp, heat,    &
                                                cool, eta, n_dt_halve,    &
                                                n_dt_halvings,            &
                                                n_steps_dt_halved)
            as_inject_op  = 0
            as_trial_reason = 0
            as_trial_op     = 0
            as_carrier_ok = .true.

            !--- Thermodynamic evolution ---!

            trial: do

            ! Save previous step solution
            ! BEFORE u_old is taken: the retry loop restarts every attempt
            ! from u_old, so a restore made after this line would be undone
            ! by the first `u = u_old` inside retry_step.
            if (upmap_n .gt. 0) call update_map_begin_step

            u_old = u
            if (attempted_step_injected_refusal(as_op_checkpoint)) then
               as_inject_op = as_op_checkpoint;  exit trial
            endif

            ! Positivity of the RK3 + HLLC step.
            !
            ! The Euler equations live on the set rho > 0, rho e > 0, and the
            ! HLLC flux leaves it through sqrt(gamma p/rho). A conservative
            ! update stays on that set only under a CFL bound -- dt(|v|+c)/dr
            ! <= 1/2 for the first-order HLL family (Einfeldt et al. 1991;
            ! Batten et al. 1997), and less than that for a high-order
            ! reconstruction -- which the run's CFL number (default 0.6) does
            ! not enforce. Where the bound is violated the update can hand the
            ! next RK stage a cell whose internal energy E - rho v^2/2 has gone
            ! negative, and the run dies on the square root, not on anything
            ! physical: measured at He/H = 1 over the molecular hot-Uranus base,
            ! where a Mach 38-224 layer forms and cell 233 goes from p = +5.6 to
            ! p = -3.6e-2 in a single step at CFL 0.6 while the same arm runs
            ! 12000 steps at CFL 0.2.
            !
            ! So each stage is tested, and the violation is answered in two
            ! steps. First locally: at the offending cells only, the
            ! high-order interface fluxes are replaced by the first-order
            ! Lax-Friedrichs flux of the neighboring cell averages and those
            ! cells' update is rebuilt (positivity_limited_fluxes). That flux
            ! preserves rho > 0, rho e > 0 for dt(|v|+c)/dr <= 1 (Perthame &
            ! Shu 1996, Numer. Math. 73, 119; Zhang & Shu 2010, J. Comput.
            ! Phys. 229, 3091), which the run's CFL number does respect, and
            ! substituting it at single interfaces instead of everywhere is
            ! the flux correction of Hu, Adams & Shu (2013, J. Comput. Phys.
            ! 242, 169), used the same way in Stone et al. (2020, ApJS 249,
            ! 4). One or two cells then cost one or two first-order
            ! interfaces instead of costing the whole grid a smaller dt.
            !
            ! Second, as the backstop: a stage still outside the set -- the
            ! source terms, gravity and geometry, are not bounded by any flux
            ! choice -- is discarded and the step retaken from u_old at half
            ! dt, the standard remedy and the only one that does not put
            ! energy into the gas that the equations did not. A run that never
            ! violates the bound enters neither branch and is the run an
            ! unguarded build gives.
            n_dt_halve    = 0
            step_accepted = .false.
            retry_step: do

            ! ONE ATTEMPT, counted whether or not it is kept (contract
            ! section 5): a step re-taken at half dt after a positivity
            ! violation is one more attempt and not one more accepted step.
            n_steps_attempted = n_steps_attempted + 1

            u = u_old
            ! The composition this attempt begins from, read out of the
            ! species vector so that an attempt retaken at half dt starts
            ! from the state the checkpoint holds and not from the
            ! composition a discarded attempt advected.
            call species_advection_begin_step(f_sp)
            n_faces_positivity_at_attempt = n_faces_flux_positivity_limited

            ! FIRST RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            if (do_profile) tp_b = omp_get_wtime()
            call reconstruction_continuation_rhs(u,WL,WR,dF,S)

            do k = 1,3
               u1(k,:) = u(k,:) - dt_loc*(dF(k,:) - S(k,:))
            enddo

            ! Apply boundary conditions
            call Apply_BC(u1)

            ! Local repair first (see the note above); the test below then
            ! decides whether it was enough.
            if (.not. positive_density_and_internal_energy(u1)) then
               call positivity_limited_fluxes(1,u,u,S,dt_loc,u1,flux_corr_ok)
               if (flux_corr_ok) call Apply_BC(u1)
            endif


            if (.not. positive_density_and_internal_energy(u1)) then
               if (n_dt_halve .eq. 0) n_steps_dt_halved = n_steps_dt_halved + 1
               n_dt_halve   = n_dt_halve + 1
               n_dt_halvings = n_dt_halvings + 1
               dt     = 0.5d0*dt
               dt_loc = 0.5d0*dt_loc
               ! Give up after n_dt_halve_max bisections: the state is then
               ! inadmissible for a reason no step size repairs, and the NaN
               ! detector below reports it with the profile.
               if (n_dt_halve .le. n_dt_halve_max) cycle retry_step
               exit retry_step
            endif

            ! THE SPECIES ROWS OF THIS STAGE, on the face mass fluxes the
            ! mass row of this same stage was built from.  Placed after the
            ! admissibility test, for two reasons: the face fluxes of an
            ! inadmissible stage are not a flow anything rides on, and a
            ! face whose flux the positivity repair replaced carries the
            ! replacement in face_flux, so the species fluxes here sum to
            ! the mass flux the density actually used.
            call species_advection_stage(1,u(1,:),u(1,:),u1(1,:),          &
                                         face_flux(1,:),dt_loc)

            !----------------------------

            ! SECOND RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call reconstruction_continuation_rhs(u1,WL,WR,dF,S)

            ! Advance in time
            do k = 1,3
               u2(k,:) = (3.0*u(k,:) + u1(k,:) - dt_loc*(dF(k,:) - S(k,:)))/4.0
            enddo

            ! Apply boundary conditions
            call Apply_BC(u2)

            if (.not. positive_density_and_internal_energy(u2)) then
               call positivity_limited_fluxes(2,u,u1,S,dt_loc,u2,flux_corr_ok)
               if (flux_corr_ok) call Apply_BC(u2)
            endif

            if (.not. positive_density_and_internal_energy(u2)) then
               if (n_dt_halve .eq. 0) n_steps_dt_halved = n_steps_dt_halved + 1
               n_dt_halve   = n_dt_halve + 1
               n_dt_halvings = n_dt_halvings + 1
               dt     = 0.5d0*dt
               dt_loc = 0.5d0*dt_loc
               ! Give up after n_dt_halve_max bisections: the state is then
               ! inadmissible for a reason no step size repairs, and the NaN
               ! detector below reports it with the profile.
               if (n_dt_halve .le. n_dt_halve_max) cycle retry_step
               exit retry_step
            endif

            call species_advection_stage(2,u(1,:),u1(1,:),u2(1,:),         &
                                         face_flux(1,:),dt_loc)

            !----------------------------

            ! THIRD RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call reconstruction_continuation_rhs(u2,WL,WR,dF,S)

            ! Advance in time
            do k = 1,3
               u(k,:) = (u(k,:) + 2.0*(u2(k,:) - dt_loc*(dF(k,:) - S(k,:))))/3.0
            enddo

            ! Apply boundary conditions
            call Apply_BC(u)

            ! u_old, not u: this stage overwrites u in place, so the state at
            ! the beginning of the step is only left in u_old.
            if (.not. positive_density_and_internal_energy(u)) then
               call positivity_limited_fluxes(3,u_old,u2,S,dt_loc,u,flux_corr_ok)
               if (flux_corr_ok) call Apply_BC(u)
            endif

            if (.not. positive_density_and_internal_energy(u)) then
               if (n_dt_halve .eq. 0) n_steps_dt_halved = n_steps_dt_halved + 1
               n_dt_halve   = n_dt_halve + 1
               n_dt_halvings = n_dt_halvings + 1
               dt     = 0.5d0*dt
               dt_loc = 0.5d0*dt_loc
               ! Give up after n_dt_halve_max bisections: the state is then
               ! inadmissible for a reason no step size repairs, and the NaN
               ! detector below reports it with the profile.
               if (n_dt_halve .le. n_dt_halve_max) cycle retry_step
               exit retry_step
            endif

            call species_advection_stage(3,u_old(1,:),u2(1,:),u(1,:),      &
                                         face_flux(1,:),dt_loc)

            ! THE REJECTION PATH, ON DEMAND. A step is refused here exactly
            ! as a positivity violation refuses one: the attempt is discarded,
            ! dt is halved and the step is retaken from u_old. It exists
            ! because the rule that a rejected trial advances no clock and no
            ! ledger cannot be tested on a run that never rejects a step, and
            ! the run's own rejections depend on the state. Unset, the probe
            ! is -1 and the branch is never entered, so a run without it is
            ! the run an unguarded build gives.
            if (count .eq. reject_step_probe .and. n_dt_halve .eq. 0) then
               n_steps_dt_halved = n_steps_dt_halved + 1
               n_dt_halve    = n_dt_halve + 1
               n_dt_halvings = n_dt_halvings + 1
               dt     = 0.5d0*dt
               dt_loc = 0.5d0*dt_loc
               write(*,'(A,I0,A)') '   EXHALE_REJECT_STEP: attempt at '//  &
                    'step ', count, ' refused; retaken at half dt'
               cycle retry_step
            endif

            ! The step is accepted here: all three stages held rho > 0,
            ! rho e > 0, so the repairs made since this attempt began belong
            ! to the trajectory the run keeps.
            n_faces_flux_positivity_limited_accepted =                     &
               n_faces_flux_positivity_limited_accepted                    &
               + (n_faces_flux_positivity_limited                          &
                  - n_faces_positivity_at_attempt)

            step_accepted = .true.
            exit retry_step
            enddo retry_step

            ! ROW 2 IS DONE. Its outcome is DEMOTED: step_accepted no
            ! longer means the step is accepted, it means the three
            ! Runge-Kutta stages produced an admissible state. The outer
            ! controller decides at the adoption boundary below.
            !
            ! THE OPERATOR THE STAGES INTEGRATED, read back rather than
            ! recomputed. The first stage formed u1 = u_old - dt_loc*(dF-S)
            ! with dF and S evaluated at u_old, so dt_loc*L(u_old) is
            ! u_old - u1 and no statement has to be placed inside the
            ! Runge-Kutta block to obtain it. That matters: a statement
            ! there changes how the compiler contracts the update into an
            ! FMA, which was measured at one ulp per step growing to 1e-11
            ! over a converging run (code site records, 2026-09-03).
            do k_as = 1, 3
               do j = 1, N
                  if (dt_loc(j) .gt. 0.0d0) then
                     L_hyd(k_as,j) = (u_old(k_as,j) - u1(k_as,j))/dt_loc(j)
                  else
                     L_hyd(k_as,j) = 0.0d0
                  endif
               enddo
            enddo
            ! THE STATE ROW 2 RETURNED, kept for the time-discrete
            ! hydrodynamic row of the acceptance predicate. That row is a
            ! verdict on the Runge-Kutta update, so it is measured on the
            ! state that update produced. Measuring it on the state at the
            ! adoption boundary instead would compare the whole operator
            ! split against the hydrodynamic operator alone, which is O(1)
            ! wherever the sources dominate and is a statement about the
            ! split, not about row 2.
            u_hyd = u

            if (.not. step_accepted) then
               ! The nested positivity retry gave up after n_dt_halve_max
               ! bisections. That is an outer rejection with the operation
               ! named, not a state the step controller may adopt.
               as_inject_op = -as_op_hydro
               exit trial
            endif
            if (attempted_step_injected_refusal(as_op_hydro)) then
               as_inject_op = as_op_hydro;  exit trial
            endif

            ! Update-map snapshots sit BETWEEN operator-split stages, never
            ! between the statements of one. A statement placed inside the RK3
            ! block changes how gfortran contracts u - dt*(dF - S) into an FMA
            ! (-ffp-contract=fast is the default), which is one ulp per step and
            ! grows to 1e-11 over a converging run: measured 2026-09-03, it broke
            ! the bitwise regression on wasp_full while the diagnostic was off.
            ! Anything measured inside a stage must be measured out of line.
            if (upmap_n .gt. 0) u_umA = u

		!------------------------------------------------!
 		
 		!---- Ionization Equilibrium ----!
 		
            ! Extract primitive variables
            if (do_profile) tp_hydro = tp_hydro + (omp_get_wtime() - tp_b)
            call U_to_W(u,W)
            rho = W(1,:)
            v   = W(2,:)
            p   = W(3,:)
            
            ! Evaluate species densities, ne and n_tot (single policy point)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                       nheiii,nheiTR,nm,ne,n_tot)

            ! Temperature profile
            call comp_T_from_p(p,n_tot,ne,T)

            if (attempted_step_injected_refusal(as_op_primitives)) then
               as_inject_op = as_op_primitives;  exit trial
            endif

            ! He/H diffusive separation, in two halves.  The ADVECTIVE half
            ! rode on the face mass fluxes inside the Runge-Kutta stages
            ! above and is written into f_sp here, at the operator-split
            ! point the element transport occupies; the DIFFUSIVE half is the
            ! implicit step that follows, taken at the advected composition.
            ! The ionization equilibrium below then re-solves the stages,
            ! conserving the new element amounts.  No-op unless he_diffusion
            ! is set (byte-identical when off).
            if (species_advection_active())                           &
               call species_advection_project(f_sp)
            if (he_diffusion)                                         &
               call element_diffusion_step(rho,v,T,f_sp,dt_loc)

            if (attempted_step_injected_refusal(as_op_diffusion)) then
               as_inject_op = as_op_diffusion;  exit trial
            endif

            ! Vertical transport of the molecular carriers (H2, OH, H2O,
            ! CO) with their chemistry, in the same operator-split slot the
            ! element diffusion occupies and for the same reason: after the
            ! hydro update, before the ionization solve, which then solves
            ! the remaining stages against the transported partition rather
            ! than recomputing it. No-op unless the oxygen chemistry is on
            ! with its transport (byte-identical when off).
            !
            ! THE INTERVAL IS CALLED DIRECTLY, not through the wrapper.
            ! carrier_transport_interval reports an uncovered interval
            ! instead of stopping, and the step controller turns that into
            ! a rejection of the whole attempted step with the operation
            ! named (design section 4.4). A3's own substep retries are
            ! unchanged and stay inside the operator, which is the right
            ! granularity for them: forcing them up to the outer boundary
            ! would discard rows 1 to 4 for a failure the operator fixes at
            ! a tenth of the cost.
            if (do_profile) tp_b = omp_get_wtime()
            as_carrier_ok = .true.
            call photochemical_transport_step(rho,v,f_sp,dt_loc,          &
                                              as_carrier_status)
            if (do_profile) tp_trans = tp_trans + (omp_get_wtime() - tp_b)

            ! AN EXHAUSTED CARRIER INTERVAL, and what each mode does with
            ! it. The operator's own stop is suppressed for the whole run
            ! (carrier_transport_stop_on_failure = .false., set at startup)
            ! because THE STEP CONTROLLER OWNS THE DECISION: in phys mode an
            ! uncovered interval refuses the whole attempted step, which is
            ! restored and retaken at half dt, and the module stopping first
            ! would take that decision away. In init mode the state is a
            ! numerical iterate, so A3b's own handling stands: the interval
            ! is recorded in the attempts ledger, the march goes on from the
            ! entry carriers, and the carrier history is marked as not
            ! certifiable. The controller does not reject there.
            if (as_carrier_status .eq. carrier_interval_exhausted) then
               as_carrier_ok = .false.
               if (run_mode .eq. run_mode_phys) exit trial
               as_carrier_ok = .true.
            endif
            if (attempted_step_injected_refusal(as_op_carriers)) then
               as_inject_op = as_op_carriers;  exit trial
            endif

            ! Refresh the lagged H(n=2) Balmer source + heating from
            ! the current state before the ionization/energy solve.
            if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)

            if (attempted_step_injected_refusal(as_op_excited_H)) then
               as_inject_op = as_op_excited_H;  exit trial
            endif

            !------------------------------------------------!
            ! THE COUPLED SOURCE STEP (rows 7 to 9 of the enumerated step,
            ! docs/b1_target_system_20260906.md T1.4 to T1.9 and T2.1).
            !
            ! ONE local source step per cell at fixed volume, with the
            ! temperature and the composition as its unknowns:
            !
            !    (a) the composition equilibrium at the returned T, with a
            !        transported carrier entering as its imposed value;
            !    (b) u_th(T,c) - u_th_old + sum_s eps_s (n_s - n_s_old)
            !          = dt [ heat(c,T) - cool(c,T) ] .
            !
            ! WHAT THIS REPLACES.  The old rows 7 to 9 were: solve the
            ! composition at the cell's T, rewrite rho from the composition,
            ! rebuild the pressure at the same T from the NEW particle count
            ! and push it into the energy row (the C1 projection: a change of
            ! the conserved energy with no source), then solve the
            ! temperature at the new composition against the energy the
            ! projection had just written.  The projection is gone: no
            ! comp_p_from_T, no W_to_U, and the energy row is anchored on
            ! u_th_old, the thermal energy the cell actually had, with the
            ! composition change entering through its formation energy.  The
            ! density is gone as well: ioniz_eq takes it intent(in) (T2.1)
            ! and returns its own mass sum as a check (T2.2).
            !
            ! FIXED POINT, NOT NEWTON, AND WHY.  The two solves are
            ! whole-grid calls with different shapes: eval_cool integrates
            ! the fine-structure line transfer over the columns before its
            ! cell-local part, and the composition sweep rebuilds the
            ! radiation field of the whole column.  A Newton iteration on the
            ! coupled (T, c) system would need the derivative of the
            ! composition with respect to T, which neither solve exposes and
            ! which is not a cell-local quantity once the field is included.
            ! What IS available is that each of the two solves is exact for
            ! its own unknown at the other's value, so the composite map is a
            ! fixed-point iteration whose error contracts by
            ! (d u_th/d c)(d c/d T) / c_v per pass -- small wherever the
            ! composition is not itself thermally driven, which is where the
            ! coupling is weak, and slow exactly where the two are strongly
            ! coupled.  The pass count is bounded and its exhaustion is a
            ! reported status, never a silently accepted state.
            !
            ! THE PAIR THAT IS RETURNED.  The loop exits when the last pass
            ! moved the temperature by at most csm_T_tol in relative terms
            ! AND moved the composition by at most csm_comp_tol.  The state
            ! that leaves is then (T, c) with c the equilibrium composition
            ! of a temperature that differs from T by less than csm_T_tol,
            ! and with the energy row of (b) satisfied to the energy solve's
            ! own residual tolerance at exactly that c.  It is also the state
            ! whose heating the run writes: at the fixed point the sweep's
            ! entry composition and its returned composition agree to
            ! csm_comp_tol, so the one-sweep rate lag that made the written
            ! heat column and Heating_breakdown.txt two different numbers
            ! closes with it.
            u_th_old_csm = u(3,:) - 0.5d0*rho*v**2.0d0
            ! The thermal energy the coupled source step starts from, for
            ! the T1.5 closure of that step alone: everything the attempted
            ! step does outside this loop moves the thermal energy too, and
            ! taking the pair here and after the loop separates the two.
            call attempted_step_thermal_energy_before_sources(u)
            call formation_energy_density(rho, f_sp, u_form_csm)
            u_form_old_csm = u_form_csm
            csm_passes  = 0
            csm_dT      = 0.0d0
            csm_dcomp   = 0.0d0
            csm_ok      = .false.
            csm_energy_status = 0
            if (csm_debug) csm_at_rest = .false.
            ! The increment history belongs to one step's pass sequence
            ! and to no other, and so does the extrapolation's own state.
            csm_geom_have_prev  = .false.
            csm_geom_seq_broken = .false.
            csm_err_ok = .false.
            csm_err_T  = 0.0d0
            csm_err_c  = 0.0d0
            csm_x_pending   = .false.
            csm_x_blocked   = .false.
            csm_x_used      = 0
            csm_probe_active = .false.
            csm_probe_replay = .false.

            coupled_source: do i_csm = 1, csm_pass_cap
               csm_passes = i_csm
               f_sp_csm = f_sp

               ! (a) the composition at the current temperature.
               ! The ledger of THIS sweep, so that the marching stop can say
               ! whether the state it is about to call steady rests on a cell
               ! whose composition is not a root of the chemical network
               ! (acceptance classes 4 and 6). Without it the stop makes no
               ! statement
               ! about the chemistry and steady_gates_met rests on its other
               ! three conditions.
               if (do_profile) tp_a = omp_get_wtime()
               call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep,          &
                             rho_recon=rho_rec_csm)
               if (do_profile) tp_ion = tp_ion + (omp_get_wtime() - tp_a)

               ! Evaluate partial densities, ne and n_tot (single policy point)
               call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                          nheiii,nheiTR,nm,ne,n_tot)

               ! T2.2: the species mass sum against the density the step was
               ! given.  Under the declared convention (the electron mass is
               ! carried with its ion) an ionization moves no mass between
               ! species, so this is a statement about the element budget of
               ! the sweep and not about the ionization state.  It is
               ! measured, reported and never applied.
               csm_mass_resid = maxval(abs(rho_rec_csm(1:N) - rho(1:N))    &
                                       /max(rho(1:N), 1.0d-99))
               if (csm_mass_resid .gt. csm_mass_resid_run) then
                  csm_mass_resid_run  = csm_mass_resid
                  csm_mass_resid_step = count
               endif

               ! The formation/excitation reservoir of the composition just
               ! returned, and its change over the step.  p0 = n0 kB T0 is
               ! the code's energy-density unit.
               call formation_energy_density(rho, f_sp, u_form_csm)
               du_form_csm = (u_form_csm - u_form_old_csm)/p0
               csm_uform_frac = maxval(abs(du_form_csm(1:N)))              &
                  /max(maxval(abs(dt_loc(1:N)*(heat(1:N)-cool(1:N)))),     &
                       1.0d-99)
               ! Step 0 is excluded: there the composition jumps from the
               ! initial condition to its equilibrium, which is not a step of
               ! the trajectory and would set this record for the whole run.
               if (count .gt. 0 .and. csm_uform_frac .gt. csm_uform_frac_run)&
                  csm_uform_frac_run = csm_uform_frac

               ! (b) the temperature that closes the energy row at this
               ! composition.
               !
               ! WHY du_form IS NOT PASSED HERE, AND WHAT IT WOULD TAKE.
               ! The energy row of b1 T1.6 reads
               !   u_th(T,c) - u_th_old + sum_s eps_s dn_s = dt [Q_ext] ,
               ! and T1.6 states that Q_ext carries ONLY the exchange with
               ! the radiation field, every collisional reaction heat being
               ! in the reservoir term instead.  This code's `heat` and
               ! `cool` are not that quantity.  They are the NET THERMAL
               ! source, and each of them already carries the reservoir
               ! transfers of its own processes:
               !   - photoionization deposits h nu - E_th, not h nu, so the
               !     ionization potential the photon paid never passes
               !     through the thermal pool at all (PH_heat_HHe, and the
               !     recipient tables of h2_photo_channels);
               !   - collisional ionization takes E_th out of the electrons
               !     and it is the `coio` term of eval_cool
               !     (util_ion_eq.f90:1799, e_th_HI_erg*a_ion_HI*nhi*ne);
               !   - radiative recombination returns E_th as a photon that
               !     leaves the cell, and the electron's own kinetic energy
               !     is the `reco` term;
               !   - the collisional molecular network, the two Penning
               !     branches, the associative branch, the Lyman-Werner
               !     fragment and fluorescence channels and the oxygen
               !     channels are each deposited as a heat above, from the
               !     same species formation-energy table this reservoir is
               !     built from.
               ! Adding sum_s eps_s dn_s to that row therefore counts every
               ! collisional transfer twice and charges the gas for the
               ! ionization energy that the photons paid.  MEASURED: with
               ! the term added, wasp_full fails at step 0 with 368 cells
               ! driven to the equation-of-state floor, the base cell asking
               ! for 23.6 K in place of 2358 K, because the composition jump
               ! from the cold isothermal initial condition to its
               ! equilibrium is charged to the thermal pool; and on the
               ! relaxed state of the B3a report the term is 28 percent of
               ! the step's thermal budget (d_u_form = -2.31e-4 against
               ! d_u_th = -8.16e-4), so it is not a small effect that could
               ! be carried and reported.
               !
               ! What T1.6 asks for is reachable, and it is a redefinition
               ! of `heat` and `cool` and not an addition to this row: the
               ! photoheating would have to deposit the whole absorbed
               ! photon, `coio` would have to leave the cooling sum, every
               ! collisional reaction heat above would have to be deleted,
               ! and the recombination photon's potential energy would have
               ! to appear as a negative external exchange.  Those live in
               ! util_ion_eq.f90, Cool_coeff.f90 and
               ! electron_energy_degradation.f90.  Until that redefinition
               ! is made, the reservoir term this code owes the row is
               ! identically what `heat` and `cool` already carry, so what
               ! is passed is the anchor u_th_old alone -- which is the part
               ! of T1.6 that IS a defect of this row: the energy the step
               ! starts from used to be rebuilt at the NEW composition (the
               ! C1 projection) instead of being the energy the cell had.
               ! The reservoir itself is measured every step and reported.
               if (use_semi_implicit_energy) then
                  if (do_profile) tp_b = omp_get_wtime()
                  call solve_energy_semi_implicit(u,W,dt_loc,heat,cool,    &
                          f_sp,count,status=csm_energy_status,             &
                          T_start=T, u_th_old=u_th_old_csm)
                  if (do_profile)                                          &
                     tp_energy = tp_energy + (omp_get_wtime() - tp_b)
                  T_csm_prev = T
                  call comp_T_from_p(W(3,:),n_tot,ne,T)
                  csm_dT_cell(1:N) = abs(T(1:N) - T_csm_prev(1:N))         &
                                     /max(abs(T(1:N)), 1.0d-99)
                  csm_dT = maxval(csm_dT_cell(1:N))
               else
                  ! The explicit forward-Euler update
                  ! ("Energy solver: Explicit"), on the same energy row and
                  ! with the same anchor u_th_old, so the projection is gone
                  ! from this path too. ONE PASS, and that is what the option
                  ! means: an explicit update evaluates its source at the
                  ! state the step starts from, and iterating it to a fixed
                  ! point would make it implicit and no longer the update the
                  ! key asks for. The pair it returns is therefore not a
                  ! fixed point, which is the accuracy the option chooses.
                  u(3,:) = u_th_old_csm + 0.5d0*rho*v**2.0d0               &
                         + dt_loc*(heat - cool)
                  T_csm_prev = T
                  call U_to_W(u,W)
                  call comp_T_from_p(W(3,:),n_tot,ne,T)
                  csm_dT_cell(1:N) = abs(T(1:N) - T_csm_prev(1:N))         &
                                     /max(abs(T(1:N)), 1.0d-99)
                  csm_dT = maxval(csm_dT_cell(1:N))
                  csm_ok = .true.
                  exit coupled_source
               endif

               do j = 1, N
                  csm_dc_cell(j) = maxval(abs(f_sp(j,:) - f_sp_csm(j,:)))
               enddo
               csm_dcomp = maxval(csm_dc_cell(1:N))

               ! (c) THE ERROR THE PAIR STILL CARRIES, from the geometric
               ! decay of this step's own pass sequence.  It is measured
               ! before the stopping test because the test is a statement
               ! about that error and not about the last step towards it.
               call coupled_pair_geometric_error_estimate

               if (csm_debug) then
                  ! How many cells are still moving, and where the worst
                  ! one is: the pass count of the whole grid is set by that
                  ! cell alone, so a global maximum on its own does not say
                  ! what the iteration is waiting for.
                  ! (the intrinsic COUNT is shadowed here: `count` is the
                  ! marching step index of this program)
                  csm_n_moving = sum(merge(1, 0,                           &
                                     csm_dT_cell(1:N) .gt. csm_T_tol       &
                                .or. csm_dc_cell(1:N) .gt. csm_comp_tol))
                  ! Wake-ups: cells that were at rest after the previous
                  ! pass and moved again on this one. They are the cells a
                  ! cell-local skip would have frozen at a state that is not
                  ! the fixed point.
                  csm_n_woke = sum(merge(1, 0, csm_at_rest(1:N) .and.      &
                                    (csm_dT_cell(1:N) .gt. csm_T_tol       &
                                .or. csm_dc_cell(1:N) .gt. csm_comp_tol)))
                  csm_at_rest(1:N) = .not.                                 &
                                    (csm_dT_cell(1:N) .gt. csm_T_tol       &
                                .or. csm_dc_cell(1:N) .gt. csm_comp_tol)
                  ! errT and errC are what the stopping test is taken
                  ! on: the estimated remaining error where the sequence
                  ! gives one (est=T) and the increment itself where it
                  ! does not.  The moving count is still counted against
                  ! the increment, because it says which cells are still
                  ! moving and not which are still in error.
                  write(*,'(a,i0,a,i0,a,es10.3,a,es10.3,a,es10.3,a,'//    &
                        'es10.3,a,l1,a,i0,a,i0,a,i0,a,i0)')                &
                     ' CSMDBG step ', count, ' pass ', i_csm, ' dT= ',     &
                     csm_dT, ' dC= ', csm_dcomp, ' errT= ', csm_err_T,     &
                     ' errC= ', csm_err_c, ' est= ', csm_err_ok,           &
                     ' moving= ', csm_n_moving, ' worstT= ',               &
                     maxloc(csm_dT_cell(1:N), 1), ' worstC= ',             &
                     maxloc(csm_dc_cell(1:N), 1), ' woke= ', csm_n_woke
               endif
               ! (d) THE STOPPING TEST, on the error the pair carries in
               ! each of the two norms: the estimate where the sequence
               ! gives one, the last increment where it does not.  The two
               ! detours are the probe's: a probe pass runs the same
               ! sequence on to csm_probe_tol so that the accepted state's
               ! distance to the fixed point can be MEASURED, and the
               ! replay pass after it retakes the accepted pass from the
               ! state that pass started from, so the trajectory is the
               ! trajectory of a run without the probe.
               if (csm_probe_replay) then
                  csm_stop = .true.
               elseif (csm_probe_active) then
                  csm_stop = (csm_dT .le. csm_probe_tol .and.             &
                              csm_dcomp .le. csm_probe_tol)               &
                             .or. i_csm .ge. csm_pass_cap - 1
               else
                  csm_stop = csm_err_T .le. csm_T_tol .and.               &
                             csm_err_c .le. csm_comp_tol
               endif

               if (csm_stop) then
                  if (csm_probe_replay) then
                     ! The accepted pass again, from the pair it was taken
                     ! from, so that everything it wrote is what leaves the
                     ! loop.  The pass count reported is the accepted one:
                     ! the probe's passes are diagnostic and belong to no
                     ! trajectory.
                     csm_probe_replay = .false.
                     csm_passes       = csm_probe_pass_accepted
                     csm_ok           = .true.
                     exit coupled_source
                  elseif (csm_probe_active) then
                     if (csm_dT .gt. csm_probe_tol .or.                   &
                         csm_dcomp .gt. csm_probe_tol)                    &
                        n_csm_probe_stall = n_csm_probe_stall + 1
                     call coupled_pair_fixed_point_distance
                     csm_probe_active = .false.
                     csm_probe_replay = .true.
                     T    = csm_probe_T_entry
                     f_sp = csm_probe_f_entry
                     cycle coupled_source
                  endif
                  ! From here the pass that met the test is a pass of the
                  ! trajectory, so its exit goes in the run's ledger:
                  ! which of the two bounds it was taken on.  An exit on
                  ! the estimate is a statement about the error; an exit
                  ! on the increment is the stricter statement the loop
                  ! falls back to.  It is counted here and not at the
                  ! probe's replay exit, which retakes this same pass.
                  if (csm_err_ok) then
                     n_csm_exit_est = n_csm_exit_est + 1
                  else
                     n_csm_exit_incr = n_csm_exit_incr + 1
                  endif
                  if (csm_err_probe .and.                                 &
                      i_csm .le. csm_pass_cap - 3) then
                     csm_probe_T_acc   = T
                     csm_probe_f_acc   = f_sp
                     csm_probe_T_entry = T_csm_prev
                     csm_probe_f_entry = f_sp_csm
                     csm_probe_pass_accepted = i_csm
                     csm_probe_active  = .true.
                     cycle coupled_source
                  endif
                  csm_ok = .true.
                  ! An extrapolated state that the stopping test met on the
                  ! very next pass: the passes the plain iteration would
                  ! have spent walking down to the tolerance are the ones
                  ! this saved, and counting these is the only honest
                  ! measure of what the extrapolation bought.
                  if (csm_x_pending) n_csm_x_exit = n_csm_x_exit + 1
                  exit coupled_source
               endif

               ! THE PASS BUDGET OF THE TRAJECTORY IS csm_max_pass
               ! WHETHER OR NOT THE PROBE IS ON.  The probe's own passes
               ! run past that budget on purpose, and letting them spend it
               ! would change which state a capped step accepts, so the
               ! exhaustion is taken here and not from the DO bound.
               if (.not. csm_probe_active .and. .not. csm_probe_replay     &
                   .and. i_csm .ge. csm_max_pass) exit coupled_source

               ! (e) THE GEOMETRIC LIMIT OF THE SEQUENCE AS THE NEXT
               ! ITERATE, taken only on a pass that did NOT satisfy the
               ! test above, so that the state which leaves the loop is
               ! always the state of a pass that met the test.  A candidate
               ! taken on the previous pass has its verdict read first,
               ! from the increment the pass just measured.  A probe pass
               ! never takes one: its sequence is being walked to the end
               ! on purpose.
               if (csm_extrap_on .and. .not. csm_probe_active .and.       &
                   .not. csm_probe_replay) then
                  if (csm_x_pending) then
                     call coupled_pair_extrapolation_verdict
                     csm_x_pending = .false.
                  endif
                  if (.not. csm_x_blocked) call coupled_pair_geometric_limit
               endif

               ! (f) This pass's increment becomes the history the next
               ! pass's estimate is taken from, unless the sequence has
               ! been broken deliberately -- the extrapolation moved the
               ! pair, or undid such a move -- in which case a fresh pair
               ! has to be collected before another estimate can be taken.
               if (csm_geom_seq_broken) then
                  csm_geom_have_prev  = .false.
                  csm_geom_seq_broken = .false.
               else
                  csm_geom_eT_prev(1:N)   = csm_geom_eT(1:N)
                  csm_geom_ec_prev(1:N,:) = csm_geom_ec(1:N,:)
                  csm_geom_have_prev = .true.
               endif

            enddo coupled_source
            ! The other end of that pair, with the source the loop settled
            ! on.  dt_loc and not dt: a local-dt run advances each cell by
            ! its own interval.
            call attempted_step_thermal_energy_after_sources(u, dt_loc,    &
                                                             heat, cool)

            ! WHAT THIS PASS'S NONLINEAR SOLVES LEFT BEHIND, for the
            ! integration-error estimate to judge itself against: the
            ! coupled pair's remaining composition and temperature errors
            ! in their own norms, and the relative row imbalance the
            ! carrier system stopped at. Recorded on every step; read only
            ! inside a sampled macrostep.
            call carrier_transport_diagnostics(ct_steps, ct_resid,        &
                                               ct_nlim, ct_worst)
            ! CELL BY CELL. csm_err_c and csm_err_T are the pair's error in
            ! the max norm of the whole grid; the estimate compares a
            ! difference at ONE cell, so the error handed over is the same
            ! estimate formed at each cell: the safety-factored geometric
            ! factor the loop measured, on that cell's own increment. The
            ! factor is read back from the two numbers the loop already
            ! carries rather than recomputed, so the two cannot disagree.
            if (csm_dcomp .gt. 0.0d0) then
               csm_inner_fac = csm_err_c/csm_dcomp
            else
               csm_inner_fac = 1.0d0
            endif
            csm_inner_c = csm_inner_fac*csm_dc_cell
            if (csm_dT .gt. 0.0d0) then
               csm_inner_fac = csm_err_T/csm_dT
            else
               csm_inner_fac = 1.0d0
            endif
            csm_inner_T = csm_inner_fac*csm_dT_cell
            call attempted_step_note_inner_error(csm_inner_c, csm_inner_T,&
                                                 ct_resid)

            n_csm_passes_total = n_csm_passes_total + csm_passes
            n_csm_calls = n_csm_calls + 1
            if (csm_passes .gt. n_csm_passes_worst) then
               n_csm_passes_worst = csm_passes
               n_csm_worst_step   = count
            endif
            if (.not. csm_ok) then
               n_csm_iter_cap = n_csm_iter_cap + 1
               if (n_csm_iter_cap .le. 5)                                  &
                  write(*,'(a,i0,a,i0,a,es10.3,a,es10.3)')                 &
                    ' (coupled source step) pass cap reached at step ',    &
                    count, ' after ', csm_passes,                          &
                    ' passes: dT/T= ', csm_dT, '  dcomp= ', csm_dcomp
            endif

            ! The injection points of the attempted-step controller. They sit
            ! after the coupled step because the loop calls each of the two
            ! solves more than once and a refusal is a statement about the
            ! operation, not about one of its passes. Row 8 no longer holds
            ! the composition projection -- it is the point at which the
            ! composition and the energy have become one state.
            if (attempted_step_injected_refusal(as_op_ioniz_eq)) then
               as_inject_op = as_op_ioniz_eq;  exit trial
            endif
            if (attempted_step_injected_refusal(as_op_composition)) then
               as_inject_op = as_op_composition;  exit trial
            endif
            if (attempted_step_injected_refusal(as_op_energy)) then
               as_inject_op = as_op_energy;  exit trial
            endif

            ! THE ENERGY UPDATE REFUSES THE ATTEMPT IN BOTH MODES. A
            ! status other than ENERGY_UPDATE_OK means the solve found no
            ! temperature at all for at least one cell: nothing was
            ! assembled from it, so there is no state to march on. The
            ! leniency the initialization mode grants elsewhere
            ! (docs/a0_run_mode_contract_20260906.md section 2) is about an
            ! iterate that exists, and a temperature that does not exist is
            ! not one. The controller restores the checkpoint and retakes
            ! the step at half dt.
            if (csm_energy_status .ne. ENERGY_UPDATE_OK) then
               as_trial_reason = as_reject_energy
               as_trial_op     = as_op_energy
               exit trial
            endif

            ! A coupled step that did not reach its fixed point is a state
            ! whose composition and temperature do not belong to each other.
            ! In phys mode that refuses the attempt, and the controller
            ! retakes it at half dt; in init mode the state is a numerical
            ! iterate (contract section 2), the exhaustion is recorded above
            ! and the march goes on. Unlike the energy failure above, the
            ! pair the loop returns IS a state: both halves were solved, and
            ! what the cap says is that they have not stopped moving.
            if (.not. csm_ok .and. run_mode .eq. run_mode_phys) then
               as_trial_reason = as_reject_source_fixed_point
               as_trial_op     = as_op_ioniz_eq
               exit trial
            endif

            if (upmap_n .gt. 0) then
               u_umB = u
               upmap_dfsp = maxval(abs(f_sp(1:N,:) - f_sp_umSave(1:N,:)))
            endif

		call Apply_BC(u)

            if (attempted_step_injected_refusal(as_op_bc)) then
               as_inject_op = as_op_bc;  exit trial
            endif

            if (upmap_n .gt. 0) u_umC = u

            ! Molecular transport (viscous momentum diffusion + its
            ! dissipation + heat conduction), Crank-Nicolson, as an
            ! operator-split stage. The spatial operators are the same ones
            ! assemble_residual subtracts, so the fixed point of the marching
            ! loop is the zero of the steady residual the Newton solver drives
            ! down. No-op unless "Viscosity:"/"Conduction:" were given, so a
            ! run without those keys is byte-identical to the inviscid code.
            ! Serial: the tridiagonal solves recur along r.
            if (transport_active()) then
               call U_to_W(u,W)
               rho = W(1,:);  v = W(2,:);  p = W(3,:)
               call comp_T_from_p(p,n_tot,ne,T)
               call viscous_conduction_step(u,W,T,n_tot+ne,dt_loc,count)
               call Apply_BC(u)
            endif

            if (attempted_step_injected_refusal(as_op_conduction)) then
               as_inject_op = as_op_conduction;  exit trial
            endif

            if (upmap_n .gt. 0) u_umD = u

            ! Periodic Shapiro low-pass filter to damp the gravity-unbalanced
            ! sound waves (base breathing), as in CETIMB (Koskinen et al. 2013a).
            if (shapiro_eps .gt. 0.0d0 .and.                             &
                mod(count, shapiro_every) .eq. 0) then
               call shapiro_filter(u)
               call Apply_BC(u)
               ! B6 CATEGORY 4, RECORDED (advisor decision 8). The filter
               ! alters the adopted state with no source term behind it, so
               ! an active filter is a validity state of every step it runs
               ! on. It is carried inside the trial until B3b budgets the
               ! energy it removes.
               n_shapiro_applied = n_shapiro_applied + 1
            endif

            if (attempted_step_injected_refusal(as_op_shapiro)) then
               as_inject_op = as_op_shapiro;  exit trial
            endif

            exit trial
            enddo trial

            !------------------------------------------------!
            ! ============ THE ADOPTION BOUNDARY ============
            ! The first point at which the physical state of the step is
            ! complete: row 12 is the last operation that writes u, and
            ! everything below only reads it. The acceptance predicate of
            ! design section 3.3 is evaluated here, on the trial state.
            ! ONE evaluation of the primitive state, not two. The
            ! post-adoption block below reads exactly these arrays: calling
            ! U_to_W, get_species_densities and comp_T_from_p a second time
            ! there would put the caloric composition through
            ! caloric_state_from_composition twice in one step, and the two
            ! passes are not the same arithmetic when the mixture is
            ! molecular.
            call U_to_W(u,W)
            rho = W(1,:)
            v   = W(2,:)
            p   = W(3,:)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                       nheiii,nheiTR,nm,ne,n_tot)
            call comp_T_from_p(p,n_tot,ne,T)

            if (as_trial_reason .ne. 0) then
               ! The trial left early at an operation that knows its own
               ! reason: the coupled source step of rows 7 to 9.
               as_verd%accepted  = .false.
               as_verd%reason    = as_trial_reason
               as_verd%operation = as_trial_op
               if (as_trial_reason .eq. as_reject_energy) then
                  as_verd%what   = 'the energy update found no'//          &
                                   ' temperature for at least one cell'
               else
                  as_verd%what   = 'the temperature and the composition'// &
                                   ' were still moving at the pass cap'
               endif
            else if (as_inject_op .ne. 0) then
               ! The trial left early: either a refusal was injected after a
               ! named operation, or the nested positivity retry exhausted
               ! its own budget (as_inject_op < 0).
               as_verd%accepted  = .false.
               if (as_inject_op .gt. 0) then
                  as_verd%reason    = as_reject_injected
                  as_verd%operation = as_inject_op
                  as_verd%what      = 'a refusal was injected after this'//&
                                      ' operation'
               else
                  as_verd%reason    = as_reject_hydro_stage
                  as_verd%operation = -as_inject_op
                  as_verd%what      = 'the nested positivity retry'//      &
                                      ' exhausted its bisections'
               endif
            else if (.not. as_carrier_ok) then
               as_verd%accepted  = .false.
               as_verd%reason    = as_reject_carrier
               as_verd%operation = as_op_carriers
               as_verd%what      = 'the carrier transport interval was'//  &
                                   ' not covered'
            else
               call attempted_step_evaluate(as_chk, u, W, T, f_sp, heat,  &
                        cool, dt, u_hyd, L_hyd,                           &
                        n_cells_without_chemical_root(last_sweep%acc_n),  &
                        as_carrier_ok, as_verd)
            endif

            ! THE ENERGY IDENTITY T1.5 (design section 9 decision 7): a
            ! sum over cells of the thermal energy change of the COUPLED
            ! SOURCE STEP ALONE against the source that step carries,
            ! printed every step in phys mode and a GATE there -- a failed
            ! closure refuses the attempt with as_reject_source_energy at
            ! operation 9.  The residual is the closure of that source step
            ! and no longer the size of the composition projection, which
            ! row 8 no longer makes.  d_u_th_transport is the thermal energy
            ! every other operation of the attempted step moved, measured
            ! between the two energies the loop is bracketed by rather than
            ! modeled, and it is on the right of the identity because T1.5
            ! constrains a local source step at fixed volume and fixed mass.
            ! d_u_form is the formation reservoir's own change, reported
            ! beside the row and not a term of it (the heating and cooling
            ! already carry every collisional transfer).
            if (as_verd%identity_evaluated .and.                          &
                run_mode .eq. run_mode_phys)                              &
               write(*,'(A,I0,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5,'//     &
                     'A,ES12.5)')                                         &
                    '   energy-identity: step=', count,                   &
                    ' d_u_th=', as_verd%d_u_th,                           &
                    ' d_u_form=', as_verd%d_u_form,                       &
                    ' Q_ext_dt=', as_verd%q_ext_dt,                       &
                    ' d_u_th_transport=', as_verd%d_u_th_transport,       &
                    ' residual/scale=', as_verd%identity_residual         &
                                        /as_verd%identity_scale
            if (as_verd%identity_evaluated .and.                          &
                run_mode .eq. run_mode_phys)                              &
               write(*,'(A,I0,A,ES12.5,A,I0,A,I0,A,L1)')                  &
                    '   hydro-time-row: step=', count,                    &
                    ' max=', as_verd%hydro_row_max,                       &
                    ' cell=', as_verd%hydro_jworst,                       &
                    ' row=', as_verd%hydro_kworst,                        &
                    ' above_tol=', as_verd%hydro_row_above

            if (as_verd%accepted) exit outer_attempt

            ! ---- A REJECTED ATTEMPT ----
            call attempted_step_note_step(as_verd)
            ! WHICH TRANSACTION CARRIES THE RETRY BUDGET. Outside a sampled
            ! macrostep the unit is the single step, and the retry shortens
            ! its dt in place. Inside one the unit is the macrostep: the
            ! three passes must share the endpoint t + H, so a shorter
            ! interval in one pass is not admissible and the whole
            ! transaction restarts at a shorter H.
            if (n_err_passes .eq. 3) then
               n_retry_now = n_macro_retry
            else
               n_retry_now = n_step_retry
            endif
            if (n_retry_now .eq. 0)                                       &
               n_steps_with_rejection = n_steps_with_rejection + 1
            ! The history this transaction is building: which interval
            ! was tried, and what refused it. The interval recorded is the
            ! one the retry budget is spent on, so inside a sampled
            ! macrostep it is H and not the half a pass covered. Read only
            ! by the exhaustion report.
            if (n_retry_now .ge. 0 .and.                                  &
                n_retry_now .le. n_step_retry_max) then
               if (n_err_passes .eq. 3) then
                  as_hist_dt(n_retry_now)  = dt_macro
               else
                  as_hist_dt(n_retry_now)  = dt
               endif
               as_hist_reason(n_retry_now) = as_verd%reason
               as_hist_op(n_retry_now)     = as_verd%operation
            endif
            write(*,'(A,I0,A,A,A,A)') '   step ', count,                  &
                 ' REFUSED at ', trim(attempted_step_operation_text(      &
                 as_verd%operation)), ': ',                               &
                 trim(attempted_step_reason_text(as_verd%reason))
            if (as_verd%jworst .gt. 0)                                    &
               write(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)')                   &
                    '     cell ', as_verd%jworst, ', row ',               &
                    as_verd%kworst, ', measure ', as_verd%measure,        &
                    ', tolerance ', as_verd%tol

            if (n_retry_now .ge. n_step_retry_max)                       &
               call stop_on_exhausted_retry_budget

            if (n_err_passes .eq. 3) then
               ! A REFUSAL INSIDE A SAMPLED MACROSTEP REJECTS THE WHOLE
               ! TRANSACTION. Shortening this pass would leave the three
               ! passes covering different intervals, and their difference
               ! would then not be a temporal error at a common time, so
               ! the state the macrostep began at is restored, H is
               ! shortened, and both comparison trajectories are recomputed
               ! from that state. No physical accumulation and no counter
               ! of the clock has moved: the adoption boundary is below the
               ! macrostep loop.
               call attempted_step_checkpoint_restore(as_entry, u, f_sp,  &
                                                      heat, cool, eta)
               if (.not. attempted_step_checkpoint_matches(as_entry, u,   &
                        f_sp, heat, cool, eta, as_bad_item)) then
                  write(*,'(A,A)') ' ATTEMPTED-STEP RESTORE INCOMPLETE:'//&
                       ' first item that did not come back: ',            &
                       trim(as_bad_item)
                  error stop 1
               endif
               call rebuild_state_from_checkpoint
               dt_reduced = attempted_step_reduced_dt(dt_macro,           &
                                                      as_verd%reason)
               ! In seconds, the unit the clock is printed in, so that
               ! the interval refused and the interval the clock later
               ! accepts are the same kind of number.
               write(*,'(A,I0,A,A,A,A,A,ES23.16,A,ES23.16)')              &
                    '   macrostep-reject: step=', count, ' pass=',        &
                    trim(err_pass_name), ' reason=',                      &
                    trim(attempted_step_reason_text(as_verd%reason)),     &
                    ' H=', dt_macro*R0/v0,                                &
                    ' H_new=', dt_reduced*R0/v0
               if (dt_macro .gt. 0.0d0)                                   &
                  dtloc_macro = dtloc_macro*(dt_reduced/dt_macro)
               dt_macro = dt_reduced
               n_macro_retry = n_macro_retry + 1
               dt_bounded_by_rejection = .true.
               cycle macrostep
            endif

            call attempted_step_checkpoint_restore(as_chk, u, f_sp, heat, &
                                                   cool, eta)
            ! THE ROUND TRIP IS NOT ASSUMED. A restore that missed an item
            ! is a silent corruption of the next attempt, so it is checked
            ! here and it is loud (b1 AT-7 (i)).
            if (.not. attempted_step_checkpoint_matches(as_chk, u, f_sp,  &
                     heat, cool, eta, as_bad_item)) then
               write(*,'(A,A)') ' ATTEMPTED-STEP RESTORE INCOMPLETE:'//   &
                    ' first item that did not come back: ',               &
                    trim(as_bad_item)
               error stop 1
            endif
            call rebuild_state_from_checkpoint
            dt_reduced = attempted_step_reduced_dt(dt, as_verd%reason)
            if (dt .gt. 0.0d0) dt_loc = dt_loc*(dt_reduced/dt)
            dt = dt_reduced
            n_step_retry = n_step_retry + 1
            dt_bounded_by_rejection = .true.

            enddo outer_attempt
            enddo err_pass

            if (n_err_passes .eq. 3) then
               ! e = the local error of the FULL step this macrostep
               ! keeps, over its tolerance: the difference between the step
               ! at H and the two at H/2, divided by (1 - 2**-p) for the
               ! retained trajectory and by the row's own scale. p is the
               ! order of the INTEGRATOR THAT CARRIES THE ROW'S CLASS,
               ! three for the conservative rows and for the electron
               ! fraction and the temperature, one for the transported
               ! fractions (attempted_step.f90, err_order_p); the class
               ! whose lower bound the gate reads also sets the exponent of
               ! the reduction a refusal applies.
               call attempted_step_error_estimate(u_err_full,             &
                         f_sp_err_full, T_err_full, ne_err_full,          &
                         u, f_sp, T, ne, as_est)
               e_err  = as_est%e_resolved
               je_err = as_est%jworst
               ke_err = as_est%kworst
               write(*,'(A,I0,A,ES12.5,A,ES12.5,A,ES12.5,A,A,A,I0,'//    &
                     'A,I0,A,L1,A,A,A,I0)')                               &
                    '   integration-error estimate: step=', count,        &
                    ' e=', as_est%e, ' e_lower=', as_est%e_lower,         &
                    ' e_upper=', as_est%e_upper,                          &
                    ' class=', trim(err_class_name(max(1,                 &
                                    as_est%worst_class))),                &
                    ' cell=', as_est%jworst, ' row=', as_est%kworst,      &
                    ' resolved=', as_est%resolved,                        &
                    ' deciding=', trim(err_class_name(max(1,              &
                                    as_est%deciding_class))),             &
                    ' p=', err_order_p(max(1, as_est%deciding_class))
               ! ONE LINE PER ROW CLASS, so that a rejection can be read
               ! back to the physics that caused it and a class that
               ! carries no difference is visibly absent.
               if (trace_step_clock) then
                  do i_ecl = 1, n_err_class
                     if (.not. as_est%class_present(i_ecl)) cycle
                     write(*,'(A,A,A,ES12.5,A,ES12.5,A,ES12.5,'//         &
                           'A,I0,A,I0,A,L1)')                             &
                          '     err-class: ',                             &
                          trim(err_class_name(i_ecl)), ' e=',             &
                          as_est%e_class(i_ecl), ' diff=',                &
                          as_est%d_class(i_ecl), ' scale=',               &
                          as_est%s_class(i_ecl), ' cell=',                &
                          as_est%j_class(i_ecl), ' row=',                 &
                          as_est%k_class(i_ecl), ' resolved=',            &
                          as_est%class_resolved(i_ecl)
                  enddo
               endif
               ! AN UNRESOLVED ESTIMATE JUDGES NOTHING. Where the inner
               ! nonlinear error of the passes is not below
               ! err_inner_fraction of the difference measured, the
               ! difference is the solvers' own error and not a temporal
               ! discretization error, so refusing on it would be a false
               ! rejection and accepting on it a false acceptance. The
               ! macrostep is adopted, the outcome is named, and the run
               ! reports how many of its steps were adopted this way: a
               ! trajectory with unresolved estimates in it is not an
               ! error-controlled trajectory.
               if (.not. as_est%resolved) then
                  write(*,'(A,I0,A,A,A,ES12.5,A,ES12.5,A,ES12.5)')        &
                       '   integration-error UNRESOLVED: step=', count,   &
                       ' first unresolved class=',                        &
                       trim(err_class_name(max(1,                         &
                            as_est%unresolved_class))),                   &
                       ' inner/difference=',                              &
                       as_est%inner_over_difference,                      &
                       ' required below ', err_inner_fraction,            &
                       ' e_upper=', as_est%e_upper
               endif
               if (.not. ieee_is_finite(e_err) .or. e_err .gt. 1.0d0) then
                  ! THE ESTIMATE IS A REJECTION CRITERION. A scaled error
                  ! above one says the macrostep did not resolve the
                  ! trajectory to the tolerance it claims, so the state the
                  ! macrostep began at is restored and the whole
                  ! transaction is retaken at a shorter H. A NONFINITE
                  ! estimate is the same rejection and not a stop: the two
                  ! trajectories differing by a nonfinite amount is a
                  ! statement about this interval, not about the state,
                  ! and a state from which no interval works reaches the
                  ! named exhaustion above through the same retry budget.
                  ! The reduction policy uses the estimate only when it is
                  ! finite (attempted_step_reduced_dt).
                  as_verd%accepted  = .false.
                  as_verd%reason    = as_reject_int_error
                  as_verd%operation = as_op_boundary
                  as_verd%what      = 'the two comparison trajectories'// &
                                      ' of this macrostep differ by more'&
                                      //' than its tolerance'
                  err_pass_name = 'the error estimate'
                  call attempted_step_note_step(as_verd)
                  n_retry_now = n_macro_retry
                  if (n_retry_now .eq. 0)                                 &
                     n_steps_with_rejection = n_steps_with_rejection + 1
                  if (n_retry_now .ge. 0 .and.                            &
                      n_retry_now .le. n_step_retry_max) then
                     as_hist_dt(n_retry_now)     = dt_macro
                     as_hist_reason(n_retry_now) = as_reject_int_error
                     as_hist_op(n_retry_now)     = as_op_boundary
                  endif
                  if (n_retry_now .ge. n_step_retry_max)                  &
                     call stop_on_exhausted_retry_budget
                  call attempted_step_checkpoint_restore(as_entry, u,     &
                                                         f_sp, heat,      &
                                                         cool, eta)
                  ! THE ROUND TRIP IS NOT ASSUMED here either: the next
                  ! macrostep is computed from this state.
                  if (.not. attempted_step_checkpoint_matches(as_entry,   &
                           u, f_sp, heat, cool, eta, as_bad_item)) then
                     write(*,'(A,A)') ' ATTEMPTED-STEP RESTORE'//         &
                          ' INCOMPLETE: first item that did not come'//   &
                          ' back: ', trim(as_bad_item)
                     error stop 1
                  endif
                  call rebuild_state_from_checkpoint
                  dt_reduced = attempted_step_reduced_dt(dt_macro,        &
                                                    as_reject_int_error)
                  write(*,'(A,I0,A,A,A,A,A,ES23.16,A,ES23.16)')           &
                       '   macrostep-reject: step=', count, ' pass=',     &
                       trim(err_pass_name), ' reason=',                   &
                       trim(attempted_step_reason_text(                   &
                            as_reject_int_error)),                        &
                       ' H=', dt_macro*R0/v0,                             &
                       ' H_new=', dt_reduced*R0/v0
                  if (dt_macro .gt. 0.0d0)                                &
                     dtloc_macro = dtloc_macro*(dt_reduced/dt_macro)
                  dt_macro = dt_reduced
                  n_macro_retry = n_macro_retry + 1
                  dt_bounded_by_rejection = .true.
                  cycle macrostep
               endif
               ! THE TRAJECTORY THE RUN KEEPS IS THE ONE AT H, which is
               ! the state the acceptance predicate was evaluated on and
               ! the state the estimate above refers to; the two half steps
               ! were the measurement, not the step.
               call attempted_step_checkpoint_restore(as_full, u, f_sp,   &
                                                      heat, cool, eta)
               call rebuild_state_from_checkpoint
               dt = dt_macro;  dt_loc = dtloc_macro
               error_estimate_seconds = error_estimate_seconds            &
                                        + (omp_get_wtime() - t_err0)
            endif

            exit macrostep
            enddo macrostep

            ! THE STEP IS ADOPTED HERE, and this is the only place the
            ! physical clock moves (docs/a0_run_mode_contract_20260906.md
            ! section 5). Physical time advances by the GLOBAL dt of the
            ! step that was accepted -- never by the minimum pseudo-step,
            ! which belongs to no single trajectory, and never by an
            ! attempt the controller threw away, whose state the run does
            ! not keep.
            n_steps_accepted = n_steps_accepted + 1
            dt_last_accepted = dt
            ! The code's time unit is R0/v0 (dr_j is in R0, v in v0), so
            ! the increment is stated in seconds here and the header needs
            ! no unit of its own.
            if (run_mode .eq. run_mode_phys) t_phys = t_phys + dt*R0/v0
            ! The state hash is the same functional the rollback test
            ! reads, so an exhaustion's restored state and the last state
            ! the clock accepted are comparable numbers.
            if (trace_step_clock)                                         &
               write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES23.16,A,ES23.16,'//      &
                     'A,ES23.16)')                                        &
                    '   step-clock: count=', count,                       &
                    ' outer_attempts=', n_outer_attempts,                 &
                    ' attempted=', n_steps_attempted,                     &
                    ' accepted=', n_steps_accepted,                       &
                    ' dt=', dt*R0/v0, ' t_phys=', t_phys,                 &
                    ' state=', attempted_step_checkpoint_checksum(u,      &
                                        f_sp, heat, cool, eta)

            if (upmap_n .gt. 0) call update_map_end_step

            ! The primitive state is the one the adoption boundary above
            ! evaluated on the state it adopted; only the energy row is
            ! taken here.
            E   = u(3,:)

            ! Is the lower boundary still admitting SUBSONIC inflow? Measured
            ! here because W now holds the state every Apply_BC of this step
            ! has finished writing. Diagnostic only: it reads W and changes
            ! nothing, so a run that never crosses is the run an unguarded
            ! build produces.
            call check_base_inflow_is_subsonic(W,count)
            
            ! Base-cell startup diagnostic (first 500 steps): trace r, n, v, T,
            ! heat, cool for the lowest cells to expose IC-startup transients.
            if (diag_base .and. count .le. 500) then
               do j = 1,6
                  write(778,'(I7,I4,1X,F10.6,5(1X,ES13.6))') count, j, r(j), &
                       W(1,j)*n0, W(2,j)*v0, T(j)*T0, heat(j), cool(j)
               enddo
            endif

		! Evaluate momentum
		mom = rho*v*r*r

            !---------------------------------------------------!
        
            !--- Loop counters and escape condition ---!
            
            ! Update counter
            count = count + 1
            
            ! The convergence measure of the marching: the radial spread of
            ! the mass flux rho v r^2 over the window r >= r_esc, evaluated
            ! ONCE per step by the one functional (mass_flux_spread, below)
            ! that the reported flux spread also reads.
            du = mass_flux_spread(mom)

            ! Arm the du triggers the first time the flux spread is seen at or
            ! above their thresholds; from then on only a descending crossing
            ! can fire them (see the declarations of du_*_armed).
            if (du .ge. du_th .and. .not.du_stop_armed) then
               du_stop_armed = .true.
               write(*,'(A,ES12.4,A,I0)') '    -> du stop armed at du =',   &
                                          du, ', step ', count
            endif
            if (du .ge. newton_du_switch .and. .not.du_newton_armed) then
               du_newton_armed = .true.
               write(*,'(A,ES12.4,A,I0)') '    -> Newton hand-off armed '// &
                                          'at du =', du, ', step ', count
            endif

            ! RELATIVE CHANGE OF THE CONSERVED STATE OVER THE STEP, as an
            ! infinity norm over the convergence window:
            !
            !     dtu = max over rows k and cells j in [j_min,N] of
            !           | 1 - u(k,j)/u_old(k,j) | .
            !
            ! The conserved state is written by the operator split alone --
            ! the Riemann fluxes, the geometric and gravitational sources,
            ! the semi-implicit energy update, the transport and diffusion
            ! steps and Apply_BC. A diagnostic evaluated afterwards reads u;
            ! it must not add to it, or the momentum row carries an increment
            ! no flux and no source produced, which then enters the next
            ! step's hydro, the steady residual and the mass flux the
            ! convergence test measures.
            !
            ! The denominator is guarded where it belongs, at the point of
            ! use: a row that was EXACTLY zero at the start of the step has
            ! no relative change defined, and the state did change, so the
            ! window is taken as unbounded rather than dropped. Dropping it
            ! would let the "state stopped changing" stop fire on a cell
            ! whose momentum moved away from zero.
            dtu = 0.0d0
            do krow = 1,3
               if (any(u_old(krow,j_min:N) .eq. 0.0d0)) then
                  dtu = huge(1.0d0)
               else
                  dtu = max(dtu, maxval(abs(1.0d0 - u(krow,j_min:N)         &
                                                   /u_old(krow,j_min:N))))
               endif
            enddo

            ! Periodic steady-residual monitor (Resid tol option). R = du/dt
            ! evaluated fresh on the current state (one Reconstruct + RK_rhs;
            ! reuses this step's heat/cool, which near convergence are
            ! consistent with the final u). du/dtu can be small at a premature
            ! operator-split balance while R is large (premature WASP golden:
            ! du,dtu tiny but R_energy ~ 30), so R is the trustworthy gate.
            ! Reference diagnostics every N_resid steps. Per Caldiroli (2021,
            ! ATES) and Koskinen (2013a, CETIMB) the CONVERGENCE DECISION is
            ! flux-based: ATES stops at d(Mdot)/Mdot < 1e-3, which here is the
            ! du<du_th test (du is the radial spread of rho*v*r^2); CETIMB asks
            ! that rho*v*r^2 = F_c be constant with altitude. The steady residual
            ! ||R|| (the maximum over cells of a cell's own scaled
            ! reported FOR REFERENCE ONLY; it gates the stop solely when the user
            ! explicitly requests it via "Resid tol:" (resid_th>0).
            !
            ! WHEN IT DOES GATE, IT GATES ON THE SAME TEST THE STEADY SOLVERS
            ! USE. steady_gates_met (steady_residual.f90) is the one definition
            ! of an accepted steady state -- ||R|| < resid_th AND the mass flux
            ! flat over r >= r_flux -- so a marched state and a JFNK solution
            ! are accepted by the same statement. Before this was shared, the
            ! marching loop accepted a WASP-121b state at a flux spread of
            ! 9.9e-2, twenty times the default threshold the JFNK was
            ! rejecting 1.9e-2 on.
            if (mod(count, N_resid) .eq. 0) then
               call assemble_residual(u, n_tot + ne, heat, cool, Rres)
               call residual_norms(Rres, u, resid_c)
               resid_max   = maxval(resid_c)
               ! The spread the stop measured on this state, not a second
               ! functional of the same window: one quantity, one number.
               flux_spread = du
               ! carrier_is_unknown = .false.: the marching loop's carriers
               ! are moved by the operator-split transport step, not solved
               ! for, so this state has no carrier residual of a solve to
               ! gate on. A coupled-carrier run reaches its answer through
               ! the JFNK finish, which does state one.
               n_no_root_now = n_cells_without_chemical_root(               &
                                                last_sweep%acc_n)
               gates_met_now = steady_gates_met(resid_max, u, resid_th,      &
                                                fspread_gate,               &
                                                .false., 0.0d0,             &
                                                n_no_root_now)
               ! THE CERTIFICATION OF THE STATE THIS ROUTE WOULD DECLARE
               ! SOLVED. It measures every equation the configuration makes
               ! independent, on the state Rres was just assembled from, and
               ! decides nothing about the marching: a state that refuses
               ! certification is still marched on and still written, and the
               ! '# coupling:' header of the file says which it was. Reported
               ! in full the first time and whenever the verdict changes, so
               ! a long run states its verdict without printing one block per
               ! monitor step.
               call certification_evaluate(cert_context_stationary, u,      &
                        Rres, f_sp, resid_th, n_no_root_now, .true.,        &
                        cert_now)
               if (cert_first .or.                                          &
                   (cert_now%certified .neqv. cert_was_certified)) then
                  write(cert_label,'(A,I0)') 'marching state at step ',      &
                        count
                  call certification_report_write(cert_now, cert_label)
                  cert_first         = .false.
                  cert_was_certified = cert_now%certified
               endif
               ! Reporting only, and it does not enter gates_met_now above:
               ! the gate window cannot show how far the non-flatness reaches.
               call flux_spread_above_radius(u, 1.03d0, fspread_103, dum)
               call flux_spread_above_radius(u, 1.10d0, fspread_110, dum)
               write(*,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
                    '   [diag] step ', count,                                 &
                    '  flux rho*v*r^2 spread=', flux_spread,                  &
                    ' (gate window ', fspread_gate,                           &
                    '; r>=1.03 ', fspread_103,                                &
                    ', r>=1.10 ', fspread_110,                                &
                    ')   ||R||(ref)=', resid_max
            endif

            ! P54/(AD) time series of the FACE mass flux at fixed radii
            ! (env gate; assemble_residual first so face_flux belongs to
            ! the current u).
            if (p54_ts_every .gt. 0) then
               if (mod(count, p54_ts_every) .eq. 0) then
                  call assemble_residual(u, n_tot + ne, heat, cool, Rres)
                  call p54_flux_timeseries_append
               endif
            endif

            ! Mass-flux LEVEL stability: du is the spread of rho*v*r^2 and is
            ! blind to a uniform drift of the whole profile; require the mean
            ! level over the wind to be unchanged across the last N_stall
            ! steps before any converged/stalled stop is accepted.
            if (lev_th .gt. 0.0d0) then
               lev_now = sum(abs(mom(j_min:N)))/dble(N - j_min + 1)
               if (lev_count .ge. N_stall) then
                  lev_rel = abs(lev_now - lev_hist(mod(lev_count, N_stall))) &
                            /max(lev_now, 1.0d-30)
                  is_level_stable = (lev_rel .lt. lev_th)
               endif
               lev_hist(mod(lev_count, N_stall)) = lev_now
               lev_count = lev_count + 1
            endif

            ! ---- PLM -> WENO3 continuation: advance lambda ----
            ! Placed before the stage logic below so that the hand-off step
            ! itself only starts the ramp; lambda is advanced from the step
            ! after it, once one step of the new operator has been measured.
            if (in_lambda_ramp) then
               lam_step_count = lam_step_count + 1
               if (recon_lambda_adaptive) then
                  if (dtu .le. recon_lambda_dtu_tol*max(dtu_ref,1.0d-300)) then
                     lam_good     = lam_good + 1
                     recon_lambda = min(1.0d0,                            &
                                        recon_lambda + recon_lambda_step)
                     if (lam_good .ge. 5) then
                        recon_lambda_step = min(recon_lambda_step0,       &
                                                2.0d0*recon_lambda_step)
                        lam_good = 0
                     endif
                  else
                     ! Too big a step in the state: halve the step in lambda
                     ! and do not advance it this time.
                     lam_good          = 0
                     recon_lambda_step = max(1.0d-6,                      &
                                             0.5d0*recon_lambda_step)
                     n_lambda_backoff  = n_lambda_backoff + 1
                  endif
               else
                  recon_lambda = min(1.0d0, recon_lambda + recon_lambda_step0)
               endif
               if (recon_lambda .ge. 1.0d0) then
                  ! Done: disarm, so that the JFNK finish and every gate see
                  ! the pure WENO3 operator with no blending left in it.
                  recon_lambda    = 1.0d0
                  in_lambda_ramp  = .false.
                  recon_lambda_on = .false.
                  rec_method      = 'WENO3'
                  use_plm         = .false.
                  use_weno3       = .true.
                  write(*,'(A,I0,A,I0,A,I0,A,ES12.4)')                    &
                     '    -> PLM -> WENO3 continuation reached lambda = '//&
                     '1 at step ', count, ' (', lam_step_count,           &
                     ' ramp steps, ', n_lambda_backoff,                   &
                     ' back-offs), du =', du
               else if (mod(lam_step_count, 25) .eq. 0) then
                  write(*,'(A,I0,A,F8.5,A,ES10.3,A,ES10.3)')              &
                     '    [lambda] step ', count, '  lambda=',            &
                     recon_lambda, '  du=', du, '  dtu=', dtu
               endif
            endif

            ! Adjust logicals for the loop (two-stage aware)
            if (in_plm_stage) then
               ! Stage 1 (PLM): do not converge; hand off to WENO3 once du has
               ! dropped below du_th_plm OR the PLM du has plateaued (stall).
               if (count .gt. 1) then
                  ! THE SENTINEL IS KEPT OUT OF THE ARITHMETIC. du_prev
                  ! starts at huge(1.0d0) and is re-set to it at each stage
                  ! switch, so the first evaluation after either forms
                  ! huge/du, which overflows to +Inf. Unmasked that is
                  ! harmless (+Inf is not below stall_tol, so the count
                  ! resets, which is what the code wants), but it raises the
                  ! overflow trap once per stage of every run. The nested
                  ! form takes the same branch without forming the value.
                  ! Not rewritten as abs(du-du_prev) < stall_tol*max(du,..):
                  ! multiplying instead of dividing can move the comparison
                  ! by an ulp at the threshold, and this form cannot.
                  if (du_prev .lt. huge(1.0d0)) then
                     if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                        stall_count = stall_count + 1
                     else
                        stall_count = 0
                     endif
                  else
                     stall_count = 0
                  endif
               endif
               du_prev = du
               if (du .lt. du_th_plm .or. stall_count .ge. N_stall) then
                  ! What the change of operator does to the equations, on the
                  ! state it is about to be applied to, before any step is
                  ! taken with it.
                  if (opdiff_mode .gt. 0)                                  &
                     call write_operator_difference(count, du)
                  if (recon_lambda_step0 .gt. 0.0d0) then
                     ! Walk the hand-off. The flags stay on PLM:
                     ! reconstruction_continuation_rhs supplies both
                     ! operators and combines them, and Apply_BC_W blends the
                     ! outer ghost, until lambda reaches 1 and the controller
                     ! above disarms both. dlambda = 1 puts lambda at 1 here,
                     ! which is the one-step switch this replaces.
                     in_lambda_ramp    = .true.
                     in_plm_stage      = .false.
                     recon_lambda_on   = .true.
                     recon_lambda      = min(1.0d0, recon_lambda_step0)
                     recon_lambda_step = recon_lambda_step0
                     lam_step_count    = 0
                     lam_good          = 0
                     n_lambda_backoff  = 0
                     dtu_ref           = dtu
                     stall_count       = 0
                     du_prev           = huge(1.0d0)
                     write(*,'(A,F8.5,A,ES12.4,A,I0)')                     &
                        '    -> started PLM -> WENO3 continuation at '//   &
                        'lambda =', recon_lambda, ', du =', du,            &
                        ', step ', count
                  else
                     rec_method   = 'WENO3'
                     ! Flip the discretization flags too (Source/Num_Fluxes/
                     ! RK_rhs/Apply_BC read these, not rec_method); without
                     ! this the WENO3 stage ran with PLM-style pressure
                     ! flux/source terms.
                     use_plm      = .false.
                     use_weno3    = .true.
                     in_plm_stage = .false.
                     stall_count  = 0
                     du_prev      = huge(1.0d0)
                     write(*,'(A,ES12.4,A,I0)') '    -> switched PLM -> WENO3 at du =', &
                                                du, ', step ', count
                  endif
                  if (opdiff_mode .ge. 2) then
                     write(*,*) '(EXHALE_main) EXHALE_OPDIFF>=2: '//       &
                                'operator difference written, stopping.'
                     stop
                  endif
               endif
               mass_flux_converged    = .false.
               step_change_converged  = .false.
               du_plateaued           = .false.
               steady_gates_converged = .false.
            else if (use_newton_solver) then
               ! Newton warm-up: never stop on du/dtu/stall (du dips are
               ! premature while the energy residual is still large); the
               ! R-based switch below hands over to the JFNK finish.
               ! DO track the du plateau: runs whose flux metric stalls just
               ! above newton_du_switch (seen with He_diffusion: du frozen at
               ! ~1.08e-2 vs the 1e-2 switch for 1e6 steps) hand off on stall.
               if (count .gt. 1) then
                  ! THE SENTINEL IS KEPT OUT OF THE ARITHMETIC. du_prev
                  ! starts at huge(1.0d0) and is re-set to it at each stage
                  ! switch, so the first evaluation after either forms
                  ! huge/du, which overflows to +Inf. Unmasked that is
                  ! harmless (+Inf is not below stall_tol, so the count
                  ! resets, which is what the code wants), but it raises the
                  ! overflow trap once per stage of every run. The nested
                  ! form takes the same branch without forming the value.
                  ! Not rewritten as abs(du-du_prev) < stall_tol*max(du,..):
                  ! multiplying instead of dividing can move the comparison
                  ! by an ulp at the threshold, and this form cannot.
                  if (du_prev .lt. huge(1.0d0)) then
                     if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                        stall_count = stall_count + 1
                     else
                        stall_count = 0
                     endif
                  else
                     stall_count = 0
                  endif
               endif
               du_prev      = du
               mass_flux_converged    = .false.
               step_change_converged  = .false.
               du_plateaued           = .false.
               steady_gates_converged = .false.
            else if (resid_th .gt. 0.0d0) then
               ! Stage 2, residual-based convergence: the STEADY GATES decide
               ! convergence in place of du (which is blind to the
               ! operator-split energy imbalance). du/dtu kept for logging.
               ! The gates are steady_gates_met's pair, evaluated above -- the
               ! residual AND the flatness of rho v r^2 over the wind -- so
               ! this stop and the steady solvers' acceptance are the same
               ! statement about the same state.
               !
               ! It does NOT replace the stall net. A residual gate is one
               ! more stop condition, not a licence to switch the others off:
               ! a run whose du has stopped moving will not reach the gate by
               ! marching further, and with every other stop disabled it ran
               ! to count_max = 1e6. Measured: jfnk_cold reaches du = 1.5e-4
               ! and ptc_warm du = 5.5e-5 -- both flux-converged by any
               ! ordinary standard -- and neither crossed its residual gate in
               ! 40,000 steps. So the plateau detector stays live here, and
               ! the stop it produces is reported as what it is: a stall with
               ! the gate target unmet, never a converged residual.
               mass_flux_converged    = .false.
               step_change_converged  = .false.
               steady_gates_converged = gates_met_now
               if (count .gt. 1) then
                  ! THE SENTINEL IS KEPT OUT OF THE ARITHMETIC. du_prev
                  ! starts at huge(1.0d0) and is re-set to it at each stage
                  ! switch, so the first evaluation after either forms
                  ! huge/du, which overflows to +Inf. Unmasked that is
                  ! harmless (+Inf is not below stall_tol, so the count
                  ! resets, which is what the code wants), but it raises the
                  ! overflow trap once per stage of every run. The nested
                  ! form takes the same branch without forming the value.
                  ! Not rewritten as abs(du-du_prev) < stall_tol*max(du,..):
                  ! multiplying instead of dividing can move the comparison
                  ! by an ulp at the threshold, and this form cannot.
                  if (du_prev .lt. huge(1.0d0)) then
                     if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                        stall_count = stall_count + 1
                     else
                        stall_count = 0
                     endif
                  else
                     stall_count = 0
                  endif
               endif
               du_prev    = du
               du_plateaued = (stall_count .ge. N_stall) .and. is_level_stable
            else
               ! Stage 2 (WENO3) or single-stage: normal convergence + stall,
               ! each additionally gated on mass-flux level stability.
               mass_flux_converged   = (du .lt. du_th) .and. du_stop_armed  &
                                                       .and. is_level_stable
               step_change_converged = (dtu .lt. dtu_th) .and. is_level_stable
               steady_gates_converged = .false.
               ! Stall detection: du settled on a plateau
               if (count .gt. 1) then
                  ! THE SENTINEL IS KEPT OUT OF THE ARITHMETIC. du_prev
                  ! starts at huge(1.0d0) and is re-set to it at each stage
                  ! switch, so the first evaluation after either forms
                  ! huge/du, which overflows to +Inf. Unmasked that is
                  ! harmless (+Inf is not below stall_tol, so the count
                  ! resets, which is what the code wants), but it raises the
                  ! overflow trap once per stage of every run. The nested
                  ! form takes the same branch without forming the value.
                  ! Not rewritten as abs(du-du_prev) < stall_tol*max(du,..):
                  ! multiplying instead of dividing can move the comparison
                  ! by an ulp at the threshold, and this form cannot.
                  if (du_prev .lt. huge(1.0d0)) then
                     if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                        stall_count = stall_count + 1
                     else
                        stall_count = 0
                     endif
                  else
                     stall_count = 0
                  endif
               endif
               du_prev    = du
               du_plateaued = (stall_count .ge. N_stall) .and. is_level_stable
            endif

            ! Staged secondary ionization: the SvS85 coupling is applied only
            ! after the wind has first converged without it -- from a cold IC the
            ! secondary-ionization base feedback amplifies the startup transient
            ! into a runaway (NaN); from a converged state the coupling is benign
            ! and the run re-converges with the full physics included.
            if (use_sec_ion .and. .not.sec_ion_active .and.                   &
                (mass_flux_converged .or. step_change_converged .or.          &
                 du_plateaued .or. steady_gates_converged)) then
               sec_ion_active = .true.
               sec_flip_step  = count
               sec_ion_armed_step = count
               mass_flux_converged = .false.;  step_change_converged = .false.
               du_plateaued = .false.;  steady_gates_converged = .false.
               stall_count  = 0
               du_prev      = huge(1.0d0)
               ! restart the mass-flux level-stability history
               lev_count    = 0
               is_level_stable = (lev_th .le. 0.0d0)   ! gate disabled => always pass
               write(*,'(A,I0,A,ES10.2)') '    -> secondary ionization activated at step ', &
                                          count, ', du =', du
            else if (sec_flip_step .ge. 0 .and.                               &
                     count - sec_flip_step .lt. N_stall) then
               ! Hold any stop for N_stall steps after the flip: du reacts only
               ! once the base adjustment wave driven by the new coupling has
               ! formed, so an immediate stop would freeze a state that has not
               ! yet incorporated the secondary-ionization physics.
               mass_flux_converged = .false.;  step_change_converged = .false.
               du_plateaued = .false.;  steady_gates_converged = .false.
            endif

            ! Production Newton finish ("Solver: Newton"): once the cheap
            ! marching warm-up has flattened the wind to du < newton_du_switch
            ! (the FLUX metric -- the radial spread of rho*v*r^2 -- consistent
            ! with the flux-based convergence decision), hand over to the JFNK
            ! steady solver to polish the rest of the way, refresh the
            ! thermodynamic state from the solved u, and exit the loop (standard
            ! final outputs and post-processing follow). Newton's quadratic local
            ! convergence tightens du from ~1e-2 to <1e-3 in far fewer steps than
            ! continued marching.
            if (use_newton_solver .and. .not.in_plm_stage .and.          &
                du_newton_armed .and.                                     &
                (newton_fail_step .lt. 0 .or.                             &
                 count - newton_fail_step .ge. N_stall) .and.             &
                (du .lt. newton_du_switch .or.                            &
                 (stall_count .ge. N_stall .and.                          &
                  du .lt. 5.0d0*newton_du_switch))) then
               if (use_sec_ion .and. .not.sec_ion_active) then
               ! Newton may only be asked to solve the system the marched
               ! state approximately satisfies. Flipping the secondary-
               ! ionization coupling on in the same breath as starting JFNK
               ! hands the solver a state that was relaxed WITHOUT the
               ! coupling -- from there GMRES stagnates (info=2, no residual
               ! reduction; measured on HD 209458 b, 2026-08-10). So flip
               ! first, re-arm the marching stops exactly as the staged
               ! activation above does, and let the wind re-relax under the
               ! full physics; the JFNK hand-off re-fires once du crosses
               ! the switch again after the N_stall hold.
                  sec_ion_active = .true.
                  sec_flip_step  = count
                  sec_ion_armed_step = count
                  mass_flux_converged = .false.;  step_change_converged = .false.
                  du_plateaued = .false.;  steady_gates_converged = .false.
                  stall_count  = 0
                  du_prev      = huge(1.0d0)
                  lev_count    = 0
                  is_level_stable = (lev_th .le. 0.0d0)
                  write(*,'(A,I0,A,ES10.2)') '    -> secondary ionization '// &
                       'activated ahead of the Newton finish at step ',       &
                       count, ', du =', du
               else if (sec_flip_step .lt. 0 .or.                         &
                        count - sec_flip_step .ge. N_stall) then
               if (du .ge. newton_du_switch)                              &
                  write(*,'(A,ES10.2)') ' (EXHALE_main) du plateaued '//  &
                       'near the hand-off threshold; engaging JFNK at '// &
                       'du =', du
               resid_max = resid_th
               if (resid_max .le. 0.0d0) resid_max = 1.0d-5
               write(*,'(A,I0,A,ES10.2)') ' (EXHALE_main) Newton finish '// &
                    'at step ', count, ', target ||R|| <', resid_max
               ! A STATIONARY SOLVE IS CONTINUATION, WHATEVER MODE THE RUN
               ! IS IN (contract section 4). Its trials and iterates are
               ! numerical states, not states the flow passed through, so
               ! what they leave in the ledgers belongs to the
               ! initialization family and not to the history of accepted
               ! steps. The clock does not move across it either: no global
               ! dt is taken, so there is nothing to add.
               ledger_family_held = ledger_family
               ledger_family      = ledger_family_init
               call steady_wind_with_element_diffusion(                  &
                                     jfnk_outer_iterations_default,       &
                                     1.0d0, .true., j)
               ledger_family = ledger_family_held
               if (j .eq. 0) then
                  ! newton_finished is itself a stop condition of the loop,
                  ! so no other flag is set: the run reports the JFNK finish
                  ! and not a du criterion it never had to satisfy.
                  du_plateaued    = .false.
                  newton_finished = .true.
               else if (j .eq. 2) then
                  ! info = 2 IS NOT A SOLVER FAILURE. The solve converged;
                  ! what it says is that the state it hands back does not
                  ! meet the gate the loop top met, and it names the row
                  ! that refuses (B5pre). Marching on from it and retrying
                  ! makes the answer WORSE, not better: MEASURED by B5pre
                  ! on wasp_full_newton, the retry ends on a du plateau at
                  ! step 15886 with all three hydrodynamic rows refusing
                  ! (mass 1.078e-03, momentum 2.387e-01, energy 2.356e-01)
                  ! where the JFNK state refuses on the energy row alone
                  ! (4.221e-04 at cell 376).
                  !
                  ! So the JFNK state IS the state the run writes. It made a
                  ! stationarity claim and the certification refused it, and
                  ! the A2 rule for that is to write the state with
                  ! certified=F, name the refusing row, and exit 2. The
                  ! marching stop conditions are left alone; newton_finished
                  ! ends the loop.
                  du_plateaued    = .false.
                  newton_finished = .true.
                  newton_returned_uncertified = .true.
                  write(*,'(A)') ' (EXHALE_main) JFNK returned info = 2:'//&
                       ' the solve converged, and the state it hands back'
                  write(*,'(A)') '   does not meet the gate its last'//    &
                       ' iterate met. The refusing row is named above.'
                  write(*,'(A)') '   That state is written as the final'// &
                       ' state, uncertified: retrying from it produces a'
                  write(*,'(A)') '   state the certification measures as'//&
                       ' worse, not better.'
               else
                  ! A GENUINE SOLVER FAILURE (info = 1: no convergence at
                  ! all). It returned its best iterate. Do NOT accept
                  ! an unconverged state as the answer -- go back to plain
                  ! time-marching. The failure is a property of the STATE the
                  ! hand-off was given, not of the run, so the solver is not
                  ! disabled for the rest of the run: marching moves the wind,
                  ! and after N_stall further steps the hand-off is offered the
                  ! state it has reached by then. Only after n_newton_attempt_max
                  ! such attempts is the finish given up, because by then the
                  ! marching is not producing states the solver can use.
                  newton_attempts = newton_attempts + 1
                  newton_fail_step = count
                  use_newton_solver = (newton_attempts .lt. n_newton_attempt_max)
                  if (use_newton_solver) then
                     write(*,'(A,I0,A,I0,A,I0,A)') ' (EXHALE_main) JFNK '//  &
                          'failed (info=', j, '); attempt ', newton_attempts, &
                          ' of ', n_newton_attempt_max,                       &
                          ' -- marching on and retrying'
                  else
                     write(*,'(A,I0,A,I0,A)') ' (EXHALE_main) JFNK failed '// &
                          '(info=', j, ') on all ', n_newton_attempt_max,     &
                          ' attempts; resuming time-marching with du-based'// &
                          ' stops'
                  endif
               endif
               endif   ! sec-ion pending / hold / engage
            endif

            ! Write to standard output (4th column: relative mass-flux level
            ! change over the last N_stall steps; -1 until the window fills)
            write(*,*) count,du,dtu,lev_rel
				
            !---------------------------------------------------!
            
            ! Detect NaNs
            do j = 1-Ng,N+Ng
            	do k = 1,3
            		dum = u(k,j)
            		if (dum.ne.dum) then
                              write(*,*)
            			write(*,'(A20,F8.6)') 'NaN detected at r = ',r(j)
                              write(*,'(A6,E13.6)') 'rho = ',W(1,j)*n0
                              write(*,'(A4,E13.6)') 'v = ',W(2,j)*v0/1.0e5
                              write(*,'(A4,E13.6)') 'p = ',W(3,j)*p0
                              write(*,'(A4,F8.1)')  'T = ',T(j)*T0
                              write(*,'(A6,E13.6)') 'nhi = ',nhi(j)*n0
                              write(*,'(A7,E13.6)') 'nhii = ',nhii(j)*n0
                              write(*,'(A7,E13.6)') 'nhei = ',nhei(j)*n0
                              write(*,'(A8,E13.6)') 'nheii = ',nheii(j)*n0
                              write(*,'(A9,E13.6)') 'nheiii = ',nheiii(j)*n0
                              write(*,'(A9,E13.6)') 'nheiTR = ',nheiTR(j)*n0
                              l_isnan = .true.
                        endif
            	enddo
         	enddo
            if(l_isnan) exit
            
            !---------------------------------------------------!
            
            !--- Write to files every 1000th iteration---! 
            if (mod(count,1000).eq.1) then
                  
                  ! Write thermodynamic and ionization profiles
                  ! The molecular columns are read from the sweep arrays: refresh
                  ! them from THIS state first, so one row is one state (the lower
                  ! ghosts' rho moves after the sweep; see the routine).
                  call molecular_carrier_densities_from_state(rho,f_sp)
                  call write_output(rho,v,p,T,heat,cool,eta,             &
                                    nhi,nhii,nhei,nheii,nheiii,          &
                                    nheiTR,nm,'eq')
                                   
            endif     
            
            ! Exit from temporal loop if only post processing has to be done
	      if (do_only_pp) exit

            ! Force continue for the first 1000 loops if force_start is enabled
            if (force_start) force_start = count .le. 1000

            ! Deterministic step cap for serial-vs-parallel verification.
            if (max_steps .gt. 0 .and. count .ge. max_steps) then
               hit_max_steps = .true.
               exit
            endif

            if (do_profile) then
               tp_tot = tp_tot + (omp_get_wtime() - tp_step0)
               if (mod(count, 500) .eq. 0 .and. tp_tot .gt. 0.0d0)          &
                  write(*,'(A,I7,A,F6.1,A)') ' (profile) count=', count,     &
                     '  ioniz_eq fraction=', 100.0d0*tp_ion/tp_tot, ' %'
            endif

      !---------------------------------------------------!

      ! End of temporal while loop
      enddo
      write(*,*) '(EXHALE_main.f90) Time integration done.'

      ! ONE STATE PER OUTPUT FILE. Everything written from here on -- the
      ! profiles, the heating and cooling breakdowns, the '# coupling:'
      ! header and the advection-corrected post-process -- describes the
      ! state the loop reached, under the physics the loop reached it under.
      ! The secondary-ionization coupling is therefore left exactly as the
      ! loop left it: a run whose staged activation never fired writes a
      ! state that was relaxed WITHOUT the coupling, and says so
      ! (write_coupling_state_header, sec_ion=F).
      !
      ! The alternative, arming the coupling here so that the rates "carry
      ! the full physics", produces two states in one set of files: the heat
      ! and cool columns of Hydro_ioniz.txt come from the marching sweeps
      ! and have the coupling off, while Heating_breakdown.txt and the _adv
      ! profiles are rebuilt with it on. Measured on the mol_base_handoff
      ! snapshot, the two names for the heating differed by a median 72 per
      ! cent. Neither number is wrong about its own state; the file is
      ! wrong about which state it holds.
      !
      ! It also makes the restart consistent: the header a run writes is
      ! read back at EXHALE_main.f90 about line 862, and a state relaxed
      ! without the coupling is one the staged activation must be allowed to
      ! stage in again, exactly as it would from any other cold start.
      ! (Development plan rev 3, section 10.2 item 2.)
      if (do_profile) then
         write(*,'(A)')          ' (profile) phase breakdown over the run:'
         write(*,'(A,F10.3,A)')  '   ioniz_eq total :', tp_ion, ' s'
         write(*,'(A,F10.3,A)')  '   full step total:', tp_tot, ' s'
         write(*,'(A,F10.3,A)')  '   hydro (3 RK stages):', tp_hydro, ' s'
         write(*,'(A,F10.3,A)')  '   carrier transport  :', tp_trans, ' s'
         write(*,'(A,F10.3,A)')  '   energy semi-impl.  :', tp_energy, ' s'
         write(*,'(A,F10.3,A)')  '   other (difference) :',                &
              tp_tot - tp_ion - tp_hydro - tp_trans - tp_energy, ' s'
         write(*,'(A)') '   eval_cool sub-blocks (all calls, all callers):'
         do jj_ec = 1, 5
            write(*,'(A,A,F10.3,A)') '     ', ec_name(jj_ec), ec_t(jj_ec), ' s'
         enddo
         if (tp_tot .gt. 0.0d0) write(*,'(A,F6.1,A)')                        &
            '   ioniz_eq fraction:', 100.0d0*tp_ion/tp_tot, ' %'
      endif

      ! Report which criterion stopped the loop. The tests below are the
      ! stop FLAGS themselves, in the order the loop leaves on them, so the
      ! line printed names the condition that actually fired -- never a
      ! different criterion that happens to share a variable. (Before this,
      ! a "Resid tol:" stop set the du flag and was printed as "momentum
      ! constant (du < du_th)" at du = 3.2e-3 with du_th = 1e-3.)
      if (l_isnan) then
         ! A NaN invalidates any convergence the flags may also carry, so it
         ! is tested first.
         write(*,'(A)') '     -> stopped: NaN in the conserved state'//      &
              ' (see the cell dump above)'
      else if (newton_finished .and. newton_returned_uncertified) then
         write(*,'(A,ES10.3,A)')                                              &
              '     -> stopped: the JFNK solve converged and the state it'//  &
              ' hands back does not meet ||R|| <', resid_max,                 &
              '; it is written UNCERTIFIED'
      else if (newton_finished) then
         ! resid_max is the target the finish was actually GIVEN: it is
         ! resid_th where the run set one and the code's own 1e-5 where it
         ! did not. Printing resid_th here reported -1 for every run that
         ! never set "Resid tol".
         write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
              '     -> converged: JFNK steady solution (||R|| <',             &
              resid_max, ' AND flux spread <', flux_spread_th, ')'
      else if (steady_gates_converged) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
              '     -> converged: both steady gates met while marching '//    &
              '(||R||=', resid_max, ' , flux spread=', fspread_gate, ')'
         write(*,'(A,ES10.3,A,ES10.3)')                                       &
              '        targets: "Resid tol"', resid_th,                       &
              ' , "Flux spread tol"', flux_spread_th
      else if (mass_flux_converged) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
              '     -> stopped: mass flux rho*v*r^2 constant over the '//     &
              'wind (du=', du, ' < du_th', du_th, ')'
         call report_marching_stop
      else if (step_change_converged) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
              '     -> stopped: the state stopped changing (dtu=', dtu,       &
              ' < dtu_th', dtu_th, ')'
         call report_marching_stop
      else if (du_plateaued) then
         write(*,'(A,I0,A)') '     -> stopped: du plateau detected '//       &
              '(stalled ', stall_count, ' steps)'
         call report_marching_stop
      else if (do_only_pp) then
         write(*,'(A)') '     -> stopped: "Do only PP" -- no time'//         &
              ' integration was requested'
      else if (hit_max_steps) then
         write(*,'(A,I0,A)') '     -> stopped: reached EXHALE_MAXSTEPS = ',  &
              max_steps, ' steps WITHOUT convergence'
         call report_marching_stop
      else if (count .ge. count_max) then
         write(*,'(A,I0,A)') '     -> stopped: reached count_max = ',          &
                             count_max, ' iterations WITHOUT convergence'
         call report_marching_stop
      endif
      write(*,'(A,I0,A,ES11.4,A,ES11.4)') '     final: count=', count,         &
                             '  du=', du, '  dtu=', dtu
      ! The value of the convergence functional the run actually used at its
      ! last step, with the window it was taken over, so that it can be
      ! recomputed from the state this run writes
      ! (src/tests/grid_and_gates/mass_flux_spread_functional.py).
      write(*,'(A,I0,A,I0,A,ES23.15)')                                        &
           '     mass_flux_spread: window cells ', j_min, '..', N,            &
           '  value=', du

      call write_run_counter_report

      !---------------------------------------------------!

      ! THE FINAL STATE IS THE ONE THE LAST ACCEPTED STEP PRODUCED, and no
      ! part of it is recomputed here (b1 T1.8: one consistent instantaneous
      ! evaluation feeds certification, the endpoint checks and every output
      ! file).  There used to be an excited_H_update on this line, refreshing
      ! the H(n=2) populations, the Balmer proton source and heat_balmer from
      ! the written state.  It made the written state two states: the `heat`
      ! column had already been assembled with the previous H(n=2) closure
      ! and was not rebuilt, so the heating this run reports and the channel
      ! sum that Heating_breakdown.txt reconstructs on the same state were
      ! two different numbers, which is what
      ! src/tests/grid_and_gates/output_state_consistency.sh measures.
      !
      ! What is left is the documented lag of the closure itself (b1 T1.7:
      ! H(n=2) is a 2x2 statistical equilibrium solved one outer pass before
      ! the state it closes), and that lag is a property of the accepted
      ! state, carried by every file that describes it.  Refreshing it after
      ! the adoption boundary does not remove the lag; it moves it into the
      ! gap between two files.

      !---------------------------------------------------!

      ! Write final thermodynamic and ionization profiles
      write(*,*) '(EXHALE_main.f90) Writing final results to file..'
      call element_census_reservoir('output (final state)', rho, f_sp)

      ! DID THIS RUN CLAIM A STATIONARY STATE AT ALL? Three stops declare
      ! one: the du threshold (the flux functional flat to du_th), the
      ! residual gate met while marching, and the steady finish returning
      ! info = 0. A stall, a du plateau, a step cap, a NaN and "Do only PP"
      ! declare nothing, and the state such a run writes is a relaxation
      ! snapshot (A0 run-mode contract section 2): it is written certified=F
      ! with that as the reason, and the run exits 0, because there was no
      ! claim to refuse.
      call certification_note_stationarity_claim(newton_finished .or.     &
                     steady_gates_converged .or. mass_flux_converged)

      ! THE CERTIFICATION OF THE STATE THAT IS ABOUT TO BE WRITTEN. Every
      ! earlier certification was made on a state the run then marched away
      ! from; this one is made on the state the files carry, so the
      ! 'certified=' field of the '# coupling:' header is a statement about
      ! the file it stands in. The residual is assembled here from that same
      ! state, and heat and cool are the ones the state carries.
      call assemble_residual(u, n_tot + ne, heat, cool, Rres)
      call certification_evaluate(cert_context_stationary, u, Rres, f_sp,  &
               resid_th, n_cells_without_chemical_root(last_sweep%acc_n), &
               .true., cert_now)
      call certification_report_write(cert_now, 'final state, as written')
      ! The molecular columns are read from the sweep arrays: refresh
      ! them from THIS state first, so one row is one state (the lower
      ! ghosts' rho moves after the sweep; see the routine).
      call molecular_carrier_densities_from_state(rho,f_sp)
      call write_output(rho,v,p,T,heat,cool,eta,                         &
                        nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')

      call assert_written_state_is_the_accepted_one

      ! Diagnostic: radiative cooling in each channel vs radius
      ! (reuses eval_cool's coefficients; see Huang et al. 2023 Fig. 10).
      call write_cool_breakdown_eq(T,rho,f_sp)

      ! Diagnostic: volumetric heating in each channel vs radius
      ! (photoionization per absorber + Balmer + He-recomb + Penning).
      call write_heat_breakdown_eq(T,rho,f_sp)

      ! Diagnostic: H(n=2) populations, Balmer proton source, and
      ! photoelectric/de-excitation heating vs radius (Huang Figs. 11/27/10/26).
      if (use_excited_H) call write_excited_H

      ! Rewrite the resolved-configuration record now that the wind exists:
      ! the configuration in it is the same one written before the run, and
      ! the elemental fluxes over the overlap window only exist once a
      ! diffusion step has measured them.
      call write_resolved_config

      !---------------------------------------------------!

      ! Post processing to include advection
      write(*,*) '(EXHALE_main.f90) Starting the post processing routine..'

      call post_process_adv(rho,v,p,T,heat,cool,eta,    &
                            nhi,nhii,nhei,nheii,nheiii,nheiTR,nm)
      
      write(*,*) '(EXHALE_main.f90) Post processing routine done.'

      ! Cells at which the ground singlet n(1^1S) = n(He I) - n(2^3S) came out
      ! negative and was floored at zero (he_ground_singlet_density in
      ! composition.f90). Reported here rather than with the counters above
      ! because the post-process is one of the consumers, so this is the first
      ! point at which the count is complete. Silent for a run whose two
      ! helium columns never crossed -- which is every run in which the
      ! metastable stays well below the summed He I, and the statement that
      ! such a run is the one an unguarded build produces.
      if (n_cells_he_singlet_clamped .gt. 0) then
         write(*,'(A,I0,A)')                                                  &
            '     helium ground singlet: ', n_cells_he_singlet_clamped,       &
            ' cell evaluation(s) floored at zero'
      endif
      
      !---------------------------------------------------!                            
                                    
      ! Get CPU time
      finish = omp_get_wtime()
      
      ! Write final execution time
      exec_time = finish - start
      n_hrs = floor(exec_time/3600.0)
      n_min = floor((exec_time - 3600.0*n_hrs)/60.0)
      n_sec = exec_time - 3600.0*n_hrs - 60.0*n_min

      write(*,*) ' '
      write(*,100) 'Execution Time = ', n_hrs,' h ', &
                                        n_min,' m ', &
                                        n_sec,' s'
100   format  (A17,I2,A3,I2,A3,F8.5,A2) 
      
      !---------------------------------------------------!
      
      ! Choose index to evaluate Mdot far enough from the top boundary.
      ! The offset is 20 cells; on a grid with fewer than 21 physical cells
      ! that index would land in the base ghosts, where the boundary closure
      ! and not the wind sets rho and v, so it is held inside [1,N].
      j = max(N - 20, 1)
      
      ! Evaluate steady state Mdot = 4 pi rho v r^2 [g/s]. A logarithm exists
      ! only where the flow at r(j) is directed outward: rho v <= 0 there is
      ! an inflow, which is not a mass-loss rate, and log10 of it is not a
      ! number to report.
      Mdot_cgs = 4.0*pi*rho(j)*v(j)*r(j)*r(j)*n0*v0*mu*R0*R0
      outflow_at_mdot_cell = (Mdot_cgs .gt. 0.0d0)
      if (outflow_at_mdot_cell) then
         Mdot = log10(Mdot_cgs)
      
         ! Correct for the 2D approximation used
         if (appx_mth.eq.'Rate/2 + Mdot/2') Mdot = Mdot - log10(2.0)
         if (appx_mth.eq.'Mdot/4') Mdot = Mdot - log10(4.0)
      else
         Mdot = 0.0d0
      endif
      
      
      ! Write Mdot in output
      write(*,*) ' '
      write(*,*) '----- Results -----'
      write(*,*) ' '
      write(*,*) '---> 2D approximate method: ', appx_mth
      if (outflow_at_mdot_cell) then
         write(*,101) ' ---> Log10 of steady-state Mdot = ', Mdot, ' g/s'
      else
         write(*,'(A,ES11.3,A)') ' ---> no outward mass flux at the '//      &
              'evaluation radius: rho*v*r^2 = ', Mdot_cgs,                   &
              ' g/s, no steady-state Mdot'
      endif
      
      ! Write Mdot to report file 
      open(unit = outfile, file = 'EXHALE_setup.out', access = 'append' )
      	write(outfile,101) ' '
      	if (outflow_at_mdot_cell) then
      	   write(outfile,102) ' - Log10 of steady-state Mdot = ', Mdot, ' g/s'
      	else
      	   write(outfile,'(A,ES11.3,A)') ' - no outward mass flux at the '// &
      	        'evaluation radius: rho*v*r^2 = ', Mdot_cgs, ' g/s'
      	endif
      close(unit = outfile)
      
101   format (A35,F5.2,A4)
102   format (A32,F5.2,A4)

      !---------------------------------------------------!

      ! THE RUN'S EXIT STATUS SAYS WHETHER WHAT IT WROTE IS CERTIFIED.
      ! Everything is on disk by now -- the profiles, the breakdowns, the
      ! advection-corrected post-process, the resolved-configuration record
      ! and the mass-loss rate above -- so the status is a verdict about a
      ! complete set of files and never a reason for an incomplete one.
      !
      ! THE STATUS FOLLOWS THE CLAIM, NOT THE RUN MODE. A run that DECLARED a
      ! stationary state -- the du stop, the marching residual gate, or the
      ! JFNK finish -- and then wrote a state the certification refused exits
      ! 2 and names the refusing entries, in EVERY run mode. A run that ended
      ! on a step cap, a stall, a du plateau or any other bound made no such
      ! claim and exits 0, in every run mode.
      ! The `init` exemption of the run-mode contract (section 2) covers the
      ! PHYSICAL-step claims -- the clock, the budgets, the histories -- which
      ! an initialization run does not make. Stationarity is not one of them:
      ! a du-stopped initialization run asserts that its written state is
      ! stationary, and the status is what tells a caller whether that
      ! assertion was certified. Which stops count as a claim is recorded by
      ! certification_note_stationarity_claim above, and the routine below
      ! exits 0 by itself for a run that claimed nothing.
      call certification_stop_uncertified

      !---------------------------------------------------!

      contains

      ! ------------------------------------------------------!

      real*8 function mass_flux_spread(mom_in)
      ! THE RADIAL SPREAD OF THE MASS FLUX rho v r^2 over the convergence
      ! window r >= "Escape radius [R_p]" (cells j_min..N):
      !
      !     spread = ( max_j F_j - min_j F_j ) / | mean_j F_j | ,
      !     F_j    = rho_j v_j r_j^2 .
      !
      ! A steady spherical wind carries one mass flux at every radius, so
      ! this is the statement that the marching has reached one: ATES
      ! (Caldiroli et al. 2021) stops on the constancy of rho v r^2 and
      ! CETIMB (Koskinen et al. 2013a) asks that F_c = rho v r^2 be constant
      ! with altitude. It is evaluated once per step and read by the stop,
      ! by the arming of the du triggers, by the secondary-ionization stage
      ! flip, by the JFNK hand-off and by the flux spread the run reports,
      ! so all of them speak about one number.
      !
      ! THE SIGN IS KEPT AND THE NORMALIZATION IS THE MEAN. Taking |F|
      ! before the extrema reads a window whose flux reverses -- an inflow
      ! at one radius and an outflow at another, which is not a steady wind
      ! at all -- as flat, and normalizing by the minimum measures the
      ! spread against the least representative point of the window instead
      ! of its level. The floor on the denominator keeps a window whose mean
      ! flux vanishes (equal inflow and outflow) from returning Inf; that
      ! state is as far from steady as the window can be, and the floor
      ! makes the measure large rather than undefined.
      ! (Development plan rev 3, section 10.2 item 5, decision 15.)
      real*8, dimension(1-Ng:N+Ng), intent(in) :: mom_in
      mass_flux_spread = (maxval(mom_in(j_min:N))                         &
                          - minval(mom_in(j_min:N)))                      &
           /max(abs(sum(mom_in(j_min:N))/dble(N - j_min + 1)), 1.0d-30)
      end function mass_flux_spread

      ! ------------------------------------------------------!

      subroutine write_operator_difference(step_no, du_now)
      ! FIXED-STATE OPERATOR DIFFERENCE AT THE PLM -> WENO3 HAND-OFF.
      !
      ! The two-stage recipe replaces one discretization by another between
      ! two time steps. What that costs is a property of the STATE and the
      ! two operators alone, and it can be measured before any step is taken
      ! with the new one: assemble the steady residual R = du/dt of the very
      ! same interior cells twice, once with each scheme -- each with its own
      ! ghost fill, its own reconstructed-face boundary states and its own
      ! momentum pressure split -- and write both, row by row and cell by
      ! cell, with the state, the local time step and the scale each row is
      ! measured against.
      !
      ! R_WENO3 - R_PLM is the instantaneous change the hand-off makes to the
      ! right-hand side; multiplied by the step it is the state change the
      ! first step after it inherits, which is what the transient that
      ! follows has to be compared with. Measured on the hot-Uranus molecular
      ! case (docs/f_plm_weno_continuation.md): dt |R|/|u| over the dtu window
      ! predicts the observed dtu of the last PLM step to 0.5 percent and of
      ! the first WENO3 step to 8 percent, so the transient there IS the
      ! operator difference and not an instability of WENO3.
      integer, intent(in) :: step_no
      real*8,  intent(in) :: du_now
      integer :: jj, uu
      logical :: sav_plm, sav_weno, sav_lam_on
      character(len=:), allocatable :: sav_method
      real*8  :: sc1, sc2, sc3

      sav_plm    = use_plm
      sav_weno   = use_weno3
      sav_method = rec_method
      sav_lam_on = recon_lambda_on
      recon_lambda_on = .false.       ! measure the two pure operators

      rec_method = 'PLM';    use_plm = .true.;  use_weno3 = .false.
      u_probe = u
      call Apply_BC(u_probe)
      call assemble_residual(u_probe, n_tot + ne, heat, cool, R_plm)

      rec_method = 'WENO3';  use_plm = .false.; use_weno3 = .true.
      u_probe = u
      call Apply_BC(u_probe)
      call assemble_residual(u_probe, n_tot + ne, heat, cool, R_weno)

      rec_method      = sav_method
      use_plm         = sav_plm
      use_weno3       = sav_weno
      recon_lambda_on = sav_lam_on

      open(newunit=uu, file='output/opdiff_profile.txt',                 &
           status='replace', action='write')
      write(uu,'(A,I0,A,ES12.4)') '# PLM -> WENO3 operator difference at '//&
           'step ', step_no, ', du = ', du_now
      write(uu,'(A)') '# same interior state, each operator with its own '//&
           'ghost fill and pressure split'
      write(uu,'(A,ES16.8)') '# dt = ', dt
      write(uu,'(A)') '# columns: r[Rp] n[cm-3] v[cm/s] T[K] '//          &
           'R_PLM_mass R_PLM_mom R_PLM_ene '//                            &
           'R_WENO3_mass R_WENO3_mom R_WENO3_ene '//                      &
           'scale_mass scale_mom scale_ene u_mass u_mom u_ene dt_loc'
      do jj = 1, N
         sc1 = residual_row_scale(1, jj, u)
         sc2 = residual_row_scale(2, jj, u)
         sc3 = residual_row_scale(3, jj, u)
         write(uu,'(1X,17(ES16.8,1X))') r(jj), W(1,jj)*n0, W(2,jj)*v0,   &
              T(jj)*T0,                                                  &
              R_plm(1,jj), R_plm(2,jj), R_plm(3,jj),                     &
              R_weno(1,jj), R_weno(2,jj), R_weno(3,jj), sc1, sc2, sc3,   &
              u(1,jj), u(2,jj), u(3,jj), dt_loc(jj)
      enddo
      close(uu)
      write(*,'(A,I0)') ' (EXHALE_main) wrote output/opdiff_profile.txt '//&
           'at step ', step_no

      end subroutine write_operator_difference

      ! ------------------------------------------------------!

      subroutine write_run_counter_report
      ! Everything a finished run has to say about ITSELF: how often each
      ! guard, limiter and solver branch was taken over the whole run.
      ! It is one routine because there is one such report, and BOTH routes
      ! that end a run call it -- the marching loop below its stop line, and
      ! the direct steady route (EXHALE_PTC=1) before it stops. The direct
      ! route used to stop without printing any of it, so a PTC solve
      ! reported no ionization-solver ledger, no positivity-limiter count and
      ! no base-inflow Mach at all.
      ! A ledger nothing was recorded in stays silent, so an ordinary
      ! marching run prints exactly what it printed before.
      ! Domain record of the one-sided CO destruction model, cumulative over
      ! the run: cell visits out of domain, above the shielding table's
      ! excitation-temperature limit, and with He+ + CO the leading He+
      ! loss; and the two worst ratios of the ordering with the radius of
      ! the first.
      integer :: co_dom_out, co_dom_hot, co_dom_hep
      real*8  :: co_dom_ratio, co_dom_r, co_dom_form

      ! THE THREE COUNTERS OF THE RUN MODE (A0 contract section 5), which are
      ! not one number: the steps tried, the steps kept, and, for a physical
      ! integration, the elapsed time those kept steps add up to. The gap
      ! between the first two is what a run spends in bisections; it is not
      ! part of any trajectory.
      write(*,'(A,I0,A,I0,A)') '     steps: ', n_steps_accepted,          &
           ' accepted of ', n_steps_attempted, ' attempted'
      if (run_mode .eq. run_mode_phys) then
         write(*,'(A,ES14.7,A)') '     physical elapsed time: ', t_phys,  &
              ' s (sum of the accepted global dt)'
      else
         write(*,'(A)') '     initialization / continuation: no'//        &
              ' physical elapsed time is claimed'
      endif
      call attempted_step_report(n_dt_halvings)
      ! The rows the carrier solve left at the arithmetic round-off of
      ! their own full terms: accepted and flagged, informational, never a
      ! statement that a state is wrong (A1scale3; review 2 section 5.3).
      if (thereis_mol .and. carrier_transport) then
         call carrier_roundoff_limited_record(nro_last, nro_total, nro_sub)
         if (nro_total .gt. 0)                                            &
            write(*,'(A,I0,A,I0,A)') '     carrier rows accepted at'//    &
                 ' their own round-off: ', nro_total, ' over ', nro_sub,  &
                 ' substep(s) (informational, never decisive)'
         ! Intervals the carrier retry could not cover. An attempt is what
         ! happened, so this ledger is never rolled back; a marked history
         ! is a statement about the states the run passed through and is
         ! reported as such.
         call carrier_exhausted_record(nex_n, nex_first, nex_last,        &
                  nex_j, nex_ic, nex_ratio, nex_phys, nex_full)
         if (nex_n .gt. 0) then
            write(*,'(A,I0,A,I0,A,I0)') '     carrier intervals NOT'//    &
                 ' covered: ', nex_n, ', first at step ', nex_first,      &
                 ', last at step ', nex_last
            write(*,'(A,I0,A,I0,A,ES12.5)') '       worst row: cell ',    &
                 nex_j, ', carrier ', nex_ic, ', |res|/physical terms ',  &
                 nex_ratio
            write(*,'(A,ES12.5,A,ES12.5)') '       physical terms ',      &
                 nex_phys, ', full terms ', nex_full
            if (.not. carrier_history_certifiable())                      &
               write(*,'(A)') '       THE CARRIER HISTORY OF THIS RUN'//  &
                    ' IS NOT CERTIFIABLE: an interval was left uncovered'
         endif
      endif

      ! Task 2: report how often the analytic-Jacobian Newton handled the
      ! ionization-equilibrium solve vs. fell back to MINPACK hybrd1.
      if (nt_calls .gt. 0) then
         write(*,'(A,I0,A,I0,A,F6.2,A)')                                       &
            '     ioniz-eq solver: Newton ', nt_calls - nt_fallback,           &
            ' / ', nt_calls, ' solves (',                                      &
            100.0d0*dble(nt_calls-nt_fallback)/dble(nt_calls),                 &
            '%); rest fell back to hybrd1'
      endif

      ! Acceptance ledger of the equilibrium solves, printed once per state
      ! kind: the states the marching loop held, the steady solver's own
      ! iterates, and the states the steady solver only probed to build its
      ! Newton model. Root validation, the molecular clamp and hybrd1 exit
      ! codes, the acceptance classes of section 113 with their largest
      ! normalized reaction residual and the longest non-root streak, the cost
      ! of the constrained element-conserving solve, the residual-decade
      ! histograms, and the admissibility signal of section 121 all belong to
      ! one ledger each. A ledger nothing was recorded in stays silent, so an
      ! ordinary marching run prints exactly what it printed before.
      call write_ioniz_eq_acceptance_report


      ! Cells whose semi-implicit energy update asked for a temperature below
      ! the range where the equation of state and the chemical network are
      ! defined (energy_semi_implicit). Such a cell has NOT converged to a
      ! physical temperature: the energy balance has no root inside the
      ! bracket, and at 0.01 T0 almost any composition satisfies its reaction
      ! balance, so a state built there is constrained by nothing. No state is
      ! assembled from one: the counters below are ATTEMPTS. Each activation
      ! refuses the attempted step it happened in, so a run can carry many of
      ! them and still end on accepted states -- what they count is how often
      ! the controller had to shorten a step for this reason. Silent for a run
      ! whose energy update stayed inside its bracket.
      if (n_energy_floor_hits .gt. 0) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                                  &
            '     energy floor WARNING: ', n_energy_floor_hits,              &
            ' failed attempt(s) in ', n_energy_floor_cells(),                &
            ' distinct cell(s), steps ', energy_floor_first_step, ' to ',    &
            energy_floor_last_step, ''
      endif

      ! Cells whose Crank-Nicolson transport stage solved for a temperature
      ! below the same floor (viscous_conduction). The floor is the lower end
      ! of the range in which the equation of state and the transport
      ! coefficients are defined, and it is not a reservoir: returning it in
      ! place of the solved temperature would put c_v (T_floor - T_solved) per
      ! volume into the state that no term of the equation supplied. So no
      ! state is built there either, the counters below are ATTEMPTS, and in
      ! production the run stops at the first of them; more than one can only
      ! be reached by a test driver that suppresses the stop
      ! (conduction_stop_on_failure). Silent for a run whose transport stage
      ! stayed above the floor, and for every run without the
      ! "Viscosity:"/"Conduction:" keys, which never enters the stage.
      if (n_conduction_floor_hits .gt. 0) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                                  &
            '     conduction floor WARNING: ', n_conduction_floor_hits,      &
            ' failed attempt(s) in ', n_conduction_floor_cells(),            &
            ' distinct cell(s), steps ', conduction_floor_first_step,        &
            ' to ', conduction_floor_last_step, ''
      endif

      ! Faces at which the high-order reconstruction produced a non-positive
      ! density or pressure and was scaled back toward the cell average
      ! (positivity_limited_faces). Silent for a run that never needed it --
      ! which is every run whose solution the reconstruction resolves, and the
      ! statement that such a run is the one an unguarded build produces.
      if (n_faces_positivity_limited .gt. 0) then
         write(*,'(A,I0,A)')                                                  &
            '     reconstruction: ', n_faces_positivity_limited,              &
            ' face state(s) scaled toward the cell average for positivity'
      endif

      ! Boundary ghost cells whose linear extrapolation left rho > 0, p > 0
      ! and was dropped to the zero-gradient state (the outer free-outflow
      ! ghosts under WENO3). Silent
      ! for a run whose boundary extrapolations stayed admissible, which is
      ! the run an unguarded build produces.
      if (n_ghost_cells_positivity_limited .gt. 0) then
         write(*,'(A,I0,A)')                                                  &
            '     boundary: ', n_ghost_cells_positivity_limited,              &
            ' ghost cell state(s) limited to zero gradient for positivity'
      endif

      ! The Mach number of the inflow the lower boundary admitted. Always
      ! reported for a run that had a base inflow at all: the maximum is a
      ! statement about how close the boundary came to over-specification,
      ! and it is worth having in the log of a run that never crossed as
      ! well as one that did.
      if (base_inflow_mach_max .gt. 0.0d0) then
         if (base_supersonic_inflow_first_step .ge. 0) then
            write(*,'(A,I0,A,I0,A,ES10.3)')                                   &
               '     base boundary WARNING: supersonic inflow on ',           &
               n_base_supersonic_inflow_steps, ' step(s), first at ',         &
               base_supersonic_inflow_first_step, ', max Mach ',              &
               base_inflow_mach_max
            write(*,*) '       The interior downstream of a supersonic'//     &
                       ' inflow is not a solution of'
            write(*,*) '       the stated problem; see section 117 of'//      &
                       ' docs/Update_EXHALE.'
         else
            write(*,'(A,ES10.3)')                                             &
               '     base boundary: inflow stayed subsonic, max Mach ',       &
               base_inflow_mach_max
         endif
      endif

      ! Cells in which the He recombination photons ionized H I faster than
      ! the stellar radiation field did, by three decades or more. The
      ! coupling is a correction to the direct photoionization integral, so
      ! this says the on-the-spot photon budget of those cells has run away.
      ! Silent for a run in which it never happened, which is every run whose
      ! cells keep an absorber for those photons.
      if (n_cells_he_rec_photoionization_dominant .gt. 0) then
         write(*,'(A,I0,A,ES10.2,A,ES8.1,A)')                                 &
            '     He recombination coupling WARNING: ',                      &
            n_cells_he_rec_photoionization_dominant,                         &
            ' cell-step(s) with dP_HI(He rec)/P_HI up to ',                  &
            he_rec_photoionization_ratio_max, ' (threshold ',                &
            he_rec_dominant_ratio, ')'
      endif

      ! The first-order flux correction that holds an RK stage inside
      ! rho > 0, rho e > 0 (positivity_limited_fluxes). Printed even when
      ! every count is zero: zero says that no stage of the run ever left the
      ! admissible set, so no interface carried a first-order flux and the
      ! trajectory is the one an unguarded build produces, and that is a
      ! statement about the run worth reading. The three numbers are
      ! different quantities:
      !   calls      = stages offered the repair, successful or not;
      !   faces      = interfaces repaired in the calls that succeeded;
      !   in accepted steps = those of them belonging to steps the marching
      !                loop kept, the repairs of an attempt discarded by the
      !                dt bisection being counted in "faces" but not here.
      write(*,'(A,I0,A,I0,A,I0,A)')                                           &
         '     flux correction: ', n_calls_flux_positivity_repair,            &
         ' positivity repair call(s), ', n_faces_flux_positivity_limited,     &
         ' face flux(es) dropped to first order, ',                           &
         n_faces_flux_positivity_limited_accepted, ' of them in accepted steps'

      ! Faces evaluated with the Lax-Friedrichs flux because the run selected
      ! it. An LLF run prints it even when it is zero, for the same reason the
      ! Roe branch below prints its counter; no other flux reports it, the
      ! repair faces above being counted separately.
      if (flux .eq. 'LLF') then
         write(*,'(A,I0)')                                                    &
            '     LLF faces: ', n_faces_llf
      endif

      ! Faces at which the Roe flux was replaced by the HLLE flux because
      ! the estimated star state does not exist (separating flow leaves
      ! vacuum) or is not admissible (docs/a2_roe_interface.md section 3).
      ! A Roe run prints it even when it is zero: zero is the statement that
      ! every face kept the Roe flux, and that is worth reading. The counter
      ! belongs to the Roe branch, so no other flux reports it.
      if (flux .eq. 'ROE') then
         write(*,'(A,I0)')                                                    &
            '     ROE faces on the HLLE branch: ', n_faces_roe_hlle
      endif

      ! Steps retaken at half dt because an RK stage left rho > 0, rho e > 0.
      ! Silent for a run that stayed inside the positivity bound of its CFL
      ! number, which is the run an unguarded build produces.
      if (n_steps_dt_halved .gt. 0) then
         write(*,'(A,I0,A,I0,A)')                                             &
            '     positivity step control: ', n_steps_dt_halved,              &
            ' step(s) retaken, ', n_dt_halvings, ' dt bisection(s)'
      endif

      ! The artificial stress is a numerical dissipation, so it is admissible
      ! only where it is negligible against the physical fluxes. Record on the
      ! final state how big it actually was and where its gate was still open.
      if (low_mach_damping_active()) then
         call contact_mode_dissipation_magnitude(u, lowmach_ratio,           &
                                     lowmach_r_peak, lowmach_r_gate)
         if (lowmach_r_gate .gt. 0.0d0) then
            write(*,'(A,ES9.2,A,F7.4,A,F7.4,A)')                             &
               '     low-Mach damping: peak |D_p|/|rho v^2+p| =',            &
               lowmach_ratio, ' at r =', lowmach_r_peak,                     &
               ' R_p; Mach gate open out to r =', lowmach_r_gate, ' R_p'
         else
            write(*,'(A)') '     low-Mach damping: Mach gate closed'//       &
               ' everywhere on the final state (term identically zero)'
         endif
      endif

      ! THE DOMAIN OF THE ONE-SIDED CO DESTRUCTION MODEL, OVER THE WHOLE
      ! RUN.  The CO row destroys CO and never forms it; that is legitimate
      ! in a cell where tau_dest << tau_res << tau_form, and both
      ! inequalities are measured on the state each accepted carrier step
      ! hands on.  The record is informational: the rates are on
      ! everywhere, and where the ordering fails the omitted formation
      ! fails with it, so the transported value stands.  The line is
      ! printed even when nothing was out of domain, because a zero in a
      ! cumulative record IS the statement that the model stayed in domain.
      if (carrier_transport .and. thereis_oxychem) then
         call carrier_co_domain_record(co_dom_out, co_dom_hot, co_dom_hep,&
                                       co_dom_ratio, co_dom_r, co_dom_form)
         write(*,'(A,ES9.2,A,I0,A,ES11.4,A,F8.4,A,ES11.4)')               &
            '     CO destruction domain: f_dom =',                       &
            carrier_co_domain_f_dom(), ', ', co_dom_out,                 &
            ' cell visits with tau_dest > f_dom tau_res, worst ratio ',  &
            co_dom_ratio, ' at r =', co_dom_r,                           &
            ' R_p; worst tau_res/tau_form ', co_dom_form
         write(*,'(A,I0,A,I0,A)')                                        &
            '       of the same visits, ', co_dom_hot,                   &
            ' above the 512 K excitation-temperature limit of the CO'//  &
            ' shielding table and ', co_dom_hep,                         &
            ' with He+ + CO the leading He+ loss of the cell'
      endif

      ! THE COUPLED SOURCE STEP, OVER THE WHOLE RUN (b1 T1.4 to T1.9, T2.1).
      ! Three attempt statistics and one measured invariant. The pass count
      ! is the cost of the coupling; the cap count is the number of steps
      ! whose composition and temperature did not reach their fixed point,
      ! and it is a property of attempts, not of the adopted state. The mass
      ! residual is T2.2: the species mass sum against the density the step
      ! was given, which the chemistry may not write. Under the declared
      ! convention (the electron mass is carried with its ion) the bound
      ! implied by what is left out is n_e m_e / rho, at most m_e/m_H = 5.4e-4
      ! in a fully ionized hydrogen gas and smaller everywhere else.
      if (count .gt. 0) then
         write(*,'(A,F6.2,A,I0,A,I0,A,I0,A)')                            &
            '     coupled source step: ', dble(n_csm_passes_total)/       &
            dble(max(n_csm_calls,1)), ' passes per coupled step on'//     &
            ' average over ', n_csm_calls, ' of them, worst ',            &
            n_csm_passes_worst, ' (step ', n_csm_worst_step, ')'
         if (n_csm_iter_cap .gt. 0) then
            write(*,'(A,I0,A,I0,A)')                                     &
               '     coupled source step: ', n_csm_iter_cap,             &
               ' step(s) reached the pass cap of ', csm_max_pass,        &
               ' without a fixed point'
         else
            write(*,'(A)') '     coupled source step: every step reached'&
               //' its fixed point within the pass cap'
         endif
         ! WHICH BOUND THE EXITS WERE TAKEN ON.  The stopping test is a
         ! statement about the error of the state it accepts: on a pass
         ! whose sequence gives an estimate the bound is
         ! |theta/(1-theta)| times the increment, and on a pass that gives
         ! none it is the increment itself, which is the stricter of the
         ! two wherever the guards admit a decaying mode.
         write(*,'(A,I0,A,I0,A)')                                        &
            '     coupled source step: ', n_csm_exit_est,                &
            ' exit(s) on the estimated error, ', n_csm_exit_incr,        &
            ' on the last increment'
         ! WHAT THE ACCEPTED STATE IS WORTH, MEASURED against the fixed
         ! point of the same sequence (csm_err_probe).  The two worst
         ! distances are in the two norms of the stopping test and belong
         ! beside its two tolerances; the counts say on how many steps the
         ! accepted state stood further away than the tolerance it was
         ! accepted under.
         if (csm_err_probe) then
            write(*,'(A,I0,A,ES10.3,A,ES10.3)')                          &
               '     accepted state against its fixed point: ',          &
               n_csm_probe, ' measured, worst distance dT/T ',           &
               csm_probe_worst_T, ', worst dcomp ', csm_probe_worst_c
            write(*,'(A,I0,A,I0,A,I0,A)')                                &
               '     accepted state against its fixed point: ',          &
               n_csm_probe_over_T, ' over csm_T_tol, ',                  &
               n_csm_probe_over_c, ' over csm_comp_tol, ',               &
               n_csm_probe_stall,                                        &
               ' probe(s) short of csm_probe_tol at the pass cap'
         endif
         write(*,'(A,ES10.3,A,I0)')                                      &
            '     chemistry mass closure: worst |sum n_s m_s - rho|/rho =',&
            csm_mass_resid_run, ' at step ', csm_mass_resid_step
         ! The size of the formation/excitation reservoir change against the
         ! step's own thermal source. It is NOT a term of the energy row
         ! under this code's definition of heat and cool (see the coupled
         ! source step), and this is what it would be worth if they were
         ! redefined as the gross exchange with the radiation field.
         write(*,'(A,ES10.3)')                                            &
            '     formation reservoir change, worst |du_form| /'//        &
            ' |dt (heat-cool)| =', csm_uform_frac_run
         ! The extrapolation of the pass sequence's geometric limit. The
         ! last line is the admissibility assertion: both records are
         ! identically zero when the ratio test that builds the candidate
         ! is right, so a nonzero value says an inadmissible state was
         ! handed to the next pass.
         if (csm_extrap_on) then
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)')                  &
               '     geometric limit of the pass sequence: ',             &
               n_csm_x_attempt, ' attempt(s), ', n_csm_x_guard,           &
               ' refused by the guards, ', n_csm_x_applied, ' taken, ',   &
               n_csm_x_kept, ' beat the pass they replaced, ',            &
               n_csm_x_slack, ' fell short of it and ', n_csm_x_refused,  &
               ' undone'
            write(*,'(A,I0,A,I0,A,F6.3)')                                 &
               '     geometric limit: ', n_csm_x_exit, ' of the ',        &
               n_csm_x_applied, ' taken were met by the stopping test'//  &
               ' on the very next pass; smallest damping of an'//         &
               ' extrapolation step ', csm_x_damp_min_seen
            write(*,'(A,ES10.3,A,ES10.3)')                                &
               '     geometric limit admissibility: worst species'//      &
               ' fraction outside [0,1] ', csm_x_worst_fneg,              &
               ', worst temperature outside its bracket ',                &
               csm_x_worst_T_out
         endif
      endif

      end subroutine write_run_counter_report

      ! ------------------------------------------------------!

      subroutine coupled_pair_geometric_error_estimate
      ! THE ERROR THE COUPLED PAIR STILL CARRIES, FROM THREE TERMS OF ITS
      ! OWN PASS SEQUENCE.  Called on every pass, after the pass's two
      ! increments have been measured and before the stopping test is
      ! taken.  It is the ONE place the sequence's ratio is formed: the
      ! stopping test reads the error it returns, and the extrapolation
      ! reads the same ratio to place its candidate at the limit.
      !
      ! A sequence of ratio theta stands theta/(1 - theta) times its last
      ! increment from its limit, so
      !     csm_err = csm_err_safety * |theta/(1 - theta)|
      !               * (this pass's increment)
      ! in each of the two norms, the factor carrying the L2 quantity
      ! theta is taken from into the max norm the test is taken in (see
      ! csm_err_safety, where the measurement that sets it is written).
      ! theta is the Rayleigh quotient of the whole increment VECTOR, one
      ! scalar for the temperature block and one for the composition
      ! block; see the constants' declaration for why a ratio taken cell
      ! by cell cannot work and was measured not to.
      !
      ! WHERE THE SEQUENCE GIVES NO ESTIMATE the increment itself is
      ! returned and csm_err_ok is false.  That is not a guess at the
      ! error: it is the only bound the loop can state from what it has,
      ! and it is the stricter of the two on this loop, since
      ! csm_err_safety*|theta/(1 - theta)| is 0.2 to 0.7 at the |theta|
      ! the guards admit here.
      real*8  :: p11, p12, p22, q11, q12, q22, mT, mc
      integer :: j, isp
      logical :: geometric

      ! The increments of THIS pass, in the stopping test's own units:
      ! relative for the temperature, absolute mass fraction for the
      ! composition.  Over 1:N, because that is the range the test is taken
      ! over and the norm has to be the test's norm.
      do j = 1, N
         csm_geom_eT(j) = (T(j) - T_csm_prev(j))/max(abs(T(j)), 1.0d-99)
      enddo
      do isp = 1, n_species
         do j = 1, N
            csm_geom_ec(j,isp) = f_sp(j,isp) - f_sp_csm(j,isp)
         enddo
      enddo

      csm_err_ok       = .false.
      csm_err_T        = csm_dT
      csm_err_c        = csm_dcomp
      csm_geom_theta_T = 0.0d0
      csm_geom_theta_c = 0.0d0
      if (csm_geom_refuse_all) return
      if (.not. csm_geom_have_prev) return
      if (i_csm .lt. csm_geom_first_pass) return

      p22 = sum(csm_geom_eT_prev(1:N)**2)
      p12 = sum(csm_geom_eT(1:N)*csm_geom_eT_prev(1:N))
      p11 = sum(csm_geom_eT(1:N)**2)
      q22 = sum(csm_geom_ec_prev(1:N,:)**2)
      q12 = sum(csm_geom_ec(1:N,:)*csm_geom_ec_prev(1:N,:))
      q11 = sum(csm_geom_ec(1:N,:)**2)

      geometric = (p22 .gt. 0.0d0) .and. (q22 .gt. 0.0d0)
      if (geometric) then
         csm_geom_theta_T = p12/p22
         csm_geom_theta_c = q12/q22
         ! A decaying mode of EITHER sign whose amplification is bounded,
         ! and an increment that has actually shrunk.  The sign is not a
         ! free choice: MEASURED on mol_base_handoff, the cosine between
         ! consecutive increments is -0.85 to -1.00, so the mode of this
         ! loop ALTERNATES and its ratio is negative (-0.06 to -0.16 in
         ! the L2 norm of the pair).  That is what an alternating solve of
         ! two blocks does, each of the two overshooting what the other
         ! asks for; the magnitudes COST5 measured (0.072 and 0.171) are
         ! |theta| and carry no sign.  The geometric sum holds for any
         ! |theta| < 1, and for theta negative it places the limit BETWEEN
         ! the last two iterates, at |theta|/(1+|theta|) of the way back.
         geometric = abs(csm_geom_theta_T) .gt. 0.0d0                    &
               .and. abs(csm_geom_theta_T) .lt. csm_geom_theta_max       &
               .and. abs(csm_geom_theta_c) .gt. 0.0d0                    &
               .and. abs(csm_geom_theta_c) .lt. csm_geom_theta_max       &
               .and. p11 .lt. p22 .and. q11 .lt. q22
      endif
      ! WHAT THE SEQUENCE LOOKS LIKE, for the measurement the guards rest
      ! on.  The magnitude ratio is || d_k || / || d_{k-1} ||, the cosine
      ! is the angle between the two increments, and the model residual
      ! below is what a single decaying mode leaves unexplained.  A
      ! magnitude ratio that decays while the cosine is small is a
      ! sequence whose ERROR NORM contracts geometrically without its
      ! error DIRECTION persisting, and no single-mode statement acts on
      ! it.
      if (csm_geom_debug .and. p22 .gt. 0.0d0 .and. q22 .gt. 0.0d0)     &
         write(*,'(a,i0,a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')             &
            ' CSMXSEQ step ', count, ' pass ', i_csm, ' ratioT= ',      &
            sqrt(p11/p22), ' cosT= ', p12/sqrt(max(p11*p22,1.0d-99)),   &
            ' ratioC= ', sqrt(q11/q22), ' cosC= ',                      &
            q12/sqrt(max(q11*q22,1.0d-99))
      if (geometric) then
         ! The geometric model's own residual, || d_k - theta d_{k-1} ||
         ! against || d_k ||.  This is the question a scalar ratio taken
         ! component by component cannot be asked, and for a
         ! Rayleigh-quotient theta it is exactly |sin| of the angle
         ! between the two increments.
         mT = max(p11 - 2.0d0*csm_geom_theta_T*p12                       &
                      + csm_geom_theta_T**2*p22, 0.0d0)
         mc = max(q11 - 2.0d0*csm_geom_theta_c*q12                       &
                      + csm_geom_theta_c**2*q22, 0.0d0)
         geometric = mT .le. csm_geom_model_tol**2*p11                    &
               .and. mc .le. csm_geom_model_tol**2*q11
      endif

      if (geometric) then
         csm_err_ok = .true.
         ! csm_err_safety carries the norm the estimate is taken in over
         ! to the norm the test is taken in; see its declaration for the
         ! measurement that sets it.
         csm_err_T  = csm_err_safety                                      &
                      *abs(csm_geom_theta_T/(1.0d0 - csm_geom_theta_T))   &
                      *csm_dT
         csm_err_c  = csm_err_safety                                      &
                      *abs(csm_geom_theta_c/(1.0d0 - csm_geom_theta_c))   &
                      *csm_dcomp
      endif

      end subroutine coupled_pair_geometric_error_estimate

      ! ------------------------------------------------------!

      subroutine coupled_pair_fixed_point_distance
      ! HOW FAR THE ACCEPTED STATE STOOD FROM THE FIXED POINT, MEASURED.
      ! Called on the probe pass whose increment reached csm_probe_tol,
      ! four decades below the stopping tolerance, so the pair the loop
      ! stands at is the fixed point at this precision and the distance to
      ! the state the stopping test accepted is the error that test
      ! admitted.  Both norms are the test's own: relative max for the
      ! temperature, absolute mass fraction max for the composition.  One
      ! line per accepted state, so that the whole distribution can be read
      ! from a run and not just its worst entry.
      real*8  :: dist_T, dist_c
      integer :: j, isp

      dist_T = 0.0d0
      do j = 1, N
         dist_T = max(dist_T, abs(T(j) - csm_probe_T_acc(j))              &
                              /max(abs(T(j)), 1.0d-99))
      enddo
      dist_c = 0.0d0
      do isp = 1, n_species
         do j = 1, N
            dist_c = max(dist_c, abs(f_sp(j,isp) - csm_probe_f_acc(j,isp)))
         enddo
      enddo

      n_csm_probe = n_csm_probe + 1
      if (dist_T .gt. csm_T_tol)    n_csm_probe_over_T = n_csm_probe_over_T + 1
      if (dist_c .gt. csm_comp_tol) n_csm_probe_over_c = n_csm_probe_over_c + 1
      csm_probe_worst_T = max(csm_probe_worst_T, dist_T)
      csm_probe_worst_c = max(csm_probe_worst_c, dist_c)
      write(*,'(a,i0,a,i0,a,i0,a,es11.4,a,es11.4,a,es11.4,a,es11.4)')     &
         ' CSMERR step ', count, ' accepted at pass ',                    &
         csm_probe_pass_accepted, ' fixed point at pass ', i_csm,         &
         ' distT= ', dist_T, ' distC= ', dist_c, ' dT= ', csm_dT,        &
         ' dC= ', csm_dcomp

      end subroutine coupled_pair_fixed_point_distance

      ! ------------------------------------------------------!

      subroutine coupled_pair_geometric_limit
      ! THE LIMIT OF THE COUPLED PASS SEQUENCE AS THE NEXT ITERATE, from
      ! the ratio coupled_pair_geometric_error_estimate has just measured.
      ! Called at the end of a pass that did not meet the stopping test.
      ! Where the estimate exists and the error it reports is already
      ! within reach of the tolerance, the pair is moved to the limit of
      ! the geometric sequence it belongs to, damped to the largest
      ! multiple that keeps the state admissible.  Taking the move breaks
      ! the sequence deliberately, and the loop then collects a fresh pair
      ! before another estimate can be taken.
      real*8  :: aT, ac, sdamp, bnd, stp, gj, gstp
      real*8  :: t_floor_j, t_ceil_j, fx
      integer :: j, isp
      logical :: geometric

      if (csm_x_used .ge. csm_extrap_max) return
      ! No estimate, no limit: the two are the same statement about the
      ! same three terms.
      if (.not. csm_err_ok) return

      n_csm_x_attempt = n_csm_x_attempt + 1
      geometric = .not. csm_extrap_refuse_all

      if (geometric) then
         ! THE REACH OF THE MODEL, and this guard is what decides whether
         ! the extrapolation pays at all.  What the candidate removes is
         ! the modeled part of the error; what it LEAVES is the part the
         ! model does not explain, and MEASURED that part is of the same
         ! order as the modeled part itself unless the sequence is
         ! already close to its limit (see csm_x_reach).  So the candidate
         ! is taken only where the estimated error is already inside the
         ! ball the loop is aiming for, loosened by csm_x_reach: that is
         ! the stopping test itself at a looser tolerance, and it is why
         ! the window left to the extrapolation is now only the band
         ! between the tolerance and csm_x_reach times it.
         geometric = csm_err_T .le. csm_x_reach*csm_T_tol                 &
               .and. csm_err_c .le. csm_x_reach*csm_comp_tol
      endif

      if (geometric) then
         aT = csm_geom_theta_T/(1.0d0 - csm_geom_theta_T)
         ac = csm_geom_theta_c/(1.0d0 - csm_geom_theta_c)
         ! THE RATIO TEST.  sdamp is the largest multiple of the
         ! extrapolation step for which every constraint still holds.
         sdamp = 1.0d0
         do j = 1, N
            ! The temperature stays inside the bracket the energy solve
            ! solves in (energy_semi_implicit, T_floor / T_ceil).
            t_floor_j = max(T_floor_code_min, T_eos_floor_K/T0)
            t_ceil_j  = T_ceiling_atomic_K/T0
            if (caloric_mixture_active .and.                             &
                allocated(molecular_cell)) then
               if (molecular_cell(j)) t_ceil_j = T_ceiling_mol_K/T0
            endif
            t_ceil_j = max(t_ceil_j, t_floor_j)
            stp = aT*(T(j) - T_csm_prev(j))
            if (stp .lt. 0.0d0) then
               bnd = (t_floor_j - T(j))/stp
               sdamp = min(sdamp, max(bnd, 0.0d0))
            elseif (stp .gt. 0.0d0) then
               bnd = (t_ceil_j - T(j))/stp
               sdamp = min(sdamp, max(bnd, 0.0d0))
            endif
            ! Every species fraction stays in [0,1].
            do isp = 1, n_species
               stp = ac*csm_geom_ec(j,isp)
               if (stp .lt. 0.0d0) then
                  bnd = -f_sp(j,isp)/stp
                  sdamp = min(sdamp, max(bnd, 0.0d0))
               elseif (stp .gt. 0.0d0) then
                  bnd = (1.0d0 - f_sp(j,isp))/stp
                  sdamp = min(sdamp, max(bnd, 0.0d0))
               endif
            enddo
            ! The He 2^3S column is a SUB-POPULATION of the He I column
            ! (species_table, bsp_is_excited_level), so the singlet
            ! population is their difference and it has to stay
            ! non-negative.  An affine step keeps each column
            ! non-negative on its own but not their difference.
            if (thereis_HeITR) then
               gj   = f_sp(j,isp_HeI) - f_sp(j,isp_HeTR)
               gstp = ac*(csm_geom_ec(j,isp_HeI)                          &
                          - csm_geom_ec(j,isp_HeTR))
               if (gstp .lt. 0.0d0) then
                  bnd = -gj/gstp
                  sdamp = min(sdamp, max(bnd, 0.0d0))
               endif
            endif
         enddo
         sdamp = min(1.0d0, csm_x_backoff*sdamp)
         geometric = sdamp .ge. csm_x_damp_min
      endif

      if (geometric) then
         csm_x_damp = sdamp
         csm_x_damp_min_seen = min(csm_x_damp_min_seen, sdamp)
         ! The iterate as it stands, so that a candidate whose pass turns
         ! out worse than the pass it replaced can be undone.
         csm_x_T_save(1:N)   = T(1:N)
         csm_x_f_save(1:N,:) = f_sp(1:N,:)
         ! What the unextrapolated pass would have moved, from the two
         ! ratios just measured.  The verdict one pass later is against
         ! this and against nothing else.
         csm_x_pred_dT = abs(csm_geom_theta_T)*csm_dT
         csm_x_pred_dc = abs(csm_geom_theta_c)*csm_dcomp
         ! The increments the iterate itself stands at, which is what the
         ! extrapolated pass has to improve on for its state to be worth
         ! keeping at all.
         csm_x_at_dT = csm_dT
         csm_x_at_dc = csm_dcomp
         do j = 1, N
            T(j) = T(j) + sdamp*aT*(T(j) - T_csm_prev(j))
         enddo
         do isp = 1, n_species
            do j = 1, N
               f_sp(j,isp) = f_sp(j,isp) + sdamp*ac*csm_geom_ec(j,isp)
            enddo
         enddo
         ! THE ADMISSIBILITY ASSERTION.  Both records are identically zero
         ! if the ratio test above is right; a nonzero value is the
         ! statement that an inadmissible state was handed on, and the run
         ! reports it.
         do j = 1, N
            fx = minval(f_sp(j,1:n_species))
            csm_x_worst_fneg = max(csm_x_worst_fneg, -fx)
            fx = maxval(f_sp(j,1:n_species))
            csm_x_worst_fneg = max(csm_x_worst_fneg, fx - 1.0d0)
            if (thereis_HeITR) csm_x_worst_fneg = max(                    &
               csm_x_worst_fneg,                                          &
               f_sp(j,isp_HeTR) - f_sp(j,isp_HeI))
            t_floor_j = max(T_floor_code_min, T_eos_floor_K/T0)
            t_ceil_j  = T_ceiling_atomic_K/T0
            if (caloric_mixture_active .and.                              &
                allocated(molecular_cell)) then
               if (molecular_cell(j)) t_ceil_j = T_ceiling_mol_K/T0
            endif
            t_ceil_j = max(t_ceil_j, t_floor_j)
            csm_x_worst_T_out = max(csm_x_worst_T_out,                     &
                                    t_floor_j - T(j), T(j) - t_ceil_j)
         enddo
         n_csm_x_applied     = n_csm_x_applied + 1
         csm_x_used          = csm_x_used + 1
         csm_x_pending       = .true.
         csm_geom_seq_broken = .true.
         if (csm_geom_debug)                                              &
            write(*,'(a,i0,a,i0,a,f6.3,a,f6.3,a,f6.3,a,es10.3,a,'//       &
                  'es10.3)')                                              &
               ' CSMXTR step ', count, ' pass ', i_csm, ' thetaT= ',      &
               csm_geom_theta_T, ' thetaC= ', csm_geom_theta_c,           &
               ' damp= ', sdamp, ' predT= ', csm_x_pred_dT, ' predC= ',   &
               csm_x_pred_dc
      else
         n_csm_x_guard = n_csm_x_guard + 1
      endif

      end subroutine coupled_pair_geometric_limit

      ! ------------------------------------------------------!

      subroutine coupled_pair_extrapolation_verdict
      ! WHETHER THE EXTRAPOLATED STATE WAS WORTH THE PASS, read from the
      ! increment the pass taken FROM it just measured.  The pass the plain
      ! iteration would have taken instead would have moved the pair by
      ! theta times the previous increment, in each of the two norms, and
      ! that is what the candidate is held to.  A candidate that does not
      ! beat it is undone: the pair (T, f_sp) is the whole iterate of this
      ! loop, so restoring it puts the iteration back exactly where the
      ! extrapolated pass started and the next pass reproduces the pass the
      ! plain iteration would have taken.  Nothing is lost but that one
      ! pass, and no further limit is taken in this step.
      if (csm_dT .le. csm_x_pred_dT .and. csm_dcomp .le. csm_x_pred_dc)   &
      then
         n_csm_x_kept = n_csm_x_kept + 1
         if (csm_geom_debug)                                            &
            write(*,'(a,i0,a,i0,a,es10.3,a,es10.3)')                      &
               ' CSMXTR step ', count, ' pass ', i_csm, ' KEPT dT= ',     &
               csm_dT, ' dC= ', csm_dcomp
      elseif (.not. csm_extrap_strict .and.                               &
              csm_dT .le. csm_x_at_dT .and.                               &
              csm_dcomp .le. csm_x_at_dc) then
         ! The candidate's own pass did not beat the pass it replaced, but
         ! it did contract: the state it left is closer to the fixed point
         ! than the iterate the candidate was taken from, so undoing it
         ! would throw away a pass of real progress to buy back a pass that
         ! was already spent.  It is kept and no further candidate is taken
         ! in this step, because the model has been shown not to hold here.
         ! MEASURED on mol_metals at 400 steps: undoing these instead
         ! reads 7.29 passes a step against 6.94 for keeping them and 6.99
         ! for not extrapolating at all, so the undone passes outweigh
         ! everything the accepted candidates win on that case.
         n_csm_x_slack = n_csm_x_slack + 1
         csm_x_blocked = .true.
         if (csm_geom_debug)                                            &
            write(*,'(a,i0,a,i0,a,es10.3,a,es10.3)')                      &
               ' CSMXTR step ', count, ' pass ', i_csm, ' SLACK dT= ',    &
               csm_dT, ' dC= ', csm_dcomp
      else
         ! The candidate's own pass moved the pair further than the
         ! candidate itself stood from its predecessor: the extrapolation
         ! did not act on this sequence at all.  The pair (T, f_sp) is the
         ! whole iterate of this loop, so restoring it puts the iteration
         ! back exactly where the extrapolated pass started.
         n_csm_x_refused = n_csm_x_refused + 1
         T(1:N)      = csm_x_T_save(1:N)
         f_sp(1:N,:) = csm_x_f_save(1:N,:)
         csm_x_blocked = .true.
         ! The restored pair is the one the extrapolated pass started
         ! from, so this pass's increment does not continue the sequence
         ! the previous one belongs to and a fresh pair has to be
         ! collected before the next estimate.
         csm_geom_seq_broken = .true.
         if (csm_geom_debug)                                            &
            write(*,'(a,i0,a,i0,a,es10.3,a,es10.3)')                      &
               ' CSMXTR step ', count, ' pass ', i_csm, ' UNDONE dT= ',   &
               csm_dT, ' dC= ', csm_dcomp
      endif
      end subroutine coupled_pair_extrapolation_verdict

      ! ------------------------------------------------------!

      subroutine stationary_state_of_the_loaded_restart
      ! THE STATIONARY RESIDUAL AND CERTIFICATION OF A LOADED STATE, MEASURED
      ! ON THE STATE AS LOADED, AND THE STATIONARY SOLVE ENTERED FROM IT
      ! (docs/restart_contract_design_20260909.md section 2,
      ! "Restart intent: stationary").
      !
      ! WHY THE STATE AS LOADED. A stationary state is a state the equations
      ! hold on; whether they hold is a question about THAT state and about
      ! no other. Every path into the solver until now went through the
      ! marching loop, which takes at least two CFL steps before it hands
      ! off, and a step is a different state: the O(dt) kick of a reload
      ! moves the state before anything measures it, so what was certified
      ! was never the state the file carried. This route takes no step.
      !
      ! WHAT IS REBUILT, AND IN WHICH ORDER. Nothing of the conserved state
      ! is touched; what is rebuilt is what is DERIVED from it, in the one
      ! order in which each consumer reads a quantity that already belongs to
      ! this state: the composition first, because the pressure map and the
      ! boundary read it; the boundary second, so the face states are this
      ! state's; then the primitive variables and the temperature. That is
      ! the order rebuild_state_from_checkpoint uses, for the same reason.
      !
      ! THE LOADED COMPOSITION IS NOT PUT ON ITS OWN FIXED POINT FIRST. The
      ! residual contains one equilibrium sweep from the composition it is
      ! given (eval_residual does exactly this), so the sweep below is part
      ! of the residual's definition and not a change of the state. What
      ! equilibrate_loaded_composition would do instead is iterate that
      ! sweep to a fixed point and rebuild the conserved state from the
      ! result, which moves u: a state whose composition is not the sweep's
      ! own fixed point is a FINDING about the file, reported here as the
      ! particle-count change of the single sweep, and hiding it would make
      ! the measurement answer for a state the file does not carry. It is
      ! available, and reported when taken, as the second word "equilibrate".
      real*8, dimension(1-Ng:N+Ng) :: npart_entry
      real*8  :: dnp
      integer :: jj, info_jfnk

      write(*,'(A)') ' (EXHALE_main) Restart intent: stationary -- the'//   &
           ' loaded state is measured as it stands; no CFL step is taken.'

      ! ONE DISCRETIZATION for the measurement and the solve, and it is the
      ! one a stationary state is a state of: a PLM+WENO3 run starts its
      ! marching in PLM and switches later, so without this the residual of
      ! a WENO3 solution would be assembled by the other operator.
      rec_method   = 'WENO3'
      use_weno3    = .true.
      use_plm      = .false.
      in_plm_stage = .false.

      if (stationary_equilibrate_loaded) then
         write(*,'(A)') '   "equilibrate" was asked for: the loaded'//      &
              ' composition is put on its own fixed point BEFORE the'
         write(*,'(A)') '   evaluation, so what is measured below is that'//&
              ' state and not the one the file carries.'
         call equilibrate_loaded_composition
      endif

      rho = u(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,             &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call Apply_BC(u)
      call U_to_W(u,W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,             &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call comp_T_from_p(p,n_tot,ne,T)
      npart_entry = n_tot + ne
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
      call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,             &
                                 nheiii,nheiTR,nm,ne,n_tot)
      ! T is the temperature of the composition beside it (the same line the
      ! steady outer pass carries, and for the same reason).
      call comp_T_from_p(p,n_tot,ne,T)
      call molecular_carrier_densities_from_state(rho,f_sp)

      ! HOW FAR THE FILE'S COMPOSITION IS FROM THE SWEEP'S OWN ROOT AT THIS
      ! (rho, T). At a state the sweep reproduces, this is the round-off of
      ! the file's own digits; anything larger is a property of the state
      ! that was written and is reported as such, never removed here.
      dnp = 0.0d0
      do jj = 1, N
         dnp = max(dnp, abs((n_tot(jj) + ne(jj) - npart_entry(jj))          &
                            /max(npart_entry(jj), 1.0d-99)))
      enddo
      write(*,'(A,ES10.3)') '   the loaded composition against the'//       &
           ' sweep''s own root: max |d(n_tot+n_e)|/(n_tot+n_e) =', dnp

      resid_max = resid_th
      if (resid_max .le. 0.0d0) resid_max = 1.0d-5

      call assemble_residual(u, n_tot + ne, heat, cool, Rres)
      call residual_norms(Rres, u, resid_c)
      write(*,'(A)') '   stationary residual of the state AS LOADED'//      &
           ' (max over cells of |R|/scale):'
      do jj = 1,3
         write(*,'(A,I2,4X,ES16.6)') '     k=', jj, resid_c(jj)
      enddo
      write(*,'(A,ES12.4)') '     ||R|| = max_k : ', maxval(resid_c)
      call write_residual_breakdown(Rres, u, heat, cool,                    &
                                    'state as loaded (stationary intent)')

      ! THE CLAIM THIS EVALUATION ANSWERS IS THE FILE'S OWN. A state written
      ! as certified is a stationary claim about itself, and re-measuring it
      ! either confirms the claim or refuses it, which is what the exit
      ! status says (contract section 8). A relaxation snapshot claimed
      ! nothing, and measuring it refuses nothing.
      call certification_note_stationarity_claim(ic_certified)
      call certification_evaluate(cert_context_stationary, u, Rres, f_sp,   &
               resid_th, n_cells_without_chemical_root(last_sweep%acc_n),   &
               .true., cert_now)
      call certification_report_write(cert_now,                             &
           'state as loaded (Restart intent: stationary)')

      if (stationary_evaluate_only) then
         ! The state is returned unchanged, with the certification just made
         ! on it: the files this run writes carry the state the files it read
         ! carried, and a 'certified=' field that is a statement about THAT
         ! state (set_state_certified, called by the evaluation above).
         call write_run_counter_report
         call element_census_reservoir('output (stationary evaluation)',    &
                                       rho, f_sp)
         call write_output(rho,v,p,T,heat,cool,eta,                         &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
         call assert_written_state_is_the_accepted_one
         write(*,'(A)') ' (EXHALE_main) Restart intent: stationary'//       &
              ' evaluate -- the state was measured and written back'//      &
              ' unchanged; no step and no solve were taken.'
         call certification_stop_uncertified
         stop
      endif

      ! THE SOLVE, ENTERED AT ONCE. dt here is not a step: eval_dt gives the
      ! CFL interval, which is the pseudo-time the continuation starts from
      ! (the direct steady route uses the same number for the same purpose).
      ! No hydro stage is taken with it.
      call eval_dt(W, dt, dt_loc)
      ! THE CONTINUATION OF A RELOADED STATE STARTS WHERE THE MARCHING
      ! HAND-OFF STARTS IT, at a pseudo-time of 1.0, not at the CFL interval:
      ! a state written by a run is already relaxed, and a continuation
      ! started at the CFL dt from it does not reach the root.  MEASURED
      ! (2026-09-11, docs/Update_EXHALE.md section 8, 8 threads): the
      ! hydrodynamic solve of the hot-Uranus carrier reload stands at
      ! ||R|| 1.13 after 40 iterations from the CFL start and at 1.5e-9 in
      ! 12 iterations from 1.0; the HD 209458 b element reload stagnates at
      ! 1.49 from the CFL start and reaches 1.2e-8 from 1.0.  The direct
      ! steady route keeps the CFL start, which is the pseudo-transient
      ! continuation it was designed as.  EXHALE_PTC_DTAU0=<val> replaces
      ! the start for measurement.
      dt = 1.0d0
      call get_environment_variable('EXHALE_PTC_DTAU0', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) dt
      write(*,'(A,ES10.2)') ' (EXHALE_main) stationary restart: dtau0 = ', dt
      call steady_wind_with_element_diffusion(                            &
              jfnk_outer_iterations_default, dt, .true., info_jfnk)
      call write_run_counter_report
      call element_census_reservoir('output (stationary restart)',          &
                                    rho, f_sp)
      call molecular_carrier_densities_from_state(rho,f_sp)
      call certification_note_stationarity_claim(info_jfnk .eq. 0)
      call assemble_residual(u, n_tot + ne, heat, cool, Rres)
      call certification_evaluate(cert_context_stationary, u, Rres, f_sp,   &
               resid_th, n_cells_without_chemical_root(last_sweep%acc_n),   &
               .true., cert_now)
      call certification_report_write(cert_now, 'final state, as written')
      call write_output(rho,v,p,T,heat,cool,eta,                            &
                        nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
      call assert_written_state_is_the_accepted_one
      write(*,'(A,I0,A)') ' (EXHALE_main) Restart intent: stationary --'//  &
           ' the stationary solve returned info = ', info_jfnk,             &
           '; output written.'
      call certification_stop_uncertified
      stop
      end subroutine stationary_state_of_the_loaded_restart

      ! ------------------------------------------------------!

      subroutine equilibrate_loaded_composition
      ! PUT A STATE READ FROM A FILE ON ITS OWN EQUILIBRIUM COMPOSITION,
      ! ONCE, BEFORE THE FIRST HYDRO STEP.
      !
      ! WHY. Hydro_ioniz_IC.txt and Ion_species_IC.txt carry the state to the
      ! precision the writer emits, and the species fractions they carry are
      ! therefore NOT the equilibrium fractions of the (rho, T) they sit
      ! beside. The marching step is
      !
      !    hydro RK3  ->  T = p/(n_tot+n_e)  ->  ioniz_eq  ->  p = (n_tot+n_e) T
      !
      ! and its last stage rewrites the pressure at the particle count the
      ! sweep just produced. On the first step of a restart that stage is not
      ! a rate but a FINITE JUMP: the particle count moves by the whole
      ! truncation error of the file, and the energy moves with it. Measured
      ! on the two states of docs/p55_base_mode.md, the jump is the same size
      ! at CFL 0.6, 0.06, 0.006 and 0.0006 -- it does not scale with dt, which
      ! is what makes it a projection and not a term of the equation -- and it
      ! reaches 6.2e-3 of the cell energy in one step at r = 1.10 on
      ! WASP-121b, twelve times that step's entire hydro contribution.
      !
      ! WHAT IS HELD FIXED. rho, v and p -- the primitive state the file
      ! carries -- and the composition is iterated to the fixed point of
      !
      !    f_sp  <-  ioniz_eq( rho, T )   with   T = p/(n_tot(f_sp) + n_e(f_sp))
      !
      ! so that the temperature is a derived, self-consistent quantity rather
      ! than the file's second, independently truncated number. At that fixed
      ! point the marching step's own sweep returns the composition it was
      ! given, comp_p_from_T returns the pressure it started from, and the
      ! projection is the identity to round-off. It is the SAME sweep the
      ! marching loop and the steady residual call, with the same secondary-
      ! ionization arming and the same excited-hydrogen refresh, called from
      ! the one place where the two disagree.
      !
      ! It is a Picard iteration and it is not guaranteed to contract, so it
      ! is capped and the achieved change in the particle count is reported;
      ! a single sweep is its first iterate, and on the states measured the
      ! second iterate already moves n_tot + n_e by less than 1e-12.
      !
      ! COST: a few ionization sweeps once per run, against 1e4 or more in the
      ! loop. The alternative -- writing the restart files with more digits --
      ! is cheaper still but only pushes the jump down by the digits added and
      ! leaves T and p in the file independently rounded; this makes the state
      ! self-consistent whatever the file's precision.
      integer, parameter :: nsweep_max = 64
      ! Floor of the iteration, not a target accuracy: the cell solve stops at
      ! its own xtol = sqrt(machine epsilon), so the particle count it returns
      ! cannot be reproduced from one sweep to the next below about 1e-11.
      ! Measured on both states of docs/p55_base_mode.md; the sweep history is
      ! printed when the cap is reached so a state that is NOT contracting is
      ! visible rather than silently accepted.
      real*8,  parameter :: dnpart_tol = 1.0d-10
      real*8, dimension(1-Ng:N+Ng) :: npart_prev
      real*8, dimension(64) :: dnp_hist
      integer :: isw
      real*8  :: dnp

      call U_to_W(u,W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,             &
                                 nheiii,nheiTR,nm,ne,n_tot)
      ! A1: the restart equilibration re-solves the composition of a loaded
      ! state, so the element ratios it carries must survive it.
      call element_census_take('restart equilibration', rho, f_sp, cen_main)
      reload_eq_sweeps = 0
      reload_eq_dnpart = huge(1.0d0)
      do isw = 1, nsweep_max
         npart_prev = n_tot + ne
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
         reload_eq_sweeps = isw
         dnp = 0.0d0
         do j = 1, N
            dnp = max(dnp, abs((n_tot(j) + ne(j) - npart_prev(j))          &
                               /max(npart_prev(j), 1.0d-99)))
         enddo
         reload_eq_dnpart = dnp
         dnp_hist(isw)    = dnp
         if (dnp .lt. dnpart_tol) exit
      enddo
      ! rho is an output of the sweep (calc_rho on the equilibrium species),
      ! and p is the file's pressure, held: the conserved state is rebuilt
      ! from the two so that u and the composition describe one state.
      call comp_T_from_p(p,n_tot,ne,T)
      W(1,:) = rho
      W(2,:) = v
      W(3,:) = p
      call W_to_U(W,u)
      call Apply_BC(u)
      write(*,'(A,I0,A,ES10.3,A,ES10.3)') ' (EXHALE_main) restart'//        &
           ' composition equilibrated in ', reload_eq_sweeps,               &
           ' sweeps; max |d(n_tot+n_e)|/(n_tot+n_e) removed =',             &
           dnp_hist(1), ', left =', reload_eq_dnpart
      if (reload_eq_dnpart .ge. dnpart_tol) then
         write(*,'(A)') '     WARNING: the composition did not reach the'// &
              ' fixed point within the sweep cap; the first hydro step'//   &
              ' still carries part of the restart projection.'
         write(*,'(A,64(1X,ES9.2))') '     sweep history:',                 &
              (dnp_hist(isw), isw = 1, reload_eq_sweeps)
      endif
      call element_census_verify(cen_main, rho, f_sp)
      end subroutine equilibrate_loaded_composition

      subroutine physical_handoff_check
      ! THE STATE A PHYSICAL INTEGRATION STARTS FROM, TESTED BEFORE IT STARTS
      ! (docs/a0_run_mode_contract_20260906.md section 3, "Initialization to
      ! physical integration").
      !
      ! A trajectory is a sequence of states the equations describe. If the
      ! first one is not such a state, every later one inherits that: the
      ! elapsed time then labels a sequence of numerical iterates, which is
      ! the confusion the run modes exist to remove. So four things are
      ! required of the handoff state, and a failure stops the run rather
      ! than repairing anything.
      !
      !  1. ADMISSIBLE. Positive density, pressure and temperature in every
      !     physical cell, and a positive internal energy in the conserved
      !     state: outside that set the Euler equations have no sound speed
      !     and the equation of state no temperature.
      !  2. SOURCE-CONSISTENT CHEMISTRY. Every cell's composition is a ROOT
      !     of the chemical network it is solved from. The sweep's acceptance
      !     classes carry this: classes 1, 2, 3 and 5 each passed the same
      !     two tests, a state inside the element simplex and a normalized
      !     reaction residual at or below ieq_res_tol, and classes 4 and 6
      !     are the two that are not roots -- the state accepted under the
      !     relaxation amnesty, and the one whose candidate was above the
      !     amnesty's residual cap so the cell kept the composition it
      !     entered the sweep with. The amnesty is a device of
      !     initialization and has no place in a trajectory, and a kept
      !     entry composition solved an earlier cell state and not this
      !     one. The element and charge constraints come with it:
      !     the simplex test IS the element budget, and the electron density
      !     is formed from the stage charges of that same composition.
      !  3. CARRIER BACKGROUNDS ESTABLISHED. Where the molecular carriers are
      !     transported, their frozen background must have been filled by a
      !     sweep; a transport step before that runs on zeros.
      !  4. CACHES REBUILT FROM THIS STATE. The sweep below is taken on the
      !     handoff state itself, so the rates, the photon columns, the
      !     excited populations and the molecular columns that the first step
      !     reads belong to it and not to the last initialization iterate.
      !
      ! Steps 1 to 4 of the routine are that list in order.
      integer :: jj, n_bad_chem, j_bad, n_nonpos
      logical :: ok

      write(*,'(A)') ' (EXHALE_main) physical integration: testing the'//  &
           ' handoff state (A0 contract section 3)'

      ! 4 first: the caches are rebuilt, and the tests below then measure the
      ! state as the first step will see it. f_sp is this state's here for
      ! the same reason as in the residual diagnostic above: init and
      ! equilibrate_loaded_composition, the only two paths that reach this
      ! call, each end with get_species_densities, so the caloric arrays
      ! U_to_W reads already hold this composition.
      call Apply_BC(u)
      call U_to_W(u,W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,            &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call comp_T_from_p(p,n_tot,ne,T)
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
      call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,            &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call comp_T_from_p(p,n_tot,ne,T)
      call molecular_carrier_densities_from_state(rho,f_sp)

      ok = .true.

      ! 1. Admissibility.
      n_nonpos = 0
      j_bad    = 0
      do jj = 1, N
         if (.not. (rho(jj) .gt. 0.0d0 .and. p(jj) .gt. 0.0d0             &
                    .and. T(jj) .gt. 0.0d0)) then
            n_nonpos = n_nonpos + 1
            if (j_bad .eq. 0) j_bad = jj
         endif
      enddo
      if (n_nonpos .gt. 0) then
         ok = .false.
         write(*,'(A,I0,A,I0,A,F9.5)') '   REFUSED: ', n_nonpos,          &
              ' cell(s) are not admissible; first is cell ', j_bad,       &
              ' at r =', r(j_bad)
         write(*,'(A,3(1X,ES12.5))') '     rho, p, T (code units):',      &
              rho(j_bad), p(j_bad), T(j_bad)
      endif
      if (.not. positive_density_and_internal_energy(u)) then
         ok = .false.
         write(*,'(A)') '   REFUSED: the conserved state does not have a'//&
              ' positive density and internal energy in every cell.'
      endif

      ! 2. Source-consistent chemistry.
      n_bad_chem = n_cells_without_chemical_root(last_sweep%acc_n)
      if (n_bad_chem .gt. 0) then
         ok = .false.
         j_bad = 0
         if (allocated(ieq_nonroot_streak)) then
            do jj = 1, N
               if (ieq_nonroot_streak(jj) .gt. 0) then
                  j_bad = jj;  exit
               endif
            enddo
         endif
         write(*,'(A,I0,A)') '   REFUSED: ', n_bad_chem, ' cell(s) of the'//&
              ' handoff sweep were accepted as class 4 or class 6, the two'
         write(*,'(A)') '     acceptance classes that are NOT roots of'//  &
              ' the chemical network (the relaxation amnesty, and the'
         write(*,'(A)') '     non-root above its residual cap whose cell'//&
              ' kept the composition it entered the sweep with).'
         if (j_bad .gt. 0) then
            write(*,'(A,I0,A,F9.5,A,ES11.4)') '     first such cell: ',   &
                 j_bad, ' at r =', r(j_bad), ' R_p, T =', T(j_bad)*T0
         else
            write(*,'(A)') '     the sweep did not leave a cell index'//  &
                 ' with a standing non-root streak.'
         endif
      endif

      ! 3. Carrier backgrounds.
      if (thereis_mol .and. carrier_transport .and. .not. bg_ready) then
         ok = .false.
         write(*,'(A)') '   REFUSED: the molecular carriers are'//         &
              ' transported but their frozen background is not'
         write(*,'(A)') '     established, so the first transport step'//  &
              ' would run on zeros.'
      endif

      if (.not. ok) then
         write(*,'(A)') '   The state is not one the equations describe,'//&
              ' so no trajectory starts from it. Relax it first with'
         write(*,'(A)') '   "Run mode: init" and restart the physical'//   &
              ' integration from the state that run writes. Aborting.'
         error stop 1
      endif

      ! The handoff is made. From here the ledgers of accepted steps are the
      ! ones being written (contract section 5): what was counted while the
      ! state was being reached is a diagnostic of a relaxation and is kept
      ! under the initialization family.
      ledger_family = ledger_family_phys
      write(*,'(A,ES14.7,A)') '   handoff accepted; physical time'//       &
           ' origin t_phys =', t_phys, ' s'
      end subroutine physical_handoff_check

      ! ------------------------------------------------------!

      subroutine parse_update_map_request(spec)
      ! EXHALE_UPDATE_MAP=<f1,f2,...>, the time step sizes as FRACTIONS of
      ! the CFL step, at most 8 of them; "1" alone means the default triple
      ! 1, 0.1, 0.01 the review asks for ("at least three decreasing time
      ! steps").
      character(len=*), intent(in) :: spec
      integer :: i, i0, ln, ios
      character(len=64) :: tok
      ln = len_trim(spec)
      if (trim(spec) .eq. '1') then
         upmap_n = 3
         upmap_f(1) = 1.0d0;  upmap_f(2) = 0.1d0;  upmap_f(3) = 0.01d0
         return
      endif
      upmap_n = 0;  i0 = 1
      do i = 1, ln+1
         if (i .gt. ln .or. spec(i:i) .eq. ',') then
            if (i .gt. i0) then
               tok = spec(i0:i-1)
               if (upmap_n .lt. 8) then
                  read(tok,*,iostat=ios) upmap_f(upmap_n+1)
                  if (ios .eq. 0) upmap_n = upmap_n + 1
               endif
            endif
            i0 = i + 1
         endif
      enddo
      end subroutine parse_update_map_request

      ! ------------------------------------------------!

      subroutine update_map_begin_step
      ! Restore the state, and (once) measure the steady residual of it.
      !
      ! The comparison the review's B2 asks for is between the STEADY
      ! RESIDUAL of a state and the PRODUCTION UPDATE MAP applied to the same
      ! state. The production map here is the marching loop's own body -- the
      ! step that follows this call is the step the run takes, with stop
      ! logic and file output suppressed by the caps set at startup -- and NOT
      ! a copy of it, which is the only way the comparison can be trusted.
      if (upmap_done .eq. 0) then
         u_umSave    = u
         f_sp_umSave = f_sp
         if (thereis_mol .and. carrier_transport) then
            upmap_carrier = .true.
            allocate(nH2_um0(1:N), R_carr_um(1:N), T_carr_um(1:N))
            allocate(res_um_all(1:N,4), terms_um_all(1:N,4))
            allocate(y_um0(1:N), nrho_um0(1:N))
         endif
         ! The residual of the state, by the same route eval_residual takes.
         ! The state here is the handoff state, or a state a refused
         ! attempt restored, and both leave f_sp and the caloric arrays
         ! agreeing (rebuild_state_from_checkpoint refreshes the
         ! composition first), so the pressure map below reads this
         ! composition and not the previous one.
         !
         ! THE INVARIANT Apply_BC CARRIES is measured on the way past, not
         ! assumed: it writes the ghost cells and hands the interior back
         ! unchanged, bit for bit, so the number below is required to be
         ! EXACTLY ZERO and anything else is that invariant broken. It is
         ! not a caloric round trip: no interior cell is converted to
         ! primitive variables and back.
         u_umA = u
         call Apply_BC(u)
         upmap_bc_interior = 0.0d0
         do j = 1, N
            do k = 1, 3
               if (u_umA(k,j) .ne. 0.0d0)                                  &
                  upmap_bc_interior = max(upmap_bc_interior,               &
                       abs(u(k,j)-u_umA(k,j))/abs(u_umA(k,j)))
            enddo
         enddo
         write(*,'(A,ES10.3,A)') ' [update map] Apply_BC moves the'//     &
              ' INTERIOR by, at most, a relative', upmap_bc_interior,     &
              ' (required exactly zero: it writes only the ghosts)'
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call assemble_residual(u, n_tot + ne, heat, cool, R_um)
         ! The carrier row's steady residual of the SAME state, per cell,
         ! with the row's own terms beside it.
         if (upmap_carrier) then
            call carrier_steady_residual(rho, v, f_sp, upmap_crmax,     &
                 upmap_cj, upmap_cic, res_out=res_um_all,                  &
                 terms_out=terms_um_all)
            R_carr_um = res_um_all(1:N,1)
            T_carr_um = terms_um_all(1:N,1)
         endif
         u_umSave    = u
         f_sp_umSave = f_sp
      else
         ! COMPOSITION BEFORE PRESSURE. This branch INSTALLS a saved state,
         ! so the caloric arrays of caloric_eos still hold the composition
         ! of whatever state the code last refreshed them from, and
         ! U_to_W's energy-to-pressure map reads them. rho is u(1,:)
         ! exactly, which is the assignment U_to_W itself makes, so the
         ! density the composition needs is available before the pressure
         ! map is built.
         !
         ! THE LOWER BOUNDARY WITH IT (the order of
         ! rebuild_state_from_checkpoint: composition, boundary, primitive
         ! state). The installed u carries the ghost cells and the base face
         ! state of whatever state the code last applied the boundary to,
         ! and Rec_BC reads that face state in EVERY Reconstruct, so the
         ! first Reconstruct of the step that follows would build the base
         ! flux of this step from another state's boundary. Apply_BC derives
         ! the ghosts from the interior alone and writes nothing else, so on
         ! the installed interior it restates this state's own boundary.
         u    = u_umSave
         f_sp = f_sp_umSave
         rho  = u(1,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call Apply_BC(u)
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:)
         call comp_T_from_p(p,n_tot,ne,T)
      endif
      u_um0 = u
      if (upmap_carrier) then
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
         do j = 1, N
            nH2_um0(j)  = f_sp(j,isp_H2)*rho(j)*n0
            y_um0(j)    = f_sp(j,isp_H2)
            nrho_um0(j) = rho(j)*n0
         enddo
      endif
      end subroutine update_map_begin_step

      ! ------------------------------------------------!

      subroutine stop_on_exhausted_retry_budget
      ! THE RETRY BUDGET OF THE TRANSACTION IS SPENT. Reached from the
      ! rejection of a pass and from the rejection of the error estimate,
      ! so that one exhaustion is reported one way whichever refused it;
      ! the caller has set as_verd, err_pass_name, n_retry_now and the
      ! retry history before calling.
      !
      ! IN BOTH MODES AND IN EVERY PASS. The step has
      ! been refused at every interval the policy offered, so the
      ! code cannot produce an acceptable state from the one it
      ! holds, and that is the same statement in a physical run
      ! and in a relaxation: neither has anywhere to go. Inside a
      ! sampled macrostep the exhausted budget is the macrostep's
      ! and the state restored is the one the macrostep began at,
      ! which is the last accepted state; outside one it is the
      ! step's own and the checkpoint of the attempt is that same
      ! state. Either way the run says what refused it, in which
      ! pass, at which interval, and how the refusals changed as
      ! the interval was narrowed; then it exits with status 2,
      ! the status certification.f90 uses for a refused
      ! stationary claim, so that one status means one thing to a
      ! caller reading it from the shell. No fabricated state, no
      ! floor clamp, no continuation.
      n_step_exhaustions = n_step_exhaustions + 1
      if (n_err_passes .eq. 3) then
         call attempted_step_checkpoint_restore(as_entry, u,     &
                                                f_sp, heat,      &
                                                cool, eta)
      else
         call attempted_step_checkpoint_restore(as_chk, u, f_sp, &
                                                heat, cool, eta)
      endif
      call rebuild_state_from_checkpoint
      write(*,*)
      write(*,'(A,I0,A,I0,A)') ' ATTEMPTED STEP EXHAUSTED at'//  &
           ' step ', count, ' after ', n_retry_now + 1,          &
           ' attempts.'
      if (n_err_passes .eq. 3)                                   &
         write(*,'(A,A,A)') '   the pass that refused it: ',     &
              trim(err_pass_name), ' of the error-control'//     &
              ' macrostep'
      write(*,'(A,A)') '   the operation that refused it: ',     &
           trim(attempted_step_operation_text(as_verd%operation))
      write(*,'(A,A)') '   the reason: ',                        &
           trim(attempted_step_reason_text(as_verd%reason))
      write(*,'(A,A)') '   ', trim(as_verd%what)
      if (as_verd%jworst .gt. 0)                                 &
         write(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)')                &
              '   cell ', as_verd%jworst, ', row ',              &
              as_verd%kworst, ', measure ', as_verd%measure,     &
              ', tolerance ', as_verd%tol
      if (n_err_passes .eq. 3) then
         write(*,'(A,ES12.5,A,I0,A)') '   the macro-interval'//  &
              ' reached ', dt_macro*R0/v0, ' s after ',          &
              n_step_retry_max, ' reductions.'
      else
         write(*,'(A,ES12.5,A,I0,A)') '   the step reached ',    &
              dt, ' (code units), the floor of dt/2**',          &
              n_step_retry_max, '.'
      endif
      write(*,'(A)') '   The state restored is the LAST ACCEPTED'&
           //' one and carries its own elapsed time.'
      write(*,'(A,ES23.16)') '   restored-state: state=',       &
           attempted_step_checkpoint_checksum(u, f_sp, heat,    &
                                              cool, eta)
      write(*,'(A)') '   the retry history of this step:'
      do k_as = 0, n_retry_now
         write(*,'(A,I2,A,ES12.5,A,A,A,A)') '     attempt ',     &
              k_as + 1, '  dt = ', as_hist_dt(k_as),             &
              '  refused at ',                                   &
              trim(attempted_step_operation_text(               &
                   as_hist_op(k_as))), ': ',                     &
              trim(attempted_step_reason_text(                   &
                   as_hist_reason(k_as)))
      enddo
      write(*,'(A)') '   The code could not produce an'//        &
           ' acceptable state, so the run exits with status 2.'
      flush(6)
      stop 2
      end subroutine stop_on_exhausted_retry_budget

      subroutine rebuild_state_from_checkpoint
      ! After a restore, the quantities operation 3 derives from u and f_sp
      ! are put back on the restored state. They are NOT in the checkpoint,
      ! deliberately: the membership rule is that an item is saved only if
      ! something reads it before the step rewrites it, and these are
      ! rewritten by operation 3 of the next attempt. Rebuilding them here
      ! rather than saving them is what makes that rule a checked property:
      ! the rollback test hashes them, so a read of one of them before its
      ! rewrite would show up as a hash that does not come back.
      !
      ! THE ORDER BELOW IS THE PHYSICS, NOT A CONVENIENCE. Three maps are
      ! rebuilt and each needs the previous one to be the restored state's:
      !
      !  (1) The composition first. The pressure of a molecular cell is the
      !      inverse of the caloric equation of state, so
      !      pressure_from_energy_density -- and through it U_to_W -- is a
      !      function of the mixture (caloric_eos: nk_per_mass, x_h2,
      !      molecular_cell, caloric_mixture_active). Those arrays are
      !      refreshed only inside get_species_densities, so a U_to_W that
      !      preceded it would map the restored energy density through the
      !      composition of some other state, and the p and T of the
      !      restored state would not be the restored state's. rho is
      !      u(1,:) exactly, which is what U_to_W assigns, so the density
      !      the composition needs is available before the map is built.
      !
      !  (2) The lower boundary next. base_boundary derives the state AT the
      !      base face, the ghost cell averages and the state at the next
      !      face down from the first interior cell average, and Rec_BC
      !      reads that face state in EVERY Reconstruct: it IS the left
      !      state of the base face and therefore one input of the base
      !      HLLC flux. The first Reconstruct after a restore precedes the
      !      pass's own Apply_BC, so the boundary is restated here, on the
      !      restored interior, and the face the next attempt sees belongs
      !      to the state it starts from. Apply_BC writes the ghosts and
      !      nothing else and derives them from the interior alone, so on a
      !      state whose ghosts already state this boundary it returns the
      !      same doubles; what it always does is put base_face_W and
      !      base_face_lower_W back on the restored interior.
      !
      !  (3) The primitive state and the temperature last, from the state
      !      the two steps above have made consistent.
      rho = u(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,           &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call Apply_BC(u)
      call U_to_W(u,W)
      rho = W(1,:)
      v   = W(2,:)
      p   = W(3,:)
      call comp_T_from_p(p,n_tot,ne,T)
      end subroutine rebuild_state_from_checkpoint

      ! ------------------------------------------------------------!

      subroutine update_map_end_step
      ! G_dt = [Phi_dt(q) - q]/dt against the steady residual R(q), by row,
      ! by cell, and by operator-split stage.
      !
      !   hydro      the SSP-RK3 of (dF - S) with Apply_BC after each stage
      !   chem       the composition sweep and the pressure re-derived from T
      !              at the post-sweep particle count (this stage also carries
      !              whatever the element-diffusion and carrier-transport
      !              steps did to f_sp, since those change no conserved
      !              variable directly)
      !   energy     the semi-implicit heating-cooling step
      !   transport  viscosity and conduction (no-op unless asked for)
      !   filter     the Shapiro filter (no-op unless asked for)
      !
      ! The steady residual contains the hydro and heating-cooling terms and
      ! nothing else, so a non-zero chem, transport or filter column is a term
      ! of the production map that the residual does not have -- and, for the
      ! filter, one that the review says cannot be expected to share a fixed
      ! point with it.
      integer :: j, k
      real*8  :: dtl, g, gmax(3), rmax(3), dmax(3)
      real*8  :: smax(5,3)
      character(len=9), parameter :: sname(5) = (/'hydro    ', 'chem     ', &
           'energy   ', 'transport', 'filter   '/)
      character(len=8), parameter :: rname(3) = (/'mass    ', 'momentum',   &
           'energy  '/)
      if (.not. upmap_open()) return
      dtl = dt
      gmax = 0.0d0;  rmax = 0.0d0;  dmax = 0.0d0;  smax = 0.0d0
      do j = 1, N
         do k = 1, 3
            g = (u(k,j) - u_um0(k,j))/dt_loc(j)
            gmax(k) = max(gmax(k), abs(g))
            rmax(k) = max(rmax(k), abs(R_um(k,j)))
            dmax(k) = max(dmax(k), abs(g + R_um(k,j)))
            smax(1,k) = max(smax(1,k), abs(u_umA(k,j)-u_um0(k,j))/dt_loc(j))
            smax(2,k) = max(smax(2,k), abs(u_umB(k,j)-u_umA(k,j))/dt_loc(j))
            smax(3,k) = max(smax(3,k), abs(u_umC(k,j)-u_umB(k,j))/dt_loc(j))
            smax(4,k) = max(smax(4,k), abs(u_umD(k,j)-u_umC(k,j))/dt_loc(j))
            smax(5,k) = max(smax(5,k), abs(u(k,j)   -u_umD(k,j))/dt_loc(j))
            write(upmap_unit,'(1X,ES12.5,1X,I5,1X,I2,1X,10(ES16.8,1X))')   &
                 dtl, j, k, r(j), R_um(k,j), g,                            &
                 (u_umA(k,j)-u_um0(k,j))/dt_loc(j),                        &
                 (u_umB(k,j)-u_umA(k,j))/dt_loc(j),                        &
                 (u_umC(k,j)-u_umB(k,j))/dt_loc(j),                        &
                 (u_umD(k,j)-u_umC(k,j))/dt_loc(j),                        &
                 (u(k,j)   -u_umD(k,j))/dt_loc(j)
         enddo
      enddo
      flush(upmap_unit)
      write(*,'(A,ES11.4,A,ES9.2,A,ES11.4)') ' [update map] dt =', dtl,   &
           '  (', upmap_scale, ' x CFL)   dt_loc(1) =', dt_loc(1)
      write(*,'(A)') '   max over cells, code units 1/t_s'
      write(*,'(A)') '   row       |R|        |G_dt|     |G_dt+R|   '//    &
           'hydro      chem       energy     transport  filter     '//     &
           'Apply_BC interior movement (relative, measured once)'
      do k = 1, 3
         write(*,'(A,A,A,9(ES11.3))') '   ', rname(k), ' ', rmax(k),       &
              gmax(k), dmax(k), smax(1,k), smax(2,k), smax(3,k),           &
              smax(4,k), smax(5,k),                                        &
              upmap_bc_interior
      enddo
      write(*,'(A,ES10.3)') '   largest change of any species fraction'//  &
           ' in the step:', upmap_dfsp
      if (upmap_carrier) call update_map_carrier_row
      upmap_done = upmap_done + 1
      if (upmap_done .ge. upmap_n) max_steps = count
      end subroutine update_map_end_step

      ! ------------------------------------------------!

      subroutine update_map_carrier_row
      ! G_c = [n(H2)_after - n(H2)_before]/dt against -R_carrier of the state
      ! the step started from, cell by cell, in cm^-3 s^-1.  n(H2) is read out
      ! of the state the step leaves, so this carries every stage of the
      ! production map -- the hydro's change of rho at fixed f_sp included --
      ! and not only the transport operator.
      integer :: j, jw, jr
      real*8  :: gc, gmx, rmx, dmx, drel
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,             &
                                 nheiii,nheiTR,nm,ne,n_tot)
      gmx = 0.0d0;  rmx = 0.0d0;  dmx = 0.0d0;  drel = 0.0d0
      jw  = 0;  jr = 0
      do j = 1, N
         ! THE MAP IS MEASURED IN THE ROW'S OWN UNKNOWN -- the carrier
         ! fraction per unit mass, which is what f_sp stores and what the
         ! hydrodynamic stage leaves alone.  The rate of change of the
         ! carrier DENSITY is written beside it: it differs by
         ! n_c d(ln rho)/dt, a hydrodynamic rate that no carrier row
         ! contains and that does not fall with the step (sec. 158).
         gc  = nrho_um0(j)*(f_sp(j,isp_H2) - y_um0(j))/(dt_loc(j)*t_s)
         gmx = max(gmx, abs(gc))
         rmx = max(rmx, abs(R_carr_um(j)))
         if (abs(gc + R_carr_um(j)) .gt. dmx) then
            dmx = abs(gc + R_carr_um(j));  jw = j
         endif
         if (abs(gc + R_carr_um(j))/max(T_carr_um(j),1.0d-300)             &
             .gt. drel) then
            drel = abs(gc + R_carr_um(j))/max(T_carr_um(j),1.0d-300)
            jr   = j
         endif
         write(upmap_unit,'(1X,ES12.5,1X,I5,1X,I2,1X,6(ES16.8,1X))')       &
              dt, j, 4, r(j), R_carr_um(j), gc, gc + R_carr_um(j),         &
              T_carr_um(j),                                                &
              (f_sp(j,isp_H2)*rho(j)*n0 - nH2_um0(j))/(dt_loc(j)*t_s)
      enddo
      write(*,'(A,3(ES11.3),A,I4,A,ES11.3,A,I4,A,F8.4)')                   &
           '   carrier   ', rmx, gmx, dmx, '   worst cell j=', jw,         &
           '   |G+R|/row terms, max =', drel, ' at j=', jr,                &
           ' r=', r(max(jr,1))
      end subroutine update_map_carrier_row

      logical function upmap_open()
      ! Open output/update_map.txt on first use.
      if (.not. upmap_file_open) then
         open(newunit=upmap_unit, file='output/update_map.txt',            &
              status='replace', action='write')
         upmap_file_open = .true.
         write(upmap_unit,'(A)') '# dt j k r  R  G_dt  hydro chem'//       &
              ' energy transport filter    (k = 1 mass, 2 momentum,'//     &
              ' 3 energy; all rates in code units 1/t_s)'
      endif
      upmap_open = .true.
      end function upmap_open

      subroutine report_marching_stop
      ! A DU, DTU OR PLATEAU STOP IS A MARCHING STOP, NOT AN ACCEPTED STEADY
      ! STATE, and the run now says so and prints the three gates measured on
      ! the state it is stopping at.
      !
      ! The three are different concepts and the external review of
      ! 2026-09-03 (section 5, item 7) is explicit that they should not be
      ! forced onto one threshold: du terminates the marching, the du
      ! hand-off threshold decides when the JFNK is offered a state, and the
      ! gates of steady_gates_met decide what may be called steady. This
      ! routine only makes the third visible wherever the first fires; the
      ! du window is deliberately NOT unified with the flux gate's window.
      !
      ! The gates are MEASURED here rather than read from the in-loop
      ! monitor, which is not evaluated at all unless "Resid tol:" was given
      ! and is otherwise stale.
      real*8, dimension(3) :: rc_now
      real*8 :: fsp_now, fmean_now
      call assemble_residual(u, n_tot + ne, heat, cool, Rres)
      call residual_norms(Rres, u, rc_now)
      call flux_spread_of_state(u, fsp_now, fmean_now)
      write(*,'(A)') '        THIS IS A MARCHING STOP, NOT AN ACCEPTED'//    &
           ' STEADY STATE. The three gates on this state:'
      if (resid_th .gt. 0.0d0) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                    &
              '          residual (max over cells of |R|/scale) ',           &
              maxval(rc_now), '   target "Resid tol" ', resid_th,            &
              merge('  MET    ', '  NOT MET', maxval(rc_now) .lt. resid_th)
      else
         write(*,'(A,ES10.3,A)')                                             &
              '          residual (max over cells of |R|/scale) ',           &
              maxval(rc_now), '   no "Resid tol" was set: this run never'//  &
              ' tested it'
      endif
      if (flux_spread_th .gt. 0.0d0) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                    &
              '          flux spread over r >= r_flux             ',         &
              fsp_now, '   target "Flux spread tol" ', flux_spread_th,       &
              merge('  MET    ', '  NOT MET', fsp_now .lt. flux_spread_th)
      else
         write(*,'(A,ES10.3,A)')                                             &
              '          flux spread over r >= r_flux             ',         &
              fsp_now, '   the flux gate is disabled in this run'
      endif
      if (thereis_mol .and. carrier_transport) then
         ! Measured here, on this state. The marching's carriers are moved by
         ! the operator-split transport step and are not solved for, so no
         ! carrier residual OF A SOLVE exists for a marched state; this is the
         ! row evaluated on the state, which is what a reader needs.
         call carrier_steady_residual(rho, v, f_sp, crc_max,              &
                                      crc_j, crc_ic, rvol=crc_vol)
         write(*,'(A,ES10.3,A,ES10.3)')                                      &
              '          carrier row (evaluated here)             ',         &
              crc_vol, '   target "Carrier resid tol" ',                     &
              carrier_resid_th
      endif
      write(*,'(A)') '        (a state that has to be quoted as steady'//    &
           ' must pass all three; see Update_EXHALE_stage1.md section 145)'
      call write_residual_breakdown(Rres, u, heat, cool, 'marching stop')
      end subroutine report_marching_stop


      subroutine assert_written_state_is_the_accepted_one
      ! THE STATE THE GATES ACCEPTED IS THE STATE THAT WAS JUST WRITTEN.
      !
      ! Nothing between solve_steady_jfnk's return and write_output touches
      ! u. Inside steady_wind_with_element_diffusion the refresh recomputes
      ! rho, v, p FROM u (U_to_W) and moves only the composition, heat, cool
      ! and T; the rest of the marching-loop body and the block between the
      ! loop and the final write read u without assigning it. So rho*v*r^2
      ! cannot move, and the flux gate's number is a property of the file the
      ! user reads.
      !
      ! That is a property of the call graph, not of any one run, so it is
      ! ASSERTED on every run rather than left as a claim: the flux spread is
      ! re-measured on the state as written and printed beside the value the
      ! solve accepted, and a difference is a loud warning. It costs one pass
      ! over the wind window.
      !
      ! (What DID differ, and is not a state change at all, is reading the
      ! output FILE without dropping its ghost rows -- see the row header
      ! write_output now emits and section 133.6.)
      real*8 :: fspread_now, fmean_now
      if (gate_fspread_accepted .lt. 0.0d0) return   ! no steady solve ran
      call flux_spread_of_state(u, fspread_now, fmean_now)
      write(*,'(A,ES11.4,A,ES11.4)') ' (EXHALE_main) flux gate: '//       &
           'accepted', gate_fspread_accepted, '   as written', fspread_now
      call flux_spread_above_radius(u, 1.03d0, fspread_103, fmean_now)
      call flux_spread_above_radius(u, 1.10d0, fspread_110, fmean_now)
      write(*,'(A,ES11.4,A,ES11.4,A)') ' (EXHALE_main) flux spread by '//  &
           'window: r>=1.03', fspread_103, '   r>=1.10', fspread_110,      &
           '  (reporting only; the gate is the r>=r_flux window above)'
      ! The comparison carries an ABSOLUTE floor beside the relative one:
      ! the spread is a difference of O(1) numbers (rho v r^2 over the
      ! window) divided by their mean, so two evaluations of one state
      ! agree only to the rounding of that difference, of order N epsilon
      ! ~ 1e-13. MEASURED (2026-09-11, the partitioned hot-Uranus carrier
      ! reload): a state flat to fourteen digits reads 3.0982e-14 as written
      ! against 2.7487e-14 accepted, and a relative test alone called that
      ! a change of u. 1e-12 stands a decade above that rounding and four
      ! decades below the smallest spread a gate has ever accepted (1.5e-15
      ! is the well-balanced wasp_full_newton's, and that state is flat to
      ! the last bit); a real change of u moves the spread by far more.
      if (abs(fspread_now - gate_fspread_accepted) .gt.                   &
          max(1.0d-10*abs(gate_fspread_accepted), 1.0d-12))               &
         write(*,'(A)') ' (EXHALE_main) WARNING: the written state is '// &
              'NOT the state the gates accepted -- something between the'//&
              ' steady solve and write_output changed u.'
      end subroutine assert_written_state_is_the_accepted_one

      ! ------------------------------------------------!

      subroutine steady_wind_with_element_diffusion(maxit, dtau0,       &
                                                   use_jfnk, jfnk_info)
      ! THE STATIONARY STATE OF THE WIND AND OF EVERY SPECIES IT CARRIES,
      ! ACCEPTED BY ONE JOINT TEST.
      !
      ! WHICH EQUATIONS ARE ALTERNATED WITH THE WIND. The steady residual
      ! (steady_residual.f90) carries the three hydrodynamic rows and, where
      ! the run registers them, the transported species rows. A species the
      ! run moves by an operator-split transport step instead -- the
      ! diffusing elements always, the molecular carriers where
      ! "Coupled carrier solve" is off -- has a stationary balance that no
      ! hydrodynamic solve touches, so a single steady solve would freeze it
      ! at whatever state it was handed. Those balances are alternated with
      ! the wind: solve the wind at fixed composition, relax the composition
      ! to ITS steady state in that wind, solve again.
      !
      ! WHAT IS ACCEPTED, AND BY WHICH EVALUATION. Only a state on which the
      ! certification of the FULL set of active equations passes: the three
      ! hydrodynamic rows, the transported carrier balances, the elemental
      ! transport balances, the level populations, the eliminated-species
      ! closure and the conservation records, each against its own tolerance
      ! (certification.f90). It is evaluated on the REFRESHED state, the one
      ! this routine would hand back, from a residual assembled here. A
      ! hydrodynamic solve returning info = 0 while an alternated species row
      ! refuses is NOT an accepted state: that flag is a statement about the
      ! rows in the Newton registry, and the alternated rows stand outside it
      ! by construction.
      !
      ! WHAT ENDS THE LOOP, and what jfnk_info then means:
      !   the joint certification passes             -> 0, the state accepted
      !   a row of the state is not finite           -> 1, named
      !   the element composition update is refused  -> 1, named
      !   the joint measure stops falling            -> 1, named
      !   the pass budget ends uncertified           -> 1, named
      ! No caller separates 1 from 2, and every ending but the first is a
      ! refusal, so one nonzero flag carries them and the line printed beside
      ! it says which. The state handed back on a refusal is the refreshed
      ! state of the last hydrodynamic solve: a stationary wind at a
      ! composition whose own balance still refuses.
      !
      ! A NONZERO HYDRODYNAMIC FLAG DOES NOT END THE LOOP BY ITSELF. info = 1
      ! is the iteration budget, info = 2 a returned state a row refuses, and
      ! both hand back a valid state on which the composition update still
      ! makes progress: MEASURED on the atomic reload, the worst elemental
      ! wind row stood at 2.61e-3 after the first budget-ended solve and at
      ! 2.78e-4 after the update that followed it
      ! (docs/solver_partition_experiment_20260911.md sec. 5.1). Ending the
      ! alternation on that flag discards what the pass achieved. The flag is
      ! printed every pass.
      !
      ! WHAT CONTROLS THE PASSES. The quantity the joint acceptance waits on
      ! is the worst gated species row of the certification: the elemental
      ! and carrier balances are judged in the wind, r >= cert_regime_wind_r,
      ! and driving that measure down is the whole purpose of the
      ! alternation. Where it fails to fall over a pass both updates are
      ! shortened -- the element relaxation's under-relaxation omega is
      ! halved (floor 0.125) and the carrier movement bound is halved (floor
      ! carrier_trust_floor) -- and outer_no_fall_max consecutive passes
      ! without a fall end the loop rather than spend the budget on a fixed
      ! point the alternation is not approaching.
      !
      ! THE COMPOSITION DRIFT AND THE VOLUME-WEIGHTED CARRIER RESIDUAL ARE
      ! REPORTED AND ARE NOT TESTS. The transport operator's smallest step,
      ! one cell crossing time, already moves the composition by 1.1e-2 on
      ! the He/H = 0.0793 hot Uranus (READ from
      ! docs/p50_carrier_wind_alternation.md), so a drift gate below that
      ! stands under anything a pass can produce; and an average over a
      ! column says nothing about the cell in which the equation is worst
      ! satisfied, which is what every row of the certification is measured
      ! by. Neither can be an acceptance.
      !
      ! EVERY ENDING OF THIS ROUTINE HANDS BACK A STATE WHOSE COMPOSITION,
      ! PARTICLE DENSITIES, TEMPERATURE AND RESIDUAL BELONG TOGETHER; A PASS
      ! THAT ENDS THE ITERATION TAKES NO COMPOSITION UPDATE. The caller
      ! certifies and writes the arrays this routine leaves behind
      ! (stationary_state_of_the_loaded_restart reads f_sp for the
      ! certification and n_HI ... n_m and T for the output files), so a
      ! composition that had moved past the particle densities, the
      ! temperature and the residual beside it would be certified as one
      ! state and written as another. This is why the progress control
      ! stands AHEAD of the composition update: the pass it ends takes no
      ! update, and what is handed back is the very state the certification
      ! of that pass was evaluated on. The other endings hold the same way:
      ! a certified or non-finite pass reaches no update, the refused
      ! element update is the entry composition restored, and the pass
      ! budget's last pass takes none for the reason below.
      !
      ! THE LAST PASS TAKES NO COMPOSITION UPDATE. An update exists to be
      ! consumed by the next hydrodynamic solve; with no pass left there is
      ! none, and taking it would replace a stationary wind by a state whose
      ! hydrodynamic rows the update itself spoiled (MEASURED on the
      ! molecular reload: mass row 4.99e-12 after a hydrodynamic solve,
      ! 1.54e-1 after the update that followed it, same document sec. 3).
      !
      ! WITH NO ALTERNATED SPECIES the body is one steady solve followed by
      ! the state refresh and the flag is the solve's own; no joint
      ! evaluation is made and none of the lines below is printed. That is
      ! what an atomic run and a coupled carrier solve both do.
      !
      ! maxit / dtau0 / use_jfnk select the steady solver and its budget:
      ! every route passes jfnk_outer_iterations_default, the marching
      ! hand-off starts the continuation from dtau0 = 1 and the two direct
      ! routes from the CFL dt, and either route may use JFNK or
      ! pseudo-transient continuation.
      integer, intent(in)  :: maxit
      real*8,  intent(in)  :: dtau0
      logical, intent(in)  :: use_jfnk
      integer, intent(out) :: jfnk_info

      ! The named endings of the outer iteration.
      integer, parameter :: outer_running                = 0
      integer, parameter :: outer_certified              = 1
      integer, parameter :: outer_state_not_finite       = 2
      integer, parameter :: outer_element_update_refused = 3
      integer, parameter :: outer_no_progress            = 4
      integer, parameter :: outer_pass_budget            = 5
      ! How many consecutive passes may leave the joint measure not falling
      ! before the loop ends. A single rise is the alternation's own
      ! feedback and not stagnation: MEASURED on the atomic reload the worst
      ! elemental wind row went 2.78e-4 -> 3.81e-3 across a hydrodynamic
      ! solve and 3.81e-3 -> 7.57e-5 across the update that followed it.
      integer, parameter :: outer_no_fall_max = 3
      ! The floor of the carrier movement bound. The bound is enforced on
      ! the accepted step (diffusive_photochemistry.f90) and a shorter step
      ! is always admissible, so halving it cannot make a pass refuse; below
      ! this floor a pass would move the composition by less than a tenth of
      ! what the front covers in one cell crossing time.
      real*8, parameter  :: carrier_trust_floor = 1.0d-3

      integer :: hydro_info, elem_status, carrier_outcome, outer_ending
      integer :: n_no_fall, icert, sp_cell, ref_cell, cell_here, pass_cap
      integer :: n_unjudged
      real*8  :: sp_worst, sp_worst_prev, sp_meas, sp_tol, trust_pass
      real*8  :: row_mass, row_mom, row_ene, t_pass0
      real*8  :: ref_worst, ref_meas, ref_tol, row_here
      logical :: species_alternated, pass_certified, rows_finite
      logical :: update_taken
      character(len=52) :: sp_name, ref_name, unjudged_name
      character(len=76) :: unjudged_why
      character(len=64) :: elem_text

      ! THE CAP NAMED FOR A MEASUREMENT IS THE RUN'S, NOT EACH SOLVE'S (see
      ! jfnk_run_cap_named at the declarations). Read once, at the first
      ! entry; after the first solve under a named cap the run does not
      ! enter the solver again, and says so where it would have.
      if (n_stationary_solves_run .eq. 0) then
         call get_environment_variable('EXHALE_JFNK_MAXIT', jfnk_cap_env)
         jfnk_run_cap_named = (len_trim(jfnk_cap_env) .gt. 0)
      endif
      if (jfnk_run_cap_named .and. n_stationary_solves_run .ge. 1) then
         write(*,'(A,A,A)') ' (EXHALE_main) EXHALE_JFNK_MAXIT = ',        &
              trim(jfnk_cap_env), ' caps the stationary solving of this'// &
              ' RUN, and the run has taken its solve: this entry is not'
         write(*,'(A)') '   taken, so the measurement is of one solve'//   &
              ' from one state.'
         jfnk_info = 1
         return
      endif
      n_stationary_solves_run = n_stationary_solves_run + 1

      comp_omega      = 0.5d0
      ! Diagnostic override (EXHALE_DIFF_OMEGA=<val>): pins the starting
      ! under-relaxation factor, so the undamped loop (1.0) can be compared
      ! against the damped one without a rebuild.
      call get_environment_variable('EXHALE_DIFF_OMEGA', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) comp_omega
      comp_drift      = 0.0d0
      elem_drift      = 0.0d0

      call get_environment_variable('EXHALE_CARRIER_TRUST', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) carrier_trust
      call get_environment_variable('EXHALE_OUTER_PASSES', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) outer_pass_cap
      ! The bound one pass may move the carriers by is shortened by the
      ! progress control below; the run's own setting is left where a later
      ! entry can read it.
      trust_pass = carrier_trust

      species_alternated = he_diffusion .or. (thereis_mol .and.           &
                           carrier_transport .and. .not. carrier_in_newton)
      pass_cap      = merge(outer_pass_cap, 1, species_alternated)
      outer_ending  = outer_running
      hydro_info    = 0
      n_no_fall     = 0
      sp_worst      = -1.0d0
      sp_worst_prev = huge(1.0d0)
      sp_meas       = 0.0d0
      sp_tol        = 0.0d0
      sp_cell       = 0
      sp_name       = 'no gated species row in this configuration'
      ref_cell      = 0
      ref_name      = 'none'
      ref_meas      = 0.0d0
      ref_tol       = 0.0d0
      n_unjudged    = 0
      unjudged_name = 'none'
      unjudged_why  = ''
      row_mass      = 0.0d0
      row_mom       = 0.0d0
      row_ene       = 0.0d0

      do it_diff = 1, pass_cap
         t_pass0 = omp_get_wtime()
         if (use_jfnk) then
            ! THE COUPLED ROUTE (section 139). With the carrier row among
            ! the unknowns the outer loop is not an alternation at all: one
            ! solve returns a wind and a carrier partition that are steady
            ! states of each other, and the loop below runs once.
            ! WHICH BALANCES BECOME ROWS is the registry's own question,
            ! answered from the configuration: the carriers where they are
            ! transported, the elements where they diffuse.  What is asked
            ! here is only whether the transported balances are to be solved
            ! WITH the wind instead of alternated with it, which is one
            ! choice for all of them.
            call set_transported_species_rows(carrier_in_newton)
            call solve_steady_jfnk(u, f_sp, resid_max, maxit, dtau0,    &
                                   40, hydro_info)
         else
            call solve_steady_ptc(u, f_sp, resid_max, maxit, dtau0,     &
                                  hydro_info)
         endif
         call set_transported_species_rows(.false.)
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:);  E = u(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,       &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,       &
                                    nheiii,nheiTR,nm,ne,n_tot)
         ! T IS THE TEMPERATURE OF THE COMPOSITION BESIDE IT. Without this
         ! line T is still p/(n_tot + n_e) evaluated at the particle count
         ! BEFORE the sweep above, so the T column of the file written from
         ! here belongs to a composition one sweep behind its own species
         ! columns. Measured on the hot Uranus steady solution: the two
         ! disagree by up to 1.8e-5, by ~1e-7 elsewhere, and nothing but the
         ! T column moves. p is not touched -- it is the conserved state's
         ! own pressure -- so this makes (p, T, f_sp) a consistent triple.
         call comp_T_from_p(p,n_tot,ne,T)
         ! A1: the state this pass would hand on as accepted, against the
         ! reservoirs the run resolved.
         call element_census_reservoir('steady outer pass (accepted '//   &
              'state)', rho, f_sp)
         ! WITH NOTHING ALTERNATED there is one solve and the flag is its
         ! own: the transported balances, where the run has any, were rows
         ! of that solve.
         if (.not. species_alternated) then
            jfnk_info = hydro_info
            return
         endif

         ! ---- THE JOINT TEST, ON THE STATE THIS PASS WOULD HAND BACK ----
         ! The residual is assembled here, so that every row measured below
         ! -- and the face mass flux the elemental and carrier balances ride
         ! on, which store_row_terms leaves behind for them -- belongs to
         ! the refreshed state and not to the last state the solver happened
         ! to evaluate inside its iteration.
         call assemble_residual(u, n_tot + ne, heat, cool, Rres)
         call certification_evaluate(cert_context_stationary, u, Rres,     &
                  f_sp, resid_th,                                         &
                  n_cells_without_chemical_root(last_sweep%acc_n),        &
                  .true., cert_now)
         pass_certified = cert_now%certified
         ! The worst GATED species row of this state is the measure the
         ! joint acceptance waits on, ranked by distance from its own
         ! tolerance because the elemental and the carrier rows carry
         ! tolerances of their own. The worst REFUSING entry of the whole
         ! inventory is read in the same pass over it, so that a refusal
         ! names an equation and a cell and not a count.
         rows_finite = .true.
         sp_worst = -1.0d0;  sp_meas = 0.0d0;  sp_tol = 0.0d0
         sp_cell  = 0
         sp_name  = 'no gated species row in this configuration'
         ref_worst = -1.0d0;  ref_meas = 0.0d0;  ref_tol = 0.0d0
         ref_cell  = 0
         ref_name  = 'none'
         n_unjudged = 0
         unjudged_name = 'none';  unjudged_why = ''
         do icert = 1, cert_now%n
            if (cert_now%e(icert)%status .eq. cert_unavailable) then
               n_unjudged = n_unjudged + 1
               if (n_unjudged .eq. 1) then
                  unjudged_name = cert_now%e(icert)%name
                  unjudged_why  = cert_now%e(icert)%reason
               endif
               cycle
            endif
            if (cert_now%e(icert)%status .ne. cert_evaluated) cycle
            if (.not. cert_now%e(icert)%finite) rows_finite = .false.
            if (cert_now%e(icert)%regime_gated) then
               row_here  = cert_now%e(icert)%row_max_gate
               cell_here = cert_now%e(icert)%jworst_gate
            else
               row_here  = cert_now%e(icert)%row_max
               cell_here = cert_now%e(icert)%jworst
            endif
            if (cert_now%e(icert)%tol .le. 0.0d0) cycle
            if (cert_now%e(icert)%regime_gated .and.                      &
                row_here/cert_now%e(icert)%tol .gt. sp_worst) then
               sp_worst = row_here/cert_now%e(icert)%tol
               sp_meas  = row_here
               sp_tol   = cert_now%e(icert)%tol
               sp_cell  = cell_here
               sp_name  = cert_now%e(icert)%name
            endif
            if (.not. cert_now%e(icert)%within_tol .and.                  &
                row_here/cert_now%e(icert)%tol .gt. ref_worst) then
               ref_worst = row_here/cert_now%e(icert)%tol
               ref_meas  = row_here
               ref_tol   = cert_now%e(icert)%tol
               ref_cell  = cell_here
               ref_name  = cert_now%e(icert)%name
            endif
         enddo
         icert = certification_entry_index(cert_now,                      &
                                           'hydrodynamic mass row')
         if (icert .gt. 0) row_mass = cert_now%e(icert)%row_max
         icert = certification_entry_index(cert_now,                      &
                                           'hydrodynamic momentum row')
         if (icert .gt. 0) row_mom = cert_now%e(icert)%row_max
         icert = certification_entry_index(cert_now,                      &
                                           'hydrodynamic energy row')
         if (icert .gt. 0) row_ene = cert_now%e(icert)%row_max
         if (pass_certified) then
            outer_ending = outer_certified
         else if (.not. rows_finite) then
            outer_ending = outer_state_not_finite
         endif

         ! ---- THE PROGRESS CONTROL, AHEAD OF ANY UPDATE OF THIS PASS ----
         ! The joint measure is what the passes are spent on, so it is what
         ! the step lengths are chosen by. omega damps the element
         ! relaxation and trust_pass bounds the carrier pass; both are
         ! shortened together, because which of the two failed to move the
         ! measure is not a question this loop can answer. It is decided
         ! here, before the update, so that the pass which ends the
         ! iteration takes no update and hands back the state its own
         ! certification above was evaluated on (the invariant in the
         ! header); the ending itself is announced with the others, below
         ! the summary line of the pass.
         if (outer_ending .eq. outer_running) then
            if (sp_cell .gt. 0 .and. sp_worst .ge. sp_worst_prev) then
               n_no_fall = n_no_fall + 1
               if (he_diffusion .and. comp_omega .gt. 0.125d0) then
                  comp_omega = max(0.5d0*comp_omega, 0.125d0)
                  write(*,'(A,F6.3)') '    -> the worst gated species'//  &
                       ' row did not fall; under-relaxation omega =',     &
                       comp_omega
               endif
               if (thereis_mol .and. carrier_transport .and.              &
                   .not. carrier_in_newton .and.                          &
                   trust_pass .gt. carrier_trust_floor) then
                  trust_pass = max(0.5d0*trust_pass, carrier_trust_floor)
                  write(*,'(A,ES9.2)') '    -> the worst gated species'//  &
                       ' row did not fall; carrier movement bound =',     &
                       trust_pass
               endif
               if (n_no_fall .ge. outer_no_fall_max)                      &
                  outer_ending = outer_no_progress
            else if (sp_cell .gt. 0) then
               n_no_fall = 0
            endif
            sp_worst_prev = sp_worst
         endif

         ! ---- THE COMPOSITION UPDATE THE NEXT SOLVE WILL CONSUME ----
         update_taken = .false.
         comp_drift   = 0.0d0
         kd           = 0
         if (outer_ending .eq. outer_running .and. it_diff .lt. pass_cap) &
             then
            update_taken = .true.
            ! Element composition relaxed to its steady state at the fixed
            ! wind. The wind is the face mass flux of THIS state, read from
            ! the mass row this state assembled and handed to the operator:
            ! the element module is given its advecting flux and does not
            ! reach into the steady residual for one. A run with no helium
            ! relaxes nothing, and no flux is fetched for it.
            if (he_diffusion) then
               Frho_elem = 0.0d0
               if (thereis_He) call face_mass_flux_of_state(rho, Frho_elem)
               call relax_element_composition(rho,v,T,f_sp,Frho_elem,     &
                                              comp_omega,elem_drift,kd,   &
                                              status=elem_status)
               select case (elem_status)
               case (element_relaxation_converged)
                  elem_text = 'the fixed point of the element operator'
               case (element_relaxation_step_budget)
                  elem_text = 'the step budget, the fixed point not'//    &
                              ' reached'
               case default
                  elem_text = 'NO ADMISSIBLE ADVANCE: the entry'//        &
                              ' composition was restored'
               end select
               if (elem_status .eq. element_relaxation_failed) then
                  ! The composition this routine holds is the one the
                  ! relaxation was entered with, restored by the operator,
                  ! so the state handed back is the hydrodynamic solve's
                  ! own and no chemical refresh is due on a composition
                  ! that was never adopted.
                  outer_ending = outer_element_update_refused
               else
                  comp_drift = elem_drift
                  ! THE SWEEP IS ENTERED AT THE RELAXED COMPOSITION'S OWN
                  ! PARTICLE COUNT AND TEMPERATURE. The relaxation moved
                  ! f_sp, so n_tot, n_e and T = p/(n_tot + n_e) still
                  ! belong to the composition it was entered with, and a
                  ! sweep taken at that T equilibrates the new composition
                  ! at another state's temperature. p is the conserved
                  ! state's own pressure and is not touched. The carrier
                  ! update below has the same shape for the same reason.
                  call get_species_densities(rho,f_sp,nhi,nhii,nhei,      &
                                       nheii,nheiii,nheiTR,nm,ne,n_tot)
                  call comp_T_from_p(p,n_tot,ne,T)
                  call ioniz_eq(T,rho,f_sp,heat,cool,eta,last_sweep)
               endif
            endif
            ! The molecular carriers are a THIRD participant in the same
            ! Picard iteration, relaxed towards their own steady state at
            ! the same fixed wind, bounded by trust_pass rather than
            ! damped. The header records why the damping exists at all --
            ! nothing in a Picard iteration of two solves keeps them from
            ! chasing each other -- and a third participant makes that risk
            ! larger.
            !
            ! THE PASS CARRIES THE CHEMISTRY WITH IT. It equilibrates the
            ! eliminated species, the electron density, the temperature and
            ! the rate coefficients on every step it keeps, at this pass's
            ! fixed (rho, v, p), so the state it hands back is already a
            ! state the chemistry has closed on and the sweep that used to
            ! stand here would re-solve its own answer. What is left is to
            ! fill the density columns this routine keeps beside the
            ! composition.
            if (outer_ending .eq. outer_running .and. thereis_mol .and.   &
                carrier_transport .and. .not. carrier_in_newton) then
               call relax_photochemical_composition(rho,v,p,T,f_sp,       &
                                              heat,cool,eta,              &
                                              trust_pass, carrier_drift,  &
                                              kc, outcome=carrier_outcome)
               call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,   &
                                          nheiii,nheiTR,nm,ne,n_tot)
               call comp_T_from_p(p,n_tot,ne,T)
               comp_drift = max(comp_drift, carrier_drift)
               kd = max(kd, kc)
            endif
         endif

         ! ---- ONE LINE PER PASS ----
         write(*,'(A,I0,A,I0,A,ES9.2,A,ES8.1,A,I0,A,A,A,ES9.2,A,ES9.2,A,ES9.2,A,F6.3,A,ES8.1,A,F8.2,A)') &
              ' (EXHALE_main) outer pass ', it_diff, ': hydro info=',     &
              hydro_info, ', worst gated species row ', sp_meas, ' of ',  &
              sp_tol, ' at cell ', sp_cell, ' (', trim(sp_name),          &
              '), mass ', row_mass, ', momentum ', row_mom, ', energy ',  &
              row_ene, ', omega', comp_omega, ', trust', trust_pass,      &
              ',', omp_get_wtime() - t_pass0, ' s'
         if (update_taken .and. he_diffusion)                             &
            write(*,'(A,A,A,ES10.2,A,I0,A)') '    element relaxation'//   &
                 ' ended on ', trim(elem_text), '; drift', elem_drift,    &
                 ' in ', kd, ' steps'
         if (update_taken .and. outer_ending .ne.                         &
             outer_element_update_refused .and. thereis_mol .and.         &
             carrier_transport .and. .not. carrier_in_newton)             &
            write(*,'(A,A,A,ES10.2,A,I0,A)') '    carrier relaxation'//   &
                 ' ended on ', trim(carrier_relax_outcome_text(           &
                 carrier_outcome)), '; drift', carrier_drift, ' in ', kc, &
                 ' transport steps'
         if (.not. update_taken) then
            if (outer_ending .eq. outer_running) then
               write(*,'(A)') '    no composition update was taken:'//    &
                    ' this is the last pass and no solve is left to'//    &
                    ' consume one'
            else
               write(*,'(A)') '    no composition update was taken:'//    &
                    ' this pass ends the iteration, so the state handed'//&
                    ' back is the state certified above'
            endif
         endif

         ! ---- THE DIAGNOSTICS OF THE PASS ----
         ! Is the state a steady state of the CARRIER equation, and where is
         ! it least so? The gated carrier row of the certification is the
         ! acceptance; these numbers say how the refusal is distributed over
         ! the column and where the front now stands.
         if (thereis_mol .and. carrier_transport .and. rows_finite) then
            call carrier_steady_residual(rho, v, f_sp, crc_max,       &
                                         crc_j, crc_ic, rvol=crc_vol,    &
                                         rlegacy=crc_leg)
            if (crc_j .gt. 0)                                            &
               write(*,'(A,ES10.2,A,ES10.2,A,I4,A,F7.3,A,A)')            &
                    '    carrier steady residual: volume-weighted',      &
                    crc_vol, ', worst cell', crc_max, ' at ', crc_j,     &
                    ' r=', r(crc_j), ' carrier ',                        &
                    trim(carrier_name(crc_ic))
            if (crc_j .gt. 0)                                            &
               write(*,'(A,ES10.2)') '    (on the retired '//            &
                    'n_H(|v|+c_s)/dr scale the same residual reads', crc_leg
            ! Where the H2 front stands after this pass: the innermost cell
            ! with 2 n(H2)/n_H below 0.5 and below 1e-2.  The alternation is
            ! a continuation in the front's position, so the position is the
            ! quantity to read per pass.
            block
              integer :: jf, j50, j01
              real*8  :: x2
              j50 = 0;  j01 = 0
              do jf = 1, N
                 x2 = 2.0d0*f_sp(jf,isp_H2)/max(f_sp(jf,1) + f_sp(jf,2)  &
                      + 2.0d0*f_sp(jf,isp_H2), 1.0d-300)
                 if (j50 .eq. 0 .and. x2 .lt. 0.5d0)  j50 = jf
                 if (j01 .eq. 0 .and. x2 .lt. 1.0d-2) j01 = jf
              enddo
              write(*,'(A,F8.4,A,F8.4,A,I0,A,ES9.2,A)')                  &
                   '    H2 front: x2=0.5 at r=', r(max(j50,1)),          &
                   '  x2=1e-2 at r=', r(max(j01,1)), '  (pass ', it_diff, &
                   ', trust', trust_pass, ')'
            end block
            call carrier_drift_location(cdl_j, cdl_ic)
            if (cdl_j .gt. 0)                                            &
               write(*,'(A,I4,A,F7.3,A,A)') '    carrier drift worst '// &
                    'at cell ', cdl_j, ' r=', r(cdl_j), ' carrier ',     &
                    trim(carrier_name(cdl_ic))
         endif
         call get_environment_variable('EXHALE_DIFFUSION_CHECK', diag_env)
         if (trim(diag_env) .eq. '1')                                    &
            call write_diffusion_pass_profile(it_diff)

         ! ---- THE ENDING OF THIS PASS, NAMED ----
         if (outer_ending .eq. outer_certified) then
            jfnk_info = 0
            write(*,'(A,I0,A)') ' (EXHALE_main) outer pass ', it_diff,   &
                 ': ACCEPTED -- every active equation of this state is'// &
                 ' within its own tolerance.'
            call certification_report_write(cert_now,                     &
                 'accepted state of the stationary outer iteration')
            return
         else if (outer_ending .eq. outer_state_not_finite) then
            jfnk_info = 1
            write(*,'(A,I0,A)') ' (EXHALE_main) outer pass ', it_diff,   &
                 ': REFUSED -- a row of this state is not finite, so it'// &
                 ' is not an iterate the alternation can continue from.'
            return
         else if (outer_ending .eq. outer_element_update_refused) then
            jfnk_info = 1
            write(*,'(A,I0,A)') ' (EXHALE_main) outer pass ', it_diff,   &
                 ': REFUSED -- the element composition relaxation found'// &
                 ' no admissible advance and its entry composition was'//  &
                 ' restored, so no further pass could differ from this'//  &
                 ' one.'
            return
         else if (outer_ending .eq. outer_no_progress) then
            jfnk_info = 1
            write(*,'(A,I0,A,I0,A)') ' (EXHALE_main) outer pass ',        &
                 it_diff, ': REFUSED -- the worst gated species row'//     &
                 ' has not fallen in ', n_no_fall, ' consecutive'//        &
                 ' passes; the alternation is not approaching a joint'//   &
                 ' fixed point at these step lengths.'
            write(*,'(A,ES10.3,A,ES10.3,A,A,A,I0)') '    it stands'//     &
                 ' at ', sp_meas, ' against ', sp_tol, ' (',              &
                 trim(sp_name), ') at cell ', sp_cell
            return
         endif
      enddo

      ! ---- THE PASS BUDGET ENDED WITHOUT AN ACCEPTED STATE ----
      ! The state handed back is the refreshed state of the last
      ! hydrodynamic solve. It is not accepted, and the entry that refuses
      ! it is named, so that the next budget, omega or movement bound is
      ! chosen against an equation and a cell.
      outer_ending = outer_pass_budget
      jfnk_info    = 1
      write(*,'(A,I0,A,I0,A)') ' (EXHALE_main) the stationary outer'//    &
           ' iteration spent its budget of ', pass_cap, ' passes'//       &
           ' without an accepted state: ', cert_now%n_failing,            &
           ' active equations refuse it.'
      if (ref_cell .gt. 0) then
         write(*,'(A,A,A,ES10.3,A,ES10.3,A,I0,A,F8.4)') '    worst'//     &
              ' refusing entry: ', trim(ref_name), ', measure ',          &
              ref_meas, ' against ', ref_tol, ' at cell ', ref_cell,      &
              ', r =', r(ref_cell)
      else
         write(*,'(A)') '    no measured row carries the refusal.'
      endif
      if (n_unjudged .gt. 0)                                              &
         write(*,'(A,I0,A,A,A,A)') '    and ', n_unjudged, ' active'//    &
              ' equations could not be judged on this state, the first'// &
              ' of them ', trim(unjudged_name), ': ', trim(unjudged_why)
      if (cert_now%n_no_chem_root .ne. 0)                                 &
         write(*,'(A,I0,A)') '    and ', cert_now%n_no_chem_root,         &
              ' cells hold a composition that is not a root of the'//     &
              ' chemical network.'
      call certification_report_write(cert_now,                           &
           'state of the last outer pass, NOT accepted')

      end subroutine steady_wind_with_element_diffusion

      ! ------------------------------------------------!

      subroutine write_diffusion_pass_profile(pass)
      ! Diagnostic (EXHALE_DIFFUSION_CHECK=1): append the element ratio and
      ! the mass flux of one outer pass to output/diffusion_pass_profiles.txt.
      ! Written so that successive passes can be compared directly -- the
      ! question a non-converging outer loop raises is whether it is cycling
      ! between two states or wandering, and that is a question about the
      ! profiles, not about the scalar drift.
      integer, intent(in) :: pass
      real*8, dimension(1-Ng:N+Ng) :: heh_l
      real*8 :: nh_l, x2_l
      logical, save :: dpp_opened = .false.
      integer :: j, uu
      logical :: first

      heh_l = element_ratio_HeH(f_sp)
      uu    = 773
      first = (pass .eq. 1) .and. (.not. dpp_opened)
      dpp_opened = .true.
      if (first) then
         open(unit=uu, file='output/diffusion_pass_profiles.txt',          &
              status='replace')
         write(uu,'(A)') '# outer-pass composition profiles '//            &
              '(EXHALE_DIFFUSION_CHECK=1)'
         write(uu,'(A)') '# columns: pass  r[Rp]  (He/H)/HeH  '//          &
              'r^2 rho v [n0 mH cm/s Rp^2]  T[K]  x_H2[2nH2/nH]  '//      &
              'n_H3p[cm^-3]  n_HeII[cm^-3]  n_e[cm^-3]'
      else
         open(unit=uu, file='output/diffusion_pass_profiles.txt',          &
              status='old', position='append')
      endif
      do j = 1, N
         ! The carrier partition belongs here too: with the molecular
         ! carriers transported it is the composition the outer loop is
         ! iterating on, and a scalar drift cannot say whether the passes
         ! cycle or wander.
         nh_l  = rho(j)*n0*(f_sp(j,1) + f_sp(j,2)                         &
               + 2.0d0*(f_sp(j,isp_H2) + f_sp(j,isp_H2p))                 &
               + 3.0d0*f_sp(j,isp_H3p) + f_sp(j,isp_HeHp))
         x2_l  = 2.0d0*rho(j)*n0*f_sp(j,isp_H2)/max(nh_l, 1.0d-99)
         write(uu,'(I5,8ES16.7)') pass, r(j), heh_l(j)/HeH,               &
              r(j)**2*rho(j)*v(j), T(j)*T0, x2_l,                         &
              rho(j)*n0*f_sp(j,isp_H3p), nheii(j)*n0, ne(j)*n0
      enddo
      close(uu)

      end subroutine write_diffusion_pass_profile

      ! ------------------------------------------------!

      subroutine allocate_state_vectors
      ! Allocate the grid-sized state, flux and diagnostic vectors of the
      ! marching loop, now that the number of computational cells N is known.
      ! The bounds 1-Ng:N+Ng (and the leading 3 of the conservative /
      ! primitive vectors) are the ones the declarations used to carry, and
      ! the zero start reproduces the static storage they came from.

      allocate(Rres(3,1-Ng:N+Ng))
      allocate(R_plm(3,1-Ng:N+Ng), R_weno(3,1-Ng:N+Ng),                   &
               u_probe(3,1-Ng:N+Ng))
      allocate(f_sp_test(1-Ng:N+Ng,n_species))
      allocate(heat0(1-Ng:N+Ng), cool0(1-Ng:N+Ng), npart0(1-Ng:N+Ng))
      allocate(mom(1-Ng:N+Ng))
      allocate(dt_loc(1-Ng:N+Ng))
      allocate(rho(1-Ng:N+Ng), v(1-Ng:N+Ng), E(1-Ng:N+Ng),                &
               p(1-Ng:N+Ng), T(1-Ng:N+Ng), cs(1-Ng:N+Ng))
      allocate(heat(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      allocate(Frho_elem(1-Ng:N+Ng))
      allocate(eta(1-Ng:N+Ng))
      allocate(nhi(1-Ng:N+Ng), nhii(1-Ng:N+Ng))
      allocate(nhei(1-Ng:N+Ng), nheii(1-Ng:N+Ng), nheiii(1-Ng:N+Ng))
      allocate(nheiTR(1-Ng:N+Ng))
      allocate(nm(1-Ng:N+Ng,n_mion))
      allocate(ne(1-Ng:N+Ng), n_tot(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species))
      allocate(u(3,1-Ng:N+Ng), u1(3,1-Ng:N+Ng), u2(3,1-Ng:N+Ng),          &
               u_old(3,1-Ng:N+Ng))
      allocate(W(3,1-Ng:N+Ng), WL(3,1-Ng:N+Ng), WR(3,1-Ng:N+Ng))
      allocate(u_th_old_csm(1-Ng:N+Ng), u_form_csm(1-Ng:N+Ng),            &
               u_form_old_csm(1-Ng:N+Ng), du_form_csm(1-Ng:N+Ng),         &
               T_csm_prev(1-Ng:N+Ng), rho_rec_csm(1-Ng:N+Ng))
      allocate(f_sp_csm(1-Ng:N+Ng,n_species))
      allocate(csm_inner_c(1-Ng:N+Ng), csm_inner_T(1-Ng:N+Ng))
      csm_inner_c = 0.0d0;  csm_inner_T = 0.0d0
      allocate(csm_dT_cell(1-Ng:N+Ng), csm_dc_cell(1-Ng:N+Ng),            &
               csm_at_rest(1-Ng:N+Ng))
      csm_dT_cell = 0.0d0;  csm_dc_cell = 0.0d0;  csm_at_rest = .false.
      call get_environment_variable('EXHALE_CSM_DEBUG', csm_dbg_env)
      csm_debug = (len_trim(csm_dbg_env) .gt. 0)
      ! The geometric decay of the pass sequence: the two increments the
      ! estimate is taken from, and the arrays the extrapolation saves.
      ! Only the cells the stopping test is taken over carry an increment
      ! (1:N); the save arrays hold the same range, because that is what
      ! the candidate moves.
      allocate(csm_geom_eT(1:N), csm_geom_eT_prev(1:N), csm_x_T_save(1:N))
      allocate(csm_geom_ec(1:N,n_species),                                &
               csm_geom_ec_prev(1:N,n_species),                           &
               csm_x_f_save(1:N,n_species))
      csm_geom_eT = 0.0d0;  csm_geom_eT_prev = 0.0d0
      csm_geom_ec = 0.0d0;  csm_geom_ec_prev = 0.0d0
      csm_x_T_save = 0.0d0; csm_x_f_save = 0.0d0
      csm_geom_have_prev  = .false.;  csm_geom_seq_broken = .false.
      csm_geom_theta_T = 0.0d0;  csm_geom_theta_c = 0.0d0
      csm_err_ok = .false.;  csm_err_T = 0.0d0;  csm_err_c = 0.0d0
      csm_x_pending = .false.
      csm_x_blocked = .false.;  csm_x_used    = 0
      csm_x_pred_dT = 0.0d0;  csm_x_pred_dc = 0.0d0
      csm_x_at_dT   = 0.0d0;  csm_x_at_dc   = 0.0d0
      csm_x_damp    = 0.0d0
      call get_environment_variable('EXHALE_CSM_GEOM_PASS', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_geom_first_pass
      call get_environment_variable('EXHALE_CSM_GEOM_MODEL', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_geom_model_tol
      call get_environment_variable('EXHALE_CSM_GEOM_THETA', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_geom_theta_max
      call get_environment_variable('EXHALE_CSM_GEOM_SAFETY', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_err_safety
      call get_environment_variable('EXHALE_CSM_GEOM_REFUSE', csm_x_env)
      csm_geom_refuse_all = (len_trim(csm_x_env) .gt. 0)
      call get_environment_variable('EXHALE_CSM_GEOM_DEBUG', csm_x_env)
      csm_geom_debug = (len_trim(csm_x_env) .gt. 0)
      ! The probe that MEASURES the accepted state against the fixed
      ! point. Its four arrays hold whole states, ghost cells included,
      ! because the pair it restores is the pair the accepted pass was
      ! taken from and the ghost cells are part of it.
      call get_environment_variable('EXHALE_CSM_ERR_PROBE', csm_x_env)
      csm_err_probe = (len_trim(csm_x_env) .gt. 0)
      csm_pass_cap = csm_max_pass
      if (csm_err_probe) then
         csm_pass_cap = csm_probe_pass_cap
         allocate(csm_probe_T_acc(1-Ng:N+Ng),                             &
                  csm_probe_T_entry(1-Ng:N+Ng))
         allocate(csm_probe_f_acc(1-Ng:N+Ng,n_species),                   &
                  csm_probe_f_entry(1-Ng:N+Ng,n_species))
         csm_probe_T_acc = 0.0d0;  csm_probe_T_entry = 0.0d0
         csm_probe_f_acc = 0.0d0;  csm_probe_f_entry = 0.0d0
      endif
      csm_probe_active = .false.;  csm_probe_replay = .false.
      csm_probe_pass_accepted = 0
      call get_environment_variable('EXHALE_CSM_EXTRAP', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0)                                     &
         csm_extrap_on = (trim(csm_x_env) .ne. '0')
      call get_environment_variable('EXHALE_CSM_EXTRAP_MAX', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_extrap_max
      call get_environment_variable('EXHALE_CSM_EXTRAP_REFUSE', csm_x_env)
      csm_extrap_refuse_all = (len_trim(csm_x_env) .gt. 0)
      call get_environment_variable('EXHALE_CSM_EXTRAP_REACH', csm_x_env)
      if (len_trim(csm_x_env) .gt. 0) read(csm_x_env,*) csm_x_reach
      call get_environment_variable('EXHALE_CSM_EXTRAP_STRICT', csm_x_env)
      csm_extrap_strict = (len_trim(csm_x_env) .gt. 0)
      ! The attempted-step controller's two working arrays: the operator the
      ! Runge-Kutta stages integrated, read back from the first stage, and
      ! the full-step result held aside by the step-doubling estimate.
      allocate(L_hyd(3,1-Ng:N+Ng), u_err_full(3,1-Ng:N+Ng))
      allocate(u_hyd(3,1-Ng:N+Ng));  u_hyd = 0.0d0
      allocate(dtloc_macro(1-Ng:N+Ng))
      allocate(f_sp_err_full(1-Ng:N+Ng,n_species))
      allocate(T_err_full(1-Ng:N+Ng), ne_err_full(1-Ng:N+Ng))
      L_hyd = 0.0d0;  u_err_full = 0.0d0;  dtloc_macro = 0.0d0
      f_sp_err_full = 0.0d0;  T_err_full = 0.0d0;  ne_err_full = 0.0d0
      allocate(dF(3,1-Ng:N+Ng), S(3,1-Ng:N+Ng))
      if (upmap_n .gt. 0) then
         allocate(u_um0(3,1-Ng:N+Ng), u_umA(3,1-Ng:N+Ng),                  &
                  u_umB(3,1-Ng:N+Ng), u_umC(3,1-Ng:N+Ng),                  &
                  u_umD(3,1-Ng:N+Ng), u_umSave(3,1-Ng:N+Ng),               &
                  R_um(3,1-Ng:N+Ng))
         allocate(f_sp_umSave(1-Ng:N+Ng,n_species))
         u_um0 = 0.0d0;  u_umA = 0.0d0;  u_umB = 0.0d0
         u_umC = 0.0d0;  u_umD = 0.0d0;  u_umSave = 0.0d0
         R_um  = 0.0d0;  f_sp_umSave = 0.0d0
      endif

      Rres      = 0.0d0
      f_sp_test = 0.0d0
      heat0     = 0.0d0
      cool0     = 0.0d0
      npart0    = 0.0d0
      mom       = 0.0d0
      dt_loc    = 0.0d0
      rho       = 0.0d0
      v         = 0.0d0
      E         = 0.0d0
      p         = 0.0d0
      T         = 0.0d0
      cs        = 0.0d0
      heat      = 0.0d0
      cool      = 0.0d0
      eta       = 0.0d0
      nhi       = 0.0d0
      nhii      = 0.0d0
      nhei      = 0.0d0
      nheii     = 0.0d0
      nheiii    = 0.0d0
      nheiTR    = 0.0d0
      nm        = 0.0d0
      ne        = 0.0d0
      n_tot     = 0.0d0
      f_sp      = 0.0d0
      u         = 0.0d0
      u1        = 0.0d0
      u2        = 0.0d0
      u_old     = 0.0d0
      W         = 0.0d0
      WL        = 0.0d0
      WR        = 0.0d0
      dF        = 0.0d0
      S         = 0.0d0

      end subroutine allocate_state_vectors

      ! End of program
      subroutine p54_flux_timeseries_append
      ! DIAGNOSTIC ONLY (P54/(AD), env EXHALE_P54_TS).  See the declaration
      ! block for what each column is.  face_flux is RK_integration's cached
      ! Riemann flux; the caller refreshed it with assemble_residual.
      integer :: i, jj
      real*8  :: ff(4), fc(4), sp103, sp110
      if (.not. p54_ts_open) then
         do i = 1,4
            p54_ts_face(i) = N
            do jj = 1, N
               if (r_edg(jj) .ge. p54_ts_r(i)) then
                  p54_ts_face(i) = jj
                  exit
               endif
            enddo
         enddo
         open(newunit=p54_ts_unit, file='output/p54_flux_ts.txt',          &
              status='replace', action='write')
         write(p54_ts_unit,'(A)') '# step  F_face(1.005 1.03 1.10 1.20)'// &
              '  F_cell(same cells)  spread_cell(r>=1.03) spread_cell(r>=1.10)'
         p54_ts_open = .true.
      endif
      do i = 1,4
         jj    = p54_ts_face(i)
         ff(i) = face_flux(1,jj)*r_edg(jj)*r_edg(jj)
         fc(i) = u(2,jj)*r(jj)*r(jj)
      enddo
      call p54_cell_spread(p54_ts_face(2), sp103)
      call p54_cell_spread(p54_ts_face(3), sp110)
      write(p54_ts_unit,'(1X,I8,1X,10(ES20.12,1X))') count,                &
           ff(1), ff(2), ff(3), ff(4), fc(1), fc(2), fc(3), fc(4),         &
           sp103, sp110
      flush(p54_ts_unit)
      end subroutine p54_flux_timeseries_append

      subroutine p54_cell_spread(ja, sp)
      ! Spread of the CELL-CENTRED rho v r^2 over [ja:N], the functional
      ! item (AD) quotes, evaluated from an arbitrary inner cell.
      integer, intent(in)  :: ja
      real*8,  intent(out) :: sp
      integer :: jj
      real*8  :: f, fmx, fmn, fsum
      sp = 0.0d0
      if (ja .gt. N) return
      fmx = -huge(1.0d0);  fmn = huge(1.0d0);  fsum = 0.0d0
      do jj = ja, N
         f    = u(2,jj)*r(jj)*r(jj)
         fmx  = max(fmx, f);  fmn = min(fmn, f);  fsum = fsum + f
      enddo
      sp = (fmx - fmn)/max(abs(fsum/dble(N - ja + 1)), 1.0d-30)
      end subroutine p54_cell_spread

      end program Hydro_ioniz

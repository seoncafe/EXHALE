      module steady_newton
      ! Steady-state residual in vector form, F(Y) = 0, for a direct
      ! (Newton / pseudo-transient-continuation) solve. Increment (ii)-2:
      ! pack/unpack + newton_residual. The Jacobian and the PTC driver are
      ! added in later increments.
      !
      ! Unknowns: the hydro conserved variables on the PHYSICAL cells
      ! j = 1..N (neq = 3*N). The ghost cells (1-Ng..0 lower, N+1..N+Ng
      ! upper) are NOT unknowns: Apply_BC fills them from the interior each
      ! residual evaluation (lower = fixed rho_bc / p = 1+dp_bc with the
      ! v-valve max(v_1,0); upper = zero-gradient / WENO3 extrapolation).
      !
      ! Ionization/temperature are eliminated LOCALLY inside the residual:
      ! given Y, ioniz_eq re-solves the equilibrium fractions and returns the
      ! consistent heat/cool, so F(Y) is a pure-hydro residual of size 3*N
      ! (the "method A" choice from ATES_sundials_solver_plan.md).

      use global_parameters
      use species_table,  only: n_mion, isp_H2, isp_H2p, isp_H3p,   &
                                isp_HeHp, isp_HI, isp_HII,           &
                                isp_OH, isp_H2O, isp_CO, iel_O, iel_C
      use Conversion,      only: U_to_W, U_to_W_interior
      use composition,     only: get_species_densities, comp_T_from_p
      use BC_Apply,        only: Apply_BC
      use ionization_equilibrium, only: ioniz_eq, ioniz_eq_ledger,      &
                 set_ioniz_eq_sweep_state_kind, ieq_state_marching,      &
                 ieq_state_steady_iterate, ieq_state_steady_candidate,   &
                 finite_real, keep_background_of_adopted_state,          &
                 keep_background_of_best_iterate,                        &
                 adopt_background_of_best_iterate,                       &
                 install_background_of_adopted_state,                  &
                 ieq_sweep_ledger_last, xuv_self_field_passes,        &
                 ieq_report_cell,                                      &
                 molecular_decay_rate_of_the_cell,                     &
                 set_molecular_decay_cells
      use excited_hydrogen,       only: excited_H_update
      ! The face states and interface fluxes the flux assembly last built,
      ! and the count of reconstructions the positivity limiter scaled back.
      ! Read by the jump scan only, to say WHICH operator of the pipeline
      ! state -> ghosts -> face states -> flux -> row is the first that
      ! jumps; the counter is saved and restored around the scan's own
      ! reconstruction so a measurement never enters the run's ledger.
      use RK_integration,         only: face_flux, face_p
      use Reconstruction_step,    only: n_faces_positivity_limited
      ! The two counters that say whether a DISCRETE BRANCH of the flux
      ! assembly fired: the faces where no admissible star state existed
      ! and the Roe flux gave way to HLLE, and the face states the
      ! positivity limiter scaled back. A branch that flips between two
      ! nearby states is a step of the map; one that never fires is not.
      use Numerical_Fluxes,       only: n_faces_roe_hlle
      ! The census of the composition's inner cell solve: how many the
      ! analytic-Jacobian Newton handled and how many fell through to
      ! MINPACK hybrd1. Read by the additivity hook only, to say what a
      ! stopping tolerance costs the inner solve.
      use newton_solver,          only: nt_calls, nt_fallback
      use steady_residual_mod,    only: assemble_residual,             &
                                        reconstruction_continuation_rhs, &
                                        face_mass_flux_of_state,       &
                                        residual_row_scale,            &
                                        mass_row_rounding_floor,       &
                                        state_scales_of_cell,          &
                                        flux_spread_of_state,          &
                                        flux_spread_above_radius,      &
                                        steady_gates_met,              &
                                        n_cells_without_chemical_root, &
                                        write_residual_breakdown
      use binary_element_diffusion, only: element_transport_residual,  &
                                        element_mass_fractions,         &
                                        project_element_mass_fractions, &
                                        mixture_mass_sum,               &
                                        trace_row_terms_diag,           &
                                        trace_row_terms_from
      use species_table,          only: n_melem, melem_name
      ! THE ONE STOICHIOMETRIC MAP OF THE NUCLEI OF A CELL (item N4a). The
      ! shared element constraint of the step, its derivatives and the
      ! element totals the write-back restores are all read from here, so
      ! that the step, the feasibility screen and the composition it writes
      ! cannot disagree about how many nuclei a species holds.
      use element_inventory,      only: element_inventory_cell,          &
                                        element_inventory_of_cell,       &
                                        element_constraint_derivatives,  &
                                        nuclei_per_particle,             &
                                        inventory_feasible_slack,        &
                                        n_inv_element, ien_H, ien_He,    &
                                        ien_O, ien_C, inv_element_name
      use diffusive_photochemistry, only: carrier_steady_residual,     &
                                        carrier_column_scale,          &
                                        carrier_row_term_scale,        &
                                        carrier_headroom,              &
                                        carrier_headroom_known,        &
                                        carrier_species_index,         &
                                        carrier_solved, carrier_name,  &
                                        n_carrier, n_carrier_max, ic_H2, &
                                        ic_OH, ic_H2O, ic_CO, ic_Hp,    &
                                        carrier_module_state,          &
                                        save_carrier_module_state,     &
                                        restore_carrier_module_state,  &
                                        carrier_state,                 &
                                        carrier_source,                &
                                        carrier_mass_amu,              &
                                        carrier_advective_divergence,  &
                                        carrier_diffusion_coefficient
      ! The acceptance contexts of docs/a2_certification_contract_20260906.md
      ! section 3. The probe and the trial decisions are taken there, on the
      ! facts one residual evaluation records, so that this solver and the
      ! marching stop read ONE statement of what makes a state usable.
      use certification,          only: cert_evaluation_facts,          &
                                        trial_state_is_admissible,      &
                                        probe_direction_is_usable,      &
                                        cert_report, cert_context_stationary, &
                                        certification_evaluate,         &
                                        certification_report_write,     &
                                        certification_last_report,      &
                                        cert_evaluated, cert_unavailable, &
                                        certification_row_measure,       &
                                        cert_scale_floor,                &
                                        cert_tol_element_at,             &
                                        cert_tol_carrier_at,             &
                                        mass_row_cell_verdict,           &
                                        cert_tol_momentum,               &
                                        cert_tol_energy,                 &
                                        cert_regime_wind_r,              &
                                        certification_species_row_gate

      implicit none
      private
      public :: neq_newton, pack_U, unpack_U, newton_residual,         &
                pack_species_rows, write_species_rows_into_composition, &
                cell_state_scales, cell_row_scales,                     &
                merit_row_scale_from_certification,                     &
                eval_residual, frozen_residual, build_banded_jac,      &
                band_matvec, band_matvec_transpose, kl_jac, ku_jac,    &
                ncolor_jac, solve_steady_ptc,                          &
                jv_product, solve_steady_jfnk, set_base_fix,           &
                gate_rnorm_accepted, gate_fspread_accepted,            &
                set_transported_species_rows,                          &
                carrier_unknown_on, n_species_rows,                    &
                species_row_carrier_index,                             &
                distance_from_certification, resid_tol_of_solve,        &
                hydrodynamic_distance_from_certification_by_cell,        &
                increment_of_what_the_residual_reads,                   &
                eq_channel_particle_count, eq_channel_h2_caloric,       &
                eq_channel_radiative,                                   &
                species_unknowns_outside_their_bounds,                   &
                hold_the_step_where_it_leaves_the_species_box,            &
                species_box_set_for_test, species_box_clear_for_test,     &
                probe_length_of_the_jacobian_action,                      &
                freeze_species_unknown_box,                              &
                read_species_unknown_space_controls, nvar_jac,           &
                pgmres, dogleg_step,                                      &
                cauchy_length_along_the_banded_gradient,                  &
                predicted_model_decrease,                                 &
                linear_operator_set_for_test,                             &
                linear_operator_clear_for_test, lin_test_dense,           &
                lin_test_dense_preconditioner_fails,                      &
                lin_test_dense_action_drifts, lin_test_drift_after,       &
                gm_tolerance_reached, gm_subspace_exhausted,              &
                gm_no_direction_sampled, gm_jacobian_action_unusable,     &
                gm_preconditioner_failed, gm_reduced_system_singular,     &
                gm_nonfinite, gm_outcome_text,                            &
                initial_trust_region_radius,                              &
                tr_radius_from_the_krylov_step,                           &
                tr_radius_from_the_gradient,                              &
                tr_radius_no_feasible_descent,                            &
                ray_slope_verdict, ray_factor, ray_floor_factor,          &
                ratio_sound_band, eta_accept,                             &
                trust_region_radius_absolute_cap,                         &
                grow_the_trust_region_ceiling,                            &
                cauchy_leg_min_fraction, tr_step_taken,                   &
                cauchy_leg_image_share_max, cauchy_leg_step_share_max,    &
                cauchy_leg_shares_of_step_and_image,                      &
                cauchy_leg_is_image_without_step,                         &
                cauchy_leg_is_dropped_on_its_length,                      &
                cauchy_leg_admitted_by_image_alone,                       &
                tr_ratio_below_acceptance,                                &
                tr_ray_probe_inadmissible,                                &
                tr_ray_slopes_disagree_sign,                              &
                tr_ray_slopes_disagree_size,                              &
                tr_ray_below_the_noise_floor,                             &
                jfnk_outer_iteration_cap,                                 &
                stationary_rows_of_the_returned_state,                    &
                row_maxima_of_the_certification, replay_distance,          &
                certification_row_class_name,                             &
                front_of_the_steepest_fraction,                           &
                resolved_gradient_length,                                 &
                row_condition_number, row_rounding_floor_in_the_measure,  &
                species_row_base_equation,                                &
                element_unknown_column_scale,                             &
                ritz_values_of_a_matrix_free_operator,                    &
                column_split_against_the_band,                            &
                unknown_mask_of_the_row_class,                            &
                row_class_shares_of_a_vector,                             &
                precond_spectrum_on, krylov_size_scan_on,                 &
                band_difference_on, gm_residual_history_on,               &
                front_row_on, the_binding_species_row,                    &
                terms_of_a_transported_species_row,                       &
                species_row_term_set, species_row_from_its_terms,         &
                attributed_jacobian_entry,                                &
                value_of_the_species_unknown,                             &
                n_ritz_step,                                              &
                element_column_scale_floor,                               &
                base_row_control_volume_balance,                          &
                base_row_reservoir_condition,                             &
                measure_the_closure_map_at_the_base,                      &
                outer_residual_closure_target,                            &
                closure_amplification_at_the_base,                        &
                equilibrate_the_scaled_columns,                           &
                column_equilibration_on,                                  &
                equilibrate_the_preconditioner,                           &
                preconditioner_equilibrated,                              &
                precon_row_scale, precon_col_scale,                       &
                banded_preconditioner_solve,                              &
                freeze_element_constraint_rows,                           &
                element_constraint_rows_known,                            &
                element_constraint_of_the_cell,                           &
                element_constraint_budget_of_the_cell,                    &
                element_constraint_gradient_of_the_cell,                  &
                project_out_of_the_active_element_constraints,            &
                largest_step_inside_the_element_constraints,              &
                restore_the_element_budget,                               &
                n_element_constraint_rows,                                &
                element_constraint_restoration_trigger,                   &
                element_budget_violation_of_composition,                  &
                largest_step_inside_the_species_box,                      &
                zero_the_blocked_components_of,                           &
                fix_active_species_bounds, n_active_element_rows,          &
                release_active_species_bounds,                            &
                read_trust_region_restart_controls,                        &
                trust_region_restart_is_due,                               &
                restart_the_trust_region_state, tr_reset_name,             &
                tr_reset_wanted, n_tr_reset_kind,                          &
                tr_reset_radius_and_ceiling, tr_reset_pseudo_transient,    &
                tr_reset_closure_map, tr_reset_acceptance_memory,          &
                tr_reset_back_to_best, tr_restart_forced_at,               &
                tr_restart_without_improvement, n_tr_restart,              &
                gm_restart_cycles, gm_restart_cycles_wanted,               &
                gm_restart_cycles_from, pgmres_with_restarts,              &
                gm_stopped_on_the_trust_ball,                              &
                krylov_truncated_on_the_trust_ball,                        &
                model_row_equilibration_on, model_row_equilibration,       &
                model_rows_equilibrated,                                   &
                unit_infinity_norm_row_scaling_of_the_band,                &
                gm_step_by_its_true_residual, gm_true_residual_turned,     &
                gm_true_residual_first, gm_true_residual_stride,           &
                gm_reorthogonalize, gm_measure_orthogonality,              &
                gm_orthogonality_loss_cycle,                               &
                gm_reorthogonalization_loss_level,                         &
                gm_verify_the_arnoldi_image,                               &
                jv_additivity_on, jv_additivity_here,                      &
                resid_jump_scan_on, resid_jump_scan_here,                  &
                scan_the_residual_along_a_direction,                       &
                median_of_the_magnitudes, steps_in_a_sampled_row,          &
                jv_probe_on_the_column_scales, state_column_scale,         &
                jv_probe_arc_scale,                                        &
                jacobian_action_of_direction,                              &
                probe_step_on_the_column_scales,                           &
                jv_probe_length_factor,                                    &
                measure_the_additivity_of_the_jacobian_action,             &
                the_point_where_the_step_leaves_the_ball,                  &
                name_of_unknown

      ! Banded-Jacobian geometry. WENO3 reconstruction reaches +/-2 cells,
      ! and within a cell all 3 hydro variables couple, so in the flat
      ! ordering Y(3*(j-1)+k) a column couples to rows within
      ! |row-col| <= 3*2 + 2 = 8. Half-bandwidths kl = ku = 8; a graph
      ! coloring with stride ncolor = kl+ku+1 = 17 makes same-color columns'
      ! row supports disjoint, so the FROZEN-radiation residual (strictly
      ! local) is captured exactly by 17 finite-difference probes.
      ! The two gate values of the state the LAST steady solve returned, kept
      ! so that the caller can re-measure them on the state it finally writes
      ! and prove the two are the same (EXHALE_main, after write_output).
      ! Nothing between the solve's return and the write touches u -- the
      ! refresh at EXHALE_main lines 1599-1607 recomputes rho, v, p FROM u and
      ! only moves the composition -- so the check is an assertion, not a
      ! correction; it exists because "the state the user reads" and "the state
      ! the gate accepted" must be the same object and nothing was saying so
      ! (section 133.6).
      real*8 :: gate_rnorm_accepted   = -1.0d0
      real*8 :: gate_fspread_accepted = -1.0d0

      ! THE SPECIES ROWS OF THE STATIONARY SYSTEM (B5).
      !
      ! Every transported balance the configuration activates is a ROW of
      ! this system and its transported quantity is an UNKNOWN of the cell:
      ! the stationary equation div(F_s + Phi_s) = P_s - L_s, on the
      ! independent space of that configuration. The local equilibrium
      ! partition of an element among its stages stays a closure of the cell
      ! solve; what is carried here is the balance OUTSIDE that solve.
      !
      ! The registry below is the only place the mapping row -> physics is
      ! written. Slot i of cell j is nvar_jac*(j-1) + 3 + i.
      ! srow_carrier       a molecular or proton carrier; srow_idx is its
      !                    carrier index, srow_isp its f_sp column, and the
      !                    unknown is n_c/n0 = f_sp*rho.
      ! srow_element_he    the helium/hydrogen partition; the unknown is the
      !                    helium MASS FRACTION X = rho_He/rho, already
      !                    dimensionless and O(0.1).
      ! srow_element_trace one transported trace element; srow_idx is its
      !                    index in the metal element table and the unknown
      !                    is that element's mass fraction.
      !
      ! An element unknown is not one column of f_sp: it is an element
      ! TOTAL, and the map from it back to a species vector is
      ! project_element_mass_fractions, the same map the advective stages
      ! use.  The element rows are therefore registered BEFORE the carriers,
      ! and written back before them, so that a carrier the Newton has set
      ! is not rescaled afterwards by the element projection.
      integer, parameter :: srow_carrier       = 1
      integer, parameter :: srow_element_he    = 2
      integer, parameter :: srow_element_trace = 3
      integer, parameter :: n_species_row_max = n_carrier_max + 1 + n_melem
      integer :: nspec_row = 0
      integer :: srow_kind(n_species_row_max) = 0
      integer :: srow_idx (n_species_row_max) = 0
      integer :: srow_isp (n_species_row_max) = 0
      ! The rows and the scales of the elemental balances of the state the
      ! last residual evaluation formed, kept for the row scaling exactly as
      ! diffusive_photochemistry keeps the carriers'.
      real*8, allocatable :: erow_he(:), escale_he(:)
      real*8, allocatable :: erow_tr(:,:), escale_tr(:,:)
      ! Which trace elements the state carries a balance for.  An element
      ! the run holds none of has no equation: its row is the statement that
      ! it stays at none, which is a unit diagonal and not an empty column.
      logical :: erow_tr_carried(n_melem) = .false.
      ! THE CARRIER ROWS OF THE STATE THE LAST EVALUATION FORMED and the sum
      ! of the magnitudes of each row's own terms, kept beside the elemental
      ! pair above and for the same reason: a row and the sum of its terms
      ! are what say how much of the row survives the cancellation of those
      ! terms (row_condition_number), and that sum is also the scale the
      ! certification divides the row by, so the two must come from ONE
      ! assembly of ONE state.
      real*8, allocatable :: crow_res_last(:,:), crow_terms_last(:,:)
      ! THE DIRICHLET RESERVOIR OF THE ELEMENT OPERATOR, one value per
      ! element row, recorded from the state the unknown vector is packed
      ! from.
      !
      ! Cell 1 carries no transport balance in either arm: it is the
      ! reservoir the operator holds at the input composition
      ! (`solve_mass_fraction` pins Xhe(1-Ng:1) = X_base, and the trace solve
      ! pins the base mixing ratio at its own entry value), so
      ! element_transport_residual fills its rows over j = 2..N and reports
      ! cell 1 as zero against the floor.  As a DIAGNOSTIC that is the right
      ! statement.  As a ROW of a stationary system it is an empty row: the
      ! residual of the base element unknown would be zero for every Y, the
      ! unknown would be unconstrained and the banded model singular -- the
      ! same defect the trace branch of eval_residual already names for an
      ! element with no reservoir.  The equation of that unknown is therefore
      ! the boundary condition the operator states, X(1) = X_reservoir, and
      ! this is the value it is anchored at.  It is the state's own base
      ! value because that state has been relaxed by the operator that pins
      ! it, so the two are the same number and the solve is told not to move
      ! the base composition, which is what the operator says.
      real*8  :: srow_base_value(n_species_row_max) = 0.0d0
      ! WHICH EQUATION A SPECIES ROW CARRIES AT THE FIRST PHYSICAL CELL.
      ! The two possibilities are named because the residual and the row
      ! scaling both have to know, and a row whose equation is one thing in
      ! the residual and another in the scaling is nondimensionalized by a
      ! quantity its equation does not contain.
      integer, parameter :: base_row_control_volume_balance = 1
      integer, parameter :: base_row_reservoir_condition    = 2
      ! Whether the emptiness of the banded model has been stated for this
      ! registry (banded_model_conditioning): a property of the system posed,
      ! so it is said once and not once per outer iteration.
      logical :: band_emptiness_stated = .false.

      ! THE ORDERED DIAGNOSTIC OF ONE ELEMENT-ROW STEP, off unless
      ! EXHALE_ELEM_DIAG is set. At the selected outer iterations it states,
      ! in this order: the extremes of the two scalings with the quantity
      ! each belongs to; the columns the banded model does not sample; the
      ! status, the pivots and a condition estimate of the band
      ! factorization; the linear outcome of the Krylov leg with the TRUE
      ! residual of the step it returned; the unknowns the step holds at a
      ! bound beside the samples that were refused; the two legs of the
      ! dogleg with their predicted decreases before any trial; and the
      ! share each row kind holds of the convergence measure and of the
      ! merit.
      !
      ! IT DECIDES NOTHING. The only arithmetic it adds is one Jacobian
      ! product for the true residual of the Krylov leg, taken with the
      ! solve's refusal statistics held aside (solve_refusal_statistics), so
      ! that the step the solve takes is the step it takes with the hook
      ! off. TEST ONLY.
      !
      ! EXHALE_ELEM_DIAG=1 selects outer iterations 1, 2 and 60; an integer
      ! above 1 replaces the third of those.
      logical :: elem_diag_on    = .false.
      integer :: elem_diag_third = 60
      ! Whether THIS outer iteration is one of the selected ones.
      logical :: elem_diag_here  = .false.
      ! WHETHER EVERY TRIAL OF THIS OUTER ITERATION STATES ITS OWN MODEL.
      ! The predicted decrease is linear in the step, so on the SAME iterate
      ! a step cut by a factor must carry a decrease cut by the same factor
      ! up to the quadratic term; where it does not, the step, its image and
      ! the two terms of the model have to be read side by side, at each
      ! trial and not once per iteration. True from elem_diag_third onward,
      ! which is where the steps of a stalling solve are short enough for
      ! the question to arise. Costs one Jacobian product at each trial.
      logical :: elem_diag_every_trial = .false.

      ! WHERE THE ELEMENT-ROW TERM DIAGNOSTIC STARTS, and whether the helium
      ! row prints its own terms.  The trace hook of the element operator
      ! prints from a cell the caller names (trace_row_terms_from); the
      ! helium row is formed by the binary composition equation, which has
      ! no such hook, so its terms are split and printed here
      ! (write_element_row_terms).  EXHALE_ELEM_DIAG_FROM names the first
      ! cell of the window; with nothing set it is the outermost six cells,
      ! which is where the diagnostic was first needed.  Read only by a
      ! print.
      integer :: elem_diag_cell_from = 0
      ! Armed for ONE evaluation of the iterate, by the same block that arms
      ! the trace hook, and disarmed again, so no probe state prints.
      logical :: element_row_terms_diag = .false.

      ! HOW MANY OUTER ITERATIONS THE STATIONARY SOLVE MAY TAKE, when the
      ! environment names a number instead of the caller. The cap is the
      ! caller's argument; EXHALE_JFNK_MAXIT overrides it so that a bounded
      ! diagnostic of a solve that does not converge stops where it is asked
      ! to instead of being cut off by a wall clock, which leaves the last
      ! iteration of the log half written and the statistics over a length
      ! nobody chose. Zero or absent means the caller's value. TEST ONLY:
      ! nothing about the step depends on it.
      integer :: jfnk_maxit_wanted = 0

      ! WHETHER THE FINITE-DIFFERENCE ACTION IS A LINEAR MAP OF THE
      ! DIRECTION, off unless EXHALE_JV_ADDITIVITY is set. The Arnoldi
      ! relation, the predicted decrease read from it and the least-squares
      ! problem GMRES solves are all relations of a LINEAR operator, while
      ! the action jv_product returns is a secant of the residual over an
      ! arc of fixed length: homogeneous in the direction by construction
      ! but not additive across directions unless the residual is affine on
      ! the arc and the arc is the same one on both sides. The hook
      ! measures the departure and separates its three possible sources:
      ! the curvature of a smooth residual over the arc (the defect falls
      ! with the arc as the arc), rounding or a non-smooth residual (it
      ! rises as the reciprocal of the arc), and a direction-dependent
      ! feasible-set decision inside jv_product itself (it does neither).
      ! TEST ONLY: it decides nothing and every product it takes is thrown
      ! away.
      ! WHETHER THE PROBE DISPLACES EVERY UNKNOWN BY THE SAME FRACTION OF
      ! ITS OWN SCALE (probe_step_on_the_column_scales) instead of fixing
      ! the LENGTH of the displacement in the unknowns themselves
      ! (probe_length_of_the_jacobian_action). Default off, arm
      ! EXHALE_JV_COLUMN_SCALE=1.
      logical :: jv_probe_on_the_column_scales = .false.
      logical :: jv_additivity_on   = .false.
      ! Whether THIS outer iteration is one the hook speaks at.
      logical :: jv_additivity_here = .false.
      ! A MULTIPLIER OF THE PROBE ARC, exactly one everywhere but inside the
      ! hook's length scan, where the same three products are taken at a
      ! tenth and at ten times the standard arc. At one the arithmetic of
      ! probe_length_of_the_jacobian_action is untouched.
      real*8  :: jv_probe_length_factor = 1.0d0
      ! A STANDING MULTIPLIER OF THE PROBE ARC (EXHALE_JV_PROBE_ARC=<x>,
      ! default 1). The arc sqrt(epsilon)(1 + ||Y||) balances the rounding
      ! of the residual against the curvature of it only where the
      ! residual's spread is its own last bit. MEASURED (item N31, atomic
      ! element reload): along the PRECONDITIONED Krylov directions the
      ! residual carries a floor five to seven decades above that, the
      ! additivity defect is that floor divided by the increment the probe
      ! produces, and it falls as the reciprocal of the arc over two
      ! decades either side of the standard one with no crossover, so a
      ! longer arc buys accuracy in proportion until the curvature term
      ! meets it. This is the arm that measures where that is.
      real*8  :: jv_probe_arc_scale = 1.0d0
      ! WHAT THE LAST PROBE OF jv_product ACTUALLY DID, so that a defect can
      ! be attributed to a decision of the feasible set rather than to the
      ! residual. Written by every product; read only by a print.
      real*8  :: jv_probe_step_nominal_last = 0.0d0
      real*8  :: jv_probe_step_last         = 0.0d0
      integer :: jv_probe_blocked_last      = 0
      logical :: jv_probe_backward_last     = .false.

      ! IS THE RESIDUAL CONTINUOUS ALONG THE PROBE ARC, AND IF NOT, WHICH
      ! OPERATOR BREAKS IT (EXHALE_RESID_JUMP_SCAN=1). A second difference
      ! that does not fall with the square of the spacing but sits on a
      ! floor says F itself steps somewhere inside the arc; where the step
      ! is, and which piece of the pipeline
      !     unknowns -> ghosts -> face states -> interface flux -> row
      ! takes it, is a matter of sampling F densely along one direction and
      ! reading the pipeline on either side of the step. Default off; the
      ! scan takes about seventy residual evaluations at the one outer
      ! iteration it speaks at and decides nothing.
      logical :: resid_jump_scan_on   = .false.
      ! Whether THIS outer iteration is one the scan speaks at.
      logical :: resid_jump_scan_here = .false.
      ! THE STATE THE FLUX ASSEMBLY SAW, kept by eval_residual while the
      ! scan asks for it: the conserved variables INCLUDING the ghost cells
      ! Apply_BC wrote, the temperature the sweep last formed and the
      ! particle count n_tot + n_e the residual divides the pressure by.
      ! The ghosts are the point of it: they are not unknowns, they are a
      ! construction of the boundary condition, and no caller of
      ! eval_residual can otherwise see what that construction produced.
      logical :: resid_capture_operator_state = .false.
      real*8, dimension(:,:), allocatable :: captured_state_with_ghosts
      real*8, dimension(:),   allocatable :: captured_temperature
      real*8, dimension(:),   allocatable :: captured_particle_count

      ! WHAT ONE RESIDUAL EVALUATION LEAVES IN THE MODULES IT CALLS, AND A
      ! LATER READER CONSUMES.
      !
      ! F(Y) is a function of Y. The evaluation nevertheless writes the
      ! elemental rows and their scales above, the carrier measures below,
      ! and the carrier arrays of diffusive_photochemistry, because those
      ! are the products of the operators the species rows come from; the
      ! Newton ROW SCALE, the acceptance gate and the certification then
      ! read them back. An evaluation at a point the caller DISCARDS -- a
      ! Jacobian column, a Krylov product, a point on a trust-region ray, a
      ! self-test -- must therefore put them back, or the state the solve
      ! does hold is scaled and judged by the numbers of a state nobody
      ! keeps.
      !
      ! MEASURED, which is why this exists: with EXHALE_SPECIES_JAC_TEST=1
      ! the hot Uranus element solve agreed with the untested run to every
      ! printed digit for five outer iterations and was on a different path
      ! by iteration 13, because the self-test's probe evaluations at
      ! Y + eps r left that state's element row scales behind for the first
      ! outer iteration to scale itself by. It is the same statement
      ! save_carrier_module_state makes for the certification: measuring a
      ! balance may not be an operation on the state.
      type :: residual_evaluation_products
         real*8, allocatable :: erow_he(:), escale_he(:)
         real*8, allocatable :: erow_tr(:,:), escale_tr(:,:)
         real*8, allocatable :: crow_res(:,:), crow_terms(:,:)
         logical :: erow_tr_carried(n_melem) = .false.
         real*8  :: carrier_relnorm  = 0.0d0
         real*8  :: carrier_cellmax  = 0.0d0
         integer :: carrier_cell     = 0
         integer :: n_headroom_last    = 0
         real*8  :: v_headroom_iterate = 0.0d0
         integer :: n_eq_sweeps_last   = 0
         integer :: n_no_chem_root_last = 0
         type(carrier_module_state) :: carrier
      end type residual_evaluation_products

      ! THE REFUSAL AND COST STATISTICS OF A SOLVE, held aside around a
      ! DIAGNOSTIC.
      !
      ! These belong to the solve and a probe's refusal belongs in them: a
      ! Newton column the solve was refused is a column it did not get,
      ! whoever asked for it. What may not write them is a SELF-TEST, which
      ! is not part of the solve at all; a run reporting different refusal
      ! counts because a diagnostic was switched on is reporting the
      ! diagnostic.
      type :: solve_refusal_statistics
         integer :: n_trial_headroom      = 0
         integer :: n_trial_no_chem_root  = 0
         integer :: no_chem_root_cell_max = 0
         real*8  :: no_chem_root_res_max  = 0.0d0
         integer :: n_eval_refusal(0:6)   = 0
      end type solve_refusal_statistics

      ! WITH THE SPECIES ROWS THE GEOMETRY CHANGES, so these are set once per
      ! solve by set_transported_species_rows rather than fixed at compile
      ! time.
      !
      ! WENO3 reaches +/-2 cells and within a cell every variable couples, so
      ! in the flat ordering Y(nvar*(j-1)+k) a column couples to rows within
      !   |row - col| <= 2*nvar + (nvar - 1) = 3*nvar - 1 ,
      ! and a graph coloring with stride ncolor = kl + ku + 1 = 6*nvar - 1
      ! makes same-color columns' row supports disjoint, so the
      ! FROZEN-radiation residual (strictly local) is captured exactly by
      ! that many finite-difference probes. The formula reproduces the two
      ! geometries the code carried before it: nvar = 3 -> 8 and 17,
      ! nvar = 4 -> 11 and 23.
      !
      ! A species row itself reaches only +/-1: its face fluxes are formed
      ! from the two neighboring cells and the second-order part of the
      ! reconstruction is not differentiated, so the hydro rows still set
      ! the bandwidth and each extra variable only widens the stride.
      integer :: nvar_jac = 3
      integer :: kl_jac = 8, ku_jac = 8
      integer :: ncolor_jac = 17
      ! J^T J of a band of half-width kl_jac has half-width 2*kl_jac, and it
      ! is what the damped Gauss-Newton (Levenberg-Marquardt) step below
      ! factors when the Newton direction ascends the merit.
      integer :: kl_jac_normal = 16
      integer :: ku_jac_normal = 16
      ! LAPACK general-band leading dimension for that matrix.
      integer :: ldab_normal   = 49

      ! HOW A PROBE STATE THAT CANNOT BE DESCRIBED IS ANSWERED
      ! (docs/Update_EXHALE_stage1.md section 121).
      !
      ! The Jacobian columns and the Krylov products are directional
      ! derivatives OF THE RESIDUAL AT Y. The point Y + h*v at which they are
      ! sampled is not physics: h is ours to choose, and v is a Krylov
      ! direction with no physical meaning at all. So when the sample point
      ! leaves the set of states the code can describe -- a balance row of the
      ! chemical network goes to NaN or overflows, and the equilibrium sweep
      ! reports it -- the answer is to sample CLOSER TO Y, not to accept the
      ! non-finite residual and not to declare the derivative undefined. The
      ! step is halved and the probe repeated, at most this many times; over
      ! six halvings the finite-difference step falls by 64 and stays far
      ! above the round-off floor of the residual.
      !
      ! If even the smallest of those steps is inadmissible, the boundary of
      ! the describable states runs through Y itself along v, and there is no
      ! derivative to sample. Then the derivative is NOT invented: the Krylov
      ! cycle stops at the directions it did build (a shorter GMRES cycle is
      ! still a valid approximate solve, and the line search judges the step
      ! it produces), and a Jacobian color that cannot be probed leaves its
      ! columns at zero, which degrades the PRECONDITIONER of those unknowns
      ! to plain relaxation and costs iterations, not correctness. Both
      ! events are counted and printed.
      integer, parameter :: n_probe_step_halvings = 6

      ! A SECOND WAY A PROBE STATE CANNOT BE DESCRIBED, and the same answer.
      ! The residual's energy row contains heat - cool, and both are evaluated
      ! on the composition the equilibrium sweep accepted. When that
      ! composition is not a root of the network -- acceptance class 4, kept
      ! only by the relaxation amnesty, or class 6, whose candidate was so far
      ! from a root that the cell kept the composition it entered the sweep
      ! with -- the residual is a number computed on
      ! a state the chemistry does not describe, and its value carries no
      ! information about how far the hydro is from steady.
      !
      ! WHICH CLASSES ARE ROOTS is the producer's definition, and there is one
      ! reading of it: n_cells_without_chemical_root (steady_residual.f90)
      ! counts classes 4 and 6 and nothing else. Classes 1, 2, 3 and 5 are all
      ! of the requested equations, certified by the same residual test; a
      ! class-5 state is a root the constrained element-conserving
      ! continuation FOUND, and the route to a root is not a weaker
      ! certificate of it.
      !
      ! WHAT THIS COUNT IS USED FOR HERE is the LOCAL VALIDITY OF THE
      ! ELIMINATED CHEMISTRY of a trial, not stationarity: the algebraic
      ! closure the energy row rests on either holds in every cell or it does
      ! not. The test is therefore a COMPARISON, not a threshold: a trial is
      ! refused when it leaves MORE cells without a chemical root than the
      ! iterate it is being measured against. That is the statement the line
      ! search actually needs -- it compares a trial's merit with the
      ! iterate's, and the comparison means something only while the trial's
      ! chemistry is no worse -- and at an iterate whose every cell has a root
      ! it is the same thing as refusing any such trial at all. No stationary
      ! tolerance enters here; that belongs to the certification of a final
      ! state (steady_gates_met).
      !
      ! WHY NOT THE ABSOLUTE FORM. Refusing every trial with a rootless
      ! cell was written first and measured: on the hot Uranus hand-off the
      ! iterate ITSELF carries 11 such cells, so all 17 Jacobian colors and
      ! every Krylov direction were refused, the solve had no Newton model at
      ! all (`||grad merit|| = 0`) and aborted at iteration 1. A rule that
      ! rejects the neighborhood of the point it is standing on cannot be
      ! used to leave that point.
      !
      ! It is applied only where `admissible` is asked for, so the marching
      ! path and the current iterate are untouched.
      !
      ! Grounds (docs/Update_EXHALE_stage1.md section 138; measurement P49). On the
      ! hot Uranus hand-off the solve spent 66 iterations driving the layer at
      ! 1.19-1.23 R_p from 3.5e-14 to 7.7e-16 g cm^-3 and from 900 K to
      ! 65,000 K while ||R|| improved from 8.2e-2 to 2.8e-2. What held that
      ! state up was cell 259, accepted class 4 with a reaction residual of
      ! 2.4e-3, whose H3+ infrared cooling -- 93% of the blocking energy row --
      ! came from an H3+ population that cannot exist at 65,000 K. The
      ! residual was measuring a composition, not a wind.

      ! Frozen-base option: the residual of the first nfix_base physical
      ! cells is replaced by the anchor row F_j = Y_j - Yfix_j, pinning them
      ! to the warm-start state, so that Newton solves for the WIND on top of
      ! a given base state. Set via set_base_fix before a solve.
      ! NB: it was previously assumed that the cell-1 momentum row cannot be
      ! satisfied at all, the lower BC (hard-pinned rho_bc/p_bc ghosts + the
      ! v-valve max(v,0)) leaving it non-zero at any steady interior solution.
      ! Direct measurement on the HD 209458 b hand-off state does not support
      ! that: R(2,1) is a four-order cancellation of gravity against the
      ! pressure gradient whose remainder is 6.5e-5 of the gravity term, and a
      ! 2.9 ppm change of the ghost pressure drives it to zero
      ! (docs/newton_scaling_and_base_wall.md). The option is kept as an
      ! experiment, not as a remedy for a base that cannot be solved.
      ! Trial states this solve refused because a cell of their chemistry was
      ! not a root of the network (see eval_residual). Counted here
      ! rather than in an ioniz_eq
      ! ledger, so that what a probe state did stays out of the run's own
      ! acceptance statistics. Reset at the start of each solve.
      integer :: n_trial_no_chem_root   = 0
      integer :: no_chem_root_cell_max  = 0
      real*8  :: no_chem_root_res_max   = 0.0d0
      ! Cells the LAST equilibrium sweep accepted on a state that is NOT a
      ! root of the chemical network (acceptance class 4 or 6), the same count for
      ! the state the solve currently holds, and for the best iterate kept.
      ! The trial test is the comparison of the first two (eval_residual); the
      ! state count is what the final gate is given, so that the certification
      ! refers to the state handed back and not to the last trial evaluated.
      integer :: n_no_chem_root_last    = 0
      integer :: n_no_chem_root_state   = 0
      integer :: n_no_chem_root_best    = 0

      integer :: nfix_base = 0
      real*8, allocatable :: Yfix_base(:)

      ! THE CARRIER ROWS OF THE LAST FULL RESIDUAL EVALUATION (section 139).
      ! With a transported balance among the unknowns the state has a row per
      ! cell that is not a hydrodynamic one: it is assembled by
      ! diffusive_photochemistry, in cm^-3 s^-1, on the frozen background the
      ! same evaluation's equilibrium sweep left. These hold the worst of
      ! them, its volume-weighted measure and its worst cell for the state
      ! eval_residual last saw, and again for the BEST iterate -- kept
      ! separately because the last evaluation of a solve is a REJECTED
      ! trial, not the state the solve returns. The row of EACH balance
      ! separately is what the certification reports; these are the solve's
      ! own running measure.
      real*8  :: carrier_relnorm_last = 0.0d0
      real*8  :: carrier_cellmax_last = 0.0d0
      integer :: carrier_cell_worst   = 0
      real*8  :: carrier_relnorm_best = 0.0d0
      real*8  :: carrier_cellmax_best = 0.0d0
      integer :: carrier_cell_best    = 0
      ! WHICH SCREEN REFUSED THE LAST RESIDUAL EVALUATION, AND WHERE.
      ! "no Krylov direction could be sampled" and "no admissible state on
      ! the ray" are outcomes, not diagnoses: the sample was refused by ONE
      ! of four screens, and which one, at which cell and which unknown,
      ! says whether the obstruction is a bound the probe crossed, a
      ! chemistry that did not close, or an overflowed row. Written by every
      ! evaluation that sets `admissible`, read by the exits that report a
      ! refusal. It describes the evaluation just made, so a probe leaves it
      ! behind deliberately and it is not in residual_evaluation_products.
      integer, parameter :: refuse_none                     = 0
      integer, parameter :: refuse_species_unknown_negative = 1
      integer, parameter :: refuse_element_fraction_above_1 = 2
      integer, parameter :: refuse_element_budget_breach    = 3
      integer, parameter :: refuse_sweep_nonfinite          = 4
      integer, parameter :: refuse_row_nonfinite            = 5
      integer, parameter :: refuse_no_chemical_root         = 6
      integer :: eval_refusal        = refuse_none
      integer :: eval_refusal_cell   = 0
      integer :: eval_refusal_row    = 0
      real*8  :: eval_refusal_value  = 0.0d0
      ! Species unknowns of the evaluated state that sit at or below zero,
      ! i.e. on the lower bound of their own variable. A column with many of
      ! them is a state at a corner of the admissible set, where a
      ! finite-difference probe along a direction with any negative
      ! component leaves that set whatever its length.
      integer :: n_species_at_lower_bound = 0
      ! How many refusals of this solve each screen accounts for. Reset at
      ! the top of a solve; a probe's refusal counts, because a refused
      ! probe is a Newton column or a Krylov direction the solve did not
      ! get.
      integer :: n_eval_refusal(0:6) = 0
      ! Jacobian-vector products this solve had to take on the BACKWARD side
      ! of the direction because no forward step along it landed on a
      ! describable state (jv_product). Zero over a solve is the statement
      ! that every Krylov direction was sampled forward.
      integer :: n_jv_backward_sample = 0
      ! AND THE SAME MEASURE OF THE STATE THE SOLVE HOLDS -- its current
      ! iterate, or the trial it has just adopted. That is what the
      ! acceptance gate and the best-iterate ledger have to read, and it is
      ! not the last evaluation's: between the evaluation of an iterate and
      ! the gate that judges it the solve evaluates Jacobian columns, Krylov
      ! products, points on a trust-region ray and trials it then refuses,
      ! and every one of them overwrites carrier_relnorm_last. Refreshed
      ! wherever n_no_chem_root_state is, and for the same reason.
      real*8  :: carrier_relnorm_state = 0.0d0
      real*8  :: carrier_cellmax_state = 0.0d0
      integer :: carrier_cell_state    = 0
      ! Trial and probe states refused because a cell's transported density
      ! left the element budget of its own element, went negative, or an
      ! element mass fraction left [0,1]. The budget is a CONSTRAINT, not a
      ! row: clamping it would put a kink in the residual (the defect of
      ! section 138), and a state outside it is not one the write-back can
      ! describe. It is the SUM over those screens; n_eval_refusal is the
      ! census that separates them, and it had to be written because the
      ! single count was read as an element-budget count for a solve whose
      ! 514 refusals were all the positivity of one carrier unknown
      ! (report B5g section 6.3).
      integer :: n_trial_headroom  = 0
      integer :: n_headroom_last    = 0
      ! WHAT A CANDIDATE'S BREACH IS COMPARED WITH, and it is a MAGNITUDE
      ! (item N4b): the worst relative violation of a shared element row
      ! that the iterate the rows were frozen at itself carries.
      !
      ! It is set where the rows are frozen and nowhere else
      ! (freeze_element_constraint_rows), and it comes out at zero, because
      ! what a row allows is already raised to the iterate's own demand
      ! where the iterate stands above its frozen budget
      ! (elem_row_target). The comparison is then "no candidate may stand
      ! outside what the row allows by more than the representability bound
      ! of the map", which is the continuation device of a frozen budget
      ! made quantitative rather than a threshold on the budget itself.
      real*8  :: v_headroom_iterate = 0.0d0
      ! Full residual evaluations of a solve, the unit its cost is counted
      ! in: each one is an equilibrium sweep over the column.
      integer :: n_resid_eval      = 0
      ! Equilibrium sweeps of the LAST residual evaluation, and the total
      ! over a solve: the composition is eliminated at the state being
      ! evaluated (eval_residual), so a residual evaluation is no longer one
      ! sweep and the two counts are what the cost of that elimination is
      ! read from.
      integer :: n_eq_sweeps_last  = 0
      integer :: n_eq_sweeps_tot   = 0
      ! Passes the finite-difference probes of the Newton model take. TWO
      ! requirements, and together they say it must be a FIXED number and
      ! must not be 1.
      !
      ! (a) A probe must sample ONE map: F(Y + h) and the base F it is
      !     differenced against have to be the same number of passes from
      !     the same seed, or the quotient is the gap between two maps.
      ! (b) It must be enough passes for the composition to respond. The
      !     eliminated residual is R(Y) = L(Y) + S(Y, c*(Y)); its Jacobian
      !     carries dS/dc . dc*/dY, and dc*/dY is the SUM of the Picard
      !     series, (I - M)^-1 times the one-pass response. One pass returns
      !     only the first term. The convergence test cannot supply the rest
      !     at a probe point: the probe displaces Y by the square root of
      !     machine epsilon, so the composition moves by about that much per
      !     pass and the test would stop after one pass while the series is
      !     still summing. Nothing about the probe point says how far it is
      !     from the fixed point.
      !
      ! MEASURED, which is why this is written as it is: with the probes on
      ! the iterate's own count -- 1, because at the loop top the seed IS the
      ! composition of that state and the elimination is already converged --
      ! `wasp_full_newton` had the Newton direction ascending from iteration
      ! 1, GMRES taking two iterations instead of one, the damped
      ! Gauss-Newton escape firing at nearly every iteration and ||R|| still
      ! 1.8 after 18 outer iterations, against 6.06e-06 after 26 with the
      ! lagged residual. Two builds, one with the probes at the iterate's
      ! count and one with them pinned to 1, produced iteration lines that
      ! agree digit for digit: the count WAS 1 in both.
      !
      ! So the count is taken from the elimination of a state whose seed is
      ! NOT its own composition -- the marching hand-off, and afterwards the
      ! accepted trials, each seeded from the previous iterate -- and that is
      ! n_eq_sweeps_model below.
      !
      ! WITH IT, MEASURED on the same case: 15 Picard terms, ||R|| falling
      ! 1.995 -> 4.40e-07 in 9 outer iterations, info = 0 and every row of the
      ! certification within tolerance, against 26 iterations to a state that
      ! refused its own energy row at 4.81e-05 with the lagged residual.
      integer :: n_eq_sweeps_model = 1
      ! Change of the quantities the residual reads the composition FOR,
      ! below which the fixed-state sweep iteration has reached its fixed
      ! point. EXHALE_RESID_EQ_TOL overrides it.
      ! increment_of_what_the_residual_reads is the measure.
      real*8  :: eq_sweep_reltol   = 1.0d-8
      ! WHAT THE OUTER RESIDUAL MAY CARRY OF THE INNER CLOSURE'S ERROR, and
      ! the measured factor between the two.
      !
      ! eq_sweep_reltol above is what the SWEEP is asked for; this is what
      ! the outer residual is allowed to keep of it. The two are the same
      ! number only where the boundary does not amplify a composition error,
      ! and at the base it does: eval_residual refreshes the particle count
      ! of the first cell before the boundary condition is applied, because
      ! the ghost temperature depends on the composition, so a closure error
      ! there reaches the energy row through the ghost as well as through
      ! the equation of state. `eq_sweep_reltol * merit` is therefore not a
      ! bound on the outer residual's closure error, and the ray test's
      ! floor used it as one.
      !
      ! closure_amplification_at_the_base is that factor, MEASURED once per
      ! solve on the first iterate (measure_the_closure_map_at_the_base) as
      ! the change of the base cell's scaled rows between two closures
      ! divided by the closure increment the sweep certified. The sweep is
      ! then asked for outer_residual_closure_target divided by it, so that
      ! the target IS a bound on what the outer residual carries, and the
      ! ray floor reads the target.
      real*8  :: outer_residual_closure_target = 1.0d-8
      real*8  :: closure_amplification_at_the_base = 1.0d0
      ! The closure increment the last adaptive sweep exited on, so that the
      ! error the sweep certified can be read from outside it. Zero after a
      ! fixed-count evaluation, which takes no such test.
      real*8  :: eq_sweep_increment_last = 0.0d0
      ! The measure this replaced, kept switchable for measurement
      ! (EXHALE_EQ_MEASURE=species): the largest relative change of any
      ! SPECIES FRACTION of any cell, with an absolute floor of 1e-12 under
      ! its denominator.
      !
      ! WHY IT IS NOT THE MEASURE ANY MORE. In the fully ionized wind the
      ! He I and H2+ fractions are ~1e-20, eight decades below that floor,
      ! so their denominator is the floor and the quotient reports their
      ! round-off rather than a change of the state. MEASURED on the
      ! molecular hot Uranus (report B5h section 2): the increment plateaus
      ! at 6e-11 wandering among cells 439-461 in He I and H2+ and stays
      ! there for sixty passes at a fixed count, so `eq_sweep_reltol` below
      ! about 1e-10 is unreachable by construction and the passes a tighter
      ! setting buys are spent on species that cannot move the residual.
      ! It is kept because the two measures have to be comparable on one
      ! binary, not because either reading is optional.
      logical :: eq_measure_is_species_fraction = .false.
      real*8, parameter :: eq_sweep_floor = 1.0d-12
      ! Print the measure of every pass, with the channel and the cell that
      ! set it (EXHALE_EQ_MEASURE_TRACE=1). Diagnostic only.
      logical :: eq_measure_trace = .false.
      ! WHICH OF THE THREE CHANNELS SET THE MEASURE. The composition enters
      ! assemble_residual through these and through nothing else once u is
      ! fixed, which is what makes them the measure.
      integer, parameter :: eq_channel_particle_count = 1
      integer, parameter :: eq_channel_h2_caloric     = 2
      integer, parameter :: eq_channel_radiative      = 3
      ! Cap on that iteration from the environment (EXHALE_RESID_SC_MAX);
      ! negative means the parameter file's n_selfconsistent_max. 1 is the
      ! lagged residual of section 155, kept so the ladder can be measured
      ! without a rebuild.
      integer :: n_eq_sweeps_cap   = -1

      ! WHY A TRUST-REGION ITERATION TOOK NO STEP.  One code per exit of
      ! trust_region_step, because the exits ask for different repairs and a
      ! single counter cannot be read: an inadmissible state on the ray wants
      ! the ray shortened to the admissible boundary, a model that
      ! over-predicts wants the radius rule, and a finite difference taken
      ! below its own cancellation floor wants neither -- it carries no
      ! verdict at all.  Until this ledger the solve printed "model refused
      ! by the ray test" for every one of them.
      integer, parameter :: tr_step_taken                 = 0
      integer, parameter :: tr_no_krylov_direction        = 1
      integer, parameter :: tr_cauchy_probe_inadmissible  = 2
      integer, parameter :: tr_dogleg_step_is_zero        = 3
      integer, parameter :: tr_no_admissible_state_on_ray = 4
      integer, parameter :: tr_model_promises_no_drop     = 5
      integer, parameter :: tr_ray_probe_inadmissible     = 6
      integer, parameter :: tr_ray_slopes_disagree_sign   = 7
      integer, parameter :: tr_ray_slopes_disagree_size   = 8
      integer, parameter :: tr_trial_no_longer_admissible = 9
      integer, parameter :: tr_ratio_below_acceptance     = 10
      integer, parameter :: tr_ray_below_the_noise_floor  = 11
      integer, parameter :: tr_projected_step_has_no_model = 12
      integer, parameter :: tr_no_direction_reduces_model  = 13
      integer, parameter :: n_tr_exit_reason              = 13

      ! WHAT THE RAY TEST COMPARES AND WHAT IT NEEDS BEFORE IT MAY SPEAK
      ! (ray_slope_verdict).
      !
      ! ray_factor: how far the model's directional derivative and a finite
      ! difference of the true merit along the same ray may differ in size
      ! before the model is refused. A factor 4 is loose enough that
      ! ordinary finite-difference noise does not reject a good model.
      real*8,  parameter :: ray_factor = 4.0d0
      ! ray_floor_factor: how far the difference f2p - f2m must stand above
      ! the MEASURED reproducibility of the merit before the quotient built
      ! from it carries a verdict. The merit of a state is a function of Y
      ! only to the accuracy with which the composition is eliminated, and
      ! two evaluations of it at one state differ by that much; a difference
      ! of a few times that is a difference of noise. Ten, because a single
      ! pair of evaluations estimates a spread whose own scatter is of its
      ! own size, and one decade of margin is what separates "above the
      ! noise" from "of the order of the noise".
      real*8,  parameter :: ray_floor_factor = 1.0d1
      ! ratio_sound_band: how close the reduction ratio has to be to unity
      ! for the ray test to lose its power of refusal. The ratio is the
      ! actual reduction over the predicted one AT THE FULL STEP, with no
      ! finite difference in it; inside this band the quadratic model
      ! reproduces the true reduction of the step to a percent, and a slope
      ! quotient taken over a thousandth of the same step cannot contradict
      ! evidence taken over the whole of it. The over-prediction guard is
      ! the ratio itself (tr_ratio_below_acceptance) and it is untouched.
      real*8,  parameter :: ratio_sound_band = 1.0d-2
      ! eta_accept: the acceptance threshold on the reduction ratio, the
      ! standard "accept anything that reduces the merit by a tenth of what
      ! the model promised". It is at module scope because the ray test
      ! reads it: a trial the acceptance test refuses on its own is named by
      ! the ratio and not by a slope quotient (ray_slope_verdict).
      real*8,  parameter :: eta_accept = 1.0d-1
      ! cauchy_leg_min_fraction: below this fraction of the radius the
      ! approximate-gradient leg is not part of the step and is dropped
      ! from it. The leg's direction is the transpose of the BANDED operator
      ! applied to r0, and its length tau ||g|| = (r0.Ag)||g||/||Ag||^2
      ! collapses where the band's row scales make ||Ag|| stand decades
      ! above r0 . A g: MEASURED on the atomic element reload, ||sU||/delta
      ! was 3.3e-4 at entry and 1.7e-6 at outer iteration 60 (N7 report,
      ! proposal P3). A leg that short moves the step by less than the
      ! radius resolves, while A sU is not small at all and dominates the
      ! model slope r0 . A s -- so the slope described a leg that was not in
      ! the step, and the ray test rightly refused it (N7b). Dropping the
      ! leg makes the model slope the slope of the step actually taken; the
      ! dogleg is then the Krylov leg cut to the ball, which does predict a
      ! reduction whenever the cycle reduced the model residual.
      ! A thousandth: three decades below the radius is below the accuracy
      ! of the radius's own control, which moves it by factors of 4 and 2.
      real*8,  parameter :: cauchy_leg_min_fraction = 1.0d-3
      ! cauchy_leg_image_share_max, cauchy_leg_step_share_max: the leg is
      ! dropped as well where it is nearly none of the step and nearly all
      ! of the model's image of it. The step is s = aU sU + aN sN and, by
      ! linearity, A s = aU A sU + aN A sN, so each leg owns a share of the
      ! length and a share of the image and the two are independent
      ! quantities. Only the image share enters the model slope r0 . A s,
      ! so a leg that supplies the image while supplying no length makes
      ! that slope a statement about a direction the step does not contain,
      ! and the predicted decrease 0.5 (r0.Ag)^2/||Ag||^2 it then carries
      ! has no radius in it at all (N7b; N22 report, noticed item 1).
      ! MEASURED at outer iteration 145 of the atomic element reload:
      ! ||sU||/delta was 6.856e-3, which the length test above keeps
      ! because it stands above cauchy_leg_min_fraction, and the step being
      ! on the ball that is also its share of the step; while ||A s||
      ! equalled tau_c ||A g|| to four digits, an image share of 1.000,
      ! because ||A g|| = 7.923e14 stands four decades above
      ! r0 . A g = 8.262e13.
      ! Half the image is the threshold because a leg that owns more of the
      ! image than the whole of the rest of the step does is what sets the
      ! slope; a tenth of the step is the threshold below which the leg is
      ! not the step, one decade below the factors of 4 and 2 the radius
      ! control moves the region by.
      real*8,  parameter :: cauchy_leg_image_share_max = 5.0d-1
      real*8,  parameter :: cauchy_leg_step_share_max  = 1.0d-1
      ! cauchy_leg_admitted_by_image_alone: which of the two rules above
      ! decides. The length test asks whether the leg is small beside the
      ! RADIUS; the question the model asks is whether the leg is small
      ! beside the STEP while owning the image, and the two agree only
      ! where the step stands on the ball. With this true the length test
      ! is not taken and the image test alone decides, so a leg that is
      ! neither the step nor the image is KEPT however short it is, which
      ! costs nothing (both images are formed already) and leaves the model
      ! slope a statement about the step, the leg's image share being
      ! small.
      ! OFF BY DEFAULT BECAUSE IT IS MEASURED TO BE WORSE, not because it
      ! is untried (N24). The two rules select different legs from the
      ! seventh outer iteration of the atomic element reload onward, and
      ! there the image rule alone hands back ||R|| 2.296e-2 against the
      ! length rule's 1.356e-2 (71 legs dropped on their image, none on
      ! their length); on the carrier reload it hands back 3.244e-1
      ! against 8.993e-2 and runs to its iteration cap where the length
      ! rule stopped on the stagnation detector.
      logical :: cauchy_leg_admitted_by_image_alone = .false.

      ! WHAT ONE CYCLE OF THE LINEAR SOLVE RETURNED. One name per outcome,
      ! because the repair each of them asks for is a different one and a
      ! single "truncated" flag cannot separate them (R3, F5): reaching the
      ! subspace budget with too little reduction, failing to sample a
      ! direction at all, a preconditioner that is not a map, a reduced
      ! problem that has lost rank, and arithmetic that left the reals are
      ! five different statements about the step that comes back.
      !
      ! A VANISHING ARNOLDI REMAINDER IS NEITHER SUCCESS NOR FAILURE BY
      ! ITSELF. The Krylov space is then invariant, so no further direction
      ! exists; whether the step in it solves the system is decided by the
      ! rank of the reduced problem and confirmed by the TRUE residual, not
      ! by the norm of the remainder.
      integer, parameter :: gm_tolerance_reached        = 0
      integer, parameter :: gm_subspace_exhausted       = 1
      ! No direction at all: the Jacobian action along the first
      ! preconditioned direction has no admissible sample, which on a
      ! bound-constrained unknown space is a statement about the feasible
      ! set and not about the step length (jv_product,
      ! largest_step_inside_the_species_box). The step is zero.
      integer, parameter :: gm_no_direction_sampled     = 2
      ! The action was refused, or was not a finite number, after at least
      ! one direction had been built: the step is the minimizer over the
      ! subspace that WAS constructed.
      integer, parameter :: gm_jacobian_action_unusable = 3
      integer, parameter :: gm_preconditioner_failed    = 4
      integer, parameter :: gm_reduced_system_singular  = 5
      integer, parameter :: gm_nonfinite                = 6
      ! The iterate left the trust ball and the point on its boundary is
      ! what comes back: a cycle stopped by the region, not by its subspace
      ! or by its arithmetic (krylov_truncated_on_the_trust_ball).
      integer, parameter :: gm_stopped_on_the_trust_ball = 7
      ! The cycle was stopped because the TRUE residual of the step it would
      ! return has risen at two consecutive checks
      ! (gm_step_by_its_true_residual): the reduced least-squares residual
      ! keeps falling, the residual the step reaches against the operator
      ! does not, and the step handed back is the best true residual the
      ! cycle saw. Only this arm can name it.
      integer, parameter :: gm_true_residual_turned      = 8
      integer, parameter :: n_gm_outcome_kind           = 8
      ! The reduced problem's rank, and an Arnoldi breakdown, are decided
      ! RELATIVE to the size of the numbers they are read from: a triangular
      ! diagonal entry this far below the largest one carries no information
      ! about its unknown, and a remainder this far below the column it
      ! belongs to adds no direction. 1e2 times the double-precision
      ! epsilon: the Arnoldi recursion loses orthogonality at a rate set by
      ! the conditioning, so a threshold at epsilon itself would call a
      ! rounding-level entry usable.
      real*8,  parameter :: gm_rank_relative = 1.0d2*epsilon(1.0d0)
      ! WHICH DIRECTION THE TRUST REGION'S FIRST RADIUS WAS MEASURED ON,
      ! and whether there was one. A radius is a length, so it can only be
      ! set from a direction the step may actually take: an arbitrary
      ! positive number does not create a direction that is missing (R1).
      integer, parameter :: tr_radius_from_the_krylov_step = 1
      integer, parameter :: tr_radius_from_the_gradient    = 2
      integer, parameter :: tr_radius_no_feasible_descent  = 3
      ! The first radius is this fraction of the direction it is measured
      ! on: the ray scans of section 142 put the whole of the merit's
      ! descent inside t = 1e-2 of the Newton step and the rise dominant by
      ! t = 0.1, so the region starts at a hundredth of that length and is
      ! then controlled by the reduction ratio.
      real*8,  parameter :: tr_first_radius_fraction = 1.0d-2
      ! The ceiling starts a hundredfold above the first radius and then
      ! follows it (see the radius update in solve_steady_jfnk).
      real*8,  parameter :: tr_radius_ceiling_factor = 1.0d2
      ! A POSITIVE FLOOR, in the SCALED coordinates the step lives in. There
      ! D is already divided out, so the floor is a dimensionless length and
      ! the only thing it may depend on is the number of unknowns: a
      ! Euclidean norm over neq components of a state that is O(1) in each
      ! is O(sqrt(neq)), so a floor of 1e-6 of sqrt(neq) is six decades
      ! below the scale of the state and cannot mask a genuinely tiny valid
      ! direction, while it does stop a radius of exactly zero from being an
      ! absorbing state. MEASURED against the radius the two species-row
      ! fixtures actually take (item N1 report), which is far above it.
      real*8,  parameter :: tr_delta_min_dimensionless = 1.0d-6

      ! --- THE STATE A STATIONARY SOLVE CARRIES BESIDE ITS ITERATE, AND
      !     WHAT A RESTART OF IT RESETS ---
      !
      ! A trust-region solve carries five things that are not the iterate:
      ! the radius and its ceiling, the pseudo-transient shift the banded
      ! preconditioner is built with, the tolerance the composition
      ! elimination is asked for (measured once from the closure map at the
      ! base), the Grippo window and monotone flag the acceptance compares
      ! against, and the ledger of the best iterate seen. Entering a solve
      ! afresh takes all five from the state it is handed; continuing an
      ! iteration carries them from the iteration before.
      !
      ! The band itself is NOT among them: it is rebuilt from the iterate
      ! and refactored by dgbtrf at the top of every outer iteration, so a
      ! re-entry changes the preconditioner only through the shift.
      ! Likewise the Krylov cycle starts from x = 0 at every outer
      ! iteration, so a re-entry resets nothing there; a second cycle of
      ! the same subspace size is a change to the linear solve and is
      ! counted separately (pgmres_with_restarts).
      integer, parameter :: tr_reset_radius_and_ceiling = 1
      integer, parameter :: tr_reset_pseudo_transient   = 2
      integer, parameter :: tr_reset_closure_map        = 3
      integer, parameter :: tr_reset_acceptance_memory  = 4
      integer, parameter :: tr_reset_back_to_best       = 5
      integer, parameter :: n_tr_reset_kind             = 5
      ! WHICH OF THE FIVE A RESTART OF THIS SOLVE TAKES
      ! (read_trust_region_restart_controls). All five is what a re-entry
      ! does; a single one is how the arms of item N21 separate them.
      logical :: tr_reset_wanted(n_tr_reset_kind) = .true.
      ! The outer iteration at which one restart is forced whatever the
      ! trigger says, 0 for none (EXHALE_TR_RESTART_AT). It exists so that
      ! the restart can be measured at one named iterate, on the state a
      ! capped solve hands back.
      integer :: tr_restart_forced_at = 0
      ! THE TRIGGER, and the statistic it reads: the number of consecutive
      ! outer iterations in which the best iterate did not improve. That is
      ! the stagnation detector's own statistic (n_since_best against
      ! n_stall_best in solve_steady_jfnk), read at a fraction of the count
      ! that stops the solve, so a restart is attempted while the solve
      ! still has iterations left rather than at the point it gives up.
      ! Zero disarms it.
      integer :: tr_restart_without_improvement = 0
      ! THE DEFAULT OF THAT TRIGGER. It is 0 -- disarmed -- until the
      ! restart is measured to help on both species-row fixtures; the
      ! measurement of item N21 is what may change it, and the arms are
      ! reached through EXHALE_TR_RESTART_STALL.
      integer, parameter :: tr_restart_stall_default = 0
      ! How many restarts this solve has taken, and of which kind.
      integer :: n_tr_restart = 0
      ! HOW MANY CYCLES OF THE SUBSPACE THE KRYLOV SOLVE MAY SPEND, the
      ! subspace size itself unchanged (EXHALE_GM_CYCLES, pgmres_with_
      ! restarts). One is a single cycle, which is what GMRES(m) without
      ! restarts is.
      integer :: gm_restart_cycles = 1
      ! What EXHALE_GM_CYCLES asked for, and the first outer iteration at
      ! which it is spent (EXHALE_GM_CYCLES_FROM, 1 by default). The two
      ! exist so that a restarted cycle can be armed part way into a solve,
      ! on the same prefix of iterations a restart of the trust-region
      ! state is armed on, and the two arms then differ in one thing.
      integer :: gm_restart_cycles_wanted = 1
      integer :: gm_restart_cycles_from   = 1

      ! WHETHER THE ARNOLDI BASIS IS STILL ORTHOGONAL, measured on request
      ! (EXHALE_GM_ORTHO=1) as max |V_i . V_j| over i < j of the basis each
      ! cycle builds, and kept as the worst of the solve. Modified
      ! Gram-Schmidt loses orthogonality at a rate set by the conditioning
      ! of the operator; whether reorthogonalization is worth its products
      ! is decided by this number and not assumed. Off by default: the
      ! measurement costs the same inner products a second orthogonalization
      ! pass would.
      logical :: gm_measure_orthogonality   = .false.
      real*8  :: gm_orthogonality_loss_max  = 0.0d0
      ! The same quantity for the LAST cycle alone, so that the basis of one
      ! named outer iteration can be read rather than the worst of a solve.
      real*8  :: gm_orthogonality_loss_cycle = 0.0d0
      ! WHETHER THE ARNOLDI RECURSION ORTHOGONALIZES TWICE
      ! (EXHALE_GM_REORTHO=1). Modified Gram-Schmidt once leaves
      ! max |V_i . V_j| growing with the conditioning of the operator; a
      ! second pass over the same basis removes the components along the
      ! earlier basis vectors that the first pass leaves, at the cost of the
      ! inner products alone and of no further product of the operator. The
      ! condition for taking it is a MEASURED loss above
      ! gm_reorthogonalization_loss_level; below that level the basis is
      ! orthogonal to the precision the model is read in and the second pass
      ! buys nothing.
      logical :: gm_reorthogonalize = .false.
      ! The level the measured loss has to exceed for the second pass to be
      ! worth its inner products: the model residual is read to a relative
      ! 1e-1 and the predicted decrease to a few digits, so a basis whose
      ! worst inner product is below 1e-8 cannot be what makes either of
      ! them wrong (the condition item N3 stated, where the carrier fixture
      ! measured 6.9e-13).
      real*8,  parameter :: gm_reorthogonalization_loss_level = 1.0d-8
      ! WHETHER THE IMAGE THE ARNOLDI RELATION RETURNS IS CHECKED AGAINST
      ! THE OPERATOR, basis vector by basis vector (EXHALE_GM_IMAGE_CHECK=1,
      ! off by default because each check costs one product of the
      ! operator). The trust region's predicted decrease is read from
      ! A_z M^-1 V_k = V_{k+1} Hbar_k; where that relation has drifted, the
      ! model is not the operator's and the drift is a property of one
      ! column of the recursion, which a single check of the assembled step
      ! cannot separate. The check also samples one direction twice, so that
      ! a matrix-free action which is not reproducible is told apart from a
      ! basis that has lost orthogonality.
      logical :: gm_verify_the_arnoldi_image = .false.
      ! WHETHER THE KRYLOV CYCLE PRINTS ITS RESIDUAL HISTORY
      ! (gm_residual_history_on, EXHALE_GM_HISTORY=1). The reduced
      ! least-squares residual after every product, and every
      ! gm_history_stride products the TRUE residual ||b - A_z x_j||/||b||
      ! of the iterate that subspace gives, which costs one product of the
      ! operator and one triangular solve. The two separate a cycle whose
      ! reduced problem is converging while the operator's residual is not
      ! (a basis that has lost orthogonality, or a rank-deficient reduced
      ! problem) from a cycle that is simply stagnating.
      logical :: gm_residual_history_on   = .false.
      ! Armed only at the outer iteration the element diagnostic speaks at,
      ! and disarmed inside the subspace-size scan so that the scan's own
      ! cycles do not print two histories at once.
      logical :: gm_residual_history_here = .false.
      integer, parameter :: gm_history_stride = 20
      ! WHETHER THE CYCLE IS REPEATED AT SEVERAL SUBSPACE SIZES AT THE SAME
      ! ITERATE (krylov_size_scan_on, EXHALE_KRYLOV_SIZE_SCAN=1). A cycle
      ! that reaches its tolerance once the subspace is large enough is a
      ! conditioning that GMRES beats by size; one that stagnates at every
      ! size has a part of the operator the preconditioner does not touch.
      ! The scan does not adopt any of its steps: the solve's own cycle has
      ! already run and the trajectory is unchanged.
      logical :: krylov_size_scan_on = .false.
      ! THE SUBSPACE SIZES THE SCAN RUNS, and whether the last two are
      ! repeated with the basis orthogonalized twice.
      integer, parameter :: n_krylov_scan_size = 4
      integer, parameter, dimension(n_krylov_scan_size) ::               &
                          krylov_scan_size = (/ 40, 80, 160, 320 /)
      ! WHETHER THE RITZ VALUES OF THE PRECONDITIONED OPERATOR ARE MEASURED
      ! (precond_spectrum_on, EXHALE_PRECOND_SPECTRUM=1). An Arnoldi
      ! recursion of n_ritz_step products from a deterministic start, whose
      ! Hessenberg is diagonalized densely: the Ritz values approximate the
      ! spectrum of A_z M^-1, and a cluster of them near zero is the part of
      ! the operator the band leaves behind. Run on the whole operator and
      ! on its compressions onto the species rows and onto the hydrodynamic
      ! rows, so that the cluster can be attributed to a row class.
      logical :: precond_spectrum_on = .false.
      integer, parameter :: n_ritz_step = 200
      ! The Ritz vectors of the three smallest-magnitude Ritz values of the
      ! whole operator, kept so that the band-difference measurement below
      ! can be taken on the directions the spectrum names. Empty until the
      ! spectrum hook has run at this iterate.
      real*8, allocatable :: ritz_vector_of_the_smallest(:,:)
      real*8, allocatable :: ritz_value_of_the_smallest(:)
      ! WHETHER THE BAND IS COMPARED WITH THE FULL JACOBIAN ACTION
      ! (band_difference_on, EXHALE_BAND_DIFFERENCE=1). On the directions
      ! the spectrum named and on the first two Arnoldi directions of the
      ! solve's own right-hand side: the difference (A - A_band) v split by
      ! row class and by cell, and, on the columns those directions live on,
      ! the split of the true column into the part inside the band and the
      ! part outside it, with the radiation frozen and with it live.
      logical :: band_difference_on = .false.
      ! WHETHER THE ROW THAT BINDS IS READ TERM BY TERM AND ITS JACOBIAN
      ! ENTRIES ATTRIBUTED TO THOSE TERMS (front_row_on,
      ! EXHALE_FRONT_ROW=1, off by default).
      !
      ! The species row that carries the largest scaled residual at the
      ! selected outer iteration is a control-volume balance of a
      ! transported species: an advective divergence of the face mass flux,
      ! a diffusive or settling divergence, and, for a carrier, a chemical
      ! source and sink. Its diagonal in the banded model is its response
      ! to its OWN unknown, and at a front that diagonal was measured three
      ! decades below the row's response to the hydrodynamic unknowns of the
      ! same stencil. Which TERM of the row carries each of those responses
      ! is not readable off the assembled row, so it is measured: every
      ! term of the row is re-formed at the state and at a state displaced
      ! along one column, and the difference quotient of each term is
      ! reported beside the band's entry for that column, which their sum
      ! reproduces.
      !
      ! The measurement costs one residual evaluation for each column
      ! probed, four unknowns in each of five cells, and adopts nothing.
      logical :: front_row_on = .false.
      ! The half width in cells of the column stencil the attribution
      ! probes, around the binding cell: the band reaches one cell either
      ! side of a row, and the two neighbor rows reported beside the
      ! binding one reach one cell further.
      integer, parameter :: n_front_row_halfwidth = 2
      ! EVERY TERM OF ONE TRANSPORTED SPECIES ROW, in the row's own
      ! physical units, together with the two face quantities its advective
      ! term is built from.  A transported species row is a control-volume
      ! balance,
      !
      !     div(J_diffusive) + div(F_rho Y^face) - S_chemical = 0 ,
      !
      ! so the three terms below sum to the row and the two face
      ! contributions sum to the advective one.  to_code is the factor that
      ! writes the row per code time in code density units, the one the
      ! residual vector carries.
      type species_row_term_set
         real*8  :: row          = 0.0d0
         real*8  :: diffusive    = 0.0d0
         real*8  :: adv_total    = 0.0d0
         real*8  :: adv_left     = 0.0d0
         real*8  :: adv_right    = 0.0d0
         ! The same advective divergence with both faces carrying the donor
         ! cell's own average instead of the reconstructed face value: the
         ! first-order upwind term the limited reconstruction sits on.
         real*8  :: adv_donor    = 0.0d0
         real*8  :: chemical     = 0.0d0
         real*8  :: row_scale    = 0.0d0
         real*8  :: frho_left    = 0.0d0
         real*8  :: frho_right   = 0.0d0
         real*8  :: yface_left   = 0.0d0
         real*8  :: yface_right  = 0.0d0
         real*8  :: ydonor_left  = 0.0d0
         real*8  :: ydonor_right = 0.0d0
         ! The unknown of the row in the code's units, and the physical
         ! density it names.
         real*8  :: unknown      = 0.0d0
         real*8  :: density      = 0.0d0
         real*8  :: to_code      = 1.0d0
         logical :: ok           = .false.
      end type species_row_term_set
      ! WHETHER THE CYCLE STOPS ON THE TRUST BALL (Steihaug-Toint
      ! truncation, EXHALE_KRYLOV_ON_THE_BALL=1, off by default).
      !
      ! The cycle asked for a relative residual 1e-1 returns the minimizer
      ! over its whole subspace, and the trust region then uses the fraction
      ! of it that fits inside the radius: MEASURED on the atomic element
      ! reload, ||sN|| 8.647e+02 against a radius 5.918e-02, so 7e-05 of a
      ! direction that spent 40 products. A truncated cycle stops where the
      ! iterate leaves the ball and returns the point on its boundary, so
      ! the products are spent on the step the region can take.
      !
      ! The GMRES iterate is not monotone in norm, so "the boundary" is
      ! taken as the FIRST iterate whose norm exceeds the radius, backtracked
      ! along the last increment to the boundary: x = x_{k-1} + t (x_k -
      ! x_{k-1}) with ||x|| = delta. Both endpoints carry their own image
      ! from the Arnoldi relation and the image of the returned point is the
      ! same combination of the two, so the model of the step is the model of
      ! the point returned and costs no further product.
      logical :: krylov_truncated_on_the_trust_ball = .false.

      ! ROW EQUILIBRATION OF THE LINEAR MODEL ONLY
      ! (model_row_equilibration_on, EXHALE_MODEL_ROW_EQUIL=1, off by
      ! default).
      !
      ! WHAT IS FREE AND WHAT IS NOT. The merit, the gate, the trust
      ! region's predicted decrease and the certification all read the
      ! residual on the certification row scales (decision 20 a,
      ! merit_row_scale_from_certification), and that fixes r0 = F/Drow.
      ! It does not fix the row scaling of the LINEAR MODEL: solving
      ! (E A) s = E r0 with E diagonal and positive is the same linear
      ! system, and its exact solution is the same step. What E changes is
      ! the norm the truncated Krylov cycle minimizes, and with it the
      ! meaning of "relative residual 1e-1" and the subspace the cycle
      ! selects. THIS ARM CHANGES THE LINEAR MODEL'S INNER PRODUCT AND
      ! NOTHING THE SOLVE IS JUDGED BY: the step comes back in the same
      ! coordinates, its image A s is mapped back to the certification
      ! scales before it leaves the cycle, and the relative residual the
      ! caller reads is the certification one.
      !
      ! E IS UNIT INFINITY NORM PER ROW OF THE BANDED MODEL, not unit
      ! diagonal. The row that binds on the carrier fixture carries a
      ! diagonal 3.3e3 below its coupling to the hydrodynamic unknowns of
      ! its own stencil (section N35), so a scaling to unit diagonal would
      ! multiply that row by three decades against the rest of the system
      ! and equilibrate nothing; and a row scaling of either kind leaves
      ! the diagonal-to-off-diagonal ratio WITHIN a row exactly where it
      ! was, which is the reason neither cures that row.
      !
      ! The factor is a power of two, so the rescaling of the band and of
      ! the right-hand side is exact in binary floating point.
      logical :: model_row_equilibration_on = .false.
      ! E, and whether the array describes the band currently factored.
      real*8, allocatable :: model_row_equilibration(:)
      logical :: model_rows_equilibrated = .false.
      ! The extremes of E at the last outer iteration, for the diagnostic.
      real*8  :: model_row_equil_min = 1.0d0, model_row_equil_max = 1.0d0

      ! THE STEP THE CYCLE RETURNS CHOSEN BY ITS TRUE RESIDUAL
      ! (gm_step_by_its_true_residual, EXHALE_GM_TRUE_RESIDUAL=1, off by
      ! default).
      !
      ! The reduced least-squares residual |gg(j+1)| is the residual of the
      ! step only while the Arnoldi relation is a statement about the
      ! operator. The action here is the secant of a NONLINEAR residual
      ! over a fixed probe arc, and MEASURED at the binding iterate of both
      ! fixtures (section N35) the reduced residual falls monotonically
      ! with the subspace size while the residual the returned step reaches
      ! against the operator stops falling at 60 to 80 products and rises.
      ! With this arm the cycle forms its candidate step every
      ! gm_true_residual_stride products from gm_true_residual_first on,
      ! and at every point where the cycle would otherwise end -- its last
      ! product, an Arnoldi breakdown, the reduced problem reaching the
      ! tolerance -- measures ||b - A x_k||/||b|| against the operator at
      ! one product per check, keeps the best it has seen, and stops when
      ! that residual has risen at gm_true_residual_rises consecutive
      ! checks or has reached the tolerance asked for IN THE TRUE NORM.
      ! The reduced residual then ends no cycle by itself. What comes back
      ! is that best step, its own image from the product the check took,
      ! and the true relative residual as gm_resid_rel, so the step control
      ! reads a residual and a model image that mean what they say.
      logical :: gm_step_by_its_true_residual = .false.
      integer, parameter :: gm_true_residual_first  = 20
      integer, parameter :: gm_true_residual_stride = 10
      integer, parameter :: gm_true_residual_rises  = 2

      ! THE OPERATOR THE LINEAR SOLVE IS MEASURED ON.
      !
      ! src/tests/krylov_and_dogleg drives the production Krylov cycle with
      ! a dense operator and an identity preconditioner that the test states
      ! itself, so that each named outcome of the cycle is read from the
      ! routine the solver runs and not from a copy of its text. With no
      ! operator installed -- every run of the code -- the cycle samples the
      ! true Jacobian action (jv_product) and the banded factorization
      ! (dgbtrs), and the two dispatch points below are that call unchanged.
      ! The same convention as carrier_headroom_set_for_test.
      !
      ! An operator whose action is not a finite number is stated by putting
      ! a NaN in the dense matrix itself, which is the map the cycle then
      ! samples; there is no separate mode for it.
      integer, parameter :: lin_test_absent = 0
      integer, parameter :: lin_test_dense  = 1
      ! The dense operator with a preconditioner that reports a LAPACK
      ! failure: the vector is returned untouched and the status is nonzero.
      integer, parameter :: lin_test_dense_preconditioner_fails = 2
      ! The dense operator whose action STOPS BEING THE SAME LINEAR MAP
      ! after lin_test_drift_after products: from there on the action is
      ! (A + c (k - lin_test_drift_after) I) v with k the number of
      ! products taken. It states, in the small, what the matrix-free
      ! action does in the large -- the secant of a nonlinear residual is
      ! not one operator over a whole cycle (section N35) -- so that a
      ! cycle whose reduced residual keeps falling while the residual its
      ! step reaches against the current action rises can be exercised
      ! without a reload.
      integer, parameter :: lin_test_dense_action_drifts = 3
      integer, parameter :: lin_test_drift_after = 40
      real*8,  parameter :: lin_test_drift_per_product = 1.0d0
      integer :: lin_test_kind = lin_test_absent
      integer :: lin_test_n_product = 0
      real*8, allocatable :: lin_test_A(:,:)

      ! WHETHER A TRIAL IS WRITTEN ONTO THE FACES OF THE SPECIES BOX before
      ! it is evaluated. A species density and a species mass fraction are
      ! bounded unknowns and the trust region's radius cannot express that;
      ! why the step is projected rather than shortened, and what shortening
      ! it instead was measured to cost, is at
      ! species_unknowns_outside_their_bounds.
      ! EXHALE_SPECIES_BOUND_PROJECT=0 restores the pure step shortening for
      ! measurement. Inert on the three-unknown route, which has no species
      ! unknown to bound.
      logical :: project_species_trial = .true.
      ! Whether an unknown at a bound is held whichever way the model wants
      ! to move it (EXHALE_SPECIES_BOUND_HOLD=all), instead of only when the
      ! model pushes it out. A measurement switch; see
      ! fix_active_species_bounds.
      logical :: hold_every_bound = .false.
      ! How many trials of a solve had to be projected, and how many of them
      ! then had no describable sample of the model. Reported with the
      ! carrier row: a solve spending its iterations on the faces of the box
      ! is a different thing from one descending its interior.
      integer :: n_projected_trials = 0
      ! The smallest species unknown, as the DENSITY or the mass fraction it
      ! names, that any state this solve ADOPTED carried. A species density
      ! and a species mass fraction are non-negative, so this may not fall
      ! below zero: the screen of eval_residual refuses a trial that leaves
      ! that set and the trust region writes its trial onto the face of it,
      ! and this is the number that says the two together held. Reported at
      ! the end of a solve that carried a species row.
      real*8  :: species_unknown_min_adopted = 0.0d0
      ! WHICH UNKNOWNS THE STEP OF THIS OUTER ITERATION MAY NOT MOVE: the
      ! species unknowns already ON a bound of their own box whose model
      ! gradient points out of it. Fixed once per step, so that the operator
      ! the Krylov cycle samples is LINEAR: a projection that depended on
      ! the sign of the direction it is applied to would not be.
      ! .false. everywhere outside a trust-region step, where nothing is
      ! held.
      logical, allocatable :: species_bound_is_active(:)
      ! The most unknowns any step of a solve held on their bounds.
      integer :: n_active_bounds_max = 0
      ! How many steps of a solve were the Cauchy point alone, the Krylov leg
      ! being unavailable. A trust region always has that step; before it was
      ! taken those iterations were lost.
      integer :: n_cauchy_only_steps = 0
      ! How many dogleg trials of a solve left the approximate-gradient leg
      ! out of the step because it was shorter than cauchy_leg_min_fraction
      ! of the radius. Counted so that the leg's cost -- one Jacobian
      ! product per outer iteration -- can be read against what it buys.
      integer :: n_cauchy_leg_short = 0
      ! How many dogleg trials of a solve left the approximate-gradient leg
      ! out of the step because it carried nearly the whole model image of
      ! the step while being nearly none of its length
      ! (cauchy_leg_image_share_max, cauchy_leg_step_share_max).
      integer :: n_cauchy_leg_all_image = 0
      ! THE COLUMN SCALES cell_state_scales LAST FORMED, kept aside from the
      ! column scaling the linear algebra works in. The two are the same
      ! array only while the equilibration below is off: the state scale is
      ! the characteristic magnitude of an unknown and it is what the
      ! finite-difference probe step of the banded Jacobian is taken as
      ! (build_banded_jac_full), so a length in unknown space -- the width
      ! of the epsilon-active set of the species box, for instance -- is
      ! measured on THIS array and not on the equilibrated one.
      real*8, allocatable :: state_column_scale(:)
      ! COLUMN EQUILIBRATION OF THE SCALED OPERATOR (EXHALE_COLUMN_EQUIL=0
      ! turns it off, and then the column scaling is the state scaling
      ! exactly). Species-row route only; a three-unknown solve never
      ! reaches it.
      !
      ! The row scaling of the species-row system is fixed by the gate: the
      ! merit is read on the certification scales (cell_row_scales), so
      ! r0 = F/Drow is not free. The one freedom left in the scaled operator
      ! A = Drow^-1 J D is the column scaling D, and equilibrating it is
      ! what keeps a band whose rows span eleven decades factorizable.
      logical :: column_equilibration_on = .false.
      ! TWO-SIDED EQUILIBRATION OF THE BAND THAT IS FACTORIZED
      ! (EXHALE_PRECON_EQUIL=0 turns it off). Species-row route only.
      !
      ! It is the PRECONDITIONER that is rescaled and not the model: ab, the
      ! merit, the trust-region coordinates and the radius are all left
      ! where the row scaling of the gate and the state scales of the
      ! unknowns put them, and only the matrix handed to dgbtrf is
      ! equilibrated. A preconditioner is an approximate inverse and any
      ! invertible rescaling of it is another one, so this is free of the
      ! model; what it buys is a factorization that is not numerically
      ! singular.
      logical :: preconditioner_equilibration_on = .true.
      ! Whether the factors below describe the band currently factored, and
      ! the two diagonal rescalings, as powers of two. M^-1 v is then
      ! C * dgbtrs(R v) (banded_preconditioner_solve).
      logical :: preconditioner_equilibrated = .false.
      real*8, allocatable :: precon_row_scale(:), precon_col_scale(:)
      ! The smallest column scale an element unknown is given, as a fraction
      ! of the reservoir its column is anchored to. The arithmetic that
      ! fixes it is at the use (cell_state_scales).
      real*8, parameter :: element_column_scale_floor = 1.0d-4
      ! The extremes of the equilibration factor of the last outer
      ! iteration, and the column norms of the scaled operator before and
      ! after it, for the diagnostic.
      real*8  :: colequil_min = 1.0d0, colequil_max = 1.0d0
      real*8  :: colnorm_before_min = 0.0d0, colnorm_before_max = 0.0d0
      real*8  :: colnorm_after_min  = 0.0d0, colnorm_after_max  = 0.0d0
      ! The largest entry of each column of the scaled operator BEFORE the
      ! equilibration, so that the census of unsampled columns keeps its
      ! meaning: after the equilibration every column is of size one by
      ! construction and the census would read zero whatever the operator.
      real*8, allocatable :: scaled_column_norm_before(:)
      ! A CARRIER UNKNOWN IS ITS LOGARITHM, ln n. That is the unknown space
      ! of the coupled route (Coupled carrier solve: True);
      ! EXHALE_CARRIER_LOG_UNKNOWN=0 restores the density unknown, which is
      ! there to be measured against.
      !
      ! WHY. A carrier density is positive and the stationary system has no
      ! way of knowing it. With n itself the unknown the Newton proposes
      ! negative densities, so every trial, probe and Krylov direction has to
      ! be refused or projected at a face of a box, and an unknown that has
      ! REACHED zero cannot be sampled at all: a two-sided finite difference
      ! along a direction with an outward component leaves the admissible set
      ! on both sides. In ln n positivity is a property of the
      ! parametrization -- exp of any real is positive -- so the face is not
      ! there to be negotiated with. MEASURED on `mol_carrier` (item B5j),
      ! the same case and the same tolerance: ||R|| 1.415 with the density
      ! unknown and its bound-aware region against 0.9972 with ln n, the
      ! merit 7.51 against 0.924, 28 of 64 steps the Cauchy point alone
      ! against none, and the smallest carrier density any adopted iterate
      ! carried exactly 0 against 2.2e-11.
      !
      ! WHAT IT COSTS, in three places.
      !   * THE JACOBIAN COLUMN CARRIES THE FACTOR n: dR/d ln n = n dR/dn.
      !     The column scale of the unknown is 1 (a logarithm is its own
      !     scale, cell_state_scales), so the banded model, which
      !     differentiates the residual column by column, and the
      !     matrix-free action jv_product, which differentiates it along one
      !     direction, both step ln n by a multiple of sqrt(eps_mach) and so
      !     both probe a RELATIVE change of the density. Where the carrier
      !     vanishes that column vanishes with it, and whether the model then
      !     has an empty column is a measurement rather than an assumption:
      !     MEASURED on `mol_carrier`, structurally empty rows 0, empty
      !     columns 0, smallest |U(j,j)| of the banded LU 1.83.
      !   * THE TRUST REGION IS ASYMMETRIC IN THE DENSITY. Its radius bounds
      !     a norm of the scaled step, so a radius Delta permits
      !     n -> n exp(+/-Delta): the same radius multiplies and divides the
      !     carrier by the same factor and can never carry it to zero. The
      !     region is symmetric in the ratio and not in the density, which is
      !     the geometry a quantity spanning decades asks for, and it is why
      !     the dogleg's cut-backs against a bound disappear here (MEASURED:
      !     largest cut 2^-4 against 2^-49).
      !   * ln 0 IS NOT A NUMBER, and zero is what the outermost cell of an
      !     ionized wind holds, so the space needs a smallest representable
      !     density -- carrier_log_floor, below.
      !
      ! WHEN THE SPACE IS ARMED. The floor is the element budget the CARRIER
      ! OPERATOR froze at the iterate, and that budget does not exist until
      ! the operator has assembled once, which first happens inside the
      ! solve's own opening residual evaluation. So the opening state is
      ! packed as a density, evaluated, and the carrier slots are then
      ! rewritten as logarithms with the floor the evaluation produced
      ! (arm_carrier_log_unknown). The residual is a function of the STATE
      ! and not of the coordinates it is named in, so that opening residual
      ! stands for the armed vector as well.
      logical :: carrier_log_unknown_wanted = .true.
      logical :: carrier_unknown_is_logarithmic = .false.
      ! THE FRACTION OF ITS ELEMENT BELOW WHICH A CARRIER DENSITY IS
      ! UNOBSERVABLE, which is not a number chosen here: it is
      ! carrier_residual's own absolute floor on the row scale (1e-20 of the
      ! free density of the carrier's element, diffusive_photochemistry.f90),
      ! where the reason is written -- "a density that small cannot change
      ! any observable, so a row below it IS converged".
      real*8, parameter :: carrier_unobservable_fraction_of_element = 1.0d-20
      ! THE SMALLEST CARRIER DENSITY THE LOGARITHMIC UNKNOWN CAN NAME, cell
      ! by cell and carrier by carrier, in the code's density units.
      !
      ! WHAT IT IS: the element budget of that carrier in that cell -- the
      ! largest density the carrier can reach there, which is the free
      ! density of its element less the nuclei sitting in the stages the
      ! carrier step holds frozen, divided by the nuclei one carrier holds
      ! (carrier_headroom, hydrogen_available_to_carriers) -- times the
      ! fraction above. The row floor's reference is that same element
      ! density without either the subtraction or the division, so this
      ! floor is at or below the row's, which is the safe direction: the
      ! space never forbids a density the row measure can still resolve.
      !
      ! WHY IT IS TIED TO THE STATE. A carrier density has no physical lower
      ! bound other than zero, so the floor of the unknown space is
      ! numerical; what it may not be is a constant, because the quantity it
      ! bounds spans the column. The floor of a cell holding 1e12 cm^-3 of
      ! hydrogen is not the floor of one holding 1e-2, and a fixed number
      ! would be a hard bound in one and unreachable noise in the other.
      ! MEASURED on `mol_carrier`, N = 500: the floor runs from 5.6e-28 to
      ! 4.6e-21 in code density units, seven decades across the column,
      ! which is the spread a constant would have had to stand in for.
      !
      ! WHAT HAPPENS AT IT: it is the lower bound of the unknown in log
      ! coordinates, so a trial below it is written onto it and the
      ! active-set rule holds it there while the model still points down
      ! (species_unknowns_outside_their_bounds, fix_active_species_bounds) --
      ! the treatment an element fraction gets at 0 and at 1. Whether a solve
      ! ever reaches it is reported and not assumed: see
      ! carrier_over_its_floor_min.
      !
      ! It is fixed once per solve from the state the solve begins at: an
      ! unknown space that moved with the iterate would make the residual a
      ! function of the iteration and not of the state. It is formed for the
      ! density unknown as well, so that the two arms are measured against
      ! one yardstick.
      real*8, allocatable :: carrier_log_floor(:,:)
      ! THE BOX THE SPECIES UNKNOWNS OF THIS STEP LIVE IN, one lower and one
      ! upper bound per unknown of the flat vector, in the coordinates the
      ! unknown is carried in (so ln n for a logarithmic carrier).  It is
      ! the ONE statement of what the admissible set is: the projection of a
      ! trial, the active set of the model and the fraction-to-the-boundary
      ! rule of the finite-difference probe all read it
      ! (species_unknowns_outside_their_bounds, fix_active_species_bounds,
      ! largest_step_inside_the_species_box).  Written by
      ! freeze_species_unknown_box, where each face is stated.
      !
      ! WHY IT IS FROZEN FOR THE STEP.  Two of its faces depend on the
      ! state: a carrier cannot exceed the cell it sits in, and it cannot
      ! exceed the element it is made of.  A box that moved with the trial
      ! would make the projection a nonlinear map of the trial, so the
      ! active set the Krylov cycle is sampled under and the set the trial
      ! is written onto would be two different sets -- and the three readers
      ! above disagreed about which state the faces came from before this
      ! array existed: the projection read the trial's own density while the
      ! other two read the iterate's.
      !
      ! Unallocated means no box, which is the three-unknown route: it has
      ! no species unknown and every reader returns before asking.
      real*8, allocatable :: species_box_lo(:), species_box_hi(:)
      ! WHETHER THE SHARED ELEMENT BUDGET IS A ROW OF THE STEP OR A FACE
      ! OF THIS BOX (item N4b, decision 14 route (i)).
      !
      ! A BUDGET IS NOT A COORDINATE BOUND. Two carriers of one element
      ! spend one budget, so the feasible set of the carriers of a cell is
      !
      !     sum_i nuclei_per_particle(i,ie) n_i  <=  budget(ie)
      !
      ! a half-space of the cell's unknowns and not a product of intervals.
      ! Each carrier's own ceiling budget/nuclei is one CORNER of it, and
      ! the conjunction of those corners admits a cell holding two hydrogen
      ! nuclei where the cell has one (item N4a, measured breach 1.0 of the
      ! available hydrogen). The corners are also what cost the step whole
      ! directions: a component with no room between its value and its own
      ! ceiling makes the fraction-to-the-boundary rule return zero for the
      ! whole direction, however much room the other components have.
      !
      ! So with element_constraint_rows_on the box carries only the faces
      ! that ARE coordinate bounds -- a density is non-negative and a
      ! carrier cannot be more of the cell than the cell is -- and the
      ! shared budget is carried as a linearized row of the step
      ! (freeze_element_constraint_rows): the step is projected onto the
      ! rows that are active, so a blocked component is moved along the
      ! constraint instead of being zeroed. Nothing is lost by removing the
      ! corners: with every carrier density non-negative the shared row
      ! implies each of them.
      !
      ! EXHALE_ELEMENT_CONSTRAINT_ROWS=0 is the measurement arm: the budget
      ! goes back under a coordinate face and no row is formed, which is the
      ! box the B5j to B5l measurements were made in.
      logical :: element_constraint_rows_on = .true.
      ! WHETHER THE WRITE-BACK RESTORES THE ELEMENT TOTALS OF THE CELL IT
      ! WRITES (item N4b deliverable 4). It is a conservation statement and
      ! not an option: an operator that moves one element of a cell without
      ! the other changes the composition the run was given.
      ! EXHALE_ELEMENT_WRITE_BACK=0 is the arm it is measured against, in
      ! which the carrier column is written and nothing puts the element
      ! total back, and the He/H of the state then leaves its reservoir
      ! along the solve (item N4a section 6b).
      logical :: element_write_back_conserves = .true.
      ! ONE FLOATING-POINT STEP INSIDE THE ELEMENT BUDGET, and it is a
      ! representability guard rather than a margin.  The face is stored in
      ! log coordinates for a logarithmic carrier, so the density the screen
      ! of eval_residual compares against the budget is exp(ln(budget)),
      ! which may land one unit in the last place ABOVE the budget; a cell
      ! written exactly onto the face would then be counted as outside it and
      ! the projected trial refused.  Eight epsilons cover the rounding of
      ! log and exp together.
      real*8, parameter :: budget_face_inside = 1.0d0 - 8.0d0*epsilon(1.0d0)
      ! Cells whose carrier already sits ABOVE its own element budget at the
      ! iterate the box was frozen at, summed over the outer iterations of a
      ! solve.  The budget face is raised to the iterate in such a cell (see
      ! freeze_species_unknown_box), so this is how often that had to happen.
      integer :: n_budget_face_at_the_iterate = 0

      ! ---- THE SHARED ELEMENT CONSTRAINT AS ROWS OF THE STEP (N4b) ----
      ! One row for each element a carrier of the cell holds: hydrogen,
      ! oxygen and carbon. Helium is held by no carrier and has no row.
      integer, parameter :: n_element_constraint_rows = 3
      integer, parameter ::                                              &
         element_constraint_element(n_element_constraint_rows) =          &
            [ ien_H, ien_O, ien_C ]
      ! c of each row at the iterate [cm^-3], and the budget it is measured
      ! against [cm^-3], frozen with the box.
      real*8, allocatable :: elem_row_c(:,:)
      real*8, allocatable :: elem_row_budget(:,:)
      ! WHAT A CANDIDATE MAY DEMAND: the budget, or the iterate's own demand
      ! where the iterate already stands above it. The budget depends on the
      ! very unknowns it constrains and is frozen at the iterate, so the
      ! iterate itself can stand outside it; a target that cut through the
      ! point the step starts from would leave the restoration nothing to
      ! restore to. It is the statement freeze_species_unknown_box makes
      ! about a raised face, made about a shared row.
      real*8, allocatable :: elem_row_target(:,:)
      ! dc/dY of each row over the unknowns of its own cell
      ! [cm^-3 per unit unknown], and the same row in the coordinates the
      ! linear algebra works in (D*dc/dY), which is where the projection
      ! acts because that is where the trust region measures length.
      real*8, allocatable :: elem_row_grad(:,:,:)
      real*8, allocatable :: elem_row_grad_scaled(:,:,:)
      ! The orthonormal basis of the ACTIVE rows of each cell, in the scaled
      ! coordinates and with the held coordinate bounds already removed, so
      ! that removing each vector in turn is the exact projection onto the
      ! intersection of the active constraint hyperplanes.
      real*8, allocatable :: elem_row_basis(:,:,:)
      integer, allocatable :: n_elem_row_basis(:)
      ! The nuclei of each element one particle of each f_sp species holds,
      ! read once from the map so that the restoration of one cell is a
      ! table lookup and not a search of the species table.
      integer :: elem_nu_of_species(n_species,n_inv_element) = 0
      logical :: elem_nu_table_known = .false.
      logical :: elem_species_is_carried(n_species) = .false.
      ! WHICH SPECIES MAY TAKE BACK THE REMAINDER OF AN ELEMENT: the ones
      ! that hold nuclei of THAT element and of no other element a carrier
      ! constrains, and that the solve does not carry. A species holding two
      ! constrained elements at once cannot be rescaled onto one of them
      ! without moving the other, so its nuclei are held fixed with the
      ! locked hydrogen of HeH+ and the remainder is shared over the rest.
      logical :: elem_species_absorbs(n_species,n_inv_element) = .false.
      logical :: elem_rows_frozen = .false.
      logical :: elem_rows_active_set_frozen = .false.
      ! HOW NEARLY ACTIVE IS ACTIVE. A row within this fraction of its own
      ! budget of being tight is treated as active, for the reason the
      ! epsilon-active set of a coordinate bound exists: a model built from
      ! finite differences cannot tell a point on a constraint from a point
      ! one probe step inside it. The probe step of a carrier unknown is
      ! sqrt(eps) of its own scale, so the band is that same order.
      real*8, parameter :: element_constraint_active_band = 1.0d-8
      ! WHEN A CANDIDATE ENTERS A RESTORATION STEP: a violation above this
      ! fraction of what the row allows. Below it the violation is inside
      ! the rounding of the demand sum and the budget sum, which are two
      ! sums of the same nuclei, and restoring would move a state by noise.
      real*8, parameter :: element_constraint_restoration_trigger         &
                         = 1.0d-10
      ! A row whose direction is already spanned by the rows of its own cell
      ! carries no new constraint, and a basis built from it would divide by
      ! its round-off. Counted over the solve: it is the degeneracy the
      ! shared H and O budgets of OH and H2O can produce.
      integer :: n_elem_row_dropped_for_rank = 0
      ! Rows whose iterate already stands outside its own frozen budget, and
      ! restoration steps taken, summed over the solve.
      integer :: n_elem_row_target_at_the_iterate = 0
      integer :: n_elem_restoration_cells = 0
      integer :: n_elem_restoration_steps = 0
      real*8  :: elem_restoration_worst_before = 0.0d0
      real*8  :: elem_restoration_worst_after  = 0.0d0
      ! The write-back's own statement: cells in which the carriers the
      ! Newton set demanded more nuclei of an element than the cell holds,
      ! so that the closure species could not absorb the remainder and the
      ! element total of the written cell is short of the total it was
      ! formed from. It is the one place a write-back can create or destroy
      ! nuclei, and it is counted rather than silent.
      ! Components a finite-difference probe could not move at any step
      ! length, summed over the probes of a solve: they are the ones sitting
      ! exactly on a face with the direction pointing out of it.
      integer :: n_probe_blocked_components = 0
      integer :: n_write_back_element_deficit = 0
      real*8  :: write_back_worst_element_deficit = 0.0d0
      ! THE SMALLEST RATIO OF AN ADOPTED CARRIER DENSITY TO ITS OWN FLOOR
      ! over every state the solve adopted, and how many (state, unknown)
      ! pairs sat at or below that floor. Above one, and a count of zero, is
      ! the statement that the floor was never reached and the unknown space
      ! never had to bound anything.
      real*8  :: carrier_over_its_floor_min = huge(1.0d0)
      integer :: n_carrier_at_its_floor = 0
      ! THE LARGEST RATIO OF AN ADOPTED CARRIER DENSITY TO ITS OWN ELEMENT
      ! BUDGET, and how many (state, unknown) pairs stood above that budget.
      ! A carrier cannot hold more nuclei of its element than the cell has
      ! free, so at or below one, and a count of zero, is the statement that
      ! the element budget held over every state the solve adopted. The
      ! reference is the budget itself, not the face of the box.
      real*8  :: carrier_over_its_budget_max = 0.0d0
      integer :: n_carrier_above_its_budget = 0
      ! The residual tolerance of the solve now running, which is the number
      ! the ACCEPTANCE GATE reads |R| against (steady_gates_met). It is the
      ! reference of the plain row measure that rmax(1) of
      ! row_maxima_of_the_certification reports, and of nothing else: the
      ! judged distance divides each hydrodynamic row by the tolerance that
      ! row carries at the cell instead
      ! (hydrodynamic_distance_from_certification_by_cell). Set at the top of a
      ! solve; negative means no solve is running.
      real*8  :: resid_tol_of_solve = -1.0d0

      ! Explicit interfaces for the external LAPACK banded-LU routines used by
      ! the direct/preconditioned Newton solves below (double-precision,
      ! general band form). The array dummies are assumed-size so that the
      ! call sites, which pass the right-hand side as a rank-1 vector
      ! (nrhs = 1), match without any argument change; this only gives the
      ! compiler kind/rank/intent information and does not alter the call.
      interface
         subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            integer,          intent(in)    :: m, n, kl, ku, ldab
            real*8,           intent(inout) :: ab(ldab,*)
            integer,          intent(out)   :: ipiv(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrf

         subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            character(1),     intent(in)    :: trans
            integer,          intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real*8,           intent(in)    :: ab(ldab,*)
            integer,          intent(in)    :: ipiv(*)
            real*8,           intent(inout) :: b(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrs

         ! The reciprocal 1-norm condition number of a factored general
         ! band, from Hager's estimator: it needs the 1-norm of the
         ! UNFACTORED matrix beside the factors, so a caller has to form
         ! that norm itself. Read only by the element-row diagnostic.
         subroutine dgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond,  &
                           work, iwork, info)
            character(1),     intent(in)    :: norm
            integer,          intent(in)    :: n, kl, ku, ldab
            real*8,           intent(in)    :: ab(ldab,*)
            integer,          intent(in)    :: ipiv(*)
            real*8,           intent(in)    :: anorm
            real*8,           intent(out)   :: rcond
            real*8,           intent(out)   :: work(*)
            integer,          intent(out)   :: iwork(*)
            integer,          intent(out)   :: info
         end subroutine dgbcon
         ! Dense eigenvalues (and right eigenvectors) of the small
         ! Hessenberg the Arnoldi recursion of the spectrum hook leaves;
         ! the Ritz values of the preconditioned operator are exactly its
         ! eigenvalues.
         subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl,     &
                          vr, ldvr, work, lwork, info)
         character*1, intent(in) :: jobvl, jobvr
         integer, intent(in)     :: n, lda, ldvl, ldvr, lwork
         double precision, intent(inout) :: a(lda,*)
         double precision, intent(out)   :: wr(*), wi(*)
         double precision, intent(out)   :: vl(ldvl,*), vr(ldvr,*)
         double precision, intent(out)   :: work(*)
         integer, intent(out)    :: info
         end subroutine dgeev
      end interface

      contains

      ! ------------------------------------------------------!

      subroutine set_transported_species_rows(on)
      ! Register a row and an unknown for every transported balance this
      ! configuration activates, and set the band geometry that goes with
      ! them. Called once before a solve; with `on` false the registry is
      ! empty and the 3-unknown path is untouched, which is what keeps every
      ! run without a transported balance bit-for-bit unchanged.
      !
      ! The transported set is the one carrier_set_init fixed once the keys
      ! were parsed (H2 always; OH, H2O and CO under the oxygen chemistry;
      ! H+ under the ionization transport), read here through
      ! carrier_solved so that the stationary system and the transport
      ! operator cannot disagree about which balances exist.
      logical, intent(in) :: on
      integer :: ic, im
      nspec_row = 0
      srow_kind = 0;  srow_idx = 0;  srow_isp = 0
      if (on) then
         ! THE ELEMENTS FIRST.  Their write-back is a projection of the
         ! whole species vector onto new element totals, so it has to run
         ! before a carrier slot is written or the carrier the Newton set
         ! would be rescaled by it.
         if (he_diffusion .and. thereis_He) then
            nspec_row = nspec_row + 1
            srow_kind(nspec_row) = srow_element_he
            srow_idx (nspec_row) = 1
         endif
         if (he_diffusion .and. he_metal_diffusion .and. thereis_metals) &
         then
            do im = 1, n_melem
               ! AN ELEMENT WITH NO RESERVOIR IS NOT A DEGREE OF FREEDOM OF
               ! THE ATMOSPHERE. Its mass fraction is identically zero in
               ! every cell, its transported balance is 0 = 0, and its
               ! unknown sits exactly on the lower bound of its own box with
               ! no room below it -- which is not a harmless extra column:
               ! the fraction-to-the-boundary rule returns zero for ANY
               ! direction whose component there points below zero
               ! (largest_step_inside_the_species_box), so the whole Krylov
               ! direction is refused and the solve has no step at all.
               ! MEASURED on the atomic element reload of
               ! backup/regression/atomic_elem_newton, whose metals.inp
               ! carries seven of the ten elements: the three absent ones
               ! put 1500 of the unknowns on a bound, the first Krylov cycle
               ! managed 0 products of 40, the initial radius came out zero
               ! and eleven outer iterations reported a zero dogleg step at
               ! ||R|| 1.888 (report B5m; item N1). The same test that the
               ! energy grid and the restart reader already use for "this
               ! element is not in this atmosphere" (set_energy_vectors,
               ! load_IC).
               if (allocated(melem_ab)) then
                  if (melem_ab(im) .le. 0.0d0) cycle
               endif
               nspec_row = nspec_row + 1
               srow_kind(nspec_row) = srow_element_trace
               srow_idx (nspec_row) = im
            enddo
         endif
         do ic = 1, merge(n_carrier, 0, thereis_mol .and.               &
                                          carrier_transport)
            if (.not. carrier_solved(ic)) cycle
            nspec_row = nspec_row + 1
            srow_kind(nspec_row) = srow_carrier
            srow_idx (nspec_row) = ic
            srow_isp (nspec_row) = carrier_species_index(ic)
         enddo
      endif
      if (.not. allocated(erow_he))                                      &
         allocate(erow_he(1:N), escale_he(1:N),                          &
                  erow_tr(1:N,n_melem), escale_tr(1:N,n_melem))
      erow_he = 0.0d0;  escale_he = 1.0d0
      erow_tr = 0.0d0;  escale_tr = 1.0d0
      erow_tr_carried = .false.
      srow_base_value = 0.0d0
      band_emptiness_stated = .false.
      nvar_jac      = 3 + nspec_row
      kl_jac        = 3*nvar_jac - 1
      ku_jac        = kl_jac
      ncolor_jac    = kl_jac + ku_jac + 1
      kl_jac_normal = 2*kl_jac
      ku_jac_normal = 2*kl_jac
      ldab_normal   = 2*kl_jac_normal + ku_jac_normal + 1
      end subroutine set_transported_species_rows

      ! ------------------------------------------------------!

      ! ------------------------------------------------------!

      logical function carrier_unknown_on()
      carrier_unknown_on = (nspec_row .gt. 0)
      end function carrier_unknown_on

      ! ------------------------------------------------------!

      integer function n_species_rows()
      n_species_rows = nspec_row
      end function n_species_rows

      ! ------------------------------------------------------!

      integer function species_row_carrier_index(i)
      ! The carrier this species row carries, or 0 if the row is not a
      ! carrier row.
      integer, intent(in) :: i
      species_row_carrier_index = 0
      if (i .lt. 1 .or. i .gt. nspec_row) return
      if (srow_kind(i) .eq. srow_carrier)                                &
         species_row_carrier_index = srow_idx(i)
      end function species_row_carrier_index

      ! ------------------------------------------------------!

      subroutine set_base_fix(Y, nfix)
      ! Anchor the first nfix physical cells at their current values.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      integer, intent(in) :: nfix
      nfix_base = max(0, min(nfix, N-10))
      if (.not. allocated(Yfix_base)) allocate(Yfix_base(nvar_jac*N))
      Yfix_base = Y
      end subroutine set_base_fix

      ! ------------------------------------------------------!

      double precision function element_unknown_column_scale(x, reservoir)&
                       result(d)
      ! THE COLUMN SCALE OF ONE ELEMENT UNKNOWN: a mass fraction, so its own
      ! magnitude where the cell holds any of the element, and the reservoir
      ! the column is anchored to where it holds none.
      !
      ! THE FLOOR SETS THE PROBE STEP OF THAT UNKNOWN'S JACOBIAN COLUMN.
      ! build_banded_jac_full perturbs the unknown by sqrt(eps) times this
      ! scale; the transport row of a cell responds with about
      ! escale * h / X_reservoir, since escale/X_reservoir is the rate at
      ! which the row moves that element; and the row is assembled to the
      ! relative accuracy of double precision. So the column is resolvable
      ! only while h stands well above eps*X_reservoir. With the relative
      ! floor below, h = sqrt(eps)*1e-4*X_reservoir = 1.5e-12 X_reservoir,
      ! about seven thousand times that round-off, which leaves three to
      ! four digits in the column entry.
      !
      ! MEASURED with the absolute floor of 1e-20 this replaced, on the
      ! atomic element reload at outer iteration 11: the carbon unknown of
      ! cell 500 stood at that floor, its probe step was 1.5e-28, its
      ! Jacobian column came back EXACTLY ZERO, the band's smallest pivot
      ! was 3.7e-19 at that same unknown, and dgbcon reported a reciprocal
      ! condition estimate of 4.0e-27. The Krylov cycle then returned a step
      ! of norm 2.5e-15 inside a ball of radius 19.6 and reported a relative
      ! residual of 0.076 for it, and the solve stopped with no descent. An
      ! empty column is an unknown no equation responds to, and a floor that
      ! is not a scale of the problem is how it was made empty.
      !
      ! THE FLOOR IS GENEROUS where a cell holds a little of the element: it
      ! overrides the cell's own value below 1e-4 of the reservoir. That
      ! costs nothing in conditioning, because the columns are equilibrated
      ! after the row scaling (equilibrate_the_scaled_columns). What the
      ! floor has to do is make the column non-empty; what makes the columns
      ! comparable is the equilibration.
      !
      ! The absolute 1e-20 remains underneath, for a row with no reservoir
      ! recorded yet: it keeps the step finite and nothing more.
      real*8, intent(in) :: x, reservoir
      d = max(abs(x), element_column_scale_floor*abs(reservoir), 1.0d-20)
      end function element_unknown_column_scale

      ! ------------------------------------------------------!

      integer function species_row_base_equation(i) result(eqn)
      ! THE EQUATION SPECIES ROW i CARRIES AT THE FIRST PHYSICAL CELL, and
      ! the one definition of it: eval_residual assembles this equation
      ! there and cell_row_scales nondimensionalizes THAT equation, so the
      ! two cannot drift apart.
      !
      ! AN ELEMENT ROW IS A DIRICHLET RESERVOIR CONDITION at the base. The
      ! element operator poses no transport balance in cell 1 -- it is the
      ! reservoir the column is anchored to, and element_transport_residual
      ! reports the cell's row as zero against the floor (READ,
      ! src/modules/functions/binary_element_diffusion.f90, the header of
      ! that routine) -- so the equation of that unknown is
      ! X(1) = X_reservoir and the scale that makes it dimensionless is the
      ! unknown's own.
      !
      ! A CARRIER ROW IS A CONTROL-VOLUME BALANCE at the base, as in every
      ! other cell. The inner ghost is data and not an unknown: it carries
      ! the composition of the handoff, and where the base face flows
      ! inward the advective term it contributes is a constant of the row
      ! (READ, src/modules/lower_atmosphere/diffusive_photochemistry.f90,
      ! the j == 1 branch of the advective block of carrier_residual). The
      ! cell balances transport against chemistry there, so the scale is
      ! the row's own terms, the same quantity as in every other cell.
      integer, intent(in) :: i
      if (i .ge. 1 .and. i .le. nspec_row) then
         if (srow_kind(i) .eq. srow_carrier) then
            eqn = base_row_control_volume_balance
         else
            eqn = base_row_reservoir_condition
         endif
      else
         eqn = base_row_control_volume_balance
      endif
      end function species_row_base_equation

      ! ------------------------------------------------------!

      subroutine cell_state_scales(Y, D)
      ! Diagonal scaling of the Newton system and of the line-search merit:
      ! the characteristic scale of each conserved unknown in its own cell.
      ! The three numbers come from state_scales_of_cell (steady_residual.f90).
      ! Since section 143 they are NOT what the convergence measure divides by:
      ! each residual row is divided by the largest term that row itself
      ! contains (residual_row_scale), so the merit and the acceptance test are
      ! no longer one expression apart. That is deliberate and it was measured
      ! -- scaling this system by those quantities as well, the mass row's
      ! varying by a factor 9 across the first two cells, left the molecular
      ! hot-Uranus solve with no descent direction after 179 iterations against
      ! 9 for the measure-only build (docs/p54_base_layer_mass_flux.md section
      ! 10.4). The unknowns keep the scales below; the rows do not.
      !
      ! It sets the finite-difference step sizes, the scaled Newton system
      ! D^-1 J D, and the line-search merit ||D^-1 F||_2.
      real*8, dimension(nvar_jac*N), intent(in)  :: Y
      real*8, dimension(nvar_jac*N), intent(out) :: D
      real*8  :: vv, cs
      integer :: j, i1, i2, i3, i, is
      do j = 1, N
         i1 = nvar_jac*(j-1) + 1;  i2 = i1 + 1;  i3 = i1 + 2
         call state_scales_of_cell(j, Y(i1), Y(i2), Y(i3),               &
                                   D(i1), D(i2), D(i3), vv, cs)
         ! A SPECIES UNKNOWN'S COLUMN SCALE is its own magnitude with the
         ! floor diffusive_photochemistry builds for it -- the same
         ! construction as the momentum unknown's |rho v| + rho c_s, and for
         ! the same reason. Why it is not the ROW scale, and what happens if
         ! it is, is recorded where the two are defined.
         do i = 1, nspec_row
            is = i3 + i
            if (srow_kind(i) .eq. srow_carrier) then
               if (carrier_unknown_is_logarithmic) then
                  ! A LOGARITHM IS ITS OWN SCALE: a unit change of ln n is a
                  ! factor e in n at every density, so the column needs no
                  ! scale from the state. This 1 is what makes the banded
                  ! Jacobian's column step and the Krylov probe step
                  ! RELATIVE changes of the density, and it is what makes the
                  ! trust-region radius a factor on n rather than an amount
                  ! of it; both are stated at the declaration of
                  ! carrier_unknown_is_logarithmic, together with the factor
                  ! n the Jacobian column then carries.
                  D(is) = 1.0d0
               else
                  D(is) = max(carrier_column_scale(j, srow_idx(i))/n0,   &
                              1.0d-30*max(abs(Y(is)), 1.0d0))
               endif
            else
               D(is) = element_unknown_column_scale(Y(is),               &
                                                    srow_base_value(i))
            endif
         enddo
      enddo
      ! KEPT ASIDE AS THE STATE SCALE. The column scaling the linear algebra
      ! works in may be equilibrated away from this array
      ! (equilibrate_the_scaled_columns); everything that needs the
      ! characteristic magnitude of an unknown reads this one.
      if (.not. allocated(state_column_scale))                           &
         allocate(state_column_scale(nvar_jac*N))
      if (size(state_column_scale) .ne. nvar_jac*N) then
         deallocate(state_column_scale)
         allocate(state_column_scale(nvar_jac*N))
      endif
      state_column_scale = D
      end subroutine cell_state_scales

      ! ------------------------------------------------------!

      double precision function merit_row_scale_from_certification(      &
                                 cert_scale, to_code_units) result(s)
      ! THE DENOMINATOR ONE ROW IS READ THROUGH, AND THE ONLY DEFINITION OF
      ! IT ON THE SPECIES-ROW ROUTE.
      !
      ! cert_scale is the scale the certification divides that row by --
      ! residual_row_scale for a hydrodynamic row, escale_he or escale_tr
      ! for an element row, the row's own physical terms for a carrier row
      ! -- and to_code_units is the factor that writes the row per code time
      ! in code density units, so that the quotient is the SAME number
      ! certification_row_measure forms. The floor is the certification's
      ! own (cert_scale_floor), so a cell whose terms all vanish reads zero
      ! instead of dividing by zero, exactly as it does in the measure.
      !
      ! Row scaling for the linear algebra and certification normalization
      ! are ONE thing here: the merit the step control descends is then the
      ! 2-norm of the quantities the gate takes a maximum of, and descending
      ! the merit is descending the gate.
      real*8, intent(in) :: cert_scale, to_code_units
      s = max(cert_scale, cert_scale_floor)*to_code_units
      end function merit_row_scale_from_certification

      ! ------------------------------------------------------!

      subroutine cell_row_scales(Y, D, Dr, u)
      ! ROW scaling of the Newton system, beside the column scaling D.
      !
      ! ON THE SPECIES-ROW ROUTE EVERY ROW IS READ ON THE SCALE THE
      ! CERTIFICATION READS IT ON (merit_row_scale_from_certification), so
      ! that the merit || F/Dr ||_2 is the 2-norm of the very quantities the
      ! gate takes a maximum of. Until this was so the two were different
      ! functionals on different scalings and a rule that bounded the first
      ! bounded nothing about the second: MEASURED on the atomic element
      ! reload at outer iteration 60, the eight element rows held 0.437 of
      ! the squared merit and nothing of ||R||, and accepted steps swung
      ! ||R|| by 40 percent while the merit moved 1e-4 of itself
      ! (N7 report, item 7).
      !
      ! THE PRICE IS THAT THE SCALED OPERATOR CHANGES WITH IT. Dr divides
      ! the rows of the banded Jacobian and of the matrix-free product as
      ! well as the residual, so the preconditioner, its condition estimate
      ! and the trust-region model are all built on this scaling. MEASURED,
      ! item N20: the carrier reload's best ||R|| falls 1.151 -> 0.446, and
      ! the ATOMIC ELEMENT reload stops with no descent at ||R|| = 1.777
      ! against the 7.4e-2 the state scaling reaches, its row scales
      ! spanning 3.7e-10 to 92.7 and dgbcon 7.3e-8 -> 8.9e-10 -- the failure
      ! docs/p54_base_layer_mass_flux.md section 10.4 measured for this
      ! scaling, with a corrected trust region. The freedom left is the
      ! COLUMN scaling D, which is not this routine's.
      !
      ! THE THREE-UNKNOWN ROUTE IS UNTOUCHED, and by construction rather
      ! than by branch: Dr is a copy of D and the routine returns before it
      ! reads anything else whenever no species row is registered. Every
      ! site that divides by the row scale keeps the original expression,
      ! character for character, in the three-unknown branch.
      !
      ! u IS THE CONSERVED STATE Y DESCRIBES, and the hydrodynamic rows can
      ! be put on their certification scale only with it: residual_row_scale
      ! is the largest term the row itself holds, which is a property of an
      ! assembled state and not of a vector of unknowns. A caller that
      ! cannot supply it (the replay driver of the main program, and the
      ! beyond-band column probe) gets the state scales in those three
      ! slots, which is what the row scaling was before this change.
      !
      ! IT IS A LOCAL OF THE SOLVE AND NOT A MODULE ARRAY, so that it has the
      ! same aliasing status as D and the solve owns its own copy. That alone
      ! does NOT make a three-unknown solve bit-identical to one that divides
      ! by D throughout -- measured, the two WASP-121b Newton cases moved by
      ! 1.2e-9 and 4.0e-8 -- so every site that divides by the row scale
      ! branches on nvar_jac and keeps the original expression, character for
      ! character, in the three-unknown branch. The run-time band geometry is
      ! NOT the cause of that difference and that was tested: a build with
      ! kl_jac and the rest back at compile time is bit-identical to this one.
      real*8, dimension(nvar_jac*N),  intent(in)  :: Y, D
      real*8, dimension(nvar_jac*N),  intent(out) :: Dr
      real*8, dimension(3,1-Ng:N+Ng), intent(in), optional :: u
      real*8  :: tscale_code
      integer :: j, i, is, k
      Dr = D
      if (nspec_row .le. 0) return
      tscale_code = R0/(v0*n0)
      if (present(u)) then
         do j = 1, N
            do k = 1, 3
               Dr(nvar_jac*(j-1)+k) =                                     &
                  merit_row_scale_from_certification(                     &
                       residual_row_scale(k, j, u), 1.0d0)
            enddo
         enddo
      endif
      do j = 1, N
         do i = 1, nspec_row
            is = nvar_jac*(j-1) + 3 + i
            select case (srow_kind(i))
            case (srow_carrier)
               ! The row's own largest terms without the time term, which is
               ! the quantity carrier_steady_residual divides by, converted
               ! from cm^-3 s^-1 to the code's clock and density. There is
               ! no separate base cell here, and that is a statement and not
               ! an omission: the carrier row of cell 1 is a control-volume
               ! balance (species_row_base_equation), so it is
               ! nondimensionalized by its own terms like every other cell.
               if (allocated(crow_terms_last)) then
                  Dr(is) = merit_row_scale_from_certification(            &
                              crow_terms_last(j,srow_idx(i)), tscale_code)
               else
                  Dr(is) = merit_row_scale_from_certification(            &
                              carrier_row_term_scale(j, srow_idx(i))/n0,  &
                              1.0d0)
               endif
            case (srow_element_he)
               if (j .eq. 1 .and. species_row_base_equation(i) .eq.       &
                                  base_row_reservoir_condition) then
                  ! THE BASE ROW IS THE RESERVOIR STATEMENT and not a
                  ! transport balance: the operator poses no balance there
                  ! and the certification measures none, so the row is
                  ! nondimensionalized by the unknown's own scale.
                  Dr(is) = D(is)
               else
                  Dr(is) = merit_row_scale_from_certification(            &
                              escale_he(j), elem_he_to_code())
               endif
            case default
               if (.not. erow_tr_carried(srow_idx(i))) then
                  ! No reservoir, no equation: the row says the element
                  ! stays absent and it is already dimensionless.
                  Dr(is) = 1.0d0
               else if (j .eq. 1 .and. species_row_base_equation(i) .eq.  &
                                       base_row_reservoir_condition) then
                  Dr(is) = D(is)
               else
                  Dr(is) = merit_row_scale_from_certification(            &
                              escale_tr(j,srow_idx(i)), elem_tr_to_code())
               endif
            end select
         enddo
      enddo
      end subroutine cell_row_scales

      ! ------------------------------------------------------!

      subroutine equilibrate_the_scaled_columns(ab, D)
      ! EVERY COLUMN OF THE SCALED OPERATOR MADE THE SAME SIZE.
      !
      ! ab enters holding A = Drow^-1 J D, the operator the banded
      ! factorization, the Krylov space and the trust-region model are all
      ! built on. Each column is multiplied by a factor e that brings its
      ! largest entry to within a factor sqrt(2) of one, and D is multiplied
      ! by the same factor, so that ab leaves holding Drow^-1 J (D e) and the
      ! caller's D is the column scaling of what ab holds. This is the
      ! column half of the two-sided equilibration LAPACK's dgeequ forms;
      ! the row half is not free here, because the row scaling is the
      ! certification's (cell_row_scales).
      !
      ! WHY THE COLUMNS AND NOT THE ROWS. The residual of a row is divided
      ! by the scale the gate reads it on, and that is what makes the merit
      ! the step control descends the functional the state is judged by. The
      ! spread of those scales is a property of the physics -- MEASURED on
      ! the atomic element reload, 3.7e-10 to 92.7, eleven decades -- and it
      ! is carried into the operator whether or not anything is done about
      ! it. What can absorb it is the column scaling, which appears in the
      ! model only through the change of variables s -> e^-1 s.
      !
      ! IT CHANGES THE LINEAR ALGEBRA AND NOT THE FUNCTIONAL. The merit is
      ! 0.5 ||r0 + A s||^2 with r0 = F/Drow, which holds no D at all, and
      ! the step in unknown space is dY = D s: replacing D by D e and s by
      ! e^-1 s leaves both the residual and the trial state unchanged. What
      ! moves is the shape of the trust-region ball, the Krylov subspace and
      ! the pivoting of the band.
      !
      ! THE FACTOR IS A POWER OF TWO, so D*e, the scaled columns and the
      ! division by e are exact in binary floating point and the
      ! equilibration adds no rounding of its own to the operator.
      !
      ! A column with no finite positive entry keeps e = 1: an unresolved
      ! color leaves its column at exactly zero (build_banded_jac_full) and
      ! there is nothing to equilibrate.
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(inout) :: ab
      real*8, dimension(nvar_jac*N),                   intent(inout) :: D
      integer :: neq, jcol, irow, ilo, ihi
      real*8  :: colmax, e, aij
      neq = nvar_jac*N
      if (.not. allocated(scaled_column_norm_before)) then
         allocate(scaled_column_norm_before(neq))
      else if (size(scaled_column_norm_before) .ne. neq) then
         deallocate(scaled_column_norm_before)
         allocate(scaled_column_norm_before(neq))
      endif
      colequil_min = huge(1.0d0);  colequil_max = 0.0d0
      colnorm_before_min = huge(1.0d0);  colnorm_before_max = 0.0d0
      colnorm_after_min  = huge(1.0d0);  colnorm_after_max  = 0.0d0
      do jcol = 1, neq
         ilo = max(1,   jcol - ku_jac)
         ihi = min(neq, jcol + kl_jac)
         colmax = 0.0d0
         do irow = ilo, ihi
            aij = ab(kl_jac+ku_jac+1 + irow - jcol, jcol)
            if (finite_real(aij)) colmax = max(colmax, abs(aij))
         enddo
         scaled_column_norm_before(jcol) = colmax
         colnorm_before_min = min(colnorm_before_min, colmax)
         colnorm_before_max = max(colnorm_before_max, colmax)
         if (colmax .gt. 0.0d0) then
            e = 2.0d0**(-nint(log(colmax)/log(2.0d0)))
         else
            e = 1.0d0
         endif
         colequil_min = min(colequil_min, e)
         colequil_max = max(colequil_max, e)
         if (e .ne. 1.0d0) then
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                 &
                  ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*e
            enddo
            D(jcol) = D(jcol)*e
         endif
         colnorm_after_min = min(colnorm_after_min, colmax*e)
         colnorm_after_max = max(colnorm_after_max, colmax*e)
      enddo
      end subroutine equilibrate_the_scaled_columns

      ! ------------------------------------------------------!

      subroutine equilibrate_the_preconditioner(abf)
      ! THE BAND THAT IS FACTORIZED, EQUILIBRATED ON BOTH SIDES.
      !
      ! abf enters holding the shifted scaled system the Krylov cycle is
      ! preconditioned with. Its rows carry the certification scaling, which
      ! is fixed by the gate, and its columns the state scales of the
      ! unknowns; both spans are properties of the problem and neither is
      ! free in the MODEL. In the PRECONDITIONER they are: an approximate
      ! inverse rescaled on either side is another approximate inverse, so
      ! abf is replaced by R abf C with R and C diagonal, and the solve
      ! becomes C dgbtrs(R v) (banded_preconditioner_solve). The model ab,
      ! the merit, the trust-region coordinates and the radius are
      ! untouched.
      !
      ! R AND C ARE THE TWO STEPS OF ONE PASS: rows first, to unit largest
      ! entry, then columns of the row-scaled band, also to unit largest
      ! entry. One pass of each is what LAPACK's dgeequ does and it is
      ! enough here; both factors are powers of two, so the rescaling adds
      ! no rounding of its own to the band.
      !
      ! WHY IT IS NEEDED, MEASURED on the atomic element reload at outer
      ! iteration 11: without it the band's smallest pivot stood at 3.7e-19
      ! against a largest of 1.05e7 and dgbcon reported a reciprocal
      ! condition estimate of 4.0e-27, so the triangular solve was a
      ! division by round-off. The Krylov cycle then reported a relative
      ! residual of 0.076 for a step of norm 2.5e-15 -- a number the true
      ! operator, whose largest scaled row is 2.8e7, cannot produce for a
      ! step that short -- and the solve stopped with no descent step.
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(inout) :: abf
      integer :: neq, jcol, irow, ilo, ihi
      real*8  :: rmax(nvar_jac*N), cmax, aij
      neq = nvar_jac*N
      if (.not. allocated(precon_row_scale)) then
         allocate(precon_row_scale(neq), precon_col_scale(neq))
      else if (size(precon_row_scale) .ne. neq) then
         deallocate(precon_row_scale, precon_col_scale)
         allocate(precon_row_scale(neq), precon_col_scale(neq))
      endif
      rmax = 0.0d0
      do jcol = 1, neq
         ilo = max(1,   jcol - ku_jac)
         ihi = min(neq, jcol + kl_jac)
         do irow = ilo, ihi
            aij = abf(kl_jac+ku_jac+1 + irow - jcol, jcol)
            if (finite_real(aij)) rmax(irow) = max(rmax(irow), abs(aij))
         enddo
      enddo
      do irow = 1, neq
         if (rmax(irow) .gt. 0.0d0) then
            precon_row_scale(irow) =                                     &
               2.0d0**(-nint(log(rmax(irow))/log(2.0d0)))
         else
            precon_row_scale(irow) = 1.0d0
         endif
      enddo
      do jcol = 1, neq
         ilo = max(1,   jcol - ku_jac)
         ihi = min(neq, jcol + kl_jac)
         cmax = 0.0d0
         do irow = ilo, ihi
            abf(kl_jac+ku_jac+1 + irow - jcol, jcol) =                   &
               abf(kl_jac+ku_jac+1 + irow - jcol, jcol)                  &
               *precon_row_scale(irow)
            aij = abf(kl_jac+ku_jac+1 + irow - jcol, jcol)
            if (finite_real(aij)) cmax = max(cmax, abs(aij))
         enddo
         if (cmax .gt. 0.0d0) then
            precon_col_scale(jcol) = 2.0d0**(-nint(log(cmax)/log(2.0d0)))
         else
            precon_col_scale(jcol) = 1.0d0
         endif
         do irow = ilo, ihi
            abf(kl_jac+ku_jac+1 + irow - jcol, jcol) =                   &
               abf(kl_jac+ku_jac+1 + irow - jcol, jcol)                  &
               *precon_col_scale(jcol)
         enddo
      enddo
      preconditioner_equilibrated = .true.
      end subroutine equilibrate_the_preconditioner

      ! ------------------------------------------------------!

      subroutine unit_infinity_norm_row_scaling_of_the_band(ab)
      ! THE DIAGONAL E THAT BRINGS EVERY ROW OF THE BANDED MODEL TO UNIT
      ! INFINITY NORM, as a power of two, formed and kept but NOT applied
      ! to ab.
      !
      ! ab is left where it is because the merit's gradient is read from
      ! it (band_matvec_transpose on r0 = F/Drow) and that gradient belongs
      ! to the certification scales. What E is applied to is the band the
      ! preconditioner is factored from and, inside the Krylov cycle, the
      ! right-hand side and the operator's action: the cycle then solves
      ! (E A) s = E r0, the same linear system in a different norm.
      !
      ! A row with no finite positive entry keeps E = 1: an unresolved
      ! color leaves its column at exactly zero (build_banded_jac_full) and
      ! a row of zeros has no size to equilibrate.
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab
      integer :: neq, jcol, irow, ilo, ihi
      real*8  :: aij
      real*8, dimension(nvar_jac*N) :: rmax
      neq = nvar_jac*N
      if (.not. allocated(model_row_equilibration)) then
         allocate(model_row_equilibration(neq))
      else if (size(model_row_equilibration) .ne. neq) then
         deallocate(model_row_equilibration)
         allocate(model_row_equilibration(neq))
      endif
      rmax = 0.0d0
      do jcol = 1, neq
         ilo = max(1,   jcol - ku_jac)
         ihi = min(neq, jcol + kl_jac)
         do irow = ilo, ihi
            aij = ab(kl_jac+ku_jac+1 + irow - jcol, jcol)
            if (finite_real(aij)) rmax(irow) = max(rmax(irow), abs(aij))
         enddo
      enddo
      model_row_equil_min = huge(1.0d0);  model_row_equil_max = 0.0d0
      do irow = 1, neq
         if (rmax(irow) .gt. 0.0d0) then
            model_row_equilibration(irow) =                              &
               2.0d0**(-nint(log(rmax(irow))/log(2.0d0)))
         else
            model_row_equilibration(irow) = 1.0d0
         endif
         model_row_equil_min = min(model_row_equil_min,                  &
                                   model_row_equilibration(irow))
         model_row_equil_max = max(model_row_equil_max,                  &
                                   model_row_equilibration(irow))
      enddo
      model_rows_equilibrated = .true.
      end subroutine unit_infinity_norm_row_scaling_of_the_band

      ! ------------------------------------------------------!

      logical function the_linear_model_rows_are_equilibrated(neq)       &
               result(yes)
      ! Whether the Krylov cycle is to work in the equilibrated rows: the
      ! arm is on, a scaling was formed for the band currently factored,
      ! and it is a scaling of this many rows.
      integer, intent(in) :: neq
      yes = model_row_equilibration_on .and. model_rows_equilibrated
      if (yes) yes = allocated(model_row_equilibration)
      if (yes) yes = (size(model_row_equilibration) .eq. neq)
      end function the_linear_model_rows_are_equilibrated

      ! ------------------------------------------------------!

      subroutine hold_residual_evaluation_products(p)
      ! Everything a residual evaluation writes and a LATER reader consumes,
      ! held aside; see residual_evaluation_products for what is in the set
      ! and why. Allocation status is carried with the values, so a module
      ! array that was not allocated before the evaluation comes back
      ! unallocated.
      type(residual_evaluation_products), intent(out) :: p
      if (allocated(erow_he)) then
         p%erow_he   = erow_he
         p%escale_he = escale_he
         p%erow_tr   = erow_tr
         p%escale_tr = escale_tr
      endif
      if (allocated(crow_res_last)) then
         p%crow_res   = crow_res_last
         p%crow_terms = crow_terms_last
      endif
      p%erow_tr_carried    = erow_tr_carried
      p%carrier_relnorm    = carrier_relnorm_last
      p%carrier_cellmax    = carrier_cellmax_last
      p%carrier_cell       = carrier_cell_worst
      p%v_headroom_iterate = v_headroom_iterate
      p%n_headroom_last    = n_headroom_last
      p%n_eq_sweeps_last   = n_eq_sweeps_last
      p%n_no_chem_root_last = n_no_chem_root_last
      call save_carrier_module_state(p%carrier)
      end subroutine hold_residual_evaluation_products

      ! ------------------------------------------------------!

      subroutine hold_solve_refusal_statistics(t)
      ! See solve_refusal_statistics.
      type(solve_refusal_statistics), intent(out) :: t
      t%n_trial_headroom      = n_trial_headroom
      t%n_trial_no_chem_root  = n_trial_no_chem_root
      t%no_chem_root_cell_max = no_chem_root_cell_max
      t%no_chem_root_res_max  = no_chem_root_res_max
      t%n_eval_refusal        = n_eval_refusal
      end subroutine hold_solve_refusal_statistics

      ! ------------------------------------------------------!

      subroutine put_back_solve_refusal_statistics(t)
      type(solve_refusal_statistics), intent(in) :: t
      n_trial_headroom      = t%n_trial_headroom
      n_trial_no_chem_root  = t%n_trial_no_chem_root
      no_chem_root_cell_max = t%no_chem_root_cell_max
      no_chem_root_res_max  = t%no_chem_root_res_max
      n_eval_refusal        = t%n_eval_refusal
      end subroutine put_back_solve_refusal_statistics

      ! ------------------------------------------------------!

      subroutine put_back_residual_evaluation_products(p)
      ! The other half of hold_residual_evaluation_products.
      type(residual_evaluation_products), intent(in) :: p
      if (allocated(erow_he) .and. allocated(p%erow_he)) then
         erow_he   = p%erow_he
         escale_he = p%escale_he
         erow_tr   = p%erow_tr
         escale_tr = p%escale_tr
      endif
      if (allocated(crow_res_last) .and. allocated(p%crow_res)) then
         crow_res_last   = p%crow_res
         crow_terms_last = p%crow_terms
      endif
      erow_tr_carried      = p%erow_tr_carried
      carrier_relnorm_last = p%carrier_relnorm
      carrier_cellmax_last = p%carrier_cellmax
      carrier_cell_worst   = p%carrier_cell
      v_headroom_iterate   = p%v_headroom_iterate
      n_headroom_last      = p%n_headroom_last
      n_eq_sweeps_last     = p%n_eq_sweeps_last
      n_no_chem_root_last  = p%n_no_chem_root_last
      call restore_carrier_module_state(p%carrier)
      end subroutine put_back_residual_evaluation_products

      ! ------------------------------------------------------!

      subroutine apply_base_fix(Y, Fvec)
      ! Replace the residual rows of the frozen base cells by anchor rows.
      real*8, dimension(nvar_jac*N), intent(in)    :: Y
      real*8, dimension(nvar_jac*N), intent(inout) :: Fvec
      integer :: j, k
      do j = 1, nfix_base
         do k = 1, nvar_jac
            Fvec(nvar_jac*(j-1)+k) = Y(nvar_jac*(j-1)+k) - Yfix_base(nvar_jac*(j-1)+k)
         enddo
      enddo
      end subroutine apply_base_fix

      ! ------------------------------------------------------!

      integer function neq_newton()
      ! Number of Newton unknowns = 3 hydro variables over physical cells.
      neq_newton = nvar_jac*N
      end function neq_newton

      ! ------------------------------------------------------!

      subroutine pack_U(u, Y)
      ! Physical-cell conserved variables -> flat unknown vector.
      ! Y(3*(j-1)+k) = u(j,k), j = 1..N, k = 1..3.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(nvar_jac*N),         intent(out) :: Y
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Y(nvar_jac*(j-1)+k) = u(k,j)
         enddo
      enddo
      end subroutine pack_U

      ! ------------------------------------------------------!

      subroutine unpack_U(Y, u)
      ! Flat unknown vector -> physical-cell conserved variables.
      ! Ghost cells are left untouched (Apply_BC sets them).
      real*8, dimension(nvar_jac*N),         intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            u(k,j) = Y(nvar_jac*(j-1)+k)
         enddo
      enddo
      end subroutine unpack_U

      ! ------------------------------------------------------!

      subroutine pack_species_rows(u, f_sp, Y)
      ! Fill the species slots of the unknown vector from the state's own
      ! composition: Y(nvar*(j-1)+3+i) = n_i(j)/n0, IN THE CODE'S OWN DENSITY
      ! UNITS and not in cm^-3. The Newton takes norms of the whole vector --
      ! the Krylov step size is sqrt(eps)(1 + ||Y||)/||v|| -- so a component
      ! thirteen decades from the hydrodynamic ones would set that step by
      ! itself.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      real*8, dimension(nvar_jac*N),          intent(inout) :: Y
      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Yel
      integer :: j, i
      if (nspec_row .le. 0) return
      if (element_rows_registered()) call element_mass_fractions(f_sp, Yel)
      do j = 1, N
         do i = 1, nspec_row
            select case (srow_kind(i))
            case (srow_carrier)
               Y(nvar_jac*(j-1)+3+i) =                                    &
                  carrier_unknown_from_density(f_sp(j,srow_isp(i))*u(1,j), &
                                               j, i)
            case (srow_element_he)
               Y(nvar_jac*(j-1)+3+i) = Yel(j,1)
            case default
               Y(nvar_jac*(j-1)+3+i) = Yel(j,1+srow_idx(i))
            end select
         enddo
      enddo
      ! The reservoir composition this solve is to hold the base at, taken
      ! from the state the vector is packed from (see srow_base_value).
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier)                              &
            srow_base_value(i) = Y(3+i)
      enddo
      end subroutine pack_species_rows

      ! ------------------------------------------------------!

      logical function element_rows_registered()
      ! Whether any elemental balance is a row of this solve.
      integer :: i
      element_rows_registered = .false.
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier)                             &
            element_rows_registered = .true.
      enddo
      end function element_rows_registered

      ! ------------------------------------------------------!

      subroutine write_species_rows_into_composition(Y, u, f_sp)
      ! The species unknowns of Y, written into the species vector: the
      ! element totals first, through the one map an element mass fraction
      ! has back into a species vector, and the carrier columns after, so a
      ! carrier the Newton set is not rescaled by the element projection.
      !
      ! AND THE ELEMENT TOTALS OF THE WRITTEN CELL ARE THE TOTALS THE
      ! ELEMENT ROWS CARRY (item N4b deliverable 4). A carrier holds nuclei
      ! of an element that other species of the same cell hold as well, so
      ! writing a carrier column moves that element's total; with the
      ! carrier written after the element projection nothing put it back,
      ! and with helium diffusion off no elemental row would have. MEASURED
      ! inside one coupled carrier JFNK: the He/H nucleus ratio of the state
      ! left its reservoir by 1.4e-11 at the first iterate and 1.1e-1 at the
      ! thirteenth, monotonically from the third on, while every carrier
      ! budget stood at breach zero (item N4a section 6b) -- the carriers
      ! took hydrogen nuclei that no species gave up, so the cell held more
      ! hydrogen than it was formed with and its helium did not follow.
      !
      ! The remainder is restored WITHIN the cell and over that element's
      ! own species: the nuclei the carriers no longer hold go back to the
      ! species that close the element (H I, and the H II, H2+ and H3+ the
      ! step does not re-solve; the oxygen and carbon stages), in proportion
      ! to what they held, which is carrier_write_back's rule on the
      ! marching path. The reservoir ratio is never rescaled: helium is held
      ! by no carrier, so restoring the hydrogen total restores the ratio.
      ! Where the carriers demand MORE of an element than the cell holds the
      ! closure species have nothing left to give, and that cell is counted
      ! (n_write_back_element_deficit) and refused by the feasibility screen
      ! rather than written with nuclei nobody supplied.
      real*8, dimension(nvar_jac*N),          intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Yel
      real*8  :: absorb_before(n_inv_element)
      real*8  :: carried_before(n_inv_element)
      real*8  :: carried_after(n_inv_element)
      real*8  :: absorb_target, rscale, nd_cell, deficit
      integer :: j, i, ie, isp, nu
      logical :: restore_the_totals
      if (nspec_row .le. 0) return
      if (element_rows_registered()) then
         ! The LOWER ghosts keep the composition they were given: the base is
         ! the Dirichlet reservoir of the element operator (srow_base_value).
         ! The UPPER ghosts are the zero-gradient continuation of the
         ! outermost physical cell, written below.  Only the physical cells
         ! are unknowns.
         call element_mass_fractions(f_sp, Yel)
         do j = 1, N
            do i = 1, nspec_row
               ! A MASS FRACTION IS WRITTEN ONTO [0,1] HERE, and the trial
               ! that left it is refused elsewhere.  The projection divides
               ! the element totals by masses and a fraction outside the
               ! range describes a cell with negative hydrogen, which is not
               ! a state the sweep below can be asked about; the positivity
               ! screen of eval_residual has already recorded the departure,
               ! so nothing is hidden by evaluating the clamped state.
               if (srow_kind(i) .eq. srow_element_he) then
                  Yel(j,1) = min(max(Y(nvar_jac*(j-1)+3+i), 0.0d0), 1.0d0)
               else if (srow_kind(i) .eq. srow_element_trace) then
                  Yel(j,1+srow_idx(i)) =                                 &
                     min(max(Y(nvar_jac*(j-1)+3+i), 0.0d0), 1.0d0)
               endif
            enddo
         enddo
         ! THE UPPER GHOST OF A TRANSPORTED ELEMENT IS THE OUTERMOST CELL'S,
         ! AND IT IS THIS ITERATE'S.  The element operator poses a
         ! zero-gradient outflow condition at the top of the column: nothing
         ! enters the domain from outside through an element, so the ghost
         ! that closes the outer face of cell N carries the mixing ratio cell
         ! N carries, and the condition is then a condition on the state
         ! being evaluated.  It is the same rule the marching path leaves the
         ! ghost at (Yetr(N+1:N+Ng) = Yetr(N) in the advective stage of
         ! binary_element_diffusion.f90, projected by this same map).
         !
         ! With a zero-gradient element mass fraction the element MIXING
         ! RATIO is zero-gradient too: msum/n_H = m_1 + m_He (n_He/n_H) is
         ! fixed by the helium mass fraction, so equal helium and trace mass
         ! fractions at cell N and its ghost mean equal n_He/n_H and
         ! n_X/n_H, which is what the trace row reads (fX in
         ! element_transport_residual).
         do i = 1, nspec_row
            if (srow_kind(i) .eq. srow_element_he) then
               do j = N+1, N+Ng
                  Yel(j,1) = Yel(N,1)
               enddo
            else if (srow_kind(i) .eq. srow_element_trace) then
               do j = N+1, N+Ng
                  Yel(j,1+srow_idx(i)) = Yel(N,1+srow_idx(i))
               enddo
            endif
         enddo
         call project_element_mass_fractions(f_sp, Yel)
      endif
      restore_the_totals = carrier_rows_registered() .and.               &
                           element_write_back_conserves
      if (restore_the_totals) call build_the_nuclei_table
      if (restore_the_totals) call mark_the_carried_species
      do j = 1, N
         ! THE ELEMENT TOTALS THE CELL IS FORMED WITH, read before a carrier
         ! column is written, so that what the carriers take can be given
         ! back by the species that close the element.
         if (restore_the_totals) then
            nd_cell = max(u(1,j), 0.0d0)*n0
            carried_before = 0.0d0
            absorb_before  = 0.0d0
            do isp = 1, n_species
               do ie = 1, n_inv_element
                  nu = elem_nu_of_species(isp,ie)
                  if (nu .le. 0) cycle
                  if (elem_species_is_carried(isp)) then
                     carried_before(ie) = carried_before(ie)              &
                        + dble(nu)*f_sp(j,isp)*nd_cell
                  else if (elem_species_absorbs(isp,ie)) then
                     absorb_before(ie) = absorb_before(ie)                &
                        + dble(nu)*f_sp(j,isp)*nd_cell
                  endif
               enddo
            enddo
         endif
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            ! A CARRIER MASS FRACTION IS WRITTEN ONTO [0,1] HERE, for the
            ! same reason the element fractions above are and by the same
            ! rule: it is a mass fraction, its bounds are physical, and the
            ! trial that left it outside them has already been recorded by
            ! the positivity screen of eval_residual, so nothing is hidden
            ! by evaluating the clamped state. Every consumer of the carrier
            ! fraction is linear in it -- the H2 opacity and self-shielding
            ! column, the H2 share of the caloric mixture, the reaction rates
            ! of the molecular network -- so a negative value is not an
            ! instability but a negative column density and a negative
            ! reaction rate, both unphysical, and the sweep cannot be asked
            ! about that state at all. It is the enforcement of the bound,
            ! not an approximation (the same statement
            ! he_ground_singlet_density_cell makes for n(1^1S)).
            !
            ! MEASURED, why this was not left as it was: on the coupled
            ! `mol_carrier` reload the Krylov probe puts n(H2) of cell 500
            ! at -7.4e-19 of the code density unit, and the fraction formed
            ! from it reached f_sp unclamped while the element fractions
            ! beside it were clamped -- one kind of bound enforced and the
            ! other not, for the same kind of quantity (report B5g
            ! section 8).
            f_sp(j,srow_isp(i)) =                                        &
               min(max(carrier_density_from_unknown(                     &
                          Y(nvar_jac*(j-1)+3+i), j), 0.0d0), u(1,j))     &
               /max(u(1,j), 1.0d-99)
         enddo
         if (.not. restore_the_totals) cycle
         ! WHAT THE CARRIERS NOW HOLD, and the remainder the closure species
         ! of the same element have to take back.
         carried_after = 0.0d0
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            isp = srow_isp(i)
            do ie = 1, n_inv_element
               carried_after(ie) = carried_after(ie)                      &
                  + dble(elem_nu_of_species(isp,ie))*f_sp(j,isp)*nd_cell
            enddo
         enddo
         do ie = 1, n_inv_element
            if (ie .eq. ien_He) cycle
            if (carried_after(ie) .eq. carried_before(ie)) cycle
            absorb_target = absorb_before(ie) + carried_before(ie)         &
                          - carried_after(ie)
            if (absorb_before(ie) .le. 0.0d0) then
               ! Nothing of this element is left in a species that could
               ! take the remainder, so the carriers cannot be given one.
               if (carried_after(ie) .gt. carried_before(ie)) then
                  n_write_back_element_deficit =                          &
                     n_write_back_element_deficit + 1
                  write_back_worst_element_deficit =                      &
                     max(write_back_worst_element_deficit,                &
                         (carried_after(ie) - carried_before(ie))         &
                         /max(carried_after(ie), 1.0d-300))
               endif
               cycle
            endif
            if (absorb_target .lt. 0.0d0) then
               ! The carriers demand more of the element than the cell
               ! holds: the closure species go to zero and the written cell
               ! is short of the total it was formed with by this much.
               deficit = -absorb_target                                   &
                        /max(absorb_before(ie) + carried_before(ie),      &
                             1.0d-300)
               n_write_back_element_deficit =                             &
                  n_write_back_element_deficit + 1
               write_back_worst_element_deficit =                         &
                  max(write_back_worst_element_deficit, deficit)
               absorb_target = 0.0d0
            endif
            rscale = absorb_target/absorb_before(ie)
            do isp = 1, n_species
               if (.not. elem_species_absorbs(isp,ie)) cycle
               f_sp(j,isp) = f_sp(j,isp)*rscale
            enddo
         enddo
      enddo
      call set_upper_ghosts_of_the_carrier_columns(f_sp, restore_the_totals)
      end subroutine write_species_rows_into_composition

      ! ------------------------------------------------------!

      subroutine set_upper_ghosts_of_the_carrier_columns(f_sp,            &
                                                         restore_totals)
      ! THE UPPER GHOST OF A TRANSPORTED CARRIER IS THE OUTERMOST CELL'S
      ! PARTITION, AND IT IS THIS ITERATE'S.  The carrier operator poses a
      ! zero-gradient outflow condition at the top of the column -- the same
      ! fc(N+1:N+Ng) = fc(N) it writes before it measures its own residual
      ! (diffusive_photochemistry.f90) -- so the ghost that closes the outer
      ! face of cell N carries the fraction cell N carries in the state being
      ! evaluated, and no value from outside the iterate enters the row of
      ! the outermost cell.  The lower ghosts are the inflow reservoir and
      ! are not touched here.
      !
      ! The element totals of the ghost are then restored exactly as they are
      ! for a physical cell: a carrier is a partition WITHIN an element, and
      ! moving it may not move the element the ghost holds, which the
      ! zero-gradient element condition above has just written.  Every term
      ! of that restoration is proportional to the cell's density and only
      ! their ratio is used, so the ghost's own density -- which Apply_BC has
      ! not written yet at this point of a residual evaluation -- is not
      ! needed and the counts below are per unit density.
      !
      ! A GHOST DEFICIT IS NOT A REFUSAL.  A ghost carries its own density,
      ! so the same fraction is a different number of nuclei there and can
      ! ask for more of an element than the ghost holds; the carrier operator
      ! clamps that at its own outer ghosts rather than refusing the state
      ! (limit_to_element_budget of diffusive_photochemistry.f90), and the
      ! clamp below is the same statement.  A physical cell is an unknown and
      ! its deficit IS refused, above.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      logical,                                intent(in)    :: restore_totals
      real*8  :: absorb_before(n_inv_element)
      real*8  :: carried_before(n_inv_element)
      real*8  :: carried_after(n_inv_element)
      real*8  :: absorb_target, rscale
      integer :: j, i, ie, isp, nu
      if (nspec_row .le. 0) return
      if (.not. carrier_rows_registered()) return
      do j = N+1, N+Ng
         if (restore_totals) then
            carried_before = 0.0d0
            absorb_before  = 0.0d0
            do isp = 1, n_species
               do ie = 1, n_inv_element
                  nu = elem_nu_of_species(isp,ie)
                  if (nu .le. 0) cycle
                  if (elem_species_is_carried(isp)) then
                     carried_before(ie) = carried_before(ie)              &
                        + dble(nu)*f_sp(j,isp)
                  else if (elem_species_absorbs(isp,ie)) then
                     absorb_before(ie) = absorb_before(ie)                &
                        + dble(nu)*f_sp(j,isp)
                  endif
               enddo
            enddo
         endif
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            f_sp(j,srow_isp(i)) = f_sp(N,srow_isp(i))
         enddo
         if (.not. restore_totals) cycle
         carried_after = 0.0d0
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            isp = srow_isp(i)
            do ie = 1, n_inv_element
               carried_after(ie) = carried_after(ie)                      &
                  + dble(elem_nu_of_species(isp,ie))*f_sp(j,isp)
            enddo
         enddo
         do ie = 1, n_inv_element
            ! Helium is held by no carrier, so its total cannot move here.
            if (ie .eq. ien_He) cycle
            if (carried_after(ie) .eq. carried_before(ie)) cycle
            if (absorb_before(ie) .le. 0.0d0) cycle
            absorb_target = absorb_before(ie) + carried_before(ie)         &
                          - carried_after(ie)
            if (absorb_target .lt. 0.0d0) absorb_target = 0.0d0
            rscale = absorb_target/absorb_before(ie)
            do isp = 1, n_species
               if (.not. elem_species_absorbs(isp,ie)) cycle
               f_sp(j,isp) = f_sp(j,isp)*rscale
            enddo
         enddo
      enddo
      end subroutine set_upper_ghosts_of_the_carrier_columns

      ! ------------------------------------------------------!

      double precision function elem_he_to_code() result(c)
      ! The helium row is a mass rate [g cm^-3 s^-1]; this is the factor
      ! that writes it per code time in code density units, so the merit
      ! compares one clock and one density with the hydrodynamic rows.
      c = R0/(v0*n0*mu)
      end function elem_he_to_code

      ! ------------------------------------------------------!

      double precision function elem_tr_to_code() result(c)
      ! A trace element row is a mixing ratio per second; this writes it per
      ! code time.
      c = R0/v0
      end function elem_tr_to_code

      ! ------------------------------------------------------!

      subroutine unpack_species_rows(Y, u, f_sp)
      ! Inverse of pack_species_rows, onto the physical cells only.
      real*8, dimension(nvar_jac*N),          intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      call write_species_rows_into_composition(Y, u, f_sp)
      end subroutine unpack_species_rows

      ! ------------------------------------------------------!

      subroutine eval_residual(Y, f_sp_seed, f_sp, Fvec, heat, cool, n_part, &
                               admissible, may_be_adopted, n_eq_sweeps_fixed,&
                               rows_judged, resid_relnorm_judged,           &
                               state_is_discarded,                          &
                               rowmax_judged, cells_judged, which_judged)
      ! Full steady residual F(Y) with local ionization-equilibrium
      ! elimination, AND the heat/cool it used (so the caller can FREEZE the
      ! radiation when building the banded Jacobian).
      !   unpack Y -> u(1:N); Apply_BC fills ghosts;
      !   (rho,v,p) -> densities, T; refresh excited-H; ioniz_eq -> heat,cool;
      !   Apply_BC again, so the ghosts carry THAT composition;
      !   R = assemble_residual(u, n_part, heat, cool); pack R(1:N) -> Fvec.
      ! State, then the composition of that state, then the ghosts, then the
      ! fluxes -- the order the marching loop has, so that the two are the same
      ! discrete operator (sections 141.6 and 144.1).
      ! THE SEED IS AN ARGUMENT AND THE RESULT IS ANOTHER ONE.  f_sp_seed is
      ! the composition the equilibrium sweep starts from -- the cell's stored
      ! composition, or whatever initial guess the caller means -- and f_sp is
      ! where the sweep's answer is written.  They are separate dummies, and
      ! the FIRST thing this routine does is overwrite f_sp from f_sp_seed, so
      ! whatever f_sp happened to hold on entry is never read.
      !
      ! WHY THAT IS THE INTERFACE.  The sweep is a Newton solve per cell and it
      ! starts from a guess; when the guess came from whatever the caller's
      ! workspace last held, F was not a function of Y.  MEASURED on the
      ! converged 1 microbar hot Uranus (section 146.3): evaluate F at a state,
      ! at two others, then at the first again through one shared workspace and
      ! 1406 of 1500 entries differ, by 3.0e-3 of the row scale against an
      ! acceptance tolerance of 1e-5; hold the seed and the three evaluations
      ! are bitwise identical.  Every call site used to defend against that by
      ! copying the adopted composition into a workspace first, five times, by
      ! convention.  The copy is now here, once, and a caller cannot get it
      ! wrong: it has to name the seed.
      !
      ! It is the same rule the chemistry of a single cell already keeps -- a
      ! cell's root may not depend on the cell solved before it, which is what
      ! seed_species_of_cell states in constrained_chemical_equilibrium.f90 --
      ! raised one level, from the cell to the sweep.
      !
      ! n_part = n_tot + n_e is returned too, so a caller that later builds a
      ! FROZEN-radiation residual can hand the same particle count back and
      ! keep T = p/n_part (hence the transport coefficients) consistent.
      real*8, dimension(nvar_jac*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)    :: f_sp_seed
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(out)   :: f_sp
      real*8, dimension(nvar_jac*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng),            intent(out)   :: heat, cool
      real*8, dimension(1-Ng:N+Ng), optional,  intent(out)   :: n_part
      ! Whether the state Y can be described at all: .false. when the
      ! equilibrium sweep met a cell whose reaction residual is not a finite
      ! number, when any cell's composition was accepted WITHOUT being a
      ! certified root of the network (acceptance class 4, the relaxation
      ! amnesty, or class 6, the candidate that was not adopted at all), or when a row of
      ! the assembled residual is not finite. A
      ! caller that is PROBING -- a Krylov direction, a Jacobian column, a
      ! line-search trial -- must not put such a residual into its Newton
      ! model; the marching loop never asks and is unaffected
      ! (docs/Update_EXHALE_stage1.md section 121).
      logical, optional,                       intent(out)   :: admissible
      ! Whether the caller could ADOPT this state. A line-search or damped
      ! Gauss-Newton candidate could; a Jacobian color or a Krylov
      ! directional-derivative sample could not -- it is a point where the
      ! residual is read, never a state anyone keeps. The rootless-
      ! composition refusal below applies only to the first kind: a derivative
      ! sampled through a cell whose chemistry has no root is
      ! a worse derivative, but refusing it leaves the solve with no Newton
      ! model at all (13 of 17 colors zeroed, measured on the hot Uranus
      ! hand-off), and a model built from an imperfect derivative is still a
      ! model. Defaults to .true., the strict reading.
      logical, optional,                       intent(in)    :: may_be_adopted
      ! Run the composition elimination for EXACTLY this many passes, with no
      ! convergence test. A finite-difference probe must sample one map: see
      ! the elimination block below.
      integer, optional,                       intent(in)    :: n_eq_sweeps_fixed
      ! The measures the state would be JUDGED by if it were handed back
      ! (certified_row_measures). Asked for by the step controls, which have
      ! to decide whether a trial may be adopted and cannot read that off the
      ! merit. Formed here because this is where the state's own u, with the
      ! ghosts and the row terms of this evaluation, exists.
      real*8, dimension(3), optional,          intent(out)   :: rows_judged
      ! |R| of the same evaluation: the plain maximum over the three
      ! hydrodynamic rows, which is what the acceptance gate reads. The
      ! judged hydrodynamic slot above divides each of the three rows by its
      ! OWN certification tolerance first and takes the maximum of those, so
      ! the two are different numbers and need not be held by the same row;
      ! a caller that reports both asks for both.
      real*8, optional,                        intent(out)   :: resid_relnorm_judged
      ! WHETHER THE CALLER KEEPS THE STATE IT IS ASKING ABOUT. A Jacobian
      ! column, a Krylov product, a point on a trust-region ray and a
      ! self-test all read the residual at a point and then throw the point
      ! away; such an evaluation puts back the products it overwrote
      ! (residual_evaluation_products), so that the state the solve does
      ! hold keeps its own element row scales, carrier measure and carrier
      ! background. Defaults to .false., the reading a caller that adopts or
      ! may adopt the state needs.
      logical, optional,                       intent(in)    :: state_is_discarded
      ! The same three row classes WITH THE CELL each maximum sits in and
      ! the row it names (row_maxima_of_the_certification). A maximum over
      ! cells does not carry where it is, and a comparison of two
      ! evaluations of ONE state needs that, so it is formed here, beside
      ! the products of this evaluation, and not from the caller's own copy
      ! of a state whose products have since been put back.
      real*8,  dimension(3), optional,         intent(out)   :: rowmax_judged
      integer, dimension(3), optional,         intent(out)   :: cells_judged
      integer, dimension(3), optional,         intent(out)   :: which_judged

      real*8, dimension(3)           :: rmx_j
      integer, dimension(3)          :: jcl_j, iwh_j
      real*8, dimension(3,1-Ng:N+Ng) :: u, W, R
      real*8, dimension(1-Ng:N+Ng)   :: rho, v, p, T
      real*8, dimension(1-Ng:N+Ng)   :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng)   :: ne, n_tot, eta
      real*8 :: rel_change
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_sp_prev
      real*8  :: df_eq
      integer :: it_eq, n_sweep_max
      logical :: eq_fixed_count
      ! WHAT THE RESIDUAL READS THE COMPOSITION FOR, as the pass now running
      ! ENTERED: the particle count of the thermal equation of state, the H2
      ! share of it the caloric map weights its rovibrational ladder by, and
      ! the two radiative terms of the energy row. The increment test at the
      ! end of the pass is the change of these three and of nothing else
      ! (increment_of_what_the_residual_reads).
      real*8, dimension(1:N) :: npart_pass, xh2_pass, heat_pass, cool_pass
      real*8, dimension(1:N) :: npart_now, xh2_now
      logical :: radiative_field_of_previous_pass
      integer :: jworst_eq, kworst_eq
      character(len=40) :: eq_channel_txt
      type(ioniz_eq_ledger) :: sweep
      integer :: irow
      logical :: gate_chem_root, headroom_ok
      ! The shared element constraint of this candidate: the worst relative
      ! violation over the grid, and the cell and the element row that carry
      ! it (element_budget_violation_of_composition).
      real*8  :: v_headroom_last
      integer :: jheadroom, kheadroom, n_write_back_deficit_entry
      type(cert_evaluation_facts) :: cfacts
      real*8, dimension(1:N,n_carrier_max) :: cres
      ! The face mass flux the elemental rows ride on, read from the mass row
      ! this evaluation just assembled and handed to the element operator.
      real*8, dimension(1-Ng:N+Ng) :: Frho_elem
      real*8  :: crc_max, crc_vol, tscale_code, nspec_j
      integer :: jc4, jj4, ii4, crc_j, crc_ic
      logical :: ehe_ok, etr_ok
      logical :: put_products_back
      type(residual_evaluation_products) :: products_at_entry

      put_products_back = .false.
      if (present(state_is_discarded)) put_products_back = state_is_discarded
      if (put_products_back)                                             &
         call hold_residual_evaluation_products(products_at_entry)

      n_resid_eval = n_resid_eval + 1
      ! The sweep starts from the seed the caller named, and from nothing else.
      f_sp = f_sp_seed
      u = 0.0d0
      call unpack_U(Y, u)
      ! THE SPECIES UNKNOWNS ENTER HERE, and they enter as the composition
      ! the equilibrium sweep is handed. With the carrier transport on,
      ! ioniz_eq already holds each transported partition fixed at the value
      ! it is given (x_h2_fix, x_oh_fix, x_h2o_fix, x_hp_fix) and solves
      ! everything else around it -- that is what made the Picard splitting
      ! possible in the first place. The only change the coupled solve makes
      ! is WHERE those numbers come from: not from the last transport step,
      ! but from the Newton unknowns. So the sweep, the particle count it
      ! returns, the temperature p/(n_tot+n_e) that follows from it, and
      ! every heating and cooling term are all functions of the species
      ! slots as well as of the hydro rows, and the Jacobian sees them.
      !
      ! The LOWER ghost is not set from Y: its partition is incoming data
      ! (section 117 pins it from the handoff), and for an element row it is
      ! the Dirichlet reservoir the operator states. Only the physical cells
      ! are unknowns.
      !
      ! THE UPPER GHOST IS THE ITERATE'S, AND IT IS ZERO GRADIENT. The
      ! element and carrier operators pose an outflow condition at the top of
      ! the column, so the mixing ratio of every transported element and the
      ! fraction of every transported carrier in cells N+1..N+Ng are the ones
      ! cell N carries in the state being evaluated
      ! (write_species_rows_into_composition and
      ! set_upper_ghosts_of_the_carrier_columns). The outer face of cell N is
      ! then reconstructed against the iterate alone, which is what makes the
      ! outermost row an equation in its own unknown.
      headroom_ok = .true.
      n_headroom_last = 0
      v_headroom_last = 0.0d0
      jheadroom = 0;  kheadroom = 0
      n_write_back_deficit_entry = n_write_back_element_deficit
      eval_refusal = refuse_none
      eval_refusal_cell = 0;  eval_refusal_row = 0
      eval_refusal_value = 0.0d0
      n_species_at_lower_bound = 0
      if (nspec_row .gt. 0) then
         do jj4 = 1, N
         do ii4 = 1, nspec_row
            jc4  = nvar_jac*(jj4-1) + 3 + ii4
            ! THE SCREEN IS ABOUT THE DENSITY, not about the unknown that
            ! names it: with the carrier carried in logarithm the two are
            ! not the same number, and the bound belongs to the density.
            nspec_j = Y(jc4)
            if (srow_kind(ii4) .eq. srow_carrier)                         &
               nspec_j = carrier_density_from_unknown(Y(jc4), jj4)
            if (nspec_j .le. 0.0d0)                                       &
               n_species_at_lower_bound = n_species_at_lower_bound + 1
            if (nspec_j .lt. 0.0d0) then
               headroom_ok = .false.
               ! The FIRST offender is recorded: the cells are swept from
               ! the base outwards, so it is the innermost one, and a
               ! direction is refused by its innermost violation.
               if (eval_refusal .eq. refuse_none) then
                  eval_refusal       = refuse_species_unknown_negative
                  eval_refusal_cell  = jj4
                  eval_refusal_row   = ii4
                  eval_refusal_value = nspec_j
               endif
            endif
            ! THE HEADROOM IS A CARRIER CONSTRAINT.  It asks whether a
            ! carrier has taken more nuclei than its element holds, which is
            ! a question about the partition of an element among its
            ! species; an element unknown IS the element and has no such
            ! ceiling above it.  Its own bound, a mass fraction in [0,1], is
            ! the positivity screen above and the projection below.
            !
            ! IT IS ALSO A FACE OF THE UNKNOWN BOX now
            ! (freeze_species_unknown_box), so a projected trial cannot
            ! leave it and this screen cannot refuse one.  It remains the
            ! backstop of the paths no projection covers: the arm that
            ! shortens the step instead (EXHALE_SPECIES_BOUND_PROJECT=0),
            ! the ray probes of the model test, and the marching route's own
            ! evaluations.
            if (srow_kind(ii4) .eq. srow_carrier) then
               ! A carrier's own ceiling is one corner of the shared
               ! constraint and no longer a screen of its own: the shared
               ! magnitude below is the verdict.
               continue
            else if (nspec_j .gt. 1.0d0) then
               headroom_ok = .false.
               if (eval_refusal .eq. refuse_none) then
                  eval_refusal       = refuse_element_fraction_above_1
                  eval_refusal_cell  = jj4
                  eval_refusal_row   = ii4
                  eval_refusal_value = nspec_j
               endif
            endif
         enddo
         enddo
         call write_species_rows_into_composition(Y, u, f_sp)
         ! THE SHARED ELEMENT CONSTRAINT AT THE CANDIDATE, BY MAGNITUDE
         ! (item N4b deliverable 3). What decides feasibility is how far a
         ! candidate stands outside the budget of a cell, not how many cells
         ! stand outside it: one cell outside by a factor of two and twenty
         ! cells outside by a part in 1e12 are not the same state, and a
         ! count cannot tell them apart. The magnitude is the verdict and
         ! the count is kept beside it as the diagnostic it is.
         !
         ! A COMPARISON, NOT A THRESHOLD, for the reason section 138 gives
         ! for the chemistry gate: the budget is frozen at the iterate and
         ! depends on the very unknowns it constrains, so the iterate itself
         ! can sit outside it, and a rule that rejects the neighborhood of
         ! the point it stands on cannot be used to leave that point. The
         ! iterate's own magnitude is therefore the reference, and the slack
         ! is the representability bound of the map -- eight units in the
         ! last place of two summation orders of the same nuclei -- and not
         ! a tolerance chosen to let a state through.  A negative density
         ! has no such excuse and is refused outright.
         !
         ! It is evaluated on the composition the write-back WROTE, so it
         ! carries the restored element totals and the frozen stages, which
         ! is the state the sweep below is about to be asked about.
         call element_budget_violation_of_composition(u, f_sp,             &
                              v_headroom_last, jheadroom, kheadroom,      &
                              n_headroom_last)
         if (v_headroom_last .gt. v_headroom_iterate                       &
                                  + inventory_feasible_slack) then
            headroom_ok = .false.
            if (eval_refusal .eq. refuse_none) then
               eval_refusal       = refuse_element_budget_breach
               eval_refusal_cell  = jheadroom
               eval_refusal_value = v_headroom_last
            endif
         endif
         ! AND THE WRITE-BACK'S OWN STATEMENT: a cell whose carriers took
         ! nuclei that no species of the cell gave up is not a state,
         ! whatever its budgets say, and it is refused here rather than
         ! evaluated.
         if (n_write_back_element_deficit .gt. n_write_back_deficit_entry) &
         then
            headroom_ok = .false.
            if (eval_refusal .eq. refuse_none) then
               eval_refusal       = refuse_element_budget_breach
               eval_refusal_value = write_back_worst_element_deficit
            endif
         endif
      endif
      ! Composition of the INTERIOR of Y, evaluated before the ghosts are
      ! filled: Apply_BC reads n_part_cell1 for the continuous-temperature base
      ! ghost (T(1) = p(1)/n_part_cell1) and that global is written by
      ! get_species_densities. Without this call it still holds the value left
      ! by the PREVIOUS residual evaluation, so F would depend on the previous Y
      ! as well as on Y -- and a finite-difference Jacobian column would then
      ! mix two states. Cost is one array pass; ioniz_eq below dominates.
      ! u(1,:) is already the density, so no U_to_W here: the ghosts are still
      ! zero at this point and U_to_W would divide by them.
      rho = u(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call Apply_BC(u)           ! fill ghosts from the interior

      ! THE COMPOSITION IS AN ELIMINATED VARIABLE AND IT IS ELIMINATED AT
      ! THIS STATE: the residual is R(Y) = L(Y) + S(Y, c*(Y)) with c*(Y) the
      ! equilibrium composition OF Y, so the sweep below is repeated at fixed
      ! Y until the composition stops moving.
      !
      ! ONE SWEEP IS NOT c*(Y). The optical depths, the excited-hydrogen
      ! populations and the electron density that set the rate coefficients
      ! are formed from the composition the sweep was SEEDED with, so a
      ! single sweep is one Picard step of that nonlocal coupling and its
      ! answer still carries the seed. Evaluating S with the seed's
      ! composition makes the Newton drive a LAGGED system to zero, and at
      ! the fixed point of the lagged system the composition sweep still
      ! moves the state: that fixed point is not a root of the coupled
      ! hydrodynamic and composition system.
      !
      ! MEASURED on `wasp_full_newton` (section 155 and the ladder recorded
      ! with it): the energy row of the state handed back, re-measured after
      ! 0, 1, 2, 3 and 8 passes at fixed Y, was 5.466e-06, 1.927e-04,
      ! 3.219e-04, 3.850e-04 and 4.221e-04 -- a factor 77 between the number
      ! the loop-top gate read and the number that state carries.
      !
      ! THE STOPPING TEST IS ON THE COMPOSITION, the quantity being
      ! converged, and not on the norm of the assembled residual: it is
      ! available without the flux assembly, so a pass costs a sweep and
      ! nothing else.
      !
      ! A FIXED COUNT, WHEN THE CALLER ASKS. A finite-difference probe has to
      ! sample ONE map: F(Y + eps v) and F(Y) must be the same number of
      ! passes from the same seed, or their difference divided by eps is
      ! dominated by the difference between two maps rather than by a
      ! directional derivative. n_eq_sweeps_fixed runs exactly that many
      ! passes with no test; an early-exit loop that stopped after k passes
      ! and a fixed count of k are the same map, the test being taken after
      ! the pass.
      n_sweep_max = 1
      if (resid_at_own_composition)                                      &
         n_sweep_max = max(1, n_selfconsistent_max)
      if (n_eq_sweeps_cap .ge. 0) n_sweep_max = max(1, n_eq_sweeps_cap)
      eq_fixed_count = .false.
      if (present(n_eq_sweeps_fixed)) then
         n_sweep_max    = max(1, n_eq_sweeps_fixed)
         eq_fixed_count = .true.
      endif
      n_eq_sweeps_last = 0
      eq_sweep_increment_last = 0.0d0
      radiative_field_of_previous_pass = .false.
      npart_pass = 0.0d0;  xh2_pass = 0.0d0
      heat_pass  = 0.0d0;  cool_pass = 0.0d0
      do it_eq = 1, n_sweep_max
      f_sp_prev = f_sp
      call U_to_W(u, W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      ! The three channels as this pass enters. The radiative pair is the one
      ! the PREVIOUS pass left in the dummy arrays, which is why the first
      ! pass has none; the two equation-of-state channels are formed from the
      ! composition the pass is entered with, by the same expressions
      ! caloric_state_from_composition uses.
      npart_pass(1:N) = n_tot(1:N) + ne(1:N)
      xh2_pass = 0.0d0
      if (thereis_mol) then
         where (npart_pass .gt. 0.0d0)                                    &
            xh2_pass = rho(1:N)*f_sp(1:N,isp_H2)/npart_pass
      endif
      if (radiative_field_of_previous_pass) then
         heat_pass(1:N) = heat(1:N)
         cool_pass(1:N) = cool(1:N)
      endif
      call comp_T_from_p(p, n_tot, ne, T)
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,rel_change)
      call ioniz_eq(T,rho,f_sp,heat,cool,eta,sweep)
      ! Refresh the particle count from the equilibrium fractions, exactly as
      ! the marching loop does before its transport stage, so the residual's
      ! T = p/(n_tot+n_e) is the same temperature the marching step diffuses.
      ! Only the transport terms read it; with them off nothing changes.
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      ! THE GHOSTS ARE WRITTEN AGAIN, WITH THE COMPOSITION THE FLUXES WILL
      ! USE. The Apply_BC above ran before the sweep, because U_to_W needs a
      ! ghost density to divide by and ioniz_eq needs a ghost pressure; the
      ! conservative ghost state it left therefore carried the composition of
      ! the PREVIOUS evaluation, while Reconstruct and RK_rhs below read it
      ! alongside interior cells carrying this one. Under a constant adiabatic
      ! index the primitive-to-conservative round trip returns the pressure the
      ! boundary condition asked for and the two cannot disagree; under the
      ! caloric EOS of section 141 they can, by as much as the ghost
      ! composition moved in one sweep -- largest during relaxation, which is
      ! when the solver is trying to make progress (section 141.6, item (AB)).
      !
      ! This is now the order the MARCHING LOOP has always had, and that is the
      ! point: there Apply_BC follows the composition refresh, so the ghosts a
      ! step's fluxes read were written with that step's composition. The
      ! residual and the marching right-hand side have to be the same discrete
      ! operator or they cannot have the same fixed point (section 144.1).
      !
      ! It costs nothing but the call: Apply_BC writes the ghosts and returns
      ! the interior unchanged, bit for bit, which is what made a second call
      ! possible at all.
      !
      ! WHAT IS STILL ONE SWEEP BEHIND, and deliberately. The ghost cells' own
      ! composition was solved by the sweep above at the ghost pressure the
      ! FIRST Apply_BC wrote, so the EOS used here is that composition and not
      ! one re-solved at the pressure just written. Iterating the pair to a
      ! ghost fixed point would make this residual something the marching loop
      ! is not, and the requirement above is the stronger one.
      call Apply_BC(u)
      n_eq_sweeps_last = it_eq
      if (.not. eq_fixed_count) then
         if (eq_measure_is_species_fraction) then
            df_eq = maxval(abs(f_sp(1:N,:) - f_sp_prev(1:N,:))            &
                           /max(abs(f_sp(1:N,:)), eq_sweep_floor))
            jworst_eq = 0;  kworst_eq = 0
         else
            npart_now(1:N) = n_tot(1:N) + ne(1:N)
            xh2_now = 0.0d0
            if (thereis_mol) then
               where (npart_now .gt. 0.0d0)                               &
                  xh2_now = rho(1:N)*f_sp(1:N,isp_H2)/npart_now
            endif
            call increment_of_what_the_residual_reads(npart_now,          &
                    npart_pass, xh2_now, xh2_pass, heat(1:N), cool(1:N),  &
                    heat_pass, cool_pass, radiative_field_of_previous_pass,&
                    df_eq, jworst_eq, kworst_eq)
         endif
         if (eq_measure_trace) then
            eq_channel_txt = 'species fraction'
            if (.not. eq_measure_is_species_fraction)                     &
               call eq_channel_name(kworst_eq, eq_channel_txt)
            write(*,'(A,I4,A,ES10.3,A,I5,A,A)') ' (eq_measure) pass',     &
                 it_eq, '  increment ', df_eq, '  cell ', jworst_eq,      &
                 '  channel ', trim(eq_channel_txt)
         endif
         eq_sweep_increment_last = df_eq
         if (df_eq .le. eq_sweep_reltol) exit
      endif
      radiative_field_of_previous_pass = .true.
      enddo
      n_eq_sweeps_tot = n_eq_sweeps_tot + n_eq_sweeps_last

      if (resid_capture_operator_state) then
         ! The pipeline's own inputs, kept for the jump scan before the
         ! flux assembly consumes them.
         if (.not. allocated(captured_state_with_ghosts))                &
            allocate(captured_state_with_ghosts(3,1-Ng:N+Ng),            &
                     captured_temperature(1-Ng:N+Ng),                    &
                     captured_particle_count(1-Ng:N+Ng))
         captured_state_with_ghosts = u
         captured_temperature       = T
         captured_particle_count    = n_tot + ne
      endif
      call assemble_residual(u, n_tot + ne, heat, cool, R)
      call pack_R(R, Fvec)
      ! The carrier row of the same state: the steady continuity equation of
      ! n(H2), assembled by the module that owns it, on the background this
      ! evaluation's own sweep just left. Its natural units are cm^-3 s^-1
      ! while the hydrodynamic rows are per code time, so it is converted by
      ! the code's own time scale -- otherwise the merit would be comparing
      ! two different clocks.
      if (nspec_row .gt. 0) then
         crc_max = 0.0d0;  crc_vol = 0.0d0;  crc_j = 0;  crc_ic = 1
         cres = 0.0d0
         if (.not. allocated(crow_res_last))                             &
            allocate(crow_res_last(1:N,n_carrier_max),                   &
                     crow_terms_last(1:N,n_carrier_max))
         crow_res_last   = 0.0d0
         crow_terms_last = 0.0d0
         if (thereis_mol .and. carrier_transport)                        &
            call carrier_steady_residual(rho, v, f_sp, crc_max, crc_j,   &
                                         crc_ic, rvol=crc_vol,           &
                                         res_out=cres,                   &
                                         terms_out=crow_terms_last)
         crow_res_last = cres
         if (element_row_terms_diag)                                     &
            call write_carrier_row_terms(f_sp, trace_row_terms_from, N)
         ! THE ELEMENTAL ROWS OF THE SAME STATE, assembled by the operator
         ! that owns them, on this evaluation's own composition. Their
         ! natural units are a mass rate and a mixing ratio per second while
         ! the hydrodynamic rows are per code time, so each is converted by
         ! its own factor -- otherwise the merit would compare two clocks.
         if (element_rows_registered()) then
            call comp_T_from_p(p, n_tot, ne, T)
            ! The wind those rows ride on: the face mass flux of THIS state,
            ! the one the mass row assembled just above differences. Read
            ! here and handed over, because the element operator stands
            ! below this solver and cannot reach it.
            call face_mass_flux_of_state(rho, Frho_elem)
            call element_transport_residual(rho, T, f_sp, Frho_elem,      &
                     erow_he, escale_he, ehe_ok, erow_tr, escale_tr,      &
                     etr_ok, tr_carried = erow_tr_carried)
            if (.not. ehe_ok) then
               erow_he = 0.0d0;  escale_he = 1.0d0
            endif
            if (.not. etr_ok) then
               erow_tr = 0.0d0;  escale_tr = 1.0d0
            endif
            ! The terms of the HELIUM row at the cells a caller asked for,
            ! on this evaluation's own composition and on the wind the rows
            ! above were formed with. Armed for one evaluation of the
            ! iterate and disarmed by the caller.
            if (element_row_terms_diag .and. ehe_ok)                      &
               call write_element_row_terms(rho, T, f_sp, Frho_elem,      &
                        erow_he, escale_he, trace_row_terms_from, N)
         endif
         tscale_code = R0/(v0*n0)
         do jj4 = 1, N
         do ii4 = 1, nspec_row
            select case (srow_kind(ii4))
            case (srow_carrier)
               Fvec(nvar_jac*(jj4-1)+3+ii4) =                            &
                  cres(jj4,srow_idx(ii4))*tscale_code
            case (srow_element_he)
               if (jj4 .eq. 1 .and. species_row_base_equation(ii4) .eq.  &
                                    base_row_reservoir_condition) then
                  ! THE BASE CELL IS THE DIRICHLET RESERVOIR of the element
                  ! operator, so its equation is that statement and not the
                  ! transport balance, which the operator does not pose
                  ! there (species_row_base_equation, srow_base_value).
                  Fvec(nvar_jac*(jj4-1)+3+ii4) =                         &
                     Y(nvar_jac*(jj4-1)+3+ii4) - srow_base_value(ii4)
               else
                  Fvec(nvar_jac*(jj4-1)+3+ii4) =                         &
                     erow_he(jj4)*elem_he_to_code()
               endif
            case default
               if (.not. erow_tr_carried(srow_idx(ii4))) then
                  ! No reservoir, no equation: the row says the element
                  ! stays absent.  A zero row would leave the column of
                  ! this unknown empty and the band factorization singular.
                  Fvec(nvar_jac*(jj4-1)+3+ii4) = Y(nvar_jac*(jj4-1)+3+ii4)
               else if (jj4 .eq. 1 .and.                                 &
                        species_row_base_equation(ii4) .eq.              &
                        base_row_reservoir_condition) then
                  Fvec(nvar_jac*(jj4-1)+3+ii4) =                         &
                     Y(nvar_jac*(jj4-1)+3+ii4) - srow_base_value(ii4)
               else
                  Fvec(nvar_jac*(jj4-1)+3+ii4) =                         &
                     erow_tr(jj4,srow_idx(ii4))*elem_tr_to_code()
               endif
            end select
         enddo
         enddo
         carrier_relnorm_last = crc_vol
         carrier_cellmax_last = crc_max
         carrier_cell_worst   = crc_j
      endif
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      if (present(n_part)) n_part = n_tot + ne
      ! Recorded on EVERY evaluation, so that a caller which does not ask for
      ! admissibility -- the one that evaluates the current iterate -- can
      ! still read off what that iterate's own chemistry cost.
      n_no_chem_root_last = n_cells_without_chemical_root(sweep%acc_n)

      if (present(admissible)) then
         ! Two independent statements, both required. The sweep's own count
         ! catches a state that broke a balance row of the chemical network,
         ! which is where a probe state fails first and where the diagnosis
         ! lies; the row scan catches anything non-finite that reached the
         ! residual by another route (an overflowed flux, a cooling term).
         gate_chem_root = .true.
         if (present(may_be_adopted)) gate_chem_root = may_be_adopted
         ! The element budget is a CONSTRAINT on the state, so it is tested
         ! here and not imposed by a clamp inside the row. A trial whose
         ! n(H2) exceeds the hydrogen its cell has, or goes negative, is not
         ! a state the write-back could describe; refusing it is the same
         ! treatment the positivity of rho and p already gets in the line
         ! search.
         if (.not. headroom_ok) n_trial_headroom = n_trial_headroom + 1
         ! WHAT THIS EVALUATION FOUND, handed to the context that asked for
         ! it. The decision itself lives in certification.f90 so that the
         ! probe and the trial are two named readings of one record and not
         ! two copies of a predicate. The row scan is taken only where the
         ! rest of the record already admits the state, exactly as before:
         ! it is the expensive half.
         cfacts%n_sweep_nonfinite    = sweep%n_nonfinite
         cfacts%n_no_chem_root_trial = n_no_chem_root_last
         cfacts%n_no_chem_root_state = n_no_chem_root_state
         cfacts%chem_root_gate       = gate_chem_root
         cfacts%headroom_ok          = headroom_ok
         cfacts%rows_finite          = .true.
         if (gate_chem_root .and.                                       &
             n_no_chem_root_last .gt. n_no_chem_root_state) then
            n_trial_no_chem_root  = n_trial_no_chem_root + 1
            no_chem_root_cell_max = max(no_chem_root_cell_max,              &
                                       n_no_chem_root_last)
            ! The worst reaction residual of a cell the chemistry does not
            ! describe, over BOTH non-root classes: 4, the state kept under
            ! the relaxation amnesty, and 6, the candidate above the amnesty
            ! cap for which the cell kept its entry composition. The count
            ! this residual belongs to (n_cells_without_chemical_root) is the
            ! sum of the two, so class 4 alone reported zero for a sweep
            ! whose non-root cells were all class 6.
            no_chem_root_res_max  = max(no_chem_root_res_max,               &
                                        max(sweep%acc_resmax(4),            &
                                            sweep%acc_resmax(6)))
         endif
         if (trial_state_is_admissible(cfacts)) then
            do irow = 1, nvar_jac*N
               if (.not. finite_real(Fvec(irow))) then
                  cfacts%rows_finite = .false.
                  if (eval_refusal .eq. refuse_none) then
                     eval_refusal      = refuse_row_nonfinite
                     eval_refusal_cell = (irow - 1)/nvar_jac + 1
                     eval_refusal_row  = irow - nvar_jac*(eval_refusal_cell-1)
                  endif
                  exit
               endif
            enddo
         endif
         if (cfacts%n_sweep_nonfinite .gt. 0 .and.                       &
             eval_refusal .eq. refuse_none) then
            eval_refusal       = refuse_sweep_nonfinite
            eval_refusal_value = dble(cfacts%n_sweep_nonfinite)
         endif
         if (gate_chem_root) then
            admissible = trial_state_is_admissible(cfacts)
            if (.not. admissible .and. eval_refusal .eq. refuse_none) then
               eval_refusal       = refuse_no_chemical_root
               eval_refusal_value = dble(n_no_chem_root_last              &
                                         - n_no_chem_root_state)
            endif
         else
            admissible = probe_direction_is_usable(cfacts)
         endif
         if (admissible) eval_refusal = refuse_none
         n_eval_refusal(eval_refusal) = n_eval_refusal(eval_refusal) + 1
      endif

      if (present(rows_judged)) call certified_row_measures(Fvec, u,     &
                                        rows_judged, resid_relnorm_judged)
      if (present(rowmax_judged) .or. present(cells_judged) .or.         &
          present(which_judged)) then
         call row_maxima_of_the_certification(Fvec, u, rmx_j, jcl_j, iwh_j)
         if (present(rowmax_judged)) rowmax_judged = rmx_j
         if (present(cells_judged))  cells_judged  = jcl_j
         if (present(which_judged))  which_judged  = iwh_j
      endif

      ! The residual and the measures of THIS state are in the caller's
      ! arguments; the modules go back to the state this evaluation found
      ! them in.
      if (put_products_back)                                             &
         call put_back_residual_evaluation_products(products_at_entry)

      end subroutine eval_residual

      ! ------------------------------------------------------!

      subroutine newton_residual(Y, f_sp_seed, f_sp, Fvec)
      ! Thin wrapper: full residual, discarding the heat/cool diagnostics.
      ! Seed in, composition out, for the reason eval_residual gives.
      real*8, dimension(nvar_jac*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)    :: f_sp_seed
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(out)   :: f_sp
      real*8, dimension(nvar_jac*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      call eval_residual(Y, f_sp_seed, f_sp, Fvec, heat, cool)
      end subroutine newton_residual

      ! ------------------------------------------------------!

      subroutine frozen_residual(Y, n_part, heat, cool, Fvec)
      ! FROZEN-radiation residual: the hydro flux/gravity residual at Y with
      ! heat/cool held FIXED (no ioniz_eq, no column-density recompute). This
      ! is strictly local (WENO3 stencil) and hence exactly banded, so the
      ! colored finite-difference Jacobian of THIS residual is exact. It is
      ! the approximate Jacobian / preconditioner for the inexact Newton: the
      ! weakly non-local radiation response is left to the outer iteration.
      real*8, dimension(nvar_jac*N),       intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: n_part, heat, cool
      real*8, dimension(nvar_jac*N),       intent(out) :: Fvec
      real*8, dimension(3,1-Ng:N+Ng) :: u, R
      integer :: j, i
      u = 0.0d0
      call unpack_U(Y, u)
      call Apply_BC(u)
      call assemble_residual(u, n_part, heat, cool, R)
      call pack_R(R, Fvec)
      ! THE SPECIES SLOTS CARRY NO EQUATION IN THIS RESIDUAL, and they are
      ! written rather than left as they were found.  pack_R fills three
      ! rows of each cell; with a species row registered the vector has
      ! more, and an output argument that is only partly written is read by
      ! the caller as whatever its memory held -- two calls then differ in
      ! those slots by an amount that is not a property of Y, and a finite
      ! difference of the two divides that by eps.  MEASURED before this
      ! line: the directional check of a two-species system reported
      ! ||J r - dFD||/||J r|| = 1.0, entirely from those slots.
      !
      ! Zero is the right value and not a placeholder: this residual holds
      ! the radiation frozen and re-evaluates no chemistry, so it states no
      ! species balance, and the banded Jacobian built from it is a
      ! hydrodynamic preconditioner whose species block is empty.  The
      ! Jacobian the steady solve actually uses is build_banded_jac_full,
      ! on eval_residual, which does carry the species rows.
      do j = 1, N
         do i = 1, nspec_row
            Fvec(nvar_jac*(j-1)+3+i) = 0.0d0
         enddo
      enddo
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      end subroutine frozen_residual

      ! ------------------------------------------------------!

      subroutine pack_R(R, Fvec)
      ! Pack residual array (physical cells) into the flat vector.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: R
      real*8, dimension(nvar_jac*N),         intent(out) :: Fvec
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Fvec(nvar_jac*(j-1)+k) = R(k,j)
         enddo
      enddo
      end subroutine pack_R

      ! ------------------------------------------------------!

      subroutine build_banded_jac(Y, n_part, heat, cool, ab)
      ! Colored finite-difference banded Jacobian of the FROZEN residual,
      ! stored in LAPACK general-band form for dgbtrf/dgbtrs:
      !   ab(kl+ku+1 + i - j, j) = J(i,j),  i in [j-ku, j+kl]
      ! ldab = 2*kl+ku+1. Columns sharing a color (stride ncolor = kl+ku+1)
      ! have disjoint row supports, so one probe per color recovers their
      ! band entries with no cross-contamination.
      real*8, dimension(nvar_jac*N),               intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng),         intent(in)  :: n_part, heat, cool
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(out) :: ab
      real*8, dimension(nvar_jac*N) :: F0, Fp, Yp, dYc
      integer :: neq, color, jcol, irow, ilo, ihi
      real*8  :: sqeps

      neq   = nvar_jac*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0

      call frozen_residual(Y, n_part, heat, cool, F0)

      do color = 1, ncolor_jac
         Yp  = Y
         dYc = 0.0d0
         do jcol = color, neq, ncolor_jac
            dYc(jcol) = sqeps*max(abs(Y(jcol)), 1.0d0)
            Yp(jcol)  = Y(jcol) + dYc(jcol)
         enddo
         call frozen_residual(Yp, n_part, heat, cool, Fp)
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac

      ! ------------------------------------------------------!

      subroutine build_banded_jac_full(Y, f_sp_base, ab, n_color_unresolved)
      ! Colored-FD banded Jacobian of the FULL residual (eval_residual, i.e.
      ! including the local ionization-equilibrium + heat/cool response). The
      ! frozen version omits the dominant local d(heat-cool)/dE coupling, so
      ! Newton has no descent direction even near the solution; this version
      ! captures it. The weakly non-local column-density coupling slightly
      ! contaminates same-color columns (treated as preconditioner error).
      ! Each probe restores f_sp from the base copy (ioniz_eq mutates it).
      real*8, dimension(nvar_jac*N),                  intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(out) :: ab
      ! Colors whose probe stayed inadmissible down to the smallest step, so
      ! their columns are left at zero (see n_probe_step_halvings above).
      integer, intent(out) :: n_color_unresolved
      real*8, dimension(nvar_jac*N) :: F0, Fp, Yp, dYc, Dsc
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      integer :: neq, color, jcol, irow, ilo, ihi, ish
      real*8  :: sqeps, hstep
      logical :: probe_ok

      neq   = nvar_jac*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0
      n_color_unresolved = 0

      ! FD column steps RELATIVE to the scale for each unknown (a flat
      ! max(|Y|,1) floor gives wind cells, ~1e-7 of the base in code
      ! units, order-unity relative kicks and garbage columns).
      call cell_state_scales(Y, Dsc)

      ! THE BASE POINT AND EVERY COLUMN ARE THE SAME MAP: a fixed number of
      ! composition passes, not a converged one. A converged elimination
      ! would take a different number of passes at Y and at Y + h, and the
      ! difference of two maps divided by h is not a derivative.
      call eval_residual(Y, f_sp_base, fwork, F0, heat, cool,             &
                         admissible=probe_ok,                             &
                         may_be_adopted=.false.,                          &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      if (.not. probe_ok) then
         ! The base point of the differences is itself indescribable. No
         ! column of this Jacobian means anything, so none is built and the
         ! caller sees every color unresolved.
         n_color_unresolved = ncolor_jac
         return
      endif

      do color = 1, ncolor_jac
         hstep    = sqeps
         probe_ok = .false.
         do ish = 0, n_probe_step_halvings
            Yp  = Y;  dYc = 0.0d0
            do jcol = color, neq, ncolor_jac
               dYc(jcol) = hstep*Dsc(jcol)
               Yp(jcol)  = Y(jcol) + dYc(jcol)
            enddo
            call eval_residual(Yp, f_sp_base, fwork, Fp, heat, cool,      &
                               admissible=probe_ok,                       &
                               may_be_adopted=.false.,                    &
                               state_is_discarded=.true.,                 &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            if (probe_ok) exit
            hstep = 0.5d0*hstep
         enddo
         if (.not. probe_ok) then
            n_color_unresolved = n_color_unresolved + 1
            cycle
         endif
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac_full

      ! ------------------------------------------------------!

      subroutine band_matvec(ab, x, y)
      ! y = J*x for the band-stored Jacobian ab (same layout as above).
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8, dimension(nvar_jac*N),                    intent(in)  :: x
      real*8, dimension(nvar_jac*N),                    intent(out) :: y
      integer :: neq, irow, jcol, jlo, jhi
      neq = nvar_jac*N
      do irow = 1, neq
         y(irow) = 0.0d0
         jlo = max(1,   irow - kl_jac)
         jhi = min(neq, irow + ku_jac)
         do jcol = jlo, jhi
            y(irow) = y(irow) + ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*x(jcol)
         enddo
      enddo
      end subroutine band_matvec

      ! ------------------------------------------------------!

      subroutine band_matvec_transpose(ab, x, y)
      ! y = J^T*x for the band-stored Jacobian ab (same layout as above).
      ! The transpose is needed because the gradient of the least-squares
      ! merit 1/2*||F||^2 is J^T F, and the Newton/PTC step does not use it.
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8, dimension(nvar_jac*N),                    intent(in)  :: x
      real*8, dimension(nvar_jac*N),                    intent(out) :: y
      integer :: neq, irow, jcol, jlo, jhi
      neq = nvar_jac*N
      y = 0.0d0
      do irow = 1, neq
         jlo = max(1,   irow - kl_jac)
         jhi = min(neq, irow + ku_jac)
         do jcol = jlo, jhi
            y(jcol) = y(jcol) + ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*x(irow)
         enddo
      enddo
      end subroutine band_matvec_transpose

      ! ------------------------------------------------------!

      subroutine name_of_unknown(i, txt)
      ! The cell and the physical quantity the flat index i belongs to, for
      ! the diagnostics below: an index alone cannot say which equation is
      ! being reported.
      integer,          intent(in)  :: i
      character(len=*), intent(out) :: txt
      integer :: j, k, is
      character(len=24) :: what
      j = (i - 1)/nvar_jac + 1
      k = i - nvar_jac*(j - 1)
      select case (k)
      case (1);  what = 'mass'
      case (2);  what = 'momentum'
      case (3);  what = 'energy'
      case default
         is = k - 3
         if (is .lt. 1 .or. is .gt. nspec_row) then
            what = 'unregistered species'
         else if (srow_kind(is) .eq. srow_carrier) then
            what = 'carrier '//trim(carrier_name(srow_idx(is)))
         else if (srow_kind(is) .eq. srow_element_he) then
            what = 'element He'
         else
            what = 'element '//trim(melem_name(srow_idx(is)))
         endif
      end select
      write(txt,'(A,A,I0)') trim(what), ' of cell ', j
      end subroutine name_of_unknown

      ! ------------------------------------------------------!

      subroutine tr_exit_text(code, txt)
      ! The sentence that belongs to one trust-region exit code.
      integer,          intent(in)  :: code
      character(len=*), intent(out) :: txt
      select case (code)
      case (tr_step_taken)
         txt = 'a step was taken'
      case (tr_no_krylov_direction)
         txt = 'no Krylov direction could be sampled'
      case (tr_cauchy_probe_inadmissible)
         txt = 'the Cauchy probe has no admissible residual'
      case (tr_dogleg_step_is_zero)
         txt = 'the dogleg step is zero'
      case (tr_no_admissible_state_on_ray)
         txt = 'no admissible state anywhere on the dogleg'
      case (tr_model_promises_no_drop)
         txt = 'the model promises no reduction'
      case (tr_ray_probe_inadmissible)
         txt = 'a ray probe has no admissible residual'
      case (tr_ray_slopes_disagree_sign)
         txt = 'the model and the true slope differ in sign'
      case (tr_ray_slopes_disagree_size)
         txt = 'the model and the true slope differ in size'
      case (tr_trial_no_longer_admissible)
         txt = 'the accepted trial is no longer admissible'
      case (tr_ratio_below_acceptance)
         txt = 'the reduction ratio is below eta'
      case (tr_ray_below_the_noise_floor)
         txt = 'the ray difference is below its cancellation floor'
      case (tr_projected_step_has_no_model)
         txt = 'the step projected onto the species bounds has no'//      &
               ' describable sample'
      case (tr_no_direction_reduces_model)
         txt = 'neither leg is known to reduce the model'
      case default
         txt = 'unnamed'
      end select
      end subroutine tr_exit_text

      ! ------------------------------------------------------!

      subroutine write_last_evaluation_refusal(tag)
      ! WHICH SCREEN REFUSED THE LAST RESIDUAL SAMPLE, at which cell and
      ! which unknown. Printed where a step control reports that it had no
      ! sample to work with, so that "no Krylov direction could be sampled"
      ! carries its cause. Silent when the last evaluation was admissible.
      character(len=*), intent(in) :: tag
      character(len=48) :: txt
      integer :: i
      if (eval_refusal .eq. refuse_none) return
      txt = ''
      if (eval_refusal_cell .gt. 0 .and. eval_refusal_row .gt. 0) then
         i = nvar_jac*(eval_refusal_cell-1) + 3 + eval_refusal_row
         if (eval_refusal .eq. refuse_row_nonfinite)                      &
            i = nvar_jac*(eval_refusal_cell-1) + eval_refusal_row
         call name_of_unknown(i, txt)
      endif
      select case (eval_refusal)
      case (refuse_species_unknown_negative)
         write(*,'(A,A,A,A,ES10.3,A,I0,A)') tag,                          &
              '      refused: ', trim(txt), ' is negative,',              &
              eval_refusal_value, '; ', n_species_at_lower_bound,         &
              ' species unknown(s) of this state sit on their lower bound'
      case (refuse_element_fraction_above_1)
         write(*,'(A,A,A,A,ES10.3)') tag, '      refused: ', trim(txt),   &
              ' is a mass fraction above one,', eval_refusal_value
      case (refuse_element_budget_breach)
         write(*,'(A,A,ES10.3,A,I0,A,I0,A)') tag, '      refused: the'//  &
              ' shared element budget is breached by ',                   &
              eval_refusal_value, ' of what the row allows, at cell ',    &
              eval_refusal_cell, ' (', n_headroom_last,                   &
              ' cell-rows breaching)'
      case (refuse_sweep_nonfinite)
         write(*,'(A,A,I0,A)') tag, '      refused: ',                    &
              int(eval_refusal_value), ' cell(s) whose reaction'//        &
              ' residual is not a finite number'
      case (refuse_row_nonfinite)
         write(*,'(A,A,A,A)') tag, '      refused: the residual row of ', &
              trim(txt), ' is not a finite number'
      case (refuse_no_chemical_root)
         write(*,'(A,A,I0,A)') tag, '      refused: ',                    &
              int(eval_refusal_value), ' cell(s) more than the iterate'// &
              ' carry a composition that is not a root of the network'
      end select
      end subroutine write_last_evaluation_refusal

      ! ------------------------------------------------------!

      subroutine banded_model_conditioning(ab, abf, D, Drow, tag)
      ! WHAT THE BANDED MODEL STATES, AND WHAT IT LEAVES UNSTATED.
      !
      ! The band is the local part of the operator and the preconditioner of
      ! the Krylov solve, and the damped Gauss-Newton escape factors J^T J of
      ! it. All three rest on the band being a NON-SINGULAR statement about
      ! every unknown. A row with no entry inside the band says nothing about
      ! any unknown, and a column with none says that no equation of the model
      ! responds to that unknown; either makes the banded system singular, so
      ! its LU leaves a pivot at round-off and the preconditioner it defines
      ! amplifies that direction without bound. This reports both, and the
      ! smallest pivot of the factorization beside them, with the cell and the
      ! quantity named.
      !
      ! THE EMPTY-ROW AND EMPTY-COLUMN STATEMENT IS UNCONDITIONAL, and it
      ! prints only when it fires: it is an assertion that the system posed
      ! is not singular by construction, and it costs one pass over the band
      ! against the 6*nvar - 1 residual evaluations the band itself cost.
      ! The pivot and largest-entry lines, which print on every iteration,
      ! are behind EXHALE_BAND_REPORT=1.
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab, abf
      ! The two scalings the band was formed under, so that the largest
      ! entry can be decomposed into the derivative and the scales: an entry
      ! of 1e22 is either a real sensitivity or a row scale that has fallen
      ! below the size on which its own row responds, and only the
      ! decomposition says which.
      real*8, dimension(nvar_jac*N), intent(in) :: D, Drow
      character(len=*), intent(in) :: tag
      integer :: neq, irow, jcol, n_row0, n_col0, i_row0, i_col0, i_piv
      integer :: i_amax, j_amax
      real*8  :: rowmax, colmax, upiv, umin, umax, amax
      character(len=48) :: txt, txt2
      character(len=8)  :: env
      logical :: verbose
      call get_environment_variable('EXHALE_BAND_REPORT', env)
      verbose = (trim(env) .eq. '1')
      neq = nvar_jac*N
      n_row0 = 0;  n_col0 = 0;  i_row0 = 0;  i_col0 = 0
      amax = 0.0d0;  i_amax = 1;  j_amax = 1
      do irow = 1, neq
         rowmax = 0.0d0
         do jcol = max(1, irow-kl_jac), min(neq, irow+ku_jac)
            rowmax = max(rowmax, abs(ab(kl_jac+ku_jac+1+irow-jcol, jcol)))
         enddo
         if (rowmax .gt. amax) then
            amax = rowmax
            i_amax = irow
            do jcol = max(1, irow-kl_jac), min(neq, irow+ku_jac)
               if (abs(ab(kl_jac+ku_jac+1+irow-jcol, jcol)) .eq. rowmax)  &
                  j_amax = jcol
            enddo
         endif
         if (rowmax .le. 0.0d0) then
            n_row0 = n_row0 + 1
            if (i_row0 .eq. 0) i_row0 = irow
         endif
      enddo
      do jcol = 1, neq
         colmax = 0.0d0
         do irow = max(1, jcol-ku_jac), min(neq, jcol+kl_jac)
            colmax = max(colmax, abs(ab(kl_jac+ku_jac+1+irow-jcol, jcol)))
         enddo
         if (colmax .le. 0.0d0) then
            n_col0 = n_col0 + 1
            if (i_col0 .eq. 0) i_col0 = jcol
         endif
      enddo
      ! The assertion, stated once per solve: an equation that responds to no
      ! unknown, or an unknown no equation responds to. It is a property of
      ! the system posed and not of the iterate, so repeating it every outer
      ! iteration would say nothing further.
      if (band_emptiness_stated) then
         if (.not. verbose) return
      endif
      band_emptiness_stated = .true.
      if (i_row0 .gt. 0) then
         call name_of_unknown(i_row0, txt)
         write(*,'(A,A,I0,A,A)') tag, ' the banded model has a'//         &
              ' structurally empty row: ', n_row0, ' in all, the first'// &
              ' the equation of ', trim(txt)
      endif
      if (i_col0 .gt. 0) then
         call name_of_unknown(i_col0, txt)
         write(*,'(A,A,I0,A,A)') tag, ' the banded model has a'//         &
              ' structurally empty column: ', n_col0, ' in all, the'//    &
              ' first the unknown ', trim(txt)
      endif
      if (.not. verbose) return
      ! U(j,j) of the banded LU sits on the same diagonal of the factored
      ! array as the matrix diagonal did.
      umin = huge(1.0d0);  umax = 0.0d0;  i_piv = 1
      do jcol = 1, neq
         upiv = abs(abf(kl_jac+ku_jac+1, jcol))
         if (upiv .lt. umin) then
            umin = upiv;  i_piv = jcol
         endif
         umax = max(umax, upiv)
      enddo
      write(*,'(A,A,I0,A,I0)') tag, ' [band] structurally empty rows ',   &
           n_row0, ', empty columns ', n_col0
      call name_of_unknown(i_piv, txt)
      write(*,'(A,A,ES11.3,A,ES11.3,A,A)') tag, ' [band] |U(j,j)| min ',  &
           umin, ' max ', umax, ' ; smallest at ', trim(txt)
      call name_of_unknown(i_amax, txt)
      call name_of_unknown(j_amax, txt2)
      write(*,'(A,A,ES11.3,A,A,A,A)') tag, ' [band] largest entry ',      &
           amax, ' in the row of ', trim(txt), ' , column ', trim(txt2)
      write(*,'(A,A,ES11.3,A,ES11.3,A,ES11.3)') tag,                      &
           ' [band]   = dF/dY ', amax*Drow(i_amax)/max(D(j_amax),1.0d-300),&
           '  times the column scale ', D(j_amax),                        &
           '  over the row scale ', Drow(i_amax)
      ! The row scales of the species rows, over the column: a row whose
      ! scale has collapsed onto its floor is what makes an entry like the
      ! one above, and the spread says whether that has happened.
      if (nspec_row .gt. 0) then
         do irow = 1, nspec_row
            umin = huge(1.0d0);  umax = 0.0d0;  i_piv = 3 + irow
            do jcol = 1, N
               upiv = Drow(nvar_jac*(jcol-1) + 3 + irow)
               if (upiv .lt. umin) then
                  umin = upiv;  i_piv = nvar_jac*(jcol-1) + 3 + irow
               endif
               umax = max(umax, upiv)
            enddo
            call name_of_unknown(i_piv, txt)
            write(*,'(A,A,A,A,ES11.3,A,ES11.3,A,A)') tag,                 &
                 ' [band] row scale of the ', trim(txt(1:index(txt,' of')-1)),&
                 ' rows: min ', umin, ' max ', umax, ' ; smallest at ',   &
                 trim(txt)
         enddo
      endif
      end subroutine banded_model_conditioning

      ! ------------------------------------------------------!

      subroutine write_scales_columns_and_pivots(D, Drow, ab, abf,        &
                                                 ipiv, idtau, lu_info,    &
                                                 n_unresolved, tag)
      ! WHAT THE TWO SCALINGS AND THE BAND SAY ABOUT THIS ITERATE, in the
      ! order a stalled step has to be read in.
      !
      ! The band is the preconditioner of the Krylov cycle and the model the
      ! predicted reduction is taken in, so three properties of it decide
      ! whether a step exists at all: whether every unknown is sampled,
      ! whether the factorization is a map, and how far the two scalings
      ! spread. A column the coloring never resolved is left at exactly zero
      ! (build_banded_jac_full), which is the same thing as an unknown no
      ! equation responds to; a column whose largest entry stands at the
      ! round-off of the largest column in the band is sampled in name only,
      ! because the preconditioner divides by it.
      !
      ! The condition estimate is dgbcon's reciprocal 1-norm estimate of the
      ! SHIFTED band, the matrix abf holds the factors of, so the 1-norm it
      ! needs is formed here from ab and the pseudo-transient shift. The
      ! pivot ratio is printed beside it: the two answer different
      ! questions, the ratio being a property of one elimination and the
      ! estimate a bound on the solution it defines.
      !
      ! Diagnostic only (elem_diag_on). It costs one pass over the band and
      ! one Hager estimate, and no residual evaluation.
      real*8, dimension(nvar_jac*N), intent(in) :: D, Drow
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab, abf
      integer, dimension(nvar_jac*N), intent(in) :: ipiv
      ! The pseudo-transient inverse time the factored band carries on its
      ! diagonal, in the same form the caller added it.
      real*8,  intent(in) :: idtau
      ! The status dgbtrf returned for abf, and the colors whose probe was
      ! refused at every step length so their columns were left at zero.
      integer, intent(in) :: lu_info, n_unresolved
      character(len=*), intent(in) :: tag
      integer :: neq, ldb, jcol, irow, i, nz, nnf, nsmall, nshown
      integer :: i_dmin, i_dmax, i_rmin, i_rmax, i_piv
      real*8  :: colmax, colbig, anorm, s, aij, rcond, upiv, umin, umax
      real*8  :: dmin, dmax, rmin, rmax
      real*8, allocatable :: work(:), colscan(:), rowscan(:)
      integer, allocatable :: iwork(:)
      logical :: equilibrated
      character(len=48) :: txt
      ! A column is sampled in name only when its largest entry stands this
      ! far below the largest column of the band: the preconditioner solves
      ! with it, so a column at the round-off of the band is a division by
      ! round-off.
      real*8, parameter :: column_round_off = 1.0d-14
      neq = nvar_jac*N
      ldb = 2*kl_jac + ku_jac + 1
      allocate(colscan(neq), rowscan(neq))
      ! --- 1. the extremes of the two scalings, with their quantities ---
      i_dmin = 1;  i_dmax = 1;  i_rmin = 1;  i_rmax = 1
      dmin = huge(1.0d0);  dmax = -1.0d0
      rmin = huge(1.0d0);  rmax = -1.0d0
      do i = 1, neq
         if (D(i) .lt. dmin) then
            dmin = D(i);  i_dmin = i
         endif
         if (D(i) .gt. dmax) then
            dmax = D(i);  i_dmax = i
         endif
         if (Drow(i) .lt. rmin) then
            rmin = Drow(i);  i_rmin = i
         endif
         if (Drow(i) .gt. rmax) then
            rmax = Drow(i);  i_rmax = i
         endif
      enddo
      call name_of_unknown(i_dmin, txt)
      write(*,'(A,A,ES11.3,A,A)') tag, ' [diag 1] column scale min ',     &
           dmin, ' at the ', trim(txt)
      call name_of_unknown(i_dmax, txt)
      write(*,'(A,A,ES11.3,A,A)') tag, ' [diag 1] column scale max ',     &
           dmax, ' at the ', trim(txt)
      call name_of_unknown(i_rmin, txt)
      write(*,'(A,A,ES11.3,A,A)') tag, ' [diag 1] row scale    min ',     &
           rmin, ' at the equation of ', trim(txt)
      call name_of_unknown(i_rmax, txt)
      write(*,'(A,A,ES11.3,A,A)') tag, ' [diag 1] row scale    max ',     &
           rmax, ' at the equation of ', trim(txt)
      ! --- 2. the columns the banded model does not sample ---
      colbig = 0.0d0
      do jcol = 1, neq
         colmax = 0.0d0
         do irow = max(1, jcol-ku_jac), min(neq, jcol+kl_jac)
            aij = ab(kl_jac+ku_jac+1+irow-jcol, jcol)
            if (finite_real(aij)) then
               colmax = max(colmax, abs(aij))
            else
               colmax = -1.0d0
               exit
            endif
         enddo
         colscan(jcol) = colmax
         if (colmax .gt. colbig) colbig = colmax
      enddo
      ! THE CENSUS BELOW IS OF THE OPERATOR THE COLUMN SCALING WAS CHOSEN
      ! FOR, not of the equilibrated one: after the equilibration every
      ! column is of size one and a column that is sampled in name only can
      ! no longer be told from one that is sampled.
      equilibrated = (column_equilibration_on .and. nspec_row .gt. 0 .and. &
                      allocated(scaled_column_norm_before))
      if (equilibrated) then
         colbig = 0.0d0
         do jcol = 1, neq
            if (colscan(jcol) .lt. 0.0d0) cycle
            colscan(jcol) = scaled_column_norm_before(jcol)
            colbig = max(colbig, colscan(jcol))
         enddo
      endif
      nz = 0;  nnf = 0;  nsmall = 0
      do jcol = 1, neq
         if (colscan(jcol) .lt. 0.0d0) then
            nnf = nnf + 1
         else if (colscan(jcol) .le. 0.0d0) then
            nz = nz + 1
         else if (colscan(jcol) .le. column_round_off*colbig) then
            nsmall = nsmall + 1
         endif
      enddo
      write(*,'(A,A,I0,A,I0,A,I0,A,I0)') tag, ' [diag 2] columns: zero ', &
           nz, ', nonfinite ', nnf, ', at the round-off of the band ',    &
           nsmall, ', unresolved colors ', n_unresolved
      ! --- 2b. how far the scaled operator's rows and columns spread ---
      ! ab holds A = Drow^-1 J D at this point, so these are the norms the
      ! band factorization pivots on and the Krylov space is built in. A
      ! column spread of many decades is what a badly conditioned band is
      ! made of, and it is the one the column scaling can remove.
      i_dmin = 1;  i_dmax = 1
      dmin = huge(1.0d0);  dmax = 0.0d0
      do jcol = 1, neq
         if (colscan(jcol) .lt. 0.0d0) cycle
         if (colscan(jcol) .lt. dmin) then
            dmin = colscan(jcol);  i_dmin = jcol
         endif
         if (colscan(jcol) .gt. dmax) then
            dmax = colscan(jcol);  i_dmax = jcol
         endif
      enddo
      call name_of_unknown(i_dmin, txt)
      write(*,'(A,A,ES11.3,A,A)') tag,                                    &
           ' [diag 2] scaled column norm min ', dmin,                     &
           ' before any equilibration, at the ', trim(txt)
      call name_of_unknown(i_dmax, txt)
      write(*,'(A,A,ES11.3,A,A)') tag,                                    &
           ' [diag 2] scaled column norm max ', dmax,                     &
           ' before any equilibration, at the ', trim(txt)
      rowscan = 0.0d0
      do jcol = 1, neq
         do irow = max(1, jcol-ku_jac), min(neq, jcol+kl_jac)
            aij = ab(kl_jac+ku_jac+1+irow-jcol, jcol)
            if (finite_real(aij))                                         &
               rowscan(irow) = max(rowscan(irow), abs(aij))
         enddo
      enddo
      i_rmin = 1;  i_rmax = 1
      rmin = huge(1.0d0);  rmax = 0.0d0
      do irow = 1, neq
         if (rowscan(irow) .lt. rmin) then
            rmin = rowscan(irow);  i_rmin = irow
         endif
         if (rowscan(irow) .gt. rmax) then
            rmax = rowscan(irow);  i_rmax = irow
         endif
      enddo
      call name_of_unknown(i_rmin, txt)
      write(*,'(A,A,ES11.3,A,A)') tag,                                    &
           ' [diag 2] scaled row    norm min ', rmin,                     &
           ' as factorized, at the equation of ', trim(txt)
      call name_of_unknown(i_rmax, txt)
      write(*,'(A,A,ES11.3,A,A)') tag,                                    &
           ' [diag 2] scaled row    norm max ', rmax,                     &
           ' as factorized, at the equation of ', trim(txt)
      if (equilibrated) then
         write(*,'(A,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') tag,          &
              ' [diag 2] column equilibration: factor ', colequil_min,    &
              ' .. ', colequil_max, ', column norms after ',              &
              colnorm_after_min, ' .. ', colnorm_after_max
      endif
      nshown = 0
      do jcol = 1, neq
         if (nshown .ge. 3) exit
         if (colscan(jcol) .le. column_round_off*colbig) then
            call name_of_unknown(jcol, txt)
            write(*,'(A,A,A,A,ES11.3)') tag, ' [diag 2]   the ',          &
                 trim(txt), ': largest entry ', colscan(jcol)
            nshown = nshown + 1
         endif
      enddo
      ! --- 3. the band factorization ---
      ! The 1-norm dgbcon needs is that of the matrix whose factors abf
      ! holds, so where the preconditioner was equilibrated the same two
      ! rescalings are applied here.
      anorm = 0.0d0
      do jcol = 1, neq
         s = 0.0d0
         do irow = max(1, jcol-ku_jac), min(neq, jcol+kl_jac)
            aij = ab(kl_jac+ku_jac+1+irow-jcol, jcol)
            if (irow .eq. jcol) then
               if (nspec_row .gt. 0) then
                  aij = aij + idtau*D(jcol)/Drow(jcol)
               else
                  aij = aij + idtau
               endif
            endif
            if (preconditioner_equilibrated)                              &
               aij = aij*precon_row_scale(irow)*precon_col_scale(jcol)
            s = s + abs(aij)
         enddo
         anorm = max(anorm, s)
      enddo
      umin = huge(1.0d0);  umax = 0.0d0;  i_piv = 1
      do jcol = 1, neq
         upiv = abs(abf(kl_jac+ku_jac+1, jcol))
         if (upiv .lt. umin) then
            umin = upiv;  i_piv = jcol
         endif
         umax = max(umax, upiv)
      enddo
      rcond = 0.0d0
      i = 0
      if (lu_info .eq. 0 .and. anorm .gt. 0.0d0 .and.                     &
          finite_real(anorm)) then
         allocate(work(3*neq), iwork(neq))
         call dgbcon('1', neq, kl_jac, ku_jac, abf, ldb, ipiv, anorm,     &
                     rcond, work, iwork, i)
         deallocate(work, iwork)
      endif
      call name_of_unknown(i_piv, txt)
      write(*,'(A,A,I0,A,ES11.3,A,ES11.3,A,A)') tag,                      &
           ' [diag 3] dgbtrf info ', lu_info, ', |U(j,j)| min ', umin,    &
           ' max ', umax, ', smallest at the ', trim(txt)
      write(*,'(A,A,ES11.3,A,ES11.3,A,ES11.3,A,I0)') tag,                 &
           ' [diag 3] 1-norm of the shifted band ', anorm,                &
           ', dgbcon rcond ', rcond, ', pivot ratio ',                    &
           umin/max(umax,1.0d-300), ', dgbcon info ', i
      deallocate(colscan, rowscan)
      end subroutine write_scales_columns_and_pivots

      ! ------------------------------------------------------!

      subroutine name_of_row_kind(k, what)
      ! The kind of unknown slot k within a cell: the three hydrodynamic
      ! rows, then the species rows in registry order.
      integer,          intent(in)  :: k
      character(len=*), intent(out) :: what
      integer :: is
      if (k .le. 3) then
         select case (k)
         case (1);      what = 'mass'
         case (2);      what = 'momentum'
         case default;  what = 'energy'
         end select
         return
      endif
      is = k - 3
      if (srow_kind(is) .eq. srow_carrier) then
         what = 'carrier '//trim(carrier_name(srow_idx(is)))
      else if (srow_kind(is) .eq. srow_element_he) then
         what = 'element He'
      else
         what = 'element '//trim(melem_name(srow_idx(is)))
      endif
      end subroutine name_of_row_kind

      ! ------------------------------------------------------!

      subroutine write_row_kind_shares(F, u, Drow, tag)
      ! WHICH ROW KIND HOLDS THE CONVERGENCE MEASURE AND WHICH HOLDS THE
      ! MERIT, side by side on one state.
      !
      ! The two are different functionals and the difference is the point of
      ! printing them together (certified_row_measures says why): ||R|| is a
      ! MAXIMUM over cells of each hydrodynamic row against its own largest
      ! term, and the elemental rows are not in it at all; the merit is the
      ! 2-norm of every row on the Newton's own row scales. A step control
      ! descends the second and the state is judged by the first, so a kind
      ! that holds almost all of the merit while another holds ||R|| is a
      ! statement about the scaling and not about the physics.
      !
      ! The elemental row measures are the products of the LAST residual
      ! evaluation (residual_evaluation_products), so this must be called
      ! while those describe the state F belongs to.
      !
      ! Diagnostic only (elem_diag_on); no residual evaluation.
      real*8, dimension(nvar_jac*N),  intent(in) :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8, dimension(nvar_jac*N),  intent(in) :: Drow
      character(len=*), intent(in) :: tag
      real*8  :: rc(3), rnorm, share(3+n_species_row_max)
      real*8  :: worst(3+n_species_row_max), q, total
      integer :: cell_worst(3+n_species_row_max)
      integer :: j, k, i, is, jw, k_merit, k_gate
      logical :: fin
      character(len=24) :: what, what_gate
      share = 0.0d0;  worst = 0.0d0;  cell_worst = 0
      do j = 1, N
         do k = 1, min(nvar_jac, 3 + n_species_row_max)
            i = nvar_jac*(j-1) + k
            q = abs(F(i))/max(Drow(i), 1.0d-300)
            share(k) = share(k) + q*q
            if (q .gt. worst(k)) then
               worst(k) = q;  cell_worst(k) = j
            endif
         enddo
      enddo
      total = sum(share(1:min(nvar_jac, 3 + n_species_row_max)))
      call resid_relnorm(F, u, rc, rnorm)
      write(*,'(A,A,ES11.3,A,ES11.3)') tag, ' [diag 7] ||R|| ', rnorm,    &
           ', squared merit ', total
      ! WHICH KIND HOLDS EACH OF THE TWO, in one line. Since the merit reads
      ! every row on its certification scale (cell_row_scales) the largest
      ! scaled row of the merit and the row ||R|| reports are the same row
      ! whenever the largest scaled row is a hydrodynamic one; where an
      ! element or carrier row is larger still, it is a row ||R|| does not
      ! carry at all and the gate reads it through cert_tol_element_at
      ! or cert_tol_carrier_at instead (distance_from_certification).
      k_merit = 1;  q = -1.0d0
      do k = 1, min(nvar_jac, 3 + n_species_row_max)
         if (worst(k) .gt. q) then
            q = worst(k);  k_merit = k
         endif
      enddo
      k_gate = 1
      do k = 2, 3
         if (rc(k) .gt. rc(k_gate)) k_gate = k
      enddo
      call name_of_row_kind(k_merit, what)
      call name_of_row_kind(k_gate, what_gate)
      write(*,'(A,A,A,A,ES11.3,A,A,A,ES11.3,A,L1)') tag,                  &
           ' [diag 7]   largest scaled row of the merit: ', trim(what),   &
           ' at ', worst(k_merit), '; ||R|| is the ', trim(what_gate),    &
           ' row at ', rnorm, '; same row ', (k_merit .eq. k_gate)
      do k = 1, min(nvar_jac, 3 + n_species_row_max)
         call name_of_row_kind(k, what)
         if (k .le. 3) then
            ! The cell belongs to the worst SCALED row, which is the
            ! merit's own offender; rows_over_window returns the maximum of
            ! the certification measure over two windows and not the cell it
            ! sits in, so none is claimed for it.
            write(*,'(A,A,A,A,F8.4,A,ES11.3,A,I0,A,F8.4,A,ES11.3)') tag,  &
                 ' [diag 7]   ', trim(what), ': merit share ',            &
                 share(k)/max(total,1.0d-300), ', worst scaled row ',     &
                 worst(k), ' at cell ', cell_worst(k),                    &
                 ' -> ||R|| share ', rc(k)/max(rnorm,1.0d-300),           &
                 ', its own measure ', rc(k)
         else
            is = k - 3
            ! The row's own certification measure, which is what the state
            ! is judged by; it is not in ||R||.
            q = 0.0d0;  jw = 0
            if (srow_kind(is) .eq. srow_element_he) then
               if (allocated(escale_he))                                  &
                  call certification_row_measure(N, j_min, erow_he(1:N),  &
                                                 escale_he(1:N), q, jw,   &
                                                 fin)
            else if (srow_kind(is) .ne. srow_carrier) then
               if (allocated(escale_tr) .and.                             &
                   erow_tr_carried(srow_idx(is)))                         &
                  call certification_row_measure(N, j_min,                &
                       erow_tr(1:N,srow_idx(is)),                         &
                       escale_tr(1:N,srow_idx(is)), q, jw, fin)
            else
               q = carrier_cellmax_last;  jw = carrier_cell_worst
            endif
            write(*,'(A,A,A,A,F8.4,A,ES11.3,A,I0,A,ES11.3,A,I0)') tag,    &
                 ' [diag 7]   ', trim(what), ': merit share ',            &
                 share(k)/max(total,1.0d-300), ', worst scaled row ',     &
                 worst(k), ' at cell ', cell_worst(k),                    &
                 ', its own measure ', q, ' at cell ', jw
         endif
      enddo
      end subroutine write_row_kind_shares

      ! ------------------------------------------------------!

      subroutine write_element_row_terms(rho, T, f_sp, Frho, res_he,      &
                                         sc_he, jlo, jhi)
      ! WHAT THE HELIUM ELEMENT ROW OF A CELL BALANCES, term by term.
      !
      ! The row is the stationary limit of the binary He/H composition
      ! equation (element_transport_residual): a diffusive divergence, which
      ! carries the concentration gradient and the settling drift on one face
      ! flux, against the divergence of the face element mass flux the wind
      ! carries.  THERE IS NO SOURCE: no reaction makes or destroys a helium
      ! nucleus, so the row is transport against transport and nothing else.
      !
      ! The operator returns the two terms summed.  The split is taken here
      ! by evaluating the SAME row on the SAME state with the face mass flux
      ! set to zero: the face element flux is F_rho times the face fraction
      ! and its divergence is exactly zero when F_rho is (species_face_flux,
      ! species_flux_divergence), so the second evaluation is the diffusive
      ! half alone and the difference of the two is the advective half.
      ! Nothing is approximated and no second copy of the operator exists.
      !
      ! The concentration and the settling halves of the diffusive term are
      ! inside ONE face flux (element_face_flux) and are not separable
      ! through the operator's interface; what is printed is their sum.
      !
      ! Also printed: the helium mass fraction the row is written in, its
      ! logarithmic radial gradient, the row's own scale (the sum of the
      ! magnitudes of the terms it balances) and its certification measure
      ! |res|/scale, which is what the judged distance divides by
      ! cert_tol_element_at.
      !
      ! Diagnostic only (element_row_terms_diag).  One extra evaluation of
      ! the element operator, which writes no persistent state (its own
      ! header states so), and no evaluation of the residual: the step the
      ! solve takes is the step it takes with the hook off.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho, T, Frho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1:N),                 intent(in) :: res_he, sc_he
      integer,                                intent(in) :: jlo, jhi

      real*8, dimension(1-Ng:N+Ng)           :: Frho_still
      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Yel
      real*8, dimension(1:N)         :: rhe_d, she_d
      real*8, dimension(1:N,n_melem) :: rtr_d, str_d
      logical, dimension(n_melem)    :: carried_d
      logical :: ok_he_d, ok_tr_d
      real*8  :: gradln, xhe, xhem, xhep
      integer :: j, j1, j2

      j1 = max(2, jlo);  j2 = min(N, jhi)
      if (j2 .lt. j1) return

      Frho_still = 0.0d0
      call element_transport_residual(rho, T, f_sp, Frho_still,           &
               rhe_d, she_d, ok_he_d, rtr_d, str_d, ok_tr_d,              &
               tr_carried = carried_d)
      if (.not. ok_he_d) return
      call element_mass_fractions(f_sp, Yel)

      write(*,'(A,I0,A,I0)') ' (JFNK) [diag 11] the helium element row,'  &
           //' cells ', j1, ' to ', j2
      do j = j1, j2
         xhe  = Yel(j,1)
         xhem = Yel(j-1,1)
         xhep = Yel(min(j+1,N+Ng),1)
         gradln = 0.0d0
         if (xhem .gt. 0.0d0 .and. xhep .gt. 0.0d0)                       &
            gradln = (log(xhep) - log(xhem))                              &
                    /max(log(r(min(j+1,N+Ng))) - log(r(j-1)), 1.0d-300)
         write(*,'(A,I4,A,F9.6,A,ES13.6,A,ES11.3,A,ES12.5,A,ES12.5,A,'    &
                 //'ES12.5,A,ES12.5,A,ES11.3)')                           &
              ' (JFNK) [diag 11] cell ', j, ': r ', r(j),                 &
              ', X_He ', xhe, ', dlnX_He/dlnr ', gradln,                  &
              ', diffusive divergence ', rhe_d(j),                        &
              ', advective divergence ', res_he(j) - rhe_d(j),            &
              ', row ', res_he(j), ', row scale ', sc_he(j),              &
              ', measure ', abs(res_he(j))/max(sc_he(j), 1.0d-300)
      enddo

      ! The wind the advective half rides on, at the same cells: a row whose
      ! advective term is the larger one is a statement about this flux and
      ! about the composition on its faces, not about the cell's own helium.
      do j = j1, j2
         write(*,'(A,I4,A,ES12.5,A,ES12.5)')                              &
              ' (JFNK) [diag 11] cell ', j, ': face mass flux in ',       &
              Frho(j-1), ', out ', Frho(j)
      enddo

      end subroutine write_element_row_terms

      ! ------------------------------------------------------!

      subroutine write_carrier_row_terms(f_sp, jlo, jhi)
      ! WHAT THE CARRIER ROW OF A CELL BALANCES, and the composition it is a
      ! statement about.
      !
      ! The carrier balance is transport against reaction
      ! (carrier_steady_residual), so unlike the element row it HAS a source,
      ! and the operator returns the row together with the sum of the
      ! magnitudes of the terms it balances (row_terms_phys), which is the
      ! scale the measure divides by.  The two are the products of the
      ! evaluation this is called from, so what is printed is the row of the
      ! state F belongs to.
      !
      ! Beside them the partition the row is about: the carrier fraction of
      ! the cell and its logarithmic radial gradient, and the neutral
      ! hydrogen fraction, because the front the carrier row binds at is also
      ! where the H I partition of the eliminated closure has more than one
      ! root (item N5).
      !
      ! Diagnostic only (element_row_terms_diag).  It evaluates nothing.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      integer,                                intent(in) :: jlo, jhi

      real*8  :: gradln, fc, fcm, fcp
      integer :: j, j1, j2, i, ic, isp

      if (.not. allocated(crow_res_last)) return
      j1 = max(1, jlo);  j2 = min(N, jhi)
      if (j2 .lt. j1) return

      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         ic  = srow_idx(i)
         isp = srow_isp(i)
         write(*,'(A,A,A,I0,A,I0)') ' (JFNK) [diag 12] the carrier ',     &
              trim(carrier_name(ic)), ' row, cells ', j1, ' to ', j2
         do j = j1, j2
            fc  = f_sp(j,isp)
            fcm = f_sp(j-1,isp)
            fcp = f_sp(min(j+1,N+Ng),isp)
            gradln = 0.0d0
            if (fcm .gt. 0.0d0 .and. fcp .gt. 0.0d0)                      &
               gradln = (log(fcp) - log(fcm))                             &
                       /max(log(r(min(j+1,N+Ng))) - log(r(j-1)), 1.0d-300)
            write(*,'(A,I4,A,F9.6,A,ES13.6,A,ES11.3,A,ES13.6,A,ES12.5,A,' &
                    //'ES12.5,A,ES11.3)')                                 &
                 ' (JFNK) [diag 12] cell ', j, ': r ', r(j),              &
                 ', carrier fraction ', fc, ', dln/dlnr ', gradln,        &
                 ', f(H I) ', f_sp(j,isp_HI),                             &
                 ', row ', crow_res_last(j,ic),                           &
                 ', row scale ', crow_terms_last(j,ic),                   &
                 ', measure ', abs(crow_res_last(j,ic))                   &
                              /max(crow_terms_last(j,ic), 1.0d-300)
         enddo
      enddo

      end subroutine write_carrier_row_terms

      ! ------------------------------------------------------!

      subroutine normal_equations_band(ab, bn)
      ! bn = J^T J in band form, from the band-stored J. A band of
      ! half-width kl_jac squared has half-width 2*kl_jac, and the product is
      ! symmetric, so bn is stored in the same LAPACK general-band layout
      ! with kl = ku = kl_jac_normal.
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N),           intent(in)  :: ab
      real*8, dimension(ldab_normal, nvar_jac*N),           intent(out) :: bn
      integer :: neq, i, j, k, klo, khi
      real*8  :: acc
      neq = nvar_jac*N
      bn  = 0.0d0
      do j = 1, neq
         do i = max(1, j-kl_jac_normal), min(neq, j+ku_jac_normal)
            acc = 0.0d0
            klo = max(1,   max(i-ku_jac, j-ku_jac))
            khi = min(neq, min(i+kl_jac, j+kl_jac))
            do k = klo, khi
               acc = acc + ab(kl_jac+ku_jac+1+k-i, i)*ab(kl_jac+ku_jac+1+k-j, j)
            enddo
            bn(kl_jac_normal+ku_jac_normal+1+i-j, j) = acc
         enddo
      enddo
      end subroutine normal_equations_band

      ! ------------------------------------------------------!

      subroutine levenberg_marquardt_step(bn, grad_merit, mu_damp, dZ, ok)
      ! One damped Gauss-Newton (Levenberg-Marquardt) step of the scaled
      ! model:  (J^T J + mu_damp*diag(J^T J)) dZ = -grad_merit,
      ! grad_merit = J^T s.  The matrix is symmetric positive definite for
      ! every mu_damp > 0 whenever diag(J^T J) > 0, so
      ! grad^T dZ = -grad_merit^T (.)^-1 grad_merit < 0: the step is a
      ! descent direction of the merit for ANY damping. mu_damp is the trust
      ! parameter -- mu_damp -> infinity recovers the scaled steepest-descent
      ! direction, mu_damp -> 0 the Gauss-Newton step. The Marquardt scaling
      ! (mu_damp times the diagonal, not mu_damp times the identity) makes
      ! the damping independent of the units of the individual unknowns.
      real*8, dimension(ldab_normal, nvar_jac*N), intent(in)  :: bn
      real*8, dimension(nvar_jac*N),                    intent(in)  :: grad_merit
      real*8,                                    intent(in)  :: mu_damp
      real*8, dimension(nvar_jac*N),                    intent(out) :: dZ
      logical,                                   intent(out) :: ok
      real*8,  allocatable :: bf(:,:)
      integer, allocatable :: ip(:)
      integer :: neq, ld, j, lpinfo
      neq = nvar_jac*N;  ld = ldab_normal
      allocate(bf(ld,neq), ip(neq))
      bf = bn
      do j = 1, neq
         bf(kl_jac_normal+ku_jac_normal+1, j) =                          &
              bf(kl_jac_normal+ku_jac_normal+1, j)                        &
            + mu_damp*abs(bn(kl_jac_normal+ku_jac_normal+1, j))
      enddo
      call dgbtrf(neq, neq, kl_jac_normal, ku_jac_normal, bf, ld, ip, lpinfo)
      if (lpinfo .ne. 0) then
         dZ = 0.0d0;  ok = .false.
         deallocate(bf, ip);  return
      endif
      dZ = -grad_merit
      call dgbtrs('N', neq, kl_jac_normal, ku_jac_normal, 1, bf, ld, ip,  &
                  dZ, neq, lpinfo)
      ok = (lpinfo .eq. 0)
      if (.not. ok) dZ = 0.0d0
      deallocate(bf, ip)
      end subroutine levenberg_marquardt_step

      ! ------------------------------------------------------!

      subroutine levenberg_marquardt_descent(Y, F, f_sp_base, D, Drow,    &
                                             ab, f2,                      &
                                             dY, f2_new, mu_used, gnorm,  &
                                             found)
      ! Find a step that LOWERS the merit when the Newton/PTC direction
      ! raises it.
      !
      ! The merit the solve minimizes is m(x) = 1/2*||s||^2 with s = F/D and
      ! x the scaled unknowns (Y = D*x); its gradient is A^T s with
      ! A = D^-1 J D, which the band-stored ab already holds. The Newton/PTC
      ! direction solves (I/dtau + A) dZ = -s -- an implicit Euler step, not a
      ! minimizer of m -- and on the He/H = 0.3 hand-off state it points
      ! UPHILL at every damping, including the dtau -> 0 limit dZ = -s, which
      ! is uphill exactly when s^T A s < 0 (docs/Update_EXHALE_stage1.md section
      ! 126). The damped Gauss-Newton family above descends for every mu, so
      ! what is left is choosing mu, and that is a one-dimensional search
      ! along the Levenberg-Marquardt curve: walk mu down by decades from a
      ! near-steepest-descent value while the TRUE merit keeps improving, and
      ! return the best step seen.
      !
      ! The trials are measured with weno_mode = 0, the same true merit the
      ! backtracking line search uses, so the two are directly comparable.
      real*8, dimension(nvar_jac*N),                    intent(in)  :: Y, F, D
      real*8, dimension(nvar_jac*N),                    intent(in)  :: Drow
      real*8, dimension(1-Ng:N+Ng,n_species),    intent(in)  :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8,                                    intent(in)  :: f2
      real*8, dimension(nvar_jac*N),                    intent(out) :: dY
      real*8,                                    intent(out) :: f2_new
      real*8,                                    intent(out) :: mu_used
      ! Norm of the merit gradient in the scaled space, ||A^T s||. Reported so
      ! that a search that finds nothing can be told apart from a stationary
      ! point, where there is nothing to find.
      real*8,                                    intent(out) :: gnorm
      logical,                                   intent(out) :: found
      real*8, allocatable :: bn(:,:)
      real*8, dimension(nvar_jac*N) :: s, grad_merit, dZ, dYtry, Ytry, Ftry
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng)   :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng) :: utry, Wtry
      real*8  :: mu_damp, f2t, f2prev
      integer :: im, saved_mode
      logical :: okstep, okres
      ! Damping sweep: mu_start is far enough on the steepest-descent side
      ! that the first step is short and safe, and n_mu decades reach the
      ! Gauss-Newton end. Measured on the He/H = 0.3 hand-off state, the merit
      ! falls monotonically from mu_damp = 1e2 down to 1e-4 (0.514 -> 0.260)
      ! and rises again at 1e-6.
      real*8,  parameter :: mu_start = 1.0d2
      integer, parameter :: n_mu     = 10
      dY = 0.0d0;  f2_new = f2;  mu_used = 0.0d0;  found = .false.
      saved_mode = weno_mode
      allocate(bn(ldab_normal, nvar_jac*N))
      if (nspec_row .gt. 0) then
         s = F/Drow
      else
         s = F/D
      endif
      call band_matvec_transpose(ab, s, grad_merit)
      gnorm = sqrt(sum(grad_merit*grad_merit))
      call normal_equations_band(ab, bn)
      f2prev = huge(1.0d0)
      mu_damp = mu_start
      do im = 1, n_mu
         call levenberg_marquardt_step(bn, grad_merit, mu_damp, dZ, okstep)
         if (okstep) then
            dYtry = D*dZ
            Ytry  = Y + dYtry
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                     &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               weno_mode = 0
               ! A POINT ON THE LEVENBERG-MARQUARDT CURVE. This routine
               ! returns a step and no state; the caller re-evaluates the
               ! one it decides to take.
               call eval_residual(Ytry, f_sp_base, fwork, Ftry, heat,     &
                                  cool,                                   &
                                  admissible=okres,                       &
                                  state_is_discarded=.true.)
               if (nspec_row .gt. 0) then
                  f2t = sqrt(sum((Ftry/Drow)**2))
               else
                  f2t = sqrt(sum((Ftry/D)**2))
               endif
               if (okres .and. f2t .lt. f2_new) then
                  dY = dYtry;  f2_new = f2t;  mu_used = mu_damp;  found = .true.
               endif
               ! The merit along the LM curve falls, reaches a minimum and
               ! rises again; once it has turned, smaller mu only overshoots
               ! further, so the sweep stops at the turn.
               if (okres .and. f2t .gt. f2prev) exit
               if (okres) f2prev = f2t
            endif
         endif
         mu_damp = 0.1d0*mu_damp
      enddo
      weno_mode = saved_mode
      deallocate(bn)
      end subroutine levenberg_marquardt_descent

      ! ------------------------------------------------------!


      ! ------------------------------------------------------!

      subroutine certified_row_measures(Fvec, u, rows, rnorm)
      ! THE MEASURES THE RETURNED STATE IS JUDGED BY, formed on one state
      ! (Fvec, u) exactly as the certification forms them:
      !
      !   rows(1)  the three hydrodynamic rows on residual_row_scale,
      !            CELL BY CELL over the whole physical column, each cell's
      !            measure over the tolerance that row carries AT THAT CELL
      !            and the largest of them taken
      !            (hydrodynamic_distance_from_certification_by_cell):
      !            momentum over 1e-8 and energy over 1e-6, one number each
      !            for the column, and the continuity row over the tolerance
      !            of its own cell, 3e-12 where the cell's arithmetic
      !            resolves it and ten times the cell's own rounding floor
      !            where it does not (cert_tol_mass_at, certification.f90).
      !            Those are the numbers hydro_row_entry gates with, taken
      !            the way mass_row_verdict takes them.
      !            It is already a distance, and is the slot that says
      !            WHICH hydrodynamic row refuses; the plain maximum of the
      !            same measures, which is |R| and the number the
      !            acceptance gate reads, is returned in rnorm beside it;
      !   rows(2)  the largest elemental transport row the registry carries,
      !            on the operator's own scale, which is the number
      !            transport_row_entry gates with cert_tol_element_at;
      !   rows(3)  the largest carrier row over the same cells, which
      !            carrier_row_entry gates with cert_tol_carrier_at.
      !
      ! THE TWO SPECIES ROWS ARE THE GATED PART OF THEIR COLUMN, r >=
      ! cert_regime_wind_r, because that is the part their tolerance is a
      ! statement about (decision 22 (a); the declaration of the two radii
      ! in certification.f90 says why the rest is reported and does not
      ! gate). The hydrodynamic row is the whole column, which is what its
      ! own tolerances are anchored on.
      !
      ! NONE OF THEM IS THE MERIT, and that is why this exists. The merit the
      ! step control minimizes is || F/Drow ||_2, a 2-norm on the Newton's
      ! own column scales; a state is judged by a MAXIMUM over cells of each
      ! row against its own largest term. The two are different functionals
      ! on different scalings, so a rule that bounds the first bounds nothing
      ! about the second: MEASURED on the hot Uranus reload with an element
      ! row, the default route drove the merit from 1.44e+02 to 3.86e+01
      ! while ||R|| went 1.72, 6.44, 1.82, 25.2, 1.93 (B5d section 2.7).
      real*8, dimension(nvar_jac*N),  intent(in)  :: Fvec
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rows
      ! |R| itself, the plain maximum over the three hydrodynamic rows and
      ! the two windows. The gate reads this number and the ledger reads
      ! rows(1), and the two are different functionals, so a caller that
      ! prints both asks for both here rather than re-deriving either.
      real*8, optional,               intent(out) :: rnorm
      real*8  :: rc(3), q, dg, rn
      ! The three hydrodynamic rows of EVERY cell, and the continuity row's
      ! rounding floor at each of them: the tolerance of that row is a
      ! function of the cell, so the distance is not a functional of the
      ! three row maxima and the cells are kept.
      real*8  :: qhyd(3,N), floor_mass(N)
      integer :: i, j, k, jw, ic
      rows = 0.0d0
      ! |R_kj|/residual_row_scale(k,j), the measure the certification takes
      ! of a cell (the expression mass_row_verdict takes cell by cell), over
      ! the whole physical column [1:N]: that is the union of the wind and
      ! the layer resid_relnorm norms separately, and a maximum over a union
      ! is the larger of the two window maxima.
      do j = 1, N
         do k = 1, 3
            qhyd(k,j) = abs(Fvec(nvar_jac*(j-1)+k))                      &
                        /max(residual_row_scale(k, j, u), cert_scale_floor)
         enddo
         floor_mass(j) = mass_row_rounding_floor(j, u)
      enddo
      rows(1) = hydrodynamic_distance_from_certification_by_cell(N,       &
                                                     qhyd, floor_mass)
      ! |R| of the same evaluation, when it is asked for: the plain maximum
      ! of the same measures with no tolerance in it, of which resid_relnorm
      ! stays the definition.
      if (present(rnorm)) then
         call resid_relnorm(Fvec, u, rc, rn)
         rnorm = rn
      endif
      if (nspec_row .le. 0) return
      do i = 1, nspec_row
         select case (srow_kind(i))
         case (srow_carrier)
            if (.not. allocated(crow_res_last)) cycle
            ic = srow_idx(i)
            call certification_species_row_gate(N, .true.,               &
                     crow_res_last(1:N,ic), crow_terms_last(1:N,ic),     &
                     dg, q, jw)
            if (jw .gt. 0) rows(3) = max(rows(3), q)
         case (srow_element_he)
            call certification_species_row_gate(N, .false., erow_he(1:N), &
                     escale_he(1:N), dg, q, jw)
            if (jw .gt. 0) rows(2) = max(rows(2), q)
         case default
            if (erow_tr_carried(srow_idx(i))) then
               call certification_species_row_gate(N, .false.,           &
                        erow_tr(1:N,srow_idx(i)),                        &
                        escale_tr(1:N,srow_idx(i)), dg, q, jw)
               if (jw .gt. 0) rows(2) = max(rows(2), q)
            endif
         end select
      enddo
      end subroutine certified_row_measures

      ! ------------------------------------------------------!

      pure real*8 function                                               &
             hydrodynamic_distance_from_certification_by_cell            &
                      (nc, q, floor_mass) result(d)
      ! HOW FAR THE THREE HYDRODYNAMIC ROWS OF A COLUMN ARE FROM BEING
      ! CERTIFIED, and the ONE expression of that distance: the largest,
      ! over cells AND over the three rows, of a cell's row measure over the
      ! tolerance that row carries at that cell.
      !
      !    d = max_j max_k  q(k,j) / tol_k(j)
      !
      ! q(k,j) = |R_kj|/residual_row_scale(k,j) is the measure the
      ! certification takes of a cell. tol_2 and tol_3 are one number each
      ! for the column (cert_tol_momentum 1e-8, cert_tol_energy 1e-6);
      ! tol_1(j) is the continuity row's tolerance AT CELL j, 3e-12 where
      ! the cell's arithmetic resolves it and ten times the cell's own
      ! rounding floor where it does not, and it is read from
      ! mass_row_cell_verdict, which is the one statement of that rule and
      ! returns the distance itself rather than a tolerance to divide by.
      !
      ! ONE TOLERANCE PER ROW, NOT ONE FOR THE THREE. The three tolerances
      ! are anchored one by one on a converged state, where the rows
      ! themselves stand decades apart on residual_row_scale
      ! (certification.f90: mass 3.0e-12, momentum 1.0e-8, energy 1.0e-6),
      ! and hydro_row_entry gates each row with its own. The plain maximum
      ! of the three measures, |R|, what the acceptance gate reads against
      ! one number, therefore says nothing about which row refuses: MEASURED
      ! on the HD 209458 b element reload, every solve of the partitioned
      ! route ends with |R| at 1e-8 held by the energy row while the mass row
      ! stands at 1e-9 against 3.0e-12, so |R|/tol is about 1 for every
      ! iterate and the mass row, hundreds of times outside the fixed
      ! tolerance, is invisible in it.
      !
      ! WHY THE ROWS CANNOT BE MAXIMIZED OVER CELLS FIRST. The maximum of a
      ! ratio is not the ratio of the maxima once the divisor moves with the
      ! cell: a base cell whose mass row is 1e-9 against its own rounding
      ! floor stands well inside its tolerance while a wind cell at 1e-11
      ! against 3e-12 refuses, and the row maximum is the base cell's. A
      ! triple of row maxima divided by cert_tol_mass, which is the FLOOR of
      ! the continuity tolerance and the tolerance of no particular cell,
      ! reports 333 where the condition the acceptance applies is 3.3; it
      ! stands at or above this distance, by the ratio of the cell's floor to
      ! 3e-12, MEASURED at up to 3.2e+03 on the HD 209458 b base layer
      ! (docs/certification_tolerance_anchoring_20260910.md, anchor 6).
      !
      ! THE TOLERANCE IS A FUNCTION OF THE STATE, through the floor, so two
      ! iterates are each ranked against the arithmetic of their own base
      ! layer. That is the same rule the acceptance applies to the state it
      ! is handed (mass_row_verdict), which is what makes d < 1 exactly the
      ! hydrodynamic part of stationary_rows_of_the_returned_state.
      integer,                 intent(in) :: nc
      real*8, dimension(3,nc), intent(in) :: q
      real*8, dimension(nc),   intent(in) :: floor_mass
      real*8  :: tol, dist
      logical :: within, anchored
      integer :: j
      d = 0.0d0
      do j = 1, nc
         call mass_row_cell_verdict(q(1,j), floor_mass(j), tol, dist,    &
                                    within, anchored)
         d = max(d, dist)
         d = max(d, q(2,j)/cert_tol_momentum)
         d = max(d, q(3,j)/cert_tol_energy)
      enddo
      end function hydrodynamic_distance_from_certification_by_cell

      ! ------------------------------------------------------!

      real*8 function distance_from_certification(rows) result(d)
      ! HOW FAR A STATE IS FROM BEING CERTIFIED, in one number: each judged
      ! row divided by the tolerance that row is judged against, and the
      ! largest of those. d < 1 is exactly the condition
      ! stationary_rows_of_the_returned_state states on the rows the
      ! registry carries.
      !
      ! The division by the tolerances is what makes the slots commensurate:
      ! a hydrodynamic row is measured against a scale that bounds its own
      ! largest term and a converged state sits at 1e-6 to 1e-13 on it row
      ! by row, while a carrier row measured against its own terms sits near
      ! 1e-3 (resid_relnorm says why they may not simply be maximized
      ! together). Divided by their own tolerances they are one quantity.
      !
      ! THE HYDRODYNAMIC SLOT ARRIVES ALREADY DIVIDED, cell by cell and by
      ! the tolerance each of its three rows carries at that cell rather
      ! than by one number
      ! (hydrodynamic_distance_from_certification_by_cell, which says why
      ! one number cannot do it), so it is taken as it stands.
      real*8, dimension(3), intent(in) :: rows
      d = rows(1)
      d = max(d, rows(2)/cert_tol_element_at(cert_regime_wind_r))
      d = max(d, rows(3)/cert_tol_carrier_at(cert_regime_wind_r))
      end function distance_from_certification

      ! ------------------------------------------------------!

      subroutine resid_relnorm(F, u, rc, rnorm)
      ! THE CONVERGENCE MEASURE of the steady solve: the relative residual of
      ! each row over the WHOLE physical column [1:N], every cell divided by
      ! residual_row_scale (steady_residual.f90) -- |u| for the mass and
      ! energy rows, |rho v| for momentum in the wind and the gravitational
      ! force density for momentum below the escape radius, where the layer
      ! is quasi-hydrostatic and |rho v| carries no information.
      !
      ! The wind [j_min:N] and the layer [1:j_min-1] are normed SEPARATELY
      ! and combined by the LARGER of the two. Until section 127 this measure
      ! looked at [j_min:N] alone, i.e. at r > r_esc, and the marching du it
      ! takes over from still does; the column between the base and r_esc was
      ! then controlled by nothing, and two solves of one hand-off state both
      ! returned info = 0 on winds differing by 0.22 dex in the rate, the
      ! WORSE of them scoring the better ||R||. Combining by the maximum,
      ! rather than by a volume- or cell-count-weighted sum, is what keeps
      ! that from recurring: the outer region carries almost all of
      ! sum r^2 dr, so any sum lets it average the inner column away.
      real*8, dimension(nvar_jac*N),         intent(in)  :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8,                         intent(out) :: rnorm
      real*8, dimension(3,1-Ng:N+Ng) :: Rres
      real*8, dimension(3)           :: rc_wind, rc_layer
      integer :: j, k
      Rres = 0.0d0
      do j = 1, N
         do k = 1, 3
            Rres(k,j) = F(nvar_jac*(j-1)+k)
         enddo
      enddo
      ! ONE DEFINITION OF THE ROW MEASURE. The gate and the certification
      ! used to form max_j |R_kj|/residual_row_scale(k,j) from two separate
      ! expressions, and only a run-time assertion kept them together (see
      ! gate_equals_certification below, which stays as the guard). The
      ! expression is certification_row_measure's, and this is the one call
      ! site that used the other copy; relnorm_over_cells keeps its own
      ! callers in steady_residual.
      call rows_over_window(Rres, u, j_min, N,       rc_wind)
      call rows_over_window(Rres, u, 1,     j_min-1, rc_layer)
      rc    = max(rc_wind, rc_layer)
      rnorm = maxval(rc)
      ! THE CARRIER ROW IS NOT FOLDED IN HERE, and that is the same
      ! judgment section 133 made about the flux gate: one number cannot say
      ! two things. `Resid tol` is calibrated against a row scale that BOUNDS
      ! each hydrodynamic row's largest term, and converged states sit at
      ! 1e-6 on it; the carrier row is measured against the terms
      ! themselves, where a converged row sits near 1e-3. Taking the maximum
      ! would hold the carrier row to a tolerance meant for a different
      ! scale. It is a third gate instead (steady_gates_met).
      end subroutine resid_relnorm

      ! ------------------------------------------------------!

      subroutine rows_over_window(Res, u, ja, jb, rc)
      ! The three hydrodynamic rows of Res over the cells [ja:jb], each
      ! measured by certification_row_measure on its own row scale. This
      ! gate is the maximum over cells, so the volume-weighted companion is
      ! not asked for and is not formed: it is optional there, and a JFNK
      ! solve evaluates this measure at every residual.  jsplit is set past
      ! the end because with no companion the split has nothing to divide.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: Res, u
      integer,                        intent(in)  :: ja, jb
      real*8, dimension(3),           intent(out) :: rc
      real*8, dimension(N) :: rr, ss
      integer :: k, j, nc, jworst
      logical :: ok_finite
      rc = 0.0d0
      if (jb .lt. ja) return
      nc = jb - ja + 1
      do k = 1, 3
         do j = 1, nc
            rr(j) = Res(k, ja + j - 1)
            ss(j) = residual_row_scale(k, ja + j - 1, u)
         enddo
         call certification_row_measure(nc, nc + 1, rr(1:nc), ss(1:nc),   &
                                        rc(k), jworst, ok_finite)
      enddo
      end subroutine rows_over_window

      ! ------------------------------------------------------!

      subroutine resid_relnorm_below_escape(F, u, rc, j_worst, r_worst, &
                                            k_worst)
      ! The relative residual of the layer [1:j_min-1] BELOW the escape
      ! radius on its own, and the cell in it that carries the largest
      ! cell-wise scaled residual.
      !
      ! Since section 127 this layer is one of the two halves of the
      ! convergence measure (resid_relnorm returns the larger of it and the
      ! wind), so what this adds to that number is WHERE in the layer the
      ! residual sits -- which the volume-weighted sum averages away -- and
      ! how the layer stands against the wind at the same iterate.
      !
      ! The row scales are residual_row_scale's: |u| for mass and energy,
      ! and for momentum the gravitational force density
      !    s_grav(j) = |u(1,j)| * |Gphi_i(j) - Gphi_i(j-1)| / dr_j(j)
      ! (the same discrete potential difference the source term uses;
      ! Source.f90), so rc(2) is the fractional violation of hydrostatic
      ! balance in the layer.
      real*8, dimension(nvar_jac*N),         intent(in)  :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8,                         intent(out) :: r_worst
      integer,                        intent(out) :: j_worst, k_worst
      real*8, dimension(3,1-Ng:N+Ng) :: Rres
      integer :: j, k
      real*8  :: cellrel, amx, uscale

      rc = 0.0d0;  j_worst = 0;  k_worst = 0;  r_worst = 0.0d0
      if (j_min .le. 1) return          ! no cells below the escape radius

      Rres = 0.0d0
      do j = 1, j_min-1
         do k = 1, 3
            Rres(k,j) = F(nvar_jac*(j-1)+k)
         enddo
      enddo
      call rows_over_window(Rres, u, 1, j_min-1, rc)

      amx = -1.0d0
      do j = 1, j_min-1
         do k = 1, 3
            uscale  = residual_row_scale(k, j, u)
            cellrel = abs(F(nvar_jac*(j-1)+k))/max(uscale, 1.0d-30)
            if (cellrel .gt. amx) then
               amx = cellrel;  j_worst = j;  k_worst = k
            endif
         enddo
      enddo
      if (j_worst .gt. 0) r_worst = r(j_worst)

      end subroutine resid_relnorm_below_escape

      ! ------------------------------------------------------!

      subroutine write_resid_below_escape(tag, F, u)
      ! One line for the layer below the escape radius, which is one of the
      ! two halves of the convergence measure (see resid_relnorm_below_escape
      ! and resid_relnorm). Silent when the layer is empty. (F, u) must be a
      ! consistent pair: the caller passes the residual of the state
      ! currently held in u. The three numbers are the volume-weighted row
      ! values: mass and energy relative to |u| (rates), momentum relative to
      ! the gravitational force density (fractional violation of hydrostatic
      ! balance in the layer). The trailing cell is where the largest
      ! cell-wise scaled residual of the layer sits.
      character(*),                   intent(in) :: tag
      real*8, dimension(nvar_jac*N),         intent(in) :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8  :: rc_below(3), r_worst
      integer :: j_worst, k_worst
      if (j_min .le. 1) return
      call resid_relnorm_below_escape(F, u, rc_below, j_worst, r_worst,   &
                                      k_worst)
      write(*,'(A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,F8.4,A,I0)')     &
           tag, ' below r_esc [1:',j_min-1,                               &
           '] |R|: mass=',rc_below(1), ' mom/grav=',rc_below(2),          &
           ' energy=',rc_below(3),                                        &
           '  max cell j=',j_worst,' r=',r_worst,' k=',k_worst
      end subroutine write_resid_below_escape

      ! ------------------------------------------------------!

      subroutine reset_no_chem_root_trial_count
      n_trial_no_chem_root  = 0
      no_chem_root_cell_max = 0
      no_chem_root_res_max  = 0.0d0
      n_no_chem_root_state = 0
      end subroutine reset_no_chem_root_trial_count

      ! ------------------------------------------------------!

      subroutine write_no_chem_root_trial_count(tag)
      ! What the solve refused, and how far outside a root those states were.
      ! Silent when every probe and trial landed on a root of the network.
      character(*), intent(in) :: tag
      if (n_trial_no_chem_root .le. 0) return
      write(*,'(A,A,I0,A,I0,A,ES10.3)') tag,                             &
           ' trial states refused for a composition without a chemical'//  &
           ' root: ',                                                     &
           n_trial_no_chem_root, '; worst cell count ', no_chem_root_cell_max,&
           ', worst reaction residual ', no_chem_root_res_max
      end subroutine write_no_chem_root_trial_count

      ! ------------------------------------------------------!

      ! ------------------------------------------------------!

      subroutine certify_returned_state(F, u, f_sp, resid_tol,            &
                                        n_no_chem_root, label)
      ! THE STATIONARY CERTIFICATION CONTEXT at a point where this solver
      ! declares a state solved (contract section 3, row 4). It measures
      ! every equation the configuration makes independent on the state
      ! about to be handed back, and it decides nothing about the solve: the
      ! merit rule, the step control and info are untouched. A state that
      ! refuses certification is still written; what it may not do is be
      ! called certified.
      real*8, dimension(:),                   intent(in) :: F
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,  intent(in) :: resid_tol
      integer, intent(in) :: n_no_chem_root
      character(len=*), intent(in) :: label
      real*8, dimension(3,1-Ng:N+Ng) :: Rrows
      type(cert_report) :: rep
      integer :: jj, kk
      Rrows = 0.0d0
      do jj = 1, N
         do kk = 1, 3
            Rrows(kk,jj) = F(nvar_jac*(jj-1)+kk)
         enddo
      enddo
      call certification_evaluate(cert_context_stationary, u, Rrows,      &
                                  f_sp, resid_tol, n_no_chem_root,        &
                                  .true., rep)
      call certification_report_write(rep, label)
      end subroutine certify_returned_state

      ! ------------------------------------------------------!

      subroutine measure_the_closure_map_at_the_base(Y, f_sp_base, Drow)
      ! THE MAP FROM THE INNER CLOSURE'S ERROR TO THE OUTER RESIDUAL AT THE
      ! BASE, measured on the iterate and not assumed.
      !
      ! WHY IT IS NOT ONE. eval_residual eliminates the composition by
      ! sweeps and stops when the change of what the residual reads --
      ! particle count, H2 share of it, heating and cooling -- falls below
      ! the tolerance the sweep is asked for. That is an error in the
      ! COMPOSITION. At the base it reaches the energy row twice: once
      ! through the equation of state of the cell, and once through the
      ! ghost, whose temperature is built from the particle count refreshed
      ! just before the boundary condition is applied. So the base row can
      ! carry many times the closure error, and a bound of the form
      ! "tolerance times the merit" on the outer residual is not a bound.
      !
      ! WHAT IS MEASURED, in four evaluations of one state.
      !
      ! The first settles the background the modules carry, so that what
      ! follows is not the transient of the state the caller left. The
      ! second is the production map, and the increment its sweep exited on
      ! is its remaining closure error in the sweep's own measure. The third
      ! repeats it, and the change between the two is the REPRODUCIBILITY:
      ! what two evaluations of one state differ by for reasons that are not
      ! the closure. The fourth runs three passes MORE than the second took,
      ! at a fixed count, which is a tighter closure whatever the tolerance
      ! is doing, and the change from the second is the closure error's
      ! footprint on the outer rows.
      !
      ! THE REFERENCE IS A PASS COUNT AND NOT A TOLERANCE, and it has to be:
      ! MEASURED on the atomic element reload the sweep exits after ONE pass
      ! at an increment of 9.0e-12, four decades below the 1e-8 it is asked
      ! for, so asking it for less does not buy a pass and a tolerance-based
      ! reference measures nothing.
      !
      ! It is taken ONCE per solve, on the first iterate. Its probes belong
      ! to the candidate ledger and not to the run's acceptance statistics,
      ! and the products they leave in the modules are put back, as every
      ! measurement in this module does.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heatw, coolw
      real*8, allocatable :: Fa(:), Fb(:), Fc(:)
      real*8  :: closure_error, base_change, all_change
      real*8  :: base_repro, amplification, q
      integer :: neq, k, i, na
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry
      if (nspec_row .le. 0) return
      neq = nvar_jac*N
      allocate(Fa(neq), Fb(neq), Fc(neq))
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      call eval_residual(Y, f_sp_base, fwork, Fb, heatw, coolw,           &
                         may_be_adopted=.false.,                          &
                         state_is_discarded=.true.)
      call eval_residual(Y, f_sp_base, fwork, Fa, heatw, coolw,           &
                         may_be_adopted=.false.,                          &
                         state_is_discarded=.true.)
      closure_error = eq_sweep_increment_last
      na            = n_eq_sweeps_last
      call eval_residual(Y, f_sp_base, fwork, Fc, heatw, coolw,           &
                         may_be_adopted=.false.,                          &
                         state_is_discarded=.true.)
      call eval_residual(Y, f_sp_base, fwork, Fb, heatw, coolw,           &
                         may_be_adopted=.false.,                          &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=na+3)
      base_change = 0.0d0
      base_repro  = 0.0d0
      do k = 1, min(3, nvar_jac)
         base_change = max(base_change,                                   &
                           abs(Fa(k) - Fb(k))/max(Drow(k), 1.0d-300))
         base_repro  = max(base_repro,                                    &
                           abs(Fa(k) - Fc(k))/max(Drow(k), 1.0d-300))
      enddo
      all_change = 0.0d0
      do i = 1, neq
         all_change = max(all_change,                                     &
                          abs(Fa(i) - Fb(i))/max(Drow(i), 1.0d-300))
      enddo
      ! A CHANGE AT OR BELOW THE REPRODUCIBILITY IS NOT A CLOSURE ERROR: the
      ! extra passes then moved the rows no further than a second evaluation
      ! of the same state does, and there is nothing to amplify. The
      ! amplification stays at one and the line below says so.
      q             = 0.0d0
      amplification = 1.0d0
      if (closure_error .gt. 0.0d0 .and. base_change .gt. base_repro) then
         q = base_change/closure_error
         ! A QUOTIENT BELOW ONE IS AN ATTENUATION AND NOT AN AMPLIFICATION,
         ! and then the target is already a bound on what the outer residual
         ! keeps: the sweep is asked for the target itself. The quotient is
         ! printed either way, because which of the two it is, is the
         ! measurement.
         if (finite_real(q)) amplification = max(1.0d0, q)
      endif
      closure_amplification_at_the_base = amplification
      eq_sweep_reltol = outer_residual_closure_target/amplification
      write(*,'(A,ES11.3,A,I0,A,I0,A)') ' (JFNK) closure map: the sweep'// &
           ' exited at increment ', closure_error, ' after ', na,          &
           ' passes; the reference ran ', na+3, ' passes'
      write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)') ' (JFNK) closure map:'//     &
           ' scaled row change at the base ', base_change,                 &
           ', its reproducibility ', base_repro,                           &
           ', over the whole column ', all_change
      write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)') ' (JFNK) closure map:'//     &
           ' the base carries ', q, ' of the closure error, so the'//      &
           ' amplification is ', amplification,                            &
           ' and the sweep is asked for ', eq_sweep_reltol
      call put_back_residual_evaluation_products(products_at_entry)
      call put_back_solve_refusal_statistics(statistics_at_entry)
      deallocate(Fa, Fb, Fc)
      end subroutine measure_the_closure_map_at_the_base

      ! ------------------------------------------------------!

      subroutine species_row_jacobian_action_check(Y, f_sp_base)
      ! J r AGAINST [F(Y + eps r) - F(Y)]/eps ON THE SYSTEM THIS SOLVE
      ! CARRIES, to the standard the run-level self-test states: a correct
      ! band layout and coloring match to about 1e-6, a wrong one is out by
      ! O(1).
      !
      ! IT IS HERE AND NOT AT THE TOP OF THE RUN because the row registry is
      ! empty until a solve is entered: the run-level self-test therefore
      ! measures the three-unknown system in every configuration and can say
      ! nothing about a species row, its column, or the widened coloring
      ! stride.
      !
      ! IT TESTS THE JACOBIAN THE SOLVE USES, build_banded_jac_full on
      ! eval_residual, and not the frozen one: the frozen residual states no
      ! species balance at all (its species block is empty by construction),
      ! so a match there would say nothing about these rows. The comparison
      ! is exact in the band only up to the weakly non-local radiation
      ! coupling the full residual carries, which is why the reference is
      ! the directional difference of the SAME map -- the same fixed number
      ! of composition passes at both points, as the columns use.
      !
      ! The direction is scaled by the unknowns' own column scales, so the
      ! perturbation of a species slot is a perturbation of that species and
      ! not of the hydrodynamic magnitude that happens to share the vector.
      !
      ! Off unless EXHALE_SPECIES_JAC_TEST=1: it costs one Jacobian build
      ! and two full residual evaluations.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, allocatable :: abj(:,:), F0(:), Fp(:), rdir(:), Jr(:), Yp(:)
      real*8, allocatable :: Dsc(:)
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heatw, coolw
      real*8            :: eps_fd, err, jrmax, errsp, jrmaxsp
      integer           :: neq, i, j, k, nunres, saved_wmode, istep
      integer           :: jcell_test, i0_test, ios_test
      logical           :: okp
      character(len=32) :: env
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry
      call get_environment_variable('EXHALE_SPECIES_JAC_TEST', env)
      if (trim(env) .ne. '1') return
      ! A SELF-TEST MAY NOT CHANGE THE RUN IT MEASURES. Its three residual
      ! evaluations are probes, so their sweeps belong to the candidate
      ! ledger and not to the run's own acceptance statistics, and the
      ! products they leave in the modules are put back at the end. MEASURED
      ! before this was so: the tested and the untested run agreed to every
      ! printed digit for five outer iterations and were on different paths
      ! by iteration 13, one converging and one stalling.
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
      neq = nvar_jac*N
      allocate(abj(2*kl_jac+ku_jac+1, neq))
      allocate(F0(neq), Fp(neq), rdir(neq), Jr(neq), Yp(neq), Dsc(neq))
      call cell_state_scales(Y, Dsc)
      ! THE ITERATE'S OWN REFERENCES FIRST, in the mode that declares this
      ! state to be the iterate: the element budget a probe is compared
      ! against is the iterate's, and without this evaluation the reference
      ! is whatever the caller left, so every color would be refused for
      ! standing where the iterate stands.
      saved_wmode = weno_mode
      weno_mode   = 1
      call eval_residual(Y, f_sp_base, fwork, F0, heatw, coolw)
      weno_mode   = 2
      call build_banded_jac_full(Y, f_sp_base, abj, nunres)
      call eval_residual(Y, f_sp_base, fwork, F0, heatw, coolw,           &
                         admissible=okp, may_be_adopted=.false.,          &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      ! A deterministic direction that excites every row of every cell.
      do i = 1, neq
         rdir(i) = sin(0.1d0*dble(i))
      enddo
      rdir   = rdir/maxval(abs(rdir))
      rdir   = rdir*Dsc
      call band_matvec(abj, rdir, Jr)
      eps_fd = sqrt(epsilon(1.0d0))
      Yp     = Y + eps_fd*rdir
      call eval_residual(Yp, f_sp_base, fwork, Fp, heatw, coolw,          &
                         admissible=okp, may_be_adopted=.false.,          &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      Fp     = (Fp - F0)/eps_fd
      jrmax  = maxval(abs(Jr))
      err    = maxval(abs(Fp - Jr))/max(jrmax, 1.0d-30)
      ! The species rows alone, so that a mismatch confined to them cannot
      ! be averaged away by the hydrodynamic rows that dominate |J r|.
      errsp = 0.0d0;  jrmaxsp = 0.0d0
      do j = 1, N
         do k = 1, nspec_row
            i = nvar_jac*(j-1) + 3 + k
            jrmaxsp = max(jrmaxsp, abs(Jr(i)))
            errsp   = max(errsp,   abs(Fp(i) - Jr(i)))
         enddo
      enddo
      if (jrmaxsp .gt. 0.0d0) errsp = errsp/jrmaxsp
      write(*,'(A,I0,A,I0,A,I0,A,I0)') ' (species_jac_test) nvar = ',     &
           nvar_jac, ', species rows = ', nspec_row, ', kl = ku = ',      &
           kl_jac, ', colors = ', ncolor_jac
      write(*,'(A,ES12.4,A,ES12.4,A,I0)') ' (species_jac_test) '//        &
           '||J r - dFD||_inf / ||J r||_inf: all rows ', err,             &
           ', species rows ', errsp, ', unresolved colors ', nunres
      ! THE SAME DERIVATIVE AT ONE CELL, AGAINST THE PERTURBATION SIZE.
      ! One number at one step length cannot say whether a mismatch is the
      ! Jacobian or the difference: a correct pair falls like the step until
      ! the closure's error divided by the step takes over, and the least
      ! error is where the two meet. The direction is supported on ONE
      ! physical cell alone, so the rows read are that cell's.
      !
      ! The default cell is the first physical one, the coupled boundary the
      ! composition, the particle count, the caloric map and the ghost
      ! reconstruction all reach at once. EXHALE_JAC_TEST_CELL=j moves it to
      ! cell j, which is how the curve is read at the cell that BINDS a
      ! solve: the least error there is the smallest row measure a Newton
      ! step can resolve at that cell, and no tolerance below it is
      ! reachable by this method whatever the state
      ! (docs/certification_tolerance_anchoring_20260910.md).
      jcell_test = 1
      call get_environment_variable('EXHALE_JAC_TEST_CELL', env)
      if (len_trim(env) .gt. 0) then
         read(env, *, iostat = ios_test) jcell_test
         if (ios_test .ne. 0) jcell_test = 1
         jcell_test = max(1, min(N, jcell_test))
      endif
      i0_test = nvar_jac*(jcell_test - 1)
      do i = 1, neq
         rdir(i) = 0.0d0
      enddo
      do k = 1, nvar_jac
         rdir(i0_test + k) = Dsc(i0_test + k)
      enddo
      call band_matvec(abj, rdir, Jr)
      jrmax = 0.0d0
      do k = 1, nvar_jac
         jrmax = max(jrmax, abs(Jr(i0_test + k)))
      enddo
      write(*,'(A,I0,A,ES12.5,A,I0,A)') ' (species_jac_test) cell ',      &
           jcell_test, ' at r = ', r(jcell_test), ' alone, ', nvar_jac,   &
           ' rows: step, error of the directional difference'
      do istep = 1, 7
         eps_fd = sqrt(epsilon(1.0d0))*(1.0d2**(3-istep))
         Yp = Y + eps_fd*rdir
         call eval_residual(Yp, f_sp_base, fwork, Fp, heatw, coolw,       &
                            admissible=okp, may_be_adopted=.false.,       &
                            state_is_discarded=.true.,                    &
                            n_eq_sweeps_fixed=n_eq_sweeps_model)
         errsp = 0.0d0
         do k = 1, nvar_jac
            errsp = max(errsp, abs((Fp(i0_test + k) - F0(i0_test + k))    &
                                   /eps_fd - Jr(i0_test + k)))
         enddo
         write(*,'(A,ES10.3,A,ES12.4,A,ES12.4)')                          &
              ' (species_jac_test)   step ', eps_fd,                      &
              '  error ', errsp, '  relative ', errsp/max(jrmax, 1.0d-30)
      enddo
      write(*,'(A,ES11.3,A,ES11.3)') ' (species_jac_test) the closure'//  &
           ' amplification the base carries ',                            &
           closure_amplification_at_the_base,                             &
           ', so the sweep is asked for ', eq_sweep_reltol
      weno_mode = saved_wmode
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      call put_back_residual_evaluation_products(products_at_entry)
      call put_back_solve_refusal_statistics(statistics_at_entry)
      deallocate(abj, F0, Fp, rdir, Jr, Yp, Dsc)
      end subroutine species_row_jacobian_action_check

      ! ------------------------------------------------------!

      subroutine jacobian_column_reach_beyond_the_band(Y, f_sp_base)
      ! HOW MUCH OF AN UNKNOWN'S TRUE JACOBIAN COLUMN LIES OUTSIDE THE BAND.
      !
      ! The band holds |row - col| <= 3*nvar - 1, which is the stencil of the
      ! LOCAL operator. A hydrodynamic unknown is local to that stencil up to
      ! the column-density coupling of the radiation. An element unknown need
      ! not be: the helium mass fraction of a cell sets that cell's opacity
      ! and with it the field, hence the heating, of every cell INSIDE it, so
      ! its column can carry entries in the energy rows of cells far below,
      ! which the band cannot hold and the preconditioner therefore does not
      ! see.
      !
      ! Measured directly: one full-residual difference per unknown, no
      ! coloring, and the fraction of the column's own 2-norm outside the
      ! band. The columns are taken in the SCALED space the solve works in,
      ! A = Drow^-1 J D, because that is the matrix the preconditioner and the
      ! damped Gauss-Newton model are built from.
      !
      ! Off unless EXHALE_JAC_COLUMN_REACH=1: it costs nvar residual
      ! evaluations per probe cell.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, allocatable :: F0(:), Fp(:), Yp(:), Dsc(:), Drw(:), col(:)
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heatw, coolw
      integer, parameter :: n_probe = 5
      integer            :: jprobe(n_probe)
      real*8  :: hstep, s_in, s_out, cmax_out, cnorm
      integer :: neq, ip, k, jc, i, irow, i_out, saved_wmode
      logical :: okp
      character(len=48) :: txt, txt2
      character(len=8)  :: env
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry
      call get_environment_variable('EXHALE_JAC_COLUMN_REACH', env)
      if (trim(env) .ne. '1') return
      ! As in species_row_jacobian_action_check: a measurement of the system
      ! is not an operation on the state.
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
      neq = nvar_jac*N
      allocate(F0(neq), Fp(neq), Yp(neq), Dsc(neq), Drw(neq), col(neq))
      call cell_state_scales(Y, Dsc)
      Drw = Dsc
      if (nspec_row .gt. 0) call cell_row_scales(Y, Dsc, Drw)
      jprobe = (/ 1, max(2,N/5), max(3,2*N/5), max(4,3*N/5), max(5,4*N/5) /)
      saved_wmode = weno_mode
      weno_mode   = 1
      call eval_residual(Y, f_sp_base, fwork, F0, heatw, coolw)
      weno_mode   = 2
      call eval_residual(Y, f_sp_base, fwork, F0, heatw, coolw,           &
                         admissible=okp, may_be_adopted=.false.,          &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      hstep = sqrt(epsilon(1.0d0))
      write(*,'(A,I0,A)') ' (jac_column_reach) the fraction of a scaled'// &
           ' Jacobian column that lies outside the band, kl = ku = ',      &
           kl_jac, ':'
      do ip = 1, n_probe
         jc = jprobe(ip)
         if (jc .lt. 1 .or. jc .gt. N) cycle
         do k = 1, nvar_jac
            i  = nvar_jac*(jc-1) + k
            Yp = Y
            Yp(i) = Y(i) + hstep*Dsc(i)
            call eval_residual(Yp, f_sp_base, fwork, Fp, heatw, coolw,     &
                               admissible=okp, may_be_adopted=.false.,     &
                               state_is_discarded=.true.,                  &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            if (.not. okp) then
               call name_of_unknown(i, txt)
               write(*,'(A,A,A)') ' (jac_column_reach)   ', trim(txt),     &
                    ': the probe state is not describable, no column'
               cycle
            endif
            ! The scaled column A(:,i) = Drow^-1 (dF/dY_i) D_i.
            col = (Fp - F0)/hstep/Drw
            s_in = 0.0d0;  s_out = 0.0d0;  cmax_out = 0.0d0;  i_out = i
            do irow = 1, neq
               if (abs(irow - i) .le. kl_jac) then
                  s_in  = s_in  + col(irow)*col(irow)
               else
                  s_out = s_out + col(irow)*col(irow)
                  if (abs(col(irow)) .gt. cmax_out) then
                     cmax_out = abs(col(irow));  i_out = irow
                  endif
               endif
            enddo
            cnorm = sqrt(s_in + s_out)
            call name_of_unknown(i, txt)
            call name_of_unknown(i_out, txt2)
            if (cnorm .le. 0.0d0) then
               write(*,'(A,A,A)') ' (jac_column_reach)   ', trim(txt),     &
                    ': the column is identically zero'
            else
               write(*,'(A,A,A,F8.4,A,ES10.3,A,A)')                        &
                    ' (jac_column_reach)   ', trim(txt),                   &
                    ': outside/total ', sqrt(s_out)/cnorm,                 &
                    ' , ||col|| ', cnorm, ' , largest outside entry in ',  &
                    trim(txt2)
            endif
         enddo
      enddo
      weno_mode = saved_wmode
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      call put_back_residual_evaluation_products(products_at_entry)
      call put_back_solve_refusal_statistics(statistics_at_entry)
      deallocate(F0, Fp, Yp, Dsc, Drw, col)
      end subroutine jacobian_column_reach_beyond_the_band

      ! ------------------------------------------------------!

      subroutine stationary_rows_of_the_returned_state(rep, ok, iworst,   &
                                                       none_measured,     &
                                                       n_rows_read)
      ! THE EQUATIONS THIS SOLVER DRIVES TO ZERO, read off the certification
      ! report of the state it is asked about -- the same record the
      ! certification block printed beside it, not a second evaluation and
      ! not a second formula. `ok` is true when every row of
      ! THE SYSTEM THIS SOLVE CARRIED was evaluated on that state, is finite
      ! and is within its own tolerance. `n_rows_read`, when asked for, is
      ! how many rows that was: a caller reporting that the rows are met
      ! says over how many equations it is speaking.
      ! `iworst` points at the refusing row with the largest measure;
      ! `none_measured` says a row could not be measured on this state at
      ! all, which refuses on its own and carries no number.
      !
      ! WHICH ROWS THOSE ARE. The three hydrodynamic rows always, and one
      ! species row for every transported balance the registry carried
      ! (B5, B5b: the carriers and the elemental partitions): a solve whose
      ! unknown vector held n(H2), n(H+) or an element fraction and whose
      ! residual drove that balance to zero has not delivered a stationary
      ! state until that balance is within cert_tol_carrier_at, and reading
      ! only the hydrodynamic rows would call such a state converged on
      ! three of its four equations. A balance the system did NOT carry is
      ! measured and reported by the certification but does not decide this
      ! flag: it is an equation the solve was never given.
      !
      ! The row measure is max_j |R_kj|/residual_row_scale(k,j), which is
      ! what resid_relnorm forms too (relnorm_over_cells over the two
      ! windows, combined by the larger). So this reads the same rows the
      ! residual gate reads, and the number printed by a refusal is the
      ! number the certification block prints.
      type(cert_report), intent(in)  :: rep
      logical,           intent(out) :: ok
      integer,           intent(out) :: iworst
      logical,           intent(out) :: none_measured
      integer, optional, intent(out) :: n_rows_read
      integer :: i, isr, n_read
      real*8  :: worst
      logical :: is_row
      ok = .true.;  iworst = 0;  none_measured = .false.;  worst = -1.0d0
      n_read = 0
      do i = 1, rep%n
         is_row = (index(rep%e(i)%name, 'hydrodynamic') .eq. 1)
         do isr = 1, nspec_row
            select case (srow_kind(isr))
            case (srow_carrier)
               if (trim(rep%e(i)%name) .eq.                              &
                   'carrier balance '//trim(carrier_name(srow_idx(isr))))&
                  is_row = .true.
            case (srow_element_he)
               if (trim(rep%e(i)%name) .eq.                              &
                   'elemental transport He/H partition')                 &
                  is_row = .true.
            case default
               if (trim(rep%e(i)%name) .eq. 'elemental transport '//     &
                   trim(melem_name(srow_idx(isr))))                      &
                  is_row = .true.
            end select
         enddo
         if (.not. is_row) cycle
         n_read = n_read + 1
         if (rep%e(i)%status .eq. cert_unavailable) then
            ok = .false.;  none_measured = .true.
         else if (rep%e(i)%status .eq. cert_evaluated) then
            if (.not. rep%e(i)%within_tol .or. .not. rep%e(i)%finite) then
               ok = .false.
               if (rep%e(i)%row_max .gt. worst) then
                  worst = rep%e(i)%row_max;  iworst = i
               endif
            endif
         endif
      enddo
      if (present(n_rows_read)) n_rows_read = n_read
      end subroutine stationary_rows_of_the_returned_state

      ! ------------------------------------------------------!

      subroutine certified_rows_of_the_state(F, u, f_sp, resid_tol,       &
                                             n_no_chem_root, ok, iworst,  &
                                             rep, n_rows_read)
      ! THE ROWS AN ACCEPTED STATE MUST SATISFY, MEASURED ON THE STATE
      ! GIVEN AND NOT REPORTED.
      !
      ! The stationary certification is the whole statement of what makes a
      ! state a solution of the equations the run carries, and each of its
      ! rows is judged against its own tolerance (mass 3e-12, momentum
      ! 1e-8, energy 1e-6, a transported species row 1e-5 in the wind),
      ! while the acceptance gate reads one number against the run's
      ! "Resid tol". A STOPPING TEST BUILT ON THE GATE ALONE THEREFORE
      ! STOPS FOR SUCCESS AT STATES THE ACCEPTANCE REFUSES: MEASURED on the
      ! molecular partitioned run, three hydrodynamic solves stopped at
      ! ||R|| below 1e-8 and were refused at return with mass rows
      ! 5.95e-12, 8.00e-12 and 4.99e-12 against 3e-12
      ! (docs/solver_partition_experiment_20260911.md section 7.3).
      !
      ! So this is the same evaluation and the same formula the acceptance
      ! at return makes: certification_evaluate builds the report of the
      ! state, and stationary_rows_of_the_returned_state reads out of it
      ! the rows of the system this solve carries. Nothing is written here;
      ! the caller states in one line which row refuses.
      real*8, dimension(:),                   intent(in) :: F
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,  intent(in)  :: resid_tol
      integer, intent(in)  :: n_no_chem_root
      logical, intent(out) :: ok
      integer, intent(out) :: iworst
      type(cert_report),  intent(out) :: rep
      integer,  optional, intent(out) :: n_rows_read
      real*8, dimension(3,1-Ng:N+Ng) :: Rrows
      logical :: none_measured
      integer :: jj, kk, n_read
      Rrows = 0.0d0
      do jj = 1, N
         do kk = 1, 3
            Rrows(kk,jj) = F(nvar_jac*(jj-1)+kk)
         enddo
      enddo
      call certification_evaluate(cert_context_stationary, u, Rrows,      &
                                  f_sp, resid_tol, n_no_chem_root,        &
                                  .true., rep)
      call stationary_rows_of_the_returned_state(rep, ok, iworst,         &
                                                 none_measured, n_read)
      if (present(n_rows_read)) n_rows_read = n_read
      end subroutine certified_rows_of_the_state

      ! ------------------------------------------------------!

      subroutine write_gate_met_while_a_row_refuses(tag, rep, iworst,     &
                                                    iter)
      ! WHY AN ITERATION THAT MEETS THE ACCEPTANCE GATE GOES ON. The gate
      ! is met and a row the acceptance of this state rests on is not, so
      ! the state is not a solution of the equations the solve was given
      ! and the remaining budget is spent on it.
      character(len=*),  intent(in) :: tag
      type(cert_report), intent(in) :: rep
      integer,           intent(in) :: iworst, iter
      ! THE CELL THAT REFUSED, WITH ITS OWN MEASURE AND ITS OWN TOLERANCE.
      ! Where a row's tolerance is a function of the cell, its verdict is
      ! taken at the cell standing furthest outside its own tolerance and
      ! that need not be the cell of the largest measure (mass_row_verdict,
      ! certification.f90, which records the binding cell as jbind); a pair
      ! made of the largest measure and the binding cell's tolerance would
      ! be a ratio of two different cells. jbind is zero on a row held to
      ! one number for the whole column, and there the largest measure is
      ! the one that refuses.
      real*8  :: row_refusing
      integer :: cell_refusing
      if (iworst .gt. 0) then
         if (rep%e(iworst)%jbind .gt. 0) then
            row_refusing  = rep%e(iworst)%row_at_bind
            cell_refusing = rep%e(iworst)%jbind
         else
            row_refusing  = rep%e(iworst)%row_max
            cell_refusing = rep%e(iworst)%jworst
         endif
         write(*,'(A,A,I0,A,A,A,ES10.3,A,ES8.1,A,I0,A)') tag,             &
              ' it ', iter, ': the acceptance gate is met and the'//      &
              ' certified ', trim(rep%e(iworst)%name), ' is ',            &
              row_refusing, ', above its own ',                           &
              rep%e(iworst)%tol, ' at cell ', cell_refusing,              &
              '; the iteration goes on within the remaining budget'
      else
         write(*,'(A,A,I0,A)') tag, ' it ', iter, ': the acceptance gate'//&
              ' is met and a row the acceptance rests on could not be'//  &
              ' measured on this iterate; the iteration goes on within'// &
              ' the remaining budget'
      endif
      end subroutine write_gate_met_while_a_row_refuses

      ! ------------------------------------------------------!

      subroutine gate_equals_certification(rep, rnorm_gate, tag, agree)
      ! THE ACCEPTANCE NUMBER AND THE CERTIFIED NUMBER ARE ONE MEASUREMENT
      ! OF ONE STATE, asserted here.
      !
      ! `rnorm_gate` is what steady_gates_met is given: max over the three
      ! hydrodynamic rows of the volume-independent row measure, formed by
      ! resid_relnorm as the larger of that maximum over the layer and over
      ! the wind. The certification forms the same rows independently
      ! (certification_row_measure) from the same F, so the maximum of its
      ! evaluated hydrodynamic entries must be the same number. Since both
      ! are maxima of the same finite set of ratios, agreement is exact in
      ! IEEE arithmetic and the 1e-12 below is a margin, not a tolerance on
      ! a physical quantity.
      !
      ! A DISAGREEMENT IS A DEFECT, not a property of the state: the two
      ! expressions have drifted apart, and neither number can then be said
      ! to describe the state handed back. The caller refuses the solve.
      type(cert_report), intent(in)  :: rep
      real*8,            intent(in)  :: rnorm_gate
      character(len=*),  intent(in)  :: tag
      logical,           intent(out) :: agree
      integer :: i, n_hydro
      real*8  :: rmax, rel
      rmax = 0.0d0;  n_hydro = 0
      do i = 1, rep%n
         if (index(rep%e(i)%name, 'hydrodynamic') .ne. 1) cycle
         if (rep%e(i)%status .ne. cert_evaluated) cycle
         n_hydro = n_hydro + 1
         rmax = max(rmax, rep%e(i)%row_max)
      enddo
      if (n_hydro .eq. 0) then
         agree = .false.
         write(*,'(A,A)') tag, ' the certification of the state handed'//  &
              ' back evaluated no hydrodynamic row, so the acceptance'//   &
              ' number cannot be checked against it'
         return
      endif
      rel = abs(rmax - rnorm_gate)/max(abs(rnorm_gate), 1.0d-300)
      agree = (rel .le. 1.0d-12)
      if (.not. agree)                                                    &
         write(*,'(A,A,ES11.3,A,ES11.3,A,ES9.2)') tag,                    &
              ' the acceptance measure and the certified rows of ONE'//    &
              ' state disagree: gate ', rnorm_gate, ', certified rows ',   &
              rmax, ', relative ', rel
      end subroutine gate_equals_certification

      ! ------------------------------------------------------!

      subroutine fix_active_species_bounds(Y, g, nactive)
      ! WHICH SPECIES UNKNOWNS THIS STEP MAY NOT MOVE, and why a step control
      ! has to decide that before it builds its model rather than after.
      !
      ! An unknown sitting ON a bound of its own box, whose model gradient
      ! points OUT of the box, is at an ACTIVE constraint: the model's own
      ! answer for it is outside the admissible set, so the best admissible
      ! value is the bound itself and the step's business is with the other
      ! unknowns. Leaving it free costs the whole iteration and not just that
      ! component: every finite-difference sample along a direction with a
      ! component pointing out of the box lands on a state the code cannot
      ! describe, at every step length, so GMRES returns no direction at all
      ! and the trust region has no model to step on. MEASURED on the coupled
      ! `mol_carrier` reload with the trial projected but the model left
      ! free: 12 of 20 outer iterations ended on "no Krylov direction could
      ! be sampled" and 189 residual samples were refused for a negative
      ! species unknown.
      !
      ! THE SET IS FIXED FOR THE STEP, from Y and the model gradient alone,
      ! and that is not a convenience: a projection that depended on the sign
      ! of the direction it is applied to would make the sampled operator
      ! nonlinear in its argument, and the Arnoldi relation the predicted
      ! reduction is taken from would no longer hold.
      !
      ! THE GRADIENT TEST IS WHAT KEEPS IT FROM STICKING. Freezing every
      ! unknown that merely touches a bound would freeze one that the model
      ! wants to bring back INTO the set, for ever. The rule here is the
      ! published one for a bound-constrained Newton (Bertsekas 1982,
      ! "Projected Newton methods for optimization problems with simple
      ! constraints", SIAM J. Control Optim. 20, 221): the active set is
      ! the indices at a bound whose steepest-descent direction leaves the
      ! feasible set, so an unknown the model pushes inward stays free and
      ! moves in the same iteration.
      !
      ! g is the model gradient of the merit, A^T r0 from the banded A. It is
      ! zeroed on the active set here, so that the Cauchy leg built from it
      ! has no component the step may not take.
      real*8, dimension(nvar_jac*N), intent(in)    :: Y
      real*8, dimension(nvar_jac*N), intent(inout) :: g
      integer,                       intent(out)   :: nactive
      real*8  :: lo, hi, val, band
      integer :: j, i, is
      nactive = 0
      if (.not. allocated(species_bound_is_active))                       &
         allocate(species_bound_is_active(nvar_jac*N))
      species_bound_is_active = .false.
      if (nspec_row .le. 0) return
      ! THE WIDTH OF THE ACTIVE SET IS A LENGTH IN UNKNOWN SPACE, so it is
      ! taken on the state scales and not on the coordinates the linear
      ! algebra works in (state_column_scale; the two differ once the
      ! columns are equilibrated). Without a state scale there is no probe
      ! step to measure the width against, and that can only happen before
      ! the solve has scaled its first iterate. The SHARED rows below do not
      ! need one -- their width is a fraction of a budget -- so they are
      ! fixed whether or not the coordinate bounds could be.
      if (project_species_trial .and. allocated(species_box_lo) .and.     &
          allocated(state_column_scale)) then
      do j = 1, N
         do i = 1, nspec_row
            is  = nvar_jac*(j-1) + 3 + i
            val = Y(is)
            ! The faces are the step's own, frozen at this iterate
            ! (freeze_species_unknown_box), so the set held here and the set
            ! the trial is written onto are one set.
            lo  = species_box_lo(is)
            hi  = species_box_hi(is)
            ! WITHIN ONE FINITE-DIFFERENCE STEP OF A BOUND IS ON IT, as
            ! far as a model built from finite differences is concerned. The
            ! probe step at this unknown is the square root of machine
            ! epsilon times its own column scale, so an unknown closer to
            ! the bound than that is carried across it by the sample
            ! whatever the sample is for -- and the test that asked whether
            ! the value was exactly at the bound left it free. MEASURED on
            ! the coupled `mol_carrier` reload with the exact test: the
            ! merit fell from 79.9 to 2.52 in nine outer iterations and then
            ! iteration 10 ended on "no Krylov direction could be sampled"
            ! at cell 215, whose carrier had reached 1e-30 of its own scale
            ! without ever being zero. This is the epsilon-active set of the
            ! published rule (Bertsekas 1982, section 3).
            band = sqrt(epsilon(1.0d0))*max(state_column_scale(is), 0.0d0)
            ! -g is the steepest-descent direction of the model, so a
            ! positive gradient at a lower bound points out of the box.
            ! WITHOUT THE GRADIENT TEST (EXHALE_SPECIES_BOUND_HOLD=all)
            ! every unknown at a bound is held, whichever way the model
            ! wants to move it. It is a measurement and not the default: an
            ! unknown the model wants to bring back INTO the set is then
            ! held there for every subsequent step as well, because nothing
            ! else can move it, and what that costs is item B5j's rung (a5).
            ! AN UNKNOWN WHOSE BOX IS A POINT CANNOT MOVE AT ALL, whichever
            ! way the model wants to move it, so it is held without asking
            ! the gradient. It happens where a carrier's element has no free
            ! density in the cell: the budget is then zero and the only
            ! admissible carrier density is zero. Left free, its component
            ! of any direction would give the fraction-to-the-boundary rule
            ! no room and cost the whole Krylov cycle, not just that
            ! component.
            if ((hi .le. lo) .or.                                          &
                (val .le. lo + band .and.                                  &
                 (g(is) .gt. 0.0d0 .or. hold_every_bound)) .or.            &
                (val .ge. hi - band .and.                                  &
                 (g(is) .lt. 0.0d0 .or. hold_every_bound))) then
               species_bound_is_active(is) = .true.
               g(is) = 0.0d0
               nactive = nactive + 1
            endif
         enddo
      enddo
      endif
      ! AND THE SHARED ROWS THE STEP MOVES ALONG RATHER THAN ACROSS. They
      ! are fixed here, from the same iterate and the same model gradient,
      ! and with the coordinate bounds just held removed from them, so that
      ! the projection and the hold are one declaration about one step
      ! (fix_active_element_constraints).
      call fix_active_element_constraints(g)
      end subroutine fix_active_species_bounds

      ! ------------------------------------------------------!

      subroutine release_active_species_bounds
      ! No constraint is held outside the step that fixed it: every other
      ! caller of jv_product samples the unrestricted operator. The held
      ! coordinate bounds and the active shared rows are released together,
      ! being one declaration about one step.
      if (allocated(species_bound_is_active))                             &
         species_bound_is_active = .false.
      call release_the_active_element_constraints
      end subroutine release_active_species_bounds

      ! ------------------------------------------------------!

      subroutine hold_the_active_bounds_of(v)
      ! Zero the components of a direction that the step may not move
      ! (species_bound_is_active). Applied to every direction the model is
      ! sampled along and to both legs of the dogleg, so that the operator
      ! the Krylov cycle sees and the step the trial takes are one map.
      real*8, dimension(nvar_jac*N), intent(inout) :: v
      if (.not. allocated(species_bound_is_active)) return
      where (species_bound_is_active) v = 0.0d0
      end subroutine hold_the_active_bounds_of

      ! ------------------------------------------------------!

      double precision function largest_step_inside_the_species_box(Y, v, &
                               n_blocked) result(tmax)
      ! THE LARGEST MULTIPLE OF A DIRECTION THAT KEEPS EVERY SPECIES UNKNOWN
      ! INSIDE ITS OWN BOX -- the fraction-to-the-boundary rule.
      !
      ! WHY A FINITE-DIFFERENCE PROBE NEEDS IT. The probe step is
      ! sqrt(eps_mach)(1 + ||Y||)/||v|| times v, a length set by the whole
      ! vector; the component of v at one species unknown is set by the
      ! banded preconditioner and can be enormous relative to that unknown's
      ! own distance from its bound. MEASURED on the coupled `mol_carrier`
      ! reload with the active set held but the step length unchanged: the
      ! probe put the carrier of cell 217 at -2.4e-16 with 26 species
      ! unknowns of the state within round-off of their bounds, the sample
      ! was refused at every halving, and 12 of 20 outer iterations ended
      ! with no Krylov direction at all.
      !
      ! Halving the step until it fits is a search for this number; this
      ! computes it. It is the same rule an interior-point step uses, and it
      ! keeps the quotient a directional derivative: a shorter finite
      ! difference of the same function along the same direction.
      !
      ! ONE BLOCKED COMPONENT NO LONGER COSTS THE WHOLE DIRECTION. A
      ! component sitting ON a face with the direction pointing out of it
      ! has no room at any step length, and returning zero for the whole
      ! direction threw away the motion of every other component with it:
      ! MEASURED on `mol_diffusion`, 18 Krylov cycles truncated for that
      ! reason. Such a component is COUNTED and skipped, and the length
      ! returned is the largest one the components with room allow; the
      ! caller zeroes the counted components before it probes
      ! (jv_product), so the difference quotient it takes is a directional
      ! derivative of the residual composed with the projection onto the
      ! box -- which is the map the trial itself evaluates, the trial being
      ! projected onto the same box, and the model of a projected step is
      ! recomputed with its own product (trust_region_step). The components
      ! this reaches are exactly those the active set left free because the
      ! model wanted them the other way (fix_active_species_bounds).
      !
      ! huge(1) when no component of v points out of the box, which is every
      ! case on the three-unknown route -- there are no species unknowns and
      ! nothing here applies.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      real*8, dimension(nvar_jac*N), intent(in) :: v
      integer, optional, intent(out) :: n_blocked
      real*8  :: lo, hi, val, room
      integer :: j, i, is
      tmax = huge(1.0d0)
      if (present(n_blocked)) n_blocked = 0
      if (nspec_row .le. 0) return
      if (.not. allocated(species_box_lo)) return
      do j = 1, N
         do i = 1, nspec_row
            is  = nvar_jac*(j-1) + 3 + i
            if (v(is) .eq. 0.0d0) cycle
            val = Y(is)
            ! The same faces the trial is projected onto and the active set
            ! is taken from (freeze_species_unknown_box).
            lo  = species_box_lo(is)
            hi  = species_box_hi(is)
            if (v(is) .lt. 0.0d0) then
               if (lo .le. -huge(1.0d0)) cycle
               room = val - lo
            else
               room = hi - val
            endif
            if (room .le. 0.0d0) then
               if (present(n_blocked)) n_blocked = n_blocked + 1
               cycle
            endif
            tmax = min(tmax, room/abs(v(is)))
         enddo
      enddo
      end function largest_step_inside_the_species_box

      ! ------------------------------------------------------!

      subroutine zero_the_blocked_components_of(Y, v)
      ! The components of a direction that have no room at any step length,
      ! set to zero. It is the same test largest_step_inside_the_species_box
      ! counts, applied to the direction it was counted on, so that the
      ! probe steps along the direction whose length that function returned.
      real*8, dimension(nvar_jac*N), intent(in)    :: Y
      real*8, dimension(nvar_jac*N), intent(inout) :: v
      real*8  :: lo, hi, room
      integer :: j, i, is
      if (nspec_row .le. 0) return
      if (.not. allocated(species_box_lo)) return
      do j = 1, N
         do i = 1, nspec_row
            is = nvar_jac*(j-1) + 3 + i
            if (v(is) .eq. 0.0d0) cycle
            lo = species_box_lo(is)
            hi = species_box_hi(is)
            if (v(is) .lt. 0.0d0) then
               if (lo .le. -huge(1.0d0)) cycle
               room = Y(is) - lo
            else
               room = hi - Y(is)
            endif
            if (room .le. 0.0d0) v(is) = 0.0d0
         enddo
      enddo
      end subroutine zero_the_blocked_components_of

      ! ------------------------------------------------------!

      subroutine note_adopted_species_unknowns(Y)
      ! Record, over every state the solve adopts: the smallest species
      ! density or mass fraction the state names
      ! (species_unknown_min_adopted), and for a carrier its distance from
      ! the floor of its own unknown space as the ratio n/floor
      ! (carrier_over_its_floor_min) with the count of the unknowns at or
      ! below it (n_carrier_at_its_floor), and its distance from the element
      ! budget above it as the ratio n/budget with the count of the unknowns
      ! above THAT (carrier_over_its_budget_max,
      ! n_carrier_above_its_budget).
      !
      ! THE RATIO IS FORMED IN BOTH ARMS. With the density unknown the floor
      ! is not a bound of the space and a carrier driven to zero reports a
      ! ratio of zero, which is the whole difference between the two arms in
      ! one number; with ln n the density cannot be zero and the ratio says
      ! how far above the floor the solve stayed.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      real*8  :: val, fl, nbu
      integer :: j, i, is
      if (nspec_row .le. 0) return
      do j = 1, N
         do i = 1, nspec_row
            is  = nvar_jac*(j-1) + 3 + i
            val = Y(is)
            if (srow_kind(i) .eq. srow_carrier) then
               val = carrier_density_from_unknown(Y(is), j)
               fl  = carrier_log_floor_at(j, i)
               carrier_over_its_floor_min =                               &
                  min(carrier_over_its_floor_min, val/fl)
               if (val .le. fl)                                           &
                  n_carrier_at_its_floor = n_carrier_at_its_floor + 1
               ! AND THE SAME QUESTION AT THE OTHER FACE: no adopted state
               ! may carry more of a carrier than the element it is made of
               ! has free in that cell. The reference is the element budget
               ! itself and not the face the box uses, which is one
               ! floating-point step inside it and is raised to the iterate
               ! where the iterate stands outside, so this states the
               ! physics rather than the step control.
               if (carrier_headroom_known()) then
                  nbu = carrier_headroom(j, srow_idx(i))/n0
                  if (nbu .gt. 0.0d0) then
                     carrier_over_its_budget_max =                        &
                        max(carrier_over_its_budget_max, val/nbu)
                     if (val .gt. nbu) n_carrier_above_its_budget =       &
                        n_carrier_above_its_budget + 1
                  endif
               endif
            endif
            species_unknown_min_adopted =                                 &
               min(species_unknown_min_adopted, val)
         enddo
      enddo
      end subroutine note_adopted_species_unknowns

      ! ------------------------------------------------------!

      function species_unknown_space_text() result(txt)
      ! What the species unknowns are, and how their bounds are kept, in one
      ! line of the solve's own report.
      character(len=100) :: txt
      ! A SOLVE WITH NO CARRIER ROW SAYS NOTHING ABOUT THE CARRIER UNKNOWN.
      ! The coupled route is entered by a configuration whose only
      ! transported balance is an element as well (He_diffusion in an atomic
      ! gas), and the unknown space of that solve is the element fractions
      ! alone.
      if (.not. carrier_rows_registered()) then
         if (project_species_trial) then
            txt = 'an element is its mass fraction; trials projected'//   &
                  ' onto their bounds'
         else
            txt = 'an element is its mass fraction; steps shortened at'// &
                  ' a bound'
         endif
      else if (carrier_unknown_is_logarithmic                             &
               .and. .not. element_constraint_rows_on) then
         txt = 'a carrier is ln n, an element is its mass fraction;'//    &
               ' the element budget is a face of the box'
      else if (carrier_unknown_is_logarithmic) then
         txt = 'a carrier is ln n, an element is its mass fraction'
      else if (project_species_trial) then
         txt = 'a carrier is its density, an element its mass fraction;'// &
               ' trials projected onto their bounds'
      else
         txt = 'a carrier is its density, an element its mass fraction;'// &
               ' steps shortened at a bound'
      endif
      end function species_unknown_space_text

      ! ------------------------------------------------------!

      subroutine read_species_unknown_space_controls
      ! WHAT THE SPECIES UNKNOWNS ARE AND HOW THEIR BOUNDS ARE KEPT, read
      ! once per solve.
      !   EXHALE_SPECIES_BOUND_PROJECT=0   shorten the step against a bound
      !                                    instead of writing the trial onto
      !                                    the face of the box.
      !   EXHALE_ELEMENT_WRITE_BACK=0      write the carrier column without
      !                                    restoring the element totals of
      !                                    the cell, which is the arm the
      !                                    conservation is measured against
      !                                    (item N4b).
      !   EXHALE_ELEMENT_CONSTRAINT_ROWS=0 carry the shared element budget
      !                                    as a coordinate face of the box
      !                                    instead of as a row of the step,
      !                                    which is the box the B5j to B5l
      !                                    measurements were made in and is
      !                                    the arm the rows are measured
      !                                    against (item N4b).
      !   EXHALE_CARRIER_LOG_UNKNOWN=0     carry the carrier DENSITY as the
      !                                    unknown instead of ln n, which is
      !                                    the arm the logarithm was measured
      !                                    against (item B5j) and is not the
      !                                    default (user decision,
      !                                    2026-09-08).
      character(len=32) :: env
      integer :: iselect
      real*8  :: xselect
      project_species_trial = .true.
      call get_environment_variable('EXHALE_SPECIES_BOUND_PROJECT', env)
      if (trim(env) .eq. '0') project_species_trial = .false.
      element_constraint_rows_on = .true.
      call get_environment_variable('EXHALE_ELEMENT_CONSTRAINT_ROWS', env)
      if (trim(env) .eq. '0') element_constraint_rows_on = .false.
      element_write_back_conserves = .true.
      call get_environment_variable('EXHALE_ELEMENT_WRITE_BACK', env)
      if (trim(env) .eq. '0') element_write_back_conserves = .false.
      hold_every_bound = .false.
      call get_environment_variable('EXHALE_SPECIES_BOUND_HOLD', env)
      if (trim(env) .eq. 'all') hold_every_bound = .true.
      carrier_log_unknown_wanted = .true.
      call get_environment_variable('EXHALE_CARRIER_LOG_UNKNOWN', env)
      if (trim(env) .eq. '0') carrier_log_unknown_wanted = .false.
      !   EXHALE_CAUCHY_LEG_BY_IMAGE=1     admit the approximate-gradient
      !                                    leg by its share of the model
      !                                    image alone, the length test
      !                                    against the radius not taken
      !                                    (cauchy_leg_admitted_by_image_alone).
      cauchy_leg_admitted_by_image_alone = .false.
      call get_environment_variable('EXHALE_CAUCHY_LEG_BY_IMAGE', env)
      if (trim(env) .eq. '1') cauchy_leg_admitted_by_image_alone = .true.
      !   EXHALE_GM_ORTHO=1              measure the loss of orthogonality
      !                                  of the Arnoldi basis
      !                                  (gm_orthogonality_loss_max).
      gm_measure_orthogonality  = .false.
      call get_environment_variable('EXHALE_GM_ORTHO', env)
      if (trim(env) .eq. '1') gm_measure_orthogonality = .true.
      gm_orthogonality_loss_max   = 0.0d0
      gm_orthogonality_loss_cycle = 0.0d0
      !   EXHALE_GM_REORTHO=1            orthogonalize the Arnoldi basis
      !                                  twice (gm_reorthogonalize). The
      !                                  measured loss of orthogonality is
      !                                  what decides whether it is worth
      !                                  its inner products
      !                                  (gm_reorthogonalization_loss_level).
      gm_reorthogonalize = .false.
      call get_environment_variable('EXHALE_GM_REORTHO', env)
      if (trim(env) .eq. '1') gm_reorthogonalize = .true.
      !   EXHALE_GM_IMAGE_CHECK=1        check the image the Arnoldi
      !                                  relation returns against the
      !                                  operator, column by column, and
      !                                  sample one direction twice
      !                                  (gm_verify_the_arnoldi_image). One
      !                                  product of the operator per column.
      gm_verify_the_arnoldi_image = .false.
      call get_environment_variable('EXHALE_GM_IMAGE_CHECK', env)
      if (trim(env) .eq. '1') gm_verify_the_arnoldi_image = .true.
      !   EXHALE_GM_HISTORY=1            print the Krylov cycle's residual
      !                                  history: the reduced problem's
      !                                  residual after every product and
      !                                  the true one every twenty
      !                                  (gm_residual_history_on).
      gm_residual_history_on = .false.
      call get_environment_variable('EXHALE_GM_HISTORY', env)
      if (trim(env) .eq. '1') gm_residual_history_on = .true.
      gm_residual_history_here = .false.
      !   EXHALE_KRYLOV_SIZE_SCAN=1      repeat the Krylov cycle of the
      !                                  selected outer iteration at 40,
      !                                  80, 160 and 320 products, and the
      !                                  last two again with the basis
      !                                  orthogonalized twice, adopting
      !                                  none of the steps
      !                                  (krylov_size_scan_on).
      krylov_size_scan_on = .false.
      call get_environment_variable('EXHALE_KRYLOV_SIZE_SCAN', env)
      if (trim(env) .eq. '1') krylov_size_scan_on = .true.
      !   EXHALE_PRECOND_SPECTRUM=1      the Ritz values of the
      !                                  preconditioned operator at the
      !                                  selected outer iteration, on the
      !                                  whole operator and on its
      !                                  compressions onto the species and
      !                                  the hydrodynamic rows
      !                                  (precond_spectrum_on).
      precond_spectrum_on = .false.
      call get_environment_variable('EXHALE_PRECOND_SPECTRUM', env)
      if (trim(env) .eq. '1') precond_spectrum_on = .true.
      if (allocated(ritz_vector_of_the_smallest))                        &
         deallocate(ritz_vector_of_the_smallest)
      if (allocated(ritz_value_of_the_smallest))                         &
         deallocate(ritz_value_of_the_smallest)
      !   EXHALE_BAND_DIFFERENCE=1       the difference between the full
      !                                  Jacobian action and the banded
      !                                  one, on the directions the
      !                                  spectrum named and on the first
      !                                  two Arnoldi directions
      !                                  (band_difference_on).
      band_difference_on = .false.
      call get_environment_variable('EXHALE_BAND_DIFFERENCE', env)
      if (trim(env) .eq. '1') band_difference_on = .true.
      !   EXHALE_FRONT_ROW=1             the species row that binds at the
      !                                  selected outer iteration and its
      !                                  two neighbors, term by term, with
      !                                  the band's entries attributed to
      !                                  those terms (front_row_on).
      front_row_on = .false.
      call get_environment_variable('EXHALE_FRONT_ROW', env)
      if (trim(env) .eq. '1') front_row_on = .true.
      !   EXHALE_JV_ADDITIVITY=1         measure whether the
      !                                  finite-difference action is
      !                                  additive across directions, and
      !                                  where and why it is not
      !                                  (jv_additivity_on). Costs about
      !                                  twenty residual evaluations at
      !                                  each of three outer iterations.
      jv_additivity_on = .false.
      call get_environment_variable('EXHALE_JV_ADDITIVITY', env)
      if (trim(env) .eq. '1') jv_additivity_on = .true.
      !   EXHALE_JV_COLUMN_SCALE=1       take the probe step of the
      !                                  matrix-free action on the column
      !                                  scales, so that every unknown is
      !                                  displaced by the same fraction of
      !                                  its own magnitude
      !                                  (probe_step_on_the_column_scales).
      jv_probe_on_the_column_scales = .false.
      call get_environment_variable('EXHALE_JV_COLUMN_SCALE', env)
      if (trim(env) .eq. '1') jv_probe_on_the_column_scales = .true.
      !   EXHALE_JV_PROBE_ARC=<x>        multiply the probe arc of the
      !                                  matrix-free action by x
      !                                  (jv_probe_arc_scale). A
      !                                  non-positive or unreadable value
      !                                  leaves it at one.
      jv_probe_arc_scale = 1.0d0
      call get_environment_variable('EXHALE_JV_PROBE_ARC', env)
      if (len_trim(env) .gt. 0) then
         xselect = 0.0d0
         read(env,*,err=193,end=193) xselect
  193    continue
         if (xselect .gt. 0.0d0) jv_probe_arc_scale = xselect
      endif
      !   EXHALE_RESID_JUMP_SCAN=1       sample the residual densely along
      !                                  one direction and say where it
      !                                  JUMPS, in which term of the row
      !                                  and in which operator of the flux
      !                                  pipeline (resid_jump_scan_on).
      !                                  About seventy residual
      !                                  evaluations at one outer
      !                                  iteration.
      resid_jump_scan_on = .false.
      call get_environment_variable('EXHALE_RESID_JUMP_SCAN', env)
      if (trim(env) .eq. '1') resid_jump_scan_on = .true.
      resid_jump_scan_here   = .false.
      jv_additivity_here     = .false.
      jv_probe_length_factor = 1.0d0
      !   EXHALE_KRYLOV_ON_THE_BALL=1    stop the Krylov cycle where the
      !                                  iterate leaves the trust ball and
      !                                  return the point on its boundary
      !                                  (krylov_truncated_on_the_trust_ball).
      krylov_truncated_on_the_trust_ball = .false.
      call get_environment_variable('EXHALE_KRYLOV_ON_THE_BALL', env)
      if (trim(env) .eq. '1') krylov_truncated_on_the_trust_ball = .true.
      !   EXHALE_MODEL_ROW_EQUIL=1       equilibrate the rows of the
      !                                  LINEAR MODEL to unit infinity
      !                                  norm for the Krylov solve alone,
      !                                  the merit, the gate, the trust
      !                                  region and the certification
      !                                  keeping the certification row
      !                                  scales (model_row_equilibration_on).
      model_row_equilibration_on = .false.
      call get_environment_variable('EXHALE_MODEL_ROW_EQUIL', env)
      if (trim(env) .eq. '1') model_row_equilibration_on = .true.
      model_rows_equilibrated = .false.
      !   EXHALE_GM_TRUE_RESIDUAL=1      let the Krylov cycle return the
      !                                  step its TRUE residual chooses
      !                                  and stop where that residual
      !                                  turns upward
      !                                  (gm_step_by_its_true_residual).
      !                                  One product of the operator every
      !                                  gm_true_residual_stride products
      !                                  past gm_true_residual_first.
      gm_step_by_its_true_residual = .false.
      call get_environment_variable('EXHALE_GM_TRUE_RESIDUAL', env)
      if (trim(env) .eq. '1') gm_step_by_its_true_residual = .true.
      !   EXHALE_ELEM_DIAG=1             the ordered element-row diagnostic
      !                                  at outer iterations 1, 2 and 60;
      !                                  an integer above 1 replaces 60.
      elem_diag_on    = .false.
      elem_diag_third = 60
      elem_diag_here  = .false.
      call get_environment_variable('EXHALE_ELEM_DIAG', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=91,end=91) iselect
   91    continue
         elem_diag_on = (iselect .ge. 1)
         if (iselect .gt. 1) elem_diag_third = iselect
      endif
      !   EXHALE_ELEM_DIAG_FROM=<j>      the first cell the element-row term
      !                                  diagnostic prints, for both the
      !                                  trace hook and the helium row.
      !                                  Absent: the outermost six cells.
      elem_diag_cell_from = 0
      call get_environment_variable('EXHALE_ELEM_DIAG_FROM', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=96,end=96) iselect
   96    continue
         elem_diag_cell_from = max(0, iselect)
      endif
      !   EXHALE_COLUMN_EQUIL=0          leave the column scaling at the
      !                                  state scales, so that the scaled
      !                                  operator is the one the row
      !                                  scaling alone produces.
      !   EXHALE_COLUMN_EQUIL=1          equilibrate the columns of the
      !                                  MODEL as well, which moves the
      !                                  trust-region coordinates with them
      !                                  (measured, off by default).
      !   EXHALE_PRECON_EQUIL=0          factor the band as the scalings
      !                                  leave it, with no equilibration of
      !                                  the preconditioner.
      column_equilibration_on = .false.
      call get_environment_variable('EXHALE_COLUMN_EQUIL', env)
      if (trim(env) .eq. '1') column_equilibration_on = .true.
      preconditioner_equilibration_on = .true.
      call get_environment_variable('EXHALE_PRECON_EQUIL', env)
      if (trim(env) .eq. '0') preconditioner_equilibration_on = .false.
      preconditioner_equilibrated = .false.
      !   EXHALE_JFNK_MAXIT=<n>          take at most n outer iterations of
      !                                  the stationary solve, in place of
      !                                  the number the caller passed.
      jfnk_maxit_wanted = 0
      call get_environment_variable('EXHALE_JFNK_MAXIT', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=92,end=92) iselect
   92    continue
         if (iselect .gt. 0) jfnk_maxit_wanted = iselect
      endif
      ! NOT ARMED HERE: the opening state is packed and evaluated as a
      ! density, and arm_carrier_log_unknown rewrites it once the floor the
      ! logarithm needs exists.
      carrier_unknown_is_logarithmic = .false.
      if (allocated(carrier_log_floor)) deallocate(carrier_log_floor)
      ! NOR IS THE BOX FORMED HERE, for the same reason: two of its faces
      ! are the element budget and the cell density of the iterate, and the
      ! solve has not evaluated one yet. Dropped so that no reader of it can
      ! see the box of the previous solve.
      if (allocated(species_box_lo)) deallocate(species_box_lo)
      if (allocated(species_box_hi)) deallocate(species_box_hi)
      n_projected_trials = 0
      n_active_bounds_max = 0
      n_cauchy_only_steps = 0
      n_cauchy_leg_short  = 0
      n_cauchy_leg_all_image = 0
      species_unknown_min_adopted = huge(1.0d0)
      carrier_over_its_floor_min = huge(1.0d0)
      n_carrier_at_its_floor = 0
      carrier_over_its_budget_max = 0.0d0
      n_carrier_above_its_budget = 0
      call release_active_species_bounds
      end subroutine read_species_unknown_space_controls

      ! ------------------------------------------------------!

      integer function jfnk_outer_iteration_cap(maxit_of_the_caller)     &
                       result(ncap)
      ! HOW MANY OUTER ITERATIONS THIS SOLVE MAY TAKE: the caller's number,
      ! or the one EXHALE_JFNK_MAXIT named if it named one
      ! (read_species_unknown_space_controls reads it once per solve).
      integer, intent(in) :: maxit_of_the_caller
      ncap = maxit_of_the_caller
      if (jfnk_maxit_wanted .gt. 0) ncap = jfnk_maxit_wanted
      end function jfnk_outer_iteration_cap

      ! ------------------------------------------------------!

      subroutine read_trust_region_restart_controls
      ! WHETHER THIS SOLVE RESTARTS ITS TRUST-REGION STATE, WHEN, AND WHICH
      ! PART OF IT, read once per solve.
      !
      !   EXHALE_TR_RESTART_AT=<n>     force one restart after outer
      !                                iteration n, whatever the trigger
      !                                says. The measurement arm: n is the
      !                                iteration a capped solve hands its
      !                                state back at, so the restart acts on
      !                                the same iterate a re-entry would.
      !   EXHALE_TR_RESTART_WHAT=<s>   which of the five resets a restart
      !                                takes, as letters: r the radius and
      !                                its ceiling, t the pseudo-transient
      !                                shift, c the closure map at the base,
      !                                a the acceptance memory, b the return
      !                                to the best iterate. The default is
      !                                all five, which is what entering a
      !                                solve afresh does.
      !   EXHALE_TR_RESTART_STALL=<n>  restart whenever the best iterate has
      !                                not improved for n consecutive outer
      !                                iterations. 0 disarms the trigger.
      !   EXHALE_GM_CYCLES=<n>         let the Krylov solve spend n cycles
      !                                of its subspace instead of one, the
      !                                subspace size unchanged.
      character(len=32) :: env
      integer :: iselect, k
      tr_restart_forced_at        = 0
      tr_restart_without_improvement = tr_restart_stall_default
      n_tr_restart                = 0
      tr_reset_wanted             = .true.
      gm_restart_cycles           = 1
      call get_environment_variable('EXHALE_TR_RESTART_AT', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=81,end=81) iselect
   81    continue
         if (iselect .gt. 0) tr_restart_forced_at = iselect
      endif
      call get_environment_variable('EXHALE_TR_RESTART_STALL', env)
      if (len_trim(env) .gt. 0) then
         iselect = -1
         read(env,*,err=82,end=82) iselect
   82    continue
         if (iselect .ge. 0) tr_restart_without_improvement = iselect
      endif
      call get_environment_variable('EXHALE_TR_RESTART_WHAT', env)
      if (len_trim(env) .gt. 0) then
         tr_reset_wanted = .false.
         do k = 1, len_trim(env)
            select case (env(k:k))
            case ('r');  tr_reset_wanted(tr_reset_radius_and_ceiling) = .true.
            case ('t');  tr_reset_wanted(tr_reset_pseudo_transient)   = .true.
            case ('c');  tr_reset_wanted(tr_reset_closure_map)        = .true.
            case ('a');  tr_reset_wanted(tr_reset_acceptance_memory)  = .true.
            case ('b');  tr_reset_wanted(tr_reset_back_to_best)       = .true.
            end select
         enddo
      endif
      gm_restart_cycles_wanted = 1
      gm_restart_cycles_from   = 1
      call get_environment_variable('EXHALE_GM_CYCLES', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=83,end=83) iselect
   83    continue
         if (iselect .gt. 0) gm_restart_cycles_wanted = iselect
      endif
      call get_environment_variable('EXHALE_GM_CYCLES_FROM', env)
      if (len_trim(env) .gt. 0) then
         iselect = 0
         read(env,*,err=84,end=84) iselect
   84    continue
         if (iselect .gt. 0) gm_restart_cycles_from = iselect
      endif
      gm_restart_cycles = gm_restart_cycles_wanted
      if (gm_restart_cycles_from .gt. 1) gm_restart_cycles = 1
      end subroutine read_trust_region_restart_controls

      ! ------------------------------------------------------!

      subroutine tr_reset_name(kind, txt)
      ! The sentence that belongs to one of the five things a restart resets.
      integer,          intent(in)  :: kind
      character(len=*), intent(out) :: txt
      select case (kind)
      case (tr_reset_radius_and_ceiling)
         txt = 'the radius and its ceiling'
      case (tr_reset_pseudo_transient)
         txt = 'the pseudo-transient shift of the preconditioner'
      case (tr_reset_closure_map)
         txt = 'the closure map at the base'
      case (tr_reset_acceptance_memory)
         txt = 'the acceptance memory'
      case (tr_reset_back_to_best)
         txt = 'the iterate: back to the best '
      case default
         txt = 'nothing named'
      end select
      end subroutine tr_reset_name

      ! ------------------------------------------------------!

      logical function trust_region_restart_is_due(iter,                  &
                              n_without_improvement) result(due)
      ! IS A RESTART OF THE TRUST-REGION STATE DUE AT THE END OF THIS OUTER
      ! ITERATION?
      !
      ! Two ways for it to be due, and they are read in this order.
      !
      ! The forced iteration is the measurement arm: it names one iterate,
      ! the one a solve capped there hands back, so that a restart on that
      ! state can be compared with continuing from it.
      !
      ! The trigger is the stagnation detector's own statistic taken
      ! earlier. That detector stops the solve when the best iterate has not
      ! improved for 2*n_stall_best consecutive iterations, and goes back to
      ! the best iterate with a monotone acceptance at n_stall_best; the
      ! same count read at tr_restart_without_improvement is a precursor of
      ! both. It is a count of iterations without improvement and NOT a
      ! count of refusals: a solve whose steps are all accepted and all
      ! sideways is exactly the circling the detector was written for
      ! (n_stall_best), and a run of refusals that still lowers the best
      ! iterate is progress.
      integer, intent(in) :: iter, n_without_improvement
      due = .false.
      if (tr_restart_forced_at .gt. 0 .and.                               &
          iter .eq. tr_restart_forced_at) then
         due = .true.
         return
      endif
      if (tr_restart_without_improvement .le. 0) return
      due = (n_without_improvement .ge. tr_restart_without_improvement)
      end function trust_region_restart_is_due

      ! ------------------------------------------------------!

      subroutine restart_the_trust_region_state(tr_delta, tr_dmax, dtau,  &
                              dtau0, n_since_best, n_no_descent,          &
                              closure_map_due, acceptance_window_due)
      ! RESET THE STATE THIS SOLVE CARRIES BESIDE ITS ITERATE, and nothing
      ! else. The iterate is not an argument: a restart is a statement about
      ! the step control and not about the state, so the composition, the
      ! conserved variables and the residual are untouched here. The one
      ! reset that does move the iterate -- the return to the best one seen
      ! -- is taken by the caller, which holds that ledger.
      !
      ! A negative radius and a negative ceiling are the absence of a
      ! radius, which is the condition solve_steady_jfnk re-initializes on
      ! (initial_trust_region_radius): the next iteration measures a fresh
      ! radius on a direction the step may take, as the first iteration of a
      ! solve does.
      !
      ! TWO OF THE FIVE ARE MEASURED ON THE NEXT ITERATE AND ARE THEREFORE
      ! ARMED HERE RATHER THAN TAKEN. The closure map at the base is
      ! measured on the row scales of the iterate it describes, and the
      ! Grippo window is collapsed onto the merit and the judged distance of
      ! the state the acceptance will compare against; both of those exist
      ! only after that iterate's own evaluation, and a value carried over
      ! from the state just left would be a window about a different point.
      real*8,  intent(inout) :: tr_delta, tr_dmax, dtau
      real*8,  intent(in)    :: dtau0
      integer, intent(inout) :: n_since_best, n_no_descent
      logical, intent(inout) :: closure_map_due, acceptance_window_due
      if (tr_reset_wanted(tr_reset_radius_and_ceiling)) then
         tr_delta = -1.0d0
         tr_dmax  = -1.0d0
      endif
      if (tr_reset_wanted(tr_reset_pseudo_transient)) dtau = dtau0
      if (tr_reset_wanted(tr_reset_closure_map)) closure_map_due = .true.
      if (tr_reset_wanted(tr_reset_acceptance_memory)) then
         acceptance_window_due = .true.
         n_since_best          = 0
         n_no_descent          = 0
      endif
      n_tr_restart = n_tr_restart + 1
      end subroutine restart_the_trust_region_state


      ! ------------------------------------------------------!

      subroutine increment_of_what_the_residual_reads(npart, npart_prev,  &
                    xh2, xh2_prev, heat, cool, heat_prev, cool_prev,      &
                    have_previous_radiation, dmax, jmax, kmax)
      ! HOW MUCH ONE PASS OF THE COMPOSITION ELIMINATION MOVED THE RESIDUAL,
      ! measured on the quantities the residual reads the composition for
      ! and on no others.
      !
      ! WITH THE STATE u FIXED -- which it is throughout the elimination --
      ! the composition reaches assemble_residual through exactly three
      ! channels, and this is the whole list (read from assemble_residual,
      ! steady_residual.f90, and caloric_state_from_composition,
      ! caloric_eos.f90):
      !
      !   n_part = n_tot + n_e   the thermal equation of state, T = p/n_part,
      !                          and the caloric map's n_part/rho;
      !   x_H2 = n(H2)/n_part    the share of the particle count whose
      !                          internal energy follows the H2
      !                          rovibrational ladder instead of 3/2 kT;
      !   heat, cool             the two radiative terms the energy row is
      !                          the difference of.
      !
      ! So a pass that leaves all three where they were leaves the residual
      ! where it was, whatever it did to the individual fractions, and a pass
      ! that moves any of them moves the residual.
      !
      ! HOW EACH CHANNEL IS NORMALIZED, and each by the quantity it enters
      ! the residual through:
      !
      !   n_part   relative to itself. It is positive definite, so no floor
      !            is needed, and T is inversely proportional to it, so the
      !            relative change of n_part IS the relative change of the
      !            temperature the rows are assembled at.
      !   x_H2     absolutely. It is a weight in [0,1] multiplying an O(1)
      !            correction to the internal energy per particle
      !            (internal_energy_of_mixture), so its absolute change is
      !            the relative change of the caloric map up to that factor.
      !            A relative measure would report the round-off of an x_H2
      !            of 1e-30, which changes no energy.
      !   heat-cool  relative to the larger of the two terms it is the
      !            difference of, which is the pair energy_row_scale divides
      !            the energy row by. Since that scale also holds the flux
      !            terms, this quotient is an UPPER bound on the change of
      !            the energy row in the units the gate reads it in: below
      !            the tolerance here, the row moved by less than the
      !            tolerance of its own scale.
      !
      ! THE FLOOR UNDER THE RADIATIVE PAIR is 1e-20 of the largest radiative
      ! term anywhere in the column, the same construction and the same
      ! number the carrier row's absolute floor uses (carrier_residual,
      ! diffusive_photochemistry.f90): a heating rate that small cannot
      ! change any observable, so a cell below it is converged whatever its
      ! own two terms do relative to each other. Without it a shielded cell
      ! whose heating and cooling are both 1e-40 of the column's would
      ! report an O(1) change for ever.
      !
      ! THE FIRST PASS HAS NO PREVIOUS RADIATIVE FIELD -- heat and cool come
      ! from the sweep, so there is none until a sweep has run -- and its
      ! measure is the two equation-of-state channels alone. The second pass
      ! supplies the first radiative comparison.
      real*8, dimension(1:N), intent(in)  :: npart, npart_prev
      real*8, dimension(1:N), intent(in)  :: xh2, xh2_prev
      real*8, dimension(1:N), intent(in)  :: heat, cool
      real*8, dimension(1:N), intent(in)  :: heat_prev, cool_prev
      logical,                intent(in)  :: have_previous_radiation
      real*8,                 intent(out) :: dmax
      integer,                intent(out) :: jmax, kmax
      real*8  :: d, rad_col, rad_scale
      integer :: j
      dmax = 0.0d0;  jmax = 0;  kmax = eq_channel_particle_count
      rad_col = 0.0d0
      if (have_previous_radiation) then
         do j = 1, N
            rad_col = max(rad_col, abs(heat(j)), abs(cool(j)))
         enddo
      endif
      do j = 1, N
         if (npart(j) .gt. 0.0d0) then
            d = abs(npart(j) - npart_prev(j))/npart(j)
            if (d .gt. dmax) then
               dmax = d;  jmax = j;  kmax = eq_channel_particle_count
            endif
         endif
         d = abs(xh2(j) - xh2_prev(j))
         if (d .gt. dmax) then
            dmax = d;  jmax = j;  kmax = eq_channel_h2_caloric
         endif
         if (have_previous_radiation) then
            rad_scale = max(abs(heat(j)), abs(cool(j)),                   &
                            abs(heat_prev(j)), abs(cool_prev(j)),         &
                            1.0d-20*rad_col)
            if (rad_scale .gt. 0.0d0) then
               d = abs((heat(j) - cool(j))                                &
                       - (heat_prev(j) - cool_prev(j)))/rad_scale
               if (d .gt. dmax) then
                  dmax = d;  jmax = j;  kmax = eq_channel_radiative
               endif
            endif
         endif
      enddo
      end subroutine increment_of_what_the_residual_reads

      ! ------------------------------------------------------!

      subroutine eq_channel_name(k, txt)
      ! The channel one increment measure was set by.
      integer,          intent(in)  :: k
      character(len=*), intent(out) :: txt
      select case (k)
      case (eq_channel_particle_count);  txt = 'particle count n_tot+n_e'
      case (eq_channel_h2_caloric);      txt = 'H2 share of the particle count'
      case (eq_channel_radiative);       txt = 'heat - cool'
      case default;                      txt = 'unnamed'
      end select
      end subroutine eq_channel_name

      ! ------------------------------------------------------!

      subroutine read_composition_elimination_controls
      ! The two numbers that govern the composition elimination inside a
      ! residual evaluation, read once per solve so that no evaluation pays
      ! for an environment lookup.
      !   EXHALE_RESID_SC_MAX   passes allowed at fixed Y; 1 is the lagged
      !                         residual of section 155, kept so that the
      !                         ladder of that section can be re-measured
      !                         without a rebuild.
      !   EXHALE_RESID_EQ_TOL   change of what the residual reads below which
      !                         the elimination has reached its fixed point.
      !   EXHALE_EQ_MEASURE     `species` restores the species-fraction
      !                         increment this replaced, for measurement.
      character(len=32) :: env
      integer, parameter :: ncell_row_max = 8
      integer :: dcells(ncell_row_max), ndcell
      n_eq_sweeps_cap = -1
      call get_environment_variable('EXHALE_RESID_SC_MAX', env)
      if (len_trim(env) .gt. 0) read(env,*) n_eq_sweeps_cap
      outer_residual_closure_target = 1.0d-8
      call get_environment_variable('EXHALE_RESID_EQ_TOL', env)
      if (len_trim(env) .gt. 0) read(env,*) outer_residual_closure_target
      ! Until the amplification has been measured the sweep is asked for the
      ! target itself, which is what it was asked for before the map was
      ! measured at all.
      closure_amplification_at_the_base = 1.0d0
      eq_sweep_reltol = outer_residual_closure_target
      call get_environment_variable('EXHALE_EQ_MEASURE', env)
      eq_measure_is_species_fraction = (trim(env) .eq. 'species')
      call get_environment_variable('EXHALE_EQ_MEASURE_TRACE', env)
      eq_measure_trace = (trim(env) .eq. '1')
      call get_environment_variable('EXHALE_XUV_SELF_PASSES', env)
      if (len_trim(env) .gt. 0) read(env,*) xuv_self_field_passes
      call get_environment_variable('EXHALE_IEQ_REPORT_CELL', env)
      if (len_trim(env) .gt. 0) read(env,*) ieq_report_cell
      ! The cells the chemical-decay diagnostic measures at are armed here,
      ! before any sweep of this solve runs, because the reaction Jacobian
      ! of a cell can be formed only where that cell's rate coefficients are
      ! live. Nothing is armed unless the diagnostic is asked for, so an
      ! ordinary solve pays no eigenvalue problem.
      ndcell = 0
      call get_environment_variable('EXHALE_RESID_ROW_CONDITION', env)
      if (trim(env) .eq. '1') then
         call cells_of_the_row_report(dcells, ndcell, ncell_row_max)
         call set_molecular_decay_cells(dcells, ndcell)
      endif
      n_eq_sweeps_tot   = 0
      n_eq_sweeps_model = 1
      end subroutine read_composition_elimination_controls

      ! ------------------------------------------------------!

      subroutine trace_composition_elimination(Y, f_sp, u)
      ! WHERE THE ELIMINATION FAILS TO REACH ONE COMPOSITION, pass by pass.
      ! R(Y) = L(Y) + S(Y, c*(Y)) is a function of Y only if the Picard
      ! iteration that eliminates the composition reaches the same c*(Y)
      ! from any admissible seed. This traces that iteration from TWO seeds
      ! in lockstep, one pass at a time, and reports for each pass how far
      ! the composition moved and where, the acceptance census of the sweep
      ! (a cell accepted under the relaxation amnesty or above its cap has no
      ! root at this state and keeps whatever it was seeded with), the two
      ! residuals, and the distance between them in the acceptance measure.
      ! A stall at a floor, a limit cycle between two states and a slow
      ! approach are three different pass histories and cannot be told apart
      ! from the endpoint alone.
      !
      ! A chain of single-pass evaluations is the same map as one multi-pass
      ! evaluation: a pass ends with the composition refresh and Apply_BC
      ! that the next evaluation's entry repeats, and both are functions of
      ! (Y, f_sp) alone.
      !
      ! Enabled by EXHALE_RESID_SC_TRACE=<passes>, which STOPS the run when
      ! the trace is written.
      real*8, dimension(nvar_jac*N),           intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species) :: comp_a, comp_b, comp_a_prev, comp_b_prev
      real*8, dimension(nvar_jac*N) :: Fvec_a, Fvec_b, Fvec_d
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool
      real*8, dimension(1-Ng:N+Ng)  :: heat_a_prof, cool_a_prof
      real*8  :: rna, rnb, rnd, rc(3), dfa, dfb, dab, seed_scale
      integer :: k, npass, ja, sa, jb, sb, jab, sab, acc4a, acc6a, acc4b, acc6b
      character(len=16) :: env

      call get_environment_variable('EXHALE_RESID_SC_TRACE', env)
      if (len_trim(env) .eq. 0) return
      npass = 0
      read(env,*) npass
      if (npass .le. 0) return

      comp_a = f_sp
      comp_b = f_sp
      ! Seed B is either the column's own composition shifted by five cells
      ! (a seed no cell is the equilibrium of, but one the network can
      ! describe) or the same composition with the molecular hydrogen
      ! fraction scaled, which is the small-perturbation direction the basin
      ! scan uses.
      call get_environment_variable('EXHALE_RESID_SC_TRACE_SCALE', env)
      if (len_trim(env) .gt. 0) then
         read(env,*) seed_scale
         comp_b(1:N,isp_H2) = f_sp(1:N,isp_H2)*(1.0d0 + seed_scale)
      else
         do k = 1, N
            comp_b(k,:) = f_sp(min(N, k+5),:)
         enddo
      endif

      write(*,'(A)') ' (sc_trace) pass | seed A: d(comp) cell sp  amnesty'// &
           '(4/6) | seed B: d(comp) cell sp  amnesty(4/6) |'//               &
           ' ||R||_A  ||R||_B  ||R_B-R_A||  max|fB-fA| cell sp'
      do k = 1, npass
         comp_a_prev = comp_a;  comp_b_prev = comp_b
         call eval_residual(Y, comp_a_prev, comp_a, Fvec_a, heat, cool,                  &
                            n_eq_sweeps_fixed = 1, state_is_discarded=.true.)
         acc4a = ieq_sweep_ledger_last%acc_n(4)
         acc6a = ieq_sweep_ledger_last%acc_n(6)
         call resid_relnorm(Fvec_a, u, rc, rna)
         call worst_composition_change(comp_a, comp_a_prev, dfa, ja, sa)
         call eval_residual(Y, comp_b_prev, comp_b, Fvec_b, heat, cool,                  &
                            n_eq_sweeps_fixed = 1, state_is_discarded=.true.)
         acc4b = ieq_sweep_ledger_last%acc_n(4)
         acc6b = ieq_sweep_ledger_last%acc_n(6)
         call resid_relnorm(Fvec_b, u, rc, rnb)
         call worst_composition_change(comp_b, comp_b_prev, dfb, jb, sb)
         Fvec_d = Fvec_b - Fvec_a
         call resid_relnorm(Fvec_d, u, rc, rnd)
         call worst_composition_change(comp_b, comp_a, dab, jab, sab)
         write(*,'(A,I4,2(ES10.2,I5,I4,I6,A,I0),3ES11.3,ES10.2,I5,I4)')      &
              ' (sc_trace) ', k,                                             &
              dfa, ja, sa, acc4a, '/', acc6a,                                &
              dfb, jb, sb, acc4b, '/', acc6b,                                &
              rna, rnb, rnd, dab, jab, sab
      enddo

      ! WHERE the two fixed points differ, cell by cell: the profile is what
      ! says whether the difference is one front cell or a whole region,
      ! which species carries it and which rows of the residual it reaches.
      open(unit=77, file='sc_trace_profile.txt', status='replace')
      write(77,'(A)') '# columns: j r  f_A(1:n_species)  f_B(1:n_species)'//&
           '  heat_A cool_A heat_B cool_B  Rmass_A Rmom_A Rene_A'//         &
           '  Rmass_B Rmom_B Rene_B'
      call eval_residual(Y, comp_a, comp_a_prev, Fvec_a, heat, cool,        &
                         n_eq_sweeps_fixed = 1, state_is_discarded=.true.)
      heat_a_prof = heat;  cool_a_prof = cool
      call eval_residual(Y, comp_b, comp_b_prev, Fvec_b, heat, cool,        &
                         n_eq_sweeps_fixed = 1, state_is_discarded=.true.)
      do k = 1, N
         write(77,'(I5,ES16.8)', advance='no') k, r(k)
         write(77,'(1000ES23.15)', advance='no') comp_a_prev(k,:)
         write(77,'(1000ES23.15)', advance='no') comp_b_prev(k,:)
         write(77,'(10ES23.15)') heat_a_prof(k), cool_a_prof(k),            &
            heat(k), cool(k),                                              &
            Fvec_a(nvar_jac*(k-1)+1), Fvec_a(nvar_jac*(k-1)+2),             &
            Fvec_a(nvar_jac*(k-1)+3),                                      &
            Fvec_b(nvar_jac*(k-1)+1), Fvec_b(nvar_jac*(k-1)+2),            &
            Fvec_b(nvar_jac*(k-1)+3)
      enddo
      close(77)
      write(*,'(A)') ' (sc_trace) profile of the two fixed points written'//&
           ' to sc_trace_profile.txt'
      stop
      end subroutine trace_composition_elimination

      ! ------------------------------------------------------!

      subroutine scan_composition_elimination_basin(Y, f_sp, u)
      ! HOW FAR THE SEED HAS TO MOVE for the elimination to reach a
      ! DIFFERENT composition. The elimination of this state has more than
      ! one fixed point (trace_composition_elimination above), so c*(Y) is
      ! defined only once a rule says which one; whether that matters to a
      ! Newton depends on where the basin boundary sits relative to the
      ! seeds the solve actually hands over. Each seed here is the state's
      ! own composition with the molecular hydrogen fraction scaled by
      ! (1 + delta), the elimination is run to its fixed point, and the
      ! residual is compared with the one the unscaled seed reaches.
      !
      ! Enabled by EXHALE_RESID_SC_BASIN=<passes>, which STOPS the run.
      real*8, dimension(nvar_jac*N),           intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species) :: comp_0, comp_d, comp_prev
      real*8, dimension(nvar_jac*N) :: Fvec_0, Fvec_d, Fvec_diff
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool
      real*8  :: rn0, rnd, rc(3), delta, h2_214
      integer :: k, id, npass, jab, sab
      real*8  :: dab
      real*8, parameter :: dlist(9) = (/ 1.0d-8, 1.0d-6, 1.0d-4, 1.0d-3,   &
                                         1.0d-2, 1.0d-1, -1.0d-2,          &
                                        -1.0d-1, -5.0d-1 /)
      character(len=16) :: env

      call get_environment_variable('EXHALE_RESID_SC_BASIN', env)
      if (len_trim(env) .eq. 0) return
      npass = 0
      read(env,*) npass
      if (npass .le. 0) return

      comp_0 = f_sp
      do k = 1, npass
         comp_prev = comp_0
         call eval_residual(Y, comp_prev, comp_0, Fvec_0, heat, cool,      &
                            n_eq_sweeps_fixed = 1, state_is_discarded=.true.)
      enddo
      call resid_relnorm(Fvec_0, u, rc, rn0)
      write(*,'(A,ES11.3,A,ES12.5)') ' (sc_basin) unscaled seed: ||R|| ',  &
           rn0, ', x(H2) of cell 214 ', comp_0(214,isp_H2)

      do id = 1, size(dlist)
         delta  = dlist(id)
         comp_d = f_sp
         comp_d(1:N,isp_H2) = f_sp(1:N,isp_H2)*(1.0d0 + delta)
         do k = 1, npass
            comp_prev = comp_d
            call eval_residual(Y, comp_prev, comp_d, Fvec_d, heat, cool,   &
                               n_eq_sweeps_fixed = 1,                      &
                               state_is_discarded=.true.)
         enddo
         Fvec_diff = Fvec_d - Fvec_0
         call resid_relnorm(Fvec_diff, u, rc, rnd)
         call worst_composition_change(comp_d, comp_0, dab, jab, sab)
         h2_214 = comp_d(214,isp_H2)
         write(*,'(A,ES10.2,A,ES11.3,A,ES12.5,A,ES10.2,I5,I4)')            &
              ' (sc_basin) seed x(H2)*(1+', delta, '): ||R_d-R_0|| ',      &
              rnd, ', x(H2) of cell 214 ', h2_214,                         &
              ', max|f_d-f_0| ', dab, jab, sab
      enddo
      stop
      end subroutine scan_composition_elimination_basin

      ! ------------------------------------------------------!

      subroutine worst_composition_change(f1, f2, dmax, jmax, smax)
      ! The largest relative difference between two compositions over the
      ! physical cells, and where it sits. Same measure and same floor as
      ! the elimination's own stopping test, so the two numbers compare.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f1, f2
      real*8,                                 intent(out) :: dmax
      integer,                                intent(out) :: jmax, smax
      real*8  :: d
      integer :: j, i
      dmax = 0.0d0;  jmax = 0;  smax = 0
      do i = 1, n_species
      do j = 1, N
         d = abs(f1(j,i) - f2(j,i))/max(abs(f1(j,i)), eq_sweep_floor)
         if (d .gt. dmax) then
            dmax = d;  jmax = j;  smax = i
         endif
      enddo
      enddo
      end subroutine worst_composition_change

      ! ------------------------------------------------------!

      subroutine probe_eliminated_composition_is_seed_independent(Y, f_sp, u)
      ! THE INVARIANT THIS ITEM RESTS ON: the composition is an ELIMINATED
      ! variable, so the residual of a state does not depend on the
      ! composition the caller happened to hand the elimination. Evaluate F
      ! at one Y from the run's own composition, then from a composition that
      ! is NOT the equilibrium of that state, and the two must be the same
      ! residual.
      !
      ! THE DISPLACED SEED is the column's own composition shifted by five
      ! cells, so every cell is given the composition of a neighbor a few
      ! scale-heights away: positive, normalized, and not the equilibrium of
      ! the cell it is given to. A seed invented from arithmetic could be
      ! outside the network's domain, and then the probe would be measuring
      ! the sweep's failure and not its convergence.
      !
      ! WHAT IT MEASURES: the two residuals differenced and passed through
      ! the acceptance measure itself (resid_relnorm), so the number is in
      ! the units `Resid tol` is written in and can be compared with it
      ! directly. With the elimination lagged (EXHALE_RESID_SC_MAX=1) this is
      ! the size of the lag; with the elimination converged it is the
      ! accuracy of the fixed point.
      !
      ! Enabled by EXHALE_RESID_SC_PROBE=1, which STOPS the run once the
      ! two evaluations are reported: the probe is about one residual
      ! evaluation and there is nothing after it to measure.
      real*8, dimension(nvar_jac*N),           intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_a, f_b, f_seed
      real*8, dimension(nvar_jac*N) :: Fa, Fb, Fd
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool
      real*8  :: rna, rnb, rnd, rc(3), dseed
      integer :: na, nb, j, jsrc
      integer, parameter :: nshift = 5
      character(len=8) :: env

      call get_environment_variable('EXHALE_RESID_SC_PROBE', env)
      if (trim(env) .ne. '1') return

      call eval_residual(Y, f_sp, f_a, Fa, heat, cool)
      na = n_eq_sweeps_last
      call resid_relnorm(Fa, u, rc, rna)

      f_seed = f_a
      do j = 1, N
         jsrc = min(N, j + nshift)
         f_seed(j,:) = f_a(jsrc,:)
      enddo
      dseed = maxval(abs(f_seed(1:N,:) - f_a(1:N,:))                      &
                     /max(abs(f_a(1:N,:)), eq_sweep_floor))

      call eval_residual(Y, f_seed, f_b, Fb, heat, cool)
      nb = n_eq_sweeps_last
      call resid_relnorm(Fb, u, rc, rnb)
      Fd = Fb - Fa
      call resid_relnorm(Fd, u, rc, rnd)

      write(*,'(A,ES9.2,A,I0,A,I0)') ' (resid_sc_probe) displaced seed,'//&
           ' relative composition change ', dseed, '; passes from the'//   &
           ' run''s composition ', na, ', from the displaced seed ', nb
      write(*,'(A,ES11.3,A,ES11.3)') ' (resid_sc_probe) ||R|| from the'// &
           ' run''s composition ', rna, ', from the displaced seed ', rnb
      write(*,'(A,ES11.3)') ' (resid_sc_probe) seed dependence of the'//  &
           ' residual, in the units of Resid tol: ', rnd
      stop
      end subroutine probe_eliminated_composition_is_seed_independent

      ! ------------------------------------------------------!

      subroutine row_maxima_of_the_certification(Fvec, u, rmax, jcell,   &
                                                 iwhich)
      ! THE THREE ROW CLASSES OF THE STATIONARY SYSTEM, each as the largest
      ! value the certification measures over the column, AND THE CELL that
      ! carries it:
      !   1  hydrodynamic, |R_kj|/residual_row_scale, the PLAIN measure
      !      that the acceptance gate reads against resid_tol_of_solve --
      !      not the judged distance, which divides each of the three rows
      !      by the tolerance it carries at each cell
      !      (hydrodynamic_distance_from_certification_by_cell);
      !   2  elemental transport, on the operator's own scale, against
      !      cert_tol_element_at;
      !   3  carrier, on the row's own terms, against cert_tol_carrier_at.
      !
      ! These three are the WHOLE physical column, gated cells and
      ! reported cells together, because what they answer is where the
      ! largest row of each class sits. What decides is the gated part
      ! alone (certified_row_measures).
      !
      ! certified_row_measures forms the same three numbers and remains the
      ! definition of them; what this adds is WHERE each one sits, which a
      ! maximum over cells does not carry and which is what says whether two
      ! evaluations of one state differ in one front cell or over a region.
      ! The hydrodynamic maximum is over the whole physical column, which is
      ! the number resid_relnorm returns: that gate combines the layer and
      ! the wind by the LARGER of the two, and a maximum over the union of
      ! two windows is the larger of the two window maxima.
      !
      ! iwhich names the row inside its class: the hydrodynamic row index
      ! 1..3, the element (0 helium, otherwise the trace element's index in
      ! the metal table) and the carrier index. Zero where the class carries
      ! no row.
      !
      ! The element and carrier entries are read from the products of the
      ! LAST residual evaluation (erow_*, escale_*, carrier_cellmax_last),
      ! so the caller evaluates and measures with nothing in between, which
      ! is what certified_row_measures requires of its caller for the same
      ! reason.
      real*8, dimension(nvar_jac*N),  intent(in)  :: Fvec
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rmax
      integer, dimension(3),          intent(out) :: jcell, iwhich
      real*8, dimension(N) :: rr, ss
      real*8  :: q
      integer :: i, j, k, jw
      logical :: fin
      rmax   = 0.0d0
      jcell  = 0
      iwhich = 0
      do k = 1, 3
         do j = 1, N
            rr(j) = Fvec(nvar_jac*(j-1)+k)
            ss(j) = residual_row_scale(k, j, u)
         enddo
         call certification_row_measure(N, N+1, rr, ss, q, jw, fin)
         if (q .gt. rmax(1)) then
            rmax(1) = q;  jcell(1) = jw;  iwhich(1) = k
         endif
      enddo
      if (nspec_row .le. 0) return
      do i = 1, nspec_row
         select case (srow_kind(i))
         case (srow_carrier)
            if (carrier_cellmax_last .gt. rmax(3)) then
               rmax(3)   = carrier_cellmax_last
               jcell(3)  = carrier_cell_worst
               iwhich(3) = srow_idx(i)
            endif
         case (srow_element_he)
            if (.not. allocated(erow_he)) cycle
            call certification_row_measure(N, j_min, erow_he(1:N),        &
                                           escale_he(1:N), q, jw, fin)
            if (q .gt. rmax(2)) then
               rmax(2) = q;  jcell(2) = jw;  iwhich(2) = 0
            endif
         case default
            if (.not. allocated(erow_tr)) cycle
            if (.not. erow_tr_carried(srow_idx(i))) cycle
            call certification_row_measure(N, j_min,                      &
                     erow_tr(1:N,srow_idx(i)), escale_tr(1:N,srow_idx(i)),&
                     q, jw, fin)
            if (q .gt. rmax(2)) then
               rmax(2) = q;  jcell(2) = jw;  iwhich(2) = srow_idx(i)
            endif
         end select
      enddo
      end subroutine row_maxima_of_the_certification

      ! ------------------------------------------------------!

      subroutine replay_distance(Fa, Fb, Drow, ndiff, dcl, jcl, kcl)
      ! THE DISTANCE BETWEEN TWO EVALUATIONS OF ONE STATE, by row class.
      !
      ! Contract 1 of the residual is that F(Y) comes back the same bits
      ! after any number of discarded evaluations of other states, so the
      ! first number is the count of entries that differ at all. A count
      ! alone does not say whether what moved matters, so the entries are
      ! also measured in the units the Newton scales them in -- the row
      ! scale Drow, which is the scale the merit and the step control
      ! divide by -- and reported separately for the three hydrodynamic
      ! rows and for the species rows, with the cell and the slot that
      ! carry each maximum. The two classes are kept apart because their
      ! scales are different quantities (cell_row_scales says why) and one
      ! number over both would be whichever class happens to be larger.
      !
      !   dcl(1), jcl(1), kcl(1)   hydrodynamic rows, the cell and the row
      !   dcl(2), jcl(2), kcl(2)   species rows, the cell and the slot
      real*8, dimension(nvar_jac*N), intent(in)  :: Fa, Fb, Drow
      integer,                       intent(out) :: ndiff
      real*8,  dimension(2),         intent(out) :: dcl
      integer, dimension(2),         intent(out) :: jcl, kcl
      real*8  :: q
      integer :: j, k, is, icl
      ndiff = 0
      dcl = 0.0d0;  jcl = 0;  kcl = 0
      do j = 1, N
         do k = 1, nvar_jac
            is = nvar_jac*(j-1) + k
            if (Fa(is) .ne. Fb(is)) ndiff = ndiff + 1
            icl = 1
            if (k .gt. 3) icl = 2
            q = abs(Fa(is) - Fb(is))/max(Drow(is), 1.0d-300)
            if (q .gt. dcl(icl)) then
               dcl(icl) = q;  jcl(icl) = j;  kcl(icl) = k
            endif
         enddo
      enddo
      end subroutine replay_distance

      ! ------------------------------------------------------!

      subroutine certification_row_class_name(k, txt)
      ! The three row classes of certified_row_measures by name, in one
      ! place, so that a report and a test read the same word.
      integer,          intent(in)  :: k
      character(len=*), intent(out) :: txt
      select case (k)
      case (1);       txt = 'hydrodynamic'
      case (2);       txt = 'element'
      case (3);       txt = 'carrier'
      case default;   txt = 'unnamed'
      end select
      end subroutine certification_row_class_name

      ! ------------------------------------------------------!

      integer function front_of_the_steepest_fraction(f) result(jf)
      ! THE CELL OF A FRONT in a species fraction: the largest relative
      ! change between neighboring cells. It is a front of the PROFILE and
      ! not of the grid -- the quantity is |f(j+1)-f(j)| divided by the
      ! larger of the two, so a fraction falling by decades over a few cells
      ! is found wherever the cells are and a smooth region scores nothing.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: f
      real*8  :: d, dmx
      integer :: j
      jf = 1;  dmx = -1.0d0
      do j = 1, N-1
         d = abs(f(j+1) - f(j))                                          &
             /max(abs(f(j+1)), abs(f(j)), eq_sweep_floor)
         if (d .gt. dmx) then
            dmx = d;  jf = j
         endif
      enddo
      end function front_of_the_steepest_fraction

      ! ------------------------------------------------------!

      subroutine branch_and_closure_of_the_seed_family(Y, f_sp, u)
      ! WHAT THE ELIMINATED CLOSURE DOES TO THE OUTER RESIDUAL, AND WHICH
      ! ROOT OF THE FRONT EACH SEED REACHES, on one state Y.
      !
      ! The residual of the stationary system is R(Y) = L(Y) + S(Y, c*(Y)),
      ! and c*(Y) is a fixed point of the composition elimination. Two
      ! statements about it have to be kept apart, because the physics
      ! allows one to fail while the other holds:
      !
      !  * CONDITIONAL CLOSURE CONVERGENCE. Seeds that reach the SAME root
      !    must give the same outer residual to within a budget, and the
      !    budget is on the ROW MAXIMA of the certification: the window
      !    integrals of the assembled residual are not Lipschitz in the
      !    state in the odd-even band above the base -- one round trip at
      !    4.9e-16 in the composition moves the 1.03 to 1.10 mass window by
      !    16 percent -- while the row maxima and their cells stay put, so a
      !    sum over cells is not a quantity a seed budget can be stated on.
      !  * BRANCH DISCOVERY. A seed that reaches a DIFFERENT root is a
      !    statement about the state and not an error: the ionization front
      !    of this configuration has two exact roots, both with a reaction
      !    residual 1e-25 against an acceptance of 1e-6, differing by 8.9x
      !    in x(H2) at the H2 front. It is RECORDED with its composition and
      !    its residual. Nothing is averaged and no seed is forced onto a
      !    preferred root: which root a stationary residual may land on is a
      !    physics rule the code does not have.
      !
      ! THE BRANCH DATUM is the composition the sweep accepted at the two
      ! fronts of the reference seed, and the grouping is by the largest
      ! relative composition difference over the column. Two seeds are on
      ! one branch when that difference is below distinct_root_separation;
      ! the elimination reaches its fixed point to eq_sweep_reltol = 1e-8,
      ! so a separation six decades above it is not the same root under any
      ! reading of the tolerance, and every measured distance is printed
      ! beside the grouping so that the grouping can be checked.
      !
      ! THE SEED FAMILY is the state's own composition with the H2 partition
      ! scaled by (1 + delta): the eliminated direction the closure was
      ! measured in and the one a Newton step moves. Zero is in the family
      ! and is the reference.
      !
      ! Enabled by EXHALE_RESID_BRANCH_REPORT=1, which STOPS the run: the
      ! family is one measurement on one state and there is nothing after it
      ! to measure.
      real*8, dimension(nvar_jac*N),           intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      ! THE FAMILY. Seven seeds scale one ELIMINATED partition of the
      ! state's own composition by (1 + delta) -- the direction the closure
      ! was measured in, and the one a Newton step moves -- and three
      ! displace the whole composition by 1, 5 and 10 cells, which is the
      ! construction that found the second root of the front: every cell is
      ! then given the composition of a neighbor a few scale heights away,
      ! positive, normalized and not the equilibrium of the cell it is given
      ! to. Seed 1 is the state's own composition and is the reference.
      integer, parameter :: nscale = 7, nshift = 3, nseed = nscale + nshift
      real*8,  parameter :: dlist(nscale) = (/ 0.0d0, 1.0d-6, -1.0d-6,   &
                                               1.0d-4, -1.0d-4,          &
                                               1.0d-2, -1.0d-2 /)
      integer, parameter :: slist(nshift) = (/ 1, 5, 10 /)
      ! WHAT COUNTS AS A DIFFERENT ROOT. The two exact roots this
      ! configuration is known to have differ by 8.9x in x(H2) at the H2
      ! front and by 19 percent in x(H I) at the ionization front
      ! (ISSUES_20260909 3.1 c), so a rule at a tenth of the quantity
      ! separates them; and the largest scaled seed of the family is 1e-2,
      ! so a seed perturbation RETAINED IN FULL does not open a root of its
      ! own, which a rule at the perturbation's own size would do. Seven
      ! decades above the fixed point the elimination reaches
      ! (eq_sweep_reltol = 1e-8), and every measured distance is printed
      ! beside the grouping so that the grouping can be checked against
      ! another rule without another run.
      real*8,  parameter :: distinct_root_separation = 1.0d-1
      real*8, allocatable :: comp(:,:,:)
      real*8, dimension(1-Ng:N+Ng,n_species) :: seed, cwork
      real*8, dimension(nvar_jac*N) :: Fvec
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool, npart_d
      real*8  :: rmax(3,nseed), sep(nseed), tol(3), rmx(3)
      integer :: jcl(3), iwh(3)
      integer :: jcell(3,nseed), iwhich(3,nseed), isp_seed
      integer :: nacc4(nseed), nacc6(nseed)
      real*8  :: datum(nseed)
      integer :: root_of(nseed), jsep(nseed), ssep(nseed)
      real*8  :: spread_row, lo, hi, dmx, d
      integer :: is, k, j, iroot, nroot, ncount, jh2, jhi, isc
      character(len=8) :: cfam
      ! `count` is a global integer of this build (global_parameters, the
      ! marching step counter), so the intrinsic of that name is not
      ! available here and a seed tally is written out.
      character(len=8)  :: env
      character(len=16) :: cname

      call get_environment_variable('EXHALE_RESID_BRANCH_REPORT', env)
      if (trim(env) .ne. '1') return

      allocate(comp(1-Ng:N+Ng,n_species,nseed))
      ! THE REFERENCE OF EACH ROW CLASS. The hydrodynamic entry of
      ! row_maxima_of_the_certification is the PLAIN row measure, the
      ! largest of the three rows with no tolerance in it, which is the
      ! number the acceptance gate reads; so its reference here is the
      ! gate's own target and not one of the three certification tolerances
      ! (which row of the three carries the maximum is printed beside it as
      ! iwhich, and dividing a mixed-row maximum by one row's tolerance
      ! would be a ratio of two different rows). The judged distance, each
      ! row over the tolerance it carries at each cell, is
      ! hydrodynamic_distance_from_certification_by_cell and is what the
      ! ledger ranks on.
      tol(1) = max(resid_tol_of_solve, 1.0d-300)
      tol(2) = cert_tol_element_at(cert_regime_wind_r)
      tol(3) = cert_tol_carrier_at(cert_regime_wind_r)

      ! WHICH PARTITION THE SCALED SEEDS MOVE. It has to be an ELIMINATED
      ! one that the configuration HOLDS. Where a partition is transported
      ! the sweep replaces the seeded value by the transported one (the
      ! x_h2_fix row of the network), and where the run carries no molecules
      ! the H2 column is identically zero, so scaling either measures
      ! nothing. H2 is the direction the seed dependence was found in and is
      ! used wherever the run is molecular and H2 is eliminated; otherwise
      ! the family moves the H+ partition, which is the partition of the
      ! ionization front.
      isp_seed = isp_HII
      if (thereis_mol) then
         isp_seed = isp_H2
         do k = 1, nspec_row
            if (srow_kind(k) .eq. srow_carrier .and.                      &
                srow_isp(k) .eq. isp_H2) isp_seed = isp_HII
         enddo
      endif
      do is = 1, nseed
         seed = f_sp
         if (is .le. nscale) then
            do j = 1, N
               seed(j,isp_seed) = min(1.0d0, max(0.0d0,                   &
                    f_sp(j,isp_seed)*(1.0d0 + dlist(is))))
            enddo
         else
            do j = 1, N
               seed(j,:) = f_sp(min(N, j + slist(is-nscale)),:)
            enddo
         endif
         ! The state is not adopted: the family is a probe of the closure
         ! and every one of its evaluations is discarded.
         ! The state is not adopted, so the evaluation puts back what it
         ! overwrote and the next seed of the family starts from the same
         ! frozen background as this one; the row maxima are therefore
         ! taken INSIDE the evaluation, where the rows of this seed exist.
         call eval_residual(Y, seed, cwork, Fvec, heat, cool,             &
                            state_is_discarded=.true.,                    &
                            rowmax_judged=rmx, cells_judged=jcl,          &
                            which_judged=iwh, n_part=npart_d)
         comp(:,:,is)   = cwork
         ! The two NON-ROOT acceptance classes of the last sweep of this
         ! evaluation: a cell in class 4 or 6 has no root of the network at
         ! this state and keeps what it was seeded with, so a composition
         ! difference sitting on such a cell is a cell the chemistry could
         ! not solve and not a second root of it. The census is not among
         ! the products a discarded evaluation puts back, so it is the one
         ! this seed left.
         nacc4(is) = ieq_sweep_ledger_last%acc_n(4)
         nacc6(is) = ieq_sweep_ledger_last%acc_n(6)
         rmax(:,is)     = rmx
         jcell(:,is)    = jcl
         iwhich(:,is)   = iwh
      enddo

      ! The two fronts of the reference seed, from its own composition: the
      ! cell of the steepest relative change of the H2 partition and of the
      ! neutral hydrogen fraction. A front is where the closure has two
      ! roots to have, so this is where the branch datum is read.
      jh2 = front_of_the_steepest_fraction(comp(:,isp_H2,1))
      jhi = front_of_the_steepest_fraction(comp(:,1,1))

      do is = 1, nseed
         call worst_composition_change(comp(:,:,is), comp(:,:,1), dmx,    &
                                       jsep(is), ssep(is))
         sep(is)   = dmx
         datum(is) = branch_datum_distance(comp(:,:,is), comp(:,:,1),     &
                                           jh2, jhi)
      enddo

      ! GROUPING BY THE BRANCH DATUM, not averaging: each seed joins the
      ! first root whose composition AT THE FRONTS it agrees with below the
      ! separation, and opens a new root otherwise. The datum is the front
      ! and not the whole column because that is what a branch of this
      ! closure is: away from a front the closure has one root and what
      ! differs there is a trace species the network left where it was
      ! seeded, which the non-root census beside each seed names.
      root_of = 0
      nroot   = 0
      do is = 1, nseed
         do iroot = 1, nroot
            do k = 1, nseed
               if (root_of(k) .ne. iroot) cycle
               d = branch_datum_distance(comp(:,:,is), comp(:,:,k),       &
                                         jh2, jhi)
               if (d .le. distinct_root_separation) root_of(is) = iroot
               exit
            enddo
            if (root_of(is) .gt. 0) exit
         enddo
         if (root_of(is) .eq. 0) then
            nroot = nroot + 1
            root_of(is) = nroot
         endif
      enddo

      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)')                               &
           ' (resid_branch) seed family of ', nseed,                      &
           ': ', nscale, ' scalings of species column ', isp_seed,        &
           ' and ', nshift, ' displacements; H2 front cell ', jh2
      write(*,'(A,I0)') ' (resid_branch) H I front cell ', jhi
      write(*,'(A)') ' (resid_branch) seed  scale/shift  root |'//        &
           ' hydro (cell,row)  element (cell,elem)  carrier (cell,ic) |'//&
           '  branch datum |  max|df/f| against the reference'//          &
           ' (cell, species) |  cells with no root (class 4 / class 6)'
      do is = 1, nseed
         if (is .le. nscale) then
            write(*,'(A,I3,A,ES10.2,I5,3(ES12.4,I5,I4),2ES11.3,I5,I4,'//  &
                 '2I7)')                                                  &
                 ' (resid_branch) ', is, ' scale ', dlist(is),            &
                 root_of(is),                                             &
                 rmax(1,is), jcell(1,is), iwhich(1,is),                   &
                 rmax(2,is), jcell(2,is), iwhich(2,is),                   &
                 rmax(3,is), jcell(3,is), iwhich(3,is),                   &
                 datum(is), sep(is), jsep(is), ssep(is),                  &
                 nacc4(is), nacc6(is)
         else
            write(*,'(A,I3,A,I10,I5,3(ES12.4,I5,I4),2ES11.3,I5,I4,'//     &
                 '2I7)')                                                  &
                 ' (resid_branch) ', is, ' shift ', slist(is-nscale),     &
                 root_of(is),                                             &
                 rmax(1,is), jcell(1,is), iwhich(1,is),                   &
                 rmax(2,is), jcell(2,is), iwhich(2,is),                   &
                 rmax(3,is), jcell(3,is), iwhich(3,is),                   &
                 datum(is), sep(is), jsep(is), ssep(is),                  &
                 nacc4(is), nacc6(is)
         endif
      enddo

      ! CONTRACT 2, on the seeds of ONE branch. The spread is taken inside
      ! each root, because a spread across roots is not a closure error.
      !
      ! It is reported twice. The SCALED seeds are the closure-convergence
      ! family: small admissible perturbations of one eliminated partition,
      ! which is the question contract 2 asks, and their spread is the
      ! budget to compare with the tolerance. ALL the seeds of the root
      ! include the displaced compositions, which are the branch-discovery
      ! family: a seed given a neighbor's whole composition is not a small
      ! perturbation of this state and its spread is reported, not budgeted.
      do isc = 1, 2
      do iroot = 1, nroot
         ncount = 0
         do is = 1, nseed
            if (root_of(is) .ne. iroot) cycle
            if (isc .eq. 1 .and. is .gt. nscale) cycle
            ncount = ncount + 1
         enddo
         if (ncount .le. 0) cycle
         cfam = 'scaled'
         if (isc .eq. 2) cfam = 'all'
         do k = 1, 3
            lo =  huge(1.0d0);  hi = -huge(1.0d0)
            do is = 1, nseed
               if (root_of(is) .ne. iroot) cycle
               if (isc .eq. 1 .and. is .gt. nscale) cycle
               lo = min(lo, rmax(k,is));  hi = max(hi, rmax(k,is))
            enddo
            spread_row = hi - lo
            call certification_row_class_name(k, cname)
            write(*,'(A,A,A,I2,A,I0,A,A,A,ES11.3,A,ES11.3,A,ES11.3)')     &
                 ' (resid_branch) ', trim(cfam), ' seeds of root ',       &
                 iroot, ', ', ncount,                                     &
                 ' of them, row class ', trim(cname), ': maximum ', hi,   &
                 ', spread ', spread_row, ', spread over tolerance ',     &
                 spread_row/tol(k)
         enddo
      enddo
      enddo

      ! CONTRACT 3. Each root named by the cells the branch datum is read
      ! at, its composition there and its residual.
      ! THE BUDGET PER PERTURBATION SIZE, which is the sensitivity contract
      ! 2 is really about: how far the row maximum of the SAME state moves
      ! when the seed of its closure moves by a stated amount, in the units
      ! of the tolerance that row is judged with. A pair of seeds at the
      ! same |delta| bounds it from both sides.
      do k = 1, 3
         call certification_row_class_name(k, cname)
         do is = 2, nscale, 2
            lo = abs(rmax(k,is)   - rmax(k,1))
            hi = abs(rmax(k,is+1) - rmax(k,1))
            write(*,'(A,A,A,ES9.2,A,ES11.3,A,ES11.3)')                    &
                 ' (resid_branch) closure budget, row class ',            &
                 trim(cname), ', seed ', abs(dlist(is)),                  &
                 ': largest row-maximum move ', max(lo, hi),              &
                 ', over the tolerance ', max(lo, hi)/tol(k)
         enddo
      enddo
      ! THE RETENTION of a scaled seed at the front: how much of the seed
      ! perturbation the closure leaves in the accepted composition there.
      ! One is a perturbation retained in full, zero a closure that forgets
      ! its seed; above one the front amplifies it. It is the quantity B5h
      ! measured on the H2 partition of the eliminated-H2 configuration
      ! (0.77 at the base) and is measured here on the partition this
      ! configuration eliminates.
      do is = 2, nscale
         write(*,'(A,ES10.2,A,ES11.3)') ' (resid_branch) scaled seed ',   &
              dlist(is), ': retention of the seed at the fronts ',        &
              datum(is)/max(abs(dlist(is)), 1.0d-300)
      enddo
      dmx = 0.0d0
      do is = 1, nseed
         dmx = max(dmx, datum(is))
      enddo
      write(*,'(A,ES11.3,A,ES11.3)')                                      &
           ' (resid_branch) largest branch-datum separation in the'//     &
           ' family ', dmx, ', separation rule ', distinct_root_separation
      write(*,'(A,I0,A,I0,A)') ' (resid_branch) roots found ', nroot,     &
           ' among ', nseed, ' seeds'
      do iroot = 1, nroot
         ncount = 0
         do is = 1, nseed
            if (root_of(is) .eq. iroot) ncount = ncount + 1
         enddo
         do is = 1, nseed
            if (root_of(is) .ne. iroot) cycle
            write(*,'(A,I2,A,I0,A,I0,A,I5,A,ES22.15,A,ES22.15,A,I5,A,'//  &
                 'ES22.15,A,ES12.4)')                                     &
                 ' (resid_branch) root ', iroot, ' holds ', ncount,       &
                 ' of ', nseed, ' seeds; cell ', jh2, ' x(H2) ',          &
                 comp(jh2,isp_H2,is), ' x(H I) ', comp(jh2,1,is),         &
                 '; cell ', jhi, ' x(H I) ', comp(jhi,1,is),              &
                 '; hydrodynamic row maximum ', rmax(1,is)
            exit
         enddo
      enddo
      if (nroot .gt. 1) then
         write(*,'(A)') ' (resid_branch) VERDICT: the seed family'//      &
              ' reaches more than one root of the closure; every root'//  &
              ' is recorded above and none is adopted.'
      else
         write(*,'(A)') ' (resid_branch) VERDICT: every seed of the'//    &
              ' family reaches one root of the closure.'
      endif
      deallocate(comp)
      write(*,*) '(steady_newton) EXHALE_RESID_BRANCH_REPORT=1: stopping.'
      stop
      end subroutine branch_and_closure_of_the_seed_family

      ! ------------------------------------------------------!

      double precision function branch_datum_distance(fa, fb, jh2, jhi)   &
                               result(d)
      ! THE DISTANCE BETWEEN TWO SEEDS' BRANCH DATA: the largest relative
      ! difference of the accepted composition AT THE FRONT CELLS, which is
      ! where a closure with more than one root has them. Its measure and
      ! its floor are the elimination's own stopping test, so this number
      ! and the fixed-point tolerance compare directly.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: fa, fb
      integer, intent(in) :: jh2, jhi
      real*8  :: q
      integer :: i, j, k
      d = 0.0d0
      do k = 1, 2
         j = jh2
         if (k .eq. 2) j = jhi
         if (j .lt. 1 .or. j .gt. N) cycle
         do i = 1, n_species
            q = abs(fa(j,i) - fb(j,i))                                    &
                /max(abs(fa(j,i)), eq_sweep_floor)
            d = max(d, q)
         enddo
      enddo
      end function branch_datum_distance

      ! ------------------------------------------------------!

      subroutine cells_of_the_row_report(cells, ncell, ncell_max)
      ! The cells the row diagnostics report at, from
      ! EXHALE_RESID_ROW_CELLS as a comma-separated list, and by default the
      ! base cell, the cell that carries the largest carrier row of the
      ! state just evaluated and the outermost cell -- the base, the front
      ! and the wind of this column.
      integer, intent(in)  :: ncell_max
      integer, intent(out) :: cells(ncell_max), ncell
      character(len=64) :: env
      integer :: i, k, j
      cells = 0
      ncell = 0
      call get_environment_variable('EXHALE_RESID_ROW_CELLS', env)
      if (len_trim(env) .gt. 0) then
         do i = 1, len_trim(env)
            if (env(i:i) .eq. ',') env(i:i) = ' '
         enddo
         read(env,*,iostat=k) (cells(i), i = 1, ncell_max)
         do i = 1, ncell_max
            if (cells(i) .ge. 1 .and. cells(i) .le. N) ncell = i
         enddo
         if (ncell .gt. 0) return
      endif
      j = carrier_cell_worst
      if (j .lt. 1 .or. j .gt. N) j = max(1, N/2)
      cells(1) = 1;  cells(2) = j;  cells(3) = N
      ncell = 3
      end subroutine cells_of_the_row_report

      ! ------------------------------------------------------!

      subroutine row_cancellation_and_the_rounding_floor(f_sp, u)
      ! HOW MUCH OF EACH SPECIES ROW SURVIVES THE CANCELLATION OF ITS OWN
      ! TERMS, and where double precision puts a floor under the row because
      ! of it.
      !
      !   C     = sum |terms| / |sum terms|
      !   floor = epsilon * sum |terms|          [the row's own units]
      !
      ! C is the amplification a relative error in one term suffers on its
      ! way into the row: a row whose terms cancel to one part in C carries
      ! C times the relative error of its largest term. The floor is the
      ! smallest row value a summation in double precision separates from
      ! zero, and what decides whether it matters is the floor DIVIDED BY
      ! THE ROW'S CERTIFICATION SCALE, which is what the tolerance is
      ! written in: certification_row_measure judges
      ! |row|/max(scale, cert_scale_floor), and the scale of both operators
      ! is the same sum of the magnitudes of the row's own terms
      ! (element_transport_residual and carrier_steady_residual state so at
      ! their headers). It is reported against the tolerance of the row's
      ! class; no compensated summation is introduced here.
      !
      ! The three hydrodynamic rows are not in the table: residual_row_scale
      ! bounds each of those rows' largest term rather than summing every
      ! term, so sum|terms| is not a quantity that operator returns and
      ! inventing one here would report a number no acceptance test reads.
      !
      ! Enabled by EXHALE_RESID_ROW_CONDITION=1, which STOPS the run.
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      integer, parameter :: ncell_max = 8
      integer :: cells(ncell_max), ncell, i, j, ic, ie
      real*8  :: res, sc, cnd, flr, judged, tolr, eps
      character(len=64) :: env
      character(len=16) :: cname
      call get_environment_variable('EXHALE_RESID_ROW_CONDITION', env)
      if (trim(env) .ne. '1') return
      eps = epsilon(1.0d0)
      call cells_of_the_row_report(cells, ncell, ncell_max)
      write(*,'(A,ES13.4E3,A,ES13.4E3)') ' (row_condition) epsilon ',    &
           eps,                                                            &
           ', cert_scale_floor ', cert_scale_floor
      write(*,'(A)') ' (row_condition) row  cell  |row|  sum|terms| '//   &
           ' C  floor=eps*sum|terms|  floor in the judged measure '//     &
           ' tolerance (huge where the cell is reported and not gated)'
      do i = 1, ncell
         j = cells(i)
         if (j .lt. 1 .or. j .gt. N) cycle
         do ie = 1, nspec_row
            select case (srow_kind(ie))
            case (srow_carrier)
               if (.not. allocated(crow_res_last)) cycle
               ic    = srow_idx(ie)
               res   = crow_res_last(j,ic)
               sc    = crow_terms_last(j,ic)
               cname = carrier_name(ic)
               tolr  = cert_tol_carrier_at(r(j))
            case (srow_element_he)
               if (.not. allocated(erow_he)) cycle
               res   = erow_he(j);  sc = escale_he(j)
               cname = 'He/H'
               tolr  = cert_tol_element_at(r(j))
            case default
               if (.not. allocated(erow_tr)) cycle
               if (.not. erow_tr_carried(srow_idx(ie))) cycle
               res   = erow_tr(j,srow_idx(ie))
               sc    = escale_tr(j,srow_idx(ie))
               cname = melem_name(srow_idx(ie))
               tolr  = cert_tol_element_at(r(j))
            end select
            cnd    = row_condition_number(res, sc)
            flr    = eps*sc
            judged = row_rounding_floor_in_the_measure(sc)
            write(*,'(A,A6,I6,6ES15.4E3)') ' (row_condition) ',           &
                 trim(cname), j, abs(res), sc, cnd, flr, judged, tolr
         enddo
      enddo
      call chemical_decay_against_the_transport_times(f_sp, u)
      write(*,*) '(steady_newton) EXHALE_RESID_ROW_CONDITION=1: stopping.'
      stop
      end subroutine row_cancellation_and_the_rounding_floor

      ! ------------------------------------------------------!

      subroutine chemical_decay_against_the_transport_times(f_sp, u)
      ! THE DECISION THIS DIAGNOSTIC IS FOR: is a molecular carrier's local
      ! chemistry fast against the transport of the same cell, so that it
      ! may stay an eliminated closure, or is it a slow mode that belongs in
      ! the unknown vector?
      !
      ! The chemical side is measured by ionization_equilibrium at the
      ! ACCEPTED state of the cell, on the reaction Jacobian of the local
      ! system with the conservation modes removed
      ! (molecular_decay_rate_of_the_cell says how, and why the code's own
      ! fraction parameterization already carries the elemental
      ! constraints). Two numbers of it are read: the row's own decay rate,
      ! which is attributable to one molecule, and the slowest decay of the
      ! whole coupled spectrum, which is the one that decides for the system.
      !
      ! The transport side is measured on the RESOLVED GRADIENT LENGTH of
      ! the molecule's own profile, L = f/|df/dr| by a central difference,
      ! and NOT on the cell width: the cell width is a property of the grid,
      ! and a smooth profile on a fine grid would report an arbitrarily
      ! short transport time.
      !   t_adv  = L/|v|          t_diff = L^2/D
      ! with D the carrier's own diffusivity. A ratio t*lambda far above one
      ! is a mode the closure may eliminate; a ratio at or below one is a
      ! mode transport sets, and the outcome of the comparison is what
      ! decides promotion into Y or a tighter closure.
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      integer, parameter :: ncell_max = 8
      integer, parameter :: nmolrow = 4
      ! The four molecular rows of the local network, in the layout order of
      ! System_HeH_mol, beside the f_sp column of the same molecule.
      integer, parameter :: molsp(nmolrow) =                              &
                 (/ isp_H2, isp_H2p, isp_H3p, isp_HeHp /)
      integer :: cells(ncell_max), ncell, i, j, k, ic
      real*8  :: lam_row, lam_slow, lam_fast, L, tadv, tdiff, D, vcm
      real*8, dimension(3,1-Ng:N+Ng) :: W
      character(len=8) :: mname
      logical :: have
      if (.not. thereis_mol) return
      call cells_of_the_row_report(cells, ncell, ncell_max)
      call U_to_W(u, W)
      write(*,'(A)') ' (decay_vs_transport) molecule cell  L[cm]  '//     &
           'D[cm2/s]  1/|lambda_row|[s]  1/lambda_slow[s]  '//            &
           '1/lambda_fast[s]  t_adv[s]  t_diff[s]  t_adv*lambda_row  '//  &
           't_diff*lambda_row   (a zero is a quantity this path does '//  &
           'not define; see the routine)'
      do i = 1, ncell
         j = cells(i)
         if (j .lt. 1 .or. j .gt. N) cycle
         do k = 1, nmolrow
            call molecular_decay_rate_of_the_cell(j, k, mname, lam_row,   &
                                                  lam_slow, lam_fast, have)
            if (.not. have) cycle
            L    = resolved_gradient_length(f_sp, molsp(k), j)
            vcm  = abs(W(2,j))*v0
            ic   = carrier_index_of_species(molsp(k))
            D    = 0.0d0
            if (ic .gt. 0) D = carrier_diffusion_coefficient(j, ic)
            ! A time or a ratio that is not defined is written as zero, and
            ! the header says so: the only diffusivity this path forms is
            ! that of a TRANSPORTED carrier, so a molecule the run does not
            ! transport has no diffusive time here, and a cell at rest has
            ! no advective one.
            tadv  = 0.0d0
            tdiff = 0.0d0
            if (vcm .gt. 0.0d0) tadv  = L/vcm
            if (D   .gt. 0.0d0) tdiff = L*L/D
            write(*,'(A,A6,I6,9ES14.4E3)') ' (decay_vs_transport) ',      &
                 trim(mname), j, L, D,                                    &
                 1.0d0/max(abs(lam_row), 1.0d-300),                       &
                 1.0d0/max(lam_slow,1.0d-300),                            &
                 1.0d0/max(lam_fast,1.0d-300), tadv, tdiff,               &
                 tadv*lam_row, tdiff*lam_row
         enddo
      enddo
      end subroutine chemical_decay_against_the_transport_times

      ! ------------------------------------------------------!

      integer function carrier_index_of_species(isp) result(ic)
      ! Which transported carrier names species isp, and 0 where the run
      ! transports none of it. The map is the registry's own
      ! (carrier_species_index) read the other way round.
      integer, intent(in) :: isp
      integer :: i
      ic = 0
      do i = 1, n_carrier
         if (carrier_species_index(i) .eq. isp) ic = i
      enddo
      end function carrier_index_of_species

      ! ------------------------------------------------------!

      double precision function resolved_gradient_length(f_sp, isp, j)    &
                               result(L)
      ! The length over which species isp's fraction changes by itself,
      ! L = f/|df/dr| from a central difference of the profile the state
      ! carries. It is the physical gradient scale of the state, bounded
      ! below by half the two-cell width because a profile cannot be
      ! resolved on less than that, and bounded above by the domain, which
      ! is what a uniform profile returns. [cm]
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      integer, intent(in) :: isp, j
      real*8  :: fm, fp, f0, dfdr, dr2
      L = (r(N) - r(1))*R0
      if (j .lt. 2 .or. j .gt. N-1) return
      fm = f_sp(j-1,isp);  f0 = f_sp(j,isp);  fp = f_sp(j+1,isp)
      dr2 = (r(j+1) - r(j-1))*R0
      if (dr2 .le. 0.0d0) return
      dfdr = (fp - fm)/dr2
      if (abs(dfdr) .le. 0.0d0) return
      L = min(L, max(abs(f0/dfdr), 0.5d0*dr2))
      end function resolved_gradient_length

      ! ------------------------------------------------------!

      double precision function row_condition_number(row, term_sum)      &
                               result(C)
      ! HOW MUCH OF A ROW SURVIVES THE CANCELLATION OF ITS OWN TERMS:
      !   C = sum |terms| / |sum terms| ,
      ! the amplification a relative error in one term suffers on its way
      ! into the row. A row whose terms cancel to one part in C carries C
      ! times the relative error of its largest term. C >= 1 always, and
      ! C = 1 is a row with no cancellation at all. A row that is exactly
      ! zero has no finite C and the largest representable number is
      ! returned for it, which is the statement that no relative accuracy
      ! survives.
      real*8, intent(in) :: row, term_sum
      if (abs(row) .le. 0.0d0) then
         C = huge(1.0d0)
      else
         C = term_sum/abs(row)
      endif
      end function row_condition_number

      ! ------------------------------------------------------!

      double precision function row_rounding_floor_in_the_measure(       &
                               term_sum) result(f)
      ! THE SMALLEST ROW VALUE DOUBLE PRECISION SEPARATES FROM ZERO, in the
      ! units the acceptance tolerance is written in.
      !
      ! Summing the terms of a row in double precision leaves a rounding
      ! residue of order epsilon times the sum of their magnitudes; the
      ! certification then judges |row|/max(sum|terms|, cert_scale_floor)
      ! (certification_row_measure), so in the judged measure the floor is
      !   epsilon * sum|terms| / max(sum|terms|, cert_scale_floor)
      ! which is epsilon wherever the scale stands above the floor of the
      ! measure and smaller where it does not. It is the number a row
      ! tolerance has to stand above to be a statement about the physics
      ! rather than about the arithmetic.
      real*8, intent(in) :: term_sum
      f = epsilon(1.0d0)*term_sum/max(term_sum, cert_scale_floor)
      end function row_rounding_floor_in_the_measure

      ! ------------------------------------------------------!

      double precision function carrier_density_from_unknown(yv, j)       &
                               result(nc)
      ! The carrier density in the code's units that the unknown of cell j
      ! names: the unknown itself, or its exponential where the carrier is
      ! carried in logarithm. The exponent is capped so that a Newton step
      ! into a state no gas can be cannot produce an infinity; the headroom
      ! comparison and the upper bound refuse it as they refuse any other
      ! such state.
      real*8,  intent(in) :: yv
      integer, intent(in) :: j
      if (carrier_unknown_is_logarithmic) then
         nc = exp(min(yv, 1.0d2))
      else
         nc = yv
      endif
      end function carrier_density_from_unknown

      ! ------------------------------------------------------!

      double precision function carrier_unknown_from_density(nc, j, irow)  &
                               result(yv)
      ! The inverse: the unknown that names a carrier density. In logarithm
      ! a density at or below the floor of this cell and this carrier is
      ! named by the floor, because ln 0 is not a number and n = 0 is what
      ! the outermost cell of an ionized wind holds
      ! (carrier_log_floor states which density that is and why).
      real*8,  intent(in) :: nc
      integer, intent(in) :: j, irow
      if (carrier_unknown_is_logarithmic) then
         yv = log(max(nc, carrier_log_floor_at(j, irow)))
      else
         yv = nc
      endif
      end function carrier_unknown_from_density

      ! ------------------------------------------------------!

      logical function carrier_rows_registered()
      ! Whether any transported CARRIER is a row of this solve. A solve whose
      ! only transported balance is an element has no carrier unknown, so
      ! nothing about the carrier unknown space is reported for it.
      integer :: i
      carrier_rows_registered = .false.
      do i = 1, nspec_row
         if (srow_kind(i) .eq. srow_carrier) carrier_rows_registered = .true.
      enddo
      end function carrier_rows_registered

      ! ------------------------------------------------------!

      double precision function carrier_log_floor_at(j, irow) result(fl)
      ! The floor under the carrier of cell j named by species row irow, in
      ! the code's density units, and 1e-300 where the solve has not formed
      ! one (which is every solve without a carrier row). The guard is a
      ! representability guard and not a choice: log of it is finite.
      integer, intent(in) :: j, irow
      fl = 1.0d-300
      if (.not. allocated(carrier_log_floor)) return
      if (j .lt. 1 .or. j .gt. N) return
      if (irow .lt. 1 .or. irow .gt. nspec_row) return
      fl = max(carrier_log_floor(j,irow), 1.0d-300)
      end function carrier_log_floor_at

      ! ------------------------------------------------------!

      subroutine arm_carrier_log_unknown(Y)
      ! THE FLOOR OF THE LOGARITHMIC CARRIER UNKNOWN, one number per cell
      ! and per carrier row, formed once per solve from the state the solve
      ! starts at: the element budget of that carrier in that cell times the
      ! fraction of its element below which a carrier density changes no
      ! observable. Both are stated at the declarations
      ! (carrier_log_floor, carrier_unobservable_fraction_of_element).
      ! With the floor in hand the carrier slots of Y, which were packed and
      ! evaluated as densities, are rewritten as their logarithms and the
      ! space is armed for the rest of the solve.
      !
      ! IT IS FORMED FOR THE DENSITY UNKNOWN AS WELL, so that a run of either
      ! arm is measured against the same yardstick and the reported ratio of
      ! a carrier to its own floor means the same thing in both.
      !
      ! WHERE THE STATE DOES NOT DEFINE IT the logarithm is not available:
      ! carrier_headroom is the element budget the carrier operator itself
      ! froze at the iterate, and a cell whose budget is not a positive
      ! finite density has none to read. The solve then says so and carries
      ! the density unknown, rather than substituting a number of its own.
      real*8, dimension(nvar_jac*N), intent(inout) :: Y
      integer :: j, i, is
      real*8  :: nb
      logical :: budget_known
      if (nspec_row .le. 0) return
      if (allocated(carrier_log_floor)) deallocate(carrier_log_floor)
      allocate(carrier_log_floor(1:N,1:nspec_row))
      carrier_log_floor = 1.0d-300
      budget_known = .true.
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         do j = 1, N
            ! [cm^-3]: the largest density this carrier can reach in this
            ! cell, from the one element table the carrier rows read.
            nb = carrier_headroom(j, srow_idx(i))
            if (.not. (nb .gt. 0.0d0 .and. nb .lt. 1.0d100)) then
               budget_known = .false.
            else
               carrier_log_floor(j,i) =                                   &
                  max(carrier_unobservable_fraction_of_element*nb/n0,     &
                      1.0d-300)
            endif
         enddo
      enddo
      if (.not. budget_known) then
         write(*,'(A)') ' (JFNK) the element budget of a carrier cell is'//&
              ' not available, so ln n has no floor the state defines:'//  &
              ' this solve carries the carrier density itself'
         ! Nothing formed from an unavailable budget is reported as a floor.
         carrier_log_floor = 1.0d-300
      endif
      carrier_unknown_is_logarithmic =                                    &
         carrier_log_unknown_wanted .and. budget_known
      if (carrier_unknown_is_logarithmic) then
         do j = 1, N
            do i = 1, nspec_row
               if (srow_kind(i) .ne. srow_carrier) cycle
               is = nvar_jac*(j-1) + 3 + i
               ! Y holds the density here, by construction of the pack that
               ! ran before the space was armed.
               Y(is) = log(max(Y(is), carrier_log_floor(j,i)))
            enddo
         enddo
      endif
      call report_carrier_unknown_space
      end subroutine arm_carrier_log_unknown

      ! ------------------------------------------------------!

      subroutine report_carrier_unknown_space
      ! What the species unknowns of THIS solve are, and the floor the
      ! logarithm rests on, printed where the solve begins so that a run
      ! stopped before its end still says which space it was in.
      real*8  :: flmin, flmax
      integer :: j, i, ncar
      if (nspec_row .le. 0) return
      write(*,'(A,A)') ' (JFNK) species unknown space: ',                 &
           trim(species_unknown_space_text())
      if (.not. allocated(carrier_log_floor)) return
      flmin = huge(1.0d0);  flmax = 0.0d0;  ncar = 0
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         ncar = ncar + 1
         do j = 1, N
            flmin = min(flmin, carrier_log_floor(j,i))
            flmax = max(flmax, carrier_log_floor(j,i))
         enddo
      enddo
      if (ncar .le. 0) return
      write(*,'(A,ES11.3,A,ES11.3,A)') ' (JFNK) carrier floor ('//        &
           'code density units): smallest ', flmin, ', largest ', flmax,  &
           ' -- 1e-20 of the cell'//"'"//'s own element budget, the'//     &
           ' fraction the carrier row scale calls unobservable'
      end subroutine report_carrier_unknown_space

      ! ------------------------------------------------------!

      subroutine freeze_species_unknown_box(Y)
      ! THE FACES OF THE BOX THE SPECIES UNKNOWNS OF THIS STEP LIVE IN,
      ! formed once from the iterate Y and from the element budget the
      ! carrier operator froze at that same iterate.
      !
      ! AN ELEMENT UNKNOWN is a mass fraction: [0, 1], and there is nothing
      ! state-dependent about either face.
      !
      ! A CARRIER UNKNOWN has three faces and each is a physical statement.
      !   * BELOW, the smallest density the space can name. In the density
      !     unknown that is zero, a density being non-negative; in ln n it
      !     is carrier_log_floor, stated at its declaration.
      !   * ABOVE, the cell's own density: a carrier cannot be more of the
      !     cell than the cell is.
      ! THE ELEMENT BUDGET IS NOT A FACE OF THIS BOX. A carrier cannot hold
      ! more nuclei of its element than the cell has free, but that ceiling
      ! is SHARED with every other carrier of the same element, so it is a
      ! half-space of the cell's unknowns and not a coordinate bound; it is
      ! carried as a row of the step instead
      ! (freeze_element_constraint_rows, element_constraint_rows_on). With
      ! every carrier density non-negative the shared row implies each
      ! carrier's own ceiling, so nothing is lost, and what is gained is
      ! that a direction blocked by the budget is projected ALONG the
      ! constraint rather than losing the whole step: a component with no
      ! room between its value and its own ceiling made the
      ! fraction-to-the-boundary rule return zero for the whole direction.
      !
      ! In the measurement arm EXHALE_ELEMENT_CONSTRAINT_ROWS=0 the budget
      ! comes back under a coordinate face, which is the box the B5j to B5l
      ! measurements were made in, and it is then RAISED to the iterate
      ! wherever the iterate already stands outside it: the budget is frozen
      ! at the iterate and depends on the very unknown it constrains, so a
      ! face cutting through the point the step starts from would leave that
      ! unknown no room on either side. How often that happens is counted
      ! (n_budget_face_at_the_iterate).
      !
      ! The cell-density face is NOT raised: it is a hard statement about
      ! the state, not a frozen by-product, and a carrier above it is a
      ! state to be pulled back in.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      real*8  :: nb, hi, lo, ycar
      integer :: j, i, is
      if (allocated(species_box_lo)) deallocate(species_box_lo)
      if (allocated(species_box_hi)) deallocate(species_box_hi)
      if (nspec_row .le. 0) return
      allocate(species_box_lo(nvar_jac*N), species_box_hi(nvar_jac*N))
      ! The hydrodynamic slots are unbounded here: density and pressure are
      ! screened by U_to_W_interior and the momentum is signed.
      species_box_lo = -huge(1.0d0)
      species_box_hi =  huge(1.0d0)
      do j = 1, N
         do i = 1, nspec_row
            is = nvar_jac*(j-1) + 3 + i
            if (srow_kind(i) .ne. srow_carrier) then
               species_box_lo(is) = 0.0d0
               species_box_hi(is) = 1.0d0
               cycle
            endif
            ! [code density units] the cell's own density, and the element
            ! budget where the carrier operator has one.
            lo = 0.0d0
            hi = max(Y(nvar_jac*(j-1)+1), 0.0d0)
            if (carrier_headroom_known() .and.                          &
                .not. element_constraint_rows_on) then
               nb = carrier_headroom(j, srow_idx(i))/n0
               nb = nb*budget_face_inside
               ycar = carrier_density_from_unknown(Y(is), j)
               if (ycar .gt. nb) then
                  nb = ycar
                  n_budget_face_at_the_iterate =                          &
                     n_budget_face_at_the_iterate + 1
               endif
               hi = min(hi, nb)
            endif
            if (carrier_unknown_is_logarithmic) then
               lo = log(carrier_log_floor_at(j, i))
               hi = log(max(hi, 1.0d-300))
            endif
            species_box_lo(is) = lo
            species_box_hi(is) = max(hi, lo)
         enddo
      enddo
      end subroutine freeze_species_unknown_box

      ! ------------------------------------------------------!

      subroutine build_the_nuclei_table
      ! The nuclei of each element one particle of each f_sp column holds,
      ! read once from the map. It is a table so that the write-back and the
      ! restoration can classify a species in one lookup; the numbers are
      ! nuclei_per_particle's and nothing here is a second stoichiometry.
      integer :: isp, ie
      if (elem_nu_table_known) return
      elem_nu_of_species = 0
      do isp = 1, n_species
         do ie = 1, n_inv_element
            elem_nu_of_species(isp,ie) = nuclei_per_particle(isp, ie)
         enddo
      enddo
      elem_nu_table_known = .true.
      end subroutine build_the_nuclei_table

      ! ------------------------------------------------------!

      subroutine mark_the_carried_species
      ! The f_sp columns THIS solve carries as carrier unknowns. The shared
      ! budget constrains exactly those: a species the solve does not carry
      ! is data of the step and its nuclei are on the other side of the
      ! constraint.
      integer :: i, isp, k, ie, nel
      elem_species_is_carried = .false.
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         if (srow_isp(i) .ge. 1 .and. srow_isp(i) .le. n_species)         &
            elem_species_is_carried(srow_isp(i)) = .true.
      enddo
      elem_species_absorbs = .false.
      do isp = 1, n_species
         if (isp .eq. isp_HeHp) cycle
         if (elem_species_is_carried(isp)) cycle
         nel = 0
         do k = 1, n_element_constraint_rows
            if (elem_nu_of_species(isp,element_constraint_element(k))     &
                .gt. 0) nel = nel + 1
         enddo
         if (nel .ne. 1) cycle
         do k = 1, n_element_constraint_rows
            ie = element_constraint_element(k)
            if (elem_nu_of_species(isp,ie) .gt. 0)                        &
               elem_species_absorbs(isp,ie) = .true.
         enddo
      enddo
      end subroutine mark_the_carried_species

      ! ------------------------------------------------------!

      logical function element_constraint_rows_known() result(known)
      ! Whether this step carries the shared element budget as rows. False
      ! on the three-unknown route, on a solve with no carrier unknown, in
      ! the measurement arm, and before the first residual evaluation has
      ! made the carrier operator freeze a budget.
      known = elem_rows_frozen
      end function element_constraint_rows_known

      ! ------------------------------------------------------!

      double precision function element_constraint_of_the_cell(j, k)      &
                               result(c)
      ! c of row k of cell j at the iterate the rows were frozen at
      ! [cm^-3]: negative inside the shared budget, positive outside it.
      integer, intent(in) :: j, k
      c = 0.0d0
      if (.not. elem_rows_frozen) return
      if (j .lt. 1 .or. j .gt. N) return
      if (k .lt. 1 .or. k .gt. n_element_constraint_rows) return
      c = elem_row_c(j,k)
      end function element_constraint_of_the_cell

      ! ------------------------------------------------------!

      double precision function element_constraint_budget_of_the_cell(j,  &
                               k) result(b)
      ! The budget row k of cell j is measured against [cm^-3], the
      ! conditional budget the carrier operator froze at the iterate.
      integer, intent(in) :: j, k
      b = 0.0d0
      if (.not. elem_rows_frozen) return
      if (j .lt. 1 .or. j .gt. N) return
      if (k .lt. 1 .or. k .gt. n_element_constraint_rows) return
      b = elem_row_budget(j,k)
      end function element_constraint_budget_of_the_cell

      ! ------------------------------------------------------!

      double precision function                                          &
               element_constraint_gradient_of_the_cell(j, k, islot)      &
               result(a)
      ! dc/dY of row k of cell j at slot islot of that cell (1 to nvar_jac),
      ! in the unknown coordinates [cm^-3 per unit unknown].
      integer, intent(in) :: j, k, islot
      a = 0.0d0
      if (.not. elem_rows_frozen) return
      if (j .lt. 1 .or. j .gt. N) return
      if (k .lt. 1 .or. k .gt. n_element_constraint_rows) return
      if (islot .lt. 1 .or. islot .gt. nvar_jac) return
      a = elem_row_grad(islot,j,k)
      end function element_constraint_gradient_of_the_cell

      ! ------------------------------------------------------!

      subroutine freeze_element_constraint_rows(Y, f_sp, D)
      ! THE SHARED ELEMENT BUDGET OF EVERY CELL AS A LINEARIZED ROW OF THIS
      ! STEP, formed once from the iterate and from the budget the carrier
      ! operator froze at that same iterate (item N4b, decision 14 route
      ! (i); docs/constrained_step_design_20260909.md section 2).
      !
      ! WHAT THE ROW IS. For each element a carrier of the cell holds,
      !
      !     c_k(Y) = sum_i nuclei_per_particle(i,ie) n_i(Y) - budget(ie)
      !
      ! with the sum over the carriers THIS solve transports, and c_k <= 0
      ! the feasible half-space. The value and the three derivatives come
      ! from element_constraint_derivatives, so the constraint the step
      ! moves along and the constraint a candidate is judged by are one
      ! expression of one stoichiometry.
      !
      ! WHERE THE BUDGET COMES FROM. The carrier operator freezes one number
      ! per cell and carrier, the largest density that carrier could reach
      ! if it held the whole budget alone (carrier_headroom, which is
      ! element_box_side of the map). Multiplying it back by the nuclei that
      ! carrier holds returns the budget itself, exactly: the hydrogen
      ! available to the carriers from the H2 side, the free oxygen from the
      ! OH side, and the smaller of the free oxygen and the free carbon from
      ! the CO side. The carbon row therefore states n(CO) <= nC_free
      ! wherever carbon is the binding element and repeats the oxygen row
      ! where oxygen is, which is a redundant row and not a tighter one.
      !
      ! WHICH UNKNOWNS THE ROW REACHES.
      !   * the carrier unknowns of the cell, through dn/d ln n = n where
      !     the carrier is carried in logarithm;
      !   * the cell's own density, whose column is c/nd and therefore
      !     VANISHES on the constraint surface: the constraint is
      !     homogeneous of degree one in the density at fixed composition,
      !     so the density direction is tangent to the surface and an active
      !     row does not push the hydrodynamic mass unknown;
      !   * an element unknown of the cell, through the budget. The element
      !     mass fractions and the element densities are related by the
      !     two-component mass closure the element operator states,
      !     m_1 n_H + m_He n_He = rho/n0 with the trace masses inside m_1
      !     (element_mass_fractions), so at fixed density
      !     n_H is proportional to 1 - X_He and a trace element's nucleus
      !     density is proportional to its own mass fraction. The nuclei the
      !     frozen stages hold are held fixed with the budget, which is what
      !     makes this a linearization and not an identity; the nonlinear
      !     constraint is re-evaluated at every candidate
      !     (element_budget_violation_of_composition) and a candidate that
      !     violates it is restored (restore_the_element_budget), so the
      !     linearization supplies a direction and never a guarantee.
      !
      ! It is frozen for the step for the reason the box is: a constraint
      ! that moved with the trial would make the projection a nonlinear map
      ! of the trial, and the operator the Krylov cycle is sampled under
      ! and the step the trial takes would be two different maps.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(nvar_jac*N),          intent(in) :: D
      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Yel
      real*8  :: budget(n_inv_element), c(n_inv_element)
      real*8  :: dc_dn(n_species,n_inv_element)
      real*8  :: dc_dnd(n_inv_element), dc_dbudget(n_inv_element)
      real*8  :: nd, nc, dn_dY, dbudget_dY, xhe, xtr
      integer :: j, i, k, ie, is, islot
      logical :: have_budget

      if (allocated(elem_row_c)) deallocate(elem_row_c)
      if (allocated(elem_row_budget)) deallocate(elem_row_budget)
      if (allocated(elem_row_target)) deallocate(elem_row_target)
      if (allocated(elem_row_grad)) deallocate(elem_row_grad)
      if (allocated(elem_row_grad_scaled)) deallocate(elem_row_grad_scaled)
      if (allocated(elem_row_basis)) deallocate(elem_row_basis)
      if (allocated(n_elem_row_basis)) deallocate(n_elem_row_basis)
      elem_rows_frozen = .false.
      elem_rows_active_set_frozen = .false.
      if (nspec_row .le. 0) return
      if (.not. element_constraint_rows_on) return
      if (.not. carrier_rows_registered()) return
      if (.not. carrier_headroom_known()) return
      call build_the_nuclei_table
      call mark_the_carried_species
      allocate(elem_row_c(1:N,n_element_constraint_rows),                 &
               elem_row_budget(1:N,n_element_constraint_rows),            &
               elem_row_target(1:N,n_element_constraint_rows),            &
               elem_row_grad(nvar_jac,1:N,n_element_constraint_rows),     &
               elem_row_grad_scaled(nvar_jac,1:N,                         &
                                    n_element_constraint_rows),           &
               elem_row_basis(nvar_jac,n_element_constraint_rows,1:N),    &
               n_elem_row_basis(1:N))
      elem_row_c = 0.0d0;  elem_row_budget = 0.0d0
      elem_row_target = 0.0d0
      elem_row_grad = 0.0d0;  elem_row_grad_scaled = 0.0d0
      elem_row_basis = 0.0d0; n_elem_row_basis = 0
      ! The element mass fractions of the iterate, read through the one map
      ! an element fraction has into a species vector, so that the element
      ! column of the row below and the unknown the element row carries are
      ! the same quantity.
      Yel = 0.0d0
      if (element_rows_registered()) call element_mass_fractions(f_sp, Yel)

      do j = 1, N
         nd = max(Y(nvar_jac*(j-1)+1), 0.0d0)*n0
         ! The budget of each element, from the side the operator froze.
         budget = 0.0d0
         have_budget = .true.
         budget(ien_H) = carrier_headroom(j, ic_H2)                       &
                        *dble(nuclei_per_particle(isp_H2, ien_H))
         budget(ien_O) = carrier_headroom(j, ic_OH)                       &
                        *dble(nuclei_per_particle(isp_OH, ien_O))
         budget(ien_C) = carrier_headroom(j, ic_CO)                       &
                        *dble(nuclei_per_particle(isp_CO, ien_C))
         do k = 1, n_element_constraint_rows
            ie = element_constraint_element(k)
            if (budget(ie) .ge. 0.5d0*huge(1.0d0) .or.                    &
                .not. finite_real(budget(ie))) have_budget = .false.
         enddo
         if (.not. have_budget) then
            deallocate(elem_row_c, elem_row_budget, elem_row_target,      &
                       elem_row_grad, elem_row_grad_scaled,               &
                       elem_row_basis, n_elem_row_basis)
            return
         endif
         call element_constraint_derivatives(nd, f_sp(j,:),               &
                                elem_species_is_carried, budget, c,       &
                                dc_dn, dc_dnd, dc_dbudget)
         do k = 1, n_element_constraint_rows
            ie = element_constraint_element(k)
            elem_row_c(j,k)      = c(ie)
            elem_row_budget(j,k) = budget(ie)
            elem_row_target(j,k) = budget(ie) + max(c(ie), 0.0d0)
            if (c(ie) .gt. 0.0d0) n_elem_row_target_at_the_iterate =      &
               n_elem_row_target_at_the_iterate + 1
            ! the cell's own density
            elem_row_grad(1,j,k) = dc_dnd(ie)*n0
            do i = 1, nspec_row
               islot = 3 + i
               is    = nvar_jac*(j-1) + islot
               if (srow_kind(i) .eq. srow_carrier) then
                  nc    = carrier_density_from_unknown(Y(is), j)
                  dn_dY = n0
                  if (carrier_unknown_is_logarithmic) dn_dY = n0*nc
                  elem_row_grad(islot,j,k) = dc_dn(srow_isp(i),ie)*dn_dY
               else if (srow_kind(i) .eq. srow_element_he) then
                  ! The hydrogen nuclei of a cell are proportional to
                  ! 1 - X_He at fixed density; the nuclei the frozen stages
                  ! hold are frozen with the budget.
                  dbudget_dY = 0.0d0
                  if (ie .eq. ien_H) then
                     xhe = min(max(Yel(j,1), 0.0d0), 1.0d0)
                     if (1.0d0 - xhe .gt. 1.0d-30)                        &
                        dbudget_dY = -budget(ien_H)/(1.0d0 - xhe)
                  endif
                  elem_row_grad(islot,j,k) = dc_dbudget(ie)*dbudget_dY
               else if (srow_kind(i) .eq. srow_element_trace) then
                  ! A trace element's nucleus density is proportional to its
                  ! own mass fraction at fixed density.
                  dbudget_dY = 0.0d0
                  if ((ie .eq. ien_O .and. srow_idx(i) .eq. iel_O) .or.   &
                      (ie .eq. ien_C .and. srow_idx(i) .eq. iel_C)) then
                     xtr = Yel(j,1+srow_idx(i))
                     if (xtr .gt. 1.0d-30) dbudget_dY = budget(ie)/xtr
                  endif
                  elem_row_grad(islot,j,k) = dc_dbudget(ie)*dbudget_dY
               endif
            enddo
            do islot = 1, nvar_jac
               elem_row_grad_scaled(islot,j,k) =                          &
                  elem_row_grad(islot,j,k)*D(nvar_jac*(j-1)+islot)
            enddo
         enddo
      enddo
      elem_rows_frozen = .true.
      ! AND THE MAGNITUDE A CANDIDATE IS COMPARED WITH, taken from the state
      ! these rows were frozen at so that the reference and the rows belong
      ! to one iterate. What a row allows is the budget raised to the
      ! iterate's own demand, so this is zero by construction and the screen
      ! of eval_residual is then a statement about the candidate alone; it
      ! is computed rather than assumed so that a state whose write-back
      ! moved its own carriers would be seen.
      v_headroom_iterate = 0.0d0
      do j = 1, N
         do k = 1, n_element_constraint_rows
            v_headroom_iterate = max(v_headroom_iterate,                  &
               relative_violation_of_the_row(elem_row_c(j,k)              &
                                             + elem_row_budget(j,k),      &
                                             elem_row_target(j,k)))
         enddo
      enddo
      end subroutine freeze_element_constraint_rows

      ! ------------------------------------------------------!

      subroutine fix_active_element_constraints(g)
      ! WHICH SHARED ROWS THIS STEP MOVES ALONG RATHER THAN ACROSS, and the
      ! orthonormal basis of them, taken once from the iterate and the model
      ! gradient as the active set of the coordinate bounds is.
      !
      ! A row is ACTIVE when the iterate is within
      ! element_constraint_active_band of its own budget AND the model's
      ! steepest-descent direction -g increases c, i.e. would leave the
      ! feasible half-space. That is the published rule for a bound written
      ! for a general linear constraint (Bertsekas 1982, "Projected Newton
      ! methods for optimization problems with simple constraints", SIAM J.
      ! Control Optim. 20, 221, section 3): a constraint the model wants to
      ! move AWAY from stays free and the step moves off it in the same
      ! iteration.
      !
      ! THE HELD COORDINATE BOUNDS ARE REMOVED FROM THE ROW FIRST, so that
      ! the projection cannot give motion back to an unknown the step may
      ! not move: on those slots the row's coefficient is dropped, and the
      ! projection then acts on the free slots alone and commutes with the
      ! hold.
      !
      ! A row whose direction is already spanned by the rows of its own cell
      ! carries no constraint the step does not already respect, and the
      ! Gram-Schmidt below drops it rather than dividing by its round-off.
      ! That is the degeneracy a shared H and O budget can produce, OH and
      ! H2O sitting in both, and it is counted.
      real*8, dimension(nvar_jac*N), intent(inout) :: g
      real*8  :: a(nvar_jac), agv, an, an0, gdotc
      integer :: j, k, q, i, i0, nb
      elem_rows_active_set_frozen = .false.
      if (.not. elem_rows_frozen) return
      n_elem_row_basis = 0
      elem_row_basis   = 0.0d0
      do j = 1, N
         i0 = nvar_jac*(j-1)
         nb = 0
         do k = 1, n_element_constraint_rows
            ! nearly active?
            if (elem_row_c(j,k) .lt.                                      &
                -element_constraint_active_band                           &
                *max(elem_row_budget(j,k), 0.0d0)) cycle
            ! the row over the free slots of this cell
            a = elem_row_grad_scaled(:,j,k)
            if (allocated(species_bound_is_active)) then
               do i = 1, nvar_jac
                  if (species_bound_is_active(i0+i)) a(i) = 0.0d0
               enddo
            endif
            an0 = sqrt(sum(a*a))
            if (an0 .le. 0.0d0) cycle
            ! does the model want to leave the half-space?
            gdotc = 0.0d0
            do i = 1, nvar_jac
               gdotc = gdotc + a(i)*g(i0+i)
            enddo
            if (gdotc .ge. 0.0d0) cycle
            ! orthogonalize against the rows already held for this cell
            do q = 1, nb
               agv = sum(a*elem_row_basis(:,q,j))
               a   = a - agv*elem_row_basis(:,q,j)
            enddo
            an = sqrt(sum(a*a))
            if (an .le. 1.0d-8*an0) then
               n_elem_row_dropped_for_rank =                              &
                  n_elem_row_dropped_for_rank + 1
               cycle
            endif
            nb = nb + 1
            elem_row_basis(:,nb,j) = a/an
         enddo
         n_elem_row_basis(j) = nb
      enddo
      elem_rows_active_set_frozen = .true.
      call project_out_of_the_active_element_constraints(g)
      end subroutine fix_active_element_constraints

      ! ------------------------------------------------------!

      subroutine project_out_of_the_active_element_constraints(s)
      ! THE COMPONENT OF A DIRECTION THAT LEAVES AN ACTIVE SHARED ROW,
      ! removed. What is left moves ALONG the constraint: the whole
      ! direction is kept except the one combination of components that
      ! would spend nuclei the cell does not have, where zeroing a blocked
      ! component would have thrown away the direction's motion in every
      ! other carrier of the same cell as well.
      !
      ! s is in the coordinates the linear algebra works in (the step taken
      ! is D*s), which is where the basis was formed. The basis is
      ! orthonormal and fixed for the step, so this is a fixed linear
      ! idempotent map: the Arnoldi relation the predicted reduction is
      ! taken from is a relation of A composed with this projection, and it
      ! holds because the composition is linear.
      real*8, dimension(nvar_jac*N), intent(inout) :: s
      real*8  :: d
      integer :: j, q, i, i0
      if (.not. elem_rows_active_set_frozen) return
      do j = 1, N
         if (n_elem_row_basis(j) .le. 0) cycle
         i0 = nvar_jac*(j-1)
         do q = 1, n_elem_row_basis(j)
            d = 0.0d0
            do i = 1, nvar_jac
               d = d + elem_row_basis(i,q,j)*s(i0+i)
            enddo
            do i = 1, nvar_jac
               s(i0+i) = s(i0+i) - d*elem_row_basis(i,q,j)
            enddo
         enddo
      enddo
      end subroutine project_out_of_the_active_element_constraints

      ! ------------------------------------------------------!

      subroutine report_element_constraint_rows
      ! WHAT THE SHARED ELEMENT BUDGET IS, printed once where the solve
      ! begins: the range of the budget over the grid and of how tightly the
      ! iterate runs against it. A solve that spends its iterations on the
      ! element budget says what that budget is before it does.
      real*8  :: occmin, occmax, occ
      integer :: j, k, nrow
      if (nspec_row .le. 0) return
      if (.not. carrier_rows_registered()) return
      if (.not. elem_rows_frozen) then
         write(*,'(A)') ' (JFNK) the shared element budget is a face of'//&
              ' the unknown box, not a row of the step'//                 &
              ' (EXHALE_ELEMENT_CONSTRAINT_ROWS=0)'
         return
      endif
      occmin = huge(1.0d0);  occmax = 0.0d0;  nrow = 0
      do j = 1, N
         do k = 1, n_element_constraint_rows
            if (elem_row_budget(j,k) .le. 0.0d0) cycle
            occ = (elem_row_c(j,k) + elem_row_budget(j,k))                &
                 /elem_row_budget(j,k)
            occmin = min(occmin, occ);  occmax = max(occmax, occ)
            nrow = nrow + 1
         enddo
      enddo
      if (nrow .le. 0) return
      write(*,'(A,I0,A,ES10.3,A,ES10.3,A)') ' (JFNK) shared element'//    &
           ' rows of the step: ', nrow, ', occupancy of the budget from ', &
           occmin, ' to ', occmax,                                        &
           ' (the nuclei the carriers of a cell hold over the nuclei'//   &
           ' that element leaves them)'
      end subroutine report_element_constraint_rows

      ! ------------------------------------------------------!

      integer function n_active_element_rows() result(nr)
      ! How many shared rows this step is projected on, over the whole grid.
      ! A step with none of them is a step the shared budget does not reach.
      integer :: j
      nr = 0
      if (.not. elem_rows_active_set_frozen) return
      if (.not. allocated(n_elem_row_basis)) return
      do j = 1, N
         nr = nr + n_elem_row_basis(j)
      enddo
      end function n_active_element_rows

      ! ------------------------------------------------------!

      subroutine release_the_active_element_constraints
      ! No constraint is held outside the step that fixed it.
      elem_rows_active_set_frozen = .false.
      if (allocated(n_elem_row_basis)) n_elem_row_basis = 0
      end subroutine release_the_active_element_constraints

      ! ------------------------------------------------------!

      double precision function                                          &
               largest_step_inside_the_element_constraints(Y, v)         &
               result(tmax)
      ! THE LARGEST MULTIPLE OF A DIRECTION THAT KEEPS EVERY SHARED ROW
      ! SATISFIED, to first order in the step: the fraction-to-the-boundary
      ! rule of the constraint rows, beside the one the coordinate box has.
      !
      ! c is affine in the carrier densities, so c(Y + t v) = c + t (a . v)
      ! is exact in the density unknown and first order in the logarithmic
      ! one, and the step to the surface is -c/(a . v) where a . v is
      ! positive.
      !
      ! A ROW WHOSE ITERATE ALREADY STANDS OUTSIDE ITS OWN BUDGET LIMITS
      ! NOTHING. The budget is frozen at the iterate and depends on the very
      ! unknowns it constrains, so the iterate itself can stand above it;
      ! a rule that returned zero there would refuse every neighborhood of
      ! the point the step starts from. It is the raised target of
      ! elem_row_target, applied to a step length.
      !
      ! huge(1) where no row is formed, which is every solve without a
      ! carrier unknown.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      real*8, dimension(nvar_jac*N), intent(in) :: v
      real*8  :: av, room
      integer :: j, k, i, i0
      tmax = huge(1.0d0)
      if (.not. elem_rows_frozen) return
      do j = 1, N
         i0 = nvar_jac*(j-1)
         do k = 1, n_element_constraint_rows
            av = 0.0d0
            do i = 1, nvar_jac
               av = av + elem_row_grad(i,j,k)*v(i0+i)
            enddo
            if (av .le. 0.0d0) cycle
            room = elem_row_target(j,k) - elem_row_budget(j,k)            &
                 - elem_row_c(j,k)
            room = max(room, 0.0d0)
            if (room .le. 0.0d0) cycle
            tmax = min(tmax, room/av)
         enddo
      enddo
      end function largest_step_inside_the_element_constraints

      ! ------------------------------------------------------!

      subroutine restore_the_element_budget(Yv, project, ncell, nrow,     &
                                            worst_before, worst_after)
      ! A CANDIDATE THAT SPENDS NUCLEI THE CELL DOES NOT HAVE IS STEPPED
      ! BACK ONTO THE CONSTRAINT BEFORE ITS MERIT IS JUDGED (item N4b
      ! deliverable 3): the restoration phase of the constrained step.
      !
      ! The measure is the map's own violation magnitude, relative to what
      ! the row allows, and the step is the one the marching path's limiter
      ! already takes on the same constraint (limit_to_element_budget): the
      ! carriers holding the element are scaled by the one factor that puts
      ! their demand on the surface. It is a step toward feasibility whose
      ! measure falls to zero, and it is one rule for the two paths.
      !
      ! THE ORDER IS CARBON, OXYGEN, HYDROGEN, and it is not arbitrary:
      ! scaling CO down lowers the oxygen and the hydrogen demand, scaling
      ! OH and H2O down lowers the hydrogen demand, and scaling H2 and H+
      ! lowers neither of the others. One pass in that order therefore
      ! leaves every row satisfied, where the reverse order would undo
      ! itself.
      !
      ! WHAT IT RESTORES TO is elem_row_target, the budget or the iterate's
      ! own demand where the iterate stands above it, for the reason stated
      ! at that array: the point the step starts from stays feasible.
      !
      ! With project false nothing is written and the violation is only
      ! measured, which is what the minus ray of the model test needs.
      real*8, dimension(nvar_jac*N), intent(inout) :: Yv
      logical, intent(in)  :: project
      integer, intent(out) :: ncell, nrow
      real*8,  intent(out) :: worst_before, worst_after
      real*8  :: demand(n_element_constraint_rows)
      real*8  :: nc(n_species_row_max), r, v
      integer :: j, k, i, ie, is, korder(n_element_constraint_rows)
      logical :: moved
      ncell = 0;  nrow = 0
      worst_before = 0.0d0;  worst_after = 0.0d0
      if (.not. elem_rows_frozen) return
      ! carbon, oxygen, hydrogen
      korder = [ 3, 2, 1 ]
      do j = 1, N
         ! the carrier densities of this cell, in cm^-3
         nc = 0.0d0
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            is = nvar_jac*(j-1) + 3 + i
            nc(i) = carrier_density_from_unknown(Yv(is), j)*n0
         enddo
         call element_demand_of_the_carriers(nc, demand)
         moved = .false.
         do k = 1, n_element_constraint_rows
            v = relative_violation_of_the_row(demand(k),                  &
                                              elem_row_target(j,k))
            worst_before = max(worst_before, v)
         enddo
         do i = 1, n_element_constraint_rows
            k  = korder(i)
            ie = element_constraint_element(k)
            v  = relative_violation_of_the_row(demand(k),                 &
                                               elem_row_target(j,k))
            if (v .le. element_constraint_restoration_trigger) cycle
            nrow = nrow + 1
            if (.not. project) cycle
            r = elem_row_target(j,k)/max(demand(k), 1.0d-300)
            r = min(max(r, 0.0d0), 1.0d0)
            call scale_the_carriers_of_the_element(j, ie, r, nc, Yv)
            call element_demand_of_the_carriers(nc, demand)
            moved = .true.
         enddo
         if (moved) ncell = ncell + 1
         do k = 1, n_element_constraint_rows
            v = relative_violation_of_the_row(demand(k),                  &
                                              elem_row_target(j,k))
            worst_after = max(worst_after, v)
         enddo
      enddo
      if (project .and. ncell .gt. 0) then
         n_elem_restoration_steps = n_elem_restoration_steps + 1
         n_elem_restoration_cells = n_elem_restoration_cells + ncell
         elem_restoration_worst_before =                                  &
            max(elem_restoration_worst_before, worst_before)
         elem_restoration_worst_after =                                   &
            max(elem_restoration_worst_after, worst_after)
      endif
      end subroutine restore_the_element_budget

      ! ------------------------------------------------------!

      subroutine element_demand_of_the_carriers(nc, demand)
      ! The nuclei of each constrained element the carrier unknowns of one
      ! cell hold [cm^-3], from the densities nc of the carrier rows and the
      ! nuclei table of the map.
      real*8, intent(in)  :: nc(n_species_row_max)
      real*8, intent(out) :: demand(n_element_constraint_rows)
      integer :: i, k, ie
      demand = 0.0d0
      do k = 1, n_element_constraint_rows
         ie = element_constraint_element(k)
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            demand(k) = demand(k)                                         &
               + dble(elem_nu_of_species(srow_isp(i),ie))*nc(i)
         enddo
      enddo
      end subroutine element_demand_of_the_carriers

      ! ------------------------------------------------------!

      double precision function relative_violation_of_the_row(demand,     &
                               allowed) result(v)
      ! How far a demand stands above what a row allows, as a fraction of
      ! it. A demand for nuclei of an element the cell has none of is
      ! infeasible and has no finite fraction, which is the map's own
      ! convention (relative_breach).
      real*8, intent(in) :: demand, allowed
      v = 0.0d0
      if (allowed .gt. 0.0d0) then
         v = max(demand - allowed, 0.0d0)/allowed
      else if (demand .gt. 0.0d0) then
         v = huge(1.0d0)
      endif
      end function relative_violation_of_the_row

      ! ------------------------------------------------------!

      subroutine scale_the_carriers_of_the_element(j, ie, r, nc, Yv)
      ! Every carrier of cell j that holds element ie scaled by r, in the
      ! unknown and in the density array together, so that the caller's
      ! demand and the vector it will evaluate are one state. In the
      ! logarithmic unknown a scaling is a shift by ln r, and a factor of
      ! zero is the floor of that unknown's own space.
      integer, intent(in)    :: j, ie
      real*8,  intent(in)    :: r
      real*8,  intent(inout) :: nc(n_species_row_max)
      real*8, dimension(nvar_jac*N), intent(inout) :: Yv
      integer :: i, is
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         if (elem_nu_of_species(srow_isp(i),ie) .le. 0) cycle
         is    = nvar_jac*(j-1) + 3 + i
         nc(i) = nc(i)*r
         Yv(is) = carrier_unknown_from_density(nc(i)/n0, j, i)
      enddo
      end subroutine scale_the_carriers_of_the_element

      ! ------------------------------------------------------!

      subroutine element_budget_violation_of_composition(u, f_sp, worst,  &
                                        jworst, kworst, ncell)
      ! THE NONLINEAR SHARED CONSTRAINT AT A CANDIDATE, evaluated on the
      ! composition that candidate actually formed and not on the linearized
      ! row (item N4b deliverable 3).
      !
      ! The demand is the map's own sum over the species vector the
      ! write-back wrote, so it carries the write-back's rescaling of the
      ! closure species and the sweep's repartition of the frozen stages; the
      ! budget is the one the step was built on. The magnitude is the
      ! verdict and the cell count is reported beside it, never instead of
      ! it: a state one cell outside a budget by a factor of two and a state
      ! twenty cells outside it by a part in 1e12 are not the same state.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)  :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8,  intent(out) :: worst
      integer, intent(out) :: jworst, kworst, ncell
      type(element_inventory_cell) :: inv
      real*8  :: v
      integer :: j, k, ie
      worst = 0.0d0;  jworst = 0;  kworst = 0;  ncell = 0
      if (.not. elem_rows_frozen) return
      do j = 1, N
         call element_inventory_of_cell(max(u(1,j), 0.0d0)*n0, f_sp(j,:), &
                                        elem_species_is_carried, inv)
         do k = 1, n_element_constraint_rows
            ie = element_constraint_element(k)
            v  = relative_violation_of_the_row(inv%n_carried(ie),         &
                                               elem_row_target(j,k))
            if (v .gt. inventory_feasible_slack) ncell = ncell + 1
            if (v .gt. worst) then
               worst = v;  jworst = j;  kworst = k
            endif
         enddo
      enddo
      end subroutine element_budget_violation_of_composition

      ! ------------------------------------------------------!

      subroutine report_species_unknown_box
      ! THE UPPER FACE OF THE CARRIER BOX AS A DENSITY, printed once where
      ! the solve begins: the box is what a trial is written onto, so a run
      ! says what that face is before it spends its iterations on it.  The
      ! range over the column is what is printed, the face being one number
      ! per cell and carrier.  The SHARED element budget is not this face;
      ! it is a row of the step and report_element_constraint_rows says what
      ! it is.
      real*8  :: himin, himax, nc
      integer :: j, i, is, ncar
      if (nspec_row .le. 0) return
      if (.not. allocated(species_box_hi)) return
      if (.not. carrier_rows_registered()) return
      himin = huge(1.0d0);  himax = 0.0d0;  ncar = 0
      do j = 1, N
         do i = 1, nspec_row
            if (srow_kind(i) .ne. srow_carrier) cycle
            is = nvar_jac*(j-1) + 3 + i
            nc = carrier_density_from_unknown(species_box_hi(is), j)
            himin = min(himin, nc);  himax = max(himax, nc)
            ncar  = ncar + 1
         enddo
      enddo
      if (ncar .le. 0) return
      if (element_constraint_rows_on) then
         write(*,'(A,ES11.3,A,ES11.3,A)') ' (JFNK) carrier ceiling ('//   &
              'code density units): smallest ', himin, ', largest ',      &
              himax, ' -- the density of the cell the carrier sits in'
      else
         write(*,'(A,ES11.3,A,ES11.3,A)') ' (JFNK) carrier ceiling ('//   &
              'code density units): smallest ', himin, ', largest ',      &
              himax, ' -- the element budget of the cell, or its own'//   &
              ' density where that is smaller'
      endif
      end subroutine report_species_unknown_box

      ! ------------------------------------------------------!

      subroutine species_unknowns_outside_their_bounds(Yv, nout, project)
      ! THE ADMISSIBLE SET OF A SPECIES UNKNOWN IS A BOX, NOT A BALL, and a
      ! trust region cannot express it: a radius bounds a norm and this
      ! obstruction is a face. The faces are freeze_species_unknown_box's,
      ! fixed for the step; this counts the unknowns of Yv outside them and,
      ! when asked, writes them onto the nearest one.
      !
      ! WHY THE STEP IS PROJECTED RATHER THAN SHORTENED. An unknown sitting
      ! ON its bound is driven outside by every step whose component there
      ! points out of the box, however short: shortening the WHOLE step by
      ! theta divides every component by theta and leaves that one component
      ! outside. MEASURED on the coupled `mol_carrier` reload before this
      ! projection existed: n(H2) of cell 500 is zero, and the dogleg cut
      ! its trial back by 2^-25 to 2^-49 against that one unknown at
      ! iterations 18 to 23, taking steps of ||s|| = 2.4e-07 against a
      ! radius of 5.0e-04, so 49 iterations moved ||R|| from 1.873 to 1.873
      ! (report B5g section 6.5). Projecting keeps the other 1999 components
      ! of the step at their full length and gives up only the component
      ! that has nowhere to go.
      !
      ! The projection is the same statement write_species_rows_into_
      ! composition makes about the composition it forms, made about the
      ! unknown itself, so the two cannot disagree about which state is
      ! being evaluated.
      real*8, dimension(nvar_jac*N), intent(inout) :: Yv
      integer,                       intent(out)   :: nout
      logical,                       intent(in)    :: project
      real*8  :: lo, hi, val
      integer :: j, i, is
      nout = 0
      if (nspec_row .le. 0) return
      if (.not. allocated(species_box_lo)) return
      do j = 1, N
         do i = 1, nspec_row
            is  = nvar_jac*(j-1) + 3 + i
            val = Yv(is)
            lo  = species_box_lo(is)
            hi  = species_box_hi(is)
            if (val .lt. lo .or. val .gt. hi) then
               nout = nout + 1
               if (project) Yv(is) = min(max(val, lo), hi)
            endif
         enddo
      enddo
      end subroutine species_unknowns_outside_their_bounds

      ! ------------------------------------------------------!

      subroutine hold_the_step_where_it_leaves_the_species_box(Y, Yv,     &
                                                               nheld)
      ! WHERE THE STEP CARRIES A SPECIES UNKNOWN OUT OF ITS BOX, THE STEP
      ! AT THAT UNKNOWN IS ZERO: the unknown keeps the value the iterate
      ! has, which is inside the box because the iterate is.
      !
      ! WHY NOT THE FACE. A component ON a face has no admissible sample of
      ! the operator on either side, so every product the model is built
      ! from drops it (jv_product: hold_the_active_bounds_of and
      ! zero_the_blocked_components_of). A trial written onto the face
      ! therefore carries a move of that unknown that NO model of the step
      ! contains, and the reduction ratio then compares the promise made
      ! about one step with the reduction of another.
      !
      ! MEASURED on the atomic element reload (item N22), outer iterations
      ! 140 to 160 of the entry text: exactly the trials with one unknown
      ! written onto its face are refused (11 of 11) and exactly the trials
      ! with none are accepted at a reduction ratio of 1.000 (10 of 10).
      ! On the refused trials the projection lengthens the step by four in
      ! the SCALED norm and leaves the step in the unknowns themselves,
      ! ||D s||, unchanged to 0.4 percent (5.0146e-10 against 4.9952e-10 at
      ! iterations 150 and 151), so the move is entirely in unknowns whose
      ! column scale is negligible: the image A s and the predicted decrease
      ! do not move at all (relative gap to a fresh product exactly zero)
      ! while the true merit RISES by a fixed 2.9e-3 to 3.6e-3. Holding the
      ! step there makes the trial the step the model is of, and the
      ! iterate stays feasible for the same reason the face did.
      !
      ! The unknown is given up for this step either way; what changes is
      ! that it is given up in the model as well as in the state.
      real*8, dimension(nvar_jac*N), intent(in)    :: Y
      real*8, dimension(nvar_jac*N), intent(inout) :: Yv
      integer,                       intent(out)   :: nheld
      real*8  :: lo, hi, val
      integer :: j, i, is
      nheld = 0
      if (nspec_row .le. 0) return
      if (.not. allocated(species_box_lo)) return
      do j = 1, N
         do i = 1, nspec_row
            is  = nvar_jac*(j-1) + 3 + i
            val = Yv(is)
            lo  = species_box_lo(is)
            hi  = species_box_hi(is)
            if (val .lt. lo .or. val .gt. hi) then
               nheld = nheld + 1
               Yv(is) = Y(is)
            endif
         enddo
      enddo
      end subroutine hold_the_step_where_it_leaves_the_species_box

      ! ------------------------------------------------------!

      subroutine species_box_set_for_test(nrow, lo, hi)
      ! TEST ONLY: install a species-unknown box of the shape the suite
      ! states, so that the box rules can be driven without a state to
      ! freeze one from (freeze_species_unknown_box).
      integer,                       intent(in) :: nrow
      real*8, dimension(nvar_jac*N), intent(in) :: lo, hi
      nspec_row = nrow
      if (allocated(species_box_lo)) deallocate(species_box_lo)
      if (allocated(species_box_hi)) deallocate(species_box_hi)
      allocate(species_box_lo(nvar_jac*N), species_box_hi(nvar_jac*N))
      species_box_lo = lo
      species_box_hi = hi
      end subroutine species_box_set_for_test

      ! ------------------------------------------------------!

      subroutine species_box_clear_for_test
      ! TEST ONLY: the empty registry the three-unknown route runs with.
      nspec_row = 0
      if (allocated(species_box_lo)) deallocate(species_box_lo)
      if (allocated(species_box_hi)) deallocate(species_box_hi)
      end subroutine species_box_clear_for_test

      ! ------------------------------------------------------!

      subroutine scan_carrier_row_across_its_lower_bound(Y, f_sp, u, iter)
      ! WHERE THE ROOT OF A CARRIER ROW SITS RELATIVE TO THE BOUND OF ITS
      ! OWN UNKNOWN, measured at one cell.
      !
      ! A carrier density is a non-negative unknown. At the outermost cell of
      ! a wind whose molecules are long gone its value is zero, i.e. on that
      ! bound, and then every finite-difference sample along a direction with
      ! a negative component there leaves the admissible set however short
      ! the step. Two different things can put it there and they ask for
      ! different repairs: the row's own root can lie at a NEGATIVE value, in
      ! which case the discretization is not positivity-preserving at that
      ! cell and no step control can help; or the root can lie at a small
      ! POSITIVE value that the step keeps overshooting, in which case the
      ! unknown space and the step control are what have to change.
      !
      ! WHAT IS MEASURED. The carrier row of the cell is the entry of F at
      ! that unknown's own slot, so the scan sets the slot to a ladder of
      ! values across the bound, evaluates the full residual at each, and
      ! reports that entry raw and divided by the row's own scale, together
      ! with the merit and whether the state was describable at all. The row
      ! is linear in its own unknown wherever the destruction is first order,
      ! so two values on one side give its root.
      !
      ! EXHALE_CARRIER_BOUND_SCAN=<cell> selects the cell and STOPS the run
      ! once the ladder is reported: it is a measurement of one row and there
      ! is nothing after it to measure.
      real*8, dimension(nvar_jac*N),           intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),          intent(in) :: u
      ! Which outer iterate this call sits at; 0 is the hand-off state,
      ! before the loop. EXHALE_CARRIER_BOUND_SCAN_IT selects one, so that
      ! the line can be measured at an iterate that has ALREADY driven the
      ! unknown onto its bound and not only at the state the solve starts
      ! from.
      integer,                                 intent(in) :: iter
      real*8, dimension(nvar_jac*N) :: Ys, Fs, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool
      real*8, dimension(nvar_jac*N) :: F0
      real*8  :: step, val, f2, rowsc
      integer :: jcell, i, islot, k, ic, iter_want
      logical :: okres
      character(len=32) :: env
      ! The ladder, as MASS FRACTIONS of the cell's own density: a carrier
      ! density in a cell is between the whole of it and nothing, so the
      ! fraction is the variable that spans the range with no scale from the
      ! iterate in it. Thirteen decades down to the bound and one step past
      ! it, plus the iterate itself as the first row.
      !
      ! BELOW THE BOUND the composition the sweep is handed is the clamped
      ! one (write_species_rows_into_composition), so the row there reads
      ! what it reads at the bound; the rows are printed so that this is
      ! visible rather than assumed.
      real*8, parameter :: ladder(12) = (/ 1.0d-2, 1.0d-4, 1.0d-6,       &
              1.0d-8, 1.0d-10, 1.0d-12, 1.0d-14, 1.0d-16, 1.0d-18,       &
              0.0d0, -1.0d-16, -1.0d-6 /)

      call get_environment_variable('EXHALE_CARRIER_BOUND_SCAN', env)
      if (len_trim(env) .le. 0) return
      read(env,*) jcell
      if (jcell .lt. 1 .or. jcell .gt. N) return
      if (nspec_row .le. 0) return
      iter_want = 0
      call get_environment_variable('EXHALE_CARRIER_BOUND_SCAN_IT', env)
      if (len_trim(env) .gt. 0) read(env,*) iter_want
      if (iter .ne. iter_want) return
      write(*,'(A,I0)') ' (carrier_bound_scan) at outer iterate ', iter

      ! THE SCALES COME FROM AN ASSEMBLY, so one has to have happened: the
      ! carrier column and row scales are published by the operator that
      ! forms the rows, and before the first evaluation of a solve they do
      ! not exist and read as their floors.
      call eval_residual(Y, f_sp, fwork, F0, heat, cool)
      call cell_state_scales(Y, D)
      call cell_row_scales(Y, D, Drow, u)
      do i = 1, nspec_row
         if (srow_kind(i) .ne. srow_carrier) cycle
         ic    = srow_idx(i)
         islot = nvar_jac*(jcell-1) + 3 + i
         step  = max(u(1,jcell), 0.0d0)
         rowsc = Drow(islot)
         write(*,'(A,A,A,I0,A,ES11.3,A,ES11.3,A,ES11.3)')                &
              ' (carrier_bound_scan) ', trim(carrier_name(ic)),           &
              ' of cell ', jcell, ': iterate ', Y(islot),                 &
              ', cell density ', step, ', row scale ', rowsc
         write(*,'(A)') ' (carrier_bound_scan)     n(code)'//             &
              '     row F        F/rowscale    merit      describable'
         do k = 0, size(ladder)
            Ys       = Y
            val      = Y(islot)
            if (k .ge. 1) val = ladder(k)*step
            Ys(islot) = val
            call eval_residual(Ys, f_sp, fwork, Fs, heat, cool,           &
                               admissible=okres, may_be_adopted=.false.,  &
                               state_is_discarded=.true.)
            f2 = sqrt(sum((Fs/Drow)**2))
            write(*,'(A,ES13.5,ES13.5,ES14.5,ES12.4,L4)')                 &
                 ' (carrier_bound_scan) ', val, Fs(islot),                &
                 Fs(islot)/max(rowsc, 1.0d-300), f2, okres
         enddo
      enddo
      stop
      end subroutine scan_carrier_row_across_its_lower_bound

      ! ------------------------------------------------------!

      subroutine write_flag_is_about_the_returned_state(tag, rep, iworst, &
                                                        none_measured)
      ! Why a solve whose stop test was met at the loop top does not end
      ! info = 0. The stop test reads the acceptance gate AND the certified
      ! rows of the ITERATE, whose composition is the one the previous
      ! accepted trial was evaluated with; the state handed back carries the
      ! composition its own equilibrium sweep left, and its rows there are
      ! the ones below. That change of composition is the only difference
      ! between the two measurements, and it is what this flag reports.
      character(len=*),  intent(in) :: tag
      type(cert_report), intent(in) :: rep
      integer,           intent(in) :: iworst
      logical,           intent(in) :: none_measured
      write(*,'(A,A)') tag, ' info=2: the stop test -- the acceptance'//  &
           ' gate and every certified row -- was met at the'
      write(*,'(A,A)') tag, '         loop top, on the iterate under the'//&
           ' composition of the previous evaluation; the state'
      write(*,'(A,A)') tag, '         handed back, measured at its own'// &
           ' composition, does not meet it. The state is handed'
      write(*,'(A,A)') tag, '         back and is NOT certified.'
      if (iworst .gt. 0)                                                  &
         write(*,'(A,A,A,A,ES10.3,A,ES8.1,A,I0)') tag, '         refusing'//&
              ' row: ', trim(rep%e(iworst)%name), ' measure ',            &
              rep%e(iworst)%row_max, ' above ', rep%e(iworst)%tol,        &
              ' at cell ', rep%e(iworst)%jworst
      if (none_measured)                                                  &
         write(*,'(A,A)') tag, '         a hydrodynamic row could not be'//&
              ' measured on the state handed back'
      end subroutine write_flag_is_about_the_returned_state

      ! ------------------------------------------------------!

      subroutine solve_steady_ptc(u, f_sp, resid_tol, maxit, dtau0, info)
      ! NOT AVAILABLE WITH THE CARRIER UNKNOWN. This route builds its
      ! Jacobian from frozen_residual, which by construction does no
      ! equilibrium sweep and therefore has no carrier row; giving it one
      ! would mean assembling that row from a background nothing in this
      ! routine refreshes. The coupled solve is the JFNK route, and the
      ! caller is told so rather than being given a silently 3-row answer.
      ! Pseudo-transient-continuation inexact Newton solve of the steady
      ! residual F(Y)=0. Per iteration:
      !   F = eval_residual(Y)                       (full, nonlocal radiation)
      !   J = build_banded_jac (frozen radiation)    (banded approx / precond)
      !   (I/dtau + J) dY = -F   via dgbtrf/dgbtrs    (LAPACK banded LU)
      !   Y <- Y + lam*dY        (backtracking line search + positivity)
      !   dtau <- dtau * rnorm_old/rnorm_new (SER ramp; cut on failure)
      ! As dtau->inf this is Newton; small dtau behaves like explicit
      ! relaxation, giving the robust startup PTC is designed for.
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j
      real*8  :: rnorm, rc(3), dtau, lam, f2, f2_try
      real*8  :: fspread
      integer :: n_jac_unresolved
      logical :: ok, try_ok
      ! The certification of the state handed back, read back so that the
      ! completion flag and the printed refusal are the same record.
      type(cert_report) :: cert_out
      logical :: rows_ok, rows_unmeasured, gate_is_cert
      integer :: i_row
      ! How many rows of the certification the stop test read at the iterate
      ! it stopped on, and the outer iterations that met the acceptance gate
      ! while one of those rows refused.
      integer :: n_rows_at_stop, n_gate_met_row_refused, i_row_refused_last

      if (nspec_row .gt. 0) then
         write(*,'(A)') ' (PTC) the carrier unknown is a JFNK-only'//&
              ' route; this solve is not attempted'
         info = 2;  return
      endif
      neq  = nvar_jac*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      dtau = dtau0
      info = 1
      n_rows_at_stop = 0;  n_gate_met_row_refused = 0
      i_row_refused_last = -1
      call reset_no_chem_root_trial_count
      call read_composition_elimination_controls
      call read_species_unknown_space_controls

      call pack_U(u, Y)
      ! The sweeps below evaluate states the SOLVER invents, not states the
      ! run holds; tag them so that the acceptance ledgers, the non-root
      ! streak and its stop keep the two apart (section 121).
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      ! Seed from the adopted composition, write the swept one into the work
      ! array, then ADOPT it. The copy is an adoption and not a seeding
      ! convention: the seed is named in the call.
      call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0)
      f_sp = f_sp_j
      n_eq_sweeps_model = n_eq_sweeps_last
      n_no_chem_root_state = n_no_chem_root_last
      call keep_background_of_adopted_state
      call resid_relnorm(F, u, rc, rnorm)
      f2 = sqrt(sum(F*F))                ! smooth line-search merit
      write(*,'(A,ES11.3,A,3ES10.2,A,ES10.2)') ' (PTC) start ||R||=',rnorm, &
           '  R(m,p,E)=',rc,'  ||F||2=',f2
      call write_resid_below_escape(' (PTC)', F, u)

      do iter = 1, maxit
         ! THE STOP TEST OF THE ITERATION, WHICH IS THE ACCEPTANCE TEST OF
         ! THE STATE: the acceptance gate -- residual, flux spread and the
         ! chemistry, n_no_chem_root_state belonging to the iterate's own
         ! sweep, so a state resting on a cell without a chemical root
         ! cannot be declared solved here -- AND every row the
         ! certification judges this system by, each against its own
         ! tolerance. The gate alone is not the acceptance (D3,
         ! docs/solver_partition_experiment_20260911.md section 7.3): it
         ! reads one number against the run's "Resid tol" while the mass
         ! row is certified against 3e-12. When the gate is met and a
         ! certified row is not, the state is not a solution of the
         ! equations this solve was given and the iteration goes on.
         if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nspec_row .gt. 0, carrier_relnorm_last,   &
                          n_no_chem_root_state)) then
            call certified_rows_of_the_state(F, u, f_sp, resid_tol,      &
                     n_no_chem_root_state, rows_ok, i_row, cert_out,     &
                     n_rows_at_stop)
            if (rows_ok) then
               info = 0
               ! The report of the state declared solved, which is the
               ! record the acceptance below reads: this route hands back
               ! the iterate it stopped on, so the two are one state.
               call certification_report_write(cert_out,                 &
                    '(PTC) state declared solved and handed back')
               write(*,'(A,I0,A)') ' (PTC) loop-top stop: the'//         &
                    ' acceptance gate is met and every one of the ',     &
                    n_rows_at_stop, ' certified row(s) of the system'//  &
                    ' this solve carries is within its own tolerance'
               exit
            endif
            n_gate_met_row_refused = n_gate_met_row_refused + 1
            if (i_row .ne. i_row_refused_last) then
               call write_gate_met_while_a_row_refuses(' (PTC)',         &
                        cert_out, i_row, iter)
               i_row_refused_last = i_row
            endif
         endif

         ! Full-residual banded Jacobian (includes the local source
         ! derivative d(heat-cool)/dE), then the PTC system (I/dtau + J).
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
         call build_banded_jac_full(Y, f_sp, ab, n_jac_unresolved)
         if (n_jac_unresolved .gt. 0)                                    &
            write(*,'(A,I0,A,I0,A)') ' (PTC) Jacobian: ',                &
                 n_jac_unresolved, ' of ', ncolor_jac,                   &
                 ' color(s) had no admissible probe; their columns are'//&
                 ' zero'
         abf = ab
         do jc = 1, neq
            abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + 1.0d0/dtau
         enddo
         dY = -F
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         if (lpinfo .ne. 0) then
            write(*,'(A,I0,A)') ' (PTC) dgbtrf info=',lpinfo,            &
                 ' (singular); cutting dtau'
            dtau = max(dtau*0.25d0, dtau0*1.0d-3);  cycle
         endif
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, ldab, ipiv,      &
                     dY, neq, lpinfo)

         ! Backtracking line search on the SMOOTH merit ||F||_2 (the
         ! Newton direction reduces this; the component-wise max-relative
         ! rnorm is non-smooth and rejected valid steps). Armijo condition
         ! with positivity (rho>0, p>0).
         lam = 1.0d0;  ok = .false.
         do ls = 1, 20
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)   ! ghosts of utry are not set
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok)
               f2_try = sqrt(sum(Ftry*Ftry))
               ! A trial the code cannot describe is not a descent
               ! candidate whatever its merit evaluates to; it is rejected
               ! and the step is halved, exactly as a negative density is.
               if (try_ok .and. f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2) then
                  ok = .true.;  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            n_no_chem_root_state = n_no_chem_root_last
            n_eq_sweeps_model = max(n_eq_sweeps_model, n_eq_sweeps_last)
            ! The last residual evaluation was this accepted trial, so the
            ! frozen background now describes the state just adopted.
            call keep_background_of_adopted_state
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            ! SER ramp on the merit ratio, scaled by the accepted step.
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            ! No descent found: shrink the pseudo-time step (more
            ! relaxation-like) but keep it >= the explicit-stable dtau0.
            dtau = max(dtau*0.25d0, dtau0)
         endif

         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES10.2,A,ES9.2)') ' (PTC) it', &
              iter,'  ||R||=',rnorm,'  ||F||2=',f2,'  dtau=',dtau,'  lam=',lam
      enddo

      call unpack_U(Y, u);  call Apply_BC(u)
      call resid_relnorm(F, u, rc, rnorm)
      call set_ioniz_eq_sweep_state_kind(ieq_state_marching)
      call install_background_of_adopted_state
      ! THE COMPLETION FLAG IS A STATEMENT ABOUT THE STATE HANDED BACK. The
      ! gate that set info = 0 was taken at the loop top, on the iterate;
      ! this asks the same question of the state that leaves the routine,
      ! with its rows read off the certification made of it. A route that
      ! meets the gate at the loop top and not here has not produced a
      ! stationary state and says info = 2.
      if (info .eq. 0) then
         cert_out = certification_last_report()
         call stationary_rows_of_the_returned_state(cert_out, rows_ok,   &
                                                 i_row, rows_unmeasured)
         call gate_equals_certification(cert_out, rnorm, ' (PTC)',       &
                                        gate_is_cert)
         if (.not. (rows_ok .and. gate_is_cert .and.                     &
             steady_gates_met(rnorm, u, resid_tol, fspread, .false.,     &
                              0.0d0, n_no_chem_root_state))) then
            info = 2
            call write_flag_is_about_the_returned_state(' (PTC)',        &
                              cert_out, i_row, rows_unmeasured)
         endif
      endif
      write(*,'(A,I0,A,ES11.3)') ' (PTC) done info=',info,' ||R||=',rnorm
      ! Iterates at which the acceptance gate was met while a row of the
      ! certification refused: the iterations this solve spent on a row the
      ! gate does not see.
      if (n_gate_met_row_refused .gt. 0)                                  &
         write(*,'(A,I0)') ' (PTC) outer iterations that met the'//       &
              ' acceptance gate while a certified row refused: ',         &
              n_gate_met_row_refused
      if (n_no_chem_root_state .gt. 0)                                    &
         write(*,'(A,I0,A)') ' (PTC) gate NOT met: ',                     &
              n_no_chem_root_state, ' cell(s) of the state handed back'// &
              ' carry a composition that is not a root of the chemical'// &
              ' network (acceptance class 4 or 6); the state is written'// &
              ' but'// &
              ' NOT certified'
      call write_no_chem_root_trial_count(' (PTC)')
      call write_resid_below_escape(' (PTC)', F, u)

      gate_rnorm_accepted = rnorm;  gate_fspread_accepted = fspread
      deallocate(Y,F,Ftry,dY,Ytry,ab,abf,ipiv)
      end subroutine solve_steady_ptc

      ! ------------------------------------------------------!

      pure real*8 function probe_length_of_the_jacobian_action(Y)       &
                   result(probe_length)
      ! HOW FAR THE MATRIX-FREE PRODUCT MOVES THE STATE, as a length in the
      ! unknowns themselves. The forward difference is taken at
      ! Y + eps0 v with eps0 = probe_length/||v||, so the DISPLACEMENT has
      ! this norm whatever the direction and whatever its length: the
      ! quotient is homogeneous of degree one in v, which is what makes the
      ! action a linear map, and the arc the residual is sampled over is a
      ! property of the iterate alone.
      !
      ! sqrt(epsilon) (1 + ||Y||) is the standard scaled interval: the
      ! rounding of the residual, of relative size epsilon, and the
      ! curvature of the residual over the interval contribute equally to
      ! the error of the quotient at this length.
      !
      ! CONSEQUENCE, MEASURED (item N22, atomic element reload): the arc is
      ! 1.845e-7 with ||Y|| = 11.4, nine decades above the reproducibility
      ! of the residual, so the action is not noise; and at outer iteration
      ! 151 the step itself is ||D s|| = 5.0e-10, which is 2.7e-3 of the
      ! arc, so the model of that step is a secant over an arc the step
      ! does not reach. Shortening the step cannot shorten the arc.
      !
      ! THE ARC IS SCALED BY jv_probe_length_factor, which is exactly one
      ! outside the additivity hook's length scan, so the standard arc is
      ! the expression above unchanged.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      probe_length = sqrt(epsilon(1.0d0))*(1.0d0 + sqrt(sum(Y*Y)))
      if (jv_probe_arc_scale .ne. 1.0d0)                                  &
         probe_length = probe_length*jv_probe_arc_scale
      if (jv_probe_length_factor .ne. 1.0d0)                              &
         probe_length = probe_length*jv_probe_length_factor
      end function probe_length_of_the_jacobian_action

      ! ------------------------------------------------------!

      real*8 function probe_step_on_the_column_scales(Y, v) result(eps0)
      ! THE PROBE STEP THAT DISPLACES EVERY UNKNOWN BY THE SAME FRACTION OF
      ! ITS OWN SCALE.
      !
      !   eps0 = sqrt(epsilon) (1 + ||Dc^-1 Y||) / ||Dc^-1 v||,
      !
      ! Dc the state column scale (state_column_scale, the characteristic
      ! magnitude of each unknown). Written in the coordinates Dc^-1 Y, in
      ! which every unknown is of order one, this is the same rule
      ! probe_length_of_the_jacobian_action states in the unknowns
      ! themselves: the displacement Dc^-1 (eps0 v) has the norm
      ! sqrt(epsilon)(1 + ||Dc^-1 Y||), so a direction spread over the n
      ! unknowns moves each of them by about sqrt(epsilon) of its own
      ! scale, which is the interval at which the truncation of a forward
      ! difference and the rounding of the residual contribute equally.
      !
      ! WHY THE UNSCALED RULE IS NOT THAT. It fixes the length of the
      ! displacement at sqrt(epsilon)(1 + ||Y||) with ||Y|| the norm of the
      ! MIXED vector, so the length is set by whichever unknowns are
      ! largest and every unknown below them is displaced by a smaller
      ! fraction of itself. MEASURED on the atomic element reload at outer
      ! iteration 1 (item N31, EXHALE_JV_ADDITIVITY): ||Y|| = 11.4 is
      ! carried by the innermost cells, whose mass unknown is then
      ! displaced by 2.2e-8 of itself, about sqrt(epsilon); the unknowns of
      ! the outer cells are displaced by 2e-10 to 6e-10 of themselves, one
      ! to two decades short of it, and the difference quotient along the
      ! preconditioned Krylov directions is then rounding and not a
      ! derivative (the additivity defect rises as the reciprocal of the
      ! arc, exponent -1.03, and the second difference of the residual sits
      ! on a floor instead of scaling as the square of the step).
      !
      ! VALIDITY. It needs a column scale that is a magnitude of the
      ! unknown; where none has been formed yet the unscaled rule stands.
      ! It says nothing about a residual that is non-smooth for a reason
      ! other than a displacement too short to resolve: a longer step
      ! raises the truncation error of the forward difference in
      ! proportion, so this trades one error for the other at their
      ! balance point and no further.
      real*8, dimension(nvar_jac*N), intent(in) :: Y, v
      real*8  :: yn, vn
      ! NOT n: that name is the grid size N of global_parameters, and
      ! Fortran does not distinguish the case.
      integer :: nunk
      nunk = nvar_jac*N
      if (.not. allocated(state_column_scale)) then
         eps0 = probe_length_of_the_jacobian_action(Y)                    &
                /max(sqrt(sum(v*v)), 1.0d-300)
         return
      endif
      if (size(state_column_scale) .ne. nunk) then
         eps0 = probe_length_of_the_jacobian_action(Y)                    &
                /max(sqrt(sum(v*v)), 1.0d-300)
         return
      endif
      yn = sqrt(sum((Y/state_column_scale)**2))
      vn = sqrt(sum((v/state_column_scale)**2))
      eps0 = sqrt(epsilon(1.0d0))*(1.0d0 + yn)/max(vn, 1.0d-300)
      if (jv_probe_arc_scale .ne. 1.0d0)                                  &
         eps0 = eps0*jv_probe_arc_scale
      if (jv_probe_length_factor .ne. 1.0d0)                              &
         eps0 = eps0*jv_probe_length_factor
      end function probe_step_on_the_column_scales

      ! ------------------------------------------------------!

      subroutine jv_product(Y, F0, f_sp_base, v, Jv, ok)
      ! Matrix-free Jacobian-vector product J*v by a forward directional
      ! finite difference of the FULL residual (captures the non-local
      ! radiation coupling the banded Jacobian omits):
      !   J*v ~= ( F(Y + eps*v) - F0 ) / eps,   F0 = F(Y).
      ! eps is the standard scaled step. f_sp is restored from the base
      ! copy (eval_residual mutates it).
      real*8, dimension(nvar_jac*N),                  intent(in)  :: Y, F0, v
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)  :: f_sp_base
      real*8, dimension(nvar_jac*N),                  intent(out) :: Jv
      ! .false. when no step along v, down to the smallest of
      ! n_probe_step_halvings, lands on a state the code can describe: then
      ! there is no directional derivative to return and Jv is meaningless
      ! (it is set to zero so that nothing non-finite can escape) -- the
      ! caller must stop using v, not use the value.
      logical,                                 intent(out) :: ok
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(nvar_jac*N) :: Fp
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      real*8, dimension(nvar_jac*N) :: vh, vb
      real*8  :: vn, eps, eps0, step_fwd, step_bwd
      integer :: ish, n_blk_fwd, n_blk_bwd
      logical :: take_backward
      ! How much of the room to the face of the box a probe step may use. A
      ! step to the face itself would put an unknown ON its bound, where the
      ! sweep's own reaction network has no interior to linearize about.
      real*8, parameter :: box_safety = 0.9d0
      Jv = 0.0d0
      ok = .true.
      ! THE UNKNOWNS THE STEP MAY NOT MOVE ARE NOT SAMPLED EITHER. A
      ! direction with a component pointing out of the species box crosses
      ! that face at every step length, so the sample does not exist; the
      ! step control has already declared those components held
      ! (fix_active_species_bounds) and this is the same declaration applied
      ! to the operator it samples, which is what makes the two one map.
      ! Outside a bound-constrained step nothing is held and this is the
      ! identity.
      vh = v
      call hold_the_active_bounds_of(vh)
      vn = sqrt(sum(vh*vh))
      if (vn .le. 0.0d0) return
      if (jv_probe_on_the_column_scales) then
         eps0 = probe_step_on_the_column_scales(Y, vh)
      else
         eps0 = probe_length_of_the_jacobian_action(Y)/vn
      endif
      ! WHAT THIS PROBE DOES, for the additivity hook to read afterwards.
      ! Written by every product and read by nothing that decides.
      jv_probe_step_nominal_last = eps0
      jv_probe_step_last         = eps0
      jv_probe_blocked_last      = 0
      jv_probe_backward_last     = .false.
      ! THE STEP IS CUT TO THE BOX BEFORE IT IS TAKEN, on whichever side of
      ! the direction has more room inside it.
      !
      ! The length above is set by the whole vector; the component of the
      ! direction at one species unknown is set by the banded preconditioner
      ! and can be many decades larger than that unknown's own distance from
      ! its bound. Halving the step until it fits is a SEARCH for the number
      ! largest_step_inside_the_species_box computes directly, and a search
      ! that runs out of halvings loses the whole Krylov cycle: MEASURED on
      ! the coupled `mol_carrier` reload, 12 of 20 outer iterations ended on
      ! "no Krylov direction could be sampled" with the probe putting the
      ! carrier of cell 217 at -2.4e-16.
      !
      ! Which side is a question about the box and not about the physics, so
      ! it is decided by the room: the quotient is first-order accurate in
      ! eps on either side and the longer step is the better difference. This
      ! subsumes the backward sample, which was the same choice made after
      ! every forward halving had already failed.
      ! AND TO THE SHARED ELEMENT ROWS, which is the other half of the
      ! feasible set: a probe that spends nuclei the cell does not have
      ! lands on a state the chemistry has no root in, exactly as a negative
      ! density does (largest_step_inside_the_element_constraints).
      ! BOTH SIDES ARE CAPPED AT THE PROBE STEP BEFORE THEY ARE COMPARED:
      ! where neither side is the binding one the two are the same number
      ! and the side is the forward one, which is what keeps the choice from
      ! turning on a room neither side uses.
      step_fwd = min(eps0,                                                &
             box_safety*largest_step_inside_the_species_box(Y, vh,         &
                                                           n_blk_fwd),    &
             box_safety*largest_step_inside_the_element_constraints(Y, vh))
      vb = -vh
      step_bwd = min(eps0,                                                &
             box_safety*largest_step_inside_the_species_box(Y, vb,         &
                                                           n_blk_bwd),    &
             box_safety*largest_step_inside_the_element_constraints(Y, vb))
      ! WHICH SIDE, and the WHOLE direction is sampled wherever a side
      ! exists along which every component has room. A quotient taken with
      ! one component zeroed is a derivative of the residual composed with
      ! the projection onto the box, which is the map the trial evaluates
      ! but is not linear in the direction, so it is used only where the
      ! alternative is no sample at all: with one component on a face one of
      ! the two sides is normally free, and the choice between two free
      ! sides is the longer step, which is the better difference.
      !
      ! WHERE BOTH SIDES ARE BLOCKED the entry text lost the whole Krylov
      ! cycle -- MEASURED on `mol_diffusion`, 18 cycles truncated -- and the
      ! blocked components of the better side are then zeroed and the rest
      ! of the direction is sampled. The step of a projected trial has its
      ! own Jacobian product either way (trust_region_step).
      if (n_blk_fwd .eq. 0 .or. n_blk_bwd .eq. 0) then
         if (n_blk_fwd .gt. 0) then
            take_backward = .true.
         else if (n_blk_bwd .gt. 0) then
            take_backward = .false.
         else
            take_backward = (step_bwd .gt. step_fwd)
         endif
         if (take_backward) then
            vh   = -vb
            eps0 = min(eps0, step_bwd)
         else
            eps0 = min(eps0, step_fwd)
         endif
      else
         take_backward = (n_blk_bwd .lt. n_blk_fwd) .or.                  &
                         (n_blk_bwd .eq. n_blk_fwd .and.                  &
                          step_bwd .gt. step_fwd)
         if (take_backward) then
            call zero_the_blocked_components_of(Y, vb)
            vh   = -vb
            eps0 = min(eps0, step_bwd)
            n_probe_blocked_components =                                  &
               n_probe_blocked_components + n_blk_bwd
            jv_probe_blocked_last = n_blk_bwd
         else
            call zero_the_blocked_components_of(Y, vh)
            eps0 = min(eps0, step_fwd)
            n_probe_blocked_components =                                  &
               n_probe_blocked_components + n_blk_fwd
            jv_probe_blocked_last = n_blk_fwd
         endif
      endif
      vn = sqrt(sum(vh*vh))
      jv_probe_step_last     = eps0
      jv_probe_backward_last = take_backward
      if (eps0 .le. 0.0d0 .or. vn .le. 0.0d0) then
         ok = .false.;  return
      endif
      if (take_backward) then
         call backward_directional_difference(Y, F0, f_sp_base, vh, eps0,  &
                                              Jv, ok)
         return
      endif
      eps  = eps0
      do ish = 0, n_probe_step_halvings
         ! THE SAME MAP AS F0. F0 was formed at Y from f_sp_base with
         ! n_eq_sweeps_model composition passes, and this sample takes
         ! exactly that many from the same seed, so the difference below is
         ! a directional derivative of one function and not the gap between
         ! a converged elimination and a truncated one.
         call eval_residual(Y + eps*vh, f_sp_base, fwork, Fp, heat, cool,  &
                            admissible=ok,                                 &
                         may_be_adopted=.false.,                           &
                         state_is_discarded=.true.,                        &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
         if (ok) exit
         eps = 0.5d0*eps
      enddo
      if (ok) then
         jv_probe_step_last = eps
         Jv = (Fp - F0)/eps
         return
      endif
      ! THE OTHER SIDE OF THE SAME DERIVATIVE, when no step FORWARD along v
      ! lands on a describable state.
      !
      ! A species unknown sitting AT its lower bound is driven negative by
      ! every positive eps along a direction whose component there is
      ! negative, however small eps is, so shortening the step cannot
      ! recover the sample and the whole Krylov cycle is lost with it.
      ! MEASURED on the coupled `mol_carrier` reload: n(H2) of cell 500 is
      ! zero, the probe puts it at -7.444e-19 of the code density unit, that
      ! is the ONLY species unknown of the state on its bound, the sample is
      ! refused at every halving, GMRES returns no direction, and 18 of the
      ! 49 outer iterations took no step for that reason (report B5g
      ! section 6.3).
      !
      ! [F(Y) - F(Y - eps v)]/eps is a sample of the SAME directional
      ! derivative, first-order accurate in eps exactly as the forward
      ! quotient is, and it steps away from the bound the forward step
      ! crossed. With one unknown on the bound one of the two sides is
      ! always available; with unknowns on their bounds in both directions
      ! neither is, and that is a bound-constrained unknown space rather
      ! than a step-length question.
      !
      ! It is tried ONLY where the forward sample was refused at every step
      ! length, so no case in which the forward sample exists can move.
      call backward_directional_difference(Y, F0, f_sp_base, vh, eps0, Jv,  &
                                           ok)
      end subroutine jv_product

      ! ------------------------------------------------------!

      subroutine backward_directional_difference(Y, F0, f_sp_base, v, eps0, &
                                                 Jv, ok)
      ! [F(Y) - F(Y - eps v)]/eps: a sample of the SAME directional
      ! derivative, first-order accurate in eps exactly as the forward
      ! quotient is, taken on the other side of the direction.
      !
      ! Two callers, one statement. The step control takes this side when
      ! the box leaves more room there (jv_product), and it takes it as the
      ! last resort when no forward step of any length lands on a
      ! describable state -- a species unknown AT its bound is carried out
      ! of the set by every positive step along a direction whose component
      ! there points out, so shortening cannot recover that sample.
      real*8, dimension(nvar_jac*N),          intent(in)  :: Y, F0, v
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp_base
      real*8,                                 intent(in)  :: eps0
      real*8, dimension(nvar_jac*N),          intent(out) :: Jv
      logical,                                intent(out) :: ok
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(nvar_jac*N) :: Fm
      real*8, dimension(1-Ng:N+Ng)  :: heat, cool
      real*8  :: eps
      integer :: ish
      Jv = 0.0d0
      ok = .false.
      if (eps0 .le. 0.0d0) return
      eps = eps0
      do ish = 0, n_probe_step_halvings
         call eval_residual(Y - eps*v, f_sp_base, fwork, Fm, heat, cool,   &
                            admissible=ok,                                 &
                         may_be_adopted=.false.,                           &
                         state_is_discarded=.true.,                        &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
         if (ok) exit
         eps = 0.5d0*eps
      enddo
      if (.not. ok) return
      n_jv_backward_sample = n_jv_backward_sample + 1
      Jv = (F0 - Fm)/eps
      end subroutine backward_directional_difference

      ! ------------------------------------------------------!

      subroutine gm_outcome_text(code, txt)
      ! The sentence that belongs to one outcome of the linear solve.
      integer,          intent(in)  :: code
      character(len=*), intent(out) :: txt
      select case (code)
      case (gm_tolerance_reached)
         txt = 'the requested linear tolerance was reached'
      case (gm_subspace_exhausted)
         txt = 'the subspace was exhausted with a finite approximate step'
      case (gm_no_direction_sampled)
         txt = 'not one direction of the operator could be sampled'
      case (gm_jacobian_action_unusable)
         txt = 'the Jacobian action stopped being usable inside the cycle'
      case (gm_preconditioner_failed)
         txt = 'the preconditioner returned a failure'
      case (gm_reduced_system_singular)
         txt = 'the reduced problem is singular or inconsistent'
      case (gm_nonfinite)
         txt = 'the arithmetic of the cycle left the finite reals'
      case (gm_stopped_on_the_trust_ball)
         txt = 'the iterate left the trust ball and its boundary point'//  &
               ' was returned'
      case (gm_true_residual_turned)
         txt = 'the true residual of the step turned upward and the'//     &
               ' best one seen was returned'
      case default
         txt = 'unnamed'
      end select
      end subroutine gm_outcome_text

      ! ------------------------------------------------------!

      pure logical function every_component_is_finite(v) result(ok)
      ! Whether a vector is an ordinary vector of reals: no NaN and no
      ! infinity anywhere in it. Written out rather than taken from
      ! ieee_arithmetic for the reason given at finite_real.
      real*8, dimension(:), intent(in) :: v
      ok = all(v .eq. v) .and. all(abs(v) .le. huge(1.0d0))
      end function every_component_is_finite

      ! ------------------------------------------------------!

      pure real*8 function scaled_two_norm(v) result(nrm)
      ! ||v||_2 taken about the largest component, so that it is the norm
      ! even where the sum of squares would underflow or overflow. Used
      ! where the value of a small residual is the answer -- the true
      ! residual verified at a breakdown -- and NOT on the Arnoldi
      ! recursion, whose remainder norm is left as it was.
      real*8, dimension(:), intent(in) :: v
      real*8 :: sc
      sc = maxval(abs(v))
      if (sc .le. 0.0d0) then
         nrm = 0.0d0
      else
         nrm = sc*sqrt(sum((v/sc)**2))
      endif
      end function scaled_two_norm

      ! ------------------------------------------------------!

      subroutine linear_operator_set_for_test(A, kind)
      ! Install the dense operator the Krylov tests state, in place of the
      ! Jacobian action and the banded preconditioner.
      real*8, dimension(:,:), intent(in) :: A
      integer,                intent(in) :: kind
      if (allocated(lin_test_A)) deallocate(lin_test_A)
      allocate(lin_test_A(size(A,1),size(A,2)))
      lin_test_A    = A
      lin_test_kind = kind
      lin_test_n_product = 0
      end subroutine linear_operator_set_for_test

      ! ------------------------------------------------------!

      subroutine linear_operator_clear_for_test
      ! Back to the true Jacobian action and the banded factorization.
      if (allocated(lin_test_A)) deallocate(lin_test_A)
      lin_test_kind = lin_test_absent
      end subroutine linear_operator_clear_for_test

      ! ------------------------------------------------------!

      subroutine banded_preconditioner_solve(abf, ipiv, neq, z, info)
      ! M^-1 z with M the factored banded scaled system (I/dtau + D^-1 J_b D),
      ! or the identity a test operator states in its place. info is the
      ! LAPACK status of the triangular solve.
      integer,                                  intent(in)    :: neq
      real*8, dimension(2*kl_jac+ku_jac+1,neq), intent(in)    :: abf
      integer, dimension(neq),                  intent(in)    :: ipiv
      real*8, dimension(neq),                   intent(inout) :: z
      integer,                                  intent(out)   :: info
      if (lin_test_kind .eq. lin_test_absent) then
         ! WHERE THE FACTORED BAND WAS EQUILIBRATED, R abf C is what was
         ! factored, so abf^-1 v = C (R abf C)^-1 R v and the two diagonal
         ! rescalings are undone here. Both factors are powers of two, so
         ! this costs two elementwise multiplies and no rounding.
         if (preconditioner_equilibrated) then
            z = precon_row_scale*z
            call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf,                &
                        2*kl_jac+ku_jac+1, ipiv, z, neq, info)
            z = precon_col_scale*z
            return
         endif
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, 2*kl_jac+ku_jac+1, &
                     ipiv, z, neq, info)
         return
      endif
      info = 0
      if (lin_test_kind .eq. lin_test_dense_preconditioner_fails) info = 1
      end subroutine banded_preconditioner_solve

      ! ------------------------------------------------------!

      subroutine jacobian_action_of_direction(Y, F0, f_sp_base, v, Jv, ok)
      ! J v, the residual's directional derivative along v: the matrix-free
      ! difference of the true residual, or the dense operator a test states
      ! in its place. ok is false where no admissible sample of the true
      ! residual exists along v (jv_product), so that the caller stops using
      ! v rather than using the value.
      real*8, dimension(nvar_jac*N),          intent(in)  :: Y, F0, v
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp_base
      real*8, dimension(nvar_jac*N),          intent(out) :: Jv
      logical,                                intent(out) :: ok
      if (lin_test_kind .eq. lin_test_absent) then
         call jv_product(Y, F0, f_sp_base, v, Jv, ok)
         return
      endif
      lin_test_n_product = lin_test_n_product + 1
      Jv = matmul(lin_test_A, v)
      if (lin_test_kind .eq. lin_test_dense_action_drifts .and.          &
          lin_test_n_product .gt. lin_test_drift_after)                  &
         Jv = Jv + lin_test_drift_per_product                            &
                  *dble(lin_test_n_product - lin_test_drift_after)*v
      ok = .true.
      end subroutine jacobian_action_of_direction

      ! ------------------------------------------------------!

      subroutine deterministic_unit_direction(seed, v)
      ! A REPRODUCIBLE DIRECTION OF UNIT LENGTH. The additivity measurement
      ! needs a pair of directions that are not the Krylov basis's, and it
      ! has to measure the same pair in two runs of one binary, so the
      ! entries come from a linear congruential recurrence of the seed
      ! (the Numerical Recipes ranqd1 multiplier and increment, reduced to
      ! 32 bits) mapped onto [-1,1) and normalized. Nothing here depends on
      ! the quality of the sequence; it is a direction.
      integer,              intent(in)  :: seed
      real*8, dimension(:), intent(out) :: v
      integer*8 :: st
      integer   :: i
      real*8    :: vn
      st = int(seed, 8)
      do i = 1, size(v)
         st   = iand(1664525_8*st + 1013904223_8, 4294967295_8)
         v(i) = 2.0d0*(real(st, 8)/4294967296.0d0) - 1.0d0
      enddo
      vn = sqrt(sum(v*v))
      if (vn .gt. 0.0d0) v = v/vn
      end subroutine deterministic_unit_direction

      ! ------------------------------------------------------!

      subroutine measure_the_additivity_of_the_jacobian_action(Y, F0,     &
                     f_sp_base, D, u1, u2, pair_name)
      ! IS THE FINITE-DIFFERENCE ACTION A LINEAR MAP OF THE DIRECTION?
      !
      ! GMRES, the Arnoldi relation and the predicted decrease the trust
      ! region reads from it are relations of a LINEAR operator. What
      ! jv_product returns is
      !     A(v) = [ F(Y + eps(v) v) - F(Y) ] / eps(v),
      !     eps(v) = L / ||P(v) v||,   L = probe_length_of_the_jacobian_action,
      ! P the hold of the active bounds and eps further cut to the species
      ! box and the element constraints. A(v) is homogeneous of degree one
      ! in v by construction, but additivity is a separate property and it
      ! can fail in exactly three ways:
      !
      !   CURVATURE OF A SMOOTH RESIDUAL. Expanding to second order,
      !   A(v) = J v + (L/2) (v^T H v)/||v|| + O(L^2), so the defect
      !   A(v1+v2) - A(v1) - A(v2) is proportional to L: it FALLS WITH THE
      !   ARC as the arc.
      !
      !   ROUNDING, or a residual that is not smooth on the arc (a branch, a
      !   clamp, a max, a table lookup, an inner iteration stopped on a
      !   tolerance). The quotient carries the residual's own irreducible
      !   spread divided by eps, so the defect RISES AS THE RECIPROCAL OF
      !   THE ARC.
      !
      !   A DECISION OF THE FEASIBLE SET INSIDE jv_product. The side taken
      !   (forward or backward), the components zeroed where both sides are
      !   blocked, and the cut of eps to the box are all functions of the
      !   DIRECTION, so the map sampled along v1 + v2 need not be the map
      !   sampled along v1 and along v2. A forward and a backward quotient
      !   of one direction differ in the SIGN of the curvature term, and a
      !   zeroed component removes a column of the operator outright. This
      !   source is discrete: the defect then does NEITHER of the above.
      !
      ! The three prints below are that separation, and the fourth reads
      ! the residual's own smoothness on the arc directly, by second
      ! differences along the sum direction: a smooth residual gives a
      ! second difference proportional to the square of the step, and an
      ! inner solve stopped on a tolerance gives a floor instead.
      !
      ! TEST ONLY. Every product taken here is discarded, the probe-length
      ! multiplier is restored before returning, and nothing outside the
      ! prints reads any of it.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y, F0, D
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      ! The two directions, in the unknowns themselves (what jv_product is
      ! handed, that is D times a direction of the scaled space).
      real*8, dimension(nvar_jac*N),          intent(in) :: u1, u2
      character(len=*),                       intent(in) :: pair_name
      real*8, dimension(:), allocatable :: a1, a2, as, g, vs, Fh1, Fh2, Fh4
      ! The defect and the two actions AT THE STANDARD ARC, kept while the
      ! length scan overwrites the work vectors with the other two arcs.
      real*8, dimension(:), allocatable :: gs, sum_of_the_two_actions
      real*8, dimension(:), allocatable :: cellsum
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      real*8  :: gtot, rel(3), fac(3), q, epsz, vn, cum
      real*8  :: d2a, d2b, d2c, dnom, room_lo, room_hi
      real*8  :: step_side_fwd, step_side_bwd, disp
      integer :: neq, i, j, k, ic, nc90, jhi, jlo, itop(3), ntop, ib
      integer :: n_blk_a, n_blk_b
      logical :: ok1, ok2, oks, okr, bwd1, bwd2
      real*8  :: f_on_those_rows, rounding_of_those_rows
      ! The energy source term at the four probe spacings and at the state,
      ! and the two terms it is the difference of, kept so that the floor
      ! can be read against the CANCELLATION of the row and not against
      ! the row itself.
      real*8, dimension(1-Ng:N+Ng,0:4) :: src
      real*8, dimension(1-Ng:N+Ng) :: heat_at_the_state, cool_at_the_state
      real*8  :: s2a, s2b, s2c, term_largest, cancellation_of_the_row
      character(len=40) :: what
      ! WHAT THE PROBE EVALUATIONS OF THIS HOOK COST, and how the inner
      ! composition solve behaved inside them: the number of full residual
      ! evaluations the hook spends, the wall time of one of them, how many
      ! cell solves the analytic-Jacobian Newton handled and how many fell
      ! through to MINPACK hybrd1. A stopping tolerance too tight for the
      ! Newton shows up as a fallback fraction, and a fallback is where a
      ! hybrd1 exit code other than 1 can appear.
      integer :: n_eval_at_entry, nt_calls_at_entry, nt_fall_at_entry
      integer :: n_eval_here, nt_calls_here, nt_fall_here
      integer(kind=8) :: clock_at_entry, clock_now, clock_rate
      neq = nvar_jac*N
      n_eval_at_entry   = n_resid_eval
      nt_calls_at_entry = nt_calls
      nt_fall_at_entry  = nt_fallback
      call system_clock(clock_at_entry, clock_rate)
      allocate(a1(neq), a2(neq), as(neq), g(neq), vs(neq))
      allocate(gs(neq), sum_of_the_two_actions(neq))
      allocate(Fh1(neq), Fh2(neq), Fh4(neq), cellsum(N))
      gs = 0.0d0;  sum_of_the_two_actions = 0.0d0;  gtot = 0.0d0
      vs = u1 + u2
      fac = (/ 0.1d0, 1.0d0, 10.0d0 /)
      rel = -1.0d0

      ! --- (b) the defect against the length of the arc ---
      do ib = 1, 3
         jv_probe_length_factor = fac(ib)
         call jacobian_action_of_direction(Y, F0, f_sp_base, u1, a1, ok1)
         call jacobian_action_of_direction(Y, F0, f_sp_base, u2, a2, ok2)
         call jacobian_action_of_direction(Y, F0, f_sp_base, vs, as, oks)
         if (.not. (ok1 .and. ok2 .and. oks)) cycle
         g    = as - a1 - a2
         dnom = sqrt(sum((a1 + a2)**2))
         rel(ib) = sqrt(sum(g*g))/max(dnom, 1.0d-300)
         if (ib .eq. 2) then
            ! the standard arc: this is the defect the Arnoldi relation
            ! carries, and the vector localized below.
            gs   = g
            sum_of_the_two_actions = a1 + a2
            gtot = sum(gs*gs)
         endif
      enddo
      jv_probe_length_factor = 1.0d0
      if (rel(2) .lt. 0.0d0) then
         write(*,'(A,A,A)') ' (JFNK) [diag 13] ', trim(pair_name),        &
              ': no admissible sample of one of the three directions'
         deallocate(a1,a2,as,g,gs,sum_of_the_two_actions,vs,              &
                    Fh1,Fh2,Fh4,cellsum)
         return
      endif

      ! --- (a) the defect at the standard arc, and where it sits ---
      write(*,'(A,A,A,ES11.3,A,ES11.3,A,ES11.3)')                         &
           ' (JFNK) [diag 13] ', trim(pair_name),                         &
           ': additivity defect (relative) ', rel(2),                     &
           ', ||J v1 + J v2|| ', sqrt(sum(sum_of_the_two_actions**2)),    &
           ', ||defect|| ', sqrt(gtot)
      ! by the kind of unknown the row belongs to
      do k = 1, nvar_jac
         q = 0.0d0
         do j = 1, N
            q = q + gs(nvar_jac*(j-1)+k)**2
         enddo
         if (q .le. 0.0d0) cycle
         call name_of_row_kind(k, what)
         write(*,'(A,A,A,F8.4)') ' (JFNK) [diag 13]   ', trim(what),      &
              ' rows hold a share of the squared defect of ',             &
              q/max(gtot, 1.0d-300)
      enddo
      ! and by cell: how few cells carry nine tenths of it
      do j = 1, N
         q = 0.0d0
         do k = 1, nvar_jac
            q = q + gs(nvar_jac*(j-1)+k)**2
         enddo
         cellsum(j) = q
      enddo
      cum = 0.0d0;  nc90 = 0;  jlo = N + 1;  jhi = 0;  ntop = 0
      itop = 0
      do while (cum .lt. 0.9d0*gtot .and. nc90 .lt. N)
         ic = 1
         do j = 2, N
            if (cellsum(j) .gt. cellsum(ic)) ic = j
         enddo
         if (cellsum(ic) .le. 0.0d0) exit
         cum  = cum + cellsum(ic)
         nc90 = nc90 + 1
         jlo  = min(jlo, ic);  jhi = max(jhi, ic)
         if (ntop .lt. 3) then
            ntop = ntop + 1;  itop(ntop) = ic
         endif
         cellsum(ic) = 0.0d0
      enddo
      write(*,'(A,I0,A,I0,A,I0,A,I0)') ' (JFNK) [diag 13]   nine tenths'//&
           ' of the squared defect sits in ', nc90, ' of ', N,            &
           ' cells, between cell ', jlo, ' and cell ', jhi
      do i = 1, ntop
         j = itop(i)
         ! the largest row of that cell
         k = 1
         do ic = 2, nvar_jac
            if (abs(gs(nvar_jac*(j-1)+ic)) .gt.                           &
                abs(gs(nvar_jac*(j-1)+k))) k = ic
         enddo
         call name_of_row_kind(k, what)
         write(*,'(A,I0,A,A,A,ES11.3)') ' (JFNK) [diag 13]   cell ', j,   &
              ': largest defect row is ', trim(what), ' at ',             &
              gs(nvar_jac*(j-1)+k)
      enddo

      ! --- (b) printed: the exponent of the arc ---
      if (rel(1) .gt. 0.0d0 .and. rel(3) .gt. 0.0d0) then
         q = log10(rel(3)/rel(1))/2.0d0
      else
         q = 0.0d0
      endif
      write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                             &
           ' (JFNK) [diag 13]   defect at a tenth of the arc ', rel(1),   &
           ', at the arc ', rel(2), ', at ten times it ', rel(3)
      write(*,'(A,F7.3,A)') ' (JFNK) [diag 13]   the defect scales as'//  &
           ' the arc to the power ', q,                                   &
           '  (+1 curvature, -1 rounding or a non-smooth residual,'//     &
           ' 0 a decision of the feasible set)'
      ! WHAT THE PROBE OF THE SUM DIRECTION ACTUALLY DID, which is the
      ! measurement that names the third source outright.
      call jacobian_action_of_direction(Y, F0, f_sp_base, vs, as, oks)
      write(*,'(A,ES11.3,A,ES11.3,A,I0,A,L1)')                            &
           ' (JFNK) [diag 13]   the probe of the sum: nominal step ',     &
           jv_probe_step_nominal_last, ', step taken ',                   &
           jv_probe_step_last, ', components blocked ',                   &
           jv_probe_blocked_last, ', backward side ',                     &
           jv_probe_backward_last
      call jacobian_action_of_direction(Y, F0, f_sp_base, u1, a1, ok1)
      step_side_fwd = jv_probe_step_last
      n_blk_a       = jv_probe_blocked_last
      bwd1          = jv_probe_backward_last
      call jacobian_action_of_direction(Y, F0, f_sp_base, u2, a2, ok2)
      step_side_bwd = jv_probe_step_last
      n_blk_b       = jv_probe_blocked_last
      bwd2          = jv_probe_backward_last
      write(*,'(A,ES11.3,A,I0,A,L1,A,ES11.3,A,I0,A,L1)')                  &
           ' (JFNK) [diag 13]   the probe of v1: step ', step_side_fwd,   &
           ', blocked ', n_blk_a, ', backward ', bwd1,                    &
           ';  of v2: step ', step_side_bwd, ', blocked ', n_blk_b,       &
           ', backward ', bwd2

      ! --- (c) the displacement at the defect-holding unknowns ---
      ! The arc is a length of the WHOLE vector; what a single unknown is
      ! displaced by is eps times that unknown's own component, and a
      ! displacement that is a large fraction of the unknown itself, or of
      ! its room inside the box, is not a derivative of anything.
      vn = sqrt(sum(vs*vs))
      epsz = 0.0d0
      if (vn .gt. 0.0d0) then
         if (jv_probe_on_the_column_scales) then
            epsz = probe_step_on_the_column_scales(Y, vs)
         else
            epsz = probe_length_of_the_jacobian_action(Y)/vn
         endif
      endif
      do i = 1, ntop
         j = itop(i)
         do k = 1, nvar_jac
            ic   = nvar_jac*(j-1) + k
            disp = epsz*abs(vs(ic))
            if (disp .le. 0.0d0) cycle
            room_lo = huge(1.0d0);  room_hi = huge(1.0d0)
            if (allocated(species_box_lo) .and. k .gt. 3) then
               room_lo = Y(ic) - species_box_lo(ic)
               room_hi = species_box_hi(ic) - Y(ic)
            endif
            call name_of_row_kind(k, what)
            write(*,'(A,I0,A,A,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')     &
                 ' (JFNK) [diag 13]   cell ', j, ' ', trim(what),         &
                 ': displacement ', disp, ', unknown ', Y(ic),            &
                 ', column scale ', D(ic),                                &
                 ', displacement over the column scale ',                 &
                 disp/max(abs(D(ic)), 1.0d-300)
            if (room_lo .lt. huge(1.0d0))                                 &
               write(*,'(A,ES11.3,A,ES11.3)')                             &
                    ' (JFNK) [diag 13]     room below ', room_lo,         &
                    ', room above ', room_hi
         enddo
      enddo

      ! --- (d) is the residual smooth on the arc? ---
      ! One-sided second differences along the sum direction at three
      ! spacings, read on the rows that hold the defect. F(Y) - 2 F(Y+h v)
      ! + F(Y+2h v) is h^2 v^T H v + O(h^3) for a smooth residual, so the
      ! three norms below stand in the ratio 1 : 4 : 16. A ratio near one
      ! is a floor: the residual is not smooth at this spacing.
      if (epsz .gt. 0.0d0) then
         ! THE SOURCE TERM OF THE ENERGY ROW AT THE SAME FOUR SPACINGS.
         ! An energy row is heat - cool plus a flux divergence, and heat
         ! and cool are each far larger than their difference in the
         ! photoionized layer, so the row's OWN rounding is epsilon times
         ! the largest term it cancels, not epsilon times the row. The
         ! captures below let the floor be compared against that.
         call eval_residual(Y, f_sp_base, fwork, Fh1,                     &
                            heat, cool, admissible=okr,                   &
                            may_be_adopted=.false.,                       &
                            state_is_discarded=.true.,                    &
                            n_eq_sweeps_fixed=n_eq_sweeps_model)
         src(:,0) = heat - cool
         heat_at_the_state = heat
         cool_at_the_state = cool
         call eval_residual(Y + 0.5d0*epsz*vs, f_sp_base, fwork, Fh1,     &
                            heat, cool, admissible=okr,                   &
                            may_be_adopted=.false.,                       &
                            state_is_discarded=.true.,                    &
                            n_eq_sweeps_fixed=n_eq_sweeps_model)
         src(:,1) = heat - cool
         if (okr) then
            call eval_residual(Y + epsz*vs, f_sp_base, fwork, Fh2,        &
                               heat, cool, admissible=okr,                &
                               may_be_adopted=.false.,                    &
                               state_is_discarded=.true.,                 &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            src(:,2) = heat - cool
         endif
         if (okr) then
            call eval_residual(Y + 2.0d0*epsz*vs, f_sp_base, fwork, as,   &
                               heat, cool, admissible=okr,                &
                               may_be_adopted=.false.,                    &
                               state_is_discarded=.true.,                 &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            src(:,3) = heat - cool
         endif
         if (okr) then
            call eval_residual(Y + 4.0d0*epsz*vs, f_sp_base, fwork, Fh4,  &
                               heat, cool, admissible=okr,                &
                               may_be_adopted=.false.,                    &
                               state_is_discarded=.true.,                 &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            src(:,4) = heat - cool
         endif
         if (okr) then
            d2a = 0.0d0;  d2b = 0.0d0;  d2c = 0.0d0
            do i = 1, ntop
               j = itop(i)
               do k = 1, nvar_jac
                  ic  = nvar_jac*(j-1) + k
                  d2a = d2a + (F0(ic) - 2.0d0*Fh1(ic) + Fh2(ic))**2
                  d2b = d2b + (F0(ic) - 2.0d0*Fh2(ic) + as(ic))**2
                  d2c = d2c + (F0(ic) - 2.0d0*as(ic)  + Fh4(ic))**2
               enddo
            enddo
            d2a = sqrt(d2a);  d2b = sqrt(d2b);  d2c = sqrt(d2c)
            write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                       &
                 ' (JFNK) [diag 13]   second difference of the residual'//&
                 ' on those rows: at half the arc ', d2a,                 &
                 ', at the arc ', d2b, ', at twice it ', d2c
            write(*,'(A,F9.3,A,F9.3,A)')                                  &
                 ' (JFNK) [diag 13]     ratios ', d2b/max(d2a,1.0d-300),  &
                 ' and ', d2c/max(d2b,1.0d-300),                          &
                 ' (4 and 4 for a residual smooth on the arc, 1 and 1'//  &
                 ' for a floor)'
            ! WHAT SIZE OF FLOOR THE ASSEMBLY OF THE RESIDUAL ITSELF CAN
            ! EXPLAIN. Double precision leaves the rows of F reproducible
            ! to about epsilon times their own size; a floor far above
            ! that is not the arithmetic of the assembly but something
            ! inside the evaluation that stops on a tolerance of its own.
            f_on_those_rows = 0.0d0
            do i = 1, ntop
               j = itop(i)
               do k = 1, nvar_jac
                  ic = nvar_jac*(j-1) + k
                  f_on_those_rows = f_on_those_rows + F0(ic)**2
               enddo
            enddo
            f_on_those_rows        = sqrt(f_on_those_rows)
            rounding_of_those_rows = epsilon(1.0d0)*f_on_those_rows
            write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                       &
                 ' (JFNK) [diag 13]     ||F|| on those rows ',            &
                 f_on_those_rows, ', epsilon times it ',                  &
                 rounding_of_those_rows,                                  &
                 ', the floor over that rounding ',                       &
                 d2b/max(rounding_of_those_rows, 1.0d-300)
            ! THE SAME FLOOR READ AGAINST THE CANCELLATION OF THE ROW.
            ! The second difference of the energy source term alone, and
            ! the largest of heat and cool on the same cells: a row that
            ! is the difference of two terms 1/x of its own size carries
            ! a rounding floor x times its own epsilon.
            s2a = 0.0d0;  s2b = 0.0d0;  s2c = 0.0d0
            term_largest = 0.0d0
            do i = 1, ntop
               j = itop(i)
               s2a = s2a + (src(j,0) - 2.0d0*src(j,1) + src(j,2))**2
               s2b = s2b + (src(j,0) - 2.0d0*src(j,2) + src(j,3))**2
               s2c = s2c + (src(j,0) - 2.0d0*src(j,3) + src(j,4))**2
               term_largest = max(term_largest, abs(heat_at_the_state(j)),&
                                  abs(cool_at_the_state(j)))
            enddo
            s2a = sqrt(s2a);  s2b = sqrt(s2b);  s2c = sqrt(s2c)
            cancellation_of_the_row = 0.0d0
            if (f_on_those_rows .gt. 0.0d0)                               &
               cancellation_of_the_row = term_largest/f_on_those_rows
            write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                       &
                 ' (JFNK) [diag 13]     second difference of the energy'//&
                 ' source alone: at half the arc ', s2a,                  &
                 ', at the arc ', s2b, ', at twice it ', s2c
            write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                       &
                 ' (JFNK) [diag 13]     largest of heat and cool on'//    &
                 ' those cells ', term_largest,                           &
                 ', epsilon times it ', epsilon(1.0d0)*term_largest,      &
                 ', the floor over THAT rounding ',                       &
                 d2b/max(epsilon(1.0d0)*term_largest, 1.0d-300)
            write(*,'(A,ES11.3)')                                         &
                 ' (JFNK) [diag 13]     the row cancels its largest'//    &
                 ' term by a factor ', cancellation_of_the_row
            ! WHICH ROW CLASS AND WHICH CELL HOLD THE FLOOR. The source
            ! is one term of the energy row only; a floor that sits in the
            ! mass row is in the flux assembly and cannot be the
            ! chemistry.
            do k = 1, nvar_jac
               q = 0.0d0
               do j = 1, N
                  ic = nvar_jac*(j-1) + k
                  q  = q + (F0(ic) - 2.0d0*Fh2(ic) + as(ic))**2
               enddo
               call name_of_row_kind(k, what)
               write(*,'(A,A,A,ES11.3)')                                  &
                    ' (JFNK) [diag 13]     second difference over the'//  &
                    ' whole column, ', trim(what), ' rows ', sqrt(q)
            enddo
            ! AND WHERE ALONG THE COLUMN IT SITS. A floor confined to the
            ! first cells is a boundary construction; one spread over the
            ! column is the interior scheme.
            do j = 1, N
               q = 0.0d0
               do k = 1, nvar_jac
                  ic = nvar_jac*(j-1) + k
                  q  = q + (F0(ic) - 2.0d0*Fh2(ic) + as(ic))**2
               enddo
               cellsum(j) = q
            enddo
            q   = sum(cellsum)
            cum = 0.0d0
            nc90 = 0;  jlo = N;  jhi = 1
            do i = 1, N
               k = maxloc(cellsum(1:N), 1)
               cum = cum + cellsum(k)
               nc90 = nc90 + 1
               jlo = min(jlo, k);  jhi = max(jhi, k)
               cellsum(k) = -1.0d0
               if (cum .ge. 0.9d0*q) exit
            enddo
            write(*,'(A,I0,A,I0,A,I0,A,I0)')                              &
                 ' (JFNK) [diag 13]     nine tenths of the squared'//     &
                 ' second difference sits in ', nc90, ' of ', N,          &
                 ' cells, between cell ', jlo, ' and cell ', jhi
         else
            write(*,'(A)') ' (JFNK) [diag 13]   the residual is not'//    &
                 ' admissible at one of the four spacings'
         endif
      endif
      call system_clock(clock_now, clock_rate)
      n_eval_here   = n_resid_eval - n_eval_at_entry
      nt_calls_here = nt_calls     - nt_calls_at_entry
      nt_fall_here  = nt_fallback  - nt_fall_at_entry
      if (n_eval_here .gt. 0 .and. clock_rate .gt. 0)                     &
         write(*,'(A,I0,A,ES11.3,A,I0,A,I0,A,F7.3,A)')                    &
              ' (JFNK) [diag 13]   cost: ', n_eval_here,                  &
              ' residual evaluations at ',                                &
              real(clock_now - clock_at_entry,8)                          &
              /real(clock_rate,8)/real(n_eval_here,8),                    &
              ' s each; inner cell solves ', nt_calls,                    &
              ', of them fell back to hybrd1 ', nt_fall_here, ' (',       &
              100.0d0*real(nt_fall_here,8)                                &
              /real(max(nt_calls_here,1),8), '%)'
      deallocate(a1,a2,as,g,gs,sum_of_the_two_actions,vs,              &
                 Fh1,Fh2,Fh4,cellsum)
      end subroutine measure_the_additivity_of_the_jacobian_action

      ! ------------------------------------------------------!

      subroutine pipeline_of_the_residual_at_a_point(Y, f_sp_base, jcell, &
                     q, qname, nq, ok)
      ! EVERY QUANTITY THE ROW OF ONE CELL IS BUILT FROM, in the order the
      ! evaluation builds them:
      !
      !   unknowns -> ghost cells (Apply_BC and the base boundary)
      !            -> temperature and particle count of the sweep
      !            -> reconstructed face states WL, WR of the cell's two
      !               faces
      !            -> interface flux and face pressure there
      !            -> the flux difference dF, the geometric source S and
      !               the radiative source heat - cool of the row.
      !
      ! Read in that order, the FIRST entry that steps between two nearby
      ! states names the operator that is not continuous; everything after
      ! it is downstream of the step and steps with it.
      !
      ! The reconstruction is re-run here on the captured state, which is
      ! the same call assemble_residual makes on the same input, so WL, WR
      ! and the interface fluxes are the ones the row was assembled from.
      ! The positivity-limiter ledger is restored afterwards: this is a
      ! measurement of a state, not a step of the run.
      real*8, dimension(nvar_jac*N),          intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp_base
      integer,                                intent(in)  :: jcell
      real*8, dimension(:),                   intent(out) :: q
      character(len=32), dimension(:),        intent(out) :: qname
      integer,                                intent(out) :: nq
      logical,                                intent(out) :: ok
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(nvar_jac*N)  :: Fv
      real*8, dimension(1-Ng:N+Ng)   :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dFa, Sa
      integer :: g, i, n_limited_before, n_roe_hlle_before
      integer :: n_limited_at_entry, n_roe_hlle_at_entry
      character(len=8) :: cname(3)
      cname = (/ 'rho     ', 'rho v   ', 'E       ' /)
      n_limited_at_entry  = n_faces_positivity_limited
      n_roe_hlle_at_entry = n_faces_roe_hlle
      resid_capture_operator_state = .true.
      call eval_residual(Y, f_sp_base, fwork, Fv, heat, cool,             &
                         admissible=ok, may_be_adopted=.false.,           &
                         state_is_discarded=.true.,                       &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      resid_capture_operator_state = .false.
      if (.not. ok) then
         nq = 0
         return
      endif
      n_limited_before  = n_faces_positivity_limited
      n_roe_hlle_before = n_faces_roe_hlle
      call reconstruction_continuation_rhs(captured_state_with_ghosts,    &
                                           WL, WR, dFa, Sa)
      nq = 0
      ! the ghost cells, lower then upper
      do g = 1, Ng
         do i = 1, 3
            nq = nq + 1
            q(nq) = captured_state_with_ghosts(i,1-g)
            write(qname(nq),'(A,I0,A,A)') 'ghost cell ', 1-g, ' ',        &
                 trim(cname(i))
         enddo
      enddo
      do g = 1, Ng
         do i = 1, 3
            nq = nq + 1
            q(nq) = captured_state_with_ghosts(i,N+g)
            write(qname(nq),'(A,I0,A,A)') 'ghost cell ', N+g, ' ',        &
                 trim(cname(i))
         enddo
      enddo
      nq = nq + 1;  q(nq) = captured_temperature(jcell)
      qname(nq) = 'T of the cell'
      nq = nq + 1;  q(nq) = captured_particle_count(jcell)
      qname(nq) = 'n_tot + n_e of the cell'
      ! the two faces of the cell: face jcell-1 below, face jcell above
      do i = 1, 3
         nq = nq + 1;  q(nq) = WL(i,jcell-1)
         qname(nq) = 'WL at the lower face, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = WR(i,jcell-1)
         qname(nq) = 'WR at the lower face, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = WL(i,jcell)
         qname(nq) = 'WL at the upper face, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = WR(i,jcell)
         qname(nq) = 'WR at the upper face, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = face_flux(i,jcell-1)
         qname(nq) = 'flux at the lower face, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = face_flux(i,jcell)
         qname(nq) = 'flux at the upper face, '//trim(cname(i))
      enddo
      nq = nq + 1;  q(nq) = face_p(jcell-1)
      qname(nq) = 'face pressure below'
      nq = nq + 1;  q(nq) = face_p(jcell)
      qname(nq) = 'face pressure above'
      do i = 1, 3
         nq = nq + 1;  q(nq) = dFa(i,jcell)
         qname(nq) = 'dF of the row, '//trim(cname(i))
      enddo
      do i = 1, 3
         nq = nq + 1;  q(nq) = Sa(i,jcell)
         qname(nq) = 'S of the row, '//trim(cname(i))
      enddo
      nq = nq + 1;  q(nq) = heat(jcell)
      qname(nq) = 'heat of the cell'
      nq = nq + 1;  q(nq) = cool(jcell)
      qname(nq) = 'cool of the cell'
      ! WHETHER A DISCRETE BRANCH OF THE ASSEMBLY FIRED, as a count. A
      ! counter that changes between two nearby states names the branch;
      ! one that stands still says the assembly took the same code path
      ! on both sides of the step and the step is arithmetic.
      nq = nq + 1;  q(nq) = real(n_faces_roe_hlle - n_roe_hlle_before, 8)
      qname(nq) = 'faces that fell back to HLLE'
      nq = nq + 1
      q(nq) = real(n_faces_positivity_limited - n_limited_before, 8)
      qname(nq) = 'face states the limiter scaled'
      ! THE GEOMETRY THAT AMPLIFIES WHATEVER THE FACES CARRY. The row is
      ! a flux difference divided by the cell volume, so an error of the
      ! interface flux enters the row multiplied by r^2/dV, which is
      ! about one over the cell width.
      nq = nq + 1;  q(nq) = dr_j(jcell)
      qname(nq) = 'width of the cell'
      nq = nq + 1
      q(nq) = r_edg(jcell)**2                                            &
              /((r_edg(jcell)**3 - r_edg(jcell-1)**3)/3.0d0)
      qname(nq) = 'r^2 over the cell volume'
      ! AND THE CANCELLATION AT THE FACE. The Roe flux is built from the
      ! DIFFERENCES of the two reconstructed states, and in a nearly
      ! hydrostatic layer those differences are decades below the states
      ! themselves, so the flux carries the rounding of the states and
      ! not of its own value.
      nq = nq + 1;  q(nq) = WR(3,jcell-1) - WL(3,jcell-1)
      qname(nq) = 'pressure jump, lower face'
      nq = nq + 1;  q(nq) = WR(1,jcell-1) - WL(1,jcell-1)
      qname(nq) = 'density jump, lower face'
      nq = nq + 1;  q(nq) = WR(3,jcell) - WL(3,jcell)
      qname(nq) = 'pressure jump, upper face'
      nq = nq + 1;  q(nq) = WR(1,jcell) - WL(1,jcell)
      qname(nq) = 'density jump, upper face'
      ! The ledgers are restored: this is a measurement of a state, not a
      ! step of the run, and the counts the run reports must stay the
      ! counts of its own updates.
      n_faces_positivity_limited = n_limited_at_entry
      n_faces_roe_hlle           = n_roe_hlle_at_entry
      end subroutine pipeline_of_the_residual_at_a_point

      ! ------------------------------------------------------!

      real*8 function median_of_the_magnitudes(x)
      ! Median of |x|, by insertion sort of a copy. The scan needs a
      ! measure of how large a difference of F between two neighboring
      ! samples USUALLY is, and the mean is not it: one step of the map
      ! inside the sampled interval would raise the mean it is compared
      ! against. Half the samples lie below the median whatever the other
      ! half do.
      real*8, dimension(:), intent(in) :: x
      real*8, dimension(size(x)) :: a
      real*8  :: t
      integer :: i, k, m
      m = size(x)
      a = abs(x)
      do i = 2, m
         t = a(i)
         k = i - 1
         do while (k .ge. 1)
            if (a(k) .le. t) exit
            a(k+1) = a(k)
            k = k - 1
         enddo
         a(k+1) = t
      enddo
      if (mod(m,2) .eq. 0) then
         median_of_the_magnitudes = 0.5d0*(a(m/2) + a(m/2+1))
      else
         median_of_the_magnitudes = a((m+1)/2)
      endif
      end function median_of_the_magnitudes

      ! ------------------------------------------------------!

      subroutine steps_in_a_sampled_row(f, factor, nstep, marked)
      ! WHICH INTERVALS OF A UNIFORMLY SAMPLED ROW HOLD A STEP OF THE MAP.
      !
      ! For a map that is continuous and sampled finely enough, the
      ! successive differences f(k+1) - f(k) are its derivative times the
      ! spacing and vary slowly with k: over a window of a few dozen
      ! samples they stay within a small factor of their own median. A
      ! STEP of the map lands entirely in one interval and makes that one
      ! difference stand far above the median, whatever the derivative
      ! does. The test is therefore a comparison against the row's own
      ! median difference and carries no absolute scale of its own, so it
      ! reads a row of any size and any units.
      !
      ! factor is how far above the median counts as a step. A row whose
      ! differences are all zero has no median to compare against and is
      ! reported as having none.
      real*8,  dimension(:), intent(in)  :: f
      real*8,                intent(in)  :: factor
      integer,               intent(out) :: nstep
      logical, dimension(:), intent(out) :: marked
      real*8, dimension(size(f)-1) :: dd
      real*8  :: med
      integer :: k
      nstep  = 0
      marked = .false.
      do k = 1, size(f) - 1
         dd(k) = f(k+1) - f(k)
      enddo
      med = median_of_the_magnitudes(dd)
      if (.not. (med .gt. 0.0d0)) return
      do k = 1, size(f) - 1
         if (abs(dd(k)) .gt. factor*med) then
            marked(k) = .true.
            nstep     = nstep + 1
         endif
      enddo
      end subroutine steps_in_a_sampled_row

      ! ------------------------------------------------------!

      subroutine scan_the_residual_along_a_direction(Y, F0, f_sp_base,    &
                     v, dir_name)
      ! WHERE THE RESIDUAL STEPS ALONG ONE DIRECTION, AT WHAT SCALE, AND IN
      ! WHICH OPERATOR.
      !
      ! A second difference of F that stands on a FLOOR instead of falling
      ! with the square of the spacing (N31, N32) says F is not smooth on
      ! the probe arc: a kink gives a second difference proportional to the
      ! spacing, curvature gives its square, and a floor is what a map that
      ! STEPS somewhere inside the arc gives. This routine finds the steps.
      !
      ! F is sampled at 65 equally spaced points of Y + t v over a window
      ! [-h, +h], and the window is then shrunk by a factor four, five
      ! times over, from two probe arcs down to two arcs divided by 256.
      ! ONE WIDTH IS NOT ENOUGH: a window that holds many steps has no
      ! interval standing above its neighbors, because every interval
      ! holds a step, and a window narrower than one step is smooth. Three
      ! measures separate the two, taken on each row of each window:
      !
      !   the TOTAL VARIATION of the row over the window against the NET
      !   change across it. A map with a derivative gives them equal; a
      !   dither of steps gives a total variation that stays where it is
      !   as the window shrinks while the net change falls with the width.
      !
      !   the MEDIAN successive difference. For a smooth row it falls in
      !   proportion to the spacing; for a dither of steps it stops
      !   falling at the size of one step.
      !
      !   an interval whose difference stands a hundred times above the
      !   row's own median difference: an ISOLATED step, which appears
      !   once the window is narrow enough to hold few of them.
      !
      ! At the largest isolated step the row is split by TERM and by
      ! OPERATOR (pipeline_of_the_residual_at_a_point) over three
      ! consecutive intervals, the middle one holding the step: the
      ! earliest quantity of the pipeline whose change over the stepping
      ! interval stands far above its change over the two neighboring
      ! intervals is the one that steps, and it names the operator.
      !
      ! TEST ONLY. Every evaluation is discarded and nothing outside the
      ! prints reads any of it.
      real*8, dimension(nvar_jac*N),          intent(in) :: Y, F0
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(nvar_jac*N),          intent(in) :: v
      character(len=*),                       intent(in) :: dir_name
      integer, parameter :: n_scan = 65
      integer, parameter :: n_win  = 10
      ! How far above the row's own median difference a difference has to
      ! stand to be called an isolated step.
      real*8,  parameter :: step_factor = 1.0d2
      real*8, dimension(:,:), allocatable :: Fs
      real*8, dimension(:),   allocatable :: heat, cool, Fv
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      integer, dimension(:), allocatable :: n_steps_at_k
      real*8, dimension(n_scan-1) :: dd
      logical, dimension(n_scan-1) :: marked
      real*8  :: eps0, dt, vn, big, relbig, tk, hw, med, tv, net
      real*8  :: med_track, tv_track, net_track, relmed, d2_track
      real*8  :: d2_first, geo, wmax, pjmp, pface
      real*8, dimension(N) :: cellwork
      real*8  :: q
      integer :: j_track
      real*8  :: q_at(4,200), dq(3), dmid, dnb
      character(len=32) :: qname(200)
      character(len=40) :: what
      integer :: neq, i, j, k, kk, ib, nqk, n_row_steps, n_k_with_step
      integer :: kbig, ibig, jbig, kindbig, jlo, jhi, nq_at, nzero
      integer :: n_by_kind(nvar_jac), n_bad, nstep, iw
      integer :: i_track, n_sign_changes
      integer :: kstep_at, iw_step, nq_here
      real*8  :: hw_step, dt_step
      logical :: okr, ok4(4)
      neq = nvar_jac*N
      vn  = sqrt(sum(v*v))
      if (vn .le. 0.0d0) return
      eps0 = probe_length_of_the_jacobian_action(Y)/vn
      if (.not. (eps0 .gt. 0.0d0)) return
      allocate(Fs(neq,n_scan), Fv(neq))
      allocate(heat(1-Ng:N+Ng), cool(1-Ng:N+Ng), n_steps_at_k(n_scan-1))
      write(*,'(A,A,A,ES11.3,A,I0,A,I0,A)')                               &
           ' (JFNK) [diag 14] ', trim(dir_name),                          &
           ': the probe arc along it is ', eps0, '; ', n_scan,            &
           ' samples of F on each of ', n_win,                            &
           ' windows, each a quarter of the one before'
      i_track  = 0
      kstep_at = 0;  iw_step = 0;  hw_step = 0.0d0;  dt_step = 0.0d0
      ibig = 0;  kbig = 0
      do iw = 1, n_win
         hw = 2.0d0*eps0/4.0d0**(iw-1)
         dt = 2.0d0*hw/real(n_scan-1, 8)
         n_bad = 0
         do k = 1, n_scan
            tk = -hw + real(k-1,8)*dt
            call eval_residual(Y + tk*v, f_sp_base, fwork, Fv, heat,      &
                               cool, admissible=okr,                      &
                               may_be_adopted=.false.,                    &
                               state_is_discarded=.true.,                 &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            Fs(:,k) = Fv
            if (.not. okr) n_bad = n_bad + 1
         enddo
         ! THE ROW THIS WINDOW IS READ ON. Chosen once, on the widest
         ! window, as the row whose successive differences are the largest
         ! fraction of the row itself: the row along which this direction
         ! moves the residual most, which is where a floor of the map
         ! reaches the difference quotient first.
         if (iw .eq. 1) then
            ! THE ROW THAT CARRIES THE FLOOR, which is the row this whole
            ! measurement is about: the largest second difference of F
            ! along the direction. A smooth row's second difference is
            ! the curvature times the square of the spacing and is
            ! negligible here; a row that carries a step of the map has a
            ! second difference of the size of the step.
            relmed = 0.0d0
            do i = 1, neq
               med = 0.0d0
               do k = 1, n_scan-2
                  med = max(med, abs(Fs(i,k) - 2.0d0*Fs(i,k+1)           &
                                     + Fs(i,k+2)))
               enddo
               if (med .gt. relmed) then
                  relmed  = med
                  i_track = i
               endif
            enddo
            if (i_track .le. 0) i_track = 1
            j       = (i_track-1)/nvar_jac + 1
            kindbig = mod(i_track-1, nvar_jac) + 1
            call name_of_row_kind(kindbig, what)
            write(*,'(A,I0,A,A,A,ES11.3,A,ES11.3)')                       &
                 ' (JFNK) [diag 14]   read on cell ', j, ', the ',        &
                 trim(what), ' row, whose value at the state is ',        &
                 F0(i_track), ' and whose largest second difference'//    &
                 ' over the widest window is ', relmed
         endif
         ! the tracked row's three measures on this window
         do k = 1, n_scan-1
            dd(k) = Fs(i_track,k+1) - Fs(i_track,k)
         enddo
         med_track = median_of_the_magnitudes(dd)
         d2_track  = 0.0d0
         do k = 1, n_scan-2
            d2_track = max(d2_track, abs(Fs(i_track,k)                    &
                           - 2.0d0*Fs(i_track,k+1) + Fs(i_track,k+2)))
         enddo
         if (iw .eq. 1) d2_first = d2_track
         tv_track  = sum(abs(dd))
         net_track = abs(Fs(i_track,n_scan) - Fs(i_track,1))
         nzero     = 0
         n_sign_changes = 0
         do k = 1, n_scan-1
            if (dd(k) .eq. 0.0d0) nzero = nzero + 1
            if (k .gt. 1) then
               if (dd(k)*dd(k-1) .lt. 0.0d0)                              &
                  n_sign_changes = n_sign_changes + 1
            endif
         enddo
         ! isolated steps, over every row of the column
         n_steps_at_k = 0;  n_by_kind = 0;  n_row_steps = 0
         big = 0.0d0;  relbig = 0.0d0
         jlo = N + 1;  jhi = 0
         do i = 1, neq
            call steps_in_a_sampled_row(Fs(i,:), step_factor, nstep,      &
                                        marked)
            if (nstep .eq. 0) cycle
            j   = (i-1)/nvar_jac + 1
            jlo = min(jlo, j);  jhi = max(jhi, j)
            kk  = mod(i-1, nvar_jac) + 1
            n_by_kind(kk) = n_by_kind(kk) + nstep
            n_row_steps   = n_row_steps + nstep
            do k = 1, n_scan-1
               if (.not. marked(k)) cycle
               n_steps_at_k(k) = n_steps_at_k(k) + 1
               if (abs(Fs(i,k+1) - Fs(i,k)) .gt. big) then
                  big    = abs(Fs(i,k+1) - Fs(i,k))
                  relbig = big/max(abs(F0(i)), 1.0d-300)
                  kbig   = k
                  ibig   = i
               endif
            enddo
         enddo
         n_k_with_step = 0
         do k = 1, n_scan-1
            if (n_steps_at_k(k) .gt. 0) n_k_with_step = n_k_with_step + 1
         enddo
         write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,F7.4)')          &
              ' (JFNK) [diag 14]   window +-', hw, ' (', hw/eps0,         &
              ' arcs), spacing ', dt, ': median difference ', med_track,  &
              ', total variation over the net change ',                   &
              tv_track/max(net_track, 1.0d-300)
         write(*,'(A,ES10.3,A,ES10.3,A,I0,A,I0,A,I0)')                    &
              ' (JFNK) [diag 14]     total variation ', tv_track,         &
              ', net change ', net_track, ', differences that are'//      &
              ' exactly zero ', nzero, ' of ', n_scan-1,                  &
              ', sign changes ', n_sign_changes
         write(*,'(A,ES10.3,A,ES10.3,A,I0,A,I0,A,I0)')                    &
              ' (JFNK) [diag 14]     largest second difference ',         &
              d2_track, ', over the median difference ',                  &
              d2_track/max(med_track, 1.0d-300),                          &
              ', differences that are exactly zero ', nzero,              &
              ' of ', n_scan-1, ', sign changes ', n_sign_changes
         if (n_bad .gt. 0)                                                &
            write(*,'(A,I0,A)') ' (JFNK) [diag 14]     ', n_bad,          &
                 ' of the samples are not admissible states'
         if (n_row_steps .gt. 0) then
            write(*,'(A,I0,A,I0,A,I0,A,I0)')                              &
                 ' (JFNK) [diag 14]     isolated steps: ', n_row_steps,   &
                 ' (row, interval) pairs in ', n_k_with_step, ' of ',     &
                 n_scan-1, ' intervals, in cells ', jlo
            write(*,'(A,I0,A,ES11.3,A,ES11.3)')                           &
                 ' (JFNK) [diag 14]     to cell ', jhi,                   &
                 '; the largest moves its row by ', big, ', which is ',   &
                 relbig
            do kk = 1, nvar_jac
               if (n_by_kind(kk) .eq. 0) cycle
               call name_of_row_kind(kk, what)
               write(*,'(A,A,A,I0)') ' (JFNK) [diag 14]       ',          &
                    trim(what), ' rows step ', n_by_kind(kk), ' times'
            enddo
            ! keep the FINEST window that isolates a step: that is where
            ! the pipeline on either side of it differs by the step alone.
            kstep_at = kbig;  iw_step = iw
            hw_step  = hw;    dt_step = dt
            jbig     = (ibig-1)/nvar_jac + 1
         else
            write(*,'(A)') ' (JFNK) [diag 14]     no isolated step on'//  &
                 ' this window'
         endif
      enddo
      ! --- WHAT SIZE OF FLOOR THE ARITHMETIC OF THE ASSEMBLY EXPLAINS ---
      ! The row is a flux difference divided by the cell volume, and the
      ! interface flux is built from the DIFFERENCES of two reconstructed
      ! states that a nearly hydrostatic layer makes decades smaller than
      ! the states themselves. The flux therefore carries the last bit of
      ! those O(1) states, not of its own value, and the row carries it
      ! multiplied by r^2 over the cell volume, which is about one over
      ! the cell width. That product is the floor a correctly rounded
      ! assembly cannot go below, and it is the number to compare the
      ! measured floor against; epsilon times the row itself is not.
      j_track = (i_track-1)/nvar_jac + 1
      call pipeline_of_the_residual_at_a_point(Y, f_sp_base, j_track,     &
               q_at(1,:), qname, nqk, okr)
      if (okr) then
         geo  = 0.0d0
         wmax = 0.0d0
         pjmp = 0.0d0;  pface = 0.0d0
         do i = 1, nqk
            if (qname(i) .eq. 'r^2 over the cell volume')                 &
               geo = q_at(1,i)
            if (qname(i)(1:2) .eq. 'WL' .or. qname(i)(1:2) .eq. 'WR')     &
               wmax = max(wmax, abs(q_at(1,i)))
            if (qname(i)(1:13) .eq. 'pressure jump')                      &
               pjmp = max(pjmp, abs(q_at(1,i)))
            if (qname(i) .eq. 'WL at the upper face, E')                  &
               pface = abs(q_at(1,i))
         enddo
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                          &
              ' (JFNK) [diag 14]   the row of that cell stands on'//      &
              ' r^2 over the cell volume ', geo,                          &
              ', its largest face state is ', wmax,                       &
              ', and the largest pressure jump across a face of it is ',  &
              pjmp
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                          &
              ' (JFNK) [diag 14]     the face pressure cancels by a'//    &
              ' factor ', pface/max(pjmp, 1.0d-300),                      &
              '; epsilon times the face state, divided by the cell'//     &
              ' width, is ', epsilon(1.0d0)*wmax*geo,                     &
              ', against the measured floor ', d2_first
         write(*,'(A,F9.3)')                                              &
              ' (JFNK) [diag 14]     the measured floor over that'//      &
              ' arithmetic bound ',                                       &
              d2_first/max(epsilon(1.0d0)*wmax*geo, 1.0d-300)
      ! AND THE SAME COMPARISON CELL BY CELL. If the floor is the
      ! rounding of the assembly then it follows the cell's own geometry
      ! and state and not the size of its row, so the ratio below stays
      ! near one order of magnitude over the whole column, across the
      ! decade the cell width spans. A branch, which fires in one place
      ! and not another, would not.
      if (allocated(captured_state_with_ghosts)) then
         write(*,'(A)') ' (JFNK) [diag 14]   cell by cell, the'//         &
              ' largest second difference of the cell against'//          &
              ' epsilon times its own state over its width:'
         do i = 1, N
            cellwork(i) = 0.0d0
            do k = 1, n_scan-2
               do kk = 1, nvar_jac
                  cellwork(i) = max(cellwork(i),                          &
                     abs(Fs(nvar_jac*(i-1)+kk,k)                          &
                         - 2.0d0*Fs(nvar_jac*(i-1)+kk,k+1)                &
                         + Fs(nvar_jac*(i-1)+kk,k+2)))
               enddo
            enddo
         enddo
         do ib = 1, 10
            j = 1
            do i = 2, N
               if (cellwork(i) .gt. cellwork(j)) j = i
            enddo
            if (.not. (cellwork(j) .gt. 0.0d0)) exit
            geo  = r_edg(j)**2                                            &
                   /((r_edg(j)**3 - r_edg(j-1)**3)/3.0d0)
            wmax = max(abs(captured_state_with_ghosts(1,j)),              &
                       abs(captured_state_with_ghosts(3,j)))
            write(*,'(A,I0,A,ES11.3,A,ES11.3,A,ES11.3,A,F9.3)')           &
                 ' (JFNK) [diag 14]     cell ', j,                        &
                 ': second difference ', cellwork(j), ', width ',         &
                 dr_j(j), ', arithmetic bound ',                          &
                 epsilon(1.0d0)*wmax*geo, ', the ratio ',                 &
                 cellwork(j)/max(epsilon(1.0d0)*wmax*geo, 1.0d-300)
            cellwork(j) = -1.0d0
         enddo
         ! AND ALONG THE COLUMN, where the cell width spans a decade or
         ! more: the bound is a floor and not an equality, so a cell
         ! whose residual has real curvature stands above it, but no
         ! cell may stand below it.
         do ib = 1, 7
            j = min(N, 1 + (ib-1)*N/6)
            geo  = r_edg(j)**2                                            &
                   /((r_edg(j)**3 - r_edg(j-1)**3)/3.0d0)
            wmax = max(abs(captured_state_with_ghosts(1,j)),              &
                       abs(captured_state_with_ghosts(3,j)))
            q = 0.0d0
            do k = 1, n_scan-2
               do kk = 1, nvar_jac
                  q = max(q, abs(Fs(nvar_jac*(j-1)+kk,k)                  &
                                 - 2.0d0*Fs(nvar_jac*(j-1)+kk,k+1)        &
                                 + Fs(nvar_jac*(j-1)+kk,k+2)))
               enddo
            enddo
            write(*,'(A,I0,A,ES11.3,A,ES11.3,A,ES11.3,A,ES10.3)')         &
                 ' (JFNK) [diag 14]     cell ', j,                        &
                 ' of the column: second difference ', q, ', width ',     &
                 dr_j(j), ', arithmetic bound ',                          &
                 epsilon(1.0d0)*wmax*geo, ', the ratio ',                 &
                 q/max(epsilon(1.0d0)*wmax*geo, 1.0d-300)
         enddo
      endif
      endif
      if (iw_step .eq. 0) then
         write(*,'(A)') ' (JFNK) [diag 14]   no window isolates a'//      &
              ' step; no operator split'
         deallocate(Fs, Fv, heat, cool, n_steps_at_k)
         return
      endif
      ! --- the term and the operator that step ---
      kindbig = mod(ibig-1, nvar_jac) + 1
      call name_of_row_kind(kindbig, what)
      tk = -hw_step + real(kstep_at-1,8)*dt_step
      write(*,'(A,I0,A,A,A,ES11.3,A,ES11.3,A)')                           &
           ' (JFNK) [diag 14]   splitting the step of cell ', jbig,       &
           ', the ', trim(what), ' row, at t ', tk, ' (',                 &
           tk/eps0, ' arcs from the state)'
      do ib = 1, 4
         k  = min(max(kstep_at - 1 + ib - 1, 1), n_scan)
         tk = -hw_step + real(k-1,8)*dt_step
         call pipeline_of_the_residual_at_a_point(Y + tk*v, f_sp_base,    &
                  jbig, q_at(ib,:), qname, nqk, ok4(ib))
         if (ib .eq. 1) nq_at = nqk
      enddo
      if (.not. all(ok4)) then
         write(*,'(A)') ' (JFNK) [diag 14]   one of the four states'//    &
              ' around the step is not admissible; no operator split'
         deallocate(Fs, Fv, heat, cool, n_steps_at_k)
         return
      endif
      write(*,'(A)') ' (JFNK) [diag 14]   the pipeline across the'//      &
           ' step, in the order the row is built (value, then the'//      &
           ' change over the interval before the step, over the'//        &
           ' stepping interval, and over the one after):'
      nq_here = nq_at
      do i = 1, nq_here
         dq(1) = q_at(2,i) - q_at(1,i)
         dq(2) = q_at(3,i) - q_at(2,i)
         dq(3) = q_at(4,i) - q_at(3,i)
         dmid  = abs(dq(2))
         dnb   = max(abs(dq(1)), abs(dq(3)))
         write(*,'(A,A,A,ES23.15,A,ES11.3,A,ES11.3,A,ES11.3,A,ES9.2)')    &
              ' (JFNK) [diag 14]     ', qname(i), ' = ', q_at(2,i),       &
              ', changes by ', dq(1), ' | ', dq(2), ' | ', dq(3),         &
              ', the stepping interval over its neighbors ',             &
              dmid/max(dnb, 1.0d-300)
      enddo
      deallocate(Fs, Fv, heat, cool, n_steps_at_k)
      end subroutine scan_the_residual_along_a_direction

      ! ------------------------------------------------------!

      subroutine pgmres(Y, F0, f_sp_base, D, Drow, abf, ipiv, idtau, b, x, m, &
                        rtol, gm_iters, gm_outcome, gm_resid_rel, gm_snorm, &
                        Ax, delta_ball)
      ! Right-preconditioned GMRES(m), single cycle, for the JFNK step, in
      ! the DIAGONALLY SCALED space:
      !   solve  A_z x = b,   A_z v = v*idtau + D^-1 * J*(D v),
      ! (J*(Dv) matrix-free; b and x are scaled quantities; the caller maps
      ! dY = D*x back). Preconditioner M = factored banded SCALED system
      ! (I/dtau + D^-1 J_banded D) supplied in abf.
      real*8, dimension(nvar_jac*N),  intent(in)  :: Y, F0, b, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in)  :: ipiv
      real*8,  intent(in)  :: idtau, rtol
      integer, intent(in)  :: m
      real*8, dimension(nvar_jac*N), intent(out) :: x
      ! The SCALED operator applied to the step it returns, A_z x, taken from
      ! the Arnoldi relation A_z M^-1 V_k = V_{k+1} Hbar_k rather than from a
      ! fresh product, so it costs no residual evaluation. The trust region
      ! needs it for the predicted reduction, and needs it to be the SAME
      ! model the Krylov step was built from. Hbar is the Hessenberg BEFORE
      ! the Givens rotations, which is why it is kept separately below.
      ! OPTIONAL: the three-unknown route never asks for it.
      real*8, dimension(nvar_jac*N), intent(out), optional :: Ax
      ! HOW MANY PRODUCTS OF THE OPERATOR THE CYCLE MANAGED.
      integer, intent(out) :: gm_iters
      ! WHICH OUTCOME THIS CYCLE HAD, one of the gm_* names above. It is not
      ! derivable from gm_iters: a cycle that ran m products may have
      ! converged at exactly m or been cut off there, and a cycle that ran
      ! none may have been refused by the feasible set or by its
      ! preconditioner.
      integer, intent(out) :: gm_outcome
      ! WHAT THE CYCLE REACHED AGAINST WHAT IT WAS ASKED FOR: the relative
      ! residual ||b - A_z x||/||b|| of the step returned. It is the
      ! residual of the reduced least-squares problem wherever that problem
      ! is a faithful statement about the system, and the TRUE residual,
      ! from one extra product of the operator, wherever it is not -- at an
      ! Arnoldi breakdown and on a reduced problem that has lost rank.
      real*8,  intent(out) :: gm_resid_rel
      ! The norm of the step returned, in the scaled coordinates the cycle
      ! works in.
      real*8,  intent(out) :: gm_snorm
      ! THE RADIUS OF THE TRUST BALL THE STEP HAS TO LIE IN, where the caller
      ! has one and the truncation is armed
      ! (krylov_truncated_on_the_trust_ball). The cycle then stops at the
      ! first iterate outside the ball and returns the point of the boundary
      ! on the segment to it, with the image of THAT point, so the caller's
      ! model is the model of the step it is handed. Without it, or with the
      ! truncation off, the cycle is the untruncated one and nothing here
      ! reads this argument.
      real*8, intent(in), optional :: delta_ball

      integer :: neq, i, j, kk, lpinfo, kv_set, krank, jb, ib
      real*8, allocatable :: V(:,:), Hs(:,:), Hb(:,:), gg(:), cs(:), sn(:), yy(:)
      real*8, allocatable :: z(:), w(:), u(:)
      ! The armed-truncation work space: the coefficients of the iterate of
      ! the subspace built so far, that iterate, its image, and the last
      ! pair that was still inside the ball.
      real*8, allocatable :: yb(:), xb(:), Axb(:), xb_in(:), Axb_in(:)
      ! The additivity hook's two directions that are not the basis's.
      real*8, allocatable :: rv1(:), rv2(:)
      real*8 :: beta, denom, tmp, res_now, hcol, hmax, tball, hmin_diag
      real*8 :: gap_col, gap_worst, repro
      logical :: jv_ok, arnoldi_breakdown, verify_true_residual
      logical :: on_the_ball, stopped_on_the_ball
      ! The residual history's own work space and its running quantities
      ! (gm_residual_history_on): the coefficients of the iterate of the
      ! subspace built so far, that iterate, and the true relative residual
      ! it reaches against the operator.
      real*8, allocatable :: yh(:), xh(:)
      real*8  :: res_true_now
      logical :: history_here
      ! WHETHER THIS CYCLE WORKS IN THE EQUILIBRATED ROWS OF THE LINEAR
      ! MODEL (model_row_equilibration_on), and the norm of the right-hand
      ! side on the CERTIFICATION scales, which is what the relative
      ! residual handed back is measured against whether the arm is on or
      ! off.
      logical :: row_equil_here
      real*8  :: beta_cert
      ! The true-residual arm's work space (gm_step_by_its_true_residual):
      ! the best step its residual against the operator has seen, that
      ! step's image, its residual, the residual of the check before, how
      ! many consecutive checks have risen and at which product the best
      ! was found.
      real*8, allocatable :: xbest(:), Axbest(:), yh2(:)
      real*8  :: res_true_best, res_true_last
      integer :: n_true_rise, j_true_best
      logical :: true_resid_here, have_a_best_step, check_here
      ! Whether the reduced problem reaching the tolerance ends the cycle.
      ! It always does without the true-residual arm; with it, only where
      ! the iterate could not be measured against the operator.
      logical :: reduced_tolerance_stops

      neq = nvar_jac*N
      allocate(V(neq,m+1), Hs(m+1,m), gg(m+1), cs(m), sn(m), yy(m))
      allocate(z(neq), w(neq), u(neq))
      ! Allocated last and on its own, so that every array this cycle has
      ! always had keeps the address it has always had (section 139.9).
      allocate(Hb(m+1,m))
      Hs = 0.0d0; gg = 0.0d0; cs = 0.0d0; sn = 0.0d0; x = 0.0d0
      yy = 0.0d0
      ! V is deliberately NOT zeroed: touching an array the cycle has always
      ! had is what moved two atomic cases by 1e-9 in section 139.9. The Ax
      ! sum below skips any column the Arnoldi never wrote, which is only
      ! ever a column whose coefficient is exactly zero anyway.
      Hb = 0.0d0;  kv_set = 1
      if (present(Ax)) Ax = 0.0d0
      ! THE TRUNCATION IS ARMED ONLY WHERE ITS THREE CONDITIONS HOLD: the
      ! switch, a radius from the caller, and the caller asking for the
      ! image -- without the image there is no model of the boundary point
      ! and the truncation would hand back a step no model describes.
      on_the_ball = krylov_truncated_on_the_trust_ball .and.              &
                    present(delta_ball) .and. present(Ax)
      if (on_the_ball) on_the_ball = (delta_ball .gt. 0.0d0)
      stopped_on_the_ball = .false.
      tball = 0.0d0
      if (on_the_ball) then
         allocate(yb(m), xb(neq), Axb(neq), xb_in(neq), Axb_in(neq))
         yb = 0.0d0
         ! The iterate of the empty subspace is the zero step, which is
         ! inside every ball of positive radius, and its image is zero.
         xb_in = 0.0d0;  Axb_in = 0.0d0
      endif
      history_here = gm_residual_history_on .and. gm_residual_history_here
      if (history_here) then
         allocate(yh(m), xh(neq))
         yh = 0.0d0
         write(*,'(A,I0,A,ES11.3)') ' (JFNK) [diag 16] residual history'//&
              ' of a cycle of ', m, ' products, tolerance asked ', rtol
      endif
      ! THE ROW SCALING OF THE LINEAR MODEL, and the arm that lets the
      ! cycle keep the step its true residual chooses. Both are read once
      ! here, so that a cycle works in one model from its first product to
      ! its last.
      row_equil_here  = the_linear_model_rows_are_equilibrated(neq)
      true_resid_here = gm_step_by_its_true_residual
      res_true_best   = huge(1.0d0)
      res_true_last   = huge(1.0d0)
      n_true_rise     = 0
      j_true_best     = 0
      have_a_best_step = .false.
      if (true_resid_here) then
         ! Allocated after everything the cycle has always had, for the
         ! reason V is not zeroed.
         allocate(xbest(neq), Axbest(neq), yh2(m))
         xbest = 0.0d0;  Axbest = 0.0d0;  yh2 = 0.0d0
      endif
      gm_orthogonality_loss_cycle = 0.0d0
      gm_outcome   = gm_tolerance_reached
      gm_resid_rel = 0.0d0
      gm_snorm     = 0.0d0
      arnoldi_breakdown    = .false.
      verify_true_residual = .false.

      ! r0 = b - A*0 = b, in the rows the cycle works in. beta_cert is the
      ! same quantity on the certification scales, and the two are one
      ! number wherever the row equilibration is off.
      if (row_equil_here) then
         beta = sqrt(sum((b*model_row_equilibration)**2))
      else
         beta = sqrt(sum(b*b))
      endif
      beta_cert = sqrt(sum(b*b))
      gm_iters = 0
      if (.not. every_component_is_finite(b)) then
         ! There is no linear system here to solve.
         gm_outcome   = gm_nonfinite
         gm_resid_rel = 1.0d0
         call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,     &
                                           yb,xb,Axb,xb_in,Axb_in)
         return
      endif
      if (beta .le. 0.0d0) then
         ! A ZERO RIGHT-HAND SIDE IS SOLVED BY THE ZERO STEP, and zero is
         ! not a scale to normalize a direction by. Returned before the
         ! preconditioner is applied to anything and before any direction is
         ! formed, so the outcome costs no product.
         call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,     &
                                           yb,xb,Axb,xb_in,Axb_in)
         return
      endif
      ! The residual of the empty subspace is the right-hand side itself.
      res_now = beta
      if (row_equil_here) then
         V(:,1) = b*model_row_equilibration/beta
      else
         V(:,1) = b/beta
      endif
      gg(1)  = beta

      do j = 1, m
         gm_iters = j
         ! z = M^{-1} V(:,j)   (right preconditioning)
         z = V(:,j)
         call banded_preconditioner_solve(abf, ipiv, neq, z, lpinfo)
         if (lpinfo .ne. 0 .or. .not. every_component_is_finite(z)) then
            ! THE PRECONDITIONER IS NOT A MAP HERE. A nonzero LAPACK status
            ! or a nonfinite result means the triangular solution did not
            ! happen, and the vector in z is not M^-1 V(:,j); sampling the
            ! operator along it would build the Arnoldi recursion on a
            ! direction of the preconditioner's arithmetic and not of the
            ! system.
            gm_iters   = j - 1
            gm_outcome = gm_preconditioner_failed
            exit
         endif
         ! THE DIRECTION IS PROJECTED ONTO THE ACTIVE SHARED ROWS BEFORE
         ! THE OPERATOR IS APPLIED TO IT, so that the operator this
         ! recursion builds a basis of is the one the step is taken with
         ! (project_out_of_the_active_element_constraints). The projection
         ! is a fixed linear idempotent map for the whole step, so the
         ! Arnoldi relation is a relation of the composition and the
         ! predicted reduction may be read from it.
         call project_out_of_the_active_element_constraints(z)
         ! w = A_z z = z*idtau + D^-1 * J*(D z)
         call jacobian_action_of_direction(Y, F0, f_sp_base, D*z, w, jv_ok)
         if (.not. jv_ok) then
            ! No admissible sample of the residual along this direction:
            ! stop the cycle here rather than feed the Arnoldi recursion a
            ! product that does not exist. With not one direction sampled
            ! the step is zero, and which of the two happened is the
            ! difference between a feasible set with no room along the first
            ! preconditioned direction and a cycle that ran out of them.
            gm_iters   = j - 1
            if (j .eq. 1) then
               gm_outcome = gm_no_direction_sampled
            else
               gm_outcome = gm_jacobian_action_unusable
            endif
            exit
         endif
         ! THE THREE-UNKNOWN BRANCH IS THE ORIGINAL EXPRESSION, CHARACTER
         ! FOR CHARACTER. Drow holds exactly D's values when there is no
         ! carrier row, so the two branches are the same arithmetic -- but
         ! not the same instruction sequence, and a reduction summed in a
         ! different order differs in its last bits. A change with no
         ! intended effect has to be byte-identical (section 139.9).
         if (nspec_row .gt. 0) then
            w = w/Drow + idtau*(D/Drow)*z
         else
            w = w/D + idtau*z
         endif
         ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
         ! the operator the cycle works with is E A_z and the band it is
         ! preconditioned by is E times the same band, so the two agree
         ! (unit_infinity_norm_row_scaling_of_the_band).
         if (row_equil_here) w = w*model_row_equilibration
         if (.not. every_component_is_finite(w)) then
            ! The action returned a value that is not a number. The subspace
            ! built so far is still a subspace; this direction is not.
            gm_iters   = j - 1
            gm_outcome = gm_jacobian_action_unusable
            exit
         endif
         ! Modified Gram-Schmidt
         hcol = 0.0d0
         do i = 1, j
            Hs(i,j) = sum(w*V(:,i))
            w = w - Hs(i,j)*V(:,i)
            hcol = max(hcol, abs(Hs(i,j)))
         enddo
         if (gm_reorthogonalize) then
            ! A SECOND PASS OVER THE SAME BASIS (gm_reorthogonalize). One
            ! pass leaves components along the earlier basis vectors of the
            ! size of the rounding times the conditioning of the operator;
            ! the second pass removes them and its coefficients belong in
            ! the Hessenberg column, because the column IS the coordinates
            ! of the image in the basis and the Arnoldi relation the trust
            ! region's model is read from holds only if it carries all of
            ! them. No product of the operator is taken here.
            do i = 1, j
               tmp     = sum(w*V(:,i))
               w       = w - tmp*V(:,i)
               Hs(i,j) = Hs(i,j) + tmp
               hcol    = max(hcol, abs(Hs(i,j)))
            enddo
         endif
         Hs(j+1,j) = sqrt(sum(w*w))
         if (.not. finite_real(Hs(j+1,j))) then
            gm_iters   = j - 1
            gm_outcome = gm_nonfinite
            exit
         endif
         ! THE REMAINDER CARRIES NO NEW DIRECTION. Measured against the
         ! column it belongs to, so that the test is about information and
         ! not about the units the operator happens to be scaled in. Where
         ! it holds, the Krylov space is invariant: no basis vector j+1 is
         ! written and the cycle may not build a column on one.
         arnoldi_breakdown = (Hs(j+1,j) .le.                              &
                              gm_rank_relative*max(hcol, beta))
         if (.not. arnoldi_breakdown) then
            V(:,j+1) = w/Hs(j+1,j)
            kv_set   = j + 1
            if (gm_measure_orthogonality) then
               do i = 1, j
                  gm_orthogonality_loss_max =                             &
                     max(gm_orthogonality_loss_max,                       &
                         abs(sum(V(:,i)*V(:,j+1))))
                  gm_orthogonality_loss_cycle =                           &
                     max(gm_orthogonality_loss_cycle,                     &
                         abs(sum(V(:,i)*V(:,j+1))))
               enddo
            endif
         endif
         ! The Arnoldi column as built, kept before the rotations overwrite it.
         Hb(1:j+1,j) = Hs(1:j+1,j)
         ! Apply previous Givens rotations to column j
         do i = 1, j-1
            tmp       =  cs(i)*Hs(i,j) + sn(i)*Hs(i+1,j)
            Hs(i+1,j) = -sn(i)*Hs(i,j) + cs(i)*Hs(i+1,j)
            Hs(i,j)   =  tmp
         enddo
         ! New Givens rotation to zero Hs(j+1,j). hypot is the length of
         ! the two-component vector without forming either square, so the
         ! rotation exists over the whole range in which its entries do.
         denom  = hypot(Hs(j,j), Hs(j+1,j))
         if (denom .gt. 0.0d0) then
            cs(j)  = Hs(j,j)/denom
            sn(j)  = Hs(j+1,j)/denom
         else
            ! BOTH ENTRIES OF THE COLUMN VANISH. There is no rotation that
            ! zeroes anything, and the pair (0, 0) is not one: a rotation
            ! of the zero column would put a zero on the diagonal of the
            ! triangular system and call the resulting residual zero. The
            ! identity keeps the recursion for gg defined; what may be
            ! returned is then decided by the rank test below, on the
            ! diagonal this leaves at zero.
            cs(j)  = 1.0d0
            sn(j)  = 0.0d0
         endif
         Hs(j,j)   = cs(j)*Hs(j,j) + sn(j)*Hs(j+1,j)
         Hs(j+1,j) = 0.0d0
         gg(j+1) = -sn(j)*gg(j)
         gg(j)   =  cs(j)*gg(j)
         ! |gg(j+1)| is the residual norm of the least-squares problem the
         ! rotations have reduced, i.e. of the step this subspace gives.
         res_now = abs(gg(j+1))
         ! --- THE RESIDUAL OF THIS SUBSPACE, REDUCED AND TRUE ---
         ! The reduced problem's residual after every product; every
         ! gm_history_stride products, and at the last one, the residual
         ! the iterate of this subspace actually reaches against the
         ! operator, which costs one triangular solve and one product. The
         ! two part company exactly where the recursion has stopped being a
         ! statement about the operator.
         if (history_here) then
            res_true_now = -1.0d0
            if (mod(j, gm_history_stride) .eq. 0 .or. j .eq. m) then
               hmin_diag = huge(1.0d0)
               do i = 1, j
                  hmin_diag = min(hmin_diag, abs(Hs(i,i)))
               enddo
               if (hmin_diag .gt. 0.0d0) then
                  yh(1:j) = 0.0d0
                  do i = j, 1, -1
                     tmp = gg(i)
                     do ib = i+1, j
                        tmp = tmp - Hs(i,ib)*yh(ib)
                     enddo
                     yh(i) = tmp/Hs(i,i)
                  enddo
                  xh = 0.0d0
                  do i = 1, j
                     xh = xh + yh(i)*V(:,i)
                  enddo
                  call banded_preconditioner_solve(abf, ipiv, neq, xh,   &
                                                   lpinfo)
                  if (lpinfo .eq. 0 .and.                                &
                      every_component_is_finite(xh)) then
                     call project_out_of_the_active_element_constraints( &
                                                                     xh)
                     call jacobian_action_of_direction(Y, F0, f_sp_base, &
                                                       D*xh, u, jv_ok)
                     if (jv_ok) then
                        if (nspec_row .gt. 0) then
                           u = u/Drow + idtau*(D/Drow)*xh
                        else
                           u = u/D + idtau*xh
                        endif
                        ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
                        ! the operator the cycle works with is E A_z and the band it is
                        ! preconditioned by is E times the same band, so the two agree
                        ! (unit_infinity_norm_row_scaling_of_the_band).
                        if (row_equil_here) u = u*model_row_equilibration
                        if (every_component_is_finite(u)) then
                           if (row_equil_here) then
                              res_true_now = sqrt(sum((b                 &
                                   *model_row_equilibration - u)**2))/beta
                           else
                              res_true_now = sqrt(sum((b - u)**2))/beta
                           endif
                        endif
                     endif
                  endif
               endif
            endif
            if (res_true_now .ge. 0.0d0) then
               write(*,'(A,I4,A,ES12.5,A,ES12.5)')                       &
                    ' (JFNK) [diag 16]   product ', j,                   &
                    ': reduced ', res_now/beta, ', TRUE ', res_true_now
            else
               write(*,'(A,I4,A,ES12.5)')                                &
                    ' (JFNK) [diag 16]   product ', j,                   &
                    ': reduced ', res_now/beta
            endif
         endif
         ! --- THE STEP THIS SUBSPACE GIVES, AND WHAT IT REALLY REACHES ---
         ! (gm_step_by_its_true_residual). From gm_true_residual_first
         ! products on, every gm_true_residual_stride and at the last
         ! product of the cycle: the iterate of the subspace built so far,
         ! and the residual it reaches AGAINST THE OPERATOR, at one
         ! triangular solve and one product. The best of them is kept with
         ! its own image, so the step that comes back and the model image
         ! the trust region reads on it are both the operator's and not the
         ! reduced problem's. The cycle stops where that residual has risen
         ! at gm_true_residual_rises consecutive checks, which is where the
         ! subspace has stopped buying anything the nonlinear action can
         ! deliver (section N35: the optimum lies at 60 to 80 products
         ! while the reduced residual keeps falling to 320).
         reduced_tolerance_stops = .true.
         if (true_resid_here) then
            ! THE CHECK POINTS ARE THE PERIODIC ONES AND EVERY POINT AT
            ! WHICH THE CYCLE WOULD OTHERWISE END: the last product of the
            ! subspace, an Arnoldi breakdown, and the reduced problem
            ! reaching the tolerance. With this arm the reduced residual
            ! stops nothing on its own -- it is the quantity the arm
            ! exists to stop trusting -- so where it reaches the tolerance
            ! the iterate is measured against the operator and the cycle
            ! ends only if the TRUE residual reached it too. Measuring the
            ! last iterate of every exit is also what keeps the step that
            ! comes back from being worse than the one the cycle had in
            ! its hand.
            res_true_now = -1.0d0
            check_here = ((j .ge. gm_true_residual_first) .and.          &
                          mod(j - gm_true_residual_first,                &
                              gm_true_residual_stride) .eq. 0)           &
                         .or. j .eq. m .or. arnoldi_breakdown            &
                         .or. (res_now .le. rtol*beta)
            if (check_here) then
               hmin_diag = huge(1.0d0)
               do i = 1, j
                  hmin_diag = min(hmin_diag, abs(Hs(i,i)))
               enddo
               if (hmin_diag .gt. 0.0d0) then
                  yh2(1:j) = 0.0d0
                  do i = j, 1, -1
                     tmp = gg(i)
                     do ib = i+1, j
                        tmp = tmp - Hs(i,ib)*yh2(ib)
                     enddo
                     yh2(i) = tmp/Hs(i,i)
                  enddo
                  u = 0.0d0
                  do i = 1, j
                     u = u + yh2(i)*V(:,i)
                  enddo
                  call banded_preconditioner_solve(abf, ipiv, neq, u,    &
                                                   lpinfo)
                  if (lpinfo .eq. 0 .and.                                &
                      every_component_is_finite(u)) then
                     call project_out_of_the_active_element_constraints( &
                                                                      u)
                     z = u
                     call jacobian_action_of_direction(Y, F0, f_sp_base, &
                                                       D*z, w, jv_ok)
                     if (jv_ok) then
                        if (nspec_row .gt. 0) then
                           w = w/Drow + idtau*(D/Drow)*z
                        else
                           w = w/D + idtau*z
                        endif
                        if (row_equil_here)                              &
                           w = w*model_row_equilibration
                        if (every_component_is_finite(w)) then
                           if (row_equil_here) then
                              res_true_now = sqrt(sum((b                 &
                                   *model_row_equilibration - w)**2))/beta
                           else
                              res_true_now = sqrt(sum((b - w)**2))/beta
                           endif
                           ! A CANDIDATE OUTSIDE THE TRUST BALL IS NOT AN
                           ! ADMISSIBLE STEP, and the ball's own truncation
                           ! below is what returns the boundary point.
                           if (on_the_ball) then
                              if (sqrt(sum(z*z)) .gt. delta_ball)        &
                                 res_true_now = -1.0d0
                           endif
                           if (res_true_now .ge. 0.0d0) then
                              if (res_true_now .lt. res_true_best) then
                                 res_true_best = res_true_now
                                 xbest  = z
                                 Axbest = w
                                 j_true_best = j
                                 have_a_best_step = .true.
                              endif
                              if (res_true_now .gt. res_true_last) then
                                 n_true_rise = n_true_rise + 1
                              else
                                 n_true_rise = 0
                              endif
                              res_true_last = res_true_now
                              if (res_true_now .le. rtol) then
                                 gm_outcome = gm_tolerance_reached
                                 exit
                              endif
                              if (n_true_rise .ge. gm_true_residual_rises) &
                                 then
                                 gm_outcome = gm_true_residual_turned
                                 exit
                              endif
                              ! The iterate was measured and it has not
                              ! reached the tolerance: the reduced
                              ! problem's claim that it has does not end
                              ! the cycle.
                              reduced_tolerance_stops = .false.
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
         ! --- HAS THE ITERATE LEFT THE TRUST BALL? (Steihaug-Toint) ---
         ! The iterate of the subspace built so far, and its image, are
         ! formed here from the rotated system and the Arnoldi columns: no
         ! product of the operator, one triangular solve of the
         ! preconditioner. Where the iterate is outside the ball, the point
         ! of the boundary on the segment from the last iterate that was
         ! inside is the step, and the same combination of the two images is
         ! its image. The GMRES iterate is not monotone in norm, so this is
         ! the FIRST crossing and not the last.
         if (on_the_ball) then
            hmin_diag = huge(1.0d0)
            do i = 1, j
               hmin_diag = min(hmin_diag, abs(Hs(i,i)))
            enddo
            ! A rank-deficient reduced problem has no iterate to measure
            ! against the ball; the rank test after the loop is what speaks
            ! about it, and the truncation stands aside.
            if (hmin_diag .gt. 0.0d0) then
               yb(1:j) = 0.0d0
               do i = j, 1, -1
                  tmp = gg(i)
                  do ib = i+1, j
                     tmp = tmp - Hs(i,ib)*yb(ib)
                  enddo
                  yb(i) = tmp/Hs(i,i)
               enddo
               u = 0.0d0
               do i = 1, j
                  u = u + yb(i)*V(:,i)
               enddo
               call banded_preconditioner_solve(abf, ipiv, neq, u, lpinfo)
               if (lpinfo .eq. 0 .and. every_component_is_finite(u)) then
                  xb = u
                  call project_out_of_the_active_element_constraints(xb)
                  Axb = 0.0d0
                  do jb = 1, j
                     do ib = 1, min(jb+1, kv_set)
                        if (Hb(ib,jb) .ne. 0.0d0)                         &
                           Axb = Axb + (Hb(ib,jb)*yb(jb))*V(:,ib)
                     enddo
                  enddo
                  if (sqrt(sum(xb*xb)) .gt. delta_ball) then
                     tball = the_point_where_the_step_leaves_the_ball(     &
                                  xb_in, xb, delta_ball)
                     x  = xb_in  + tball*(xb  - xb_in)
                     Ax = Axb_in + tball*(Axb - Axb_in)
                     stopped_on_the_ball = .true.
                     gm_outcome = gm_stopped_on_the_trust_ball
                     exit
                  endif
                  xb_in  = xb
                  Axb_in = Axb
               endif
            endif
         endif
         if (res_now .le. rtol*beta .and. reduced_tolerance_stops) then
            ! CONVERGENCE CLAIMED AT A BREAKDOWN IS A CLAIM ABOUT THE
            ! REDUCED PROBLEM ONLY. Verified against the operator below.
            verify_true_residual = arnoldi_breakdown
            exit
         endif
         if (arnoldi_breakdown) then
            ! The space is invariant and the tolerance was not reached: this
            ! subspace is all there is, and what its step achieves is a
            ! question for the operator.
            verify_true_residual = .true.
            exit
         endif
      enddo

      ! --- IS THE IMAGE THE ARNOLDI RELATION RETURNS THE OPERATOR'S? ---
      ! Column by column: A_z M^-1 V_k against V_{k+1} Hbar(:,k), the
      ! relation the trust region's predicted decrease is read from. A gap
      ! in ONE column says which direction of the recursion drifted, which
      ! the gap of the assembled step cannot. One direction is then sampled
      ! TWICE, so that a matrix-free action which is not reproducible is
      ! told apart from a basis that has lost orthogonality: the first
      ! disagreement is a property of the residual evaluation, the second of
      ! the recursion. One product of the operator per column; armed only by
      ! gm_verify_the_arnoldi_image.
      ! Only where the element diagnostic prints, which is the iterate the
      ! check is wanted at: one product of the operator per column is a
      ! multiple of the cycle's own cost and no solve carries it at every
      ! iteration.
      if (gm_verify_the_arnoldi_image .and. elem_diag_here .and.          &
          gm_iters .ge. 1) then
         gap_worst = 0.0d0
         do jb = 1, gm_iters
            z = V(:,jb)
            call banded_preconditioner_solve(abf, ipiv, neq, z, lpinfo)
            if (lpinfo .ne. 0) cycle
            call project_out_of_the_active_element_constraints(z)
            call jacobian_action_of_direction(Y, F0, f_sp_base, D*z, w,   &
                                              jv_ok)
            if (.not. jv_ok) cycle
            if (nspec_row .gt. 0) then
               w = w/Drow + idtau*(D/Drow)*z
            else
               w = w/D + idtau*z
            endif
            ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
            ! the operator the cycle works with is E A_z and the band it is
            ! preconditioned by is E times the same band, so the two agree
            ! (unit_infinity_norm_row_scaling_of_the_band).
            if (row_equil_here) w = w*model_row_equilibration
            u = 0.0d0
            do ib = 1, min(jb+1, kv_set)
               u = u + Hb(ib,jb)*V(:,ib)
            enddo
            gap_col = sqrt(sum((u - w)**2))                               &
                    / max(sqrt(sum(w*w)), 1.0d-300)
            gap_worst = max(gap_worst, gap_col)
            write(*,'(A,I3,A,ES11.3,A,ES11.3,A,ES11.3)')                  &
                 ' (JFNK) [image] column ', jb,                           &
                 ': relative gap of the Arnoldi image ', gap_col,         &
                 ', ||A_z M^-1 V_k|| ', sqrt(sum(w*w)),                    &
                 ', ||V Hbar_k|| ', sqrt(sum(u*u))
            if (jb .eq. 1) then
               ! The same direction, sampled a second time.
               call jacobian_action_of_direction(Y, F0, f_sp_base, D*z, u, &
                                                 jv_ok)
               if (jv_ok) then
                  if (nspec_row .gt. 0) then
                     u = u/Drow + idtau*(D/Drow)*z
                  else
                     u = u/D + idtau*z
                  endif
                  ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
                  ! the operator the cycle works with is E A_z and the band it is
                  ! preconditioned by is E times the same band, so the two agree
                  ! (unit_infinity_norm_row_scaling_of_the_band).
                  if (row_equil_here) u = u*model_row_equilibration
                  repro = sqrt(sum((u - w)**2))                           &
                        / max(sqrt(sum(w*w)), 1.0d-300)
                  write(*,'(A,ES11.3)') ' (JFNK) [image]   two products'//&
                       ' of the first direction differ by (relative) ',   &
                       repro
               endif
            endif
         enddo
         write(*,'(A,ES11.3,A,ES11.3)') ' (JFNK) [image] worst column'//  &
              ' gap ', gap_worst, ', loss of orthogonality of this'//     &
              ' basis ', gm_orthogonality_loss_cycle
         ! IS THE ACTION A LINEAR MAP OF THE DIRECTION AT ALL? The Arnoldi
         ! relation is a relation of a LINEAR operator, and the matrix-free
         ! action is a secant of the residual over an arc of fixed length
         ! (probe_length_of_the_jacobian_action): the quotient is
         ! homogeneous in the direction by construction, but it is the
         ! secant of a nonlinear residual in the direction sampled, so the
         ! action of a SUM of two directions need not be the sum of their
         ! actions. That difference is measured here on the first two
         ! directions of the basis, and it is the one property the column
         ! checks above cannot see: a column is compared with its own
         ! product, while the assembled step is a combination of them.
         if (gm_iters .ge. 2) then
            u = V(:,1) + V(:,2)
            call banded_preconditioner_solve(abf, ipiv, neq, u, lpinfo)
            if (lpinfo .eq. 0) then
               call project_out_of_the_active_element_constraints(u)
               call jacobian_action_of_direction(Y, F0, f_sp_base, D*u, w, &
                                                 jv_ok)
               if (jv_ok) then
                  if (nspec_row .gt. 0) then
                     w = w/Drow + idtau*(D/Drow)*u
                  else
                     w = w/D + idtau*u
                  endif
                  ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
                  ! the operator the cycle works with is E A_z and the band it is
                  ! preconditioned by is E times the same band, so the two agree
                  ! (unit_infinity_norm_row_scaling_of_the_band).
                  if (row_equil_here) w = w*model_row_equilibration
                  ! The sum of the two columns' own images, from the
                  ! Hessenberg the recursion wrote.
                  z = 0.0d0
                  do ib = 1, min(2, kv_set)
                     z = z + Hb(ib,1)*V(:,ib)
                  enddo
                  do ib = 1, min(3, kv_set)
                     z = z + Hb(ib,2)*V(:,ib)
                  enddo
                  write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                 &
                       ' (JFNK) [image] the action of the sum of the'//    &
                       ' first two directions departs from the sum of'//   &
                       ' their actions by (relative) ',                    &
                       sqrt(sum((w - z)**2))                               &
                       /max(sqrt(sum(w*w)), 1.0d-300),                     &
                       ', probe arc ',                                     &
                       probe_length_of_the_jacobian_action(Y),             &
                       ', ||Y|| ', sqrt(sum(Y*Y))
               endif
            endif
         endif
      endif

      ! --- IS THE ACTION ADDITIVE ACROSS DIRECTIONS? (jv_additivity_on) ---
      ! Measured on the first two directions of the basis, which is the
      ! pair the assembled step is built from, and on a pair that owes
      ! nothing to the basis or to the preconditioner, so that a defect of
      ! the operator is told from a defect the recursion selected for.
      if (jv_additivity_on .and. jv_additivity_here .and.                 &
          gm_iters .ge. 2) then
         z = V(:,1)
         call banded_preconditioner_solve(abf, ipiv, neq, z, lpinfo)
         if (lpinfo .eq. 0) then
            u = V(:,2)
            call banded_preconditioner_solve(abf, ipiv, neq, u, lpinfo)
         endif
         if (lpinfo .eq. 0) then
            call project_out_of_the_active_element_constraints(z)
            call project_out_of_the_active_element_constraints(u)
            call measure_the_additivity_of_the_jacobian_action(Y, F0,     &
                 f_sp_base, D, D*z, D*u,                                  &
                 'the first two directions of the basis')
         endif
         allocate(rv1(neq), rv2(neq))
         call deterministic_unit_direction(20260910, rv1)
         call deterministic_unit_direction(20260911, rv2)
         call measure_the_additivity_of_the_jacobian_action(Y, F0,        &
              f_sp_base, D, D*rv1, D*rv2, 'two directions off the basis')
         deallocate(rv1, rv2)
      endif

      ! --- WHERE DOES THE RESIDUAL STEP? (resid_jump_scan_on) ---
      ! Along the first direction of the basis, which is the one the
      ! Arnoldi relation is built on and the one whose probe arc is short
      ! enough for the floor to reach the quotient, and along one that
      ! owes nothing to the basis, so that the density of steps in an arc
      ! can be compared between the two.
      if (resid_jump_scan_on .and. resid_jump_scan_here .and.             &
          gm_iters .ge. 1) then
         z = V(:,1)
         call banded_preconditioner_solve(abf, ipiv, neq, z, lpinfo)
         if (lpinfo .eq. 0) then
            call project_out_of_the_active_element_constraints(z)
            call scan_the_residual_along_a_direction(Y, F0, f_sp_base,    &
                 D*z, 'the first direction of the basis')
         endif
         allocate(rv1(neq))
         call deterministic_unit_direction(20260910, rv1)
         call scan_the_residual_along_a_direction(Y, F0, f_sp_base,       &
              D*rv1, 'a direction off the basis')
         deallocate(rv1)
         ! Once for the solve: a restart cycle of the same iteration is
         ! the same state.
         resid_jump_scan_here = .false.
      endif

      ! --- the boundary point, where the cycle stopped on the ball ---
      ! The step and its image are already the point's; what is reported is
      ! the residual of THAT point, taken from the image the Arnoldi
      ! relation gives and so at no further product.
      if (stopped_on_the_ball) then
         if (every_component_is_finite(x)) then
            ! THE IMAGE LEAVES THE CYCLE ON THE CERTIFICATION SCALES, and
            ! so does the relative residual read from it, whatever rows
            ! the cycle itself worked in.
            if (row_equil_here) Ax = Ax/model_row_equilibration
            gm_resid_rel = scaled_two_norm(b - Ax)                        &
                         /max(beta_cert, 1.0d-300)
            gm_snorm     = scaled_two_norm(x)
         else
            x = 0.0d0;  Ax = 0.0d0
            gm_outcome   = gm_nonfinite
            gm_resid_rel = 1.0d0
            gm_snorm     = 0.0d0
         endif
         call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,     &
                                           yb,xb,Axb,xb_in,Axb_in)
         return
      endif

      ! --- the step the TRUE residual chose (gm_step_by_its_true_residual)
      ! What comes back is the candidate whose residual against the
      ! operator was smallest, its own image from the product that measured
      ! it, and that residual on the CERTIFICATION scales, which is the
      ! norm every reader of gm_resid_rel works in. The reduced problem's
      ! residual is not consulted: it is the quantity this arm exists to
      ! stop trusting.
      if (true_resid_here .and. have_a_best_step) then
         x = xbest
         if (row_equil_here) then
            if (present(Ax)) Ax = Axbest/model_row_equilibration
            gm_resid_rel = scaled_two_norm(b                             &
                 - Axbest/model_row_equilibration)                        &
                 /max(beta_cert, 1.0d-300)
         else
            if (present(Ax)) Ax = Axbest
            gm_resid_rel = scaled_two_norm(b - Axbest)                   &
                 /max(beta_cert, 1.0d-300)
         endif
         gm_snorm = scaled_two_norm(x)
         if (gm_outcome .ne. gm_true_residual_turned) then
            if (res_true_best .le. rtol) then
               gm_outcome = gm_tolerance_reached
            else
               gm_outcome = gm_subspace_exhausted
            endif
         endif
         if (gm_residual_history_on .and. gm_residual_history_here)      &
            write(*,'(A,I4,A,ES12.5,A,ES12.5)') ' (JFNK) [diag 16] the'//&
                 ' step of product ', j_true_best,                        &
                 ' was kept: true residual in the cycle norm ',           &
                 res_true_best, ', on the certification scales ',         &
                 gm_resid_rel
         call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,     &
                                           yb,xb,Axb,xb_in,Axb_in)
         return
      endif

      gm_resid_rel = res_now/beta
      if (gm_outcome .eq. gm_tolerance_reached .and.                      &
          res_now .gt. rtol*beta) gm_outcome = gm_subspace_exhausted

      ! --- the rank of the reduced problem, BEFORE any division by it ---
      ! The rotations have left an upper triangular matrix; a diagonal entry
      ! far below the largest one is an unknown the subspace says nothing
      ! about, and dividing by it manufactures a step out of rounding. The
      ! leading block whose diagonal stands above the threshold is
      ! nonsingular, and its solution is the least-squares step of that
      ! block with the remaining coefficients zero.
      kk   = gm_iters
      hmax = 0.0d0
      do i = 1, kk
         hmax = max(hmax, abs(Hs(i,i)))
      enddo
      krank = 0
      do i = 1, kk
         if (abs(Hs(i,i)) .le. gm_rank_relative*hmax) exit
         krank = i
      enddo
      if (krank .lt. kk) then
         if (gm_outcome .eq. gm_tolerance_reached .or.                    &
             gm_outcome .eq. gm_subspace_exhausted)                       &
            gm_outcome = gm_reduced_system_singular
         verify_true_residual = .true.
      endif
      ! Back-substitute H(1:krank,1:krank) yy = gg(1:krank)
      do i = krank, 1, -1
         tmp = gg(i)
         do j = i+1, krank
            tmp = tmp - Hs(i,j)*yy(j)
         enddo
         yy(i) = tmp/Hs(i,i)
      enddo
      ! u = V(:,1:kk) yy ;  x = M^{-1} u  (undo right preconditioning)
      u = 0.0d0
      do i = 1, kk
         u = u + yy(i)*V(:,i)
      enddo
      call banded_preconditioner_solve(abf, ipiv, neq, u, lpinfo)
      x = u
      ! THE STEP IS THE ONE THE OPERATOR WAS SAMPLED WITH. Every direction
      ! of the basis above was projected onto the active shared rows before
      ! the Jacobian was applied to it, so the image the recursion returns
      ! is the image of the PROJECTED combination; the leg handed back is
      ! projected too, and the two then describe one step. The projection is
      ! zero on every coordinate bound the step holds, so it commutes with
      ! the caller's hold.
      call project_out_of_the_active_element_constraints(x)
      if (lpinfo .ne. 0 .or. .not. every_component_is_finite(x)) then
         ! The step that comes back is not a step. Zero is, and it is what
         ! the caller is told, with the residual of the zero step.
         x = 0.0d0
         if (present(Ax)) Ax = 0.0d0
         if (lpinfo .ne. 0) then
            gm_outcome = gm_preconditioner_failed
         else
            gm_outcome = gm_nonfinite
         endif
         gm_resid_rel = 1.0d0
         gm_snorm     = 0.0d0
         call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,     &
                                           yb,xb,Axb,xb_in,Axb_in)
         return
      endif
      gm_snorm = scaled_two_norm(x)

      ! A_z x = A_z M^-1 (V_kk yy) = V_{kk+1} Hbar(1:kk+1,1:kk) yy.
      if (present(Ax)) then
         do j = 1, kk
            do i = 1, min(j+1, kv_set)
               if (Hb(i,j) .ne. 0.0d0) Ax = Ax + (Hb(i,j)*yy(j))*V(:,i)
            enddo
         enddo
      endif
      ! --- BACK TO THE CERTIFICATION SCALES, where the cycle worked in
      ! the equilibrated rows. The image is divided by E, and the relative
      ! residual handed back is read from THAT image against the
      ! certification right-hand side: the reduced problem's residual
      ! res_now/beta is a number in the cycle's own norm and no reader of
      ! gm_resid_rel works in it. Neither costs a product. Where the
      ! caller did not ask for the image, the same combination is assembled
      ! into the work vector w for the residual alone.
      if (row_equil_here) then
         if (present(Ax)) then
            Ax = Ax/model_row_equilibration
            gm_resid_rel = scaled_two_norm(b - Ax)/max(beta_cert,1.0d-300)
         else
            w = 0.0d0
            do j = 1, kk
               do i = 1, min(j+1, kv_set)
                  if (Hb(i,j) .ne. 0.0d0) w = w + (Hb(i,j)*yy(j))*V(:,i)
               enddo
            enddo
            gm_resid_rel = scaled_two_norm(b - w/model_row_equilibration) &
                         /max(beta_cert, 1.0d-300)
         endif
         if (elem_diag_here)                                             &
            write(*,'(A,ES12.5,A,ES12.5)') ' (JFNK) [row equil] the'//   &
                 ' cycle reached (equilibrated rows) ', res_now/beta,     &
                 ', on the certification scales ', gm_resid_rel
      endif

      ! --- what the step achieves against the operator itself ---
      ! One product of the scaled operator, taken only where the reduced
      ! problem is not a faithful statement about the system: the norm of
      ! the Arnoldi remainder underflows, or the reduced problem lost rank.
      ! A zero step costs nothing here, because the directional difference
      ! of a zero direction is zero without a residual evaluation
      ! (jv_product).
      if (verify_true_residual) then
         call jacobian_action_of_direction(Y, F0, f_sp_base, D*x, w, jv_ok)
         if (jv_ok) then
            if (nspec_row .gt. 0) then
               w = w/Drow + idtau*(D/Drow)*x
            else
               w = w/D + idtau*x
            endif
            ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
            ! the operator the cycle works with is E A_z and the band it is
            ! preconditioned by is E times the same band, so the two agree
            ! (unit_infinity_norm_row_scaling_of_the_band).
            if (row_equil_here) w = w*model_row_equilibration
            if (every_component_is_finite(w)) then
               if (row_equil_here) then
                  gm_resid_rel = scaled_two_norm(b                       &
                       - w/model_row_equilibration)                       &
                       /max(beta_cert, 1.0d-300)
               else
                  gm_resid_rel = scaled_two_norm(b - w)/beta
               endif
            endif
         endif
      endif

      call release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,        &
                                        yb,xb,Axb,xb_in,Axb_in)
      end subroutine pgmres

      ! ------------------------------------------------------!

      subroutine release_the_cycle_work_space(V,Hs,Hb,gg,cs,sn,yy,z,w,u,  &
                                              yb,xb,Axb,xb_in,Axb_in)
      ! The arrays one Krylov cycle holds, given back in one place: the
      ! truncated cycle carries five of them and the untruncated one does
      ! not allocate those, so every exit of the cycle releases what it has
      ! and nothing it has not.
      real*8, allocatable, intent(inout) :: V(:,:), Hs(:,:), Hb(:,:)
      real*8, allocatable, intent(inout) :: gg(:), cs(:), sn(:), yy(:)
      real*8, allocatable, intent(inout) :: z(:), w(:), u(:)
      real*8, allocatable, intent(inout) :: yb(:), xb(:), Axb(:)
      real*8, allocatable, intent(inout) :: xb_in(:), Axb_in(:)
      deallocate(V,Hs,Hb,gg,cs,sn,yy,z,w,u)
      if (allocated(yb))     deallocate(yb)
      if (allocated(xb))     deallocate(xb)
      if (allocated(Axb))    deallocate(Axb)
      if (allocated(xb_in))  deallocate(xb_in)
      if (allocated(Axb_in)) deallocate(Axb_in)
      end subroutine release_the_cycle_work_space

      ! ------------------------------------------------------!

      pure real*8 function the_point_where_the_step_leaves_the_ball(       &
                           x_in, x_out, delta) result(t)
      ! WHERE THE SEGMENT FROM AN ITERATE INSIDE THE BALL TO ONE OUTSIDE IT
      ! CROSSES THE BOUNDARY: the t of [0,1] with
      !   || x_in + t (x_out - x_in) || = delta ,
      ! the positive root of the quadratic
      !   ||d||^2 t^2 + 2 (x_in . d) t + ||x_in||^2 - delta^2 = 0 ,
      ! d = x_out - x_in. With ||x_in|| <= delta the constant term is not
      ! positive, so the root exists and is unique in [0, 1]. It is taken in
      ! the form -c/(b + sqrt(b^2 - a c)) rather than as the difference of
      ! two nearly equal numbers, so that a crossing close to the inner
      ! endpoint keeps its digits.
      real*8, dimension(:), intent(in) :: x_in, x_out
      real*8,               intent(in) :: delta
      real*8, dimension(size(x_in)) :: d
      real*8 :: aa, bb, cc, disc
      d  = x_out - x_in
      aa = sum(d*d)
      bb = sum(x_in*d)
      cc = sum(x_in*x_in) - delta*delta
      t  = 0.0d0
      if (aa .le. 0.0d0) return
      disc = bb*bb - aa*cc
      if (disc .lt. 0.0d0) disc = 0.0d0
      if (bb + sqrt(disc) .gt. 0.0d0) then
         t = -cc/(bb + sqrt(disc))
      else
         t = (-bb + sqrt(disc))/aa
      endif
      t = max(0.0d0, min(1.0d0, t))
      end function the_point_where_the_step_leaves_the_ball

      ! ------------------------------------------------------!

      subroutine pgmres_with_restarts(Y, F0, f_sp_base, D, Drow, abf,     &
                        ipiv, idtau, b, x, m, rtol, gm_iters, gm_outcome, &
                        gm_resid_rel, gm_snorm, Ax, delta_ball)
      ! GMRES(m) WITH RESTARTS: at most gm_restart_cycles cycles of the same
      ! subspace size, each started from the residual the one before it
      ! left.
      !
      ! The subspace size is untouched. A cycle of m products returns the
      ! minimizer of the residual over a Krylov space of dimension m; where
      ! that space is exhausted above the tolerance asked for, the residual
      ! it leaves is a right hand side like any other, and a second cycle of
      ! m from it costs the same m products with the same memory. The step
      ! is the sum of the cycles' steps and the operator's image of it is
      ! the sum of their images, both exactly, because the operator and the
      ! projection onto the active shared rows are linear and the same at
      ! every cycle (pgmres applies the projection to every direction it
      ! samples).
      !
      ! WHAT IS REPORTED is the relative residual of the SUM against the
      ! ORIGINAL right hand side, taken from the images the cycles hand
      ! back, so it costs no product; the cycle count is the total number of
      ! products; the outcome is the last cycle's, or "the requested linear
      ! tolerance was reached" where the sum reached it.
      !
      ! With gm_restart_cycles = 1 this is one call of pgmres and nothing
      ! else, which is the single cycle the solve has always taken.
      real*8, dimension(nvar_jac*N),  intent(in)  :: Y, F0, b, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in)  :: ipiv
      real*8,  intent(in)  :: idtau, rtol
      integer, intent(in)  :: m
      real*8, dimension(nvar_jac*N), intent(out) :: x
      real*8, dimension(nvar_jac*N), intent(out), optional :: Ax
      integer, intent(out) :: gm_iters, gm_outcome
      real*8,  intent(out) :: gm_resid_rel, gm_snorm
      ! The trust ball the step has to lie in, handed on unchanged. A cycle
      ! that stopped on the boundary has spent the region's whole length,
      ! so there is no residual to continue from inside it and the restarts
      ! do not run.
      real*8, intent(in), optional :: delta_ball

      integer :: icycle, it_cycle
      real*8  :: beta, resid_cycle, snorm_cycle
      real*8, dimension(nvar_jac*N) :: rhs, xc, Axc, Ax_sum

      call pgmres(Y, F0, f_sp_base, D, Drow, abf, ipiv, idtau, b, x, m,   &
                  rtol, gm_iters, gm_outcome, gm_resid_rel, gm_snorm, Ax, &
                  delta_ball)
      if (gm_restart_cycles .le. 1) return
      if (gm_outcome .eq. gm_stopped_on_the_trust_ball) return
      ! A cycle that returned no step, or one whose step is not a step, is
      ! not a residual to continue from: the outcome names why and a second
      ! cycle of the same operator would meet the same obstruction.
      if (gm_iters .le. 0) return
      if (gm_outcome .eq. gm_nonfinite .or.                               &
          gm_outcome .eq. gm_preconditioner_failed) return
      if (.not. present(Ax)) then
         ! Without the image of the step there is no residual to restart
         ! from that does not cost a product of its own, and the callers
         ! that do not ask for it are the ones that do not restart.
         return
      endif
      beta = sqrt(sum(b*b))
      if (.not. (beta .gt. 0.0d0)) return
      Ax_sum = Ax
      do icycle = 2, gm_restart_cycles
         if (gm_resid_rel .le. rtol) exit
         rhs = b - Ax_sum
         call pgmres(Y, F0, f_sp_base, D, Drow, abf, ipiv, idtau, rhs,    &
                     xc, m, rtol, it_cycle, gm_outcome, resid_cycle,      &
                     snorm_cycle, Axc)
         if (it_cycle .le. 0) exit
         if (.not. every_component_is_finite(xc)) exit
         x        = x + xc
         Ax_sum   = Ax_sum + Axc
         gm_iters = gm_iters + it_cycle
         gm_resid_rel = sqrt(sum((b - Ax_sum)**2))/beta
         gm_snorm     = scaled_two_norm(x)
      enddo
      Ax = Ax_sum
      if (gm_resid_rel .le. rtol) gm_outcome = gm_tolerance_reached
      end subroutine pgmres_with_restarts

      ! ------------------------------------------------------!

      subroutine krylov_cycle_over_subspace_sizes(Y, F0, f_sp_base, D,   &
                        Drow, abf, ipiv, idtau, b, rtol)
      ! IS IT THE SUBSPACE OR THE OPERATOR? The same cycle, at the same
      ! iterate, over a ladder of subspace sizes (krylov_scan_size), and
      ! the last two sizes again with the basis orthogonalized twice.
      !
      ! GMRES over a Krylov space of growing dimension reaches the exact
      ! solution in at most neq products, and it reaches a given relative
      ! residual in a number of products set by the clustering of the
      ! spectrum of the preconditioned operator: a cycle that converges once
      ! the subspace is large enough is an operator whose eigenvalues are
      ! merely spread, which the size of the space can beat. A cycle that
      ! stagnates at every size is an operator with a part the
      ! preconditioner does not touch, and no amount of subspace removes it.
      !
      ! None of these steps is adopted: the solve's own cycle has already
      ! run and returned the leg the trust region uses, and this scan is a
      ! measurement of the same linear system with a different truncation.
      ! Its residual samples are therefore held out of the solve's counts.
      real*8, dimension(nvar_jac*N),  intent(in) :: Y, F0, b, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in) :: ipiv
      real*8,  intent(in) :: idtau, rtol
      real*8, allocatable :: xs(:), w(:)
      integer :: neq, is, msize, gmit, outcome, ipass
      real*8  :: resid_rel, snorm, resid_true
      logical :: jv_ok, reortho_at_entry
      character(len=64) :: why_txt
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry

      neq = nvar_jac*N
      allocate(xs(neq), w(neq))
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      reortho_at_entry = gm_reorthogonalize
      write(*,'(A,ES11.3)') ' (JFNK) [diag 16] the same linear system'// &
           ' over a ladder of subspace sizes, ||b|| ', sqrt(sum(b*b))
      ! Two passes: the basis orthogonalized once, then twice on the two
      ! largest sizes, which is where a single Gram-Schmidt pass has had the
      ! most room to lose the orthogonality the reduced problem assumes.
      do ipass = 1, 2
         gm_reorthogonalize = (ipass .eq. 2)
         do is = 1, n_krylov_scan_size
            msize = krylov_scan_size(is)
            if (msize .gt. neq) cycle
            if (ipass .eq. 2 .and. is .le. 2) cycle
            gm_residual_history_here = .true.
            call pgmres(Y, F0, f_sp_base, D, Drow, abf, ipiv, idtau, b,  &
                        xs, msize, rtol, gmit, outcome, resid_rel, snorm)
            gm_residual_history_here = .false.
            ! WHAT THE RETURNED STEP REACHES AGAINST THE OPERATOR, not
            ! against the reduced problem: one product of the operator on
            ! the step the cycle hands back.
            resid_true = -1.0d0
            if (snorm .gt. 0.0d0) then
               call jacobian_action_of_direction(Y, F0, f_sp_base, D*xs, &
                                                 w, jv_ok)
               if (jv_ok) then
                  if (nspec_row .gt. 0) then
                     w = w/Drow + idtau*(D/Drow)*xs
                  else
                     w = w/D + idtau*xs
                  endif
                  if (every_component_is_finite(w))                      &
                     resid_true = sqrt(sum((b - w)**2))                  &
                                / max(sqrt(sum(b*b)), 1.0d-300)
               endif
            endif
            call gm_outcome_text(outcome, why_txt)
            write(*,'(A,I4,A,L1,A,I4,A,ES12.5,A,ES12.5,A,ES11.3)')       &
                 ' (JFNK) [diag 16] size ', msize,                       &
                 ', orthogonalized twice ', gm_reorthogonalize,          &
                 ': products ', gmit, ', reduced ', resid_rel,           &
                 ', TRUE ', resid_true, ', ||x|| ', snorm
            write(*,'(A,A)') ' (JFNK) [diag 16]   outcome: ',            &
                 trim(why_txt)
         enddo
      enddo
      gm_reorthogonalize = reortho_at_entry
      call put_back_solve_refusal_statistics(statistics_at_entry)
      call put_back_residual_evaluation_products(products_at_entry)
      deallocate(xs, w)
      end subroutine krylov_cycle_over_subspace_sizes

      ! ------------------------------------------------------!

      subroutine unknown_mask_of_the_row_class(kind_wanted, mask)
      ! ONE OF THE THREE SETS OF UNKNOWNS the spectrum is compressed onto:
      ! every unknown (0), the species rows, element and carrier alike (1),
      ! and the three hydrodynamic rows (2). The mask is the diagonal of an
      ! orthogonal projection, so a recursion that applies it after every
      ! product builds the spectrum of the COMPRESSION P A P restricted to
      ! the range of P, which is what "the operator on those rows" means.
      integer,                       intent(in)  :: kind_wanted
      real*8, dimension(nvar_jac*N), intent(out) :: mask
      integer :: i, j, k
      do i = 1, nvar_jac*N
         j = (i - 1)/nvar_jac + 1
         k = i - nvar_jac*(j - 1)
         select case (kind_wanted)
         case (1);  mask(i) = merge(1.0d0, 0.0d0, k .gt. 3)
         case (2);  mask(i) = merge(1.0d0, 0.0d0, k .le. 3)
         case default;  mask(i) = 1.0d0
         end select
      enddo
      end subroutine unknown_mask_of_the_row_class

      ! ------------------------------------------------------!

      subroutine row_class_shares_of_a_vector(v, share, jcell, ncell_hit)
      ! WHERE A VECTOR LIVES: the share of its squared norm carried by each
      ! of the five row classes (mass, momentum, energy, element, carrier),
      ! and the smallest set of CELLS that carries nine tenths of it, which
      ! is what says whether a direction is a localized front feature or a
      ! mode spread over the domain.
      real*8, dimension(nvar_jac*N), intent(in)  :: v
      real*8, dimension(5),          intent(out) :: share
      integer, dimension(10),        intent(out) :: jcell
      integer,                       intent(out) :: ncell_hit
      real*8, dimension(N) :: wcell
      real*8  :: total, running, wmax
      integer :: i, j, k, ic, jmax
      logical, dimension(N) :: taken
      share = 0.0d0;  wcell = 0.0d0
      do i = 1, nvar_jac*N
         j = (i - 1)/nvar_jac + 1
         k = i - nvar_jac*(j - 1)
         if (k .le. 3) then
            ic = k
         else if (srow_kind(k-3) .eq. srow_carrier) then
            ic = 5
         else
            ic = 4
         endif
         share(ic) = share(ic) + v(i)*v(i)
         wcell(j)  = wcell(j)  + v(i)*v(i)
      enddo
      total = sum(share)
      if (total .gt. 0.0d0) share = share/total
      ! The cells taken in order of the weight they carry, until nine
      ! tenths of the squared norm is accounted for.
      taken = .false.;  running = 0.0d0;  ncell_hit = 0;  jcell = 0
      do while (running .lt. 0.9d0*total .and. ncell_hit .lt. N)
         wmax = -1.0d0;  jmax = 0
         do j = 1, N
            if (taken(j)) cycle
            if (wcell(j) .gt. wmax) then
               wmax = wcell(j);  jmax = j
            endif
         enddo
         if (jmax .eq. 0 .or. wmax .le. 0.0d0) exit
         taken(jmax) = .true.
         running     = running + wmax
         ncell_hit   = ncell_hit + 1
         if (ncell_hit .le. 10) jcell(ncell_hit) = jmax
      enddo
      end subroutine row_class_shares_of_a_vector

      ! ------------------------------------------------------!

      subroutine ritz_values_of_a_matrix_free_operator(Y, F0, f_sp_base,&
                        D, Drow, abf, ipiv, idtau, mask, kmax, V, wr, wi,&
                        VR, ldvr, kdone, info)
      ! THE RITZ VALUES OF A_z M^-1 COMPRESSED ONTO THE RANGE OF A MASK.
      !
      ! An Arnoldi recursion of at most kmax products from a deterministic
      ! start inside the range of the mask, the mask applied after every
      ! product so that the recursion stays there, and the basis
      ! orthogonalized twice so that the Hessenberg is the operator's
      ! compression and not the recursion's rounding. The Ritz values are
      ! the eigenvalues of that Hessenberg (dgeev), and where kmax reaches
      ! the dimension of the space they are the eigenvalues of the operator
      ! itself.
      !
      ! The basis and the eigenvectors of the Hessenberg are handed back so
      ! that the caller can form Ritz vectors from them; info is dgeev's,
      ! and kdone the number of products the recursion managed.
      real*8, dimension(nvar_jac*N),  intent(in) :: Y, F0, D, Drow, mask
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in) :: ipiv
      real*8,  intent(in) :: idtau
      integer, intent(in) :: kmax, ldvr
      real*8, dimension(nvar_jac*N,kmax+1), intent(out) :: V
      real*8, dimension(kmax),      intent(out) :: wr, wi
      real*8, dimension(ldvr,kmax), intent(out) :: VR
      integer, intent(out) :: kdone, info
      real*8, allocatable :: Hs(:,:), Hd(:,:), z(:), w(:), work(:), VL(:,:)
      real*8  :: hnext, tmp
      integer :: neq, i, j, kk, lpinfo, lwork
      logical :: jv_ok, row_equil_here
      neq   = nvar_jac*N
      ! The spectrum measured is that of the operator the CYCLE runs, so
      ! it carries the row scaling of the linear model where that arm is
      ! on (the_linear_model_rows_are_equilibrated).
      row_equil_here = the_linear_model_rows_are_equilibrated(neq)
      kdone = 0;  info = 0
      wr = 0.0d0;  wi = 0.0d0;  VR = 0.0d0
      allocate(Hs(kmax+1,kmax), z(neq), w(neq))
      Hs = 0.0d0
      call deterministic_unit_direction(20260911, V(:,1))
      V(:,1) = mask*V(:,1)
      tmp = sqrt(sum(V(:,1)**2))
      if (tmp .le. 0.0d0) then
         deallocate(Hs, z, w)
         return
      endif
      V(:,1) = V(:,1)/tmp
      do j = 1, kmax
         z = V(:,j)
         call banded_preconditioner_solve(abf, ipiv, neq, z, lpinfo)
         if (lpinfo .ne. 0 .or. .not. every_component_is_finite(z)) exit
         call project_out_of_the_active_element_constraints(z)
         call jacobian_action_of_direction(Y, F0, f_sp_base, D*z, w,     &
                                           jv_ok)
         if (.not. jv_ok) exit
         if (nspec_row .gt. 0) then
            w = w/Drow + idtau*(D/Drow)*z
         else
            w = w/D + idtau*z
         endif
         ! AND THE ROW SCALING OF THE LINEAR MODEL, where the arm is on:
         ! the operator the cycle works with is E A_z and the band it is
         ! preconditioned by is E times the same band, so the two agree
         ! (unit_infinity_norm_row_scaling_of_the_band).
         if (row_equil_here) w = w*model_row_equilibration
         if (.not. every_component_is_finite(w)) exit
         w = mask*w
         do kk = 1, 2
            do i = 1, j
               tmp     = sum(w*V(:,i))
               w       = w - tmp*V(:,i)
               Hs(i,j) = Hs(i,j) + tmp
            enddo
         enddo
         hnext = sqrt(sum(w*w))
         kdone = j
         if (hnext .le. 1.0d2*epsilon(1.0d0)*max(maxval(abs(Hs(1:j,j))), &
                                                 1.0d-300)) exit
         if (j .eq. kmax) exit
         Hs(j+1,j) = hnext
         V(:,j+1)  = w/hnext
      enddo
      if (kdone .lt. 1) then
         deallocate(Hs, z, w)
         return
      endif
      lwork = max(8*kdone, 16)
      allocate(Hd(kdone,kdone), work(lwork), VL(1,1))
      Hd = Hs(1:kdone,1:kdone)
      call dgeev('N', 'V', kdone, Hd, kdone, wr, wi, VL, 1, VR, ldvr,    &
                 work, lwork, info)
      deallocate(Hd, work, VL, Hs, z, w)
      end subroutine ritz_values_of_a_matrix_free_operator

      ! ------------------------------------------------------!

      subroutine ritz_values_of_the_preconditioned_operator(Y, F0,       &
                        f_sp_base, D, Drow, abf, ipiv, idtau)
      ! THE SPECTRUM OF THE PRECONDITIONED OPERATOR A_z M^-1, measured.
      !
      ! GMRES over a subspace of dimension k reduces the residual by a
      ! factor bounded by the smallest value a polynomial of degree k that
      ! is one at the origin can take on the spectrum. A spectrum clustered
      ! away from the origin is beaten by a low-degree polynomial and the
      ! cycle converges in a few products; eigenvalues arbitrarily close to
      ! the origin force such a polynomial to stay near one there and no
      ! subspace of moderate size reduces anything. So the question the
      ! stalling cycle asks is answered by the eigenvalues near zero: how
      ! many there are and which unknowns their eigenvectors live on.
      !
      ! Run three times: on the whole operator, and on its compressions
      ! onto the species rows and onto the hydrodynamic rows, so that a
      ! cluster near zero can be attributed to a class of unknowns.
      real*8, dimension(nvar_jac*N),  intent(in) :: Y, F0, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in) :: ipiv
      real*8,  intent(in) :: idtau
      real*8, allocatable :: V(:,:), mask(:), yritz(:)
      real*8, allocatable :: wr(:), wi(:), VR(:,:)
      integer, allocatable :: iorder(:)
      real*8  :: tmp, amag, amin, amax
      real*8  :: share(5)
      integer :: neq, kdone, i, j, kk, ievinfo, iwhich, kmax
      integer :: n_below_2, n_below_4, n_complex, ipick, jpick
      integer :: jcell(10), ncell_hit, ic
      character(len=24) :: what
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry

      neq  = nvar_jac*N
      kmax = min(n_ritz_step, neq)
      if (kmax .lt. 2) return
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      allocate(V(neq,kmax+1), mask(neq), yritz(neq))
      allocate(wr(kmax), wi(kmax), VR(kmax,kmax), iorder(kmax))

      do iwhich = 0, 2
         call unknown_mask_of_the_row_class(iwhich, mask)
         select case (iwhich)
         case (1);  what = 'the species rows'
         case (2);  what = 'the hydrodynamic rows'
         case default;  what = 'every row'
         end select
         if (iwhich .eq. 1 .and. nspec_row .le. 0) cycle
         call ritz_values_of_a_matrix_free_operator(Y, F0, f_sp_base, D, &
                   Drow, abf, ipiv, idtau, mask, kmax, V, wr, wi, VR,    &
                   kmax, kdone, ievinfo)
         if (kdone .lt. 2 .or. ievinfo .ne. 0) then
            write(*,'(A,A,A,I0,A,I0)') ' (JFNK) [diag 15] ', trim(what), &
                 ': the recursion managed ', kdone, ' products, dgeev'// &
                 ' info ', ievinfo
            cycle
         endif
         ! The Ritz values ordered by magnitude, smallest first.
         do i = 1, kdone
            iorder(i) = i
         enddo
         do i = 1, kdone-1
            do j = i+1, kdone
               if (hypot(wr(iorder(j)), wi(iorder(j))) .lt.              &
                   hypot(wr(iorder(i)), wi(iorder(i)))) then
                  kk = iorder(i);  iorder(i) = iorder(j);  iorder(j) = kk
               endif
            enddo
         enddo
         amin = hypot(wr(iorder(1)), wi(iorder(1)))
         amax = hypot(wr(iorder(kdone)), wi(iorder(kdone)))
         n_below_2 = 0;  n_below_4 = 0;  n_complex = 0
         do i = 1, kdone
            amag = hypot(wr(i), wi(i))
            if (amag .lt. 1.0d-2) n_below_2 = n_below_2 + 1
            if (amag .lt. 1.0d-4) n_below_4 = n_below_4 + 1
            if (wi(i) .ne. 0.0d0) n_complex = n_complex + 1
         enddo
         write(*,'(A,A,A,I0,A)') ' (JFNK) [diag 15] the Ritz values'//   &
              ' of the preconditioned operator on ', trim(what), ' (',   &
              kdone, ' products):'
         write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                         &
              ' (JFNK) [diag 15]   largest magnitude ', amax,            &
              ', smallest ', amin, ', their ratio ',                     &
              amax/max(amin, 1.0d-300)
         write(*,'(A,I0,A,I0,A,I0,A,I0)')                                &
              ' (JFNK) [diag 15]   below 1e-2: ', n_below_2,             &
              ', below 1e-4: ', n_below_4, ', in complex pairs: ',       &
              n_complex, ' of ', kdone
         do i = 1, min(6, kdone)
            write(*,'(A,I2,A,ES12.5,A,ES12.5,A,ES12.5)')                 &
                 ' (JFNK) [diag 15]   smallest ', i, ': real ',          &
                 wr(iorder(i)), ', imaginary ', wi(iorder(i)),           &
                 ', magnitude ', hypot(wr(iorder(i)), wi(iorder(i)))
         enddo
         ! --- where the smallest Ritz vectors live ---
         ! Only for the whole operator: the compressions' vectors live on
         ! their own class by construction and say nothing new.
         if (iwhich .ne. 0) cycle
         if (allocated(ritz_vector_of_the_smallest))                     &
            deallocate(ritz_vector_of_the_smallest)
         if (allocated(ritz_value_of_the_smallest))                      &
            deallocate(ritz_value_of_the_smallest)
         allocate(ritz_vector_of_the_smallest(neq,3))
         allocate(ritz_value_of_the_smallest(3))
         ritz_vector_of_the_smallest = 0.0d0
         ritz_value_of_the_smallest  = 0.0d0
         do ipick = 1, min(3, kdone)
            jpick = iorder(ipick)
            ! A complex pair is stored by dgeev as two consecutive columns
            ! holding the real and the imaginary part of one eigenvector;
            ! the column indexed here is a vector of the two-dimensional
            ! invariant subspace, which is the direction wanted.
            yritz = 0.0d0
            do i = 1, kdone
               yritz = yritz + VR(i,jpick)*V(:,i)
            enddo
            tmp = sqrt(sum(yritz*yritz))
            if (tmp .gt. 0.0d0) yritz = yritz/tmp
            ritz_vector_of_the_smallest(:,ipick) = yritz
            ritz_value_of_the_smallest(ipick)    =                       &
                 hypot(wr(jpick), wi(jpick))
            call row_class_shares_of_a_vector(yritz, share, jcell,       &
                                              ncell_hit)
            write(*,'(A,I2,A,ES12.5,A)') ' (JFNK) [diag 15]   the Ritz'//&
                 ' vector of the smallest ', ipick, ' (magnitude ',      &
                 ritz_value_of_the_smallest(ipick), '):'
            do ic = 1, 5
               if (share(ic) .lt. 1.0d-3) cycle
               select case (ic)
               case (1);  what = 'mass'
               case (2);  what = 'momentum'
               case (3);  what = 'energy'
               case (4);  what = 'element'
               case default;  what = 'carrier'
               end select
               write(*,'(A,A,A,F8.4)') ' (JFNK) [diag 15]     ',         &
                    trim(what), ' rows carry a share ', share(ic)
            enddo
            write(*,'(A,I0,A,10(1X,I0))') ' (JFNK) [diag 15]     nine'// &
                 ' tenths of it sits in ', ncell_hit, ' cells, the'//    &
                 ' largest of them', (jcell(i), i = 1, min(10,           &
                 ncell_hit))
         enddo
      enddo

      deallocate(V, mask, yritz, wr, wi, VR, iorder)
      call put_back_solve_refusal_statistics(statistics_at_entry)
      call put_back_residual_evaluation_products(products_at_entry)
      end subroutine ritz_values_of_the_preconditioned_operator

      ! ------------------------------------------------------!

      pure subroutine column_split_against_the_band(col, icol, ab,       &
                        s_inside, s_band, s_outside, s_whole)
      ! ONE COLUMN OF THE SCALED JACOBIAN SPLIT AGAINST THE BAND, in
      ! squared norms: the part inside the band (s_inside), how far the
      ! band's own entries there are from it (s_band), the part outside
      ! the band, which the storage cannot hold (s_outside), and the whole
      ! difference between the column and the band's padded column
      ! (s_whole). The band is zero outside its half-width, so
      ! s_whole = s_band + s_outside identically, and reporting both is
      ! what makes the split a decomposition rather than two numbers.
      real*8, dimension(nvar_jac*N), intent(in)  :: col
      integer,                       intent(in)  :: icol
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab
      real*8, intent(out) :: s_inside, s_band, s_outside, s_whole
      integer :: irow
      real*8  :: bandv
      s_inside = 0.0d0;  s_band = 0.0d0;  s_outside = 0.0d0
      s_whole  = 0.0d0
      do irow = 1, nvar_jac*N
         if (abs(irow - icol) .le. kl_jac) then
            bandv    = ab(kl_jac+ku_jac+1 + irow - icol, icol)
            s_inside = s_inside + col(irow)*col(irow)
            s_band   = s_band   + (col(irow) - bandv)**2
            s_whole  = s_whole  + (col(irow) - bandv)**2
         else
            s_outside = s_outside + col(irow)*col(irow)
            s_whole   = s_whole   + col(irow)*col(irow)
         endif
      enddo
      end subroutine column_split_against_the_band

      ! ------------------------------------------------------!

      subroutine what_the_band_omits_of_the_jacobian(Y, F0, f_sp_base,   &
                        D, Drow, ab, abf, ipiv, idtau, r0)
      ! WHAT THE BANDED PRECONDITIONER MISSES, measured against the full
      ! finite-difference action at the same iterate.
      !
      ! The band holds |row - col| <= kl_jac, the stencil of the local
      ! operator, and its entries are a colored finite difference of the
      ! FULL residual (build_banded_jac_full), so two things separate it
      ! from the operator the Krylov cycle samples: the entries OUTSIDE the
      ! band, which the storage cannot hold at all and which are where the
      ! radiation's column coupling and the elemental budget's reach live,
      ! and the entries INSIDE it, which the coloring contaminates because
      ! columns of one color are not in fact disjoint once the residual is
      ! non-local.
      !
      ! Both are measured here. On the directions the spectrum named and on
      ! the first two Arnoldi directions of the solve's own right-hand side:
      ! the difference (A - A_band) v by row class and by cell. On the
      ! columns those directions live on: the true scaled column, split into
      ! the part inside the band, which is compared with the band's own
      ! entries, and the part outside it, which the band does not have; and
      ! the same column of the FROZEN-radiation residual, which is exactly
      ! banded by construction, so that whatever the full column carries
      ! outside the band is the non-local response and not the stencil.
      real*8, dimension(nvar_jac*N),  intent(in) :: Y, F0, D, Drow, r0
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in) :: ipiv
      real*8,  intent(in) :: idtau
      real*8, allocatable :: v(:), Av(:), Abv(:), dif(:), z(:), w(:)
      real*8, allocatable :: Fa(:), Fp(:), Yp(:), acol(:), fcol(:)
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heatw, coolw, npartw
      real*8  :: share(5), sh_dif(5)
      integer :: jcell(10), ncell_hit
      integer :: neq, i, j, k, iv, n_vec, lpinfo, icol, ncol_hit
      integer :: irow, ic, jbind, kbind
      integer, dimension(6) :: colpick
      real*8  :: tmp, nrm_full, nrm_dif, s_in, s_out, s_band, s_all
      real*8  :: s_whole
      real*8  :: hstep, worst, diag_here, off_hydro, off_elem, off_carr
      logical :: jv_ok, okp
      character(len=48) :: txt
      character(len=24) :: what
      type(residual_evaluation_products) :: products_at_entry
      type(solve_refusal_statistics)     :: statistics_at_entry

      neq = nvar_jac*N
      call hold_residual_evaluation_products(products_at_entry)
      call hold_solve_refusal_statistics(statistics_at_entry)
      allocate(v(neq), Av(neq), Abv(neq), dif(neq), z(neq), w(neq))
      allocate(Fa(neq), Fp(neq), Yp(neq), acol(neq), fcol(neq))

      ! --- the directions this is measured on ---
      ! Three Ritz vectors of the smallest Ritz values, where the spectrum
      ! hook has run at this iterate, and the first two Arnoldi directions
      ! of the cycle's own right-hand side, which are the pair the assembled
      ! step is mostly built from.
      n_vec = 2
      if (allocated(ritz_vector_of_the_smallest)) n_vec = 5

      do iv = 1, n_vec
         if (iv .le. 3 .and. allocated(ritz_vector_of_the_smallest)) then
            v = ritz_vector_of_the_smallest(:,iv)
            write(txt,'(A,I0)') 'the Ritz vector of the smallest ', iv
         else
            ! The Arnoldi directions, formed here: V1 is the normalized
            ! right-hand side, V2 the part of A_z M^-1 V1 orthogonal to it.
            tmp = sqrt(sum(r0*r0))
            if (tmp .le. 0.0d0) cycle
            v = -r0/tmp
            if (iv .eq. n_vec) then
               z = v
               call banded_preconditioner_solve(abf, ipiv, neq, z,       &
                                                lpinfo)
               if (lpinfo .ne. 0) cycle
               call project_out_of_the_active_element_constraints(z)
               call jacobian_action_of_direction(Y, F0, f_sp_base, D*z,  &
                                                 w, jv_ok)
               if (.not. jv_ok) cycle
               if (nspec_row .gt. 0) then
                  w = w/Drow + idtau*(D/Drow)*z
               else
                  w = w/D + idtau*z
               endif
               w = w - sum(w*v)*v
               tmp = sqrt(sum(w*w))
               if (tmp .le. 0.0d0) cycle
               v = w/tmp
               txt = 'the second Arnoldi direction'
            else
               txt = 'the first Arnoldi direction'
            endif
         endif
         if (sqrt(sum(v*v)) .le. 0.0d0) cycle
         ! --- (A - A_band) v ---
         call jacobian_action_of_direction(Y, F0, f_sp_base, D*v, Av,    &
                                           jv_ok)
         if (.not. jv_ok) cycle
         Av = Av/Drow
         call band_matvec(ab, v, Abv)
         dif = Av - Abv
         nrm_full = sqrt(sum(Av*Av))
         nrm_dif  = sqrt(sum(dif*dif))
         write(*,'(A,A)') ' (JFNK) [diag 17] ', trim(txt)
         call row_class_shares_of_a_vector(v, share, jcell, ncell_hit)
         call row_class_shares_of_a_vector(dif, sh_dif, jcell, ncell_hit)
         write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                         &
              ' (JFNK) [diag 17]   ||A v|| ', nrm_full,                  &
              ', ||(A - A_band) v|| ', nrm_dif, ', their ratio ',        &
              nrm_dif/max(nrm_full, 1.0d-300)
         do ic = 1, 5
            if (sh_dif(ic) .lt. 1.0d-3 .and. share(ic) .lt. 1.0d-3) cycle
            select case (ic)
            case (1);  what = 'mass'
            case (2);  what = 'momentum'
            case (3);  what = 'energy'
            case (4);  what = 'element'
            case default;  what = 'carrier'
            end select
            write(*,'(A,A,A,F8.4,A,F8.4)') ' (JFNK) [diag 17]     ',     &
                 trim(what), ' rows: share of the direction ', share(ic),&
                 ', share of the difference ', sh_dif(ic)
         enddo
         write(*,'(A,I0,A,10(1X,I0))') ' (JFNK) [diag 17]     nine'//    &
              ' tenths of the difference sits in ', ncell_hit,           &
              ' cells, the largest of them',                             &
              (jcell(i), i = 1, min(10, ncell_hit))

         ! --- the columns the direction lives on ---
         ! Up to six unknowns of largest weight in v: for each, the true
         ! scaled column against the band's own, inside and outside.
         colpick = 0;  ncol_hit = 0
         do k = 1, 6
            worst = -1.0d0;  icol = 0
            do i = 1, neq
               if (any(colpick(1:k-1) .eq. i)) cycle
               if (abs(v(i)) .gt. worst) then
                  worst = abs(v(i));  icol = i
               endif
            enddo
            if (icol .eq. 0 .or. worst .le. 0.0d0) exit
            colpick(k) = icol;  ncol_hit = k
         enddo
         if (iv .gt. 1) cycle
         ! The base point of the differences, and the heat, cool and
         ! particle count the frozen residual needs, from one evaluation of
         ! the iterate.
         call eval_residual(Y, f_sp_base, fwork, Fa, heatw, coolw,       &
                            npartw, admissible=okp,                      &
                            may_be_adopted=.false.,                      &
                            state_is_discarded=.true.,                   &
                            n_eq_sweeps_fixed=n_eq_sweeps_model)
         if (.not. okp) cycle
         call frozen_residual(Y, npartw, heatw, coolw, fcol)
         hstep = sqrt(epsilon(1.0d0))
         do k = 1, ncol_hit
            icol = colpick(k)
            Yp = Y
            Yp(icol) = Y(icol) + hstep*D(icol)
            call eval_residual(Yp, f_sp_base, fwork, Fp, heatw, coolw,   &
                               admissible=okp, may_be_adopted=.false.,   &
                               state_is_discarded=.true.,                &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            call name_of_unknown(icol, txt)
            if (.not. okp) then
               write(*,'(A,A,A)') ' (JFNK) [diag 17]   ', trim(txt),     &
                    ': the probe state is not describable, no column'
               cycle
            endif
            ! The SCALED column A(:,icol) = Drow^-1 (dF/dY) D, and the
            ! step was taken on D, so the column scale cancels the step
            ! and what is left is the difference over the step alone.
            acol = (Fp - Fa)/hstep/Drow
            ! The same column with the radiation held fixed, which is
            ! exactly banded: whatever the full column carries beyond the
            ! stencil and this one does not is the non-local response.
            call frozen_residual(Yp, npartw, heatw, coolw, Fp)
            w = (Fp - fcol)/hstep/Drow
            call column_split_against_the_band(acol, icol, ab, s_in,     &
                                               s_band, s_out, s_whole)
            s_all = s_in + s_out
            tmp = 0.0d0
            do irow = 1, neq
               if (abs(irow - icol) .gt. kl_jac) tmp = tmp + w(irow)**2
            enddo
            write(*,'(A,A)') ' (JFNK) [diag 17]   the column of ',       &
                 trim(txt)
            write(*,'(A,ES11.3,A,F9.5,A,F9.5)')                          &
                 ' (JFNK) [diag 17]     ||column|| ', sqrt(s_all),       &
                 ', outside the band ', sqrt(s_out)/max(sqrt(s_all),     &
                 1.0d-300), ', inside it, against the band''s entries ', &
                 sqrt(s_band)/max(sqrt(s_all), 1.0d-300)
            write(*,'(A,ES11.3,A,ES11.3)')                               &
                 ' (JFNK) [diag 17]     the frozen-radiation column'//   &
                 ' outside the band ', sqrt(tmp), ', of ',               &
                 sqrt(sum(w*w))
            ! THE SPLIT ACCOUNTS FOR THE WHOLE DIFFERENCE. Reported so that
            ! the two parts are read as a decomposition and not as two
            ! unrelated numbers.
            write(*,'(A,ES11.3,A,ES11.3)')                               &
                 ' (JFNK) [diag 17]     the two parts sum to ',          &
                 sqrt(s_band + s_out),                                   &
                 ' against the whole difference ', sqrt(s_whole)
         enddo
      enddo

      ! --- (c) THE BAND ROW OF THE ROW THAT BINDS ---
      ! Which unknowns the model believes the binding row responds to,
      ! grouped by class and read against the row's own diagonal: an element
      ! row whose coupling to the hydrodynamic rows of its own cell is a
      ! large multiple of its diagonal is a row the preconditioner cannot
      ! solve one unknown at a time.
      worst = -1.0d0;  irow = 1
      do i = 1, neq
         if (abs(r0(i)) .gt. worst) then
            worst = abs(r0(i));  irow = i
         endif
      enddo
      jbind = (irow - 1)/nvar_jac + 1
      kbind = irow - nvar_jac*(jbind - 1)
      call name_of_unknown(irow, txt)
      diag_here = ab(kl_jac+ku_jac+1, irow)
      off_hydro = 0.0d0;  off_elem = 0.0d0;  off_carr = 0.0d0
      do i = max(1, irow - kl_jac), min(neq, irow + ku_jac)
         if (i .eq. irow) cycle
         tmp = abs(ab(kl_jac+ku_jac+1 + irow - i, i))
         j = (i - 1)/nvar_jac + 1
         k = i - nvar_jac*(j - 1)
         if (k .le. 3) then
            off_hydro = off_hydro + tmp
         else if (srow_kind(k-3) .eq. srow_carrier) then
            off_carr = off_carr + tmp
         else
            off_elem = off_elem + tmp
         endif
      enddo
      write(*,'(A,A)') ' (JFNK) [diag 17] the band row of the row the'// &
           ' merit is largest on: ', trim(txt)
      write(*,'(A,ES12.5,A,ES12.5)') ' (JFNK) [diag 17]   diagonal ',    &
           diag_here, ', sum of the off-diagonal magnitudes ',           &
           off_hydro + off_elem + off_carr
      write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                            &
           ' (JFNK) [diag 17]   as a fraction of the diagonal:'//        &
           ' hydrodynamic ', off_hydro/max(abs(diag_here), 1.0d-300),    &
           ', element ', off_elem/max(abs(diag_here), 1.0d-300),         &
           ', carrier ', off_carr/max(abs(diag_here), 1.0d-300)

      deallocate(v, Av, Abv, dif, z, w)
      deallocate(Fa, Fp, Yp, acol, fcol)
      call put_back_solve_refusal_statistics(statistics_at_entry)
      call put_back_residual_evaluation_products(products_at_entry)
      end subroutine what_the_band_omits_of_the_jacobian

      ! ------------------------------------------------------!

      double precision function value_of_the_species_unknown(yv, isrow,   &
                                                    jcell) result(x)
      ! THE QUANTITY A SPECIES UNKNOWN NAMES, whichever variable the solve
      ! carries it in: the carrier density in code units, through the
      ! exponential where the carrier is carried in logarithm
      ! (carrier_density_from_unknown), and the element mass fraction or
      ! mixing ratio itself for an element row.
      real*8,  intent(in) :: yv
      integer, intent(in) :: isrow, jcell
      x = yv
      if (isrow .lt. 1 .or. isrow .gt. nspec_row) return
      if (srow_kind(isrow) .eq. srow_carrier)                             &
         x = carrier_density_from_unknown(yv, jcell)
      end function value_of_the_species_unknown

      ! ------------------------------------------------------!

      double precision function species_row_from_its_terms(t) result(rw)
      ! THE ROW ITS TERMS COMPOSE.  A transported species row is
      !     div(J_diffusive) + div(F_rho Y^face) - S_chemical ,
      ! so the three terms of a term set sum to the row with the chemical
      ! one carrying the minus sign.  It is stated once, here, so that the
      ! measurement and the test read one expression.
      type(species_row_term_set), intent(in) :: t
      rw = t%diffusive + t%adv_total - t%chemical
      end function species_row_from_its_terms

      ! ------------------------------------------------------!

      double precision function attributed_jacobian_entry(term_at_step,   &
                          term_at_base, hstep, to_code, dcol, drow)       &
                          result(a)
      ! ONE TERM'S CONTRIBUTION TO ONE ENTRY OF THE SCALED BANDED MODEL.
      !
      ! The band holds Drow^-1 J D, and the residual row is the operator's
      ! row written per code time in code density units (to_code), so the
      ! entry a term contributes to column jcol of row irow is the term's
      ! own one-sided difference quotient along that column, in the same
      ! two-sided scaling:
      !
      !     (term(Y + h e_jcol) - term(Y)) / h  x  to_code x D(jcol)
      !                                          / Drow(irow) .
      !
      ! With the step h that built the column (sqrt(eps) D(jcol)) the sum
      ! of the terms' quotients is the row's own, which is the band entry
      ! up to the contamination a coloured probe suffers from the other
      ! columns of its colour.
      real*8, intent(in) :: term_at_step, term_at_base, hstep, to_code
      real*8, intent(in) :: dcol, drow
      a = 0.0d0
      if (hstep .eq. 0.0d0 .or. drow .eq. 0.0d0) return
      a = (term_at_step - term_at_base)/hstep*to_code*dcol/drow
      end function attributed_jacobian_entry

      ! ------------------------------------------------------!

      subroutine terms_of_a_transported_species_row(isrow, jcell, rho,   &
                                                    Tst, f_sp, t)
      ! EVERY TERM OF ONE TRANSPORTED SPECIES ROW, on the state whose
      ! residual was evaluated last.
      !
      ! The operator that owns a row returns its terms summed.  They are
      ! separated here the way write_element_row_terms separates two of
      ! them: the SAME operator is evaluated again on the SAME state with
      ! the face mass flux MASKED, and the divergence of a face species
      ! flux is exactly zero wherever the mass flux it multiplies is
      ! (species_face_flux, species_flux_divergence).  Masking every face
      ! but one leaves that face's contribution standing alone, so the two
      ! faces of the cell are read separately and nothing is approximated.
      !
      ! The chemical source of a carrier is the same evaluation the row
      ! made, at the cell's own frozen background (carrier_source).  An
      ! element row has none: no reaction makes or destroys a nucleus, so
      ! the row is transport against transport.
      !
      ! WHAT THE CALLER MUST HAVE DONE.  The row values and their scales
      ! are the products of the last residual evaluation (crow_res_last,
      ! erow_he, erow_tr), and the face mass flux is the stored one of that
      ! same assembly, which face_mass_flux_of_state refuses to hand over
      ! for another state.  So this reads the state the caller evaluated
      ! and cannot silently read a different one.
      integer,                                intent(in)  :: isrow, jcell
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tst
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      type(species_row_term_set),             intent(out) :: t

      real*8, dimension(1-Ng:N+Ng) :: Frho, Fmask, msum
      real*8, dimension(1-Ng:N+Ng) :: ntot_c, nrho_c, wfac_c, TK_c
      real*8, dimension(1-Ng:N+Ng) :: mbar_c, nH_free, nO_free, nC_free
      real*8, dimension(1-Ng:N+Ng,n_carrier_max) :: fc
      real*8, dimension(1:N,n_carrier_max) :: adv, advmag
      real*8, dimension(1:N)         :: rhe_d, she_d
      real*8, dimension(1:N,n_melem) :: rtr_d, str_d
      logical, dimension(n_melem)    :: carried_d
      real*8  :: nc(n_carrier_max), src(n_carrier_max)
      real*8  :: fr_rp, fr_rm, fr_ap, fr_am, fr_dv, fr_tconv, fr_mc
      real*8  :: fr_cf
      logical :: ok_he_d, ok_tr_d
      integer :: ic, im, jd

      t%ok = .false.
      if (isrow .lt. 1 .or. isrow .gt. nspec_row) return
      if (jcell .lt. 1 .or. jcell .gt. N) return

      call face_mass_flux_of_state(rho, Frho)
      t%frho_left  = Frho(jcell-1)
      t%frho_right = Frho(jcell)
      fr_rp = r_edg(jcell)
      fr_rm = r_edg(jcell-1)
      fr_ap = fr_rp*fr_rp
      fr_am = fr_rm*fr_rm
      fr_dv = (fr_ap*fr_rp - fr_am*fr_rm)/3.0d0

      if (srow_kind(isrow) .eq. srow_carrier) then
         if (.not. allocated(crow_res_last)) return
         ic = srow_idx(isrow)
         call carrier_state(rho, f_sp, fc, ntot_c, nrho_c, wfac_c, TK_c, &
                            mbar_c, nH_free, nO_free, nC_free)
         call mixture_mass_sum(f_sp, msum)
         call carrier_advective_divergence(fc, msum, Frho, adv, advmag)
         t%adv_total = adv(jcell,ic)
         Fmask = 0.0d0;  Fmask(jcell) = Frho(jcell)
         call carrier_advective_divergence(fc, msum, Fmask, adv, advmag)
         t%adv_right = adv(jcell,ic)
         Fmask = 0.0d0;  Fmask(jcell-1) = Frho(jcell-1)
         call carrier_advective_divergence(fc, msum, Fmask, adv, advmag)
         t%adv_left = adv(jcell,ic)
         ! The face mass fraction each masked evaluation carried, read back
         ! out of its own divergence: the term is A F_rho Y^face / dV times
         ! the conversion, so dividing by everything but Y^face returns it.
         fr_tconv = n0*v0/R0
         fr_mc    = carrier_mass_amu(ic)
         fr_cf    = msum(jcell)/fr_mc*fr_tconv
         if (Frho(jcell) .ne. 0.0d0)                                     &
            t%yface_right = t%adv_right*fr_dv/(fr_ap*Frho(jcell)*fr_cf)
         if (Frho(jcell-1) .ne. 0.0d0)                                   &
            t%yface_left = -t%adv_left*fr_dv/(fr_am*Frho(jcell-1)*fr_cf)
         ! The donor cell of each face is the side the face MASS FLUX
         ! selects, which is the rule species_face_fraction upwinds on.
         jd = jcell;  if (Frho(jcell)   .lt. 0.0d0) jd = jcell+1
         t%ydonor_right = fr_mc*fc(min(jd,N+Ng),ic)/msum(min(jd,N+Ng))
         jd = jcell-1;  if (Frho(jcell-1) .lt. 0.0d0) jd = jcell
         t%ydonor_left  = fr_mc*fc(max(jd,1-Ng),ic)/msum(max(jd,1-Ng))
         t%adv_donor = (fr_ap*Frho(jcell)*t%ydonor_right                 &
                      - fr_am*Frho(jcell-1)*t%ydonor_left)/fr_dv*fr_cf
         nc = fc(jcell,:)*nrho_c(jcell)
         call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell),  &
                             src)
         t%chemical  = src(ic)
         t%row       = crow_res_last(jcell,ic)
         t%row_scale = crow_terms_last(jcell,ic)
         t%diffusive = t%row - t%adv_total + t%chemical
         t%density   = nc(ic)
         t%to_code   = R0/(v0*n0)
         t%ok        = .true.
         return
      endif

      ! --- an element row: transport against transport ---
      if (.not. allocated(erow_he)) return
      Fmask = 0.0d0
      call element_transport_residual(rho, Tst, f_sp, Fmask, rhe_d,      &
               she_d, ok_he_d, rtr_d, str_d, ok_tr_d, tr_carried=carried_d)
      if (.not. ok_he_d) return
      if (srow_kind(isrow) .eq. srow_element_he) then
         t%diffusive = rhe_d(jcell)
         t%row       = erow_he(jcell)
         t%row_scale = escale_he(jcell)
         t%to_code   = elem_he_to_code()
      else
         im = srow_idx(isrow)
         if (.not. ok_tr_d) return
         if (.not. erow_tr_carried(im)) return
         t%diffusive = rtr_d(jcell,im)
         t%row       = erow_tr(jcell,im)
         t%row_scale = escale_tr(jcell,im)
         t%to_code   = elem_tr_to_code()
      endif
      t%adv_total = t%row - t%diffusive
      Fmask = 0.0d0;  Fmask(jcell) = Frho(jcell)
      call element_transport_residual(rho, Tst, f_sp, Fmask, rhe_d,      &
               she_d, ok_he_d, rtr_d, str_d, ok_tr_d, tr_carried=carried_d)
      if (srow_kind(isrow) .eq. srow_element_he) then
         t%adv_right = rhe_d(jcell) - t%diffusive
      else
         t%adv_right = rtr_d(jcell,srow_idx(isrow)) - t%diffusive
      endif
      Fmask = 0.0d0;  Fmask(jcell-1) = Frho(jcell-1)
      call element_transport_residual(rho, Tst, f_sp, Fmask, rhe_d,      &
               she_d, ok_he_d, rtr_d, str_d, ok_tr_d, tr_carried=carried_d)
      if (srow_kind(isrow) .eq. srow_element_he) then
         t%adv_left = rhe_d(jcell) - t%diffusive
      else
         t%adv_left = rtr_d(jcell,srow_idx(isrow)) - t%diffusive
      endif
      t%chemical = 0.0d0
      t%ok       = .true.

      end subroutine terms_of_a_transported_species_row

      ! ------------------------------------------------------!

      subroutine the_binding_species_row(Y, f_sp, D, Drow, ab, r0)
      ! THE SPECIES ROW THAT BINDS, TERM BY TERM, AND ITS JACOBIAN ENTRIES
      ! ATTRIBUTED TO THOSE TERMS (front_row_on, EXHALE_FRONT_ROW=1).
      !
      ! The band row of the binding species row was measured to have a
      ! diagonal three decades below its coupling to the hydrodynamic
      ! unknowns of its own stencil.  A diagonal is the row's response to
      ! its own unknown; the row is a sum of terms; so the question is
      ! which term carries which response, and it is answered by the only
      ! device that can answer it, re-forming every term at a displaced
      ! state.
      !
      ! WHAT IS MEASURED.  For the binding species row and its two
      ! neighbors in the same species:
      !   (a) every term of the row on its certification scale, the two
      !       advective faces separately with the upwind side and the
      !       reconstructed face value beside the donor cell's own;
      !   (b) the band's entries for that row, to its own unknown, to the
      !       same unknown in the two neighbors and to the three
      !       hydrodynamic unknowns of the three cells, each as a number
      !       and as a multiple of the diagonal;
      !   (c) the same entries split by TERM: each term of the row is
      !       re-formed at Y and at Y displaced along one column by the
      !       step build_banded_jac_full uses for that column, and the
      !       difference quotients are reported.  Their sum is the row's
      !       own difference quotient, which is what the band entry is up
      !       to the contamination one coloured probe suffers from the
      !       other columns of its colour;
      !   (d) the local cell Peclet number of the row (its advective term
      !       against its diffusive one), the logarithmic gradient of the
      !       unknown across the cell, and, for a carrier, the ratio of
      !       the chemical time scale to the advective one.
      !
      ! COST: one residual evaluation for each column probed, four
      ! unknowns in each of 2*n_front_row_halfwidth+1 cells.  It adopts
      ! nothing: the products of every evaluation are held at entry and put
      ! back at exit, as the other measurements of this iterate do.
      real*8, dimension(nvar_jac*N),                   intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),          intent(in) :: f_sp
      real*8, dimension(nvar_jac*N),                   intent(in) :: D
      real*8, dimension(nvar_jac*N),                   intent(in) :: Drow
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: ab
      real*8, dimension(nvar_jac*N),                   intent(in) :: r0

      type(residual_evaluation_products) :: products_at_entry
      type(species_row_term_set) :: t0(3), tp(3)
      real*8, dimension(:), allocatable :: Yp, Fvec, heat, cool
      real*8, dimension(:,:), allocatable :: fwork
      real*8, dimension(:),   allocatable :: rho_s, T_s
      real*8, dimension(:,:), allocatable :: f_s
      character(len=64) :: txt, cname
      real*8  :: sqeps, hstep, band, dq_row, dq_adv, dq_dif, dq_chm
      real*8  :: dq_donor, diag_here, dscale, xm, xp, gradln
      real*8  :: tau_chem, tau_adv, peclet, worst
      integer :: fr_neq, i, k, jj, irow, ibind, jbind, isrow, ir
      integer :: jlo, jhi, jc, jcol, kk, m, nrows, jrow(3)
      logical :: cap_saved, ok_p

      if (.not. front_row_on) return
      if (nspec_row .le. 0) return
      fr_neq = nvar_jac*N

      ! --- the row that binds among the SPECIES rows ---
      ! The largest scaled row of the whole vector may be hydrodynamic
      ! (measured, the atomic reload), and this item is about the species
      ! rows, so the two are reported separately and the species one is
      ! the subject.
      worst = -1.0d0;  ibind = 0
      do i = 1, fr_neq
         jj = (i - 1)/nvar_jac + 1
         k  = i - nvar_jac*(jj - 1)
         if (k .le. 3) cycle
         if (abs(r0(i)) .gt. worst) then
            worst = abs(r0(i));  ibind = i
         endif
      enddo
      if (ibind .le. 0) return
      jbind = (ibind - 1)/nvar_jac + 1
      isrow = ibind - nvar_jac*(jbind - 1) - 3
      if (isrow .lt. 1 .or. isrow .gt. nspec_row) return
      call name_of_unknown(ibind, txt)
      worst = -1.0d0;  irow = 1
      do i = 1, fr_neq
         if (abs(r0(i)) .gt. worst) then
            worst = abs(r0(i));  irow = i
         endif
      enddo
      call name_of_unknown(irow, cname)
      write(*,'(A,A,A,ES12.5)') ' (JFNK) [diag 18] the species row that'//&
           ' binds: ', trim(txt), ', scaled residual ',                  &
           abs(r0(ibind))
      write(*,'(A,A,A,ES12.5)') ' (JFNK) [diag 18]   the largest scaled'//&
           ' row of the whole vector: ', trim(cname), ', ',              &
           abs(r0(irow))

      nrows = 0
      do m = -1, 1
         jc = jbind + m
         if (jc .lt. 1 .or. jc .gt. N) cycle
         nrows = nrows + 1
         jrow(nrows) = jc
      enddo

      allocate(Yp(fr_neq), Fvec(fr_neq))
      allocate(heat(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      allocate(fwork(1-Ng:N+Ng,n_species))
      allocate(rho_s(1-Ng:N+Ng), T_s(1-Ng:N+Ng))
      allocate(f_s(1-Ng:N+Ng,n_species))

      call hold_residual_evaluation_products(products_at_entry)
      cap_saved = resid_capture_operator_state
      resid_capture_operator_state = .true.

      ! --- (a) the terms of the three rows at the iterate ---
      call eval_residual(Y, f_sp, fwork, Fvec, heat, cool,               &
                         n_eq_sweeps_fixed=n_eq_sweeps_model)
      rho_s = captured_state_with_ghosts(1,:)
      T_s   = captured_temperature
      f_s   = fwork
      do ir = 1, nrows
         call terms_of_a_transported_species_row(isrow, jrow(ir), rho_s, &
                                                 T_s, f_s, t0(ir))
         t0(ir)%unknown = Y(nvar_jac*(jrow(ir)-1)+3+isrow)
      enddo

      do ir = 1, nrows
         jc = jrow(ir)
         i  = nvar_jac*(jc-1) + 3 + isrow
         call name_of_unknown(i, txt)
         if (.not. t0(ir)%ok) then
            write(*,'(A,A,A)') ' (JFNK) [diag 18] ', trim(txt),          &
                 ': the operator posed no row here'
            cycle
         endif
         dscale = max(t0(ir)%row_scale, cert_scale_floor)
         write(*,'(A,A,A,F9.6,A,ES13.6,A,ES12.5)')                       &
              ' (JFNK) [diag 18] ', trim(txt), ': r ', r(jc),            &
              ', unknown ', t0(ir)%unknown, ', row scale ',              &
              t0(ir)%row_scale
         write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)')                &
              ' (JFNK) [diag 18]   on the certification scale: row ',    &
              t0(ir)%row/dscale, ', diffusive ',                         &
              t0(ir)%diffusive/dscale, ', advective ',                   &
              t0(ir)%adv_total/dscale, ', chemical ',                    &
              -t0(ir)%chemical/dscale
         ! The first-order upwind form of the advective term is formed
         ! for a CARRIER row only: it is the reconstruction of the face
         ! composition that the limiter acts on, and the carrier row is
         ! the one whose limiter is in question here.
         if (srow_kind(isrow) .eq. srow_carrier) then
            write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                      &
                 ' (JFNK) [diag 18]   advective faces: inner ',          &
                 t0(ir)%adv_left/dscale, ', outer ',                     &
                 t0(ir)%adv_right/dscale, ', first-order upwind total ', &
                 t0(ir)%adv_donor/dscale
         else
            write(*,'(A,ES12.5,A,ES12.5)')                               &
                 ' (JFNK) [diag 18]   advective faces: inner ',          &
                 t0(ir)%adv_left/dscale, ', outer ',                     &
                 t0(ir)%adv_right/dscale
         endif
         write(*,'(A,ES12.5,A,ES12.5)')                                  &
              ' (JFNK) [diag 18]   face mass flux: inner ',              &
              t0(ir)%frho_left, ', outer ', t0(ir)%frho_right
         if (srow_kind(isrow) .eq. srow_carrier) then
            write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)')             &
                 ' (JFNK) [diag 18]   face fraction reconstructed:'//    &
                 ' inner ', t0(ir)%yface_left, ' against donor ',        &
                 t0(ir)%ydonor_left, ', outer ', t0(ir)%yface_right,     &
                 ' against donor ', t0(ir)%ydonor_right
            tau_chem = 0.0d0;  tau_adv = 0.0d0
            if (t0(ir)%chemical .ne. 0.0d0)                              &
               tau_chem = t0(ir)%density/abs(t0(ir)%chemical)
            if (t0(ir)%adv_total .ne. 0.0d0)                             &
               tau_adv = t0(ir)%density/abs(t0(ir)%adv_total)
            write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                      &
                 ' (JFNK) [diag 18]   time scales [s]: chemical ',       &
                 tau_chem, ', advective ', tau_adv, ', ratio ',          &
                 tau_chem/max(tau_adv, 1.0d-300)
         endif
         peclet = abs(t0(ir)%adv_total)/max(abs(t0(ir)%diffusive),       &
                                            1.0d-300)
         ! The logarithmic gradient of the unknown across the cell, taken
         ! on the DENSITY the unknown names, so that it is one number
         ! whether the carrier is carried as n or as ln n.
         xm = 0.0d0;  xp = 0.0d0;  gradln = 0.0d0
         if (jc .gt. 1)  xm = value_of_the_species_unknown(               &
                                 Y(nvar_jac*(jc-2)+3+isrow), isrow, jc-1)
         if (jc .lt. N)  xp = value_of_the_species_unknown(               &
                                 Y(nvar_jac*jc+3+isrow), isrow, jc+1)
         if (xm .gt. 0.0d0 .and. xp .gt. 0.0d0)                          &
            gradln = (log(xp) - log(xm))                                 &
                    /max(log(r(min(jc+1,N))) - log(r(max(jc-1,1))),      &
                         1.0d-300)
         write(*,'(A,ES12.5,A,ES12.5)')                                  &
              ' (JFNK) [diag 18]   cell Peclet (advective over'//        &
              ' diffusive) ', peclet, ', dln(unknown)/dlnr ', gradln
         ! The two identities the split has to satisfy on the state, and
         ! they are printed rather than asserted: the terms compose the
         ! row, and the two masked faces compose the advective term.
         write(*,'(A,ES11.3,A,ES11.3)')                                  &
              ' (JFNK) [diag 18]   the terms compose the row to ',       &
              (species_row_from_its_terms(t0(ir)) - t0(ir)%row)/dscale,  &
              ', the two faces compose the advective term to ',          &
              (t0(ir)%adv_left + t0(ir)%adv_right                        &
               - t0(ir)%adv_total)/dscale
      enddo

      ! --- (b) and (c): the band's entries and their term split ---
      jlo = max(1, jbind - n_front_row_halfwidth)
      jhi = min(N, jbind + n_front_row_halfwidth)
      sqeps = sqrt(epsilon(1.0d0))
      do jc = jlo, jhi
         do kk = 1, 4
            k = kk
            if (kk .eq. 4) k = 3 + isrow
            jcol  = nvar_jac*(jc-1) + k
            hstep = sqeps*D(jcol)
            if (hstep .eq. 0.0d0) cycle
            Yp = Y
            Yp(jcol) = Y(jcol) + hstep
            call eval_residual(Yp, f_sp, fwork, Fvec, heat, cool,        &
                               admissible=ok_p, may_be_adopted=.false.,  &
                               n_eq_sweeps_fixed=n_eq_sweeps_model)
            do ir = 1, nrows
               call terms_of_a_transported_species_row(isrow, jrow(ir),  &
                        captured_state_with_ghosts(1,:),                 &
                        captured_temperature, fwork, tp(ir))
            enddo
            call name_of_unknown(jcol, cname)
            do ir = 1, nrows
               i = nvar_jac*(jrow(ir)-1) + 3 + isrow
               if (i - jcol .gt. kl_jac .or. jcol - i .gt. ku_jac) cycle
               if (.not. (t0(ir)%ok .and. tp(ir)%ok)) cycle
               band = ab(kl_jac+ku_jac+1 + i - jcol, jcol)
               diag_here = ab(kl_jac+ku_jac+1, i)
               dq_row   = attributed_jacobian_entry(tp(ir)%row,          &
                             t0(ir)%row, hstep, t0(ir)%to_code, D(jcol), &
                             Drow(i))
               dq_dif   = attributed_jacobian_entry(tp(ir)%diffusive,    &
                             t0(ir)%diffusive, hstep, t0(ir)%to_code,    &
                             D(jcol), Drow(i))
               dq_adv   = attributed_jacobian_entry(tp(ir)%adv_total,    &
                             t0(ir)%adv_total, hstep, t0(ir)%to_code,    &
                             D(jcol), Drow(i))
               dq_chm   = attributed_jacobian_entry(tp(ir)%chemical,     &
                             t0(ir)%chemical, hstep, t0(ir)%to_code,     &
                             D(jcol), Drow(i))
               dq_donor = attributed_jacobian_entry(tp(ir)%adv_donor,    &
                             t0(ir)%adv_donor, hstep, t0(ir)%to_code,    &
                             D(jcol), Drow(i))
               call name_of_unknown(i, txt)
               write(*,'(A,A,A,A)') ' (JFNK) [diag 18] row ', trim(txt), &
                    ' against the column of ', trim(cname)
               write(*,'(A,ES12.5,A,ES11.3,A,ES12.5)')                   &
                    ' (JFNK) [diag 18]   band entry ', band,             &
                    ', as a multiple of the diagonal ',                  &
                    band/max(abs(diag_here), 1.0d-300),                  &
                    ', the row re-formed ', dq_row
               write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)')          &
                    ' (JFNK) [diag 18]   by term: diffusive ', dq_dif,   &
                    ', advective ', dq_adv, ', chemical ', -dq_chm,      &
                    ', their sum ', dq_dif + dq_adv - dq_chm
               if (srow_kind(isrow) .eq. srow_carrier)                   &
                  write(*,'(A,ES12.5,A,ES12.5)')                         &
                       ' (JFNK) [diag 18]   the advective term with'//   &
                       ' the reconstruction frozen at first order ',     &
                       dq_donor, ', against the limited one ', dq_adv
            enddo
         enddo
      enddo

      ! --- what the diagonal would be on another column scale ---
      ! The diagonal is dF/dY times D(column)/Drow(row), so a column scale
      ! raised or lowered by a decade moves it by the same decade while
      ! every entry of the row to ANOTHER column stands where it was: the
      ! ratio the preconditioner has to stand on is a choice of scale and
      ! not a property of the operator.  What the row's response to its own
      ! unknown is per unit of the PHYSICAL density is the number that does
      ! not move, and it is reported beside it.
      i = nvar_jac*(jbind-1) + 3 + isrow
      diag_here = ab(kl_jac+ku_jac+1, i)
      write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                            &
           ' (JFNK) [diag 18] the diagonal of the binding row ',         &
           diag_here, ', a decade of column scale up ', diag_here*10.0d0,&
           ', down ', diag_here*0.1d0
      write(*,'(A,ES12.5,A,ES12.5)')                                     &
           ' (JFNK) [diag 18]   its column scale ', D(i),                &
           ', its row scale ', Drow(i)
      if (D(i) .gt. 0.0d0)                                               &
         write(*,'(A,ES12.5)')                                           &
              ' (JFNK) [diag 18]   the same response per unit of the'//  &
              ' unknown in code units ', diag_here/D(i)
      ! THE COLUMN SCALE OF EVERY UNKNOWN OF THE BINDING CELL BESIDE THE
      ! UNKNOWN ITSELF.  An entry of the scaled band is the row's response
      ! to a displacement of one COLUMN SCALE of its unknown, so a column
      ! whose scale stands far above the unknown it scales is read at a
      ! displacement the state does not have: the momentum unknown's scale
      ! is |rho v| + rho c_s (state_scales_of_cell), which in a deeply
      ! subsonic layer is the sound speed and not the wind, so it exceeds
      ! the momentum itself by one over the Mach number.
      do kk = 1, 4
         k = kk
         if (kk .eq. 4) k = 3 + isrow
         jcol = nvar_jac*(jbind-1) + k
         call name_of_unknown(jcol, cname)
         write(*,'(A,A,A,ES12.5,A,ES12.5,A,ES11.3,A,ES12.5)')            &
              ' (JFNK) [diag 18]   ', trim(cname), ': unknown ',         &
              Y(jcol), ', column scale ', D(jcol),                       &
              ', scale over unknown ',                                   &
              D(jcol)/max(abs(Y(jcol)), 1.0d-300),                       &
              ', row scale ', Drow(jcol)
      enddo

      resid_capture_operator_state = cap_saved
      call put_back_residual_evaluation_products(products_at_entry)
      deallocate(Yp, Fvec, heat, cool, fwork, rho_s, T_s, f_s)

      end subroutine the_binding_species_row

      ! ------------------------------------------------------!

      subroutine initial_trust_region_radius(gm_outcome, dZ, g_held,      &
                                             tr_delta, tr_dmax, status)
      ! THE RADIUS THE TRUST REGION STARTS FROM, TAKEN FROM A DIRECTION THE
      ! STEP MAY TAKE.
      !
      ! The Krylov step is the direction of choice, because the region is
      ! sized against the arc over which the Newton model was measured to
      ! describe the merit (tr_first_radius_fraction). It may be used only
      ! where the cycle's own outcome says it is a step: a cycle that could
      ! not sample a direction returns the zero vector, and a hundredth of
      ! zero is not a small radius but the absence of one. A ball of radius
      ! zero contains a single point, so an iteration handed one can report
      ! nothing but a zero step, and the state is absorbing unless the
      ! initialization is guarded by "is this a positive finite length"
      ! rather than by a sign test (R1).
      !
      ! Where the Krylov step is not available the PROJECTED GRADIENT is:
      ! it has been zeroed on every unknown held at a bound
      ! (fix_active_species_bounds), so it points nowhere out of the
      ! feasible set and a step along it always has room inside. That is the
      ! direction the dogleg then degenerates to, and it is the direction
      ! every convergence proof of a trust-region method rests on (Nocedal
      ! and Wright, Numerical Optimization, 2nd ed., section 4.1).
      !
      ! With neither direction there is nothing to size, and saying so is
      ! the answer: an arbitrary positive radius does not create a direction
      ! that is missing.
      !
      ! The lengths are in the SCALED coordinates, where the column scales
      ! are already divided out, so no scale enters here.
      integer,                       intent(in)  :: gm_outcome
      real*8, dimension(nvar_jac*N), intent(in)  :: dZ, g_held
      real*8,                        intent(out) :: tr_delta, tr_dmax
      integer,                       intent(out) :: status
      real*8 :: nz, ng, dmin
      tr_delta = 0.0d0
      tr_dmax  = -1.0d0
      nz = scaled_two_norm(dZ)
      ng = scaled_two_norm(g_held)
      dmin = tr_delta_min_dimensionless*sqrt(dble(nvar_jac*N))
      if ((gm_outcome .eq. gm_tolerance_reached  .or.                     &
           gm_outcome .eq. gm_subspace_exhausted .or.                     &
           gm_outcome .eq. gm_reduced_system_singular) .and.              &
          nz .gt. 0.0d0 .and. finite_real(nz)) then
         tr_delta = tr_first_radius_fraction*nz
         status   = tr_radius_from_the_krylov_step
      else if (ng .gt. 0.0d0 .and. finite_real(ng)) then
         tr_delta = tr_first_radius_fraction*ng
         status   = tr_radius_from_the_gradient
      else
         status   = tr_radius_no_feasible_descent
         return
      endif
      tr_delta = max(tr_delta, dmin)
      tr_dmax  = min(tr_radius_ceiling_factor*tr_delta,                   &
                     trust_region_radius_absolute_cap())
      tr_delta = min(tr_delta, tr_dmax)
      end subroutine initial_trust_region_radius

      ! ------------------------------------------------------!

      double precision function trust_region_radius_absolute_cap()        &
                               result(dcap)
      ! THE LARGEST RADIUS THAT IS STILL A CONTROL, in the scaled
      ! coordinates the step lives in.
      !
      ! There the column scales are divided out, so a step of Euclidean
      ! length sqrt(neq) changes EVERY unknown by one unit of its own
      ! characteristic scale, which is the largest step the state's own
      ! magnitudes admit as a single move. The cap is a thousand of those:
      ! far above anything a step takes, so it binds on no ordinary
      ! iteration, and still finite, so the ceiling cannot run away.
      !
      ! WHY A CAP AT ALL. The ceiling doubles on every accepted step whose
      ! ratio exceeded three quarters and had nothing above it: MEASURED on
      ! the atomic element reload, it reached delta = 2.684e8 against a step
      ! of ||s|| = 1.018e-5, eleven decades apart, and the dogleg was
      ! returning the interior Krylov point throughout (N7b, noticed item
      ! 4). A radius eleven decades above the step it bounds bounds nothing
      ! -- the admissibility cuts and the species box are what limit the
      ! step there -- and it makes the printed radius useless as a
      ! diagnostic. The feasible set is still what bounds the step; this
      ! bounds the number that claims to.
      real*8, parameter :: tr_radius_absolute_cap_factor = 1.0d3
      dcap = tr_radius_absolute_cap_factor*sqrt(dble(nvar_jac*N))
      end function trust_region_radius_absolute_cap

      ! ------------------------------------------------------!

      subroutine grow_the_trust_region_ceiling(tr_dmax)
      ! The ceiling follows the radius: an accepted step whose actual
      ! reduction reached three quarters of what the model promised is a
      ! model good over the whole radius, so the region may double next
      ! iteration and the ceiling with it. The absolute cap is what it may
      ! not pass (trust_region_radius_absolute_cap).
      real*8, intent(inout) :: tr_dmax
      tr_dmax = min(2.0d0*tr_dmax, trust_region_radius_absolute_cap())
      end subroutine grow_the_trust_region_ceiling

      ! ------------------------------------------------------!

      subroutine tr_radius_status_text(code, txt)
      ! The sentence that belongs to one outcome of the radius
      ! initialization.
      integer,          intent(in)  :: code
      character(len=*), intent(out) :: txt
      select case (code)
      case (tr_radius_from_the_krylov_step)
         txt = 'the first radius was measured on the Krylov step'
      case (tr_radius_from_the_gradient)
         txt = 'the first radius was measured on the projected gradient'
      case (tr_radius_no_feasible_descent)
         txt = 'no feasible descent direction exists at this iterate'
      case default
         txt = 'unnamed'
      end select
      end subroutine tr_radius_status_text

      ! ------------------------------------------------------!

      pure subroutine cauchy_length_along_the_banded_gradient(rAg, nAg,    &
                                                              tau, ascends)
      ! HOW FAR ALONG -g THE MODEL 1/2||r0 + A s||^2 IS MINIMIZED, and
      ! whether it falls along -g at all.
      !
      ! g is the transpose of the BANDED A applied to r0, because a
      ! transpose is not available matrix-free, so it is not the gradient of
      ! the model it is used in and nothing makes -g a descent direction of
      ! that model. The slope of the model along -g at s = 0 is -(r0 . A g),
      ! so the model falls only while r0 . A g is positive. With the exact
      ! gradient that quantity is ||A^T r0||^2 and the question cannot
      ! arise.
      !
      ! rAg is r0 . A g and nAg is ||A g||^2.
      real*8,  intent(in)  :: rAg, nAg
      real*8,  intent(out) :: tau
      logical, intent(out) :: ascends
      tau     = 0.0d0
      ! THE DIRECTION ASCENDS THE MODEL. Shortening a step along it cannot
      ! help -- the predicted reduction is linear in the step length with a
      ! negative slope, so every radius the region tries predicts an
      ! increase and the region is shrunk again. MEASURED on the hot Uranus
      ! element solve: eleven consecutive outer iterations exiting on "the
      ! model promises no reduction", pred falling exactly by 4 at each
      ! factor-4 shrink of the radius, down to a radius of 1.1e-16 (report
      ! B5e section 3.4). The caller drops the leg rather than shortening
      ! it, and the dogleg reduces to the Krylov leg cut to the ball.
      ascends = (nAg .le. 0.0d0 .or. rAg .le. 0.0d0)
      if (ascends) return
      ! THE MINIMIZER ALONG THE DIRECTION ACTUALLY USED.
      ! d/dt [1/2 ||r0 - t A g||^2] = -(r0 . A g) + t ||A g||^2 vanishes at
      ! t = (r0 . A g)/||A g||^2, and that is the minimum of the model along
      ! -g whatever g is. The expression ||g||^2/||A g||^2 minimizes the
      ! model only where g is its exact gradient, in which case
      ! r0 . A g = ||A^T r0||^2 = ||g||^2 and the two agree; g here is the
      ! transpose of the BANDED A, the two numerators are then different
      ! numbers, and the old length does not minimize the model and need not
      ! even decrease it: with r0 = 1, A = 1 and a band of 10 the merit at
      ! the point it returns is 40.5 against 0.5 at the iterate, while the
      ! minimizer reaches 0 (src/tests/krylov_and_dogleg, row
      ! approximate_gradient_length_decreases_the_model; R2).
      !
      ! The point is therefore an APPROXIMATE DESCENT POINT of the model
      ! along a direction the band supplies, and not the Cauchy point of an
      ! exact gradient. The classical trust-region guarantee is a statement
      ! about the latter; naming this one after it is what made the ascent
      ! case look like an anomaly instead of the general situation for an
      ! inexact gradient.
      tau = rAg/nAg
      end subroutine cauchy_length_along_the_banded_gradient

      ! ------------------------------------------------------!

      pure real*8 function predicted_model_decrease(r0, As) result(pred)
      ! THE MODEL'S PROMISED DROP FOR A STEP WHOSE IMAGE IS As, WITHOUT THE
      ! CANCELLATION.
      !   1/2 ( ||r0||^2 - ||r0 + As||^2 ) = -(r0 . As) - 1/2 ||As||^2
      ! algebraically, and the right-hand form is the one to compute: the
      ! left differences two sums of the same size, so on a short step its
      ! value is the rounding of ||r0||^2 and its SIGN is noise -- and the
      ! sign is what decides whether the region shrinks again.
      real*8, dimension(nvar_jac*N), intent(in) :: r0, As
      pred = -sum(r0*As) - 0.5d0*sum(As*As)
      end function predicted_model_decrease

      ! ------------------------------------------------------!

      ! DOGLEG GEOMETRY.
      !
      ! The classical Powell dogleg: the path runs from the origin to the
      ! Cauchy point sU (the minimizer of the linear model along the steepest
      ! descent direction) and from there to the Newton point sN, and the step
      ! is where that path first meets the ball of radius delta. Returned as
      ! the two coefficients of s = aU*sU + aN*sN, because the trust region
      ! needs A*s and gets it from A*sU and A*sN by linearity rather than from
      ! another matrix-vector product.
      !
      ! Three cases, in the order they are tested:
      !   ||sN|| <= delta            -> the Newton point itself
      !   ||sU|| >= delta            -> the Cauchy direction cut to the ball
      !   otherwise                  -> the crossing of the second leg
      subroutine dogleg_step(sU, sN, delta, aU, aN, on_boundary)
      real*8, dimension(nvar_jac*N), intent(in)  :: sU, sN
      real*8,                        intent(in)  :: delta
      real*8,                        intent(out) :: aU, aN
      logical,                       intent(out) :: on_boundary
      real*8, dimension(nvar_jac*N) :: d
      real*8 :: nU, nN, a, b, c, disc, tau
      nN = sqrt(sum(sN*sN))
      nU = sqrt(sum(sU*sU))
      on_boundary = .true.
      if (nN .le. delta) then
         aU = 0.0d0;  aN = 1.0d0;  on_boundary = (nN .ge. delta)
         return
      endif
      if (nU .ge. delta) then
         ! the Cauchy leg already leaves the ball
         if (nU .gt. 0.0d0) then
            aU = delta/nU
         else
            aU = 0.0d0
         endif
         aN = 0.0d0
         return
      endif
      ! second leg: sU + tau (sN - sU), tau in [0,1], ||.|| = delta
      d = sN - sU
      a = sum(d*d)
      b = 2.0d0*sum(sU*d)
      c = sum(sU*sU) - delta*delta
      if (a .le. 0.0d0) then
         aU = 1.0d0;  aN = 0.0d0;  return
      endif
      disc = b*b - 4.0d0*a*c
      if (disc .lt. 0.0d0) disc = 0.0d0
      tau = (-b + sqrt(disc))/(2.0d0*a)
      tau = min(1.0d0, max(0.0d0, tau))
      aU = 1.0d0 - tau
      aN = tau
      end subroutine dogleg_step

      ! ------------------------------------------------------!

      ! THE TWO SHARES A DOGLEG LEG HAS: OF THE STEP, AND OF ITS IMAGE.
      ! The step is s = aU sU + aN sN and, by linearity of A, its image is
      ! A s = aU A sU + aN A sN, so the approximate-gradient leg owns a
      ! share |aU| ||sU|| of the length the two legs contribute and a share
      ! |aU| ||A sU|| of the image they contribute. The shares are taken
      ! against the SUM of the two contributions rather than against the
      ! norm of their sum, so that each lies in [0,1] whatever the angle
      ! between the legs. The model slope of the step is r0 . A s, so it is
      ! the image share that says whether the model describes the step,
      ! while the step share says how much of the step the leg is.
      pure subroutine cauchy_leg_shares_of_step_and_image(sU, sN, AsU,     &
                            AsN, aU, aN, step_share, image_share)
      real*8, dimension(nvar_jac*N), intent(in)  :: sU, sN, AsU, AsN
      real*8,                        intent(in)  :: aU, aN
      real*8,                        intent(out) :: step_share, image_share
      real*8 :: lU, lN, iU, iN
      lU = abs(aU)*sqrt(sum(sU*sU))
      lN = abs(aN)*sqrt(sum(sN*sN))
      iU = abs(aU)*sqrt(sum(AsU*AsU))
      iN = abs(aN)*sqrt(sum(AsN*AsN))
      step_share  = lU/max(lU + lN, 1.0d-300)
      image_share = iU/max(iU + iN, 1.0d-300)
      end subroutine cauchy_leg_shares_of_step_and_image

      ! ------------------------------------------------------!

      ! IS THE LEG THE IMAGE OF A STEP IT IS NOT PART OF? A leg above
      ! cauchy_leg_image_share_max of the image and below
      ! cauchy_leg_step_share_max of the length carries the model slope of a
      ! direction the step barely contains, and the model is then not a
      ! model of the step. Both conditions are required: a leg that is most
      ! of the step is allowed to be most of its image, and a leg that is
      ! little of the image may be as short as it likes.
      pure logical function cauchy_leg_is_image_without_step(step_share,   &
                                             image_share) result(drop)
      real*8, intent(in) :: step_share, image_share
      drop = (image_share .gt. cauchy_leg_image_share_max) .and.          &
             (step_share  .lt. cauchy_leg_step_share_max)
      end function cauchy_leg_is_image_without_step

      ! ------------------------------------------------------!

      ! IS THE LEG SHORTER THAN THE RADIUS RESOLVES? A leg below
      ! cauchy_leg_min_fraction of the radius moves the step by less than
      ! the radius's own control does, so it is not part of the step. The
      ! rule is a proxy for the image test above and holds only where the
      ! step stands on the ball; where the caller admits the leg by its
      ! image alone the length is not asked about at all
      ! (cauchy_leg_admitted_by_image_alone).
      pure logical function cauchy_leg_is_dropped_on_its_length(          &
                            leg_length, delta, by_image_alone) result(drop)
      real*8,  intent(in) :: leg_length, delta
      logical, intent(in) :: by_image_alone
      drop = (.not. by_image_alone) .and.                                 &
             (leg_length .lt. cauchy_leg_min_fraction*delta)
      end function cauchy_leg_is_dropped_on_its_length

      ! ------------------------------------------------------!

      subroutine ray_slope_verdict(mslope, fslope, f2p, f2m, ok_probes,   &
                                   ratio, merit_reproducibility,          &
                                   ray_floor, ray_verdict, model_ok,      &
                                   refuse)
      ! DOES THE MODEL'S SLOPE DESCRIBE THE FUNCTION, and is the finite
      ! difference that asks the question large enough to answer it.
      !
      ! mslope is the model's directional derivative along the step, fslope
      ! the difference quotient of the true merit along the same ray, f2p
      ! and f2m the two probe merits the quotient was built from, ok_probes
      ! whether both probes are states the equilibrium sweep could describe,
      ! ratio the actual reduction over the predicted one at the FULL step,
      ! and merit_reproducibility the measured spread of the merit over two
      ! evaluations of one state.
      !
      ! ray_verdict is what the comparison found, one of the tr_ray_* codes
      ! or tr_step_taken when the two slopes agree. refuse and model_ok are
      ! what the step control acts on, and they are not the same statement:
      ! a verdict against the model is not acted on when the ratio itself is
      ! sound (ratio_sound_band), so the verdict remains available to be
      ! reported while the step goes on to the acceptance test; and where
      ! the ratio is at or below eta_accept the refusal is NAMED by the
      ! ratio, because that is the informative statement about a trial the
      ! acceptance test refuses on its own.
      real*8,  intent(in)  :: mslope, fslope, f2p, f2m, ratio
      real*8,  intent(in)  :: merit_reproducibility
      logical, intent(in)  :: ok_probes
      real*8,  intent(out) :: ray_floor
      integer, intent(out) :: ray_verdict, refuse
      logical, intent(out) :: model_ok
      logical :: slopes_agree
      ! THE FLOOR OF THE DIFFERENCE ITSELF. f2p and f2m are two evaluations
      ! of a merit of size max(f2p,f2m); their difference carries
      ! information only while it stands above the noise of the two numbers
      ! it is built from, and that noise is not machine epsilon here: the
      ! composition is eliminated iteratively at each probe, so the merit of
      ! a state is a function of Y only to that accuracy. Below the floor
      ! the quotient is noise divided by a small number and says nothing
      ! about the model, which is a different verdict from "the model
      ! disagrees".
      !
      ! THE FIRST PART is outer_residual_closure_target times the size of
      ! the two numbers: what the outer residual is allowed to keep of the
      ! inner closure's error. It is a bound only because the sweep is asked
      ! for that target divided by the MEASURED amplification of the base
      ! boundary (closure_amplification_at_the_base); asking the sweep for
      ! the target itself and using the target here, as this line did,
      ! assumes an amplification of one where the base has none. THE SECOND
      ! PART is what two evaluations of one state actually differ by, with
      ! nothing assumed about what the elimination reached or about what one
      ! evaluation leaves behind for the next.
      !
      ! WHICH OF THE TWO GOVERNS IS A MEASUREMENT, and on the atomic element
      ! reload it is the assumed one (N7b, 20 outer iterations): two
      ! evaluations of one state differ by 1.78e-15 while the merit is 14.95
      ! and by 0 to 4.3e-19 while it is 2.7e-5, which is one to a few tens
      ! of units in the last place of the merit itself, four to seven
      ! decades below the assumed part. So on that operator the composition
      ! elimination reaches far more than it is asked for, the assumed floor
      ! is conservative, and a difference of 3e-12 between two probes is a
      ! real difference between two states and not noise. The measured part
      ! is kept because a floor that is asserted is not a floor: on an
      ! operator whose elimination does not close, it is the one that binds.
      ray_floor = max(outer_residual_closure_target                       &
                         *max(abs(f2p), abs(f2m)),                        &
                      ray_floor_factor*merit_reproducibility)
      slopes_agree = (mslope*fslope .gt. 0.0d0) .and.                     &
                     (abs(fslope) .le. ray_factor*abs(mslope)) .and.      &
                     (abs(mslope) .le. ray_factor*abs(fslope))
      ! A probe the equilibrium solve could not describe carries no merit
      ! value, so the difference quotient built from it is not the slope of
      ! the function along the ray and cannot confirm the model. An
      ! inadmissible probe therefore means NO verdict, which is the same
      ! answer the test gives when the two slopes disagree: shrink the
      ! region and rebuild the model.
      if (.not. ok_probes) then
         ray_verdict = tr_ray_probe_inadmissible
      else if (slopes_agree) then
         ray_verdict = tr_step_taken
      else if (abs(f2p - f2m) .le. ray_floor) then
         ray_verdict = tr_ray_below_the_noise_floor
      else if (mslope*fslope .le. 0.0d0) then
         ray_verdict = tr_ray_slopes_disagree_sign
      else
         ray_verdict = tr_ray_slopes_disagree_size
      endif
      model_ok = (ray_verdict .eq. tr_step_taken)
      refuse   = ray_verdict
      ! A SOUND RATIO IS NOT OVERRULED BY A SLOPE QUOTIENT. The ray test
      ! exists to catch a model that promises a reduction the function does
      ! not deliver, and the reduction ratio measures exactly that, at the
      ! full step, with no finite difference in it. Where the ratio stands
      ! within ratio_sound_band of unity the slope comparison adds no
      ! information and can only remove steps the model described correctly.
      ! The verdict is still returned, so the diagnostic is not lost, and a
      ! model that really over-predicts is still refused: by this test
      ! outside the band, and by the acceptance threshold below it.
      !
      ! An inadmissible probe is not covered: there the test could not be
      ! taken at all, and a neighborhood of the step that holds states the
      ! equilibrium sweep cannot describe is a statement about the step and
      ! not about the noise of a quotient.
      if (.not. model_ok .and.                                            &
          ray_verdict .ne. tr_ray_probe_inadmissible .and.                &
          abs(ratio - 1.0d0) .le. ratio_sound_band) then
         model_ok = .true.
         refuse   = tr_step_taken
      endif
      ! A TRIAL THE ACCEPTANCE TEST ALREADY REFUSES IS NAMED BY THE RATIO.
      ! The three verdicts above about the quotient -- the two slope
      ! disagreements and the noise floor -- are statements about a finite
      ! difference taken over a thousandth of the step; the ratio is a
      ! statement about the whole of it, and where it stands at or below
      ! eta_accept the step is refused whatever the quotient says. The
      ! ledger then records the informative name: MEASURED on the atomic
      ! element reload, iterations at ratio -533 and -1493 (actual reduction
      ! NEGATIVE) exited as "the model and the true slope differ in sign"
      ! (N7b, noticed item 1). model_ok stays false, so the region is cut on
      ! the same rule as before and no step changes; only the name does.
      !
      ! An inadmissible probe is not covered: there the neighborhood of the
      ! step holds states the equilibrium sweep cannot describe, which is
      ! something the ratio does not say.
      if (.not. model_ok .and.                                            &
          ray_verdict .ne. tr_ray_probe_inadmissible .and.                &
          ratio .le. eta_accept) refuse = tr_ratio_below_acceptance
      end subroutine ray_slope_verdict

      ! ------------------------------------------------------!

      ! A SCALED TRUST-REGION STEP FOR THE FOUR-UNKNOWN COUPLED SOLVE.
      !
      ! WHY. The pseudo-transient line search controls the step by a
      ! pseudo-time and then asks the merit whether the result is acceptable.
      ! Measured (docs/Update_EXHALE_stage1.md section 142), on the coupled hot
      ! Uranus that control fails in a specific way: the step it builds is an
      ! ASCENT direction of the merit in the model's own arithmetic from outer
      ! iteration 5 onward, because the I/dtau term's contribution to the
      ! merit slope carries no fixed sign under the two-sided scaling; and
      ! raising dtau until the slope turns negative does not help either,
      ! because the arc over which the raised direction descends is far
      ! shorter than the step. Both failures are failures of STEP LENGTH
      ! control, which is what a trust region is.
      !
      ! WHAT. The model is the same scaled linear model the Krylov cycle is
      ! built from,
      !
      !     m(s) = 1/2 || r0 + A s ||^2 ,   r0 = F/Drow ,   A = Drow^-1 J D ,
      !
      ! minimized over ||s|| <= delta by a dogleg between the Cauchy point and
      ! the Krylov point. A*sN comes from the Arnoldi relation for free; A*sU
      ! costs one matrix-free product; every other A*s on the dogleg path is a
      ! linear combination of those two, so the predicted reduction is exact
      ! for the model and needs no further evaluation.
      !
      ! The gradient A^T r0 is taken from the BANDED A, because a transpose
      ! is not available matrix-free. That makes the Cauchy direction an
      ! approximation while the Newton direction and every predicted reduction
      ! are not; the ray test below is what stops a bad Cauchy direction from
      ! being accepted on a model it does not describe.
      !
      ! ADMISSIBILITY IS PART OF THE STEP, not a test after it: the step is
      ! scaled back by theta until the trial state has positive rho and p, a
      ! non-negative carrier density, a carrier density inside the hydrogen
      ! its cell has, and a chemistry the equilibrium sweep can certify --
      ! which is exactly what eval_residual's "admissible" reports. The
      ! predicted reduction is recomputed for the scaled-back step, so the
      ! ratio compares like with like.
      !
      ! THE MODEL IS REJECTED, not merely disbelieved, when its directional
      ! derivative disagrees with a central finite difference of the true
      ! merit along the same ray. That is the check the ray scans of section
      ! 142 showed is needed: at the stagnating iterates the merit turns over
      ! inside a step, and a model that does not see the turn cannot be used
      ! to size one.
      subroutine trust_region_step(Y, F, F_jac, f_sp, D, Drow, ab, abf, ipiv, &
                                   idtau, gm_m, gm_rtol, f2, delta,        &
                                   dY, f2_new, pred, actual, ratio, snorm, &
                                   n_theta_cuts, gmit, gm_outcome,         &
                                   gm_resid_rel,                           &
                                   model_ok, accepted, Yacc, Facc, f_sp_acc, &
                                   f2_ref_out, refuse, rows_acc,           &
                                   cauchy_ascends)
      real*8, dimension(nvar_jac*N),          intent(in)    :: Y, F, D, Drow
      ! The base point of every finite difference of the Newton model: the
      ! residual at Y under n_eq_sweeps_model passes of the composition
      ! elimination. F is the residual at
      ! the state's own composition and is what the step must reduce.
      real*8, dimension(nvar_jac*N),          intent(in)    :: F_jac
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in)    :: ab
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(inout) :: abf
      integer, dimension(nvar_jac*N),         intent(inout) :: ipiv
      real*8,  intent(in)    :: idtau, gm_rtol, f2
      integer, intent(in)    :: gm_m
      real*8,  intent(inout) :: delta
      real*8, dimension(nvar_jac*N), intent(out) :: dY
      ! The accepted trial, handed back whole: the caller adopts Y, F and the
      ! species array together, and the LAST residual evaluation this routine
      ! makes is of exactly this state, so the frozen background the caller
      ! keeps describes the state it adopts.
      real*8, dimension(nvar_jac*N),          intent(out) :: Yacc, Facc
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp_acc
      real*8,  intent(out) :: f2_new, pred, actual, ratio, snorm
      ! The merit of the iterate measured in the TRIAL's mode. Reported beside
      ! the caller's own f2 so that the gap between the two evaluation modes
      ! is visible rather than hidden inside the ratio.
      real*8,  intent(out) :: f2_ref_out
      integer, intent(out) :: n_theta_cuts, gmit
      ! The relative residual the Krylov cycle of the Newton leg reached,
      ! against the gm_rtol it was asked for (pgmres).
      real*8,  intent(out) :: gm_resid_rel
      ! WHICH OUTCOME THE KRYLOV CYCLE OF THE NEWTON LEG HAD, one of the
      ! gm_* names (pgmres).
      integer, intent(out) :: gm_outcome
      logical, intent(out) :: model_ok, accepted
      ! WHICH EXIT THIS ITERATION TOOK, one of the tr_* codes above. It is
      ! not derivable from model_ok and accepted: those two are false
      ! together for seven different reasons, and the repair each of them
      ! asks for is a different one.
      integer, intent(out) :: refuse
      ! The judged rows of the ACCEPTED trial (certified_row_measures), so
      ! that the caller's best-iterate ledger ranks the state this routine
      ! hands it on the measure the state will be judged by. Zero when no
      ! step was taken.
      real*8, dimension(3), intent(out) :: rows_acc
      ! Whether the banded gradient ascended the model at this iterate, so
      ! that the Cauchy leg was dropped and the dogleg reduced to the Krylov
      ! leg cut to the ball. Counted by the caller: it is a property of the
      ! band as a model of the transpose, not of the state.
      logical, intent(out) :: cauchy_ascends

      real*8, dimension(nvar_jac*N) :: r0, g, Ag, sN, AsN, sU, AsU, s, As
      real*8, dimension(nvar_jac*N) :: Ytry, Ftry, Jv, Ftrial, Ystep
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng)   :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng) :: utry, Wtry
      real*8  :: aU, aN, nAg, tau_c, theta, f2t, mslope, fslope
      ! The approximate-gradient leg as the minimizer along -g makes it,
      ! before the dogleg decides whether it is long enough to be part of
      ! the step (cauchy_leg_min_fraction).
      real*8  :: nU_full
      logical :: cauchy_leg_short
      real*8  :: cauchy_step_share, cauchy_image_share
      logical :: cauchy_leg_all_image
      real*8  :: f2p, f2m, hray, gnorm, f2_ref
      integer :: it, saved_mode
      logical :: okres, on_boundary, ok_state
      ! Both ray probes describable states (see the ray test below).
      logical :: ok_ray
      ! The radius floor exists so that a solve that cannot move says so
      ! instead of grinding the radius to zero. The acceptance threshold
      ! eta_accept is at module scope, beside the constants of the ray test
      ! that reads it.
      real*8,  parameter :: delta_floor = 1.0d-10
      integer, parameter :: n_theta_max = 12
      ! Ray test: the fraction of the step the finite difference is taken
      ! over. 1e-3 of the step is inside the arc the section 142 ray scans
      ! measured as model-valid. How far the two slopes may then differ, and
      ! what their difference must stand above, are ray_factor,
      ! ray_floor_factor and ratio_sound_band at module scope.
      real*8,  parameter :: ray_h      = 1.0d-3
      character(len=8) :: tr_trace_env
      logical :: tr_trace
      real*8  :: idtau_leg
      real*8  :: delta_in, ray_floor, worst_row, rAg, gm_snorm_leg
      ! WHAT THE RAY TEST FOUND, whether or not it was acted on: the two
      ! probes and the slopes are already above, this is the verdict
      ! ray_slope_verdict returned, and the reproducibility of the merit
      ! that set the floor the verdict used.
      integer :: ray_verdict
      real*8  :: merit_reproducibility, merit_first
      ! Elimination passes the first of the two reference evaluations took,
      ! beside the second's, so that a reproducibility set by a pass the
      ! elimination did or did not take can be read as such.
      integer :: n_passes_first
      integer :: i_worst, n_active
      logical :: krylov_leg_missing
      ! Species unknowns of the trial, and of the minus ray, that sat
      ! outside their own bounds before the projection.
      integer :: n_out, n_ray_out
      ! The restoration phase's own quantities at the trial that was
      ! evaluated: cells and rows it moved, and the worst relative violation
      ! of the shared element rows before and after it
      ! (restore_the_element_budget).
      integer :: n_rest_cell, n_rest_row, n_rest_out
      real*8  :: rest_before, rest_after
      character(len=48) :: what_worst
      ! Long enough for the longest sentence tr_exit_text and gm_outcome_text
      ! can return (67 characters), so an exit is never named in part.
      character(len=80) :: why_txt
      ! The element-row diagnostic's own quantities (elem_diag_on): the
      ! model gradient BEFORE the held set was zeroed out of it, so that the
      ! unknowns the step holds can be ranked by how hard the model pushed
      ! them out of their box; the refusal statistics at entry, so that the
      ! samples this step was refused can be counted; and the dogleg at the
      ! entry radius, taken before any trial.
      real*8, allocatable :: g_unheld(:), As_true(:)
      ! THE RADIUS PAIR THE MODEL IS READ AT, and the leg pair it is read
      ! with: one step, image and trial state for a radius and for a
      ! quarter of it, with the approximate-gradient leg admitted and with
      ! it dropped, so that the promise's dependence on the radius can be
      ! read off directly ([diag 9]).
      real*8, allocatable :: dl_sU(:), dl_AsU(:), dl_s(:), dl_As(:)
      real*8  :: dl_delta, dl_aU, dl_aN, dl_pred, dl_actual, dl_f2
      ! The FULL leg's two shares, for the reading: the decision above uses
      ! them only where the length test kept the leg, but the reading wants
      ! them at every trial, so that the window in which the length test
      ! keeps a leg the image test removes can be counted.
      real*8  :: dl_step_share, dl_image_share
      integer :: k_radius, k_leg, n_dl_out
      logical :: dl_bnd, dl_ok
      real*8  :: aU_diag, aN_diag, gm_resid_true, arnoldi_gap
      real*8  :: pred_cauchy, pred_krylov, pred_dogleg, pred_fresh
      integer :: n_refused_at_entry, i_top(3), n_top, ihit
      logical :: on_boundary_diag, ok_true
      type(solve_refusal_statistics) :: statistics_at_entry

      call get_environment_variable('EXHALE_TR_TRACE', tr_trace_env)
      tr_trace = (trim(tr_trace_env) .eq. '1')
      ! The dogleg's Newton leg is solved WITHOUT the pseudo-transient shift
      ! (A sN = -r0, not (idtau I + A) sN = -r0), because the dogleg's
      ! monotone-model premise requires the leg to be the minimizer of the
      ! SAME model the predicted reduction is measured in.  With the shifted
      ! leg the model residual at full length is exactly the shift,
      ! idtau ||(D/Drow) sN||, which sat above ||r0|| and made the model
      ! predict a rise on half of all iterations, so the radius walked
      ! halve-halve-double downward and could never carry a long step
      ! (verified by trace and by a Krylov-tolerance ladder that changed
      ! nothing, section 159.3).  The banded preconditioner keeps its shift:
      ! a preconditioner need not be exact.  EXHALE_TR_SHIFTED_LEG=1 restores
      ! the shifted leg for measurement.
      call get_environment_variable('EXHALE_TR_SHIFTED_LEG', tr_trace_env)
      idtau_leg = 0.0d0
      if (trim(tr_trace_env) .eq. '1') idtau_leg = idtau
      n_rest_cell = 0;  n_rest_row = 0;  n_rest_out = 0
      rest_before = 0.0d0;  rest_after = 0.0d0
      dY = 0.0d0;  f2_new = f2;  accepted = .false.;  model_ok = .true.
      f2_ref = f2
      Yacc = Y;  Facc = F;  f_sp_acc = f_sp
      pred = 0.0d0;  actual = 0.0d0;  ratio = 0.0d0;  snorm = 0.0d0
      n_theta_cuts = 0;  gmit = 0;  gm_outcome = gm_tolerance_reached
      gm_resid_rel = 0.0d0;  gm_snorm_leg = 0.0d0
      saved_mode = weno_mode
      refuse = tr_step_taken;  delta_in = delta
      rows_acc = 0.0d0;  cauchy_ascends = .false.;  rAg = 0.0d0
      nAg = 0.0d0;  tau_c = 0.0d0;  gnorm = 0.0d0
      ! The image of the approximate-gradient direction, zero until the
      ! probe forms it: the dogleg rebuilds the Cauchy leg as -tau_c*Ag at
      ! every trial, and with no probe tau_c is zero and Ag has to be a
      ! number for the product to be one.
      Ag = 0.0d0
      krylov_leg_missing = .false.
      mslope = 0.0d0;  fslope = 0.0d0;  ray_floor = 0.0d0
      f2p = 0.0d0;  f2m = 0.0d0
      ray_verdict = tr_step_taken;  merit_reproducibility = 0.0d0
      merit_first = 0.0d0;  n_passes_first = 0
      Ftrial = F
      n_refused_at_entry = sum(n_eval_refusal(1:6))

      ! THE REFERENCE MERIT IS MEASURED IN THE SAME MODE AS THE TRIALS.
      ! The caller's f2 is the merit of this iterate with the WENO3 weights
      ! and the limited slope FROZEN at it (mode 1); every trial is measured
      ! with them recomputed (mode 0), because that is the residual the solve
      ! is driving to zero (section 126). Those are different functions, and a
      ! reduction ratio taken across them is not a ratio of anything: measured
      ! on the coupled hot Uranus the two differ by enough that a step which
      ! lowers the mode-0 merit is followed by a mode-1 re-evaluation that
      ! raises it, so the region grows on a reduction the next iteration does
      ! not see. One extra evaluation at the iterate fixes it, and both the
      ! predicted and the actual reduction are then taken about the same
      ! baseline. The MODEL is still the Krylov one; only its constant term is
      ! the residual actually being reduced.
      ! AND HOW REPRODUCIBLE THAT MERIT IS, measured here and used by the
      ! ray test below as the floor its difference quotient has to stand
      ! above. The merit of a state is not a function of Y alone: the
      ! composition is an eliminated variable solved iteratively at every
      ! evaluation, and what one evaluation leaves in the modules it calls
      ! is what the next one starts from. TWO EVALUATIONS OF THIS STATE,
      ! from the same seed f_sp, in the same reconstruction mode and under
      ! the same elimination control as every probe below, therefore differ
      ! by exactly the amount two probes of two nearby states differ by for
      ! reasons that are not the states. The second is the one the model is
      ! built about, so that the background the modules hold is the one this
      ! iterate's own last evaluation left, as everywhere else in this
      ! routine.
      !
      ! COST: one residual evaluation per outer iteration, which is the
      ! price of one Jacobian product. It buys a floor that is measured
      ! rather than asserted; whether the measured part or the assumed one
      ! then binds is stated at ray_slope_verdict, with what each of them
      ! came to on the atomic element reload.
      weno_mode = 0
      call eval_residual(Y, f_sp, fwork, Ftry, heat, cool)
      merit_first    = 0.5d0*sum((Ftry/Drow)**2)
      n_passes_first = n_eq_sweeps_last
      ! THE TERMS OF THE ELEMENT ROWS AT THE OUTERMOST CELLS, printed by the
      ! operator that forms them, on THIS evaluation of the iterate
      ! (trace_row_terms_diag). The row of the last cell balances an
      ! advective divergence against a diffusive flux whose outer face
      ! carries none, and whether it is a balance a step can restore is read
      ! from the sizes of those terms. Armed for one evaluation and
      ! disarmed again, so no probe state prints.
      if (elem_diag_here) then
         trace_row_terms_diag   = .true.
         element_row_terms_diag = .true.
         trace_row_terms_from   = max(1, N - 5)
         if (elem_diag_cell_from .gt. 0)                                  &
            trace_row_terms_from = min(max(2, elem_diag_cell_from), N)
      endif
      call eval_residual(Y, f_sp, fwork, Ftry, heat, cool)
      trace_row_terms_diag   = .false.
      element_row_terms_diag = .false.
      r0 = Ftry/Drow
      f2_ref = sqrt(sum(r0*r0))
      merit_reproducibility = abs(0.5d0*f2_ref*f2_ref - merit_first)
      write(*,'(A,ES11.3,A,I0,A,I0,A)')                                   &
           ' (JFNK) [TR]      merit reproducibility at the iterate',      &
           merit_reproducibility, ' (two evaluations, ',                   &
           n_passes_first, ' and ', n_eq_sweeps_last,                     &
           ' elimination passes)'

      ! --- the bounds this step holds, before any direction is sampled ---
      ! The model gradient decides which of them are active, so it is formed
      ! here rather than with the Cauchy leg below: every sample of the
      ! operator, the Krylov cycle's included, has to see the same held set
      ! or the model and the step are two different maps
      ! (fix_active_species_bounds).
      call band_matvec_transpose(ab, r0, g)
      if (elem_diag_here) then
         allocate(g_unheld(nvar_jac*N))
         g_unheld = g
      endif
      call fix_active_species_bounds(Y, g, n_active)
      n_active_bounds_max = max(n_active_bounds_max, n_active)
      ! WHICH UNKNOWNS THIS STEP MAY NOT MOVE, ranked by the size of the
      ! model gradient the hold removed: that is how hard the model was
      ! pushing each of them out of its own box. A held component is a
      ! direction the step does not have; a refused SAMPLE is a direction
      ! the step could not measure, and the two are counted separately
      ! because they call for different repairs.
      if (elem_diag_here) then
         i_top = 0;  n_top = 0
         if (allocated(species_bound_is_active)) then
            do ihit = 1, 3
               i_worst = 0;  worst_row = -1.0d0
               do it = 1, nvar_jac*N
                  if (.not. species_bound_is_active(it)) cycle
                  if (it .eq. i_top(1) .or. it .eq. i_top(2)) cycle
                  if (abs(g_unheld(it)) .gt. worst_row) then
                     worst_row = abs(g_unheld(it));  i_worst = it
                  endif
               enddo
               if (i_worst .eq. 0) exit
               i_top(ihit) = i_worst;  n_top = ihit
            enddo
         endif
         write(*,'(A,I0,A,I0)') ' (JFNK) [diag 5] unknowns held at a'//   &
              ' bound ', n_active, ', residual samples refused so far in'//&
              ' this solve ', n_refused_at_entry
         do ihit = 1, n_top
            call name_of_unknown(i_top(ihit), what_worst)
            write(*,'(A,A,A,ES11.3)') ' (JFNK) [diag 5]   the ',          &
                 trim(what_worst), ': model gradient removed by the hold ',&
                 g_unheld(i_top(ihit))
         enddo
      endif

      ! --- THE BAND ROW OF THE ROW THE MERIT IS LARGEST ON ---
      ! Which unknowns the model believes that row responds to, and how
      ! large each derivative is beside the scale the row is judged on. A
      ! row whose only entries are its own cell and its two neighbors is a
      ! transport balance; an entry pattern that stops at the last cell is
      ! what an outer boundary with no equation beyond it looks like from
      ! inside the model. The Jacobian is the BANDED one the preconditioner
      ! is built from (band_matvec_transpose states the storage), so this is
      ! the model's own derivative and not a re-differentiation.
      if (elem_diag_here) then
         i_worst = 1;  worst_row = -1.0d0
         do it = 1, nvar_jac*N
            if (abs(r0(it)) .gt. worst_row) then
               worst_row = abs(r0(it));  i_worst = it
            endif
         enddo
         call name_of_unknown(i_worst, what_worst)
         write(*,'(A,A,A,ES11.3,A,ES11.3)') ' (JFNK) [diag 10] the row'// &
              ' the merit is largest on: ', trim(what_worst),             &
              ', scaled residual ', worst_row, ', row scale ',            &
              Drow(i_worst)
         do it = max(1, i_worst - kl_jac), min(nvar_jac*N, i_worst + ku_jac)
            if (ab(kl_jac+ku_jac+1 + i_worst - it, it) .eq. 0.0d0) cycle
            call name_of_unknown(it, what_worst)
            write(*,'(A,A,A,ES13.5,A,ES13.5)') ' (JFNK) [diag 10]   d/d ',&
                 trim(what_worst), ': ',                                  &
                 ab(kl_jac+ku_jac+1 + i_worst - it, it), ', scaled by the'&
                 //' row ', ab(kl_jac+ku_jac+1 + i_worst - it, it)        &
                 /max(Drow(i_worst), 1.0d-300)
         enddo
      endif

      ! --- the Krylov (Newton) point, and A applied to it, for free ---
      ! THE RADIUS IS HANDED TO THE CYCLE, and read by it only where the
      ! truncation is armed: the cycle then stops on the boundary of THIS
      ! ball and the dogleg below has a leg it can take whole, instead of
      ! a direction four decades longer than the region
      ! (krylov_truncated_on_the_trust_ball). The radius is the one this
      ! iteration entered with; the trial loop shrinks it afterwards and
      ! the dogleg cuts the leg as it always has.
      call pgmres_with_restarts(Y, F_jac, f_sp, D, Drow, abf, ipiv,        &
                  idtau_leg, -r0, sN,                                     &
                  gm_m, gm_rtol, gmit, gm_outcome, gm_resid_rel,          &
                  gm_snorm_leg, AsN, delta)
      ! A KRYLOV LEG IS NOT WHAT A TRUST REGION GUARANTEES; THE CAUCHY POINT
      ! IS. Where the Krylov cycle could not sample its first direction the
      ! step used to be abandoned, and on a bound-constrained system that is
      ! a whole iteration lost to the geometry of the feasible set: MEASURED
      ! on the coupled `mol_carrier` reload, 12 of 20 outer iterations ended
      ! on "no Krylov direction could be sampled" with the merit standing at
      ! 7.93 and the carrier of cell 215 pinned at zero eleven decades below
      ! its own row's root of 2.75e-06.
      !
      ! The Cauchy direction is available in exactly that case, and not by
      ! luck: the gradient it is built from has been zeroed on every unknown
      ! held at a bound (fix_active_species_bounds), so it has no component
      ! pointing out of the box anywhere and a step along it always has room
      ! inside. So the dogleg degenerates to the Cauchy point cut to the
      ! ball -- which is the classical guarantee of a trust-region method,
      ! the step every convergence proof of one rests on (Nocedal and
      ! Wright, Numerical Optimization, 2nd ed., section 4.1) -- and the
      ! iteration takes a step instead of shrinking its radius against
      ! nothing.
      ! The leg is missing exactly where the cycle managed no product: the
      ! feasible set left no room along the first preconditioned direction
      ! (gm_no_direction_sampled) or the preconditioner failed on it
      ! (gm_preconditioner_failed).
      krylov_leg_missing = (gmit .le. 0)
      if (krylov_leg_missing) then
         sN = 0.0d0;  AsN = 0.0d0
      endif
      ! A_z = idtau*(D/Drow) + A, so remove the pseudo-transient part to get
      ! the Jacobian action the model is written in.
      AsN = AsN - idtau_leg*(D/Drow)*sN
      ! The banded preconditioner has no knowledge of the held set, so the
      ! leg it returns can carry components there; they are dropped, and
      ! A sN is unaffected because the Arnoldi relation was built from
      ! directions that were held (jv_product).
      call hold_the_active_bounds_of(sN)
      call project_out_of_the_active_element_constraints(sN)
      ! The model residual at the FULL Newton leg is, for an exact leg,
      ! r0 + A sN = -idtau (D/Drow) sN: the pseudo-transient shift, not zero.
      ! Printed against ||r0|| because their ratio decides the sign of the
      ! full-length predicted reduction.
      if (tr_trace) write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')               &
           ' (JFNK) [TR-trace] ||r0||=', sqrt(sum(r0*r0)),                &
           '  idtau*||(D/Drow)sN||=', idtau*sqrt(sum(((D/Drow)*sN)**2)),  &
           '  ||r0+AsN||=', sqrt(sum((r0+AsN)*(r0+AsN)))

      ! WHAT THE KRYLOV LEG IS, MEASURED RATHER THAN REPORTED. The relative
      ! residual pgmres hands back is the reduced least-squares problem's,
      ! and the image AsN it hands back comes from the Arnoldi relation; a
      ! subspace that lost orthogonality, or a reduced problem that lost
      ! rank, makes both of them statements about the recursion and not
      ! about the operator. One product of the operator on the leg the
      ! dogleg will actually use settles it, and the gap between that image
      ! and the Arnoldi one is the fidelity of the model the predicted
      ! reduction is taken in.
      !
      ! The product's refusals belong to the diagnostic and not to the
      ! solve, so the statistics are put back around it
      ! (solve_refusal_statistics).
      if (elem_diag_here) then
         gm_resid_true = -1.0d0;  arnoldi_gap = -1.0d0
         if (sqrt(sum(sN*sN)) .gt. 0.0d0) then
            allocate(As_true(nvar_jac*N))
            call hold_solve_refusal_statistics(statistics_at_entry)
            call jv_product(Y, F_jac, f_sp, D*sN, Jv, ok_true)
            call put_back_solve_refusal_statistics(statistics_at_entry)
            if (ok_true) then
               As_true = Jv/Drow
               gm_resid_true = sqrt(sum((r0 + As_true)**2))              &
                             / max(sqrt(sum(r0*r0)), 1.0d-300)
               arnoldi_gap = sqrt(sum((AsN - As_true)**2))               &
                           / max(sqrt(sum(As_true*As_true)), 1.0d-300)
            endif
            deallocate(As_true)
         endif
         call gm_outcome_text(gm_outcome, why_txt)
         write(*,'(A,I0,A,I0,A,ES11.3,A,ES9.2)')                          &
              ' (JFNK) [diag 4] Krylov leg: products ', gmit, ' of ',     &
              gm_m, ', relative residual returned ', gm_resid_rel,        &
              ' against ', gm_rtol
         write(*,'(A,A)') ' (JFNK) [diag 4]   outcome: ', trim(why_txt)
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                          &
              ' (JFNK) [diag 4]   TRUE relative residual of the returned'//&
              ' step ', gm_resid_true, ', ||sN|| ', sqrt(sum(sN*sN)),     &
              ', relative gap of the Arnoldi image ', arnoldi_gap
         ! The loss of orthogonality of THIS cycle's basis, where it was
         ! measured (EXHALE_GM_ORTHO=1), beside the gap it may explain.
         if (gm_measure_orthogonality)                                    &
            write(*,'(A,ES11.3,A,L1)') ' (JFNK) [diag 4]   loss of'//     &
                 ' orthogonality of this basis ',                         &
                 gm_orthogonality_loss_cycle,                             &
                 ', orthogonalized twice ', gm_reorthogonalize
      endif

      ! --- WHAT HOLDS THE CYCLE ON THIS SYSTEM ---
      ! Three measurements of the SAME linear system the leg above was
      ! taken from, at this iterate, adopting nothing: the cycle over a
      ! ladder of subspace sizes, the Ritz values of the preconditioned
      ! operator, and the difference between the full Jacobian action and
      ! the banded one on the directions the spectrum names. Each is off by
      ! default and each costs products of the operator, so all three speak
      ! only where the element diagnostic does.
      if (elem_diag_here .and. krylov_size_scan_on)                      &
         call krylov_cycle_over_subspace_sizes(Y, F_jac, f_sp, D, Drow,  &
                   abf, ipiv, idtau_leg, -r0, gm_rtol)
      if (elem_diag_here .and. precond_spectrum_on)                      &
         call ritz_values_of_the_preconditioned_operator(Y, F_jac, f_sp, &
                   D, Drow, abf, ipiv, idtau_leg)
      if (elem_diag_here .and. band_difference_on)                       &
         call what_the_band_omits_of_the_jacobian(Y, F_jac, f_sp, D,     &
                   Drow, ab, abf, ipiv, idtau_leg, r0)
      ! And the row that binds read term by term, with the band's entries
      ! for it attributed to those terms (front_row_on).
      if (elem_diag_here .and. front_row_on)                             &
         call the_binding_species_row(Y, f_sp, D, Drow, ab, r0)

      ! --- the Cauchy point, from the gradient formed above ---
      gnorm = sqrt(sum(g*g))
      tau_c = 0.0d0
      if (gnorm .le. 0.0d0) then
         sU = 0.0d0;  AsU = 0.0d0
      else
         call jv_product(Y, F_jac, f_sp, D*g, Jv, ok_state)
         if (.not. ok_state) then
            model_ok = .false.
            refuse = tr_cauchy_probe_inadmissible
            f2_ref_out = f2_ref
            weno_mode = saved_mode
            call release_active_species_bounds
            return
         endif
         Ag  = Jv/Drow
         nAg = sum(Ag*Ag)
         ! THE SLOPE OF THE MODEL ALONG THE DIRECTION ACTUALLY USED.
         ! d/dt [ 1/2 ||r0 - t A g||^2 ] at t = 0 is -(r0 . A g), so the
         ! model falls along -g only while r0 . A g is positive. With the
         ! EXACT gradient that quantity is ||A^T r0||^2 and the question
         ! cannot arise; this g is the transpose of the BANDED A, because a
         ! transpose is not available matrix-free, and nothing makes it a
         ! descent direction of the model it is used in.
         rAg = sum(r0*Ag)
         ! The length along -g, and whether that direction descends the
         ! model at all (cauchy_length_along_the_banded_gradient). Where it
         ! ascends the leg is dropped and the dogleg reduces to the Krylov
         ! leg cut to the ball, which does predict a reduction: GMRES
         ! returns ||r0 + A sN|| < ||r0||, hence r0 . A sN < -||A sN||^2/2,
         ! and the model change along t sN is positive for every t in (0,1].
         call cauchy_length_along_the_banded_gradient(rAg, nAg, tau_c,     &
                                                      cauchy_ascends)
         if (cauchy_ascends) then
            tau_c = 0.0d0
         endif
         sU  = -tau_c*g
         AsU = -tau_c*Ag
      endif
      ! HOW MUCH OF THE STEP THE LEG CAN BE, before the dogleg is built on
      ! it: its length against the radius it has to share with the Krylov
      ! leg (cauchy_leg_min_fraction).
      nU_full = sqrt(sum(sU*sU))
      ! With no Krylov leg the dogleg's second leg is the Cauchy point
      ! itself, so the path is the Cauchy direction cut to the ball. If the
      ! Cauchy leg is missing too -- an empty gradient, an unsamplable
      ! probe, or a gradient that ascends the model -- there is no direction
      ! at all and the exit says which.
      if (krylov_leg_missing) then
         if (sqrt(sum(sU*sU)) .le. 0.0d0) then
            model_ok = .false.
            refuse = tr_no_krylov_direction
            f2_ref_out = f2_ref
            weno_mode = saved_mode
            call release_active_species_bounds
            return
         endif
         sN = sU;  AsN = AsU
         n_cauchy_only_steps = n_cauchy_only_steps + 1
      endif
      ! IS EITHER LEG KNOWN TO REDUCE THE MODEL? The Cauchy leg was dropped
      ! because it ascends, and the argument that the Krylov leg descends
      ! rests on the cycle having reduced the residual of the model,
      ! ||r0 + A sN|| < ||r0||. Where the cycle's own achieved relative
      ! residual is one or more that is not true, and building a dogleg on
      ! it would be building it on an assumption the cycle has already
      ! contradicted (R2).
      if (cauchy_ascends .and. gm_resid_rel .ge. 1.0d0) then
         model_ok   = .false.
         refuse     = tr_no_direction_reduces_model
         f2_ref_out = f2_ref
         weno_mode  = saved_mode
         call release_active_species_bounds
         return
      endif

      ! WHAT THE MODEL PROMISES ALONG EACH LEG, BEFORE ANY TRIAL. The three
      ! predicted decreases are the same expression the ratio test uses
      ! (predicted_model_decrease), taken at the FULL Cauchy point, at the
      ! full Krylov point and at the dogleg point of the entry radius. They
      ! cost no product: the images of both legs are already formed. Read
      ! against the merit itself they say whether the region, the leg, or
      ! the direction is what limits the step.
      if (elem_diag_here) then
         call dogleg_step(sU, sN, delta, aU_diag, aN_diag, on_boundary_diag)
         pred_cauchy = predicted_model_decrease(r0, AsU)
         pred_krylov = predicted_model_decrease(r0, AsN)
         pred_dogleg = predicted_model_decrease(r0,                       &
                            aU_diag*AsU + aN_diag*AsN)
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')                 &
              ' (JFNK) [diag 6] ||g|| ', gnorm, ', r0 . A g ', rAg,       &
              ', ||A g|| ', sqrt(max(nAg,0.0d0)), ', tau ', tau_c
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')        &
              ' (JFNK) [diag 6]   ||sU|| ', nU_full,                      &
              ', ||sU||/radius ', nU_full/max(delta, 1.0d-300),           &
              ', ||sN|| ', sqrt(sum(sN*sN)), ', radius ', delta,          &
              ', merit ', f2_ref
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')        &
              ' (JFNK) [diag 6]   predicted decrease: Cauchy ',           &
              pred_cauchy, ', Krylov ', pred_krylov, ', dogleg ',         &
              pred_dogleg, ' at aU ', aU_diag, ', aN ', aN_diag
      endif

      ! --- the dogleg, then admissibility, then the ratio ---
      do it = 1, n_theta_max
         ! THE LEG THAT IS NOT IN THE STEP IS NOT IN THE MODEL OF IT EITHER.
         ! Where the approximate-gradient leg is shorter than
         ! cauchy_leg_min_fraction of the radius it moves the step by less
         ! than the radius resolves, while its image A sU need not be small
         ! at all: dropping it from the step and keeping it in the model
         ! slope is what made the model slope describe a direction the step
         ! did not take (N7b). It is dropped from both, and the dogleg
         ! reduces to the Krylov leg cut to the ball. The test is taken at
         ! the CURRENT radius, which shrinks inside this loop, so a leg that
         ! is short beside a large region is used again beside a small one.
         ! With no Krylov leg the Cauchy point is the whole step and is kept
         ! whatever its length: it is the direction every convergence
         ! guarantee of a trust region rests on.
         cauchy_leg_short = .false.
         if (.not. krylov_leg_missing) then
            cauchy_leg_short = cauchy_leg_is_dropped_on_its_length(       &
                    nU_full, delta, cauchy_leg_admitted_by_image_alone)
         endif
         if (cauchy_leg_short) then
            sU = 0.0d0;  AsU = 0.0d0
         else
            sU  = -tau_c*g
            AsU = -tau_c*Ag
         endif
         if (cauchy_leg_short) n_cauchy_leg_short = n_cauchy_leg_short + 1
         call dogleg_step(sU, sN, delta, aU, aN, on_boundary)
         ! AND THE LEG THAT IS ALL OF THE IMAGE AND NONE OF THE STEP IS
         ! DROPPED TOO. The test above compares the leg's LENGTH with the
         ! radius; what decides whether the model slope r0 . A s describes
         ! the step is the leg's share of the IMAGE, and the two part
         ! company wherever ||A g|| stands decades above r0 . A g: the leg
         ! is then long enough for the length test to keep it and still
         ! supplies the whole image, so the slope handed to the ray test
         ! belongs to a direction the step barely contains and the
         ! predicted decrease carries no radius
         ! (cauchy_leg_image_share_max, cauchy_leg_step_share_max). The
         ! dogleg is rebuilt without the leg, which costs no product: both
         ! images are formed already.
         cauchy_step_share   = 0.0d0
         cauchy_image_share  = 0.0d0
         cauchy_leg_all_image = .false.
         if (.not. krylov_leg_missing .and. .not. cauchy_leg_short) then
            call cauchy_leg_shares_of_step_and_image(sU, sN, AsU, AsN,    &
                     aU, aN, cauchy_step_share, cauchy_image_share)
            cauchy_leg_all_image =                                        &
                 cauchy_leg_is_image_without_step(cauchy_step_share,      &
                                                  cauchy_image_share)
         endif
         if (cauchy_leg_all_image) then
            sU = 0.0d0;  AsU = 0.0d0
            n_cauchy_leg_all_image = n_cauchy_leg_all_image + 1
            call dogleg_step(sU, sN, delta, aU, aN, on_boundary)
         endif
         ! THE MODEL'S PROMISE AGAINST THE TRUE DECREASE ALONG THE STEP,
         ! AT TWO RADII, WITH THE LEG ADMITTED AND WITH IT DROPPED.
         ! The predicted decrease -(r0 . A s) - 0.5 ||A s||^2 is
         ! homogeneous of degree one in the step through its first term, so
         ! a step cut by four must carry a promise cut by four. A leg that
         ! is the whole image of the step breaks that: the image is
         ! tau_c ||A g||, a length of its own that the radius does not
         ! enter, so the promise stands still while the step shrinks. The
         ! four readings here separate the two, in one run, at one iterate:
         ! for the radius and a quarter of it, the dogleg is built with the
         ! leg and without it, the true decrease is taken from a discarded
         ! evaluation of the trial state, and the ratio of the two is
         ! printed. The step is held where it leaves the species box, as the
         ! trial path holds it, but the shared element budget is not
         ! restored: this is a reading of the model, not a candidate state.
         if (elem_diag_every_trial .and. it .eq. 1) then
            allocate(dl_sU(nvar_jac*N), dl_AsU(nvar_jac*N),               &
                     dl_s(nvar_jac*N),  dl_As(nvar_jac*N))
            call hold_solve_refusal_statistics(statistics_at_entry)
            do k_radius = 1, 2
               dl_delta = delta/4.0d0**(k_radius - 1)
               do k_leg = 1, 2
                  if (k_leg .eq. 1) then
                     dl_sU  = -tau_c*g
                     dl_AsU = -tau_c*Ag
                  else
                     dl_sU  = 0.0d0
                     dl_AsU = 0.0d0
                  endif
                  call dogleg_step(dl_sU, sN, dl_delta, dl_aU, dl_aN,     &
                                   dl_bnd)
                  dl_s  = dl_aU*dl_sU  + dl_aN*sN
                  dl_As = dl_aU*dl_AsU + dl_aN*AsN
                  dl_pred = predicted_model_decrease(r0, dl_As)
                  Ytry    = Y + D*dl_s
                  n_dl_out = 0
                  if (project_species_trial)                              &
                     call hold_the_step_where_it_leaves_the_species_box(Y,&
                                                       Ytry, n_dl_out)
                  dl_actual = 0.0d0
                  dl_f2     = -1.0d0
                  call unpack_U(Ytry, utry)
                  call U_to_W_interior(utry, Wtry)
                  dl_ok = (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.           &
                           minval(Wtry(3,1:N)) .gt. 0.0d0)
                  if (dl_ok) then
                     weno_mode = 0
                     call eval_residual(Ytry, f_sp, fwork, Ftry, heat,    &
                                        cool, admissible=okres,           &
                                        state_is_discarded=.true.)
                     dl_ok = okres
                     if (okres) then
                        dl_f2     = sqrt(sum((Ftry/Drow)**2))
                        dl_actual = 0.5d0*(f2_ref*f2_ref - dl_f2*dl_f2)
                     endif
                  endif
                  write(*,'(A,A,A,ES11.3,A,ES14.6,A,ES14.6,A,ES11.3,'//   &
                          'A,ES11.3,A,I0,A,L1)')                           &
                       ' (JFNK) [diag 9] leg ',                            &
                       merge('admitted', 'dropped ', k_leg .eq. 1),        &
                       ', radius ', dl_delta, ', pred ', dl_pred,          &
                       ', actual ', dl_actual, ', ratio ',                 &
                       dl_actual/max(abs(dl_pred), 1.0d-300)               &
                       *sign(1.0d0, dl_pred), ', ||s|| ',                  &
                       sqrt(sum(dl_s*dl_s)), ', unknowns held ',           &
                       n_dl_out, ', describable ', dl_ok
               enddo
            enddo
            call put_back_solve_refusal_statistics(statistics_at_entry)
            deallocate(dl_sU, dl_AsU, dl_s, dl_As)
         endif
         if (tr_trace) write(*,'(A,I3,A,ES11.3,A,L1)')                    &
              ' (JFNK) [TR-trace] trial', it, '  delta_in=', delta,        &
              '  on_boundary=', on_boundary
         s  = aU*sU  + aN*sN
         As = aU*AsU + aN*AsN
         snorm = sqrt(sum(s*s))
         if (snorm .le. 0.0d0) then
            model_ok = .false.;  refuse = tr_dogleg_step_is_zero;  exit
         endif
         ! Scale back until the trial is a state the code can describe. The
         ! model is linear, so the predicted reduction of the scaled-back step
         ! is recomputed rather than scaled.
         theta = 1.0d0
         ok_state = .false.
         n_out = 0
         do while (theta .ge. 1.0d0/2.0d0**n_theta_max)
            Ytry = Y + D*(theta*s)
            ! THE STEP IS HELD AT WHATEVER LEAVES THE BOX FIRST, then
            ! shortened for whatever else the state fails: the two are
            ! different obstructions and only one of them yields to a
            ! shorter step. An unknown the step carries out of its box
            ! keeps the iterate's value, so the trial is a step the
            ! operator can be sampled along and the model is a model of it
            ! (hold_the_step_where_it_leaves_the_species_box).
            if (project_species_trial) then
               call hold_the_step_where_it_leaves_the_species_box(Y, Ytry, &
                                                                  n_out)
               ! AND THE NONLINEAR SHARED CONSTRAINT, RESTORED BEFORE THE
               ! MERIT IS JUDGED. The rows the step was projected on are a
               ! LINEARIZATION of a constraint whose budget depends on the
               ! cell's own density and element unknowns, so a candidate can
               ! violate it however carefully the direction was projected.
               ! A candidate that violates by more than
               ! element_constraint_restoration_trigger of what the row
               ! allows is stepped back onto the constraint, measured by the
               ! map's violation magnitude, and only then is its merit taken
               ! (restore_the_element_budget). The restoration moves the
               ! trial off the plane the dogleg was built in, exactly as the
               ! box projection does, so it is counted into n_out and the
               ! model of the step is recomputed below.
               call restore_the_element_budget(Ytry, .true., n_rest_cell,  &
                                    n_rest_row, rest_before, rest_after)
               if (n_rest_cell .gt. 0) then
                  call hold_the_step_where_it_leaves_the_species_box(Y,    &
                                                    Ytry, n_rest_out)
                  n_out = n_out + n_rest_cell + n_rest_out
               endif
            endif
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                     &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               weno_mode = 0
               ! A POINT ON THE DOGLEG, not a state anyone holds yet: if the
               ! ratio test then accepts it, it is evaluated again below as
               ! the adopted state.
               call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,    &
                                  admissible=okres,                       &
                                  state_is_discarded=.true.)
               if (okres) then
                  f2t = sqrt(sum((Ftry/Drow)**2))
                  ok_state = .true.
                  exit
               endif
            endif
            theta = 0.5d0*theta
            n_theta_cuts = n_theta_cuts + 1
         enddo
         if (.not. ok_state) then
            ! nothing on this dogleg is describable; shrink the region
            delta = 0.25d0*delta
            if (delta .lt. delta_floor) then
               model_ok = .false.
               refuse = tr_no_admissible_state_on_ray
               exit
            endif
            cycle
         endif
         Ftrial = Ftry
         ! THE VECTOR THAT WAS EVALUATED IS THE VECTOR THAT IS ADOPTED, and
         ! it is kept rather than rebuilt. Rebuilding it as Y + D s puts a
         ! projected component back through a division and a multiplication
         ! by the same scale, and the round trip does not close: MEASURED on
         ! the coupled `mol_carrier` reload, the acceptance re-evaluation
         ! then found the carrier of cell 258 at -4.815e-35 on a trial whose
         ! own evaluation had it at exactly zero, and six consecutive outer
         ! iterations ended on "the accepted trial is no longer admissible".
         ! For a step that was not projected this is the same bits: s already
         ! holds theta times the dogleg step and Y + D s is how Ytry was
         ! formed.
         Ystep = Ytry
         if (n_out .gt. 0) then
            ! THE STEP THE TRIAL IS, AND THE MODEL OF THAT STEP. The
            ! projected trial is not Y + D theta s any more, so the
            ! predicted reduction may not be formed from theta*As: the
            ! projected step leaves the plane the dogleg was built in and
            ! A applied to it is not a combination of A sU and A sN. One
            ! matrix-free product buys the model of the step actually taken,
            ! and without it the ratio test would compare the reduction of
            ! one step with the promise made about another.
            s     = (Ytry - Y)/D
            snorm = sqrt(sum(s*s))
            call jv_product(Y, F_jac, f_sp, D*s, Jv, ok_state)
            if (.not. ok_state) then
               model_ok = .false.
               refuse = tr_projected_step_has_no_model
               delta = 0.25d0*snorm
               if (delta .lt. delta_floor) exit
               cycle
            endif
            As = Jv/Drow
            n_projected_trials = n_projected_trials + 1
         else
            s     = theta*s
            As    = theta*As
            snorm = theta*snorm
         endif
         pred   = predicted_model_decrease(r0, As)
         actual = 0.5d0*(f2_ref*f2_ref - f2t*f2t)
         ! THE MODEL OF THIS TRIAL, TERM BY TERM, AND THE SAME MODEL FROM A
         ! FRESH IMAGE OF THE SAME STEP. The predicted decrease is
         ! -(r0 . A s) - 1/2 ||A s||^2, homogeneous of degree one in the
         ! step through its first term, so a step cut by four must carry a
         ! first term cut by four. Three things can break that and they are
         ! separated here: an image A s carried over from a step of another
         ! length (||A s|| no longer proportional to ||s||), a directional
         ! derivative that has stopped resolving the direction (the fresh
         ! image disagreeing with the one the dogleg built), and
         ! cancellation between the two terms (their sizes printed
         ! separately). The fresh product is diagnostic only: the solve's
         ! refusal statistics are held aside so the step is the step the
         ! solve takes with the hook off.
         if (elem_diag_every_trial) then
            allocate(As_true(nvar_jac*N))
            call hold_solve_refusal_statistics(statistics_at_entry)
            call jv_product(Y, F_jac, f_sp, D*s, As_true, ok_true)
            call put_back_solve_refusal_statistics(statistics_at_entry)
            if (ok_true) then
               As_true = As_true/Drow
               pred_fresh  = predicted_model_decrease(r0, As_true)
               arnoldi_gap = sqrt(sum((As - As_true)**2))                 &
                           / max(sqrt(sum(As_true*As_true)), 1.0d-300)
            else
               pred_fresh = 0.0d0;  arnoldi_gap = -1.0d0
               As_true    = 0.0d0
            endif
            write(*,'(A,I3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,'//       &
                    'A,ES11.3,A,I0)')                                      &
                 ' (JFNK) [diag 8] trial', it, '  radius ', delta,         &
                 ', ||sN|| ', sqrt(sum(sN*sN)), ', aU ', aU, ', aN ', aN, &
                 ', theta at ', theta, ', unknowns held ', n_out
            ! AND THE TWO SHARES THE APPROXIMATE-GRADIENT LEG HAS OF THE
            ! STEP AND OF ITS IMAGE, beside which of the two rules removed
            ! it. The image share is what the model slope is a statement
            ! about; the step share is what the step contains
            ! (cauchy_leg_shares_of_step_and_image). They are taken here at
            ! the FULL leg whether or not a rule removed it, so that a leg
            ! the length test drops can still be read; both are -1 with no
            ! Krylov leg to share with.
            dl_step_share  = -1.0d0
            dl_image_share = -1.0d0
            if (.not. krylov_leg_missing) then
               allocate(dl_sU(nvar_jac*N), dl_AsU(nvar_jac*N))
               dl_sU  = -tau_c*g
               dl_AsU = -tau_c*Ag
               call dogleg_step(dl_sU, sN, delta, dl_aU, dl_aN, dl_bnd)
               call cauchy_leg_shares_of_step_and_image(dl_sU, sN,        &
                        dl_AsU, AsN, dl_aU, dl_aN, dl_step_share,         &
                        dl_image_share)
               deallocate(dl_sU, dl_AsU)
            endif
            write(*,'(A,ES11.3,A,ES11.3,A,L1,A,L1)')                      &
                 ' (JFNK) [diag 8]   leg share of the step ',             &
                 dl_step_share, ', of the image ', dl_image_share,        &
                 ', dropped on its length ', cauchy_leg_short,            &
                 ', on its image ', cauchy_leg_all_image
            write(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)')              &
                 ' (JFNK) [diag 8]   ||s|| ', snorm, ', ||A s|| ',        &
                 sqrt(sum(As*As)), ', r0 . A s ', sum(r0*As),             &
                 ', -0.5 ||A s||^2 ', -0.5d0*sum(As*As)
            write(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES11.3)')     &
                 ' (JFNK) [diag 8]   pred ', pred, ', actual ', actual,   &
                 ', pred from a fresh image ', pred_fresh,                &
                 ', ||A s|| fresh ', sqrt(sum(As_true*As_true)),          &
                 ', relative gap ', arnoldi_gap
            ! AND THE LENGTH THE DIRECTIONAL DERIVATIVE IS SAMPLED OVER,
            ! BESIDE THE LENGTH OF THE STEP IT IS THE MODEL OF. The probe
            ! length of the matrix-free action is
            ! sqrt(epsilon) (1 + ||Y||) as a length in the unscaled
            ! unknowns, whatever the direction (jv_product: the step is
            ! eps0 v with eps0 inversely proportional to ||v||). It is a
            ! property of the iterate and not of the step, so once the step
            ! is shorter than it the model of the step is a secant over an
            ! arc the step does not reach, and cutting the radius further
            ! cannot make the two agree.
            write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')                       &
                 ' (JFNK) [diag 8]   ||D s|| ', sqrt(sum((D*s)**2)),      &
                 ', probe length ',                                        &
                 probe_length_of_the_jacobian_action(Y),                  &
                 ', step over probe ', sqrt(sum((D*s)**2))                &
                 / max(probe_length_of_the_jacobian_action(Y), 1.0d-300)
            deallocate(As_true)
         endif
         if (pred .le. 0.0d0) then
            ! a model that promises no reduction is not a model to step on
            model_ok = .false.
            refuse = tr_model_promises_no_drop
            if (tr_trace) write(*,'(A,I3,A,ES11.3,A,ES11.3)')             &
                 ' (JFNK) [TR-trace] trial', it, '  pred<=0  pred=', pred, &
                 '  snorm=', snorm
            delta = 0.25d0*snorm
            if (delta .lt. delta_floor) exit
            cycle
         endif
         ratio = actual/pred

         ! --- the ray test: does the model's slope describe the function? ---
         mslope = sum(r0*As)
         hray   = ray_h
         weno_mode = 0
         Ytry = Y + D*(hray*s)
         call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,          &
                            admissible=okres,                             &
                            state_is_discarded=.true.)
         ok_ray = okres
         f2p = 0.5d0*sum((Ftry/Drow)**2)
         ! A CENTRAL DIFFERENCE DOES NOT EXIST ON A FACE OF THE FEASIBLE
         ! SET. Where the step moves a species unknown UP off its lower
         ! bound the minus ray moves it below that bound, at every hray, and
         ! there is no state there to evaluate. The one-sided quotient about
         ! the iterate is first order in hray where the central one is
         ! second order, and it exists; the alternative is the verdict "no
         ! sample", which shrinks the region against an obstruction a
         ! shorter ray does not remove.
         Ytry = Y - D*(hray*s)
         n_ray_out = 0
         if (project_species_trial) then
            call species_unknowns_outside_their_bounds(Ytry, n_ray_out,    &
                                                       .false.)
            ! A ray that spends nuclei the cell does not have is as
            ! undescribable as one with a negative density, and for the same
            ! reason: the chemistry of that cell has no root. It is counted
            ! and not restored, because a restored minus ray is no longer
            ! the reflection of the step.
            call restore_the_element_budget(Ytry, .false., n_rest_cell,    &
                                 n_rest_row, rest_before, rest_after)
            n_ray_out = n_ray_out + n_rest_row
         endif
         if (n_ray_out .gt. 0) then
            f2m    = 0.5d0*f2_ref*f2_ref
            fslope = (f2p - f2m)/hray
         else
            call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,       &
                               admissible=okres,                          &
                               state_is_discarded=.true.)
            ok_ray = ok_ray .and. okres
            f2m = 0.5d0*sum((Ftry/Drow)**2)
            fslope = (f2p - f2m)/(2.0d0*hray)
         endif
         call ray_slope_verdict(mslope, fslope, f2p, f2m, ok_ray, ratio,  &
                                merit_reproducibility, ray_floor,         &
                                ray_verdict, model_ok, refuse)
         ! WHAT THE RAY TEST FOUND WHERE IT WAS NOT ACTED ON. A verdict
         ! against the model beside a ratio at unity is a statement about
         ! the PROBE, and it is the statement that says whether the finite
         ! difference is worth its two evaluations on this operator, so it
         ! is printed rather than discarded.
         if (model_ok .and. ray_verdict .ne. tr_step_taken) then
            call tr_exit_text(ray_verdict, why_txt)
            write(*,'(A,A)') ' (JFNK) [TR]      the ray test is not'//    &
                 ' acted on, the ratio being sound: ', trim(why_txt)
            write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')              &
                 ' (JFNK) [TR]      model slope ', mslope,                &
                 ', ray slope ', fslope, ', ray floor ', ray_floor,       &
                 ', |f2p - f2m| ', abs(f2p - f2m)
         endif
         if (.not. model_ok) then
            delta = 0.25d0*snorm
            exit
         endif

         ! --- radius update and acceptance ---
         if (ratio .lt. 0.25d0) then
            delta = 0.25d0*snorm
         else if (ratio .gt. 0.75d0 .and. on_boundary) then
            delta = 2.0d0*delta
         endif
         if (ratio .gt. eta_accept) then
            ! Re-evaluate the accepted trial so that it is the LAST state this
            ! routine evaluated: the caller's keep_background_of_adopted_state,
            ! carrier_relnorm_last and the acceptance ledgers all read what the
            ! previous evaluation left behind, and the ray test above has been
            ! evaluating points that are not the step.
            weno_mode = 0
            Yacc = Ystep
            call eval_residual(Yacc, f_sp, f_sp_acc, Facc, heat, cool,    &
                               admissible=okres, rows_judged=rows_acc)
            ! WHAT IS NOT DONE HERE, and it was tried. Refusing the step
            ! when distance_from_certification(rows_acc) rises above jref
            ! -- the merit's own five-iterate window applied to the judged
            ! rows -- was implemented and MEASURED on the hot Uranus
            ! element solve: the region was shrunk by 13 such refusals and
            ! the solve stalled at ||R|| = 1.92 against the 4.345e-09 it
            ! reaches without them (report B5e section 3.3). A step control
            ! bounds a smooth merit; the judged measure is a maximum over
            ! cells and over rows, and it belongs to the ledger that
            ! chooses the returned state.
            if (okres) then
               f2t = sqrt(sum((Facc/Drow)**2))
               dY = Yacc - Y;  f2_new = f2t;  accepted = .true.
            else
               ! it was admissible a moment ago and is not now; refuse it
               model_ok = .false.
               refuse = tr_trial_no_longer_admissible
               delta = 0.25d0*snorm
            endif
         else
            refuse = tr_ratio_below_acceptance
         endif
         exit
      enddo
      f2_ref_out = f2_ref
      weno_mode = saved_mode
      call release_active_species_bounds
      ! AND WHAT THE STEP ACTUALLY DID, against the promise printed above,
      ! beside the samples this one step was refused.
      if (elem_diag_here) then
         call tr_exit_text(refuse, why_txt)
         write(*,'(A,ES11.3,A,ES11.3,A,F9.4,A,ES11.3)')                   &
              ' (JFNK) [diag 6]   pred ', pred, ', actual ', actual,      &
              ', ratio ', ratio, ', ||s|| ', snorm
         write(*,'(A,A,A,I0)') ' (JFNK) [diag 6]   exit: ', trim(why_txt),&
              ', samples refused during this step ',                      &
              sum(n_eval_refusal(1:6)) - n_refused_at_entry
         ! AND WHAT THE CONSTRAINED STEP DID WITH THE SHARED ROWS: how many
         ! of them the projection was taken on, how many cells the
         ! restoration moved, and the violation magnitude it removed.
         write(*,'(A,I0,A,I0,A,ES11.3,A,ES11.3)')                         &
              ' (JFNK) [diag 6]   active shared rows ',                   &
              n_active_element_rows(), ', cells restored ', n_rest_cell,   &
              ', worst violation before ', rest_before, ', after ',       &
              rest_after
         ! AND WHAT THE RAY TEST COMPARED, which is what decides the exit
         ! whenever the ratio itself is sound: the model's slope along the
         ! step, the difference quotient of the merit along the same ray,
         ! and the cancellation floor of that quotient. A ratio at unity
         ! beside two slopes that differ by more than the factor the test
         ! allows is a statement about the PROBE and not about the model.
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)')                 &
              ' (JFNK) [diag 6]   model slope ', mslope,                  &
              ', ray slope ', fslope, ', ray floor ', ray_floor,          &
              ', |f2p - f2m| ', abs(f2p - f2m)
         if (allocated(g_unheld)) deallocate(g_unheld)
      endif

      ! WHAT THE REFUSING STATE LOOKED LIKE.  The radius alone cannot say
      ! whether the region shrank because the step left the describable set
      ! or because the model over-predicted on a perfectly ordinary state,
      ! so the trial's own worst scaled row is printed beside the two
      ! slopes the ray test compared.
      if (tr_trace .and. refuse .ne. tr_step_taken) then
         call tr_exit_text(refuse, why_txt)
         worst_row = -1.0d0;  i_worst = 1
         do it = 1, nvar_jac*N
            if (abs(Ftrial(it))/Drow(it) .gt. worst_row) then
               worst_row = abs(Ftrial(it))/Drow(it);  i_worst = it
            endif
         enddo
         call name_of_unknown(i_worst, what_worst)
         write(*,'(A,A)') ' (JFNK) [TR-why] ', trim(why_txt)
         write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,I0)')                      &
              ' (JFNK) [TR-why]   delta in', delta_in, ' out', delta,      &
              '  ||s||', snorm, '  theta cuts ', n_theta_cuts
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')                           &
              ' (JFNK) [TR-why]   pred', pred, '  actual', actual,         &
              '  model slope', mslope
         write(*,'(A,ES11.3,A,L1)')                                        &
              ' (JFNK) [TR-why]   r0 . A g', rAg,                          &
              '  Cauchy leg dropped ', cauchy_ascends
         write(*,'(A,ES11.3,A,ES11.3,A,ES11.3,A,ES10.3)')                  &
              ' (JFNK) [TR-why]   true slope', fslope, '  f2p', f2p,       &
              '  f2m', f2m, '  ray floor', ray_floor
         write(*,'(A,ES11.3,A,A)')                                         &
              ' (JFNK) [TR-why]   worst scaled row of the refusing trial', &
              worst_row, ' at ', trim(what_worst)
      endif
      end subroutine trust_region_step

      ! ------------------------------------------------------!

      subroutine solve_steady_jfnk(u, f_sp, resid_tol, maxit, dtau0,    &
                                   gm_m, info)
      ! Pseudo-transient-continuation Jacobian-free Newton-Krylov solve.
      ! Per outer iteration:
      !   F  = eval_residual(Y)                          (full residual)
      !   M  = factored banded (I/dtau + J_banded)       (preconditioner)
      !   solve (I/dtau + J) dY = -F  by right-precond. GMRES (J*v matrix-
      !                                free; captures non-local radiation)
      !   Y <- Y + lam*dY   (||F||_2 line search + positivity)
      !   dtau <- SER ramp
      !
      ! info: 0 the returned state met the acceptance gate and every judged
      ! row; 1 the iteration limit was reached; 2 the solve ended without a
      ! state that meets the gate; 3 an iterate at which no feasible descent
      ! direction exists -- no Krylov direction and a projected gradient of
      ! zero -- so there is nothing to size a step on
      ! (initial_trust_region_radius).
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit, gm_m
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc, gmit
      ! The outer iteration cap this solve runs under, the caller's argument
      ! unless EXHALE_JFNK_MAXIT named another (jfnk_outer_iteration_cap).
      integer :: maxit_used
      ! THE KRYLOV SUBSPACE SIZE THIS SOLVE USES.  gm_m is what the caller
      ! passed; EXHALE_GM_M overrides it so that the dependence of the step
      ! on the size of the subspace can be MEASURED without a rebuild, the
      ! caller's value being the default.  A restarted GMRES with a subspace
      ! too small for the spectrum returns a step that minimizes the
      ! residual over a subspace and not the model's own minimizer, and the
      ! dogleg falls monotonically along its path only when the Newton leg
      ! IS that minimizer (section 159.3), so the size is part of the step
      ! control and not only of its cost.
      integer :: gm_m_used
      character(len=32) :: gm_env
      ! The relative residual the Newton leg's cycle reached, per outer
      ! iteration, against the tolerance it was asked for.
      real*8  :: tr_gm_resid
      integer :: jj, kk, jworst, kworst, irow, ilo, ihi, isr
      real*8  :: amx
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      ! THE BASE POINT OF THE NEWTON MODEL: the ONE-PASS residual at the
      ! current iterate, about which every finite difference of the banded
      ! Jacobian and of the Krylov products is taken (n_eq_sweeps_model). F
      ! itself is the residual at the iterate's own composition and is the
      ! right hand side, the merit and the gate; the system solved is
      ! therefore J_1 dY = -F, an inexact Newton whose root is F's root.
      real*8, allocatable :: F_jac(:)
      real*8, allocatable :: dZ(:), D(:), dY_lm(:), Drow(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(1-Ng:N+Ng)            :: npart0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j, f_sp_best
      real*8  :: rnorm, rc(3), dtau, lam, f2, f2_try, idtau, amx0
      real*8  :: f2hist(5), f2ref, rnorm_best
      ! Flux gate: the spread of rho v r^2 over r >= r_flux, of the current
      ! iterate and of the best one kept, and whether each meets BOTH gates.
      real*8  :: fspread, fspread_best
      real*8  :: fspread_103, fspread_110, dum_fmean
      logical :: gates_now, gates_best
      integer :: n_no_descent, n_nonmonotone_accepts
      integer :: n_since_best
      logical :: monotone_search
      ! Stagnation limit: consecutive outer iterations in which the line
      ! search found NO acceptable step at all. Measured on the HD 209458 b
      ! hand-off state (docs/newton_scaling_and_base_wall.md): runs that go on
      ! to converge never string more than 4 such iterations together, runs
      ! that are truly stuck string 34 or more.
      integer, parameter :: n_no_descent_max = 12
      ! P51. STAGNATION OF THE BEST ITERATE, and what to do about it.
      !
      ! The non-monotone (Grippo) line search compares a trial against the
      ! WORST merit of the last five iterates, so it accepts sideways steps.
      ! That is what lets the solve leave a bad neighborhood, and it is also
      ! what lets it circle: measured on the coupled hot Uranus hand-off, the
      ! best iterate is reached at outer iteration 10 and then 490 further
      ! iterations produce no improvement of ||R|| at all while the search
      ! accepts 61 sideways steps.
      !
      ! So when the best iterate has not improved for n_stall_best
      ! iterations the solve GOES BACK TO IT and switches the acceptance to
      ! monotone Armijo -- descent against the current merit rather than
      ! against the worst of five. If that still buys nothing in another
      ! n_stall_best iterations, the solve stops and says it STAGNATED,
      ! which is a different statement from "no descent direction exists".
      !
      ! WHERE 20 COMES FROM. Over every solve of the P50 campaign that ended
      ! info = 0 on a configuration that survived measurement, the largest
      ! gap between two improvements of the best iterate is 6 outer
      ! iterations (both WASP-121b Newton cases); on the coupled hot Uranus
      ! it is 3. Two accepted solves show gaps of 68 and 85, and both are
      ! configurations that were measured and rejected. So 20 stands a
      ! factor 3.3 above anything a surviving route has needed, and 40 --
      ! the point at which the solve gives up -- stands a factor 6.7 above
      ! it and below every stall observed (tails of 278 to 729).
      integer, parameter :: n_stall_best = 20
      ! THE ACCEPTANCE BOUNDS THE QUANTITY THE STATE IS JUDGED BY, beside
      ! the merit it descends on.
      !
      ! The two are different functionals, and until B5e only the first was
      ! bounded: the merit is || F/Drow ||_2, a 2-norm on the Newton's own
      ! column scales, while info = 0 rests on a MAXIMUM over cells of each
      ! row against its own largest term (resid_relnorm and the
      ! certification's row entries). A Grippo window on the merit alone
      ! therefore places no bound at all on the number the state is judged
      ! by: MEASURED on the hot Uranus reload with an element row, the
      ! default route drove the merit from 1.44e+02 down to 3.86e+01 while
      ! ||R|| ran 1.72, 6.44, 1.82, 25.2, 1.93, 5.09 -- fifteenfold
      ! excursions of the judged quantity, every one of them an accepted
      ! step (B5d section 2.7).
      !
      ! WHERE THE BOUND IS, AND WHERE IT IS NOT. Making it a condition of
      ! the STEP -- a trial refused when it leaves the state further from
      ! certification than the worst of the last five iterates, the merit's
      ! own window applied to the judged rows -- was implemented and
      ! MEASURED: on the element solve above it shrank the trust region 13
      ! times and left the solve stalled at ||R|| = 1.92 against the
      ! 4.345e-09 the same text reaches without it (report B5e section 3.3).
      ! A step control bounds a SMOOTH merit; the judged measure is a
      ! maximum over cells and over rows, and a solve legitimately trades one
      ! row against the others on its way down.
      !
      ! So the judgement enters at the LEDGER that chooses the state to hand
      ! back, which is the place the solve makes its claim: the best iterate
      ! is ranked on distance_from_certification, the same rows the
      ! certification reads and each against its own tolerance, instead of
      ! on ||R||, which is the largest of the three hydrodynamic row
      ! measures with no tolerance in it. jhist and jref are kept for the
      ! count below, which measures how far apart the two functionals ran.
      real*8  :: jhist(5), jref, dj, dj_try, dj_best
      real*8  :: rows_j(3), rows_try(3)
      ! |R| of the iterate's own evaluation, beside the judged rows of the
      ! same evaluation: the gate reads this and the ledger reads the
      ! distance, and the iteration line prints both so that a reader can
      ! tell which of the two a solve was held by.
      real*8  :: rnorm_j
      logical :: keep_this_iterate
      integer :: n_judged_excursions
      logical :: ok, try_ok, carrier_ok
      ! --- the scaled trust region: THE STEP CONTROL OF THE COUPLED ROUTE ---
      !
      ! It is the step control wherever the system carries a species row, and
      ! the three-unknown route never reaches one statement of it.
      ! EXHALE_TRUST_REGION=0 restores the pseudo-transient line search there
      ! for measurement.
      !
      ! WHY IT IS THE DEFAULT AND THE LINE SEARCH IS NOT. MEASURED on the
      ! 1 microbar hot Uranus reloaded with He_diffusion and
      ! `Coupled carrier solve: True`, i.e. one element row and nvar = 4,
      ! from the same 300-step marching snapshot and with everything else
      ! equal:
      !
      !   line search   ||R|| 1.986 -> 1.723, aborted at outer iteration 18
      !                 on "no descent direction exists for the banded
      !                 model", the accepted steps running ||R|| 1.72, 6.44,
      !                 1.82, 25.2 (B5d section 2.7);
      !   trust region  ||R|| 1.986 -> 4.345e-09 against a target of 1e-08,
      !                 54 accepted steps in 66 outer iterations, and the
      !                 elemental transport row of the returned state
      !                 4.633e-11, at the 1e-08 the elemental transport
      !                 tolerance carried then.
      !
      ! No case of the regression matrix carries a species row, so this
      ! default moves no golden; what it changes is the route a coupled solve
      ! takes, and it changes it from one that does not converge to one that
      ! does.
      logical :: use_tr
      real*8  :: tr_delta, tr_pred, tr_actual, tr_ratio, tr_snorm, tr_dmax
      ! Relative tolerance of the Krylov solve that builds the dogleg's
      ! Newton leg.  The dogleg falls monotonically along its path only
      ! when that leg is the model's own minimizer, so a loose leg makes
      ! the model predict an increase at full step length (section 159.3).
      ! Overridable so the dependence can be measured; the default is the
      ! value the section 153 draft used.
      real*8  :: tr_gm_rtol
      real*8  :: tr_f2ref, tr_f2_end
      integer :: tr_cuts, n_tr_accept, n_tr_reject, n_tr_model_bad
      ! One counter per named exit of trust_region_step, so the summary line
      ! reports why the region shrank and not merely how often.
      integer :: tr_refuse, n_tr_exit(0:n_tr_exit_reason)
      character(len=80) :: tr_why_txt
      logical :: tr_model_ok, tr_accept, tr_cauchy_ascends
      ! Outer iterations at which the banded gradient ascended the model,
      ! so the Cauchy leg was dropped and the dogleg reduced to the Krylov
      ! leg cut to the ball.
      integer :: n_tr_cauchy_dropped
      character(len=24) :: tr_env
      ! Damped Gauss-Newton escape from an ascending Newton direction.
      logical :: lm_tried, lm_found
      real*8  :: f2_lm, mu_lm, gnorm_lm
      ! Probe states of THIS solve that could not be described: Jacobian
      ! colors left at zero, and GMRES cycles cut short because the next
      ! Krylov direction had no admissible sample (section 121).
      integer :: n_jac_unresolved, n_jac_unresolved_tot
      ! HOW THE LINEAR SOLVE OF EACH OUTER ITERATION ENDED, by name
      ! (gm_outcome_text). A solve that ends without a certified state
      ! failed for one of these reasons and the iteration counts do not say
      ! which.
      integer :: gm_outcome, n_gm_outcome(0:n_gm_outcome_kind)
      real*8  :: gm_snorm_last
      character(len=80) :: gm_why_txt
      ! The cell and the quantity of the largest scaled residual of the
      ! iteration line (name_of_unknown).
      character(len=48) :: what_worst_row
      ! The residual and the projected model gradient the FIRST radius is
      ! measured on, and which direction it came from
      ! (initial_trust_region_radius).
      real*8, dimension(nvar_jac*N) :: r0_init, g_init
      integer :: n_active_init, tr_init_status
      real*8, allocatable :: Ybest(:)
      ! The certification of the state handed back, read back so that the
      ! completion flag and the printed refusal are the same record.
      type(cert_report) :: cert_out
      logical :: rows_ok, rows_unmeasured, gate_is_cert
      integer :: i_row
      ! How many rows of the certification the stop test read at the iterate
      ! it stopped on, and the outer iterations that met the acceptance gate
      ! while one of those rows refused: the iterations the solve spends on a
      ! row the gate does not see.
      integer :: n_rows_at_stop, n_gate_met_row_refused, i_row_refused_last
      ! ||R|| of the accepted iterate, kept so that the state handed back can
      ! be compared with the number that iterate carried.
      real*8  :: rnorm_handback_prev
      ! WHETHER THE CLOSURE MAP AT THE BASE IS TO BE RE-MEASURED at the top
      ! of the next outer iteration, where that iterate's own row scales
      ! exist (a restart's third reset, restart_the_trust_region_state).
      logical :: closure_map_due
      ! And whether the Grippo window and the monotone flag are to be
      ! collapsed onto the next iterate, for the same reason: the window is
      ! a bound about the state the acceptance compares against.
      logical :: acceptance_window_due
      ! Whether a restart of the trust-region state is due at the end of
      ! this outer iteration, and the sentence naming what it reset.
      logical :: tr_restart_due
      integer :: i_reset
      character(len=64) :: reset_txt

      neq  = nvar_jac*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(F_jac(neq))
      allocate(dZ(neq), D(neq), Ybest(neq), dY_lm(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      ! Allocated last and on its own, so that every array the three-unknown
      ! solve has always had keeps the address it has always had.
      allocate(Drow(neq))
      dtau = dtau0;  info = 1
      resid_tol_of_solve = resid_tol
      gm_m_used = gm_m
      call get_environment_variable('EXHALE_GM_M', gm_env)
      if (len_trim(gm_env) .gt. 0) read(gm_env,*) gm_m_used
      gm_m_used = max(gm_m_used, 1)
      tr_gm_resid = 0.0d0

      n_jac_unresolved_tot = 0;  n_gm_outcome = 0
      n_rows_at_stop = 0;  n_gate_met_row_refused = 0
      i_row_refused_last = -1
      gm_outcome = gm_tolerance_reached;  gm_snorm_last = 0.0d0
      use_tr = (nspec_row .gt. 0)
      call get_environment_variable('EXHALE_TRUST_REGION', tr_env)
      if (trim(tr_env) .eq. '0') use_tr = .false.
      n_tr_accept = 0;  n_tr_reject = 0;  n_tr_model_bad = 0
      n_tr_exit = 0;  tr_refuse = tr_step_taken
      n_tr_cauchy_dropped = 0;  tr_cauchy_ascends = .false.
      tr_f2_end = -1.0d0
      tr_delta = -1.0d0
      tr_dmax  = -1.0d0
      tr_gm_rtol = 1.0d-1
      call get_environment_variable('EXHALE_TR_GMRTOL', tr_env)
      if (len_trim(tr_env) .gt. 0) read(tr_env,*) tr_gm_rtol
      call get_environment_variable('EXHALE_TR_DELTA0', tr_env)
      if (len_trim(tr_env) .gt. 0) read(tr_env,*) tr_delta
      call reset_no_chem_root_trial_count
      call read_composition_elimination_controls
      call read_species_unknown_space_controls
      call read_trust_region_restart_controls
      closure_map_due = .false.;  acceptance_window_due = .false.
      maxit_used = jfnk_outer_iteration_cap(maxit)

      call pack_U(u, Y)
      call pack_species_rows(u, f_sp, Y)
      n_trial_headroom = 0;  n_resid_eval = 0
      n_eval_refusal = 0;  n_jv_backward_sample = 0
      n_budget_face_at_the_iterate = 0
      n_elem_row_dropped_for_rank = 0
      n_elem_row_target_at_the_iterate = 0
      n_elem_restoration_cells = 0;  n_elem_restoration_steps = 0
      elem_restoration_worst_before = 0.0d0
      elem_restoration_worst_after  = 0.0d0
      n_write_back_element_deficit = 0
      write_back_worst_element_deficit = 0.0d0
      n_probe_blocked_components = 0
      ! From here to the end of the solve the equilibrium sweeps are the
      ! steady solver's, not the run's. The current iterate is a state the
      ! run holds; everything the solver evaluates to BUILD its Newton model
      ! -- Jacobian columns, Krylov products, line-search trials -- is not,
      ! and is tagged separately so that the acceptance ledgers, the non-root
      ! streak and its stop are not written by states no one adopted
      ! (docs/Update_EXHALE_stage1.md section 121).
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      ! The invariant the acceptance of this solve rests on, measured on the
      ! hand-off state when it is asked for. Its sweeps are the solver's and
      ! are tagged as such by the line above.
      call probe_eliminated_composition_is_seed_independent(Y, f_sp, u)
      call trace_composition_elimination(Y, f_sp, u)
      call scan_composition_elimination_basin(Y, f_sp, u)
      ! Which root of the front each admissible seed reaches, and what the
      ! spread of the closure costs the row maxima of the state, when that
      ! is asked for.
      call branch_and_closure_of_the_seed_family(Y, f_sp, u)
      ! Where the root of a carrier row sits relative to the lower bound of
      ! its own unknown, when that is asked for.
      call scan_carrier_row_across_its_lower_bound(Y, f_sp, u, 0)
      ! Seed from the adopted composition, write the swept one into the work
      ! array, then ADOPT it. The copy is an adoption and not a seeding
      ! convention: the seed is named in the call.
      call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0, n_part=npart0,  &
                         rows_judged=rows_j, resid_relnorm_judged=rnorm_j)
      f_sp = f_sp_j
      ! THE CARRIER UNKNOWNS BECOME LOGARITHMS HERE. The evaluation above is
      ! what makes the carrier operator freeze the element budget the floor
      ! is formed from, so this is the earliest point at which the space can
      ! be armed; F belongs to the state and stands for the rewritten vector
      ! as it stands for the packed one.
      call arm_carrier_log_unknown(Y)
      ! AND THE BOX THEY LIVE IN IS FORMED, from this state and from the
      ! element budget the evaluation above made the carrier operator
      ! freeze. It is re-formed at the top of every outer iteration, from
      ! that iteration's own iterate.
      call freeze_species_unknown_box(Y)
      call report_species_unknown_box
      ! The passes the elimination needed from a composition that was NOT
      ! this state's own: the marching hand-off. That is the number of terms
      ! of the Picard series the Newton model has to carry (see
      ! n_eq_sweeps_model).
      n_eq_sweeps_model = n_eq_sweeps_last
      n_no_chem_root_state = n_no_chem_root_last
      carrier_relnorm_state = carrier_relnorm_last
      carrier_cellmax_state = carrier_cellmax_last
      carrier_cell_state    = carrier_cell_worst
      call keep_background_of_adopted_state
      ! How much of each species row survives the cancellation of its own
      ! terms, where the rounding of that sum puts a floor under the row,
      ! and how the local chemical decay of the molecular carriers stands
      ! against the transport of the same cell, when those are asked for.
      ! It reads the rows the evaluation above formed, with nothing in
      ! between.
      call row_cancellation_and_the_rounding_floor(f_sp, u)
      call resid_relnorm(F, u, rc, rnorm)
      call cell_state_scales(Y, D)
      if (nspec_row .gt. 0) call cell_row_scales(Y, D, Drow, u)
      ! AND THE SHARED ELEMENT BUDGET AS ROWS OF THE FIRST STEP. They are
      ! formed after the scales because the projection acts in the
      ! coordinates the linear algebra works in, and after the iterate's own
      ! evaluation because the budget is a product of that evaluation
      ! (freeze_element_constraint_rows).
      if (nspec_row .gt. 0)                                              &
         call freeze_element_constraint_rows(Y, f_sp, D)
      call report_element_constraint_rows
      ! HOW MUCH OF THE INNER CLOSURE'S ERROR THE BASE ROW CARRIES, measured
      ! on this iterate, once. It sets the tolerance the sweep is asked for
      ! from here on, so that outer_residual_closure_target is a bound on
      ! what the outer residual keeps of that error rather than an
      ! assumption that the boundary amplifies nothing.
      if (nspec_row .gt. 0)                                              &
         call measure_the_closure_map_at_the_base(Y, f_sp, Drow)
      ! Line-search merit: the scaled 2-norm over the WHOLE domain 1..N. It has
      ! covered every physical cell since it was written, and it stays that
      ! way: restricting it to the wind -- a volume-weighted RMS of F/D over
      ! [j_min:N] -- was tried on 2026-08-11 and again on 2026-08-13, and both
      ! times it broke solves that converge with the whole-domain norm:
      ! photo_deep_secion_cont goes info = 0 at ||R|| = 7.28e-4 in 4 outer
      ! iterations against info = 2 at 3.25e-3 in 15, and examples/15
      ! (HD 209458 b, molecular) goes info = 0 at 9.54e-6 in 272 against
      ! info = 2 at 1.17e-4 in 30. The Newton step is computed from the full
      ! residual over ALL rows, so a merit that ignores most of those rows
      ! rejects the steps that step actually takes.
      !
      ! The merit and rnorm now cover the same cells, but they remain
      ! DIFFERENT FUNCTIONALS: a 2-norm of F/D, D being the cell's own signal
      ! scale rho(|v|+c_s) in the momentum row (cell_state_scales), against a
      ! maximum over rows and over the two regions of volume-weighted ratios
      ! taken on the scale each region's physics sets (residual_row_scale).
      ! So the iterate with the smallest ||R|| still need not be the last one
      ! -- handled at the end of this routine, not by changing the merit
      ! (docs/newton_scaling_and_base_wall.md section 4,
      ! docs/hd209_metal_stagnation.md).
      if (nspec_row .gt. 0) then
         f2 = sqrt(sum((F/Drow)**2))  ! merit in the SCALED space
      else
         f2 = sqrt(sum((F/D)**2))     ! merit in the SCALED space
      endif
      f2hist = f2                  ! non-monotone line-search memory
      ! The judged distance of the starting state, and the ledger's first
      ! entry. Every distance in this routine comes from the rows
      ! eval_residual measured on the state it evaluated, so the ledger
      ! compares one construction with itself.
      dj = distance_from_certification(rows_j)
      jhist = dj;  dj_best = dj
      rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp;  n_no_descent = 0
      n_since_best = 0;  monotone_search = .false.
      carrier_relnorm_best = carrier_relnorm_state
      carrier_cellmax_best = carrier_cellmax_state
      carrier_cell_best    = carrier_cell_state
      n_no_chem_root_best  = n_no_chem_root_state
      ! The best-iterate RANKING is deliberately not given the chemistry
      ! count: it chooses which state to hand back, and among states that
      ! all rest on a rootless cell the one closest to balance is still the
      ! one to return. The certification of the state handed back is made
      ! once, below, and it does read the count.
      gates_best = steady_gates_met(rnorm, u, resid_tol, fspread_best,   &
                                    nspec_row .gt. 0, carrier_relnorm_state)
      gates_now  = gates_best;  fspread = fspread_best
      call keep_background_of_best_iterate
      n_nonmonotone_accepts = 0;  n_judged_excursions = 0
      write(*,'(A,ES11.3,A,ES10.2)') ' (JFNK) start ||R||=',rnorm,        &
           '  ||Fs||2=',f2
      if (jfnk_maxit_wanted .gt. 0) then
         write(*,'(A,I0,A,I0,A)') ' (JFNK) outer iteration cap ',         &
              maxit_used, ', from EXHALE_JFNK_MAXIT in place of the ',    &
              maxit, ' the caller asked for'
      else
         write(*,'(A,I0,A)') ' (JFNK) outer iteration cap ', maxit_used,  &
              ', from the caller'
      endif
      call write_resid_below_escape(' (JFNK)', F, u)
      ! THE JACOBIAN ACTION OF THE SYSTEM THIS SOLVE ACTUALLY CARRIES,
      ! asserted on the state it starts from (EXHALE_SPECIES_JAC_TEST=1).
      ! The self-test at the top of the run tests the three-unknown system
      ! only, because the registry is empty until this routine is entered;
      ! with species rows the band layout, the coloring stride and the
      ! species columns are exactly what a wrong registry would get wrong,
      ! and an O(1) mismatch is what a wrong layout gives.
      call species_row_jacobian_action_check(Y, f_sp)
      ! AND WHERE THE COLUMNS OF THAT SYSTEM REACH, asserted on the same
      ! state (EXHALE_JAC_COLUMN_REACH=1): the band is the preconditioner and
      ! the damped Gauss-Newton model, so what it cannot hold is what those
      ! two do not see.
      call jacobian_column_reach_beyond_the_band(Y, f_sp)

      do iter = 1, maxit_used
         ! THE STOP TEST OF THE ITERATION, WHICH IS THE ACCEPTANCE TEST OF
         ! THE STATE (the same statement the PTC route makes at its loop
         ! top): the acceptance gate -- residual, flux spread, chemistry
         ! and the carrier norm -- AND every row the certification judges
         ! the system this solve carries by, each against its own
         ! tolerance. The gate alone is not the acceptance: it reads one
         ! number against the run's "Resid tol" while the certification
         ! reads the mass row against 3e-12, so a solve stopping on the
         ! gate stops at states the return then refuses (D3,
         ! docs/solver_partition_experiment_20260911.md section 7.3). When
         ! the gate is met and a certified row is not, the state is not a
         ! solution and the iteration goes on.
         !
         ! THE ROWS ARE MEASURED ONLY WHERE THEY CAN DECIDE, i.e. behind
         ! the gate: the certification evaluates the transported species
         ! balances as well, and an evaluation of it at every iteration
         ! would be paid for at every iteration for a test that can only
         ! change the outcome once the gate is already met.
         if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nspec_row .gt. 0, carrier_relnorm_state,  &
                          n_no_chem_root_state)) then
            call certified_rows_of_the_state(F, u, f_sp, resid_tol,       &
                     n_no_chem_root_state, rows_ok, i_row, cert_out,      &
                     n_rows_at_stop)
            if (rows_ok) then
               info = 0
               write(*,'(A,I0,A)') ' (JFNK) loop-top stop: the'//         &
                    ' acceptance gate is met and every one of the ',      &
                    n_rows_at_stop, ' certified row(s) of the system'//   &
                    ' this solve carries is within its own tolerance'
               exit
            endif
            n_gate_met_row_refused = n_gate_met_row_refused + 1
            if (i_row .ne. i_row_refused_last) then
               call write_gate_met_while_a_row_refuses(' (JFNK)',         &
                        cert_out, i_row, iter)
               i_row_refused_last = i_row
            endif
         endif
         ! Where the root of a carrier row sits relative to the bound of its
         ! own unknown, at THIS iterate, when that is asked for.
         call scan_carrier_row_across_its_lower_bound(Y, f_sp, u, iter)
         idtau = 1.0d0/dtau
         ! WHETHER THIS OUTER ITERATION IS ONE THE ELEMENT-ROW DIAGNOSTIC
         ! SPEAKS AT (elem_diag_on): the first two, where the model is still
         ! being built from an admitted state, and one inside whatever the
         ! solve settles into.
         elem_diag_here = elem_diag_on .and. (iter .le. 2 .or.            &
                                              iter .eq. elem_diag_third)
         ! AND WHETHER THIS ONE IS AN ITERATION THE ADDITIVITY HOOK
         ! SPEAKS AT: the first, one well inside whatever the solve
         ! settles into, and the last of a capped run, so that a defect
         ! measured at the entry state can be told from the one at the
         ! iterate the solve stops on.
         jv_additivity_here = jv_additivity_on .and.                     &
                              (iter .eq. 1 .or. iter .eq. 20 .or.        &
                               iter .eq. maxit_used)
         ! AND WHETHER THIS ONE IS THE ITERATION THE JUMP SCAN SPEAKS AT:
         ! the first, the state at which the floor was measured (N31, N32).
         if (iter .eq. 1) resid_jump_scan_here = resid_jump_scan_on
         ! AND WHETHER EVERY TRIAL OF IT STATES ITS OWN MODEL: from the
         ! selected iteration onward, so that consecutive iterations at one
         ! state can be compared step against step.
         elem_diag_every_trial = elem_diag_on .and.                       &
                                 (iter .ge. elem_diag_third)
         ! HOW MANY CYCLES THE KRYLOV SOLVE OF THIS ITERATION MAY SPEND.
         gm_restart_cycles = 1
         if (iter .ge. gm_restart_cycles_from)                            &
            gm_restart_cycles = gm_restart_cycles_wanted

         ! FREEZE the WENO3 weights at the current iterate: one mode-1
         ! residual evaluation stores the smoothness factors and yields the
         ! frozen-weights F; the inner evaluations that build the Newton model
         ! (Jacobian probes, J*v) then reuse them (mode 2), so the inner
         ! problem excludes the strongly nonlinear weight response (the
         ! standard lagged-weights remedy for FV steady solves). The line
         ! search does NOT use them -- see there.
         weno_mode = 1
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
         call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0,             &
                            rows_judged=rows_j,                          &
                            resid_relnorm_judged=rnorm_j)
         dj = distance_from_certification(rows_j)
         ! The reference the trial admissibility test compares against: how
         ! many cells THIS iterate's own chemistry left without a root.
         ! Refreshed
         ! here, at the top of every outer iteration, from the iterate itself.
         n_no_chem_root_state = n_no_chem_root_last
         carrier_relnorm_state = carrier_relnorm_last
         carrier_cellmax_state = carrier_cellmax_last
         carrier_cell_state    = carrier_cell_worst
         call keep_background_of_adopted_state
         ! Y at the top of an outer iteration IS the state the solve holds,
         ! so this covers every state it adopts. Reported here as well as at
         ! the end of the solve, so that the statement survives a run that
         ! was stopped rather than completed.
         !
         ! IT IS TAKEN AFTER THE ITERATE'S OWN EVALUATION, and that is what
         ! makes the element budget the number to compare against: the
         ! budget is a product of the carrier assembly, so before that
         ! evaluation carrier_headroom still holds the budget of the
         ! PREVIOUS iterate and the ratio would be a carrier of one state
         ! over the budget of another. MEASURED with the call at the loop
         ! top: `mol_carrier` with the budget a face reported 67 adopted
         ! carriers above their budget and a worst ratio of 1.058 at
         ! gm_m = 80, all of it the budget moving with the state. The count
         ! of species unknowns at a bound is the same statement: written by
         ! every residual evaluation, before this it described whatever
         ! trial the previous iteration evaluated last.
         call note_adopted_species_unknowns(Y)
         if (nspec_row .gt. 0 .and. carrier_rows_registered()) then
            write(*,'(A,ES11.3,A,I0,A,ES10.2,A,I0,A,ES10.2,A,I0)')        &
                 ' (JFNK) species unknowns:'//                            &
                 ' smallest adopted ', species_unknown_min_adopted,        &
                 ', unknowns on a bound ', n_species_at_lower_bound,       &
                 ', smallest carrier over its own floor ',                 &
                 carrier_over_its_floor_min,                               &
                 ', carriers at that floor ', n_carrier_at_its_floor,      &
                 ', largest carrier over its own element budget ',         &
                 carrier_over_its_budget_max,                              &
                 ', carriers above that budget ', n_carrier_above_its_budget
         else if (nspec_row .gt. 0) then
            write(*,'(A,ES11.3,A,I0)') ' (JFNK) species unknowns:'//       &
                 ' smallest adopted ', species_unknown_min_adopted,        &
                 ', unknowns on a bound ', n_species_at_lower_bound
         endif
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
         weno_mode = 2
         ! DIAGONAL SCALING FOR THIS OUTER ITERATION, and it is taken AFTER
         ! the evaluation of the iterate. A species row's scale is not a
         ! function of Y alone: cell_row_scales reads the elemental row
         ! scales element_transport_residual formed and the carrier scales
         ! carrier_steady_residual formed, which are the products of the
         ! LAST residual evaluation. The evaluation above is the iterate's
         ! own, so those products describe Y; taken before it, the scaling
         ! of this iteration is that of whatever state the previous
         ! iteration evaluated last -- a refused trial, or a point on a
         ! trust-region ray.
         call cell_state_scales(Y, D)
         if (nspec_row .gt. 0) call cell_row_scales(Y, D, Drow, u)
         ! AND THE FACES OF THE SPECIES BOX FOR THIS STEP, taken after the
         ! iterate's own evaluation for the same reason the scales are: the
         ! element budget one of the faces is is a product of that
         ! evaluation (freeze_species_unknown_box).
         if (nspec_row .gt. 0) call freeze_species_unknown_box(Y)
         if (nspec_row .gt. 0)                                            &
            call freeze_element_constraint_rows(Y, f_sp, D)
         ! AND THE CLOSURE MAP AT THE BASE AGAIN, where a restart asked for
         ! it. Taken here because the map is measured on the row scales of
         ! the iterate it describes, and those were just formed from this
         ! iterate's own evaluation.
         if (closure_map_due .and. nspec_row .gt. 0) then
            call measure_the_closure_map_at_the_base(Y, f_sp, Drow)
            closure_map_due = .false.
         endif
         if (nspec_row .gt. 0) then
            f2 = sqrt(sum((F/Drow)**2)) ! TRUE merit of the current iterate
         else
            f2 = sqrt(sum((F/D)**2))    ! TRUE merit of the current iterate
         endif
         ! WHICH ROW KIND HOLDS THE MEASURE AND WHICH HOLDS THE MERIT. Said
         ! here because the row scales were just formed from the iterate's
         ! own evaluation and nothing has been evaluated since, so the
         ! elemental row measures this reads still describe the iterate.
         if ((elem_diag_here .or. elem_diag_every_trial) .and.            &
             nspec_row .gt. 0)                                            &
            call write_row_kind_shares(F, u, Drow, ' (JFNK)')
         ! THE BASE POINT THE NEWTON MODEL IS DIFFERENCED ABOUT: the same
         ! state, the same seed, the same reconstruction mode and the same
         ! number of composition passes as the probes. Differencing two
         ! different maps is not a derivative of either.
         call eval_residual(Y, f_sp, f_sp_j, F_jac, heat0, cool0,          &
                            n_eq_sweeps_fixed=n_eq_sweeps_model)

         ! Grippo non-monotone reference: the worst true merit of the last 5
         ! outer iterates, the current one included. Recording it here, and not
         ! only when a step is accepted, is what keeps the reference from
         ! lagging behind the state the line search actually starts from.
         ! WHAT THE STEP INTO THIS STATE WAS ALLOWED TO DO, and what it did:
         ! the state's distance from certification against the window the
         ! acceptance held the step to. Printed before the window moves, so
         ! the two numbers on one line are the bound and the outcome, and
         ! `distance <= window` from the second outer iteration on is the
         ! invariant the acceptance installs.
         write(*,'(A,ES11.3,A,ES11.3,A,ES10.3,A,ES10.3,A,ES10.3,'//       &
                 'A,ES10.3)')                                             &
              ' (JFNK) judged rows: distance', dj, '  window',            &
              maxval(jhist), '   ||R||', rnorm_j, '  hydrodynamic',       &
              rows_j(1), '  element',                                     &
              rows_j(2), '  carrier', rows_j(3)
         f2hist = (/ f2hist(2:5), f2 /)
         ! And the same window on the quantity the state is judged by.
         jhist  = (/ jhist(2:5),  dj /)
         ! A RESTART COLLAPSES BOTH WINDOWS ONTO THIS STATE. Each is a bound
         ! about the point the acceptance compares against, and after a
         ! restart that point is this iterate and not the worst of four that
         ! are behind it; the monotone flag goes back to what entering a
         ! solve afresh leaves it as.
         if (acceptance_window_due) then
            f2hist = f2
            jhist  = dj
            monotone_search       = .false.
            acceptance_window_due = .false.
         endif

         ! Banded preconditioner of the SCALED system:
         ! M = I/dtau + D^-1 J_banded D, factored.
         call build_banded_jac_full(Y, f_sp, ab, n_jac_unresolved)
         n_jac_unresolved_tot = n_jac_unresolved_tot + n_jac_unresolved
         ! The pseudo-transient term is the identity of the UNSCALED
         ! system, so under a two-sided scaling it is diag(D/Drow), not the
         ! identity -- the same wherever the two scalings are, which is
         ! every row but the carrier's. The three-unknown branch is kept as
         ! the original text so that it stays byte-identical.
         if (nspec_row .gt. 0) then
            do jc = 1, neq
               ilo = max(1,   jc - ku_jac)
               ihi = min(neq, jc + kl_jac)
               do irow = ilo, ihi
                  ab(kl_jac+ku_jac+1 + irow - jc, jc) =                  &
                       ab(kl_jac+ku_jac+1 + irow - jc, jc)*D(jc)/Drow(irow)
               enddo
            enddo
            ! AND THE ONE FREEDOM LEFT IN THE SCALED OPERATOR IS TAKEN
            ! HERE: with the rows fixed by the gate, the columns are
            ! equilibrated and D is carried with them, so that everything
            ! downstream -- the factorization, the Krylov space, the
            ! trust-region ball and dY = D*s -- works in one set of
            ! coordinates. The merit is not touched: it holds no D.
            if (column_equilibration_on)                                 &
               call equilibrate_the_scaled_columns(ab, D)
            ! AND THE ROW SCALING OF THE LINEAR MODEL, which decision 20 a
            ! left free: formed from ab and applied to the band that is
            ! factored and, inside the cycle, to the right-hand side and
            ! the action. ab itself is untouched, so the merit's gradient
            ! and the trust-region model stay on the certification scales
            ! (unit_infinity_norm_row_scaling_of_the_band).
            model_rows_equilibrated = .false.
            if (model_row_equilibration_on)                              &
               call unit_infinity_norm_row_scaling_of_the_band(ab)
            abf = ab
            do jc = 1, neq
               abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc)       &
                                        + idtau*D(jc)/Drow(jc)
            enddo
            if (model_rows_equilibrated) then
               do jc = 1, neq
                  ilo = max(1,   jc - ku_jac)
                  ihi = min(neq, jc + kl_jac)
                  do irow = ilo, ihi
                     abf(kl_jac+ku_jac+1 + irow - jc, jc) =              &
                          abf(kl_jac+ku_jac+1 + irow - jc, jc)           &
                          *model_row_equilibration(irow)
                  enddo
               enddo
            endif
            ! AND THE BAND THAT IS FACTORIZED IS EQUILIBRATED, the model
            ! above it is not. The two rescalings are undone inside the
            ! preconditioner solve, so the Krylov cycle sees the same
            ! approximate inverse of the same operator, computed from a
            ! factorization that is not a division by round-off.
            preconditioner_equilibrated = .false.
            if (preconditioner_equilibration_on)                         &
               call equilibrate_the_preconditioner(abf)
         else
            do jc = 1, neq
               ilo = max(1,   jc - ku_jac)
               ihi = min(neq, jc + kl_jac)
               do irow = ilo, ihi
                  ab(kl_jac+ku_jac+1 + irow - jc, jc) =                  &
                       ab(kl_jac+ku_jac+1 + irow - jc, jc)*D(jc)/D(irow)
               enddo
            enddo
            abf = ab
            do jc = 1, neq
               abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + idtau
            enddo
            preconditioner_equilibrated = .false.
            model_rows_equilibrated = .false.
            if (model_row_equilibration_on) then
               call unit_infinity_norm_row_scaling_of_the_band(ab)
               do jc = 1, neq
                  ilo = max(1,   jc - ku_jac)
                  ihi = min(neq, jc + kl_jac)
                  do irow = ilo, ihi
                     abf(kl_jac+ku_jac+1 + irow - jc, jc) =              &
                          abf(kl_jac+ku_jac+1 + irow - jc, jc)           &
                          *model_row_equilibration(irow)
                  enddo
               enddo
            endif
         endif
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         call banded_model_conditioning(ab, abf, D, Drow, ' (JFNK)')
         if (elem_diag_here)                                              &
            call write_scales_columns_and_pivots(D, Drow, ab, abf, ipiv,  &
                      idtau, lpinfo, n_jac_unresolved, ' (JFNK)')
         if (lpinfo .ne. 0) then
            dtau = max(dtau*0.25d0, dtau0);  cycle
         endif

         ! The window on the judged quantity, formed here because both step
         ! controls below are bounded by it (see jhist).
         if (monotone_search) then
            jref = dj
         else
            jref = maxval(jhist)
         endif

         ! ---- D2: the scaled trust region, four-unknown branch only ----
         if (use_tr) then
            ! THE RADIUS IS RE-INITIALIZED WHENEVER IT IS NOT A POSITIVE
            ! FINITE LENGTH. A hundredth of a zero Krylov step is zero, and
            ! a ball of radius zero contains a single point: an iteration
            ! handed one can report nothing but a zero step, so a sign test
            ! leaves it absorbing (R1).
            if (.not. (tr_delta .gt. 0.0d0) .or.                          &
                .not. finite_real(tr_delta)) then
               ! THE BOUNDS THIS ITERATE HOLDS ARE FIXED BEFORE THE CYCLE
               ! SAMPLES ANYTHING, for the same reason trust_region_step
               ! fixes them before its own cycle: a direction with a
               ! component pointing out of the species box has no admissible
               ! sample at any step length, so a cycle run with nothing held
               ! returns no direction at all and there is nothing to measure
               ! a radius on.
               r0_init = F/Drow
               call band_matvec_transpose(ab, r0_init, g_init)
               call fix_active_species_bounds(Y, g_init, n_active_init)
               call pgmres(Y, F_jac, f_sp, D, Drow, abf, ipiv, idtau,      &
                           -r0_init,                                      &
                           dZ, gm_m_used, 1.0d-1, gmit,                   &
                           gm_outcome, tr_gm_resid, gm_snorm_last)
               call release_active_species_bounds
               call initial_trust_region_radius(gm_outcome, dZ, g_init,   &
                                                tr_delta, tr_dmax,        &
                                                tr_init_status)
               call tr_radius_status_text(tr_init_status, tr_why_txt)
               call gm_outcome_text(gm_outcome, gm_why_txt)
               write(*,'(A,ES11.3,A,A)') ' (JFNK) [TR] initial radius',   &
                    tr_delta, ': ', trim(tr_why_txt)
               write(*,'(A,ES11.3,A,ES11.3,A,A)')                          &
                    ' (JFNK) [TR]      Krylov step', gm_snorm_last,        &
                    '   projected gradient', scaled_two_norm(g_init),      &
                    '; ', trim(gm_why_txt)
               write(*,'(A,ES9.2,A,I0)') ' (JFNK) [TR] Krylov tolerance'// &
                    ' of the dogleg Newton leg:', tr_gm_rtol,               &
                    '   subspace size: ', gm_m_used
               if (tr_init_status .eq. tr_radius_no_feasible_descent) then
                  ! NOTHING TO STEP ALONG, and repeating zero-radius
                  ! iterations is not recovery. The solve stops and says
                  ! which failure it was, with the linear outcome beside it.
                  write(*,'(A,A)') ' (JFNK) [TR] stopping: ',              &
                       trim(tr_why_txt)
                  info = 3
                  exit
               endif
            endif
            ! IS THE MERIT THE SAME FUNCTION IT WAS LAST ITERATION? Drow is
            ! rebuilt from the state at the top of every outer iteration, so
            ! ||F/Drow|| is a different functional each time. A trust region
            ! compares a predicted with an actual reduction of ONE function;
            ! if the function itself moves between iterations the comparison
            ! is only valid inside an iteration. This prints the size of that
            ! move: the merit of the state just adopted, measured with the
            ! scaling it was adopted under, against the same state measured
            ! with the scaling of the iteration that follows.
            if (tr_f2_end .gt. 0.0d0)                                      &
               write(*,'(A,ES11.3,A,ES11.3,A,F9.3)')                       &
                    ' (JFNK) [TR]      merit across the iteration boundary:'//&
                    ' old scaling', tr_f2_end, '   new scaling', f2,        &
                    '   factor', f2/max(tr_f2_end,1.0d-99)
            call trust_region_step(Y, F, F_jac, f_sp, D, Drow, ab, abf,    &
                                   ipiv,                                   &
                                   idtau, gm_m_used, tr_gm_rtol, f2,       &
                                   tr_delta,                               &
                                   dY, f2_try, tr_pred, tr_actual,         &
                                   tr_ratio, tr_snorm, tr_cuts, gmit,      &
                                   gm_outcome, tr_gm_resid,                &
                                   tr_model_ok, tr_accept,                 &
                                   Ytry, Ftry, f_sp_j, tr_f2ref, tr_refuse,&
                                   rows_try, tr_cauchy_ascends)
            n_tr_exit(tr_refuse) = n_tr_exit(tr_refuse) + 1
            if (tr_cauchy_ascends) n_tr_cauchy_dropped =                  &
                                   n_tr_cauchy_dropped + 1
            ! THE CEILING IS NOT A PROPERTY OF THE FIRST STEP. An accepted
            ! step whose actual reduction reached three quarters of what the
            ! model promised is a model good over the whole radius, and the
            ! region is doubled for the next iteration; a ceiling fixed at a
            ! hundredfold of the FIRST radius would then stop the region
            ! growing and the solve would carry short steps for the rest of
            ! its life. The ceiling follows the radius by the same factor,
            ! up to the one length that is absolute
            ! (trust_region_radius_absolute_cap). The feasible set still
            ! bounds the STEP, through the projection onto the species box
            ! and the admissibility cuts; the cap bounds the NUMBER that
            ! claims to.
            if (tr_accept .and. tr_ratio .gt. 0.75d0)                     &
               call grow_the_trust_region_ceiling(tr_dmax)
            if (tr_dmax .gt. 0.0d0) tr_delta = min(tr_delta, tr_dmax)
            n_gm_outcome(gm_outcome) = n_gm_outcome(gm_outcome) + 1
            ok  = tr_accept
            ! The judged rows of the step the region accepted, for the
            ! best-iterate ledger below.
            dj_try = dj
            if (tr_accept) dj_try = distance_from_certification(rows_try)
            lam = 1.0d0
            tr_f2_end = -1.0d0
            if (tr_accept) then
               tr_f2_end   = f2_try
               n_tr_accept = n_tr_accept + 1
            else
               n_tr_reject = n_tr_reject + 1
            endif
            if (.not. tr_model_ok) n_tr_model_bad = n_tr_model_bad + 1
            ! One line per iteration, accepted or not: radius, predicted and
            ! actual reduction, their ratio, and the step actually taken.
            write(*,'(A,I4,A,ES10.3,A,ES11.3,A,ES11.3,A,F9.3,A,ES10.3,'//  &
                    'A,I2,A,L1,A,L1)')                                     &
                 ' (JFNK) [TR] it', iter, '  delta=', tr_delta,            &
                 '  pred=', tr_pred, '  actual=', tr_actual,               &
                 '  ratio=', tr_ratio, '  ||s||=', tr_snorm,               &
                 '  cuts=', tr_cuts, '  model_ok=', tr_model_ok,           &
                 '  accepted=', tr_accept
            if (tr_refuse .ne. tr_step_taken) then
               call tr_exit_text(tr_refuse, tr_why_txt)
               write(*,'(A,A)') ' (JFNK) [TR]      no step: ',             &
                    trim(tr_why_txt)
               ! AND WHICH SCREEN REFUSED THE SAMPLE. The exits that end on
               ! "no sample" name the outcome; this names the cause, from
               ! the record the refusing evaluation left. Nothing between
               ! trust_region_step's return and here evaluates a residual,
               ! so the record is still that refusal's.
               call write_last_evaluation_refusal(' (JFNK) [TR]')
            endif
            write(*,'(A,ES11.3,A,ES11.3,A,F8.3)')                          &
                 ' (JFNK) [TR]      merit of the iterate: frozen-mode',    &
                 f2, '   trial-mode', tr_f2ref, '   ratio', f2/max(tr_f2ref,1.0d-99)
            ! WHAT THE KRYLOV CYCLE OF THE NEWTON LEG REACHED AGAINST WHAT
            ! IT WAS ASKED FOR, and how much subspace it was given. A cycle
            ! that used its whole budget and is still above the tolerance
            ! returned the minimizer over a subspace and not the model's
            ! own, and the dogleg's monotonicity rests on the second
            ! (section 159.3); the two numbers side by side are what say
            ! which of the two happened.
            call gm_outcome_text(gm_outcome, gm_why_txt)
            write(*,'(A,I0,A,I0,A,ES10.3,A,ES9.2,A,A)')                    &
                 ' (JFNK) [TR]      Krylov leg: products ', gmit,          &
                 ' of ', gm_m_used, ', relative residual reached ',        &
                 tr_gm_resid, ' against ', tr_gm_rtol,                     &
                 '; ', trim(gm_why_txt)
         else
         ! GMRES solve of the scaled PTC system; unscale the step.
         if (nspec_row .gt. 0) then
            call pgmres(Y, F_jac, f_sp, D, Drow, abf, ipiv, idtau, -F/Drow, &
                        dZ, gm_m_used, 1.0d-1, gmit, gm_outcome,           &
                        tr_gm_resid, gm_snorm_last)
         else
         call pgmres(Y, F_jac, f_sp, D, Drow, abf, ipiv, idtau, -F/D, dZ, &
                     gm_m_used, 1.0d-1, gmit, gm_outcome, tr_gm_resid,    &
                     gm_snorm_last)
         endif
         n_gm_outcome(gm_outcome) = n_gm_outcome(gm_outcome) + 1
         dY = D*dZ
         endif

         ! Scaled ||F/D||_2 backtracking line search with positivity, tested
         ! against the non-monotone reference above.
         !
         ! The trial states are evaluated with weno_mode = 0, i.e. with the
         ! smoothness weights recomputed at the trial state, so what decides
         ! acceptance is the residual the solve is driving to zero. Evaluating
         ! the trials with the weights frozen at Y instead -- the quantity the
         ! inner Newton model minimizes -- is a different measure: on the
         ! HD 189733 b hand-off state the two differ by a median factor 2.2 at
         ! lam = 1, and 7 of the 14 steps the frozen test accepted RAISED the
         ! true residual, by up to 9x (docs/newton_scaling_and_base_wall.md
         ! §10). Mode 0 leaves the stored weights alone, so the frozen model of
         ! this outer iteration survives the search.
         ! Armijo after a stagnation restart: descent against the merit of
         ! the state the search starts from, not against the worst of the
         ! last five. Until the restart happens this is the original line.
         if (monotone_search) then
            f2ref = f2
            jref  = dj
         else
            f2ref = maxval(f2hist)
            jref  = maxval(jhist)
         endif
         if (.not. use_tr) then
         lam = 1.0d0;  ok = .false.;  f2_try = huge(1.0d0)
         weno_mode = 0
         ! A cycle that could sample NO Krylov direction returns the zero
         ! step. The zero step must not go through the line search: the
         ! non-monotone test compares the trial against the WORST merit of
         ! the last five iterates, so Ytry = Y passes it whenever the current
         ! iterate is not that worst one, and the iteration would be recorded
         ! as a descent step that moves nothing -- resetting the no-descent
         ! counter and letting the solve spin to maxit instead of aborting.
         ! There is no Newton model this iteration; that is a no-descent
         ! iteration, and it is counted as one.
         if (gmit .le. 0) then
            lam = 0.0d0
            write(*,'(A,I4,A)') ' (JFNK) it',iter,'  no Krylov direction'// &
                 ' could be sampled; no Newton model this iteration'
         endif
         do ls = 1, 20
            if (gmit .le. 0) exit
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)   ! ghosts of utry are not set
            ! The carrier unknown gets the same treatment rho and p get: a
            ! trial that puts a negative H2 density into a cell is rejected
            ! before the residual is asked for it.
            carrier_ok = .true.
            if (nspec_row .gt. 0) then
               do jj = 1, N
               do isr = 1, nspec_row
                  if (Ytry(nvar_jac*(jj-1)+3+isr) .lt. 0.0d0)            &
                     carrier_ok = .false.
               enddo
               enddo
            endif
            if (carrier_ok .and.                                        &
                minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok, rows_judged=rows_try)
               dj_try = distance_from_certification(rows_try)
               if (nspec_row .gt. 0) then
                  f2_try = sqrt(sum((Ftry/Drow)**2))
               else
                  f2_try = sqrt(sum((Ftry/D)**2))
               endif
               ! A trial state the code cannot describe is rejected here,
               ! not left to the merit comparison: a non-finite merit fails
               ! the test by accident of IEEE ordering, and a trial whose
               ! chemistry broke in only a few cells can still produce a
               ! finite, smaller merit and be adopted (section 121).
               !
               ! HOW OFTEN THE STEP CONTROL AND THE JUDGEMENT DISAGREE,
               ! counted and reported and NOT acted on. Making this a
               ! condition of the step was tried and MEASURED to stop a
               ! solve that otherwise converges (report B5e section 3.3):
               ! the judged measure is a maximum over cells and over rows
               ! and is far less smooth than the merit, so the merit's own
               ! five-iterate window pins a solve that is legitimately
               ! trading one row against the others. The place the
               ! judgement belongs is the ledger that chooses the state to
               ! hand back, and that is where it now is.
               if (try_ok .and. dj_try .gt. jref .and.                   &
                   f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2ref)               &
                    n_judged_excursions = n_judged_excursions + 1
               if (try_ok .and.                                          &
                   f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2ref) then
                  ok = .true.
                  if (f2_try .ge. (1.0d0 - 1.0d-4*lam)*f2)              &
                       n_nonmonotone_accepts =                    &
                            n_nonmonotone_accepts + 1
                  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         ! WHEN BACKTRACKING FINDS NOTHING AND SHORTENING THE STEP CANNOT
         ! HELP. The search above only shortens the Newton/PTC direction. If
         ! that direction points uphill, every backtrack points uphill too --
         ! the excess merit falls linearly with lam and never changes sign --
         ! and with dtau already at its floor the next outer iteration would
         ! rebuild the same preconditioner, the same Krylov direction and the
         ! same search, and repeat this iteration bit for bit until the
         ! stagnation counter fires (11 identical repetitions, measured on the
         ! He/H = 0.3 hand-off state). What is wrong there is the DIRECTION,
         ! so a direction that is guaranteed downhill is built instead
         ! (docs/Update_EXHALE_stage1.md section 126). While dtau is still above its
         ! floor the reduction below does change the next iteration, and the
         ! PTC ramp is left to do its work.
         endif   ! .not. use_tr -- the trust region replaces the search
         ! The damped Gauss-Newton escape is what the trust region REPLACES:
         ! it looks for a descent direction after the fact, where the region
         ! controls the step by model agreement before the fact. Running both
         ! would make the reduction ratio meaningless, so it is skipped.
         lm_tried = .false.;  lm_found = .false.
         if ((.not. use_tr) .and. dtau .le. dtau0*(1.0d0 + 1.0d-12)) then
            ! Two ways for the Newton/PTC direction to leave the iterate where
            ! it is, both of them repeatable bit for bit at this dtau: the
            ! backtracking search rejects every step, or it accepts one that
            ! does not lower the true merit -- the non-monotone test compares
            ! against the worst of the last five iterates, so a sideways step
            ! passes it and the solve can circle for hundreds of iterations
            ! (measured: 487 such acceptances in 500 iterations on the
            ! He/H = 0.3 hand-off state).
            lm_tried = (.not. ok)
            if (ok) lm_tried = (f2_try .ge. f2)
         endif
         if (lm_tried) then
            call levenberg_marquardt_descent(Y, F, f_sp, D, Drow, ab, f2,&
                            dY_lm, f2_lm, mu_lm, gnorm_lm, lm_found)
            if (lm_found .and. f2_lm .ge. f2_try) lm_found = .false.
            if (lm_found) then
               ! BACKTRACK ALONG THE DAMPED GAUSS-NEWTON DIRECTION, in the
               ! coupled route only. The escape used to evaluate the full
               ! step and nothing else, which is enough while the only way
               ! it can fail is to ascend -- the direction is built to
               ! descend, so a shorter one descends too. With the carrier
               ! unknown a second failure mode appears and it is not about
               ! descent: the full step lowers the merit (2.28 -> 1.50,
               ! measured) but puts 73 cells' molecular equilibrium outside
               ! the physical simplex, so the trial is refused and the solve
               ! aborts with a descent direction in hand. Shortening the
               ! step is what that asks for. The 3-unknown route keeps the
               ! single evaluation it has always had, so no atomic solve
               ! moves.
               if (nspec_row .le. 0) then
               Ytry = Y + dY_lm
               weno_mode = 0
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok, rows_judged=rows_try)
               dj_try = distance_from_certification(rows_try)
               f2_try = sqrt(sum((Ftry/D)**2))
               ok  = try_ok
               lam = 1.0d0
               else
                  lam = 1.0d0
                  do ls = 1, 12
                     Ytry = Y + lam*dY_lm
                     carrier_ok = .true.
                     do jj = 1, N
                     do isr = 1, nspec_row
                        if (Ytry(nvar_jac*(jj-1)+3+isr) .lt. 0.0d0)      &
                           carrier_ok = .false.
                     enddo
                     enddo
                     if (carrier_ok) then
                        weno_mode = 0
                        call eval_residual(Ytry, f_sp, f_sp_j, Ftry,      &
                                           heat0,                        &
                                           cool0, admissible=try_ok,      &
                                           rows_judged=rows_try)
                        dj_try = distance_from_certification(rows_try)
                        f2_try = sqrt(sum((Ftry/Drow)**2))
                     else
                        try_ok = .false.;  f2_try = huge(1.0d0)
                        dj_try = huge(1.0d0)
                     endif
                     if (try_ok .and. f2_try .lt. f2) exit
                     lam = 0.5d0*lam
                  enddo
                  ok  = try_ok .and. f2_try .lt. f2
               endif
               write(*,'(A,I4,A,ES9.2,A,ES10.3,A,ES10.3)')               &
                    ' (JFNK) it',iter,'  damped Gauss-Newton escape,'//  &
                    ' mu=',mu_lm,'  ||Fs||2',f2,' ->',f2_try
            else if (.not. ok) then
               ! Nothing downhill along either direction. With the damped
               ! Gauss-Newton family exhausted this is a stationary point of
               ! the merit to within what the banded model can see, so the
               ! gradient norm is what says whether that is true.
               write(*,'(A,I4,A,ES10.3)') ' (JFNK) it',iter,'  no descent'// &
                    ' along the Newton or the damped Gauss-Newton'//     &
                    ' direction; ||grad merit||=', gnorm_lm
               call write_last_evaluation_refusal(' (JFNK)')
            endif
         endif

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            ! The last residual evaluation was this accepted trial, so the
            ! frozen background, and the count of cells whose composition is
            ! not a root of the network, now describe the state just adopted.
            n_no_chem_root_state = n_no_chem_root_last
            carrier_relnorm_state = carrier_relnorm_last
            carrier_cellmax_state = carrier_cellmax_last
            carrier_cell_state    = carrier_cell_worst
            ! And its elimination started from the PREVIOUS iterate's
            ! composition, so its pass count is the length of the Picard
            ! series the next Newton model must carry.
            n_eq_sweeps_model    = max(n_eq_sweeps_model, n_eq_sweeps_last)
            call keep_background_of_adopted_state
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            ! The judged distance of the state just adopted, from the rows
            ! the evaluation that accepted it measured.
            dj = dj_try
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            dtau = max(dtau*0.25d0, dtau0)
         endif

         ! Keep the best iterate seen so that a failed solve returns it
         ! rather than wherever it stopped, and never prefer one that fails
         ! the flux gate over one that passes it.
         !
         ! THE RANKING IS THE MEASURE THE RETURNED STATE IS JUDGED BY, which
         ! is not ||R||, AND THAT IS TRUE OF THE THREE HYDRODYNAMIC ROWS ON
         ! THEIR OWN. ||R|| is the largest of the three row measures; the
         ! state handed back is judged row by row, each against its own
         ! tolerance (mass 3e-12, momentum 1e-8, energy 1e-6) and, where the
         ! registry carries them, by every species row as well
         ! (stationary_rows_of_the_returned_state). So a ledger ranking on
         ! ||R|| can keep either an iterate whose element or carrier balance
         ! is the worse of two, or -- with no species row at all -- an
         ! iterate whose mass row is the worse of two while the energy row
         ! that holds ||R|| is the better, and hand it back to be refused on
         ! exactly that row. distance_from_certification is the same rows the
         ! certification reads, each over the tolerance it is judged against,
         ! and it is the ranking on every route.
         !
         ! THIS IS WHERE THE JUDGED MEASURE BELONGS AND THE STEP CONTROL IS
         ! NOT. Bounding the STEP by this quantity was implemented and
         ! measured to stall the element solve (report B5e section 3.3): the
         ! merit is a smooth 2-norm and can be globalized, the judged measure
         ! is a maximum over cells and over rows and is a stopping test.
         gates_now = steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nspec_row .gt. 0, carrier_relnorm_state)
         keep_this_iterate = (dj .lt. dj_best)
         if ((gates_now .and. .not. gates_best) .or.                      &
             ((gates_now .eqv. gates_best) .and. keep_this_iterate)) then
            dj_best = dj
            rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp
            fspread_best = fspread;  gates_best = gates_now
            call keep_background_of_best_iterate
            n_since_best = 0
            carrier_relnorm_best = carrier_relnorm_state
            carrier_cellmax_best = carrier_cellmax_state
            carrier_cell_best    = carrier_cell_state
            n_no_chem_root_best  = n_no_chem_root_state
         else
            n_since_best = n_since_best + 1
         endif

         ! P51: the best iterate has stopped moving. Go back to it and stop
         ! accepting sideways steps; if that buys nothing either, stop.
         if (n_since_best .ge. n_stall_best) then
            if (.not. monotone_search) then
               Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
               n_no_chem_root_state = n_no_chem_root_best
               carrier_relnorm_state = carrier_relnorm_best
               carrier_cellmax_state = carrier_cellmax_best
               carrier_cell_state    = carrier_cell_best
               call adopt_background_of_best_iterate
               call unpack_U(Y, u)
               call unpack_species_rows(Y, u, f_sp)
               call Apply_BC(u)
               monotone_search = .true.
               n_since_best    = 0
               ! THE MEASURE THE DETECTOR READ, beside the one the gate
               ! reports. With a species row registered the ledger ranks
               ! iterates on the judged distance from certification, not on
               ! ||R||, so a stop explained by ||R|| alone quotes a
               ! functional the decision was not taken in.
               write(*,'(A,I0,A,ES11.3,A,ES11.3,A,ES11.3,A)')            &
                    ' (JFNK) best iterate has not improved in ',          &
                    n_stall_best, ' iterations; restarting from it'//      &
                    ' (judged distance ', dj, ', best ', dj_best,         &
                    ', ||R||=', rnorm_best,                               &
                    ') with a monotone line search'
            else
               info = 2
               write(*,'(A,I0,A,ES11.3,A,ES11.3,A,ES11.3)')             &
                    ' (JFNK) STAGNATED: no improvement of the best'//     &
                    ' iterate in ', 2*n_stall_best, ' iterations,'//      &
                    ' monotone search included -- stopping at a judged'// &
                    ' distance ', dj, ', best ', dj_best, ', ||R||=',     &
                    rnorm_best
               exit
            endif
         endif

         ! Stagnation: the failure mode of the non-monotone search is that it
         ! stops finding ANY acceptable step and the iterate no longer moves.
         ! Count consecutive such iterations. (A watchdog on rnorm instead
         ! cannot work here: the solver minimizes the merit, not rnorm, and
         ! the opening pseudo-transient legitimately raises rnorm for ~30
         ! iterations while it repairs the sub-sonic region -- see
         ! docs/newton_scaling_and_base_wall.md.)
         if (ok) then
            n_no_descent = 0
         else
            n_no_descent = n_no_descent + 1
         endif
         if (lm_tried .and. .not. ok) then
            ! The escape above was reached, which means repeating this
            ! iteration would reproduce it exactly. There is nothing left to
            ! try, so the solve stops here instead of counting to
            ! n_no_descent_max identical iterations.
            info = 2
            write(*,'(A)') ' (JFNK) no descent direction exists for the'// &
                 ' banded model at this state -- aborting'
            exit
         endif
         if (n_no_descent .ge. n_no_descent_max) then
            info = 2
            write(*,'(A,I0,A)') ' (JFNK) line search found no descent '// &
                 'step in ', n_no_descent_max, ' consecutive iterations'// &
                 ' -- aborting'
            exit
         endif

         ! --- A RESTART OF THE TRUST-REGION STATE ---
         !
         ! What a solve carries beside its iterate is reset here and the
         ! iterate is not touched, except by the one reset that is a return
         ! to the best iterate seen. It is taken AFTER the stagnation block
         ! above, so that block keeps its priority within an iteration; a
         ! trigger threshold below n_stall_best does however reset the
         ! count that block reads, and the two are then one mechanism with
         ! the restart's period.
         !
         ! It is reached on the coupled route only: the three-unknown route
         ! never sets use_tr and carries no trust-region state to restart.
         if (use_tr) then
            tr_restart_due = trust_region_restart_is_due(iter, n_since_best)
            if (tr_restart_due) then
               if (tr_reset_wanted(tr_reset_back_to_best)) then
                  Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
                  n_no_chem_root_state = n_no_chem_root_best
                  carrier_relnorm_state = carrier_relnorm_best
                  carrier_cellmax_state = carrier_cellmax_best
                  carrier_cell_state    = carrier_cell_best
                  call adopt_background_of_best_iterate
                  call unpack_U(Y, u)
                  call unpack_species_rows(Y, u, f_sp)
                  call Apply_BC(u)
                  dj = dj_best
               endif
               call restart_the_trust_region_state(tr_delta, tr_dmax,     &
                        dtau, dtau0, n_since_best, n_no_descent,          &
                        closure_map_due, acceptance_window_due)
               write(*,'(A,I0,A,I0,A)') ' (JFNK) [TR] restart ',          &
                    n_tr_restart, ' of the trust-region state after'//    &
                    ' outer iteration ', iter, ':'
               do i_reset = 1, n_tr_reset_kind
                  if (.not. tr_reset_wanted(i_reset)) cycle
                  call tr_reset_name(i_reset, reset_txt)
                  write(*,'(A,A)') ' (JFNK) [TR]      reset ',            &
                       trim(reset_txt)
               enddo
            endif
         endif

         ! Locate the cell carrying the largest scaled residual |F/D|, i.e.
         ! the term the merit is actually dominated by, over the whole domain
         ! 1..N. (Normalizing by max_j|u(k,j)| instead reports the base cell
         ! almost always, because the base holds the global maximum of |rho v|
         ! while carrying no wind.)
         amx = -1.0d0;  jworst = 1;  kworst = 1
         do jj = 1, N
            do kk = 1, nvar_jac
               if (nspec_row .gt. 0) then
                  amx0 = abs(F(nvar_jac*(jj-1)+kk))                      &
                       / Drow(nvar_jac*(jj-1)+kk)
               else
                  amx0 = abs(F(nvar_jac*(jj-1)+kk))                      &
                       / D(nvar_jac*(jj-1)+kk)
               endif
               if (amx0 .gt. amx) then
                  amx = amx0
                  jworst = jj;  kworst = kk
               endif
            enddo
         enddo
         ! WHICH ROW CARRIES THAT MAXIMUM, BY NAME. A species row's slot
         ! number says the cell but not the element or the carrier, and an
         ! integer field cannot hold a slot of two digits at all -- the
         ! sodium row of cell 500 printed as k=* for fifty iterations
         ! (item N24). name_of_unknown is the one place the flat index is
         ! turned into a cell and a quantity.
         call name_of_unknown(nvar_jac*(jworst-1)+kworst,              &
                              what_worst_row)
         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES9.2,A,I3,A,F7.3,A,A)')      &
              ' (JFNK) it',iter,'  ||R||=',rnorm,'  ||Fs||2=',f2,        &
              '  lam=',lam,'  gm=',gmit,'  worst r=',r(jworst),          &
              '  worst row: ',trim(what_worst_row)
      enddo

      weno_mode = 0                 ! restore default reconstruction
      call set_ioniz_eq_sweep_state_kind(ieq_state_marching)

      ! What the solve could not describe. Both are properties of the
      ! iterates this solve visited, not of the returned state, and neither
      ! is an error by itself: a zeroed Jacobian color only degrades the
      ! preconditioner, and a truncated cycle only shortens the Krylov
      ! subspace. They are printed because a solve that ends info = 2 with
      ! either of them non-zero failed for a reason the iteration counts do
      ! not show.
      if (use_tr) then
         write(*,'(A,I0,A,I0,A,I0,A,ES10.3)')                             &
              ' (JFNK) [TR] steps accepted ', n_tr_accept, ', rejected ', &
              n_tr_reject, ', model refused ',                            &
              n_tr_model_bad, '; final radius', tr_delta
         if (n_tr_cauchy_dropped .gt. 0)                                  &
            write(*,'(A,I0,A)') ' (JFNK) [TR] the banded gradient'//      &
                 ' ascended the model at ', n_tr_cauchy_dropped,          &
                 ' iteration(s); the Cauchy leg was dropped there and'//  &
                 ' the dogleg reduced to the Krylov leg'
         ! WHICH REFUSAL, not how many. An inadmissible state on the ray and
         ! a model that over-predicts are different failures with different
         ! repairs, and the single counter above cannot separate them.
         do jj = 1, n_tr_exit_reason
            if (n_tr_exit(jj) .le. 0) cycle
            call tr_exit_text(jj, tr_why_txt)
            write(*,'(A,I4,A,A)') ' (JFNK) [TR] iterations that ended on ',&
                 n_tr_exit(jj), ': ', trim(tr_why_txt)
         enddo
      endif
      if (n_jac_unresolved_tot .gt. 0)                                   &
         write(*,'(A,I0,A)') ' (JFNK) probe states without an '//        &
              'admissible residual: ', n_jac_unresolved_tot,             &
              ' Jacobian color(s) zeroed'
      ! HOW THE LINEAR SOLVE ENDED, ITERATION BY ITERATION, BY NAME. A
      ! cycle that reached its tolerance and one that could not sample a
      ! single direction are both "one linear solve"; only the names say
      ! which of them this solve was made of.
      do jj = 0, n_gm_outcome_kind
         if (n_gm_outcome(jj) .le. 0) cycle
         call gm_outcome_text(jj, gm_why_txt)
         write(*,'(A,I4,A,A)') ' (JFNK) linear solves that ended on ',    &
              n_gm_outcome(jj), ': ', trim(gm_why_txt)
      enddo
      if (gm_measure_orthogonality)                                       &
         write(*,'(A,ES10.3)') ' (JFNK) worst loss of orthogonality of'// &
              ' an Arnoldi basis of this solve:', gm_orthogonality_loss_max

      ! Residual of the layer below the escape radius, resolved into its three
      ! rows and its worst cell -- the half of ||R|| that the single number
      ! reports only when it is the larger one. Printed HERE, before the
      ! best-iterate restore below, because this is the last point at which F
      ! and u are a consistent pair: the restore replaces Y (hence u) by the
      ! best iterate without re-evaluating F, and the line is not worth a
      ! further residual evaluation. So it refers to the LAST iterate visited,
      ! which is also the returned state whenever no restore happens.
      call write_resid_below_escape(' (JFNK)', F, u)

      ! Return the best iterate SEEN, not the last one visited, and judge
      ! convergence on it. The line search minimizes a merit; the solve is
      ! accepted on ||R||. These are different functionals, so the iterate with
      ! the smallest ||R|| is not in general the last one, and a run that
      ! reached ||R|| < resid_tol at some iterate and then wandered off has
      ! nonetheless produced a state that satisfies the code's own convergence
      ! criterion -- reporting it as a failure throws that state away and sends
      ! the caller back to time-marching (measured on examples/16: a 3.958e-5
      ! iterate discarded at iteration 2, docs/hd209_metal_stagnation.md).
      call unpack_U(Y, u);  call Apply_BC(u)
      gates_now = steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nspec_row .gt. 0, carrier_relnorm_state)
      ! On the same measure the ledger ranked with, and for the same reason:
      ! the state handed back is judged row by row against each row's own
      ! tolerance, and on the species rows too.
      keep_this_iterate = (dj_best .lt. dj)
      if ((gates_best .and. .not. gates_now) .or.                         &
          ((gates_best .eqv. gates_now) .and. keep_this_iterate)) then
         Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
         fspread = fspread_best
         n_no_chem_root_state = n_no_chem_root_best
         carrier_relnorm_state = carrier_relnorm_best
         carrier_cellmax_state = carrier_cellmax_best
         carrier_cell_state    = carrier_cell_best
         call adopt_background_of_best_iterate
         write(*,'(A,ES11.3)') ' (JFNK) returning best iterate, '//      &
              '||R||=', rnorm
      endif
      call unpack_U(Y, u)
      call unpack_species_rows(Y, u, f_sp)
      call Apply_BC(u)
      ! ONE MEASUREMENT OF THE STATE ACTUALLY RETURNED, and it is the one
      ! the gate and the certification below both read.
      !
      ! WHY IT IS TAKEN AT ALL. ||R|| is carried from the iterate that set
      ! rnorm_best, and the best-iterate restore above can replace Y without
      ! re-evaluating F; on the coupled route the carrier half of that number
      ! was moreover assembled with the limited slopes frozen at a different
      ! outer iterate. This evaluation makes F, u, f_sp and rnorm one state
      ! again.
      !
      ! WHY IT IS NOW ONE EVALUATION AND NOT A SEPARATE SWEEP LOOP. Until
      ! section 155 the residual was formed at the PREVIOUS iterate's
      ! composition, so the state handed back had to be re-measured at its
      ! own composition by a loop of fixed-Y sweeps placed here; the number
      ! the gate had read and the number that state carried then differed by
      ! up to a factor 77 (wasp_full_newton, energy row 5.466e-06 against
      ! 4.221e-04). The elimination now happens inside every residual
      ! evaluation, at the state being evaluated, so this call already
      ! returns the residual at the state's own composition and the loop is
      ! gone with the lag that made it necessary.
      weno_mode = 0
      f_sp_best = f_sp
      call eval_residual(Y, f_sp_best, f_sp_j, F, heat0, cool0)
      f_sp_best = f_sp_j
      f_sp = f_sp_best
      rnorm_handback_prev = rnorm
      call resid_relnorm(F, u, rc, rnorm)
      if (nspec_row .gt. 0) then
         carrier_relnorm_best = carrier_relnorm_last
         carrier_cellmax_best = carrier_cellmax_last
         carrier_cell_best    = carrier_cell_worst
      endif
      carrier_relnorm_state = carrier_relnorm_last
      carrier_cellmax_state = carrier_cellmax_last
      carrier_cell_state    = carrier_cell_worst
      n_no_chem_root_state = n_no_chem_root_last
      ! HOW FAR THE ACCEPTED ITERATE'S OWN NUMBER WAS FROM THE STATE IT
      ! BECAME, printed because it is the measure of the lag this item
      ! removed: with the elimination inside the evaluation the two are the
      ! same state and the ratio is 1 to the accuracy of the composition
      ! fixed point, and any departure is a composition that did not reach
      ! it within n_selfconsistent_max passes.
      write(*,'(A,I0,A,ES11.3,A,ES9.2)') ' (JFNK) residual of the state'//&
           ' handed back at its own composition, ', n_eq_sweeps_last,     &
           ' sweeps: ||R||=', rnorm, '   ratio to the accepted iterate''s'//&
           ' own ||R||=', rnorm/max(rnorm_handback_prev, 1.0d-300)

      ! THE CERTIFICATION OF THE STATE HANDED BACK. Its chemistry is part of
      ! it: n_no_chem_root_state counts the cells of that state whose accepted
      ! composition is not a root of the network (acceptance class 4 or 6), and a
      ! state resting on one of those is not a solution of the equations that
      ! were asked for, whatever its residual norms say -- heat and cool in
      ! its energy row were evaluated on a composition the chemistry does not
      ! describe. The run still writes the state; what it may not do is call
      ! it certified.
      !
      ! AND info IS A STATEMENT ABOUT THAT STATE, in both directions. The
      ! flag is set to 0 at the loop top, on the iterate; the measurement
      ! above is the one of the state handed back, and if that one refuses
      ! then so does the flag. info = 2 is what a solve that did not deliver
      ! a stationary state ends on.
      !
      ! THE GATE AND THE CERTIFICATION NOW MEASURE ONE OBJECT BY
      ! CONSTRUCTION: both read the F this routine just evaluated at the
      ! returned state, at that state's own composition. The equality is
      ! asserted rather than assumed, because the two form the row measure
      ! from separate expressions -- max_j |R_kj|/residual_row_scale(k,j) in
      ! certification_row_measure, and the larger of the same maximum over
      ! the two windows in resid_relnorm -- and nothing else keeps them
      ! together.
      call certify_returned_state(F, u, f_sp, resid_tol,                 &
           n_no_chem_root_state,                                         &
           '(JFNK) state handed back')
      cert_out = certification_last_report()
      call stationary_rows_of_the_returned_state(cert_out, rows_ok,      &
                                                 i_row, rows_unmeasured)
      call gate_equals_certification(cert_out, rnorm, ' (JFNK)',          &
                                     gate_is_cert)
      if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nspec_row .gt. 0, carrier_relnorm_state,  &
                          n_no_chem_root_state) .and. rows_ok                &
          .and. gate_is_cert) then
         info = 0
      else if (info .eq. 0) then
         info = 2
         call write_flag_is_about_the_returned_state(' (JFNK)', cert_out, &
                                                 i_row, rows_unmeasured)
      endif
      ! The frozen background the carrier transport reads must describe the
      ! state handed back, not the last state this solve happened to evaluate
      ! (docs/Update_EXHALE_stage1.md section 121).
      call install_background_of_adopted_state
      gate_rnorm_accepted = rnorm;  gate_fspread_accepted = fspread
      write(*,'(A,I0,A,ES11.3,A,ES10.3,A,I0)') ' (JFNK) done info=',info, &
           ' ||R||=',rnorm,'  flux spread=',fspread,                      &
           '  non-monotone accepts=', n_nonmonotone_accepts
      ! HOW MANY ITERATIONS WERE SPENT ON A ROW THE GATE DOES NOT SEE:
      ! iterates at which the acceptance gate was met while a row of the
      ! certification refused. A solve that ends without a certified state
      ! and reports none of these was held by the gate itself; one that
      ! reports many was held by a single row below the gate's resolution.
      if (n_gate_met_row_refused .gt. 0)                                  &
         write(*,'(A,I0)') ' (JFNK) outer iterations that met the'//      &
              ' acceptance gate while a certified row refused: ',         &
              n_gate_met_row_refused
      ! HOW OFTEN THE TWO FUNCTIONALS DISAGREED: steps the merit accepted
      ! that left the judged rows above the worst of the last five iterates.
      ! Reported and not acted on -- the ledger that chooses the returned
      ! state is where the judgement enters (see there).
      write(*,'(A,I0)') ' (JFNK) accepted steps that raised the judged'// &
           ' rows above their five-iterate window: ', n_judged_excursions
      ! AND THE LEDGER'S OWN NUMBER: the smallest distance from
      ! certification any iterate of this solve reached, which is the
      ! iterate the restore above hands back. It is the minimum of every
      ! iteration's `judged rows` line by construction, and it would not be
      ! if the ledger ranked on ||R|| while the binding row was a species
      ! row, or a hydrodynamic row other than the one holding ||R||.
      write(*,'(A,ES11.3)') ' (JFNK) best judged iterate: distance',      &
           dj_best
      ! Both norms on the state handed back, with the window contributions
      ! and the worst cell's terms (section 145).
      block
        real*8, dimension(3,1-Ng:N+Ng) :: Rp145
        integer :: jp145, kp145
        Rp145 = 0.0d0
        do jp145 = 1, N
           do kp145 = 1, 3
              Rp145(kp145,jp145) = F(nvar_jac*(jp145-1)+kp145)
           enddo
        enddo
        call write_residual_breakdown(Rp145, u, heat0, cool0,             &
             'JFNK state handed back')
      end block

      ! Where the non-flatness reaches, for reading only: the acceptance test
      ! is steady_gates_met on the r >= r_flux window alone and it is unchanged.
      call flux_spread_above_radius(u, 1.03d0, fspread_103, dum_fmean)
      call flux_spread_above_radius(u, 1.10d0, fspread_110, dum_fmean)
      write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A)')                           &
           ' (JFNK) flux spread by window: r>=1.03', fspread_103,         &
           '   r>=1.10', fspread_110, '   r>=r_flux', fspread,            &
           '  (only the last is the gate)'
      if (nspec_row .gt. 0) then
         write(*,'(A,ES10.2,A,ES10.2,A,I4,A,I0)') ' (JFNK) carrier row '//&
              'of the returned state: volume-weighted',                  &
              carrier_relnorm_best, ', worst cell',                      &
              carrier_cellmax_best, ' at j=', carrier_cell_best,         &
              '; trial and probe states refused by an admissibility'//  &
              ' screen: ', n_trial_headroom
         ! WHICH SCREEN REFUSED THE SAMPLES THIS SOLVE DID NOT GET. The
         ! first entry is the samples that were accepted; the rest are the
         ! Newton columns, Krylov directions and trials the solve was
         ! refused, by cause.
         write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')                  &
              ' (JFNK) residual samples: admitted ', n_eval_refusal(0),   &
              ', refused for a negative species unknown ',                &
              n_eval_refusal(refuse_species_unknown_negative),            &
              ', an element fraction above one ',                         &
              n_eval_refusal(refuse_element_fraction_above_1),            &
              ', a breach of the shared element budget ',                 &
              n_eval_refusal(refuse_element_budget_breach),               &
              ', a non-finite sweep ',                                    &
              n_eval_refusal(refuse_sweep_nonfinite),                     &
              ', a non-finite row ',                                      &
              n_eval_refusal(refuse_row_nonfinite),                       &
              ', no chemical root ',                                      &
              n_eval_refusal(refuse_no_chemical_root)
         write(*,'(A,I0,A)') ' (JFNK) Jacobian-vector products taken on'//&
              ' the backward side of their direction: ',                  &
              n_jv_backward_sample, ' (no forward step along them landed'//&
              ' on a describable state)'
         write(*,'(A,A,A,I0,A)') ' (JFNK) species unknown space: ',       &
              trim(species_unknown_space_text()),                         &
              '; trials written onto the faces of the species box: ',     &
              n_projected_trials, ' (a solve spending its steps on those'//&
              ' faces is not descending the interior of the set)'
         write(*,'(A,ES11.3,A)') ' (JFNK) smallest species unknown any'// &
              ' adopted iterate carried: ', species_unknown_min_adopted,  &
              ' (a density and a mass fraction are non-negative)'
         ! THE SAME STATEMENT AS A RATIO, which is the one that separates the
         ! two unknown spaces: a carrier the solve drove to zero reads zero
         ! here whatever its own scale is, and one that stayed inside the
         ! space reads how many decades above its floor it stopped.
         if (carrier_rows_registered())                                   &
            write(*,'(A,ES11.3,A,I0,A)') ' (JFNK) smallest adopted'//     &
              ' carrier density over the floor of its own cell: ',        &
              carrier_over_its_floor_min,                                 &
              '; adopted carriers at or below that floor: ',              &
              n_carrier_at_its_floor,                                     &
              ' (the floor is 1e-20 of the cell'//"'"//'s element budget)'
         ! THE STATEMENT AT THE OTHER FACE, which is the one this item made
         ! a face: a carrier cannot hold more nuclei of its element than the
         ! cell has free, so at or below one and a count of zero is what an
         ! admissible solve reports.
         if (carrier_rows_registered())                                   &
            write(*,'(A,ES11.3,A,I0,A)') ' (JFNK) largest adopted'//      &
              ' carrier density over the element budget of its own cell: ',&
              carrier_over_its_budget_max,                                &
              '; adopted carriers above that budget: ',                   &
              n_carrier_above_its_budget,                                 &
              ' (the budget is the free density of the carrier'//"'"//    &
              's element divided by the nuclei one carrier holds)'
         write(*,'(A,I0,A)') ' (JFNK) most species unknowns held on an'//  &
              ' active bound by one step: ', n_active_bounds_max,          &
              ' (their rows are constrained, not rooted)'
         ! THE FACE THE ELEMENT BUDGET IS, and how often it had to be raised
         ! to the iterate because the iterate already stood outside it
         ! (freeze_species_unknown_box).
         if (carrier_rows_registered())                                   &
            write(*,'(A,I0,A)') ' (JFNK) carrier cells whose element'//   &
              ' budget face was raised to the iterate: ',                 &
              n_budget_face_at_the_iterate,                               &
              ' (summed over the outer iterations; the budget is frozen'// &
              ' at the iterate and the iterate can stand outside it)'
         ! THE SHARED ELEMENT ROWS OF THE STEP: how often the point the step
         ! started from already stood outside its own frozen budget, how
         ! many rows carried no direction the rows of their own cell did not
         ! already carry, how many candidates had to be stepped back onto
         ! the constraint and by how much, and how often the write-back
         ! could not give an element total back.
         if (carrier_rows_registered()) then
            write(*,'(A,I0,A,I0)') ' (JFNK) shared element rows: cell'//  &
              '-rows whose iterate stood outside its own budget ',        &
              n_elem_row_target_at_the_iterate,                          &
              ', rows dropped as already spanned ',                       &
              n_elem_row_dropped_for_rank
            write(*,'(A,I0,A,I0,A,ES10.3,A,ES10.3)')                     &
              ' (JFNK) restoration phase: trials restored ',             &
              n_elem_restoration_steps, ', cell-steps ',                 &
              n_elem_restoration_cells,                                  &
              ', worst violation removed from ',                         &
              elem_restoration_worst_before, ' to ',                     &
              elem_restoration_worst_after
            write(*,'(A,I0,A,ES10.3)') ' (JFNK) write-back cells whose'//&
              ' closure species could not return an element total: ',    &
              n_write_back_element_deficit, ', worst relative deficit ', &
              write_back_worst_element_deficit
            write(*,'(A,I0,A)') ' (JFNK) probe components with no room'//&
              ' at any step length: ', n_probe_blocked_components,       &
              ' (zeroed in the direction rather than costing it whole)'
         endif
         write(*,'(A,I0,A)') ' (JFNK) [TR] steps that were the Cauchy'//   &
              ' point alone: ', n_cauchy_only_steps,                       &
              ' (no Krylov direction could be sampled there)'
         write(*,'(A,I0,A)') ' (JFNK) [TR] dogleg trials that left the'// &
              ' approximate-gradient leg out of the step: ',              &
              n_cauchy_leg_short,                                          &
              ' (shorter than cauchy_leg_min_fraction of the radius;'//   &
              ' the leg costs one Jacobian product per outer iteration)'
         if (cauchy_leg_admitted_by_image_alone)                          &
            write(*,'(A)') ' (JFNK) [TR] the length test was not taken:'//&
              ' the leg was admitted by its share of the model image'//   &
              ' alone (EXHALE_CAUCHY_LEG_BY_IMAGE)'
         write(*,'(A,I0,A)') ' (JFNK) [TR] dogleg trials that left the'// &
              ' approximate-gradient leg out because it was the image'//  &
              ' and not the step: ', n_cauchy_leg_all_image,               &
              ' (above cauchy_leg_image_share_max of the image, below'//  &
              ' cauchy_leg_step_share_max of the step)'
         write(*,'(A,I0,A,I0,A,I0,A)') ' (JFNK) cost: ', nvar_jac,       &
              ' unknowns per cell, ', ncolor_jac, ' colors, ',           &
              n_resid_eval, ' residual evaluations'
      endif
      call write_no_chem_root_trial_count(' (JFNK)')
      ! WHAT THE COMPOSITION ELIMINATION COST. A residual evaluation is no
      ! longer one equilibrium sweep: it is as many as the composition needs
      ! at fixed Y, and the ratio is the price of a residual that is the
      ! residual of its own composition.
      if (n_resid_eval .gt. 0)                                            &
         write(*,'(A,I0,A,I0,A,F6.2,A,I0,A)') ' (JFNK) composition'//     &
              ' elimination: ', n_eq_sweeps_tot, ' sweeps over ',        &
              n_resid_eval, ' residual evaluations, ',                   &
              real(n_eq_sweeps_tot,8)/real(n_resid_eval,8),              &
              ' per evaluation; the Newton model carried ',              &
              n_eq_sweeps_model, ' Picard term(s)'
      ! WHICH GATES ARE LEFT -- each one asked separately, because there are
      ! three of them and "the one left" was a guess about which.
      !
      ! The chain this replaces tested the residual and the flux gate TOGETHER
      ! first and then fell through to the carrier row, so a state that met the
      ! residual gate and missed the other two printed "the CARRIER gate is the
      ! one left" while the flux gate was also unmet, by a factor 4.9.
      ! Measured on the 1 microbar hot Uranus (docs/p23_transport_on_state.md
      ! section 6): ||R|| = 1.547e-6 against 1e-5, carrier row 7.535e-3 against
      ! 1e-3, flux spread 2.447e-2 against 5e-3 -- two gates left, one named.
      ! One line per unmet gate says what is true whichever combination occurs.
      if (info .ne. 0) then
         if (rnorm .ge. resid_tol)                                        &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'residual ||R|| =', rnorm, ' >=', resid_tol
         if (flux_spread_th .gt. 0.0d0 .and. fspread .ge. flux_spread_th) &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'flux spread', fspread, ' >=', flux_spread_th
         if (nspec_row .gt. 0 .and.                                        &
             carrier_relnorm_best .ge. carrier_resid_th)                  &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'carrier row', carrier_relnorm_best,                     &
                 ' >=', carrier_resid_th
         if (n_no_chem_root_state .gt. 0)                                 &
            write(*,'(A,I0,A)') ' (JFNK) gate NOT met: ',                 &
                 n_no_chem_root_state, ' cell(s) of the state handed'//   &
                 ' back carry a composition that is not a root of the'//  &
                 ' chemical network (acceptance class 4 or 6); the'//     &
                 ' state is'// &
                 ' written but NOT certified'
      endif
      deallocate(Y,F,F_jac,Ftry,dY,Ytry,dZ,D,Ybest,dY_lm,Drow,ab,abf,ipiv)
      end subroutine solve_steady_jfnk

      ! End of module
      end module steady_newton

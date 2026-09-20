      module certification
      ! ONE EVALUATOR OF THE EQUATIONS A RUN IS ACTUALLY SOLVING, AND THE
      ! CONTEXTS THAT READ IT.
      !
      ! docs/a2_certification_contract_20260906.md is the contract; this is
      ! steps 1 to 3 of its section 7. The inventory of what a given run
      ! solves is docs/b1a_active_equation_inventory_20260906.md section 2,
      ! and the run states are docs/a0_run_mode_contract_20260906.md
      ! section 2.
      !
      ! WHAT IT IS FOR. Until now "converged" was a statement about three
      ! hydrodynamic rows, a mass-flux spread, sometimes a carrier norm and,
      ! since the acceptance-class contract, a count of cells without a
      ! chemical root. A run that also solves an elemental transport balance,
      ! a level population or a set of eliminated species reached its answer
      ! with nothing said about those equations at all. The evaluator here
      ! lists EVERY equation the configuration makes independent and returns,
      ! for each, one of three statuses:
      !
      !   not_applicable  the configuration does not carry the unknown
      !   evaluated       a residual of that equation was measured ON THIS
      !                   STATE, with its row measure, its scale and its floor
      !   unavailable     the equation is active and no state-consistent
      !                   measurement of it exists; NEVER reported as a zero
      !
      ! The third status is the point. carrier_steady_residual returns zeros
      ! when the frozen background is not ready, which at the call site is
      ! indistinguishable from a balanced carrier; an active equation with no
      ! evaluator returns nothing at all. Both now say so, and a state
      ! carrying either is written and NOT certified.
      !
      ! WHAT IS EVALUATED. Every equation of the inventory that the
      ! configuration makes independent:
      !
      !   the three hydrodynamic rows          (B1a 2.1)
      !   each transported carrier balance     (B1a 2.2)
      !   the He/H partition and the trace     (B1a 2.3)
      !     elemental transport balances
      !   the He 2^3S and H(n=2) level         (B1a 2.5)
      !     balances
      !   the eliminated-species closure of    (B1a 2.4)
      !     the coupled cell system
      !
      ! The last three kinds are the state-consistent evaluators of step 4
      ! of the contract: a stationary residual of the elemental transport
      ! (element_transport_residual), a residual of the 2s/2p level balance
      ! at the state's field and populations
      ! (excited_hydrogen_level_residual), and the reaction residual of the
      ! closure at the state's composition and the sweep's own rates
      ! (ionization_closure_residual_profile). None of them solves anything
      ! and none of them writes f_sp, rho or a ledger.
      !
      ! The closure measure IS the measure the cell sweep judges a root by,
      ! so a cell of acceptance class 1, 2, 3 or 5 reads at or below
      ! ieq_res_tol and a class-4 or class-6 cell above it, and a certification is
      ! refused when any cell is above it.
      !
      ! THE ROW MEASURE. The maximum over the physical cells of
      ! |res| / max(scale, floor), never a volume average: a front or a thin
      ! layer must not be averaged away by the wind that carries most of the
      ! volume. The volume average is reported beside it, and is never
      ! decisive (contract section 4; the carrier `rvol` is the case review
      ! R7 raised). The scales are the ones the residual code already
      ! defines: residual_row_scale for the hydrodynamic rows (the largest
      ! term the row itself contains) and row_terms for a carrier row (the
      ! sum of the magnitudes of the terms that row balances, with the
      ! absolute floor carrier_residual builds into it).
      !
      ! NO SIDE EFFECTS. The carrier balance is measured inside an isolated
      ! workspace: every module array carrier_steady_residual writes is held
      ! aside before the call and put back after it, allocation status
      ! included (save_carrier_module_state / restore_carrier_module_state,
      ! diffusive_photochemistry.f90). A certification is a measurement of a
      ! state, not an operation on it.

      use global_parameters
      use Conversion,               only: U_to_W
      use steady_residual_mod,      only: residual_row_scale,             &
                                          n_cells_without_chemical_root,  &
                                          face_mass_flux_of_state,        &
                                          mass_row_rounding_floor,        &
                                          row_terms_describe_state
      use diffusive_photochemistry, only: carrier_steady_residual,        &
                                          n_carrier, n_carrier_max,       &
                                          carrier_name, carrier_solved,   &
                                          carrier_co_domain_record,       &
                                          carrier_module_state,           &
                                          save_carrier_module_state,      &
                                          restore_carrier_module_state,   &
                                          carrier_module_state_matches,   &
                                          carrier_roundoff_limited_record,&
                                          carrier_history_certifiable,   &
                                          carrier_row_terms_on,          &
                                          carrier_row_terms_write,       &
                                          ionization_stage_sum_measure
      use ionization_equilibrium,   only: bg_ready, ieq_nonroot_streak,   &
                                          finite_real, ieq_res_tol,       &
                                          ieq_triplet_row,                &
                                          transported_rows_exist,         &
                                    ionization_closure_residual_profile
      use excited_hydrogen,         only: excited_hydrogen_level_residual
      use binary_element_diffusion, only: element_transport_residual
      use element_inventory,        only: ien_H, ien_He
      ! What identity (2) of ionization_stage_transport can stand at in
      ! floating point. The identity is algebraic, so this is the whole of
      ! what the stage sum entry can clear, and it is formed there and not
      ! here.
      use ionization_stage_transport, only:                               &
                                    ionization_stage_sum_rounding_bound
      use composition,              only: get_species_densities,          &
                                          comp_T_from_p
      use species_table,            only: n_mion, n_melem, melem_name,  &
                                          mion_fsp, isp_HI, isp_HII,    &
                                          isp_HeI, isp_HeII, isp_HeIII, &
                                          isp_H2, isp_H2p, isp_H3p,     &
                                          isp_HeHp, isp_OH, isp_H2O,    &
                                          isp_CO
      ! The H3+ cooling model's own domain records (B3b-H3): four counters
      ! of the cells at which a collider density or a temperature left the
      ! range the published fits cover. Category (2) of B1a section 4.
      use h3p_cooling,              only: h3p_n_below_collider,          &
                                          h3p_n_below_fit_T,            &
                                          h3p_n_above_fit_T,            &
                                          h3p_n_outside_nonlte_T
      ! The domain record of the Ly-alpha damping-wing closure (LYA-BETA):
      ! cell visits whose (a tau)^(1/3) fell below the published limit of
      ! the slab solution (Neufeld 1990 section Va); there the closure is
      ! carried by its thin limit. Category (2) of B1a section 4.
      use lya_rt,                   only: lya_wing_domain_record
      use energy_semi_implicit,     only: n_energy_floor_hits,            &
                                          n_energy_floor_cells,           &
                                          energy_update_last_status,      &
                                          energy_update_last_cell,        &
                                          energy_update_last_residual,    &
                                          ENERGY_UPDATE_OK
      use viscous_conduction,       only: n_conduction_floor_hits,        &
                                          n_conduction_floor_cells,       &
                                          conduction_last_status,         &
                                          conduction_last_cell,           &
                                          CONDUCTION_OK
      use utils,                    only: set_state_certified, calc_rho
      use stationary_operator,      only: face_mass_flux_budget
      ! The direction the base contact was upwinded on at its last
      ! evaluation, and what decided it (item D2b).
      use base_boundary,            only: base_face_swind_last,           &
                                          base_face_Mwind_last,           &
                                          base_face_blend_last

      implicit none
      private

      ! Evaluation status of one equation of the inventory.
      integer, parameter, public :: cert_not_applicable = 0
      integer, parameter, public :: cert_evaluated      = 1
      integer, parameter, public :: cert_unavailable    = 2

      ! The contexts of contract section 3 that this step installs.
      integer, parameter, public :: cert_context_probe      = 1
      integer, parameter, public :: cert_context_trial      = 2
      integer, parameter, public :: cert_context_stationary = 3
      ! THE PHYSICAL-STEP CONTEXT of contract section 3, row "Physical time
      ! step": the adoption boundary of the marching step (B3a). It is a
      ! context of THIS evaluator and not a second implementation, so the
      ! tolerances a physical step is judged by and the tolerances a
      ! stationary state is certified by are one set of numbers (advisor
      ! decision 5, docs/b3a_attempted_step_controller_design_20260906.md
      ! section 9). What it asks for is different from what the stationary
      ! context asks for -- the time-discrete balances of the step, not the
      ! stationary residual, which a valid finite-time state does not have
      ! to zero -- but the row measure and the floors are the same.
      integer, parameter, public :: cert_context_physical_step = 4

      ! THE NUMBER OF COUPLED CELL SYSTEMS of B1a section 2.4. Exactly one
      ! of them is active in a run; the report carries an entry for each so
      ! that it says which one was measured and which were not the run's.
      ! active_system_variant selects the active one and
      ! system_variant_name names them.
      integer, parameter, public :: n_system_variant = 7

      ! THE CAPACITY OF A REPORT IS THE SIZE OF THE INVENTORY. It is not a
      ! round number: it is the count of equations the configuration can
      ! make independent at once, summed over the kinds of section 2 above,
      ! so that a report has room for every equation of the largest
      ! configuration and for nothing the inventory does not have.
      !
      ! RAISING THE NUMBER IS NOT THE REPAIR. What an overfilled report
      ! costs is that an entry writer attaches a new equation's measure,
      ! tolerance and verdict to the name of the equation before it, and a
      ! certification then states a verdict about a set of equations it does
      ! not name. A larger constant leaves that available one equation
      ! later. The contract is add_entry's, which refuses to overfill at
      ! all, and this count, which follows the inventory when a carrier, a
      ! trace element or a cell system is added to it.
      integer, parameter, public :: cert_max_entries =                    &
             3                &  ! hydrodynamic mass, momentum, energy rows
           + n_carrier_max    &  ! one balance for each transported carrier
           + 2                &  ! the stage nucleus sums, H and He
           + 1                &  ! the He/H partition transport balance
           + n_melem          &  ! one transport balance for each trace element
           + 2                &  ! the He 2^3S and H(n=2) level balances
           + n_system_variant    ! the eliminated-species closures

      ! TOLERANCES, one per kind of row, named so that the report says which
      ! number refused a state. THE CONVERGENCE STUDY OF CONTRACT SECTION 9
      ! FIXES THEM (A2tol, 2026-09-06); none was chosen to make a present
      ! snapshot pass, and section 9.6 of that study records that most of
      ! the matrix does not pass at these values.
      !
      ! ONE NUMBER PER HYDRODYNAMIC ROW, NOT ONE FOR THE THREE. On the same
      ! converged state the three rows reach values seven decades apart
      ! (MEASURED, wasp_full_newton at "Resid tol" 1e-7: mass 2.630e-13,
      ! momentum 1.280e-09, energy 4.700e-06) and their one-ulp arithmetic
      ! floors are two decades apart (1.018e-13, 5.127e-12, 1.643e-11). A
      ! single number either certifies a mass row eight decades above its
      ! own floor or refuses an energy row no solve in this code has
      ! reached. Each value below is that measured anchor with about a
      ! decade of margin, and each stands at least a decade above its row's
      ! floor: 29x for mass, 2000x for momentum, 3e6x for energy.
      !
      ! THEY NO LONGER INHERIT "Resid tol". Section 9.2 MEASURED a
      ! structural factor of 41 to 47 between the target a JFNK solve is
      ! given and the row the state it hands back actually carries, because
      ! the completion flag was decided before the state's own residual was
      ! re-measured. An inherited resid_th was therefore a condition the
      ! solver could not meet whatever the state was.
      !
      ! RE-ANCHORED AFTER B5a (contract section 10, A2tol2, 2026-09-06).
      ! With the self-consistent residual the returned row is 0.36 to 0.93
      ! of the solve's target, and wasp_full_newton at "Resid tol" 1e-8
      ! reaches mass 2.6e-13, momentum 1.3e-9, energy 6.275e-10. Mass and
      ! momentum keep their section-9 values. The energy value moves from
      ! 5e-5, which had been anchored on the lagged-gate energy row
      ! 4.700e-06 and refused nothing, to 1e-6: a decade above 1.212e-07,
      ! the largest energy row of any state in that study whose solve met
      ! its own target, and 5.7e4x above its floor 1.765e-11.
      real*8, parameter, public :: cert_tol_mass     = 3.0d-12
      real*8, parameter, public :: cert_tol_momentum = 1.0d-8
      real*8, parameter, public :: cert_tol_energy   = 1.0d-6
      ! THE MASS ROW'S TOLERANCE IS A FUNCTION OF THE CELL, cert_tol_mass_at
      ! below, and every reader of it goes through that accessor. The
      ! constant above is the FLOOR of that function: where the arithmetic
      ! permits it, the row is held to 3e-12 and to nothing looser.
      !
      ! WHY A FIXED NUMBER ALONE REFUSES STATES THAT CARRY NO SIGNAL. The
      ! mass row is a difference of two face fluxes over a cell volume, and
      ! in a nearly hydrostatic base layer that difference is formed
      ! between quantities the Riemann solve builds from the cell's whole
      ! momentum at its own signal speed, so its rounding is about eps/Mach
      ! of the flux and not eps of it (mass_row_rounding_floor,
      ! steady_residual.f90, states the derivation; N33 measured the
      ! phenomenon). MEASURED at the base of the three fixtures of item
      ! P16: the estimated floor is 5.7E-14 on `wasp_full_newton` (Mach
      ! 3.9e-3), 9.3E-13 on the hot-Uranus carrier reload (Mach 2.4e-4) and
      ! 3.2E-11 on the HD 209458 b element reload (Mach 5.1e-6), four
      ! decades across three states of one code. A single 3e-12 therefore
      ! holds one of them to eight ulps of its own base layer and another
      ! to a tenth of one ulp.
      !
      ! c_round, THE ONE NUMBER ANCHORED, and how it was measured. One ulp
      ! is added to the density of every cell and the whole flux assembly is
      ! run again; the STEP the row takes is the non-smoothness the assembly
      ! carries, and the ratio of that step to the estimate above is what
      ! the margin has to cover. MEASURED (P16, EXHALE_MASS_FLOOR_SCAN, the
      ! largest ratio over all 500 cells of each state as loaded):
      ! `wasp_full_newton` 1.59, HD 209458 b 1.89, hot-Uranus 2.68. The
      ! value below stands 3.7 times above the largest of them, and N33's
      ! independent reading of the same floor against its own bound was 5.5
      ! to 12, so a decade is the margin the two measurements together
      ! support; nothing larger is measured and nothing larger is taken.
      real*8, parameter, public :: cert_mass_round_margin = 1.0d1
      ! THE CEILING OF THE ANCHORED TOLERANCE. |R_1|/s_1 is the fractional
      ! change of the mass flux across the cell, so a row at 1 says the
      ! flux changes by the whole of itself there. No arithmetic argument
      ! admits that, whatever the estimated floor says, and a state whose
      ! estimated floor reaches this value is one whose mass row cannot be
      ! judged at all rather than one that passes. It is not an operative
      ! branch on any state measured: it would need a base Mach number
      ! below 1e-15.
      real*8, parameter, public :: cert_tol_mass_ceiling = 1.0d0
      ! THE TWO OUTCOMES OF THE MASS ROW AT ONE CELL (review R5, item Q3).
      ! `cert_tol_mass_ceiling` above says a cell whose estimated floor
      ! reaches it cannot be judged, not that it passes; mass_row_cell_
      ! verdict below returns one of these through its optional `status`,
      ! and a cell it names UNRESOLVED never reads `within`, whatever q is.
      integer, parameter, public :: cert_mass_row_resolved   = 1
      integer, parameter, public :: cert_mass_row_unresolved = 2
      ! THE TWO SPECIES-ROW TOLERANCES ARE VALUES OF A FUNCTION OF THE
      ! RADIUS, cert_tol_carrier_at and cert_tol_element_at below, and every
      ! reader of them goes through those two accessors. The constants here
      ! are the value that function takes in the wind.
      real*8, parameter, public :: cert_tol_carrier_wind = 1.0d-5
      ! The closure and the two level balances are measured in units of the
      ! turnover of their own row, so their tolerance is the one the cell
      ! sweep already judges a root by: ieq_res_tol, the value the measured
      ! acceptance populations of the whole regression matrix fix (its
      ! declaration states the measurement). Using any other number here
      ! would mean a cell the sweep called a root and the certification
      ! called unbalanced, or the reverse.
      real*8, parameter, public :: cert_tol_closure = ieq_res_tol
      real*8, parameter, public :: cert_tol_level   = ieq_res_tol
      ! The elemental transport balance in the wind, on the composition
      ! Newton's own measure: composition_residual returns the residual
      ! relative to the size of the terms of its own row, which is the
      ! construction of row_terms, so the element and the carrier row are
      ! measured alike and take one value.
      real*8, parameter, public :: cert_tol_element_wind = 1.0d-5
      ! WHAT A CELL OUTSIDE THE GATING REGIME IS JUDGED BY: nothing. Its row
      ! is measured and reported and enters no verdict, and dividing its
      ! measure by this value is what "it does not gate" is arithmetically.
      real*8, parameter, public :: cert_tol_reported_only = huge(1.0d0)
      ! THE TIME-DISCRETE HYDRODYNAMIC ROW of the physical-step context:
      !
      !     |u_new - u_old - dt L(u)|  <=  cert_tol_hydro_step
      !                                    * (|u_old| + dt |L(u)|)
      !
      ! per row and per physical cell, with L(u) = -(dF - S) the operator
      ! the three Runge-Kutta stages integrate, evaluated at the state the
      ! step started from. The value is INHERITED, not measured, and it is
      ! stated here as a parameter so that a reader sees what it is: a
      ! third-order SSP-RK3 update of a smooth state differs from one Euler
      ! evaluation of the same operator by O(dt) times the row's own terms,
      ! which is of order one in this relative measure, so this is a LOOSE
      ! ADMISSIBILITY STATEMENT about the update and NOT a convergence
      ! statement. What carries the integration error of the step is the
      ! step-doubling estimate of the controller, not this row. A
      ! convergence study (contract section 4) replaces it; nothing was
      ! chosen to make a snapshot pass.
      real*8, parameter, public :: cert_tol_hydro_step = 1.0d0
      !
      ! WHY THIS ROW IS EVALUATED AND REPORTED AND DOES NOT GATE, in the
      ! interim. MEASURED on hydrostatic_column at its cold start
      ! (2026-09-06, this item): the measure reaches 1.92 in cell 2 and 1.70
      ! in cell 1 over the first four steps and is below 1 everywhere
      ! afterwards. Those are the base cells of a HYDROSTATIC column, where
      ! the momentum of the state is zero, so |u_old| vanishes for that row
      ! and the scale is dt|L(u_old)| alone; the measure is then the
      ! stage-to-stage variation of the operator, L(u2) against L(u_old),
      ! and an SSP-RK3 update has no bound of order one on that: its stages
      ! are convex combinations of Euler steps, so the residual against one
      ! Euler evaluation is bounded by 1 + max_stage|L|/|L(u_old)|, which is
      ! large wherever the right-hand side moves across the stages.
      !
      ! So no defensible number for this tolerance exists yet, and choosing
      ! one that admits the measured startup would be choosing a tolerance
      ! to make a snapshot pass. The row is measured every step and
      ! reported with its worst cell; it becomes a gate when the
      ! convergence study of contract section 4 fixes the number. The value
      ! above is the order-one bound a naive reading of the measure
      ! suggests, kept so that the report says what the row is being
      ! compared against.
      !
      ! What carries the integration error of the step meanwhile is the
      ! step-doubling estimate of the controller, which is a property of the
      ! whole operator split and needs no tolerance of this kind.

      ! The floor under every row scale: it exists only so that a division
      ! by an identically zero scale cannot produce a NaN, and it is far
      ! below any physical scale of either row.
      real*8, parameter, public :: cert_scale_floor = 1.0d-300

      ! THE TWO REGIMES A SPECIES ROW IS READ IN. The transport balance of
      ! an element or a carrier is a cancellation between terms whose ratio
      ! is set by the flow: in the wind the advective divergence carries the
      ! row and the diffusive one is a correction to it, while below the
      ! homopause the eddy and settling terms carry it and each is orders of
      ! magnitude larger than their difference. The precision the same
      ! equation can be resolved to is therefore not one number, and the
      ! report states the two separately instead of averaging them into one
      ! threshold (docs/certification_tolerance_anchoring_20260910.md).
      !
      ! The radii: r_layer is the top of the diffusion-dominated column of
      ! this code's planets (the K_zz term of binary_element_diffusion still
      ! exceeds the advective one there), r_wind the radius above which the
      ! flow is transonic on every case of the matrix.
      real*8, parameter, public :: cert_regime_layer_r = 1.10d0
      real*8, parameter, public :: cert_regime_wind_r  = 1.20d0

      ! WHICH CELLS A SPECIES ROW IS GATED IN, and what the number is
      ! (decision 22 (a), user 2026-09-10; the five anchors are MEASURED in
      ! docs/certification_tolerance_anchoring_20260910.md sections 2, 4
      ! and 9).
      !
      ! GATING, at r >= cert_regime_wind_r: 1e-5. It is the coarser of the
      ! two floors the element operator states for itself (about 1e-12 in a
      ! wind and about 1e-5 in a K_zz-dominated column,
      ! binary_element_diffusion), it is a decade above the discretization
      ! error a smooth manufactured column shows at the production spacing
      ! (1.52e-6 at dr = 1.84e-3 R_p, MEASURED), and it is one and a half
      ! decades below what the two candidate states carry in that window
      ! (element 3.05e-4, carrier 7.19e-2, MEASURED), so it admits no state
      ! that exists.
      !
      ! REPORTED AND NOT GATING, below cert_regime_wind_r. Two reasons, one
      ! per part of that column:
      !   r < cert_regime_layer_r  the element fluxes of the operator are
      !     not conserved there at all: their radial spread over that
      !     window is 9.8 times the median (MEASURED, after the mass-closure
      !     repair; 39 before it), which is the standing base sound wave and
      !     not a tolerance question. No threshold means anything against
      !     it, and the window gates nothing until it is repaired.
      !   between the two radii  the discretization of the row is 4.6e-4 at
      !     the radius that binds (r = 1.153, 78 percent of the row measured
      !     there, MEASURED by grid refinement), so a 1e-5 gate would ask a
      !     residual to fall one and a half decades below the error of the
      !     discrete equation it is the residual of. No separate anchor for
      !     the band exists, and inventing one would be choosing a number
      !     rather than measuring it.
      ! A state whose wind rows are within tolerance is therefore certified
      ! IN THE WIND, and the report and the state file say so in one line
      ! with the worst reported row beside it.

      ! What the run asked its stationary solver for ("Resid tol"), kept so
      ! that a report can name the solver's target beside the row measures.
      ! It is NOT a certification tolerance.
      real*8, save, public :: cert_resid_tol_of_run = -1.0d0

      type, public :: cert_entry
         character(len=52) :: name        = ''
         integer           :: status      = cert_not_applicable
         character(len=76) :: reason      = ''
         ! max over the physical cells of |res|/max(scale, floor)
         real*8            :: row_max     = 0.0d0
         ! the same rows volume weighted, reported beside row_max, never
         ! decisive
         real*8            :: row_vol     = 0.0d0
         integer           :: jworst      = 0
         logical           :: finite      = .true.
         real*8            :: tol         = 0.0d0
         logical           :: within_tol  = .true.
         character(len=76) :: units_floor = ''
         ! INFORMATIONAL, NEVER DECISIVE. Rows whose residual already sits
         ! at the arithmetic round-off of their own full terms, which is
         ! what an exact solve leaves there: they are accepted and flagged
         ! (A1scale3; review 2 section 5.3 calls them informational).
         ! Reported beside the measure so that a reader can tell a row the
         ! solve did not resolve from a row it got wrong; it can never
         ! invalidate a state.
         integer           :: n_roundoff_limited = -1
         ! CELLS IN WHICH THE SPECIES OF THIS ROW IS ABSOLUTELY ABSENT, and
         ! which therefore report and do not gate (carrier_row_entry). A
         ! carrier below carrier_absent_fraction of the free reservoir of
         ! its own element is not a species of the state: it carries no
         ! mass, no charge, no opacity and no energy any reported quantity
         ! depends on, and its balance is a statement about round-off.
         integer           :: n_absent = 0
         integer           :: j_first_absent = 0
         ! THE SAME MEASURE OVER THE TWO REGIMES, and at one probe radius.
         ! Reported, never decisive: row_max above is the whole column and
         ! remains the only number the verdict reads. These say WHERE the
         ! distance is, which is what a tolerance can be anchored on.
         real*8            :: row_max_layer = 0.0d0
         integer           :: jworst_layer  = 0
         real*8            :: row_max_wind  = 0.0d0
         integer           :: jworst_wind   = 0
         ! The measure at the cell nearest the probe radius the run names
         ! (EXHALE_CERT_ANCHOR_R), so one physical radius can be followed
         ! across two grids. jprobe = 0 when no radius was named.
         real*8            :: row_at_probe  = 0.0d0
         integer           :: jprobe        = 0
         ! WHETHER THIS ROW'S TOLERANCE IS A FUNCTION OF THE RADIUS. True
         ! for the elemental transport and carrier rows alone: their
         ! verdict is taken over the cells the accessors gate and the rest
         ! of the column is reported beside it. row_max above stays the
         ! whole column for every row, gated or not.
         logical           :: regime_gated  = .false.
         ! The gated part of such a row: the largest measure over the
         ! gated cells and the cell that carries it. jworst_gate = 0 on a
         ! gated row says no cell of the column is gated at all, which is
         ! a row the state cannot be judged on.
         real*8            :: row_max_gate  = 0.0d0
         integer           :: jworst_gate   = 0
         ! The other side of the same split: the largest measure over the
         ! cells the accessor does NOT gate, r < cert_regime_wind_r. It is
         ! the whole reported part of the column and not the layer window
         ! alone, because the band between the two radii holds the cell
         ! that binds the element row of a real state (r = 1.15).
         real*8            :: row_max_reported = 0.0d0
         integer           :: jworst_reported  = 0
         ! WHERE A TOLERANCE THAT IS A FUNCTION OF THE CELL BINDS. The
         ! continuity row's tolerance is one (cert_tol_mass_at), so the
         ! cell whose measure stands furthest outside ITS OWN tolerance is
         ! not in general the cell of the largest measure, and the verdict
         ! is taken cell by cell. jbind = 0 says the row's tolerance is one
         ! number for the whole column and jworst carries the verdict.
         integer           :: jbind        = 0
         real*8            :: row_at_bind  = 0.0d0
         real*8            :: dist_bind    = 0.0d0
         ! The tolerance that applied at the cell of the LARGEST measure,
         ! which is the other cell a reader looks for when the two differ:
         ! it says why a row_max above the fixed value was admitted.
         real*8            :: tol_at_jworst = 0.0d0
         ! Whether the tolerance at that cell came from the rounding anchor
         ! or from the fixed floor, so a reader can tell which gate
         ! refused or admitted the state.
         logical           :: tol_anchored = .false.
         ! MASS ROW ONLY (review R5, item Q3): how many cells of the column
         ! carry an estimated rounding floor that reaches the ceiling, and
         ! the first such cell in cell order, so the report can name a cell
         ! whose balance was never judged even when a different cell binds
         ! the verdict. Zero and 0 on every other row, and on a mass row
         ! with no unresolved cell.
         integer           :: n_mass_unresolved = 0
         integer           :: j_first_mass_unresolved = 0
      end type cert_entry

      type, public :: cert_report
         integer          :: context = cert_context_stationary
         integer          :: n       = 0
         type(cert_entry) :: e(cert_max_entries)
         ! Cells of the sweep that produced this state whose accepted
         ! composition is not a root of the chemical network (acceptance
         ! class 4 or 6). known = .false. says the caller made no statement.
         integer          :: n_no_chem_root       = 0
         logical          :: chem_root_known      = .false.
         integer          :: first_no_chem_cell   = 0
         ! THE MASS CLOSURE OF THE COMPOSITION, max_j |sum_i f_i A_i - 1|
         ! over the physical cells, and the cell it attains. The mass
         ! fractions are defined by f_i = n_i m_i / rho, so this sum is one
         ! identically and any departure is the arithmetic the state was
         ! built by. REPORTED AND NEVER GATED: it is a statement about the
         ! consistency of the two halves of the state, not an equation the
         ! state has to satisfy, and its natural scale (1e-16) is nine
         ! decades below every tolerance in the inventory. Item L19.
         real*8           :: mass_closure         = 0.0d0
         integer          :: j_mass_closure       = 0
         ! THE FIVE VALIDITY STATES of B1a section 4, in its order. A
         ! negative number is "not produced": no counter exists in the code
         ! for that state, and inventing one here would be reporting a
         ! measurement the run never made.
         ! 4.1 has NO PRODUCER in the code. The one item that had a
         ! cumulative record was the CO thermal ceiling, and it is gone: the
         ! CO row carries published destruction rates, so what was
         ! unvalidated physics is now a closure with a measured domain and
         ! is counted under 4.2. A count is not invented here to replace it.
         integer :: n_active_unvalidated_physics       = -1
         integer :: n_out_of_domain_closure            = -1
         integer :: n_rejected_trial                   = -1
         integer :: n_unbudgeted_accepted_correction   =  0
         integer :: n_specified_external_reservoir     = -1
         ! Cells of the energy update and of the conduction stage that
         ! reached their lower bracket end. Neither stage clamps any more:
         ! both are bracketed and residual-controlled and a floor
         ! activation is a FAILURE that stops the run, not an accepted
         ! state. So these are ATTEMPTS diagnostics and not validity
         ! states, and neither invalidates: a rejected attempt leaves no
         ! contribution to the state being judged (review 2 section 5.3).
         ! The H3+ cooling domain records, reported one by one under the
         ! out-of-domain category so that a count says WHICH edge was left.
         integer :: n_h3p_below_collider               =  0
         integer :: n_h3p_below_fit_T                  =  0
         integer :: n_h3p_above_fit_T                  =  0
         integer :: n_h3p_outside_nonlte_T             =  0
         ! The domain record of the one-sided CO destruction model, reported
         ! beside the H3+ records for the same reason: a count says WHICH
         ! edge of the model's domain the run left. Cell visits summed over
         ! the run's carrier intervals, and the worst ratio of the ordering
         ! tau_dest << tau_res that makes a destruction-only row legitimate.
         integer :: n_co_out_of_domain                 =  0
         integer :: n_co_above_shield_T                =  0
         integer :: n_co_HeII_led                      =  0
         real*8  :: co_worst_dest_over_res             =  0.0d0
         ! The Ly-alpha damping-wing domain record: cell visits below the
         ! (a tau)^(1/3) limit of the slab solution, the visits seen, and
         ! the smallest (a tau)^(1/3) of the run, with the limit itself.
         integer :: n_lya_wing_out                     =  0
         integer :: n_lya_wing_seen                    =  0
         real*8  :: lya_wing_worst                     =  0.0d0
         real*8  :: lya_wing_limit                     =  0.0d0
         ! Whether the run's carrier history is certifiable at all: false
         ! once an interval was left uncovered and the march went on from
         ! the entry carriers (A3b). It is a property of the HISTORY, not
         ! of the state's own rows, and it refuses certification on its own.
         logical :: carrier_history_ok                 = .true.
         integer :: n_energy_floor_attempts            =  0
         integer :: n_conduction_floor_attempts        =  0
         logical :: certified = .false.
         integer :: n_failing = 0
      end type cert_report

      ! WHAT ONE RESIDUAL EVALUATION FOUND ABOUT THE STATE IT WAS GIVEN, as
      ! the probe and trial contexts of contract section 3 read it. It is
      ! filled by the residual evaluation itself (eval_residual,
      ! steady_newton.f90) because that is where the sweep and the packed
      ! rows are; the two decisions below are then functions of this record
      ! alone and are not re-derived at either call site.
      type, public :: cert_evaluation_facts
         integer :: n_sweep_nonfinite     = 0
         integer :: n_no_chem_root_trial  = 0
         integer :: n_no_chem_root_state  = 0
         ! Whether the caller asks for the chemical-root condition at all: a
         ! Jacobian probe is a sample of a directional derivative and is not
         ! a state the run may adopt, so it is not held to it.
         logical :: chem_root_gate        = .true.
         logical :: headroom_ok           = .true.
         logical :: rows_finite           = .true.
      end type cert_evaluation_facts

      ! ---------------------------------------------------------------- !
      ! THE PHYSICAL-STEP VERDICT (context cert_context_physical_step).
      ! Contract section 3 asks the adoption boundary of a marching step
      ! for: the time-discrete balances of every active equation within
      ! their tolerances, admissibility, the element and charge invariants
      ! to round-off, the energy identity of the step, and the
      ! integration-error requirement. The first three are here, because
      ! they are equations and this module owns the equations and their
      ! tolerances. The energy identity and the integration-error estimate
      ! belong to the controller: the identity needs a formation-energy
      ! reservoir this module does not define, and the estimate is a
      ! property of the step-size policy, not of a state.
      !
      ! What it is NOT. This is not a stationary certification, and a
      ! nonzero stationary residual is not a failure of it: a valid
      ! finite-time state of a physical integration has one (review 2, F1).
      integer, parameter, public :: cert_step_ok            = 0
      integer, parameter, public :: cert_step_nonfinite     = 1
      integer, parameter, public :: cert_step_positivity    = 2
      integer, parameter, public :: cert_step_element       = 3
      integer, parameter, public :: cert_step_carrier       = 5
      integer, parameter, public :: cert_step_energy        = 6
      integer, parameter, public :: cert_step_conduction    = 7
      integer, parameter, public :: cert_step_chem_root     = 8
      integer, parameter, public :: cert_step_hydro_row     = 9

      type :: cert_step_verdict
         logical :: accepted  = .true.
         integer :: reason    = cert_step_ok
         ! The operation of the fourteen (b1 section 7.1) whose returned
         ! state refused, so a rejection names a row of that table.
         integer :: operation = 13
         integer :: jworst    = 0
         integer :: kworst    = 0
         real*8  :: measure   = 0.0d0
         real*8  :: tol       = 0.0d0
         ! The worst time-discrete hydrodynamic row of the step, reported
         ! whether or not it refused.
         real*8  :: hydro_row_max = 0.0d0
         integer :: hydro_jworst  = 0
         integer :: hydro_kworst  = 0
         ! Whether that row stands above cert_tol_hydro_step. Reported,
         ! never decisive in the interim: see the tolerance's declaration.
         logical :: hydro_row_above = .false.
         character(len=88) :: what = ''
      end type cert_step_verdict

      public :: certification_evaluate, certification_report_write
      public :: certification_active_equation_count
      public :: certification_row_measure
      public :: certification_row_measure_over
      public :: cert_tol_element_at, cert_tol_carrier_at
      public :: cert_tol_mass_at, mass_row_cell_verdict
      public :: mass_row_column_verdict
      public :: cert_mass_gate_name
      public :: certification_species_row_gate
      public :: certification_entry_index
      ! The stage sum entry, so that src/tests/certification/ can state its
      ! verdict on a measure it chooses rather than on one a solve happens
      ! to produce.
      public :: ionization_stage_sum_entry
      ! What decided the direction of the base contact and whether the
      ! base face flux agrees with it (base_contact_direction_agreement).
      public :: base_contact_direction_agreement
      integer, parameter, public :: contact_direction_from_face_velocity    = 0
      integer, parameter, public :: contact_direction_window_agrees         = 1
      integer, parameter, public :: contact_direction_face_flux_at_rounding = 2
      integer, parameter, public :: contact_direction_window_disagrees      = 3
      public :: certification_last_report, certification_stop_uncertified
      public :: certification_note_stationarity_claim
      public :: trial_state_is_admissible, probe_direction_is_usable
      public :: certification_evaluate_physical_step
      public :: cert_step_verdict
      public :: cert_step_reason_text

      ! The verdict of the LAST stationary certification made, so that the
      ! end of a run can act on the state it actually wrote.
      type(cert_report), save :: cert_last
      ! WHETHER THE RUN CLAIMED A STATIONARY STATE AT ALL. The run states of
      ! docs/a0_run_mode_contract_20260906.md section 2 are three, and a
      ! relaxation snapshot is not one of the other two: a run that ends on a
      ! step cap, on a stall, on a NaN or on any other bound has made no
      ! claim about stationarity, so refusing to certify it is not a finding
      ! about the state -- there was nothing to certify. Only a run that
      ! DECLARED a stationary state and then failed the certification exits
      ! nonzero. Set by the caller that knows which stop fired.
      logical, save :: stationarity_claimed = .false.

      contains

      ! ------------------------------------------------------!

      logical function trial_state_is_admissible(f) result(ok)
      ! THE STATIONARY NEWTON TRIAL CONTEXT (contract section 3, row 2): a
      ! trial is admissible if it is a state at all -- finite rows, a finite
      ! sweep, inside the element headroom -- and if its eliminated chemistry
      ! is locally valid, which here means it leaves no MORE cells without a
      ! chemical root than the iterate it is compared against. The solver's
      ! own merit-decrease and trust-region rules then decide the step; no
      ! stationary tolerance enters this decision (review 2, F1).
      type(cert_evaluation_facts), intent(in) :: f
      ok = (f%n_sweep_nonfinite .eq. 0) .and. f%headroom_ok                &
           .and. f%rows_finite
      if (f%chem_root_gate) ok = ok .and.                                  &
           (f%n_no_chem_root_trial .le. f%n_no_chem_root_state)
      end function trial_state_is_admissible

      ! ------------------------------------------------------!

      logical function probe_direction_is_usable(f) result(ok)
      ! THE PROBE CONTEXT (contract section 3, row 1): a residual sample
      ! taken to build a Jacobian column or a Jacobian-vector product is
      ! usable when every row it produced is a finite number and the sweep
      ! that produced them found none. Convergence is not asked of it, and
      ! neither is the chemical-root condition: a probe is not a state the
      ! run may adopt, and holding it to that condition would refuse the
      ! neighborhood of a point the solver is standing on.
      type(cert_evaluation_facts), intent(in) :: f
      type(cert_evaluation_facts) :: g
      g = f
      g%chem_root_gate = .false.
      ok = trial_state_is_admissible(g)
      end function probe_direction_is_usable

      ! ------------------------------------------------------!

      integer function certification_active_equation_count(rep) result(n)
      ! How many equations of the inventory the configuration makes
      ! independent, whatever their evaluation status.
      type(cert_report), intent(in) :: rep
      integer :: i
      n = 0
      do i = 1, rep%n
         if (rep%e(i)%status .ne. cert_not_applicable) n = n + 1
      enddo
      end function certification_active_equation_count

      ! ------------------------------------------------------!

      function certification_last_report() result(rep)
      ! The last stationary certification made in this run.
      type(cert_report) :: rep
      rep = cert_last
      end function certification_last_report

      ! ------------------------------------------------------!

      subroutine certification_note_stationarity_claim(claimed)
      ! Whether the stop that ended the run DECLARED a stationary state: the
      ! du-threshold stop, the residual gate met while marching, or the
      ! steady finish returning info = 0. A stall, a step cap, a NaN and
      ! "do only PP" are not such stops.
      logical, intent(in) :: claimed
      stationarity_claimed = claimed
      end subroutine certification_note_stationarity_claim

      ! ------------------------------------------------------!

      subroutine certification_stop_uncertified()
      ! A RUN WHOSE FINAL STATE IS NOT CERTIFIED EXITS NONZERO, after its
      ! outputs are written and after the failing entries are named, so that
      ! a script can tell a certified stationary solution from a state that
      ! was merely written (user decision, contract section 8). The status is
      ! 2, chosen to be distinct from the 1 that a Fortran error stop leaves.
      ! It belongs ONLY to a run that declared a stationary state and then
      ! failed the certification; a run that made no such claim exits 0.
      ! The files are complete when this runs: the regression harness
      ! compares file contents and reads its verdict from those, so the
      ! nonzero status changes what a caller sees from the shell and not what
      ! the harness compares.
      if (cert_last%certified) return
      if (.not. stationarity_claimed) then
         write(*,'(A)') ' (certification) the state written is a '//       &
              'relaxation snapshot: the run made no stationary claim, so'
         write(*,'(A)') '   it is written certified=F with the reason '//  &
              '"no stationary claim" and the run exits 0.'
         return
      endif
      call certification_report_write(cert_last,                           &
           'FINAL STATE AS WRITTEN -- this run is NOT a certified '//      &
           'stationary solution')
      write(*,'(A)') ' (certification) the state was written in full; '//  &
           'the run exits with status 2 because it is not certified.'
      flush(6)
      stop 2
      end subroutine certification_stop_uncertified

      ! ------------------------------------------------------!

      subroutine certification_row_measure(nc, jsplit, res, scale, rmax,  &
                                           jworst, ok_finite, wvol, rvol)
      ! THE MEASURE OF ONE ROW OVER ONE COLUMN OF CELLS, and the single
      ! definition of it.
      !
      !   rmax = max_j |res_j| / max(scale_j, floor)
      !
      ! the maximum over cells, never a volume average: a front or a thin
      ! layer must not be paid for by the volume of the wind. rvol is the
      ! volume-weighted companion, formed separately over the cells below
      ! jsplit and at or above it and combined by the larger of the two
      ! (residual_norms gives the reason: the wind carries almost all of
      ! sum r^2 dr, so a single sum averages the inner column away). It is
      ! REPORTED and is never the condition -- a species whose rates are
      ! small beside a balanced dominant one can be out by any factor of its
      ! own terms and move a ratio of sums by nothing at all (review R7).
      !
      ! THE COMPANION IS OPTIONAL, and so is the cell volume it needs: a
      ! caller that reads only rmax -- the convergence gate of the Newton
      ! solve is one, and it asks for this measure at every residual
      ! evaluation -- asks for neither and the two sums are not formed.
      ! rmax, jworst and ok_finite do not depend on them, so the number a
      ! caller reads is the same either way.
      integer, intent(in)  :: nc, jsplit
      real*8,  intent(in)  :: res(nc), scale(nc)
      real*8,  intent(out) :: rmax
      integer, intent(out) :: jworst
      logical, intent(out) :: ok_finite
      real*8,  intent(in),  optional :: wvol(nc)
      real*8,  intent(out), optional :: rvol
      real*8  :: q, num_l, den_l, num_w, den_w
      integer :: j
      logical :: want_vol
      want_vol = present(rvol) .and. present(wvol)
      rmax = 0.0d0;  jworst = 0;  ok_finite = .true.
      num_l = 0.0d0;  den_l = 0.0d0;  num_w = 0.0d0;  den_w = 0.0d0
      do j = 1, nc
         if (.not. finite_real(res(j)) .or. .not. finite_real(scale(j)))  &
            ok_finite = .false.
         q = abs(res(j))/max(scale(j), cert_scale_floor)
         if (q .gt. rmax) then
            rmax = q;  jworst = j
         endif
         if (.not. want_vol) cycle
         if (j .lt. jsplit) then
            num_l = num_l + abs(res(j))*wvol(j)
            den_l = den_l + scale(j)*wvol(j)
         else
            num_w = num_w + abs(res(j))*wvol(j)
            den_w = den_w + scale(j)*wvol(j)
         endif
      enddo
      ! A companion asked for without the cell volume to weight it with is
      ! not a measure of anything, so it is returned as zero rather than
      ! as an unweighted sum that would read like one.
      if (present(rvol)) then
         if (want_vol) then
            rvol = max(num_l/max(den_l, cert_scale_floor),                &
                       num_w/max(den_w, cert_scale_floor))
         else
            rvol = 0.0d0
         endif
      endif
      end subroutine certification_row_measure

      ! ------------------------------------------------------!

      double precision function cert_tol_mass_at(j, u) result(tol)
      ! THE CONTINUITY ROW'S TOLERANCE AT ONE CELL, and the only reader of
      ! cert_tol_mass: the fixed value where the cell's arithmetic permits
      ! it, and the cell's own rounding floor times the measured margin
      ! where it does not.
      !
      !   tol(j) = max( cert_tol_mass,
      !                 min( ceiling, c_round * floor(j) ) )
      !
      ! A state is therefore never asked to resolve its mass row below the
      ! rounding of its own base layer, and never allowed to do worse than
      ! 3e-12 where the rounding is below it. The declarations of
      ! cert_mass_round_margin and of mass_row_rounding_floor
      ! (steady_residual.f90) carry the derivation and the measurements.
      !
      ! WHAT IS GIVEN UP, stated because it is the cost of the anchor. Where
      ! the floor binds, the mass row of that cell is judged at eps/Mach of
      ! its own flux instead of at 3e-12: on the HD 209458 b base layer that
      ! is about 3e-10 rather than 3e-12. What the row is FOR is undamaged --
      ! the flux errors it was given this scale to catch are fractions of a
      ! percent to tens of percent of the flux (the 30 percent at 1.03 R_p
      ! of docs/p54_base_layer_mass_flux.md), seven decades above any floor
      ! measured -- and the mass-flux spread gate reads the same quantity
      ! independently.
      integer,                        intent(in) :: j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8  :: dist
      logical :: within, anchored
      call mass_row_cell_verdict(0.0d0, mass_row_rounding_floor(j, u),    &
                                 tol, dist, within, anchored)
      end function cert_tol_mass_at

      ! ------------------------------------------------------!

      pure subroutine mass_row_cell_verdict(q, floor_q, tol, dist,        &
                                            within, anchored, status)
      ! THE VERDICT ON ONE CELL OF THE CONTINUITY ROW, and the single
      ! expression of it: the tolerance that applies there, the distance
      ! from it, which of the two gates set it, and (optionally) whether
      ! the cell could be judged at all.
      !
      !   raw = c_round * floor
      !
      ! UNRESOLVED: raw >= ceiling, or floor_q is not a finite number. The
      ! comment beside cert_tol_mass_ceiling states the reason: an estimated
      ! floor that reaches the ceiling is a row whose balance the assembly
      ! cannot resolve, not one that passes. Such a cell never reads
      ! `within`, whatever q is; its distance is max(1, q/ceiling), so a
      ! reader that ranks cells by distance still sees a refusal rather than
      ! a zero.
      !
      ! RESOLVED: raw < ceiling. The rule is the one anchored by item P16,
      !
      !   tol = max( cert_tol_mass, raw )
      !
      ! and `anchored` says whether the fixed floor or the rounding anchor
      ! set it.
      !
      ! q is the cell's row measure |R_1|/s_1 and floor_q the rounding
      ! floor OF THAT SAME MEASURE (mass_row_rounding_floor). The rule is
      ! separated from the column it is applied over so that it can be
      ! stated on one (q, floor) pair and nothing re-derives it.
      real*8,  intent(in)  :: q, floor_q
      real*8,  intent(out) :: tol, dist
      logical, intent(out) :: within, anchored
      integer, intent(out), optional :: status
      real*8  :: raw
      logical :: floor_is_not_finite
      ! A NaN floor fails every ordered comparison, including .ge. against
      ! the ceiling, so it is caught by the one comparison an IEEE NaN
      ! never satisfies: equality with itself.
      floor_is_not_finite = .not. (floor_q .eq. floor_q)
      raw = cert_mass_round_margin*floor_q
      if (floor_is_not_finite .or. raw .ge. cert_tol_mass_ceiling) then
         tol      = cert_tol_mass_ceiling
         anchored = .true.
         dist     = max(1.0d0, q/tol)
         within   = .false.
         if (present(status)) status = cert_mass_row_unresolved
      else
         tol      = max(cert_tol_mass, raw)
         anchored = (tol .gt. cert_tol_mass)
         dist     = q/tol
         within   = (dist .lt. 1.0d0)
         if (present(status)) status = cert_mass_row_resolved
      endif
      end subroutine mass_row_cell_verdict

      ! ------------------------------------------------------!

      pure function cert_mass_gate_name(anchored) result(nm)
      ! WHICH OF THE TWO GATES SET THE CONTINUITY ROW'S TOLERANCE AT A
      ! CELL, in the words the report prints, so that the report and a
      ! reader of the verdict take the name from one place.
      logical, intent(in) :: anchored
      character(len=32)   :: nm
      if (anchored) then
         nm = 'the rounding anchor of that cell'
      else
         nm = 'the fixed tolerance'
      endif
      end function cert_mass_gate_name

      ! ------------------------------------------------------!

      pure function cert_tol_element_at(rr) result(tol)
      ! THE ELEMENTAL TRANSPORT TOLERANCE AT ONE RADIUS, and the only
      ! reader of cert_tol_element_wind: what is a gate in the wind is a
      ! reported measure below it, for the reasons stated at the
      ! declaration of the two radii.
      real*8, intent(in) :: rr
      real*8 :: tol
      tol = cert_tol_reported_only
      if (rr .ge. cert_regime_wind_r) tol = cert_tol_element_wind
      end function cert_tol_element_at

      ! ------------------------------------------------------!

      pure function cert_tol_carrier_at(rr) result(tol)
      ! The same for a transported carrier row. The two are separate
      ! functions because they are tolerances of two different equations
      ! and nothing but the anchoring makes them equal today.
      real*8, intent(in) :: rr
      real*8 :: tol
      tol = cert_tol_reported_only
      if (rr .ge. cert_regime_wind_r) tol = cert_tol_carrier_wind
      end function cert_tol_carrier_at

      ! ------------------------------------------------------!

      subroutine certification_species_row_gate(nc, is_carrier, res,      &
                                                scale, dgate, rgate,      &
                                                jgate, absent)
      ! THE PART OF A SPECIES ROW THAT DECIDES, in one reduction and one
      ! place: over the physical cells, the largest of
      !
      !     |res_j| / max(scale_j, floor) / tol(r_j)
      !
      ! with tol from the accessor of this row's kind. dgate is that
      ! number, which is below 1 exactly when the row is within tolerance
      ! everywhere it is gated; rgate is the raw measure of the cell that
      ! carries it and jgate is that cell. A cell the accessor does not
      ! gate divides by cert_tol_reported_only and contributes nothing.
      !
      ! jgate = 0 says NO CELL OF THE COLUMN IS GATED, which is not a
      ! satisfied row: the caller reports it as a row the state cannot be
      ! judged on rather than as one it passed.
      integer, intent(in)  :: nc
      logical, intent(in)  :: is_carrier
      real*8,  intent(in)  :: res(nc), scale(nc)
      real*8,  intent(out) :: dgate, rgate
      integer, intent(out) :: jgate
      ! THE SECOND WAY A CELL LEAVES THE GATE: the species of this row is
      ! absolutely absent there. The radius accessor above excludes a cell
      ! because the DISCRETIZATION cannot be judged there; this excludes it
      ! because there is no species to judge. Both leave the cell reported.
      logical, optional, intent(in) :: absent(nc)
      real*8  :: q, tol, d
      integer :: j
      dgate = 0.0d0;  rgate = 0.0d0;  jgate = 0
      do j = 1, nc
         if (present(absent)) then
            if (absent(j)) cycle
         endif
         if (is_carrier) then
            tol = cert_tol_carrier_at(r(j))
         else
            tol = cert_tol_element_at(r(j))
         endif
         if (tol .ge. cert_tol_reported_only) cycle
         q = abs(res(j))/max(scale(j), cert_scale_floor)
         d = q/tol
         if (jgate .eq. 0 .or. d .gt. dgate) then
            dgate = d;  rgate = q;  jgate = j
         endif
      enddo
      end subroutine certification_species_row_gate

      ! ------------------------------------------------------!

      subroutine certification_row_measure_over(nc, res, scale, rlo, rhi,  &
                                                rmax, jworst)
      ! THE SAME MEASURE OVER A RADIAL WINDOW, rlo <= r(j) < rhi.
      !
      ! Same definition as certification_row_measure, restricted to the
      ! cells of one regime; a window that holds no cell returns zero with
      ! jworst = 0, which the report prints as "no cell" rather than as a
      ! satisfied row. It is a REPORTING reduction: nothing in the verdict
      ! reads it.
      integer, intent(in)  :: nc
      real*8,  intent(in)  :: res(nc), scale(nc), rlo, rhi
      real*8,  intent(out) :: rmax
      integer, intent(out) :: jworst
      real*8  :: q
      integer :: j
      rmax = 0.0d0;  jworst = 0
      do j = 1, nc
         if (r(j) .lt. rlo .or. r(j) .ge. rhi) cycle
         q = abs(res(j))/max(scale(j), cert_scale_floor)
         if (jworst .eq. 0 .or. q .gt. rmax) then
            rmax = q;  jworst = j
         endif
      enddo
      end subroutine certification_row_measure_over

      ! ------------------------------------------------------!

      subroutine fill_regime_measures(ent, res, scale)
      ! The regime split and the probe-radius reading of one entry, taken
      ! from the same res and scale the entry's own measure was taken from.
      type(cert_entry), intent(inout) :: ent
      real*8,           intent(in)    :: res(1:N), scale(1:N)
      real*8  :: rp, dbest, dj
      integer :: j
      call certification_row_measure_over(N, res, scale, -1.0d0,           &
               cert_regime_layer_r, ent%row_max_layer, ent%jworst_layer)
      call certification_row_measure_over(N, res, scale, -1.0d0,           &
               cert_regime_wind_r, ent%row_max_reported,                   &
               ent%jworst_reported)
      call certification_row_measure_over(N, res, scale,                   &
               cert_regime_wind_r, huge(1.0d0), ent%row_max_wind,          &
               ent%jworst_wind)
      ent%row_at_probe = 0.0d0
      ent%jprobe       = 0
      rp = anchor_probe_radius()
      if (rp .le. 0.0d0) return
      dbest = huge(1.0d0)
      do j = 1, N
         dj = abs(r(j) - rp)
         if (dj .lt. dbest) then
            dbest = dj;  ent%jprobe = j
         endif
      enddo
      if (ent%jprobe .gt. 0) ent%row_at_probe =                            &
           abs(res(ent%jprobe))/max(scale(ent%jprobe), cert_scale_floor)
      end subroutine fill_regime_measures

      ! ------------------------------------------------------!

      subroutine gate_species_row(ent, is_carrier, res, scale, absent)
      ! THE VERDICT OF A ROW WHOSE TOLERANCE IS A FUNCTION OF THE RADIUS.
      ! The entry keeps its whole-column measure, and what decides is the
      ! gated part of it (certification_species_row_gate). A column with
      ! no gated cell cannot be judged on this row and is refused with
      ! that as the reason, so a domain that ends below the wind radius
      ! never reads as a satisfied balance.
      type(cert_entry), intent(inout) :: ent
      logical,          intent(in)    :: is_carrier
      real*8,           intent(in)    :: res(1:N), scale(1:N)
      logical, optional, intent(in)   :: absent(1:N)
      real*8  :: dgate
      if (present(absent)) then
         call certification_species_row_gate(N, is_carrier, res, scale,   &
                  dgate, ent%row_max_gate, ent%jworst_gate,               &
                  absent = absent)
      else
         call certification_species_row_gate(N, is_carrier, res, scale,   &
                  dgate, ent%row_max_gate, ent%jworst_gate)
      endif
      if (ent%jworst_gate .eq. 0) then
         ent%status = cert_unavailable
         if (ent%n_absent .gt. 0) then
            ent%reason = 'every cell at or above the wind radius is '//   &
                         'empty of this species'
         else
            ent%reason = 'no cell at or above the wind radius to '//      &
                         'gate this row'
         endif
         ent%within_tol = .false.
         return
      endif
      ent%within_tol = ent%finite .and. (dgate .lt. 1.0d0)
      end subroutine gate_species_row

      ! ------------------------------------------------------!

      function certification_entry_index(rep, name) result(idx)
      ! The index of the entry of that name, 0 if the inventory has none.
      ! An accessor: a caller that wants one row's numbers should not have
      ! to walk the array and match the padded name itself.
      type(cert_report), intent(in) :: rep
      character(len=*),  intent(in) :: name
      integer :: idx, i
      idx = 0
      do i = 1, rep%n
         if (trim(rep%e(i)%name) .eq. trim(name)) then
            idx = i;  return
         endif
      enddo
      end function certification_entry_index

      ! ------------------------------------------------------!

      function anchor_report_on() result(on)
      ! Whether the run asks for the anchoring block (EXHALE_CERT_ANCHOR=1):
      ! the row measures at full double precision, split by regime. It is a
      ! second reading of numbers the block above already prints, at the
      ! precision a reproducibility or a grid comparison needs, and it is
      ! off by default so that a run's report is unchanged by it.
      logical :: on
      character(len=32) :: env
      call get_environment_variable('EXHALE_CERT_ANCHOR', env)
      on = (trim(env) .eq. '1')
      end function anchor_report_on

      ! ------------------------------------------------------!

      function anchor_probe_radius() result(rp)
      ! The radius EXHALE_CERT_ANCHOR_R names, at which every row is read
      ! in addition to its maxima, so that one physical radius can be
      ! followed across two grids. Zero or unset: no probe.
      real*8 :: rp
      integer :: ios
      character(len=32) :: env
      rp = 0.0d0
      call get_environment_variable('EXHALE_CERT_ANCHOR_R', env)
      if (len_trim(env) .eq. 0) return
      read(env, *, iostat = ios) rp
      if (ios .ne. 0) rp = 0.0d0
      end function anchor_probe_radius

      ! ------------------------------------------------------!

      subroutine certification_evaluate(context, u, Res, f_sp,            &
                                        resid_tol, n_no_chem_root,        &
                                        chem_root_known, rep)
      ! THE EVALUATOR. Builds the active-equation inventory from the
      ! configuration flags, measures what can be measured on the state
      ! given, and takes the verdict.
      !
      ! The state is passed in and is never taken from a cache: u and Res are
      ! the conserved variables and the assembled residual of the SAME state,
      ! and f_sp is its composition. The density and velocity the carrier
      ! balance needs are taken from u itself rather than from the caller's
      ! copies, so the two halves of the measurement cannot describe two
      ! states. The caller supplies Res because it has just assembled it;
      ! assembling a second one here would answer for a state built by this
      ! routine rather than for the state that is about to be written.
      integer,                        intent(in)  :: context
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: Res
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,                         intent(in)  :: resid_tol
      integer,                        intent(in)  :: n_no_chem_root
      logical,                        intent(in)  :: chem_root_known
      type(cert_report),              intent(out) :: rep

      integer :: k, ic, i, j, itr, imel
      real*8, dimension(1:N,n_carrier_max) :: cres, cterms
      logical, dimension(1:N,n_carrier_max) :: cabsent
      real*8, dimension(3,1-Ng:N+Ng)       :: Wcert
      real*8, dimension(1-Ng:N+Ng) :: Tcert, cn_hi, cn_hii, cn_hei
      real*8, dimension(1-Ng:N+Ng) :: cn_heii, cn_heiii, cn_heitr
      real*8, dimension(1-Ng:N+Ng) :: cn_e, cn_tot
      real*8, dimension(1-Ng:N+Ng) :: Frho_cert
      real*8, dimension(1-Ng:N+Ng,n_mion) :: cn_m
      real*8, dimension(1:N) :: ehe_res, ehe_sc
      real*8, dimension(1:N,n_melem) :: etr_res, etr_sc
      logical, dimension(n_melem)    :: etr_carried
      real*8, dimension(1:N) :: clo_res, clo_tr, n2_res, n2_sc
      logical :: carriers_measured, he_ok, tr_ok, clo_ok, n2_ok
      ! The stage sum identity of this state, its face and whether it was
      ! measured at all (ionization_stage_sum_measure), ONE PER ELEMENT
      ! whose stages are carried.
      real*8  :: ssum_max, ssum_g
      integer :: ssum_j, iestg
      ! The index add_entry returns. The short carrier entries of the loop
      ! below write no field through it; it is taken because the index is
      ! the only thing that says which entry a call created.
      integer :: ient
      logical :: ssum_known
      character(len=76) :: clo_why, n2_why

      rep%context         = context
      rep%n               = 0
      rep%n_no_chem_root  = n_no_chem_root
      rep%chem_root_known = chem_root_known
      ! resid_tol is the run's own "Resid tol". It is the SOLVER's target
      ! and is no longer the certification's: the three hydrodynamic rows
      ! carry the fixed numbers of the convergence study (their
      ! declarations say why). The argument is kept because the caller
      ! passes what it asked its solver for and the report names both.
      cert_resid_tol_of_run = resid_tol

      ! ---- the primitive state, and the temperature it implies ----
      ! Derived from u and f_sp themselves -- the pressure of u through the
      ! caloric EOS, the particle and electron counts of f_sp -- and not
      ! from any cached array, so every measurement below describes ONE
      ! state.
      cn_hei = 0.0d0;  cn_heii = 0.0d0;  cn_heiii = 0.0d0
      cn_heitr = 0.0d0
      call U_to_W(u, Wcert)
      call get_species_densities(Wcert(1,:), f_sp, cn_hi, cn_hii, cn_hei,  &
                                 cn_heii, cn_heiii, cn_heitr, cn_m,        &
                                 cn_e, cn_tot)
      call comp_T_from_p(Wcert(3,:), cn_tot, cn_e, Tcert)

      ! ---- the mass closure of the composition (reported, item L19) ----
      ! sum_i f_i A_i - 1, from the same (u, f_sp) every row above is
      ! measured on: calc_rho weighs the species of this state with the
      ! run's own mass policy, and the density it returns is rho times that
      ! sum. Reported with its worst cell and gating nothing.
      call mass_closure_of_state(Wcert(1,:), f_sp, rep%mass_closure,      &
                                 rep%j_mass_closure)

      ! ---- the hydrodynamic rows: always present (B1a section 2.1) ----
      call hydro_row_entry(rep, 1, 'hydrodynamic mass row',      u, Res)
      call hydro_row_entry(rep, 2, 'hydrodynamic momentum row',  u, Res)
      call hydro_row_entry(rep, 3, 'hydrodynamic energy row',    u, Res)

      ! ---- the transported carrier balances (B1a section 2.2) ----
      ! The set is carrier_set_init's, fixed once after the keys are parsed;
      ! a carrier the run does not solve carries no equation.
      carriers_measured = .false.
      cabsent = .false.
      if (transported_rows_exist() .and. bg_ready) then
         call carrier_rows_of_state(Wcert(1,:), Wcert(2,:), f_sp,         &
                                    cres, cterms, cabsent)
         carriers_measured = .true.
      endif
      do ic = 1, n_carrier_max
         if (ic .gt. n_carrier .or. .not. carrier_solved(ic)               &
             .or. .not. transported_rows_exist()) then
            call add_entry(rep, 'carrier balance '//trim(carrier_name(ic)),&
                           cert_not_applicable,                           &
                           'the run does not carry this species', ient)
            cycle
         endif
         if (.not. carriers_measured) then
            call add_entry(rep, 'carrier balance '//trim(carrier_name(ic)),&
                           cert_unavailable,                              &
                           'the frozen background is not ready (bg_ready)',&
                           ient)
            cycle
         endif
         call carrier_row_entry(rep, ic, cres, cterms, cabsent)
      enddo

      ! ---- the ionization stage sum: moving charge moves no nucleus ----
      ! With the ionization state transported, the stage fluxes of an
      ! element must sum, face by face and over ALL its stages, to that
      ! element's own nucleus flux (ionization_stage_transport, equation
      ! 2). The identity is structural, so what stands in it is the
      ! rounding of the sums, and a measure above that says the
      ! construction has been broken.
      !
      ! IT GATES, at the floating-point bound of the identity: the tolerance
      ! and the verdict are ionization_stage_sum_entry's, one entry per
      ! element, and the bound is read from the one place that forms it
      ! (ionization_stage_sum_rounding_bound). Where no stage of an element
      ! is transported at all there is no entry.
      if (ionization_transport .and. carriers_measured) then
         do iestg = ien_H, ien_He
            call ionization_stage_sum_measure(iestg, ssum_max, ssum_j,    &
                                              ssum_known, ssum_g)
            call ionization_stage_sum_entry(rep, iestg, ssum_max, ssum_j, &
                                            ssum_known, ssum_g)
         enddo
      endif

      ! ---- the elemental transport balances (B1a section 2.3) ----
      ! Stationary: transport against transport, the time term removed the
      ! way the carrier balance removes it. Measured on this state's own
      ! composition inside the operator's own coefficients, with the five
      ! diagnostics of the last solve held aside and put back.
      he_ok = .false.;  tr_ok = .false.;  etr_carried = .false.
      if (he_diffusion .and. thereis_He) then
         ! The wind the elemental rows ride on: the face mass flux of the
         ! state being certified, read from the mass row the caller
         ! assembled for it and handed to the element operator, which stands
         ! below the steady residual and does not reach up to it.
         call face_mass_flux_of_state(Wcert(1,:), Frho_cert)
         call element_transport_residual(Wcert(1,:), Tcert,                &
                  f_sp, Frho_cert, ehe_res, ehe_sc, he_ok, etr_res,        &
                  etr_sc, tr_ok, tr_carried = etr_carried)
      endif
      call transport_row_entry(rep, 'elemental transport He/H partition',   &
               he_diffusion .and. thereis_He, he_ok, ehe_res, ehe_sc,       &
               'sum of the row''s own terms [g cm^-3 s^-1]; floor '//      &
               '1e-20 rho X_base/(R0/v0)')
      ! ONE ENTRY PER TRACE ELEMENT.  The balance of each element is its
      ! own equation of the inventory, and the reduction that used to
      ! report the worst element at each cell could not say which element
      ! was out, nor be a row of a stationary system.
      do imel = 1, n_melem
         call transport_row_entry(rep, 'elemental transport '//            &
                  trim(melem_name(imel)),                                  &
                  he_diffusion .and. he_metal_diffusion .and.              &
                  thereis_metals .and. etr_carried(imel),                  &
                  tr_ok, etr_res(:,imel), etr_sc(:,imel),                  &
                  'sum of the row''s own terms [1/s]; '//                 &
                  'floor 1e-20 f_base/(R0/v0)')
      enddo

      ! ---- the eliminated-species closure and the He 2^3S level row ----
      ! One evaluation answers both (B1a sections 2.4 and 2.5): the
      ! metastable is an unknown of the same coupled cell system, so its
      ! balance is one row of the vector the closure measure takes the
      ! maximum of. The closure entry is that maximum, which is exactly the
      ! measure the sweep judges a root by, so a cell of acceptance class
      ! 1, 2, 3 or 5 reads at or below ieq_res_tol and a class-4 or class-6 cell above
      ! it. The triplet row is reported separately as well, so a refusal can
      ! say whether it is the level that is out.
      call ionization_closure_residual_profile(Wcert(1,:), f_sp, clo_res,  &
                                               clo_tr, clo_ok, clo_why)
      itr = ieq_triplet_row()
      call unit_scale_entry(rep, 'level balance He 2^3S',                   &
               thereis_He .and. thereis_HeITR, clo_ok .and. (itr .gt. 0),   &
               clo_tr, cert_tol_level, clo_why,                             &
               'the row''s own turnover rate; dimensionless, scale 1')

      ! ---- the H(n=2) level balance (B1a section 2.5) ----
      call excited_hydrogen_level_residual(Tcert, Wcert(1,:), f_sp,        &
                                           n2_res, n2_sc, n2_ok, n2_why)
      call transport_row_entry(rep, 'level balance H(n=2)', use_excited_H,  &
               n2_ok, n2_res, n2_sc,                                       &
               'sum of the row''s own terms [cm^-3 s^-1], worst of '//    &
               '2s/2p; floor 1e-30 cm^-3 s^-1',                            &
               tol_whole_column = cert_tol_level)

      ! ---- the eliminated-species closure, one entry per system variant ----
      ! (B1a section 2.4). Exactly one of these is active in a run: the
      ! variant the configuration selects. Its residual is re-evaluated for
      ! acceptance INSIDE the sweep that also rewrites f_sp and rho, and
      ! there is no entry point that takes a state and returns the residual
      ! of its eliminated species without solving, so the active variant is
      ! unavailable and every run therefore stands uncertified until the
      ! evaluator of contract step 4 exists.
      do i = 1, n_system_variant
         call unit_scale_entry(rep,                                        &
              'eliminated-species closure '//trim(system_variant_name(i)),  &
              i .eq. active_system_variant(), clo_ok, clo_res,              &
              cert_tol_closure, clo_why,                                    &
              'the row''s own turnover rate; dimensionless, scale 1')
      enddo

      ! ---- the validity states (B1a section 4) ----
      call read_validity_states(rep)

      ! THE CARRIER HISTORY. An interval the carrier retry could not cover,
      ! from which the march went on with the entry carriers, means the
      ! transported composition of this run is not the one the equation
      ! asked for at that step. That is a property of the HISTORY and not
      ! of the state's own rows, so no row measure can see it and it is
      ! read here (A3b, carrier_history_certifiable).
      rep%carrier_history_ok = .true.
      if (transported_rows_exist())                                        &
         rep%carrier_history_ok = carrier_history_certifiable()

      ! ---- the verdict ----
      rep%n_failing = 0
      rep%certified = .true.
      do i = 1, rep%n
         select case (rep%e(i)%status)
            case (cert_unavailable)
               rep%certified = .false.
               rep%n_failing = rep%n_failing + 1
            case (cert_evaluated)
               if (.not. rep%e(i)%within_tol .or. .not. rep%e(i)%finite)   &
                  then
                  rep%certified = .false.
                  rep%n_failing = rep%n_failing + 1
               endif
         end select
      enddo
      if (.not. chem_root_known) rep%certified = .false.
      if (n_no_chem_root .ne. 0) rep%certified = .false.
      ! An unbudgeted accepted correction is a change of a conserved
      ! quantity with no source term behind it, and it stands in the history
      ! of the state being judged (B1a section 4.4). The energy and
      ! conduction temperature floors are the two of them that have a
      ! counter today.
      if (rep%n_unbudgeted_accepted_correction .gt. 0)                     &
         rep%certified = .false.
      ! A marked carrier history refuses certification with that as the
      ! reason: a state reached through an interval that was never covered
      ! is not a state of the equations the run states it solves.
      if (.not. rep%carrier_history_ok) rep%certified = .false.

      ! The first cell of the state whose composition is not a root, so a
      ! refusal names a cell and not only a count.
      rep%first_no_chem_cell = 0
      if (n_no_chem_root .gt. 0 .and. allocated(ieq_nonroot_streak)) then
         do j = 1, N
            if (ieq_nonroot_streak(j) .gt. 0) then
               rep%first_no_chem_cell = j;  exit
            endif
         enddo
      endif

      if (context .eq. cert_context_stationary) then
         if (rep%certified) then
            ! WHAT THE FILE SAYS ABOUT A STATE CERTIFIED WITH SPECIES
            ! ROWS: those rows were judged in the wind alone, so the
            ! header carries that qualification and a reader of the file
            ! knows it without the run log. A state with no species row
            ! is certified over its whole column and carries no token.
            call set_state_certified(.true.,                              &
                 trim(wind_certification_token(rep)))
         else if (.not. stationarity_claimed) then
            call set_state_certified(.false., 'no_stationary_claim')
         else
            call set_state_certified(.false., 'failing_entries')
         endif
         cert_last = rep
      endif

      end subroutine certification_evaluate

      ! ------------------------------------------------------!

      subroutine carrier_rows_of_state(rho, v, f_sp, cres, cterms, cabsent)
      ! The carrier balance of a state, measured INSIDE AN ISOLATED
      ! WORKSPACE. carrier_steady_residual refreshes the frozen background,
      ! the photolysis rates, the advection correction, the row terms, the
      ! two H2 scales and the element headroom; all of that is put back, so
      ! the next transport step and the next output read exactly what they
      ! would have read had the measurement not been made.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, v
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1:N,n_carrier_max),   intent(out) :: cres, cterms
      ! Where each carrier is absolutely absent, measured inside the same
      ! isolated evaluation that forms the rows, so the mark and the row
      ! describe one state (carrier_absent_fraction, carrier_residual).
      logical, dimension(1:N,n_carrier_max),  intent(out) :: cabsent
      type(carrier_module_state) :: ws
      real*8  :: rcmax, rvol
      integer :: jworst, icworst
      call save_carrier_module_state(ws)
      call carrier_steady_residual(rho, v, f_sp, rcmax, jworst, icworst,  &
                                   rvol=rvol, res_out=cres,              &
                                   terms_out=cterms, absent_out=cabsent)
      ! The terms of the rows just measured, on request
      ! (EXHALE_CARRIER_ROW_TERMS=1, default off): the record belongs to
      ! THIS evaluation, so it is written here and not left for a later
      ! assembly to overwrite.  It is a file, not module state, so it is
      ! outside the round trip asserted below.
      if (carrier_row_terms_on())                                        &
         call carrier_row_terms_write('output/carrier_row_terms.txt')
      call restore_carrier_module_state(ws)
      ! The round trip is ASSERTED, not assumed: a measurement that changed
      ! the state it measured would move the next transport step, and the
      ! next output, by an amount nothing else in the run would explain.
      if (.not. carrier_module_state_matches(ws))                          &
         write(*,'(A)') ' (certification) WARNING: the carrier module '//  &
              'state was NOT reinstated after the isolated evaluation'
      end subroutine carrier_rows_of_state

      ! ------------------------------------------------------!

      subroutine hydro_row_entry(rep, k, name, u, Res)
      ! Row k of the hydrodynamic system, measured on the state (u, Res)
      ! with the scale the residual code defines: the largest term the row
      ! itself contains (residual_row_scale, steady_residual.f90).
      !
      ! The scale is read from the terms store_row_terms left behind, so it
      ! is a number about a state; unless those terms belong to THIS u the
      ! row is unavailable and is not given a value.
      type(cert_report), intent(inout) :: rep
      integer,           intent(in)    :: k
      character(len=*),  intent(in)    :: name
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, Res
      real*8, dimension(1:N) :: rr, ss, ww
      integer :: j, idx
      if (.not. row_terms_describe_state(u)) then
         call add_entry(rep, name, cert_unavailable,                       &
              'no residual terms assembled from this state', idx)
         return
      endif
      call add_entry(rep, name, cert_evaluated, '', idx)
      select case (k)
         case (1);    rep%e(idx)%tol = cert_tol_mass
         case (2);    rep%e(idx)%tol = cert_tol_momentum
         case default;rep%e(idx)%tol = cert_tol_energy
      end select
      rep%e(idx)%units_floor = 'residual_row_scale, the row''s largest '// &
                              'term; floor 1e-300 (numerical only)'
      do j = 1, N
         rr(j) = Res(k,j)
         ss(j) = residual_row_scale(k, j, u)
         ww(j) = r(j)*r(j)*dr_j(j)
      enddo
      call certification_row_measure(N, j_min, rr, ss,                    &
                                     rep%e(idx)%row_max,                  &
                                     rep%e(idx)%jworst,                   &
                                     rep%e(idx)%finite,                   &
                                     wvol = ww, rvol = rep%e(idx)%row_vol)
      call fill_regime_measures(rep%e(idx), rr, ss)
      if (k .eq. 1) then
         call mass_row_verdict(rep%e(idx), rr, ss, u)
      else
         rep%e(idx)%within_tol = rep%e(idx)%finite .and.                   &
                                 (rep%e(idx)%row_max .lt. rep%e(idx)%tol)
      endif
      end subroutine hydro_row_entry

      ! ------------------------------------------------------!

      subroutine mass_row_verdict(ent, rr, ss, u)
      ! THE CONTINUITY ROW'S VERDICT ON A PHYSICAL STATE: the rounding floor
      ! of every cell, read off u, handed to the aggregator below. Nothing
      ! here decides the verdict; mass_row_column_verdict is the one
      ! expression of it, and is exercised directly by the test suite on
      ! synthetic floors so that a verdict can be asserted without also
      ! standing up a flux reconstruction on a physical state.
      type(cert_entry), intent(inout) :: ent
      real*8, dimension(1:N),                 intent(in) :: rr, ss
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u
      real*8, dimension(N) :: floors
      integer :: j
      do j = 1, N
         floors(j) = mass_row_rounding_floor(j, u)
      enddo
      call mass_row_column_verdict(N, rr, ss, floors, ent)
      end subroutine mass_row_verdict

      ! ------------------------------------------------------!

      pure subroutine mass_row_column_verdict(nc, rr, ss, floors, ent)
      ! THE CONTINUITY ROW'S VERDICT, TAKEN CELL BY CELL, because its
      ! tolerance is a function of the cell (cert_tol_mass_at). The cell
      ! that binds is the one whose measure stands furthest outside its own
      ! tolerance, which is not in general the cell of the largest measure:
      ! a wind cell at 1e-11 against 3e-12 refuses while a base cell at
      ! 1e-9 against its own rounding floor does not.
      !
      ! row_max and jworst are left as they are. They are the measure of
      ! the row over the column and the cell that carries it, one
      ! expression for every row of the report, and the acceptance gate is
      ! asserted against them (gate_equals_certification, steady_newton.f90).
      ! What this routine decides is the verdict and which number took it.
      !
      ! REVIEW R5 (item Q3): a cell whose floor makes the row UNRESOLVED
      ! (mass_row_cell_verdict) is counted and the first one is named
      ! (n_mass_unresolved, j_first_mass_unresolved), whether or not it is
      ! also the binding cell.
      integer,                 intent(in)    :: nc
      real*8,  dimension(nc),  intent(in)    :: rr, ss, floors
      type(cert_entry),        intent(inout) :: ent
      real*8  :: q, tt, dd
      logical :: wi, an, wi_bind
      integer :: j, st
      ent%jbind = 0;  ent%dist_bind = 0.0d0;  ent%row_at_bind = 0.0d0
      ent%n_mass_unresolved = 0;  ent%j_first_mass_unresolved = 0
      wi_bind = .true.
      do j = 1, nc
         q = abs(rr(j))/max(ss(j), cert_scale_floor)
         call mass_row_cell_verdict(q, floors(j), tt, dd, wi, an, st)
         if (st .eq. cert_mass_row_unresolved) then
            ent%n_mass_unresolved = ent%n_mass_unresolved + 1
            if (ent%n_mass_unresolved .eq. 1) ent%j_first_mass_unresolved = j
         endif
         if (ent%jbind .eq. 0 .or. dd .gt. ent%dist_bind) then
            ent%jbind        = j
            ent%dist_bind    = dd
            ent%row_at_bind  = q
            ent%tol          = tt
            ent%tol_anchored = an
            wi_bind          = wi
         endif
         if (j .eq. ent%jworst) ent%tol_at_jworst = tt
      enddo
      ent%within_tol = ent%finite .and. wi_bind
      end subroutine mass_row_column_verdict

      ! ------------------------------------------------------!

      subroutine carrier_row_entry(rep, ic, cres, cterms, cabsent)
      ! The steady balance of one transported carrier, transport minus
      ! reaction, on the row's own terms.
      !
      ! WHY PER CARRIER AND PER CELL. carrier_steady_residual also forms a
      ! ratio of SUMS over all carriers of a region; a species whose rates
      ! are small next to a balanced dominant one can be out by any factor
      ! of its own terms and move that ratio by nothing at all (review R7).
      ! The row measure below is that species' own equation.
      type(cert_report), intent(inout) :: rep
      integer,           intent(in)    :: ic
      real*8, dimension(1:N,n_carrier_max), intent(in) :: cres, cterms
      logical, dimension(1:N,n_carrier_max), intent(in) :: cabsent
      real*8, dimension(1:N) :: rr, ss, ww
      logical, dimension(1:N) :: aa
      integer :: j, idx, nro_last, nro_total, nro_sub
      call add_entry(rep, 'carrier balance '//trim(carrier_name(ic)),      &
                     cert_evaluated, '', idx)
      rep%e(idx)%regime_gated = .true.
      rep%e(idx)%tol = cert_tol_carrier_at(cert_regime_wind_r)
      ! The round-off limited rows of the transport steps this run has
      ! taken, carried as an informational field of the carrier entry.
      call carrier_roundoff_limited_record(nro_last, nro_total, nro_sub)
      rep%e(idx)%n_roundoff_limited = nro_total
      rep%e(idx)%units_floor = 'row_terms [cm^-3 s^-1], the sum of the '// &
           'row''s own terms with the 1e-20 free-element floor inside it'
      ! AND AN ABSOLUTE FLOOR UNDER THE SPECIES ITSELF, on the same footing
      ! as the element rows' 1e-20 rho X_base: a carrier below
      ! carrier_absent_fraction of the free reservoir of its own element is
      ! reported and does not gate. Without it a relative row measure gates
      ! on a species that is not there -- MEASURED on the LHS 1140 b
      ! molecular wind, n(H2) = 5.6e-27 cm^-3 against a gas of 1e6 cm^-3 at
      ! 2.7 R_p, whose row read exactly 1.000 and refused every state.
      rep%e(idx)%n_absent       = 0
      rep%e(idx)%j_first_absent = 0
      do j = 1, N
         rr(j) = cres(j,ic)
         ss(j) = cterms(j,ic)
         aa(j) = cabsent(j,ic)
         ww(j) = r(j)*r(j)*dr_j(j)
         if (aa(j)) then
            rep%e(idx)%n_absent = rep%e(idx)%n_absent + 1
            if (rep%e(idx)%j_first_absent .eq. 0)                         &
               rep%e(idx)%j_first_absent = j
         endif
      enddo
      ! THE WHOLE-COLUMN MEASURE KEEPS EVERY CELL. What the absence changes
      ! is the VERDICT, not the measurement: a reader still sees the row
      ! measure of the empty cells, and the count below says how many of
      ! them there were.
      call certification_row_measure(N, j_min, rr, ss,                    &
                                     rep%e(idx)%row_max,                  &
                                     rep%e(idx)%jworst,                   &
                                     rep%e(idx)%finite,                   &
                                     wvol = ww, rvol = rep%e(idx)%row_vol)
      call fill_regime_measures(rep%e(idx), rr, ss)
      call gate_species_row(rep%e(idx), .true., rr, ss, absent = aa)
      end subroutine carrier_row_entry

      ! ------------------------------------------------------!

      subroutine transport_row_entry(rep, name, active, measured, res,    &
                                     scale, units_floor, tol_whole_column)
      ! One equation of the inventory whose residual and scale are given in
      ! the row's own physical units: the elemental transport balances and
      ! the H(n=2) level balance. `active` says the configuration carries
      ! the unknown; `measured` says a state-consistent residual of it was
      ! obtained. Active and unmeasured is `unavailable`, never a zero.
      type(cert_report), intent(inout) :: rep
      character(len=*),  intent(in)    :: name
      logical,           intent(in)    :: active, measured
      real*8,            intent(in)    :: res(1:N), scale(1:N)
      character(len=*),  intent(in)    :: units_floor
      ! A row whose tolerance is ONE NUMBER OVER THE WHOLE COLUMN: the
      ! H(n=2) level balance is such a row. Without it the entry is an
      ! elemental transport row and is judged by cert_tol_element_at,
      ! which gates the wind cells and reports the rest.
      real*8, optional,  intent(in)    :: tol_whole_column
      real*8, dimension(1:N) :: ww
      integer :: j, idx
      if (.not. active) then
         call add_entry(rep, name, cert_not_applicable,                    &
                        'the run does not carry this unknown', idx)
         return
      endif
      if (.not. measured) then
         call add_entry(rep, name, cert_unavailable,                       &
                        'no state-consistent measurement of this equation',&
                        idx)
         return
      endif
      call add_entry(rep, name, cert_evaluated, '', idx)
      rep%e(idx)%regime_gated = .not. present(tol_whole_column)
      if (present(tol_whole_column)) then
         rep%e(idx)%tol = tol_whole_column
      else
         rep%e(idx)%tol = cert_tol_element_at(cert_regime_wind_r)
      endif
      rep%e(idx)%units_floor = units_floor
      do j = 1, N
         ww(j) = r(j)*r(j)*dr_j(j)
      enddo
      call certification_row_measure(N, j_min, res, scale,                 &
                                     rep%e(idx)%row_max,                   &
                                     rep%e(idx)%jworst,                    &
                                     rep%e(idx)%finite,                    &
                                     wvol = ww, rvol = rep%e(idx)%row_vol)
      call fill_regime_measures(rep%e(idx), res, scale)
      if (rep%e(idx)%regime_gated) then
         call gate_species_row(rep%e(idx), .false., res, scale)
      else
         rep%e(idx)%within_tol = rep%e(idx)%finite .and.                   &
                                 (rep%e(idx)%row_max .lt. rep%e(idx)%tol)
      endif
      end subroutine transport_row_entry

      ! ------------------------------------------------------!

      subroutine unit_scale_entry(rep, name, active, measured, res, tol,   &
                                  why, units_floor)
      ! One equation whose residual ALREADY carries its scale: the
      ! eliminated-species closure and the He 2^3S level row are measured in
      ! units of the turnover rate of their own row, so the row measure is
      ! the maximum of the values themselves and the scale is 1.
      type(cert_report), intent(inout) :: rep
      character(len=*),  intent(in)    :: name
      logical,           intent(in)    :: active, measured
      real*8,            intent(in)    :: res(1:N), tol
      character(len=*),  intent(in)    :: why, units_floor
      real*8, dimension(1:N) :: ww, ss
      integer :: j, idx
      if (.not. active) then
         call add_entry(rep, name, cert_not_applicable,                    &
                        'the run does not carry this unknown', idx)
         return
      endif
      if (.not. measured) then
         call add_entry(rep, name, cert_unavailable, why, idx)
         return
      endif
      call add_entry(rep, name, cert_evaluated, '', idx)
      rep%e(idx)%tol         = tol
      rep%e(idx)%units_floor = units_floor
      do j = 1, N
         ww(j) = r(j)*r(j)*dr_j(j)
         ss(j) = 1.0d0
      enddo
      call certification_row_measure(N, j_min, res, ss,                    &
                                     rep%e(idx)%row_max,                   &
                                     rep%e(idx)%jworst,                    &
                                     rep%e(idx)%finite,                    &
                                     wvol = ww, rvol = rep%e(idx)%row_vol)
      call fill_regime_measures(rep%e(idx), res, ss)
      rep%e(idx)%within_tol = rep%e(idx)%finite .and.                      &
                              (rep%e(idx)%row_max .le. rep%e(idx)%tol)
      end subroutine unit_scale_entry

      ! ------------------------------------------------------!

      subroutine mass_closure_of_state(rho, f_sp, dev, jworst)
      ! max_j |sum_i f_i A_i - 1| over the physical cells, and where it sits.
      !
      ! The mass fractions of a state are f_i = n_i m_i / rho by definition,
      ! so the sum is one and this number is zero for a state whose two
      ! halves describe one gas. It is not zero: the composition returned by
      ! a sweep is normalized by the density the sweep was given and is only
      ! that density's composition to the arithmetic of the sweep, and
      ! without the projection of item L19 that departure ratchets by about
      ! 1e-14 per outer pass of a stationary solve.
      !
      ! calc_rho is the run's own mass policy term for term (the trace-metal
      ! mass under eos_metals, the molecular and oxygen-carrier masses, He
      ! 2^3S inside the He I column), so the ratio it gives is exactly the
      ! sum above and no second mass table exists here to drift from it.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8,                                 intent(out) :: dev
      integer,                                intent(out) :: jworst
      real*8, dimension(1-Ng:N+Ng)         :: rho_of_comp
      real*8, dimension(1-Ng:N+Ng,n_mion)  :: nm_c
      real*8, dimension(1-Ng:N+Ng,4)       :: nmol_c
      real*8, dimension(1-Ng:N+Ng,3)       :: nox_c
      real*8  :: d
      integer :: j, im
      dev = 0.0d0;  jworst = 0
      do im = 1, n_mion
         nm_c(:,im) = rho*f_sp(:,mion_fsp(im))
      enddo
      nmol_c(:,1) = rho*f_sp(:,isp_H2)
      nmol_c(:,2) = rho*f_sp(:,isp_H2p)
      nmol_c(:,3) = rho*f_sp(:,isp_H3p)
      nmol_c(:,4) = rho*f_sp(:,isp_HeHp)
      nox_c(:,1)  = rho*f_sp(:,isp_OH)
      nox_c(:,2)  = rho*f_sp(:,isp_H2O)
      nox_c(:,3)  = rho*f_sp(:,isp_CO)
      call calc_rho(rho*f_sp(:,isp_HI),   rho*f_sp(:,isp_HII),            &
                    rho*f_sp(:,isp_HeI),  rho*f_sp(:,isp_HeII),           &
                    rho*f_sp(:,isp_HeIII), rho_of_comp,                   &
                    nm_c, nmol_c, nox_c)
      do j = 1, N
         if (rho(j) .le. 0.0d0) cycle
         d = abs(rho_of_comp(j) - rho(j))/rho(j)
         if (d .gt. dev) then
            dev = d;  jworst = j
         endif
      enddo
      end subroutine mass_closure_of_state

      ! ------------------------------------------------------!

      integer function base_contact_direction_agreement(s_wind, M_wind,  &
                        flux_base, window_mean, window_available)
      ! Which quantity decided the direction of the base contact at the last
      ! boundary evaluation, and, where the wind window decided it, whether
      ! the base face mass flux the Riemann solve assembled for the same
      ! state carries the flux in the same direction. The comment at the
      ! report that prints it gives the reason (item D2b). A base face flux
      ! at or below 1e-12 of the window mean is at the rounding of the flux
      ! and carries no direction.
      real*8,  intent(in) :: s_wind, M_wind, flux_base, window_mean
      logical, intent(in) :: window_available
      if (.not. (s_wind .gt. 0.0d0 .and. window_available .and.            &
                 window_mean .ne. 0.0d0)) then
         base_contact_direction_agreement =                                &
              contact_direction_from_face_velocity
      else if (abs(flux_base) .le. 1.0d-12*abs(window_mean)) then
         base_contact_direction_agreement =                                &
              contact_direction_face_flux_at_rounding
      else if ((M_wind .gt. 0.0d0) .eqv. (flux_base .gt. 0.0d0)) then
         base_contact_direction_agreement = contact_direction_window_agrees
      else
         base_contact_direction_agreement =                                &
              contact_direction_window_disagrees
      endif
      end function base_contact_direction_agreement

      ! ------------------------------------------------------!

      subroutine ionization_stage_sum_entry(rep, ien, dmax, jworst,       &
                                            known, g)
      ! THE STAGE SUM IDENTITY OF ONE ELEMENT AS AN ENTRY OF THE REPORT.
      !
      !    sum_k F_k(f) + F_close(f) = N_el(f)   at every face,
      !
      ! the sum over ALL stages of that element (ionization_stage_transport,
      ! equation 2), measured by the stationary evaluation that formed the
      ! stage fluxes and read here rather than re-formed, so the number the
      ! report carries is the one the rows were built with.
      !
      ! IT GATES, at the floating-point bound of the identity. The identity
      ! is algebraic and has no truncation error, so the whole of what a
      ! tolerance on it can clear is the rounding of the sums, which
      ! ionization_stage_sum_rounding_bound states and this entry reads from
      ! that one place. The bound is 2 (5 nk + 2) eps g, and g is the ratio
      ! of magnitude sums MEASURED at the face the measure was taken at, so
      ! the row is gated at the rounding of its own face and not at the
      ! largest rounding any face could carry. Where the operator reports no
      ! g (no face of the column carried a nonzero scale) the bound falls
      ! back on the algebraic ceiling of g, 3/2 for one carried stage and 2
      ! for two.
      !
      ! A measure above the bound is a broken construction and not a state
      ! the solver could have done better on; the constructions that break
      ! the identity on purpose stand at 1.5e-4 to 3.6e-4 (MEASURED,
      ! src/tests/ionization_stage_flux/), ten decades above it, and the
      ! states this code produces at 0.07 to 0.08 of it. The derivation, the
      ! executed validation and the production readings are anchor (8) of
      ! docs/certification_tolerance_anchoring_20260910.md.
      !
      ! ONE ENTRY PER ELEMENT. The identity is a statement about one
      ! element's own nucleus flux, so a single entry over both elements
      ! could not say which element's construction it belongs to; that is
      ! the same reason each trace element's transport balance is its own
      ! entry.
      !
      ! The number of carried stages is the element's own and not a free
      ! parameter: hydrogen carries x(H II) and helium x(He II) and
      ! x(He III), and the measure is known only when every one of them is
      ! solved (ionization_stage_nucleus_sum).
      type(cert_report), intent(inout) :: rep
      integer,           intent(in)    :: ien
      real*8,            intent(in)    :: dmax
      integer,           intent(in)    :: jworst
      logical,           intent(in)    :: known
      ! The measured ratio of magnitude sums at that face. Absent, or
      ! negative, means the operator measured none and the bound is taken at
      ! its ceiling.
      real*8, optional,  intent(in)    :: g
      character(len=2)  :: estg_name(2)
      character(len=76) :: scale_text
      integer :: nk, idx
      real*8  :: gm
      estg_name(ien_H)  = 'H '
      estg_name(ien_He) = 'He'
      if (.not. known) then
         call add_entry(rep, 'ionization stage nucleus sum '//            &
                        trim(estg_name(ien)), cert_unavailable,           &
                        'this element carries no transported'//           &
                        ' stage, or no stationary evaluation has'//       &
                        ' formed its stage fluxes yet', idx)
         return
      endif
      nk = 1
      if (ien .eq. ien_He) nk = 2
      gm = -1.0d0
      if (present(g)) gm = g
      call add_entry(rep, 'ionization stage nucleus sum '//               &
                     trim(estg_name(ien)), cert_evaluated,                &
                     'gating at the floating-point bound of an'//         &
                     ' algebraic identity (anchor 8)', idx)
      rep%e(idx)%row_max     = dmax
      rep%e(idx)%jworst      = jworst
      if (gm .gt. 0.0d0) then
         rep%e(idx)%tol = ionization_stage_sum_rounding_bound(nk, gm)
         write(scale_text,'(A,F9.5,A)')                                   &
              'relative to max(|N_el|, sum_k |F_k|) at the face; g =',    &
              gm, ' there'
      else
         rep%e(idx)%tol = ionization_stage_sum_rounding_bound(nk)
         scale_text =                                                     &
           'relative to max(|N_el|, sum_k |F_k|); g at its ceiling'
      endif
      rep%e(idx)%finite      = finite_real(dmax)
      rep%e(idx)%within_tol  = rep%e(idx)%finite .and.                    &
                               (dmax .le. rep%e(idx)%tol)
      rep%e(idx)%units_floor = scale_text
      end subroutine ionization_stage_sum_entry

      ! ------------------------------------------------------!

      subroutine add_entry(rep, name, status, reason, ient)
      ! THE ONE PLACE A REPORT GAINS AN ENTRY, AND THE ONLY PLACE AN ENTRY
      ! INDEX COMES FROM. ient is the index of the entry this call created,
      ! and every field of that entry is written through it. THE INVARIANT:
      ! an entry writer only ever addresses the entry its own call created.
      ! A writer that addressed rep%n instead would address whatever entry
      ! stood last when its own was not created, and the report would carry
      ! one equation's name over another equation's measure, tolerance and
      ! verdict.
      !
      ! AN ENTRY THAT DOES NOT FIT STOPS THE RUN, before any field is
      ! written, so the entries the report already holds stand as they were.
      ! cert_max_entries is the size of the inventory itself (its
      ! declaration says how it is formed), so no configuration reaches
      ! this; a configuration that did would have grown an equation past the
      ! count that declares the inventory, and a certification written from
      ! a report that cannot hold all of its equations is not a statement
      ! about the equations the run solves. The refusal names the equation
      ! that did not fit, the capacity, and the entry that stands last and
      ! unchanged.
      type(cert_report), intent(inout) :: rep
      character(len=*),  intent(in)    :: name, reason
      integer,           intent(in)    :: status
      integer,           intent(out)   :: ient
      if (rep%n .ge. cert_max_entries) then
         write(*,'(A,A,A)') ' (certification) REPORT FULL: the equation "',&
              trim(name), '" does not fit'
         write(*,'(A,I0,A)') ' (certification)   capacity '//             &
              'cert_max_entries = ', cert_max_entries,                    &
              ' entries, all of them in use'
         write(*,'(A,A,A,ES13.6)') ' (certification)   last entry '//     &
              'unchanged: ', trim(rep%e(cert_max_entries)%name),          &
              '  max=', rep%e(cert_max_entries)%row_max
         flush(6)
         error stop 1
      endif
      rep%n = rep%n + 1
      ient  = rep%n
      rep%e(ient)             = cert_entry()
      rep%e(ient)%name        = name
      rep%e(ient)%status      = status
      rep%e(ient)%reason      = reason
      end subroutine add_entry

      ! ------------------------------------------------------!

      integer function active_system_variant() result(iv)
      ! Which of the seven coupled cell systems (B1a section 2.4) the
      ! configuration selects. One of them is always active: every run
      ! eliminates the neutral stage of each element and the electron
      ! density by a closure.
      if (thereis_mol) then
         iv = 6;  if (thereis_metals) iv = 7
      else if (thereis_He) then
         if (thereis_HeITR) then
            iv = 3;  if (thereis_metals) iv = 5
         else
            iv = 2;  if (thereis_metals) iv = 4
         endif
      else
         iv = 1
      endif
      end function active_system_variant

      ! ------------------------------------------------------!

      function system_variant_name(iv) result(nm)
      integer, intent(in) :: iv
      character(len=20)   :: nm
      select case (iv)
         case (1);     nm = 'System_H'
         case (2);     nm = 'System_HeH'
         case (3);     nm = 'System_HeH_TR'
         case (4);     nm = 'System_HeH_metals'
         case (5);     nm = 'System_HeH_TR_metals'
         case (6);     nm = 'System_HeH_mol'
         case default; nm = 'System_HeH_mol_metals'
      end select
      end function system_variant_name

      ! ------------------------------------------------------!

      subroutine read_validity_states(rep)
      ! The five validity states of B1a section 4, read from the producers
      ! that exist. A state with no producer is left negative and reported
      ! as "not produced": the report says what the run measured and does
      ! not invent the rest.
      type(cert_report), intent(inout) :: rep
      integer  :: nco_out, nco_hot, nco_hep
      real*8   :: co_ratio, co_r, co_form
      ! 4.1 active unvalidated physics: NOT PRODUCED. The one item of that
      ! list with a cumulative record was the CO thermal ceiling, a
      ! thermodynamic bound standing in for a kinetic process; it is gone,
      ! the CO row carries published destruction rates, and the CO items are
      ! counted under 4.2 as a closure with a domain.
      rep%n_active_unvalidated_physics = -1
      ! 4.2 OUT-OF-DOMAIN CLOSURE. The H3+ cooling model is the first
      ! closure in the code that records its own domain (B3b-H3): the cells
      ! at which the collider density fell below the tabulated range, at
      ! which it rose above it, and at which the temperature was clamped to
      ! the low or the high end of the published fits.
      !
      ! THEY ARE INFORMATIONAL AND DO NOT INVALIDATE BY THEMSELVES
      ! (decision 8, 2026-09-06): the extrapolating branch below the
      ! tabulated collider range is the exact collisional limit of the same
      ! model, not a guess outside it, so a cell counted there is described
      ! by the model and not merely clamped to its edge. The count is
      ! reported so that a reader can see how much of the domain was
      ! outside the table; a closure whose out-of-domain branch is NOT a
      ! limit of its own model has to say so where it is counted.
      !
      ! THE CO DESTRUCTION MODEL IS THE SECOND SUCH CLOSURE (decision 10,
      ! 2026-09-06). Its row destroys CO and never forms it, which holds
      ! where tau_dest << tau_res << tau_form; the cells where that ordering
      ! fails are counted here. They are informational for the same reason:
      ! where the destruction is too slow to matter the omitted formation is
      ! too, so the transported value stands, which is what a transport
      ! operator should do in a quenched layer.
      rep%n_out_of_domain_closure = 0
      if (thereis_mol) rep%n_out_of_domain_closure =                      &
           h3p_n_below_collider + h3p_n_below_fit_T                        &
           + h3p_n_above_fit_T + h3p_n_outside_nonlte_T
      rep%n_h3p_below_collider   = h3p_n_below_collider
      rep%n_h3p_below_fit_T      = h3p_n_below_fit_T
      rep%n_h3p_above_fit_T      = h3p_n_above_fit_T
      rep%n_h3p_outside_nonlte_T = h3p_n_outside_nonlte_T
      ! The Ly-alpha record is zero unless jlya_escape_prob ran (the
      ! accessor reports 0 seen); it joins the out-of-domain count like
      ! the H3+ and CO records, and like them it does not invalidate.
      call lya_wing_domain_record(rep%n_lya_wing_out, rep%n_lya_wing_seen, &
                                  rep%lya_wing_worst, rep%lya_wing_limit)
      rep%n_out_of_domain_closure = rep%n_out_of_domain_closure           &
                                    + rep%n_lya_wing_out
      if (thereis_oxychem .and. carrier_transport) then
         call carrier_co_domain_record(nco_out, nco_hot, nco_hep,         &
                                       co_ratio, co_r, co_form)
         rep%n_co_out_of_domain      = nco_out
         rep%n_co_above_shield_T     = nco_hot
         rep%n_co_HeII_led           = nco_hep
         rep%co_worst_dest_over_res  = co_ratio
         rep%n_out_of_domain_closure = rep%n_out_of_domain_closure        &
                                       + nco_out
      endif
      ! 4.3 rejected numerical trial with no adopted contribution: counted,
      ! but by locals of the marching loop, not by a module a certification
      ! can read. It does not invalidate in any case.
      rep%n_rejected_trial             = -1
      ! 4.4 unbudgeted accepted correction: a change of a conserved
      ! quantity in the ADOPTED state with no source term behind it.
      ! NEITHER temperature floor is one any more. The energy update and
      ! the conduction stage are both bracketed and residual-controlled,
      ! and reaching the lower bracket end is a FAILURE that stops the run
      ! rather than a state the run adopts; a rejected attempt contributes
      ! nothing to the state being judged (review 2 section 5.3). Both
      ! counters are therefore reported as attempts beside the validity
      ! states and neither invalidates. The category itself stays: a
      ! producer of a genuine unbudgeted correction reaches it here.
      ! Row 12 of the attempted step is one of them: the Shapiro filter
      ! alters the adopted state with no source term behind the change.
      ! Advisor decision 8 (design section 9): it is carried until B3b
      ! reaches it, is recorded when it fires, and NO STATE CERTIFIES WHILE
      ! IT IS ACTIVE. That is the interim contract made visible where it
      ! decides something, rather than left in a comment. Row 8's
      ! composition projection was the other one; B3c removed its producer,
      ! so it can no longer contribute here.
      rep%n_unbudgeted_accepted_correction = n_shapiro_applied
      rep%n_energy_floor_attempts          = n_energy_floor_hits
      rep%n_conduction_floor_attempts      = n_conduction_floor_hits
      ! 4.5 specified external reservoir: informational, and reported by the
      ! resolved-configuration record rather than by a counter.
      rep%n_specified_external_reservoir = -1
      end subroutine read_validity_states

      ! ------------------------------------------------------!

      subroutine certification_report_write(rep, label, face_budget)
      ! One block per certification, with every entry of the inventory, its
      ! status, its measure and the scale that measure is taken in.
      !
      ! face_budget, when the caller measured one, is the mass flux the same
      ! operator puts through every face of the same state (stationary_
      ! operator). It is REPORTED ONLY: no entry of the inventory reads it
      ! and no verdict of this report depends on it.
      type(cert_report),           intent(in) :: rep
      character(len=*),            intent(in) :: label
      type(face_mass_flux_budget), intent(in), optional :: face_budget
      integer :: i
      character(len=14) :: st
      real*8 :: fb_scale
      write(*,'(A)') ' '
      write(*,'(A,A)') ' (certification) ', trim(label)
      write(*,'(A,I0,A,I0,A)') '   active equations ',                     &
           certification_active_equation_count(rep), ' of ', rep%n,        &
           ' in the inventory'
      ! THE TOLERANCES THIS REPORT JUDGED BY, printed once so that a reader
      ! does not have to find them in the source, and so that a refusal can
      ! be read against the number that refused it. Each row's own value is
      ! printed beside it below as tol=.
      write(*,'(A)') '   tolerances (convergence study, contract'//        &
           ' sections 9 and 10; the hydrodynamic values were'
      write(*,'(A)') '   re-anchored on the self-consistent JFNK'//        &
           ' residual after B5a):'
      write(*,'(A,ES9.2,A,ES9.2,A,ES9.2)')                                 &
           '     hydrodynamic mass ', cert_tol_mass,                       &
           ', momentum ', cert_tol_momentum,                               &
           ', energy ', cert_tol_energy
      write(*,'(A,F5.1,A)') '     the mass value is a FLOOR: where the'// &
           ' cell''s own rounding of the flux difference stands above'//   &
           ' it,',                                                        &
           cert_mass_round_margin, ' times that floor is the tolerance'
      write(*,'(A)') '     of that cell (cert_tol_mass_at), and the'//    &
           ' verdict on the mass row is taken cell by cell'
      write(*,'(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)')                         &
           '     closure ', cert_tol_closure,                              &
           ', level ', cert_tol_level,                                     &
           '; and IN THE WIND, at r >=', cert_regime_wind_r,               &
           ', carrier ', cert_tol_carrier_at(cert_regime_wind_r)
      write(*,'(A,ES9.2,A)')                                               &
           '     and elemental transport ',                                &
           cert_tol_element_at(cert_regime_wind_r),                        &
           '; the species rows of the cells below that radius are'
      write(*,'(A,ES9.2,A)') '     REPORTED AND DO NOT GATE: the'//        &
           ' layer r <', cert_regime_layer_r, ', whose element fluxes'
      write(*,'(A)') '     are not conserved (spread 9.8), and the'//      &
           ' band above it, whose discretization (4.6e-4 at the'
      write(*,'(A)') '     binding radius) stands above any tolerance'//   &
           ' one could set there'
      if (cert_resid_tol_of_run .gt. 0.0d0)                                &
         write(*,'(A,ES9.2,A)') '     the run''s own "Resid tol" is ',     &
              cert_resid_tol_of_run, ', which is the SOLVER''s target'//   &
              ' and not a certification tolerance'
      do i = 1, rep%n
         select case (rep%e(i)%status)
            case (cert_not_applicable); st = 'not_applicable'
            case (cert_evaluated);      st = 'evaluated'
            case default;               st = 'UNAVAILABLE'
         end select
         if (rep%e(i)%status .eq. cert_evaluated) then
            write(*,'(A,A52,A,A14,A,ES10.3,A,ES10.3,A,I0,A,ES8.1,A,A)')    &
                 '   ', rep%e(i)%name, ' ', st,                            &
                 '  max=', rep%e(i)%row_max,                               &
                 '  vol=', rep%e(i)%row_vol,                               &
                 '  cell=', rep%e(i)%jworst,                               &
                 '  tol=', rep%e(i)%tol,                                   &
                 '  ', merge('within  ', 'ABOVE   ',                       &
                             rep%e(i)%within_tol)
            if (.not. rep%e(i)%finite)                                     &
               write(*,'(A)') '        the row is not a finite number'
            write(*,'(A,A)') '        scale: ', trim(rep%e(i)%units_floor)
            ! WHICH PART OF A REGIME-GATED ROW DECIDED, so that max= above
            ! is never read as the number the verdict was taken on.
            if (rep%e(i)%regime_gated .and. rep%e(i)%jworst_gate .gt. 0)   &
               write(*,'(A,ES10.3,A,ES9.2,A,ES10.3,A,I0,A)')               &
                    '        gated at r >=', cert_regime_wind_r,          &
                    ' against ', rep%e(i)%tol, ': ',                       &
                    rep%e(i)%row_max_gate, ' at cell ',                    &
                    rep%e(i)%jworst_gate,                                  &
                    '; the cells below that radius are reported only'
            ! WHICH CELL AND WHICH TOLERANCE TOOK THE VERDICT on a row
            ! whose tolerance is a function of the cell, so that a reader
            ! can tell the fixed gate from the rounding anchor.
            if (rep%e(i)%jbind .gt. 0)                                     &
               write(*,'(A,I0,A,ES12.5,A,ES10.3,A,ES8.1,A,ES10.3,A,A)')    &
                    '        verdict at cell ', rep%e(i)%jbind,            &
                    ' (r=', r(rep%e(i)%jbind), '): ',                      &
                    rep%e(i)%row_at_bind, ' against ', rep%e(i)%tol,       &
                    ', distance ', rep%e(i)%dist_bind, ' -- ',             &
                    trim(cert_mass_gate_name(rep%e(i)%tol_anchored))
            if (rep%e(i)%jbind .gt. 0 .and.                                &
                rep%e(i)%jbind .ne. rep%e(i)%jworst)                       &
               write(*,'(A,I0,A,ES8.1,A)')                                 &
                    '        the largest measure sits at cell ',           &
                    rep%e(i)%jworst, ', whose own tolerance is ',          &
                    rep%e(i)%tol_at_jworst, ' there'
            ! REVIEW R5, ITEM Q3: an unresolved mass row is named even when
            ! a different, worse-refusing cell is the one that binds the
            ! verdict, so a reader is never left inferring "certified" from
            ! a within_tol=F line whose text describes a different cell.
            if (rep%e(i)%n_mass_unresolved .gt. 0)                         &
               write(*,'(A,I0,A,I0,A)')                                    &
                    '        mass row unresolved: the estimated'//         &
                    ' rounding floor of the flux difference reaches'//     &
                    ' the row itself at cell ',                            &
                    rep%e(i)%j_first_mass_unresolved,                      &
                    ' (', rep%e(i)%n_mass_unresolved,                      &
                    ' cell(s) of the column)'
            if (rep%e(i)%n_absent .gt. 0)                                  &
               write(*,'(A,I0,A,I0,A)')                                    &
                    '        the species is absent (below 1e-20 of the'//  &
                    ' free reservoir of its element) in ',                 &
                    rep%e(i)%n_absent, ' cell(s), first cell ',            &
                    rep%e(i)%j_first_absent,                               &
                    '; those cells are reported and do not gate'
            if (rep%e(i)%n_roundoff_limited .gt. 0)                        &
               write(*,'(A,I0,A)') '        round-off limited rows over'// &
                    ' the run: ', rep%e(i)%n_roundoff_limited,             &
                    ' (accepted and flagged; informational, never'//       &
                    ' decisive)'
         else if (rep%e(i)%status .eq. cert_unavailable) then
            write(*,'(A,A52,A,A14,A,A)') '   ', rep%e(i)%name, ' ', st,    &
                 '  ', trim(rep%e(i)%reason)
         endif
      enddo
      ! THE MASS FLUX THROUGH EVERY FACE, on the state the rows above were
      ! measured on and under the operator they were measured with. In
      ! spherical symmetry the mass equation transports r_f^2 (rho v)_f
      ! across face f, so on a stationary state that number is the wind's
      ! own mass flux at every face, and the base face carries it too.
      ! REPORTED ONLY: no entry of the inventory reads it, no tolerance
      ! judges it, and the verdict below does not depend on it.
      if (present(face_budget)) then
      if (face_budget%available) then
         write(*,'(A)') '   face mass flux r_f^2 (rho v)_f of this'//      &
              ' state (reported beside the rows above; it gates nothing):'
         write(*,'(A,A,A,A,A,L1)') '     operator: ',                      &
              trim(face_budget%operator_name), ', numerical flux: ',       &
              trim(face_budget%flux_name), ', well balanced: ',            &
              face_budget%well_balanced
         if (face_budget%window_available .and.                            &
             face_budget%window_mean .ne. 0.0d0) then
            fb_scale = face_budget%window_mean
            write(*,'(A,ES12.5)') '     in units of the wind-window mean'//&
                 ' of rho v r^2, which is ', fb_scale
            write(*,'(A,ES22.15,A,I0)') '       minimum over the faces ',  &
                 face_budget%flux_min/fb_scale, ' at face ',               &
                 face_budget%j_min
            write(*,'(A,ES22.15,A,I0)') '       maximum over the faces ',  &
                 face_budget%flux_max/fb_scale, ' at face ',               &
                 face_budget%j_max
            write(*,'(A,ES22.15,A,ES12.5)') '       base face ',           &
                 face_budget%flux_base/fb_scale,                           &
                 ', offset from the window mean ',                         &
                 face_budget%flux_base/fb_scale - 1.0d0
            write(*,'(A)') '       a face value and the mean of a'//       &
                 ' cell-centered product are not the same quantity, so'//  &
                 ' an offset of the size of the'
            write(*,'(A)') '       wind''s own discretization is not a'//  &
                 ' mass leak'
         else
            write(*,'(A)') '     the wind window carries no usable mean'// &
                 ' flux on this state, so the faces are reported in code'//&
                 ' units'
            write(*,'(A,ES22.15,A,I0,A,ES22.15,A,I0)')                     &
                 '       minimum ', face_budget%flux_min, ' at face ',     &
                 face_budget%j_min, ', maximum ', face_budget%flux_max,    &
                 ' at face ', face_budget%j_max
            write(*,'(A,ES22.15)') '       base face ',                    &
                 face_budget%flux_base
         endif
         ! THE DIRECTION OF THE BASE CONTACT AGAINST THE FLUX IT CARRIES.
         ! Where the wind window has standing, the base boundary upwinds
         ! the contact on the sign of the window's mass flux and not on the
         ! matched face velocity, whose sign on a converged state is set by
         ! the boundary's own pressure residual (item D2b, user decision of
         ! 2026-09-19). That choice presumes the window and the base face
         ! carry the flux in one direction; the base face mass flux the
         ! Riemann solve assembled for this same state is the local
         ! measure, so a disagreement is stated here instead of assumed
         ! absent. A base face flux below 1e-12 of the window mean is at
         ! the rounding of the flux and carries no direction.
         select case (base_contact_direction_agreement(                    &
                          base_face_swind_last, base_face_Mwind_last,      &
                          face_budget%flux_base, face_budget%window_mean,  &
                          face_budget%window_available))
         case (contact_direction_from_face_velocity)
            write(*,'(A,F4.1,A)') '     base contact direction: read'//   &
                 ' from the matched face velocity (the window has no'//    &
                 ' standing; w_rev =', base_face_blend_last, ')'
         case (contact_direction_face_flux_at_rounding)
            write(*,'(A)') '     base contact direction: read from'//     &
                 ' the wind window; the base face flux is at its'//        &
                 ' rounding and carries no direction'
         case (contact_direction_window_agrees)
            write(*,'(A,F4.1,A)') '     base contact direction:'//        &
                 ' read from the wind window, which agrees with the'//     &
                 ' base face flux (w_rev =', base_face_blend_last, ')'
         case default
            write(*,'(A,F4.1,A)') '     base contact direction:'//        &
                 ' DISAGREEMENT, the wind window and the base face'//      &
                 ' flux carry opposite signs (w_rev =',                    &
                 base_face_blend_last, ')'
            write(*,'(A)') '       the contact was upwinded on the'//     &
                 ' window; the local flux reverses, which is the case'//   &
                 ' the boundary''s direction rule does not cover'
         end select
      endif
      endif
      ! THE ANCHORING BLOCK (EXHALE_CERT_ANCHOR=1), off by default.
      ! The same measures at full double precision and split by regime, so
      ! that two evaluations of one state, two grids, or a perturbed state
      ! can be compared at the precision those comparisons need. Nothing
      ! here is a verdict: the block above holds every number the
      ! certification decides on.
      if (anchor_report_on()) then
         write(*,'(A,ES10.3,A,ES10.3,A)') '   anchoring block: regimes '// &
              'are the layer r <', cert_regime_layer_r,                    &
              ' and the wind r >=', cert_regime_wind_r, '; measures at'//  &
              ' full precision'
         if (anchor_probe_radius() .gt. 0.0d0)                             &
            write(*,'(A,ES16.9)') '     probe radius (EXHALE_CERT_'//      &
                 'ANCHOR_R): ', anchor_probe_radius()
         do i = 1, rep%n
            if (rep%e(i)%status .ne. cert_evaluated) cycle
            write(*,'(A,A52,A,ES23.16,A,I0)') '   anchor ',                &
                 rep%e(i)%name, ' all  ', rep%e(i)%row_max,                &
                 '  cell=', rep%e(i)%jworst
            if (rep%e(i)%jworst_layer .gt. 0) then
               write(*,'(A,ES23.16,A,I0,A,ES12.5)')                        &
                    '          layer ', rep%e(i)%row_max_layer,            &
                    '  cell=', rep%e(i)%jworst_layer,                      &
                    '  r=', r(rep%e(i)%jworst_layer)
            else
               write(*,'(A)') '          layer  no cell in the window'
            endif
            if (rep%e(i)%jworst_wind .gt. 0) then
               write(*,'(A,ES23.16,A,I0,A,ES12.5)')                        &
                    '          wind  ', rep%e(i)%row_max_wind,             &
                    '  cell=', rep%e(i)%jworst_wind,                       &
                    '  r=', r(rep%e(i)%jworst_wind)
            else
               write(*,'(A)') '          wind   no cell in the window'
            endif
            if (rep%e(i)%jprobe .gt. 0)                                    &
               write(*,'(A,ES23.16,A,I0,A,ES12.5)')                        &
                    '          probe ', rep%e(i)%row_at_probe,             &
                    '  cell=', rep%e(i)%jprobe,                            &
                    '  r=', r(rep%e(i)%jprobe)
         enddo
      endif
      if (rep%chem_root_known) then
         write(*,'(A,I0,A)') '   cells without a chemical root: ',         &
              rep%n_no_chem_root, ' (acceptance class 4 or 6)'
         if (rep%first_no_chem_cell .gt. 0)                                &
            write(*,'(A,I0)') '     first such cell: j = ',                &
                 rep%first_no_chem_cell
      else
         write(*,'(A)') '   cells without a chemical root: NOT STATED'//   &
              ' by this caller'
      endif
      ! The mass closure of the composition this state carries: reported
      ! beside the equations and gating nothing (item L19). A state whose
      ! two halves describe one gas reads a few units in the last place.
      if (rep%j_mass_closure .gt. 0) then
         write(*,'(A,ES10.3,A,I0,A,F10.5)') '   mass closure of the'//    &
              ' composition, max |sum_i f_i A_i - 1| = ',                 &
              rep%mass_closure, ' at cell ', rep%j_mass_closure,          &
              ', r =', r(rep%j_mass_closure)
         write(*,'(A)') '     reported only: the sum is one by the'//     &
              ' definition of the mass fractions, so this is the'//       &
              ' arithmetic of the state and no equation of it'
      endif
      write(*,'(A)') '   validity states (B1a section 4):'
      call write_validity('active unvalidated physics',                    &
           rep%n_active_unvalidated_physics)
      call write_validity('out-of-domain closure activations',             &
           rep%n_out_of_domain_closure)
      if (rep%n_out_of_domain_closure .gt. 0) then
         write(*,'(A)') '       H3+ cooling domain records'//              &
              ' (informational; they do not invalidate):'
         write(*,'(A,I0,A,I0,A,I0,A,I0)')                                  &
              '         collider below the table: ',                       &
              rep%n_h3p_below_collider,                                    &
              ', temperature below the fits: ', rep%n_h3p_below_fit_T,     &
              ', above them: ', rep%n_h3p_above_fit_T,                     &
              ', outside the non-LTE range: ', rep%n_h3p_outside_nonlte_T
         if (rep%n_co_out_of_domain .gt. 0) then
            write(*,'(A)') '       one-sided CO destruction domain'//     &
                 ' (informational; it does not invalidate):'
            write(*,'(A,I0,A,ES11.4,A,I0,A,I0)')                          &
                 '         cell visits with tau_dest > f_dom tau_res: ',  &
                 rep%n_co_out_of_domain, ', worst tau_dest/tau_res ',     &
                 rep%co_worst_dest_over_res,                              &
                 ', above the shielding table 512 K limit: ',             &
                 rep%n_co_above_shield_T,                                 &
                 ', with He+ + CO the leading He+ loss: ',                &
                 rep%n_co_HeII_led
         endif
         if (rep%n_lya_wing_out .gt. 0) then
            write(*,'(A)') '       Ly-alpha damping-wing closure domain'//&
                 ' (informational; it does not invalidate):'
            write(*,'(A,I0,A,I0,A,ES11.4,A,ES11.4)')                      &
                 '         cell visits below the (a tau)^(1/3) limit: ',  &
                 rep%n_lya_wing_out, ' of ', rep%n_lya_wing_seen,         &
                 ', smallest (a tau)^(1/3) ', rep%lya_wing_worst,         &
                 ', limit ', rep%lya_wing_limit
         endif
      endif
      call write_validity('rejected trials with no adopted contribution',  &
           rep%n_rejected_trial)
      call write_validity('unbudgeted accepted corrections',               &
           rep%n_unbudgeted_accepted_correction)
      if (n_shapiro_applied .gt. 0)                                       &
         write(*,'(A,I0,A)') '       of which the Shapiro filter of'//    &
              ' step row 12: ', n_shapiro_applied, ' (B3b removes it)'
      call write_validity('specified external reservoirs',                 &
           rep%n_specified_external_reservoir)
      ! Rejected attempts, reported and never decisive: a stage that
      ! reached its lower bracket end failed and stopped the run, so it
      ! left no contribution to any state a certification judges.
      if (.not. rep%carrier_history_ok)                                    &
         write(*,'(A)') '   CARRIER HISTORY NOT CERTIFIABLE: an interval'//&
              ' was left uncovered and the march went on from the'//       &
              ' entry carriers'
      write(*,'(A,I0,A,I0,A)') '   attempts (not validity states): ',      &
           rep%n_energy_floor_attempts, ' energy-update and ',             &
           rep%n_conduction_floor_attempts, ' conduction cell(s) '//       &
           'reached the lower bracket end'
      if (rep%certified) then
         write(*,'(A)') '   CERTIFIED: every active equation was '//       &
              'evaluated and is within its tolerance'
         call write_wind_certification_line(rep)
      else
         write(*,'(A,I0,A)') '   NOT CERTIFIED: ', rep%n_failing,          &
              ' entry/entries of the inventory refuse it'
         do i = 1, rep%n
            if (rep%e(i)%status .eq. cert_unavailable) then
               write(*,'(A,A,A,A)') '     ', trim(rep%e(i)%name),          &
                    ': ', trim(rep%e(i)%reason)
            else if (rep%e(i)%status .eq. cert_evaluated .and.             &
                     .not. rep%e(i)%within_tol) then
               if (rep%e(i)%regime_gated) then
                  write(*,'(A,A,A,ES10.3,A,ES8.1,A,I0,A)') '     ',        &
                       trim(rep%e(i)%name), ': gated row measure ',        &
                       rep%e(i)%row_max_gate, ' above ', rep%e(i)%tol,     &
                       ' at cell ', rep%e(i)%jworst_gate,                  &
                       ' (a wind cell)'
               else if (rep%e(i)%jbind .gt. 0) then
                  ! A row whose tolerance is a function of the cell names
                  ! the cell that REFUSED it, not the cell of its largest
                  ! measure: on this row the two differ, and the largest
                  ! measure may stand inside its own tolerance.
                  write(*,'(A,A,A,ES10.3,A,ES8.1,A,I0)') '     ',          &
                       trim(rep%e(i)%name), ': row measure ',              &
                       rep%e(i)%row_at_bind, ' above ', rep%e(i)%tol,      &
                       ' at cell ', rep%e(i)%jbind
               else
                  write(*,'(A,A,A,ES10.3,A,ES8.1,A,I0)') '     ',          &
                       trim(rep%e(i)%name), ': row measure ',              &
                       rep%e(i)%row_max, ' above ', rep%e(i)%tol,          &
                       ' at cell ', rep%e(i)%jworst
               endif
            endif
         enddo
         if (rep%chem_root_known .and. rep%n_no_chem_root .gt. 0)          &
            write(*,'(A,I0,A,I0)') '     chemistry: ',                     &
                 rep%n_no_chem_root, ' cell(s) without a root, first at '//&
                 'j = ', rep%first_no_chem_cell
         if (.not. rep%chem_root_known)                                    &
            write(*,'(A)') '     chemistry: the caller made no statement'//&
                 ' about it'
         if (rep%n_unbudgeted_accepted_correction .gt. 0)                  &
            write(*,'(A,I0,A)') '     validity: ',                         &
                 rep%n_unbudgeted_accepted_correction,                     &
                 ' unbudgeted accepted correction(s) in the history'
      endif
      write(*,'(A)') ' '
      end subroutine certification_report_write

      ! ------------------------------------------------------!

      function wind_certification_token(rep) result(tok)
      ! The one word the state file carries beside certified=T when the
      ! species rows of the state were gated in the wind and reported
      ! below it. Empty when the state carries no such row, so a file of
      ! the three-unknown route is unchanged.
      type(cert_report), intent(in) :: rep
      character(len=32) :: tok
      integer :: i
      tok = ''
      do i = 1, rep%n
         if (rep%e(i)%status .ne. cert_evaluated) cycle
         if (rep%e(i)%regime_gated) then
            tok = 'certified_in_wind';  return
         endif
      enddo
      end function wind_certification_token

      ! ------------------------------------------------------!

      subroutine write_wind_certification_line(rep)
      ! WHAT A CERTIFIED STATE WITH SPECIES ROWS IS CERTIFIED IN, in one
      ! line: the species rows were judged in the wind alone, so the state
      ! is certified there and its rows below the wind radius are a
      ! measurement and not a claim. The worst of those reported measures
      ! and its cell stand in the line, so that reading it is enough to
      ! know how far the ungated part of the column is.
      type(cert_report), intent(in) :: rep
      real*8  :: worst
      integer :: i, jw, ngated
      character(len=52) :: nm
      ngated = 0;  worst = -1.0d0;  jw = 0;  nm = ''
      do i = 1, rep%n
         if (rep%e(i)%status .ne. cert_evaluated) cycle
         if (.not. rep%e(i)%regime_gated) cycle
         ngated = ngated + 1
         ! The reported part of the row is the whole column with the gated
         ! cells taken out, layer and band together: the cell that binds
         ! the element row of the candidate states sits in the band, so a
         ! line that quoted the layer window alone would leave the largest
         ! ungated measure of the state out of the report.
         if (rep%e(i)%row_max_reported .gt. worst) then
            worst = rep%e(i)%row_max_reported
            jw    = rep%e(i)%jworst_reported
            nm    = rep%e(i)%name
         endif
      enddo
      if (ngated .eq. 0) return
      if (jw .gt. 0) then
         write(*,'(A,ES9.2,A,ES10.3,A,A,A,I0,A,ES10.3)')                  &
              '   certified IN THE WIND, r >=', cert_regime_wind_r,       &
              '; the element and carrier rows below it are reported'//    &
              ' and do not gate: worst ', worst, ' (', trim(nm),          &
              ') at cell ', jw, ', r =', r(jw)
      else
         write(*,'(A,ES9.2,A)')                                           &
              '   certified IN THE WIND, r >=', cert_regime_wind_r,       &
              '; the element and carrier rows below it are reported'//    &
              ' and do not gate: no cell below that radius'
      endif
      end subroutine write_wind_certification_line

      ! ------------------------------------------------------!

      subroutine write_validity(name, n)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: n
      if (n .lt. 0) then
         write(*,'(A,A,A)') '     ', name, ': not produced'
      else
         write(*,'(A,A,A,I0)') '     ', name, ': ', n
      endif
      end subroutine write_validity

      ! ------------------------------------------------------!

      subroutine certification_evaluate_physical_step(u_old, u, u_hydro,  &
                             W, T, f_sp, f_sp_old, dt_g, L_hydro,          &
                             n_no_chem_root, carrier_completed, verd)
      ! THE PHYSICAL-STEP CONTEXT. Evaluated on the TRIAL state at the
      ! adoption boundary of the marching step, in the order
      ! docs/b3a_attempted_step_controller_design_20260906.md section 3.3
      ! fixes: cheapest and most decisive first, so that a state that is
      ! not a state at all is refused before any residual is assembled.
      !
      !   1 finiteness and positivity of u, T and f_sp;
      !   2 the element and charge invariants to round-off;
      !   3 every active returned-state verdict that exists today: the
      !     carriers (A1, through the interval status A3b reports), the
      !     energy update (B2), the conduction stage (B2b), and the
      !     chemical-root count,
      !     a class-4 or class-6 cell being a non-root and therefore not a
      !     root (contract section 3);
      !   4 the time-discrete hydrodynamic residual, the measure B3a
      !     introduces, against cert_tol_hydro_step.
      !
      ! Entries 5 and 6 of that list -- the energy identity T1.5 and the
      ! integration-error estimate -- are the controller's, for the reason
      ! stated at the type above.
      !
      ! NO SIDE EFFECTS: it reads the state and the returned-state records
      ! and writes nothing.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u_old, u, W
      ! The state the Runge-Kutta stages returned. The time-discrete
      ! hydrodynamic row is a verdict on THAT update, so it is measured
      ! there and not at the adoption boundary: the operator split between
      ! the two changes the energy and momentum rows by amounts that belong
      ! to the source operators, and comparing them against the
      ! hydrodynamic operator alone would be a statement about the split.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: u_hydro
      real*8, dimension(3,1-Ng:N+Ng),         intent(in) :: L_hydro
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: T
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp, f_sp_old
      real*8,  intent(in) :: dt_g
      integer, intent(in) :: n_no_chem_root
      logical, intent(in) :: carrier_completed
      type(cert_step_verdict), intent(out) :: verd

      integer :: j, k, i, jbad, ibad
      real*8  :: res, scal, m, dmax, etol, etol2

      verd = cert_step_verdict()

      ! ---- 1. finiteness and positivity -------------------------------
      do j = 1, N
         do k = 1, 3
            if (.not. finite_real(u(k,j))) then
               call step_refuse(verd, cert_step_nonfinite, 2, j, k,        &
                    0.0d0, 0.0d0, 'a conserved variable is not finite')
               return
            endif
         enddo
         if (W(1,j) .le. 0.0d0) then
            call step_refuse(verd, cert_step_positivity, 2, j, 1,          &
                 W(1,j), 0.0d0, 'the mass density is not positive')
            return
         endif
         if (W(3,j) .le. 0.0d0) then
            call step_refuse(verd, cert_step_positivity, 2, j, 3,          &
                 W(3,j), 0.0d0, 'the pressure is not positive')
            return
         endif
         if (.not. finite_real(T(j)) .or. T(j) .le. 0.0d0) then
            call step_refuse(verd, cert_step_positivity, 9, j, 3,          &
                 T(j), 0.0d0,                                              &
                 'the temperature is not positive and finite')
            return
         endif
         do i = 1, n_species
            if (.not. finite_real(f_sp(j,i))) then
               call step_refuse(verd, cert_step_nonfinite, 7, j, i,        &
                    0.0d0, 0.0d0, 'a species fraction is not finite')
               return
            endif
         enddo
      enddo

      ! ---- 2. the element and charge invariants -----------------------
      ! WHAT IS INVARIANT ACROSS ONE STEP, and what is not. The fraction
      ! layout is per unit mass, so the nucleus count of an element in a
      ! cell is the weighted sum of its species fractions, and that sum is
      ! NOT one: it is the element's abundance, which the configuration
      ! fixes. What a step may not change is that sum ITSELF, cell by cell,
      ! because no operation of rows 1 to 12 creates or destroys a nucleus.
      ! The hydrodynamic update does not touch f_sp at all, the ionization
      ! sweep moves nuclei between stages of one element and conserves each
      ! element's total by construction, and the composition projection
      ! rewrites the pressure and not the composition.
      !
      ! The two operations that DO change a cell's element totals are the
      ! element diffusion of row 4 and the carrier transport of row 5:
      ! both move nuclei between cells. Where either is active the cell-by-cell
      ! statement is not an invariant of the step and is not made; the
      ! global statement that would replace it needs the transport across
      ! the two faces of the domain, which no operator reports today, so it
      ! is not invented here. The entry then says which operator is why.
      !
      ! Charge follows: n_e is formed from the stage charges of the same
      ! composition (calc_ne), so a composition inside its element budget
      ! carries the electron count that budget implies and there is no
      ! second statement to make.
      if (.not. he_diffusion .and.                                        &
          .not. transported_rows_exist()) then
         etol = 100.0d0*dble(n_species)*epsilon(1.0d0)
         dmax = 0.0d0;  jbad = 0;  ibad = 0
         do j = 1, N
            call element_nucleus_sums(f_sp(j,:), scal, res)
            call element_nucleus_sums(f_sp_old(j,:), m, etol2)
            dmax = max_defect(dmax, abs(scal - m), j, 1, jbad, ibad)
            dmax = max_defect(dmax, abs(res - etol2), j, 3, jbad, ibad)
         enddo
         if (dmax .gt. etol) then
            call step_refuse(verd, cert_step_element, 7, jbad, ibad,      &
                 dmax, etol,                                              &
                 'a cell element total moved: no operation of the step'// &
                 ' creates or destroys a nucleus')
            return
         endif
      endif

      ! ---- 3. the returned-state verdicts that exist ------------------
      ! THE CARRIER VERDICT IS READ THROUGH THE INTERVAL STATUS, not
      ! through carrier_last_verdict_of_run. A1's returned-state acceptance
      ! is what A3's retry controller acts on INSIDE the operator: an
      ! interval reported covered is one whose substeps were every one of
      ! them accepted, and an interval reported exhausted is one no substep
      ! length could get through. carrier_last_verdict_of_run is not the
      ! answer to the question "was this step's transport accepted": it
      ! holds the last verdict any interval produced, so a step in which the
      ! operator had nothing to do (no background yet, no carriers moved)
      ! reads the previous interval's verdict, or the default, which is a
      ! refusal.  carrier_step_verdict is the step-scoped companion, reset
      ! to "no interval" at the entry of every transport step.
      ! MEASURED: reading it here refused step 0 of mol_carrier, where the
      ! operator itself covered the interval and printed nothing.
      if (transported_rows_exist()) then
         if (.not. carrier_completed) then
            call step_refuse(verd, cert_step_carrier, 5, 0, 0, 0.0d0,      &
                 0.0d0, 'the carrier transport interval was not covered')
            return
         endif
      endif
      if (energy_update_last_status .ne. ENERGY_UPDATE_OK) then
         call step_refuse(verd, cert_step_energy, 9,                       &
              energy_update_last_cell, 3, energy_update_last_residual,     &
              0.0d0, 'the energy update returned a named failure')
         return
      endif
      if (conduction_last_status .ne. CONDUCTION_OK) then
         call step_refuse(verd, cert_step_conduction, 11,                  &
              conduction_last_cell, 3, 0.0d0, 0.0d0,                       &
              'the conduction stage returned a named failure')
         return
      endif
      ! A non-root chemistry (class 4 or 6) rejects a physical step (contract
      ! section 3). It
      ! is admissible only in initialization mode, where the evaluator
      ! still runs and its report is labelled initialization.
      if (run_mode .eq. run_mode_phys .and. n_no_chem_root .gt. 0) then
         call step_refuse(verd, cert_step_chem_root, 7, 0, 0,              &
              dble(n_no_chem_root), 0.0d0,                                 &
              'a cell carries a composition that is not a root'//         &
              ' of the network (class 4 or 6)')
         return
      endif

      ! ---- 4. the time-discrete hydrodynamic residual -----------------
      if (dt_g .gt. 0.0d0) then
         do j = 1, N
            do k = 1, 3
               res  = abs(u_hydro(k,j) - u_old(k,j) - dt_g*L_hydro(k,j))
               scal = abs(u_old(k,j)) + dt_g*abs(L_hydro(k,j))
               m    = res/max(scal, cert_scale_floor)
               if (m .gt. verd%hydro_row_max) then
                  verd%hydro_row_max = m
                  verd%hydro_jworst  = j
                  verd%hydro_kworst  = k
               endif
            enddo
         enddo
         verd%hydro_row_above = (verd%hydro_row_max .gt. cert_tol_hydro_step)
         verd%measure = verd%hydro_row_max
         verd%tol     = cert_tol_hydro_step
      endif
      end subroutine certification_evaluate_physical_step

      ! ------------------------------------------------------!

      subroutine element_nucleus_sums(fs, s_H, s_He)
      ! The hydrogen and helium nucleus totals of one cell, in the units of
      ! the fraction layout: H I and H II carry one hydrogen nucleus each,
      ! H2 and H2+ two, H3+ three, HeH+ one hydrogen and one helium, OH one
      ! and H2O two; He I (the triplet inside it), He II and He III carry
      ! one helium each. The species indices are species_table's.
      real*8, dimension(n_species), intent(in)  :: fs
      real*8,                       intent(out) :: s_H, s_He
      s_H  = fs(1) + fs(2)
      s_He = 0.0d0
      if (thereis_He) s_He = fs(3) + fs(4) + fs(5)
      if (thereis_mol) then
         s_H  = s_H + 2.0d0*fs(34) + 2.0d0*fs(35) + 3.0d0*fs(36) + fs(37)
         s_He = s_He + fs(37)
      endif
      if (thereis_oxychem) s_H = s_H + fs(38) + 2.0d0*fs(39)
      end subroutine element_nucleus_sums

      ! ------------------------------------------------------!

      real*8 function max_defect(dcur, d, j, k, jbad, ibad) result(dout)
      real*8,  intent(in)    :: dcur, d
      integer, intent(in)    :: j, k
      integer, intent(inout) :: jbad, ibad
      dout = dcur
      if (d .gt. dcur) then
         dout = d;  jbad = j;  ibad = k
      endif
      end function max_defect

      ! ------------------------------------------------------!

      subroutine step_refuse(verd, reason, op, j, k, m, tol, what)
      type(cert_step_verdict), intent(inout) :: verd
      integer, intent(in) :: reason, op, j, k
      real*8,  intent(in) :: m, tol
      character(len=*), intent(in) :: what
      verd%accepted  = .false.
      verd%reason    = reason
      verd%operation = op
      verd%jworst    = j
      verd%kworst    = k
      verd%measure   = m
      verd%tol       = tol
      verd%what      = what
      end subroutine step_refuse

      ! ------------------------------------------------------!

      function cert_step_reason_text(reason) result(s)
      integer, intent(in) :: reason
      character(len=64) :: s
      select case (reason)
         case (cert_step_ok);         s = 'accepted'
         case (cert_step_nonfinite);  s = 'a value was not finite'
         case (cert_step_positivity); s =                                  &
              'the state left the admissible set'
         case (cert_step_element);    s = 'an element invariant broke'
         case (cert_step_carrier);    s = 'the carrier verdict refused it'
         case (cert_step_energy);     s = 'the energy update failed'
         case (cert_step_conduction); s = 'the conduction stage failed'
         case (cert_step_chem_root);  s =                                  &
              'a cell was not a root of the chemical network'
         case (cert_step_hydro_row);  s =                                  &
              'the time-discrete hydrodynamic row was out of tolerance'
         case default;                s = 'unnamed'
      end select
      end function cert_step_reason_text

      ! ------------------------------------------------------!

      ! End of module
      end module certification

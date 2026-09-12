      module diffusive_photochemistry
      ! Vertical transport of the molecular carriers, solved together with
      ! the chemistry that makes and destroys them:
      !
      !   d n_i/dt + (1/r^2) d/dr [ r^2 ( n_i v + Phi_i ) ] = P_i - L_i
      !
      !   Phi_i = - n_tot (D_i + K_zz) d f_i/dr
      !           - n_i D_i [ 1/H_i - 1/H_atm ]        (f_i = n_i/n_tot)
      !
      ! for i = H2, OH, H2O, CO.  This is milestone M3 of
      ! docs/a2_oxygen_option_design.md; sections 3 and 4 of that document
      ! are the specification and this header states what was built and what
      ! was measured, not what was intended.
      !
      ! ---------------------------------------------------------------
      ! 1. WHY THE CARRIERS CANNOT BE A LOCAL STEADY STATE
      !
      ! Everywhere else in EXHALE the composition of a cell is the root of a
      ! local algebraic system: the partition of an element among its stages
      ! is recomputed from scratch each step by ioniz_eq.  That is correct
      ! wherever the chemistry is fast compared with the flow, and it is the
      ! reason the code can carry thirty ion stages cheaply.
      !
      ! The measurement that motivated this module says it is not correct
      ! for H2 at a cool base (design sec. 3.1, from the P4 block of
      ! docs/oxygen_chemistry_new_plan.md): at the HD 189733 b 1-microbar
      ! level tau_chem(H2)/tau_adv is 0.20 to 1.67 on H/v, and 16.4 when the
      ! two are compared at the same pressure, against 1e-3 on HD 209458 b.
      !
      ! WHAT M3 THEN MEASURED IN EXHALE'S OWN STRUCTURE, AND IT IS NOT THE
      ! SAME NUMBER.  That ratio was built from the PHOTOCHEMICAL model's
      ! scale height and velocity at that pressure.  In an EXHALE run of the
      ! same planet the base cell has v = 0 EXACTLY -- the lower boundary
      ! condition puts it there -- and the flow time r/|v| in the cells above
      ! it is 1e8-2e8 s against tau_chem(H2) = 178-552 s.  The diffusive time
      ! of the base cell is 8.0e5 s with K_zz = 0, i.e. 4500 chemical times,
      ! and 634 s with K_zz = 1e9, i.e. 3.6.  So in this code's structure the
      ! base partition is a LOCAL quantity, and transport is not the term
      ! that decides it; what transport does decide is the profile of a
      ! species with no chemistry of its own, which here is CO.  The M3
      ! result block of the design document carries the numbers and what
      ! follows from them.
      !
      ! That is why every run writes tau_chem, tau_adv AND the diffusive time
      ! of each cell to output/Oxygen_chemistry.txt: the regime is a
      ! measurement of the run, not a property of the option.
      !
      ! ---------------------------------------------------------------
      ! 2. WHICH SPECIES ARE TRANSPORTED, AND WHY THE OTHERS ARE NOT
      !
      ! Transported: H2, OH, H2O, CO.  Between them they carry every H
      ! nucleus that is not atomic or in a molecular ion, and every O and C
      ! nucleus that is not a free atom.  OH is in the set even though its
      ! lifetime is short: measured over the H2 layer of the HD 189733 b
      ! run, OH lives 4.0e-3 to 0.40 s against O1 plus its photolysis while
      ! H2 lives 33 to 553 s against the oxygen cycle, i.e. four decades
      ! longer.  It is carried because it holds up to 14% of the oxygen
      ! where the carriers live (against 43% for H2O and 55% for CO), so
      ! leaving it out would open a hole in the oxygen closure, and because
      ! an implicit step carries a short-lived species at no cost: it simply
      ! relaxes to its own local balance.
      !
      ! NOT transported, and closed by their local balance instead:
      !  - H2+, H3+, HeH+.  Measured on the HD 209458 b molecular example
      !    over the cells that hold H2 (r = 1.000-1.018 R_p), the three of
      !    them together carry at most 1.2e-9 of the H nuclei.  Leaving them
      !    local is therefore a closure statement at that level and not a
      !    leak: the H-nucleus budget is closed against the ELEMENT total,
      !    not against a sum of transported species, so whatever they hold
      !    is simply not available to the carriers.
      !  - O(1D).  Its lifetime is 1.6e-3 s at the HD 189733 b base and it
      !    has one sink; the design's sec. 2.3 states the case, and it is the
      !    same argument as above rather than a different one.
      !  - H I / H II and the ion stages of every element.  Their balance is
      !    local, and gate G2 requires that limit to be reproduced.  They are
      !    where the transported carriers' nuclei come from and go to: the
      !    write-back below moves nuclei between a carrier and the closure
      !    species of its element, so no element can drift.
      !
      ! The one term this drops that is not small everywhere is the AMBIPOLAR
      ! field.  All four transported species are neutral, so eq. (2) of the
      ! design has no (Z e E)/(k T) term for any of them; the field is not
      ! neglected, it is absent.  Had the molecular ions been transported it
      ! would not be, and settling_coefficient in binary_element_diffusion
      ! already builds it.
      !
      ! ---------------------------------------------------------------
      ! 3. THE COEFFICIENTS
      !
      ! MOLECULAR DIFFUSION D_i.  Blanc's law over the background carriers
      ! (H I, H II, H2, H2+, H3+, He I, He II, He III, HeH+), each pair from
      ! the rigid-sphere coefficient of Banks & Kockarts (1973) that
      ! binary_element_diffusion already carries
      ! (hard_sphere_pair_diffusion, reused rather than restated).
      !
      ! What that leaves out is the ION-NEUTRAL polarization enhancement.
      ! For an ion colliding with one of these molecules the induced-dipole
      ! (Langevin) cross section is larger than the rigid-sphere one, so the
      ! rigid-sphere D is an UPPER bound on that pair.  Two reasons it is
      ! used anyway, and both are measurements rather than preferences.
      ! First, the transported species only exist where the gas is neutral:
      ! over the cells that hold H2 on the HD 209458 b molecular example the
      ! electron fraction per H nucleus is at most 1.0e-7 (3.1e-11 at the
      ! base cell), so the ion pairs carry that share of the Blanc friction
      ! sum and the choice of coefficient for them is numerically
      ! irrelevant.  Second, for H2O the induced-dipole formula would be the
      ! wrong one in any case: a permanent dipole of 1.85 D dominates the
      ! ion-molecule capture rate, and the locked-dipole treatment that needs
      ! is not in this code.  Putting an induced polarizability there would
      ! look more careful and be less honest.
      !
      ! EDDY DIFFUSION K_zz.  Read from kzz_cell and nowhere else.  This
      ! module defines no eddy coefficient and no key for one (design
      ! sec. 3.5): kzz_cell is filled once by eddy_diffusion_on_grid, from
      ! the scalar He_Kzz or from the lower-atmosphere profile's own K_zz
      ! column, and every diffusive term in the code reads that one array.
      ! With no profile and no He_Kzz it is zero, and the transport is then
      ! pure molecular diffusion -- which on a lower atmosphere is the wrong
      ! limit, since eddy mixing is what holds the composition well mixed
      ! below the homopause.  The run says so at startup rather than
      ! refusing, so that the K_zz = 0 limit stays available as a test.
      !
      ! SETTLING.  G_i = (m_i - m_bar) g/(k T), with m_bar the mean mass per
      ! gas particle (rho/n_tot) and g from the same Dphi the hydro uses.
      ! Thermal diffusion is not carried (alpha_T = 0): there is no evaluated
      ! alpha_T for H2O, OH or CO in an H2/He bath, and the element operator
      ! makes the same choice for the same reason.
      !
      ! ---------------------------------------------------------------
      ! 4. THE CHEMISTRY IS THE SAME CHEMISTRY
      !
      ! P_i - L_i is not written here.  It is obtained by calling the very
      ! rows the local equilibrium solve uses -- mol_heh_rows row 4 for H2
      ! and oxygen_carrier_rows for OH and H2O -- at the trial densities,
      ! with the background frozen.  A rate that changes there changes here,
      ! and the two cannot drift apart the way two transcriptions of one
      ! coefficient do (the R16-R20 lesson recorded in mol_rates.f90).
      ! The Jacobian is a forward difference of those same calls.
      !
      ! CO has no chemical source: it is frozen chemically by decision D4 and
      ! only its transport moves it.  What that decision MEANS changes when
      ! CO is transported, and the change is not cosmetic.  In the local
      ! solve "inert" meant "at its own CO <-> C + O chemical equilibrium",
      ! which dissociates CO above about 4000 K.  Transported, "inert" would
      ! mean INDESTRUCTIBLE, and it was: measured on HD 189733 b before the
      ! ceiling of limit_to_element_budget was added, the wind carried CO out
      ! to 1.67 R_p and 2e4 K and it held 55% of the oxygen there, switching
      ! off the O I and C II cooling of the entire wind.  The ceiling is the
      ! statement the thermodynamics supports and the audited set has no rate
      ! for; its own limits are written at that code site.
      !
      ! WHAT IS FROZEN INSIDE ONE STEP: the ion stages, the molecular ions,
      ! n_e, n_tot, T, and the photoionization rates.  The photolysis rates
      ! are NOT frozen -- they are rebuilt from the current H2O, OH and H2
      ! columns at the top of every step, so the self-shielding of the water
      ! layer responds within the relaxation instead of only between Picard
      ! passes.
      !
      ! ---------------------------------------------------------------
      ! 5. BOUNDARY CONDITIONS (design decision D6)
      !
      ! WHETHER THE BASE FACE IS AN INFLOW is decided by the wind's mass
      ! flux and not by the base cell's velocity.  Measured on the converged
      ! hot Uranus (docs/p44_base_sawtooth.md): the cell-centred rho v r^2 of
      ! cell 1 is 160 to 200 times the wind's flux INWARD while the face flux
      ! is 0.36 to 0.98 times the wind's flux OUTWARD, and the discrepancy
      ! decays over four cells, which is the signature of a collocated
      ! odd-even velocity mode and not of infall.  The stationary evaluation
      ! reads this; the marching rows read the face mass flux itself, which
      ! carries the same statement without needing a rule.
      !
      ! ELEMENT RESERVOIRS ARE DIRICHLET; THE PARTITION AMONG CARRIERS IS
      ! NOT.  The H, He, O and C totals at the base are boundary data the
      ! handoff legitimately supplies, and they are imposed where they always
      ! were -- through melem_ab, the base density and the EOS.  The split of
      ! each element among its carriers is what this option exists to
      ! compute, so it is not imposed anywhere: the base face carries zero
      ! diffusive flux, and the advective inflow, which is now the base face
      ! mass flux times the ghost composition, carries the base cell's own
      ! partition because the ghost is given that partition (df/dr = 0 at the
      ! face), which is the same statement.
      !
      ! That is what makes the A/B gate of design sec. 6.2 a test.  Had the
      ! base H2 fraction been pinned from the same handoff the gate compares
      ! against, the gate would measure nothing.
      !
      ! OUTER: zero gradient, as the element operator does.  The carriers are
      ! negligible there by many decades in every configuration this option
      ! targets; a run in which they are not is outside its validity range
      ! and output/Oxygen_chemistry.txt says so.
      !
      ! ---------------------------------------------------------------
      ! 6. DISCRETIZATION AND THE SIMPLEX
      !
      ! Cell-centred f_i, face-centred fluxes on r_edg, one backward-Euler
      ! step per call, Newton with a backtracking line search, block
      ! tridiagonal in space (4x4 blocks: the transport is diagonal in
      ! species and the chemistry is dense within a cell).  Faces 0 and N
      ! carry zero diffusive flux.  The drift term switches from central to
      ! donor-cell on the same Peclet test the element operator uses.
      !
      ! THE MATERIAL ADVECTION OF A CARRIER IS NOT IN THESE ROWS.  It is the
      ! divergence of the hydrodynamic face mass flux, F_rho(j) Y_c(j)/m_c,
      ! and it is taken with the mass row itself inside the Runge-Kutta
      ! stages, on the same faces, areas, volumes and time step
      ! (species_face_flux.f90, docs/b4_spatial_operator_design_20260906.md
      ! T-B4.1 to T-B4.3).  A cell can then only lose the carrier that the
      ! mass row says it loses, whatever the velocity field does, and a
      ! uniform partition is preserved to round-off.  What these rows carry
      ! is the operator-split remainder: the diffusive and drift face fluxes
      ! and the chemistry.  The one place the advective term still appears is
      ! carrier_steady_residual, because the balance a stationary state
      ! satisfies is the whole equation and not the remainder of a split.
      !
      ! WHY THE ADVECTIVE HALF HAS TO BE ON THE FACE MASS FLUX, and it is a
      ! numerical statement about the front this operator produces.  A
      ! first-order donor-cell difference on the cell velocity carries a
      ! numerical diffusivity |v| dr/2.  Measured on the 500-cell grid of the
      ! hot Uranus that is 8.5, 20 and 22 times the molecular D of H2 at 1.15,
      ! 1.30 and 1.50 R_p, and a controlled grid test says it is not a
      ! bookkeeping detail: refining the H2 front by a factor two at UNCHANGED
      ! CFL step (500 -> 1000 cells, the base spacing untouched, so the global
      ! dt moves by 2%) cut the front's advance rate by a factor five and
      ! removed its re-acceleration entirely.  The front this operator
      ! produced was being carried by the scheme as much as by the flow.  The
      ! face form reconstructs the composition with the SAME limiter as the
      ! primitive variables, so it is second order where the profile is
      ! smooth and it is the hydrodynamics' own accuracy rather than a
      ! separate correction with its own limiter.
      !
      ! POSITIVITY AND THE ELEMENT SIMPLEX.  f_i >= 0 is enforced on every
      ! Newton update.  After the step the carriers are limited so that the
      ! nuclei they hold do not exceed the element's own total, in H, O and
      ! C alike; every limiter is counted and reported, never silent.  The
      ! write-back then shares the remaining nuclei out over the element's
      ! other species in proportion to what they already held.
      !
      ! WHAT THAT MAKES THIS OPERATOR, STATED PLAINLY.  Each cell's ELEMENT
      ! TOTALS are unchanged by the step: only the partition among an
      ! element's carriers moves.  So element conservation is exact by
      ! construction rather than by cancellation -- which is what gate G5
      ! measures -- and rho is untouched, because bsp_mass of a carrier is
      ! exactly the sum of its nuclei's masses.
      !
      ! It also fixes the approximation.  A water molecule that diffuses out
      ! of a cell really takes its oxygen with it; here the oxygen stays and
      ! reappears as O I.  That is the same statement EXHALE already makes
      ! about every metal when trace-metal diffusion (`He_metal_diffusion`) is off -- the
      ! element is slaved to hydrogen at melem_ab -- and it is why this
      ! operator transports SPECIATION and not elements.  Elemental
      ! separation of oxygen is he_metal_diffusion's subject, and the two
      ! are refused together until CO, which carries an oxygen AND a carbon
      ! nucleus, has a rule for following two element factors at once
      ! (input_read states the refusal).  For hydrogen the distinction does
      ! not arise: the H total is the element, so transporting the partition
      ! 2 n_H2/n_H IS transporting what the gate measures.
      !
      ! ---------------------------------------------------------------
      ! 7. VALIDITY
      !
      ! The same range the reaction set has (oxygen_rates header): this is
      ! valid where the base H2/H partition is kinetic and set by the oxygen
      ! cycle.  Transport moves the boundary of validity outward, it does not
      ! remove it -- design sec. 9 -- and the Damkohler profile in
      ! output/Oxygen_chemistry.txt is where each run states its own.
      !
      ! References: Banks & Kockarts (1973), Aeronomy, Part B, ch. 15 (the
      ! diffusive flux and the rigid-sphere coefficient); Draine & Bertoldi
      ! (1996) for the H2 band; the reaction set is
      ! docs/a2_reaction_audit.md.

      use global_parameters
      use caloric_eos, only: adiabatic_index_at_T
      use grav_func,     only: Dphi
      use species_table, only: n_bsp, bsp_fsp, bsp_mass,                 &
                               isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                               isp_HeIII, isp_HeTR,                      &
                               isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                               isp_OH, isp_H2O, isp_CO,                  &
                               n_mion, mion_fsp, melem_i0, melem_top,    &
                               iel_O, iel_C
      use binary_element_diffusion, only: hard_sphere_pair_diffusion,     &
                          advected_carrier_reset,                        &
                          advected_carrier_register,                     &
                          advected_carrier_count,                        &
                          advected_carrier_fractions,                    &
                          advected_carrier_state_ready,                  &
                          advected_carrier_state_consumed,               &
                          carrier_mass_fractions, mixture_mass_sum
      use ion_cell_state, only: ieq_cell, ion_rates
      use oxygen_rates, only: rk_D1_Hep_CO,                              &
                              rk_CO_radiative_association
      use co_self_shielding_table, only: co_shield_tex_limit_K
      use System_HeH_mol, only: set_mol_coeffs, set_oxygen_coeffs,       &
                                mol_heh_rows, oxygen_carrier_rows,       &
                                oj3, oj4, oj5, oj7
      use water_photolysis, only: n_fuv_band
      use utils_ion_eq, only: fuv_lw_photon_field
      use utils, only: calc_ne
      use ionization_equilibrium, only: bg_cell, bg_ready, finite_real,   &
                                        ioniz_eq, ioniz_eq_ledger
      use ion_cell_state, only: ion_rates
      use caloric_eos, only: pressure_from_energy_density
      use steady_residual_mod, only: carrier_row_scale,                  &
                                     face_mass_flux_of_state
      use species_advective_transport, only: species_face_fraction,      &
                                     species_face_flux,                  &
                                     species_flux_divergence
      use composition, only: base_h2_composition_imposed,               &
                             get_species_densities, comp_T_from_p
      use element_census, only: element_census_state, element_census_take,&
                               element_census_verify
      ! The one stoichiometric map of the H, He, O and C nuclei of a cell.
      ! The element densities this operator spends, the hydrogen budget its
      ! carriers may take and the box side of each carrier are quantities
      ! derived from that map, so the stoichiometry behind them is written
      ! in exactly one place.
      use element_inventory, only: carrier_element_totals,                &
                                   carrier_hydrogen_budget,              &
                                   element_box_side,                     &
                                   ien_H, ien_O,                         &
                                   element_inventory_report,             &
                                   element_inventory_report_on

      implicit none
      private

      public :: photochemical_transport_step
      public :: relax_photochemical_composition
      ! WHAT ENDED A CARRIER RELAXATION PASS, in words.  The pass advances
      ! the carriers at a fixed wind under a bound on how far the
      ! composition may move, and its five endings are physically different
      ! states of affairs and not degrees of one, so the caller reads which
      ! one instead of reading success into a stop.  The named values are
      ! the carrier_relax_* parameters below.
      public :: carrier_relax_outcome_text
      public :: carrier_transport_diagnostics
      public :: carrier_diffusion_coefficient
      public :: carrier_set_init
      public :: n_carrier_max
      ! Element-closure unit test (src/tests/element_census_tests.f90): the
      ! write-back's invariant is that the element totals it is handed come
      ! back unchanged, and the test states it on the two routines that make
      ! the pair.
      public :: carrier_state, carrier_write_back
      ! The advective term of the carrier rows, exposed so that an
      ! acceptance test measures it against the Runge-Kutta stage's own
      ! update of the same state and not against a transcription of it.
      public :: carrier_advective_divergence, carrier_mass_amu
      ! Save and restore of the module arrays carrier_steady_residual
      ! overwrites, so that the balance of a state can be MEASURED without
      ! the measurement being an operation on the state.
      public :: carrier_module_state
      public :: save_carrier_module_state, restore_carrier_module_state
      public :: carrier_module_state_matches
      public :: carrier_steady_residual
      public :: carrier_drift_location
      ! The two row scales of the last assembly side by side, so that
      ! src/tests/carrier_retry/ can state on the production arrays what
      ! each of them does to the measure when the step is shortened: the
      ! one the Newton iterates against carries the time term, the one a
      ! returned state is judged against does not.
      public :: carrier_row_scales
      ! Exposed so an acceptance test reads the mass this module gives a
      ! collision partner, and not a second copy of it.
      public :: species_mass_amu, carrier_background_mass
      ! The coupled steady solve (sec. 139): n(H2) is a Newton unknown, so
      ! the solver needs the carrier row, the scale of its unknown, and the
      ! element headroom it must stay inside.
      public :: carrier_column_scale, carrier_row_term_scale,           &
                carrier_headroom, carrier_headroom_known,               &
                carrier_element_headroom,                               &
                carrier_species_index
      ! Test only: writes the element budget of a synthetic cell, so that
      ! src/tests/steady_species_rows/ can state the projection onto the
      ! budget face of the coupled solve's unknown box.
      public :: carrier_headroom_set_for_test
      public :: n_carrier, carrier_name, ic_H2, ic_OH, ic_H2O, ic_CO
      public :: ic_Hp, carrier_solved
      ! The element table of the carriers and the perturbation the chemistry
      ! rows are differentiated with: public because src/tests/physics_probe/
      ! carrier_reference_scales.f90 measures the directional derivative the
      ! second one produces against a central difference of the first.
      public :: carrier_element_reference_density,                       &
                carrier_source_derivative_step, carrier_source
      ! The domain of the one-sided CO destruction model, cumulative over
      ! the whole run: the cells in which tau_dest << tau_res << tau_form
      ! fails, which is the ordering that makes a CO row carrying
      ! destruction and no formation legitimate.
      public :: carrier_co_domain_record, carrier_co_domain_f_dom
      ! Test only: writes distinct amounts into the two ledger families of
      ! the domain record, so that src/tests/carrier_retry/ can state on the
      ! production arrays that the record separates the families, that the
      ! run total composes them (a sum for the counts, a maximum for the
      ! worst ratio), and that a refused attempt does not roll it back.
      public :: carrier_co_domain_perturb_for_test
      ! Whether an element constraint moved a given cell in the last
      ! limiter call: src/tests/carrier_constraint_attribution/ measures
      ! that the CO thermal ceiling marks the cell it acts in, which is the
      ! cell whose OH and H2O rows it moves.
      public :: carrier_cell_is_constrained, limit_to_element_budget
      ! ACCEPTANCE OF THE STATE THE CARRIER SOLVE RETURNS, and the
      ! classification of a stalling iteration it keeps separate from it.
      ! Both are decisions taken on their arguments alone, so
      ! src/tests/carrier_returned_state_acceptance/ can state them on
      ! constructed rows, and so the retry controller of rev 3 sec. 4.1 can
      ! ask for a verdict on a trial state without running a solve.
      public :: carrier_returned_state_verdict, carrier_stall_class
      ! TWO RECORDS WITH TWO SCOPES, each named for the one it keeps.
      ! carrier_step_verdict is THIS transport step's: it is reset to
      ! "no interval" at the entry of photochemical_transport_step, so a
      ! caller asking whether this step's transport was accepted can never
      ! read a previous step's answer.  carrier_last_verdict_of_run is the
      ! last verdict any interval of the run produced and outlives the step
      ! that produced it, which is what makes it the wrong gate for a step.
      public :: carrier_verdict_reason_text
      public :: carrier_step_verdict, carrier_last_verdict_of_run
      ! THE CARRIER-LOCAL CHECKPOINT AND THE RETRY CONTROLLER (PLAN
      ! 20260906 rev 2, step A3).  The checkpoint is the carrier fractions
      ! plus every module array and counter one transport attempt writes;
      ! carrier_transport_interval is the step's body with the stop taken
      ! out of it, so src/tests/carrier_retry/ can reach the exhaustion
      ! branch and read what the controller did.
      public :: carrier_checkpoint
      public :: carrier_checkpoint_take, carrier_checkpoint_restore
      public :: carrier_checkpoint_matches
      public :: carrier_transport_interval
      ! WHETHER THE INTERVAL WAS COVERED, as a value the caller reads
      ! instead of inferring it from a stop (PLAN 20260906 rev 2, step
      ! A3b).  photochemical_transport_step returns it, so the step
      ! controller around this operator can reject the whole attempted
      ! step on an exhausted carrier interval rather than the run ending
      ! inside the operator.
      public :: carrier_interval_covered, carrier_interval_exhausted
      public :: carrier_last_interval_status
      ! THE EXHAUSTED INTERVALS OF AN INITIALIZATION RUN, and the mark
      ! that the carrier history of the run contains one.  In
      ! initialization mode a non-root iterate is a numerical guess and
      ! the march continues from the entry carriers
      ! (a0_run_mode_contract_20260906 sec. 2), so the record is what is
      ! left of the failure and the certification's carrier entry reads
      ! the mark.
      public :: carrier_exhausted_record, carrier_history_certifiable
      public :: carrier_attempt_record
      public :: carrier_last_interval_halvings
      public :: carrier_roundoff_limited_record
      public :: carrier_refused_rows, carrier_refused_row
      public :: carrier_record_refused_rows
      public :: carrier_perturb_checkpointed_state_for_test
      public :: carrier_reject_leading_attempts_for_test
      public :: carrier_history_reset_for_test
      public :: carrier_initial_substep_for_test
      public :: carrier_transport_stop_on_failure
      public :: carrier_transport_stops_suppressed

      integer, parameter :: dp = kind(1.0d0)

      ! HOW A CARRIER NEWTON ENDED.  A solve that stopped short returns a
      ! state that satisfies no equation, so the outcome is a value the
      ! caller reads and not a line in a log.
      integer, parameter, public :: carrier_solve_converged          = 0
      integer, parameter, public :: carrier_solve_iteration_cap      = 1
      integer, parameter, public :: carrier_solve_line_search_failed = 2
      ! The iteration stopped falling while its residual was still above
      ! the absolute floor.  It is a termination, never an acceptance: the
      ! state it leaves is judged by carrier_returned_state_verdict like
      ! any other.
      integer, parameter, public :: carrier_solve_stagnated          = 3
      ! The iteration drove its residual down to the round-off of the row
      ! terms and the state is STILL above what the acceptance asks of it
      ! on its physical terms.  That is not a failure of the iteration: at
      ! a substep far below the physical time scale of the row the time
      ! term is the largest term of the row, so the absolute residual an
      ! exact solve leaves, round-off times that term, is 1/dt times the
      ! physical terms the state is judged against.  Nothing further is
      ! reachable HERE, and a shorter substep raises the floor rather than
      ! lowering it (carrier_row_roundoff).  The rows in that condition are
      ! the ones the verdict counts as round-off limited and accepts.
      integer, parameter, public :: carrier_solve_unreachable_at_substep = 4

      ! HOW A STALLING ITERATION STOPPED.  The three are separated because
      ! only one of them is an accuracy statement: reaching the absolute
      ! floor newton_floor says the residual is as small as this
      ! discretization resolves, while falling by newton_drop RELATIVE to
      ! the residual the iteration started from says nothing about the size
      ! of what is left (a start of 1e8 admits a returned residual of 1e2).
      ! The second is therefore a diagnostic label and not a verdict.
      integer, parameter, public :: carrier_stall_none             = 0
      integer, parameter, public :: carrier_stall_at_floor         = 1
      integer, parameter, public :: carrier_stall_at_relative_drop = 2
      integer, parameter, public :: carrier_stall_open             = 3

      ! WHY A RETURNED CARRIER STATE WAS REFUSED.  carrier_accept_ok is the
      ! only value an accepted state carries.
      integer, parameter, public :: carrier_accept_ok           = 0
      integer, parameter, public :: carrier_reject_row_residual  = 1
      integer, parameter, public :: carrier_reject_nonfinite     = 2
      ! Refused by the test hook below and by nothing else, so a forced
      ! rejection can never be read in a log as a physical one.
      integer, parameter, public :: carrier_reject_forced_for_test = 3
      ! NO INTERVAL WAS INTEGRATED IN THIS STEP, so there is no state to
      ! pass a verdict on.  It is the value carrier_step_verdict carries in
      ! a step the operator had nothing to do in (a run with no carriers,
      ! or one whose background is not ready), and it is neither an
      ! acceptance nor a refusal: accepted stays .false., because nothing
      ! was accepted.
      integer, parameter, public :: carrier_no_interval = 4
      ! A ROW AT THE ROUND-OFF OF ITS OWN TERMS IS NOT A REFUSAL REASON.
      ! Such a row carries the residual an exact solve of it leaves, so it
      ! is no evidence of a defect in the returned state: it is the limit
      ! of double precision at that substep.  It is accepted and counted as
      ! round-off limited (n_rows_roundoff_limited), and there is therefore
      ! no reason value for it.

      ! THE VERDICT ON THE STATE THE SOLVE HANDS BACK, which is the state
      ! the write-back puts into the run.  It is decided on that state's own
      ! rows and not on how the iteration ended: the raw status and the
      ! stall class are carried along for diagnosis only.  A rejection and
      ! its retry (rev 3 sec. 4.1, A3 increment (b)) consume this type.
      type, public :: carrier_verdict
         logical  :: accepted   = .false.
         integer  :: reason     = carrier_reject_row_residual
         ! The unconstrained row that decided it, and its floored relative
         ! imbalance |R| / (sum of the row's own term magnitudes + floor).
         integer  :: jworst     = 0
         integer  :: icworst    = 1
         real(dp) :: rworst     = 0.0d0
         ! What the element constraints held, so that a rejection can never
         ! be read as a constraint failure and a constrained acceptance
         ! always says how much of the state it did not certify.
         integer  :: n_cells_constrained = 0
         integer  :: n_rows_constrained  = 0
         ! How many unconstrained rows stand above the floor AND carry more
         ! imbalance than the arithmetic that formed them explains, and how
         ! many of those sit in a cell NEXT TO a constrained one.  A clamp
         ! moves the carrier densities of its own cell, and those enter the
         ! face fluxes of both its neighbours, so a refusal whose rows are
         ! all neighbours of clamped cells is the constraint reaching one
         ! cell further, while a refusal away from every clamp is the
         ! solve's.  The two are counted rather than argued.  These are the
         ! rows a refusal rests on: n_rows_above = 0 is an acceptance.
         integer  :: n_rows_above        = 0
         integer  :: n_rows_above_beside_constrained = 0
         ! How many rows stand above the floor with a residual that is
         ! already the round-off of their own full terms, which is what an
         ! exact solve leaves there.  They are ACCEPTED: nothing in such a
         ! row says the state is wrong, only that double precision resolves
         ! it no further at this substep.  A trace carrier whose physical
         ! terms sit at the 1e-20 element floor of carrier_residual is
         ! typically here.  Informational, and reported as such.
         integer  :: n_rows_roundoff_limited = 0
         ! Diagnosis, never the decision.
         integer  :: raw_status = carrier_solve_converged
         integer  :: stall      = carrier_stall_none
      end type carrier_verdict

      ! The transported set, in the order the solve uses.
      !
      ! WHICH SPECIES ARE IN IT IS A PROPERTY OF THE RUN, NOT OF THE MODULE.
      ! H2 is carrier 1 in every configuration -- it is the one carrier
      ! whose Damkohler number is of order unity in the region the option
      ! exists for -- and the oxygen cycle adds OH, H2O and CO behind it.
      ! Measured on the converged He/H = 0.0793 hot Uranus, the H2 chemical
      ! time is 0.02 to 1.3 times the flow time r/|v| across the front,
      ! while the three molecular ions that are left local sit at 1e-9 to
      ! 5e-5 of it and hold at most 2e-8 of the H nuclei, so the split
      ! between carried and local is a measurement and not a convention
      ! (docs/supersonic_molecular_base.md sec. 13.6).
      !
      ! n_carrier_max dimensions every array; n_carrier is how many of them
      ! this run solves, set once by carrier_set_init.  The inactive columns
      ! stay at zero -- the species behind them do not exist without the
      ! oxygen chemistry -- so the whole-array statements below are exact.
      !
      ! THE PROTON IS THE FIFTH, AND IT IS LAST FOR A REASON. Its index has
      ! to be a constant -- the ionization rows, the write-back and the
      ! element budget all name it -- and the four existing indices have to
      ! keep the values they have, so the only place left is behind them.
      ! That makes the solved set non-contiguous in the one configuration
      ! that carries H2 and H+ without the oxygen cycle, which is why
      ! carrier_solved exists below: n_carrier is the LAST index solved and
      ! carrier_solved says which of 1..n_carrier actually are.
      integer, parameter :: n_carrier_max = 5
      integer, save      :: n_carrier = 4   ! set by carrier_set_init
      integer, parameter :: ic_H2 = 1, ic_OH = 2, ic_H2O = 3, ic_CO = 4
      integer, parameter :: ic_Hp = 5
      character(len=3), parameter :: carrier_name(n_carrier_max) =       &
           (/ 'H2 ', 'OH ', 'H2O', 'CO ', 'H+ ' /)
      ! Which of 1..n_carrier this run solves. A carrier that is not solved
      ! keeps the identity row nrho/dt on the diagonal and a zero residual,
      ! so its unknown does not move and the block structure the Thomas
      ! sweep needs is unchanged; it costs one 5x5 block instead of one 2x2
      ! and no chemistry evaluation. Set once by carrier_set_init.
      logical, save :: carrier_solved(n_carrier_max) = .true.

      ! Background collision partners for Blanc's law: the species whose
      ! friction a transported carrier feels, named by their f_sp column.
      ! Their MASSES are not listed here; carrier_background_mass reads each
      ! from the species table, so the reduced mass of every pair is formed
      ! from the same number calc_rho weighs the species with and the
      ! friction and the density cannot disagree about what a helium atom
      ! weighs.
      integer, parameter :: n_bkg = 9
      integer, parameter :: bkg_isp(n_bkg) =                             &
           (/ isp_HI, isp_HII, isp_H2, isp_H2p, isp_H3p,                 &
              isp_HeI, isp_HeII, isp_HeIII, isp_HeHp /)

      ! Atomic mass unit [g], the value binary_element_diffusion uses.
      ! The mass unit of the density normalization is the hydrogen ATOM mass mu
      ! of parameters.f90 (one definition); this module reads it and keeps no copy.

      ! Newton controls, mirroring solve_mass_fraction, except for the
      ! iteration budget. The damped line search crosses a stiff transient
      ! (an H2 destruction time shorter than the step, after the base
      ! ionization jumps) in steps of about one percent of the merit before
      ! the quadratic basin is reached: measured 141 iterations to 3.9e-12
      ! on the oxygen_chemistry case where 30 stopped the run, against 2 to
      ! 5 iterations on every ordinary step. The budget therefore has to
      ! cover the globalization, not only the quadratic tail.
      integer,  parameter :: newton_maxit  = 200
      real(dp), parameter :: newton_tol    = 1.0d-12
      real(dp), parameter :: newton_floor  = 1.0d-8
      real(dp), parameter :: newton_drop   = 1.0d-6
      integer,  parameter :: newton_halves = 8
      ! WHAT AN ACCEPTED RETURNED STATE MUST SATISFY, row by row.  The
      ! measure is the row residual divided by that row's own PHYSICAL
      ! scale: the sum of the magnitudes of the terms the row balances that
      ! belong to the state, the transport flux divergence and the reaction
      ! production and loss, PLUS an absolute floor of 1e-20 of the free
      ! density of the carrier's element in the row's volumetric-rate units
      ! (carrier_residual states the floor and why it is there).  The time
      ! term is deliberately NOT in it: it is 1/dt times a density, so it
      ! belongs to the integrator and not to the state, and a measure that
      ! carried it would fall like dt for one fixed physical imbalance
      ! (declaration of row_terms_phys).  So this is a
      ! relative condition with an absolute floor already inside it: a
      ! carrier at 1e-28 of its element cannot fail it on round-off, and a
      ! carrier that holds the layer cannot pass it by being compared with
      ! a bound instead of with its own terms.
      !
      ! The number is newton_floor, the absolute residual the solve already
      ! declares to be what this discretization resolves.  It is the same
      ! floor on the same measure, asked of the state that is handed back
      ! rather than of the iteration's last trial.
      real(dp), parameter :: carrier_accept_tol = newton_floor
      ! THE SMALLEST ROW RESIDUAL THIS ASSEMBLY CAN TELL FROM ZERO, as a
      ! fraction of the row's FULL terms.  A carrier row is a sum of the
      ! time term, two face fluxes, an advective term, the deferred
      ! second-order correction and the chemistry, each of them built from
      ! several products, so forming it commits an accumulated rounding of
      ! order a few tens of eps times the largest term in it, and that
      ! largest term is what row_terms measures.  A residual at or below
      ! this fraction of row_terms is what an EXACT solve of the row would
      ! leave: no further iteration reduces it, and halving the substep
      ! raises it, because the time term it is round-off of grows like
      ! 1/dt while the physical terms the state is judged against do not.
      ! MEASURED on the molecular column of src/tests/carrier_retry/: the
      ! Newton reaches 2.6e-16 of the full row terms, about one eps, so the
      ! bound is not tight against the iteration it classifies.
      !
      ! WHAT IT DECIDES.  A row whose residual is at or below this fraction
      ! of its full terms is ACCEPTED and counted as round-off limited: the
      ! residual left in it is the arithmetic of forming the row, so it is
      ! no evidence of a defect in the state.  A trace carrier whose
      ! physical terms sit at the 1e-20 element floor of carrier_residual
      ! is typically in that condition and is accepted with the flag.  For
      ! every row whose physical terms stand above their own round-off the
      ! acceptance is unchanged and stays independent of the substep: it is
      ! the ratio of the residual to the physical terms, and neither of
      ! those carries dt.
      real(dp), parameter, public :: carrier_row_roundoff =                &
                                     64.0d0*epsilon(1.0d0)
      ! Relaxation controls, mirroring relax_element_composition.
      integer,  parameter :: relax_maxstep = 400
      real(dp), parameter :: relax_tol     = 1.0d-10
      ! THE SHORTEST TRIAL STEP A RELAXATION PASS WILL TAKE, as a fraction
      ! of the cell's own transport time.  A trial that leaves the movement
      ! bound is halved and retried, and at this length it can add at most
      ! about 1e-3 of the bound: the first step at the full cell crossing
      ! time moves the largest mixing ratio by 1.1e-2 on the hot Uranus
      ! column, and the movement of a step falls with its length.  That is
      ! below the tolerance the bound is stated to, so a trial that still
      ! leaves the bound here says the bound is reached and not that the
      ! step is too long.
      real(dp), parameter :: relax_grow_min = 2.0d0**(-10)
      ! THE ENDINGS OF A RELAXATION PASS.  MOVEMENT_BOUND: the returned
      ! state lies inside the bound and the shortest admissible trial would
      ! leave it.  FIXED_POINT: a full-length step no longer moves the
      ! state.  STEP_BUDGET: relax_maxstep trials were taken without
      ! either.  INTERVAL_REFUSED: the shortest admissible trial still could
      ! not be covered by the transport step, so the state inside the bound
      ! is what stands.  NO_CARRIERS: the
      ! configuration transports no carrier, or no frozen background has
      ! been formed yet, so there is nothing to advance and the unmoved
      ! state is inside every bound.
      integer, parameter, public :: carrier_relax_movement_bound   = 0
      integer, parameter, public :: carrier_relax_fixed_point      = 1
      integer, parameter, public :: carrier_relax_step_budget      = 2
      integer, parameter, public :: carrier_relax_interval_refused = 3
      integer, parameter, public :: carrier_relax_nothing_to_advance = 4
      ! The chemistry of a kept transport step did not close (a sweep that
      ! left a cell non-finite), and the shortest admissible trial did no
      ! better: the entry state is handed back.
      integer, parameter, public :: carrier_relax_chemistry_refused = 5

      ! The background this step holds frozen: the ion stages, the molecular
      ! ions and the free electrons (sec. 4).  Filled once per step by
      ! carrier_state from the f_sp the step was handed, so the chemistry
      ! rows see a state and not a moving target.
      real(dp), dimension(:), allocatable :: cbg_nhii, cbg_nh2p
      real(dp), dimension(:), allocatable :: cbg_nh3p, cbg_nhehp
      real(dp), dimension(:), allocatable :: cbg_nhei, cbg_nheii
      real(dp), dimension(:), allocatable :: cbg_nheiii, cbg_nheiTR
      real(dp), dimension(:), allocatable :: cbg_ne
      ! Ionized stages of oxygen and carbon, frozen with the rest of the
      ! background: they are what separates the element headroom (which the
      ! limiter uses) from the FREE ATOMIC density the reaction rows need.
      real(dp), dimension(:), allocatable :: cbg_nOion, cbg_nCion
      ! Photolysis rates rebuilt from the CURRENT carrier columns at the top
      ! of each step, so the water layer shields itself inside the
      ! relaxation rather than only between Picard passes.
      real(dp), dimension(:),   allocatable :: cph_klw
      real(dp), dimension(:,:), allocatable :: cph_jh2o, cph_joh
      ! CO photodissociation rate of the current columns [s^-1], rebuilt
      ! beside cph_klw so that the CO layer shields ITSELF inside the
      ! relaxation and not only between Picard passes.  It is the cell mean,
      ! as cph_klw is.
      real(dp), dimension(:),   allocatable :: cph_kco

      ! THE BASE FACE OF AN INFLOWING CARRIER (design decision D6, amended).
      !
      ! D6 gave the carrier partition a ZERO-FLUX base: the element totals
      ! are boundary data the handoff supplies, but the split of an element
      ! among its carriers is what the oxygen option exists to COMPUTE, so
      ! imposing it at the boundary would have made that option's own A/B
      ! gate circular.  That argument holds exactly when nothing upstream
      ! states the partition.
      !
      ! When a lower-atmosphere handoff DOES state it -- base.inp's
      ! q_H2_base, or the profile's value at the matching level -- the
      ! partition is not a quantity this operator computes but incoming
      ! data, and section 117 already pins it on the ghosts of the
      ! ionization sweep.  Leaving the base cell free while the ghost below
      ! it is pinned puts the same defect section 117 removed one cell
      ! higher up: the equation of state and the species state would again
      ! describe different gas across one face.  So where the handoff speaks,
      ! the advective inflow at the base face carries ITS composition.
      !
      ! The two branches cannot both apply.  Decision D7 refuses q_H2_base
      ! while the oxygen chemistry is on, so a run either has a handoff
      ! partition (H2 only -- nothing states OH, H2O or CO) or computes its
      ! own (zero flux, as D6 wrote it).
      !
      ! ONLY THE ADVECTIVE TERM.  Face 0 still carries no diffusive flux.
      ! The handoff supplies the composition of the gas that flows in, not a
      ! mixing rate across a boundary whose gradient is set by the ghost
      ! spacing and by a K_zz that describes an unresolved region.
      !
      ! It is carried entirely by the GHOST COMPOSITION now: the advective
      ! term is the divergence of the face species fluxes, whose face
      ! composition is reconstructed from the ghosts wherever the face mass
      ! flux flows inward (carrier_advective_divergence,
      ! carrier_face_mass_fraction).  The Dirichlet record the cell-velocity
      ! form needed -- which carrier is pinned, at what value, and whether
      ! the base face is inflowing -- is gone with that form: the direction
      ! is the sign of the face mass flux and the value is the ghost's.
      ! WHETHER THE ROWS OF THIS OPERATOR CARRY THE MATERIAL ADVECTIVE TERM.
      ! In the marching loop they do not: the advection of a carrier is the
      ! divergence of the hydrodynamic face mass flux and is taken with the
      ! mass row inside the Runge-Kutta stages, so a term here would be that
      ! transport a second time.  In the FIXED-WIND RELAXATION they do: that
      ! path takes many transport steps between two hydrodynamic stages and
      ! is asked for the state at which the whole carrier equation balances
      ! at a wind that does not move, so the advection has to be in the rows
      ! it solves.  It is the same division the element operator makes.
      logical, save :: carrier_rows_advect = .false.

      ! Diagnostics of the last step, read by write_output.
      integer,  save :: pct_newton_steps = 0
      real(dp), save :: pct_newton_resid = 0.0d0
      ! The same residual measured on the Newton's last trial, BEFORE the
      ! element limiter rescaled it.  Kept so that the ordering of the two
      ! can be measured rather than asserted: where the limiter did nothing
      ! the pair is bitwise equal, and where it acted the reported residual
      ! is the one of the state that actually leaves the solve.  Printed by
      ! the EXHALE_CARRIER_DEBUG line only.
      real(dp), save :: pct_newton_resid_before_limit = 0.0d0
      integer,  save :: pct_limited      = 0
      ! WHICH CELLS THE ELEMENT CONSTRAINTS OF THE LAST CALL MOVED.  A cell
      ! an element clamp touched no longer holds the state the Newton's last
      ! residual was measured at: the clamp changes the carrier densities,
      ! and through the free-atom closures of carrier_source (n_hi and n_o0,
      ! each of them the element total less the carriers) it changes every
      ! coupled row of that cell.  So a row residual left in ANY carrier of
      ! such a cell is the constraint's and not the solver's, and the
      ! question has to be asked of the CELL.
      !
      ! A CARRIER-INDEXED FLAG CANNOT ANSWER IT, in both directions.  The
      ! clamps act cell by cell while a carrier flag is a maximum over the
      ! grid, so a clamp active anywhere leaves that carrier's entry
      ! permanently true and excuses a solver failure in any cell of the
      ! domain; and a clamp on CO alone leaves OH and H2O unflagged although
      ! it moves their rows through n_o0.
      logical, save, allocatable :: pct_cell_constrained(:)
      ! The verdict on the state the last carrier solve returned.  Kept so
      ! the caller and the retry controller read one decision rather than
      ! each forming their own from the raw status.
      type(carrier_verdict), save :: pct_verdict
      ! How many times the retry controller halved the substep of the LAST
      ! interval it was asked for.  It is a property of that interval and
      ! not of any one attempt, so it is not part of the attempt
      ! checkpoint: the stop message and the retry test read it to say
      ! whether the controller went to its floor or stopped short of it
      ! because halving could not help.
      integer,  save :: pct_retry_halvings = 0
      ! ---- THE DOMAIN RECORD OF THE ONE-SIDED CO MODEL -----------------
      !
      ! The CO row destroys CO and never forms it, and what makes that
      ! legitimate in a cell is the ordering
      !
      !     tau_dest  <<  tau_res  <<  tau_form
      !
      ! (docs/b1_target_system_20260906.md T5.2).  The first inequality says
      ! the destruction is fast enough that the transported CO reaches the
      ! balance the rates set before it leaves the cell; the second says the
      ! omitted formation is too slow to rebuild what was destroyed while
      ! the gas is there.  Neither is an assumption this code may make
      ! silently, so both are evaluated on the state and reported.
      !
      !     tau_dest = 1/(k_D1 n(He+) + k_CO)
      !     tau_res  = min( r/|v| , dr^2/(D_CO + K_zz) )
      !     tau_form = 1/(k_8597 n_O)     [record only, never a row]
      !
      ! tau_res uses the code's own two timescale definitions, the same ones
      ! output/Oxygen_chemistry.txt already writes for H2, with the CO
      ! diffusion coefficient in place of the H2 one; T5.2's gradient-scale
      ! form needs a gradient scale this code does not form and is a later
      ! refinement.  A cell is IN DOMAIN when tau_dest <= f_dom tau_res.
      !
      ! WHAT THE RECORD DOES AND DOES NOT DO.  It is informational.  The
      ! rates are physics and are evaluated in every cell; nothing is
      ! switched off where the test fails.  What fails below the helium
      ! ionization front is the OMISSION OF FORMATION, not the destruction
      ! terms, and over the run's own physical time a destruction that slow
      ! cannot remove what an equally omitted formation cannot rebuild, so
      ! the transported value stands there.  The record says where that is
      ! the case, and how far.
      !
      ! Two further out-of-domain states are counted beside it.  Cells above
      ! the 512 K excitation-temperature limit of the shielding table
      ! (co_self_shielding_table states why the table stops there and what
      ! it costs); and cells in which the He+ charge transfer removes He+
      ! faster than it is otherwise lost, where the frozen-He+ approximation
      ! of the carrier step over-states channel D1.
      !
      ! Indexed by the ledger family exactly as the counters above are.
      !
      ! THEY ARE NOT IN THE STEP CHECKPOINT, and that is the difference
      ! between a record of a CLAMP and a record of a STATE.  The ceiling
      ! counters above roll back with a refused trial, because a clamp that
      ! was applied inside a trial the run then discarded never acted on the
      ! run's history.  This record says whether the CO model was in domain
      ! in the cells of a state, and that question has an answer whether or
      ! not the step that read it was accepted -- indeed a relaxation that
      ! never covers an interval is exactly the case in which a reader needs
      ! it.  Rolling it back would erase it in that case.
      integer,  save :: co_dom_cells_out(2)   = 0
      integer,  save :: co_dom_cells_hot(2)   = 0
      integer,  save :: co_dom_cells_hep(2)   = 0
      real(dp), save :: co_dom_worst_ratio(2) = 0.0d0
      real(dp), save :: co_dom_worst_r(2)     = 0.0d0
      real(dp), save :: co_dom_form_ratio(2)  = 0.0d0
      ! The stated constant of the domain test.  One decade is the weakest
      ! reading of "much less than" that is still a statement, and the
      ! measured profile of the oxygen_chemistry case
      ! (docs/co_destruction_rates_literature_20260906.md sec. 12.5) puts
      ! every cell above r = 1.05 at tau_dest/tau_res of 1.8e-1 to 9.6e-5
      ! and every cell below r = 1.04 at 1.5 to 6.7e2, so any threshold
      ! between 1e-1 and 1 selects the same boundary on that state.  It is a
      ! parameter of the model, not a knob to tune a result with.
      real(dp), parameter :: co_domain_f_dom = 0.1d0

      integer,  save :: pct_worst_j       = 0
      integer,  save :: pct_worst_ic      = 1
      ! Where the last fixed-wind relaxation moved the composition most.
      integer,  save :: pct_drift_j     = 0
      integer,  save :: pct_drift_ic    = 1
      real(dp), save :: pct_worst_limit  = 0.0d0
      ! Molecular diffusion coefficient of each carrier on the grid, kept
      ! from the last step so the run can print the transport time scale
      ! beside the chemical one it already prints.
      real(dp), dimension(:,:), allocatable :: pct_Dco
      ! The sum of each carrier row's own term magnitudes, from the last
      ! residual assembly, THE TIME TERM INCLUDED.  This is the scale of the
      ! equation the Newton iterates on (sec. 139), and the solver's row and
      ! column scales are built from it, so there is one definition of it.
      real(dp), dimension(:,:), allocatable :: row_terms
      ! THE SAME ROW WITHOUT ITS TIME TERM, and this is the scale a state is
      ! JUDGED on.  A carrier row is
      !
      !     (n - n_old)/dt = transport + reaction,
      !
      ! and row_terms above sums the magnitudes of all three, the time term
      ! included, because that is the size of the equation the Newton is
      ! iterating on at the step it was given.  It is the wrong scale for a
      ! statement ABOUT THE RETURNED STATE.  The time term carries 1/dt, so
      ! it is not a property of the state at all: the integrator sets it and
      ! can make it arbitrarily large by shortening the step.  An imbalance
      ! the state really carries -- a flux left unbalanced by the element
      ! clamp of a neighbouring cell, say -- is a TRANSPORT term and does not
      ! scale with dt, so |res|/row_terms of one fixed physical defect falls
      ! like dt and passes any fixed relative floor once the step is short
      ! enough.  MEASURED on oxygen_chemistry: the refused row of the CO
      ! balance at cell 299, 1.04e-7 of its terms at the full step, read
      ! 7.9e-9 after ten halvings and passed a 1e-8 floor with nothing about
      ! the clamp behind it changed.
      !
      ! row_terms_phys is therefore |transport| + |reaction| with the same
      ! absolute floor, and the floor is converted to volumetric-rate units
      ! by the cell's signal-crossing rate ALONE and not by max(1/dt, rate):
      ! a floor that grows as the step shrinks reopens the same door in the
      ! cells where the floor is what the row is measured against.  This is
      ! the scale the returned-state acceptance and the stationary
      ! certification read; the Newton keeps row_terms for its own iterate.
      real(dp), dimension(:,:), allocatable :: row_terms_phys
      ! Scale of the carrier ROW, and separately of its UNKNOWN, from that
      ! same assembly.  The ROW scale is row_terms over one code time unit,
      ! so that the residual divided by it is the dimensionless imbalance
      ! the acceptance test uses.  The COLUMN scale is the unknown's own
      ! magnitude with a floor that does not collapse where it does -- the
      ! same construction as the momentum unknown's |rho v| + rho c_s.
      !
      ! THE TWO MUST DIFFER, and that was measured.  n(H2) alone cannot be
      ! either: it falls four decades across the front and is effectively
      ! zero beyond it, which is the pathology carrier_row_scale was written
      ! to avoid.  Using the ROW scale for both makes the row's merit
      ! contribution the CHEMICAL rate, about 1e3 per code time at the front
      ! against the hydrodynamic rows' 1e-6, so the line search sees nothing
      ! but the carrier and the damped Gauss-Newton escape lands on a state
      ! whose chemistry cannot be certified (merit 89.8, abort at iteration
      ! 2).  Using the COLUMN scale for both makes the Newton propose
      ! carrier steps 84 times the cell's own n(H2) and cuts the line search
      ! to lam = 4.9e-4.
      !
      ! ONE PAIR FOR EACH SOLVED CARRIER, not for H2 alone: the stationary
      ! system carries a row and an unknown for every transported balance
      ! the configuration activates, and a carrier scaled by another
      ! carrier's terms is a row the Newton cannot weigh.
      real(dp), dimension(:,:), allocatable :: col_scale_car, row_scale_car
      ! Element headroom of each cell from that same assembly [cm^-3]: the
      ! free density of the carrier's OWN element divided by the nuclei of
      ! that element the carrier holds (carrier_element_headroom).
      real(dp), dimension(:,:), allocatable :: headroom_car
      logical :: headroom_set = .false.

      ! THE MODULE ARRAYS THE CARRIER RESIDUAL OVERWRITES, HELD ASIDE AND PUT
      ! BACK, so that measuring the carrier balance of a state does not change
      ! what the next transport step or the next output reads.
      !
      ! carrier_steady_residual declares its arguments intent(in) but is not a
      ! function of them alone: it refreshes the frozen background densities,
      ! the photolysis rates, the advection correction, the row terms, the two
      ! H2 scales and (under weno_mode = 1) the element headroom, all of which
      ! are module state that later calls read. A certification of a state has
      ! to be able to evaluate that balance without being an operation ON the
      ! state, so the pair below stores exactly the arrays the call writes and
      ! reinstates them, allocation status included: an array that was not
      ! allocated before the call is deallocated again, so the module comes
      ! back to the state it was in and not merely to the same numbers.
      type, public :: carrier_module_state
         real(dp), allocatable :: cbg_nhii(:), cbg_nh2p(:), cbg_nh3p(:)
         real(dp), allocatable :: cbg_nhehp(:), cbg_nhei(:), cbg_nheii(:)
         real(dp), allocatable :: cbg_nheiii(:), cbg_nheiTR(:), cbg_ne(:)
         real(dp), allocatable :: cbg_nOion(:), cbg_nCion(:)
         real(dp), allocatable :: cph_klw(:), cph_kco(:)
         real(dp), allocatable :: cph_jh2o(:,:), cph_joh(:,:)
         real(dp), allocatable :: row_terms(:,:)
         real(dp), allocatable :: row_terms_phys(:,:)
         real(dp), allocatable :: col_scale_car(:,:), row_scale_car(:,:)
         real(dp), allocatable :: headroom_car(:,:)
         logical  :: headroom_set = .false.
         integer  :: pct_worst_j  = 0
         integer  :: pct_worst_ic = 1
      end type carrier_module_state

      ! HOW MANY TIMES ONE TRANSPORT INTERVAL MAY HALVE ITS CARRIER SUBSTEP,
      ! and therefore the smallest substep the controller will attempt,
      ! dt/2**carrier_retry_max.  The substep only ever halves inside an
      ! interval and never grows back, so the bound on the retries and the
      ! floor on the substep are one rule and not two: the rejection that
      ! arrives when the substep is already at the floor is the exhaustion,
      ! and no interval can cost more than 2**carrier_retry_max substeps.
      integer,  parameter, public :: carrier_retry_max = 8
      real(dp), parameter :: carrier_substep_floor =                     &
                             1.0d0/2.0d0**carrier_retry_max

      ! TEST HOOKS.  Both default to the production behavior, no production
      ! path writes either, and src/tests/carrier_retry/ is the only reader
      ! of what they cause.  They exist because the branches they reach --
      ! restore, halve, exhaust -- cannot otherwise be driven on a column
      ! that solves, and a column that does not solve would test the
      ! chemistry and not the controller.
      !   (a) refuse the leading attempts of every interval, whatever the
      !       returned state satisfies;
      integer,  save :: carrier_reject_leading_attempts_for_test = 0
      !   (b) start the interval at this fraction of dt instead of at 1, so
      !       that a subdivided integration OVER THE SAME FROZEN BACKGROUND
      !       is available as the reference a retried result is compared
      !       with.  Calling the transport step twice with half the step is
      !       not that reference: it rebuilds the background from the
      !       half-advanced composition.
      real(dp), save :: carrier_initial_substep_for_test = 1.0d0

      ! ATTEMPT STATISTICS of the controller, never restored: they count
      ! attempts, which is what happened, and not contributions, which is
      ! what the accepted state holds.
      integer, save :: carrier_substeps_attempted = 0
      integer, save :: carrier_substeps_rejected  = 0
      integer, save :: carrier_substeps_accepted  = 0
      integer, save :: carrier_intervals_retried  = 0

      ! THE INTERVALS THAT WERE NOT COVERED, and what the run did instead.
      ! An exhausted interval is an attempt and belongs to the attempts
      ! ledger, so it is never restored and it carries no family index.
      ! The first and the last marching step at which one occurred bound
      ! the part of the history that rests on frozen carriers, and the
      ! worst row is the largest imbalance, relative to that row's own
      ! physical terms, any of those refusals rested on.
      integer,  save :: carrier_intervals_exhausted = 0
      integer,  save :: carrier_exhausted_step_first = -1
      integer,  save :: carrier_exhausted_step_last  = -1
      integer,  save :: carrier_exhausted_worst_j    = 0
      integer,  save :: carrier_exhausted_worst_ic   = 1
      real(dp), save :: carrier_exhausted_worst_ratio = 0.0d0
      real(dp), save :: carrier_exhausted_worst_phys  = 0.0d0
      real(dp), save :: carrier_exhausted_worst_full  = 0.0d0

      ! WHETHER THE CARRIER SUBSYSTEM OF THIS RUN CAN BE CERTIFIED FROM ITS
      ! HISTORY.  False once an interval was left uncovered and the march
      ! went on from the entry carriers: those carriers solve no equation
      ! over that step, and no later step makes that step's state a
      ! solution.  It says nothing about the FINAL state, which the
      ! stationary certification measures on its own rows; the two answers
      ! are different questions and the certification reports both.
      logical, save :: carrier_history_certifiable_flag = .true.

      ! Whether an uncovered interval ends the run where the run mode says
      ! it must (physical integration).  Production is .true.; a test
      ! driver that has to read the status of the refusal writes .false.,
      ! as the energy update and the conduction update are read.
      logical, save :: carrier_transport_stop_on_failure = .true.
      ! How many times the run WOULD have ended in this operator and did
      ! not, because a test driver had suppressed the stop.  Zero in every
      ! production run, and the only way a test can state that the stop was
      ! reached.
      integer, save :: carrier_transport_stops_suppressed = 0

      ! WHAT THE TRANSPORT STEP DID WITH THE INTERVAL IT WAS ASKED FOR.
      integer, parameter :: carrier_interval_covered   = 0
      integer, parameter :: carrier_interval_exhausted = 1
      integer, save      :: carrier_interval_status_last =                &
                            carrier_interval_covered

      ! THE VERDICT OF THE STEP IN PROGRESS.  Reset at the entry of
      ! photochemical_transport_step and written only where an interval was
      ! actually integrated, so it never carries an earlier step's answer.
      type(carrier_verdict), save :: carrier_step_verdict_val =           &
           carrier_verdict(accepted = .false., reason = carrier_no_interval)

      ! ROUND-OFF LIMITED ROWS OVER THE RUN.  A row whose residual is the
      ! round-off of its own full terms is accepted (carrier_verdict), and
      ! how often that happened is informational: it says in how many
      ! accepted substeps double precision, and not the physics, set the
      ! accuracy of a row.  Zero over a run IS the statement that every
      ! accepted row was certified on its physical terms.
      integer, save :: carrier_roundoff_rows_total = 0
      integer, save :: carrier_roundoff_substeps   = 0
      integer, save :: carrier_roundoff_rows_last  = 0

      ! THE ROWS A REFUSAL RESTS ON, with the two scales each of them was
      ! measured against: the physical terms the acceptance uses and the
      ! full terms whose round-off decides whether the row could have been
      ! solved any further.  Kept here rather than in the verdict because
      ! the checkpoint restore that follows a refusal puts row_terms and
      ! row_terms_phys back to the entry assembly, so the stop could not
      ! read them afterwards.  Bounded: a stop needs the rows, not all of
      ! them.
      integer, parameter, public :: carrier_refused_rows_max = 8
      integer,  save :: n_refused_reported = 0
      integer,  save :: refused_j(carrier_refused_rows_max)  = 0
      integer,  save :: refused_ic(carrier_refused_rows_max) = 1
      real(dp), save :: refused_res(carrier_refused_rows_max)  = 0.0d0
      real(dp), save :: refused_phys(carrier_refused_rows_max) = 0.0d0
      real(dp), save :: refused_full(carrier_refused_rows_max) = 0.0d0

      ! THE CARRIER-LOCAL CHECKPOINT: everything one transport attempt
      ! writes and nothing else.
      !
      ! LABELED PARTIAL RECOVERY.  This restores the carrier subsystem --
      ! the transported fractions, the frozen background arrays, the
      ! photolysis rates, the advection correction, the row scales, the
      ! headroom, the constraint diagnostics and the last verdict.  It does
      ! NOT restore the hydrodynamic state,
      ! the ionization state or the energy of the step around it, and it
      ! makes no claim to: the enumerated checkpoint of the whole attempted
      ! step is B3a (b1_target_system_20260906_draft sec. 7, contract T7.3).
      ! It is sound on its own because the transport step reads rho, v and
      ! the background as data and writes only f_sp, and only after the
      ! interval is complete.
      !
      ! The carrier fractions are optional in the pair below because at the
      ! top of a step they do not exist yet: the entry state of the
      ! carriers is f_sp, which the step has not written.
      type, public :: carrier_checkpoint
         real(dp), allocatable :: fc(:,:)
         ! The arrays a residual assembly overwrites, in the enumeration
         ! save_carrier_module_state already keeps (A2), reused rather than
         ! written a second time so the two cannot drift apart.
         type(carrier_module_state) :: bg
         real(dp), allocatable :: pct_Dco(:,:)
         logical,  allocatable :: pct_cell_constrained(:)
         integer  :: pct_newton_steps = 0
         real(dp) :: pct_newton_resid = 0.0d0
         real(dp) :: pct_newton_resid_before_limit = 0.0d0
         integer  :: pct_limited     = 0
         real(dp) :: pct_worst_limit = 0.0d0
         integer  :: pct_drift_j     = 0
         integer  :: pct_drift_ic    = 1
         type(carrier_verdict) :: pct_verdict
      end type carrier_checkpoint


      contains

      ! Fix the transported set for this run.  Called once from input_read,
      ! after every key is parsed, so that no routine can see a half-decided
      ! carrier set.
      subroutine carrier_set_init()
      carrier_solved = .false.
      carrier_solved(ic_H2) = .true.
      if (thereis_oxychem) then
         carrier_solved(ic_OH)  = .true.
         carrier_solved(ic_H2O) = .true.
         carrier_solved(ic_CO)  = .true.
      endif
      if (ionization_transport) carrier_solved(ic_Hp) = .true.
      ! n_carrier is the last index in the solved set: the loops run
      ! 1..n_carrier and skip the ones carrier_solved leaves out.
      if (ionization_transport) then
         n_carrier = ic_Hp              ! ... H2 [, OH, H2O, CO] and H+
      else if (thereis_oxychem) then
         n_carrier = ic_CO              ! H2, OH, H2O, CO
      else
         n_carrier = ic_H2              ! H2 alone
      endif

      ! THE TRANSPORTED SET.  Each solved carrier is advected on the
      ! hydrodynamic face mass fluxes inside the Runge-Kutta stages, on the
      ! same faces, areas, volumes and time step as the density, so that a
      ! cell can only lose the carrier the mass row says it loses.  The rows
      ! this module solves then carry the diffusive and drift fluxes and the
      ! chemistry alone.
      !
      ! THE BASE INFLOW COMPOSITION travels with the declaration.  Where a
      ! handoff states a carrier's base partition -- base.inp's q_H2_base, or
      ! the profile's value at the matching level -- the inner ghosts hold it
      ! and the inflowing base face carries it in.  Where nothing states it
      ! the ghosts take the base cell's own partition, which is the
      ! zero-gradient condition of design decision D6 written on the face the
      ! gas crosses.  The proton is never stated: no key gives an ionization
      ! fraction below the base, so it keeps the zero-gradient condition in
      ! every configuration.
      call advected_carrier_reset()
      if (thereis_mol .and. carrier_transport) then
         if (carrier_solved(ic_H2))                                      &
            call advected_carrier_register(isp_H2,                       &
                                           base_h2_composition_imposed())
         if (carrier_solved(ic_OH))                                      &
            call advected_carrier_register(isp_OH,  .false.)
         if (carrier_solved(ic_H2O))                                     &
            call advected_carrier_register(isp_H2O, .false.)
         if (carrier_solved(ic_CO))                                      &
            call advected_carrier_register(isp_CO,  .false.)
         if (carrier_solved(ic_Hp))                                      &
            call advected_carrier_register(isp_HII, .false.)
      endif
      end subroutine carrier_set_init

      ! ------------------------------------------------------------- !

      ! Diagnostics of the last transport step: Newton iterations, the
      ! relative residual it stopped at, how many cells needed the element
      ! limiter, and the worst relative overshoot it removed.
      subroutine carrier_transport_diagnostics(nsteps, resid, nlim, worst)
      integer,  intent(out) :: nsteps, nlim
      real(dp), intent(out) :: resid, worst
      nsteps = pct_newton_steps
      resid  = pct_newton_resid
      nlim   = pct_limited
      worst  = pct_worst_limit
      end subroutine carrier_transport_diagnostics

      ! ------------------------------------------------------------- !

      ! How far the retry controller halved the substep of the last
      ! interval: 0 when the first attempt settled it, carrier_retry_max
      ! when it went to its floor, and anything between when it stopped
      ! because a further halving could not change the outcome.
      integer function carrier_last_interval_halvings() result(nh)
      nh = pct_retry_halvings
      end function carrier_last_interval_halvings

      ! ------------------------------------------------------------- !

      ! Did an element constraint of the last limiter call move cell j?
      ! Out of range, or before any solve has run, the answer is no: there
      ! is no constraint to attribute a residual to.
      logical function carrier_cell_is_constrained(j) result(bound)
      integer, intent(in) :: j
      bound = .false.
      if (.not. allocated(pct_cell_constrained)) return
      if (j .lt. 1 .or. j .gt. size(pct_cell_constrained)) return
      bound = pct_cell_constrained(j)
      end function carrier_cell_is_constrained

      ! ------------------------------------------------------------- !

      ! HOW A STALLING NEWTON ITERATION STOPPED, as a statement about the
      ! residual it is stopping at rather than about the caller.
      !
      ! An iteration is stalling when it has taken at least three steps and
      ! the last one did not halve the residual.  What the stall then means
      ! depends on where the residual stands:
      !
      !  - at or below newton_floor: the absolute floor of this
      !    discretization, an accuracy statement;
      !  - at or below newton_drop times the residual the iteration STARTED
      !    from: not an accuracy statement at all.  The predicate admits
      !    rstart = 1e8, rprev = 1e2, rnorm = 60, which satisfies both
      !    "did not halve" and "fell six decades", with an absolute
      !    residual of 60 (review R6).  It is returned as its own class so
      !    that it can be reported and never confused with the first;
      !  - above both: the iteration is simply stuck.
      integer function carrier_stall_class(it, rnorm, rprev, rstart)      &
                       result(cls)
      integer,  intent(in) :: it
      real(dp), intent(in) :: rnorm, rprev, rstart
      cls = carrier_stall_none
      if (it .lt. 3) return
      if (rnorm .le. 0.5d0*rprev) return
      if (rnorm .le. newton_floor) then
         cls = carrier_stall_at_floor
      else if (rnorm .le. newton_drop*rstart) then
         cls = carrier_stall_at_relative_drop
      else
         cls = carrier_stall_open
      endif
      end function carrier_stall_class

      ! ------------------------------------------------------------- !

      ! IS THE STATE THAT IS ABOUT TO BE WRITTEN BACK A SOLUTION?
      !
      ! The question is asked of the RETURNED state and of nothing else.
      ! Whether the Newton reported convergence, ran out of iterations or
      ! lost its line search is a fact about the iteration, not about the
      ! state: the element limiter runs after the iteration ends and moves
      ! whole carrier families in the cells it binds, and an iteration that
      ! stalled at a relative drop can return a residual of any size at all.
      ! So the raw status is carried in the verdict for diagnosis and takes
      ! no part in the decision.
      !
      ! EVERY UNCONSTRAINED ROW HAS TO PASS, not the worst one.  Taking the
      ! largest returned residual and asking only whether ITS cell was
      ! constrained lets a constraint acting anywhere excuse a failed row
      ! everywhere else in the column (the 3.3 defect of the review): a
      ! constrained cell says something about that cell and nothing about
      ! any other.
      !
      ! WHAT A CONSTRAINT DOES EXCUSE is the cell it moved, and the whole of
      ! it: the clamps rescale carriers cell by cell, and through the free
      ! atomic closures of carrier_source that moves every coupled row of
      ! that cell, so the attribution is per cell (measured in
      ! src/tests/carrier_constraint_attribution/).  Those rows are counted
      ! and reported as the constraint's, never certified.
      !
      ! A NON-FINITE ROW IS REFUSED WHEREVER IT IS, constrained or not: a
      ! constraint cannot make a NaN a solution, and the maximum the solve
      ! reports would not see one (a comparison with a NaN is false, so it
      ! never becomes the running maximum).
      subroutine carrier_returned_state_verdict(res, terms, constrained,  &
                                                raw_status, stall,        &
                                                verdict, terms_full,      &
                                                row_roundoff_limited)
      ! Row residuals and the matching row scales of the returned state,
      ! (cell, carrier).  The scale is the PHYSICAL one, row_terms_phys,
      ! which is dt-independent; the decision is taken on the arguments
      ! alone, so a test states it on constructed rows.
      real(dp), dimension(:,:), intent(in)  :: res, terms
      ! Which cells an element constraint moved in the call that produced
      ! this state.
      logical,  dimension(:),   intent(in)  :: constrained
      integer,                  intent(in)  :: raw_status, stall
      type(carrier_verdict),    intent(out) :: verdict
      ! The FULL terms of the same rows, the time term included.  With it
      ! the verdict can say of a row above the floor whether its residual
      ! is already the round-off of the arithmetic that formed the row.
      ! Such a row is at the limit of double precision at this substep and
      ! is ACCEPTED, counted as round-off limited: it is not evidence of a
      ! defect in the returned state, and no iteration and no shorter
      ! substep can reduce it (carrier_row_roundoff).  Absent, no row is
      ! classified and every row above the floor refuses, which is what a
      ! caller that builds rows by hand gets.
      real(dp), dimension(:,:), optional, intent(in) :: terms_full
      ! Which rows were classified round-off limited, row by row, so that
      ! the caller can report the two scales of the rows that were not.
      logical, dimension(:,:), optional, intent(out) :: row_roundoff_limited

      integer  :: j, ic, ncell, ncar
      logical  :: has_full, at_roundoff
      real(dp) :: rel

      verdict%accepted            = .true.
      verdict%reason              = carrier_accept_ok
      verdict%jworst              = 0
      verdict%icworst             = 1
      verdict%rworst              = 0.0d0
      verdict%n_cells_constrained = 0
      verdict%n_rows_constrained  = 0
      verdict%n_rows_above        = 0
      verdict%n_rows_above_beside_constrained = 0
      verdict%n_rows_roundoff_limited = 0
      verdict%raw_status          = raw_status
      verdict%stall               = stall
      if (present(row_roundoff_limited)) row_roundoff_limited = .false.

      has_full = .false.
      ncell = min(size(res,1), size(terms,1), size(constrained))
      ncar  = min(n_carrier, size(res,2), size(terms,2))
      if (present(terms_full)) then
         has_full = .true.
         ncell    = min(ncell, size(terms_full,1))
         ncar     = min(ncar,  size(terms_full,2))
      endif

      do j = 1, ncell
         if (constrained(j)) verdict%n_cells_constrained =                &
                             verdict%n_cells_constrained + 1
         do ic = 1, ncar
            ! A carrier this run does not solve carries no equation, so it
            ! has no row to accept or refuse.
            if (.not. carrier_solved(ic)) cycle
            if (.not. finite_real(res(j,ic)) .or.                         &
                .not. finite_real(terms(j,ic))) then
               verdict%accepted = .false.
               verdict%reason   = carrier_reject_nonfinite
               verdict%jworst   = j
               verdict%icworst  = ic
               verdict%rworst   = res(j,ic)
               return
            endif
            if (constrained(j)) then
               verdict%n_rows_constrained = verdict%n_rows_constrained + 1
               cycle
            endif
            rel = abs(res(j,ic))/max(terms(j,ic), 1.0d-300)
            ! A ROW AT THE ROUND-OFF OF ITS OWN FULL TERMS CARRIES NO
            ! INFORMATION ABOUT THE STATE.  Its residual is what forming
            ! the row costs in double precision, so it neither refuses the
            ! state nor decides the worst row: it is counted, flagged and
            ! left out of the decision.  The acceptance on the physical
            ! terms is unchanged for every row that is not there, and that
            ! acceptance stays independent of the substep.
            at_roundoff = .false.
            if (has_full .and. rel .gt. carrier_accept_tol)               &
               at_roundoff = (abs(res(j,ic)) .le. carrier_row_roundoff    &
                                                  *terms_full(j,ic))
            if (at_roundoff) then
               verdict%n_rows_roundoff_limited =                          &
               verdict%n_rows_roundoff_limited + 1
               if (present(row_roundoff_limited))                         &
                  row_roundoff_limited(j,ic) = .true.
               cycle
            endif
            if (rel .gt. verdict%rworst) then
               verdict%rworst  = rel
               verdict%jworst  = j
               verdict%icworst = ic
            endif
            if (rel .gt. carrier_accept_tol) then
               verdict%n_rows_above = verdict%n_rows_above + 1
               if (constrained(max(j-1,1)) .or.                           &
                   constrained(min(j+1,ncell)))                           &
                  verdict%n_rows_above_beside_constrained =               &
                  verdict%n_rows_above_beside_constrained + 1
            endif
         enddo
      enddo

      ! rworst is the worst row that CAN say something about the state, so
      ! a state whose only rows above the floor are round-off limited is
      ! accepted and carries the count of them.
      if (verdict%rworst .gt. carrier_accept_tol) then
         verdict%accepted = .false.
         verdict%reason   = carrier_reject_row_residual
      endif
      end subroutine carrier_returned_state_verdict

      ! ------------------------------------------------------------- !

      ! The verdict on the state the last carrier solve of the RUN
      ! returned.  It belongs to whichever interval produced it and to no
      ! particular step, so a step in which the operator had nothing to do
      ! leaves it standing: use carrier_step_verdict to ask about a step.
      type(carrier_verdict) function carrier_last_verdict_of_run() result(v)
      v = pct_verdict
      end function carrier_last_verdict_of_run

      ! The verdict of THIS transport step, reset at the entry of
      ! photochemical_transport_step.  A step that integrated no interval
      ! carries reason carrier_no_interval and accepted = .false.
      type(carrier_verdict) function carrier_step_verdict() result(v)
      v = carrier_step_verdict_val
      end function carrier_step_verdict

      ! Why a state was refused, in one line, for the log and the stop.
      function carrier_verdict_reason_text(v) result(txt)
      type(carrier_verdict), intent(in) :: v
      character(len=72) :: txt
      select case (v%reason)
      case (carrier_accept_ok)
         txt = 'every unconstrained row within the returned-state floor'
      case (carrier_reject_nonfinite)
         txt = 'a carrier row of the returned state is not a finite number'
      case (carrier_reject_forced_for_test)
         txt = 'the test hook refused this attempt; no physical row failed'
      case (carrier_no_interval)
         txt = 'no carrier interval was integrated in this step'
      case default
         txt = 'an unconstrained carrier row of the returned state is'//   &
               ' above the floor'
      end select
      end function carrier_verdict_reason_text

      ! How the iteration itself ended, for the same line.  It is a fact
      ! about the iteration and never an acceptance.
      function carrier_solve_status_text(status) result(txt)
      integer, intent(in) :: status
      character(len=72) :: txt
      select case (status)
      case (carrier_solve_converged)
         txt = 'the residual reached the absolute floor'
      case (carrier_solve_line_search_failed)
         txt = 'the line search exhausted its halvings without reducing'//&
               ' the residual'
      case (carrier_solve_stagnated)
         txt = 'it stagnated at a relative drop, above the absolute floor'
      case (carrier_solve_unreachable_at_substep)
         txt = 'it reached the round-off of the row terms, still above'// &
               ' the acceptance'
      case default
         txt = 'the iteration cap was reached'
      end select
      end function carrier_solve_status_text

      ! ------------------------------------------------------------- !

      ! THE DOMAIN RECORD OF THE ONE-SIDED CO MODEL, taken on the state an
      ! accepted carrier step hands on.  It reads and never writes the
      ! physics: nothing here changes a density or a rate.
      !
      ! It runs over the INTERIOR cells only, because a ghost carries no
      ! volume of the domain, and it marks the ledger family of the moment,
      ! so a cell first out of domain during a numerical relaxation does not
      ! appear in the physical history's record.
      !
      ! The three counters and the two worst-case fields are cumulative:
      ! what a reader needs is whether the model was ever out of domain
      ! during the run, and over how much of the integration, not whether
      ! it was at the last step.
      subroutine carrier_co_domain_take(fc, nrho, TK, v, rp, nO_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fc
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: nrho, TK, v, rp
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: nO_free
      real(dp) :: n_co, k_dest, tau_dest, tau_adv, tau_dif, tau_res
      real(dp) :: ratio, dr_cm, d_co, n_o0, k_form, tau_form, hep_loss
      real(dp) :: k_hep_co
      integer  :: j, ifam

      if (.not. thereis_oxychem) return
      ifam = ledger_family
      k_hep_co = rk_D1_Hep_CO()
      do j = 1, N
         n_co = max(fc(j,ic_CO)*nrho(j), 0.0d0)
         if (n_co .le. 0.0d0) cycle

         ! tau_dest: the two channels of the row, at the same rates the row
         ! used.
         k_dest = k_hep_co*cbg_nheii(j) + cph_kco(j)
         if (k_dest .le. 0.0d0) then
            tau_dest = huge(1.0d0)
         else
            tau_dest = 1.0d0/k_dest
         endif

         ! tau_res: the code's own advective and diffusive cell times, with
         ! the CO diffusion coefficient.
         if (abs(v(j)) .gt. 0.0d0) then
            tau_adv = rp(j)/abs(v(j)*v0)
         else
            tau_adv = huge(1.0d0)
         endif
         dr_cm   = dr_j(j)*R0
         d_co    = carrier_diffusion_coefficient(j, ic_CO)
         tau_dif = dr_cm*dr_cm/max(d_co + kzz_cell(j), 1.0d-30)
         tau_res = min(tau_adv, tau_dif)

         ratio = tau_dest/max(tau_res, 1.0d-300)
         if (ratio .gt. co_domain_f_dom) then
            co_dom_cells_out(ifam) = co_dom_cells_out(ifam) + 1
            if (ratio .gt. co_dom_worst_ratio(ifam)) then
               co_dom_worst_ratio(ifam) = ratio
               co_dom_worst_r(ifam)     = r(j)
            endif
         endif

         ! The second inequality of T5.2, measured rather than asserted:
         ! the radiative association C + O -> CO + photon of RATE22 entry
         ! 8597 on the free atomic oxygen of the cell.  It enters no row.
         ! Free atomic oxygen: the SAME closure carrier_source uses, so the
         ! two cannot state different densities for one cell.
         n_o0   = max(nO_free(j) - cbg_nOion(j)                           &
                      - fc(j,ic_OH)*nrho(j) - fc(j,ic_H2O)*nrho(j)        &
                      - n_co, 0.0d0)
         k_form = rk_CO_radiative_association(max(TK(j), 1.0d0))*n_o0
         if (k_form .gt. 0.0d0) then
            tau_form = 1.0d0/k_form
            co_dom_form_ratio(ifam) = max(co_dom_form_ratio(ifam),        &
                                          tau_res/tau_form)
         endif

         ! Above the excitation-temperature limit of the shielding table the
         ! table is still evaluated -- there is nothing better -- and the
         ! cell is counted.
         if (TK(j) .gt. co_shield_tex_limit_K())                          &
            co_dom_cells_hot(ifam) = co_dom_cells_hot(ifam) + 1

         ! Where D1 removes He+ faster than recombination does it is the
         ! leading He+ sink of the cell, and the He+ this step reads is one
         ! sweep behind the CO it is reading; the sweep's own He+ row
         ! carries the reaction (System_HeH_mol row 2), so what this counts
         ! is the size of the operator split's lag, not a missing sink.
         hep_loss = bg_cell(j)%rcheiiB*cbg_ne(j)
         if (k_hep_co*n_co .gt. hep_loss)                                 &
            co_dom_cells_hep(ifam) = co_dom_cells_hep(ifam) + 1
      enddo
      end subroutine carrier_co_domain_take

      ! The record, read by the setup report and by the output headers.
      ! family absent = the run total over both ledgers.
      subroutine carrier_co_domain_record(ncells_out, ncells_hot,         &
                                          ncells_hep, worst_ratio,        &
                                          worst_r, form_ratio, family)
      integer,  intent(out) :: ncells_out, ncells_hot, ncells_hep
      real(dp), intent(out) :: worst_ratio, worst_r, form_ratio
      integer, intent(in), optional :: family
      if (present(family)) then
         ncells_out  = co_dom_cells_out(family)
         ncells_hot  = co_dom_cells_hot(family)
         ncells_hep  = co_dom_cells_hep(family)
         worst_ratio = co_dom_worst_ratio(family)
         worst_r     = co_dom_worst_r(family)
         form_ratio  = co_dom_form_ratio(family)
         return
      endif
      ncells_out  = sum(co_dom_cells_out)
      ncells_hot  = sum(co_dom_cells_hot)
      ncells_hep  = sum(co_dom_cells_hep)
      worst_ratio = maxval(co_dom_worst_ratio)
      worst_r     = co_dom_worst_r(maxloc(co_dom_worst_ratio, 1))
      form_ratio  = maxval(co_dom_form_ratio)
      end subroutine carrier_co_domain_record

      ! The stated constant of the domain test, so that the reports print
      ! the value in force rather than a copy of it.
      double precision function carrier_co_domain_f_dom() result(f)
      f = co_domain_f_dom
      end function carrier_co_domain_f_dom

      ! ------------------------------------------------------------- !

      ! TEST ONLY: a distinct amount into each ledger family of the domain
      ! record, so that the composition rules and the non-restoration can be
      ! measured on the production arrays.  The counts compose as a sum over
      ! the families and the worst ratio as a maximum, which are two
      ! different rules on one record; and the record is not part of the
      ! carrier checkpoint, so a restore must leave what this wrote in
      ! place.  No production path calls this.
      subroutine carrier_co_domain_perturb_for_test(n1, n2)
      integer, intent(in) :: n1, n2
      co_dom_cells_out   = co_dom_cells_out + [n1, n2]
      co_dom_cells_hot   = co_dom_cells_hot + [n1, n2]
      co_dom_cells_hep   = co_dom_cells_hep + [n1, n2]
      co_dom_worst_ratio = co_dom_worst_ratio + [dble(n1), dble(n2)]
      co_dom_worst_r     = [dble(n1), dble(n2)]
      co_dom_form_ratio  = co_dom_form_ratio + [dble(n1), dble(n2)]
      end subroutine carrier_co_domain_perturb_for_test

      ! The controller's attempt statistics over the run: carrier substeps
      ! attempted, refused, accepted, and how many intervals needed at least
      ! one retry.
      subroutine carrier_attempt_record(nattempt, nreject, naccept,       &
                                        nretried)
      integer, intent(out) :: nattempt, nreject, naccept, nretried
      nattempt = carrier_substeps_attempted
      nreject  = carrier_substeps_rejected
      naccept  = carrier_substeps_accepted
      nretried = carrier_intervals_retried
      end subroutine carrier_attempt_record

      ! WHAT THE LAST TRANSPORT STEP DID WITH ITS INTERVAL: covered, or
      ! exhausted and not covered.  The same value the step returns, kept
      ! so a caller that did not take it can still read it.
      integer function carrier_last_interval_status() result(st)
      st = carrier_interval_status_last
      end function carrier_last_interval_status

      ! THE INTERVALS THIS RUN LEFT UNCOVERED.  n is how many there were,
      ! step_first and step_last the first and the last marching step at
      ! which one occurred (-1 when there was none), and the last four the
      ! worst row any of those refusals rested on: its cell, its carrier,
      ! its imbalance relative to the row's own physical terms, and the two
      ! scales it was measured against.  A nonzero n is the statement that
      ! the carriers of that many steps are frozen entry values rather than
      ! a solution of the transport-chemistry system.
      subroutine carrier_exhausted_record(n, step_first, step_last,       &
                                          jworst, icworst, ratio_worst,   &
                                          terms_phys, terms_full)
      integer,  intent(out) :: n, step_first, step_last, jworst, icworst
      real(dp), intent(out) :: ratio_worst, terms_phys, terms_full
      n           = carrier_intervals_exhausted
      step_first  = carrier_exhausted_step_first
      step_last   = carrier_exhausted_step_last
      jworst      = carrier_exhausted_worst_j
      icworst     = carrier_exhausted_worst_ic
      ratio_worst = carrier_exhausted_worst_ratio
      terms_phys  = carrier_exhausted_worst_phys
      terms_full  = carrier_exhausted_worst_full
      end subroutine carrier_exhausted_record

      ! WHETHER THE CARRIER HISTORY OF THIS RUN CAN BE CERTIFIED.  False
      ! once an interval was left uncovered and the march went on from the
      ! entry carriers: no step of that history solved the carrier system
      ! over that interval, and no later step repairs it.  A certification
      ! of the FINAL state is a different question, taken on that state's
      ! own rows, and both belong in the carrier entry of the report.
      logical function carrier_history_certifiable() result(ok)
      ok = carrier_history_certifiable_flag
      end function carrier_history_certifiable

      ! The history mark put back to its start, for a test driver that
      ! marks it in one section and needs an unmarked run in the next.  A
      ! run never calls this: a history once marked stays marked.
      subroutine carrier_history_reset_for_test()
      carrier_history_certifiable_flag = .true.
      end subroutine carrier_history_reset_for_test

      ! HOW OFTEN DOUBLE PRECISION, AND NOT THE PHYSICS, SET THE ACCURACY
      ! OF AN ACCEPTED ROW.  nlast is the count in the last accepted
      ! substep, ntotal the sum over the run and nsubsteps how many
      ! accepted substeps carried at least one such row.  Informational:
      ! these rows are accepted, and a nonzero count says that the state
      ! they belong to was certified on rows one of which could not be
      ! solved any further at that substep.
      subroutine carrier_roundoff_limited_record(nlast, ntotal, nsubsteps)
      integer, intent(out) :: nlast, ntotal, nsubsteps
      nlast     = carrier_roundoff_rows_last
      ntotal    = carrier_roundoff_rows_total
      nsubsteps = carrier_roundoff_substeps
      end subroutine carrier_roundoff_limited_record

      ! ------------------------------------------------------------- !

      ! The rows the last verdict refused, at most carrier_refused_rows_max
      ! of them, and the two scales each was measured against.  A row is
      ! here only if its imbalance is larger than the round-off of the
      ! arithmetic that formed it, so the pair of scales says how far the
      ! refusal is from the limit of the precision.
      integer function carrier_refused_rows() result(n)
      n = n_refused_reported
      end function carrier_refused_rows

      subroutine carrier_refused_row(k, j, ic, absres, terms_phys,        &
                                     terms_full)
      integer,  intent(in)  :: k
      integer,  intent(out) :: j, ic
      real(dp), intent(out) :: absres, terms_phys, terms_full
      j = 0;  ic = 1
      absres = 0.0d0;  terms_phys = 0.0d0;  terms_full = 0.0d0
      if (k .lt. 1 .or. k .gt. n_refused_reported) return
      j          = refused_j(k)
      ic         = refused_ic(k)
      absres     = refused_res(k)
      terms_phys = refused_phys(k)
      terms_full = refused_full(k)
      end subroutine carrier_refused_row

      ! ------------------------------------------------------------- !

      ! Keep the rows a refusal rests on, with both scales, before the
      ! restore erases the assembly they came from.  A row that the verdict
      ! classified round-off limited is not one of them: it was accepted.
      subroutine carrier_record_refused_rows(res, terms_phys, terms_full, &
                                             roundoff, constrained)
      real(dp), dimension(:,:), intent(in) :: res, terms_phys, terms_full
      logical,  dimension(:,:), intent(in) :: roundoff
      logical,  dimension(:),   intent(in) :: constrained
      integer :: j, ic, ncell, ncar
      n_refused_reported = 0
      ncell = min(size(res,1), size(constrained))
      ncar  = min(n_carrier, size(res,2))
      do j = 1, ncell
         if (constrained(j)) cycle
         do ic = 1, ncar
            if (.not. carrier_solved(ic)) cycle
            if (roundoff(j,ic)) cycle
            if (abs(res(j,ic)) .le. carrier_accept_tol                    &
                                    *max(terms_phys(j,ic), 1.0d-300))     &
               cycle
            if (n_refused_reported .ge. carrier_refused_rows_max) return
            n_refused_reported = n_refused_reported + 1
            refused_j(n_refused_reported)    = j
            refused_ic(n_refused_reported)   = ic
            refused_res(n_refused_reported)  = abs(res(j,ic))
            refused_phys(n_refused_reported) = terms_phys(j,ic)
            refused_full(n_refused_reported) = terms_full(j,ic)
         enddo
      enddo
      end subroutine carrier_record_refused_rows

      ! Molecular diffusion coefficient of carrier ic in cell j [cm^2/s],
      ! from the last transport step; zero before the first one.
      double precision function carrier_diffusion_coefficient(j, ic)      &
                                result(D)
      integer, intent(in) :: j, ic
      if (allocated(pct_Dco)) then
         D = pct_Dco(j,ic)
      else
         D = 0.0d0
      endif
      end function carrier_diffusion_coefficient

      ! ------------------------------------------------------------- !

      ! The coupled carrier transport-chemistry system advanced by backward
      ! Euler over the interval dt_code of each cell (adimensional, as the
      ! hydro's own): in one step where the returned state is accepted, and
      ! otherwise in accepted substeps that cover the same interval.  Updates
      ! f_sp in place, and only once the whole interval is covered; a no-op
      ! unless the molecular chemistry is on with its carrier transport, and
      ! unless ioniz_eq has run at least once (the frozen background comes
      ! from there).
      !
      ! WHAT AN UNCOVERED INTERVAL MEANS DEPENDS ON WHAT THE RUN IS DOING
      ! (a0_run_mode_contract_20260906 sec. 2), and the run says which in
      ! its input.
      !
      !  * PHYSICAL INTEGRATION.  Every accepted step satisfies its
      !    balances, so a carrier interval that no substep could cover is
      !    a step that cannot be accepted.  The entry state is back, the
      !    status says exhausted, and the run ends here until the step
      !    controller around this operator turns that status into a
      !    rejection of the whole attempted step.
      !  * INITIALIZATION AND CONTINUATION.  A non-root iterate is a
      !    numerical guess and nothing physical accumulates, so the march
      !    goes on from the entry carriers, which are frozen over this
      !    step: that is a guess this mode allows, and it is the same
      !    entry state a covered interval would have started from.  What
      !    is left of the failure is the record: the attempts ledger
      !    counts the interval with its step and its worst row, one line
      !    names it, and the run's carrier history is marked not
      !    certifiable, because no step of it solved the carrier system
      !    over that interval.  A certification of the FINAL state is
      !    unaffected and still reads the carrier rows of that state.
      !
      ! carrier_interval_status is the value the caller reads instead of
      ! inferring the outcome from a stop.
      subroutine photochemical_transport_step(rho, v, f_sp, dt_code,      &
                                              carrier_interval_status,   &
                                              trial)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_code
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      integer, intent(out), optional :: carrier_interval_status
      ! A TRIAL is a step the caller keeps only if it likes the result: the
      ! stationary relaxation takes one, reads the status, and undoes the
      ! composition when the interval was not covered or the state left its
      ! bound.  A trial that fails is therefore not an interval the run
      ! adopted without integrating it, and it leaves the carrier history
      ! of the run untouched; the attempts ledger still counts it.  A step
      ! the marching loop takes is never a trial.
      logical, intent(in), optional :: trial
      logical :: is_trial

      type(carrier_verdict) :: verdict
      logical  :: completed
      real(dp) :: frac_done
      integer  :: nsub
      integer  :: krow, nrow, jrow, icrow
      real(dp) :: rres, rphys, rfull
      integer  :: jw, icw
      real(dp) :: rw, pw, fw
      logical  :: inv_carried(n_species)

      ! THE STATUS IS DEFINED ON EVERY RETURN PATH, and it is set here,
      ! before anything can return, rather than on each of them.  The
      ! configuration that carries no carriers has nothing to integrate and
      ! leaves through the first return below: its interval is covered by
      ! definition, and a caller reads that instead of an undefined value
      ! it would have had to pre-set itself.
      is_trial = .false.
      if (present(trial)) is_trial = trial
      carrier_interval_status_last = carrier_interval_covered
      if (present(carrier_interval_status))                               &
           carrier_interval_status = carrier_interval_covered
      ! This step's verdict starts at "no interval": what the last interval
      ! of the run decided is not an answer about this step.
      carrier_step_verdict_val = carrier_verdict(accepted = .false.,      &
                                        reason = carrier_no_interval)

      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      ! nsub is the number of substeps attempted, so nsub = 0 is the step
      ! that integrated nothing and whose verdict stays "no interval".
      if (nsub .gt. 0) carrier_step_verdict_val = verdict
      ! THE ELEMENTAL CONSTRAINT OF THE STATE THIS STEP WROTE, reported and
      ! not enforced here (EXHALE_INVENTORY_REPORT=1, default off).  The
      ! budgets the step spent were formed at the state it entered while
      ! its carriers moved, so the shared constraint of the written state
      ! is a statement about the step and not an identity.
      if (element_inventory_report_on()) then
         call carrier_species_mask(inv_carried)
         call element_inventory_report(                                   &
              'the state the carrier transport step wrote', rho, f_sp,    &
              inv_carried)
      endif
      if (completed) return
      carrier_interval_status_last = carrier_interval_exhausted
      if (present(carrier_interval_status))                               &
           carrier_interval_status = carrier_interval_exhausted

      ! THE ATTEMPTS LEDGER TAKES THE INTERVAL IN EITHER MODE: an attempt
      ! is what happened, and it is never rolled back.  The worst row is
      ! the largest imbalance relative to its own physical terms among the
      ! rows this refusal rested on.
      call exhausted_worst_row(verdict, jw, icw, rw, pw, fw)
      carrier_intervals_exhausted = carrier_intervals_exhausted + 1
      if (carrier_exhausted_step_first .lt. 0)                            &
           carrier_exhausted_step_first = count
      carrier_exhausted_step_last = count
      if (rw .gt. carrier_exhausted_worst_ratio) then
         carrier_exhausted_worst_ratio = rw
         carrier_exhausted_worst_j     = jw
         carrier_exhausted_worst_ic    = icw
         carrier_exhausted_worst_phys  = pw
         carrier_exhausted_worst_full  = fw
      endif

      ! A trial is undone by its caller: no history mark, no message, no
      ! stop, whatever the run mode.
      if (is_trial) return
      if (run_mode .ne. run_mode_phys) then
         ! Initialization and continuation: one line, and the march goes
         ! on from the entry carriers.  The mark is what a certification
         ! of this run's carrier history reads.
         carrier_history_certifiable_flag = .false.
         write(*,'(a,i0,a,a,a,i0,a,es12.4,a,es12.4,a,es12.4,a)')          &
            ' (photochemical_transport_step) step ', count,               &
            ': the carrier interval was NOT covered; worst row ',         &
            trim(carrier_name(icw)), ' at j = ', jw, ', |res| ', rw*pw,   &
            ' against physical terms ', pw, ' and full terms ', fw,       &
            '; the carriers of this step are frozen at their entry'//     &
            ' values (initialization mode)'
         return
      endif

      ! A STATE THAT SATISFIES NEITHER ITS EQUATIONS NOR A CONSTRAINT IS
      ! NOT A STATE THIS CODE MAY CONTINUE FROM.
      !
      ! The transported partition enters the free atomic hydrogen the next
      ! ionization sweep closes on, the chemical heat of the molecular
      ! layer, and the H2 opacity; a partition that is a solution of nothing
      ! puts that into every one of them without a mark.  The controller has
      ! already refused the returned state at every substep down to the
      ! floor and has put the entry state back, so the interval is NOT
      ! covered and there is nothing left to continue from.
      !
      ! THE DECISION IS THE VERDICT'S, taken on the returned state's own
      ! rows, and the raw status below is printed for diagnosis only.
      write(*,'(a)') ' (photochemical_transport_step) the carrier'//      &
         ' transport-chemistry state was refused'
      write(*,'(a)') '   why: '//trim(carrier_verdict_reason_text(        &
         verdict))
      write(*,'(a,es12.4,a,a,a,i0)')                                      &
         '   worst unconstrained row of the returned state ',             &
         verdict%rworst, ' in carrier ',                                  &
         trim(carrier_name(verdict%icworst)), ' at cell j = ',            &
         verdict%jworst
      write(*,'(a,es9.2)') '   the returned state is accepted below ',    &
         carrier_accept_tol
      write(*,'(a,i0,a,i0,a)') '   element constraints moved ',           &
         verdict%n_cells_constrained, ' cells (',                         &
         verdict%n_rows_constrained, ' rows), whose rows are the'//       &
         ' constraint'//''''//'s and are not certified here'
      write(*,'(a,i0,a,i0,a)') '   unconstrained rows above the'//        &
         ' floor: ', verdict%n_rows_above, ', of which ',                 &
         verdict%n_rows_above_beside_constrained, ' next to a'//          &
         ' constrained cell'
      write(*,'(a,i0,a)') '   rows accepted as round-off limited: ',      &
         verdict%n_rows_roundoff_limited, ' (their residual is what'//    &
         ' forming the row costs in double precision)'
      ! THE TWO SCALES OF EVERY ROW THE REFUSAL RESTS ON.  The acceptance
      ! is taken on the physical terms; whether anything could have been
      ! done about the row is decided by its full terms, whose round-off
      ! is the residual an exact solve leaves.  Printing both is what
      ! separates a state that is out of balance from a row the arithmetic
      ! cannot resolve, without either being argued from the other.
      nrow = carrier_refused_rows()
      if (nrow .gt. 0) then
         write(*,'(a)') '   the rows it rests on, with the residual and'//&
            ' the two scales it is measured against:'
         write(*,'(a)') '      carrier    cell   |residual|   physical'// &
            ' terms    full terms   |res|/phys   |res|/full'
         do krow = 1, nrow
            call carrier_refused_row(krow, jrow, icrow, rres, rphys,      &
                                     rfull)
            write(*,'(a,a8,i7,5es14.4)') '      ',                        &
               trim(carrier_name(icrow)), jrow, rres, rphys, rfull,       &
               rres/max(rphys, 1.0d-300), rres/max(rfull, 1.0d-300)
         enddo
         write(*,'(a,es9.2,a)') '   a row is round-off limited, and'//    &
            ' accepted, at |res|/full terms below ',                      &
            carrier_row_roundoff, '; these are above it'
      endif
      write(*,'(a)') '   how the iteration ended: '//                     &
         trim(carrier_solve_status_text(verdict%raw_status))
      write(*,'(a,es12.4,a,es12.4)')                                      &
         '   row residual before the element limiter ',                   &
         pct_newton_resid_before_limit, ', of the state returned ',       &
         pct_newton_resid
      write(*,'(a,i0,a,i0,a)') '   the retry controller made ', nsub,     &
         ' attempts and halved the substep down to dt/2**',               &
         pct_retry_halvings, ', and refused the state each time'
      write(*,'(a,f12.9,a)')                                              &
         '   THE INTERVAL WAS NOT COMPLETED.  Accepted substeps had'//    &
         ' covered ', frac_done, ' of dt before the failure;'
      write(*,'(a)') '   they were discarded with the entry state, and'// &
         ' the composition f_sp was not written, so'
      write(*,'(a)') '   no partial interval has entered the run.'
      write(*,'(a)') '   See docs/PLAN_20260906_rev2.md step A3.'
      write(*,'(a)') '   Physical integration accepts only a step that'// &
         ' satisfies its balances, so this'
      write(*,'(a)') '   step cannot be accepted.'
      if (carrier_transport_stop_on_failure) error stop 1
      carrier_transport_stops_suppressed =                                &
           carrier_transport_stops_suppressed + 1

      end subroutine photochemical_transport_step

      ! ------------------------------------------------------------- !

      ! THE WORST ROW A REFUSAL RESTED ON, with both scales it was measured
      ! against: the largest |res|/physical terms among the rows kept by
      ! carrier_record_refused_rows.  Where no row was kept -- a refusal
      ! whose reason is the iteration itself and not a row -- the verdict's
      ! own worst row is returned and the two scales are zero, because the
      ! assembly they came from has been restored.
      subroutine exhausted_worst_row(verdict, j, ic, ratio, terms_phys,   &
                                     terms_full)
      type(carrier_verdict), intent(in) :: verdict
      integer,  intent(out) :: j, ic
      real(dp), intent(out) :: ratio, terms_phys, terms_full
      integer  :: k, n, jk, ick
      real(dp) :: res, phys, full, rat
      ! The verdict is the argument and not pct_verdict, which the entry
      ! restore has already put back to the assembly of the previous step.
      j = verdict%jworst;  ic = verdict%icworst
      ratio = verdict%rworst
      terms_phys = 0.0d0;  terms_full = 0.0d0
      n = carrier_refused_rows()
      if (n .le. 0) return
      ratio = 0.0d0
      do k = 1, n
         call carrier_refused_row(k, jk, ick, res, phys, full)
         rat = res/max(phys, 1.0d-300)
         if (rat .gt. ratio) then
            ratio = rat;  j = jk;  ic = ick
            terms_phys = phys;  terms_full = full
         endif
      enddo
      end subroutine exhausted_worst_row

      ! ------------------------------------------------------------- !

      ! COVER THE INTERVAL dt_code WITH ACCEPTED CARRIER SUBSTEPS, or cover
      ! none of it.  This is the body of the step above with the stop taken
      ! out, so that the caller owns the stop and a test can read what the
      ! controller did.  completed is the only statement that the interval
      ! was integrated; frac_done is how much of it the accepted substeps
      ! covered, and it is 1 exactly when completed is true.  Where the
      ! interval is not covered, frac_done is what the accepted substeps
      ! had reached before the failure, and those are discarded with the
      ! entry state, so nothing of the interval survives.
      !
      ! WHAT THE SUBSTEPS HOLD FROZEN, AND WHAT FOLLOWS THE CARRIER STATE.
      ! The retry halves the substep of the carrier system alone; it does
      ! not re-enter the operators around it.  Frozen at the state the
      ! interval was entered with, and identical in every substep:
      !   * the background composition -- the ion stages, the molecular
      !     ions and the free electrons of cbg_* -- and the free element
      !     densities nH_free, nO_free, nC_free the carriers share with it;
      !   * the temperature TK, the total density ntot, nrho, wfac and the
      !     mean mass mbar, all of them the frozen cell state's;
      !   * the radiation the layer sees: the photolysis rates cph_klw,
      !     cph_jh2o and cph_joh, and with them the shielding columns they
      !     were built from;
      !   * the transport coefficients Dco, Agrd, Bdrf and the upwind
      !     directions updrf, the geometry rp and rep, the gravity gphys
      !     and the signal rate;
      !   * the face mass flux the advective term rides on and the mixture
      !     mass it is converted with.
      ! Following the carrier state, and rebuilt at every substep from the
      ! state that substep starts at:
      !   * the carrier densities in the reaction terms of carrier_source,
      !     and through the free-atom closures every coupled row of the cell;
      !   * the advective and diffusive face fluxes, which are those
      !     densities on the frozen coefficients above;
      !   * the time term, whose old state is the previous accepted substep.
      ! So the substeps integrate the carrier transport-chemistry system at a
      ! fixed background, which is the operator one full step already
      ! applies: subdividing refines its time integration and does not change
      ! the splitting.  The splitting error against the operators around it
      ! is the error of freezing that background over dt and is unchanged by
      ! the subdivision.  src/tests/carrier_retry/ measures the
      ! time-integration error by step doubling and prints it.
      !
      ! LABELED PARTIAL RECOVERY.  The checkpoint is the carrier
      ! subsystem's.  A refused substep leaves the hydrodynamic, ionization
      ! and energy state of the surrounding step untouched because this
      ! operator does not write them; it does not restore them either, and
      ! makes no claim to.  The enumerated checkpoint of the whole attempted
      ! step is B3a (b1_target_system_20260906_draft sec. 7, contract T7.3).
      subroutine carrier_transport_interval(rho, v, f_sp, dt_code,        &
                                            completed, frac_done, nsub,   &
                                            verdict)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_code
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      logical,               intent(out) :: completed
      real(dp),              intent(out) :: frac_done
      integer,               intent(out) :: nsub
      type(carrier_verdict), intent(out) :: verdict

      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fc, fc_old, fctry
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco
      real(dp), dimension(1-Ng:N+Ng)           :: ntot, TK, gphys, mbar
      real(dp), dimension(1-Ng:N+Ng)           :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng)           :: rp, rep, dt_phys
      real(dp), dimension(1-Ng:N+Ng)           :: sigrate
      ! The mixture mass and the face mass flux the advective term of the
      ! rows would ride on.  The marching rows carry no advective term (it
      ! is taken with the mass row inside the Runge-Kutta stages), so they
      ! are formed here only to be handed on unread.
      real(dp), dimension(1-Ng:N+Ng)           :: msum, Frho
      real(dp), dimension(1-Ng:N+Ng)           :: nH_free, nO_free, nC_free
      real(dp), dimension(0:N,n_carrier_max)       :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max)       :: updrf
      type(carrier_checkpoint) :: chk_entry, chk_attempt
      real(dp) :: sub, step, frac_left
      integer  :: nhalve, naccepted
      logical  :: forced
      ! THE CONSTRAINT DIAGNOSTICS OF THE WHOLE INTERVAL, over its ACCEPTED
      ! substeps and no others.  solve_carriers resets pct_limited,
      ! pct_worst_limit and pct_cell_constrained at every call, so after
      ! several accepted substeps those describe the last one and not the
      ! step the run took.  What the element constraints did to the state
      ! the interval hands on is the sum of the cell applications, the
      ! largest overshoot removed and the cells any of them moved.
      integer  :: acc_limited
      real(dp) :: acc_worst_limit
      logical, allocatable :: acc_constrained(:)
      ! A1 element-budget assertions.  The transport step holds rho fixed and
      ! its write-back restores the ENTRY element totals cell by cell (that is
      ! what nH_free / nO_free / nC_free are), so the absolute nucleus density
      ! of every element is invariant across both the write-back alone and the
      ! whole step.  Off unless EXHALE_ELEMENT_ASSERT is set.
      type(element_census_state) :: cen_step, cen_wb

      ! Nothing to integrate is a covered interval, not a failure.
      completed = .true.
      frac_done = 1.0d0
      nsub      = 0
      if (.not. thereis_mol)       return
      if (.not. carrier_transport) return
      if (.not. bg_ready)          return
      completed = .false.
      frac_done = 0.0d0

      ! The entry checkpoint is taken before anything is built, so a
      ! permanent failure returns the module to the state the PREVIOUS step
      ! left and not to a half-built one.  The carrier fractions are not
      ! part of it: at this point they still live in f_sp, which this
      ! routine writes only once the interval is complete.
      call carrier_checkpoint_take(chk_entry)

      call element_census_take('photochemical_transport_step',            &
                               rho, f_sp, cen_step)
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK,            &
                         mbar, nH_free, nO_free, nC_free)
      fc_old = fc

      call carrier_geometry(dt_code, rp, rep, dt_phys, gphys)
      call carrier_advective_state(carrier_rows_advect, rho, f_sp,     &
                                   msum, Frho)
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      if (.not. allocated(pct_Dco)) allocate(pct_Dco(1-Ng:N+Ng,n_carrier_max))
      pct_Dco = Dco
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)
      call carrier_photolysis(rho, TK, f_sp)
      call carrier_signal_rate(v, TK, mbar, sigrate)

      ! THE DOMAIN OF THE ONE-SIDED CO MODEL, taken once per interval on the
      ! state the interval starts from -- which is the state the previous
      ! accepted step handed on -- and with the rates carrier_photolysis has
      ! just rebuilt on that state's own columns.  It is taken here and not
      ! after the substeps for two reasons.  It is then taken whether or not
      ! the interval is covered, and a relaxation that never covers one is
      ! exactly the case in which a reader needs to know whether the CO model
      ! was in domain.  And the entry checkpoint, which a failed interval
      ! restores, was taken before this routine allocated the frozen
      ! background, so after that restore there is nothing to read.
      ! It reads only; nothing here changes a density or a rate.
      call carrier_co_domain_take(fc, nrho, TK, v, rp, nO_free)

      frac_left       = 1.0d0
      sub             = max(min(carrier_initial_substep_for_test, 1.0d0), &
                            carrier_substep_floor)
      nhalve             = 0
      pct_retry_halvings = 0
      naccepted          = 0
      acc_limited     = 0
      acc_worst_limit = 0.0d0
      do
         step = min(sub, frac_left)
         call carrier_checkpoint_take(chk_attempt, fc)
         fctry = fc
         call solve_carriers(fctry, fc, nrho, wfac, step*dt_phys, rp,     &
                             rep, msum, Frho, Agrd, Bdrf, updrf, TK,      &
                             sigrate, nH_free, nO_free, nC_free, verdict)
         nsub = nsub + 1
         carrier_substeps_attempted = carrier_substeps_attempted + 1
         forced = (nsub .le. carrier_reject_leading_attempts_for_test)
         if (forced) then
            verdict%accepted = .false.
            verdict%reason   = carrier_reject_forced_for_test
            pct_verdict      = verdict
         endif
         if (verdict%accepted) then
            fc        = fctry
            frac_left = frac_left - step
            naccepted = naccepted + 1
            carrier_substeps_accepted = carrier_substeps_accepted + 1
            acc_limited     = acc_limited    + pct_limited
            acc_worst_limit = max(acc_worst_limit, pct_worst_limit)
            ! Rows this accepted substep could not be solved any further
            ! at, because their residual is the round-off of their own
            ! terms.  Counted, never a refusal.
            carrier_roundoff_rows_last = verdict%n_rows_roundoff_limited
            if (verdict%n_rows_roundoff_limited .gt. 0) then
               carrier_roundoff_rows_total = carrier_roundoff_rows_total  &
                                        + verdict%n_rows_roundoff_limited
               carrier_roundoff_substeps = carrier_roundoff_substeps + 1
            endif
            if (allocated(pct_cell_constrained)) then
               if (.not. allocated(acc_constrained)) then
                  acc_constrained = pct_cell_constrained
               else
                  acc_constrained = acc_constrained .or. pct_cell_constrained
               endif
            endif
            if (frac_left .le. 0.0d0) exit
         else
            ! Back to the state this substep started from: a refused trial
            ! leaves no contribution, in the carrier fractions or in the
            ! accepted records.
            call carrier_checkpoint_restore(chk_attempt, fc)
            carrier_substeps_rejected = carrier_substeps_rejected + 1
            if (nhalve .eq. 0) carrier_intervals_retried =                &
                                    carrier_intervals_retried + 1
            ! EVERY REFUSAL LEFT IS ONE A SHORTER SUBSTEP MAY ANSWER: a
            ! row at the round-off of its own terms no longer refuses, it
            ! is accepted and counted, so the controller halves only for
            ! rows that carry more imbalance than the arithmetic explains.
            if (nhalve .ge. carrier_retry_max) then
               ! The substep is already at dt/2**carrier_retry_max.  The
               ! entry state goes back and the interval is not covered.
               frac_done = 1.0d0 - frac_left
               call carrier_checkpoint_restore(chk_entry)
               return
            endif
            sub                = 0.5d0*sub
            nhalve             = nhalve + 1
            pct_retry_halvings = nhalve
         endif
      enddo
      completed = .true.
      frac_done = 1.0d0

      ! One substep leaves the diagnostics exactly as it wrote them, so the
      ! single-substep path is untouched.
      if (naccepted .gt. 1) then
         pct_limited     = acc_limited
         pct_worst_limit = acc_worst_limit
         if (allocated(acc_constrained))                                  &
              pct_cell_constrained = acc_constrained
      endif

      call element_census_take('carrier_write_back', rho, f_sp, cen_wb)
      call carrier_write_back(rho, f_sp, fc, nrho, nH_free, nO_free,     &
                              nC_free)
      call element_census_verify(cen_wb,   rho, f_sp, rho_is_fixed=.true.)
      call element_census_verify(cen_step, rho, f_sp, rho_is_fixed=.true.)

      if (carrier_debug_on()) then
         write(*,'(a,i4,a,es9.2,a,f6.0,a,i0,a,es9.2,a,es9.2,'//         &
              'a,es9.2,a,es9.2)')                                        &
            ' (carrier_transport) newton ', pct_newton_steps,            &
            ' resid ', pct_newton_resid, ' at cell '//                  &
            trim(carrier_name(pct_worst_ic))//' j=', dble(pct_worst_j),  &
            ' limited ', pct_limited,                                    &
            ' cells, worst ', pct_worst_limit,                           &
            '; base x_H2 ',                                              &
            2.0d0*fc(1,ic_H2)*nrho(1)/max(nH_free(1), 1.0d-99),          &
            ' <- ', 2.0d0*fc_old(1,ic_H2)*nrho(1)                        &
                    /max(nH_free(1), 1.0d-99),                           &
            ' dt ', dt_phys(1)
         write(*,'(a,es12.5,a,es12.5,a,i0)')                             &
            ' (carrier_transport) row residual before the element'//     &
            ' limiter ', pct_newton_resid_before_limit,                  &
            ', of the state returned ', pct_newton_resid,                &
            '; limited cells ', pct_limited
         ! The accepted state, said in the terms it was accepted on: the
         ! worst UNCONSTRAINED row, which is the one the verdict measures,
         ! beside the way the iteration itself ended.  A step that ended by
         ! stagnating at the relative drop is visible here and nowhere else.
         write(*,'(a,es12.4,a,a,a,i0,a,i0,a)')                           &
            ' (carrier_transport) worst unconstrained row ',             &
            verdict%rworst, ' in ',                                      &
            trim(carrier_name(verdict%icworst)), ' at j = ',             &
            verdict%jworst, '; ends: '//                                 &
            trim(carrier_solve_status_text(verdict%raw_status))//        &
            '; constrained cells ', verdict%n_cells_constrained, ''
         if (nsub .gt. 1) write(*,'(a,i0,a)')                            &
            ' (carrier_transport) the interval was covered in ', nsub,   &
            ' attempts'
      endif

      end subroutine carrier_transport_interval

      ! EXHALE_CARRIER_DEBUG=1 turns on one line per step from the transport
      ! operator: Newton iterations, the relative residual, how many cells
      ! the element limiter touched, and how far the base partition moved.
      logical function carrier_debug_on()
      character(len=8) :: env
      call get_environment_variable('EXHALE_CARRIER_DEBUG', env)
      carrier_debug_on = (trim(env) .eq. '1')
      end function carrier_debug_on

      ! ------------------------------------------------------------- !

      ! THE SPECIES THIS RUN'S CARRIER OPERATOR TRANSPORTS, as a mask over
      ! the f_sp columns, so that the elemental map is asked about the set
      ! carrier_set_init fixed and not about a wider one.  A species outside
      ! the set holds its nuclei on the frozen side of the conditional
      ! budget.
      subroutine carrier_species_mask(carried)
      logical, intent(out) :: carried(n_species)
      integer :: ic
      carried = .false.
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         carried(carrier_species_index(ic)) = .true.
      enddo
      end subroutine carrier_species_mask

      ! ------------------------------------------------------------- !

      ! The carrier mixing ratios of the current state, and the element
      ! headroom each of them is limited against.
      !
      ! nH_free is the H nuclei NOT held by the ion stages and the molecular
      ! ions, i.e. what H I and the transported carriers share; nO_free and
      ! nC_free are the same for oxygen and carbon against their ion stages.
      ! They are the element totals minus the parts this step cannot move, so
      ! a carrier that stays inside them cannot break an element budget.
      !
      ! THE TEMPERATURE IS THE FROZEN BACKGROUND'S, not a separate argument:
      ! TK, ntot and everything built from them come from bg_cell, the cell
      ! states the last ionization sweep left, so the chemistry, the
      ! diffusivities and the settling term of one transport step all describe
      ! one gas. A caller cannot hand this routine a temperature the sweep has
      ! not seen.
      subroutine carrier_state(rho, f_sp, fc, ntot, nrho, wfac,          &
                               TK, mbar, nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(out) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: ntot, TK
      ! nrho = rho*n0, the density the carrier fraction is a fraction OF --
      ! the one the hydrodynamic mass row conserves exactly.  wfac converts
      ! that fraction into the PARTICLE mixing ratio the diffusive flux is
      ! driven by: X = wfac*y, wfac = rho*n0/n_tot = mbar/m_amu.
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: mbar
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nC_free
      real(dp), dimension(1-Ng:N+Ng) :: nd
      real(dp), dimension(1-Ng:N+Ng,4) :: nmol_l
      real(dp), dimension(1-Ng:N+Ng,n_mion) :: nm_l
      ! The advected carriers, in the order they were declared, which is the
      ! order of the solved columns of fc.
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fcadv
      integer :: j, i0, k, im, ic, ncar

      if (.not. allocated(cbg_nhii)) then
         allocate(cbg_nhii(1-Ng:N+Ng),  cbg_nh2p(1-Ng:N+Ng),             &
                  cbg_nh3p(1-Ng:N+Ng),  cbg_nhehp(1-Ng:N+Ng),            &
                  cbg_nhei(1-Ng:N+Ng),  cbg_nheii(1-Ng:N+Ng),            &
                  cbg_nheiii(1-Ng:N+Ng),cbg_nheiTR(1-Ng:N+Ng),           &
                  cbg_ne(1-Ng:N+Ng),                                     &
                  cbg_nOion(1-Ng:N+Ng), cbg_nCion(1-Ng:N+Ng))
         allocate(cph_klw(1-Ng:N+Ng), cph_kco(1-Ng:N+Ng),                &
                  cph_jh2o(1-Ng:N+Ng,n_fuv_band),                        &
                  cph_joh(1-Ng:N+Ng,n_fuv_band))
         cph_klw  = 0.0d0
         cph_kco  = 0.0d0
         cph_jh2o = 0.0d0
         cph_joh  = 0.0d0
      endif

      nd = rho*n0                      ! cm^-3 per unit of f_sp

      ! The frozen background of this step.
      cbg_nhii   = f_sp(:,isp_HII)  *nd
      cbg_nh2p   = f_sp(:,isp_H2p)  *nd
      cbg_nh3p   = f_sp(:,isp_H3p)  *nd
      cbg_nhehp  = f_sp(:,isp_HeHp) *nd
      cbg_nhei   = f_sp(:,isp_HeI)  *nd
      cbg_nheii  = f_sp(:,isp_HeII) *nd
      cbg_nheiii = f_sp(:,isp_HeIII)*nd
      cbg_nheiTR = 0.0d0
      if (thereis_HeITR) cbg_nheiTR = f_sp(:,isp_HeTR)*nd
      do im = 1, n_mion
         nm_l(:,im) = f_sp(:,mion_fsp(im))
      enddo
      nmol_l(:,1) = f_sp(:,isp_H2)
      nmol_l(:,2) = f_sp(:,isp_H2p)
      nmol_l(:,3) = f_sp(:,isp_H3p)
      nmol_l(:,4) = f_sp(:,isp_HeHp)
      call calc_ne(f_sp(:,isp_HII), f_sp(:,isp_HeII), f_sp(:,isp_HeIII), &
                   cbg_ne, nm_l, nmol_l)
      cbg_ne = cbg_ne*nd
      cbg_nOion = 0.0d0
      i0 = melem_i0(iel_O)
      do k = 1, melem_top(iel_O)
         cbg_nOion = cbg_nOion + f_sp(:,mion_fsp(i0+k))*nd
      enddo
      cbg_nCion = 0.0d0
      i0 = melem_i0(iel_C)
      do k = 1, melem_top(iel_C)
         cbg_nCion = cbg_nCion + f_sp(:,mion_fsp(i0+k))*nd
      enddo

      ! Gas-particle density and the mean mass per particle.  Both come from
      ! the cell state ioniz_eq last solved, so they are the same numbers the
      ! EOS is using; freezing them inside one step is stated in sec. 4.
      do j = 1-Ng, N+Ng
         ntot(j) = bg_cell(j)%ntot
         TK(j)   = bg_cell(j)%T_K
      enddo
      do j = 1-Ng, N+Ng
         if (ntot(j) .gt. 0.0d0) then
            mbar(j) = rho(j)*n0*mu/ntot(j)
         else
            mbar(j) = mu
         endif
      enddo

      ! THE CARRIER UNKNOWN IS THE FRACTION PER UNIT MASS, not per particle.
      ! Section 158: the marching loop's hydrodynamic stage moves rho at
      ! frozen f_sp, so the carrier density it hands on changes by
      ! n_c d(ln rho)/dt.  A carrier advected as n_c/n_tot instead composes
      ! with that stage into the conservation law PLUS a spurious
      ! n_c v d(ln mbar)/dr, which lives exactly across the dissociation and
      ! ionization fronts where mbar varies -- measured at 13 times the row's
      ! own steady residual in the wind.  n_c/(rho n0) is the variable the
      ! mass row conserves, and it is the variable f_sp already stores.
      nrho = nd
      wfac = nd/max(ntot, 1.0d-99)
      fc(:,ic_H2)  = f_sp(:,isp_H2)
      fc(:,ic_OH)  = f_sp(:,isp_OH)
      fc(:,ic_H2O) = f_sp(:,isp_H2O)
      fc(:,ic_CO)  = f_sp(:,isp_CO)
      ! The proton column is filled whether or not it is solved -- it costs
      ! one assignment and it keeps fc a complete picture of the state -- but
      ! only a run with ionization_transport lets the solve move it.
      fc(:,ic_Hp)  = f_sp(:,isp_HII)

      ! THE ADVECTED CARRIERS ARE THE STATE THIS STEP STARTS FROM.  The
      ! material advection of every solved carrier was taken with the mass
      ! row inside the Runge-Kutta stages, and what stands in the species
      ! vector at this point is the composition BEFORE that advection: the
      ! carriers are not written into f_sp between the two, because moving a
      ! carrier there without moving the stages it competes with would leave
      ! its element total wrong until the write-back below restores it.  So
      ! the advected fractions are read here, and the write-back at the end
      ! of the step puts carrier and element totals back together in one
      ! place.
      ! It is a state of ONE marching step, so it is taken up once and
      ! marked spent; and the fixed-wind relaxation, whose own rows carry
      ! the advective term, reads the species vector like every other path
      ! that runs between two marching steps.
      if (advected_carrier_state_ready() .and.                           &
          .not. carrier_rows_advect) then
         ncar = advected_carrier_count()
         call advected_carrier_fractions(f_sp, fcadv(:,1:ncar))
         k = 0
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            k = k + 1
            fc(:,ic) = fcadv(:,k)
         enddo
         call advected_carrier_state_consumed()
      endif
      where (fc .lt. 0.0d0) fc = 0.0d0

      ! ---- element headroom -------------------------------------------
      ! What the transported carriers of each element may hold.  It is the
      ! ELEMENT TOTAL, not the neutral stage: transport moves nuclei, and
      ! which stage a nucleus sits in is re-solved from scratch by the very
      ! next ionization sweep.  Limiting the carriers to the neutral stage
      ! instead would cap them wherever the element happens to be ionized --
      ! measured on the HD 209458 b molecular gate, that clamped 150-240 of
      ! 503 cells by up to 2% every step, which is a boundary of the
      ! bookkeeping and not of the physics.
      !
      ! THAT ARGUMENT DOES NOT REACH AS FAR AS IT WAS TAKEN TO, and hydrogen
      ! is where it breaks.  "The next sweep re-solves the stages" is true of
      ! the sweep; it is not true of the step in between.  This step's
      ! write-back does not re-solve H I, H II, H2+ and H3+ -- it rescales
      ! them in proportion to what they already held -- so nuclei the
      ! carriers take are nuclei those stages simply lose, and a carrier
      ! holding the whole element leaves the sweep a cell with no electrons.
      ! So the ELEMENT total is the right headroom for the LIMITER's outer
      ! question (which element, and how much of it exists), and
      ! hydrogen_available_to_carriers is the right one for what the
      ! carriers may actually take.  nH_free below is the first; the limiter
      ! and the chemistry rows both go through the second.
      !
      ! The one nucleus this step genuinely cannot move is the H in HeH+,
      ! because that molecule also carries a helium nucleus and rescaling it
      ! would move helium.  It is therefore taken out of the hydrogen
      ! headroom rather than shared.
      !
      ! The three densities are the element_free of the one stoichiometric
      ! map (element_inventory), which reads the nucleus counts of the
      ! species table, so no call site of this operator carries a
      ! stoichiometric coefficient of its own.
      call carrier_element_totals(nd, f_sp, nH_free, nO_free, nC_free)

      end subroutine carrier_state

      ! ------------------------------------------------------------- !

      ! The row scales of the last carrier residual assembly: terms_full is
      ! the sum of the magnitudes of every term of the row, the time term
      ! included, and terms_physical the same sum without it (declaration of
      ! row_terms_phys).  Both come back zero where no assembly has run.
      subroutine carrier_row_scales(terms_full, terms_physical)
      real(dp), dimension(1:N,n_carrier_max), intent(out) :: terms_full,  &
                                                             terms_physical
      terms_full     = 0.0d0
      terms_physical = 0.0d0
      if (allocated(row_terms))      terms_full     = row_terms(1:N,:)
      if (allocated(row_terms_phys)) terms_physical = row_terms_phys(1:N,:)
      end subroutine carrier_row_scales

      ! ------------------------------------------------------------- !

      ! Cell and carrier that the last fixed-wind relaxation moved most.
      subroutine carrier_drift_location(j, ic)
      integer, intent(out) :: j, ic
      j  = pct_drift_j
      ic = pct_drift_ic
      end subroutine carrier_drift_location

      ! ------------------------------------------------------------- !

      subroutine save_carrier_module_state(s)
      type(carrier_module_state), intent(out) :: s
      if (allocated(cbg_nhii))     s%cbg_nhii     = cbg_nhii
      if (allocated(cbg_nh2p))     s%cbg_nh2p     = cbg_nh2p
      if (allocated(cbg_nh3p))     s%cbg_nh3p     = cbg_nh3p
      if (allocated(cbg_nhehp))    s%cbg_nhehp    = cbg_nhehp
      if (allocated(cbg_nhei))     s%cbg_nhei     = cbg_nhei
      if (allocated(cbg_nheii))    s%cbg_nheii    = cbg_nheii
      if (allocated(cbg_nheiii))   s%cbg_nheiii   = cbg_nheiii
      if (allocated(cbg_nheiTR))   s%cbg_nheiTR   = cbg_nheiTR
      if (allocated(cbg_ne))       s%cbg_ne       = cbg_ne
      if (allocated(cbg_nOion))    s%cbg_nOion    = cbg_nOion
      if (allocated(cbg_nCion))    s%cbg_nCion    = cbg_nCion
      if (allocated(cph_klw))      s%cph_klw      = cph_klw
      if (allocated(cph_kco))      s%cph_kco      = cph_kco
      if (allocated(cph_jh2o))     s%cph_jh2o     = cph_jh2o
      if (allocated(cph_joh))      s%cph_joh      = cph_joh
      if (allocated(row_terms))    s%row_terms    = row_terms
      if (allocated(row_terms_phys))                                     &
                                   s%row_terms_phys = row_terms_phys
      if (allocated(col_scale_car)) s%col_scale_car = col_scale_car
      if (allocated(row_scale_car)) s%row_scale_car = row_scale_car
      if (allocated(headroom_car))  s%headroom_car  = headroom_car
      s%headroom_set   = headroom_set
      s%pct_worst_j    = pct_worst_j
      s%pct_worst_ic   = pct_worst_ic
      end subroutine save_carrier_module_state

      subroutine restore_carrier_module_state(s)
      type(carrier_module_state), intent(in) :: s
      call put_back_1d(s%cbg_nhii,     cbg_nhii)
      call put_back_1d(s%cbg_nh2p,     cbg_nh2p)
      call put_back_1d(s%cbg_nh3p,     cbg_nh3p)
      call put_back_1d(s%cbg_nhehp,    cbg_nhehp)
      call put_back_1d(s%cbg_nhei,     cbg_nhei)
      call put_back_1d(s%cbg_nheii,    cbg_nheii)
      call put_back_1d(s%cbg_nheiii,   cbg_nheiii)
      call put_back_1d(s%cbg_nheiTR,   cbg_nheiTR)
      call put_back_1d(s%cbg_ne,       cbg_ne)
      call put_back_1d(s%cbg_nOion,    cbg_nOion)
      call put_back_1d(s%cbg_nCion,    cbg_nCion)
      call put_back_1d(s%cph_klw,      cph_klw)
      call put_back_1d(s%cph_kco,      cph_kco)
      call put_back_2d(s%cph_jh2o,     cph_jh2o)
      call put_back_2d(s%cph_joh,      cph_joh)
      call put_back_2d(s%row_terms,    row_terms)
      call put_back_2d(s%row_terms_phys, row_terms_phys)
      call put_back_2d(s%col_scale_car, col_scale_car)
      call put_back_2d(s%row_scale_car, row_scale_car)
      call put_back_2d(s%headroom_car,  headroom_car)
      headroom_set   = s%headroom_set
      pct_worst_j    = s%pct_worst_j
      pct_worst_ic   = s%pct_worst_ic
      end subroutine restore_carrier_module_state

      logical function carrier_module_state_matches(s) result(ok)
      ! Whether the module arrays now hold EXACTLY what was held aside in s,
      ! allocation status included. The round trip of an isolated evaluation
      ! is asserted with this, so that "the measurement did not change the
      ! state" is a check the run makes and not a claim a comment makes.
      type(carrier_module_state), intent(in) :: s
      ok = same_1d(s%cbg_nhii,     cbg_nhii)
      ok = ok .and. same_1d(s%cbg_nh2p,     cbg_nh2p)
      ok = ok .and. same_1d(s%cbg_nh3p,     cbg_nh3p)
      ok = ok .and. same_1d(s%cbg_nhehp,    cbg_nhehp)
      ok = ok .and. same_1d(s%cbg_nhei,     cbg_nhei)
      ok = ok .and. same_1d(s%cbg_nheii,    cbg_nheii)
      ok = ok .and. same_1d(s%cbg_nheiii,   cbg_nheiii)
      ok = ok .and. same_1d(s%cbg_nheiTR,   cbg_nheiTR)
      ok = ok .and. same_1d(s%cbg_ne,       cbg_ne)
      ok = ok .and. same_1d(s%cbg_nOion,    cbg_nOion)
      ok = ok .and. same_1d(s%cbg_nCion,    cbg_nCion)
      ok = ok .and. same_1d(s%cph_klw,      cph_klw)
      ok = ok .and. same_1d(s%cph_kco,      cph_kco)
      ok = ok .and. same_2d(s%cph_jh2o,     cph_jh2o)
      ok = ok .and. same_2d(s%cph_joh,      cph_joh)
      ok = ok .and. same_2d(s%row_terms,    row_terms)
      ok = ok .and. same_2d(s%row_terms_phys, row_terms_phys)
      ok = ok .and. same_2d(s%col_scale_car, col_scale_car)
      ok = ok .and. same_2d(s%row_scale_car, row_scale_car)
      ok = ok .and. same_2d(s%headroom_car,  headroom_car)
      ok = ok .and. (s%headroom_set .eqv. headroom_set)
      ok = ok .and. (s%pct_worst_j  .eq. pct_worst_j)
      ok = ok .and. (s%pct_worst_ic .eq. pct_worst_ic)
      end function carrier_module_state_matches

      logical function same_1d(kept, live) result(ok)
      real(dp), allocatable, intent(in) :: kept(:), live(:)
      ok = (allocated(kept) .eqv. allocated(live))
      if (.not. ok) return
      if (.not. allocated(kept)) return
      ok = (size(kept) .eq. size(live))
      if (ok) ok = all(kept .eq. live)
      end function same_1d

      logical function same_2d(kept, live) result(ok)
      real(dp), allocatable, intent(in) :: kept(:,:), live(:,:)
      ok = (allocated(kept) .eqv. allocated(live))
      if (.not. ok) return
      if (.not. allocated(kept)) return
      ok = (size(kept,1) .eq. size(live,1)) .and.                         &
           (size(kept,2) .eq. size(live,2))
      if (ok) ok = all(kept .eq. live)
      end function same_2d

      subroutine put_back_1d(kept, live)
      real(dp), allocatable, intent(in)    :: kept(:)
      real(dp), allocatable, intent(inout) :: live(:)
      if (allocated(kept)) then
         live = kept
      else if (allocated(live)) then
         deallocate(live)
      endif
      end subroutine put_back_1d

      subroutine put_back_2d(kept, live)
      real(dp), allocatable, intent(in)    :: kept(:,:)
      real(dp), allocatable, intent(inout) :: live(:,:)
      if (allocated(kept)) then
         live = kept
      else if (allocated(live)) then
         deallocate(live)
      endif
      end subroutine put_back_2d

      subroutine put_back_1d_logical(kept, live)
      logical, allocatable, intent(in)    :: kept(:)
      logical, allocatable, intent(inout) :: live(:)
      if (allocated(kept)) then
         live = kept
      else if (allocated(live)) then
         deallocate(live)
      endif
      end subroutine put_back_1d_logical

      logical function same_1d_logical(kept, live) result(ok)
      logical, allocatable, intent(in) :: kept(:), live(:)
      ok = (allocated(kept) .eqv. allocated(live))
      if (.not. ok) return
      if (.not. allocated(kept)) return
      ok = (size(kept) .eq. size(live))
      if (ok) ok = all(kept .eqv. live)
      end function same_1d_logical

      subroutine put_back_2d_logical(kept, live)
      logical, allocatable, intent(in)    :: kept(:,:)
      logical, allocatable, intent(inout) :: live(:,:)
      if (allocated(kept)) then
         live = kept
      else if (allocated(live)) then
         deallocate(live)
      endif
      end subroutine put_back_2d_logical

      logical function same_2d_logical(kept, live) result(ok)
      logical, allocatable, intent(in) :: kept(:,:), live(:,:)
      ok = (allocated(kept) .eqv. allocated(live))
      if (.not. ok) return
      if (.not. allocated(kept)) return
      ok = (size(kept,1) .eq. size(live,1))                               &
           .and. (size(kept,2) .eq. size(live,2))
      if (ok) ok = all(kept .eqv. live)
      end function same_2d_logical

      ! ------------------------------------------------------------- !

      ! HOLD ASIDE EVERYTHING ONE CARRIER TRANSPORT ATTEMPT WRITES, and put
      ! it back.  The list is the module's own state, item by item: the
      ! arrays a residual assembly overwrites (through the pair A2 wrote,
      ! so that one enumeration serves both), the diffusion coefficients of
      ! the step, and the last call's constraint diagnostics and verdict.
      ! Allocation status is part of the state: an array that was not allocated before the
      ! attempt is deallocated again, so the module returns to the state it
      ! was in and not merely to the same numbers.
      !
      ! fc is optional because the entry state of the carriers is f_sp: at
      ! the top of a step the local carrier array does not hold a state yet.
      subroutine carrier_checkpoint_take(chk, fc)
      type(carrier_checkpoint), intent(out) :: chk
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in), optional :: fc
      if (present(fc)) chk%fc = fc
      call save_carrier_module_state(chk%bg)
      if (allocated(pct_Dco))              chk%pct_Dco = pct_Dco
      if (allocated(pct_cell_constrained))                                &
                       chk%pct_cell_constrained = pct_cell_constrained
      chk%pct_newton_steps              = pct_newton_steps
      chk%pct_newton_resid              = pct_newton_resid
      chk%pct_newton_resid_before_limit = pct_newton_resid_before_limit
      chk%pct_limited            = pct_limited
      chk%pct_worst_limit        = pct_worst_limit
      chk%pct_drift_j            = pct_drift_j
      chk%pct_drift_ic           = pct_drift_ic
      chk%pct_verdict            = pct_verdict
      end subroutine carrier_checkpoint_take

      subroutine carrier_checkpoint_restore(chk, fc)
      type(carrier_checkpoint), intent(in) :: chk
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(out), optional :: fc
      if (present(fc)) then
         if (allocated(chk%fc)) fc = chk%fc
      endif
      call restore_carrier_module_state(chk%bg)
      call put_back_2d(chk%pct_Dco, pct_Dco)
      call put_back_1d_logical(chk%pct_cell_constrained,                  &
                               pct_cell_constrained)
      pct_newton_steps              = chk%pct_newton_steps
      pct_newton_resid              = chk%pct_newton_resid
      pct_newton_resid_before_limit = chk%pct_newton_resid_before_limit
      pct_limited            = chk%pct_limited
      pct_worst_limit        = chk%pct_worst_limit
      pct_drift_j            = chk%pct_drift_j
      pct_drift_ic           = chk%pct_drift_ic
      pct_verdict            = chk%pct_verdict
      end subroutine carrier_checkpoint_restore

      ! Whether the module now holds EXACTLY what the checkpoint kept,
      ! allocation status included.  Restoration is measured with this and
      ! not asserted by a comment.
      logical function carrier_checkpoint_matches(chk, fc) result(ok)
      type(carrier_checkpoint), intent(in) :: chk
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in), optional :: fc
      ok = carrier_module_state_matches(chk%bg)
      ok = ok .and. same_2d(chk%pct_Dco, pct_Dco)
      ok = ok .and. same_1d_logical(chk%pct_cell_constrained,             &
                                    pct_cell_constrained)
      ok = ok .and. (chk%pct_newton_steps .eq. pct_newton_steps)
      ok = ok .and. (chk%pct_newton_resid .eq. pct_newton_resid)
      ok = ok .and. (chk%pct_newton_resid_before_limit .eq.               &
                     pct_newton_resid_before_limit)
      ok = ok .and. (chk%pct_limited     .eq. pct_limited)
      ok = ok .and. (chk%pct_worst_limit .eq. pct_worst_limit)
      ok = ok .and. (chk%pct_drift_j     .eq. pct_drift_j)
      ok = ok .and. (chk%pct_drift_ic    .eq. pct_drift_ic)
      ok = ok .and. verdicts_agree(chk%pct_verdict, pct_verdict)
      if (present(fc) .and. ok) then
         if (allocated(chk%fc)) ok = all(chk%fc .eq. fc)
      endif
      end function carrier_checkpoint_matches

      logical function verdicts_agree(a, b) result(ok)
      type(carrier_verdict), intent(in) :: a, b
      ok = (a%accepted .eqv. b%accepted)                                  &
           .and. (a%reason  .eq. b%reason)                                &
           .and. (a%jworst  .eq. b%jworst)                                &
           .and. (a%icworst .eq. b%icworst)                               &
           .and. (a%rworst  .eq. b%rworst)                                &
           .and. (a%n_cells_constrained .eq. b%n_cells_constrained)       &
           .and. (a%n_rows_constrained  .eq. b%n_rows_constrained)        &
           .and. (a%n_rows_above        .eq. b%n_rows_above)              &
           .and. (a%n_rows_above_beside_constrained .eq.                  &
                  b%n_rows_above_beside_constrained)                      &
           .and. (a%n_rows_roundoff_limited .eq.                          &
                  b%n_rows_roundoff_limited)                              &
           .and. (a%raw_status .eq. b%raw_status)                         &
           .and. (a%stall      .eq. b%stall)
      end function verdicts_agree

      ! ------------------------------------------------------------- !

      ! WRITE A DISTINCT CHANGE INTO EVERY ITEM THE CHECKPOINT COVERS.
      ! src/tests/carrier_retry/ states restoration as a measurement over
      ! the whole list, and not over the items one particular column
      ! happens to move: an item left out of carrier_checkpoint_take or of
      ! carrier_checkpoint_restore is then a failing assertion.  An array
      ! that is not allocated is allocated here, so the restore is asked to
      ! give the module its allocation status back and not only its
      ! numbers.  Test only: no production path calls this.
      subroutine carrier_perturb_checkpointed_state_for_test(fc)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout),        &
                optional :: fc
      if (present(fc)) fc = fc + 1.0d0
      call bump_1d(cbg_nhii,   1-Ng, N+Ng)
      call bump_1d(cbg_nh2p,   1-Ng, N+Ng)
      call bump_1d(cbg_nh3p,   1-Ng, N+Ng)
      call bump_1d(cbg_nhehp,  1-Ng, N+Ng)
      call bump_1d(cbg_nhei,   1-Ng, N+Ng)
      call bump_1d(cbg_nheii,  1-Ng, N+Ng)
      call bump_1d(cbg_nheiii, 1-Ng, N+Ng)
      call bump_1d(cbg_nheiTR, 1-Ng, N+Ng)
      call bump_1d(cbg_ne,     1-Ng, N+Ng)
      call bump_1d(cbg_nOion,  1-Ng, N+Ng)
      call bump_1d(cbg_nCion,  1-Ng, N+Ng)
      call bump_1d(cph_klw,    1-Ng, N+Ng)
      call bump_1d(cph_kco,    1-Ng, N+Ng)
      call bump_2d(cph_jh2o,   1-Ng, N+Ng, 1, n_fuv_band)
      call bump_2d(cph_joh,    1-Ng, N+Ng, 1, n_fuv_band)
      call bump_2d(row_terms,  1, N, 1, n_carrier_max)
      call bump_2d(row_terms_phys, 1, N, 1, n_carrier_max)
      call bump_2d(col_scale_car, 1, N, 1, n_carrier_max)
      call bump_2d(row_scale_car, 1, N, 1, n_carrier_max)
      call bump_2d(headroom_car,  1, N, 1, n_carrier_max)
      call bump_2d(pct_Dco,    1-Ng, N+Ng, 1, n_carrier_max)
      if (.not. allocated(pct_cell_constrained))                          &
           allocate(pct_cell_constrained(1:N+Ng))
      pct_cell_constrained = .not. pct_cell_constrained
      headroom_set   = .not. headroom_set
      pct_worst_j    = pct_worst_j  + 3
      pct_worst_ic   = pct_worst_ic + 1
      pct_newton_steps              = pct_newton_steps + 5
      pct_newton_resid              = pct_newton_resid + 0.25d0
      pct_newton_resid_before_limit = pct_newton_resid_before_limit       &
                                      + 0.5d0
      pct_limited     = pct_limited    + 7
      pct_worst_limit = pct_worst_limit + 0.75d0
      pct_drift_j     = pct_drift_j  + 13
      pct_drift_ic    = pct_drift_ic + 1
      pct_verdict%accepted   = .not. pct_verdict%accepted
      pct_verdict%reason     = pct_verdict%reason  + 1
      pct_verdict%jworst     = pct_verdict%jworst  + 17
      pct_verdict%icworst    = pct_verdict%icworst + 1
      pct_verdict%rworst     = pct_verdict%rworst  + 0.375d0
      pct_verdict%n_cells_constrained = pct_verdict%n_cells_constrained + 2
      pct_verdict%n_rows_constrained  = pct_verdict%n_rows_constrained  + 3
      pct_verdict%n_rows_above        = pct_verdict%n_rows_above        + 4
      pct_verdict%n_rows_above_beside_constrained =                       &
           pct_verdict%n_rows_above_beside_constrained + 5
      pct_verdict%raw_status = pct_verdict%raw_status + 1
      pct_verdict%stall      = pct_verdict%stall      + 1
      end subroutine carrier_perturb_checkpointed_state_for_test

      subroutine bump_1d(a, lo, hi)
      real(dp), allocatable, intent(inout) :: a(:)
      integer, intent(in) :: lo, hi
      integer :: k
      if (.not. allocated(a)) then
         allocate(a(lo:hi))
         a = 0.0d0
      endif
      do k = lbound(a,1), ubound(a,1)
         a(k) = a(k) + 1.0d0 + 0.5d0*dble(k)
      enddo
      end subroutine bump_1d

      subroutine bump_2d(a, lo1, hi1, lo2, hi2)
      real(dp), allocatable, intent(inout) :: a(:,:)
      integer, intent(in) :: lo1, hi1, lo2, hi2
      integer :: k, m
      if (.not. allocated(a)) then
         allocate(a(lo1:hi1,lo2:hi2))
         a = 0.0d0
      endif
      do m = lbound(a,2), ubound(a,2)
         do k = lbound(a,1), ubound(a,1)
            a(k,m) = a(k,m) + 1.0d0 + 0.5d0*dble(k) + 0.25d0*dble(m)
         enddo
      enddo
      end subroutine bump_2d

      ! Is the state (rho, v, f_sp) a steady state OF THE CARRIER EQUATION,
      ! and where is it least so? The temperature of that state is the one
      ! the frozen background carries (see carrier_state), so it is not a
      ! separate argument.
      !
      ! The steady solver never sees this equation: it holds the carrier
      ! partition fixed and solves the wind, and an outer Picard iteration
      ! relaxes the carriers at that wind (EXHALE_main). Their common fixed
      ! point is a steady state of both, but ONLY the hydrodynamic half of it
      ! is measured by ||R||. This is the other half.
      !
      ! It is the residual the transport step already assembles, evaluated
      ! with the time term switched off (dt -> infinity with fc_old = fc
      ! makes it identically zero), divided by carrier_row_scale -- the same
      ! place, and the same principle, as every hydrodynamic row scale
      ! (steady_residual.f90). A no-op, reported as zero, when the carriers
      ! are not transported.
      subroutine carrier_steady_residual(rho, v, f_sp, rcmax,            &
                                         jworst, icworst, rvol, rlegacy, &
                                         res_out, terms_out)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho, v
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), intent(out) :: rcmax
      integer,  intent(out) :: jworst, icworst
      ! The volume-weighted measure, the one a gate should read: the larger
      ! of the layer's and the wind's ratio of sums, exactly as
      ! relnorm_over_cells forms the hydrodynamic rows.
      real(dp), intent(out), optional :: rvol
      ! The same residual on the retired n_H(|v|+c_s)/dr scale, so that a run
      ! can report both while the two are being compared.
      real(dp), intent(out), optional :: rlegacy
      real(dp), dimension(1:N,n_carrier_max), intent(out), optional ::   &
                                                 res_out, terms_out

      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fc, Dco
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar, gphys
      real(dp), dimension(1-Ng:N+Ng) :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng) :: rp, rep, dt_big, sigrate
      real(dp), dimension(1-Ng:N+Ng) :: msum, Frho
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp), dimension(0:N,n_carrier_max) :: Agrd, Bdrf, Jf, dJl, dJr
      integer,  dimension(0:N,n_carrier_max) :: updrf
      real(dp), dimension(1:N,n_carrier_max) :: res
      real(dp) :: rnorm_unused, cs, sc, rel, w, nH_avail
      real(dp) :: num_l, den_l, num_w, den_w
      integer  :: j, ic
      logical  :: inv_carried(n_species)

      rcmax = 0.0d0;  jworst = 0;  icworst = 1
      if (present(rvol))      rvol      = 0.0d0
      if (present(rlegacy))   rlegacy   = 0.0d0
      if (present(res_out))   res_out   = 0.0d0
      if (present(terms_out)) terms_out = 0.0d0
      if (.not. thereis_mol)       return
      if (.not. carrier_transport) return
      if (.not. bg_ready)          return

      if (.not. allocated(row_terms)) allocate(row_terms(1:N,n_carrier_max))
      if (.not. allocated(row_terms_phys))                               &
         allocate(row_terms_phys(1:N,n_carrier_max))
      row_terms      = 0.0d0
      row_terms_phys = 0.0d0
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK,            &
                         mbar, nH_free, nO_free, nC_free)
      call carrier_geometry(spread(1.0d0,1,N+2*Ng), rp, rep, dt_big,     &
                            gphys)
      call carrier_advective_state(.true., rho, f_sp, msum, Frho)
      dt_big = 1.0d30                 ! kills the time term exactly
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)
      call carrier_photolysis(rho, TK, f_sp)
      call carrier_signal_rate(v, TK, mbar, sigrate)
      call carrier_residual(fc, fc, nrho, wfac, dt_big, rp, rep, msum,   &
                            Frho, .true.,                                &
                            Agrd, Bdrf, updrf, nH_free, nO_free,         &
                            nC_free, sigrate, res, Jf, dJl, dJr,         &
                            rnorm_unused)

      ! THE SCALE IS THE ROW'S OWN LARGEST TERMS, and that is a change.
      ! It used to be carrier_row_scale, n_H(|v|+c_s)/dr -- the flux
      ! divergence the cell would carry if its WHOLE hydrogen inventory moved
      ! at the signal speed. That is a bound on the row, but it is not the
      ! row's own term, and measured on the He/H = 0.0793 hot Uranus it
      ! exceeds the sum of the terms the row actually balances by 10^2.7 to
      ! 10^6.4 (docs/p50_carrier_wind_alternation.md section 2.4). A state
      ! the steady solver accepted at info = 0 then reported a carrier
      ! residual of 1.6e-4 while carrying 0.60 of its own terms at the H2
      ! front, and the outer loop had no way to see that the state it was
      ! handing over was not a carrier steady state at all.
      !
      ! Section 133 already made "a bound on that row's own largest term" the
      ! principle for every hydrodynamic row. This is the carrier row on the
      ! same principle, and the quantity it needs is the one
      ! carrier_residual has always computed for its own Newton: row_terms.
      !
      ! THE SCALE HERE IS row_terms_phys, the same sum without the time
      ! term.  This routine asks whether a state is STATIONARY, a question
      ! about transport against reaction in which no step appears, so its
      ! scale may not carry one either: dividing by a scale that grows like
      ! 1/dt would report a state as more nearly stationary the shorter the
      ! step that happened to produce it.  The two agree here to the size of
      ! the time term this call leaves, dt = 1e30 (below), which is 1e-30 of
      ! the row; what the physical scale removes is the residue and the
      ! step form of the absolute floor.  The two Newton scales below are the
      ! Newton's own and keep row_terms.
      if (.not. allocated(col_scale_car))                                &
         allocate(col_scale_car(1:N,n_carrier_max),                      &
                  row_scale_car(1:N,n_carrier_max),                      &
                  headroom_car(1:N,n_carrier_max))
      num_l = 0.0d0;  den_l = 0.0d0;  num_w = 0.0d0;  den_w = 0.0d0
      do j = 1, N
         cs = sqrt(adiabatic_index_at_T(j, TK(j)/T0)*kb_erg                &
                  *max(TK(j), 1.0d0)/max(mbar(j), 1.0d-40))
         sc = carrier_row_scale(j, nH_free(j), v(j)*v0, cs)
         w  = r(j)*r(j)*dr_j(j)
         nH_avail = hydrogen_available_to_carriers(j, nH_free(j))
         do ic = 1, n_carrier_max
            row_scale_car(j,ic) = row_terms(j,ic)*R0/v0
            col_scale_car(j,ic) = fc(j,ic)*nrho(j)                       &
                            + row_terms(j,ic)/max(sigrate(j), 1.0d-300)
         enddo
         ! THE ELEMENT HEADROOM IS FROZEN AT THE OUTER ITERATE, and the
         ! reference matters. It is a constraint the solver tests trial
         ! states against, and it depends on the very unknown it constrains
         ! -- n_H_free counts H2's own two nuclei -- so evaluating it at the
         ! trial makes it move with the step. Refreshing it at every state
         ! that is not a Jacobian probe was tried and is worse than that:
         ! the last LINE-SEARCH trial then leaves its own headroom behind,
         ! and the Gauss-Newton escape that follows finds the ITERATE
         ! outside it, so every backtracked step is refused down to
         ! lam = 5e-4 and the solve aborts holding a descent direction
         ! (measured). weno_mode = 1 is the one evaluation per outer
         ! iteration that is the iterate itself.
         if (weno_mode .eq. 1 .or. .not. headroom_set) then
            do ic = 1, n_carrier_max
               headroom_car(j,ic) = carrier_element_headroom(ic, nH_avail, &
                                          nO_free(j), nC_free(j))
            enddo
         endif
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            rel = abs(res(j,ic))/max(row_terms_phys(j,ic), 1.0d-300)
            if (rel .gt. rcmax) then
               rcmax = rel;  jworst = j;  icworst = ic
            endif
            if (present(rlegacy)) rlegacy = max(rlegacy,                 &
                                     abs(res(j,ic))/max(sc, 1.0d-300))
            if (j .lt. j_min) then
               num_l = num_l + abs(res(j,ic))*w
               den_l = den_l + row_terms_phys(j,ic)*w
            else
               num_w = num_w + abs(res(j,ic))*w
               den_w = den_w + row_terms_phys(j,ic)*w
            endif
         enddo
      enddo
      ! The two regions are normed separately and combined by the larger,
      ! for the reason residual_norms gives: the wind carries almost all of
      ! sum r^2 dr, so any single sum lets it average the inner column away.
      if (present(rvol)) rvol = max(num_l/max(den_l, 1.0d-300),          &
                                    num_w/max(den_w, 1.0d-300))
      headroom_set = .true.
      ! THE ELEMENTAL CONSTRAINT OF THE ITERATE THE BUDGETS WERE FROZEN AT
      ! (EXHALE_INVENTORY_REPORT=1, default off).  weno_mode = 1 is the one
      ! evaluation of an outer iteration that is the iterate itself, so
      ! this is the state the box sides above describe, and the report is
      ! the statement of what those sides leave unbounded.
      if (element_inventory_report_on() .and. weno_mode .eq. 1) then
         call carrier_species_mask(inv_carried)
         call element_inventory_report(                                   &
              'the stationary iterate the element budget was frozen at',  &
              rho, f_sp, inv_carried)
      endif
      if (present(res_out))   res_out   = res(1:N,:)
      if (present(terms_out)) terms_out = row_terms_phys(1:N,:)
      end subroutine carrier_steady_residual

      ! ------------------------------------------------------------- !

      ! COLUMN scale of carrier ic's unknown in cell j [cm^-3]: the carrier
      ! density itself plus the amount of it a signal crossing would carry at
      ! the size of the row's own terms.  Exactly the construction of the
      ! momentum unknown's |rho v| + rho c_s -- the unknown where it means
      ! something, a floor that does not collapse where it does not.  It sets
      ! the Newton step and the finite-difference step, so it has to be the
      ! size of the variable; using the ROW scale for it instead was measured
      ! to propose carrier steps 84 times the cell's own n(H2) and to send the
      ! line search to lam = 4.9e-4.  Zero before the first assembly, which
      ! the caller reads as "no scale yet".
      double precision function carrier_column_scale(j, ic) result(D)
      integer, intent(in) :: j, ic
      D = 0.0d0
      if (allocated(col_scale_car)) D = col_scale_car(j,ic)
      end function carrier_column_scale

      ! ROW scale of carrier ic's row in cell j [cm^-3]: the row's own
      ! largest terms accumulated over one code time unit, so that the
      ! residual divided by it is the dimensionless imbalance the acceptance
      ! test uses.  It is NOT the column scale, and the two must differ:
      ! the row is stiff -- its natural rate is the chemical one, 1e3 per
      ! code time at the front -- while its unknown is small, so one vector
      ! cannot scale both without either drowning the merit or overshooting
      ! the step.
      double precision function carrier_row_term_scale(j, ic) result(D)
      integer, intent(in) :: j, ic
      D = 0.0d0
      if (allocated(row_scale_car)) D = row_scale_car(j,ic)
      end function carrier_row_term_scale

      ! The largest n(carrier ic) cell j may hold [cm^-3], from the free
      ! density of the carrier's own element and the nuclei of it the carrier
      ! holds.  A trial state above it is not a state the element budget
      ! admits.  The coupled steady solve makes it a FACE of the box its
      ! species unknowns live in and writes a trial onto it
      ! (species_unknowns_outside_their_bounds); the relaxation of this
      ! module refuses such a trial rather than clamping it, because
      ! clamping inside a Newton whose unknown IS the carrier puts a kink in
      ! the residual, which is the defect section 138 removed, and it is
      ! what made the fixed-wind relaxation return a state pinned on the
      ! ceiling in a quarter of the grid.
      !
      ! IT IS FROZEN AT THE OUTER ITERATE, not evaluated at the trial: see
      ! the assignment of headroom_car in carrier_steady_residual for why
      ! (it depends on the very unknown it constrains).  So a caller
      ! building a step reads one fixed number per cell and carrier for the
      ! whole of that step.
      double precision function carrier_headroom(j, ic) result(nmax)
      integer, intent(in) :: j, ic
      nmax = huge(1.0d0)
      if (allocated(headroom_car)) nmax = headroom_car(j,ic)
      end function carrier_headroom

      ! ------------------------------------------------------------- !

      ! WHETHER THE ELEMENT BUDGET EXISTS AT ALL, which is a question about
      ! the state and not about a magnitude.  The budget is a product of the
      ! carrier assembly, so before the first assembly there is none, and
      ! carrier_headroom then answers with huge(1) so that every comparison
      ! against it passes.  A caller that has to know the difference -- the
      ! coupled solve, which puts the budget under a face of its unknown box
      ! and cannot form a face out of huge(1) -- asks here rather than
      ! testing the returned number against a magnitude of its own.
      logical function carrier_headroom_known()
      carrier_headroom_known = allocated(headroom_car)
      end function carrier_headroom_known

      ! ------------------------------------------------------------- !

      ! Test only: no production path calls this.  The element budget is a
      ! property of a state this module has assembled, and
      ! src/tests/steady_species_rows/ states the projection onto the budget
      ! face of the coupled solve's unknown box on a SYNTHETIC cell, which
      ! has no such state.  It writes the budget [cm^-3] the face is then
      ! formed from.
      subroutine carrier_headroom_set_for_test(nmax)
      real(dp), dimension(:,:), intent(in) :: nmax
      if (allocated(headroom_car)) deallocate(headroom_car)
      allocate(headroom_car(1:size(nmax,1),1:size(nmax,2)))
      headroom_car = nmax
      headroom_set = .true.
      end subroutine carrier_headroom_set_for_test

      ! ------------------------------------------------------------- !

      ! THE LARGEST DENSITY CARRIER ic CAN REACH IN A CELL [cm^-3]: the free
      ! density of the element it is made of, divided by the nuclei of that
      ! element one carrier holds.  It is stoichiometry and not a
      ! convention, and it reads the same element table as the reference
      ! density and the derivative step (carrier_element_reference_density),
      ! so the three cannot state different elements for one carrier.
      !
      ! H2 holds two hydrogen nuclei and H+ one, and both are bounded by the
      ! hydrogen the carriers may take rather than by the element total: the
      ! stages this step holds frozen and rescales must be left something
      ! (hydrogen_available_to_carriers states what happens when they are
      ! not).  OH and H2O hold one oxygen nucleus each; CO needs one of
      ! oxygen and one of carbon, so the smaller free density bounds it.
      !
      ! IT IS ONE SIDE OF A BOX, DERIVED FROM THE MAP, AND NOT THE FEASIBLE
      ! SET.  Each side is the budget of ONE element divided by the nuclei
      ! of that element the carrier holds (element_box_side, whose
      ! multiplicities come from the species table); the element whose
      ! budget each carrier is charged to is stated here.  A side formed
      ! from one element bounds no other element the same carrier holds:
      ! the hydrogen of OH and of H2O is not bounded by their oxygen side,
      ! and the hydrogen sides of H2 and of H+ bound each carrier alone
      ! while the two spend one budget.  A state can therefore satisfy
      ! every side below and still break the shared element constraint,
      ! whose magnitude element_inventory reports.
      double precision function carrier_element_headroom(ic, nH_avail,    &
                               nO_free, nC_free) result(nmax)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nH_avail, nO_free, nC_free
      select case (ic)
      case (ic_H2)
         nmax = element_box_side(isp_H2, ien_H, nH_avail)
      case (ic_Hp)
         nmax = element_box_side(isp_HII, ien_H, nH_avail)
      case (ic_OH)
         nmax = element_box_side(isp_OH, ien_O, nO_free)
      case (ic_H2O)
         nmax = element_box_side(isp_H2O, ien_O, nO_free)
      case (ic_CO)
         nmax = element_box_side(isp_CO, ien_O, min(nO_free, nC_free))
      case default
         write(*,'(a,i0,a)') ' (carrier_element_headroom) carrier index ', &
            ic, ' belongs to no element'
         error stop 1
      end select
      end function carrier_element_headroom

      ! ------------------------------------------------------------- !

      ! The cell's signal-crossing rate (|v| + c_s)/dr [1/s], from the same
      ! background the carrier rows are assembled on.
      subroutine carrier_signal_rate(v, TK, mbar, sigrate)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: v, TK, mbar
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: sigrate
      real(dp) :: cs
      integer  :: j
      do j = 1-Ng, N+Ng
         cs = sqrt(adiabatic_index_at_T(j, TK(j)/T0)*kb_erg                &
                  *max(TK(j), 1.0d0)/max(mbar(j), 1.0d-40))
         sigrate(j) = (abs(v(j))*v0 + cs)                                &
                    / max((r_edg(min(max(j,1),N)) -                      &
                           r_edg(min(max(j,1),N)-1))*R0, 1.0d0)
      enddo
      end subroutine carrier_signal_rate

      ! ------------------------------------------------------------- !

      ! THE MATERIAL ADVECTION OF EVERY TRANSPORTED CARRIER, as the
      ! divergence of the face species fluxes the hydrodynamic stages carry.
      !
      !     F_c(j) = F_rho(j) Y_c^face(j),   Y_c = m_c f_c / msum,
      !     adv(j) = ( A_+ F_c(j) - A_- F_c(j-1) ) / dV_j  x  msum/m_c
      !
      ! with F_rho the face mass flux of the mass row of this state,
      ! Y_c^face the reconstructed and upwinded face mass fraction of
      ! species_face_fraction, and the divergence that of
      ! species_flux_divergence: the same three routines, on the same faces,
      ! areas and volumes, that species_advective_update calls inside a
      ! Runge-Kutta stage.  The factor msum/m_c returns the mass divergence
      ! to the units the carrier row is written in, n_c per unit volume and
      ! time [cm^-3 s^-1], and it is the inverse of the conversion
      ! carrier_mass_fractions applies on the way in.
      !
      ! WHY IT IS NOT A CELL-VELOCITY UPWIND DIFFERENCE ANY MORE.  Until
      ! item B5b the stationary rows and the fixed-wind relaxation carried
      ! n_tot v df/dr plus a deferred van Leer correction while the marching
      ! stages carried this divergence.  Two discretizations of one term have
      ! two fixed points, so the state a Newton converges on and the state
      ! the marching converges on were not the same object; the acceptance
      ! contract requires that they be.  The cell-velocity form also had the
      ! two failures the header of species_face_flux.f90 records: a cell
      ! whose two face fluxes straddle zero had no advective term at all,
      ! and a cell outflowing at both faces could be evacuated of a species
      ! its density does not lose.
      !
      ! The base face carries no special case: where the face mass flux
      ! flows inward the face composition is reconstructed from the ghosts,
      ! which hold the handoff partition where a handoff states one and the
      ! base cell's own where none does (carrier_mass_fractions).
      subroutine carrier_advective_divergence(fc, msum, Frho, adv, advmag)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: msum
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: Frho
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: adv
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: advmag

      real(dp), dimension(1-Ng:N+Ng) :: Y, Yf, Fs, dvF, dvM
      real(dp) :: mc, tconv
      integer  :: j, ic

      adv    = 0.0d0
      advmag = 0.0d0
      ! Code density times code velocity over code length is a rate in units
      ! of v0/R0, and f_sp counts n0 particles per unit density, so this is
      ! the factor that puts the divergence in cm^-3 s^-1 -- the same factor
      ! the time term of the row carries through nrho and dt_phys.
      tconv = n0*v0/R0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         mc = carrier_mass_amu(ic)
         call carrier_face_mass_fraction(fc, msum, ic, Y)
         call species_face_fraction(Y, Frho, Yf)
         call species_face_flux(Frho, Yf, Fs)
         call species_flux_divergence(Fs, dvF, dvM)
         do j = 1, N
            adv(j,ic)    = dvF(j)*msum(j)/mc*tconv
            advmag(j,ic) = dvM(j)*msum(j)/mc*tconv
         enddo
      enddo

      end subroutine carrier_advective_divergence

      ! ------------------------------------------------------------- !

      ! The two face coefficients of the advective term, in the row's own
      ! units per unit of the face MASS fraction:
      !
      !     advj(j) =  A_+ F_rho(j)  /dV_j ,
      !     advm(j) = -A_- F_rho(j-1)/dV_j ,
      !
      ! the same A_+, A_- and dV_j species_flux_divergence uses.  They are
      ! what the block-tridiagonal Jacobian of the fixed-wind relaxation
      ! differentiates the term with: multiplied by msum(j)/m_c they give
      ! d adv(j) / d Y(donor), and divided once more by the donor cell's own
      ! mixture mass they give d adv(j) / d f_c(donor).
      subroutine carrier_advective_face_coefficients(Frho, advj, advm)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: Frho
      real(dp), dimension(1:N),       intent(out) :: advj, advm
      real(dp) :: rpf, rmf, dAp, dAm, dV, tconv
      integer  :: j
      tconv = n0*v0/R0
      do j = 1, N
         rpf = r_edg(j)
         rmf = r_edg(j-1)
         dAp = rpf*rpf
         dAm = rmf*rmf
         dV  = (dAp*rpf - dAm*rmf)/3.0
         advj(j) =  dAp*Frho(j)  /dV*tconv
         advm(j) = -dAm*Frho(j-1)/dV*tconv
      enddo
      end subroutine carrier_advective_face_coefficients

      ! ------------------------------------------------------------- !

      ! One carrier as a mass fraction of the mixture, with the inflow
      ! composition of the inner ghosts: the handoff value where a handoff
      ! states one, the base cell's own where none does.  It is the same
      ! rule and the same expression carrier_mass_fractions applies to the
      ! declared set at the start of a marching step, written here on the
      ! carrier index this module solves in, so that the stationary rows and
      ! the stages put the same quantity on the faces.
      subroutine carrier_face_mass_fraction(fc, msum, ic, Y)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: msum
      integer,                                      intent(in)  :: ic
      real(dp), dimension(1-Ng:N+Ng),               intent(out) :: Y
      integer :: j
      Y = carrier_mass_amu(ic)*fc(:,ic)/msum
      if (.not. carrier_base_composition_imposed(ic)) then
         do j = 1-Ng, 0
            Y(j) = Y(1)
         enddo
      endif
      where (Y .lt. 0.0d0) Y = 0.0d0
      end subroutine carrier_face_mass_fraction

      ! ------------------------------------------------------------- !

      ! Whether something below the base states this carrier's partition.
      ! Only the molecular handoff does -- base.inp's q_H2_base, or the
      ! profile's value at the matching level -- and it states the H2
      ! partition alone; every other carrier keeps the zero-gradient
      ! condition of design decision D6.  It is the same predicate
      ! carrier_set_init registers the transported set with.
      logical function carrier_base_composition_imposed(ic) result(imposed)
      integer, intent(in) :: ic
      imposed = (ic .eq. ic_H2) .and. base_h2_composition_imposed()
      end function carrier_base_composition_imposed

      ! ------------------------------------------------------------- !

      ! The two state quantities the advective term of a carrier row is
      ! built from: the mixture mass per unit of f_sp, and the face mass
      ! flux of the MASS ROW of this state.  Both are frozen over the Newton
      ! that follows, exactly as the transport coefficients are: the trial
      ! states move the carrier partition and neither the wind nor the
      ! mixture mass of the frozen background moves with it.
      !
      ! Where the rows carry no advective term -- the marching operator,
      ! whose carriers were advected with the mass row inside the
      ! Runge-Kutta stages -- neither is read, and the face mass flux is not
      ! even asked for: in mid-step the stored flux belongs to the last
      ! stage and not to a state anyone is measuring.
      subroutine carrier_advective_state(advect, rho, f_sp, msum, Frho)
      logical,                                  intent(in)  :: advect
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: msum, Frho
      msum = 1.0d0
      Frho = 0.0d0
      if (.not. advect) return
      call mixture_mass_sum(f_sp, msum)
      call face_mass_flux_of_state(rho, Frho)
      end subroutine carrier_advective_state

      ! ------------------------------------------------------------- !

      ! Grid, step and gravity in physical units.
      subroutine carrier_geometry(dt_code, rp, rep, dt_phys, gphys)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: dt_code
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: rp, rep
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: dt_phys, gphys
      integer :: j
      rp      = r*R0
      rep     = r_edg*R0
      dt_phys = dt_code*R0/v0
      where (dt_phys .lt. 1.0d-30) dt_phys = 1.0d-30
      do j = 1-Ng, N+Ng
         gphys(j) = Dphi(r(j))*v0*v0/R0
      enddo
      end subroutine carrier_geometry

      ! ------------------------------------------------------------- !

      ! Molecular diffusion coefficient of each transported carrier in the
      ! background mixture, by Blanc's law over the rigid-sphere pairs
      ! (module header sec. 3).
      subroutine carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: ntot
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(out) :: Dco
      real(dp) :: y(n_bkg), ysum, fric, Dpair, nd, nc
      real(dp) :: bkg_mass(n_bkg)
      integer  :: j, ic, is
      do is = 1, n_bkg
         bkg_mass(is) = carrier_background_mass(is)
      enddo
      Dco = 0.0d0
      do j = 1-Ng, N+Ng
         nd = rho(j)*n0
         ysum = 0.0d0
         do is = 1, n_bkg
            y(is) = max(f_sp(j,bkg_isp(is))*nd, 0.0d0)
            ysum  = ysum + y(is)
         enddo
         if (ysum .le. 0.0d0) then
            Dco(j,:) = 0.0d0
            cycle
         endif
         y  = y/ysum
         nc = max(ntot(j), 1.0d0)
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            ! THE PROTON DIFFUSES AS AN ION, AND THAT IS NOT THIS OPTION.
            ! The hard-sphere pair coefficient below is the mutual diffusion
            ! of a NEUTRAL through neutrals. A proton in a mostly neutral
            ! gas is not that: it is coupled to the electrons by the
            ! ambipolar electric field, and to the other charges by Coulomb
            ! collisions whose cross section is orders of magnitude larger
            ! than the gas-kinetic one and is a function of the ionization
            ! fraction. Using the neutral coefficient for it would put a
            ! number in the equation that is not the physics of the term.
            !
            ! So the molecular diffusion of the proton is set to zero and it
            ! is carried by the bulk flow alone -- which is the term this
            ! option exists for: what Model A's Figure 15 shows carrying 43
            ! percent of the proton production at 2.0 r_base and 56
            ! percent at 2.4 is the ADVECTION, not a diffusive flux. The eddy
            ! coefficient K_zz is a separate matter and is still applied to
            ! this carrier (carrier_face_coefficients): it is a bulk mixing
            ! coefficient of the gas, blind to the charge of what it mixes.
            ! Adding the ambipolar coefficient is a physics decision of its
            ! own and is not taken here.
            if (ic .eq. ic_Hp) then
               Dco(j,ic) = 0.0d0
               cycle
            endif
            fric = 0.0d0
            do is = 1, n_bkg
               if (y(is) .le. 0.0d0) cycle
               Dpair = hard_sphere_pair_diffusion(max(TK(j), 1.0d0), nc, &
                            carrier_mass_amu(ic), bkg_mass(is))
               fric  = fric + y(is)/max(Dpair, 1.0d-99)
            enddo
            if (fric .gt. 0.0d0) then
               Dco(j,ic) = 1.0d0/fric
            else
               Dco(j,ic) = 0.0d0
            endif
         enddo
      enddo
      end subroutine carrier_diffusivities

      ! Mass of the species that occupies f_sp column isp, in units of the
      ! hydrogen atom, read from the species table so that this module and
      ! calc_rho cannot disagree about it. The row of the mass table is found
      ! through bsp_fsp, the table's own column map, so inserting a species
      ! into that table cannot silently change a mass here.
      double precision function species_mass_amu(isp) result(m)
      integer, intent(in) :: isp
      integer :: k
      ! A SPECIES WITH NO ROW IN THE MASS TABLE IS NOT A MASS OF ZERO.
      ! The mass divides the thermal speed of the hard-sphere pair
      ! diffusivity and multiplies the settling drift, so returning zero
      ! here would silently give the species an infinite mobility and no
      ! gravitational settling instead of stopping. There is no state of the
      ! species table in which a transported carrier or a collision partner
      ! of this module is absent from it, so reaching this is a defect in the
      ! table or in the species list that asked for it.
      m = -1.0d0
      do k = 1, n_bsp
         if (bsp_fsp(k) .eq. isp) m = bsp_mass(k)
      enddo
      if (m .le. 0.0d0) then
         write(*,'(a,i0,a)') ' (species_mass_amu) f_sp column ', isp,    &
            ' has no row in the species mass table'
         error stop 1
      endif
      end function species_mass_amu

      ! Mass [m_H] of background collision partner `is` of Blanc's law, the
      ! quantity that enters the reduced mass of every carrier/background
      ! pair.  It is the species table's mass of that partner.
      double precision function carrier_background_mass(is) result(m)
      integer, intent(in) :: is
      m = species_mass_amu(bkg_isp(is))
      end function carrier_background_mass

      ! THE f_sp COLUMN OF A TRANSPORTED CARRIER.  One table: the mass, the
      ! composition read-out (carrier_state) and the stationary system's
      ! unknown all name the same column through this function, so a carrier
      ! cannot be one species in one place and another somewhere else.
      integer function carrier_species_index(ic) result(isp)
      integer, intent(in) :: ic
      select case (ic)
      case (ic_H2)
         isp = isp_H2
      case (ic_OH)
         isp = isp_OH
      case (ic_H2O)
         isp = isp_H2O
      case (ic_CO)
         isp = isp_CO
      case (ic_Hp)
         isp = isp_HII
      case default
         write(*,'(a,i0,a)') ' (carrier_species_index) carrier index ',   &
            ic, ' names no species'
         error stop 1
      end select
      end function carrier_species_index

      ! Mass of a transported carrier [m_H], from the species table.
      double precision function carrier_mass_amu(ic) result(m)
      integer, intent(in) :: ic
      m = species_mass_amu(carrier_species_index(ic))
      end function carrier_mass_amu

      ! ------------------------------------------------------------- !

      ! Face gradient and drift coefficients of each carrier, and the
      ! Peclet-hybrid donor switch.  Faces 0 and N are left at zero, which
      ! IS the zero-diffusive-flux boundary of sec. 5.
      subroutine carrier_face_coefficients(ntot, TK, mbar, gphys, Dco,   &
                                           rp, Agrd, Bdrf, updrf)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: ntot, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: mbar, gphys
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: Dco
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp
      real(dp), dimension(0:N,n_carrier_max),       intent(out) :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(out) :: updrf
      real(dp), dimension(1-Ng:N+Ng) :: Gco
      real(dp) :: dr_f, ntf, Df, Gf, Kf
      integer  :: j, ic

      Agrd  = 0.0d0
      Bdrf  = 0.0d0
      updrf = 0

      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         ! Settling coefficient of this carrier, G = (m_i - m_bar) g/(k T).
         do j = 1-Ng, N+Ng
            Gco(j) = (carrier_mass_amu(ic)*mu - mbar(j))*gphys(j)   &
                     /(kb_erg*max(TK(j), 1.0d0))
         enddo
         do j = 1, N-1
            dr_f = max(rp(j+1) - rp(j), 1.0d0)
            ntf  = 0.5d0*(ntot(j) + ntot(j+1))
            Df   = 0.5d0*(Dco(j,ic) + Dco(j+1,ic))
            Gf   = 0.5d0*(Gco(j) + Gco(j+1))
            Kf   = 0.5d0*(kzz_cell(j) + kzz_cell(j+1))
            Agrd(j,ic) = ntf*(Df + Kf)/dr_f
            Bdrf(j,ic) = ntf*Df*Gf
            if (abs(Bdrf(j,ic))*dr_f .le. 2.0d0*ntf*(Df + Kf)) then
               updrf(j,ic) = 0        ! central: the drift is resolved
            else if (Bdrf(j,ic) .ge. 0.0d0) then
               updrf(j,ic) = -1       ! settling inward: donor is cell j+1
            else
               updrf(j,ic) = +1       ! rising outward: donor is cell j
            endif
         enddo
      enddo
      end subroutine carrier_face_coefficients

      ! ------------------------------------------------------------- !

      ! Rebuild the photolysis and Lyman-Werner rates from the CURRENT
      ! carrier columns and store them on the frozen cell state, so the
      ! water layer shields itself within the relaxation (sec. 4).  Uses the
      ! opa_pf weighting the last ioniz_eq left, like every other column.
      subroutine carrier_photolysis(rho, TK, f_sp)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho, TK
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), dimension(1-Ng:N+Ng) :: nH2, nH2O, nOH, nCO, nH_nuc, nHI
      real(dp), dimension(1-Ng:N+Ng) :: c1, c2, c3, c4, fsh, trl, klw
      real(dp), dimension(1-Ng:N+Ng) :: pls, pla, kco, thco
      real(dp), dimension(1-Ng:N+Ng,n_fuv_band) :: tau_b, jh2o, joh
      real(dp) :: amax, ovl
      integer  :: j
      nH2  = rho*n0*f_sp(:,isp_H2)
      nH2O = rho*n0*f_sp(:,isp_H2O)
      nOH  = rho*n0*f_sp(:,isp_OH)
      nCO  = rho*n0*f_sp(:,isp_CO)
      ! Neutral atomic hydrogen: the absorber of the Ly-alpha band, which is
      ! the H I resonance line (fuv_lw_photon_field).
      nHI  = rho*n0*f_sp(:,1)
      ! Total hydrogen NUCLEUS density, the density axis of the
      ! self-shielding table: every carrier weighted by the H nuclei it
      ! holds, the same sum ionization_equilibrium forms for nh.
      nH_nuc = rho*n0*(f_sp(:,1) + f_sp(:,2)                             &
                     + 2.0_dp*(f_sp(:,isp_H2) + f_sp(:,isp_H2p))         &
                     + 3.0_dp*f_sp(:,isp_H3p) + f_sp(:,isp_HeHp)         &
                     + f_sp(:,isp_OH) + 2.0_dp*f_sp(:,isp_H2O))
      call fuv_lw_photon_field(nH2, nH2O, nOH, nCO, nHI, TK, nH_nuc,     &
                               c1, c2, c3, c4, fsh, trl, tau_b, klw,      &
                               pls, pla, kco, thco,                       &
                               jh2o, joh, amax, ovl)
      do j = 1-Ng, N+Ng
         cph_klw(j)      = klw(j)
         cph_kco(j)      = kco(j)
         cph_jh2o(j,:)   = jh2o(j,:)
         cph_joh(j,:)    = joh(j,:)
      enddo
      end subroutine carrier_photolysis

      ! ------------------------------------------------------------- !

      ! The H nuclei cell j's transported carriers may hold [cm^-3].
      !
      ! It is NOT the hydrogen element total.  The H sitting in H II, H2+ and
      ! H3+ is frozen for the duration of this step (sec. 4), and the
      ! write-back does not re-solve those stages -- it rescales them in
      ! proportion to what they already held.  So a carrier state that used
      ! the whole element would leave the closure species nothing: the
      ! remainder handed to H I, H II, H2+ and H3+ goes to zero, the free
      ! electron density with it, and the very next ionization sweep is given
      ! a cell it cannot find a root in.  That is not hypothetical -- it is
      ! what the transported hot Uranus did in 80 of 500 cells between
      ! 1.032 and 1.156 R_p, ending in "persistent non-root chemical
      ! equilibrium" at r = 1.112 with n_e = 3.7 cm^-3.
      !
      ! This is the same budget carrier_source uses to close atomic H, and
      ! the two must not be able to drift apart, so both call this.
      double precision function hydrogen_available_to_carriers(j, nH_free) &
                                result(nH)
      ! The hydrogen nuclei the carriers of this cell may hold: the element
      ! total less the nuclei sitting in the stages this step holds FROZEN
      ! and rescales rather than re-solves.
      !
      ! WHICH STAGES THOSE ARE DEPENDS ON THE RUN, and that is the whole
      ! change the proton option makes here. Without it H+ is one of the
      ! frozen stages and its nuclei are unavailable; with it H+ is a
      ! carrier, solved by the same Newton as H2, so its nuclei are on the
      ! carrier side of the budget and are subtracted by the caller through
      ! nc(ic_Hp) instead of here. Leaving cbg_nhii in the subtraction with
      ! the proton carried would charge the proton twice -- once as frozen
      ! background and once as a carrier -- and cap the carriers below the
      ! hydrogen that exists.
      !
      ! The subtraction itself is the conditional hydrogen budget of the
      ! one stoichiometric map (element_inventory); this routine supplies
      ! the frozen background of the step it is asked about.
      integer,  intent(in) :: j
      real(dp), intent(in) :: nH_free
      nH = carrier_hydrogen_budget(nH_free, cbg_nhii(j), cbg_nh2p(j),    &
                                   cbg_nh3p(j), ionization_transport)
      end function hydrogen_available_to_carriers

      ! ------------------------------------------------------------- !

      ! THE REFERENCE DENSITY OF A CARRIER [cm^-3]: the free density of the
      ! element the carrier is made of, which is the largest density that
      ! carrier can physically reach in this cell.  It is the ONE table of
      ! the module: the Jacobian's finite-difference step and the residual
      ! row's absolute floor are both fractions of it, and both read it here
      ! so that they cannot state different elements for the same carrier.
      !
      ! WHICH ELEMENT EACH CARRIER BELONGS TO is stoichiometry, not a
      ! convention: H2 and H+ carry hydrogen nuclei only, OH and H2O are
      ! limited by the oxygen they contain (their hydrogen is abundant where
      ! their oxygen is not), and CO needs one nucleus of each of oxygen and
      ! carbon, so the smaller free density bounds it.
      !
      ! Charging the proton to the oxygen reservoir is not a small
      ! mis-scaling but a total loss of the row: in a hydrogen and helium
      ! atmosphere nO_free is identically zero, so both the step and the
      ! floor fall to their 1e-300 guard and the diagonal derivative of the
      ! proton row comes back as exactly zero
      ! (src/tests/carrier_reference_scales/).
      double precision function carrier_element_reference_density(ic,     &
                               nH_free, nO_free, nC_free) result(nref)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nH_free, nO_free, nC_free
      select case (ic)
      case (ic_H2, ic_Hp)
         nref = nH_free
      case (ic_OH, ic_H2O)
         nref = nO_free
      case (ic_CO)
         nref = min(nO_free, nC_free)
      case default
         write(*,'(a,i0,a)') ' (carrier_element_reference_density) '//    &
            'carrier index ', ic, ' belongs to no element'
         error stop 1
      end select
      end function carrier_element_reference_density

      ! ------------------------------------------------------------- !

      ! Density perturbation [cm^-3] with which the chemistry rows of one
      ! carrier are differentiated (the Jacobian block of solve_carriers).
      !
      ! THE STEP NEEDS A FLOOR TIED TO THE ELEMENT, and that was a real
      ! defect rather than a precaution.  Every source term is built from
      ! densities of order the element total, so a perturbation far below
      ! the round-off of THOSE terms returns noise for a derivative.  Where
      ! a carrier is empty -- OH high in the wind, which photolysis keeps at
      ! 1e-30 of its element -- a step proportional to the carrier itself is
      ! exactly that, and the Newton then churned to its iteration cap with
      ! a relative residual of 1 in that one cell (measured on HD 189733 b).
      ! The floor is 1e-12 of the reference density of the carrier's own
      ! element, which is resolvable against a double's 1e-16 while still
      ! being a small perturbation of the row.
      double precision function carrier_source_derivative_step(ic, nc_ic, &
                               nH_free, nO_free, nC_free) result(dn)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nc_ic, nH_free, nO_free, nC_free
      dn = max(1.0d-6*abs(nc_ic),                                        &
               1.0d-12*carrier_element_reference_density(ic, nH_free,    &
                                                         nO_free,        &
                                                         nC_free),       &
               1.0d-300)
      end function carrier_source_derivative_step

      ! ------------------------------------------------------------- !

      ! Chemical production minus loss of the four carriers in cell j, at the
      ! trial densities nc(1..4) [cm^-3].  The rows are the SAME rows the
      ! local equilibrium solve uses (sec. 4); nothing is rewritten here.
      subroutine carrier_source(j, nc, nH_free, nO_free, src)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: nc(n_carrier_max), nH_free, nO_free
      real(dp), intent(out) :: src(n_carrier_max)
      real(dp) :: fv(10)
      real(dp) :: n_hi, n_hii, n_h2p, n_h3p, n_hehp
      real(dp) :: n_hei, n_heii, n_heiii, n_heiTR, n_heiSI, n_e, n_o0

      ! THE PROTON IS EITHER FROZEN BACKGROUND OR THE TRIAL UNKNOWN.
      ! With the ionization state carried, n(H+) is what this Newton is
      ! solving for in this cell, so the row has to be evaluated at the
      ! trial value and not at the value the last sweep left behind --
      ! otherwise the recombination sink and the H closure below both refer
      ! to a different proton than the one the row is moving, and the
      ! Jacobian column for it is zero.
      if (ionization_transport) then
         n_hii = nc(ic_Hp)
      else
         n_hii = cbg_nhii(j)
      endif
      n_h2p   = cbg_nh2p(j)
      n_h3p   = cbg_nh3p(j)
      n_hehp  = cbg_nhehp(j)
      n_heii  = cbg_nheii(j)
      n_heiii = cbg_nheiii(j)
      n_heiTR = cbg_nheiTR(j)
      n_hei   = cbg_nhei(j)
      n_heiSI = max(n_hei - n_heiTR, 0.0d0)
      ! THE ELECTRONS THAT RECOMBINE WITH THE TRIAL PROTONS ARE THE TRIAL
      ! PROTONS' OWN. The rest of the electron budget -- He+, He++, the
      ! molecular ions, the metals -- is background and stays frozen with
      ! it, but holding the H+ electrons frozen while H+ itself moves would
      ! break charge neutrality inside the solve and, where H+ carries the
      ! electrons (the wind of this configuration), would evaluate the
      ! recombination sink a_hii n_e n_hii at an n_e belonging to a
      ! different ionization state. n_e is a sum over charges, so replacing
      ! one term by its trial value is exact.
      n_e     = cbg_ne(j)
      if (ionization_transport)                                              &
         n_e = max(cbg_ne(j) - cbg_nhii(j) + n_hii, 0.0d0)

      ! Atomic H closes the hydrogen budget: what the element has left after
      ! the ions, the molecular ions and the trial carriers.  This is the
      ! coupling that makes the H2 source a function of the unknown.
      ! nH_free is the element headroom, so the frozen ions come off here.
      ! With the proton carried it is one of the carriers this closure
      ! subtracts (one H nucleus each), and hydrogen_available_to_carriers
      ! has stopped subtracting it as frozen background.
      if (ionization_transport) then
         n_hi = max(hydrogen_available_to_carriers(j, nH_free)           &
                    - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O)     &
                    - nc(ic_Hp), 0.0d0)
      else
      n_hi = max(hydrogen_available_to_carriers(j, nH_free)              &
                 - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O), 0.0d0)
      endif
      ! Free atomic oxygen closes the oxygen budget the same way: the
      ! element total less its ionized stages and less the carriers.
      n_o0 = max(nO_free - cbg_nOion(j) - nc(ic_OH) - nc(ic_H2O)         &
                         - nc(ic_CO), 0.0d0)

      ! The cell state of the last sweep, with only the photolysis rates
      ! refreshed: n_ofam and n_co are carried over unchanged because the
      ! rows called below take the densities they need as arguments.
      ieq_cell       = bg_cell(j)
      ieq_cell%k_LW  = cph_klw(j)
      call set_mol_coeffs(ieq_cell%T_K, ieq_cell%ntot)
      call set_oxygen_coeffs(ieq_cell%T_K, cph_jh2o(j,:), cph_joh(j,:))

      fv = 0.0d0
      call mol_heh_rows(fv, n_hi, n_hii, nc(ic_H2), n_h2p, n_h3p,        &
                        n_hehp, n_heiSI, n_heiTR, n_heii, n_heiii, n_e,  &
                        ieq_cell%ntot,                                   &
                        ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,  &
                        ieq_cell%P_HeITR, ieq_cell%P_H2,                 &
                        ieq_cell%P_H2_di, ieq_cell%P_H2_dd,              &
                        ieq_cell%P_H2_nd, ieq_cell%k_LW,                 &
                        ieq_cell%rchiiB, ieq_cell%rcheiiB,               &
                        ieq_cell%rcheiiiB, ieq_cell%rcheiTR,             &
                        ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,           &
                        ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,       &
                        ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b,      &
                        ieq_cell%Q31, ieq_cell%A31)
      ! The photolysis channels set_oxygen_coeffs built for this cell: the
      ! transport operator evaluates the chemistry at the cell's own,
      ! unscaled radiation field.
      call oxygen_carrier_rows(fv, 9, n_hi, nc(ic_H2), nc(ic_OH),        &
                               nc(ic_H2O), n_o0, oj3, oj4, oj5, oj7)

      src(ic_H2)  = fv(4)
      src(ic_OH)  = fv(9)
      src(ic_H2O) = fv(10)
      ! ---- CO: the one-sided destruction model -----------------------
      ! Two channels and no formation term, which is what "one-sided"
      ! means (docs/b3b_co_destruction_design_20260906.md secs. 1 and 2):
      !
      !   D1  He+ + CO -> C+ + O + He, UMIST RATE22 entry 4068, measured,
      !       accuracy better than 25 percent, temperature independent at
      !       the Langevin value 1.6e-9 cm^3 s^-1 (oxygen_rates::
      !       rk_D1_Hep_CO).  n_heii is the frozen He+ of this step.
      !   D2  CO + hv -> C + O on the 912-1201 A beam, with the Visser,
      !       van Dishoeck & Black (2009) shielding function of the
      !       star-ward CO and H2 columns.  cph_kco is that rate, the cell
      !       mean, rebuilt from the CURRENT carrier columns at the top of
      !       the step, so the CO layer shields itself inside the
      !       relaxation.
      !
      ! Both are proportional to the trial nc(ic_CO), so the term is a
      ! diagonal contribution to the CO row and the Jacobian sees it.
      !
      ! WHERE THE NUCLEI GO.  The oxygen is already accounted for: n_o0
      ! above is the closure nO_free - ions - OH - H2O - CO, so an oxygen
      ! nucleus that leaves CO appears as free atomic O in the same
      ! evaluation of the same cell and the oxygen element total does not
      ! move.  The carbon returns to the carbon pool and carrier_write_back
      ! shares it over C I / C II / C III in proportion to what they held;
      ! the ionization sweep that follows re-solves that partition from its
      ! own balance.  The step moves nuclei, the sweep moves charge.
      !
      ! WHERE THE He+ IS CONSUMED.  D1 destroys a helium ion as well as a CO
      ! molecule, and the two operators split the one reaction between them:
      ! this row removes the CO at the frozen He+ of the step, and row (2)
      ! of mol_heh_rows removes the He+ at the frozen CO of the sweep, both
      ! from rk_D1_Hep_CO, so the reaction is complete over the pair and
      ! counted once in each.  What the split leaves is a lag of one
      ! operator, not a missing sink; the domain record counts the cells in
      ! which D1 is the leading He+ loss (carrier_co_domain_record), which
      ! is where that lag is worth reading.
      !
      ! WHERE THE MODEL IS OUT OF DOMAIN.  Below the helium ionization
      ! front, in the shielded molecular layer, both channels are slow
      ! against the residence time and the omitted formation is what would
      ! set the CO abundance.  The rates are still evaluated there -- they
      ! are physics, not a switch -- and the cells are counted; over the
      ! run's own physical time a destruction that slow cannot remove what
      ! an equally omitted formation cannot rebuild, so the transported
      ! value stands, which is the correct behaviour of a transport
      ! operator in a quenched layer.
      src(ic_CO)  = -(rk_D1_Hep_CO()*n_heii + cph_kco(j))*nc(ic_CO)
      ! Row (1) of mol_heh_rows is the complete proton balance and is taken
      ! whole: photoionization of H including the secondaries, the
      ! dissociative and double photoionization branches of H2, the H2+ + H
      ! and He+ + H2 channels and the He(2^3S) Penning term, less radiative
      ! recombination and less the H+ + H2 reactions R10 and R13. There is
      ! no proton chemistry written here that the local solve does not
      ! already have; the difference is only where the row is evaluated.
      src(ic_Hp)  = 0.0d0
      if (ionization_transport) src(ic_Hp) = fv(1)

      end subroutine carrier_source

      ! ------------------------------------------------------------- !

      ! Face flux of one carrier and its two derivatives.  The drift term is
      ! linear in f (not the X(1-X) form of the element operator, whose two
      ! factors exist because that variable is one half of a binary pair).
      subroutine carrier_face_flux(fl, fr, wcl, wcr, Agr, Bst, upw,      &
                                   Jf, dJl, dJr)
      ! fl, fr are the unknown -- the fraction per unit mass -- and wcl, wcr
      ! turn each into the PARTICLE mixing ratio X = wc*y whose gradient
      ! drives molecular diffusion and whose value the settling drift acts
      ! on.  The two are not the same variable wherever mbar varies.
      real(dp), intent(in)  :: fl, fr, wcl, wcr, Agr, Bst
      integer,  intent(in)  :: upw
      real(dp), intent(out) :: Jf, dJl, dJr
      real(dp) :: wl, wr
      if (upw .eq. 0) then
         wl = 0.5d0
         wr = 0.5d0
      else if (upw .lt. 0) then
         wl = 0.0d0
         wr = 1.0d0
      else
         wl = 1.0d0
         wr = 0.0d0
      endif
      Jf  = -Agr*(wcr*fr - wcl*fl) - Bst*(wl*wcl*fl + wr*wcr*fr)
      dJl =  (Agr - Bst*wl)*wcl
      dJr = (-Agr - Bst*wr)*wcr
      end subroutine carrier_face_flux

      ! ------------------------------------------------------------- !

      ! Residual of the backward-Euler transport-chemistry system, and the
      ! face-flux derivatives the Jacobian needs.  rnorm is the largest
      ! RELATIVE row imbalance, built the way composition_residual builds
      ! its own: the row divided by the sum of the magnitudes of its terms.
      subroutine carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp,   &
                                  rep, msum, Frho, advect, Agrd, Bdrf,   &
                                  updrf,                                 &
                                  nH_free, nO_free, nC_free, sigrate,    &
                                  res, Jf, dJl, dJr, rnorm, rnorm_phys)
      ! advect: whether the row carries the material advective term.  The
      ! MARCHING rows do not, because the advection of a carrier is the
      ! divergence of the hydrodynamic face mass flux and is taken with the
      ! mass row inside the Runge-Kutta stages; what is left here is the
      ! operator-split remainder, the diffusive and drift fluxes and the
      ! chemistry.  The STATIONARY evaluation does, because the balance a
      ! steady state satisfies is the whole equation and not the remainder
      ! of an operator split; so does the fixed-wind relaxation, which has
      ! no stages to ride on.  Where it is carried it is the SAME face-flux
      ! divergence the stages take (carrier_advective_divergence).
      logical, intent(in) :: advect
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp, rep
      ! Mass of the mixture per unit of f_sp, and the face mass flux of the
      ! mass row of this state: the two the advective term is built from.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: msum, Frho
      real(dp), dimension(0:N,n_carrier_max),       intent(in)  :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(in)  :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nC_free
      ! Signal-crossing rate (|v| + c_s)/dr of each cell [1/s].  It carries
      ! the absolute floor of the row scale below, which used to be carried
      ! by 1/dt -- see there.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: sigrate
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: res
      real(dp), dimension(0:N,n_carrier_max),       intent(out) :: Jf, dJl, dJr
      real(dp),                                 intent(out) :: rnorm
      ! THE SAME ROWS ON THE SCALE THE RETURNED STATE IS JUDGED ON, the
      ! physical terms alone.  The Newton needs both: the full-term measure
      ! is what its iteration is conditioned on, and this one is what its
      ! result has to satisfy, so stopping on the first alone leaves, at a
      ! short substep, a state the acceptance refuses
      ! (declaration of row_terms_phys).
      real(dp),                       optional, intent(out) :: rnorm_phys
      ! Cell and species carrying the worst relative row imbalance, for the
      ! debug line: a Newton that stops short says WHERE it stopped.
      real(dp) :: nc(n_carrier_max), src(n_carrier_max)
      real(dp) :: Kj, sL, sR, dsc, dph, rr, nfl
      real(dp) :: rnormp
      integer  :: j, ic
      ! Row scale of every row of this assembly, kept so that the worst
      ! relative imbalance can be located AFTER the cell loop rather than
      ! inside it -- see there.
      real(dp) :: rsc(1:N,n_carrier_max)
      ! The advective term of every row and its size, formed once from the
      ! face fluxes before the cell loop: the face reconstruction is an
      ! operation on the whole column and not on a cell.
      real(dp) :: adv(1:N,n_carrier_max), advmag(1:N,n_carrier_max)

      adv    = 0.0d0
      advmag = 0.0d0
      if (advect) call carrier_advective_divergence(fc, msum, Frho,      &
                                                    adv, advmag)

      res = 0.0d0
      Jf  = 0.0d0
      dJl = 0.0d0
      dJr = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N-1
            call carrier_face_flux(fc(j,ic), fc(j+1,ic), wfac(j),        &
                                   wfac(j+1), Agrd(j,ic),                &
                                   Bdrf(j,ic), updrf(j,ic),              &
                                   Jf(j,ic), dJl(j,ic), dJr(j,ic))
         enddo
      enddo

      ! CELL-LOCAL, AND THE ONE THING THAT WAS NOT.  Every term of a carrier
      ! row is built from this cell's own state and its two stored face
      ! fluxes, and carrier_source evaluates the chemistry at this cell
      ! through module state that is threadprivate (ieq_cell, the mk* rate
      ! coefficients of System_HeH_mol and the oj* photolysis channels), so
      ! each thread fills its own copy per cell.  The exception was the
      ! running maximum: rnorm together with the cell and carrier that
      ! attain it, which is a max and an argmax and cannot be taken inside a
      ! parallel loop without making the diagnostic depend on which thread
      ! got there first.
      !
      ! So the loop takes the maximum only -- a max is exact and independent
      ! of the order it is taken in, so rnorm is bitwise the sequential value
      ! at any number of threads -- and the cell that attains it is found by
      ! a scan afterwards, which is what makes pct_worst_j / pct_worst_ic
      ! deterministic.  The sequential loop updated on a STRICT improvement,
      ! so its winner was the FIRST row in (j, then ic) order attaining the
      ! maximum, and the scan below stops at exactly that row.  With no row
      ! above the initial zero it leaves both untouched, as the sequential
      ! loop did.
      rnorm  = 0.0d0
      rnormp = 0.0d0
      rsc    = 0.0d0
      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,ic,nc,src,Kj,sL,sR,dsc,dph,rr,nfl)               &
      !$omp   reduction(max:rnorm,rnormp)
      do j = 1, N
         nc = fc(j,:)*nrho(j)
         call carrier_source(j, nc, nH_free(j), nO_free(j), src)
         Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
         sL = rep(j-1)**2
         sR = rep(j)**2
         do ic = 1, n_carrier
            ! A carrier this run does not solve keeps a zero row and an
            ! identity block, so its unknown never moves and it contributes
            ! nothing to the residual norm.
            if (.not. carrier_solved(ic)) then
               res(j,ic) = 0.0d0
               rsc(j,ic) = 0.0d0
               if (allocated(row_terms)) row_terms(j,ic) = 0.0d0
               if (allocated(row_terms_phys))                              &
                  row_terms_phys(j,ic) = 0.0d0
               cycle
            endif
            rr  = nrho(j)*(fc(j,ic) - fc_old(j,ic))/dt_phys(j)
            dsc = abs(nrho(j)*fc(j,ic)/dt_phys(j))                       &
                + abs(nrho(j)*fc_old(j,ic)/dt_phys(j))
            ! dph accumulates the same terms EXCEPT the time term, and it
            ! is accumulated beside dsc rather than differenced out of it,
            ! so the Newton's own scale keeps its summation order exactly.
            dph = 0.0d0
            rr  = rr + Kj*(sR*Jf(j,ic) - sL*Jf(j-1,ic))
            dsc = dsc + Kj*(sR*abs(Jf(j,ic)) + sL*abs(Jf(j-1,ic)))
            dph = dph + Kj*(sR*abs(Jf(j,ic)) + sL*abs(Jf(j-1,ic)))
            ! THE MATERIAL ADVECTION IS THE DIVERGENCE OF THE MASS ROW'S
            ! OWN FACE FLUXES, the same expression the Runge-Kutta stages
            ! subtract from rho Y_c, so that a marched state and a Newton
            ! state are states of one equation.  The base face needs no
            ! special case: the face composition is reconstructed from the
            ! ghosts wherever the face mass flux flows inward, which is the
            ! handoff partition where a handoff states one and the base
            ! cell's own where none does.
            if (advect) then
               rr  = rr + adv(j,ic)
               dsc = dsc + advmag(j,ic)
               dph = dph + advmag(j,ic)
            endif
            rr  = rr - src(ic)
            dsc = dsc + abs(src(ic))
            dph = dph + abs(src(ic))
            ! ABSOLUTE FLOOR ON THE ROW SCALE, and it is not cosmetic. A
            ! relative residual with no floor asks a row whose species is
            ! 1e-28 of its element to balance to 1e-12 of ITSELF, and out in
            ! the wind every term of such a row is round-off: measured on
            ! the HD 209458 b example, H2O jumped four decades from cell to
            ! cell at 1e-28 and the Newton spent its whole iteration budget
            ! there at a relative residual of 1 while every cell that
            ! carries a molecule was already at 1e-13. The floor is 1e-20 of
            ! the element the carrier belongs to, converted to the row's own
            ! volumetric-rate units: a density that small cannot change any
            ! observable, so a row below it IS converged.
            !
            ! IT IS CONVERTED BY THE LARGER OF TWO RATES, and the second one
            ! is why this line changed. Dividing by the step alone made the
            ! floor vanish in the STEADY evaluation, which substitutes
            ! dt = 1e30 to kill the time term (carrier_steady_residual): the
            ! floor went with it, and the far wind was then measured with no
            ! floor at all. The cell's signal-crossing rate is the same
            ! statement written without a step, so it survives dt -> infinity
            ! and agrees with the step form whenever the step is the cell's
            ! own crossing time, which is what the relaxation uses.
            !
            ! WHICH ELEMENT the carrier belongs to is stated in exactly one
            ! place, the same one the Jacobian step reads: a floor and a
            ! differentiation step that named different elements for one
            ! carrier would measure the row against a reservoir it cannot
            ! draw on.
            nfl = carrier_element_reference_density(ic, nH_free(j),      &
                                                    nO_free(j),         &
                                                    nC_free(j))
            dsc = dsc + 1.0d-20*nfl*max(1.0d0/dt_phys(j), sigrate(j))
            ! THE PHYSICAL SCALE TAKES THE SIGNAL RATE ALONE.  The step form
            ! of the same floor grows without bound as the step shrinks, and
            ! a scale that grows like 1/dt is exactly what a statement about
            ! the returned state must not have; the cell's signal-crossing
            ! rate is the same floor written without a step (declaration of
            ! row_terms_phys).
            dph = dph + 1.0d-20*nfl*sigrate(j)
            res(j,ic)  = rr
            rsc(j,ic)  = dsc
            if (allocated(row_terms)) row_terms(j,ic) = dsc
            if (allocated(row_terms_phys)) row_terms_phys(j,ic) = dph
            if (abs(rr)/max(dsc, 1.0d-300) .gt. rnorm)                   &
               rnorm = abs(rr)/max(dsc, 1.0d-300)
            if (abs(rr)/max(dph, 1.0d-300) .gt. rnormp)                  &
               rnormp = abs(rr)/max(dph, 1.0d-300)
         enddo
      enddo
      !$omp end parallel do

      if (rnorm .gt. 0.0d0) then
         outer: do j = 1, N
            do ic = 1, n_carrier
               if (.not. carrier_solved(ic)) cycle
               if (abs(res(j,ic))/max(rsc(j,ic), 1.0d-300)               &
                   .eq. rnorm) then
                  pct_worst_j  = j
                  pct_worst_ic = ic
                  exit outer
               endif
            enddo
         enddo outer
      endif
      if (present(rnorm_phys)) rnorm_phys = rnormp
      end subroutine carrier_residual

      ! ------------------------------------------------------------- !

      ! Newton solve of the coupled system, block-tridiagonal in space with
      ! 4x4 blocks.  The chemistry Jacobian is a forward difference of the
      ! same rows the residual calls, so a change to the network reaches the
      ! Jacobian without a second edit.
      subroutine solve_carriers(fc, fc_old, nrho, wfac, dt_phys, rp,     &
                                rep, msum, Frho, Agrd, Bdrf, updrf, TK,  &
                                sigrate,                                 &
                                nH_free, nO_free, nC_free, verdict)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)    :: fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rp, rep
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: msum
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: Frho, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: sigrate
      real(dp), dimension(0:N,n_carrier_max),       intent(in)    :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(in)    :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      ! THE VERDICT ON THE STATE THIS ROUTINE RETURNS, taken on that state's
      ! own rows.  It carries how the iteration ended (one of the
      ! carrier_solve_* constants) and how it stalled, but neither of those
      ! decides it: a state that satisfies no equation puts a solution of
      ! nothing into the chemical heat, the free atomic H the ionization
      ! sweep closes on and the H2 opacity, whatever the iteration reported
      ! about itself.  The caller decides what to do with a refusal.
      type(carrier_verdict),                    intent(out)   :: verdict

      real(dp), dimension(1:N,n_carrier_max)  :: res, rhs
      ! Which rows of the returned state the verdict found at the round-off
      ! of their own full terms, so the refused rows can be reported with
      ! the two scales they were measured against.
      logical,  dimension(1:N,n_carrier_max)  :: row_roundoff
      real(dp), dimension(0:N,n_carrier_max)  :: Jf, dJl, dJr
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max) :: aa, bb, cc
      real(dp), dimension(1:N,n_carrier_max)  :: dfc
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: ftry
      real(dp) :: nc(n_carrier_max), s0(n_carrier_max), s1(n_carrier_max)
      real(dp) :: Kj, sL, sR, cadv, dn
      ! The two face coefficients of the advective term, frozen with the
      ! wind over this solve.
      real(dp), dimension(1:N) :: advj, advm
      real(dp) :: rnorm, rprev, rstart, rtry, damp
      ! The same residual on the scale the returned state is judged on.
      real(dp) :: rnormp, rtryp
      integer  :: j, ic, kc, it, ihalf
      ! Did the last Newton iteration exhaust its halvings without finding a
      ! trial state below the current residual? That is a line-search
      ! failure and is a different statement from running out of iterations:
      ! the direction was not a descent direction of the merit the solve
      ! measures.
      logical  :: line_search_failed
      ! How the last iteration stalled, one of the carrier_stall_*
      ! constants.  Kept after the loop so the termination can be named.
      integer  :: stall
      integer  :: status

      ! The row scales the verdict below is taken on are the ones this
      ! assembly writes, so the array has to exist before the first call.
      if (.not. allocated(row_terms))                                     &
         allocate(row_terms(1:N,n_carrier_max))
      if (.not. allocated(row_terms_phys))                                &
         allocate(row_terms_phys(1:N,n_carrier_max))
      if (.not. allocated(pct_cell_constrained))                          &
         allocate(pct_cell_constrained(1:N+Ng))
      advj = 0.0d0
      advm = 0.0d0
      if (carrier_rows_advect)                                            &
         call carrier_advective_face_coefficients(Frho, advj, advm)
      call carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp, rep,    &
                            msum, Frho, carrier_rows_advect, Agrd, Bdrf, &
                            updrf, nH_free, nO_free,                     &
                            nC_free, sigrate, res, Jf, dJl, dJr, rnorm,  &
                            rnormp)
      rstart = rnorm
      rprev  = rnorm
      pct_newton_steps = 0
      line_search_failed = .false.
      stall              = carrier_stall_none

      do it = 1, newton_maxit
         ! WHERE THE ITERATION MAY STOP, and it is the condition the
         ! RETURNED STATE has to satisfy and not a second, unrelated one.
         ! The acceptance measures the row against its physical terms
         ! (carrier_returned_state_verdict), so an iteration that stops at
         ! a fixed fraction of the FULL terms stops, at a substep far below
         ! the physical time scale of the row, at a state that is 1/dt
         ! times worse on the scale it is about to be judged on: the
         ! stopping test and the acceptance would measure different things.
         ! So both have to hold -- the full-term accuracy this solve has
         ! always required, and the acceptance measure itself.
         !
         ! And the iteration may also stop because there is nothing left to
         ! reduce: a residual at the round-off of the row terms is what an
         ! exact solve leaves (carrier_row_roundoff).  That is not
         ! convergence, it is named as its own outcome below, and the state
         ! is put to the acceptance exactly like any other.
         if (rnorm .le. newton_tol .and. rnormp .le. carrier_accept_tol)  &
            exit
         if (rnorm .le. carrier_row_roundoff) exit

         ! ---- assemble the block tridiagonal
         aa = 0.0d0
         bb = 0.0d0
         cc = 0.0d0
         ! CELL-LOCAL ASSEMBLY.  Row j of the block-tridiagonal system is
         ! built from cell j's own geometry, its two stored face-flux
         ! derivatives and n_carrier + 1 evaluations of carrier_source at
         ! cell j, and it writes aa/bb/cc at index j alone.  carrier_source
         ! reaches the chemistry through threadprivate module state, so each
         ! thread carries its own rate coefficients.  There is no reduction
         ! here, so the assembled blocks are bitwise the sequential ones at
         ! any number of threads.  The Thomas sweep that consumes them is
         ! sequential in space and stays serial.
         !$omp parallel do default(shared) schedule(static)              &
         !$omp   private(j,ic,kc,Kj,sL,sR,cadv,nc,s0,s1,dn)
         do j = 1, N
            Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
            sL = rep(j-1)**2
            sR = rep(j)**2
            do ic = 1, n_carrier
               bb(j,ic,ic) = nrho(j)/dt_phys(j)
               ! The identity row of a carrier this run does not solve: the
               ! diagonal above with a zero right-hand side gives it a zero
               ! update, which is what keeps it out of the system without
               ! taking it out of the block.
               if (.not. carrier_solved(ic)) cycle
               if (j .lt. N) then
                  bb(j,ic,ic) = bb(j,ic,ic) + Kj*sR*dJl(j,ic)
                  cc(j,ic,ic) = cc(j,ic,ic) + Kj*sR*dJr(j,ic)
               endif
               if (j .gt. 1) then
                  bb(j,ic,ic) = bb(j,ic,ic) - Kj*sL*dJr(j-1,ic)
                  aa(j,ic,ic) = aa(j,ic,ic) - Kj*sL*dJl(j-1,ic)
               endif
               ! THE ADVECTIVE ENTRIES ARE THE DONOR-CELL LINEARIZATION OF
               ! THE FACE-FLUX DIVERGENCE the row now carries.  With the
               ! face composition taken from the cell the face mass flux
               ! selects, the term of cell j is
               !
               !   adv(j) = [ A_+ F_rho(j) Y(don(j))
               !              - A_- F_rho(j-1) Y(don(j-1)) ] / dV_j
               !            x msum(j)/m_c ,   Y = m_c f_c/msum ,
               !
               ! so d adv(j)/d f_c(k) is the same coefficient divided by the
               ! mixture mass of the DONOR cell.  The limiter of
               ! species_face_fraction and the second-order part of the
               ! reconstruction are left out of the Jacobian, as the deferred
               ! correction of the cell-velocity form was: they are the
               ! strongly nonlinear part of the operator and they reach
               ! beyond the tridiagonal block.  The residual is the full
               ! operator either way, so this changes what the Newton
               ! converges AT, not what it converges TO.
               !
               ! The entries exist only where the rows carry the advective
               ! term, which is the fixed-wind relaxation; in the marching
               ! loop that transport is the mass row's and is not
               ! differentiated here.
               if (carrier_rows_advect) then
                  cadv = advj(j)*msum(j)/carrier_mass_amu(ic)
                  if (Frho(j) .ge. 0.0d0) then
                     bb(j,ic,ic) = bb(j,ic,ic) + cadv/msum(j)
                  else if (j .lt. N) then
                     cc(j,ic,ic) = cc(j,ic,ic) + cadv/msum(j+1)
                  endif
                  cadv = advm(j)*msum(j)/carrier_mass_amu(ic)
                  if (j .eq. 1) then
                     ! The inner ghost is data, not an unknown: where the
                     ! base face flows inward the term it carries is a
                     ! constant of the row.  Where it flows outward the
                     ! donor is cell 1 itself.
                     if (Frho(0) .lt. 0.0d0)                             &
                        bb(1,ic,ic) = bb(1,ic,ic) + cadv/msum(1)
                  else if (Frho(j-1) .ge. 0.0d0) then
                     aa(j,ic,ic) = aa(j,ic,ic) + cadv/msum(j-1)
                  else
                     bb(j,ic,ic) = bb(j,ic,ic) + cadv/msum(j)
                  endif
               endif
            enddo
            ! ---- chemistry block, by forward difference of the same rows
            !
            ! The perturbation is carrier_source_derivative_step: a fraction
            ! of the carrier itself with a floor tied to the free density of
            ! the carrier's OWN element, which is the same table the
            ! residual row's absolute floor reads.
            nc = fc(j,:)*nrho(j)
            call carrier_source(j, nc, nH_free(j), nO_free(j), s0)
            do kc = 1, n_carrier
               if (.not. carrier_solved(kc)) cycle
               dn = carrier_source_derivative_step(kc, nc(kc),           &
                                                   nH_free(j),           &
                                                   nO_free(j),           &
                                                   nC_free(j))
               nc(kc) = nc(kc) + dn
               call carrier_source(j, nc, nH_free(j), nO_free(j), s1)
               nc(kc) = nc(kc) - dn
               do ic = 1, n_carrier
                  if (.not. carrier_solved(ic)) cycle
                  bb(j,ic,kc) = bb(j,ic,kc)                              &
                              - (s1(ic) - s0(ic))/dn*nrho(j)
               enddo
            enddo
         enddo
         !$omp end parallel do

         rhs = -res
         call block_thomas(aa, bb, cc, rhs, dfc)

         ! ---- damped update, with f >= 0 enforced on every try
         damp = 1.0d0
         line_search_failed = .false.
         do ihalf = 0, newton_halves
            ftry = fc
            do ic = 1, n_carrier
               if (.not. carrier_solved(ic)) cycle
               do j = 1, N
                  ftry(j,ic) = max(fc(j,ic) + damp*dfc(j,ic), 0.0d0)
               enddo
               ! Only the OUTER ghosts mirror the interior.  The lower
               ! ghost is the inflow reservoir, not a copy of the base
               ! cell: it carries the composition the handoff states, and
               ! the advective term reconstructs the base face from it.
               ! Writing the base cell into it made the base condition
               ! degenerate into a zero-gradient one from the second
               ! relaxation pass on -- the marching loop hid that because
               ! its ionization sweep re-pins the ghost every step, and
               ! relax_photochemical_composition, which takes many
               ! transport steps between two sweeps, did not (sec. 158).
               ftry(N+1:N+Ng,ic) = ftry(N,ic)
            enddo
            call carrier_residual(ftry, fc_old, nrho, wfac, dt_phys, rp, &
                                  rep, msum, Frho, carrier_rows_advect,  &
                                  Agrd, Bdrf, updrf, nH_free,            &
                                  nO_free, nC_free, sigrate, res, Jf,    &
                                  dJl, dJr, rtry, rtryp)
            if (rtry .lt. rnorm) exit
            if (ihalf .eq. newton_halves) then
               ! Every halving down to 2^-8 of the Newton step left the
               ! residual where it was or above it: the direction is not a
               ! descent direction of the measure this solve minimizes.  The
               ! last trial is still taken -- rejection and retry are A3
               ! increment (b) and do not exist -- so the outcome has to be
               ! carried out in the status instead.
               line_search_failed = .true.
               exit
            endif
            damp = 0.5d0*damp
         enddo
         fc     = ftry
         rprev  = rnorm
         rnorm  = rtry
         rnormp = rtryp
         pct_newton_steps = it
         ! The iteration has stopped making progress. Whether that is the
         ! discretization's floor or a residual that merely fell a long way
         ! from a large start is carrier_stall_class's statement, and it is
         ! not an acceptance either way: there is nothing further this
         ! direction will do, so the loop ends and the state is judged.
         stall = carrier_stall_class(it, rnorm, rprev, rstart)
         if (stall .eq. carrier_stall_at_floor .or.                      &
             stall .eq. carrier_stall_at_relative_drop) exit
      enddo

      ! HOW THE NEWTON ITSELF ENDED, decided on the Newton's own residual
      ! and before the constraint below can move it.  The two questions are
      ! separate: this one asks whether the iteration solved the equation it
      ! was given, the verdict below asks what the state handed back
      ! satisfies.  A stall at the relative drop is named as what it is and
      ! is not counted as convergence.
      if ((rnorm .le. newton_tol .and. rnormp .le. carrier_accept_tol)    &
          .or. stall .eq. carrier_stall_at_floor) then
         status = carrier_solve_converged
      else if (rnorm .le. carrier_row_roundoff) then
         ! The residual is round-off of the row terms and the state is
         ! still above the acceptance: nothing further is reachable at this
         ! substep, and halving raises the floor rather than lowering it.
         status = carrier_solve_unreachable_at_substep
      else if (stall .eq. carrier_stall_at_relative_drop) then
         status = carrier_solve_stagnated
      else if (line_search_failed) then
         status = carrier_solve_line_search_failed
      else
         status = carrier_solve_iteration_cap
      endif

      pct_newton_resid_before_limit = rnorm

      call limit_to_element_budget(fc, nrho, TK, nH_free, nO_free,        &
                                   nC_free)

      ! THE RESIDUAL REPORTED IS THE RESIDUAL OF THE STATE RETURNED.  The
      ! Newton's last trial is not that state: limit_to_element_budget then
      ! rescales whole carrier families in any cell that asked for more
      ! nuclei than its element holds, and the rows of such a cell are not
      ! the rows the Newton last measured.  A number taken before the
      ! constraint would report a clamped state as converged, which is the
      ! one case in which the diagnostic is needed at all.  It costs one
      ! further residual assembly per solve: O(N) with one chemistry
      ! evaluation of each cell, against the several assemblies of each
      ! Newton iteration.  MEASURED on the hp_zero_seed regression case,
      ! 100 steps single-threaded, 5.78 s against 5.70 s, and the returned
      ! state is bitwise unchanged by it.
      call carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp, rep,    &
                            msum, Frho, carrier_rows_advect, Agrd,       &
                            Bdrf, updrf, nH_free, nO_free,               &
                            nC_free, sigrate, res, Jf, dJl, dJr, rnorm,  &
                            rnormp)
      pct_newton_resid = rnorm

      ! AND THE STATE IS JUDGED ON THOSE ROWS.  res is the residual of
      ! exactly the state fc that leaves here, after the constraint, and
      ! pct_cell_constrained says which cells the constraint moved.
      !
      ! THE SCALE IS row_terms_phys AND NOT row_terms.  The judgement is
      ! about the state, so it is taken against the size of the terms the
      ! row physically balances, the transport and the reaction.  row_terms
      ! carries the time term as well, which is 1/dt times a density: it
      ! belongs to the iteration the Newton just ran and not to the state,
      ! and measuring against it makes a fixed physical imbalance look
      ! smaller in proportion as the substep is shortened, so the retry
      ! controller could certify an unchanged defect by halving the step
      ! (declaration of row_terms_phys, with the measurement).
      ! The FULL terms go with it so that a row whose residual is already
      ! the round-off of the arithmetic that built it can be told from a
      ! row that is genuinely out of balance.  The first is the limit of
      ! double precision at this substep and is accepted; only the second
      ! refuses, and only the second is worth a shorter substep.
      call carrier_returned_state_verdict(res, row_terms_phys(1:N,:),     &
                                          pct_cell_constrained(1:N),      &
                                          status, stall, verdict,         &
                                          terms_full=row_terms(1:N,:),    &
                                          row_roundoff_limited=row_roundoff)
      call carrier_record_refused_rows(res, row_terms_phys(1:N,:),        &
                                       row_terms(1:N,:), row_roundoff,    &
                                       pct_cell_constrained(1:N))
      pct_verdict = verdict

      end subroutine solve_carriers

      ! ------------------------------------------------------------- !

      ! Block-tridiagonal Thomas sweep with dense n_carrier x n_carrier
      ! blocks.  Written here because the element operator's scalar Thomas
      ! sweep does not generalize: the species system is diagonal in space
      ! but dense in species through the chemistry.
      subroutine block_thomas(aa, bb, cc, dd, xx)
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max), intent(in)  :: aa,bb,cc
      real(dp), dimension(1:N,n_carrier_max),           intent(in)  :: dd
      real(dp), dimension(1:N,n_carrier_max),           intent(out) :: xx
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max) :: cp
      real(dp), dimension(1:N,n_carrier_max)           :: dp_
      real(dp) :: mm(n_carrier_max,n_carrier_max)
      real(dp) :: rhs(n_carrier_max,n_carrier_max+1)
      real(dp) :: piv, fac
      integer  :: j, i, k, l, ip

      cp  = 0.0d0
      dp_ = 0.0d0
      xx  = 0.0d0
      do j = 1, N
         ! mm = bb(j) - aa(j) cp(j-1)
         mm(1:n_carrier,1:n_carrier) = bb(j,1:n_carrier,1:n_carrier)
         if (j .gt. 1) then
            do i = 1, n_carrier
               do k = 1, n_carrier
                  do l = 1, n_carrier
                     mm(i,k) = mm(i,k) - aa(j,i,l)*cp(j-1,l,k)
                  enddo
               enddo
            enddo
         endif
         ! Solve mm * [cp(j) | dp(j)] = [cc(j) | dd(j) - aa(j) dp(j-1)]
         rhs(1:n_carrier,1:n_carrier) = cc(j,1:n_carrier,1:n_carrier)
         rhs(1:n_carrier,n_carrier+1) = dd(j,1:n_carrier)
         if (j .gt. 1) then
            do i = 1, n_carrier
               do l = 1, n_carrier
                  rhs(i,n_carrier+1) = rhs(i,n_carrier+1)                &
                                     - aa(j,i,l)*dp_(j-1,l)
               enddo
            enddo
         endif
         ! Gaussian elimination with partial pivoting on the small block.
         do i = 1, n_carrier
            ip  = i
            piv = abs(mm(i,i))
            do k = i+1, n_carrier
               if (abs(mm(k,i)) .gt. piv) then
                  piv = abs(mm(k,i))
                  ip  = k
               endif
            enddo
            if (ip .ne. i) then
               do k = 1, n_carrier
                  fac      = mm(i,k)
                  mm(i,k)  = mm(ip,k)
                  mm(ip,k) = fac
               enddo
               do k = 1, n_carrier+1
                  fac       = rhs(i,k)
                  rhs(i,k)  = rhs(ip,k)
                  rhs(ip,k) = fac
               enddo
            endif
            if (mm(i,i) .eq. 0.0d0) mm(i,i) = 1.0d-300
            do k = i+1, n_carrier
               fac = mm(k,i)/mm(i,i)
               if (fac .eq. 0.0d0) cycle
               do l = i, n_carrier
                  mm(k,l) = mm(k,l) - fac*mm(i,l)
               enddo
               do l = 1, n_carrier+1
                  rhs(k,l) = rhs(k,l) - fac*rhs(i,l)
               enddo
            enddo
         enddo
         do l = 1, n_carrier+1
            do i = n_carrier, 1, -1
               fac = rhs(i,l)
               do k = i+1, n_carrier
                  fac = fac - mm(i,k)*rhs(k,l)
               enddo
               rhs(i,l) = fac/mm(i,i)
            enddo
         enddo
         cp(j,1:n_carrier,1:n_carrier) = rhs(1:n_carrier,1:n_carrier)
         dp_(j,1:n_carrier)  = rhs(1:n_carrier,n_carrier+1)
      enddo

      xx(N,:) = dp_(N,:)
      do j = N-1, 1, -1
         do i = 1, n_carrier
            xx(j,i) = dp_(j,i)
            do k = 1, n_carrier
               xx(j,i) = xx(j,i) - cp(j,i,k)*xx(j+1,k)
            enddo
         enddo
      enddo
      end subroutine block_thomas

      ! ------------------------------------------------------------- !

      ! Hold the carriers inside the element simplex: the H nuclei they
      ! carry cannot exceed what the cell has left for them, and the same for
      ! O and C.  Every limiter is counted and the worst overshoot is
      ! reported (design sec. 3.6: no silent clamp).
      subroutine limit_to_element_budget(fc, nrho, TK, nH_free, nO_free,  &
                                         nC_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp) :: nc(n_carrier_max), got, sc, over, nH_car
      logical  :: hit
      integer  :: j, ic
      ! A carrier that holds the WHOLE of its element sits exactly on the
      ! face of the simplex -- all the carbon is CO in a molecular base --
      ! so a limiter that counted every round-off would report half the grid
      ! and mean nothing. Only an overshoot above this relative size is
      ! counted; the clamp itself is applied whatever the size.
      real(dp), parameter :: limit_report = 1.0d-10

      pct_limited     = 0
      pct_worst_limit = 0.0d0
      if (.not. allocated(pct_cell_constrained))                          &
         allocate(pct_cell_constrained(1:N+Ng))
      pct_cell_constrained = .false.
      ! The OUTER ghosts take the partition of the cell they mirror -- the
      ! zero-gradient outflow boundary -- and then go through the same
      ! limiter as an interior cell.  They HAVE to: a ghost carries its own
      ! density, so the same fraction is a different number of nuclei there
      ! and can ask for more than the ghost's element holds.  Leaving them
      ! out of the limiter was measured, and it broke the oxygen and carbon
      ! budgets by 2% in a 400-step run while every interior cell closed at
      ! round-off.  The LOWER ghosts are not mirrored and not limited: this
      ! operator does not own them (see solve_carriers).
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         fc(N+1:N+Ng,ic) = fc(N,ic)
      enddo
      do j = 1, N+Ng
         nc  = max(fc(j,:)*nrho(j), 0.0d0)
         hit = .false.
         ! ---- CO IS DESTROYED BY RATES, NOT BY A THERMODYNAMIC BOUND ---
         ! The only constraints applied here are conservation: the carriers
         ! of a cell cannot hold more nuclei than the cell has.  What sets
         ! the CO abundance is the CO row itself, which carries
         ! He+ + CO -> C+ + O + He (UMIST RATE22 4068) and the
         ! photodissociation of the 912-1201 A beam with the Visser et al.
         ! (2009) shielding function; the domain in which a
         ! destruction-only row is legitimate is measured cell by cell and
         ! reported (carrier_co_domain_take).
         ! carbon: CO is the only carbon carrier
         if (nc(ic_CO) .gt. nC_free(j)) then
            over = nc(ic_CO)/max(nC_free(j), 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            nc(ic_CO) = nC_free(j)
            pct_cell_constrained(j) = .true.
            if (over .gt. limit_report) hit = .true.
         endif
         ! oxygen: CO already fits, scale the water family into what is left
         got = nc(ic_OH) + nc(ic_H2O) + nc(ic_CO)
         if (got .gt. nO_free(j)) then
            over = got/max(nO_free(j), 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            sc = max(nO_free(j) - nc(ic_CO), 0.0d0)                      &
                 /max(nc(ic_OH) + nc(ic_H2O), 1.0d-300)
            nc(ic_OH)  = nc(ic_OH) *sc
            nc(ic_H2O) = nc(ic_H2O)*sc
            pct_cell_constrained(j) = .true.
            if (over .gt. limit_report) hit = .true.
         endif
         ! hydrogen: scale the three H-bearing carriers together, against
         ! the SAME budget the chemistry rows close atomic H with
         ! (hydrogen_available_to_carriers) -- not against the element
         ! total, which would let the carriers take the nuclei the frozen
         ! ion stages are holding and leave the write-back nothing to give
         ! back to them.
         ! The proton holds ONE H nucleus and, where it is carried, it is
         ! in this sum with the others. It is bounded the same way -- the
         ! whole H-bearing carrier set scaled by one factor when it asks for
         ! more hydrogen than the cell has -- rather than clamped on its
         ! own: a clamp on one carrier is a kink in the residual at the
         ! surface where it binds, and section 162 measured what a kink at
         ! a front does to a Newton. Scaling the set keeps the state inside
         ! the simplex along a straight line through the interior, which is
         ! the same construction the positivity limiter of the
         ! reconstruction uses.
         !
         ! WHERE IT BINDS IS THE OTHER END OF THE WIND FROM THE H2 FRONT.
         ! Below ~1.5 r_base the hydrogen is nearly all neutral and the
         ! proton asks for a small share; the budget is tight where
         ! x(H+) -> 1, high in the wind, and there H2 has long since gone,
         ! so the two carriers are rarely tight at the same place.
         nH_car = hydrogen_available_to_carriers(j, nH_free(j))
         got = 2.0d0*nc(ic_H2) + nc(ic_OH) + 2.0d0*nc(ic_H2O)
         if (carrier_solved(ic_Hp)) got = got + nc(ic_Hp)
         if (got .gt. nH_car) then
            over = got/max(nH_car, 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            sc = nH_car/max(got, 1.0d-300)
            nc(ic_H2)  = nc(ic_H2) *sc
            nc(ic_OH)  = nc(ic_OH) *sc
            nc(ic_H2O) = nc(ic_H2O)*sc
            if (carrier_solved(ic_Hp)) nc(ic_Hp) = nc(ic_Hp)*sc
            pct_cell_constrained(j) = .true.
            if (over .gt. limit_report) hit = .true.
         endif
         if (hit) pct_limited = pct_limited + 1
         fc(j,:) = nc/max(nrho(j), 1.0d-99)
      enddo
      end subroutine limit_to_element_budget

      ! ------------------------------------------------------------- !

      ! Write the solved carriers back into f_sp, moving the nucleus
      ! difference into the closure species of each element so that no
      ! element total and no mass density is changed by the step.  bsp_mass
      ! of a carrier is exactly the sum of its nuclei's masses, so the mass
      ! moved out of H I, O I and C I is the mass moved into the carriers.
      subroutine carrier_write_back(rho, f_sp, fc, nrho, nH_free,        &
                                    nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)    :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp) :: nd, nH2n, nOHn, nH2On, nCOn, nHpn, rest, held, sc
      integer  :: j, k, i0

      ! THE LOWER GHOSTS ARE NOT WRITTEN.  They hold the composition of the
      ! gas flowing in, which the ionization sweep pins from the handoff and
      ! which this operator reads as its base Dirichlet value; rewriting
      ! them with the base cell's own partition would erase the boundary
      ! condition the step is being driven by.
      do j = 1, N+Ng
         nd = rho(j)*n0
         if (nd .le. 0.0d0) cycle
         nH2n  = fc(j,ic_H2) *nrho(j)
         nOHn  = fc(j,ic_OH) *nrho(j)
         nH2On = fc(j,ic_H2O)*nrho(j)
         nCOn  = fc(j,ic_CO) *nrho(j)
         nHpn  = 0.0d0
         if (carrier_solved(ic_Hp)) nHpn = fc(j,ic_Hp)*nrho(j)
         f_sp(j,isp_H2)  = nH2n /nd
         f_sp(j,isp_OH)  = nOHn /nd
         f_sp(j,isp_H2O) = nH2On/nd
         f_sp(j,isp_CO)  = nCOn /nd
         if (carrier_solved(ic_Hp)) f_sp(j,isp_HII) = nHpn/nd

         ! The nuclei the carriers did NOT take are shared out over the
         ! element's other species IN PROPORTION TO WHAT THEY ALREADY HELD.
         ! Two reasons for the proportional rule rather than putting the
         ! whole remainder in the neutral stage: it cannot drive a stage
         ! negative, and it leaves the ionization fraction of the cell where
         ! it was, which is the right starting point for the sweep that
         ! re-solves it a moment later.  HeH+ is excluded from the hydrogen
         ! share because rescaling it would move a helium nucleus too.
         ! The nuclei the carriers did not take, shared over the stages
         ! that are still local. WHEN THE PROTON IS CARRIED IT LEAVES THAT
         ! SHARE: its density has just been written from the transported
         ! value above, and rescaling it here in proportion to what it held
         ! would immediately undo the transport. It moves to the taken side
         ! of the ledger instead, one nucleus each, exactly as H2 takes two.
         if (carrier_solved(ic_Hp)) then
            rest = nH_free(j) - 2.0d0*nH2n - nOHn - 2.0d0*nH2On - nHpn
            held = f_sp(j,isp_HI)*nd                                     &
                 + 2.0d0*f_sp(j,isp_H2p)*nd + 3.0d0*f_sp(j,isp_H3p)*nd
         else
         rest = nH_free(j) - 2.0d0*nH2n - nOHn - 2.0d0*nH2On
         held = f_sp(j,isp_HI)*nd + f_sp(j,isp_HII)*nd                   &
              + 2.0d0*f_sp(j,isp_H2p)*nd + 3.0d0*f_sp(j,isp_H3p)*nd
         endif
         if (held .gt. 0.0d0 .and. rest .gt. 0.0d0) then
            sc = rest/held
            f_sp(j,isp_HI)  = f_sp(j,isp_HI) *sc
            if (.not. carrier_solved(ic_Hp))                             &
               f_sp(j,isp_HII) = f_sp(j,isp_HII)*sc
            f_sp(j,isp_H2p) = f_sp(j,isp_H2p)*sc
            f_sp(j,isp_H3p) = f_sp(j,isp_H3p)*sc
         else if (held .le. 0.0d0) then
            f_sp(j,isp_HI) = max(rest, 0.0d0)/nd
         else
            f_sp(j,isp_HI)  = 0.0d0
            if (.not. carrier_solved(ic_Hp)) f_sp(j,isp_HII) = 0.0d0
            f_sp(j,isp_H2p) = 0.0d0
            f_sp(j,isp_H3p) = 0.0d0
         endif

         rest = nO_free(j) - nOHn - nH2On - nCOn
         i0   = melem_i0(iel_O)
         held = 0.0d0
         do k = 0, melem_top(iel_O)
            held = held + f_sp(j,mion_fsp(i0+k))*nd
         enddo
         if (held .gt. 0.0d0) then
            sc = max(rest, 0.0d0)/held
            do k = 0, melem_top(iel_O)
               f_sp(j,mion_fsp(i0+k)) = f_sp(j,mion_fsp(i0+k))*sc
            enddo
         else
            f_sp(j,mion_fsp(i0)) = max(rest, 0.0d0)/nd
         endif

         rest = nC_free(j) - nCOn
         i0   = melem_i0(iel_C)
         held = 0.0d0
         do k = 0, melem_top(iel_C)
            held = held + f_sp(j,mion_fsp(i0+k))*nd
         enddo
         if (held .gt. 0.0d0) then
            sc = max(rest, 0.0d0)/held
            do k = 0, melem_top(iel_C)
               f_sp(j,mion_fsp(i0+k)) = f_sp(j,mion_fsp(i0+k))*sc
            enddo
         else
            f_sp(j,mion_fsp(i0)) = max(rest, 0.0d0)/nd
         endif
      enddo
      end subroutine carrier_write_back

      ! ------------------------------------------------------------- !

      ! The largest change of a solved carrier column between two carrier
      ! states, over the physical cells.  ONE number for the whole grid, so
      ! that a bound written on it is a statement about the state and not
      ! about the cell that happens to hold the least of the gas.
      real(dp) function carrier_composition_displacement(fa, fb)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fa, fb
      real(dp) :: d
      integer  :: j, ic
      d = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            d = max(d, abs(fa(j,ic) - fb(j,ic)))
         enddo
      enddo
      carrier_composition_displacement = d
      end function carrier_composition_displacement

      ! ------------------------------------------------------------- !

      ! What ended a relaxation pass, in words, for the caller's log.
      function carrier_relax_outcome_text(outcome) result(txt)
      integer, intent(in) :: outcome
      character(len=64) :: txt
      select case (outcome)
      case (carrier_relax_movement_bound)
         txt = 'the movement bound, with the last step inside it'
      case (carrier_relax_fixed_point)
         txt = 'the carriers stopped moving at the fixed wind'
      case (carrier_relax_step_budget)
         txt = 'the step budget, with neither reached'
      case (carrier_relax_interval_refused)
         txt = 'the shortest admissible trial was not covered'
      case (carrier_relax_nothing_to_advance)
         txt = 'nothing to advance: no carrier, or no frozen background'
      case (carrier_relax_chemistry_refused)
         txt = 'the chemistry of the shortest admissible trial did not close'
      case default
         txt = 'undefined'
      end select
      end function carrier_relax_outcome_text

      ! ------------------------------------------------------------- !

      ! THE PRESSURE AND THE TEMPERATURE OF A COMPOSITION AT A FIXED
      ! CONSERVED STATE.  At a wind that does not move the quantities held
      ! are the conserved variables u = (rho, rho v, E): the thermal energy
      ! rho e = E - (rho v)^2/(2 rho) is fixed, and the pressure is the one
      ! the caloric equation of state of the CURRENT composition assigns
      ! to it (pressure_from_energy_density, which reads the composition
      ! get_species_densities has just refreshed, so a change of the
      ! molecular content changes the heat capacity and with it p and T at
      ! the same rho e).  A held pressure would describe a different gas
      ! from the one the hydrodynamic solve holds.
      subroutine pressure_and_temperature_at_fixed_conserved_state(u,     &
                                                     f_sp, p, T, ntot, ne)
      real(dp), dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: p, T
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: ntot, ne

      real(dp), dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii
      real(dp), dimension(1-Ng:N+Ng) :: nheiTR
      real(dp), dimension(1-Ng:N+Ng,n_mion) :: nm
      integer :: j

      nhei   = 0.0d0
      nheii  = 0.0d0
      nheiii = 0.0d0
      nheiTR = 0.0d0
      call get_species_densities(u(1,:), f_sp, nhi, nhii, nhei, nheii,    &
                                 nheiii, nheiTR, nm, ne, ntot)
      do j = 1-Ng, N+Ng
         p(j) = pressure_from_energy_density(j, u(1,j),                   &
                   u(3,j) - 0.5d0*u(2,j)*u(2,j)/u(1,j))
      enddo
      call comp_T_from_p(p, ntot, ne, T)
      end subroutine pressure_and_temperature_at_fixed_conserved_state

      ! ------------------------------------------------------------- !

      ! THE CHEMISTRY OF A COMPOSITION AT A FIXED CONSERVED STATE.
      !
      ! The carriers this operator transports are held at the values it
      ! wrote: the sweep is handed the transported partition and fixes
      ! x(H2), the oxygen carriers and, where the protons are carried, the
      ! ionized fraction, from the composition it reads.  What it returns
      ! is the chemistry OF that partition: the locally eliminated stages
      ! (H2+, H3+, HeH+, the helium and metal stages, and H+ where it is
      ! not carried), the electron density, and the rate coefficients and
      ! background of every cell.
      !
      ! The sweep changes the particle count, and with it the temperature
      ! the fixed thermal energy assigns to the gas, so one sweep at the
      ! entry temperature is not a closed thermochemical state.  The cycle
      ! (densities, p from rho e, T, sweep) is therefore repeated until T
      ! stops moving: the state handed back is one whose composition,
      ! temperature, pressure and rate coefficients belong together at the
      ! conserved variables the hydrodynamic solve holds.  The tolerance
      ! is the sweep's own reaction-residual tolerance (ieq_res_tol, 1e-6):
      ! a temperature change below it cannot move a composition the sweep
      ! has converged to that tolerance.
      !
      ! ok is false when a sweep left a cell non-finite; the caller then
      ! discards the step and the background this sweep wrote.
      subroutine equilibrate_chemistry_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, heat, cool, eta, ok,   &
                                             n_cycles)
      real(dp), dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: p, T
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: heat, cool
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: eta
      logical,                                  intent(out)   :: ok
      integer,                                  intent(out)   :: n_cycles

      integer,  parameter :: cycles_max = 5
      real(dp), parameter :: t_cycle_tol = 1.0d-6
      real(dp), dimension(1-Ng:N+Ng) :: ntot, ne, T_prev
      type(ioniz_eq_ledger) :: ledger
      integer :: k

      ok = .true.
      n_cycles = 0
      call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,     &
                                                             p, T, ntot, ne)
      do k = 1, cycles_max
         T_prev = T
         call ioniz_eq(T, u(1,:), f_sp, heat, cool, eta, ledger)
         n_cycles = k
         if (ledger%n_nonfinite .gt. 0 .or.                               &
             .not. all(f_sp(1:N,:) .eq. f_sp(1:N,:)) .or.                 &
             .not. all(abs(f_sp(1:N,:)) .le. huge(1.0d0))) then
            ok = .false.
            return
         endif
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                                             p, T, ntot, ne)
         if (maxval(abs(T(1:N) - T_prev(1:N))/max(T(1:N), 1.0d-300))      &
             .lt. t_cycle_tol) return
      enddo
      end subroutine equilibrate_chemistry_at_fixed_conserved_state

      ! ------------------------------------------------------------- !

      ! Relax the carriers towards their steady state at a fixed wind, for
      ! the Picard loop of the steady solver.  Same shape as
      ! relax_element_composition: a time scale for each cell, grown
      ! geometrically while the state keeps moving, with the steady mass
      ! flux in place of the instantaneous velocity so that a breathing
      ! base does not drive the relaxation.
      !
      ! WHAT IS HELD.  The conserved variables u of the state the
      ! hydrodynamic solve handed over.  The composition moves; rho, rho v
      ! and E do not; the pressure and the temperature beside the returned
      ! composition are the ones the caloric equation of state of THAT
      ! composition assigns to the unchanged thermal energy
      ! (pressure_and_temperature_at_fixed_conserved_state), so the
      ! conserved variables the next solve consumes and the primitive state
      ! beside them describe one gas.
      !
      ! THE BOUND is on the composition handed back, not on the transport
      ! step: every trial is taken on a copy, the chemistry of a kept step
      ! is closed on the new composition at the fixed conserved state, and
      ! only then is the displacement from the pass entry measured, on the
      ! largest H2 mixing ratio of the ENTRY state.  A trial whose returned
      ! state lies outside `trust`, whose interval the operator did not
      ! cover, whose chemistry did not close, or which left a species
      ! non-finite, is undone (composition and background) and retried at
      ! half the length down to relax_grow_min.  MEASURED (P2, hot-Uranus
      ! carrier reload): the shortest admissible trial still moves 6.5e-4 of
      ! the entry H2 maximum, so at a bound below that nothing is kept.
      !
      ! THE FIXED POINT is one at which a full-length trial no longer moves
      ! the carriers, with the chemistry closed on each side of it.
      subroutine relax_photochemical_composition(u, v, f_sp, p, T,       &
                                                 heat, cool, eta,         &
                                                 trust, drift, nstep,     &
                                                 outcome)
      real(dp), dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: v
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      ! The pressure and temperature of the composition handed back, at
      ! the conserved state u; consistent with f_sp whether or not a step
      ! was kept.
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: p, T
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: heat, cool
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: eta
      real(dp),                                 intent(in)    :: trust
      real(dp),                                 intent(out)   :: drift
      integer,                                  intent(out)   :: nstep
      integer, intent(out), optional :: outcome

      real(dp), dimension(1-Ng:N+Ng) :: rho, dt_code, ntot_e, ne_e
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fprev, fnow, fentry
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco
      real(dp), dimension(1-Ng:N+Ng,n_species) :: f_held
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar, nrho, wfac
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp), dimension(1-Ng:N+Ng) :: p_held, T_held, heat_held
      real(dp), dimension(1-Ng:N+Ng) :: cool_held, eta_held
      type(ion_rates), dimension(:), allocatable :: bg_held
      real(dp) :: drj, tdiff, tadv, grow, tscale, dmax, x_ref
      integer  :: j, k, ic, ending, step_status, refusal, n_cycles
      logical  :: refused, chem_ok
      type(element_census_state) :: cen_relax
      ! Why the last trial was undone; the ending names it when the
      ! shortest admissible trial is refused for that reason.
      integer, parameter :: refused_interval  = 1
      integer, parameter :: refused_bound     = 2
      integer, parameter :: refused_chemistry = 3

      drift  = 0.0d0
      nstep  = 0
      rho    = u(1,:)
      if (.not. thereis_mol .or. .not. carrier_transport .or.             &
          .not. bg_ready) then
         if (present(outcome)) outcome = carrier_relax_nothing_to_advance
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                                          p, T, ntot_e, ne_e)
         return
      endif
      call element_census_take('relax_photochemical_composition',         &
                               rho, f_sp, cen_relax)
      tscale = R0/v0
      call carrier_state(rho, f_sp, fprev, ntot, nrho, wfac, TK,        &
                         mbar, nH_free, nO_free, nC_free)
      fentry = fprev
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      do j = 1, N
         drj   = max((r_edg(j) - r_edg(j-1))*R0, 1.0d0)
         tdiff = drj*drj/max(maxval(Dco(j,1:n_carrier))               &
                            + kzz_cell(j), 1.0d-30)
         tadv  = drj/max(abs(v(j))*v0, 1.0d-30)
         dt_code(j) = min(tdiff, tadv)/tscale
      enddo
      dt_code(1-Ng:0)   = dt_code(1)
      dt_code(N+1:N+Ng) = dt_code(N)
      x_ref  = max(maxval(fentry(1:N,ic_H2)), 1.0d-30)
      ! The primitive state of the entry composition at u: what is handed
      ! back when no step is kept.
      call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,     &
                                                          p, T, ntot_e, ne_e)
      grow   = 1.0d0
      ending = carrier_relax_step_budget
      fnow   = fentry
      allocate(bg_held(lbound(bg_cell,1):ubound(bg_cell,1)))
      carrier_rows_advect = .true.
      do k = 1, relax_maxstep
         f_held    = f_sp
         bg_held   = bg_cell
         p_held    = p;     T_held    = T
         heat_held = heat;  cool_held = cool;  eta_held = eta
         refused   = .false.
         refusal   = 0
         call photochemical_transport_step(rho, v, f_sp, dt_code*grow,    &
                                           step_status, trial = .true.)
         if (step_status .ne. carrier_interval_covered) then
            refused = .true.;  refusal = refused_interval
         endif
         if (.not. refused) then
            call equilibrate_chemistry_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, heat, cool, eta,       &
                                             chem_ok, n_cycles)
            if (.not. chem_ok) then
               refused = .true.;  refusal = refused_chemistry
            endif
         endif
         if (.not. refused) then
            call carrier_state(rho, f_sp, fnow, ntot, nrho, wfac,        &
                               TK, mbar, nH_free, nO_free, nC_free)
            if (carrier_composition_displacement(fnow, fentry)/x_ref      &
                .gt. trust) then
               refused = .true.;  refusal = refused_bound
            endif
         endif
         if (refused) then
            ! The trial and everything it wrote are undone: composition,
            ! background, primitive state.
            f_sp = f_held
            bg_cell = bg_held
            p = p_held;  T = T_held
            heat = heat_held;  cool = cool_held;  eta = eta_held
            fnow = fprev
            if (grow .le. relax_grow_min) then
               select case (refusal)
               case (refused_interval)
                  ending = carrier_relax_interval_refused
               case (refused_chemistry)
                  ending = carrier_relax_chemistry_refused
               case default
                  ending = carrier_relax_movement_bound
               end select
               exit
            endif
            grow = max(0.5d0*grow, relax_grow_min)
            cycle
         endif
         nstep = nstep + 1
         dmax  = carrier_composition_displacement(fnow, fprev)
         if (dmax/x_ref .lt. relax_tol .and. grow .ge. 1.0d0) then
            ending = carrier_relax_fixed_point
            exit
         endif
         fprev = fnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo
      carrier_rows_advect = .false.
      deallocate(bg_held)
      call carrier_state(rho, f_sp, fnow, ntot, nrho, wfac, TK,          &
                         mbar, nH_free, nO_free, nC_free)
      drift = 0.0d0
      pct_drift_j  = 0
      pct_drift_ic = 1
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            if (abs(fnow(j,ic) - fentry(j,ic)) .gt. drift) then
               drift = abs(fnow(j,ic) - fentry(j,ic))
               pct_drift_j  = j
               pct_drift_ic = ic
            endif
         enddo
      enddo
      drift = drift/x_ref
      call element_census_verify(cen_relax, rho, f_sp, rho_is_fixed=.true.)
      if (present(outcome)) outcome = ending
      end subroutine relax_photochemical_composition

      ! End of module

      end module diffusive_photochemistry

      module diffusive_photochemistry
      ! Vertical transport of the molecular carriers, solved together with
      ! the chemistry that makes and destroys them:
      !
      !   d n_i/dt + (1/r^2) d/dr [ r^2 ( n_i v + Phi_i ) ] = P_i - L_i
      !
      !   Phi_i = - n_tot (D_i + K_zz) d f_i/dr
      !           - n_i D_i [ 1/H_i - 1/H_atm ]        (f_i = n_i/n_tot)
      !
      ! for the molecular carriers i = H2, OH, H2O, CO, and, in the same
      ! block system, the carried ionization stages of hydrogen and helium
      ! (H+, He+, He++) and the He 2^3S level, whose rows are the stage
      ! fluxes of ionization_stage_transport (section 2).  This header
      ! states what was built and what was measured, not what was intended.
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
      ! for H2 at a cool base: at the HD 189733 b 1-microbar
      ! level tau_chem(H2)/tau_adv is 0.20 to 1.67 on H/v, and 16.4 when the
      ! two are compared at the same pressure, against 1e-3 on HD 209458 b.
      !
      ! WHAT WAS THEN MEASURED IN EXHALE'S OWN STRUCTURE, AND IT IS NOT THE
      ! SAME NUMBER.  That ratio was built from the PHOTOCHEMICAL model's
      ! scale height and velocity at that pressure.  In an EXHALE run of the
      ! same planet the base cell has v = 0 EXACTLY -- the lower boundary
      ! condition puts it there -- and the flow time r/|v| in the cells above
      ! it is 1e8-2e8 s against tau_chem(H2) = 178-552 s.  The diffusive time
      ! of the base cell is 8.0e5 s with K_zz = 0, i.e. 4500 chemical times,
      ! and 634 s with K_zz = 1e9, i.e. 3.6.  So in this code's structure the
      ! base partition is a LOCAL quantity, and transport is not the term
      ! that decides it; what transport does decide is the profile of a
      ! species with no chemistry of its own, which here is CO.
      !
      ! That is why every run writes tau_chem, tau_adv AND the diffusive time
      ! of each cell to output/Oxygen_chemistry.txt: the regime is a
      ! measurement of the run, not a property of the option.
      !
      ! ---------------------------------------------------------------
      ! 2. WHICH SPECIES ARE TRANSPORTED, AND WHY THE OTHERS ARE NOT
      !
      ! THE SET IS A PROPERTY OF THE RUN (carrier_set_init), up to eight
      ! rows a cell (n_carrier_max):
      !  - H2, with "Molecular carrier transport: True" in a molecular run;
      !    OH, H2O and CO behind it when the oxygen chemistry is on;
      !  - H+, and with helium He+ and He++, with "Ionization transport:
      !    True".  They are written as fractions of their element's nuclei
      !    and ride on the element nucleus flux, with the eddy term of the
      !    stage (ionization_stage_transport, equation 1); H I and He I are
      !    the closing stages, one minus the carried ones, and not rows;
      !  - the He 2^3S level, with "He 2^3S transport: True" on top of the
      !    carried helium stages, as a fourth state of the helium partition.
      ! A configuration that solves a later row without an earlier one (the
      ! stages without the molecular carriers, say) keeps the unsolved rows
      ! as identity rows (carrier_solved).
      !
      ! THE MOLECULAR CARRIERS.  Between them H2, OH, H2O and CO carry every
      ! H nucleus that is not atomic or in a molecular ion, and every O and C
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
      !  - the ion stages of the metals, and the hydrogen and helium stages
      !    in a run without "Ionization transport".  Their balance is local,
      !    and gate G2 requires that limit to be reproduced.  They are where
      !    the transported carriers' nuclei come from and go to: the
      !    write-back below moves nuclei between a carrier and the closure
      !    species of its element, so no element can drift.
      !
      ! The one term this drops that is not small everywhere is the AMBIPOLAR
      ! field.  The four molecular carriers are neutral, so eq. (2) of the
      ! design has no (Z e E)/(k T) term for any of them; for them the field
      ! is not neglected, it is absent.  The carried ionization stages ARE
      ! charged, and for them it is neglected: they move with their
      ! element's nucleus flux plus the eddy term, with no molecular
      ! diffusion and no ambipolar drift against their own neutral
      ! (carrier_diffusivities sets their D to zero; the validity
      ! of the one-velocity approximation is measured in
      ! ionization_stage_transport).  settling_coefficient in
      ! binary_element_diffusion builds the field for the element operator.
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
      ! CO has no chemical source: it is frozen chemically, and
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
      ! 5. BOUNDARY CONDITIONS
      !
      ! WHETHER THE BASE FACE IS AN INFLOW is decided by the wind's mass
      ! flux and not by the base cell's velocity.  Measured on the converged
      ! hot Uranus: the cell-centred rho v r^2 of cell 1 is 160 to 200 times
      ! the wind's flux INWARD while the face flux
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
      ! compute, so the operator does not impose it: the base face carries
      ! zero diffusive flux, and the advective inflow, which is the base face
      ! mass flux times the ghost composition, carries the base cell's own
      ! partition because the ghost is given that partition (df/dr = 0 at the
      ! face), which is the same statement.  The one exception is data from
      ! below: where a handoff states the base H2 partition (base.inp
      ! q_H2_base, or the profile's value at the matching level), the inner
      ! ghosts hold that partition and the inflowing base face carries it in
      ! (carrier_set_init, base_h2_composition_imposed).  The carried
      ! ionization stages are never stated by a handoff and keep the
      ! zero-gradient ghost in every configuration.
      !
      ! That is what makes the A/B gate of design sec. 6.2 a test.  Had the
      ! base H2 fraction been pinned from the same handoff the gate compares
      ! against, the gate would measure nothing.
      !
      ! OUTER: no diffusive, eddy or settling flux crosses the outer face,
      ! and the carriers leave it by advection alone.  carrier_face_coefficients
      ! gives the reason, a ghost-built face flux measured to draw carriers
      ! inward; no output file flags a carrier that is not negligible there.
      !
      ! ---------------------------------------------------------------
      ! 6. DISCRETIZATION AND THE SIMPLEX
      !
      ! Cell-centred f_i, face-centred fluxes on r_edg, one backward-Euler
      ! step per call, Newton with a backtracking line search, block
      ! tridiagonal in space (n_carrier_max x n_carrier_max = 8x8 blocks,
      ! the unsolved rows identity rows: the transport is diagonal in
      ! species and the chemistry is dense within a cell).  Faces 0 and N
      ! carry zero diffusive flux.  The drift term switches from central to
      ! donor-cell on the same Peclet test the element operator uses.
      !
      ! THE MATERIAL ADVECTION OF A CARRIER IS NOT IN THESE ROWS.  It is the
      ! divergence of the hydrodynamic face mass flux, F_rho(j) Y_c(j)/m_c,
      ! and it is taken with the mass row itself inside the Runge-Kutta
      ! stages, on the same faces, areas, volumes and time step
      ! (species_face_flux.f90).  A cell can then only lose the carrier that
      ! the mass row says it loses, whatever the velocity field does, and a
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
      ! (1996) for the H2 band.

      use mol_rates, only: h2_thermochemistry_init, h2_thermochemistry_ready
      use global_parameters
      use caloric_eos, only: adiabatic_index_at_T,                        &
                             h2_rovibrational_energy_and_heat_capacity
      use grav_func,     only: Dphi
      use species_table, only: n_bsp, bsp_fsp, bsp_mass,                 &
                               bsp_nHe, bsp_is_excited_level,            &
                               isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                               isp_HeIII, isp_HeTR,                      &
                               isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                               isp_OH, isp_H2O, isp_CO,                  &
                               n_mion, mion_fsp, melem_i0, melem_top,    &
                               n_melem, iel_O, iel_C,                    &
                               bsp_charge, bsp_nH, bsp_nO, bsp_nC,       &
                               mion_elem, mion_stage
      use binary_element_diffusion, only: hard_sphere_pair_diffusion,     &
                          advected_carrier_reset,                        &
                          advected_carrier_register,                     &
                          advected_carrier_count,                        &
                          advected_carrier_fractions,                    &
                          advected_carrier_state_ready,                  &
                          advected_carrier_state_consumed,               &
                          carrier_mass_fractions, mixture_mass_sum,      &
                          element_nucleus_counts
      use ion_cell_state, only: ieq_cell, ion_rates
      use oxygen_rates, only: rk_D1_Hep_CO,                              &
                              rk_CO_radiative_association
      use co_self_shielding_table, only: co_shield_tex_limit_K
      use System_HeH_mol, only: set_mol_coeffs, set_oxygen_coeffs,       &
                                mol_heh_rows, oxygen_carrier_rows,       &
                                oj3, oj4, oj5, oj7, mk_ion_H2
      use water_photolysis, only: n_fuv_band
      use utils_ion_eq, only: fuv_lw_photon_field
      use utils, only: calc_ne
      use ionization_equilibrium, only: bg_cell, bg_ready, finite_real,   &
                                        ioniz_eq, ioniz_eq_ledger,       &
                                        set_ioniz_eq_state_may_be_refused, &
                                        transported_rows_exist,          &
                                        ieq_rate_state,                  &
                                        save_ieq_rate_state,             &
                                        restore_ieq_rate_state
      ! THE H AND He IONIZATION ROWS OF AN ATOMIC GAS, the very rows the
      ! local sweep solves there (System_HeH_TR and System_HeH reach them
      ! through ion_residual_core), so that a stage source and the sweep's
      ! own balance are one expression and cannot drift apart.
      use ion_residual_core, only: heh_tr_rows, tr_triplet_row_channels
      use molecular_reaction_heat, only: species_formation_energy
      ! He <-> H charge exchange, Huang et al. (2023) Table 4 group B. It is
      ! a term of the H+ and He+ balances in EVERY system that carries
      ! helium, so it is a term of the stage rows too.
      !
      ! The metal charge exchange of the same table -- group A, metal + H
      ! and H+, active whenever metals are present; group C, metal + He and
      ! He+, under cx_full; group D, metal + metal, which reaches no H or
      ! He stage; and the group E electron capture O2+ + H0, under its own
      ! scale -- is a term of the same two balances, and the local sweep
      ! adds it to its rows through the same reaction set. The stage
      ! sources take it in the basis their unknowns are written in, the
      ! stage densities themselves, and cx_set_cell hands the evaluator
      ! this cell's rate coefficients.
      use charge_exchange, only: he_h_cx_fvec, cx_set_cell,               &
                                 charge_exchange_stage_sources
      use ion_cell_state, only: ion_rates
      use caloric_eos, only: pressure_from_energy_density
      use steady_residual_mod, only: carrier_row_scale,                  &
                                     face_mass_flux_of_state,            &
                                     carrier_enthalpy_divergence
      use species_advective_transport, only: species_face_fraction,      &
                                     species_face_flux
      ! THE TRANSPORT OF AN IONIZATION STAGE, written on its element's own
      ! nucleus flux.  The stage rows of this operator -- H II per hydrogen
      ! nucleus today -- are the divergence of that flux and nothing is
      ! rebuilt here, which is the condition under which the stage fluxes
      ! of one element sum to the element's nucleus flux face by face
      ! (ionization_stage_transport, equation 2).
      use ionization_stage_transport, only:                              &
                                  ionization_stage_face_flux,            &
                                  ionization_stage_face_jacobian,        &
                                  stage_simplex_projection,              &
                                  stage_simplex_sum_over_limit,          &
                                  stage_fraction_under_zero,             &
                                  hydrogen_and_helium_nucleus_face_flux
      ! The one spherical geometry of this grid: face areas r_edg^2 and the
      ! exact shell volumes.  Every divergence below divides by these, so
      ! the advective and the diffusive halves of one carrier row are the
      ! divergence of one flux.
      use grid_construction, only: spherical_face_area_and_cell_volume
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
                                   ien_H, ien_He, ien_O,                 &
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
      public :: carrier_diffusivities_of_state
      public :: carrier_set_init
      public :: n_carrier_max
      ! Element closure: the write-back's invariant is that the element
      ! totals it is handed come back unchanged, on the two routines that
      ! make the pair.
      public :: carrier_state, carrier_write_back
      ! The advective term of the carrier rows, exposed so that an
      ! acceptance test measures it against the Runge-Kutta stage's own
      ! update of the same state and not against a transcription of it.
      public :: carrier_advective_divergence, carrier_mass_amu
      ! The two face coefficients of the advective term and the predicate
      ! that says whether the base ghosts are data or a copy: they carry the
      ! boundary derivative the block-tridiagonal assembly uses.
      public :: carrier_advective_face_coefficients
      public :: carrier_base_composition_imposed
      ! The diagnostic interventions, so that a test driver can state what
      ! each of them does to the operator without running the binary.
      public :: carrier_face_coefficients
      public :: l22b_setup, carrier_outflow_ghost
      public :: l22b_advective_band_derivative
      ! Save and restore of the module arrays carrier_steady_residual
      ! overwrites, so that the balance of a state can be MEASURED without
      ! the measurement being an operation on the state.
      public :: carrier_module_state
      public :: save_carrier_module_state, restore_carrier_module_state
      public :: carrier_module_state_matches
      ! The thermochemical state one trial writes, held aside and put back
      ! when the trial is undone (the type's header enumerates it), public
      ! so that a test can state the round trip on the production pair and
      ! not on a copy of it.
      public :: thermochemical_state
      public :: save_thermochemical_state, restore_thermochemical_state
      public :: carrier_steady_residual
      public :: carrier_enthalpy_active, carrier_enthalpy_divergence_of_state
      ! The two parts of that enthalpy flux, public for the tests
      ! (src/tests/carrier_enthalpy_flux, src/tests/ionization_stage_enthalpy_flux).
      public :: molecular_carrier_enthalpy_active,                       &
                molecular_carrier_enthalpy_face_flux,                    &
                ionization_stage_enthalpy_active,                        &
                ionization_stage_enthalpy_face_flux
      ! Test only: the diffusion coefficients the carrier face flux is
      ! formed from (src/tests/carrier_enthalpy_flux).
      public :: carrier_diffusivities
      ! The terms of every carrier row of the last assembly, written on
      ! request (EXHALE_CARRIER_ROW_TERMS=1) by the certification that
      ! measured them.  A record of a measurement; nothing reads it back.
      public :: carrier_row_terms_on, carrier_row_terms_write
      ! The H2 content each cell's own carrier row would settle at, for the
      ! molecular seed: the row's chemistry alone, taken from the routine
      ! that assembles the row.
      public :: carrier_h2_chemical_root
      public :: carrier_drift_location
      ! The two row scales of the last assembly side by side: the
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
      ! Test only: writes the element budget of a synthetic cell, from which
      ! the projection onto the budget face of the coupled solve's unknown
      ! box is formed.
      public :: carrier_headroom_set_for_test
      public :: n_carrier, carrier_name, ic_H2, ic_OH, ic_H2O, ic_CO
      public :: ic_Hp, ic_HeII, ic_HeIII, ic_HeTR, carrier_solved
      ! The stage sum identity of the last stationary evaluation: the
      ! certification reports it beside the rows it was measured with.
      public :: ionization_stage_sum_measure
      ! The ionization stage rows of the operator: which carriers are
      ! stages, the element nucleus density each is a fraction of, the
      ! frozen face state their fluxes ride on, one stage's face flux and
      ! its two derivatives, the identity that sums the stages of one
      ! element, and the projection onto each element's simplex.  The
      ! helium rows ride on the HELIUM nucleus flux and the two helium
      ! stages share one simplex and one closing stage.
      public :: carrier_is_ionization_stage, carrier_stage_element
      public :: carrier_is_nucleus_fraction
      ! The energy the transported He 2^3S population carries in its
      ! excitation, measured against the energy row (certification).
      public :: helium_metastable_excitation_energy_divergence
      public :: carrier_nucleus_reference, carrier_stage_face_state
      ! THE HELIUM THE TWO IONIZED STAGES MAY HOLD, and the neutral helium
      ! a trial partition of them leaves.  They are the admissible set of
      ! the helium rows and the closure the chemistry rows read, written
      ! once.
      public :: carrier_helium_available_to_stages
      public :: carrier_helium_neutral_of_partition
      public :: carrier_helium_singlet_breach
      public :: carrier_helium_background_set_for_test
      public :: carrier_stage_face_flux, ionization_stage_nucleus_sum
      public :: carrier_ionization_stage_projection
      ! The element table of the carriers and the perturbation the chemistry
      ! rows are differentiated with: the second one produces a directional
      ! derivative of the first.
      public :: carrier_element_reference_density,                       &
                carrier_source_derivative_step, carrier_source
      ! The domain of the one-sided CO destruction model, cumulative over
      ! the whole run: the cells in which tau_dest << tau_res << tau_form
      ! fails, which is the ordering that makes a CO row carrying
      ! destruction and no formation legitimate.
      public :: carrier_co_domain_record, carrier_co_domain_f_dom
      ! Test only: writes distinct amounts into the two ledger families of
      ! the domain record, which separates the families; the run total
      ! composes them (a sum for the counts, a maximum for the worst
      ! ratio), and a refused attempt does not roll it back.
      public :: carrier_co_domain_perturb_for_test
      ! Whether an element constraint moved a given cell in the last
      ! limiter call: the CO thermal ceiling marks the cell it acts in,
      ! which is the cell whose OH and H2O rows it moves.
      public :: carrier_cell_is_constrained, limit_to_element_budget
      ! ACCEPTANCE OF THE STATE THE CARRIER SOLVE RETURNS, and the
      ! classification of a stalling iteration it keeps separate from it.
      ! Both are decisions taken on their arguments alone, so the retry
      ! controller can ask for a verdict on a trial state without running
      ! a solve.
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
      ! THE CARRIER-LOCAL CHECKPOINT AND THE RETRY CONTROLLER.  The
      ! checkpoint is the carrier fractions plus every module array and
      ! counter one transport attempt writes; carrier_transport_interval is
      ! the step's body with the stop taken out of it, so the exhaustion
      ! branch can be reached and what the controller did read.
      public :: carrier_checkpoint
      public :: carrier_checkpoint_take, carrier_checkpoint_restore
      public :: carrier_checkpoint_matches
      public :: carrier_transport_interval
      ! WHETHER THE INTERVAL WAS COVERED, as a value the caller reads
      ! instead of inferring it from a stop.  photochemical_transport_step
      ! returns it, so the step
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
      ! ABSOLUTE SMALLNESS OF A CARRIER, one number for the whole module:
      ! the fraction of the free reservoir of its own element below which a
      ! carrier density is treated as absent.  Twenty decades is the same
      ! statement the elemental transport rows make about an element
      ! (1e-20 rho X_base, binary_element_diffusion).  It is a choice, not a
      ! round-off limit: a density stored as its own double stays resolved
      ! far below eps of the reservoir, and the test reads the abundance
      ! alone, not the production or the incoming flux that could raise it.
      ! It carries the row-scale floor of carrier_residual and the absolute
      ! floor under the certification's carrier row.
      real(dp), parameter, public :: carrier_absent_fraction = 1.0d-20

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
      ! its retry consume this type.
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
      ! between carried and local is a measurement and not a convention.
      !
      ! n_carrier_max dimensions every array; n_carrier is how many of them
      ! this run solves, set once by carrier_set_init.  The inactive columns
      ! stay at zero -- the species behind them do not exist without the
      ! oxygen chemistry -- so the whole-array statements below are exact.
      !
      ! THE IONIZATION STAGES ARE THE NEXT THREE, AND THAT IS FOR A REASON.
      ! Their indices have to be constants -- the ionization rows, the
      ! write-back and the element budget all name them -- and the four
      ! molecular indices have to keep the values they have, so the only
      ! place left is behind them.  The He 2^3S level ("He 2^3S
      ! transport") is the last, behind the stages it shares the helium
      ! simplex with.
      ! That makes the solved set non-contiguous in the one configuration
      ! that carries H2 and the stages without the oxygen cycle, which is why
      ! carrier_solved exists below: n_carrier is the LAST index solved and
      ! carrier_solved says which of 1..n_carrier actually are.
      !
      ! WHICH STAGES THEY ARE IS THE PARTITION OF TWO ELEMENTS, not a list
      ! of species: x(H II) per hydrogen nucleus closes against H I (and
      ! the hydrogen the molecules hold), and x(He II) with x(He III) per
      ! helium nucleus close against He I (and the helium of HeH+).  The
      ! closing stage of each element is not a row: it is one minus the
      ! carried ones, which is what makes the stage fluxes of an element
      ! sum to that element's own nucleus flux.
      !
      ! THE He 2^3S LEVEL IS A FOURTH STATE OF THE HELIUM PARTITION, not a
      ! stage: it has the charge of He I.  Carried ("He 2^3S transport"),
      ! x3 = n(2^3S)/n(He nuclei) rides on the helium nucleus flux in the
      ! same variable as the two ionized stages, and the closing state of
      ! helium is then the ground singlet, one minus the three carried
      ! ones (and the helium of HeH+).  carrier_is_nucleus_fraction is the
      ! set of carriers written that way; carrier_is_ionization_stage the
      ! three of them that are ionization stages.
      integer, parameter :: n_carrier_max = 8
      integer, save      :: n_carrier = 4   ! set by carrier_set_init
      integer, parameter :: ic_H2 = 1, ic_OH = 2, ic_H2O = 3, ic_CO = 4
      integer, parameter :: ic_Hp = 5
      integer, parameter :: ic_HeII = 6, ic_HeIII = 7
      integer, parameter :: ic_HeTR = 8
      ! Names without blanks: the row-term record and the logs write them
      ! as one whitespace-delimited field.
      character(len=6), parameter :: carrier_name(n_carrier_max) =       &
           (/ 'H2    ', 'OH    ', 'H2O   ', 'CO    ', 'H+    ',          &
              'He+   ', 'He++  ', 'He2^3S' /)
      ! Which of 1..n_carrier this run solves. A carrier that is not solved
      ! keeps the identity row nrho/dt on the diagonal and a zero residual,
      ! so its unknown does not move and the block structure the Thomas
      ! sweep needs is unchanged; it costs one n_carrier_max-square block
      ! entry and no chemistry evaluation. Set once by carrier_set_init.
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
      ! A carrier that holds the WHOLE of its element sits exactly on the
      ! face of the simplex -- all the carbon is CO in a molecular base --
      ! so a limiter that counted every round-off would report half the grid
      ! and mean nothing. Only an overshoot above this relative size is
      ! counted; the clamp itself is applied whatever the size.
      real(dp), parameter :: limit_report = 1.0d-10

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
      ! MEASURED on the molecular column: the Newton reaches 2.6e-16 of the
      ! full row terms, about one eps, so the
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
      ! state and the steady carrier residual of the present rows is at or
      ! below carrier_accept_tol.  STEP_BUDGET: relax_maxstep trials were taken without
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
      ! THE THERMOCHEMICAL CLOSURE (equilibrate_chemistry_at_fixed_
      ! conserved_state): its cycle budget, its temperature tolerance, the
      ! named reasons it returns, and the run-wide counts the relaxation
      ! report prints (the cost of the closure is measured, not guessed).
      ! chem_cycles_cap_for_test, when nonnegative,
      ! replaces the budget so that a test can force exhaustion.
      ! THE BUDGET AND THE TOLERANCE, BOTH MEASURED (on
      ! LHS 1140 b molecular_scalar_gj1132_kzz1e9/HeH2.13 with the closure
      ! traced cycle by cycle, EXHALE_CARRIER_DEBUG=1).  The cycle
      ! contracts hard for two or three cycles and then STOPS: the
      ! increments of six consecutive closures of that run are
      !
      !   7.6e-05 -> 1.1e-05 -> 2.8e-06, and from there a band
      !   1.7e-06 to 3.4e-06 that does not narrow over thirty more cycles,
      !
      ! and the cell carrying it is not the molecular layer but the far
      ! wind at r = 10 to 11 R_p, where it is the equilibrium sweep's own
      ! acceptance noise and not a state still moving.  So 1e-6 is BELOW
      ! the floor the sweep defines: the closure then converges only when
      ! the band happens to dip, which took 8, 22 and 28 cycles in three
      ! calls of that run and did not happen within 31 in a fourth.  The
      ! tolerance is set above that floor with a factor of three margin,
      ! and the budget at twice the 6 cycles the slowest of those calls
      ! needed to reach it.  With 5 and 1e-6 every call spent its budget,
      ! every trial of relax_photochemical_composition was undone for it,
      ! and the H2 carrier of that case took no transport step at all.
      integer,  parameter, public :: chem_cycles_max = 12
      real(dp), parameter, public :: chem_cycle_tol  = 1.0d-5
      integer,  parameter, public :: chem_closure_converged      = 0
      integer,  parameter, public :: chem_closure_exhausted      = 1
      integer,  parameter, public :: chem_closure_nonfinite      = 2
      ! 3 was "a cell left the element simplex", retired: the
      ! count it read is a property of the root search and not of the state
      ! the closure hands back (the closure body says why), so no closure
      ! can return it and the value is left unused rather than reassigned.
      integer,  parameter, public :: chem_closure_not_admissible = 4
      ! The sweep could not solve the lower boundary at this state: a
      ! ghost's base handoff partition or the ghost composition fixed point
      ! was left open (ieq_state_may_be_refused, ionization_equilibrium.f90).
      integer,  parameter, public :: chem_closure_boundary_open  = 5
      integer, public :: chem_cycles_cap_for_test = -1
      integer, public :: n_chem_closure_cycles   = 0
      integer, public :: n_chem_closures_reached = 0
      integer, public :: n_chem_closures_refused = 0
      integer, public :: n_chem_last_reason      = chem_closure_converged
      ! The cell, carrier and absolute change that attained the movement
      ! bound when the last trial was refused on it:
      ! a bound that refuses a pass says which cell refused it
      ! (carrier_worst_composition_change).  Diagnostic only.
      integer,  public :: bound_last_j = 0, bound_last_ic = 1
      ! The measure the bound was refused on -- a relative particle-count
      ! change by default, the retired carrier-fraction ratio under
      ! EXHALE_CARRIER_BOUND_FRACTION=1 -- and the entry value it is
      ! relative to.
      real(dp), public :: bound_last_dabs = 0.0d0
      real(dp), public :: bound_last_entry = 0.0d0
      logical,  public :: bound_last_fraction = .false.
      ! ---------------------------------------------------------------- !
      ! A DIAGNOSTIC EXPERIMENT, DEFAULT OFF.  WHAT IT MEASURES: the
      ! movement bound is one
      ! scalar over a column that holds a slow H2 front and a far wind, and
      ! the cell that attains it is the front while the cell whose carrier
      ! row refuses the certification sits decades inside the allowance.
      ! The experiment asks whether the refusing cell descends when the
      ! front's cells no longer choose the trial length.  It is NOT a
      ! candidate for production and certifies nothing: no cell is ever
      ! excluded from the certification, only from the choice of the trial
      ! interval.
      !
      ! EXHALE_L22_MASK=<file> lists, one integer per line (blank lines and
      ! lines opening with # ignored), the physical cells OMITTED FROM THE
      ! INTERVAL SELECTION.  EXHALE_L22_MASK_MODE selects what the omission
      ! means:
      !   veto  -- the particle-count bound is still evaluated on the
      !            omitted cells and still refuses the trial, so the
      !            feasible set and therefore the accepted interval are the
      !            ones of the unmasked run; only the naming of the
      !            controlling cell changes.  This is the null comparison.
      !   waive -- the bound is not evaluated on the omitted cells, so the
      !            interval is chosen by the selection set alone.  What
      !            still holds the trial then: the transport step's own
      !            interval coverage, the positivity and element and charge
      !            feasibility of the chemistry closure, the finiteness and
      !            gas conditions of thermal_state_admissible, the element
      !            census of the pass, and the safety stop below.
      ! SAFETY STOP (waive only): a trial in which any omitted cell's
      ! particle count (n_tot + n_e) changes by more than a factor two ends
      ! the relaxation at the last accepted state with
      ! carrier_relax_mask_safety_stop.  An empty selection set is not a
      ! comparison and ends the relaxation the same way with its own
      ! message.
      integer, parameter, public :: l22_mask_off   = 0
      integer, parameter, public :: l22_mask_veto  = 1
      integer, parameter, public :: l22_mask_waive = 2
      integer, public :: l22_mask_mode = l22_mask_off
      logical, save   :: l22_mask_ready = .false.
      logical, dimension(:), allocatable, save :: l22_mask_omit
      integer, public :: l22_mask_n_omit = 0, l22_mask_n_select = 0
      ! What the last trial measured on each of the two sets, for the log:
      ! the largest relative particle-count change of the selection set and
      ! of the omitted set, with the cells attaining them.
      integer,  public :: l22_last_j_select = 0, l22_last_j_omit = 0
      real(dp), public :: l22_last_d_select = 0.0d0
      real(dp), public :: l22_last_d_omit   = 0.0d0

      ! ---- Measurements of the carrier Jacobian and of what a relaxation
      ! pass hands back, each switched on by an environment key and OFF by
      ! default, so that a run with none of them set integrates the operator
      ! this module states.  The boundary of the carrier column is not among
      ! them: it is one rule, stated at carrier_face_mass_fraction and
      ! carrier_face_coefficients, and no key decides it.
      !
      ! (B) THE DEFERRED RECONSTRUCTION TERMS.  The advective entries of the
      ! block-tridiagonal Jacobian are the FIRST-ORDER donor-cell
      ! linearization (see the assembly in solve_carriers); the limited
      ! slopes of species_face_fraction are not differentiated.
      ! EXHALE_L22B_JAC_RECON=1 replaces those entries by a central
      ! difference of carrier_advective_divergence itself, restricted to the
      ! tridiagonal band, and reports the weight the band cannot hold.
      ! EXHALE_L22B_JAC_ACTION=<file> writes, at the first Jacobian assembly
      ! of the run, the action of the assembled matrix on a direction
      ! concentrated on the cells EXHALE_L22B_JAC_CELLS=lo,hi beside a
      ! central difference of the full residual along the same direction.
      !
      ! (C) EXHALE_L22B_DISPLACEMENT=1 reports what the composition the
      ! relaxation hands back did to the primitive state the wind reads:
      ! the largest relative change of p, T, the mean mass per particle and
      ! the particle count over the column, with the cells attaining them.
      logical, public :: l22b_jac_recon        = .false.
      logical, public :: l22b_jac_action       = .false.
      logical, public :: l22b_displacement     = .false.
      integer, public :: l22b_jac_lo = 0, l22b_jac_hi = 0
      character(len=512), save :: l22b_jac_file = ' '
      logical, save   :: l22b_ready = .false.
      logical, save   :: l22b_jac_action_done = .false.
      ! The largest advective Jacobian entry the tridiagonal band cannot
      ! hold, relative to the band of the same row, over the last assembly
      ! that formed the reconstruction terms.
      real(dp), public :: l22b_recon_dropped = 0.0d0

      ! The ledger counts of the last closure's last sweep (diagnostics).
      integer, public :: chem_last_offsimplex = 0, chem_last_nonfinite = 0
      integer, public :: chem_last_mol_clamped = 0
      real(dp), public :: chem_last_viol_worst = 0.0d0
      real(dp), public :: chem_last_increment    = 0.0d0
      public :: equilibrate_chemistry_at_fixed_conserved_state
      public :: chem_closure_reason_text
      integer, parameter, public :: carrier_relax_nothing_to_advance = 4
      ! The chemistry of a kept transport step did not close (a sweep that
      ! left a cell non-finite), and the shortest admissible trial did no
      ! better: the entry state is handed back.
      integer, parameter, public :: carrier_relax_chemistry_refused = 5
      ! The diagnostic interval mask waived the bound on a cell that then
      ! moved by more than a factor two in one trial, or its selection set
      ! was empty.  Only reachable with EXHALE_L22_MASK set.
      integer, parameter, public :: carrier_relax_mask_safety_stop = 6

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
      ! THE THREE IONIZATION STAGES OF EVERY METAL ELEMENT [cm^-3], frozen
      ! with the rest of the background and indexed (cell, canonical
      ! element).  They are the reactant densities of the metal charge
      ! exchange of Huang et al. (2023) Table 4, which moves H and He
      ! between their own stages and so is a term of the transported stage
      ! sources.  A metal stage is NOT a transported unknown: this operator
      ! moves the hydrogen and helium stages and the molecular carriers,
      ! and the ionization sweep that follows re-solves the metal
      ! partition from its own balance, so within one solve these are a
      ! state and not a moving target, exactly as the molecular ions are.
      real(dp), dimension(:,:), allocatable :: cbg_nm0, cbg_nm1, cbg_nm2
      ! THE NUCLEUS DENSITY OF THE ELEMENT EACH IONIZATION STAGE ROW IS A
      ! FRACTION OF [cm^-3], frozen with the rest of the background.  It is
      ! the cell quantity whose arithmetic face mean the stage flux divides
      ! by, taken from the one nucleus count of the species table
      ! (element_nucleus_counts), so the row's unknown, its time term and
      ! its face flux cannot disagree about how many nuclei a cell holds.
      real(dp), dimension(:), allocatable :: cbg_nHnuc, cbg_nHenuc
      ! THE HELIUM NUCLEI HELD OUTSIDE THE ATOMIC STAGES [cm^-3], frozen
      ! with the rest of the background: the helium of HeH+ and of every
      ! further helium-bearing molecule of the species table.  The stage
      ! rows partition the atomic helium alone -- a nucleus inside a
      ! molecule is moved by that molecule's own carrier row, and rescaling
      ! it would move a hydrogen nucleus with it -- so this density is
      ! reserved from the helium simplex, from the neutral closure of the
      ! chemistry rows and from the write-back's remainder alike.
      real(dp), dimension(:), allocatable :: cbg_nHemol
      ! THE FACE STATE THE IONIZATION STAGE ROWS RIDE ON, frozen over one
      ! interval exactly as the transport coefficients are: the element
      ! nucleus face flux of each element, its face nucleus density, the
      ! face eddy coefficient and the face spacing.  They come from
      ! hydrogen_and_helium_nucleus_face_flux, which reads the element
      ! operator's one public flux; nothing here rebuilds a face
      ! coefficient of that operator.
      real(dp), dimension(:,:), allocatable :: stg_Nel, stg_nelf
      real(dp), dimension(:),   allocatable :: stg_Kf, stg_drf
      ! The stage sum identity of the last stationary evaluation, the face
      ! that carries it, and whether any evaluation has measured it.
      ! One record per ELEMENT whose stages are carried (ien_H, ien_He):
      ! the identity is a statement about one element's own nucleus flux,
      ! so a single worst face over both elements could not say which
      ! element's construction it belongs to.
      real(dp), save :: stage_sum_max(2)    = 0.0d0
      integer,  save :: stage_sum_jworst(2) = 0
      logical,  save :: stage_sum_known(2)  = .false.
      ! The ratio of magnitude sums S/S' at that same face, which is what
      ! the rounding bound of the identity is proportional to
      ! (ionization_stage_transport, equation 5).  Negative means no face of
      ! the column carried a nonzero scale, so the identity was not measured
      ! anywhere and the bound has to fall back on its own ceiling.
      real(dp), save :: stage_sum_g(2)      = -1.0d0
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

      ! THE BASE FACE OF AN INFLOWING CARRIER.
      !
      ! The original design gave the carrier partition a ZERO-FLUX base: the
      ! element totals are boundary data the handoff supplies, but the split
      ! of an element among its carriers is what the oxygen option exists to
      ! COMPUTE, so imposing it at the boundary would have made that option's
      ! own A/B gate circular.  That argument holds exactly when nothing
      ! upstream states the partition.
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
      ! The two branches cannot both apply.  The code refuses q_H2_base
      ! while the oxygen chemistry is on, so a run either has a handoff
      ! partition (H2 only -- nothing states OH, H2O or CO) or computes its
      ! own (zero flux).
      !
      ! THE BASE BOUNDARY CONDITION OF CARRIER TRANSPORT, in one statement:
      ! ZERO DIFFUSIVE FLUX at face 0, and an advective flux whose
      ! composition is the ghost's wherever the face mass flux enters the
      ! domain.  The handoff supplies the composition of the gas that flows
      ! in, not a mixing rate across a boundary whose gradient is set by the
      ! ghost spacing and by a K_zz that describes an unresolved region, so
      ! the eddy and molecular terms are closed with no flux there and the
      ! budget entry of this boundary is the advective term alone.
      !
      ! IT IS NOT A FUNCTION OF THE ADVECTIVE CONTACT'S UPWIND CHOICE.  The
      ! direction this operator upwinds on is the sign of the FACE MASS FLUX
      ! the hydrodynamic solve returns, not the direction the base boundary
      ! condition upwinds the contact on; the two are separate statements and
      ! a change of one does not move the other.
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
      ! The first inequality says the destruction is fast enough that the
      ! transported CO reaches the
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
      ! (md/co_destruction_rates_literature_20260906.md sec. 12.5) puts
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
      ! THE LARGEST NEGATIVE GROUND-SINGLET HELIUM the write-back met in the
      ! state it was handed, relative to the helium nuclei of that cell.
      ! On the admissible set of the two ionized stages it is zero: the
      ! helium the stages did not take, less the frozen metastable, is
      ! nonnegative by construction (carrier_helium_available_to_stages).
      ! MEASURED rather than repaired, in the shape stage_fraction_under_zero
      ! has; a nonzero value says the state the step wrote holds more helium
      ! than its cells have.
      real(dp), save :: carrier_helium_singlet_under_zero = 0.0d0
      ! Molecular diffusion coefficient of each carrier on the grid [cm^2/s]
      ! at the state the run last held (carrier_diffusion_coefficient says
      ! which), so the run can print the transport time scale beside the
      ! chemical one it already prints.
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
      ! ---- THE TERMS OF EACH ROW, CELL BY CELL, FOR READING ONLY -------
      ! EXHALE_CARRIER_ROW_TERMS=1 records the individual terms of every
      ! carrier row of the last assembly, so that a row that refuses a
      ! certification can be read term by term instead of being inferred
      ! from its measure.  Default off, nothing in the solution reads these
      ! arrays, and they are not part of carrier_module_state: they are a
      ! record of a measurement, not a state the measurement may leave
      ! behind.
      real(dp), dimension(:,:), allocatable :: rowdump_dif   ! diffusive+drift divergence
      real(dp), dimension(:,:), allocatable :: rowdump_adv   ! material advective divergence
      real(dp), dimension(:,:), allocatable :: rowdump_prod  ! chemical production
      real(dp), dimension(:,:), allocatable :: rowdump_loss  ! chemical loss
      real(dp), dimension(:,:), allocatable :: rowdump_phot  ! the radiative part of the loss
      ! THE SCALE THE ACCEPTANCE DIVIDES BY, row_terms_phys of this row:
      ! the two face fluxes by their own MAGNITUDES and not by their
      ! difference, the advective magnitudes, the reaction terms and the
      ! absolute floor.  A record that rebuilt a scale from the DIVERGENCE
      ! instead understates it by the cancellation of the two faces, which
      ! in a smooth wind is two decades, so the number a reader compared
      ! with the certification would not be the certification's.
      real(dp), dimension(:,:), allocatable :: rowdump_scale
      real(dp), dimension(:,:), allocatable :: rowdump_res   ! the row itself
      real(dp), dimension(:,:), allocatable :: rowdump_nc    ! the carrier density [cm^-3]
      real(dp), dimension(:,:), allocatable :: rowdump_floor ! the absolute floor of the scale
      ! ---- THE H2 ROW, REACTION BY REACTION ---------------------------
      ! The H2 production and loss above are each a sum of a handful of
      ! reactions, and which of them carries the row is a physics question
      ! a net production and a net loss cannot answer: the base layer of a
      ! molecular case stands at a third of the lower-atmosphere handoff,
      ! and naming the channel that puts it there is what decides whether
      ! the network or the handoff is wrong.  These are the
      ! individual terms, each a volumetric rate [cm^-3 s^-1], written
      ! beside the H2 row of carrier_row_terms.txt.  The order is fixed
      ! here and mol_heh_rows fills 1 to n_h2chan_mol in it.
      integer, parameter :: n_h2chan_mol      = 13
      integer, parameter :: n_h2chan_oxy_prod = 14  ! OH + H -> H2 + O, H2O + H -> H2 + OH
      integer, parameter :: n_h2chan_oxy_loss = 15  ! H2 + O -> OH + H, H2 + OH -> H2O + H
      integer, parameter :: n_h2chan          = 15
      ! The names, in that order, for the header of the record.
      character(len=12), parameter :: h2chan_name(n_h2chan) =            &
         [ character(len=12) ::                                         &
           'R15_3body', 'R9_H2p_H', 'R6_H3p_e', 'R11_H3p_H',            &
           'P_H2_photo', 'LW_photdis', 'R10R13_Hp', 'R12_thermal',      &
           'R14_edis', 'R8_H2p_H2', 'R17R23_Hep', 'R18_HeHp',           &
           'HeI23S_H2', 'oxy_prod', 'oxy_loss' ]
      real(dp), dimension(:,:), allocatable :: rowdump_h2ch
      ! ---- THE FOUR FACE FLUXES OF EACH ROW, SEPARATELY AND SIGNED -----
      ! The divergences above answer "how much did transport move"; they
      ! cannot answer "which face, and in which direction", and that is the
      ! question a residual sitting at a boundary cell asks.  These are the
      ! diffusive and the advective flux at the INNER face (j-1) and at the
      ! OUTER face (j) of the cell, in the row's own units, so that an
      ! advective flux DIFFERENCE, a diffusive flux that CHANGES SIGN across
      ! the cell and a chemical imbalance are told apart by reading them
      ! rather than by inference.  Diagnostic only, filled with the rest of
      ! the record and written beside it.
      real(dp), dimension(:,:), allocatable :: rowdump_fdif_in
      real(dp), dimension(:,:), allocatable :: rowdump_fdif_out
      real(dp), dimension(:,:), allocatable :: rowdump_fadv_in
      real(dp), dimension(:,:), allocatable :: rowdump_fadv_out
      real(dp), dimension(:),   allocatable :: rowdump_frho_in
      real(dp), dimension(:),   allocatable :: rowdump_frho_out
      logical, save :: rowdump_on     = .false.
      logical, save :: rowdump_asked  = .false.
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
      type :: carrier_module_state
         real(dp), allocatable :: cbg_nhii(:), cbg_nh2p(:), cbg_nh3p(:)
         real(dp), allocatable :: cbg_nhehp(:), cbg_nhei(:), cbg_nheii(:)
         real(dp), allocatable :: cbg_nheiii(:), cbg_nheiTR(:), cbg_ne(:)
         real(dp), allocatable :: cbg_nOion(:), cbg_nCion(:)
         real(dp), allocatable :: cbg_nm0(:,:), cbg_nm1(:,:), cbg_nm2(:,:)
         real(dp), allocatable :: cbg_nHnuc(:), cbg_nHenuc(:)
         real(dp), allocatable :: cbg_nHemol(:)
         real(dp), allocatable :: stg_Nel(:,:), stg_nelf(:,:)
         real(dp), allocatable :: stg_Kf(:), stg_drf(:)
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
      ! path writes either.  They exist because the branches they reach --
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
      type :: carrier_checkpoint
         real(dp), allocatable :: fc(:,:)
         ! The arrays a residual assembly overwrites, in the enumeration
         ! save_carrier_module_state already keeps, reused rather than
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

      ! THE THERMOCHEMICAL STATE ONE TRIAL OF THE RELAXATION WRITES, held
      ! aside and put back when the trial is undone.
      !
      ! A trial takes a transport step and then closes the chemistry at the
      ! fixed conserved state, and the chemistry sweep writes far more than
      ! the composition it is handed: the frozen cell state the next sweep
      ! and the transport operator read, the rate state the closure
      ! evaluator measures a composition against, the composition the
      ! caloric maps carry, and the base ghost the boundary condition
      ! reads.  A trial that is refused has to put ALL of them back, or the
      ! state the caller receives is a composition of one trial beside the
      ! rates, the heat capacity and the boundary of another.  MEASURED on
      ! the carrier_retry column before this was enumerated: the rate state
      ! was left 14.87 K from the composition, so the closure residual of
      ! the returned state read 5.3e-4 where the state's own is 1.6e-16,
      ! and the energy-to-pressure map stood 1.5e-3 relative away from the
      ! entry one.
      !
      ! THE NON-ROOT PERSISTENCE COUNTER IS PART OF IT, and putting it back
      ! is what the counter means: it counts the CONSECUTIVE sweeps a cell
      ! of the accepted state has rested on a non-root acceptance, which is
      ! what separates a recovering transient from a wind built on one
      ! (ionization_equilibrium, ieq_nonroot_streak_stop).  The sweeps of a
      ! trial the caller discarded are not sweeps of the state the run
      ! carries, so they do not count towards that persistence.
      !
      ! THE CALORIC STATE IS NOT COPIED, it is rebuilt: the arrays are
      ! private to caloric_eos and its header states the rule -- a path
      ! that installs a different composition must refresh them BEFORE it
      ! maps an energy density to a pressure -- so the restore goes through
      ! get_species_densities, the one routine that turns (rho, f_sp) into
      ! number densities and refreshes them, exactly as the production
      ! paths do.
      type :: thermochemical_state
         real(dp), allocatable :: f_sp(:,:)
         real(dp), allocatable :: p(:), T(:)
         real(dp), allocatable :: heat(:), cool(:), eta(:)
         type(ion_rates), allocatable :: bg_cell(:)
         type(ieq_rate_state) :: rates
         real(dp) :: dp_bc = 0.0d0
         real(dp) :: n_part_cell1 = 0.0d0
         type(carrier_checkpoint) :: carriers
      end type thermochemical_state


      contains

      ! Fix the transported set for this run.  Called once from input_read,
      ! after every key is parsed, so that no routine can see a half-decided
      ! carrier set.
      subroutine carrier_set_init()
      carrier_solved = .false.
      ! THE MOLECULAR CARRIERS, where the network exists and its transport
      ! is selected.  An atomic gas has none of them: there is no H2 row to
      ! solve in a gas with no H2.
      if (thereis_mol .and. carrier_transport) then
         carrier_solved(ic_H2) = .true.
         if (thereis_oxychem) then
            carrier_solved(ic_OH)  = .true.
            carrier_solved(ic_H2O) = .true.
            carrier_solved(ic_CO)  = .true.
         endif
      endif
      ! THE CARRIED STAGES OF THE TWO ELEMENTS, together.  The key states
      ! that the ionization state of the gas is carried by the flow, and
      ! the ionization state of a helium/hydrogen mixture is the partition
      ! of BOTH elements: carrying the proton while helium is re-solved on
      ! its local root every sweep would hand the electron budget of the
      ! transported hydrogen to a helium partition that never saw the flow.
      ! The helium stages need helium in the mixture, which the key
      ! requires (input_read).
      if (ionization_transport) then
         carrier_solved(ic_Hp) = .true.
         if (thereis_He) then
            carrier_solved(ic_HeII)  = .true.
            carrier_solved(ic_HeIII) = .true.
         endif
      endif
      ! THE He 2^3S POPULATION, as a fourth helium state of the same
      ! partition: it rides on the helium nucleus flux the two ionized
      ! stages ride on, and its ground singlet closes against them, so it
      ! is carried only with them (input_read refuses the key otherwise).
      if (he23s_transport .and. ionization_transport .and. thereis_He    &
          .and. thereis_HeITR) carrier_solved(ic_HeTR) = .true.
      ! n_carrier is the last index in the solved set: the loops run
      ! 1..n_carrier and skip the ones carrier_solved leaves out.
      if (carrier_solved(ic_HeTR)) then
         n_carrier = ic_HeTR            ! ..., H+, He+, He++, He 2^3S
      else if (ionization_transport) then
         n_carrier = ic_HeIII           ! [H2 [, OH, H2O, CO],] H+, He+, He++
      else if (thereis_oxychem) then
         n_carrier = ic_CO              ! H2, OH, H2O, CO
      else
         n_carrier = ic_H2              ! H2 alone
      endif

      ! The H2 thermochemistry table is built here, serially, if no caller
      ! built it yet: the carrier source terms read it inside parallel
      ! regions, and a lazy build there was removed on 2026-09-13.
      ! The main program initializes it at startup; this
      ! covers a test driver that enters through the carriers.
      if (thereis_mol .and. .not. h2_thermochemistry_ready())             &
         call h2_thermochemistry_init

      ! THE ENTHALPY THE CARRIERS CARRY (the molecular carriers relative
      ! to their elements, the carried ionization stages relative to their
      ! element) enters the stationary energy row
      ! through the one procedure the residual module holds for it (that
      ! module cannot use this one, which uses it): registered here, where
      ! the carrier set is decided, and read by assemble_residual.  A run
      ! whose energy equation carries no carrier term leaves it null, and
      ! its residual is the one it had.
      if (carrier_enthalpy_active()) then
         carrier_enthalpy_divergence => carrier_enthalpy_divergence_of_state
      else
         nullify(carrier_enthalpy_divergence)
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
      ! zero-gradient condition written on the face the
      ! gas crosses.  The proton is never stated: no key gives an ionization
      ! fraction below the base, so it keeps the zero-gradient condition in
      ! every configuration.
      call advected_carrier_reset()
      if (transported_rows_exist()) then
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
         if (carrier_solved(ic_HeII))                                    &
            call advected_carrier_register(isp_HeII, .false.)
         if (carrier_solved(ic_HeIII))                                   &
            call advected_carrier_register(isp_HeIII, .false.)
         ! Registered last, so the register order stays the order of the
         ! carrier indices (carrier_state reads the advected columns in
         ! that order).  No handoff states a 2^3S fraction below the base:
         ! the inner ghosts take the base cell's own, as for the stages.
         if (carrier_solved(ic_HeTR))                                    &
            call advected_carrier_register(isp_HeTR, .false.)
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
      ! that cell, so the attribution is per cell.  Those rows are counted
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

      ! Molecular diffusion coefficient of carrier ic in cell j [cm^2/s] of
      ! the state the run last HELD, which is the last of
      !   * the entry state of an accepted carrier transport interval
      !     (carrier_transport_interval: a marching step, or a kept trial of
      !     the carrier relaxation; a refused step or trial puts the value
      !     of the state before it back through the carrier checkpoint),
      !   * the state the carrier relaxation returns
      !     (relax_photochemical_composition), and
      !   * the state the steady solver hands back (solve_steady_ptc,
      !     solve_steady_jfnk, through carrier_diffusivities_of_state), and
      !   * the state the chemical-decay diagnostic of the steady solver
      !     reports on (chemical_decay_against_the_transport_times).
      ! A residual evaluated at a trial or probe state stores nothing.  Zero
      ! before the first of them.
      double precision function carrier_diffusion_coefficient(j, ic)      &
                                result(D)
      integer, intent(in) :: j, ic
      if (allocated(pct_Dco)) then
         D = pct_Dco(j,ic)
      else
         D = 0.0d0
      endif
      end function carrier_diffusion_coefficient

      ! The molecular diffusivities of the state (rho, f_sp) the run holds,
      ! stored for carrier_diffusion_coefficient.  The temperature and the
      ! gas-particle density are those of the frozen background, which the
      ! caller has made the background of that state (for the steady solver,
      ! install_background_of_adopted_state); they are the same numbers
      ! carrier_state hands the transport operator.  Nothing is stored where
      ! no transported balance exists or no background has been formed.
      subroutine carrier_diffusivities_of_state(rho, f_sp)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), dimension(1-Ng:N+Ng) :: ntot_bg, TK_bg
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: D_state
      integer :: j
      if (.not. transported_rows_exist() .or. .not. bg_ready) return
      if (.not. allocated(bg_cell)) return
      do j = 1-Ng, N+Ng
         ntot_bg(j) = bg_cell(j)%ntot
         TK_bg(j)   = bg_cell(j)%T_K
      enddo
      call carrier_diffusivities(f_sp, rho, TK_bg, ntot_bg, D_state)
      if (.not. allocated(pct_Dco))                                      &
         allocate(pct_Dco(1-Ng:N+Ng,n_carrier_max))
      pct_Dco = D_state
      end subroutine carrier_diffusivities_of_state

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
           carrier_exhausted_step_first = marching_step
      carrier_exhausted_step_last = marching_step
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
            ' (photochemical_transport_step) step ', marching_step,               &
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
      ! the subdivision.  The time-integration error is measured by step
      ! doubling and printed.
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
      if (.not. transported_rows_exist()) return
      if (.not. bg_ready)                return
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
      call carrier_stage_face_state(rho, TK, f_sp, carrier_rows_advect,  &
                                    msum, Frho)
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
      ! EXHALE_CARRIER_ROW_TERMS=1 records the terms of every carrier row
      ! and writes them once per certification to
      ! output/carrier_row_terms.txt.  Default off; nothing in the solution
      ! reads the record, so the run it is switched on in is the run it is
      ! switched off in.
      logical function carrier_row_terms_on()
      character(len=8) :: env
      call get_environment_variable('EXHALE_CARRIER_ROW_TERMS', env)
      carrier_row_terms_on = (trim(env) .eq. '1')
      end function carrier_row_terms_on

      ! ------------------------------------------------------------- !

      ! THE CARRIER ROWS OF THE LAST ASSEMBLY, TERM BY TERM.  One line per
      ! (cell, solved carrier): the transport terms, the chemical
      ! production and loss separately, the radiative part of the loss, the
      ! row itself, and the row measured on ONE scale, rowdump_scale: the
      ! row's own physical scale that the acceptance and the certification
      ! divide by (row_terms_phys), the magnitudes of the two face fluxes and
      ! of the advective terms, the reaction terms and the 1e-20
      ! free-element floor.  The file header states the same.  For a helium
      ! stage the production and loss columns are the positive and negative
      ! parts of the net source, not gross reaction rates.
      subroutine carrier_row_terms_write(fname)
      character(len=*), intent(in) :: fname
      integer :: u, j, ic, k
      real(dp) :: dif, advj, pr, ls, fl, rr, sterm
      if (.not. rowdump_on)            return
      if (.not. allocated(rowdump_res)) return
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u,'(A)') '# carrier row terms of the last assembly '//       &
                     '(EXHALE_CARRIER_ROW_TERMS=1)'
      write(u,'(A)') '# every rate is a volumetric rate [cm^-3 s^-1]; '//&
                     'the row is transport - (production - loss)'
      write(u,'(A)') '# columns: cell r[R_p] T[K] carrier n_c[cm^-3] '// &
                     'x_c diffusive advective production loss '//        &
                     'photo_loss net_source residual floor scale measure'
      write(u,'(A)') '# scale is the row''s own physical scale, the one '//&
                     'the acceptance and the certification divide by '// &
                     '(row_terms_phys): the two face fluxes by their '// &
                     'magnitudes, the advective magnitudes, the reaction'
      write(u,'(A)') '# terms and the 1e-20 free-element floor; measure '//&
                     '= |residual|/scale.  The diffusive column beside '//&
                     'it is the DIVERGENCE, which the two faces cancel '//&
                     'into and which is therefore not that scale.'
      ! The H2 row carries its chemistry reaction by reaction as well, in
      ! the fixed order of h2chan_name; every other carrier writes zeros
      ! there.  Channels 1 to 4 are formation, 5 to 13 destruction, and the
      ! last two the oxygen chemistry's two sides when it is on.
      write(u,'(A)') '# and for the H2 row, the same production and '//  &
                     'loss reaction by reaction [cm^-3 s^-1]:'
      write(u,'(A)') '# '//trim(h2chan_header())
      write(u,'(A)') '# and, for every row, the TWO FACES separately: '// &
                     'each term already carries its area and the cell '// &
                     'volume, so fdif_in + fdif_out = diffusive and '//   &
                     'fadv_in + fadv_out = advective.'
      write(u,'(A)') '# face cell carrier fdif_in fdif_out fadv_in '//    &
                     'fadv_out Frho_in Frho_out'
      do j = 1, N
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            dif  = rowdump_dif(j,ic)
            advj = rowdump_adv(j,ic)
            pr   = rowdump_prod(j,ic)
            ls   = rowdump_loss(j,ic)
            fl   = rowdump_floor(j,ic)
            rr   = rowdump_res(j,ic)
            sterm = rowdump_scale(j,ic)
            write(u,'(I6,1X,ES14.7,1X,ES12.5,1X,A8,12(1X,ES15.7E3))')    &
                 j, r(j), bg_cell(j)%T_K, trim(carrier_name(ic)),        &
                 rowdump_nc(j,ic),                                       &
                 rowdump_nc(j,ic)/max(bg_cell(j)%ntot, 1.0d-300),        &
                 dif, advj, pr, ls, rowdump_phot(j,ic), pr - ls, rr, fl, &
                 sterm, abs(rr)/max(sterm, 1.0d-300)
            write(u,'(A,I6,1X,A8,6(1X,ES15.7E3))') '  face', j,        &
                 trim(carrier_name(ic)),                                 &
                 rowdump_fdif_in(j,ic), rowdump_fdif_out(j,ic),          &
                 rowdump_fadv_in(j,ic), rowdump_fadv_out(j,ic),          &
                 rowdump_frho_in(j), rowdump_frho_out(j)
            if (ic .eq. ic_H2) then
               write(u,'(A,I6,15(1X,ES15.7E3))') '  chan', j,            &
                    (rowdump_h2ch(j,k), k = 1, n_h2chan)
            endif
         enddo
      enddo
      close(u)
      end subroutine carrier_row_terms_write

      ! The channel names of the H2 record as one line, for its header.
      function h2chan_header() result(hdr)
      character(len=16*n_h2chan+8) :: hdr
      integer :: k
      hdr = 'chan cell'
      do k = 1, n_h2chan
         hdr = trim(hdr)//' '//trim(h2chan_name(k))
      enddo
      end function h2chan_header

      ! ------------------------------------------------------------- !

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
      ! nH_free is the cell's hydrogen less the nucleus locked in HeH+, and
      ! nO_free and nC_free are its oxygen and carbon totals whatever species
      ! holds them (carrier_element_totals in element_inventory).
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
      real(dp), dimension(1-Ng:N+Ng) :: nd, nucH, nucHe
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
         allocate(cbg_nm0(1-Ng:N+Ng,n_melem),                            &
                  cbg_nm1(1-Ng:N+Ng,n_melem),                            &
                  cbg_nm2(1-Ng:N+Ng,n_melem))
         cbg_nm0 = 0.0d0
         cbg_nm1 = 0.0d0
         cbg_nm2 = 0.0d0
         allocate(cbg_nHnuc(1-Ng:N+Ng), cbg_nHenuc(1-Ng:N+Ng))
         allocate(cbg_nHemol(1-Ng:N+Ng))
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
      ! The reactant densities of the metal charge exchange, stage by
      ! stage.  melem_i0(e) is the neutral's column in the ion table and
      ! the singly and doubly ionized stages follow it; an element carried
      ! to the singly ionized stage alone (melem_top = 1) keeps a zero
      ! doubly ionized density, and every reaction that would read it has
      ! a zero rate.
      cbg_nm0 = 0.0d0
      cbg_nm1 = 0.0d0
      cbg_nm2 = 0.0d0
      do im = 1, n_melem
         i0 = melem_i0(im)
         cbg_nm0(:,im) = f_sp(:,mion_fsp(i0))*nd
         cbg_nm1(:,im) = f_sp(:,mion_fsp(i0+1))*nd
         if (melem_top(im) .ge. 2)                                       &
            cbg_nm2(:,im) = f_sp(:,mion_fsp(i0+2))*nd
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
      ! THE ELEMENT NUCLEUS DENSITIES THE IONIZATION STAGE ROWS ARE
      ! FRACTIONS OF, from the same count whose face mean the stage flux
      ! divides by.  They are part of the frozen background of the step:
      ! transport moves a nucleus from cell to cell, and which stage it
      ! sits in is what the stage row solves, so the element total of a
      ! cell does not move while the rows are being solved.
      call element_nucleus_counts(f_sp, nucH, nucHe)
      cbg_nHnuc  = nucH *nd
      cbg_nHenuc = nucHe*nd
      ! THE HELIUM THE ATOMIC STAGES DO NOT PARTITION, from the same table:
      ! every species that carries a helium nucleus and is not one of the
      ! three atomic stages, which today is HeH+ alone.  The He 2^3S column
      ! is a level inside He I and its nucleus is already counted through
      ! He I, so it is skipped here exactly as it is in the nucleus count.
      cbg_nHemol = 0.0d0
      do k = 1, n_bsp
         if (bsp_is_excited_level(k)) cycle
         if (bsp_nHe(k) .le. 0) cycle
         if (bsp_fsp(k) .eq. isp_HeI  .or. bsp_fsp(k) .eq. isp_HeII .or. &
             bsp_fsp(k) .eq. isp_HeIII) cycle
         cbg_nHemol = cbg_nHemol + dble(bsp_nHe(k))*f_sp(:,bsp_fsp(k))*nd
      enddo
      ! The proton column is filled whether or not it is solved -- it costs
      ! one assignment and it keeps fc a complete picture of the state -- but
      ! only a run with ionization_transport lets the solve move it.  Its
      ! unknown is the fraction of the HYDROGEN NUCLEI in H II, which is
      ! the variable whose flux the element's own nucleus flux carries
      ! (carrier_is_nucleus_fraction).
      fc(:,ic_Hp)  = f_sp(:,isp_HII)*nd/max(cbg_nHnuc, 1.0d-300)
      ! The two ionized helium stages, per HELIUM nucleus, from the same
      ! count.  Filled whether or not they are solved, for the reason the
      ! proton column above is.
      fc(:,ic_HeII)  = f_sp(:,isp_HeII) *nd/max(cbg_nHenuc, 1.0d-300)
      fc(:,ic_HeIII) = f_sp(:,isp_HeIII)*nd/max(cbg_nHenuc, 1.0d-300)
      ! The He 2^3S level, per helium nucleus, from the same count.
      fc(:,ic_HeTR) = 0.0d0
      if (thereis_HeITR)                                                 &
         fc(:,ic_HeTR) = f_sp(:,isp_HeTR)*nd/max(cbg_nHenuc, 1.0d-300)

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
            ! The stages advect a fraction per unit MASS, which is what
            ! every species row of the hydrodynamic operator carries.  An
            ! ionization stage's unknown is a fraction per element NUCLEUS,
            ! so the advected value is divided by the element's nucleus
            ! count here, at the same place the entry value above was.
            !
            ! THE NUCLEUS COUNT IS THE ONE THE SPECIES VECTOR HOLDS, which
            ! at this point is the composition BEFORE the stages moved the
            ! declared set (the paragraph above).  The stages move the
            ! stage and the molecular carriers but not the neutral stage
            ! that closes the element, so the count is short by what they
            ! moved, of order the advective change of one step.  It is the
            ! same reference the chemistry rows close atomic hydrogen with
            ! (nH_free, taken from this same f_sp), so the two agree; a
            ! nucleus count taken from the advected columns alone would
            ! not, because those columns are not the whole element.
            if (carrier_is_nucleus_fraction(ic)) then
               do j = 1-Ng, N+Ng
                  fc(j,ic) = fcadv(j,k)*nd(j)                            &
                       /max(carrier_stage_nucleus_density(ic, j),        &
                            1.0d-300)
               enddo
            else
               fc(:,ic) = fcadv(:,k)
            endif
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
      if (allocated(cbg_nHnuc))    s%cbg_nHnuc    = cbg_nHnuc
      if (allocated(cbg_nHenuc))   s%cbg_nHenuc   = cbg_nHenuc
      if (allocated(cbg_nHemol))   s%cbg_nHemol   = cbg_nHemol
      if (allocated(stg_Nel))      s%stg_Nel      = stg_Nel
      if (allocated(stg_nelf))     s%stg_nelf     = stg_nelf
      if (allocated(stg_Kf))       s%stg_Kf       = stg_Kf
      if (allocated(stg_drf))      s%stg_drf      = stg_drf
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
      if (allocated(cbg_nm0))      s%cbg_nm0      = cbg_nm0
      if (allocated(cbg_nm1))      s%cbg_nm1      = cbg_nm1
      if (allocated(cbg_nm2))      s%cbg_nm2      = cbg_nm2
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
      call put_back_1d(s%cbg_nHnuc,    cbg_nHnuc)
      call put_back_1d(s%cbg_nHenuc,   cbg_nHenuc)
      call put_back_1d(s%cbg_nHemol,   cbg_nHemol)
      call put_back_2d(s%stg_Nel,      stg_Nel)
      call put_back_2d(s%stg_nelf,     stg_nelf)
      call put_back_1d(s%stg_Kf,       stg_Kf)
      call put_back_1d(s%stg_drf,      stg_drf)
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
      call put_back_2d(s%cbg_nm0,      cbg_nm0)
      call put_back_2d(s%cbg_nm1,      cbg_nm1)
      call put_back_2d(s%cbg_nm2,      cbg_nm2)
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
      ok = ok .and. same_1d(s%cbg_nHnuc,    cbg_nHnuc)
      ok = ok .and. same_1d(s%cbg_nHenuc,   cbg_nHenuc)
      ok = ok .and. same_1d(s%cbg_nHemol,   cbg_nHemol)
      ok = ok .and. same_2d(s%stg_Nel,      stg_Nel)
      ok = ok .and. same_2d(s%stg_nelf,     stg_nelf)
      ok = ok .and. same_1d(s%stg_Kf,       stg_Kf)
      ok = ok .and. same_1d(s%stg_drf,      stg_drf)
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
      ! arrays a residual assembly overwrites (through the same pair,
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

      ! ------------------------------------------------------------- !

      subroutine save_thermochemical_state(s, f_sp, p, T, heat, cool, eta)
      ! Hold aside everything one trial of the relaxation writes (the
      ! type's header says what and why).
      type(thermochemical_state), intent(out) :: s
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: p, T, heat, cool, eta
      s%f_sp = f_sp
      s%p    = p
      s%T    = T
      s%heat = heat
      s%cool = cool
      s%eta  = eta
      if (allocated(bg_cell)) s%bg_cell = bg_cell
      call save_ieq_rate_state(s%rates)
      s%dp_bc        = dp_bc
      s%n_part_cell1 = n_part_cell1
      call carrier_checkpoint_take(s%carriers)
      end subroutine save_thermochemical_state

      ! ------------------------------------------------------------- !

      subroutine restore_thermochemical_state(s, rho, f_sp, p, T, heat,   &
                                              cool, eta)
      ! Put the held state back.  The composition is reinstated first and
      ! the caloric maps are refreshed from it through the production
      ! installation path, so that the pressure the next energy-to-pressure
      ! map returns is the one of the composition standing beside it.
      type(thermochemical_state), intent(in) :: s
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: p, T, heat, cool, eta

      real(dp), dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii
      real(dp), dimension(1-Ng:N+Ng) :: nheiTR, ne, ntot
      real(dp), dimension(1-Ng:N+Ng,n_mion) :: nm

      f_sp = s%f_sp
      p    = s%p
      T    = s%T
      heat = s%heat
      cool = s%cool
      eta  = s%eta
      if (allocated(s%bg_cell)) bg_cell = s%bg_cell
      call restore_ieq_rate_state(s%rates)
      call carrier_checkpoint_restore(s%carriers)
      nhei   = 0.0d0
      nheii  = 0.0d0
      nheiii = 0.0d0
      nheiTR = 0.0d0
      call get_species_densities(rho, f_sp, nhi, nhii, nhei, nheii,       &
                                 nheiii, nheiTR, nm, ne, ntot)
      ! The base ghost is written back AFTER the refresh, because
      ! get_species_densities writes the particle count of the first cell
      ! itself: what the boundary condition read at the entry of the trial
      ! is the held pair, not the one a re-derivation happens to give.
      dp_bc        = s%dp_bc
      n_part_cell1 = s%n_part_cell1
      end subroutine restore_thermochemical_state

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
      ! Restoration is a measurement over the whole list, and not over the
      ! items one particular column happens to move: an item left out of
      ! carrier_checkpoint_restore is then a failing assertion.  An array
      ! that is not allocated is allocated here, so the restore is asked to
      ! give the module its allocation status back and not only its
      ! numbers.  Test only: no production path calls this.
      subroutine carrier_perturb_checkpointed_state_for_test(fc)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout),        &
                optional :: fc
      if (present(fc)) fc = fc + 1.0d0
      call bump_1d(cbg_nhii,   1-Ng, N+Ng)
      call bump_1d(cbg_nHnuc,  1-Ng, N+Ng)
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
                                         res_out, terms_out, absent_out)
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
      ! The cells in which a carrier stands below carrier_absent_fraction of
      ! its own element's free reservoir (carrier_residual): a row the
      ! certification reports and does not gate.
      logical,  dimension(1:N,n_carrier_max), intent(out), optional ::   &
                                                 absent_out

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
      ! A carrier a run does not solve, or a run that reaches none of this,
      ! leaves every cell judged: absence is asserted, never assumed.
      if (present(absent_out)) absent_out = .false.
      if (.not. transported_rows_exist()) return
      if (.not. bg_ready)                return

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
      call carrier_stage_face_state(rho, TK, f_sp, .true., msum, Frho)
      ! The stage sum identity of THIS state, measured on the same face
      ! fluxes the rows below are the divergence of.
      if (carrier_solved(ic_Hp)) then
         call ionization_stage_nucleus_sum(fc, ien_H,                     &
                                           stage_sum_max(ien_H),          &
                                           stage_sum_jworst(ien_H),       &
                                           stage_sum_g(ien_H))
         stage_sum_known(ien_H) = .true.
      endif
      if (carrier_solved(ic_HeII) .and. carrier_solved(ic_HeIII)) then
         call ionization_stage_nucleus_sum(fc, ien_He,                    &
                                           stage_sum_max(ien_He),         &
                                           stage_sum_jworst(ien_He),      &
                                           stage_sum_g(ien_He))
         stage_sum_known(ien_He) = .true.
      endif
      call carrier_photolysis(rho, TK, f_sp)
      call carrier_signal_rate(v, TK, mbar, sigrate)
      if (present(absent_out)) then
         call carrier_residual(fc, fc, nrho, wfac, dt_big, rp, rep, msum, &
                               Frho, .true.,                              &
                               Agrd, Bdrf, updrf, nH_free, nO_free,       &
                               nC_free, sigrate, res, Jf, dJl, dJr,       &
                               rnorm_unused, absent_out = absent_out)
      else
         call carrier_residual(fc, fc, nrho, wfac, dt_big, rp, rep, msum, &
                               Frho, .true.,                              &
                               Agrd, Bdrf, updrf, nH_free, nO_free,       &
                               nC_free, sigrate, res, Jf, dJl, dJr,       &
                               rnorm_unused)
      endif

      ! THE SCALE IS THE ROW'S OWN LARGEST TERMS, and that is a change.
      ! It used to be carrier_row_scale, n_H(|v|+c_s)/dr -- the flux
      ! divergence the cell would carry if its WHOLE hydrogen inventory moved
      ! at the signal speed. That is a bound on the row, but it is not the
      ! row's own term, and measured on the He/H = 0.0793 hot Uranus it
      ! exceeds the sum of the terms the row actually balances by 10^2.7 to
      ! 10^6.4. A state the steady solver accepted at info = 0 then reported
      ! a carrier
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
            col_scale_car(j,ic) = fc(j,ic)                               &
                            *carrier_nucleus_reference(ic, j, nrho(j))   &
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
                                          nO_free(j), nC_free(j),         &
                                          cbg_nHnuc(j), cbg_nHenuc(j))
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
      ! property of a state this module has assembled, and the projection
      ! onto the budget face of the coupled solve's unknown box is taken on
      ! a SYNTHETIC cell, which has no such state.  It writes the budget
      ! [cm^-3] the face is then
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
      !
      ! AN IONIZATION STAGE IS BOUNDED BY THE SIMPLEX AND NOT BY THIS
      ! BUDGET.  Its unknown is a fraction of its element's NUCLEI, so the
      ! largest density it can reach is the element's own nucleus density
      ! nH_nuc -- every nucleus of the element in that stage -- and not the
      ! smaller reservoir hydrogen_available_to_carriers leaves the
      ! molecular carriers.  The two differ because the molecular carriers
      ! take nuclei the frozen ion stages hold, while a stage row
      ! RE-PARTITIONS the element it belongs to and takes nothing from any
      ! other element.
      double precision function carrier_element_headroom(ic, nH_avail,    &
                               nO_free, nC_free, nH_nuc, nHe_nuc)        &
                               result(nmax)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nH_avail, nO_free, nC_free, nH_nuc, nHe_nuc
      select case (ic)
      case (ic_H2)
         nmax = element_box_side(isp_H2, ien_H, nH_avail)
      case (ic_Hp)
         nmax = element_box_side(isp_HII, ien_H, nH_nuc)
      case (ic_HeII)
         nmax = element_box_side(isp_HeII, ien_He, nHe_nuc)
      case (ic_HeIII)
         nmax = element_box_side(isp_HeIII, ien_He, nHe_nuc)
      case (ic_HeTR)
         ! One He 2^3S atom holds one helium nucleus, the nucleus of the
         ! He I atom it is a level of.  The element inventory counts that
         ! nucleus once, through He I, and gives the level none of its own
         ! (bsp_is_excited_level), so element_box_side refuses the level;
         ! its side of the box is the helium nucleus density itself.
         nmax = nHe_nuc
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
      !     adv(j) = ( A_+ F_c(j) - A_- F_c(j-1) ) / V_j  x  msum/m_c
      !
      ! with F_rho the face mass flux of the mass row of this state,
      ! Y_c^face the reconstructed and upwinded face mass fraction of
      ! species_face_fraction, and A, V the face area and the exact shell
      ! volume of spherical_face_area_and_cell_volume, which the diffusive
      ! half of the same row divides by as well.  The factor msum/m_c returns the mass divergence
      ! to the units the carrier row is written in, n_c per unit volume and
      ! time [cm^-3 s^-1], and it is the inverse of the conversion
      ! carrier_mass_fractions applies on the way in.
      !
      ! WHY IT IS NOT A CELL-VELOCITY UPWIND DIFFERENCE ANY MORE.  Earlier
      ! the stationary rows and the fixed-wind relaxation carried
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
      subroutine carrier_advective_divergence(fc, msum, Frho, adv, advmag, &
                                              advin, advout)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: msum
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: Frho
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: adv
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: advmag
      ! THE TWO FACES SEPARATELY, each already carrying its area and the
      ! cell volume, so advin + advout = adv exactly.  The divergence says
      ! how much transport moved; these say through WHICH face and in which
      ! direction, which is the question a residual sitting at a boundary
      ! cell asks.  Diagnostic and optional; nothing in the solution asks
      ! for them.
      real(dp), dimension(1:N,n_carrier_max), optional, intent(out) :: advin
      real(dp), dimension(1:N,n_carrier_max), optional, intent(out) :: advout

      real(dp), dimension(1-Ng:N+Ng) :: Y, Yf, Fs, dvF, dvM
      real(dp), dimension(0:N) :: fa
      real(dp), dimension(1:N) :: cv
      real(dp) :: mc, tconv
      integer  :: j, ic

      adv    = 0.0d0
      advmag = 0.0d0
      if (present(advin))  advin  = 0.0d0
      if (present(advout)) advout = 0.0d0
      call l22b_setup()
      call spherical_face_area_and_cell_volume(fa, cv)
      ! Code density times code velocity over code length is a rate in units
      ! of v0/R0, and f_sp counts n0 particles per unit density, so this is
      ! the factor that puts the divergence in cm^-3 s^-1 -- the same factor
      ! the time term of the row carries through nrho and dt_phys.
      tconv = n0*v0/R0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         ! An ionization stage's material advection is x_k N_el, the first
         ! term of the stage flux, and is already inside the face flux that
         ! row is the divergence of; taking it a second time on the bulk
         ! face mass flux would advect the stage twice and on the wrong
         ! flux (carrier_is_nucleus_fraction).
         if (carrier_is_nucleus_fraction(ic)) cycle
         mc = carrier_mass_amu(ic)
         call carrier_face_mass_fraction(fc, msum, ic, Y)
         call species_face_fraction(Y, Frho, Yf)
         ! Y already carries the stated ghosts at both ends
         ! (carrier_face_mass_fraction), so the reconstruction below reads
         ! one boundary rule wherever this term is evaluated.  With the
         ! outflow continuation X_ghost = X_N the forward difference
         ! entering the MC limiter of PLM is exactly zero, the limited slope
         ! of cell N is zero and the outflow face is the donor cell average
         ! by itself.
         call species_face_flux(Frho, Yf, Fs)
         do j = 1, N
            dvF(j)       = (fa(j)*Fs(j) - fa(j-1)*Fs(j-1))/cv(j)
            dvM(j)       = (fa(j)*abs(Fs(j))                             &
                            + fa(j-1)*abs(Fs(j-1)))/cv(j)
            adv(j,ic)    = dvF(j)*msum(j)/mc*tconv
            advmag(j,ic) = dvM(j)*msum(j)/mc*tconv
            if (present(advin))                                          &
               advin(j,ic)  = -fa(j-1)*Fs(j-1)/cv(j)*msum(j)/mc*tconv
            if (present(advout))                                         &
               advout(j,ic) =  fa(j)  *Fs(j)  /cv(j)*msum(j)/mc*tconv
         enddo
      enddo

      end subroutine carrier_advective_divergence

      ! ------------------------------------------------------------- !

      ! The two face coefficients of the advective term, in the row's own
      ! units per unit of the face MASS fraction:
      !
      !     advj(j) =  A_+ F_rho(j)  /V_j ,
      !     advm(j) = -A_- F_rho(j-1)/V_j ,
      !
      ! the same A_+, A_- and V_j the carrier residual divides by.  They are
      ! what the block-tridiagonal Jacobian of the fixed-wind relaxation
      ! differentiates the term with: multiplied by msum(j)/m_c they give
      ! d adv(j) / d Y(donor), and divided once more by the donor cell's own
      ! mixture mass they give d adv(j) / d f_c(donor).
      subroutine carrier_advective_face_coefficients(Frho, advj, advm)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: Frho
      real(dp), dimension(1:N),       intent(out) :: advj, advm
      real(dp), dimension(0:N) :: fa
      real(dp), dimension(1:N) :: cv
      real(dp) :: tconv
      integer  :: j
      call spherical_face_area_and_cell_volume(fa, cv)
      tconv = n0*v0/R0
      do j = 1, N
         advj(j) =  fa(j)  *Frho(j)  /cv(j)*tconv
         advm(j) = -fa(j-1)*Frho(j-1)/cv(j)*tconv
      enddo
      end subroutine carrier_advective_face_coefficients

      ! ------------------------------------------------------------- !

      ! One carrier as a mass fraction of the mixture, WITH BOTH GHOST
      ! RULES OF THE CARRIER COLUMN APPLIED HERE AND NOWHERE ELSE, so that
      ! every evaluation of the advective term -- the stationary rows the
      ! certification measures, the trials of the fixed-wind relaxation and
      ! the derivative probes -- reads one ghost and returns one row for one
      ! interior state.
      !
      ! THE INNER GHOSTS carry the inflow composition: the handoff value
      ! where a handoff states one, the base cell's own where none does.  It
      ! is the same rule and the same expression carrier_mass_fractions
      ! applies to the declared set at the start of a marching step, written
      ! here on the carrier index this module solves in, so that the
      ! stationary rows and the stages put the same quantity on the faces.
      !
      ! THE OUTER GHOSTS carry the interior continued, X_ghost = X_N.  The
      ! outer face of this domain is an outflow face of a transported
      ! scalar: its characteristics leave the domain, so the boundary states
      ! no composition of its own and the ghost is the interior's own value
      ! continued.  The equilibrium composition the ionization sweep solves
      ! in a ghost cell (ionization_equilibrium loops 1-Ng to N+Ng) is not a
      ! transported quantity and may not stand here: it is the composition a
      ! cell of that density and temperature would relax to, not the
      ! composition the wind carried out of cell N.
      !
      ! WHERE THE FACE MASS FLUX REVERSES, F_rho(N) < 0, the rule does not
      ! change, because this boundary has no reservoir to state one with:
      ! the hydrodynamic closure at the same face (free_outflow_ghost in
      ! Apply_BC) continues cell N's own state outward -- the isothermal
      ! hydrostatic density at cell N's temperature, with the velocity
      ! copied -- and holds no composition of its own.  The material an
      ! inflowing outer face returns is therefore the material the domain
      ! released, which is cell N's composition, and the carrier rows are
      ! closed with the same ghost the hydrodynamic rows are closed with.
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
      do j = N+1, N+Ng
         Y(j) = Y(N)
      enddo
      where (Y .lt. 0.0d0) Y = 0.0d0
      end subroutine carrier_face_mass_fraction

      ! ------------------------------------------------------------- !

      ! THE SAME OUTFLOW CONTINUATION ON THE CARRIER STATE ITSELF, for the
      ! places that need a ghost VALUE and not a face fraction: the element
      ! budget, whose limiter reads the ghost cells because a ghost carries
      ! its own density, and the derivative probes, which perturb a column
      ! and must leave it closed the way the residual finds it.  One rule,
      ! stated in carrier_face_mass_fraction above.
      subroutine carrier_outflow_ghost(fc)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      integer :: ic
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         fc(N+1:N+Ng,ic) = fc(N,ic)
      enddo
      end subroutine carrier_outflow_ghost

      ! ------------------------------------------------------------- !

      ! Whether something below the base states this carrier's partition.
      ! Only the molecular handoff does -- base.inp's q_H2_base, or the
      ! profile's value at the matching level -- and it states the H2
      ! partition alone; every other carrier keeps the zero-gradient
      ! condition.  It is the same predicate
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
            ! this row, inside its stage flux (ionization_stage_transport,
            ! equation 1): it is a bulk mixing coefficient of the gas, blind
            ! to the charge of what it mixes.  Adding the ambipolar
            ! coefficient is a physics decision of its own and is not taken
            ! here.
            if (carrier_is_nucleus_fraction(ic)) then
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
      case (ic_HeII)
         isp = isp_HeII
      case (ic_HeIII)
         isp = isp_HeIII
      case (ic_HeTR)
         isp = isp_HeTR
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

      ! WHETHER CARRIER ic IS AN IONIZATION STAGE OF AN ELEMENT, and
      ! therefore a fraction of that element's NUCLEI carried on that
      ! element's nucleus flux, rather than a fraction per unit mass
      ! carried on the bulk face mass flux.
      !
      ! The two are the same physical object in different variables, and
      ! the nucleus form is the one whose stage fluxes sum to the element's
      ! own nucleus flux: the element's eddy flux is already inside that
      ! flux, so a stage written per unit mass on the bulk flux and given
      ! the whole mixing-ratio eddy flux -n_tot K d(n_k/n_tot)/dr carries
      ! the element's share of the eddy transport a second time
      ! (ionization_stage_transport).
      !
      ! The molecular carriers are NOT ionization stages: H2, OH, H2O and
      ! CO carry a molecular diffusion coefficient and a settling drift of
      ! their own, which the stage flux does not have, and the nuclei they
      ! hold belong to more than one element.
      logical function carrier_is_ionization_stage(ic) result(stage)
      integer, intent(in) :: ic
      stage = (ic .eq. ic_Hp) .or. (ic .eq. ic_HeII)                     &
                               .or. (ic .eq. ic_HeIII)
      end function carrier_is_ionization_stage

      ! WHETHER CARRIER ic IS A FRACTION OF ITS ELEMENT'S NUCLEI riding on
      ! that element's nucleus flux: the three ionization stages, and the
      ! He 2^3S level, which is a state of the helium partition with the
      ! charge of He I.  Every property the transport of a stage rests on
      ! -- the variable, the face flux, the reference density, the absence
      ! of a molecular diffusion and settling drift of its own -- holds for
      ! the level in the same form, so these are the carriers that take the
      ! stage flux of ionization_stage_transport.  What separates the level
      ! from the stages is its charge, which enters only the electron
      ! budget and the stage enthalpy flux, and those name the stages.
      logical function carrier_is_nucleus_fraction(ic) result(nucfrac)
      integer, intent(in) :: ic
      nucfrac = carrier_is_ionization_stage(ic) .or. (ic .eq. ic_HeTR)
      end function carrier_is_nucleus_fraction

      ! ------------------------------------------------------------- !

      ! THE ELEMENT WHOSE NUCLEI AN IONIZATION STAGE PARTITIONS, as the
      ! index of the element inventory: hydrogen for x(H II), helium for
      ! x(He II), x(He III) and the He 2^3S level.  One table, read by the nucleus density
      ! the row is written against, by the face flux it rides on and by
      ! the identity that sums the stages of one element.
      integer function carrier_stage_element(ic) result(ien)
      integer, intent(in) :: ic
      if (ic .eq. ic_Hp) then
         ien = ien_H
      else
         ien = ien_He
      endif
      end function carrier_stage_element

      ! ------------------------------------------------------------- !

      ! The nucleus density [cm^-3] of the element an ionization stage
      ! partitions, in cell j, from the frozen background of the step
      ! (carrier_state builds cbg_nHnuc and cbg_nHenuc from the one
      ! stoichiometric map).
      double precision function carrier_stage_nucleus_density(ic, j)     &
                               result(nel)
      integer, intent(in) :: ic, j
      if (carrier_stage_element(ic) .eq. ien_H) then
         nel = cbg_nHnuc(j)
      else
         nel = cbg_nHenuc(j)
      endif
      end function carrier_stage_nucleus_density

      ! ------------------------------------------------------------- !

      ! THE DENSITY CARRIER ic's UNKNOWN IS A FRACTION OF, in cell j
      ! [cm^-3]: the element's nucleus density for an ionization stage and
      ! the density of the mass row, nrho = rho n0, for every other
      ! carrier.  It multiplies the unknown to give the species density, it
      ! is what the time term of the row is written with, and it divides
      ! the density the write-back hands back, so the three agree by
      ! construction.
      double precision function carrier_nucleus_reference(ic, j, nrho_j)  &
                               result(nref)
      integer,  intent(in) :: ic, j
      real(dp), intent(in) :: nrho_j
      if (carrier_is_nucleus_fraction(ic)) then
         nref = carrier_stage_nucleus_density(ic, j)
      else
         nref = nrho_j
      endif
      end function carrier_nucleus_reference

      ! ------------------------------------------------------------- !

      ! The species densities [cm^-3] of every carrier in cell j, each from
      ! its own reference density above.  One place, so that the chemistry
      ! rows, the limiter, the row scales and the write-back cannot read a
      ! carrier fraction against different denominators.
      subroutine carrier_cell_densities(fcj, j, nrho_j, nc)
      real(dp), intent(in)  :: fcj(n_carrier_max)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: nrho_j
      real(dp), intent(out) :: nc(n_carrier_max)
      integer :: ic
      do ic = 1, n_carrier_max
         nc(ic) = fcj(ic)*carrier_nucleus_reference(ic, j, nrho_j)
      enddo
      end subroutine carrier_cell_densities

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
      integer  :: j, ic, jlast

      Agrd  = 0.0d0
      Bdrf  = 0.0d0
      updrf = 0
      call l22b_setup()
      ! NO DIFFUSIVE, EDDY OR SETTLING FLUX CROSSES EITHER END OF THE
      ! CARRIER COLUMN.  The loop below runs over the interior faces
      ! j = 1 .. N-1 alone, so Agrd and Bdrf stay zero at r_{1/2} and at
      ! r_{N+1/2} and the column is a closed box for that flux.  The rule is
      ! one rule for the operator, for the certification that measures its
      ! rows and for the fixed-wind relaxation that drives them to zero.
      !
      ! THE OUTER FACE.  Molecular diffusion and the settling drift
      ! G = (m_c - mbar) g/(kT) are Chapman-Enskog coefficients of a
      ! COLLISIONAL mixture, and the outermost cells of these domains are
      ! transitional.  MEASURED 2026-09-23 with collisional_validity.py on
      ! the current LHS 1140 b states: Kn_bulk 0.06-0.14 at the top cell
      ! (29 R_p), the largest species value 0.41-0.59, exobase above the
      ! domain.  That is not Kn > 1, so the zero face flux is a closure
      ! chosen for the reason below, not a property of collisionless gas.
      ! The composition still leaves through that face, by advection, which
      ! is a statement about the bulk motion and not about the collisions.
      !
      ! The alternative -- forming the outer face from the ghost the way an
      ! interior face is formed -- was measured and is not defensible.  With
      ! the stated outflow ghost X_ghost = X_N the gradient part of the flux
      ! vanishes identically and only the settling drift crosses the face;
      ! H2 is heavier than the mean particle of an ionized hydrogen wind, so
      ! G > 0 and the drift is INWARD, and the domain would gain carriers
      ! through its outer boundary out of nothing.  MEASURED on the frozen
      ! LHS 1140 b state `wellmixed/HeH0.083` at cell 500, r = 29.0 R_p:
      ! that inward flux is -2.9471e-02 cm^-3 s^-1 against an interior
      ! influx of +9.7857e-04, thirty times the largest term of the row,
      ! and the row measure of the cell goes from 2.4619e-02 to 5.0629e-01.
      !
      ! THE INNER FACE is closed for the reason the base handoff states:
      ! what the layer below hands over is a COMPOSITION on the inflowing
      ! base face, carried by the base mass flux, and
      ! carrier_face_mass_fraction puts it on that face.  A diffusive flux
      ! through the same face would add a second, independent statement of
      ! the same handoff.
      jlast = N - 1

      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         ! An ionization stage carries none of these.  Its whole face flux,
         ! advective and eddy together, is the stage flux of
         ! ionization_stage_transport written on its element's nucleus
         ! flux, and a gradient coefficient left nonzero here would put the
         ! eddy transport of that stage into the row a second time.
         if (carrier_is_nucleus_fraction(ic)) cycle
         ! Settling coefficient of this carrier, G = (m_i - m_bar) g/(k T).
         do j = 1-Ng, N+Ng
            Gco(j) = (carrier_mass_amu(ic)*mu - mbar(j))*gphys(j)   &
                     /(kb_erg*max(TK(j), 1.0d0))
         enddo
         do j = 1, jlast
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

      ! THE FACE STATE THE IONIZATION STAGE ROWS RIDE ON, frozen over one
      ! interval exactly as the carrier face coefficients above are: the
      ! element nucleus face flux N_el(f), the face nucleus density
      ! n_el(f), the face eddy coefficient K(f) and the face spacing dr(f)
      ! of equation (1) of ionization_stage_transport.
      !
      ! Every one of them comes from the element operator's one public flux
      ! (element_nucleus_face_flux, through
      ! hydrogen_and_helium_nucleus_face_flux).  Nothing is rebuilt here,
      ! which is what makes the stage fluxes of an element sum to that
      ! element's own nucleus flux face by face.
      !
      ! THE FACE MASS FLUX IS HANDED OVER IN g cm^-2 s^-1.  The rows' own
      ! F_rho is the code face mass flux of the mass row, with the mixture
      ! mass per unit of f_sp still outside it; the element operator wants
      ! the physical flux, so the mixture mass enters at the face, on the
      ! arithmetic mean every other face coefficient of that operator uses.
      ! Where the rows carry no advective term -- the marching operator,
      ! whose proton rode on the Runge-Kutta stages -- it is zero, and the
      ! stage flux is then the element's diffusive half plus the stage's
      ! own eddy term, which is the operator-split remainder the molecular
      ! rows carry there too.
      subroutine carrier_stage_face_state(rho, TK, f_sp, advect, msum,    &
                                          Frho)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho, TK
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      logical,                                  intent(in) :: advect
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: msum, Frho

      real(dp), dimension(1-Ng:N+Ng) :: Frho_cgs, N_H, N_He
      real(dp), dimension(0:N) :: n_Hf, n_Hef, Kf, drf
      integer :: j

      if (.not. allocated(stg_Nel))                                      &
         allocate(stg_Nel(1-Ng:N+Ng,n_carrier_max),                      &
                  stg_nelf(0:N,n_carrier_max),                           &
                  stg_Kf(0:N), stg_drf(0:N))
      stg_Nel  = 0.0d0
      stg_nelf = 0.0d0
      stg_Kf   = 0.0d0
      stg_drf  = 1.0d0
      if (.not. ionization_transport) return

      Frho_cgs = 0.0d0
      if (advect) then
         do j = 1-Ng, N+Ng-1
            Frho_cgs(j) = Frho(j)*0.5d0*(msum(j) + msum(j+1))*n0*mu*v0
         enddo
      endif
      call hydrogen_and_helium_nucleus_face_flux(rho, TK/T0, f_sp,       &
                    Frho_cgs, N_H, N_He, n_Hf, n_Hef, Kf, drf)
      stg_Nel(:,ic_Hp)  = N_H
      stg_nelf(:,ic_Hp) = n_Hf
      ! The helium stages ride on the helium nucleus flux of the SAME
      ! evaluation of the element operator, so the two elements' fluxes
      ! cannot come from different states of the gas.
      stg_Nel(:,ic_HeII)   = N_He
      stg_Nel(:,ic_HeIII)  = N_He
      stg_nelf(:,ic_HeII)  = n_Hef
      stg_nelf(:,ic_HeIII) = n_Hef
      ! The He 2^3S level moves with the helium nucleus flux, i.e. at the
      ! element's velocity like the rest of He I; a drift of the
      ! metastable against ground-state He I is not carried
      ! (md/He23S_transport_design_20261003.md section 4).
      stg_Nel(:,ic_HeTR)  = N_He
      stg_nelf(:,ic_HeTR) = n_Hef
      stg_Kf  = Kf
      stg_drf = drf

      end subroutine carrier_stage_face_state

      ! ------------------------------------------------------------- !

      ! THE FACE FLUX OF ONE IONIZATION STAGE ROW AND ITS TWO DERIVATIVES,
      ! in the shape the carrier residual stores every face flux in, so
      ! that one divergence and one block assembly serve both kinds of row.
      !
      !    F(f) = x(f) N_el(f) - n_el(f) K(f) [x(j+1) - x(j)]/dr(f)
      !
      ! with x(f) the reconstructed, limited and upwinded face fraction --
      ! the same rule the element mass fractions take, with the donor side
      ! chosen by the sign of the ELEMENT flux -- and the eddy term zero at
      ! the two end faces, where the element itself crosses only with the
      ! gas.  Equation (1) and its boundary paragraph in
      ! ionization_stage_transport.
      !
      ! The derivatives are taken at the donor-cell face value, exactly as
      ! carrier_advective_face_coefficients takes the molecular rows' at
      ! theirs: the limiter and the second-order part of the reconstruction
      ! are the strongly nonlinear part of the operator and reach beyond
      ! the tridiagonal band.  The residual is the full operator either
      ! way, so this changes what the Newton converges AT and not what it
      ! converges TO.
      ! ONE STAGE AT A TIME IS EXACT HERE: F_k depends on x_k alone, the
      ! stages of an element coupling only through the closing stage, which
      ! is not a row (ionization_stage_face_jacobian).  The closing face
      ! fraction this call forms is therefore the one of a single-stage
      ! element and is discarded; the identity that needs all the stages of
      ! an element together is ionization_stage_nucleus_sum.
      subroutine carrier_stage_face_flux(ic, fc, Jf, dJl, dJr)
      integer,                                      intent(in)  :: ic
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      real(dp), dimension(0:N),                     intent(out) :: Jf
      real(dp), dimension(0:N),                     intent(out) :: dJl, dJr

      real(dp), dimension(1,1-Ng:N+Ng) :: xion
      real(dp), dimension(1,0:N)       :: xf, Fk
      real(dp), dimension(0:N)         :: xclose, Fclose

      xion(1,:) = fc(:,ic)
      call ionization_stage_face_flux(xion, stg_Nel(:,ic),               &
                                      stg_nelf(:,ic), stg_Kf, stg_drf,   &
                                      xf, xclose, Fk, Fclose)
      Jf = Fk(1,:)
      call ionization_stage_face_jacobian(stg_Nel(:,ic), stg_nelf(:,ic), &
                                          stg_Kf, stg_drf, dJl, dJr)

      end subroutine carrier_stage_face_flux

      ! ------------------------------------------------------------- !

      ! THE IDENTITY THAT SAYS MOVING CHARGE BETWEEN STAGES CANNOT MOVE A
      ! NUCLEUS, measured on the state in hand:
      !
      !    sum_k F_k(f) = N_el(f)   at every face,
      !
      ! the sum running over ALL stages of the element, the closing neutral
      ! stage included (ionization_stage_transport, equation 2).  It holds
      ! by construction -- one element flux multiplies every stage, the
      ! closing face fraction is one minus the others and the closing eddy
      ! term is minus the sum of the others' -- so what is left in it is
      ! the rounding of the sums, and a measure above that says the
      ! construction has been broken somewhere.
      !
      ! The measure is relative to the largest of the element flux and the
      ! sum of the magnitudes of the stage fluxes, so a face across which
      ! no nucleus moves does not read as a failure of an identity it
      ! satisfies trivially.
      ! ONE ELEMENT PER CALL, with ALL of that element's carried stages in
      ! the same array: the identity is a sum over the stages of one
      ! element, so helium's two stages have to enter it together, sharing
      ! the one closing stage that the neutral helium is.
      !
      ! g_out is S/S' AT THE FACE THAT CARRIES dmax, with
      ! S  = |N| + |F_close| + sum_k |F_k| + sum_k |E_k| and S' the same sum
      ! without the eddy terms (ionization_stage_transport, equations 4 and
      ! 5).  It is what the rounding bound of the identity is proportional
      ! to, so the certification can gate this measure at the bound of the
      ! face it was taken at rather than at the ceiling of every face.  The
      ! eddy terms come back from the face-flux routine itself, so there is
      ! no second spelling of them here.  Negative on return means no face
      ! carried a nonzero scale.
      subroutine ionization_stage_nucleus_sum(fc, ien, dmax, jworst, g_out)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      integer,                                      intent(in)  :: ien
      real(dp),                                     intent(out) :: dmax
      integer,                                      intent(out) :: jworst
      real(dp), optional,                           intent(out) :: g_out

      real(dp), dimension(:,:), allocatable :: xion, xf, Fk, Ek
      real(dp), dimension(0:N)              :: xclose, Fclose
      real(dp) :: d, sc, Fsum, sprime, eddy, gworst
      integer  :: j, k, nk, ics(3), icf

      dmax   = 0.0d0
      jworst = 0
      gworst = -1.0d0
      if (present(g_out)) g_out = gworst
      if (.not. allocated(stg_Nel)) return
      if (ien .eq. ien_H) then
         if (.not. carrier_solved(ic_Hp)) return
         nk     = 1
         ics(1) = ic_Hp
      else
         if (.not. carrier_solved(ic_HeII))  return
         if (.not. carrier_solved(ic_HeIII)) return
         nk     = 2
         ics(1) = ic_HeII
         ics(2) = ic_HeIII
         ! The carried He 2^3S level is a third carried state of the
         ! helium partition; its flux is part of the helium nucleus flux,
         ! and the closing state is then the ground singlet.
         if (carrier_solved(ic_HeTR)) then
            nk     = 3
            ics(3) = ic_HeTR
         endif
      endif
      icf = ics(1)
      allocate(xion(nk,1-Ng:N+Ng), xf(nk,0:N), Fk(nk,0:N), Ek(nk,0:N))
      do k = 1, nk
         xion(k,:) = fc(:,ics(k))
      enddo
      call ionization_stage_face_flux(xion, stg_Nel(:,icf),              &
                                      stg_nelf(:,icf), stg_Kf,           &
                                      stg_drf, xf, xclose, Fk, Fclose,   &
                                      Ek_out = Ek)
      do j = 0, N
         Fsum = Fclose(j)
         sc   = abs(Fclose(j))
         eddy = 0.0d0
         do k = 1, nk
            Fsum = Fsum + Fk(k,j)
            sc   = sc + abs(Fk(k,j))
            eddy = eddy + abs(Ek(k,j))
         enddo
         sprime = abs(stg_Nel(j,icf)) + sc
         sc = max(abs(stg_Nel(j,icf)), sc)
         if (sc .le. 0.0d0) cycle
         d = abs(Fsum - stg_Nel(j,icf))/sc
         ! The ratio is kept for the face that carries the maximum, and for
         ! the first measured face, so that a column whose identity is
         ! exactly zero everywhere still reports the ratio of a real face
         ! instead of none.  dmax and jworst are what they were.
         if (d .gt. dmax .or. gworst .lt. 0.0d0)                          &
            gworst = (sprime + eddy)/max(sprime, 1.0d-300)
         if (d .gt. dmax) then
            dmax   = d
            jworst = j
         endif
      enddo
      if (present(g_out)) g_out = gworst
      deallocate(xion, xf, Fk, Ek)

      end subroutine ionization_stage_nucleus_sum

      ! ------------------------------------------------------------- !

      ! The identity above as the last stationary evaluation measured it,
      ! and whether an evaluation has measured it at all.  The
      ! certification reads it here rather than re-forming the stage
      ! fluxes, so the number it reports is the one the rows were built
      ! with.
      ! g is the ratio of magnitude sums at the face the measure was taken
      ! at, which the certification needs to gate the measure at the
      ! rounding bound of THAT face; negative means no face of the column
      ! carried a nonzero scale.
      subroutine ionization_stage_sum_measure(ien, dmax, jworst, known, g)
      integer,  intent(in)  :: ien
      real(dp), intent(out) :: dmax
      integer,  intent(out) :: jworst
      logical,  intent(out) :: known
      real(dp), optional, intent(out) :: g
      dmax   = stage_sum_max(ien)
      jworst = stage_sum_jworst(ien)
      known  = stage_sum_known(ien)
      if (present(g)) g = stage_sum_g(ien)
      end subroutine ionization_stage_sum_measure

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

      ! Test only: install the frozen helium background the projection's
      ! bound and the chemistry rows' neutral closure are formed from --
      ! the helium nuclei of a cell, the helium of its molecules and the
      ! population of its metastable level, one value for the whole column
      ! -- so that a suite can project constructed fractions without
      ! building a column and a sweep.  Production installs the same three
      ! from the species table in carrier_state.
      subroutine carrier_helium_background_set_for_test(nHe_nuc, nHe_mol, &
                                                        nHe_TR)
      real(dp), intent(in) :: nHe_nuc, nHe_mol, nHe_TR
      if (.not. allocated(cbg_nHenuc)) allocate(cbg_nHenuc(1-Ng:N+Ng))
      if (.not. allocated(cbg_nHemol)) allocate(cbg_nHemol(1-Ng:N+Ng))
      if (.not. allocated(cbg_nheiTR)) allocate(cbg_nheiTR(1-Ng:N+Ng))
      cbg_nHenuc = nHe_nuc
      cbg_nHemol = nHe_mol
      cbg_nheiTR = nHe_TR
      end subroutine carrier_helium_background_set_for_test

      ! ------------------------------------------------------------- !

      ! THE HELIUM NUCLEI THE TWO IONIZED STAGES OF CELL j MAY HOLD
      ! [cm^-3].
      !
      ! It is NOT the helium element total.  Two parts of that total are
      ! not the stages' to take, and both are frozen for the duration of
      ! the step:
      !
      !   the helium bound in HeH+ and in any further helium-bearing
      !   molecule (cbg_nHemol).  That nucleus is moved by the molecule's
      !   own row, and the write-back leaves the molecule where it is
      !   because rescaling it would move a hydrogen nucleus with it.  A
      !   simplex written without it admits x(He II) + x(He III) close to
      !   one beside a molecule holding the rest, and the write-back then
      !   hands back more helium than the cell has: at HeH+ holding 0.10
      !   of the inventory, x(He II) = 0.60 and x(He III) = 0.35 pass a
      !   bound of one and leave an inventory of 1.05.
      !
      !   the population of the He 2^3S metastable (cbg_nheiTR), where
      !   the step holds it frozen.  The level is a sub-population of He I
      !   and not a nucleus of its own: the chemistry rows read the ground
      !   singlet as n(He I) - n(He 2^3S), so a partition leaving less
      !   neutral helium than the frozen level holds has a negative
      !   singlet, which is not a state the rows may be evaluated at.
      !   Where "He 2^3S transport" carries the level it is not frozen and
      !   not reserved: it is the third carried state of the same simplex
      !   (the paragraph above carrier_helium_available_to_stages).
      !
      ! Reserving both makes the neutral singlet of a trial the difference
      ! of this density and the two ionized stages, which is nonnegative on
      ! the admissible set and needs no clip.  It is the helium
      ! counterpart of hydrogen_available_to_carriers, and the simplex, the
      ! chemistry rows and the write-back all go through it.
      ! HOW FAR OUTSIDE THE HELIUM INVENTORY THE LAST WRITE-BACK WAS ASKED
      ! TO GO, as a fraction of the helium nuclei of the cell that carried
      ! it: zero says every state the step wrote held its two ionized
      ! stages, its molecules and its frozen metastable inside the helium
      ! the cells have.  It is a measurement of the state, never a repair
      ! of it.
      double precision function carrier_helium_singlet_breach()          &
                                result(breach)
      breach = carrier_helium_singlet_under_zero
      end function carrier_helium_singlet_breach

      ! ------------------------------------------------------------- !

      ! WITH THE He 2^3S LEVEL CARRIED ("He 2^3S transport") it is not
      ! reserved: it is the third carried state of the same simplex, and
      ! reserving its frozen population as well would count it twice.  The
      ! helium of HeH+ stays reserved in either case.
      double precision function carrier_helium_available_to_stages(j)    &
                                result(nHe)
      integer, intent(in) :: j
      nHe = cbg_nHenuc(j) - cbg_nHemol(j)
      if (.not. carrier_solved(ic_HeTR)) nHe = nHe - cbg_nheiTR(j)
      if (nHe .lt. 0.0d0) nHe = 0.0d0
      end function carrier_helium_available_to_stages

      ! ------------------------------------------------------------- !

      ! THE NEUTRAL HELIUM OF A TRIAL STAGE PARTITION in cell j [cm^-3]:
      ! the ground singlet that closes the element, and the whole neutral
      ! stage, which is that singlet plus the frozen metastable.
      !
      ! NO CLIP.  On the admissible set of the stages the singlet is
      ! nonnegative by construction (carrier_helium_available_to_stages),
      ! so a negative value here is an inadmissible trial and not a state
      ! to be repaired: clipping it to zero while the two ionized stages
      ! stand where they are leaves the cell holding more helium than it
      ! has, and the rows are then evaluated at a partition no cell can be
      ! in.  How far the STATE THE STEP WROTE stood outside the set is
      ! MEASURED in carrier_write_back
      ! (carrier_helium_singlet_under_zero); this routine is called from
      ! inside the threaded Jacobian assembly and keeps no record of its
      ! own.
      !
      ! n_heiTR is the He 2^3S population of the trial: the carried
      ! unknown where "He 2^3S transport" carries it, the frozen one of
      ! the step otherwise.  In the first case the available helium holds
      ! the level (carrier_helium_available_to_stages) and the singlet is
      ! what the three carried states leave; in the second the level is
      ! reserved there and the singlet is what the two ionized stages
      ! leave.  Either way n_hei = n_heiSI + n_heiTR is the whole neutral
      ! stage.
      subroutine carrier_helium_neutral_of_partition(j, n_heii, n_heiii, &
                                                     n_heiTR, n_hei,     &
                                                     n_heiSI)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: n_heii, n_heiii, n_heiTR
      real(dp), intent(out) :: n_hei, n_heiSI
      n_heiSI = carrier_helium_available_to_stages(j) - n_heii - n_heiii
      if (carrier_solved(ic_HeTR)) n_heiSI = n_heiSI - n_heiTR
      n_hei   = n_heiSI + n_heiTR
      end subroutine carrier_helium_neutral_of_partition

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
      ! proton row comes back as exactly zero.
      double precision function carrier_element_reference_density(ic,     &
                               nH_free, nO_free, nC_free, nHe_nuc)       &
                               result(nref)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nH_free, nO_free, nC_free, nHe_nuc
      select case (ic)
      case (ic_H2, ic_Hp)
         nref = nH_free
      case (ic_HeII, ic_HeIII, ic_HeTR)
         ! THE HELIUM NUCLEUS COUNT OF THE CELL, which bounds a helium
         ! stage from above.  What a stage can actually reach is less by
         ! the helium held in HeH+ and in the frozen He 2^3S level
         ! (carrier_helium_available_to_stages), and that smaller density
         ! is what the admissible set is written with.  The count is the
         ! one used HERE because this table sets a SCALE -- the floor of
         ! the Jacobian's difference step and the absolute floor of the
         ! row -- and a scale that can fall to zero where an element is
         ! wholly bound in a molecule would leave the row's diagonal
         ! derivative at the 1e-300 guard.  The two
         ! agree to the molecular and metastable fractions of the helium,
         ! which are trace.
         nref = nHe_nuc
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
                               nH_free, nO_free, nC_free, nHe_nuc)       &
                               result(dn)
      integer,  intent(in) :: ic
      real(dp), intent(in) :: nc_ic, nH_free, nO_free, nC_free, nHe_nuc
      dn = max(1.0d-6*abs(nc_ic),                                        &
               1.0d-12*carrier_element_reference_density(ic, nH_free,    &
                                                         nO_free,        &
                                                         nC_free,        &
                                                         nHe_nuc),       &
               1.0d-300)
      end function carrier_source_derivative_step

      ! ------------------------------------------------------------- !

      ! Chemical production minus loss of the four carriers in cell j, at the
      ! trial densities nc(1..4) [cm^-3].  The rows are the SAME rows the
      ! local equilibrium solve uses (sec. 4); nothing is rewritten here.
      subroutine carrier_source(j, nc, nH_free, nO_free, src,             &
                                sprod, sloss, sphot, sh2chan)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: nc(n_carrier_max), nH_free, nO_free
      real(dp), intent(out) :: src(n_carrier_max)
      ! THE TERMS OF EACH ROW, NOT THEIR NET.  src is production minus loss;
      ! sprod and sloss are the two separately, both positive, and sphot is
      ! the radiative part of the loss.  Where the rows hand out no split --
      ! the helium stages of either gas, and the proton of an atomic gas --
      ! the entries are the POSITIVE AND NEGATIVE PARTS of the net source
      ! and not gross channels, and the comments at each say so.  In every
      ! gas and for every carrier sprod - sloss is the assembled src, after
      ! every correction made to a row once it had returned.  Where the chemistry is fast the two
      ! cancel to many digits and |src| stands orders below either, so a row
      ! that reads as one term is not a row with one term; the record of
      ! EXHALE_CARRIER_ROW_TERMS exists to tell those apart.  The physical
      ! scale of each row (row_terms_phys, carrier_residual), which the
      ! acceptance and the certification divide by, takes sprod + sloss as
      ! its reaction terms; the solution itself, and the Newton's own scale,
      ! read src alone.
      real(dp), optional, intent(out) :: sprod(n_carrier_max)
      real(dp), optional, intent(out) :: sloss(n_carrier_max)
      real(dp), optional, intent(out) :: sphot(n_carrier_max)
      ! THE H2 ROW CHANNEL BY CHANNEL, in the order of the n_h2chan_*
      ! parameters of this module: the individual reactions the H2
      ! production and loss are sums of, each a volumetric rate
      ! [cm^-3 s^-1].  Diagnostic only, like sprod/sloss/sphot.
      real(dp), optional, intent(out) :: sh2chan(n_h2chan)
      real(dp) :: fv(10)
      real(dp) :: pHp, lHp, pH2, lH2, lH2ph
      ! What the He <-> H charge-exchange pair adds to the proton row,
      ! taken as the change it makes to that row rather than rebuilt from
      ! the two rate expressions, which live in one place
      ! (charge_exchange).
      real(dp) :: cx_dHp
      real(dp) :: pOH, lOH, pH2O, lH2O, pH2ox, lH2ox
      real(dp) :: n_hi, n_hii, n_h2p, n_h3p, n_hehp
      real(dp) :: n_hei, n_heii, n_heiii, n_heiTR, n_heiSI, n_e, n_o0
      real(dp) :: h2chan(n_h2chan)
      ! The two halves of the He 2^3S row (tr_triplet_row_channels).
      real(dp) :: hetr_prod, hetr_loss
      ! The metal charge exchange as stoichiometric stage sources
      ! [cm^-3 s^-1]: the net source of each hydrogen and helium stage and
      ! the gross production and loss the same reactions give them.  The
      ! stage densities are what this operator's unknowns are, so these
      ! enter the rows with no conversion; the local systems take the same
      ! reactions in the boundary-flow basis their rows are written in
      ! (charge_exchange, module header).
      real(dp) :: cx_sH(0:1), cx_sHe(0:2)
      real(dp) :: cx_pH(0:1), cx_lH(0:1), cx_pHe(0:2), cx_lHe(0:2)

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
      ! THE He 2^3S POPULATION IS FROZEN BACKGROUND OR THE TRIAL UNKNOWN,
      ! for the reason the proton is: carried ("He 2^3S transport"), its
      ! row is solved here, and the triplet photoionization, collisional
      ! ionization and Penning terms of the He II, H II and singlet rows
      ! below then see the population the row is moving.
      n_heiTR = cbg_nheiTR(j)
      if (carrier_solved(ic_HeTR)) n_heiTR = nc(ic_HeTR)
      ! THE TWO IONIZED HELIUM STAGES ARE EITHER FROZEN BACKGROUND OR THE
      ! TRIAL UNKNOWNS, for the reason the proton is: with the helium
      ! ionization state carried, n(He+) and n(He++) are what this Newton
      ! is solving for in this cell, and rows evaluated at the last sweep's
      ! values would give the Jacobian columns of those rows a zero
      ! derivative.  The neutral helium then closes the element budget, the
      ! same closure the hydrogen has below: the helium AVAILABLE TO THE
      ! STAGES less the two ionized stages is the ground singlet, and the
      ! neutral stage is that singlet plus the frozen metastable.  The
      ! helium held in HeH+ and the population of the He 2^3S level are
      ! RESERVED from that budget rather than subtracted and clipped
      ! afterwards (carrier_helium_available_to_stages): on the admissible
      ! set of the two stages both differences are nonnegative, and the
      ! projection that defines that set is the one the limiter and every
      ! Newton trial go through.
      if (carrier_solved(ic_HeII)) then
         n_heii  = nc(ic_HeII)
         n_heiii = nc(ic_HeIII)
         call carrier_helium_neutral_of_partition(j, n_heii, n_heiii,    &
                                                  n_heiTR, n_hei,        &
                                                  n_heiSI)
      else
         ! The frozen background's own partition.  Its singlet is the
         ! floored difference of the two columns the sweep wrote, the
         ! floor and its reason stated at he_ground_singlet_density
         ! (a restart whose two columns come from different solves can put
         ! the difference on the wrong side of zero); it is written out
         ! here rather than called because this routine runs inside the
         ! threaded Jacobian assembly and that function keeps a count.
         n_heii  = cbg_nheii(j)
         n_heiii = cbg_nheiii(j)
         n_hei   = cbg_nhei(j)
         n_heiSI = max(n_hei - n_heiTR, 0.0d0)
      endif
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
      ! The helium electrons of the same sum, He+ once and He++ twice, moved
      ! to their trial values for the same reason.  n_e is a sum over
      ! charges, so replacing one term by its trial value is exact.
      if (carrier_solved(ic_HeII))                                           &
         n_e = max(n_e - cbg_nheii(j)  + n_heii                              &
                       - 2.0d0*cbg_nheiii(j) + 2.0d0*n_heiii, 0.0d0)

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

      fv     = 0.0d0
      h2chan = 0.0d0
      cx_dHp = 0.0d0
      pHp    = 0.0d0;  lHp   = 0.0d0
      pH2    = 0.0d0;  lH2   = 0.0d0;  lH2ph = 0.0d0
      pOH    = 0.0d0;  lOH   = 0.0d0
      pH2O   = 0.0d0;  lH2O  = 0.0d0
      pH2ox  = 0.0d0;  lH2ox = 0.0d0
      src    = 0.0d0

      ! WHICH BALANCE WRITES THE H AND He IONIZATION ROWS IS THE GAS, and
      ! in each case it is the balance the LOCAL SWEEP of that same gas
      ! solves.  A stage row of this operator and the sweep's row for the
      ! same stage are then one expression evaluated in two places, so a
      ! state stationary for one is stationary for the other; two
      ! transcriptions of "the same" chemistry could not be.
      if (thereis_mol) then
      ieq_cell%k_LW  = cph_klw(j)
      call set_mol_coeffs(ieq_cell%T_K, ieq_cell%ntot)
      call set_oxygen_coeffs(ieq_cell%T_K, cph_jh2o(j,:), cph_joh(j,:))
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
                        ieq_cell%q13, ieq_cell%q31g, ieq_cell%q31a,      &
                        ieq_cell%q31b, ieq_cell%Q31, ieq_cell%A31,       &
                        p_Hp = pHp, l_Hp = lHp, p_H2 = pH2,              &
                        l_H2 = lH2, l_H2_phot = lH2ph,               &
                        h2_chan = h2chan)
      ! He <-> H charge exchange, the same pair and the same orientation
      ! System_HeH_mol applies to the rows this operator has just taken
      ! from it: rows (1) and (2) are both written production positive, so
      ! he_row_sign = +1, and the CX reservoir is the ground singlet.
      !
      ! THE PAIR MOVES THE PROTON ROW AFTER mol_heh_rows HAS SPLIT IT, so
      ! the split is completed here.  He+ + H0 -> He0 + H+ makes a proton
      ! and He0 + H+ -> He+ + H0 removes one; their net is the change this
      ! call makes to row (1), and it is carried on the side of its own
      ! sign so that sprod - sloss is the assembled source of the row and
      ! not the source before the pair was added.  The two directions are
      ! not separated because the routine returns their difference, which
      ! is the same statement the helium rows below carry.
      ! Row (2) is the He II stage source here, which He2+ + H0 -> He+ +
      ! H+ feeds (heii_row_is_stage_source).
      cx_dHp = fv(1)
      call he_h_cx_fvec(fv, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,    &
                        ieq_cell%kcx_Hepp_H0,                            &
                        n_hi, n_hii, n_heiSI, n_heii, n_heiii, 1.0d0,     &
                        .true.)
      cx_dHp = fv(1) - cx_dHp
      ! The photolysis channels set_oxygen_coeffs built for this cell: the
      ! transport operator evaluates the chemistry at the cell's own,
      ! unscaled radiation field.
      call oxygen_carrier_rows(fv, 9, n_hi, nc(ic_H2), nc(ic_OH),        &
                               nc(ic_H2O), n_o0, oj3, oj4, oj5, oj7,     &
                               p_OH = pOH, l_OH = lOH, p_H2O = pH2O,     &
                               l_H2O = lH2O, p_H2_oxy = pH2ox,           &
                               l_H2_oxy = lH2ox)
      src(ic_H2)  = fv(4)
      src(ic_OH)  = fv(9)
      src(ic_H2O) = fv(10)
      else
      ! AN ATOMIC GAS: the H/He/He 2^3S balances of ion_residual_core, which
      ! is what System_HeH_TR and System_HeH solve in this configuration.
      ! The triplet coefficients are zero where the level is not tracked
      ! (ionization_equilibrium fills the cell state either way), and the
      ! four rows then reduce term by term to the H/He balances of
      ! System_HeH.  Row (4) is the He 2^3S balance (tr_triplet_row); it
      ! is read where "He 2^3S transport" carries the level, and left to
      ! the sweep's local root otherwise.
      call heh_tr_rows(fv, n_hi, n_hii, n_heiSI, n_heiTR, n_heii,        &
                       n_heiii, n_e,                                     &
                       ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,   &
                       ieq_cell%P_HeITR,                                 &
                       ieq_cell%rchiiB, ieq_cell%rcheiiB,                &
                       ieq_cell%rcheiiiB, ieq_cell%rcheiTR,              &
                       ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,            &
                       ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,        &
                       ieq_cell%q13, ieq_cell%q31g, ieq_cell%q31a,       &
                       ieq_cell%q31b, ieq_cell%Q31, ieq_cell%A31,        &
                       transport_operator = .true.)
      ! The same pair in the orientation those rows are written in: row (2)
      ! is the summed He I balance, He-I-gain positive, so he_row_sign =
      ! -1.  The reactant of He + H+ -> He+ + H is the GROUND SINGLET
      ! He(1^1S): the rate carries the barrier exp(-12.75/T4), and
      ! 12.75e4 K = 10.99 eV is the ionization-potential difference
      ! 24.587 - 13.598 eV of ground-state helium against hydrogen, while
      ! He(2^3S) lies 19.82 eV above the singlet and the same collision is
      ! exothermic for it.  System_HeH_TR passes the singlet for that
      ! reason, and this row is that row.
      call he_h_cx_fvec(fv, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,    &
                        ieq_cell%kcx_Hepp_H0,                            &
                        n_hi, n_hii, n_heiSI, n_heii, n_heiii, -1.0d0,    &
                        .false.)
      endif

      ! THE METAL CHARGE EXCHANGE OF THE SAME TWO BALANCES.  Group A of
      ! Huang et al. (2023) Table 4 moves an electron between a metal and
      ! hydrogen and is active whenever metals are present; group C moves
      ! one between a metal and helium under cx_full; the group E electron
      ! capture O2+ + H0 -> O+ + H+ is active under its own scale.  All of
      ! them are terms of the H+ and He+ balances the local sweep solves
      ! (cx_add_to_fvec, after the same rows), so a transported stage that
      ! did not carry them would answer a different balance of the same
      ! cell from the sweep that consumes it.  The metal stages are the
      ! frozen background of the step, and the He reactant of group C is
      ! the ground singlet, for the reason the pair above names.
      !
      ! The rate coefficients belong to a cell, so this cell's are loaded
      ! first; the evaluator refuses an assembly whose stored cell is not
      ! the caller's, and this routine runs on whatever thread the
      ! parallel loop gave the cell.
      cx_sH  = 0.0d0
      cx_sHe = 0.0d0
      cx_pH  = 0.0d0
      cx_lH  = 0.0d0
      cx_pHe = 0.0d0
      cx_lHe = 0.0d0
      if (thereis_metals .and.                                          &
          (ionization_transport .or. carrier_solved(ic_HeII))) then
         ! at the background electron density the sweep formed the
         ! cell's coefficients with, so that the rates are the sweep's
         call cx_set_cell(ieq_cell%T_K, cbg_ne(j))
         call charge_exchange_stage_sources(cbg_nm0(j,:), cbg_nm1(j,:), &
                        cbg_nm2(j,:), n_hi, n_hii, n_heiSI, n_heii,     &
                        n_heiii, ieq_cell%T_K, cx_sH, cx_sHe,           &
                        p_H = cx_pH, l_H = cx_lH,                       &
                        p_He = cx_pHe, l_He = cx_lHe)
      endif
      ! ---- CO: the one-sided destruction model -----------------------
      ! Two channels and no formation term, which is what "one-sided"
      ! means:
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
      if (ionization_transport) src(ic_Hp) = fv(1) + cx_sH(1)
      ! Rows (2) and (3) of the same block are the complete He+ and He++
      ! balances -- photoionization and collisional ionization of He I and
      ! of He+, the He(2^3S) channels, radiative and dielectronic
      ! recombination, the molecular sinks R17/R23 and the He+ + CO
      ! channel -- taken whole, exactly as row (1) is.  There is no helium
      ! chemistry written here that the local solve does not already have;
      ! the difference is only where the row is evaluated.
      src(ic_HeII)  = 0.0d0
      src(ic_HeIII) = 0.0d0
      if (carrier_solved(ic_HeII)) then
         if (thereis_mol) then
            src(ic_HeII)  = fv(2) + cx_sHe(1)
            src(ic_HeIII) = fv(3) + cx_sHe(2)
         else
            ! heh_tr_rows writes the SUMMED He I balance in row (2), He I
            ! gain positive, and the He++ balance in row (3).  A helium
            ! nucleus that leaves He I enters He+, and one that leaves He+
            ! enters He++, so the He+ balance is minus the He I balance
            ! less the He++ balance.  The three stages of the element then
            ! sum to zero source, which is what makes the transport of the
            ! partition conserve helium nuclei.
            src(ic_HeII)  = -fv(2) - fv(3) + cx_sHe(1)
            src(ic_HeIII) =  fv(3)         + cx_sHe(2)
         endif
      endif
      ! THE He 2^3S ROW IS THE LEVEL BALANCE THE SWEEP SOLVES, taken
      ! whole: tr_triplet_row (ion_residual_core), row (4) of the atomic
      ! rows and row (8) of the molecular ones, the latter with the
      ! He(2^3S) + H2 quench.  Production positive, as the stage rows.
      ! A transition between 2^3S and the singlet moves no nucleus out of
      ! He I, so the summed He I balance of row (2) and hence the He II
      ! source above are the same whichever of the two levels the neutral
      ! sits in; what the level changes there is which of its own channels
      ! (triplet photoionization, collisional ionization) runs, at the
      ! trial n_heiTR.  The metal charge exchange has no 2^3S channel
      ! (its helium reactant is the ground singlet).
      src(ic_HeTR) = 0.0d0
      if (carrier_solved(ic_HeTR)) then
         if (thereis_mol) then
            src(ic_HeTR) = fv(8)
         else
            src(ic_HeTR) = fv(4)
         endif
      endif

      ! THE SAME ROWS AS SUMS OF MAGNITUDES.  Every term of a row is
      ! sign-definite in the trial densities, so the sum of the magnitudes
      ! of the individual terms is the production plus the loss; grouping
      ! the like-signed terms therefore loses nothing.  CO has no formation
      ! channel in the one-sided model above, so its two loss channels are
      ! its whole scale.
      !
      ! THE He 2^3S ROW CARRIES ITS TWO GROSS HALVES.  Unlike the
      ! stage rows its balance is one closed expression of sign-definite
      ! channels (tr_triplet_row), so the halves are formed from the
      ! same terms by tr_triplet_row_channels, with the He(2^3S) + H2
      ! quench of the molecular row on the loss side; sprod - sloss is
      ! src to the rounding of the two sums.  Where the level is in
      ! near-local balance the net is orders below either half, which is
      ! what this record has to show.
      if (carrier_solved(ic_HeTR) .and.                            &
          (present(sprod) .or. present(sloss))) then
         call tr_triplet_row_channels(hetr_prod, hetr_loss,            &
              n_hi, n_heiSI, n_heiTR, n_heii, n_e,                     &
              ieq_cell%P_HeITR, ieq_cell%rcheiTR, ieq_cell%q13,        &
              ieq_cell%q31g, ieq_cell%q31a, ieq_cell%q31b,             &
              ieq_cell%Q31, ieq_cell%A31, ieq_cell%a_ion_HeITR)
         if (thereis_mol)                                             &
            hetr_loss = hetr_loss + mk_ion_H2*n_heiTR*nc(ic_H2)
      endif
      if (present(sprod)) then
         sprod        = 0.0d0
         sprod(ic_H2) = pH2 + pH2ox
         sprod(ic_OH) = pOH
         sprod(ic_H2O)= pH2O
         sprod(ic_CO) = 0.0d0
         ! THE PROTON CARRIES ITS PRODUCTION AND LOSS WHERE THE ROWS HAND
         ! THEM OUT, and its NET where they do not.  mol_heh_rows splits
         ! row (1); heh_tr_rows returns the four atomic rows and not their
         ! halves, so in an atomic gas the record states the balance and
         ! not two magnitudes it would have to invent.  Putting it on the
         ! side of its own sign keeps sprod - sloss equal to src for every
         ! row, which is what a reader of the record needs.
         !
         ! THE CORRECTIONS MADE AFTER THE ROWS RETURNED ARE IN IT.  The
         ! He <-> H charge-exchange pair moves row (1) after mol_heh_rows
         ! has split it, so its net is added here on the side of its own
         ! sign; the record of a molecular cell is otherwise the row
         ! before that pair, and sprod - sloss is not the source the
         ! residual was assembled from.
         !
         ! THE METAL CHARGE EXCHANGE IS IN IT AS ITS OWN TWO MAGNITUDES.
         ! Those reactions hand out the gross production and the gross
         ! loss of the stage separately, so the record states them and
         ! does not fold them into a net.
         if (ionization_transport) then
            if (thereis_mol) then
               sprod(ic_Hp) = pHp + max(cx_dHp, 0.0d0) + cx_pH(1)
            else
               sprod(ic_Hp) = max(src(ic_Hp), 0.0d0)
            endif
         endif
         ! THE HELIUM STAGES CARRY THE NET, ON THE SIDE IT FALLS, AND THE
         ! TWO ENTRIES ARE PARTS OF THAT NET AND NOT GROSS CHANNELS.
         ! mol_heh_rows hands out the production/loss split for rows (1)
         ! and (4) alone; rows (2) and (3) are returned as one balance, so
         ! what this record can state for them is that balance and not its
         ! two halves.  Putting it on the side of its own sign keeps
         ! sprod - sloss equal to src for every row, which is what a reader
         ! of the record needs; the two separate magnitudes of a helium row
         ! are not available and the record does not invent them.  A reader
         ! that needs the gross channels of a helium stage has to take them
         ! from the rows themselves.
         if (carrier_solved(ic_HeII)) then
            sprod(ic_HeII)  = max(src(ic_HeII),  0.0d0)
            sprod(ic_HeIII) = max(src(ic_HeIII), 0.0d0)
         endif
         if (carrier_solved(ic_HeTR)) sprod(ic_HeTR) = hetr_prod
      endif
      if (present(sloss)) then
         sloss        = 0.0d0
         sloss(ic_H2) = lH2 + lH2ox
         sloss(ic_OH) = lOH
         sloss(ic_H2O)= lH2O
         sloss(ic_CO) = abs(src(ic_CO))
         if (ionization_transport) then
            if (thereis_mol) then
               sloss(ic_Hp) = lHp + max(-cx_dHp, 0.0d0) + cx_lH(1)
            else
               sloss(ic_Hp) = max(-src(ic_Hp), 0.0d0)
            endif
         endif
         if (carrier_solved(ic_HeII)) then
            sloss(ic_HeII)  = max(-src(ic_HeII),  0.0d0)
            sloss(ic_HeIII) = max(-src(ic_HeIII), 0.0d0)
         endif
         if (carrier_solved(ic_HeTR)) sloss(ic_HeTR) = hetr_loss
      endif
      if (present(sh2chan)) then
         sh2chan = 0.0d0
         sh2chan(1:n_h2chan_mol) = h2chan(1:n_h2chan_mol)
         sh2chan(n_h2chan_oxy_prod) = pH2ox
         sh2chan(n_h2chan_oxy_loss) = lH2ox
      endif
      if (present(sphot)) then
         sphot        = 0.0d0
         sphot(ic_H2) = lH2ph
         sphot(ic_OH) = oj7*nc(ic_OH)
         sphot(ic_H2O)= (oj3 + oj4 + oj5)*nc(ic_H2O)
         sphot(ic_CO) = cph_kco(j)*nc(ic_CO)
         if (ionization_transport) sphot(ic_Hp) = 0.0d0
      endif

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

      ! WHETHER THE ENERGY EQUATION CARRIES THE ENTHALPY OF THE RELATIVE
      ! TRANSPORT OF THIS OPERATOR'S CARRIERS: the molecular carriers
      ! moving relative to the rest of their elements, and the ionization
      ! stages that "Ionization transport" carries moving relative to their
      ! own element.  Both are removed, with the element term, by
      ! "Interdiffusion enthalpy flux: False" (the one key for the enthalpy
      ! of relative transport, parameters.f90).
      logical function carrier_enthalpy_active() result(on)
      on = molecular_carrier_enthalpy_active() .or.                      &
           ionization_stage_enthalpy_active()
      end function carrier_enthalpy_active

      ! The molecular part: where a molecular carrier (H2, OH, H2O, CO) is
      ! transported relative to the gas.
      logical function molecular_carrier_enthalpy_active() result(on)
      integer :: ic
      on = .false.
      if (.not. (interdiffusion_enthalpy_flux .and. thereis_mol .and.     &
                 carrier_transport)) return
      do ic = ic_H2, ic_CO
         if (carrier_solved(ic)) on = .true.
      enddo
      end function molecular_carrier_enthalpy_active

      ! The stage part: where the ionized stages are carried, each with its
      ! own eddy term relative to its element (ionization_stage_transport,
      ! equation 1).
      logical function ionization_stage_enthalpy_active() result(on)
      on = interdiffusion_enthalpy_flux .and. ionization_transport .and.  &
           carrier_solved(ic_Hp)
      end function ionization_stage_enthalpy_active

      ! ------------------------------------------------------------- !

      ! Which species of the table are the transported molecular carriers
      ! of this run.
      subroutine molecular_carrier_species(car_bsp)
      logical, intent(out) :: car_bsp(n_bsp)
      integer :: ic, ib, isp
      car_bsp = .false.
      do ic = ic_H2, ic_CO
         if (.not. carrier_solved(ic)) cycle
         isp = carrier_species_index(ic)
         do ib = 1, n_bsp
            if (bsp_fsp(ib) .eq. isp) car_bsp(ib) = .true.
         enddo
      enddo
      end subroutine molecular_carrier_species

      ! ------------------------------------------------------------- !

      subroutine recoiling_nucleus_enthalpies(f_sp, j, kTe, urv, car_bsp, &
                                              hbH, hbO, hbC)
      ! THE SENSIBLE ENTHALPY PER NUCLEUS OF THE SPECIES THAT TAKE UP A
      ! RECOIL in one cell [erg per nucleus], for hydrogen, oxygen and
      ! carbon.  When a carrier (a molecule, or a carried ionization stage)
      ! moves nuclei across a face that its element's flux does not move,
      ! the same nuclei move the other way in the species the composition
      ! update rescales to keep the element total (carrier_write_back), and
      ! in no others: the energy row must move the enthalpy of the nuclei
      ! the composition rows move.
      !
      ! Hydrogen recoils in the species made of hydrogen nuclei alone that
      ! the flow does not carry: H I, H2+, H3+, and H+ unless "Ionization
      ! transport" carries it (carrier_solved(ic_Hp)), in which case the
      ! write-back sets H+ from its own transported row and H+ takes no
      ! part of the recoil.  The molecular carriers (car_bsp) take none:
      ! their densities are their own rows'.  HeH+ never recoils: rescaling
      ! it would move a helium nucleus with the hydrogen one, so the
      ! write-back leaves it where it is.  O and C recoil in their atomic
      ! stages where the metals are in the equation of state.  Per particle
      ! h_s = e_s + kT with e_s of the caloric EOS: gamma_ad/(gamma_ad - 1)
      ! kT (1 + Z_s), each particle with the electrons its charge gave up,
      ! H2 adding its rovibrational energy k u_rv (urv, zero under "Caloric
      ! EOS: monatomic").  An element with no recoiling species in the cell
      ! recoils with its neutral ground atom, gamma_ad/(gamma_ad - 1) kT.
      ! j is the cell and kTe its kT [erg].
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: kTe, urv
      logical,  intent(in)  :: car_bsp(n_bsp)
      real(dp), intent(out) :: hbH, hbO, hbC
      real(dp) :: cp_part, hs, dens
      real(dp) :: eH, nHnuc, eO, nOnuc, eC, nCnuc
      integer  :: ib, im

      cp_part = gamma_ad/(gamma_ad - 1.0d0)
      eH = 0.0d0;  nHnuc = 0.0d0
      eO = 0.0d0;  nOnuc = 0.0d0
      eC = 0.0d0;  nCnuc = 0.0d0
      ! Each recoiling particle of nH nuclei moves nH nuclei and its own
      ! enthalpy hs.
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib) .or. car_bsp(ib)) cycle
         if (bsp_nH(ib) .le. 0) cycle
         if (bsp_nHe(ib) + bsp_nO(ib) + bsp_nC(ib) .gt. 0) cycle
         if (bsp_fsp(ib) .eq. isp_HII .and. carrier_solved(ic_Hp)) cycle
         dens = f_sp(j,bsp_fsp(ib))
         if (dens .le. 0.0d0) cycle
         hs = cp_part*kTe*(1.0d0 + dble(bsp_charge(ib)))
         if (bsp_fsp(ib) .eq. isp_H2) hs = hs + kb_erg*urv
         eH = eH + dens*hs
         nHnuc = nHnuc + dens*dble(bsp_nH(ib))
      enddo
      if (eos_include_metals .and. thereis_metals) then
         do im = 1, n_mion
            dens = f_sp(j,mion_fsp(im))
            if (dens .le. 0.0d0) cycle
            hs = cp_part*kTe*(1.0d0 + dble(mion_stage(im)))
            if (mion_elem(im) .eq. iel_O) then
               eO = eO + dens*hs;  nOnuc = nOnuc + dens
            else if (mion_elem(im) .eq. iel_C) then
               eC = eC + dens*hs;  nCnuc = nCnuc + dens
            endif
         enddo
      endif
      hbH = cp_part*kTe;  hbO = cp_part*kTe;  hbC = cp_part*kTe
      if (nHnuc .gt. 0.0d0) hbH = eH/nHnuc
      if (nOnuc .gt. 0.0d0) hbO = eO/nOnuc
      if (nCnuc .gt. 0.0d0) hbC = eC/nCnuc
      end subroutine recoiling_nucleus_enthalpies

      ! ------------------------------------------------------------- !

      subroutine molecular_carrier_enthalpy_face_flux(rho, Tcode, f_sp, q)
      ! THE SENSIBLE ENTHALPY FLUX OF THE MOLECULAR CARRIERS at the faces
      ! f = 0 ... N [erg cm^-2 s^-1], positive outward (code audit of
      ! 2026-09-29, md/CODE_AUDIT_20260929.md F1).
      !
      ! The energy equation of a mixture whose species move relative to the
      ! gas carries q_d = sum_s h_s J_s (Cook 2009, Phys. Fluids 21, 055109,
      ! eqs. 11-13; interdiffusion_enthalpy_face_flux states the physics and
      ! carries the ELEMENT part, in which every species of an element moves
      ! with its element).  A carrier moves relative to the rest of its own
      ! elements as well: its diffusive, eddy and settling flux Phi_c
      ! (carrier_face_flux, the flux the carrier rows are the divergence of)
      ! moves nuclei that the element fluxes do not, so the nuclei of the
      ! same elements held by the species that are NOT carriers move the
      ! other way, nucleus for nucleus, and the element fluxes are left as
      ! they are.  The carrier term is therefore
      !
      !    q_car = sum_c Phi_c [ h_c - sum_el nu_{c,el} hbar_el ] ,
      !
      ! h_c the enthalpy of one carrier particle and hbar_el the enthalpy
      ! per nucleus of el of the species that take up the recoil, each
      ! particle with the electrons its charge gave up
      ! (recoiling_nucleus_enthalpies states the recoil set: the species
      ! carrier_write_back rescales, and no others).  Per particle
      ! h_s = e_s + kT with e_s of the caloric EOS: gamma_ad/(gamma_ad - 1)
      ! kT (1 + Z_s), H2 adding its rovibrational energy k u_rv (zero under
      ! "Caloric EOS: monatomic").  So an H2 molecule rising through atomic
      ! hydrogen at zero hydrogen flux carries h_H2 - 2 h_H = k (u_rv -
      ! 5 T/2) upward per molecule.  Formation energy is not in these
      ! enthalpies, for the reason component_specific_enthalpies gives.
      ! The ionization stages that "Ionization transport" carries are the
      ! other carriers of this operator; their term is
      ! ionization_stage_enthalpy_face_flux, and the two share the recoil
      ! set, so no nucleus is counted by both.
      !
      ! THE FLUX IS THE CARRIER OPERATOR'S OWN, for the state in hand: the
      ! coefficients of carrier_face_coefficients from carrier_diffusivities
      ! and carrier_geometry, and carrier_face_flux on the unknown f_sp and
      ! the mixing-ratio factor n rho/n_tot.  The temperature and the
      ! particle count are those of the state passed (Tcode and the particle
      ! count of its composition, calc_ntot's), not the frozen background of
      ! a carrier step, so the term follows the iterate of a steady solve.
      ! Only the faces 1 ... N-1 carry a carrier flux (both ends of the
      ! carrier column are closed, carrier_face_coefficients), so q(0) =
      ! q(N) = 0 and the term moves energy within the column only.  Face
      ! values of h are the arithmetic means of their two cells.  Writes no
      ! module state (l22b_setup reads its options once).
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(0:N),                 intent(out) :: q

      real(dp), dimension(1-Ng:N+Ng) :: nd, ntot, TK, mbar, wfac
      real(dp), dimension(1-Ng:N+Ng) :: rp, rep, dtd, dtp, gphys
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco, dh
      real(dp), dimension(0:N,n_carrier_max) :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max) :: updrf
      logical  :: car_bsp(n_bsp)
      real(dp) :: cp_part, kTe, urv, crv
      real(dp) :: hbH, hbO, hbC, hc
      real(dp) :: Jfc, dJl, dJr
      integer  :: j, ic, ib, isp, im

      q = 0.0d0
      if (.not. molecular_carrier_enthalpy_active()) return

      call molecular_carrier_species(car_bsp)

      cp_part = gamma_ad/(gamma_ad - 1.0d0)
      nd = rho*n0
      do j = 1-Ng, N+Ng
         ! The gas-particle count of calc_ntot: one particle for every
         ! species of the table but the excited level, and the metal stages
         ! where the equation of state counts them.
         ntot(j) = 0.0d0
         do ib = 1, n_bsp
            if (bsp_is_excited_level(ib)) cycle
            ntot(j) = ntot(j) + f_sp(j,bsp_fsp(ib))*nd(j)
         enddo
         if (eos_include_metals .and. thereis_metals) then
            do im = 1, n_mion
               ntot(j) = ntot(j) + f_sp(j,mion_fsp(im))*nd(j)
            enddo
         endif
         TK(j) = max(Tcode(j)*T0, 1.0d0)
         if (ntot(j) .gt. 0.0d0) then
            mbar(j) = rho(j)*n0*mu/ntot(j)
         else
            mbar(j) = mu
         endif
         wfac(j) = nd(j)/max(ntot(j), 1.0d-99)
      enddo
      dtd = 1.0d0
      call carrier_geometry(dtd, rp, rep, dtp, gphys)
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)

      ! The enthalpy difference of every carrier against the nuclei that
      ! recoil, cell by cell [erg per carrier particle].
      dh = 0.0d0
      do j = 1-Ng, N+Ng
         kTe = kb_erg*TK(j)
         urv = 0.0d0
         call h2_rovibrational_energy_and_heat_capacity(TK(j), urv, crv)
         call recoiling_nucleus_enthalpies(f_sp, j, kTe, urv, car_bsp,   &
                                           hbH, hbO, hbC)
         do ib = 1, n_bsp
            if (.not. car_bsp(ib)) cycle
            hc = cp_part*kTe*(1.0d0 + dble(bsp_charge(ib)))
            if (bsp_fsp(ib) .eq. isp_H2) hc = hc + kb_erg*urv
            do ic = ic_H2, ic_CO
               if (.not. carrier_solved(ic)) cycle
               if (carrier_species_index(ic) .ne. bsp_fsp(ib)) cycle
               dh(j,ic) = hc - (dble(bsp_nH(ib))*hbH + dble(bsp_nO(ib))*hbO &
                                + dble(bsp_nC(ib))*hbC)
            enddo
         enddo
      enddo

      do ic = ic_H2, ic_CO
         if (.not. carrier_solved(ic)) cycle
         isp = carrier_species_index(ic)
         do j = 1, N-1
            call carrier_face_flux(f_sp(j,isp), f_sp(j+1,isp), wfac(j),   &
                                   wfac(j+1), Agrd(j,ic), Bdrf(j,ic),     &
                                   updrf(j,ic), Jfc, dJl, dJr)
            q(j) = q(j) + Jfc*0.5d0*(dh(j,ic) + dh(j+1,ic))
         enddo
      enddo
      end subroutine molecular_carrier_enthalpy_face_flux

      ! ------------------------------------------------------------- !

      subroutine ionization_stage_enthalpy_face_flux(rho, Tcode, f_sp, q)
      ! THE SENSIBLE ENTHALPY THE IONIZATION STAGES CARRY BY MOVING
      ! RELATIVE TO THEIR OWN ELEMENT, at the faces f = 0 ... N
      ! [erg cm^-2 s^-1], positive outward.
      !
      ! With "Ionization transport" stage k of an element crosses a face
      ! with the flux (ionization_stage_transport, equation 1)
      !
      !    F_k = x_k(f) N_el(f) + E_k ,
      !    E_k = - n_el(f) K(f) [x_k(j+1) - x_k(j)]/dr(f) ,
      !
      ! the first term moving the stage with its element and E_k, the
      ! stage's own eddy term, moving it relative to the element; the
      ! closing (neutral) stage carries - sum_k E_k, so the E_k move no
      ! nucleus.  The enthalpy of the element's own motion is carried by the
      ! hydrodynamic flux (E + p) u and, relative to the mass-weighted
      ! velocity, by the element term (interdiffusion_enthalpy_face_flux,
      ! whose h_He averages over the stages of the cell); what neither
      ! carries is the enthalpy of the E_k.  By q_d = sum_s h_s J_s (Cook
      ! 2009, Phys. Fluids 21, 055109, eqs. 11-13) it is
      !
      !    q_st = sum_el sum_k E_k [ h_k - hbar_el ] ,
      !
      ! h_k the sensible enthalpy of one particle of stage k with the Z_k
      ! electrons its charge gave up and hbar_el the enthalpy per nucleus of
      ! the species that take up the recoil, those whose densities the
      ! composition update rescales (carrier_write_back;
      ! recoiling_nucleus_enthalpies): hydrogen recoils in H I, H2+ and H3+
      ! (H+ is the carried stage itself and H2 a carrier with its own flux
      ! and its own term, molecular_carrier_enthalpy_face_flux, so no
      ! hydrogen nucleus is counted by both), helium in He I with the He
      ! 2^3S level inside it (HeH+ is reserved from the stages and never
      ! rescaled).  Where "He 2^3S transport" carries the level it is a
      ! carried state with charge 0 and the sensible enthalpy of a neutral
      ! atom, so h_3 - hbar_He = 0 and its eddy term adds nothing here; the
      ! E_HeII and E_HeIII below do not depend on whether the level is
      ! carried (equation 1 of ionization_stage_transport: E_k is a
      ! function of x_k alone).  Its excitation energy is chemical energy
      ! (helium_metastable_excitation_energy_divergence).
      !
      ! WHAT A STAGE EXCHANGE CARRIES.  Per particle h = e + kT with e of
      ! the caloric EOS (component_specific_enthalpies): gamma_ad/(gamma_ad
      ! - 1) kT for every atom, ion and electron, so the heavy particle of
      ! every stage of an element carries the same enthalpy, and the stage
      ! of charge Z_k adds the Z_k electrons the zero-current closure makes
      ! follow it.  A stage exchange at fixed element flux therefore moves
      ! the electrons' enthalpy alone: with the recoil in the neutral atom,
      !
      !    q_st = h_e sum_el sum_k Z_k E_k ,  h_e = gamma_ad/(gamma_ad - 1) kT,
      !
      ! (5/2) kT for each electron of the net electron flux sum_k Z_k E_k
      ! relative to the nuclei.  For helium that is exact here (He II and
      ! He III against He I: h_e E_HeII + 2 h_e E_HeIII); for hydrogen the
      ! recoil also takes H2+ (one electron per two nuclei) and H3+ (one per
      ! three), whose share lowers hbar_H, and that is carried as it stands.
      !
      ! FORMATION (IONIZATION) ENERGY IS NOT IN q_st, for the reason
      ! component_specific_enthalpies gives: the conserved energy of this
      ! code is kinetic plus sensible, and the stationary energy row books
      ! chemical energy only as the heating and cooling of the reactions
      ! where the ions are made and where they recombine.
      !
      ! VALIDITY.  (i) Zero current: the electrons move with the ions,
      ! n_e w_e = sum_s Z_s n_s w_s, the ambipolar closure of every enthalpy
      ! flux of this code.  (ii) No stage drift beyond the eddy term: a
      ! stage moves relative to its element by E_k alone, the one-velocity
      ! approximation of ionization_stage_transport, whose header measures
      ! the ion-neutral drift it leaves out at 0.02 to 0.6 of the advective
      ! stage flux outward of about 1.2 R_p on LHS 1140 b; the enthalpy of
      ! that drift is absent with the drift.  (iii) One temperature for
      ! electrons, ions and neutrals: h_e is that of the gas T.  (iv) As in
      ! Cook (2009), the Dufour flux and the kinetic energy of the relative
      ! motion are left out.
      !
      ! THE FLUX IS THE STAGE ROWS' OWN.  E_k is returned by
      ! ionization_stage_face_flux (Ek_out) on the face nucleus density, eddy
      ! coefficient and spacing of hydrogen_and_helium_nucleus_face_flux and
      ! the fractions x_k = n_k/n_el of the state passed, which is how the
      ! stage rows form it (carrier_state, carrier_stage_face_state); E_k
      ! does not depend on the element flux N_el, so no face mass flux is
      ! handed over.  E_k is zero at faces 0 and N (the element crosses them
      ! only with the gas), so q(0) = q(N) = 0 and the term moves energy
      ! within the column only.  Face values of h_k - hbar_el are the
      ! arithmetic means of their two cells, as in the molecular term.
      ! Writes no module state.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(0:N),                 intent(out) :: q

      real(dp), dimension(1-Ng:N+Ng) :: Fnone, N_H, N_He, nHc, nHec, nd
      real(dp), dimension(1-Ng:N+Ng) :: hgas, dhHp
      real(dp), dimension(0:N)       :: n_Hf, n_Hef, Kf, drf
      real(dp), dimension(0:N)       :: xclose, Fclose
      real(dp), dimension(1,1-Ng:N+Ng) :: xH
      real(dp), dimension(2,1-Ng:N+Ng) :: xHe
      real(dp), dimension(1,0:N)     :: xfH, FkH, EH
      real(dp), dimension(2,0:N)     :: xfHe, FkHe, EHe
      logical  :: car_bsp(n_bsp)
      real(dp) :: cp_part, TKc, kTe, urv, crv, hbH, hbO, hbC
      integer  :: j

      q = 0.0d0
      if (.not. ionization_stage_enthalpy_active()) return

      Fnone = 0.0d0
      call hydrogen_and_helium_nucleus_face_flux(rho, Tcode, f_sp, Fnone, &
                    N_H, N_He, n_Hf, n_Hef, Kf, drf,                      &
                    n_H_cell = nHc, n_He_cell = nHec)
      nd = rho*n0
      xH(1,:) = f_sp(:,isp_HII)*nd/max(nHc, 1.0d-300)
      call ionization_stage_face_flux(xH, N_H, n_Hf, Kf, drf, xfH,        &
                                      xclose, FkH, Fclose, Ek_out = EH)
      EHe = 0.0d0
      if (carrier_solved(ic_HeII) .and. carrier_solved(ic_HeIII)) then
         xHe(1,:) = f_sp(:,isp_HeII) *nd/max(nHec, 1.0d-300)
         xHe(2,:) = f_sp(:,isp_HeIII)*nd/max(nHec, 1.0d-300)
         call ionization_stage_face_flux(xHe, N_He, n_Hef, Kf, drf, xfHe, &
                                         xclose, FkHe, Fclose,            &
                                         Ek_out = EHe)
      endif

      ! h_k - hbar_el of every carried stage, cell by cell [erg]: H II
      ! against the hydrogen recoil, He II and He III against He I, whose
      ! enthalpy per nucleus is that of a neutral atom.
      cp_part = gamma_ad/(gamma_ad - 1.0d0)
      call molecular_carrier_species(car_bsp)
      do j = 1-Ng, N+Ng
         TKc = max(Tcode(j)*T0, 1.0d0)
         kTe = kb_erg*TKc
         urv = 0.0d0
         if (thereis_mol)                                                &
            call h2_rovibrational_energy_and_heat_capacity(TKc, urv, crv)
         call recoiling_nucleus_enthalpies(f_sp, j, kTe, urv, car_bsp,   &
                                           hbH, hbO, hbC)
         hgas(j) = cp_part*kTe
         dhHp(j) = 2.0d0*hgas(j) - hbH
      enddo
      do j = 0, N
         q(j) = EH(1,j)*0.5d0*(dhHp(j) + dhHp(j+1))                      &
              + (EHe(1,j) + 2.0d0*EHe(2,j))*0.5d0*(hgas(j) + hgas(j+1))
      enddo
      end subroutine ionization_stage_enthalpy_face_flux

      ! ------------------------------------------------------------- !

      subroutine carrier_enthalpy_divergence_of_state(rho, Tcode, f_sp,  &
                                                      divq, q_out)
      ! The divergence of the enthalpy flux of this operator's carriers,
      ! molecular_carrier_enthalpy_face_flux plus
      ! ionization_stage_enthalpy_face_flux, in every cell 1 ... N of the
      ! carrier column, in the code units of the energy row
      ! (q0 = n0 mu v0^3/R0), on the carrier rows' own geometry
      ! (spherical_face_area_and_cell_volume):
      !
      !    divq(j) = [ A_j q(j) - A_{j-1} q(j-1) ] / (V_j R0 q0) ,  j = 1 ... N.
      !
      ! Cell 1 is included: the carrier column, unlike the element column
      ! (whose cell 1 is the prescribed reservoir), carries a balance in
      ! every cell 1 ... N, stage rows included.  With q(0) = q(N) = 0,
      ! sum_j V_j divq(j) = 0 to rounding.  The ghosts are returned as zero.
      ! Zero, with no arithmetic, unless carrier_enthalpy_active().  q_out
      ! (optional) returns the face flux [erg cm^-2 s^-1].
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: divq
      real(dp), dimension(0:N), optional,       intent(out) :: q_out
      real(dp), dimension(0:N) :: q, qst, fa
      real(dp), dimension(1:N) :: cv
      integer :: j

      divq = 0.0d0
      if (present(q_out)) q_out = 0.0d0
      if (.not. carrier_enthalpy_active()) return
      call molecular_carrier_enthalpy_face_flux(rho, Tcode, f_sp, q)
      if (ionization_stage_enthalpy_active()) then
         call ionization_stage_enthalpy_face_flux(rho, Tcode, f_sp, qst)
         q = q + qst
      endif
      call spherical_face_area_and_cell_volume(fa, cv)
      do j = 1, N
         divq(j) = (fa(j)*q(j) - fa(j-1)*q(j-1))/(cv(j)*R0*q0)
      enddo
      if (present(q_out)) q_out = q
      end subroutine carrier_enthalpy_divergence_of_state

      ! ------------------------------------------------------------- !

      subroutine helium_metastable_excitation_energy_divergence(rho,     &
                                        Tcode, f_sp, divq, q_out)
      ! THE EXCITATION ENERGY THE CARRIED He 2^3S POPULATION MOVES BY
      ! MOVING RELATIVE TO He I, at the faces and as a divergence in the
      ! code units of the energy row (q0 = n0 mu v0^3/R0), on the carrier
      ! rows' geometry (the formula of carrier_enthalpy_divergence_of_state).
      !
      ! With "He 2^3S transport" the level crosses a face with the flux of
      ! ionization_stage_transport equation (1), F_3 = x3(f) N_He(f) + E_3,
      ! and He I (the singlet with the level inside it) with
      ! F_I = x_I(f) N_He(f) + E_I, E_I = -(E_HeII + E_HeIII) the eddy
      ! term of the whole neutral stage.  The advective halves move the
      ! level at the velocity of He I, so the flux of the level RELATIVE to
      ! He I is the eddy part alone,
      !
      !    q_x(f) = eps_3 [ E_3(f) - (x3/x_I)(f) E_I(f) ] ,
      !
      ! eps_3 = 19.82 eV the energy of 2^3S above the ground singlet
      ! (species_formation_energy), face fractions the arithmetic means of
      ! their two cells.  q_x is zero at faces 0 and N (no eddy term there).
      !
      ! IT IS NOT A TERM OF THE ENERGY EQUATION OF THIS CODE, and this
      ! routine is a measurement.  The conserved energy is kinetic plus
      ! sensible; the energy of the 2^3S level is chemical energy, booked
      ! where the reactions that make and destroy the level run: the
      ! 1^1S -> 2^3S excitation cools with q13 and the de-excitation heats
      ! with q31g at the transported population (radiative_cooling_of_cell;
      ! lambda_coex_HeI_except_23S leaves that transition out so its
      ! 19.82 eV is counted once), and the Penning and associative quench
      ! deposit their heat at the transported population too
      ! (util_ion_eq, heat_chan 11).  A population carried from where it was
      ! made to where it is destroyed therefore carries its excitation
      ! energy in the reaction terms, and q_x added to the energy row
      ! would count it twice.  It is the same statement
      ! ionization_stage_enthalpy_face_flux makes of the ionization energy
      ! of the stages.  Its size, MEASURED on the certified LHS 1140 b
      ! validation state (energy-coupled He/H 2.09 closure, K_zz = 1e9
      ! cm^2 s^-1; md/He23S_transport_design_20261003.md section 11):
      ! max |div q_x| over the energy row's largest term is 3.0e-7 over the
      ! column (cell 198, r = 1.14 R_p) and 6.6e-9 in the gated wind
      ! r >= 1.2 R_p (cell 263, r = 1.45 R_p), so it would be negligible
      ! even as a term.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: divq
      real(dp), dimension(0:N), optional,       intent(out) :: q_out

      real(dp), dimension(1-Ng:N+Ng) :: Fnone, N_H, N_He, nHc, nHec, nd
      real(dp), dimension(1-Ng:N+Ng) :: xHeI_cell, xmeta_cell
      real(dp), dimension(0:N)       :: n_Hf, n_Hef, Kf, drf, q, fa
      real(dp), dimension(0:N)       :: xclose, Fclose
      real(dp), dimension(1:N)       :: cv
      real(dp), dimension(3,1-Ng:N+Ng) :: xHe
      real(dp), dimension(3,0:N)     :: xfHe, FkHe, EHe
      real(dp) :: eps3_erg, x3f, xIf
      integer  :: j

      divq = 0.0d0
      q    = 0.0d0
      if (present(q_out)) q_out = 0.0d0
      if (.not. carrier_solved(ic_HeTR)) return

      Fnone = 0.0d0
      call hydrogen_and_helium_nucleus_face_flux(rho, Tcode, f_sp, Fnone, &
                    N_H, N_He, n_Hf, n_Hef, Kf, drf,                      &
                    n_H_cell = nHc, n_He_cell = nHec)
      nd = rho*n0
      xHe(1,:) = f_sp(:,isp_HeII) *nd/max(nHec, 1.0d-300)
      xHe(2,:) = f_sp(:,isp_HeIII)*nd/max(nHec, 1.0d-300)
      xHe(3,:) = f_sp(:,isp_HeTR) *nd/max(nHec, 1.0d-300)
      call ionization_stage_face_flux(xHe, N_He, n_Hef, Kf, drf, xfHe,    &
                                      xclose, FkHe, Fclose, Ek_out = EHe)
      xmeta_cell = xHe(3,:)
      xHeI_cell  = f_sp(:,isp_HeI)*nd/max(nHec, 1.0d-300)
      eps3_erg = species_formation_energy(isp_HeTR)/erg2eV
      do j = 0, N
         x3f = 0.5d0*(xmeta_cell(j) + xmeta_cell(j+1))
         xIf = 0.5d0*(xHeI_cell(j) + xHeI_cell(j+1))
         if (xIf .le. 0.0d0) cycle
         q(j) = eps3_erg*(EHe(3,j) + (x3f/xIf)*(EHe(1,j) + EHe(2,j)))
      enddo
      call spherical_face_area_and_cell_volume(fa, cv)
      do j = 1, N
         divq(j) = (fa(j)*q(j) - fa(j-1)*q(j-1))/(cv(j)*R0*q0)
      enddo
      if (present(q_out)) q_out = q
      end subroutine helium_metastable_excitation_energy_divergence

      ! ------------------------------------------------------------- !

      ! Residual of the backward-Euler transport-chemistry system, and the
      ! face-flux derivatives the Jacobian needs.  rnorm is the largest
      ! RELATIVE row imbalance, built the way composition_residual builds
      ! its own: the row divided by the sum of the magnitudes of its terms.
      subroutine carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp,   &
                                  rep, msum, Frho, advect, Agrd, Bdrf,   &
                                  updrf,                                 &
                                  nH_free, nO_free, nC_free, sigrate,    &
                                  res, Jf, dJl, dJr, rnorm, rnorm_phys,   &
                                  absent_out)
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
      ! WHERE THE CARRIER IS ABSENT: the cells in which this carrier stands
      ! below carrier_absent_fraction of the free reservoir of its own
      ! element.  The row is still assembled and still measured there -- it
      ! is the GATE that reads this, and only in the certification, which
      ! reports such a row and does not let it decide a verdict.  Formed
      ! here because this is the one place that holds the carrier density
      ! and the reference density of its element together.
      logical, dimension(1:N,n_carrier_max), optional, intent(out) ::     &
                                                             absent_out
      ! Cell and species carrying the worst relative row imbalance, for the
      ! debug line: a Newton that stops short says WHERE it stopped.
      real(dp) :: nc(n_carrier_max), src(n_carrier_max)
      ! The production and the loss of each row separately, and the
      ! radiative part of the loss: the scale below is their sum and the
      ! dump reports them one by one.
      real(dp) :: srcp(n_carrier_max), srcl(n_carrier_max)
      real(dp) :: srcph(n_carrier_max)
      real(dp) :: srch2c(n_h2chan)
      real(dp) :: Kj, sL, sR, dsc, dph, rr, nfl, nref
      real(dp) :: rnormp
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real(dp), dimension(0:N) :: fa
      real(dp), dimension(1:N) :: cv
      real(dp) :: R0sq, R0cb
      integer  :: j, ic
      ! Row scale of every row of this assembly, kept so that the worst
      ! relative imbalance can be located AFTER the cell loop rather than
      ! inside it -- see there.
      real(dp) :: rsc(1:N,n_carrier_max)
      ! The advective term of every row and its size, formed once from the
      ! face fluxes before the cell loop: the face reconstruction is an
      ! operation on the whole column and not on a cell.
      real(dp) :: adv(1:N,n_carrier_max), advmag(1:N,n_carrier_max)
      ! The two advective face terms, filled only for the record.
      real(dp) :: advfin(1:N,n_carrier_max), advfout(1:N,n_carrier_max)

      call l22b_setup()
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
      adv     = 0.0d0
      advmag  = 0.0d0
      advfin  = 0.0d0
      advfout = 0.0d0
      if (advect) call carrier_advective_divergence(fc, msum, Frho,      &
                                                    adv, advmag,         &
                                                    advin = advfin,      &
                                                    advout = advfout)

      ! The row dump of this assembly (EXHALE_CARRIER_ROW_TERMS=1, default
      ! off).  The environment is read once for the run.
      if (.not. rowdump_asked) then
         rowdump_on    = carrier_row_terms_on()
         rowdump_asked = .true.
      endif
      if (rowdump_on .and. .not. allocated(rowdump_res)) then
         allocate(rowdump_dif(1:N,n_carrier_max),                        &
                  rowdump_adv(1:N,n_carrier_max),                        &
                  rowdump_prod(1:N,n_carrier_max),                       &
                  rowdump_loss(1:N,n_carrier_max),                       &
                  rowdump_phot(1:N,n_carrier_max),                       &
                  rowdump_scale(1:N,n_carrier_max),                      &
                  rowdump_res(1:N,n_carrier_max),                        &
                  rowdump_nc(1:N,n_carrier_max),                         &
                  rowdump_floor(1:N,n_carrier_max),                      &
                  rowdump_h2ch(1:N,n_h2chan),                            &
                  rowdump_fdif_in(1:N,n_carrier_max),                    &
                  rowdump_fdif_out(1:N,n_carrier_max),                   &
                  rowdump_fadv_in(1:N,n_carrier_max),                    &
                  rowdump_fadv_out(1:N,n_carrier_max),                   &
                  rowdump_frho_in(1:N), rowdump_frho_out(1:N))
      endif
      if (rowdump_on) then
         rowdump_dif  = 0.0d0
         rowdump_adv  = 0.0d0
         rowdump_prod = 0.0d0
         rowdump_loss = 0.0d0
         rowdump_phot = 0.0d0
         rowdump_scale = 0.0d0
         rowdump_res  = 0.0d0
         rowdump_nc   = 0.0d0
         rowdump_floor= 0.0d0
         rowdump_h2ch = 0.0d0
         rowdump_fdif_in  = 0.0d0
         rowdump_fdif_out = 0.0d0
         rowdump_fadv_in  = 0.0d0
         rowdump_fadv_out = 0.0d0
         rowdump_frho_in  = 0.0d0
         rowdump_frho_out = 0.0d0
      endif

      res = 0.0d0
      Jf  = 0.0d0
      dJl = 0.0d0
      dJr = 0.0d0
      if (present(absent_out)) absent_out = .false.
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         ! AN IONIZATION STAGE CARRIES ITS WHOLE FACE FLUX HERE, advective
         ! and eddy together, on its element's nucleus flux; a molecular
         ! carrier carries its diffusive, eddy and settling flux here and
         ! its material advection in adv(j,ic) above.  Both end in the same
         ! array, so the divergence below and the block assembly of
         ! solve_carriers are written once for the two kinds of row.
         if (carrier_is_nucleus_fraction(ic)) then
            call carrier_stage_face_flux(ic, fc, Jf(:,ic), dJl(:,ic),    &
                                         dJr(:,ic))
            cycle
         endif
         do j = 1, N-1
            call carrier_face_flux(fc(j,ic), fc(j+1,ic), wfac(j),        &
                                   wfac(j+1), Agrd(j,ic),                &
                                   Bdrf(j,ic), updrf(j,ic),              &
                                   Jf(j,ic), dJl(j,ic), dJr(j,ic))
         enddo
         ! Jf(0) and Jf(N) stay at their initialized zero: the two end
         ! faces of the column carry no diffusive, eddy or settling flux
         ! (carrier_face_coefficients states why).
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
      !$omp   private(j,ic,nc,src,srcp,srcl,srcph,srch2c,Kj,sL,sR,dsc,   &
      !$omp           dph,rr,nfl,nref)                                  &
      !$omp   reduction(max:rnorm,rnormp)
      do j = 1, N
         call carrier_cell_densities(fc(j,:), j, nrho(j), nc)
         call carrier_source(j, nc, nH_free(j), nO_free(j), src,         &
                             sprod = srcp, sloss = srcl, sphot = srcph,  &
                             sh2chan = srch2c)
         if (rowdump_on) rowdump_h2ch(j,:) = srch2c
         Kj = 1.0d0/(cv(j)*R0cb)
         sL = fa(j-1)*R0sq
         sR = fa(j)*R0sq
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
            ! THE TIME TERM IS WRITTEN WITH THE DENSITY THE UNKNOWN IS A
            ! FRACTION OF: rho n0 for a fraction per unit mass and the
            ! element's nucleus density for an ionization stage.  With the
            ! element continuity dn_el/dt = -div N_el holding, the
            ! nonconservative form n_el dx/dt + div F_stage is the
            ! conservative d n_stage/dt + div F_stage exactly, which is the
            ! same statement the element row is written on.
            nref = carrier_nucleus_reference(ic, j, nrho(j))
            rr  = nref*(fc(j,ic) - fc_old(j,ic))/dt_phys(j)
            dsc = abs(nref*fc(j,ic)/dt_phys(j))                          &
                + abs(nref*fc_old(j,ic)/dt_phys(j))
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
            ! THE REACTION TERMS OF THE PHYSICAL SCALE ARE THE ROW'S OWN
            ! PRODUCTION AND LOSS, srcp + srcl, as the declaration of
            ! carrier_accept_tol states the measure ("the transport flux
            ! divergence and the reaction production and loss").  Where a
            ! row hands out only its net (the helium stages, the proton of an
            ! atomic gas) srcp + srcl is |src| bit for bit and nothing moves.
            ! Where it hands out its channels the net is not one of the
            ! row's terms but their difference: the He 2^3S level at
            ! 1.23 R_p on the 2000-cell LHS 1140 b He/H 0.65 state has
            ! production and loss 41.25 cm^-3 s^-1 each, net 8.8e-7 and
            ! transport 8.8e-7 (MEASURED, md/Update_EXHALE_stage3.md
            ! section 94), so a scale of |J| + |net| stood 5e-8 below the
            ! row's terms and asked the coupled state for 4e-13 of its
            ! channel rates; any change of the background at the
            ! convergence level of the wind solve moved the measure across
            ! 1e-5 and back from pass to pass.  The Newton's own scale dsc
            ! keeps |src|: it sets the step and the stop of the iteration
            ! on the cell's own row, not the statement about the state.
            dph = dph + srcp(ic) + srcl(ic)
            ! ABSOLUTE FLOOR ON THE ROW SCALE, and it is not cosmetic. A
            ! relative residual with no floor asks a row whose species is
            ! 1e-28 of its element to balance to 1e-12 of ITSELF: measured on
            ! the HD 209458 b example, H2O jumped four decades from cell to
            ! cell at 1e-28 and the Newton spent its whole iteration budget
            ! there at a relative residual of 1 while every cell that
            ! carries a molecule was already at 1e-13. The floor is 1e-20 of
            ! the element the carrier belongs to, converted to the row's own
            ! volumetric-rate units: the rate that would move the carrier by
            ! 1e-20 of its element in one signal-crossing time.  A row is
            ! judged against it; it does not declare the row closed.
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
                                                    nC_free(j),         &
                                                    cbg_nHenuc(j))
            ! AND THE SAME NUMBER DECIDES WHICH CELLS COUNT AS ABSENT.  The
            ! floor above bounds the row's SCALE; this test reads the
            ! abundance alone.  A carrier at 1e-27 cm^-3 in a gas of 1e6 can
            ! read a relative imbalance of 1 while its production is a real,
            ! resolved rate, so the test does not make the row round-off: it
            ! removes that cell from the carrier gate (certification keeps
            ! its measure in the report) without asking whether production
            ! or incoming flux could raise the carrier.  That is a choice.
            if (present(absent_out))                                     &
               absent_out(j,ic) = (nc(ic) .lt. carrier_absent_fraction*nfl)
            dsc = dsc + carrier_absent_fraction*nfl                      &
                       *max(1.0d0/dt_phys(j), sigrate(j))
            ! THE PHYSICAL SCALE TAKES THE SIGNAL RATE ALONE.  The step form
            ! of the same floor grows without bound as the step shrinks, and
            ! a scale that grows like 1/dt is exactly what a statement about
            ! the returned state must not have; the cell's signal-crossing
            ! rate is the same floor written without a step (declaration of
            ! row_terms_phys).
            dph = dph + carrier_absent_fraction*nfl*sigrate(j)
            if (rowdump_on) then
               rowdump_dif(j,ic)  = Kj*(sR*Jf(j,ic) - sL*Jf(j-1,ic))
               rowdump_adv(j,ic)  = adv(j,ic)
               ! The two faces separately, each already carrying its area
               ! and the cell volume, so the two sum to the divergence
               ! above and a reader can see WHICH face and which sign.
               rowdump_fdif_in(j,ic)  = -Kj*sL*Jf(j-1,ic)
               rowdump_fdif_out(j,ic) =  Kj*sR*Jf(j,ic)
               rowdump_fadv_in(j,ic)  = advfin(j,ic)
               rowdump_fadv_out(j,ic) = advfout(j,ic)
               rowdump_frho_in(j)     = Frho(j-1)
               rowdump_frho_out(j)    = Frho(j)
               rowdump_prod(j,ic) = srcp(ic)
               rowdump_loss(j,ic) = srcl(ic)
               rowdump_phot(j,ic) = srcph(ic)
               rowdump_res(j,ic)  = rr
               rowdump_nc(j,ic)   = nc(ic)
               rowdump_floor(j,ic)= carrier_absent_fraction*nfl          &
                                   *sigrate(j)
               rowdump_scale(j,ic) = dph
            endif
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
      real(dp) :: Kj, sL, sR, cadv, dn, dnref
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real(dp), dimension(0:N) :: fa
      real(dp), dimension(1:N) :: cv
      real(dp) :: R0sq, R0cb
      ! The two face coefficients of the advective term, frozen with the
      ! wind over this solve.
      real(dp), dimension(1:N) :: advj, advm
      ! The banded derivative of the whole advective term,
      ! formed only when EXHALE_L22B_JAC_RECON is set.
      real(dp), dimension(1:N,n_carrier_max,-1:1) :: dadv
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
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
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
         if (carrier_rows_advect .and. l22b_jac_recon)                    &
            call l22b_advective_band_derivative(fc, msum, Frho, dadv)
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
         !$omp   private(j,ic,kc,Kj,sL,sR,cadv,nc,s0,s1,dn,dnref)
         do j = 1, N
            Kj = 1.0d0/(cv(j)*R0cb)
            sL = fa(j-1)*R0sq
            sR = fa(j)*R0sq
            do ic = 1, n_carrier
               bb(j,ic,ic) = carrier_nucleus_reference(ic, j, nrho(j))   &
                            /dt_phys(j)
               ! The identity row of a carrier this run does not solve: the
               ! diagonal above with a zero right-hand side gives it a zero
               ! update, which is what keeps it out of the system without
               ! taking it out of the block.
               if (.not. carrier_solved(ic)) cycle
               ! AN IONIZATION STAGE HAS A LIVE FLUX AT BOTH END FACES, so
               ! the two guards below are not its boundary rule.  Its face
               ! flux carries the material advection, which crosses the
               ! base face and the outer face; what vanishes there is its
               ! EDDY term alone, and ionization_stage_face_jacobian has
               ! already set that half to zero at f = 0 and f = N.  So the
               ! derivative of the end face is written and the boundary
               ! rule stays where the flux states it.
               !
               ! For a molecular carrier the guards ARE the rule: its face
               ! flux is diffusive, eddy and settling only, and those are
               ! zero at both ends (carrier_face_coefficients), so the
               ! entries they would write are zero anyway.
               !
               ! At j = 1 the base face's donor derivative lands on the
               ! diagonal because nothing below the base states an
               ! ionization fraction: the inner ghosts carry cell 1's own
               ! partition, which makes them a copy of the unknown and not
               ! data (the boundary paragraph of ionization_stage_transport
               ! and carrier_base_composition_imposed, which is false for
               ! every stage).
               if (j .lt. N) then
                  bb(j,ic,ic) = bb(j,ic,ic) + Kj*sR*dJl(j,ic)
                  cc(j,ic,ic) = cc(j,ic,ic) + Kj*sR*dJr(j,ic)
               else if (carrier_is_nucleus_fraction(ic)) then
                  ! The outer ghost is cell N continued (carrier_outflow_
                  ! ghost), so both derivatives of the outer face are
                  ! derivatives of the same unknown and land together on
                  ! the diagonal.
                  bb(N,ic,ic) = bb(N,ic,ic)                              &
                              + Kj*sR*(dJl(N,ic) + dJr(N,ic))
               endif
               ! j = N carries no diffusive entry: the outer face of the
               ! column carries no such flux (carrier_face_coefficients).
               if (j .gt. 1) then
                  bb(j,ic,ic) = bb(j,ic,ic) - Kj*sL*dJr(j-1,ic)
                  aa(j,ic,ic) = aa(j,ic,ic) - Kj*sL*dJl(j-1,ic)
               else if (carrier_is_nucleus_fraction(ic)) then
                  bb(1,ic,ic) = bb(1,ic,ic)                              &
                              - Kj*sL*(dJl(0,ic) + dJr(0,ic))
               endif
               ! An ionization stage's advection is inside the face flux
               ! above; the entries below are the molecular rows' material
               ! advection on the bulk face mass flux and are not its.
               if (carrier_is_nucleus_fraction(ic)) cycle
               ! THE ADVECTIVE ENTRIES ARE THE DONOR-CELL LINEARIZATION OF
               ! THE FACE-FLUX DIVERGENCE the row now carries.  With the
               ! face composition taken from the cell the face mass flux
               ! selects, the term of cell j is
               !
               !   adv(j) = [ A_+ F_rho(j) Y(don(j))
               !              - A_- F_rho(j-1) Y(don(j-1)) ] / dV_j
               !            x msum(j)/m_c ,   Y = m_c f_c/msum ,
               !
               ! so d adv(j)/d f_c(k) is the face coefficient times
               ! msum(j)/msum(don): **the carrier mass cancels**, between
               ! the msum/m_c the row carries and the m_c/msum of Y.  It was
               ! being divided by m_c here as well, which made every
               ! advective entry a factor m_c too small -- 2 for H2.  The
               ! expression is now the one that agrees with a
               ! central difference of carrier_advective_divergence itself,
               ! to 1e-8 on a state whose mass fraction is constant so that
               ! the reconstruction contributes nothing.  The limiter of
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
               ! WHICH GHOSTS ARE DATA AND WHICH ARE COPIES, and the
               ! boundary derivatives that follow.  A ghost the residual
               ! IMPOSES is a constant of the row and contributes no
               ! derivative; a ghost the residual COPIES from an interior
               ! cell contributes the derivative of that cell.  The rule is
               ! read off carrier_face_mass_fraction, which is what the
               ! residual uses, so the two cannot disagree:
               !
               !   outer ghosts   fc(N+1:) = fc(N)          -- a COPY, always
               !                  (set beside every trial, and by the
               !                   returned-state fill), so where the outer
               !                  face flows INWARD its donor is that copy
               !                  and d adv(N)/d f_c(N) is not zero.  It was
               !                  being dropped: the old test was
               !                  "j .lt. N", which wrote the neighbour entry
               !                  for every interior face and nothing at all
               !                  at j = N.
               !   inner ghosts   Y(1-Ng:0) = Y(1) unless the handoff states
               !                  this carrier's partition -- a COPY in the
               !                  first case and DATA in the second.  So an
               !                  INFLOWING base face contributes
               !                  d adv(1)/d f_c(1) exactly when the
               !                  composition is NOT imposed.  It was being
               !                  dropped in that case too.
               !
               ! Both entries land on the diagonal, which is where a copy of
               ! the cell's own unknown belongs.  The limiter of
               ! species_face_fraction and the second-order reconstruction
               ! stay out of the Jacobian as before; these are the FIRST
               ! ORDER donor-cell terms and nothing else.
               if (carrier_rows_advect .and. l22b_jac_recon) then
                  ! The reconstructed banded Jacobian: the advective entries are
                  ! a central difference of the face-flux divergence itself,
                  ! limiter and reconstruction included, restricted to the
                  ! band a block-tridiagonal matrix can hold.  What falls
                  ! outside the band is measured in l22b_recon_dropped.
                  if (j .gt. 1) aa(j,ic,ic) = aa(j,ic,ic) + dadv(j,ic,-1)
                  bb(j,ic,ic) = bb(j,ic,ic) + dadv(j,ic,0)
                  if (j .lt. N) cc(j,ic,ic) = cc(j,ic,ic) + dadv(j,ic,1)
               else if (carrier_rows_advect) then
                  cadv = advj(j)*msum(j)
                  if (Frho(j) .ge. 0.0d0) then
                     bb(j,ic,ic) = bb(j,ic,ic) + cadv/msum(j)
                  else if (j .lt. N) then
                     cc(j,ic,ic) = cc(j,ic,ic) + cadv/msum(j+1)
                  else
                     ! j = N and the outer face flows inward: the donor is
                     ! the outer ghost, which is a copy of cell N.
                     bb(N,ic,ic) = bb(N,ic,ic) + cadv/msum(N+1)
                  endif
                  cadv = advm(j)*msum(j)
                  if (j .eq. 1) then
                     if (Frho(0) .lt. 0.0d0) then
                        ! The base face flows outward: the donor is cell 1.
                        bb(1,ic,ic) = bb(1,ic,ic) + cadv/msum(1)
                     else if (.not. carrier_base_composition_imposed(ic)) &
                        then
                        ! It flows inward and nothing states the partition,
                        ! so the ghost carries cell 1's own mass fraction.
                        bb(1,ic,ic) = bb(1,ic,ic) + cadv/msum(1)
                     endif
                     ! (the remaining case -- inflow with the handoff
                     !  imposing the partition -- is genuine data and has no
                     !  derivative)
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
            call carrier_cell_densities(fc(j,:), j, nrho(j), nc)
            call carrier_source(j, nc, nH_free(j), nO_free(j), s0)
            do kc = 1, n_carrier
               if (.not. carrier_solved(kc)) cycle
               dn = carrier_source_derivative_step(kc, nc(kc),           &
                                                   nH_free(j),           &
                                                   nO_free(j),           &
                                                   nC_free(j),           &
                                                   cbg_nHenuc(j))
               ! THE PERTURBATION STAYS INSIDE THE ADMISSIBLE SET.  A
               ! helium stage at the sum rule cannot be raised: the two
               ! ionized stages together hold at most the helium available
               ! to them, and a forward difference taken past that bound
               ! differentiates the closure's behaviour outside the set
               ! and not the chemistry.  The difference is taken inward
               ! there, which is the same derivative of the same row, and
               ! only while the stage itself stays nonnegative.
               ! With the He 2^3S level carried it is the third state of
               ! the same simplex and enters the same sum.
               if (kc .eq. ic_HeII .or. kc .eq. ic_HeIII .or.            &
                   kc .eq. ic_HeTR) then
                  if (carrier_solved(ic_HeTR)) then
                     if (nc(ic_HeII) + nc(ic_HeIII) + nc(ic_HeTR) + dn   &
                         .gt. carrier_helium_available_to_stages(j)      &
                         .and. nc(kc) .gt. dn) dn = -dn
                  else if (nc(ic_HeII) + nc(ic_HeIII) + dn .gt.          &
                      carrier_helium_available_to_stages(j) .and.        &
                      nc(kc) .gt. dn) then
                     dn = -dn
                  endif
               endif
               nc(kc) = nc(kc) + dn
               call carrier_source(j, nc, nH_free(j), nO_free(j), s1)
               nc(kc) = nc(kc) - dn
               ! The chain rule closes on the COLUMN's own reference
               ! density: d src/d f_c(kc) = (d src/d n(kc)) times the
               ! density that column's fraction is taken against, which is
               ! the element nucleus density for an ionization stage.
               dnref = carrier_nucleus_reference(kc, j, nrho(j))
               do ic = 1, n_carrier
                  if (.not. carrier_solved(ic)) cycle
                  bb(j,ic,kc) = bb(j,ic,kc)                              &
                              - (s1(ic) - s0(ic))/dn*dnref
               enddo
            enddo
         enddo
         !$omp end parallel do

         ! The action of the assembled matrix against a central
         ! difference of the full residual, once per run, on the direction
         ! the keys name.  It reassembles the residual into arrays of its
         ! own, so res and the step below are the ones this iteration built.
         if (l22b_jac_action .and. .not. l22b_jac_action_done) then
            l22b_jac_action_done = .true.
            call l22b_jacobian_action_report(fc, fc_old, nrho, wfac,      &
                       dt_phys, rp, rep, msum, Frho, Agrd, Bdrf, updrf,   &
                       nH_free, nO_free, nC_free, sigrate, aa, bb, cc)
         endif

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
            enddo
            ! Only the OUTER ghosts mirror the interior, by the stated
            ! outflow continuation.  The lower ghost is the inflow
            ! reservoir, not a copy of the base cell: it carries the
            ! composition the handoff states, and the advective term
            ! reconstructs the base face from it.  Writing the base cell
            ! into it made the base condition degenerate into a
            ! zero-gradient one from the second relaxation pass on -- the
            ! marching loop hid that because its ionization sweep re-pins
            ! the ghost every step, and relax_photochemical_composition,
            ! which takes many transport steps between two sweeps, did not
            ! (sec. 158).
            call carrier_outflow_ghost(ftry)
            ! THE TRIAL IS RETURNED TO THE ADMISSIBLE SET BEFORE ITS ROWS
            ! ARE EVALUATED, by the same projection the limiter and the
            ! write-back go through.  The nonnegativity above is one face
            ! of that set; the other is the stage sum rule, and a trial
            ! outside it has a negative neutral remainder in its element.
            ! The chemistry rows would then be evaluated at a partition no
            ! cell can be in -- with the reservation of HeH+ and of the
            ! frozen He 2^3S level, at a negative ground singlet -- and
            ! what the old closure returned there was the clip and not the
            ! chemistry.  Scaling the stages toward the interior is the
            ! same straight-line construction the nonnegativity clip above
            ! and the element limiter of the returned state both use.
            call carrier_ionization_stage_projection(ftry)
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
               ! last trial is still taken -- rejection and retry do not
               ! exist -- so the outcome has to be
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
      call carrier_outflow_ghost(fc)
      ! THE ADMISSIBLE SET OF AN IONIZATION STAGE IS THE SIMPLEX, NOT A
      ! NUCLEUS BUDGET.  Its unknown is the fraction of its element's
      ! nuclei in that stage, so what bounds it is x_k >= 0 and
      ! sum_k x_k <= the share of the element the atomic stages partition,
      ! the remainder being the neutral stage; the free nucleus density
      ! that bounds a molecular carrier is already divided out of it.  That
      ! share is one for hydrogen, whose molecules are carriers of their
      ! own, and less than one for helium wherever HeH+ or the frozen
      ! He 2^3S level holds part of the element.  What the projection had
      ! to move is measured and reported (stage_simplex_sum_over_limit,
      ! stage_fraction_under_zero) rather than asserted.
      call carrier_ionization_stage_projection(fc)
      do j = 1, N+Ng
         call carrier_cell_densities(fc(j,:), j, nrho(j), nc)
         nc  = max(nc, 0.0d0)
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
         do ic = 1, n_carrier_max
            fc(j,ic) = nc(ic)                                            &
                      /max(carrier_nucleus_reference(ic, j, nrho(j)),    &
                           1.0d-99)
         enddo
      enddo
      end subroutine limit_to_element_budget

      ! ------------------------------------------------------------- !

      ! RETURN THE IONIZATION STAGE FRACTIONS OF EVERY CELL AND GHOST TO
      ! THE SIMPLEX, element by element, through the one projection of
      ! ionization_stage_transport.  Hydrogen carries x(H II) alone today,
      ! so the simplex is the interval [0,1] and the projection is the
      ! statement that a stage cannot hold more nuclei than its element
      ! has; with the helium stages carried it is the two-dimensional face
      ! x(He II) + x(He III) <= 1 - (n_HeH+ + n_He(2^3S))/n_He,nuc and the
      ! same routine returns it, the bound handed to it cell by cell.
      ! What the projection moves is recorded on the SAME measure the
      ! element-budget limiter records its own clamps on: the amount by
      ! which the stage sum exceeds the nuclei available to it, which for
      ! one carried stage of a wholly atomic element is x - 1, the relative
      ! overshoot of that stage above its element's nucleus density.  So a stage clamp and a carrier clamp are one
      ! number in the diagnostics and the cell is marked in the same array,
      ! which is what keeps the acceptance from judging a constrained row
      ! as an unconstrained one.
      subroutine carrier_ionization_stage_projection(fc)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      real(dp), dimension(1,1-Ng:N+Ng) :: x1
      real(dp), dimension(2,1-Ng:N+Ng) :: x2
      real(dp), dimension(3,1-Ng:N+Ng) :: x3s
      real(dp), dimension(1-Ng:N+Ng)   :: xmax
      real(dp) :: over
      integer  :: j
      ! The cell record this routine marks belongs to the limiter that
      ! normally calls it; a caller that projects on its own still has to
      ! find it there.
      if (.not. allocated(pct_cell_constrained)) then
         allocate(pct_cell_constrained(1:N+Ng))
         pct_cell_constrained = .false.
      endif
      ! Hydrogen: one carried stage, so the simplex is the interval [0,1].
      if (carrier_solved(ic_Hp)) then
         do j = 1, N+Ng
            over = fc(j,ic_Hp) - 1.0d0
            if (over .le. 0.0d0) cycle
            pct_worst_limit = max(pct_worst_limit, over)
            pct_cell_constrained(j) = .true.
            if (over .gt. limit_report) pct_limited = pct_limited + 1
         enddo
         x1(1,:) = fc(:,ic_Hp)
         call stage_simplex_projection(x1)
         fc(:,ic_Hp) = x1(1,:)
      endif
      ! Helium: two carried stages sharing ONE simplex, and the bound on
      ! their sum is the share of the element they may hold,
      !
      !    x(He II) + x(He III) <= 1 - (n_HeH+ + n_He(2^3S))/n_He,nuc,
      !
      ! the remainder being the neutral helium.  They are projected
      ! together because the constraint is on their sum: taking each one to
      ! [0,1] on its own would admit a cell whose two ionized stages hold
      ! more helium than the cell has.  The bound is below one wherever
      ! helium sits outside the atomic stages -- in HeH+, or in the frozen
      ! metastable level of He I -- and reserving it is what keeps the
      ! neutral closure of the chemistry rows and the remainder of the
      ! write-back nonnegative without either clipping
      ! (carrier_helium_available_to_stages).
      if (carrier_solved(ic_HeII) .and. carrier_solved(ic_HeIII)) then
         ! The bound is a statement about the composition of the cell, so
         ! the frozen background of the step has to stand before the
         ! stages can be projected.  Without it the helium held outside
         ! the atomic stages is unknown and the only bound available is
         ! one, which is the defect this reservation exists to remove; so
         ! the absence is named rather than filled in.
         if (.not. allocated(cbg_nHenuc)) then
            write(*,'(a)') ' (carrier_ionization_stage_projection) the'//&
               ' helium nucleus count of the step is not built; call'//  &
               ' carrier_state (or, in a test,'//                        &
               ' carrier_helium_background_set_for_test) first'
            error stop 1
         endif
         do j = 1-Ng, N+Ng
            xmax(j) = carrier_helium_available_to_stages(j)              &
                     /max(cbg_nHenuc(j), 1.0d-300)
         enddo
         ! WITH THE He 2^3S LEVEL CARRIED the simplex is three-dimensional,
         !    x(He II) + x(He III) + x3 <= 1 - n_HeH+/n_He,nuc,
         ! the level no longer reserved in the bound (it is one of the
         ! states being bounded) and the ground singlet the remainder.
         if (carrier_solved(ic_HeTR)) then
            do j = 1, N+Ng
               over = fc(j,ic_HeII) + fc(j,ic_HeIII) + fc(j,ic_HeTR)     &
                    - xmax(j)
               if (over .le. 0.0d0) cycle
               pct_worst_limit = max(pct_worst_limit, over)
               pct_cell_constrained(j) = .true.
               if (over .gt. limit_report) pct_limited = pct_limited + 1
            enddo
            x3s(1,:) = fc(:,ic_HeII)
            x3s(2,:) = fc(:,ic_HeIII)
            x3s(3,:) = fc(:,ic_HeTR)
            call stage_simplex_projection(x3s, xmax)
            fc(:,ic_HeII)  = x3s(1,:)
            fc(:,ic_HeIII) = x3s(2,:)
            fc(:,ic_HeTR)  = x3s(3,:)
         else
         do j = 1, N+Ng
            over = fc(j,ic_HeII) + fc(j,ic_HeIII) - xmax(j)
            if (over .le. 0.0d0) cycle
            pct_worst_limit = max(pct_worst_limit, over)
            pct_cell_constrained(j) = .true.
            if (over .gt. limit_report) pct_limited = pct_limited + 1
         enddo
         x2(1,:) = fc(:,ic_HeII)
         x2(2,:) = fc(:,ic_HeIII)
         call stage_simplex_projection(x2, xmax)
         fc(:,ic_HeII)  = x2(1,:)
         fc(:,ic_HeIII) = x2(2,:)
         endif
      endif
      end subroutine carrier_ionization_stage_projection

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
      real(dp) :: nHeIIn, nHeIIIn, nHeTRn
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
         ! The proton's unknown is a fraction per hydrogen NUCLEUS, so it
         ! is multiplied by the element's nucleus density and not by the
         ! density of the mass row (carrier_nucleus_reference).
         if (carrier_solved(ic_Hp))                                      &
            nHpn = fc(j,ic_Hp)                                           &
                  *carrier_nucleus_reference(ic_Hp, j, nrho(j))
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

         ! ---- helium, where its two ionized stages are carried -------
         ! The stages are written from the transported fractions, and the
         ! helium nuclei they did not take go to the neutral stage.  The
         ! helium sitting in HeH+ and in any further helium-bearing
         ! molecule is left where it is, for the reason the hydrogen share
         ! leaves it: rescaling that molecule would move a nucleus of the
         ! OTHER element with it.  It is therefore taken out of the
         ! remainder, from the same frozen count the simplex reserved it
         ! with, so the element total of the cell is exactly what it was.
         ! He 2^3S is a level inside He I.  Left to the sweep, it is
         ! rescaled by the same factor and keeps its share of the neutral
         ! helium, which is the right starting point for the sweep that
         ! re-solves the level a moment later; carried ("He 2^3S
         ! transport"), it is written from its own row below.
         !
         ! THE REMAINDER IS NONNEGATIVE ON THE ADMISSIBLE SET and is not
         ! clipped: the state handed here has been through the projection,
         ! which bounds the two stages by the helium available to them, so
         ! what is left is at least the frozen metastable's own population.
         ! A negative remainder would mean the cell had been made to hold
         ! more helium than it has, so it is MEASURED into
         ! carrier_helium_singlet_under_zero and left in the state rather
         ! than hidden by a floor while the two ionized stages stand.
         if (carrier_solved(ic_HeII)) then
            nHeIIn  = fc(j,ic_HeII)                                      &
                     *carrier_nucleus_reference(ic_HeII, j, nrho(j))
            nHeIIIn = fc(j,ic_HeIII)                                     &
                     *carrier_nucleus_reference(ic_HeIII, j, nrho(j))
            f_sp(j,isp_HeII)  = nHeIIn /nd
            f_sp(j,isp_HeIII) = nHeIIIn/nd
            rest = carrier_nucleus_reference(ic_HeII, j, nrho(j))        &
                 - nHeIIn - nHeIIIn - cbg_nHemol(j)
            ! THE CARRIED He 2^3S LEVEL IS WRITTEN FROM ITS OWN ROW, and
            ! He I, which holds it, takes the whole remainder: the level is
            ! not rescaled with He I (that would undo its transport), and
            ! the singlet is rest - n(2^3S), nonnegative on the
            ! three-dimensional simplex the projection returned.
            if (carrier_solved(ic_HeTR)) then
               nHeTRn = fc(j,ic_HeTR)                                    &
                       *carrier_nucleus_reference(ic_HeTR, j, nrho(j))
               f_sp(j,isp_HeTR) = nHeTRn/nd
               if (rest - nHeTRn .lt. 0.0d0)                             &
                  carrier_helium_singlet_under_zero =                    &
                     max(carrier_helium_singlet_under_zero,              &
                         (nHeTRn - rest)/max(cbg_nHenuc(j), 1.0d-300))
               f_sp(j,isp_HeI) = rest/nd
            else
            if (rest - cbg_nheiTR(j) .lt. 0.0d0)                         &
               carrier_helium_singlet_under_zero =                       &
                  max(carrier_helium_singlet_under_zero,                 &
                      (cbg_nheiTR(j) - rest)                             &
                      /max(cbg_nHenuc(j), 1.0d-300))
            held = f_sp(j,isp_HeI)*nd
            if (held .gt. 0.0d0 .and. rest .gt. 0.0d0) then
               sc = rest/held
               f_sp(j,isp_HeI) = f_sp(j,isp_HeI)*sc
               if (thereis_HeITR)                                        &
                  f_sp(j,isp_HeTR) = f_sp(j,isp_HeTR)*sc
            else
               ! The cell holds no neutral helium to share the remainder
               ! over: it becomes the whole neutral stage, and the level
               ! inside it starts the next sweep empty.
               f_sp(j,isp_HeI) = rest/nd
               if (thereis_HeITR) f_sp(j,isp_HeTR) = 0.0d0
            endif
            endif
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

      ! THE H2 CONTENT EACH CELL'S OWN ROW WOULD SETTLE AT, at this state,
      ! with transport left out: the density at which the row's chemical
      ! production equals its chemical loss, both taken from carrier_source,
      ! so this is the row's own chemistry and not a second statement of it
      ! -- photodissociation, photoionization, the ion channels and the
      ! thermal pair all included.
      !
      ! THE ROOT IS SOLVED FOR AND NOT DIVIDED OUT, and that is a physics
      ! statement about this row, not a numerical preference.  The loss of
      ! H2 is proportional to n(H2); its PRODUCTION is not independent of
      ! n(H2), because the dominant formation channel is the three-body
      ! association R15, k15 n(H I)^2, and carrier_source closes atomic
      ! hydrogen out of the element budget as n(H I) = n_H,avail - 2 n(H2)
      ! - ... .  The production therefore FALLS as n(H2) rises and vanishes
      ! where every hydrogen nucleus is already bound, so the row is
      ! quadratic in the unknown and P/(L/n) evaluated at the trial is not
      ! its root.  Evaluated at a fully molecular trial -- which is exactly
      ! what the thermochemical fit hands this routine in a cold base layer
      ! -- P/(L/n) collapses with n(H I)^2 and returns a value decades
      ! below the root.  MEASURED on LHS 1140 b: at the base of
      ! molecular_scalar_gj1132_kzz1e9/HeH2.13, P/(L/n) read x2 = 0.089 to
      ! 0.14 where the row's own root is x2 = 0.9 or above, and the whole
      ! "the handoff is 7 to 21 times the network's root" finding was that
      ! division and not a disagreement between the two
      ! chemistries.
      !
      ! HOW.  g(n2) = production(n2) - loss(n2) is evaluated through
      ! carrier_source itself, so the chemistry is stated once.  It is
      ! strictly decreasing -- every production term is non-increasing in
      ! n(H2) and every loss term is proportional to it -- and it changes
      ! sign between n2 = 0 and the element ceiling, so a bisection on that
      ! bracket converges to the one root.  The state is otherwise frozen:
      ! the ion stages, the molecular ions and the temperature are the
      ! sweep's, so this is the root of the carrier half of the
      ! alternation and not of the fully coupled system.
      !
      ! It is evaluated at the state it is handed, so the molecular ions
      ! that carry part of the production must already be in that state; an
      ! atomic state, whose H2+ and HeH+ are zero, returns the
      ! thermochemical association alone.
      !
      ! WHAT IT IS FOR: the molecular seed.  The thermochemical
      ! fit q_H2(p, T) is the balance of the three-body association against
      ! the thermal dissociation and knows nothing of the radiation field or
      ! of the ionized gas, so in an irradiated outer wind it asserts an
      ! equilibrium that does not hold there; MEASURED on LHS 1140 b it put
      ! 13 to 21 percent of the gas into H2 from 5.9 R_p out to 29 R_p, 80
      ! to 360 times this root, with atomic H four decades BELOW the
      ! molecule it is made from.  The seed takes the smaller of the two.
      !
      ! Zero is returned where the loss rate is zero, where the carrier is
      ! not solved, and on every cell if no sweep has filled the background.
      subroutine carrier_h2_chemical_root(rho, f_sp, n_root, row_imbalance)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: n_root
      ! HOW WELL THE RETURNED DENSITY BALANCES THE ROW, |P - L|/(P + L) at
      ! it: zero says the number handed back IS the root of the chemistry
      ! and not an estimate of it.  It is the one measurement that tells a
      ! solved root from a divided one, so the seed report carries it.
      real(dp), dimension(1-Ng:N+Ng), optional, intent(out) :: row_imbalance
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fc
      real(dp), dimension(1-Ng:N+Ng) :: ntot, nrho, wfac, TK, mbar
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp) :: nc(n_carrier_max), src(n_carrier_max)
      real(dp) :: sprod(n_carrier_max), sloss(n_carrier_max)
      real(dp) :: nlo, nhi, nmid, glo, ghi, gmid, ncap
      integer  :: j, it
      ! Bisection budget.  The bracket closes on a RELATIVE width, because
      ! a root far below the element ceiling -- the far wind, where x(H2)
      ! reaches 1e-15 of the hydrogen -- is not resolved at all by an
      ! absolute fraction of that ceiling.  Reaching 1e-12 of a root 1e-16
      ! of the ceiling takes about 93 halvings, so the count below is well
      ! clear of it and the loop always exits on the width.
      integer, parameter :: n_bisect = 200
      n_root = 0.0d0
      if (present(row_imbalance)) row_imbalance = 0.0d0
      if (.not. thereis_mol) return
      if (.not. bg_ready)    return
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK, mbar,       &
                         nH_free, nO_free, nC_free)
      call carrier_photolysis(rho, TK, f_sp)
      do j = 1, N
         call carrier_cell_densities(fc(j,:), j, nrho(j), nc)
         ! The element ceiling of this cell's H2: every hydrogen nucleus
         ! the carriers may hold, less the nuclei the other carriers hold,
         ! divided by the two nuclei an H2 molecule takes.  It is the same
         ! budget carrier_source closes atomic hydrogen with, so the
         ! production is exactly zero at the ceiling.
         ncap = hydrogen_available_to_carriers(j, nH_free(j))             &
                - nc(ic_OH) - 2.0d0*nc(ic_H2O)
         ! With the ionization state carried the proton is a carrier too and
         ! holds one nucleus each, subtracted by the caller rather than by
         ! the budget above -- the same subtraction carrier_source makes.
         if (ionization_transport) ncap = ncap - nc(ic_Hp)
         ncap = 0.5d0*max(ncap, 0.0d0)
         if (ncap .le. 0.0d0) cycle
         nlo = 0.0d0
         nhi = ncap
         nc(ic_H2) = nlo
         call carrier_source(j, nc, nH_free(j), nO_free(j), src,          &
                             sprod = sprod, sloss = sloss)
         glo = sprod(ic_H2) - sloss(ic_H2)
         ! No production at all with every nucleus free: the row has no
         ! root above zero and the cell carries no H2 chemically.
         if (glo .le. 0.0d0) cycle
         nc(ic_H2) = nhi
         call carrier_source(j, nc, nH_free(j), nO_free(j), src,          &
                             sprod = sprod, sloss = sloss)
         ghi = sprod(ic_H2) - sloss(ic_H2)
         ! The ceiling itself balances or still produces: the row wants
         ! every nucleus it can have, which is the ceiling.
         if (ghi .ge. 0.0d0) then
            n_root(j) = ncap
            if (present(row_imbalance))                                   &
               row_imbalance(j) = abs(sprod(ic_H2) - sloss(ic_H2))        &
                                 /max(sprod(ic_H2) + sloss(ic_H2),        &
                                      1.0d-300)
            cycle
         endif
         do it = 1, n_bisect
            nmid = 0.5d0*(nlo + nhi)
            nc(ic_H2) = nmid
            call carrier_source(j, nc, nH_free(j), nO_free(j), src,       &
                                sprod = sprod, sloss = sloss)
            gmid = sprod(ic_H2) - sloss(ic_H2)
            if (gmid .gt. 0.0d0) then
               nlo = nmid
               glo = gmid
            else
               nhi = nmid
               ghi = gmid
            endif
            if (nhi - nlo .le. 1.0d-12*max(nhi, 1.0d-300)) exit
         enddo
         n_root(j) = 0.5d0*(nlo + nhi)
         if (present(row_imbalance)) then
            nc(ic_H2) = n_root(j)
            call carrier_source(j, nc, nH_free(j), nO_free(j), src,       &
                                sprod = sprod, sloss = sloss)
            row_imbalance(j) = abs(sprod(ic_H2) - sloss(ic_H2))           &
                              /max(sprod(ic_H2) + sloss(ic_H2), 1.0d-300)
         endif
      enddo
      n_root(1-Ng:0)   = n_root(1)
      n_root(N+1:N+Ng) = n_root(N)
      if (present(row_imbalance)) then
         row_imbalance(1-Ng:0)   = row_imbalance(1)
         row_imbalance(N+1:N+Ng) = row_imbalance(N)
      endif
      end subroutine carrier_h2_chemical_root

      ! ------------------------------------------------------------- !

      ! WHAT THE FIXED WIND FEELS OF A TRIAL, cell by cell: the relative
      ! change of the cell's PARTICLE COUNT, (n_tot + n_e), which at a fixed
      ! conserved state is the same statement as the relative change of
      ! 1/mu.  This is the quantity the movement bound is written on.
      !
      ! WHY THIS QUANTITY AND NOT THE CARRIER'S OWN FRACTION.  The bound
      ! exists because the wind is held fixed while the carriers relax, so
      ! one pass may move the composition only as far as the wind's response
      ! to it stays linear.  What the wind responds to is the particle count
      ! and the mean molecular mass of the cell -- they are what set the
      ! pressure at the conserved thermal energy -- and not how much a
      ! carrier changed relative to some other cell's abundance.  An H2 at
      ! 1e-08 of the gas moves neither and is nothing the wind can feel; an
      ! H2 at tens of percent moves both.  A bound on (n_tot + n_e) is
      ! therefore loose exactly where the carrier is negligible and tight
      ! exactly where it is not, with no threshold beyond `trust` itself.
      !
      ! WHAT THE OLD MEASURE DID, MEASURED.  It compared the
      ! largest ABSOLUTE change of any carrier column against the largest H2
      ! mixing ratio of the entry state, one number for the whole grid.  On
      ! LHS 1140 b that reference is the base value, so the allowance was an
      ! absolute 1.9e-03 anywhere, and the bound was refused on every pass at
      ! a cell of the base or the front and NEVER at the cells the
      ! certification is about: cell 3 (r = 1.0006, change 7.4e-05 of an
      ! entry 6.8e-03) on the seed whose base is the row's own root, cell 143
      ! (r = 1.0543, change 5.97e-04 of an entry 4.73e-03) on the seed whose
      ! base is the handoff's.  The wind cells that carry the refusing row
      ! changed by 1.4e-08 a pass, four decades inside the allowance, and
      ! were passengers on whatever trial length the saturating cell left.
      ! EXHALE_CARRIER_BOUND_FRACTION=1 restores that measure.
      subroutine carrier_particle_count_change(nsum_now, nsum_entry, jw, dw)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: nsum_now, nsum_entry
      integer,  intent(out) :: jw
      real(dp), intent(out) :: dw
      real(dp) :: d
      integer  :: j
      dw = 0.0d0
      jw = 0
      do j = 1, N
         if (nsum_entry(j) .le. 0.0d0) cycle
         d = abs(nsum_now(j) - nsum_entry(j))/nsum_entry(j)
         if (d .gt. dw) then
            dw = d
            jw = j
         endif
      enddo
      end subroutine carrier_particle_count_change

      ! ------------------------------------------------------------- !

      ! EXHALE_CARRIER_BOUND_FRACTION=1 puts the movement bound back on the
      ! carrier fraction against the column's largest H2 mixing ratio, the
      ! measure used earlier.  Default off; the reasoning is at
      ! carrier_particle_count_change.
      logical function carrier_bound_on_fraction()
      character(len=8) :: env
      call get_environment_variable('EXHALE_CARRIER_BOUND_FRACTION', env)
      carrier_bound_on_fraction = (trim(env) .eq. '1')
      end function carrier_bound_on_fraction

      ! ------------------------------------------------------------- !

      ! READ THE DIAGNOSTIC INTERVAL MASK ONCE, AND SAY WHAT IT IS.
      ! The declaration of l22_mask_mode carries what this mask measures and
      ! what it does not certify.  Called at the head of every relaxation;
      ! the file is read on the first call only, so the mask is frozen for
      ! the whole comparison, which is the condition the experiment is read
      ! under.  With EXHALE_L22_MASK unset the mask is off and every path
      ! below is the unmasked one.
      subroutine l22_interval_mask_setup()
      character(len=512) :: fname
      character(len=16)  :: mode
      character(len=256) :: line
      integer :: u_msk, ios, jc, j
      if (l22_mask_ready) return
      l22_mask_ready = .true.
      l22_mask_mode  = l22_mask_off
      call get_environment_variable('EXHALE_L22_MASK', fname)
      if (len_trim(fname) .eq. 0) return
      call get_environment_variable('EXHALE_L22_MASK_MODE', mode)
      select case (trim(mode))
      case ('veto');  l22_mask_mode = l22_mask_veto
      case ('waive'); l22_mask_mode = l22_mask_waive
      case default
         write(*,'(A)') '    (L22 interval mask) EXHALE_L22_MASK_MODE'//  &
              ' must be veto or waive; the mask is off'
         return
      end select
      if (.not. allocated(l22_mask_omit)) allocate(l22_mask_omit(1:N))
      l22_mask_omit = .false.
      open(newunit = u_msk, file = trim(fname), status = 'old',           &
           action = 'read', iostat = ios)
      if (ios .ne. 0) then
         write(*,'(A,A)') '    (L22 interval mask) cannot read ',         &
              trim(fname)
         l22_mask_mode = l22_mask_off
         return
      endif
      do
         read(u_msk,'(A)',iostat = ios) line
         if (ios .ne. 0) exit
         line = adjustl(line)
         if (len_trim(line) .eq. 0) cycle
         if (line(1:1) .eq. '#') cycle
         read(line,*,iostat = ios) jc
         if (ios .ne. 0) cycle
         if (jc .ge. 1 .and. jc .le. N) l22_mask_omit(jc) = .true.
      enddo
      close(u_msk)
      l22_mask_n_omit = 0
      do j = 1, N
         if (l22_mask_omit(j)) l22_mask_n_omit = l22_mask_n_omit + 1
      enddo
      l22_mask_n_select = N - l22_mask_n_omit
      write(*,'(A)') '    (L22 interval mask) THIS RUN IS A DIAGNOSTIC'// &
           ' EXPERIMENT: the trial interval of the carrier relaxation is'
      write(*,'(A,A)') '      chosen on a subset of the column.  Mask '// &
           'file ', trim(fname)
      write(*,'(A,A,A,I0,A,I0,A)') '      mode ', trim(mode),             &
           ', omitted from the interval selection ', l22_mask_n_omit,     &
           ' cell(s), selection set ', l22_mask_n_select, ' cell(s)'
      if (l22_mask_mode .eq. l22_mask_veto)                               &
         write(*,'(A)') '      veto: the particle-count bound is still'// &
              ' checked on the omitted cells, so the accepted interval'// &
              ' is the unmasked one'
      if (l22_mask_mode .eq. l22_mask_waive)                              &
         write(*,'(A)') '      waive: the particle-count bound is NOT'//  &
              ' checked on the omitted cells; positivity, element and'//  &
              ' charge feasibility, finite thermodynamics and the'//      &
              ' factor-two safety stop hold the trial'
      write(*,'(A)') '      no cell is excluded from the certification'
      end subroutine l22_interval_mask_setup

      ! ------------------------------------------------------------- !

      ! THE RELATIVE PARTICLE-COUNT CHANGE OF A TRIAL ON THE TWO SETS THE
      ! DIAGNOSTIC MASK DEFINES: the selection set, which chooses the
      ! interval, and the omitted set, which does not.  Same quantity as
      ! carrier_particle_count_change, split by the mask.
      subroutine l22_particle_count_change_by_set(nsum_now, nsum_entry,   &
                                        jsel, dsel, jomit, domit)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: nsum_now, nsum_entry
      integer,  intent(out) :: jsel, jomit
      real(dp), intent(out) :: dsel, domit
      real(dp) :: d
      integer  :: j
      jsel  = 0;  dsel  = 0.0d0
      jomit = 0;  domit = 0.0d0
      do j = 1, N
         if (nsum_entry(j) .le. 0.0d0) cycle
         d = abs(nsum_now(j) - nsum_entry(j))/nsum_entry(j)
         if (l22_mask_omit(j)) then
            if (d .gt. domit) then
               domit = d;  jomit = j
            endif
         else
            if (d .gt. dsel) then
               dsel = d;  jsel = j
            endif
         endif
      enddo
      end subroutine l22_particle_count_change_by_set

      ! ------------------------------------------------------------- !

      ! The diagnostic keys, read once for the run.  With none of them set
      ! every switch below stays at the value it is declared with and the
      ! operator is unchanged; the declarations carry what each one states.
      subroutine l22b_setup()
      character(len=512) :: env
      integer :: ios, lo, hi
      if (l22b_ready) return
      l22b_ready = .true.
      call get_environment_variable('EXHALE_L22B_JAC_RECON', env)
      l22b_jac_recon = (trim(env) .eq. '1')
      call get_environment_variable('EXHALE_L22B_DISPLACEMENT', env)
      l22b_displacement = (trim(env) .eq. '1')
      call get_environment_variable('EXHALE_L22B_JAC_ACTION', env)
      l22b_jac_file = trim(env)
      l22b_jac_action = (len_trim(env) .gt. 0)
      if (l22b_jac_action) then
         l22b_jac_lo = 1
         l22b_jac_hi = N
         call get_environment_variable('EXHALE_L22B_JAC_CELLS', env)
         if (len_trim(env) .gt. 0) then
            read(env,*,iostat=ios) lo, hi
            if (ios .eq. 0) then
               l22b_jac_lo = max(1, min(N, lo))
               l22b_jac_hi = max(l22b_jac_lo, min(N, hi))
            endif
         endif
      endif
      if (l22b_jac_recon .or. l22b_jac_action) then
         write(*,'(A)') '    (L22b) THIS RUN IS A DIAGNOSTIC'//           &
              ' EXPERIMENT: the carrier operator or its Jacobian is not'
         write(*,'(A)') '      the one the tree integrates.'
         if (l22b_jac_recon)                                              &
            write(*,'(A)') '      advective Jacobian: a central'//        &
                 ' difference of the full face-flux divergence,'//        &
                 ' restricted to the tridiagonal band'
         if (l22b_jac_action)                                             &
            write(*,'(A,I0,A,I0,A,A)') '      Jacobian action probe on'// &
                 ' cells ', l22b_jac_lo, ' to ', l22b_jac_hi, ' -> ',     &
                 trim(l22b_jac_file)
      endif
      end subroutine l22b_setup

      ! ------------------------------------------------------------- !

      ! d adv(j,ic) / d f_c(k,ic) for k = j-1, j, j+1, by a central
      ! difference of carrier_advective_divergence, which is the FULL
      ! operator: the limited slopes of species_face_fraction and the
      ! composition bound it imposes are inside it.  The advective term of
      ! carrier ic reads only that carrier's own column, so one perturbed
      ! column gives every carrier's entries at once.
      !
      ! The reconstruction stencil reaches cell j-2 as well, and a
      ! block-tridiagonal matrix cannot hold that entry.  Its largest size
      ! relative to the band of the same row is recorded in
      ! l22b_recon_dropped, so the omission is measured and not assumed
      ! small.  The outer ghosts follow cell N in every perturbed column,
      ! which is the rule the residual uses.
      subroutine l22b_advective_band_derivative(fc, msum, Frho, dadv)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: msum
      real(dp), dimension(1-Ng:N+Ng),               intent(in)  :: Frho
      real(dp), dimension(1:N,n_carrier_max,-1:1),  intent(out) :: dadv
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fp, fm
      real(dp), dimension(1:N,n_carrier_max) :: ap, am, dmag
      ! The sum of the magnitudes of the entries of each row that fall
      ! inside the tridiagonal band and outside it, so that what the band
      ! cannot hold is measured against what it does hold.
      real(dp), dimension(1:N,n_carrier_max) :: bandsum, outsum
      real(dp), dimension(n_carrier_max) :: h, fref
      real(dp) :: d
      integer  :: j, k, ic, jd
      dadv = 0.0d0
      bandsum = 0.0d0
      outsum  = 0.0d0
      l22b_recon_dropped = 0.0d0
      do ic = 1, n_carrier
         fref(ic) = 0.0d0
         if (.not. carrier_solved(ic)) cycle
         fref(ic) = maxval(fc(1:N,ic))
      enddo
      do k = 1, N
         fp = fc
         fm = fc
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            h(ic) = 1.0d-6*max(fc(k,ic), 1.0d-12*fref(ic))
            if (h(ic) .le. 0.0d0) h(ic) = 1.0d-30
            fp(k,ic) = fc(k,ic) + h(ic)
            fm(k,ic) = max(fc(k,ic) - h(ic), 0.0d0)
         enddo
         call carrier_outflow_ghost(fp)
         call carrier_outflow_ghost(fm)
         call carrier_advective_divergence(fp, msum, Frho, ap, dmag)
         call carrier_advective_divergence(fm, msum, Frho, am, dmag)
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            do j = 1, N
               d = (ap(j,ic) - am(j,ic))                                  &
                  /max(fp(k,ic) - fm(k,ic), 1.0d-300)
               jd = k - j
               if (jd .ge. -1 .and. jd .le. 1) then
                  dadv(j,ic,jd) = d
                  bandsum(j,ic) = bandsum(j,ic) + abs(d)
               else
                  outsum(j,ic) = outsum(j,ic) + abs(d)
               endif
            enddo
         enddo
      enddo
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            if (bandsum(j,ic) .gt. 0.0d0) l22b_recon_dropped =            &
               max(l22b_recon_dropped, outsum(j,ic)/bandsum(j,ic))
         enddo
      enddo
      end subroutine l22b_advective_band_derivative

      ! ------------------------------------------------------------- !

      ! THE ACTION OF THE ASSEMBLED JACOBIAN AGAINST THE ACTION OF THE
      ! OPERATOR, row by row, on a direction concentrated on the cells the
      ! keys name.  The direction is each cell's own carrier fraction there
      ! and zero elsewhere, and the outer ghosts follow cell N, which is
      ! how the line search fills them, so the two actions are taken on one
      ! rule.  Nothing here changes the solve: the report writes a file and
      ! returns.  The row scales row_terms are left describing the probe's
      ! own evaluation rather than the iterate's, which is one reason the
      ! key is a diagnostic one.
      subroutine l22b_jacobian_action_report(fc, fc_old, nrho, wfac,      &
                 dt_phys, rp, rep, msum, Frho, Agrd, Bdrf, updrf,         &
                 nH_free, nO_free, nC_free, sigrate, aa, bb, cc)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: nrho, wfac, dt_phys
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: rp, rep, msum, Frho
      real(dp), dimension(0:N,n_carrier_max), intent(in) :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max), intent(in) :: updrf
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: nH_free, nO_free
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: nC_free, sigrate
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max), intent(in) :: &
                                                             aa, bb, cc
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: vdir, fp, fm
      real(dp), dimension(1:N,n_carrier_max) :: resp, resm, Jv, Fd
      real(dp), dimension(0:N,n_carrier_max)  :: Jfl, dJll, dJrl
      real(dp) :: rl, rlp, hstep, num, den, relerr, worst
      integer  :: j, ic, kc, u, jworst, icworst
      vdir = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = l22b_jac_lo, l22b_jac_hi
            vdir(j,ic) = fc(j,ic)
         enddo
         vdir(N+1:N+Ng,ic) = vdir(N,ic)
      enddo
      hstep = 1.0d-6
      fp = fc
      fm = fc
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            fp(j,ic) = fc(j,ic) + hstep*vdir(j,ic)
            fm(j,ic) = fc(j,ic) - hstep*vdir(j,ic)
         enddo
      enddo
      call carrier_outflow_ghost(fp)
      call carrier_outflow_ghost(fm)
      call carrier_residual(fp, fc_old, nrho, wfac, dt_phys, rp, rep,     &
                            msum, Frho, carrier_rows_advect, Agrd, Bdrf,  &
                            updrf, nH_free, nO_free, nC_free, sigrate,    &
                            resp, Jfl, dJll, dJrl, rl, rlp)
      call carrier_residual(fm, fc_old, nrho, wfac, dt_phys, rp, rep,     &
                            msum, Frho, carrier_rows_advect, Agrd, Bdrf,  &
                            updrf, nH_free, nO_free, nC_free, sigrate,    &
                            resm, Jfl, dJll, dJrl, rl, rlp)
      Fd = (resp - resm)/(2.0d0*hstep)
      Jv = 0.0d0
      do j = 1, N
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            do kc = 1, n_carrier
               if (.not. carrier_solved(kc)) cycle
               ! The ghost's own derivative is already on the diagonal of
               ! row N, so the superdiagonal is read only for j < N.
               if (j .gt. 1) Jv(j,ic) = Jv(j,ic)                          &
                                      + aa(j,ic,kc)*vdir(j-1,kc)
               Jv(j,ic) = Jv(j,ic) + bb(j,ic,kc)*vdir(j,kc)
               if (j .lt. N) Jv(j,ic) = Jv(j,ic)                          &
                                      + cc(j,ic,kc)*vdir(j+1,kc)
            enddo
         enddo
      enddo
      worst = 0.0d0;  jworst = 0;  icworst = 1
      open(newunit=u, file=trim(l22b_jac_file), status='replace',         &
           action='write')
      write(u,'(A)') '# L22 step 2b: the assembled Jacobian action'//     &
                     ' against a central difference of the full residual'
      write(u,'(A,I0,A,I0,A,ES12.5)') '# direction: the carrier'//        &
           ' fraction of cells ', l22b_jac_lo, ' to ', l22b_jac_hi,       &
           ', step ', hstep
      write(u,'(A,I0)') '# advective Jacobian with the reconstruction'//  &
           ' terms: ', merge(1, 0, l22b_jac_recon)
      write(u,'(A)') '# columns: cell carrier J_action fd_action'//       &
                     ' abs_diff rel_error'
      do j = 1, N
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            num = abs(Jv(j,ic) - Fd(j,ic))
            den = max(abs(Fd(j,ic)), abs(Jv(j,ic)))
            relerr = 0.0d0
            if (den .gt. 0.0d0) relerr = num/den
            write(u,'(I6,1X,A8,4(1X,ES14.7))') j,                         &
                 trim(carrier_name(ic)), Jv(j,ic), Fd(j,ic), num, relerr
            if (den .gt. 0.0d0 .and. relerr .gt. worst) then
               worst = relerr;  jworst = j;  icworst = ic
            endif
         enddo
      enddo
      close(u)
      write(*,'(A,ES10.3,A,I0,A,A)') '    (L22b) Jacobian action:'//      &
           ' worst relative error ', worst, ' at cell ', jworst,          &
           ' carrier ', trim(carrier_name(icworst))
      if (l22b_jac_recon) write(*,'(A,ES10.3)') '    (L22b) largest'//    &
           ' advective entry outside the tridiagonal band, relative to'// &
           ' the band of its row: ', l22b_recon_dropped
      end subroutine l22b_jacobian_action_report

      ! ------------------------------------------------------------- !

      ! THE LARGEST COMPOSITION CHANGE OF A TRIAL AND WHERE IT STANDS.
      ! It is the movement bound of the relaxation only under
      ! EXHALE_CARRIER_BOUND_FRACTION=1 (carrier_bound_on_fraction): a trial
      ! is then refused when this change, divided by the largest H2 mixing
      ! ratio of the entry state, exceeds trust.  The DEFAULT bound is the
      ! relative change of the particle count n_tot + n_e of each cell
      ! (carrier_particle_count_change, relax_photochemical_composition).
      ! The reasoning for the absolute form follows.  What the wind responds to
      ! is the ABSOLUTE composition change -- it is what moves the mean
      ! molecular mass, the particle count and the equation of state -- and
      ! an absolute bound is therefore already vacuous for a cell whose
      ! carrier is decades below the column maximum: such a cell is not what
      ! the bound holds.  (Admitting every cell an e-fold of its own value
      ! on top of it was tried and MEASURED: the extra clause
      ! is looser exactly where x_j is LARGE, so it freed the base and the
      ! front, which is the opposite of what a bound on the wind's linear
      ! response is for, and it made trust inert over three decades on the
      ! carrier_retry column.  It was removed.)
      !
      ! WHAT IS RETURNED BESIDE IT: the cell and the carrier that attain the
      ! change, so that a pass refused on the bound can say which cell
      ! refused it instead of only that one did.
      subroutine carrier_worst_composition_change(fa, fb, jw, icw, dmax)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fa, fb
      integer,  intent(out) :: jw, icw
      real(dp), intent(out) :: dmax
      real(dp) :: d
      integer  :: j, ic
      dmax = 0.0d0
      jw   = 0
      icw  = 1
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            d = abs(fa(j,ic) - fb(j,ic))
            if (d .gt. dmax) then
               dmax = d
               jw   = j
               icw  = ic
            endif
         enddo
      enddo
      end subroutine carrier_worst_composition_change

      ! ------------------------------------------------------------- !

      ! The largest change of a solved carrier between two carrier states,
      ! each carrier against its OWN scale xref (the largest fraction that
      ! carrier holds), over the physical cells: max over carriers of
      ! max_j |fa - fb|/xref.  The fixed-point test of a relaxation pass;
      ! relax_photochemical_composition states why each carrier is taken on
      ! its own scale.
      real(dp) function carrier_displacement_on_own_scale(fa, fb, xref)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fa, fb
      real(dp), dimension(n_carrier_max), intent(in) :: xref
      real(dp) :: d
      integer  :: j, ic
      d = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            d = max(d, abs(fa(j,ic) - fb(j,ic))/xref(ic))
         enddo
      enddo
      carrier_displacement_on_own_scale = d
      end function carrier_displacement_on_own_scale

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
      case (carrier_relax_mask_safety_stop)
         txt = 'the diagnostic interval mask hit its safety stop'
      case default
         txt = 'undefined'
      end select
      end function carrier_relax_outcome_text

      ! ------------------------------------------------------------- !

      ! Why one trial of the relaxation was undone, in words, for the
      ! trial-by-trial log of a pass (EXHALE_CARRIER_DEBUG=1).  The codes
      ! are the local refused_* parameters of
      ! relax_photochemical_composition.
      function carrier_trial_refusal_text(refusal) result(txt)
      integer, intent(in) :: refusal
      character(len=40) :: txt
      select case (refusal)
      case (1);    txt = 'the interval was not covered'
      case (2);    txt = 'the movement bound'
      case (3);    txt = 'the chemistry did not close'
      case default; txt = 'unrecognized'
      end select
      end function carrier_trial_refusal_text

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
      ! (densities, p from rho e, T, sweep) is therefore repeated until the
      ! relative temperature increment falls below chem_cycle_tol = 1e-5,
      ! which is the one tolerance the exit tests.
      !
      ! WHAT THE STATE HANDED BACK IS, exactly.  The composition is the one
      ! the last sweep produced at the temperature of the cycle BEFORE the
      ! exit; the temperature beside it is the one that composition and the
      ! unchanged thermal energy give AFTER it.  The two therefore stand one
      ! increment apart, and the contract bounds that increment by
      ! chem_cycle_tol alone, i.e. by ten times ieq_res_tol.  MEASURED on
      ! the twelve-cell molecular column: the realized
      ! gap there is 2.4e-7 relative, well inside what the contract admits.
      !
      ! WHAT THAT GAP DOES TO THE CHEMISTRY, measured and not assumed.  A
      ! temperature change below the sweep's reaction-residual tolerance
      ! does NOT bound the abundance movement of every species by that
      ! tolerance: on the same column a displacement of ieq_res_tol = 1e-6
      ! moves a trace carrier (HeH+, 1.3e-12 of the gas) by 4.8e-6, five
      ! times the tolerance, and one of chem_cycle_tol moves it by 4.5e-5,
      ! 45 times (table 6).  What the gap does bound is the quantity the
      ! tolerance measures: the reaction residual of the returned
      ! composition AT THE RETURNED TEMPERATURE is 2.6e-14 over the column
      ! and 1.6e-16 at the binding cell, eleven decades below ieq_res_tol,
      ! so the composition handed back is a root of the network at the
      ! temperature it is handed back with.  The two readings differ
      ! because the normalized residual of a row is proportional to the
      ! abundance that row is written on: at the HeH+ row of that cell the
      ! sensitivity of the residual to its own unknown is |R|/eps =
      ! 5.5e-13, so a residual of ieq_res_tol admits a RELATIVE HeH+ error
      ! of 1.8e+6 (table 5).  A composition bound is therefore taken from
      ! the residual only with that conditioning beside it.
      !
      ! ok is false when a sweep left a cell non-finite; the caller then
      ! discards the step and the background this sweep wrote.
      subroutine equilibrate_chemistry_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, heat, cool, eta, ok,   &
                                             n_cycles, reason, increment)
      ! THE THERMOCHEMICAL CLOSURE AT A FIXED CONSERVED STATE, and its
      ! contract. Cycles of (p, T of the composition at the unchanged
      ! thermal energy; one equilibrium sweep; p, T again) are taken until
      ! the temperature stops moving. The closure is reported as reached,
      ! ok = .true., ONLY when all of the following hold on the last cycle:
      ! the relative temperature increment is below chem_cycle_tol; the
      ! sweep's ledger reports no nonfinite cell and no lower-boundary solve
      ! left open; the composition is a set
      ! of numbers; and p, T, heat, cool and eta are finite on the physical
      ! cells with p > 0 and T > 0. Anything else -- the cycle budget spent,
      ! a nonfinite composition, an unsolved lower boundary, a state that
      ! is not admissible -- is
      ! ok = .false. with the reason named, and the caller restores the
      ! trial. EVERY ONE OF THOSE IS A STATEMENT ABOUT THE STATE HANDED
      ! BACK, which is what the closure is asked about; the count of cells
      ! whose ROOT SEARCH left the element simplex is reported and decides
      ! nothing. Until 2026-09-13 ok started true and
      ! the exhausted loop fell through with it (the review of
      ! 2026-09-12): an algebraically consistent (p, T) beside a chemistry
      ! evaluated at an earlier temperature was handed back as closed.
      real(dp), dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out)   :: p, T
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: heat, cool
      real(dp), dimension(1-Ng:N+Ng),           intent(inout) :: eta
      logical,                                  intent(out)   :: ok
      integer,                                  intent(out)   :: n_cycles
      integer,  intent(out), optional :: reason
      real(dp), intent(out), optional :: increment

      real(dp), dimension(1-Ng:N+Ng) :: ntot, ne, T_prev
      type(ioniz_eq_ledger) :: ledger
      real(dp) :: dT
      integer  :: k, why, cap, j_dT

      ok       = .false.
      n_cycles = 0
      why      = chem_closure_exhausted
      dT       = huge(1.0d0)
      cap      = chem_cycles_max
      if (chem_cycles_cap_for_test .ge. 0) cap = chem_cycles_cap_for_test
      call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,     &
                                                             p, T, ntot, ne)
      do k = 1, cap
         T_prev = T
         ! The caller restores the trial whenever this closure is not
         ! reached, so the sweep may hand back a lower boundary it could not
         ! solve with that counted instead of stopping the run.
         call set_ioniz_eq_state_may_be_refused(.true.)
         call ioniz_eq(T, u(1,:), f_sp, heat, cool, eta, ledger)
         call set_ioniz_eq_state_may_be_refused(.false.)
         n_cycles = k
         chem_last_offsimplex  = ledger%n_offsimplex
         chem_last_nonfinite   = ledger%n_nonfinite
         chem_last_mol_clamped = ledger%n_mol_clamped
         chem_last_viol_worst  = ledger%viol_worst
         if (ledger%n_nonfinite .gt. 0 .or.                               &
             .not. all(f_sp(1:N,:) .eq. f_sp(1:N,:)) .or.                 &
             .not. all(abs(f_sp(1:N,:)) .le. huge(1.0d0))) then
            why = chem_closure_nonfinite;  exit
         endif
         ! A LOWER BOUNDARY THE SWEEP COULD NOT SOLVE is a statement about
         ! the state handed back, like the two above: the ghost cells the
         ! boundary is built on are an iterate, not the solution of the
         ! base handoff and of the ghost's own ionization balance.
         if (ledger%n_ghost_open .gt. 0) then
            why = chem_closure_boundary_open;  exit
         endif
         ! n_offsimplex IS A PROPERTY OF THE ROOT SEARCH AND NOT OF THE
         ! STATE HANDED BACK, so it is reported and does not refuse the
         ! closure.  It counts the cells at which no starting point stayed
         ! inside the element simplex, and BOTH branches that raise it end
         ! in an admissible composition: the molecular branch projects the
         ! closest root onto the element budget and rechecks its reaction
         ! residual, the atomic branch hands back the uncoupled ionization
         ! balance, "admissible by construction", and rechecks the same
         ! residual (ionization_equilibrium.f90, the two clamp sites).  A
         ! cell whose root sits ON A FACE of the simplex raises it by
         ! construction, and a fully dissociated H2 in a 6000 K wind cell
         ! is exactly such a root -- so reading the count as a closure
         ! failure refuses the carrier operator in the regime it exists
         ! for.  MEASURED on LHS 1140 b
         ! molecular_scalar_gj1132_kzz1e9/HeH2.13: 12 cells, all
         ! of them molecular clamps, worst element-budget excursion
         ! 6.5e-04, refused every trial of relax_photochemical_composition
         ! down to the shortest admissible one on every outer pass, so the
         ! H2 carrier never took a single transport step and the state the
         ! route certified was its seed.
         !
         ! WHAT STILL REFUSES: a cell whose accepted state broke a balance
         ! row out of the reals (n_nonfinite, above), a composition that is
         ! not a number, a thermal state that is not a gas (below), and a
         ! temperature that has not stopped moving.  Those are statements
         ! about the state this routine hands back, which is what the
         ! closure is asked about.
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                                             p, T, ntot, ne)
         if (.not. thermal_state_admissible(p, T, heat, cool, eta)) then
            why = chem_closure_not_admissible;  exit
         endif
         dT = maxval(abs(T(1:N) - T_prev(1:N))/max(T(1:N), 1.0d-300))
         ! The increment of every cycle, so that a closure that spends its
         ! budget can be read as contracting or as stalled: one number per
         ! cycle says which, and the exhausted verdict alone does not
         ! (EXHALE_CARRIER_DEBUG=1, default off).
         if (carrier_debug_on()) then
            j_dT = maxloc(abs(T(1:N) - T_prev(1:N))/max(T(1:N), 1.0d-300),&
                          dim = 1)
            write(*,'(A,I0,A,ES11.4,A,I0,A,ES10.3,A,I0,A,I0)')            &
                 '    (chemistry closure) cycle ', k, ': max |dT|/T ', dT, &
                 ' at cell ', j_dT, ', r ', r(j_dT),                       &
                 ', cells off the element simplex ', ledger%n_offsimplex,  &
                 ', nonfinite ', ledger%n_nonfinite
         endif
         if (dT .lt. chem_cycle_tol) then
            ok = .true.;  why = chem_closure_converged;  exit
         endif
      enddo
      n_chem_closure_cycles = n_chem_closure_cycles + n_cycles
      if (ok) then
         n_chem_closures_reached = n_chem_closures_reached + 1
      else
         n_chem_closures_refused = n_chem_closures_refused + 1
      endif
      if (present(reason))    reason    = why
      if (present(increment)) increment = dT
      end subroutine equilibrate_chemistry_at_fixed_conserved_state

      ! ------------------------------------------------------------- !

      logical function thermal_state_admissible(p, T, heat, cool, eta)   &
                       result(adm)
      ! Finite p, T, heat, cool and eta on the physical cells, with p > 0
      ! and T > 0: the state the closure hands back must be a gas.
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: p, T, heat, cool, eta
      adm = all(p(1:N) .gt. 0.0d0) .and. all(p(1:N) .le. huge(1.0d0))    &
            .and. all(T(1:N) .gt. 0.0d0) .and. all(T(1:N) .le. huge(1.0d0)) &
            .and. all(heat(1:N) .eq. heat(1:N))                           &
            .and. all(abs(heat(1:N)) .le. huge(1.0d0))                    &
            .and. all(cool(1:N) .eq. cool(1:N))                           &
            .and. all(abs(cool(1:N)) .le. huge(1.0d0))                    &
            .and. all(eta(1:N) .eq. eta(1:N))                             &
            .and. all(abs(eta(1:N)) .le. huge(1.0d0))
      end function thermal_state_admissible

      ! ------------------------------------------------------------- !

      function chem_closure_reason_text(why) result(txt)
      integer, intent(in) :: why
      character(len=40) :: txt
      select case (why)
      case (chem_closure_converged);      txt = 'converged'
      case (chem_closure_exhausted);      txt = 'cycle budget spent'
      case (chem_closure_nonfinite);      txt = 'nonfinite composition'
      case (chem_closure_not_admissible); txt = 'p, T, heat, cool or eta not admissible'
      case (chem_closure_boundary_open);  txt = 'lower boundary left open'
      case default;                       txt = 'unrecognized reason'
      end select
      end function chem_closure_reason_text

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
      ! only then is the displacement from the pass entry measured: the
      ! relative change of each cell's particle count n_tot + n_e
      ! (carrier_particle_count_change), or, where the run asks for the
      ! retired measure, the largest carrier change over x_ref below.  A
      ! trial whose returned state lies outside `trust`, whose interval the
      ! operator did not cover, whose chemistry did not close, or which left
      ! a species non-finite, is undone (composition and background) and
      ! retried at half the length down to relax_grow_min.  MEASURED
      ! (hot-Uranus carrier reload, on the retired measure): the shortest
      ! admissible trial still moves 6.5e-4 of the entry H2 maximum, so at a
      ! bound below that nothing is kept.
      !
      ! THE FIXED POINT is one at which a full-length trial no longer moves
      ! the carriers, with the chemistry closed on each side of it: for
      ! EVERY solved carrier, the largest change of that carrier over one
      ! step against xref_carrier, the largest fraction THAT carrier holds
      ! at the pass entry (carrier_displacement_on_own_scale).  Each
      ! carrier is judged on its own scale because the carriers span ten
      ! decades: the He 2^3S level is carried as n(2^3S)/n(He nuclei),
      ! 5e-11 at 1.2 R_p on the LHS 1140 b states, while the helium stages
      ! hold fractions of order one.  Against the largest fraction of all
      ! carriers (the test this replaces) a step that moved the level by
      ! 20 % of itself read 1e-11 and passed relax_tol = 1e-10, so a pass
      ! declared the carrier block closed after ONE transport step with the
      ! level row still at 9e-5 of its own terms (MEASURED, the 2000-cell
      ! He/H 0.65 state, phase6/grid2000/continued, every pass), and the
      ! row then sat at 2e-5 to 2e-3 for 300 passes instead of closing.
      ! A displacement test alone does not close the row either (the step
      ! residual is bounded only by n/dt; the loop below states why), so the
      ! fixed point also asks the steady residual of the rows.
      ! x_ref, the largest fraction of all solved carriers, is still what
      ! the reported drift and the retired fraction bound are written on.
      ! It is taken over the SOLVED carriers and not over H2 alone: a run
      ! with transported ionization stages and no molecular network holds
      ! no H2, its H2 maximum is zero, and a 1e-30 floor standing in for it
      ! made a step change of 1e-40 the test (MEASURED,
      ! backup/regression/iontrans_metals run.log).
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
      ! The particle count (n_tot + n_e) of the pass entry and of a trial,
      ! which is what the movement bound is written on
      ! (carrier_particle_count_change).
      real(dp), dimension(1-Ng:N+Ng) :: nsum0, nsum1
      real(dp), dimension(1-Ng:N+Ng) :: p_try, T_try, ntot_t, ne_t
      logical  :: bound_on_fraction
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fprev, fnow, fentry
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar, nrho, wfac
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp), dimension(1-Ng:N+Ng) :: p_entry, T_entry, mbar_entry
      real(dp), dimension(1-Ng:N+Ng) :: nsum_end, ntot_x, ne_x
      real(dp) :: dp_max, dT_max, dmb_max, dns_max
      integer  :: jp_max, jT_max, jmb_max, jns_max
      real(dp) :: drel
      ! EVERYTHING THE TRIAL WRITES, held aside at the top of each trial
      ! and reinstated when the trial is undone.  The enumeration is the
      ! type's, in one place, so that the composition, the rates the
      ! closure evaluator reads, the caloric maps and the base ghost come
      ! back together or not at all.
      type(thermochemical_state) :: held
      real(dp) :: drj, tdiff, tadv, grow, tscale, dmax, x_ref
      ! The largest fraction each solved carrier holds at the pass entry,
      ! the scale its own displacement is judged against at the fixed point.
      real(dp) :: xref_carrier(n_carrier_max)
      ! The steady residual of the carrier rows at the fixed wind, the
      ! measure of the certification (|res|/row_terms_phys over the cells
      ! where the carrier is present), read when a step no longer moves the
      ! carriers.
      real(dp) :: rsteady_rlx, steady_rlx
      integer  :: jsteady_rlx, icsteady_rlx
      real(dp), dimension(1:N,n_carrier_max) :: res_rlx, terms_rlx
      logical,  dimension(1:N,n_carrier_max) :: absent_rlx
      ! The cell that stands furthest outside the movement bound, and by
      ! how much, when a trial is refused on it.
      integer  :: jbnd, icbnd
      real(dp) :: dbnd
      ! HOW FAR THE PASS GOT, AND IN WHAT TIME (EXHALE_CARRIER_DEBUG=1).
      ! The physical interval the kept steps of this pass cover, cell by
      ! cell, so that the advance of a wind cell can be read against its own
      ! chemical time instead of against a step count.
      real(dp), dimension(1-Ng:N+Ng) :: dt_kept
      integer  :: jp, kp, jprobe(3)
      real(dp), parameter :: r_probe(3) = [1.20d0, 1.36d0, 1.60d0]
      integer  :: j, k, ic, ending, step_status, refusal, n_cycles
      logical  :: refused, chem_ok
      ! The diagnostic interval mask's safety stop (default unreachable:
      ! only a waived omitted set can set it).
      logical  :: l22_safety_tripped
      type(element_census_state) :: cen_relax
      ! Why the last trial was undone; the ending names it when the
      ! shortest admissible trial is refused for that reason.
      integer, parameter :: refused_interval  = 1
      integer, parameter :: refused_bound     = 2
      integer, parameter :: refused_chemistry = 3

      drift  = 0.0d0
      nstep  = 0
      rho    = u(1,:)
      if (.not. transported_rows_exist() .or. .not. bg_ready) then
         if (present(outcome)) outcome = carrier_relax_nothing_to_advance
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                                          p, T, ntot_e, ne_e)
         return
      endif
      call element_census_take('relax_photochemical_composition',         &
                               rho, f_sp, cen_relax)
      call l22_interval_mask_setup()
      l22_safety_tripped = .false.
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
      x_ref  = 1.0d-30
      xref_carrier = 1.0d-30
      do ic = 1, n_carrier
         if (carrier_solved(ic)) then
            x_ref = max(x_ref, maxval(fentry(1:N,ic)))
            xref_carrier(ic) = max(1.0d-30, maxval(fentry(1:N,ic)))
         endif
      enddo
      ! The primitive state of the entry composition at u: what is handed
      ! back when no step is kept.
      call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,     &
                                                          p, T, ntot_e, ne_e)
      ! The entry particle count, the quantity the movement bound compares a
      ! trial against, and which measure this run uses.
      nsum0             = ntot_e + ne_e
      ! The primitive state the wind reads at the entry
      ! composition, kept so that what the pass hands back can be reported
      ! against it.  mbar = rho n0 mu / n_tot is the mean mass per particle,
      ! the quantity of the composition the momentum and energy rows see.
      p_entry    = p
      T_entry    = T
      mbar_entry = mbar
      bound_on_fraction = carrier_bound_on_fraction()
      ! An empty selection set is not a comparison: no cell is left to
      ! choose the interval, so the experiment has no reading and the
      ! relaxation says so instead of advancing on an undefined multiplier.
      if (l22_mask_mode .ne. l22_mask_off .and. l22_mask_n_select .eq. 0) &
      then
         write(*,'(A)') '    (L22 interval mask) the selection set is'//  &
              ' empty: no cell is left to choose the interval'
         if (present(outcome)) outcome = carrier_relax_mask_safety_stop
         call element_census_verify(cen_relax, rho, f_sp,                 &
                                    rho_is_fixed = .true.)
         return
      endif
      grow    = 1.0d0
      ending  = carrier_relax_step_budget
      fnow    = fentry
      dt_kept = 0.0d0
      carrier_rows_advect = .true.
      do k = 1, relax_maxstep
         call save_thermochemical_state(held, f_sp, p, T, heat, cool, eta)
         refused   = .false.
         refusal   = 0
         ! The bound measure of THIS trial: zero until the bound block runs,
         ! so a trial refused before it reports no bound measure.
         dbnd      = 0.0d0
         jbnd      = 0
         l22_last_d_select = 0.0d0;  l22_last_j_select = 0
         l22_last_d_omit   = 0.0d0;  l22_last_j_omit   = 0
         call photochemical_transport_step(rho, v, f_sp, dt_code*grow,    &
                                           step_status, trial = .true.)
         if (step_status .ne. carrier_interval_covered) then
            refused = .true.;  refusal = refused_interval
         endif
         if (.not. refused) then
            call equilibrate_chemistry_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, heat, cool, eta,       &
                                             chem_ok, n_cycles,           &
                                             n_chem_last_reason,          &
                                             chem_last_increment)
            if (.not. chem_ok) then
               refused = .true.;  refusal = refused_chemistry
            endif
         endif
         if (.not. refused) then
            call carrier_state(rho, f_sp, fnow, ntot, nrho, wfac,        &
                               TK, mbar, nH_free, nO_free, nC_free)
            ! THE MOVEMENT BOUND, on the particle count of each cell --
            ! what the fixed wind responds to -- unless the run asks for the
            ! retired carrier-fraction measure
            ! (carrier_particle_count_change states both and why).
            if (bound_on_fraction) then
               call carrier_worst_composition_change(fnow, fentry, jbnd,  &
                                                     icbnd, dbnd)
               dbnd = dbnd/x_ref
            else
               call pressure_and_temperature_at_fixed_conserved_state(u,  &
                          f_sp, p_try, T_try, ntot_t, ne_t)
               nsum1 = ntot_t + ne_t
               if (l22_mask_mode .eq. l22_mask_off) then
                  call carrier_particle_count_change(nsum1, nsum0, jbnd,  &
                                                     dbnd)
               else
                  ! The diagnostic interval mask: the selection set chooses
                  ! the interval, the omitted set only vetoes it in veto
                  ! mode.  Both measures are kept for the log.
                  call l22_particle_count_change_by_set(nsum1, nsum0,     &
                       l22_last_j_select, l22_last_d_select,             &
                       l22_last_j_omit,   l22_last_d_omit)
                  if (l22_mask_mode .eq. l22_mask_veto .and.             &
                      l22_last_d_omit .gt. l22_last_d_select) then
                     jbnd = l22_last_j_omit;    dbnd = l22_last_d_omit
                  else
                     jbnd = l22_last_j_select;  dbnd = l22_last_d_select
                  endif
                  ! The safety stop of the waived set: a factor two in the
                  ! particle count of one cell in one trial is outside any
                  ! regime this experiment can be read in.
                  if (l22_mask_mode .eq. l22_mask_waive .and.            &
                      l22_last_d_omit .gt. 1.0d0) l22_safety_tripped =   &
                      .true.
               endif
               icbnd = ic_H2
            endif
            if (l22_safety_tripped) then
               refused = .true.;  refusal = refused_bound
               write(*,'(A,I0,A,ES10.3)') '    (L22 interval mask)'//     &
                    ' SAFETY STOP: omitted cell ', l22_last_j_omit,       &
                    ' moved its particle count by a relative ',           &
                    l22_last_d_omit
            else if (dbnd .gt. trust) then
               refused = .true.;  refusal = refused_bound
               bound_last_j        = jbnd
               bound_last_ic       = icbnd
               bound_last_dabs     = dbnd
               bound_last_fraction = bound_on_fraction
               if (jbnd .ge. 1) then
                  if (bound_on_fraction) then
                     bound_last_entry = fentry(jbnd,icbnd)
                  else
                     bound_last_entry = nsum0(jbnd)
                  endif
               endif
            endif
         endif
         if (refused) then
            ! THE TRIAL AND EVERYTHING IT WROTE ARE UNDONE.  Not the
            ! composition and the frozen cell state alone: the rate state
            ! the closure evaluator measures a composition against, the
            ! caloric maps, the base ghost and the carrier subsystem are
            ! put back with them, so that what the next trial and the
            ! caller read is one state and not the composition of the
            ! kept steps beside the rates of a discarded trial.
            call restore_thermochemical_state(held, rho, f_sp, p, T,      &
                                              heat, cool, eta)
            fnow = fprev
            ! EVERY REJECTED TRIAL, WITH ITS REASON AND THE CELL THAT
            ! CARRIED IT (EXHALE_CARRIER_DEBUG=1, default off): the mask
            ! experiment is read off the sequence of trials of one pass, not
            ! off the last one alone.
            if (carrier_debug_on()) then
               write(*,'(A,I0,A,ES10.3,A,A,A,ES10.3,A,I0)')               &
                    '      trial ', k, ': grow ', grow, ', REFUSED on ',  &
                    trim(carrier_trial_refusal_text(refusal)),            &
                    ', measure ', dbnd, ' at cell ', jbnd
               if (l22_mask_mode .ne. l22_mask_off)                       &
                  write(*,'(A,ES10.3,A,I0,A,ES10.3,A,I0)')                &
                       '        selection set ', l22_last_d_select,       &
                       ' at cell ', l22_last_j_select,                    &
                       ', omitted set ', l22_last_d_omit, ' at cell ',    &
                       l22_last_j_omit
            endif
            if (l22_safety_tripped) then
               ending = carrier_relax_mask_safety_stop
               exit
            endif
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
         nstep   = nstep + 1
         dt_kept = dt_kept + dt_code*grow
         ! EVERY ACCEPTED TRIAL (EXHALE_CARRIER_DEBUG=1, default off): the
         ! multiplier it was accepted at, the bound measure it stood at, and
         ! the physical interval the shortest cell of the column then
         ! covered.
         if (carrier_debug_on()) then
            write(*,'(A,I0,A,ES10.3,A,ES10.3,A,I0,A,ES10.3,A)')           &
                 '      trial ', k, ': grow ', grow, ' ACCEPTED, '//      &
                 'measure ', dbnd, ' at cell ', jbnd,                     &
                 ', shortest cell interval ',                             &
                 minval(dt_code(1:N))*grow*tscale, ' s'
            if (l22_mask_mode .ne. l22_mask_off)                          &
               write(*,'(A,ES10.3,A,I0,A,ES10.3,A,I0)')                   &
                    '        selection set ', l22_last_d_select,          &
                    ' at cell ', l22_last_j_select,                       &
                    ', omitted set ', l22_last_d_omit, ' at cell ',       &
                    l22_last_j_omit
         endif
         dmax  = carrier_displacement_on_own_scale(fnow, fprev,          &
                                                   xref_carrier)
         ! A STEP THAT NO LONGER MOVES THE CARRIERS IS NOT YET A CLOSED
         ! BLOCK.  An implicit step of length dt leaves the steady residual
         ! r = (n_new - n_old)/dt, so a displacement below relax_tol bounds
         ! r only by relax_tol n/dt, and in a row whose own terms are far
         ! below n/dt -- a level in near-local balance whose transport
         ! flux changes sign, as He 2^3S at 1.2-1.25 R_p on the LHS 1140 b
         ! states, row terms 2e-6 cm^-3 s^-1 against n/dt ~ 3e-2 -- that
         ! bound says nothing: the pass ended after one step with the row at
         ! 5e-5 to 9e-5 of its terms (MEASURED, phase6/grid2000).  The
         ! fixed point is therefore also asked of the residual itself: the
         ! steady carrier residual on the certification's measure, over the
         ! cells where the carrier is present, at or below
         ! carrier_accept_tol, the residual the step's own Newton accepts.
         ! Until it is, the steps go on growing (grow below), which is the
         ! pseudo-transient continuation to the steady block.
         if (dmax .lt. relax_tol .and. grow .ge. 1.0d0) then
            call carrier_steady_residual(rho, v, f_sp, rsteady_rlx,       &
                                         jsteady_rlx, icsteady_rlx,       &
                                         res_out = res_rlx,               &
                                         terms_out = terms_rlx,           &
                                         absent_out = absent_rlx)
            steady_rlx = 0.0d0
            do ic = 1, n_carrier
               if (.not. carrier_solved(ic)) cycle
               do j = 1, N
                  if (absent_rlx(j,ic)) cycle
                  steady_rlx = max(steady_rlx, abs(res_rlx(j,ic))         &
                                   /max(terms_rlx(j,ic), 1.0d-300))
               enddo
            enddo
            if (carrier_debug_on())                                       &
               write(*,'(A,I0,A,ES10.3)') '      trial ', k,              &
                    ': steady carrier residual of the present rows ',     &
                    steady_rlx
            if (steady_rlx .le. carrier_accept_tol) then
               ending = carrier_relax_fixed_point
               exit
            endif
         endif
         fprev = fnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo
      carrier_rows_advect = .false.
      call carrier_state(rho, f_sp, fnow, ntot, nrho, wfac, TK,          &
                         mbar, nH_free, nO_free, nC_free)
      ! The diffusivities of the composition this pass returns, at the
      ! background of its last kept chemistry (or of the entry, when no
      ! trial was kept): the state the run holds from here.
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      if (.not. allocated(pct_Dco))                                      &
         allocate(pct_Dco(1-Ng:N+Ng,n_carrier_max))
      pct_Dco = Dco
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
      ! WHAT THE PASS DID TO THE CELLS THAT CARRY THE REFUSING ROW
      ! (EXHALE_CARRIER_DEBUG=1, default off).  Three probe radii in the
      ! wind, each reported as the H2 fraction it entered the pass with, the
      ! one it leaves with, their ratio, and the physical interval the kept
      ! steps covered there -- the last being what a reader compares against
      ! the cell's own chemical time to say whether the advance is limited
      ! by the interval or by the bound.
      if (carrier_debug_on()) then
         do kp = 1, 3
            jprobe(kp) = 1
            do j = 1, N
               if (abs(r(j) - r_probe(kp)) .lt.                          &
                   abs(r(jprobe(kp)) - r_probe(kp))) jprobe(kp) = j
            enddo
         enddo
         write(*,'(A,I0,A,ES10.3)') '    (carrier relaxation) kept ',     &
              nstep, ' step(s); drift ', drift
         do kp = 1, 3
            jp = jprobe(kp)
            write(*,'(A,I0,A,F7.4,A,ES12.5,A,ES12.5,A,ES10.3,A,ES10.3,A)')&
                 '      probe cell ', jp, ' r ', r(jp), ': x(H2) ',        &
                 fentry(jp,ic_H2), ' -> ', fnow(jp,ic_H2), ', ratio ',     &
                 fnow(jp,ic_H2)/max(fentry(jp,ic_H2), 1.0d-300),           &
                 ', kept interval ', dt_kept(jp)*tscale, ' s'
         enddo
      endif
      ! WHAT THE RETURNED COMPOSITION DID TO THE STATE THE
      ! WIND READS (EXHALE_L22B_DISPLACEMENT=1, default off).  A carrier
      ! relaxation at a FIXED wind may close on its own rows and still hand
      ! back a composition whose pressure, temperature and mean mass per
      ! particle are far from the ones the hydrodynamic rows were last
      ! solved at; those three are what the momentum and energy rows carry,
      ! so their displacement is what says whether the coupled solve that
      ! follows starts inside a regime it can follow.
      if (l22b_displacement) then
         call pressure_and_temperature_at_fixed_conserved_state(u, f_sp,  &
                                             p, T, ntot_x, ne_x)
         nsum_end = ntot_x + ne_x
         dp_max = 0.0d0;  jp_max = 0
         dT_max = 0.0d0;  jT_max = 0
         dmb_max = 0.0d0; jmb_max = 0
         dns_max = 0.0d0; jns_max = 0
         do j = 1, N
            if (p_entry(j) .gt. 0.0d0) then
               drel = abs(p(j) - p_entry(j))/p_entry(j)
               if (drel .gt. dp_max) then
                  dp_max = drel;  jp_max = j
               endif
            endif
            if (T_entry(j) .gt. 0.0d0) then
               drel = abs(T(j) - T_entry(j))/T_entry(j)
               if (drel .gt. dT_max) then
                  dT_max = drel;  jT_max = j
               endif
            endif
            if (mbar_entry(j) .gt. 0.0d0) then
               drel = abs(mbar(j) - mbar_entry(j))/mbar_entry(j)
               if (drel .gt. dmb_max) then
                  dmb_max = drel;  jmb_max = j
               endif
            endif
            if (nsum0(j) .gt. 0.0d0) then
               drel = abs(nsum_end(j) - nsum0(j))/nsum0(j)
               if (drel .gt. dns_max) then
                  dns_max = drel;  jns_max = j
               endif
            endif
         enddo
         write(*,'(A)') '    (L22b) displacement of the primitive state'//&
              ' between the pass entry and the composition returned:'
         write(*,'(A,ES10.3,A,I0,A,F8.4,A,ES10.3,A,I0,A,F8.4)')           &
              '      |dp|/p  ', dp_max, ' at cell ', jp_max, ' r ',       &
              r(max(jp_max,1)), ';  |dT|/T  ', dT_max, ' at cell ',       &
              jT_max, ' r ', r(max(jT_max,1))
         write(*,'(A,ES10.3,A,I0,A,F8.4,A,ES10.3,A,I0,A,F8.4)')           &
              '      |dmbar|/mbar  ', dmb_max, ' at cell ', jmb_max,      &
              ' r ', r(max(jmb_max,1)), ';  |d(n_tot+n_e)|/(n_tot+n_e) ', &
              dns_max, ' at cell ', jns_max, ' r ', r(max(jns_max,1))
      endif
      call element_census_verify(cen_relax, rho, f_sp, rho_is_fixed=.true.)
      if (present(outcome)) outcome = ending
      end subroutine relax_photochemical_composition

      ! End of module

      end module diffusive_photochemistry

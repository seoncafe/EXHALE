      program steady_species_rows_tests
      ! THE UNKNOWN SPACE OF THE STATIONARY SYSTEM, ONE ROW PER TRANSPORTED
      ! BALANCE (docs/PLAN_20260906_rev2.md B5 target 1;
      ! docs/b1_target_system_20260906.md T1.7).
      !
      ! Every transported balance the configuration activates is a row of the
      ! stationary system and its transported quantity an unknown of the
      ! cell. The registry of steady_newton is the only place that mapping is
      ! written, and these rows assert it against the transported set the
      ! carrier operator fixed (carrier_set_init), configuration by
      ! configuration.
      !
      ! WHAT EACH ROW IS FOR.
      !   * the row count per configuration: the balances of T1.7's table,
      !     and no others. A registry that missed OH, H2O and CO left three
      !     of the oxygen cycle's four balances outside the unknown space
      !     while the solve reported convergence on the fourth.
      !   * the band geometry: |row - col| <= 3*nvar - 1 from the WENO3
      !     reach of two cells plus the nvar - 1 variables of the cell
      !     itself, and a coloring stride of kl + ku + 1. A wrong stride
      !     gives same-color columns overlapping row supports, and the
      !     colored finite-difference Jacobian is then contaminated rather
      !     than exact. The two geometries the code carried as constants,
      !     nvar = 3 -> (8, 17) and nvar = 4 -> (11, 23), are rows here.
      !   * the species each row names: read back through
      !     species_row_carrier_index, so a row cannot carry one balance and
      !     write another species.
      !
      ! Every row prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line. Exit status is nonzero if any row fails.

      use global_parameters
      use species_table
      use element_inventory
      use element_census, only: element_ratio_gate
      use diffusive_photochemistry, only: carrier_set_init, n_carrier,   &
                                          carrier_solved, ic_H2, ic_OH,  &
                                          ic_H2O, ic_CO, ic_Hp,          &
                                          carrier_species_index,         &
                                          carrier_element_headroom,      &
                                          carrier_headroom_known,        &
                                          carrier_headroom_set_for_test, &
                                          carrier_advective_divergence,  &
                                          carrier_mass_amu, n_carrier_max
      use species_advective_transport, only: species_advective_update
      use composition, only: mass_per_H_nucleus_without_He
      use test_columns, only: column_carrying_its_own_density,            &
                              column_mass_closure
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use binary_element_diffusion, only: element_diffusion_step,         &
                                          element_transport_residual,     &
                                          element_mass_fractions,         &
                                          relax_element_composition
      use certification, only: cert_tol_element_at, cert_tol_carrier_at, &
                               cert_tol_reported_only,                   &
                               cert_regime_wind_r, cert_regime_layer_r
      use steady_newton, only: set_transported_species_rows,             &
                               pack_species_rows,                        &
                               write_species_rows_into_composition,      &
                               neq_newton, kl_jac, ku_jac, ncolor_jac,   &
                               n_species_rows, species_row_carrier_index,&
                               carrier_unknown_on, cell_row_scales,      &
                               distance_from_certification,              &
                               hydrodynamic_distance_from_certification_by_cell, &
                               increment_of_what_the_residual_reads,     &
                               species_unknowns_outside_their_bounds,    &
                               freeze_species_unknown_box,               &
                               read_species_unknown_space_controls,      &
                               nvar_jac,                                 &
                               eq_channel_particle_count,                &
                               eq_channel_h2_caloric,                    &
                               eq_channel_radiative,                     &
                               species_row_base_equation,                &
                               base_row_reservoir_condition,             &
                               base_row_control_volume_balance,          &
                               freeze_element_constraint_rows,           &
                               element_constraint_rows_known,            &
                               element_constraint_of_the_cell,           &
                               element_constraint_budget_of_the_cell,    &
                               element_constraint_gradient_of_the_cell,  &
                    project_out_of_the_active_element_constraints,       &
                    largest_step_inside_the_element_constraints,         &
                               largest_step_inside_the_species_box,      &
                               zero_the_blocked_components_of,           &
                               fix_active_species_bounds,                &
                               release_active_species_bounds,            &
                               n_active_element_rows,                    &
                               restore_the_element_budget,               &
                               n_element_constraint_rows

      implicit none
      integer :: nfail
      nfail = 0

      N  = 7

      ! ---- configuration 3 of T1.7: molecular, no transported balance ----
      thereis_mol         = .true.
      carrier_transport   = .false.
      thereis_oxychem     = .false.
      ionization_transport = .false.
      call carrier_set_init()
      call set_transported_species_rows(.false.)
      call int_row('rows_no_transport', n_species_rows(), 0, nfail)
      call int_row('nvar_no_transport', neq_newton()/N,   3, nfail)
      call int_row('band_no_transport', kl_jac,           8, nfail)
      call int_row('colors_no_transport', ncolor_jac,    17, nfail)
      call log_row('carrier_unknown_off', carrier_unknown_on(), .false., &
                   nfail)

      ! ---- configuration 4: the H2 balance alone ----
      carrier_transport = .true.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_carriers_H2',   n_species_rows(), 1, nfail)
      call int_row('nvar_carriers_H2',   neq_newton()/N,   4, nfail)
      call int_row('band_carriers_H2',   kl_jac,          11, nfail)
      call int_row('colors_carriers_H2', ncolor_jac,      23, nfail)
      call int_row('row1_is_H2', species_row_carrier_index(1), ic_H2,    &
                   nfail)
      call log_row('carrier_unknown_on', carrier_unknown_on(), .true.,   &
                   nfail)
      ! A CARRIER ROW IS A CONTROL-VOLUME BALANCE AT THE BASE, declared and
      ! not inferred (item N8a, deliverable 3): the inner ghost is the
      ! inflow reservoir and the base face flux it carries is a constant of
      ! the row, so cell 1 balances transport against chemistry there like
      ! every other cell and is nondimensionalized by its own terms.
      call int_row('carrier_base_equation_is_a_control_volume_balance',   &
                   species_row_base_equation(1),                          &
                   base_row_control_volume_balance, nfail)

      ! ---- configuration 5: the four balances of the oxygen cycle ----
      thereis_oxychem = .true.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_oxygen',   n_species_rows(), 4, nfail)
      call int_row('nvar_oxygen',   neq_newton()/N,   7, nfail)
      call int_row('band_oxygen',   kl_jac,          20, nfail)
      call int_row('colors_oxygen', ncolor_jac,      41, nfail)
      call int_row('oxygen_row2_is_OH',  species_row_carrier_index(2),   &
                   ic_OH,  nfail)
      call int_row('oxygen_row3_is_H2O', species_row_carrier_index(3),   &
                   ic_H2O, nfail)
      call int_row('oxygen_row4_is_CO',  species_row_carrier_index(4),   &
                   ic_CO,  nfail)

      ! ---- configuration 8: the H2 and H+ balances ----
      thereis_oxychem      = .false.
      ionization_transport = .true.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_hp',   n_species_rows(), 2, nfail)
      call int_row('nvar_hp',   neq_newton()/N,   5, nfail)
      call int_row('band_hp',   kl_jac,          14, nfail)
      call int_row('colors_hp', ncolor_jac,      29, nfail)
      call int_row('hp_row2_is_Hp', species_row_carrier_index(2), ic_Hp, &
                   nfail)
      call int_row('band_symmetric', ku_jac, kl_jac, nfail)

      ! ---- the whole oxygen and proton set together ----
      thereis_oxychem = .true.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_oxygen_and_hp',   n_species_rows(), 5, nfail)
      call int_row('nvar_oxygen_and_hp',   neq_newton()/N,   8, nfail)
      call int_row('band_oxygen_and_hp',   kl_jac,          23, nfail)
      call int_row('colors_oxygen_and_hp', ncolor_jac,      47, nfail)

      ! ---- configuration 6 of T1.7: the He/H partition alone ----
      ! The element rows come FIRST in the registry, because their
      ! write-back is a projection of the whole species vector and a carrier
      ! written before it would be rescaled by it.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      carrier_transport    = .false.
      thereis_mol          = .false.
      thereis_He           = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_he_diffusion',   n_species_rows(), 1, nfail)
      call int_row('nvar_he_diffusion',   neq_newton()/N,   4, nfail)
      call int_row('band_he_diffusion',   kl_jac,          11, nfail)
      call int_row('colors_he_diffusion', ncolor_jac,      23, nfail)
      call int_row('he_row_is_not_a_carrier',                            &
                   species_row_carrier_index(1), 0, nfail)

      ! ---- configuration 7: the He/H partition and every trace element ----
      he_metal_diffusion = .true.
      thereis_metals     = .true.
      call set_transported_species_rows(.true.)
      call int_row('rows_element_profile', n_species_rows(),             &
                   1 + n_melem, nfail)
      call int_row('nvar_element_profile', neq_newton()/N,               &
                   4 + n_melem, nfail)
      call int_row('band_element_profile', kl_jac,                       &
                   3*(4 + n_melem) - 1, nfail)
      call int_row('colors_element_profile', ncolor_jac,                 &
                   6*(4 + n_melem) - 1, nfail)
      call int_row('trace_rows_are_not_carriers',                        &
                   species_row_carrier_index(2), 0, nfail)

      ! ---- elements and carriers together: the elements are registered
      !      first, so the carrier rows follow them ----
      thereis_mol        = .true.
      carrier_transport  = .true.
      he_metal_diffusion = .false.
      thereis_metals     = .false.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      call int_row('rows_he_and_H2', n_species_rows(), 2, nfail)
      call int_row('he_row_first', species_row_carrier_index(1), 0, nfail)
      call int_row('carrier_row_after_the_element',                      &
                   species_row_carrier_index(2), ic_H2, nfail)

      ! ---- the registry is emptied again, and the geometry with it ----
      thereis_He         = .false.
      he_diffusion       = .false.
      call set_transported_species_rows(.false.)
      call int_row('rows_reset',   n_species_rows(), 0, nfail)
      call int_row('nvar_reset',   neq_newton()/N,   3, nfail)
      call int_row('colors_reset', ncolor_jac,      17, nfail)

      ! ---- the f_sp column each carrier names ----
      call int_row('species_of_H2',  carrier_species_index(ic_H2),  isp_H2,  nfail)
      call int_row('species_of_OH',  carrier_species_index(ic_OH),  isp_OH,  nfail)
      call int_row('species_of_H2O', carrier_species_index(ic_H2O), isp_H2O, nfail)
      call int_row('species_of_CO',  carrier_species_index(ic_CO),  isp_CO,  nfail)
      call int_row('species_of_Hp',  carrier_species_index(ic_Hp),  isp_HII, nfail)

      ! ---- the element headroom of each carrier is its own stoichiometry --
      ! H2 holds two hydrogen nuclei and H+ one, OH and H2O one oxygen each,
      ! CO one of oxygen and one of carbon so the smaller free density bounds
      ! it. Charging a carrier to the wrong element is not a small
      ! mis-scaling: in a hydrogen and helium atmosphere the oxygen free
      ! density is zero.
      call real_row('headroom_H2_is_half_H',                             &
                    carrier_element_headroom(ic_H2, 8.0d0, 3.0d0, 5.0d0), &
                    4.0d0, nfail)
      call real_row('headroom_Hp_is_H',                                  &
                    carrier_element_headroom(ic_Hp, 8.0d0, 3.0d0, 5.0d0), &
                    8.0d0, nfail)
      call real_row('headroom_OH_is_O',                                  &
                    carrier_element_headroom(ic_OH, 8.0d0, 3.0d0, 5.0d0), &
                    3.0d0, nfail)
      call real_row('headroom_H2O_is_O',                                 &
                    carrier_element_headroom(ic_H2O, 8.0d0, 3.0d0, 5.0d0),&
                    3.0d0, nfail)
      call real_row('headroom_CO_is_min_O_C',                            &
                    carrier_element_headroom(ic_CO, 8.0d0, 5.0d0, 3.0d0), &
                    3.0d0, nfail)

      call element_row_scale_is_the_certification_scale(nfail)
      call judged_distance_is_the_worst_row_over_its_tolerance(nfail)
      call advective_term_is_the_stage_operator(nfail)
      call elemental_relaxation_fixed_point_is_the_row(nfail, .false.)
      call elemental_relaxation_fixed_point_is_the_row(nfail, .true.)
      call relaxation_drift_covers_every_element(nfail)
      call eq_measure_reads_what_the_residual_reads(nfail)
      call species_box_is_a_box_and_not_a_ball(nfail)
      ! AFTER the row above, which states the box of a solve whose carrier
      ! operator has not assembled and therefore has no element budget: the
      ! subroutine below writes one, and it stays written.
      call element_budget_is_a_row_of_the_step_and_not_a_face(nfail)
      call a_direction_along_a_shared_constraint_is_not_discarded(nfail)
      call an_infeasible_start_is_restored_by_measured_violation(nfail)
      call element_constraint_derivatives_against_differences(nfail)
      ! THE ELEMENTAL CONSTRAINT THE BOX SIDES DO NOT STATE (item N4a).
      call map_and_operator_totals_are_one_number(nfail)
      call shared_element_constraint_is_the_feasible_set(nfail)
      call reference_states_are_classified_under_the_map(nfail)
      ! THE UPPER GHOST OF A TRANSPORTED ROW (item N26).
      call the_upper_ghost_of_a_transported_element_is_the_iterate(nfail)
      call the_upper_ghost_of_a_transported_carrier_is_the_iterate(nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') '  ', nfail, ' row(s) failed'
         error stop 1
      endif

      contains

      ! ================================================================= !

      subroutine element_row_scale_is_the_certification_scale(nf)
      ! AN ELEMENT ROW IS READ ON THE SCALE THE CERTIFICATION READS IT ON
      ! (decision 20 a, item N20).
      !
      ! cell_row_scales divides every row of the species-row route by
      ! merit_row_scale_from_certification of that row's certification
      ! scale, so the merit the step control descends is the 2-norm of the
      ! quantities the gate takes a maximum of. For an element row that
      ! scale is the operator's own escale_he or escale_tr written per code
      ! time, with no floor of the solver's own on top of it: a floor that
      ! the certification does not apply would make the two functionals
      ! different again, which is the defect this item removes. Until N20
      ! the row was floored at the flow-time rate of the quantity it
      ! transports, and the eight element rows then held 0.437 of the
      ! squared merit while holding nothing of ||R|| (N7 report, item 7).
      !
      ! WHAT THE FLOOR WAS FOR IS NOW A MEASUREMENT AND NOT AN ASSUMPTION.
      ! The row scale divides the row's DERIVATIVES as well as its residual,
      ! so a scale far below the size on which the row responds turns an
      ! ordinary derivative into an enormous entry of the banded model:
      ! MEASURED on the hot Uranus reload with the operator's own floor
      ! (1e-20 of the flow-time rate), the helium row scales spanned
      ! 1.277e-29 to 2.360 over the column, a dF/dY of 46.7 became a scaled
      ! band entry of 4.890e+22, and the damped Gauss-Newton escape reported
      ! ||grad merit|| = 1.104e+44 and no descent at outer iteration 8. What
      ! the conditioning of the band does under the certification scaling is
      ! measured on both fixtures by dgbcon and reported with the item, not
      ! guarded by a floor the gate does not share.
      !
      ! The rows below are EXACT: with the code's own scales set to unity
      ! the operator's term scale is 1, so every cell that carries a
      ! transport balance reads 1 whatever its rho X. The base cell is the
      ! Dirichlet reservoir and takes the unknown's own scale, because its
      ! equation is X(1) = X_reservoir and not a transport balance and the
      ! certification measures no balance there.
      !
      ! u IS NOT PASSED HERE, which is the branch a caller that holds no
      ! assembled state takes: residual_row_scale is the largest term the
      ! row itself holds and exists only for a state, so the three
      ! hydrodynamic slots keep the column scales.
      integer, intent(inout) :: nf
      real*8, allocatable :: Yv(:), Dv(:), Drv(:)
      integer :: nv, j, ihe
      ! The unit conversion of the helium row is R0/(v0 n0 mu); mu is the
      ! hydrogen atom mass and a parameter, so n0 carries its reciprocal and
      ! the conversion is one.  The term scale the registry initializes
      ! (escale_he = 1) then reaches cell_row_scales as 1.
      R0 = 1.0d0;  v0 = 1.0d0;  n0 = 1.0d0/mu
      thereis_oxychem      = .false.
      ionization_transport = .false.
      carrier_transport    = .false.
      thereis_mol          = .false.
      thereis_He           = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      nv = neq_newton()/N
      call int_row('certification_scale_nvar', nv, 4, nf)
      allocate(Yv(nv*N), Dv(nv*N), Drv(nv*N))
      Yv = 1.0d0;  Dv = 1.0d0
      do j = 1, N
         ihe = nv*(j-1) + 4
         ! rho X = 4 in cells 1..3 and 0.25 from cell 4 on: under the old
         ! rule the first group took the flow-time floor and the second kept
         ! the term scale, and the two answers differed by a factor 4.
         if (j .le. 3) then
            Dv(nv*(j-1)+1) = 2.0d0;  Dv(ihe) = 2.0d0
         else
            Dv(nv*(j-1)+1) = 0.5d0;  Dv(ihe) = 0.5d0
         endif
      enddo
      call cell_row_scales(Yv, Dv, Drv)
      ! WHICH EQUATION THE ROW CARRIES AT THE BASE IS DECLARED, and the
      ! scale below follows the declaration rather than restating it (item
      ! N8a, deliverable 3): an element row is the Dirichlet reservoir
      ! condition there, because element_transport_residual poses no
      ! transport balance in cell 1.
      call int_row('element_base_equation_is_the_reservoir_condition',    &
                   species_row_base_equation(1),                          &
                   base_row_reservoir_condition, nf)
      ! The base cell is the reservoir statement: the unknown's own scale.
      call real_row('element_base_row_scale_is_the_unknown',              &
                    Drv(4), Dv(4), nf)
      ! Every other cell is the certification scale, whatever its rho X.
      call real_row('element_row_scale_is_the_operator_term_scale',       &
                    Drv(nv*1 + 4), 1.0d0, nf)
      call real_row('element_row_scale_takes_no_floor_of_its_own',        &
                    Drv(nv*3 + 4), 1.0d0, nf)
      ! With no state to read a row scale from, the hydrodynamic rows are
      ! the column scales.
      call real_row('hydro_row_scale_without_a_state_is_the_column_scale',&
                    Drv(nv*3 + 1), Dv(nv*3 + 1), nf)
      deallocate(Yv, Dv, Drv)
      end subroutine element_row_scale_is_the_certification_scale

      ! ================================================================= !

      subroutine judged_distance_is_the_worst_row_over_its_tolerance(nf)
      ! THE ACCEPTANCE BOUNDS THE QUANTITY THE STATE IS JUDGED BY (item B5e).
      !
      ! The merit the step control descends on is || F/Drow ||_2, a 2-norm on
      ! the Newton's own column scales; a state is judged by a MAXIMUM over
      ! cells of each row against its own largest term, one row at a time and
      ! each against its own tolerance. distance_from_certification is what
      ! makes those rows one comparable quantity: each divided by the
      ! tolerance it is judged against, and the largest taken, so that d < 1
      ! is exactly the condition the certification states.
      !
      ! The rows are EXACT equalities. What a wrong expression gives is a
      ! distance that follows one row and ignores another, which is the
      ! defect the whole item is about, so there is nothing here to put a
      ! tolerance on. The ratios are powers of two so that the division is
      ! exact in binary and the equality can be asked for as one.
      !
      ! THE TWO SPECIES ROWS ARE MEASURED IN THE WIND, where their
      ! tolerance is a statement (decision 22 (a)): the rows passed here
      ! are the gated part of their column and the divisor is the value
      ! the accessor takes there, 1e-5, so the ratios below are that
      ! tolerance times a power of two.
      !
      ! THE HYDRODYNAMIC SLOT ARRIVES AS A DISTANCE: each row of each cell
      ! over the certification tolerance that row carries there (momentum
      ! 1e-8 and energy 1e-6 for the column, the continuity row 3e-12 or
      ! ten times the cell's own rounding floor), the largest taken
      ! (hydrodynamic_distance_from_certification_by_cell), so that the
      ! ledger ranks iterates by the row that the acceptance refuses and
      ! not by the run's "Resid tol", which is the solver's target and not
      ! a certification tolerance.
      integer, intent(inout) :: nf
      real*8 :: rows(3), tel, tca
      real*8 :: q(3,1), fl(1)
      tel = cert_tol_element_at(cert_regime_wind_r)
      tca = cert_tol_carrier_at(cert_regime_wind_r)
      call real_row('the_element_tolerance_the_distance_divides_by',      &
                    tel, 1.0d-5, nf)
      call real_row('the_carrier_tolerance_the_distance_divides_by',      &
                    tca, 1.0d-5, nf)
      call log_row('and_a_layer_cell_is_not_gated_by_either',             &
                   (cert_tol_element_at(cert_regime_layer_r) .ge.         &
                    cert_tol_reported_only) .and.                         &
                   (cert_tol_carrier_at(cert_regime_layer_r) .ge.         &
                    cert_tol_reported_only), .true., nf)
      ! Each row alone, in the units of its own tolerance.
      rows = (/ 2.0d0, 0.0d0, 0.0d0 /)
      call real_row('judged_distance_reads_the_hydrodynamic_rows',        &
                    distance_from_certification(rows), 2.0d0, nf)
      rows = (/ 0.0d0, 4.0d0*tel, 0.0d0 /)
      call real_row('judged_distance_reads_the_element_row',              &
                    distance_from_certification(rows), 4.0d0, nf)
      rows = (/ 0.0d0, 0.0d0, 8.0d0*tca /)
      call real_row('judged_distance_reads_the_carrier_row',              &
                    distance_from_certification(rows), 8.0d0, nf)
      ! And together it is the WORST of them, not a sum and not an average:
      ! a state is certified when every row is within its own tolerance.
      rows = (/ 0.5d0, 4.0d0*tel, 2.0d0*tca /)
      call real_row('judged_distance_is_the_worst_row',                   &
                    distance_from_certification(rows), 4.0d0, nf)
      ! d < 1 is the certification condition itself.
      rows = (/ 0.5d0, 0.5d0*tel, 0.5d0*tca /)
      call log_row('judged_distance_below_one_is_certified',              &
                   distance_from_certification(rows) .lt. 1.0d0,          &
                   .true., nf)
      rows = (/ 0.5d0, 2.0d0*tel, 0.5d0*tca /)
      call log_row('judged_distance_above_one_is_not_certified',          &
                   distance_from_certification(rows) .lt. 1.0d0,          &
                   .false., nf)
      ! THE SUPERSEDED TOLERANCE, transcribed: an element row of 4e-8 was
      ! four times its tolerance and refused a state; on the anchored
      ! value it is 4e-3 of it, which is what the measurement of the
      ! discretization (4.6e-4 at the binding radius) says such a row is.
      rows = (/ 0.0d0, 4.0d-8, 0.0d0 /)
      call log_row('superseded_tolerance_refused_a_row_the_grid_'//       &
                   'cannot_resolve',                                      &
                   4.0d-8/1.0d-8 .lt. 1.0d0, .false., nf)
      call log_row('and_the_anchored_one_does_not',                       &
                   distance_from_certification(rows) .lt. 1.0d0,          &
                   .true., nf)
      ! The hydrodynamic distance reads each row against ITS OWN
      ! certification tolerance: one cell whose mass row is 1e-9 beside a
      ! momentum row of 1e-13 and an energy row of 1e-8 is 1e-9/3e-12 of
      ! the way, not 1e-8/1e-8 (the "Resid tol" reading this replaces,
      ! under which the HD 209458 b element reload's mass row at 1e-9 read
      ! as certified by the ledger while the acceptance refused it,
      ! Update_EXHALE_stage2.md section 8, P10). The cell's rounding floor is
      ! 1e-30, ten times which is far below 3e-12, so the tolerance that
      ! applies to its continuity row is the fixed one.
      q(1,1) = 1.0d-9;  q(2,1) = 1.0d-13;  q(3,1) = 1.0d-8
      fl(1)  = 1.0d-30
      call real_row('hydrodynamic_distance_reads_each_row_against_its_'// &
                    'own_tolerance',                                      &
                    hydrodynamic_distance_from_certification_by_cell(     &
                       1, q, fl),                                         &
                    1.0d-9/3.0d-12, nf)
      end subroutine judged_distance_is_the_worst_row_over_its_tolerance

      ! ================================================================= !

      subroutine advective_term_is_the_stage_operator(nf)
      ! THE ADVECTIVE TERM OF A CARRIER ROW IS THE OPERATOR THE RUNGE-KUTTA
      ! STAGES APPLY (item B5b).
      !
      ! One synthetic column, one carrier, no chemistry.  The rate of change
      ! a stage-1 species update produces from a given face mass flux is
      ! compared with the advective term carrier_advective_divergence forms
      ! on the same state: they are the same term of the same equation and
      ! must agree to the bit.
      !
      ! THE COLUMN IS CHOSEN WHERE THE TWO FORMS CANNOT AGREE.  The face
      ! mass flux changes sign inside the column, so there is a cell whose
      ! own velocity is zero while its two faces carry gas in opposite
      ! directions.  The cell-velocity upwind difference the stationary rows
      ! carried until this item, n_tot v (f_j - f_j-1)/dr, is identically
      ! zero in that cell whatever its neighbors hold; the flux divergence
      ! is not.  Both numbers are reported, so the row states the size of
      ! the disagreement it removed and not only that it is gone.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:), allocatable :: Frho, msum, rho_c, dt_loc, Yin
      real*8, dimension(:,:), allocatable :: Ycar, Ynew, fc
      real*8, dimension(:,:), allocatable :: adv, advmag
      real*8  :: dtv, mc, worst, worst_cell_velocity, rate, cadv, ntv
      integer :: j, jz

      N  = nc
      n0 = 1.0d0
      v0 = 1.0d0
      R0 = 1.0d0
      use_plm = .true.;  use_weno3 = .false.;  rec_method = 'PLM'
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_He           = .false.
      he_diffusion         = .false.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()

      if (allocated(r)) deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j)) deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r_edg(j) = 1.0d0 + 0.05d0*dble(j)
      enddo
      do j = 1-Ng, N+Ng
         r(j)    = r_edg(j) - 0.025d0
         dr_j(j) = 0.05d0
      enddo

      allocate(Frho(1-Ng:N+Ng), msum(1-Ng:N+Ng), rho_c(1-Ng:N+Ng),       &
               dt_loc(1-Ng:N+Ng), Yin(1-Ng:N+Ng))
      allocate(Ycar(1-Ng:N+Ng,1), Ynew(1-Ng:N+Ng,1))
      allocate(fc(1-Ng:N+Ng,n_carrier_max))
      allocate(adv(1:N,n_carrier_max), advmag(1:N,n_carrier_max))

      mc  = carrier_mass_amu(ic_H2)
      dtv = 1.0d-2
      jz  = nc/2
      do j = 1-Ng, N+Ng
         ! A face mass flux that changes sign at face jz, so cell jz has an
         ! inflow at both faces and no velocity of its own.
         Frho(j)   = 0.3d0*(dble(jz) - dble(j) + 0.5d0)
         msum(j)   = 1.3d0
         rho_c(j)  = 1.0d0
         dt_loc(j) = dtv
         fc(j,:)   = 0.0d0
         fc(j,ic_H2) = 0.05d0 + 0.02d0*sin(0.3d0*dble(j))
      enddo
      Ycar(:,1) = mc*fc(:,ic_H2)/msum
      ! The inner ghosts of a carrier no handoff states carry the base
      ! cell's own partition, which is what carrier_mass_fractions writes
      ! for the stages and what carrier_face_mass_fraction writes for the
      ! row.  The outer ghosts carry the outflow continuation of cell N,
      ! which is the boundary rule carrier_face_mass_fraction states for
      ! the row.  Writing both here is what makes the two comparable: what
      ! this row compares is the two forms of the TERM, on one state with
      ! one boundary rule, and not the boundary rule itself.
      do j = 1-Ng, 0
         Ycar(j,1) = Ycar(1,1)
      enddo
      do j = N+1, N+Ng
         Ycar(j,1) = Ycar(N,1)
      enddo

      call carrier_advective_divergence(fc, msum, Frho, adv, advmag)
      call species_advective_update(1, 1, 0, Ycar, Ycar, Frho, rho_c,    &
                                    rho_c, rho_c, dt_loc, 1, Ynew)

      ! The stage's rate of change of the carrier, in the row's own units.
      worst = 0.0d0
      worst_cell_velocity = 0.0d0
      do j = 1, N
         rate = -(Ynew(j,1) - Ycar(j,1))/dtv*msum(j)/mc
         worst = max(worst, abs(rate - adv(j,ic_H2))                     &
                            /max(abs(rate), 1.0d-300))
         ! The form the stationary rows carried before this item, on the
         ! same state: the cell velocity is the face flux averaged onto the
         ! cell, which is what n_tot v was.
         ntv = 0.5d0*(Frho(j) + Frho(j-1))
         if (ntv .ge. 0.0d0) then
            cadv = ntv/dr_j(j)*(fc(j,ic_H2) - fc(j-1,ic_H2))
         else if (j .lt. N) then
            cadv = ntv/dr_j(j)*(fc(j+1,ic_H2) - fc(j,ic_H2))
         else
            cadv = 0.0d0
         endif
         worst_cell_velocity = max(worst_cell_velocity,                  &
                                   abs(cadv - rate)                      &
                                   /max(abs(rate), 1.0d-300))
      enddo
      ! The bound is not zero because the reference is formed by
      ! DIFFERENCING a state: (Y_new - Y_in)/dt loses the digits the two
      ! states share, which at this step is about three.  What is asserted
      ! is that nothing but that cancellation separates the two.
      call bound_row('advective_term_is_the_stage_update', worst,        &
                     1.0d-13, nf)
      write(*,'(A,ES11.4)') '  the cell-velocity form this replaced'//    &
           ' differs from the stage by, at most, a relative ',           &
           worst_cell_velocity

      deallocate(Frho, msum, rho_c, dt_loc, Yin, Ycar, Ynew, fc,         &
                 adv, advmag)
      end subroutine advective_term_is_the_stage_operator

      ! ================================================================= !

      subroutine elemental_relaxation_fixed_point_is_the_row(nf, metals)
      ! THE FIXED-WIND ELEMENTAL RELAXATION AND THE STATIONARY ELEMENTAL ROW
      ! ARE ONE OPERATOR (item B5c).
      !
      ! The outer pass of the direct-steady route alternates a wind solved at
      ! a fixed composition with a composition relaxed at a fixed wind
      ! (relax_element_composition, which drives element_diffusion_step at a
      ! growing step until the composition stops moving), and then measures
      ! the state with the stationary elemental row
      ! (element_transport_residual). If those two carry different
      ! discretizations of the material advection they have different fixed
      ! points, and a converged alternation lands on a state the row does not
      ! read as stationary however long it runs.
      !
      ! What is asserted here is that they do not: relaxed at a GIVEN face
      ! mass flux until the helium mass fraction stops moving, the state
      ! makes the row of that same flux read its own convergence level. The
      ! step is driven here at a schedule written out in the open, so that
      ! what this row asserts is the operator and not a step schedule;
      ! relax_element_composition, which drives the same step at the
      ! composition time scales of the column, is asserted by
      ! relaxation_drift_covers_every_element below. The step, the flux it
      ! is driven with and the row that measures it are the production ones.
      !
      ! The column is stratified so the settling has something to work
      ! against, and the second row of each pair measures the same relaxed
      ! state against a row with NO advective term: the difference of the two
      ! row SCALES is the size of the advective term itself, which is what
      ! says the assertion is not vacuous.
      !
      ! metals runs the same column with trace-element diffusion on
      ! (He_metal_diffusion), which is the only path that exercises the
      ! deferred correction of solve_trace_element_in_hydrogen; its rows
      ! carry the worst trace element of the column.
      integer, intent(inout) :: nf
      logical, intent(in)    :: metals
      integer, parameter :: nc = 60
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, dt_c, Frho
      real*8, dimension(:),   allocatable :: Fzero, Xnow, Xprev
      real*8, dimension(:,:), allocatable :: f_c, Yel
      real*8, dimension(:),   allocatable :: res_he, sc_he, sc_wind, adv_new
      real*8, dimension(:,:), allocatable :: res_tr, sc_tr
      real*8  :: dr_u, m_1, mpH, X_base, grow, row, adv_share, heh_j
      real*8  :: famp, gap, rhov_old, old_term, row_tr
      character(len=16) :: tag
      logical :: ok_he, ok_tr
      integer :: j, k, nrel, ie, im

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      q0 = n0*mu*v0*v0*v0/R0
      b0 = 5.0d0
      spherical_domain   = .true.
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_mol        = .false.
      carrier_transport  = .false.
      thereis_metals     = metals
      eos_include_metals = .true.
      he_diffusion       = .true.
      he_metal_diffusion = metals
      he_ambipolar       = .false.
      he_alphaT          = 0.0d0
      he_kzz             = 1.0d11
      HeH                = 0.0833333333333333d0
      use_plm = .true.;  use_weno3 = .false.;  rec_method = 'PLM'
      recon_lambda_on = .false.
      call carrier_set_init()

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr_u = 1.0d0/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      if (metals) melem_ab = 1.0d-4

      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               dt_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng), Fzero(1-Ng:N+Ng),        &
               Xnow(1-Ng:N+Ng), Xprev(1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), Yel(1-Ng:N+Ng,1+n_melem))
      allocate(res_he(1:N), sc_he(1:N), sc_wind(1:N), adv_new(1:N),       &
               res_tr(1:N,n_melem), sc_tr(1:N,n_melem))

      m_1 = mass_per_H_nucleus_without_He()
      f_c = 0.0d0
      ! A face mass flux with r^2 F_rho constant, which is what a steady wind
      ! carries. It is not a free choice: div(F_rho X) = X div(F_rho)
      ! + F_rho grad X, so on a flux that does NOT conserve mass the term
      ! X div(F_rho) is a source of the element with no sink, and the
      ! conservative balance has no bounded steady composition at all (the
      ! relaxation drives X onto 0 and 1 and stays there). That is a property
      ! of the equation and not of the operator, and it is why the wind the
      ! outer pass hands over is the object this row is measured on.
      ! The amplitude is chosen so the advective term is a comparable part of
      ! the balance to the settling; the fraction reported below measures it.
      famp = 1.0d-4
      do j = 1-Ng, N+Ng
         ! A composition that is not the reservoir one anywhere but the base,
         ! so the relaxation has a distance to travel.
         heh_j = HeH*(1.0d0 + 0.5d0*sin(6.0d0*(r(j) - 1.0d0)))
         mpH   = m_1 + m_He_over_m_H*heh_j
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = heh_j/mpH
         ! The metals ride on hydrogen at the reservoir metal/H, in their
         ! neutral stage, which is the state set_IC and the write-back build
         ! an absent element from. m_1 already carries their mass.
         if (metals) then
            do ie = 1, n_melem
               f_c(j,mion_fsp(melem_i0(ie))) = melem_ab(ie)*f_c(j,isp_HI)
            enddo
         endif
         rho_c(j)  = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)    = 0.0d0
         T_c(j)    = 1.0d0
         Frho(j)   = famp/(r_edg(j)*r_edg(j))
      enddo
      Fzero = 0.0d0

      X_base = reservoir_helium_mass_fraction(m_1)
      call element_mass_fractions(f_c, Yel)
      Xprev = Yel(:,1)

      ! The relaxation: the composition time scale, grown geometrically, so
      ! the late steps are direct steady solves. Same schedule as
      ! relax_element_composition.
      dt_c = 1.0d-6
      grow = 1.0d0
      nrel = 0
      do k = 1, 400
         call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c*grow,     &
                                     Frho_in = Frho)
         call element_mass_fractions(f_c, Yel)
         Xnow = Yel(:,1)
         nrel = k
         if (maxval(abs(Xnow(1:N) - Xprev(1:N)))/X_base .lt. 1.0d-12) exit
         Xprev = Xnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo

      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he,       &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      row = 0.0d0
      do j = 2, N
         row = max(row, abs(res_he(j))/max(sc_he(j), 1.0d-300))
      enddo
      sc_wind = sc_he
      ! The worst trace element of the column, over every element and cell.
      ! The worst trace element of the column, over every element and cell.
      ! THE COLUMN IS MIXED (K_zz = 1e11 against D_12 ~ 5e9 here) SO THAT NO
      ! ELEMENT VANISHES. Driven harder, the heaviest elements settle out of
      ! the top cells entirely, and there the write-back of
      ! element_diffusion_step takes its vanished-element branch
      ! (nXold <= 1e-25 n_H) instead of scaling the solved mixing ratio into
      ! f_sp, so the state measured is not the one the solve returned:
      ! MEASURED on the same column at K_zz = 0, Ca, K and Fe reach exactly
      ! zero in the outermost cells and their rows there read 1e-4 to 2e-2
      ! while every element that is present reads 1e-15. That is the write-
      ! back's own guard and not this operator.
      row_tr = 0.0d0
      if (metals) then
         do im = 1, n_melem
            do j = 2, N
               row_tr = max(row_tr, abs(res_tr(j,im))                     &
                                    /max(sc_tr(j,im), 1.0d-300))
            enddo
         enddo
      endif

      ! The same row on the same state with NO advective term. Its scale is
      ! the diffusive terms alone, so the difference of the two scales IS the
      ! size of the advective term, cell by cell: what the fraction below
      ! says is how much of the balance the row above reads to 1e-13 the
      ! advection is carrying. Without it the assertion would hold for a
      ! column the wind never touched.
      call element_transport_residual(rho_c, T_c, f_c, Fzero, res_he,      &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      adv_share = 0.0d0
      do j = 2, N
         adv_share = max(adv_share, (sc_wind(j) - sc_he(j))               &
                                    /max(sc_wind(j), 1.0d-300))
         ! The advective term the new form carries at this state: the row
         ! with no advection is the diffusive divergence alone, and the full
         ! row is zero to the bound asserted below, so the advective term is
         ! minus that.
         adv_new(j) = -res_he(j)
      enddo

      ! WHAT THE FORM THIS REPLACED WOULD READ ON THE SAME STATE. The
      ! relaxation used to advect as rho v dX/dr, one-sided upwind on the
      ! smoothed steady mass flux mdot/(4 pi r^2) -- which for this column IS
      ! the wind, since r^2 F_rho is constant, so the two differ by the
      ! discretization alone. Written out here because the row above is a
      ! bound and this is the size of the disagreement it removed: it is the
      ! residual the OLD operator leaves at the fixed point of the new one,
      ! on the new one's own row scale.
      gap = 0.0d0
      do j = 2, N
         rhov_old = famp/(r(j)*r(j))*n0*mu*v0
         old_term = rhov_old*(Xnow(j) - Xnow(j-1))                        &
                    /max((r(j) - r(j-1))*R0, 1.0d0)
         ! res_he is zero at this state to the bound above, so the advective
         ! term of the new form is minus the diffusive divergence, and the
         ! old form's residual there is its own advective term minus that.
         gap = max(gap, abs(old_term - adv_new(j))                        &
                        /max(sc_wind(j), 1.0d-300))
      enddo

      tag = '_helium_only'
      if (metals) tag = '_with_metals'

      write(*,'(A,I0,A,ES11.4,A,ES11.4)')                                 &
           '  the relaxation converged in ', nrel,                        &
           ' steps; X runs from ', minval(Xnow(1:N)), ' to ',             &
           maxval(Xnow(1:N))
      write(*,'(A,ES11.4)') '  the cell-velocity form on the smoothed'//   &
           ' steady flux leaves this state at a relative ', gap
      ! The bound is not the arithmetic floor: the relaxation stops on the
      ! composition and the row is measured after it, so what is left is the
      ! step's own Newton tolerance. It is decades below what a relaxation
      ! carrying a SECOND discretization of this term leaves behind, which is
      ! the number printed above.
      call bound_row('elemental_relaxation_fixed_point_is_the_row'//      &
                     trim(tag), row, 1.0d-6, nf)
      if (metals)                                                         &
         call bound_row('trace_relaxation_fixed_point_is_the_row', row_tr,&
                        1.0d-6, nf)
      if (adv_share .gt. 0.1d0) then
         write(*,'(A,A,A,A,ES12.5,A)') 'PASS ',                           &
              'elemental_row_carries_the_advection', trim(tag),           &
              ' measured=', adv_share, ' reference=>0.1 tol=0'
      else
         write(*,'(A,A,A,A,ES12.5,A)') 'FAIL ',                           &
              'elemental_row_carries_the_advection', trim(tag),           &
              ' measured=', adv_share, ' reference=>0.1 tol=0'
         nf = nf + 1
      endif

      deallocate(rho_c, v_c, T_c, dt_c, Frho, Fzero, Xnow, Xprev, f_c,    &
                 Yel, res_he, sc_he, sc_wind, adv_new, res_tr, sc_tr)
      end subroutine elemental_relaxation_fixed_point_is_the_row

      ! ================================================================= !

      subroutine relaxation_drift_covers_every_element(nf)
      ! THE FIXED-WIND RELAXATION REPORTS THE DISTANCE EVERY ELEMENT STILL
      ! HAS TO TRAVEL, NOT HELIUM'S ALONE (item DIFT-LINK, from B5c 6.3).
      !
      ! relax_element_composition returns one number, and the outer Picard
      ! loop of the direct-steady route stops, damps and reports on it. If
      ! that number is helium's drift alone, a column whose trace metals
      ! settle on a longer time scale than its helium is declared converged
      ! while the metals are still moving: the loop exits with an elemental
      ! composition that is not a fixed point of the operator that produced
      ! it, and the stationary trace rows measured afterwards are rows of a
      ! state nothing claimed was stationary.
      !
      ! The column here puts helium at its reservoir value everywhere, so it
      ! has almost nothing to travel, and starts every trace element at a
      ! mixing ratio that swings +/-80% about its own reservoir, so the
      ! metals have a great deal. What is asserted is that the returned
      ! drift bounds every element's own change over the relaxation, each
      ! element measured in its own reservoir value; the helium-only number
      ! is printed beside it, and it is the one that does not.
      !
      ! The swing is taken OUT OF THE HYDROGEN it is measured against, so
      ! the column's species reconstruct the density it is handed (asserted
      ! before the relaxation, `drift_column_carries_its_own_density`).
      ! Added on top of a fixed hydrogen fraction it left the species short
      ! of that density by 1.7e-2, and a drift measured on ratios that do
      ! not describe the transported mass is a drift of no state.
      !
      ! omega = 1 so that the state handed back IS the relaxed one and the
      ! drift is exactly the distance between the two compositions the test
      ! measures.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 60
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c, Ybef, Yaft
      real*8, dimension(:),   allocatable :: Yres, dchange
      real*8  :: dr_u, m_1, X_base, famp, drift, worst, drift_he
      real*8  :: worst_trace, closure
      integer :: j, ie, nrel, iworst

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      q0 = n0*mu*v0*v0*v0/R0
      b0 = 5.0d0
      spherical_domain   = .true.
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_mol        = .false.
      carrier_transport  = .false.
      thereis_metals     = .true.
      eos_include_metals = .true.
      he_diffusion       = .true.
      he_metal_diffusion = .true.
      he_ambipolar       = .false.
      he_alphaT          = 0.0d0
      he_kzz             = 1.0d11
      HeH                = 0.0833333333333333d0
      use_plm = .true.;  use_weno3 = .false.;  rec_method = 'PLM'
      recon_lambda_on = .false.
      call carrier_set_init()

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr_u = 1.0d0/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4

      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), Ybef(1-Ng:N+Ng,1+n_melem),       &
               Yaft(1-Ng:N+Ng,1+n_melem))
      allocate(Yres(1+n_melem), dchange(1+n_melem))

      m_1  = mass_per_H_nucleus_without_He()
      famp = 1.0d-4
      ! The species as mass fractions that sum to one, so the column
      ! carries its own density exactly: helium at the reservoir value
      ! everywhere, where the base pins it and it has almost no distance to
      ! travel, and every metal swinging +/-80% about its own reservoir
      ! mixing ratio in its neutral stage, the swing taken off the hydrogen
      ! each element's mixing ratio is measured against. Cell 1 is the
      ! Dirichlet reservoir and stays at the reservoir value, which is what
      ! the relaxation measures each element against.
      call column_carrying_its_own_density(f_c, 0.0d0, .false., 0.8d0)
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         ! A face mass flux with r^2 F_rho constant, which is what a steady
         ! wind carries; on a flux that does not conserve mass the
         ! conservative balance has no bounded steady composition at all.
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo

      ! WHAT THE COLUMN WEIGHS AGAINST THE DENSITY IT IS HANDED, by the
      ! code's own mass policy (calc_rho, metals in the budget). The metal
      ! swing above is taken out of the hydrogen it is measured against, so
      ! the species reconstruct rho to round-off; a swing added on top of a
      ! fixed hydrogen fraction would make the element ratios the
      ! relaxation reads and the density it transports two different
      ! atmospheres, and the drift it returns a drift of neither.
      closure = column_mass_closure(rho_c, f_c)
      call bound_row('drift_column_carries_its_own_density', closure,     &
                     1.0d-14, nf)

      X_base = reservoir_helium_mass_fraction(m_1)
      call element_mass_fractions(f_c, Ybef)
      call relax_element_composition(rho_c, v_c, T_c, f_c, Frho, 1.0d0,   &
                                     drift, nrel)
      call element_mass_fractions(f_c, Yaft)

      ! Each element's own change over the relaxation, in its own reservoir
      ! value: helium in the input He/H, a metal in its base cell.
      Yres    = 0.0d0
      Yres(1) = X_base
      do ie = 1, n_melem
         Yres(1+ie) = Ybef(1,1+ie)
      enddo
      dchange = 0.0d0
      do ie = 1, 1 + n_melem
         if (Yres(ie) .le. 0.0d0) cycle
         dchange(ie) = maxval(abs(Yaft(1:N,ie) - Ybef(1:N,ie)))/Yres(ie)
      enddo
      drift_he    = dchange(1)
      worst_trace = 0.0d0
      iworst      = 0
      do ie = 1, n_melem
         if (dchange(1+ie) .gt. worst_trace) then
            worst_trace = dchange(1+ie)
            iworst      = ie
         endif
      enddo

      ! The bound: no element moved further than the number the relaxation
      ! handed the outer loop. Measured as the largest ratio of an element's
      ! own change to that number, which is 1 for the element that set it.
      worst = 0.0d0
      do ie = 1, 1 + n_melem
         if (Yres(ie) .le. 0.0d0) cycle
         worst = max(worst, dchange(ie)/max(drift, 1.0d-300))
      enddo

      write(*,'(A,I0,A,ES11.4)') '  the relaxation took ', nrel,          &
           ' steps and returned a drift of ', drift
      if (iworst .gt. 0)                                                  &
         write(*,'(A,A,A,ES11.4,A,ES11.4)') '  the worst trace element (', &
              trim(melem_name(iworst)), ') moved by ', worst_trace,       &
              ' while helium moved by ', drift_he
      call bound_row('relaxation_drift_bounds_every_element',             &
                     max(0.0d0, worst - 1.0d0), 1.0d-12, nf)
      ! THE GUARD AGAINST VACUITY, and it has no free parameter: the row
      ! above is satisfied by a helium-only drift on any column whose helium
      ! is the element that moves furthest. What states that this column is
      ! not one of those is that a trace element outran helium, which is
      ! exactly the condition under which the measure this replaced returns
      ! a number that does not bound every element.
      if (worst_trace .gt. drift_he) then
         write(*,'(A,A,ES12.5,A)') 'PASS ',                               &
              'relaxation_helium_only_drift_would_not_bound measured=',   &
              worst_trace/max(drift_he, 1.0d-300), ' reference=>1 tol=0'
      else
         write(*,'(A,A,ES12.5,A)') 'FAIL ',                               &
              'relaxation_helium_only_drift_would_not_bound measured=',   &
              worst_trace/max(drift_he, 1.0d-300), ' reference=>1 tol=0'
         nf = nf + 1
      endif

      deallocate(rho_c, v_c, T_c, Frho, f_c, Ybef, Yaft, Yres, dchange)
      end subroutine relaxation_drift_covers_every_element

      ! ================================================================= !

      subroutine eq_measure_reads_what_the_residual_reads(nf)
      ! THE COMPOSITION ELIMINATION'S STOPPING TEST MEASURES WHAT THE
      ! RESIDUAL READS THE COMPOSITION FOR (item B5j part 2).
      !
      ! The test this replaced was the largest relative change of any
      ! SPECIES FRACTION of any cell, with an absolute floor of 1e-12 under
      ! its denominator. In the fully ionized wind the He I and H2+
      ! fractions are ~1e-20, so their denominator is the floor and the
      ! quotient reports the round-off of a species eight decades below it.
      ! MEASURED on the molecular hot Uranus (report B5h section 2): the
      ! increment plateaus at 6e-11 for sixty passes at a fixed count, so
      ! `eq_sweep_reltol` below about 1e-10 is unreachable by construction.
      !
      ! THE SYNTHETIC COLUMN below is that situation with nothing else in
      ! it: a species whose fraction is 1e-20 changes by 6e-23, i.e. by
      ! 0.6 percent of itself, and nothing else moves. Its contribution to
      ! the particle count is 1e-20 of it, so the particle count moves by
      ! 6e-23 relative -- which is what the residual can see, and it is
      ! eleven decades below the tightest tolerance anyone would ask for.
      ! The expression this replaced reads the same perturbation as
      ! 6e-23/1e-12 = 6e-11, above 1e-12, which is the whole difference and
      ! is asserted here as well so that the two numbers stand side by side.
      !
      ! The three channels are then exercised one at a time, because a
      ! measure that saw only one of them would pass the first row and still
      ! be blind to a pass that redistributed the composition at fixed
      ! particle count.
      integer, intent(inout) :: nf
      real*8, allocatable :: npart(:), npart_prev(:), xh2(:), xh2_prev(:)
      real*8, allocatable :: heat(:), cool(:), heat_prev(:), cool_prev(:)
      real*8  :: d, d_old, f_tiny, df_tiny
      integer :: jm, km, nsave

      nsave = N
      N = 5
      allocate(npart(1:N), npart_prev(1:N), xh2(1:N), xh2_prev(1:N))
      allocate(heat(1:N), cool(1:N), heat_prev(1:N), cool_prev(1:N))

      ! A species at 1e-20 of the mass, changed by 0.6 percent of itself.
      f_tiny  = 1.0d-20
      df_tiny = 6.0d-23
      npart      = 1.0d0
      npart_prev = 1.0d0
      npart_prev(3) = 1.0d0 - df_tiny
      xh2 = 0.0d0;  xh2_prev = 0.0d0
      heat = 1.0d0;  cool = 1.0d0
      heat_prev = 1.0d0;  cool_prev = 1.0d0
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .true.,    &
              d, jm, km)
      call bound_row('eq_measure_ignores_a_species_the_residual_cannot'// &
                     '_see', d, 1.0d-12, nf)
      d_old = df_tiny/max(f_tiny, 1.0d-12)
      if (d_old .gt. 1.0d-12) then
         write(*,'(A,A,ES12.5,A)') 'PASS ', 'eq_measure_species'//        &
              '_fraction_expression_would_not measured=', d_old,          &
              ' reference=>1e-12 tol=0'
      else
         write(*,'(A,A,ES12.5,A)') 'FAIL ', 'eq_measure_species'//        &
              '_fraction_expression_would_not measured=', d_old,          &
              ' reference=>1e-12 tol=0'
         nf = nf + 1
      endif

      ! EVERY PERTURBATION BELOW IS A NEGATIVE POWER OF TWO of a number
      ! that is itself exact in binary, so the measure's value is exact and
      ! the rows carry no tolerance: a perturbation of 1e-6 subtracted from
      ! one would arrive 3e-11 out relative and the row would then be
      ! measuring the construction of its own input.
      !
      ! Channel 1: the particle count, which sets T = p/n_part.
      npart_prev = 1.0d0
      npart_prev(2) = 1.0d0 - 1.0d0/1048576.0d0            ! 1 - 2^-20
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .true.,    &
              d, jm, km)
      call real_row('eq_measure_particle_count_value', d,                 &
                    1.0d0/1048576.0d0, nf)
      call int_row('eq_measure_particle_count_channel', km,               &
                   eq_channel_particle_count, nf)
      call int_row('eq_measure_particle_count_cell', jm, 2, nf)

      ! Channel 2: the H2 share of the particle count, at FIXED particle
      ! count -- a redistribution the first channel cannot see.
      npart_prev = 1.0d0
      xh2_prev = 0.5d0;  xh2 = 0.5d0
      xh2(4) = 0.5d0 + 1.0d0/4194304.0d0                   ! 0.5 + 2^-22
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .true.,    &
              d, jm, km)
      call real_row('eq_measure_h2_share_value', d,                       &
                    1.0d0/4194304.0d0, nf)
      call int_row('eq_measure_h2_share_channel', km,                     &
                   eq_channel_h2_caloric, nf)

      ! Channel 3: the radiative pair, relative to the larger of the two
      ! terms the energy row is the difference of -- here heat = 2, so the
      ! change of the difference is halved by its own scale.
      xh2 = 0.5d0;  xh2_prev = 0.5d0
      heat = 2.0d0;  heat_prev = 2.0d0
      cool = 1.0d0;  cool_prev = 1.0d0
      cool(5) = 1.0d0 + 1.0d0/268435456.0d0                ! 1 + 2^-28
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .true.,    &
              d, jm, km)
      call real_row('eq_measure_radiative_value', d,                      &
                    0.5d0/268435456.0d0, nf)
      call int_row('eq_measure_radiative_channel', km,                    &
                   eq_channel_radiative, nf)

      ! A cell whose two radiative terms are 1e-40 of the column's largest
      ! cannot move the energy row whatever they do relative to each other,
      ! and the floor of the measure says so.
      xh2 = 0.5d0;  xh2_prev = 0.5d0
      heat = 1.0d0;  cool = 1.0d0
      heat_prev = 1.0d0;  cool_prev = 1.0d0
      heat(5) = 1.0d-40;  cool(5) = 1.0d-40
      heat_prev(5) = 1.0d-40;  cool_prev(5) = 5.0d-41
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .true.,    &
              d, jm, km)
      call bound_row('eq_measure_floors_a_radiatively_dead_cell', d,      &
                     1.0d-12, nf)

      ! THE FIRST PASS HAS NO PREVIOUS RADIATIVE FIELD, and the measure then
      ! reports the two equation-of-state channels alone rather than the
      ! difference against an array nothing wrote.
      heat = 1.0d0;  cool = 3.0d0
      heat_prev = 0.0d0;  cool_prev = 0.0d0
      call increment_of_what_the_residual_reads(npart, npart_prev,        &
              xh2, xh2_prev, heat, cool, heat_prev, cool_prev, .false.,   &
              d, jm, km)
      call real_row('eq_measure_first_pass_has_no_radiative_channel', d,  &
                    0.0d0, nf)

      deallocate(npart, npart_prev, xh2, xh2_prev)
      deallocate(heat, cool, heat_prev, cool_prev)
      N = nsave
      end subroutine eq_measure_reads_what_the_residual_reads

      ! ================================================================= !

      subroutine species_box_is_a_box_and_not_a_ball(nf)
      ! A SPECIES UNKNOWN IS BOUNDED AND THE TRIAL IS WRITTEN ONTO ITS FACE
      ! (item B5j part 1).
      !
      ! A carrier unknown is a density in the code's units, so it lies
      ! between nothing and the cell's own density; an element unknown is a
      ! mass fraction and lies in [0,1]. Shortening the whole step against
      ! such a bound cannot recover a component that is already ON it, which
      ! is what the trust region was doing: MEASURED on the coupled
      ! `mol_carrier` reload, the dogleg cut its trial back by 2^-25 to
      ! 2^-49 against n(H2) of cell 500 and 49 iterations moved ||R|| from
      ! 1.873 to 1.873 (report B5g section 6.5).
      !
      ! These rows state the projection itself: which unknowns are counted
      ! as outside, that the projected vector is inside, and that the
      ! hydrodynamic slots are never touched by it.
      integer, intent(inout) :: nf
      real*8, allocatable :: Yv(:), Ykeep(:)
      integer :: nout, nsave, j, is

      nsave = N
      N = 4
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_He           = .false.
      he_diffusion         = .false.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call set_transported_species_rows(.true.)
      allocate(Yv(nvar_jac*N), Ykeep(nvar_jac*N))
      ! Densities 1, 2, 3, 4 in the mass slots; the carrier of every cell
      ! outside its own box in a different way.
      do j = 1, N
         Yv(nvar_jac*(j-1)+1) = dble(j)
         Yv(nvar_jac*(j-1)+2) = -7.0d0     ! momentum: signed, never bounded
         Yv(nvar_jac*(j-1)+3) = dble(j)
      enddo
      Yv(nvar_jac*0+4) = -1.0d-19          ! below the bound
      Yv(nvar_jac*1+4) =  1.0d0            ! inside
      Yv(nvar_jac*2+4) =  9.0d0            ! above the cell's own density 3
      Yv(nvar_jac*3+4) =  0.0d0            ! on the bound
      Ykeep = Yv
      ! The faces are frozen from the state before any of them is read
      ! (freeze_species_unknown_box).  Nothing here has assembled a carrier
      ! state, so there is no element budget and the upper face is the
      ! cell's own density, which is the box these rows were written for.
      call freeze_species_unknown_box(Yv)
      call log_row('species_box_has_no_element_budget_yet',              &
                   carrier_headroom_known(), .false., nf)
      call species_unknowns_outside_their_bounds(Yv, nout, .false.)
      call int_row('species_box_counts_the_two_outside', nout, 2, nf)
      call real_row('species_box_counting_changes_nothing',              &
                    maxval(abs(Yv - Ykeep)), 0.0d0, nf)
      call species_unknowns_outside_their_bounds(Yv, nout, .true.)
      call real_row('species_box_projects_below_onto_zero',              &
                    Yv(nvar_jac*0+4), 0.0d0, nf)
      call real_row('species_box_projects_above_onto_the_density',       &
                    Yv(nvar_jac*2+4), 3.0d0, nf)
      call real_row('species_box_leaves_the_inside_alone',               &
                    Yv(nvar_jac*1+4), 1.0d0, nf)
      call real_row('species_box_leaves_the_bound_alone',                &
                    Yv(nvar_jac*3+4), 0.0d0, nf)
      call species_unknowns_outside_their_bounds(Yv, nout, .false.)
      call int_row('species_box_projected_vector_is_inside', nout, 0, nf)
      ! The hydrodynamic slots, including the signed momentum, are untouched.
      do j = 1, N
         do is = 1, 3
            if (Yv(nvar_jac*(j-1)+is) .ne. Ykeep(nvar_jac*(j-1)+is)) then
               nf = nf + 1
               write(*,'(A,I0,A,I0)') 'FAIL species_box_touched_a_'//    &
                    'hydrodynamic_slot measured=cell ', j, ' row ', is
            endif
         enddo
      enddo
      write(*,'(A)') 'PASS species_box_leaves_the_hydrodynamic_slots'//  &
           ' measured=0 reference=0 tol=0'
      deallocate(Yv, Ykeep)
      call set_transported_species_rows(.false.)
      N = nsave
      end subroutine species_box_is_a_box_and_not_a_ball

      ! ================================================================= !

      subroutine element_budget_is_a_row_of_the_step_and_not_a_face(nf)
      ! A CARRIER CANNOT EXCEED THE ELEMENT IT IS MADE OF, AND THAT CEILING
      ! IS SHARED, SO IT IS A ROW OF THE STEP AND NOT A FACE OF THE UNKNOWN
      ! BOX (item N4b, decision 14 route (i);
      ! docs/constrained_step_design_20260909.md section 2).
      !
      ! Two carriers of one element spend one budget, so what bounds them is
      !
      !     sum_i nuclei_per_particle(i,ie) n_i  <=  budget(ie)
      !
      ! a half-space of the cell's unknowns. Each carrier's own ceiling
      ! budget/nuclei is a CORNER of that half-space, and the conjunction of
      ! the corners is a strictly larger set: item N4a measures a breach of
      ! 1.0 of the available hydrogen on a state sitting on two of them at
      ! once. The corners also cost the step whole directions, a component
      ! with no room between its value and its own ceiling making the
      ! fraction-to-the-boundary rule return zero for every component.
      !
      ! So the box keeps the faces that ARE coordinate bounds -- a carrier
      ! density is non-negative and cannot be more of the cell than the cell
      ! is -- and the shared budget is a row
      ! (freeze_element_constraint_rows). These rows state the box that is
      ! left, and the row that took the budget over.
      !
      ! THE BUDGET IS A PROPERTY OF AN ASSEMBLED CARRIER STATE, which a
      ! synthetic cell has none of, so it is written here
      ! (carrier_headroom_set_for_test). The rows below are in the DENSITY
      ! coordinate; in ln n it is the same row read through the logarithm,
      ! and that coordinate is measured on the coupled run itself.
      !
      ! RED AND GREEN, ONE ENVIRONMENT VARIABLE APART.
      ! EXHALE_ELEMENT_CONSTRAINT_ROWS=0 puts the budget back under a
      ! coordinate face and forms no row, which is the box the B5j to B5l
      ! measurements were made in; the rows below then read the budget where
      ! they read the cell's own density, and no row exists to read. The
      ! controls are read here, as a solve reads them, so that the switch
      ! reaches this subroutine.
      integer, intent(inout) :: nf
      real*8, allocatable :: Yv(:), Ykeep(:), Yit(:), Dsc(:)
      real*8, allocatable :: budget(:,:), f_c(:,:)
      real*8  :: nsave_n0, nH_avail
      integer :: nout, nsave, j, kH

      nsave = N
      nsave_n0 = n0
      N  = 4
      ! The density scale the carrier unknown is measured in: the budget is
      ! a density in cm^-3 and the unknown is that density over n0.
      n0 = 1.0d10
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_He           = .false.
      he_diffusion         = .false.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call read_species_unknown_space_controls
      call set_transported_species_rows(.true.)
      allocate(Yv(nvar_jac*N), Ykeep(nvar_jac*N), Yit(nvar_jac*N))
      allocate(Dsc(nvar_jac*N))
      allocate(budget(N,n_carrier_max), f_c(1-Ng:N+Ng,n_species))
      ! Half the hydrogen of a cell holding 1e10 cm^-3 of it: the side each
      ! carrier would have alone is 0.5 in the code's density units, and the
      ! shared hydrogen budget it comes from is twice that, H2 holding two
      ! nuclei.
      budget = 0.5d0*1.0d10
      call carrier_headroom_set_for_test(budget)
      call log_row('species_box_element_budget_is_known',                &
                   carrier_headroom_known(), .true., nf)
      ! The iterate: every cell's carrier well inside its budget.
      Dsc = 1.0d0
      f_c = 0.0d0
      do j = 1, N
         Yit(nvar_jac*(j-1)+1) = dble(j)     ! density 1, 2, 3, 4
         Yit(nvar_jac*(j-1)+2) = -7.0d0
         Yit(nvar_jac*(j-1)+3) = dble(j)
         Yit(nvar_jac*(j-1)+4) = 0.1d0       ! carrier, inside the budget
         f_c(j,isp_H2) = 0.1d0/dble(j)       ! the same carrier as a fraction
      enddo
      call freeze_species_unknown_box(Yit)
      call freeze_element_constraint_rows(Yit, f_c, Dsc)
      call log_row('element_budget_is_a_row_of_the_step',                &
                   element_constraint_rows_known(), .true., nf)
      ! ---- what the box is left with -------------------------------
      ! The trial: cell 1 and cell 3 far above their own SHARED budget but
      ! inside the cell they sit in, cell 2 inside everything, cell 4 below
      ! zero. Only the last leaves the box: the budget is not a face of it.
      Yv = Yit
      Yv(nvar_jac*0+4) = 0.9d0
      Yv(nvar_jac*1+4) = 0.2d0
      Yv(nvar_jac*2+4) = 0.9d0
      Yv(nvar_jac*3+4) = -1.0d-20
      Ykeep = Yv
      call species_unknowns_outside_their_bounds(Yv, nout, .false.)
      call int_row('species_box_counts_only_the_negative_carrier',       &
                   nout, 1, nf)
      call real_row('species_box_budget_counting_changes_nothing',       &
                    maxval(abs(Yv - Ykeep)), 0.0d0, nf)
      call species_unknowns_outside_their_bounds(Yv, nout, .true.)
      call real_row('species_box_leaves_a_carrier_above_its_own_side',   &
                    Yv(nvar_jac*0+4), 0.9d0, nf)
      call real_row('species_box_leaves_a_carrier_inside_the_budget',    &
                    Yv(nvar_jac*1+4), 0.2d0, nf)
      call real_row('species_box_projects_a_negative_carrier_onto_zero', &
                    Yv(nvar_jac*3+4), 0.0d0, nf)
      ! The face that IS a coordinate bound: a carrier cannot be more of the
      ! cell than the cell is. Cell 1 holds 1 in code density units.
      Yv(nvar_jac*0+4) = 4.0d0
      call species_unknowns_outside_their_bounds(Yv, nout, .true.)
      call real_row('species_box_upper_face_is_the_cell_density',        &
                    Yv(nvar_jac*0+4), 1.0d0, nf)
      ! ---- and what the row says -----------------------------------
      ! The hydrogen row is the first of the three (H, O, C).
      kH = 1
      nH_avail = 2.0d0*0.5d0*1.0d10
      call real_row('element_row_budget_is_the_operator_budget',         &
                    element_constraint_budget_of_the_cell(1, kH),        &
                    nH_avail, nf)
      ! c is the demand less the budget: cell 1 holds n(H2) = 0.1*n0 and H2
      ! holds two hydrogen nuclei, so the demand is 0.2*n0.
      call real_row('element_row_c_is_the_demand_less_the_budget',       &
                    element_constraint_of_the_cell(1, kH),               &
                    2.0d0*0.1d0*1.0d10 - nH_avail, nf)
      ! The carrier column of the row is the stoichiometry itself, in
      ! nuclei per unit unknown: two nuclei per H2 particle, and the unknown
      ! is a density in units of n0.
      call real_row('element_row_gradient_is_the_stoichiometry',         &
                    element_constraint_gradient_of_the_cell(1, kH, 4),   &
                    2.0d0*1.0d10, nf)
      ! The density column is c/nd, so it vanishes as the row becomes
      ! tight: the constraint is homogeneous of degree one in the density
      ! and the density direction is tangent to its surface.
      call real_row('element_row_density_column_is_c_over_nd',           &
                    element_constraint_gradient_of_the_cell(1, kH, 1),   &
                    (2.0d0*0.1d0*1.0d10 - nH_avail)/(1.0d0*1.0d10)       &
                    *1.0d10, nf)
      deallocate(Yv, Ykeep, Yit, Dsc, budget, f_c)
      call set_transported_species_rows(.false.)
      N  = nsave
      n0 = nsave_n0
      end subroutine element_budget_is_a_row_of_the_step_and_not_a_face

      ! ================================================================= !

      subroutine a_direction_along_a_shared_constraint_is_not_discarded(nf)
      ! ONE BLOCKED COMPONENT MAY NOT COST THE WHOLE DIRECTION (item N4b
      ! deliverable 2).
      !
      ! Two failures were measured and they are different failures.
      !   * The fraction-to-the-boundary rule returned ZERO for a whole
      !     direction as soon as ONE component had no room between its value
      !     and its own face, so the finite-difference probe had no step and
      !     the Krylov cycle no direction: MEASURED on `mol_diffusion`, 18
      !     cycles truncated for that reason. The rule now counts such a
      !     component, skips it, and returns the length the components with
      !     room allow; the caller zeroes the counted components before it
      !     probes.
      !   * A component blocked by a SHARED budget was zeroed, which threw
      !     away the direction's motion in every other carrier of the same
      !     element as well. The step is now projected onto the constraint's
      !     hyperplane, so what is removed is the one combination that
      !     spends nuclei the cell does not have and everything orthogonal
      !     to it survives.
      integer, intent(inout) :: nf
      real*8, allocatable :: Yit(:), v(:), vk(:), g(:), Dsc(:)
      real*8, allocatable :: budget(:,:), f_c(:,:)
      real*8  :: nsave_n0, tmax, len_in, len_out, comp_in, slack
      integer :: nsave, j, nblk, nact

      nsave = N;  nsave_n0 = n0
      N  = 3
      n0 = 1.0d10
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .true.
      ionization_transport = .false.
      thereis_He           = .false.
      he_diffusion         = .false.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call read_species_unknown_space_controls
      call set_transported_species_rows(.true.)
      allocate(Yit(nvar_jac*N), v(nvar_jac*N), vk(nvar_jac*N),           &
               g(nvar_jac*N), Dsc(nvar_jac*N))
      allocate(budget(N,n_carrier_max), f_c(1-Ng:N+Ng,n_species))
      ! The hydrogen budget is 1e10 cm^-3 (H2 holding two nuclei of the
      ! 0.5e10 side), and the oxygen and carbon ones are set far above what
      ! the cell below demands, so that the HYDROGEN row is the one that is
      ! tight and the rows are one and not three.
      budget = 0.5d0*1.0d10
      budget(:,ic_OH)  = 2.0d0*1.0d10
      budget(:,ic_H2O) = 2.0d0*1.0d10
      budget(:,ic_CO)  = 2.0d0*1.0d10
      call carrier_headroom_set_for_test(budget)
      Dsc = 1.0d0
      f_c = 0.0d0
      ! Cell 1's H2 sits exactly ON the cell's own density, which IS a
      ! coordinate face: no room above it at any step length. Cell 2's is
      ! well inside everything.
      do j = 1, N
         Yit(nvar_jac*(j-1)+1) = 1.0d0
         Yit(nvar_jac*(j-1)+2) = -7.0d0
         Yit(nvar_jac*(j-1)+3) = 1.0d0
         Yit(nvar_jac*(j-1)+4) = 0.1d0
         Yit(nvar_jac*(j-1)+5) = 0.0d0
         Yit(nvar_jac*(j-1)+6) = 0.0d0
         Yit(nvar_jac*(j-1)+7) = 0.0d0
         f_c(j,isp_H2) = 0.1d0
      enddo
      Yit(nvar_jac*0+4) = 1.0d0
      f_c(1,isp_H2)     = 1.0d0
      call freeze_species_unknown_box(Yit)
      ! A direction pointing out of that face at cell 1 and inward at cell
      ! 2, where there is room for a step of (1 - 0.1)/1 = 0.9.
      v = 0.0d0
      v(nvar_jac*0+4) = 1.0d0
      v(nvar_jac*1+4) = 1.0d0
      tmax = largest_step_inside_the_species_box(Yit, v, nblk)
      call int_row('one_component_with_no_room_is_counted', nblk, 1, nf)
      call real_row('the_direction_keeps_the_length_the_rest_allows',    &
                    tmax, 0.9d0, nf)
      vk = v
      call zero_the_blocked_components_of(Yit, vk)
      call real_row('the_blocked_component_is_zeroed_in_the_direction',  &
                    vk(nvar_jac*0+4), 0.0d0, nf)
      call real_row('every_other_component_keeps_its_length',            &
                    vk(nvar_jac*1+4), 1.0d0, nf)
      ! ---- the shared row: the projection, not the zero -------------
      ! Cell 2 carries the four oxygen-cycle balances, and OH and H2O both
      ! hold hydrogen, so they share one row. Put the cell one part in 1e9
      ! inside its hydrogen budget, which is inside the band that counts as
      ! active, so the row is active AND still has a little room: a
      ! direction pointing out of it is then cut to that room, and the same
      ! direction projected is not cut by it at all.
      slack = 1.0d-9
      Yit(nvar_jac*1+4) = 0.0d0
      Yit(nvar_jac*1+5) = 0.5d0*(1.0d0 - slack)
      Yit(nvar_jac*1+6) = 0.25d0*(1.0d0 - slack)
      Yit(nvar_jac*1+7) = 0.0d0
      f_c = 0.0d0
      f_c(1,isp_H2)  = 1.0d0
      f_c(2,isp_OH)  = 0.5d0*(1.0d0 - slack)
      f_c(2,isp_H2O) = 0.25d0*(1.0d0 - slack)
      call freeze_species_unknown_box(Yit)
      call freeze_element_constraint_rows(Yit, f_c, Dsc)
      ! n(OH) + 2 n(H2O) = (1 - slack) of the whole hydrogen budget.
      call bound_row('the_shared_hydrogen_row_is_nearly_tight',          &
                     abs(element_constraint_of_the_cell(2, 1)            &
                         /element_constraint_budget_of_the_cell(2, 1)    &
                         + slack), 1.0d-14, nf)
      ! The model gradient points out of the half-space, so the row is
      ! active: -g raises c, i.e. a . g < 0.
      g = 0.0d0
      g(nvar_jac*1+5) = -1.0d0
      g(nvar_jac*1+6) = -1.0d0
      call fix_active_species_bounds(Yit, g, nact)
      call int_row('the_shared_row_of_the_cell_is_active',              &
                   n_active_element_rows(), 1, nf)
      ! A direction that raises OH alone, which points straight out of the
      ! shared hydrogen row.
      v = 0.0d0
      v(nvar_jac*1+5) = 1.0d0
      len_in  = sqrt(sum(v*v))
      comp_in = v(nvar_jac*1+5)
      ! Unprojected, the row cuts it to the room the cell has left, which is
      ! the slack above: one part in 1e9 of the step.
      call bound_row('the_unprojected_direction_is_cut_at_the_row',      &
                largest_step_inside_the_element_constraints(Yit, v)      &
                /slack - 1.0d0, 1.0d-9, nf)
      vk = v
      call project_out_of_the_active_element_constraints(vk)
      len_out = sqrt(sum(vk*vk))
      ! What is removed is ONE COMBINATION of the carriers of the row, not
      ! the components: the three species that hold hydrogen here are H2 and
      ! H2O with two nuclei each and OH with one, so the row's direction is
      ! (2, 1, 2) over them, and a direction that raised OH alone comes back
      ! as 1 - 1/9 of that motion with the two others taking up the rest
      ! with the opposite sign. That is what moving ALONG a shared budget
      ! is, and it is what zeroing the blocked component threw away.
      call real_row('the_projected_direction_keeps_the_blocked_carrier', &
                    vk(nvar_jac*1+5), 1.0d0 - 1.0d0/9.0d0, nf)
      call real_row('the_projected_direction_moves_the_other_carriers',  &
                    vk(nvar_jac*1+6), -2.0d0/9.0d0, nf)
      call real_row('the_projected_direction_moves_the_H2_of_the_row',   &
                    vk(nvar_jac*1+4), -2.0d0/9.0d0, nf)
      call bound_row('the_projection_keeps_most_of_the_length',          &
                     len_in/max(len_out, 1.0d-300), 1.2d0, nf)
      ! And the combination that is removed is exactly the outward one: the
      ! projected direction no longer changes c to first order.
      call bound_row('the_projected_direction_moves_along_the_row',      &
                     abs(2.0d0*vk(nvar_jac*1+4)                          &
                         + 1.0d0*vk(nvar_jac*1+5)                        &
                         + 2.0d0*vk(nvar_jac*1+6))                       &
                     /max(abs(comp_in), 1.0d-300), 1.0d-15, nf)
      ! THE CELL'S OWN DENSITY IS BARELY MOVED, and that is the homogeneity
      ! of the constraint and not a choice: its column is c/nd, which is the
      ! slack itself here, so the density direction is tangent to the
      ! surface and an active row does not push the hydrodynamic mass
      ! unknown.
      call bound_row('the_active_row_does_not_push_the_cell_density',    &
                     abs(vk(nvar_jac*1+1))                               &
                     /max(abs(vk(nvar_jac*1+5)), 1.0d-300), 1.0d-8, nf)
      ! And it is no longer the hydrogen row that cuts the step: what is
      ! left is the loose oxygen row, nine decades further out.
      call bound_row('the_projected_direction_is_not_cut_by_the_row',    &
                slack                                                    &
                /largest_step_inside_the_element_constraints(Yit, vk),   &
                1.0d-6, nf)
      call release_active_species_bounds
      deallocate(Yit, v, vk, g, Dsc, budget, f_c)
      call set_transported_species_rows(.false.)
      N  = nsave;  n0 = nsave_n0
      end subroutine a_direction_along_a_shared_constraint_is_not_discarded

      ! ================================================================= !

      subroutine an_infeasible_start_is_restored_by_measured_violation(nf)
      ! A CANDIDATE THAT SPENDS NUCLEI THE CELL DOES NOT HAVE IS STEPPED
      ! BACK ONTO THE CONSTRAINT BEFORE ITS MERIT IS JUDGED (item N4b
      ! deliverable 3).
      !
      ! The measure is the map's own violation magnitude, relative to what
      ! the row allows, and the step is the one the marching path's limiter
      ! takes on the same constraint: the carriers holding the element are
      ! scaled by the one factor that puts their demand on the surface.
      !
      ! RED with EXHALE_ELEMENT_CONSTRAINT_ROWS=0, in which the budget is a
      ! coordinate face and no restoration exists: the violation is neither
      ! measured nor removed and no cell is moved.
      integer, intent(inout) :: nf
      real*8, allocatable :: Yv(:), Dsc(:), budget(:,:), f_c(:,:)
      real*8  :: nsave_n0, before, after
      integer :: nsave, j, ncell, nrow

      nsave = N;  nsave_n0 = n0
      N  = 4
      n0 = 1.0d10
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_He           = .false.
      he_diffusion         = .false.
      he_metal_diffusion   = .false.
      thereis_metals       = .false.
      call carrier_set_init()
      call read_species_unknown_space_controls
      call set_transported_species_rows(.true.)
      allocate(Yv(nvar_jac*N), Dsc(nvar_jac*N))
      allocate(budget(N,n_carrier_max), f_c(1-Ng:N+Ng,n_species))
      budget = 0.5d0*1.0d10
      call carrier_headroom_set_for_test(budget)
      Dsc = 1.0d0
      f_c = 0.0d0
      ! The iterate: every carrier at a quarter of its shared budget, so
      ! every row is strictly inside and what a row allows IS its budget.
      do j = 1, N
         Yv(nvar_jac*(j-1)+1) = 1.0d0
         Yv(nvar_jac*(j-1)+2) = -7.0d0
         Yv(nvar_jac*(j-1)+3) = 1.0d0
         Yv(nvar_jac*(j-1)+4) = 0.25d0
         f_c(j,isp_H2) = 0.25d0
      enddo
      call freeze_species_unknown_box(Yv)
      call freeze_element_constraint_rows(Yv, f_c, Dsc)
      ! A feasible candidate is not moved.
      call restore_the_element_budget(Yv, .true., ncell, nrow, before,    &
                                      after)
      call int_row('a_feasible_candidate_is_not_restored', ncell, 0, nf)
      call real_row('a_feasible_candidate_has_no_violation', before,      &
                    0.0d0, nf)
      ! An infeasible start: every carrier at twice the whole budget, so
      ! the demand is four times what the row allows and the violation is
      ! 3.0 of it.
      do j = 1, N
         Yv(nvar_jac*(j-1)+4) = 2.0d0
      enddo
      call restore_the_element_budget(Yv, .false., ncell, nrow, before,   &
                                      after)
      call real_row('restoration_measures_the_violation_before', before,  &
                    3.0d0, nf)
      call int_row('restoration_counts_the_rows_it_would_move', nrow, 4,  &
                   nf)
      call int_row('a_measured_candidate_is_not_moved', ncell, 0, nf)
      call restore_the_element_budget(Yv, .true., ncell, nrow, before,    &
                                     after)
      call int_row('restoration_moves_every_infeasible_cell', ncell, 4,   &
                   nf)
      call real_row('restoration_removes_the_violation', after, 0.0d0, nf)
      ! And the state it leaves is the one on the surface: the shared
      ! hydrogen budget is 1.0 in code density units, H2 holds two nuclei,
      ! so the carrier lands on 0.5.
      call real_row('restoration_lands_on_the_constraint_surface',        &
                    Yv(nvar_jac*0+4), 0.5d0, nf)
      deallocate(Yv, Dsc, budget, f_c)
      call set_transported_species_rows(.false.)
      N  = nsave;  n0 = nsave_n0
      end subroutine an_infeasible_start_is_restored_by_measured_violation

      ! ================================================================= !

      subroutine element_constraint_derivatives_against_differences(nf)
      ! THE DERIVATIVES OF THE CONSTRAINT MAP AGAINST FINITE DIFFERENCES OF
      ! IT (item N4b deliverable 1).
      !
      ! The constraint is affine in the carrier densities and in the frozen
      ! budget and homogeneous of degree one in the cell's own particle
      ! density, so the difference quotients below are EXACT and the rows
      ! carry no tolerance beyond the rounding of the two sums the value is.
      ! A row that needed a tolerance would be saying the map and its
      ! derivative are two expressions.
      integer, intent(inout) :: nf
      real*8  :: f_cell(n_species), budget(n_inv_element)
      real*8  :: c0(n_inv_element), c1(n_inv_element)
      real*8  :: dn0(n_species,n_inv_element), dn1(n_species,n_inv_element)
      real*8  :: dnd0(n_inv_element), dnd1(n_inv_element)
      real*8  :: db0(n_inv_element), db1(n_inv_element)
      logical :: carried(n_species)
      real*8  :: h, nd, q
      integer :: isp

      call carrier_species_that_hold(.true., carried)
      ! Densities in place of fractions with nd = 1, so the map's products
      ! are exact; the budget is the conditional one of a step.
      f_cell = 0.0d0
      f_cell(isp_H2)  = 0.125d0
      f_cell(isp_HII) = 0.25d0
      f_cell(isp_HI)  = 0.5d0
      f_cell(isp_OH)  = 0.0625d0
      f_cell(isp_H2O) = 0.03125d0
      f_cell(isp_CO)  = 0.015625d0
      budget = 0.0d0
      budget(ien_H) = 1.0d0
      budget(ien_O) = 0.5d0
      budget(ien_C) = 0.25d0
      nd = 1.0d0
      call element_constraint_derivatives(nd, f_cell, carried, budget,    &
                                          c0, dn0, dnd0, db0)
      ! ---- the species columns -------------------------------------
      ! Every species the map counts, perturbed one at a time.
      h = 0.0078125d0
      q = 0.0d0
      do isp = 1, n_species
         f_cell(isp) = f_cell(isp) + h
         call element_constraint_derivatives(nd, f_cell, carried, budget,  &
                                             c1, dn1, dnd1, db1)
         f_cell(isp) = f_cell(isp) - h
         q = max(q, abs((c1(ien_H) - c0(ien_H))/h - dn0(isp,ien_H)))
         q = max(q, abs((c1(ien_O) - c0(ien_O))/h - dn0(isp,ien_O)))
         q = max(q, abs((c1(ien_C) - c0(ien_C))/h - dn0(isp,ien_C)))
      enddo
      call real_row('constraint_species_derivative_is_the_difference',    &
                    q, 0.0d0, nf)
      ! The columns themselves are the stoichiometry: two hydrogen nuclei
      ! for H2, one for H+, one oxygen for OH and one carbon for CO.
      call real_row('constraint_column_of_H2_is_two_hydrogen_nuclei',     &
                    dn0(isp_H2,ien_H), 2.0d0, nf)
      call real_row('constraint_column_of_H2O_is_two_hydrogen_nuclei',    &
                    dn0(isp_H2O,ien_H), 2.0d0, nf)
      call real_row('constraint_column_of_H2O_is_one_oxygen_nucleus',     &
                    dn0(isp_H2O,ien_O), 1.0d0, nf)
      call real_row('constraint_column_of_CO_is_one_carbon_nucleus',      &
                    dn0(isp_CO,ien_C), 1.0d0, nf)
      ! A species the calculation does not carry is not in the demand.
      call real_row('constraint_column_of_a_species_not_carried_is_zero', &
                    dn0(isp_HI,ien_H), 0.0d0, nf)
      ! ---- the cell's own density ----------------------------------
      ! Homogeneous of degree one at fixed mass fractions and with the
      ! budget read as a fixed number of nuclei per unit mass, so the budget
      ! is scaled with the density and the quotient is exact.
      h = 0.0625d0
      call element_constraint_derivatives(nd + h, f_cell, carried,         &
                                    budget*(nd + h)/nd, c1, dn1, dnd1, db1)
      q = 0.0d0
      q = max(q, abs((c1(ien_H) - c0(ien_H))/h - dnd0(ien_H)))
      q = max(q, abs((c1(ien_O) - c0(ien_O))/h - dnd0(ien_O)))
      q = max(q, abs((c1(ien_C) - c0(ien_C))/h - dnd0(ien_C)))
      call real_row('constraint_density_derivative_is_the_difference',     &
                    q, 0.0d0, nf)
      call real_row('constraint_density_column_is_c_over_nd',              &
                    dnd0(ien_H), c0(ien_H)/nd, nf)
      ! ---- the frozen budget --------------------------------------
      h = 0.03125d0
      budget(ien_H) = budget(ien_H) + h
      call element_constraint_derivatives(nd, f_cell, carried, budget,     &
                                          c1, dn1, dnd1, db1)
      budget(ien_H) = budget(ien_H) - h
      call real_row('constraint_budget_derivative_is_the_difference',      &
                    (c1(ien_H) - c0(ien_H))/h, db0(ien_H), nf)
      call real_row('constraint_budget_column_is_minus_one',               &
                    db0(ien_H), -1.0d0, nf)
      ! ---- helium is held by no carrier ---------------------------
      call real_row('constraint_of_helium_is_identically_zero',            &
                    c0(ien_He), 0.0d0, nf)
      call real_row('constraint_budget_column_of_helium_is_zero',          &
                    db0(ien_He), 0.0d0, nf)
      end subroutine element_constraint_derivatives_against_differences

      ! ================================================================= !

      subroutine shared_element_constraint_is_the_feasible_set(nf)
      ! THE ADMISSIBLE SET OF THE SPECIES UNKNOWNS IS A SIMPLEX AND NOT A
      ! PRODUCT OF INTERVALS (item N4a; docs/ISSUES_20260909_review.md R4;
      ! docs/element_inventory_contexts_20260909.md).
      !
      ! A carrier holds nuclei of an element that other species of the same
      ! cell hold as well, so the constraint the carriers of one element
      ! share is
      !
      !     2 n(H2) + n(H+) + n(OH) + 2 n(H2O)  <=  available hydrogen
      !
      ! and the box side of each carrier is one face of a set that is
      ! larger than it.  The rows below evaluate the PRODUCTION map
      ! (element_inventory_of_candidate) on states the separate sides admit
      ! and the shared constraint does not, and they read the box sides
      ! from the production carrier_element_headroom, so nothing here is a
      ! second copy of either arithmetic.
      !
      ! The budget is the frozen data of the calculation: the element
      ! densities of the state the sides were formed at, with the carriers
      ! taken from the candidate.  That is what a trial of the coupled
      ! solve is, and it is why the constraint is not an identity.
      !
      ! RED, MEASURED with EXHALE_INVENTORY_SHARED=0, which judges by the
      ! separate sides: all three infeasible states below are called
      ! feasible.  GREEN with the shared constraint: they breach their
      ! budget by 1.0, by 0.25 and by 0.125 of it.
      integer, intent(inout) :: nf
      real*8  :: f_cell(n_species), free(n_inv_element)
      logical :: carried(n_species)
      type(element_inventory_cell) :: inv
      integer :: ie

      ! Densities in place of fractions, with nd = 1: the map's products are
      ! then exact and the numbers are the review's normalized ones.
      call carrier_species_that_hold(.true., carried)

      ! ---- both carriers exactly on their own ceilings -----------------
      ! One hydrogen nucleus is available.  n(H2) = 0.5 is its whole side
      ! and n(H+) = 1 is its whole side, and together they hold two.
      f_cell = 0.0d0
      free   = 0.0d0
      free(ien_H)     = 1.0d0
      f_cell(isp_H2)  = 0.5d0
      f_cell(isp_HII) = 1.0d0
      call element_inventory_of_candidate(1.0d0, f_cell, carried, free,   &
                                          inv)
      call real_row('inventory_two_nuclei_for_one',                       &
                    inv%n_carried(ien_H), 2.0d0, nf)
      call real_row('inventory_H2_is_on_its_own_side',                    &
                    f_cell(isp_H2),                                       &
                    carrier_element_headroom(ic_H2, free(ien_H), 0.0d0,   &
                                             0.0d0), nf)
      call real_row('inventory_Hp_is_on_its_own_side',                    &
                    f_cell(isp_HII),                                      &
                    carrier_element_headroom(ic_Hp, free(ien_H), 0.0d0,   &
                                             0.0d0), nf)
      call real_row('inventory_separate_ceilings_breach_is_one',          &
                    inv%violation(ien_H), 1.0d0, nf)
      call log_row('inventory_separate_ceilings_are_infeasible',          &
                   inventory_is_feasible(inv), .false., nf)
      call int_row('inventory_separate_ceilings_worst_element',           &
                   worst_element(inv), ien_H, nf)

      ! ---- every side satisfied, the hydrogen remainder negative -------
      ! Each carrier strictly inside its own side, and the closure species
      ! are left less than nothing.
      f_cell = 0.0d0
      free   = 0.0d0
      free(ien_H)     = 1.0d0
      f_cell(isp_H2)  = 0.375d0
      f_cell(isp_HII) = 0.5d0
      call element_inventory_of_candidate(1.0d0, f_cell, carried, free,   &
                                          inv)
      call bound_row('inventory_H2_is_inside_its_own_side',               &
                     f_cell(isp_H2)                                       &
                     - carrier_element_headroom(ic_H2, free(ien_H),       &
                                                0.0d0, 0.0d0), 0.0d0, nf)
      call bound_row('inventory_Hp_is_inside_its_own_side',               &
                     f_cell(isp_HII)                                      &
                     - carrier_element_headroom(ic_Hp, free(ien_H),       &
                                                0.0d0, 0.0d0), 0.0d0, nf)
      call bound_row('inventory_hydrogen_remainder_is_negative',          &
                     inv%remainder(ien_H), -1.0d-30, nf)
      call log_row('inventory_negative_remainder_is_infeasible',          &
                   inventory_is_feasible(inv), .false., nf)
      call real_row('inventory_negative_remainder_breach',                &
                    inv%violation(ien_H), 0.25d0, nf)

      ! ---- the oxygen carriers share one budget too --------------------
      ! OH, H2O and CO each inside the free oxygen, and together above it.
      f_cell = 0.0d0
      free   = 0.0d0
      free(ien_O)     = 1.0d0
      free(ien_C)     = 1.0d0
      f_cell(isp_OH)  = 0.375d0
      f_cell(isp_H2O) = 0.375d0
      f_cell(isp_CO)  = 0.375d0
      call element_inventory_of_candidate(1.0d0, f_cell, carried, free,   &
                                          inv)
      call bound_row('inventory_OH_is_inside_its_own_side',               &
                     f_cell(isp_OH)                                       &
                     - carrier_element_headroom(ic_OH, 0.0d0,             &
                                       free(ien_O), free(ien_C)),         &
                     0.0d0, nf)
      call bound_row('inventory_CO_is_inside_its_own_side',               &
                     f_cell(isp_CO)                                       &
                     - carrier_element_headroom(ic_CO, 0.0d0,             &
                                       free(ien_O), free(ien_C)),         &
                     0.0d0, nf)
      call real_row('inventory_oxygen_carriers_hold_more_than_the_element',&
                    inv%n_carried(ien_O), 1.125d0, nf)
      call log_row('inventory_shared_oxygen_is_infeasible',               &
                   inventory_is_feasible(inv), .false., nf)

      ! ---- a state inside the shared constraint is feasible ------------
      ! The control row: the same carriers scaled into the budget.
      f_cell = 0.0d0
      free   = 0.0d0
      free(ien_H)     = 1.0d0
      f_cell(isp_H2)  = 0.125d0
      f_cell(isp_HII) = 0.5d0
      f_cell(isp_HI)  = 0.25d0
      call element_inventory_of_candidate(1.0d0, f_cell, carried, free,   &
                                          inv)
      call log_row('inventory_inside_the_shared_constraint_is_feasible',  &
                   inventory_is_feasible(inv), .true., nf)
      do ie = 1, n_inv_element
         if (inv%violation(ie) .gt. 0.0d0) then
            call real_row('inventory_feasible_state_has_no_breach',       &
                          inv%violation(ie), 0.0d0, nf)
            exit
         endif
      enddo

      end subroutine shared_element_constraint_is_the_feasible_set

      ! ================================================================= !

      integer function worst_element(inv) result(ie)
      ! The element index whose breach is the worst of the cell, which the
      ! map returns beside the magnitude so that a refusal can name the
      ! element it rests on.
      type(element_inventory_cell), intent(in) :: inv
      real*8 :: w
      w = inventory_worst_violation(inv, ie)
      end function worst_element

      ! ================================================================= !

      subroutine map_and_operator_totals_are_one_number(nf)
      ! THE TWO ENTRY POINTS OF THE MAP RETURN THE SAME NUMBERS.
      !
      ! carrier_element_totals is what the carrier operator calls over the
      ! grid; element_inventory_of_cell is what a classification calls on
      ! one cell.  Both read the nucleus counts of the species table, and
      ! the row is exact: one stoichiometric source means one number, not
      ! two numbers within a tolerance.
      integer, intent(inout) :: nf
      real*8, allocatable :: nd(:), f_sp(:,:)
      real*8, allocatable :: nH_free(:), nO_free(:), nC_free(:)
      type(element_inventory_cell) :: inv
      logical :: carried(n_species)
      integer :: j, im
      real*8  :: dh, do_, dc

      N = 5
      allocate(nd(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(nH_free(1-Ng:N+Ng), nO_free(1-Ng:N+Ng), nC_free(1-Ng:N+Ng))
      nd   = 1.0d0
      f_sp = 0.0d0
      ! A state with every element of the map in more than one species.
      do j = 1-Ng, N+Ng
         f_sp(j,isp_HI)   = 1.0d0 + 0.1d0*dble(j)
         f_sp(j,isp_HII)  = 0.3d0
         f_sp(j,isp_H2)   = 0.7d0
         f_sp(j,isp_H2p)  = 0.05d0
         f_sp(j,isp_H3p)  = 0.02d0
         f_sp(j,isp_HeHp) = 0.01d0
         f_sp(j,isp_HeI)  = 0.9d0
         f_sp(j,isp_OH)   = 0.11d0
         f_sp(j,isp_H2O)  = 0.13d0
         f_sp(j,isp_CO)   = 0.17d0
         do im = 1, n_mion
            f_sp(j,mion_fsp(im)) = 0.01d0*dble(im)
         enddo
      enddo
      call carrier_element_totals(nd, f_sp, nH_free, nO_free, nC_free)
      call carrier_species_that_hold(.false., carried)
      dh = 0.0d0;  do_ = 0.0d0;  dc = 0.0d0
      do j = 1, N
         call element_inventory_of_cell(nd(j), f_sp(j,:), carried, inv)
         dh  = max(dh,  abs(inv%element_free(ien_H) - nH_free(j)))
         do_ = max(do_, abs(inv%element_free(ien_O) - nO_free(j)))
         dc  = max(dc,  abs(inv%element_free(ien_C) - nC_free(j)))
      enddo
      call real_row('inventory_hydrogen_total_is_the_operators', dh,      &
                    0.0d0, nf)
      call real_row('inventory_oxygen_total_is_the_operators',   do_,     &
                    0.0d0, nf)
      call real_row('inventory_carbon_total_is_the_operators',   dc,      &
                    0.0d0, nf)
      ! And the groups add up to the total, which is what makes the map an
      ! inventory rather than four independent sums.  The bound is four
      ! units in the last place of the total: the groups and the total are
      ! the same nuclei added in two orders, so what separates them is the
      ! rounding of those orders and nothing else.
      call element_inventory_of_cell(nd(1), f_sp(1,:), carried, inv)
      call bound_row('inventory_groups_close_the_hydrogen_total',         &
                     abs(inv%n_carried(ien_H) + inv%n_locked(ien_H)       &
                         + inv%n_reserved(ien_H) + inv%n_closure(ien_H)   &
                         - inv%n_total(ien_H)),                           &
                     4.0d0*epsilon(1.0d0)*inv%n_total(ien_H), nf)
      call bound_row('inventory_groups_close_the_oxygen_total',           &
                     abs(inv%n_carried(ien_O) + inv%n_locked(ien_O)       &
                         + inv%n_reserved(ien_O) + inv%n_closure(ien_O)   &
                         - inv%n_total(ien_O)),                           &
                     4.0d0*epsilon(1.0d0)*inv%n_total(ien_O), nf)
      call real_row('inventory_helium_of_the_hydride_is_locked',          &
                    inv%n_locked(ien_He), f_sp(1,isp_HeHp), nf)
      deallocate(nd, f_sp, nH_free, nO_free, nC_free)
      end subroutine map_and_operator_totals_are_one_number

      ! ================================================================= !

      subroutine reference_states_are_classified_under_the_map(nf)
      ! EVERY SELECTED REFERENCE STATE EVALUATED AND CLASSIFIED UNDER THE
      ! MAP (item N4a acceptance; docs/PLAN_20260909_review.md F7).
      !
      ! The reference states are the golden Ion_species.txt of the
      ! molecular regression cases, which carry the species densities of
      ! every cell in cm^-3, so the map is called on them with nd = 1 and
      ! nothing is converted.  The carrier set used is the widest one a
      ! configuration can transport (H2, OH, H2O, CO), which puts the most
      ! nuclei on the carriers' side of the budget and is therefore the
      ! strongest statement the map can make about a state.
      !
      ! WHAT THE ROW STATES AND WHAT IT DOES NOT.  Evaluated at a state
      ! that came from its own element closure the hydrogen constraint
      ! reduces to n(H I) >= 0, so a breach here is a departure of that
      ! state from its own inventory and not a statement about the step
      ! that produced it.  The measured worst breach and the tightest
      ! occupancy of each budget are printed for every case; a golden that
      ! breached the map would expose a defect of that state, never a
      ! reason to weaken the map.
      integer, intent(inout) :: nf
      integer, parameter :: ncase = 9
      character(len=20) :: cases(ncase)
      ! Which of them TRANSPORTS helium, READ from "He_diffusion: True" in
      ! the case's own input.inp: with the helium mass fraction a
      ! transported unknown a He/H gradient is the solution, so the ratio of
      ! such a case is measured and not gated.
      logical :: helium_transported(ncase)
      character(len=256) :: root, path
      real*8  :: worst, occ, dhe
      integer :: ic, jworst, ie, jhe
      logical :: ok
      ! The element census's own gate on an exact bookkeeping identity
      ! (element_census.f90): three decades above the accumulated round-off
      ! of a marching state and a decade below the cell solver's xtol.  It
      ! is not a percentage of abundance.
      ! The gate on an element ratio, from the map: element_census's own
      ! number for an exact bookkeeping identity, not a bound chosen here.
      real*8, parameter :: census_gate = element_ratio_gate
      cases = [ character(len=20) ::                                      &
                'mol_base_handoff', 'mol_metals', 'mol_lyman_werner',     &
                'mol_diffusion', 'mol_ir_bands', 'mol_sec_ion',           &
                'mol_carrier', 'lower_profile', 'oxygen_chemistry' ]
      helium_transported = [ .false., .false., .false.,                   &
                             .true.,  .false., .false.,                   &
                             .false., .true.,  .false. ]
      call get_environment_variable('EXHALE_INVENTORY_GOLDEN_DIR', root)
      if (len_trim(root) .eq. 0) root = 'backup/regression/golden'
      do ic = 1, ncase
         path = trim(root)//'/'//trim(cases(ic))//'/Ion_species.txt'
         call classify_one_reference_state(path, worst, occ, jworst, ie,   &
                                           dhe, jhe, ok)
         if (.not. ok) then
            write(*,'(A,A,A,A,A)') 'FAIL inventory_reference_',           &
               trim(cases(ic)), ' measured=unreadable reference=',        &
               trim(path), ' tol=0'
            nf = nf + 1
            cycle
         endif
         write(*,'(A,A,A,ES12.5,A,I0,A,A,A,ES12.5)')                      &
            '  inventory_reference_', trim(cases(ic)),                    &
            ': worst breach ', worst, ' at cell ', jworst,                &
            ' in ', trim(inv_element_name(max(ie,1))),                    &
            ', tightest occupancy ', occ
         write(*,'(A,A,A,ES12.5,A,I0)')                                    &
            '  inventory_reference_', trim(cases(ic)),                    &
            ': He/H departure from the base reservoir ', dhe,             &
            ' at cell ', jhe
         call bound_row('inventory_reference_'//trim(cases(ic))//          &
                        '_is_feasible', worst,                            &
                        inventory_feasible_slack, nf)
         ! AND THE STATE STANDS IN THE ELEMENT RATIO ITS RUN WAS GIVEN.
         ! No operator of a run without element transport may move the
         ! ratio of helium to hydrogen nuclei: the ionization sweep
         ! repartitions stages, the carrier step moves nuclei inside ONE
         ! element, and the write-back restores the element totals it was
         ! handed.  The two cases that transport helium are measured
         ! instead, their gradient being the transported solution.
         if (.not. helium_transported(ic))                                &
            call bound_row('inventory_reference_'//trim(cases(ic))//       &
                           '_keeps_its_He_over_H', dhe, census_gate, nf)
      enddo
      call coupled_carrier_states_are_classified(nf)
      end subroutine reference_states_are_classified_under_the_map

      ! ================================================================= !

      subroutine coupled_carrier_states_are_classified(nf)
      ! THE STATES OF A COUPLED CARRIER SOLVE, CLASSIFIED UNDER THE SAME MAP.
      !
      ! EXHALE_INVENTORY_COUPLED_STATES is a colon-separated list of run
      ! directories whose output/Ion_species.txt is a state a JFNK with
      ! carrier rows returned.  Without it the rows are skipped rather than
      ! asserted against a state that is not there.
      !
      ! MEASURED (2026-09-09, item N9): every coupled-carrier JFNK state of
      ! items B5j and B5k departs from its input He/H by 4 to 8 percent with
      ! the carrier unknown a density and by 41 to 43 percent with it in
      ! ln n, while 12000 marching steps of the same case keep the ratio to
      ! 1e-11.  READ in the source: the carrier column is written AFTER the
      ! element projection, so the Newton's n(H2) is not rescaled onto an
      ! element total, and with helium diffusion off nothing re-imposes one.
      !
      ! The row below states that the map SEES that departure.  It is a
      ! defect of those states and not a property to keep: when the carrier
      ! row restores the element total, the state will stand in its
      ! reservoir's ratio and this row is replaced by the conserving one.
      integer, intent(inout) :: nf
      character(len=1024) :: list
      character(len=256)  :: dir, path
      real*8  :: worst, occ, dhe
      integer :: jworst, ie, jhe, i1, i2, n
      logical :: ok
      call get_environment_variable('EXHALE_INVENTORY_COUPLED_STATES',    &
                                    list)
      if (len_trim(list) .eq. 0) then
         write(*,'(A)') '  coupled carrier state rows skipped: set'//     &
            ' EXHALE_INVENTORY_COUPLED_STATES to the run'
         write(*,'(A)') '  directories of a JFNK solve that carried'//    &
            ' carrier rows'
         return
      endif
      n  = len_trim(list)
      i1 = 1
      do while (i1 .le. n)
         i2 = index(list(i1:n), ':')
         if (i2 .eq. 0) then
            dir = list(i1:n)
            i1  = n + 1
         else
            dir = list(i1:i1+i2-2)
            i1  = i1 + i2
         endif
         if (len_trim(dir) .eq. 0) cycle
         path = trim(dir)//'/Ion_species.txt'
         call classify_one_reference_state(path, worst, occ, jworst, ie,   &
                                           dhe, jhe, ok)
         if (.not. ok) then
            write(*,'(A,A,A)') 'FAIL inventory_coupled_state'//           &
               ' measured=unreadable reference=', trim(path), ' tol=0'
            nf = nf + 1
            cycle
         endif
         write(*,'(A,A,A,ES12.5,A,ES12.5,A,I0)')                          &
            '  inventory_coupled_state ', trim(dir),                      &
            ': worst carrier breach ', worst,                             &
            ', He/H departure ', dhe, ' at cell ', jhe
         ! The carrier budgets are satisfied by such a state -- each
         ! carrier is inside the element it is charged to -- and the
         ! element ratio is not, which is why the two are separate
         ! statements of the map and not one number.
         call bound_row('inventory_coupled_state_carrier_budget_holds',   &
                        worst, inventory_feasible_slack, nf)
         if (dhe .gt. element_ratio_gate) then
            write(*,'(A,ES12.5,A,ES12.5,A)')                              &
               'PASS inventory_coupled_state_He_over_H_is_seen'//         &
               ' measured=', dhe, ' reference=>', element_ratio_gate,    &
               ' tol=0'
         else
            write(*,'(A,ES12.5,A,ES12.5,A)')                              &
               'FAIL inventory_coupled_state_He_over_H_is_seen'//         &
               ' measured=', dhe, ' reference=>', element_ratio_gate,    &
               ' tol=0'
            nf = nf + 1
         endif
      enddo
      end subroutine coupled_carrier_states_are_classified

      ! ================================================================= !

      subroutine classify_one_reference_state(path, worst, occ, jworst,    &
                                              ie_worst, dhe, jhe, ok)
      ! One reference state read from its own schema header and classified
      ! cell by cell under the map.  The header names the species of each
      ! column and the rows the physical cells occupy, so the reader
      ! follows the file rather than a column count written here.
      character(len=*), intent(in)  :: path
      real*8,           intent(out) :: worst, occ
      integer,          intent(out) :: jworst, ie_worst
      ! The departure of the He/H nucleus ratio from the base cell's, which
      ! is the reservoir the run pins, and the cell that carries it.
      real*8,           intent(out) :: dhe
      integer,          intent(out) :: jhe
      logical,          intent(out) :: ok
      integer, parameter :: mxcol = 64
      character(len=8192) :: line
      character(len=32)   :: tok(mxcol)
      integer :: colsp(mxcol)
      real*8  :: val(mxcol)
      real*8  :: f_cell(n_species), v
      real*8, allocatable :: yhe(:)
      logical :: carried(n_species)
      type(element_inventory_cell) :: inv
      integer :: u, ios, ncol, ntok, i, j, nrow, jcell, ie
      ok       = .false.
      worst    = 0.0d0
      occ      = 0.0d0
      jworst   = 0
      ie_worst = 0
      dhe      = 0.0d0
      jhe      = 0
      call carrier_species_that_hold(.false., carried)
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios .ne. 0) return
      ncol = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .ne. '#') exit
         if (index(line, '# columns') .eq. 1) then
            call split_tokens(line, tok, ntok, mxcol)
            ! tok(1) = '#', tok(2) = 'columns', tok(3) = 'r[Rp]'
            ncol = ntok - 2
            do i = 1, ncol
               colsp(i) = species_column_of_name(trim(tok(i+2)))
            enddo
         endif
      enddo
      if (ncol .le. 1) then
         close(u)
         return
      endif
      ! The first data line has already been read; count the rest.
      nrow = 1
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0) cycle
         nrow = nrow + 1
      enddo
      rewind(u)
      allocate(yhe(max(nrow - 2*Ng, 1)))
      yhe = 0.0d0
      j = 0
      jcell = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#') cycle
         if (len_trim(line) .eq. 0) cycle
         j = j + 1
         ! Ng ghost cells at each end, as the schema header states; the
         ! physical cells are the ones a constraint is posed on.
         if (j .le. Ng .or. j .gt. nrow - Ng) cycle
         jcell = jcell + 1
         read(line,*,iostat=ios) val(1:ncol)
         if (ios .ne. 0) then
            close(u)
            return
         endif
         f_cell = 0.0d0
         do i = 2, ncol
            if (colsp(i) .gt. 0) f_cell(colsp(i)) = val(i)
         enddo
         call element_inventory_of_cell(1.0d0, f_cell, carried, inv)
         yhe(jcell) = helium_to_hydrogen(inv)
         v = inventory_worst_violation(inv, ie)
         if (v .gt. worst) then
            worst    = v
            jworst   = jcell
            ie_worst = ie
         endif
         do ie = 1, n_inv_element
            if (inv%occupancy(ie) .lt. huge(1.0d0))                       &
               occ = max(occ, inv%occupancy(ie))
         enddo
      enddo
      close(u)
      call ratio_departure_from_the_base(yhe, jcell, dhe, jhe)
      deallocate(yhe)
      ok = (jcell .gt. 0)
      end subroutine classify_one_reference_state

      ! ================================================================= !

      integer function species_column_of_name(name) result(isp)
      ! The f_sp column of an output-schema species name, from the same
      ! label tables the writer used.  Zero for a column that is not a
      ! species density.
      character(len=*), intent(in) :: name
      integer :: im
      isp = 0
      select case (trim(name))
      case ('HI');    isp = isp_HI
      case ('HII');   isp = isp_HII
      case ('HeI');   isp = isp_HeI
      case ('HeII');  isp = isp_HeII
      case ('HeIII'); isp = isp_HeIII
      case ('HeITR'); isp = isp_HeTR
      case ('H2');    isp = isp_H2
      case ('H2p');   isp = isp_H2p
      case ('H3p');   isp = isp_H3p
      case ('HeHp');  isp = isp_HeHp
      case ('OH');    isp = isp_OH
      case ('H2O');   isp = isp_H2O
      case ('CO');    isp = isp_CO
      case default
         do im = 1, n_mion
            if (trim(name) .eq. trim(mion_name(im))) then
               isp = mion_fsp(im)
               return
            endif
         enddo
      end select
      end function species_column_of_name

      ! ================================================================= !

      subroutine split_tokens(line, tok, ntok, mx)
      character(len=*), intent(in)  :: line
      character(len=32), intent(out) :: tok(*)
      integer,          intent(out) :: ntok
      integer,          intent(in)  :: mx
      integer :: i, i1, n
      ntok = 0
      i    = 1
      n    = len_trim(line)
      do while (i .le. n)
         do while (i .le. n)
            if (line(i:i) .ne. ' ') exit
            i = i + 1
         enddo
         if (i .gt. n) exit
         i1 = i
         do while (i .le. n)
            if (line(i:i) .eq. ' ') exit
            i = i + 1
         enddo
         if (ntok .ge. mx) return
         ntok = ntok + 1
         tok(ntok) = line(i1:min(i-1, i1+31))
      enddo
      end subroutine split_tokens

      ! ================================================================= !

      real*8 function reservoir_helium_mass_fraction(m_1)
      ! The reservoir helium mass fraction of the input He/H, the value the
      ! element operator pins its base to.
      real*8, intent(in) :: m_1
      reservoir_helium_mass_fraction = m_He_over_m_H*HeH/(m_1 + m_He_over_m_H*HeH)
      end function reservoir_helium_mass_fraction

      ! ================================================================= !


      ! ================================================================= !

      subroutine the_upper_ghost_of_a_transported_element_is_the_iterate(nf)
      ! THE UPPER GHOST OF A TRANSPORTED ELEMENT IS THE OUTERMOST CELL'S,
      ! AND IT IS THE ITERATE'S (item N26).
      !
      ! The element operator poses a zero-gradient outflow condition at the
      ! top of the column: nothing enters the domain from outside through an
      ! element, so the ghost that closes the outer face of cell N carries
      ! the mixing ratio cell N carries. The composition write-back of the
      ! stationary solve writes the physical cells from the unknown vector,
      ! and unless it writes that ghost too the outer face is reconstructed
      ! against a number the Newton cannot move. MEASURED on the atomic
      ! element reload (item N25): the sodium mixing ratio of cell N+1 stood
      ! at the pre-solve reservoir value 1.73111e-06 through the whole solve
      ! while cell N fell to 1.04e-11, and the row of cell N then measured
      ! the flux that value produces, 5.4e-02 of the row's own scale and the
      ! whole judged distance of the solve.
      !
      ! Three rows here. After the write-back the upper ghost equals the
      ! outermost physical cell, in the element MASS fraction the operator's
      ! projection carries and in the element MIXING RATIO the trace row
      ! reads. And a column with no interior gradient leaves the row of the
      ! outermost cell at zero, which it cannot do while the ghost supplies a
      ! face flux; the last pair of rows puts the stale ghost back and shows
      ! the row is not zero by construction.
      !
      ! The reconstruction is WENO3, the scheme the stationary solve runs on.
      ! Gravity is off, so the settling drift vanishes and a flat mixing
      ! ratio is an exact stationary state of the operator: what is left in
      ! the outermost row is the boundary and nothing else.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho, Y
      real*8, dimension(:),   allocatable :: res_he, sc_he
      real*8, dimension(:,:), allocatable :: res_tr, sc_tr, f_c, Yel, u
      real*8  :: m_1, mpH, dr_u, famp, dghost, dmix, row_flat, row_stale
      real*8  :: fmix_ghost, fmix_cell
      logical :: ok_he, ok_tr
      integer :: j, ie, im

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      b0 = 0.0d0
      spherical_domain     = .true.
      thereis_He           = .true.
      thereis_HeITR        = .false.
      thereis_mol          = .false.
      carrier_transport    = .false.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_metals       = .true.
      eos_include_metals   = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .true.
      he_ambipolar         = .false.
      he_alphaT            = 0.0d0
      he_kzz               = 1.0d11
      HeH                  = 0.0833333333333333d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
      recon_lambda_on = .false.
      call carrier_set_init()
      call set_transported_species_rows(.true.)

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr_u = 1.0d0/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4

      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), u(3,1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), Yel(1-Ng:N+Ng,1+n_melem))
      allocate(res_he(1:N), sc_he(1:N), res_tr(1:N,n_melem),              &
               sc_tr(1:N,n_melem), Y(nvar_jac*N))

      ! r^2 F_rho constant, the wind a steady mass row carries: with a flat
      ! composition the material divergence is then exactly zero cell by
      ! cell, and the balance the row measures is empty.
      famp = 1.0d-4
      m_1  = mass_per_H_nucleus_without_He()
      mpH  = m_1 + m_He_over_m_H*HeH
      f_c  = 0.0d0
      do j = 1-Ng, N+Ng
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = HeH/mpH
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = melem_ab(ie)*f_c(j,isp_HI)
         enddo
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo
      ! THE GHOST THE SOLVE WAS HANDED, five decades from the cell it
      ! continues, which is the distance the reload measured.
      do j = N+1, N+Ng
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = 3.0d0*HeH/mpH
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = 1.0d5*melem_ab(ie)/mpH
         enddo
      enddo
      u = 0.0d0
      u(1,:) = rho_c

      Y = 0.0d0
      call pack_species_rows(u, f_c, Y)
      call write_species_rows_into_composition(Y, u, f_c)

      call element_mass_fractions(f_c, Yel)
      dghost = 0.0d0
      do im = 1, 1 + n_melem
         do j = N+1, N+Ng
            dghost = max(dghost, abs(Yel(j,im) - Yel(N,im))               &
                                 /max(abs(Yel(N,im)), 1.0d-300))
         enddo
      enddo
      ! The bound is eight units in the last place: the ghost's fraction is
      ! written through the mass sum of the projection and read back through
      ! the same sum, so what separates it from the cell it copies is the
      ! rounding of that round trip and nothing else.
      call bound_row('upper_ghost_element_mass_fraction_is_the_cell',     &
                     dghost, 8.0d0*epsilon(1.0d0), nf)

      ! The mixing ratio the trace row reads, n_X/n_H: hydrogen sits in H I
      ! alone in this column, so it is one ratio of two entries of f_sp.
      dmix = 0.0d0
      do ie = 1, n_melem
         fmix_cell = f_c(N,mion_fsp(melem_i0(ie)))/f_c(N,isp_HI)
         do j = N+1, N+Ng
            fmix_ghost = f_c(j,mion_fsp(melem_i0(ie)))/f_c(j,isp_HI)
            dmix = max(dmix, abs(fmix_ghost - fmix_cell)/melem_ab(ie))
         enddo
      enddo
      call bound_row('upper_ghost_element_mixing_ratio_is_the_cell',      &
                     dmix, 1.0d-13, nf)

      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he,      &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      row_flat = abs(res_he(N))/max(sc_he(N), 1.0d-300)
      do im = 1, n_melem
         row_flat = max(row_flat, abs(res_tr(N,im))                       &
                                  /max(sc_tr(N,im), 1.0d-300))
      enddo
      call bound_row('outermost_element_row_of_a_flat_column_is_zero',    &
                     row_flat, 1.0d-12, nf)

      ! WHY THE WRITE-BACK HAS TO SET THE GHOST: the element operator writes
      ! no ghost of its own, so the outer boundary of the balance is the
      ! composition ghost the caller handed over.  The same interior with a
      ! ghost five decades away from its cell makes the outermost row read a
      ! face flux the column does not carry, and that is what the three rows
      ! above are worth: they state that the write-back sets the ghost the
      ! row is then closed with.
      do j = N+1, N+Ng
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = 3.0d0*HeH/mpH
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = 1.0d5*melem_ab(ie)/mpH
         enddo
      enddo
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he,      &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      row_stale = abs(res_he(N))/max(sc_he(N), 1.0d-300)
      do im = 1, n_melem
         row_stale = max(row_stale, abs(res_tr(N,im))                     &
                                    /max(sc_tr(N,im), 1.0d-300))
      enddo
      if (row_stale .gt. 1.0d-3) then
         write(*,'(A,A,A,ES12.5,A)') 'PASS ',                             &
              'a_stale_upper_ghost_is_read_by_the_outermost_row',         &
              ' measured=', row_stale, ' reference=>1e-3 tol=0'
      else
         write(*,'(A,A,A,ES12.5,A)') 'FAIL ',                             &
              'a_stale_upper_ghost_is_read_by_the_outermost_row',         &
              ' measured=', row_stale, ' reference=>1e-3 tol=0'
         nf = nf + 1
      endif

      deallocate(rho_c, v_c, T_c, Frho, u, f_c, Yel, res_he, sc_he,       &
                 res_tr, sc_tr, Y)
      end subroutine the_upper_ghost_of_a_transported_element_is_the_iterate

      ! ================================================================= !

      subroutine the_upper_ghost_of_a_transported_carrier_is_the_iterate(nf)
      ! THE SAME CONDITION FOR A CARRIER, AND THE ELEMENT IT SITS IN (N26).
      !
      ! A carrier is a partition WITHIN an element. Its upper ghost is the
      ! outermost cell's fraction, the zero-gradient outflow condition the
      ! carrier operator writes at its own ghosts before it measures its
      ! residual; and writing it may not move the element total the ghost
      ! holds, which the element condition has just set.
      !
      ! The third row is the three-unknown route: with no species row
      ! registered the write-back writes nothing at all, ghosts included.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 12
      real*8, dimension(:),   allocatable :: rho_c, Y, Y3
      real*8, dimension(:,:), allocatable :: f_c, f_keep, Yel, u
      real*8  :: m_1, mpH, dr_u, dcar, dhe, dnone
      integer :: j

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      b0 = 0.0d0
      spherical_domain     = .true.
      thereis_He           = .true.
      thereis_HeITR        = .false.
      thereis_mol          = .true.
      carrier_transport    = .true.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_metals       = .false.
      he_diffusion         = .true.
      he_metal_diffusion   = .false.
      HeH                  = 0.0833333333333333d0
      call carrier_set_init()
      call set_transported_species_rows(.true.)

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr_u = 1.0d0/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u

      allocate(rho_c(1-Ng:N+Ng), u(3,1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), f_keep(1-Ng:N+Ng,n_species),     &
               Yel(1-Ng:N+Ng,1+n_melem))
      allocate(Y(nvar_jac*N))

      m_1  = mass_per_H_nucleus_without_He()
      mpH  = m_1 + m_He_over_m_H*HeH
      f_c  = 0.0d0
      do j = 1-Ng, N+Ng
         ! A quarter of the hydrogen nuclei in H2, and a column that is not
         ! flat, so the ghost is not the cell by accident.
         f_c(j,isp_H2)  = 0.125d0/mpH*(1.0d0 + 0.02d0*dble(j))
         f_c(j,isp_HI)  = 1.0d0/mpH - 2.0d0*f_c(j,isp_H2)
         f_c(j,isp_HeI) = HeH/mpH
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
      enddo
      ! The ghost the solve was handed: a partition of its own.
      do j = N+1, N+Ng
         f_c(j,isp_H2) = 0.2d0*f_c(j,isp_H2)
         f_c(j,isp_HI) = 1.0d0/mpH - 2.0d0*f_c(j,isp_H2)
      enddo
      u = 0.0d0
      u(1,:) = rho_c

      Y = 0.0d0
      call pack_species_rows(u, f_c, Y)
      call write_species_rows_into_composition(Y, u, f_c)

      dcar = 0.0d0
      do j = N+1, N+Ng
         dcar = max(dcar, abs(f_c(j,isp_H2) - f_c(N,isp_H2))              &
                          /max(abs(f_c(N,isp_H2)), 1.0d-300))
      enddo
      call real_row('upper_ghost_carrier_fraction_is_the_cell', dcar,     &
                    0.0d0, nf)

      call element_mass_fractions(f_c, Yel)
      dhe = 0.0d0
      do j = N+1, N+Ng
         dhe = max(dhe, abs(Yel(j,1) - Yel(N,1))                          &
                        /max(abs(Yel(N,1)), 1.0d-300))
      enddo
      call bound_row('upper_ghost_helium_survives_the_carrier_copy', dhe, &
                     1.0d-14, nf)

      ! The three-unknown route: no species row, nothing written.
      call set_transported_species_rows(.false.)
      f_keep = f_c
      allocate(Y3(nvar_jac*N))
      Y3 = 1.0d0
      call write_species_rows_into_composition(Y3, u, f_c)
      dnone = maxval(abs(f_c - f_keep))
      call real_row('no_species_row_writes_no_composition', dnone, 0.0d0, &
                    nf)

      deallocate(rho_c, u, f_c, f_keep, Yel, Y, Y3)
      end subroutine the_upper_ghost_of_a_transported_carrier_is_the_iterate

      ! ================================================================= !

      subroutine bound_row(name, got, tol, nf)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: got, tol
      integer,          intent(inout) :: nf
      if (got .le. tol) then
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'PASS ', name,               &
              ' measured=', got, ' reference=0 tol=', tol
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'FAIL ', name,               &
              ' measured=', got, ' reference=0 tol=', tol
         nf = nf + 1
      endif
      end subroutine bound_row

      subroutine int_row(name, got, want, nf)
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: got, want
      integer,          intent(inout) :: nf
      if (got .eq. want) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=', got,  &
              ' reference=', want, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=', got,  &
              ' reference=', want, ' tol=0'
         nf = nf + 1
      endif
      end subroutine int_row

      subroutine real_row(name, got, want, nf)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: got, want
      integer,          intent(inout) :: nf
      if (got .eq. want) then
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'PASS ', name,             &
              ' measured=', got, ' reference=', want, ' tol=0'
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'FAIL ', name,             &
              ' measured=', got, ' reference=', want, ' tol=0'
         nf = nf + 1
      endif
      end subroutine real_row

      subroutine log_row(name, got, want, nf)
      character(len=*), intent(in)    :: name
      logical,          intent(in)    :: got, want
      integer,          intent(inout) :: nf
      if (got .eqv. want) then
         write(*,'(A,A,A,L1,A,L1,A)') 'PASS ', name, ' measured=', got,  &
              ' reference=', want, ' tol=0'
      else
         write(*,'(A,A,A,L1,A,L1,A)') 'FAIL ', name, ' measured=', got,  &
              ' reference=', want, ' tol=0'
         nf = nf + 1
      endif
      end subroutine log_row

      end program steady_species_rows_tests

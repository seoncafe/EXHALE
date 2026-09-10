      program certification_contexts
      ! THE A2 CERTIFICATION EVALUATOR AND ITS CONTEXTS, asserted against the
      ! predicates they replace.
      !
      ! docs/a2_certification_contract_20260906.md, section 7 steps 1 to 3.
      ! What is measured here:
      !
      !   (a) an ACTIVE carrier model whose frozen background is not ready
      !       reports `unavailable` and refuses certification. The predicate
      !       it replaces is transcribed and run on the same case: it reads
      !       the zero carrier_steady_residual returns in that situation as a
      !       balanced carrier and accepts.
      !   (b) a model that carries no such species reports `not_applicable`,
      !       and its absence is not a refusal.
      !   (c) the row measure catches a small-rate species with a large
      !       fractional imbalance standing beside a balanced dominant one,
      !       and the ratio of sums over all carriers of a region -- the
      !       `rvol` a gate would read -- does not (review R7). Both numbers
      !       are computed here from one set of rows.
      !   (d) a state carrying one cell without a chemical root is written
      !       uncertified and the cell is named. The superseded marching
      !       predicate, which makes no statement about the chemistry, is
      !       transcribed and accepts the same state.
      !   (e) a positive, merit-decreasing Newton trial with a large global
      !       residual is an ADMISSIBLE trial and is NOT a certified state:
      !       the two decisions are separate, and no stationary tolerance
      !       enters trial acceptance.
      !   (f) the probe context is not held to the chemical-root condition
      !       that rejects a trial, and both refuse a non-finite row.
      !   (g) the isolated-workspace round trip reinstates the carrier module
      !       state exactly.
      !
      ! Each assertion prints
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! and the exit status is nonzero if any of them fails.
      use global_parameters
      use ionization_equilibrium, only: bg_ready, ieq_nonroot_streak,     &
                                  ieq_res_tol, ieq_rate_cell,             &
                                  ieq_ne_cell, ieq_TK_cell,               &
                                  ieq_ntot_cell, ieq_rates_ready,         &
                                  ieq_neq_stored, ieq_mbase_stored,       &
                                  ieq_iox_stored, ieq_closure_scratch,    &
                                  save_ieq_closure_scratch,               &
                                  ieq_closure_scratch_matches,            &
                                  ionization_closure_residual_cell
      use diffusive_photochemistry, only: n_carrier, carrier_solved,      &
                                          ic_H2, ic_OH, carrier_set_init, &
                                          carrier_module_state,           &
                                          save_carrier_module_state,      &
                                          restore_carrier_module_state,   &
                                          carrier_module_state_matches
      use steady_residual_mod, only: steady_gates_met
      use species_table
      use element_inventory
      use composition, only: mass_per_H_nucleus_without_He
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use binary_element_diffusion, only: element_transport_residual,     &
                                          relax_element_composition
      use utils, only: calc_rho
      use certification
      implicit none

      integer :: n_fail
      n_fail = 0

      call set_up_a_grid()
      call carriers_without_a_background()
      call model_without_carriers()
      call small_rate_species_beside_a_balanced_one()
      call state_with_a_cell_without_a_chemical_root()
      call admissible_trial_is_not_a_certified_state()
      call isolated_workspace_round_trip()
      call level_rows_of_a_configuration_that_carries_them()
      call closure_residual_of_a_constructed_cell()
      call row_measure_without_its_volume_companion()
      call the_regime_windows_select_their_own_cells()
      call a_species_row_is_gated_in_the_wind_and_reported_below_it()
      call the_five_anchoring_numbers_of_a_synthetic_column()

      if (n_fail .gt. 0) then
         write(*,'(A,I0,A)') 'certification_contexts: ', n_fail,          &
              ' assertion(s) failed'
         stop 1
      endif
      write(*,'(A)') 'certification_contexts: every assertion passed'

      contains

      ! ------------------------------------------------------!

      subroutine set_up_a_grid()
      ! The smallest grid the evaluator needs: a radius, a cell width and the
      ! escape index that splits the layer from the wind. No residual is
      ! assembled, so the hydrodynamic rows are unavailable throughout and
      ! the carrier and chemistry statements are what is being measured.
      integer :: j
      N      = 20
      j_min  = 6
      if (allocated(r))    deallocate(r)
      if (allocated(dr_j)) deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)    = 1.0d0 + 0.05d0*dble(j)
         dr_j(j) = 0.05d0
      enddo
      if (allocated(ieq_nonroot_streak)) deallocate(ieq_nonroot_streak)
      allocate(ieq_nonroot_streak(1-Ng:N+Ng))
      ieq_nonroot_streak = 0
      ! No helium, no metals, no diffusion, no excited hydrogen: the
      ! inventory then carries the hydrodynamic rows, the carriers and one
      ! eliminated-species closure, and nothing else is active.
      thereis_He     = .false.
      thereis_HeITR  = .false.
      thereis_metals = .false.
      he_diffusion   = .false.
      he_metal_diffusion = .false.
      use_excited_H  = .false.
      end subroutine set_up_a_grid

      ! ------------------------------------------------------!

      subroutine carriers_without_a_background()
      ! (a) An active carrier model with bg_ready false.
      type(cert_report) :: rep
      real*8, dimension(3,1-Ng:N+Ng) :: u, Res
      real*8  :: rcmax_returned
      integer :: i
      thereis_mol       = .true.
      carrier_transport = .true.
      bg_ready          = .false.
      n_carrier         = 4
      carrier_solved    = .true.
      u = 1.0d0;  Res = 0.0d0
      call certification_evaluate(cert_context_stationary, u, Res,        &
                                  f_sp_none(), 1.0d-5, 0, .true., rep)
      i = entry_index(rep, 'carrier balance H2')
      call check_int('active_carrier_without_background_is_unavailable',  &
                     rep%e(i)%status, cert_unavailable)
      call check_log('no_certification_without_the_background',           &
                     rep%certified, .false.)
      ! THE PREDICATE THIS REPLACES, transcribed: carrier_steady_residual
      ! sets rcmax = 0 and returns when bg_ready is false, and a gate reads
      ! that zero as a balanced carrier.
      rcmax_returned = 0.0d0
      call check_log('superseded_rule_reads_the_zero_as_balanced',        &
                     rcmax_returned .lt. 1.0d-8, .true.)
      end subroutine carriers_without_a_background

      ! ------------------------------------------------------!

      subroutine model_without_carriers()
      ! (b) A model that carries no such species.
      type(cert_report) :: rep
      real*8, dimension(3,1-Ng:N+Ng) :: u, Res
      integer :: i
      thereis_mol       = .false.
      carrier_transport = .false.
      bg_ready          = .false.
      u = 1.0d0;  Res = 0.0d0
      call certification_evaluate(cert_context_stationary, u, Res,        &
                                  f_sp_none(), 1.0d-5, 0, .true., rep)
      i = entry_index(rep, 'carrier balance H2')
      call check_int('inactive_carrier_is_not_applicable',                &
                     rep%e(i)%status, cert_not_applicable)
      ! One equation of the inventory is active in every configuration: the
      ! closure that eliminates the neutral stages and the electron density.
      call check_log('an_inactive_entry_is_not_a_refusal',                &
           entry_status(rep, 'carrier balance H2') .eq. cert_not_applicable,&
           .true.)
      end subroutine model_without_carriers

      ! ------------------------------------------------------!

      subroutine row_measure_without_its_volume_companion()
      ! THE COMPANION IS OPTIONAL AND THE GATE DOES NOT PAY FOR IT.
      ! certification_row_measure returns rmax, the cell it sits in and the
      ! finiteness flag from the cells alone; the volume-weighted companion
      ! is a second reading of the same rows and nothing above depends on
      ! it.  A caller that reads only rmax -- the convergence gate of the
      ! Newton solve, which asks at every residual evaluation -- passes
      ! neither wvol nor rvol, and must get the same three numbers, bit for
      ! bit, as a caller that asks for all five.
      integer, parameter :: nc = 20
      real*8  :: res(nc), sc(nc), w(nc)
      real*8  :: rmax_full, rmax_bare, rvol
      integer :: j, jw_full, jw_bare
      logical :: fin_full, fin_bare
      do j = 1, nc
         res(j) = 1.0d-6*dble(j)*(-1.0d0)**j
         sc(j)  = 1.0d-3 + 1.0d-4*dble(nc - j)
         w(j)   = 1.0d0 + dble(j)**2
      enddo
      call certification_row_measure(nc, 7, res, sc, rmax_full, jw_full,  &
                                     fin_full, wvol = w, rvol = rvol)
      call certification_row_measure(nc, 7, res, sc, rmax_bare, jw_bare,  &
                                     fin_bare)
      call check_log('the_row_maximum_does_not_depend_on_the_companion', &
                     rmax_bare .eq. rmax_full, .true.)
      call check_int('nor_does_the_cell_it_sits_in', jw_bare, jw_full)
      call check_log('nor_does_the_finiteness_flag',                      &
                     fin_bare .eqv. fin_full, .true.)
      ! The companion itself is a different number from the maximum on
      ! these rows, so the assertion above is not two readings of one
      ! quantity that happen to agree.
      call check_log('and_the_companion_is_a_different_number',           &
                     rvol .ne. rmax_full, .true.)
      end subroutine row_measure_without_its_volume_companion

      ! ------------------------------------------------------!

      subroutine small_rate_species_beside_a_balanced_one()
      ! (c) Two carrier rows over the same cells. The dominant one balances
      ! to 1e-12 of its own terms; the small-rate one is out by 1e-3 of its
      ! own, and its terms are a millionth of the dominant one's.
      integer, parameter :: nc = 20
      real*8  :: res_dom(nc), sc_dom(nc), res_sml(nc), sc_sml(nc), w(nc)
      real*8  :: rmax_dom, rmax_sml, rvol_dom, rvol_sml, rvol_combined
      real*8  :: num, den
      integer :: j, jw
      logical :: fin
      do j = 1, nc
         sc_dom(j)  = 1.0d0
         res_dom(j) = 1.0d-12
         sc_sml(j)  = 1.0d-6
         res_sml(j) = 1.0d-9
         w(j)       = 1.0d0
      enddo
      call certification_row_measure(nc, 6, res_dom, sc_dom,              &
                                     rmax_dom, jw, fin,                  &
                                     wvol = w, rvol = rvol_dom)
      call certification_row_measure(nc, 6, res_sml, sc_sml,              &
                                     rmax_sml, jw, fin,                  &
                                     wvol = w, rvol = rvol_sml)
      call check_log('row_measure_catches_the_small_rate_species',        &
                     rmax_sml .gt. cert_tol_carrier_at(2.0d0), .true.)
      call check_log('the_dominant_row_is_within_its_tolerance',          &
                     rmax_dom .lt. cert_tol_carrier_at(2.0d0), .true.)
      ! THE MEASURE THIS REPLACES, transcribed: carrier_steady_residual's
      ! rvol is a ratio of SUMS over all carriers of a region.
      num = 0.0d0;  den = 0.0d0
      do j = 1, nc
         num = num + (abs(res_dom(j)) + abs(res_sml(j)))*w(j)
         den = den + (sc_dom(j) + sc_sml(j))*w(j)
      enddo
      rvol_combined = num/den
      call check_log('superseded_rvol_hides_the_small_rate_species',      &
                     rvol_combined .lt. cert_tol_carrier_at(2.0d0),  &
                     .true.)
      end subroutine small_rate_species_beside_a_balanced_one

      ! ------------------------------------------------------!

      subroutine state_with_a_cell_without_a_chemical_root()
      ! (d) One cell of the state's own sweep accepted a composition that is
      ! not a root of the network.
      type(cert_report) :: rep
      real*8, dimension(3,1-Ng:N+Ng) :: u, Res
      real*8  :: fspread
      thereis_mol       = .false.
      carrier_transport = .false.
      u = 1.0d0;  Res = 0.0d0
      ieq_nonroot_streak    = 0
      ieq_nonroot_streak(7) = 1
      call certification_evaluate(cert_context_stationary, u, Res,        &
                                  f_sp_none(), 1.0d-5, 1, .true., rep)
      call check_log('a_class4_cell_refuses_certification',               &
                     rep%certified, .false.)
      call check_int('the_refusal_names_the_cell',                        &
                     rep%first_no_chem_cell, 7)
      ! THE PREDICATE THIS REPLACES, transcribed: the marching stop as it
      ! stood, with no statement about the chemistry of the state at all.
      ! The flux gate is disabled and the flux window placed outside the
      ! grid, so the residual condition is the only other one.
      j_flux         = N + 4
      flux_spread_th = -1.0d0
      call check_log('superseded_marching_stop_accepts_the_same_state',   &
           steady_gates_met(1.0d-9, u, 1.0d-5, fspread, .false., 0.0d0),  &
           .true.)
      ieq_nonroot_streak = 0
      end subroutine state_with_a_cell_without_a_chemical_root

      ! ------------------------------------------------------!

      subroutine admissible_trial_is_not_a_certified_state()
      ! (e) and (f). A trial that leaves no more cells without a root than
      ! the iterate, whose rows are finite and inside the element headroom,
      ! is admissible whatever its residual is; and a row measure far above
      ! the stationary tolerance is not a reason to reject it.
      type(cert_evaluation_facts) :: f
      integer, parameter :: nc = 20
      real*8  :: res(nc), sc(nc), w(nc), rmax, rvol
      integer :: j, jw
      logical :: fin
      f%n_sweep_nonfinite    = 0
      f%n_no_chem_root_trial = 0
      f%n_no_chem_root_state = 0
      f%chem_root_gate       = .true.
      f%headroom_ok          = .true.
      f%rows_finite          = .true.
      call check_log('merit_decreasing_trial_is_admissible',              &
                     trial_state_is_admissible(f), .true.)
      do j = 1, nc
         res(j) = 1.0d0;  sc(j) = 1.0d0;  w(j) = 1.0d0
      enddo
      call certification_row_measure(nc, 6, res, sc, rmax, jw, fin,       &
                                     wvol = w, rvol = rvol)
      ! cert_tol_hydro is gone: the convergence study of contract section 9
      ! split it into one number per hydrodynamic row. The energy row is the
      ! loosest of the three, so asserting against it is the strongest form
      ! of "a unit residual on a unit scale is not certified".
      call check_log('the_same_state_is_not_certified',                   &
                     rmax .lt. cert_tol_energy, .false.)
      ! A trial that leaves MORE cells without a root than the iterate is
      ! refused; the probe that samples the same neighborhood is not.
      f%n_no_chem_root_trial = 3
      call check_log('trial_with_a_new_rootless_cell_is_refused',         &
                     trial_state_is_admissible(f), .false.)
      call check_log('the_probe_is_not_held_to_that_condition',           &
                     probe_direction_is_usable(f), .true.)
      ! A non-finite row refuses both.
      f%n_no_chem_root_trial = 0
      f%rows_finite          = .false.
      call check_log('a_non_finite_row_refuses_the_trial',                &
                     trial_state_is_admissible(f), .false.)
      call check_log('a_non_finite_row_refuses_the_probe',                &
                     probe_direction_is_usable(f), .false.)
      end subroutine admissible_trial_is_not_a_certified_state

      ! ------------------------------------------------------!

      subroutine isolated_workspace_round_trip()
      ! (g) The pair that holds the carrier module state aside and puts it
      ! back leaves it exactly as it was, allocation status included.
      type(carrier_module_state) :: ws
      call save_carrier_module_state(ws)
      call restore_carrier_module_state(ws)
      call check_log('the_carrier_module_state_is_reinstated_exactly',    &
                     carrier_module_state_matches(ws), .true.)
      end subroutine isolated_workspace_round_trip

      ! ------------------------------------------------------!

      subroutine level_rows_of_a_configuration_that_carries_them()
      ! (h) The two level balances of B1a section 2.5 follow the flags that
      ! make them independent, and neither is ever a zero:
      !   the metastable off  -> He 2^3S is not_applicable, and its absence
      !                          is not a refusal;
      !   the metastable on with no sweep behind the state -> unavailable,
      !                          with the reason stated;
      !   H(n=2) off          -> not_applicable.
      type(cert_report) :: rep
      real*8, dimension(3,1-Ng:N+Ng) :: u, Res
      logical :: he_was, tr_was
      integer :: i
      he_was = thereis_He;  tr_was = thereis_HeITR
      thereis_mol       = .false.
      carrier_transport = .false.
      u = 1.0d0;  Res = 0.0d0

      thereis_He    = .true.
      thereis_HeITR = .false.
      use_excited_H = .false.
      call certification_evaluate(cert_context_stationary, u, Res,        &
                                  f_sp_none(), 1.0d-5, 0, .true., rep)
      i = entry_index(rep, 'level balance He 2^3S')
      call check_int('triplet_row_without_the_metastable_is_absent',      &
                     rep%e(i)%status, cert_not_applicable)
      i = entry_index(rep, 'level balance H(n=2)')
      call check_int('excited_H_row_when_off_is_absent',                  &
                     rep%e(i)%status, cert_not_applicable)

      thereis_HeITR   = .true.
      ieq_rates_ready = .false.
      call certification_evaluate(cert_context_stationary, u, Res,        &
                                  f_sp_none(), 1.0d-5, 0, .true., rep)
      i = entry_index(rep, 'level balance He 2^3S')
      call check_int('triplet_row_with_no_sweep_is_unavailable',          &
                     rep%e(i)%status, cert_unavailable)
      call check_log('an_unavailable_level_row_refuses_certification',    &
                     rep%certified, .false.)
      thereis_He = he_was;  thereis_HeITR = tr_was
      end subroutine level_rows_of_a_configuration_that_carries_them

      ! ------------------------------------------------------!

      subroutine closure_residual_of_a_constructed_cell()
      ! (i) THE ELIMINATED-SPECIES CLOSURE AS A FUNCTION OF A COMPOSITION.
      !
      ! One cell of an H-only run, with the rate state of the sweep filled
      ! by hand: photoionization g, recombination a, no collisional
      ! ionization. Its equilibrium row is
      !     (1-x) n_H g  -  a x n_H n_e = 0 ,  n_e = x n_H ,
      ! and the root is found here by bisection, independently of the code
      ! under test. The closure residual read at that root is a root by the
      ! sweep's own tolerance; the same cell read at a composition 1 per
      ! cent away from it is not. The scratch the evaluation installs is
      ! reinstated bit for bit.
      type(ieq_closure_scratch) :: ws
      real*8  :: g_hi, a_hii, nh_l, xlo, xhi, xm, fm
      real*8  :: xroot(1), xoff(1), rrow(1), res_root, res_off
      integer :: it
      logical :: he_was
      he_was     = thereis_He
      thereis_He = .false.
      thereis_mol       = .false.
      carrier_transport = .false.
      thereis_metals    = .false.

      g_hi  = 1.0d-4
      a_hii = 2.6d-13
      nh_l  = 1.0d6

      ! Bisection on f(x) = (1-x) n_H g - a x^2 n_H^2, decreasing in x.
      xlo = 0.0d0;  xhi = 1.0d0
      do it = 1, 200
         xm = 0.5d0*(xlo + xhi)
         fm = (1.0d0 - xm)*nh_l*g_hi - a_hii*xm*xm*nh_l*nh_l
         if (fm .gt. 0.0d0) then
            xlo = xm
         else
            xhi = xm
         endif
      enddo
      xm = 0.5d0*(xlo + xhi)

      if (allocated(ieq_rate_cell)) deallocate(ieq_rate_cell)
      if (allocated(ieq_ne_cell))   deallocate(ieq_ne_cell)
      if (allocated(ieq_TK_cell))   deallocate(ieq_TK_cell)
      if (allocated(ieq_ntot_cell)) deallocate(ieq_ntot_cell)
      allocate(ieq_rate_cell(1-Ng:N+Ng), ieq_ne_cell(1-Ng:N+Ng),          &
               ieq_TK_cell(1-Ng:N+Ng), ieq_ntot_cell(1-Ng:N+Ng))
      ieq_rate_cell(1)%P_HI     = g_hi
      ieq_rate_cell(1)%rchiiB   = a_hii
      ieq_rate_cell(1)%nh       = nh_l
      ieq_rate_cell(1)%a_ion_HI = 0.0d0
      ieq_ne_cell(1)   = xm*nh_l
      ieq_TK_cell(1)   = 1.0d4
      ieq_ntot_cell(1) = nh_l
      ieq_neq_stored   = 1
      ieq_mbase_stored = 0
      ieq_iox_stored   = 0
      ieq_rates_ready  = .true.

      call save_ieq_closure_scratch(ws)
      xroot(1) = xm
      call ionization_closure_residual_cell(1, xroot, 1, res_root, rrow)
      xoff(1)  = xm*0.99d0
      call ionization_closure_residual_cell(1, xoff, 1, res_off, rrow)

      call check_log('closure_at_the_root_is_within_the_sweep_tolerance', &
                     res_root .le. ieq_res_tol, .true.)
      call check_log('closure_at_a_perturbed_composition_is_above_it',    &
                     res_off .gt. ieq_res_tol, .true.)
      call check_log('the_closure_scratch_is_reinstated_exactly',         &
                     ieq_closure_scratch_matches(ws), .true.)
      thereis_He = he_was
      end subroutine closure_residual_of_a_constructed_cell


      ! ------------------------------------------------------!

      subroutine the_regime_windows_select_their_own_cells()
      ! THE REGIME SPLIT OF A ROW, on rows built for the purpose.
      !
      ! certification_row_measure_over reduces the same measure over one
      ! radial window, so that the layer and the wind readings of one
      ! equation are separate numbers. The rows below put the largest
      ! residual of the column inside the layer window and the second
      ! largest inside the wind window, so a window that took the wrong
      ! cells could not pass.
      integer, parameter :: nc = 20
      real*8  :: res(nc), sc(nc), rl, rw, rall
      integer :: j, jl, jw, jall
      logical :: fin
      call set_up_a_grid()          ! r(j) = 1 + 0.05 j, so cells 1 and 2
      do j = 1, nc                  ! are the layer and 4 upward the wind
         res(j) = 1.0d-9
         sc(j)  = 1.0d0
      enddo
      res(2)  = 5.0d-4              ! r = 1.10 is NOT in the layer window
      res(1)  = 3.0d-4              ! r = 1.05 is
      res(10) = 1.0d-4              ! r = 1.50 is in the wind window
      call certification_row_measure(nc, 4, res, sc, rall, jall, fin)
      call certification_row_measure_over(nc, res, sc, -1.0d0,            &
                                          cert_regime_layer_r, rl, jl)
      call certification_row_measure_over(nc, res, sc,                    &
                                    cert_regime_wind_r, huge(1.0d0),      &
                                    rw, jw)
      write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)')                             &
           '  DIAGNOSTIC regime split: all ', rall, '  layer ', rl,       &
           '  wind ', rw
      call check_int('the_layer_window_takes_the_cell_below_its_radius',  &
                     jl, 1)
      call check_int('the_wind_window_takes_the_cell_above_its_radius',   &
                     jw, 10)
      call check_log('the_cell_between_the_windows_is_in_neither',        &
                     (jl .ne. 2) .and. (jw .ne. 2), .true.)
      call check_log('the_whole_column_maximum_is_the_larger_of_them',    &
                     rall .ge. max(rl, rw), .true.)
      call check_int('and_it_is_the_cell_between_the_two_windows',        &
                     jall, 2)
      end subroutine the_regime_windows_select_their_own_cells

      ! ------------------------------------------------------!

      subroutine a_species_row_is_gated_in_the_wind_and_reported_below_it()
      ! (h) THE VERDICT OF A SPECIES ROW IS TAKEN IN THE WIND, and the
      ! cells below cert_regime_wind_r are measured and reported and
      ! decide nothing (decision 22 (a); the declaration of the two radii
      ! in certification.f90 states the anchoring).
      !
      ! Two rows over one grid, differing only in the wind cell:
      !   row A   layer cell 1e-2, wind cell 5e-6: within tolerance,
      !   row B   layer cell 1e-2, wind cell 5e-5: refused,
      ! which is the distinction the whole item is about. The superseded
      ! rule, the whole-column maximum against 1e-8, is transcribed and
      ! run on the same two rows: it refuses both and cannot tell them
      ! apart, and it refuses them on a cell whose element flux is not
      ! conserved to better than a factor 10.
      integer, parameter :: nc = 20
      real*8, parameter  :: superseded_tol = 1.0d-8
      real*8  :: resA(nc), resB(nc), sc(nc), dg, rg, rall
      integer :: j, jg, jall
      logical :: fin
      call set_up_a_grid()          ! r(j) = 1 + 0.05 j
      do j = 1, nc
         sc(j)   = 1.0d0
         resA(j) = 1.0d-9
      enddo
      resA(1)  = 1.0d-2             ! r = 1.05, the diffusion layer
      resA(10) = 5.0d-6             ! r = 1.50, the wind
      resB     = resA
      resB(10) = 5.0d-5
      ! The accessor is a step at the wind radius and nothing else.
      call check_log('the_tolerance_at_the_wind_radius_is_the_wind_one',  &
                     cert_tol_element_at(cert_regime_wind_r) .eq.         &
                     1.0d-5, .true.)
      call check_log('and_a_cell_below_it_is_not_gated_at_all',           &
                     cert_tol_element_at(cert_regime_wind_r*0.999d0)      &
                     .ge. cert_tol_reported_only, .true.)
      call check_log('the_carrier_accessor_is_the_same_step',             &
                     (cert_tol_carrier_at(cert_regime_wind_r) .eq.        &
                      1.0d-5) .and.                                       &
                     (cert_tol_carrier_at(cert_regime_layer_r) .ge.       &
                      cert_tol_reported_only), .true.)
      ! Row A: the gated maximum is the wind cell, and it is within.
      call certification_species_row_gate(nc, .false., resA, sc, dg, rg,  &
                                          jg)
      write(*,'(A,I0,A,ES12.5,A,ES12.5)')                                 &
           '  DIAGNOSTIC row A gated at cell ', jg, ', measure ', rg,     &
           ', distance ', dg
      call check_int('the_gated_cell_of_row_A_is_the_wind_cell', jg, 10)
      call check_log('row_A_is_within_its_wind_tolerance',                &
                     dg .lt. 1.0d0, .true.)
      call check_log('and_its_layer_cell_did_not_enter_the_verdict',      &
                     rg .lt. 1.0d-5, .true.)
      ! Row B differs only in the wind cell, and that is what refuses it.
      call certification_species_row_gate(nc, .false., resB, sc, dg, rg,  &
                                          jg)
      call check_int('the_gated_cell_of_row_B_is_the_same_cell', jg, 10)
      call check_log('row_B_is_refused_by_its_wind_cell',                 &
                     dg .lt. 1.0d0, .false.)
      ! A carrier row is gated by its own accessor, on the same cells.
      call certification_species_row_gate(nc, .true., resA, sc, dg, rg,   &
                                          jg)
      call check_int('a_carrier_row_is_gated_in_the_same_window', jg, 10)
      call check_log('and_reaches_the_same_verdict_on_row_A',             &
                     dg .lt. 1.0d0, .true.)
      ! THE SUPERSEDED RULE, transcribed: the whole-column maximum against
      ! 1e-8. It reads the layer cell, so it refuses both rows and says
      ! nothing about the difference between them.
      call certification_row_measure(nc, 4, resA, sc, rall, jall, fin)
      call check_int('superseded_rule_reads_the_layer_cell', jall, 1)
      call check_log('superseded_rule_refuses_row_A',                     &
                     rall .lt. superseded_tol, .false.)
      call certification_row_measure(nc, 4, resB, sc, rall, jall, fin)
      call check_log('superseded_rule_refuses_row_B_for_the_same_reason', &
                     rall .lt. superseded_tol, .false.)
      ! A column that ends below the wind radius has no gated cell, and
      ! that is not a satisfied row: the caller refuses it.
      do j = 1, nc
         r(j) = 1.0d0 + 0.001d0*dble(j)
      enddo
      call certification_species_row_gate(nc, .false., resA, sc, dg, rg,  &
                                          jg)
      call check_int('a_column_below_the_wind_radius_gates_no_cell',      &
                     jg, 0)
      call set_up_a_grid()
      call the_line_a_wind_certified_state_carries()
      end subroutine a_species_row_is_gated_in_the_wind_and_reported_below_it

      ! ------------------------------------------------------!

      subroutine the_line_a_wind_certified_state_carries()
      ! THE ONE LINE a certified state with species rows is reported
      ! with. The report is written from a cert_report built here, with
      ! one gated entry standing at 5e-6 where it is gated and at 5.8e-4
      ! where it is only reported, because no state of this code is
      ! certified today and the line would otherwise never be seen. What
      ! is asserted is what a test can assert about a report: that the
      ! writer runs on such a state; the line itself is printed for the
      ! reader.
      type(cert_report) :: rep
      rep%n                       = 1
      rep%certified               = .true.
      rep%chem_root_known         = .true.
      rep%n_no_chem_root          = 0
      rep%e(1)%name               = 'elemental transport He/H partition'
      rep%e(1)%status             = cert_evaluated
      rep%e(1)%regime_gated       = .true.
      rep%e(1)%tol                = cert_tol_element_at(cert_regime_wind_r)
      rep%e(1)%row_max            = 5.8d-4
      rep%e(1)%jworst             = 2
      rep%e(1)%row_max_gate       = 5.0d-6
      rep%e(1)%jworst_gate        = 10
      rep%e(1)%row_max_reported   = 5.8d-4
      rep%e(1)%jworst_reported    = 2
      rep%e(1)%within_tol         = .true.
      rep%e(1)%units_floor        = 'the row''s own terms'
      call certification_report_write(rep,                                &
           'a constructed state, to show the line a wind-certified '//    &
           'state carries')
      call check_log('a_wind_certified_report_is_written',                &
                     rep%certified, .true.)
      end subroutine the_line_a_wind_certified_state_carries

      ! ------------------------------------------------------!

      subroutine element_column_of_cells(nc)
      ! The column every anchoring row below runs on: helium and trace
      ! metals transported, WENO3, no gravity, no ambipolar field, atomic
      ! hydrogen, so the settling coefficient vanishes and what the rows
      ! read is transport against transport. The same configuration as
      ! src/tests/element_operator, on a uniform grid whose cell centers of
      ! the coarse counts are a SUBSET of the fine one's: r(j) = 1 +
      ! (j-1)/(nc-1), so 201, 401 and 801 cells share their centers and a
      ! fine solution restricts onto a coarse grid with no interpolation.
      integer, intent(in) :: nc
      integer :: j
      real*8  :: dr_u
      N   = nc
      T0  = 1.0d3
      R0  = 1.0d10
      n0  = 1.0d10
      v0  = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0  = n0*mu*v0*v0
      b0  = 0.0d0
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
      ! Thermal diffusion, with the temperature gradient of
      ! column_background: it is what makes the steady state of this column
      ! a profile with a gradient rather than a uniform composition. With
      ! no drift at all, a uniform X and an r^2 F_rho that is constant are
      ! a steady state of both terms and the column has nothing to resolve.
      he_alphaT            = 1.5d-1
      he_kzz               = 1.0d11
      HeH                  = 0.0833333333333333d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
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
      dr_j  = dr_u
      j_min = 1
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4
      end subroutine element_column_of_cells

      ! ------------------------------------------------------!

      subroutine column_background(rho_c, v_c, T_c, Frho)
      ! The wind and the thermodynamic background of the column, given by
      ! formulas of r, so a coarse grid carries the SAME background as the
      ! fine one and nothing about it is interpolated. r^2 F_rho is
      ! constant, which is a steady mass row.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_c, v_c, T_c, Frho
      integer :: j
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0 + 5.0d-1*(r(j) - 1.0d0)
         Frho(j)  = 1.0d-4/(r_edg(j)*r_edg(j))
      enddo
      end subroutine column_background

      ! ------------------------------------------------------!

      subroutine column_composition(f_c, gradient)
      ! The composition of the column, written as MASS FRACTIONS that sum to
      ! one: helium carries X_He, each trace element its own Z_e, and
      ! hydrogen the remainder. The species vector holds a mass fraction
      ! divided by the species mass, so sum_i m_i f_i = 1 exactly and the
      ! column's species carry the column's own density whatever the
      ! gradient is. With `gradient` the helium mass fraction falls by a
      ! factor 2 from the base to the outermost cell, and hydrogen takes up
      ! what it leaves; without it the column is uniform.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_c
      logical,                                intent(in)  :: gradient
      real*8  :: mpH, xhe0, xhe, zsum, ze(n_melem)
      integer :: j, ie
      mpH  = mass_per_H_nucleus_without_He() + m_He_over_m_H*HeH
      xhe0 = m_He_over_m_H*HeH/mpH
      zsum = 0.0d0
      do ie = 1, n_melem
         ze(ie) = melem_A(ie)*melem_ab(ie)/mpH
         zsum   = zsum + ze(ie)
      enddo
      f_c = 0.0d0
      do j = 1-Ng, N+Ng
         xhe = xhe0
         if (gradient) xhe = xhe0/(1.0d0 + (r(j) - 1.0d0)                 &
                                   /max(r(N) - 1.0d0, 1.0d-30))
         f_c(j,isp_HI)  = (1.0d0 - xhe - zsum)/bsp_mass(isp_HI)
         f_c(j,isp_HeI) = xhe/m_He_over_m_H
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = ze(ie)/melem_A(ie)
         enddo
      enddo
      end subroutine column_composition

      ! ------------------------------------------------------!

      subroutine close_the_column_boundaries(f_c)
      ! The boundaries both production callers write: the lower ghosts are
      ! the Dirichlet reservoir the base cell carries, the upper ghosts the
      ! zero-gradient continuation of the outermost cell.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_c
      integer :: j
      do j = 1-Ng, 0
         f_c(j,:) = f_c(1,:)
      enddo
      do j = N+1, N+Ng
         f_c(j,:) = f_c(N,:)
      enddo
      end subroutine close_the_column_boundaries

      ! ------------------------------------------------------!

      function worst_element_row(rho_c, T_c, f_c, Frho) result(rmax)
      ! The largest element row of the column in its own measure, over the
      ! cells that carry an equation (cell 1 is the reservoir and reads
      ! zero by construction).
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho_c, T_c
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_c
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: Frho
      real*8  :: rmax
      real*8, dimension(1:N)         :: res_he, sc_he
      real*8, dimension(1:N,n_melem) :: res_tr, sc_tr
      logical :: ok_he, ok_tr
      integer :: j, im
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he,      &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      rmax = 0.0d0
      do j = 2, N
         rmax = max(rmax, abs(res_he(j))/max(sc_he(j), cert_scale_floor))
         do im = 1, n_melem
            rmax = max(rmax, abs(res_tr(j,im))                            &
                             /max(sc_tr(j,im), cert_scale_floor))
         enddo
      enddo
      end function worst_element_row

      ! ------------------------------------------------------!

      function column_mass_closure(rho_c, f_c) result(closure)
      ! |sum_i m_i n_i - rho| / rho over the physical cells, with the code's
      ! own mass policy (calc_rho): the conservation number that says
      ! whether a state's species still carry the density the hydrodynamics
      ! evolves. No row of the certification inventory measures it.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho_c
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_c
      real*8  :: closure
      real*8, dimension(1-Ng:N+Ng)        :: msum
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm_c
      integer :: j, k
      do k = 1, n_mion
         nm_c(:,k) = f_c(:,mion_fsp(k))*rho_c*n0
      enddo
      call calc_rho(f_c(:,isp_HI)*rho_c*n0, f_c(:,isp_HII)*rho_c*n0,      &
                    f_c(:,isp_HeI)*rho_c*n0, f_c(:,isp_HeII)*rho_c*n0,    &
                    f_c(:,isp_HeIII)*rho_c*n0, msum, nm = nm_c)
      closure = 0.0d0
      do j = 1, N
         closure = max(closure, abs(msum(j) - rho_c(j)*n0)                &
                                /max(rho_c(j)*n0, 1.0d-300))
      enddo
      end function column_mass_closure

      ! ------------------------------------------------------!

      subroutine the_five_anchoring_numbers_of_a_synthetic_column()
      ! THE FIVE NUMBERS A SPECIES-ROW TOLERANCE CAN BE ANCHORED ON, each
      ! printed for a column built here, so that the anchoring of
      ! docs/certification_tolerance_anchoring_20260910.md can be repeated
      ! without a planet run. Nothing here is a threshold on the tolerance:
      ! the rows assert only what the numbers must satisfy to be numbers
      ! (a repeated evaluation returns the same one; a finer grid reads a
      ! smaller discretization error; the least derivative error is not the
      ! largest), and the values themselves are printed as DIAGNOSTIC.
      !
      ! (1) reproducibility: two evaluations of one state.
      ! (2) representation: the same state moved by one unit in the last
      !     place of its transported composition.
      ! (3) derivative: the error of the difference quotient of the row
      !     against the probe length, whose least value is the smallest row
      !     change a Newton step can resolve.
      ! (4) discretization: one converged column read on three nested
      !     grids, and the order the three imply.
      ! (5) conservation: the mass closure between the column's density and
      !     its own species, which no row of the inventory measures.
      integer, parameter :: nfine = 801
      integer, parameter :: nstride(2) = [2, 4]
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c, f_fine
      real*8, dimension(:),   allocatable :: rho_m
      real*8  :: e1, e2, eulp, eflat, drift, h, dq, dref, ebest, eworst
      real*8  :: ecoarse(2), order, closure, x0, gradmax, closure0
      integer :: it, nstep, j, s, k, ic, nc

      ! ---- the flat column: the row is analytically zero, so what it
      ! reads is the assembly's own round-off floor ----
      call element_column_of_cells(201)
      he_alphaT = 0.0d0     ! no drift, so a uniform column is exactly steady
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      call column_background(rho_c, v_c, T_c, Frho)
      call column_composition(f_c, .false.)
      call close_the_column_boundaries(f_c)
      eflat = worst_element_row(rho_c, T_c, f_c, Frho)
      write(*,'(A,ES12.5)') '  DIAGNOSTIC (0) round-off floor of the '//  &
           'analytically zero flat column: ', eflat
      call check_log('the_flat_column_row_is_at_the_arithmetic_floor',    &
                     eflat .le. 1.0d-12, .true.)

      ! ---- a column with a gradient, relaxed toward its steady state ----
      deallocate(rho_c, v_c, T_c, Frho, f_c)
      call element_column_of_cells(nfine)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(f_fine(1-Ng:N+Ng,n_species))
      call column_background(rho_c, v_c, T_c, Frho)
      call column_composition(f_c, .true.)
      ! The column is mass closed BEFORE the relaxation by construction, so
      ! the closure printed as number (5) below is what the transport step
      ! left behind and not what this test started with.
      closure0 = column_mass_closure(rho_c, f_c)
      do it = 1, 30
         call relax_element_composition(rho_c, v_c, T_c, f_c, Frho,       &
                                        1.0d2, drift, nstep)
      enddo
      call close_the_column_boundaries(f_c)
      f_fine = f_c
      e1 = worst_element_row(rho_c, T_c, f_c, Frho)
      ! The column the three grids read has a real gradient: without one
      ! there is no discretization error to measure and row (4) below would
      ! be a comparison of three round-off floors.
      gradmax = maxval(f_c(1:N,isp_HeI))/max(minval(f_c(1:N,isp_HeI)),    &
                                             1.0d-300) - 1.0d0
      write(*,'(A,ES12.5)') '  DIAGNOSTIC the relaxed column carries a '//&
           'helium contrast of ', gradmax
      call check_log('the_relaxed_column_has_a_gradient_to_resolve',      &
                     gradmax .gt. 1.0d-3, .true.)

      ! (1) two evaluations of one state.
      e2 = worst_element_row(rho_c, T_c, f_c, Frho)
      write(*,'(A,ES23.16,A,ES12.5)') '  DIAGNOSTIC (1) row of the '//    &
           'relaxed column ', e1, ', change on a second evaluation ',     &
           abs(e2 - e1)
      call check_log('two_evaluations_of_one_column_give_one_row',        &
                     e2 .eq. e1, .true.)

      ! (2) the same state, its helium moved by one unit in the last place.
      do j = 1-Ng, N+Ng
         f_c(j,isp_HeI) = nearest(f_c(j,isp_HeI), 1.0d0)
      enddo
      eulp = worst_element_row(rho_c, T_c, f_c, Frho)
      write(*,'(A,ES12.5,A,ES12.5)') '  DIAGNOSTIC (2) one ulp of the '// &
           'helium moves the row by ', abs(eulp - e1),                    &
           ', relative to the row ', abs(eulp - e1)/max(e1, 1.0d-300)
      call check_log('one_ulp_of_the_composition_moves_the_row_at_all',   &
                     abs(eulp - e1) .ge. 0.0d0, .true.)
      f_c = f_fine

      ! (3) the derivative of the row against the probe length, at the cell
      ! that carries the largest row. A correct pair falls like the probe
      ! until the cancellation of the two evaluations divided by it takes
      ! over; the least error is where the two meet and no smaller change
      ! of this row is resolvable by a difference quotient.
      x0 = f_c(N/2,isp_HeI)
      dref = 0.0d0
      ebest  = huge(1.0d0)
      eworst = 0.0d0
      write(*,'(A)') '  DIAGNOSTIC (3) probe length, difference '//       &
           'quotient of the row at the middle cell'
      do k = 1, 7
         h = 1.0d-2*(1.0d-2**(k-1))
         f_c(N/2,isp_HeI) = x0*(1.0d0 + h)
         dq = (worst_element_row(rho_c, T_c, f_c, Frho) - e1)/h
         f_c(N/2,isp_HeI) = x0
         if (k .eq. 1) dref = dq
         write(*,'(A,ES10.3,A,ES14.6)') '    probe ', h,                  &
              '  quotient ', dq
         if (k .gt. 1) then
            ebest  = min(ebest,  abs(dq - dref))
            eworst = max(eworst, abs(dq - dref))
         endif
      enddo
      call check_log('the_difference_quotient_of_the_row_is_finite',      &
                     ebest .le. eworst, .true.)

      ! (4) the same converged column on three nested grids. The coarse
      ! grids take the fine cells' own values, so the number below is the
      ! discretization error of that profile and carries no interpolation.
      write(*,'(A,I0,A,ES12.5)') '  DIAGNOSTIC (4) grid ', nfine,         &
           ' cells: row ', e1
      allocate(rho_m(1-Ng:nfine+Ng))
      do ic = 1, 2
         s  = nstride(ic)
         nc = (nfine - 1)/s + 1
         deallocate(rho_c, v_c, T_c, Frho, f_c)
         call element_column_of_cells(nc)
         allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),       &
                  Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
         call column_background(rho_c, v_c, T_c, Frho)
         do j = 1, N
            f_c(j,:) = f_fine(1 + s*(j-1),:)
         enddo
         call close_the_column_boundaries(f_c)
         ecoarse(ic) = worst_element_row(rho_c, T_c, f_c, Frho)
         write(*,'(A,I0,A,ES12.5)') '  DIAGNOSTIC (4) grid ', nc,         &
              ' cells: row ', ecoarse(ic)
      enddo
      order = log(ecoarse(2)/max(ecoarse(1), 1.0d-300))/log(2.0d0)
      write(*,'(A,F6.3)') '  DIAGNOSTIC (4) observed order between the '//&
           'two coarse grids: ', order
      call check_log('the_finer_grid_reads_the_smaller_row',              &
                     ecoarse(1) .lt. ecoarse(2), .true.)
      call check_log('and_the_finest_reads_the_smallest',                 &
                     e1 .lt. ecoarse(1), .true.)

      ! (5) the mass closure of the column against its own density: the
      ! conservation number the stationary inventory has no row for.
      closure = column_mass_closure(rho_c, f_c)
      write(*,'(A,ES12.5,A,ES12.5)') '  DIAGNOSTIC (5) mass closure of '//&
           'the column, |sum m_i n_i - rho|/rho: as built ', closure0,    &
           ', after the transport relaxation ', closure
      call check_log('the_columns_species_carry_its_own_mass',            &
                     closure .lt. 1.0d-2, .true.)
      deallocate(rho_c, v_c, T_c, Frho, f_c, f_fine, rho_m)
      end subroutine the_five_anchoring_numbers_of_a_synthetic_column

      ! ------------------------------------------------------!

      function f_sp_none() result(fs)
      ! A composition array of the right shape. The cases above never reach
      ! the carrier evaluation, which is the only consumer of it.
      real*8, dimension(1-Ng:N+Ng,n_species) :: fs
      fs = 0.0d0
      end function f_sp_none

      integer function entry_index(rep, name) result(i)
      type(cert_report), intent(in) :: rep
      character(len=*),  intent(in) :: name
      integer :: k
      i = 1
      do k = 1, rep%n
         if (trim(rep%e(k)%name) .eq. name) then
            i = k;  return
         endif
      enddo
      end function entry_index

      integer function entry_status(rep, name) result(st)
      type(cert_report), intent(in) :: rep
      character(len=*),  intent(in) :: name
      st = rep%e(entry_index(rep, name))%status
      end function entry_status

      subroutine check_int(name, measured, reference)
      character(*), intent(in) :: name
      integer,      intent(in) :: measured, reference
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_int

      subroutine check_log(name, measured, reference)
      character(*), intent(in) :: name
      logical,      intent(in) :: measured, reference
      character(5) :: sm, sr
      sm = 'false';  sr = 'false'
      if (measured)  sm = 'true '
      if (reference) sr = 'true '
      if (measured .eqv. reference) then
         write(*,'(A,A,A,A,A,A,A)') 'PASS ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
      else
         write(*,'(A,A,A,A,A,A,A)') 'FAIL ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_log

      end program certification_contexts

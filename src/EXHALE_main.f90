      program Hydro_ioniz

      use global_parameters
      use species_table, only: n_mion, mion_fsp
      use Read_input
      use Initialization
      ! What the restart file said it was produced under (parsed by load_IC).
      use IC_load, only: ic_coupling_present, ic_sec_ion_active,          &
                        ic_rec_method
      use setup_report
      use eval_time_step
      use energy_semi_implicit
      use utils
      use composition, only: get_species_densities, comp_T_from_p,         &
                             comp_p_from_T, element_ratio_HeH,             &
                             n_cells_he_singlet_clamped
      use binary_element_diffusion, only: element_diffusion_step,          &
                                          relax_element_composition
      use diffusive_photochemistry, only: photochemical_transport_step,   &
                                          relax_photochemical_composition,&
                                          carrier_transport_diagnostics, &
                                          carrier_steady_residual,      &
                                          carrier_name,                 &
                                          carrier_drift_location
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
                                     residual_row_scale
      use viscous_conduction, only: transport_active, viscous_conduction_step
      use steady_newton, only: neq_newton, pack_U, unpack_U, newton_residual, &
                               eval_residual, frozen_residual,              &
                               build_banded_jac, band_matvec,               &
                               kl_jac, ku_jac, solve_steady_ptc,           &
                               solve_steady_jfnk, set_base_fix,          &
                               set_carrier_unknown,                      &
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
      real*8  :: upmap_bc_roundtrip = 0.0d0
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
      real*8, dimension(:,:), allocatable :: f_sp_w
      integer :: nbitdiff_b, nbitdiff_c, nvj_det, nbitdiff_ctrl
      real*8  :: dmax_b, dmax_c, det_sc_b, det_sc_c, det_rs, det_sc_ctrl
      ! Banded-Jacobian self-test scratch (EXHALE_JAC_TEST hook)
      real*8, allocatable :: abjac(:,:), rdir(:), Jr(:), dFD(:), F0f(:)
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
      real*8 :: mom_max,mom_min

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

      ! Newton-diffusion co-convergence (Solver: Newton + He_diffusion):
      ! outer iteration alternating the JFNK steady solve with diffusion
      ! relaxation of the He/H field at the converged wind.
      real*8  :: comp_drift, comp_drift_prev, comp_omega
      ! Drift and step count of the carrier relaxation, the third
      ! participant of the same damped Picard iteration.
      real*8  :: carrier_drift
      real*8  :: crc_max, crc_vol, crc_leg
      ! How far one fixed-wind carrier pass may move the composition, as a
      ! fraction of the largest H2 mixing ratio on the grid. Measured
      ! boundary of the steady solve's tolerance on the hot Uranus hand-off:
      ! it still returns info = 0 at 0.2 and stops at 0.4, so this stands a
      ! factor 20 inside it (docs/p50_carrier_wind_alternation.md).
      ! The bound on one carrier pass, as a fraction of the largest H2 mixing
      ! ratio on the grid (relax_photochemical_composition).  1e-2 is the
      ! measured safe value (docs/p50_carrier_wind_alternation.md: the wind
      ! solve still accepts the handed-over state at 0.2 and stops accepting
      ! at 0.4).  EXHALE_CARRIER_TRUST overrides it and EXHALE_OUTER_PASSES
      ! the pass cap, so the alternation can be run as a CONTINUATION that
      ! walks the H2 front to where the carrier equation puts it (section
      ! 163); both are diagnostics and neither is an input key.
      real*8  :: carrier_trust = 1.0d-2
      integer :: outer_pass_cap = 20

      integer :: crc_j, crc_ic, cdl_j, cdl_ic
      integer :: it_diff, kd, kc

      ! Admissibility record of the low-Mach contact-mode dissipation: how big
      ! the artificial stress got against the physical momentum flux it was
      ! added to, and how far out its Mach gate stayed open. Only written when
      ! "Low-Mach damping" is on.
      real*8  :: lowmach_ratio, lowmach_r_peak, lowmach_r_gate

      ! Maximum eigenvalue      
      real*8 :: alpha

      ! Temporal step (global) and cell-by-cell pseudo-time steps
      real*8 :: dt
      real*8, dimension(:), allocatable :: dt_loc

      ! Positivity step control (see the retry_step loop): bisections taken on
      ! the current step, bisections over the run, and steps that needed at
      ! least one. n_dt_halve_max bounds the bisection: 20 halvings is a factor
      ! 1e-6 on dt, past which no step size repairs the state.
      integer :: n_dt_halve, n_dt_halvings, n_steps_dt_halved
      integer, parameter :: n_dt_halve_max = 20

      ! Did the local first-order flux correction restore the admissible set
      ! for the stage it was applied to? False hands the stage to the dt
      ! bisection above.
      logical :: flux_corr_ok

      ! Mdot value
      real*8 :: Mdot
      
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
      
      ! Open output report file 
      open(unit = outfile, file = 'EXHALE_setup.out')
      
      !------------------------------------------------!
      
      ! Read planetary parameters from file
      call input_read

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
         call Apply_BC(u)
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
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
      ! (env EXHALE_RESID_DETERMINISM=1.)  Evaluate F at the loaded state, at
      ! two other states, then at the loaded state again, REUSING ONE OUTPUT
      ! WORKSPACE throughout so that it carries each evaluation's composition
      ! into the next call.  The seed is named in every call and is the same
      ! every time, so a residual that is a state function returns the same
      ! bits; anything else means the workspace's leftover was read.
      !
      ! The interface is what makes this pass: eval_residual takes the seed as
      ! one argument and writes the composition to another, and the two may not
      ! be the same array, so a caller CANNOT hand the sweep its own previous
      ! answer even by accident.  Before that separation the same test read
      ! 1406 of 1500 entries differing by 3.0e-3 of the row scale
      ! (Update_EXHALE.md section 146.3).
      call get_environment_variable('EXHALE_RESID_DETERMINISM', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion
         if (sec_ion_active) sec_ion_armed_step = 0
         allocate(Yvec(neq_newton()), Fa(neq_newton()), Fb(neq_newton()),   &
                  Fc(neq_newton()), Ytry1(neq_newton()), Ytry2(neq_newton()))
         allocate(f_sp_w(1-Ng:N+Ng,n_species))
         call pack_U(u, Yvec)
         ! Two trial states of the size a line search takes.
         nvj_det = neq_newton()/N
         Ytry1 = Yvec;  Ytry2 = Yvec
         do j = 1, N
            Ytry1(nvj_det*(j-1)+3) = Ytry1(nvj_det*(j-1)+3)*(1.0d0+1.0d-6)
            Ytry2(nvj_det*(j-1)+2) = Ytry2(nvj_det*(j-1)+2)*(1.0d0+1.0d-4)
         enddo
         ! CONTROL: the same state twice in a row, nothing in between. Any
         ! difference here is state the evaluation leaves behind for itself,
         ! not something a trial injected.
         call newton_residual(Yvec, f_sp, f_sp_w, Fa)
         call newton_residual(Yvec, f_sp, f_sp_w, Fb)
         nbitdiff_ctrl = 0
         det_sc_ctrl   = 0.0d0
         do j = 1, neq_newton()
            if (Fa(j) .ne. Fb(j)) nbitdiff_ctrl = nbitdiff_ctrl + 1
         enddo
         do j = 1, N
            do k = 1, 3
               det_rs = residual_row_scale(k, j, u)
               det_sc_ctrl = max(det_sc_ctrl,                               &
                    abs(Fa(nvj_det*(j-1)+k)-Fb(nvj_det*(j-1)+k))/det_rs)
            enddo
         enddo
         write(*,'(A,I0,A,I0,A,ES11.3)') ' (resid_determinism) control,'//  &
              ' F(Y0) twice with nothing between: ', nbitdiff_ctrl,         &
              ' of ', neq_newton(), ' entries differ, row-scale max ',      &
              det_sc_ctrl
         ! One workspace, reused; one seed, named every time.
         call newton_residual(Yvec,  f_sp, f_sp_w, Fa)
         call newton_residual(Ytry1, f_sp, f_sp_w, Fc)
         call newton_residual(Yvec,  f_sp, f_sp_w, Fb)
         call newton_residual(Ytry2, f_sp, f_sp_w, Fc)
         call newton_residual(Yvec,  f_sp, f_sp_w, Fc)
         nbitdiff_b = 0;  nbitdiff_c = 0
         dmax_b = 0.0d0;  dmax_c = 0.0d0
         det_sc_b = 0.0d0;  det_sc_c = 0.0d0
         do j = 1, neq_newton()
            if (Fa(j) .ne. Fb(j)) nbitdiff_b = nbitdiff_b + 1
            if (Fa(j) .ne. Fc(j)) nbitdiff_c = nbitdiff_c + 1
            dmax_b = max(dmax_b, abs(Fa(j)-Fb(j))/max(abs(Fa(j)),1.0d-300))
            dmax_c = max(dmax_c, abs(Fa(j)-Fc(j))/max(abs(Fa(j)),1.0d-300))
         enddo
         do j = 1, N
            do k = 1, 3
               det_rs = residual_row_scale(k, j, u)
               det_sc_b = max(det_sc_b,                                     &
                    abs(Fa(nvj_det*(j-1)+k)-Fb(nvj_det*(j-1)+k))/det_rs)
               det_sc_c = max(det_sc_c,                                     &
                    abs(Fa(nvj_det*(j-1)+k)-Fc(nvj_det*(j-1)+k))/det_rs)
            enddo
         enddo
         write(*,'(A)') ' (resid_determinism) F(Y0) re-evaluated after a'//  &
              ' trial, one workspace reused throughout:'
         write(*,'(A,I0,A,I0,A,ES11.3)') '   after trial 1: ', nbitdiff_b,   &
              ' of ', neq_newton(), ' entries differ, max relative ', dmax_b
         write(*,'(A,I0,A,I0,A,ES11.3)') '   after trial 2: ', nbitdiff_c,   &
              ' of ', neq_newton(), ' entries differ, max relative ', dmax_c
         write(*,'(A,ES11.3,A,ES11.3)') '   in row-scale units: after 1 ',   &
              det_sc_b, ', after 2 ', det_sc_c
         if (nbitdiff_b .gt. 0) then
            write(*,'(A)') '   cells and rows that moved (first 12):'
            im = 0
            do j = 1, N
               do k = 1, 3
                  if (Fa(nvj_det*(j-1)+k) .ne. Fb(nvj_det*(j-1)+k)) then
                     im = im + 1
                     if (im .le. 12) write(*,'(A,I5,A,F9.5,A,I2,A,ES11.3)') &
                          '     j=', j, ' r=', r(j), ' row=', k,            &
                          '  |dF|/s = ',                                    &
                          abs(Fa(nvj_det*(j-1)+k)-Fb(nvj_det*(j-1)+k))      &
                          /residual_row_scale(k, j, u)
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
         det_rs = max(det_sc_ctrl, 1.0d-9)
         if (det_sc_b .le. det_rs .and. det_sc_c .le. det_rs) then
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
         k = neq_newton()
         allocate(Yvec(k), Fvec(k), F0f(k), rdir(k), Jr(k), dFD(k))
         allocate(abjac(2*kl_jac+ku_jac+1, k))
         call pack_U(u, Yvec)
         f_sp_test = f_sp
         ! F0 + frozen heat/cool at Y, then the banded Jacobian
         call eval_residual(Yvec, f_sp, f_sp_test, Fvec, heat0, cool0, npart0)
         call frozen_residual(Yvec, npart0, heat0, cool0, F0f)
         call build_banded_jac(Yvec, npart0, heat0, cool0, abjac)
         ! Deterministic probe direction (max|r| = 1)
         do j = 1, k
            rdir(j) = sin(0.1d0*dble(j))
         enddo
         rdir = rdir/maxval(abs(rdir))
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
         write(*,*) '  transport: True". The steady PTC/JFNK route holds'
         write(*,*) '  the composition at its own local root and would'
         write(*,*) '  undo the transported ionization state. Use the'
         write(*,*) '  marching route for this option. Aborting.'
         error stop 1
      endif
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         if (sec_ion_active) sec_ion_armed_step = 0
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
         call steady_wind_with_element_diffusion(3000, dt,             &
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

      ! A state read from a file starts with the composition the file could
      ! hold, not the composition that state has. Put it on the equilibrium
      ! the first hydro step would otherwise jump to (see the routine).
      if ((do_load_IC .or. ic_mode .eq. 4) .and. reload_equilibrate)       &
         call equilibrate_loaded_composition

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
            if (do_profile) tp_step0 = omp_get_wtime()

            !-------------------------------------------------!

            !--- Thermodynamic evolution ---!

            ! Save previous step solution
            ! BEFORE u_old is taken: the retry loop restarts every attempt
            ! from u_old, so a restore made after this line would be undone
            ! by the first `u = u_old` inside retry_step.
            if (upmap_n .gt. 0) call update_map_begin_step

            u_old = u

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
            n_dt_halve = 0
            retry_step: do

            u = u_old

            ! FIRST RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            if (do_profile) tp_b = omp_get_wtime()
            call reconstruction_continuation_rhs(u,alpha,WL,WR,dF,S)

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

            !----------------------------

            ! SECOND RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call reconstruction_continuation_rhs(u1,alpha,WL,WR,dF,S)

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

            !----------------------------

            ! THIRD RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call reconstruction_continuation_rhs(u2,alpha,WL,WR,dF,S)

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

            exit retry_step
            enddo retry_step

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

            ! He/H diffusive separation: advect + diffuse the He element
            ! ratio (updates the He/H split in f_sp; ionization equilibrium below
            ! then re-solves the stages, conserving the new element amounts).
            ! No-op unless he_diffusion is set (byte-identical when off).
            if (he_diffusion)                                         &
               call element_diffusion_step(rho,v,T,f_sp,dt_loc)

            ! Vertical transport of the molecular carriers (H2, OH, H2O,
            ! CO) with their chemistry, in the same operator-split slot the
            ! element diffusion occupies and for the same reason: after the
            ! hydro update, before the ionization solve, which then solves
            ! the remaining stages against the transported partition rather
            ! than recomputing it. No-op unless the oxygen chemistry is on
            ! with its transport (byte-identical when off).
            if (do_profile) tp_b = omp_get_wtime()
            call photochemical_transport_step(rho,v,T,f_sp,dt_loc)
            if (do_profile) tp_trans = tp_trans + (omp_get_wtime() - tp_b)

            ! Refresh the lagged H(n=2) Balmer source + heating from
            ! the current state before the ionization/energy solve.
            if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)

            ! Evaluate ionization equilibrium
            if (do_profile) tp_a = omp_get_wtime()
            call ioniz_eq(T,rho,f_sp,heat,cool,eta)
            if (do_profile) tp_ion = tp_ion + (omp_get_wtime() - tp_a)
            ! Evaluate partial densities, ne and n_tot (single policy point)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                       nheiii,nheiTR,nm,ne,n_tot)

            ! Evaluate updated pressure
            call comp_p_from_T(T,n_tot,ne,p)
            
            ! Convert to primitive profiles
            W(1,:) = rho
            W(2,:) = v
            W(3,:) = p
            
            ! Revert to conservative 
            call W_to_U(W,u)

            if (upmap_n .gt. 0) then
               u_umB = u
               upmap_dfsp = maxval(abs(f_sp(1:N,:) - f_sp_umSave(1:N,:)))
            endif
  	
  	      !------------------------------------------------!
  	      
  	      ! Evolve solution in time due to source terms
            !     using a semi-implicit step (default) or the original
            !     explicit forward-Euler update ("Energy solver: Explicit")

            if (use_semi_implicit_energy) then
               if (do_profile) tp_b = omp_get_wtime()
               call solve_energy_semi_implicit(u,W,dt_loc,heat,cool,f_sp,count)
               if (do_profile) tp_energy = tp_energy + (omp_get_wtime() - tp_b)
            else
               u(3,:) = u(3,:) + dt_loc*(heat - cool)
            endif

		call Apply_BC(u)

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
               call viscous_conduction_step(u,W,T,n_tot+ne,dt_loc)
               call Apply_BC(u)
            endif

            if (upmap_n .gt. 0) u_umD = u

            ! Periodic Shapiro low-pass filter to damp the gravity-unbalanced
            ! sound waves (base breathing), as in CETIMB (Koskinen et al. 2013a).
            if (shapiro_eps .gt. 0.0d0 .and.                             &
                mod(count, shapiro_every) .eq. 0) then
               call shapiro_filter(u)
               call Apply_BC(u)
            endif

            !------------------------------------------------!

            if (upmap_n .gt. 0) call update_map_end_step

            ! Convert to physical variables and extract profiles
            call U_to_W(u,W)
            rho = W(1,:)
            v   = W(2,:)
            p   = W(3,:)
            E   = u(3,:)

            ! Is the lower boundary still admitting SUBSONIC inflow? Measured
            ! here because W now holds the state every Apply_BC of this step
            ! has finished writing. Diagnostic only: it reads W and changes
            ! nothing, so a run that never crosses is the run an unguarded
            ! build produces.
            call check_base_inflow_is_subsonic(W,count)
            
            ! Evaluate ionized densities, ne and n_tot (single policy point)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                       nheiii,nheiTR,nm,ne,n_tot)

  	      ! Temperature profile
  	      call comp_T_from_p(p,n_tot,ne,T)

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
            
            ! Maximum and minimum value for momentum
            mom_max = maxval(abs(mom(j_min:N)))
            mom_min = minval(abs(mom(j_min:N)))

            ! Evaluate relative momentum variation. Guard mom_min against an exact
            ! zero (a transient zero-momentum point can appear during relaxation of
            ! a deep-RLOF wind whose base sits close to L1); a zero here would give
            ! du = Inf and a spurious "converged" exit on the first step.
            du = abs((mom_max-mom_min)/max(mom_min, 1.0d-30))

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

            ! Evaluate variation of time derivative          
      	u(2,:) = u(2,:) + 1.0e-16 ! To avoid division by zero          
	      
	      ! --- Infty-norm
	      dtu = max(maxval(abs(1.0-u(1,j_min:N)/u_old(1,j_min:N))), &
	      	    maxval(abs(1.0-u(2,j_min:N)/u_old(2,j_min:N))))
		dtu = max(maxval(abs(1.0-u(3,j_min:N)/u_old(3,j_min:N))), &
	      	    dtu)

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
               flux_spread = (maxval(mom(j_min:N)) - minval(mom(j_min:N)))   &
                    /max(abs(sum(mom(j_min:N))/dble(N - j_min + 1)), 1.0d-30)
               ! carrier_is_unknown = .false.: the marching loop's carriers
               ! are moved by the operator-split transport step, not solved
               ! for, so this state has no carrier residual of a solve to
               ! gate on. A coupled-carrier run reaches its answer through
               ! the JFNK finish, which does state one.
               gates_met_now = steady_gates_met(resid_max, u, resid_th,      &
                                                fspread_gate,               &
                                                .false., 0.0d0)
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
                  if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                     stall_count = stall_count + 1
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
                  if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                     stall_count = stall_count + 1
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
                  if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                     stall_count = stall_count + 1
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
                  if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                     stall_count = stall_count + 1
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
               call steady_wind_with_element_diffusion(500, 1.0d0,      &
                                                       .true., j)
               if (j .eq. 0) then
                  ! newton_finished is itself a stop condition of the loop,
                  ! so no other flag is set: the run reports the JFNK finish
                  ! and not a du criterion it never had to satisfy.
                  du_plateaued    = .false.
                  newton_finished = .true.
               else
                  ! JFNK failed; it returned its best iterate. Do NOT accept
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

      ! Ensure the post-processing and final-output rates carry the full physics
      ! even if the loop hit the iteration cap before the stage flip fired.
      if (use_sec_ion) then
         if (.not. sec_ion_active) sec_ion_armed_step = count
         sec_ion_active = .true.
      endif
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
      else if (newton_finished) then
         write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
              '     -> converged: JFNK steady solution (||R|| <',             &
              resid_th, ' AND flux spread <', flux_spread_th, ')'
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

      call write_run_counter_report

      !---------------------------------------------------!

      ! Refresh the H(n=2) diagnostics from the converged state for
      ! write_excited_H below (the in-loop Balmer source was lagged one step).
      ! The Balmer photoelectric heating is ~1e-4 of the total, so we do NOT
      ! re-solve ioniz_eq here: that would desync T/f_sp from the converged
      ! hydro state for a negligible heat correction.
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)

      !---------------------------------------------------!

      ! Write final thermodynamic and ionization profiles
      write(*,*) '(EXHALE_main.f90) Writing final results to file..'
      call element_census_reservoir('output (final state)', rho, f_sp)
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
      
      ! Choose index to evaluate Mdot far enough from the top boundary  
      j = N - 20
      
      ! Evaluate steady state log of Mdot
      Mdot = log10(4.0*pi*rho(j)*v(j)*r(j)*r(j)*n0*v0*mu*R0*R0)
      
      ! Correct for the 2D approximation used
      if (appx_mth.eq.'Rate/2 + Mdot/2') Mdot = Mdot - log10(2.0)
      if (appx_mth.eq.'Mdot/4') Mdot = Mdot - log10(4.0)
      
      
      ! Write Mdot in output
      write(*,*) ' '
      write(*,*) '----- Results -----'
      write(*,*) ' '
      write(*,*) '---> 2D approximate method: ', appx_mth
      write(*,101) ' ---> Log10 of steady-state Mdot = ', Mdot, ' g/s'
      
      ! Write Mdot to report file 
      open(unit = outfile, file = 'EXHALE_setup.out', access = 'append' )
      	write(outfile,101) ' '
      	write(outfile,102) ' - Log10 of steady-state Mdot = ', Mdot, ' g/s'
      close(unit = outfile)
      
101   format (A35,F5.2,A4)
102   format (A32,F5.2,A4)

      !---------------------------------------------------!

      contains

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


      ! Cells whose semi-implicit energy update landed on the temperature
      ! floor (energy_semi_implicit). A cell on the floor has not converged
      ! to a physical temperature -- the update overshot and the clamp
      ! absorbed it -- and at 0.01 T0 the chemical network is frozen, so its
      ! composition is not constrained by its reaction balance either. Silent
      ! for a run that never reaches the floor, which is every matrix case.
      if (n_energy_floor_hits .gt. 0) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                                  &
            '     energy floor WARNING: ', n_energy_floor_hits,              &
            ' activation(s) in ', n_energy_floor_cells(),      &
            ' distinct cell(s), steps ', energy_floor_first_step, ' to ',    &
            energy_floor_last_step, ''
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

      ! Interfaces whose flux was dropped to first order to hold an RK stage
      ! inside rho > 0, rho e > 0 (positivity_limited_fluxes). Silent for a
      ! run that stayed inside the positivity bound of its CFL number, which
      ! is the run an unguarded build produces.
      if (n_faces_flux_positivity_limited .gt. 0) then
         write(*,'(A,I0,A)')                                                  &
            '     flux correction: ', n_faces_flux_positivity_limited,        &
            ' face flux(es) dropped to first order for positivity'
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

      end subroutine write_run_counter_report

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
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
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
         ! WHAT Apply_BC DOES TO THE INTERIOR is measured on the way past:
         ! it converts every cell to primitive variables, writes the ghosts,
         ! and converts back, and with a caloric EOS that round trip is not
         ! the identity. The review (B1) asks for the interior bytes to be
         ! preserved; this says how far they are not.
         u_umA = u
         call Apply_BC(u)
         upmap_bc_roundtrip = 0.0d0
         do j = 1, N
            do k = 1, 3
               if (u_umA(k,j) .ne. 0.0d0)                                  &
                  upmap_bc_roundtrip = max(upmap_bc_roundtrip,             &
                       abs(u(k,j)-u_umA(k,j))/abs(u_umA(k,j)))
            enddo
         enddo
         write(*,'(A,ES10.3)') ' [update map] Apply_BC round trip moves'// &
              ' the INTERIOR by, at most, a relative', upmap_bc_roundtrip
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
            call carrier_steady_residual(rho, v, T, f_sp, upmap_crmax,     &
                 upmap_cj, upmap_cic, res_out=res_um_all,                  &
                 terms_out=terms_um_all)
            R_carr_um = res_um_all(1:N,1)
            T_carr_um = terms_um_all(1:N,1)
         endif
         u_umSave    = u
         f_sp_umSave = f_sp
      else
         u    = u_umSave
         f_sp = f_sp_umSave
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,          &
                                    nheiii,nheiTR,nm,ne,n_tot)
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
           'Apply_BC round trip (relative, measured once)'
      do k = 1, 3
         write(*,'(A,A,A,9(ES11.3))') '   ', rname(k), ' ', rmax(k),       &
              gmax(k), dmax(k), smax(1,k), smax(2,k), smax(3,k),           &
              smax(4,k), smax(5,k),                                        &
              upmap_bc_roundtrip
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
         call carrier_steady_residual(rho, v, T, f_sp, crc_max,              &
                                      crc_j, crc_ic, rvol=crc_vol)
         write(*,'(A,ES10.3,A,ES10.3)')                                      &
              '          carrier row (evaluated here)             ',         &
              crc_vol, '   target "Carrier resid tol" ',                     &
              carrier_resid_th
      endif
      write(*,'(A)') '        (a state that has to be quoted as steady'//    &
           ' must pass all three; see Update_EXHALE.md section 145)'
      call write_residual_breakdown(Rres, u, heat, cool, 'marching stop')
      end subroutine report_marching_stop

      subroutine report_unmet_steady_gates
      ! Which of the two gates of steady_gates_met is the one still open, on
      ! the state the marching loop is stopping at. The steady solvers print
      ! the same statement when they return unaccepted; a marched state that
      ! stops short of the gates now says the same thing in the same terms.
      if (resid_max .ge. resid_th) then
         write(*,'(A,ES10.3,A,ES10.3)')                                      &
              '        the RESIDUAL gate is the one left: ||R||',            &
              resid_max, ' >=', resid_th
      else
         write(*,'(A,ES10.3,A,ES10.3)')                                      &
              '        the residual gate is met: ||R||',                     &
              resid_max, ' <', resid_th
      endif
      if (flux_spread_th .gt. 0.0d0) then
         if (fspread_gate .ge. flux_spread_th) then
            write(*,'(A,ES10.3,A,ES10.3)')                                   &
                 '        the FLUX gate is the one left: spread',            &
                 fspread_gate, ' >=', flux_spread_th
         else
            write(*,'(A,ES10.3,A,ES10.3)')                                   &
                 '        the flux gate is met: spread',                     &
                 fspread_gate, ' <', flux_spread_th
         endif
      endif
      end subroutine report_unmet_steady_gates


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
      if (abs(fspread_now - gate_fspread_accepted) .gt.                   &
          1.0d-10*max(abs(gate_fspread_accepted), 1.0d-30))               &
         write(*,'(A)') ' (EXHALE_main) WARNING: the written state is '// &
              'NOT the state the gates accepted -- something between the'//&
              ' steady solve and write_output changed u.'
      end subroutine assert_written_state_is_the_accepted_one

      ! ------------------------------------------------!

      subroutine steady_wind_with_element_diffusion(maxit, dtau0,       &
                                                   use_jfnk, jfnk_info)
      ! Steady wind that is self-consistent with the diffused element
      ! composition it carries.
      !
      ! The steady residual (steady_newton.f90) contains no diffusion: the
      ! composition is moved by the operator-split element_diffusion_step,
      ! so a single steady solve would freeze it at whatever state
      ! it was handed.  The two are therefore co-converged: solve the wind,
      ! relax the element composition to ITS steady state at that wind
      ! (relax_element_composition, which steps on the composition time scale
      ! and not on the hydro CFL step), repeat until the composition stops
      ! moving between passes.  The measure is absolute -- the largest change
      ! of the helium mass fraction over the pass divided by the base value --
      ! because the relative measure it replaces is meaningless in a cell the
      ! transport has emptied.  At most 20 outer passes.
      !
      ! The loop is damped.  Nothing in a Picard iteration of two solves --
      ! wind at fixed composition, composition at fixed wind -- keeps the two
      ! from chasing each other, and on the HD 209458 b Kzz = 0 wind they do:
      ! the drift settles into a limit cycle instead of falling.  The
      ! composition update is therefore under-relaxed,
      ! X <- X_old + omega (X_relaxed - X_old), starting at omega = 0.5 and
      ! halved (floor 0.125) on any pass whose drift failed to fall.  The
      ! drift tested here is the UNDAMPED distance to the fixed point, so a
      ! small omega cannot buy a false convergence.
      !
      ! With He_diffusion off the composition never moves and the body is
      ! one steady solve followed by the state refresh, which is what both
      ! call sites did before.
      !
      ! maxit / dtau0 / use_jfnk select the steady solver and its budget:
      ! the marching hand-off gives it 500 iterations from dtau0 = 1, the
      ! direct-steady route 3000 from the CFL dt, and either route may use
      ! JFNK or pseudo-transient continuation.
      integer, intent(in)  :: maxit
      real*8,  intent(in)  :: dtau0
      logical, intent(in)  :: use_jfnk
      integer, intent(out) :: jfnk_info

      comp_omega      = 0.5d0
      ! Diagnostic override (EXHALE_DIFF_OMEGA=<val>): pins the starting
      ! under-relaxation factor, so the undamped loop (1.0) can be compared
      ! against the damped one without a rebuild.
      call get_environment_variable('EXHALE_DIFF_OMEGA', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) comp_omega
      comp_drift      = 0.0d0
      comp_drift_prev = huge(1.0d0)

      call get_environment_variable('EXHALE_CARRIER_TRUST', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) carrier_trust
      call get_environment_variable('EXHALE_OUTER_PASSES', diag_env)
      if (len_trim(diag_env) .gt. 0) read(diag_env,*) outer_pass_cap
      do it_diff = 1, merge(outer_pass_cap, 1, he_diffusion .or.         &
                                  (thereis_mol .and. carrier_transport))
         if (use_jfnk) then
            ! THE COUPLED ROUTE (section 139). With the carrier row among
            ! the unknowns the outer loop is not an alternation at all: one
            ! solve returns a wind and a carrier partition that are steady
            ! states of each other, and the loop below runs once.
            call set_carrier_unknown(thereis_mol .and. carrier_transport &
                                     .and. carrier_in_newton)
            call solve_steady_jfnk(u, f_sp, resid_max, maxit, dtau0,    &
                                   40, jfnk_info)
         else
            call solve_steady_ptc(u, f_sp, resid_max, maxit, dtau0,     &
                                  jfnk_info)
         endif
         call set_carrier_unknown(.false.)
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:);  E = u(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,       &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
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
         if (jfnk_info .ne. 0) exit             ! steady solve failed
         if (.not. he_diffusion .and. .not. (thereis_mol .and.            &
             carrier_transport)) exit
         ! The coupled solve has already made the carriers part of the
         ! answer, so there is no carrier pass to take and no alternation to
         ! iterate. A run that also diffuses helium still needs its own
         ! relaxation, and falls through to it.
         if (thereis_mol .and. carrier_transport .and. carrier_in_newton  &
             .and. .not. he_diffusion) exit
         comp_drift = 0.0d0
         kd = 0
         ! Element composition relaxed to its steady state at the fixed wind
         if (he_diffusion) then
            call relax_element_composition(rho,v,T,f_sp,comp_omega,       &
                                           comp_drift,kd)
            call ioniz_eq(T,rho,f_sp,heat,cool,eta)
         endif
         ! The molecular carriers are a THIRD participant in the same Picard
         ! iteration, relaxed to their own steady state at the same fixed
         ! wind and under the same damping. The header of this routine
         ! records why the damping exists -- nothing in a Picard iteration
         ! of two solves keeps them from chasing each other -- and a third
         ! makes that risk larger, not smaller, which is why the drift
         ! reported below is the worst of the two composition drifts.
         if (thereis_mol .and. carrier_transport .and.                    &
             .not. carrier_in_newton) then
            call relax_photochemical_composition(rho,v,T,f_sp,           &
                                                 carrier_trust,          &
                                                 carrier_drift,kc)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,     &
                                       nheiii,nheiTR,nm,ne,n_tot)
            call comp_T_from_p(p,n_tot,ne,T)
            call ioniz_eq(T,rho,f_sp,heat,cool,eta)
            comp_drift = max(comp_drift, carrier_drift)
            kd = max(kd, kc)
         endif
         write(*,'(A,I0,A,I0,A,F6.3,A,ES10.2)') ' (EXHALE_main) '//      &
              'steady-wind diffusion outer pass ', it_diff, ': ', kd,   &
              ' relaxation steps, omega =', comp_omega,                 &
              ', composition drift =', comp_drift
         ! Is the state a steady state of the CARRIER equation, and where is
         ! it least so? ||R|| measures only the hydrodynamic half of this
         ! Picard fixed point; this is the other half, on the row scale that
         ! lives beside every other row scale (steady_residual.f90).
         if (thereis_mol .and. carrier_transport) then
            call carrier_steady_residual(rho, v, T, f_sp, crc_max,       &
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
              write(*,'(A,F8.4,A,F8.4,A,I0,A,ES9.2)')                    &
                   '    H2 front: x2=0.5 at r=', r(max(j50,1)),          &
                   '  x2=1e-2 at r=', r(max(j01,1)), '  (pass ', it_diff, &
                   ', trust', carrier_trust, ')'
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
         ! WHAT ENDS THE LOOP.  The composition drift is reported and is
         ! NOT the test, because it cannot be one: the transport operator's
         ! smallest step -- one cell crossing time -- already moves the
         ! composition by 1.1e-2 on the He/H = 0.0793 hot Uranus, so a 1e-3
         ! drift gate stands below anything the pass can produce and can
         ! only fire once the front has stopped moving. A gate that cannot
         ! fire reads as a control to the next person.
         !
         ! The test is the same one the wind is held to: is each row of the
         ! state steady on the scale of its own largest terms. For the
         ! carriers that is the volume-weighted carrier residual above. The
         ! helium element relaxation keeps its own drift test, which is a
         ! different measure of a different operator.
         if (he_diffusion .and. .not. (thereis_mol .and.                 &
             carrier_transport)) then
            if (comp_drift .lt. 1.0d-3) exit
         else if (thereis_mol .and. carrier_transport) then
            if (crc_vol .lt. carrier_resid_th .and.                     &
                (.not. he_diffusion .or. comp_drift .lt. 1.0d-3)) exit
         endif
         ! Damp harder whenever the drift failed to fall over the pass.
         ! Helium only: the carrier pass is bounded rather than damped.
         if (he_diffusion .and. comp_drift .ge. comp_drift_prev .and.    &
             comp_omega .gt. 0.125d0) then
            comp_omega = max(0.5d0*comp_omega, 0.125d0)
            write(*,'(A,F6.3)') '    -> composition drift did not '//    &
                 'fall; under-relaxation omega =', comp_omega
         endif
         comp_drift_prev = comp_drift
      enddo

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

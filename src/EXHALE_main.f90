      program Hydro_ioniz

      use global_parameters
      use species_table, only: n_mion, mion_fsp
      use Read_input
      use Initialization
      use setup_report
      use eval_time_step
      use energy_semi_implicit
      use utils
      use composition, only: get_species_densities, comp_T_from_p,         &
                             comp_p_from_T, element_ratio_HeH
      use binary_element_diffusion, only: element_diffusion_step,          &
                                          relax_element_composition
      use lower_column, only: lower_column_solve
      use steady_residual_mod, only: assemble_residual, residual_norms, residual_norms_vol
      use viscous_conduction, only: transport_active, viscous_conduction_step
      use steady_newton, only: neq_newton, pack_U, unpack_U, newton_residual, &
                               eval_residual, frozen_residual,              &
                               build_banded_jac, band_matvec,               &
                               kl_jac, ku_jac, solve_steady_ptc,           &
                               solve_steady_jfnk, set_base_fix
      use Conversion
      use ionization_equilibrium
      use utils_ion_eq, only: write_cool_breakdown_eq, write_heat_breakdown_eq
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
      
      ! Logical variables
      logical :: l_isnan = .false.
      logical :: in_plm_stage = .false.   ! true while in the stage-1 (PLM) phase of a two-stage PLM->WENO3 run

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
      logical :: du_newton_armed = .false.   ! JFNK hand-off, and the secondary-
                                             ! ionization flip sharing its test

      ! Residual-based convergence monitor (Resid tol option)
      real*8, dimension(:,:), allocatable :: Rres
      real*8  :: resid_c(3), resid_cv(3), resid_max, flux_spread
      logical :: is_resid_ok

      ! Newton-residual self-test scratch (EXHALE_NEWTON_TEST hook)
      real*8, allocatable :: Yvec(:), Fvec(:)
      real*8, dimension(:,:), allocatable :: f_sp_test
      ! Banded-Jacobian self-test scratch (EXHALE_JAC_TEST hook)
      real*8, allocatable :: abjac(:,:), rdir(:), Jr(:), dFD(:), F0f(:)
      real*8, dimension(:), allocatable :: heat0, cool0, npart0
      real*8 :: jac_eps, jac_err
      ! --- lightweight phase profiler (gated by env EXHALE_PROFILE=1) ---
      logical :: do_profile = .false.
      real*8  :: tp_step0, tp_a, tp_ion = 0.0d0, tp_hyd = 0.0d0, tp_tot = 0.0d0
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
      integer :: it_diff, kd

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
      call allocate_state_vectors

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
         rho = W(1,:)
         v   = W(2,:)
         p   = W(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
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
      ! max_j |R(j,k)| / max|u(:,k)| over the wind region [j_min:N]. Used to
      ! rank candidate "converged" states (premature-dip vs true steady).
      ! Evaluated in WENO3 (the production scheme the states converged under).
      call get_environment_variable('EXHALE_RESIDUAL', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
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
         call residual_norms_vol(Rres, u, resid_cv)
         write(*,'(A)') ' (EXHALE_main) EXHALE_RESIDUAL=1 steady residual ||R||:'
         write(*,'(A)') '   component    L-inf: max|R|/max|u|    vol-wt: '// &
                        'sum|R|V/sum|u|V   [1/t_s]'
         do k = 1,3
            write(*,'(A,I2,4X,ES16.6,4X,ES16.6)') '   k=', k,               &
                 resid_c(k), resid_cv(k)
         enddo
         write(*,'(A,ES12.4,A,ES12.4)') '   ||R|| = max_k :  L-inf =',      &
              maxval(resid_c), '   vol-weighted =', maxval(resid_cv)
         write(*,*) '(EXHALE_main) EXHALE_RESIDUAL=1: residual reported, stopping.'
         stop
      endif

      ! Optional Newton-residual self-test (env EXHALE_NEWTON_TEST=1): verify
      ! the vector residual F(Y) used by the steady solver reproduces the
      ! diagnostic residual. (a) pack/unpack are exact inverses; (b) the
      ! the max|F|/max|u| for each component over [j_min:N] equals the EXHALE_RESIDUAL
      ! values. Validates increment (ii)-2 before the Jacobian/PTC driver.
      call get_environment_variable('EXHALE_NEWTON_TEST', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         allocate(Yvec(neq_newton()), Fvec(neq_newton()))
         call pack_U(u, Yvec)
         W = 0.0d0                                  ! reuse W as unpack target
         call unpack_U(Yvec, W)
         write(*,'(A,ES12.3)') ' (newton_test) max|unpack(pack(u))-u| (cells 1..N) = ', &
              maxval(abs(W(:,1:N) - u(:,1:N)))
         f_sp_test = f_sp                           ! copy: newton_residual mutates it
         call newton_residual(Yvec, f_sp_test, Fvec)
         write(*,'(A)') ' (newton_test) max|F|/max|u| over [j_min:N] (cf. EXHALE_RESIDUAL):'
         do k = 1,3
            mom = 0.0d0
            do j = 1,N
               mom(j) = Fvec(3*(j-1)+k)
            enddo
            write(*,'(A,I2,4X,ES16.6)') '   k=', k,                    &
               maxval(abs(mom(j_min:N)))/max(maxval(abs(u(k,j_min:N))),1.0d-30)
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
         k = neq_newton()
         allocate(Yvec(k), Fvec(k), F0f(k), rdir(k), Jr(k), dFD(k))
         allocate(abjac(2*kl_jac+ku_jac+1, k))
         call pack_U(u, Yvec)
         f_sp_test = f_sp
         ! F0 + frozen heat/cool at Y, then the banded Jacobian
         call eval_residual(Yvec, f_sp_test, Fvec, heat0, cool0, npart0)
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
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         sec_ion_active = use_sec_ion   ! bypass mode skips the time loop
         call U_to_W(u,W)
         call eval_dt(W, dt, dt_loc)            ! CFL dt = default PTC dtau0
         resid_max = resid_th
         if (resid_max .le. 0.0d0) resid_max = 1.0d-3
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
         call write_output(rho,v,p,T,heat,cool,eta,                    &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
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

      is_mom_const = .false.
      is_zero_dt   = .false.
      is_stalled   = .false.
      du_prev      = huge(1.0d0)
      stall_count  = 0

      ! Staged secondary ionization: the SvS85 coupling starts applied only if
      ! "Secondary_ionization: Immediate" was set; otherwise it is switched on
      ! after the wind first converges without it (see the stage flip below).
      sec_ion_active = (use_sec_ion .and. sec_ion_immediate)
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
      resid_c   = huge(1.0d0)
      resid_max = huge(1.0d0)
      is_resid_ok = .false.

      ! Activate base-cell startup diagnostic if requested
      call get_environment_variable('EXHALE_DIAG_BASE', diag_env)
      if (trim(diag_env) .eq. '1') then
         diag_base = .true.
         open(unit=778, file='output/base_diag.txt', status='replace')
         write(778,'(A)') '# count  j  r[Rp]  n[cm-3]  v[cm/s]  T[K]  heat  cool'
      endif

      !------ Main temporal loop ------!
      ! Stop on any of: momentum constant (du<du_th), steady state
      ! (dtu<dtu_th), du plateau (stall), or the hard iteration cap.
      do while( ( .not.is_mom_const .and. .not.is_zero_dt .and.            &
                  .not.is_stalled   .and. count < count_max ) .or.         &
                force_start )

            !---- Time step evaluation ----!
            ! dt_loc = cell-by-cell pseudo-dt ("Time stepping: Local"), or
            ! uniformly the global dt (default; bit-identical updates).
            call eval_dt(W,dt,dt_loc)
            if (do_profile) tp_step0 = omp_get_wtime()

            !-------------------------------------------------!

            !--- Thermodynamic evolution ---!

            ! Save previous step solution
            u_old = u
            
            ! FIRST RK STEP
            
            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u,WL,WR) 
            
            ! Evaluate flux difference and source terms
            call RK_rhs(u,WL,WR,alpha,dF,S)
            
            do k = 1,3
               u1(k,:) = u(k,:) - dt_loc*(dF(k,:) - S(k,:))
            enddo
             
            ! Apply boundary conditions
            call Apply_BC(u1)
                      
            !----------------------------
            
            ! SECOND RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u1,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u1,WL,WR,alpha,dF,S)
                    
            ! Advance in time
            do k = 1,3
               u2(k,:) = (3.0*u(k,:) + u1(k,:) - dt_loc*(dF(k,:) - S(k,:)))/4.0
            enddo

            ! Apply boundary conditions            
            call Apply_BC(u2)
            
            !----------------------------
            
            ! THIRD RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u2,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u2,WL,WR,alpha,dF,S)
                  
            ! Advance in time
            do k = 1,3
               u(k,:) = (u(k,:) + 2.0*(u2(k,:) - dt_loc*(dF(k,:) - S(k,:))))/3.0
            enddo

            ! Apply boundary conditions
            call Apply_BC(u)
 		
		!------------------------------------------------!
 		
 		!---- Ionization Equilibrium ----!
 		
            ! Extract primitive variables
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
  	
  	      !------------------------------------------------!
  	      
  	      ! Evolve solution in time due to source terms
            !     using a semi-implicit step (default) or the original
            !     explicit forward-Euler update ("Energy solver: Explicit")

            if (use_semi_implicit_energy) then
               call solve_energy_semi_implicit(u,W,dt_loc,heat,cool,f_sp)
            else
               u(3,:) = u(3,:) + dt_loc*(heat - cool)
            endif

		call Apply_BC(u)

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

            ! Periodic Shapiro low-pass filter to damp the gravity-unbalanced
            ! sound waves (base breathing), as in CETIMB (Koskinen et al. 2013a).
            if (shapiro_eps .gt. 0.0d0 .and.                             &
                mod(count, shapiro_every) .eq. 0) then
               call shapiro_filter(u)
               call Apply_BC(u)
            endif

            !------------------------------------------------!

            ! Convert to physical variables and extract profiles
            call U_to_W(u,W)
            rho = W(1,:)
            v   = W(2,:)
            p   = W(3,:)
            E   = u(3,:)
            
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

            ! CETIMB-style base velocity: update the mass-flux constant F_c from
            ! the [j_min:N] constant-momentum (escape) region -- NOT the base,
            ! where rho*v*r^2 is not yet flat. A slow exponential moving average
            ! makes v0 track the STABLE flux constant rather than the step-by-step
            ! transient noise (the instantaneous value feeds the breathing back
            ! into the base velocity). Used by the next step's base BC.
            if (base_v_massflux .and. j_min .le. N) then
               dum = sum(mom(j_min:N))/dble(N - j_min + 1)
               if (base_flux_const .le. 0.0d0) then
                  base_flux_const = dum
               else
                  base_flux_const = 0.99d0*base_flux_const + 0.01d0*dum
               endif
            endif
		
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
            ! ||R|| (volume-weighted by default, resid_vol) is computed and
            ! reported FOR REFERENCE ONLY; it gates the stop solely when the user
            ! explicitly requests it via "Resid tol:" (resid_th>0).
            if (mod(count, N_resid) .eq. 0) then
               call assemble_residual(u, n_tot + ne, heat, cool, Rres)
               if (resid_vol) then
                  call residual_norms_vol(Rres, u, resid_c)
               else
                  call residual_norms(Rres, u, resid_c)
               endif
               resid_max   = maxval(resid_c)
               flux_spread = (maxval(mom(j_min:N)) - minval(mom(j_min:N)))   &
                    /max(abs(sum(mom(j_min:N))/dble(N - j_min + 1)), 1.0d-30)
               write(*,'(A,I0,A,ES10.3,A,ES10.3)') '   [diag] step ', count, &
                    '  flux rho*v*r^2 spread=', flux_spread,                  &
                    '   ||R||(ref)=', resid_max
            endif
            is_resid_ok = (resid_max .lt. resid_th)

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
                  rec_method   = 'WENO3'
                  ! Flip the discretization flags too (Source/Num_Fluxes/RK_rhs/
                  ! Apply_BC read these, not rec_method); without this the
                  ! WENO3 stage ran with PLM-style pressure flux/source terms.
                  use_plm      = .false.
                  use_weno3    = .true.
                  in_plm_stage = .false.
                  stall_count  = 0
                  du_prev      = huge(1.0d0)
                  write(*,'(A,ES12.4,A,I0)') '    -> switched PLM -> WENO3 at du =', &
                                             du, ', step ', count
               endif
               is_mom_const = .false.
               is_zero_dt   = .false.
               is_stalled   = .false.
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
               is_mom_const = .false.
               is_zero_dt   = .false.
               is_stalled   = .false.
            else if (resid_th .gt. 0.0d0) then
               ! Stage 2, residual-based convergence: trust the steady
               ! residual ||R|| instead of du (which is blind to the
               ! operator-split energy imbalance). du/dtu kept for logging.
               is_mom_const = is_resid_ok
               is_zero_dt   = .false.
               is_stalled   = .false.
            else
               ! Stage 2 (WENO3) or single-stage: normal convergence + stall,
               ! each additionally gated on mass-flux level stability.
               is_mom_const = (du .lt. du_th)   .and. du_stop_armed          &
                                                .and. is_level_stable
               is_zero_dt   = (dtu .lt. dtu_th) .and. is_level_stable
               ! Stall detection: du settled on a plateau
               if (count .gt. 1) then
                  if (abs(du - du_prev)/max(du,1.0d-30) .lt. stall_tol) then
                     stall_count = stall_count + 1
                  else
                     stall_count = 0
                  endif
               endif
               du_prev    = du
               is_stalled = (stall_count .ge. N_stall) .and. is_level_stable
            endif

            ! Staged secondary ionization: the SvS85 coupling is applied only
            ! after the wind has first converged without it -- from a cold IC the
            ! secondary-ionization base feedback amplifies the startup transient
            ! into a runaway (NaN); from a converged state the coupling is benign
            ! and the run re-converges with the full physics included.
            if (use_sec_ion .and. .not.sec_ion_active .and.                   &
                (is_mom_const .or. is_zero_dt .or. is_stalled)) then
               sec_ion_active = .true.
               sec_flip_step  = count
               is_mom_const = .false.;  is_zero_dt = .false.;  is_stalled = .false.
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
               is_mom_const = .false.;  is_zero_dt = .false.;  is_stalled = .false.
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
                  is_mom_const = .false.;  is_zero_dt = .false.
                  is_stalled   = .false.
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
               if (valve_eps .le. 0.0d0) then
                  valve_eps = 1.0d-4
                  write(*,'(A)') ' (EXHALE_main) Solver: Newton -> '//     &
                       'adopting smooth base valve, eps = 1.0e-4'
               endif
               resid_max = resid_th
               if (resid_max .le. 0.0d0) resid_max = 1.0d-3
               write(*,'(A,I0,A,ES10.2)') ' (EXHALE_main) Newton finish '// &
                    'at step ', count, ', target ||R|| <', resid_max
               call steady_wind_with_element_diffusion(500, 1.0d0,      &
                                                       .true., j)
               if (j .eq. 0) then
                  is_mom_const = .true.       ! exit the marching loop
                  is_stalled   = .false.
               else
                  ! JFNK failed; it returned its best iterate. Do NOT
                  ! accept an unconverged state as the answer -- fall back
                  ! to plain time-marching (du-based stops re-armed).
                  use_newton_solver = .false.
                  write(*,'(A,I0,A)') ' (EXHALE_main) JFNK failed (info=', &
                       j, '); resuming time-marching with du-based stops'
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
                  call write_output(rho,v,p,T,heat,cool,eta,             &
                                    nhi,nhii,nhei,nheii,nheiii,          &
                                    nheiTR,nm,'eq')
                                   
            endif     
            
            ! Exit from temporal loop if only post processing has to be done
	      if (do_only_pp) exit

            ! Force continue for the first 1000 loops if force_start is enabled
            if (force_start) force_start = count .le. 1000

            ! Deterministic step cap for serial-vs-parallel verification.
            if (max_steps .gt. 0 .and. count .ge. max_steps) exit

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
      if (use_sec_ion) sec_ion_active = .true.
      if (do_profile) then
         write(*,'(A)')          ' (profile) phase breakdown over the run:'
         write(*,'(A,F10.3,A)')  '   ioniz_eq total :', tp_ion, ' s'
         write(*,'(A,F10.3,A)')  '   full step total:', tp_tot, ' s'
         if (tp_tot .gt. 0.0d0) write(*,'(A,F6.1,A)')                        &
            '   ioniz_eq fraction:', 100.0d0*tp_ion/tp_tot, ' %'
      endif

      ! Report which criterion stopped the loop
      if (is_mom_const) then
         write(*,*) '    -> converged: momentum constant (du < du_th)'
      else if (is_zero_dt) then
         write(*,*) '    -> converged: steady state reached (dtu < dtu_th)'
      else if (is_stalled) then
         write(*,'(A,I0,A)') '     -> stopped: du plateau detected (stalled ', &
                             stall_count, ' steps); steady state assumed'
      else if (count .ge. count_max) then
         write(*,'(A,I0,A)') '     -> stopped: reached count_max = ',          &
                             count_max, ' iterations WITHOUT convergence'
      endif
      write(*,'(A,I0,A,ES11.4,A,ES11.4)') '     final: count=', count,         &
                             '  du=', du, '  dtu=', dtu

      ! Task 2: report how often the analytic-Jacobian Newton handled the
      ! ionization-equilibrium solve vs. fell back to MINPACK hybrd1.
      if (nt_calls .gt. 0) then
         write(*,'(A,I0,A,I0,A,F6.2,A)')                                       &
            '     ioniz-eq solver: Newton ', nt_calls - nt_fallback,           &
            ' / ', nt_calls, ' solves (',                                      &
            100.0d0*dble(nt_calls-nt_fallback)/dble(nt_calls),                 &
            '%); rest fell back to hybrd1'
      endif

      ! Root validation of the atomic ionization solves: how often a cell was
      ! re-solved from another starting point, how often the first root lay
      ! outside the physical simplex, how often a stored state was rejected as
      ! a starting point, and how often no starting point produced an
      ! admissible root. Silent for a run that never leaves the simplex.
      if (ieq_n_retry + ieq_n_unphys + ieq_n_reseed + ieq_n_noroot .gt. 0) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                                    &
            '     ioniz-eq roots: ', ieq_n_retry, ' cell solve(s) restarted, ',&
            ieq_n_unphys, ' root(s) outside the physical simplex, ',           &
            ieq_n_reseed, ' stored state(s) rejected, ', ieq_n_noroot,         &
            ' cell(s) left on the ionization balance'
      endif

      ! Molecular cells whose equilibrium roots all left the physical simplex,
      ! so the closest one was clamped onto the element budget. Silent for a
      ! run whose molecular solve stays inside the simplex, and for an
      ! atomic run.
      if (ieq_n_mol_clamped .gt. 0) then
         write(*,'(A,I0,A)')                                                  &
            '     ioniz-eq molecular: ', ieq_n_mol_clamped,                   &
            ' cell(s) clamped onto the element budget'
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
      call write_output(rho,v,p,T,heat,cool,eta,                         &
                        nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')

      ! Diagnostic: radiative cooling in each channel vs radius
      ! (reuses eval_cool's coefficients; see Huang et al. 2023 Fig. 10).
      call write_cool_breakdown_eq(T,rho,f_sp)

      ! Diagnostic: volumetric heating in each channel vs radius
      ! (photoionization per absorber + Balmer + He-recomb + Penning).
      call write_heat_breakdown_eq(T,rho,f_sp)

      ! Diagnostic: H(n=2) populations, Balmer proton source, and
      ! photoelectric/de-excitation heating vs radius (Huang Figs. 11/27/10/26).
      if (use_excited_H) call write_excited_H

      !---------------------------------------------------!

      ! Post processing to include advection
      write(*,*) '(EXHALE_main.f90) Starting the post processing routine..'

      call post_process_adv(rho,v,p,T,heat,cool,eta,    &
                            nhi,nhii,nhei,nheii,nheiii,nheiTR,nm)
      
      write(*,*) '(EXHALE_main.f90) Post processing routine done.'                           
      
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

      do it_diff = 1, merge(20, 1, he_diffusion)
         if (use_jfnk) then
            call solve_steady_jfnk(u, f_sp, resid_max, maxit, dtau0,    &
                                   40, jfnk_info)
         else
            call solve_steady_ptc(u, f_sp, resid_max, maxit, dtau0,     &
                                  jfnk_info)
         endif
         call U_to_W(u,W)
         rho = W(1,:);  v = W(2,:);  p = W(3,:);  E = u(3,:)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,       &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,       &
                                    nheiii,nheiTR,nm,ne,n_tot)
         if (jfnk_info .ne. 0) exit             ! steady solve failed
         if (.not. he_diffusion) exit
         ! Element composition relaxed to its steady state at the fixed wind
         call relax_element_composition(rho,v,T,f_sp,comp_omega,          &
                                        comp_drift,kd)
         call ioniz_eq(T,rho,f_sp,heat,cool,eta)
         write(*,'(A,I0,A,I0,A,F6.3,A,ES10.2)') ' (EXHALE_main) '//      &
              'steady-wind diffusion outer pass ', it_diff, ': ', kd,   &
              ' relaxation steps, omega =', comp_omega,                 &
              ', composition drift =', comp_drift
         call get_environment_variable('EXHALE_DIFFUSION_CHECK', diag_env)
         if (trim(diag_env) .eq. '1')                                    &
            call write_diffusion_pass_profile(it_diff)
         if (comp_drift .lt. 1.0d-3) exit
         ! Damp harder whenever the drift failed to fall over the pass.
         if (comp_drift .ge. comp_drift_prev .and.                       &
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
      integer :: j, uu
      logical :: first

      heh_l = element_ratio_HeH(f_sp)
      uu    = 773
      first = (pass .eq. 1)
      if (first) then
         open(unit=uu, file='output/diffusion_pass_profiles.txt',          &
              status='replace')
         write(uu,'(A)') '# outer-pass composition profiles '//            &
              '(EXHALE_DIFFUSION_CHECK=1)'
         write(uu,'(A)') '# columns: pass  r[Rp]  (He/H)/HeH  '//          &
              'r^2 rho v [n0 mH cm/s Rp^2]  T[K]'
      else
         open(unit=uu, file='output/diffusion_pass_profiles.txt',          &
              status='old', position='append')
      endif
      do j = 1, N
         write(uu,'(I5,4ES16.7)') pass, r(j), heh_l(j)/HeH,                &
              r(j)**2*rho(j)*v(j), T(j)*T0
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
      end program Hydro_ioniz

      program Hydro_ioniz

      use global_parameters
      use species_table, only: n_mion, mion_fsp
      use Read_input
      use Initialization
      use setup_report
      use eval_time_step
      use energy_semi_implicit
      use utils
      use composition, only: get_species_densities, comp_T_from_p, comp_p_from_T
      use steady_residual_mod, only: assemble_residual, residual_norms
      use steady_newton, only: neq_newton, pack_U, unpack_U, newton_residual, &
                               eval_residual, frozen_residual,              &
                               build_banded_jac, band_matvec,               &
                               kl_jac, ku_jac, solve_steady_ptc,           &
                               solve_steady_jfnk, set_base_fix
      use Conversion
      use ionization_equilibrium
      use utils_ion_eq, only: write_cool_breakdown_eq
      use excited_hydrogen, only: excited_H_update, write_excited_H
      use Reconstruction_step
      use RK_integration
      use BC_Apply
      use output_write
      use post_processing
      use ionization_equilibrium
      use newton_solver, only: nt_calls, nt_fallback   ! Task 2 usage counters
      
      implicit none
      
      ! Logical variables
      logical :: l_isnan = .false.
      logical :: in_plm_stage = .false.   ! true while in the stage-1 (PLM) phase of a two-stage PLM->WENO3 run

      ! Mass-flux level-stability tracking (ring buffer over N_stall steps)
      real*8, allocatable :: lev_hist(:)
      real*8  :: lev_now, lev_rel
      integer :: lev_count
      logical :: is_level_stable

      ! Residual-based convergence monitor (Resid tol option)
      real*8, dimension(1-Ng:N+Ng,3) :: Rres
      real*8  :: resid_c(3), resid_max
      logical :: is_resid_ok

      ! Newton-residual self-test scratch (ATES_NEWTON_TEST hook)
      real*8, allocatable :: Yvec(:), Fvec(:)
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_sp_test
      ! Banded-Jacobian self-test scratch (ATES_JAC_TEST hook)
      real*8, allocatable :: abjac(:,:), rdir(:), Jr(:), dFD(:), F0f(:)
      real*8, dimension(1-Ng:N+Ng) :: heat0, cool0
      real*8 :: jac_eps, jac_err
      
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
      real*8, dimension(1-Ng:N+Ng) :: mom
      real*8 :: mom_max,mom_min

      ! Convergence stall detection (ported from ATES_extended)
      real*8  :: du_prev
      integer :: stall_count

      ! Optional base-cell startup diagnostic (env ATES_DIAG_BASE=1): dumps the
      ! first few cells (r, n, v, T, heat, cool, dC/dT sign) for the first steps,
      ! to trace IC-startup transients (e.g. warm-Parker base breakdown). Off by
      ! default => zero effect on normal runs.
      logical            :: diag_base = .false.
      character(len=64)  :: diag_env

      ! Phase 3a decoupled outer iteration over the excited-H (H n=2) feedback
      real*8  :: exc_rel

      ! Maximum eigenvalue      
      real*8 :: alpha

      ! Temporal step (global) and per-cell pseudo-time steps
      real*8 :: dt
      real*8, dimension(1-Ng:N+Ng) :: dt_loc
      
      ! Mdot value
      real*8 :: Mdot
      
      ! Vectors of thermodynamical variables
      real*8, dimension(1-Ng:N+Ng) :: rho,v,E,p,T,cs
      real*8, dimension(1-Ng:N+Ng) :: heat,cool
      real*8, dimension(1-Ng:N+Ng) :: eta 
      real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii
      real*8, dimension(1-Ng:N+Ng) :: nheiTR
      ! Metal ion densities, 2D: column i = ion i of the species_table
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng) :: ne,n_tot
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_sp  ! 1-6: H/He(+HeITR), 7-9: CI/II/III, 10-12: OI/II/III, 13-15: NI/II/III, 16-18: MgI/II/III
       
      ! Conservative and primitive vectors
      real*8, dimension(1-Ng:N+Ng,3) :: u,u1,u2,u_old
      real*8, dimension(1-Ng:N+Ng,3) :: W,WL,WR
          
      ! Flux and source vectors
      real*8, dimension(1-Ng:N+Ng,3) :: dF,S
      
      !------------------------------------------------! 
      
      ! Open output report file 
      open(unit = outfile, file = 'ATES.out')
      
      !------------------------------------------------!
      
      ! Read planetary parameters from file
      call input_read
      
      !------------------------------------------------!
      
      ! Initialize simulations
      call init(W,u,f_sp)

      ! Optional IC-dump hook (env ATES_DUMP_IC=1): write the state exactly
      ! as initialized/loaded and stop. Lets the restart round-trip test
      ! inspect what load_IC restored BEFORE the first ionization-
      ! equilibrium solve re-equilibrates the species fractions.
      call get_environment_variable('ATES_DUMP_IC', diag_env)
      if (trim(diag_env) .eq. '1') then
         rho = W(:,1)
         v   = W(:,2)
         p   = W(:,3)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         heat = 0.0d0; cool = 0.0d0; eta = 0.0d0
         call write_output(rho,v,p,T,heat,cool,eta,                    &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
         write(*,*) '(ATES_main) ATES_DUMP_IC=1: IC state written, stopping.'
         stop
      endif

      ! Optional steady-residual diagnostic (env ATES_RESIDUAL=1): evaluate
      ! the finite-volume steady residual R = du/dt of the loaded state and
      ! stop. R = dF - S for mass/momentum and dF_E - S_E - (heat - cool) for
      ! energy; at a true steady state R = 0, INDEPENDENT of how the run was
      ! stopped (du/dtu/stall). Reported as the max relative residual rate
      ! max_j |R(j,k)| / max|u(:,k)| over the wind region [j_min:N]. Used to
      ! rank candidate "converged" states (premature-dip vs true steady).
      ! Evaluated in WENO3 (the production scheme the states converged under).
      call get_environment_variable('ATES_RESIDUAL', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         ! Apply_BC first so the residual depends only on the interior state
         ! (ghosts set by BC), matching the Newton residual F(Y) definition.
         call Apply_BC(u,u)
         call U_to_W(u,W)
         rho = W(:,1);  v = W(:,2);  p = W(:,3)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
         call assemble_residual(u, heat, cool, Rres)
         call residual_norms(Rres, u, resid_c)
         write(*,'(A)') ' (ATES_main) ATES_RESIDUAL=1 steady residual ||R||:'
         write(*,'(A)') '   component   max|R|/max|u| [1/t_s]'
         do k = 1,3
            write(*,'(A,I2,4X,ES16.6)') '   k=', k, resid_c(k)
         enddo
         write(*,*) '(ATES_main) ATES_RESIDUAL=1: residual reported, stopping.'
         stop
      endif

      ! Optional Newton-residual self-test (env ATES_NEWTON_TEST=1): verify
      ! the vector residual F(Y) used by the steady solver reproduces the
      ! diagnostic residual. (a) pack/unpack are exact inverses; (b) the
      ! per-component max|F|/max|u| over [j_min:N] equals the ATES_RESIDUAL
      ! values. Validates increment (ii)-2 before the Jacobian/PTC driver.
      call get_environment_variable('ATES_NEWTON_TEST', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         allocate(Yvec(neq_newton()), Fvec(neq_newton()))
         call pack_U(u, Yvec)
         W = 0.0d0                                  ! reuse W as unpack target
         call unpack_U(Yvec, W)
         write(*,'(A,ES12.3)') ' (newton_test) max|unpack(pack(u))-u| (cells 1..N) = ', &
              maxval(abs(W(1:N,:) - u(1:N,:)))
         f_sp_test = f_sp                           ! copy: newton_residual mutates it
         call newton_residual(Yvec, f_sp_test, Fvec)
         write(*,'(A)') ' (newton_test) max|F|/max|u| over [j_min:N] (cf. ATES_RESIDUAL):'
         do k = 1,3
            mom = 0.0d0
            do j = 1,N
               mom(j) = Fvec(3*(j-1)+k)
            enddo
            write(*,'(A,I2,4X,ES16.6)') '   k=', k,                    &
               maxval(abs(mom(j_min:N)))/max(maxval(abs(u(j_min:N,k))),1.0d-30)
         enddo
         write(*,*) '(ATES_main) ATES_NEWTON_TEST=1: done, stopping.'
         stop
      endif

      ! Optional banded-Jacobian self-test (env ATES_JAC_TEST=1): verify the
      ! colored-FD banded Jacobian of the frozen residual reproduces a
      ! directional finite difference, J*r ~= (F(Y+eps r)-F(Y))/eps. A wrong
      ! band layout gives an O(1) mismatch; a correct one matches to ~1e-6.
      call get_environment_variable('ATES_JAC_TEST', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         k = neq_newton()
         allocate(Yvec(k), Fvec(k), F0f(k), rdir(k), Jr(k), dFD(k))
         allocate(abjac(2*kl_jac+ku_jac+1, k))
         call pack_U(u, Yvec)
         f_sp_test = f_sp
         ! F0 + frozen heat/cool at Y, then the banded Jacobian
         call eval_residual(Yvec, f_sp_test, Fvec, heat0, cool0)
         call frozen_residual(Yvec, heat0, cool0, F0f)
         call build_banded_jac(Yvec, heat0, cool0, abjac)
         ! Deterministic probe direction (max|r| = 1)
         do j = 1, k
            rdir(j) = sin(0.1d0*dble(j))
         enddo
         rdir = rdir/maxval(abs(rdir))
         call band_matvec(abjac, rdir, Jr)
         jac_eps = 1.0d-7
         Yvec = Yvec + jac_eps*rdir
         call frozen_residual(Yvec, heat0, cool0, dFD)
         dFD = (dFD - F0f)/jac_eps
         jac_err = maxval(abs(dFD - Jr))/max(maxval(abs(Jr)), 1.0d-30)
         write(*,'(A,ES12.4)') ' (jac_test) max|F_frozen(Y)| = ', maxval(abs(F0f))
         write(*,'(A,ES12.4)') ' (jac_test) ||J*r - dFD||_inf / ||J*r||_inf = ', jac_err
         write(*,'(A)')        '   (correct band layout => ~1e-6; wrong => O(1))'
         write(*,*) '(ATES_main) ATES_JAC_TEST=1: done, stopping.'
         stop
      endif

      ! Optional steady-state PTC-Newton solve (env ATES_PTC=1): solve
      ! F(Y)=0 directly from the current IC, write the converged profiles,
      ! and stop. Validation of increment (ii)-4 against the marching
      ! reference (WASP-121b ~13.71).
      call get_environment_variable('ATES_PTC', diag_env)
      if (trim(diag_env) .eq. '1') then
         rec_method = 'WENO3';  use_weno3 = .true.;  use_plm = .false.
         call U_to_W(u,W)
         call eval_dt(W, dt, dt_loc)            ! CFL dt = default PTC dtau0
         resid_max = resid_th
         if (resid_max .le. 0.0d0) resid_max = 1.0d-3
         ! Optional dtau0 override for experimentation (ATES_PTC_DTAU0=<val>):
         ! a larger start probes the Newton regime directly.
         call get_environment_variable('ATES_PTC_DTAU0', diag_env)
         if (len_trim(diag_env) .gt. 0) read(diag_env,*) dt
         write(*,'(A,ES10.2)') ' (ATES_main) PTC dtau0 = ', dt
         ! Optional frozen-base experiment (ATES_PTC_NFIX=<n>): anchor the
         ! first n physical cells (non-smooth base-BC region) and solve for
         ! the wind on top of them.
         call get_environment_variable('ATES_PTC_NFIX', diag_env)
         if (len_trim(diag_env) .gt. 0) then
            read(diag_env,*) k
            allocate(Yvec(neq_newton()))
            call pack_U(u, Yvec)
            call set_base_fix(Yvec, k)
            deallocate(Yvec)
            write(*,'(A,I0,A)') ' (ATES_main) frozen-base: first ', k, ' cells anchored'
         endif
         call get_environment_variable('ATES_PTC_JFNK', diag_env)
         if (trim(diag_env) .eq. '1') then
            call solve_steady_jfnk(u, f_sp, resid_max, 3000, dt, 40, j)
         else
            call solve_steady_ptc(u, f_sp, resid_max, 3000, dt, j)
         endif
         ! Final consistent state + output
         call U_to_W(u,W);  rho = W(:,1);  v = W(:,2);  p = W(:,3)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call comp_T_from_p(p,n_tot,ne,T)
         if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
         call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
         call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,      &
                                    nheiii,nheiTR,nm,ne,n_tot)
         call write_output(rho,v,p,T,heat,cool,eta,                    &
                           nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')
         write(*,*) '(ATES_main) ATES_PTC=1: solver done, output written, stopping.'
         stop
      endif

      !------------------------------------------------!

      ! Generate report of the current setup
      call write_setup_report

	!---------------------------------------------------!

	! Close outfile
	close(unit = outfile)

      !------------------------------------------------!

      !---- Start computation ----!
      
      ! Get starting time
      write(*,*) '(ATES_main.f90) Starting time integration..'
      start = omp_get_wtime()

      ! Phase 3a excited-H feedback is a decoupled (lagged-explicit) source:
      ! at the top of every timestep the H(n=2) Balmer proton source + photo-
      ! electric heating are recomputed from the current state and frozen into
      ! the global arrays that ioniz_eq reads (never inside its Newton solve).
      ! The single hydro+ionization relaxation then converges hydro, ionization
      ! and the Balmer feedback together. Off => arrays stay zero, so the run is
      ! byte-identical to the Phase-2 result.
      exc_rel = 1.0d0

      is_mom_const = .false.
      is_zero_dt   = .false.
      is_stalled   = .false.
      du_prev      = huge(1.0d0)
      stall_count  = 0

      ! Two-stage reconstruction setup: if du_th_plm > du_th (both from input.inp
      ! "du_th [PLM,WENO3]:"), run PLM until du < du_th_plm, then switch
      ! reconstruction to WENO3 and converge at du < du_th. Otherwise the run is
      ! single-stage with the reconstruction set by "Reconstruction scheme:".
      in_plm_stage = (du_th_plm .gt. du_th)
      if (in_plm_stage) then
         rec_method = 'PLM'
         ! Keep the discretization flags consistent with rec_method: Source,
         ! Num_Fluxes, RK_rhs and Apply_BC branch on use_plm/use_weno3, not on
         ! the rec_method string.
         use_plm    = .true.
         use_weno3  = .false.
         write(*,'(A,ES9.2,A,ES9.2)') ' (ATES_main) Two-stage reconstruction: PLM until du <', &
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
      call get_environment_variable('ATES_DIAG_BASE', diag_env)
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
            ! dt_loc = per-cell pseudo-dt ("Time stepping: Local"), or
            ! uniformly the global dt (default; bit-identical updates).
            call eval_dt(W,dt,dt_loc)
            
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
               u1(:,k) = u(:,k) - dt_loc*(dF(:,k) - S(:,k))
            enddo
             
            ! Apply boundary conditions
            call Apply_BC(u1,u1)
                      
            !----------------------------
            
            ! SECOND RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u1,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u1,WL,WR,alpha,dF,S)
                    
            ! Advance in time
            do k = 1,3
               u2(:,k) = (3.0*u(:,k) + u1(:,k) - dt_loc*(dF(:,k) - S(:,k)))/4.0
            enddo

            ! Apply boundary conditions            
            call Apply_BC(u2,u2)
            
            !----------------------------
            
            ! THIRD RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u2,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u2,WL,WR,alpha,dF,S)
                  
            ! Advance in time
            do k = 1,3
               u(:,k) = (u(:,k) + 2.0*(u2(:,k) - dt_loc*(dF(:,k) - S(:,k))))/3.0
            enddo

            ! Apply boundary conditions
            call Apply_BC(u,u)
 		
		!------------------------------------------------!
 		
 		!---- Ionization Equilibrium ----!
 		
            ! Extract primitive variables
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            
            ! Evaluate species densities, ne and n_tot (single policy point)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                       nheiii,nheiTR,nm,ne,n_tot)

            ! Temperature profile
            call comp_T_from_p(p,n_tot,ne,T)

            ! Phase 3a: refresh the lagged H(n=2) Balmer source + heating from
            ! the current state before the ionization/energy solve.
            if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)

            ! Evaluate ionization equilibrium
            call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
            
            ! Evaluate partial densities, ne and n_tot (single policy point)
            call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,  &
                                       nheiii,nheiTR,nm,ne,n_tot)

            ! Evaluate updated pressure
            call comp_p_from_T(T,n_tot,ne,p)
            
            ! Convert to primitive profiles
            W(:,1) = rho
            W(:,2) = v
            W(:,3) = p
            
            ! Revert to conservative 
            call W_to_U(W,u)
  	
  	      !------------------------------------------------!
  	      
  	      ! Evolve solution in time due to source terms
            !     using a semi-implicit step (default) or the original
            !     explicit forward-Euler update ("Energy solver: Explicit")

            if (use_semi_implicit_energy) then
               call solve_energy_semi_implicit(u,W,dt_loc,heat,cool,f_sp)
            else
               u(:,3) = u(:,3) + dt_loc*(heat - cool)
            endif

		call Apply_BC(u,u)

            !------------------------------------------------!
            
            ! Convert to physical variables and extract profiles
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            E   = u(:,3)
            
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
                       W(j,1)*n0, W(j,2)*v0, T(j)*T0, heat(j), cool(j)
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
            
            ! Evaluate variation of time derivative          
      	u(:,2) = u(:,2) + 1.0e-16 ! To avoid division by zero          
	      
	      ! --- Infty-norm
	      dtu = max(maxval(abs(1.0-u(j_min:N,1)/u_old(j_min:N,1))), &
	      	    maxval(abs(1.0-u(j_min:N,2)/u_old(j_min:N,2))))
		dtu = max(maxval(abs(1.0-u(j_min:N,3)/u_old(j_min:N,3))), &
	      	    dtu)

            ! Periodic steady-residual monitor (Resid tol option). R = du/dt
            ! evaluated fresh on the current state (one Reconstruct + RK_rhs;
            ! reuses this step's heat/cool, which near convergence are
            ! consistent with the final u). du/dtu can be small at a premature
            ! operator-split balance while R is large (premature WASP golden:
            ! du,dtu tiny but R_energy ~ 30), so R is the trustworthy gate.
            if ((resid_th .gt. 0.0d0 .or. use_newton_solver) .and.       &
                mod(count, N_resid) .eq. 0) then
               call assemble_residual(u, heat, cool, Rres)
               call residual_norms(Rres, u, resid_c)
               resid_max = maxval(resid_c)
               write(*,'(A,I0,A,3ES11.3,A,ES11.3)') '   [resid] step ',     &
                    count,'  R(m,p,E)=',resid_c,'  max=',resid_max
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
               is_mom_const = (du .lt. du_th)   .and. is_level_stable
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

            ! Production Newton finish ("Solver: Newton"): once the marching
            ! warm-up has brought the steady residual below newton_R_switch
            ! (residual-based, NOT du-based: du dips transiently while the
            ! energy residual is still O(1-30)), hand over to the JFNK
            ! steady solver, refresh the thermodynamic state from the solved
            ! u, and exit the loop (standard final outputs and
            ! post-processing follow).
            if (use_newton_solver .and. .not.in_plm_stage .and.          &
                resid_max .lt. newton_R_switch) then
               if (valve_eps .le. 0.0d0) then
                  valve_eps = 1.0d-4
                  write(*,'(A)') ' (ATES_main) Solver: Newton -> '//     &
                       'adopting smooth base valve, eps = 1.0e-4'
               endif
               resid_max = resid_th
               if (resid_max .le. 0.0d0) resid_max = 1.0d-3
               write(*,'(A,I0,A,ES10.2)') ' (ATES_main) Newton finish '// &
                    'at step ', count, ', target ||R|| <', resid_max
               call solve_steady_jfnk(u, f_sp, resid_max, 500, 1.0d0, 40, j)
               call U_to_W(u,W)
               rho = W(:,1);  v = W(:,2);  p = W(:,3);  E = u(:,3)
               call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii, &
                                          nheiii,nheiTR,nm,ne,n_tot)
               call comp_T_from_p(p,n_tot,ne,T)
               if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)
               call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
               call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii, &
                                          nheiii,nheiTR,nm,ne,n_tot)
               if (j .eq. 0) then
                  is_mom_const = .true.       ! exit the marching loop
                  is_stalled   = .false.
               else
                  ! JFNK failed; it returned its best iterate. Do NOT
                  ! accept an unconverged state as the answer -- fall back
                  ! to plain time-marching (du-based stops re-armed).
                  use_newton_solver = .false.
                  write(*,'(A,I0,A)') ' (ATES_main) JFNK failed (info=', &
                       j, '); resuming time-marching with du-based stops'
               endif
            endif

            ! Write to standard output (4th column: relative mass-flux level
            ! change over the last N_stall steps; -1 until the window fills)
            write(*,*) count,du,dtu,lev_rel
				
            !---------------------------------------------------!
            
            ! Detect NaNs
            do j = 1-Ng,N+Ng
            	do k = 1,3
            		dum = u(j,k)
            		if (dum.ne.dum) then
                              write(*,*)
            			write(*,'(A20,F8.6)') 'NaN detected at r = ',r(j)
                              write(*,'(A6,E13.6)') 'rho = ',W(j,1)*n0
                              write(*,'(A4,E13.6)') 'v = ',W(j,2)*v0/1.0e5
                              write(*,'(A4,E13.6)') 'p = ',W(j,3)*p0
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
            
      !---------------------------------------------------!
            
      ! End of temporal while loop
      enddo
      write(*,*) '(ATES_main.f90) Time integration done.'

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

      !---------------------------------------------------!

      ! Phase-3a: refresh the H(n=2) diagnostics from the converged state for
      ! write_excited_H below (the in-loop Balmer source was lagged one step).
      ! The Balmer photoelectric heating is ~1e-4 of the total, so we do NOT
      ! re-solve ioniz_eq here: that would desync T/f_sp from the converged
      ! hydro state for a negligible heat correction.
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,exc_rel)

      !---------------------------------------------------!

      ! Write final thermodynamic and ionization profiles
      write(*,*) '(ATES_main.f90) Writing final results to file..'
      call write_output(rho,v,p,T,heat,cool,eta,                         &
                        nhi,nhii,nhei,nheii,nheiii,nheiTR,nm,'eq')

      ! Phase-2 diagnostic: per-channel radiative cooling vs radius
      ! (reuses eval_cool's coefficients; see Huang et al. 2023 Fig. 10).
      call write_cool_breakdown_eq(T,rho,f_sp)

      ! Phase-3a diagnostic: H(n=2) populations, Balmer proton source, and
      ! photoelectric/de-excitation heating vs radius (Huang Figs. 11/27/10/26).
      if (use_excited_H) call write_excited_H

      !---------------------------------------------------!

      ! Post processing to include advection
      write(*,*) '(ATES_main.f90) Starting the post processing routine..'

      call post_process_adv(rho,v,p,T,heat,cool,eta,    &
                            nhi,nhii,nhei,nheii,nheiii,nheiTR,nm)
      
      write(*,*) '(ATES_main.f90) Post processing routine done.'                           
      
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
      open(unit = outfile, file = 'ATES.out', access = 'append' )
      	write(outfile,101) ' '
      	write(outfile,102) ' - Log10 of steady-state Mdot = ', Mdot, ' g/s'
      close(unit = outfile)
      
101   format (A35,F5.2,A4)
102   format (A32,F5.2,A4)

      !---------------------------------------------------!
      
      ! End of program
      end program Hydro_ioniz

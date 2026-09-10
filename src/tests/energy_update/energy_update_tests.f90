      program energy_update_tests
      ! Tests of the residual-controlled temperature solve of the implicit
      ! energy source step (docs/PLAN_20260906_rev2.md step B2).
      !
      ! The solve is driven through the three routines of
      ! energy_semi_implicit that own the iteration -- energy_balance_init,
      ! energy_balance_update, energy_balance_finish -- so that the cooling
      ! the solver sees is an analytic law with a closed form or an
      ! independently integrable root, instead of eval_cool. The production
      ! wrapper solve_energy_semi_implicit differs from what runs here only
      ! in where the cooling comes from and in what it does with the answer.
      !
      ! Each assertion prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line. Lines beginning with two spaces or with DIAGNOSTIC are context.
      ! Exit status is nonzero if any assertion fails.
      use global_parameters
      use energy_semi_implicit
      use energy_update_cooling_models
      implicit none

      integer :: n_fail

      n_fail = 0
      ! One physical cell plus its ghosts. Every entry is given the same
      ! problem and cell 1 is the one asserted on.
      N = 1

      call test_linear_cooling(n_fail)
      call test_falling_branch(n_fail)
      call test_frozen_heating(n_fail)
      call test_no_bracket(n_fail)
      call test_temperature_floor(n_fail)
      call test_failed_update_assembles_nothing(n_fail)

      write(*,*)
      if (n_fail .gt. 0) then
         write(*,'(a,i0,a)') 'energy_update: ', n_fail, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'energy_update: every assertion passed'

      contains

      ! ------------------------------------------------------!

      subroutine check(name, measured, reference, tol, nf)
      ! Relative comparison, absolute when the reference is zero.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      real*8 :: d
      if (abs(reference) .gt. 0.0d0) then
         d = abs(measured - reference)/abs(reference)
      else
         d = abs(measured - reference)
      endif
      if (d .le. tol) then
         write(*,'(a,a,a,es14.7,a,es14.7,a,es9.2)') 'PASS ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(a,a,a,es14.7,a,es14.7,a,es9.2)') 'FAIL ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nf = nf + 1
      endif
      end subroutine check

      subroutine check_int(name, measured, reference, nf)
      character(len=*), intent(in) :: name
      integer, intent(in) :: measured, reference
      integer, intent(inout) :: nf
      if (measured .eq. reference) then
         write(*,'(a,a,a,i0,a,i0,a)') 'PASS ', name, ' measured=',        &
            measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(a,a,a,i0,a,i0,a)') 'FAIL ', name, ' measured=',        &
            measured, ' reference=', reference, ' tol=0'
         nf = nf + 1
      endif
      end subroutine check_int

      subroutine check_above(name, measured, threshold, nf)
      ! The measured quantity must EXCEED the threshold. Used for the
      ! residual the previous two-iteration update leaves behind: the test
      ! is red for that algorithm by construction and would stop being a
      ! meaningful reference if it ever fell below the acceptance tolerance.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, threshold
      integer, intent(inout) :: nf
      if (measured .gt. threshold) then
         write(*,'(a,a,a,es14.7,a,es14.7,a)') 'PASS ', name,              &
            ' measured=', measured, ' reference=>', threshold, ' tol=0'
      else
         write(*,'(a,a,a,es14.7,a,es14.7,a)') 'FAIL ', name,              &
            ' measured=', measured, ' reference=>', threshold, ' tol=0'
         nf = nf + 1
      endif
      end subroutine check_above

      ! ------------------------------------------------------!

      subroutine run_solve(model, T_old, a, heat, cool0, T_floor, T_ceil, s)
      ! Drive the bracketed solve with one of the analytic cooling laws.
      ! model: 1 linear, 2 falling, 3 constant.
      integer, intent(in) :: model
      real*8, dimension(1-Ng:N+Ng), intent(in) :: T_old, a, heat, cool0
      real*8, dimension(1-Ng:N+Ng), intent(in) :: T_floor, T_ceil
      type(energy_balance_state), intent(out) :: s
      real*8, dimension(1-Ng:N+Ng) :: C
      call energy_balance_init(s, T_old, a, heat, cool0, T_floor, T_ceil)
      do while (s%n_active .gt. 0 .and. s%iters .lt. energy_max_iter)
         select case (model)
         case (1)
            call cool_linear(s%T_eval, C)
         case (2)
            call cool_falling(s%T_eval, C)
         case default
            call cool_constant(s%T_eval, C)
         end select
         call energy_balance_update(s, C)
      enddo
      call energy_balance_finish(s)
      end subroutine run_solve

      ! ------------------------------------------------------!

      real*8 function residual_of(T, C, T_old, a, heat)
      ! R = u(T) - u(T_old) - a (heat - cool), monatomic branch, the same
      ! quantity the solver tests.
      real*8, intent(in) :: T, C, T_old, a, heat
      residual_of = T/(gamma_ad - 1.0d0) - T_old/(gamma_ad - 1.0d0)        &
                    - a*(heat - C)
      end function residual_of

      real*8 function scale_of(T_old, a, heat, C)
      real*8, intent(in) :: T_old, a, heat, C
      scale_of = abs(T_old/(gamma_ad - 1.0d0)) + a*(abs(heat) + abs(C))
      end function scale_of

      ! ------------------------------------------------------!

      subroutine test_linear_cooling(nf)
      ! (a) cool = c (T - T_eq), heat = 0. The implicit step
      !     u(T) - u(T_old) + a c (T - T_eq) = 0
      ! is linear in T with the closed-form root below.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor, T_ceil
      type(energy_balance_state) :: s
      real*8 :: cv, T_exact, C_ret

      write(*,*)
      write(*,'(a)') '---- linear_cooling ----'
      lin_c   = 3.0d0
      lin_Teq = 0.4d0
      T_old  = 1.0d0
      a      = 0.25d0
      heat   = 0.0d0
      T_floor = 0.01d0
      T_ceil  = 1.0d3
      cool0  = lin_c*(T_old - lin_Teq)

      cv = 1.0d0/(gamma_ad - 1.0d0)
      T_exact = (cv*T_old(1) + a(1)*lin_c*lin_Teq)/(cv + a(1)*lin_c)

      call run_solve(1, T_old, a, heat, cool0, T_floor, T_ceil, s)

      call check_int('linear_cooling_status', s%status(1),                &
                     ENERGY_UPDATE_OK, nf)
      call check('linear_cooling_closed_form', s%T_new(1), T_exact,       &
                 1.0d-12, nf)
      C_ret = lin_c*(s%T_new(1) - lin_Teq)
      call check('linear_cooling_returned_cool_at_returned_T',            &
                 s%cool_new(1), C_ret, 1.0d-14, nf)
      call check('linear_cooling_residual', s%res_scaled(1), 0.0d0,       &
                 energy_res_tol, nf)
      write(*,'(a,i0)') '  cooling evaluations: ', s%iters
      end subroutine test_linear_cooling

      ! ------------------------------------------------------!

      subroutine test_falling_branch(nf)
      ! (b) and (f). cool = Cp exp(-(T - Tp)/w) DECREASES with temperature.
      ! The reference root is an independent bisection of the same equation
      ! to the precision of the arithmetic; the previous two-iteration
      ! update is run on the same problem and its residual measured.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor, T_ceil
      type(energy_balance_state) :: s
      real*8 :: T_ref, rel_err, T_old2i, C_2i, R_2i, sc_2i, R_new, sc_new
      real*8 :: lo, hi, mid, f_mid, C_ret

      write(*,*)
      write(*,'(a)') '---- falling_branch ----'
      fall_Cp = 1.0d0
      fall_Tp = 1.0d0
      fall_w  = 0.2d0
      T_old   = 1.0d0
      a       = 0.2d0
      heat    = 2.0d0
      T_floor = 0.01d0
      T_ceil  = 1.0d2
      cool0   = cool_falling_scalar(T_old(1))

      ! Independent reference: bisection on [T_old, 2], the interval the
      ! solver's upward search covers, to 200 halvings.
      lo = T_old(1)
      hi = 2.0d0
      do while (hi - lo .gt. 1.0d-15*hi)
         mid = 0.5d0*(lo + hi)
         f_mid = residual_of(mid, cool_falling_scalar(mid), T_old(1),     &
                             a(1), heat(1))
         if (f_mid .lt. 0.0d0) then
            lo = mid
         else
            hi = mid
         endif
      enddo
      T_ref = 0.5d0*(lo + hi)

      call run_solve(2, T_old, a, heat, cool0, T_floor, T_ceil, s)

      call check_int('falling_branch_status', s%status(1),                &
                     ENERGY_UPDATE_OK, nf)
      call check('falling_branch_root', s%T_new(1), T_ref, 1.0d-11, nf)
      C_ret = cool_falling_scalar(s%T_new(1))
      call check('falling_branch_cool_at_returned_T', s%cool_new(1),      &
                 C_ret, 1.0d-14, nf)
      R_new  = residual_of(s%T_new(1), s%cool_new(1), T_old(1), a(1),     &
                           heat(1))
      sc_new = scale_of(T_old(1), a(1), heat(1), s%cool_new(1))
      call check('falling_branch_residual', abs(R_new)/sc_new, 0.0d0,     &
                 energy_res_tol, nf)
      write(*,'(a,i0)') '  cooling evaluations: ', s%iters

      ! (f) The previous update on the same problem: two iterations with a
      ! frozen cooling derivative entering as 1 + c |dC/dT|, the cooling
      ! refreshed once, and no residual test.
      call two_iteration_update(T_old(1), a(1), heat(1), cool0(1),        &
                                T_old2i, C_2i)
      R_2i  = residual_of(T_old2i, C_2i, T_old(1), a(1), heat(1))
      sc_2i = scale_of(T_old(1), a(1), heat(1), C_2i)
      write(*,'(a,es14.7,a,es14.7)') '  two-iteration T = ', T_old2i,     &
                                     '  converged T = ', T_ref
      call check_above('falling_branch_two_iteration_residual_red',       &
                       abs(R_2i)/sc_2i, energy_res_tol, nf)
      rel_err = abs(T_old2i - T_ref)/T_ref
      write(*,'(a,es14.7)')                                               &
         '  DIAGNOSTIC two-iteration relative temperature error: ',       &
         rel_err
      end subroutine test_falling_branch

      ! ------------------------------------------------------!

      subroutine two_iteration_update(T_old, a, heat, cool0, T_ret, C_ret)
      ! The update this step replaces, reproduced for the red measurement of
      ! test (f): exactly two Newton-like iterations on the temperature form
      ! with the cooling derivative frozen at T_old and entering through its
      ! absolute value, the cooling refreshed after the first iteration only,
      ! a clamp at 0.01, and no test on what comes out. The returned cooling
      ! belongs to the previous iterate.
      real*8, intent(in)  :: T_old, a, heat, cool0
      real*8, intent(out) :: T_ret, C_ret
      real*8 :: T_trial, cool_trial, delta_T, dC_dT, c_factor, F, dF_dT
      integer :: iter

      T_trial    = T_old
      cool_trial = cool0
      delta_T    = max(1.0d-5, 1.0d-5*T_trial)
      dC_dT      = (cool_falling_scalar(T_trial + delta_T) - cool_trial)  &
                   /delta_T
      c_factor   = a*(gamma_ad - 1.0d0)
      do iter = 1, 2
         F     = T_trial - T_old - c_factor*(heat - cool_trial)
         dF_dT = 1.0d0 + c_factor*abs(dC_dT)
         T_trial = T_trial - F/dF_dT
         if (T_trial .lt. 0.01d0) T_trial = 0.01d0
         if (iter .eq. 1) cool_trial = cool_falling_scalar(T_trial)
      enddo
      T_ret = T_trial
      C_ret = cool_trial
      end subroutine two_iteration_update

      ! ------------------------------------------------------!

      subroutine test_frozen_heating(nf)
      ! (c) The heating is frozen at the state the step starts from, by
      ! construction: the photoionization rates do not follow T within one
      ! step. The test states what that means for the returned residual --
      ! the balance the solve closes is the one with the frozen heating, and
      ! it closes to the tolerance -- and measures how far the same state is
      ! from the balance with a strongly temperature-dependent heating
      ! evaluated at the returned temperature.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor, T_ceil
      type(energy_balance_state) :: s
      real*8 :: H0, w_h, heat_at_T_new, R_frozen, R_unfrozen, sc

      write(*,*)
      write(*,'(a)') '---- frozen_heating ----'
      ! A heating that halves over a temperature change of w_h.
      H0  = 4.0d0
      w_h = 0.3d0
      lin_c   = 2.0d0
      lin_Teq = 0.5d0
      T_old   = 1.0d0
      a       = 0.3d0
      T_floor = 0.01d0
      T_ceil  = 1.0d3
      heat    = H0                       ! frozen at T_old
      cool0   = lin_c*(T_old - lin_Teq)

      call run_solve(1, T_old, a, heat, cool0, T_floor, T_ceil, s)

      call check_int('frozen_heating_status', s%status(1),                &
                     ENERGY_UPDATE_OK, nf)
      R_frozen = residual_of(s%T_new(1), s%cool_new(1), T_old(1), a(1),   &
                             heat(1))
      sc = scale_of(T_old(1), a(1), heat(1), s%cool_new(1))
      call check('frozen_heating_residual_with_frozen_heat',              &
                 abs(R_frozen)/sc, 0.0d0, energy_res_tol, nf)

      heat_at_T_new = H0*exp(-(s%T_new(1) - T_old(1))/w_h)
      R_unfrozen = residual_of(s%T_new(1), s%cool_new(1), T_old(1), a(1), &
                               heat_at_T_new)
      write(*,'(a,es14.7)')                                               &
         '  DIAGNOSTIC |R|/scale if the heating followed T: ',            &
         abs(R_unfrozen)/sc
      write(*,'(a)') '  The solve closes the balance it is given. A' //   &
                     ' heating that follows T'
      write(*,'(a)') '  belongs to the coupled source step of B3c, not' //&
                     ' to this solve.'
      end subroutine test_frozen_heating

      ! ------------------------------------------------------!

      subroutine test_no_bracket(nf)
      ! (d) The heating is so large that the balance is not reached anywhere
      ! below the ceiling of the physical range: R(T_ceil) < 0.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor, T_ceil
      type(energy_balance_state) :: s
      real*8 :: R_ceiling

      write(*,*)
      write(*,'(a)') '---- no_bracket ----'
      const_C0 = 1.0d0
      T_old    = 1.0d0
      a        = 1.0d0
      heat     = 1.0d6
      cool0    = const_C0
      T_floor  = 0.01d0
      T_ceil   = 1.0d1
      R_ceiling = residual_of(T_ceil(1), const_C0, T_old(1), a(1), heat(1))
      write(*,'(a,es14.7)') '  R at the ceiling: ', R_ceiling

      call run_solve(3, T_old, a, heat, cool0, T_floor, T_ceil, s)

      call check_int('no_bracket_status', s%status(1),                    &
                     ENERGY_UPDATE_NO_BRACKET, nf)
      call check_int('no_bracket_reason', s%reason(1),                    &
                     ENERGY_REASON_CEILING, nf)
      call check('no_bracket_last_iterate_at_ceiling', s%T_new(1),        &
                 T_ceil(1), 1.0d-14, nf)
      end subroutine test_no_bracket

      ! ------------------------------------------------------!

      subroutine test_temperature_floor(nf)
      ! (e) A constant cooling large enough that the balance asks for a
      ! temperature below the floor. The update must fail, must not return a
      ! state sitting on the floor, and must report the energy that would
      ! have to be supplied for the floor to be a solution.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor, T_ceil
      type(energy_balance_state) :: s
      real*8 :: R_floor_exact, missing_source_exact
      integer :: hits_before

      write(*,*)
      write(*,'(a)') '---- temperature_floor ----'
      const_C0 = 5.0d0
      T_old    = 1.0d0
      a        = 1.0d0
      heat     = 0.0d0
      cool0    = const_C0
      T_floor  = 0.01d0
      T_ceil   = 1.0d3

      R_floor_exact = residual_of(T_floor(1), const_C0, T_old(1), a(1),   &
                                  heat(1))
      missing_source_exact = R_floor_exact/a(1)

      hits_before = n_energy_floor_hits
      call run_solve(3, T_old, a, heat, cool0, T_floor, T_ceil, s)

      call check_int('floor_status', s%status(1), ENERGY_UPDATE_FLOOR, nf)
      call check_int('floor_reason', s%reason(1),                         &
                     ENERGY_REASON_FLOOR_REACHED, nf)
      call check('floor_last_iterate_is_the_floor', s%T_new(1),           &
                 T_floor(1), 1.0d-14, nf)
      call check('floor_residual_energy', s%res_energy(1), R_floor_exact, &
                 1.0d-12, nf)
      call check('floor_missing_source', s%res_energy(1)/a(1),            &
                 missing_source_exact, 1.0d-12, nf)
      write(*,'(a)') '  The floor is reported, not adopted: the solve' //  &
                     ' returns a failure and'
      write(*,'(a)') '  the caller assembles no state from it.'
      ! The counters of the module are maintained by the production wrapper,
      ! not by the core; this test asserts only that the core did not touch
      ! them behind the wrapper's back.
      call check_int('floor_counters_owned_by_the_wrapper',               &
                     n_energy_floor_hits, hits_before, nf)
      end subroutine test_temperature_floor

      ! ------------------------------------------------------!

      subroutine test_failed_update_assembles_nothing(nf)
      ! (f) THE PRODUCTION WRAPPER ON A FAILURE. A cell with no particles at
      ! all has no source coefficient dt/(n_tot + n_e), so the admissibility
      ! test of energy_balance_init refuses it before any cooling is asked
      ! for. What is asserted is what the caller depends on: the wrapper
      ! returns a status other than ENERGY_UPDATE_OK, it does not stop the
      ! program, and it writes NOTHING -- the pressure row of W, the energy
      ! row of u and the cooling column come back exactly as they were
      ! handed in, so the attempted step the caller then refuses restores a
      ! state this routine did not touch.
      integer, intent(inout) :: nf
      real*8, dimension(3,1-Ng:N+Ng) :: u, W, u_in, W_in
      real*8, dimension(1-Ng:N+Ng)   :: dt_a, heat, cool, cool_in, T_start
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_sp
      integer :: st
      real*8  :: dmax

      write(*,*)
      write(*,'(a)') '---- failed_update_assembles_nothing ----'
      n0 = 1.0d10
      T0 = 1.0d4
      p0 = n0*kb_erg*T0
      ! No mass and no species: n_tot + n_e = 0 in every cell.
      W = 0.0d0
      W(3,:) = 1.0d0
      u = 0.0d0
      u(3,:) = 1.5d0
      f_sp = 0.0d0
      dt_a = 1.0d0
      heat = 0.0d0
      cool = 7.0d0
      T_start = 1.0d0
      u_in = u;  W_in = W;  cool_in = cool

      call solve_energy_semi_implicit(u, W, dt_a, heat, cool, f_sp, 0,    &
                                      status = st, T_start = T_start)

      call check_int('failed_update_status_is_not_ok',                    &
                     merge(1, 0, st .ne. ENERGY_UPDATE_OK), 1, nf)
      call check_int('failed_update_reports_the_status_through_the_module',&
                     energy_update_last_status, st, nf)
      dmax = maxval(abs(W(3,:) - W_in(3,:)))
      call check('failed_update_left_the_pressure_row', dmax, 0.0d0,      &
                 0.0d0, nf)
      dmax = maxval(abs(u(3,:) - u_in(3,:)))
      call check('failed_update_left_the_energy_row', dmax, 0.0d0,        &
                 0.0d0, nf)
      dmax = maxval(abs(cool - cool_in))
      call check('failed_update_left_the_cooling', dmax, 0.0d0, 0.0d0, nf)
      end subroutine test_failed_update_assembles_nothing

      end program energy_update_tests

      program base_boundary_continuity_probe
      ! IS THE BASE FACE DENSITY A CONTINUOUS FUNCTION OF THE STATE, AND
      ! WHERE IS IT NOT (item L26 of docs/PLAN_20260916_rev3.md section 6).
      !
      ! The base face state is
      !
      !     rho_b = (1 - w_rev) rho_res + w_rev rho_rev ,
      !     w_rev = s_wind w_wind(M_wind) + (1 - s_wind) w_i ,
      !     s_wind = A(M_wind) C(d_window) ,
      !     d_window = sqrt(<F^2> - <F>^2)/|<F>| over r >= r_flux ,
      !
      ! with w_i the local branch weight of the first interior cell and
      ! F = rho v r^2.  A is the amplitude weight, 0 at M_wind = 0 and 1 at
      ! |M_wind| >= base_wind_window_amplitude; C is the shape weight, 1 at
      ! d_window <= base_wind_window_spread and 0 at twice it.
      !
      ! THE ZERO WINDOW IS THE POINT THIS PROGRAM IS ABOUT.  C is scale free,
      ! so on its own it has no limit there: along a uniform flux eps F it is
      ! 1 at every eps and along a flux whose spread is above the gate it is
      ! 0 at every eps, and the two paths reach the same state.  A vanishes
      ! quadratically with the window flux, so s_wind and its derivative
      ! vanish along every direction and the two limits coincide.  The rows
      ! below MEASURE that, by calling the production boundary; this program
      ! re-implements no part of the condition.
      !
      ! WHAT IT REPORTS, at every sampled point of every path below:
      ! rho_b, the two ghost densities the same call returns, w_i, s_wind,
      ! w_rev, du_window, have_F, and rho_res and rho_rev separately, since a
      ! change of branch moves rho_b only in so far as the two isentropes
      ! differ.  Each path also gets one-sided directional derivatives of
      ! rho_b by finite differences at three step sizes, so that a jump (the
      ! estimates grow like 1/h), a kink (the two sides converge to different
      ! finite numbers) and a smooth point (the estimates agree) are told
      ! apart by measurement and not by inspection.
      !
      ! THE PATHS, in the order of section 6 of the plan:
      !   (i)   a uniform window flux approaching zero from both signs;
      !   (ii)  nonuniform window shapes approaching zero, three with a
      !         relative spread above twice the threshold and one below it;
      !   (iii) zero window mean with nonzero variance;
      !   (iv)  the first cell and the window disagreeing in sign;
      !   (v)   the cell that carries the window extremum changing;
      !   (vi)  a base flux the window does not carry;
      !   (vii) the fully static state, and a static base under a live wind.
      ! Plus the repeated-evaluation and call-order identities of R48: the
      ! boundary must be a function of its argument, so ten evaluations of one
      ! state, and an evaluation with a full flux evaluation in between, must
      ! return the same bits.  Those two are assertions; everything else is a
      ! measurement printed for the memo.
      !
      ! The state is whatever the run directory this program is started in
      ! holds: it calls input_read and init exactly as EXHALE_main does, so
      ! the grid, the reservoir, the composition of cell 1 and every boundary
      ! option are the ones that case runs with.
      use global_parameters
      use Read_input,     only: input_read
      use Initialization, only: init
      use lya_rt,         only: lya_rt_allocate_arrays
      use excited_hydrogen, only: excited_H_allocate_arrays
      use ionization_equilibrium, only: ioniz_eq_allocate_arrays
      use mol_rates,      only: h2_thermochemistry_init
      use Conversion,     only: W_to_U
      use Reconstruction_step, only: Reconstruct
      use RK_integration, only: RK_rhs, face_flux
      use base_boundary,  only: base_boundary_states,                      &
                                read_base_branch_options,                  &
                                wind_window_mass_flux,                     &
                                characteristic_branch_weight,              &
                                base_face_mach_blend,                      &
                                base_wind_window_spread,                   &
                                base_face_blend_last,                      &
                                base_face_Mi_last,                         &
                                base_face_Mwind_last,                      &
                                base_face_swind_last,                      &
                                base_face_d_window_last,                  &
                                base_face_rho_res_last,                    &
                                base_face_rho_rev_last
      implicit none

      real*8, allocatable :: W0(:,:), Wq(:,:), u(:,:), f_sp(:,:)
      real*8, allocatable :: WLr(:,:), WRr(:,:), dFr(:,:), Sr(:,:)
      real*8  :: Wface(3), Wface_lower(3)
      real*8, allocatable :: Wghost(:,:)
      real*8  :: F0, du0
      logical :: have0
      integer :: j, j0, nw, n_fail
      character(len=64) :: label

      ! The two cells whose window flux is bumped in path (v), and the bump
      ! that puts du_window inside the handover band (threshold to twice it).
      integer :: jA, jB
      ! The fixed bump of path (v), set so that the window's own relative
      ! spread sits INSIDE the handover band, where the weight is strictly
      ! between 0 and 1 and a change of the extremal cell can move rho_b at
      ! all.  Two bumped cells out of nw contribute
      ! d_window = alpha sqrt(2/nw) to the relative standard deviation, so
      ! alpha is set from the threshold and the window's size and not fixed
      ! in advance: a bump chosen for one window would saturate another.
      real*8 :: alpha_extremum

      ! The three difference quotients the last call to `derivatives` made,
      ! and the value it took them at, so that the verdict rows can be taken
      ! on the same numbers the measurement rows print.
      real*8 :: d_last(3), rho_last, d_above(3), h_last(3)

      n_fail = 0
      label = 'state'
      if (command_argument_count() .ge. 1) call get_command_argument(1,label)

      call input_read
      open(unit = outfile, file = 'EXHALE_setup.out')
      call lya_rt_allocate_arrays
      call excited_H_allocate_arrays
      call ioniz_eq_allocate_arrays
      call h2_thermochemistry_init

      allocate(W0(3,1-Ng:N+Ng), Wq(3,1-Ng:N+Ng), u(3,1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species))
      allocate(Wghost(3,1-Ng:0))
      allocate(WLr(3,1-Ng:N+Ng), WRr(3,1-Ng:N+Ng),                        &
               dFr(3,1-Ng:N+Ng), Sr(3,1-Ng:N+Ng))
      call init(W0,u,f_sp)

      j0 = max(j_flux,1)
      nw = N - j0 + 1
      jA = min(N, j0 + 10)
      jB = min(N, j0 + 20)

      ! The options before the thresholds are read off them: a control
      ! experiment that overrides the threshold has to move this path's bump
      ! with it, or the sweep leaves the handover band it exists to probe.
      call read_base_branch_options()
      alpha_extremum = 1.5d0*base_wind_window_spread                      &
                       *sqrt(0.5d0*dble(nw))

      call wind_window_mass_flux(W0, F0, du0, have0)

      write(*,'(A)') ''
      write(*,'(A)') '================ base boundary continuity probe ====='
      write(*,'(A,A)')    ' state label            : ', trim(label)
      write(*,'(A,I6,A,I6,A,ES16.9)') ' window: j_flux =', j0,            &
           '  cells =', nw, '  r(j_flux) =', r(j0)
      write(*,'(A,ES16.9,A,ES16.9)') ' window flux of the state F0 =',    &
           F0, '  its spread du =', du0
      write(*,'(A,L1)') ' have_F of the state    : ', have0
      write(*,'(A,ES12.5,A,ES12.5)') ' base_face_mach_blend =',           &
           base_face_mach_blend, '  base_wind_window_spread =',           &
           base_wind_window_spread
      write(*,'(A,ES16.9,A,ES16.9)') ' cell 1: rho =', W0(1,1),           &
           '  v =', W0(2,1)*v0
      write(*,'(A)') ''
      write(*,'(A)') ' Columns of every ROW line:'
      write(*,'(A)') '   path  t  rho_b  rho_gh0  rho_gh1  M_i  M_wind'// &
                     '  w_i  s_wind  w_rev  du_window  have_F'//        &
                     '  rho_res  rho_rev'
      write(*,'(A)') '   rho_b, the ghost densities, rho_res and rho_rev'//&
                     ' are in code units (n0).'
      write(*,'(A)') ''

      ! ---------------------------------------------------------------- !
      ! the unperturbed state
      ! ---------------------------------------------------------------- !
      call sample('state ', 0, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! (i) a uniform window flux approaching zero from both signs
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- (i) uniform window flux eps*F0, eps -> 0 ----'
      call decade_sweep('unif+', 1, +1.0d0)
      call decade_sweep('unif-', 1, -1.0d0)
      call sample('unif0 ', 1, 0.0d0)
      call derivatives('unif+', 1, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_uniform_plus_derivative')
      call derivatives('unif-', 1, 0.0d0, -1.0d0)
      call derivative_converges('base_continuity_uniform_minus_derivative')
      call jump_report('uniform window vs the zero window', 1)
      call zero_window_limit('base_continuity_uniform_plus_limit', 1,      &
                             +1.0d0)
      call zero_window_limit('base_continuity_uniform_minus_limit', 1,     &
                             -1.0d0)

      ! ---------------------------------------------------------------- !
      ! (ii) nonuniform window shapes approaching zero
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- (ii) nonuniform window shapes, eps -> 0 ----'
      write(*,'(A)') '   shape 21: 1 + 0.5 sin, spread ~ 1.0'
      call decade_sweep('sin  ', 21, +1.0d0)
      call sample('sin0  ', 21, 0.0d0)
      call derivatives('sin  ', 21, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_sin_derivative')
      call zero_window_limit('base_continuity_sin_limit', 21, +1.0d0)
      write(*,'(A)') '   shape 22: alternating 1, 2; spread ~ 0.67'
      call decade_sweep('alt12', 22, +1.0d0)
      call sample('alt120', 22, 0.0d0)
      call derivatives('alt12', 22, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_alt12_derivative')
      call zero_window_limit('base_continuity_alt12_limit', 22, +1.0d0)
      write(*,'(A)') '   shape 23: linear ramp 0.5 .. 1.5, spread ~ 1.0'
      call decade_sweep('ramp ', 23, +1.0d0)
      call sample('ramp0 ', 23, 0.0d0)
      call derivatives('ramp ', 23, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_ramp_derivative')
      call zero_window_limit('base_continuity_ramp_limit', 23, +1.0d0)
      write(*,'(A)') '   shape 24: 1 + 0.002 sin, spread ~ 0.004 '//      &
                     '(below the threshold)'
      call decade_sweep('sinsm', 24, +1.0d0)
      call sample('sinsm0', 24, 0.0d0)
      call derivatives('sinsm', 24, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_small_sin_derivative')
      call zero_window_limit('base_continuity_small_sin_limit', 24, +1.0d0)

      ! ---------------------------------------------------------------- !
      ! (iii) zero window mean with nonzero variance
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- (iii) zero mean, nonzero variance ----'
      call decade_sweep('zmean', 3, +1.0d0)
      call sample('zmean0', 3, 0.0d0)
      call derivatives('zmean', 3, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_zero_mean_derivative')
      call zero_window_limit('base_continuity_zero_mean_limit', 3, +1.0d0)

      ! ---------------------------------------------------------------- !
      ! (iv) the first cell and the window disagreeing in sign
      ! ---------------------------------------------------------------- !
      ! The window is held at a uniform +F0 (gas enters) or -F0 (gas leaves)
      ! and the first cell's velocity is swept through zero, so the local and
      ! the window discriminant disagree over half of each sweep.
      write(*,'(A)') ''
      write(*,'(A)') '---- (iv) local and window signs disagreeing ----'
      write(*,'(A)') '   window at +F0, t = multiplier on v(1)'
      call sign_sweep('sgn+ ', 41)
      call derivatives('sgn+ ', 41,  0.0d0, +1.0d0)
      call derivative_converges('base_continuity_sign_sweep_plus_derivative')
      call derivatives('sgn+ ', 41,  0.0d0, -1.0d0)
      call derivative_converges(                                          &
           'base_continuity_sign_sweep_plus_derivative_from_below')
      write(*,'(A)') '   window at -F0, t = multiplier on v(1)'
      call sign_sweep('sgn- ', 42)
      call derivatives('sgn- ', 42,  0.0d0, +1.0d0)
      call derivative_converges('base_continuity_sign_sweep_minus_derivative')
      call derivatives('sgn- ', 42,  0.0d0, -1.0d0)
      call derivative_converges(                                          &
           'base_continuity_sign_sweep_minus_derivative_from_below')

      ! ---------------------------------------------------------------- !
      ! (v) the cell carrying the window extremum changes
      ! ---------------------------------------------------------------- !
      ! Two cells of an otherwise uniform window are bumped, one by a fixed
      ! alpha and one by t.  The window's max is carried by the first while
      ! t < alpha and by the second while t > alpha, and alpha is chosen so
      ! that du_window sits inside the handover band, where s_wind is
      ! strictly between 0 and 1 and the switch can move rho_b.
      write(*,'(A)') ''
      write(*,'(A)') '---- (v) the extremal cell of the window changes ----'
      write(*,'(A,I5,A,I5,A,ES12.5)') '   bumped cells ', jA, ' and ',    &
           jB, ',  fixed bump alpha =', alpha_extremum
      call extremum_sweep('extr ', 5)
      call derivatives('extr ', 5, alpha_extremum, +1.0d0)
      call derivative_converges('base_continuity_extremum_derivative_above')
      d_above = d_last
      call derivatives('extr ', 5, alpha_extremum, -1.0d0)
      call derivative_converges('base_continuity_extremum_derivative_below')
      ! THE SLOPE DOES NOT JUMP WHERE THE EXTREMAL CELL CHANGES.  With the
      ! range functional the two sides converged to -242.2 and +0.0128 on
      ! this state, a factor 1.9e4; a functional with no argmax has one
      ! slope.
      call verdict_rel('base_continuity_extremum_no_slope_jump',           &
                       abs(d_above(3) - d_last(3))                         &
                       /max(abs(d_above(3)), abs(d_last(3)), 1.0d-300),    &
                       0.0d0, 1.0d-3)

      ! ---------------------------------------------------------------- !
      ! (vi) a base flux the window does not carry
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- (vi) base cells 1..6 scaled, window unchanged --'
      call sample('base05', 6, 0.5d0)
      call sample('base10', 6, 1.0d0)
      call sample('base20', 6, 2.0d0)
      call sample('base00', 6, 0.0d0)
      call sample('baseM1', 6, -1.0d0)
      call derivatives('base ', 6, 1.0d0, +1.0d0)
      call derivative_converges('base_continuity_base_scaling_derivative')
      call derivatives('base ', 6, 1.0d0, -1.0d0)
      call derivative_converges(                                          &
           'base_continuity_base_scaling_derivative_from_below')

      ! ---------------------------------------------------------------- !
      ! (vii) the fully static state, and a static base under a live wind
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- (vii) static states ----'
      write(*,'(A)') '   v = 0 in every cell'
      call sample('stat  ', 7, 0.0d0)
      write(*,'(A)') '   v = 0 in every cell, then a uniform window '//   &
                     'flux t*F0'
      call decade_sweep('statW', 71, +1.0d0)
      call derivatives('statW', 71, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_static_window_derivative')
      call derivatives('statW', 71, 0.0d0, -1.0d0)
      call derivative_converges(                                          &
           'base_continuity_static_window_derivative_from_below')
      call zero_window_limit('base_continuity_static_window_limit', 71,    &
                             +1.0d0)
      write(*,'(A)') '   v = 0 in every cell, then a nonuniform window '// &
                     'flux t*F0 (spread above the gate)'
      call decade_sweep('statS', 73, +1.0d0)
      call sample('statS0', 73, 0.0d0)
      call derivatives('statS', 73, 0.0d0, +1.0d0)
      call derivative_converges('base_continuity_static_closed_gate_derivative')
      call zero_window_limit('base_continuity_static_closed_gate_limit',   &
                             73, +1.0d0)
      write(*,'(A)') '   v = 0 in cells 1..6, the window as the state '// &
                     'carries it'
      call sample('statB ', 72, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! the local readings of the base flux, against the face flux
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- local readings of the base mass flux ----'
      call local_flux_candidates

      ! ---------------------------------------------------------------- !
      ! R48: the boundary is a function of its argument
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- R48 repeated evaluation and call order ----'
      call repeated_and_reordered

      write(*,'(A)') ''
      if (n_fail .gt. 0) then
         write(*,'(A,I0)') 'base_boundary_continuity_probe FAILURES: ',   &
              n_fail
         stop 1
      endif
      stop 0

      contains

      !------------------------------------------!

      subroutine build_state(path, t)
      ! Wq is the probed state: a copy of the loaded state with the
      ! velocities of the window, or of the base cells, replaced.  Density
      ! and pressure are never touched, so the composition, the reservoir and
      ! the interior isentrope are those of the loaded state throughout and
      ! only the discriminant of the branch moves.
      integer, intent(in) :: path
      real*8,  intent(in) :: t
      integer :: jj, k
      real*8  :: a, x, pi

      pi = 4.0d0*atan(1.0d0)
      Wq = W0

      select case (path)
      case (0)
         ! the loaded state itself
      case (1, 21, 22, 23, 24, 3)
         do jj = j0, N
            k = jj - j0
            select case (path)
            case (1)
               a = 1.0d0
            case (21)
               a = 1.0d0 + 0.5d0*sin(2.0d0*pi*dble(k)/dble(nw))
            case (22)
               a = 1.0d0
               if (mod(k,2) .eq. 1) a = 2.0d0
            case (23)
               a = 0.5d0 + dble(k)/dble(max(nw-1,1))
            case (24)
               a = 1.0d0 + 2.0d-3*sin(2.0d0*pi*dble(k)/dble(nw))
            case (3)
               ! alternating signs; an odd cell out is set to zero so that
               ! the window mean is zero by construction and not by
               ! cancellation of unequal counts
               a = 1.0d0
               if (mod(k,2) .eq. 1) a = -1.0d0
               if (mod(nw,2) .eq. 1 .and. jj .eq. N) a = 0.0d0
            end select
            Wq(2,jj) = t*F0*a/(Wq(1,jj)*r(jj)*r(jj))
         enddo
      case (41, 42)
         x = 1.0d0
         if (path .eq. 42) x = -1.0d0
         do jj = j0, N
            Wq(2,jj) = x*F0/(Wq(1,jj)*r(jj)*r(jj))
         enddo
         Wq(2,1) = t*abs(W0(2,1))
      case (5)
         do jj = j0, N
            a = 1.0d0
            if (jj .eq. jA) a = 1.0d0 + alpha_extremum
            if (jj .eq. jB) a = 1.0d0 + t
            Wq(2,jj) = F0*a/(Wq(1,jj)*r(jj)*r(jj))
         enddo
      case (6)
         do jj = 1, min(6,N)
            Wq(2,jj) = t*W0(2,jj)
         enddo
      case (7)
         Wq(2,:) = 0.0d0
      case (71)
         Wq(2,:) = 0.0d0
         do jj = j0, N
            Wq(2,jj) = t*F0/(Wq(1,jj)*r(jj)*r(jj))
         enddo
      case (73)
         ! the same static base under a window whose spread is above twice
         ! the gate, so that the static state is approached along a closed
         ! gate as well as along an open one
         Wq(2,:) = 0.0d0
         do jj = j0, N
            k = jj - j0
            a = 1.0d0 + 0.5d0*sin(2.0d0*pi*dble(k)/dble(nw))
            Wq(2,jj) = t*F0*a/(Wq(1,jj)*r(jj)*r(jj))
         enddo
      case (72)
         do jj = 1, min(6,N)
            Wq(2,jj) = 0.0d0
         enddo
      end select
      end subroutine build_state

      !------------------------------------------!

      real*8 function rho_b_of(path, t)
      ! The base face density the production boundary returns on the state
      ! this path and parameter build.
      integer, intent(in) :: path
      real*8,  intent(in) :: t
      call build_state(path, t)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      rho_b_of = Wface(1)
      end function rho_b_of

      !------------------------------------------!

      subroutine sample(tag, path, t)
      ! One sampled point, with everything the branch was decided by.
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: path
      real*8,           intent(in) :: t
      real*8  :: Fw, duw, w_i
      logical :: hv
      character(len=*), parameter :: fmt_row =                            &
           '(A,A,A,ES13.6,A,ES22.15,A,ES22.15,A,ES22.15,'             //  &
           'A,ES13.6,A,ES13.6,'                                     //  &
           'A,F10.7,A,F10.7,A,F10.7,A,ES13.6,A,L1,A,ES16.9,A,ES16.9)'
      call build_state(path, t)
      call wind_window_mass_flux(Wq, Fw, duw, hv)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      w_i = characteristic_branch_weight(base_face_Mi_last                &
                                         /base_face_mach_blend)
      write(*,fmt_row)                                                    &
           ' ROW ', tag, ' t=', t,                                        &
           ' rho_b=', Wface(1),                                           &
           ' gh0=', Wghost(1,0),                                          &
           ' gh1=', Wghost(1,1-Ng),                                       &
           ' M_i=', base_face_Mi_last,                                    &
           ' M_wind=', base_face_Mwind_last,                               &
           ' w_i=', w_i,                                                  &
           ' s=', base_face_swind_last,                                   &
           ' w_rev=', base_face_blend_last,                               &
           ' du=', duw,                                                   &
           ' haveF=', hv,                                                 &
           ' rho_res=', base_face_rho_res_last,                           &
           ' rho_rev=', base_face_rho_rev_last
      end subroutine sample

      !------------------------------------------!

      subroutine decade_sweep(tag, path, sgn)
      ! eps = sgn * 1e-2, 1e-4, ... 1e-14 times the window flux of the state.
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: path
      real*8,           intent(in) :: sgn
      integer :: k
      do k = 2, 14, 2
         call sample(tag, path, sgn*1.0d1**(-dble(k)))
      enddo
      end subroutine decade_sweep

      !------------------------------------------!

      subroutine sign_sweep(tag, path)
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: path
      integer :: k
      real*8  :: t
      do k = -2, 2
         t = 0.5d0*dble(k)
         call sample(tag, path, t)
      enddo
      end subroutine sign_sweep

      !------------------------------------------!

      subroutine extremum_sweep(tag, path)
      ! t on both sides of the fixed bump, so the cell that carries the
      ! window maximum changes inside the sweep.
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: path
      integer :: k
      do k = -5, 5
         call sample(tag, path, alpha_extremum + dble(k)*1.0d-3)
      enddo
      end subroutine extremum_sweep

      !------------------------------------------!

      subroutine derivatives(tag, path, t0, sgn)
      ! The one-sided directional derivative of rho_b at t0 in the direction
      ! sgn, by forward differences at three step sizes a decade apart.  A
      ! jump at t0 makes the estimates grow like 1/h (each is ten times the
      ! last); a kink makes the two directions converge to different finite
      ! numbers; a smooth point makes all three agree.
      character(len=*), intent(in) :: tag
      integer,          intent(in) :: path
      real*8,           intent(in) :: t0, sgn
      real*8  :: h(3), d(3), r0
      integer :: k
      h = (/ 1.0d-6, 1.0d-7, 1.0d-8 /)
      if (path .eq. 1 .or. path .eq. 21 .or. path .eq. 22 .or.            &
          path .eq. 23 .or. path .eq. 24 .or. path .eq. 3  .or.           &
          path .eq. 71 .or. path .eq. 73)                                &
         h = (/ 1.0d-8, 1.0d-10, 1.0d-12 /)
      r0 = rho_b_of(path, t0)
      do k = 1, 3
         d(k) = (rho_b_of(path, t0 + sgn*h(k)) - r0)/(sgn*h(k))
      enddo
      write(*,'(A,A,A,ES12.5,A,F5.1,3(A,ES12.5,A,ES13.6))')              &
           ' DERIV ', tag, ' t0=', t0, ' dir=', sgn,                      &
           '   h=', h(1), ' D=', d(1),                                    &
           '   h=', h(2), ' D=', d(2),                                    &
           '   h=', h(3), ' D=', d(3)
      write(*,'(A,ES22.15,A,ES12.5,A,ES12.5)')                            &
           '   rho_b(t0)=', r0, '  |D2/D1|=',                             &
           abs(d(2))/max(abs(d(1)),1.0d-300), '  |D3/D2|=',               &
           abs(d(3))/max(abs(d(2)),1.0d-300)
      d_last = d
      h_last = h
      rho_last = r0
      end subroutine derivatives

      !------------------------------------------!

      subroutine derivative_converges(name)
      ! A DIRECTIONAL DERIVATIVE THAT EXISTS.  The three difference quotients
      ! of the last `derivatives` call are taken at steps two decades apart,
      ! and the row asks whether the quotient GROWS as the step shrinks.
      !
      ! That is the one thing a jump does and nothing else does: across a gap
      ! g the quotient is g/h, so it is multiplied by a hundred at every step
      ! of this sweep.  A derivative that exists leaves the quotient bounded:
      ! where it is nonzero the three agree, and where it is zero they fall
      ! like h, which is a C1 point with a large second derivative and not a
      ! discontinuity.  Testing the three for RELATIVE agreement would refuse
      ! that second case however smooth the face density is, so the verdict is
      ! taken on the growth and the relative agreement is printed beside it.
      !
      ! A quotient below the rounding floor of its own difference,
      ! 4 eps rho_b / h, is one unit in the last place of the face density
      ! divided by the step and is not a measurement of anything; it is read
      ! as zero, or the last step of a face density that does not move at all
      ! would read as unbounded growth.
      character(len=*), intent(in) :: name
      real*8 :: g1, g2, growth, de(3)
      integer :: k
      do k = 1, 3
         de(k) = abs(d_last(k))
         if (de(k) .le. 4.0d0*epsilon(1.0d0)*abs(rho_last)/h_last(k))     &
            de(k) = 0.0d0
      enddo
      g1 = de(2)/(de(1) + 1.0d-300)
      g2 = de(3)/(de(2) + 1.0d-300)
      growth = max(g1, g2)
      if (de(1) + de(2) + de(3) .eq. 0.0d0) growth = 0.0d0
      call verdict_le(name, growth, 1.0d0 + 1.0d-3)
      end subroutine derivative_converges

      !------------------------------------------!

      subroutine zero_window_limit(name, path, sgn)
      ! THE LIMIT AT THE ZERO WINDOW IS THE VALUE AT IT.  The face density at
      ! the smallest sampled amplitude of a path against the face density at
      ! amplitude zero exactly, relative to the face density itself.  This is
      ! the row that read 37.4 per cent of rho_b before item L26's repair.
      character(len=*), intent(in) :: name
      integer,          intent(in) :: path
      real*8,           intent(in) :: sgn
      real*8 :: r_lim, r_zero
      r_lim  = rho_b_of(path, sgn*1.0d-14)
      r_zero = rho_b_of(path, 0.0d0)
      call verdict_rel(name, abs(r_lim - r_zero)/max(abs(r_lim),1.0d-99), &
                       0.0d0, 1.0d-12)
      end subroutine zero_window_limit

      !------------------------------------------!

      subroutine jump_report(what, path)
      ! The size of the gap between the limit of a path at t -> 0 and the
      ! value the boundary returns at t = 0 exactly, in the units the memo
      ! quotes: relative to rho_b itself, and against the bound
      ! |w_i - 1/2| |rho_res - rho_rev|, which is what the gap has to be if
      ! the whole of it comes from the branch weight jumping from 1/2 to w_i.
      character(len=*), intent(in) :: what
      integer,          intent(in) :: path
      real*8 :: r_at_zero, r_limit, w_i, bound
      r_limit   = rho_b_of(path, 1.0d-14)
      r_at_zero = rho_b_of(path, 0.0d0)
      w_i   = characteristic_branch_weight(base_face_Mi_last              &
                                           /base_face_mach_blend)
      bound = abs(w_i - 0.5d0)*abs(base_face_rho_res_last                 &
                                   - base_face_rho_rev_last)
      write(*,'(A,A)') ' JUMP ', what
      write(*,'(A,ES22.15,A,ES22.15)') '   rho_b(t -> 0) =', r_limit,     &
           '   rho_b(t = 0) =', r_at_zero
      write(*,'(A,ES16.9,A,F12.8,A)') '   gap =',                         &
           r_limit - r_at_zero, '  =', 1.0d2*abs(r_limit - r_at_zero)     &
           /max(abs(r_limit),1.0d-99), ' per cent of rho_b(t -> 0)'
      write(*,'(A,ES16.9,A,ES13.6)') '   |w_i - 1/2| |rho_res - rho_rev|' &
           //' =', bound, '   w_i =', w_i
      end subroutine jump_report

      !------------------------------------------!

      subroutine local_flux_candidates
      ! CAN THE DIRECTION AT THE BASE FACE BE READ LOCALLY?  Every reading
      ! that is a function of the base cells alone, printed in units of the
      ! wind window's mean flux, beside the mass flux the Riemann solve
      ! actually puts through the base face on the same state.
      !
      ! The readings are (a) the cell-centred product of cell 1, which is the
      ! discriminant w_i is built from; (b) its mean over cells 1..k, the
      ! local average that annihilates a Nyquist component of the velocity;
      ! (c) a least-squares linear extrapolation of the same product to
      ! r_edg(0) over cells 1..n; (d) the face state's own rho_b v_b r_b^2.
      ! A reading is admissible as a discriminant only if its SIGN is the
      ! sign of the flux the face carries.
      real*8  :: Fbase, ssum, sx, sy, sxx, sxy, aa, bb, xx
      integer :: j, k
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      call W_to_U(Wq, u)
      call Reconstruct(u, WLr, WRr)
      call RK_rhs(u, WLr, WRr, dFr, Sr)
      Fbase = face_flux(1,0)*r_edg(0)*r_edg(0)
      write(*,'(A,ES16.9)') '   window mean flux F0 =', F0
      write(*,'(A,F16.6)') '   CAND face_flux_at_the_base_face /F0 =',    &
           Fbase/F0
      do j = 0, 6
         write(*,'(A,I2,A,F16.6)') '   face ', j,                          &
              '  Riemann flux /F0 =', face_flux(1,j)*r_edg(j)*r_edg(j)/F0
      enddo
      do j = 1, 6
         write(*,'(A,I2,A,F16.6,A,ES13.6)') '   cell ', j,                 &
              '  rho v r^2 /F0 =',                                         &
              Wq(1,j)*Wq(2,j)*r(j)*r(j)/F0, '   v [cm/s] =', Wq(2,j)*v0
      enddo
      do k = 1, 6
         ssum = 0.0d0
         do j = 1, k
            ssum = ssum + Wq(1,j)*Wq(2,j)*r(j)*r(j)
         enddo
         write(*,'(A,I2,A,F16.6)') '   CAND mean_of_cells_1_to_', k,       &
              ' /F0 =', ssum/dble(k)/F0
      enddo
      do k = 2, 6
         sx = 0.0d0; sy = 0.0d0; sxx = 0.0d0; sxy = 0.0d0
         do j = 1, k
            xx = r(j) - r_edg(0)
            sx = sx + xx;  sy = sy + Wq(1,j)*Wq(2,j)*r(j)*r(j)
            sxx = sxx + xx*xx
            sxy = sxy + xx*Wq(1,j)*Wq(2,j)*r(j)*r(j)
         enddo
         bb = (dble(k)*sxy - sx*sy)/(dble(k)*sxx - sx*sx)
         aa = (sy - bb*sx)/dble(k)
         write(*,'(A,I2,A,F16.6)') '   CAND lsq_extrapolation_1_to_', k,   &
              ' /F0 =', aa/F0
      enddo
      write(*,'(A,F16.6)') '   CAND face_state_rho_b_v_b_rb2      /F0 =', &
           Wface(1)*Wface(2)*r_edg(0)*r_edg(0)/F0
      end subroutine local_flux_candidates

      !------------------------------------------!

      subroutine repeated_and_reordered
      ! THE BOUNDARY IS A FUNCTION OF ITS ARGUMENT, so evaluating it ten
      ! times on one state, evaluating it on another state in between, and
      ! evaluating it with a complete flux evaluation of the state in
      ! between, must all return the same bits.  A stale face flux or any
      ! other saved quantity that had become an argument of the condition
      ! would break one of these.
      real*8  :: rho_ref, rho_k, rho_after, rho_other
      integer :: k, n_diff

      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      rho_ref = Wface(1)

      n_diff = 0
      do k = 1, 10
         call build_state(0, 0.0d0)
         call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
         if (Wface(1) .ne. rho_ref) n_diff = n_diff + 1
      enddo
      call verdict_int('base_boundary_repeated_evaluation_identical',     &
                       n_diff, 0)

      ! an evaluation on another state in between
      rho_other = rho_b_of(1, 1.0d-6)
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      rho_after = Wface(1)
      write(*,'(A,ES22.15)') '   the intervening state gave rho_b=',      &
           rho_other
      call verdict_bit('base_boundary_after_another_state_identical',     &
                       rho_after, rho_ref)

      ! a complete flux evaluation of the same state in between: the
      ! reconstruction and the Riemann solve this boundary feeds
      call build_state(0, 0.0d0)
      call W_to_U(Wq, u)
      call Reconstruct(u, WLr, WRr)
      call RK_rhs(u, WLr, WRr, dFr, Sr)
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      rho_k = Wface(1)
      call verdict_bit('base_boundary_after_flux_evaluation_identical',   &
                       rho_k, rho_ref)

      ! and the boundary evaluated BEFORE that flux evaluation, in a fresh
      ! order: boundary, flux, boundary, flux, boundary
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      rho_k = Wface(1)
      call W_to_U(Wq, u)
      call Reconstruct(u, WLr, WRr)
      call RK_rhs(u, WLr, WRr, dFr, Sr)
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      call W_to_U(Wq, u)
      call Reconstruct(u, WLr, WRr)
      call RK_rhs(u, WLr, WRr, dFr, Sr)
      call build_state(0, 0.0d0)
      call base_boundary_states(Wq, Wface, Wghost, Wface_lower)
      call verdict_bit('base_boundary_call_order_identical',              &
                       Wface(1), rho_k)
      end subroutine repeated_and_reordered

      !------------------------------------------!

      subroutine verdict_le(name, measured, bound)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, bound
      if (measured .le. bound) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'PASS ', name,              &
              ' measured=', measured, ' reference=', bound,               &
              ' tol=upper_bound'
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'FAIL ', name,              &
              ' measured=', measured, ' reference=', bound,               &
              ' tol=upper_bound'
         n_fail = n_fail + 1
      endif
      end subroutine verdict_le

      !------------------------------------------!

      subroutine verdict_rel(name, measured, reference, tol)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, reference, tol
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES10.3)') 'PASS ', name,       &
              ' measured=', measured, ' reference=', reference, ' tol=',  &
              tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES10.3)') 'FAIL ', name,       &
              ' measured=', measured, ' reference=', reference, ' tol=',  &
              tol
         n_fail = n_fail + 1
      endif
      end subroutine verdict_rel

      !------------------------------------------!

      subroutine verdict_bit(name, measured, reference)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, reference
      if (measured .eq. reference) then
         write(*,'(A,A,A,ES22.15,A,ES22.15,A)') 'PASS ', name,            &
              ' measured=', measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,ES22.15,A,ES22.15,A)') 'FAIL ', name,            &
              ' measured=', measured, ' reference=', reference, ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine verdict_bit

      !------------------------------------------!

      subroutine verdict_int(name, measured, reference)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: measured, reference
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine verdict_int

      end program base_boundary_continuity_probe

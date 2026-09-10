      program positivity_limiter_scaling
      ! The positivity scaling of the reconstructed face states, against the
      ! limiter it cites: Zhang and Shu (2010, J. Comput. Phys., 229, 8918;
      ! doi:10.1016/j.jcp.2010.08.016), Section 2.2, their Eqs. (2.4), (2.5),
      ! (2.10), (2.11) and the quadratic (2.12).
      !
      ! The routine under test is the production one,
      ! Reconstruction_step's positivity_limited_faces, together with
      ! positivity_scaling and positive_variable_scaling.
      !
      ! WHAT THE PAPER DOES.  Its limiter acts on the CONSERVED vector
      ! w = (rho, m, E)^T.  Step one scales the density polynomial about the
      ! cell average (their Eq. 2.4) by
      !     theta_1 = min{ (rho_j^n - eps)/(rho_j^n - rho_min), 1 }   (2.5),
      ! step two scales the whole vector about the cell average (their
      ! Eq. 2.10) by theta_2 = min_a |s_eps^a - w_j^n| / |qhat_j^a - w_j^n|
      ! (their Eq. 2.11), which their implementation section states is
      ! mathematically equivalent to theta_2 = min_a t_eps^a with t_eps^a the
      ! root of the quadratic
      !     p[ (1 - t_eps^a) w_j^n + t_eps^a qhat_j(xhat_j^a) ] = eps  (2.12).
      ! Their floor is absolute: "For example, we can take eps = 10^-13 in
      ! the computation", tightened in the flowchart to
      ! eps = min_j { 10^-13, rho_j^n, p(w_j^n) }.
      !
      ! WHAT THIS CODE DOES, AND WHERE IT DEPARTS.  It scales the PRIMITIVE
      ! vector W = (rho, v, p)^T about the cell average, in ONE step with
      ! theta = min over rho and p, and with the relative floor
      ! eps = epsilon(1)*q_avg rather than an absolute 10^-13.  The three
      ! consequences measured below are:
      !
      !   (1) The failure mode a conserved-space limiter has to solve a
      !       quadratic to avoid does not exist here.  Pressure is itself a
      !       scaled variable, so the limited face carries the pressure the
      !       scaling put there and its internal energy p/(gamma-1) has the
      !       same sign; there is no path on which the scaled state has a
      !       negative internal energy.  Assertions
      !       primitive_scaling_pressure_positive and
      !       primitive_scaling_internal_energy.
      !
      !       What the floor does NOT survive is round-off.  eps is one unit
      !       in the last place of the cell average, so the scaled value is a
      !       cancellation of the same size as the floor and lands on EXACTLY
      !       ZERO for about a part in 200 of the states that fire the
      !       limiter; a zero face density then divides in v = m/rho and in
      !       the sound speed sqrt(gamma p/rho) of Num_Fluxes.f90.  Section
      !       (2b) sweeps it, and shows that clamping the scaled value to the
      !       floor removes it.  Section (5) is the same fact seen from the
      !       conserved variables: at Mach 60 a floor of one ulp of p is
      !       nineteen decades below E, so the paper's admissible set of
      !       their Eq. (2.6) cannot be tested at it at all.
      !
      !   (2) The two paths are NOT the same state.  Along the paper's path
      !       rho, m and E are linear in t, so
      !       p(t) = (gamma-1)(E(t) - m(t)^2/(2 rho(t))) is CONCAVE in t
      !       (m^2/rho is the perspective function of m^2 and is jointly
      !       convex), whence p(t) >= (1-t) p_avg + t p_rec, which is exactly
      !       the pressure this code interpolates.  So the code's theta is
      !       never larger than the paper's t_eps of Eq. (2.12): at equal eps
      !       this limiter is the more restrictive of the two, never the
      !       less safe one.  Assertion primitive_theta_at_most_conserved,
      !       with both numbers printed.
      !
      !   (3) Conservativity, the third property their Section 2.2 requires
      !       of the limited polynomial, is lost.  Their single theta_2 per
      !       CELL leaves the cell average of the limited polynomial equal to
      !       w_j^n; this code takes a theta per FACE STATE, so a cell whose
      !       two ends are scaled by different factors no longer reconstructs
      !       to its own average.  Assertion
      !       endpoint_mean_departs_from_average measures the departure, and
      !       cell_theta_restores_average shows that the paper's min over the
      !       cell's two endpoint values puts it back at round-off.
      !
      ! WHAT OF THEIR THEOREM CARRIES OVER.  The scaling is a convex
      ! combination with an admissible anchor, so positivity of the face
      ! values carries over (that is Lemma 2.5, and it needs no quadrature).
      ! Their Theorem 2.1 does not: it concludes that the NEXT cell average
      ! lies in G, and its hypotheses are the N-point Legendre Gauss-Lobatto
      ! point set of their Eq. (1.7) with 2N - 3 >= k, and the CFL condition
      ! of their Eq. (2.1), lambda*max(|u| + c) <= what_1*alpha_0.  A
      ! finite-volume reconstruction to two face values supplies only the
      ! N = 2 point set, admissible for k = 1 (PLM) and not for the k = 2 of
      ! WENO3, which their Remark 2.6 names as the open implementation
      ! question for finite-volume WENO; and this code imposes no such CFL
      ! bound on the flux it then evaluates.
      !
      ! Reference state.  A hypersonic layer, the state the limiter fires in
      ! (Mach 60): the thermal pressure is a part in 3e3 of the total energy
      ! density, so a reconstruction across a jump tips it negative.

      use global_parameters, only: N, Ng, gamma_ad
      use caloric_eos, only: caloric_mixture_active
      use Reconstruction_step, only: positivity_limited_faces,            &
                                     positivity_scaling,                  &
                                     positive_variable_scaling,           &
                                     n_faces_positivity_limited
      use Conversion, only: W_to_U, U_to_W
      use assertion_report
      implicit none

      real*8, allocatable :: u_avg(:,:), w_avg(:,:)
      real*8, allocatable :: wl(:,:), wr(:,:), wl0(:,:), wr0(:,:)
      real*8 :: w_a(3), w_rec(3), u_a(3), u_rec(3)
      real*8 :: th_code, t_eps, th_left, th_right, th_cell
      real*8 :: eps_cmp, th_cmp
      real*8 :: rho_lim, v_lim, p_lim, e_int, p_round_trip
      real*8 :: mean_rho, mean_p, dep_face, dep_cell
      real*8 :: wl_end(3), wr_end(3)
      integer :: j, j0
      integer :: n_zero, n_neg, n_swept, n_zero_clamped

      ! The atomic, single-gamma configuration: the pressure of a face is
      ! (gamma_ad - 1) times its internal energy density, so the conserved
      ! path below is the one the code's own conversion would take.
      caloric_mixture_active = .false.
      N = 5
      j0 = 3

      allocate(u_avg(3,1-Ng:N+Ng), w_avg(3,1-Ng:N+Ng))
      allocate(wl(3,1-Ng:N+Ng), wr(3,1-Ng:N+Ng))
      allocate(wl0(3,1-Ng:N+Ng), wr0(3,1-Ng:N+Ng))

      ! Cell averages: uniform hypersonic layer, admissible everywhere.
      w_a = (/ 1.0d-4, 6.0d6, 6.0d5 /)
      do j = 1-Ng, N+Ng
         w_avg(:,j) = w_a
      enddo
      call W_to_U(w_avg, u_avg)

      ! --- (0) an admissible reconstruction is left bitwise alone ------ !
      do j = 1-Ng, N+Ng
         wl(:,j) = w_a*(/ 1.05d0, 1.01d0, 1.10d0 /)
         wr(:,j) = w_a*(/ 0.95d0, 0.99d0, 0.90d0 /)
      enddo
      wl0 = wl
      wr0 = wr
      n_faces_positivity_limited = 0
      call positivity_limited_faces(u_avg, wl, wr)
      call check_absolute('admissible_reconstruction_untouched',          &
                          maxval(abs(wl - wl0)) + maxval(abs(wr - wr0)),  &
                          0.0d0, 0.0d0)
      call check_absolute('admissible_reconstruction_count',              &
                          dble(n_faces_positivity_limited), 0.0d0, 0.0d0)

      ! --- (1) a negative reconstructed pressure ---------------------- !
      ! The left state of face j0 belongs to cell j0 and is anchored on that
      ! cell's average; the right state of face j0 belongs to cell j0+1.
      w_rec = (/ 1.2d-4, 6.2d6, -1.0d5 /)
      do j = 1-Ng, N+Ng
         wl(:,j) = w_a
         wr(:,j) = w_a
      enddo
      wl(:,j0) = w_rec
      n_faces_positivity_limited = 0
      call positivity_limited_faces(u_avg, wl, wr)

      rho_lim = wl(1,j0)
      v_lim   = wl(2,j0)
      p_lim   = wl(3,j0)
      e_int   = p_lim/(gamma_ad - 1.0d0)

      call check_positive('primitive_scaling_pressure_positive', p_lim)
      call check_positive('primitive_scaling_density_positive', rho_lim)
      call check_positive('primitive_scaling_internal_energy', e_int)
      call check_absolute('limited_face_count',                          &
                          dble(n_faces_positivity_limited), 1.0d0, 0.0d0)

      ! The scaling the code chose, and its closed form: pressure is the
      ! binding variable here, so theta = (p_a - eps)/(p_a - p_rec).
      th_code = positivity_scaling(w_rec(1), w_rec(3), w_a(1), w_a(3))
      call check_relative('theta_is_pressure_limited', th_code,           &
                          (w_a(3) - epsilon(1.0d0)*w_a(3))               &
                          /(w_a(3) - w_rec(3)), 1.0d-14)
      call check_relative('limited_state_is_convex_combination',          &
                          rho_lim,                                        &
                          w_a(1) + th_code*(w_rec(1) - w_a(1)), 1.0d-14)

      ! --- (2) the paper's t_eps of Eq. (2.12) ------------------------ !
      ! The floor used for the comparison is NOT the code's one ulp: along
      ! the conserved path E and m^2/(2 rho) are each 2e9 here and their
      ! difference at the root would be 1e-10, nineteen decades below them,
      ! so p(t) is not resolvable there in double precision and neither
      ! Eq. (2.12) nor a bisection on it means anything at that floor.  The
      ! comparison is therefore made at a floor the state can carry.  The
      ! cancellation in p(t) is of order epsilon(1)*E = 5e-7 in pressure
      ! units here, so eps = 1e-3 p_avg is six decades clear of it and is
      ! still deep enough that the two limiters differ.
      eps_cmp = 1.0d-3*w_a(3)
      call W_to_U_state(w_a,   u_a)
      call W_to_U_state(w_rec, u_rec)
      t_eps   = conserved_path_pressure_root(u_a, u_rec, eps_cmp)
      th_cmp  = (w_a(3) - eps_cmp)/(w_a(3) - w_rec(3))
      call check_relative('conserved_path_root_is_a_root',                &
                          conserved_path_pressure(u_a, u_rec, t_eps),     &
                          eps_cmp, 1.0d-8)
      call check_at_most('primitive_theta_at_most_conserved', th_cmp,     &
                         t_eps)
      write(*,'(a,es23.15,a,es23.15)')                                    &
           'MEASURED theta_primitive_scaling=', th_cmp,                   &
           ' t_eps_Zhang_Shu_2010_Eq_2_12=', t_eps

      ! --- (2b) the floor the scaling actually delivers --------------- !
      ! WHY THE PRODUCTION UPDATE CLAMPS.  theta puts the scaled variable on
      ! the floor eps = epsilon(1)*q_avg only in exact arithmetic.  eps is
      ! one unit in the last place of q_avg, so the bare update
      ! q_avg + theta (q_rec - q_avg) is a cancellation whose result is
      ! decided by the rounding of the product, not by theta: it scatters
      ! over roughly +-1 eps, never negative in the sweep below, but EXACTLY
      ! ZERO for a part in 197 of the reconstructions that fire the limiter,
      ! the value the routine exists to avoid (a zero face density divides
      ! in v = m/rho and in the sound speed sqrt(gamma p/rho) of
      ! Num_Fluxes.f90).  positivity_limited_faces therefore clamps the
      ! scaled density and pressure to eps after the update (2026-09-07);
      ! the sweep asserts both halves: the bare update does reach zero, and
      ! the production routine, driven over the same reconstructed values,
      ! never delivers a face density or pressure at or below zero.
      call floor_delivery_sweep(w_a(3), n_swept, n_zero, n_neg,           &
                                n_zero_clamped)
      call check_at_least('bare_update_never_negative',                   &
                          dble(n_swept - n_neg), dble(n_swept))
      call check_at_least('bare_update_reaches_zero',                     &
                          dble(n_zero), 1.0d0)
      call check_absolute('clamped_scaling_delivers_the_floor',           &
                          dble(n_zero_clamped), 0.0d0, 0.0d0)
      write(*,'(a,i0,a,i0,a,i0)')                                         &
           'MEASURED floor_sweep_states=', n_swept,                       &
           ' bare_update_exactly_zero=', n_zero, ' bare_update_negative=', n_neg
      call production_floor_sweep(u_avg, w_avg, wl, wr, n_swept, n_zero)
      call check_absolute('production_limiter_never_delivers_zero',       &
                          dble(n_zero), 0.0d0, 0.0d0)
      write(*,'(a,i0,a,i0)')                                              &
           'MEASURED production_sweep_face_states=', n_swept,             &
           ' at_or_below_zero=', n_zero

      ! --- (3) conservativity of the limited reconstruction ----------- !
      ! Cell j0 with BOTH ends reconstructed: its left end is the right
      ! state of face j0-1, its right end is the left state of face j0.  For
      ! a k = 1 reconstruction the two-point Gauss-Lobatto rule has weights
      ! 1/2, 1/2, so the unlimited pair averages to the cell average, and
      ! the paper's limited pair must do so too.
      wl_end = w_a + (/ 0.2d-4, 0.2d6, -7.0d5 /)
      wr_end = w_a - (/ 0.2d-4, 0.2d6, -7.0d5 /)
      th_right = positivity_scaling(wl_end(1), wl_end(3), w_a(1), w_a(3))
      th_left  = positivity_scaling(wr_end(1), wr_end(3), w_a(1), w_a(3))
      th_cell  = min(th_left, th_right)

      ! As the code scales, one theta per face state:
      mean_rho = 0.5d0*( (w_a(1) + th_left *(wr_end(1) - w_a(1)))         &
                       + (w_a(1) + th_right*(wl_end(1) - w_a(1))) )
      mean_p   = 0.5d0*( (w_a(3) + th_left *(wr_end(3) - w_a(3)))         &
                       + (w_a(3) + th_right*(wl_end(3) - w_a(3))) )
      dep_face = max(abs(mean_rho/w_a(1) - 1.0d0),                        &
                     abs(mean_p  /w_a(3) - 1.0d0))

      ! As Eq. (2.11) scales, one theta per cell:
      mean_rho = 0.5d0*( (w_a(1) + th_cell*(wr_end(1) - w_a(1)))          &
                       + (w_a(1) + th_cell*(wl_end(1) - w_a(1))) )
      mean_p   = 0.5d0*( (w_a(3) + th_cell*(wr_end(3) - w_a(3)))          &
                       + (w_a(3) + th_cell*(wl_end(3) - w_a(3))) )
      dep_cell = max(abs(mean_rho/w_a(1) - 1.0d0),                        &
                     abs(mean_p  /w_a(3) - 1.0d0))

      call check_at_least('endpoint_mean_departs_from_average', dep_face, &
                          1.0d-2)
      call check_at_most('cell_theta_restores_average', dep_cell, 1.0d-14)
      write(*,'(a,es23.15,a,es23.15)')                                    &
           'MEASURED endpoint_mean_departure_face_theta=', dep_face,      &
           ' cell_theta=', dep_cell

      ! --- (4) a NaN face state falls back on the cell average -------- !
      call check_absolute('nan_face_scales_to_the_average',               &
                          positive_variable_scaling(nan_value(), w_a(1)), &
                          0.0d0, 0.0d0)

      ! --- (5) the ulp floor in the paper's conserved variables ------- !
      ! Their admissible set G_eps (their Eq. 2.6) is stated on the CONSERVED
      ! vector, and a scheme that carries the state as (rho, m, E) tests
      ! p = (gamma-1)(E - m^2/(2 rho)) there.  At the code's floor of one ulp
      ! of the cell average that test cannot be performed: mapping the
      ! limited face to conserved variables and back returns a pressure that
      ! bears no relation to the one that was limited, because E and
      ! m^2/(2 rho) agree to nineteen decimals.  This is not a defect of the
      ! Riemann path, which takes the face primitives directly, but it is the
      ! reason a one-ulp floor is not a floor in the sense the paper means.
      call W_to_U_state((/ rho_lim, v_lim, p_lim /), u_rec)
      p_round_trip = (gamma_ad - 1.0d0)                                   &
                     *(u_rec(3) - 0.5d0*u_rec(2)**2/u_rec(1))
      call check_at_least('ulp_floor_lost_in_conserved_round_trip',       &
                          abs(p_round_trip - p_lim), p_lim)
      write(*,'(a,es23.15,a,es23.15)')                                    &
           'MEASURED limited_face_pressure=', p_lim,                      &
           ' recovered_from_conserved=', p_round_trip

      deallocate(u_avg, w_avg, wl, wr, wl0, wr0)

      if (assertion_failures .gt. 0) error stop 1

      contains

      !--------------!

      subroutine W_to_U_state(w, u)
      ! One state's primitive to conserved map, gamma_ad, as in the
      ! atomic branch of Conversion's W_to_U.
      real*8, intent(in)  :: w(3)
      real*8, intent(out) :: u(3)
      u(1) = w(1)
      u(2) = w(1)*w(2)
      u(3) = w(3)/(gamma_ad - 1.0d0) + 0.5d0*w(1)*w(2)**2
      end subroutine W_to_U_state

      !--------------!

      real*8 function conserved_path_pressure(ua, ub, t)
      ! p of (1 - t) ua + t ub, the pressure of their Eq. (2.6) along the
      ! straight line s_a(t) of their Eq. (2.8).
      real*8, intent(in) :: ua(3), ub(3), t
      real*8 :: u(3)
      u = (1.0d0 - t)*ua + t*ub
      conserved_path_pressure =                                           &
           (gamma_ad - 1.0d0)*(u(3) - 0.5d0*u(2)**2/u(1))
      end function conserved_path_pressure

      !--------------!

      real*8 function conserved_path_pressure_root(ua, ub, eps) result(t)
      ! The t_eps of their Eq. (2.12), by bisection.  p(t) is concave with
      ! p(0) > eps and p(1) < eps, so the root is unique in (0,1) and
      ! bisection is as good as the quadratic formula and better
      ! conditioned when E is the difference of two nearly equal numbers,
      ! which at Mach 60 it is.
      real*8, intent(in) :: ua(3), ub(3), eps
      real*8 :: lo, hi, mid
      integer :: it
      lo = 0.0d0
      hi = 1.0d0
      do it = 1, 200
         mid = 0.5d0*(lo + hi)
         if (conserved_path_pressure(ua, ub, mid) .ge. eps) then
            lo = mid
         else
            hi = mid
         endif
      enddo
      t = lo
      end function conserved_path_pressure_root

      !--------------!

      subroutine floor_delivery_sweep(q_avg, n_swept, n_zero, n_neg,      &
                                      n_zero_clamped)
      ! The value the production scaling delivers for a range of
      ! reconstructed values that all cross zero, and the value the same
      ! theta delivers when the scaled variable is clamped to the floor.
      ! theta comes from the production positive_variable_scaling; the
      ! update is the one positivity_limited_faces applies.
      real*8,  intent(in)  :: q_avg
      integer, intent(out) :: n_swept, n_zero, n_neg, n_zero_clamped
      real*8  :: q_rec, th, q_lim
      integer :: i
      n_swept = 0;  n_zero = 0;  n_neg = 0;  n_zero_clamped = 0
      do i = 1, 200000
         q_rec = -q_avg*dble(i)*1.0d-3
         th    = positive_variable_scaling(q_rec, q_avg)
         q_lim = q_avg + th*(q_rec - q_avg)
         n_swept = n_swept + 1
         if (q_lim .lt. 0.0d0) n_neg  = n_neg  + 1
         if (q_lim .eq. 0.0d0) n_zero = n_zero + 1
         if (scaled_value_clamped(q_avg, q_rec, th) .le. 0.0d0)           &
              n_zero_clamped = n_zero_clamped + 1
      enddo
      end subroutine floor_delivery_sweep

      !--------------!

      subroutine production_floor_sweep(u_avg, w_avg, wl, wr, n_swept,     &
                                        n_zero)
      ! The production routine driven over reconstructed densities and
      ! pressures that cross zero at every face state of the grid, counting
      ! the face values it hands back at or below zero.
      real*8, intent(in)    :: u_avg(:,:), w_avg(:,:)
      real*8, intent(inout) :: wl(:,:), wr(:,:)
      integer, intent(out)  :: n_swept, n_zero
      real*8, allocatable :: ul(:,:)
      integer :: i, j, k
      allocate(ul(size(u_avg,1), size(u_avg,2)))
      ul = u_avg
      n_swept = 0;  n_zero = 0
      do i = 1, 4000
         wl = w_avg;  wr = w_avg
         do j = 1, size(wl,2)
            wl(1,j) = -w_avg(1,j)*dble(i)*5.0d-2
            wr(1,j) = -w_avg(1,j)*dble(i)*5.0d-2*1.3d0
            wl(3,j) = -w_avg(3,j)*dble(i)*5.0d-2*0.7d0
            wr(3,j) = -w_avg(3,j)*dble(i)*5.0d-2
         enddo
         call positivity_limited_faces(ul, wl, wr)
         do j = 1, size(wl,2)
            do k = 1, 3, 2
               n_swept = n_swept + 2
               if (wl(k,j) .le. 0.0d0) n_zero = n_zero + 1
               if (wr(k,j) .le. 0.0d0) n_zero = n_zero + 1
            enddo
         enddo
      enddo
      deallocate(ul)
      end subroutine production_floor_sweep

      !--------------!

      real*8 function scaled_value_clamped(q_avg, q_rec, th)
      ! The proposed form: the scaled variable, clamped to the floor the
      ! scaling was solved for.  In exact arithmetic it is the same number;
      ! in floating point it is the one that cannot land on zero.
      real*8, intent(in) :: q_avg, q_rec, th
      scaled_value_clamped = max(q_avg + th*(q_rec - q_avg),              &
                                 epsilon(1.0d0)*q_avg)
      end function scaled_value_clamped

      !--------------!

      real*8 function nan_value()
      ! A quiet NaN without a compile-time 0/0.
      real*8 :: z
      z = 0.0d0
      nan_value = z/z
      end function nan_value

      end program positivity_limiter_scaling

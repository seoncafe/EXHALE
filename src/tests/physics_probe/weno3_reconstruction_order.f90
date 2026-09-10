      program weno3_reconstruction_order
      ! Convergence order of the face reconstruction, measured against the
      ! exact face value of a smooth function whose cell averages are the
      ! input.
      !
      ! WHAT IS BEING TESTED. The finite-volume update of RK_rhs.f90 divides
      ! the difference of the two face fluxes of cell j by the SHELL VOLUME
      !     dV_j = ( r_{j+1/2}^3 - r_{j-1/2}^3 ) / 3 ,
      ! so the cell state it advances is the volume average
      !     Q_j = ( int f r^2 dr ) / ( int r^2 dr )   over cell j,
      ! and the state the Riemann solver wants at a face is the POINT value
      ! f(r_{j+1/2}).  A reconstruction is p-th order when the difference
      ! between the two, over a smooth f, falls as h^p.  This program builds
      ! the volume averages by 8-point Gauss-Legendre quadrature (exact for
      ! polynomials of degree 15 in r, so the input carries no error of its
      ! own at the level measured here) and reports L_inf and L_1 of the face
      ! error together with the rate between successive resolutions.
      !
      ! WHY THE VOLUME COORDINATE IS THE RIGHT FRAME. Writing V = r^3/3, the
      ! volume average of cell j is the plain average of f over the interval
      ! [V_{j-1/2}, V_{j+1/2}] of width dV_j, and a face in r is a face in V.
      ! A polynomial reconstruction in V from cell averages in V is therefore
      ! the exact analogue of the Cartesian one, and the geometry
      ! coefficients of weno3_geometry_coefficients are built from exactly
      ! these dV.  Everything below is derived in that frame.
      !
      ! THE ESWENO3 FORM THE CODE USES (Reconstruction.f90).  With
      !     dWp = Q_{j+1} - Q_j ,  dWm = Q_j - Q_{j-1} ,
      ! the right face state of cell j is written as a convex combination of
      ! the two two-cell (linear) candidate stencils,
      !     q^L_{j+1/2} = Q_j + w0 c_p dWp + w1 c_m dWm ,
      !     w0 = S0/(S0 + D1 S1) ,   w1 = D1 S1/(S0 + D1 S1) ,
      ! and the left face state of cell j as
      !     q^R_{j-1/2} = Q_j - w0' c_p' dWp - w1' c_m' dWm ,
      !     w0' = D2 S0/(D2 S0 + S1) ,  w1' = S1/(D2 S0 + S1) .
      ! S0, S1 are the ESWENO smoothness factors 1 + tau^2/beta_k of
      ! Yamaleev and Carpenter (2009), with tau = dWp - dWm and
      ! beta_k = (jump)^2 + dr_j^2.  They tend to 1 as O(h^2) on smooth data,
      ! so in a smooth region the scheme reduces to its LINEAR weights, and
      ! the order in a smooth region is the order of that linear scheme.
      ! D1 and D2 therefore carry the ideal (optimal) weight ratios, and
      ! c_p, c_m the two candidate stencils' own face coefficients.
      !
      ! THE COEFFICIENTS A THIRD-ORDER SCHEME NEEDS.  Let a = dV_{j-1},
      ! b = dV_j, c = dV_{j+1} and put the centre of cell j at the origin of
      ! the V axis, so cell j-1 is centred at m1 = -(a+b)/2 and cell j+1 at
      ! m2 = (b+c)/2.  The two linear candidates interpolate to the right
      ! face of cell j as
      !     Q_j + [ b/(b+c) ] dWp        and     Q_j + [ b/(a+b) ] dWm ,
      ! and to its left face as
      !     Q_j - [ b/(b+c) ] dWp        and     Q_j - [ b/(a+b) ] dWm ,
      ! the same two magnitudes at both faces because both faces are b/2 from
      ! the centre.  The parabola through the three cell averages has
      !     p1 = ( beta_g dWp + alpha_g dWm ) / Det ,
      !     p2 = ( -m1 dWp - m2 dWm ) / Det ,
      !     alpha_g = m2^2 + (c^2 - b^2)/12 ,
      !     beta_g  = m1^2 + (a^2 - b^2)/12 ,
      !     Det     = m2 beta_g - m1 alpha_g ,
      ! and its face values are Q_j + p1 b/2 + p2 b^2/6 (right) and
      ! Q_j - p1 b/2 + p2 b^2/6 (left).  Matching term by term gives the
      ! ideal weights; on a grid uniform in V they are 2/3 and 1/3 at the
      ! right face and 1/3 and 2/3 at the left, so
      !     D1 = w1/w0 = 1/2   and   D2 = w0'/w1' = 1/2 .
      !
      ! WHAT THE CODE USES.  weno3_geometry_coefficients returns
      !     C1(j) = c/(b+c) ,  C2(j) = b/(b+c) ,
      !     D1(j) = c/(a+b) ,  D2(j) = a/(b+c) ,
      ! and Reconstruct takes c_p = C2(j), c_m = C1(j-1), c_p' = C2(j),
      ! c_m' = C1(j-1) since 2026-09-07.  Until that fix it took c_p = C1(j)
      ! and c_m' = C2(j-1), which is what this probe was written to measure.
      !
      ! The two WEIGHT RATIOS are admissible.  D1 and D2 are both 1/2 on a
      ! grid uniform in V and 1/2 + O(h) on a smooth stretched one, and an
      ! O(h) error in a weight ratio costs nothing: the two candidate face
      ! values differ from each other by O(h^2), so moving weight between
      ! them by O(h) moves the face value by O(h^3), which is the order the
      ! scheme is meant to have anyway.
      !
      ! The two STENCIL COEFFICIENTS c_p and c_m' were not, before the fix.
      ! The candidate that interpolates from cells j and j+1 carries the
      ! weight b/(b+c) at either face of cell j; the code applied
      ! C1(j) = c/(b+c) at the right face, which is the OTHER cell's volume
      ! share.  Likewise the candidate from cells j-1 and j carries b/(a+b),
      ! and the code applied C2(j-1) = a/(a+b) at the left face.  Both were
      ! 1/2 + O(h) instead of 1/2 - O(h): an O(h) error in a coefficient that
      ! multiplies a jump of size O(h) is an O(h^2) error in the face value,
      ! so the scheme was second order (MEASURED 2026-09-07: rate 1.997 on a
      ! grid uniform in r, 1.988 on the production grid; 3.000 and 2.991
      ! with the shares corrected).  The other two, c_m = C1(j-1) = b/(a+b)
      ! and c_p' = C2(j) = b/(b+c), were already the right ones.  The swap
      ! was inherited verbatim from ATES v2.0.
      !
      ! The prediction that identified it is kept as a test: on a grid
      ! uniform in the VOLUME coordinate, where c = b = a exactly, C1 = C2
      ! and the swap is invisible, so even the old scheme was third order
      ! there.
      !
      ! The variants below separate the coefficients: weno3_face_states_local
      ! is a copy of the production formula that takes its coefficients as
      ! arguments, checked against the production routine on the production
      ! coefficients before anything is concluded from it; the "corrected"
      ! variant now coincides with production and is kept as the statement of
      ! what the coefficients must be.
      !
      ! References:
      !   N. K. Yamaleev and M. H. Carpenter, "Third-order Energy Stable
      !   WENO scheme", J. Comput. Phys. 228, 3025 (2009),
      !   doi:10.1016/j.jcp.2009.01.011.  The published version was read.
      !   Its Eqs. (18), (21) and (22) are the smoothness factors the
      !   production routine computes, with their tau already a square,
      !   tau = (u_{j+1} - 2 u_j + u_{j-1})^2, and their ideal weights
      !   d_0 = 2/3, d_1 = 1/3 (Eq. 20) entering only as the ratio 1/2 that
      !   D1 and D2 carry.  Their Eq. (24), w_r = d_r + O(Dx^2) on smooth
      !   data, is why the smoothness factors cannot set the order and the
      !   stencil coefficients can, which is what the tables below measure.
      !   Their floor is eps = O(Dxi^2) with the solution scale in front
      !   (Eqs. 62, 64, 65); the eps scan below is the production floor
      !   eps = dr_j^2 scaled by kappa, which spans that prescription.
      !   Their scheme is a uniform-grid finite-difference reconstruction of
      !   the FLUX, and its energy stability rests on the extra dissipation
      !   term of Eq. (48), which the production routine does not compute,
      !   so the volume-coordinate finite-volume coefficients derived above
      !   are this code's own and only the accuracy result carries over.
      !   The attribution in Reconstruction.f90 states this in full.
      !   The companion paper, "A systematic methodology for constructing
      !   high-order energy stable WENO schemes", J. Comput. Phys. 228,
      !   4248 (2009), doi:10.1016/j.jcp.2009.03.002, carries the same
      !   weight form (its Eq. 58) and the same solution-scaled floor (its
      !   Eq. 79) for design orders four and above.

      use global_parameters
      use grid_construction, only: define_grid
      use Reconstruction_step, only: Reconstruct, Reconstruct_scalar,     &
                                     weno3_geometry_coefficients
      use PLM_reconstruction, only: PLM_rec_scalar
      use Conversion, only: W_to_U_comp
      use assertion_report, only: check_absolute, assertion_failures

      implicit none

      integer, parameter :: n_res = 5
      integer :: nn(n_res) = (/ 50, 100, 200, 400, 800 /)

      real*8  :: e_inf(n_res), e_one(n_res)
      real*8  :: rate_prod_sin_u, rate_fix_sin_u
      real*8  :: rate_prod_gau_u, rate_fix_gau_u
      real*8  :: rate_plm_sin_u
      real*8  :: rate_prod_sin_m, rate_fix_sin_m, rate_plm_sin_m
      real*8  :: rate_Conly_sin_m, rate_prod_sin_v
      real*8  :: rate_linD_sin_u, rate_Donly_sin_u, rate_Conly_sin_u
      real*8  :: dev_local, dev_vector, exact_err
      real*8  :: kap(4) = (/ 1.0d-2, 1.0d-1, 1.0d0, 1.0d1 /)
      real*8  :: rate_eps(4)
      integer :: i

      write(*,'(a)') '  DIAGNOSTIC weno3_reconstruction_order: face error'// &
                     ' against the exact face value'

      ! ---- 1. the local copy reproduces the production routine ----
      call local_copy_deviation(200, 1, dev_local)
      call check_absolute('weno3_local_copy_matches_production',           &
                          dev_local, 0.0d0, 1.0d-14)

      ! ---- 2. the vector reconstruction carries the same formula ----
      call vector_scalar_deviation(200, dev_vector)
      call check_absolute('weno3_vector_matches_scalar_faces',             &
                          dev_vector, 0.0d0, 1.0d-13)

      ! ---- 3. order on a grid uniform in r ----
      write(*,'(a)') '  DIAGNOSTIC grid: uniform in r on [1,3],'//         &
                     ' profile: sine'
      call order_table('production   ', 0, 1, 0, e_inf, e_one,             &
                       rate_prod_sin_u)
      call order_table('optimal C,D  ', 0, 1, 3, e_inf, e_one,             &
                       rate_fix_sin_u)
      call order_table('optimal D, production C', 0, 1, 1, e_inf, e_one,   &
                       rate_Donly_sin_u)
      call order_table('corrected C, production D', 0, 1, 2, e_inf, e_one, &
                       rate_Conly_sin_u)
      call order_table('production, S=1', 0, 1, 4, e_inf, e_one,           &
                       rate_linD_sin_u)
      call order_table('PLM          ', 0, 1, -1, e_inf, e_one,            &
                       rate_plm_sin_u)

      write(*,'(a)') '  DIAGNOSTIC grid: uniform in r on [1,3],'//         &
                     ' profile: Gaussian'
      call order_table('production   ', 0, 2, 0, e_inf, e_one,             &
                       rate_prod_gau_u)
      call order_table('optimal C,D  ', 0, 2, 3, e_inf, e_one,             &
                       rate_fix_gau_u)

      ! ---- 4. order on the production radial grid ----
      write(*,'(a)') '  DIAGNOSTIC grid: define_grid Mixed'//             &
                     ' (base cells and geometric stretching), profile: sine'
      call order_table('production   ', 1, 1, 0, e_inf, e_one,             &
                       rate_prod_sin_m)
      call order_table('optimal C,D  ', 1, 1, 3, e_inf, e_one,             &
                       rate_fix_sin_m)
      call order_table('corrected C, production D', 1, 1, 2, e_inf, e_one, &
                       rate_Conly_sin_m)
      call order_table('PLM          ', 1, 1, -1, e_inf, e_one,            &
                       rate_plm_sin_m)

      write(*,'(a)') '  DIAGNOSTIC grid: uniform in the volume'//          &
                     ' coordinate on [1,3], profile: sine'
      call order_table('production   ', 2, 1, 0, e_inf, e_one,             &
                       rate_prod_sin_v)

      ! ---- 5. the smoothness-factor floor is not the cause ----
      write(*,'(a)') '  DIAGNOSTIC ESWENO floor eps = (kappa dr_j)^2,'//   &
                     ' optimal C,D, uniform grid, sine'
      do i = 1,4
         call order_table_eps(kap(i), rate_eps(i))
         write(*,'(a,es9.2,a,f7.3)') '  DIAGNOSTIC kappa = ', kap(i),      &
                                     '  rate = ', rate_eps(i)
      enddo

      ! ---- 6. the optimal coefficients are exact on a volume quadratic ----
      call volume_quadratic_error(200, exact_err)
      call check_absolute('weno3_optimal_weights_exact_on_volume_quadratic',&
                          exact_err, 0.0d0, 1.0d-12)

      ! ---- verdicts ----
      call check_at_least('weno3_optimal_weights_third_order_uniform',     &
                          rate_fix_sin_u, 2.8d0)
      call check_at_least('weno3_optimal_weights_third_order_mixed_grid',  &
                          rate_fix_sin_m, 2.8d0)
      call check_at_least('weno3_optimal_weights_third_order_gaussian',    &
                          rate_fix_gau_u, 2.8d0)
      ! The proposed fix is the two stencil coefficients alone; these two
      ! assert that it is sufficient, on both grids, with the production D.
      call check_at_least('weno3_corrected_stencil_coefficients_uniform',  &
                          rate_Conly_sin_u, 2.8d0)
      call check_at_least('weno3_corrected_stencil_coefficients_mixed',    &
                          rate_Conly_sin_m, 2.8d0)
      ! The control: on a grid uniform in the volume coordinate C1 = C2 and
      ! the production scheme has no coefficient mismatch left, so it must
      ! reach its design order untouched.  This is what identifies the
      ! cause rather than merely exhibiting a cure.
      call check_at_least('weno3_production_third_order_on_uniform_volume',&
                          rate_prod_sin_v, 2.8d0)

      write(*,'(a,f7.3)') '  DIAGNOSTIC production WENO3 rate,'//          &
           ' uniform grid, sine     = ', rate_prod_sin_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC production WENO3 rate,'//          &
           ' uniform grid, Gaussian = ', rate_prod_gau_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC production WENO3 rate,'//          &
           ' Mixed grid, sine       = ', rate_prod_sin_m
      write(*,'(a,f7.3)') '  DIAGNOSTIC PLM rate,'//                       &
           ' uniform grid, sine                = ', rate_plm_sin_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC PLM rate,'//                       &
           ' Mixed grid, sine                  = ', rate_plm_sin_m
      write(*,'(a,f7.3)') '  DIAGNOSTIC linear weights of the code'//      &
           ' (S = 1), uniform, sine  = ', rate_linD_sin_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC optimal D, production C,'//        &
           ' uniform, sine        = ', rate_Donly_sin_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC corrected C, production D,'//      &
           ' uniform, sine      = ', rate_Conly_sin_u
      write(*,'(a,f7.3)') '  DIAGNOSTIC corrected C, production D,'//      &
           ' Mixed grid, sine   = ', rate_Conly_sin_m
      write(*,'(a,f7.3)') '  DIAGNOSTIC production WENO3,'//               &
           ' uniform in volume, sine        = ', rate_prod_sin_v

      ! The PRODUCTION reconstruction is third order on a general grid:
      ! Reconstruct applies C2(j) to dWp in its WL expression and C1(j-1)
      ! to dWm in its WR expression since 2026-09-07 (RED before that fix:
      ! shortfall 0.803 and 0.812, i.e. rates 1.997 and 1.988).
      call check_at_least('weno3_production_third_order',                  &
                          rate_prod_sin_u, 2.8d0)
      call check_at_least('weno3_production_third_order_mixed',            &
                          rate_prod_sin_m, 2.8d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'weno3_reconstruction_order: ',               &
                             assertion_failures, ' assertion(s) failed'
         error stop 1
      endif

      contains

      !--------------!

      subroutine check_at_least(name, measured, floor_value)
      ! A one-sided verdict in the format of assertion_report: what is
      ! asserted is that the rate falls no lower than the floor, so the
      ! quantity compared is the SHORTFALL max(0, floor - rate), whose
      ! reference is zero at zero tolerance.  A NaN rate gives a NaN
      ! shortfall and fails.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, floor_value
      write(*,'(a,a,a,f7.3,a,f7.3)') '  DIAGNOSTIC ', trim(name),          &
           ': rate = ', measured, ', floor = ', floor_value
      call check_absolute(name, max(0.0d0, floor_value - measured),        &
                          0.0d0, 0.0d0)
      end subroutine check_at_least

      !--------------!

      real*8 function profile(x, which)
      ! The smooth test profiles.  Both are analytic and monotone over most
      ! of [1,3]; the Gaussian adds a curvature scale that is not the domain
      ! size, so a rate measured on it is not an artefact of one shape.
      real*8, intent(in) :: x
      integer, intent(in) :: which
      real*8, parameter :: pi = 3.14159265358979323846d0
      select case (which)
      case (1)
         profile = 2.0d0 + sin(0.5d0*pi*(x - 1.0d0))
      case (2)
         profile = 1.0d0 + exp(-((x - 1.4d0)/0.35d0)**2)
      case (3)
         ! Quadratic in the volume coordinate V = r^3/3: the exact
         ! reconstruction reproduces it, whatever the grid.
         profile = 0.7d0 + 0.3d0*(x**3/3.0d0)                              &
                         + 0.11d0*(x**3/3.0d0)**2
      case default
         profile = 0.0d0
      end select
      end function profile

      !--------------!

      subroutine set_uniform_grid(nc)
      ! Faces equally spaced in r on [1,3]; cell centres at their midpoints,
      ! widths from the faces, exactly the identities define_grid stores.
      integer, intent(in) :: nc
      real*8 :: h
      integer :: j
      N = nc
      call allocate_grid()
      h = 2.0d0/dble(N)
      do j = 1-Ng, N+Ng
         r_edg(j) = 1.0d0 + dble(j)*h
         r(j)     = 1.0d0 + (dble(j) - 0.5d0)*h
         dr_j(j)  = h
      enddo
      end subroutine set_uniform_grid

      !--------------!

      subroutine set_volume_uniform_grid(nc)
      ! Faces equally spaced in the volume coordinate V = r^3/3 on
      ! [1,3], so every dV is the same.  This is the grid on which the
      ! production coefficients C1 and C2 coincide, and it is the control
      ! for the diagnosis: if the coefficient mismatch is what costs the
      ! order, the production scheme is third order here.
      integer, intent(in) :: nc
      real*8 :: v0, v1, hv, vv
      integer :: j
      N = nc
      call allocate_grid()
      v0 = 1.0d0
      v1 = 27.0d0
      hv = (v1 - v0)/dble(N)
      do j = 1-Ng, N+Ng
         vv = v0 + dble(j)*hv
         r_edg(j) = vv**(1.0d0/3.0d0)
      enddo
      do j = 2-Ng, N+Ng
         r(j)    = 0.5d0*(r_edg(j-1) + r_edg(j))
         dr_j(j) = r_edg(j) - r_edg(j-1)
      enddo
      dr_j(1-Ng) = dr_j(2-Ng)
      r(1-Ng)    = r_edg(1-Ng) - 0.5d0*dr_j(1-Ng)
      end subroutine set_volume_uniform_grid

      !--------------!

      subroutine set_mixed_grid(nc)
      ! The production grid.  Refining it means refining EVERY cell, so the
      ! uniform base region keeps its physical extent (dr_base*N_low_cells)
      ! while both of its factors follow the resolution; otherwise added
      ! cells would refine the stretched part alone and no rate would be
      ! measurable.
      integer, intent(in) :: nc
      N = nc
      call allocate_grid()
      grid_type   = 'Mixed'
      r_max       = 3.0d0
      N_low_cells = nc/5
      dr_base     = 0.5d0/dble(N_low_cells)
      call define_grid
      end subroutine set_mixed_grid

      !--------------!

      subroutine allocate_grid()
      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      r = 0.0d0; r_edg = 0.0d0; dr_j = 0.0d0
      end subroutine allocate_grid

      !--------------!

      subroutine cell_averages(which, q)
      ! Volume averages of the profile over every cell, by 8-point
      ! Gauss-Legendre quadrature of int f r^2 dr divided by
      ! (r_p^3 - r_m^3)/3.  Ghost cells are filled from the same analytic
      ! profile, so no boundary treatment enters the interior faces.
      integer, intent(in) :: which
      real*8, intent(out) :: q(1-Ng:N+Ng)
      real*8 :: xg(8), wg(8), rm, rp, mid, hw, x, s, v
      integer :: j, k
      xg = (/ -0.9602898564975363d0, -0.7966664774136267d0,               &
              -0.5255324099163290d0, -0.1834346424956498d0,               &
               0.1834346424956498d0,  0.5255324099163290d0,               &
               0.7966664774136267d0,  0.9602898564975363d0 /)
      wg = (/  0.1012285362903763d0,  0.2223810344533745d0,               &
               0.3137066458778873d0,  0.3626837833783620d0,               &
               0.3626837833783620d0,  0.3137066458778873d0,               &
               0.2223810344533745d0,  0.1012285362903763d0 /)
      do j = 1-Ng, N+Ng
         if (j .eq. 1-Ng) then
            rm = r(j) - 0.5d0*dr_j(j)
         else
            rm = r_edg(j-1)
         endif
         rp = rm + dr_j(j)
         mid = 0.5d0*(rm + rp)
         hw  = 0.5d0*(rp - rm)
         s = 0.0d0
         v = 0.0d0
         do k = 1,8
            x = mid + hw*xg(k)
            s = s + wg(k)*profile(x, which)*x*x
            v = v + wg(k)*x*x
         enddo
         q(j) = s/v
      enddo
      end subroutine cell_averages

      !--------------!

      subroutine volume_widths(dV)
      ! The shell volumes weno3_geometry_coefficients builds its weights
      ! from, in the same index range and with the same two copied ends.
      real*8, intent(out) :: dV(1-Ng:N+Ng)
      integer :: j
      do j = 0, N+1
         dV(j) = r_edg(j)**3 - r_edg(j-1)**3
      enddo
      dV(1-Ng) = dV(2-Ng)
      dV(N+Ng) = dV(N+Ng-1)
      end subroutine volume_widths

      !--------------!

      subroutine optimal_coefficients(cLp, cLm, cRp, cRm, D1o, D2o)
      ! The four stencil coefficients and the two ideal weight ratios of the
      ! third-order reconstruction, from the derivation in the header.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: cLp, cLm, cRp, cRm,     &
                                                   D1o, D2o
      real*8, dimension(1-Ng:N+Ng) :: dV
      real*8 :: a, b, c, m1, m2, al, be, det, kp, km, g0, g1
      integer :: j
      call volume_widths(dV)
      cLp = 0.0d0; cLm = 0.0d0; cRp = 0.0d0; cRm = 0.0d0
      D1o = 1.0d0; D2o = 1.0d0
      do j = 1-Ng+1, N+Ng-1
         a = dV(j-1);  b = dV(j);  c = dV(j+1)
         cLp(j) = b/(b + c)
         cLm(j) = b/(a + b)
         cRp(j) = b/(b + c)
         cRm(j) = b/(a + b)
         m1 = -0.5d0*(a + b)
         m2 =  0.5d0*(b + c)
         al = m2*m2 + (c*c - b*b)/12.0d0
         be = m1*m1 + (a*a - b*b)/12.0d0
         det = m2*be - m1*al
         ! Right face: coefficients of dWp and dWm in p(b/2) - Q_j.
         kp = (be*b/2.0d0 - m1*b*b/6.0d0)/det
         km = (al*b/2.0d0 - m2*b*b/6.0d0)/det
         g0 = kp/cLp(j)
         g1 = km/cLm(j)
         D1o(j) = g1/g0
         ! Left face: coefficients of dWp and dWm in p(-b/2) - Q_j.
         kp = (-be*b/2.0d0 - m1*b*b/6.0d0)/det
         km = (-al*b/2.0d0 - m2*b*b/6.0d0)/det
         g0 = kp/(-cRp(j))
         g1 = km/(-cRm(j))
         D2o(j) = g0/g1
      enddo
      end subroutine optimal_coefficients

      !--------------!

      subroutine weno3_face_states_local(q, cLp, cLm, cRp, cRm, D1u, D2u,  &
                                         eps_scale, use_S, qL, qR)
      ! The production ESWENO3 face formula of Reconstruction.f90, with its
      ! coefficients supplied instead of taken from
      ! weno3_geometry_coefficients, so that one coefficient at a time can be
      ! replaced.  With the production coefficients, use_S true and
      ! eps_scale 1 it is the production routine, statement for statement.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: q, cLp, cLm, cRp, cRm,  &
                                                   D1u, D2u
      real*8, intent(in) :: eps_scale
      logical, intent(in) :: use_S
      real*8, dimension(1-Ng:N+Ng), intent(out) :: qL, qR
      real*8, dimension(1-Ng:N+Ng) :: dq
      real*8 :: dqp, dqm, b0, b1, tau, S0, S1, ep
      integer :: j
      dq(1-Ng:N+Ng-1) = q(2-Ng:N+Ng) - q(1-Ng:N+Ng-1)
      dq(N+Ng) = 0.0d0
      do j = 1-Ng, N+Ng
         qL(j) = q(j)
         qR(j) = q(min(j+1, N+Ng))
      enddo
      do j = 0, N+1
         dqp = dq(j)
         dqm = dq(j-1)
         if (use_S) then
            ep  = (eps_scale*dr_j(j))**2
            b0  = dqp*dqp + ep
            b1  = dqm*dqm + ep
            tau = dqp - dqm
            S0  = 1.0d0 + tau*tau/b0
            S1  = 1.0d0 + tau*tau/b1
         else
            S0 = 1.0d0
            S1 = 1.0d0
         endif
         qL(j)   = q(j)                                                    &
            + (S0*cLp(j)*dqp + D1u(j)*S1*cLm(j)*dqm)/(S0 + D1u(j)*S1)
         qR(j-1) = q(j)                                                    &
            - (D2u(j)*S0*cRp(j)*dqp + S1*cRm(j)*dqm)/(D2u(j)*S0 + S1)
      enddo
      end subroutine weno3_face_states_local

      !--------------!

      subroutine production_coefficients(cLp, cLm, cRp, cRm, D1p, D2p)
      ! The coefficients Reconstruct actually applies, read back from the
      ! production weno3_geometry_coefficients and put in the index
      ! convention weno3_face_states_local expects: cell j's own volume
      ! share on each jump at either face, C2(j) on the jump to j+1 and
      ! C1(j-1) on the jump to j-1 (the mapping Reconstruct applies since
      ! 2026-09-07; before that fix it applied C1(j) and C2(j-1), the
      ! neighbour's shares, and this probe measured the second order that
      ! produced).
      real*8, dimension(1-Ng:N+Ng), intent(out) :: cLp, cLm, cRp, cRm,     &
                                                   D1p, D2p
      real*8, dimension(1-Ng:N+Ng) :: C1, C2, D1, D2
      integer :: j
      call weno3_geometry_coefficients(C1, C2, D1, D2)
      cLp = 0.0d0; cLm = 0.0d0; cRp = 0.0d0; cRm = 0.0d0
      D1p = 1.0d0; D2p = 1.0d0
      do j = 1-Ng+1, N+Ng-1
         cLp(j) = C2(j)
         cLm(j) = C1(j-1)
         cRp(j) = C2(j)
         cRm(j) = C1(j-1)
         D1p(j) = D1(j)
         D2p(j) = D2(j)
      enddo
      end subroutine production_coefficients

      !--------------!

      subroutine face_errors(which, variant, eps_scale, einf, eL1)
      ! L_inf and L_1 of the face error over the interior faces 1..N-1,
      ! where both adjoining cells carry a complete three-cell stencil of
      ! interior data.  Both face states of a face are counted: the left
      ! state comes from cell j and the right state from cell j+1, and a
      ! third-order scheme has to place both on the exact value.
      !
      ! variant  0 production coefficients (identical to Reconstruct_scalar)
      !          1 optimal D, production C
      !          2 corrected stencil coefficients, production D
      !          3 corrected stencil coefficients and optimal D
      !          4 production coefficients with the smoothness factors off
      !         -1 PLM
      integer, intent(in) :: which, variant
      real*8, intent(in)  :: eps_scale
      real*8, intent(out) :: einf, eL1
      real*8, dimension(1-Ng:N+Ng) :: q, qL, qR
      real*8, dimension(1-Ng:N+Ng) :: cLp, cLm, cRp, cRm, D1u, D2u
      real*8, dimension(1-Ng:N+Ng) :: oLp, oLm, oRp, oRm, D1o, D2o
      real*8 :: ex, e
      integer :: j, nf
      call cell_averages(which, q)
      if (variant .eq. -1) then
         rec_method = 'PLM'
         call PLM_rec_scalar(q, qL, qR)
      else
         rec_method = 'WENO3'
         call production_coefficients(cLp, cLm, cRp, cRm, D1u, D2u)
         call optimal_coefficients(oLp, oLm, oRp, oRm, D1o, D2o)
         select case (variant)
         case (1)
            D1u = D1o;  D2u = D2o
         case (2)
            cLp = oLp;  cLm = oLm;  cRp = oRp;  cRm = oRm
         case (3)
            cLp = oLp;  cLm = oLm;  cRp = oRp;  cRm = oRm
            D1u = D1o;  D2u = D2o
         end select
         call weno3_face_states_local(q, cLp, cLm, cRp, cRm, D1u, D2u,     &
                                      eps_scale, variant .ne. 4, qL, qR)
      endif
      einf = 0.0d0
      eL1  = 0.0d0
      nf   = 0
      do j = 1, N-1
         ex = profile(r_edg(j), which)
         e  = abs(qL(j) - ex)
         einf = max(einf, e);  eL1 = eL1 + e;  nf = nf + 1
         e  = abs(qR(j) - ex)
         einf = max(einf, e);  eL1 = eL1 + e;  nf = nf + 1
      enddo
      eL1 = eL1/dble(nf)
      end subroutine face_errors

      !--------------!

      subroutine order_table(label, grid_kind, which, variant, einf, eL1,  &
                             rate_out)
      ! One row per resolution and the L_inf rate between the two finest.
      character(len=*), intent(in) :: label
      integer, intent(in) :: grid_kind, which, variant
      real*8, intent(out) :: einf(n_res), eL1(n_res), rate_out
      real*8 :: rinf, r1
      integer :: i
      write(*,'(a,a)') '  DIAGNOSTIC   scheme = ', trim(label)
      write(*,'(a)') '  DIAGNOSTIC      N      L_inf       rate '//        &
                     '      L_1        rate'
      do i = 1, n_res
         select case (grid_kind)
         case (0)
            call set_uniform_grid(nn(i))
         case (1)
            call set_mixed_grid(nn(i))
         case (2)
            call set_volume_uniform_grid(nn(i))
         end select
         call face_errors(which, variant, 1.0d0, einf(i), eL1(i))
         if (i .eq. 1) then
            write(*,'(a,i6,es12.4,a,es12.4,a)') '  DIAGNOSTIC ', nn(i),    &
                 einf(i), '        -  ', eL1(i), '        -'
         else
            rinf = log(einf(i-1)/einf(i))/log(2.0d0)
            r1   = log(eL1(i-1)/eL1(i))/log(2.0d0)
            write(*,'(a,i6,es12.4,f10.3,es12.4,f10.3)') '  DIAGNOSTIC ',   &
                 nn(i), einf(i), rinf, eL1(i), r1
         endif
      enddo
      rate_out = log(einf(n_res-1)/einf(n_res))/log(2.0d0)
      end subroutine order_table

      !--------------!

      subroutine order_table_eps(kappa, rate_out)
      ! The rate with the optimal coefficients and the ESWENO floor scaled by
      ! kappa: eps = (kappa dr_j)^2 in place of dr_j^2.  If the floor were
      ! what costs the order, the rate would depend on kappa.
      real*8, intent(in)  :: kappa
      real*8, intent(out) :: rate_out
      real*8, dimension(1-Ng:N+Ng) :: q, qL, qR
      real*8, dimension(1-Ng:N+Ng) :: oLp, oLm, oRp, oRm, D1o, D2o
      real*8 :: ei(2), ex, e
      integer :: i, j
      do i = 1, 2
         call set_uniform_grid(nn(n_res-2+i))
         rec_method = 'WENO3'
         call cell_averages(1, q)
         call optimal_coefficients(oLp, oLm, oRp, oRm, D1o, D2o)
         call weno3_face_states_local(q, oLp, oLm, oRp, oRm, D1o, D2o,     &
                                      kappa, .true., qL, qR)
         ei(i) = 0.0d0
         do j = 1, N-1
            ex = profile(r_edg(j), 1)
            e = abs(qL(j) - ex);  ei(i) = max(ei(i), e)
            e = abs(qR(j) - ex);  ei(i) = max(ei(i), e)
         enddo
      enddo
      rate_out = log(ei(1)/ei(2))/log(2.0d0)
      end subroutine order_table_eps

      !--------------!

      subroutine local_copy_deviation(nc, which, dev)
      ! The local formula on the production coefficients against the
      ! production Reconstruct_scalar, on the Mixed grid so that every
      ! coefficient differs from every other.
      integer, intent(in) :: nc, which
      real*8, intent(out) :: dev
      real*8, dimension(1-Ng:N+Ng) :: q, qL, qR, pL, pR
      real*8, dimension(1-Ng:N+Ng) :: cLp, cLm, cRp, cRm, D1p, D2p
      real*8 :: sc
      integer :: j
      call set_mixed_grid(nc)
      rec_method = 'WENO3'
      call cell_averages(which, q)
      call Reconstruct_scalar(q, pL, pR)
      call production_coefficients(cLp, cLm, cRp, cRm, D1p, D2p)
      call weno3_face_states_local(q, cLp, cLm, cRp, cRm, D1p, D2p,        &
                                   1.0d0, .true., qL, qR)
      sc  = maxval(abs(q(1:N)))
      dev = 0.0d0
      do j = 1, N-1
         dev = max(dev, abs(qL(j) - pL(j))/sc)
         dev = max(dev, abs(qR(j) - pR(j))/sc)
      enddo
      end subroutine local_copy_deviation

      !--------------!

      subroutine vector_scalar_deviation(nc, dev)
      ! Reconstruct (the primitive vector) and Reconstruct_scalar have to
      ! place the same face value on the same profile: the density component
      ! of a state at rest is reconstructed by the same formula.  Interior
      ! faces only, since Rec_BC rewrites the boundary ones.
      integer, intent(in) :: nc
      real*8, intent(out) :: dev
      real*8, dimension(1-Ng:N+Ng) :: q, sL, sR
      real*8, dimension(3,1-Ng:N+Ng) :: u, WL, WR
      real*8 :: Wc(3), Uc(3), sc
      integer :: j
      call set_mixed_grid(nc)
      rec_method = 'WENO3'
      use_weno3  = .true.
      call cell_averages(1, q)
      do j = 1-Ng, N+Ng
         Wc(1) = q(j)
         Wc(2) = 0.0d0
         Wc(3) = 1.0d0
         call W_to_U_comp(Wc, Uc, j)
         u(:,j) = Uc
      enddo
      call Reconstruct(u, WL, WR)
      call Reconstruct_scalar(q, sL, sR)
      sc  = maxval(abs(q(1:N)))
      dev = 0.0d0
      do j = 1, N-1
         dev = max(dev, abs(WL(1,j) - sL(j))/sc)
         dev = max(dev, abs(WR(1,j) - sR(j))/sc)
      enddo
      end subroutine vector_scalar_deviation

      !--------------!

      subroutine volume_quadratic_error(nc, err)
      ! A profile that is a quadratic in V = r^3/3 is reproduced exactly by
      ! the third-order reconstruction, on any grid: this is the direct
      ! check of the derived coefficients, independent of any rate.
      integer, intent(in) :: nc
      real*8, intent(out) :: err
      real*8, dimension(1-Ng:N+Ng) :: q, qL, qR
      real*8, dimension(1-Ng:N+Ng) :: oLp, oLm, oRp, oRm, D1o, D2o
      real*8 :: ex, sc
      integer :: j
      call set_mixed_grid(nc)
      rec_method = 'WENO3'
      call cell_averages(3, q)
      call optimal_coefficients(oLp, oLm, oRp, oRm, D1o, D2o)
      ! The smoothness factors are switched off here: they are the nonlinear
      ! part, and what is asserted is that the LINEAR scheme reproduces a
      ! quadratic.
      call weno3_face_states_local(q, oLp, oLm, oRp, oRm, D1o, D2o,        &
                                   1.0d0, .false., qL, qR)
      sc  = maxval(abs(q(1:N)))
      err = 0.0d0
      do j = 1, N-1
         ex  = profile(r_edg(j), 3)
         err = max(err, abs(qL(j) - ex)/sc)
         err = max(err, abs(qR(j) - ex)/sc)
      enddo
      end subroutine volume_quadratic_error

      end program weno3_reconstruction_order

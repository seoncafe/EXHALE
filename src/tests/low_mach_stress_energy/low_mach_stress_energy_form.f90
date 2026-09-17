      program low_mach_stress_energy_form
      ! THE DISCRETE KINETIC-ENERGY QUADRATIC FORM OF THE GATED
      ! FOURTH-DIFFERENCE MOMENTUM STRESS, on the grid and with the
      ! coefficients the code actually uses.
      !
      ! The stress of src/modules/flux/low_mach_dissipation.f90 adds, at the
      ! interior face j between cells j and j+1,
      !
      !   D_j = eps4 g_j rho_f,j lambda_f,j
      !         ( v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1} )              (1)
      !
      ! to the numerical momentum flux and v_f,j D_j to the total-energy
      ! flux.  RK_rhs divides the flux difference by the spherical volume,
      !
      !   d(rho v)_j/dt = - ( A_j D_j - A_{j-1} D_{j-1} ) / dV_j ,
      !   A_j = r_edg(j)^2 ,  dV_j = ( r_edg(j)^3 - r_edg(j-1)^3 ) / 3 ,
      !
      ! and the mass flux is untouched, so at frozen density the discrete
      ! kinetic energy changes at
      !
      !   E' = sum_j dV_j v_j d(rho v)_j/dt
      !      = - sum_j v_j ( A_j D_j - A_{j-1} D_{j-1} )
      !      =   sum_{faces} A_j D_j ( v_{j+1} - v_j ) ,               (2)
      !
      ! exactly, with no leftover boundary term because D vanishes at the
      ! faces j = 0 and j = N.  Writing d_j = v_{j+1} - v_j and w_j =
      ! A_j eps4 g_j rho_f,j lambda_f,j >= 0, the third difference in (1) is
      ! the plain second difference of d, so
      !
      !   E' = sum_j w_j d_j ( d_{j+1} - 2 d_j + d_{j-1} ) .           (3)
      !
      ! For constant w this is -w sum (d_{j+1} - d_j)^2 <= 0.  For varying w
      ! it is d^T W L d with L the second difference, and the symmetric part
      ! of W L is not negative semidefinite in general: the sign claim of the
      ! module header, which integrates by parts as though the coefficient
      ! were constant, is what this program measures.
      !
      ! The correctly signed alternative carries the coefficient BETWEEN the
      ! two differences instead of outside both.  With the cell-centered
      ! second difference q_j = 2 v_j - v_{j+1} - v_{j-1} and a cell-centered
      ! weight k_j >= 0 that vanishes in the first and last cell,
      !
      !   A_j D_j = k_j q_j - k_{j+1} q_{j+1} ,                        (4)
      !
      ! substitution into (2) and one summation by parts give
      !
      !   E' = - sum_j k_j q_j^2 <= 0                                  (5)
      !
      ! for every v and every nonuniform grid, with no boundary term and no
      ! ghost value entering.  In matrix form (4) is dv/dt = -M^-1 B^T K B v
      ! with M = diag(rho_j dV_j), B the cell-centered second difference and
      ! K = diag(k_j), so that d(v^T M v/2)/dt = -(Bv)^T K (Bv).  For
      ! constant k, (4) reduces algebraically to (1): the two forms differ
      ! only through the variation of the coefficient.
      !
      ! WHAT THIS PROGRAM PRINTS.  For each configuration it assembles the
      ! matrix A of E' = -v^T A v over the N physical cells, symmetrizes it,
      ! and calls LAPACK dsyev.  Nonnegative dissipation is lambda_min >= 0.
      ! It also evaluates E' directly on random fields and on the alternating
      ! base mode v_j = (-1)^j a in cells 1-6.
      !
      ! Verdict lines are "PASS|FAIL <name> measured=<v> reference=<r>
      ! tol=<t>"; lines beginning with DIAGNOSTIC or two spaces are context.

      implicit none

      integer, parameter :: dp = kind(1.0d0)
      real(dp), parameter :: mH_cgs = 1.6735575d-24   ! g, m(H)
      real(dp), parameter :: gamma_atomic = 5.0d0/3.0d0

      ! Defaults of parameters.f90 for the damping option.
      real(dp) :: eps4  = 0.02d0
      real(dp) :: m_th  = 1.0d-3

      integer :: nfail
      character(len=512) :: state_file
      integer :: nargs

      nfail = 0
      call get_command_argument(1, state_file)
      nargs = command_argument_count()

      ! 1. Uniform grid, constant coefficients: the case in which the
      !    header's integration by parts is exact.
      call uniform_grid_case(64,  nfail)
      call uniform_grid_case(500, nfail)

      ! 2. A uniform grid with a varying coefficient only: isolates the
      !    coefficient variation from the grid.
      call varying_coefficient_case(500, nfail)

      ! 3. The gate's own step in the coefficient.
      call gate_edge_case(200, nfail)

      ! 4. The actual grid and state, when one is given.
      if (nargs .ge. 1) call state_case(trim(state_file), nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'FAILURES: ', nfail, ' assertion(s) failed'
         error stop 1
      endif

      contains

      ! ------------------------------------------------------!

      subroutine uniform_grid_case(n, nfail)
      ! Uniform dr, constant rho, constant sound speed, gate wide open.  The
      ! area weight A_j = r_edg(j)^2 is set to 1 here, so that w is strictly
      ! constant and (3) is the textbook constant-coefficient form.
      integer, intent(in)    :: n
      integer, intent(inout) :: nfail
      real(dp), allocatable :: w(:), acell(:), kcell(:)
      real(dp) :: lmin_jst, lmax_jst, lmin_bkb, lmax_bkb
      integer  :: j
      character(len=32) :: tag

      allocate(w(1:n-1), acell(0:n), kcell(1:n))
      w     = 1.0d0
      acell = 1.0d0
      do j = 1, n
         kcell(j) = 1.0d0
      enddo
      kcell(1) = 0.0d0;  kcell(n) = 0.0d0

      write(tag,'(A,I0)') 'uniform_n', n
      call spectrum_jst(n, w, lmin_jst, lmax_jst)
      call spectrum_bkb(n, kcell, lmin_bkb, lmax_bkb)

      call verdict('jst_'//trim(tag), lmin_jst, lmax_jst, nfail)
      call verdict('bkb_'//trim(tag), lmin_bkb, lmax_bkb, nfail)

      deallocate(w, acell, kcell)
      end subroutine uniform_grid_case

      ! ------------------------------------------------------!

      subroutine varying_coefficient_case(n, nfail)
      ! Uniform grid, but the face weight falls by four decades over the
      ! domain the way rho does in a stratified base, smoothly.  The only
      ! difference from uniform_grid_case is the coefficient variation.
      integer, intent(in)    :: n
      integer, intent(inout) :: nfail
      real(dp), allocatable :: w(:), kcell(:)
      real(dp) :: lmin_jst, lmax_jst, lmin_bkb, lmax_bkb, x
      integer  :: j

      allocate(w(1:n-1), kcell(1:n))
      do j = 1, n-1
         x    = dble(j-1)/dble(n-2)
         w(j) = 10.0d0**(-4.0d0*x)
      enddo
      do j = 1, n
         x        = dble(j-1)/dble(n-1)
         kcell(j) = 10.0d0**(-4.0d0*x)
      enddo
      kcell(1) = 0.0d0;  kcell(n) = 0.0d0

      call spectrum_jst(n, w, lmin_jst, lmax_jst)
      call spectrum_bkb(n, kcell, lmin_bkb, lmax_bkb)

      call verdict_indefinite('jst_smooth_coefficient_indefinite',        &
                              lmin_jst, lmax_jst, nfail)
      call verdict('bkb_smooth_coefficient', lmin_bkb, lmax_bkb, nfail)

      deallocate(w, kcell)
      end subroutine varying_coefficient_case

      ! ------------------------------------------------------!

      subroutine gate_edge_case(n, nfail)
      ! Uniform grid, uniform coefficient, and the gate (5) closed over the
      ! outer half of the domain: the coefficient steps from w to 0 at one
      ! face, which is what the gate does wherever the Mach number crosses
      ! M_th.  This is the sharpest coefficient variation the production
      ! term can produce, and it is produced by the term's own gate.
      integer, intent(in)    :: n
      integer, intent(inout) :: nfail
      real(dp), allocatable :: w(:), kcell(:), vmin(:)
      real(dp) :: lmin_jst, lmax_jst, lmin_bkb, lmax_bkb, e_direct, e_mat
      real(dp), allocatable :: amat(:,:)
      integer  :: j, lo, hi

      allocate(w(1:n-1), kcell(1:n), vmin(1:n), amat(n,n))
      do j = 1, n-1
         if (j .le. n/2) then
            w(j) = 1.0d0
         else
            w(j) = 0.0d0
         endif
      enddo
      do j = 1, n
         if (j .le. n/2) then
            kcell(j) = 1.0d0
         else
            kcell(j) = 0.0d0
         endif
      enddo
      kcell(1) = 0.0d0;  kcell(n) = 0.0d0

      call spectrum_jst(n, w, lmin_jst, lmax_jst, vmin)
      call spectrum_bkb(n, kcell, lmin_bkb, lmax_bkb)
      call verdict_indefinite('jst_gate_edge_indefinite', lmin_jst,      &
                              lmax_jst, nfail)
      call verdict('bkb_gate_edge', lmin_bkb, lmax_bkb, nfail)

      ! The negative mode, evaluated once through the assembled matrix and
      ! once face by face from the stress itself.
      call assemble_jst(n, w, amat)
      e_mat    = -quad(n, amat, vmin)
      e_direct = energy_rate_jst(n, w, vmin)
      call support_window(n, vmin, lo, hi, 0.99d0)
      write(*,'(A,ES13.6,A,ES13.6)')                                     &
         '  negative mode: E_kin dot (matrix) =', e_mat,                 &
         '  (face sum) =', e_direct
      write(*,'(A,I0,A,I0,A,I0)')                                        &
         '  negative mode lives in cells ', lo, ' to ', hi,              &
         '; the gate closes at face ', n/2

      deallocate(w, kcell, vmin, amat)
      end subroutine gate_edge_case

      ! ------------------------------------------------------!

      subroutine state_case(fname, nfail)
      ! The actual grid and state of an EXHALE Hydro_ioniz.txt: r [R_p],
      ! rho [m_H/cm^3], v [cm/s], p [cgs], T [K], with two ghost cells at
      ! each end.  Face radii are reconstructed as the midpoints of the
      ! written cell centers, since r_edg is not in the output file; the
      ! nonuniformity of the grid is therefore the state's own.
      ! The sound speed uses gamma = 5/3 (the run is atomic).
      ! This configuration is a DIAGNOSTIC, not a verdict: what it reports is
      ! the sign of the form the production coefficients give.
      character(len=*), intent(in) :: fname
      integer, intent(inout)       :: nfail
      real(dp), allocatable :: rc(:), rho(:), vel(:), pre(:)
      real(dp), allocatable :: redg(:), acell(:), w(:), kcell(:), mach2(:)
      real(dp), allocatable :: amat(:,:), vtest(:), vmin(:)
      real(dp) :: lmin_jst, lmax_jst, lmin_bkb, lmax_bkb
      real(dp) :: lmin_int, lmax_int, cs2, gate, e_ran, e_alt, e_alt_bkb
      real(dp) :: r0_cm, gmax, gmin, e_mat, e_direct
      integer  :: n, j, ng, iseed, nopen, lo, hi

      call read_state(fname, n, ng, r0_cm, rc, rho, vel, pre)
      if (n .le. 0) then
         write(*,'(A)') 'FAIL state_read measured=0 reference=N>0 tol=0'
         nfail = nfail + 1
         return
      endif
      write(*,'(A,I0,A,ES12.5,A)') 'DIAGNOSTIC state cells=', n,         &
         ' R0[cm]=', r0_cm, '  file='//trim(fname)

      allocate(redg(0:n), acell(0:n+1), w(1:n-1), kcell(1:n), mach2(0:n+1))
      allocate(vmin(1:n))
      ! Cell centers carry indices 0..n+1 (one ghost each side is enough for
      ! the four-cell stencil of the interior faces).
      do j = 0, n
         redg(j) = 0.5d0*(rc(j) + rc(j+1))*r0_cm
      enddo
      do j = 0, n+1
         acell(j) = (rc(j)*r0_cm)**2
         cs2      = gamma_atomic*pre(j)/(rho(j)*mH_cgs)
         mach2(j) = vel(j)*vel(j)/cs2
      enddo

      ! Face weights of (1), exactly as contact_mode_dissipation_flux forms
      ! them, times the area A_j of the divergence.
      nopen = 0
      do j = 1, n-1
         gate = 1.0d0 - 0.5d0*(mach2(j) + mach2(j+1))/(m_th*m_th)
         if (gate .le. 0.0d0) then
            w(j) = 0.0d0
            cycle
         endif
         gate = gate*gate
         nopen = nopen + 1
         w(j) = redg(j)**2 * eps4 * gate                                 &
              * 0.5d0*(rho(j) + rho(j+1))*mH_cgs                         &
              * 0.5d0*( abs(vel(j))   + sqrt(gamma_atomic*pre(j)         &
                                             /(rho(j)*mH_cgs))           &
                      + abs(vel(j+1)) + sqrt(gamma_atomic*pre(j+1)       &
                                             /(rho(j+1)*mH_cgs)) )
      enddo
      write(*,'(A,I0,A,I0,A)') 'DIAGNOSTIC gate open on ', nopen,        &
         ' of ', n-1, ' interior faces'
      gmax = maxval(w);  gmin = huge(1.0d0)
      do j = 1, n-1
         if (w(j) .gt. 0.0d0 .and. w(j) .lt. gmin) gmin = w(j)
      enddo
      write(*,'(A,ES12.5,A,ES12.5)') 'DIAGNOSTIC face weight w: max=',   &
         gmax, '  min over open faces=', gmin

      ! Cell weights of (4): the same coefficient, cell-centered, zero in
      ! the first and last cell so that (5) has no boundary term.
      do j = 1, n
         gate = 1.0d0 - mach2(j)/(m_th*m_th)
         if (gate .le. 0.0d0) then
            kcell(j) = 0.0d0
            cycle
         endif
         kcell(j) = 0.5d0*(acell(j) + acell(j)) * eps4 * gate*gate       &
                  * rho(j)*mH_cgs                                        &
                  * ( abs(vel(j)) + sqrt(gamma_atomic*pre(j)             &
                                         /(rho(j)*mH_cgs)) )
      enddo
      kcell(1) = 0.0d0;  kcell(n) = 0.0d0

      call spectrum_jst(n, w, lmin_jst, lmax_jst, vmin)
      call spectrum_bkb(n, kcell, lmin_bkb, lmax_bkb)
      call spectrum_jst_interior(n, w, lmin_int, lmax_int)

      write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)')                            &
         'DIAGNOSTIC jst_state          lambda_min=', lmin_jst,          &
         ' lambda_max=', lmax_jst, ' ratio=',                            &
         lmin_jst/max(abs(lmax_jst),tiny(1.0d0))
      write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)')                            &
         'DIAGNOSTIC jst_state_interior lambda_min=', lmin_int,          &
         ' lambda_max=', lmax_int, ' ratio=',                            &
         lmin_int/max(abs(lmax_int),tiny(1.0d0))
      write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)')                            &
         'DIAGNOSTIC bkb_state          lambda_min=', lmin_bkb,          &
         ' lambda_max=', lmax_bkb, ' ratio=',                            &
         lmin_bkb/max(abs(lmax_bkb),tiny(1.0d0))

      ! The negative mode: where it lives, and the same number evaluated
      ! once through the assembled matrix and once face by face.
      allocate(amat(1:n,1:n), vtest(1:n))
      call assemble_jst(n, w, amat)
      e_mat    = -quad(n, amat, vmin)
      e_direct = energy_rate_jst(n, w, vmin)
      call support_window(n, vmin, lo, hi, 0.99d0)
      write(*,'(A,ES13.6,A,ES13.6)')                                     &
         'DIAGNOSTIC jst_state negative mode: E_kin dot (matrix) =',     &
         e_mat, '  (face sum) =', e_direct
      write(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)')                           &
         'DIAGNOSTIC jst_state negative mode lives in cells ', lo,       &
         ' to ', hi, ', r = ', rc(lo), ' to ', rc(hi)

      ! The energy form on explicit fields.
      do iseed = 1, 5
         call random_field(n, iseed, vtest)
         e_ran = -quad(n, amat, vtest)
         write(*,'(A,I0,A,ES13.6)') 'DIAGNOSTIC jst_state random field ',&
            iseed, ' E_kin dot =', e_ran
      enddo
      vtest = 0.0d0
      do j = 1, 6
         vtest(j) = dble((-1)**j)
      enddo
      e_alt = -quad(n, amat, vtest)
      call assemble_bkb(n, kcell, amat)
      e_alt_bkb = -quad(n, amat, vtest)
      write(*,'(A,ES13.6)')                                              &
         'DIAGNOSTIC jst_state alternating base mode (cells 1-6) E_kin dot =',&
         e_alt
      write(*,'(A,ES13.6)')                                              &
         'DIAGNOSTIC bkb_state alternating base mode (cells 1-6) E_kin dot =',&
         e_alt_bkb

      ! The verdicts for this configuration: the existing stress is NOT
      ! negative semidefinite here, and the corrected form is.
      call verdict_indefinite('jst_state_indefinite', lmin_jst, lmax_jst,&
                              nfail)
      call verdict('bkb_state', lmin_bkb, lmax_bkb, nfail)

      deallocate(rc, rho, vel, pre, redg, acell, w, kcell, mach2,        &
                 amat, vtest, vmin)
      end subroutine state_case

      ! ------------------------------------------------------!

      subroutine verdict_indefinite(name, lmin, lmax, nfail)
      ! The assertion that the form is NOT negative semidefinite, i.e. that
      ! there is a velocity field on which the stress puts kinetic energy IN.
      ! The reference is the accuracy of the symmetric eigenvalue problem
      ! itself, 100 eps lambda_max: a negative eigenvalue below it would be
      ! the rounding of dsyev and not a property of the operator.
      character(len=*), intent(in) :: name
      real(dp), intent(in)         :: lmin, lmax
      integer, intent(inout)       :: nfail
      real(dp) :: anchor
      anchor = -100.0d0*epsilon(1.0d0)*abs(lmax)
      if (lmin .le. anchor) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'PASS ', name,             &
            ' measured=', lmin, ' reference=', anchor, ' tol=0'
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'FAIL ', name,             &
            ' measured=', lmin, ' reference=', anchor, ' tol=0'
         nfail = nfail + 1
      endif
      write(*,'(A,ES13.6,A,ES13.6)') '  lambda_min/lambda_max=',         &
         lmin/max(abs(lmax),tiny(1.0d0)), '  lambda_max=', lmax
      end subroutine verdict_indefinite

      ! ------------------------------------------------------!

      subroutine verdict(name, lmin, lmax, nfail)
      ! Nonnegative dissipation.  The reference is the same eigenvalue-solver
      ! accuracy, 100 eps lambda_max: below it a negative eigenvalue is
      ! rounding and not the operator.
      character(len=*), intent(in) :: name
      real(dp), intent(in)         :: lmin, lmax
      integer, intent(inout)       :: nfail
      real(dp) :: anchor
      anchor = -100.0d0*epsilon(1.0d0)*abs(lmax)
      if (lmin .ge. anchor) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'PASS ', name,             &
            ' measured=', lmin, ' reference=', anchor, ' tol=0'
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A)') 'FAIL ', name,             &
            ' measured=', lmin, ' reference=', anchor, ' tol=0'
         nfail = nfail + 1
      endif
      write(*,'(A,ES13.6,A,ES13.6)') '  lambda_min/lambda_max=',         &
         lmin/max(abs(lmax),tiny(1.0d0)), '  lambda_max=', lmax
      end subroutine verdict

      ! ------------------------------------------------------!

      subroutine assemble_jst(n, w, amat)
      ! A of E' = -v^T A v for the stress (1).  The ghost values v_0 and
      ! v_{n+1} the four-cell stencil of the first and last interior face
      ! reaches are closed by zero gradient (v_0 = v_1, v_{n+1} = v_n), the
      ! closure a reflecting/outflow ghost gives; the closure-independent
      ! statement is the interior submatrix below.
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: w(1:n-1)
      real(dp), intent(out) :: amat(1:n,1:n)
      real(dp) :: a(-1:2), b(-1:2)
      integer  :: j, ka, kb, ia, ib

      amat = 0.0d0
      do j = 1, n-1
         if (w(j) .eq. 0.0d0) cycle
         ! d_j = v_{j+1} - v_j
         a(-1) = 0.0d0;  a(0) = -1.0d0;  a(1) = 1.0d0;  a(2) = 0.0d0
         ! third difference v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1}
         b(-1) = -1.0d0; b(0) = 3.0d0;  b(1) = -3.0d0; b(2) = 1.0d0
         do ka = -1, 2
            if (a(ka) .eq. 0.0d0) cycle
            ia = closure(j+ka, n)
            do kb = -1, 2
               if (b(kb) .eq. 0.0d0) cycle
               ib = closure(j+kb, n)
               amat(ia,ib) = amat(ia,ib) - w(j)*a(ka)*b(kb)
            enddo
         enddo
      enddo
      end subroutine assemble_jst

      ! ------------------------------------------------------!

      subroutine assemble_jst_interior(n, w, amat, lo, hi)
      ! The same operator restricted to test fields supported in cells
      ! lo..hi, so that no ghost value enters and the sign statement is
      ! independent of the boundary closure.
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: w(1:n-1)
      real(dp), intent(out) :: amat(1:n,1:n)
      integer,  intent(in)  :: lo, hi
      real(dp) :: a(-1:2), b(-1:2)
      integer  :: j, ka, kb, ia, ib

      amat = 0.0d0
      do j = 1, n-1
         if (w(j) .eq. 0.0d0) cycle
         a(-1) = 0.0d0;  a(0) = -1.0d0;  a(1) = 1.0d0;  a(2) = 0.0d0
         b(-1) = -1.0d0; b(0) = 3.0d0;  b(1) = -3.0d0; b(2) = 1.0d0
         do ka = -1, 2
            if (a(ka) .eq. 0.0d0) cycle
            ia = j+ka
            if (ia .lt. lo .or. ia .gt. hi) cycle
            do kb = -1, 2
               if (b(kb) .eq. 0.0d0) cycle
               ib = j+kb
               if (ib .lt. lo .or. ib .gt. hi) cycle
               amat(ia,ib) = amat(ia,ib) - w(j)*a(ka)*b(kb)
            enddo
         enddo
      enddo
      end subroutine assemble_jst_interior

      ! ------------------------------------------------------!

      subroutine assemble_bkb(n, kcell, amat)
      ! A of E' = -v^T A v for the flux (4).  No ghost value enters, because
      ! k_1 = k_n = 0 removes q_1 and q_n from the sum.  The area A_j of the
      ! divergence cancels against the A_j of (4) and so does not appear:
      ! k_j carries the whole coefficient.
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: kcell(1:n)
      real(dp), intent(out) :: amat(1:n,1:n)
      real(dp) :: a(-1:2), b(-1:2)
      integer  :: j, ka, kb, ia, ib

      amat = 0.0d0
      do j = 1, n-1
         ! A_j D_j = k_j q_j - k_{j+1} q_{j+1}
         !         = k_j  (2 v_j - v_{j+1} - v_{j-1})
         !         - k_{j+1}(2 v_{j+1} - v_{j+2} - v_j)
         b(-1) = -kcell(j)
         b(0)  =  2.0d0*kcell(j) + kcell(j+1)
         b(1)  = -kcell(j) - 2.0d0*kcell(j+1)
         b(2)  =  kcell(j+1)
         if (all(b .eq. 0.0d0)) cycle
         a(-1) = 0.0d0;  a(0) = -1.0d0;  a(1) = 1.0d0;  a(2) = 0.0d0
         do ka = -1, 2
            if (a(ka) .eq. 0.0d0) cycle
            ia = closure(j+ka, n)
            do kb = -1, 2
               if (b(kb) .eq. 0.0d0) cycle
               ib = closure(j+kb, n)
               amat(ia,ib) = amat(ia,ib) - a(ka)*b(kb)
            enddo
         enddo
      enddo
      end subroutine assemble_bkb

      ! ------------------------------------------------------!

      integer function closure(idx, n)
      ! Zero-gradient ghost closure.
      integer, intent(in) :: idx, n
      closure = idx
      if (idx .lt. 1) closure = 1
      if (idx .gt. n) closure = n
      end function closure

      ! ------------------------------------------------------!

      subroutine spectrum_jst(n, w, lmin, lmax, vmin)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: w(1:n-1)
      real(dp), intent(out) :: lmin, lmax
      real(dp), intent(out), optional :: vmin(n)
      real(dp), allocatable :: amat(:,:)
      allocate(amat(n,n))
      call assemble_jst(n, w, amat)
      if (present(vmin)) then
         call sym_spectrum(n, amat, lmin, lmax, vmin)
      else
         call sym_spectrum(n, amat, lmin, lmax)
      endif
      deallocate(amat)
      end subroutine spectrum_jst

      ! ------------------------------------------------------!

      subroutine spectrum_jst_interior(n, w, lmin, lmax)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: w(1:n-1)
      real(dp), intent(out) :: lmin, lmax
      real(dp), allocatable :: amat(:,:), sub(:,:)
      integer :: lo, hi, m, i, j
      lo = 3;  hi = n-2;  m = hi - lo + 1
      allocate(amat(n,n), sub(m,m))
      call assemble_jst_interior(n, w, amat, lo, hi)
      do j = 1, m
         do i = 1, m
            sub(i,j) = amat(lo+i-1, lo+j-1)
         enddo
      enddo
      call sym_spectrum(m, sub, lmin, lmax)
      ! lo, hi are 3 and n-2: no ghost value enters this block.
      deallocate(amat, sub)
      end subroutine spectrum_jst_interior

      ! ------------------------------------------------------!

      subroutine spectrum_bkb(n, kcell, lmin, lmax)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: kcell(1:n)
      real(dp), intent(out) :: lmin, lmax
      real(dp), allocatable :: amat(:,:)
      allocate(amat(n,n))
      call assemble_bkb(n, kcell, amat)
      call sym_spectrum(n, amat, lmin, lmax)
      deallocate(amat)
      end subroutine spectrum_bkb

      ! ------------------------------------------------------!

      subroutine sym_spectrum(n, amat, lmin, lmax, vmin)
      ! Eigenvalues of the symmetric part (A + A^T)/2, which is the only
      ! part a quadratic form sees.  vmin, if asked for, returns the
      ! eigenvector of the smallest eigenvalue.
      integer,  intent(in)    :: n
      real(dp), intent(in)    :: amat(n,n)
      real(dp), intent(out)   :: lmin, lmax
      real(dp), intent(out), optional :: vmin(n)
      real(dp), allocatable :: s(:,:), evals(:), work(:)
      integer :: i, j, lwork, info
      character(len=1) :: jobz
      allocate(s(n,n), evals(n))
      do j = 1, n
         do i = 1, n
            s(i,j) = 0.5d0*(amat(i,j) + amat(j,i))
         enddo
      enddo
      lwork = 8*n + 64
      allocate(work(lwork))
      jobz = 'N'
      if (present(vmin)) jobz = 'V'
      call dsyev(jobz, 'U', n, s, n, evals, work, lwork, info)
      if (info .ne. 0) then
         write(*,'(A,I0)') 'FAIL dsyev measured=info reference=0 tol=0 info=', info
         lmin = 0.0d0;  lmax = 0.0d0
         if (present(vmin)) vmin = 0.0d0
      else
         lmin = evals(1);  lmax = evals(n)
         if (present(vmin)) vmin(1:n) = s(1:n,1)
      endif
      deallocate(s, evals, work)
      end subroutine sym_spectrum

      ! ------------------------------------------------------!

      real(dp) function quad(n, amat, v)
      integer,  intent(in) :: n
      real(dp), intent(in) :: amat(n,n), v(n)
      integer :: i, j
      quad = 0.0d0
      do j = 1, n
         do i = 1, n
            quad = quad + v(i)*amat(i,j)*v(j)
         enddo
      enddo
      end function quad

      ! ------------------------------------------------------!

      real(dp) function energy_rate_jst(n, w, v)
      ! E' of equation (2) evaluated face by face from the stress (1), with
      ! the zero-gradient ghost closure.  This is an independent evaluation
      ! of the same number the assembled matrix gives, so the two together
      ! check the assembly.
      integer,  intent(in) :: n
      real(dp), intent(in) :: w(1:n-1), v(1:n)
      integer :: j
      real(dp) :: vm1, vp2
      energy_rate_jst = 0.0d0
      do j = 1, n-1
         if (w(j) .eq. 0.0d0) cycle
         vm1 = v(closure(j-1,n))
         vp2 = v(closure(j+2,n))
         energy_rate_jst = energy_rate_jst                               &
            + w(j)*(v(j+1) - v(j))                                       &
                  *(vp2 - 3.0d0*v(j+1) + 3.0d0*v(j) - vm1)
      enddo
      end function energy_rate_jst

      ! ------------------------------------------------------!

      subroutine support_window(n, v, lo, hi, frac)
      ! The shortest window of cells carrying the fraction frac of the
      ! squared norm of v, reported so that a negative mode can be located.
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: v(n)
      integer,  intent(out) :: lo, hi
      real(dp), intent(in)  :: frac
      real(dp) :: tot, run
      integer  :: j
      tot = sum(v*v)
      run = 0.0d0
      lo = 1;  hi = n
      do j = 1, n
         run = run + v(j)*v(j)
         if (run .ge. 0.5d0*(1.0d0-frac)*tot) then
            lo = j
            exit
         endif
      enddo
      run = 0.0d0
      do j = n, 1, -1
         run = run + v(j)*v(j)
         if (run .ge. 0.5d0*(1.0d0-frac)*tot) then
            hi = j
            exit
         endif
      enddo
      end subroutine support_window

      ! ------------------------------------------------------!

      subroutine random_field(n, iseed, v)
      ! A reproducible pseudo-random field: a linear congruential sequence,
      ! so that the numbers do not depend on the compiler's random_number.
      integer,  intent(in)  :: n, iseed
      real(dp), intent(out) :: v(n)
      integer(kind=8) :: s
      integer :: j
      s = int(iseed, 8)*2654435761_8 + 12345_8
      do j = 1, n
         s = mod(6364136223846793005_8*s + 1442695040888963407_8,        &
                 4611686018427387904_8)
         v(j) = 2.0d0*(dble(abs(s))/4.611686018427387904d18) - 1.0d0
      enddo
      end subroutine random_field

      ! ------------------------------------------------------!

      subroutine read_state(fname, n, ng, r0_cm, rc, rho, vel, pre)
      ! Reads r [R_p], rho [m_H/cm^3], v [cm/s], p [cgs] of an EXHALE
      ! Hydro_ioniz.txt.  The "# grid" header line states N and R0[cm] and
      ! the "# rows" line states how many ghost cells precede the physical
      ! ones.  One ghost on each side is kept, which is the stencil the
      ! interior faces of the stress reach.
      character(len=*), intent(in) :: fname
      integer,  intent(out) :: n, ng
      real(dp), intent(out) :: r0_cm
      real(dp), allocatable, intent(out) :: rc(:), rho(:), vel(:), pre(:)
      character(len=8192) :: line
      integer :: u, ios, i, k, nrow
      real(dp) :: a, b, c, d

      n = 0;  ng = 2;  r0_cm = 0.0d0
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios .ne. 0) return

      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .ne. '#') then
            backspace(u)
            exit
         endif
         k = index(line, '# grid N ')
         if (k .gt. 0) read(line(k+9:),*,iostat=ios) n
         k = index(line, 'R0[cm] ')
         if (k .gt. 0) read(line(k+7:),*,iostat=ios) r0_cm
         k = index(line, '# rows ')
         if (k .gt. 0) then
            i = index(line, ':')
            if (i .gt. 0) read(line(i+1:),*,iostat=ios) ng
         endif
      enddo
      if (n .le. 0 .or. r0_cm .le. 0.0d0) then
         close(u)
         n = 0
         return
      endif

      nrow = n + 2*ng
      allocate(rc(0:n+1), rho(0:n+1), vel(0:n+1), pre(0:n+1))
      do i = 1, nrow
         read(u,*,iostat=ios) a, b, c, d
         if (ios .ne. 0) then
            close(u)
            n = 0
            return
         endif
         ! Rows ng+1 .. ng+n are the physical cells; keep one ghost each side.
         k = i - ng
         if (k .ge. 0 .and. k .le. n+1) then
            rc(k)  = a;  rho(k) = b;  vel(k) = c;  pre(k) = d
         endif
      enddo
      close(u)
      end subroutine read_state

      end program low_mach_stress_energy_form

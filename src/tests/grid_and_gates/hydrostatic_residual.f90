      program hydrostatic_residual
      ! THE OPERATOR IDENTITY UNDER TEST
      !   On a state that is hydrostatic (v = 0, dp/dr = -rho dphi/dr), the
      !   momentum row of the right-hand side must return zero:
      !
      !     R_mom(j) = dF(2,j) - S(2,j)
      !              = [ r^2 F(2) ]_{j+1/2} - [ r^2 F(2) ]_{j-1/2}  over dV
      !                (+ the WENO3 face-pressure term (p_R - p_L)/dr_j)
      !                + 0.5(rho_L + rho_R)(Gphi_i(j) - Gphi_i(j-1))/dr_j(j)
      !                (- the PLM geometric pressure term, which lives in S)
      !
      !   The marching stage is u = u - dt(dF - S) (EXHALE_main), so R_mom is
      !   exactly the momentum the operator adds to a cell that should stay at
      !   rest, and dv = -R_mom dt/rho is the velocity one CFL step produces.
      !
      ! WHAT IS MEASURED, AND ON WHAT
      !   The analytic isothermal hydrostatic column of the mechanical
      !   `hydrostatic_column` case (backup/regression/hydrostatic_column),
      !   evaluated as CELL VOLUME AVERAGES on the run's own grid, at
      !   Grid cells = 250, 500, 1000, 2000, for both reconstructions (PLM,
      !   ESWENO3) and both Riemann solvers (ROE, HLLC).  This is the first
      !   act of increment B4-5 of docs/b4_spatial_operator_design_20260906.md
      !   section 2.6: the marched case has no discrete hydrostatic state to
      !   measure a residual against (its README records that the column
      !   drains), while the analytic column needs no converged run and
      !   answers the question the refinement ladder presupposes -- is the
      !   pressure/gravity pair itself the source of the flow, or is it the
      !   boundary treatment.
      !
      !   Everything the residual is built from is PRODUCTION: define_grid,
      !   set_gravity_grid, Apply_BC (the characteristic base condition and
      !   the free-outflow top), Reconstruct (PLM / ESWENO3 + Rec_BC + the
      !   face positivity limiter), Num_flux (ROE / HLLC) and RK_rhs with
      !   Source.  Nothing here re-implements a flux, a source or a grid; the
      !   driver builds the state, calls the right-hand side ONCE, and reads
      !   dF and S back out.
      !
      ! THE COLUMN, in code units (lengths R0, velocities v0 = sqrt(k T0/m_H),
      ! densities n0 m_H, pressures n0 k T0, potential phi = -b0/r):
      !
      !     T(r)   = 1                     (isothermal at T0)
      !     p(r)   = nhat rho(r)           (p = n_part k T, n_part = nhat rho)
      !     rho(r) = exp[ -(phi(r) - phi(1))/nhat ]
      !     v(r)   = 0
      !
      !   with nhat the particle count per unit mass of the neutral He/H
      !   mixture of the case.  Substituting into dp/dr = -rho dphi/dr with
      !   the code's own phi reproduces the exponential exactly, so the state
      !   is hydrostatic in the CONTINUUM; whether the discrete operator sees
      !   it as such is the measurement.
      !
      ! NORMALIZATION
      !   R_mom is reported divided by the local weight rho(j) Dphi(r(j)),
      !   with Dphi the production gravitational gradient: the term the
      !   pressure divergence has to cancel.  1 means the operator leaves the
      !   whole weight of the cell unbalanced.
      !
      ! ONE GRID PER PROCESS, AND WHY
      !   `define_grid` gives the Mixed grid's stretch factor a Newton
      !   iteration whose convergence flag, `real*8 :: tol = 1.0`, carries an
      !   initializer and is therefore SAVEd: on any call after the first the
      !   loop is skipped and the stretch parameter keeps its initial guess
      !   1.01, so a second grid built in the same process is not the grid its
      !   N asks for.  Production calls `define_grid` exactly once (init.f90),
      !   so nothing in a run is affected, but a refinement ladder cannot be
      !   walked inside one process.  This driver therefore measures ONE grid
      !   per invocation and appends its numbers to `hydrostatic_ladder.dat`
      !   in the working directory; invoked with no argument it reads that
      !   file back, prints the ladder with its observed orders, and gives the
      !   verdicts.
      !
      ! ASSERTIONS (only statements that are identities today)
      !   1. flat_state_zero_residual: with b0 = 0 and a uniform state the
      !      pressure divergence and the geometric source are the same number
      !      formed two ways, so R_mom is zero to round-off.  An identity of
      !      the discretization at any grid and for either scheme.
      !   2. interior_residual_falls_with_N: the largest normalized residual
      !      over the interior cells (the physical cells whose stencil reads
      !      no ghost) is smaller at N = 2000 than at N = 250.  A consistent
      !      operator's truncation error goes to zero with dr; the ORDER at
      !      which it does is reported, not asserted.
      !   3. outer_cell_residual_is_consistent: the observed order of the
      !      OUTERMOST physical cell's residual in 1/N is at least 1.  A
      !      boundary closure that is consistent with the equations leaves a
      !      truncation error there like any other cell; one that states
      !      something the interior does not satisfy leaves a residual that
      !      does not converge at all.  The bound is the lowest order a
      !      consistent closure can have and is not a tolerance: the
      !      zero-gradient outflow ghost this replaced measured 0.006 (PLM),
      !      the isothermal hydrostatic continuation measures 1.80.
      !   The base cell is reported and not asserted: which closure it carries
      !   is the characteristic condition of base_boundary, a different object
      !   from the interior pair and from the outer condition.

      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use grav_func,                  only: phi, Dphi
      use base_boundary,              only: set_base_reservoir
      use BC_Apply,                   only: Apply_BC
      use Reconstruction_step,        only: Reconstruct
      use RK_integration,             only: RK_rhs
      use eval_time_step,             only: eval_dt
      use Conversion,                 only: U_to_W

      implicit none

      ! ---- the mechanical column of backup/regression/hydrostatic_column,
      !      READ from its input.inp ----
      real*8, parameter :: Rp_RJ   = 0.49d0     ! Planet radius [R_J]
      real*8, parameter :: Mp_MJ   = 0.0457d0   ! Planet mass [M_J]
      real*8, parameter :: Teq     = 1140.0d0   ! Equilibrium temperature [K]
      real*8, parameter :: rmax_c  = 3.0d0      ! Outer radius [R_p]
      real*8, parameter :: resc_c  = 2.0d0      ! Escape radius [R_p]
      real*8, parameter :: he_h    = 0.0793d0   ! He/H number ratio
      ! He mass in units of the hydrogen atom mass (species_table.f90).
      real*8, parameter :: m_He    = 3.9715d0

      integer, parameter :: n_scheme = 4
      character(len=12), dimension(n_scheme) :: scheme_name
      character(len=*), parameter :: ladder_file = 'hydrostatic_ladder.dat'

      ! Particles per unit mass of the neutral atomic He/H mixture, in units
      ! of one particle per hydrogen atom mass: (1 + He/H)/(1 + m_He He/H).
      ! p = n_part k T then reads p = nhat rho T in code units, so the
      ! isothermal sound speed squared of the column is nhat.
      real*8 :: nhat

      character(len=32) :: arg
      integer :: Nrun, is, nfail, ios, nlow_arg
      real*8  :: drb_arg

      scheme_name(1) = 'PLM_ROE'
      scheme_name(2) = 'PLM_HLLC'
      scheme_name(3) = 'WENO3_ROE'
      scheme_name(4) = 'WENO3_HLLC'

      nhat  = (1.0d0 + he_h)/(1.0d0 + m_He*he_h)
      nfail = 0

      ! Fixed configuration of the case, independent of the grid.
      spherical_domain = .true.        ! "Domain mode: Spherical"
      dr_base          = 2.0e-4        ! "Base grid [dr,cells]:" defaults,
      N_low_cells      = 50            !   the same default-real literals
      r_max            = rmax_c        !   input_read uses
      r_esc            = resc_c
      r_flux           = 1.20d0
      grid_type        = 'Mixed'
      CFL              = 0.6d0

      ! Usage:  hydrostatic_residual [N [dr_base N_low]]
      ! With N: measure that grid and append its records. The optional pair
      ! states "Base grid [dr,cells]:"; without it the case's own defaults
      ! stand, which hold the base region at 50 cells of 2.0e-4 R_p while
      ! N refines the stretched region alone. With no argument at all: read
      ! the records back, print the ladder and give the verdicts.
      call get_command_argument(1, arg)

      if (len_trim(arg) .gt. 0) then

         read(arg,*,iostat=ios) Nrun
         if (ios .ne. 0 .or. Nrun .lt. 100) then
            write(*,'(A)') 'FAIL hydrostatic_residual_argument '//        &
                 'measured=bad reference=integer_N tol=0'
            call exit(1)
         endif
         call get_command_argument(2, arg)
         if (len_trim(arg) .gt. 0) then
            read(arg,*,iostat=ios) drb_arg
            if (ios .eq. 0) dr_base = drb_arg
            call get_command_argument(3, arg)
            if (len_trim(arg) .gt. 0) then
               read(arg,*,iostat=ios) nlow_arg
               if (ios .eq. 0) N_low_cells = nlow_arg
            endif
         endif
         call measure_one_grid(Nrun)

      else

         call report_ladder(nfail)
         if (nfail .gt. 0) then
            write(*,'(A,I0,A)') 'hydrostatic_residual: ', nfail,          &
                                ' assertion(s) FAILED'
            call exit(1)
         endif
         write(*,'(A)') 'hydrostatic_residual: all assertions PASSED'

      endif

      contains

      ! ------------------------------------------------------------------ !

      subroutine set_scheme(is)
      ! The flags the operator dispatches on. use_plm puts the pressure
      ! inside the momentum flux (Phys_flux) and the geometric term in the
      ! source (Source); use_weno3 leaves the pressure out of the flux and
      ! adds the face-pressure difference in RK_rhs. They are the pair
      ! input_read sets from "Reconstruction scheme:" and "Numerical flux:".
      integer, intent(in) :: is
      select case (is)
      case (1)
         rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
         flux = 'ROE'
      case (2)
         rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
         flux = 'HLLC'
      case (3)
         rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
         flux = 'ROE'
      case (4)
         rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
         flux = 'HLLC'
      end select
      end subroutine set_scheme

      ! ------------------------------------------------------------------ !

      real*8 function cell_average(f, ja)
      ! Volume average of a radial function over cell ja,
      !   <f> = int f(r) r^2 dr / int r^2 dr  over [r_edg(ja-1), r_edg(ja)],
      ! by the 8-point Gauss-Legendre rule, exact for polynomials of degree
      ! 15. Its error on the exponential column is far below the
      ! reconstruction's own, so the initial data is the exact discrete
      ! representation of the atmosphere and what is measured is the operator
      ! and not the sampling.
      interface
         real*8 function f(x)
         real*8, intent(in) :: x
         end function f
      end interface
      integer, intent(in) :: ja
      real*8 :: a, b, xm, xh, x, num, den
      integer :: k
      ! Nodes and weights on [-1,1] (Abramowitz & Stegun table 25.4).
      real*8, dimension(4), parameter :: gl_x = (/                        &
         0.1834346424956498d0, 0.5255324099163290d0,                      &
         0.7966664774136267d0, 0.9602898564975363d0 /)
      real*8, dimension(4), parameter :: gl_w = (/                        &
         0.3626837833783620d0, 0.3137066458778873d0,                      &
         0.2223810344533745d0, 0.1012285362903763d0 /)

      a  = r_edg(ja-1)
      b  = r_edg(ja)
      xm = 0.5d0*(a + b)
      xh = 0.5d0*(b - a)
      num = 0.0d0
      den = 0.0d0
      do k = 1,4
         x   = xm + xh*gl_x(k)
         num = num + gl_w(k)*f(x)*x*x
         den = den + gl_w(k)*x*x
         x   = xm - xh*gl_x(k)
         num = num + gl_w(k)*f(x)*x*x
         den = den + gl_w(k)*x*x
      enddo
      cell_average = num/den

      end function cell_average

      ! ------------------------------------------------------------------ !

      real*8 function rho_column(x)
      ! The analytic isothermal hydrostatic density, with the code's own phi.
      real*8, intent(in) :: x
      rho_column = exp(-(phi(x) - phi(1.0d0))/nhat)
      end function rho_column

      real*8 function p_column(x)
      ! p = n_part k T = nhat rho T with T = T0 (code T = 1).
      real*8, intent(in) :: x
      p_column = nhat*rho_column(x)
      end function p_column

      ! ------------------------------------------------------------------ !

      subroutine measure_one_grid(Ncells)
      ! One point of the ladder. The grid is built once (see the header note
      ! on define_grid), then the flat-state identity and the analytic column
      ! are measured on it for each scheme, and one record per scheme is
      ! appended to the ladder file.
      integer, intent(in) :: Ncells

      real*8, dimension(:,:), allocatable :: u, W, WL, WR, dF, S
      real*8, dimension(:),   allocatable :: dt_loc, rho_a, p_a, Rn
      real*8  :: b0_case, dt, cs, dv, wgt, acc, sc
      real*8  :: flat_err, rmx, rrms, rint, dvcs, rmx_r, rint_r, dvcs_r
      real*8  :: Rbase1, Rbase2, Rtop2, Rtop1, Rmid, drmid
      integer :: j, jmx, jint, jmid, jdv, iu

      ! --- the grid, through the production routines, with the
      !     normalizations of input_read ---
      N  = Ncells
      T0 = Teq
      R0 = Rp_RJ*RJ
      Mp = Mp_MJ*MJ
      v0 = sqrt(kb_erg*T0/mu)
      b0_case = (Gc*Mp*mu)/(kb_erg*T0*R0)

      call allocate_grid_arrays
      call define_grid

      allocate(u(3,1-Ng:N+Ng), W(3,1-Ng:N+Ng), WL(3,1-Ng:N+Ng),           &
               WR(3,1-Ng:N+Ng), dF(3,1-Ng:N+Ng), S(3,1-Ng:N+Ng),          &
               dt_loc(1-Ng:N+Ng), rho_a(1-Ng:N+Ng), p_a(1-Ng:N+Ng),       &
               Rn(1-Ng:N+Ng))

      jmid = 1
      do j = 1,N
         if (r(j) .ge. 2.0d0) then
            jmid = j
            exit
         endif
      enddo
      drmid = dr_j(jmid)

      write(*,'(A,I5,A,ES12.5,A,ES12.5,A,F10.5)')                         &
         '  N=', N, '  dr_j(1)=', dr_j(1), '  dr_j(r=2)=', drmid,         &
         '  b0=', b0_case

      open(newunit=iu, file=ladder_file, position='append',               &
           status='unknown')

      do is = 1,n_scheme

         call set_scheme(is)

         ! ---- the flat state: uniform, no gravity ----
         b0 = 0.0d0
         call set_gravity_grid
         do j = 1-Ng,N+Ng
            u(1,j) = 1.0d0
            u(2,j) = 0.0d0
            u(3,j) = 1.0d0/(gamma_ad - 1.0d0)
         enddo
         n_part_cell1 = 1.0d0
         call set_base_reservoir(1.0d0, 1.0d0, 1.0d0, 1.0d0)

         call Apply_BC(u)
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)

         ! Normalized by the size of the two pressure terms themselves,
         ! (r_+^2 + r_-^2) p/dV, so the verdict is scale free.
         flat_err = 0.0d0
         do j = 1,N
            sc = (r_edg(j)*r_edg(j) + r_edg(j-1)*r_edg(j-1))              &
                 /((r_edg(j)**3 - r_edg(j-1)**3)/3.0d0)
            flat_err = max(flat_err, abs(dF(2,j) - S(2,j))/sc)
         enddo

         ! ---- the analytic hydrostatic column ----
         b0 = b0_case
         call set_gravity_grid

         do j = 1-Ng,N+Ng
            rho_a(j) = cell_average(rho_column, j)
            p_a(j)   = cell_average(p_column,   j)
            u(1,j)   = rho_a(j)
            u(2,j)   = 0.0d0
            u(3,j)   = p_a(j)/(gamma_ad - 1.0d0)
         enddo

         ! What the lower boundary states: the column's own (p, s) at r = 1,
         ! and the particle count of cell 1 the characteristic condition
         ! needs in order to read T(1) = p(1)/n_part(1).
         call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
         n_part_cell1 = nhat*rho_a(1)

         call Apply_BC(u)
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)

         call U_to_W(u, W)
         call eval_dt(W, dt, dt_loc)

         rmx  = 0.0d0
         rint = 0.0d0
         jmx  = 1
         jint = 3
         jdv  = 1
         acc  = 0.0d0
         dvcs = 0.0d0
         do j = 1,N
            wgt   = rho_a(j)*Dphi(r(j))
            Rn(j) = (dF(2,j) - S(2,j))/wgt
            acc   = acc + Rn(j)*Rn(j)
            if (abs(Rn(j)) .gt. rmx) then
               rmx = abs(Rn(j))
               jmx = j
            endif
            ! The interior: the physical cells whose reconstruction stencil
            ! and whose two faces read no ghost state.
            if (j .ge. 3 .and. j .le. N-2) then
               if (abs(Rn(j)) .gt. rint) then
                  rint = abs(Rn(j))
                  jint = j
               endif
            endif
            cs = sqrt(gamma_ad*p_a(j)/rho_a(j))
            dv = abs(dF(2,j) - S(2,j))*dt/rho_a(j)
            if (dv/cs .gt. dvcs) then
               dvcs = dv/cs
               jdv  = j
            endif
         enddo
         rrms   = sqrt(acc/real(N,8))
         rmx_r  = r(jmx)
         rint_r = r(jint)
         dvcs_r = r(jdv)
         Rbase1     = Rn(1)
         Rbase2     = Rn(2)
         Rtop2   = Rn(N-1)
         Rtop1     = Rn(N)
         Rmid      = Rn(jmid)

         write(*,'(A,A10,A,ES11.4,A,ES11.4,A,I6,A,F8.4,A,ES11.4,A,ES11.4)')&
            '     ', trim(scheme_name(is)),                               &
            '  max=', rmx, ' rms=', rrms, ' at j=', jmx, ' r=', rmx_r,    &
            '  interior max=', rint, '  R(r=2)=', Rmid

         write(iu,'(I8,I3,15ES24.16)') N, is, dr_j(1), drmid, flat_err,   &
            rmx, rrms, rint, rmx_r, rint_r, Rbase1, Rbase2, Rtop2, Rtop1,  &
            Rmid, dvcs, dvcs_r

      enddo

      close(iu)
      deallocate(u, W, WL, WR, dF, S, dt_loc, rho_a, p_a, Rn)

      end subroutine measure_one_grid

      ! ------------------------------------------------------------------ !

      subroutine report_ladder(nfail)
      ! Read every record written by the per-grid invocations, print the
      ! ladder with its observed orders, and give the verdicts.
      integer, intent(inout) :: nfail

      integer, parameter :: mx = 64
      integer :: Nv(mx), sv(mx)
      real*8  :: dr1(mx), drm(mx), fe(mx), rmx(mx), rrms(mx), rint(mx)
      real*8  :: rmxr(mx), rintr(mx), Rbase1(mx), Rbase2(mx), Rtop2(mx), Rtop1(mx)
      real*8  :: Rmid(mx), dvcs(mx), dvcsr(mx), hN(mx)
      integer :: nrec, iu, ios, i, k, is, kprev
      real*8  :: ord_all, ord_int, ord_dv, ord_b1, ord_mid, ord_top, femax

      nrec = 0
      open(newunit=iu, file=ladder_file, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') 'FAIL hydrostatic_residual_ladder_file '//        &
              'measured=missing reference='//ladder_file//' tol=0'
         nfail = nfail + 1
         return
      endif
      do
         if (nrec .ge. mx) exit
         read(iu,*,iostat=ios) Nv(nrec+1), sv(nrec+1), dr1(nrec+1),       &
              drm(nrec+1), fe(nrec+1), rmx(nrec+1), rrms(nrec+1),         &
              rint(nrec+1), rmxr(nrec+1), rintr(nrec+1), Rbase1(nrec+1),      &
              Rbase2(nrec+1), Rtop2(nrec+1), Rtop1(nrec+1),           &
              Rmid(nrec+1), dvcs(nrec+1), dvcsr(nrec+1)
         if (ios .ne. 0) exit
         nrec = nrec + 1
      enddo
      close(iu)

      if (nrec .lt. 2*n_scheme) then
         write(*,'(A,I0,A)') 'FAIL hydrostatic_residual_ladder_file '//   &
              'measured=', nrec, ' reference=at_least_two_grids tol=0'
         nfail = nfail + 1
         return
      endif

      do is = 1,n_scheme

         write(*,'(A)') ''
         write(*,'(A,A)') '  DIAGNOSTIC ladder ', trim(scheme_name(is))
         write(*,'(A)') '         N       dr_j(1)     dr_j(r=2)  '//      &
              '  max|R|/rho g     rms        interior max  r(int max)'//  &
              '   R(1)         R(2)        R(r=2)       R(N-1)       '//  &
              'R(N)       max dv/cs   r(dv max)'

         do i = 1,nrec
            if (sv(i) .ne. is) cycle
            write(*,'(A,I8,12ES13.4,F12.4)')                              &
               '     ', Nv(i), dr1(i), drm(i), rmx(i), rrms(i), rint(i),  &
               rintr(i), Rbase1(i), Rbase2(i), Rmid(i), Rtop2(i),         &
               Rtop1(i), dvcs(i), dvcsr(i)
         enddo

         ! Two orders, because the Mixed grid refines its two regions by
         ! different factors and no single length describes both. The
         ! quantities that live where dr shrinks with N are quoted in the
         ! LOCAL width (R at r = 2, in dr_j(r=2)); the ones that are maxima
         ! over the whole domain are quoted in 1/N, the one refinement
         ! parameter every cell shares.
         do i = 1,nrec
            hN(i) = 1.0d0/real(Nv(i),8)
         enddo
         call ladder_order(is, nrec, sv, hN,  rmx,    ord_all)
         call ladder_order(is, nrec, sv, hN,  rint,   ord_int)
         call ladder_order(is, nrec, sv, hN,  dvcs,   ord_dv)
         call ladder_order(is, nrec, sv, dr1, Rbase1, ord_b1)
         call ladder_order(is, nrec, sv, drm, Rmid,   ord_mid)
         write(*,'(A,3(A,F8.4))') '     DIAGNOSTIC order in 1/N       ',  &
            '  max=', ord_all, '  interior max=', ord_int,                &
            '  max dv/cs=', ord_dv
         write(*,'(A,2(A,F8.4))') '     DIAGNOSTIC order in local dr  ',  &
            '  R(r=2) in dr_j(r=2)=', ord_mid,                            &
            '  R(1) in dr_j(1)=', ord_b1

      enddo

      write(*,'(A)') ''
      do is = 1,n_scheme
         femax = 0.0d0
         do i = 1,nrec
            if (sv(i) .eq. is) femax = max(femax, fe(i))
         enddo
         call verdict_below('flat_state_zero_residual['//                 &
                            trim(scheme_name(is))//']', femax, 1.0d-14,   &
                            nfail)
      enddo

      do is = 1,n_scheme
         k     = 0
         kprev = 0
         do i = 1,nrec
            if (sv(i) .ne. is) cycle
            if (kprev .eq. 0) kprev = i
            k = i
         enddo
         call verdict_below('interior_residual_falls_with_N['//           &
                            trim(scheme_name(is))//']',                   &
                            rint(k)/rint(kprev), 1.0d0, nfail)
      enddo

      ! The outermost physical cell, the cell the free-outflow ghost is read
      ! into. Its residual must converge, and at least at first order; the
      ! order is the assertion because the LEVEL of a boundary residual is a
      ! property of the grid and the level alone cannot tell a converging
      ! closure from a stalled one.
      do is = 1,n_scheme
         call ladder_order(is, nrec, sv, hN, Rtop1, ord_top)
         call verdict_above('outer_cell_residual_is_consistent['//        &
                            trim(scheme_name(is))//']', ord_top, 1.0d0,   &
                            nfail)
      enddo

      end subroutine report_ladder

      ! ------------------------------------------------------------------ !

      subroutine ladder_order(is, nrec, sv, h, q, ord)
      ! Least-squares slope of log|q| against log h over the records of one
      ! scheme: the observed order of q in the length h.
      integer, intent(in)  :: is, nrec
      integer, intent(in)  :: sv(*)
      real*8,  intent(in)  :: h(*), q(*)
      real*8,  intent(out) :: ord
      real*8  :: sx, sy, sxx, sxy, x, y
      integer :: i, m
      sx = 0.0d0; sy = 0.0d0; sxx = 0.0d0; sxy = 0.0d0; m = 0
      do i = 1,nrec
         if (sv(i) .ne. is) cycle
         if (.not. (abs(q(i)) .gt. 0.0d0)) cycle
         if (.not. (h(i) .gt. 0.0d0)) cycle
         x = log(h(i))
         y = log(abs(q(i)))
         sx = sx + x;  sy = sy + y
         sxx = sxx + x*x;  sxy = sxy + x*y
         m = m + 1
      enddo
      ! Fewer than two points, or a ladder in which this length does not
      ! move at all (a base refinement holds N fixed, and the reverse), has
      ! no slope to report.
      if (m .lt. 2 .or. .not. (abs(m*sxx - sx*sx) .gt. 0.0d0)) then
         ord = 0.0d0
         return
      endif
      ord = (m*sxy - sx*sy)/(m*sxx - sx*sx)
      end subroutine ladder_order

      ! ------------------------------------------------------------------ !

      subroutine verdict_below(name, measured, bound, nfail)
      ! One-sided verdict: the measured quantity must be strictly below the
      ! stated bound. Written as the negation of "below", so a NaN fails.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, bound
      integer,          intent(inout) :: nfail
      if (measured .lt. bound) then
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'PASS ', name,              &
              ' measured=', measured, ' reference=', bound, ' tol=0'
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'FAIL ', name,              &
              ' measured=', measured, ' reference=', bound, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_below

      ! ------------------------------------------------------------------ !

      subroutine verdict_above(name, measured, bound, nfail)
      ! One-sided verdict: the measured quantity must be at or above the
      ! stated bound. Written as the negation of "at or above", so a NaN
      ! fails.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, bound
      integer,          intent(inout) :: nfail
      if (measured .ge. bound) then
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'PASS ', name,              &
              ' measured=', measured, ' reference=', bound, ' tol=0'
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5,A)') 'FAIL ', name,              &
              ' measured=', measured, ' reference=', bound, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_above

      end program hydrostatic_residual

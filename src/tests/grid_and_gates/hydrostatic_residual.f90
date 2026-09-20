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
      !
      ! THE MOMENTUM ROW'S REFERENCE SCALE UNDER THE WELL-BALANCED KEY
      !   The three assertions of measure_momentum_row_scale below are of the
      !   NORMALIZED row, R_2/s_2 with s_2 = residual_row_scale(2,...), the
      !   one expression the stationary solver's acceptance test reads: on
      !   the scheme's own discrete equilibrium it must be at the rounding
      !   level and not at one; on a perturbed state it must be below one and
      !   proportional to the perturbation; and with no gravity it must be
      !   the scale the base scheme uses.  The key cancels the equilibrium
      !   pressure force and the gravitational source against each other
      !   before the row is formed, so a scale built from the row's remaining
      !   terms alone is the row itself and reads one however small the
      !   imbalance is (Update_EXHALE N37).
      !
      ! THE MOMENTUM ROW'S REFERENCE SCALE AGAINST THE TERMS OF THE EQUATION
      !   The five assertions of measure_momentum_physical_terms below ask
      !   what the scale is on four states whose momentum terms are known in
      !   closed form: no term at all (no gravity, uniform pressure, at
      !   rest), the weight alone (the scheme's own discrete hydrostatic
      !   equilibrium, with the key off and on), an imbalance proportional to a
      !   perturbation, and a pure ram divergence (supersonic uniform flow).
      !   A scale made of the discretization's pieces, max(|dF_2|, |S_2|),
      !   fails the first and the last under PLM, where the pressure is
      !   split between the momentum flux and the geometric source and the
      !   two O(2 p/r) halves cancel between them.

      use global_parameters
      use grid_construction,          only: define_grid,                 &
                                            spherical_cell_volume
      use gravity_grid_construction,  only: set_gravity_grid
      use grav_func,                  only: phi, Dphi
      use base_boundary,              only: set_base_reservoir
      use BC_Apply,                   only: Apply_BC
      use Reconstruction_step,        only: Reconstruct
      use RK_integration,             only: RK_rhs, face_flux,            &
                                            face_q_up, face_q_dn
      use steady_residual_mod,        only: residual_row_scale
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
      character(len=*), parameter :: wb_file = 'hydrostatic_wellbalanced.dat'
      character(len=*), parameter :: sc_file = 'hydrostatic_momentum_scale.dat'
      character(len=*), parameter :: pt_file = 'hydrostatic_momentum_terms.dat'

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
                 /spherical_cell_volume(j)
            flat_err = max(flat_err, abs(dF(2,j) - S(2,j))/sc)
         enddo

         ! ---- the analytic hydrostatic column ----
         b0 = b0_case
         call set_gravity_grid

         ! The lowest ghost has no lower face (r_edg does not reach below
         ! r_edg(1-Ng)), so its volume average is not defined and is not
         ! taken: Apply_BC writes that column, and RK_rhs forms no row for
         ! it.  It carries the neighboring average here so that nothing
         ! reads an undefined double.
         do j = 2-Ng,N+Ng
            rho_a(j) = cell_average(rho_column, j)
            p_a(j)   = cell_average(p_column,   j)
         enddo
         rho_a(1-Ng) = rho_a(2-Ng)
         p_a(1-Ng)   = p_a(2-Ng)
         do j = 1-Ng,N+Ng
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

      ! ---- the well-balanced key, on the same grid ----
      call measure_well_balanced(u, WL, WR, dF, S, rho_a, p_a, Rn, b0_case)
      call measure_momentum_row_scale(u, WL, WR, dF, S, rho_a, b0_case)
      call measure_momentum_physical_terms(u, WL, WR, dF, S, rho_a,       &
                                           b0_case)

      deallocate(u, W, WL, WR, dF, S, dt_loc, rho_a, p_a, Rn)

      end subroutine measure_one_grid

      ! ------------------------------------------------------------------ !

      subroutine measure_well_balanced(u, WL, WR, dF, S, rho_a, p_a, Rn,  &
                                       b0_case)
      ! WHAT IS MEASURED, AND WHY TWO COLUMNS.  The well-balanced key
      ! ("Well balanced:", default off) preserves the equilibrium ITS OWN
      ! discretization defines, which is the state whose two neighboring
      ! equilibrium extrapolations agree at every shared face,
      !
      !   p_j+1 + rho_j+1 (phi_c(j+1) - phi_i(j))
      !         = p_j - rho_j (phi_i(j) - phi_c(j)),
      !
      ! the spherical form of Kaeppeli and Mishra (2016, A&A 587, A94,
      ! their eq. 18).  That is the column measured here as `discrete`, built
      ! by integrating the relation outward from the innermost ghost with the
      ! analytic density averages; on it the momentum, mass and energy rows
      ! must fall to the rounding level at every N and for both
      ! reconstructions.  The ANALYTIC column of the ladder above is measured
      ! as well (`analytic`): it satisfies the CONTINUUM balance and departs
      ! from the discrete one by the truncation error of the equilibrium
      ! extrapolation, so on it the key is second order like the base scheme
      ! and the two columns together say which part of the residual is the
      ! discretization of gravity and which is the floating-point assembly.
      ! Their paper measures the same pair: their initial data is the
      ! DISCRETE equilibrium solved for by a Newton iteration (their eq. 42),
      ! and their fig. 2 is the rounding-level line that follows.
      !
      ! The rows are normalized by the size of the terms they are built from:
      ! the momentum row by max_j rho_j |dphi/dr|, the weight the pressure
      ! gradient has to carry; the mass row by max_j (A+ + A-) rho_j c_j/dV
      ! and the energy row by the same with c_j^3, the size of an interface
      ! flux over a cell.  The window is the interior, cells 3 to N-2, whose
      ! stencils and faces read no ghost: the base face carries the
      ! characteristic boundary condition and the top the free-outflow
      ! extrapolation, neither of which is the interior pair.
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u, WL, WR, dF, S
      real*8, dimension(1-Ng:N+Ng),   intent(inout) :: rho_a, p_a, Rn
      real*8, intent(in) :: b0_case

      real*8, dimension(1-Ng:N+Ng) :: p_d, rho_d
      real*8  :: mom_ref, mass_ref, ene_ref, sc, cs, dVj
      real*8  :: mom_mx, mass_mx, ene_mx
      integer :: j, iw, icol, iwb, iu2
      logical :: sav_wb

      sav_wb = well_balanced

      b0 = b0_case
      call set_gravity_grid

      do j = 2-Ng,N+Ng
         rho_a(j) = cell_average(rho_column, j)
         p_a(j)   = cell_average(p_column,   j)
      enddo
      rho_a(1-Ng) = rho_a(2-Ng)
      p_a(1-Ng)   = p_a(2-Ng)

      ! The discrete equilibrium of the scheme, at the column's own
      ! temperature T = 1: with p = nhat rho T the face-matching condition
      !   p_j+1 + rho_j+1 (phi_c(j+1) - phi_i(j))
      !         = p_j - rho_j (phi_i(j) - phi_c(j))
      ! is a two-term recursion for the density,
      !   rho_j+1 = rho_j (nhat - a_j)/(nhat + b_j),
      ! with a_j and b_j the two half-cell potential differences.  It is
      ! solved outward from the innermost ghost of the analytic column, which
      ! is the isothermal case of the Newton solve Kaeppeli and Mishra use to
      ! initialize their exactness test (2016, A&A 587, A94, their eq. 42).
      ! Integrating the PRESSURE instead with the analytic densities is not
      ! the same state and leaves the outer cells at negative pressure: that
      ! density profile has no positive discrete equilibrium on this grid.
      call discrete_equilibrium_density(rho_a(1), rho_d)
      p_d = nhat*rho_d

      open(newunit=iu2, file=wb_file, position='append', status='unknown')

      do is = 1,n_scheme
      do iwb = 0,1
      do icol = 1,2

         call set_scheme(is)
         well_balanced = (iwb .eq. 1)

         do j = 1-Ng,N+Ng
            u(2,j) = 0.0d0
            if (icol .eq. 1) then
               u(1,j) = rho_a(j)
               u(3,j) = p_a(j)/(gamma_ad - 1.0d0)
            else
               u(1,j) = rho_d(j)
               u(3,j) = p_d(j)/(gamma_ad - 1.0d0)
            endif
         enddo

         call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
         n_part_cell1 = nhat*rho_a(1)

         call Apply_BC(u)
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)

         mom_ref  = 0.0d0
         mass_ref = 0.0d0
         ene_ref  = 0.0d0
         mom_mx   = 0.0d0
         mass_mx  = 0.0d0
         ene_mx   = 0.0d0
         do j = 3,N-2
            dVj = spherical_cell_volume(j)
            sc  = (r_edg(j)*r_edg(j) + r_edg(j-1)*r_edg(j-1))/dVj
            cs  = sqrt(gamma_ad*u(3,j)*(gamma_ad - 1.0d0)/u(1,j))
            mom_ref  = max(mom_ref,  u(1,j)*Dphi(r(j)))
            mass_ref = max(mass_ref, sc*u(1,j)*cs)
            ene_ref  = max(ene_ref,  sc*u(1,j)*cs*cs*cs)
            mom_mx   = max(mom_mx,  abs(dF(2,j) - S(2,j)))
            mass_mx  = max(mass_mx, abs(dF(1,j) - S(1,j)))
            ene_mx   = max(ene_mx,  abs(dF(3,j) - S(3,j)))
         enddo

         write(iu2,'(I8,3I3,3ES24.16)') N, is, iwb, icol,                 &
            mom_mx/mom_ref, mass_mx/mass_ref, ene_mx/ene_ref

      enddo
      enddo
      enddo

      close(iu2)
      well_balanced = sav_wb

      end subroutine measure_well_balanced

      ! ------------------------------------------------------------------ !

      subroutine discrete_equilibrium_density(rho_base, rho_d)
      ! THE DISCRETE EQUILIBRIUM OF THE SCHEME at the column's own
      ! temperature T = 1.  With p = nhat rho T the face-matching condition
      !   p_j+1 + rho_j+1 (phi_c(j+1) - phi_i(j))
      !         = p_j - rho_j (phi_i(j) - phi_c(j))
      ! is a two-term recursion for the density,
      !   rho_j+1 = rho_j (nhat - a_j)/(nhat + b_j),
      ! with a_j and b_j the two half-cell potential differences.  It is
      ! solved outward and inward from the density of cell 1, which is the
      ! isothermal case of the Newton solve Kaeppeli and Mishra use to
      ! initialize their exactness test (2016, A&A 587, A94, their eq. 42).
      ! Integrating the PRESSURE instead with the analytic densities is not
      ! the same state and leaves the outer cells at negative pressure:
      ! that density profile has no positive discrete equilibrium on this
      ! grid.
      real*8, intent(in)  :: rho_base
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_d
      integer :: j
      rho_d(1) = rho_base
      do j = 1,N+Ng-1
         rho_d(j+1) = rho_d(j)                                           &
            *(nhat - (Gphi_i(j)   - Gphi_c(j)))                          &
            /(nhat + (Gphi_c(j+1) - Gphi_i(j)))
      enddo
      do j = 1,2-Ng,-1
         rho_d(j-1) = rho_d(j)                                           &
            *(nhat + (Gphi_c(j)   - Gphi_i(j-1)))                        &
            /(nhat - (Gphi_i(j-1) - Gphi_c(j-1)))
      enddo
      end subroutine discrete_equilibrium_density

      ! ------------------------------------------------------------------ !

      subroutine measure_momentum_row_scale(u, WL, WR, dF, S, rho_a,     &
                                            b0_case)
      ! WHAT IS MEASURED: the momentum row DIVIDED BY ITS OWN REFERENCE
      ! SCALE, R_2(j)/residual_row_scale(2,j,u), the quantity the stationary
      ! solver's certification compares against its tolerance.  Three states
      ! and the largest value over the interior cells 3..N-2 on each:
      !
      !   equilibrium  the discrete equilibrium of the scheme, on which the
      !                key's row is at the rounding level.  The scaled row
      !                must be there too; a scale built from the row's
      !                remaining terms alone returns the row itself and
      !                reads exactly one.
      !   perturbed    the same column with the pressure multiplied by
      !                1 + eps sin(...), at eps and eps/2.  The density and
      !                the potential are untouched, so the reference scale
      !                is the same number on both, and the scaled row must
      !                stay below one and halve with the perturbation: the
      !                equilibrium departure and every face quantity built
      !                from it are homogeneous of degree one in eps (the
      !                limited slope is positively homogeneous and the WENO3
      !                weights are scale free), so the ratio is 2 to the
      !                accuracy of the O(eps^2) terms.
      !   no gravity   b0 = 0 on a state with structure in it, where the
      !                equilibrium pressure force is zero and the key's
      !                scale must therefore be the two terms the row still
      !                holds, max(|ram|, |dp/dr|), exactly.  The departure
      !                of the key's scale from the BASE SCHEME's on the same
      !                state is reported and not asserted: with no gravity
      !                each side's own hydrostatic equilibrium is its own
      !                cell pressure, so the key's pressure gradient is the
      !                base scheme's, term for term, and the two scales are
      !                the same number up to the rounding of the two flux
      !                assemblies (MEASURED below).
      !
      ! The scale is a function of the state, so the SAME state is never
      ! asked for it with the key off and on without a state in between:
      ! refresh_row_terms answers a repeated request for one state from what
      ! it already holds.
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u, WL, WR, dF, S
      real*8, dimension(1-Ng:N+Ng),   intent(inout) :: rho_a
      real*8, intent(in) :: b0_case

      real*8, parameter :: eps_p = 1.0d-6
      real*8, dimension(1-Ng:N+Ng) :: rho_d, p_d, s_on, s_off, s_row
      real*8  :: s_eq, s_p1, s_p2, s_zero, s_own, shape_p, twopi, span
      real*8  :: x_on, s_break
      real*8  :: dAp, dAm, dV, ram, pgr
      integer :: j, iu3, ip
      logical :: sav_wb

      sav_wb = well_balanced
      twopi  = 8.0d0*atan(1.0d0)

      b0 = b0_case
      call set_gravity_grid
      do j = 2-Ng,N+Ng
         rho_a(j) = cell_average(rho_column, j)
      enddo
      rho_a(1-Ng) = rho_a(2-Ng)
      call discrete_equilibrium_density(rho_a(1), rho_d)
      p_d  = nhat*rho_d
      span = r(N) - r(1)

      open(newunit=iu3, file=sc_file, position='append',                 &
           status='unknown')

      do is = 1,n_scheme

         call set_scheme(is)
         well_balanced = .true.

         ! ---- the discrete equilibrium, and the two perturbations of it ----
         s_eq = 0.0d0
         s_p1 = 0.0d0
         s_p2 = 0.0d0
         do ip = 0,2

            do j = 1-Ng,N+Ng
               shape_p = 1.0d0
               if (ip .gt. 0) shape_p = 1.0d0                            &
                  + eps_p/dble(2**(ip-1))                                &
                    *sin(4.0d0*twopi*(r(j) - r(1))/span)
               u(1,j) = rho_d(j)
               u(2,j) = 0.0d0
               u(3,j) = shape_p*p_d(j)/(gamma_ad - 1.0d0)
            enddo

            call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
            n_part_cell1 = nhat*rho_d(1)

            call Apply_BC(u)
            call Reconstruct(u, WL, WR)
            call RK_rhs(u, WL, WR, dF, S)

            ! The row terms are cached against the state they belong to,
            ! and this state repeats across the schemes; asking for
            ! another state's scale first is what makes the request below
            ! form the terms of THIS scheme.
            s_break = residual_row_scale(2, 3, 1.5d0*u)

            x_on = 0.0d0
            do j = 3,N-2
               x_on = max(x_on, abs(dF(2,j) - S(2,j))                    &
                                /residual_row_scale(2, j, u))
            enddo
            if (ip .eq. 0) s_eq = x_on
            if (ip .eq. 1) s_p1 = x_on
            if (ip .eq. 2) s_p2 = x_on

         enddo

         ! ---- no gravity: the key's scale is the base scheme's ----
         b0 = 0.0d0
         call set_gravity_grid
         do j = 1-Ng,N+Ng
            u(1,j) = 1.0d0/r(j)**3
            u(2,j) = u(1,j)*0.3d0*r(j)
            u(3,j) = 0.5d0*u(2,j)*u(2,j)/u(1,j)                          &
                   + (0.8d0/r(j)**4)/(gamma_ad - 1.0d0)
         enddo
         call set_base_reservoir(0.8d0, 1.0d0, nhat, 1.0d0)
         n_part_cell1 = nhat*u(1,1)
         call Apply_BC(u)

         well_balanced = .true.
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)
         ! THE ROW'S OWN TERMS ARE THE TERMS OF THE MOMENTUM EQUATION, the
         ! ram divergence and the pressure gradient, each formed here from
         ! the face data this evaluation stored: |dF_2| is their SUM and is
         ! smaller than either wherever they cancel, which is what the scale
         ! must not be.  The key's pressure gradient is the face pressure
         ! measured from each side's own equilibrium, and with no gravity
         ! that equilibrium is the cell pressure itself, so these are also
         ! the base scheme's two terms.
         do j = 3,N-2
            dAp = r_edg(j)*r_edg(j)
            dAm = r_edg(j-1)*r_edg(j-1)
            dV  = spherical_cell_volume(j)
            ram = (dAp*face_flux(2,j) - dAm*face_flux(2,j-1))/dV
            if (use_plm) then
               pgr = (dAp*face_q_up(j) - dAm*face_q_dn(j-1))/dV
            else
               pgr = (face_q_up(j) - face_q_dn(j-1))/dr_j(j)
            endif
            s_row(j) = max(abs(ram), abs(pgr))
         enddo
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         do j = 3,N-2
            s_on(j) = residual_row_scale(2, j, u)
         enddo
         ! Another state in between, so that the scale of u is formed
         ! again with the key the other way and not answered from what this one
         ! left.
         well_balanced = .false.
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         do j = 3,N-2
            s_off(j) = residual_row_scale(2, j, u)
         enddo

         s_zero = 0.0d0
         s_own  = 0.0d0
         do j = 3,N-2
            if (s_off(j) .gt. 0.0d0)                                     &
               s_zero = max(s_zero, abs(s_on(j) - s_off(j))/s_off(j))
            if (s_row(j) .gt. 0.0d0)                                     &
               s_own = max(s_own, abs(s_on(j) - s_row(j))/s_row(j))
         enddo

         well_balanced = .true.
         b0 = b0_case
         call set_gravity_grid

         write(iu3,'(I8,I3,5ES24.16)') N, is, s_eq, s_p1, s_p2, s_zero,  &
            s_own

      enddo

      close(iu3)
      well_balanced = sav_wb

      end subroutine measure_momentum_row_scale

      ! ------------------------------------------------------------------ !

      subroutine measure_momentum_physical_terms(u, WL, WR, dF, S,       &
                                                 rho_a, b0_case)
      ! WHAT IS MEASURED: the momentum row's REFERENCE SCALE,
      ! residual_row_scale(2,j,u), against the terms of the momentum
      ! equation
      !
      !   d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc ,
      !
      ! on four states whose terms are known in closed form.  The scale is
      ! what the stationary solver's certification divides the row by, so
      ! what it measures decides which states are accepted; a scale built
      ! from the DISCRETIZATION's pieces, max(|dF_2|, |S_2|), is not those
      ! terms, because under PLM the pressure sits partly in the momentum
      ! flux (Phys_flux adds p) and partly in the source (the geometric term
      ! (A+ - A-) p_c/dV) with the two O(2 p/r) parts cancelling between
      ! them.
      !
      !   zero gravity, uniform pressure, at rest.  Every term of the
      !        equation is zero, so the scale must be zero to the rounding
      !        of its own assembly.  It is reported as a FRACTION OF THE
      !        GEOMETRIC PRESSURE TERM (A+ - A-) p_j/dV, which is what a
      !        scale made of the discretization's pieces returns here and
      !        is 2 p/r as dr goes to zero.  A fully zero row is divided by
      !        the tiny floor of momentum_row_scale, the only place the
      !        denominator is not a term of the row; the two schemes reach
      !        that state at the rounding of two different assemblies, eps p
      !        over dr under WENO3 against eps p over r under PLM, so the
      !        statement is made against the physical term and not between
      !        the two schemes, whose measured values are reported.
      !   the discrete hydrostatic equilibrium of the scheme, at rest.  The
      !        pressure gradient and the weight are equal and opposite and
      !        are the only terms, so the scale is the weight, and the
      !        weight of the scheme the row was assembled by: the half-sum of
      !        the two face densities times the interface potential
      !        difference for the base scheme, as Source.f90 forms it, and
      !        the pressure force of the cell's own equilibrium under
      !        "Well balanced:", where that cancellation is algebraic.  The
      !        two differ by the O(dr/r) of their geometric weighting and
      !        are not interchangeable references.  Under the key that
      !        weight is the scale exactly; for the base scheme what is
      !        asserted is that the scale is not below it.
      !   the same equilibrium with the pressure perturbed by eps and
      !        eps/2.  The scale is the same number on both (the density and
      !        the potential are untouched) and the imbalance is
      !        homogeneous of degree one in eps, so the scaled row halves.
      !   supersonic uniform flow, no gravity, uniform pressure, v = 3 c_s.
      !        The ram divergence is the only term: (A+ - A-) rho v v/dV
      !        exactly, the faces of a uniform state being the cell values.
      !        A scale that keeps the pressure inside the momentum flux
      !        reads (A+ - A-)(rho v v + p)/dV instead, high by 1/(gamma M M).
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u, WL, WR, dF, S
      real*8, dimension(1-Ng:N+Ng),   intent(inout) :: rho_a
      real*8, intent(in) :: b0_case

      real*8, parameter :: eps_p = 1.0d-6
      real*8, parameter :: mach  = 3.0d0
      real*8, dimension(1-Ng:N+Ng) :: rho_d, p_d, R0
      real*8  :: z_frac, z_scale, e_base, b_base, e_wb, p_rat, u_ram
      real*8  :: x1, x2, sc, wgt, epf, ramr, cs_u, v_u
      real*8  :: shape_p, twopi, span, geo, dAp, dAm, dV
      real*8  :: s_break
      integer :: j, iu4, ip
      logical :: sav_wb

      sav_wb = well_balanced
      twopi  = 8.0d0*atan(1.0d0)

      b0 = b0_case
      call set_gravity_grid
      do j = 2-Ng,N+Ng
         rho_a(j) = cell_average(rho_column, j)
      enddo
      rho_a(1-Ng) = rho_a(2-Ng)
      call discrete_equilibrium_density(rho_a(1), rho_d)
      p_d  = nhat*rho_d
      span = r(N) - r(1)

      open(newunit=iu4, file=pt_file, position='append',                 &
           status='unknown')

      do is = 1,n_scheme

         call set_scheme(is)
         well_balanced = .false.

         ! ---- zero gravity, uniform pressure, at rest ----
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
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         z_frac  = 0.0d0
         z_scale = 0.0d0
         do j = 3,N-2
            dAp = r_edg(j)*r_edg(j)
            dAm = r_edg(j-1)*r_edg(j-1)
            dV  = spherical_cell_volume(j)
            geo = (dAp - dAm)*1.0d0/dV
            sc  = residual_row_scale(2, j, u)
            z_frac  = max(z_frac,  sc/geo)
            z_scale = max(z_scale, sc)
         enddo

         ! ---- the discrete equilibrium: the scale is the weight ----
         b0 = b0_case
         call set_gravity_grid
         do j = 1-Ng,N+Ng
            u(1,j) = rho_d(j)
            u(2,j) = 0.0d0
            u(3,j) = p_d(j)/(gamma_ad - 1.0d0)
         enddo
         call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
         n_part_cell1 = nhat*rho_d(1)
         call Apply_BC(u)
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)

         ! The weight of the base scheme, as `source` forms it from the two
         ! reconstructed face densities of this cell.
         !
         ! WHAT IS ASSERTED IS ONE SIDED: the scale must not fall BELOW the
         ! weight, (w - s)/w <= 0, which says that the weight is one of the
         ! terms the max runs over.  A scale made of the discretization's
         ! pieces has no such term under PLM, where S_2 is the weight minus
         ! the geometric pressure term, and it comes out 28 percent below
         ! the weight on this state (MEASURED on the entry text).  The
         ! departure |s - w|/w is REPORTED and not asserted, because it is
         ! not round-off on this state for the base scheme: the state is the
         ! equilibrium of the WELL-BALANCED reconstruction, on which the base
         ! scheme's own row is its truncation error, and the HLLC face
         ! pressure is one side's own value rather than an average, so the
         ! momentum flux carries a dissipative part of the size dr/H times
         ! the weight at rest.  Both fall with the grid; a bound on that
         ! quantity would be a tolerance chosen for them.
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         e_base  = 0.0d0
         b_base  = -1.0d0
         do j = 3,N-2
            wgt = abs(0.5d0*(WR(1,j-1) + WL(1,j))                        &
                      *(Gphi_i(j) - Gphi_i(j-1))/dr_j(j))
            sc  = residual_row_scale(2, j, u)
            if (wgt .gt. 0.0d0) then
               e_base = max(e_base, abs(sc - wgt)/wgt)
               b_base = max(b_base, (wgt - sc)/wgt)
            endif
         enddo

         ! The same state under the key, whose weight is the pressure force
         ! of the cell's own equilibrium (the right-hand side of the
         ! cancellation identity, in the form the assembled row belongs to).
         well_balanced = .true.
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         e_wb   = 0.0d0
         do j = 3,N-2
            dAp = r_edg(j)*r_edg(j)
            dAm = r_edg(j-1)*r_edg(j-1)
            dV  = spherical_cell_volume(j)
            if (use_plm) then
               epf = abs(u(1,j)*(dAp*(Gphi_i(j)   - Gphi_c(j))           &
                               + dAm*(Gphi_c(j)   - Gphi_i(j-1))))/dV
            else
               epf = abs(u(1,j)*(Gphi_i(j) - Gphi_i(j-1)))/dr_j(j)
            endif
            sc = residual_row_scale(2, j, u)
            if (epf .gt. 0.0d0) e_wb = max(e_wb,                       &
               (abs(sc - epf) - abs(dF(2,j) - S(2,j)))/epf)
         enddo
         well_balanced = .false.

         ! ---- the perturbed equilibrium, base scheme: the scaled
         !      imbalance the perturbation adds halves with it ----
         ! WHAT IS DIFFERENCED, AND WHY.  The base scheme's row on the
         ! well-balanced discrete equilibrium is not zero: that state is the
         ! equilibrium of that reconstruction, and the base scheme's
         ! own truncation error on it is eps-independent and larger than the
         ! perturbation's imbalance here (MEASURED on the entry text: the
         ! undifferenced ratio is 1.0012 at N = 250 and 1.586 at N = 2000,
         ! the truncation falling with the grid).  The perturbation's own
         ! contribution to each cell's row is what has to halve, so the
         ! unperturbed row of the same cell is subtracted; it and every face
         ! quantity built from it are homogeneous of degree one in the
         ! perturbation, so the ratio is 2 up to O(eps).  A scale that is
         ! the row itself, which is what the well-balanced key had before
         ! the weight was put back, returns 1 here and not 2.
         do j = 1-Ng,N+Ng
            u(1,j) = rho_d(j)
            u(2,j) = 0.0d0
            u(3,j) = p_d(j)/(gamma_ad - 1.0d0)
         enddo
         call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
         n_part_cell1 = nhat*rho_d(1)
         call Apply_BC(u)
         call Reconstruct(u, WL, WR)
         call RK_rhs(u, WL, WR, dF, S)
         do j = 1-Ng,N+Ng
            R0(j) = dF(2,j) - S(2,j)
         enddo

         x1 = 0.0d0
         x2 = 0.0d0
         do ip = 1,2
            do j = 1-Ng,N+Ng
               shape_p = 1.0d0 + eps_p/dble(2**(ip-1))                   &
                         *sin(4.0d0*twopi*(r(j) - r(1))/span)
               u(1,j) = rho_d(j)
               u(2,j) = 0.0d0
               u(3,j) = shape_p*p_d(j)/(gamma_ad - 1.0d0)
            enddo
            call set_base_reservoir(p_column(1.0d0), 1.0d0, nhat, 1.0d0)
            n_part_cell1 = nhat*rho_d(1)
            call Apply_BC(u)
            call Reconstruct(u, WL, WR)
            call RK_rhs(u, WL, WR, dF, S)
            s_break = residual_row_scale(2, 3, 1.5d0*u)
            sc = 0.0d0
            do j = 3,N-2
               sc = max(sc, abs(dF(2,j) - S(2,j) - R0(j))                &
                            /residual_row_scale(2, j, u))
            enddo
            if (ip .eq. 1) x1 = sc
            if (ip .eq. 2) x2 = sc
         enddo
         p_rat = 0.0d0
         if (x2 .gt. 0.0d0) p_rat = x1/x2

         ! ---- supersonic uniform flow, no gravity, uniform pressure ----
         b0 = 0.0d0
         call set_gravity_grid
         cs_u = sqrt(gamma_ad*1.0d0/1.0d0)
         v_u  = mach*cs_u
         do j = 1-Ng,N+Ng
            u(1,j) = 1.0d0
            u(2,j) = v_u
            u(3,j) = 0.5d0*v_u*v_u + 1.0d0/(gamma_ad - 1.0d0)
         enddo
         n_part_cell1 = 1.0d0
         call set_base_reservoir(1.0d0, 1.0d0, 1.0d0, 1.0d0)
         call Apply_BC(u)
         s_break = residual_row_scale(2, 3, 1.5d0*u)
         u_ram = 0.0d0
         do j = 3,N-2
            dAp  = r_edg(j)*r_edg(j)
            dAm  = r_edg(j-1)*r_edg(j-1)
            dV   = spherical_cell_volume(j)
            ramr = (dAp - dAm)*v_u*v_u/dV
            sc   = residual_row_scale(2, j, u)
            u_ram = max(u_ram, abs(sc - ramr)/ramr)
         enddo

         b0 = b0_case
         call set_gravity_grid

         write(iu4,'(I8,I3,7ES24.16)') N, is, z_frac, z_scale, e_base,   &
            b_base, e_wb, p_rat, u_ram

      enddo

      close(iu4)
      well_balanced = sav_wb

      end subroutine measure_momentum_physical_terms

      ! ------------------------------------------------------------------ !

      subroutine report_momentum_physical_terms(nfail)
      ! The momentum row's reference scale against the terms of the momentum
      ! equation, on the four states measure_momentum_physical_terms builds.
      integer, intent(inout) :: nfail

      integer, parameter :: mxs = 256
      integer :: Nv(mxs), sv(mxs)
      real*8  :: zf(mxs), zs(mxs), eb(mxs), bb(mxs), ea(mxs), pr(mxs)
      real*8  :: ur(mxs)
      integer :: nrec, iu, ios, i, is
      real*8  :: wzf, wzs, web, wbb, wea, wpr, wur

      nrec = 0
      open(newunit=iu, file=pt_file, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') 'FAIL hydrostatic_momentum_terms_file '//        &
              'measured=missing reference='//pt_file//' tol=0'
         nfail = nfail + 1
         return
      endif
      do
         if (nrec .ge. mxs) exit
         read(iu,*,iostat=ios) Nv(nrec+1), sv(nrec+1), zf(nrec+1),       &
              zs(nrec+1), eb(nrec+1), bb(nrec+1), ea(nrec+1), pr(nrec+1), &
              ur(nrec+1)
         if (ios .ne. 0) exit
         nrec = nrec + 1
      enddo
      close(iu)

      write(*,'(A)') ''
      write(*,'(A)') '  DIAGNOSTIC momentum row scale against the terms'//&
           ' of the momentum equation, interior cells 3..N-2'
      write(*,'(A)') '        N  scheme      zero-g/(2p/r)'//            &
           '     zero-g scale   base |s-w|/w'//                          &
           '    base (w-s)/w      WB |s-w|/w'//                          &
           '   halving ratio     supersonic ram'
      do i = 1,nrec
         write(*,'(A,I8,2X,A10,2X,5ES16.5,F15.6,ES17.5)')                &
            '     ', Nv(i), trim(scheme_name(sv(i))), zf(i), zs(i),      &
            eb(i), bb(i), ea(i), pr(i), ur(i)
      enddo

      ! THE BOUNDS.  The two identity statements -- no force at all, and a
      ! pure ram divergence -- are asserted at 1e-8.  That is eight decades
      ! below what a scale built from the discretization's pieces returns on
      ! those states (1.000 of the geometric pressure term, and 1/(gamma M M)
      ! = 6.7e-2 above the ram term at M = 3), and above the rounding the
      ! scale's own assembly leaves: under WENO3 the reconstruction of a
      ! uniform pressure is exact and the scale is the tiny floor itself
      ! (MEASURED 2.22507E-308, which is tiny(1.0d0)), while under PLM the
      ! geometric term and the flux pressure cancel to the rounding of the
      ! cell pressure.  On the discrete equilibrium the base scheme's
      ! statement is one sided, (w - s)/w <= 0, which says the weight is one
      ! of the terms the max runs over and is an inequality of the assembly
      ! and not a tolerance; the key's is two sided and exact, the scale
      ! being the equilibrium pressure force itself (MEASURED 0.00000E+00
      ! at every N and every pair).  The departure |s - w|/w of the base
      ! scheme is reported, not asserted: on this state, which is the
      ! equilibrium of the well-balanced reconstruction, it carries the base
      ! scheme's truncation error and, with HLLC, the dissipative part of a
      ! momentum flux whose face pressure is one side's own value.  The halving
      ! ratio is asserted against two to a tenth, which is an algebraic
      ! property and not a tolerance: the imbalance the perturbation adds
      ! and every face quantity built from it are homogeneous of degree one
      ! in it (the limited slope is positively homogeneous, the WENO3
      ! weights are scale free), so the ratio is 2 up to O(eps).
      do is = 1,n_scheme
         wzf = 0.0d0
         wzs = 0.0d0
         web = 0.0d0
         wbb = -1.0d0
         wea = 0.0d0
         wur = 0.0d0
         wpr = 0.0d0
         do i = 1,nrec
            if (sv(i) .ne. is) cycle
            wzf = max(wzf, zf(i))
            wzs = max(wzs, zs(i))
            web = max(web, eb(i))
            wbb = max(wbb, bb(i))
            wea = max(wea, ea(i))
            wur = max(wur, ur(i))
            wpr = max(wpr, abs(pr(i) - 2.0d0))
         enddo
         call verdict_below('zero_gravity_uniform_state_momentum_'//     &
                            'scale_vanishes['//                          &
                            trim(scheme_name(is))//']', wzf, 1.0d-8,     &
                            nfail)
         call verdict_below('momentum_scale_is_not_below_the_weight_'// &
                            'on_the_discrete_equilibrium['//             &
                            trim(scheme_name(is))//']', wbb, 1.0d-12,    &
                            nfail)
         call verdict_below('well_balanced_momentum_scale_is_the_'//     &
                            'equilibrium_pressure_force['//              &
                            trim(scheme_name(is))//']', wea, 1.0d-12,    &
                            nfail)
         call verdict_below('momentum_scaled_row_tracks_the_'//          &
                            'imbalance['//                              &
                            trim(scheme_name(is))//']', wpr, 1.0d-1,     &
                            nfail)
         call verdict_below('supersonic_uniform_flow_momentum_scale_'//  &
                            'is_the_ram_divergence['//                   &
                            trim(scheme_name(is))//']', wur, 1.0d-8,     &
                            nfail)
         write(*,'(A,A,A,ES12.5,A,ES12.5)') '  DIAGNOSTIC zero_gravity'//&
              '_scale[', trim(scheme_name(is)), '] = ', wzs,             &
              '   equilibrium |s-w|/w = ', web
      enddo

      end subroutine report_momentum_physical_terms

      ! ------------------------------------------------------------------ !

      subroutine report_momentum_row_scale(nfail)
      ! The normalized momentum row of the well-balanced key, on the three
      ! states measure_momentum_row_scale builds.
      integer, intent(inout) :: nfail

      integer, parameter :: mxs = 256
      integer :: Nv(mxs), sv(mxs)
      real*8  :: seq(mxs), sp1(mxs), sp2(mxs), szr(mxs), sow(mxs)
      integer :: nrec, iu, ios, i, is
      real*8  :: weq, wp, wzr, wown, wrat, rat

      nrec = 0
      open(newunit=iu, file=sc_file, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') 'FAIL hydrostatic_momentum_scale_file '//        &
              'measured=missing reference='//sc_file//' tol=0'
         nfail = nfail + 1
         return
      endif
      do
         if (nrec .ge. mxs) exit
         read(iu,*,iostat=ios) Nv(nrec+1), sv(nrec+1), seq(nrec+1),      &
              sp1(nrec+1), sp2(nrec+1), szr(nrec+1), sow(nrec+1)
         if (ios .ne. 0) exit
         nrec = nrec + 1
      enddo
      close(iu)

      write(*,'(A)') ''
      write(*,'(A)') '  DIAGNOSTIC momentum row over its own reference'// &
           ' scale, well-balanced key ON, interior cells 3..N-2'
      write(*,'(A)') '        N  scheme        equilibrium'//            &
           '    perturbed(eps)  perturbed(eps/2)   ratio'//              &
           '   zero-g vs base  zero-g vs own terms'
      do i = 1,nrec
         rat = 0.0d0
         if (sp2(i) .gt. 0.0d0) rat = sp1(i)/sp2(i)
         write(*,'(A,I8,2X,A10,2X,3ES17.5,F10.4,2ES17.5)')               &
            '     ', Nv(i), trim(scheme_name(sv(i))), seq(i), sp1(i),    &
            sp2(i), rat, szr(i), sow(i)
      enddo

      ! THE BOUNDS.  On the discrete equilibrium the key's momentum row is
      ! at the rounding level of a row assembled from cell pressures,
      ! epsilon x p x r^2/dV, which is 1e-12 of the cell's own weight at the
      ! base of this grid and far below it above; 1e-10 is above that bound
      ! and twelve decades below the one a degenerate scale returns.  The
      ! perturbed state is asserted against one, which is what a row that is
      ! its own scale reads, and its two perturbations against the factor
      ! two between them, to a tenth: the departure is homogeneous of degree
      ! one in the perturbation, so the ratio is 2 up to O(eps).  The
      ! zero-gravity statement is an identity and not a bound: with no
      ! gravity the equilibrium pressure force is zero, so the scale is the
      ! max of the same two numbers the row's own terms give and the
      ! departure is exactly zero.  It is exactly zero only while this
      ! fixture divides by the shell volume the rows divide by, which is why
      ! it reads spherical_cell_volume and does not spell the volume out:
      ! the two spellings of it differ by up to 6e-13 at the base of a
      ! refined grid, which this row reads as a broken identity.
      do is = 1,n_scheme
         weq  = 0.0d0
         wp   = 0.0d0
         wzr  = 0.0d0
         wown = 0.0d0
         wrat = 0.0d0
         do i = 1,nrec
            if (sv(i) .ne. is) cycle
            weq = max(weq, seq(i))
            wp  = max(wp,  max(sp1(i), sp2(i)))
            wzr = max(wzr, szr(i))
            wown = max(wown, sow(i))
            rat = 0.0d0
            if (sp2(i) .gt. 0.0d0) rat = sp1(i)/sp2(i)
            wrat = max(wrat, abs(rat - 2.0d0))
         enddo
         call verdict_below('well_balanced_scaled_momentum_row_on_its'// &
                            '_own_equilibrium['//                        &
                            trim(scheme_name(is))//']', weq, 1.0d-10,    &
                            nfail)
         call verdict_below('well_balanced_scaled_momentum_row_below'//  &
                            '_one_when_perturbed['//                     &
                            trim(scheme_name(is))//']', wp, 1.0d0,       &
                            nfail)
         call verdict_below('well_balanced_scaled_momentum_row_tracks'// &
                            '_the_imbalance['//                          &
                            trim(scheme_name(is))//']', wrat, 1.0d-1,    &
                            nfail)
         call verdict_below('zero_gravity_momentum_scale_is_the_rows'//  &
                            '_own_terms['//                              &
                            trim(scheme_name(is))//']', wown, 1.0d-15,   &
                            nfail)
         write(*,'(A,A,A,ES12.5)') '  DIAGNOSTIC zero_gravity_scale_'//  &
              'against_the_base_scheme[', trim(scheme_name(is)),         &
              '] = ', wzr
      enddo

      end subroutine report_momentum_row_scale

      ! ------------------------------------------------------------------ !

      subroutine report_ladder(nfail)
      ! Read every record written by the invocations that measured one grid
      ! each, print the ladder with its observed orders, and give the
      ! verdicts.
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

      call report_well_balanced(nfail)
      call report_momentum_row_scale(nfail)
      call report_momentum_physical_terms(nfail)
      call report_cfl_restriction(nfail)

      end subroutine report_ladder

      ! ------------------------------------------------------------------ !

      subroutine report_cfl_restriction(nfail)
      ! WHICH CELLS AND WHICH FACES THE EXPLICIT-STABLE INTERVAL IS TAKEN
      ! OVER, on the analytic column of this program at rest.
      !
      ! The restriction eval_dt returns belongs to the EVOLVED cells: each
      ! one is bounded by its own shell volume, its two face areas and the
      ! fastest signal of each face's Riemann problem,
      !
      !   dt_j = CFL V_j / max( A_{j-1/2} S_{j-1/2}, A_{j+1/2} S_{j+1/2} ),
      !   S_{j+1/2} = max(|v| + c) over the two states bounding the face.
      !
      ! Three things are asserted, and each one separates that statement
      ! from a plausible neighbor:
      !
      !   the value itself, against the expression above formed here, which
      !   pins the geometry (a cell-centered dr_j/(|v| + c) misses it by
      !   dr/r, 2e-4 on this base grid);
      !
      !   that a cell which bounds NO face of an evolved cell cannot
      !   restrict the interval, asserted by giving the outermost ghost a
      !   sound speed far above every other cell's and requiring the
      !   interval not to move at all;
      !
      !   that the base face's own Riemann problem DOES restrict cell 1,
      !   asserted by giving the base ghost that sound speed and requiring
      !   the interval to become CFL V_1/(A_1/2 S_1/2) with it. This is the
      !   half a "take the minimum over 1..N" would drop: no ghost cell is
      !   evolved, but the reservoir can carry the fastest wave into the
      !   first cell.
      integer, intent(inout) :: nfail

      real*8, dimension(:,:), allocatable :: W
      real*8, dimension(:),   allocatable :: dt_loc
      real*8  :: dt, dt_expect, dt_ghost, dt_base, s_lo, s_hi, s_big
      real*8  :: dtj, b0_case
      integer :: j, jbind

      N  = 500
      T0 = Teq
      R0 = Rp_RJ*RJ
      Mp = Mp_MJ*MJ
      v0 = sqrt(kb_erg*T0/mu)
      b0_case = (Gc*Mp*mu)/(kb_erg*T0*R0)

      call allocate_grid_arrays
      call define_grid
      b0 = b0_case
      call set_gravity_grid

      allocate(W(3,1-Ng:N+Ng), dt_loc(1-Ng:N+Ng))
      do j = 1-Ng,N+Ng
         W(1,j) = rho_column(r(j))
         W(2,j) = 0.0d0
         W(3,j) = p_column(r(j))
      enddo

      ! The column is at rest and isothermal, so every signal speed is the
      ! adiabatic sound speed of the same temperature and the binding cell
      ! is the one with the smallest V_j/A_{j+1/2}: cell 1 on this grid.
      call eval_dt(W, dt, dt_loc)

      dt_expect = huge(1.0d0)
      jbind     = 1
      do j = 1,N
         s_lo = max(signal_speed_of(W(:,j-1)), signal_speed_of(W(:,j)))
         s_hi = max(signal_speed_of(W(:,j)),   signal_speed_of(W(:,j+1)))
         dtj  = CFL*spherical_cell_volume(j)                              &
                /max(r_edg(j-1)*r_edg(j-1)*s_lo, r_edg(j)*r_edg(j)*s_hi)
         if (dtj .lt. dt_expect) then
            dt_expect = dtj
            jbind     = j
         endif
      enddo
      write(*,'(A,I5,A,ES13.6,A,ES13.6)')                                 &
           '  cfl: binding cell ', jbind, '  dt=', dt,                     &
           '  cell-centred CFL dr/(|v|+c) of that cell=',                  &
           CFL*dr_j(jbind)/signal_speed_of(W(:,jbind))
      call verdict_below('cfl_is_the_face_restriction_of_its_own_volume',  &
           abs(dt - dt_expect)/dt_expect, 1.0d-14, nfail)

      ! A cell beyond every face of an evolved cell: the outermost ghost.
      ! Ng is at least 2 here, so r_edg(N) is bounded by cells N and N+1
      ! and cell N+Ng bounds no face of cell N.
      s_big  = 1.0d3*signal_speed_of(W(:,1))
      W(3,N+Ng) = W(1,N+Ng)*s_big*s_big/gamma_ad
      call eval_dt(W, dt_ghost, dt_loc)
      call verdict_below('cfl_ignores_a_cell_bounding_no_evolved_face',    &
           abs(dt_ghost - dt)/dt, 1.0d-15, nfail)
      W(3,N+Ng) = p_column(r(N+Ng))

      ! The base face. Its Riemann problem is between the lowest ghost and
      ! cell 1, and its wave restricts cell 1 through the face area of that
      ! face and the volume of that cell.
      W(3,0) = W(1,0)*s_big*s_big/gamma_ad
      call eval_dt(W, dt_base, dt_loc)
      dt_expect = CFL*spherical_cell_volume(1)                            &
                  /(r_edg(0)*r_edg(0)*max(signal_speed_of(W(:,0)),         &
                                          signal_speed_of(W(:,1))))
      call verdict_below('cfl_reads_the_base_face_riemann_problem',        &
           abs(dt_base - dt_expect)/dt_expect, 1.0d-14, nfail)
      W(3,0) = p_column(r(0))

      deallocate(W, dt_loc)

      end subroutine report_cfl_restriction

      ! ------------------------------------------------------------------ !

      real*8 function signal_speed_of(Wj) result(c)
      ! |v| + c of one primitive state, the bound on the fastest wave of a
      ! Riemann problem it takes part in, with the adiabatic sound speed of
      ! the atomic column this program builds (no caloric mixture is
      ! installed here, so gamma_ad is the index eval_dt uses too).
      real*8, intent(in) :: Wj(3)
      c = abs(Wj(2)) + sqrt(gamma_ad*Wj(3)/Wj(1))
      end function signal_speed_of

      ! ------------------------------------------------------------------ !

      subroutine report_well_balanced(nfail)
      ! The well-balanced table, both ways: the discrete equilibrium of the
      ! scheme, on which the key must return the rounding level, and the
      ! analytic column, on which the key and the base scheme both carry the
      ! truncation error of the discretization of gravity.
      integer, intent(inout) :: nfail

      integer, parameter :: mxw = 512
      integer :: Nv(mxw), sv(mxw), wv(mxw), cv(mxw)
      real*8  :: mom(mxw), mas(mxw), ene(mxw)
      integer :: nrec, iu, ios, i, is
      real*8  :: wmom, wmas, wene
      character(len=10) :: colname(2)

      colname(1) = 'analytic'
      colname(2) = 'discrete'

      nrec = 0
      open(newunit=iu, file=wb_file, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') 'FAIL hydrostatic_wellbalanced_file '//           &
              'measured=missing reference='//wb_file//' tol=0'
         nfail = nfail + 1
         return
      endif
      do
         if (nrec .ge. mxw) exit
         read(iu,*,iostat=ios) Nv(nrec+1), sv(nrec+1), wv(nrec+1),        &
              cv(nrec+1), mom(nrec+1), mas(nrec+1), ene(nrec+1)
         if (ios .ne. 0) exit
         nrec = nrec + 1
      enddo
      close(iu)

      write(*,'(A)') ''
      write(*,'(A)') '  DIAGNOSTIC well-balanced rows, normalized by the'//&
           ' size of their own terms, interior cells 3..N-2'
      write(*,'(A)') '        N  scheme        column    Well balanced'// &
           '     momentum         mass          energy'
      do i = 1,nrec
         write(*,'(A,I8,2X,A10,2X,A8,4X,L1,7X,3ES15.4)')                  &
            '     ', Nv(i), trim(scheme_name(sv(i))), trim(colname(cv(i))),&
            wv(i) .eq. 1, mom(i), mas(i), ene(i)
      enddo

      ! THE ASSERTION.  On the discrete equilibrium of the scheme the key
      ! must leave the rounding level in every row, at every N and for both
      ! reconstructions.  1e-13 is not a tolerance chosen to pass: the
      ! arithmetic bound of a row assembled from cell pressures is
      ! epsilon x p x r^2/dV, which on this grid is 1e-12 of the momentum
      ! weight at the base and far below it in the interior.
      do is = 1,n_scheme
         wmom = 0.0d0
         wmas = 0.0d0
         wene = 0.0d0
         do i = 1,nrec
            if (sv(i) .ne. is .or. wv(i) .ne. 1 .or. cv(i) .ne. 2) cycle
            wmom = max(wmom, mom(i))
            wmas = max(wmas, mas(i))
            wene = max(wene, ene(i))
         enddo
         call verdict_below('well_balanced_momentum_row_at_rounding['//   &
                            trim(scheme_name(is))//']', wmom, 1.0d-13,    &
                            nfail)
         call verdict_below('well_balanced_mass_row_at_rounding['//       &
                            trim(scheme_name(is))//']', wmas, 1.0d-13,    &
                            nfail)
         call verdict_below('well_balanced_energy_row_at_rounding['//     &
                            trim(scheme_name(is))//']', wene, 1.0d-13,    &
                            nfail)
      enddo

      end subroutine report_well_balanced

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

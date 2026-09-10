      module hydrodynamic_rows_double
      ! The hydrodynamic rows of the stationary residual in DOUBLE
      ! precision, from the kind-generic text of
      ! hydrodynamic_rows_body.inc.  This instantiation exists to be
      ! compared: it must reproduce Reconstruction_step, Numerical_Fluxes,
      ! source_func and RK_integration bit for bit on any state, which is
      ! what makes the quadruple instantiation below a measurement of the
      ! ARITHMETIC and not of a second discretization.

      use global_parameters
      use caloric_eos, only: caloric_cell_mixture,                       &
                             h2_rovibrational_table_grid,                &
                             h2_rovibrational_table_node
      use BC_Apply, only: base_face_W, base_face_lower_W
      use Reconstruction_step, only: S0sav, S1sav
      use S_estimate_ROE, only: ROE_STAR_OK, ROE_STAR_VACUUM,            &
                                ROE_STAR_INADMISSIBLE

      implicit none
      private
      public :: hydrodynamic_flux_difference_and_source

      integer, parameter :: fv = kind(1.0d0)
      ! Stopping rule of the caloric energy -> temperature inverse, in the
      ! units of this kind: the root is resolved to the precision the rows
      ! are evaluated in, and no tighter, so the Newton terminates.
      real(fv), parameter :: fv_newton_tol = 1.0e-15_fv
      integer,  parameter :: fv_newton_maxit = 80

      contains

      include 'hydrodynamic_rows_body.inc'

      end module hydrodynamic_rows_double

      !-----------------------------------------------------------!

      module hydrodynamic_rows_quadruple
      ! The same rows in QUADRUPLE precision, from the same text.  The
      ! composition, the heating and the cooling remain double: what is
      ! evaluated here is the arithmetic of the finite-volume operator on a
      ! given composition, which is where the residual's non-smoothness
      ! floor is (docs/Update_EXHALE.md item N33: the Roe flux of a face in
      ! the near-hydrostatic base layer is built from state jumps 3e5 to
      ! 6e6 times smaller than the states, so it carries the last bit of
      ! O(1) quantities, and the row divides the flux difference by the cell
      ! volume, multiplying that by r^2/dV = 5.1e3 at the base).  Eighteen
      ! decades of extra mantissa put that floor below every other quantity
      ! of the row.

      use global_parameters
      use caloric_eos, only: caloric_cell_mixture,                       &
                             h2_rovibrational_table_grid,                &
                             h2_rovibrational_table_node
      use BC_Apply, only: base_face_W, base_face_lower_W
      use Reconstruction_step, only: S0sav, S1sav
      use S_estimate_ROE, only: ROE_STAR_OK, ROE_STAR_VACUUM,            &
                                ROE_STAR_INADMISSIBLE

      implicit none
      private
      public :: hydrodynamic_flux_difference_and_source

      integer, parameter :: fv = selected_real_kind(30, 300)
      real(fv), parameter :: fv_newton_tol = 1.0e-31_fv
      integer,  parameter :: fv_newton_maxit = 200

      contains

      include 'hydrodynamic_rows_body.inc'

      end module hydrodynamic_rows_quadruple

      !-----------------------------------------------------------!

      module hydrodynamic_rows
      ! WHICH ARITHMETIC THE STATIONARY RESIDUAL'S HYDRODYNAMIC ROWS ARE
      ! ASSEMBLED IN, and the arm that changes it.
      !
      ! Default: the production routines (Reconstruct, RK_rhs, Num_flux,
      ! source), in double, reached through
      ! steady_residual_mod.reconstruction_continuation_rhs.  Nothing here
      ! is called and nothing is different.
      !
      ! EXHALE_RESID_QUAD=1 arms a CONTROL EXPERIMENT.  It replaces those
      ! four routines, for the stationary evaluations only, by the
      ! quadruple-precision instantiation of one kind-generic text, from the
      ! conserved state through the reconstruction (the PLM limiter, the
      ! WENO3 weights), the Riemann flux (Roe with its entropy fix, the HLLE
      ! fallback), the geometric and gravitational sources and the division
      ! by the cell volume, rounded to double once, at the row.
      !
      ! WHAT IT IS FOR.  The stationary residual carries a non-smoothness
      ! floor of 1.0e-11 to 1.6e-11 along a preconditioned Krylov direction,
      ! flat over 5.4 decades of sampling spacing, and that floor is the
      ! ROUNDING of the flux assembly (item N33): the base-layer face states
      ! cancel to 3e5 to 6e6, so the interface flux carries the last bit of
      ! O(1) quantities, and the row multiplies it by r^2/dV.  Divided by
      ! the small increment a preconditioned direction produces, that floor
      ! IS the additivity defect of the finite-difference Jacobian action
      ! (item N31), which is why the Krylov cycle spends every product at a
      ! true relative residual near one.  If the same residual with the
      ! floor eighteen decades lower lets the Krylov cycle converge, the
      ! diagnosis is confirmed before any discretization is changed.
      !
      ! IT IS NOT A DISCRETIZATION AND NOTHING IS ADOPTED FROM IT.  The
      ! equations, the scheme, the boundary construction and the
      ! composition are the same; only the width of the mantissa the
      ! arithmetic runs in is different, and the arm is off unless the
      ! environment names it.  The physically right treatment of the floor
      ! is a well-balanced flux difference at the base, not extra precision.
      !
      ! Restrictions of the arm, each of which stops the run rather than
      ! silently evaluating something else: the low-Mach contact-mode
      ! dissipation flux is not part of the generic text, and the
      ! reconstruction flags must be the consistent pair the schemes are
      ! selected by.

      use global_parameters
      use RK_integration, only: face_flux, face_p
      use Reconstruction_step, only: n_faces_positivity_limited
      use Numerical_Fluxes, only: n_faces_roe_hlle, n_faces_llf
      use low_mach_dissipation, only: low_mach_damping_active
      use hydrodynamic_rows_double,                                      &
             only: rows_in_double => hydrodynamic_flux_difference_and_source
      use hydrodynamic_rows_quadruple,                                   &
             only: rows_in_quad   => hydrodynamic_flux_difference_and_source

      implicit none
      private
      public :: generic_precision_rows_arm
      public :: ARM_OFF, ARM_QUADRUPLE, ARM_GENERIC_DOUBLE
      public :: hydrodynamic_rows_in_quadruple_precision
      public :: hydrodynamic_rows_in_double_precision
      public :: quad_rows_calls, quad_rows_seconds

      ! What EXHALE_RESID_QUAD selects.
      integer, parameter :: ARM_OFF            = 0  ! the production routines
      integer, parameter :: ARM_QUADRUPLE      = 1  ! the generic text at real(16)
      integer, parameter :: ARM_GENERIC_DOUBLE = 2  ! the generic text at real(8)

      ! How many times the quadruple pipeline ran and how long it spent
      ! there, so the cost of the arm is a measured number and not an
      ! estimate.
      integer :: quad_rows_calls   = 0
      real*8  :: quad_rows_seconds = 0.0d0

      contains

      ! ------------------------------------------------------!

      integer function generic_precision_rows_arm() result(mode)
      ! EXHALE_RESID_QUAD: absent or anything else ARM_OFF, 1
      ! ARM_QUADRUPLE, 2 ARM_GENERIC_DOUBLE.  Read at every call rather than
      ! cached, so that a test can arm and disarm it inside one process; the
      ! read costs microseconds against a residual evaluation of
      ! milliseconds.
      !
      ! Mode 2 is the control OF the control: the same generic text at the
      ! double kind, which must reproduce the production routines bit for
      ! bit.  A whole run under it that is bitwise the production run is the
      ! strongest form of that assertion, because it exercises every state
      ! the solve visits and not one stated by a test.
      character(len=32) :: env
      integer :: st, ln
      mode = ARM_OFF
      env = ' '
      call get_environment_variable('EXHALE_RESID_QUAD', env,            &
                                    length=ln, status=st)
      if (st .ne. 0 .or. ln .le. 0) return
      select case (trim(adjustl(env)))
      case ('1')
         mode = ARM_QUADRUPLE
      case ('2')
         mode = ARM_GENERIC_DOUBLE
      end select
      end function generic_precision_rows_arm

      ! ------------------------------------------------------!

      subroutine hydrodynamic_rows_in_double_precision(u,WL,WR,dF,S)
      ! The generic text at kind(1.0d0).  Not on any run's path: it is the
      ! control of the arm, compared against the production routines.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL,WR,dF,S
      real*8, dimension(3,1-Ng:N+Ng) :: ff
      real*8, dimension(1-Ng:N+Ng)   :: fp
      integer :: n_scaled, n_hlle, n_llf
      call the_generic_text_is_usable()
      call rows_in_double(u,WL,WR,dF,S,ff,fp,n_scaled,n_hlle,n_llf)
      call store_the_interface_fluxes(ff,fp,n_scaled,n_hlle,n_llf)
      end subroutine hydrodynamic_rows_in_double_precision

      ! ------------------------------------------------------!

      subroutine hydrodynamic_rows_in_quadruple_precision(u,WL,WR,dF,S)
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL,WR,dF,S
      real*8, dimension(3,1-Ng:N+Ng) :: ff
      real*8, dimension(1-Ng:N+Ng)   :: fp
      integer :: n_scaled, n_hlle, n_llf
      real*8  :: t_a, t_b
      call the_generic_text_is_usable()
      call cpu_seconds(t_a)
      call rows_in_quad(u,WL,WR,dF,S,ff,fp,n_scaled,n_hlle,n_llf)
      call cpu_seconds(t_b)
      quad_rows_calls   = quad_rows_calls + 1
      quad_rows_seconds = quad_rows_seconds + (t_b - t_a)
      call store_the_interface_fluxes(ff,fp,n_scaled,n_hlle,n_llf)
      end subroutine hydrodynamic_rows_in_quadruple_precision

      ! ------------------------------------------------------!

      subroutine store_the_interface_fluxes(ff,fp,n_scaled,n_hlle,n_llf)
      ! The interface flux and face pressure the rows were differenced from,
      ! in the arrays every later reader of them uses (the species face
      ! flux, the flux gate, the jump scan), and the three counters of the
      ! run summary.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: ff
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: fp
      integer, intent(in) :: n_scaled, n_hlle, n_llf
      if (.not. allocated(face_flux)) then
         allocate(face_flux(3,1-Ng:N+Ng), face_p(1-Ng:N+Ng))
      endif
      face_flux = ff
      face_p    = fp
      n_faces_positivity_limited = n_faces_positivity_limited + n_scaled
      n_faces_roe_hlle           = n_faces_roe_hlle + n_hlle
      n_faces_llf                = n_faces_llf + n_llf
      end subroutine store_the_interface_fluxes

      ! ------------------------------------------------------!

      subroutine the_generic_text_is_usable()
      ! The two configurations the generic text does not carry.  Both stop
      ! the run: an arm that quietly evaluated a different operator would
      ! answer a question nobody asked.
      if (low_mach_damping_active()) then
         write(*,*) 'ERROR: EXHALE_RESID_QUAD does not carry the',       &
                    ' low-Mach contact-mode dissipation flux.'
         error stop 1
      endif
      if ((use_plm .eqv. use_weno3) .or.                                 &
          (use_plm .and. rec_method .ne. 'PLM') .or.                     &
          (use_weno3 .and. rec_method .ne. 'WENO3')) then
         write(*,*) 'ERROR: EXHALE_RESID_QUAD needs a consistent',       &
                    ' reconstruction selection; rec_method = ',          &
                    trim(rec_method)
         error stop 1
      endif
      end subroutine the_generic_text_is_usable

      ! ------------------------------------------------------!

      subroutine cpu_seconds(t)
      ! Wall-clock seconds from the system clock, for the cost measurement.
      real*8, intent(out) :: t
      integer(kind=8) :: c, cr
      call system_clock(count=c, count_rate=cr)
      if (cr .gt. 0_8) then
         t = dble(c)/dble(cr)
      else
         t = 0.0d0
      endif
      end subroutine cpu_seconds

      end module hydrodynamic_rows

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
      ! floor is (docs/Update_EXHALE_stage2.pdf item N33: the Roe flux of a face in
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
      ! ASSEMBLED IN, and what changes it.
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
      ! arithmetic runs in is different, and it is off unless the
      ! environment names it.  The physically right treatment of the floor
      ! is a well-balanced flux difference at the base, not extra precision.
      !
      ! Restrictions of this assembly, each of which stops the run rather than
      ! silently evaluating something else: the low-Mach contact-mode
      ! dissipation flux is not part of the generic text, the
      ! reconstruction flags must be the consistent pair the schemes are
      ! selected by, and a reconstruction continuation strictly between the
      ! two schemes is refused because the momentum row's terms are read
      ! off by a single-scheme rule (the_generic_text_is_usable).
      !
      ! WHAT IT HANDS BACK BESIDE THE ROWS.  The face flux, the face
      ! pressure and, under the well-balanced key, the face departures the
      ! momentum row was assembled from, all in the module arrays of
      ! RK_integration that every later reader of them uses
      ! (store_the_interface_fluxes).  That is what lets the momentum row's
      ! terms and its reference scale describe the state just
      ! assembled here and not the last state RK_rhs saw.

      use global_parameters
      use RK_integration, only: face_flux, face_p, face_q_up, face_q_dn
      use Reconstruction_step, only: n_faces_positivity_limited
      use Numerical_Fluxes, only: n_faces_roe_hlle, n_faces_llf
      use low_mach_dissipation, only: low_mach_damping_active
      use hydrodynamic_rows_double,                                      &
             only: rows_in_double => hydrodynamic_flux_difference_and_source
      use hydrodynamic_rows_quadruple,                                   &
             only: rows_in_quad   => hydrodynamic_flux_difference_and_source

      implicit none
      private
      public :: generic_precision_rows_selected
      public :: ROWS_PRODUCTION, ROWS_QUADRUPLE, ROWS_GENERIC_DOUBLE
      public :: hydrodynamic_rows_in_quadruple_precision
      public :: hydrodynamic_rows_in_double_precision
      public :: quad_rows_calls, quad_rows_seconds

      ! What EXHALE_RESID_QUAD selects.
      integer, parameter :: ROWS_PRODUCTION     = 0  ! the production routines
      integer, parameter :: ROWS_QUADRUPLE      = 1  ! the generic text at real(16)
      integer, parameter :: ROWS_GENERIC_DOUBLE = 2  ! the generic text at real(8)

      ! How many times the quadruple assembly ran and how long it spent
      ! there, so its cost is a measured number and not an estimate.
      integer :: quad_rows_calls   = 0
      real*8  :: quad_rows_seconds = 0.0d0

      contains

      ! ------------------------------------------------------!

      integer function generic_precision_rows_selected() result(mode)
      ! EXHALE_RESID_QUAD: absent or anything else ROWS_PRODUCTION, 1
      ! ROWS_QUADRUPLE, 2 ROWS_GENERIC_DOUBLE.  Read at every call rather than
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
      mode = ROWS_PRODUCTION
      env = ' '
      call get_environment_variable('EXHALE_RESID_QUAD', env,            &
                                    length=ln, status=st)
      if (st .ne. 0 .or. ln .le. 0) return
      select case (trim(adjustl(env)))
      case ('1')
         mode = ROWS_QUADRUPLE
      case ('2')
         mode = ROWS_GENERIC_DOUBLE
      end select
      end function generic_precision_rows_selected

      ! ------------------------------------------------------!

      subroutine hydrodynamic_rows_in_double_precision(u,WL,WR,dF,S)
      ! The generic text at kind(1.0d0).  Not on any run's path: it is the
      ! control of that experiment, compared against the production routines.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL,WR,dF,S
      real*8, dimension(3,1-Ng:N+Ng) :: ff
      real*8, dimension(1-Ng:N+Ng)   :: fp
      real*8, dimension(1-Ng:N+Ng)   :: fqu, fqd
      integer :: n_scaled, n_hlle, n_llf
      call the_generic_text_is_usable()
      call rows_in_double(u,WL,WR,dF,S,ff,fp,n_scaled,n_hlle,n_llf,     &
                          fq_up8=fqu, fq_dn8=fqd)
      call store_the_interface_fluxes(ff,fp,fqu,fqd,                     &
                                      n_scaled,n_hlle,n_llf)
      end subroutine hydrodynamic_rows_in_double_precision

      ! ------------------------------------------------------!

      subroutine hydrodynamic_rows_in_quadruple_precision(u,WL,WR,dF,S)
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL,WR,dF,S
      real*8, dimension(3,1-Ng:N+Ng) :: ff
      real*8, dimension(1-Ng:N+Ng)   :: fp
      real*8, dimension(1-Ng:N+Ng)   :: fqu, fqd
      integer :: n_scaled, n_hlle, n_llf
      real*8  :: t_a, t_b
      call the_generic_text_is_usable()
      call cpu_seconds(t_a)
      call rows_in_quad(u,WL,WR,dF,S,ff,fp,n_scaled,n_hlle,n_llf,       &
                        fq_up8=fqu, fq_dn8=fqd)
      call cpu_seconds(t_b)
      quad_rows_calls   = quad_rows_calls + 1
      quad_rows_seconds = quad_rows_seconds + (t_b - t_a)
      call store_the_interface_fluxes(ff,fp,fqu,fqd,                     &
                                      n_scaled,n_hlle,n_llf)
      end subroutine hydrodynamic_rows_in_quadruple_precision

      ! ------------------------------------------------------!

      subroutine store_the_interface_fluxes(ff,fp,fqu,fqd,               &
                                            n_scaled,n_hlle,n_llf)
      ! The face data the rows were differenced from, in the arrays every
      ! later reader of them uses (the species face flux, the flux gate,
      ! the jump scan, the momentum row's terms), and the three counters of
      ! the run summary.
      !
      ! THE FACE DEPARTURES BELONG TO THE STATE THAT WAS JUST ASSEMBLED.
      ! Under the well-balanced option the momentum row's pressure-gradient
      ! term is the gradient of the departure, (A+ q_up - A- q_dn)/dV under
      ! PLM and (q_up - q_dn)/dr under WENO3
      ! (momentum_row_terms_of_cell, RK_rhs.f90), so a reader handed the
      ! row without them would form that term from whatever face data was
      ! last left in the module: the previous evaluation's, or the zeros of
      ! the first allocation.  With the key off the row carries the cell's
      ! own pressure instead, nothing reads these two arrays, and they are
      ! left as they stand.
      !
      ! There is one pair of departures to store because this assembly evaluates
      ! ONE reconstruction: the_generic_text_is_usable refuses a
      ! reconstruction continuation with 0 < lambda < 1, where the blended
      ! row is not the row of any single pair.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: ff
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: fp
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: fqu, fqd
      integer, intent(in) :: n_scaled, n_hlle, n_llf
      if (.not. allocated(face_flux)) then
         allocate(face_flux(3,1-Ng:N+Ng), face_p(1-Ng:N+Ng))
      endif
      face_flux = ff
      face_p    = fp
      if (well_balanced) then
         if (.not. allocated(face_q_up)) then
            allocate(face_q_up(1-Ng:N+Ng), face_q_dn(1-Ng:N+Ng))
         endif
         face_q_up = fqu
         face_q_dn = fqd
      endif
      n_faces_positivity_limited = n_faces_positivity_limited + n_scaled
      n_faces_roe_hlle           = n_faces_roe_hlle + n_hlle
      n_faces_llf                = n_faces_llf + n_llf
      end subroutine store_the_interface_fluxes

      ! ------------------------------------------------------!

      subroutine the_generic_text_is_usable()
      ! The three configurations the generic text does not carry.  All of
      ! them stop the run: a control experiment that quietly evaluated a different
      ! operator would answer a question nobody asked.
      !
      ! A RECONSTRUCTION CONTINUATION STRICTLY INSIDE (0,1) IS ONE OF THEM.
      ! The generic text does blend the two schemes' rows there, as
      ! reconstruction_continuation_rhs does, but the momentum row's terms
      ! are then not recoverable: the face-pressure difference carries
      ! weight one under WENO3 and zero under PLM and
      ! momentum_row_terms_of_state reads that weight off the scheme flags,
      ! a single-scheme rule, while the two schemes give each interface
      ! different well-balanced departures so that no one pair satisfies
      ! the blended row's identity.  The endpoints lambda <= 0 and
      ! lambda >= 1 are a single scheme and are carried.
      if (recon_lambda_on) then
         if (recon_lambda .gt. 0.0d0 .and. recon_lambda .lt. 1.0d0) then
            write(*,*) 'ERROR: EXHALE_RESID_QUAD does not carry the',    &
                       ' PLM/WENO3 reconstruction continuation;',        &
                       ' recon_lambda = ', recon_lambda
            error stop 1
         endif
      endif
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

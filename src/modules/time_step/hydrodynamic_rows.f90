      module hydrodynamic_rows_double
      ! Double-precision finite-volume rows from hydrodynamic_rows_body.inc.
      ! This instantiation must match Reconstruction_step, Numerical_Fluxes,
      ! source_func, and RK_integration bitwise for the precision comparison.

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


      module hydrodynamic_rows_quadruple
      ! Quadruple-precision finite-volume arithmetic on double-precision
      ! composition, heating, and cooling. Near the hydrostatic base, Roe
      ! face jumps are 3e5-6e6 times smaller than their states, and r^2/dV
      ! reaches 5.1e3 (docs/Update_EXHALE_stage2.pdf). The wider mantissa
      ! lowers the measured residual roundoff floor by about 18 decades.

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


      module hydrodynamic_rows
      ! Stationary hydrodynamic rows use the production double-precision
      ! routines by default. EXHALE_RESID_QUAD=1 selects the same finite-volume
      ! reconstruction, Roe/HLLE flux, and sources in quadruple precision;
      ! rows, face fluxes, pressures, and sources are returned in double.
      !
      ! This diagnostic tests the measured 1.0e-11 to 1.6e-11 residual
      ! roundoff floor along a preconditioned Krylov direction, flat over 5.4
      ! decades of sampling spacing (docs/Update_EXHALE_stage2.pdf). Near
      ! the base, face-state cancellation reaches 3e5-6e6 and the row
      ! amplifies flux roundoff by r^2/dV. A well-balanced flux difference
      ! addresses that floor physically; extra precision diagnoses it.
      !
      ! The generic rows refuse low-Mach contact dissipation, inconsistent
      ! reconstruction flags, and 0 < recon_lambda < 1, whose blended
      ! momentum terms cannot be recovered from one pair of face departures.
      ! They store the current face flux, pressure, and well-balanced
      ! departures in RK_integration arrays for later row diagnostics.

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


      integer function generic_precision_rows_selected() result(mode)
      ! EXHALE_RESID_QUAD: 1 selects quadruple generic rows; 2 selects
      ! double generic rows for bitwise comparison with production; all
      ! other values select production. Read on each call so a test can
      ! change modes within one process.
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


      subroutine store_the_interface_fluxes(ff,fp,fqu,fqd,               &
                                            n_scaled,n_hlle,n_llf)
      ! Store the current face fluxes, pressure, and run counters.
      ! With well_balanced, store the departures used in the momentum row:
      ! (A+ q_up - A- q_dn)/dV for PLM, (q_up - q_dn)/dr for WENO3
      ! (momentum_row_terms_of_cell, RK_rhs.f90). Without it, the row uses
      ! cell pressure and the departure arrays are not read. A single pair
      ! suffices because the_generic_text_is_usable rejects 0 < lambda < 1.
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


      subroutine the_generic_text_is_usable()
      ! Reject configurations the generic rows cannot reproduce: low-Mach
      ! contact dissipation, inconsistent reconstruction flags, and
      ! 0 < recon_lambda < 1. In the blended case, the WENO3 face-pressure
      ! term has weight one and the PLM term zero, while
      ! momentum_row_terms_of_state reads a single scheme flag; one pair of
      ! well-balanced departures cannot reconstruct the blended row.
      ! The endpoints select a single scheme.
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

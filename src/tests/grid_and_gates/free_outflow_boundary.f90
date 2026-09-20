      program free_outflow_boundary
      ! WHAT THIS TEST IS ABOUT
      !   The outer boundary of these runs is a free outflow.  Two statements
      !   about it are checked here on a SUPERSONIC wind, and one on the ghost
      !   the boundary writes in any state.
      !
      !   1. A supersonic outflow carries no information inward THROUGH THE
      !      RIEMANN SOLUTION.  With every wave speed of the face positive,
      !      HLLC takes the branch SL >= 0 and returns the physical flux of
      !      the LEFT state alone (Num_Fluxes.f90; SL is clamped at zero by
      !      SL = min(0, ...), so the branch is entered with SL exactly 0).
      !      Replacing the right state of that face by any other admissible
      !      gas state must therefore leave the numerical flux bitwise
      !      unchanged.  Measured here by calling the production `Num_flux`
      !      twice at the same face with two different right states.
      !
      !   2. That statement is about the SOLVER and not about the stencil.
      !      The same perturbation at a SUBSONIC face changes the flux, and
      !      that half is asserted too: without it the first assertion could
      !      be passed by a solver that ignored its right state everywhere.
      !
      !   3. The ghost `Apply_BC` writes at an outflow must be an admissible
      !      gas state, must not reverse the flow (a ghost velocity of the
      !      opposite sign would make the boundary an INFLOW, which carries
      !      characteristics into the domain that only a reservoir may state,
      !      and this closure has none), and must leave the temperature of the
      !      last physical cell alone, since the ionization sweep and the
      !      column integral read the ghosts as cells, and `eval_dt` reads
      !      the outer ghost as one of the two states bounding the last face
      !      of cell N.
      !
      !   4. The outermost physical cell must carry a pressure slope of its
      !      own.  A zero-gradient outflow ghost makes the forward difference
      !      entering the MC limiter of `PLM_rec` exactly zero, so the limited
      !      slope of cell N is exactly zero and that cell reconstructs flat:
      !      it then carries no pressure gradient against its own weight, which
      !      is the O(1) hydrostatic imbalance measured by
      !      `hydrostatic_residual`.  Asserted as the two face states of cell
      !      N differing.
      !
      ! THE STATE
      !   A smooth isothermal radial outflow on the production grid, built to
      !   cross the sound speed inside the domain so that one face of it is
      !   subsonic and one is supersonic:
      !
      !     v(r)   = v_in + (v_out - v_in) (r - 1)/(r_max - 1)
      !     rho(r) = rho_1 v(1)/(v(r) r^2)          (steady continuity)
      !     p(r)   = nhat rho(r)                    (isothermal, T = T0)
      !
      !   with v_in = 0.2 c_s and v_out = 2.5 c_s, c_s = sqrt(gamma nhat) in
      !   code units.  Nothing here is a solution of the momentum equation and
      !   nothing needs to be: what is tested is the boundary closure and the
      !   Riemann solver's upwinding, both of which are statements about the
      !   local state.
      !
      !   Everything is PRODUCTION: define_grid, set_gravity_grid, Apply_BC,
      !   Reconstruct and Num_flux.  The driver builds the state and reads
      !   the results back.

      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use base_boundary,              only: set_base_reservoir
      use BC_Apply,                   only: Apply_BC
      use Reconstruction_step,        only: Reconstruct
      use Numerical_Fluxes,           only: Num_flux
      use Conversion,                 only: U_to_W, W_to_U_comp

      implicit none

      ! ---- the mechanical column of backup/regression/hydrostatic_column,
      !      READ from its input.inp: any bound planet serves, and this one
      !      is the geometry hydrostatic_residual already measures on ----
      real*8, parameter :: Rp_RJ   = 0.49d0     ! Planet radius [R_J]
      real*8, parameter :: Mp_MJ   = 0.0457d0   ! Planet mass [M_J]
      real*8, parameter :: Teq     = 1140.0d0   ! Equilibrium temperature [K]
      real*8, parameter :: rmax_c  = 3.0d0      ! Outer radius [R_p]
      real*8, parameter :: he_h    = 0.0793d0   ! He/H number ratio
      real*8, parameter :: m_He    = 3.9715d0   ! He mass in m_H

      integer, parameter :: n_scheme = 2
      character(len=8), dimension(n_scheme) :: scheme_name

      real*8  :: nhat, cs
      integer :: nfail, is

      scheme_name(1) = 'PLM'
      scheme_name(2) = 'WENO3'

      nhat  = (1.0d0 + he_h)/(1.0d0 + m_He*he_h)
      nfail = 0

      spherical_domain = .true.
      dr_base          = 2.0e-4
      N_low_cells      = 50
      r_max            = rmax_c
      r_esc            = 2.0d0
      r_flux           = 1.20d0
      grid_type        = 'Mixed'
      CFL              = 0.6d0
      N                = 500
      T0               = Teq
      R0               = Rp_RJ*RJ
      Mp               = Mp_MJ*MJ
      v0               = sqrt(kb_erg*T0/mu)
      b0               = (Gc*Mp*mu)/(kb_erg*T0*R0)
      cs               = sqrt(gamma_ad*nhat)

      call allocate_grid_arrays
      call define_grid
      call set_gravity_grid

      do is = 1,n_scheme
         call measure_one_scheme(is)
      enddo

      write(*,'(A)') ''
      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'free_outflow_boundary: ', nfail,            &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'free_outflow_boundary: all assertions PASSED'

      contains

      ! ------------------------------------------------------------------ !

      subroutine measure_one_scheme(is)
      integer, intent(in) :: is

      real*8, dimension(:,:), allocatable :: u, W, WL, WR
      real*8  :: v_in, v_out, rho1, vloc, Wpert(3)
      real*8  :: NF(3), NFp(3), p_out, p_outp, dNF, mach_N, mach_sub
      real*8  :: dslope, dtemp, tN, tg
      integer :: j, k, j_sub

      allocate(u(3,1-Ng:N+Ng), W(3,1-Ng:N+Ng), WL(3,1-Ng:N+Ng),           &
               WR(3,1-Ng:N+Ng))

      select case (is)
      case (1)
         rec_method = 'PLM';   use_plm = .true.;  use_weno3 = .false.
      case (2)
         rec_method = 'WENO3'; use_plm = .false.; use_weno3 = .true.
      end select
      flux = 'HLLC'

      ! ---- the wind ----
      v_in  = 0.2d0*cs
      v_out = 2.5d0*cs
      rho1  = 1.0d0
      do j = 1-Ng,N+Ng
         vloc   = v_in + (v_out - v_in)*(r(j) - 1.0d0)/(r_max - 1.0d0)
         W(1,j) = rho1*v_in/(vloc*r(j)*r(j))
         W(2,j) = vloc
         W(3,j) = nhat*W(1,j)
         call W_to_U_comp(W(:,j), u(:,j), j)
      enddo

      ! What the lower boundary states. The base is not the subject here; it
      ! is given the state of cell 1 so that it writes something admissible.
      n_part_cell1 = nhat*W(1,1)
      call set_base_reservoir(W(3,1), 1.0d0, nhat, 1.0d0)

      call Apply_BC(u)
      call U_to_W(u, W)
      call Reconstruct(u, WL, WR)

      ! ---- 3. the ghost the boundary wrote ----
      do k = 1,Ng
         call verdict_true('outflow_ghost_is_admissible['//               &
              trim(scheme_name(is))//']',                                 &
              W(1,N+k) .gt. 0.0d0 .and. W(3,N+k) .gt. 0.0d0, nfail)
         call verdict_true('outflow_ghost_is_not_an_inflow['//            &
              trim(scheme_name(is))//']',                                 &
              W(2,N+k)*W(2,N) .ge. 0.0d0, nfail)
      enddo
      tN = W(3,N)/W(1,N)
      tg = W(3,N+1)/W(1,N+1)
      dtemp = abs(tg - tN)/tN
      call verdict_below('outflow_ghost_keeps_the_temperature['//         &
           trim(scheme_name(is))//']', dtemp, 1.0d-14, nfail)

      ! ---- 4. the outermost cell's own slope ----
      ! WL(:,N) is cell N reconstructed to its outer face and WR(:,N-1) is
      ! the same cell reconstructed to its inner face; they differ by the
      ! cell's limited slope and are equal only if that slope is zero.
      dslope = abs(WL(3,N) - WR(3,N-1))/W(3,N)
      call verdict_above('outer_cell_carries_a_pressure_slope['//         &
           trim(scheme_name(is))//']', dslope, 1.0d-12, nfail)

      ! ---- 1. and 2. the Riemann solver's upwinding ----
      ! The perturbation of the right state: an admissible gas state well
      ! away from the one the reconstruction produced.
      mach_N = W(2,N)/sqrt(gamma_ad*W(3,N)/W(1,N))

      call Num_flux(WL(:,N), WR(:,N), NF, p_out, N, N+1)
      Wpert(1) = 3.0d0*WR(1,N)
      Wpert(2) = 0.5d0*WR(2,N)
      Wpert(3) = 3.0d0*WR(3,N)
      call Num_flux(WL(:,N), Wpert, NFp, p_outp, N, N+1)
      dNF = maxval(abs(NFp - NF))/maxval(abs(NF))
      write(*,'(A,A,A,F8.4,A,ES12.5)') '  DIAGNOSTIC ',                   &
         trim(scheme_name(is)), '  outer-face Mach = ', mach_N,           &
         '  relative flux change under the perturbed ghost = ', dNF
      call verdict_true('outer_face_is_supersonic['//                     &
           trim(scheme_name(is))//']', mach_N .gt. 1.0d0, nfail)
      ! BITWISE, not "small": the HLLC branch returns Phys_flux(WL) and never
      ! forms an expression in the right state at all.
      call verdict_true('supersonic_face_flux_ignores_the_ghost['//       &
           trim(scheme_name(is))//']', dNF .eq. 0.0d0, nfail)

      ! The same perturbation at the innermost SUBSONIC face, which must move
      ! the flux: otherwise the assertion above would hold for a solver that
      ! never read its right state.
      j_sub = 1
      do j = 1,N-1
         if (W(2,j)/sqrt(gamma_ad*W(3,j)/W(1,j)) .lt. 1.0d0) then
            j_sub = j
            exit
         endif
      enddo
      mach_sub = W(2,j_sub)/sqrt(gamma_ad*W(3,j_sub)/W(1,j_sub))
      call Num_flux(WL(:,j_sub), WR(:,j_sub), NF, p_out, j_sub, j_sub+1)
      Wpert(1) = 3.0d0*WR(1,j_sub)
      Wpert(2) = 0.5d0*WR(2,j_sub)
      Wpert(3) = 3.0d0*WR(3,j_sub)
      call Num_flux(WL(:,j_sub), Wpert, NFp, p_outp, j_sub, j_sub+1)
      dNF = maxval(abs(NFp - NF))/maxval(abs(NF))
      write(*,'(A,A,A,I5,A,F8.4)') '  DIAGNOSTIC ', trim(scheme_name(is)),&
         '  subsonic face j = ', j_sub, '  Mach = ', mach_sub
      call verdict_above('subsonic_face_flux_reads_the_right_state['//    &
           trim(scheme_name(is))//']', dNF, 1.0d-12, nfail)

      deallocate(u, W, WL, WR)

      end subroutine measure_one_scheme

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

      ! ------------------------------------------------------------------ !

      subroutine verdict_true(name, ok, nfail)
      ! A verdict with no number to report: the condition holds or it does
      ! not. 1 and 0 are printed so the line has the shape of every other.
      character(len=*), intent(in)    :: name
      logical,          intent(in)    :: ok
      integer,          intent(inout) :: nfail
      if (ok) then
         write(*,'(A,A,A)') 'PASS ', name,                                &
              ' measured= 1.00000E+00 reference= 1.00000E+00 tol=0'
      else
         write(*,'(A,A,A)') 'FAIL ', name,                                &
              ' measured= 0.00000E+00 reference= 1.00000E+00 tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_true

      end program free_outflow_boundary

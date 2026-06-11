      module steady_residual_mod
      ! Finite-volume STEADY residual  R = du/dt  (zero at a true steady
      ! state), shared by the ATES_RESIDUAL diagnostic, the in-loop
      ! convergence monitor, and the steady-state Newton/PTC solver.
      !
      !   R(:,1) = dF - S                 (mass)
      !   R(:,2) = dF - S                 (momentum)
      !   R(:,3) = dF_E - S_E - (heat-cool)   (energy: flux/gravity minus
      !                                        the radiative net source)
      !
      ! dF, S come from the existing Reconstruct + RK_rhs (HLLC; alpha is
      ! unused by the HLLC flux so 0 is passed). heat/cool are supplied by
      ! the caller, which decides how they were obtained: the in-loop
      ! monitor reuses the current step's values; the diagnostic and the
      ! Newton residual first call ioniz_eq (local ionization-equilibrium
      ! elimination) to get heat/cool consistent with u.
      !
      ! The reconstruction scheme is whatever use_plm/use_weno3/rec_method
      ! currently select; callers set WENO3 for a production residual.

      use global_parameters
      use Reconstruction_step
      use RK_integration

      implicit none
      private
      public :: assemble_residual, residual_norms

      contains

      ! ------------------------------------------------------!

      subroutine assemble_residual(u, heat, cool, R)
      real*8, dimension(1-Ng:N+Ng,3), intent(in)  :: u
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: heat, cool
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: R
      ! Local scratch so callers' own WL/WR/dF/S are untouched
      real*8, dimension(1-Ng:N+Ng,3) :: WL, WR, dF, S

      call Reconstruct(u, WL, WR)
      call RK_rhs(u, WL, WR, 0.0d0, dF, S)
      R(:,1) = dF(:,1) - S(:,1)
      R(:,2) = dF(:,2) - S(:,2)
      R(:,3) = dF(:,3) - S(:,3) - (heat - cool)

      end subroutine assemble_residual

      ! ------------------------------------------------------!

      subroutine residual_norms(R, u, rc)
      ! Per-component max relative residual rate over the wind [j_min:N]:
      !   rc(k) = max_j |R(j,k)| / max|u(:,k)|     [units 1/t_s]
      real*8, dimension(1-Ng:N+Ng,3), intent(in)  :: R
      real*8, dimension(1-Ng:N+Ng,3), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      integer :: k
      do k = 1,3
         rc(k) = maxval(abs(R(j_min:N,k)))                                &
                 /max(maxval(abs(u(j_min:N,k))), 1.0d-30)
      enddo
      end subroutine residual_norms

      ! End of module
      end module steady_residual_mod

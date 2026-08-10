      module steady_residual_mod
      ! Finite-volume STEADY residual  R = du/dt  (zero at a true steady
      ! state), shared by the EXHALE_RESIDUAL diagnostic, the in-loop
      ! convergence monitor, and the steady-state Newton/PTC solver.
      !
      !   R(:,1) = dF - S                          (mass)
      !   R(:,2) = dF - S - F_mu                   (momentum)
      !   R(:,3) = dF_E - S_E - (heat-cool)
      !                       - (w F_mu + q_mu + conduction)   (energy)
      !
      ! dF, S come from the existing Reconstruct + RK_rhs (HLLC; alpha is
      ! unused by the HLLC flux so 0 is passed). heat/cool are supplied by
      ! the caller, which decides how they were obtained: the in-loop
      ! monitor reuses the current step's values; the diagnostic and the
      ! Newton residual first call ioniz_eq (local ionization-equilibrium
      ! elimination) to get heat/cool consistent with u.
      !
      ! The molecular-transport terms are the operator-split stage that
      ! viscous_conduction_step relaxes in the marching loop, evaluated by
      ! the SAME routine (viscous_conduction_sources) from the same
      ! tridiagonal operator, so the Newton solver solves exactly the system
      ! the marching relaxes. They vanish identically unless "Viscosity:" /
      ! "Conduction:" are set. The caller passes n_part = n_tot + n_e (the
      ! adimensional particle count) rather than T itself, so that T = p/n_part
      ! stays a function of the unknowns and the temperature dependence of
      ! the conduction operator is picked up by the residual's linearization.
      !
      ! The reconstruction scheme is whatever use_plm/use_weno3/rec_method
      ! currently select; callers set WENO3 for a production residual.

      use global_parameters
      use Conversion, only: U_to_W
      use Reconstruction_step
      use RK_integration
      use viscous_conduction, only: transport_active,                   &
                                    viscous_conduction_sources

      implicit none
      private
      public :: assemble_residual, residual_norms, residual_norms_vol

      contains

      ! ------------------------------------------------------!

      subroutine assemble_residual(u, n_part, heat, cool, R)
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: n_part
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: R
      ! Local scratch so callers' own WL/WR/dF/S are untouched
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dF, S, W
      real*8, dimension(1-Ng:N+Ng)   :: Tc, Smom, Sene

      call Reconstruct(u, WL, WR)
      call RK_rhs(u, WL, WR, 0.0d0, dF, S)
      R(1,:) = dF(1,:) - S(1,:)
      R(2,:) = dF(2,:) - S(2,:)
      R(3,:) = dF(3,:) - S(3,:) - (heat - cool)

      if (transport_active()) then
         call U_to_W(u, W)
         Tc = W(3,:)/n_part
         call viscous_conduction_sources(W(2,:), Tc, Smom, Sene)
         R(2,:) = R(2,:) - Smom
         R(3,:) = R(3,:) - Sene
      endif

      end subroutine assemble_residual

      ! ------------------------------------------------------!

      subroutine residual_norms(R, u, rc)
      ! The max relative residual rate for each component over the wind [j_min:N]:
      !   rc(k) = max_j |R(j,k)| / max|u(:,k)|     [units 1/t_s]
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: R
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      integer :: k
      do k = 1,3
         rc(k) = maxval(abs(R(k,j_min:N)))                                &
                 /max(maxval(abs(u(k,j_min:N))), 1.0d-30)
      enddo
      end subroutine residual_norms

      ! ------------------------------------------------------!

      subroutine residual_norms_vol(Res, u, rc)
      ! VOLUME-WEIGHTED relative residual over the wind [j_min:N]:
      !   rc(k) = sum_j |R(j,k)| V_j / sum_j |u(j,k)| V_j,   V_j = r_j^2 dr_j
      ! (the spherical cell volume up to the 4*pi factor, which cancels).
      ! Physical meaning: the numerator is d/dt of the volume-integrated
      ! conserved quantity k, so rc(k) is the fractional drift rate of the
      ! GLOBAL mass/momentum/energy budget. Unlike the L-inf residual_norms
      ! (max over cells, ~1/dr so dominated by the smallest cells on the
      ! non-uniform grid), this weights each cell by its volume and is not
      ! biased by an isolated small near-base cell.
      ! NB: the residual argument is named Res (not R) because Fortran is
      ! case-insensitive and the grid radius array is r.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: Res, u
      real*8, dimension(3),           intent(out) :: rc
      integer :: k, j
      real*8  :: num, den, w
      do k = 1,3
         num = 0.0d0;  den = 0.0d0
         do j = j_min, N
            w   = r(j)*r(j)*dr_j(j)
            num = num + abs(Res(k,j))*w
            den = den + abs(u(k,j))*w
         enddo
         rc(k) = num/max(den, 1.0d-30)
      enddo
      end subroutine residual_norms_vol

      ! End of module
      end module steady_residual_mod

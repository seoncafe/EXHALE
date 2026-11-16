      module low_mach_dissipation
      ! Mach-gated numerical stress for the 2 dr contact/entropy mode left
      ! weakly damped by HLLC as v approaches zero. Jameson, Schmidt &
      ! Turkel (1981) motivate the fourth-difference momentum flux:
      ! D_p(j+1/2) = eps4*g*rho_f*lambda_f
      !              *(v(j+2)-3*v(j+1)+3*v(j)-v(j-1)),
      ! lambda_f = mean(|v|+c_s), D_E = v_f*D_p, and D_rho = 0.
      ! D_p and D_E vanish at faces 0 and N, conserving total momentum and
      ! total energy. The paired stress work transfers energy between kinetic
      ! and internal forms; its sign is not guaranteed where the face
      ! coefficient varies. The actual variable-coefficient form is
      ! sign-indefinite at the Mach gate; the nonpositive alternative
      ! -M^-1 B^T K B has not been adopted. Do not interpret this option as
      ! physical viscosity or as universally dissipative.
      ! g = max(0,1-M_f^2/M_th^2)^2, M_f^2 = mean(M_j^2,M_(j+1)^2).
      ! It is C1 at the threshold and zero for M_f >= M_th. On a resolved
      ! smooth profile the term vanishes as dr^3; the 2 dr mode decays at
      ! 16*eps4*g*lambda/dr. Explicit stability requires
      ! eps4 < 1/(16*CFL), checked by input_read.
      ! Valid only where the stress is negligible beside physical fluxes;
      ! report its size in any run that enables it. The stress and work term
      ! enter RK_rhs, so marching and the stationary residual use one flux.
      ! "Low-Mach damping: <eps4> [<M_th>]" is off by default (eps4 <= 0).
      ! Variable-coefficient energy analysis and measurements:
      ! ../EXHALE_v1.00/md/lhs1140b_stationary_L25_20260916.md section 1.
      ! THE SIGN OF THAT EXCHANGE IS INDEFINITE for the stress (1).  It is
      ! a dissipation only where the face coefficient
      !
      !   w_j = r_edg(j)^2 eps4 g_j rho_f,j lambda_f,j >= 0
      !
      ! is constant.  At frozen density, on the spherical weights RK_rhs
      ! forms the cell update with, the discrete kinetic-energy rate of (1)
      ! is the quadratic form
      !
      !   E' = sum_{j=1}^{N-1} w_j d_j (L d)_j = d^T W L d ,            (E)
      !   d_j = v_{j+1} - v_j ,   (L d)_j = d_{j+1} - 2 d_j + d_{j-1} ,
      !
      ! and it carries no boundary term, because D_p is set to zero at the
      ! faces j = 0 and j = N.  L is negative semidefinite, so E' <= 0 when
      ! W = w I; for a varying W the form is governed by the symmetric part
      ! of W L, which is not sign definite. The smooth gate (5) makes W
      ! vary across the Mach threshold. Where E' > 0
      ! the pair takes energy out of the internal energy and puts it into
      ! the velocity field.
      !
      ! It is written in M^2 = v^2/c_s^2 rather than in |v|/c_s because |v|
      ! has a kink at v = 0 and the stagnant layer is exactly where v changes
      ! sign; M^2 is a smooth function of the conserved variables everywhere.

      use global_parameters
      use caloric_eos, only: pressure_from_energy_density,             &
                            adiabatic_index_from_state

      implicit none
      private
      public :: low_mach_damping_active, contact_mode_dissipation_flux,  &
                contact_mode_dissipation_magnitude

      contains


      logical function low_mach_damping_active()
      low_mach_damping_active = (lowmach_damp_eps .gt. 0.0d0)
      end function low_mach_damping_active


      subroutine contact_mode_dissipation_flux(u, Dflux, gate_face)
      ! Face-centered dissipative flux, Eqs. (1) and (4), in code units.
      ! Dflux(:,j) belongs to the face r_edg(j) between cells j and j+1 and
      ! is nonzero only on the interior faces j = 1..N-1.  Component 1 (mass)
      ! is identically zero -- see section 2.  gate_face, if asked for, returns
      ! the gate (5) at the same faces, so that the one definition of the gate
      ! serves both the flux and the admissibility report below.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: Dflux
      real*8, dimension(1-Ng:N+Ng), intent(out), optional :: gate_face
      ! Cell-centered fields over the stencil range of the interior faces
      real*8, dimension(0:N+1) :: dens, vel, lam, mach2
      real*8  :: pres, cs2, gate, m2_face, lam_face, rho_face, v_face
      real*8  :: m2_th
      real*8, parameter :: cs2_floor = 1.0d-30    ! code units, (v/v0)^2
      integer :: j

      Dflux = 0.0d0
      if (present(gate_face)) gate_face = 0.0d0
      if (lowmach_damp_eps .le. 0.0d0) return
      m2_th = lowmach_damp_mach_th*lowmach_damp_mach_th

      ! Floor the squared sound speed before dividing. This does not reject
      ! nonpositive pressure: a zero-velocity cell still has Mach number zero.
      do j = 0, N+1
         dens(j) = max(u(1,j), tiny(1.0d0))
         vel(j)  = u(2,j)/dens(j)
         pres    = pressure_from_energy_density(j, dens(j),             &
                      u(3,j) - 0.5d0*u(2,j)*vel(j))
         cs2     = max(adiabatic_index_from_state(j, dens(j), pres)     &
                       *pres/dens(j), cs2_floor)
         lam(j)   = abs(vel(j)) + sqrt(cs2)
         mach2(j) = vel(j)*vel(j)/cs2
      enddo

      do j = 1, N-1
         m2_face = 0.5d0*(mach2(j) + mach2(j+1))
         gate    = 1.0d0 - m2_face/m2_th
         if (gate .le. 0.0d0) cycle
         gate     = gate*gate
         if (present(gate_face)) gate_face(j) = gate
         lam_face = 0.5d0*(lam(j)  + lam(j+1))
         rho_face = 0.5d0*(dens(j) + dens(j+1))
         v_face   = 0.5d0*(vel(j)  + vel(j+1))
         Dflux(2,j) = lowmach_damp_eps*gate*rho_face*lam_face             &
                    *(vel(j+2) - 3.0d0*vel(j+1)                           &
                               + 3.0d0*vel(j) - vel(j-1))
         Dflux(3,j) = v_face*Dflux(2,j)
      enddo

      end subroutine contact_mode_dissipation_flux


      subroutine contact_mode_dissipation_magnitude(u, ratio_max, r_peak,  &
                                                    r_gate_out)
      ! Size of the artificial stress (1) against the PHYSICAL momentum flux
      ! rho v^2 + p it is added to, at the same faces, plus the outermost
      ! radius at which the gate (5) is still open.  This is the
      ! admissibility check section 4 asks for: a numerical dissipation is
      ! admissible only where it is negligible against the physical fluxes,
      ! and neither the fourth difference nor the gate is a bound on that
      ! ratio by itself.  Returns ratio_max = 0 and r_peak = r_gate_out = -1
      ! when the term is off or the gate is closed everywhere.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, intent(out) :: ratio_max      ! max |D_p| / |rho v^2 + p|
      real*8, intent(out) :: r_peak          ! radius of that face [R_p]
      real*8, intent(out) :: r_gate_out     ! outermost open-gate face [R_p]
      real*8, dimension(3,1-Ng:N+Ng) :: Dflux
      real*8, dimension(1-Ng:N+Ng)   :: gate_face
      real*8 :: rho_f, v_f, p_f, phys, ratio
      integer :: j

      ratio_max = 0.0d0;  r_peak = -1.0d0;  r_gate_out = -1.0d0
      if (lowmach_damp_eps .le. 0.0d0) return

      call contact_mode_dissipation_flux(u, Dflux, gate_face)

      do j = 1, N-1
         if (gate_face(j) .le. 0.0d0) cycle
         r_gate_out = r_edg(j)
         rho_f = 0.5d0*(u(1,j) + u(1,j+1))
         v_f   = 0.5d0*(u(2,j)/max(u(1,j),  tiny(1.0d0))                 &
                      + u(2,j+1)/max(u(1,j+1),tiny(1.0d0)))
         p_f   = 0.5d0*( pressure_from_energy_density(j,                &
                            max(u(1,j),  tiny(1.0d0)),                  &
                            u(3,j)   - 0.5d0*u(2,j)**2                  &
                                       /max(u(1,j),  tiny(1.0d0)))      &
                       + pressure_from_energy_density(j+1,              &
                            max(u(1,j+1),tiny(1.0d0)),                  &
                            u(3,j+1) - 0.5d0*u(2,j+1)**2                &
                                       /max(u(1,j+1),tiny(1.0d0))) )
         phys  = abs(rho_f*v_f*v_f + p_f)
         if (phys .le. 0.0d0) cycle
         ratio = abs(Dflux(2,j))/phys
         if (ratio .gt. ratio_max) then
            ratio_max = ratio
            r_peak     = r_edg(j)
         endif
      enddo

      end subroutine contact_mode_dissipation_magnitude

      end module low_mach_dissipation

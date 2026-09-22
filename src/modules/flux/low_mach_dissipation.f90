      module low_mach_dissipation
      ! An artificial viscous stress that damps the 2 dr velocity mode a
      ! contact-resolving upwind flux leaves undamped where the flow stalls.
      !
      ! -------------------------------------------------------------------
      ! 1. What loses its damping, and where
      ! -------------------------------------------------------------------
      ! The HLLC flux resolves the middle (contact/entropy) wave exactly: the
      ! jump carried by that wave moves at the contact speed S* ~ v, and the
      ! dissipation the solver applies to it is proportional to |S*|.  It
      ! therefore VANISHES as the flow stagnates, v -> 0.  The acoustic
      ! families keep their damping at |v| +- c_s, but the entropy field and
      ! the velocity field that the near-hydrostatic force balance ties to it
      ! acquire a 2 dr pattern that nothing in the scheme removes.
      !
      ! In HD 209458 b with molecular chemistry AND solar C, N, O
      ! (examples/16_molecular_metals) the metal cooling exceeds the
      ! photoionization heating by a factor ~10 at r ~ 1.02 R_p and stops the
      ! flow: the mean |v| over r = 1.015-1.030 falls to 6-24 cm/s against
      ! ~207 cm/s in the metals-off twin, i.e. M ~ 3e-5 against a sound speed
      ! of ~4.5 km/s, and the alternating component of v grows to 4.7-10
      ! times the local mean |v| where the metals-off solution carries 0.001.
      ! Every failed steady solve of that configuration puts its worst scaled
      ! residual in that shell.
      !
      ! -------------------------------------------------------------------
      ! 2. The stress
      ! -------------------------------------------------------------------
      ! A fourth-difference (Jameson, Schmidt & Turkel 1981) stress is added
      ! to the numerical MOMENTUM flux at each interior face,
      !
      !   D_p{j+1/2} = eps4 g(M) rho_f lambda_f
      !                ( v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1} ) ,      (1)
      !
      ! with rho_f, lambda_f, v_f the arithmetic face averages and
      !
      !   lambda_f = 1/2 [ (|v| + c_s)_j + (|v| + c_s)_{j+1} ]          (2)
      !
      ! the face spectral radius -- the same signal speed the upwind flux
      ! uses, which in the stagnant layer is just c_s.  The undivided third
      ! difference in (1) is dr^3 d^3v/dr^3 + O(dr^4), so (1) is the
      ! hyper-viscous stress
      !
      !   tau_art = - mu_4 d^3 v/dr^3 ,   mu_4 = eps4 g rho c_s dr^3 ,   (3)
      !
      ! i.e. the acoustic momentum diffusivity rho c_s dr that an upwind
      ! scheme already carries for the acoustic waves, reapplied at fourth
      ! order to the family that lost it.
      !
      ! A stress does work, so the TOTAL energy flux carries with it
      !
      !   D_E{j+1/2} = v_f D_p{j+1/2} ,                                 (4)
      !
      ! exactly as viscous_conduction pairs w F_mu + q_mu with F_mu for the
      ! physical Navier-Stokes stress.  Without (4) the kinetic energy the
      ! stress removes would be taken out of the INTERNAL energy instead of
      ! being converted into it, which is a spurious heat source of
      ! indefinite sign.  With (4) the pair conserves total energy exactly,
      ! so the internal-energy change is exactly minus the kinetic-energy
      ! change.
      !
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
      ! of W L, which is not sign definite, and the gate (5) below makes W
      ! vary by a STEP wherever the Mach number crosses M_th.  Where E' > 0
      ! the pair takes energy out of the internal energy and puts it into
      ! the velocity field.
      !
      ! MEASURED (the smallest eigenvalue
      ! of the symmetric part of the operator of E' over its largest, against
      ! the eigensolver's own rounding): +7.2e-17 on a uniform grid at
      ! constant coefficient, -1.1e-12 for a coefficient falling four decades
      ! smoothly over 500 cells, -6.9e-03 with the gate closing at one face,
      ! and -1.1e-10 on an LHS 1140 b wind, where lambda_min = -4.147e+02
      ! against a rounding of 8.0e-02 and the negative mode sits on the cells
      ! astride the gate edge at r = 4.0-5.1 R_p.  The same lambda_min comes
      ! out of the interior block that no ghost value reaches, so the sign
      ! failure is interior and no boundary closure produces it.
      !
      ! THE FORM THAT IS DISSIPATIVE ON EVERY GRID writes the same stencil as
      ! a conservative flux of the CELL second difference,
      !
      !   q_j = 2 v_j - v_{j+1} - v_{j-1} ,
      !   A_j D_p{j+1/2} = k_j q_j - k_{j+1} q_{j+1} ,  k_1 = k_N = 0 ,  (G)
      !
      ! with k_j = eps4 g(M_j^2) rho_j lambda_j A^c_j >= 0 read at the cell,
      ! which is where M_j^2 is defined.  In matrix form (G) is
      ! dv/dt = -M^-1 B^T K B v with M = diag(rho_j dV_j) the mass matrix,
      ! (B v)_j = q_j and K = diag(k_j), so
      ! d(v^T M v/2)/dt = -(B v)^T K (B v) <= 0 for every velocity field and
      ! every nonuniform grid, and no ghost value enters (G) at all.  For
      ! constant k, (G) IS (1), so the order of accuracy, the 2 dr decay rate
      ! and the explicit-stability bound of sections 4 and 6 carry over
      ! unchanged.  MEASURED with the same driver: -4.5e-17 with the gate
      ! closing at one face and -5.9e-18 on the LHS 1140 b wind.
      !
      ! (G) IS NOT WHAT THIS MODULE APPLIES.  The executable form is (1) with
      ! the pair (4); (G) is not adopted here, so a run that turns this
      ! option on carries the indefinite form.
      !
      ! THE MASS FLUX IS NOT TOUCHED, and neither is the energy flux other
      ! than through (4).  MEASURED on this configuration, a fourth-difference
      ! flux on the conserved variables themselves -- the plain JST form --
      ! is inadmissible here: the physical mass flux rho v and energy flux
      ! v(E+p) both vanish with the velocity, while dr^3 d^3(rho)/dr^3 and
      ! dr^3 d^3(E)/dr^3 do not, because rho and E are stratified over a
      ! barely resolved scale height.  Their dissipative divergences come out
      ! 1e4-1e5 times the physical mass-flux divergence and ~5e3 times the
      ! radiative source, i.e. they would rewrite the base velocity and
      ! thermal structure rather than damp an oscillation.  The momentum flux
      ! rho v^2 + p does NOT vanish at stagnation -- it tends to p -- which is
      ! why (1) stays small there.
      !
      ! -------------------------------------------------------------------
      ! 3. The gate
      ! -------------------------------------------------------------------
      !   g = [ max(0, 1 - M_f^2/M_th^2) ]^2 ,
      !   M_f^2 = 1/2 (M_j^2 + M_{j+1}^2)                               (5)
      !
      ! is 1 at M = 0, exactly 0 for M_f >= M_th, and C^1 at M_f = M_th so
      ! that the steady residual stays differentiable for the Newton solve.
      ! It is written in M^2 = v^2/c_s^2 rather than in |v|/c_s because |v|
      ! has a kink at v = 0 and the stagnant layer is exactly where v changes
      ! sign; M^2 is a smooth function of the conserved variables everywhere.
      ! With M_th = 1e-3 the term is dead above r ~ 1.10 R_p in the
      ! HD 209458 b models -- far below the r_esc = 2 R_p escape radius
      ! (M ~ 0.10) over which convergence is measured and the r = N-20 cell
      ! at which Mdot is evaluated, so neither can be affected by it.
      !
      ! -------------------------------------------------------------------
      ! 4. Why a fourth difference, and why it is bounded
      ! -------------------------------------------------------------------
      ! Written as a flux difference, (1) contributes -eps4 g lambda
      ! (delta^4 v)_j/dr to dv_j/dt.  On the 2 dr mode v_j = (-1)^j a one has
      ! (delta^4 v)_j = 16 a, so the mode decays at 16 eps4 g lambda/dr; on a
      ! resolved profile the same expression is -eps4 g lambda dr^3 d^4v/dr^4,
      ! i.e. O(dr^3), far below the scheme's own truncation error.  A 1-2-1
      ! Shapiro filter, by contrast, is a SECOND difference: it damps the
      ! 2 dr mode at the same order at which it smooths everything else.
      !
      ! Where the pair (1)+(4) removes kinetic energy, what it can move is
      ! bounded by the KINETIC energy of the oscillation it removes: at Mach
      ! number M the kinetic energy density is g(g-1)/2 M^2 times the
      ! internal one, and it is consumed as the oscillation dies.  MEASURED
      ! over the region where the gate is open in the HD 209458 b molecular
      ! runs, M = 2e-5 to 5e-4, so that ratio is 2e-10 to 1e-7: whatever the
      ! pair moves that way, it cannot be a thermally significant amount of
      ! energy.  On the modes for which E' > 0 (section 2) the transfer runs
      ! the other way and this bound does not apply to it, which is why the
      ! size of the term is reported for any run that uses it (VALIDITY).
      !
      ! -------------------------------------------------------------------
      ! 5. Discrete consistency and conservation
      ! -------------------------------------------------------------------
      ! (1) and (4) are added to the numerical flux inside RK_rhs, the SAME
      ! routine the marching loop and assemble_residual (hence the JFNK
      ! steady solver) both call.  The equation the Newton residual measures
      ! is therefore the equation the marching loop relaxes, and the fixed
      ! point of one is the zero of the other.  This is the property the
      ! Shapiro filter lacks: it is applied to the marching state only, so
      ! the Newton residual never sees it, and the mode it suppresses during
      ! marching is still an undamped mode of the system Newton solves.
      !
      ! They are FLUXES, so momentum and total energy are conserved to
      ! round-off.  Both are set to zero at the base face (j = 0) and at the
      ! outer face (j = N), so nothing is injected or removed through the
      ! boundaries; those are also the faces whose four-cell stencil would
      ! reach outside the ghost layer.  Cell j then depends on cells
      ! j-2 .. j+2, which is the stencil width the WENO3 residual already
      ! has, so the banded Jacobian bandwidth (kl_jac = ku_jac = 8) is
      ! unchanged.
      !
      ! -------------------------------------------------------------------
      ! 6. Choosing eps4
      ! -------------------------------------------------------------------
      ! From section 4 the 2 dr mode decays at 16 eps4 lambda/dr, while the
      ! explicit marching step is dt = CFL dr/lambda, so the mode is damped
      ! by a factor 16 eps4 CFL per step and explicit stability requires
      !
      !   eps4 < 1/(16 CFL)                                             (6)
      !
      ! (0.10 at the default CFL = 0.6); input_read warns when (6) is
      ! violated.  Values around 0.01-0.03 -- the classical JST range 1/64
      ! to 1/32 -- damp the mode by 10-30 per cent per step.
      !
      ! VALIDITY.  This is a numerical dissipation, not a transport
      ! coefficient: it has no physical counterpart and is admissible only
      ! where it is negligible against the physical fluxes.  What makes it so
      ! is structural rather than tuned -- the fourth difference vanishes as
      ! dr^3 on resolved structure and the gate (5) vanishes identically
      ! outside the stagnant layer -- but the size of the term relative to
      ! the physical fluxes must be reported for any run that uses it.
      ! Physical damping of the same layer belongs to viscous_conduction
      ! (Navier-Stokes viscosity and conduction); MEASURED on this
      ! configuration that route made the steady solve worse, which is why a
      ! numerical dissipation is offered separately instead of being folded
      ! into it.
      !
      ! OFF BY DEFAULT ("Low-Mach damping: <eps4> [<M_th>]", eps4 <= 0
      ! disables), and the call site skips it entirely when off, so a run
      ! without the key is byte-identical to the code without this module.

      use global_parameters
      use caloric_eos, only: pressure_from_energy_density,             &
                            adiabatic_index_from_state

      implicit none
      private
      public :: low_mach_damping_active, contact_mode_dissipation_flux,  &
                contact_mode_dissipation_magnitude

      contains

      ! ------------------------------------------------------!

      logical function low_mach_damping_active()
      low_mach_damping_active = (lowmach_damp_eps .gt. 0.0d0)
      end function low_mach_damping_active

      ! ------------------------------------------------------!

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

      ! cs2_floor keeps the Mach number finite in a cell whose reconstructed
      ! pressure has gone non-positive: the gate then reads M^2 = huge, which
      ! switches the term off in that cell rather than dividing by zero.
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

      ! ------------------------------------------------------!

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

      ! End of module
      end module low_mach_dissipation

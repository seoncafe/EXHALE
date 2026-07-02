      module lower_column
      ! Tier-1 analytic lower/middle atmosphere (Koskinen et al. 2022, Sec 2.2).
      !
      ! Integrates the hypsometric relation from the 1-bar level (radius
      ! R_1bar, e.g. the broadband transit radius) up to the escape-model base
      ! pressure (default 1 microbar), isothermal at T_eq, with the mean
      ! molecular weight mu(p,T) following the chemical-equilibrium H2/H/He
      ! partition of Visscher et al. (2006) as fitted in Koskinen et al. (2022,
      ! Eqs. 11-13):
      !
      !   q_H2 = [1.9845 + 10^u - sqrt(10^u (3.9690 + 10^u))]/2.3670,
      !   u    = -23672/T - log10(p[bar]) + 6.2645,
      !
      ! with q_He from the fixed elemental He/H and q_H = 1 - q_H2 - q_He.
      ! (Validated against NASA CEA at 1000-2500 K in that paper.)
      !
      ! Deliverables: the radius of the base-pressure level r0 (the quantity
      ! the user currently guesses as "Planet radius"), and the H2/H/He volume
      ! mixing ratios + mean molecular weight at the base.  NOTE: in CHEMICAL
      ! EQUILIBRIUM q_H2(1 ubar) stays large up to T ~ 2000 K (fully atomic
      ! only above ~2400 K), so for Teq ~ 1000-2000 K hot Jupiters the fit
      ! reports a molecular base even though photochemistry + the (hotter)
      ! actual base temperature dissociate H2 there (Koskinen 2013a; Moses
      ! 2011).  The solver therefore also returns r_base_atomic (fully atomic
      ! mu) to BRACKET the base radius; a large q_H2 on a genuinely cool
      ! planet signals that Tier-2 molecular physics is required.
      !
      ! Approximations (documented): isothermal column at T_eq; spherical
      ! gravity g = GM/r^2 (the Roche/tidal correction to the column is a
      ! refinement -- Koskinen 2022 use equipotentials, which matters near
      ! RLO); chemical equilibrium only (photochemistry UNDERESTIMATES H
      ! dissociation at Teff ~ 1000-2000 K -- the Tier-3 hook).
      !
      ! Validation gate (docs/lower_atmosphere_coupling.*): Koskinen 2022
      ! Model A -- Uranus-mass planet (M = 8.6813e28 g, R_1bar = 2.5559e9 cm)
      ! around a Sun-like star at 0.05 au (T_eq = 1140 K, A_B = 0.3):
      ! r0(1 ubar) = 1.34 R_1bar with q_H2(base) ~ 0.84.

      implicit none
      private
      public :: lower_column_solve, q_h2_equilibrium

      real*8, parameter :: kbol   = 1.380649d-16   ! erg/K
      real*8, parameter :: mh_g   = 1.6726d-24     ! H mass [g]
      real*8, parameter :: Gcgs   = 6.67259d-8     ! CGS gravitational constant

      contains

      ! ------------------------------------------------------------------ !

      ! Chemical-equilibrium H2 volume mixing ratio q_H2(p,T) for an H/He gas
      ! (Koskinen et al. 2022 Eq. 11, from Visscher et al. 2006).
      ! p in bar, T in K.  Clamped to [0,1] for numerical safety far outside
      ! the fit range (very cold -> 1, very hot -> 0, as physically expected).
      double precision function q_h2_equilibrium(p_bar, T) result(qh2)
      real*8, intent(in) :: p_bar, T
      real*8 :: u, tenu
      u = -23672.0d0/T - log10(max(p_bar, 1.0d-30)) + 6.2645d0
      if (u .gt. 30.0d0) then
         qh2 = 1.0d0                                          ! -> fully H2
      else if (u .lt. -30.0d0) then
         qh2 = 0.0d0                                          ! fully atomic
      else
         tenu = 10.0d0**u
         qh2  = (1.9845d0 + tenu - sqrt(tenu*(3.9690d0 + tenu)))/2.3670d0
         if (qh2 .lt. 0.0d0) qh2 = 0.0d0
         if (qh2 .gt. 1.0d0) qh2 = 1.0d0
      endif
      end function q_h2_equilibrium

      ! ------------------------------------------------------------------ !

      ! Mean molecular weight [amu] of the H2/H/He mixture at (p,T), given the
      ! ELEMENTAL helium abundance fhe = n_He/n_H(nuclei).  Per unit H nucleus:
      ! H2 binds a fraction x2 of H nuclei into x2/2 molecules, so
      !   particles = (1 - x2) + x2/2 + fhe,   mass = 1 + 4 fhe   [amu].
      double precision function mu_mixture(p_bar, T, fhe) result(mu)
      real*8, intent(in) :: p_bar, T, fhe
      real*8 :: q, x2, npart
      ! The Visscher/Koskinen fit returns the MIXTURE volume mixing ratio
      ! q_H2 = n_H2/(n_H2+n_H+n_He) (protosolar He baked into the fit).
      ! Invert for the H-nucleus fraction bound in H2:
      !   x2 = 2 q (1+fhe) / (1+q)   (per H nucleus: n_H2=x2/2, n_H=1-x2,
      !   n_He=fhe; verified against Koskinen 2022 Model A base composition
      !   q_H2=0.84, q_H=0.026 -> x2=0.985).
      q  = q_h2_equilibrium(p_bar, T)
      x2 = 2.0d0*q*(1.0d0 + fhe)/(1.0d0 + q)
      if (x2 .gt. 1.0d0) x2 = 1.0d0
      npart = (1.0d0 - x2) + 0.5d0*x2 + fhe
      mu = (1.0d0 + 4.0d0*fhe)/npart
      end function mu_mixture

      ! ------------------------------------------------------------------ !

      ! Integrate the isothermal hypsometric relation
      !   d ln p / dr = - mu(p,T) m_H g(r) / (k T),   g = G M / r^2,
      ! from (p = p_deep at r = R_deep) up to p = p_base, with 4th-order RK in
      ! x = ln p (adaptive-free fixed step; the integrand is smooth).
      !
      ! Inputs (cgs):  Mp [g], R_deep [cm] (radius of the p_deep level),
      !                T [K] (isothermal), fhe (elemental He/H),
      !                p_deep, p_base [bar]
      ! Outputs:       r_base [cm], q_H2/q_H/q_He (volume fractions) and
      !                mu [amu] at the base.
      ! r_base_atomic brackets the equilibrium answer from the other side:
      ! the same column integrated with a FULLY ATOMIC mu = (1+4fhe)/(1+fhe)
      ! (chemical equilibrium is known to UNDERESTIMATE H dissociation at
      ! Teff ~ 1000-2000 K because photochemistry is ignored; Koskinen 2022).
      subroutine lower_column_solve(Mp, R_deep, T, fhe, p_deep, p_base,   &
                                    r_base, qh2_b, qh_b, qhe_b, mu_b,     &
                                    r_base_atomic)
      real*8, intent(in)  :: Mp, R_deep, T, fhe, p_deep, p_base
      real*8, intent(out) :: r_base, qh2_b, qh_b, qhe_b, mu_b
      real*8, intent(out) :: r_base_atomic

      integer, parameter :: nstep = 2000
      real*8 :: x, dx, r, k1, k2, k3, k4
      real*8 :: q, x2, ntot_per_H
      integer :: i

      ! integrate dr/dx = -kT / (mu(p) m_H g(r)),  x = ln p  (x decreasing)
      x  = log(p_deep)
      dx = (log(p_base) - log(p_deep))/dble(nstep)
      r  = R_deep
      do i = 1, nstep
         k1 = drdx(r,              exp(x))
         k2 = drdx(r + 0.5d0*dx*k1, exp(x + 0.5d0*dx))
         k3 = drdx(r + 0.5d0*dx*k2, exp(x + 0.5d0*dx))
         k4 = drdx(r + dx*k3,       exp(x + dx))
         r  = r + dx*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
         x  = x + dx
      enddo
      r_base = r

      ! atomic bracket: same integration with mu fixed at the atomic value
      x  = log(p_deep)
      r  = R_deep
      do i = 1, nstep
         k1 = drdx_atomic(r)
         k2 = drdx_atomic(r + 0.5d0*dx*k1)
         k3 = drdx_atomic(r + 0.5d0*dx*k2)
         k4 = drdx_atomic(r + dx*k3)
         r  = r + dx*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
         x  = x + dx
      enddo
      r_base_atomic = r

      ! base composition (volume fractions of the full H2/H/He mixture)
      q  = q_h2_equilibrium(p_base, T)
      x2 = 2.0d0*q*(1.0d0 + fhe)/(1.0d0 + q) ! H nuclei bound in H2
      if (x2 .gt. 1.0d0) x2 = 1.0d0
      ntot_per_H = (1.0d0 - x2) + 0.5d0*x2 + fhe
      qh2_b = 0.5d0*x2/ntot_per_H
      qh_b  = (1.0d0 - x2)/ntot_per_H
      qhe_b = fhe/ntot_per_H
      mu_b  = (1.0d0 + 4.0d0*fhe)/ntot_per_H

      contains

         double precision function drdx(rr, p_bar)
         real*8, intent(in) :: rr, p_bar
         real*8 :: g, mu
         g    = Gcgs*Mp/(rr*rr)
         mu   = mu_mixture(p_bar, T, fhe)
         drdx = -kbol*T/(mu*mh_g*g)
         end function drdx

         double precision function drdx_atomic(rr)
         real*8, intent(in) :: rr
         real*8 :: g, mu
         g    = Gcgs*Mp/(rr*rr)
         mu   = (1.0d0 + 4.0d0*fhe)/(1.0d0 + fhe)
         drdx_atomic = -kbol*T/(mu*mh_g*g)
         end function drdx_atomic

      end subroutine lower_column_solve

      ! End of module
      end module lower_column

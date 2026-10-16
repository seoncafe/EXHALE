      module lower_column
      ! analytic lower/middle atmosphere (Koskinen et al. 2022, Sec 2.2).
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
      ! planet signals that molecular physics is required.
      !
      ! Approximations (documented): isothermal column at T_eq; spherical
      ! gravity g = GM/r^2 (the Roche/tidal correction to the column is a
      ! refinement -- Koskinen 2022 use equipotentials, which matters near
      ! RLO); chemical equilibrium only (photochemistry UNDERESTIMATES H
      ! dissociation at Teff ~ 1000-2000 K -- the external-handoff hook);
      ! mu is the H2/H/He mean molecular weight only (the trace-metal mass is
      ! neglected in the column mu, a << 1% effect at solar metallicity).
      !
      ! Validation gate (docs/lower_atmosphere_coupling.*): Koskinen 2022
      ! Model A -- Uranus-mass planet (M = 8.6813e28 g, R_1bar = 2.5559e9 cm)
      ! around a Sun-like star at 0.05 au (T_eq = 1140 K, A_B = 0.3):
      ! r0(1 ubar) = 1.34 R_1bar with q_H2(base) ~ 0.84.

      ! Every physical constant this column uses is the one of
      ! global_parameters, renamed only where the local name would collide
      ! with a variable of this module (mu is the mean molecular weight
      ! here, so the hydrogen atom mass keeps the name mh_g).  In
      ! particular m_He_over_m_H, the weight the column's mean molecular
      ! weight gives one helium nucleus, is the ratio of the two measured
      ! atom masses formed there, so this column and the species table
      ! cannot weigh a helium atom differently.
      ! The hydrogen ATOM mass is the right mass here because it multiplies
      ! a dimensionless mean molecular weight, so it has to be the mass the
      ! rest of the code normalizes to and not the proton mass, which is
      ! 5.6e-4 lighter and would shorten the scale height by as much.
      use global_parameters, only: m_He_over_m_H,                       &
                                   kbol => kb_erg,                      &
                                   mh_g => mu,                          &
                                   Gcgs => Gc

      implicit none
      private
      public :: lower_column_solve, q_h2_equilibrium
      ! Exposed so an acceptance test reads the helium weight this column
      ! uses, and not a second copy of it.
      public :: m_He_over_m_H

      contains

      ! ------------------------------------------------------------------ !

      ! Chemical-equilibrium H2 volume mixing ratio q_H2(p,T) for an H/He gas
      ! (Koskinen et al. 2022, ApJ 929, 52, Eq. 11, quoting Visscher et al.
      ! 2006; read from the published paper). p in bar, T in K.
      !
      ! THE LIMITS OF THE FIT ITSELF, which the expression must reach on its
      ! own. u = -23672/T - log10(p) + 6.2645 runs to -infinity as the gas
      ! gets COLD and to +infinity as it gets HOT, and the expression
      !     q = (1.9845 + 10^u - sqrt(10^u (3.9690 + 10^u)))/2.3670
      ! therefore tends to
      !     cold, u -> -inf :  1.9845/2.3670 = 0.83840   (fully molecular)
      !     hot,  u -> +inf :  0                         (fully atomic)
      ! The cold limit is the fit's approximation of the fully molecular
      ! solar-composition value 0.5/(0.5 + He/H) = 0.863 at He/H = 0.0793;
      ! it is a volume mixing ratio n_H2/(n_H2 + n_H + n_He), so it cannot
      ! reach 1 in a gas that contains helium.
      !
      ! VALIDITY IN He/H: THIS FIT IS FOR SOLAR COMPOSITION AND MUST NOT BE
      ! USED IN A GAS RICHER IN HELIUM THAN ABOUT He/H = 0.1.
      ! The cold asymptote 1.9845/2.3670 = 0.83840 is a CONSTANT, while the
      ! attainable ceiling 0.5/(0.5 + He/H) falls as helium is added. They
      ! cross at
      !     He/H = 0.5/0.83840 - 0.5 = 0.0964,
      ! and above that the fit returns a mixing ratio the element ratio
      ! cannot supply. Solar composition clears the ceiling by only 3%
      ! (0.8384 against 0.8631); at He/H = 1 the fit returns 0.8384 against a
      ! ceiling of 0.3333, 2.5x over, and at He/H = 10 it is 17.6x over. The
      ! whole of the He/H = 1-1000 diagnostic ladder is outside this fit.
      ! A molecular base there has to state q_H2 explicitly through a
      ! lower-atmosphere handoff; input_read refuses the fit's value rather
      ! than capping it (section 117 of docs/Update_EXHALE.*, which is where
      ! the silent cap that used to hide this was removed).
      !
      ! AN ASYMPTOTIC BRANCH ON EITHER SIDE IS A DISCONTINUITY AT ITS OWN
      ! THRESHOLD AND THERE IS NONE HERE. A cold cut at u < -30 (T <~ 529 K at
      ! 1e-6 bar, the reachable side) would hand the gas the fully ATOMIC
      ! value where the fit gives 0.8384, a step of 0.84 in mixing ratio, and
      ! it would hand it to the two places that choose which basin a cold cell
      ! is solved from: the seed of the constrained continuation solve and the
      ! molecular-basin retry of ioniz_eq. A hot cut at u > 30 (p < 1e-24 bar)
      ! would return 0 where the fit gives 8.3e-31.
      !
      ! The cold side is left to the expression, which reaches the limit on
      ! its own: 10^u underflows to zero far below any threshold and the
      ! formula returns 1.9845/2.3670 exactly.
      !
      ! THE FORM THE FIT IS EVALUATED IN. Written as the fit states it,
      !     q = (A + t - sqrt(t (2A + t)))/B,   A = 1.9845, B = 2.3670,
      ! the two differenced terms both grow like t while their difference
      ! falls like A^2/(2t), so the difference drops below ulp(t) = t eps at
      !     t = A/sqrt(2 eps) = 9.4e7,
      ! and above that the value returned is the rounding of t, quantized in
      ! units of ulp(t)/B. On the LHS 1140 b column (T 418 to 5838 K,
      ! p 1.6e-15 to 9.5e-7 bar) that form is wrong by more than 1 percent in
      ! 186 of 500 cells, returns exactly zero in 164 of them where the fit is
      ! positive, and is 1.3e4 times the fit's value at its worst (MEASURED).
      !
      ! Multiplying by the conjugate,
      !     A + t - sqrt(t (2A + t)) = A^2/((A + t) + sqrt(t (2A + t)))
      ! exactly, and the right-hand side sums two positive terms instead of
      ! differencing two large ones. It is the same function: the cold limit
      ! is A/B = 0.83840304182509... to the last bit, it agrees with the form
      ! above wherever that form is not cancelling, and it carries the fit's
      ! own A^2/(2Bt) hot tail down to underflow. No asymptotic branch is
      ! needed on either side. The square root is taken as sqrt(t) sqrt(2A+t)
      ! so that the product under it cannot overflow at the largest t the
      ! exponent allows, and min(u, 300) keeps t itself finite; q there is
      ! 6.6e-301 and the fit is long since zero to any use it is put to.
      ! The form is positive and falls monotonically from A/B = 0.8384, so
      ! the ceiling below is a statement of the range and never fires, and no
      ! floor at 0 is needed at all. The ceiling the ELEMENT ratio imposes,
      ! 0.5/(0.5 + He/H), is a different statement and is the callers'.
      double precision function q_h2_equilibrium(p_bar, T) result(qh2)
      real*8, intent(in) :: p_bar, T
      real*8 :: u, tenu
      u    = -23672.0d0/T - log10(max(p_bar, 1.0d-30)) + 6.2645d0
      tenu = 10.0d0**min(u, 300.0d0)
      qh2  = 1.9845d0**2                                                 &
           / (2.3670d0*((1.9845d0 + tenu)                                &
              + sqrt(tenu)*sqrt(3.9690d0 + tenu)))
      if (qh2 .gt. 1.0d0) qh2 = 1.0d0
      end function q_h2_equilibrium

      ! ------------------------------------------------------------------ !

      ! Mean molecular weight of the H2/H/He mixture at (p,T), in units of the
      ! hydrogen ATOM (mh_g above, the unit it is multiplied by here and the
      ! one parameters.f90 normalizes to), NOT of the atomic mass unit u,
      ! given the ELEMENTAL helium abundance fhe = n_He/n_H(nuclei).
      ! Per unit H nucleus:
      ! H2 binds a fraction x2 of H nuclei into x2/2 molecules, so
      !   particles = (1 - x2) + x2/2 + fhe,
      !   mass      = 1 + (m_He/m_H) fhe   [hydrogen atoms].
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
      mu = (1.0d0 + m_He_over_m_H*fhe)/npart
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
      !                mu [hydrogen atoms] at the base.
      ! r_base_atomic brackets the equilibrium answer from the other side:
      ! the same column integrated with a FULLY ATOMIC
      ! mu = (1 + (m_He/m_H) fhe)/(1+fhe)
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
      mu_b  = (1.0d0 + m_He_over_m_H*fhe)/ntot_per_H

      contains

         double precision function drdx(rr, p_bar)
         real*8, intent(in) :: rr, p_bar
         real*8 :: grav_accel, mu_mean
         grav_accel    = Gcgs*Mp/(rr*rr)
         mu_mean   = mu_mixture(p_bar, T, fhe)
         drdx = -kbol*T/(mu_mean*mh_g*grav_accel)
         end function drdx

         double precision function drdx_atomic(rr)
         real*8, intent(in) :: rr
         real*8 :: grav_accel, mu_mean
         grav_accel    = Gcgs*Mp/(rr*rr)
         mu_mean   = (1.0d0 + m_He_over_m_H*fhe)/(1.0d0 + fhe)
         drdx_atomic = -kbol*T/(mu_mean*mh_g*grav_accel)
         end function drdx_atomic

      end subroutine lower_column_solve

      ! End of module
      end module lower_column

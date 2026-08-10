      module viscous_conduction
      ! Molecular transport in the radial (1-D, spherically symmetric) bulk
      ! flow: the Navier-Stokes viscous momentum force, its dissipation, and
      ! thermal conduction.  This is the Navier-Stokes level that CETIMB
      ! (Koskinen et al. 2013a, Icarus 226, 1678; Koskinen et al. 2022, ApJ
      ! 929, 52, Appendix B) solves and that EXHALE's inviscid HLLC scheme
      ! lacks; it supplies the physical damping of the near-base momentum
      ! imbalance ("breathing base").
      !
      ! -------------------------------------------------------------------
      ! 1. Equations
      ! -------------------------------------------------------------------
      ! Koskinen et al. (2022) write the bulk momentum and energy equations
      ! as their Eqs. (B2) and (B3),
      !
      !   d(rho w)/dt + (1/r^2) d/dr(r^2 rho w^2 + p) = -rho dU/dr + F_mu
      !   d(rho u)/dt + (1/r^2) d/dr(r^2 rho u w)
      !        = rho q - p (1/r^2) d/dr(r^2 w)
      !          + (1/r^2) d/dr(r^2 kappa dT/dr) + q_mu ,
      !
      ! with u = c_v T, and give the viscous momentum force F_mu and the
      ! dissipation functional q_mu as their Eqs. (B5) and (B6).  As PRINTED
      ! those two read
      !
      !   F_mu = (4/3)(1/r^2) d/dr(r^2 mu dw/dr) - (dmu/dr)(dw/dr)
      !          - (16/3) mu w/r^2                                    (B5)
      !   q_mu = (10/3) mu (dw/dr)^2 - (8/3) mu (w/r)(dw/dr)
      !          + (16/3) mu (w/r)^2                                  (B6)
      !
      ! and neither expression is the Navier-Stokes result for a radial flow.
      ! For the Newtonian stress
      !   tau_ij = mu (d_i w_j + d_j w_i) - (2/3) mu delta_ij (div w)
      ! in spherical symmetry the only non-zero components are
      !   tau_rr = (4/3) mu (w' - w/r),   tau_tt = tau_pp = -(1/2) tau_rr ,
      ! so that, exactly,
      !
      !   F_mu = (div tau)_r = (1/r^2) d/dr(r^2 tau_rr) + tau_rr/r
      !        = (4/3)(1/r^2) d/dr(r^2 mu dw/dr) - (4/3)(dmu/dr)(w/r)
      !          - (8/3) mu w/r^2                                      (1)
      !   q_mu = tau : grad w = (4/3) mu (w' - w/r)^2 .                (2)
      !
      ! (1) and (B5) agree in the leading constant-mu diffusion terms
      ! mu w'' and (8/3) mu w'/r but differ in the mu-gradient term
      ! ((4/3) mu' w/r versus mu' w') and in the geometric term (8/3 versus
      ! 16/3).  (2) and (B6) share the middle term but not the other two;
      ! (2) is manifestly non-negative, as a dissipation must be, and is the
      ! exact partner of (1): w F_mu + q_mu = div(tau . w), so the pair
      ! conserves total energy.  The printed (B5)/(B6) pair does not satisfy
      ! that identity.  The same appendix prints the pressure term inside the
      ! spherical divergence in (B2), which likewise cannot be meant
      ! literally.  We therefore implement (1) and (2) -- physical
      ! correctness is the acceptance criterion -- and record the
      ! discrepancy in docs/viscosity_conduction.md.
      !
      ! EXHALE evolves the TOTAL energy E = rho w^2/2 + p/(g-1), so the
      ! energy source that belongs with the momentum source F_mu is
      !
      !   S_E = w F_mu + q_mu + (1/r^2) d/dr(r^2 kappa dT/dr) .        (3)
      !
      ! -------------------------------------------------------------------
      ! 2. Transport coefficients
      ! -------------------------------------------------------------------
      ! Neither Koskinen et al. (2013a) nor (2022) prints mu(T) or kappa(T);
      ! O'Neill & Chorlton (1989) is cited there only as the textbook source
      ! of the dissipation functional.  We therefore use the atomic-hydrogen
      ! heat conduction coefficient of Watson, Donahue & Walker (1981), the
      ! standard choice in this literature (quoted e.g. as Eq. (A6) of Erkaev
      ! et al. 2016, MNRAS 460, 1300):
      !
      !   kappa(T) = 4.45e4 (T/1000 K)^0.7   erg cm^-1 s^-1 K^-1       (4)
      !
      ! and fix mu(T) to it through the Chapman-Enskog/Eucken relation for a
      ! monatomic gas, kappa = (15/4)(k_B/m) mu, i.e. Prandtl number 2/3:
      !
      !   mu(T) = (4/15)(m_H/k_B) kappa(T)
      !         = 1.44e-4 (T/1000 K)^0.7     g cm^-1 s^-1 .            (5)
      !
      ! VALIDITY.  (4) and (5) are the NEUTRAL atomic-hydrogen values.  Above
      ! the ionization front the electron (Spitzer) conductivity, ~ T^{5/2},
      ! is much larger, and Coulomb collisions raise the ion viscosity; that
      ! regime is NOT covered here.  The terms are included for the dense,
      ! largely neutral base region (r <~ 1.1 R_p), which is where the
      ! momentum imbalance they are meant to damp lives; in the tenuous
      ! ionized wind they are numerically negligible either way (Koskinen et
      ! al. 2013a note that conduction and viscosity are not important in the
      ! thermosphere of HD 209458 b).
      !
      ! -------------------------------------------------------------------
      ! 3. Code units
      ! -------------------------------------------------------------------
      ! Lengths in R0, velocities in v0 = sqrt(k T0/m_H), densities in
      ! rho0 = n0 m_H, temperatures in T0, so a force density scales with
      ! rho0 v0^2/R0 and an energy-density rate with q0 = rho0 v0^3/R0.
      ! Substituting into (1)-(3) leaves the expressions form-invariant with
      !
      !   mu_code    = mu_cgs   /(rho0 v0 R0)          (= 1/Reynolds)
      !   kappa_code = kappa_cgs T0/(rho0 v0^3 R0) = (15/4) mu_code .
      !
      ! -------------------------------------------------------------------
      ! 4. Discretization and boundary conditions
      ! -------------------------------------------------------------------
      ! Finite volume on the existing grid, with the same face areas and cell
      ! volume RK_rhs uses: A_j+1/2 = r_edg(j)^2, dV_j = (r_edg(j)^3 -
      ! r_edg(j-1)^3)/3.  Both operators are assembled ONCE, as tridiagonal
      ! coefficient triplets, and the SAME triplets are used (a) to evaluate
      ! the source that enters the steady residual and (b) as the matrix of
      ! the Crank-Nicolson update in the marching loop.  The marching scheme
      ! and the Newton residual therefore cannot disagree about the operator.
      !
      ! Inner boundary: the flux through the base face is formed from the
      ! ghost cell j = 0, i.e. from the state Apply_BC anchors there
      ! (rho_bc, the valved/mass-flux base velocity, the BC pressure).  The
      ! base cell then feels a viscous stress and a heat flux relative to
      ! that anchored state -- this is the whole point of the term: the
      ! breathing base cell can now dissipate against its anchor instead of
      ! only exchanging momentum with the cell above it.  During an implicit
      ! step the ghost is held fixed (Dirichlet), and its contribution moves
      ! to the right-hand side.
      !
      ! Outer boundary: ZERO diffusive flux through the last face.  The upper
      ! ghosts are a zero-gradient copy (PLM) or a linear extrapolation
      ! (WENO3) of the interior; the extrapolated ghost carries a spurious
      ! gradient that would inject an unphysical viscous stress and heat flux
      ! into the outflow.  A vanishing diffusive flux is the correct
      ! condition for a supersonic outflow boundary.
      !
      ! -------------------------------------------------------------------
      ! Runtime switches (both default OFF, so a run without the keys is
      ! byte-identical to the inviscid code):
      !   "Viscosity: True"     -> visc_on, calibrated mu(T) of (5) + q_mu
      !   "Viscosity: <mu0> [<s>]" -> visc_mu0 > 0, diagnostic power-law
      !                            override mu = visc_mu0 T^visc_s in CODE
      !                            units (unchanged legacy meaning)
      !   "Conduction: True"    -> cond_on, calibrated kappa(T) of (4)

      use global_parameters
      use Conversion, only: U_to_W

      implicit none
      private
      public :: viscosity_active, conduction_active, transport_active
      public :: dynamic_viscosity, thermal_conductivity
      public :: viscous_momentum_source, viscous_dissipation
      public :: thermal_conduction_source, viscous_conduction_sources
      public :: viscous_conduction_step

      ! Watson, Donahue & Walker (1981) atomic-hydrogen heat conduction,
      ! kappa = kappa_1000K (T/1000 K)^kappa_expo  [erg cm^-1 s^-1 K^-1].
      real*8, parameter :: kappa_1000K = 4.45d4
      real*8, parameter :: kappa_expo  = 0.7d0
      ! Chapman-Enskog/Eucken ratio for a monatomic gas: kappa = eucken k_B mu/m.
      real*8, parameter :: eucken      = 3.75d0        ! = 15/4

      contains

      ! ------------------------------------------------------!

      logical function viscosity_active()
      viscosity_active = (visc_on .or. visc_mu0 .gt. 0.0d0)
      end function viscosity_active

      ! ------------------------------------------------------!

      logical function conduction_active()
      conduction_active = cond_on
      end function conduction_active

      ! ------------------------------------------------------!

      logical function transport_active()
      transport_active = (viscosity_active() .or. conduction_active())
      end function transport_active

      ! ------------------------------------------------------!

      subroutine dynamic_viscosity(Tcell, visc)
      ! Dynamic viscosity in CODE units, Eq. (5) (or the diagnostic
      ! power-law override, which is already expressed in code units).
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: visc
      real*8 :: cgs_to_code, pref
      integer :: j
      if (.not. viscosity_active()) then
         visc = 0.0d0
         return
      endif
      if (visc_mu0 .gt. 0.0d0) then
         do j = 1-Ng, N+Ng
            visc(j) = visc_mu0*max(Tcell(j), 1.0d-8)**visc_s
         enddo
         return
      endif
      ! mu_cgs = (4/15)(m_H/k_B) kappa_cgs ; then to code units.
      cgs_to_code = 1.0d0/(n0*mu*v0*R0)
      pref = (mu/(eucken*kb_erg))*kappa_1000K*cgs_to_code                &
             *(T0/1.0d3)**kappa_expo
      do j = 1-Ng, N+Ng
         visc(j) = pref*max(Tcell(j), 1.0d-8)**kappa_expo
      enddo
      end subroutine dynamic_viscosity

      ! ------------------------------------------------------!

      subroutine thermal_conductivity(Tcell, kap)
      ! Heat conduction coefficient in CODE units, Eq. (4).
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: kap
      real*8 :: pref
      integer :: j
      if (.not. conduction_active()) then
         kap = 0.0d0
         return
      endif
      ! kappa_code = kappa_cgs T0/(rho0 v0^3 R0)
      pref = kappa_1000K*(T0/1.0d3)**kappa_expo*T0/(n0*mu*v0**3*R0)
      do j = 1-Ng, N+Ng
         kap(j) = pref*max(Tcell(j), 1.0d-8)**kappa_expo
      enddo
      end subroutine thermal_conductivity

      ! ------------------------------------------------------!

      subroutine viscous_momentum_coeffs(Tcell, alo, adi, aup, dcl, dcd, dcu)
      ! Tridiagonal coefficients of the radial viscous force (1),
      !   F_mu(j) = alo(j) w(j-1) + adi(j) w(j) + aup(j) w(j+1) ,
      ! and of the strain rate D = dw/dr - w/r that (2) squares,
      !   D(j)    = dcl(j) w(j-1) + dcd(j) w(j) + dcu(j) w(j+1) .
      ! j runs over the physical cells 1..N; w(0) is the base ghost.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(N),         intent(out) :: alo, adi, aup
      real*8, dimension(N),         intent(out) :: dcl, dcd, dcu
      real*8, dimension(1-Ng:N+Ng) :: visc
      real*8 :: Ap, Am, dV, rp, rm, vp, vm, wp, wm, dp, dm, c43
      real*8 :: pdi, pup, mlo, mdi, mup
      integer :: j

      call dynamic_viscosity(Tcell, visc)
      c43 = 4.0d0/3.0d0
      alo = 0.0d0;  adi = 0.0d0;  aup = 0.0d0
      dcl = 0.0d0;  dcd = 0.0d0;  dcu = 0.0d0

      do j = 1, N
         rp = r_edg(j);    rm = r_edg(j-1)
         Ap = rp*rp;       Am = rm*rm
         dV = (Ap*rp - Am*rm)/3.0d0

         ! ---- strain rate D = dw/dr - w/r at the cell center ----
         ! (one-sided in the outermost cell, where w(N+1) is an
         !  extrapolated ghost we do not want to lean on)
         if (j .lt. N) then
            dm = 1.0d0/(r(j+1) - r(j-1))
            dcl(j) = -dm
            dcu(j) =  dm
            dcd(j) = -1.0d0/r(j)
         else
            dm = 1.0d0/(r(j) - r(j-1))
            dcl(j) = -dm
            dcd(j) =  dm - 1.0d0/r(j)
            dcu(j) =  0.0d0
         endif

         ! ---- face stresses tau_rr = (4/3) mu (dw/dr - w/r) ----
         ! right face (zero diffusive flux at the outer boundary)
         pdi = 0.0d0;  pup = 0.0d0
         if (j .lt. N) then
            dp = 1.0d0/(r(j+1) - r(j))
            wp = (rp - r(j))*dp                     ! interpolation weight
            vp = visc(j) + wp*(visc(j+1) - visc(j))
            pdi = c43*vp*(-dp - (1.0d0 - wp)/rp)
            pup = c43*vp*( dp - wp/rp)
         endif
         ! left face (base face at j = 1 uses the anchored ghost)
         dm = 1.0d0/(r(j) - r(j-1))
         wm = (rm - r(j-1))*dm
         vm = visc(j-1) + wm*(visc(j) - visc(j-1))
         mlo = c43*vm*(-dm - (1.0d0 - wm)/rm)
         mdi = c43*vm*( dm - wm/rm)
         mup = 0.0d0

         ! ---- F_mu = (A_p tau_p - A_m tau_m)/dV + tau_c/r ----
         alo(j) = (       - Am*mlo)/dV + c43*visc(j)*dcl(j)/r(j)
         adi(j) = (Ap*pdi - Am*mdi)/dV + c43*visc(j)*dcd(j)/r(j)
         aup(j) = (Ap*pup - Am*mup)/dV + c43*visc(j)*dcu(j)/r(j)
      enddo

      end subroutine viscous_momentum_coeffs

      ! ------------------------------------------------------!

      subroutine thermal_conduction_coeffs(Tcell, blo, bdi, bup)
      ! Tridiagonal coefficients of (1/r^2) d/dr(r^2 kappa dT/dr):
      !   Q(j) = blo(j) T(j-1) + bdi(j) T(j) + bup(j) T(j+1) .
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(N),         intent(out) :: blo, bdi, bup
      real*8, dimension(1-Ng:N+Ng) :: kap
      real*8 :: Ap, Am, dV, rp, rm, kp, km, wp, wm, dp, dm
      integer :: j

      call thermal_conductivity(Tcell, kap)
      blo = 0.0d0;  bdi = 0.0d0;  bup = 0.0d0

      do j = 1, N
         rp = r_edg(j);    rm = r_edg(j-1)
         Ap = rp*rp;       Am = rm*rm
         dV = (Ap*rp - Am*rm)/3.0d0
         ! right face (zero heat flux at the outer boundary)
         if (j .lt. N) then
            dp = 1.0d0/(r(j+1) - r(j))
            wp = (rp - r(j))*dp
            kp = kap(j) + wp*(kap(j+1) - kap(j))
            bdi(j) = bdi(j) - Ap*kp*dp/dV
            bup(j) = bup(j) + Ap*kp*dp/dV
         endif
         ! left face (base face at j = 1 uses the anchored ghost)
         dm = 1.0d0/(r(j) - r(j-1))
         wm = (rm - r(j-1))*dm
         km = kap(j-1) + wm*(kap(j) - kap(j-1))
         bdi(j) = bdi(j) - Am*km*dm/dV
         blo(j) = blo(j) + Am*km*dm/dV
      enddo

      end subroutine thermal_conduction_coeffs

      ! ------------------------------------------------------!

      subroutine viscous_momentum_source(vel, Tcell, Fv)
      ! Radial viscous force density (1), code units, on cells 1..N.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: vel, Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Fv
      real*8, dimension(N) :: alo, adi, aup, dcl, dcd, dcu
      integer :: j
      Fv = 0.0d0
      if (.not. viscosity_active()) return
      call viscous_momentum_coeffs(Tcell, alo, adi, aup, dcl, dcd, dcu)
      do j = 1, N
         Fv(j) = alo(j)*vel(j-1) + adi(j)*vel(j) + aup(j)*vel(j+1)
      enddo
      end subroutine viscous_momentum_source

      ! ------------------------------------------------------!

      subroutine viscous_dissipation(vel, Tcell, qv)
      ! Viscous dissipation (2), q_mu = (4/3) mu (dw/dr - w/r)^2 >= 0.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: vel, Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: qv
      real*8, dimension(N) :: alo, adi, aup, dcl, dcd, dcu
      real*8, dimension(1-Ng:N+Ng) :: visc
      real*8  :: strain
      integer :: j
      qv = 0.0d0
      if (.not. viscosity_active()) return
      call viscous_momentum_coeffs(Tcell, alo, adi, aup, dcl, dcd, dcu)
      call dynamic_viscosity(Tcell, visc)
      do j = 1, N
         strain = dcl(j)*vel(j-1) + dcd(j)*vel(j) + dcu(j)*vel(j+1)
         qv(j)  = (4.0d0/3.0d0)*visc(j)*strain*strain
      enddo
      end subroutine viscous_dissipation

      ! ------------------------------------------------------!

      subroutine thermal_conduction_source(Tcell, Qc)
      ! Heat conduction source (1/r^2) d/dr(r^2 kappa dT/dr), code units.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Qc
      real*8, dimension(N) :: blo, bdi, bup
      integer :: j
      Qc = 0.0d0
      if (.not. conduction_active()) return
      call thermal_conduction_coeffs(Tcell, blo, bdi, bup)
      do j = 1, N
         Qc(j) = blo(j)*Tcell(j-1) + bdi(j)*Tcell(j) + bup(j)*Tcell(j+1)
      enddo
      end subroutine thermal_conduction_source

      ! ------------------------------------------------------!

      subroutine viscous_conduction_sources(vel, Tcell, Smom, Sene)
      ! The pair that enters the conserved-variable equations:
      !   Smom = F_mu                          (momentum, Eq. 1)
      !   Sene = w F_mu + q_mu + conduction    (TOTAL energy, Eq. 3)
      ! This is the single definition used by BOTH the steady residual and
      ! the marching update, so the two cannot describe different systems.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: vel, Tcell
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Smom, Sene
      real*8, dimension(1-Ng:N+Ng) :: qv, Qc
      Smom = 0.0d0;  Sene = 0.0d0
      if (.not. transport_active()) return
      call viscous_momentum_source(vel, Tcell, Smom)
      call viscous_dissipation(vel, Tcell, qv)
      call thermal_conduction_source(Tcell, Qc)
      Sene = vel*Smom + qv + Qc
      end subroutine viscous_conduction_sources

      ! ------------------------------------------------------!

      subroutine solve_tridiagonal(dlo, ddi, dup, rhs, sol)
      ! Thomas algorithm for the N x N tridiagonal system
      !   dlo(j) x(j-1) + ddi(j) x(j) + dup(j) x(j+1) = rhs(j).
      ! Radial recursion: inherently serial in j, so this must stay OUTSIDE
      ! any cell-parallel region.
      real*8, dimension(N), intent(in)  :: dlo, ddi, dup, rhs
      real*8, dimension(N), intent(out) :: sol
      real*8, dimension(N) :: cp, dp
      real*8  :: den
      integer :: j
      den   = ddi(1)
      cp(1) = dup(1)/den
      dp(1) = rhs(1)/den
      do j = 2, N
         den   = ddi(j) - dlo(j)*cp(j-1)
         cp(j) = dup(j)/den
         dp(j) = (rhs(j) - dlo(j)*dp(j-1))/den
      enddo
      sol(N) = dp(N)
      do j = N-1, 1, -1
         sol(j) = dp(j) - cp(j)*sol(j+1)
      enddo
      end subroutine solve_tridiagonal

      ! ------------------------------------------------------!

      subroutine viscous_conduction_step(u, W, Tcell, n_part, dt)
      ! Crank-Nicolson update of the two transport operators, applied as an
      ! operator-split stage of the marching loop, as CETIMB integrates the
      ! same terms.  Both are diffusive, so an explicit update would be bound
      ! by dt < rho dr^2/mu (and C dr^2/kappa); the implicit form removes that
      ! bound unconditionally.  MEASURED at the WASP-121 b resolution, though,
      ! the diffusive limit is ~1 t_s against an advective step of ~9e-5 t_s,
      ! so these terms are NOT stiff there and the implicit treatment is
      ! insurance for a denser base, a coarser grid or a larger coefficient
      ! rather than a necessity (docs/viscosity_conduction.md Sec. 4).
      !
      !   (i)  momentum   rho (w* - w)/dt = (1/2)[F_mu(w) + F_mu(w*)]
      !        The total energy is advanced by the viscous WORK that goes
      !        with it, so the internal energy is untouched by this stage:
      !        d(rho w^2/2) = w_avg rho (w* - w) = dt w_avg F_mu_avg exactly,
      !        which is what holding p fixed and rebuilding E accomplishes.
      !   (ii) temperature  C (T* - T)/dt = (1/2)[Q(T) + Q(T*)] + q_mu,
      !        C = n_part/(g-1) the internal energy per unit temperature,
      !        q_mu evaluated at the updated velocity.
      !
      ! The spatial operators are the SAME tridiagonal triplets that
      ! viscous_conduction_sources evaluates for the steady residual, so the
      ! fixed point of this update is exactly the zero of that residual.
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u, W
      real*8, dimension(1-Ng:N+Ng),   intent(in)    :: Tcell, n_part, dt
      real*8, dimension(1-Ng:N+Ng)   :: vnew, Tnew, qv
      real*8, dimension(N) :: alo, adi, aup, dcl, dcd, dcu
      real*8, dimension(N) :: blo, bdi, bup
      real*8, dimension(N) :: dlo, ddi, dup, rhs, sol
      real*8  :: half, Lold, cap
      integer :: j

      if (.not. transport_active()) return
      half = 0.5d0
      vnew = W(2,:)
      Tnew = Tcell

      ! ---- (i) viscous momentum diffusion ----
      if (viscosity_active()) then
         call viscous_momentum_coeffs(Tcell, alo, adi, aup, dcl, dcd, dcu)
         do j = 1, N
            Lold   = alo(j)*W(2,j-1) + adi(j)*W(2,j) + aup(j)*W(2,j+1)
            cap    = W(1,j)/dt(j)
            dlo(j) = -half*alo(j)
            ddi(j) = cap - half*adi(j)
            dup(j) = -half*aup(j)
            rhs(j) = cap*W(2,j) + half*Lold
         enddo
         ! The base ghost is Dirichlet during the step: its column moves to
         ! the right-hand side.  There is no upper-ghost column (zero
         ! diffusive flux at the outer face, and the outermost strain rate
         ! is one-sided), so nothing is dropped at j = N.
         rhs(1) = rhs(1) - dlo(1)*W(2,0)
         dlo(1) = 0.0d0
         dup(N) = 0.0d0
         call solve_tridiagonal(dlo, ddi, dup, rhs, sol)
         do j = 1, N
            vnew(j) = sol(j)
         enddo
         ! Rebuild the conserved state at FIXED pressure: the kinetic energy
         ! change is precisely the viscous work, so no internal energy is
         ! created or destroyed by the momentum stage.
         do j = 1, N
            W(2,j) = vnew(j)
            u(2,j) = W(1,j)*vnew(j)
            u(3,j) = 0.5d0*W(1,j)*vnew(j)**2 + W(3,j)/(g - 1.0d0)
         enddo
      endif

      ! ---- (ii) heat conduction + viscous dissipation ----
      if (conduction_active() .or. viscosity_active()) then
         call viscous_dissipation(W(2,:), Tcell, qv)
         if (conduction_active()) then
            call thermal_conduction_coeffs(Tcell, blo, bdi, bup)
         else
            blo = 0.0d0;  bdi = 0.0d0;  bup = 0.0d0
         endif
         do j = 1, N
            Lold   = blo(j)*Tcell(j-1) + bdi(j)*Tcell(j) + bup(j)*Tcell(j+1)
            cap    = n_part(j)/((g - 1.0d0)*dt(j))
            dlo(j) = -half*blo(j)
            ddi(j) = cap - half*bdi(j)
            dup(j) = -half*bup(j)
            rhs(j) = cap*Tcell(j) + half*Lold + qv(j)
         enddo
         rhs(1) = rhs(1) - dlo(1)*Tcell(0)
         dlo(1) = 0.0d0
         dup(N) = 0.0d0
         call solve_tridiagonal(dlo, ddi, dup, rhs, sol)
         do j = 1, N
            Tnew(j) = max(sol(j), 1.0d-2)
            W(3,j)  = n_part(j)*Tnew(j)
            u(3,j)  = 0.5d0*W(1,j)*W(2,j)**2 + W(3,j)/(g - 1.0d0)
         enddo
      endif

      end subroutine viscous_conduction_step

      ! End of module
      end module viscous_conduction

      module binary_element_diffusion
      ! Binary (two-component) H/He element transport for the single-fluid
      ! EXHALE wind.  Formulation: docs/binary_diffusion_design.md (sections
      ! 2-4); this module implements its atomic-region core (milestone M2).
      !
      ! COMPONENTS.  The gas is treated as two components moving through each
      ! other with a single mass-averaged velocity v (the hydro's velocity):
      !   component 1 : hydrogen carriers together with the trace metals
      !                 slaved to them at fixed metal/H.  Mass per H nucleus
      !                 m_1 = mass_per_H_nucleus_without_He() [m_H units];
      !   component He: helium, all stages.
      ! With that split rho_1 + rho_He = rho exactly (the metal mass is in rho
      ! under the same eos_metals policy that puts it into m_1), so the two
      ! diffusive mass fluxes close, J_1 = -J_He, and the composition update
      ! creates or destroys no mass.
      !
      ! STATE.  The transported variable is the helium MASS FRACTION
      !
      !    X = rho_He/rho = m_He n_He / (m_1 n_H + m_He n_He),   0 <= X <= 1
      !
      ! with n_H, n_He the ELEMENT (nucleus) densities counted over every
      ! species through the bsp_nH / bsp_nHe weights of species_table, so the
      ! count is the same quantity in the atomic, triplet-helium and molecular
      ! regions.  X is bounded at both ends of the composition axis, unlike the
      ! ratio f = n_He/n_H of the earlier trace kernel, which diverges exactly
      ! in the helium-rich limit this module exists for.
      !
      ! TRANSPORT.  Equation (1) of the memo,
      !
      !    d(rho X)/dt + (1/r^2) d/dr [ r^2 ( rho X v + J ) ] = 0
      !
      ! solved here in the ADVECTIVE form (memo eq. 5)
      !
      !    dX/dt + v dX/dr = -(1/(rho r^2)) d/dr [ r^2 J ]
      !
      ! because this routine receives only the new rho, v, T and dt -- not the
      ! hydro's Runge-Kutta face mass fluxes -- so an independent conservative
      ! advection of rho X could not keep a uniform X uniform in a compressing
      ! or expanding flow.  In the advective form a uniform X is preserved
      ! exactly whatever the hydro does (test T6), the composition and the
      ! conservative form agree exactly at a steady state (where rho v r^2 is
      ! constant), and in a transient the helium mass is conserved only to the
      ! order of the hydro's own truncation.  The advective term is a ONE-SIDED
      ! UPWIND difference taken with the CELL velocity v_j,
      !
      !    rho_j v_j (X_j - X_{j-1})/(r_j - r_{j-1})    for v_j >= 0
      !    rho_j v_j (X_{j+1} - X_j)/(r_{j+1} - r_j)    for v_j <  0
      !
      ! so a uniform X is exact, the row sum of the advective part is zero and
      ! the donor neighbour is the only off-diagonal it creates.  It is NOT
      ! built from face-averaged velocities: with v_f = (v_j + v_{j+1})/2 the
      ! upwind selection can pick the OUTWARD neighbour at both faces of a cell
      ! whose two face velocities straddle zero (v_f(j-1) < 0 < v_f(j)), and the
      ! advective term of that cell then vanishes identically even though its
      ! own v_j is large.  The cell is left with nothing but the molecular
      ! diffusion time dr^2/D_12 -- ~10^6 s at the base against a ~1 s hydro
      ! step -- so whatever composition it holds is frozen there.  That is what
      ! produced the isolated helium hole in the first free cell above the base
      ! of the breathing HD 209458 b wind, where the base sound wave alternates
      ! the sign of v from cell to cell.
      !
      ! DIFFUSIVE MASS FLUX (memo eqs. 2-4, 6).  For a binary mixture only one
      ! diffusive flux is independent; the helium one is
      !
      !    J = -rho (D_12 + K_zz) dX/dr  -  rho D_12 G X (1-X)
      !    G = [ (m_He - m_1) m_H g - (Zbar_He - Zbar_1) eE ] / (k T)
      !        + alpha_T dlnT/dr                                      [1/cm]
      !
      ! and the hydrogen-component flux is its negative.  The gradient
      ! coefficient of the memo, A = (m_1 m_He n/mbar) D_12 (dx/dX), collapses:
      ! with x = n_He/(n_1+n_He) one has x = m_1 X/(m_He - (m_He-m_1)X), hence
      ! dx/dX = mbar^2/(m_1 m_He) with mbar = rho/n, so A = n mbar D_12 = rho
      ! D_12.  Nothing divides by a vanishing element density anywhere -- that
      ! division is what made the trace kernel's self-consistent coupling
      ! diverge.  The settling prefactor rho X (1-X) switches the settling off
      ! as either element is exhausted, while the gradient flux does not vanish
      ! there, so a helium-free cell next to a helium-bearing one still
      ! receives helium.
      !
      ! Sign check (memo 2.3): neutral gas at uniform composition has G > 0 and
      ! therefore J < 0 -- helium drifts inward, i.e. settles.
      !
      ! AMBIPOLAR FIELD.  eE is computed, not assumed:
      !
      !    eE = -(1/n_e) d(n_e k T)/dr = -k T dln(n_e T)/dr        (memo 3a)
      !
      ! (Koskinen et al. 2013, section 2.1), evaluated from the solved electron
      ! density and temperature of the cell.  The logarithmic form is used
      ! because it is exact for the barometric stratification of the base while
      ! being the same second-order central difference elsewhere.  Limits: a
      ! neutral gas gives Zbar_He - Zbar_1 = 0 (relative settling mass 3 m_H),
      ! an H+ plasma eE = m_H g/2 (2.5), a fully ionized He++ plasma
      ! eE = 4/3 m_H g (5/3).  he_ambipolar = .false. sets eE = 0.
      !
      ! DIFFUSION COEFFICIENT.  D_12 = 1.52e18 (1/m_H + 1/m_He)^(1/2)
      ! T^(1/2) / n, n = n_H + n_He nuclei (Banks & Kockarts 1973).  This is a
      ! NEUTRAL hard-sphere binary coefficient, symmetric in the two species;
      ! the Coulomb correction for the ionized wind is NOT applied.  Validity:
      ! it is the right coefficient in the weakly ionized base and lower wind,
      ! and it sets WHERE separation turns on rather than the form of the
      ! transport equation in the ionized wind (memo, decision D3).  The metal
      ! nuclei are left out of n and of the reduced mass as a trace
      ! approximation (their number fraction is ~1e-4 at solar abundance).
      !
      ! DISCRETIZATION (memo section 3).  Finite volume on the existing grid,
      ! faces carrying r^2 areas, one implicit (backward-Euler, tridiagonal)
      ! solve per call with the coefficients frozen within the step and one
      ! Picard sweep on the X-dependent settling prefactor.  The settling term
      ! keeps the Peclet-based central/upwind hybrid: central where
      ! |W| dr < 2 (A + E), upwind toward the settling direction otherwise.
      !
      ! M-MATRIX CONDITION.  Row j of the backward-Euler matrix has
      !   off-diagonals   aa = -K r^2 Lf(j-1) <= 0,  cc = K r^2 Rf(j) <= 0
      !   diagonal        bb = rho/dt + (outflow face terms) > 0
      ! because (i) the gradient coefficient A + E >= 0 contributes +A/dr to
      ! Lf and -A/dr to Rf, (ii) the settling hybrid never lets its
      ! contribution flip the sign of an off-diagonal (that is exactly the
      ! Peclet switch: in the central branch |W|/2 <= A/dr), and (iii) the
      ! upwind advection contributes +rho|v|/dr to the diagonal and the same
      ! amount, negated, to the donor neighbour and to nothing else.  The row sum is
      ! rho/dt > 0, so the matrix is a strictly diagonally dominant M-matrix:
      ! its inverse is nonnegative and X^new >= 0 follows from X^old >= 0 and a
      ! nonnegative boundary value.  The same operator applied to 1-X (the
      ! hydrogen equation, which is this equation with J -> -J and hence the
      ! settling direction reversed) gives 1 - X^new >= 0 for a lagged
      ! prefactor in [0,1].  Hence the bounds hold per step WITHOUT clipping;
      ! the round-off clip below is an assertion, not a limiter.  There is no
      ! cap f <= HeH and no base pile-up limiter: a pile-up, if one appears, is
      ! physics for the base boundary condition to answer.
      !
      ! Gated on he_diffusion (default .false.); when off the module is never
      ! entered, so flag-off runs are byte-identical to before.

      use global_parameters
      use grav_func,     only: Dphi
      use species_table, only: n_bsp, bsp_fsp, bsp_nH, bsp_nHe,           &
                               bsp_charge, isp_HI, isp_HeI,               &
                               n_mion, mion_fsp, mion_stage,              &
                               n_melem, melem_i0, melem_top, melem_A
      use composition,   only: mass_per_H_nucleus_without_He

      implicit none
      private
      public :: element_diffusion_step, relax_element_composition
      public :: relative_settling_mass

      ! Species masses in the m_H units the code counts f_sp in (species_table
      ! bsp_mass literals), and the mass of that H = 1 unit in grams (the
      ! hydrogen ATOM as in parameters.f90 -- not the atomic mass unit u).
      real*8, parameter :: m_He_amu = 4.0d0
      real*8, parameter :: m_H_amu  = 1.0d0
      real*8, parameter :: m_amu_g  = 1.67353284d-24
      ! D_12 = D12_pref sqrt(T)/n,  D12_pref = 1.52e18 (1/m_H + 1/m_He)^(1/2)
      real*8, parameter :: D12_pref = 1.52d18*1.118033989d0       ! = 1.699e18
      ! Composition relaxation at a fixed wind (relax_element_composition):
      ! step budget, and the absolute convergence measure max|dX|/X_base.
      integer, parameter :: he_relax_maxstep = 400
      real*8,  parameter :: he_relax_tol     = 1.0d-12

      contains

      ! ------------------------------------------------------------------ !

      subroutine element_diffusion_step(rho, v, Tcode, f_sp, dt_code,     &
                                        closed_base, Jface_out, rhov_in)
      ! Advance the helium mass fraction X one relaxation step and project the
      ! new element totals back into f_sp.  rho, v, Tcode are the current
      ! adimensional primitives, dt_code the adimensional relaxation timestep;
      ! f_sp is modified in place.
      !
      ! closed_base (optional, default .false.) replaces the Dirichlet
      ! reservoir base by a zero-flux inner boundary, which is what the closed
      ! column of tests T1a/T4/T6 needs; production runs never set it.
      ! Jface_out (optional) returns the diffusive helium mass flux
      ! [g cm^-2 s^-1] at the faces r_edg(0:N) evaluated with the coefficients
      ! the step actually used and the NEW X, so that a discrete elemental
      ! budget closes exactly against it (test T1b).
      ! rhov_in (optional) is the advecting momentum density rho v
      ! [g cm^-2 s^-1] the composition is carried by.  Absent, it is the
      ! cell's own rho_j v_j, which is what the marching path wants; the
      ! relaxation at a fixed wind supplies the STEADY mass flux mdot/(4 pi
      ! r^2) instead -- see relax_element_composition for why.

      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho, v, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: dt_code
      logical, optional,                      intent(in)    :: closed_base
      real*8, dimension(0:N), optional,       intent(out)   :: Jface_out
      real*8, dimension(1-Ng:N+Ng), optional, intent(in)    :: rhov_in

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, Xhe, Xold
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, TK, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng) :: rp, rep, nH_phys, ntot_phys, dmeff
      real*8, dimension(1-Ng:N+Ng) :: nX, nXold, DcoX, GcoX, zbX, zb1, rhov
      real*8, dimension(0:N)       :: PdL, PdR, Jf
      ! NB: local scalars are checked against global_parameters case-
      ! insensitively.  In particular the time scale is tscale, NOT t0 (a local
      ! t0 would alias the global temperature normalization T0), and nothing
      ! here is named N, Ng, r, g, mu, info, count or du.
      real*8 :: tscale, X_base, m_1, mX, DprefX, fXbase, rXsc
      integer :: j, jlo, im, i0m, top, k
      logical :: shut_base

      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      shut_base = .false.
      if (present(closed_base)) shut_base = closed_base
      jlo = 2
      if (shut_base) jlo = 1

      m_1 = mass_per_H_nucleus_without_He()

      ! --- element nucleus counts per unit mass, and the helium mass fraction
      call element_nucleus_counts(f_sp, nucH, nucHe)
      msum = m_1*nucH + m_He_amu*nucHe        ! = 1 to round-off (see header)
      where (msum .lt. 1.0d-30) msum = 1.0d-30
      Xhe  = m_He_amu*nucHe/msum
      Xold = Xhe

      ! --- dimensional fields
      TK       = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      rp       = r*R0
      rep      = r_edg*R0
      rho_phys = rho*n0*m_amu_g*msum                       ! [g/cm^3]
      ntot_phys = (nucH + nucHe)*rho*n0                    ! [cm^-3] nuclei
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      Dco = D12_pref*sqrt(TK)/ntot_phys                    ! D_12 [cm^2/s]
      tscale  = R0/v0
      dt_phys = dt_code*tscale
      where (dt_phys .lt. 1.0d-30) dt_phys = 1.0d-30

      ! --- advecting momentum density rho v [g cm^-2 s^-1]
      if (present(rhov_in)) then
         rhov = rhov_in
      else
         rhov = rho_phys*v*v0
      endif

      ! --- settling coefficient G [1/cm] (gravity + ambipolar field + thermal)
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1)

      ! --- base reservoir composition (Dirichlet), X at the input He/H
      X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)

      ! --- implicit solve, with one Picard sweep on the settling prefactor
      call face_coefficients(Xhe, rho_phys, Dco, Gco, rp, PdL, PdR)
      call solve_mass_fraction(Xold, Xhe, rho_phys, dt_phys, rp, rep, rhov,&
                               PdL, PdR, X_base, jlo, shut_base)
      call face_coefficients(Xhe, rho_phys, Dco, Gco, rp, PdL, PdR)
      call solve_mass_fraction(Xold, Xhe, rho_phys, dt_phys, rp, rep, rhov,&
                               PdL, PdR, X_base, jlo, shut_base)

      ! Round-off assertion (NOT a limiter): the M-matrix argument in the
      ! header gives 0 <= X <= 1 per step; only round-off can leave the range.
      where (Xhe .lt. 0.0d0) Xhe = 0.0d0
      where (Xhe .gt. 1.0d0) Xhe = 1.0d0
      if (.not. shut_base) Xhe(1-Ng:1) = X_base    ! base + inner ghosts
      Xhe(N+1:N+Ng) = Xhe(N)                       ! zero-gradient outer ghost

      ! --- diffusive face flux actually carried by the step (diagnostic)
      if (present(Jface_out) .or. diffusion_check_on()) then
         do j = 0, N
            Jf(j) = PdL(j)*Xhe(j) + PdR(j)*Xhe(j+1)
         enddo
         if (present(Jface_out)) Jface_out = Jf
      endif
      if (diffusion_check_on()) then
         call write_element_flux_profile(rep, Jf, Xhe, rho_phys, v, dmeff)
      endif

      ! --- project the new element totals back into the species vector
      call project_elements(f_sp, Xhe, msum, nucH, nucHe, m_1)

      ! --- trace metals: each element diffuses against the (post-projection)
      ! hydrogen background with its own mass and mean charge.  Carried over
      ! unchanged from the validated trace kernel (it is a trace treatment by
      ! construction: a metal element cannot feed back on the background it
      ! diffuses through); only the background definition is the element count
      ! above rather than HI+HII.  Default OFF.
      if (he_metal_diffusion .and. thereis_metals) then
         call element_nucleus_counts(f_sp, nucH, nucHe)
         nH_phys = nucH*rho*n0
         where (nH_phys .lt. 1.0d-30) nH_phys = 1.0d-30
         do im = 1, n_melem
            i0m = melem_i0(im)
            top = melem_top(im)
            mX  = melem_A(im)
            nX  = 0.0d0
            zbX = 0.0d0
            do k = 0, top
               nX  = nX  + f_sp(:,mion_fsp(i0m+k))*rho*n0
               zbX = zbX + dble(k)*f_sp(:,mion_fsp(i0m+k))*rho*n0
            enddo
            where (nX .gt. 1.0d-30)
               zbX = zbX/nX
            elsewhere
               zbX = 0.0d0
            end where
            nXold = nX
            DprefX = 1.52d18*sqrt(1.0d0 + 1.0d0/mX)   ! Banks&Kockarts, X-in-H
            DcoX = DprefX*sqrt(TK)/ntot_phys
            if (he_ambipolar) then
               dmeff = (mX - m_H_amu) - 0.5d0*(zbX - zb1)
            else
               dmeff = mX - m_H_amu
            endif
            do j = 1-Ng, N+Ng
               GcoX(j) = dmeff(j)*m_amu_g*(Dphi(r(j))*v0*v0/R0)           &
                         /(kb_erg*TK(j))
            enddo
            if (he_alphaT .ne. 0.0d0) then
               do j = 2-Ng, N+Ng-1
                  GcoX(j) = GcoX(j) + he_alphaT*(log(TK(j+1))-log(TK(j-1)))&
                            / max((r(j+1)-r(j-1))*R0, 1.0d0)
               enddo
            endif
            fXbase = nXold(1)/nH_phys(1)                  ! reservoir metal/H
            call solve_trace_element_in_hydrogen(nX, nH_phys, DcoX, GcoX,  &
                                                 fXbase, dt_phys, rp, rep,  &
                                                 rhov/rho_phys)
            do j = 1-Ng, N+Ng
               ! Target density from the solved mixing ratio.  There is no
               ! cap at the reservoir ratio: settling piles an element up as
               ! readily as it depletes one, and clipping the pile-up is the
               ! same limiter the memo (section 1) rejects for helium.
               rXsc = nX(j)
               if (rXsc .lt. 0.0d0) rXsc = 0.0d0
               if (nXold(j) .gt. 1.0d-25*nH_phys(j)) then
                  rXsc = rXsc/nXold(j)                    ! scale factor
                  do k = 0, top
                     f_sp(j,mion_fsp(i0m+k)) = f_sp(j,mion_fsp(i0m+k))*rXsc
                  enddo
               else if (rXsc .gt. 1.0d-25*nH_phys(j)) then
                  ! element returned to an exhausted cell: re-seed (neutral)
                  f_sp(j,mion_fsp(i0m)) = rXsc/max(rho(j)*n0, 1.0d-30)
                  do k = 1, top
                     f_sp(j,mion_fsp(i0m+k)) = 0.0d0
                  enddo
               endif
            enddo
         enddo
      endif

      if (diffusion_check_on()) then
         call element_nucleus_counts(f_sp, nucH, nucHe)
         call report_step(Xhe, Jf, rho_phys, rp, rep,                     &
              maxval(abs(m_1*nucH(1:N) + m_He_amu*nucHe(1:N) - msum(1:N)) &
                     /msum(1:N)))
      endif

      end subroutine element_diffusion_step

      ! ------------------------------------------------------------------ !

      subroutine relax_element_composition(rho, v, Tcode, f_sp, omega,     &
                                           drift, nstep)
      ! Relax the element composition to its steady state in a FIXED wind.
      !
      ! The step size is the COMPOSITION time scale, not the hydro CFL step.
      ! The implicit operator of element_diffusion_step is unconditionally
      ! stable, so nothing ties its dt to a sound-crossing time, while the
      ! quantities it relaxes evolve on
      !
      !    dr^2/(D_12 + K_zz)   (diffusion)   and   dr/max(|v|, w_s)
      !
      ! with w_s = D_12 |G| the settling drift speed.  At the base of
      ! HD 209458 b the first of these is ~10^6 s against a ~1 s CFL step, so
      ! a relaxation driven at the CFL step moves the composition by ~10^-6 of
      ! the way per step and no practical number of steps reaches the steady
      ! state -- which is what left the outer loop drifting at ~0.75 per pass
      ! and hitting its pass limit.  Here dt starts at the smallest of those
      ! scales and is grown geometrically, so the late steps are direct steady
      ! solves and the coefficients are Picard-updated along the way.
      !
      ! The advecting flow is the STEADY mass flux, rho v = mdot/(4 pi r^2)
      ! with mdot the wind's own mass-loss rate (the median of 4 pi r^2 rho v
      ! over the escape window [j_min:N], which is the region the solver
      ! declares steady), NOT the cell's rho_j v_j.  The composition equation
      !
      !    rho (dX/dt + v dX/dr) = -(1/r^2) d(r^2 J)/dr
      !
      ! is the conservative equation only where rho and v satisfy continuity.
      ! Below ~1.02 R_p the converged states do not: measured on HD 209458 b
      ! and LHS 1140 b, the spread of r^2 rho v there is 10^2-10^4 times its
      ! own median, because the base carries a standing sound wave whose sign
      ! alternates from cell to cell.  Relaxed on that field, the advective
      ! form converges to the composition of a flow that neither conserves
      ! mass nor exists: on the HD 209458 b Kzz = 0 wind it put a 5-fold cliff
      ! at the fourth cell and a plateau at 0.10 of the reservoir ratio, where
      ! the barometric He-H separation scale is 24 cells.  With the steady
      ! flux the same wind relaxes to a base flat to 1% and a profile that
      ! leaves it on the barometric scale (0.97 at 1.05 R_p, 0.85 at 1.2 R_p).
      ! The marching path keeps rho_j v_j: there the wind is genuinely
      ! transient and the cell velocity is the consistent choice.
      !
      ! The stopping measure is ABSOLUTE: the largest change of the helium mass
      ! fraction over a step, divided by the reservoir value X_base.  The
      ! relative measure it replaces, max|dX|/X, is meaningless in a cell the
      ! transport has emptied.
      !
      ! UNDER-RELAXATION.  omega damps the composition update handed back to
      ! the outer Picard loop,
      !
      !    X <- X_old + omega (X_relaxed - X_old),      0 < omega <= 1
      !
      ! because that loop -- solve the wind at fixed composition, relax the
      ! composition at fixed wind -- has no damping of its own and the two can
      ! chase each other instead of converging (measured on the HD 209458 b
      ! Kzz = 0 wind: a limit cycle at 0.15-0.17 over twenty passes).  The
      ! caller owns the schedule.  The returned drift is the UNDAMPED distance
      ! max|X_relaxed - X_old|/X_base, not the damped step actually applied:
      ! it is the distance to the fixed point, so omega cannot buy a false
      ! convergence by making the applied step small.
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho, v, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,                                 intent(in)    :: omega
      real*8,                                 intent(out)   :: drift
      integer,                                intent(out)   :: nstep

      real*8, dimension(1-Ng:N+Ng) :: dt_code, nucH, nucHe, Xnow, Xpass
      real*8, dimension(1-Ng:N+Ng) :: msum, Xmix
      real*8, dimension(1-Ng:N+Ng) :: Xprev, TK, ntot, Dcl, Gcl, dmcl, zbcl
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, rhov, mflx
      real*8 :: m_1, X_base, tscale, drj, tdiff, tadv, wset, grow, dt_ref
      real*8 :: mflx_steady
      integer :: j, k

      drift = 0.0d0
      nstep = 0
      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      m_1    = mass_per_H_nucleus_without_He()
      X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)
      tscale = R0/v0

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      call element_nucleus_counts(f_sp, nucH, nucHe)
      ntot = (nucH + nucHe)*rho*n0
      where (ntot .lt. 1.0d0) ntot = 1.0d0
      Dcl = D12_pref*sqrt(TK)/ntot
      call settling_coefficient(rho, Tcode, f_sp, Gcl, dmcl, zbcl)

      ! Steady mass flux r^2 rho v [g cm^-1 s^-1], taken as the median over
      ! the escape window and imposed on the whole column as mflx_steady/r^2.
      rho_phys = rho*n0*m_amu_g*(m_1*nucH + m_He_amu*nucHe)
      mflx     = (r*R0)**2*rho_phys*v*v0
      mflx_steady = median_of(mflx(j_min:N))
      do j = 1-Ng, N+Ng
         rhov(j) = mflx_steady/((r(j)*R0)**2)
      enddo

      do j = 1, N
         drj   = max((r_edg(j) - r_edg(j-1))*R0, 1.0d0)
         tdiff = drj*drj/max(Dcl(j) + he_kzz, 1.0d-30)
         wset  = Dcl(j)*abs(Gcl(j))
         tadv  = drj/max(abs(rhov(j))/max(rho_phys(j), 1.0d-30),           &
                         wset, 1.0d-30)
         dt_code(j) = min(tdiff, tadv)/tscale
      enddo
      dt_code(1-Ng:0)   = dt_code(1)
      dt_code(N+1:N+Ng) = dt_code(N)
      dt_ref = minval(dt_code(1:N))

      call helium_mass_fraction(f_sp, m_1, Xpass)
      Xprev = Xpass
      grow  = 1.0d0
      do k = 1, he_relax_maxstep
         call element_diffusion_step(rho, v, Tcode, f_sp, dt_code*grow,    &
                                     rhov_in = rhov)
         call helium_mass_fraction(f_sp, m_1, Xnow)
         nstep = k
         if (maxval(abs(Xnow(1:N) - Xprev(1:N)))/X_base                   &
             .lt. he_relax_tol) exit
         Xprev = Xnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo
      drift = maxval(abs(Xnow(1:N) - Xpass(1:N)))/X_base

      ! Damped update: blend the relaxed composition with the one this pass
      ! started from and project the blend back into the species vector.  The
      ! projection is the same one element_diffusion_step ends on, applied to
      ! the CURRENT f_sp (which carries X_relaxed), so both element totals are
      ! met exactly and no species can go negative.
      if (omega .lt. 1.0d0) then
         Xmix = Xpass + omega*(Xnow - Xpass)
         where (Xmix .lt. 0.0d0) Xmix = 0.0d0
         where (Xmix .gt. 1.0d0) Xmix = 1.0d0
         call element_nucleus_counts(f_sp, nucH, nucHe)
         msum = m_1*nucH + m_He_amu*nucHe
         where (msum .lt. 1.0d-30) msum = 1.0d-30
         call project_elements(f_sp, Xmix, msum, nucH, nucHe, m_1)
      endif

      ! dt_ref is the shortest composition time scale of the grid; report it
      ! so the log shows what the relaxation was measured against.
      if (diffusion_check_on())                                           &
         write(0,'(A,ES10.3,A,ES12.5,A,I0)')                              &
              ' (diffusion) relaxation dt_0 = ', dt_ref*tscale,           &
              ' s, steady 4pi r^2 rho v = ', 4.0d0*pi*mflx_steady,        &
              ' g/s, steps = ', nstep

      end subroutine relax_element_composition

      ! ------------------------------------------------------------------ !

      real*8 function median_of(a)
      ! Median of a, by insertion sort of a local copy (the arrays here are
      ! one grid long and this runs once per relaxation).
      real*8, dimension(:), intent(in) :: a
      real*8, dimension(size(a)) :: b
      real*8  :: t
      integer :: i, k, m
      b = a
      do i = 2, size(b)
         t = b(i)
         k = i - 1
         do while (k .ge. 1)
            if (b(k) .le. t) exit
            b(k+1) = b(k)
            k = k - 1
         enddo
         b(k+1) = t
      enddo
      m = size(b)
      if (mod(m,2) .eq. 0) then
         median_of = 0.5d0*(b(m/2) + b(m/2+1))
      else
         median_of = b((m+1)/2)
      endif
      end function median_of

      ! ------------------------------------------------------------------ !

      subroutine helium_mass_fraction(f_sp, m_1, Xhe)
      ! Helium mass fraction X = m_He n_He/(m_1 n_H + m_He n_He) of the species
      ! vector -- the same definition element_diffusion_step transports.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8,                                 intent(in)  :: m_1
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Xhe
      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum

      call element_nucleus_counts(f_sp, nucH, nucHe)
      msum = m_1*nucH + m_He_amu*nucHe
      where (msum .lt. 1.0d-30) msum = 1.0d-30
      Xhe = m_He_amu*nucHe/msum

      end subroutine helium_mass_fraction

      ! ------------------------------------------------------------------ !

      subroutine element_nucleus_counts(f_sp, nucH, nucHe)
      ! Element (nucleus) counts per unit mass, over every species that
      ! carries them: n_H/rho = f_HI + f_HII + 2 f_H2 + 2 f_H2+ + 3 f_H3+
      ! + f_HeH+, n_He/rho = f_HeI + f_HeII + f_HeIII + f_HeTR + f_HeH+.
      ! Columns a run does not carry are zero, so no flag has to be tested.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucH, nucHe
      integer :: ib

      nucH  = 0.0d0
      nucHe = 0.0d0
      do ib = 1, n_bsp
         if (bsp_nH(ib)  .gt. 0)                                          &
            nucH  = nucH  + dble(bsp_nH(ib)) *f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0)                                          &
            nucHe = nucHe + dble(bsp_nHe(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      where (nucH  .lt. 0.0d0) nucH  = 0.0d0
      where (nucHe .lt. 0.0d0) nucHe = 0.0d0

      end subroutine element_nucleus_counts

      ! ------------------------------------------------------------------ !

      subroutine mean_charges_and_electrons(f_sp, nucH, nucHe, zb1, zbHe,  &
                                            ne_rel)
      ! Mean charge of each component per nucleus and the free-electron
      ! density per unit mass.  The component charges use the species that
      ! carry ONE element only (the atomic region: HI/HII against
      ! HeI/HeII/HeIII/HeTR); HeH+, which carries both, is left out of both
      ! mean charges as a trace approximation and the metal charge is left out
      ! of component 1 (metal/H ~ 1e-4).  The electron density, by contrast,
      ! counts every charge under the same policy as utils/calc_ne, because it
      ! is the physical n_e the ambipolar field is built from.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: nucH, nucHe
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: zb1, zbHe
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: ne_rel
      integer :: ib, im

      zb1    = 0.0d0
      zbHe   = 0.0d0
      ne_rel = 0.0d0
      do ib = 1, n_bsp
         if (bsp_charge(ib) .eq. 0) cycle
         if (bsp_nH(ib) .gt. 0 .and. bsp_nHe(ib) .eq. 0)                  &
            zb1  = zb1  + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0 .and. bsp_nH(ib) .eq. 0)                  &
            zbHe = zbHe + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
         ne_rel = ne_rel + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      if (eos_include_metals .and. thereis_metals) then
         do im = 1, n_mion
            if (mion_stage(im) .gt. 0)                                    &
               ne_rel = ne_rel + dble(mion_stage(im))*f_sp(:,mion_fsp(im))
         enddo
      endif
      zb1  = zb1 /max(nucH,  1.0d-30)
      zbHe = zbHe/max(nucHe, 1.0d-30)

      end subroutine mean_charges_and_electrons

      ! ------------------------------------------------------------------ !

      subroutine settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1)
      ! The settling coefficient G of the header, and the gravity-normalized
      ! relative settling mass it comes from,
      !
      !   dmeff = (m_He - m_1) - (Zbar_He - Zbar_1) eE/(m_H g)   [m_H units]
      !   G     = dmeff m_H g/(k T) + alpha_T dlnT/dr            [1/cm]
      !
      ! with eE = -k T dln(n_e T)/dr (central difference, one-sided at the
      ! ends).  dmeff is returned so the ambipolar limits can be tested
      ! against 3 (neutral), 2.5 (H+ plasma) and 5/3 (He++ plasma) without a
      ! second definition of the field anywhere.  Where the local gravity
      ! vanishes dmeff is reported as zero (it is a ratio to g); G itself is
      ! always built from the forces, never from dmeff.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Gco, dmeff, zb1

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, zbHe, ne_rel, ne_phys
      real*8, dimension(1-Ng:N+Ng) :: TK, gphys, eEf, lnpe
      real*8 :: m_1, dr2
      integer :: j

      m_1 = mass_per_H_nucleus_without_He()
      call element_nucleus_counts(f_sp, nucH, nucHe)
      call mean_charges_and_electrons(f_sp, nucH, nucHe, zb1, zbHe, ne_rel)

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      ne_phys = ne_rel*rho*n0
      do j = 1-Ng, N+Ng
         gphys(j) = Dphi(r(j))*v0*v0/R0                    ! [cm/s^2], inward
      enddo

      ! Ambipolar field eE = -k T dln(n_e T)/dr [erg/cm].
      eEf = 0.0d0
      if (he_ambipolar) then
         lnpe = log(max(ne_phys*TK, 1.0d-300))
         do j = 2-Ng, N+Ng-1
            dr2 = (r(j+1) - r(j-1))*R0
            eEf(j) = -kb_erg*TK(j)*(lnpe(j+1) - lnpe(j-1))/dr2
         enddo
         eEf(1-Ng) = eEf(2-Ng)
         eEf(N+Ng) = eEf(N+Ng-1)
      endif

      do j = 1-Ng, N+Ng
         Gco(j) = ((m_He_amu - m_1)*m_amu_g*gphys(j)                      &
                   - (zbHe(j) - zb1(j))*eEf(j))/(kb_erg*TK(j))
         if (abs(gphys(j)) .gt. 0.0d0) then
            dmeff(j) = (m_He_amu - m_1)                                   &
                       - (zbHe(j) - zb1(j))*eEf(j)/(m_amu_g*gphys(j))
         else
            dmeff(j) = 0.0d0
         endif
      enddo

      ! Thermal diffusion: alpha_T dlnT/dr, default he_alphaT = 0 (no-op).
      if (he_alphaT .ne. 0.0d0) then
         do j = 2-Ng, N+Ng-1
            Gco(j) = Gco(j) + he_alphaT*(log(TK(j+1)) - log(TK(j-1)))     &
                     / max((r(j+1)-r(j-1))*R0, 1.0d0)
         enddo
      endif

      end subroutine settling_coefficient

      ! ------------------------------------------------------------------ !

      function relative_settling_mass(rho, Tcode, f_sp) result(dmeff)
      ! Diagnostic: the gravity-normalized relative settling mass of helium
      ! against component 1, in m_H units (3 neutral, 2.5 in an H+ plasma,
      ! 5/3 in a fully ionized He++ plasma).  Same definition the operator
      ! uses -- there is no second copy of the ambipolar field.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng) :: dmeff
      real*8, dimension(1-Ng:N+Ng) :: Gco, zb1

      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1)

      end function relative_settling_mass

      ! ------------------------------------------------------------------ !

      subroutine face_coefficients(Xhe, rho_phys, Dco, Gco, rp, PdL, PdR)
      ! Face coefficients of the diffusive flux at r_edg(0:N).  The diffusive
      ! flux at face j is
      !     J(j) = PdL(j) X(j) + PdR(j) X(j+1)
      ! Faces 0 and N are the two boundaries and carry zero diffusive flux.
      ! The advection is not a face quantity here -- it is the cell-velocity
      ! upwind difference assembled in solve_mass_fraction (see the header) --
      ! and the outer ghost carries the top cell's own composition, so an
      ! inflowing top boundary brings in gas of the same composition instead of
      ! a silent zero flux.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Xhe, rho_phys, Dco, Gco
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rp
      real*8, dimension(0:N),       intent(out) :: PdL, PdR

      real*8 :: dr_f, rhof, Df, Gf, Xf, Agrad, Wset
      integer :: j

      PdL = 0.0d0 ; PdR = 0.0d0

      do j = 1, N-1
         dr_f = max(rp(j+1)-rp(j), 1.0d0)
         rhof = 0.5d0*(rho_phys(j)+rho_phys(j+1))
         Df   = 0.5d0*(Dco(j)+Dco(j+1))
         Gf   = 0.5d0*(Gco(j)+Gco(j+1))
         Xf   = 0.5d0*(Xhe(j)+Xhe(j+1))
         if (Xf .lt. 0.0d0) Xf = 0.0d0
         if (Xf .gt. 1.0d0) Xf = 1.0d0
         Agrad = rhof*(Df + he_kzz)              ! gradient + eddy, A + E >= 0
         Wset  = rhof*Df*Gf*(1.0d0 - Xf)         ! settling, flux = -Wset*X_f
         ! Peclet hybrid: central where the settling drift is resolved,
         ! upwind toward the settling direction otherwise.  Both branches keep
         ! PdL >= 0 and PdR <= 0, which is the M-matrix condition.
         if (abs(Wset)*dr_f .lt. 2.0d0*Agrad) then
            PdL(j) =  Agrad/dr_f - 0.5d0*Wset
            PdR(j) = -Agrad/dr_f - 0.5d0*Wset
         else if (Wset .ge. 0.0d0) then         ! settles inward: donor is j+1
            PdL(j) =  Agrad/dr_f
            PdR(j) = -Agrad/dr_f - Wset
         else                                   ! rises: donor is j
            PdL(j) =  Agrad/dr_f - Wset
            PdR(j) = -Agrad/dr_f
         endif
      enddo

      end subroutine face_coefficients

      ! ------------------------------------------------------------------ !

      subroutine solve_mass_fraction(Xold, Xhe, rho_phys, dt_phys, rp, rep, &
                                     rhov, PdL, PdR, X_base, jlo, shut_base)
      ! One backward-Euler tridiagonal solve of
      !   rho (X^new - X^old)/dt + div(r^2 J)/r^2 + rho v dX/dr = 0
      ! for rows jlo..N, the advective term entering as the cell-velocity
      ! upwind difference of the header.  Xold is the state at the start of the
      ! step (so the routine may be called twice with the same Xold for the
      ! Picard sweep) and Xhe carries the iterate in and the solution out.
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: Xold, rho_phys, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: rp, rep, rhov
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: Xhe
      real*8, dimension(0:N),       intent(in)    :: PdL, PdR
      real*8,                       intent(in)    :: X_base
      integer,                      intent(in)    :: jlo
      logical,                      intent(in)    :: shut_base

      real*8, dimension(1-Ng:N+Ng) :: aa, bb, cc, dd, cpv, dpv
      real*8 :: Kj, mden, sL, sR, cadv
      integer :: j

      do j = jlo, N
         Kj    = 1.0d0/(rp(j)**2*max(rep(j)-rep(j-1), 1.0d0))
         sL    = rep(j-1)**2
         sR    = rep(j)**2
         aa(j) = -Kj*sL*PdL(j-1)
         bb(j) =  rho_phys(j)/dt_phys(j)                                   &
                + Kj*(sR*PdL(j) - sL*PdR(j-1))
         cc(j) =  Kj*sR*PdR(j)
         dd(j) =  rho_phys(j)*Xold(j)/dt_phys(j)

         ! Advective term rho v dX/dr, one-sided upwind on the CELL velocity:
         ! +cadv on the diagonal and -cadv on the donor neighbour, so the
         ! advective coefficients of the row sum to zero (uniform X exact) and
         ! the only off-diagonal it creates is <= 0 (M-matrix preserved).
         if (rhov(j) .ge. 0.0d0) then
            ! Dropped at a closed inner boundary, which is zero-gradient.
            if (j .gt. jlo .or. .not. shut_base) then
               cadv  = rhov(j)/max(rp(j) - rp(j-1), 1.0d0)
               bb(j) = bb(j) + cadv
               aa(j) = aa(j) - cadv
            endif
         else
            ! Dropped at the top, whose zero-gradient ghost makes it vanish.
            if (j .lt. N) then
               cadv  = rhov(j)/max(rp(j+1) - rp(j), 1.0d0)
               bb(j) = bb(j) - cadv
               cc(j) = cc(j) + cadv
            endif
         endif
      enddo
      cc(N) = 0.0d0

      if (jlo .eq. 2) then
         dd(2) = dd(2) - aa(2)*X_base
         aa(2) = 0.0d0
      endif

      cpv(jlo) = cc(jlo)/bb(jlo)
      dpv(jlo) = dd(jlo)/bb(jlo)
      do j = jlo+1, N
         mden   = bb(j) - aa(j)*cpv(j-1)
         cpv(j) = cc(j)/mden
         dpv(j) = (dd(j) - aa(j)*dpv(j-1))/mden
      enddo
      Xhe(N) = dpv(N)
      do j = N-1, jlo, -1
         Xhe(j) = dpv(j) - cpv(j)*Xhe(j+1)
      enddo
      if (jlo .eq. 2) Xhe(1-Ng:1) = X_base
      Xhe(N+1:N+Ng) = Xhe(N)

      end subroutine solve_mass_fraction

      ! ------------------------------------------------------------------ !

      subroutine project_elements(f_sp, Xhe, msum, nucH_old, nucHe_old, m_1)
      ! Nonnegative projection of the species vector onto the new element
      ! totals (memo section 3).  Helium-only species are scaled by
      ! r_He = n_He^new/n_He^old, hydrogen-only species (H2, H2+, H3+ included)
      ! and the slaved metals by r_H, and a species carrying both elements
      ! (HeH+) by min(r_H, r_He); the nuclei this under-counts for the element
      ! with the larger factor are deposited into that element's neutral
      ! ground species.  Every operation is a multiplication by a nonnegative
      ! factor or an addition, so no negative intermediate can arise, and both
      ! element totals are met exactly.  The result is the initial guess for
      ! ioniz_eq, which owns the split WITHIN an element; this step owns the
      ! element totals.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: Xhe, msum
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: nucH_old
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: nucHe_old
      real*8,                                 intent(in)    :: m_1

      real*8  :: nucH_new, nucHe_new, rH, rHe, rBoth, gotH, gotHe
      integer :: j, ib, im

      do j = 1-Ng, N+Ng
         nucHe_new = Xhe(j)*msum(j)/m_He_amu
         nucH_new  = (1.0d0 - Xhe(j))*msum(j)/m_1
         if (nucH_old(j) .gt. 1.0d-30) then
            rH = nucH_new/nucH_old(j)
         else
            rH = 0.0d0
         endif
         if (nucHe_old(j) .gt. 1.0d-30) then
            rHe = nucHe_new/nucHe_old(j)
         else
            rHe = 0.0d0
         endif
         if (rH .lt. 0.0d0)  rH  = 0.0d0
         if (rHe .lt. 0.0d0) rHe = 0.0d0
         rBoth = min(rH, rHe)

         do ib = 1, n_bsp
            if (bsp_nH(ib) .gt. 0 .and. bsp_nHe(ib) .gt. 0) then
               f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rBoth
            else if (bsp_nH(ib) .gt. 0) then
               f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rH
            else if (bsp_nHe(ib) .gt. 0) then
               f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rHe
            endif
         enddo
         do im = 1, n_mion                       ! metals slaved to hydrogen
            f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
         enddo

         ! deposit the shortfall of each element into its neutral ground stage
         gotH  = 0.0d0
         gotHe = 0.0d0
         do ib = 1, n_bsp
            if (bsp_nH(ib)  .gt. 0)                                       &
               gotH  = gotH  + dble(bsp_nH(ib)) *f_sp(j,bsp_fsp(ib))
            if (bsp_nHe(ib) .gt. 0)                                       &
               gotHe = gotHe + dble(bsp_nHe(ib))*f_sp(j,bsp_fsp(ib))
         enddo
         if (nucH_new  .gt. gotH)                                         &
            f_sp(j,isp_HI)  = f_sp(j,isp_HI)  + (nucH_new  - gotH)
         if (nucHe_new .gt. gotHe)                                        &
            f_sp(j,isp_HeI) = f_sp(j,isp_HeI) + (nucHe_new - gotHe)
      enddo

      end subroutine project_elements

      ! ------------------------------------------------------------------ !

      subroutine solve_trace_element_in_hydrogen(nX, nHl, Dco, Gco, fbase,  &
                                                 dt_phys, rp, rep, vadv)
      ! One backward-Euler tridiagonal solve for a TRACE element diffusing
      ! relative to a hydrogen background nHl, with element diffusion
      ! coefficient Dco, settling coefficient Gco and fixed reservoir base
      ! ratio fbase = nX/nHl, advected by vadv [cm/s] (the same flow as the
      ! helium equation, divided by rho).  Valid only for an element whose own
      ! mass does
      ! not shape the background it moves through (metal/H ~ 1e-4); helium is
      ! NOT solved this way -- that is what the binary operator above exists
      ! for.
      !
      ! The transported variable is the MIXING RATIO fX = nX/n_H, in the same
      ! advective form as the helium equation,
      !
      !    dfX/dt + v dfX/dr = (1/(n_H r^2)) d/dr [ r^2 n_H
      !                        ( (D_1X + K_zz) dfX/dr + D_1X G fX ) ]
      !
      ! and the advection is the cell-velocity upwind difference of the header.
      ! Solving the DENSITY conservatively with face-averaged velocities, which
      ! is what this routine did before, evacuates any cell whose two face
      ! velocities point outward: such a cell has an outflow term at both faces
      ! and no inflow term at either, and the drain has no counterpart in the
      ! hydro's own rho, whose face fluxes come from the Riemann solver and
      ! carry no such divergence.  At the breathing base of HD 209458 b that
      ! emptied the first free cell of every metal by ~10^3 (measured
      ! 2026-08-25) and with it the metal-line cooling of that cell.  In the
      ! mixing-ratio advective form a uniform fX is preserved for any velocity
      ! field, so no spurious divergence can create or destroy the element.
      ! nX is intent(inout): supply the current density, receive the solved one.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: nX
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: nHl, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: rp, rep, vadv
      real*8,                       intent(in)    :: fbase

      real*8, dimension(0:N)       :: PL, PR
      real*8, dimension(1-Ng:N+Ng) :: aa, bb, cc, dd, cpv, dpv, fX
      real*8 :: dr_f, nHf, Df, Gf, DK, Kj, mden, sL, sR, cadv
      integer :: j

      fX = nX/max(nHl, 1.0d-30)

      PL = 0.0d0
      PR = 0.0d0
      do j = 1, N-1
         dr_f = max(rp(j+1)-rp(j), 1.0d0)
         nHf  = 0.5d0*(nHl(j)+nHl(j+1))
         Df   = 0.5d0*(Dco(j)+Dco(j+1))
         Gf   = 0.5d0*(Gco(j)+Gco(j+1))
         DK   = Df + he_kzz
         ! Peclet hybrid, as in face_coefficients: central where the settling
         ! drift is resolved, upwind toward the settling direction otherwise.
         if (abs(Df*Gf)*dr_f .lt. 2.0d0*DK) then
            PL(j) =  nHf*(DK/dr_f - 0.5d0*Df*Gf)
            PR(j) = -nHf*(DK/dr_f + 0.5d0*Df*Gf)
         else if (Gf .ge. 0.0d0) then          ! settles inward: donor is j+1
            PL(j) =  nHf*(DK/dr_f)
            PR(j) = -nHf*(DK/dr_f + Df*Gf)
         else                                  ! rises: donor is j
            PL(j) =  nHf*(DK/dr_f - Df*Gf)
            PR(j) = -nHf*(DK/dr_f)
         endif
      enddo

      do j = 2, N
         Kj    = 1.0d0/(max(nHl(j), 1.0d-30)*rp(j)**2                      &
                        *max(rep(j)-rep(j-1), 1.0d0))
         sL    = rep(j-1)**2
         sR    = rep(j)**2
         aa(j) = -Kj*sL*PL(j-1)
         bb(j) =  1.0d0/dt_phys(j) + Kj*(sR*PL(j) - sL*PR(j-1))
         cc(j) =  Kj*sR*PR(j)
         dd(j) =  fX(j)/dt_phys(j)
         if (vadv(j) .ge. 0.0d0) then
            cadv  = vadv(j)/max(rp(j) - rp(j-1), 1.0d0)
            bb(j) = bb(j) + cadv
            aa(j) = aa(j) - cadv
         else if (j .lt. N) then
            cadv  = vadv(j)/max(rp(j+1) - rp(j), 1.0d0)
            bb(j) = bb(j) - cadv
            cc(j) = cc(j) + cadv
         endif
      enddo
      cc(N) = 0.0d0

      dd(2) = dd(2) - aa(2)*fbase                         ! base Dirichlet
      aa(2) = 0.0d0

      cpv(2) = cc(2)/bb(2)
      dpv(2) = dd(2)/bb(2)
      do j = 3, N
         mden   = bb(j) - aa(j)*cpv(j-1)
         cpv(j) = cc(j)/mden
         dpv(j) = (dd(j) - aa(j)*dpv(j-1))/mden
      enddo
      fX(N) = dpv(N)
      do j = N-1, 2, -1
         fX(j) = dpv(j) - cpv(j)*fX(j+1)
      enddo
      fX(1-Ng:1) = fbase
      fX(N+1:N+Ng) = fX(N)
      ! Round-off assertion, not a limiter: the M-matrix keeps fX >= 0.
      where (fX .lt. 0.0d0) fX = 0.0d0
      nX = fX*nHl

      end subroutine solve_trace_element_in_hydrogen

      ! ------------------------------------------------------------------ !

      logical function diffusion_check_on()
      ! EXHALE_DIFFUSION_CHECK=1 turns on the stderr diagnostic written each step.
      character(len=8) :: envv
      call get_environment_variable('EXHALE_DIFFUSION_CHECK', envv)
      diffusion_check_on = (trim(envv) .eq. '1')
      end function diffusion_check_on

      ! ------------------------------------------------------------------ !

      subroutine write_element_flux_profile(rep, Jf, Xhe, rho_phys, v,     &
                                            dmeff)
      ! Radial profile of the elemental helium face flux carried by the step
      ! (EXHALE_DIFFUSION_CHECK=1), written to ./diffusion_faceflux.txt and
      ! replaced at every call, so that after a run the file holds the state
      ! the run ended on.  This is the T8 diagnostic of
      ! docs/binary_diffusion_design.md section 6:
      !
      !   F_He(r_f) = 4 pi r_f^2 ( rho X v + J )   [g/s]
      !
      ! with the advective part rho_f v_f X_upwind (the physical face flux --
      ! the operator carries the advection as a cell-velocity upwind difference,
      ! which has no face representation) and J from the diffusive face
      ! coefficients with the new X, so the diffusive column is the flux the
      ! step actually used and not a re-derived one.  The total mass flux
      ! 4 pi r_f^2 rho_f v_f is written beside it: at a steady state both are
      ! constant with radius, and the comparison of their radial spreads is
      ! the pass criterion.  dmeff is the face-averaged relative settling
      ! mass of the ambipolar diagnostic (3 neutral, 2.5 H+ plasma, 5/3 He++
      ! plasma).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rep, Xhe, rho_phys, v
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dmeff
      real*8, dimension(0:N),       intent(in) :: Jf
      real*8  :: area, adv, rhof, vf, dmf
      integer :: j, uu

      uu = 771
      open(unit=uu, file='diffusion_faceflux.txt', status='replace')
      write(uu,'(A)') '# face flux of the helium element carried by '//    &
         'binary_element_diffusion (EXHALE_DIFFUSION_CHECK=1)'
      write(uu,'(A)') '# columns: j  r_face[R_p]  X_face  J_diff'//        &
         '[g/cm2/s]  F_adv[g/cm2/s]  F_He=4pi r^2 (F_adv+J)[g/s]'//        &
         '  Mdot_face=4pi r^2 rho v[g/s]  dmeff_face'
      do j = 1, N-1
         area = 4.0d0*pi*rep(j)**2
         rhof = 0.5d0*(rho_phys(j) + rho_phys(j+1))
         vf   = 0.5d0*(v(j) + v(j+1))*v0
         if (vf .ge. 0.0d0) then
            adv = rhof*vf*Xhe(j)
         else
            adv = rhof*vf*Xhe(j+1)
         endif
         dmf  = 0.5d0*(dmeff(j) + dmeff(j+1))
         write(uu,'(I6,7ES16.7)') j, rep(j)/R0, 0.5d0*(Xhe(j)+Xhe(j+1)),  &
              Jf(j), adv, area*(adv + Jf(j)), area*rhof*vf, dmf
      enddo
      close(uu)

      end subroutine write_element_flux_profile

      ! ------------------------------------------------------------------ !

      subroutine report_step(Xhe, Jf, rho_phys, rp, rep, mass_resid)
      ! Diagnostic written each step (EXHALE_DIFFUSION_CHECK=1): the range of X, the
      ! largest |J_He + J_1| over the faces, the largest relative error of the
      ! two-component mass closure m_1 n_H + m_He n_He = rho after the
      ! write-back, and the total helium mass in the domain.
      !
      ! A binary mixture has ONE independent diffusive flux: the operator
      ! computes J_He and the hydrogen component carries -J_He, so the printed
      ! |J_He + J_1| is zero by construction and is reported for completeness.
      ! The quantity that actually guards the discretization is the mass
      ! closure: it is what keeps rho -- owned by the hydro -- consistent with
      ! the composition the operator hands back, and it is nonzero the moment
      ! the write-back stops meeting both element totals exactly.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: Xhe, rho_phys, rp, rep
      real*8, dimension(0:N),       intent(in) :: Jf
      real*8,                       intent(in) :: mass_resid
      real*8  :: mHe_tot, sumJ
      integer :: j

      mHe_tot = 0.0d0
      do j = 1, N
         mHe_tot = mHe_tot                                                &
                 + rho_phys(j)*Xhe(j)*rp(j)**2*(rep(j)-rep(j-1))
      enddo
      sumJ = 0.0d0
      do j = 0, N
         sumJ = max(sumJ, abs(Jf(j) - Jf(j)))
      enddo
      write(0,'(A,ES13.6,A,ES13.6,A,ES9.2,A,ES9.2,A,ES16.9)')             &
         ' (diffusion) X min ', minval(Xhe(1:N)), ' max ',                &
         maxval(Xhe(1:N)), ' max|J_He+J_1| ', sumJ,                       &
         ' mass closure ', mass_resid, ' He mass ', 4.0d0*pi*mHe_tot

      end subroutine report_step

      ! End of module
      end module binary_element_diffusion

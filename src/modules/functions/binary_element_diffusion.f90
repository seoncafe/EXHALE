      module binary_element_diffusion
      ! Binary (two-component) H/He element transport for the single-fluid
      ! EXHALE wind.  Formulation: docs/binary_diffusion_design.md (sections
      ! 2-5), including the molecular-region closure of its section 5
      ! (milestone M4).
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
      !    G = [ (m_He - m_c1) m_H g - (Zbar_He - Zbar_c1) eE ] / (k T)
      !        + alpha_T dlnT/dr  -  dln(psi)/dr                      [1/cm]
      !
      ! and the hydrogen-component flux is its negative.  The masses and
      ! charges in G are CARRIER quantities: m_c1 is the mass of component 1
      ! per collision partner and Zbar_c1 its mean charge per collision
      ! partner (see CARRIERS below), which in the atomic region are the mass
      ! and charge per nucleus, m_1 and Zbar_1.
      !
      ! The gradient coefficient of the memo,
      ! A = (m_c1 m_He n_car/mbar_car) D_12 (dx/dX), collapses to rho D_12,
      ! and it does so with the carriers as well as with the nuclei.  Write
      ! n_car = n_1c + n_He for the collision partners, psi = n_1c/n_H for the
      ! carriers each hydrogen nucleus is spread over (1 atomic, 1/2 fully
      ! H2), m_c1 = m_1/psi and mbar_car = rho/n_car.  Then
      ! x = n_He/n_car obeys x = (X/m_He)/(X/m_He + psi (1-X)/m_1), so at
      ! frozen psi -- the chemistry is lagged over the step, as every other
      ! coefficient is -- dx/dX = psi mbar_car^2/(m_1 m_He), while
      ! rho_1 rho_He/(rho x (1-x)) = m_c1 m_He n_car/mbar_car.  The product is
      ! rho D_12 (m_c1 psi/m_1) = rho D_12 exactly, for any psi.  So the
      ! gradient term of the molecular region is the same Fickian rho D dX/dr
      ! as the atomic one.
      !
      ! What the frozen-psi derivative leaves out is the piece of dx/dr the
      ! CHEMISTRY carries, and that piece is a genuine driver: the
      ! Chapman-Cowling force is the gradient of the MOLE fraction x, and
      !
      !    logit(x) = logit(X) + ln(m_1/(m_He psi)),
      !
      ! so d logit(x)/dr = d logit(X)/dr - dln(psi)/dr.  Where the hydrogen
      ! turns molecular going down, each nucleus is spread over fewer
      ! collision partners, helium's mole fraction rises, and helium diffuses
      ! down that gradient -- outward across the molecular front -- even at a
      ! uniform mass fraction.  The term therefore rides in G as
      ! -dln(psi)/dr, which leaves the discretization, the Peclet hybrid and
      ! the M-matrix argument untouched and vanishes identically wherever psi
      ! is constant, i.e. over the whole atomic region.  The eddy term does
      ! not carry it: eddy mixing transports the mixture as a whole and has no
      ! preferred species, so it acts on dX/dr alone.
      !
      ! Nothing divides by a vanishing element density anywhere -- that
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
      ! DIFFUSION COEFFICIENT -- STAGE-RESOLVED (memo 2.6, decision D3).  The
      ! friction between the two elements is built pair by pair over their
      ! IONIZATION STAGES, because a helium ion moving through a proton gas is
      ! held by Coulomb collisions whose momentum-transfer cross section at
      ! T ~ 1e4 K is ~1e-12 cm^2 against the ~1e-15 cm^2 of a hard sphere.  A
      ! single neutral coefficient lets helium and the metals settle in the
      ! ionized wind by two to three orders of magnitude too fast.  Each pair
      ! coefficient is the Chapman-Enskog first approximation
      !
      !    D_st = 3 k T / (16 n mu_st Omega_st^(1,1))
      !
      ! with n the total CARRIER density (the Chapman-Enskog n is a density of
      ! colliding particles; in the atomic region it is the nucleus density,
      ! and where the hydrogen is molecular it is smaller by the nuclei bound
      ! into each molecule) and mu_st the reduced mass, evaluated
      ! in the limit the pair belongs to (routines below):
      !   neutral-neutral   hard_sphere_pair_diffusion    ~ T^(1/2)/n
      !   ion-neutral       ion_neutral_pair_diffusion    (polarization and
      !                                                    rigid core added)
      !   ion-ion           coulomb_pair_diffusion        ~ T^(5/2)/n
      ! and the element-element coefficient is the stage-fraction weighted
      ! harmonic sum (frictions add; Blanc's law across stages)
      !
      !    1/D_eff = sum_{s in He} sum_{t in comp.1} y_s y_t / D_st .
      !
      ! Limits: an all-neutral gas leaves one term, the Banks & Kockarts
      ! hard-sphere coefficient the earlier milestones used everywhere; a
      ! fully ionized gas leaves the He++/H+ Coulomb term alone (test T12).
      ! Resonant charge exchange H+ + H is NOT a term here: both partners
      ! carry the same element, and a binary diffusion coefficient is driven
      ! by the friction BETWEEN the elements, so every ion-neutral pair that
      ! appears is non-resonant.  HeH+ carries both elements and is left out
      ! of both carrier lists, and the metal nuclei are left out of n and of
      ! the reduced mass, both as trace approximations.
      !
      ! CARRIERS, AND THE MOLECULAR REGION (memo section 5, milestone M4).
      ! Equation (1) transports ELEMENTS, and the chemistry moves hydrogen
      ! between H, H+, H2, H2+, H3+ and helium between its stages without
      ! changing either element's mass, so the transport equation is the same
      ! above and below the molecular front.  What the species set changes is
      ! the CLOSURE, because the collision partners of helium are then a
      ! mixture of carriers:
      !
      !   friction   the stage sum above already runs over the carrier lists
      !              hcar_isp = {HI, HII, H2, H2+, H3+} and hecar_isp, each
      !              carrier entering with its own mass (an H2 is one partner
      !              of mass 2 m_H), charge and polarizability.  Blanc's law
      !              across carriers is the same harmonic sum as across
      !              stages, so it reduces to D(He,H) in the atomic region and
      !              to D(He,H2) in a fully molecular one, and interpolates
      !              between (test T7);
      !   forces     gravity and the ambipolar field act on the particles that
      !              actually collide, so G uses the mean carrier mass m_c1
      !              and mean carrier charge Zbar_c1 of the hydrogen component
      !              (per nucleus the weight of H2 equals that of H; per
      !              collision partner it is twice), and the trace-metal loop
      !              settles a metal against the same mean carrier mass;
      !   density    n in every pair coefficient is the carrier density.
      !
      ! The molecular columns of f_sp are zero unless the molecular chemistry
      ! is on, so no flag is tested anywhere: with it off every carrier list,
      ! mean carrier mass and carrier density collapses onto the nuclei and
      ! the atomic-region behavior is recovered identically.
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
                               bsp_is_excited_level,                      &
                               bsp_charge, isp_HI, isp_HeI,               &
                               isp_HII, isp_HeII, isp_HeIII, isp_HeTR,    &
                               isp_H2, isp_H2p, isp_H3p,                  &
                               n_mion, mion_fsp, mion_stage,              &
                               n_melem, melem_i0, melem_top, melem_A
      use composition,   only: mass_per_H_nucleus_without_He
      use lower_atmosphere_profile, only: lap_in_use, lap_r_top_RJ,       &
                          lap_flux_measured, lap_flux_window_empty,       &
                          lap_flux_nface, lap_FH_median, lap_FH_spread,   &
                          lap_FHe_median, lap_FHe_spread,                 &
                          lap_Mdot_median, lap_Mdot_spread,               &
                          lap_flux_r_lo_Rp, lap_flux_r_hi_Rp,             &
                          lap_flux_r_lo_measured,                         &
                          lap_steady_r_lo_Rp, lap_steady_nface,           &
                          lap_steady_FH_median,  lap_steady_FH_spread,    &
                          lap_steady_FHe_median, lap_steady_FHe_spread,   &
                          lap_steady_Mdot_median, lap_steady_Mdot_spread

      implicit none
      private
      public :: element_diffusion_step, relax_element_composition
      public :: relative_settling_mass
      ! Exposed so the acceptance tests read the same coefficients the
      ! operator uses -- there is no second copy of the friction anywhere.
      public :: helium_hydrogen_diffusion
      public :: hard_sphere_pair_diffusion, polarization_pair_diffusion
      public :: ion_neutral_pair_diffusion
      public :: coulomb_pair_diffusion, coulomb_logarithm
      public :: alpha_HI, alpha_HeI

      ! Species masses in the m_H units the code counts f_sp in (species_table
      ! bsp_mass literals), and the mass of that H = 1 unit in grams (the
      ! hydrogen ATOM as in parameters.f90 -- not the atomic mass unit u).
      real*8, parameter :: m_He_amu = 4.0d0
      real*8, parameter :: m_H_amu  = 1.0d0
      real*8, parameter :: m_amu_g  = 1.67353284d-24
      ! Elementary charge in electrostatic units, from the CODATA 2018 exact
      ! coulomb value: e = 1.602176634e-19 C / (10/c) = 4.803204713e-10 esu.
      real*8, parameter :: e_esu    = 4.803204713d-10
      ! Hard-sphere prefactor of Banks & Kockarts (1973): the neutral binary
      ! coefficient is D = 1.52e18 (1/A_s + 1/A_t)^(1/2) sqrt(T)/n.
      real*8, parameter :: bk_hs_pref = 1.52d18

      ! Static dipole polarizabilities of the NEUTRAL collision partners, the
      ! only atomic datum the ion-neutral (polarization) coefficient needs.
      ! Quoted in atomic units a_0^3 and converted with
      ! a_0^3 = 1.481847e-25 cm^3.  alpha(H) = 4.5 a_0^3 is the exact
      ! nonrelativistic ground-state value; the rest are the recommended
      ! static dipole polarizabilities of Schwerdtfeger & Nagle (2019,
      ! Mol. Phys. 117, 1200).  As a cross-check H and He come out at 0.667
      ! and 0.205 Angstrom^3, the values tabulated for aeronomy by
      ! Schunk & Nagy (Ionospheres, Table 4.1).
      real*8, parameter :: a0cub_cm3 = 1.481847d-25
      real*8, parameter :: alpha_HI  =   4.500d0*a0cub_cm3
      real*8, parameter :: alpha_HeI =   1.383d0*a0cub_cm3
      real*8, parameter :: alpha_H2  =   5.315d0*a0cub_cm3
      ! neutral metals, in the melem order C O N Mg Si Ca Na K S Fe
      real*8, parameter :: alpha_melem(n_melem) = a0cub_cm3*                &
           [  11.3d0,   5.3d0,   7.4d0,  71.2d0,  37.3d0,                   &
             160.8d0, 162.7d0, 289.7d0,  19.4d0,  62.0d0 ]

      ! Carriers of component 1 (hydrogen).  H2/H2+/H3+ are zero unless the
      ! molecular chemistry is on, so no flag has to be tested; with it off
      ! the list collapses to HI/HII.  Masses are per CARRIER (an H2 is one
      ! collision partner of mass 2), which is what the friction counts.
      integer, parameter :: n_hcar = 5
      integer, parameter :: hcar_isp(n_hcar) =                              &
           [ isp_HI, isp_HII, isp_H2, isp_H2p, isp_H3p ]
      real*8,  parameter :: hcar_Z(n_hcar) =                                &
           [ 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0 ]
      real*8,  parameter :: hcar_m(n_hcar) =                                &
           [ 1.0d0, 1.0d0, 2.0d0, 2.0d0, 3.0d0 ]
      real*8,  parameter :: hcar_alpha(n_hcar) =                            &
           [ alpha_HI, 0.0d0, alpha_H2, 0.0d0, 0.0d0 ]
      character(len=5), parameter :: hcar_name(n_hcar) =                    &
           [ 'HI   ', 'HII  ', 'H2   ', 'H2+  ', 'H3+  ' ]

      ! Carriers of the helium component.  The 2^3S triplet is a neutral
      ! helium atom for the friction -- same mass, same charge, same
      ! polarizability as a ground-state He I atom -- and it is an EXCITED
      ! LEVEL of He I, so the HeI column already counts it
      ! (bsp_is_excited_level).  It is therefore not a list entry of its own:
      ! a separate HeTR row made every triplet atom two collision partners.
      integer, parameter :: n_hecar = 3
      integer, parameter :: hecar_isp(n_hecar) =                            &
           [ isp_HeI, isp_HeII, isp_HeIII ]
      real*8,  parameter :: hecar_Z(n_hecar) =                              &
           [ 0.0d0, 1.0d0, 2.0d0 ]
      real*8,  parameter :: hecar_m(n_hecar) =                              &
           [ 4.0d0, 4.0d0, 4.0d0 ]
      real*8,  parameter :: hecar_alpha(n_hecar) =                          &
           [ alpha_HeI, 0.0d0, 0.0d0 ]
      character(len=5), parameter :: hecar_name(n_hecar) =                  &
           [ 'HeI  ', 'HeII ', 'HeIII' ]

      ! Largest number of ionization stages carried for one metal element.
      integer, parameter :: n_mstage = 3
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
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, dlnpsi
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, TK, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng) :: rp, rep, nH_phys, ntot_phys, dmeff
      real*8, dimension(1-Ng:N+Ng) :: nX, nXold, DcoX, GcoX, zbX, zb1, rhov
      real*8, dimension(1-Ng:N+Ng) :: Dneut, eEf, ne_phys, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng,n_hcar)    :: yH
      real*8, dimension(1-Ng:N+Ng,n_mstage)  :: yX
      real*8, dimension(n_mstage)  :: ZXs, mXs, alXs
      integer, dimension(1-Ng:N+Ng):: idom
      integer, dimension(n_mstage) :: ispX
      real*8, dimension(0:N)       :: PdL, PdR, Jf
      ! NB: local scalars are checked against global_parameters case-
      ! insensitively.  In particular the time scale is tscale, NOT t0 (a local
      ! t0 would alias the global temperature normalization T0), and nothing
      ! here is named N, Ng, r, g, mu, info, count or du.
      real*8 :: tscale, X_base, m_1, mX, fXbase, rXsc
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
      ! Carrier (collision-partner) density, which is what the Chapman-Enskog
      ! coefficients divide by.  Equal to the nucleus density in the atomic
      ! region; smaller where the hydrogen is bound into molecules.
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ntot_phys = (carH + carHe)*rho*n0                    ! [cm^-3] carriers
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      ! D_12 [cm^2/s], stage-resolved (module header / memo 2.6)
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dco, Dneut, idom)
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
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)

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
      if (present(Jface_out) .or. diffusion_check_on() .or. lap_in_use) then
         do j = 0, N
            Jf(j) = PdL(j)*Xhe(j) + PdR(j)*Xhe(j+1)
         enddo
         if (present(Jface_out)) Jface_out = Jf
      endif
      ! Written on demand (EXHALE_DIFFUSION_CHECK=1) and unconditionally
      ! whenever a lower-atmosphere profile is in use: there the elemental
      ! fluxes are not a diagnostic but the quantity the two models have to
      ! agree on, so the run must always leave them behind.
      if (diffusion_check_on() .or. lap_in_use) then
         call write_element_flux_profile(rep, Jf, Xhe, rho_phys, v, dmeff, &
                                         Dco, Dneut, idom)
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
         call carrier_counts(f_sp, carH, carHe, mcarH)
         ntot_phys = (carH + carHe)*rho*n0
         where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
         call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
         ne_phys = max(ne_rel*rho*n0, 1.0d0)
         call carrier_fractions(f_sp, n_hcar, hcar_isp, yH)
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
            ! Stage-resolved friction against the hydrogen carriers, exactly
            ! as for helium: a metal ion in the ionized wind is held to the
            ! protons by the Coulomb coefficient, which is what keeps it from
            ! settling out (Koskinen et al. 2013, section 3.2.2).  Stages
            ! beyond this element's top are given zero weight.
            ispX = mion_fsp(i0m)
            ZXs  = 0.0d0
            mXs  = mX
            alXs = 0.0d0
            alXs(1) = alpha_melem(im)
            do k = 0, top
               ispX(k+1) = mion_fsp(i0m+k)
               ZXs(k+1)  = dble(k)
            enddo
            call carrier_fractions(f_sp, top+1, ispX(1:top+1),             &
                                   yX(:,1:top+1))
            if (top+1 .lt. n_mstage) yX(:,top+2:n_mstage) = 0.0d0
            call stage_mixture_diffusion(TK, ntot_phys, ne_phys,           &
                                         top+1, yX(:,1:top+1),             &
                                         ZXs(1:top+1), mXs(1:top+1),       &
                                         alXs(1:top+1),                    &
                                         n_hcar, yH, hcar_Z, hcar_m,       &
                                         hcar_alpha, DcoX)
            ! Same ambipolar field as helium, computed from the electron
            ! pressure gradient (memo 3a) -- NOT the hydrogen-plasma constant
            ! eE = m_H g/2 this loop used to assume, which is wrong wherever
            ! helium or the metals carry a significant share of the electrons.
            ! The partner mass is the mean HYDROGEN mass per carrier: a metal
            ! atom settles against the particles it collides with, which are H
            ! atoms and ions above the molecular front and H2/H3+ below it.
            ! mcarH is 1 in the atomic region, so nothing changes there.
            ! -dln(psi)/dr enters here for the same reason it enters the
            ! helium equation: the driver is the metal's MOLE fraction among
            ! the hydrogen carriers, n_X/n_1c, while the transported variable
            ! is the mixing ratio n_X/n_H (module header).
            do j = 1-Ng, N+Ng
               GcoX(j) = ((mX - mcarH(j))*m_amu_g*(Dphi(r(j))*v0*v0/R0)   &
                          - (zbX(j) - zb1(j))*eEf(j))/(kb_erg*TK(j))      &
                         - dlnpsi(j)
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
      real*8, dimension(1-Ng:N+Ng) :: Xprev, TK, Dcl, Gcl, dmcl, zbcl
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
      ! Only the composition TIME SCALES below are built from this; the
      ! operator recomputes its own D_12 each step from the current state.
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dcl)
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
         tdiff = drj*drj/max(Dcl(j) + kzz_cell(j), 1.0d-30)
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
      ! + f_HeH+, n_He/rho = f_HeI + f_HeII + f_HeIII + f_HeH+.  Columns a
      ! run does not carry are zero, so no flag has to be tested.  He 2^3S
      ! is skipped: it is an excited level of He I and its nucleus is
      ! already inside the HeI column (bsp_is_excited_level).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucH, nucHe
      integer :: ib

      nucH  = 0.0d0
      nucHe = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib)  .gt. 0)                                          &
            nucH  = nucH  + dble(bsp_nH(ib)) *f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0)                                          &
            nucHe = nucHe + dble(bsp_nHe(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      where (nucH  .lt. 0.0d0) nucH  = 0.0d0
      where (nucHe .lt. 0.0d0) nucHe = 0.0d0

      end subroutine element_nucleus_counts

      ! ------------------------------------------------------------------ !

      subroutine carrier_counts(f_sp, carH, carHe, mcarH)
      ! Collision-partner (CARRIER) counts per unit mass of the two
      ! components, and the mean hydrogen mass per carrier of component 1:
      !
      !   carH  = sum_c f_c            over hcar_isp  = HI HII H2 H2+ H3+
      !   carHe = sum_s f_s            over hecar_isp = HeI HeII HeIII
      !   mcarH = sum_c m_c f_c / carH   [m_H per carrier]
      !
      ! The friction and the force terms of the diffusion equation count
      ! PARTICLES, not nuclei: an H2 molecule is one collision partner and
      ! carries twice the weight of an H atom (memo section 5).  With the
      ! molecular chemistry off the molecular columns are zero, so carH and
      ! carHe are the nucleus counts and mcarH = 1 exactly, which is why no
      ! flag is tested here and why the atomic region is untouched.
      !
      ! HeH+ carries a nucleus of each element and is in neither list, the
      ! same trace approximation that leaves it out of the mean charges.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: carH, carHe
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: mcarH

      real*8, dimension(1-Ng:N+Ng) :: mHsum
      integer :: ic

      carH  = 0.0d0
      carHe = 0.0d0
      mHsum = 0.0d0
      do ic = 1, n_hcar
         carH  = carH  + max(f_sp(:,hcar_isp(ic)), 0.0d0)
         mHsum = mHsum + hcar_m(ic)*max(f_sp(:,hcar_isp(ic)), 0.0d0)
      enddo
      do ic = 1, n_hecar
         carHe = carHe + max(f_sp(:,hecar_isp(ic)), 0.0d0)
      enddo
      if (present(mcarH)) then
         where (carH .gt. 0.0d0)
            mcarH = mHsum/carH
         elsewhere
            mcarH = m_H_amu
         end where
      endif

      end subroutine carrier_counts

      ! ------------------------------------------------------------------ !
      !  Pair diffusion coefficients (memo 2.6, eqs. 8-10)
      !
      !  All three are the Chapman-Enskog first approximation
      !      D_st = 3 k T / (16 n mu_st Omega_st^(1,1)),
      !      Omega^(1,1) = sqrt(kT/2 pi mu) int_0^inf e^-gam^2 gam^5 Q^(1) dgam,
      !      gam^2 = mu g^2/(2 k T),
      !  evaluated with the momentum-transfer cross section Q^(1) of the
      !  interaction the pair actually has.  n is the total nucleus density in
      !  every case, so the three are directly comparable and the mixture rule
      !  below can add their inverses.
      ! ------------------------------------------------------------------ !

      real*8 function hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
      ! NEUTRAL-NEUTRAL.  Banks & Kockarts (1973), Aeronomy Part B: the
      ! binary coefficient of two neutrals of atomic weights A_s, A_t is
      !
      !    D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2) / n     [cm^2/s]
      !
      ! symmetric in the two species, so one number serves He-in-H and
      ! H-in-He.  Consistency check of the framework: the rigid-sphere
      ! collision integral Omega^(1,1) = pi d^2 (kT/2 pi mu)^(1/2) put into
      ! the Chapman-Enskog D above returns exactly this expression with the
      ! collision diameter d = 2.7 Angstrom, so (8), (9) and (10) are the same
      ! approximation and not three unrelated fits.
      real*8, intent(in) :: TK, ntot, A_s, A_t
      hard_sphere_pair_diffusion = bk_hs_pref*sqrt(1.0d0/A_s + 1.0d0/A_t)  &
                                   *sqrt(TK)/ntot
      end function hard_sphere_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
      ! ION-NEUTRAL, NON-RESONANT.  The ion polarizes the neutral, giving the
      ! Langevin interaction V = -alpha_n e^2/(2 r^4).  For that potential the
      ! product g Q^(1)(g) is independent of the relative speed,
      !
      !    g Q^(1) = 2.21 pi e (alpha_n/mu)^(1/2),
      !
      ! the constant behind the non-resonant ion-neutral momentum-transfer
      ! collision frequency of Schunk & Nagy (Ionospheres, eq. 4.88),
      ! nu_in = 2.21 pi (n_n m_n/(m_i+m_n)) (gamma_n e^2/mu_in)^(1/2).  The
      ! textbook is not in references/, so the constant below was re-derived
      ! from the Chapman-Enskog integral rather than copied: with
      ! Q^(1) = C/g the integral gives Omega^(1,1) = 3 C/16, and the 3/16 of
      ! the Chapman-Enskog D cancels it exactly, leaving
      !
      !    D = k T / (2.21 pi e n (alpha_n mu)^(1/2))        [cm^2/s]
      !
      ! alpha_n [cm^3] is the neutral partner's static dipole polarizability,
      ! e is in esu and mu_g the reduced mass in grams.  Note D ~ T, not
      ! T^(1/2): the polarization coupling weakens as the pair gets faster.
      ! Valid for NON-RESONANT pairs only; a resonant pair (H+ + H, He+ + He)
      ! is charge-exchange dominated and never appears here, because both of
      ! its partners carry the same element (see the module header).
      !
      ! TEMPERATURE RANGE.  This is the LOW-energy limit of the ion-neutral
      ! interaction: it keeps only the induced-dipole attraction and no
      ! repulsive core, so its friction weakens as T^-1 where a rigid core
      ! would hold it at T^-1/2.  It is therefore NOT called on its own --
      ! ion_neutral_pair_diffusion below adds it to the rigid-core channel,
      ! which is what keeps the friction from vanishing at high temperature.
      real*8, intent(in) :: TK, ntot, alpha_n, mu_g
      polarization_pair_diffusion = kb_erg*TK                              &
           /(2.21d0*pi*e_esu*ntot*sqrt(max(alpha_n*mu_g, 1.0d-60)))
      end function polarization_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function coulomb_pair_diffusion(TK, ntot, ne, Z_s, Z_t, mu_g)
      ! ION-ION.  Screened Coulomb, momentum-transfer cross section
      ! Q^(1) = 4 pi b_90^2 ln(Lambda) with b_90 = Z_s Z_t e^2/(mu g^2) the
      ! impact parameter of a 90 degree deflection.  Carrying it through the
      ! Chapman-Enskog Omega^(1,1) integral (int gam e^-gam^2 dgam = 1/2)
      ! gives
      !
      !    D = 3 (k T)^(5/2)
      !        / [ 4 (2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln(Lambda) ]
      !
      ! the standard T^(5/2)/(n Z^2 Z^2 ln Lambda mu^(1/2)) Chapman-Enskog
      ! Coulomb result (Paquette et al. 1986, ApJS 61, 177, their section II;
      ! Schunk & Nagy eq. 4.142 give the equivalent collision frequency).
      ! Neither is in references/, so the numerical constant 3/(4 sqrt(2 pi))
      ! was re-derived here from the collision integral, not copied.
      real*8, intent(in) :: TK, ntot, ne, Z_s, Z_t, mu_g
      real*8 :: zz
      zz = abs(Z_s*Z_t)
      coulomb_pair_diffusion = 3.0d0*(kb_erg*TK)**2.5d0                    &
           /(4.0d0*sqrt(2.0d0*pi*mu_g)*ntot*(zz*e_esu*e_esu)**2            &
             *coulomb_logarithm(TK, ne, zz))
      end function coulomb_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function coulomb_logarithm(TK, ne, zz)
      ! ln(Lambda) = ln(lambda_D/b_90) with the Debye length screened by the
      ! electrons, lambda_D = (k T/(4 pi n_e e^2))^(1/2), and the 90 degree
      ! impact parameter evaluated at the thermal relative energy
      ! mu g^2 = 3 k T:
      !
      !    ln(Lambda) = ln( 3 k T lambda_D / (Z_s Z_t e^2) )
      !
      ! (Spitzer 1962, Physics of Fully Ionized Gases, section 5.2).  Floored
      ! at 1: the expansion the Coulomb cross section comes from needs
      ! Lambda >> 1, and a cell where it is not is one where the ion-ion pair
      ! carries no weight anyway.
      real*8, intent(in) :: TK, ne, zz
      real*8 :: lam_D
      lam_D = sqrt(kb_erg*TK/(4.0d0*pi*max(ne, 1.0d0)*e_esu*e_esu))
      coulomb_logarithm = max(log(3.0d0*kb_erg*TK*lam_D                    &
                                  /(max(zz, 1.0d0)*e_esu*e_esu)), 1.0d0)
      end function coulomb_logarithm

      ! ------------------------------------------------------------------ !

      real*8 function ion_neutral_pair_diffusion(TK, ntot, alpha_n,        &
                                                 A_s, A_t, mu_g)
      ! ION-NEUTRAL, NON-RESONANT: the induced-dipole (polarization) channel
      ! and the rigid-core channel taken TOGETHER,
      !
      !    1/D_in = 1/D_polarization + 1/D_hard sphere .
      !
      ! Reason: the two are momentum-transfer cross sections of the same
      ! encounter -- the long-range attraction and the short-range repulsion
      ! of one interaction potential -- and to first order their Q^(1) add,
      ! so their collision integrals add, so their FRICTIONS add.  Adding
      ! frictions is adding inverse diffusion coefficients (D = 3kT/(16 n mu
      ! Omega^(1,1)) is linear in 1/Omega), which is the same rule
      ! stage_mixture_diffusion uses across stages.  The physical content is
      ! that opening a second channel cannot make a pair MORE mobile: D_in is
      ! bounded above by the weaker-friction channel and can never exceed
      ! either limit.
      !
      ! Both limits are reproduced exactly.  D_pol ~ T and D_hs ~ T^(1/2), so
      ! their ratio crosses unity near 1.5e3 K for an H/He pair: below that
      ! the polarization term dominates the friction and D_in -> D_pol; above
      ! it the rigid core does and D_in -> D_hs.  At the 1e4 K of an ionized
      ! wind D_pol is 2.6 to 4.7 times D_hs, so the combined value sits within
      ! 20-30% of the hard sphere -- taking the polarization value alone there
      ! would understate the friction by a factor of a few and let an element
      ! settle correspondingly faster through the partially ionized layer,
      ! which is where ion-neutral pairs carry the friction at all.
      !
      ! The rigid-core channel uses the same Banks & Kockarts collision
      ! diameter as the neutral-neutral pair: no separate ion-neutral core
      ! radius is available, and the ionic radius differs from the atomic one
      ! by much less than the factor the combination is settling.
      real*8, intent(in) :: TK, ntot, alpha_n, A_s, A_t, mu_g
      real*8 :: Dpol, Dhs

      Dpol = polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
      Dhs  = hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
      ion_neutral_pair_diffusion = 1.0d0                                   &
           /(1.0d0/max(Dpol, 1.0d-99) + 1.0d0/max(Dhs, 1.0d-99))

      end function ion_neutral_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function stage_pair_diffusion(TK, ntot, ne, Z_s, m_s, al_s,   &
                                                       Z_t, m_t, al_t)
      ! The pair coefficient of two ionization stages, in whichever of the
      ! three limits above the pair belongs to.  Masses in m_H units,
      ! polarizabilities in cm^3 (used only when that partner is the neutral
      ! one).
      real*8, intent(in) :: TK, ntot, ne, Z_s, m_s, al_s, Z_t, m_t, al_t
      real*8 :: mu_g

      mu_g = m_s*m_t/(m_s + m_t)*m_amu_g
      if (Z_s .eq. 0.0d0 .and. Z_t .eq. 0.0d0) then
         stage_pair_diffusion = hard_sphere_pair_diffusion(TK, ntot,       &
                                                           m_s, m_t)
      else if (Z_s .ne. 0.0d0 .and. Z_t .ne. 0.0d0) then
         stage_pair_diffusion = coulomb_pair_diffusion(TK, ntot, ne,       &
                                                       Z_s, Z_t, mu_g)
      else if (Z_s .eq. 0.0d0) then
         stage_pair_diffusion = ion_neutral_pair_diffusion(TK, ntot,       &
                                                 al_s, m_s, m_t, mu_g)
      else
         stage_pair_diffusion = ion_neutral_pair_diffusion(TK, ntot,       &
                                                 al_t, m_s, m_t, mu_g)
      endif

      end function stage_pair_diffusion

      ! ------------------------------------------------------------------ !

      subroutine stage_mixture_diffusion(TK, ntot, ne, nsA, yA, ZA, mA,    &
                                         alA, nsB, yB, ZB, mB, alB,        &
                                         Deff, idom)
      ! Effective binary coefficient between two ELEMENTS whose ionization
      ! stages are resolved (memo eq. 11):
      !
      !    1/D_eff = sum_s sum_t y_s y_t / D_st
      !
      ! The element moves as one body (one velocity per element, not per
      ! stage), so its friction against the other element is the sum of the
      ! frictions its stages feel.  Frictions add, hence the INVERSE
      ! coefficients are averaged, weighted by the carrier fractions y.  This
      ! is Blanc's law, exact for a trace species in a mixture and the
      ! standard approximation otherwise; it is the same rule section 5 of
      ! the memo adopts across molecular carriers, applied across stages too.
      ! With one stage populated in each element it returns that pair exactly.
      !
      ! idom (optional) returns the packed index (s-1)*nsB + t of the pair
      ! carrying the largest share of the friction, for the diagnostic
      ! profile.
      integer,                         intent(in)  :: nsA, nsB
      real*8, dimension(1-Ng:N+Ng),    intent(in)  :: TK, ntot, ne
      real*8, dimension(1-Ng:N+Ng,nsA),intent(in)  :: yA
      real*8, dimension(1-Ng:N+Ng,nsB),intent(in)  :: yB
      real*8, dimension(nsA),          intent(in)  :: ZA, mA, alA
      real*8, dimension(nsB),          intent(in)  :: ZB, mB, alB
      real*8, dimension(1-Ng:N+Ng),    intent(out) :: Deff
      integer, dimension(1-Ng:N+Ng), optional, intent(out) :: idom

      real*8  :: fric, fpair, wgt, Dst, fbest
      integer :: j, is, it, kbest

      do j = 1-Ng, N+Ng
         fric  = 0.0d0
         fbest = 0.0d0
         kbest = 1
         do is = 1, nsA
            if (yA(j,is) .le. 0.0d0) cycle
            do it = 1, nsB
               if (yB(j,it) .le. 0.0d0) cycle
               wgt = yA(j,is)*yB(j,it)
               Dst = stage_pair_diffusion(TK(j), ntot(j), ne(j),           &
                                          ZA(is), mA(is), alA(is),         &
                                          ZB(it), mB(it), alB(it))
               fpair = wgt/max(Dst, 1.0d-99)
               fric  = fric + fpair
               if (fpair .gt. fbest) then
                  fbest = fpair
                  kbest = (is-1)*nsB + it
               endif
            enddo
         enddo
         if (fric .gt. 0.0d0) then
            Deff(j) = 1.0d0/fric
         else
            ! Neither element has a populated stage in this cell: the
            ! settling prefactor vanishes there anyway, and the gradient term
            ! only needs a finite coefficient.  Take the neutral pair.
            Deff(j) = hard_sphere_pair_diffusion(TK(j), ntot(j),           &
                                                 mA(1), mB(1))
         endif
         if (present(idom)) idom(j) = kbest
      enddo

      end subroutine stage_mixture_diffusion

      ! ------------------------------------------------------------------ !

      subroutine carrier_fractions(f_sp, ns, isp_list, yc)
      ! Fraction of an element's collision partners in each of its stages,
      ! sum_s y_s = 1.  f_sp holds particle counts per unit mass, so a ratio
      ! of f_sp entries is a ratio of particles -- which is what the friction
      ! weights, not a ratio of nuclei (an H2 is one collision partner).
      ! A cell with no carrier at all is given to the first (neutral) stage,
      ! so y is always normalized and D_eff always finite.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      integer,                                intent(in)  :: ns
      integer, dimension(ns),                 intent(in)  :: isp_list
      real*8, dimension(1-Ng:N+Ng,ns),        intent(out) :: yc

      real*8, dimension(1-Ng:N+Ng) :: tot
      integer :: is

      tot = 0.0d0
      do is = 1, ns
         yc(:,is) = max(f_sp(:,isp_list(is)), 0.0d0)
         tot      = tot + yc(:,is)
      enddo
      do is = 1, ns
         where (tot .gt. 0.0d0)
            yc(:,is) = yc(:,is)/tot
         elsewhere
            yc(:,is) = 0.0d0
         end where
      enddo
      where (tot .le. 0.0d0) yc(:,1) = 1.0d0

      end subroutine carrier_fractions

      ! ------------------------------------------------------------------ !

      subroutine helium_hydrogen_diffusion(rho, Tcode, f_sp, Deff, Dneut,  &
                                           idom)
      ! The stage-resolved binary coefficient between helium and component 1,
      ! which is what the transport operator carries as D_12.  Dneut is the
      ! all-neutral coefficient of the same cell, returned so the diagnostic
      ! and the reports can quote the suppression D_eff/D_neutral directly.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Deff
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: Dneut
      integer, dimension(1-Ng:N+Ng), optional, intent(out) :: idom

      real*8, dimension(1-Ng:N+Ng)        :: TK, ntot, carH, carHe
      real*8, dimension(1-Ng:N+Ng)        :: zb1, zbHe, ne_rel, ne_phys
      real*8, dimension(1-Ng:N+Ng,n_hecar):: yHe
      real*8, dimension(1-Ng:N+Ng,n_hcar) :: yH
      integer :: j

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      call carrier_counts(f_sp, carH, carHe)
      ntot = (carH + carHe)*rho*n0
      where (ntot .lt. 1.0d0) ntot = 1.0d0
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ne_phys = max(ne_rel*rho*n0, 1.0d0)

      call carrier_fractions(f_sp, n_hecar, hecar_isp, yHe)
      call carrier_fractions(f_sp, n_hcar,  hcar_isp,  yH)
      call stage_mixture_diffusion(TK, ntot, ne_phys,                      &
                                   n_hecar, yHe, hecar_Z, hecar_m,         &
                                   hecar_alpha,                            &
                                   n_hcar,  yH,  hcar_Z,  hcar_m,          &
                                   hcar_alpha, Deff, idom)

      if (present(Dneut)) then
         do j = 1-Ng, N+Ng
            Dneut(j) = hard_sphere_pair_diffusion(TK(j), ntot(j),          &
                                                  m_He_amu, m_H_amu)
         enddo
      endif

      end subroutine helium_hydrogen_diffusion

      ! ------------------------------------------------------------------ !

      subroutine mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ! Mean charge of each component per CARRIER and the free-electron
      ! density per unit mass.  The charge that enters the ambipolar force of
      ! the diffusion equation is the charge of a colliding particle, so the
      ! normalization is the carrier count, not the nucleus count: an H3+
      ! carries one charge for three hydrogen nuclei.  In the atomic region
      ! carriers and nuclei are the same particles and this is the charge per
      ! nucleus, as before.
      !
      ! The component charges use the species that carry ONE element only
      ! (HI/HII/H2/H2+/H3+ against HeI/HeII/HeIII/HeTR); HeH+, which carries
      ! both, is left out of both mean charges as a trace approximation and
      ! the metal charge is left out of component 1 (metal/H ~ 1e-4).  The
      ! electron density, by contrast, counts every charge under the same
      ! policy as utils/calc_ne, because it is the physical n_e the ambipolar
      ! field is built from.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: zb1, zbHe
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: ne_rel
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe
      integer :: ib, im

      zb1    = 0.0d0
      zbHe   = 0.0d0
      ne_rel = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
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
      call carrier_counts(f_sp, carH, carHe)
      zb1  = zb1 /max(carH,  1.0d-30)
      zbHe = zbHe/max(carHe, 1.0d-30)

      end subroutine mean_charges_and_electrons

      ! ------------------------------------------------------------------ !

      subroutine settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1,   &
                                      eEout, dlnpsi_out)
      ! The settling coefficient G of the header, and the gravity-normalized
      ! relative settling mass it comes from,
      !
      !   dmeff = (m_He - m_c1) - (Zbar_He - Zbar_c1) eE/(m_H g) [m_H units]
      !   G     = dmeff m_H g/(k T) + alpha_T dlnT/dr            [1/cm]
      !
      ! with CARRIER quantities: m_c1 = m_1 mcarH is the mass of component 1
      ! per collision partner (m_1 is its mass per hydrogen nucleus and mcarH
      ! the hydrogen nuclei bound into one carrier), and Zbar_c1, Zbar_He are
      ! the mean charges per carrier.  Gravity and the ambipolar field act on
      ! the particles that collide, so an H2 molecule weighs twice an H atom
      ! here even though it carries the same mass per nucleus (memo section
      ! 5).  In the atomic region mcarH = 1 and this is the nucleus form.
      !
      ! with eE = -k T dln(n_e T)/dr (central difference, one-sided at the
      ! ends).  dmeff is returned so the ambipolar limits can be tested
      ! against 3 (neutral), 2.5 (H+ plasma) and 5/3 (He++ plasma) without a
      ! second definition of the field anywhere.  Where the local gravity
      ! vanishes dmeff is reported as zero (it is a ratio to g); G itself is
      ! always built from the forces, never from dmeff.  dmeff is the FORCE
      ! ratio only: the chemistry term -dln(psi)/dr that G also carries in the
      ! molecular region is not a mass and does not belong in it.
      !
      ! psi = n_1c/n_H is the number of collision partners component 1 spreads
      ! one hydrogen nucleus over: 1 in the atomic region, 1/2 where the
      ! hydrogen is H2.  Its gradient is the mole-fraction driver of the
      ! module header, and it is returned (dlnpsi_out) so that the trace-metal
      ! loop uses the one definition of it.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Gco, dmeff, zb1
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: eEout
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: dlnpsi_out

      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng) :: ne_phys, mc1, lnmcar, dlnpsi
      real*8, dimension(1-Ng:N+Ng) :: TK, gphys, eEf, lnpe
      real*8 :: m_1, dr2
      integer :: j

      m_1 = mass_per_H_nucleus_without_He()
      call carrier_counts(f_sp, carH, carHe, mcarH)
      mc1 = m_1*mcarH
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)

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

      ! Mole-fraction driver of the molecular region, -dln(psi)/dr with
      ! psi = n_1c/n_H the collision partners per hydrogen nucleus (module
      ! header).  psi = 1 identically wherever the hydrogen is atomic, so the
      ! central difference below is exactly zero there.
      dlnpsi = 0.0d0
      lnmcar = log(max(mcarH, 1.0d-30))      ! ln psi = -ln mcarH
      do j = 2-Ng, N+Ng-1
         dlnpsi(j) = -(lnmcar(j+1) - lnmcar(j-1))                         &
                     /max((r(j+1) - r(j-1))*R0, 1.0d0)
      enddo
      dlnpsi(1-Ng) = dlnpsi(2-Ng)
      dlnpsi(N+Ng) = dlnpsi(N+Ng-1)

      do j = 1-Ng, N+Ng
         Gco(j) = ((m_He_amu - mc1(j))*m_amu_g*gphys(j)                   &
                   - (zbHe(j) - zb1(j))*eEf(j))/(kb_erg*TK(j))            &
                  - dlnpsi(j)
         if (abs(gphys(j)) .gt. 0.0d0) then
            dmeff(j) = (m_He_amu - mc1(j))                                &
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

      ! The field itself, so that every element settling against component 1
      ! -- helium here, the trace metals below -- uses ONE ambipolar field
      ! and there is no second, cruder copy of it anywhere in the module.
      if (present(eEout)) eEout = eEf
      if (present(dlnpsi_out)) dlnpsi_out = dlnpsi

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

      real*8 :: dr_f, rhof, Df, Gf, Kf, Xf, Agrad, Wset
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
         Kf   = 0.5d0*(kzz_cell(j)+kzz_cell(j+1))
         Agrad = rhof*(Df + Kf)                  ! gradient + eddy, A + E >= 0
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
            if (bsp_is_excited_level(ib)) cycle
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
         DK   = Df + 0.5d0*(kzz_cell(j)+kzz_cell(j+1))
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
                                            dmeff, Dco, Dneut, idom)
      ! Radial profile of the elemental face fluxes carried by the step,
      ! written to ./output/element_flux_profile.txt and replaced at every
      ! call, so that after a run the file holds the state the run ended on.
      ! This is the T8 diagnostic of docs/binary_diffusion_design.md
      ! section 6, and the measurement the Phase-E flux closure reads
      ! (docs/phase_e_flux_closure_design.md section 3.4):
      !
      !   F_He(r_f) = 4 pi r_f^2 ( rho X v + J )         [g/s]
      !   F_H (r_f) = 4 pi r_f^2 ( rho (1-X) v - J )     [g/s]
      !
      ! A binary mixture has ONE independent diffusive flux, so the hydrogen
      ! element carries -J against the helium element's +J; the two elemental
      ! fluxes are written from the same X and the same J the step used, and
      ! neither is reconstructed anywhere else.  The advective part is
      ! rho_f v_f X_upwind (the physical face flux -- the operator carries the
      ! advection as a cell-velocity upwind difference, which has no face
      ! representation).  The total mass flux 4 pi r_f^2 rho_f v_f is written
      ! beside them: at a steady state all three are constant with radius, and
      ! the comparison of their radial spreads is the pass criterion.  dmeff
      ! is the face-averaged relative settling mass of the ambipolar
      ! diagnostic (3 neutral, 2.5 H+ plasma, 5/3 He++ plasma).
      !
      ! With a lower-atmosphere profile in use the routine also reduces each
      ! elemental flux to ONE number -- its median and its relative radial
      ! spread -- over each of TWO radial windows, for write_resolved_config
      ! and the flux closure of docs/phase_e_flux_closure_design.md section 6.
      !
      ! (a) The OVERLAP window, the interval both models describe.  Its upper
      !     edge is the radius the profile reaches, r(p_top).  Its lower edge
      !     is MEASURED rather than assumed: the first face at which the
      !     5-face moving spread of r^2 rho v falls below 10 per cent, which
      !     is where the standing base sound wave stops dominating
      !     (docs/binary_diffusion_design.md section 7.3, where that spread is
      !     10^2-10^4 times its own median).  When no face qualifies the edge
      !     falls back to 1.02 R_p and lap_flux_r_lo_measured says so.  An
      !     empty window is recorded as empty and never replaced by a number
      !     from outside it.
      !
      ! (b) The STEADY-FLUX window, r >= r_esc, i.e. the [j_min:N] escape
      !     region the solver uses to declare the wind steady.  At a steady
      !     state the elemental flux 4 pi r^2 (rho X v + J) is independent of
      !     radius, so the flux measured there IS the flux through the
      !     matching level.  That identity is why this window is a legitimate
      !     statement about the handoff and not a different quantity, and it
      !     is the only window wide enough on a planet whose lower-atmosphere
      !     column spans ~0.01 R_p.
      !
      ! The last three columns are the stage-resolved friction of memo 2.6:
      ! the effective coefficient D_eff the step used, the all-neutral
      ! hard-sphere coefficient of the same cell, their ratio (which is the
      ! Coulomb suppression where both elements are ionized), and the name of
      ! the stage pair carrying the largest share of the friction.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rep, Xhe, rho_phys, v
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dmeff, Dco, Dneut
      integer, dimension(1-Ng:N+Ng),intent(in) :: idom
      real*8, dimension(0:N),       intent(in) :: Jf
      real*8  :: area, adv, adv_H, rhof, vf, dmf, Df, Dnf, Xf
      real*8, dimension(1:N) :: FH_win, FHe_win, M_win, r_win
      integer :: j, uu, is, it, nwin

      nwin = 0

      uu = 771
      open(unit=uu, file='./output/element_flux_profile.txt',             &
           status='replace')
      write(uu,'(A)') '# elemental face fluxes carried by '//             &
         'binary_element_diffusion (T8 / Phase-E closure measurement)'
      write(uu,'(A)') '#   F_He   = 4 pi r^2 (rho X v + J)      [g/s]'
      write(uu,'(A)') '#   F_H    = 4 pi r^2 (rho (1-X) v - J)  [g/s]'
      write(uu,'(A)') '#   Mdot_f = 4 pi r^2 rho v              [g/s]'//  &
         '   (= F_H + F_He identically)'
      ! One whitespace-free token per column, so the schema line can be read
      ! the way every other EXHALE product's is.
      write(uu,'(A)') '# columns: j r_face[R_p] X_face J_diff'//           &
         '[g/cm2/s] F_adv[g/cm2/s] F_He[g/s] F_H[g/s]'//                  &
         ' Mdot_face[g/s] dmeff_face'//                                   &
         ' D_eff[cm2/s] D_neutral[cm2/s] D_eff/D_neutral dominant_pair'
      do j = 1, N-1
         area = 4.0d0*pi*rep(j)**2
         rhof = 0.5d0*(rho_phys(j) + rho_phys(j+1))
         vf   = 0.5d0*(v(j) + v(j+1))*v0
         if (vf .ge. 0.0d0) then
            Xf = Xhe(j)
         else
            Xf = Xhe(j+1)
         endif
         adv   = rhof*vf*Xf
         adv_H = rhof*vf*(1.0d0 - Xf)
         dmf  = 0.5d0*(dmeff(j) + dmeff(j+1))
         Df   = 0.5d0*(Dco(j)   + Dco(j+1))
         Dnf  = 0.5d0*(Dneut(j) + Dneut(j+1))
         is   = (idom(j) - 1)/n_hcar + 1
         it   = idom(j) - (is-1)*n_hcar
         write(uu,'(I6,11ES16.7,2X,A)') j, rep(j)/R0,                     &
              0.5d0*(Xhe(j)+Xhe(j+1)),                                    &
              Jf(j), adv, area*(adv + Jf(j)),                             &
              area*(adv_H - Jf(j)), area*rhof*vf, dmf,                    &
              Df, Dnf, Df/max(Dnf, 1.0d-99),                              &
              trim(hecar_name(is))//'-'//trim(hcar_name(it))
         if (lap_in_use) then
            nwin = nwin + 1
            r_win(nwin)   = rep(j)
            FHe_win(nwin) = area*(adv + Jf(j))
            FH_win(nwin)  = area*(adv_H - Jf(j))
            M_win(nwin)   = area*rhof*vf
         endif
      enddo
      close(uu)

      if (lap_in_use) call reduce_element_flux_windows(nwin, r_win, FH_win,&
                                                       FHe_win, M_win)

      end subroutine write_element_flux_profile

      ! ------------------------------------------------------------------ !

      subroutine reduce_element_flux_windows(nf, r_f, FH_f, FHe_f, M_f)
      ! Reduce the face-resolved elemental fluxes to one median and one
      ! relative spread over each of the two windows described above, and
      ! leave them in the lower_atmosphere_profile module for
      ! write_resolved_config.  Nothing here recomputes a flux: it only
      ! selects faces and takes order statistics of what the step carried.
      integer, intent(in) :: nf
      real*8, dimension(nf), intent(in) :: r_f, FH_f, FHe_f, M_f
      real*8, dimension(nf) :: a
      real*8  :: r_lo, r_hi, sprd, med
      integer :: i, i1, i2, n
      ! Half-width of the moving window the base-wave test is taken over, and
      ! the spread below which r^2 rho v is called free of the standing wave.
      integer, parameter :: nhalf = 2
      real*8,  parameter :: spread_flat = 0.10d0

      lap_flux_measured = .true.

      ! ---- (a) overlap window ------------------------------------------- !
      ! Lower edge: the first face whose 5-face moving spread of the mass flux
      ! is below spread_flat.  r^2 rho v is written here as M_f/(4 pi), and a
      ! ratio of a max-minus-min to a median is insensitive to that constant,
      ! so M_f is used directly.
      r_lo = -1.0d0
      lap_flux_r_lo_measured = .false.
      do i = 1 + nhalf, nf - nhalf
         med  = median_of(M_f(i-nhalf:i+nhalf))
         sprd = (maxval(M_f(i-nhalf:i+nhalf))                             &
                 - minval(M_f(i-nhalf:i+nhalf)))                          &
                /max(abs(med), 1.0d-99)
         if (sprd .lt. spread_flat) then
            r_lo = r_f(i)
            lap_flux_r_lo_measured = .true.
            exit
         endif
      enddo
      if (.not. lap_flux_r_lo_measured) r_lo = 1.02d0*R0

      r_hi = -1.0d0
      if (lap_r_top_RJ .gt. 0.0d0) r_hi = lap_r_top_RJ*RJ

      lap_flux_r_lo_Rp = r_lo/R0
      lap_flux_r_hi_Rp = r_hi/R0

      i1 = 0; i2 = -1
      do i = 1, nf
         if (r_f(i) .ge. r_lo .and. r_f(i) .le. r_hi) then
            if (i1 .eq. 0) i1 = i
            i2 = i
         endif
      enddo
      n = 0
      if (i1 .gt. 0) n = i2 - i1 + 1
      lap_flux_nface        = n
      lap_flux_window_empty = (n .lt. 1)
      if (n .ge. 1) then
         a(1:n) = FH_f(i1:i2)
         call median_and_spread(a(1:n), lap_FH_median,  lap_FH_spread)
         a(1:n) = FHe_f(i1:i2)
         call median_and_spread(a(1:n), lap_FHe_median, lap_FHe_spread)
         a(1:n) = M_f(i1:i2)
         call median_and_spread(a(1:n), lap_Mdot_median, lap_Mdot_spread)
      endif

      ! ---- (b) steady-flux window, r >= r_esc --------------------------- !
      i1 = 0; i2 = -1
      do i = 1, nf
         if (r_f(i) .ge. r_esc*R0) then
            if (i1 .eq. 0) i1 = i
            i2 = i
         endif
      enddo
      n = 0
      if (i1 .gt. 0) n = i2 - i1 + 1
      lap_steady_nface   = n
      lap_steady_r_lo_Rp = r_esc
      if (n .ge. 1) then
         a(1:n) = FH_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_FH_median,             &
                                        lap_steady_FH_spread)
         a(1:n) = FHe_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_FHe_median,            &
                                        lap_steady_FHe_spread)
         a(1:n) = M_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_Mdot_median,           &
                                        lap_steady_Mdot_spread)
      endif

      end subroutine reduce_element_flux_windows

      ! ------------------------------------------------------------------ !

      subroutine median_and_spread(a, med, sprd)
      ! Median and relative radial spread (max - min)/|median| of one window.
      real*8, dimension(:), intent(in)  :: a
      real*8,               intent(out) :: med, sprd
      med  = median_of(a)
      sprd = (maxval(a) - minval(a))/max(abs(med), 1.0d-99)
      end subroutine median_and_spread

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

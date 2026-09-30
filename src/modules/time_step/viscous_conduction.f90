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
      ! discrepancy in md/viscosity_conduction.md.
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
      ! THE CONDUCTIVITY OF THE MIXTURE (2026-09-30; replaces (4) in the
      ! heat conduction term, (5) is kept for the viscosity).  The heat
      ! conduction coefficient is the number-weighted sum of the
      ! Banks & Kockarts (1973) coefficients of the gas constituents, the
      ! form Salz et al. (2015, A&A 576, A21, Eqs. 24-26) use for an
      ! electron-proton plasma and neutral hydrogen and Sutton et al.
      ! (2015, JGR 120, 6884, Eq. 5) for a neutral mixture including
      ! helium:
      !
      !   kappa = (1/n) [ n_e kappa_ei + n_HI kappa_H + n_HeI kappa_He ]
      !   kappa_ei = 1.2e-6 T^(5/2)   (Salz 2015 Eq. 24, Spitzer 1978)
      !   kappa_H  = 379 T^0.69       (Salz 2015 Eq. 25)
      !   kappa_He = 299 T^0.69       (Sutton 2015 Eq. 5)
      !                                        erg cm^-1 s^-1 K^-1 ,   (4')
      !
      ! n the number density of all particles (electrons, atoms and ions,
      ! the metals included), one temperature for all of them.  The ions
      ! contribute nothing of their own (their coefficient, Yelle 2004
      ! Eq. 9, is 7.4e-8 T^(5/2)/sqrt(mu), below 1e-2 of kappa_H at 1e4 K);
      ! they enter only through n.  kappa_H is (4) to within 2 per cent
      ! from 300 to 7000 K (0.06 per cent at 1000 K).  At He/H = 9.7 the
      ! helium term decides kappa (0.79 kappa_H at 1000 K).  The helium
      ! coefficient agrees with the tabulated helium conductivity of
      ! Incropera et al. (2007, Fundamentals of Heat and Mass Transfer, 6th
      ! ed., Table A.4) to 1 per cent at 300 and 1000 K.
      !
      ! MOLECULAR HYDROGEN.  None of the three sources gives an H2
      ! coefficient; H2 enters (4') as n_H2 kappa_H2 with
      !
      !   kappa_H2 = A T^s [ 1 + c E(theta_v/T) ],  E(x) = x^2 e^x/(e^x - 1)^2,
      !   A = 272.523, s = 0.735652, c = 0.370446, theta_v = 5987 K
      !                                        erg cm^-1 s^-1 K^-1,   (4'')
      !
      ! the Eucken form (a viscosity T^s times a heat capacity whose
      ! vibrational part is the Einstein function of the v = 0 -> 1 spacing
      ! of H2, 5987.000 K in the Roueff et al. 2019 ladder this code carries,
      ! h2_vibrational_relaxation / molecular_infrared_data) fitted to the
      ! tabulated conductivity of Incropera et al. (2007), Table A.4, over
      ! 200-2000 K: largest deviation 2.5 per cent (at 200 K), rms 1.05 per
      ! cent (h2_conductivity).  The table is not a power law (local
      ! exponent 0.74 from 300 to 1000 K, 0.98 from 1900 to 2000 K); a
      ! single power law misses it by 6 per cent.  VALIDITY: 200-2000 K.
      ! Below 200 K the rotational heat capacity freezes out and (4'')
      ! overestimates (+7.6 per cent at 150 K, +20 per cent at 100 K); above
      ! 2000 K it is an EXTRAPOLATION of the fitted form (no tabulated
      ! value).  E is evaluated as x^2 e^-x/(1 - e^-x)^2, which cannot
      ! overflow at large x (low T), and by its Taylor series in x below
      ! x = 0.05, where 1 - e^-x loses digits (einstein_heat_capacity).  The
      ! molecular ions (H2+, H3+, HeH+) are trace particles and enter only
      ! through n and n_e.
      !
      ! VALIDITY.  (5) is the NEUTRAL atomic-hydrogen viscosity; Coulomb
      ! collisions raise the ion viscosity above the ionization front and
      ! that regime is NOT covered by (5).  The conductivity is the mixture
      ! value (4') since 2026-09-30, which carries the electron (Spitzer)
      ! term, the neutral H and He terms and the H2 term of (4''), so the
      ! ionized wind IS covered by it.  The terms are included for the dense,
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
      !   kappa_code = kappa_cgs T0/(rho0 v0^3 R0)
      !              (= (15/4) mu_code only for the Watson value (4) of a
      !              neutral atomic gas; not for the mixture value (4')).
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
      !   "Conduction: True"    -> cond_on, the mixture kappa of (4'), with
      !                            EXHALE_CONDUCTION_SCALE as its
      !                            continuation factor
      !
      ! -------------------------------------------------------------------
      ! 5. What the transport stage returns, and what it refuses to return
      ! -------------------------------------------------------------------
      ! The Crank-Nicolson stage returns a state only when every cell of it
      ! is admissible: the tridiagonal eliminations completed, every solved
      ! value is finite, and every solved temperature is at or above the
      ! floor of the range in which the equation of state and the transport
      ! coefficients are defined (T_floor, the same lower end the implicit
      ! energy source step brackets its root in, energy_semi_implicit).
      ! Otherwise the stage FAILS with a named status and reason, the cell,
      ! the value it solved for, and the energy the floor would have had to
      ! inject for the floor value to be a solution.
      !
      ! THE FLOOR IS NOT A RESERVOIR.  A temperature below T_floor means the
      ! diffusion operator asks this cell to give up more heat over dt than
      ! it holds down to the validity limit.  Returning the floor instead
      ! would add c_v (T_floor - T_solved) per volume to the state with no
      ! source term behind it, so the returned column would not be a zero of
      ! the steady residual that viscous_conduction_sources assembles from
      ! the same triplets, and the state would carry an unbudgeted accepted
      ! correction.
      !
      ! INTERIM ACTION ON FAILURE.  This module stops the run (error stop 1)
      ! after printing the diagnostics.  Step B3a replaces that stop by a
      ! REJECTION of the whole attempted step: the controller restores the
      ! checkpoint and retries with a shorter step, and this routine will
      ! then only report the status.  Until it exists no recovery is
      ! claimed here and stopping is the stop-safe action.
      !
      ! CONSERVATION.  The two operators are assembled in flux form: the
      ! coefficient of the face between cells j and j+1 is built from the
      ! same interpolated coefficient, the same face area and the same
      ! spacing on both sides, so sum_j dV_j Q_j telescopes to the net flux
      ! through the two ends of the column.  The outer end carries zero
      ! diffusive flux by construction, so a column whose base face also
      ! carries none conserves its thermal energy exactly under the
      ! conduction operator.  As boundary-conditioned in a RUN the base face
      ! does not carry none: the ghost is a Dirichlet anchor, and the heat
      ! that crosses that face is precisely the exchange with the anchored
      ! lower atmosphere this term exists to represent.  The stage is
      ! therefore conservative up to that named boundary flux, and closed
      ! only when the base flux vanishes.  Both statements have been verified.

      use global_parameters
      use Conversion, only: U_to_W
      use caloric_eos, only: caloric_mixture_active,                    &
                             energy_density_from_pressure,              &
                             heat_capacity_per_particle
      ! The lower end of the range in which the equation of state and the
      ! rate fits are defined, declared once for the whole code by the
      ! implicit energy source step and used here as the same floor.
      use energy_semi_implicit, only: T_eos_floor_K, T_floor_code_min
      ! The lower atmosphere's own temperature at a radius, which is the
      ! base boundary condition of this operator (conduction_base_level_T).
      use base_boundary, only: base_reservoir_temperature_at
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,      &
                               isp_HeIII, isp_H2, isp_H2p, isp_H3p,       &
                               isp_HeHp, isp_OH, isp_H2O, isp_CO,         &
                               n_mion, mion_fsp
      use utils, only: calc_ne, calc_ntot

      implicit none
      private
      public :: viscosity_active, conduction_active, transport_active
      public :: dynamic_viscosity, thermal_conductivity
      public :: viscous_momentum_source, viscous_dissipation
      public :: thermal_conduction_source, viscous_conduction_sources
      public :: thermal_conduction_coeffs
      public :: viscous_conduction_step
      public :: conduction_base_level_T, conduction_base_heat_flux
      public :: conduction_scale, h2_conductivity, einstein_heat_capacity

      ! THE BUDGET ENTRY OF THE BASE BOUNDARY CONDITION of this operator:
      ! the conductive heat flux through the base face at the last
      ! evaluation of thermal_conduction_source, code units, positive
      ! inward.  Zero bulk velocity does not make it zero, which is why it
      ! is measured and not assumed.
      real*8 :: conduction_base_heat_flux = 0.0d0
      public :: n_conduction_floor_hits, n_conduction_floor_cells
      public :: n_conduction_floor_hits_family
      public :: conduction_floor_first_step, conduction_floor_last_step
      ! THE CELL-BY-CELL FLOOR RECORD, readable the way
      ! energy_semi_implicit's energy_floor_cell_hits is, so that the two
      ! temperature floors of the marching loop report their attempts in
      ! the same shape and a reader of one needs no second convention for
      ! the other.  It is a diagnostic extremum over the run: unallocated
      ! until the floor is first reached, and one count per activation of
      ! the cell.
      public :: conduction_floor_cell_hits
      public :: CONDUCTION_OK, CONDUCTION_FLOOR, CONDUCTION_NONFINITE
      public :: CONDUCTION_SOLVE_FAILED
      public :: CONDUCTION_REASON_NONE, CONDUCTION_REASON_FLOOR_REACHED
      public :: CONDUCTION_REASON_NONFINITE_INPUT
      public :: CONDUCTION_REASON_NONFINITE_SOLUTION
      public :: CONDUCTION_REASON_PIVOT_BREAKDOWN
      public :: conduction_last_status, conduction_last_reason
      public :: conduction_last_cell, conduction_last_step
      public :: conduction_last_T_solution, conduction_last_T_floor
      public :: conduction_last_energy_deficit, conduction_last_nfail
      public :: conduction_stop_on_failure
      public :: conduction_temperature_floor

      ! Status of the Crank-Nicolson transport stage, reported through the
      ! optional `status` argument of viscous_conduction_step and left in
      ! conduction_last_status, so the physical-step context and the B3a
      ! controller read a verdict instead of inferring one from the state.
      integer, parameter :: CONDUCTION_OK           = 0
      integer, parameter :: CONDUCTION_FLOOR        = 1
      integer, parameter :: CONDUCTION_NONFINITE    = 2
      integer, parameter :: CONDUCTION_SOLVE_FAILED = 3

      ! Reasons, one per way the stage can fail. NONFINITE_SOLUTION and
      ! PIVOT_BREAKDOWN share the SOLVE_FAILED status and say different
      ! things about the matrix: a zero or non-finite pivot means the
      ! elimination broke down, a non-finite solution means it completed on
      ! a matrix that was not diagonally usable.
      integer, parameter :: CONDUCTION_REASON_NONE               = 0
      integer, parameter :: CONDUCTION_REASON_FLOOR_REACHED      = 1
      integer, parameter :: CONDUCTION_REASON_NONFINITE_INPUT    = 2
      integer, parameter :: CONDUCTION_REASON_NONFINITE_SOLUTION = 3
      integer, parameter :: CONDUCTION_REASON_PIVOT_BREAKDOWN    = 4

      ! Set .false. ONLY by the test drivers that exercise the failure
      ! statuses on purpose. Production runs stop.
      logical, save :: conduction_stop_on_failure = .true.

      ! Verdict of the last call, for a caller that reads it after the fact
      ! (EXHALE_main does not pass `status` today).
      !   T_solution        the temperature the operator solved for in the
      !                     failing cell, code units, unclamped
      !   T_floor           the floor it fell below, code units
      !   energy_deficit    c_v (T_floor - T_solution) per volume, erg/cm3:
      !                     the energy a floor would have injected into that
      !                     cell with no source term behind it
      integer, save :: conduction_last_status = CONDUCTION_OK
      integer, save :: conduction_last_reason = CONDUCTION_REASON_NONE
      integer, save :: conduction_last_cell   = 0
      integer, save :: conduction_last_step   = -1
      integer, save :: conduction_last_nfail  = 0
      real*8,  save :: conduction_last_T_solution     = 0.0d0
      real*8,  save :: conduction_last_T_floor        = 0.0d0
      real*8,  save :: conduction_last_energy_deficit = 0.0d0

      ! ATTEMPTS at the temperature floor of the Crank-Nicolson transport
      ! stage. A cell whose solved temperature falls below the floor is a
      ! FAILURE of the stage (see section 5 of the header): no state is
      ! built from the floor value, so these count attempted steps that were
      ! refused, never corrections carried by an adopted state. They are
      ! attempt statistics and are kept, not restored, on a rejected step. The
      ! implicit energy source step counts its own floor in the same terms.
      !   hits        total activations over the run
      !   first/last  first and last marching step on which the floor was hit
      !   cell_hits   activations of each cell, so the number of DISTINCT
      !               cells that ever reached the floor can be reported
      integer, save :: n_conduction_floor_hits = 0
      ! The same activations split by ledger family: index
      ! ledger_family_init counts the activations taken while the run was
      ! reaching a state, index ledger_family_phys those taken inside
      ! accepted physical steps. The total above is their sum.
      integer, save :: n_conduction_floor_hits_family(2) = 0
      integer, save :: conduction_floor_first_step = -1
      integer, save :: conduction_floor_last_step  = -1
      integer, allocatable, save :: conduction_floor_cell_hits(:)

      ! Watson, Donahue & Walker (1981) atomic-hydrogen heat conduction,
      ! kappa = kappa_1000K (T/1000 K)^kappa_expo  [erg cm^-1 s^-1 K^-1].
      real*8, parameter :: kappa_1000K = 4.45d4
      real*8, parameter :: kappa_expo  = 0.7d0
      ! Banks & Kockarts (1973) coefficients of (4'), erg cm^-1 s^-1 K^-1
      ! with T in K: electron-ion (Salz et al. 2015 Eq. 24), neutral
      ! hydrogen (Eq. 25) and neutral helium (Sutton et al. 2015 Eq. 5).
      real*8, parameter :: kappa_ei_coef = 1.2d-6, kappa_ei_expo = 2.5d0
      real*8, parameter :: kappa_H_coef  = 379.0d0
      real*8, parameter :: kappa_He_coef = 299.0d0
      real*8, parameter :: kappa_n_expo  = 0.69d0
      ! Molecular hydrogen, Eq. (4''): the fit to Incropera et al. (2007)
      ! Table A.4 over 200-2000 K (section 2), theta_v the v = 0 -> 1
      ! spacing of H2 (Roueff et al. 2019).
      real*8, parameter :: kappa_H2_A = 272.523d0, kappa_H2_s = 0.735652d0
      real*8, parameter :: kappa_H2_c = 0.370446d0, theta_v_H2 = 5987.0d0
      ! Below this x = theta_v/T the Einstein function is its Taylor series
      ! (einstein_heat_capacity).
      real*8, parameter :: einstein_series_x = 0.05d0
      ! Chapman-Enskog/Eucken ratio for a monatomic gas: kappa = eucken k_B mu/m.
      real*8, parameter :: eucken      = 3.75d0        ! = 15/4

      contains

      ! ------------------------------------------------------!

      ! Number of DISTINCT cells that have reached the temperature floor of
      ! the transport stage at least once.
      integer function n_conduction_floor_cells()
      integer :: j
      n_conduction_floor_cells = 0
      if (.not. allocated(conduction_floor_cell_hits)) return
      do j = lbound(conduction_floor_cell_hits,1),                       &
             ubound(conduction_floor_cell_hits,1)
         if (conduction_floor_cell_hits(j) .gt. 0)                       &
            n_conduction_floor_cells = n_conduction_floor_cells + 1
      enddo
      end function n_conduction_floor_cells

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

      real*8 function conduction_scale()
      ! A CONTINUATION FACTOR, NOT PHYSICS. Measurement key
      ! EXHALE_CONDUCTION_SCALE = s, 0 <= s <= 1 (1 when unset, which is
      ! the equation): the conductivity of (4') is multiplied by s, in the
      ! stationary row and the marching update alike. It exists to reach a
      ! state with conduction from one solved without it in steps: on the
      ! LHS 1140 b states the omitted term is 0.5 to 100 times the heating
      ! (md/Update_EXHALE_stage3.md section 40), and a solve that switches
      ! it on at once does not converge. A state solved at 0 < s < 1 is a
      ! step of that continuation and not a prediction of the model. The
      ! factor is part of the state identity: the 'cond' token of the
      ! options field carries its value wherever s /= 1 (load_IC,
      ! opt_value).
      real*8,  save :: scale = 1.0d0
      logical, save :: scale_read = .false.
      character(len=32) :: env
      integer :: st
      if (.not. scale_read) then
         scale_read = .true.
         call get_environment_variable('EXHALE_CONDUCTION_SCALE', env,     &
                                       status=st)
         if (st .eq. 0 .and. len_trim(env) .gt. 0) then
            read(env,*,iostat=st) scale
            if (st .ne. 0 .or. .not. (scale .ge. 0.0d0 .and.               &
                scale .le. 1.0d0)) then
               write(*,*) '(conduction) ERROR: EXHALE_CONDUCTION_SCALE'//  &
                          ' takes a number from 0 to 1, not "'//           &
                          trim(env)//'".'
               error stop 1
            endif
            write(*,'(A,ES13.6,A)') ' (conduction) heat conduction scaled'//&
                 ' by', scale, ' (continuation factor'//                   &
                 ' EXHALE_CONDUCTION_SCALE; 1 is the equation)'
         endif
      endif
      conduction_scale = scale
      end function conduction_scale

      ! ------------------------------------------------------!

      subroutine thermal_conductivity(Tcell, f_sp, kap)
      ! Heat conduction coefficient of the mixture, Eq. (4'), in CODE
      ! units, times the continuation factor conduction_scale. The number
      ! fractions come from the composition f_sp of the same cells (rho
      ! f_sp is a number density, so rho cancels in every ratio); the
      ! electron and particle counts are the equation of state's own
      ! (calc_ne, calc_ntot).
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng), intent(out) :: kap
      real*8, dimension(1-Ng:N+Ng) :: xe, xH, xHe, zero, nhe0, nhe1, nhe2
      real*8, dimension(1-Ng:N+Ng) :: ne_w, nt_w, nh2
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm_w
      real*8, dimension(1-Ng:N+Ng,4) :: nmol_w
      real*8, dimension(1-Ng:N+Ng,3) :: nox_w
      real*8 :: to_code, TK, kap_H2
      integer :: j, im
      if (.not. conduction_active()) then
         kap = 0.0d0
         return
      endif
      zero = 0.0d0
      if (thereis_He) then
         nhe0 = f_sp(:,isp_HeI);  nhe1 = f_sp(:,isp_HeII)
         nhe2 = f_sp(:,isp_HeIII)
      else
         nhe0 = 0.0d0;  nhe1 = 0.0d0;  nhe2 = 0.0d0
      endif
      do im = 1, n_mion
         nm_w(:,im) = f_sp(:,mion_fsp(im))
      enddo
      if (thereis_mol) then
         nmol_w(:,1) = f_sp(:,isp_H2);   nmol_w(:,2) = f_sp(:,isp_H2p)
         nmol_w(:,3) = f_sp(:,isp_H3p);  nmol_w(:,4) = f_sp(:,isp_HeHp)
         call calc_ne(f_sp(:,isp_HII), nhe1, nhe2, ne_w, nm_w, nmol_w)
         if (thereis_oxychem) then
            nox_w(:,1) = f_sp(:,isp_OH);  nox_w(:,2) = f_sp(:,isp_H2O)
            nox_w(:,3) = f_sp(:,isp_CO)
            call calc_ntot(f_sp(:,isp_HI), f_sp(:,isp_HII), nhe0, nhe1,  &
                           nhe2, nt_w, nm_w, nmol_w, nox_w)
         else
            call calc_ntot(f_sp(:,isp_HI), f_sp(:,isp_HII), nhe0, nhe1,  &
                           nhe2, nt_w, nm_w, nmol_w)
         endif
         nh2 = f_sp(:,isp_H2)
      else
         call calc_ne(f_sp(:,isp_HII), nhe1, nhe2, ne_w, nm_w)
         call calc_ntot(f_sp(:,isp_HI), f_sp(:,isp_HII), nhe0, nhe1,      &
                        nhe2, nt_w, nm_w)
         nh2 = 0.0d0
      endif
      ! cgs -> code: kappa_code = kappa_cgs T0/(rho0 v0^3 R0)
      to_code = conduction_scale()*T0/(n0*mu*v0**3*R0)
      do j = 1-Ng, N+Ng
         TK = max(Tcell(j), 1.0d-8)*T0
         xe(j)  = ne_w(j)/(ne_w(j) + nt_w(j))
         xH(j)  = f_sp(j,isp_HI)/(ne_w(j) + nt_w(j))
         xHe(j) = nhe0(j)/(ne_w(j) + nt_w(j))
         ! The H2 term only where there is H2, so that its coefficient is
         ! never evaluated for a mixture that does not carry it.
         kap_H2 = 0.0d0
         if (nh2(j) .gt. 0.0d0)                                           &
            kap_H2 = nh2(j)/(ne_w(j) + nt_w(j))*h2_conductivity(TK)
         kap(j) = to_code*( xe(j)*kappa_ei_coef*TK**kappa_ei_expo         &
                          + xH(j)*kappa_H_coef*TK**kappa_n_expo           &
                          + xHe(j)*kappa_He_coef*TK**kappa_n_expo         &
                          + kap_H2)
      enddo
      end subroutine thermal_conductivity

      ! ------------------------------------------------------!

      real*8 function h2_conductivity(TK)
      ! Thermal conductivity of molecular hydrogen [erg cm^-1 s^-1 K^-1] at
      ! TK [K], Eq. (4''): valid 200-2000 K, an extrapolation outside
      ! (section 2 of the header).
      real*8, intent(in) :: TK
      if (.not. (TK .gt. 0.0d0 .and. TK .le. huge(TK))) then
         write(*,'(A,ES24.16)') ' (conduction) ERROR: H2 conductivity'//   &
              ' asked for at a temperature that is not positive and'//     &
              ' finite, TK [K] =', TK
         error stop 1
      endif
      h2_conductivity = kappa_H2_A*TK**kappa_H2_s                          &
                        *(1.0d0 + kappa_H2_c                               &
                                  *einstein_heat_capacity(theta_v_H2/TK))
      end function h2_conductivity

      ! ------------------------------------------------------!

      pure real*8 function einstein_heat_capacity(x)
      ! Einstein function E(x) = x^2 e^x/(e^x - 1)^2 of Eq. (4''), the
      ! vibrational heat capacity of one oscillator in units of k_B, at
      ! x = theta_v/T >= 0.  Written with e^-x, x^2 e^-x/(1 - e^-x)^2, so
      ! that no factor overflows.  In the e^x form (e^x - 1)^2 overflows
      ! above x ~ 355 (T < 16.9 K for H2), where that form returned 0 in
      ! place of E < 1e-148, and x^2 e^x above x = 696.69, where it
      ! returned NaN (Inf/Inf) up to x = 700, i.e. for
      ! 8.553 K <= T < 8.593 K.  Beyond x = 700 E < 1e-298 and is set to
      ! zero, which keeps e^-x out of the subnormal range.  Below
      ! x = einstein_series_x = 0.05 the difference 1 - e^-x cancels
      ! (relative rounding ~ eps/x) and E is its Taylor series
      !   E = 1 - x^2/12 + x^4/240 - x^6/6048 + O(x^8/172800),
      ! whose truncation error there is below 2.3e-16 relative.  The two
      ! forms agree at the switch to 6.7e-16 relative (gfortran and ifx;
      ! test einstein_series_continuity, src/tests/energy_update).
      real*8, intent(in) :: x
      real*8 :: x2, em
      if (x .lt. einstein_series_x) then
         x2 = x*x
         einstein_heat_capacity = 1.0d0 - x2*(1.0d0/12.0d0                 &
                                  - x2*(1.0d0/240.0d0 - x2/6048.0d0))
      else if (x .gt. 700.0d0) then
         einstein_heat_capacity = 0.0d0
      else
         em = exp(-x)
         einstein_heat_capacity = x*x*em/(1.0d0 - em)**2
      endif
      end function einstein_heat_capacity

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

      real*8 function conduction_base_level_T(Tcell) result(Tb)
      ! THE BASE BOUNDARY CONDITION OF THERMAL CONDUCTION, stated here and
      ! not inherited from the advective ghost: a prescribed temperature at
      ! the base ghost cell, equal to the lower atmosphere's own temperature
      ! at that radius on the reservoir's hydrostatic isentrope.
      !
      ! WHY THE RESERVOIR AND NOT THE GHOST THE ADVECTION LEAVES.  The lower
      ! atmosphere below the base level is a heat bath on the time scales of
      ! a stationary solution (base_boundary.f90, "WHO OWNS THE LEVEL"), and
      ! a bath does not stop existing when the gas above it drains: the
      ! temperature it holds at the level is the same whichever way the
      ! contact is upwinded.  Reading the advective ghost instead would make
      ! this operator's boundary condition a silent function of the
      ! direction of the flow, and a reverse flow would leave the base
      ! adiabatic, which is not a statement anything in this model makes.
      !
      ! WHAT IT COSTS.  The ghost's own p/n_part and the reservoir's T at the
      ! level are not the same number, and the distance is small: the
      ! prescribed particles per unit mass and the count the composition
      ! sweep solves in the ghost differ by 2.7e-08 on the hot-Uranus
      ! carrier state and 1.3e-07 on the LHS 1140 b molecular state, and the
      ! ghost's temperature stands 9.1e-07 and 1.3e-04 from the level's at
      ! the level's own radius (MEASURED, with every quantity located at its
      ! own radius; the 4.7 per cent this comment once carried compared the
      ! count prescribed AT THE LEVEL with the count solved one cell below
      ! it, a difference of location that follows the grid).  Both
      ! counts are reported at base_boundary's report_base_boundary_model.
      ! This operator takes the reservoir's, because it is the bath's temperature
      ! that drives the heat flux across the level.
      !
      ! FALLBACK: a reservoir that returns no positive temperature leaves the
      ! ghost's own value standing, so a configuration without a stated base
      ! level conducts against the state it holds and not against a zero.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: Tcell
      Tb = base_reservoir_temperature_at(r(0))
      if (.not. (Tb .gt. 0.0d0)) Tb = Tcell(0)
      end function conduction_base_level_T

      ! ------------------------------------------------------!

      subroutine thermal_conduction_coeffs(Tcell, f_sp, blo, bdi, bup)
      ! Tridiagonal coefficients of (1/r^2) d/dr(r^2 kappa dT/dr):
      !   Q(j) = blo(j) T(j-1) + bdi(j) T(j) + bup(j) T(j+1) .
      !
      ! The conductivity at the base ghost is evaluated at the base level's
      ! own temperature (conduction_base_level_T), the same value the
      ! callers multiply blo(1) by, so the coefficient and the temperature
      ! it acts on belong to one boundary condition.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(N),         intent(out) :: blo, bdi, bup
      real*8, dimension(1-Ng:N+Ng) :: kap, Tb
      real*8 :: Ap, Am, dV, rp, rm, kp, km, wp, wm, dp, dm
      integer :: j

      Tb = Tcell
      Tb(1-Ng:0) = conduction_base_level_T(Tcell)
      call thermal_conductivity(Tb, f_sp, kap)
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

      subroutine thermal_conduction_source(Tcell, f_sp, Qc)
      ! Heat conduction source (1/r^2) d/dr(r^2 kappa dT/dr), code units.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Qc
      real*8, dimension(N) :: blo, bdi, bup
      real*8  :: Tb
      integer :: j
      Qc = 0.0d0
      conduction_base_heat_flux = 0.0d0
      if (.not. conduction_active()) return
      call thermal_conduction_coeffs(Tcell, f_sp, blo, bdi, bup)
      Tb = conduction_base_level_T(Tcell)
      do j = 1, N
         if (j .eq. 1) then
            Qc(j) = blo(j)*Tb + bdi(j)*Tcell(j) + bup(j)*Tcell(j+1)
         else
            Qc(j) = blo(j)*Tcell(j-1) + bdi(j)*Tcell(j) + bup(j)*Tcell(j+1)
         endif
      enddo
      ! The heat the lower atmosphere puts through the base face, per unit
      ! area and time, in code units: the budget entry of this operator's
      ! own boundary condition.  Positive is heat entering the domain.  It
      ! is the base face's own term of the row above, blo(1) (T_base - T_1),
      ! multiplied back by the cell volume and divided by the face area, so
      ! it is exactly the flux the row differences and not a second estimate
      ! of it.
      conduction_base_heat_flux = blo(1)*(Tb - Tcell(1))                  &
           *((r_edg(1)**3 - r_edg(0)**3)/3.0d0)/(r_edg(0)*r_edg(0))
      end subroutine thermal_conduction_source

      ! ------------------------------------------------------!

      subroutine viscous_conduction_sources(vel, Tcell, f_sp, Smom, Sene)
      ! The pair that enters the conserved-variable equations:
      !   Smom = F_mu                          (momentum, Eq. 1)
      !   Sene = w F_mu + q_mu + conduction    (TOTAL energy, Eq. 3)
      ! This is the single definition used by BOTH the steady residual and
      ! the marching update, so the two cannot describe different systems.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: vel, Tcell
      ! The composition of the same cells (the conductivity of the mixture).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Smom, Sene
      real*8, dimension(1-Ng:N+Ng) :: qv, Qc
      Smom = 0.0d0;  Sene = 0.0d0
      if (.not. transport_active()) return
      call viscous_momentum_source(vel, Tcell, Smom)
      call viscous_dissipation(vel, Tcell, qv)
      call thermal_conduction_source(Tcell, f_sp, Qc)
      Sene = vel*Smom + qv + Qc
      end subroutine viscous_conduction_sources

      ! ------------------------------------------------------!

      logical function is_finite(x)
      ! .true. for a normal or subnormal number: neither a NaN (x /= x) nor
      ! an infinity.
      real*8, intent(in) :: x
      is_finite = (x .eq. x) .and. (abs(x) .le. huge(x))
      end function is_finite

      ! ------------------------------------------------------!

      subroutine solve_tridiagonal(dlo, ddi, dup, rhs, sol, ok)
      ! Thomas algorithm for the N x N tridiagonal system
      !   dlo(j) x(j-1) + ddi(j) x(j) + dup(j) x(j+1) = rhs(j).
      ! Radial recursion: inherently serial in j, so this must stay OUTSIDE
      ! any cell-parallel region.
      ! ok is .false. if an elimination pivot is zero or not finite, which
      ! is the breakdown of the algorithm: the caller must not read sol.
      real*8, dimension(N), intent(in)  :: dlo, ddi, dup, rhs
      real*8, dimension(N), intent(out) :: sol
      logical,              intent(out) :: ok
      real*8, dimension(N) :: cp, dp
      real*8  :: den
      integer :: j
      ok = .true.
      den   = ddi(1)
      if (den .eq. 0.0d0 .or. .not. is_finite(den)) then
         ok = .false.;  sol = 0.0d0;  return
      endif
      cp(1) = dup(1)/den
      dp(1) = rhs(1)/den
      do j = 2, N
         den   = ddi(j) - dlo(j)*cp(j-1)
         if (den .eq. 0.0d0 .or. .not. is_finite(den)) then
            ok = .false.;  sol = 0.0d0;  return
         endif
         cp(j) = dup(j)/den
         dp(j) = (rhs(j) - dlo(j)*dp(j-1))/den
      enddo
      sol(N) = dp(N)
      do j = N-1, 1, -1
         sol(j) = dp(j) - cp(j)*sol(j+1)
      enddo
      end subroutine solve_tridiagonal

      ! ------------------------------------------------------!

      real*8 function conduction_temperature_floor()
      ! Lower end of the range in which the equation of state, the cooling
      ! and the transport coefficients are defined, in CODE units. The same
      ! value the implicit energy source step brackets its root in: 0.01 T0,
      ! raised to the 1 K lower node of the H2 rovibrational table if T0 is
      ! low enough for 0.01 T0 to fall below it.
      conduction_temperature_floor = max(T_floor_code_min, T_eos_floor_K/T0)
      end function conduction_temperature_floor

      ! ------------------------------------------------------!

      subroutine conduction_reason_text(reason_in, text)
      integer, intent(in) :: reason_in
      character(len=*), intent(out) :: text
      select case (reason_in)
      case (CONDUCTION_REASON_FLOOR_REACHED)
         text = 'transport step asks for a temperature below the floor'
      case (CONDUCTION_REASON_NONFINITE_INPUT)
         text = 'non-finite state, particle count or step on entry'
      case (CONDUCTION_REASON_NONFINITE_SOLUTION)
         text = 'non-finite solution of the Crank-Nicolson system'
      case (CONDUCTION_REASON_PIVOT_BREAKDOWN)
         text = 'zero or non-finite pivot in the tridiagonal elimination'
      case default
         text = 'admissible'
      end select
      end subroutine conduction_reason_text

      ! ------------------------------------------------------!

      subroutine conduction_report(stat, reason, j_bad, nfail, step,      &
                                   T_sol, T_flr, deficit, status)
      ! Record the verdict of one call, hand it to the caller, and take the
      ! interim action on a failure (print the diagnostics and stop). B3a
      ! replaces the stop by a rejection of the whole attempted step.
      integer, intent(in) :: stat, reason, j_bad, nfail, step
      real*8,  intent(in) :: T_sol, T_flr, deficit
      integer, intent(out), optional :: status
      character(len=72) :: reason_text

      conduction_last_status         = stat
      conduction_last_reason         = reason
      conduction_last_cell           = j_bad
      conduction_last_step           = step
      conduction_last_nfail          = nfail
      conduction_last_T_solution     = T_sol
      conduction_last_T_floor        = T_flr
      conduction_last_energy_deficit = deficit
      if (present(status)) status = stat
      if (stat .eq. CONDUCTION_OK) return

      call conduction_reason_text(reason, reason_text)
      write(*,*)
      write(*,'(a)')      ' transport stage FAILURE (viscosity/conduction)'
      write(*,'(a,i0,a,i0,a)') '   step ', step, ', ', nfail,             &
                               ' cell(s) unacceptable'
      write(*,'(a,i0)')   '   first failing cell   j = ', j_bad
      write(*,'(a,a)')    '   reason               ', trim(reason_text)
      write(*,'(a,i0)')   '   status               ', stat
      if (stat .eq. CONDUCTION_FLOOR) then
         write(*,'(a,es13.6)') '   solved T         [K] ', T_sol*T0
         write(*,'(a,es13.6)') '   floor            [K] ', T_flr*T0
         write(*,'(a,es13.6)') '   energy a floor would inject [erg/cm3] ', &
                               deficit
         write(*,'(a)') '   That energy has no source term behind it, so'
         write(*,'(a)') '   the floor value is not returned as a state.'
      endif
      write(*,'(a)') '   The transport stage returns no state for this step.'
      if (conduction_stop_on_failure) then
         write(*,'(a)') '   Stopping: B3a will reject and retry the step'
         write(*,'(a)') '   instead.'
         error stop 1
      endif
      end subroutine conduction_report

      ! ------------------------------------------------------!

      subroutine viscous_conduction_step(u, W, Tcell, n_part, f_sp, dt,  &
                                         step, status)
      ! Crank-Nicolson update of the two transport operators, applied as an
      ! operator-split stage of the marching loop, as CETIMB integrates the
      ! same terms.  Both are diffusive, so an explicit update would be bound
      ! by dt < rho dr^2/mu (and C dr^2/kappa); the implicit form removes that
      ! bound unconditionally.  MEASURED at the WASP-121 b resolution, though,
      ! the diffusive limit is ~1 t_s against an advective step of ~9e-5 t_s,
      ! so these terms are NOT stiff there and the implicit treatment is
      ! insurance for a denser base, a coarser grid or a larger coefficient
      ! rather than a necessity (md/viscosity_conduction.md Sec. 4).
      !
      !   (i)  momentum   rho (w* - w)/dt = (1/2)[F_mu(w) + F_mu(w*)]
      !        The total energy is advanced by the viscous WORK that goes
      !        with it, so the internal energy is untouched by this stage:
      !        d(rho w^2/2) = w_avg rho (w* - w) = dt w_avg F_mu_avg exactly,
      !        which is what holding p fixed and rebuilding E accomplishes.
      !   (ii) temperature  C (T* - T)/dt = (1/2)[Q(T) + Q(T*)] + q_mu,
      !        C = n_part c_v/k the internal energy per unit temperature,
      !        q_mu evaluated at the updated velocity.
      !
      ! The spatial operators are the SAME tridiagonal triplets that
      ! viscous_conduction_sources evaluates for the steady residual, so the
      ! fixed point of this update is exactly the zero of that residual.
      !
      ! The state is assembled ONLY when both solves completed and every
      ! solved temperature is admissible; a cell below the floor is a
      ! failure, not a clamp (section 5 of the module header).
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u, W
      real*8, dimension(1-Ng:N+Ng),   intent(in)    :: Tcell, n_part, dt
      ! The composition of the same cells (the conductivity of the mixture).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      ! Marching step index, for the failure report and the floor counters.
      integer,                        intent(in)    :: step
      ! Verdict of the stage. Optional, so the call site is unchanged; the
      ! same value is left in conduction_last_status.
      integer, intent(out), optional :: status
      real*8, dimension(1-Ng:N+Ng)   :: vnew, Tnew, qv
      real*8, dimension(N) :: alo, adi, aup, dcl, dcd, dcu
      real*8, dimension(N) :: blo, bdi, bup
      real*8, dimension(N) :: dlo, ddi, dup, rhs, sol
      real*8  :: half, Lold, cap, T_flr, c_T, deficit
      integer :: j, stat, reason, j_bad, nfail
      logical :: ok

      stat    = CONDUCTION_OK
      reason  = CONDUCTION_REASON_NONE
      j_bad   = 0
      nfail   = 0
      deficit = 0.0d0
      T_flr   = conduction_temperature_floor()
      if (.not. transport_active()) then
         call conduction_report(CONDUCTION_OK, CONDUCTION_REASON_NONE, 0, &
                                0, step, 0.0d0, T_flr, 0.0d0, status)
         return
      endif
      half = 0.5d0
      vnew = W(2,:)
      Tnew = Tcell

      ! ---- entry admissibility ----
      ! The stencils reach one cell past each end of 1..N; nothing beyond
      ! that is read, so nothing beyond that is tested.
      do j = 0, N+1
         if (is_finite(Tcell(j)) .and. is_finite(W(1,j))                  &
             .and. is_finite(W(2,j))) cycle
         nfail = nfail + 1
         if (j_bad .eq. 0) then
            j_bad  = j
            stat   = CONDUCTION_NONFINITE
            reason = CONDUCTION_REASON_NONFINITE_INPUT
         endif
      enddo
      do j = 1, N
         if (is_finite(n_part(j)) .and. is_finite(dt(j))                  &
             .and. dt(j) .gt. 0.0d0) cycle
         nfail = nfail + 1
         if (j_bad .eq. 0) then
            j_bad  = j
            stat   = CONDUCTION_NONFINITE
            reason = CONDUCTION_REASON_NONFINITE_INPUT
         endif
      enddo
      if (stat .ne. CONDUCTION_OK) then
         call conduction_report(stat, reason, j_bad, nfail, step,         &
                                0.0d0, T_flr, 0.0d0, status)
         return
      endif

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
         call solve_tridiagonal(dlo, ddi, dup, rhs, sol, ok)
         if (.not. ok) then
            call conduction_report(CONDUCTION_SOLVE_FAILED,               &
                                   CONDUCTION_REASON_PIVOT_BREAKDOWN, 1,  &
                                   1, step, 0.0d0, T_flr, 0.0d0, status)
            return
         endif
         do j = 1, N
            if (is_finite(sol(j))) cycle
            nfail = nfail + 1
            if (j_bad .eq. 0) then
               j_bad  = j
               stat   = CONDUCTION_SOLVE_FAILED
               reason = CONDUCTION_REASON_NONFINITE_SOLUTION
            endif
         enddo
         if (stat .ne. CONDUCTION_OK) then
            call conduction_report(stat, reason, j_bad, nfail, step,      &
                                   0.0d0, T_flr, 0.0d0, status)
            return
         endif
         do j = 1, N
            vnew(j) = sol(j)
         enddo
      endif

      ! ---- (ii) heat conduction + viscous dissipation ----
      if (conduction_active() .or. viscosity_active()) then
         ! The dissipation belongs to the UPDATED velocity field, which is
         ! vnew: it equals the incoming W(2,:) wherever the momentum stage
         ! did not run, and in the ghosts, which it never touches.
         call viscous_dissipation(vnew, Tcell, qv)
         if (conduction_active()) then
            call thermal_conduction_coeffs(Tcell, f_sp, blo, bdi, bup)
         else
            blo = 0.0d0;  bdi = 0.0d0;  bup = 0.0d0
         endif
         do j = 1, N
            if (j .eq. 1) then
               Lold = blo(j)*conduction_base_level_T(Tcell)              &
                      + bdi(j)*Tcell(j) + bup(j)*Tcell(j+1)
            else
               Lold = blo(j)*Tcell(j-1) + bdi(j)*Tcell(j)                &
                      + bup(j)*Tcell(j+1)
            endif
            ! Heat capacity of the cell, frozen at the stage temperature
            ! exactly as the conductivity above it is.  This is a
            ! linearization of the implicit step, not of the answer: at the
            ! fixed point Tnew = Tcell the cap terms cancel and what remains
            ! is Lold + qv = 0, which is the zero of
            ! viscous_conduction_sources whatever cap was.
            if (caloric_mixture_active) then
               cap = n_part(j)*heat_capacity_per_particle(j, Tcell(j))   &
                     /dt(j)
            else
            cap    = n_part(j)/((gamma_ad - 1.0d0)*dt(j))
            endif
            dlo(j) = -half*blo(j)
            ddi(j) = cap - half*bdi(j)
            dup(j) = -half*bup(j)
            rhs(j) = cap*Tcell(j) + half*Lold + qv(j)
         enddo
         ! The base row's Dirichlet value is this operator's own boundary
         ! condition (conduction_base_level_T), not the advective ghost.
         rhs(1) = rhs(1) - dlo(1)*conduction_base_level_T(Tcell)
         dlo(1) = 0.0d0
         dup(N) = 0.0d0
         call solve_tridiagonal(dlo, ddi, dup, rhs, sol, ok)
         if (.not. ok) then
            call conduction_report(CONDUCTION_SOLVE_FAILED,               &
                                   CONDUCTION_REASON_PIVOT_BREAKDOWN, 1,  &
                                   1, step, 0.0d0, T_flr, 0.0d0, status)
            return
         endif
         ! Admissibility of the solved temperatures. A cell at or below the
         ! floor is refused: the energy that would have to be injected to
         ! put it there, c_v (T_floor - T_solved) per volume, has no source
         ! term behind it. The code energy density unit is n0 k_B T0.
         do j = 1, N
            if (.not. is_finite(sol(j))) then
               nfail = nfail + 1
               if (j_bad .eq. 0) then
                  j_bad  = j
                  stat   = CONDUCTION_SOLVE_FAILED
                  reason = CONDUCTION_REASON_NONFINITE_SOLUTION
               endif
               cycle
            endif
            if (sol(j) .gt. T_flr) cycle
            nfail = nfail + 1
            if (.not. allocated(conduction_floor_cell_hits)) then
               allocate(conduction_floor_cell_hits(1-Ng:N+Ng))
               conduction_floor_cell_hits = 0
            endif
            n_conduction_floor_hits = n_conduction_floor_hits + 1
            n_conduction_floor_hits_family(ledger_family) =              &
                 n_conduction_floor_hits_family(ledger_family) + 1
            conduction_floor_cell_hits(j) =                              &
                 conduction_floor_cell_hits(j) + 1
            if (conduction_floor_first_step .lt. 0)                      &
               conduction_floor_first_step = step
            conduction_floor_last_step = step
            if (j_bad .eq. 0 .or. stat .eq. CONDUCTION_OK) then
               if (caloric_mixture_active) then
                  c_T = n_part(j)*heat_capacity_per_particle(j, Tcell(j))
               else
                  c_T = n_part(j)/(gamma_ad - 1.0d0)
               endif
               j_bad   = j
               stat    = CONDUCTION_FLOOR
               reason  = CONDUCTION_REASON_FLOOR_REACHED
               deficit = c_T*(T_flr - sol(j))*n0*kb_erg*T0
               conduction_last_T_solution = sol(j)
            endif
         enddo
         if (stat .ne. CONDUCTION_OK) then
            call conduction_report(stat, reason, j_bad, nfail, step,      &
                                   conduction_last_T_solution, T_flr,     &
                                   deficit, status)
            return
         endif
         do j = 1, N
            Tnew(j) = sol(j)
         enddo
      endif

      ! ---- assemble the returned state ----
      ! Reached only when every cell of both solves is admissible. The
      ! momentum stage rebuilds the conserved state at FIXED pressure: the
      ! kinetic energy change is precisely the viscous work, so no internal
      ! energy is created or destroyed by it, and the pressure row is then
      ! set by the temperature the second stage solved for.
      do j = 1, N
         if (viscosity_active()) then
            W(2,j) = vnew(j)
            u(2,j) = W(1,j)*vnew(j)
         endif
         W(3,j)  = n_part(j)*Tnew(j)
         u(3,j)  = 0.5d0*W(1,j)*W(2,j)**2                                &
                 + energy_density_from_pressure(j, W(1,j), W(3,j))
      enddo

      call conduction_report(CONDUCTION_OK, CONDUCTION_REASON_NONE, 0, 0, &
                             step, 0.0d0, T_flr, 0.0d0, status)

      end subroutine viscous_conduction_step

      ! End of module
      end module viscous_conduction

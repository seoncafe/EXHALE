      module caloric_eos
      ! Caloric equation of state of the gas: the map between the internal
      ! energy density, the pressure and the temperature.
      !
      ! WHAT IS WRONG WITH ONE CONSTANT gamma.  The thermal EOS of this code
      ! is already composition exact -- composition.f90 carries
      ! p = (n_tot + n_e) k T with the species count of the cell -- but the
      ! CALORIC half, e = p/(gamma - 1), assumed a monatomic gamma = 5/3
      ! everywhere.  That says every particle stores (3/2) kT.  Molecular
      ! hydrogen does not: at the 400-1300 K of the molecular layer its
      ! rotational ladder is fully excited (the J = 1 level sits at 170 K), so
      ! an H2 molecule stores about (5/2) kT, and the vibrational ladder adds
      ! to that above ~800 K.  Where the gas is 80 percent H2 the specific
      ! heat is therefore wrong by a factor 1.67.
      !
      ! WHAT THIS MODULE HOLDS.  For a mixture in which every atom, ion and
      ! free electron carries (3/2) k and H2 carries its own rovibrational
      ! heat capacity as well, the internal energy density is
      !
      !     e = (3/2) (n_tot + n_e) k T  +  n(H2) k u_rv(T)
      !
      ! with u_rv the mean rovibrational energy of an H2 molecule in units of
      ! k (so u_rv has the dimension of a temperature).  Because
      ! p = (n_tot + n_e) k T is unchanged, the exponent that makes
      ! e = p/(gamma_eff - 1) exact is
      !
      !     gamma_eff = 1 + (n_tot + n_e)/C_V_over_k ,
      !     C_V_over_k = (3/2)(n_tot + n_e) + n(H2) c_rv(T) ,
      !
      ! and gamma_eff is a function of the composition AND of T.  The module
      ! provides the four maps the hydrodynamics needs -- (rho, e) -> T,
      ! (rho, e) -> p, (rho, p) -> e and (rho, p) -> gamma_eff -- plus
      ! gamma_eff itself for the wave-speed estimates.
      !
      ! ONE H2 MOLECULE IS THE WHOLE MOLECULAR CORRECTION.  H2+, H3+ and HeH+
      ! are counted with (3/2) k each, i.e. as if monatomic, and so are the
      ! oxygen carriers OH, H2O and CO.  This is deliberate and it is an
      ! approximation: those species do have internal degrees of freedom, but
      ! in this code's molecular layer they are trace.  Measured on the
      ! converged He/H = 0.0793 hot Uranus, over every cell where H2 is more
      ! than 1e-3 of the particle count (r = 1.000 to 1.101 r_base), the
      ! largest ratios to H2 are n(H2+)/n(H2) = 1.7e-7, n(HeH+)/n(H2) = 1.2e-7
      ! and n(H3+)/n(H2) = 3.2e-8.  Giving all three the full H2 heat capacity
      ! instead of the monatomic one would move C_V by less than 1e-6, six
      ! orders below the H2 term this module exists for.  If a configuration
      ! ever makes one of them a major carrier, it belongs here and the
      ! assumption above stops holding.
      !
      ! WHERE THE LEVEL DATA COMES FROM.  u_rv and c_rv are evaluated from the
      ! COMPLETE bound rovibrational ladder of H2 X^1 Sigma_g^+ -- 302 levels
      ! up to 51966 K, i.e. up to the 4.478 eV dissociation limit -- which the
      ! code already carries in molecular_infrared_data (Roueff et al. 2019,
      ! A&A 630, A58, table 2, via VizieR J/A+A/630/A58): h2_lev_T is the term
      ! energy in K and h2_lev_g the weight g_I (2J+1) with the nuclear spin
      ! factor in it.  Using the observed ladder rather than a rigid rotor
      ! plus a harmonic oscillator is not a refinement for its own sake.
      ! Against a rigid rotor with theta_rot = 85.35 K and a harmonic
      ! oscillator with theta_vib = 5987 K, the ladder gives c_v(H2)/k larger
      ! by 1.4 percent at 200-900 K, 2.3 percent at 1300 K and 7.3 percent at
      ! 3000 K -- and the levels are already here, so there is no second copy
      ! of the H2 spectroscopic constants to drift out of step with the
      ! first.
      !
      ! ONE POTENTIAL FOR THE ENERGY AND FOR THE CHEMISTRY.  The Boltzmann
      ! sum over this ladder, h2_partition_function, is also the internal
      ! partition function the H + H <-> H2 equilibrium constant of mol_rates
      ! is built from, and its zero, the v = 0, J = 0 level, is the level the
      ! dissociation energy D0 there is measured from.  u_rv and ln K_eq are
      ! therefore the first moment and the free energy of the same Q, and are
      ! thermodynamically consistent by construction: no second level model
      ! exists to disagree with this one.
      !
      ! Because the sum runs over the single equilibrium partition function
      ! with the spin weights inside it, this is the ortho/para EQUILIBRIUM
      ! mixture.  A frozen 3:1 mixture gives the same c_rv to four decimals at
      ! and above 300 K and differs only below ~250 K, which is colder than
      ! the molecular layer reaches, so the distinction does not arise here.
      !
      ! RANGE OF VALIDITY.  The ladder is the electronic ground state only, so
      ! this is the heat capacity of BOUND H2.  Dissociation is not a heat
      ! capacity and is not in here: it is a chemical source term and is
      ! accounted separately (mol_rates / lyman_werner).  Above ~5000 K
      ! essentially no H2 survives, and the table's own upper end is 50000 K.
      !
      ! BYTE IDENTITY.  A run with no molecules must reproduce the constant
      ! gamma_ad arithmetic to the bit, so the branch is STRUCTURAL, never a
      ! limit of the formula: with n(H2) = 0 the mixture expression gives
      ! 1 + 1/1.5, which is not the same double as 5/3.  Two gates:
      !   * caloric_mixture_active is .false. for a run with molecular
      !     chemistry off, and every routine then evaluates the legacy
      !     expression verbatim;
      !   * cell by cell, molecular_cell(j) is .false. where n(H2) is exactly
      !     zero, and those cells take the same legacy expression.

      use global_parameters
      use species_table, only: isp_H2
      use molecular_infrared_data, only: n_h2_lev, h2_lev_T, h2_lev_g

      implicit none
      private

      public :: caloric_state_from_composition
      public :: caloric_mixture_active, molecular_cell
      public :: pressure_from_energy_density, energy_density_from_pressure
      public :: adiabatic_index_from_state, adiabatic_index_at_T
      public :: heat_capacity_per_particle, internal_energy_per_particle
      public :: temperature_from_energy_per_particle
      ! The same three maps addressed by the H2 particle fraction instead of
      ! by a cell index, for the post-process energy equation, which walks a
      ! saved profile rather than the live grid.
      public :: heat_capacity_of_mixture, internal_energy_of_mixture
      public :: temperature_of_mixture
      public :: h2_particle_fraction
      public :: h2_rovibrational_energy_and_heat_capacity
      ! The direct level sum and the partition function it is a moment of.
      ! Public because the chemistry (mol_rates) builds its equilibrium
      ! constants from the same Q, and because the thermodynamic identity
      ! u_rv = T^2 d ln Q/dT is asserted against the exact sum.
      public :: h2_rovibrational_sum, h2_partition_function
      public :: caloric_eos_state_line
      ! The composition ratios and the tabulated H2 ladder, for an evaluation
      ! of the same maps in another real kind (hydrodynamic_rows_body.inc).
      ! Read only: they carry no arithmetic of their own.
      public :: caloric_cell_mixture
      public :: h2_rovibrational_table_grid, h2_rovibrational_table_node
      public :: caloric_mixture_record
      public :: hold_caloric_mixture, restore_caloric_mixture

      ! THE COMPOSITION THE CALORIC MAPS READ, AS ONE RECORD, so that a
      ! caller that needs the maps of another composition for a moment (the
      ! pass snapshot of the stationary outer iteration, which evaluates the
      ! state the next pass starts from) can put the module back exactly as
      ! it found it and leave the trajectory of the run untouched.
      type :: caloric_mixture_record
         logical :: held   = .false.
         logical :: active = .false.
         real*8,  allocatable :: nk_per_mass(:), x_h2(:)
         logical, allocatable :: molecular_cell(:)
      end type caloric_mixture_record

      ! Is any cell of this run molecular?  Set by
      ! caloric_state_from_composition; .false. keeps every map on the
      ! constant gamma_ad arithmetic of an atomic gas.
      logical :: caloric_mixture_active = .false.

      ! The composition of each cell that the maps need.  Both are RATIOS of number
      ! densities, so both are functions of the mass fractions alone and are
      ! unchanged by a change of density: a reconstructed interface state may
      ! use its owning cell's values without any lag in rho.  They change only
      ! when the composition does, which is exactly when
      ! caloric_state_from_composition is called (from get_species_densities,
      ! the single point that turns (rho, f_sp) into number densities).
      !
      ! COROLLARY, AND IT BINDS EVERY CALLER: a path that installs a
      ! different f_sp must refresh these arrays BEFORE it maps an energy
      ! density to a pressure, or the pressure it gets belongs to the
      ! previous composition. The energy-to-pressure map is an inverse of
      ! the caloric equation of state, not an arithmetic identity, so the
      ! error does not cancel anywhere downstream. This is why
      ! rebuild_state_from_checkpoint calls get_species_densities before
      ! U_to_W on a restored state.
      !
      ! ONE-SWEEP LAG AT THE GHOSTS, stated because it is real.  The steady
      ! residual fills the ghosts (Apply_BC) BEFORE its ionization sweep runs,
      ! so the ghost conservative state it writes carries the composition of
      ! the previous sweep while the fluxes assembled afterwards carry this
      ! one.  Under a constant gamma the two could not disagree.  This is the
      ! same lag the base-ghost particle count n_part_cell1 already has, and
      ! for the same reason: moving the boundary condition after the sweep
      ! would add a second u -> W -> u round trip, which is not the identity
      ! bitwise and would move every atomic golden.  The marching loop has no
      ! such lag -- there Apply_BC follows the composition refresh.
      real*8,  allocatable :: nk_per_mass(:)   ! (n_tot + n_e)/rho
      real*8,  allocatable :: x_h2(:)          ! n(H2)/(n_tot + n_e)
      logical, allocatable :: molecular_cell(:)

      ! Tabulated mean rovibrational energy of H2, u_rv(T) [K], and its
      ! derivative T du/dT on a logarithmic temperature grid.  The table is
      ! interpolated by a CUBIC HERMITE in ln T built on exact nodal values of
      ! both u_rv and du_rv/dlnT, so the interpolant is C^1 and its own
      ! derivative -- which is what the Newton inverse below differentiates --
      ! is exact at the nodes and smooth between them.  A JFNK residual
      ! differentiates this map, so a piecewise-constant heat capacity (what
      ! linear interpolation of u_rv would give) is not acceptable.
      !
      ! Why a table at all: the direct 302-level sum costs one exponential per
      ! level per evaluation, and the energy -> temperature inverse evaluates
      ! it several times per cell per conversion.  The table reduces that to
      ! one bracket and one cubic.
      integer, parameter :: n_utab   = 4096
      real*8,  parameter :: T_utab_lo = 1.0d0        ! [K]
      real*8,  parameter :: T_utab_hi = 5.0d4        ! [K]
      real*8  :: lnT_utab_lo, dlnT_utab
      real*8  :: u_utab(n_utab)       ! u_rv at the node [K]
      real*8  :: dudlnT_utab(n_utab)  ! T c_rv at the node [K]
      logical :: utab_built = .false.

      contains

      ! ------------------------------------------------------!

      subroutine caloric_state_from_composition(rho, f_sp, ne, n_tot)
      ! Refresh the composition of each cell that the caloric maps read.  Called from
      ! get_species_densities with the same (rho, f_sp, n_e, n_tot) that set
      ! the thermal EOS, so the two halves of the equation of state can never
      ! describe different gas.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho, ne, n_tot
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8  :: nk, nh2
      integer :: j

      if (.not. allocated(nk_per_mass)) then
         allocate(nk_per_mass(1-Ng:N+Ng), x_h2(1-Ng:N+Ng),                &
                  molecular_cell(1-Ng:N+Ng))
         nk_per_mass    = 0.0d0
         x_h2           = 0.0d0
         molecular_cell = .false.
      endif

      if (.not. thereis_mol) then
         caloric_mixture_active = .false.
         molecular_cell         = .false.
         return
      endif

      if (.not. utab_built) call build_h2_rovibrational_table()

      caloric_mixture_active = .false.
      do j = 1-Ng, N+Ng
         nk  = n_tot(j) + ne(j)
         nh2 = rho(j)*f_sp(j,isp_H2)
         ! A cell qualifies only if it has both a positive particle count to
         ! divide by and H2 in it.  The test is the negation of "greater than
         ! zero", so a NaN density leaves the cell on the legacy branch and
         ! the NaN propagates to the run's own detector instead of being
         ! turned into a temperature here.
         if (.not. (nk .gt. 0.0d0) .or. .not. (nh2 .gt. 0.0d0)) then
            molecular_cell(j) = .false.
            nk_per_mass(j)    = 0.0d0
            x_h2(j)           = 0.0d0
         else
            molecular_cell(j)      = .true.
            nk_per_mass(j)         = nk/rho(j)
            x_h2(j)                = nh2/nk
            caloric_mixture_active = .true.
         endif
      enddo

      end subroutine caloric_state_from_composition

      ! ------------------------------------------------------!

      subroutine hold_caloric_mixture(rec)
      ! Copy the composition the caloric maps read into rec.
      type(caloric_mixture_record), intent(out) :: rec
      rec%active = caloric_mixture_active
      rec%held   = allocated(nk_per_mass)
      if (rec%held) then
         rec%nk_per_mass    = nk_per_mass
         rec%x_h2           = x_h2
         rec%molecular_cell = molecular_cell
      endif
      end subroutine hold_caloric_mixture

      ! ------------------------------------------------------!

      subroutine restore_caloric_mixture(rec)
      ! Put back the composition hold_caloric_mixture copied into rec. The
      ! arrays are allocated by the first refresh and never deallocated, so
      ! a record taken before that first refresh restores the flag alone.
      type(caloric_mixture_record), intent(in) :: rec
      caloric_mixture_active = rec%active
      if (rec%held .and. allocated(nk_per_mass)) then
         nk_per_mass    = rec%nk_per_mass
         x_h2           = rec%x_h2
         molecular_cell = rec%molecular_cell
      endif
      end subroutine restore_caloric_mixture

      ! ------------------------------------------------------!

      function caloric_eos_state_line() result(line)
      ! One line for the startup configuration report.
      character(len=96) :: line
      if (caloric_eos_monatomic) then
         line = 'monatomic gamma = 5/3 by request ("Caloric EOS:'//     &
                ' monatomic", a comparison option; H2 stores 3/2 kT)'
      else if (caloric_mixture_active) then
         line = 'composition-dependent: H2 rovibrational ladder'//        &
                ' (Roueff+2019), atoms/ions/electrons 3/2 k'
      else if (thereis_mol) then
         line = 'monatomic gamma = 5/3 (molecular chemistry on,'//        &
                ' no H2 in any cell yet)'
      else
         line = 'monatomic gamma = 5/3'
      endif
      end function caloric_eos_state_line

      subroutine caloric_cell_mixture(j, nk_per_rho, x2, is_molecular)
      ! The two composition ratios the caloric maps of cell j are evaluated
      ! at: (n_tot + n_e)/rho and n(H2)/(n_tot + n_e).  Both are ratios of
      ! number densities and so are unchanged by a change of density, which
      ! is why a reconstructed face state may use its owning cell's values.
      ! is_molecular is .false. wherever the maps take the constant gamma_ad
      ! branch, and the two ratios are then not defined.
      integer, intent(in)  :: j
      real*8,  intent(out) :: nk_per_rho, x2
      logical, intent(out) :: is_molecular

      nk_per_rho  = 0.0d0
      x2          = 0.0d0
      is_molecular = .false.
      if (.not. caloric_mixture_active) return
      if (.not. allocated(molecular_cell)) return
      if (.not. molecular_cell(j)) return
      nk_per_rho  = nk_per_mass(j)
      x2          = x_h2(j)
      is_molecular = .true.

      end subroutine caloric_cell_mixture

      ! ------------------------------------------------------!

      subroutine h2_rovibrational_table_grid(n_nodes, lnT_lo, dlnT,      &
                                             T_lo, T_hi)
      ! The logarithmic temperature grid the Hermite table below sits on.
      ! Builds the table if it is not there yet, so a caller that reaches the
      ! nodes without going through caloric_state_from_composition still gets
      ! a table.
      integer, intent(out) :: n_nodes
      real*8,  intent(out) :: lnT_lo, dlnT, T_lo, T_hi

      if (.not. utab_built) call build_h2_rovibrational_table()
      n_nodes = n_utab
      lnT_lo  = lnT_utab_lo
      dlnT    = dlnT_utab
      T_lo    = T_utab_lo
      T_hi    = T_utab_hi

      end subroutine h2_rovibrational_table_grid

      ! ------------------------------------------------------!

      subroutine h2_rovibrational_table_node(i, u_node, dudlnT_node)
      ! Node i of the table: u_rv [K] and T c_rv [K] at that temperature.
      integer, intent(in)  :: i
      real*8,  intent(out) :: u_node, dudlnT_node

      if (.not. utab_built) call build_h2_rovibrational_table()
      u_node      = u_utab(i)
      dudlnT_node = dudlnT_utab(i)

      end subroutine h2_rovibrational_table_node

      ! ------------------------------------------------------!

      ! ------------------------------------------------------!
      ! The H2 rovibrational ladder
      ! ------------------------------------------------------!

      subroutine build_h2_rovibrational_table()
      ! Tabulate u_rv and T c_rv on the logarithmic grid, from the exact sum.
      integer :: i
      real*8  :: TK, u, c

      lnT_utab_lo = log(T_utab_lo)
      dlnT_utab   = (log(T_utab_hi) - lnT_utab_lo)/dble(n_utab - 1)
      do i = 1, n_utab
         TK = exp(lnT_utab_lo + dlnT_utab*dble(i-1))
         call h2_rovibrational_sum(TK, u, c)
         u_utab(i)      = u
         dudlnT_utab(i) = TK*c
      enddo
      utab_built = .true.

      end subroutine build_h2_rovibrational_table

      ! ------------------------------------------------------!

      subroutine h2_rovibrational_sum(TK, u_rv, c_rv, Q)
      ! Internal partition function Q, mean rovibrational energy
      ! u_rv = <E>/k [K] and heat capacity c_rv = du_rv/dT [dimensionless] of
      ! one H2 molecule at TK [K], summed directly over the bound ladder.
      ! c_rv is the energy variance, (<E^2> - <E>^2)/(k T)^2, which is the
      ! exact derivative of the mean.
      !
      ! ONE POTENTIAL.  Q, u_rv and c_rv are the zeroth, first and second
      ! moments of ONE Boltzmann sum over ONE level set, so they satisfy
      ! u_rv = T^2 d ln Q/dT and c_rv = du_rv/dT identically rather than by
      ! agreement between two models.  Q is what the H + H <-> H2
      ! equilibrium constant of mol_rates is built from and u_rv is what the
      ! energy equation carries, so the chemistry and the equation of state
      ! hold the same H2 molecule.
      !
      ! The zero of the sum is the v = 0, J = 0 level (the zero of
      ! h2_lev_T), which is also the level the dissociation energy D0 of
      ! mol_rates is measured from and the level the H2 entry of the
      ! formation reservoir uses; hence Q -> 1 as T -> 0.
      !
      ! The exponentials are taken relative to the ground level, so nothing
      ! overflows; at the cold end every excited term underflows to zero and
      ! the sum correctly returns Q = 1, u_rv = 0.
      real*8, intent(in)  :: TK
      real*8, intent(out) :: u_rv, c_rv
      real*8, intent(out), optional :: Q
      real*8  :: Z, S1, S2, w, E, Tsafe
      integer :: i

      Tsafe = max(TK, 1.0d-30)
      Z  = 0.0d0;  S1 = 0.0d0;  S2 = 0.0d0
      do i = 1, n_h2_lev
         E = h2_lev_T(i)
         w = h2_lev_g(i)*exp(-E/Tsafe)
         Z  = Z  + w
         S1 = S1 + w*E
         S2 = S2 + w*E*E
      enddo
      u_rv = S1/Z
      c_rv = (S2/Z - u_rv*u_rv)/(Tsafe*Tsafe)
      if (present(Q)) Q = Z

      end subroutine h2_rovibrational_sum

      ! ------------------------------------------------------!

      double precision function h2_partition_function(TK) result(Q)
      ! Internal (rovibrational) partition function of H2 X^1 Sigma_g^+ at
      ! TK [K], measured from v = 0, J = 0, over the observed bound ladder of
      ! Roueff et al. (2019), A&A 630, A58, table 2.
      !
      ! THE NUCLEAR SPIN WEIGHTS ARE INSIDE THE SUM.  h2_lev_g is
      ! g = g_I (2J+1) with g_I = 1 for even J (para) and 3 for odd J
      ! (ortho), so this is the single equilibrium partition function of the
      ! ortho/para mixture, not a para-only sum and not a frozen 3:1 mixture.
      ! Any equilibrium constant built from it must count the nuclear spin of
      ! the free atoms on the same convention (2 spin states each), which is
      ! what mol_rates does.
      !
      ! This is the ONLY H2 internal partition function in the code: the
      ! equilibrium constants of mol_rates and the caloric equation of state
      ! read it, so one potential generates both.
      real*8, intent(in) :: TK
      real*8 :: u_rv, c_rv

      call h2_rovibrational_sum(TK, u_rv, c_rv, Q)

      end function h2_partition_function

      ! ------------------------------------------------------!

      subroutine h2_rovibrational_energy_and_heat_capacity(TK, u_rv, c_rv)
      ! The same two quantities from the cubic Hermite table.  Outside the
      ! tabulated range the interpolant is continued LINEARLY in T with the
      ! end-point heat capacity, which keeps u_rv monotone in T -- the
      ! property the energy -> temperature inverse below relies on.
      real*8, intent(in)  :: TK
      real*8, intent(out) :: u_rv, c_rv
      real*8  :: x, s, s2, s3, h00, h10, h01, h11, dh00, dh10, dh01, dh11
      integer :: i

      ! Lazy build.  Every run with H2 in it reaches this through
      ! caloric_state_from_composition first, which is serial, so the
      ! table is already there by the time any parallel region can ask;
      ! the test is kept so that a caller that arrives another way (a
      ! unit test, the post-process energy equation) still gets a table.
      if (.not. utab_built) call build_h2_rovibrational_table()

      ! Monatomic comparison option: H2 stores no rotational or vibrational
      ! energy, so u_rv = c_rv = 0 and gamma_eff = 5/3 for any composition.
      if (caloric_eos_monatomic) then
         u_rv = 0.0d0; c_rv = 0.0d0
         return
      endif

      if (TK .le. T_utab_lo) then
         c_rv = dudlnT_utab(1)/T_utab_lo
         u_rv = u_utab(1) + c_rv*(TK - T_utab_lo)
         if (u_rv .lt. 0.0d0) u_rv = 0.0d0
         return
      endif
      if (TK .ge. T_utab_hi) then
         c_rv = dudlnT_utab(n_utab)/T_utab_hi
         u_rv = u_utab(n_utab) + c_rv*(TK - T_utab_hi)
         return
      endif

      x = (log(TK) - lnT_utab_lo)/dlnT_utab
      i = 1 + int(x)
      if (i .lt. 1)          i = 1
      if (i .gt. n_utab - 1) i = n_utab - 1
      s  = x - dble(i-1)
      s2 = s*s;  s3 = s2*s

      ! Cubic Hermite basis on [0,1] and its derivative with respect to s.
      h00 =  2.0d0*s3 - 3.0d0*s2 + 1.0d0
      h10 =        s3 - 2.0d0*s2 + s
      h01 = -2.0d0*s3 + 3.0d0*s2
      h11 =        s3 -       s2
      dh00 =  6.0d0*s2 - 6.0d0*s
      dh10 =  3.0d0*s2 - 4.0d0*s + 1.0d0
      dh01 = -6.0d0*s2 + 6.0d0*s
      dh11 =  3.0d0*s2 - 2.0d0*s

      u_rv = h00*u_utab(i)      + h10*dlnT_utab*dudlnT_utab(i)            &
           + h01*u_utab(i+1)    + h11*dlnT_utab*dudlnT_utab(i+1)
      ! du/dT = (du/ds)/(dlnT_utab * T)
      c_rv = ( dh00*u_utab(i)   + dh10*dlnT_utab*dudlnT_utab(i)           &
             + dh01*u_utab(i+1) + dh11*dlnT_utab*dudlnT_utab(i+1) )       &
             /(dlnT_utab*TK)

      end subroutine h2_rovibrational_energy_and_heat_capacity

      ! ------------------------------------------------------!
      ! The maps
      ! ------------------------------------------------------!

      real*8 function internal_energy_of_mixture(x2, T)
      ! Internal energy per (n_tot + n_e) of a gas in which a fraction x2 of
      ! the particle-plus-electron count is H2, in units of k T0, at the code
      ! temperature T.  x2 = 0 is the monatomic value written exactly as
      ! T/(gamma_ad - 1) so the legacy arithmetic survives to the bit.
      real*8, intent(in) :: x2, T
      real*8 :: u_rv, c_rv

      if (x2 .le. 0.0d0) then
         internal_energy_of_mixture = T/(gamma_ad - 1.0d0);  return
      endif
      call h2_rovibrational_energy_and_heat_capacity(T*T0, u_rv, c_rv)
      internal_energy_of_mixture = 1.5d0*T + x2*u_rv/T0

      end function internal_energy_of_mixture

      ! ------------------------------------------------------!

      real*8 function heat_capacity_of_mixture(x2, T)
      ! C_V/((n_tot + n_e) k) = 1/(gamma_eff - 1) of the same mixture.
      real*8, intent(in) :: x2, T
      real*8 :: u_rv, c_rv

      if (x2 .le. 0.0d0) then
         heat_capacity_of_mixture = 1.0d0/(gamma_ad - 1.0d0);  return
      endif
      call h2_rovibrational_energy_and_heat_capacity(T*T0, u_rv, c_rv)
      heat_capacity_of_mixture = 1.5d0 + x2*c_rv

      end function heat_capacity_of_mixture

      ! ------------------------------------------------------!

      real*8 function h2_particle_fraction(j)
      ! n(H2)/(n_tot + n_e) of cell j; zero where the cell has no molecules.
      integer, intent(in) :: j
      if (caloric_mixture_active) then
         if (molecular_cell(j)) then
            h2_particle_fraction = x_h2(j)
         else
            h2_particle_fraction = 0.0d0
         endif
      else
         h2_particle_fraction = 0.0d0
      endif
      end function h2_particle_fraction

      ! ------------------------------------------------------!

      real*8 function internal_energy_per_particle(j, T)
      ! The mixture energy above at the composition of cell j.
      integer, intent(in) :: j
      real*8,  intent(in) :: T

      if (.not. caloric_mixture_active) then
         internal_energy_per_particle = T/(gamma_ad - 1.0d0);  return
      endif
      if (.not. molecular_cell(j)) then
         internal_energy_per_particle = T/(gamma_ad - 1.0d0);  return
      endif
      internal_energy_per_particle = internal_energy_of_mixture(x_h2(j), T)

      end function internal_energy_per_particle

      ! ------------------------------------------------------!

      real*8 function heat_capacity_per_particle(j, T)
      ! The mixture heat capacity above at the composition of cell j.
      integer, intent(in) :: j
      real*8,  intent(in) :: T

      if (.not. caloric_mixture_active) then
         heat_capacity_per_particle = 1.0d0/(gamma_ad - 1.0d0);  return
      endif
      if (.not. molecular_cell(j)) then
         heat_capacity_per_particle = 1.0d0/(gamma_ad - 1.0d0);  return
      endif
      heat_capacity_per_particle = heat_capacity_of_mixture(x_h2(j), T)

      end function heat_capacity_per_particle

      ! ------------------------------------------------------!

      real*8 function adiabatic_index_at_T(j, T)
      ! gamma_eff = 1 + (n_tot + n_e) k/C_V of this cell at code temperature T.
      integer, intent(in) :: j
      real*8,  intent(in) :: T

      if (.not. caloric_mixture_active) then
         adiabatic_index_at_T = gamma_ad;  return
      endif
      if (.not. molecular_cell(j)) then
         adiabatic_index_at_T = gamma_ad;  return
      endif
      adiabatic_index_at_T = 1.0d0 + 1.0d0/heat_capacity_per_particle(j, T)

      end function adiabatic_index_at_T

      ! ------------------------------------------------------!

      real*8 function temperature_of_mixture(x2, u_in)
      ! Invert internal_energy_of_mixture for T.  With the vibrational ladder
      ! in c_v this is a genuinely nonlinear scalar equation, and it is solved
      ! HERE, to full precision, rather than lagged at a previous temperature:
      ! a residual that depended on the previous state would blur the fixed
      ! point and the finite-difference Jacobian of the JFNK solver alike.
      !
      ! u_rv is non-negative and strictly increasing in T, so
      !     f(T) = 1.5 T + x2 u_rv(T T0)/T0 - u_in
      ! is strictly increasing and 0 < T <= u_in/1.5 brackets its root.
      ! Newton from the upper end is safeguarded by bisection on that bracket,
      ! so it cannot leave it whatever the curvature does.
      real*8, intent(in) :: x2, u_in
      real*8  :: Tlo, Thi, T, Tnew, f, df, u_rv, c_rv
      integer :: it
      integer, parameter :: max_it = 80

      ! No molecules, or not a physical energy (or a NaN): the monatomic value,
      ! written exactly as the legacy expression, and a NaN passes through to
      ! the caller's own positivity guard.
      if (x2 .le. 0.0d0 .or. .not. (u_in .gt. 0.0d0)) then
         temperature_of_mixture = u_in*(gamma_ad - 1.0d0)
         return
      endif

      Tlo = 0.0d0
      Thi = u_in/1.5d0
      T   = Thi
      do it = 1, max_it
         call h2_rovibrational_energy_and_heat_capacity(T*T0, u_rv, c_rv)
         f  = 1.5d0*T + x2*u_rv/T0 - u_in
         df = 1.5d0 + x2*c_rv
         ! An exact zero IS the root; the bracket update below would file it as
         ! an endpoint and the safeguard would then bisect it away.  The
         ! equality test on a real is deliberate and is not a tolerance in
         ! disguise (gfortran -Wcompare-reals flags it): the question asked is
         ! literally "did this evaluation return the floating-point zero", and
         ! any tolerance around it would stop the iteration early at a T that
         ! is not the root.
         if (f .eq. 0.0d0) exit
         if (f .gt. 0.0d0) then
            Thi = T
         else
            Tlo = T
         endif
         Tnew = T - f/df
         ! A Newton step that leaves the bracket is replaced by a bisection of
         ! it.  The test is on the STEP THE ITERATION ACTUALLY TAKES, not on
         ! the Newton step that may have been discarded: a converged Newton
         ! step lands on a bracket endpoint, and testing the discarded step
         ! would leave the loop running while the safeguard threw the root
         ! away and returned the midpoint of the bracket instead.
         if (Tnew .le. Tlo .or. Tnew .ge. Thi) Tnew = 0.5d0*(Tlo + Thi)
         if (abs(Tnew - T) .le. 1.0d-15*T) then
            T = Tnew
            exit
         endif
         T = Tnew
      enddo
      temperature_of_mixture = T

      end function temperature_of_mixture

      ! ------------------------------------------------------!

      real*8 function temperature_from_energy_per_particle(j, u_in)
      ! The inverse above at the composition of cell j.
      integer, intent(in) :: j
      real*8,  intent(in) :: u_in

      if (.not. caloric_mixture_active) then
         temperature_from_energy_per_particle = u_in*(gamma_ad - 1.0d0)
         return
      endif
      if (.not. molecular_cell(j)) then
         temperature_from_energy_per_particle = u_in*(gamma_ad - 1.0d0)
         return
      endif
      temperature_from_energy_per_particle =                              &
         temperature_of_mixture(x_h2(j), u_in)

      end function temperature_from_energy_per_particle

      ! ------------------------------------------------------!

      real*8 function pressure_from_energy_density(j, rho_c, rho_e)
      ! p from the internal energy density rho e of a state whose composition
      ! is that of cell j.  rho_c is the density of THAT state, which for a
      ! reconstructed face is not the cell average; the composition ratios
      ! this module stores are density independent, so the same cell entry is
      ! the right one at the face.
      integer, intent(in) :: j
      real*8,  intent(in) :: rho_c, rho_e
      real*8 :: nk, T

      if (.not. caloric_mixture_active) then
         pressure_from_energy_density = (gamma_ad - 1.0d0)*rho_e;  return
      endif
      if (.not. molecular_cell(j) .or. .not. (rho_c .gt. 0.0d0)) then
         pressure_from_energy_density = (gamma_ad - 1.0d0)*rho_e;  return
      endif
      nk = rho_c*nk_per_mass(j)
      T  = temperature_from_energy_per_particle(j, rho_e/nk)
      pressure_from_energy_density = nk*T

      end function pressure_from_energy_density

      ! ------------------------------------------------------!

      real*8 function energy_density_from_pressure(j, rho_c, p)
      ! The inverse map, rho e from p.  Direct: T = p/(n_tot + n_e) needs no
      ! iteration, and the energy follows from it.
      integer, intent(in) :: j
      real*8,  intent(in) :: rho_c, p
      real*8 :: nk

      if (.not. caloric_mixture_active) then
         energy_density_from_pressure = p/(gamma_ad - 1.0d0);  return
      endif
      if (.not. molecular_cell(j) .or. .not. (rho_c .gt. 0.0d0)) then
         energy_density_from_pressure = p/(gamma_ad - 1.0d0);  return
      endif
      nk = rho_c*nk_per_mass(j)
      energy_density_from_pressure = nk*internal_energy_per_particle(j, p/nk)

      end function energy_density_from_pressure

      ! ------------------------------------------------------!

      real*8 function adiabatic_index_from_state(j, rho_c, p)
      ! gamma_eff of a state (rho_c, p) at the composition of cell j.  Used by
      ! the sound speeds: c = sqrt(gamma_eff p/rho).
      integer, intent(in) :: j
      real*8,  intent(in) :: rho_c, p

      if (.not. caloric_mixture_active) then
         adiabatic_index_from_state = gamma_ad;  return
      endif
      if (.not. molecular_cell(j) .or. .not. (rho_c .gt. 0.0d0)) then
         adiabatic_index_from_state = gamma_ad;  return
      endif
      adiabatic_index_from_state =                                        &
         adiabatic_index_at_T(j, p/(rho_c*nk_per_mass(j)))

      end function adiabatic_index_from_state

      ! End of module
      end module caloric_eos

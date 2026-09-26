      module mol_rates
      ! molecular (H2 / H2+ / H3+ / HeH+) reaction-rate coefficients, on
      ! the reaction list of Koskinen et al. (2022), ApJ 929:52, Table 1
      ! (page 19; verified against the PDF in references/).  All rates in
      ! cgs (cm^3 s^-1; the three-body rates R13/R15 are returned as
      ! two-body equivalents, the cm^6 s^-1 coefficient times the third-body
      ! density, as in the table's "n" factor).
      ! T = heavy-particle temperature [K]; Te = electron temperature [K]
      ! (EXHALE uses a common T).
      !
      ! R1-R11, R13, R14, R16, R18, R19 and R23 are the Table-1 expressions
      ! as printed.  The entries that are NOT, each argued at its own
      ! function:
      !   R12  the thermal dissociation of H2, built by detailed balance
      !        from R15 instead of from the Table-1 fit, because the two
      !        Table-1 entries are the same channel in opposite directions
      !        and their source measurements do not overlap in temperature;
      !   R15  the Cohen & Westberg (1983) recommended coefficient with the
      !        third bodies resolved, in place of the Ham et al. (1970)
      !        room-temperature value Table 1 prints;
      !   R17  the measured two-body total of Boehringer & Arnold (1986)
      !        less the radiative branch R23, plus the Table-1 Arrhenius
      !        term;
      !   R20  not carried: its cited source bounds the channel far below
      !        the Table-1 value (the retired-R20 block below);
      !   R21, R22  not in this module: the He <-> H charge-exchange pair
      !        is charge_exchange::he_h_cx_rates, applied in every system
      !        that carries He.
      ! One reaction is not in Table 1 at all: H2+ + He -> HeH+ + H (Black
      ! 1978; rk_H2p_He_HeHp), the HeH+ source Table 1 omits.
      !
      ! Photo-rates P1-P5 (H, He, H2 photoionization; H2 dissociative and
      ! double photoionization) are "SC" in the paper -- computed from cross
      ! sections x stellar flux x column densities.  In EXHALE they are the
      ! photoionization rates util_ion_eq builds from the cross sections of
      ! cross_sec.f90 (P_H2, P_H2_di, P_H2_dd; the double-ionization part
      ! is zero unless a model for it is selected); they are NOT part of
      ! this module.
      !
      ! THIRD BODY M (R12, R13, R15).  Koskinen et al. (2022) write the
      ! three-body rates as a two-body coefficient times "n" and do not say
      ! which particles that n counts; their models are H2/H-dominated
      ! Jovian and Neptunian envelopes, where M is H2 and H.
      !
      ! THE COLLIDERS OF R12/R15 ARE RESOLVED, and the two reactions take
      ! the H2-EQUIVALENT third-body density of
      ! h2_association_collider_density below, not the total heavy-particle
      ! density.  Cohen & Westberg (1983), J. Phys. Chem. Ref. Data 12, 531,
      ! p. 559, recommend a separate coefficient for each third body on
      ! their own data sheets (READ): k1(H2) = 2.8e-31 T^-0.6 over
      ! 50-5000 K, k1(H) = 8.8e-33 temperature independent over 50-5000 K,
      ! k1(Ar) = 1.9e-30 T^-1 over 77-5000 K, all cm^6 molecule^-2 s^-1.  At
      ! 1000 K atomic hydrogen is 2.0 times as efficient a third body as H2
      ! and a monatomic inert gas 0.43 times as efficient, so one
      ! coefficient applied to the total density overstates the monatomic
      ! colliders by about a factor 2 and understates the atomic-hydrogen
      ! one by about the same factor.
      !
      ! R13 is NOT collider-resolved: its source (Miller et al. 1968) is a
      ! single coefficient and no evaluation in hand splits it by third
      ! body, so it keeps the total density.
      !
      ! Helium has no coefficient of its own anywhere in hand.  The only
      ! published helium third-body calculation, Paolini, Ohlinger & Forrey
      ! (2011), Phys. Rev. A 83, 042713, resolves the product states but
      ! tabulates no rate: its total is shown in figures over 0-350 K only
      ! (READ, their Figs. 6 and 7), which neither reaches this layer's
      ! 200-1600 K nor gives a temperature dependence.  Helium is therefore
      ! given Cohen & Westberg's argon coefficient, the nearest tabulated
      ! monatomic inert third body of the same evaluation.  That is an
      ! ARGON-BASED ESTIMATE FOR HELIUM and not a bound on the helium
      ! efficiency; its range (77-5000 K), its stated uncertainty (+/-0.3 in
      ! log) and the reason the two atoms are not ordered by mass and
      ! polarizability are at the coefficient.  The substitution is stated
      ! there and is not a fit.
      !
      ! HE-DOMINATED LIMIT.  Two further limits of this network are noted
      ! rather than patched here, because both would replace a
      ! Koskinen Table-1 entry with a rate from another compilation:
      !   * HeH+ formation.  REPAIRED 2026-09-15, and the
      !     paragraph is kept because the repair is the point.  Table 1
      !     forms HeH+ only through R20, He+ + H2 -> HeH+ + H at 4.2e-13,
      !     citing Schauer et al. (1989) -- a measurement that neither made
      !     that channel nor supports that number (see the retired-R20 block
      !     below).  R20 is therefore GONE from this module, and the HeH+
      !     source the network now carries is the one Table 1 omitted,
      !     H2+ + He -> HeH+ + H, from the primary source (Black 1978;
      !     rk_H2p_He_HeHp).  It needs no He+ and is therefore the route
      !     that grows with the helium fraction.  The ~10% associative
      !     branch of He(2^3S) + H / + H2 (Garcia Munoz 2025, Appendix) is
      !     the other HeH+ source and the molecular system carries it
      !     (row 7 of mol_heh_rows).
      !   * HeH+ destruction. Against the Garcia Munoz (2025) Table A.6
      !     values, R16 is 3.4-8.6x smaller, R18 1.2x larger and R19
      !     1.4-2.6x smaller.  Both compilations are published and the
      !     differences are of order unity to ten; no entry is overruled.
      !   * He+ + H2 BELOW ~650 K.  R17 is the dissociative charge transfer,
      !     and the Koskinen Table-1 form for it, 1e-9 exp(-5700/T) (Moses
      !     & Bass 2000), is an extrapolation there: three thermal
      !     measurements put that channel at 1.1e-13 to 1.5e-13 over
      !     15-300 K where the fit gives 5.6e-18 at 300 K.  The numbers and
      !     the sources are at R17.  Not overruled here, because replacing
      !     it needs the two measurements' own papers, which were not
      !     obtained.
      !
      ! NOTE (Koskinen 2022 baseline): neutral H2 photodissociation through
      ! the Lyman-Werner bands is NOT part of their Table 1; their
      ! sensitivity test with a Backx et al. (1976) cross section and
      ! dissociation probability 0.125 changed Mdot by <= 1.4x.  EXHALE adds
      ! it separately and opt-in, in
      ! src/modules/lower_atmosphere/lyman_werner.f90 (the dissociation
      ! cross section of the line-by-line self-shielding table
      ! h2_self_shielding_table.f90; key "Stellar LW flux"), so it is
      ! deliberately absent here.

      ! Physical constants and the hydrogen-atom mass, from the one place
      ! that owns them.  Nothing else of global_parameters is in scope here,
      ! so the local names of this module are checked against these six and
      ! against the module's own declarations only.
      use global_parameters, only: kb_erg, kb_eV, hp_erg, c_light, mu, pi
      ! The internal partition function of H2 and the level ladder behind it
      ! live with the caloric equation of state, which is the other consumer
      ! of the same Boltzmann sum.  See the H2 THERMOCHEMISTRY block below.
      use caloric_eos, only: h2_partition_function
      ! The vibrational heat fraction is built from the same H2 level ladder
      ! as the equilibrium constant, and its all-level radiative rate is
      ! reduced from the line list by an initializer this module's own
      ! serial prologue calls.  Nothing else of that module is used here.
      use h2_vibrational_relaxation, only: h2_vibrational_relaxation_init

      implicit none
      private
      public :: rk_R1_Hp_rec,   rk_R2_Hep_rec,  rk_R3_H_cion,            &
                rk_R4_He_cion,  rk_R5_H2p_dr,   rk_R6_H3p_dr_H2,         &
                rk_R7_H3p_dr_3H, rk_R8_H2p_H2,  rk_R9_H2p_H,             &
                rk_R10_Hp_H2v4, rk_R11_H3p_H,   rk_R12_H2_thdis,         &
                rk_R13_Hp_H2_M, rk_R14_H2_edis, rk_R15_3body_H2,         &
                rk_R16_HeHp_dr, rk_R17_Hep_H2_diss, rk_R18_HeHp_H2,      &
                rk_R19_HeHp_H,  rk_H2p_He_HeHp,    rk_R23_H2_Hep_cx
      ! The photon R23 emits, read by the reaction-heat ledger so that a
      ! radiative channel does not deposit its own photon as heat, and the
      ! two-body/three-body domain guard of R17.
      public :: rct_photon_energy_eV, hep_h2_three_body_domain
      ! H2 thermochemistry: the three-body recombination coefficient that R15
      ! and R12 share, the H + H <-> H2 equilibrium constant that turns one
      ! into the other, and the table builder.  The internal partition
      ! function the equilibrium constant is built from is not owned here:
      ! it is h2_partition_function of caloric_eos, the one Boltzmann sum
      ! over the one H2 level ladder.  See the H2 THERMOCHEMISTRY block
      ! below.
      public :: k3b_H_H_to_H2, keq_H_H_to_H2, h2_thermochemistry_ready,                             &
                h2_thermochemistry_init, h2_dissociation_energy_eV
      ! The third-body efficiencies of H + H + M -> H2 + M, and the
      ! H2-equivalent collider density they build (see THIRD BODY M above).
      public :: k3b_H_H_to_H2_atomic_H, k3b_H_H_to_H2_monatomic,         &
                h2_association_collider_density

      ! ------------------------------------------------------------------ !
      ! H2 THERMOCHEMISTRY
      !
      ! ONE LEVEL SET.  The internal partition function of H2 that the
      ! equilibrium constant below is built from is h2_partition_function of
      ! caloric_eos: the Boltzmann sum over the observed bound rovibrational
      ! ladder of H2 X^1 Sigma_g^+, 302 levels to 36118 cm^-1, of Roueff et
      ! al. (2019), A&A 630, A58, table 2.  It is the same sum whose first
      ! moment is the rovibrational internal energy the caloric equation of
      ! state carries, and its zero is the v = 0, J = 0 level, which is the
      ! level D0 below is measured from and the level the H2 entry of the
      ! formation reservoir uses.  One potential therefore generates the
      ! chemistry and the energy equation, and the reverse rate R12, the
      ! recombination R15 and the reaction enthalpy cannot disagree about how
      ! much energy a bound H2 molecule holds.
      !
      ! THE NUCLEAR SPIN CONVENTION MATCHES ON BOTH SIDES.  The ladder's
      ! weights are g_I (2J+1) with g_I = 3 for odd J and 1 for even J, so
      ! the molecular sum counts nuclear spin; the free H atom is given
      ! 2 (electronic doublet) x 2 (nuclear spin I = 1/2) = 4 in K_eq below,
      ! which is the matching count, and the spin factors cancel out of the
      ! equilibrium.
      !
      ! THE DISSOCIATION ENERGY FROM v = 0, J = 0, in cm^-1.  It is the one
      ! bond energy of this code.
      !
      ! D0 IS CHECKED AGAINST THE THERMOCHEMISTRY OF THE SAME REACTION.
      ! 36118.11 cm^-1 = 432.05 kJ/mol; adding the 298.15 K enthalpy
      ! functions of two H atoms (2 x 5/2 RT) and removing that of H2
      ! (7/2 RT) gives dH(298.15 K) = 435.8 kJ/mol, against the 436 kJ/mol
      ! tabulated on the H2 (+M) <-> H + H (+M) data sheet of Baulch et al.
      ! (1992) J. Phys. Chem. Ref. Data 21, 411, p. 550.  It also agrees with
      ! the top of the Roueff ladder, whose highest bound level stands at
      ! 51965.8 K = 36118.04 cm^-1: the ladder and D0 share their zero AND
      ! their limit, to 0.07 cm^-1.
      real*8, parameter :: D0_H2_cm    = 36118.11d0   ! D0(v=0,J=0)

      ! Tabulation of ln K_eq on a log-spaced temperature grid.  The level
      ! sum costs 302 exponentials, far too many to evaluate at every cell
      ! inside the Newton solve of the molecular system (~1e9 over a run), so
      ! it is evaluated once here and interpolated linearly in
      ! (ln T, ln K_eq).
      !
      ! GRID.  100-20000 K, 2001 points, i.e. 2.65e-3 in ln T.  ln K_eq is
      ! dominated by D0 hc / kT, whose curvature in ln T is D0 hc / kT and
      ! largest at the cold end; the linear-interpolation error is therefore
      ! bounded by (1/8) (D0 hc/kT) (dlnT)^2 = 4.6e-4 in ln K_eq at 100 K,
      ! i.e. 0.046 % in K_eq, and falls as 1/T from there.  Measured against
      ! the direct sum at the midpoint of every interval, the largest error
      ! is 4.55e-4 (0.046 %) in the first interval and 1.42e-6 (1.4e-4 %) in
      ! the last.
      !
      ! ENDS.  T is clamped into [100, 20000] K before the lookup.  Below
      ! 100 K the clamp overestimates K_eq's fall and so overestimates the
      ! dissociation rate, but K_eq(100 K) = 1.2e202 cm^3 already makes that
      ! rate 1e-234 cm^3 s^-1 -- and 100 K is also as cold as K_eq may be
      ! evaluated at all before exp(D0 hc/kT) overflows double precision.
      ! Above 20000 K the clamp underestimates the dissociation rate, at a
      ! temperature where it is already ~1e-10 cm^3 s^-1 times the total
      ! density and H2 cannot survive a single time step either way.
      !
      ! The table is written once, by h2_thermochemistry_init, and only read
      ! afterwards, so it is NOT threadprivate: every thread of the parallel
      ! cell sweep shares one copy.  EXHALE_main calls the initializer before
      ! any parallel region opens, exactly as it does for the infrared band
      ! tables; the guard inside keq_H_H_to_H2 is a safety net for any other
      ! entry point and never fires in the cell sweep.
      ! ---- the He+ + H2 two-body-only domain guard (see R17) ----------
      ! The ratio above which the omitted three-body channel stops being a
      ! rounding error: k3 n / k2 = 0.1 is a ten percent rate error, which
      ! is the size of the two-body measurement's own scatter.
      real*8,  parameter :: three_body_warn_ratio = 0.1d0
      integer, save :: n_hep_h2_evaluations = 0
      integer, save :: n_hep_h2_three_body  = 0
      real*8,  save :: hep_h2_worst_ratio   = 0.0d0
      logical, save :: hep_h2_warned        = .false.
      integer, parameter :: n_keq    = 2001
      real*8,  parameter :: T_keq_lo = 1.0d2
      real*8,  parameter :: T_keq_hi = 2.0d4
      real*8,  save      :: keq_lnT(n_keq), keq_lnK(n_keq)
      logical, save      :: keq_table_ready = .false.
      ! ------------------------------------------------------------------ !

      contains

      logical function h2_thermochemistry_ready()
      ! Whether h2_thermochemistry_init has been called; the ionization
      ! sweep's serial prologue calls the initializer when it has not, so
      ! that no parallel region ever builds the table.
      h2_thermochemistry_ready = keq_table_ready
      end function h2_thermochemistry_ready


      ! D0(H2) in eV.  The bond energy has ONE definition in this code, the
      ! spectroscopic D0_H2_cm above; anything that needs it as an energy
      ! reads it through here rather than writing 4.478 of its own.
      double precision function h2_dissociation_energy_eV() result(e)
      e = D0_H2_cm*1.239841984d-4
      end function

      ! R1-R4 are the atomic H/He reactions of Table 1.  They are also the
      ! set the input key "Atomic rate set: Koskinen2022" puts in front of
      ! EXHALE's own atomic rates: Cool_coeff calls these four functions in
      ! that mode, so the published coefficients have one definition, here,
      ! whether the run is molecular or atomic.  They are elemental for that
      ! caller, whose accessors are elemental themselves.

      ! R1: H+ + e -> H + hv                    (Storey & Hummer 1995)
      elemental double precision function rk_R1_Hp_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.0d-12*(300.0d0/Te)**0.64d0
      end function

      ! R2: He+ + e -> He + hv                  (Storey & Hummer 1995)
      elemental double precision function rk_R2_Hep_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.6d-12*(300.0d0/Te)**0.64d0
      end function

      ! R3: H + e -> H+ + 2e                    (Voronov 1997)
      ! U = 13.6 eV / E_e;  E_e = kB Te in eV.  The eV Boltzmann constant is
      ! the CODATA value owned by global_parameters, not a rounded local
      ! copy: it is not part of the published fit.  The 13.6 eV IS: it is
      ! Voronov's Table I entry, fitted together with the four coefficients
      ! below, so it stays as published rather than reading the measured
      ! threshold e_th_HI.  Same for the 24.6 eV of R4.
      elemental double precision function rk_R3_H_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 13.6d0/(kb_eV*Te)
      k = 2.91d-8*U**0.39d0*exp(-U)/(0.232d0 + U)
      end function

      ! R4: He + e -> He+ + 2e                  (Voronov 1997)
      elemental double precision function rk_R4_He_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 24.6d0/(kb_eV*Te)
      k = 1.75d-8*U**0.35d0*exp(-U)/(0.180d0 + U)
      end function

      ! R5: H2+ + e -> H + H                    (Auerbach et al. 1977)
      double precision function rk_R5_H2p_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.3d-8*(300.0d0/Te)**0.4d0
      end function

      ! R6 and R7: the two product channels of H3+ dissociative
      ! recombination, H3+ + e -> H2 + H and H3+ + e -> H + H + H.
      !
      ! HOW THE TWO CONSTANTS ARE BUILT, since neither is printed anywhere
      ! as it stands here.  Larsson, McCall & Orel (2008), Chem. Phys.
      ! Lett. 462, 145, give one total and one branching, both on p. 149:
      ! the thermal rate constant alpha(300 K) = (7.2 +- 1.1)e-8 cm^3 s^-1
      ! of the Kokoouline and Greene calculation, "in good agreement with
      ! the new storage ring results", and the three-body branching ratio
      ! 0.70 +- 0.07, "in very good agreement with the CRYRING storage ring
      ! results".  The two constants below are that total split by that
      ! branching: 7.2e-8 * 0.30 = 2.16e-8 into H2 + H, 7.2e-8 * 0.70 =
      ! 5.04e-8 into three atoms.
      !
      ! THE EXPONENT IS NOT LARSSON'S.  0.65 is the temperature index of
      ! the earlier storage-ring fits of Sundstrom et al. (1994) and Datz
      ! et al. (1995) as Yelle (2004) Table 1 prints them (his R16a, R16b).
      ! Larsson quotes a single 300 K number and no fitted temperature
      ! dependence, so the shape is inherited and only the 300 K value and
      ! the branching are his.  Above about 1000 K the (300/T)^0.65
      ! extrapolation is doing the work, and nothing in either source
      ! measures it there.
      !
      ! WHY NOT THE OLDER PAIR.  Yelle's 2.9e-8 and 8.6e-8, total 1.15e-7,
      ! are 1.60x this total at every temperature, and the excess is
      ! accounted for in the source: Larsson's p. 149 states that "the
      ! early results obtained at CRYRING [23,24] and ASTRID [27], which
      ! gave results just above or at 1e-7 cm^3 s^-1, were slightly too
      ! high because of rotational excitations", and his Ref. 24 is
      ! Sundstrom et al.  The rotationally cold storage-ring measurements
      ! and the calculation that reproduces them supersede that pair, so
      ! the later value is carried.
      double precision function rk_R6_H3p_dr_H2(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.16d-8*(300.0d0/Te)**0.65d0
      end function

      double precision function rk_R7_H3p_dr_3H(Te) result(k)
      real*8, intent(in) :: Te
      k = 5.04d-8*(300.0d0/Te)**0.65d0
      end function

      ! R8: H2+ + H2 -> H3+ + H                 (Theard & Huntress 1974)
      double precision function rk_R8_H2p_H2() result(k)
      k = 2.0d-9
      end function

      ! R9: H2+ + H -> H+ + H2                  (Karpas et al. 1979)
      double precision function rk_R9_H2p_H() result(k)
      k = 6.4d-10
      end function

      ! R10: H+ + H2(v>=4) -> H2+ + H           (Yelle 2004)
      double precision function rk_R10_Hp_H2v4(T) result(k)
      real*8, intent(in) :: T
      k = 1.0d-9*exp(-21900.0d0/T)
      end function

      ! R11: H3+ + H -> H2+ + H2                (Harada et al. 2010)
      double precision function rk_R11_H3p_H(T) result(k)
      real*8, intent(in) :: T
      k = 2.1d-9*exp(-20000.0d0/T)
      end function

      ! R12: H2 + M -> H + H + M   (detailed balance of R15; see below)
      ! Two-body-equivalent: multiply by the H2-EQUIVALENT third-body
      ! density h2_association_collider_density(T, n_H2, n_HI, n_He) (see
      ! the THIRD BODY M note in the module header).
      !
      ! WHAT THIS RATE IS.  Not an independent fit.  R12 and R15 are the same
      ! channel run in opposite directions, so detailed balance genuinely
      ! applies to the pair:
      !     k_diss(T) = k_rec(T) / K_eq(T),   K_eq(T) = n_H2 / n_H^2  [cm^3]
      ! with k_rec = k3b_H_H_to_H2 (cm^6 s^-1, the coefficient of R15 WITHOUT
      ! its third-body density factor) and K_eq = keq_H_H_to_H2.  cm^6 s^-1
      ! divided by cm^3 leaves the cm^3 s^-1 this function returns.
      !
      ! (2) WHY DETAILED BALANCE RATHER THAN THE PUBLISHED DISSOCIATION FIT.
      ! The primary measurement behind the Table-1 dissociation fit was
      ! itself interpreted through detailed balance.  Breshears & Bird (1973)
      ! analyze their shock-tube data with the rate law
      !     -d(A2)/dt = k_d (M) [ (A2) - (A)^2 / K_eq ]
      ! (their Eq. 2) under the relation k_d / k_r = K_eq (their Eq. 3),
      ! introduced with the words "is assumed".  Deriving one direction of
      ! this channel from the other through K_eq is therefore not an
      ! assumption added here; it is the relation the measured coefficient
      ! already carries.
      !
      ! What that leaves is which direction to anchor on, and the Koskinen
      ! Table-1 pair does not overlap in temperature:
      !   * recombination -- Ham, Trainor & Kaufman (1970) J. Chem. Phys. 53,
      !     4395: (8.3 +/- 0.4)e-33 cm^6 s^-1 at 298 K following T^-0.6,
      !     measured over 77-300 K;
      !   * dissociation -- Baulch et al. (1992) J. Phys. Chem. Ref. Data 21,
      !     411, Table 2 and the data sheet on p. 550: 1.5e-9 exp(-48350/T)
      !     cm^3 s^-1 recommended over 2500-8000 K.
      ! The molecular layer modelled here sits at ~1000-2000 K, BETWEEN the
      ! two ranges, and nothing makes the pair thermodynamically consistent
      ! where neither was measured: taken at face value they imply K_eq
      ! values differing from the thermochemical one by up to a factor 10 at
      ! 1000 K.
      !
      ! The recombination is the direction anchored on, because a
      ! recombination coefficient exists whose published range covers the
      ! layer while no dissociation fit does.  R15 therefore carries Cohen &
      ! Westberg's (1983) recommended k1(H2) = 2.8e-31 T^-0.6 over
      ! 50-5000 K rather than Ham's 77-300 K measurement, which it
      ! reproduces in shape and exceeds by 14 % in amplitude (see R15 for the
      ! full argument), and R12 is its thermodynamic reverse.  Reversing the
      ! recombination is also what both reference networks in the workspace
      ! do with the near-identical Baulch p. 495 coefficient (VULCAN reaction
      ! 531, Photochem zahnle_earth "H + H (+ M) <=> H2 (+ M)", both
      ! 2.70e-31 T^-0.6).
      !
      ! WHICH VIBRATIONAL STATE THE COEFFICIENTS REFER TO is answered by the
      ! source rather than inferred.  Breshears & Bird read their initial
      ! rate off the minimum of the postshock density-gradient "dip", i.e.
      ! "the rate corresponding to conditions of vibrational equilibrium in
      ! the thermal sense but with negligible extent of dissociation".  So
      ! the dissociation coefficient is a vibrationally equilibrated rate,
      ! which is the state K_eq assumes on both sides; and at Ham's 77-300 K
      ! H2 is in v = 0 in any case.  The same paper bounds where its rate law
      ! applies: the "observation of induction periods for diatomic
      ! dissociation clearly indicates that the assumptions implicit in Eqs.
      ! (2)-(4) are not valid during the very early course of the reaction;
      ! nevertheless, they may be applicable during the 'steady state' which
      ! follows the induction period."  A hydrodynamic model that resolves
      ! nothing faster than a vibrational relaxation time is in that steady
      ! state throughout.
      !
      ! (1) THE THIRD BODY IS H2.  The anchored coefficient is Cohen &
      ! Westberg's k1(H2), labelled M = H2 on their data sheet and separate
      ! from their k1(H); the studies behind it are H2-as-third-body flow
      ! tubes, Ham et al.'s among them.  Every fit compared against here is
      ! the M = H2 one too: Baulch's "H2 + H2 -> 2H + H2" row, not the
      ! "H2 + Ar" row printed next to it, and Breshears & Bird's k_d^H2.
      ! (Breshears & Bird did not shock pure H2: they measured 5, 10 and
      ! 20 % H2 in Ar and in Xe and extracted k(M = H2) through the linear
      ! composition rule k_d = sum_M X_M k_d^M of their Eq. 4, whose validity
      ! is one of the things their precision was meant to test.)  What is
      ! derived here is therefore k(M = H2), and what the caller multiplies
      ! it by is the H2-EQUIVALENT collider density of
      ! h2_association_collider_density, in which atomic hydrogen counts
      ! k1(H)/k1(H2) and helium k1(Ar)/k1(H2) (see the THIRD BODY M note in
      ! the module header).  Dividing one k1(M) by K_eq and multiplying the
      ! same k1(M) by the density leaves the ratio of the two directions
      ! equal to K_eq for every collider separately, so the pair stays an
      ! exact detailed balance whatever the mixture is.
      !
      ! (3) VERIFICATION.  Measured with this implementation:
      !   k_diss(detailed balance) / k_diss(Baulch 1992 fit) =
      !     0.90 at 8000 K, 1.07 at 5000 K, 0.92 at 3000 K, 0.79 at 2500 K,
      !     0.59 at 2000 K, 0.35 at 1500 K, 0.11 at 1000 K.
      !   The agreement INSIDE Baulch's stated 2500-8000 K range, and the
      !   departure only below it, is the evidence that the construction is
      !   right: where the fit was measured the two agree, where it is an
      !   extrapolation they part company.  Against the PRIMARY measurement
      !   rather than the evaluation -- Breshears & Bird's own M = H2
      !   coefficient 5.48e-9 exp(-52989/T), their abstract's 3.30e15
      !   exp(-105300/RT) cc/mole/s converted -- the same ratio is 1.04 at
      !   3500 K, the bottom of their measured 3500-8000 K range, 0.92 at
      !   4000 K, 0.74 at 5000 K and 0.44 at 8000 K.  Evaluation and
      !   measurement are not the same curve: BB73 divided by Baulch is 0.035
      !   at 1000 K, 0.36 at 2000 K, 0.78 at 3000 K, 1.15-1.45 over
      !   4000-5000 K and 2.05 at 8000 K, matching in the middle and parting
      !   at both ends as a synthesis of several data sets does.
      ! K_eq itself has three independent confirmations:
      !   * the chemical-equilibrium H2 fit this code already carries,
      !     q_h2_equilibrium (lower_column.f90, Koskinen et al. 2022 Eq. 11
      !     quoting Visscher et al. 2006), evaluated in its DILUTE limit
      !     where n_H2 << n_H and q -> 0.831856 / 10^u, with the He/H = 0.0793
      !     the fit is for turning its mixture mixing ratio into n_H2/n_H^2:
      !     that K_eq divided by this one is 0.958 at 2000 K, 0.983 at
      !     2500 K, 1.024 at 3000 K.
      !   * the NIST-JANAF Shomate enthalpies and entropies of H and H2 that
      !     oxygen_rates.f90 already carries: K_eq from dG(T) divided by this
      !     one is 0.9997 at 600 K, 0.9998 at 1000 K, 1.0000 at 2000 K,
      !     1.0001 at 3000 K, 1.0008 at 5000 K.  Two independent routes to
      !     the free energy of H2 -- an observed level ladder summed here, a
      !     calorimetric evaluation there -- agree to better than 0.1 %
      !     wherever the Shomate table is valid.  This is the check that
      !     identified the model ladder retired on 2026-09-06: it deviated
      !     the same way, 1.010 at 600 K rising to 1.101 at 5000 K, which is
      !     the ratio of the two partition functions and nothing else.  An
      !     assertion holds this agreement.
      !   * Cohen & Westberg (1983) recommend BOTH directions on the same
      !     data sheet, so their ratio k1(H2)/k-1(H2) is an equilibrium
      !     constant from an evaluation that used none of the above:
      !     1.958e-25 T^0.1 exp(52530/T) cm^3 over 600-5000 K, which divided
      !     by this K_eq is 1.009 at 600 K, 0.949 at 1000 K, 1.026 at
      !     2000 K, 1.059 at 3000 K and 1.007 at 5000 K.  (It is also the
      !     reason the k_diss above is not simply k-1(H2): the two agree to
      !     those same 12 %, and taking K_eq from the physics instead of from
      !     a two-parameter fit is what keeps R12 and R15 exactly reversible
      !     inside the solver.)
      !
      ! (4) A COMPARISON AGAINST ANOTHER NETWORK IS MEANINGLESS UNLESS THE
      ! THIRD BODY MATCHES.  Breshears & Bird measured all three colliders in
      ! one experiment with one analysis, so their DISSOCIATION coefficients
      ! are the clean statement of how much the third body matters:
      !     M = Ar  1.55e-10 exp(-44736/T)
      !     M = H2  5.48e-9  exp(-52989/T)
      !     M = H   3.52e-9  exp(-43881/T)   [cm^3 s^-1]
      ! Inside their measured 3500-8000 K range k(M = H)/k(M = H2) is 8.7 at
      ! 3500 K and 4.0 at 5000 K, and Cohen & Westberg's independent
      ! RECOMBINATION pair -- k1(H) = 8.8e-33, temperature-independent over
      ! 50-5000 K, against k1(H2) = 2.8e-31 T^-0.6 -- gives 4.2 and 5.2 at
      ! the same two temperatures.  Detailed balance forces the two ratios to
      ! be equal, and where both sources are inside their ranges they are.
      ! They are NOT equal where the sources are extrapolated: at 1000 K the
      ! Breshears ratio reaches 5.8e3 (its two colliders have activation
      ! temperatures 52989 and 43881 K, so the ratio runs away below the
      ! measured range) while Cohen & Westberg give 2.0.  Take 2 to 5 as the
      ! collider spread where H2 exists, and read the 5.8e3 as what
      ! extrapolating an Arrhenius fit costs.  Note also that C&W's M = H
      ! recommendation is itself uncertain by +/-0.5 in log throughout.
      ! A second axis has to match as well: rates quoted as the n -> 0
      ! collisional-dissociation limit, with H2 held in v = 0 (Glover &
      ! Jappsen 2007, ApJ 666, 1, R9 and R10), are ~4e9 below this rate at
      ! 1000 K for the SAME third body, because this one is the vibrationally
      ! equilibrated regime of item (2).
      !
      ! (5) SOURCES, all read from the published papers in references/:
      !   Ham, D. O., Trainor, D. W. & Kaufman, F. 1970, J. Chem. Phys. 53,
      !     4395 (Ham_1970JCP.pdf)
      !   Baulch, D. L. et al. 1992, J. Phys. Chem. Ref. Data 21, 411,
      !     p. 550 and p. 495 (Baulch_1992JPCRD.pdf)
      !   Breshears, W. D. & Bird, P. F. 1973, 14th Symposium (International)
      !     on Combustion, 211 (Breshears_1973.pdf) -- the shock-tube
      !     measurement behind the Baulch dissociation recommendation
      !   Cohen, N. & Westberg, K. R. 1983, J. Phys. Chem. Ref. Data 12, 531,
      !     p. 559 (Cohen_1983JPCRD_12_531.pdf) -- the recombination
      !     coefficient R15 carries, and both directions of the pair
      double precision function rk_R12_H2_thdis(T) result(k)
      real*8, intent(in) :: T
      k = k3b_H_H_to_H2(T)/keq_H_H_to_H2(T)
      end function

      ! R13: H+ + H2 + M -> H3+ + M             (Miller et al. 1968)
      ! Returns the 2-body-equivalent rate: 3.2e-29 * n  [cm^3 s^-1], with n
      ! the third-body density (see the THIRD BODY M note in the header).
      double precision function rk_R13_Hp_H2_M(n) result(k)
      real*8, intent(in) :: n
      k = 3.2d-29*n
      end function

      ! R14: H2 + e -> H + H + e                (Stibbe & Tennyson 1999)
      double precision function rk_R14_H2_edis(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.33d-6*(300.0d0/Te)**0.91d0*exp(-55800.0d0/Te)
      end function

      ! Three-body recombination coefficient of H + H + H2 -> H2 + H2
      ! [cm^6 s^-1], Cohen & Westberg (1983) J. Phys. Chem. Ref. Data 12,
      ! 531, p. 559, recommended k1(H2) over 50-5000 K.  Defined here once
      ! because BOTH directions of this channel read it: R15 multiplies it by
      ! the third-body density, R12 divides it by the equilibrium constant.
      ! Provenance and validity are on R15 below; the reason R12 is built
      ! from it is on R12 above.
      double precision function k3b_H_H_to_H2(T) result(k)
      real*8, intent(in) :: T
      k = 2.8d-31*T**(-0.6d0)
      end function

      ! Three-body recombination coefficient of H + H + H -> H2 + H
      ! [cm^6 s^-1], Cohen & Westberg (1983) J. Phys. Chem. Ref. Data 12,
      ! 531, p. 559, recommended k1(H) over 50-5000 K: 8.8e-33, temperature
      ! independent, k(298) = 8.8e-33 (READ from their recommended-rate
      ! sheet, the same sheet k1(H2) is printed on).  Their stated
      ! uncertainty is +/-0.5 in log k1(H) throughout the range, five times
      ! the +/-0.2 to +/-0.4 of k1(H2).
      !
      ! Atomic hydrogen is 2.0 times as efficient as H2 at 1000 K and 1.5
      ! times at 200 K on these two recommendations; this is the same
      ! evaluation and the same experiments for both, so the ratio is what
      ! the source itself supports and is not built from two sources.
      double precision function k3b_H_H_to_H2_atomic_H() result(k)
      k = 8.8d-33
      end function

      ! Three-body recombination coefficient of H + H + M -> H2 + M for a
      ! MONATOMIC INERT third body [cm^6 s^-1], Cohen & Westberg (1983),
      ! p. 559, recommended k1(Ar) over 77-5000 K: 1.9e-30 T^-1,
      ! k(298) = 6.4e-33, uncertainty +/-0.3 in log throughout (READ, the
      ! argon data sheet, printed there in cm^6 molecule^-2 s^-1).
      !
      ! WHAT IT IS USED FOR HERE: AN ARGON-BASED ESTIMATE FOR HELIUM.
      ! Helium is the dominant third body of a helium-rich base and no
      ! evaluation in hand gives a helium coefficient: the only published
      ! helium calculation, Paolini, Ohlinger & Forrey (2011), Phys. Rev.
      ! A 83, 042713, plots its total over 0-350 K and tabulates nothing
      ! (READ).  Argon is the nearest tabulated monatomic inert third body
      ! of the same evaluation, and its coefficient is used for helium.
      !
      ! IT IS AN ESTIMATE AND NOT A BOUND.  Mass and polarizability order
      ! the two atoms, but they do not order thermally averaged quantum
      ! three-body recombination rates over 77-5000 K: Paolini et al. treat
      ! He and Ar with different interaction potentials and different
      ! resonance and continuum contributions, discuss the sensitivity to
      ! the potential energy surface, distinguish equilibrium from
      ! steady-state populations of the intermediate complex, and find an
      ! additional exchange contribution that matters for Ar at low
      ! temperature and not for He (READ).  No claim k1(He) <= k1(Ar) is
      ! made here.
      !
      ! RANGE AND UNCERTAINTY OF THE NUMBER ITSELF: 77-5000 K, +/-0.3 in
      ! log k1(Ar) throughout, which is Cohen & Westberg's own statement for
      ! argon and is the only uncertainty this coefficient carries; the
      ! distance from argon to helium is not inside it.  Against the H2
      ! coefficient this is 0.43 at 1000 K and 0.86 at 200 K, i.e. it
      ! removes the factor ~2 by which applying k1(H2) to helium overstated
      ! the monatomic colliders.
      !
      ! THE SAME CHOICE IS MADE IN BOTH DIRECTIONS.  R15 and R12 multiply
      ! one collider sum built from this coefficient
      ! (h2_association_collider_density), so whatever the helium estimate
      ! is, the pair stays an exact detailed balance collider by collider
      ! and the estimate cannot displace the equilibrium constant.
      !
      ! IT IS NOT CLAMPED to its stated 77-5000 K range, and neither is
      ! k1(H2), because the grid reaches both ends of the wind and clamping
      ! one direction of the R12/R15 pair without the other would break that
      ! detailed balance.  WHAT THE OUT-OF-RANGE EXTRAPOLATION COSTS IS
      ! MEASURED AND NOT ARGUED FROM THE H2 ABUNDANCE: R15 forms H2 out of
      ! ATOMIC hydrogen at n(H)^2 sum_M k1(M) n_M, so the disappearance of
      ! H2 above the front does not make it vanish.  Measured on the
      ! certified molecular base of LHS 1140 b and on the hottest atomic
      ! wind of the catalog, the association source at every cell with T
      ! outside 77-5000 K stands far below the other H2 formation channels
      ! there.
      double precision function k3b_H_H_to_H2_monatomic(T) result(k)
      real*8, intent(in) :: T
      k = 1.9d-30/max(T, 1.0d0)
      end function

      ! The H2-EQUIVALENT third-body density of H + H + M -> H2 + M and of
      ! its reverse [cm^-3]: the collider sum divided by the H2
      ! coefficient, so that
      !
      !   sum_M k1(M) n_M  =  k3b_H_H_to_H2(T) * h2_association_collider_density
      !
      ! exactly, and R12 and R15 remain an exact detailed-balance pair
      ! collider by collider when both are multiplied by this one density:
      ! the pair is k_rec = k1(M) and k_diss = k1(M)/K_eq with the SAME
      ! k1(M), so the ratio is K_eq whatever the collider mix is.
      !
      ! Three colliders are carried, the three Cohen & Westberg (1983)
      ! recommend a coefficient for on the sheets of this reaction: H2, H
      ! and a monatomic inert atom standing for helium (see THIRD BODY M in
      ! the module header for the substitution and its uncertainty).  Every
      ! other heavy particle is left out rather than counted at the H2
      ! efficiency: no evaluation in hand gives a coefficient for H+, for
      ! the molecular ions or for a metal atom, and in the molecular layer
      ! they are together below 1e-4 of this sum.  Above the H2 front the
      ! ionized colliders are not a small share of the heavy particles, and
      ! leaving them out then understates the R15 association rate there;
      ! what that association rate is worth against the other H2 sources
      ! above the front is measured, not argued from the H2 abundance,
      ! because R15 runs on atomic hydrogen.
      double precision function h2_association_collider_density           &
                                  (T, n_H2, n_HI, n_He) result(n_eff)
      real*8, intent(in) :: T, n_H2, n_HI, n_He
      real*8 :: k_h2
      k_h2  = k3b_H_H_to_H2(T)
      n_eff = max(n_H2, 0.0d0)                                            &
            + max(n_HI, 0.0d0)*k3b_H_H_to_H2_atomic_H()/k_h2              &
            + max(n_He, 0.0d0)*k3b_H_H_to_H2_monatomic(T)/k_h2
      end function

      ! R15: H + H + M -> H2 + M           (Cohen & Westberg 1983)
      ! Returns the 2-body-equivalent rate: 2.8e-31 T^-0.6 * n, with n the
      ! H2-EQUIVALENT third-body density h2_association_collider_density(T,
      ! n_H2, n_HI, n_He), so that the product is the collider sum
      ! k1(H2) n(H2) + k1(H) n(H) + k1(Ar) n(He) (see the THIRD BODY M note
      ! in the header).
      ! Table 1 prints the temperature of this fit as Te; the reaction is a
      ! neutral three-body recombination with no electron in it, so the
      ! heavy-particle T is used.  EXHALE carries a single T, so the two
      ! readings coincide numerically.
      !
      ! WHICH COEFFICIENT, AND WHY NOT KOSKINEN'S.  Koskinen et al. (2022)
      ! Table 1 gives this rate as 8e-33 (300/T)^0.6, the room-temperature
      ! measurement of Ham, Trainor & Kaufman (1970) J. Chem. Phys. 53, 4395:
      ! k_R = (8.3 +/- 0.4)e-33 cm^6 s^-1 at 298 K, following T^-0.6 "fairly
      ! closely", over 77-300 K.  Written as A T^-0.6 that is 2.45e-31, and
      ! the molecular layer here (~1000-2000 K) sits three to seven times
      ! ABOVE the measured range.  What is used instead is the recommended
      ! k1(H2) of Cohen & Westberg (1983) J. Phys. Chem. Ref. Data 12, 531,
      ! p. 559: 2.8e-31 T^-0.6 cm^6 s^-1 over 50-5000 K, k(298) = 9.0e-33.
      ! Its range covers the layer outright, and its exponent is Ham's own
      ! regime carried across that range -- "the low temperature (77-300 K)
      ! data for M = H2 are in reasonable agreement and suggest a temperature
      ! dependence of T^-0.6 for k1(H2)" -- with Ham et al. one of the four
      ! studies the recommendation rests on.  The change of source is 14 % in
      ! amplitude and nothing in shape.  A third, independent evaluation
      ! agrees: Baulch et al. (1992) J. Phys. Chem. Ref. Data 21, 411,
      ! Table 3 and p. 495 recommend ko = 2.7e-31 T^-0.6 over 100-5000 K,
      ! 3.6 % below this one, and it is the coefficient both reference
      ! networks in the workspace carry (VULCAN reaction 531, Photochem
      ! zahnle_earth, both 2.70e-31 T^-0.6).
      !
      ! Stated uncertainty on the value used: +/-0.2 in log k1(H2) at 300 K,
      ! rising to +/-0.4 at 5000 K.
      !
      ! WHY NO HIGH-PRESSURE LIMIT IS CARRIED.  A falloff term would be
      ! needed where k_0 n approaches the high-pressure limit k_inf. The two
      ! reference networks in the workspace disagree about k_inf by three and
      ! a half orders (VULCAN 3.31e-6 T^-1, Photochem 1.0e-12 cm^3 s^-1), so
      ! the onset density n_crit = k_inf/k_0 was evaluated with BOTH.
      ! Measured on the five molecular regression cases and both diagnostic
      ! runs, the largest molecular-layer particle density is 5.4e13 cm^-3
      ! (matrix) and 7.2e14 cm^-3 (He/H = 10 run), against n_crit >=
      ! 4.4e19 cm^-3 from the STRICTER of the two at this k_0. Every
      ! configuration is therefore at least four and a half orders below the
      ! onset, and the two references' disagreement never enters. The
      ! low-pressure form is exact here; no k_inf is carried.
      double precision function rk_R15_3body_H2(T, n) result(k)
      real*8, intent(in) :: T, n
      k = k3b_H_H_to_H2(T)*n
      end function

      ! ------------------------------------------------------------------ !
      ! H + H <-> H2 equilibrium constant and its rovibrational partition
      ! function. See the H2 THERMOCHEMISTRY block above the contains
      ! statement for the constants, the table and its error, and R12 for
      ! what the equilibrium constant is used for.
      ! ------------------------------------------------------------------ !

      ! The internal partition function of H2 is NOT defined here.  It is
      ! h2_partition_function of caloric_eos, the Boltzmann sum over the
      ! observed Roueff et al. (2019) ladder, and this module reads it: one
      ! level set generates both the equilibrium constant below and the
      ! rovibrational internal energy of the equation of state.
      !
      ! WHAT THE LADDER REPLACED, AND BY HOW MUCH.  Until 2026-09-06 this
      ! module summed its own model ladder from the Huber and Herzberg
      ! (1979) constants (we = 4401.213, we xe = 121.336, Be = 60.853,
      ! alpha_e = 3.062 cm^-1), a rigid rotor with a vibration-rotation
      ! coupling and no centrifugal distortion, truncated where a level
      ! passed D0.  Against the observed ladder that model places the levels
      ! too HIGH -- already in v = 0 the J = 1, 2, 3, 5 levels stand +0.21,
      ! +2.21, +9.11 and +56.8 K above the observed ones, the J^4 signature
      ! of the missing distortion term -- so it underpopulates the ladder.
      ! MEASURED, mean internal energy of the model ladder against the
      ! observed one: 0.63 % low at 300 K, 1.8 % at 1000 K, 3.1 % at 2000 K,
      ! 5.1 % at 4000 K, 8.3 % at 8000 K; partition function 0.46 % low at
      ! 300 K, 1.76 % at 1000 K, 3.65 % at 2000 K, 7.7 % at 4000 K.  Of the
      ! 8000 K gap only 3.1 % comes from the 55 levels the model ladder does
      ! not reach; the rest is placement.  K_eq rises by the partition
      ! function ratio, and the dissociation rate R12 falls by it.
      !
      ! Build the ln K_eq table. Called from EXHALE_main before any parallel
      ! region opens; idempotent, so a second call costs nothing.
      !
      !     K_eq(T) = Lam(2 m_H) q_int(T) exp(D0/kT) / [ Lam(m_H) x 4 ]^2
      !     Lam(m)  = (2 pi m kB T / h^2)^(3/2)                    [cm^-3]
      !
      ! q_int is h2_partition_function, measured from v = 0, J = 0, which is
      ! the level the D0 in the exponent is measured from: the zero of the
      ! ladder and the zero of the bond energy are the same level, so the
      ! exponential and the sum do not double count.  The 4 on the H atom is
      ! 2 (electronic doublet) x 2 (nuclear spin I = 1/2); the g_I weights
      ! inside the molecular sum are the matching nuclear-spin counting, so
      ! both sides of the equilibrium are on the same convention and the spin
      ! factors cancel.
      subroutine h2_thermochemistry_init
      integer :: i
      real*8  :: ln_lo, d_ln, t_k, lam_h, lam_h2, hc_k
      ! The largest total spontaneous decay rate over the H2 bound ladder,
      ! which the vibrational heat fraction of the SAME level set needs and
      ! which is reduced from the line list rather than written down.  It
      ! is set here because this is the serial prologue every molecular
      ! configuration already passes through before a parallel region opens
      ! and because both quantities are properties of one level ladder.
      ! Idempotent, like the table below.
      call h2_vibrational_relaxation_init
      if (keq_table_ready) return
      hc_k  = hp_erg*c_light/kb_erg
      ln_lo = log(T_keq_lo)
      d_ln  = (log(T_keq_hi) - ln_lo)/dble(n_keq - 1)
      do i = 1, n_keq
         keq_lnT(i) = ln_lo + d_ln*dble(i - 1)
         t_k        = exp(keq_lnT(i))
         lam_h      = (2.0d0*pi*mu*kb_erg*t_k/hp_erg**2)**1.5d0
         lam_h2     = (2.0d0*pi*2.0d0*mu*kb_erg*t_k/hp_erg**2)**1.5d0
         keq_lnK(i) = log(lam_h2) + log(h2_partition_function(t_k))       &
                    + D0_H2_cm*hc_k/t_k - 2.0d0*log(4.0d0*lam_h)
      enddo
      keq_table_ready = .true.
      end subroutine h2_thermochemistry_init

      ! K_eq(T) = n_H2 / n_H^2 [cm^3] at chemical equilibrium, interpolated
      ! from the table. T is clamped into [100, 20000] K (see the ENDS note
      ! above the contains statement).
      double precision function keq_H_H_to_H2(T) result(keq)
      real*8, intent(in) :: T
      real*8  :: ln_t, frac
      integer :: i
      ! The table is built ONCE, serially, by h2_thermochemistry_init (the
      ! main program calls it before any parallel sweep; a test driver that
      ! uses this function calls it first). A lazy build here was removed on
      ! 2026-09-13: the readiness flag was read outside the
      ! critical region that built the table, which is a data race for any
      ! caller that reaches this function concurrently before the build.
      if (.not. keq_table_ready) then
         write(*,*) '(mol_rates) ERROR: keq_H_H_to_H2 called before'//   &
                    ' h2_thermochemistry_init. Call the initializer'//    &
                    ' serially first. Aborting.'
         error stop 1
      endif
      ln_t = log(min(max(T, T_keq_lo), T_keq_hi))
      i    = int((ln_t - keq_lnT(1))/(keq_lnT(2) - keq_lnT(1))) + 1
      i    = min(max(i, 1), n_keq - 1)
      frac = (ln_t - keq_lnT(i))/(keq_lnT(i+1) - keq_lnT(i))
      keq  = exp(keq_lnK(i) + frac*(keq_lnK(i+1) - keq_lnK(i)))
      end function keq_H_H_to_H2

      ! R16: HeH+ + e -> He + H                 (Yousif & Mitchell 1989)
      ! 3.4x below the Florescu-Mitchell & Mitchell (2006) values tabulated
      ! by Garcia Munoz (2025) Table A.6 at 500 K and 8.6x below at 10^4 K:
      ! both fall with T, this one faster.  See the HE-DOMINATED LIMIT note
      ! in the header.
      double precision function rk_R16_HeHp_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.0d-8*(300.0d0/Te)**0.6d0
      end function

      ! R17: He+ + H2 -> H+ + H + He
      !      (measured two-body total: Boehringer & Arnold 1986;
      !       rising branch: the Moses & Bass 2000 ESTIMATE -- see below)
      !
      ! THE DISSOCIATIVE CHARGE TRANSFER, AND IT HAS TWO MECHANISMS.  The
      ! Koskinen Table-1 entry for it, 1e-9 exp(-5700/T) (Moses & Bass 2000),
      ! is an Arrhenius form with a 0.49 eV barrier taken from a model of
      ! SATURN's ionosphere (Moses, J. I. & Bass, S. F. 2000, J. Geophys.
      ! Res. 105, 7013; bibliographic record READ through ADS, the paper
      ! itself not obtained), and BY ITSELF it is
      ! wrong at every temperature a molecular base layer reaches: it gives
      ! 5.6e-18 at 300 K where the reaction has been measured three times at
      ! ~1.1e-13, and 6.5e-16 at 400 K against a measured ~1.35e-13.  The
      ! reason the measurements do not switch off is stated by the
      ! experimenters: the reactants and products do not correlate
      ! adiabatically (Mahan 1971), and what carries the reaction at thermal
      ! energy is TUNNELLING through a small barrier out of a long-lived
      ! He+-H2 complex (Preston et al. 1978; Boehringer & Arnold 1986,
      ! Discussion).  That channel has a weak, slightly NEGATIVE temperature
      ! dependence.  Over the barrier, above ~400 K, an Arrhenius branch
      ! takes over.  The two are different mechanisms of one reaction and
      ! they ADD.
      !
      ! WHAT IS MEASURED, read from the published papers in references/:
      !
      !   Boehringer, H. & Arnold, F. 1986, J. Chem. Phys. 84, 1459
      !     (Bohringer_1986JCP_84_1459.pdf), selected-ion drift tube,
      !     18 <= T <= 408 K, He+ in pure H2, two-body and three-body
      !     separated by the pressure dependence:
      !         k2 = 1.1e-13 (300/T)^0.24(+/-0.04)  cm^3 s^-1
      !         k3 = 1.6e-30 (100/T)^1.27(+/-0.4)   cm^6 s^-1
      !     with "At temperatures below 80 K the rate coefficient is no
      !     longer increasing and an upper limit of k2 < 2e-13 cm^3 s^-1 can
      !     be given for the studied temperature range."
      !
      !   Johnsen, R., Chen, A. & Biondi, M. A. 1980, J. Chem. Phys. 72,
      !     3085 (Johnsen_1980JCP_72_3085.pdf), drift tube, 78 to 330 K,
      !     their Table I:
      !         330 K   k = (1.1 +/- 0.1)e-13,  K = (4.4 +/- 2.0)e-31
      !          78 K   k = (1.5 +/- 0.15)e-13, K = (1.8 +/- 0.4)e-30
      !     Their Fig. 3, whose title IS this reaction, carries the curve
      !     above 330 K as an EFFECTIVE temperature from their earlier
      !     elevated-ion-energy data: READ from it, ~1.45e-13 at 100 K, a
      !     minimum ~1.05e-13 at 300 K, ~1.35e-13 at 400 K, ~2.0e-13 at
      !     500 K, ~2.5e-13 at 600 K and ~3.0e-13 at 700 K.
      !
      !   Schauer et al. 1989 (see the retired-R20 block) extend the same
      !     reaction to 15-40 K, 3.0e-14 to 4.9e-14, and their Fig. 3 is
      !     where the two papers above were first seen together.
      !
      ! THE FORM USED, and why it is a sum and not a choice:
      !
      !     k(R17) = [ k2(T) - k(R23) ]  +  1e-9 exp(-5700/T)
      !
      ! The first bracket is the measured two-body TOTAL of He+ + H2 less
      ! the radiative branch this network carries separately as R23, so the
      ! coded channels sum to the measured total instead of each being
      ! quoted independently; the branching is then stated once, here, and
      ! cannot drift from R23.  (Johnsen et al. attribute at least 80
      ! percent of the total to this dissociative channel, and Boehringer &
      ! Arnold decline to fix the branching more tightly than that, so
      ! "total less the radiative branch" is the assignment the measurements
      ! support.)  The second term is the Koskinen Table-1 Arrhenius, kept
      ! because it is the over-barrier branch and it is what the data above
      ! 400 K show.
      !
      ! VERIFICATION against Johnsen et al.'s Fig. 3, the only measurement
      ! of the rising branch: the sum gives 1.03e-13 at 300 K (read 1.05e-13),
      ! 9.6e-14 at 400 K (1.35e-13), 1.01e-13 at 500 K (2.0e-13), 1.60e-13 at
      ! 600 K (2.5e-13) and 3.7e-13 at 700 K (3.0e-13) -- inside a factor 2
      ! everywhere, against the 200x and 2e4x the Arrhenius alone was out by
      ! at 400 and 300 K.
      !
      ! VALIDITY RANGE, and what is extrapolation.  The measured bracket is
      ! 18-408 K (Boehringer & Arnold) and 78-330 K (Johnsen et al.); above
      ! 408 K it is carried up with its own weak exponent, which is mild
      ! (a factor 0.82 from 300 K to 1800 K) and, above ~625 K where the
      ! Arrhenius branch passes it, no longer controls the rate.  Above
      ! ~700 K NOTHING is measured, and the Arrhenius branch then decides the
      ! rate over the whole hot part of a molecular layer while being an
      ! extrapolation of a fit whose own support is effective-temperature
      ! drift data below 700 K.  That is a MODEL UNCERTAINTY of this network
      ! and not an established rate: a result that depends on the He+ + H2
      ! destruction of H2 above 700 K inherits it.  Carried as published
      ! because no measurement exists to put in its place.  Below 80 K the measured rate stops
      ! rising and k2 < 2e-13 is an upper limit; this form gives 2.1e-13 at
      ! 18 K, i.e. it sits at that limit, and no molecular layer of this
      ! code is that cold.
      !
      ! THE THREE-BODY CHANNEL IS NOT CARRIED, AND THE OMISSION IS GUARDED
      ! CELL BY CELL RATHER THAN ARGUED FROM A REPRESENTATIVE DENSITY.  Both
      ! papers measure one, He+ + 2H2 -> products, with
      ! k3 = 1.6e-30 (100/T)^1.27 cm^6 s^-1 (Boehringer & Arnold 1986); at
      ! the 1e16-1e17 cm^-3 of their drift tubes it is what the pressure
      ! dependence they report IS.  Whether it may be dropped is a statement
      ! about k3 n / k2 in the cell being integrated and nowhere else, so
      ! that ratio is formed here whenever the caller supplies the density,
      ! and a cell above THREE_BODY_WARN_RATIO is counted and reported
      ! through hep_h2_three_body_domain rather than passing silently.
      !
      ! Garcia Munoz (2025) Table A.7 gives this channel as 3e-14 and
      ! T-independent, which is the Schauer et al. 15-40 K value carried
      ! flat; it is below both drift-tube measurements at every temperature
      ! above 40 K.
      double precision function rk_R17_Hep_H2_diss(T, n_M) result(k)
      real*8, intent(in) :: T
      ! Third-body density of the cell [cm^-3], supplied by the callers that
      ! have it.  It changes no rate: it is read only by the domain guard
      ! below, which is why it is optional and why a diagnostic caller that
      ! has no density may leave it out.
      real*8, intent(in), optional :: n_M
      real*8 :: k_two_body_total, k_tunnelling, ratio
      ! The measured two-body total of He+ + H2 (Boehringer & Arnold 1986).
      k_two_body_total = 1.1d-13*(300.0d0/T)**0.24d0
      ! Less the radiative branch the network carries as R23, so that the
      ! coded channels sum to that total.  The subtraction cannot go
      ! negative for any temperature this code reaches (the total is
      ! 3.1e-14 even at 1e5 K against R23's 7.2e-15), and the guard is
      ! there so that a future change of R23 fails loudly in the comparison
      ! above rather than quietly as a negative rate.
      k_tunnelling = max(k_two_body_total - rk_R23_H2_Hep_cx(), 0.0d0)
      k = k_tunnelling + 1.0d-9*exp(-5700.0d0/T)
      ! THE DOMAIN GUARD, in this cell's own n and T.
      if (present(n_M)) then
         ratio = 1.6d-30*(100.0d0/T)**1.27d0*n_M/max(k, 1.0d-300)
         !$omp critical (hep_h2_three_body_guard)
         n_hep_h2_evaluations = n_hep_h2_evaluations + 1
         if (ratio .gt. hep_h2_worst_ratio) hep_h2_worst_ratio = ratio
         if (ratio .gt. three_body_warn_ratio) then
            n_hep_h2_three_body = n_hep_h2_three_body + 1
            if (.not. hep_h2_warned) then
               hep_h2_warned = .true.
               write(*,'(A,ES10.3,A,F8.1,A,ES10.3)')                      &
                  ' (mol_rates) WARNING: the He+ + H2 three-body channel'//&
                  ' is not negligible here: k3 n / k2 =', ratio,          &
                  ' at T =', T, ' K, n =', n_M
               write(*,'(A)') '   R17 carries the two-body rate only'//   &
                  ' (Boehringer & Arnold 1986 measure both); above a'//   &
                  ' ratio of 0.1 the omission is a modelling error and'// &
                  ' not a rounding one.'
            endif
         endif
         !$omp end critical (hep_h2_three_body_guard)
      endif
      end function

      ! How the He+ + H2 two-body-only assumption fared over the run: how
      ! many evaluations were made, how many exceeded the warn ratio, and
      ! the worst k3 n / k2 seen.  A counter of a stated approximation, so
      ! that "it was negligible" is a measurement and not an assertion.
      subroutine hep_h2_three_body_domain(n_eval, n_above, worst)
      integer, intent(out) :: n_eval, n_above
      real*8,  intent(out) :: worst
      n_eval  = n_hep_h2_evaluations
      n_above = n_hep_h2_three_body
      worst   = hep_h2_worst_ratio
      end subroutine

      ! R18: HeH+ + H2 -> H3+ + He              (Bohme et al. 1980)
      double precision function rk_R18_HeHp_H2() result(k)
      k = 1.5d-9
      end function

      ! R19: HeH+ + H -> H2+ + He               (Karpas et al. 1979)
      double precision function rk_R19_HeHp_H() result(k)
      k = 9.1d-10
      end function

      ! ------------------------------------------------------------------ !
      ! R20 OF KOSKINEN TABLE 1 IS NOT IN THIS MODULE, and the reason is the
      ! paper it cites.  Table 1 gives He+ + H2 -> HeH+ + H as 4.2e-13,
      ! T-independent, citing Schauer, Jefferts, Barlow & Dunn (1989),
      ! J. Chem. Phys. 91, 4593 (references/Schauer_1989JCP_91_4593.pdf,
      ! READ).  That paper measures the RADIATIVE (He+ + H2 -> He + H2+ + hv,
      ! this network's R23) and the DISSOCIATIVE (He+ + H2 -> He + H+ + H,
      ! R17) charge transfer over 15 < T < 40 K.  It does not measure an
      ! HeH+ channel, and it rules that channel out at thermal energies:
      !
      !   "Reaction (3) apparently does not become allowed until the
      !   collision energy approaches 9 eV."   (Sec. I, reaction (3) being
      !   He+ + H2 -> H + HeH+; the reference given is Jones, Wu, Hughes &
      !   Hopper 1980, J. Chem. Phys. 73, 5631.)
      !
      !   "Both H2+ and HeH+ react rapidly with H2 to form H3+, so
      !   observation of this product serves to monitor the sum of
      !   reactions (2) and (3). Because reaction (3) is apparently ruled
      !   out, as just discussed, we nominally attribute observation of H3+
      !   to reaction (2), and refer to it as radiative charge transfer
      !   (RCT)."
      !
      ! The paper's H3+ signal is therefore the SUM of the radiative charge
      ! transfer and the HeH+ channel, and its Fig. 4 gives that sum as
      ! 7.5e-15 to 1.1e-14 over 15-40 K, with an upper limit of 1.0e-14 from
      ! 70 to 320 K (Johnsen et al. 1980, drawn in the same figure).  The
      ! measurement therefore BOUNDS the HeH+ channel at <= 1.0e-14, i.e. 42
      ! times below the 4.2e-13 Table 1 assigns to it, and the reason the
      ! paper gives for the bound is structural: the reactants and products
      ! of He+ + H2 -> HeH+ + H do not correlate (their Fig. 1, adapted from
      ! Hopper 1980) and the channel opens near 9 eV of collision energy,
      ! which is 1e5 K in relative kinetic energy -- decades above anything
      ! this code integrates.  The channel is closed here, so it is not
      ! carried; a coefficient that is identically zero over the whole
      ! domain is not a rate.
      !
      ! WHAT REMOVING IT DOES, measured rather than argued
      ! (molecular_scalar_gj1132_kzz1e9 at He/H = 2.13, first cell, with the
      ! ionization re-solved on both sides).  At the SHARE level R20 looked
      ! decisive: it carried 25 to 30 percent of the H2 loss and, with R18
      ! consuming the HeH+ it made, about half of it.  With the balance
      ! actually re-solved it is not, and the reason is worth stating:
      ! helium ions in that layer are SINK-limited, so removing one He+ sink
      ! raises n(He+) until the others carry the same flux.
      !
      !   at 1023 K   n(He+) x 1.11, n(HeH+) / 487, n(H2+) x 0.71,
      !               total H2 loss  1.916e+03 -> 1.904e+03  (-0.7 %)
      !   at  629 K   n(He+) x 4.4,  n(HeH+) / 8.0e4, n(H+) x 3.9,
      !               total H2 loss  4.455e+03 -> 4.443e+03  (-0.3 %)
      !
      ! So the repair rewires the HeH+ and H+ budgets of the base layer by
      ! four to five decades and by a factor 4, and leaves the NET H2
      ! destruction there where it was.  It is made because the channel is
      ! not supported by its source, not because of what it moves.
      ! ------------------------------------------------------------------ !

      ! H2+ + He -> HeH+ + H                             (Black 1978)
      !
      ! THE HeH+ SOURCE THE NETWORK NEEDS, and the one Koskinen Table 1 does
      ! not carry.  With the unsupported R20 gone, Table 1 has no HeH+
      ! formation channel left at all, and a species with destruction
      ! channels and no formation is not a model.  This is the published
      ! route that needs no He+ and is therefore the one that matters in a
      ! helium-rich gas.  It is ON by default for that reason: the
      ! alternative is a network whose HeH+ row has no source term.
      !
      ! SOURCE, read from the published paper: Black, J. H. 1978, ApJ 222,
      ! 125, "Molecules in planetary nebulae" (references/
      ! Black_1978ApJ_222_125.pdf), his reaction (12) and the rate table of
      ! that paper, his reference 13 being the Chupka, Berkowitz & Russell
      ! (1969) cross sections.  The same coefficient is carried by Garcia
      ! Munoz (2025), A&A 698, A199, Table A.6 and network row 197, which is
      ! where the omission was noticed.
      !
      ! VALIDITY, in the source's own words (p. 126, verbatim):
      !
      !   "is also a significant source of HeH+. This reaction is
      !   substantially endothermic for ground-state H2+; however, it is
      !   rapid for H2+ with v > 3, states which may be well populated at
      !   T > 5000 K. Cross sections for reaction (12) have been determined
      !   by Chupka, Berkowitz, and Russell (1969) for specific v. When
      !   these cross section data are integrated over a Maxwellian velocity
      !   distribution and averaged over a thermal distribution of
      !   vibrational populations of H2+, we find an approximate
      !   representation of the rate coefficient
      !       k12 ~ 3 x 10^-10 exp(-6717/T) cm^3 s^-1.
      !   This result is in harmony with the rate measured by Neynaber and
      !   Magnuson (1973)."
      !
      ! and, on the assumption that carries it:
      !
      !   "In what follows, the populations of vibrational states of H2+ are
      !   assumed to be thermalized at the kinetic temperature."
      !
      ! WHICH WAY THE ASSUMPTION ERRS, in the source's own words rather than
      ! by argument here (p. 126):
      !
      !   "If, however, the rate of reaction (10) is as large as suggested
      !   recently by Bottcher (1976), then some excited vibrational states
      !   may have nonthermal populations, and the predicted abundances of
      !   H2+ and HeH+ may be overestimates."
      !
      ! So this coefficient is an APPROXIMATION VALID UNDER AN ASSUMED
      ! THERMAL H2+ VIBRATIONAL DISTRIBUTION -- not a bound in either
      ! direction -- and the sentence above states only which way it would
      ! err if that distribution is not thermal.  WHETHER IT IS thermal here
      ! is a comparison of time scales, and this layer's is not favourable:
      ! at the base of the LHS 1140 b molecular case the H2+ lifetime
      ! against R5, R8 and R9 is 1/(k5 n_e + k8 n(H2) + k9 n(H I)), of order
      ! 1e-3 s, while the radiative vibrational relaxation of H2+ is of
      ! order 1 s and the collisional one is carried by the same collisions
      ! that destroy the ion.  The ion therefore reacts long before its
      ! vibrational population is set by T, and that population is whatever
      ! the process that made it left behind.  The coefficient is used
      ! OUTSIDE the assumption it was derived under; the error may run
      ! either way and its size is not known without a v-resolved H2+,
      ! which this code does not carry.
      ! Values: 4.4e-16 at 500 K, 1.0e-11 at 2000 K, 7.8e-11 at 5000 K,
      ! 1.5e-10 at 10^4 K.
      !
      ! THE OTHER TWO HeH+ SOURCES OF THAT PAPER, and why only one is here.
      ! Black's rate table lists three:
      !   (12) H2+ + He -> HeH+ + H          3.0e-10 exp(-6717/T)  -- this
      !   (13) He* + H2 -> HeH+ + H + e      1e-9                  -- carried
      !        already, as the associative branch (1 - f_penning_HeI23S) of
      !        the He(2^3S) + H2 ionization in row 7 of mol_heh_rows, with
      !        the Garcia Munoz (2025) coefficient and branching rather than
      !        Black's single 1e-9
      !   (14) H+ + He -> HeH+* -> HeH+ + hv 1e-18                 -- NOT
      !        carried.  It is not a measurement: Black introduces it as a
      !        PROPOSAL ("Dabrowski and Herzberg (1977) have proposed that
      !        there is a substantial probability of forming HeH+ by
      !        vibrational inverse predissociation") argued from a count of
      !        quasi-bound levels ("The number of quasi-bound levels of HeH+
      !        and its large dipole moment suggest a large rate coefficient"),
      !        and adopting an unmeasured 1978 estimate into a channel is
      !        what put R20 in this network in the first place.
      !        WHAT IT WOULD DO, measured so the omission is bounded and not
      !        merely noted (first cell at 1023 K): 1e-18 n(H+)
      !        n(He) = 55 cm^-3 s^-1 against 0.15 for reaction (12), so it
      !        would be the dominant HeH+ source and would raise n(HeH+)
      !        from 7.7e-05 to ~3e-02 cm^-3.  It changes nothing observable:
      !        HeH+ would still be 1e-13 of the helium, and the H2 loss it
      !        drives through R18 would be 0.8 percent of the total instead
      !        of 0.002.  Resolving it needs a computed radiative-association
      !        rate for H+ + He, which exists (e.g. the He+ + H channel of
      !        Courtney et al. 2021, ApJ 919, 70) and was not obtained here.
      double precision function rk_H2p_He_HeHp(T) result(k)
      real*8, intent(in) :: T
      k = 3.0d-10*exp(-6717.0d0/T)
      end function

      ! R23: H2 + He+ -> H2+ + He + hv          (Barlow 1984)
      !
      ! THE PRODUCTS INCLUDE A PHOTON, and the name of the channel says so:
      ! this is RADIATIVE charge transfer.  It agrees with the direct
      ! measurement: Schauer, Jefferts, Barlow & Dunn (1989), J. Chem. Phys.
      ! 91, 4593, Fig. 4 (READ) gives 7.5e-15 to 1.1e-14 over 15-40 K with
      ! no temperature dependence, and an upper limit of 1.0e-14 over
      ! 70-320 K from Johnsen et al. (1980) drawn in the same figure.  The
      ! coded 7.2e-15 is that measurement's own group (Barlow's thesis is
      ! reference 1 of the 1989 paper) and sits inside it.  NOTE what the
      ! measurement actually bounds: the 1989 experiment detects H3+, which
      ! both H2+ and HeH+ make, so its number is the sum of this channel and
      ! the HeH+ channel of the retired R20 -- which is why that channel
      ! cannot be 4.2e-13 (see the retired-R20 block above).
      !
      ! WHERE THE 9.16 eV GOES, and why the heat ledger may not have it all.
      ! Boehringer & Arnold (1986), J. Chem. Phys. 84, 1459, p. 1461,
      ! describing Hopper's mechanism for this channel (READ, verbatim):
      !
      !   "Hopper reconsidered the reaction dynamics in more detail and
      !   suggested that also a reaction channel leading to H2+ should be
      !   possible. This process involves a radiative transition from the
      !   first formed excited state of the collision complex (He+.H2) to
      !   the ground state (He.H2+) which then decays into the products
      !   (1a). The wavelength of the emitted photon should be about 153 nm
      !   and H2+ should preferentially be produced in a vibrationally
      !   excited state (v = 2)."
      !
      ! 153 nm is 8.103 eV of the reaction's 9.161 eV, and it LEAVES: at
      ! 1530 A it is longward of the Lyman-Werner bands (912-1201 A) and of
      ! the Lyman continuum, so neither H2 nor H I absorbs it where it is
      ! made.  What stays with the gas is the remaining 1.058 eV, and the
      ! paper says most of that is INTERNAL -- the H2+ is born in v = 2,
      ! about 0.55 eV on the H2+ ladder -- so only what is collisionally
      ! relaxed before the ion reacts becomes heat.
      !
      ! THE APPROXIMATION MADE, with its range: the heat ledger deposits
      ! q(R23) - E_photon and nothing is re-absorbed, because this code
      ! carries no transfer for a 1530 A line and no v-resolved H2+ to relax.
      ! That OVERSTATES the deposit by at most the 0.55 eV of vibration
      ! (half of what is left), and it understates nothing.  A model that
      ! followed the H2+ vibration would place that 0.55 eV in whichever
      ! reaction destroys the ion instead.
      double precision function rct_photon_energy_eV() result(e)
      ! The photon R23 emits, from the wavelength Boehringer & Arnold quote
      ! for Hopper's mechanism.  Read by the reaction-heat ledger, which is
      ! the only consumer, so the 153 nm is written down once.
      e = 1.239841984d3/153.0d0
      end function
      double precision function rk_R23_H2_Hep_cx() result(k)
      k = 7.2d-15
      end function

      ! End of module
      end module mol_rates

      module oxygen_rates
      ! Rate coefficients and thermodynamic data for the oxygen-hydrogen
      ! chemistry of the A2 option (docs/a2_oxygen_option_design.md).
      !
      ! THIS MODULE IS NOT WIRED INTO THE CODE.  It is the milestone-M1
      ! deliverable: coefficients, their sources and their validity ranges,
      ! plus the thermodynamic data that generates the reverse rates.  It is
      ! deliberately absent from SRC in the Makefile; M2 connects it.  The
      ! standalone driver that checks it is src/tests/a2_m1/.
      !
      ! The audit that fixed every number below, with the two-network
      ! comparison table and the measured share each reaction carries at the
      ! HD 189733 b base, is docs/a2_reaction_audit.md.
      !
      ! Units are cgs throughout: two-body rates in cm^3 s^-1, three-body
      ! rates in cm^6 s^-1 (the low-pressure limit k0), T in K.  EXHALE
      ! carries a single temperature, so no distinction is made between the
      ! heavy-particle and the electron temperature; none of these reactions
      ! has an electron in it.
      !
      ! ---------------------------------------------------------------
      ! SOURCES AND THE TWO-NETWORK RULE
      !
      ! Every coefficient is transcribed from a named publication, and every
      ! transcription is checked against two independent ports of the same
      ! literature that are on this machine:
      !
      !   Photochem  photochem/data/reaction_mechanisms/zahnle_earth.yaml
      !              (Kevin Zahnle's earth network, maintained by Wogan;
      !              k = A T^b exp(-Ea/T), Ea in K; provenance in the ref:
      !              keys, resolved through photochem/data/bib.bib)
      !   VULCAN     VULCAN/thermo/NCHO_photo_network.txt
      !              (k = A T^B exp(-C/T); odd ids forward, even ids the
      !              thermodynamic reverse)
      !
      ! Where the two disagree the disagreement is recorded here with its
      ! size, in the way mol_rates.f90 records the R16-R20 discrepancies.
      !
      ! ---------------------------------------------------------------
      ! REVERSE RATES: THERMODYNAMIC REVERSAL, NOT A SECOND TRANSCRIPTION
      !
      ! O1, O2, O8, O9 and O10 all run close to cancellation on a hot base
      ! (the HD 209458 b gross rates are three decades above the net), and an
      ! independently transcribed reverse rate then produces an arbitrary net
      ! rather than a small error -- which is exactly the defect the R21/R22
      ! note in mol_rates.f90 records for the He charge-exchange pair.  Both
      ! reference networks avoid it the same way: neither carries an explicit
      ! reverse for any of these reactions, and each computes
      !
      !     k_rev = k_fwd / K_c ,   K_c = exp(-dG/RT) (P0/kB T)^(dn) ,
      !
      ! with P0 = 1 bar = 1e6 dyn cm^-2 and dn = n_products - n_reactants.
      ! That is rate_from_detailed_balance below.  Detailed balance then
      ! holds by construction and the hot limit reduces to chemical
      ! equilibrium.
      !
      ! ---------------------------------------------------------------
      ! THERMODYNAMIC DATA: NIST SHOMATE
      !
      ! Shomate coefficients A-G for H, H2, O, OH, H2O, CO and C,
      ! transcribed from the species: block of zahnle_earth.yaml and checked
      ! digit by digit against the NIST Chemistry WebBook (SRD 69) gas-phase
      ! thermochemistry pages, whose source is Chase (1998), NIST-JANAF
      ! Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9.
      !
      ! Atomic C was added at M2, for the CO <-> C + O equilibrium that sets
      ! the carbon monoxide reservoir (co_equilibrium_density below).  Like
      ! atomic O it has no Shomate table on the WebBook, so it is checked
      ! against the WebBook's NIST-JANAF gas-phase entry instead: the
      ! coefficients return dfH(298.15 K) = 716.672 kJ/mol against the
      ! tabulated 716.68 and S(298.15 K) = 158.102 J/mol/K against 158.100.
      ! The eighth NIST coefficient H is not carried: NIST sets
      ! H = dfH(298.15 K), so dropping it makes the Shomate enthalpy the
      ! absolute enthalpy including the enthalpy of formation, which is what
      ! a reaction dG needs.
      !
      ! CAVEAT -- THE OH ENTHALPY OF FORMATION.  The NIST-JANAF value carried
      ! here is dfH(OH, 298.15 K) = 38.99 kJ/mol.  Burcat's NASA-9 set, which
      ! VULCAN uses (VULCAN/thermo/NASA9/OH.txt), gives 37.28 kJ/mol.  The
      ! 1.71 kJ/mol offset is a pure enthalpy-of-formation difference (the
      ! entropies agree to 0.03 J/mol/K and every other species in this set
      ! agrees to better than 0.02 kJ/mol up to 1500 K; above that the two
      ! H2O fits part company by 0.25 kJ/mol at 2000 K and 0.62 at 2500 K,
      ! which is a fitting difference in the top interval, not a difference
      ! of thermochemistry), and it multiplies any reverse
      ! rate in which OH appears by exp(1.72 kJ/mol / RT) per OH: 1.99 at
      ! 300 K, 1.23 at 1000 K, 1.11 at 2000 K.  The JANAF value is kept
      ! because the rate evaluations these coefficients reverse (Baulch
      ! et al. 1992; Lifshitz & Michael 1991) were themselves referred to
      ! JANAF thermochemistry, and the measured test of the reversal in
      ! docs/a2_reaction_audit.md sec. 5 confirms it: reversing O10 with the
      ! Shomate data reproduces VULCAN's independently transcribed
      ! O + H2O -> OH + OH to within 0.82-1.36 over 300-2500 K, against
      ! 1.60-3.27 with the NASA-9 data.  A run that needs the modern
      ! (Active Thermochemical Tables) OH enthalpy must change the OH F
      ! coefficient below by -1.71 and say so.
      !
      ! ---------------------------------------------------------------
      ! THIRD BODY M (O8, O9)
      !
      ! Both three-body reactions are given as the low-pressure limit k0 with
      ! a Lindemann falloff onto k_inf; the reference networks apply no
      ! collider efficiencies to them, so every heavy particle counts equally
      ! and the same caveat mol_rates.f90 states for R12/R13/R15 applies: a
      ! monatomic third body has no internal modes to take up the released
      ! energy, so He is in general less efficient than H2 and these rates
      ! with He as M are upper bounds.  At the base densities of interest
      ! both reactions are deeply in the low-pressure limit -- the measured
      ! reduced pressure at the HD 189733 b base (n = 7.87e12 cm^-3, 864 K)
      ! is 1.2e-8 for O8 and 2.6e-9 for O9 -- so k_inf and the choice of
      ! falloff function are numerically irrelevant there.  They are carried
      ! anyway because the option is meant to run down to the 1e-4 bar level.
      !
      ! ---------------------------------------------------------------
      ! VALIDITY RANGE OF THE SET AS A WHOLE
      !
      ! This set is valid where the base H2/H partition is kinetic and set by
      ! the oxygen cycle: the regime measured on HD 189733 b (864 K, oxygen
      ! family 96-99.6% of the net H2 destruction).  Where the base is hot
      ! enough that H2 + M -> H + H + M runs the net (HD 209458 b, 2331 K,
      ! 70-97%), these reactions add a few percent and the R12/R15 pair of
      ! mol_rates.f90 is the physics.  Where the base is sulfur-rich, the
      ! set is missing S + H2 <-> H + HS, measured at 19% of the net at the
      ! HD 209458 b 1-microbar level; sulfur is excluded by decision D8 of
      ! the design.
      !
      ! ---------------------------------------------------------------
      ! REACTIONS NOT IN THIS MODULE, AND WHY
      !
      ! O8, O9, O10 and O12 below are transcribed but are NOT part of the
      ! minimal set the measured budget supports; they are here so that the
      ! audit records a value and a verdict rather than an omission, in the
      ! way mol_rates.f90 carries R21/R22.  Their measured shares of the
      ! dominant channel O1 at the HD 189733 b base are 7.8e-9 (O8),
      ! 9.7e-8 (O9), 5.5e-7 (O10) and 2.8e-5 (O12).
      !
      ! O11 of the design -- physical quenching O(1D) + M -> O + M -- has no
      ! entry here, and the reason is now published rather than assumed.
      ! The IUPAC evaluation that supplies O6 quotes a single rate for
      ! O(1D) + H2 that is the SUM of the reactive and the quenching
      ! channel, and states that the reactive channel is above 95% of it
      ! (Atkinson et al. 2004, data sheet I.A2.18).  So quenching by H2 is
      ! inside O6's own coefficient, is under 5% of it, and has no
      ! separately evaluated rate to transcribe.  Neither reference network
      ! carries a quenching channel for H2 or for He, the two colliders that
      ! make up 85% of this gas; and O6 alone gives O(1D) a lifetime of
      ! 1.6e-3 s at the HD 189733 b base, which is
      ! shorter than every other time scale in the problem by many decades.
      ! The local steady state for O(1D) is therefore closed by O6.
      !
      ! No ion chemistry of the water family (H3O+, OH+, H2O+, O+ + H2): the
      ! trace-metal ionization solver assumes coronal-like rates and no
      ! three-body reactions, neither of which holds for molecules, and the
      ! measured oxygen at these levels is 100% neutral.

      implicit none
      private

      public :: rk_O1_OH_H2_water,    rk_O2_O_H2_hydroxyl,               &
                rk_O6_O1D_H2_hydroxyl, rk_O8_O_H_assoc,                  &
                rk_O9_H_OH_assoc,     rk_O10_OH_OH_water,                &
                rk_O12_O1D_H2O_hydroxyl,                                 &
                lindemann_rate
      public :: enthalpy_shomate, entropy_shomate, gibbs_energy_shomate,  &
                equilibrium_constant_conc, rate_from_detailed_balance,    &
                co_equilibrium_density, oxygen_chemical_equilibrium_fractions
      public :: ith_H, ith_H2, ith_O, ith_OH, ith_H2O, ith_CO, ith_C,     &
                n_ox_sp, ox_species_name
      public :: n_fuv_band, fuv_band_lo_A, fuv_band_hi_A, fuv_band_name,  &
                qy_H2O_OH_H, qy_H2O_H2_O1D, qy_H2O_O_H_H, qy_OH_O_H
      public :: photolysis_threshold_erg, e_excite_O1D_erg,               &
                ich_H2O_OH_H, ich_H2O_H2_O, ich_H2O_O_H_H, ich_OH_O_H

      integer, parameter :: dp = kind(1.0d0)

      ! Universal gas constant [J mol^-1 K^-1], CODATA 2018 exact value
      ! R = N_A k_B with N_A = 6.02214076e23 and k_B = 1.380649e-23 J/K.
      real(dp), parameter :: R_gas_J = 8.31446261815324d0
      ! Boltzmann constant [erg K^-1], CODATA 2018 exact
      real(dp), parameter :: kb_erg  = 1.380649d-16
      ! Standard-state pressure, 1 bar [dyn cm^-2]
      real(dp), parameter :: p_std_cgs = 1.0d6

      ! ---------------------------------------------------------------
      ! Species indices for the THERMODYNAMIC TABLE of this module.  The
      ! prefix is ith_ and not isp_ on purpose: species_table already owns
      ! isp_H2 / isp_OH / isp_H2O / isp_CO as f_sp COLUMN numbers (34, 38,
      ! 39, 40), and Fortran is case-insensitive, so a module using both
      ! tables under the same prefix would either fail to compile or, with
      ! a narrow only: list, silently index one table with the other's
      ! numbers.  Two different meanings, two different names.
      integer, parameter :: ith_H = 1, ith_H2 = 2, ith_O = 3,             &
                            ith_OH = 4, ith_H2O = 5, ith_CO = 6,          &
                            ith_C = 7
      integer, parameter :: n_ox_sp = 7
      character(len=3), parameter :: ox_species_name(n_ox_sp) =           &
           (/ 'H  ', 'H2 ', 'O  ', 'OH ', 'H2O', 'CO ', 'C  ' /)

      ! Shomate temperature-range boundaries.  Each species has at most
      ! n_shom_max intervals; n_shom(i) gives how many it has, and
      ! t_shom(0:n_shom(i), i) their edges [K].  The lowest interval of each
      ! species is the 10-298 K fit that Photochem added below the NIST
      ! tables; EXHALE never reaches those temperatures, but it is carried so
      ! that the table is the same object as the reference file.
      integer, parameter :: n_shom_max = 4
      integer, parameter :: n_shom(n_ox_sp) = (/ 1, 4, 2, 3, 3, 3, 2 /)
      real(dp), parameter :: t_shom(0:n_shom_max, n_ox_sp) = reshape(     &
        (/ &
        !   H : NIST gives one interval, 298-6000 K
           0.0d0,  6000.0d0,     0.0d0,    0.0d0,    0.0d0,               &
        !   H2: 10-298 (Photochem fit), then NIST 298-1000, 1000-2500,
        !       2500-6000
          10.0d0,   298.0d0,  1000.0d0, 2500.0d0, 6000.0d0,               &
        !   O : 10-298 (Photochem fit), then 298-6000
          10.0d0,   298.0d0,  6000.0d0,    0.0d0,    0.0d0,               &
        !   OH: 10-298, then NIST 298-1300, 1300-6000
          10.0d0,   298.0d0,  1300.0d0, 6000.0d0,    0.0d0,               &
        !   H2O: 10-298, then NIST 500-1700 (used from 298), 1700-6000
          10.0d0,   298.0d0,  1700.0d0, 6000.0d0,    0.0d0,               &
        !   CO: 10-298, then NIST 298-1300, 1300-6000
          10.0d0,   298.0d0,  1300.0d0, 6000.0d0,    0.0d0,               &
        !   C : 10-298 (Photochem fit), then 298-6000
          10.0d0,   298.0d0,  6000.0d0,    0.0d0,    0.0d0                &
        /), (/ n_shom_max+1, n_ox_sp /) )

      ! Shomate coefficients A,B,C,D,E,F,G (NIST's H is dropped; see header).
      ! Cp  [J/mol/K] = A + B t + C t^2 + D t^3 + E/t^2
      ! H   [kJ/mol]  = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F
      ! S   [J/mol/K] = A ln t + B t + C t^2/2 + D t^3/3 - E/(2 t^2) + G
      ! with t = T/1000.
      real(dp), parameter :: c_shom(7, n_shom_max, n_ox_sp) = reshape(    &
        (/ &
        ! ---- H  (NIST WebBook C12385136, 298-6000 K, Chase 1998) -------
          20.78603d0, 4.85d-10, -1.58d-10, 1.53d-11, 3.2d-11,             &
          211.8d0, 139.87d0,                                              &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
        ! ---- H2 (NIST WebBook C1333740, Chase 1998) --------------------
          28.93235d0, 6.885775d0, -71.92969d0, 160.2033d0, 1.129446d-05,  &
          -8.613061d0, 165.4217d0,                                        &
          33.066178d0, -11.36342d0, 11.432816d0, -2.772874d0,             &
          -0.158558d0, -9.980797d0, 172.708d0,                            &
          18.563083d0, 12.257357d0, -2.859786d0, 0.268238d0, 1.97799d0,   &
          -1.147438d0, 156.2881d0,                                        &
          43.41356d0, -4.293079d0, 1.272428d0, -0.096876d0, -20.53386d0,  &
          -38.51515d0, 162.0814d0,                                        &
        ! ---- O  (no Shomate table on the WebBook; see header note) -----
          21.11945d0, -13.50082d0, 137.6414d0, -280.234d0, -2.248628d-05, &
          242.8134d0, 186.9937d0,                                         &
          21.1861d0, -0.502314d0, 0.168694d0, -0.008962d0, 0.075664d0,    &
          243.1306d0, 187.2591d0,                                         &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
        ! ---- OH (NIST WebBook C3352576, Chase 1998) --------------------
          29.35267d0, -10.20092d0, 103.7819d0, -209.5712d0,               &
          -1.701198d-05, 30.18574d0, 219.509d0,                           &
          32.27768d0, -11.36291d0, 13.60545d0, -3.846486d0, -0.001335d0,  &
          29.75113d0, 225.5783d0,                                         &
          28.74701d0, 4.714489d0, -0.814725d0, 0.054748d0, -2.747829d0,   &
          26.41439d0, 214.1166d0,                                         &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
        ! ---- H2O (NIST WebBook C7732185, Chase 1998) -------------------
          33.25764d0, 0.08576948d0, -2.790822d0, 20.93015d0,              &
          -5.017802d-08, -251.7621d0, 228.9958d0,                         &
          30.092d0, 6.832514d0, 6.793435d0, -2.53448d0, 0.082139d0,       &
          -250.881d0, 223.3967d0,                                         &
          41.96426d0, 8.622053d0, -1.49978d0, 0.098119d0, -11.15764d0,    &
          -272.1797d0, 219.7809d0,                                        &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
        ! ---- CO (NIST WebBook C630080, Chase 1998) ---------------------
          29.14634d0, -1.87642d0, 19.73517d0, -44.9939d0, -3.064321d-06,  &
          -119.2159d0, 233.0144d0,                                        &
          25.56759d0, 6.09613d0, 4.054656d0, -2.671301d0, 0.131021d0,     &
          -118.0089d0, 227.3665d0,                                        &
          35.1507d0, 1.300095d0, -0.205921d0, 0.01355d0, -3.28278d0,      &
          -127.8375d0, 231.712d0,                                         &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
        ! ---- C  (no Shomate table on the WebBook; see header note) -----
          20.78655d0, -0.007461516d0, -0.1354649d0, 1.958757d0,           &
          -3.374423d-08, 710.4722d0, 183.2452d0,                          &
          21.1751d0, -0.81243d0, 0.448537d0, -0.043256d0, -0.013103d0,    &
          710.347d0, 183.8734d0,                                          &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,                      &
          0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0                       &
        /), (/ 7, n_shom_max, n_ox_sp /) )

      ! ---------------------------------------------------------------
      ! FUV bands for the photolysis channels O3-O5 (H2O) and O7 (OH).
      !
      ! The band edges are fixed by the branching-ratio intervals and by the
      ! absorbers, not chosen: the three-body branch H2O + hv -> O + H + H
      ! opens only in the interval that contains Ly-alpha.  912 A is the H
      ! Lyman edge, below which atomic hydrogen absorbs everything; 2304 A is
      ! the long-wavelength end of the H2O photodissociation cross section
      ! (230.413 nm, set by the measurement of Ranjan et al. 2020).
      !
      ! 1110 A IS A BAND EDGE BECAUSE THE H2 LYMAN-WERNER SYSTEM ENDS THERE.
      ! Draine & Bertoldi (1996) take 912-1110 A for the Solomon process
      ! (their footnote 4: essentially all H2 pumping out of v = 0 happens
      ! shortward of 1110 A), and lyman_werner.f90 uses that interval.  Over
      ! it the SAME photons are absorbed by H2 in lines and by H2O and OH in
      ! a continuum, so 912-1110 A must be ONE band with ONE incident flux
      ! and one beam that the three absorbers share; otherwise the energy of
      ! the interval is counted twice, once in "Stellar LW flux" and once
      ! inside a band that reaches across it.  The first band below is
      ! therefore the Lyman-Werner interval itself, carrying the flux of the
      ! existing "Stellar LW flux" key, and B1 begins where it ends.  The
      ! shared beam is built in water_photolysis.f90 section 3.
      !
      ! PROVENANCE OF THE QUANTUM YIELDS.  The table below is the reading of
      ! the photodissociation-qy dataset of photochem/data/xsections/H2O.h5,
      ! whose nodes sit at 105, 120.1, 120.2, 123.0, 123.1, 145.0, 145.1 and
      ! 185 nm and are constantly extrapolated outside them.  The numbers
      ! themselves come from two measurements, by way of the summary in
      ! JPL Publication 19-5 (Burkholder et al. 2020) entry B2, pp. 4-37 to
      ! 4-42:
      !   * Stief, Payne & Klemm (1975), J. Chem. Phys. 62, 4000, give
      !     0.89 / 0.11 over 105-145 nm and >= 0.99 / <= 0.01 over
      !     145-185 nm, band-averaged over a lithium-fluoride flash lamp;
      !     they give no single-wavelength Ly-alpha value.
      !   * Slanger & Black (1982), J. Chem. Phys. 77, 2432, give
      !     0.78 / 0.10 / 0.12 at 121.567 nm.
      !
      ! FOUR CAVEATS, EVERY ONE OF THEM PHYSICS AND NOT BOOKKEEPING.
      !
      ! 1. JPL 19-5 does NOT recommend these yields.  Within entry B2 the
      !    word recommend attaches only to the absorption cross section.
      !    The yields are quoted as literature values, and at Ly-alpha the
      !    entry sets out a SECOND, disagreeing set from Mordaunt et al.
      !    (1994) -- 0.64 / 0.11 / 0.11 with a fourth channel
      !    OH(A 2Sigma+) + H at 0.14 -- without choosing between them.
      !
      ! 2. The 0.78 at Ly-alpha is Phi2 + Phi4, the sum of the ground-state
      !    OH(X 2Pi) + H channel and the electronically excited
      !    OH(A 2Sigma+) + H channel.  Writing it as one OH + H channel
      !    buries the OH(A) production, which Mordaunt et al. split off at
      !    0.14.  A network that wants the OH(A) emission separately cannot
      !    take this table as it stands.
      !
      ! 3. The 1.00 / 0.00 of band B4 is a rounding of Stief's >= 0.99 and
      !    <= 0.01.
      !
      ! 4. VULCAN's H2O_branch.csv is NOT the same table as a function of
      !    wavelength, although its four ratio triplets are the same
      !    numbers.  Its Ly-alpha node sits at 121.0 nm, 0.567 nm short of
      !    the line, and VULCAN interpolates linearly between nodes, so at
      !    the actual line center it uses 0.837 / 0.105 / 0.058 -- the
      !    three-body branch is 0.058 there, not 0.12.  It also has no flat
      !    1231-1450 A interval, ramping instead from 0.89 to 1.00.  A third
      !    source, the PHIDRATES database behind Huebner & Mukherjee (2015),
      !    gives 0.75 / 0.10 / 0.15 flat over 984-1304 A and attributes it
      !    to the same Slanger & Black paper.  The spread on the Ly-alpha
      !    three-body branch across these sources is 0.058 to 0.15.
      !
      ! OH PHOTOLYSIS IS A MERGED CHANNEL, NOT A MEASURED UNIT YIELD.  The
      ! qy_OH_O_H = 1 below is the Leiden default: Heays, Bosman & van
      ! Dishoeck (2017), A&A 602, A105, sec. 3.1 state that they generally
      ! do not divide the photoabsorption cross section into decay channels,
      ! and OH is not one of their exceptions.  Their Table 1 gives OH a
      ! dissociation threshold of 279 nm and an ab initio cross section from
      ! van Dishoeck & Dalgarno; the 264.9 nm that the data file ends at is
      ! a file endpoint, not a statement in the paper.  VULCAN's
      ! OH_branch.csv says outright in its header that it combines the
      ! O(1D) branch below 150 nm with the O(3P) branch above it, and
      ! PHIDRATES puts the O(1D) branch at 0.88 at Ly-alpha.  This matters
      ! for A2 specifically: O(1D) made by OH photolysis feeds straight back
      ! into O6, and a merged unit yield erases that.  M2 must either split
      ! the OH channel or state that it does not.
      !
      ! CROSS SECTIONS, for M2.  The H2O grid of H2O.h5 runs 0.1-230.413 nm
      ! with a photodissociation peak of 2.75515e-17 cm^2 at 111.5 nm; the
      ! OH grid runs 0.06-264.9 nm with a photodissociation peak of
      ! 1.36397e-17 cm^2 at 100 nm.  VULCAN's H2O_cross.csv covers only
      ! 6.2-230.413 nm and peaks 5.7% lower, at 2.598e-17 cm^2.  The
      ! wavelength ranges are concatenations made by Photochem, not
      ! statements of any one paper: Huebner & Mukherjee (2015) below
      ! 6.3 nm, Heays et al. (2017) to 192.056 nm (their own Table 7 runs to
      ! 193.9 nm), Ranjan et al. (2020) from there to 230.413 nm, the last
      ! being a direct 292 K measurement over 186-230 nm.
      integer, parameter :: n_fuv_band = 5
      real(dp), parameter :: fuv_band_lo_A(n_fuv_band) =                  &
           (/  912.0d0, 1110.0d0, 1202.0d0, 1231.0d0, 1451.0d0 /)
      real(dp), parameter :: fuv_band_hi_A(n_fuv_band) =                  &
           (/ 1110.0d0, 1201.0d0, 1230.0d0, 1450.0d0, 2304.0d0 /)
      character(len=2), parameter :: fuv_band_name(n_fuv_band) =          &
           (/ 'LW', 'B1', 'B2', 'B3', 'B4' /)

      ! H2O photodissociation quantum yields on those bands.  The yield
      ! table of H2O.h5 is flat below its 120.1 nm node, so the LW band and
      ! B1 -- the two halves of what was one 912-1201 A interval -- carry
      ! the same triplet.  B2 contains Ly-alpha (1215.67 A) and is the only
      ! band where the three-body branch is open.
      real(dp), parameter :: qy_H2O_OH_H(n_fuv_band)   =                  &
           (/ 0.89d0, 0.89d0, 0.78d0, 0.89d0, 1.00d0 /)
      real(dp), parameter :: qy_H2O_H2_O1D(n_fuv_band) =                  &
           (/ 0.11d0, 0.11d0, 0.10d0, 0.11d0, 0.00d0 /)
      real(dp), parameter :: qy_H2O_O_H_H(n_fuv_band)  =                  &
           (/ 0.00d0, 0.00d0, 0.12d0, 0.00d0, 0.00d0 /)
      ! OH: one merged dissociation channel, unit yield.  See the caveat
      ! above -- this is a Leiden default, not an evaluated branching.
      real(dp), parameter :: qy_OH_O_H(n_fuv_band) =                      &
           (/ 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.00d0 /)

      ! Photolysis channel identifiers, for photolysis_threshold_erg.
      integer, parameter :: ich_H2O_OH_H  = 1, ich_H2O_H2_O = 2,          &
                            ich_H2O_O_H_H = 3, ich_OH_O_H   = 4
      ! Avogadro constant [mol^-1], CODATA 2018 exact
      real(dp), parameter :: N_avog = 6.02214076d23

      ! Excitation energy of O(1D) above the ground O(3P) term [erg].
      ! Taken as the difference of the F Shomate coefficients that
      ! zahnle_earth.yaml carries for O1D and O, 432.6334 - 243.1306 =
      ! 189.5028 kJ/mol = 1.9641 eV.  The file's own note on that entry is
      ! "Estimated from thermodynamic data at 298 K and species O", so this
      ! is a thermal 298 K value referred to the statistically averaged 3P
      ! term, not a spectroscopic 1D2 - 3P2 term difference; the two differ
      ! by the 3P fine-structure population and the difference is at the
      ! 0.2% level.  Only the photolysis energetics of O4 use it, and that
      ! is M2's subject, not M1's.
      real(dp), parameter :: e_excite_O1D_erg = 3.14676d-12

      contains

! =====================================================================
!  Rate coefficients
! =====================================================================

      ! O1: OH + H2 -> H2O + H
      !     k = 3.6e-16 T^1.52 exp(-1740/T) cm^3 molecule^-1 s^-1,
      !     250-2500 K, dlog k = +-0.1 at 250 K rising to +-0.3 at 2500 K.
      !     Baulch et al. (2005), J. Phys. Chem. Ref. Data 34, 757,
      !     Table 4.1 (OH radical reactions) and the data sheet on p. 1029.
      !
      ! The dominant channel of the whole network: 95.3% of the net H2 loss
      ! at the HD 189733 b base and 20.6% at the HD 209458 b base.  Its
      ! reverse, H2O + H -> OH + H2, is 38.1% of the H2O loss there and is
      ! generated by rate_from_detailed_balance, not transcribed.
      !
      ! WHY THE 2005 SUPPLEMENT AND NOT THE 1992 EVALUATION.  The earlier
      ! evaluation, Baulch et al. (1992), JPCRD 21, 411, Table 1 p. 419 and
      ! data sheet p. 552, recommends 1.7e-16 T^1.6 exp(-1660/T) over
      ! 300-2500 K, adopted there from Zellner (1979).  That is the value
      ! zahnle_earth.yaml carries (as 1.740708e-16, 2.4% higher).  VULCAN's
      ! id 1 carries 3.57e-16 T^1.52 exp(-1740/T) from Oldenborg & Loge
      ! (1992), which is the 2005 recommendation to 0.8%.  So the two
      ! reference networks are not two readings of one evaluation; they are
      ! the two successive CEC evaluations of the same reaction, and the
      ! later one is adopted here.  Measured, the 1992 form is 2.7% below
      ! the 2005 one at 300 K, 11% below at 864 K and 8.8% below at 2500 K,
      ! all far inside the +-0.1 to +-0.3 dex the evaluations quote.
      double precision function rk_O1_OH_H2_water(T) result(k)
      real(dp), intent(in) :: T
      k = 3.6d-16*T**1.52d0*exp(-1740.0d0/T)
      end function

      ! O2: O + H2 -> OH + H
      !     k = [6.34e-12 exp(-4000/T) + 1.46e-9 exp(-9650/T)]
      !     cm^3 molecule^-1 s^-1, 298-3300 K, dlog k = +-0.2 over the whole
      !     range.  Baulch et al. (2005), Table 4.1 (O atom reactions) and
      !     the data sheet on p. 804.
      !
      ! 3.7% of the net at the HD 189733 b base, 7.8% at the HD 209458 b
      ! base.  Its reverse, OH + H -> O + H2, is 2.4% of the OH loss.
      !
      ! WHY THE 2005 SUPPLEMENT.  Both reference networks carry the 1992
      ! form, 8.5e-20 T^2.67 exp(-3163/T), 300-2500 K, dlog k = +-0.5 at
      ! 300 K falling to +-0.2 above 500 K (Baulch et al. 1992, Table 1
      ! p. 416, data sheet p. 430) -- and both carry it with exp(-3160/T),
      ! a 3 K slip from the published 3163 that is worth 0.1% at 1000 K.
      ! Measured, the 2005 two-term form agrees with the 1992 one to 10%
      ! at 300 K and to 12% at 2500 K, but is 1.84x SMALLER at 864 K and
      ! 1.75x smaller at 1000 K -- a bigger change than either evaluation's
      ! quoted uncertainty allows for, and the largest single change this
      ! audit makes to a rate that is actually in the minimal set.  The
      ! later evaluation is adopted; the consequence for the measured
      ! budget is in docs/a2_reaction_audit.md.
      double precision function rk_O2_O_H2_hydroxyl(T) result(k)
      real(dp), intent(in) :: T
      k = 6.34d-12*exp(-4000.0d0/T) + 1.46d-9*exp(-9650.0d0/T)
      end function

      ! O6: O(1D) + H2 -> OH + H
      !     k = 1.1e-10 cm^3 molecule^-1 s^-1, independent of temperature
      !     over 200-350 K, dlog k = +-0.1 at 298 K.  Atkinson et al.
      !     (2004), Atmos. Chem. Phys. 4, 1461, data sheet I.A2.18,
      !     pp. 1508-1509 (IUPAC subcommittee evaluation).
      !
      ! The reason O(1D) is carried at all, and the reason its steady state
      ! closes: this is its dominant sink in an H2-rich gas.
      !
      ! THE LARGEST OPEN UNCERTAINTY IN THIS MODULE, on three counts.
      !   * Range.  The evaluation covers 200-350 K and the base is at
      !     560-2400 K.  The value is used outside its stated range because
      !     the reaction is a barrierless insertion whose rate the same
      !     evaluation finds temperature-independent where it was measured;
      !     that is an extrapolation and it is stated as one.
      !   * Competing values.  VULCAN's id 615 gives 2.87e-10 over
      !     100-2100 K, cited to Tully (1975), J. Chem. Phys. 62, 1893 --
      !     2.6x larger, and the NIST Chemical Kinetics Database record for
      !     that paper labels it an RRK(M) extrapolation, not a measurement,
      !     so it disagrees with the four independent room-temperature
      !     measurements the IUPAC evaluation averages.  Photochem's
      !     zahnle_earth.yaml gives 1.5e-10, 36% above the IUPAC value,
      !     with the key Ba92; that attribution does not hold, because
      !     Baulch et al. (1992), Baulch et al. (2005) and Tsang & Hampson
      !     (1986) are combustion evaluations and none of them carries any
      !     O(1D) chemistry.  The Photochem value has no source this audit
      !     could trace.
      !   * Channel.  The 1.1e-10 is the TOTAL k1 + k2, where channel 2 is
      !     physical quenching O(1D) + H2 -> O(3P) + H2.  The evaluation's
      !     comment is that channel 1 is the dominant pathway, above 95%
      !     (Wine & Ravishankara 1982).  Using the total for the reactive
      !     channel therefore overstates it by at most 5%, and it is the
      !     right number for the O(1D) lifetime either way.  This is also
      !     the published answer to O11: quenching by H2 exists, it is under
      !     5% of the collisions, and it has no separately evaluated rate.
      !
      ! The consequence of the 2.6x is measured: with 2.87e-10 this channel
      ! is 5.5% of O1 at the HD 189733 b base and with 1.1e-10 it is 2.1%,
      ! against 5.3% for ground-state O through O2.  The design's statement
      ! that O(1D) beats ground-state O there is a correct report of a run
      ! made with the VULCAN rate set and does not survive the IUPAC value.
      double precision function rk_O6_O1D_H2_hydroxyl() result(k)
      k = 1.1d-10
      end function

      ! O8: O + H + M -> OH + M
      !     k0 = 1.3e-29 T^-1 cm^6 molecule^-2 s^-1, no temperature range
      !     stated, uncertainty a factor of 10.  Tsang & Hampson (1986),
      !     J. Phys. Chem. Ref. Data 15, 1087, entry 5,4, data sheet
      !     p. 1111 and summary table p. 1091.
      !
      ! Returns the low-pressure limit k0; combine with k_inf through
      ! lindemann_rate.  The oxygen analogue of R15 (H + H + M -> H2 + M)
      ! in mol_rates.f90.
      !
      ! WHAT THE SOURCE ACTUALLY SAYS.  The data sheet states that there are
      ! no definitive measurements of this rate, that the numbers in the
      ! literature are rough estimates spanning 1e-33 to 1e-30 cm^6 s^-1,
      ! that the uncertainty is a factor of 10, and that the reaction is not
      ! very important under combustion conditions.  This is an estimate,
      ! not an evaluated measurement, and no temperature range is attached
      ! to it.  Both reference networks carry it (Photochem as 1.29e-29,
      ! 0.8% low; VULCAN's id 657 as 1.30e-29 over 300-2500 K, a range the
      ! source does not give).  Their high-pressure limits differ by a
      ! factor of 10 (1.0e-12 against 1.0e-11) and neither is in Tsang &
      ! Hampson: Photochem's 1.0e-12 is the same placeholder it puts on
      ! H + H (+M) and O + O (+M).  The VULCAN value is used below where a
      ! k_inf is needed, and it is numerically irrelevant -- see the THIRD
      ! BODY M note in the header.
      !
      ! NOT IN THE MINIMAL SET: measured at 8.5e-9 of O1 at the
      ! HD 189733 b base.
      double precision function rk_O8_O_H_assoc(T) result(k0)
      real(dp), intent(in) :: T
      k0 = 1.3d-29*T**(-1.0d0)
      end function

      ! O9: H + OH + M -> H2O + M
      !     k0(N2)  = 6.1e-26 T^-2.0 cm^6 molecule^-2 s^-1, 300-3000 K,
      !               dlog k = +-0.5
      !     k0(Ar)  = 2.3e-26 T^-2.0, 300-3000 K, dlog k = +-0.3
      !     k0(H2O) = 3.9e-25 T^-2.0, 300-3000 K, dlog k = +-0.5
      !     Baulch et al. (1992), JPCRD 21, 411, Table 3 p. 428 and data
      !     sheets pp. 496-498; carried unchanged into Baulch et al. (2005),
      !     Table 4.2 and data sheet p. 913.
      !
      ! The function returns the N2 value.  N2 is the diatomic collider of
      ! the three and is the closest published analogue to an H2 bath; Ar is
      ! the monatomic analogue for He, a factor 2.7 below it; H2O is 6.4x
      ! above N2 and is the wrong collider for a hydrogen-helium atmosphere.
      !
      ! THIS RESOLVES THE APPARENT TWO-DECADE DISAGREEMENT between the
      ! reference networks, which is not a disagreement about the reaction.
      ! VULCAN's id 659, k0 = 3.89e-25 T^-2, IS Baulch's H2O-collider value.
      ! Photochem's 1.050748e-26 T^-2.1 is Javoy et al. (2003), Exp. Therm.
      ! Fluid Sci. 27, 371, section 3.4, read here.  It is their ARGON
      ! value, k = 3.75e21 T^-2.1 cm^6 mol^-2 s^-1 = 1.034e-26 cm^6
      ! molecule^-2 s^-1 after division by N_A^2 = 3.6266e47, over
      ! 2790-3200 K at about 250 kPa total in 1200-4500 ppm H2O diluted in
      ! Ar, quoted at +-25%.  It is not a direct measurement of the
      ! association: they measured the H2O dissociation rate behind
      ! reflected shock waves and inverted it through the equilibrium
      ! constant.  The coefficient Photochem carries is 1.6% above that
      ! conversion, far inside the paper's own +-25%, so the NIST record
      ! 2003JAV/NAU371-377:10 is confirmed; what the record does not say is
      ! that the collider is Ar.
      !
      ! The two networks therefore quote DIFFERENT COLLIDERS, VULCAN H2O
      ! and Photochem Ar, on top of a source disagreement.  Like for like
      ! against Baulch's own Ar value, Javoy is 3.9x low at 300 K and 5.0x
      ! low at 3000 K; against the N2 value adopted here it is 10.4x low at
      ! 300 K and 13.1x low at 3000 K.  The 66-81x between the networks
      ! over 300-2500 K is about 17x of collider (Baulch's own H2O/Ar ratio
      ! is 17.0, Javoy's assumed efficiency 18) times 3.9-4.9x of source.
      ! Neither network is using a value appropriate to an H2/He bath, and
      ! Javoy's validity range begins well above the molecular layer, so
      ! the Baulch N2 value stands.
      !
      ! Photochem's high-pressure limit, 2.7e-10 exp(-75/T), is Cobos &
      ! Troe (1985), J. Chem. Phys. 83, 1010, Table I entry (14)
      ! H + OH -> H2O, read here.  The paper prints no Arrhenius form: it
      ! tabulates k_rec,inf = 2.1e-10 at 300 K and 2.6e-10 at 2100 K in
      ! cm^3 molecule^-1 s^-1, and 2.7e-10 exp(-75/T) is the two-point fit
      ! through them (2.10e-10, 2.61e-10) -- which is where the 300-2100 K
      ! of the NIST record 1985COB/TRO1010-1015:15 comes from.  The method
      ! is the simplified statistical adiabatic channel model of that
      ! paper's Part I with the looseness parameter fitted (alpha = 1.0
      ! A^-1, beta = 2.1 A^-1, alpha/beta = 0.48), not transition-state
      ! theory as the record labels it, and its 300 K experimental anchor
      ! is the isotope exchange OH + D -> OD + H, not a direct H + OH
      ! association measurement.  Baulch recommends no k_inf for this
      ! reaction.
      !
      ! NOT IN THE MINIMAL SET: measured at 1.1e-7 of O1 at the HD 189733 b
      ! base with the N2 value adopted here.
      double precision function rk_O9_H_OH_assoc(T) result(k0)
      real(dp), intent(in) :: T
      k0 = 6.1d-26*T**(-2.0d0)
      end function

      ! O10: OH + OH -> H2O + O
      !     k = 5.56e-20 T^2.42 exp(+970/T) cm^3 molecule^-1 s^-1,
      !     250-2400 K, dlog k = +-0.15.  Baulch et al. (2005), Table 4.1
      !     and data sheet p. 1032.
      !
      ! NOT IN THE MINIMAL SET: measured at 1.5e-6 of O1 at the HD 189733 b
      ! base.  It is transcribed because it is the one reaction in this
      ! module whose reverse is written out explicitly in the other network
      ! (VULCAN id 5, O + H2O -> OH + OH, 8.20e-14 T^0.95 exp(-8570/T),
      ! 250-2400 K), so the pair is an end-to-end test of
      ! rate_from_detailed_balance against two independent transcriptions of
      ! opposite directions of the same reaction.
      !
      ! WHY THE 2005 SUPPLEMENT, AND A MIS-ATTRIBUTION FIXED.  The 1992
      ! evaluation recommends 2.5e-15 T^1.14 exp(-50/T) over 250-2500 K,
      ! dlog k = +-0.2 (Table 1 p. 419, data sheet p. 555), taken there from
      ! Ernst, Wagner & Zellner; that is what zahnle_earth.yaml carries, as
      ! 2.549944e-15.  The yaml labels it "Ba92, Li91", but Lifshitz &
      ! Michael (1991) is a study of the REVERSE direction and it is not in
      ! the 1992 data sheet's reference list at all; it is cited by the 2005
      ! evaluation, which is the one that revised the recommendation.
      ! Measured, the two forms agree to 1.3% at 300 K but the 1992 one is
      ! 2.4x larger at 864 K and 1.3x larger at 2500 K.
      double precision function rk_O10_OH_OH_water(T) result(k)
      real(dp), intent(in) :: T
      k = 5.56d-20*T**2.42d0*exp(970.0d0/T)
      end function

      ! O12: O(1D) + H2O -> OH + OH
      !     k = 2.2e-10 cm^3 molecule^-1 s^-1, independent of temperature
      !     over 200-350 K, dlog k = +-0.1 at 298 K.  Atkinson et al.
      !     (2004), Atmos. Chem. Phys. 4, 1461, data sheet I.A2.19, p. 1510
      !     (total of the three channels; the 2 HO channel dominates).
      !
      ! NOT IN THE MINIMAL SET: the second O(1D) sink after O6, measured at
      ! 1.2e-5 of O1 at the HD 189733 b base and giving O(1D) a lifetime of
      ! 1.1 s there against 1.6e-3 s from O6 -- 700x slower, so it changes
      ! the O(1D) steady state by 0.14%.  It is transcribed so that the
      ! O(1D) closure records a value for its only competitor rather than an
      ! omission.  This is the value zahnle_earth.yaml carries; VULCAN's
      ! id 619 gives 1.62e-10 exp(+65/T) over 235-370 K, 21% below the
      ! evaluation at 300 K.
      double precision function rk_O12_O1D_H2O_hydroxyl() result(k)
      k = 2.2d-10
      end function

      ! Lindemann combination of a low-pressure limit k0 [cm^6 s^-1] and a
      ! high-pressure limit k_inf [cm^3 s^-1] at third-body density n_M
      ! [cm^-3], returning the two-body-equivalent rate [cm^3 s^-1]:
      !     k = k_inf Pr/(1+Pr),   Pr = k0 n_M / k_inf .
      ! This is what both reference codes reduce to when no Troe parameters
      ! are supplied, which is the case for O8 and O9.
      double precision function lindemann_rate(k0, kinf, n_M) result(k)
      real(dp), intent(in) :: k0, kinf, n_M
      real(dp) :: Pr
      Pr = k0*n_M/kinf
      k  = kinf*Pr/(1.0d0 + Pr)
      end function
! =====================================================================
!  Thermodynamics and detailed balance
! =====================================================================

      ! Temperature at which the Shomate polynomials are evaluated: T
      ! clamped to the tabulated range of the species.
      !
      ! WHY A CLAMP AND NOT AN EXTRAPOLATION.  A Shomate fit is a quartic in
      ! t = T/1000, so evaluating it outside its range does not degrade
      ! gracefully -- it diverges.  At 1e6 K, a temperature an unconverged
      ! cold-start wind passes through, the t^4 term reaches 1e9 kJ/mol and
      ! the equilibrium constant built from it is meaningless; measured, the
      ! unclamped table kept CO fully associated at 9e5 K, which is the
      ! opposite of the physics.  Holding the polynomial at its endpoint is
      ! not thermodynamically exact either, but it is the correct LIMIT for
      ! every use this module has: the table tops out at 6000 K, and every
      ! molecule in this set (CO, H2O, OH, H2) is already dissociated by
      ! many decades there, so the clamped constant says "fully dissociated"
      ! and goes on saying it.  A run that needs real thermochemistry above
      ! 6000 K needs a different table, not this one extrapolated.
      double precision function t_shomate(isp, T) result(Tc)
      integer,  intent(in) :: isp
      real(dp), intent(in) :: T
      Tc = min(max(T, t_shom(0, isp)), t_shom(n_shom(isp), isp))
      end function

      ! Index of the Shomate interval that contains T, clamped to the end
      ! intervals outside the tabulated range.
      integer function shomate_interval(isp, T) result(k)
      integer,  intent(in) :: isp
      real(dp), intent(in) :: T
      integer :: i
      k = n_shom(isp)
      do i = 1, n_shom(isp)
        if (T < t_shom(i, isp)) then
          k = i
          exit
        end if
      end do
      if (T < t_shom(0, isp)) k = 1
      end function

      ! Absolute molar enthalpy [kJ mol^-1], including the enthalpy of
      ! formation (NIST's H coefficient is not subtracted; see header).
      double precision function enthalpy_shomate(isp, T) result(H)
      integer,  intent(in) :: isp
      real(dp), intent(in) :: T
      real(dp) :: c(7), tred
      integer  :: k
      ! tred = T/1000, the reduced temperature the Shomate form uses.  It is
      ! deliberately not named t: Fortran is case-insensitive and t would
      ! shadow the dummy T.
      k = shomate_interval(isp, T)
      c = c_shom(:, k, isp)
      tred = t_shomate(isp, T)/1000.0d0
      H = c(1)*tred + c(2)*tred**2/2.0d0 + c(3)*tred**3/3.0d0            &
        + c(4)*tred**4/4.0d0 - c(5)/tred + c(6)
      end function

      ! Standard molar entropy at 1 bar [J mol^-1 K^-1].
      double precision function entropy_shomate(isp, T) result(S)
      integer,  intent(in) :: isp
      real(dp), intent(in) :: T
      real(dp) :: c(7), tred
      integer  :: k
      k = shomate_interval(isp, T)
      c = c_shom(:, k, isp)
      tred = t_shomate(isp, T)/1000.0d0
      S = c(1)*log(tred) + c(2)*tred + c(3)*tred**2/2.0d0                &
        + c(4)*tred**3/3.0d0 - c(5)/(2.0d0*tred**2) + c(7)
      end function

      ! Standard molar Gibbs energy at 1 bar [J mol^-1].
      double precision function gibbs_energy_shomate(isp, T) result(G)
      integer,  intent(in) :: isp
      real(dp), intent(in) :: T
      ! The T that multiplies the entropy is the clamped one too: G has to
      ! be the Gibbs energy of one state, not a mixture of an endpoint
      ! enthalpy with an out-of-range temperature.
      G = 1000.0d0*enthalpy_shomate(isp, T)                               &
        - t_shomate(isp, T)*entropy_shomate(isp, T)
      end function

      ! Equilibrium constant in number-density units for
      !     sum(reactants) -> sum(products) ,
      ! i.e. K_c = prod(n_products)/prod(n_reactants) at equilibrium,
      ! in cm^(3*dn) with dn = n_prod - n_react:
      !     K_c = exp(-dG/RT) * (P0/(kB T))^dn ,   P0 = 1 bar.
      ! ireact and iprod hold species indices (ith_* above); repeat an index
      ! to give it a stoichiometric coefficient of 2.
      double precision function equilibrium_constant_conc(ireact, iprod, T)&
                                result(Kc)
      integer,  intent(in) :: ireact(:), iprod(:)
      real(dp), intent(in) :: T
      real(dp) :: dG, n_std, Teq
      integer  :: i, dn
      dG = 0.0d0
      do i = 1, size(iprod)
        dG = dG + gibbs_energy_shomate(iprod(i), T)
      end do
      do i = 1, size(ireact)
        dG = dG - gibbs_energy_shomate(ireact(i), T)
      end do
      dn    = size(iprod) - size(ireact)
      ! Evaluated at the same clamped temperature the Gibbs energies above
      ! were, so that K_c is the equilibrium constant of one state (see
      ! t_shomate). Every species of this table shares the same 6000 K
      ! ceiling, so one clamp covers the whole reaction.
      Teq   = t_shomate(ith_H2, T)
      n_std = p_std_cgs/(kb_erg*Teq)
      Kc    = exp(-dG/(R_gas_J*Teq))*n_std**dn
      end function

      ! Threshold energy of a photolysis channel [erg per event], as the
      ! reaction enthalpy at 298.15 K taken from the Shomate table above.
      ! A photon of energy hv deposits hv - threshold as fragment kinetic
      ! energy, the way e_lw_fragment_erg does in lyman_werner.f90.
      !
      ! APPROXIMATION.  The exact threshold is the 0 K dissociation energy
      ! D0, and dH(298.15) exceeds it by the 0-298 K enthalpy content of the
      ! fragments minus that of the parent -- of order a few kJ/mol against
      ! reaction enthalpies of 428-927 kJ/mol, i.e. under 1%.  The Shomate
      ! fits do not reach 0 K, so dH(298.15) is used and the approximation is
      ! stated here rather than hidden.  Channel ich_H2O_H2_O returns the
      ! ground-state O(3P) threshold; the O(1D) channel O4 of the design adds
      ! e_excite_O1D_erg to it.
      double precision function photolysis_threshold_erg(ichan) result(e)
      integer, intent(in) :: ichan
      real(dp), parameter :: T_ref = 298.15d0
      real(dp) :: dH_kJmol
      select case (ichan)
      case (ich_H2O_OH_H)     ! H2O + hv -> OH + H
        dH_kJmol = enthalpy_shomate(ith_OH, T_ref)                        &
                 + enthalpy_shomate(ith_H,  T_ref)                        &
                 - enthalpy_shomate(ith_H2O, T_ref)
      case (ich_H2O_H2_O)     ! H2O + hv -> H2 + O(3P)
        dH_kJmol = enthalpy_shomate(ith_H2, T_ref)                        &
                 + enthalpy_shomate(ith_O,  T_ref)                        &
                 - enthalpy_shomate(ith_H2O, T_ref)
      case (ich_H2O_O_H_H)    ! H2O + hv -> O + H + H
        dH_kJmol = enthalpy_shomate(ith_O, T_ref)                         &
                 + 2.0d0*enthalpy_shomate(ith_H, T_ref)                   &
                 - enthalpy_shomate(ith_H2O, T_ref)
      case (ich_OH_O_H)       ! OH + hv -> O + H
        dH_kJmol = enthalpy_shomate(ith_O, T_ref)                         &
                 + enthalpy_shomate(ith_H, T_ref)                         &
                 - enthalpy_shomate(ith_OH, T_ref)
      case default
        dH_kJmol = 0.0d0
      end select
      e = dH_kJmol*1.0d3/N_avog*1.0d7      ! kJ/mol -> erg per event
      end function

      ! Carbon monoxide density [cm^-3] of a gas holding n_C_tot carbon and
      ! n_O_tot oxygen nuclei per cm^3 at temperature T, from the
      !
      !     CO  <->  C + O
      !
      ! chemical equilibrium of this module's own thermodynamic table.  With
      ! K_c = n_C n_O / n_CO from equilibrium_constant_conc, conservation of
      ! the two elements gives
      !
      !     (n_C_tot - n_CO)(n_O_tot - n_CO) = K_c n_CO ,
      !
      ! a quadratic whose physical root is the smaller one; it is evaluated
      ! in the form that does not cancel when K_c is small (n_CO -> the
      ! limiting element) or large (n_CO -> n_C n_O / K_c).
      !
      ! WHY THIS IS THE CARBON RESERVOIR OF THE A2 NETWORK.  Decision D4 of
      ! docs/a2_oxygen_option_design.md carries CO as an unreactive oxygen
      ! reservoir -- it holds 45-46% of the oxygen at every level of every
      ! measured arm, so a network that gives the whole oxygen abundance to
      ! the water family over-supplies the OH cycle by about a factor two --
      ! and takes its abundance from a chemical-equilibrium C/O partition
      ! (that document's section 3.7).  The C=O bond is 11.1 eV, four times
      ! the H2O OH bond, so CO survives to temperatures at which the water
      ! family is long gone: with this table the partition turns over near
      ! 4000-5000 K, above the H2 -> H front and below the wind, which is
      ! where the reservoir has to disappear if the carbon and the oxygen of
      ! the wind are to be atomic.
      !
      ! THREE APPROXIMATIONS, EACH STATED RATHER THAN HIDDEN.
      !  1. The equilibrium is written against the FREE ATOMIC C and O of
      !     the gas, which is the only alternative carbon carrier the code
      !     has (the metal block solves C I / C II / C III).  Competition
      !     for the oxygen from H2O and OH is not in the quadratic: the
      !     arguments are the element totals.  Where CO is associated
      !     (T <~ 3000 K) K_c is many decades below either density and the
      !     root is min(n_C, n_O) whatever the other carriers hold, and
      !     where it dissociates the water family is already empty, so the
      !     neglected competition moves the partition only inside the
      !     turnover.
      !  2. Equilibrium is not kinetics.  CO in an irradiated upper
      !     atmosphere is quenched -- it survives above the level at which
      !     equilibrium would destroy it -- so this is a lower bound on the
      !     reservoir there.  Carrying the kinetics would need CO
      !     photodissociation, which predissociates in lines and
      !     self-shields (van Dishoeck & Black 1988), i.e. a second
      !     shielding function and a second band; it is not in the audited
      !     set of docs/a2_reaction_audit.md and is not attempted here.
      !  3. Ionization is not in the balance.  Where CO is associated the
      !     gas is shielded and both elements are neutral, so this costs
      !     nothing there; above the turnover n_CO is negligible.
      double precision function co_equilibrium_density(n_C_tot, n_O_tot, T)&
                                result(n_CO)
      real(dp), intent(in) :: n_C_tot, n_O_tot, T
      real(dp) :: Kc, ssum, prod, disc
      prod = n_C_tot*n_O_tot
      if (prod .le. 0.0d0) then
        n_CO = 0.0d0
        return
      end if
      Kc   = equilibrium_constant_conc((/ ith_CO /), (/ ith_C, ith_O /), T)
      ssum = n_C_tot + n_O_tot + Kc
      disc = ssum*ssum - 4.0d0*prod
      if (disc .lt. 0.0d0) disc = 0.0d0
      n_CO = 2.0d0*prod/(ssum + sqrt(disc))
      ! The physical root cannot exceed the limiting element, and in the
      ! associated limit it IS that element, so the expression above returns
      ! it to within round-off -- occasionally a few parts in 1e16 ABOVE it.
      ! The caller subtracts n_CO from both element totals, so that
      ! round-off would leave one of them negative and put a negative
      ! density into the cooling and the electron sum. Capped here, where
      ! the bound is a statement about the reaction rather than a guard
      ! bolted onto each consumer.
      n_CO = min(n_CO, min(n_C_tot, n_O_tot))
      end function

      ! Chemical-equilibrium partition of the free oxygen family among atomic
      ! O, OH and H2O at the local (T, n_H2, n_HI), returned as fractions of
      ! the family. The two equilibria are the audited reactions themselves,
      !
      !   O  + H2 <-> OH  + H   ->  n_OH /n_O  = K_c(O2) n_H2/n_H
      !   OH + H2 <-> H2O + H   ->  n_H2O/n_OH = K_c(O1) n_H2/n_H
      !
      ! both dimensionless (equal particle counts on the two sides), so the
      ! partition depends on the gas only through the ratio n_H2/n_H.
      !
      ! Two uses, and they are different in kind. (a) It is the molecular-basin
      ! starting point of the cell solve and the seed of a restart that carries
      ! no oxygen columns -- zero is a valid root of the water cycle and hybrd1
      ! is known to be bistable from a zero molecular seed, so the seed must
      ! not be zero. (b) It is what the SOLVED partition must reduce to when
      ! the photolysis rates go to zero, which is gate G3 of the design and the
      ! test that makes the thermodynamically reversed pairs trustworthy.
      !
      ! Evaluated through logarithms because K_c(O1) K_c(O2) (n_H2/n_H)^2
      ! overflows a double at the cool, molecular end.
      subroutine oxygen_chemical_equilibrium_fractions(T, n_h2, n_hi,       &
                                                       f_oh, f_h2o)
      real*8, intent(in)  :: T, n_h2, n_hi
      real*8, intent(out) :: f_oh, f_h2o
      real*8 :: kc1, kc2, ratio, lw(3), lwmax, w(3), wsum
      integer :: i

      f_oh  = 0.0d0
      f_h2o = 0.0d0
      if (n_h2 .le. 0.0d0 .or. n_hi .le. 0.0d0) return

      kc2 = equilibrium_constant_conc((/ ith_O, ith_H2 /),                  &
                                     (/ ith_OH, ith_H /), T)
      kc1 = equilibrium_constant_conc((/ ith_OH, ith_H2 /),                 &
                                     (/ ith_H2O, ith_H /), T)
      if (kc1 .le. 0.0d0 .or. kc2 .le. 0.0d0) return
      ratio = n_h2/n_hi

      lw(1) = 0.0d0                                        ! atomic O
      lw(2) = log(kc2) + log(ratio)                        ! OH
      lw(3) = lw(2) + log(kc1) + log(ratio)                ! H2O
      lwmax = max(lw(1), max(lw(2), lw(3)))
      wsum  = 0.0d0
      do i = 1,3
      	w(i) = exp(max(lw(i) - lwmax, -700.0d0))
      	wsum = wsum + w(i)
      enddo
      f_oh  = w(2)/wsum
      f_h2o = w(3)/wsum

      end subroutine oxygen_chemical_equilibrium_fractions

      ! Reverse rate coefficient of a reaction whose forward coefficient is
      ! k_fwd, by detailed balance: k_rev = k_fwd/K_c.  Units follow from
      ! the stoichiometry, so a two-body forward gives a two-body reverse and
      ! the low-pressure limit of a three-body forward gives the
      ! two-body-equivalent unimolecular reverse.
      double precision function rate_from_detailed_balance(k_fwd, ireact,  &
                                iprod, T) result(k_rev)
      real(dp), intent(in) :: k_fwd
      integer,  intent(in) :: ireact(:), iprod(:)
      real(dp), intent(in) :: T
      k_rev = k_fwd/equilibrium_constant_conc(ireact, iprod, T)
      end function

      end module oxygen_rates

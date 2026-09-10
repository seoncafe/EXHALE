      module molecular_reaction_heat
      ! Net chemical heat of the H2/He network: the energy the COLLISIONAL
      ! reactions of the molecular chemistry release into, or take out of, the
      ! gas.
      !
      ! WHY IT IS NEEDED, AND WHY IT IS ON BY DEFAULT.  In an atomic gas the
      ! ionization energy a photon spends is not a heating channel: the
      ! photoelectron carries away hv - I(H), which PH_heat_HHe deposits, and
      ! the I(H) itself comes back as a Lyman photon when H+ recombines
      ! radiatively and leaves the gas.  In a MOLECULAR gas that return path
      ! is different.  Follow one closed cycle:
      !
      !     H2 + hv  -> H2+ + e       photon pays I(H2) = 15.43 eV; the
      !                               photoelectron keeps hv - I(H2) and
      !                               PH_heat_HHe already deposits that
      !     H2+ + H2 -> H3+ + H       + 1.70 eV
      !     H3+ + e  -> H2 + H        + 9.25 eV   (DISSOCIATIVE: the energy
      !                               goes into the fragments, not a photon)
      !
      ! The cycle returns to H2 having turned one H2 into H + H, and it has
      ! left I(H2) - D0(H2) = 15.43 - 4.48 = 10.95 eV in the gas.  There is no
      ! photon anywhere in it to carry that away.  A code that deposits only
      ! the photoelectron therefore loses 10.95 eV per cycle, and on the
      ! converged He/H = 0.0793 hot Uranus that loss is 81 per cent of the
      ! total heating rate at 1.02 r_base.  Omitting it is not an
      ! approximation of the molecular layer's energy budget; it is most of
      ! it.  Hence the default is ON, and `Molecular reaction heat: False`
      ! exists to reproduce the state the code was in before it.
      !
      ! An atomic run cannot reach this module -- the caller gates on
      ! thereis_mol -- so every non-molecular result is unchanged to the bit
      ! either way.
      !
      ! HOW LARGE.  Measured on the converged He/H = 0.0793 hot Uranus
      ! (docs/p53_h2_adiabatic_index.md and section 141 of
      ! docs/Update_EXHALE_stage1.md), the sum below is 2.4 times the whole
      ! photoelectric heating rate over 1.00-1.05 r_base, falling through 1 at
      ! 1.053 and to 1 percent by 1.10.  The argument above makes that ratio
      ! I(H2)/(hv - I(H2)), so 2.4 corresponds to a mean absorbed photon of
      ! 22 eV: a consistency check on the sum, not an independent derivation
      ! of it.
      !
      ! HOW IT IS BUILT.  Not as a list of reaction enthalpies -- a list can
      ! disagree with itself -- but from ONE table of species formation
      ! energies, species_formation_energy below, which is the eps_s of the
      ! energy ledger (docs/b1_target_system_20260906.md T1.2) and covers
      ! every species the code carries, not only the molecular ones.
      ! eps(X) is the energy needed to build X out of neutral ground-state
      ! free atoms at rest and free electrons at rest, so
      !
      !     heat released by  A + B -> C + D
      !         =  eps(A) + eps(B) - eps(C) - eps(D) ,
      !
      ! and every reaction of the network, in either direction, is consistent
      ! with every other by construction.  A cycle that returns to its
      ! starting species releases exactly zero, which is the property a
      ! reaction-by-reaction list cannot be trusted to have.  The
      ! stoichiometry of each reaction is a table too (mreac_react /
      ! mreac_prod), so that the heat is formed from the reaction rather
      ! than transcribed beside it, and so that a test can check the
      ! nucleus and charge balance of every channel.
      !
      ! WHAT IS DELIBERATELY NOT IN THE SUM, so that nothing is counted twice:
      !
      !   * PHOTON-DRIVEN reactions -- H2 photoionization (P_H2), its
      !     dissociative branch (P_H2_di), Lyman-Werner photodissociation
      !     (k_LW) and the H/He photoionizations.  Their enthalpy is paid by
      !     the absorbed photon, and the share of it that stays in the gas is
      !     already deposited by PH_heat_HHe and by the Lyman-Werner fragment
      !     and fluorescence terms in ionization_equilibrium.
      !   * RADIATIVE recombination of H+ and He+ (R1, R2): the energy leaves
      !     as a photon, and the electron's own kinetic energy is already the
      !     `reco` term of eval_cool.
      !   * COLLISIONAL ionization of H and He (R3, R4): already the `coio`
      !     sink of eval_cool.
      !   * PENNING ionization of H and of H2 by He(2^3S): already deposited
      !     in ionization_equilibrium.  The ASSOCIATIVE branch of the same
      !     collision, He(2^3S) + H -> HeH+ + e, is a reaction of the same
      !     table and its heat is the difference of the reservoir contents,
      !     oxygen_reaction_energy_eV(ir_assoc_HeTR); it is deposited at the
      !     Penning site of ionization_equilibrium, which owns the collision
      !     rate Q31 and the branching f_penning_HeI23S.  It is in the
      !     oreac table below rather than in this sum because that is where
      !     its rate lives, not because it is a different kind of energy.
      !   * H <-> He charge exchange (R21, R22): handled by charge_exchange.
      !
      ! WHERE THE ENERGY IS DEPOSITED.  A reaction that leaves H2 in a
      ! vibrationally excited level can radiate that part away instead of
      ! thermalizing it.  The three-body association H + H + M -> H2 + M is
      ! the channel where the whole 4.48 eV goes into one new molecule, so it
      ! carries Hollenbach & McKee's (1979, ApJS 41, 555, p. 586) collisional
      ! branching (1 + n_cr/n)^-1, evaluated by the module that owns n_cr
      ! (h2_vibrational_relaxation).  In this layer that factor is 1 to a part
      ! in 1e3 -- n(H2) = 4e9 to 5e13 cm^-3 against n_cr ~ 1e6 -- so it is
      ! carried for correctness at the top of a thinner layer, not because it
      ! changes the number here.  The other channels deposit in full.

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,      &
                               isp_HeIII, isp_HeTR, isp_H2, isp_H2p,    &
                               isp_H3p, isp_HeHp, isp_OH, isp_H2O,      &
                               isp_CO, n_mion, mion_elem, mion_ethr,    &
                               melem_i0, mion_fsp, im_OI, im_CII
      use oxygen_rates, only: oxygen_formation_energy_eV, e_excite_O1D_erg, &
                              ith_OH, ith_H2O, ith_CO, ith_H, ith_H2,     &
                              ith_O,                                      &
                              rk_O1_OH_H2_water, rk_O2_O_H2_hydroxyl,     &
                              rate_from_detailed_balance
      use mol_rates, only: rk_R5_H2p_dr,   rk_R6_H3p_dr_H2,             &
                           rk_R7_H3p_dr_3H, rk_R8_H2p_H2,               &
                           rk_R9_H2p_H,    rk_R10_Hp_H2v4,              &
                           rk_R11_H3p_H,   rk_R12_H2_thdis,             &
                           rk_R13_Hp_H2_M, rk_R14_H2_edis,              &
                           rk_R15_3body_H2, rk_R16_HeHp_dr,             &
                           rk_R17_Hep_H2_diss, rk_R18_HeHp_H2,          &
                           rk_R19_HeHp_H,  rk_R20_Hep_H2_HeHp,          &
                           rk_R23_H2_Hep_cx, h2_dissociation_energy_eV
      use h2_vibrational_relaxation, only: h2_vibrational_heat_fraction

      implicit none
      private
      public :: molecular_chemical_heating, formation_energy_report
      public :: species_formation_energy, molecular_reaction_energy_eV
      public :: isp_eps_electron, isp_eps_H_n2, isp_eps_O1D
      public :: n_mreac, mreac_react, mreac_prod, mreac_name
      public :: n_oreac, oreac_react, oreac_prod, oreac_name
      public :: oxygen_reaction_energy_eV, oxygen_chemical_heating
      public :: ir_O1, ir_O1r, ir_O2, ir_O2r, ir_O6, ir_assoc_HeTR,   &
                ir_D1
      public :: formation_energy_density

      ! ------------------------------------------------------------------ !
      ! THE SPECIES FORMATION-ENERGY TABLE eps_s [eV].
      !
      ! REFERENCE STATE, one for the whole code
      ! (docs/b1_target_system_20260906.md T1.2): every element a neutral,
      ! ground-state, free atom at rest, at the zero of its internal level
      ! ladder; a free electron carries zero, so the ionization energy of an
      ! ion sits on the ion.  A metal stage therefore carries the CUMULATIVE
      ! ionization potentials from its neutral ground state, never a zero
      ! placed at whichever stage happens to dominate the base.
      !
      ! NON-OVERLAP WITH THE EQUATION OF STATE (T1.3).  Each internal ladder
      ! has exactly one owner.  eps(H2) = -D0 is measured from v = 0, J = 0,
      ! which is the same zero the rovibrational energy u_rv of caloric_eos
      ! sums from, so the bound ladder is the equation of state's and the
      ! bond is the reservoir's, and neither holds the other's energy.  The
      ! oxygen and carbon molecules follow the same rule through the 0 K
      ! formation energies of oxygen_rates.
      !
      ! PROVENANCE, one line per constant, and what was actually checked.
      !
      ! The ionization energies are the e_th_* of global_parameters, which
      ! carry the full NIST values (13.598434599, 24.587389, 54.417765 and
      ! 15.425927 eV) and not a rounded threshold: one definition of each
      ! potential in the code, and the same energy the photoionization
      ! channels are charged, so a photoionization followed by a
      ! recombination closes exactly.
      !
      !   IP(H), IP(He)   NIST Atomic Spectra Database, ionization energies.
      !                   READ FROM THE DATABASE: 13.598434599702(12) eV and
      !                   24.587389011(25) eV.
      !   IP(H2)          NIST Chemistry WebBook, evaluated ionization energy
      !                   of H2.  READ FROM THE DATABASE: 15.42593(5) eV.
      !   D0(H2)          not repeated here: read through
      !                   mol_rates::h2_dissociation_energy_eV, which is the
      !                   code's one definition of the bond energy
      !                   (D0_H2_cm = 36118.11 cm^-1, Huber & Herzberg 1979).
      !   D0(H3+)         Mizus, Polyansky, McKemmish & Tennyson (2019),
      !                   Mol. Phys. 117, 1663, READ IN THE PUBLISHED TEXT
      !                   (references/Mizus_2019MP_117_1663.pdf).  Their
      !                   abstract: "An improved dissociation energy for H3+
      !                   is derived as D0 = 35, 076 +/- 2 cm-1"; the same
      !                   value with its energy unit in the section that
      !                   follows the potential construction, "and
      !                   D0 = 35, 076 +/- 2 cm-1 = 4.3489 +/- 0.0002 eV".
      !                   The channel is theirs too: "The new BOPES75K can be
      !                   used to give the dissociation limit into
      !                   H3+ -> H2 + H+", which is the reaction the enthalpy
      !                   below is written for.  Their own comparison: this
      !                   is above the earlier theoretical 4.337 +/- 0.002 eV
      !                   of Lie and Frye and below the experimental
      !                   4.381 +/- 0.021 eV of Cosby and Helm.  NOTE this is
      !                   a 0 K dissociation energy and NOT the proton
      !                   affinity of H2, which NIST gives as 422.3 kJ/mol =
      !                   4.377 eV at 298 K -- the wrong quantity for a 0 K
      !                   enthalpy ledger.
      !   D0(HeH+)        DERIVED FROM Coxon & Hajigeorgiou (1999),
      !                   J. Mol. Spectrosc. 193, 306, read in the publisher
      !                   PDF (references/Coxon_1999JMS_193_306.pdf):
      !                       Table 9  De(4HeH+)  = 16448.84   cm^-1
      !                       Table 8  G_0(4HeH+) =  1566.6764 cm^-1
      !                       D0       = De - G_0 = 14882.1636 cm^-1
      !                   Their Table 9 is headed "Dissociation Energies and
      !                   Be Values (cm^-1) for Isotopomers of HeH+" and its
      !                   footnote reads "D_e values are derived from the
      !                   results of the global fit in Table 4"; Table 8 gives
      !                   the vibrational term values G_v of the effective
      !                   potential of each isotopomer.  The paper states the
      !                   channel this D0 belongs to: "Since the dissociation
      !                   limit of the ground state of HeH+ is He + H+",
      !                   which is the reaction the enthalpy below is written
      !                   for.
      !
      !                   THE PAPER NEVER PRINTS A D0.  What it determines is
      !                   De -- 16456.24 +/- 0.1 cm^-1 for the
      !                   Born-Oppenheimer potential, 16448.84 for the 4HeH+
      !                   isotopomer -- so a citation of it for "D0" is a
      !                   citation for the wrong quantity, and the value
      !                   14874.215 cm^-1 that the secondary literature
      !                   attributes to it appears nowhere in it.  That value
      !                   is 7.95 cm^-1 (9.9e-4 eV, 0.053 per cent) below the
      !                   De - G_0 above; where the difference comes from was
      !                   not established, so the arithmetic that IS supported
      !                   by a primary text is what is used here.
      ! ------------------------------------------------------------------ !
      real*8, parameter :: cm_to_eV  = 1.239841984d-4
      real*8, parameter :: eV_to_erg = 1.602176634d-12

      ! The three ionization potentials this module returns to the gas are
      ! the global thresholds, the same energies the photoionization channels
      ! are charged: one definition of each potential in the code.
      real*8, parameter :: IP_H   = e_th_HI   ! H  -> H+  + e-  [eV]
      real*8, parameter :: IP_He  = e_th_HeI  ! He -> He+ + e-  [eV]
      real*8, parameter :: IP_H2  = e_th_H2   ! H2 -> H2+ + e-  [eV]
      ! D0(H3+ -> H2 + H+); Mizus et al. (2019), published text.
      real*8, parameter :: D0_H3p_eV  = 35076.0d0*cm_to_eV
      ! D0(HeH+ -> He + H+) = De - G_0 = 16448.84 - 1566.6764 cm^-1,
      ! Coxon & Hajigeorgiou (1999) Tables 9 and 8; see above.
      real*8, parameter :: D0_HeHp_eV = (16448.84d0 - 1566.6764d0)*cm_to_eV

      ! ------------------------------------------------------------------ !
      ! KEYS OF THE TABLE.  species_formation_energy is indexed by the f_sp
      ! species column of species_table (isp_HI = 1 to isp_CO = 40, the metal
      ! stages filling 7 to 33), so one integer names a species everywhere.
      ! Three participants of the ledger have no f_sp column of their own and
      ! are given keys past the end of that layout: the free electron, which
      ! is an identity of the charge balance rather than a solved column, and
      ! the two excited states whose populations live in arrays of their own.
      ! ------------------------------------------------------------------ !
      integer, parameter :: isp_eps_electron = 41
      integer, parameter :: isp_eps_H_n2     = 42   ! H(n=2), 2s and 2p
      integer, parameter :: isp_eps_O1D      = 43   ! O(1D)

      ! ------------------------------------------------------------------ !
      ! THE REACTIONS OF THE NETWORK, as stoichiometry rather than as
      ! energies.  Column ir holds up to three reactant and three product
      ! species keys, 0 meaning an unused slot; the third body M of R12, R13
      ! and R15 is on both sides and is left out of both.  The heat of ir is
      ! molecular_reaction_energy_eV(ir), formed from the table above, so no
      ! reaction energy is written down anywhere in this file.
      ! ------------------------------------------------------------------ !
      integer, parameter :: n_mreac = 17
      integer, parameter :: ir_R5  =  1, ir_R6  =  2, ir_R7  =  3,        &
                            ir_R8  =  4, ir_R9  =  5, ir_R10 =  6,        &
                            ir_R11 =  7, ir_R12 =  8, ir_R13 =  9,        &
                            ir_R14 = 10, ir_R15 = 11, ir_R16 = 12,        &
                            ir_R17 = 13, ir_R18 = 14, ir_R19 = 15,        &
                            ir_R20 = 16, ir_R23 = 17
      character(len=26), parameter :: mreac_name(n_mreac) = (/            &
        'R5  H2+ + e  -> H + H     ', 'R6  H3+ + e  -> H2 + H    ',       &
        'R7  H3+ + e  -> H + H + H ', 'R8  H2+ + H2 -> H3+ + H   ',       &
        'R9  H2+ + H  -> H+ + H2   ', 'R10 H+ + H2  -> H2+ + H   ',       &
        'R11 H3+ + H  -> H2+ + H2  ', 'R12 H2 + M   -> H + H + M ',       &
        'R13 H+ +H2+M -> H3+ + M   ', 'R14 H2 + e   -> H + H + e ',       &
        'R15 H + H +M -> H2 + M    ', 'R16 HeH+ + e -> He + H    ',       &
        'R17 He+ + H2 -> H+ + H +He', 'R18 HeH+ +H2 -> H3+ + He  ',       &
        'R19 HeH+ + H -> H2+ + He  ', 'R20 He+ + H2 -> HeH+ + H  ',       &
        'R23 H2 + He+ -> H2+ + He  ' /)
      integer, parameter :: mreac_react(3, n_mreac) = reshape( (/         &
        isp_H2p,  isp_eps_electron, 0,                                    &
        isp_H3p,  isp_eps_electron, 0,                                    &
        isp_H3p,  isp_eps_electron, 0,                                    &
        isp_H2p,  isp_H2,           0,                                    &
        isp_H2p,  isp_HI,           0,                                    &
        isp_HII,  isp_H2,           0,                                    &
        isp_H3p,  isp_HI,           0,                                    &
        isp_H2,   0,                0,                                    &
        isp_HII,  isp_H2,           0,                                    &
        isp_H2,   isp_eps_electron, 0,                                    &
        isp_HI,   isp_HI,           0,                                    &
        isp_HeHp, isp_eps_electron, 0,                                    &
        isp_HeII, isp_H2,           0,                                    &
        isp_HeHp, isp_H2,           0,                                    &
        isp_HeHp, isp_HI,           0,                                    &
        isp_HeII, isp_H2,           0,                                    &
        isp_H2,   isp_HeII,         0 /), (/ 3, n_mreac /) )
      integer, parameter :: mreac_prod(3, n_mreac) = reshape( (/          &
        isp_HI,   isp_HI,           0,                                    &
        isp_H2,   isp_HI,           0,                                    &
        isp_HI,   isp_HI,           isp_HI,                               &
        isp_H3p,  isp_HI,           0,                                    &
        isp_HII,  isp_H2,           0,                                    &
        isp_H2p,  isp_HI,           0,                                    &
        isp_H2p,  isp_H2,           0,                                    &
        isp_HI,   isp_HI,           0,                                    &
        isp_H3p,  0,                0,                                    &
        isp_HI,   isp_HI,           isp_eps_electron,                     &
        isp_H2,   0,                0,                                    &
        isp_HeI,  isp_HI,           0,                                    &
        isp_HII,  isp_HI,           isp_HeI,                              &
        isp_H3p,  isp_HeI,          0,                                    &
        isp_H2p,  isp_HeI,          0,                                    &
        isp_HeHp, isp_HI,           0,                                    &
        isp_H2p,  isp_HeI,          0 /), (/ 3, n_mreac /) )

      ! ------------------------------------------------------------------ !
      ! THE COLLISIONAL OXYGEN CHANNELS AND THE ASSOCIATIVE He(2^3S) BRANCH.
      !
      ! These are reactions of the SAME species table; they are kept in a
      ! table of their own only because their rate coefficients come from
      ! oxygen_rates and from the metastable collision rate of the ionization
      ! sweep, not from mol_rates.  Their energies are formed the same way,
      ! as a difference of reservoir contents, so no reaction enthalpy is
      ! written down here either.
      !
      ! O1, O1r, O2 and O2r are the four collisional channels of the oxygen
      ! network (docs/a2_reaction_audit.md section 2); the photolysis channels
      ! O3, O4, O5 and O7 are paid by the absorbed photon and are deposited by
      ! water_photolysis through heat_per_water_dissociation and
      ! heat_per_hydroxyl_dissociation, so they are absent here.
      !
      ! O6 is the sink of the eliminated O(1D).  The elimination transfers the
      ! reservoir with the nuclei (docs/b1_target_system_20260906.md T1.9
      ! item 1): the flux through O6 equals the O(1D) production oj4 n(H2O),
      ! and it carries eps(O(1D)) to the products.  Writing the channel with
      ! eps(O(1D)) among the reactants is what makes the O4 + O6 pair deposit
      ! exactly as much as the direct route H2O + hv -> OH + H at the same
      ! photon energy, which is acceptance test AT-1d (ii).
      !
      ! The associative branch He(2^3S) + H -> HeH+ + e is the other channel
      ! whose energy had nowhere to go (D0 C25).  Its heat is the difference
      ! eps(He 2^3S) + eps(H) - eps(HeH+) - eps(e), and nothing about it is
      ! transcribed.
      !
      ! D1 is the CO destruction channel He+ + CO -> C+ + O + He (UMIST
      ! RATE22 entry 4068; the rate coefficient is oxygen_rates::
      ! rk_D1_Hep_CO, which states its provenance and its extrapolation).
      ! Its heat is eps(He+) + eps(CO) - eps(C+) - eps(O) - eps(He), which
      ! the table makes +2.2117 eV: the ionization potential of helium, less
      ! that of carbon, less the C=O bond energy.  The reaction is
      ! exothermic, its products leave with that energy as kinetic energy,
      ! and it is deposited as heat.
      !
      ! WHAT THE ROW DOES NOT SAY, and it belongs here rather than in a
      ! reader's head: the carrier step that runs this channel holds the
      ! ion stages frozen, so the row moves NUCLEI and not charge.  The
      ! carbon nucleus released here returns to the carbon element pool and
      ! the ionization sweep that follows re-solves C I / C II / C III from
      ! its own balance; the helium ion consumed here is not taken out of
      ! n(He+) inside the step.  The domain record of
      ! diffusive_photochemistry counts the cells in which that second
      ! omission would matter, i.e. where this channel removes He+ faster
      ! than recombination does.
      ! ------------------------------------------------------------------ !
      integer, parameter :: isp_OI_col  = mion_fsp(im_OI)
      integer, parameter :: isp_CII_col = mion_fsp(im_CII)
      integer, parameter :: n_oreac = 7
      integer, parameter :: ir_O1  = 1, ir_O1r = 2, ir_O2  = 3,           &
                            ir_O2r = 4, ir_O6  = 5, ir_assoc_HeTR = 6,    &
                            ir_D1  = 7
      character(len=26), parameter :: oreac_name(n_oreac) = (/            &
        'O1  OH + H2  -> H2O + H   ', 'O1r H2O + H  -> OH + H2   ',       &
        'O2  O + H2   -> OH + H    ', 'O2r OH + H   -> O + H2    ',       &
        'O6  O(1D)+H2 -> OH + H    ', 'A1  He3S + H -> HeH+ + e  ',       &
        'D1  He+ + CO -> C+ +O +He ' /)
      integer, parameter :: oreac_react(3, n_oreac) = reshape( (/         &
        isp_OH,       isp_H2, 0,                                          &
        isp_H2O,      isp_HI, 0,                                          &
        isp_OI_col,   isp_H2, 0,                                          &
        isp_OH,       isp_HI, 0,                                          &
        isp_eps_O1D,  isp_H2, 0,                                          &
        isp_HeTR,     isp_HI, 0,                                          &
        isp_HeII,     isp_CO, 0 /), (/ 3, n_oreac /) )
      integer, parameter :: oreac_prod(3, n_oreac) = reshape( (/          &
        isp_H2O,      isp_HI,           0,                                &
        isp_OH,       isp_H2,           0,                                &
        isp_OH,       isp_HI,           0,                                &
        isp_OI_col,   isp_H2,           0,                                &
        isp_OH,       isp_HI,           0,                                &
        isp_HeHp,     isp_eps_electron, 0,                                &
        isp_CII_col,  isp_OI_col,       isp_HeI /), (/ 3, n_oreac /) )

      contains

      ! ------------------------------------------------------!

      subroutine molecular_chemical_heating(T_K, nhi, nhii, nheii, nmol,   &
                                            ne, ntot, gamma_chem)
      ! Volumetric chemical heating rate [erg cm^-3 s^-1] of every cell.
      ! nmol(:,1:4) are n(H2), n(H2+), n(H3+), n(HeH+) in cm^-3, the same
      ! columns ionization_equilibrium already carries; everything else is
      ! cgs as well.  Positive is heating.
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: T_K, nhi, nhii, nheii
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: ne, ntot
      real*8, dimension(1-Ng:N+Ng,4), intent(in)  :: nmol
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: gamma_chem
      real*8 :: q(n_mreac)
      real*8 :: T, nh2, nh2p, nh3p, nhehp, f_heat
      real*8 :: k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19
      real*8 :: k20,k23
      real*8 :: g_eV
      integer :: j, ir

      ! Heat released by one event of each reaction [eV], from the one
      ! formation-energy table.  Temperature-independent, so it is formed
      ! once for the whole grid.
      do ir = 1, n_mreac
         q(ir) = molecular_reaction_energy_eV(ir)
      enddo

      do j = 1-Ng, N+Ng
         T     = max(T_K(j), 1.0d0)
         nh2   = nmol(j,1);  nh2p  = nmol(j,2)
         nh3p  = nmol(j,3);  nhehp = nmol(j,4)
         k5  = rk_R5_H2p_dr(T);        k6  = rk_R6_H3p_dr_H2(T)
         k7  = rk_R7_H3p_dr_3H(T);     k8  = rk_R8_H2p_H2()
         k9  = rk_R9_H2p_H();          k10 = rk_R10_Hp_H2v4(T)
         k11 = rk_R11_H3p_H(T);        k12 = rk_R12_H2_thdis(T)
         k13 = rk_R13_Hp_H2_M(ntot(j));k14 = rk_R14_H2_edis(T)
         k15 = rk_R15_3body_H2(T, ntot(j))
         k16 = rk_R16_HeHp_dr(T);      k17 = rk_R17_Hep_H2_diss(T)
         k18 = rk_R18_HeHp_H2();       k19 = rk_R19_HeHp_H()
         k20 = rk_R20_Hep_H2_HeHp();   k23 = rk_R23_H2_Hep_cx()
         ! Share of the association energy that thermalizes rather than being
         ! radiated in the infrared quadrupole lines.
         f_heat = h2_vibrational_heat_fraction(T, nhi(j), nh2)
         ! The heat of each channel is q(ir) above; what is written here is
         ! only which densities and which rate coefficient multiply it.
         ! R10 is endothermic and consumes vibrational energy that came from
         ! the thermal bath, so its sink is the ground-state difference like
         ! every other channel; R15 deposits only the share f_heat that
         ! thermalizes rather than leaving in the infrared quadrupole lines.
         g_eV =                                                            &
              k5 *nh2p*ne(j)          *q(ir_R5)                            &
            + k6 *nh3p*ne(j)          *q(ir_R6)                            &
            + k7 *nh3p*ne(j)          *q(ir_R7)                            &
            + k8 *nh2p*nh2            *q(ir_R8)                            &
            + k9 *nh2p*nhi(j)         *q(ir_R9)                            &
            + k10*nhii(j)*nh2         *q(ir_R10)                           &
            + k11*nh3p*nhi(j)         *q(ir_R11)                           &
            + k12*ntot(j)*nh2         *q(ir_R12)                           &
            + k13*nhii(j)*nh2         *q(ir_R13)                           &
            + k14*ne(j)*nh2           *q(ir_R14)                           &
            + k15*nhi(j)*nhi(j)       *q(ir_R15)*f_heat                    &
            + k16*nhehp*ne(j)         *q(ir_R16)                           &
            + k17*nheii(j)*nh2        *q(ir_R17)                           &
            + k18*nhehp*nh2           *q(ir_R18)                           &
            + k19*nhehp*nhi(j)        *q(ir_R19)                           &
            + k20*nheii(j)*nh2        *q(ir_R20)                           &
            + k23*nh2*nheii(j)        *q(ir_R23)
         gamma_chem(j) = g_eV*eV_to_erg
      enddo

      end subroutine molecular_chemical_heating

      ! ------------------------------------------------------!

      double precision function species_formation_energy(isp) result(eps)
      ! Formation plus excitation energy [eV] of one particle of species
      ! isp, measured from the single reference state at the top of this
      ! module: every element a neutral, ground-state, free atom at rest,
      ! the free electron at zero.  isp is the f_sp species column of
      ! species_table, or one of the three keys isp_eps_* for the
      ! participants that have no column.
      !
      ! This is the eps_s of the energy ledger, and it is the ONLY place any
      ! of these energies is written down: a reaction heat, a photoevent
      ! partition and the reservoir density sum_s n_s eps_s all read it, so
      ! none of them can disagree with the others about how much energy a
      ! particle holds.
      integer, intent(in) :: isp
      integer :: im, k
      real*8  :: D0_H2

      select case (isp)
      ! ---- the reference and the free electron ----
      case (isp_HI, isp_HeI, isp_eps_electron)
         eps = 0.0d0
      ! ---- H and He ionization, cumulative over the stages ----
      case (isp_HII)
         eps = IP_H
      case (isp_HeII)
         eps = IP_He
      case (isp_HeIII)
         eps = IP_He + e_th_HeII
      ! ---- excited states: excitation energy only, the parent's nucleus
      !      and translational energy being counted through the parent ----
      case (isp_HeTR)                 ! He 2^3S above He I ground
         eps = IP_He - e_th_HeTR
      case (isp_eps_H_n2)             ! H(n=2), 2s and 2p, above H(1s)
         eps = IP_H*(1.0d0 - 0.25d0)
      case (isp_eps_O1D)              ! O(1D) above the O(3P) ground term
         eps = e_excite_O1D_erg/eV_to_erg
      ! ---- the molecular species of the H2/He network ----
      case (isp_H2)
         eps = -h2_dissociation_energy_eV()
      case (isp_H2p)                  ! H + H -> H2+ + e
         D0_H2 = h2_dissociation_energy_eV()
         eps = IP_H2 - D0_H2
      case (isp_H3p)                  ! H2 + H+ -> H3+
         D0_H2 = h2_dissociation_energy_eV()
         eps = -D0_H2 + IP_H - D0_H3p_eV
      case (isp_HeHp)                 ! He + H+ -> HeH+
         eps = IP_H - D0_HeHp_eV
      ! ---- the oxygen and carbon molecules, from the 0 K formation
      !      energies of the thermodynamic table that owns them ----
      case (isp_OH)
         eps = oxygen_formation_energy_eV(ith_OH)
      case (isp_H2O)
         eps = oxygen_formation_energy_eV(ith_H2O)
      case (isp_CO)
         eps = oxygen_formation_energy_eV(ith_CO)
      case default
         ! ---- the metal stages, f_sp columns 7 to 33 ----
         ! Cumulative ionization potentials from the neutral ground state,
         ! summed over the stages of the ion's own element, so the neutral
         ! is at zero for every element and no zero is placed at whichever
         ! stage dominates the base.
         im = isp - 6
         if (im .ge. 1 .and. im .le. n_mion) then
            eps = 0.0d0
            do k = melem_i0(mion_elem(im)), im - 1
               eps = eps + mion_ethr(k)
            enddo
         else
            write(*,*) 'species_formation_energy: no entry for species ',  &
                       isp
            stop
         endif
      end select

      end function species_formation_energy

      ! ------------------------------------------------------!

      double precision function molecular_reaction_energy_eV(ir) result(q)
      ! Heat [eV] released by one event of reaction ir of the H2/He network,
      ! as the formation energy of its reactants minus that of its products.
      ! Positive is exothermic.  Nothing is transcribed: the stoichiometry
      ! table and the species table together fix every channel, so a forward
      ! and a reverse channel are exact negatives and a closed cycle sums to
      ! zero.
      integer, intent(in) :: ir
      integer :: i

      q = 0.0d0
      do i = 1, 3
         if (mreac_react(i, ir) .ne. 0)                                    &
            q = q + species_formation_energy(mreac_react(i, ir))
      enddo
      do i = 1, 3
         if (mreac_prod(i, ir) .ne. 0)                                     &
            q = q - species_formation_energy(mreac_prod(i, ir))
      enddo

      end function molecular_reaction_energy_eV

      ! ------------------------------------------------------!

      double precision function oxygen_reaction_energy_eV(ir) result(q)
      ! Heat [eV] released by one event of channel ir of the collisional
      ! oxygen network, or of the associative He(2^3S) branch, as the
      ! formation energy of its reactants minus that of its products.  Same
      ! construction as molecular_reaction_energy_eV and the same table, so
      ! O1 and O1r are exact negatives and the O4 + O6 pair deposits what the
      ! direct dissociation of H2O deposits.
      integer, intent(in) :: ir
      integer :: i

      q = 0.0d0
      do i = 1, 3
         if (oreac_react(i, ir) .ne. 0)                                    &
            q = q + species_formation_energy(oreac_react(i, ir))
      enddo
      do i = 1, 3
         if (oreac_prod(i, ir) .ne. 0)                                     &
            q = q - species_formation_energy(oreac_prod(i, ir))
      enddo

      end function oxygen_reaction_energy_eV

      ! ------------------------------------------------------!

      subroutine oxygen_chemical_heating(T_K, nhi, nh2, n_o0, nox,        &
                                         flux_O6, gamma_ox)
      ! Volumetric heating rate [erg cm^-3 s^-1] of the four COLLISIONAL
      ! oxygen channels O1, O1r, O2, O2r and of the O(1D) sink O6.  Positive
      ! is heating.  The photolysis channels are not here: their enthalpy is
      ! paid by the absorbed photon and water_photolysis deposits the excess.
      !
      ! n_o0 is the FREE atomic oxygen of the cell, the O I population the
      ! oxygen carriers and CO have been taken out of, which is the density
      ! the O2 channel runs on in System_HeH_mol::oxygen_carrier_rows.
      ! flux_O6 is the O(1D) production rate oj4 n(H2O) [cm^-3 s^-1]: O(1D)
      ! is eliminated by its local steady state, so the flux through its
      ! single sink O6 equals its production whatever n(H2) is, and the
      ! caller supplies that product because it owns the band photolysis
      ! rates.
      !
      ! The rate coefficients are the same functions of oxygen_rates that
      ! System_HeH_mol::set_oxygen_coeffs evaluates, with the same arguments
      ! and the same detailed-balance reverses, so there is one definition of
      ! each coefficient in the code and this is a second evaluation of it,
      ! not a second transcription.
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: T_K, nhi, nh2, n_o0
      real*8, dimension(1-Ng:N+Ng,3), intent(in)  :: nox
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: flux_O6
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: gamma_ox
      real*8  :: q(n_oreac)
      real*8  :: T, k1, k1r, k2, k2r, g_eV
      integer :: j, ir

      do ir = 1, n_oreac
         q(ir) = oxygen_reaction_energy_eV(ir)
      enddo

      do j = 1-Ng, N+Ng
         T   = max(T_K(j), 1.0d0)
         k1  = rk_O1_OH_H2_water(T)
         k1r = rate_from_detailed_balance(k1, (/ ith_OH, ith_H2 /),        &
                                              (/ ith_H2O, ith_H /), T)
         k2  = rk_O2_O_H2_hydroxyl(T)
         k2r = rate_from_detailed_balance(k2, (/ ith_O, ith_H2 /),         &
                                              (/ ith_OH, ith_H /), T)
         g_eV = k1 *nox(j,1)*nh2(j)   *q(ir_O1)                            &
              + k1r*nox(j,2)*nhi(j)   *q(ir_O1r)                           &
              + k2 *n_o0(j) *nh2(j)   *q(ir_O2)                            &
              + k2r*nox(j,1)*nhi(j)   *q(ir_O2r)                           &
              + flux_O6(j)            *q(ir_O6)
         gamma_ox(j) = g_eV*eV_to_erg
      enddo

      end subroutine oxygen_chemical_heating

      ! ------------------------------------------------------!

      subroutine formation_energy_density(rho, f_sp, u_form, n_H_n2, n_O1D)
      ! The formation/excitation reservoir u_form = sum_s n_s eps_s
      ! [erg cm^-3] of every cell (docs/b1_target_system_20260906.md T1.2).
      !
      ! rho is the adimensional mass density and f_sp the species fractions,
      ! so n_s = rho f_sp(s) n0 in cm^-3, the same conversion
      ! get_species_densities makes.
      !
      ! THE TWO EXCITED POPULATIONS ARE ADDED ON TOP OF THEIR PARENT COLUMN,
      ! not instead of it.  The He I column of f_sp contains the 2^3S
      ! metastable and the H I column contains H(n=2), and both parents carry
      ! eps = 0, so a metastable atom contributes exactly its excitation
      ! energy and its nucleus is counted once.  The same holds for O(1D)
      ! inside the O I column.  Optional because a run without the excited
      ! hydrogen, or without the oxygen chemistry, has no such population.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: u_form
      real*8, dimension(1-Ng:N+Ng), intent(in), optional  :: n_H_n2, n_O1D
      real*8  :: eps(n_species)
      integer :: isp

      do isp = 1, n_species
         eps(isp) = species_formation_energy(isp)*eV_to_erg
      enddo

      u_form = 0.0d0
      do isp = 1, n_species
         if (eps(isp) .eq. 0.0d0) cycle
         u_form = u_form + eps(isp)*rho*f_sp(:,isp)*n0
      enddo

      if (present(n_H_n2)) u_form = u_form                                 &
         + species_formation_energy(isp_eps_H_n2)*eV_to_erg*n_H_n2
      if (present(n_O1D))  u_form = u_form                                 &
         + species_formation_energy(isp_eps_O1D) *eV_to_erg*n_O1D

      end subroutine formation_energy_density

      ! ------------------------------------------------------!

      subroutine formation_energy_report(iu)
      ! Write the table and the reaction energies it implies, so that a run
      ! records the ledger it used rather than leaving it to be re-derived.
      integer, intent(in) :: iu
      integer :: ir

      write(iu,'(A)') '# Species formation energies eps_s [eV], zero at'
      write(iu,'(A)') '# neutral ground-state free atoms and free e- at rest'
      write(iu,'(A,F10.4)') '#   eps(H+)   = ',                            &
                            species_formation_energy(isp_HII)
      write(iu,'(A,F10.4)') '#   eps(He+)  = ',                            &
                            species_formation_energy(isp_HeII)
      write(iu,'(A,F10.4)') '#   eps(He++) = ',                            &
                            species_formation_energy(isp_HeIII)
      write(iu,'(A,F10.4)') '#   eps(He2^3S) = ',                          &
                            species_formation_energy(isp_HeTR)
      write(iu,'(A,F10.4)') '#   eps(H n=2) = ',                           &
                            species_formation_energy(isp_eps_H_n2)
      write(iu,'(A,F10.4)') '#   eps(H2)   = ',                            &
                            species_formation_energy(isp_H2)
      write(iu,'(A,F10.4)') '#   eps(H2+)  = ',                            &
                            species_formation_energy(isp_H2p)
      write(iu,'(A,F10.4)') '#   eps(H3+)  = ',                            &
                            species_formation_energy(isp_H3p)
      write(iu,'(A,F10.4)') '#   eps(HeH+) = ',                            &
                            species_formation_energy(isp_HeHp)
      write(iu,'(A,F10.4)') '#   eps(OH)   = ',                            &
                            species_formation_energy(isp_OH)
      write(iu,'(A,F10.4)') '#   eps(H2O)  = ',                            &
                            species_formation_energy(isp_H2O)
      write(iu,'(A,F10.4)') '#   eps(CO)   = ',                            &
                            species_formation_energy(isp_CO)
      write(iu,'(A,F10.4)') '#   eps(O 1D) = ',                            &
                            species_formation_energy(isp_eps_O1D)
      write(iu,'(A)') '# Heat released [eV] per reaction:'
      do ir = 1, n_mreac
         write(iu,'(A,A,A,F9.4)') '#   ', mreac_name(ir), ' : ',           &
                                  molecular_reaction_energy_eV(ir)
      enddo
      do ir = 1, n_oreac
         write(iu,'(A,A,A,F9.4)') '#   ', oreac_name(ir), ' : ',           &
                                  oxygen_reaction_energy_eV(ir)
      enddo

      end subroutine formation_energy_report

      ! ------------------------------------------------------!

      ! End of module
      end module molecular_reaction_heat

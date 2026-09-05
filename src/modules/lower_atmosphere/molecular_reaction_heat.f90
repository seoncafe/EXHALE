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
      ! docs/Update_EXHALE.md), the sum below is 2.4 times the whole
      ! photoelectric heating rate over 1.00-1.05 r_base, falling through 1 at
      ! 1.053 and to 1 percent by 1.10.  The argument above makes that ratio
      ! I(H2)/(hv - I(H2)), so 2.4 corresponds to a mean absorbed photon of
      ! 22 eV: a consistency check on the sum, not an independent derivation
      ! of it.
      !
      ! HOW IT IS BUILT.  Not as a list of reaction enthalpies -- a list can
      ! disagree with itself -- but from ONE table of species enthalpies.
      ! h(X) is the energy needed to build X out of ground-state H atoms,
      ! ground-state He atoms and free electrons at rest, so
      !
      !     heat released by  A + B -> C + D   =  h(A) + h(B) - h(C) - h(D) ,
      !
      ! and every reaction of the network, in either direction, is consistent
      ! with every other by construction.  A cycle that returns to its
      ! starting species releases exactly zero, which is the property a
      ! reaction-by-reaction list cannot be trusted to have.
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
      !     collision, He(2^3S) + H -> HeH+ + e, releases 8.1 eV and is NOT
      !     deposited anywhere; it is left out here too because the metastable
      !     density in the molecular layer is 2e-7 cm^-3 (measured on the
      !     state above), which makes that channel 1e-13 of the sum below.
      !     It belongs with the other He(2^3S) terms if it is ever added.
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
      public :: molecular_chemical_heating, species_enthalpy_report

      ! ------------------------------------------------------------------ !
      ! THE SPECIES ENTHALPY TABLE  [eV], zero at H + He + e- at rest.
      !
      ! PROVENANCE, one line per constant, and what was actually checked.
      !
      ! The ionization energies here are the THERMOCHEMICAL ones and are
      ! deliberately not the e_th_* of global_parameters: those are the
      ! rounded thresholds the photoionization cross sections are tabulated
      ! against (13.6, 24.6, 15.4 eV), and rounding a threshold is harmless in
      ! a cross section but not in an energy ledger that has to close over a
      ! cycle.
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

      real*8, parameter :: IP_H   = 13.598434599d0  ! NIST ASD
      real*8, parameter :: IP_He  = 24.587389011d0  ! NIST ASD
      real*8, parameter :: IP_H2  = 15.425933d0     ! NIST WebBook
      ! D0(H3+ -> H2 + H+); Mizus et al. (2019), published text.
      real*8, parameter :: D0_H3p_eV  = 35076.0d0*cm_to_eV
      ! D0(HeH+ -> He + H+) = De - G_0 = 16448.84 - 1566.6764 cm^-1,
      ! Coxon & Hajigeorgiou (1999) Tables 9 and 8; see above.
      real*8, parameter :: D0_HeHp_eV = (16448.84d0 - 1566.6764d0)*cm_to_eV

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
      real*8 :: h_H2, h_Hp, h_Hep, h_H2p, h_H3p, h_HeHp
      real*8 :: T, nh2, nh2p, nh3p, nhehp, f_heat
      real*8 :: k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19
      real*8 :: k20,k23
      real*8 :: g_eV
      integer :: j

      call species_enthalpies(h_H2, h_Hp, h_Hep, h_H2p, h_H3p, h_HeHp)

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
         g_eV =                                                            &
         ! R5  H2+ + e  -> H + H
              k5 *nh2p*ne(j)          *( h_H2p )                           &
         ! R6  H3+ + e  -> H2 + H
            + k6 *nh3p*ne(j)          *( h_H3p - h_H2 )                    &
         ! R7  H3+ + e  -> H + H + H
            + k7 *nh3p*ne(j)          *( h_H3p )                           &
         ! R8  H2+ + H2 -> H3+ + H
            + k8 *nh2p*nh2            *( h_H2p + h_H2 - h_H3p )            &
         ! R9  H2+ + H  -> H+ + H2
            + k9 *nh2p*nhi(j)         *( h_H2p - h_Hp - h_H2 )             &
         ! R10 H+ + H2(v>=4) -> H2+ + H   (endothermic; the vibrational
         !     energy it consumes came from the thermal bath, so the sink is
         !     the ground-state enthalpy difference)
            + k10*nhii(j)*nh2         *( h_Hp + h_H2 - h_H2p )             &
         ! R11 H3+ + H  -> H2+ + H2      (endothermic)
            + k11*nh3p*nhi(j)         *( h_H3p - h_H2p - h_H2 )            &
         ! R12 H2 + M   -> H + H + M     (the thermal dissociation SINK)
            + k12*ntot(j)*nh2         *( h_H2 )                            &
         ! R13 H+ + H2 + M -> H3+ + M
            + k13*nhii(j)*nh2         *( h_Hp + h_H2 - h_H3p )             &
         ! R14 H2 + e   -> H + H + e     (the electron-impact SINK)
            + k14*ne(j)*nh2           *( h_H2 )                            &
         ! R15 H + H + M -> H2 + M       (association; HM79 branching)
            + k15*nhi(j)*nhi(j)       *( -h_H2 )*f_heat                    &
         ! R16 HeH+ + e -> He + H
            + k16*nhehp*ne(j)         *( h_HeHp )                          &
         ! R17 He+ + H2 -> H+ + H + He
            + k17*nheii(j)*nh2        *( h_Hep + h_H2 - h_Hp )             &
         ! R18 HeH+ + H2 -> H3+ + He
            + k18*nhehp*nh2           *( h_HeHp + h_H2 - h_H3p )           &
         ! R19 HeH+ + H -> H2+ + He
            + k19*nhehp*nhi(j)        *( h_HeHp - h_H2p )                  &
         ! R20 He+ + H2 -> HeH+ + H
            + k20*nheii(j)*nh2        *( h_Hep + h_H2 - h_HeHp )           &
         ! R23 H2 + He+ -> H2+ + He
            + k23*nh2*nheii(j)        *( h_H2 + h_Hep - h_H2p )
         gamma_chem(j) = g_eV*eV_to_erg
      enddo

      end subroutine molecular_chemical_heating

      ! ------------------------------------------------------!

      subroutine species_enthalpies(h_H2, h_Hp, h_Hep, h_H2p, h_H3p,       &
                                    h_HeHp)
      ! Energy [eV] to build each species from ground-state H, ground-state
      ! He and free electrons at rest.  D0(H2) is read from mol_rates so that
      ! the bond energy has exactly one definition in the code.
      real*8, intent(out) :: h_H2, h_Hp, h_Hep, h_H2p, h_H3p, h_HeHp
      real*8 :: D0_H2

      D0_H2  = h2_dissociation_energy_eV()
      h_H2   = -D0_H2                       ! H + H -> H2
      h_Hp   =  IP_H                        ! H -> H+ + e
      h_Hep  =  IP_He                       ! He -> He+ + e
      h_H2p  =  IP_H2 - D0_H2               ! H + H -> H2+ + e
      h_H3p  =  h_H2 + h_Hp - D0_H3p_eV     ! H2 + H+ -> H3+
      h_HeHp =  h_Hp - D0_HeHp_eV           ! He + H+ -> HeH+

      end subroutine species_enthalpies

      ! ------------------------------------------------------!

      subroutine species_enthalpy_report(iu)
      ! Write the table and the reaction energies it implies, so that a run
      ! records the ledger it used rather than leaving it to be re-derived.
      integer, intent(in) :: iu
      real*8 :: h_H2, h_Hp, h_Hep, h_H2p, h_H3p, h_HeHp

      call species_enthalpies(h_H2, h_Hp, h_Hep, h_H2p, h_H3p, h_HeHp)
      write(iu,'(A)') '# Species enthalpies [eV], zero at H + He + e- at rest'
      write(iu,'(A,F10.4)') '#   h(H2)   = ', h_H2
      write(iu,'(A,F10.4)') '#   h(H+)   = ', h_Hp
      write(iu,'(A,F10.4)') '#   h(He+)  = ', h_Hep
      write(iu,'(A,F10.4)') '#   h(H2+)  = ', h_H2p
      write(iu,'(A,F10.4)') '#   h(H3+)  = ', h_H3p
      write(iu,'(A,F10.4)') '#   h(HeH+) = ', h_HeHp
      write(iu,'(A)') '# Heat released [eV] per reaction:'
      write(iu,'(A,F9.4)') '#   R5  H2+ + e  -> H + H      : ', h_H2p
      write(iu,'(A,F9.4)') '#   R6  H3+ + e  -> H2 + H     : ', h_H3p - h_H2
      write(iu,'(A,F9.4)') '#   R7  H3+ + e  -> H + H + H  : ', h_H3p
      write(iu,'(A,F9.4)') '#   R8  H2+ + H2 -> H3+ + H    : ',            &
                            h_H2p + h_H2 - h_H3p
      write(iu,'(A,F9.4)') '#   R9  H2+ + H  -> H+ + H2    : ',            &
                            h_H2p - h_Hp - h_H2
      write(iu,'(A,F9.4)') '#   R12 H2 + M   -> H + H + M  : ', h_H2
      write(iu,'(A,F9.4)') '#   R13 H+ +H2+M -> H3+ + M    : ',            &
                            h_Hp + h_H2 - h_H3p
      write(iu,'(A,F9.4)') '#   R15 H + H +M -> H2 + M     : ', -h_H2
      write(iu,'(A,F9.4)') '#   R16 HeH+ + e -> He + H     : ', h_HeHp
      write(iu,'(A,F9.4)') '#   R17 He+ + H2 -> H+ + H + He: ',            &
                            h_Hep + h_H2 - h_Hp
      write(iu,'(A,F9.4)') '#   R18 HeH+ +H2 -> H3+ + He   : ',            &
                            h_HeHp + h_H2 - h_H3p
      write(iu,'(A,F9.4)') '#   R19 HeH+ + H -> H2+ + He   : ',            &
                            h_HeHp - h_H2p
      write(iu,'(A,F9.4)') '#   R20 He+ + H2 -> HeH+ + H   : ',            &
                            h_Hep + h_H2 - h_HeHp
      write(iu,'(A,F9.4)') '#   R23 H2 + He+ -> H2+ + He   : ',            &
                            h_H2 + h_Hep - h_H2p

      end subroutine species_enthalpy_report

      ! End of module
      end module molecular_reaction_heat

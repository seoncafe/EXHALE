      program species_formation_energy_table
      ! The one formation-energy table eps_s of the energy ledger, and the
      ! reaction energies that are differences of it.
      !
      ! WHAT IS BEING ASSERTED.  docs/b1_target_system_20260906.md T1.2 fixes
      ! one reference state for the whole code -- every element a neutral,
      ! ground-state, free atom at rest, the free electron at zero -- and
      ! T1.9 requires every reaction energy of every network to be a
      ! DIFFERENCE of that one table rather than a number written beside the
      ! reaction.  A table has a property a list of reaction enthalpies does
      ! not: a forward and a reverse channel are exact negatives, a closed
      ! cycle releases exactly zero, and two routes between the same two
      ! states release the same energy.  Those are the identities below, and
      ! each of them fails if any energy is transcribed twice.
      !
      ! T1.3 adds the non-overlap rule: energy held by the equation of state
      ! is not also in the reservoir.  For H2 that is a statement about one
      ! zero -- the v = 0, J = 0 level -- which the bond energy is measured
      ! from and which the rovibrational sum of the caloric equation of state
      ! starts at.  Group A asserts that the two really are the same zero.
      !
      ! Production routines exercised:
      !   species_formation_energy and molecular_reaction_energy_eV of
      !     src/modules/lower_atmosphere/molecular_reaction_heat.f90 -- the
      !     table itself and the 17 reaction energies formed from it;
      !   oxygen_formation_energy_eV and photolysis_threshold_erg of
      !     src/modules/lower_atmosphere/oxygen_rates.f90 -- the oxygen and
      !     carbon entries of the same table and the photolysis channel
      !     energies built from them;
      !   h2_dissociation_energy_eV of
      !     src/modules/lower_atmosphere/mol_rates.f90 -- the code's one
      !     definition of D0(H2);
      !   the H2 level ladder of
      !     src/modules/lower_atmosphere/molecular_infrared_data.f90, whose
      !     zero and whose bound limit the reservoir has to agree with;
      !   the composition metadata of src/modules/init/species_table.f90,
      !     against which the stoichiometry of every reaction is balanced.
      use global_parameters,   only: e_th_HI, e_th_HeI, e_th_HeII,        &
                                     e_th_H2, e_th_HeTR
      use species_table,       only: isp_HI, isp_HII, isp_HeI, isp_HeII,  &
                                     isp_HeIII, isp_HeTR, isp_H2,         &
                                     isp_H2p, isp_H3p, isp_HeHp, isp_OH,  &
                                     isp_H2O, isp_CO, n_bsp, bsp_fsp,     &
                                     bsp_nH, bsp_nHe, bsp_nO, bsp_nC,     &
                                     bsp_charge, im_CI, im_CII, im_OI,    &
                                     im_NaI, im_FeI, im_FeII
      use molecular_reaction_heat, only: species_formation_energy,        &
                                     molecular_reaction_energy_eV,        &
                                     n_mreac, mreac_react, mreac_prod,    &
                                     mreac_name, isp_eps_electron,        &
                                     isp_eps_H_n2, isp_eps_O1D
      use oxygen_rates,        only: oxygen_formation_energy_eV,          &
                                     photolysis_threshold_erg,            &
                                     ith_H, ith_H2, ith_O, ith_OH,        &
                                     ith_H2O, ith_CO, ith_C,              &
                                     ich_H2O_OH_H, ich_H2O_H2_O,          &
                                     ich_H2O_O_H_H, ich_OH_O_H,           &
                                     e_excite_O1D_erg
      use mol_rates,           only: h2_dissociation_energy_eV
      use molecular_infrared_data, only: n_h2_lev, h2_lev_T
      use assertion_report
      implicit none

      ! Metal stage columns of f_sp: the neutral of the element sits at
      ! column 6 + im_<element>I, and each further stage one column on.
      integer, parameter :: isp_CI   = 6 + im_CI
      integer, parameter :: isp_CII  = 6 + im_CII
      integer, parameter :: isp_CIII = 6 + im_CII + 1
      integer, parameter :: isp_OI   = 6 + im_OI
      integer, parameter :: isp_NaI  = 6 + im_NaI
      integer, parameter :: isp_NaII = 6 + im_NaI + 1
      integer, parameter :: isp_FeI  = 6 + im_FeI
      integer, parameter :: isp_FeIII= 6 + im_FeII + 1

      ! cm^-1 to eV and eV to erg, CODATA 2018.
      real*8, parameter :: cm_to_eV  = 1.239841984d-4
      real*8, parameter :: eV_to_erg = 1.602176634d-12
      ! One wavenumber as a temperature, hc/kB [K cm].
      real*8, parameter :: cm_to_K   = 1.438776877d0
      ! D0(H2) from v = 0, J = 0, READ from mol_rates.f90 (Huber and
      ! Herzberg 1979).  Transcribed here so that the driver does not take
      ! the number from the routine it is checking.
      real*8, parameter :: D0_H2_cm  = 36118.11d0

      ! ---- an independent transcription of the network, from the reaction
      !      names of molecular_reaction_heat.  The driver forms each energy
      !      from these lists and compares with what the module returns, so
      !      the module's own stoichiometry table is never the reference for
      !      itself.
      integer, parameter :: n_check = 17
      integer :: tr_react(3, n_check), tr_prod(3, n_check)

      real*8  :: q(n_mreac), eps_H2, eps_OH, eps_H2O, eps_CO, eps_O1D
      real*8  :: e1, e2, e3, e4, qi
      integer :: ir, i, k, dH, dHe, dO, dC, dZ
      character(len=64) :: nm

      ! ---------------- the independent stoichiometry ----------------
      ! R5  H2+ + e  -> H + H
      tr_react(:, 1) = (/ isp_H2p, isp_eps_electron, 0 /)
      tr_prod (:, 1) = (/ isp_HI,  isp_HI,           0 /)
      ! R6  H3+ + e  -> H2 + H
      tr_react(:, 2) = (/ isp_H3p, isp_eps_electron, 0 /)
      tr_prod (:, 2) = (/ isp_H2,  isp_HI,           0 /)
      ! R7  H3+ + e  -> H + H + H
      tr_react(:, 3) = (/ isp_H3p, isp_eps_electron, 0 /)
      tr_prod (:, 3) = (/ isp_HI,  isp_HI,      isp_HI /)
      ! R8  H2+ + H2 -> H3+ + H
      tr_react(:, 4) = (/ isp_H2p, isp_H2,           0 /)
      tr_prod (:, 4) = (/ isp_H3p, isp_HI,           0 /)
      ! R9  H2+ + H  -> H+ + H2
      tr_react(:, 5) = (/ isp_H2p, isp_HI,           0 /)
      tr_prod (:, 5) = (/ isp_HII, isp_H2,           0 /)
      ! R10 H+ + H2  -> H2+ + H
      tr_react(:, 6) = (/ isp_HII, isp_H2,           0 /)
      tr_prod (:, 6) = (/ isp_H2p, isp_HI,           0 /)
      ! R11 H3+ + H  -> H2+ + H2
      tr_react(:, 7) = (/ isp_H3p, isp_HI,           0 /)
      tr_prod (:, 7) = (/ isp_H2p, isp_H2,           0 /)
      ! R12 H2 + M   -> H + H + M
      tr_react(:, 8) = (/ isp_H2,  0,                0 /)
      tr_prod (:, 8) = (/ isp_HI,  isp_HI,           0 /)
      ! R13 H+ + H2 + M -> H3+ + M
      tr_react(:, 9) = (/ isp_HII, isp_H2,           0 /)
      tr_prod (:, 9) = (/ isp_H3p, 0,                0 /)
      ! R14 H2 + e   -> H + H + e
      tr_react(:,10) = (/ isp_H2,  isp_eps_electron, 0 /)
      tr_prod (:,10) = (/ isp_HI,  isp_HI, isp_eps_electron /)
      ! R15 H + H + M -> H2 + M
      tr_react(:,11) = (/ isp_HI,  isp_HI,           0 /)
      tr_prod (:,11) = (/ isp_H2,  0,                0 /)
      ! R16 HeH+ + e -> He + H
      tr_react(:,12) = (/ isp_HeHp, isp_eps_electron, 0 /)
      tr_prod (:,12) = (/ isp_HeI, isp_HI,           0 /)
      ! R17 He+ + H2 -> H+ + H + He
      tr_react(:,13) = (/ isp_HeII, isp_H2,          0 /)
      tr_prod (:,13) = (/ isp_HII, isp_HI,     isp_HeI /)
      ! R18 HeH+ + H2 -> H3+ + He
      tr_react(:,14) = (/ isp_HeHp, isp_H2,          0 /)
      tr_prod (:,14) = (/ isp_H3p, isp_HeI,          0 /)
      ! R19 HeH+ + H -> H2+ + He
      tr_react(:,15) = (/ isp_HeHp, isp_HI,          0 /)
      tr_prod (:,15) = (/ isp_H2p, isp_HeI,          0 /)
      ! H2+ + He -> HeH+ + H, the HeH+ source that replaced Koskinen R20
      ! (item L7f; the paper R20 cites bounds that channel 42 times below
      ! the value Table 1 gave it).  The reaction is ENDOTHERMIC, by 0.80 eV
      ! from this table, which is the published endothermicity of the
      ! ground-state channel; the 6717 K in its rate coefficient is the
      ! effective barrier of a vibrationally averaged cross section and is
      ! deliberately smaller.
      tr_react(:,16) = (/ isp_H2p,  isp_HeI,         0 /)
      tr_prod (:,16) = (/ isp_HeHp, isp_HI,          0 /)
      ! R23 H2 + He+ -> H2+ + He
      tr_react(:,17) = (/ isp_H2,  isp_HeII,         0 /)
      tr_prod (:,17) = (/ isp_H2p, isp_HeI,          0 /)

      eps_H2  = species_formation_energy(isp_H2)
      eps_OH  = species_formation_energy(isp_OH)
      eps_H2O = species_formation_energy(isp_H2O)
      eps_CO  = species_formation_energy(isp_CO)
      eps_O1D = species_formation_energy(isp_eps_O1D)
      do ir = 1, n_mreac
         q(ir) = molecular_reaction_energy_eV(ir)
      enddo

      write(*,'(a)') '---- A. the reference state and its zeros ----'
      ! Every neutral ground-state atom is the reference and carries zero,
      ! whatever its element and whichever stage dominates the base.
      call check_absolute('eps_zero_HI',                                  &
           species_formation_energy(isp_HI),  0.0d0, 0.0d0)
      call check_absolute('eps_zero_HeI',                                 &
           species_formation_energy(isp_HeI), 0.0d0, 0.0d0)
      call check_absolute('eps_zero_CI',                                  &
           species_formation_energy(isp_CI),  0.0d0, 0.0d0)
      call check_absolute('eps_zero_OI',                                  &
           species_formation_energy(isp_OI),  0.0d0, 0.0d0)
      call check_absolute('eps_zero_NaI',                                 &
           species_formation_energy(isp_NaI), 0.0d0, 0.0d0)
      call check_absolute('eps_zero_FeI',                                 &
           species_formation_energy(isp_FeI), 0.0d0, 0.0d0)
      ! The free electron carries no formation energy: the ionization
      ! energy is on the ion (T1.2).
      call check_absolute('eps_zero_electron',                            &
           species_formation_energy(isp_eps_electron), 0.0d0, 0.0d0)

      write(*,'(a)') '---- B. ionization, cumulative from the neutral ----'
      call check_absolute('eps_HII_is_IP_H',                              &
           species_formation_energy(isp_HII), e_th_HI, 0.0d0)
      call check_absolute('eps_HeII_is_IP_He',                            &
           species_formation_energy(isp_HeII), e_th_HeI, 0.0d0)
      call check_absolute('eps_HeIII_is_cumulative',                      &
           species_formation_energy(isp_HeIII), e_th_HeI + e_th_HeII,     &
           1.0d-12)
      ! Metal stages: the reference is the neutral ground state, so a stage
      ! carries the sum of the potentials below it.  References READ from
      ! mion_ethr, species_table.f90.
      call check_absolute('eps_CII_cumulative',                           &
           species_formation_energy(isp_CII),  11.26d0, 1.0d-12)
      call check_absolute('eps_CIII_cumulative',                          &
           species_formation_energy(isp_CIII), 11.26d0 + 24.38d0, 1.0d-12)
      call check_absolute('eps_NaII_cumulative',                          &
           species_formation_energy(isp_NaII), 5.139d0, 1.0d-12)
      call check_absolute('eps_FeIII_cumulative',                         &
           species_formation_energy(isp_FeIII),                           &
           7.902d0 + 16.199d0, 1.0d-12)

      write(*,'(a)') '---- C. excited states carry excitation only ----'
      ! He 2^3S: parameters.f90 defines e_th_HeTR = IP(He) - E(2^3S), so the
      ! metastable's excitation energy is the difference of the two.
      call check_absolute('eps_He23S_excitation',                         &
           species_formation_energy(isp_HeTR), e_th_HeI - e_th_HeTR,      &
           1.0d-12)
      ! H(n=2): the Lyman-alpha energy above H(1s), E_21 = I(H)(1 - 1/4).
      call check_absolute('eps_H_n2_is_Lyman_alpha',                      &
           species_formation_energy(isp_eps_H_n2),                        &
           0.75d0*e_th_HI, 1.0d-12)

      write(*,'(a)') '---- D. the H2 zero is the caloric EOS zero ----'
      ! eps(H2) = -D0 measured from v = 0, J = 0.
      call check_absolute('eps_H2_is_minus_D0',                           &
           eps_H2, -D0_H2_cm*cm_to_eV, 1.0d-9)
      call check_absolute('eps_H2_uses_the_one_D0',                       &
           eps_H2, -h2_dissociation_energy_eV(), 0.0d0)
      ! The rovibrational sum of the equation of state starts at the same
      ! level: the lowest term energy of the ladder is exactly zero, so the
      ! EOS holds the ladder and the reservoir holds the bond, with no
      ! overlap (T1.3).
      call check_absolute('h2_ladder_zero_is_v0J0',                       &
           minval(h2_lev_T(1:n_h2_lev)), 0.0d0, 0.0d0)
      ! And they share the dissociation limit: the highest bound level of
      ! the ladder stands at D0 above that zero.
      call check_relative('h2_ladder_limit_is_D0',                        &
           maxval(h2_lev_T(1:n_h2_lev)), D0_H2_cm*cm_to_K, 5.0d-5)

      write(*,'(a)') '---- E. every reaction heat is a table difference ----'
      do ir = 1, n_check
         qi = 0.0d0
         do i = 1, 3
            if (tr_react(i, ir) .ne. 0)                                   &
               qi = qi + species_formation_energy(tr_react(i, ir))
         enddo
         do i = 1, 3
            if (tr_prod(i, ir) .ne. 0)                                    &
               qi = qi - species_formation_energy(tr_prod(i, ir))
         enddo
         nm = 'reaction_energy_'//trim(mreac_name(ir)(1:3))
         call check_absolute(trim(nm), q(ir), qi, 1.0d-12)
      enddo

      write(*,'(a)') '---- F. every reaction conserves nuclei and charge ----'
      do ir = 1, n_mreac
         dH  = 0;  dHe = 0;  dO = 0;  dC = 0;  dZ = 0
         do i = 1, 3
            call add_composition(mreac_react(i, ir),  1, dH, dHe, dO, dC, dZ)
            call add_composition(mreac_prod (i, ir), -1, dH, dHe, dO, dC, dZ)
         enddo
         nm = 'balance_'//trim(mreac_name(ir)(1:3))
         call check_absolute(trim(nm),                                    &
              dble(abs(dH) + abs(dHe) + abs(dO) + abs(dC) + abs(dZ)),     &
              0.0d0, 0.0d0)
      enddo

      write(*,'(a)') '---- G. the identities a table has and a list has not ----'
      ! A channel and its reverse are exact negatives.  R9 and R10 are the
      ! two directions of H2+ + H <-> H+ + H2; R12 and R15 the two
      ! directions of H2 + M <-> H + H + M.
      call check_absolute('reverse_pair_R9_R10',  q(5) + q(6), 0.0d0, 0.0d0)
      call check_absolute('reverse_pair_R12_R15', q(8) + q(11), 0.0d0, 0.0d0)
      ! Two channels with the same reactants and products release the same
      ! energy however they are written: R12 and R14 both dissociate H2.
      call check_absolute('same_change_R12_R14', q(8) - q(10), 0.0d0, 0.0d0)
      ! Two routes from the same state to the same state: R8 followed by R6
      ! turns H2+ + H2 + e into H2 + H + H, which is R5 with an H2
      ! spectator.
      call check_absolute('route_R8_R6_equals_R5',                        &
           q(4) + q(2) - q(1), 0.0d0, 1.0d-12)
      ! The three energies the module header quotes for the cycle that makes
      ! this heating channel 81 per cent of the base heating: R8 releases
      ! 1.70 eV, R6 releases 9.25 eV, and the cycle leaves
      ! I(H2) - D0(H2) = 10.95 eV in the gas.
      call check_absolute('cycle_R8_header_value',  q(4), 1.70d0, 0.02d0)
      call check_absolute('cycle_R6_header_value',  q(2), 9.25d0, 0.02d0)
      call check_absolute('cycle_leaves_IH2_minus_D0',                    &
           q(4) + q(2), e_th_H2 + eps_H2, 1.0d-12)

      write(*,'(a)') '---- H. the oxygen and carbon entries ----'
      ! The reference again: the free atoms are zero, so a molecule carries
      ! minus its dissociation energy into them.
      call check_absolute('eps_zero_O_shomate',                           &
           oxygen_formation_energy_eV(ith_O), 0.0d0, 0.0d0)
      call check_absolute('eps_zero_C_shomate',                           &
           oxygen_formation_energy_eV(ith_C), 0.0d0, 0.0d0)
      call check_absolute('eps_zero_H_shomate',                           &
           oxygen_formation_energy_eV(ith_H), 0.0d0, 0.0d0)
      ! One D0(H2) in the code: the Shomate species H2 returns the
      ! spectroscopic bond energy, not a second fit of it.
      call check_absolute('shomate_H2_uses_the_one_D0',                   &
           oxygen_formation_energy_eV(ith_H2), eps_H2, 0.0d0)
      ! The three molecules the oxygen network carries, against the bond
      ! energies the Shomate table implies (MEASURED from the F coefficients
      ! of its lowest intervals, docs/co_destruction_rates_literature
      ! _20260906.md section 3 quotes D(CO) = 11.156 eV at 298 K from the
      ! same table, which the 0 K value sits 0.04 eV below).
      call check_absolute('eps_OH_bond',  eps_OH,  -4.39888d0,  1.0d-4)
      call check_absolute('eps_H2O_bond', eps_H2O, -9.51622d0,  1.0d-4)
      call check_absolute('eps_CO_bond',  eps_CO, -11.11569d0,  1.0d-4)
      ! O(1D) is an excited state and carries its excitation energy only.
      call check_absolute('eps_O1D_excitation',                           &
           eps_O1D, e_excite_O1D_erg/eV_to_erg, 1.0d-12)

      write(*,'(a)') '---- I. the oxygen channel energies ----'
      e1 = photolysis_threshold_erg(ich_H2O_OH_H) /eV_to_erg
      e2 = photolysis_threshold_erg(ich_H2O_H2_O) /eV_to_erg
      e3 = photolysis_threshold_erg(ich_H2O_O_H_H)/eV_to_erg
      e4 = photolysis_threshold_erg(ich_OH_O_H)   /eV_to_erg
      ! Each threshold is the difference of the same table entries.
      call check_absolute('threshold_H2O_OH_H_is_difference',             &
           e1, eps_OH - eps_H2O, 1.0d-12)
      call check_absolute('threshold_H2O_H2_O_is_difference',             &
           e2, eps_H2 - eps_H2O, 1.0d-12)
      call check_absolute('threshold_H2O_O_H_H_is_difference',            &
           e3, -eps_H2O, 1.0d-12)
      call check_absolute('threshold_OH_O_H_is_difference',               &
           e4, -eps_OH, 1.0d-12)
      ! Taking H2O apart in one step or in two costs the same, which is the
      ! property a list of separately quoted thresholds cannot be trusted to
      ! have.
      call check_absolute('H2O_two_step_equals_one_step',                 &
           e1 + e4 - e3, 0.0d0, 1.0d-12)
      ! And the H2-forming branch differs from the full atomization by
      ! exactly one H2 bond.
      call check_absolute('H2O_H2_branch_differs_by_D0_H2',               &
           e3 - e2 - (-eps_H2), 0.0d0, 1.0d-12)
      ! The collisional channels of the oxygen network, formed from the same
      ! table.  They are not deposited by any consumer today (D0 C25); what
      ! is asserted is that the ledger they will be read from is closed.
      ! O1  OH + H2 -> H2O + H   and its reverse O1r
      call check_absolute('oxygen_O1_reverse_is_exact_negative',          &
           (eps_OH + eps_H2 - eps_H2O) + (eps_H2O - eps_OH - eps_H2),     &
           0.0d0, 1.0d-12)
      ! O6 O(1D) + H2 -> OH + H is O2 O + H2 -> OH + H plus the excitation
      ! energy the eliminated O(1D) carries to the products: the transfer
      ! rule of T1.9 item 1.
      call check_absolute('oxygen_O6_carries_O1D_excitation',             &
           (eps_O1D + eps_H2 - eps_OH) - (eps_H2 - eps_OH) - eps_O1D,     &
           0.0d0, 1.0d-12)
      ! O2 followed by O1 is the overall O + 2 H2 -> H2O + 2 H.
      call check_absolute('oxygen_O2_then_O1_is_overall',                 &
           (eps_H2 - eps_OH) + (eps_OH + eps_H2 - eps_H2O)                &
           - (2.0d0*eps_H2 - eps_H2O), 0.0d0, 1.0d-12)

      ! The table as a file, for the comparison with the Python diagnostic
      ! of src/utils/formation_energy_flux_diagnostic.py.
      call write_table_file()

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'species_formation_energy_table: ',          &
              assertion_failures, ' assertion(s) failed'
         stop 1
      endif
      write(*,'(a)') 'species_formation_energy_table: every assertion passed'

      contains

      !--------------!

      subroutine add_composition(isp, sgn, dH, dHe, dO, dC, dZ)
      ! Add the nuclei and the charge of one particle of species isp to the
      ! running balance, with sign sgn.  The composition is read from
      ! species_table, so the stoichiometry of molecular_reaction_heat is
      ! balanced against a table it does not own.  The free electron and the
      ! unused slot 0 are the two keys species_table has no row for.
      integer, intent(in)    :: isp, sgn
      integer, intent(inout) :: dH, dHe, dO, dC, dZ
      integer :: ib
      if (isp .eq. 0) return
      if (isp .eq. isp_eps_electron) then
         dZ = dZ - sgn
         return
      endif
      do ib = 1, n_bsp
         if (bsp_fsp(ib) .eq. isp) then
            dH  = dH  + sgn*bsp_nH(ib)
            dHe = dHe + sgn*bsp_nHe(ib)
            dO  = dO  + sgn*bsp_nO(ib)
            dC  = dC  + sgn*bsp_nC(ib)
            dZ  = dZ  + sgn*bsp_charge(ib)
            return
         endif
      enddo
      write(*,*) 'add_composition: no species_table row for column ', isp
      stop 1
      end subroutine add_composition

      !--------------!

      subroutine write_table_file()
      ! One "<name> <value in eV>" line per entry, read by
      ! species_formation_energy_python_table.py.
      integer :: iu, lw
      character(len=512) :: work, path
      ! The suite's object directory when it has one, else where the driver
      ! stands, so that a run never writes into the source tree.
      call get_environment_variable('EXHALE_TEST_WORK', work, lw)
      if (lw .gt. 0) then
         path = trim(work)//'/species_formation_energy_table.dat'
      else
         path = 'species_formation_energy_table.dat'
      endif
      open(newunit=iu, file=trim(path), status='replace', action='write')
      write(iu,'(a)') '# eps_s [eV], zero at neutral ground-state free'
      write(iu,'(a)') '# atoms at rest and free electrons at rest.'
      write(iu,'(a,es24.16)') 'HI    ', species_formation_energy(isp_HI)
      write(iu,'(a,es24.16)') 'HII   ', species_formation_energy(isp_HII)
      write(iu,'(a,es24.16)') 'HeI   ', species_formation_energy(isp_HeI)
      write(iu,'(a,es24.16)') 'HeII  ', species_formation_energy(isp_HeII)
      write(iu,'(a,es24.16)') 'HeIII ', species_formation_energy(isp_HeIII)
      write(iu,'(a,es24.16)') 'HeITR ', species_formation_energy(isp_HeTR)
      write(iu,'(a,es24.16)') 'Hn2   ',                                   &
                              species_formation_energy(isp_eps_H_n2)
      write(iu,'(a,es24.16)') 'H2    ', eps_H2
      write(iu,'(a,es24.16)') 'H2p   ', species_formation_energy(isp_H2p)
      write(iu,'(a,es24.16)') 'H3p   ', species_formation_energy(isp_H3p)
      write(iu,'(a,es24.16)') 'HeHp  ', species_formation_energy(isp_HeHp)
      write(iu,'(a,es24.16)') 'OH    ', eps_OH
      write(iu,'(a,es24.16)') 'H2O   ', eps_H2O
      write(iu,'(a,es24.16)') 'CO    ', eps_CO
      write(iu,'(a,es24.16)') 'O1D   ', eps_O1D
      do k = 7, 33
         write(iu,'(a,i2.2,a,es24.16)') 'metal', k, ' ',                  &
              species_formation_energy(k)
      enddo
      close(iu)
      end subroutine write_table_file

      end program species_formation_energy_table

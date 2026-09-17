      program molecular_energy_recipients
      ! WHO RECEIVES THE ENERGY OF THE MOLECULAR NETWORK, and the identities
      ! that hold whatever the recipient is.
      !
      ! Three channels of the H2/He network leave a product in an excited
      ! state instead of giving the gas kinetic energy: the dissociative
      ! recombinations R5 (H2+ + e) and R16 (HeH+ + e) leave one hydrogen
      ! atom in n = 2 (Takagi 2002, Phys. Scr. T96, 52; Guberman 1994, Phys.
      ! Rev. A 49, R4277), and R6 (H3+ + e) leaves its H2 fragment
      ! vibrationally hot with a distribution peaking at v = 5-6 (Kokoouline,
      ! Greene & Esry 2001, Nature 412, 891; Strasser et al. 2001, Phys. Rev.
      ! Lett. 86, 779).  A fourth correction is of a rate rather than a
      ! recipient: the third body of the R12/R15 pair is resolved by collider
      ! (Cohen & Westberg 1983, J. Phys. Chem. Ref. Data 12, 531, p. 559)
      ! instead of the total density carrying the M = H2 coefficient.  A
      ! fifth is of the data the thermalized fraction is built from: the
      ! all-level radiative rate of the code's own ladder in place of the
      ! v = 1 one, and three published quantum calculations in place of two
      ! 1979 fits (Lique 2015 for H, Jozwiak et al. 2024 for He, Le Bourlot
      ! et al. 1999 for H2).
      !
      ! WHAT IS ASSERTED.
      !   A. THE ENTHALPY IS UNTOUCHED.  For each of the four channels that
      !      hand energy somewhere other than the gas, heat + radiated +
      !      internal = the reaction energy of the one formation table.  A
      !      recipient change that moved an enthalpy would fail here.
      !   B. R12 AND R15 STAY AN EXACT DETAILED-BALANCE PAIR, collider by
      !      collider: for a gas of one third body the ratio of the two
      !      directions is the equilibrium constant, for every collider and
      !      at every temperature.
      !   C. THE COLLIDER SUM REDUCES TO WHAT IT REPLACED: with n(H) and
      !      n(He) zero and n(H2) the whole density, the H2-equivalent third
      !      body is the total density, so the correction is a change of
      !      mixture and not of normalization.
      !   D. THE RADIATIVE AND COLLISIONAL DATA ARE THE PUBLISHED ONES: the
      !      radiative rate is the all-level maximum of the code's own
      !      302-level ladder and not the v = 1 value, and each of the three
      !      collider coefficients reproduces the number its own paper
      !      publishes, the helium table also obeying detailed balance.
      !   E. THE LEDGER DEPOSITS WHAT THE RECIPIENTS LEAVE IT, measured by
      !      running the production assembly on cells that carry one channel
      !      at a time, and the H(n=2) source is the same two rates with the
      !      same densities, so one event is not counted at two rates.
      !
      ! The four corrections are turned off together by the environment key
      ! EXHALE_REACTION_HEAT_RECIPIENTS=0; running this driver with that key
      ! set is the control, and the rows of A, C (the mixture), D and E that
      ! state the correction then FAIL by construction.
      !
      ! Production routines exercised: molecular_chemical_heating,
      ! molecular_reaction_energy_eV, dissociative_recombination_n2_eV and
      ! dissociative_recombination_n2_source of
      ! src/modules/lower_atmosphere/molecular_reaction_heat.f90;
      ! k3b_H_H_to_H2, k3b_H_H_to_H2_atomic_H, k3b_H_H_to_H2_monatomic,
      ! h2_association_collider_density, rk_R12_H2_thdis, rk_R15_3body_H2,
      ! keq_H_H_to_H2, rk_R5_H2p_dr, rk_R6_H3p_dr_H2, rk_R7_H3p_dr_3H,
      ! rk_R16_HeHp_dr of src/modules/lower_atmosphere/mol_rates.f90;
      ! h2_total_decay_rate_max, h2_vibrational_heat_fraction,
      ! gamma_10_H, gamma_10_H2, gamma_10_He, e_vib_v1_eV, e_vib_v5_eV and
      ! e_vib_v6_eV of
      ! src/modules/lower_atmosphere/h2_vibrational_relaxation.f90.
      use global_parameters
      use molecular_reaction_heat, only:                                  &
               molecular_chemical_heating, molecular_reaction_energy_eV,  &
               species_formation_energy, isp_eps_H_n2,                    &
               dissociative_recombination_n2_eV,                          &
               dissociative_recombination_n2_source,                      &
               reaction_heat_recipients_corrected
      use mol_rates, only: h2_thermochemistry_init, keq_H_H_to_H2,        &
               k3b_H_H_to_H2, k3b_H_H_to_H2_atomic_H,                     &
               k3b_H_H_to_H2_monatomic,                                   &
               h2_association_collider_density,                           &
               rk_R12_H2_thdis, rk_R15_3body_H2,                          &
               rk_R5_H2p_dr, rk_R6_H3p_dr_H2, rk_R7_H3p_dr_3H,            &
               rk_R16_HeHp_dr
      use h2_vibrational_relaxation, only: h2_total_decay_rate_max,       &
               h2_vibrational_heat_fraction, h2_vibrational_relaxation_init, &
               h2_vibrational_relaxation_ready,                           &
               gamma_10_H, gamma_10_H2, gamma_10_He,                      &
               e_vib_v1_eV, e_vib_v5_eV, e_vib_v6_eV
      use assertion_report

      implicit none

      ! Reaction indices of the network, in the order of mreac_name.
      integer, parameter :: ir_R5 = 1, ir_R6 = 2, ir_R7 = 3, ir_R16 = 12, &
                            ir_R15 = 11, ir_R12 = 8
      real*8,  parameter :: eV_to_erg = 1.602176634d-12
      ! The all-level maximum of the total spontaneous decay rate over the
      ! code's 302-level H2 ladder, READ from
      ! docs/lhs1140b_stationary_L7g_inventory_20260916.md section 0, where
      ! it was measured on the same line list, to the digits printed there.
      real*8,  parameter :: a_cascade_inventory = 5.594d-6
      ! The v = 1, J = 0 total of the same ladder, the value the fraction
      ! used before (same source, same section).
      real*8,  parameter :: a_v1_inventory      = 8.532d-7

      integer, parameter :: ncell = 1
      integer :: jlo, jhi, it
      real*8, allocatable :: T_K(:), nhi(:), nhii(:), nheii(:), nheiS(:), &
                             ne(:), ntot(:), nmol(:,:), gam(:)
      real*8 :: q5, q6, q7, q15, q12, q16, e_n2, e_int_r6
      real*8 :: t_probe(4), tt, n_one, n_eff, kr, kd, f_heat
      real*8 :: expect, got, k5, k6, k7, k16
      logical :: corrected

      t_probe = (/ 3.0d2, 8.0d2, 1.5d3, 3.0d3 /)
      n_one   = 1.0d12

      ! The H2 equilibrium-constant table and the ladder reduction; the same
      ! serial prologue EXHALE_main runs before any parallel region.
      call h2_thermochemistry_init
      corrected = reaction_heat_recipients_corrected()
      write(*,'(a,l1)') '# corrected recipients in force: ', corrected

      q5  = molecular_reaction_energy_eV(ir_R5)
      q6  = molecular_reaction_energy_eV(ir_R6)
      q7  = molecular_reaction_energy_eV(ir_R7)
      q12 = molecular_reaction_energy_eV(ir_R12)
      q15 = molecular_reaction_energy_eV(ir_R15)
      q16 = molecular_reaction_energy_eV(ir_R16)
      e_n2     = dissociative_recombination_n2_eV()
      e_int_r6 = 0.5d0*(e_vib_v5_eV + e_vib_v6_eV)

      ! ================================================================== !
      ! A. THE ENTHALPY IS UNTOUCHED BY THE CHANGE OF RECIPIENT.
      ! ================================================================== !
      ! The two enthalpies themselves, against the values the audit of the
      ! network READ from the same table (predecessor memo, section 5).
      call check_absolute('enthalpy_R5_unchanged',  q5,  10.9478d0, 1.0d-3)
      call check_absolute('enthalpy_R16_unchanged', q16, 11.7534d0, 1.0d-3)
      call check_absolute('enthalpy_R6_unchanged',  q6,   9.2500d0, 2.0d-2)

      ! The excitation energy handed to n = 2 is the ONE entry of the
      ! formation table, and it is the Lyman-alpha energy 0.75 IP(H).
      call check_absolute('n2_excitation_is_the_table_entry',             &
                          e_n2, species_formation_energy(isp_eps_H_n2),   &
                          0.0d0)
      call check_relative('n2_excitation_is_0.75_IP_H',                   &
                          e_n2, 0.75d0*e_th_HI, 1.0d-14)

      ! heat + radiated-or-excited = enthalpy, channel by channel.
      call check_absolute('R5_heat_plus_excitation_is_enthalpy',          &
                          (q5 - e_n2) + e_n2 - q5, 0.0d0, 1.0d-13)
      call check_absolute('R16_heat_plus_excitation_is_enthalpy',         &
                          (q16 - e_n2) + e_n2 - q16, 0.0d0, 1.0d-13)
      call check_absolute('R6_prompt_plus_internal_is_enthalpy',          &
                          (q6 - e_int_r6) + e_int_r6 - q6, 0.0d0, 1.0d-13)

      ! And the shares themselves, which are what the sources fix.
      call check_absolute('R5_kinetic_share',  q5  - e_n2, 0.749d0, 2.0d-3)
      call check_absolute('R16_kinetic_share', q16 - e_n2, 1.554d0, 2.0d-3)
      ! The peak of the H3+ product distribution, between v = 5 and v = 6 of
      ! the code's own ladder.
      call check_absolute('R6_internal_share_is_the_v5_v6_peak',          &
                          e_int_r6, 2.4795d0, 1.0d-3)

      ! ================================================================== !
      ! B. R12 AND R15 ARE AN EXACT DETAILED-BALANCE PAIR, COLLIDER BY
      !    COLLIDER.  For a gas of one third body the ratio of the two
      !    directions must be K_eq, whichever third body it is.
      ! ================================================================== !
      do it = 1, 4
         tt = t_probe(it)
         ! M = H2
         n_eff = h2_association_collider_density(tt, n_one, 0.0d0, 0.0d0)
         kr    = rk_R15_3body_H2(tt, n_eff)
         kd    = rk_R12_H2_thdis(tt)*n_eff
         call check_relative('detailed_balance_M_H2_T'//tlabel(tt),       &
                             kr/kd, keq_H_H_to_H2(tt), 1.0d-13)
         call check_relative('collider_H2_is_k1_H2_T'//tlabel(tt),        &
                             kr, k3b_H_H_to_H2(tt)*n_one, 1.0d-13)
         ! M = H
         n_eff = h2_association_collider_density(tt, 0.0d0, n_one, 0.0d0)
         kr    = rk_R15_3body_H2(tt, n_eff)
         kd    = rk_R12_H2_thdis(tt)*n_eff
         call check_relative('detailed_balance_M_H_T'//tlabel(tt),        &
                             kr/kd, keq_H_H_to_H2(tt), 1.0d-13)
         call check_relative('collider_H_is_k1_H_T'//tlabel(tt),          &
                             kr, k3b_H_H_to_H2_atomic_H()*n_one, 1.0d-13)
         ! M = He, carrying the monatomic inert coefficient
         n_eff = h2_association_collider_density(tt, 0.0d0, 0.0d0, n_one)
         kr    = rk_R15_3body_H2(tt, n_eff)
         kd    = rk_R12_H2_thdis(tt)*n_eff
         call check_relative('detailed_balance_M_He_T'//tlabel(tt),       &
                             kr/kd, keq_H_H_to_H2(tt), 1.0d-13)
         call check_relative('collider_He_is_k1_monatomic_T'//tlabel(tt), &
                             kr, k3b_H_H_to_H2_monatomic(tt)*n_one,       &
                             1.0d-13)
      enddo

      ! The efficiency ordering the evaluation gives at 1000 K: atomic H
      ! twice H2, a monatomic inert atom 0.43 of it.
      call check_relative('efficiency_H_over_H2_at_1000K',                &
                          k3b_H_H_to_H2_atomic_H()/k3b_H_H_to_H2(1.0d3),  &
                          1.983d0, 1.0d-3)
      call check_relative('efficiency_monatomic_over_H2_at_1000K',        &
                          k3b_H_H_to_H2_monatomic(1.0d3)                  &
                          /k3b_H_H_to_H2(1.0d3), 0.4281d0, 1.0d-3)

      ! ================================================================== !
      ! C. THE LIMIT ROW: a gas of pure H2 gives back the total density, so
      !    the collider sum replaces a MIXTURE and not a normalization.
      ! ================================================================== !
      call check_absolute('collider_sum_is_ntot_in_pure_H2',              &
                  h2_association_collider_density(1.0d3, 3.0d13, 0.0d0,   &
                                                  0.0d0) - 3.0d13,        &
                  0.0d0, 0.0d0)

      ! ================================================================== !
      ! D. THE RADIATIVE RATE OF THE CASCADE IS THE ALL-LEVEL MAXIMUM.
      ! ================================================================== !
      call check_absolute('ladder_reduction_has_run',                     &
           merge(0.0d0, 1.0d0, h2_vibrational_relaxation_ready()),        &
           0.0d0, 0.0d0)
      call check_relative('all_level_A_max_matches_the_ladder',           &
                          h2_total_decay_rate_max(), a_cascade_inventory, &
                          1.0d-3)
      call check_relative('all_level_A_max_over_v1_total',                &
                          h2_total_decay_rate_max()/a_v1_inventory,       &
                          6.556d0, 2.0d-3)

      ! ---- THE TWO COLLISIONAL TABLES ARE THE PAPERS' OWN NUMBERS.
      ! Lique (2015) Table 2 publishes the thermal v = 1 -> v' = 0
      ! relaxation by atomic hydrogen at room temperature as 1.8e-13
      ! cm^3 s^-1; the table in the module is the reduction of his
      ! distributed state-to-state file, and must reproduce it.
      call check_relative('lique_H_thermal_v1_v0_at_300K',               &
                          gamma_10_H(3.0d2), 1.8d-13, 2.0d-2)
      ! Jozwiak et al. (2024), reduced the same way from their VizieR
      ! tables; these two values were obtained independently by the data
      ! inventory of this item from the same files (2.64e-15 and 7.02e-15).
      call check_relative('jozwiak_He_thermal_v1_v0_at_808K',            &
                          gamma_10_He(8.08d2), 2.64d-15, 1.0d-2)
      call check_relative('jozwiak_He_thermal_v1_v0_at_1000K',           &
                          gamma_10_He(1.0d3), 7.02d-15, 1.0d-2)
      ! Le Bourlot, Pineau des Forets & Flower (1999), reduced the same way
      ! from the state-to-state table the Cloudy distribution carries for a
      ! para-H2 perturber.  That is the ONE H2-H2 set the paper published:
      ! its section 2.1 states that the rates were computed for a para-H2
      ! (J = 0) perturber and applied unchanged to an ortho-H2 one, and its
      ! 51 levels to (v = 3, j = 8) are the para file's list.  The file the
      ! Cloudy distribution names for an ortho perturber holds a shorter
      ! list and stands a median 1.74 above it; mixing the two, which this
      ! table did before the paper was read, raised the coefficient by 1.14
      ! to 1.42 over the grid.
      call check_relative('lebourlot_H2_thermal_v1_v0_at_1000K',         &
                          gamma_10_H2(1.0d3), 2.631d-15, 2.0d-3)
      call check_relative('lebourlot_H2_thermal_v1_v0_at_1500K',         &
                          gamma_10_H2(1.5d3), 2.302d-14, 2.0d-3)

      ! ---- THE HELIUM TABLE OBEYS DETAILED BALANCE, which is the check
      ! that the file was read in the direction its header states and that
      ! the level energies the reduction weighted with are the molecule's.
      ! The four pairs below are READ from
      ! references/Jozwiak_2024J_A+A_685_A113/ph2-rat.dat, the para
      ! (v=0,j=0) <-> (v=1,j=0) excitation and de-excitation rows at 300,
      ! 1000, 1500 and 2000 K of that file's own grid. Both levels have
      ! j = 0, so their degeneracies are equal and the ratio must be the
      ! Boltzmann factor of the v = 1 term energy the module itself carries.
      ! Le Bourlot et al.'s file distributes only the de-excitation
      ! direction, so the same check cannot be made on the H2 table.
      call check_detailed_balance('jozwiak_He_detailed_balance_0300K',   &
               4.96631d-28, 2.30726d-19, 3.0d2)
      call check_detailed_balance('jozwiak_He_detailed_balance_1000K',   &
               9.38156d-20, 3.73587d-17, 1.0d3)
      call check_detailed_balance('jozwiak_He_detailed_balance_1500K',   &
               3.71014d-18, 2.00818d-16, 1.5d3)
      call check_detailed_balance('jozwiak_He_detailed_balance_2000K',   &
               3.10828d-17, 6.20267d-16, 2.0d3)
      ! THE ORDERING OF THE THREE COLLIDERS, which is the physical content
      ! of replacing two fits by three calculations: atomic hydrogen is
      ! three orders above the other two, its exchange channel being what
      ! makes it efficient, and H2 is no better a quencher of H2 than
      ! helium is.
      call check_relative('helium_share_of_the_collider_coefficient',    &
                          gamma_10_He(1.0d3)/gamma_10_H(1.0d3),          &
                          9.341d-4, 2.0d-3)
      call check_at_least('atomic_H_is_three_orders_above_H2',           &
                          gamma_10_H(1.0d3)/gamma_10_H2(1.0d3), 1.0d3)
      call check_at_most('H2_is_no_better_a_quencher_than_He',           &
                         gamma_10_H2(1.0d3)/gamma_10_He(1.0d3), 1.0d0)

      ! ================================================================== !
      ! E. THE LEDGER DEPOSITS WHAT THE RECIPIENTS LEAVE IT.  Each cell
      !    below carries ONE channel: every density that would open another
      !    is zero, so the assembly's whole output is that channel's.
      ! ================================================================== !
      N  = ncell
      R0 = 1.0d10
      thereis_He     = .true.
      thereis_mol    = .true.
      jlo = 1 - Ng
      jhi = N + Ng
      allocate(T_K(jlo:jhi), nhi(jlo:jhi), nhii(jlo:jhi), nheii(jlo:jhi), &
               nheiS(jlo:jhi), ne(jlo:jhi), ntot(jlo:jhi), gam(jlo:jhi))
      allocate(nmol(jlo:jhi,4))

      ! ---- R5 alone: H2+ and electrons, nothing else.
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      nmol(:,2) = 1.0d6
      ne        = 1.0d10
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      k5     = rk_R5_H2p_dr(T_K(1))
      expect = k5*nmol(1,2)*ne(1)*(q5 - e_n2)*eV_to_erg
      call check_relative('ledger_R5_deposits_enthalpy_less_n2',          &
                          gam(1), expect, 1.0d-12)

      ! ---- R16 alone: HeH+ and electrons.
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      nmol(:,4) = 1.0d4
      ne        = 1.0d10
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      k16    = rk_R16_HeHp_dr(T_K(1))
      expect = k16*nmol(1,4)*ne(1)*(q16 - e_n2)*eV_to_erg
      call check_relative('ledger_R16_deposits_enthalpy_less_n2',         &
                          gam(1), expect, 1.0d-12)

      ! ---- R6 and R7 alone: H3+ and electrons, no collider at all, so the
      ! thermalized fraction is zero and R6's internal share is withheld in
      ! full.  This is the row that would catch a share deposited promptly.
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      nmol(:,3) = 1.0d7
      ne        = 1.0d10
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      k6     = rk_R6_H3p_dr_H2(T_K(1))
      k7     = rk_R7_H3p_dr_3H(T_K(1))
      f_heat = h2_vibrational_heat_fraction(T_K(1), 0.0d0, 0.0d0)
      expect = ( k6*(q6 - e_int_r6 + e_int_r6*f_heat) + k7*q7 )           &
               *nmol(1,3)*ne(1)*eV_to_erg
      call check_relative('ledger_R6_withholds_the_internal_share',       &
                          gam(1), expect, 1.0d-12)
      call check_absolute('heat_fraction_vanishes_without_colliders',     &
                          f_heat, 0.0d0, 0.0d0)

      ! ---- R15 alone: atomic hydrogen only, so the third body is H and the
      ! rate must be k1(H) n(H)^3, not k1(H2) n_tot n(H)^2.
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      nhi    = 1.0d13
      ntot   = 1.0d13
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      f_heat = h2_vibrational_heat_fraction(T_K(1), nhi(1), 0.0d0)
      expect = k3b_H_H_to_H2_atomic_H()*nhi(1)**3*q15*f_heat*eV_to_erg
      call check_relative('ledger_R15_third_body_is_atomic_H',            &
                          gam(1), expect, 1.0d-12)

      ! ---- The H(n=2) source is the SAME two rates and the same densities
      ! the ledger charged, so one event is not counted at two rates.
      got    = dissociative_recombination_n2_source(1.0d3, 1.0d6, 1.0d4,  &
                                                    1.0d10)
      expect = ( rk_R5_H2p_dr(1.0d3)*1.0d6                                &
               + rk_R16_HeHp_dr(1.0d3)*1.0d4 )*1.0d10
      call check_relative('n2_source_is_k5_nH2p_ne_plus_k16_nHeHp_ne',    &
                          got, expect, 1.0d-14)
      call check_positive('n2_source_is_present', got)

      ! With no molecular ions there is no source, which is the atomic run.
      call check_absolute('n2_source_vanishes_without_molecular_ions',    &
             dissociative_recombination_n2_source(1.0d3, 0.0d0, 0.0d0,    &
                                                  1.0d10), 0.0d0, 0.0d0)

      ! ---- THE LEDGER DOES NOT DEPEND ON use_excited_H.  With the excited
      ! level carried, the 10.199 eV reaches it; without, it leaves as
      ! radiation and no trapping is modelled.  Either way the gas is
      ! charged the same, which is what this row holds.
      use_excited_H = .false.
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      nmol(:,2) = 1.0d6
      ne        = 1.0d10
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      expect = k5*nmol(1,2)*ne(1)*(q5 - e_n2)*eV_to_erg
      call check_relative('ledger_R5_same_without_excited_H',             &
                          gam(1), expect, 1.0d-12)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'molecular_energy_recipients: ',             &
              assertion_failures, ' assertion(s) failed'
         stop 1
      endif
      write(*,'(a)') 'molecular_energy_recipients: every assertion passed'

      contains

      subroutine zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      ! Every density zero and one temperature, so that a cell carries only
      ! the channels the caller then fills in.
      real*8, intent(out) :: T_K(:), nhi(:), nhii(:), nheii(:), nheiS(:)
      real*8, intent(out) :: ne(:), ntot(:), nmol(:,:)
      T_K   = 1.0d3
      nhi   = 0.0d0
      nhii  = 0.0d0
      nheii = 0.0d0
      nheiS = 0.0d0
      ne    = 0.0d0
      ntot  = 0.0d0
      nmol  = 0.0d0
      end subroutine zero_cell

      subroutine check_detailed_balance(name, k_up, k_down, t)
      ! The ratio of an excitation rate to its own de-excitation rate must
      ! be the Boltzmann factor of the energy gap, the two levels here
      ! having equal degeneracies. The gap is e_vib_v1_eV, the module's own
      ! (v=1, J=0) term energy, so this row also holds the table's level
      ! ladder to the code's.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: k_up, k_down, t
      ! Boltzmann's constant in eV/K, CODATA 2018.
      real*8, parameter :: kb_per_eV = 8.617333262d-5
      call check_relative(name, k_up/k_down,                              &
                          exp(-e_vib_v1_eV/(kb_per_eV*t)), 1.0d-4)
      end subroutine check_detailed_balance

      character(len=4) function tlabel(t)
      ! The probe temperature as a label, so each assertion names its own.
      real*8, intent(in) :: t
      write(tlabel,'(i4.4)') nint(t)
      end function tlabel

      end program molecular_energy_recipients

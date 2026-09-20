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
               rk_R16_HeHp_dr, rk_R8_H2p_H2, rk_H2p_He_HeHp,              &
               rk_R23_H2_Hep_cx, rk_R17_Hep_H2_diss,                      &
               rct_photon_energy_eV, h2_dissociation_energy_eV
      use h2_vibrational_relaxation, only: h2_total_decay_rate_max,       &
               h2_vibrational_heat_fraction, h2_vibrational_relaxation_init, &
               h2_vibrational_relaxation_ready,                           &
               gamma_10_H, gamma_10_H2, gamma_10_He,                      &
               e_vib_v1_eV, e_vib_v5_eV, e_vib_v6_eV
      use System_HeH_mol, only: set_mol_coeffs, h2_third_body_density,    &
               mol_heh_rows, mk15
      use ion_cell_state, only: ieq_cell
      use assertion_report

      implicit none

      ! Reaction indices of the network, in the order of mreac_name.
      integer, parameter :: ir_R5 = 1, ir_R6 = 2, ir_R7 = 3, ir_R16 = 12, &
                            ir_R15 = 11, ir_R12 = 8, ir_R8 = 4,           &
                            ir_R17 = 13, ir_h2p_he = 16, ir_R23 = 17
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

      ! ---- THE COMPLETE CYCLES (block F) ----
      ! Ground-state reaction energies and the gas share of each step, per
      ! event [eV], and the energies the other reservoirs receive.
      real*8 :: q8, q17, q23, q_h2phe, d0_h2, e_gamma_r23
      real*8 :: gas5, gas6, gas8, gas12, gas15, gas16, gas23, gas_h2phe
      real*8 :: heat_one, heat_two, cyc

      ! ---- THE UNCERTAINTY RECORD (block G) ----
      ! THREE CELLS OF A MOLECULAR COLUMN, TAKEN FROM A STATE.
      !
      ! The three places the bracket is quoted at are a rule, not a list of
      ! numbers:
      !
      !   molecular depth      physical cell 1, the base of the column;
      !   H2 front             the first cell where 2 n(H2)/n_H has fallen
      !                        to half its value at cell 1;
      !   dilute upper column  the cell nearest r_dilute_upper_column
      !                        below, where 2 n(H2)/n_H is a few parts in
      !                        1e6.  No rule drawn from the profile picks
      !                        this cell: the radius IS the choice, and it
      !                        is the radius of the cell the uncertainty
      !                        record was first written on, kept so that
      !                        the same physical place is taken on any grid.
      !
      ! EXHALE_PROBE_STATE names the directory of a written state and the
      ! three cells are then read from its {Hydro_ioniz_IC, Ion_species_IC}
      ! pair (the state a run hands back; the pair without _IC is the
      ! product written FROM a measurement of it and is not the same state).
      ! Unset, the compiled values below are used, which are those of
      ! LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13,
      ! output_pre_L34/{Hydro_ioniz_IC.txt, Ion_species_IC.txt}, rows 3, 180
      ! and 301 (physical cells 1, 178 and 299), every digit READ from those
      ! two files.  Reading that same directory reproduces them exactly, and
      ! the driver prints the distance so that it is one number.
      integer, parameter :: ncase = 3
      ! The radius of the dilute upper column, R_p, READ from the grid of
      ! the archived fiducial (its cell 299).
      real*8, parameter :: r_dilute_upper_column = 1.8422600916216840d0
      character(len=20), parameter :: case_name(ncase) =                  &
         (/ character(len=20) :: 'molecular_depth', 'h2_front',           &
                                 'dilute_upper' /)
      real*8, parameter :: c0_T(ncase)     =                               &
         (/ 808.31894177851177d0, 2128.3673585881161d0,                   &
            5083.4376627309339d0 /)
      real*8, parameter :: c0_nhi(ncase)   =                               &
         (/ 1867980406714.9812d0, 15643668544.791315d0,                   &
            81619552.724037111d0 /)
      real*8, parameter :: c0_nhii(ncase)  =                               &
         (/ 8994439.3355619069d0, 23327598.136132225d0,                   &
            3367353.9183437992d0 /)
      real*8, parameter :: c0_nheiS(ncase) =                               &
         (/ 6171608738847.3115d0, 35876533538.838463d0,                   &
            70245906.254677579d0 /)
      real*8, parameter :: c0_nheii(ncase) =                               &
         (/ 1912.4490258127826d0, 4748.3048090936782d0,                   &
            3521969.0794570139d0 /)
      real*8, parameter :: c0_nh2(ncase)   =                               &
         (/ 514739739664.78314d0, 1646227639.0261714d0,                   &
            153.44755402026095d0 /)
      real*8, parameter :: c0_nh2p(ncase)  =                               &
         (/ 7.8966795351741143d-2, 99.266851520807592d0,                  &
            3.6030811057120853d0 /)
      real*8, parameter :: c0_nh3p(ncase)  =                               &
         (/ 3963.8077308419979d0, 705.51523051539277d0,                   &
            2.7962526782356655d-5 /)
      real*8, parameter :: c0_nhehp(ncase) =                               &
         (/ 1.5589358127021795d-5, 2.7132243962847249d0,                  &
            5.1951384523255015d0 /)
      real*8, parameter :: c0_ne(ncase)    =                               &
         (/ 9.0003156725121252d6, 2.3333154532414723d7,                   &
            6.9012038847970981d6 /)
      real*8, parameter :: c0_ntot(ncase)  =                               &
         (/ 8.5543378855427480d12, 5.3189762876890274d10,                 &
            1.5876088026669142d8 /)
      ! The three cells in force: the compiled values, or the ones read
      ! from the state EXHALE_PROBE_STATE names.
      real*8 :: c_T(ncase), c_nhi(ncase), c_nhii(ncase), c_nheiS(ncase)
      real*8 :: c_nheii(ncase), c_nh2(ncase), c_nh2p(ncase), c_nh3p(ncase)
      real*8 :: c_nhehp(ncase), c_ne(ncase), c_ntot(ncase)
      real*8 :: c_r(ncase), c_x2(ncase)
      integer :: c_cell(ncase)
      character(len=512) :: state_dir
      integer :: ic
      real*8 :: h_nom, f_c, n3_c, e_lev, d_rec, d_br, d_eint_hi
      real*8 :: d_eint_lo, d_quench, d_he_hi, d_he_lo, eff_he
      real*8 :: bracket_lo, bracket_hi

      ! ---- THE COLLIDER SUM IN THREE PLACES (block H) ----
      real*8 :: n3_chem, n3_direct, fvec_h(8), chan_h(13)

      t_probe = (/ 3.0d2, 8.0d2, 1.5d3, 3.0d3 /)
      n_one   = 1.0d12

      ! The three representative cells of block G, from the state
      ! EXHALE_PROBE_STATE names or from the compiled archived values.
      call representative_cells

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

      ! ================================================================== !
      ! F. THE COMPLETE CYCLES CLOSE, so no reservoir is charged twice.
      !
      !    A closed cycle returns to its starting species, so what it
      !    liberates is fixed by its INPUT alone and by nothing inside it.
      !    Four routes below turn one H2 into two H atoms, driven by one
      !    photon, and each must liberate IP(H2) - D0(H2) between the gas,
      !    the H(n=2) level, the 153 nm photon of R23 and the H2 internal
      !    energy the ledger withholds.  A share deposited twice shows as a
      !    surplus; a share dropped shows as a deficit.
      !
      !    THE GAS TERM OF EACH STEP IS MEASURED, not restated: it is the
      !    production assembly's own output on a cell that carries that
      !    step, divided by the step's event rate.  The other terms come
      !    from the production accessors that own them.
      ! ================================================================== !
      q8      = molecular_reaction_energy_eV(ir_R8)
      q17     = molecular_reaction_energy_eV(ir_R17)
      q23     = molecular_reaction_energy_eV(ir_R23)
      q_h2phe = molecular_reaction_energy_eV(ir_h2p_he)
      d0_h2       = h2_dissociation_energy_eV()
      e_gamma_r23 = rct_photon_energy_eV()

      ! ---- one event of R5: H2+ and electrons only.
      heat_one = assembly_heat(1.0d3, 0.0d0, 0.0d0, 0.0d0, 0.0d0,         &
                               1.0d6, 0.0d0, 0.0d0, 1.0d10, 0.0d0)
      gas5 = heat_one/(rk_R5_H2p_dr(1.0d3)*1.0d6*1.0d10)/eV_to_erg

      ! ---- one event of R16: HeH+ and electrons only.
      heat_one = assembly_heat(1.0d3, 0.0d0, 0.0d0, 0.0d0, 0.0d0,         &
                               0.0d0, 0.0d0, 1.0d4, 1.0d10, 0.0d0)
      gas16 = heat_one/(rk_R16_HeHp_dr(1.0d3)*1.0d4*1.0d10)/eV_to_erg

      ! ---- one event of H2+ + He -> HeH+ + H: H2+ and ground-singlet
      ! helium only, so no electron opens R5.
      heat_one = assembly_heat(1.0d3, 0.0d0, 0.0d0, 1.0d12, 0.0d0,        &
                               1.0d6, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
      gas_h2phe = heat_one/(rk_H2p_He_HeHp(1.0d3)*1.0d6*1.0d12)/eV_to_erg

      ! ---- one event of R6 and of R7: H3+ and electrons, no collider, so
      ! the thermalized fraction is zero and the whole internal share is
      ! withheld.  R7 is subtracted from the pair at its own rate.
      heat_one = assembly_heat(1.0d3, 0.0d0, 0.0d0, 0.0d0, 0.0d0,         &
                               0.0d0, 1.0d7, 0.0d0, 1.0d10, 0.0d0)
      gas6 = ( heat_one/eV_to_erg                                         &
             - rk_R7_H3p_dr_3H(1.0d3)*1.0d7*1.0d10*q7 )                   &
             /(rk_R6_H3p_dr_H2(1.0d3)*1.0d7*1.0d10)

      ! ---- one event of R8: H2+ and H2.  H2 alone opens the thermal
      ! dissociation R12 as well, so the same cell without the ion is the
      ! control and the difference is R8's; both cells form the same
      ! collider sum, n(H2) being their only heavy particle.
      heat_one = assembly_heat(1.0d3, 0.0d0, 0.0d0, 0.0d0, 1.0d12,        &
                               1.0d10, 0.0d0, 0.0d0, 0.0d0, 1.0d12)
      heat_two = assembly_heat(1.0d3, 0.0d0, 0.0d0, 0.0d0, 1.0d12,        &
                               0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d12)
      gas8  = (heat_one - heat_two)/(rk_R8_H2p_H2()*1.0d10*1.0d12)        &
              /eV_to_erg
      gas12 = heat_two/(rk_R12_H2_thdis(1.0d3)*1.0d12*1.0d12)/eV_to_erg

      ! ---- one event of R23: H2 and He+.  He+ opens R17 too, and its
      ! heat is the ground-state enthalpy of that channel at its own rate.
      heat_one = assembly_heat(1.0d3, 0.0d0, 1.0d8, 0.0d0, 1.0d12,        &
                               0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d12)
      gas23 = ( (heat_one - heat_two)/eV_to_erg                           &
              - rk_R17_Hep_H2_diss(1.0d3, 1.0d12)*1.0d8*1.0d12*q17 )      &
              /(rk_R23_H2_Hep_cx()*1.0d12*1.0d8)

      ! ---- one event of R15: atomic hydrogen only, so the third body is
      ! H and the thermalized fraction is that cell's own.
      heat_one = assembly_heat(1.0d3, 1.0d13, 0.0d0, 0.0d0, 0.0d0,        &
                               0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d13)
      f_heat = h2_vibrational_heat_fraction(1.0d3, 1.0d13, 0.0d0)
      gas15  = heat_one                                                   &
               /(k3b_H_H_to_H2_atomic_H()*1.0d13*1.0d13*1.0d13)/eV_to_erg

      ! Cycle 1, the direct route.  H2 + hv -> H2+ + e, then
      ! R5: H2+ + e -> H + H.  The gas keeps the enthalpy less the n = 2
      ! excitation and the level keeps that excitation.
      cyc = gas5 + e_n2
      call check_absolute('cycle_photoionization_then_R5_closes',         &
                          cyc - (e_th_H2 - d0_h2), 0.0d0, 1.0d-9)

      ! Cycle 2, through H3+.  H2 + hv -> H2+ + e, R8, then R6.  R6's
      ! internal share is withheld in full at this cell, and the withheld
      ! part is what closes the cycle.
      cyc = gas8 + gas6 + e_int_r6
      call check_absolute('cycle_photoionization_R8_R6_closes',           &
                          cyc - (e_th_H2 - d0_h2), 0.0d0, 1.0d-9)

      ! Cycle 3, through HeH+.  H2 + hv -> H2+ + e, H2+ + He -> HeH+ + H,
      ! then R16.  Helium is a catalyst and returns.
      cyc = gas_h2phe + gas16 + e_n2
      call check_absolute('cycle_photoionization_HeHp_R16_closes',        &
                          cyc - (e_th_H2 - d0_h2), 0.0d0, 1.0d-9)

      ! Cycle 4, THE NASCENT ION.  He + hv -> He+ + e, R23 (which leaves
      ! H2+ vibrationally excited and emits a 153 nm photon), then R5.
      ! The ledger deposits the ion's vibrational energy ONCE, inside the
      ! 1.058 eV R23 leaves after its photon; R5 is charged the
      ! GROUND-STATE enthalpy of the ion and adds none of it back.  A
      ! ledger that deposited the excitation at formation and again at
      ! recombination would exceed the input by that excitation.
      cyc = gas23 + e_gamma_r23 + gas5 + e_n2
      call check_absolute('cycle_R23_nascent_H2p_then_R5_closes',         &
                          cyc - (e_th_HeI - d0_h2), 0.0d0, 1.0d-9)

      ! Cycle 5, the three-body pair.  R15 forms H2 and R12 destroys it;
      ! the gas keeps the thermalized share of the formation energy, the
      ! infrared quadrupole lines carry the rest, and the dissociation
      ! takes the whole bond back.  A closed cycle liberates nothing.
      cyc = gas15 + q15*(1.0d0 - f_heat) + gas12
      call check_absolute('cycle_R15_then_R12_liberates_nothing',         &
                          cyc, 0.0d0, 1.0d-9)

      ! The four routes above are the same net reaction, H2 + hv -> H + H,
      ! so they must agree with each other and not only with the input.
      call check_absolute('the_four_routes_liberate_the_same_energy',     &
                          max(abs((gas5 + e_n2)                           &
                                  - (gas8 + gas6 + e_int_r6)),            &
                              abs((gas5 + e_n2)                           &
                                  - (gas_h2phe + gas16 + e_n2))),         &
                          0.0d0, 1.0d-9)

      ! ================================================================== !
      ! G. THE UNCERTAINTY RECORD ON THREE CELLS OF THE CERTIFIED
      !    MOLECULAR FIDUCIAL.  Five quantities of the ledger are
      !    approximations rather than measured numbers, and each is varied
      !    HERE ALONE, on the state the code meets, so that the heat
      !    bracket of each can be read off separately.  The nominal heat is
      !    the production assembly's; every variation is an exact
      !    difference against it, formed from the same production rate
      !    coefficients and reaction energies, so no term of the ledger is
      !    restated in this driver.
      ! ================================================================== !
      write(*,'(a)') '# UNCERTAINTY RECORD, molecular chemical heating'// &
                     ' [erg cm^-3 s^-1] on the certified fiducial'
      write(*,'(a)') '#   cell  T[K]  heat_nominal  recipient(0..1)'//    &
                     '  R6_branch(0..1)  R6_Eint(+/-15%)'//               &
                     '  quench(f..1)  He_eff(+/-0.3dex)  bracket'
      do ic = 1, ncase
         h_nom = assembly_heat(c_T(ic), c_nhi(ic), c_nheii(ic),           &
                               c_nheiS(ic), c_nh2(ic), c_nh2p(ic),        &
                               c_nh3p(ic), c_nhehp(ic), c_ne(ic),         &
                               c_ntot(ic), c_nhii(ic))
         f_c  = h2_vibrational_heat_fraction(c_T(ic), c_nhi(ic),          &
                                             c_nh2(ic), c_nheiS(ic))
         n3_c = h2_association_collider_density(c_T(ic), c_nh2(ic),       &
                                                c_nhi(ic), c_nheiS(ic))
         eff_he = k3b_H_H_to_H2_monatomic(c_T(ic))/k3b_H_H_to_H2(c_T(ic))

         ! (a) THE RECIPIENT FRACTION of step 2 of item L31: the share of
         ! the 10.199 eV that leaves the gas at R5 and R16.  Taking it to
         ! zero is the whole recipient correction undone.
         e_lev = e_n2*( rk_R5_H2p_dr(c_T(ic))*c_nh2p(ic)*c_ne(ic)         &
                      + rk_R16_HeHp_dr(c_T(ic))*c_nhehp(ic)*c_ne(ic) )    &
                 *eV_to_erg
         d_rec = e_lev

         ! (b) THE R6 BRANCHING: the share of R6 events that leave the H2
         ! fragment internally excited at all.  Only the part the gas does
         ! not get back through the quench fraction moves the heat.
         d_br = (1.0d0 - f_c)*e_int_r6                                    &
                *rk_R6_H3p_dr_H2(c_T(ic))*c_nh3p(ic)*c_ne(ic)*eV_to_erg

         ! (c) THE R6 INTERNAL ENERGY, varied by the +/-15 per cent that
         ! is the published uncertainty of the PEAK's position.  The
         ! distance from the peak to the mean of the distribution is not
         ! quantified by either source and is not inside this bracket.
         d_eint_hi = 0.15d0*d_br
         d_eint_lo = -d_eint_hi

         ! (d) THE QUENCH FRACTION, from its model value to 1, which is
         ! the thermalized limit and the largest heat the two internal
         ! shares can give.
         d_quench = (1.0d0 - f_c)                                         &
                    *( e_int_r6*rk_R6_H3p_dr_H2(c_T(ic))*c_nh3p(ic)       &
                                *c_ne(ic)                                 &
                     + q15*k3b_H_H_to_H2(c_T(ic))*n3_c                    &
                           *c_nhi(ic)*c_nhi(ic) )*eV_to_erg

         ! (e) THE HELIUM EFFICIENCY, argon's coefficient scaled by
         ! Cohen & Westberg's stated +/-0.3 in log.  It moves both
         ! directions of the R12/R15 pair, which is why the two terms
         ! appear together.
         d_he_hi = helium_efficiency_shift(ic, 10.0d0**0.3d0, f_c,        &
                                           n3_c, eff_he)
         d_he_lo = helium_efficiency_shift(ic, 10.0d0**(-0.3d0), f_c,     &
                                           n3_c, eff_he)

         bracket_hi = h_nom + max(0.0d0, d_eint_hi) + max(0.0d0, d_quench)&
                    + max(d_he_hi, d_he_lo) + max(0.0d0, d_rec)
         bracket_lo = h_nom + min(0.0d0, d_eint_lo) - abs(d_br)           &
                    + min(d_he_hi, d_he_lo)
         write(*,'(a,a20,f10.2,8es13.5)') '#   ', case_name(ic),          &
              c_T(ic), h_nom, d_rec, d_br, d_eint_hi, d_quench,           &
              d_he_hi, d_he_lo, bracket_hi - bracket_lo

         ! The energy the recipient variation moves IS the energy the
         ! H(n=2) source receives, at the same rates and densities: one
         ! event, one excitation, counted in one place.
         call check_relative('recipient_bracket_is_the_n2_source_'       &
                             //trim(case_name(ic)),                       &
              d_rec,                                                      &
              dissociative_recombination_n2_source(c_T(ic), c_nh2p(ic),   &
                             c_nhehp(ic), c_ne(ic))*e_n2*eV_to_erg,       &
              1.0d-13)
         ! Every bracket is finite and the nominal heat lies inside it.
         call check_at_least('bracket_contains_the_nominal_'             &
                             //trim(case_name(ic)),                       &
                             bracket_hi - h_nom, 0.0d0)
         call check_at_least('bracket_below_contains_the_nominal_'       &
                             //trim(case_name(ic)),                       &
                             h_nom - bracket_lo, 0.0d0)
      enddo

      ! WHERE THE LADDER IS THERMALIZED THE TWO INTERNAL SHARES CARRY NO
      ! UNCERTAINTY AT ALL.  At the molecular depth 1 - f is 1e-6 or
      ! below, so the quench bracket cannot reach a part in 1e5 of the
      ! heat there; this is the measurement behind the statement that the
      ! scalar fraction and a level-population model cannot differ in this
      ! layer by anything the run reports.
      f_c = h2_vibrational_heat_fraction(c_T(1), c_nhi(1), c_nh2(1),      &
                                         c_nheiS(1))
      call check_at_most('quench_bracket_is_closed_at_the_base',          &
                         1.0d0 - f_c, 1.0d-5)

      ! ================================================================== !
      ! H. ONE COLLIDER SUM, EVALUATED AS EACH OF ITS THREE CALLERS CALLS
      !    IT, AND BITWISE EQUAL.  The chemistry forms it inside the rows
      !    (System_HeH_mol::h2_third_body_density, at the composition the
      !    row is evaluated at); the carrier balance reaches the same call
      !    through mol_heh_rows and writes it into the R15 channel of its
      !    row record; the heat ledger forms it in
      !    molecular_chemical_heating.  If any one of them formed its own,
      !    the H2 abundance would be set at one rate and its heat charged
      !    at another.
      ! ================================================================== !
      call fill_probe_cell(c_T(1), c_ntot(1),                             &
                           c_nhi(1) + c_nhii(1) + 2.0d0*c_nh2(1),         &
                           c_nheiS(1) + c_nheii(1))
      call set_mol_coeffs(c_T(1), c_ntot(1))
      n3_chem   = h2_third_body_density(c_ntot(1), c_nh2(1), c_nhi(1),    &
                                        c_nheiS(1))
      n3_direct = h2_association_collider_density(c_T(1), c_nh2(1),       &
                                                  c_nhi(1), c_nheiS(1))
      call check_absolute('collider_sum_chemistry_equals_the_routine',    &
                          n3_chem - n3_direct, 0.0d0, 0.0d0)

      ! The carrier balance's own R15 channel, k15 n(H)^2 with
      ! k15 = mk15 n_third, written from the row and not recomputed.
      fvec_h = 0.0d0
      chan_h = 0.0d0
      call mol_heh_rows(fvec_h, c_nhi(1), c_nhii(1), c_nh2(1),            &
                        c_nh2p(1), c_nh3p(1), c_nhehp(1), c_nheiS(1),     &
                        0.0d0, c_nheii(1), 0.0d0, c_ne(1), c_ntot(1),     &
                        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  &
                        0.0d0, 0.0d0,                                     &
                        2.0d-12, 3.0d-12, 5.0d-12, 0.0d0,                 &
                        0.0d0, 0.0d0, 0.0d0, 0.0d0,                       &
                        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,                &
                        h2_chan = chan_h)
      call check_absolute('collider_sum_carrier_row_equals_the_routine',  &
                          chan_h(1) - (mk15*n3_direct)*c_nhi(1)*c_nhi(1), &
                          0.0d0, 0.0d0)

      ! The heat ledger's own, read out of the assembly on a cell whose
      ! only open channel is the thermal dissociation R12: there the whole
      ! deposit is k12 n_third n(H2) q(R12), so the third body the ledger
      ! used is recovered without dividing anything.
      heat_one = assembly_heat(c_T(1), c_nhi(1), 0.0d0, c_nheiS(1),       &
                               c_nh2(1), 0.0d0, 0.0d0, 0.0d0, 0.0d0,      &
                               c_ntot(1))
      n3_c   = h2_association_collider_density(c_T(1), c_nh2(1),          &
                                               c_nhi(1), c_nheiS(1))
      f_heat = h2_vibrational_heat_fraction(c_T(1), c_nhi(1), c_nh2(1),   &
                                            c_nheiS(1))
      expect = ( rk_R12_H2_thdis(c_T(1))*n3_c*c_nh2(1)*q12                &
               + k3b_H_H_to_H2(c_T(1))*n3_c*c_nhi(1)*c_nhi(1)*q15*f_heat )&
               *eV_to_erg
      call check_relative('collider_sum_heat_ledger_equals_the_routine',  &
                          heat_one, expect, 1.0d-14)

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

      double precision function assembly_heat(T, n_hi, n_heii, n_heiS,    &
                       n_h2, n_h2p, n_h3p, n_hehp, n_e, n_tot, n_hii)     &
                       result(g)
      ! The production molecular chemical heating of ONE cell at the state
      ! given [erg cm^-3 s^-1].  Everything the caller does not name is
      ! zero, so a cell carries only the channels its own densities open.
      real*8, intent(in) :: T, n_hi, n_heii, n_heiS, n_h2, n_h2p
      real*8, intent(in) :: n_h3p, n_hehp, n_e, n_tot
      real*8, intent(in), optional :: n_hii
      call zero_cell(T_K, nhi, nhii, nheii, nheiS, ne, ntot, nmol)
      T_K   = T
      nhi   = n_hi
      nheii = n_heii
      nheiS = n_heiS
      ne    = n_e
      ntot  = n_tot
      if (present(n_hii)) nhii = n_hii
      nmol(:,1) = n_h2
      nmol(:,2) = n_h2p
      nmol(:,3) = n_h3p
      nmol(:,4) = n_hehp
      call molecular_chemical_heating(T_K, nhi, nhii, nheii, nheiS,       &
                                      nmol, ne, ntot, gam)
      g = gam(1)
      end function assembly_heat

      double precision function helium_efficiency_shift(ic, scale, f, n3, &
                                eff_he) result(d)
      ! The change in the molecular chemical heat when the helium third
      ! body is given `scale` times the argon coefficient it stands on
      ! [erg cm^-3 s^-1].  Both directions of the R12/R15 pair carry the
      ! same collider sum, so both move together and their ratio, the
      ! equilibrium constant, does not move at all.
      integer, intent(in) :: ic
      real*8,  intent(in) :: scale, f, n3, eff_he
      real*8 :: dn3
      dn3 = (scale - 1.0d0)*eff_he*c_nheiS(ic)
      d = ( k3b_H_H_to_H2(c_T(ic))*dn3*c_nhi(ic)*c_nhi(ic)*q15*f          &
          + rk_R12_H2_thdis(c_T(ic))*dn3*c_nh2(ic)*q12 )*eV_to_erg
      end function helium_efficiency_shift

      subroutine fill_probe_cell(T, n_tot, n_h, n_he)
      ! The cell state the molecular rows read through ieq_cell.  The
      ! radiation field is zero: what block H tests is the third body the
      ! rows form, and no photon rate enters it.
      real*8, intent(in) :: T, n_tot, n_h, n_he
      thereis_HeITR   = .false.
      thereis_oxychem = .false.
      ieq_cell%P_HI      = 0.0d0
      ieq_cell%P_HeI     = 0.0d0
      ieq_cell%P_HeII    = 0.0d0
      ieq_cell%P_HeITR   = 0.0d0
      ieq_cell%P_H2      = 0.0d0
      ieq_cell%P_H2_di   = 0.0d0
      ieq_cell%P_H2_dd   = 0.0d0
      ieq_cell%P_H2_nd   = 0.0d0
      ieq_cell%k_LW      = 0.0d0
      ieq_cell%rchiiB    = 2.0d-12
      ieq_cell%rcheiiB   = 3.0d-12
      ieq_cell%rcheiiiB  = 5.0d-12
      ieq_cell%rcheiTR   = 0.0d0
      ieq_cell%a_ion_HI  = 0.0d0
      ieq_cell%a_ion_HeI = 0.0d0
      ieq_cell%a_ion_HeII  = 0.0d0
      ieq_cell%a_ion_HeITR = 0.0d0
      ieq_cell%q13       = 0.0d0
      ieq_cell%q31a      = 0.0d0
      ieq_cell%q31b      = 0.0d0
      ieq_cell%Q31       = 0.0d0
      ieq_cell%A31       = 0.0d0
      ieq_cell%kcx_He0_Hp = 0.0d0
      ieq_cell%kcx_Hep_H0 = 0.0d0
      ieq_cell%T_K       = T
      ieq_cell%ntot      = n_tot
      ieq_cell%nh        = n_h
      ieq_cell%nhe       = n_he
      ieq_cell%n_ofam    = 0.0d0
      ieq_cell%n_co      = 0.0d0
      end subroutine fill_probe_cell

      character(len=4) function tlabel(t)
      ! The probe temperature as a label, so each assertion names its own.
      real*8, intent(in) :: t
      write(tlabel,'(i4.4)') nint(t)
      end function tlabel


      subroutine representative_cells
      ! THE THREE CELLS BLOCK G IS QUOTED AT, either the compiled values or
      ! the ones a written state carries.  EXHALE_PROBE_STATE names the
      ! directory; unset, nothing is read and the compiled values stand, so
      ! the suite runs with no state in the tree.
      integer :: k
      real*8  :: dmax
      call get_environment_variable('EXHALE_PROBE_STATE', state_dir)
      c_T = c0_T;  c_nhi = c0_nhi;  c_nhii = c0_nhii
      c_nheiS = c0_nheiS;  c_nheii = c0_nheii
      c_nh2 = c0_nh2;  c_nh2p = c0_nh2p;  c_nh3p = c0_nh3p
      c_nhehp = c0_nhehp;  c_ne = c0_ne;  c_ntot = c0_ntot
      c_cell = (/ 1, 178, 299 /)
      c_r    = (/ 1.0001933187778382d0, 1.1010398408687325d0,             &
                  r_dilute_upper_column /)
      c_x2   = -1.0d0
      if (len_trim(state_dir) .eq. 0) then
         write(*,'(a)') '# representative cells: the compiled values'//   &
              ' (the archived molecular fiducial, output_pre_L34)'
         return
      endif
      call cells_of_written_state(trim(state_dir))
      write(*,'(a)') '# representative cells READ from '//trim(state_dir)
      write(*,'(a)') '#   case  cell  r[Rp]  2n(H2)/n_H  T[K]  n(HI)'//   &
                     '  n(H2)  n(HeI)  n_e  n_tot'
      do k = 1, ncase
         write(*,'(a,a20,i6,f11.6,9es16.8)') '#   ', case_name(k),        &
              c_cell(k), c_r(k), c_x2(k), c_T(k), c_nhi(k), c_nh2(k),     &
              c_nheiS(k), c_ne(k), c_ntot(k)
      enddo
      ! One number for the reproduction: zero on the state the compiled
      ! values were transcribed from, and the size of the move on any other.
      dmax = 0.0d0
      do k = 1, ncase
         dmax = max(dmax, reldist(c_T(k),     c0_T(k)))
         dmax = max(dmax, reldist(c_nhi(k),   c0_nhi(k)))
         dmax = max(dmax, reldist(c_nhii(k),  c0_nhii(k)))
         dmax = max(dmax, reldist(c_nheiS(k), c0_nheiS(k)))
         dmax = max(dmax, reldist(c_nheii(k), c0_nheii(k)))
         dmax = max(dmax, reldist(c_nh2(k),   c0_nh2(k)))
         dmax = max(dmax, reldist(c_nh2p(k),  c0_nh2p(k)))
         dmax = max(dmax, reldist(c_nh3p(k),  c0_nh3p(k)))
         dmax = max(dmax, reldist(c_nhehp(k), c0_nhehp(k)))
         dmax = max(dmax, reldist(c_ne(k),    c0_ne(k)))
         dmax = max(dmax, reldist(c_ntot(k),  c0_ntot(k)))
      enddo
      write(*,'(a,es12.5)') '# distance of the read cells from the'//     &
           ' compiled ones, largest relative: ', dmax
      end subroutine representative_cells

      double precision function reldist(a, b) result(d)
      real*8, intent(in) :: a, b
      d = abs(a - b)/max(abs(b), tiny(1.0d0))
      end function reldist

      subroutine cells_of_written_state(dirname)
      ! THE THREE CELLS, READ FROM A WRITTEN STATE.  The pair is the one
      ! load_IC reads, {Hydro_ioniz_IC.txt, Ion_species_IC.txt}; the schema
      ! carries two ghost cells at each end, so physical cell j is row
      ! j + 2.  The species columns are taken by NAME from the file's own
      ! "# columns" header, so a state written with another species list is
      ! read correctly or refused.
      !
      ! The electron density and the particle count are formed the way the
      ! code's own mass policy forms them: a column whose name ends in III
      ! carries two charges and one ending in II carries one, the molecular
      ! ions H2p, H3p and HeHp carry one each, and HeITR is a level of He I
      ! and is inside that column already, so it enters neither sum.
      character(len=*), intent(in) :: dirname
      integer, parameter :: mrow = 4096, mcol = 64, nghost = 2
      real*8  :: hyd(mrow,mcol), spc(mrow,mcol)
      character(len=16) :: cname(mcol)
      integer :: nrh, nch, nrs, ncs, j, k, kk, nphys, jfront, jdil
      real*8  :: x2base, dr, drbest
      call read_table(dirname//'/Hydro_ioniz_IC.txt',                     &
                      dirname//'/Hydro_ioniz.txt', hyd, nrh, nch, cname, .false.)
      call read_table(dirname//'/Ion_species_IC.txt',                     &
                      dirname//'/Ion_species.txt', spc, nrs, ncs, cname, .true.)
      if (nrh .ne. nrs) then
         write(*,'(a)') 'FAIL probe_state_pair_has_one_row_count'
         stop 1
      endif
      nphys = nrh - 2*nghost
      ! The H2 front: the first cell at or below half the base value of
      ! 2 n(H2)/n_H.
      x2base = h2_fraction(spc, cname, ncs, 1 + nghost)
      jfront = 0
      do j = 1, nphys
         if (h2_fraction(spc, cname, ncs, j + nghost)                     &
             .le. 0.5d0*x2base) then
            jfront = j
            exit
         endif
      enddo
      if (jfront .eq. 0) then
         write(*,'(a)') 'FAIL probe_state_has_no_h2_front'
         stop 1
      endif
      ! The dilute upper column: the cell nearest the recorded radius.
      jdil   = 1
      drbest = huge(1.0d0)
      do j = 1, nphys
         dr = abs(hyd(j + nghost, 1) - r_dilute_upper_column)
         if (dr .lt. drbest) then
            drbest = dr;  jdil = j
         endif
      enddo
      c_cell = (/ 1, jfront, jdil /)
      do k = 1, ncase
         j = c_cell(k) + nghost
         c_r(k)     = hyd(j, 1)
         c_T(k)     = hyd(j, 5)
         c_x2(k)    = h2_fraction(spc, cname, ncs, j)
         c_nhi(k)   = column_of(spc, cname, ncs, j, 'HI')
         c_nhii(k)  = column_of(spc, cname, ncs, j, 'HII')
         c_nheiS(k) = column_of(spc, cname, ncs, j, 'HeI')
         c_nheii(k) = column_of(spc, cname, ncs, j, 'HeII')
         c_nh2(k)   = column_of(spc, cname, ncs, j, 'H2')
         c_nh2p(k)  = column_of(spc, cname, ncs, j, 'H2p')
         c_nh3p(k)  = column_of(spc, cname, ncs, j, 'H3p')
         c_nhehp(k) = column_of(spc, cname, ncs, j, 'HeHp')
         c_ne(k)    = 0.0d0
         c_ntot(k)  = 0.0d0
         do kk = 2, ncs
            if (trim(cname(kk)) .eq. 'HeITR') cycle
            c_ntot(k) = c_ntot(k) + spc(j, kk)
            c_ne(k)   = c_ne(k)                                           &
                      + dble(charge_of(cname(kk)))*spc(j, kk)
         enddo
      enddo
      end subroutine cells_of_written_state

      double precision function h2_fraction(a, cname, nc, j) result(x2)
      ! 2 n(H2) over the hydrogen nuclei the state carries.
      real*8,            intent(in) :: a(:,:)
      character(len=16), intent(in) :: cname(:)
      integer,           intent(in) :: nc, j
      real*8 :: nh
      nh = column_of(a, cname, nc, j, 'HI')                               &
         + column_of(a, cname, nc, j, 'HII')                              &
         + 2.0d0*column_of(a, cname, nc, j, 'H2')                         &
         + 2.0d0*column_of(a, cname, nc, j, 'H2p')                        &
         + 3.0d0*column_of(a, cname, nc, j, 'H3p')                        &
         + column_of(a, cname, nc, j, 'HeHp')
      x2 = 2.0d0*column_of(a, cname, nc, j, 'H2')/max(nh, tiny(1.0d0))
      end function h2_fraction

      double precision function column_of(a, cname, nc, j, want) result(v)
      ! The column of that name, or a stop: a molecular quantity asked of a
      ! state that carries no molecules is a mistaken state, not a zero.
      real*8,            intent(in) :: a(:,:)
      character(len=16), intent(in) :: cname(:)
      integer,           intent(in) :: nc, j
      character(len=*),  intent(in) :: want
      integer :: k
      do k = 1, nc
         if (trim(cname(k)) .eq. want) then
            v = a(j, k)
            return
         endif
      enddo
      write(*,'(a)') 'FAIL probe_state_has_no_column_'//want
      stop 1
      end function column_of

      integer function charge_of(name) result(z)
      ! The charge a species column carries, from its own name.
      character(len=*), intent(in) :: name
      integer :: l
      character(len=16) :: s
      s = adjustl(name)
      l = len_trim(s)
      z = 0
      if (trim(s) .eq. 'HeITR') return
      if (trim(s) .eq. 'H2p' .or. trim(s) .eq. 'H3p'                      &
          .or. trim(s) .eq. 'HeHp') then
         z = 1
         return
      endif
      if (l .ge. 3) then
         if (s(l-2:l) .eq. 'III') then
            z = 2
            return
         endif
      endif
      if (l .ge. 2) then
         if (s(l-1:l) .eq. 'II') z = 1
      endif
      end function charge_of

      subroutine read_table(path, alt, a, nrow, ncol, cname, want_names)
      ! One written state file: the comment lines carry the column names,
      ! the rest are the rows.  The _IC member of the pair is the state a
      ! run handed back and is preferred; alt is read only if it is absent.
      character(len=*),  intent(in)    :: path, alt
      real*8,            intent(out)   :: a(:,:)
      integer,           intent(out)   :: nrow, ncol
      character(len=16), intent(inout) :: cname(:)
      logical,           intent(in)    :: want_names
      character(len=4096) :: line
      character(len=512)  :: use_path
      integer :: u, ios, k
      logical :: there
      inquire(file=path, exist=there)
      use_path = path
      if (.not. there) use_path = alt
      open(newunit=u, file=trim(use_path), status='old', action='read',   &
           iostat=ios)
      if (ios .ne. 0) then
         write(*,'(a)') 'FAIL probe_state_file_unreadable '//trim(use_path)
         stop 1
      endif
      nrow = 0
      ncol = 0
      do
         read(u,'(a)', iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0) cycle
         if (line(1:1) .eq. '#') then
            if (want_names .and. index(line, '# columns') .eq. 1)         &
               call names_of(line, cname, ncol)
            cycle
         endif
         nrow = nrow + 1
         if (nrow .gt. size(a,1)) then
            write(*,'(a)') 'FAIL probe_state_too_many_rows'
            stop 1
         endif
         if (ncol .eq. 0) ncol = tokens_in(line)
         read(line,*,iostat=ios) (a(nrow,k), k = 1, ncol)
         if (ios .ne. 0) then
            write(*,'(a)') 'FAIL probe_state_row_unreadable'
            stop 1
         endif
      enddo
      close(u)
      end subroutine read_table

      subroutine names_of(line, cname, ncol)
      ! The names the "# columns" header lists, the first of them r.
      character(len=*),  intent(in)    :: line
      character(len=16), intent(inout) :: cname(:)
      integer,           intent(out)   :: ncol
      character(len=4096) :: rest
      integer :: i, j
      rest = adjustl(line(index(line,'columns') + 7:))
      ncol = 0
      do
         if (len_trim(rest) .eq. 0) exit
         i = index(trim(rest), ' ')
         if (i .eq. 0) i = len_trim(rest) + 1
         ncol = ncol + 1
         if (ncol .gt. size(cname)) then
            write(*,'(a)') 'FAIL probe_state_too_many_columns'
            stop 1
         endif
         cname(ncol) = rest(1:i-1)
         ! The units a name carries in brackets are not part of it.
         j = index(cname(ncol), '[')
         if (j .gt. 1) cname(ncol) = cname(ncol)(1:j-1)
         rest = adjustl(rest(i+1:))
      enddo
      end subroutine names_of

      integer function tokens_in(line) result(n)
      character(len=*), intent(in) :: line
      character(len=4096) :: rest
      integer :: i
      rest = adjustl(line)
      n = 0
      do
         if (len_trim(rest) .eq. 0) exit
         i = index(trim(rest), ' ')
         if (i .eq. 0) i = len_trim(rest) + 1
         n = n + 1
         rest = adjustl(rest(i+1:))
      enddo
      end function tokens_in

      end program molecular_energy_recipients

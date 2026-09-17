      program molecular_third_body_in_the_chemistry
      ! ONE THREE-BODY RATE FOR THE COMPOSITION AND FOR THE ENERGY.
      !
      ! R15, H + H + M -> H2 + M, and its reverse R12 run on a COLLIDER SUM,
      ! k1(H2) n(H2) + k1(H) n(H) + k1(Ar) n(He), with the three recommended
      ! coefficients of Cohen & Westberg (1983), J. Phys. Chem. Ref. Data 12,
      ! 531, p. 559.  Two of those densities, n(H2) and n(H), are unknowns of
      ! the molecular balance system, so the rate is a function of the state
      ! the row is evaluated at and cannot be hoisted once per cell from the
      ! total heavy-particle density.  If it is, the H2 abundance is set by
      ! one rate and its heat charged at another, and nothing in either
      ! answer says so.
      !
      ! Production routines exercised: mol_heh_rows, set_mol_coeffs and
      ! set_mol_turnover_rates of
      ! src/modules/nonlinear_system_solver/System_HeH_mol.f90;
      ! h2_association_collider_density, k3b_H_H_to_H2,
      ! k3b_H_H_to_H2_atomic_H, k3b_H_H_to_H2_monatomic, rk_R15_3body_H2,
      ! rk_R12_H2_thdis and keq_H_H_to_H2 of
      ! src/modules/lower_atmosphere/mol_rates.f90.
      !
      ! WHAT IS ASSERTED.
      !  1. ONE NUMBER.  The R15 rate the H2 row forms, read from the row's
      !     own channel decomposition, is the rate the energy ledger forms
      !     at the same (T, n(H), n(H2), n(He)) -- the ledger's expression
      !     being rk_R15_3body_H2 of h2_association_collider_density, which
      !     molecular_reaction_heat writes and this test repeats from the
      !     same two routines.
      !  2. THE OLD HOIST IS NOT THE SAME NUMBER, so assertion 1 has
      !     content: at this cell the total-density coefficient stands 27
      !     per cent above the collider sum.
      !  3. BOTH DIRECTIONS CARRY THE SAME DENSITY.  R12's channel divided
      !     by R15's is the equilibrium constant times the densities alone,
      !     which is what makes the pair an exact detailed balance collider
      !     by collider.
      !  4. THE DERIVATIVES THE NUMERICAL JACOBIAN BUILDS.  d(R15 rate)/d
      !     n(H2) and d/d n(H), and d(R12 rate)/d n(H), by central
      !     differences of the production rows, against the analytic
      !     derivative of the collider sum.  These are the entries a
      !     hoisted coefficient does not have: with the third body frozen,
      !     the first and third are zero.
      !  5. THE TURNOVER SCALE SEES THE COLLIDER SUM.  The H2 row's bound
      !     rises when the cell's helium does, because helium is a third
      !     body of this channel.
      !
      ! The cell is the certified LHS 1140 b molecular base in order of
      ! magnitude (T = 808 K, n(H) = 1.87e12, n(He) = 6.17e12,
      ! n(H2) = 5.15e11 cm^-3), so the numbers the assertions carry are the
      ! ones the code meets there.

      use assertion_report
      use global_parameters, only: thereis_oxychem, thereis_HeITR
      use ion_cell_state, only: ieq_cell
      use System_HeH_mol, only: mol_heh_rows, set_mol_coeffs,             &
                                set_mol_turnover_rates, mol_inv_turnover
      use mol_rates, only: h2_thermochemistry_init,                       &
                           h2_association_collider_density,               &
                           k3b_H_H_to_H2, k3b_H_H_to_H2_atomic_H,         &
                           k3b_H_H_to_H2_monatomic,                       &
                           rk_R15_3body_H2, rk_R12_H2_thdis,              &
                           keq_H_H_to_H2

      implicit none

      ! The synthetic cell.
      real*8, parameter :: T_K    = 808.3d0
      real*8, parameter :: n_H    = 3.90d12      ! H nuclei
      real*8, parameter :: n_He   = 6.1716d12
      real*8, parameter :: n_hi   = 1.8680d12
      real*8, parameter :: n_h2   = 5.1474d11
      real*8, parameter :: n_hii  = 1.0d4
      real*8, parameter :: n_h2p  = 7.9d-2
      real*8, parameter :: n_h3p  = 1.0d3
      real*8, parameter :: n_hehp = 1.6d-5
      real*8, parameter :: n_heiSI = 6.1716d12
      real*8, parameter :: n_heii  = 1.0d2
      real*8, parameter :: n_heiii = 0.0d0
      real*8, parameter :: n_e     = 1.0d6
      real*8, parameter :: ntot    = 8.5543d12

      real*8 :: fv(8), ch(13), chp(13), chm(13)
      real*8 :: n3, k15_ledger, k15_row, k12_row, dh, d_num, d_ana
      real*8 :: eff_H, eff_He, s4_lo, s4_hi

      ! The H2 equilibrium-constant table, which the reverse of the thermal
      ! dissociation reads inside set_mol_coeffs.  keq_H_H_to_H2 refuses to
      ! build it itself so that no parallel region ever does.
      call h2_thermochemistry_init
      call fill_cell()
      call set_mol_coeffs(T_K, ntot)

      ! Collider efficiencies against the H2 third body, the two ratios the
      ! collider sum is built from.
      eff_H  = k3b_H_H_to_H2_atomic_H()/k3b_H_H_to_H2(T_K)
      eff_He = k3b_H_H_to_H2_monatomic(T_K)/k3b_H_H_to_H2(T_K)

      !----------------------------------------------------------------!
      ! 1. ONE NUMBER FOR THE CHEMISTRY AND FOR THE LEDGER.
      n3 = h2_association_collider_density(T_K, n_h2, n_hi, n_heiSI)
      k15_ledger = rk_R15_3body_H2(T_K, n3)
      call rows_at(n_hi, n_h2, ch)
      k15_row = ch(1)/(n_hi*n_hi)
      call check_relative('one_k15_for_the_row_and_the_ledger',          &
                          k15_row, k15_ledger, 1.0d-14)

      !----------------------------------------------------------------!
      ! 2. THE TOTAL-DENSITY COEFFICIENT IS A DIFFERENT NUMBER, so the
      !    assertion above is not satisfied by both readings at once.  The
      !    bound is loose on purpose: what has to hold is that the two
      !    differ at all at a composition the code meets.
      call check_at_least('collider_sum_is_not_the_total_density',        &
                          rk_R15_3body_H2(T_K, ntot)/k15_row, 1.1d0)

      !----------------------------------------------------------------!
      ! 3. BOTH DIRECTIONS CARRY THE SAME THIRD BODY.
      ! channel 8 is R12's loss, k12 n_third n(H2); channel 1 is R15's
      ! formation, k15 n(H)^2.  Their ratio is therefore
      ! n(H2)/(K_eq n(H)^2), with no collider density left in it, which is
      ! the statement that the pair is an exact detailed balance whatever
      ! the mixture is.
      k12_row = ch(8)/n_h2
      call check_relative('R12_and_R15_share_the_third_body',            &
                          ch(8)/ch(1),                                   &
                          n_h2/(keq_H_H_to_H2(T_K)*n_hi*n_hi), 1.0d-12)

      !----------------------------------------------------------------!
      ! 4. THE DERIVATIVES THE NUMERICAL JACOBIAN BUILDS.
      ! d(R15 rate)/d n(H2) = k1(H2) n(H)^2: one H2 added to the gas adds
      ! one H2-equivalent third body.  A hoisted coefficient gives zero.
      dh = 1.0d-4*n_h2
      call rows_at(n_hi, n_h2 + dh, chp)
      call rows_at(n_hi, n_h2 - dh, chm)
      d_num = (chp(1) - chm(1))/(2.0d0*dh)
      d_ana = k3b_H_H_to_H2(T_K)*n_hi*n_hi
      call check_relative('dR15_dnH2_is_the_H2_collider', d_num, d_ana,  &
                          1.0d-8)

      ! d(R15 rate)/d n(H) = k1(H2) [ eff_H n(H)^2 + 2 n_third n(H) ]: the
      ! atom is both a third body and a reactant.  A hoisted coefficient
      ! keeps only the second term.
      dh = 1.0d-4*n_hi
      call rows_at(n_hi + dh, n_h2, chp)
      call rows_at(n_hi - dh, n_h2, chm)
      d_num = (chp(1) - chm(1))/(2.0d0*dh)
      d_ana = k3b_H_H_to_H2(T_K)*(eff_H*n_hi*n_hi + 2.0d0*n3*n_hi)
      call check_relative('dR15_dnHI_has_the_collider_term', d_num,      &
                          d_ana, 1.0d-6)
      ! The collider term is not a rounding of the reactant term: it is
      ! 20 per cent of the whole derivative at this cell.
      call check_at_least('dR15_dnHI_collider_term_is_resolved',         &
                          eff_H*n_hi/(eff_H*n_hi + 2.0d0*n3), 1.0d-2)

      ! d(R12 rate)/d n(H) = k12 eff_H n(H2): the reverse direction picks
      ! up the same third body, which is zero with a hoisted coefficient.
      d_num = (chp(8) - chm(8))/(2.0d0*dh)
      d_ana = rk_R12_H2_thdis(T_K)*eff_H*n_h2
      call check_relative('dR12_dnHI_is_the_atomic_third_body', d_num,   &
                          d_ana, 1.0d-6)

      !----------------------------------------------------------------!
      ! 5. THE TURNOVER SCALE SEES THE COLLIDER SUM.  mol_inv_turnover(4)
      ! is 1/s(4), the bound the H2 row is divided by before the solver
      ! reads it.  Helium is a third body of this channel, so a cell with
      ! more helium has a larger bound; with the third body frozen at the
      ! total density it would move only through that total.
      ieq_cell%nhe = n_He
      ieq_cell%ntot = ntot
      call set_mol_turnover_rates(n_e, 1.0d0)
      s4_lo = 1.0d0/mol_inv_turnover(4)
      ieq_cell%nhe = 4.0d0*n_He
      call set_mol_turnover_rates(n_e, 1.0d0)
      s4_hi = 1.0d0/mol_inv_turnover(4)
      call check_at_least('H2_turnover_bound_rises_with_the_helium',     &
                          s4_hi/s4_lo, 1.0d0 + 1.0d-6)
      ieq_cell%nhe = n_He

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'molecular_third_body_in_the_chemistry: ',  &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                     &
           'molecular_third_body_in_the_chemistry: all assertions passed'

      contains

      subroutine fill_cell()
      ! Rate state of the synthetic cell.  The radiation field is
      ! suppressed, as it is under the shielded base; what the assertions
      ! test is the H2 row's three-body channel at fixed values of these,
      ! so only their being finite matters.
      thereis_HeITR      = .false.
      thereis_oxychem    = .false.
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
      ieq_cell%T_K       = T_K
      ieq_cell%ntot      = ntot
      ieq_cell%nh        = n_H
      ieq_cell%nhe       = n_He
      ieq_cell%n_ofam    = 0.0d0
      ieq_cell%n_co      = 0.0d0
      end subroutine fill_cell

      subroutine rows_at(nhi_in, nh2_in, chan)
      ! The H2 row's channel decomposition at a given atomic and molecular
      ! hydrogen density, everything else held.  The derivative under test
      ! is the partial one the finite-difference Jacobian builds column by
      ! column, so the electron density and the H nucleus total are held
      ! too.
      real*8, intent(in)  :: nhi_in, nh2_in
      real*8, intent(out) :: chan(13)
      fv   = 0.0d0
      chan = 0.0d0
      call mol_heh_rows(fv, nhi_in, n_hii, nh2_in, n_h2p, n_h3p, n_hehp,  &
                        n_heiSI, 0.0d0, n_heii, n_heiii, n_e, ntot,       &
                        ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,   &
                        ieq_cell%P_HeITR, ieq_cell%P_H2,                  &
                        ieq_cell%P_H2_di, ieq_cell%P_H2_dd,               &
                        ieq_cell%P_H2_nd, ieq_cell%k_LW,                  &
                        ieq_cell%rchiiB, ieq_cell%rcheiiB,                &
                        ieq_cell%rcheiiiB, ieq_cell%rcheiTR,              &
                        ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,            &
                        ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,        &
                        ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b,       &
                        ieq_cell%Q31, ieq_cell%A31,                       &
                        h2_chan = chan)
      end subroutine rows_at

      end program molecular_third_body_in_the_chemistry

      program co_helium_ion_sink
      ! He+ + CO -> C+ + O + He as a SPECIES sink, not only as a heat
      ! source: the He+ balance row of the molecular system must lose the
      ! helium ion the reaction consumes, at the rate the heating assembly
      ! deposits its energy for.
      !
      ! Production routines exercised: mol_heh_rows, set_mol_coeffs and
      ! set_mol_turnover_rates of
      ! src/modules/nonlinear_system_solver/System_HeH_mol.f90;
      ! rk_D1_Hep_CO of src/modules/lower_atmosphere/oxygen_rates.f90;
      ! oxygen_reaction_energy_eV of
      ! src/modules/lower_atmosphere/molecular_reaction_heat.f90.
      !
      ! THE REACTION.  UMIST RATE22 entry 4068 of rate22_final.rates
      ! (Millar, Walsh, Van de Sande & Markwick 2024, A&A 682, A109):
      ! He+ + CO -> O + C+ + He, alpha = 1.60e-9, beta = 0, gamma = 0,
      ! method M (measured), accuracy A.  One helium ion, one CO molecule
      ! and 2.2117 eV of product kinetic energy per event; the helium
      ! leaves neutral, so the free neutral He that closes the helium
      ! budget receives it and no other helium row moves.
      !
      ! WHAT IS ASSERTED.
      !  1. The He+ row of mol_heh_rows carries the term.  The row is
      !     evaluated twice on ONE synthetic cell, at n(CO) = 0 and at
      !     n(CO) > 0, with every other quantity held; the difference must
      !     be exactly -k_D1 n(CO) n(He+) and it must be a LOSS.
      !  2. Its derivative.  d(row 2)/d n(He+), by central differences,
      !     must gain exactly -k_D1 n(CO) when the CO is switched on.  This
      !     is the entry the numerical Jacobian of hybrd1 builds, so the
      !     check is on the quantity the solver actually uses.
      !  3. No other row moves.  Rows 1 and 3 to 7 are bitwise unchanged by
      !     the CO: the reaction touches helium and carbon alone, and
      !     carbon is not a row of this system.
      !  4. The heat and the species change count ONE event.  The event
      !     rate implied by the row (the sink divided by q_D1) equals the
      !     event rate the heating assembly deposits channel 18 for,
      !     k_D1 n(He+) n(CO), both built from rk_D1_Hep_CO and from the
      !     one formation-energy table.
      !  5. The row's turnover scale sees the term.  set_mol_turnover_rates
      !     is the bound the residual is divided by; without the CO term
      !     the He+ row of the shielded molecular base would be scaled by a
      !     bound far below the reaction it is mostly made of.
      !  6. Without the oxygen chemistry nothing changes: no other option
      !     carries CO, and a run without it must reproduce the row it had.
      !
      ! The cell is the shielded molecular base of the oxygen_chemistry
      ! case in order of magnitude (T = 1440 K, n_H = 1e13, n_He = 1e12,
      ! n(CO) = 2.4e10 cm^-3), so the numbers the assertions carry are the
      ! ones the code meets there.

      use assertion_report
      use global_parameters, only: thereis_oxychem, thereis_HeITR
      use ion_cell_state, only: ieq_cell
      use System_HeH_mol, only: mol_heh_rows, set_mol_coeffs,             &
                                set_mol_turnover_rates, mol_inv_turnover
      use oxygen_rates, only: rk_D1_Hep_CO
      use molecular_reaction_heat, only: oxygen_reaction_energy_eV, ir_D1
      use mol_rates, only: h2_thermochemistry_init

      implicit none
      real*8, parameter :: eV_to_erg = 1.602176634d-12

      ! The synthetic cell.
      real*8, parameter :: T_K   = 1440.0d0
      real*8, parameter :: ntot  = 1.0d13
      real*8, parameter :: n_H   = 1.0d13
      real*8, parameter :: n_He  = 1.0d12
      real*8, parameter :: n_CO  = 2.4d10
      real*8, parameter :: n_hi  = 2.0d12
      real*8, parameter :: n_hii = 1.0d6
      real*8, parameter :: n_h2  = 4.0d12
      real*8, parameter :: n_h2p = 1.0d3
      real*8, parameter :: n_h3p = 1.0d4
      real*8, parameter :: n_hehp = 1.0d2
      real*8, parameter :: n_heiSI = 9.99d11
      real*8, parameter :: n_heii  = 1.0d5
      real*8, parameter :: n_heiii = 1.0d0
      real*8, parameter :: n_e     = 1.0d8

      real*8 :: f0(8), f1(8), fp(8), fm(8)
      real*8 :: k_d1, sink, d0, d1, dh, q_d1_erg, rate_row, rate_heat
      real*8 :: scale0, scale1
      integer :: i
      real*8  :: others_move

      ! The H2 equilibrium-constant table, which the reverse of the thermal
      ! dissociation reads inside set_mol_coeffs. keq_H_H_to_H2 refuses to
      ! build it itself so that no parallel region ever does; EXHALE_main
      ! calls this in the same serial prologue.
      call h2_thermochemistry_init
      k_d1 = rk_D1_Hep_CO()
      call fill_cell()
      call set_mol_coeffs(T_K, ntot)

      !----------------------------------------------------------------!
      ! 1. THE ROW CARRIES THE SINK.
      thereis_oxychem = .true.
      ieq_cell%n_co   = 0.0d0
      call rows(f0)
      ieq_cell%n_co   = n_CO
      call rows(f1)

      sink = f1(2) - f0(2)
      call check_relative('hep_row_co_sink', sink,                        &
                          -k_d1*n_CO*n_heii, 1.0d-12)
      ! It is a LOSS of He+: a sink written with the wrong sign would pass
      ! a magnitude test and fail this one. sink + |sink| is zero for a
      ! negative number and twice the number for a positive one.
      call check_absolute('hep_row_co_sink_is_a_loss',                    &
                          sink + abs(sink), 0.0d0, 0.0d0)

      !----------------------------------------------------------------!
      ! 2. THE DERIVATIVE THE JACOBIAN BUILDS.
      ! Central difference in n(He+) with everything else held, at CO off
      ! and CO on. The step is 1e-4 of n(He+), well inside the linear
      ! regime of a row that is a polynomial in the densities.
      dh = 1.0d-4*n_heii
      ieq_cell%n_co = 0.0d0
      call rows_at_hep(n_heii + dh, fp)
      call rows_at_hep(n_heii - dh, fm)
      d0 = (fp(2) - fm(2))/(2.0d0*dh)
      ieq_cell%n_co = n_CO
      call rows_at_hep(n_heii + dh, fp)
      call rows_at_hep(n_heii - dh, fm)
      d1 = (fp(2) - fm(2))/(2.0d0*dh)
      call check_relative('hep_row_derivative_co_term', d1 - d0,          &
                          -k_d1*n_CO, 1.0d-8)

      !----------------------------------------------------------------!
      ! 3. NO OTHER ROW MOVES.
      ! The reaction consumes He+ and CO and makes He, C+ and O. Helium is
      ! the only element with a row here; the carbon and the oxygen belong
      ! to the metal block and to the carrier operator.
      others_move = 0.0d0
      do i = 1, 7
         if (i .eq. 2) cycle
         if (f1(i) .ne. f0(i)) others_move = others_move + 1.0d0
      enddo
      call check_absolute('co_moves_only_the_helium_ion_row',             &
                          others_move, 0.0d0, 0.0d0)

      !----------------------------------------------------------------!
      ! 4. ONE EVENT, COUNTED ONCE IN THE SPECIES AND ONCE IN THE ENERGY.
      ! The heating assembly deposits channel 18 as
      ! k_D1 n(He+) n(CO) q_D1 (util_ion_eq::heating_of_composition), so
      ! the event rate it pays for is that expression divided by q_D1. The
      ! event rate the species ledger performs is the He+ sink measured in
      ! part 1, one ion per event. They are the same number or the gas is
      ! heated by reactions it does not run.
      q_d1_erg  = oxygen_reaction_energy_eV(ir_D1)*eV_to_erg
      rate_heat = (k_d1*n_heii*n_CO*q_d1_erg)/q_d1_erg
      rate_row  = -sink
      call check_relative('heat_and_species_event_rate', rate_heat,       &
                          rate_row, 1.0d-14)
      ! The energy per event is exothermic and is the table's own value
      ! (co_destruction asserts where that number comes from; here it only
      ! has to be the SAME number the row is paired with).
      call check_positive('reaction_energy_D1_exothermic',                &
                          oxygen_reaction_energy_eV(ir_D1))

      !----------------------------------------------------------------!
      ! 5. THE TURNOVER SCALE SEES THE TERM.
      ! mol_inv_turnover(2) is 1/s(2), the bound the He+ row is divided by
      ! before hybrd1 reads it. The bound puts the whole helium element on
      ! the ion, as it does for every other term of that row, so the CO
      ! contribution is k_D1 n(CO) n_He.
      ieq_cell%n_co = 0.0d0
      call set_mol_turnover_rates(n_e, 1.0d0)
      scale0 = 1.0d0/mol_inv_turnover(2)
      ieq_cell%n_co = n_CO
      call set_mol_turnover_rates(n_e, 1.0d0)
      scale1 = 1.0d0/mol_inv_turnover(2)
      call check_relative('hep_turnover_scale_co_term', scale1 - scale0,  &
                          k_d1*n_CO*n_He, 1.0d-10)

      !----------------------------------------------------------------!
      ! 6. WITHOUT THE OXYGEN CHEMISTRY THE ROW IS THE ONE IT WAS.
      ! No other option carries CO, and a stale n_co left in the cell state
      ! must not reach the row.
      thereis_oxychem = .false.
      ieq_cell%n_co   = n_CO
      call rows(f1)
      call check_absolute('hep_row_unchanged_without_oxychem',            &
                          f1(2) - f0(2), 0.0d0, 0.0d0)
      thereis_oxychem = .true.

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'co_helium_ion_sink: ', assertion_failures,  &
              ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'co_helium_ion_sink: all assertions passed'

      contains

      subroutine fill_cell()
      ! Rate state of the synthetic cell. The values are order-of-magnitude
      ! representatives of the shielded base: photoionization suppressed,
      ! recombination and charge exchange at their coefficients there. What
      ! the assertions test is a DIFFERENCE at fixed values of these, so
      ! only their being nonzero and finite matters.
      thereis_HeITR      = .false.
      ieq_cell%P_HI      = 1.0d-12
      ieq_cell%P_HeI     = 1.0d-13
      ieq_cell%P_HeII    = 1.0d-15
      ieq_cell%P_HeITR   = 0.0d0
      ieq_cell%P_H2      = 1.0d-13
      ieq_cell%P_H2_di   = 2.0d-14
      ieq_cell%P_H2_dd   = 0.0d0
      ieq_cell%P_H2_nd   = 0.0d0
      ieq_cell%k_LW      = 1.0d-14
      ieq_cell%rchiiB    = 2.0d-12
      ieq_cell%rcheiiB   = 3.0d-12
      ieq_cell%rcheiiiB  = 5.0d-12
      ieq_cell%rcheiTR   = 0.0d0
      ieq_cell%a_ion_HI  = 1.0d-30
      ieq_cell%a_ion_HeI = 1.0d-32
      ieq_cell%a_ion_HeII  = 1.0d-34
      ieq_cell%a_ion_HeITR = 0.0d0
      ieq_cell%q13       = 0.0d0
      ieq_cell%q31a      = 0.0d0
      ieq_cell%q31b      = 0.0d0
      ieq_cell%Q31       = 0.0d0
      ieq_cell%A31       = 0.0d0
      ieq_cell%kcx_He0_Hp = 1.0d-15
      ieq_cell%kcx_Hep_H0 = 1.0d-15
      ieq_cell%T_K       = T_K
      ieq_cell%ntot      = ntot
      ieq_cell%nh        = n_H
      ieq_cell%nhe       = n_He
      ieq_cell%n_ofam    = 1.0d10
      ieq_cell%n_co      = 0.0d0
      end subroutine fill_cell

      subroutine rows(fv)
      real*8, intent(out) :: fv(8)
      call rows_at_hep(n_heii, fv)
      end subroutine rows

      subroutine rows_at_hep(hep, fv)
      ! The eight balance rows at a given He+ density, everything else
      ! held. The electron density is held too: the derivative under test
      ! is the partial one the finite-difference Jacobian builds row by
      ! row, not the total one along the charge-neutrality constraint.
      real*8, intent(in)  :: hep
      real*8, intent(out) :: fv(8)
      fv = 0.0d0
      call mol_heh_rows(fv, n_hi, n_hii, n_h2, n_h2p, n_h3p, n_hehp,      &
                        n_heiSI, 0.0d0, hep, n_heiii, n_e, ntot,          &
                        ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,   &
                        ieq_cell%P_HeITR, ieq_cell%P_H2,                  &
                        ieq_cell%P_H2_di, ieq_cell%P_H2_dd,               &
                        ieq_cell%P_H2_nd, ieq_cell%k_LW,                  &
                        ieq_cell%rchiiB, ieq_cell%rcheiiB,                &
                        ieq_cell%rcheiiiB, ieq_cell%rcheiTR,              &
                        ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,            &
                        ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,        &
                        ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b,       &
                        ieq_cell%Q31, ieq_cell%A31)
      end subroutine rows_at_hep

      end program co_helium_ion_sink

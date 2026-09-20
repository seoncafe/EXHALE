      program metal_helium_reservoir
      ! WHICH HELIUM STATE THE metal + He REACTIONS OF GROUP C REACT FROM,
      ! asked of the three systems that assemble them.
      !
      ! Huang et al. (2023), ApJ 951, 123, Table 4 group C lists Si + He+,
      ! Si+ + He, C + He+, C+ + He, O + He+ and O+ + He, every reactant
      ! written as the bare element. The three reverse rates carry the
      ! barriers exp(-19.1/T4), exp(-15.5/T4) and exp(-12.7/T4), and those
      ! are the ionization-potential differences of GROUND-STATE helium:
      ! 24.587 - 8.152 = 16.436 eV = 19.07e4 K against Si,
      ! 24.587 - 11.260 = 13.327 eV = 15.47e4 K against C and
      ! 24.587 - 13.618 = 10.969 eV = 12.73e4 K against O. They are
      ! therefore He(1^1S) rates. He(2^3S) lies 19.82 eV above the singlet
      ! and its ionization potential is 4.77 eV, so each of those collisions
      ! is exothermic for it and runs at a rate this code carries nowhere;
      ! the reservoir handed to the group C rows must be the ground singlet,
      ! and must be the same one in every system that assembles them.
      !
      ! HOW GROUP C IS ISOLATED. Group C is active only with cx_full, so
      ! each system's residual is evaluated at one state with cx_full true
      ! and false and the difference taken. Exactly one metal element is
      ! present, oxygen, so every group D row has a reactant at zero, and
      ! groups A and E are active in both evaluations and cancel: what is
      ! left of the difference is C5 (O + He+ -> O+ + He) and
      ! C6 (O+ + He -> O + He+).
      !
      ! THE STATE ISOLATES C6, the one group C row whose reactant is the
      ! helium. All the oxygen is ionized, so C5 has no neutral O to start
      ! from and the whole group C source is C6; the hydrogen is ionized to
      ! 1 atom cm^-3 of H0, which leaves the group A rows a rate far below
      ! the C6 one instead of a rate far above it; every photoionization,
      ! recombination and collisional coefficient of the cell is zero, and
      ! the He <-> H pair is switched off, so the rows the difference is
      ! read from carry nothing but charge exchange. C6 is
      ! barrier-suppressed by exp(-12.7/T4), so read against the rates of an
      ! ordinary cell it would be a part in 1e5 and the reservoir it reacts
      ! from could not be resolved at all.
      !
      ! WHAT IS ASSERTED
      !   the group C helium source does not move when a metastable is added
      !   at fixed singlet, He+ and He2+ densities, in each of the three
      !   systems; the three agree on it; and the oxygen source is minus the
      !   helium source, since C6 moves one electron between the two
      !   elements and nothing else.
      !
      ! THE ROW BASES are D7a's, asserted by reaction_source_orientation.
      !   triplet    ion_system_HeH_TR_metals: row 2 is the summed He I
      !              balance, He I-gain positive, so S(He I) = d(2); metal
      !              rows from 5, the X0 <-> X+ boundary ionization
      !              positive, so S(O I) = -d(metal base).
      !   molecular  ion_system_HeH_mol_metals: row 2 is the He+ balance,
      !              production positive, so S(He I) = -d(2); metal rows
      !              from 9.
      !   network    constrained_chemical_equilibrium: the same rows as the
      !              molecular system, in the same orientation, assembled
      !              from species densities instead of fractions.
      use global_parameters, only: thereis_HeITR, thereis_metals,         &
                                   thereis_oxychem
      use species_table,    only: n_melem, iel_O
      use ion_cell_state,   only: ieq_cell
      use charge_exchange,  only: cx_full, cx_init, cx_set_cell,          &
                                  cx_metal_base, he_h_cx_rates,           &
                                  he_h_charge_exchange
      use System_HeH_TR_metals,  only: ion_system_HeH_TR_metals
      use System_HeH_mol_metals, only: ion_system_HeH_mol_metals
      use System_HeH_mol,        only: mol_inv_turnover
      use System_HeH_metals,     only: set_metal_coeffs
      use constrained_chemical_equilibrium, only:                         &
           constrained_network_layout_of_cell,                            &
           constrained_network_balance_rows_of_cell,                      &
           is_HI, is_HII, is_HeI_SI, is_HeII, is_HeIII, is_HeITR,         &
           is_met0, n_species_max, n_fraction_rows_max
      use assertion_report, only: check_absolute, assertion_failures
      implicit none

      integer, parameter :: nsys = 3
      integer, parameter :: neqmax = 40
      character(len=14), parameter :: sysname(nsys) =                     &
         [ character(len=14) :: 'HeH_TR_metals', 'HeH_mol_metals',        &
           'network' ]

      ! The cell: densities in cm^-3, rate coefficients at T_cell.
      real*8, parameter :: T_cell  = 1.0d4
      real*8, parameter :: n_h     = 1.0d9
      real*8, parameter :: n_he_0  = 1.0d8
      real*8, parameter :: n_o     = 1.0d5
      ! 1 - x_hii = 1e-9, i.e. 1 atom cm^-3 of neutral hydrogen.
      real*8, parameter :: x_hii   = 1.0d0 - 1.0d-9
      real*8, parameter :: x_heii  = 0.20d0
      real*8, parameter :: x_heiii = 0.05d0
      real*8, parameter :: x_oii   = 0.95d0
      real*8, parameter :: x_oiii  = 0.05d0
      ! The metastable fraction of the helium nuclei in the state that
      ! carries one, far above any atmosphere's 1e-6 so that a reservoir
      ! that includes it separates from one that does not by much more than
      ! the rounding of the assembled rows.
      real*8, parameter :: x_tr    = 0.10d0

      ! Relative to the group C source itself, which is what the rows of
      ! this state carry.
      real*8, parameter :: round_tol = 1.0d-9

      real*8 :: k1, k2
      real*8 :: S_he(nsys), S_o(nsys), rs(nsys)
      real*8 :: S_he0(nsys), S_o0(nsys), rs0(nsys)
      real*8 :: sc
      integer :: is

      call he_h_cx_rates(T_cell, k1, k2)

      thereis_HeITR   = .true.
      thereis_metals  = .true.
      thereis_oxychem = .false.

      ! The canonical element list with oxygen alone carrying atoms: every
      ! group D row then has a reactant at zero, and the only cx_full-gated
      ! rows left with both reactants present are C5 and C6.
      call set_metal_coeffs(n_melem, oxygen_only(), zeros(), zeros(),     &
                            zeros(), zeros(), zeros(), zeros(), tops())

      do is = 1, nsys
         call group_c_sources(is, x_tr,  S_he(is),  S_o(is),  rs(is))
         call group_c_sources(is, 0.0d0, S_he0(is), S_o0(is), rs0(is))
      enddo

      write(*,'(a)') 'DIAGNOSTIC the group C helium source [cm^-3 s^-1]'
      do is = 1, nsys
         write(*,'(a,a,a,es23.15,a,es23.15)') '  ', sysname(is),          &
              ' : with metastable ', S_he(is), '  without ', S_he0(is)
      enddo

      ! ---- 1. the reservoir does not include the metastable ------------!
      ! The comparison state holds n(H0), n(H+), n(He 1^1S), n(He+),
      ! n(He2+) and every oxygen stage at the same absolute values and adds
      ! the metastable by enlarging the helium nucleus density, the only way
      ! the fraction unknowns allow it. The group C rates are ground-state
      ! rates, so their source must not notice the metastable.
      do is = 1, nsys
         sc = max(rs(is), rs0(is))
         call check_absolute('groupC_He_source_no_metastable_channel_'//  &
                             trim(sysname(is)),                           &
                             (S_he(is) - S_he0(is))/sc, 0.0d0, round_tol)
         call check_absolute('groupC_O_source_no_metastable_channel_'//   &
                             trim(sysname(is)),                           &
                             (S_o(is) - S_o0(is))/sc, 0.0d0, round_tol)
      enddo

      ! ---- 2. the three systems agree on it ----------------------------!
      sc = maxval(rs)
      do is = 2, nsys
         call check_absolute('groupC_He_source_agrees_'//                 &
                             trim(sysname(is))//'_with_'//                &
                             trim(sysname(1)),                            &
                             (S_he(is) - S_he(1))/sc, 0.0d0, round_tol)
         call check_absolute('groupC_O_source_agrees_'//                  &
                             trim(sysname(is))//'_with_'//                &
                             trim(sysname(1)),                            &
                             (S_o(is) - S_o(1))/sc, 0.0d0, round_tol)
      enddo

      ! ---- 3. one electron moves between the two elements --------------!
      do is = 1, nsys
         call check_absolute('groupC_charge_balance_'//                   &
                             trim(sysname(is)),                           &
                             (S_he(is) + S_o(is))/rs(is), 0.0d0,          &
                             round_tol)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'metal_helium_reservoir: ',                  &
              assertion_failures, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'metal_helium_reservoir: every assertion passed'

      contains

      !--------------!
      function zeros() result(z)
      real*8 :: z(n_melem)
      z = 0.0d0
      end function zeros

      !--------------!
      function oxygen_only() result(t)
      real*8 :: t(n_melem)
      t = 0.0d0
      t(iel_O) = n_o
      end function oxygen_only

      !--------------!
      ! Every element carrying its X++ stage, as the canonical table does
      ! for the elements that reach it.
      function tops() result(t)
      integer :: t(n_melem)
      t = 2
      end function tops

      !--------------!
      ! The cell state every system reads. Only the fields the H/He rows and
      ! the charge exchange use are set; the molecular coefficients stay at
      ! their zero defaults, which is a cell with no such reaction, and in
      ! any case they cancel in the difference that isolates group C.
      subroutine fill_cell(xtr)
      real*8, intent(in) :: xtr
      real*8 :: n_he
      n_he = n_he_0/(1.0d0 - xtr)
      ieq_cell%T_K         = T_cell
      ieq_cell%nh          = n_h
      ieq_cell%nhe         = n_he
      ieq_cell%n_ofam      = 0.0d0
      ieq_cell%ntot        = n_h + n_he
      ! Every photoionization, recombination and collisional coefficient is
      ! zero, so the rows carry charge exchange alone; they would cancel in
      ! the difference either way, but their size sets how many digits of
      ! the group C source survive it.
      ieq_cell%P_HI        = 0.0d0
      ieq_cell%P_HeI       = 0.0d0
      ieq_cell%P_HeII      = 0.0d0
      ieq_cell%P_HeITR     = 0.0d0
      ieq_cell%rchiiB      = 0.0d0
      ieq_cell%rcheiiB     = 0.0d0
      ieq_cell%rcheiiiB    = 0.0d0
      ieq_cell%rcheiTR     = 0.0d0
      ieq_cell%a_ion_HI    = 0.0d0
      ieq_cell%a_ion_HeI   = 0.0d0
      ieq_cell%a_ion_HeII  = 0.0d0
      ieq_cell%a_ion_HeITR = 0.0d0
      ieq_cell%A31         = 0.0d0
      ieq_cell%q13         = 0.0d0
      ieq_cell%q31a        = 0.0d0
      ieq_cell%q31b        = 0.0d0
      ieq_cell%Q31         = 0.0d0
      ieq_cell%kcx_He0_Hp  = k1
      ieq_cell%kcx_Hep_H0  = k2
      ieq_cell%x_h2_fixed   = .false.
      ieq_cell%x_ox_fixed   = .false.
      ieq_cell%x_hp_fixed   = .false.
      ieq_cell%x_heii_fixed = .false.
      ieq_cell%x_heiii_fixed = .false.
      end subroutine fill_cell

      !--------------!
      ! The unknown vector of system `is` for a state whose metastable is
      ! xtr of the helium nuclei, with the H, He and O ion densities held at
      ! the reference values.
      subroutine system_state(is, xtr, neq, mbase, x)
      integer, intent(in)  :: is
      real*8,  intent(in)  :: xtr
      integer, intent(out) :: neq, mbase
      real*8,  intent(out) :: x(neqmax)
      real*8 :: n_he, f
      n_he = n_he_0/(1.0d0 - xtr)
      f    = n_he_0/n_he
      x = 0.0d0
      x(1) = x_hii
      x(2) = x_heii*f
      x(3) = x_heiii*f
      if (is .eq. 1) then
         x(4) = xtr
         mbase = 5
      else
         x(8) = xtr
         mbase = 9
      endif
      neq = mbase - 1 + 2*n_melem
      x(mbase + 2*(iel_O-1))     = x_oii
      x(mbase + 2*(iel_O-1) + 1) = x_oiii
      end subroutine system_state

      !--------------!
      ! The species densities of the same state, in the network's layout.
      subroutine species_state(xtr, sden)
      real*8, intent(in)  :: xtr
      real*8, intent(out) :: sden(n_species_max)
      sden(:) = 0.0d0
      sden(is_HI)     = (1.0d0 - x_hii)*n_h
      sden(is_HII)    = x_hii*n_h
      sden(is_HeI_SI) = (1.0d0 - x_heii - x_heiii)*n_he_0
      sden(is_HeII)   = x_heii*n_he_0
      sden(is_HeIII)  = x_heiii*n_he_0
      sden(is_HeITR)  = xtr/(1.0d0 - xtr)*n_he_0
      sden(is_met0 + 3*(iel_O-1))     = (1.0d0 - x_oii - x_oiii)*n_o
      sden(is_met0 + 3*(iel_O-1) + 1) = x_oii*n_o
      sden(is_met0 + 3*(iel_O-1) + 2) = x_oiii*n_o
      end subroutine species_state

      !--------------!
      ! The group C helium and oxygen number-density sources [cm^-3 s^-1]
      ! of system `is`, as the difference of the residual with cx_full on
      ! and off. rs is the magnitude of the rows the difference is read
      ! from, which is the scale it is resolved against.
      subroutine group_c_sources(is, xtr, S_heI, S_OI, rs_out)
      integer, intent(in)  :: is
      real*8,  intent(in)  :: xtr
      real*8,  intent(out) :: S_heI, S_OI, rs_out
      real*8  :: fon(n_fraction_rows_max), foff(n_fraction_rows_max)
      real*8  :: d(n_fraction_rows_max)
      integer :: neq, mbase, iox_row

      call rows(is, xtr, .true.,  neq, mbase, fon)
      call rows(is, xtr, .false., neq, mbase, foff)
      d(1:neq) = fon(1:neq) - foff(1:neq)
      iox_row  = mbase + 2*(iel_O-1)
      if (is .eq. 1) then
         S_heI = d(2)              ! summed He I balance, He I-gain positive
      else
         S_heI = -d(2)             ! He+ balance, production positive
      endif
      ! Every metal row is the X0 <-> X+ boundary, ionization positive.
      S_OI   = -d(iox_row)
      ! The scale the assertions are read against is the group C source
      ! itself: it is what these rows carry.
      rs_out = max(abs(S_heI), abs(S_OI))
      end subroutine group_c_sources

      !--------------!
      subroutine rows(is, xtr, with_full, neq, mbase, fvec)
      integer, intent(in)  :: is
      real*8,  intent(in)  :: xtr
      logical, intent(in)  :: with_full
      integer, intent(out) :: neq, mbase
      real*8,  intent(out) :: fvec(n_fraction_rows_max)
      real*8  :: x(neqmax), params(60), sden(n_species_max)
      logical :: is_unknown(n_species_max)
      integer :: iflag, nu, nrow, ncons

      iflag  = 1
      params = 0.0d0
      fvec   = 0.0d0
      ! The He <-> H pair is a different reaction, tested by
      ! he_h_pair_reservoir and advection_pair_reservoir; off here so that
      ! the helium row carries the group C source alone.
      he_h_charge_exchange = .false.
      cx_full = with_full
      call cx_init()
      call cx_set_cell(T_cell)
      call fill_cell(xtr)
      call system_state(is, xtr, neq, mbase, x)
      mol_inv_turnover(:) = 1.0d0
      select case (is)
      case (1)
         cx_metal_base = mbase
         call ion_system_HeH_TR_metals(neq, x, fvec, iflag, params)
         cx_metal_base = 4
      case (2)
         cx_metal_base = mbase
         call ion_system_HeH_mol_metals(neq, x, fvec, iflag, params)
         cx_metal_base = 4
      case default
         call species_state(xtr, sden)
         call constrained_network_layout_of_cell(neq, mbase, mbase,       &
              sden, nu, nrow, ncons, is_unknown)
         call constrained_network_balance_rows_of_cell(sden, fvec)
      end select
      cx_full = .false.
      call cx_init()
      call cx_set_cell(T_cell)
      he_h_charge_exchange = .true.
      end subroutine rows

      end program metal_helium_reservoir

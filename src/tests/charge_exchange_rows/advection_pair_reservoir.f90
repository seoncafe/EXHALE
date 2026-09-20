      program advection_pair_reservoir
      ! WHICH HELIUM STATE THE He <-> H CHARGE-EXCHANGE PAIR REACTS FROM IN
      ! THE ADVECTION-CORRECTED SYSTEM, and whether it is the same state the
      ! equilibrium system reacts from.
      !
      ! Huang et al. (2023), ApJ 951, 123, Table 4 lists He + H+ -> He+ + H
      ! from Glover & Jappsen (2007) with the barrier exp(-12.75/T4):
      ! 12.75e4 K is 10.99 eV, the ionization-potential difference
      ! 24.587 - 13.598 eV of GROUND-STATE helium against hydrogen, so the
      ! reactant is He(1^1S). He(2^3S) lies 19.82 eV above the singlet, its
      ! ionization potential is 4.77 eV, and the same collision is
      ! exothermic for it: a different reaction, at a different rate,
      ! carried nowhere in this code. System_implicit_adv_HeH_TR solves the
      ! singlet and the metastable as separate unknowns and its row 2 is the
      ! ground-singlet balance, so the reservoir it hands the pair must be
      ! the singlet x(2), and must be the one the equilibrium system uses.
      !
      ! HOW THE PAIR IS ISOLATED. The residual is evaluated twice at one
      ! state, with he_h_charge_exchange true and false, and the difference
      ! is the pair's contribution to each row and nothing else: every other
      ! term, and every coefficient the system reads, is identical between
      ! the two evaluations.
      !
      ! THE ROW BASES.
      !   advection  System_implicit_adv_HeH_TR: fvec(1) tracks n_HI per
      !              n_H and fvec(2) the ground singlet per n_He, both
      !              neutral-gain positive and both carrying the factor
      !              c1 = dr/v. The number-density sources are therefore
      !              S(H I) = d(1) n_H/c1 and S(He I) = d(2) n_He/c1.
      !   triplet    System_HeH_TR: row 1 is H I -> H II positive and row 2
      !              the summed He I balance, He I-gain positive, both in
      !              cm^-3 s^-1, so S(H I) = -d(1) and S(He I) = d(2).
      use global_parameters, only: thereis_HeITR, thereis_metals,         &
                                   thereis_oxychem
      use ion_cell_state,   only: ieq_cell, adv_cell
      use charge_exchange,  only: he_h_charge_exchange, he_h_cx_rates
      use System_implicit_adv_HeH_TR, only: adv_implicit_HeH_TR
      use System_HeH_TR,              only: ion_system_HeH_TR
      use assertion_report, only: check_absolute, assertion_failures
      implicit none

      ! The cell: densities in cm^-3, rate coefficients at T_cell.
      real*8, parameter :: T_cell  = 1.0d4
      real*8, parameter :: n_h     = 1.0d9
      real*8, parameter :: n_he_0  = 1.0d8
      real*8, parameter :: x_hii   = 0.30d0
      real*8, parameter :: x_heii  = 0.20d0
      real*8, parameter :: x_heiii = 0.05d0
      ! The metastable fraction of the helium nuclei in the state that
      ! carries one. Far above any atmosphere's 1e-6, so that a pair
      ! charged to the summed He I separates from one charged to the
      ! singlet by much more than the rounding of the assembled rows.
      real*8, parameter :: x_tr    = 0.10d0
      ! dr/v of the upwind integration, s. Chosen so that c1 times the
      ! cell's rates is of order one, which is where the advection system
      ! actually runs.
      real*8, parameter :: c1      = 1.0d4

      ! The rounding floor of a source read out of an assembled row.
      real*8, parameter :: round_tol = 1.0d-11

      real*8 :: k1, k2, R1, R2
      real*8 :: S_he_adv, S_h_adv, sc_adv
      real*8 :: S_he_adv0, S_h_adv0, sc_adv0
      real*8 :: S_he_eq, S_h_eq, sc_eq
      real*8 :: sc

      call he_h_cx_rates(T_cell, k1, k2)

      ! The closed form the pair must produce at the reference state. The
      ! ground singlet is n_he_0 (1 - x_HeII - x_HeIII) in every state
      ! below: the metastable is added on top of it.
      R1 = k1*n_he_0*(1.0d0 - x_heii - x_heiii)*x_hii*n_h   ! He(1^1S) + H+
      R2 = k2*x_heii*n_he_0*(1.0d0 - x_hii)*n_h             ! He+ + H0

      write(*,'(a)') 'DIAGNOSTIC the two Table 4 group B rate coefficients'
      write(*,'(a,es23.15)') '  k(He + H+)  [cm^3 s^-1]  = ', k1
      write(*,'(a,es23.15)') '  k(He+ + H)  [cm^3 s^-1]  = ', k2
      write(*,'(a,es23.15)') '  R1 = k1 n(He 1^1S) n(H+) = ', R1
      write(*,'(a,es23.15)') '  R2 = k2 n(He+) n(H0)     = ', R2
      write(*,'(a)')         '  k1 n(He 2^3S) n(H+), the term a summed'
      write(*,'(a,es23.15)') '    He I reservoir would add to R1  = ',    &
           k1*x_tr/(1.0d0 - x_tr)*n_he_0*x_hii*n_h

      thereis_HeITR   = .true.
      thereis_metals  = .false.
      thereis_oxychem = .false.

      ! ---- 1. the reservoir is the ground singlet ----------------------!
      call advection_pair_sources(x_tr, S_he_adv, S_h_adv, sc_adv)
      call check_absolute('adv_pair_He_source_singlet',                   &
                          (S_he_adv - (R2 - R1))/sc_adv, 0.0d0, round_tol)
      call check_absolute('adv_pair_H_source_singlet',                    &
                          (S_h_adv - (R1 - R2))/sc_adv, 0.0d0, round_tol)

      ! ---- 2. adding a metastable at fixed singlet moves nothing -------!
      ! The comparison state holds n(H0), n(H+), n(He 1^1S), n(He+) and
      ! n(He2+) at the same absolute values and adds the metastable by
      ! enlarging the helium nucleus density, the only way the fraction
      ! unknowns allow it. The He(2^3S) + H+ collision has its own rate and
      ! is carried nowhere, so the pair must not notice the metastable.
      call advection_pair_sources(0.0d0, S_he_adv0, S_h_adv0, sc_adv0)
      sc = max(sc_adv, sc_adv0)
      call check_absolute('adv_pair_He_source_no_metastable_channel',     &
                          (S_he_adv - S_he_adv0)/sc, 0.0d0, round_tol)
      call check_absolute('adv_pair_H_source_no_metastable_channel',      &
                          (S_h_adv - S_h_adv0)/sc, 0.0d0, round_tol)

      ! ---- 3. the same state as the equilibrium system reacts from -----!
      call equilibrium_pair_sources(x_tr, S_he_eq, S_h_eq, sc_eq)
      sc = max(sc_adv, sc_eq)
      call check_absolute('adv_pair_He_source_agrees_with_HeH_TR',        &
                          (S_he_adv - S_he_eq)/sc, 0.0d0, round_tol)
      call check_absolute('adv_pair_H_source_agrees_with_HeH_TR',         &
                          (S_h_adv - S_h_eq)/sc, 0.0d0, round_tol)

      write(*,'(a)') 'DIAGNOSTIC the pair sources [cm^-3 s^-1]'
      write(*,'(a,es23.15)') '  S(He I) advection-corrected = ', S_he_adv
      write(*,'(a,es23.15)') '  S(He I) equilibrium         = ', S_he_eq
      write(*,'(a,es23.15)') '  S(H I)  advection-corrected = ', S_h_adv
      write(*,'(a,es23.15)') '  S(H I)  equilibrium         = ', S_h_eq

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'advection_pair_reservoir: ',                &
              assertion_failures, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'advection_pair_reservoir: every assertion passed'

      contains

      !--------------!
      ! The advection cell, for a state whose metastable is xtr of the
      ! helium nuclei. The helium nucleus density grows by exactly the
      ! metastable density, so n(He 1^1S), n(He+) and n(He2+) are the
      ! reference ones at every xtr. Only the fields the H/He rows and the
      ! pair read are set; the rest stay at zero, which is a cell with no
      ! such reaction and in any case cancels in the difference.
      subroutine fill_adv_cell(xtr)
      real*8, intent(in) :: xtr
      real*8 :: n_he
      n_he = n_he_0/(1.0d0 - xtr)
      adv_cell%c1          = c1
      adv_cell%nh          = n_h
      adv_cell%heh_loc     = n_he/n_h
      adv_cell%xhi_old     = 1.0d0 - x_hii
      adv_cell%xheiS_old   = (1.0d0 - x_heii - x_heiii)*n_he_0/n_he
      adv_cell%xheiii_old  = x_heiii*n_he_0/n_he
      adv_cell%xheiTR_old  = xtr
      adv_cell%P_HI        = 1.0d-4
      adv_cell%P_HeI       = 3.0d-5
      adv_cell%P_HeII      = 1.0d-6
      adv_cell%P_HeITR     = 2.0d-3
      adv_cell%rchiiB      = 2.6d-13
      adv_cell%rcheiiB     = 4.3d-13
      adv_cell%rcheiiiB    = 1.5d-12
      adv_cell%rcheiTR     = 1.5d-13
      adv_cell%a_ion_HI    = 1.0d-11
      adv_cell%a_ion_HeI   = 1.0d-13
      adv_cell%a_ion_HeII  = 1.0d-15
      adv_cell%a_ion_HeITR = 1.0d-8
      adv_cell%A31         = 1.272d-4
      adv_cell%q13         = 2.1d-8
      adv_cell%q31a        = 5.0d-10
      adv_cell%q31b        = 2.0d-9
      adv_cell%Q31         = 5.0d-10
      adv_cell%xe_metal    = 0.0d0
      adv_cell%kcx_He0_Hp  = k1
      adv_cell%kcx_Hep_H0  = k2
      end subroutine fill_adv_cell

      !--------------!
      ! The equilibrium cell at the same physical state.
      subroutine fill_ieq_cell(xtr)
      real*8, intent(in) :: xtr
      real*8 :: n_he
      n_he = n_he_0/(1.0d0 - xtr)
      ieq_cell%T_K         = T_cell
      ieq_cell%nh          = n_h
      ieq_cell%nhe         = n_he
      ieq_cell%ntot        = n_h + n_he
      ieq_cell%P_HI        = 1.0d-4
      ieq_cell%P_HeI       = 3.0d-5
      ieq_cell%P_HeII      = 1.0d-6
      ieq_cell%P_HeITR     = 2.0d-3
      ieq_cell%rchiiB      = 2.6d-13
      ieq_cell%rcheiiB     = 4.3d-13
      ieq_cell%rcheiiiB    = 1.5d-12
      ieq_cell%rcheiTR     = 1.5d-13
      ieq_cell%a_ion_HI    = 1.0d-11
      ieq_cell%a_ion_HeI   = 1.0d-13
      ieq_cell%a_ion_HeII  = 1.0d-15
      ieq_cell%a_ion_HeITR = 1.0d-8
      ieq_cell%A31         = 1.272d-4
      ieq_cell%q13         = 2.1d-8
      ieq_cell%q31a        = 5.0d-10
      ieq_cell%q31b        = 2.0d-9
      ieq_cell%Q31         = 5.0d-10
      ieq_cell%kcx_He0_Hp  = k1
      ieq_cell%kcx_Hep_H0  = k2
      end subroutine fill_ieq_cell

      !--------------!
      ! The unknown vector of one system at the reference physical state.
      ! The two layouts share x(3) = n_HeIII/n_He and x(4) = n(He 2^3S)/n_He
      ! and differ in the first two entries: the advection system solves for
      ! the neutral fractions n_HI/n_H and n(He 1^1S)/n_He, the equilibrium
      ! system for the ionized ones n_HII/n_H and n_HeII/n_He.
      !
      ! The helium nucleus density grows by exactly the metastable density,
      ! so n(He 1^1S) = (1 - x_HeII - x_HeIII) n_he_0, n(He+) = x_HeII
      ! n_he_0 and n(He2+) = x_HeIII n_he_0 at every xtr.
      subroutine system_state(advection, xtr, x)
      logical, intent(in)  :: advection
      real*8,  intent(in)  :: xtr
      real*8,  intent(out) :: x(4)
      real*8 :: f
      f = 1.0d0 - xtr                     ! = n_he_0/n_he
      if (advection) then
         x(1) = 1.0d0 - x_hii
         x(2) = (1.0d0 - x_heii - x_heiii)*f
      else
         x(1) = x_hii
         x(2) = x_heii*f
      endif
      x(3) = x_heiii*f
      x(4) = xtr
      end subroutine system_state

      !--------------!
      ! The pair's physical He I and H I number-density sources [cm^-3
      ! s^-1] in the advection-corrected system, as the difference of two
      ! residual evaluations. rs is the magnitude of the two rows it writes
      ! into, converted the same way, which is the scale the difference is
      ! resolved against.
      subroutine advection_pair_sources(xtr, S_heI, S_hI, rs)
      real*8, intent(in)  :: xtr
      real*8, intent(out) :: S_heI, S_hI, rs
      real*8  :: x(4), fon(4), foff(4), params(25)
      real*8  :: n_he
      integer :: iflag
      n_he   = n_he_0/(1.0d0 - xtr)
      iflag  = 1
      params = 0.0d0
      call fill_adv_cell(xtr)
      call system_state(.true., xtr, x)
      he_h_charge_exchange = .true.
      call adv_implicit_HeH_TR(4, x, fon, iflag, params)
      he_h_charge_exchange = .false.
      call adv_implicit_HeH_TR(4, x, foff, iflag, params)
      he_h_charge_exchange = .true.
      S_hI  = (fon(1) - foff(1))*n_h/c1
      S_heI = (fon(2) - foff(2))*n_he/c1
      rs    = max(abs(fon(1)), abs(foff(1)))*n_h/c1
      rs    = max(rs, max(abs(fon(2)), abs(foff(2)))*n_he/c1)
      end subroutine advection_pair_sources

      !--------------!
      ! The same two sources in System_HeH_TR, whose rows are already in
      ! cm^-3 s^-1.
      subroutine equilibrium_pair_sources(xtr, S_heI, S_hI, rs)
      real*8, intent(in)  :: xtr
      real*8, intent(out) :: S_heI, S_hI, rs
      real*8  :: x(4), fon(4), foff(4), params(40)
      integer :: iflag
      iflag  = 1
      params = 0.0d0
      call fill_ieq_cell(xtr)
      call system_state(.false., xtr, x)
      he_h_charge_exchange = .true.
      call ion_system_HeH_TR(4, x, fon, iflag, params)
      he_h_charge_exchange = .false.
      call ion_system_HeH_TR(4, x, foff, iflag, params)
      he_h_charge_exchange = .true.
      S_hI  = -(fon(1) - foff(1))
      S_heI =  (fon(2) - foff(2))
      rs    = max(abs(fon(1)), abs(foff(1)), abs(fon(2)), abs(foff(2)))
      end subroutine equilibrium_pair_sources

      end program advection_pair_reservoir

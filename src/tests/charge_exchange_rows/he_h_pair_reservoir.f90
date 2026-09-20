      program he_h_pair_reservoir
      ! WHICH HELIUM STATE THE He <-> H CHARGE-EXCHANGE PAIR REACTS FROM,
      ! asked of the four ionization systems that carry helium and a
      ! metastable, by running their residuals and isolating the pair.
      !
      ! Huang et al. (2023), ApJ 951, 123, Table 4 lists He + H+ -> He+ + H
      ! and He+ + H -> He + H+ from Glover & Jappsen (2007). The forward
      ! rate carries the barrier exp(-12.75/T4): 12.75e4 K is 10.99 eV, the
      ! ionization-potential difference 24.587 - 13.598 eV of GROUND-STATE
      ! helium against hydrogen, so the reactant is He(1^1S) and not the sum
      ! over He I. He(2^3S) lies 19.82 eV above the singlet, its ionization
      ! potential is 4.77 eV, and the same collision is exothermic for it:
      ! a different reaction with a different rate, carried nowhere in this
      ! code. The reservoir the systems hand the pair must therefore be the
      ! ground singlet, and must be the same one in every system.
      !
      ! HOW THE PAIR IS ISOLATED. Each system's residual is evaluated twice
      ! at one state, with he_h_charge_exchange true and false, and the
      ! difference is the pair's contribution to each row and nothing else:
      ! every other term of the system, and every module coefficient it
      ! reads, is identical between the two evaluations.
      !
      ! THE ROW BASES are the ones D7a established and the suite's
      ! reaction_source_orientation test asserts. The pair touches only the
      ! H row and the He I <-> He II boundary row, so the physical sources
      ! it makes are, with u the net upward flow across a boundary,
      !   triplet   (System_HeH_TR, System_HeH_TR_metals): row 2 is the
      !             summed He I balance, He I-gain positive, u(He,0) =
      !             -fvec(2); row 1 is H I -> H II positive, u(H,0) = fvec(1).
      !   molecular (System_HeH_mol, System_HeH_mol_metals): row 2 is the
      !             He II balance, production positive, u(He,0) = fvec(2)
      !             for a reaction that does not touch He2+; row 1 is the
      !             H II balance, production positive, u(H,0) = fvec(1).
      ! S(He I) = -u(He,0) and S(H I) = -u(H,0) in both.
      use global_parameters, only: thereis_HeITR, thereis_metals,         &
                                   thereis_oxychem
      use species_table,    only: n_melem
      use ion_cell_state,   only: ieq_cell
      use charge_exchange,  only: he_h_charge_exchange, he_h_cx_rates,    &
                                  cx_init, cx_set_cell, cx_metal_base
      use System_HeH_TR,        only: ion_system_HeH_TR
      use System_HeH_TR_metals, only: ion_system_HeH_TR_metals
      use System_HeH_mol,       only: ion_system_HeH_mol, mol_inv_turnover
      use System_HeH_mol_metals, only: ion_system_HeH_mol_metals
      use System_HeH_metals,    only: set_metal_coeffs
      use assertion_report, only: check_absolute, check_relative,         &
                                  assertion_failures
      implicit none

      integer, parameter :: nsys = 4
      integer, parameter :: neqmax = 40
      character(len=14), parameter :: sysname(nsys) =                     &
         [ character(len=14) :: 'HeH_TR', 'HeH_TR_metals', 'HeH_mol',     &
           'HeH_mol_metals' ]

      ! The cell: densities in cm^-3, rate coefficients at T_cell.
      real*8, parameter :: T_cell  = 1.0d4
      real*8, parameter :: n_h     = 1.0d9
      real*8, parameter :: x_hii   = 0.30d0
      real*8, parameter :: x_heii  = 0.20d0
      real*8, parameter :: x_heiii = 0.05d0
      ! The metastable fraction of the helium nuclei in the state that
      ! carries one. Far above any atmosphere's 1e-6, so that a pair
      ! charged to the summed He I separates from one charged to the
      ! singlet by much more than the rounding of the assembled rows.
      real*8, parameter :: x_tr    = 0.10d0
      real*8, parameter :: n_he_0  = 1.0d8

      ! The rounding floor of a source read out of an assembled row: the
      ! row carries terms of order rowscale, so a difference of two
      ! evaluations of it is meaningful down to eps times that.
      real*8, parameter :: round_tol = 1.0d-13

      real*8 :: k1, k2
      real*8 :: S_he(nsys), S_h(nsys), rowscale(nsys)
      real*8 :: S_he_0, S_h_0, rs0, R1, R2, sc
      integer :: is

      call he_h_cx_rates(T_cell, k1, k2)
      call cx_init()
      call cx_set_cell(T_cell)
      ! The canonical element list with no atoms in any of it: metal_rows
      ! then writes the identity rows of an absent element, no element
      ! contributes an electron or a charge-exchange reaction, and the two
      ! metal systems still run the block they differ from the others by.
      ! The list has to be the canonical one and not a shorter stand-in:
      ! the generic charge-exchange assembly addresses elements by their
      ! canonical index, so it reads all n_melem of the density arrays.
      call set_metal_coeffs(n_melem, zeros(n_melem), zeros(n_melem),      &
                            zeros(n_melem), zeros(n_melem),               &
                            zeros(n_melem), zeros(n_melem),               &
                            zeros(n_melem), tops(n_melem))
      call fill_cell(n_he_0)

      ! The closed form the pair must produce at the reference state, built
      ! here from the paper's two reactions and this state's densities. The
      ! ground singlet is n_He_0 (1 - x_HeII - x_HeIII) in every state
      ! below: the metastable is added on top of it.
      R1 = k1*n_he_0*(1.0d0 - x_heii - x_heiii)*x_hii*n_h   ! He(1^1S) + H+
      R2 = k2*x_heii*n_he_0*(1.0d0 - x_hii)*n_h             ! He+ + H0

      write(*,'(a)') 'DIAGNOSTIC the two Table 4 group B rate coefficients'
      write(*,'(a,es23.15)') '  k(He + H+)  [cm^3 s^-1]  = ', k1
      write(*,'(a,es23.15)') '  k(He+ + H)  [cm^3 s^-1]  = ', k2
      write(*,'(a,es23.15)') '  R1 = k1 n(He 1^1S) n(H+) = ', R1
      write(*,'(a,es23.15)') '  R2 = k2 n(He+) n(H0)     = ', R2
      write(*,'(a,es23.15)') '  k1 n(He 2^3S) n(H+), the term a summed'
      write(*,'(a,es23.15)') '    He I reservoir would add to R1  = ',    &
           k1*x_tr/(1.0d0 - x_tr)*n_he_0*x_hii*n_h

      ! ---- 1. the reservoir is the ground singlet, in every system -----!
      ! The state carries a metastable population of 10 percent of the
      ! helium nuclei, so a pair charged to the summed He I raises R1 by
      ! 1/9 of itself.
      do is = 1, nsys
         call pair_sources(is, x_tr, S_he(is), S_h(is), rowscale(is))
         call check_absolute('pair_He_source_singlet_'//                  &
                             trim(sysname(is)),                           &
                             (S_he(is) - (R2 - R1))/rowscale(is),         &
                             0.0d0, round_tol)
         call check_absolute('pair_H_source_singlet_'//                   &
                             trim(sysname(is)),                           &
                             (S_h(is) - (R1 - R2))/rowscale(is),          &
                             0.0d0, round_tol)
      enddo

      ! ---- 2. the four systems agree on the pair, at one state ---------!
      sc = maxval(rowscale)
      do is = 2, nsys
         call check_absolute('pair_He_source_agrees_'//                   &
                             trim(sysname(is))//'_with_'//                &
                             trim(sysname(1)),                            &
                             (S_he(is) - S_he(1))/sc, 0.0d0, round_tol)
         call check_absolute('pair_H_source_agrees_'//                    &
                             trim(sysname(is))//'_with_'//                &
                             trim(sysname(1)),                            &
                             (S_h(is) - S_h(1))/sc, 0.0d0, round_tol)
      enddo

      ! ---- 3. adding a metastable at fixed singlet moves nothing -------!
      ! The comparison state holds n(H0), n(H+), n(He 1^1S), n(He+) and
      ! n(He2+) at the same absolute values and adds the metastable by
      ! enlarging the helium nucleus density, which is the only way the
      ! fraction unknowns allow it: n(He 1^1S) = (1 - x2 - x3 - x4) n_He is
      ! unchanged when n_He grows by exactly the metastable density. The
      ! He 2^3S + H+ collision has its own rate and is carried nowhere, so
      ! the pair's sources must not notice the metastable at all.
      do is = 1, nsys
         call pair_sources(is, 0.0d0, S_he_0, S_h_0, rs0)
         call check_absolute('pair_He_source_no_metastable_channel_'//    &
                             trim(sysname(is)),                           &
                             (S_he(is) - S_he_0)/max(rowscale(is), rs0),  &
                             0.0d0, round_tol)
         call check_absolute('pair_H_source_no_metastable_channel_'//     &
                             trim(sysname(is)),                           &
                             (S_h(is) - S_h_0)/max(rowscale(is), rs0),    &
                             0.0d0, round_tol)
      enddo

      ! ---- 4. the pair's derivatives in the triplet basis --------------!
      call pair_jacobian_vs_central_difference()

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'he_h_pair_reservoir: ',                     &
              assertion_failures, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'he_h_pair_reservoir: every assertion passed'

      contains

      !--------------!
      function zeros(n) result(z)
      integer, intent(in) :: n
      real*8 :: z(n)
      z = 0.0d0
      end function zeros

      !--------------!
      ! Every element carrying its X++ stage, as the canonical table does
      ! for the elements that reach it; with no atoms present the stage is
      ! an identity row either way.
      function tops(n) result(t)
      integer, intent(in) :: n
      integer :: t(n)
      t = 2
      end function tops

      !--------------!
      ! The cell state every system reads. Only the fields the H/He rows
      ! and the pair use are set; the molecular and metal coefficients stay
      ! at their zero defaults, which is a cell with no such reaction, and
      ! in any case they cancel in the difference that isolates the pair.
      subroutine fill_cell(n_he)
      real*8, intent(in) :: n_he
      ieq_cell%T_K        = T_cell
      ieq_cell%nh         = n_h
      ieq_cell%nhe        = n_he
      ieq_cell%ntot       = n_h + n_he
      ieq_cell%P_HI       = 1.0d-4
      ieq_cell%P_HeI      = 3.0d-5
      ieq_cell%P_HeII     = 1.0d-6
      ieq_cell%P_HeITR    = 2.0d-3
      ieq_cell%rchiiB     = 2.6d-13
      ieq_cell%rcheiiB    = 4.3d-13
      ieq_cell%rcheiiiB   = 1.5d-12
      ieq_cell%rcheiTR    = 1.5d-13
      ieq_cell%a_ion_HI   = 1.0d-11
      ieq_cell%a_ion_HeI  = 1.0d-13
      ieq_cell%a_ion_HeII = 1.0d-15
      ieq_cell%a_ion_HeITR = 1.0d-8
      ieq_cell%A31        = 1.272d-4
      ieq_cell%q13        = 2.1d-8
      ieq_cell%q31a       = 5.0d-10
      ieq_cell%q31b       = 2.0d-9
      ieq_cell%Q31        = 5.0d-10
      ieq_cell%kcx_He0_Hp = k1
      ieq_cell%kcx_Hep_H0 = k2
      end subroutine fill_cell

      !--------------!
      ! The unknown vector of system `is` for a state whose metastable is
      ! xtr of the helium nuclei and whose H and He ion densities are the
      ! reference ones.
      subroutine system_state(is, xtr, neq, x)
      integer, intent(in)  :: is
      real*8,  intent(in)  :: xtr
      integer, intent(out) :: neq
      real*8,  intent(out) :: x(neqmax)
      real*8 :: n_he
      n_he = n_he_0/(1.0d0 - xtr)
      x = 0.0d0
      x(1) = x_hii
      x(2) = x_heii*n_he_0/n_he
      x(3) = x_heiii*n_he_0/n_he
      ! The metal block adds two rows per canonical element.
      select case (is)
      case (1)                         ! the triplet system: x(4) = He 2^3S
         neq  = 4
         x(4) = xtr
      case (2)                         ! metals from row 5
         neq  = 4 + 2*n_melem
         x(4) = xtr
      case (3)                         ! the molecular systems: x(8) = He 2^3S
         neq  = 8
         x(4) = 0.0d0                  ! H2, as a fraction of the H nuclei
         x(8) = xtr
      case default                     ! metals from row 9
         neq  = 8 + 2*n_melem
         x(4) = 0.0d0
         x(8) = xtr
      end select
      end subroutine system_state

      !--------------!
      ! The pair's physical He I and H I number-density sources [cm^-3
      ! s^-1] in system `is`, at the state of system_state, obtained as the
      ! difference of two residual evaluations. rs is the magnitude of the
      ! two rows the pair writes into, which is the scale the difference is
      ! resolved against.
      subroutine pair_sources(is, xtr, S_heI, S_hI, rs)
      integer, intent(in)  :: is
      real*8,  intent(in)  :: xtr
      real*8,  intent(out) :: S_heI, S_hI, rs
      real*8  :: x(neqmax), fon(neqmax), foff(neqmax), d(neqmax)
      integer :: neq
      ! The metastable is added on top of the reference helium: the
      ! nucleus density grows by exactly its density, so the ground
      ! singlet (1 - x2 - x3 - x4) n_He is the reference one either way.
      call system_state(is, xtr, neq, x)
      call fill_cell(n_he_0/(1.0d0 - xtr))
      call residual(is, neq, x, .true.,  fon)
      call residual(is, neq, x, .false., foff)
      d(1:neq) = fon(1:neq) - foff(1:neq)
      rs = max(abs(fon(1)), abs(fon(2)), abs(foff(1)), abs(foff(2)))
      select case (is)
      case (1, 2)            ! row 2 He I-gain positive: u(He,0) = -fvec(2)
         S_heI = d(2)
      case default           ! row 2 the He II balance: u(He,0) = fvec(2)
         S_heI = -d(2)
      end select
      S_hI = -d(1)
      end subroutine pair_sources

      !--------------!
      subroutine residual(is, neq, x, with_pair, fvec)
      integer, intent(in)  :: is, neq
      real*8,  intent(in)  :: x(neqmax)
      logical, intent(in)  :: with_pair
      real*8,  intent(out) :: fvec(neqmax)
      real*8  :: xloc(neqmax), params(60)
      integer :: iflag
      iflag  = 1
      params = 0.0d0
      xloc   = x
      fvec   = 0.0d0
      he_h_charge_exchange = with_pair
      thereis_HeITR   = .true.
      thereis_oxychem = .false.
      mol_inv_turnover(:) = 1.0d0
      select case (is)
      case (1)
         thereis_metals = .false.
         call ion_system_HeH_TR(neq, xloc, fvec, iflag, params(1:40))
      case (2)
         thereis_metals = .true.
         cx_metal_base  = 5
         call ion_system_HeH_TR_metals(neq, xloc, fvec, iflag, params)
         cx_metal_base  = 4
      case (3)
         thereis_metals = .false.
         call ion_system_HeH_mol(neq, xloc, fvec, iflag, params(1:40))
      case (4)
         thereis_metals = .true.
         cx_metal_base  = 9
         call ion_system_HeH_mol_metals(neq, xloc, fvec, iflag, params)
         cx_metal_base  = 4
      end select
      he_h_charge_exchange = .true.
      end subroutine residual

      !--------------!
      ! The pair's derivative with respect to each unknown of the triplet
      ! layout, against a central difference of the rows it is isolated
      ! from. The rate is bilinear in densities linear in the unknowns, so
      ! the central difference is exact to rounding at any step.
      !
      ! Both rows of the triplet basis receive R2 - R1, with
      !   R1 = k1 n(He 1^1S) n(H+),  n(He 1^1S) = (1 - x2 - x3 - x4) n_He
      !   R2 = k2 n(He+)     n(H0),  n(He+) = x2 n_He, n(H0) = (1 - x1) n_H
      ! so the x4 column is nonzero exactly because the reservoir is the
      ! singlet: a metastable made at fixed He+ and He2+ is a singlet
      ! removed, and the pair's rate follows it.
      subroutine pair_jacobian_vs_central_difference()
      real*8  :: x(neqmax), xp(neqmax), fp(neqmax), fm(neqmax)
      real*8  :: an(4), fdk, n_he, n_hii, n_hi, n_heii, n_heiSI, h, rs
      integer :: neq, k, is
      character(len=1) :: kc

      n_he = n_he_0/(1.0d0 - x_tr)
      call fill_cell(n_he)
      call system_state(1, x_tr, neq, x)
      n_hii   = x(1)*n_h
      n_hi    = (1.0d0 - x(1))*n_h
      n_heii  = x(2)*n_he
      n_heiSI = (1.0d0 - x(2) - x(3) - x(4))*n_he

      ! d(R2 - R1)/dx_k, the increment both rows of the triplet basis take
      an(1) = -k2*n_heii*n_h - k1*n_heiSI*n_h
      an(2) =  k2*n_he*n_hi  + k1*n_he*n_hii
      an(3) =  k1*n_he*n_hii
      an(4) =  k1*n_he*n_hii

      ! A step large enough that the difference of two assembled rows is
      ! far above their rounding, and exact all the same.
      h = 1.0d-2
      do is = 1, 2
         call system_state(is, x_tr, neq, x)
         do k = 1, 4
            xp = x;  xp(k) = x(k) + h
            call pair_rows(is, neq, xp, fp, rs)
            xp = x;  xp(k) = x(k) - h
            call pair_rows(is, neq, xp, fm, rs)
            fdk = (fp(2) - fm(2))/(2.0d0*h)
            write(kc,'(i1)') k
            call check_relative('pair_jacobian_column'//kc//'_'//         &
                                trim(sysname(is)), fdk, an(k), 1.0d-5)
         enddo
      enddo
      end subroutine pair_jacobian_vs_central_difference

      !--------------!
      ! The pair's contribution to the rows of system `is`, at x.
      subroutine pair_rows(is, neq, x, d, rs)
      integer, intent(in)  :: is, neq
      real*8,  intent(in)  :: x(neqmax)
      real*8,  intent(out) :: d(neqmax), rs
      real*8  :: fon(neqmax), foff(neqmax)
      call residual(is, neq, x, .true.,  fon)
      call residual(is, neq, x, .false., foff)
      d = 0.0d0
      d(1:neq) = fon(1:neq) - foff(1:neq)
      rs = max(abs(fon(1)), abs(fon(2)))
      ! The triplet basis writes the same increment into both rows, so the
      ! row read below is row 2 and row 1 carries the same number.
      end subroutine pair_rows

      end program he_h_pair_reservoir

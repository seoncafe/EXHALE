      program transported_stage_sources
      ! THE CHARGE-EXCHANGE SOURCE OF A TRANSPORTED IONIZATION STAGE,
      ! against the same reactions read out of the local systems' rows
      ! (PLAN_20260918_rev2 item D7b).
      !
      ! A transported stage carries the stage DENSITY as its unknown, so
      ! its source is the stoichiometric number-density source of that
      ! stage; a local ionization system carries BOUNDARY FLOWS as its
      ! rows.  charge_exchange_stage_sources returns the first and
      ! cx_add_to_fvec assembles the second, both from one reaction table
      ! and one rate vector.  Every row below compares the two on the
      ! PHYSICAL sources they represent, never on raw residual arrays,
      ! whose orientation is exactly what differs between the bases.
      !
      ! WHAT IS ASSERTED
      !   <rxn>_<basis>_stage_sources_match
      !       one active reaction at a time, in each of the three bases the
      !       metal systems assemble: the stage sources equal the sources
      !       the basis's own conversion recovers from the rows, for every
      !       stage of every element.
      !   <rxn>_gross_split, <rxn>_gross_nonnegative
      !       the gross production and loss the stage-source interface
      !       hands out are nonnegative and their difference is the net
      !       source, stage by stage.
      !   <rxn>_stoichiometry, <rxn>_nucleus_balance, <rxn>_charge_source
      !       a single electron transfer moves exactly -R, +R, -R, +R over
      !       four stages, creates and destroys no nucleus, and makes no
      !       net charge, read off the stage sources directly.
      !   production_state_<basis>_stage_sources_match
      !       the same equality with every reaction of the full table
      !       active at once, at densities of the launch region rather
      !       than at one reactant pair.
      !   group_C_absent_without_cx_full, group_D_absent_without_cx_full
      !   group_E_active_without_cx_full, group_E_off_at_zero_scale
      !   group_B_not_in_stage_sources
      !       the active set of the stage sources is the active set of the
      !       residual assembly, reaction group by reaction group, and the
      !       He <-> H pair is counted once because it is in neither.
      !   cell_order_independence, thread_count_independence
      !       eight cells at different temperatures, evaluated in
      !       increasing and then in decreasing order and then on several
      !       OpenMP threads: each cell's stage sources are bitwise the
      !       same, because the rates belong to the cell the caller names
      !       and not to the cell a worker saw before it.
      !
      ! The reactions and their directions are Huang et al. (2023), ApJ
      ! 951, 123, Table 4 (groups A, C and D), read from the journal PDF,
      ! plus the O2+ + H0 electron capture of Barragan et al. (2006) that
      ! the code carries as group E outside Table 4.  They are written out
      ! here from the paper, not read from the module's descriptor tables.
      use species_table,  only: n_melem, iel_C, iel_O, iel_Si
      use charge_exchange, only: cx_full, cx_init, cx_set_cell,           &
                                 cx_metal_base, cx_add_to_fvec,           &
                                 cx_o2p_h_scale, he_h_cx_rates,           &
                                 charge_exchange_stage_sources
      use assertion_report, only: check_absolute, check_positive,         &
                                  assertion_failures
      implicit none

      integer, parameter :: ncode = 12      ! 10 metals, then H = 11, He = 12
      integer, parameter :: cx_H = 11, cx_He = 12
      integer, parameter :: neqmax = 32
      integer, parameter :: nbasis = 3
      character(len=9), parameter :: basis_name(nbasis) =                 &
         [ character(len=9) :: 'standard', 'triplet', 'molecular' ]
      integer, parameter :: basis_base(nbasis) = [ 4, 5, 9 ]
      real*8,  parameter :: basis_sign(nbasis) = [ 1.0d0, -1.0d0, 1.0d0 ]

      ! The reactions, as Table 4 and the group-E source print them.
      ! C1..C6 are metal + He / He+, D1 is metal + metal (no H or He
      ! source at all), A9 is metal + H+, E1 the electron capture.
      integer, parameter :: nreac = 9
      character(len=16), parameter :: rname(nreac) =                      &
         [ character(len=16) :: 'C1_Si_Hep', 'C2_Sip_He', 'C3_C_Hep',     &
           'C4_Cp_He', 'C5_O_Hep', 'C6_Op_He', 'A9_Si_Hp', 'E1_O2p_H',    &
           'D1_C_Sip' ]
      integer :: don_el(nreac), don_stg(nreac)
      integer :: acc_el(nreac), acc_stg(nreac)

      real*8  :: T_cell
      integer :: ir, ib

      T_cell = 1.0d4

      ! donor: element, stage that loses the electron (stage -> stage+1)
      ! acceptor: element, stage that gains it (stage -> stage-1)
      don_el  = [ iel_Si, cx_He,  iel_C,  cx_He, iel_O,  cx_He,           &
                  iel_Si, cx_H,   iel_C ]
      don_stg = [ 0,      0,      0,      0,     0,      0,               &
                  0,      0,      0 ]
      acc_el  = [ cx_He,  iel_Si, cx_He,  iel_C, cx_He,  iel_O,           &
                  cx_H,   iel_O,  iel_Si ]
      acc_stg = [ 1,      1,      1,      1,     1,      1,               &
                  1,      2,      1 ]

      ! ---- one reaction at a time, in each basis ----------------------!
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      do ir = 1, nreac
         call stage_source_stoichiometry(ir)
         do ib = 1, nbasis
            call one_reaction_against_rows(ir, ib)
         enddo
      enddo

      ! ---- the whole set at once, at launch-region densities -----------!
      do ib = 1, nbasis
         call production_state_against_rows(ib)
      enddo

      ! ---- which reaction groups the stage sources carry ---------------!
      call reaction_group_membership()

      ! ---- the cell the rates belong to, and the thread count ----------!
      call cell_order_independence()
      call thread_count_independence()

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'transported_stage_sources: ',               &
              assertion_failures, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'transported_stage_sources: every assertion passed'

      contains

      !--------------!
      ! The physical number-density sources [cm^-3 s^-1] that a basis's
      ! rows represent.  u(el,lower) is the net upward flow across the
      ! boundary between stage `lower` and stage lower+1; the He I <-> He II
      ! boundary is the one the bases disagree about.
      subroutine sources_from_rows(ib, fvec, S)
      integer, intent(in)  :: ib
      real*8,  intent(in)  :: fvec(neqmax)
      real*8,  intent(out) :: S(ncode,0:2)
      real*8  :: u(ncode,0:1)
      integer :: e, base
      base = basis_base(ib)
      u = 0.0d0
      u(cx_H,0) = fvec(1)
      select case (trim(basis_name(ib)))
      case ('standard')      ! row 2 is the He I -> He II boundary flow
         u(cx_He,0) = fvec(2)
         u(cx_He,1) = fvec(3)
      case ('triplet')       ! row 2 is the He I source, row 3 the boundary
         u(cx_He,0) = -fvec(2)
         u(cx_He,1) = fvec(3)
      case ('molecular')     ! row 2 is the He II source, row 3 the He III one
         u(cx_He,1) = fvec(3)
         u(cx_He,0) = fvec(2) + fvec(3)
      end select
      do e = 1, n_melem
         u(e,0) = fvec(base + 2*(e-1))
         u(e,1) = fvec(base + 2*(e-1) + 1)
      enddo
      do e = 1, ncode
         S(e,0) = -u(e,0)
         S(e,1) =  u(e,0) - u(e,1)
         S(e,2) =  u(e,1)
      enddo
      end subroutine sources_from_rows

      !--------------!
      subroutine reactant_pair_densities(ir, nm0, nm1, nm2,               &
                                         n_hi, n_hii, n_hei, n_heii, n_heiii)
      integer, intent(in)  :: ir
      real*8,  intent(out) :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8,  intent(out) :: n_hi, n_hii, n_hei, n_heii, n_heiii
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      n_hi = 0.0d0; n_hii = 0.0d0
      n_hei = 0.0d0; n_heii = 0.0d0; n_heiii = 0.0d0
      call put_one(don_el(ir), don_stg(ir), nm0, nm1, nm2,                &
                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      call put_one(acc_el(ir), acc_stg(ir), nm0, nm1, nm2,                &
                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      end subroutine reactant_pair_densities

      subroutine put_one(el, stg, nm0, nm1, nm2,                          &
                         n_hi, n_hii, n_hei, n_heii, n_heiii)
      integer, intent(in)    :: el, stg
      real*8,  intent(inout) :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8,  intent(inout) :: n_hi, n_hii, n_hei, n_heii, n_heiii
      if (el .eq. cx_H) then
         if (stg .eq. 0) n_hi  = 1.0d0
         if (stg .eq. 1) n_hii = 1.0d0
      else if (el .eq. cx_He) then
         if (stg .eq. 0) n_hei   = 1.0d0
         if (stg .eq. 1) n_heii  = 1.0d0
         if (stg .eq. 2) n_heiii = 1.0d0
      else
         if (stg .eq. 0) nm0(el) = 1.0d0
         if (stg .eq. 1) nm1(el) = 1.0d0
         if (stg .eq. 2) nm2(el) = 1.0d0
      endif
      end subroutine put_one

      !--------------!
      ! The stage sources of one reaction, judged on themselves: the four
      ! stages a single-electron transfer moves, the element balance and
      ! the charge balance, and the gross production and loss.
      subroutine stage_source_stoichiometry(ir)
      integer, intent(in) :: ir
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      real*8  :: pH(0:1), lH(0:1), pHe(0:2), lHe(0:2)
      real*8  :: S(ncode,0:2), Sexp(ncode,0:2)
      real*8  :: R, dev, bal, chg, gdev, gmin
      integer :: e, k

      call reactant_pair_densities(ir, nm0, nm1, nm2,                     &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM,       &
                       p_H = pH, l_H = lH, p_He = pHe, l_He = lHe)
      call gather(sH, sHe, sM, S)

      R = -S(don_el(ir), don_stg(ir))
      call check_positive(trim(rname(ir))//'_stage_rate', R)
      if (R .le. 0.0d0) return

      Sexp = 0.0d0
      Sexp(don_el(ir), don_stg(ir))   = Sexp(don_el(ir), don_stg(ir))   - R
      Sexp(don_el(ir), don_stg(ir)+1) = Sexp(don_el(ir), don_stg(ir)+1) + R
      Sexp(acc_el(ir), acc_stg(ir))   = Sexp(acc_el(ir), acc_stg(ir))   - R
      Sexp(acc_el(ir), acc_stg(ir)-1) = Sexp(acc_el(ir), acc_stg(ir)-1) + R

      dev = 0.0d0
      do e = 1, ncode
         do k = 0, 2
            dev = max(dev, abs(S(e,k) - Sexp(e,k)))
         enddo
      enddo
      call check_absolute(trim(rname(ir))//'_stoichiometry', dev/R,       &
                          0.0d0, 1.0d-13)

      bal = 0.0d0
      do e = 1, ncode
         bal = max(bal, abs(S(e,0) + S(e,1) + S(e,2)))
      enddo
      call check_absolute(trim(rname(ir))//'_nucleus_balance', bal/R,     &
                          0.0d0, 1.0d-13)

      chg = 0.0d0
      do e = 1, ncode
         chg = chg + S(e,1) + 2.0d0*S(e,2)
      enddo
      call check_absolute(trim(rname(ir))//'_charge_source', chg/R,       &
                          0.0d0, 1.0d-13)

      ! The gross production and loss of the two transported elements.
      gdev = 0.0d0
      gmin = 0.0d0
      do k = 0, 1
         gdev = max(gdev, abs((pH(k) - lH(k)) - sH(k)))
         gmin = min(gmin, pH(k), lH(k))
      enddo
      do k = 0, 2
         gdev = max(gdev, abs((pHe(k) - lHe(k)) - sHe(k)))
         gmin = min(gmin, pHe(k), lHe(k))
      enddo
      call check_absolute(trim(rname(ir))//'_gross_split', gdev/R,        &
                          0.0d0, 1.0d-13)
      call check_absolute(trim(rname(ir))//'_gross_nonnegative', gmin,    &
                          0.0d0, 0.0d0)
      end subroutine stage_source_stoichiometry

      !--------------!
      subroutine gather(sH, sHe, sM, S)
      real*8, intent(in)  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      real*8, intent(out) :: S(ncode,0:2)
      integer :: e, k
      S = 0.0d0
      S(cx_H,0) = sH(0)
      S(cx_H,1) = sH(1)
      do k = 0, 2
         S(cx_He,k) = sHe(k)
      enddo
      do e = 1, n_melem
         do k = 0, 2
            S(e,k) = sM(e,k)
         enddo
      enddo
      end subroutine gather

      !--------------!
      ! One reaction, the stage sources against the rows of one basis.
      subroutine one_reaction_against_rows(ir, ib)
      integer, intent(in) :: ir, ib
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: fvec(neqmax), Srow(ncode,0:2), Sstg(ncode,0:2)
      real*8  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      real*8  :: R, dev
      integer :: e, k

      call reactant_pair_densities(ir, nm0, nm1, nm2,                     &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      cx_metal_base = basis_base(ib)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          basis_sign(ib), T_cell)
      call sources_from_rows(ib, fvec, Srow)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, Sstg)

      R = -Sstg(don_el(ir), don_stg(ir))
      if (R .le. 0.0d0) return
      dev = 0.0d0
      do e = 1, ncode
         do k = 0, 2
            dev = max(dev, abs(Srow(e,k) - Sstg(e,k)))
         enddo
      enddo
      call check_absolute(trim(rname(ir))//'_'//trim(basis_name(ib))//    &
                          '_stage_sources_match', dev/R, 0.0d0, 1.0d-13)
      end subroutine one_reaction_against_rows

      !--------------!
      ! Every reaction of the full table active at once, at densities of
      ! the launch region: hydrogen and helium partly ionized, every metal
      ! present and spread over its stages.  With 64 reactions summing into
      ! the same rows the comparison is of two sums whose terms are added
      ! in different orders, so the tolerance is the relative one of the
      ! largest stage source and not zero.
      subroutine production_state_against_rows(ib)
      integer, intent(in) :: ib
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: fvec(neqmax), Srow(ncode,0:2), Sstg(ncode,0:2)
      real*8  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      real*8  :: dev, scale, abund(n_melem), n_X
      integer :: e, k

      ! Solar abundances relative to hydrogen, the order of magnitude of
      ! the mixture the metal runs carry (C, O, N, Mg, Si, Ca, Na, K, S,
      ! Fe), used here only to make no element negligible.
      abund = [ 2.7d-4, 4.9d-4, 6.8d-5, 4.0d-5, 3.2d-5,                   &
                2.2d-6, 1.7d-6, 1.2d-7, 1.3d-5, 3.2d-5 ]
      n_hi    = 7.0d8
      n_hii   = 3.0d8
      n_hei   = 8.0d7
      n_heii  = 1.6d7
      n_heiii = 4.0d6
      do e = 1, n_melem
         n_X = 1.0d9*abund(e)
         nm0(e) = 0.55d0*n_X
         nm1(e) = 0.35d0*n_X
         nm2(e) = 0.10d0*n_X
      enddo

      cx_metal_base = basis_base(ib)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          basis_sign(ib), T_cell)
      call sources_from_rows(ib, fvec, Srow)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, Sstg)

      dev = 0.0d0;  scale = 0.0d0
      do e = 1, ncode
         do k = 0, 2
            dev   = max(dev, abs(Srow(e,k) - Sstg(e,k)))
            scale = max(scale, abs(Sstg(e,k)))
         enddo
      enddo
      call check_positive('production_state_'//trim(basis_name(ib))//     &
                          '_scale', scale)
      call check_absolute('production_state_'//trim(basis_name(ib))//     &
                          '_stage_sources_match', dev/scale, 0.0d0,       &
                          1.0d-13)
      end subroutine production_state_against_rows

      !--------------!
      ! The active set of the stage sources is the active set of the
      ! residual assembly, group by group.
      subroutine reaction_group_membership()
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      real*8  :: S(ncode,0:2), k1, k2
      real*8  :: scale_save

      ! Group C (metal + He) and group D (metal + metal) are off unless
      ! cx_full is set.
      cx_full = .false.
      call cx_init()
      call cx_set_cell(T_cell)
      call reactant_pair_densities(1, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_absolute('group_C_absent_without_cx_full',               &
                          maxval(abs(S)), 0.0d0, 0.0d0)
      call reactant_pair_densities(9, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_absolute('group_D_absent_without_cx_full',               &
                          maxval(abs(S)), 0.0d0, 0.0d0)

      ! Group E is not in Table 4 and is gated by its own scale, not by
      ! cx_full: on at the default scale with cx_full off, and gone when
      ! the scale is zeroed.
      call reactant_pair_densities(8, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_positive('group_E_active_without_cx_full', -S(cx_H,0))
      scale_save = cx_o2p_h_scale
      cx_o2p_h_scale = 0.0d0
      call cx_init()
      call cx_set_cell(T_cell)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_hei, n_heii, n_heiii, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_absolute('group_E_off_at_zero_scale', maxval(abs(S)),    &
                          0.0d0, 0.0d0)
      cx_o2p_h_scale = scale_save

      ! Group B, the He <-> H pair, is in neither set: it is applied by
      ! he_h_cx_fvec in every system with helium, so the stage sources
      ! return nothing for it and the pair is counted once.
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      call charge_exchange_stage_sources(nm0, nm1, nm2, 0.0d0, 1.0d0,     &
                       1.0d0, 0.0d0, 0.0d0, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_absolute('group_B_He0_Hp_not_in_stage_sources',          &
                          maxval(abs(S)), 0.0d0, 0.0d0)
      call charge_exchange_stage_sources(nm0, nm1, nm2, 1.0d0, 0.0d0,     &
                       0.0d0, 1.0d0, 0.0d0, T_cell, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      call check_absolute('group_B_Hep_H0_not_in_stage_sources',          &
                          maxval(abs(S)), 0.0d0, 0.0d0)
      ! The pair itself is not zero at that state, so the two rows above
      ! say that the generic set omits it and not that the state is empty.
      call he_h_cx_rates(T_cell, k1, k2)
      call check_positive('group_B_rate_is_nonzero_at_this_state', k1)
      end subroutine reaction_group_membership

      !--------------!
      ! The rates belong to the cell the caller names.  Eight cells are
      ! evaluated in increasing and then in decreasing temperature order;
      ! each cell's stage sources must be bitwise the same both times.
      subroutine cell_order_independence()
      integer, parameter :: ncell = 8
      real*8  :: Tc(ncell)
      real*8  :: up(ncode,0:2,ncell), dn(ncode,0:2,ncell)
      real*8  :: S(ncode,0:2), dev, scale
      integer :: j, e, k
      cx_full = .true.
      call cx_init()
      do j = 1, ncell
         Tc(j) = 2.0d3 + 2.0d3*dble(j-1)
      enddo
      do j = 1, ncell
         call sources_at_cell(Tc(j), S)
         up(:,:,j) = S
      enddo
      do j = ncell, 1, -1
         call sources_at_cell(Tc(j), S)
         dn(:,:,j) = S
      enddo
      dev = 0.0d0;  scale = 0.0d0
      do j = 1, ncell
         do e = 1, ncode
            do k = 0, 2
               dev   = max(dev, abs(up(e,k,j) - dn(e,k,j)))
               scale = max(scale, abs(up(e,k,j)))
            enddo
         enddo
      enddo
      call check_positive('cell_order_scale', scale)
      call check_absolute('cell_order_independence', dev, 0.0d0, 0.0d0)
      end subroutine cell_order_independence

      !--------------!
      ! The same eight cells shared out over several OpenMP threads.  The
      ! rate coefficients are thread-local and each cell loads its own, so
      ! the source of a cell cannot depend on how many workers there are or
      ! on which of them took it.
      subroutine thread_count_independence()
      integer, parameter :: ncell = 8
      real*8  :: Tc(ncell)
      real*8  :: one(ncode,0:2,ncell), many(ncode,0:2,ncell)
      real*8  :: S(ncode,0:2), dev, scale
      integer :: j, e, k
      cx_full = .true.
      call cx_init()
      do j = 1, ncell
         Tc(j) = 2.0d3 + 2.0d3*dble(j-1)
      enddo
      do j = 1, ncell
         call sources_at_cell(Tc(j), S)
         one(:,:,j) = S
      enddo
      many = 0.0d0
      !$omp parallel do default(shared) private(j,S) schedule(dynamic,1)  &
      !$omp   num_threads(4)
      do j = 1, ncell
         call sources_at_cell(Tc(j), S)
         many(:,:,j) = S
      enddo
      !$omp end parallel do
      dev = 0.0d0;  scale = 0.0d0
      do j = 1, ncell
         do e = 1, ncode
            do k = 0, 2
               dev   = max(dev, abs(one(e,k,j) - many(e,k,j)))
               scale = max(scale, abs(one(e,k,j)))
            enddo
         enddo
      enddo
      call check_positive('thread_count_scale', scale)
      call check_absolute('thread_count_independence', dev, 0.0d0, 0.0d0)
      end subroutine thread_count_independence

      ! One cell: load its rate coefficients on this thread, then evaluate
      ! the stage sources of a fixed composition at its temperature.
      subroutine sources_at_cell(T, S)
      real*8, intent(in)  :: T
      real*8, intent(out) :: S(ncode,0:2)
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: sH(0:1), sHe(0:2), sM(n_melem,0:2)
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      nm0(iel_Si) = 1.0d3
      nm0(iel_C)  = 2.0d3
      nm1(iel_O)  = 5.0d2
      call cx_set_cell(T)
      call charge_exchange_stage_sources(nm0, nm1, nm2, 1.0d6, 5.0d5,     &
                       2.0d5, 1.0d4, 1.0d2, T, sH, sHe, sM)
      call gather(sH, sHe, sM, S)
      end subroutine sources_at_cell

      end program transported_stage_sources

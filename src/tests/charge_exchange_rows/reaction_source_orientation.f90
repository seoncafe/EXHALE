      program reaction_source_orientation
      ! One charge-exchange reaction at a time, in each of the three
      ! equation bases the metal systems assemble, judged on the PHYSICAL
      ! number-density sources the basis's own conversion recovers from the
      ! rows -- never on the raw rows, whose orientation is what differs.
      !
      ! WHAT IS ASSERTED, for every reaction and every basis
      !   stoichiometry   the four stages a single-electron transfer moves
      !                   carry exactly -R, +R, -R, +R (donor lower, donor
      !                   upper, acceptor upper, acceptor lower) and every
      !                   other stage of every element carries nothing;
      !   nucleus balance sum over the stages of an element is zero, so no
      !                   element is created or destroyed;
      !   charge          sum over all species of stage times source is
      !                   zero: an exchange moves one electron and makes no
      !                   net charge.
      ! The reactions and their directions are the ones of Huang et al.
      ! (2023), ApJ 951, 123, Table 4 (groups A and C), read from the
      ! journal PDF, plus the O2+ + H0 electron capture of Barragan et al.
      ! (2006) that the code carries as group E outside Table 4. They are
      ! written out here from the paper, not read from the module's own
      ! descriptor tables, so a sign error in those tables is visible.
      !
      ! THE BASES
      !   standard  ion_system_HeH_metals: row 1 the H I -> H II boundary
      !             flow, row 2 the He I -> He II boundary flow, row 3 the
      !             He II -> He III boundary flow, metals from row 4. He row
      !             ionization positive, he_row_sign = +1.
      !   triplet   ion_system_HeH_TR_metals: row 2 is the SUMMED He I
      !             balance of heh_tr_rows, written He I-gain positive, so
      !             he_row_sign = -1 and the He I -> He II boundary flow is
      !             minus that row. Metals from row 5.
      !   molecular ion_system_HeH_mol_metals: row 2 is the He II balance
      !             and row 3 the He III balance of mol_heh_rows, both
      !             production positive, so he_row_sign = +1. Metals from
      !             row 9 here (metal_row_base with the oxygen carriers on).
      !
      ! ENTRY-TEXT BEHAVIOR, printed as DIAGNOSTIC and not asserted: before
      ! the row-orientation argument existed, the triplet system applied the
      ! generic assembly in the standard orientation. The diagnostic repeats
      ! that by asking for he_row_sign = +1 in the triplet basis, where it
      ! gives the He II source the wrong sign and a nonzero charge source.
      use species_table,  only: n_melem, iel_C, iel_O, iel_Si
      use charge_exchange, only: cx_full, cx_init, cx_set_cell,           &
                                 cx_metal_base, cx_add_to_fvec,           &
                                 cx_add_to_jac, cx_o2p_h_scale,           &
                                 he_h_cx_rates, he_h_cx_fvec
      use assertion_report, only: check_absolute, check_positive,         &
                                  check_relative, assertion_failures
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
      integer, parameter :: nreac = 8
      character(len=16), parameter :: rname(nreac) =                      &
         [ character(len=16) :: 'C1_Si_Hep', 'C2_Sip_He', 'C3_C_Hep',     &
           'C4_Cp_He', 'C5_O_Hep', 'C6_Op_He', 'A9_Si_Hp', 'E1_O2p_H' ]
      integer :: don_el(nreac), don_stg(nreac)
      integer :: acc_el(nreac), acc_stg(nreac)
      logical :: needs_full(nreac)

      real*8  :: T_cell
      integer :: ir, ib

      T_cell = 1.0d4

      ! donor: element, stage that loses the electron (stage -> stage+1)
      ! acceptor: element, stage that gains it (stage -> stage-1)
      don_el  = [ iel_Si, cx_He, iel_C,  cx_He, iel_O,  cx_He, iel_Si, cx_H ]
      don_stg = [ 0,      0,     0,      0,     0,      0,     0,      0    ]
      acc_el  = [ cx_He,  iel_Si, cx_He, iel_C, cx_He,  iel_O, cx_H,   iel_O]
      acc_stg = [ 1,      1,     1,      1,     1,      1,     1,      2    ]
      ! groups C1..C6 are Table 4 metal + He / He+ and need cx_full; A9 is
      ! metal + H+, active by default; E1 is not in Table 4 and is gated by
      ! cx_o2p_h_scale alone.
      needs_full = [ .true., .true., .true., .true., .true., .true.,      &
                     .false., .false. ]

      ! ---- the full Table-4 set, one reaction at a time, in each basis --!
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      do ir = 1, nreac
         do ib = 1, nbasis
            call one_reaction(ir, ib)
         enddo
      enddo

      ! ---- the entry-text orientation in the triplet basis ------------!
      call entry_text_triplet_diagnostic()

      ! ---- group C is off, and group E on, without cx_full -------------!
      call default_set_membership()

      ! ---- the He <-> H pair is counted once ---------------------------!
      call group_b_counted_once()

      ! ---- zero abundances ---------------------------------------------!
      call zero_abundance_limit()

      ! ---- analytic Jacobian against a central difference ---------------!
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      do ib = 1, nbasis
         call jacobian_vs_central_difference(ib)
      enddo

      ! ---- the stored cell is the caller's, whatever the cell order -----!
      call cell_order_independence()

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'reaction_source_orientation: ',             &
              assertion_failures, ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'reaction_source_orientation: every assertion passed'

      contains

      !--------------!
      ! The physical number-density sources [cm^-3 s^-1] that a basis's
      ! rows represent. u(el,lower) is the net upward flow across the
      ! boundary between stage `lower` and stage lower+1, which every basis
      ! writes the same way for H and for the metals; the He I <-> He II
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
      ! Densities with the two reactants of reaction ir at 1 cm^-3 and every
      ! other species absent, so that reaction is the only event with a
      ! nonzero rate.
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
      subroutine one_reaction(ir, ib)
      integer, intent(in) :: ir, ib
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: fvec(neqmax), S(ncode,0:2), Sexp(ncode,0:2)
      real*8  :: R, dev, bal, chg
      integer :: e, k
      character(len=64) :: tag

      call reactant_pair_densities(ir, nm0, nm1, nm2,                     &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      cx_metal_base = basis_base(ib)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          basis_sign(ib), T_cell)
      call sources_from_rows(ib, fvec, S)

      ! The rate the assembly produced, read off the donor's own stage.
      R = -S(don_el(ir), don_stg(ir))
      tag = trim(rname(ir))//'_'//trim(basis_name(ib))
      if (ib .eq. 1) call check_positive(trim(rname(ir))//'_rate', R)
      if (R .le. 0.0d0) return

      ! One electron moves: the donor rises a stage, the acceptor falls one.
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
      call check_absolute(trim(tag)//'_stoichiometry', dev/R, 0.0d0, 1.0d-13)

      bal = 0.0d0
      do e = 1, ncode
         bal = max(bal, abs(S(e,0) + S(e,1) + S(e,2)))
      enddo
      call check_absolute(trim(tag)//'_nucleus_balance', bal/R, 0.0d0,    &
                          1.0d-13)

      chg = 0.0d0
      do e = 1, ncode
         chg = chg + S(e,1) + 2.0d0*S(e,2)
      enddo
      call check_absolute(trim(tag)//'_charge_source', chg/R, 0.0d0,      &
                          1.0d-13)
      end subroutine one_reaction

      !--------------!
      ! What the triplet system produced before the orientation argument
      ! existed: the generic assembly in the standard orientation, applied
      ! to rows whose He I balance is He I-gain positive.
      subroutine entry_text_triplet_diagnostic()
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: fvec(neqmax), S(ncode,0:2), R, chg
      integer :: e

      call reactant_pair_densities(1, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      cx_metal_base = 5
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          1.0d0, T_cell)
      call sources_from_rows(2, fvec, S)
      R = -S(iel_Si,0)
      chg = 0.0d0
      do e = 1, ncode
         chg = chg + S(e,1) + 2.0d0*S(e,2)
      enddo
      write(*,'(a)') 'DIAGNOSTIC entry-text orientation, Si I + He II ->'
      write(*,'(a)') 'DIAGNOSTIC   Si II + He I read in the triplet basis'
      write(*,'(a,es23.15)') '  Si II source   ', S(iel_Si,1)
      write(*,'(a,es23.15)') '  He II source   ', S(cx_He,1)
      write(*,'(a,es23.15)') '  He I source    ', S(cx_He,0)
      write(*,'(a,es23.15)') '  charge source  ', chg
      write(*,'(a,es23.15)') '  required He II source ', -R
      write(*,'(a,es23.15)') '  required charge source ', 0.0d0
      end subroutine entry_text_triplet_diagnostic

      !--------------!
      ! With cx_full off, the metal + He group is not assembled at all,
      ! while the group-E electron capture, which is not in Table 4 and is
      ! gated by its own scale factor, still is.
      subroutine default_set_membership()
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: fvec(neqmax), S(ncode,0:2)
      cx_full = .false.
      call cx_init()
      call cx_set_cell(T_cell)
      cx_metal_base = 4

      call reactant_pair_densities(1, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          1.0d0, T_cell)
      call sources_from_rows(1, fvec, S)
      call check_absolute('group_C_absent_without_cx_full',               &
                          maxval(abs(S)), 0.0d0, 0.0d0)

      call reactant_pair_densities(8, nm0, nm1, nm2,                      &
                                   n_hi, n_hii, n_hei, n_heii, n_heiii)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          n_hi, n_hii, n_hei, n_heii, n_heiii,            &
                          1.0d0, T_cell)
      call sources_from_rows(1, fvec, S)
      call check_positive('group_E_active_without_cx_full', -S(cx_H,0))
      call check_absolute('group_E_scale_is_its_own_gate',                &
                          cx_o2p_h_scale, 1.0d0, 0.0d0)
      end subroutine default_set_membership

      !--------------!
      ! The He <-> H pair is excluded from the generic set, so the generic
      ! assembly returns nothing for it in either direction and the pair
      ! reaches a row only through he_h_cx_fvec: one reaction, counted once.
      subroutine group_b_counted_once()
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: fvec(neqmax), k1, k2, expect
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      cx_metal_base = 4
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0

      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0,              &
                          1.0d0, T_cell)
      call check_absolute('group_B_He0_Hp_not_in_generic_set',            &
                          maxval(abs(fvec)), 0.0d0, 0.0d0)
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          1.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,              &
                          1.0d0, T_cell)
      call check_absolute('group_B_Hep_H0_not_in_generic_set',            &
                          maxval(abs(fvec)), 0.0d0, 0.0d0)

      ! The dedicated pair, He0 + H+ alone: it ionizes He and neutralizes H.
      call he_h_cx_rates(T_cell, k1, k2)
      fvec = 0.0d0
      call he_h_cx_fvec(fvec, k1, k2, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0)
      expect = k1
      call check_relative('group_B_He0_Hp_ionizes_He', fvec(2), expect,   &
                          1.0d-14)
      call check_relative('group_B_He0_Hp_neutralizes_H', fvec(1),        &
                          -expect, 1.0d-14)
      ! The same reaction in the He I-gain basis of the triplet systems.
      fvec = 0.0d0
      call he_h_cx_fvec(fvec, k1, k2, 0.0d0, 1.0d0, 1.0d0, 0.0d0, -1.0d0)
      call check_relative('group_B_He0_Hp_triplet_basis', fvec(2),        &
                          -expect, 1.0d-14)
      end subroutine group_b_counted_once

      !--------------!
      subroutine zero_abundance_limit()
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: fvec(neqmax)
      integer :: ib
      cx_full = .true.
      call cx_init()
      call cx_set_cell(T_cell)
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      do ib = 1, nbasis
         cx_metal_base = basis_base(ib)
         fvec = 0.0d0
         call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                 &
                             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,           &
                             basis_sign(ib), T_cell)
         call check_absolute('zero_abundance_'//trim(basis_name(ib)),     &
                             maxval(abs(fvec)), 0.0d0, 0.0d0)
      enddo
      end subroutine zero_abundance_limit

      !--------------!
      ! The analytic Jacobian is the derivative of the rows the residual
      ! assembles, in whatever orientation those rows carry. The rate is
      ! bilinear in densities that are linear in the unknowns, so a central
      ! difference of the residual is exact to rounding and the comparison
      ! is a strict one.
      subroutine jacobian_vs_central_difference(ib)
      integer, intent(in) :: ib
      integer, parameter  :: nx = neqmax
      real*8  :: x(nx), xp(nx), fp(neqmax), fm(neqmax)
      real*8  :: fjac(neqmax,neqmax), fd(neqmax,neqmax)
      real*8  :: n_X(n_melem), n_h, n_he, h, dev, scale
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      integer :: e, k, i, ix, nxu

      n_h  = 1.0d8
      n_he = 1.0d7
      do e = 1, n_melem
         n_X(e) = 1.0d4*dble(e)
      enddo
      cx_metal_base = basis_base(ib)
      nxu = basis_base(ib) + 2*n_melem - 1

      x = 0.0d0
      x(1) = 0.30d0                  ! H II per hydrogen nucleus
      x(2) = 0.20d0                  ! He II per helium nucleus
      x(3) = 0.05d0                  ! He III per helium nucleus
      do e = 1, n_melem
         ix = basis_base(ib) + 2*(e-1)
         x(ix)   = 0.25d0
         x(ix+1) = 0.10d0
      enddo

      call densities_from_x(ib, x, n_X, n_h, n_he, nm0, nm1, nm2,         &
                            n_hi, n_hii, n_hei, n_heii, n_heiii)
      fjac = 0.0d0
      call cx_add_to_jac(neqmax, fjac, nm0, nm1, nm2,                     &
                         n_hi, n_hii, n_hei, n_heii, n_heiii,             &
                         n_X, n_h, n_he, basis_sign(ib), T_cell)

      fd = 0.0d0
      h  = 1.0d-4
      do k = 1, nxu
         xp = x;  xp(k) = x(k) + h
         call densities_from_x(ib, xp, n_X, n_h, n_he, nm0, nm1, nm2,     &
                               n_hi, n_hii, n_hei, n_heii, n_heiii)
         fp = 0.0d0
         call cx_add_to_fvec(neqmax, fp, nm0, nm1, nm2,                   &
                             n_hi, n_hii, n_hei, n_heii, n_heiii,         &
                             basis_sign(ib), T_cell)
         xp = x;  xp(k) = x(k) - h
         call densities_from_x(ib, xp, n_X, n_h, n_he, nm0, nm1, nm2,     &
                               n_hi, n_hii, n_hei, n_heii, n_heiii)
         fm = 0.0d0
         call cx_add_to_fvec(neqmax, fm, nm0, nm1, nm2,                   &
                             n_hi, n_hii, n_hei, n_heii, n_heiii,         &
                             basis_sign(ib), T_cell)
         do i = 1, nxu
            fd(i,k) = (fp(i) - fm(i))/(2.0d0*h)
         enddo
      enddo

      dev = 0.0d0;  scale = 0.0d0
      do k = 1, nxu
         do i = 1, nxu
            dev   = max(dev, abs(fjac(i,k) - fd(i,k)))
            scale = max(scale, abs(fd(i,k)))
         enddo
      enddo
      call check_positive('jacobian_scale_'//trim(basis_name(ib)), scale)
      call check_absolute('jacobian_vs_fd_'//trim(basis_name(ib)),        &
                          dev/scale, 0.0d0, 1.0d-9)
      end subroutine jacobian_vs_central_difference

      subroutine densities_from_x(ib, x, n_X, n_h, n_he, nm0, nm1, nm2,   &
                                  n_hi, n_hii, n_hei, n_heii, n_heiii)
      ! The unknown layout of cx_dens_lin: x1 = n(H II)/n_H,
      ! x2 = n(He II)/n_He, x3 = n(He III)/n_He, and for metal e the two
      ! ionized fractions at cx_metal_base + 2*(e-1).
      integer, intent(in)  :: ib
      real*8,  intent(in)  :: x(neqmax), n_X(n_melem), n_h, n_he
      real*8,  intent(out) :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8,  intent(out) :: n_hi, n_hii, n_hei, n_heii, n_heiii
      integer :: e, ix
      n_hii   = n_h*x(1)
      n_hi    = n_h - n_hii
      n_heii  = n_he*x(2)
      n_heiii = n_he*x(3)
      n_hei   = n_he - n_heii - n_heiii
      do e = 1, n_melem
         ix = basis_base(ib) + 2*(e-1)
         nm1(e) = n_X(e)*x(ix)
         nm2(e) = n_X(e)*x(ix+1)
         nm0(e) = n_X(e) - nm1(e) - nm2(e)
      enddo
      end subroutine densities_from_x

      !--------------!
      ! The rates belong to the cell they were filled for, not to the cell
      ! the thread happened to see before it. Eight cells are assembled in
      ! increasing and then in decreasing temperature order, and each cell's
      ! rows must be the same both times.
      subroutine cell_order_independence()
      integer, parameter :: ncell = 8
      real*8  :: Tc(ncell), f_up(neqmax,ncell), f_dn(neqmax,ncell)
      real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8  :: fvec(neqmax), dev, scale
      integer :: j, i
      cx_full = .true.
      call cx_init()
      nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
      nm0(iel_Si) = 1.0d0
      do j = 1, ncell
         Tc(j) = 2.0d3 + 2.0d3*dble(j-1)
      enddo
      do j = 1, ncell
         call assemble_cell(Tc(j), nm0, nm1, nm2, fvec)
         f_up(:,j) = fvec
      enddo
      do j = ncell, 1, -1
         call assemble_cell(Tc(j), nm0, nm1, nm2, fvec)
         f_dn(:,j) = fvec
      enddo
      dev = 0.0d0;  scale = 0.0d0
      do j = 1, ncell
         do i = 1, neqmax
            dev   = max(dev, abs(f_up(i,j) - f_dn(i,j)))
            scale = max(scale, abs(f_up(i,j)))
         enddo
      enddo
      call check_positive('cell_order_scale', scale)
      call check_absolute('cell_order_independence', dev, 0.0d0, 0.0d0)
      end subroutine cell_order_independence

      subroutine assemble_cell(T, nm0, nm1, nm2, fvec)
      real*8, intent(in)  :: T, nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real*8, intent(out) :: fvec(neqmax)
      call cx_set_cell(T)
      cx_metal_base = 5
      fvec = 0.0d0
      call cx_add_to_fvec(neqmax, fvec, nm0, nm1, nm2,                    &
                          0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,              &
                          -1.0d0, T)
      end subroutine assemble_cell

      end program reaction_source_orientation

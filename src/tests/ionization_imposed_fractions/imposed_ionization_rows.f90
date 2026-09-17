      program imposed_ionization_rows
      ! The substitution that hands an ionization residual the fraction the
      ! flow carries, measured in the seven systems that can reach the rows
      ! it writes.
      !
      ! Production routine exercised:
      ! impose_transported_ionization_fractions
      ! (src/modules/nonlinear_system_solver/ion_residual_core.f90), called
      ! by ion_system_H, ion_system_HeH, ion_system_HeH_TR,
      ! ion_system_HeH_metals, ion_system_HeH_TR_metals, ion_system_HeH_mol
      ! and ion_system_HeH_mol_metals.
      !
      ! REFERENCE. Where an ionization stage is transported, its fraction in
      ! a cell is not a root of that cell's photoionization balance, and the
      ! balance row that would have computed it is replaced by the carried
      ! value:
      !
      !       fvec(1) = x(1) - x_hp_fix ,   x(1) = n(H II)/n(H nuclei).
      !
      ! Three statements follow, and each is measured here for every system:
      !
      !   (a) with the fraction imposed, row 1 IS that expression, to the
      !       last bit -- an identity row whose Jacobian row is the identity;
      !   (b) with nothing imposed, the residual does not depend on the
      !       carried value at all, so a cell the operator does not own is
      !       solved by the same rows as before the substitution existed;
      !   (c) the imposed row is one number, the same in every system, since
      !       the substitution is written once and shared.
      !
      ! Rows 2 and 3 are He II and He III per He nucleus. The cell state
      ! states no imposed helium fraction, so the routine writes neither,
      ! and (a) covers row 1 alone.
      !
      ! The tolerances are exact: an imposed row is an assignment, not an
      ! arithmetic result, and a row the substitution does not touch is the
      ! same floating-point number it was.
      use global_parameters,  only: thereis_HeITR, thereis_metals,        &
                                    thereis_oxychem
      use ion_cell_state,     only: ieq_cell
      use ion_residual_core,  only: impose_transported_ionization_fractions
      use System_H,           only: ion_system_H
      use System_HeH,         only: ion_system_HeH
      use System_HeH_TR,      only: ion_system_HeH_TR
      use System_HeH_metals,  only: ion_system_HeH_metals, set_metal_coeffs
      use System_HeH_TR_metals, only: ion_system_HeH_TR_metals
      use System_HeH_mol,     only: ion_system_HeH_mol, set_mol_coeffs,   &
                                    set_mol_turnover_rates
      use System_HeH_mol_metals, only: ion_system_HeH_mol_metals,         &
                                    set_mol_metal_turnover_rates
      use mol_rates,          only: h2_thermochemistry_init
      use species_table,      only: n_melem
      use assertion_report
      implicit none

      ! The metal element list is the canonical one of species_table: the
      ! charge-exchange assembly indexes the metal densities by the same
      ! element codes, so a shorter list than n_melem is not a smaller
      ! problem but an out-of-bounds read. Two elements carry an abundance
      ! and the rest carry none, which is the production shape of a trace
      ! metal run.
      integer, parameter :: nelem = n_melem
      integer, parameter :: nmax  = 7 + 2*nelem
      ! The carried H+ fraction, and a second value of it: (b) compares the
      ! residual at the two with nothing imposed, where neither may be read.
      real*8, parameter  :: x_fix_a = 0.77d0
      real*8, parameter  :: x_fix_b = 0.11d0

      integer, parameter :: nsys = 7
      character(len=24)  :: sysname(nsys)
      integer :: nrow(nsys)
      real*8  :: row1_imposed(nsys)
      real*8  :: x(nmax), fvec_off(nmax), fvec_alt(nmax), fvec_on(nmax)
      real*8  :: mtot(nelem), mg0(nelem), mg1(nelem)
      real*8  :: mb0(nelem), mb1(nelem), ma1(nelem), ma2(nelem)
      integer :: mtop(nelem)
      real*8  :: probe(nmax), probe_in(nmax)
      real*8  :: dmax, constraint
      integer :: is, i

      sysname(1) = 'H'
      sysname(2) = 'HeH'
      sysname(3) = 'HeH_TR'
      sysname(4) = 'HeH_metals'
      sysname(5) = 'HeH_TR_metals'
      sysname(6) = 'HeH_mol'
      sysname(7) = 'HeH_mol_metals'
      ! Unknowns of each system: H alone; H/He; H/He + the metastable;
      ! H/He + 2 per metal element; the same with the metastable; the
      ! molecular network (H+, He+, He++, H2, H2+, H3+, HeH+); and that
      ! network with the metals above it.
      nrow(1) = 1
      nrow(2) = 3
      nrow(3) = 4
      nrow(4) = 3 + 2*nelem
      nrow(5) = 4 + 2*nelem
      nrow(6) = 7
      nrow(7) = 7 + 2*nelem

      thereis_HeITR   = .false.
      thereis_metals  = .true.
      thereis_oxychem = .false.

      call set_cell_rates()

      ! Two elements present at trace abundance, three stages each; the
      ! rest absent, which metal_rows writes as identity rows.
      mtot(:) = 0.0d0
      mg0(:)  = 3.0d-6
      mg1(:)  = 4.0d-7
      mb0(:)  = 2.0d-11
      mb1(:)  = 5.0d-12
      ma1(:)  = 3.0d-13
      ma2(:)  = 9.0d-13
      mtop(:) = 2
      mtot(1) = 1.0d4
      mtot(4) = 3.0d3
      call set_metal_coeffs(nelem, mtot, mg0, mg1, mb0, mb1, ma1, ma2, mtop)

      ! The molecular rate coefficients and the turnover scale the molecular
      ! rows are divided by, both as ioniz_eq sets them for a cell. The H2
      ! equilibrium-constant table is built serially first, as the sweep
      ! does before it enters its parallel region.
      call h2_thermochemistry_init()
      call set_mol_coeffs(ieq_cell%T_K, ieq_cell%ntot)
      call set_mol_turnover_rates(3.0d7, 1.0d0)
      call set_mol_metal_turnover_rates(3.0d7, 1.0d0)

      ! A composition inside every element budget: the H nuclei held by H+,
      ! H2, H2+, H3+ and HeH+ sum to less than one, and so do the helium
      ! stages and the metal stages.
      x(:)  = 0.0d0
      x(1)  = 0.30d0
      x(2)  = 0.20d0
      x(3)  = 0.05d0
      x(4)  = 0.10d0
      x(5)  = 1.0d-4
      x(6)  = 1.0d-5
      x(7)  = 1.0d-6
      do i = 8,nmax
         x(i) = 0.30d0
         if (mod(i,2) .eq. 1) x(i) = 0.10d0
      enddo
      constraint = x(1) - x_fix_a

      ! --- the routine itself, with nothing imposed ---
      ieq_cell%x_hp_fixed = .false.
      ieq_cell%x_hp_fix   = x_fix_a
      do i = 1,nmax
         probe_in(i) = dble(i)*1.5d0 - 4.0d0
      enddo
      probe = probe_in
      call impose_transported_ionization_fractions(ieq_cell, x, probe)
      call check_absolute('routine_writes_nothing_when_unset',            &
           spread_of(probe, probe_in, 1, nmax), 0.0d0, 0.0d0)

      ! --- the routine itself, with the H+ fraction imposed ---
      ieq_cell%x_hp_fixed = .true.
      probe = probe_in
      call impose_transported_ionization_fractions(ieq_cell, x, probe)
      call check_absolute('routine_writes_the_constraint_in_row_one',     &
           probe(1) - constraint, 0.0d0, 0.0d0)
      call check_absolute('routine_leaves_the_helium_rows_alone',         &
           spread_of(probe, probe_in, 2, nmax), 0.0d0, 0.0d0)

      ! --- the same three statements inside each system ---
      do is = 1,nsys
         fvec_off = 0.0d0
         fvec_alt = 0.0d0
         fvec_on  = 0.0d0

         ieq_cell%x_hp_fixed = .false.
         ieq_cell%x_hp_fix   = x_fix_a
         call residual_of(is, nrow(is), x, fvec_off)

         ieq_cell%x_hp_fix   = x_fix_b
         call residual_of(is, nrow(is), x, fvec_alt)

         ieq_cell%x_hp_fixed = .true.
         ieq_cell%x_hp_fix   = x_fix_a
         call residual_of(is, nrow(is), x, fvec_on)

         row1_imposed(is) = fvec_on(1)

         ! (a) the imposed row is the constraint, exactly.
         call check_absolute(trim(sysname(is))//'_imposed_row_is_the_'    &
              //'constraint', fvec_on(1) - constraint, 0.0d0, 0.0d0)

         ! (a) and it is the ONLY row the substitution moved.
         if (nrow(is) .ge. 2) then
            call check_absolute(trim(sysname(is))//'_imposed_moves_no_'   &
                 //'other_row', spread_of(fvec_on, fvec_off, 2, nrow(is)),&
                 0.0d0, 0.0d0)
         endif

         ! (b) with nothing imposed the carried value is not read: the
         !     residual at two different carried values is the same
         !     floating-point number in every row.
         call check_absolute(trim(sysname(is))//'_unimposed_ignores_the_' &
              //'carried_value', spread_of(fvec_off, fvec_alt, 1,         &
              nrow(is)), 0.0d0, 0.0d0)

         ! (b) and row 1 is still the photoionization balance, not the
         !     constraint: the two differ by the whole balance.
         dmax = abs(fvec_off(1) - constraint)
         call check_positive(trim(sysname(is))//'_unimposed_row_is_a_'    &
              //'balance', dmax)
      enddo

      ! (c) one substitution, written once: the imposed row carries the
      !     same number in all seven systems.
      dmax = 0.0d0
      do is = 2,nsys
         dmax = max(dmax, abs(row1_imposed(is) - row1_imposed(1)))
      enddo
      call check_absolute('imposed_row_agrees_across_the_seven_systems',  &
           dmax, 0.0d0, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'imposed_ionization_rows: ',                 &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'imposed_ionization_rows: all assertions passed'

      contains

      subroutine set_cell_rates()
      ! One cell of a warm, partly ionized, partly molecular gas: every rate
      ! the seven residuals read, none of them zero, so that no row is
      ! trivially satisfied and the balance of row 1 is far from the
      ! constraint.
      ieq_cell%jcell       = 1
      ieq_cell%nh          = 1.0d8
      ieq_cell%nhe         = 8.0d6
      ieq_cell%T_K         = 8.0d3
      ieq_cell%ntot        = 1.1d8
      ieq_cell%P_HI        = 1.0d-4
      ieq_cell%P_HeI       = 5.0d-5
      ieq_cell%P_HeII      = 1.0d-6
      ieq_cell%rchiiB      = 2.6d-13
      ieq_cell%rcheiiB     = 4.5d-13
      ieq_cell%rcheiiiB    = 1.3d-12
      ieq_cell%a_ion_HI    = 1.0d-11
      ieq_cell%a_ion_HeI   = 3.0d-12
      ieq_cell%a_ion_HeII  = 1.0d-14
      ieq_cell%a_ion_HeITR = 1.0d-9
      ieq_cell%rcheiTR     = 1.0d-13
      ieq_cell%A31         = 1.272d-4
      ieq_cell%P_HeITR     = 1.0d-5
      ieq_cell%q13         = 1.0d-10
      ieq_cell%q31a        = 5.0d-10
      ieq_cell%q31b        = 1.0d-9
      ieq_cell%Q31         = 5.0d-10
      ieq_cell%kcx_He0_Hp  = 1.0d-14
      ieq_cell%kcx_Hep_H0  = 1.0d-14
      ieq_cell%P_H2        = 1.0d-8
      ieq_cell%P_H2_di     = 1.0d-9
      ieq_cell%P_H2_dd     = 0.0d0
      ieq_cell%P_H2_nd     = 0.0d0
      ieq_cell%k_LW        = 1.0d-10
      ieq_cell%n_ofam      = 0.0d0
      ieq_cell%n_co        = 0.0d0
      ieq_cell%x_h2_fixed  = .false.
      ieq_cell%x_ox_fixed  = .false.
      ieq_cell%x_h2_fix    = 0.0d0
      ieq_cell%x_oh_fix    = 0.0d0
      ieq_cell%x_h2o_fix   = 0.0d0
      end subroutine set_cell_rates

      subroutine residual_of(isys, n, xin, fout)
      ! The production residual of one system, on the cell standing in
      ! ieq_cell. params is the MINPACK transport argument and is unread by
      ! every one of them.
      integer, intent(in)  :: isys, n
      real*8,  intent(in)  :: xin(nmax)
      real*8,  intent(out) :: fout(nmax)
      real*8  :: xl(n), fl(n), par40(40), par60(60)
      integer :: iflag
      iflag    = 1
      par40(:) = 0.0d0
      par60(:) = 0.0d0
      xl(1:n)  = xin(1:n)
      fl(:)    = 0.0d0
      select case (isys)
      case (1)
         call ion_system_H(n, xl, fl, iflag, par40)
      case (2)
         call ion_system_HeH(n, xl, fl, iflag, par40)
      case (3)
         call ion_system_HeH_TR(n, xl, fl, iflag, par40)
      case (4)
         call ion_system_HeH_metals(n, xl, fl, iflag, par60)
      case (5)
         call ion_system_HeH_TR_metals(n, xl, fl, iflag, par60)
      case (6)
         call ion_system_HeH_mol(n, xl, fl, iflag, par40)
      case (7)
         call ion_system_HeH_mol_metals(n, xl, fl, iflag, par60)
      end select
      fout(:)   = 0.0d0
      fout(1:n) = fl(1:n)
      end subroutine residual_of

      double precision function spread_of(a, b, i1, i2)
      ! The largest absolute difference between two residuals over a range
      ! of rows. Zero means every row of the range is the same
      ! floating-point number in both.
      real*8, intent(in)  :: a(nmax), b(nmax)
      integer, intent(in) :: i1, i2
      integer :: k
      spread_of = 0.0d0
      do k = i1,i2
         spread_of = max(spread_of, abs(a(k) - b(k)))
      enddo
      end function spread_of

      end program imposed_ionization_rows

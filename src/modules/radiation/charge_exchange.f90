      module charge_exchange
      ! Charge-exchange (charge-transfer) reactions from Huang et al. (2023),
      ! ApJ 951, 123, Table 4. Each reaction transfers a single electron; the
      ! two directions of a pair are tabulated and fitted independently (the
      ! reverse rates were obtained by Huang via microscopic balance where no
      ! direct data exist). The image-verified source transcription lives in
      ! docs/charge_exchange_table4.md.
      !
      ! Active set:
      !   cx_full = .false. (default): the 23 metal + H/H+ reactions (Group A),
      !             which couple each metal's ionization balance to the H+
      !             fraction. This is the default.
      !   cx_full = .true. : the full Table 4 -- adds He + H (Group B),
      !             metal + He (Group C) and metal + metal (Group D).
      !
      ! NOTE on C/N/O: Table 4 contains only the neutral<->singly-ionized
      ! charge transfer for C, N, O (rows A13-A18). Re-sourcing C/N/O to Huang
      ! therefore DROPS the doubly-ionized terms (C2+/N2+/O2+ + H) that the
      ! previous Kingdon & Ferland implementation carried. This is an
      ! intentional behavior change to match the paper and is documented for
      ! validation.
      !
      ! Species are addressed by an element code: the canonical metal indices
      ! iel_C..iel_Fe (= 1..10 from species_table) plus cx_H = 11 and
      ! cx_He = 12. A reaction is a donor (loses an electron, stage s -> s+1)
      ! plus an acceptor (gains one, stage s -> s-1). The residual
      ! contribution is assembled generically in cx_add_to_fvec, so the
      ! zero-abundance guard is preserved automatically: the rate
      ! R = k * n(donor) * n(acceptor) vanishes when either reactant is absent.

      use species_table, only: iel_C, iel_O, iel_N, iel_Mg, iel_Si,        &
                               iel_Ca, iel_Na, iel_K, iel_S, iel_Fe,       &
                               n_melem

      implicit none
      private

      public :: cx_init, cx_set_cell, cx_add_to_fvec, cx_add_to_jac, cx_full
      public :: cx_metal_base

      ! Pseudo-element codes for H and He (metals use iel_* = 1..10).
      integer, parameter :: cx_H  = 11
      integer, parameter :: cx_He = 12

      ! Total number of Table-4 rows carried (23 A + 2 B + 6 C + 32 D).
      integer, parameter :: n_cxreac = 63

      ! Reaction descriptors (row r = 1..n_cxreac, in Table-4 order A,B,C,D).
      ! One electron passes from donor to acceptor:
      !   donor    : element cx_don_el(r) at stage cx_don_stg(r) -> stage+1
      !   acceptor : element cx_acc_el(r) at stage cx_acc_stg(r) -> stage-1
      integer, parameter :: cx_don_el(n_cxreac) = [                        &
      ! -- Group A: metal + H / H+ (A1..A23) --
        iel_Mg, cx_H,   iel_Mg, cx_H,   iel_Fe, cx_H,   iel_Fe, cx_H,      &
        iel_Si, cx_H,   iel_Si, cx_H,   iel_O,  cx_H,   iel_C,  cx_H,      &
        iel_N,  cx_H,   iel_S,  cx_H,   iel_Na, cx_H,   iel_K,             &
      ! -- Group B: He + H / H+ (B1..B2) --
        cx_He,  cx_H,                                                      &
      ! -- Group C: metal + He / He+ (C1..C6) --
        iel_Si, cx_He,  iel_C,  cx_He,  iel_O,  cx_He,                     &
      ! -- Group D: metal + metal (D1..D32) --
        iel_C,  iel_Si, iel_Mg, iel_Si, iel_Fe, iel_Si, iel_Na, iel_Mg,    &
        iel_Mg, iel_C,  iel_Mg, iel_N,  iel_Na, iel_Fe, iel_Fe, iel_C,     &
        iel_Fe, iel_O,  iel_Mg, iel_S,  iel_Na, iel_S,  iel_Fe, iel_S,     &
        iel_S,  iel_C,  iel_Fe, iel_N,  iel_Na, iel_Ca, iel_K,  iel_Na ]

      integer, parameter :: cx_don_stg(n_cxreac) = [                       &
      ! -- A --
        0, 0, 1, 0,  0, 0, 1, 0,  0, 0, 1, 0,  0, 0,  0, 0,  0, 0,         &
        0, 0,  0, 0,  0,                                                   &
      ! -- B --
        0, 0,                                                              &
      ! -- C --
        0, 0,  0, 0,  0, 0,                                                &
      ! -- D --
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,                 &
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0 ]

      integer, parameter :: cx_acc_el(n_cxreac) = [                        &
      ! -- Group A --
        cx_H,   iel_Mg, cx_H,   iel_Mg, cx_H,   iel_Fe, cx_H,   iel_Fe,    &
        cx_H,   iel_Si, cx_H,   iel_Si, cx_H,   iel_O,  cx_H,   iel_C,     &
        cx_H,   iel_N,  cx_H,   iel_S,  cx_H,   iel_Na, cx_H,              &
      ! -- Group B --
        cx_H,   cx_He,                                                    &
      ! -- Group C --
        cx_He,  iel_Si, cx_He,  iel_C,  cx_He,  iel_O,                     &
      ! -- Group D --
        iel_Si, iel_C,  iel_Si, iel_Mg, iel_Si, iel_Fe, iel_Mg, iel_Na,   &
        iel_C,  iel_Mg, iel_N,  iel_Mg, iel_Fe, iel_Na, iel_C,  iel_Fe,   &
        iel_O,  iel_Fe, iel_S,  iel_Mg, iel_S,  iel_Na, iel_S,  iel_Fe,   &
        iel_C,  iel_S,  iel_N,  iel_Fe, iel_Ca, iel_Na, iel_Na, iel_K ]

      integer, parameter :: cx_acc_stg(n_cxreac) = [                       &
      ! -- A --
        1, 1, 1, 2,  1, 1, 1, 2,  1, 1, 1, 2,  1, 1,  1, 1,  1, 1,         &
        1, 1,  1, 1,  1,                                                   &
      ! -- B --
        1, 1,                                                              &
      ! -- C --
        1, 1,  1, 1,  1, 1,                                                &
      ! -- D --
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,                 &
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1 ]

      ! .true. for the Group-A metal + H reactions (the default active set);
      ! .false. for Groups B/C/D (enabled only when cx_full = .true.).
      logical, parameter :: cx_default(n_cxreac) = [                       &
      ! -- A (23) --
        .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,      &
        .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,      &
        .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,      &
        .true.,  .true.,                                                   &
      ! -- B (2) --
        .false., .false.,                                                  &
      ! -- C (6) --
        .false., .false., .false., .false., .false., .false.,             &
      ! -- D (32) --
        .false., .false., .false., .false., .false., .false., .false.,    &
        .false., .false., .false., .false., .false., .false., .false.,    &
        .false., .false., .false., .false., .false., .false., .false.,    &
        .false., .false., .false., .false., .false., .false., .false.,    &
        .false., .false., .false., .false. ]

      ! Switch: .false. = Group A only (metal + H); .true. = full Table 4.
      logical, save :: cx_full = .false.

      ! fvec row index of the first metal element's neutral<->+ boundary.
      ! Default 4 matches ion_system_HeH_metals (x1=HII, x2=HeII, x3=HeIII,
      ! then metals at 4+2*(e-1)). The merged HeITR+metals system inserts the
      ! He-triplet unknown at row 4, pushing the metals to row 5; that solver
      ! sets cx_metal_base = 5 around its solve and resets to 4 afterward.
      integer, save :: cx_metal_base = 4

      ! Active reaction list (row indices) and per-cell rate coefficients.
      ! cx_act / cx_nact are built ONCE by cx_init (read-only during the sweep)
      ! and stay shared; cx_kc and cx_metal_base are PER-CELL state, so they are
      ! threadprivate now that the ionization cell sweep runs OpenMP-parallel.
      ! cx_kc is lazily allocated per thread in cx_set_cell; cx_metal_base is
      ! broadcast to each thread (copyin) at the parallel region and toggled
      ! 4<->5 per cell within a thread.
      integer, save :: cx_nact = 0
      integer, allocatable, save :: cx_act(:)
      real*8,  allocatable, save :: cx_kc(:)
      !$omp threadprivate(cx_kc, cx_metal_base)

      ! Upper clamp on any evaluated rate [cm^3 s^-1]; guards the lnT-
      ! polynomial fits (N/S/Na/K), which diverge as T -> 1 K (never reached
      ! in the model). The lower clamp is 0 -- a few KF96 bracket fits dip
      ! slightly negative below their validity floor (e.g. Mg+ + H+ at
      ! T < ~1500 K), where the true rate is negligible anyway.
      real*8, parameter :: cx_kc_max = 1.0d-6

      contains

      ! Build the active reaction index list for the current mode (cx_full,
      ! set beforehand from metals.inp). Called once by the driver after the
      ! input is read and before the ionization sweep.
      subroutine cx_init
      integer :: r, k
      cx_nact = 0
      do r = 1, n_cxreac
         if (cx_default(r) .or. cx_full) cx_nact = cx_nact + 1
      enddo
      if (allocated(cx_act)) deallocate(cx_act)
      if (allocated(cx_kc))  deallocate(cx_kc)
      allocate(cx_act(cx_nact), cx_kc(cx_nact))
      k = 0
      do r = 1, n_cxreac
         if (cx_default(r) .or. cx_full) then
            k = k + 1
            cx_act(k) = r
         endif
      enddo
      cx_kc = 0.0d0
      end subroutine cx_init

      ! Evaluate the active-reaction rate coefficients at temperature T [K].
      ! Called per cell before the hybrd1 solve, like set_metal_coeffs.
      subroutine cx_set_cell(T)
      real*8, intent(in) :: T
      integer :: i
      real*8  :: kc
      ! cx_kc is threadprivate: cx_init allocated only the master thread's copy,
      ! so each worker thread allocates its own on first use here (cx_nact is the
      ! shared, setup-once active-reaction count).
      if (.not. allocated(cx_kc)) allocate(cx_kc(cx_nact))
      do i = 1, cx_nact
         kc = cx_rate(cx_act(i), T)
         cx_kc(i) = min(max(kc, 0.0d0), cx_kc_max)
      enddo
      end subroutine cx_set_cell

      ! Add the charge-exchange source terms to an already-built residual
      ! vector. dens(el,stage) is assembled from the per-species densities;
      ! each reaction moves R from the donor's lower boundary (+R) to the
      ! acceptor's lower boundary (-R), reproducing the previous hard-coded
      ! C/N/O terms exactly while scaling to all metals.
      subroutine cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                 &
                                n_hi, n_hii, n_hei, n_heii, n_heiii)
      integer, intent(in)    :: N_eq
      real*8,  intent(inout) :: fvec(N_eq)
      real*8,  intent(in)    :: nm0(:), nm1(:), nm2(:)
      real*8,  intent(in)    :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8  :: dens(12,0:2), rrate
      integer :: i, r, de, ds, ae, as, idon, iacc

      dens = 0.0d0
      do i = 1, n_melem
         dens(i,0) = nm0(i)
         dens(i,1) = nm1(i)
         dens(i,2) = nm2(i)
      enddo
      dens(cx_H,0)  = n_hi
      dens(cx_H,1)  = n_hii
      dens(cx_He,0) = n_hei
      dens(cx_He,1) = n_heii
      dens(cx_He,2) = n_heiii

      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r);  ds = cx_don_stg(r)
         ae = cx_acc_el(r);  as = cx_acc_stg(r)
         rrate = cx_kc(i)*dens(de,ds)*dens(ae,as)
         idon  = cx_fvidx(de, ds)       ! donor's lower-stage boundary
         iacc  = cx_fvidx(ae, as-1)     ! acceptor's lower-stage boundary
         fvec(idon) = fvec(idon) + rrate
         fvec(iacc) = fvec(iacc) - rrate
      enddo
      end subroutine cx_add_to_fvec

      ! Map (element code, lower stage of a boundary) to the fvec row index:
      !   H        -> 1                              (HI <-> HII)
      !   He       -> 2 + lower                      (HeI<->HeII, HeII<->HeIII)
      !   metal e  -> cx_metal_base + 2*(e-1) + lower (X0<->X+,  X+<->X++)
      ! cx_metal_base = 4 for ion_system_HeH_metals, 5 for the merged
      ! HeITR+metals system (He-triplet unknown occupies row 4).
      integer function cx_fvidx(el, lower)
      integer, intent(in) :: el, lower
      if (el .eq. cx_H) then
         cx_fvidx = 1
      else if (el .eq. cx_He) then
         cx_fvidx = 2 + lower
      else
         cx_fvidx = cx_metal_base + 2*(el-1) + lower
      endif
      end function cx_fvidx

      ! Rate coefficient [cm^3 s^-1] for Table-4 row id at temperature T [K].
      ! Coefficients are transcribed verbatim from docs/charge_exchange_table4.md.
      real*8 function cx_rate(id, T)
      integer, intent(in) :: id
      real*8,  intent(in) :: T
      real*8 :: t4, tr, lnT, lnT2, lnT3, lnT4

      t4   = T/1.0d4
      tr   = T/300.0d0
      lnT  = log(T)
      lnT2 = lnT*lnT
      lnT3 = lnT2*lnT
      lnT4 = lnT3*lnT

      select case (id)
      ! ---------------- Group A: metal + H / H+ ----------------
      case (1);  cx_rate = kf(9.76d-12, 3.14d0,   55.54d0, -1.12d0,   0.0d0)   ! A1  Mg + H+
      case (2);  cx_rate = kf(2.95d-12, 3.28d0,   55.54d0, -1.12d0,   6.91d0)  ! A2  Mg+ + H
      case (3);  cx_rate = kf(7.6d-14,  0.0d0,    -1.97d0, -4.32d0,   1.67d0)  ! A3  Mg+ + H+
      case (4);  cx_rate = kf(8.58d-14, 2.49d-3,  0.0293d0,-4.33d0,   0.0d0)   ! A4  Mg2+ + H
      case (5);  cx_rate = 4.0d-9                                              ! A5  Fe + H+
      case (6);  cx_rate = gp(1.16d-9,  0.072d0,  6.61d0)                      ! A6  Fe+ + H
      case (7);  cx_rate = gp(2.3d-9,   0.0d0,    3.0d0)                       ! A7  Fe+ + H+
      case (8);  cx_rate = kf(1.26d-9,  0.0772d0, -0.41d0, -7.31d0,   0.0d0)   ! A8  Fe2+ + H
      case (9);  cx_rate = gp(7.41d-11, 0.85d0,   0.0d0)                       ! A9  Si + H+
      case (10); cx_rate = gp(4.71d-11, 0.95d0,   6.32d0)                      ! A10 Si+ + H
      case (11); cx_rate = kf(4.1d-10,  0.24d0,   3.17d0, -4.18d-3,   3.18d0)  ! A11 Si+ + H+
      case (12); cx_rate = kf(1.26d-9,  0.24d0,   3.17d0, -4.18d-3,   0.0d0)   ! A12 Si2+ + H
      case (13); cx_rate = 2.08d-9*t4**0.405d0 + 1.11d-11*t4**(-0.458d0)       ! A13 O + H+
      case (14); cx_rate = (1.26d-9*t4**0.517d0 + 4.25d-10*t4**6.69d-3)       &
                           *exp(-227.0d0/T)                                    ! A14 O+ + H
      case (15); cx_rate = gp(1.31d-15, 0.213d0,  0.0d0)                       ! A15 C + H+
      case (16); cx_rate = gp(6.3d-17,  1.96d0,   17.0d0)                      ! A16 C+ + H
      case (17); cx_rate = p4(-35.4d0,  1.94d0,  -0.154d0, -6.3d-3, -1.16d-3, 0.0d0)  ! A17 N + H+
      case (18); cx_rate = p4(-40.1d0,  6.4d0,   -1.75d0,   0.18d0, -5.96d-3, 0.0d0)  ! A18 N+ + H
      case (19); cx_rate = p4(-50.0d0,  13.3d0,  -2.77d0,   0.243d0,-7.24d-3, 0.0d0)  ! A19 S + H+
      case (20); cx_rate = p4(-50.14d0, 13.3d0,  -2.77d0,   0.243d0,-7.24d-3, 3.76d0) ! A20 S+ + H
      case (21); cx_rate = p4(48.2d0,  -43.0d0,   8.74d0,  -0.77d0, -0.0256d0, 0.0d0) ! A21 Na + H+
      case (22); cx_rate = p4(46.6d0,  -41.6d0,   8.02d0,  -0.65d0,  0.0197d-3, 9.817d0) ! A22 Na+ + H
      case (23); cx_rate = p4(-27.8d0,  0.125d0,  0.0663d0,-0.0237d0,-1.36d-3, 0.0d0)  ! A23 K + H+
      ! ---------------- Group B: He + H / H+ ----------------
      case (24); cx_rate = 1.75d-11*tr**(-0.75d0)*exp(-12.75d0/t4)             ! B1 He + H+
      case (25); cx_rate = 1.25d-15*tr**0.25d0                                 ! B2 He+ + H
      ! ---------------- Group C: metal + He / He+ ----------------
      case (26); cx_rate = fsi(T)                                             ! C1 Si + He+
      case (27); cx_rate = 1.415d0*fsi(T)*T**0.103d0*exp(-19.1d0/t4)          ! C2 Si+ + He
      case (28); cx_rate = gp(2.5d-15,  1.597d0, 0.0d0)                        ! C3 C + He+
      case (29); cx_rate = gp(6.75d-15, 1.654d0, 15.5d0)                       ! C4 C+ + He
      case (30); cx_rate = 4.99d-15*t4**0.379d0                              &
                           + 2.78d-15*t4**(-0.216d0)*exp(t4/81.97d0)           ! C5 O + He+
      case (31); cx_rate = 3.2d0*(5.0d-15*t4**0.38d0                         &
                           + 2.78d-15*t4**(-0.22d0)*exp(t4/81.97d0))         &
                           *T**0.0377d0*exp(-12.7d0/t4)                        ! C6 O+ + He
      ! ---------------- Group D: metal + metal ----------------
      case (32); cx_rate = (0.724d0*T**0.0463d0*exp(-3.61d0/t4))/dcsi(T)      ! D1  C + Si+
      case (33); cx_rate = 1.0d0/dcsi(T)                                      ! D2  C+ + Si
      case (34); cx_rate = 2.9d-9                                             ! D3  Mg + Si+
      case (35); cx_rate = gp(9.8d-10,  -0.0264d0, 0.59d0)                    ! D4  Mg+ + Si
      case (36); cx_rate = 1.9d-9                                             ! D5  Fe + Si+
      case (37); cx_rate = gp(1.4d-9,   -0.236d0,  0.299d0)                   ! D6  Fe+ + Si
      case (38); cx_rate = 1.0d-11                                            ! D7  Na + Mg+
      case (39); cx_rate = gp(3.76d-11, -0.0302d0, 2.909d0)                   ! D8  Na+ + Mg
      case (40); cx_rate = 1.1d-9                                             ! D9  Mg + C+
      case (41); cx_rate = gp(3.01d-10, -0.0832d0, 4.19d0)                    ! D10 Mg+ + C
      case (42); cx_rate = 1.2d-9                                             ! D11 Mg + N+
      case (43); cx_rate = gp(9.48d-10,  0.139d0,  7.99d0)                    ! D12 Mg+ + N
      case (44); cx_rate = 1.0d-11                                            ! D13 Na + Fe+
      case (45); cx_rate = gp(2.4d-11,  -0.0987d0, 3.21d0)                    ! D14 Na+ + Fe
      case (46); cx_rate = 2.6d-9                                             ! D15 Fe + C+
      case (47); cx_rate = gp(1.12d-9,   0.0147d0, 3.90d0)                    ! D16 Fe+ + C
      case (48); cx_rate = 1.71d-9                                            ! D17 Fe + O+
      case (49); cx_rate = gp(5.0d-10,  -0.0337d0, 6.63d0)                    ! D18 Fe+ + O
      case (50); cx_rate = 2.8d-10                                            ! D19 Mg + S+
      case (51); cx_rate = gp(5.35d-11,  0.121d0,  3.15d0)                    ! D20 Mg+ + S
      case (52); cx_rate = 2.6d-10                                            ! D21 Na + S+
      case (53); cx_rate = gp(1.86d-10,  0.151d0,  6.06d0)                    ! D22 Na+ + S
      case (54); cx_rate = 1.8d-10                                            ! D23 Fe + S+
      case (55); cx_rate = gp(5.38d-11,  0.052d0,  2.85d0)                    ! D24 Fe+ + S
      case (56); cx_rate = 3.26d-11*tr**0.289d0*exp(-338.0d0/T)               ! D25 S + C+
      case (57); cx_rate = gp(6.97d-11,  0.356d0,  1.06d0)                    ! D26 S+ + C
      case (58); cx_rate = 9.46d-10                                           ! D27 Fe + N+
      case (59); cx_rate = gp(1.17d-9,   0.07d0,   7.696d0)                   ! D28 Fe+ + N
      case (60); cx_rate = 3.0d-9                                             ! D29 Na + Ca+
      case (61); cx_rate = gp(1.7d-8,   -0.168d0,  1.13d0)                    ! D30 Na+ + Ca
      case (62); cx_rate = 1.0d-11                                            ! D31 K + Na+
      case (63); cx_rate = gp(6.58d-12, -0.203d0,  0.927d0)                   ! D32 K+ + Na
      case default; cx_rate = 0.0d0
      end select

      contains

         ! KF96 bracket form: a * T4^b * (1 + c*exp(d*T4)) * exp(-Ek/T4).
         ! d carries the table's multiplicative-T4 coefficient (negative);
         ! Ek is the trailing activation barrier in units of 1e4 K (Ek=0 if
         ! no trailing Boltzmann factor).
         real*8 function kf(a, b, c, d, Ek)
         real*8, intent(in) :: a, b, c, d, Ek
         kf = a * t4**b * (1.0d0 + c*exp(d*t4)) * exp(-Ek/t4)
         end function kf

         ! Power-law form: a * (T/300)^b * exp(-Ek/T4).
         real*8 function gp(a, b, Ek)
         real*8, intent(in) :: a, b, Ek
         gp = a * tr**b * exp(-Ek/t4)
         end function gp

         ! exp[ c0 + c1 lnT + c2 lnT^2 + c3 lnT^3 + c4 lnT^4 - Ek/T4 ],
         ! capped at cx_kc_max. These radiative-CT fits (Lin/Zhao/Dutta/
         ! Watanabe) are tiny over the model T-range and diverge only as
         ! T -> 1 K, which never occurs.
         real*8 function p4(c0, c1, c2, c3, c4, Ek)
         real*8, intent(in) :: c0, c1, c2, c3, c4, Ek
         real*8 :: arg
         arg = c0 + c1*lnT + c2*lnT2 + c3*lnT3 + c4*lnT4 - Ek/t4
         p4  = exp(min(arg, log(cx_kc_max)))
         end function p4

         ! f_Si(T) used by the Si + He charge-transfer rows (Satta 2013).
         real*8 function fsi(Tk)
         real*8, intent(in) :: Tk
         fsi = 3.32d-13*sqrt(Tk) + 1.2d-16*Tk + 4.2d-9/sqrt(Tk) - 7.9d-13
         end function fsi

         ! Denominator D_CSi(T) used by the C + Si charge-transfer rows
         ! (Satta 2013).
         real*8 function dcsi(Tk)
         real*8, intent(in) :: Tk
         dcsi = 1.87d8 + 5.09d10*Tk**(-0.527d0)
         end function dcsi

      end function cx_rate

      ! Analytic Jacobian counterpart of cx_add_to_fvec (Task 2): add
      ! d(CX source)/dx to the dense Jacobian fjac. Each per-species density is
      ! a LINEAR function of the unknowns, so for a reaction rate
      ! R = kc*D*A (donor density D, acceptor density A),
      !   dR/dx_k = kc*(dD/dx_k * A + D * dA/dx_k),
      ! added to row idon (+) and row iacc (-), exactly mirroring the
      ! fvec(idon)+=R, fvec(iacc)-=R bookkeeping. n_X (= met_ntot), n_h, n_he
      ! supply the linear-derivative coefficients.
      subroutine cx_add_to_jac(N_eq, fjac, nm0, nm1, nm2,                 &
                               n_hi, n_hii, n_hei, n_heii, n_heiii,       &
                               n_X, n_h, n_he)
      integer, intent(in)    :: N_eq
      real*8,  intent(inout) :: fjac(N_eq,N_eq)
      real*8,  intent(in)    :: nm0(:), nm1(:), nm2(:)
      real*8,  intent(in)    :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8,  intent(in)    :: n_X(:), n_h, n_he
      real*8  :: dens(12,0:2), D, A, kc, g
      integer :: i, r, de, ds, ae, as, idon, iacc, m
      integer :: nD, nA, kD(2), kA(2)
      real*8  :: cD(2), cA(2)

      dens = 0.0d0
      do i = 1, n_melem
         dens(i,0) = nm0(i); dens(i,1) = nm1(i); dens(i,2) = nm2(i)
      enddo
      dens(cx_H,0)  = n_hi
      dens(cx_H,1)  = n_hii
      dens(cx_He,0) = n_hei
      dens(cx_He,1) = n_heii
      dens(cx_He,2) = n_heiii

      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r); ds = cx_don_stg(r)
         ae = cx_acc_el(r); as = cx_acc_stg(r)
         D  = dens(de,ds);  A  = dens(ae,as);  kc = cx_kc(i)
         idon = cx_fvidx(de, ds)
         iacc = cx_fvidx(ae, as-1)
         ! donor-density derivatives: dR = kc * (dD) * A
         call cx_dens_lin(de, ds, n_X, n_h, n_he, nD, kD, cD)
         do m = 1, nD
            g = kc*cD(m)*A
            fjac(idon,kD(m)) = fjac(idon,kD(m)) + g
            fjac(iacc,kD(m)) = fjac(iacc,kD(m)) - g
         enddo
         ! acceptor-density derivatives: dR = kc * D * (dA)
         call cx_dens_lin(ae, as, n_X, n_h, n_he, nA, kA, cA)
         do m = 1, nA
            g = kc*D*cA(m)
            fjac(idon,kA(m)) = fjac(idon,kA(m)) + g
            fjac(iacc,kA(m)) = fjac(iacc,kA(m)) - g
         enddo
      enddo
      end subroutine cx_add_to_jac

      ! Linear decomposition of dens(el,stg) in the unknowns: returns the
      ! up-to-two unknown indices kidx(1:nidx) with d dens(el,stg)/d x = coef.
      ! Matches the x layout of ion_system_HeH_metals (x1=HII, x2=HeII,
      ! x3=HeIII; metal e at ix=4+2*(e-1): x_ix=X+, x_{ix+1}=X++).
      subroutine cx_dens_lin(el, stg, n_X, n_h, n_he, nidx, kidx, coef)
      integer, intent(in)  :: el, stg
      real*8,  intent(in)  :: n_X(:), n_h, n_he
      integer, intent(out) :: nidx, kidx(2)
      real*8,  intent(out) :: coef(2)
      integer :: ix
      nidx = 0; kidx = 0; coef = 0.0d0
      if (el .eq. cx_H) then
         nidx = 1; kidx(1) = 1
         if (stg .eq. 0) then
            coef(1) = -n_h
         else
            coef(1) =  n_h
         endif
      else if (el .eq. cx_He) then
         if (stg .eq. 0) then
            nidx = 2; kidx(1) = 2; kidx(2) = 3
            coef(1) = -n_he; coef(2) = -n_he
         else if (stg .eq. 1) then
            nidx = 1; kidx(1) = 2; coef(1) = n_he
         else
            nidx = 1; kidx(1) = 3; coef(1) = n_he
         endif
      else
         ix = cx_metal_base + 2*(el-1)
         if (stg .eq. 0) then
            nidx = 2; kidx(1) = ix; kidx(2) = ix+1
            coef(1) = -n_X(el); coef(2) = -n_X(el)
         else if (stg .eq. 1) then
            nidx = 1; kidx(1) = ix; coef(1) = n_X(el)
         else
            nidx = 1; kidx(1) = ix+1; coef(1) = n_X(el)
         endif
      endif
      end subroutine cx_dens_lin

      end module charge_exchange

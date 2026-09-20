      module charge_exchange
      ! Charge-exchange (charge-transfer) reactions from Huang et al. (2023),
      ! ApJ 951, 123, Table 4. Each reaction transfers a single electron; the
      ! two directions of a pair are tabulated and fitted independently (the
      ! reverse rates were obtained by Huang via microscopic balance where no
      ! direct data exist). The image-verified source transcription lives in
      ! docs/charge_exchange_table4.md.
      !
      ! Active set (generic cx_act assembly, cx_add_to_fvec/cx_add_to_jac):
      !   cx_full = .false. (default): the 23 metal + H/H+ reactions (Group A),
      !             which couple each metal's ionization balance to the H+
      !             fraction. This is the default.
      !   cx_full = .true. : adds metal + He (Group C) and metal + metal
      !             (Group D). It NO LONGER adds He + H (Group B): the He<->H
      !             pair is now applied by the dedicated he_h_cx_* routines
      !             below and is excluded from cx_act by cx_init (rows 24,25),
      !             so it is never counted twice.
      !
      ! He <-> H charge exchange (Group B, He0+H+ and He++H0) is handled
      ! separately so it is available in EVERY ionization system containing He
      ! -- not only the metals systems that call cx_add_to_fvec. It is gated by
      ! he_h_charge_exchange (default .true.), independent of cx_full, and the
      ! same two rate coefficients feed every system via he_h_cx_rates.
      !
      ! NOTE on C/N/O: Table 4 contains only the neutral<->singly-ionized
      ! charge transfer for C, N, O (rows A13-A18). Re-sourcing C/N/O to Huang
      ! therefore DROPS the doubly-ionized terms (C2+/N2+/O2+ + H) that the
      ! previous Kingdon & Ferland implementation carried. This is an
      ! intentional behavior change to match the paper and is documented for
      ! validation.
      !
      ! Group E is not from Table 4. It carries the one reaction whose absence
      ! from that table leaves the O III profile wrong by decades in the region
      ! a transit probes, O2+ + H0 -> O+ + H+. It is ON by default at the
      ! published rate (cx_o2p_h_scale = 1); "cx_O2p_H <scale>" in metals.inp
      ! rescales it, and 0 restores the Table-4-only reaction set exactly.
      ! cx_full does NOT control it -- cx_full means "all of Table 4", and this
      ! row is not in Table 4. See electron_capture_O2p_from_H below.
      !
      ! Species are addressed by an element code: the canonical metal indices
      ! iel_C..iel_Fe (= 1..10 from species_table) plus cx_H = 11 and
      ! cx_He = 12. A reaction is a donor (loses an electron, stage s -> s+1)
      ! plus an acceptor (gains one, stage s -> s-1). The residual
      ! contribution is assembled generically in cx_add_to_fvec, so the
      ! zero-abundance guard is preserved automatically: the rate
      ! R = k * n(donor) * n(acceptor) vanishes when either reactant is absent.
      !
      ! THE EQUATION BASIS OF THE HELIUM ROW IS THE CALLER'S, NOT THIS
      ! MODULE'S. cx_add_to_fvec and cx_add_to_jac write the ionization-
      ! positive form: a reaction is +R on the donor's lower-stage boundary
      ! row and -R on the acceptor's. The one row whose orientation differs
      ! between the solvers is the He I <-> He II boundary (row 2): it is
      ! written He I-gain positive in the summed He I balance of the
      ! He 2^3S systems (ion_residual_core, heh_tr_rows) and ionization
      ! positive everywhere else. The caller therefore passes he_row_sign,
      ! +1 for an ionization-positive He I <-> He II row and -1 for a
      ! He I-gain one, and it multiplies the increments of THAT row alone;
      ! the H row, the He II <-> He III row and the metal rows are
      ! ionization positive in every system and never carry it. The only
      ! reactions that reach the helium row are the metal + He / He+ pairs
      ! of Huang et al. (2023) Table 4 (group C, Si/C/O with He and He+),
      ! active only when cx_full = .true.; none of them touches He2+.
      !
      ! THE HELIUM REACTANT OF GROUP C IS THE GROUND SINGLET He(1^1S), NOT
      ! THE SUM OVER He I, so n_hei of cx_add_to_fvec and cx_add_to_jac is
      ! the singlet density. Table 4 writes those reactants as the bare
      ! element, and their barriers say which state it is: the reverse rates
      ! Si+ + He, C+ + He and O+ + He carry exp(-19.1/T4), exp(-15.5/T4) and
      ! exp(-12.7/T4), against ionization-potential differences of
      ! ground-state helium of 24.587 - 8.152 = 16.436 eV = 19.07e4 K,
      ! 24.587 - 11.260 = 13.327 eV = 15.47e4 K and
      ! 24.587 - 13.618 = 10.969 eV = 12.73e4 K. He(2^3S) lies 19.82 eV
      ! above the singlet and its ionization potential is 4.77 eV, so every
      ! one of those collisions is exothermic for it and runs at a rate this
      ! code carries nowhere; charging the metastable to a ground-state rate
      ! would be a reaction that is not in the set. The three metal systems
      ! and the constrained equilibrium therefore all pass the singlet.
      !
      ! TWO EQUATION BASES, ONE RATE. A reaction is one event with one
      ! volumetric rate R = k n(donor) n(acceptor) [cm^-3 s^-1], and the
      ! solvers differ only in what they write it into. The local
      ! ionization systems carry BOUNDARY FLOWS as their rows (the net
      ! upward flow across the boundary between two stages of an element),
      ! and cx_add_to_fvec adds R to them. The transported ionization
      ! stages carry the STAGE DENSITIES themselves, and
      ! charge_exchange_stage_sources returns the stoichiometric source of
      ! each stage, -R on the donor's own stage and +R on the stage above
      ! it, -R on the acceptor's and +R on the stage below. Both read the
      ! rates from cx_reaction_rates, so the rate and the reaction table
      ! are written once and a source assembled in either basis describes
      ! the same events. The two are equivalent term by term: for an
      ! element whose boundary flows are u_0 and u_1, the stage sources are
      ! S_0 = -u_0, S_1 = u_0 - u_1 and S_2 = u_1, which is the conversion
      ! src/tests/charge_exchange_rows asserts reaction by reaction.
      !
      ! WHO OWNS THE RATE COEFFICIENTS. cx_kc and cx_metal_base are
      ! thread-local and hold ONE cell: cx_set_cell(T) fills cx_kc for that
      ! cell's temperature and records it in cx_cell_T. Every caller of
      ! cx_add_to_fvec, cx_add_to_jac or cx_add_to_turnover must have called
      ! cx_set_cell for the cell it is assembling, on the same thread. The
      ! ionization cell sweep does (ionization_equilibrium), and so must any
      ! other evaluator of these rows, for instance a transported stage
      ! source. cx_add_to_fvec, cx_add_to_jac and cx_add_to_turnover take
      ! the cell temperature and stop the run if it is not the one cx_kc was
      ! filled at, so a stale set of rates inherited from whatever cell the
      ! thread handled last cannot be used silently. The turnover's cell is
      ! optional in the interface while two of its call sites still do not
      ! state one; a call that states none is refused when the thread holds
      ! no cell at all, since then there is no rate set to bound a row
      ! with.

      use species_table, only: iel_C, iel_O, iel_N, iel_Mg, iel_Si,        &
                               iel_Ca, iel_Na, iel_K, iel_S, iel_Fe,       &
                               n_melem

      implicit none
      private

      public :: cx_init, cx_set_cell, cx_add_to_fvec, cx_add_to_jac, cx_full
      public :: cx_metal_base, cx_add_to_turnover
      ! The reaction set as stoichiometric number-density sources, for a
      ! solver whose unknowns ARE stage densities rather than boundary
      ! flows (the transported ionization stages).
      public :: charge_exchange_stage_sources
      ! Dedicated He <-> H charge-exchange pair (Group B), available in every
      ! system with He, independent of cx_full.
      public :: he_h_charge_exchange, he_h_cx_rates
      public :: he_h_cx_fvec, he_h_cx_jac, he_h_cx_fvec_adv
      ! Group E: the omitted O2+ + H electron capture, as a sensitivity dial.
      public :: cx_o2p_h_scale, electron_capture_O2p_from_H

      ! Pseudo-element codes for H and He (metals use iel_* = 1..10).
      integer, parameter :: cx_H  = 11
      integer, parameter :: cx_He = 12

      ! Table-4 row ids of the He <-> H pair (Group B). Excluded from the
      ! generic cx_act set (cx_init) because the pair is applied by he_h_cx_*.
      integer, parameter :: cx_B1_He0_Hp = 24   ! He0 + H+ -> He+ + H0
      integer, parameter :: cx_B2_Hep_H0 = 25   ! He+ + H0 -> He0 + H+

      ! Master switch for the He <-> H charge-exchange pair (default on). Set
      ! from input.inp ("He_H_charge_exchange: False") by input_read.
      logical, save :: he_h_charge_exchange = .true.

      ! Reaction rows carried: 63 from Table 4 (23 A + 2 B + 6 C + 32 D) plus
      ! the one Group-E row that is NOT in Table 4 (E1, O2+ + H0).
      integer, parameter :: n_cxreac = 64
      ! Row id of the Group-E reaction, and its scale factor.
      !   Source  : Barragan, Errea, Mendez, Rabadan & Riera (2006), ApJ 636,
      !             544, Table 1 column k1 -- molecular close-coupling on MELD
      !             multireference-CI wavefunctions; the calculation CHIANTI
      !             v11 (Dufresne et al. 2024) adopts and Cloudy fits.
      !   Validity: 1e2 - 1e5 K, the calculated range. Outside it the rate is
      !             held at the endpoint, not extrapolated (see
      !             electron_capture_O2p_from_H).
      !   Default : 1.0 = the published rate, made the default on 2026-08-30
      !             (docs/Update_EXHALE_stage1.md section 107, measurement in section
      !             103). Charge transfer beats O III recombination wherever
      !             n(H0)/n_e exceeds a few 1e-3 -- k_CT/alpha_rec is
      !             122/280/617 at 5e3/1e4/2e4 K -- which is the whole launch
      !             region, so omitting the row left x(O III) too high by 4.7
      !             decades at 1.05 Rp on HD 209458 b.
      ! Setting "cx_O2p_H 0" in metals.inp reproduces the Table-4-only reaction
      ! set exactly: cx_init then leaves the row out of cx_act altogether, so
      ! no rate is evaluated and no residual term is assembled.
      integer, parameter :: cx_E1_O2p_H0 = 64
      real*8, save :: cx_o2p_h_scale = 1.0d0

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
        iel_S,  iel_C,  iel_Fe, iel_N,  iel_Na, iel_Ca, iel_K,  iel_Na,    &
      ! -- Group E: NOT in Table 4 (E1, O2+ + H0 -> O+ + H+) --
        cx_H ]

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
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,                 &
      ! -- E --
        0 ]

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
        iel_C,  iel_S,  iel_N,  iel_Fe, iel_Ca, iel_Na, iel_Na, iel_K,    &
      ! -- Group E --
        iel_O ]

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
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,                 &
      ! -- E: the acceptor is O at stage 2, i.e. O III -> O II --
        2 ]

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
        .false., .false., .false., .false.,                                &
      ! -- E (1): not in Table 4; gated on cx_o2p_h_scale, not on cx_full --
        .false. ]

      ! Switch: .false. = Group A only (metal + H); .true. = full Table 4.
      logical, save :: cx_full = .false.

      ! fvec row index of the first metal element's neutral<->+ boundary.
      ! Default 4 matches ion_system_HeH_metals (x1=HII, x2=HeII, x3=HeIII,
      ! then metals at 4+2*(e-1)). The merged HeITR+metals system inserts the
      ! He-triplet unknown at row 4, pushing the metals to row 5; that solver
      ! sets cx_metal_base = 5 around its solve and resets to 4 afterward.
      integer, save :: cx_metal_base = 4

      ! Active reaction list (row indices) and cell-by-cell rate coefficients.
      ! cx_act / cx_nact are built ONCE by cx_init (read-only during the sweep)
      ! and stay shared; cx_kc and cx_metal_base are CELL-BY-CELL state, so they are
      ! threadprivate now that the ionization cell sweep runs OpenMP-parallel.
      ! cx_kc is lazily allocated per thread in cx_set_cell; cx_metal_base is
      ! broadcast to each thread (copyin) at the parallel region and toggled
      ! 4<->5 per cell within a thread.
      integer, save :: cx_nact = 0
      integer, allocatable, save :: cx_act(:)
      real*8,  allocatable, save :: cx_kc(:)
      ! Temperature [K] the thread's cx_kc was filled at by cx_set_cell. The
      ! negative initial value is "this thread holds no cell", which no
      ! temperature can match, so an assembly that skipped cx_set_cell is
      ! refused rather than run on another cell's rates. It is threadprivate
      ! with cx_kc because it labels cx_kc, and it is deliberately NOT in the
      ! copyin list of the cell sweep: a worker entering the region holds no
      ! cell until it loads one.
      real*8, save :: cx_cell_T = -1.0d0
      !$omp threadprivate(cx_kc, cx_metal_base, cx_cell_T)

      ! Upper clamp on any evaluated rate [cm^3 s^-1]; guards the lnT-
      ! polynomial fits (N/S/Na/K), which diverge as T -> 1 K (never reached
      ! in the model). The lower clamp is 0 -- a few KF96 bracket fits dip
      ! slightly negative below their validity floor (e.g. Mg+ + H+ at
      ! T < ~1500 K), where the true rate is negligible anyway.
      real*8, parameter :: cx_kc_max = 1.0d-6

      contains

      ! Electron capture by doubly ionized oxygen from atomic hydrogen,
      !     O2+(2s^2 2p^2 3P) + H(1s) -> O+ + H+ ,
      ! rate coefficient in cm^3 s^-1. This is reaction (1) of Barragan,
      ! Errea, Mendez, Rabadan & Riera (2006), ApJ 636, 544, and the values
      ! below are their Table 1 column k1 as published (units 1e-9 cm^3 s^-1),
      ! read from the journal PDF. Their method is a molecular close-coupling
      ! expansion on MELD multireference-CI wavefunctions, quantal below
      ! 250 eV/amu, over 1e2 - 1e5 K; it is the calculation CHIANTI v11
      ! (Dufresne et al. 2024, ApJ 974, 71, appendix A.3.3) and Cloudy both
      ! adopt for O III -> O II charge transfer, superseding the
      ! Landau-Zener-era value of Kingdon & Ferland (1996, ApJS 106, 205),
      ! which is 1.6x higher at 5000 K and 1.2x at 1e4 K.
      !
      ! Interpolated linearly in log10(k) against log10(T) between the
      ! published points -- the table itself, not a fit to it. T is clamped to
      ! [1e2, 1e5] K, the range the calculation covers, so outside it the rate
      ! is held at the endpoint rather than extrapolated (the same convention
      ! as fine_structure_upsilon in Cool_coeff.f90). The authors' own caveat
      ! is that radiative electron capture, which they do not include,
      ! "probably limits the application of our rate coefficients at low
      ! temperatures".
      !
      ! For scale: the classical Langevin capture bound for this pair,
      ! k_L = 2 pi q sqrt(alpha/mu) with alpha(H) = 4.5 a0^3, is
      ! 3.9e-9 cm^3 s^-1, so the published rate is 16-32% of the rate at which
      ! the two particles reach small separation at all -- physically sensible.
      ! The reverse direction O+ + H+ -> O2+ + H0 is NOT carried: it is
      ! endothermic by IP(O II) - IP(H I) = 35.12112 - 13.5984 = 21.52 eV
      ! (NIST ASD), so detailed balance suppresses it by exp(-21.52 eV/kT),
      ! i.e. by 1e-11 at 1e4 K and 4e-6 at 2e4 K, relative to the forward rate.
      pure double precision function electron_capture_O2p_from_H(T)
      real*8, intent(in) :: T
      integer, parameter :: n_o2p = 26
      real*8, parameter :: T_o2p(n_o2p) = [                                &
         1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 1.0d3, 1.5d3, 2.0d3,          &
         2.5d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 7.0d3, 8.0d3, 9.0d3,          &
         1.0d4, 2.0d4, 3.0d4, 4.0d4, 5.0d4, 6.0d4, 7.0d4, 8.0d4,          &
         9.0d4, 1.0d5 ]
      real*8, parameter :: k_o2p(n_o2p) = [                                &
         0.53362d-9, 0.52153d-9, 0.49852d-9, 0.48016d-9, 0.46613d-9,      &
         0.43607d-9, 0.43975d-9, 0.45679d-9, 0.47931d-9, 0.50418d-9,      &
         0.55659d-9, 0.60995d-9, 0.66288d-9, 0.71464d-9, 0.76481d-9,      &
         0.81319d-9, 0.85973d-9, 1.24359d-9, 1.53854d-9, 1.78217d-9,      &
         1.98734d-9, 2.16168d-9, 2.31143d-9, 2.44169d-9, 2.55653d-9,      &
         2.65907d-9 ]
      real*8  :: x, w
      integer :: i
      x = log10(min(max(T, T_o2p(1)), T_o2p(n_o2p)))
      i = 1
      do while (i .lt. n_o2p - 1 .and. x .gt. log10(T_o2p(i+1)))
         i = i + 1
      enddo
      w = (x - log10(T_o2p(i)))/(log10(T_o2p(i+1)) - log10(T_o2p(i)))
      electron_capture_O2p_from_H =                                       &
         10.0d0**(log10(k_o2p(i)) + w*(log10(k_o2p(i+1))                  &
                                       - log10(k_o2p(i))))
      end function electron_capture_O2p_from_H

      ! Build the active reaction index list for the current mode (cx_full,
      ! set beforehand from metals.inp). Called once by the driver after the
      ! input is read and before the ionization sweep.
      subroutine cx_init
      integer :: r, k
      cx_nact = 0
      do r = 1, n_cxreac
         if (cx_row_is_active(r)) cx_nact = cx_nact + 1
      enddo
      if (allocated(cx_act)) deallocate(cx_act)
      if (allocated(cx_kc))  deallocate(cx_kc)
      allocate(cx_act(cx_nact), cx_kc(cx_nact))
      k = 0
      do r = 1, n_cxreac
         if (cx_row_is_active(r)) then
            k = k + 1
            cx_act(k) = r
         endif
      enddo
      cx_kc = 0.0d0
      end subroutine cx_init

      ! Membership test for the generic active set, written once so cx_init's
      ! two passes cannot drift apart.
      !   Group B (He<->H) is applied by he_h_cx_*, so it is never here.
      !   Group E is not in Table 4, so cx_full does not reach it; it is
      !   gated on its own scale factor, positive by default and zeroed by
      !   "cx_O2p_H 0".
      logical function cx_row_is_active(r)
      integer, intent(in) :: r
      if (r .eq. cx_B1_He0_Hp .or. r .eq. cx_B2_Hep_H0) then
         cx_row_is_active = .false.
      else if (r .eq. cx_E1_O2p_H0) then
         cx_row_is_active = (cx_o2p_h_scale .gt. 0.0d0)
      else
         cx_row_is_active = (cx_default(r) .or. cx_full)
      endif
      end function cx_row_is_active

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
      ! The rates now belong to this temperature; the assembly routines
      ! check that the cell they are handed is this one.
      cx_cell_T = T
      end subroutine cx_set_cell

      ! Refuse an assembly whose rate coefficients belong to another cell.
      ! cx_kc holds one cell at a time on each thread, so a caller that did
      ! not run cx_set_cell for the cell it is assembling would use whichever
      ! temperature that thread happened to load last, and the source it
      ! built would depend on the order the cells were handed out.
      subroutine cx_require_cell(T_cell, caller)
      real*8,           intent(in) :: T_cell
      character(len=*), intent(in) :: caller
      if (cx_cell_T .ne. T_cell) then
         write(*,*) '(charge_exchange) ERROR: ', trim(caller),            &
                    ' was called for a cell at T =', T_cell, ' K,'
         write(*,*) '  but the charge-exchange rate coefficients of this', &
                    ' thread were last filled at T =', cx_cell_T, ' K.'
         write(*,*) '  cx_kc holds one cell; call cx_set_cell(T) for the', &
                    ' cell being assembled, on this thread, first.'
         write(*,*) '  Aborting.'
         error stop 1
      endif
      end subroutine cx_require_cell

      ! Add the charge-exchange source terms to an already-built residual
      ! vector. dens(el,stage) is assembled from each species's densities;
      ! each reaction moves R from the donor's lower boundary (+R) to the
      ! acceptor's lower boundary (-R).
      !
      ! he_row_sign is the orientation of the caller's He I <-> He II row
      ! (row 2): +1 when that row is written He I -> He II positive, -1 when
      ! it is written He I-gain positive, as the summed He I balance of the
      ! He 2^3S systems is. It multiplies the increments of that row alone,
      ! so a reaction that destroys He II and makes He I enters a He I-gain
      ! row with the positive sign its physics requires. Every other row of
      ! every caller is ionization positive and keeps the sign written here.
      ! T_cell is the temperature of the cell being assembled, checked
      ! against the one cx_set_cell filled cx_kc at.
      subroutine cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                 &
                                n_hi, n_hii, n_hei, n_heii, n_heiii,       &
                                he_row_sign, T_cell)
      integer, intent(in)    :: N_eq
      real*8,  intent(inout) :: fvec(N_eq)
      real*8,  intent(in)    :: nm0(:), nm1(:), nm2(:)
      real*8,  intent(in)    :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8,  intent(in)    :: he_row_sign, T_cell
      real*8  :: dens(12,0:2), rrate(cx_nact), sdon, sacc
      integer :: i, r, de, ds, ae, as, idon, iacc

      call cx_reactant_densities(nm0, nm1, nm2, n_hi, n_hii, n_hei,       &
                                 n_heii, n_heiii, dens)
      call cx_reaction_rates(dens, T_cell, 'cx_add_to_fvec', rrate)

      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r);  ds = cx_don_stg(r)
         ae = cx_acc_el(r);  as = cx_acc_stg(r)
         idon  = cx_fvidx(de, ds)       ! donor's lower-stage boundary
         iacc  = cx_fvidx(ae, as-1)     ! acceptor's lower-stage boundary
         ! The He I <-> He II boundary is the one row whose orientation is
         ! the caller's; everything else is ionization positive.
         sdon = cx_he_i_boundary_sign(de, ds,   he_row_sign)
         sacc = cx_he_i_boundary_sign(ae, as-1, he_row_sign)
         fvec(idon) = fvec(idon) + sdon*rrate(i)
         fvec(iacc) = fvec(iacc) - sacc*rrate(i)
      enddo
      end subroutine cx_add_to_fvec

      ! The reactant densities of one cell, by element code and ionization
      ! stage: canonical metals 1..n_melem, then cx_H = 11 and cx_He = 12.
      ! The helium entry at stage 0 is the GROUND SINGLET He(1^1S), which
      ! is the reactant of the group C reactions (module header).
      subroutine cx_reactant_densities(nm0, nm1, nm2, n_hi, n_hii,        &
                                       n_hei, n_heii, n_heiii, dens)
      real*8, intent(in)  :: nm0(:), nm1(:), nm2(:)
      real*8, intent(in)  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8, intent(out) :: dens(12,0:2)
      integer :: i
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
      end subroutine cx_reactant_densities

      ! The volumetric rate of every active reaction, R = k n(donor)
      ! n(acceptor) [cm^-3 s^-1], at the cell whose rate coefficients this
      ! thread holds. THE ONE PLACE A CHARGE-EXCHANGE RATE IS FORMED: both
      ! equation bases read this array (module header), so a reaction
      ! cannot acquire a second spelling in a second solver. The rate
      ! vanishes when either reactant is absent, which is the zero-
      ! abundance guard of the whole set.
      subroutine cx_reaction_rates(dens, T_cell, caller, rrate)
      real*8,           intent(in)  :: dens(12,0:2)
      real*8,           intent(in)  :: T_cell
      character(len=*), intent(in)  :: caller
      real*8,           intent(out) :: rrate(cx_nact)
      integer :: i, r
      call cx_require_cell(T_cell, caller)
      do i = 1, cx_nact
         r = cx_act(i)
         rrate(i) = cx_kc(i)*dens(cx_don_el(r), cx_don_stg(r))            &
                            *dens(cx_acc_el(r), cx_acc_stg(r))
      enddo
      end subroutine cx_reaction_rates

      ! THE REACTION SET AS STOICHIOMETRIC NUMBER-DENSITY SOURCES
      ! [cm^-3 s^-1], stage by stage: one electron leaves the donor, so the
      ! donor's stage loses R and the stage above it gains R, and one
      ! electron arrives at the acceptor, so the acceptor's stage loses R
      ! and the stage below it gains R. This is the form a solver whose
      ! unknowns are the stage densities themselves needs; the local
      ! ionization systems, whose rows are boundary flows, take the same
      ! rates through cx_add_to_fvec (module header states the conversion).
      !
      ! s_H and s_He are the sources of the hydrogen and helium stages, and
      ! s_metal those of the metals, indexed by the canonical element order.
      ! The optional p_* and l_* are the GROSS production and loss of each
      ! stage, both positive, whose difference is the source: a stage whose
      ! charge exchange is fast carries production and loss orders above
      ! their difference, and a reader of the net cannot recover them.
      !
      ! Group B, the He <-> H pair, is NOT in this set (cx_init leaves it
      ! out): it is applied by he_h_cx_fvec in every system with helium,
      ! so counting it here would count it twice. Groups C and D enter only
      ! with cx_full, group E under cx_o2p_h_scale, exactly as they do in
      ! the residual assembly, because both read one active set.
      !
      ! T_cell is the temperature of the cell being evaluated, checked
      ! against the one cx_set_cell filled the rate coefficients at.
      subroutine charge_exchange_stage_sources(nm0, nm1, nm2,             &
                                n_hi, n_hii, n_hei, n_heii, n_heiii,      &
                                T_cell, s_H, s_He, s_metal,               &
                                p_H, l_H, p_He, l_He)
      real*8, intent(in)  :: nm0(:), nm1(:), nm2(:)
      real*8, intent(in)  :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8, intent(in)  :: T_cell
      real*8, intent(out) :: s_H(0:1), s_He(0:2)
      real*8, optional, intent(out) :: s_metal(n_melem,0:2)
      real*8, optional, intent(out) :: p_H(0:1), l_H(0:1)
      real*8, optional, intent(out) :: p_He(0:2), l_He(0:2)
      ! The working arrays carry one stage above the highest the reaction
      ! table reaches. Every donor of Table 4 and of group E is a NEUTRAL
      ! and every acceptor is singly or doubly ionized, so a reaction
      ! writes stages 0 to 2 only; the extra row exists so that a future
      ! table row promoting an already ionized donor is caught here by
      ! name instead of writing past the metal arrays.
      real*8  :: dens(12,0:2), rrate(cx_nact)
      real*8  :: ssrc(12,0:3), gprod(12,0:3), gloss(12,0:3)
      integer :: i, r, de, ds, ae, as, e

      call cx_reactant_densities(nm0, nm1, nm2, n_hi, n_hii, n_hei,       &
                                 n_heii, n_heiii, dens)
      call cx_reaction_rates(dens, T_cell,                                &
                             'charge_exchange_stage_sources', rrate)

      ssrc  = 0.0d0
      gprod = 0.0d0
      gloss = 0.0d0
      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r);  ds = cx_don_stg(r)
         ae = cx_acc_el(r);  as = cx_acc_stg(r)
         ssrc(de,ds)    = ssrc(de,ds)    - rrate(i)
         ssrc(de,ds+1)  = ssrc(de,ds+1)  + rrate(i)
         ssrc(ae,as)    = ssrc(ae,as)    - rrate(i)
         ssrc(ae,as-1)  = ssrc(ae,as-1)  + rrate(i)
         gloss(de,ds)   = gloss(de,ds)   + rrate(i)
         gprod(de,ds+1) = gprod(de,ds+1) + rrate(i)
         gloss(ae,as)   = gloss(ae,as)   + rrate(i)
         gprod(ae,as-1) = gprod(ae,as-1) + rrate(i)
      enddo
      if (maxval(abs(ssrc(:,3))) .gt. 0.0d0) then
         write(*,*) '(charge_exchange) ERROR: a reaction of the active',  &
                    ' set promotes a donor above ionization stage 2,'
         write(*,*) '  which the three-stage element layout of this',     &
                    ' module and of the metal systems cannot hold.'
         write(*,*) '  Aborting.'
         error stop 1
      endif

      s_H(0:1)  = ssrc(cx_H,0:1)
      s_He(0:2) = ssrc(cx_He,0:2)
      if (present(s_metal)) then
         do e = 1, n_melem
            s_metal(e,0:2) = ssrc(e,0:2)
         enddo
      endif
      if (present(p_H))  p_H(0:1)  = gprod(cx_H,0:1)
      if (present(l_H))  l_H(0:1)  = gloss(cx_H,0:1)
      if (present(p_He)) p_He(0:2) = gprod(cx_He,0:2)
      if (present(l_He)) l_He(0:2) = gloss(cx_He,0:2)
      end subroutine charge_exchange_stage_sources

      ! Orientation factor of the row that (element, lower stage) addresses:
      ! he_row_sign for the He I <-> He II boundary, whose equation basis
      ! differs between the solvers, and 1 for every other row.
      real*8 function cx_he_i_boundary_sign(el, lower, he_row_sign)
      integer, intent(in) :: el, lower
      real*8,  intent(in) :: he_row_sign
      if (el .eq. cx_He .and. lower .eq. 0) then
         cx_he_i_boundary_sign = he_row_sign
      else
         cx_he_i_boundary_sign = 1.0d0
      endif
      end function cx_he_i_boundary_sign

      ! Upper bound of the charge-exchange rate each residual row can carry,
      ! added onto the turnover scale of each row (the acceptance normalization of
      ! ioniz_eq). The convention is that of set_mol_turnover_rates
      ! (System_HeH_mol): every reactant is set to the WHOLE of its element,
      ! so each active reaction contributes kc * N(donor element) *
      ! N(acceptor element) to both rows it moves. el_tot is indexed by the
      ! same element codes as dens in cx_add_to_fvec: canonical metals
      ! 1..n_melem, then cx_H = 11 and cx_He = 12; an absent element carries
      ! el_tot = 0 and contributes nothing, exactly as its rate does. The
      ! rows are addressed through cx_fvidx, so cx_metal_base must hold the
      ! layout of the system being judged, as it must for cx_add_to_fvec.
      ! T_cell is the temperature of the cell whose turnover is being
      ! bounded, checked against the one cx_set_cell filled the rate
      ! coefficients at: the bound is kc(T) times two element counts, so
      ! rates belonging to another cell would normalize this cell's residual
      ! by another cell's turnover, and the acceptance verdict a caller then
      ! reads would depend on the order the cells were handed out.
      subroutine cx_add_to_turnover(srow, el_tot, T_cell)
      real*8, intent(inout) :: srow(*)
      real*8, intent(in)    :: el_tot(12)
      real*8, optional, intent(in) :: T_cell
      integer :: i, r, idon, iacc
      real*8  :: bound

      if (present(T_cell)) then
         call cx_require_cell(T_cell, 'cx_add_to_turnover')
      else if (cx_cell_T .lt. 0.0d0) then
         write(*,*) '(charge_exchange) ERROR: cx_add_to_turnover was',    &
                    ' called on a thread that holds no cell.'
         write(*,*) '  cx_kc is filled cell by cell; call cx_set_cell(T)', &
                    ' for the cell being judged, on this thread, first.'
         write(*,*) '  Aborting.'
         error stop 1
      endif

      do i = 1, cx_nact
         r     = cx_act(i)
         bound = cx_kc(i)*el_tot(cx_don_el(r))*el_tot(cx_acc_el(r))
         idon  = cx_fvidx(cx_don_el(r), cx_don_stg(r))
         iacc  = cx_fvidx(cx_acc_el(r), cx_acc_stg(r)-1)
         srow(idon) = srow(idon) + bound
         srow(iacc) = srow(iacc) + bound
      enddo
      end subroutine cx_add_to_turnover

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
      ! Oxygen O <-> H+ near-resonant charge exchange. Huang et al. (2023)
      ! Table 4 prints the two rate coefficients with the reactant labels
      ! exchanged: the exp(-227/T) Boltzmann factor is printed on the O+ + H0
      ! row, but IP(O I)=13.6181 eV > IP(H I)=13.5984 eV, so the O0 + H+ -> O+
      ! + H0 ionizing channel is the endothermic one (dE/k = 227.7 K) and must
      ! carry the barrier; O+ + H0 -> O0 + H+ is exothermic and must not.
      ! Assigned here to the physically correct rows (detailed balance then
      ! holds to ~9%, and each direction matches Cloudy c25.00 to 1-3%; as
      ! printed it violates detailed balance by 1.47x at 8000 K). See
      ! docs/HUANG2023_TABLE4_OXYGEN_ERRATUM.md.
      case (13); cx_rate = (1.26d-9*t4**0.517d0 + 4.25d-10*t4**6.69d-3)       &
                           *exp(-227.0d0/T)                                    ! A13 O + H+  (endothermic, ionizing)
      case (14); cx_rate = 2.08d-9*t4**0.405d0 + 1.11d-11*t4**(-0.458d0)       ! A14 O+ + H  (exothermic, recombining)
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
      ! THE TWO DIRECTIONS ARE DIFFERENT CHANNELS, and no detailed-balance
      ! relation is expected between them. Both reach the code through
      ! Huang et al. (2023) Table 4 (footnote b -> Glover & Jappsen 2007),
      ! which prints them as a forward/reverse pair, but their primary
      ! sources are two separate calculations of two separate processes:
      !
      !   B1  He(1s^2) + H+ -> He+ + H         NON-RADIATIVE (collisional)
      !       Kimura, Lane, Dalgarno & Dixson (1993), ApJ 405, 801:
      !       molecular close coupling on full-CI potentials. Their Table 3
      !       tabulates the rate coefficient from 6000 K to 1e5 K ONLY, and
      !       their charge-transfer cross sections start at 20 eV; the
      !       11 eV threshold makes this direction endothermic, so below
      !       6000 K the expression here is pure Arrhenius extrapolation
      !       with no computed point behind it. It is also negligible
      !       there: 1.7e-60 cm^3 s^-1 at 1140 K, measured.
      !
      !   B2  He+ + H -> He + H+ + photon      RADIATIVE charge transfer
      !       Stancil, Lepp & Dalgarno (1998), ApJ 509, 1, Table 1 row (19),
      !       which prints the photon in the exit channel and whose fit
      !       carries the 0.25 approach-probability factor of Stancil &
      !       Zygelman (1996). The underlying calculation is Zygelman,
      !       Dalgarno, Kimura & Lane (1989), Phys. Rev. A 40, 2340, read
      !       here from the published pages. Their Table I, captioned "Rate
      !       coefficients for direct radiative charge transfer and for
      !       total radiative decay in units of 10^-15 cm^3 sec^-1", gives
      !       the DIRECT column 5.36, 4.21, 4.41, 4.50, 4.83, 5.99 at
      !       T = 1, 10, 100, 200, 400, 1000 K. The 0.25 factor multiplies
      !       that DIRECT column -- their process (2), radiative charge
      !       transfer -- and not the "Total" column, which is the total
      !       radiative decay including radiative association to HeH+:
      !       interpolating Direct to 300 K gives 4.67e-15, and 0.25 times
      !       that is 1.17e-15, which is the 1.25e-15 of the fit to the
      !       rounding.
      !       PUBLISHED RANGE 1 - 1000 K: their introduction says they
      !       present the rate coefficient of the radiative processes "for
      !       temperatures up to 1000 K", and Table I stops there. EXHALE
      !       evaluates B2 at the 1140 K molecular base and up into the
      !       wind, so it is extrapolated above its tabulation everywhere
      !       it matters here.
      !       WHAT THE EXTRAPOLATION CAVEAT DOES AND DOES NOT COVER.
      !       Courtney, Forrey, McArdle, Stancil & Babb (2021), ApJ 919, 70
      !       recomputed this exact process -- their Table 1 reaction 9,
      !       He+ + H -> He + H+ + h-nu -- and report that "The smooth
      !       backgrounds for both processes agree very well with the
      !       previous results of Zygelman et al. (1989) and Zygelman &
      !       Dalgarno (1990)", so the basis of this fit is supported by an
      !       independent modern calculation. Two corrections that could
      !       have applied to an extrapolation out of a 1000 K table are
      !       ruled out by them for our regime: the narrow J' = 8 and 11
      !       resonances that separate their LTE and NLTE-ZDL rates become
      !       "negligible above 1000 K", and "Stimulated processes are
      !       shown to only be important above about 50,000 K". EXHALE
      !       evaluates B2 from 1140 K upward, so neither is needed.
      !       Their Figure 3 also puts the recomputed rate coefficient near
      !       1.6e-15 cm^3 s^-1 at 1e3 K and 2.5e-15 at 1e4 K (read from
      !       the figure), against 1.69e-15 and 3.00e-15 from the fit
      !       below: agreement to 6% and 20% across the range this code
      !       uses. The paper gives NO closed-form fit for reaction 9 --
      !       its Table 2 coefficients cover reactions 5, 6 and 3 only, and
      !       reaction 9 lives in that figure and in the UGA Molecular
      !       Opacity Project database -- so there is nothing to swap this
      !       expression for, and it is kept.
      !       What remains a caveat is therefore only the extrapolation
      !       beyond the published tabulation itself, not the resonance or
      !       stimulated-emission physics.
      !
      ! A photon-emitting channel has no collisional reverse, so the ratio
      ! k(B1)/k(B2) is under no obligation to equal 4 exp(-dE/kT); the
      ! factor 1.05e6/T by which it departs from that (921x at 1140 K,
      ! 210x at 5000 K, 105x at 1e4 K, evaluated from the two expressions
      ! below) is a consequence
      ! of comparing two channels, not evidence of a defect in either fit.
      ! The label in Table 4, and in this file before 2026-08-31, dropped
      ! the photon and with it the reason. See
      ! docs/molecular_chemistry_audit_he_rich.md section 4.3.
      !
      ! B1 HAS TWO BRANCHES. Glover & Jappsen (2007), ApJ 666, 1, Table 1
      ! reaction R27, read from the published paper, fit the Kimura table as
      !     k27 = 1.26e-9 T^-0.75 exp(-127,500/T)   for T <= 10,000 K
      !         = 4.0e-37 T^4.74                    for T >  10,000 K
      ! and only the low branch was carried here before 2026-08-31. The two
      ! sources do NOT disagree on that branch: the form below, written with
      ! tr = T/300, expands to 1.2614e-9 T^-0.75 exp(-127,500/T), which is
      ! Glover & Jappsen's low branch to the rounding of their 1.26e-9, so
      ! the Huang et al. (2023) Table 4 route and the original agree and
      ! there is nothing to adjudicate. It is kept verbatim so that the
      ! branch below 1e4 K is unchanged.
      ! The branches join cleanly: 3.661e-18 against 3.648e-18 cm^3 s^-1 at
      ! the switch, 0.4% apart (measured). Above it the low branch runs high
      ! -- 1.9x at 1.1e4 K, 13.1x at 2e4 K, 7.1x at 4e4 K -- so carrying the
      ! published high branch LOWERS B1 there.
      case (24)
         if (T .le. 1.0d4) then
            cx_rate = 1.75d-11*tr**(-0.75d0)*exp(-12.75d0/t4)                  ! B1 He + H+   (non-radiative, Kimura 1993)
         else
            cx_rate = 4.0d-37*T**4.74d0                                        ! B1 high branch, Glover & Jappsen (2007) R27
         endif
      case (25); cx_rate = 1.25d-15*tr**0.25d0                                 ! B2 He+ + H   (RADIATIVE, Stancil 1998)
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
      ! ---------------- Group E: outside Huang Table 4 ----------------
      case (64); cx_rate = cx_o2p_h_scale                                   &
                           *electron_capture_O2p_from_H(T)                   ! E1 O2+ + H0
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
      ! d(CX source)/dx to the dense Jacobian fjac. Each species's density is
      ! a LINEAR function of the unknowns, so for a reaction rate
      ! R = kc*D*A (donor density D, acceptor density A),
      !   dR/dx_k = kc*(dD/dx_k * A + D * dA/dx_k),
      ! added to row idon (+) and row iacc (-), exactly mirroring the
      ! fvec(idon)+=R, fvec(iacc)-=R bookkeeping. n_X (= met_ntot), n_h, n_he
      ! supply the linear-derivative coefficients. he_row_sign and T_cell
      ! carry the same meaning as in cx_add_to_fvec, and the derivative rows
      ! are oriented exactly as the residual rows are, so this stays the
      ! derivative of what that routine assembles.
      subroutine cx_add_to_jac(N_eq, fjac, nm0, nm1, nm2,                 &
                               n_hi, n_hii, n_hei, n_heii, n_heiii,       &
                               n_X, n_h, n_he, he_row_sign, T_cell)
      integer, intent(in)    :: N_eq
      real*8,  intent(inout) :: fjac(N_eq,N_eq)
      real*8,  intent(in)    :: nm0(:), nm1(:), nm2(:)
      real*8,  intent(in)    :: n_hi, n_hii, n_hei, n_heii, n_heiii
      real*8,  intent(in)    :: n_X(:), n_h, n_he
      real*8,  intent(in)    :: he_row_sign, T_cell
      real*8  :: dens(12,0:2), D, A, kc, dR_dn, sdon, sacc
      integer :: i, r, de, ds, ae, as, idon, iacc, m
      integer :: nD, nA, kD(2), kA(2)
      real*8  :: cD(2), cA(2)

      call cx_require_cell(T_cell, 'cx_add_to_jac')

      call cx_reactant_densities(nm0, nm1, nm2, n_hi, n_hii, n_hei,       &
                                 n_heii, n_heiii, dens)

      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r); ds = cx_don_stg(r)
         ae = cx_acc_el(r); as = cx_acc_stg(r)
         D  = dens(de,ds);  A  = dens(ae,as);  kc = cx_kc(i)
         idon = cx_fvidx(de, ds)
         iacc = cx_fvidx(ae, as-1)
         sdon = cx_he_i_boundary_sign(de, ds,   he_row_sign)
         sacc = cx_he_i_boundary_sign(ae, as-1, he_row_sign)
         ! donor-density derivatives: dR = kc * (dD) * A
         call cx_dens_lin(de, ds, n_X, n_h, n_he, nD, kD, cD)
         do m = 1, nD
            dR_dn = kc*cD(m)*A
            fjac(idon,kD(m)) = fjac(idon,kD(m)) + sdon*dR_dn
            fjac(iacc,kD(m)) = fjac(iacc,kD(m)) - sacc*dR_dn
         enddo
         ! acceptor-density derivatives: dR = kc * D * (dA)
         call cx_dens_lin(ae, as, n_X, n_h, n_he, nA, kA, cA)
         do m = 1, nA
            dR_dn = kc*D*cA(m)
            fjac(idon,kA(m)) = fjac(idon,kA(m)) + sdon*dR_dn
            fjac(iacc,kA(m)) = fjac(iacc,kA(m)) - sacc*dR_dn
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

      ! ================================================================
      ! Dedicated He <-> H charge-exchange pair (Huang 2023, Table 4,
      ! group B). Applied by every ionization system that contains He, so it
      ! is not tied to the metal charge-exchange assembly (cx_add_to_fvec).
      !   B1: He0 + H+ -> He+ + H0   (k1, endothermic, exp(-12.75/T4))
      !   B2: He+ + H0 -> He0 + H+   (k2, exothermic)
      ! B1 ionizes He (HeI->HeII) and recombines H (HII->HI); B2 ionizes H
      ! (HI->HII) and recombines He (HeII->HeI). Both directions share the
      ! two rate coefficients from he_h_cx_rates. Gated by
      ! he_h_charge_exchange (default .true.): each routine adds nothing when
      ! the switch is off.
      ! ================================================================

      ! Rate coefficients [cm^3 s^-1] of the He<->H pair at temperature T [K].
      ! k_He0_Hp is B1 (He0+H+), k_Hep_H0 is B2 (He++H0). The formulas live in
      ! one place (cx_rate rows 24,25), so this reuses them.
      ! Non-radiative electron capture by singly ionized helium from atomic
      ! hydrogen,
      !     He+ + H(1s) -> He + H+ ,
      ! rate coefficient in cm^3 s^-1. This is the He+ row of Table 1 of
      ! Kingdon & Ferland (1996), ApJS 106, 205, read from the journal PDF:
      ! parameters a = 7.47e-6 (in units of 1e-9 cm^3 s^-1), b = 2.06,
      ! c = 9.93, d = -3.89 over 6e3 - 1e5 K, in their fitting form (their
      ! eq. 7)
      !     alpha(t4) = a t4^b [1 + c exp(d t4)] ,   t4 = T/1e4 ,
      ! from the calculation of Zygelman et al. (1989).
      !
      ! WHY THIS IS A SEPARATE CHANNEL FROM B2, AND NOT A DUPLICATE OF IT.
      ! Confirmed against the common source, Zygelman, Dalgarno, Kimura &
      ! Lane (1989), read from the published pages. That paper computes both
      ! channels and keeps them apart throughout: its title is "Radiative
      ! and nonradiative charge transfer in He+ + H collisions at low
      ! energy", its Fig. 8 compares "nonradiative (direct) charge transfer
      ! cross sections (circles) and radiative charge transfer cross
      ! sections (triangles)", and its Table I tabulates RATE COEFFICIENTS
      ! for the RADIATIVE processes only, from 1 to 1000 K. So:
      !   - the radiative rate coefficient (Table I, Direct column) is the
      !     one Stancil et al. (1998) fit and the code carries as B2;
      !   - the non-radiative channel exists in that paper only as CROSS
      !     SECTIONS over 1 - 100 eV (their Fig. 8), which is the collision
      !     energy range that maps onto the 6e3 - 1e5 K of the Kingdon &
      !     Ferland fit adopted here.
      ! The two therefore come from different parts of the calculation and
      ! are additive; the total He+ + H removal rate is their sum. Kingdon &
      ! Ferland's own footnote on this row, that radiative charge transfer
      ! may be the most rapid mechanism, is consistent with their tabulated
      ! value being the non-radiative one.
      ! Zygelman et al. also give the physical reason this matters in a
      ! wind: the radiative processes dominate at low collision energy,
      ! while above a few eV the direct (non-radiative) process becomes the
      ! faster removal mechanism, and in the 5-8 eV range the two cross
      ! sections are comparable at about 1e-20 cm^2.
      ! Measured sizes of the channel added here, relative to B2: 0.36x at
      ! 1140 K, 1.16x at 3020 K and 2.99x at 1e4 K, so the total He+ + H
      ! removal rate rises by 1.36x, 2.16x and 3.99x respectively.
      !
      ! VALIDITY AND EXTRAPOLATION. Kingdon & Ferland state, of every fit in
      ! their Table 1, that "the fits to these data are only valid within
      ! the given temperature range; extrapolation significantly outside
      ! this range is not recommended". The range here is 6e3 - 1e5 K, and
      ! the molecular base of interest sits near 1140 K, a factor 5 below
      ! it. The two sides are treated differently, on the shape of the fit:
      !   - BELOW the range the expression is evaluated as fitted. It decays
      !     as t4^2.06 toward low temperature, which is the behavior a
      !     curve-crossing channel should have, so downward extrapolation
      !     under-weights the channel and cannot make it spuriously large.
      !     Clamping to the 6000 K value instead was rejected: that would
      !     freeze the rate at 5.1e-15 and overstate this channel at the
      !     base by about a factor 8.
      !   - ABOVE the range the fit GROWS as t4^2.06 without bound, so the
      !     temperature is capped at the published ceiling. This follows the
      !     treatment already used for the Group E oxygen row in this file.
      ! The hottest cell of the regression matrix is about 1.1e4 K, so the
      ! upper cap is not reached in any current run.
      !
      ! PUBLISHED PRECEDENT for handling a Kingdon & Ferland fit outside its
      ! stated range. Ziegler, U. (2018), A&A 620, A81 section 2.1 states
      ! that an implementation "has to ensure that the implementation does
      ! not suffer from unphysical behavior or unboundedness of the rate
      ! coefficients in the asymptotic limits T, n -> 0, infinity. All rate
      ! coefficients are checked in this respect and, if necessary, are
      ! reasonably modified", and his worked example is a Kingdon & Ferland
      ! (1996) coefficient from the same table as this one: their C+ + H
      ! rate, range [5.5e3 K, 1e5 K], which "Below T = 5390 K ... becomes
      ! negative, which is unphysical. Therefore, k is set to zero below
      ! this threshold. Likewise, an upper floor is set at T = 1e9 K."
      ! So bounding these fits outside their published range on physical
      ! grounds, asymmetrically and by inspection of the fit's own shape, is
      ! established practice rather than something invented here.
      ! The precedent is for the PRACTICE, not for the specific action: his
      ! C+ + H fit turns negative below threshold and is therefore zeroed,
      ! whereas the He+ fit here stays positive and merely decays, so it is
      ! evaluated rather than zeroed. His ceiling (1e9 K) and ours (1e5 K,
      ! the published range boundary) likewise differ.
      double precision function nonradiative_electron_capture_Hep_from_H(T) &
                                result(k)
      real*8, intent(in) :: T
      real*8, parameter  :: a = 7.47d-6, b = 2.06d0
      real*8, parameter  :: c = 9.93d0,  d = -3.89d0
      real*8, parameter  :: T_max = 1.0d5     ! published ceiling
      real*8 :: t4
      t4 = min(T, T_max)/1.0d4
      if (t4 .le. 0.0d0) then
         k = 0.0d0
      else
         k = a*t4**b*(1.0d0 + c*exp(d*t4))*1.0d-9
      endif
      end function nonradiative_electron_capture_Hep_from_H

      subroutine he_h_cx_rates(T, k_He0_Hp, k_Hep_H0)
      ! The two H <-> He rate coefficients, in the same admissible band
      ! [0, cx_kc_max] that cx_set_cell imposes on every other row. Group B
      ! is excluded from cx_act and reaches the systems only through here,
      ! so before 2026-08-31 it was the one pair that skipped that band.
      ! Both values sit far inside it (1e-60 to 4e-15 over 300 - 2e4 K,
      ! measured), so this is a consistency repair and not a numerical
      ! change; one kind of row gets one rule.
      real*8, intent(in)  :: T
      real*8, intent(out) :: k_He0_Hp, k_Hep_H0
      ! He+ + H removal is the SUM of the two channels that do it: the
      ! radiative one of row B2 and the non-radiative one above. Carrying
      ! only B2, as the code did before 2026-08-31, left out a channel of
      ! the same order (see the function's header for the sizes).
      k_He0_Hp = min(max(cx_rate(cx_B1_He0_Hp, T), 0.0d0), cx_kc_max)
      k_Hep_H0 = min(max(cx_rate(cx_B2_Hep_H0, T)                        &
                         + nonradiative_electron_capture_Hep_from_H(T),  &
                         0.0d0), cx_kc_max)
      end subroutine he_h_cx_rates

      ! Residual terms for the systems whose H row (fvec(1)) is written
      ! HI->HII (ionization) positive and He row (fvec(2)) carries the
      ! HeI<->HeII balance.
      ! THE HELIUM REACTANT n_hei IS THE GROUND SINGLET He(1^1S), NOT THE
      ! SUM OVER He I. The Table 4 rate for He + H+ -> He+ + H (Glover &
      ! Jappsen 2007) carries the barrier exp(-12.75/T4), and 12.75e4 K =
      ! 10.99 eV is the ionization-potential difference 24.587 - 13.598 eV
      ! of ground-state helium against hydrogen. He(2^3S) sits 19.82 eV
      ! above the singlet, so the same collision is exothermic for it and
      ! runs at a different rate; that reaction is not carried anywhere in
      ! the code, and passing the summed He I here would charge it to this
      ! rate. Every caller therefore passes the singlet: the systems
      ! without a metastable have no other helium, and the two He 2^3S
      ! systems and the molecular ones pass n_heiSI.
      ! he_row_sign = +1 when that He row is written
      ! HeI->HeII (ionization) positive (System_HeH / System_HeH_metals, and
      ! the He+-production row of System_HeH_mol); he_row_sign = -1 when it is
      ! written HeI-gain positive (the summed HeI row of the TR systems). The
      ! H row is HI->HII positive in all of them. k1 = He0+H+, k2 = He++H0.
      subroutine he_h_cx_fvec(fvec, k1, k2, n_hi, n_hii, n_hei, n_heii,   &
                              he_row_sign)
      real*8              :: fvec(*)
      real*8, intent(in)  :: k1, k2, n_hi, n_hii, n_hei, n_heii, he_row_sign
      real*8 :: R1, R2
      if (.not. he_h_charge_exchange) return
      R1 = k1*n_hei*n_hii       ! He0 + H+ -> He+ + H0
      R2 = k2*n_heii*n_hi       ! He+ + H0 -> He0 + H+
      ! H row (HI->HII positive): B2 ionizes H0 (+R2), B1 recombines H+ (-R1).
      fvec(1) = fvec(1) + R2 - R1
      ! He row: B1 ionizes HeI->HeII (+R1), B2 recombines HeII->HeI (-R2), in
      ! the ionization-positive orientation; he_row_sign flips it for the
      ! HeI-gain (TR) orientation.
      fvec(2) = fvec(2) + he_row_sign*(R1 - R2)
      end subroutine he_h_cx_fvec

      ! Analytic-Jacobian counterpart of he_h_cx_fvec for the systems that
      ! carry a written Jacobian (System_HeH, System_HeH_metals), whose He row
      ! is HeI->HeII positive (he_row_sign = +1). Unknown layout
      ! x1 = n_HII/n_H, x2 = n_HeII/n_He, x3 = n_HeIII/n_He, so the pair
      ! touches rows 1,2 and columns 1,2,3 only. nsz = system size.
      subroutine he_h_cx_jac(nsz, fjac, k1, k2, n_h, n_he,                &
                             n_hi, n_hii, n_hei, n_heii)
      integer, intent(in) :: nsz
      real*8              :: fjac(nsz,nsz)
      real*8, intent(in)  :: k1, k2, n_h, n_he, n_hi, n_hii, n_hei, n_heii
      real*8 :: j1, j2, j3
      if (.not. he_h_charge_exchange) return
      ! j_k = d(R2 - R1)/dx_k (the H-row derivative), with the linear density
      ! derivatives dn_hi/dx1=-n_h, dn_hii/dx1=n_h, dn_heii/dx2=n_he,
      ! dn_hei/dx2 = dn_hei/dx3 = -n_he. The He row (he_row_sign=+1) is minus
      ! the H row, so fjac(2,k) += -j_k.
      j1 = -k2*n_heii*n_h - k1*n_hei*n_h           ! d/dx1
      j2 =  k2*n_he*n_hi  + k1*n_he*n_hii          ! d/dx2
      j3 =  k1*n_he*n_hii                          ! d/dx3
      fjac(1,1) = fjac(1,1) + j1
      fjac(1,2) = fjac(1,2) + j2
      fjac(1,3) = fjac(1,3) + j3
      fjac(2,1) = fjac(2,1) - j1
      fjac(2,2) = fjac(2,2) - j2
      fjac(2,3) = fjac(2,3) - j3
      end subroutine he_h_cx_jac

      ! Residual terms for the advection-correction systems
      ! (System_implicit_adv_HeH / _HeH_TR). Those rows are fraction-
      ! normalized and written neutral-gain positive: fvec(1) tracks n_HI
      ! (per n_H), fvec(2) tracks n_HeI (per n_He = heh_loc*n_H). c1 = dr/v.
      ! xhi/xhii = HI/HII fractions of H; xhei/xheii = HeI/HeII fractions of
      ! He. k1 = He0+H+, k2 = He++H0.
      ! THE HELIUM REACTANT xhei IS THE GROUND SINGLET He(1^1S), as it is in
      ! he_h_cx_fvec and for the same reason: the Table 4 rate carries the
      ! barrier exp(-12.75/T4) = 10.99 eV, the ionization-potential
      ! difference of ground-state helium against hydrogen, and the
      ! metastable's own charge exchange is a different reaction that this
      ! code does not carry. System_implicit_adv_HeH has no metastable and
      ! its neutral helium is the singlet; System_implicit_adv_HeH_TR passes
      ! x(2), the singlet unknown of its row 2.
      subroutine he_h_cx_fvec_adv(fvec, c1, xhi, xhii, xhei, xheii,       &
                                  heh_loc, n_h, k1, k2)
      real*8              :: fvec(*)
      real*8, intent(in)  :: c1, xhi, xhii, xhei, xheii, heh_loc, n_h, k1, k2
      real*8 :: rr1, rr2
      if (.not. he_h_charge_exchange) return
      ! rr1 = R1/n_he = k1*n_hei*n_hii/n_he, rr2 = R2/n_he. Then R/n_h = rr*heh_loc.
      rr1 = k1*xhei*xhii*n_h        ! He0 + H+ -> He+ + H0, per n_he
      rr2 = k2*xheii*xhi*n_h        ! He+ + H0 -> He0 + H+, per n_he
      ! H row (HI-gain positive): B1 makes HI (+R1), B2 destroys HI (-R2); /n_h.
      fvec(1) = fvec(1) + c1*heh_loc*(rr1 - rr2)
      ! He row (HeI-gain positive): B2 makes HeI (+R2), B1 destroys HeI (-R1); /n_he.
      fvec(2) = fvec(2) + c1*(rr2 - rr1)
      end subroutine he_h_cx_fvec_adv

      end module charge_exchange

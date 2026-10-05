      module charge_exchange
      ! Charge-exchange (charge-transfer) reactions from Huang et al. (2023),
      ! ApJ 951, 123, Table 4. Each reaction transfers a single electron.
      !
      ! ONE RULE FOR THE TWO DIRECTIONS OF A PAIR. A rate coefficient that
      ! a source calculated or measured is carried as that source gives it
      ! (published_rate_coefficient). Where the source gives one direction
      ! only, Huang et al. obtained the other "assuming microscopic balance"
      ! (their Section 2.4.2 and Eqs. 1-2) and printed a power-law fit of
      ! it; those fitted reverses are NOT carried. The reverse is computed
      ! here from the forward rate at the same temperature,
      !
      !   k_r = k_f [Q(X^s) Q(Y^(t+1))] / [Q(X^(s+1)) Q(Y^t)]
      !         x exp[(IP(X^s) - IP(Y^t))/kT]
      !
      ! for the forward X^s + Y^(t+1) -> X^(s+1) + Y^t (exothermic, so the
      ! exponent is negative and the product cannot overflow), with Q the
      ! internal partition functions of internal_partition_function and IP
      ! the NIST ionization energies (detailed_balance_ratio; the rows are
      ! listed in cx_reverse_of), times the factor
      ! product_population_factor(f, T, n_e) that corrects it for level
      ! populations of the forward's products that are not Boltzmann (the
      ! note at q_level_energy): 1 at the collisional limit and wherever
      ! the forward's product levels are not known, and set by the local
      ! electron density where the source names them. A forward channel that emits a photon
      ! (radiative charge transfer), or that ends in an excited state which
      ! decays by an allowed line, has no collisional reverse from the
      ! populated states, and its Table 4 reverse is not carried
      ! (cx_A2_Mgp_H). src/tests/charge_exchange_detailed_balance checks
      ! every pair against this relation with the same partition functions.
      !
      ! THE ENERGY OF EACH REACTION. A transfer X^s + Y^(t+1) -> X^(s+1) +
      ! Y^t changes the ionization energy stored in the gas by
      ! dE = IP(Y^t) - IP(X^s). What the translational energy of the gas
      ! gains is dE less the excitation of the product state that leaves as
      ! radiation (an allowed decay, or the photon of a radiative transfer);
      ! an endothermic reaction takes its dE from the gas.
      ! charge_exchange_heating forms that heat for the reactions this
      ! module applies, with the product states of cx_radiated_energy_eV
      ! and, for a product formed in a metastable level the gas does not
      ! hold populated, the forbidden-line loss of
      ! product_manifold_radiation_eV.
      !
      ! THE THRESHOLD LAW every carried rate obeys. A Maxwellian average of
      ! a cross section sigma(E) >= 0 over reactants with thermal internal
      ! populations has d ln(k Q_reactants)/d(-1/kT) = <E_tot>_reacting -
      ! (3/2) kT (Tolman), and <E_tot>_reacting >= max(dE_endothermic, 0),
      ! so k Q_reactants T^(3/2) exp(max(dE_endo,0)/kT) can never decrease
      ! with T. src/tests/charge_exchange_detailed_balance checks it on
      ! every row.
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
      ! He <-> H charge exchange (Group B, He0+H+ and He++H0, plus the
      ! radiative electron capture He2+ + H0 -> He+ + H+ that Table 4 does
      ! not carry) is handled separately so it is available in EVERY
      ! ionization system containing He -- not only the metals systems that
      ! call cx_add_to_fvec. It is gated by he_h_charge_exchange (default
      ! .true.), independent of cx_full, and the same three rate
      ! coefficients feed every system via he_h_cx_rates.
      !
      ! NOTE on C/N/O: Table 4 contains only the neutral<->singly-ionized
      ! charge transfer for C, N, O (rows A13-A18). The captures by the
      ! doubly ionized stages from atomic hydrogen that the earlier Kingdon
      ! & Ferland implementation carried are back for O and N (group E,
      ! below, from Barragan et al. 2006); C2+ + H0 is not carried.
      !
      ! Group E is not from Table 4. It carries the electron captures by
      ! doubly ionized oxygen and nitrogen from atomic hydrogen, O2+ + H0 ->
      ! O+ + H+ (E1), whose absence from that table leaves the O III profile
      ! wrong by decades in the region a transit probes, and N2+ + H0 -> N+
      ! + H+ (E2), both from the same calculation (Barragan et al. 2006).
      ! Each is ON by default at the published rate (cx_o2p_h_scale = 1,
      ! cx_n2p_h_scale = 1); "cx_O2p_H <scale>" and "cx_N2p_H <scale>" in
      ! metals.inp rescale them, and 0 leaves the row out exactly. cx_full
      ! does NOT control them -- cx_full means "all of Table 4", and these
      ! rows are not in Table 4. See electron_capture_O2p_from_H and
      ! electron_capture_N2p_from_H below.
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
      ! the singlet density. The one group C row carried whose reactant is
      ! the helium atom, C2 (Si+ + He -> Si + He+, as Huang et al. print it),
      ! carries exp(-19.1/T4), and 19.07e4 K = 16.436 eV = 24.587 - 8.152 eV
      ! is the ionization-potential difference of GROUND-STATE helium
      ! against silicon. He(2^3S) lies 19.82 eV above the singlet and its
      ! ionization potential is 4.77 eV, so the same collision is
      ! exothermic for it and runs at a rate this code carries nowhere;
      ! charging the metastable to a ground-state rate would be a reaction
      ! that is not in the set. The three metal systems and the constrained
      ! equilibrium therefore all pass the singlet.
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
      ! S_0 = -u_0, S_1 = u_0 - u_1 and S_2 = u_1.
      !
      ! WHO OWNS THE RATE COEFFICIENTS. cx_kc and cx_metal_base are
      ! thread-local and hold ONE cell: cx_set_cell(T, n_e) fills cx_kc for
      ! that cell's temperature and electron density (the reverses depend on
      ! the level populations the electrons maintain) and records the
      ! temperature in cx_cell_T. Every caller of
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
                               n_melem, n_mion, mion_elem, mion_stage
      use global_parameters, only: kb_eV, kb_erg, hp_erg, c_light,        &
                                   e_th_HI, e_th_HeI, e_th_HeII, erg2eV
      ! The Maxwell-averaged collision-strength constant, defined once.
      use Cooling_Coefficients, only: coll_rate_prefactor

      implicit none
      private

      public :: cx_init, cx_set_cell, cx_add_to_fvec, cx_add_to_jac, cx_full
      public :: cx_metal_base, cx_add_to_turnover
      ! The volumetric heat the applied reactions deliver to the gas.
      public :: charge_exchange_heating
      ! The reaction set as stoichiometric number-density sources, for a
      ! solver whose unknowns ARE stage densities rather than boundary
      ! flows (the transported ionization stages).
      public :: charge_exchange_stage_sources
      ! Dedicated He <-> H charge-exchange pair (Group B), available in every
      ! system with He, independent of cx_full.
      public :: he_h_charge_exchange, he_h_cx_rates
      public :: he_h_cx_fvec, he_h_cx_jac, he_h_cx_fvec_adv
      ! Group E: the O2+ + H and N2+ + H electron captures that Table 4
      ! omits, each with its metals.inp scale.
      public :: cx_o2p_h_scale, electron_capture_O2p_from_H
      public :: cx_n2p_h_scale, electron_capture_N2p_from_H
      ! The non-radiative He+ + H0 channel and He2+ + H0 -> He+ + H+ +
      ! photon, applied with the He <-> H pair (read by the gates).
      public :: nonradiative_electron_capture_Hep_from_H
      public :: radiative_electron_capture_Hepp_from_H
      ! The five-level populations behind the density dependence of the
      ! derived reverses (read by the detailed-balance gate).
      public :: metastable_level_populations, metastable_excess_radiation
      public :: lp_NI, lp_SII, lp_CI

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

      ! He2+ + H0 -> He+(1s) + H+ + photon: rate coefficient from the
      ! function radiative_electron_capture_Hepp_from_H (source, range and
      ! uncertainty written there).

      ! Reaction rows carried: 63 from Table 4 (23 A + 2 B + 6 C + 32 D) plus
      ! the two Group-E rows that are NOT in Table 4 (E1, O2+ + H0; E2,
      ! N2+ + H0).
      integer, parameter :: n_cxreac = 65
      ! Row id of the Group-E reaction, and its scale factor.
      !   Source  : Barragan, Errea, Mendez, Rabadan & Riera (2006), ApJ 636,
      !             544, Table 1 column k1 -- molecular close-coupling on MELD
      !             multireference-CI wavefunctions; the calculation CHIANTI
      !             v11 (Dufresne et al. 2024) adopts and Cloudy fits.
      !   Validity: 1e2 - 1e5 K, the calculated range. Outside it the rate is
      !             held at the endpoint, not extrapolated (see
      !             electron_capture_O2p_from_H).
      !   Default : 1.0 = the published rate, made the default on 2026-08-30
      !             (docs/Update_EXHALE_stage1.pdf section 107, measurement in section
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
      ! Row id of E2, N2+ + H0 -> N+ + H+, and its scale factor.
      !   Source  : Barragan et al. (2006), the same calculation as E1, their
      !             Table 2 column k4 (electron_capture_N2p_from_H).
      !   Validity: 1e2 - 1e5 K, held at the endpoints outside.
      !   Default : 1.0 = the published rate. "cx_N2p_H 0" in metals.inp
      !             leaves the row out of cx_act altogether.
      integer, parameter :: cx_E2_N2p_H0 = 65
      real*8, save :: cx_n2p_h_scale = 1.0d0

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
      ! -- Group E: NOT in Table 4 (E1 O2+ + H0, E2 N2+ + H0) --
        cx_H,   cx_H ]

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
        0, 0 ]

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
        iel_O,  iel_N ]

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
      ! -- E: the acceptors are O and N at stage 2 (O III -> O II,
      !    N III -> N II) --
        2, 2 ]

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
      ! -- E (2): not in Table 4; gated on cx_o2p_h_scale and
      !    cx_n2p_h_scale, not on cx_full --
        .false., .false. ]

      ! Switch: .false. = Group A only (metal + H); .true. = full Table 4.
      logical, save :: cx_full = .false.

      ! THE DETAILED-BALANCE REVERSES. cx_reverse_of(r) is the forward row
      ! whose rate row r is computed from (module header); 0 marks a row
      ! carried as its source gives it. Each forward named here is the
      ! exothermic direction of a non-radiative channel that its source
      ! (read) calculated, measured or estimated and that ends in the
      ! populated manifold of the product (the note at q_level_energy), and
      ! each reverse replaces a Huang et al. (2023) power-law fit of
      ! microscopic balance (or, for A3 and A7, Kingdon & Ferland's 1996
      ! refit of the same relation):
      !   A3 <- A4, A7 <- A8, A11 <- A12, A17 <- A18 (Kingdon & Ferland
      !   1996, Butler & Dalgarno 1979), the metal + metal pairs of
      !   Prasad & Huntress (1980) and Lavvas et al. (2014),
      !   D4 <- D3, ..., D24 <- D23 (not D18) and D32 <- D31, A20 <- A19
      !   (Zhao et al. 2005: S + H+ ends in S+ 4S, 2D and 2P, all in the
      !   populated manifold of S+), and D26 <- D25 (Chenel et al. 2010
      !   computed C+(2P) + S(3P) -> C(3P) + S+(4S); the 0.90 eV defect is
      !   below the first excited term of either product, C 1D 1.26 eV and
      !   S+ 2D 1.84 eV, NIST, so the channel ends in the ground terms).
      ! Carried as their sources give both directions: A13/A14 (Stancil et
      ! al. 1999) and A15/A16 (Stancil et al. 1998). Carried as Huang et
      ! al. print them, because the original of the forward is unknown:
      ! D1/D2, which Table 4 attributes to Satta et al. (2013), a paper
      ! (read) that computes Si + He+ only. Group B has its own routine,
      ! he_h_cx_rates, which applies the same relation.
      integer, parameter :: cx_reverse_of(n_cxreac) = [                    &
      ! -- A1..A23 --
        0,  0,  4,  0,  0,  0,  8,  0,  0,  0, 12,  0,                      &
        0,  0,  0,  0, 18,  0,  0, 19,  0,  0,  0,                          &
      ! -- B1, B2 (he_h_cx_rates) --
        0,  0,                                                              &
      ! -- C1..C6 --
        0,  0,  0,  0,  0,  0,                                              &
      ! -- D1..D32 --
        0,  0,  0, 34,  0, 36,  0, 38,  0, 40,  0, 42,  0, 44,  0, 46,      &
        0,  0,  0, 50,  0, 52,  0, 54,  0, 56,  0,  0,  0,  0,  0, 62,      &
      ! -- E1, E2 --
        0,  0 ]

      ! TEN TABLE 4 ROWS ARE NOT CARRIED. Each is Huang et al.'s
      ! microscopic-balance reverse of a forward whose source says it has
      ! no collisional reverse from the populated states:
      !   A2, Mg+ + H -> Mg + H+. Its forward A1 is the calculation of
      !       Allan, Clegg, Dickinson & Flower (1988), MNRAS 235, 1245, in
      !       which "The Mg+ ions are formed in the 3p 2Po state and 4s 2S
      !       state, which cascades to 3p 2Po", both emptied by allowed
      !       lines; the reverse of those channels starts from states that
      !       are not populated, and Kingdon & Ferland (1996) give no
      !       reverse for this reaction because it is "negligible due to
      !       the energetics".
      !   A10, Si+ + H -> Si + H+. Its forward A9 is the channel
      !       Si(3P) + H+ -> Si+(4P) + H of Kimura et al. (1996), ApJ 473,
      !       1114; the reverse of that channel starts from Si+ 3s3p2 4P,
      !       5.31 eV up and emptied by the Si II] 2335 A lines, and the
      !       channel to the ground term, whose reverse would start from
      !       the populated Si+, is "much smaller" in their calculation.
      !   C4, C+ + He -> C + He+. Its forward C3 is the non-radiative
      !       calculation of Kimura, Dalgarno, Chantranupong, Li, Hirsch &
      !       Buenker (1993), ApJ 417, 812, whose charge-transfer channel
      !       is "[He(1S) + C+(2D)]": the product is C+ 2s2p2 2D, 9.29 eV
      !       up, which decays to the ground term by an allowed line, so the
      !       reverse of that channel starts from a state that is not
      !       populated.
      !   C6, O+ + He -> O + He+. Its forward C5 is RADIATIVE charge
      !       transfer, O(3P) + He+ -> O+ + He + photon (Zhao et al. 2004,
      !       ApJ 615, 1063, Table 3 caption), which has no collisional
      !       reverse.
      !   A22, Na+ + H -> Na + H+. Its forward A21 is, for the part that
      !       is not radiative, the channel "Na(3s) + H+ -> Na+ + H (n=2)
      !       - 1.682 eV" of Dutta, Nordlander, Kimura & Dalgarno (2001),
      !       Phys. Rev. A 63, 022709 (read): the product is H(n=2), which
      !       decays by Lyman alpha, so the reverse of that channel starts
      !       from H(n=2), not from the ground-state H this row would be
      !       applied to; the radiative part (Watanabe et al. 2002) has no
      !       collisional reverse.
      !   A6, Fe+ + H -> Fe + H+, and D18, Fe+ + O -> Fe + O+. Their
      !       forwards A5 and D17 are the measurements of Rutherford & Vroom
      !       (1972), J. Chem. Phys. 57, 3091 (read), whose cross sections
      !       are large (about 1e-15 cm^2) from 500 eV down to their lowest
      !       energies. Their text: "only for cases of very near resonance
      !       are the cross sections large at low energies", "Iron, with its
      !       multitude of excited ion states allows a near resonance
      !       situation to exist for all the reactions shown in Figs. 1 and
      !       2", and "the state of excitation of the products could not be
      !       determined". With defects of 5.70 and 5.72 eV the near-resonant
      !       Fe+ terms are 3d6(5D)4p z4F and z4D (5.48 - 5.62 eV, NIST),
      !       which decay by allowed lines, so the reverse of the measured
      !       channel starts from states that are not populated, and Huang
      !       et al.'s microscopic-balance reverses of a ground-state channel
      !       are not carried (cx_radiated_energy_eV gives the heat).
      !   D28, Fe+ + N -> Fe + N+. The same paper's N+ + Fe (D27) falls
      !       under the same statement of near resonance, which for its
      !       6.63 eV defect is an odd-parity Fe+ term near 6.6 eV that
      !       decays by allowed lines; not carried for the same reason.
      !   C2, Si+ + He -> Si + He+. Its forward C1 is the calculation of
      !       Satta, Grassi, Gianturco, Yakovleva & Belyaev (2013), MNRAS
      !       436, 2722 (read): "the only exit channels which could have
      !       non-adiabatic crossings with the initial PEC correspond to
      !       electronically excited Si+ atoms", Rydberg-like Si+* at
      !       about 1.19e5 - 1.32e5 cm^-1 (their Figs. 1 and 2) that "could
      !       then become stabilized by radiative emission into the final
      !       ground electronic states of the Si cation"; the reverse of
      !       that channel starts from Si+* and is not carried.
      !   D30, Na+ + Ca -> Na + Ca+. Its forward D29 is the measurement of
      !       Kwolek, Goodman, Slayton, Bluemel, Wells, Narducci & Smith
      !       (2019), Phys. Rev. A 99, 052703 (read). For the ground-state
      !       entrance channel Na[S] + Ca+[S] "there exists only one
      !       lower-energy exit channel, Ca[S] + Na+[S]", reachable only by
      !       radiative charge transfer, predicted at about 1e-16 cm^3 s^-1,
      !       whereas they measure 1e-11 to 1e-9 with a collision-energy
      !       threshold, and "we predict that the dominant process is
      !       endothermic", into "Ca[4s4p 3P0] + Na+[2p6 1S]" (symbols as
      !       typeset). Ca 4s4p 3P decays by the 657 nm intercombination
      !       line, so the reverse of the measured channel starts from a
      !       state that is not populated. (As the reverse of a ground-state
      !       channel it also exceeded the Langevin capture rate of Na+ +
      !       Ca above about 7e3 K.)
      integer, parameter :: cx_A2_Mgp_H  = 2
      integer, parameter :: cx_A6_Fep_H  = 6
      integer, parameter :: cx_A10_Sip_H = 10
      integer, parameter :: cx_A22_Nap_H = 22
      integer, parameter :: cx_C4_Cp_He  = 29
      integer, parameter :: cx_C6_Op_He  = 31
      integer, parameter :: cx_D18_Fep_O = 49
      integer, parameter :: cx_C2_Sip_He = 27
      integer, parameter :: cx_D28_Fep_N = 59
      integer, parameter :: cx_D30_Nap_Ca = 61

      ! IONIZATION ENERGIES [eV] of stage s -> s+1, NIST Atomic Spectra
      ! Database, "Ionization Energies" form (retrieved 2026-09-24), by the
      ! element codes of this module (canonical metals 1..10, then H, He).
      ! The H I, He I and He II thresholds are the named constants of
      ! global_parameters, where each is defined once. Hydrogen has no
      ! second entry (0 is never read).
      real*8, parameter :: cx_ionization_energy(12,0:1) = reshape( [      &
      ! stage 0: C, O, N, Mg, Si, Ca, Na, K, S, Fe, H, He
        11.2602880d0, 13.618055d0, 14.53413d0, 7.646236d0, 8.15168d0,     &
        6.1131549210d0, 5.13907696d0, 4.34066373d0, 10.3600167d0,         &
        7.9024681d0, e_th_HI, e_th_HeI,                                   &
      ! stage 1
        24.383143d0, 35.12112d0, 29.60125d0, 15.035271d0, 16.34585d0,     &
        11.871719d0, 47.28636d0, 31.62500d0, 23.33788d0, 16.19921d0,      &
        0.0d0, e_th_HeII ], [12,2] )

      ! THE LEVELS OF THE INTERNAL PARTITION FUNCTIONS. For each species,
      ! the NIST Atomic Spectra Database levels (energies in cm^-1 as NIST
      ! prints them, retrieved 2026-09-24; weight 2J+1) OF THE GROUND
      ! PARITY BELOW THE LOWEST LEVEL OF OPPOSITE PARITY: the ground term,
      ! its fine structure and the metastable terms, i.e. every level that
      ! has no electric-dipole decay channel. The cut-off is physical, not
      ! numerical. A level with an allowed decay (A ~ 1e7-1e9 s^-1) is
      ! emptied radiatively far faster than collisions at the densities of
      ! these atmospheres refill it, whereas the metastable terms (C I 1D,
      ! Fe II a4F, Ca II 3d 2D, ...) are held near their Boltzmann
      ! populations by collisions; the partition function is therefore the
      ! sum over the populated manifold.
      !
      ! WHAT DETAILED BALANCE REQUIRES WHEN THE LEVELS ARE NOT THERMAL.
      ! Microscopic reversibility holds level by level, g_a k(a->c) =
      ! g_c k(c->a) exp(-(E_c - E_a)/kT) for any two levels a, c of the
      ! reactant and product pairs, whatever the level populations are. A
      ! thermal rate coefficient is a population-weighted sum of these, so
      ! a function of T alone exists only when the populations are
      ! themselves functions of T. Write the forward as the code applies
      ! it (the source's rate for every level of the reactant) with a
      ! branching b_c into the product levels c, and let p_c be the actual
      ! fractional populations of the product species. The reverse then is
      !   k_r = k_f [Q(X^s) Q(Y^(t+1)) / Q(Y^t)] exp(-dE/kT)
      !         x sum_c p_c b_c exp(E_c/kT)/g_c ,
      ! so the product's partition function enters as the inverse of that
      ! sum. With Boltzmann populations p_c = g_c exp(-E_c/kT)/Q it is Q
      ! for ANY branching; with thermal branching it is Q for ANY
      ! populations. Only when the metastables are empty AND the forward
      ! lands in the ground term alone does it become the ground-term
      ! weight G (the low-density limit, where it is G/p_ground).
      !
      ! THE REVERSE AT ANY DENSITY. Write p_c = beta_c g_c exp(-E_c/kT)/Q,
      ! beta_c the departure coefficient of product level c from Boltzmann.
      ! The sum above is then (1/Q) sum_c b_c beta_c, and
      !   k_r(T, n_e) = k_r(T)|thermal x F,   F = sum_c b_c beta_c ,
      ! with k_r|thermal the relation of the module header: the reverse per
      ! product particle is the thermal one times the branching-weighted
      ! departure coefficient of the product levels the forward reaches
      ! (product_population_factor). F = 1 at the collisional limit
      ! (beta = 1), and F = 1 at every density when the branching is
      ! thermal (b_c = g_c exp(-E_c/kT)/Q), which is the statement that
      ! every level of the product manifold reacts back with one rate
      ! coefficient; this is carried for every forward whose source does
      ! not name its product levels (A4, A8, A12, the group D estimates).
      ! Where the source names them, F follows the populations:
      !   A18 N+ + H -> N(4S) + H+ (the only N term the 0.94 eV reaches):
      !       F = beta(N I 4S), which rises to Q(N I)/4 = 1.18 at 1e4 K as
      !       n_e -> 0;
      !   A19 S + H+ -> S+(4S, 2D, 2P) + H, the product terms of Zhao et
      !       al. (2005) Table III, dominated by S+ 2P (3.04 eV up):
      !       F = sum_tau b_tau beta_tau(S II), which falls to
      !       b_4S Q(S II)/4, 7e-5 at 1e4 K, as n_e -> 0, because an S+ ion
      !       the gas does not hold in 2D or 2P can only return through the
      !       small S+(4S) channel;
      !   D25 S + C+ -> S+(4S) + C(3P) (ground terms by energetics):
      !       F = beta(C I 3P) beta(S II 4S).
      ! The departure coefficients are those of
      ! metastable_level_populations, the five-level statistical
      ! equilibrium of N I, S II and C I (exactly the levels of this table)
      ! at the cell's (T, n_e), with the CHIANTI 11.0.2 A-values and
      ! electron collision strengths. MEASURED against the full CHIANTI
      ! model (all levels, electron and proton collisions, ChiantiPy 0.15.2)
      ! the term departures agree within 0.1% at 3e3 and 1e4 K for 1e2 <=
      ! n_e <= 1e10 cm^-3, and within 2% at 2e4 K, where the levels above
      ! the five cascade into them. CHIANTI carries no collisions with
      ! neutral hydrogen, so in a neutral gas of low n_e the populations err
      ! toward the low-density limit. VALIDITY: every density; the
      ! remaining uncertainty is the branching assumed where the source
      ! gives none. The ground-branching limit of those pairs, MEASURED on
      ! CHIANTI's populations at the densities of the wasp_full wind (n_e =
      ! 3e8 - 4e9 cm^-3), is 0.99 of this Q for Fe II and 1.00 for Fe III,
      ! S II, N II, O I and C I at 1e4 K and n_e = 1e9, and 0.94 (Fe II) at
      ! n_e = 1e7; at 2e4 K it is 0.90 (Fe II), 0.94 (Fe III) (ratios on
      ! CHIANTI's own level set). The exception is Ca II, whose 3d 2D
      ! CHIANTI leaves below Boltzmann there (G/p_ground = 0.81 Q at 1e4 K,
      ! n_e = 1e9; 0.53 Q at 2e4 K); no carried pair is derived through it.
      !
      ! The populated-manifold Q is also finite and independent
      ! of any Rydberg or pressure cut-off, which an LTE sum over all bound
      ! levels is not: for Na I, K I and Ca I the allowed levels would
      ! multiply Q by 1.4-3.6 at 1e4 K and by 3-17 at 2e4 K and diverge
      ! with the Rydberg series. Left out as well are the levels above the
      ! cut whose only E1 decays are intercombination lines (Mg I 3s3p 3Po,
      ! Ca I 4s4p 3Po and 3d4s 1D, the Fe I septets and the levels that
      ! decay only to them): their population lies between zero and
      ! Boltzmann depending on the density, and thermally populated they
      ! would raise Q(Mg I), Q(Ca I) and Q(Fe I) by 39%, 122% and 32% at
      ! 1e4 K. The cut level of each species is
      ! named in its comment (term and energy in cm^-1). Species that no
      ! pair reads are absent (q_first = 0) and are refused by
      ! internal_partition_function.
      integer, parameter :: n_q_level = 206
      real*8, parameter :: q_level_energy(n_q_level) = [                  &
      ! C I, 3 terms: 3P 1D 1S; cut 5So 33735.121
         0.0000000d0, 16.4167130d0, 43.4134567d0, 10192.657d0, 21648.030d0, &
      ! C II, 1 term: 2Po; cut 4P 43002.8
         0.000000d0, 63.395090d0, &
      ! O I, 3 terms: 3P 1D 1S; cut 5So 73768.200
         0.000d0, 158.265d0, 226.977d0, 15867.862d0, 33792.583d0, &
      ! O II, 3 terms: 4So 2Do 2Po; cut 4P 119837.21
         0.00d0, 26810.55d0, 26830.57d0, 40468.01d0, 40470.00d0, &
      ! O III, 3 terms: 3P 1D 1S; cut 5So 60324.79
         0.000d0, 113.178d0, 306.174d0, 20273.27d0, 43185.74d0, &
      ! N I, 3 terms: 4So 2Do 2Po; cut 4P 83284.070
         0.000d0, 19224.464d0, 19233.177d0, 28838.920d0, 28839.306d0, &
      ! N II, 3 terms: 3P 1D 1S; cut 5So 46784.6
         0.0d0, 48.7d0, 130.8d0, 15316.2d0, 32688.8d0, &
      ! Mg I, 1 term: 1S; cut 3Po 21850.405
         0.000d0, &
      ! Mg II, 1 term: 2S; cut 2Po 35669.31
         0.00d0, &
      ! Mg III, 1 term: 1S; cut 3Po 425640.3
         0.0d0, &
      ! Si I, 3 terms: 3P 1D 1S; cut 5So 33326.053
         0.000d0, 77.115d0, 223.157d0, 6298.850d0, 15394.370d0, &
      ! Si II, 1 term: 2Po; cut 4P 42824.29
         0.00d0, 287.24d0, &
      ! Si III, 1 term: 1S; cut 3Po 52724.69
         0.00d0, &
      ! Ca I, 1 term: 1S; cut 3Po 15157.901
         0.000d0, &
      ! Ca II, 2 terms: 2S 2D; cut 2Po 25191.51
         0.00d0, 13650.19d0, 13710.88d0, &
      ! Na I, 1 term: 2S; cut 2Po 16956.17025
         0.00000d0, &
      ! Na II, 1 term: 1S; cut 3Po 264924.32
         0.00d0, &
      ! K I, 1 term: 2S; cut 2Po 12985.185724
         0.0000d0, &
      ! K II, 1 term: 1S; cut 3Po 162502.7
         0.0d0, &
      ! S I, 3 terms: 3P 1D 1S; cut 5So 52623.605
         0.0000d0, 396.05648d0, 573.59573d0, 9238.6090d0, 22179.9542d0, &
      ! S II, 3 terms: 4So 2Do 2Po; cut 4P 79395.52
         0.00d0, 14852.98d0, 14884.75d0, 24524.87d0, 24571.56d0, &
      ! Fe I, 5 terms: a5D a5F a3F a5P a3P2; cut z7Do 19350.891
         0.000d0, 415.933d0, 704.007d0, 888.132d0, 978.074d0, &
         6928.268d0, 7376.764d0, 7728.060d0, 7985.785d0, 8154.714d0, &
         11976.239d0, 12560.934d0, 12968.554d0, 17550.181d0, 17726.988d0, &
         17927.382d0, 18378.186d0, &
      ! Fe II, 24 terms: a6D a4F a4D a4P ...; cut z6Do 38458.9934
         0.0000d0, 384.7872d0, 667.6829d0, 862.6118d0, 977.0498d0, &
         1872.5998d0, 2430.1369d0, 2837.9807d0, 3117.4877d0, 7955.3186d0, &
         8391.9554d0, 8680.4706d0, 8846.7837d0, 13474.4474d0, 13673.2045d0, &
         13904.8604d0, 15844.6485d0, 16369.4098d0, 18360.6399d0, 18886.773d0, &
         20340.2461d0, 20516.9534d0, 20805.7632d0, 20830.5534d0, 21251.5833d0, &
         21307.999d0, 21430.3564d0, 21581.6151d0, 21711.8963d0, 21812.0454d0, &
         22409.8178d0, 22637.1950d0, 22810.3459d0, 22939.3512d0, 23031.2829d0, &
         23317.6351d0, 25428.7893d0, 25787.5816d0, 25805.3270d0, 25981.6451d0, &
         26055.4120d0, 26170.1810d0, 26352.7671d0, 26932.7348d0, 27314.9183d0, &
         27620.4033d0, 30388.5451d0, 30764.4737d0, 31364.4554d0, 31368.453d0, &
         31387.9788d0, 31483.1979d0, 31811.814d0, 31999.049d0, 32875.6409d0, &
         32909.883d0, 33466.4943d0, 33501.2957d0, 36126.412d0, 36252.931d0, &
         37227.277d0, 38164.220d0, 38214.503d0, &
      ! Fe III, 21 terms: 5D a3P 3H a3F ...; cut 7Po 82002.12
         0.00d0, 436.19d0, 738.87d0, 932.52d0, 1027.26d0, &
         19404.42d0, 20051.36d0, 20300.59d0, 20481.72d0, 20688.12d0, &
         21208.12d0, 21461.92d0, 21699.79d0, 21857.08d0, 24558.44d0, &
         24941.12d0, 25142.12d0, 30089.69d0, 30355.88d0, 30715.98d0, &
         30725.46d0, 30857.69d0, 30886.19d0, 34811.9d0, 35803.01d0, &
         41000.31d0, 42896.27d0, 49149.58d0, 49577.27d0, 50184.97d0, &
         50276.19d0, 50295.10d0, 50411.95d0, 57221.20d0, 63425.74d0, &
         63466.85d0, 63487.32d0, 63494.73d0, 63495.00d0, 66465.04d0, &
         66523.30d0, 66592.08d0, 69696.16d0, 69747.95d0, 69788.52d0, &
         69837.15d0, 69838.15d0, 70694.40d0, 70725.46d0, 70729.12d0, &
         73728.05d0, 73849.35d0, 73936.34d0, 76957.02d0, 77044.67d0, &
         77075.66d0, 77102.66d0, 79840.34d0, 79845.04d0, 79860.75d0, &
      ! H I, 1 term: 2S; cut 2Po 82258.9191133
         0.0000000000d0, &
      ! H II: one state
         0.0d0, &
      ! He I: 1s2 1S only (He 2^3S is carried as its own species)
         0.0000d0, &
      ! He II, 1 term: 2S; cut 2Po 329179.2939406
         0.0000d0, &
      ! He III: one state
         0.0d0 ]
      integer, parameter :: q_level_weight(n_q_level) = [                 &
         1, 3, 5, 5, 1, &
         2, 4, &
         5, 3, 1, 5, 1, &
         4, 6, 4, 4, 2, &
         1, 3, 5, 5, 1, &
         4, 6, 4, 2, 4, &
         1, 3, 5, 5, 1, &
         1, &
         2, &
         1, &
         1, 3, 5, 5, 1, &
         2, 4, &
         1, &
         1, &
         2, 4, 6, &
         2, &
         1, &
         2, &
         1, &
         5, 3, 1, 5, 1, &
         4, 4, 6, 2, 4, &
         9, 7, 5, 3, 1, 11, 9, 7, 5, 3, 9, 7, 5, 7, 5, 3, &
         5, &
         10, 8, 6, 4, 2, 10, 8, 6, 4, 8, 6, 4, 2, 6, 4, 2, &
         10, 8, 4, 2, 12, 6, 10, 6, 14, 4, 12, 10, 8, 4, 2, 10, &
         8, 6, 4, 6, 12, 4, 10, 8, 6, 12, 10, 2, 8, 6, 10, 8, &
         4, 2, 6, 8, 6, 8, 14, 12, 10, 8, 4, 6, 2, 6, 4, &
         9, 7, 5, 3, 1, 5, 13, 11, 9, 3, 1, 9, 7, 5, 11, 9, &
         7, 7, 13, 5, 3, 7, 9, 1, 5, 5, 7, 1, 3, 5, 9, 7, &
         5, 9, 13, 11, 9, 7, 5, 7, 5, 3, 9, 1, 3, 7, 5, 11, &
         7, 9, 5, 3, 1, 7, 5, 3, 5, 15, 13, 11, &
         2, &
         1, &
         1, &
         2, &
         1 ]
      integer, parameter :: q_first(12,0:2) = reshape( [                   &
           1,   8,  23,  33,  36,  44,  48,  50,  52,  62, 202, 204, &
           6,  13,  28,  34,  41,  45,  49,  51,  57,  79, 203, 205, &
           0,  18,   0,  35,  43,   0,   0,   0,   0, 142,   0, 206 ], [12,3] )
      integer, parameter :: q_last(12,0:2) = reshape( [                   &
           5,  12,  27,  33,  40,  44,  48,  50,  56,  78, 202, 204, &
           7,  17,  32,  34,  42,  47,  49,  51,  61, 141, 203, 205, &
          -1,  22,  -1,  35,  43,  -1,  -1,  -1,  -1, 201,  -1, 206 ], [12,3] )

      ! THE FIVE-LEVEL MODELS OF N I, S II AND C I (the levels of
      ! q_level_energy for these species, in the CHIANTI 11.0.2 level
      ! order). Data:
      !   energies: the NIST levels of q_level_energy, so that the Boltzmann
      !     limit of these models is exactly that table's partition function;
      !   A-values [s^-1]: CHIANTI 11.0.2 .wgfa files, whose sources are
      !     N I Tachiev & Froese Fischer (2002), A&A 385, 716; S II Tayal &
      !     Zatsarinny (2010), ApJS 188, 32 (E2 and M1 from Zatsarinny);
      !     C I Hibbert et al. (1993), A&AS 99, 179 (forbidden lines);
      !   effective collision strengths (electrons): the CHIANTI 11.0.2
      !     .scups files (N I Tayal 2006, ApJS 163, 207; S II Tayal &
      !     Zatsarinny 2010; C I Wang et al. 2013, Phys. Rev. A 87, 012704,
      !     with Zatsarinny's fine-structure split), descaled by ChiantiPy
      !     0.15.2 onto the temperatures lp_T (DERIVED; the tables are
      !     interpolated linearly in ln T and held at the ends).
      ! Pairs are ordered (1,2) (1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (3,4)
      ! (3,5) (4,5) in CHIANTI level numbers; the upper level of a pair is
      ! the one of higher energy (in N I level 2, 2D3/2, lies above level
      ! 3, 2D5/2). A pair with no radiative line in the file has A = 0.
      integer, parameter :: lp_NI = 1, lp_SII = 2, lp_CI = 3
      integer, parameter :: n_lp_T = 12, n_lp_pair = 10
      real*8, parameter :: lp_T(n_lp_T) = [ 1.0d3, 1.5d3, 2.0d3, 3.0d3,  &
         5.0d3, 7.0d3, 1.0d4, 1.5d4, 2.0d4, 3.0d4, 5.0d4, 1.0d5 ]
      integer, parameter :: lp_lo(n_lp_pair) = [1,1,1,1,2,2,2,3,3,4]
      integer, parameter :: lp_hi(n_lp_pair) = [2,3,4,5,3,4,5,4,5,5]
      ! Level energies [cm^-1] and statistical weights, CHIANTI order.
      real*8, parameter :: lp_E(5,3) = reshape( [                        &
      ! N I 4S3/2, 2D3/2, 2D5/2, 2P1/2, 2P3/2
         0.000d0, 19233.177d0, 19224.464d0, 28838.920d0, 28839.306d0,     &
      ! S II 4S3/2, 2D3/2, 2D5/2, 2P1/2, 2P3/2
         0.00d0, 14852.98d0, 14884.75d0, 24524.87d0, 24571.56d0,          &
      ! C I 3P0, 3P1, 3P2, 1D2, 1S0
         0.0000000d0, 16.4167130d0, 43.4134567d0, 10192.657d0,            &
         21648.030d0 ], [5,3] )
      integer, parameter :: lp_g(5,3) = reshape( [                       &
         4, 4, 6, 2, 4,   4, 4, 6, 2, 4,   1, 3, 5, 5, 1 ], [5,3] )
      ! The term each level belongs to: 1 the ground term, 2 and 3 the
      ! metastable terms (N I, S II: 4S, 2D, 2P; C I: 3P, 1D, 1S).
      integer, parameter :: lp_term(5,3) = reshape( [                    &
         1, 2, 2, 3, 3,   1, 2, 2, 3, 3,   1, 1, 1, 2, 3 ], [5,3] )
      real*8, parameter :: lp_A(n_lp_pair,3) = reshape( [                &
      ! N I
         2.0600d-05, 7.7970d-06, 2.6110d-03, 6.5180d-03, 1.0690d-08,      &
         5.0550d-02, 2.6370d-02, 3.3030d-02, 5.8810d-02, 0.0000d+00,      &
      ! S II
         6.3200d-04, 2.1950d-04, 7.6400d-02, 1.9000d-01, 1.7100d-07,      &
         1.4740d-01, 1.1650d-01, 7.1600d-02, 1.6050d-01, 0.0000d+00,      &
      ! C I
         7.9600d-08, 0.0000d+00, 9.3400d-08, 0.0000d+00, 2.6600d-07,      &
         6.1700d-05, 2.1100d-03, 1.8100d-04, 1.9800d-05, 6.3800d-01 ],    &
         [n_lp_pair,3] )
      real*8, parameter :: lp_ups(n_lp_T,n_lp_pair,3) = reshape( [       &
      ! N I (1,2) ... (4,5)
         2.29951d-02, 3.69951d-02, 5.10006d-02, 7.80091d-02, 1.27994d-01, &
         1.70322d-01, 2.23952d-01, 2.94963d-01, 3.48987d-01, 4.27037d-01, &
         5.19017d-01, 5.62465d-01,                                        &
         3.49843d-02, 5.52179d-02, 7.59829d-02, 1.16988d-01, 1.90955d-01, &
         2.55111d-01, 3.36876d-01, 4.41885d-01, 5.22921d-01, 6.39989d-01, &
         7.77983d-01, 8.43211d-01,                                        &
         5.00070d-03, 7.86974d-03, 1.10007d-02, 1.80010d-02, 2.99978d-02, &
         4.04243d-02, 5.50064d-02, 7.30019d-02, 8.69925d-02, 1.06012d-01, &
         1.27988d-01, 1.38571d-01,                                        &
         1.00015d-02, 1.64279d-02, 2.30014d-02, 3.60019d-02, 5.99952d-02, &
         8.20335d-02, 1.09012d-01, 1.45004d-01, 1.72985d-01, 2.12023d-01, &
         2.54978d-01, 2.75140d-01,                                        &
         1.29980d-02, 2.29814d-02, 3.50040d-02, 6.30158d-02, 1.22003d-01, &
         1.79320d-01, 2.56945d-01, 3.67959d-01, 4.59999d-01, 6.06103d-01, &
         8.01056d-01, 8.97115d-01,                                        &
         5.39956d-02, 6.29104d-02, 6.99998d-02, 8.30034d-02, 1.02997d-01, &
         1.20195d-01, 1.40980d-01, 1.71981d-01, 1.98992d-01, 2.44022d-01, &
         3.11011d-01, 3.45175d-01,                                        &
         5.90030d-02, 7.14794d-02, 8.20120d-02, 1.01021d-01, 1.33018d-01, &
         1.60803d-01, 1.94999d-01, 2.44012d-01, 2.85031d-01, 3.52092d-01, &
         4.46057d-01, 4.93061d-01,                                        &
         4.00018d-02, 4.93599d-02, 5.70078d-02, 7.00146d-02, 9.40116d-02, &
         1.13521d-01, 1.38998d-01, 1.75006d-01, 2.04020d-01, 2.52064d-01, &
         3.20040d-01, 3.54120d-01,                                        &
         1.28999d-01, 1.53031d-01, 1.72012d-01, 2.05025d-01, 2.60014d-01, &
         3.06346d-01, 3.65976d-01, 4.49991d-01, 5.21024d-01, 6.41124d-01, &
         8.16076d-01, 9.04589d-01,                                        &
         2.29986d-02, 2.90875d-02, 3.50018d-02, 4.70066d-02, 7.10011d-02, &
         9.37331d-02, 1.22980d-01, 1.63985d-01, 1.98000d-01, 2.50036d-01, &
         3.17019d-01, 3.49673d-01,                                        &
      ! S II (1,2) ... (4,5)
         3.22098d+00, 3.17714d+00, 3.10933d+00, 2.95191d+00, 2.69986d+00, &
         2.60000d+00, 2.60000d+00, 2.49992d+00, 2.39991d+00, 2.29999d+00, &
         2.03441d+00, 1.36115d+00,                                        &
         4.32471d+00, 4.27251d+00, 4.22388d+00, 4.13708d+00, 3.99993d+00, &
         3.90003d+00, 3.80001d+00, 3.69994d+00, 3.59997d+00, 3.39992d+00, &
         3.19348d+00, 2.48136d+00,                                        &
         7.02548d-01, 7.03083d-01, 7.01857d-01, 6.97342d-01, 6.89997d-01, &
         6.90000d-01, 7.00004d-01, 7.20017d-01, 7.30000d-01, 7.40011d-01, &
         6.95298d-01, 5.03257d-01,                                        &
         1.42262d+00, 1.42359d+00, 1.42138d+00, 1.41324d+00, 1.40000d+00, &
         1.40000d+00, 1.40000d+00, 1.40006d+00, 1.50001d+00, 1.50000d+00, &
         1.46568d+00, 1.14345d+00,                                        &
         7.96414d+00, 7.86020d+00, 7.76147d+00, 7.58265d+00, 7.29996d+00, &
         7.10016d+00, 6.90011d+00, 6.59996d+00, 6.30005d+00, 6.00001d+00, &
         5.62009d+00, 4.40226d+00,                                        &
         1.49995d+00, 1.49995d+00, 1.49995d+00, 1.49997d+00, 1.50000d+00, &
         1.50000d+00, 1.50000d+00, 1.50000d+00, 1.50000d+00, 1.50000d+00, &
         1.46784d+00, 1.15069d+00,                                        &
         2.39992d+00, 2.39991d+00, 2.39992d+00, 2.39995d+00, 2.40000d+00, &
         2.40000d+00, 2.40000d+00, 2.40000d+00, 2.40000d+00, 2.40001d+00, &
         2.34856d+00, 1.84114d+00,                                        &
         1.80036d+00, 1.80037d+00, 1.80034d+00, 1.80021d+00, 1.80000d+00, &
         1.80000d+00, 1.80000d+00, 1.80000d+00, 1.80000d+00, 1.79999d+00, &
         1.55738d+00, 1.01830d+00,                                        &
         4.09986d+00, 4.09985d+00, 4.09986d+00, 4.09992d+00, 4.10000d+00, &
         4.10000d+00, 4.10000d+00, 4.10000d+00, 4.10000d+00, 4.10001d+00, &
         4.01212d+00, 3.14528d+00,                                        &
         1.79045d+00, 1.79004d+00, 1.79098d+00, 1.79441d+00, 1.80000d+00, &
         1.80000d+00, 1.80000d+00, 1.79996d+00, 1.80008d+00, 1.89994d+00, &
         1.66904d+00, 1.14154d+00,                                        &
      ! C I (1,2) ... (4,5)
         3.57089d-02, 7.84514d-02, 1.25218d-01, 2.09305d-01, 3.18990d-01, &
         3.78504d-01, 4.23984d-01, 4.57496d-01, 4.73005d-01, 4.87996d-01, &
         5.00002d-01, 5.05994d-01,                                        &
         3.15178d-02, 6.67698d-02, 1.03423d-01, 1.61662d-01, 2.18009d-01, &
         2.42091d-01, 2.57000d-01, 2.65259d-01, 2.69003d-01, 2.76001d-01, &
         2.89010d-01, 3.09998d-01,                                        &
         2.65021d-02, 4.06038d-02, 5.44323d-02, 7.98910d-02, 1.21991d-01, &
         1.56010d-01, 1.95977d-01, 2.42758d-01, 2.75006d-01, 3.14986d-01, &
         3.49006d-01, 3.59000d-01,                                        &
         2.95040d-03, 4.57854d-03, 6.21070d-03, 9.31876d-03, 1.46995d-02, &
         1.91381d-02, 2.44977d-02, 3.10326d-02, 3.55016d-02, 4.06988d-02, &
         4.43006d-02, 4.34008d-02,                                        &
         1.14989d-01, 2.47914d-01, 3.89174d-01, 6.26129d-01, 8.90907d-01, &
         1.01916d+00, 1.10994d+00, 1.17338d+00, 1.20000d+00, 1.22997d+00, &
         1.28000d+00, 1.32997d+00,                                        &
         7.95051d-02, 1.21760d-01, 1.63293d-01, 2.40236d-01, 3.67970d-01, &
         4.69864d-01, 5.88928d-01, 7.29086d-01, 8.25014d-01, 9.43952d-01, &
         1.05002d+00, 1.08000d+00,                                        &
         8.76005d-03, 1.36037d-02, 1.84629d-02, 2.77225d-02, 4.36942d-02, &
         5.67333d-02, 7.24876d-02, 9.18910d-02, 1.04999d-01, 1.19993d-01, &
         1.30001d-01, 1.28003d-01,                                        &
         1.32021d-01, 2.02826d-01, 2.72495d-01, 4.01264d-01, 6.13992d-01, &
         7.83593d-01, 9.81930d-01, 1.21648d+00, 1.38007d+00, 1.57996d+00, &
         1.75005d+00, 1.80000d+00,                                        &
         1.46012d-02, 2.26858d-02, 3.07972d-02, 4.62465d-02, 7.28943d-02, &
         9.46890d-02, 1.20985d-01, 1.53217d-01, 1.75004d-01, 1.99991d-01, &
         2.18002d-01, 2.14005d-01,                                        &
         7.47032d-02, 9.15621d-02, 1.07632d-01, 1.34251d-01, 1.70994d-01, &
         2.00905d-01, 2.39979d-01, 2.94237d-01, 3.39012d-01, 4.09972d-01, &
         5.07034d-01, 6.55957d-01 ], [n_lp_T,n_lp_pair,3] )


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

      ! Upper clamp on any evaluated rate [cm^3 s^-1], far above any capture
      ! rate, so that a fit evaluated outside its range cannot run away. The
      ! lower clamp is 0: a fit that dips negative below its validity floor
      ! has a negligible true rate there.
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
      ! as upsilon_quartic_logT in Cool_coeff.f90). The authors' own caveat
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

      ! Electron capture by doubly ionized nitrogen from atomic hydrogen,
      !     N2+(2s^2 2p 2Po) + H(1s) -> N+(2s 2p^3 3Do) + H+ ,
      ! rate coefficient in cm^3 s^-1: reaction (4) of Barragan, Errea,
      ! Mendez, Rabadan & Riera (2006), ApJ 636, 544, the calculation that
      ! row E1 carries for O2+ + H, their Table 2 column k4 as published
      ! (units 1e-9 cm^3 s^-1, read from the page image). The product term is
      ! the one their text gives: the rates "agree with ... [those] evaluated
      ! from the cross sections for populating N+ (2s2p3 3Do) of Bienstock
      ! et al. (1986)" (symbols as typeset), and Butler, Heil & Dalgarno
      ! (1980), ApJ 241, 442, Table 2 (read) list the 1084 A decay of that
      ! term for N+2 + H. N+ 2s2p3 3Do (92244.49 cm^-1, the 2J+1 mean of the
      ! NIST levels) decays to the ground term by the allowed N II 1084-1086
      ! A multiplet, so the reverse of this channel starts from a state that
      ! is not populated and is not carried; it is endothermic by
      ! IP(N II) - IP(H I) = 16.00 eV besides. The same paper's reaction (5),
      ! from the metastable N2+(2s 2p^2 4P) 7.09 eV up, is not carried: that
      ! term holds 5e-4 of N III at 1e4 K even at the collisional limit
      ! (DERIVED from the NIST levels).
      ! Interpolated linearly in log10(k) against log10(T) between the
      ! published points; T is held inside the calculated 1e2 - 1e5 K. As
      ! for E1, the authors expect radiative electron capture, not included,
      ! to matter at the lowest temperatures.
      pure double precision function electron_capture_N2p_from_H(T)     &
                                     result(k)
      real*8, intent(in) :: T
      real*8, parameter :: T_n2p(26) = [                                   &
         1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 1.0d3, 1.5d3, 2.0d3,          &
         2.5d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 7.0d3, 8.0d3, 9.0d3,          &
         1.0d4, 2.0d4, 3.0d4, 4.0d4, 5.0d4, 6.0d4, 7.0d4, 8.0d4,          &
         9.0d4, 1.0d5 ]
      real*8, parameter :: k_n2p(26) = [                                   &
         0.46312d-9, 0.59312d-9, 0.65966d-9, 0.70483d-9, 0.73830d-9,      &
         0.82782d-9, 0.86587d-9, 0.88640d-9, 0.89963d-9, 0.90935d-9,      &
         0.92378d-9, 0.93484d-9, 0.94395d-9, 0.95174d-9, 0.95857d-9,      &
         0.96469d-9, 0.97027d-9, 1.01334d-9, 1.06250d-9, 1.13596d-9,      &
         1.22873d-9, 1.33023d-9, 1.43343d-9, 1.53502d-9, 1.63390d-9,      &
         1.73000d-9 ]
      k = log_log_table(T, T_n2p, k_n2p)
      end function electron_capture_N2p_from_H

      ! ---------------------------------------------------------------
      ! THE FIVE-LEVEL STATISTICAL EQUILIBRIUM OF N I, S II AND C I.
      ! The levels are the populated manifold of q_level_energy (data at
      ! lp_E). Processes: spontaneous radiative decay (lp_A) and electron
      ! collisions,
      !   q_ul = C Upsilon/(g_u sqrt T),  q_lu = q_ul (g_u/g_l) exp(-dE/kT),
      ! with C = coll_rate_prefactor (Cooling_Coefficients) and Upsilon the
      ! CHIANTI effective collision strength. Left out, and why: the
      ! levels above the five (allowed decays; their cascade feeding
      ! changes the term departures by less than 0.1% at <= 1e4 K and by
      ! up to 2% at 2e4 K against the full CHIANTI model, MEASURED), proton
      ! collisions (fine structure only), collisions with neutral hydrogen
      ! (not in CHIANTI), and the escape probability of the forbidden lines
      ! (optically thin). R(i,j) is the rate [s^-1] from level i to level j.
      pure subroutine metastable_transition_rates(isp, T, n_e, R)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T, n_e
      real*8,  intent(out) :: R(5,5)
      real*8, parameter :: hc_over_k = hp_erg*c_light/kb_erg
      real*8  :: x, w, ups, qd, qu
      integer :: it, k, u, l
      x = log(min(max(T, lp_T(1)), lp_T(n_lp_T)))
      it = 1
      do while (it .lt. n_lp_T - 1 .and. x .gt. log(lp_T(it+1)))
         it = it + 1
      enddo
      w = (x - log(lp_T(it)))/(log(lp_T(it+1)) - log(lp_T(it)))
      R = 0.0d0
      do k = 1, n_lp_pair
         if (lp_E(lp_hi(k),isp) .ge. lp_E(lp_lo(k),isp)) then
            u = lp_hi(k);  l = lp_lo(k)
         else
            u = lp_lo(k);  l = lp_hi(k)
         endif
         ups = (1.0d0 - w)*lp_ups(it,k,isp) + w*lp_ups(it+1,k,isp)
         qd  = coll_rate_prefactor*ups/(dble(lp_g(u,isp))*sqrt(T))
         qu  = qd*dble(lp_g(u,isp))/dble(lp_g(l,isp))                     &
               *exp(-hc_over_k*(lp_E(u,isp) - lp_E(l,isp))/T)
         R(u,l) = R(u,l) + lp_A(k,isp) + n_e*qd
         R(l,u) = R(l,u) + n_e*qu
      enddo
      end subroutine metastable_transition_rates

      ! The rate matrix M of dp/dt = M p, M(k,j) = R(j,k) for j /= k and
      ! M(k,k) = -sum_j R(k,j); its columns sum to zero (the atom is
      ! conserved), so one row is replaced by the normalization.
      pure subroutine metastable_rate_matrix(isp, T, n_e, M, R)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T, n_e
      real*8,  intent(out) :: M(5,5), R(5,5)
      integer :: k, j
      call metastable_transition_rates(isp, T, n_e, R)
      do k = 1, 5
         do j = 1, 5
            M(k,j) = R(j,k)
         enddo
         M(k,k) = -(sum(R(k,:)) - R(k,k))
      enddo
      end subroutine metastable_rate_matrix

      ! Fractional level populations p(1:5) (sum 1) at T [K] and electron
      ! density n_e [cm^-3]. At n_e -> 0 every level above the ground
      ! decays and p -> (1,0,0,0,0); at large n_e p -> the Boltzmann
      ! distribution over these levels, whose sum is exactly the partition
      ! function of internal_partition_function for the species.
      pure subroutine metastable_level_populations(isp, T, n_e, p)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T, n_e
      real*8,  intent(out) :: p(5)
      real*8  :: M(5,5), R(5,5), b(5)
      call metastable_rate_matrix(isp, T, n_e, M, R)
      M(1,:) = 1.0d0
      b    = 0.0d0
      b(1) = 1.0d0
      call gauss_solve_5(M, b, p)
      end subroutine metastable_level_populations

      ! Departure coefficient of each TERM (1 ground, 2 and 3 metastable,
      ! lp_term) from Boltzmann: beta_tau = p_tau / p_tau(Boltzmann), the
      ! J levels summed within the term. Below about 60 K the Boltzmann
      ! weight of the upper metastable terms underflows (exp(-hc E/kT) <
      ! 1e-308 for N I 2P, S II 2P, C I 1S), and the population the
      ! collisions give such a term underflows with it; the term is then
      ! given zero departure instead of 0/0. Nothing physical is lost: every
      ! reverse that reads a departure carries exp(-dE/kT) < 1e-60 there.
      pure subroutine metastable_term_departures(isp, T, n_e, beta)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T, n_e
      real*8,  intent(out) :: beta(3)
      real*8  :: p(5), wb(5), pt(3), wt(3)
      integer :: c
      call metastable_level_populations(isp, T, n_e, p)
      call metastable_boltzmann_weights(isp, T, wb)
      pt = 0.0d0
      wt = 0.0d0
      do c = 1, 5
         pt(lp_term(c,isp)) = pt(lp_term(c,isp)) + p(c)
         wt(lp_term(c,isp)) = wt(lp_term(c,isp)) + wb(c)
      enddo
      beta = 0.0d0
      do c = 1, 3
         if (wt(c) .gt. 0.0d0) beta(c) = pt(c)*sum(wb)/wt(c)
      enddo
      end subroutine metastable_term_departures

      ! g_c exp(-E_c/kT) of the five levels.
      pure subroutine metastable_boltzmann_weights(isp, T, wb)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T
      real*8,  intent(out) :: wb(5)
      real*8, parameter :: hc_over_k = hp_erg*c_light/kb_erg
      integer :: c
      do c = 1, 5
         wb(c) = dble(lp_g(c,isp))*exp(-hc_over_k*lp_E(c,isp)/T)
      enddo
      end subroutine metastable_boltzmann_weights

      ! THE EXCESS RADIATED ENERGY X_u [eV] OF AN ION PLACED IN LEVEL u:
      ! the energy the ion radiates in forbidden lines beyond what an ion
      ! of the statistical-equilibrium distribution p radiates, over the
      ! time it takes to forget where it started,
      !   X_u = sum_k rad_k y_k,   y = integral_0^inf (P(t|u) - p) dt,
      !   M y = -(e_u - p),  sum_k y_k = 0,
      ! with rad_k = sum_l A(k->l) (E_k - E_l) the power radiated from
      ! level k. This is the excitation a product formed in level u does
      ! NOT return to the gas: at n_e -> 0 X_u = E_u (the level radiates
      ! down), at large n_e X_u -> 0 (collisions return it), which are the
      ! two limits the populated-manifold heat of charge_exchange_heating
      ! lies between. (The first-passage radiation to the ground level is
      ! not this quantity: it also counts radiation that an ion of the
      ! equilibrium distribution would emit in the same time.)
      ! p (optional) returns the equilibrium populations the excess is
      ! measured from (metastable_level_populations), formed here anyway.
      pure subroutine metastable_excess_radiation(isp, T, n_e, X, p_out)
      integer, intent(in)  :: isp
      real*8,  intent(in)  :: T, n_e
      real*8,  intent(out) :: X(5)
      real*8,  intent(out), optional :: p_out(5)
      real*8, parameter :: cm_to_eV = hp_erg*c_light*erg2eV
      real*8  :: M(5,5), R(5,5), Mn(5,5), p(5), b(5), y(5), rad(5)
      integer :: u, k, l
      call metastable_rate_matrix(isp, T, n_e, M, R)
      Mn = M
      Mn(1,:) = 1.0d0
      b    = 0.0d0
      b(1) = 1.0d0
      call gauss_solve_5(Mn, b, p)
      if (present(p_out)) p_out = p
      ! Radiated power of each level [cm^-1 s^-1]: the radiative part of
      ! R, which is lp_A on the downward entries.
      rad = 0.0d0
      do k = 1, n_lp_pair
         if (lp_E(lp_hi(k),isp) .ge. lp_E(lp_lo(k),isp)) then
            u = lp_hi(k);  l = lp_lo(k)
         else
            u = lp_lo(k);  l = lp_hi(k)
         endif
         rad(u) = rad(u) + lp_A(k,isp)*(lp_E(u,isp) - lp_E(l,isp))
      enddo
      do u = 1, 5
         b = p
         b(u) = b(u) - 1.0d0
         b(1) = 0.0d0
         call gauss_solve_5(Mn, b, y)
         X(u) = sum(rad*y)*cm_to_eV
      enddo
      end subroutine metastable_excess_radiation

      ! x = A^-1 b for a 5 x 5 system, Gaussian elimination with partial
      ! pivoting (the level systems span twenty decades between their
      ! radiative and collisional entries).
      pure subroutine gauss_solve_5(A_in, b_in, x)
      real*8, intent(in)  :: A_in(5,5), b_in(5)
      real*8, intent(out) :: x(5)
      real*8  :: A(5,5), b(5), f, row(5), t
      integer :: i, j, k, ip
      A = A_in
      b = b_in
      do k = 1, 4
         ip = k
         do i = k+1, 5
            if (abs(A(i,k)) .gt. abs(A(ip,k))) ip = i
         enddo
         if (ip .ne. k) then
            row = A(k,:);  A(k,:) = A(ip,:);  A(ip,:) = row
            t = b(k);      b(k) = b(ip);      b(ip) = t
         endif
         do i = k+1, 5
            f = A(i,k)/A(k,k)
            A(i,k:5) = A(i,k:5) - f*A(k,k:5)
            b(i) = b(i) - f*b(k)
         enddo
      enddo
      do i = 5, 1, -1
         t = b(i)
         do j = i+1, 5
            t = t - A(i,j)*x(j)
         enddo
         x(i) = t/A(i,i)
      enddo
      end subroutine gauss_solve_5

      ! The level-population model of the species at (element code,
      ! stage), 0 when there is none.
      pure integer function metastable_model_of(el, stg) result(isp)
      integer, intent(in) :: el, stg
      isp = 0
      if (el .eq. iel_N .and. stg .eq. 0) isp = lp_NI
      if (el .eq. iel_S .and. stg .eq. 1) isp = lp_SII
      if (el .eq. iel_C .and. stg .eq. 0) isp = lp_CI
      end function metastable_model_of

      ! F of the note at q_level_energy for the forward row f: the factor
      ! by which the populations of f's product levels at (T, n_e) move the
      ! detailed-balance reverse away from its thermal value. 1 for every
      ! forward whose source does not name its product levels (thermal
      ! branching, for which F = 1 at every density).
      double precision function product_population_factor(f, T, n_e)    &
                                result(pf)
      integer, intent(in) :: f
      real*8,  intent(in) :: T, n_e
      real*8  :: beta(3), beta2(3), bs(3)
      select case (f)
      case (18)
         ! A18: N(4S) alone.
         call metastable_term_departures(lp_NI, T, n_e, beta)
         pf = beta(1)
      case (19)
         ! A19: S+(4S, 2D, 2P) with the branching of Zhao et al. (2005).
         call sulfur_product_term_branching(T, bs)
         call metastable_term_departures(lp_SII, T, n_e, beta)
         pf = sum(bs*beta)
      case (56)
         ! D25: C(3P) + S+(4S).
         call metastable_term_departures(lp_CI,  T, n_e, beta)
         call metastable_term_departures(lp_SII, T, n_e, beta2)
         pf = beta(1)*beta2(1)
      case default
         pf = 1.0d0
      end select
      end function product_population_factor

      ! The fraction of S + H+ -> S+ + H that ends in each S+ term (4S, 2D,
      ! 2P), for S I in the thermal distribution over 3P and 1D that
      ! sulfur_charge_transfer_ionization weights the rate with: Zhao et
      ! al. (2005), Phys. Rev. A 71, 062713, Table III (read from the page
      ! image), the exit-channel columns, interpolated in log-log and held
      ! at the ends of 20 K - 1e6 K. S(1D) + H+ has no S+(4S) exit column
      ! (the quartet is spin-forbidden from the singlet). The fractions are
      ! formed from the channel columns; the "Total" columns, which the
      ! rate itself uses, equal their sums to the rounding of the table.
      ! At 1e4 K: 5.1e-5, 0.068, 0.932 (DERIVED).
      subroutine sulfur_product_term_branching(T, bs)
      real*8, intent(in)  :: T
      real*8, intent(out) :: bs(3)
      real*8, parameter :: hc_over_k = hp_erg*c_light/kb_erg
      integer, parameter :: n_z = 25
      real*8, parameter :: T_z(n_z) = [ 2.0d1, 4.0d1, 6.0d1, 8.0d1,     &
         1.0d2, 2.0d2, 4.0d2, 6.0d2, 8.0d2, 1.0d3, 2.0d3, 4.0d3, 6.0d3,  &
         8.0d3, 1.0d4, 2.0d4, 4.0d4, 6.0d4, 8.0d4, 1.0d5, 2.0d5, 4.0d5,  &
         6.0d5, 8.0d5, 1.0d6 ]
      ! S(3P) + H+ -> S+(4S), S+(2D), S+(2P) + H
      real*8, parameter :: k_3p_4s(n_z) = [ 3.10d-16, 3.90d-16,         &
         4.09d-16, 4.21d-16, 4.12d-16, 3.71d-16, 3.28d-16, 3.07d-16,     &
         2.96d-16, 2.93d-16, 3.11d-16, 4.05d-16, 5.35d-16, 6.97d-16,     &
         8.90d-16, 8.10d-15, 7.94d-13, 4.34d-12, 1.03d-11, 1.77d-11,     &
         6.30d-11, 1.59d-10, 2.48d-10, 3.40d-10, 4.39d-10 ]
      real*8, parameter :: k_3p_2d(n_z) = [ 1.87d-15, 4.33d-15,         &
         7.32d-15, 1.02d-14, 1.28d-14, 2.21d-14, 3.03d-14, 3.24d-14,     &
         3.33d-14, 3.43d-14, 4.48d-14, 1.23d-13, 3.22d-13, 6.64d-13,     &
         1.19d-12, 8.65d-12, 5.95d-11, 1.52d-10, 2.69d-10, 4.02d-10,     &
         1.13d-9, 2.53d-9, 3.76d-9, 4.83d-9, 5.81d-9 ]
      real*8, parameter :: k_3p_2p(n_z) = [ 3.16d-13, 8.51d-13,         &
         1.20d-12, 1.45d-12, 1.66d-12, 2.44d-12, 3.28d-12, 3.65d-12,     &
         3.85d-12, 3.98d-12, 4.39d-12, 5.34d-12, 7.52d-12, 1.11d-11,     &
         1.63d-11, 6.92d-11, 3.13d-10, 6.97d-10, 1.17d-9, 1.68d-9,       &
         4.40d-9, 9.34d-9, 1.35d-8, 1.70d-8, 2.01d-8 ]
      ! S(1D) + H+ -> S+(2D), S+(2P) + H
      real*8, parameter :: k_1d_2d(n_z) = [ 1.17d-16, 1.40d-16,         &
         1.68d-16, 1.99d-16, 2.14d-16, 2.69d-16, 3.05d-16, 3.16d-16,     &
         3.27d-16, 3.42d-16, 5.38d-16, 1.76d-15, 4.25d-15, 8.78d-15,     &
         1.87d-14, 4.45d-13, 7.64d-12, 2.68d-11, 5.51d-11, 8.83d-11,     &
         2.65d-10, 6.16d-10, 1.04d-9, 1.56d-9, 2.19d-9 ]
      real*8, parameter :: k_1d_2p(n_z) = [ 1.57d-18, 3.21d-18,         &
         4.46d-18, 5.84d-18, 7.28d-18, 1.34d-17, 2.16d-17, 3.15d-17,     &
         4.97d-17, 8.39d-17, 7.25d-16, 6.62d-15, 1.84d-14, 3.56d-14,     &
         6.83d-14, 1.64d-12, 1.63d-11, 3.54d-11, 5.25d-11, 6.76d-11,     &
         1.52d-10, 5.35d-10, 1.24d-9, 2.17d-9, 3.21d-9 ]
      real*8 :: w3p, w1d
      ! S I 3P2, 3P1, 3P0 and 1D2 (NIST, the levels of q_level_energy),
      ! the weights of sulfur_charge_transfer_ionization.
      w3p = 5.0d0 + 3.0d0*exp(-hc_over_k*396.05648d0/T)                   &
          + 1.0d0*exp(-hc_over_k*573.59573d0/T)
      w1d = 5.0d0*exp(-hc_over_k*9238.6090d0/T)
      bs(1) = w3p*log_log_table(T, T_z, k_3p_4s)
      bs(2) = w3p*log_log_table(T, T_z, k_3p_2d)                          &
            + w1d*log_log_table(T, T_z, k_1d_2d)
      bs(3) = w3p*log_log_table(T, T_z, k_3p_2p)                          &
            + w1d*log_log_table(T, T_z, k_1d_2p)
      bs = bs/sum(bs)
      end subroutine sulfur_product_term_branching

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
      call cx_check_reverse_pairs()
      end subroutine cx_init

      ! The reverse table is data, so check once that it describes what
      ! detailed_balance_ratio assumes: a reverse row is its forward read
      ! backwards (donor and acceptor exchanged, each one stage apart), the
      ! forward is exothermic, so the exponent of the ratio is negative and
      ! cannot overflow, and the forward is itself a published row.
      subroutine cx_check_reverse_pairs()
      integer :: r, f
      logical :: ok
      do r = 1, n_cxreac
         f = cx_reverse_of(r)
         if (f .eq. 0) cycle
         ok = (cx_reverse_of(f) .eq. 0)                                   &
              .and. cx_don_el(r) .eq. cx_acc_el(f)                         &
              .and. cx_acc_el(r) .eq. cx_don_el(f)                         &
              .and. cx_don_stg(r) .eq. cx_acc_stg(f) - 1                   &
              .and. cx_acc_stg(r) .eq. cx_don_stg(f) + 1
         if (ok) ok = (cx_ionization_energy(cx_acc_el(f),                 &
                                            cx_acc_stg(f)-1)              &
                       .gt. cx_ionization_energy(cx_don_el(f),            &
                                                 cx_don_stg(f)))
         if (.not. ok) then
            write(*,*) '(charge_exchange) ERROR: row', r, ' is declared', &
                       ' the detailed-balance reverse of row', f, ','
            write(*,*) '  but the two are not one exothermic reaction', &
                       ' read in two directions. Aborting.'
            error stop 1
         endif
      enddo
      end subroutine cx_check_reverse_pairs

      ! Membership test for the generic active set, written once so cx_init's
      ! two passes cannot drift apart.
      !   Group B (He<->H) is applied by he_h_cx_*, so it is never here.
      !   Group E is not in Table 4, so cx_full does not reach it; each row
      !   is gated on its own scale factor, positive by default and zeroed
      !   by "cx_O2p_H 0" or "cx_N2p_H 0".
      !   A2, A6, A10, A22, C2, C4, C6, D18, D28 and D30 are not carried
      !   (the notes at cx_A2_Mgp_H).
      logical function cx_row_is_active(r)
      integer, intent(in) :: r
      if (r .eq. cx_B1_He0_Hp .or. r .eq. cx_B2_Hep_H0) then
         cx_row_is_active = .false.
      else if (r .eq. cx_A2_Mgp_H .or. r .eq. cx_A10_Sip_H             &
               .or. r .eq. cx_A6_Fep_H .or. r .eq. cx_A22_Nap_H        &
               .or. r .eq. cx_D18_Fep_O .or. r .eq. cx_D28_Fep_N       &
               .or. r .eq. cx_C2_Sip_He .or. r .eq. cx_D30_Nap_Ca      &
               .or. r .eq. cx_C4_Cp_He .or. r .eq. cx_C6_Op_He) then
         cx_row_is_active = .false.
      else if (r .eq. cx_E1_O2p_H0) then
         cx_row_is_active = (cx_o2p_h_scale .gt. 0.0d0)
      else if (r .eq. cx_E2_N2p_H0) then
         cx_row_is_active = (cx_n2p_h_scale .gt. 0.0d0)
      else
         cx_row_is_active = (cx_default(r) .or. cx_full)
      endif
      end function cx_row_is_active

      ! Evaluate the active-reaction rate coefficients at temperature T [K]
      ! and electron density n_e [cm^-3] of one cell (n_e sets the level
      ! populations the derived reverses depend on, the note at
      ! q_level_energy). Called per cell before the hybrd1 solve, like
      ! set_metal_coeffs, with the electron density the cell's other
      ! coefficients are formed with.
      subroutine cx_set_cell(T, n_e)
      real*8, intent(in) :: T, n_e
      integer :: i
      real*8  :: kc
      ! cx_kc is threadprivate: cx_init allocated only the master thread's copy,
      ! so each worker thread allocates its own on first use here (cx_nact is the
      ! shared, setup-once active-reaction count).
      if (.not. allocated(cx_kc)) allocate(cx_kc(cx_nact))
      do i = 1, cx_nact
         kc = cx_rate(cx_act(i), T, n_e)
         cx_kc(i) = min(max(kc, 0.0d0), cx_kc_max)
      enddo
      ! The rates now belong to this cell; the assembly routines check
      ! that the cell they are handed is this one by its temperature.
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
                                he_row_sign, T_cell, gross)
      integer, intent(in)    :: N_eq
      real*8,  intent(inout) :: fvec(N_eq)
      ! The gross rate of each row [cm^-3 s^-1], to which every reaction
      ! adds its magnitude on both rows it reaches (the normalization of
      ! the molecular systems, System_HeH_mol).
      real*8,  optional, intent(inout) :: gross(N_eq)
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
         if (present(gross)) then
            gross(idon) = gross(idon) + abs(rrate(i))
            gross(iacc) = gross(iacc) + abs(rrate(i))
         endif
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
      ! table reaches. The donors of Table 4 and of group E are neutral
      ! or, in rows A3 (Mg+), A7 (Fe+) and A11 (Si+), singly ionized, and
      ! every acceptor is singly or doubly ionized, so a reaction writes
      ! stages 0 to 2 only; the extra row exists so that a future table
      ! row with a doubly ionized donor, which would write stage 3, is
      ! caught here by name instead of writing past the metal arrays.
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

      ! Rate coefficient [cm^3 s^-1] of row id at temperature T [K] and
      ! electron density n_e [cm^-3]: the source's own expression for a
      ! row carried as published, and for a row listed in cx_reverse_of the
      ! detailed-balance reverse of its forward row, formed from that
      ! forward at the same temperature and corrected for the populations
      ! of the forward's product levels (product_population_factor).
      real*8 function cx_rate(id, T, n_e)
      integer, intent(in) :: id
      real*8,  intent(in) :: T, n_e
      integer :: f
      f = cx_reverse_of(id)
      if (f .gt. 0) then
         cx_rate = published_rate_coefficient(f, T)                        &
                   *detailed_balance_ratio(f, T)                           &
                   *product_population_factor(f, T, n_e)
      else
         cx_rate = published_rate_coefficient(id, T)
      endif
      end function cx_rate

      ! k_reverse/k_forward of the forward row f, X^s + Y^(t+1) ->
      ! X^(s+1) + Y^t (donor X^s, acceptor Y^(t+1)), in thermal
      ! equilibrium at T:
      !   [Q(X^s) Q(Y^(t+1))] / [Q(X^(s+1)) Q(Y^t)] exp(-dE/kT),
      !   dE = IP(Y^t) - IP(X^s) > 0 (cx_check_reverse_pairs),
      ! so the exponent is negative and the ratio cannot overflow; at low T
      ! it underflows to zero, which is the value. The reduced-mass factor
      ! (mu_f/mu_r)^(3/2) differs from 1 by electron masses only, a few
      ! parts in 1e4, and is left out.
      real*8 function detailed_balance_ratio(f, T)
      integer, intent(in) :: f
      real*8,  intent(in) :: T
      integer :: xe, xs, ye, yt
      real*8  :: dE
      xe = cx_don_el(f);  xs = cx_don_stg(f)
      ye = cx_acc_el(f);  yt = cx_acc_stg(f) - 1
      dE = cx_ionization_energy(ye, yt) - cx_ionization_energy(xe, xs)
      detailed_balance_ratio =                                            &
           internal_partition_function(xe, xs,   T)                       &
         * internal_partition_function(ye, yt+1, T)                       &
         / ( internal_partition_function(xe, xs+1, T)                     &
           * internal_partition_function(ye, yt,   T) )                   &
         * exp(-dE/(kb_eV*T))
      end function detailed_balance_ratio

      ! Internal partition function of element code el (canonical metals
      ! 1..n_melem, then 11 = H, 12 = He) at ionization stage stg and
      ! temperature T [K]: the sum of (2J+1) exp(-E/kT) over the levels of
      ! q_level_energy, the populated manifold of the species (the note at
      ! that table). hc/k = hp_erg c_light / kb_erg [cm K].
      double precision function internal_partition_function(el, stg, T)  &
                                result(q)
      integer, intent(in) :: el, stg
      real*8,  intent(in) :: T
      real*8, parameter :: hc_over_k = hp_erg*c_light/kb_erg
      integer :: i
      if (el .lt. 1 .or. el .gt. 12 .or. stg .lt. 0 .or. stg .gt. 2) then
         q = -1.0d0
      else if (q_first(el,stg) .eq. 0) then
         q = -1.0d0
      else
         q = 0.0d0
         do i = q_first(el,stg), q_last(el,stg)
            q = q + dble(q_level_weight(i))                                &
                    *exp(-hc_over_k*q_level_energy(i)/T)
         enddo
      endif
      if (q .le. 0.0d0) then
         write(*,*) '(charge_exchange) ERROR: no partition function is', &
                    ' tabulated for element code', el, ' stage', stg
         write(*,*) '  Aborting.'
         error stop 1
      endif
      end function internal_partition_function

      ! C + H+ -> C+ + H (+ photon): the SUM of the non-radiative and the
      ! radiative channel, reactions (2) + (3) of Stancil, Havener, Krstic,
      ! Schultz, Kimura, Gu, Hirsch, Buenker & Zygelman (1998), ApJ 502,
      ! 1006, read from the published pages. The values are their Table 1,
      ! column "Reactions (2) + (3)", 10 - 1e7 K, interpolated linearly in
      ! log10(k) against log10(T); outside that range the endpoint is held.
      ! The table, not their eq. (4) fit, is carried: with the Table 1
      ! parameters (a2 = 1.28e-15, b2 = 3.08, c2 = 2.19e6 K) the printed
      ! form exp(-c/T) gives 3.2e-15 at 2e4 K and 3.9e-15 at 5e4 K against
      ! the tabulated 1.16e-14 and 1.67e-13, so the fit as printed does not
      ! represent the calculation above 1e4 K; its first term alone is what
      ! Huang et al. (2023) Table 4 and Glover & Jappsen (2007) R39 print.
      ! Below 1e4 K the radiative channel dominates ("for lower temperatures
      ! the reaction of H+ with C proceeds primarily by radiative charge
      ! transfer"), so this row has no collisional reverse there; its
      ! non-radiative reverse is row A16.
      pure double precision function carbon_charge_transfer_ionization(T) &
                                     result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_c = 19
      real*8, parameter :: T_c(n_c) = [                                  &
         1.0d1, 2.0d1, 5.0d1, 1.0d2, 2.0d2, 5.0d2, 1.0d3, 2.0d3, 5.0d3,   &
         1.0d4, 2.0d4, 5.0d4, 1.0d5, 2.0d5, 5.0d5, 1.0d6, 2.0d6, 5.0d6,   &
         1.0d7 ]
      real*8, parameter :: k_c(n_c) = [                                  &
         1.27d-15, 5.03d-16, 9.75d-16, 1.64d-15, 1.86d-15, 1.66d-15,      &
         1.62d-15, 1.81d-15, 2.51d-15, 3.90d-15, 1.16d-14, 1.67d-13,      &
         1.50d-12, 1.27d-11, 2.13d-10, 1.29d-9, 5.14d-9, 1.92d-8,         &
         3.95d-8 ]
      real*8  :: x, w
      integer :: i
      x = log10(min(max(T, T_c(1)), T_c(n_c)))
      i = 1
      do while (i .lt. n_c - 1 .and. x .gt. log10(T_c(i+1)))
         i = i + 1
      enddo
      w = (x - log10(T_c(i)))/(log10(T_c(i+1)) - log10(T_c(i)))
      k = 10.0d0**(log10(k_c(i)) + w*(log10(k_c(i+1)) - log10(k_c(i))))
      end function carbon_charge_transfer_ionization

      ! Linear interpolation of log10(k) in log10(T) through a table,
      ! held at the end points outside it.
      pure double precision function log_log_table(T, Tt, kt) result(k)
      real*8, intent(in) :: T, Tt(:), kt(:)
      real*8  :: x, w
      integer :: i, n
      n = size(Tt)
      x = log10(min(max(T, Tt(1)), Tt(n)))
      i = 1
      do while (i .lt. n - 1 .and. x .gt. log10(Tt(i+1)))
         i = i + 1
      enddo
      w = (x - log10(Tt(i)))/(log10(Tt(i+1)) - log10(Tt(i)))
      k = 10.0d0**(log10(kt(i)) + w*(log10(kt(i+1)) - log10(kt(i))))
      end function log_log_table

      ! Charge transfer of H+ (ion = 1), O+ (2) and N+ (3) with neutral iron,
      ! X+ + Fe -> X + Fe+, rate coefficient [cm^3 s^-1]: Rutherford &
      ! Vroom (1972), J. Chem. Phys. 57, 3091, Table I (read), "Calculated
      ! cross sections and rate coefficients" at 300, 600 and 1200 K:
      !   H+ + Fe 73.5, 53.4, 39.7; O+ + Fe 29.4, 22.2, 17.1;
      !   N+ + Fe 14.7, 11.4, 9.46   (units 1e-10 cm^3 s^-1),
      ! interpolated in log-log and held at the end points. The table is an
      ! extrapolation of cross sections measured from 2 to 500 eV (H+ from
      ! 8 eV) that "operates on the assumption that the relative abundance
      ! of the various products emerging from the reaction stays constant
      ! over the total energy range from thermal to 500 eV"; the
      ! maximum error of the measurements is "+-50% at high impact energies
      ! increasing to as much as a factor of 2.2 at the lowest energies"
      ! (their plus-minus sign written as +-). Huang et al. (2023) Table 4
      ! print the 1200 K values (4.0e-9, 1.71e-9, 9.46e-10) as constants,
      ! which also sets them at 300 - 1200 K, where the paper gives up to
      ! 1.85 times more.
      ! ABOVE 1200 K THE 1200 K VALUE IS HELD, and this is why. The measured
      ! H+ + Fe cross section (Fig. 1, read from the page image) runs from
      ! 8 eV (laboratory; 7.86 eV in the center of mass) to 500 eV and
      ! is nearly flat, 9.8e-16 cm^2 at 8 eV falling to 7.8e-16 at 500 eV
      ! (log-log slope -0.045); the paper says the cross sections "show an
      ! increase in magnitude with decreasing ion energy". A Maxwellian at
      ! 1e4 K draws all of its rate from below 7.9 eV, so no temperature of
      ! these winds is covered by the measurement itself: every thermal
      ! rate is an extrapolation below the lowest measured energy. DERIVED
      ! Maxwellian averages at 3e3, 5e3, 1e4, 2e4 K:
      !   (a) the measured cross section held at its 8 eV value below it:
      !       0.79, 1.01, 1.43, 2.02 x 1e-9 (its measured slope continued
      !       down instead: 0.88, 1.12, 1.53, 2.10);
      !   (b) the measured cross section joined, log-log, to the paper's own
      !       low-energy extrapolation (the Table I cross sections, which
      !       times the rms speed give the Table I rates, placed at
      !       E = 3kT/2): 4.21, 4.06, 3.96, 3.91 x 1e-9.
      ! (b) is the source's statement: the Table I extrapolation was
      ! accepted by the authors only where it fits the lowest measured
      ! points to +-10%, and the join between 0.16 and 7.9 eV has slope
      ! -0.51, the E^-1/2 of capture by the ion-induced dipole attraction,
      ! which is what a cross section rising toward low energy, as the
      ! measured one does, follows there. (a) is a lower bound that assumes
      ! the rise stops at the lowest measured point. The held 3.97e-9 is
      ! (b) to 6% over 3e3 - 2e4 K and continuous with Table I at 1200 K,
      ! so it is carried; the band (a) - (b), 1.4e-9 - 4.0e-9 at 1e4 K, is
      ! within the paper's "factor of 2.2 at the lowest energies" and is the
      ! uncertainty of this rate in a wind. As a secondary cross-check only,
      ! Pequignot & Aldrovandi (1986), A&A 161, 169, p. 172 (read), quote
      ! "beta(H+, Fe) = 2.7 10^-9 cm^3 s^-1 inferred from Rutherford and
      ! Vroom (1972) around 10^4 K" (symbols as typeset), inside the band.
      ! VALIDITY: 300 - 1200 K as tabulated; above 1200 K the held value
      ! is (b) within 6% up to 2e4 K, with the band above. For O+ and N+
      ! the same construction (b) gives 1.6e-9, 1.3e-9, 1.1e-9 and 1.05e-9,
      ! 1.12e-9, 1.23e-9 at 3e3, 1e4, 2e4 K against the held 1.71e-9 and
      ! 9.46e-10.
      pure double precision function iron_charge_transfer_rate(ion, T)   &
                                     result(k)
      integer, intent(in) :: ion
      real*8,  intent(in) :: T
      real*8, parameter :: T_rv(3) = [ 3.0d2, 6.0d2, 1.2d3 ]
      real*8, parameter :: k_rv(3,3) = reshape( [                        &
         73.5d-10, 53.4d-10, 39.7d-10,                                   &
         29.4d-10, 22.2d-10, 17.1d-10,                                   &
         14.7d-10, 11.4d-10, 9.46d-10 ], [3,3] )
      k = log_log_table(T, T_rv, k_rv(:,ion))
      end function iron_charge_transfer_rate

      ! Na(3s) + H+ -> Na+ + H(n=2) - 1.682 eV, the NON-RADIATIVE channel
      ! of Dutta, Nordlander, Kimura & Dalgarno (2001), Phys. Rev. A 63,
      ! 022709 (read): "Cross sections for nonradiative charge transfer in
      ! H+ + Na(3s) collisions at energies less than or equal to 40 eV",
      ! a six-state close-coupling calculation, the process written
      ! "Na(3s) + H+ -> Na+ + H (n=2) - 1.682 eV" (symbols as typeset).
      ! Rate coefficient [cm^3 s^-1] = the Maxwellian average of their total
      ! cross section (Fig. 5, "Present calculation", open circles read from
      ! the page image: 0.110, 0.124, 0.237, 0.340, 0.511, 0.723, 0.976
      ! x 1e-16 cm^2 at 2.26, 2.83, 5.65, 8.51, 14.2, 22.7, 34.0 eV),
      ! log-log interpolated, held at 0.110e-16 cm^2 from their lowest
      ! point down to the 1.682 eV threshold and zero below it (DERIVED,
      ! tabulated below; with zero below 2.26 eV instead the rate is 0.69,
      ! 0.36, 0.15 of these at 1e4, 5e3, 3e3 K).
      ! The paper's own rate coefficients (its Fig. 7) are not carried:
      ! below about 5000 K they exceed what the paper's cross sections
      ! give (5e-13 at 1e3 K against at most 3.5e-19 from Fig. 5, for a
      ! channel that is endothermic by 1.68 eV) and break the threshold law
      ! (module header), which a cross section that vanishes below
      ! threshold cannot do; at 2e4 K the two agree (2.8e-11 against
      ! 2.87e-11). VALIDITY: cross sections to 40 eV, rates quoted "at
      ! temperatures below 20 000 K"; the table runs to 5e4 K and is held
      ! above; below 500 K the rate is below 1e-26 and is set to zero.
      pure double precision function sodium_nonradiative_ionization(T)   &
                                     result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_d = 15
      real*8, parameter :: T_d(n_d) = [ 5.0d2, 7.0d2, 1.0d3, 1.5d3,     &
         2.0d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 8.0d3, 1.0d4, 1.5d4, 2.0d4,  &
         3.0d4, 5.0d4 ]
      real*8, parameter :: k_d(n_d) = [ 1.6235d-27, 9.6747d-23,         &
         3.5282d-19, 1.9760d-16, 4.5396d-15, 1.0113d-13, 4.7200d-13,     &
         1.1917d-12, 2.2245d-12, 4.9678d-12, 8.2860d-12, 1.7913d-11,     &
         2.8686d-11, 5.2491d-11, 1.0601d-10 ]
      if (T .lt. T_d(1)) then
         k = 0.0d0
      else
         k = log_log_table(T, T_d, k_d)
      endif
      end function sodium_nonradiative_ionization

      ! Na(3s) + H+ -> Na+ + H(1s) + photon and K(4s) + H+ -> K+ + H(1s)
      ! + photon, the RADIATIVE channels of Watanabe, Dutta, Nordlander,
      ! Kimura & Dalgarno (2002), Phys. Rev. A 66, 044701 (read): an
      ! optical-potential calculation, "the radiative charge-transfer rate
      ! coefficients were obtained for temperatures below 10 000 K". Their
      ! Fig. 5 (circles Na, triangles K; linear axis in 1e-12 cm^3 s^-1),
      ! the computed points read from the page image, 1e2 (Na) or 4e2 (K)
      ! to 1e4 K; log-log interpolated and held at the end points (the
      ! rates fall as about T^-0.3 near 1e4 K, so holding overstates them a
      ! little above it; they are 1e-2 of the non-radiative Na rate there).
      pure double precision function sodium_radiative_charge_transfer(T) &
                                     result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_w = 15
      real*8, parameter :: T_w(n_w) = [ 1.0d2, 2.0d2, 4.0d2, 6.0d2,     &
         8.0d2, 1.0d3, 2.0d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 7.0d3, 8.0d3,  &
         9.0d3, 1.0d4 ]
      real*8, parameter :: k_w(n_w) = [ 0.638d-12, 0.540d-12,           &
         0.419d-12, 0.354d-12, 0.311d-12, 0.280d-12, 0.198d-12,          &
         0.161d-12, 0.139d-12, 0.123d-12, 0.112d-12, 0.102d-12,          &
         0.097d-12, 0.091d-12, 0.086d-12 ]
      k = log_log_table(T, T_w, k_w)
      end function sodium_radiative_charge_transfer

      pure double precision function potassium_radiative_charge_transfer(T) &
                                     result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_w = 13
      real*8, parameter :: T_w(n_w) = [ 4.0d2, 6.0d2, 8.0d2, 1.0d3,     &
         2.0d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 7.0d3, 8.0d3, 9.0d3, 1.0d4 ]
      real*8, parameter :: k_w(n_w) = [ 0.674d-12, 0.546d-12,           &
         0.465d-12, 0.411d-12, 0.278d-12, 0.222d-12, 0.190d-12,          &
         0.167d-12, 0.150d-12, 0.136d-12, 0.127d-12, 0.120d-12,          &
         0.114d-12 ]
      k = log_log_table(T, T_w, k_w)
      end function potassium_radiative_charge_transfer

      ! S + H+ -> S+ + H, rate coefficient [cm^3 s^-1] for the S I of the
      ! populated manifold: Zhao, Stancil, Gu, Liebermann, Funke, Buenker,
      ! Zygelman, Kimura & Dalgarno (2005), Phys. Rev. A 71, 062713 (read),
      ! a quantum-mechanical molecular-orbital close-coupling calculation,
      ! Table III, the "Total" columns (summed over S+ 4S, 2D and 2P) for
      ! S(3P) + H+ and S(1D) + H+, 20 K - 1e6 K ("temperatures between
      ! 10 K and 2.0x10^6 K"), each interpolated in log-log and held at the
      ! ends. The two are weighted with the Boltzmann populations of 3P
      ! and 1D within Q(S I); the 1S term (2.75 eV, 0.4% of S I at 1e4 K)
      ! was not calculated and contributes nothing. Every product term lies
      ! in the populated manifold of S+, so the reverse A20 is the
      ! detailed-balance one. Huang et al. (2023) print an exp[sum c_i
      ! (ln T)^i] fit (1.83e-11 at 1e4 K, against 1.75e-11 for S(3P) in
      ! the table); Kingdon & Ferland's (1996) 1e-14 assumed a radiative
      ! process, which Zhao et al. "do not support". Their dominant channel
      ! is "S(3P) + H+ -> S+(2P0) + H" (as typeset).
      double precision function sulfur_charge_transfer_ionization(T)     &
                                result(k)
      real*8, intent(in) :: T
      real*8, parameter :: hc_over_k = hp_erg*c_light/kb_erg
      integer, parameter :: n_z = 25
      real*8, parameter :: T_z(n_z) = [ 2.0d1, 4.0d1, 6.0d1, 8.0d1,     &
         1.0d2, 2.0d2, 4.0d2, 6.0d2, 8.0d2, 1.0d3, 2.0d3, 4.0d3, 6.0d3,  &
         8.0d3, 1.0d4, 2.0d4, 4.0d4, 6.0d4, 8.0d4, 1.0d5, 2.0d5, 4.0d5,  &
         6.0d5, 8.0d5, 1.0d6 ]
      real*8, parameter :: k_3p(n_z) = [ 3.18d-13, 8.56d-13, 1.21d-12,  &
         1.46d-12, 1.67d-12, 2.47d-12, 3.31d-12, 3.68d-12, 3.88d-12,     &
         4.01d-12, 4.44d-12, 5.46d-12, 7.84d-12, 1.18d-11, 1.75d-11,     &
         7.79d-11, 3.74d-10, 8.53d-10, 1.45d-9, 2.10d-9, 5.59d-9,        &
         1.20d-8, 1.75d-8, 2.22d-8, 2.64d-8 ]
      real*8, parameter :: k_1d(n_z) = [ 1.19d-16, 1.43d-16, 1.72d-16,  &
         2.05d-16, 2.21d-16, 2.83d-16, 3.26d-16, 3.48d-16, 3.77d-16,     &
         4.26d-16, 1.26d-15, 8.38d-15, 2.27d-14, 4.44d-14, 8.70d-14,     &
         2.09d-12, 2.39d-11, 6.23d-11, 1.08d-10, 1.56d-10, 4.17d-10,     &
         1.15d-9, 2.28d-9, 3.73d-9, 5.40d-9 ]
      real*8 :: w3p, w1d
      ! S I 3P2, 3P1, 3P0 and 1D2 (NIST, the levels of q_level_energy)
      w3p = 5.0d0 + 3.0d0*exp(-hc_over_k*396.05648d0/T)                   &
          + 1.0d0*exp(-hc_over_k*573.59573d0/T)
      w1d = 5.0d0*exp(-hc_over_k*9238.6090d0/T)
      k = (w3p*log_log_table(T, T_z, k_3p) + w1d*log_log_table(T, T_z, k_1d)) &
          /internal_partition_function(iel_S, 0, T)
      end function sulfur_charge_transfer_ionization

      ! Na + H+ -> Na+ + H, the sum of the two channels above (row A21).
      pure double precision function sodium_charge_transfer_ionization(T) &
                                     result(k)
      real*8, intent(in) :: T
      k = sodium_nonradiative_ionization(T)                              &
        + sodium_radiative_charge_transfer(T)
      end function sodium_charge_transfer_ionization

      ! The rate coefficient [cm^3 s^-1] a source gives for row id at T [K].
      ! Every row reached here is carried as published; the reverses built
      ! by detailed balance (cx_reverse_of) and the rows that are not
      ! carried (A2, A6, A10, A22, C2, C4, C6, D18, D28, D30) never reach
      ! it, and a call for one stops the run.
      real*8 function published_rate_coefficient(id, T)
      integer, intent(in) :: id
      real*8,  intent(in) :: T
      real*8 :: t4, tr, rk, tc

      t4   = T/1.0d4
      tr   = T/300.0d0

      select case (id)
      ! ---------------- Group A: metal + H / H+ ----------------
      ! A1 Mg + H+ -> Mg+(3p, 4s) + H: Kingdon & Ferland (1996) Table 3
      ! row Mg0, a fit to Allan et al. (1988) Table 1 (read; within 1% of
      ! it at 1e4 K), fitted 5e3 - 3e4 K, evaluated as fitted outside. Its
      ! reverse A2 is not carried (the note at cx_A2_Mgp_H).
      case (1);  rk = kf(9.76d-12, 3.14d0,   55.54d0, -1.12d0,   0.0d0)
      ! A4 Mg2+ + H -> Mg+ + H+: Kingdon & Ferland (1996) Table 1 row
      ! Mg+2, a fit to the Landau-Zener rate of Butler & Dalgarno (1980),
      ! ApJ 241, 838, Table 1A (read: 8.6e-14 - 8.7e-14 at 1e3 - 3e4 K, no
      ! emission line listed for the product), fitted 1e3 - 3e4 K and
      ! evaluated as fitted outside. Its reverse A3 is detailed balance
      ! (and so shares that range); Butler & Dalgarno's own Table 2
      ! ionization rate, 7.6e-14 exp(-1.44 eV/kT), is that relation to 12%.
      case (4);  rk = kf(8.58d-14, 2.49d-3,  0.0293d0,-4.33d0,   0.0d0)
      ! A5 Fe + H+ -> Fe+ + H: Rutherford & Vroom (1972), J. Chem. Phys.
      ! 57, 3091, Table I (read), iron_charge_transfer_rate. Its reverse
      ! A6 is not carried (the note at cx_A2_Mgp_H).
      case (5);  rk = iron_charge_transfer_rate(1, T)                       ! A5  Fe + H+
      ! A8 Fe2+ + H -> Fe+ + H+: Kingdon & Ferland (1996) Table 1 row
      ! Fe+2 (Neufeld & Dalgarno 1987), fitted 1e3 - 1e5 K and evaluated
      ! as fitted outside. Its reverse A7 is detailed balance, as Kingdon &
      ! Ferland's Table 3 row Fe+ is.
      case (8);  rk = kf(1.26d-9,  0.0772d0, -0.41d0, -7.31d0,   0.0d0)
      ! A9 Si(3P) + H+ -> Si+(4P) + H: Glover & Jappsen (2007) Table 1
      ! R42, both branches as printed (5.88e-13 T^0.848 up to 1e4 K,
      ! 1.45e-13 T above; Huang et al. carry the first alone). R42 cites
      ! Kingdon & Ferland (1996), whose published tables carry no Si0 row;
      ! the fit reproduces within 9% the calculation of Kimura, Sannigrahi,
      ! Gu, Hirsch, Buenker & Shimamura (1996), ApJ 473, 1114, Table 1
      ! (read), 200 - 2e5 K, and is evaluated as fitted outside that range.
      ! Their channel ends in the metastable Si+ 3s3p2 4P term, 5.31 eV up,
      ! which decays to the ground term by the Si II] 2335 A
      ! intercombination lines; Si+(2P) formation is "much smaller". The
      ! reverse A10 is therefore NOT carried: detailed balance relates this
      ! rate to Si+(4P) + H, not to the ground-term Si+ it would be applied
      ! to (row A10 is excluded in cx_row_is_active).
      case (9)
         if (T .le. 1.0d4) then
            rk = 5.88d-13*T**0.848d0
         else
            rk = 1.45d-13*T
         endif
      ! A12 Si2+ + H -> Si+ + H+: Kingdon & Ferland (1996) Table 1 row
      ! Si+2, a = 1.23 (1e-9 cm^3 s^-1), b = 0.24, c = 3.17, d = +4.18e-3,
      ! fitted 1e1 - 1e6 K to Gargaud, McCarroll & Valiron (1982), A&A
      ! 106, 197, Table 3 (read; the fit is within 3% of it at 5e3, 1e4
      ! and 2e4 K). Glover & Jappsen (2007) R46 print the same. Huang et
      ! al. (2023) Table 4 prints 1.26e-9 and exp(-4.18e-3 T4). Gargaud et
      ! al. obtain their ionization rate from this one by detailed balance
      ! with the Si+ 2Po weight 6; the reverse A11 is formed that way.
      case (12); rk = kf(1.23d-9,  0.24d0,   3.17d0,   4.18d-3,  0.0d0)
      ! Oxygen O <-> H+ near-resonant charge exchange. Huang et al. (2023)
      ! Table 4 prints the two rate coefficients with the reactant labels
      ! exchanged: the exp(-227/T) Boltzmann factor is printed on the O+ + H0
      ! row, but IP(O I)=13.6181 eV > IP(H I)=13.5984 eV, so the O0 + H+ -> O+
      ! + H0 ionizing channel is the endothermic one (dE/k = 227.7 K) and must
      ! carry the barrier; O+ + H0 -> O0 + H+ is exothermic and must not.
      ! Assigned here to the physically correct rows (detailed balance then
      ! holds to ~9%, and each direction matches Cloudy c25.00 to 1-3%; as
      ! printed it violates detailed balance by 1.47x at 8000 K). See
      ! md/HUANG2023_TABLE4_OXYGEN_ERRATUM.md.
      case (13); rk = (1.26d-9*t4**0.517d0 + 4.25d-10*t4**6.69d-3)         &
                      *exp(-227.0d0/T)                                      ! A13 O + H+  (endothermic, ionizing)
      case (14); rk = 2.08d-9*t4**0.405d0 + 1.11d-11*t4**(-0.458d0)          ! A14 O+ + H  (exothermic, recombining)
      ! A15 C + H+ -> C+ + H (+ photon): carbon_charge_transfer_ionization.
      case (15); rk = carbon_charge_transfer_ionization(T)
      ! A16 C+ + H -> C + H+: reaction (1) of Stancil et al. (1998), their
      ! eq. (4) with the Table 1 parameters a1 = 6.08e-14, b1 = 1.96,
      ! c1 = 1.7e5 K, which Huang et al. (2023) print as 6.3e-17
      ! (T/300)^1.96 exp(-17/T4). It is the non-radiative channel, whose
      ! reverse is the non-radiative part of A15; Stancil et al. obtained
      ! that part from this one by detailed balance at the cross-section
      ! level, so the pair is theirs and consistent. The calculated range
      ! is 3e4 - 1e7 K (their Table 1 starts at 5e4 K); BELOW 3e4 K THE FIT
      ! IS EXTRAPOLATED, with no calculated point behind it. There it lies
      ! far below the detailed-balance reverse of the whole of A15 (by
      ! 1e-5 at 1e4 K), which is the bound it must respect, since A15 is
      ! mostly radiative there.
      case (16); rk = 6.08d-14*t4**1.96d0*exp(-1.7d5/T)
      ! A18 N+(3P) + H -> N(4S) + H+ + 0.94 eV: Kingdon & Ferland (1996)
      ! Table 1 row N+1, a = 1.01(-3) (1e-9 cm^3 s^-1), b = -0.29,
      ! c = -0.92, d = -8.38, fitted 1e2 - 5e4 K (read from the ADS scan),
      ! a fit to the distorted-wave calculation of Butler & Dalgarno
      ! (1979), ApJ 234, 765 (read): spin-orbit coupling of the X 2Pi and
      ! a 4Sigma- states of NH+, "1.0 x 10^-12 cm3 s^-1 at a temperature
      ! of 10^4 K" and 1.2e-12 at 1e3 K, their Fig. 2 spanning log T = 2 -
      ! 4.8, uncertain by "no more than 30%". The fit reproduces both
      ! values (1.01e-12, 1.19e-12). Evaluated as fitted outside
      ! 1e2 - 5e4 K (it falls as T^-0.29 there). N(4S) is the only N term
      ! the 0.94 eV can reach (N 2D lies at 2.38 eV, NIST), so the product
      ! is the ground term and the reverse A17 is detailed balance, which
      ! is also Kingdon & Ferland's own Table 3 row N0 (their 4.55(-3) is
      ! this 1.01(-3) times the ground-term weight ratio 9 x 2 / 4 = 4.5,
      ! with dE/k = 1.086e4 K).
      ! Huang et al. (2023) Table 4 print both directions as
      ! exp[sum c_i (ln T)^i] fits citing Lin, Stancil, Gu, Buenker &
      ! Kimura (2005), Phys. Rev. A 71, 062708, which could not be
      ! obtained. As printed the endothermic N + H+ fell from 4.1e-18 at
      ! 5e3 K to 8.1e-22 at 2e4 K, faster than T^(-3/2), which no cross
      ! section gives (the threshold law, module header), and at 1e4 K both rows
      ! were 2e7 and 70 times below the rates above; they are not carried.
      case (18); rk = kf(1.01d-12, -0.29d0,  -0.92d0,  -8.38d0,   0.0d0)
      ! A19 S + H+ -> S+ + H: sulfur_charge_transfer_ionization (Zhao et
      ! al. 2005, Table III). Its reverse A20 is detailed balance.
      case (19); rk = sulfur_charge_transfer_ionization(T)                 ! A19 S + H+
      ! A21 Na + H+ -> Na+ + H: sodium_charge_transfer_ionization, the sum
      ! of the non-radiative channel to H(n=2) (Dutta et al. 2001) and the
      ! radiative channel to H(1s) (Watanabe et al. 2002). Its reverse A22
      ! is not carried (the note at cx_A2_Mgp_H).
      case (21); rk = sodium_charge_transfer_ionization(T)                 ! A21 Na + H+
      ! A23 K + H+ -> K+ + H(1s) + photon: Watanabe et al. (2002), Fig. 5
      ! (potassium_radiative_charge_transfer). Radiative, so it has no
      ! collisional reverse and none is carried. The non-radiative
      ! K(4s) + H+ -> K+ + H(n=2), endothermic by 0.94 eV (NIST), was not
      ! calculated in either paper and is not carried.
      case (23); rk = potassium_radiative_charge_transfer(T)               ! A23 K + H+
      ! ---------------- Group B: He + H / H+ ----------------
      ! THE TWO ROWS ARE DIFFERENT CHANNELS. B2 is radiative and has no
      ! collisional reverse. B1 is the reverse of the non-radiative
      ! He+ + H channel that he_h_cx_rates adds to B2
      ! (nonradiative_electron_capture_Hep_from_H), and is formed from it.
      !
      !   B1  He(1s^2) + H+ -> He+ + H         NON-RADIATIVE (collisional)
      !       The detailed-balance reverse of the non-radiative channel
      !       (nonradiative_electron_capture_Hep_from_H), with the
      !       descriptors of row B2 (same reactants and products):
      !       k(B1) = k_nonrad 4 exp(-10.989 eV/kT) for the ground singlet.
      !       k_nonrad is the reciprocity transform of the H+ + He(1s^2)
      !       cross section of Loreau, Ryabchenko & Vaeck (2014), so k(B1)
      !       is the Maxwellian average of that computed cross section
      !       itself (with the near-threshold form of the note at the
      !       function); the
      !       detailed-balance form only moves the tabulation to the
      !       exothermic direction, where the table interpolates well.
      !       The rate table of an earlier calculation of this direction,
      !       Kimura, Lane, Dalgarno & Dixson (1993), ApJ 405, 801, Table 3
      !       (read), which Huang et al. (2023) Table 4 carry through the
      !       low branch of the Glover & Jappsen (2007) R27 fit, is not
      !       carried: it is 5100x the carried B1 at 1e4 K, 570x at 2e4 K
      !       and 560x at 1e5 K, and it is not the Maxwellian average of
      !       Kimura et al.'s own Table 1 cross sections either (that
      !       average is 6e-8 of Table 3 at 1e4 K and 1/220 at 1e5 K;
      !       DERIVED, the Table 1 energies read as c.m.). Below 2e4 K
      !       Kimura et al.'s rate is
      !       set by the cross section within a few eV of the 11 eV
      !       threshold, which the paper does not tabulate (its Table 1
      !       starts at 20 eV), and its table implies activation energies
      !       of 15 eV between 6000 and 8000 K and 4.3 eV between 8000 and
      !       1e4 K, neither the 11 eV of the threshold; the exothermic
      !       channel at those temperatures samples collision energies of a
      !       few eV, where Zygelman et al. computed the cross sections
      !       directly (their Fig. 8, 1 - 100 eV). The reverse of the
      !       better-determined direction is therefore the one carried.
      !       THE THRESHOLD LAW DECIDES IT. For an endothermic rate the
      !       Tolman activation energy E_a = -d ln k/d(1/kT) equals the mean
      !       energy of the reacting collisions less (3/2) kT, and every
      !       reacting collision carries at least the 10.989 eV, so
      !       E_a >= 10.989 eV - (3/2) kT, and a chord between two
      !       temperatures can be no lower than that bound at the upper one.
      !       Kimura et al.'s Table 3 (read from the ADS scan: 4.7e-22,
      !       7.4e-19, 2.6e-18, 2.7e-16, 3.7e-15 at 6e3, 8e3, 1e4, 2e4,
      !       4e4 K) gives 4.33 eV between 8e3 and 1e4 K against a bound of
      !       9.70 eV, and 8.00 against 8.40 between 1e4 and 2e4 K; its
      !       2-digit values cannot account for that (the first needs a
      !       factor 16.7 over the chord, the table has 3.5). The table is
      !       not consistent with its own Table 1 either: the charge-
      !       transfer cross section there is 3.7e-23 cm^2 at 20 eV and
      !       falls toward threshold (their Fig. 1), and holding it at
      !       3.7e-23 all the way down to 10.99 eV gives at most 2.4e-21,
      !       2.4e-17 and 7.4e-16 cm^3 s^-1 at 1e4, 4e4 and 1e5 K, 1e3,
      !       150 and 200 times below the Table 3 values (DERIVED). The
      !       computed cross sections of Loreau et al. (2014) and Kimura
      !       et al.'s Table 1 agree within 25% at 20 eV (c.m.) and fall
      !       steeply toward the threshold, as the law requires.
      !
      !   B2  He+ + H -> He + H+ + photon      RADIATIVE charge transfer
      !       Stancil, Lepp & Dalgarno (1998), ApJ 509, 1, Table 1 row (19),
      !       1.20E-15 (T/300)^0.25 (read from the page image), which prints
      !       the photon in the exit channel and whose note (10) reads
      !       "Zygelman et al. 1989, multiplied by 0.25 to account for
      !       approach probability factor (see Stancil & Zygelman 1996)";
      !       Stancil & Zygelman (1996), ApJ 472, 102, footnote 1: the
      !       radiative rate coefficients of Zygelman et al. (1989) "must be
      !       multiplied by a factor of 1/4 to properly account for the
      !       statistical weight of the A 1Sigma+ channel formed during the
      !       initial approach of the atoms". The 1.25e-15 carried before
      !       2026-09-26 is the value Glover & Jappsen (2007), ApJ 666, 1,
      !       print for their R26 citing Zygelman et al. (1989), 4% above
      !       the Stancil et al. fit. The underlying calculation is Zygelman,
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
      !       that is 1.17e-15, against 1.20e-15 from the fit (3%).
      !       The radiative cross sections of their Fig. 8 (triangles,
      !       about 3.5e-4 a0^2 from 1 to 20 eV, DERIVED from a reading of
      !       the figure, held below 1 eV) times 1/4 give 2.1e-15, 4.0e-15
      !       and 7.0e-15 at 3e3, 1e4 and 3e4 K against 2.1e-15, 2.9e-15
      !       and 3.8e-15 from the fit: 1.4x and 1.9x the fit above 1e4 K,
      !       where the recomputation of Courtney et al. (below) is 0.85x
      !       of it. Band carried above 1e4 K: 0.85x - 1.9x.
      !       BELOW 1000 K a digitization of Courtney et al. (2021) Fig. 3
      !       (the NLTE zero-density curve without stimulated processes,
      !       about 3% a point, made in the MoCHII tree on 2026-09-26)
      !       lies 1.9x above this fit at 100 K and 2.8x at 10 K, while it
      !       agrees within 2% over 1e3 - 1e4 K. A fit to their tabulated
      !       rate (University of Georgia Molecular Opacity Project) would
      !       replace this power law; the site no longer exists (2026-09-27). The same
      !       row is carried by MoCHII (agreed 2026-09-26).
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
      !       the figure), against 1.62e-15 and 2.88e-15 from the fit
      !       below: agreement to 1% and 15% across the range this code
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
      ! k(B1)/k(B2) is under no obligation to equal 4 exp(-dE/kT). The
      ! label in Table 4 drops the photon and with it the reason.
      case (24); rk = nonradiative_electron_capture_Hep_from_H(T)          &
                      *detailed_balance_ratio(cx_B2_Hep_H0, T)
      case (25); rk = 1.20d-15*tr**0.25d0                                    ! B2 He+ + H   (RADIATIVE, Stancil 1998)
      ! ---------------- Group C: metal + He / He+ ----------------
      ! C1 Si(3P) + He+ -> Si+* + He: Satta, Grassi, Gianturco, Yakovleva
      ! & Belyaev (2013), MNRAS 436, 2722 (read), a multiple-crossing
      ! Landau-Zener calculation, their eq. (8) k = c0 T^0.5 + c1 T +
      ! c2 T^-0.5 + c3 with the Table 1 coefficients of k1, "the single
      ! ionization channel": c0 = 3.226 28(-13), c1 = 1.327 26(-16),
      ! c2 = 4.1.533(-9) as printed (read here as 4.1533e-9; the other
      ! reading 4.11533e-9 moves k1 by 0.5% at 1e4 K), c3 = -9.162 43(-13),
      ! over the 10 - 1e5 K of their figures. Huang et al. (2023) print
      ! 3.32e-13, 1.2e-16, 4.2e-9 and -7.9e-13 (2% above at 1e4 K); the
      ! paper's values are carried. The product is an excited Si+* that
      ! radiates (the note at cx_A2_Mgp_H); C2 is not carried. Their k2,
      ! Si + He+ -> Si++ + He + e (electron shake-off), 1.5e-12 at 1e4 K
      ! against 7.4e-11 for k1, moves two electrons and has no place in the
      ! one-electron layout of this table; it is not carried.
      case (26); rk = 3.22628d-13*sqrt(T) + 1.32726d-16*T                  &
                    + 4.1533d-9/sqrt(T) - 9.16243d-13                        ! C1 Si + He+
      ! C3 C + He+ -> C+(2D) + He: Glover & Jappsen (2007) R41, their three
      ! branches (Huang et al. 2023 carry the T > 2000 K one alone, 1.4x
      ! low at 1140 K), a fit to the non-radiative calculation of Kimura,
      ! Dalgarno, Chantranupong, Li, Hirsch & Buenker (1993), ApJ 417, 812,
      ! Table 1, 10 - 1e4 K (read; the branches reproduce it at 1e3 and
      ! 1e4 K). Above 1e4 K the last branch is extrapolated.
      case (28)
         if (T .le. 2.0d2) then
            rk = 8.58d-17*T**0.757d0
         else if (T .le. 2.0d3) then
            rk = 3.25d-17*T**0.968d0
         else
            rk = 2.77d-19*T**1.597d0
         endif
      ! C5 O + He+ -> O+ + He + photon: Zhao et al. (2004), ApJ 615, 1063,
      ! Table 3 and eq. (19), sum over i of a_i (T/1e4)^b_i exp(-T/c_i),
      ! a1 = 4.991e-15, b1 = 0.3794, c1 = 1.121e6 K, a2 = 2.780e-15,
      ! b2 = -0.2163, c2 = -8.158e5 K (read), 10 - 1e6 K, within 6% of
      ! their table. Huang et al. (2023) drop the first exp and print
      ! 81.97 for c2/1e4. RADIATIVE, so its reverse C6 is not carried.
      ! The temperature is held inside the calculated 10 - 1e6 K: above
      ! 1e6 K the second term, exp(+T/8.158e5), would grow without bound.
      case (30)
         tc = min(max(T, 1.0d1), 1.0d6)
         rk = 4.991d-15*(tc/1.0d4)**0.3794d0*exp(-tc/1.121d6)             &
            + 2.780d-15*(tc/1.0d4)**(-0.2163d0)*exp(tc/8.158d5)
      ! ---------------- Group D: metal + metal ----------------
      ! D1/D2 C + Si+ <-> C+ + Si: Huang et al. (2023) Table 4 cite Satta
      ! et al. (2013), which (read) computes Si + He+ only and gives no
      ! C + Si+ rate; the origin of these rows is unknown and they are
      ! carried as printed.
      case (32); rk = (0.724d0*T**0.0463d0*exp(-3.61d0/t4))/dcsi(T)          ! D1  C + Si+
      case (33); rk = 1.0d0/dcsi(T)                                          ! D2  C+ + Si
      ! The metal + metal forwards of Prasad & Huntress (1980), ApJS 43, 1,
      ! Table 2 (read), rate A (T/300)^alpha exp(-beta/T) with alpha =
      ! beta = 0: reactions 243 (D3), 241 (D5), 274 (D7), 138 (D9), 188
      ! (D11), 273 (D13), 136 (D15), 254 (D19), 253 (D21), 252 (D23). Their
      ! text says the rates of the ion-molecule reactions not studied in
      ! the laboratory "were estimated by the use of the guidelines
      ! described by Huntress (1977)"; 188 is from the review of Ferguson
      ! (1973). They are temperature-independent values for interstellar
      ! clouds, evaluated here at 1e3 - 2e4 K without a calculation or a
      ! measurement there. Each reverse (the even row after it) is
      ! detailed balance.
      case (34); rk = 2.9d-9                                                 ! D3  Mg + Si+
      case (36); rk = 1.9d-9                                                 ! D5  Fe + Si+
      case (38); rk = 1.0d-11                                                ! D7  Na + Mg+
      case (40); rk = 1.1d-9                                                 ! D9  Mg + C+
      case (42); rk = 1.2d-9                                                 ! D11 Mg + N+
      case (44); rk = 1.0d-11                                                ! D13 Na + Fe+
      case (46); rk = 2.6d-9                                                 ! D15 Fe + C+
      case (50); rk = 2.8d-10                                                ! D19 Mg + S+
      case (52); rk = 2.6d-10                                                ! D21 Na + S+
      case (54); rk = 1.8d-10                                                ! D23 Fe + S+
      ! D17 Fe + O+ and D27 Fe + N+: Rutherford & Vroom (1972) Table I
      ! (read), iron_charge_transfer_rate. D18 and D28 are not carried
      ! (the note at cx_A2_Mgp_H).
      case (48); rk = iron_charge_transfer_rate(2, T)                       ! D17 Fe + O+
      case (58); rk = iron_charge_transfer_rate(3, T)                       ! D27 Fe + N+
      ! D25 S(3P) + C+(2P) -> S+(4S) + C(3P): Chenel, Mangaud, Justum,
      ! Talbi, Bacchus-Montabonel & Desouter-Lecomte (2010), J. Phys. B 43,
      ! 245701 (read), Table 2, "Wave packet" column, the Maxwellian average
      ! of cross sections from quantum wave-packet dynamics over 0.5 - 10 eV
      ! completed by semiclassical values above ("the most precise
      ! evaluation" in their words): 1.3e-11, 3.7e-11, 8.7e-11, 1.0e-10,
      ! 1.5e-10, 2.2e-10 at 500, 1e3, 5e3, 1e4, 5e4, 1e5 K, interpolated in
      ! log-log and held outside. Semiclassical and wave-packet values
      ! differ by "about 10-20% for most of the temperatures with a higher
      ! error of 37% at 10 000 K". Huang et al. (2023) print 3.26e-11
      ! (T/300)^0.289 exp(-338/T), 0.89 of the table at 1e3 K and 0.87 at
      ! 1e4 K (its exp(-338/T) is a fit parameter, not a threshold of this
      ! exothermic reaction). The reverse D26 is detailed balance
      ! (cx_reverse_of), as Chenel et al. obtain theirs ("deduced from the
      ! symmetry properties of the S-matrix").
      case (56); rk = log_log_table(T, [5.0d2, 1.0d3, 5.0d3, 1.0d4,        &
                                        5.0d4, 1.0d5],                       &
                                    [1.3d-11, 3.7d-11, 8.7d-11, 1.0d-10,     &
                                     1.5d-10, 2.2d-10])                      ! D25 S + C+
      ! D29 Na + Ca+ -> Na+ + Ca: Kwolek et al. (2019), Phys. Rev. A 99,
      ! 052703 (read), a hybrid-trap measurement at ion energies of 0.04 -
      ! 0.24 eV (their Fig. 6; relative collision energies 0.37 of that,
      ! about 0.015 - 0.09 eV, i.e. 3kT/2 at 100 - 700 K). The ground-state
      ! channel Na[S] + Ca+[S] rises from 4e-11 to 2e-9 and 4e-9 cm^3 s^-1
      ! over that range and the Ca+[D] channels are 2.5e-9 - 3.5e-9; all
      ! lie at the Langevin capture rate they give, "kL(g) ~ 2.995 x 10^-9
      ! cm3/s" for Na[3S] (their Eq. (15); Na[P] + Ca+[D] is above it).
      ! The printed 3.0e-9 is that Langevin value and is carried: a capture
      ! rate independent of energy, the classical upper bound, which the
      ! measurement reaches at its highest energies. Above the measured
      ! energies it is an assumption (the channel they infer is
      ! endothermic by 0.9 eV asymptotically, open at wind temperatures).
      ! The product is Ca 4s4p 3P (the note at cx_A2_Mgp_H); D30 is not
      ! carried.
      case (60); rk = 3.0d-9                                                 ! D29 Na + Ca+
      ! D31 K + Na+ -> K+ + Na: Lavvas, Koskinen & Yelle (2014), ApJ 796,
      ! 15, Table 2 R16 (read), k = 1e-11, marked "Est" (an estimate). Its
      ! reverse D32 is detailed balance; Lavvas et al.'s own reaction
      ! constant for it, 0.7844 T^0.032 exp(-9244/T), is that relation to
      ! within 5% at 1e4 K.
      case (62); rk = 1.0d-11                                                ! D31 K + Na+
      ! ---------------- Group E: outside Huang Table 4 ----------------
      case (64); rk = cx_o2p_h_scale*electron_capture_O2p_from_H(T)         ! E1 O2+ + H0
      case (65); rk = cx_n2p_h_scale*electron_capture_N2p_from_H(T)         ! E2 N2+ + H0
      case default
         write(*,*) '(charge_exchange) ERROR: row', id, ' has no',        &
                    ' published expression: it is a detailed-balance',    &
                    ' reverse or a row that is not carried.'
         write(*,*) '  Aborting.'
         error stop 1
      end select
      published_rate_coefficient = rk

      contains

         ! KF96 bracket form: a * T4^b * (1 + c*exp(d*T4)) * exp(-Ek/T4).
         ! d carries the table's multiplicative-T4 coefficient with its
         ! printed sign; Ek is the trailing activation barrier in units of
         ! 1e4 K (Ek=0 if no trailing Boltzmann factor).
         real*8 function kf(a, b, c, d, Ek)
         real*8, intent(in) :: a, b, c, d, Ek
         kf = a * t4**b * (1.0d0 + c*exp(d*t4)) * exp(-Ek/t4)
         end function kf

         ! Power-law form: a * (T/300)^b * exp(-Ek/T4).
         real*8 function gp(a, b, Ek)
         real*8, intent(in) :: a, b, Ek
         gp = a * tr**b * exp(-Ek/t4)
         end function gp

         ! Denominator D_CSi(T) used by the C + Si charge-transfer rows
         ! (Satta 2013).
         real*8 function dcsi(Tk)
         real*8, intent(in) :: Tk
         dcsi = 1.87d8 + 5.09d10*Tk**(-0.527d0)
         end function dcsi

      end function published_rate_coefficient

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
      ! Dedicated He <-> H charge exchange (Huang 2023, Table 4, group B,
      ! and the He2+ + H capture Table 4 does not carry). Applied by every
      ! ionization system that contains He, so it is not tied to the metal
      ! charge-exchange assembly (cx_add_to_fvec).
      !   B1: He0 + H+ -> He+ + H0   (k1, endothermic: the detailed-balance
      !       reverse of the non-radiative part of k2, 4 exp(-10.989 eV/kT))
      !   B2: He+ + H0 -> He0 + H+   (k2, exothermic: radiative plus
      !       non-radiative)
      !   B3: He2+ + H0 -> He+ + H+ + photon  (k3, radiative, no reverse;
      !       radiative_electron_capture_Hepp_from_H)
      ! B1 ionizes He (HeI->HeII) and recombines H (HII->HI); B2 ionizes H
      ! (HI->HII) and recombines He (HeII->HeI); B3 ionizes H and
      ! recombines He (HeIII->HeII). The three rate coefficients come from
      ! he_h_cx_rates. Gated by he_h_charge_exchange (default .true.): each
      ! routine adds nothing when the switch is off.
      ! ================================================================

      ! Rate coefficients [cm^3 s^-1] of the He<->H pair at temperature T [K].
      ! k_He0_Hp is B1 (He0+H+), k_Hep_H0 is B2 (He++H0). The formulas live in
      ! one place (cx_rate rows 24,25 and the function below), so this
      ! reuses them.
      ! Non-radiative electron capture by singly ionized helium from atomic
      ! hydrogen,
      !     He+(1s) + H(1s) -> He(1s^2 1S) + H+ + 10.989 eV ,
      ! rate coefficient [cm^3 s^-1] at T [K]. he_h_cx_rates adds it to the
      ! radiative channel (row B2); row B1 is its detailed-balance reverse.
      !
      ! SOURCE. Loreau, Ryabchenko & Vaeck (2014), J. Phys. B 47, 135204,
      ! read from the published pages. They compute the reverse reaction,
      ! H+ + He(1s^2) -> H(1s) + He+(1s) (their reaction (1)), with a
      ! time-dependent and a time-independent quantal method on ab initio
      ! HeH+ potentials and couplings ("excellent agreement" between the
      ! two, and with the low-energy measurements of Latypov & Shaporenko),
      ! and recommend the fit of their eq. (3), log10 sigma[1e-16 cm^2] as
      ! a Fourier series in log10 E[eV/u], Table 2 coefficients, "valid for
      ! collision energies from 10 eV/u up to 10 MeV/u", reproducing their
      ! cross sections "with an uncertainty of less than 25% between 100
      ! and 400 eV/u ... and of less than 5% over the rest". The forward
      ! cross section follows by microscopic reversibility,
      !     sigma_f(E_f) = (1/4) (mu_r E_r) / (mu_f E_f) sigma_r(E_r) ,
      !     E_r = E_f + 10.989 eV (c.m.), statistical weights 4 for
      !     He+(1s) + H(1s) and 1 for He(1^1S) + H+ ,
      ! with E[eV/u] fixing the relative speed, E_cm = mu[u] E, so that the
      ! reverse threshold lies at x_th = 13.655 eV/u. The table below is
      ! the Maxwellian average of this sigma_f (DERIVED, mu_f = 0.8051 u;
      ! generator cooling_data/heii_h_charge_exchange_rate.py). Their
      ! calculation starts at 15 eV/u (their Fig. 1) and the fit does not
      ! follow the threshold (it stays finite under it, 1.3e-23 cm^2 at
      ! 10 eV/u, where the channel is closed), so between x_th and 15 eV/u
      ! the reverse cross section is taken as
      !     sigma_r(x) = sigma_r(15) sqrt((x - x_th)/(15 - x_th)),
      ! zero below x_th: by reciprocity, sigma_f v constant at low forward
      ! energy, capture on the ion-induced dipole potential with a transfer
      ! probability independent of energy (an assumption over forward
      ! energies 0 - 1.08 eV that no calculation tests). The same form is
      ! carried by MoCHII (2026-09-26), so the two codes agree.
      ! Interpolated linearly in ln k against 1/T; the first chord is
      ! continued below 300 K and the 1e5 K value is held above 1e5 K.
      !
      ! UNCERTAINTY, as factors on the carried value (DERIVED from the
      ! published cross sections):
      !   - Near threshold, T < 3e4 K. The two treatments of 13.655 -
      !     15 eV/u that bracket the carried one: the fit carried through
      !     the threshold (finite there, so sigma_f ~ 1/E_f and a rate
      !     growing without bound as T falls) gives 1.32x at 1e4 K, 2.1x at
      !     3e3 K and 3.6x at 1e3 K; sigma_r = 0 below 15 eV/u, carried
      !     here before 2026-09-26, closes a channel that is open and gives
      !     0.46x at 1e4 K, 0.039x at 3e3 K and 1.5e-5x at 1e3 K. At the
      !     carried value the channel is 3.3% of the radiative one at 1e4 K
      !     and 5.9% at 1e3 K.
      !   - Above 3e4 K, where the calculations disagree. Kimura, Lane,
      !     Dalgarno & Dixson (1993), ApJ 405, 801, Table 1 (the reverse
      !     cross section at 20 - 5000 eV, energies read as c.m.), turned
      !     forward the same way, gives 0.4x at 4e4 K, 1.0x at 6e4 K and
      !     2.5x at 1e5 K (its cross sections are 5-8x Loreau's at 50 -
      !     100 eV and agree within 25% at 20 eV). Band carried: 0.4x -
      !     2.5x. The non-radiative cross sections of Zygelman, Dalgarno,
      !     Kimura & Lane (1989), Phys. Rev. A 40, 2340, Fig. 8 (circles,
      !     5 - 51 eV, DERIVED from a reading of the figure), computed on
      !     the A 1Sigma+ entrance channel alone and so conditional on a
      !     singlet approach, times the 1/4 of that approach, give 0.9x at
      !     1e4 K, 26x at 2e4 K, 100x at 3e4 K, 350x at 5e4 K and 790x at
      !     1e5 K. They are not carried: at the cross-section level they
      !     exceed the reciprocity transform of Loreau et al. by 12x at
      !     5.5 eV rising to 930x at 51 eV, and the reverse cross sections
      !     of Kimura et al. (1993, the same group) by 230 - 370x at 20
      !     and 50 eV, while the two reverse calculations agree with each
      !     other to within a factor 8.
      ! Relative to the radiative channel B2 the carried rate is 3.3% at
      ! 1e4 K, 2.6% at 3e4 K and 4.7% at 1e5 K: He+ + H is removed almost
      ! entirely by radiative charge transfer in this whole range.
      !
      ! Kingdon & Ferland (1996), ApJS 106, 205, Table 1 row He+ (a =
      ! 7.47e-6, b = 2.06, c = 9.93, d = -3.89, 6e3 - 1e5 K), carried here
      ! before 2026-09-26, is a fit to the Zygelman et al. values WITHOUT
      ! the 1/4: it lies within 13 - 30% of the singlet-conditional
      ! Maxwellian average at 5e4 - 1e5 K and 55x above it at 1e4 K, which
      ! puts it 93x above the carried value at 1e4 K and 3500x at 1e5 K
      ! (DERIVED).
      double precision function nonradiative_electron_capture_Hep_from_H(T) &
                                result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_l = 23
      real*8, parameter  :: T_l(n_l) = [ 3.0d2, 5.0d2, 7.0d2, 1.0d3,     &
         1.5d3, 2.0d3, 3.0d3, 4.0d3, 5.0d3, 6.0d3, 7.0d3, 8.0d3, 1.0d4,  &
         1.2d4, 1.5d4, 2.0d4, 2.5d4, 3.0d4, 4.0d4, 5.0d4, 6.0d4, 8.0d4,  &
         1.0d5 ]
      real*8, parameter  :: k_l(n_l) = [ 9.484d-17, 9.506d-17,           &
         9.528d-17, 9.562d-17, 9.617d-17, 9.670d-17, 9.754d-17,          &
         9.797d-17, 9.805d-17, 9.789d-17, 9.758d-17, 9.721d-17,          &
         9.644d-17, 9.579d-17, 9.522d-17, 9.541d-17, 9.693d-17,          &
         9.962d-17, 1.081d-16, 1.205d-16, 1.368d-16, 1.811d-16,          &
         2.418d-16 ]
      real*8  :: x, w, lnk
      integer :: i
      if (T .le. 0.0d0) then
         k = 0.0d0
         return
      endif
      x = 1.0d0/min(T, T_l(n_l))
      i = 1
      do while (i .lt. n_l - 1 .and. min(T, T_l(n_l)) .gt. T_l(i+1))
         i = i + 1
      enddo
      w = (x - 1.0d0/T_l(i))/(1.0d0/T_l(i+1) - 1.0d0/T_l(i))
      lnk = log(k_l(i)) + w*(log(k_l(i+1)) - log(k_l(i)))
      k = exp(max(lnk, -700.0d0))
      end function nonradiative_electron_capture_Hep_from_H

      ! Radiative electron capture by doubly ionized helium from atomic
      ! hydrogen,
      !     He2+ + H(1s) -> He+(1s) + H+ + photon (40.81 eV defect) ,
      ! rate coefficient [cm^3 s^-1] at T [K] (row B3).
      !
      ! SOURCE. West, Lane & Cohen (1982), Phys. Rev. A 26, 3164, read from
      ! the published pages: the cross section for radiative transfer to
      ! He+(1s) over 0.1 - 1000 eV (c.m.), their Fig. 4 (linear axis) and
      ! Fig. 6 (log axes), with the low-energy structure in Fig. 5 (0.01 -
      ! 0.3 eV; predissociation resonances, of which "only the widest" are
      ! resolved). The paper gives no rate coefficient. The table below is
      ! the Maxwellian average of that cross section (DERIVED): Fig. 5(a)
      ! from 0.0125 eV, Fig. 5(b) above 0.1 eV, Fig. 6 above 0.3 eV and
      ! Fig. 4 above 4.5 eV, digitized (Figs. 4 and 6 agree within 3% over
      ! 0.15 - 3 eV), the resolved resonances included (they raise the
      ! energy average of the cross section by 6% over 0.01 - 0.1 eV and
      ! 12% over 0.1 - 0.3 eV, and the rate by 8% at 1e3 K, 1.5% at
      ! 1e4 K), power laws continued from the neighboring decade outside
      ! 0.0125 - 60 eV. Log-log interpolated, held at the end points,
      ! 100 K - 1e5 K.
      !
      ! RANGE AND UNCERTAINTY. Digitization about 3%; resonances narrower
      ! than those resolved are not in Fig. 5 and could add a few percent
      ! below 3e3 K; the part of the average below the figure's lowest
      ! energy (a power-law continuation) is 3% of the rate at 1e3 K, 17%
      ! at 300 K and 56% at 100 K, so the table is a calculation above
      ! 1e3 K and an extrapolation below. The non-radiative channel to
      ! He+(n = 2), which dominates above a few hundred eV, is not carried:
      ! West et al. expect it to "decrease monotonically at low energies"
      ! and the radiative process to dominate below about 10 eV (c.m.);
      ! the Coulomb-trajectory values they plot in Fig. 6 (Winter & Lane;
      ! 0.034 a0^2 at 100 eV, falling steeply below) give below 1e-15 at
      ! 1e5 K, 0.5% of this rate, if taken as zero below 100 eV (DERIVED).
      ! Over 1e3 - 3e4 K the rate is 1.6 - 2.0e-13, 16 - 20 times the
      ! constant 1.0e-14 of Kingdon & Ferland (1996) Table 1 row He+2
      ! carried before 2026-09-26. That constant is not an average of these
      ! cross sections: their Section 2 sets every reaction without a
      ! significant non-radiative channel to "1 x 10^-14 (see Butler,
      ! Guberman, & Dalgarno 1977)", and the He+2 row carries the same
      ! 1.00(-5) and footnote as Li+3, Be+2, B+3 and Ne+2.
      ! The photon carries the whole defect: no collisional reverse, and no
      ! heat to the gas (charge_exchange_heating).
      double precision function radiative_electron_capture_Hepp_from_H(T) &
                                result(k)
      real*8, intent(in) :: T
      integer, parameter :: n_w = 13
      real*8, parameter  :: T_w(n_w) = [ 1.0d2, 3.0d2, 1.0d3, 2.0d3,     &
         3.0d3, 5.0d3, 8.0d3, 1.0d4, 1.5d4, 2.0d4, 3.0d4, 5.0d4, 1.0d5 ]
      real*8, parameter  :: k_w(n_w) = [ 1.684d-13, 1.856d-13,           &
         2.015d-13, 1.976d-13, 1.900d-13, 1.785d-13, 1.688d-13,          &
         1.651d-13, 1.600d-13, 1.579d-13, 1.572d-13, 1.600d-13,          &
         1.706d-13 ]
      if (T .le. 0.0d0) then
         k = k_w(1)
      else
         k = log_log_table(T, T_w, k_w)
      endif
      end function radiative_electron_capture_Hepp_from_H

      subroutine he_h_cx_rates(T, k_He0_Hp, k_Hep_H0, k_Hepp_H0)
      ! The three H <-> He rate coefficients, in the same admissible band
      ! [0, cx_kc_max] that cx_set_cell imposes on every other row: Group B
      ! is excluded from cx_act and reaches the systems only through here.
      ! None depends on the electron density: the products of the forwards
      ! are He(1^1S), H+, He+(1s) and H(1s), none of which has a metastable
      ! level in the manifolds of this module.
      real*8, intent(in)  :: T
      real*8, intent(out) :: k_He0_Hp, k_Hep_H0, k_Hepp_H0
      ! He+ + H removal is the SUM of the two channels that do it: the
      ! radiative one of row B2 and the non-radiative one above. He + H+ is
      ! the reverse of the non-radiative one alone (row B1).
      k_He0_Hp = min(max(published_rate_coefficient(cx_B1_He0_Hp, T),   &
                         0.0d0), cx_kc_max)
      k_Hep_H0 = min(max(published_rate_coefficient(cx_B2_Hep_H0, T)     &
                         + nonradiative_electron_capture_Hep_from_H(T),  &
                         0.0d0), cx_kc_max)
      k_Hepp_H0 = min(max(radiative_electron_capture_Hepp_from_H(T),    &
                          0.0d0), cx_kc_max)
      end subroutine he_h_cx_rates

      ! Residual terms for the systems whose H row (fvec(1)) is written
      ! HI->HII (ionization) positive and He row (fvec(2)) carries the
      ! HeI<->HeII balance.
      ! THE HELIUM REACTANT n_hei IS THE GROUND SINGLET He(1^1S), NOT THE
      ! SUM OVER He I. The rate of He + H+ -> He+ + H (he_h_cx_rates) is the
      ! detailed-balance reverse of He+ + H -> He(1^1S) + H+ and carries
      ! the barrier exp(-10.989 eV/kT), the ionization-potential difference
      ! 24.587 - 13.598 eV of ground-state helium against hydrogen, with the
      ! weight of the singlet alone. He(2^3S) sits 19.82 eV
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
      ! H row is HI->HII positive in all of them. k1 = He0+H+, k2 = He++H0,
      ! k3 = He2+ + H0.
      ! THE He2+ + H CAPTURE WRITES ROW 3, AND ROW 2 DEPENDS ON THE BASIS.
      ! Row 3 is He II -> He III positive in every system (a boundary flow,
      ! which for the top stage is also the He III source), so the capture
      ! enters it as -R3. Row 2 is a BOUNDARY FLOW (He I <-> He II) in the
      ! atomic systems, which the capture does not cross, and the He II
      ! STAGE SOURCE in the molecular ones (System_HeH_mol and its metals
      ! form, where src(He II) = fvec(2)), where it enters as +R3;
      ! heii_row_is_stage_source says which. The H row gains +R3 in both.
      subroutine he_h_cx_fvec(fvec, k1, k2, k3, n_hi, n_hii, n_hei,       &
                              n_heii, n_heiii, he_row_sign,               &
                              heii_row_is_stage_source, gross)
      real*8              :: fvec(*)
      ! The gross rate of rows 1-3, to which each reaction adds its
      ! magnitude on the rows it reaches (System_HeH_mol).
      real*8, optional, intent(inout) :: gross(*)
      real*8, intent(in)  :: k1, k2, k3, n_hi, n_hii, n_hei, n_heii
      real*8, intent(in)  :: n_heiii, he_row_sign
      logical, intent(in) :: heii_row_is_stage_source
      real*8 :: R1, R2, R3
      if (.not. he_h_charge_exchange) return
      R1 = k1*n_hei*n_hii       ! He0 + H+ -> He+ + H0
      R2 = k2*n_heii*n_hi       ! He+ + H0 -> He0 + H+
      R3 = k3*n_heiii*n_hi      ! He2+ + H0 -> He+ + H+ + photon
      ! H row (HI->HII positive): B2 and B3 ionize H0 (+R2, +R3), B1
      ! recombines H+ (-R1).
      fvec(1) = fvec(1) + R2 - R1 + R3
      ! He row: B1 ionizes HeI->HeII (+R1), B2 recombines HeII->HeI (-R2), in
      ! the ionization-positive orientation; he_row_sign flips it for the
      ! HeI-gain (TR) orientation.
      fvec(2) = fvec(2) + he_row_sign*(R1 - R2)
      if (heii_row_is_stage_source) fvec(2) = fvec(2) + R3
      fvec(3) = fvec(3) - R3
      if (present(gross)) then
         gross(1) = gross(1) + abs(R1) + abs(R2) + abs(R3)
         gross(2) = gross(2) + abs(R1) + abs(R2)
         if (heii_row_is_stage_source) gross(2) = gross(2) + abs(R3)
         gross(3) = gross(3) + abs(R3)
      endif
      end subroutine he_h_cx_fvec

      ! Analytic-Jacobian counterpart of he_h_cx_fvec for the systems that
      ! carry a written Jacobian (System_HeH, System_HeH_metals), whose He row
      ! is HeI->HeII positive (he_row_sign = +1). Unknown layout
      ! x1 = n_HII/n_H, x2 = n_HeII/n_He, x3 = n_HeIII/n_He, so the pair
      ! touches rows 1,2 and columns 1,2,3 only. nsz = system size.
      subroutine he_h_cx_jac(nsz, fjac, k1, k2, k3, n_h, n_he,            &
                             n_hi, n_hii, n_hei, n_heii, n_heiii)
      integer, intent(in) :: nsz
      real*8              :: fjac(nsz,nsz)
      real*8, intent(in)  :: k1, k2, k3, n_h, n_he, n_hi, n_hii, n_hei
      real*8, intent(in)  :: n_heii, n_heiii
      real*8 :: j1, j2, j3, d1, d3
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
      ! He2+ + H0: R3 = k3 n_heiii n_hi with dn_hi/dx1 = -n_h and
      ! dn_heiii/dx3 = n_he; +R3 in row 1, -R3 in row 3 (boundary basis,
      ! the only one these two systems use).
      d1 = -k3*n_heiii*n_h                         ! dR3/dx1
      d3 =  k3*n_he*n_hi                           ! dR3/dx3
      fjac(1,1) = fjac(1,1) + d1
      fjac(1,3) = fjac(1,3) + d3
      fjac(3,1) = fjac(3,1) - d1
      fjac(3,3) = fjac(3,3) - d3
      end subroutine he_h_cx_jac

      ! Residual terms for the advection-correction systems
      ! (System_implicit_adv_HeH / _HeH_TR). Those rows are fraction-
      ! normalized and written neutral-gain positive: fvec(1) tracks n_HI
      ! (per n_H), fvec(2) tracks n_HeI (per n_He = heh_loc*n_H). c1 = dr/v.
      ! xhi/xhii = HI/HII fractions of H; xhei/xheii = HeI/HeII fractions of
      ! He. k1 = He0+H+, k2 = He++H0.
      ! THE HELIUM REACTANT xhei IS THE GROUND SINGLET He(1^1S), as it is in
      ! he_h_cx_fvec and for the same reason: the rate carries the barrier
      ! exp(-10.989 eV/kT), the ionization-potential difference of
      ! ground-state helium against hydrogen, and the
      ! metastable's own charge exchange is a different reaction that this
      ! code does not carry. System_implicit_adv_HeH has no metastable and
      ! its neutral helium is the singlet; System_implicit_adv_HeH_TR passes
      ! x(2), the singlet unknown of its row 2.
      ! Row 3 of both systems is the He III fraction balance, written He
      ! III-gain positive (per n_he), which He2+ + H0 -> He+ + H+ empties.
      subroutine he_h_cx_fvec_adv(fvec, c1, c1_he, xhi, xhii, xhei,       &
                                  xheii, xheiii, heh_loc, n_h, k1, k2, k3)
      ! c1 and c1_he are the rate weights of the hydrogen and helium rows
      ! (g h_j over the velocity of that element's nuclei, ion_cell_state).
      real*8              :: fvec(*)
      real*8, intent(in)  :: c1, c1_he, xhi, xhii, xhei, xheii, xheiii
      real*8, intent(in)  :: heh_loc
      real*8, intent(in)  :: n_h, k1, k2, k3
      real*8 :: rr1, rr2, rr3
      if (.not. he_h_charge_exchange) return
      ! rr1 = R1/n_he = k1*n_hei*n_hii/n_he, rr2 = R2/n_he. Then R/n_h = rr*heh_loc.
      rr1 = k1*xhei*xhii*n_h        ! He0 + H+ -> He+ + H0, per n_he
      rr2 = k2*xheii*xhi*n_h        ! He+ + H0 -> He0 + H+, per n_he
      rr3 = k3*xheiii*xhi*n_h       ! He2+ + H0 -> He+ + H+, per n_he
      ! H row (HI-gain positive): B1 makes HI (+R1), B2 and B3 destroy HI
      ! (-R2, -R3); /n_h.
      fvec(1) = fvec(1) + c1*heh_loc*(rr1 - rr2 - rr3)
      ! He row (HeI-gain positive): B2 makes HeI (+R2), B1 destroys HeI (-R1); /n_he.
      fvec(2) = fvec(2) + c1_he*(rr2 - rr1)
      ! He III row (HeIII-gain positive): B3 destroys HeIII (-R3); /n_he.
      fvec(3) = fvec(3) - c1_he*rr3
      end subroutine he_h_cx_fvec_adv

      ! ================================================================
      ! CHARGE-EXCHANGE HEATING
      ! ================================================================

      ! Volumetric heat [erg cm^-3 s^-1] that the charge-exchange reactions
      ! this module applies deliver to the translational energy of the gas,
      ! cell by cell: the sum over the active reactions of R (dE - E_rad),
      ! R = k n(donor) n(acceptor), dE the energy defect (positive
      ! exothermic, cx_energy_defect_eV) and E_rad the excitation of the
      ! product state that leaves as radiation (cx_radiated_energy_eV),
      ! plus the He <-> H pair of he_h_cx_rates. The ionization energy the
      ! reactions move between species is not an energy of the gas in
      ! this code (photoionization deposits h nu - IP, recombination
      ! radiates IP), so this is the whole energy change of the gas and the
      ! ledger closes reaction by reaction: heat + E_rad = dE.
      !
      ! The reactions and rates are those of the balance rows: the active
      ! set cx_act of cx_init (empty when cx_init was not called, as in a
      ! run without metals), each coefficient cx_rate clamped into
      ! [0, cx_kc_max] as cx_set_cell clamps it, and the He <-> H pair
      ! under he_h_charge_exchange with the coefficients of he_h_cx_rates,
      ! whose He + H+ reactant is the ground singlet. It evaluates the
      ! rates itself and leaves the thread's cx_kc untouched, so it may be
      ! called anywhere. For a pair formed by detailed balance the reverse
      ! takes back exactly what the forward gives at the collisional limit
      ! (both use -dE and +dE, the forward ending in the populated manifold
      ! of its product, the note at q_level_energy), so a pair in balance
      ! deposits nothing there; below the critical densities the product's
      ! metastable excitation partly leaves in its forbidden line
      ! (product_manifold_radiation_eV), and a pair in balance then cools
      ! the gas by that line.
      !
      ! nm(:, im) are the metal ion stage densities in the order of
      ! species_table (mion_elem, mion_stage); n_hei is the ground-singlet
      ! He I, the reactant the He rows carry; n_e is the electron density
      ! the balance rows' rate coefficients are formed with (cx_set_cell),
      ! which also sets how much of a metastable product's excitation the
      ! gas keeps (product_manifold_radiation_eV).
      subroutine charge_exchange_heating(T_K, n_e, n_hi, n_hii, n_hei,    &
                                         n_heii, n_heiii, nm, heat)
      real*8, intent(in)  :: T_K(:), n_e(:), n_hi(:), n_hii(:), n_hei(:)
      real*8, intent(in)  :: n_heii(:), n_heiii(:)
      real*8, intent(in)  :: nm(:,:)
      real*8, intent(out) :: heat(:)
      real*8  :: dens(12,0:2), T, q, n_pair, kc, k_nr, k_b1, k_b2, k_b3
      real*8  :: lvl_p(5,3), lvl_X(5,3)
      logical :: lvl_done(3)
      integer :: j, i, r, im
      real*8, parameter :: dE_He_H = e_th_HeI - e_th_HI

      do j = 1, size(T_K)
         T = T_K(j)
         dens = 0.0d0
         do im = 1, n_mion
            dens(mion_elem(im), mion_stage(im)) = nm(j,im)
         enddo
         dens(cx_H,0)  = n_hi(j)
         dens(cx_H,1)  = n_hii(j)
         dens(cx_He,0) = n_hei(j)
         dens(cx_He,1) = n_heii(j)
         dens(cx_He,2) = n_heiii(j)
         q = 0.0d0
         lvl_done = .false.
         do i = 1, cx_nact
            r = cx_act(i)
            n_pair = dens(cx_don_el(r), cx_don_stg(r))                     &
                    *dens(cx_acc_el(r), cx_acc_stg(r))
            if (n_pair .le. 0.0d0) cycle
            kc = min(max(cx_rate(r, T, n_e(j)), 0.0d0), cx_kc_max)
            q  = q + kc*n_pair*(cx_energy_defect_eV(r)                   &
                                 - cx_radiated_energy_eV(r, T)            &
                                 - product_manifold_radiation_eV(r, T,    &
                                        n_e(j), lvl_p, lvl_X, lvl_done))
         enddo
         ! He <-> H: the non-radiative He+ + H -> He(1^1S) + H+ ends in the
         ! only He state the 10.989 eV can reach (He 2^3S lies at 19.82 eV)
         ! and gives it all to the gas; the radiative channel B2 gives it to
         ! its photon; B1, the reverse of the non-radiative channel, takes
         ! it back. He2+ + H0 -> He+(1s) + H+ is radiative: its photon
         ! carries the whole 40.81 eV defect and the gas receives nothing
         ! (the kinetic share of a radiative transfer is not given by the
         ! source and is taken as zero, as for B2), so k_b3 adds no term.
         if (he_h_charge_exchange) then
            call he_h_cx_rates(T, k_b1, k_b2, k_b3)
            k_nr = min(nonradiative_electron_capture_Hep_from_H(T), k_b2)
            q = q + (k_nr*n_heii(j)*n_hi(j) - k_b1*n_hei(j)*n_hii(j))      &
                   *dE_He_H
         endif
         heat(j) = q/erg2eV
      enddo
      end subroutine charge_exchange_heating

      ! Energy defect [eV] of row r read forward, IP(acceptor's lower
      ! stage) - IP(donor): positive when the reaction releases energy.
      double precision function cx_energy_defect_eV(r) result(dE)
      integer, intent(in) :: r
      dE = cx_ionization_energy(cx_acc_el(r), cx_acc_stg(r)-1)            &
         - cx_ionization_energy(cx_don_el(r), cx_don_stg(r))
      end function cx_energy_defect_eV

      ! Energy [eV] that one event of row r at T [K] sends out as radiation:
      ! the excitation of the product state when that state decays by an
      ! allowed or intercombination line (it lies above the populated
      ! manifold of q_level_energy), or the photon of a radiative transfer.
      ! A product in the populated manifold (ground term or a collisionally
      ! populated metastable term) keeps its excitation in the gas, within
      ! the validity of that manifold (the note at q_level_energy).
      ! Term energies are NIST ASD levels weighted by 2J+1, the product
      ! J levels taken as statistically populated.
      !
      ! PRODUCT STATES KNOWN FROM THE SOURCES READ:
      !   A1  Mg + H+: Mg+ 3p 2Po (35730.36 cm^-1) and 4s 2S (69804.95),
      !       the channels into the 2 1Sigma and 4 1Sigma states of Allan,
      !       Clegg, Dickinson & Flower (1988), MNRAS 235, 1245, Table 1
      !       (read; "The Mg+ ions are formed in the 3p 2Po state and 4s 2S
      !       state, which cascades to 3p 2Po", Fig. 3 naming the 2 1Sigma
      !       curve the 3p channel). The 4s channel is endothermic by 2.70 eV
      !       and takes that from the gas. The fraction into 3p is their
      !       Table 1 ratio, interpolated in log T over 5e3 - 3.2e4 K and
      !       held at the ends.
      !   A9  Si + H+: Si+ 3s3p2 4P (43002.21 cm^-1), Kimura et al. (1996)
      !       "H+ + Si(3P) -> H + Si+(4P) + 1112 cm^-1"; 0.115 eV to the gas.
      !   A15 C + H+: the non-radiative part of reaction (2) of Stancil et
      !       al. (1998) ends in C+ 2Po and gives 2.338 eV; the rest of the
      !       tabulated (2) + (3) is radiative. The non-radiative part is
      !       the detailed-balance reverse of A16 (their reaction (1)), at
      !       most the whole of A15.
      !   C3  C + He+: C+ 2s2p2 2D (74931.09 cm^-1), Kimura et al. (1993),
      !       ApJ 417, 812 ("[He(1S) + C+(2D)]"); 4.037 eV to the gas.
      !   E1  O2+ + H: O+ 2s2p4 4P (119932.56 cm^-1), the "Decay 833 A"
      !       that Butler, Heil & Dalgarno (1980), ApJ 241, 442, Table 2
      !       (read) give for O+2 + H, the calculation Barragan et al.
      !       (2006), whose rate is carried, agree with below 1e4 K (their
      !       paper names no product state); 6.653 eV to the gas.
      !   E2  N2+ + H: N+ 2s2p3 3Do (the 2J+1 mean of the NIST levels
      !       92237.2, 92250.3, 92251.8 cm^-1), the term Barragan et
      !       al. (2006) name and whose 1084 A decay Butler, Heil & Dalgarno
      !       (1980) Table 2 list; it decays to the ground term by the
      !       allowed N II 1084-1086 A multiplet, so 16.003 - 11.437 =
      !       4.566 eV goes to the gas.
      !   A4, A12, A18, D25 and their detailed-balance reverses: the
      !       energy defect reaches no level above the populated manifold
      !       of the product (Mg+ 3p 4.43 eV > 1.44; Si+ 4P 5.31 > 2.75;
      !       N 2D 2.38 > 0.94; the notes at cx_reverse_of), so the product
      !       is in it by energetics.
      !   A19 S + H+ and its reverse: S+ 4S, 2D, 2P (Zhao et al. 2005), all
      !       in the populated manifold; the whole defect is heat.
      !   C1 Si + He+: Si+* Rydberg-like terms at about 1.19e5 - 1.32e5
      !       cm^-1 that radiate (Satta et al. 2013, Figs. 1 and 2, read
      !       from the page image; the crossing curves they draw lie at
      !       1.19e5 - 1.28e5); E_rad is taken at the middle, 1.24e5 cm^-1
      !       (15.37 eV), so 1.06 eV of the 16.44 eV defect is heat, within
      !       0.5 - 1.7 eV over the range drawn.
      !   D29 Na + Ca+: Ca 4s4p 3P (15263.09 cm^-1, 2J+1 mean, NIST), the
      !       exit channel Kwolek et al. (2019) infer, which radiates the
      !       657 nm intercombination line: the event takes 0.918 eV from
      !       the gas.
      !   A8 (2.60 eV) and A19 (3.24 eV): every level of Fe+ and S+ below
      !       the defect lies in the populated manifold (Fe II up to
      !       4.77 eV, S II up to 9.84 eV), so the whole defect is heat
      !       whatever the product level.
      !   A21 Na + H+: the non-radiative part ends in H(n=2) (Dutta et al.
      !       2001), whose 10.199 eV (NIST, 82259.11 cm^-1, the 2J+1 mean of
      !       2s and 2p) leaves as Lyman
      !       alpha, so it takes 10.199 - 8.459 = 1.740 eV from the gas (the
      !       paper writes the defect -1.682 eV); the radiative part (to
      !       H(1s) + photon, Watanabe et al. 2002) gives the gas nothing.
      !   A5 Fe + H+ and D17 Fe + O+ (defects 5.696 and 5.716 eV):
      !       Rutherford & Vroom (1972) could not determine the product
      !       state and attribute the large cross sections to near
      !       resonance (quoted at cx_A2_Mgp_H). The near-resonant Fe+
      !       terms are 3d6(5D)4p z4F and z4D (mean 44726.42 cm^-1 =
      !       5.545 eV, NIST, 2J+1 weighted), 0.15 - 0.17 eV below the
      !       defect, which is the kinetic energy released. Their
      !       spin-allowed one-electron decay 4p -> 4s leads to 3d6(5D)4s a4D
      !       (mean 8320.49 cm^-1 = 1.032 eV), a metastable term of the
      !       populated manifold whose energy returns to the gas, so
      !       E_rad = 5.545 - 1.032 = 4.514 eV and the heat is 1.18 eV
      !       (A5) and 1.20 eV (D17) (DERIVED; the decay branching is taken
      !       from the selection rules, not from a calculation read). The
      !       limits are 0.15 eV (no energy returned) and the full defect
      !       (ground-state product, which the large cross sections argue
      !       against). Pequignot & Aldrovandi (1986), A&A 161, 169, Table 3
      !       (read) list H+ + Fe among the reactions for which "several
      !       energetically favorable channels exist" within their window
      !       of energy defects, -0.2 to 1.1 eV (their Table 3 values "are
      !       nothing but guesses"): an independent statement of the same
      !       near-resonance principle, not a measurement of the product.
      !   RADIATIVE, E_rad = dE: C5 (O + He+, Zhao et al. 2004), A23 (K +
      !       H+, Watanabe et al. 2002) and the radiative He+ + H (B2). The
      !       photon carries the defect; the kinetic-energy change dE - h nu
      !       of a radiative transfer is not given by those sources and is
      !       taken as zero.
      ! ASSUMED IN THE POPULATED MANIFOLD (E_rad = 0) because the source
      ! names no product: A13/A14 (O, 0.02 eV), D1/D2 (source unknown),
      ! D27 (N+ + Fe, 6.63 eV: Rutherford & Vroom's near resonance puts the
      ! product in odd Fe+ terms near 6.6 eV whose decays, to terms not
      ! identified here, radiate most of it, so this is an upper bound),
      ! the Prasad & Huntress (1980) D-group forwards and D31.
      double precision function cx_radiated_energy_eV(r, T) result(e)
      integer, intent(in) :: r
      real*8,  intent(in) :: T
      real*8, parameter :: cm_to_eV = hp_erg*c_light*erg2eV
      real*8, parameter :: E_Mgp_3p = 35730.36d0*cm_to_eV
      real*8, parameter :: E_Mgp_4s = 69804.95d0*cm_to_eV
      real*8, parameter :: E_Sip_4P = 43002.21d0*cm_to_eV
      real*8, parameter :: E_Cp_2D  = 74931.09d0*cm_to_eV
      real*8, parameter :: E_Op_4P  = 119932.56d0*cm_to_eV
      real*8, parameter :: E_H_n2   = 82259.11d0*cm_to_eV
      real*8, parameter :: E_Fep_z4 = (44726.42d0 - 8320.49d0)*cm_to_eV
      real*8, parameter :: E_Sip_Ryd = 1.24d5*cm_to_eV
      real*8, parameter :: E_Ca_3P  = 15263.09d0*cm_to_eV
      real*8, parameter :: E_Np_3Do = (7.0d0*92237.2d0 + 5.0d0*92250.3d0 &
                                     + 3.0d0*92251.8d0)/15.0d0*cm_to_eV
      real*8 :: f3p, ratio, k_nr, k15, k_na
      select case (r)
      case (1)
         f3p = magnesium_3p_fraction(T)
         e   = f3p*E_Mgp_3p + (1.0d0 - f3p)*E_Mgp_4s
      case (5, 48)
         e = E_Fep_z4
      case (9)
         e = E_Sip_4P
      case (21)
         k_na = sodium_charge_transfer_ionization(T)
         if (k_na .gt. 0.0d0) then
            e = (sodium_nonradiative_ionization(T)*E_H_n2                 &
                 + sodium_radiative_charge_transfer(T)                     &
                   *cx_energy_defect_eV(21))/k_na
         else
            e = cx_energy_defect_eV(21)
         endif
      case (23)
         e = cx_energy_defect_eV(23)
      case (26)
         e = E_Sip_Ryd
      case (60)
         e = E_Ca_3P
      case (15)
         k15   = published_rate_coefficient(15, T)
         ratio = detailed_balance_ratio(15, T)
         k_nr  = 0.0d0
         if (ratio .gt. 0.0d0 .and. k15 .gt. 0.0d0)                        &
            k_nr = min(published_rate_coefficient(16, T)/ratio, k15)
         if (k15 .gt. 0.0d0) then
            e = cx_energy_defect_eV(15)*(1.0d0 - k_nr/k15)
         else
            e = cx_energy_defect_eV(15)
         endif
      case (28)
         e = E_Cp_2D
      case (30)
         e = cx_energy_defect_eV(30)
      case (64)
         e = E_Op_4P
      case (65)
         e = E_Np_3Do
      case default
         e = 0.0d0
      end select
      end function cx_radiated_energy_eV

      ! Energy [eV] that one event of row r at (T, n_e) loses to forbidden
      ! lines because the populated manifold of a product, or of a
      ! reactant, is not held at its equilibrium populations. The cooling
      ! of this code charges every ion of a species with the radiation of
      ! the equilibrium distribution p of metastable_level_populations; an
      ! ion that a reaction places in level c instead radiates X_c more
      ! (metastable_excess_radiation) before it has forgotten where it
      ! started, and an ion that a reaction removes from level c takes
      ! X_c of the charged radiation with it. So
      !   E_manifold = sum_c b_c X_c (products)  -  sum_c w_c X_c (reactants),
      ! over the species that have a level model (N I, S II, C I;
      ! metastable_model_of), with b_c the fraction of events that leave the
      ! product in level c and w_c the fraction that take the reactant from
      ! level c. Limits: X_c -> 0 at the collisional limit, where the
      ! populated-manifold heat of cx_radiated_energy_eV is exact, and
      ! X_c -> E_c as n_e -> 0, where the excitation of a metastable product
      ! leaves as its forbidden line (S+ 2P formed by S + H+ radiates 3.04 of
      ! the 3.24 eV defect).
      ! The product branching b_c is the source's where it names the
      ! product levels (A16 C(3P), A18 N(4S), A19 S+(4S, 2D, 2P) of Zhao et
      ! al. 2005 Table III, D25 C(3P) + S+(4S)), with the J levels of a term
      ! populated in proportion to g_J exp(-E_J/kT), and thermal over the
      ! five levels otherwise, which is the branching the reverses assume
      ! too (the note at q_level_energy). The reactant weights follow from
      ! the same rule: a forward reacts from every level with the one rate
      ! the source gives, so w = p and sum_c p_c X_c = 0 (the equilibrium
      ! ion radiates no excess); a derived reverse of a forward that names
      ! its product levels reacts from level c in proportion to
      ! beta_c b_c (the note at q_level_energy), so w_c = beta_c b_c / F.
      ! The internal energy an ion borrows from the gas when it is formed
      ! and returns when it is destroyed (the mean of E_c over p) is not a
      ! heat of either reaction and is not counted, as everywhere in this
      ! ledger. Products and reactants without a level model keep the
      ! populated-manifold heat, valid above the critical densities of their
      ! metastable terms (n_e >~ 1e7 cm^-3 at 1e4 K for Fe II, Fe III, O I,
      ! O II, N II, Si I, S I; lower for most).
      ! lvl_p, lvl_X and lvl_done hold, for the cell being evaluated, the
      ! equilibrium populations and excess radiation of each level model,
      ! formed on first use (a cell reads each model for several rows).
      double precision function product_manifold_radiation_eV(r, T, n_e, &
                                lvl_p, lvl_X, lvl_done) result(e)
      integer, intent(in)    :: r
      real*8,  intent(in)    :: T, n_e
      real*8,  intent(inout) :: lvl_p(5,3), lvl_X(5,3)
      logical, intent(inout) :: lvl_done(3)
      integer :: side, el, stg, isp, f, rb, c
      real*8  :: X(5), b(5), p(5), wb(5), w(5)
      e = 0.0d0
      ! the products of r, formed with the branching of r
      do side = 1, 2
         call product_species_of_row(r, side, el, stg)
         if (el .gt. n_melem) cycle
         isp = metastable_model_of(el, stg)
         if (isp .eq. 0) cycle
         call product_level_branching(r, isp, T, b)
         call level_model_of_cell(isp)
         e = e + sum(b*lvl_X(:,isp))
      enddo
      ! the reactants of a derived reverse whose forward names its product
      ! levels: they are that forward's products, removed with weights
      ! beta_c b_c
      f = cx_reverse_of(r)
      if (f .eq. 18 .or. f .eq. 19 .or. f .eq. 56) then
         do side = 1, 2
            call product_species_of_row(f, side, el, stg)
            if (el .gt. n_melem) cycle
            isp = metastable_model_of(el, stg)
            if (isp .eq. 0) cycle
            rb = f
            call product_level_branching(rb, isp, T, b)
            call level_model_of_cell(isp)
            p = lvl_p(:,isp)
            X = lvl_X(:,isp)
            call metastable_boltzmann_weights(isp, T, wb)
            ! w_c = beta_c b_c, zero where the Boltzmann weight of the
            ! level underflows (metastable_term_departures)
            w = 0.0d0
            do c = 1, 5
               if (wb(c) .gt. 0.0d0) w(c) = b(c)*p(c)*sum(wb)/wb(c)
            enddo
            if (sum(w) .gt. 0.0d0) then
               w = w/sum(w)
               e = e - sum(w*X)
            endif
         enddo
      endif

      contains

         subroutine level_model_of_cell(i)
         integer, intent(in) :: i
         if (.not. lvl_done(i)) then
            call metastable_excess_radiation(i, T, n_e, lvl_X(:,i),       &
                                             lvl_p(:,i))
            lvl_done(i) = .true.
         endif
         end subroutine level_model_of_cell

      end function product_manifold_radiation_eV

      ! Product species (element code, stage) of row r on one side: side 1
      ! the donor after the transfer, side 2 the acceptor after it.
      pure subroutine product_species_of_row(r, side, el, stg)
      integer, intent(in)  :: r, side
      integer, intent(out) :: el, stg
      if (side .eq. 1) then
         el = cx_don_el(r);  stg = cx_don_stg(r) + 1
      else
         el = cx_acc_el(r);  stg = cx_acc_stg(r) - 1
      endif
      end subroutine product_species_of_row

      ! Fraction b_c of the events of row r that leave the product species
      ! isp in level c (the rule stated at product_manifold_radiation_eV).
      subroutine product_level_branching(r, isp, T, b)
      integer, intent(in)  :: r, isp
      real*8,  intent(in)  :: T
      real*8,  intent(out) :: b(5)
      real*8  :: wb(5), wt(3), bs(3)
      integer :: c
      call metastable_boltzmann_weights(isp, T, wb)
      wt = 0.0d0
      do c = 1, 5
         wt(lp_term(c,isp)) = wt(lp_term(c,isp)) + wb(c)
      enddo
      bs = wt/sum(wt)
      select case (r)
      case (16, 18, 56)
         bs = [1.0d0, 0.0d0, 0.0d0]
      case (19)
         call sulfur_product_term_branching(T, bs)
      end select
      b = 0.0d0
      do c = 1, 5
         if (wt(lp_term(c,isp)) .gt. 0.0d0)                                &
            b(c) = bs(lp_term(c,isp))*wb(c)/wt(lp_term(c,isp))
      enddo
      end subroutine product_level_branching

      ! Fraction of Mg + H+ -> Mg+ + H that ends in Mg+ 3p: Allan, Clegg,
      ! Dickinson & Flower (1988), MNRAS 235, 1245, Table 1 (read), the
      ! "into 2 1Sigma" (3p) column over the sum of it and the "into
      ! 4 1Sigma" (4s) column, units 1e-10 cm^3 s^-1, 5e3 - 3.2e4 K. The
      ! sum, not the printed "total", is used: they agree to the last digit
      ! except at 2.6e4 K, where the total is printed 7.954 against the sum
      ! 7.854 of 1.504 and 6.350, and the neighboring totals (7.237,
      ! 8.439) place it near the sum. Interpolated linearly in log T and
      ! held at the end points.
      pure double precision function magnesium_3p_fraction(T) result(f)
      real*8, intent(in) :: T
      integer, parameter :: n_a = 16
      real*8, parameter :: T_a(n_a) = [ 5.0d3, 6.0d3, 8.0d3, 1.0d4,      &
         1.2d4, 1.4d4, 1.5d4, 1.6d4, 1.8d4, 2.0d4, 2.2d4, 2.4d4, 2.6d4,   &
         2.8d4, 3.0d4, 3.2d4 ]
      real*8, parameter :: k_3p(n_a) = [ 0.2549d0, 0.3430d0, 0.5055d0,   &
         0.6419d0, 0.7582d0, 0.8632d0, 0.9138d0, 0.9641d0, 1.066d0,       &
         1.170d0, 1.278d0, 1.390d0, 1.504d0, 1.620d0, 1.736d0, 1.850d0 ]
      real*8, parameter :: k_4s(n_a) = [ 0.07057d0, 0.1856d0, 0.6125d0,  &
         1.236d0, 1.955d0, 2.696d0, 3.060d0, 3.415d0, 4.093d0, 4.723d0,   &
         5.307d0, 5.847d0, 6.350d0, 6.819d0, 7.259d0, 7.673d0 ]
      real*8  :: x, w, f1, f2
      integer :: i
      x = log(min(max(T, T_a(1)), T_a(n_a)))
      i = 1
      do while (i .lt. n_a - 1 .and. x .gt. log(T_a(i+1)))
         i = i + 1
      enddo
      w  = (x - log(T_a(i)))/(log(T_a(i+1)) - log(T_a(i)))
      f1 = k_3p(i)/(k_3p(i) + k_4s(i))
      f2 = k_3p(i+1)/(k_3p(i+1) + k_4s(i+1))
      f  = (1.0d0 - w)*f1 + w*f2
      end function magnesium_3p_fraction

      end module charge_exchange

   module h2_vibrational_relaxation
   ! What becomes of the vibrational energy an H2 molecule is left with, and
   ! how much energy that is per event. Two drivers put energy into the same
   ! place and the constants belong in one file so they cannot drift apart:
   !
   !   * a fast photoelectron, which excites v = 1 and v = 2 directly and
   !     also excites the B and C states, most of which fluoresce back into
   !     bound, vibrationally excited levels of the ground state
   !     (electron_energy_degradation.f90);
   !   * a stellar Lyman or Werner band photon, which excites the same B and
   !     C states with the same outcome (lyman_werner.f90).
   !
   ! The vibrationally excited molecule then either radiates the energy away
   ! in the infrared quadrupole lines, or is collisionally de-excited and the
   ! energy becomes heat. Hollenbach & McKee (1979, ApJS 41, 555), p. 586:
   ! "The term (1 + n_cr/n)^-1 expresses the fraction of the excitation
   ! energy going into heat via collisional de-excitation rather than into
   ! radiation." Burton, Hollenbach & Tielens (1990, ApJ 365, 620), p. 623,
   ! say the same thing from the other side: "when n >> n_crit and
   ! n/G_0 >> 1, collisional deexcitation dominates and nearly every UV pump
   ! into excited v states (typically ~2 eV above the ground state) leads to
   ! the transformation of the vibrational energy into heat."
   !
   ! ---------------------------------------------------------------------
   !  VALIDITY
   ! ---------------------------------------------------------------------
   !  * COLLISION PARTNERS. Three are carried, and none of them is a fit any
   !    longer: atomic hydrogen from Lique (2015), helium from Jozwiak et
   !    al. (2024) and H2 from Le Bourlot, Pineau des Forets & Flower
   !    (1999), each reduced to the thermal v = 1 -> v' = 0 relaxation and
   !    tabulated at its own coefficient below. Helium used to be absent for
   !    want of a rate and is now present with its own.
   !    Electrons are left out: their rate coefficient is about two
   !    orders above the H2 one of Hollenbach & McKee (Burton et al.
   !    Table 7), but the molecular layer runs at an electron fraction of
   !    1e-16 to 5e-5, so they contribute a part in 1e3 or less of the
   !    collider sum there.
   !  * THE ORDERING OF THE THREE, MEASURED AT 1000 K: atomic hydrogen
   !    7.51e-12, helium 7.01e-15, H2 2.63e-15 cm^3 s^-1. Atomic hydrogen is
   !    three orders above the other two, its exchange channel being what
   !    makes it efficient, and H2 is a poor quencher of H2. The Hollenbach
   !    & McKee (1979) eq. (6.29) fits this module carried before had the
   !    atomic rate 1.14 to 3.47 times too large over 200-1000 K and the H2
   !    rate 27 to 500 times too large over 300-1527 K, i.e. they overstated
   !    the thermalized share, the H2 one grossly.
   !  * TEMPERATURE. The three tables cover 100-5000 K (H), 20-8000 K (He)
   !    and 100-6000 K (H2), each held at its end values rather than
   !    extrapolated. The molecular layer runs 200-1600 K, inside all three.
   !  * EVERY LEVEL OF THE CASCADE, NOT ONLY v = 1. The radiative rate in
   !    n_cr is the LARGEST total spontaneous decay rate over the whole
   !    bound X ladder, reduced at initialization from the same line list
   !    the infrared cooling uses (Roueff et al. 2019, A&A 630, A58,
   !    table 2, through molecular_infrared_data): 302 levels, 1833
   !    electric-quadrupole and magnetic-dipole lines. So the RADIATIVE side
   !    of the branching is the whole ladder's and not v = 1's. THE
   !    COLLISIONAL SIDE IS STILL v = 1's, one coefficient per collider, so
   !    the fraction is a first-event branching model and not a bound over
   !    the ladder; what it approximates, and what brackets it, are at
   !    h2_vibrational_heat_fraction below.
   !    MEASURED on that ladder, the maximum is 5.594e-06 s^-1, 6.6 times
   !    the v = 1, J = 0 total of 8.532e-07 s^-1 (itself within 3 per cent
   !    of the 8.3e-07 s^-1 Hollenbach & McKee adopt). Burton et al. instead
   !    use a pseudolevel at v = 6 with A_(6->0) = 2e-7 s^-1; the all-level
   !    maximum is 28 times that, i.e. the most conservative of the three.
   !    CROSS-CHECKED AGAINST THE COMPLETE QUADRUPOLE SET. The tree's line
   !    list leaves fifteen near-dissociation levels with no downward
   !    transition, so the maximum it gives could in principle be the
   !    maximum of an incomplete network. It is not: the complete set of
   !    Wolniewicz, Simbotin & Dalgarno (1998), ApJS 115, 293, which
   !    connects ALL bound rovibrational levels (4661 lines, read from the
   !    copy the Cloudy c25.00 distribution carries) gives 5.5945e-06 s^-1,
   !    agreeing to 1.1e-04 relative, and locates it at the same level,
   !    (v = 1, J = 28) -- a high-J rotational level of low v, not a
   !    near-dissociation one. Against the complete set only three levels of
   !    the code's ladder have no downward transition at all, (0,0), (0,1)
   !    and (14,4), and none of them can raise a maximum.
   !
   !  The distinction only matters where n is comparable to n_cr. Carrying
   !  the all-level maximum through the collisional coefficients below,
   !  n_cr is 1.27e6 cm^-3 at 808 K and 2.97e5 cm^-3 at 1527 K against the
   !  atomic hydrogen that dominates this layer's collider sum. MEASURED on
   !  the certified LHS 1140 b molecular base, that sum is 99.80 per cent
   !  atomic hydrogen, 0.20 per cent helium and 0.006 per cent H2 at cell 1,
   !  and 1 - f is 6.8e-07 there and 3.8e-07 at cell 25. It is
   !  written out rather than assumed because a shallower base, or the top
   !  of a thinner layer, would not be: against an H2 collider ALONE n_cr is
   !  2.1e9 cm^-3 at 1000 K, and at 1527 K with n(H2) = 4e9 cm^-3 and no
   !  atomic hydrogen the fraction would be 0.947, not 1. That is why the H2
   !  coefficient has to be the calculation and not the fit, which would
   !  have claimed 0.999 there.
   !
   !  Energies in eV where marked, otherwise cgs.

   ! The H2 line list, the same table the infrared cooling sums: the
   ! all-level radiative rate is reduced from it here rather than written
   ! down, so the two cannot describe different molecules. Only the lines
   ! are read, never the level ladder, which has one consumer
   ! (h2_partition_function of caloric_eos).
   use molecular_infrared_data, only: n_h2_line, h2_line_A, h2_line_Tu

   implicit none
   private

   public :: h2_energy_per_bound_fluorescence_eV,                         &
             h2_energy_per_bound_fluorescence_erg
   public :: e_vib_v1_eV, e_vib_v2_eV, e_vib_v5_eV, e_vib_v6_eV
   public :: h2_vibrational_heat_fraction, h2_vibrational_critical_density
   public :: h2_vibrational_relaxation_init, h2_total_decay_rate_max,     &
             h2_vibrational_relaxation_ready
   ! The three collider coefficients, for the assertions that hold each
   ! table to the value its own paper publishes.
   public :: gamma_10_H, gamma_10_H2, gamma_10_He

   ! ----- Energy left in the ground state per bound fluorescence ----- !
   ! An H2 molecule lifted into B or C either dissociates into the X
   ! continuum or fluoresces back into a bound level X(v'', J''), and the
   ! internal energy of that landing level is what collisions can turn into
   ! heat. The quantity below is the MEAN INTERNAL ENERGY OF THE LANDING
   ! LEVEL, i.e. per excitation that ends BOUND. It carries no branching
   ! ratio of its own: the caller multiplies by however many of its
   ! excitations end bound, which is (1 - p_diss) with the p_diss of the
   ! agent that did the exciting -- the photon branching for a Lyman-Werner
   ! pump, the state-resolved electron-impact branching for a photoelectron.
   !
   ! WHAT IT REPLACES. Burton, Hollenbach & Tielens (1990), Appendix A,
   ! p. 635: "The process is approximated as per Paper I, except that the
   ! effective heating per pump, E_*, has been taken as 2.0 eV." That is an
   ! ADOPTED value -- they revise it from their Paper I without showing a
   ! derivation, and Black & van Dishoeck (1987) pp. 412-413, the paper
   ! usually cited nearby, supply only the branching (a dissociation
   ! probability of 0.10, so 90% of the fluorescent transitions "populate
   ! various bound excited vibration-rotation levels of the ground state"),
   ! not a mean energy.
   !
   ! IT IS NOW COMPUTED. Summing the Abgrall, Roueff & Drira (2000)
   ! transition probabilities over every decay of every Lyman-Werner
   ! transition that starts on a thermally populated X(v, J) inside
   ! 912-1110 A, weighted by the absorption oscillator strength of the line
   ! and by the X(v'', J'') energy the decay lands on, gives
   !
   !     E_bound = 2.060 eV at  700 K   ->   2.128 eV at 3200 K
   !
   ! i.e. Burton et al.'s adopted 2.0 eV is 3 to 6 per cent low. Per PUMP
   ! rather than per bound fluorescence the same sum is 1.769 to 1.799 eV,
   ! which is the number to compare with 2.0 x (1 - p_diss); the decays that
   ! land back on v'' = 0 and leave only rotational energy are 16.3 to 16.5
   ! per cent of the bound fluorescences and are ALREADY INSIDE these means,
   ! so they must not be discounted a second time.
   ! Source and arithmetic: md/p39_lw_cross_section_sources.md sec. 5.2.
   !
   ! VALIDITY.
   !  * 700-3200 K, the grid of the level-resolved calculation; the fit is
   !    clamped outside it rather than extrapolated. Over that range a
   !    quadratic reproduces the nine computed points to 0.03 per cent (a
   !    straight line would do to 0.31 per cent).
   !  * UNSHIELDED WEIGHTS. The average is over the line mix an unattenuated
   !    band pumps. Deeper in the layer self-shielding removes the strongest
   !    lines first and the surviving mix lands slightly higher: the same sum
   !    with shielded weights gives 2.17-2.25 eV at N(H2) = 1e18-1e21 cm^-2,
   !    +4 to +9 per cent. The face value is used, so this term is a lower
   !    bound inside the shielded layer, by about as much as the correction
   !    it applies.
   !  * ELECTRON IMPACT SHARES IT, as the adopted constant did. The upper
   !    states and the fluorescence cascade are the same, but the
   !    distribution over the upper vibrational levels a photoelectron
   !    produces is not the one a band photon produces, so using the same
   !    mean landing energy there is an assumption and is marked as one.

   ! ----- Vibrational term energies of the X state ----- !
   ! Roueff et al. (2019), A&A 630, A58, table 2 -- the same line list
   ! molecular_infrared_data.f90 carries: (v=1,J=0) lies 5987.000 K and
   ! (v=2,J=0) 11635.400 K above (v=0,J=0). NOT the harmonic 2 x 0.516 eV:
   ! anharmonicity puts v = 2 at 1.0027 eV, 2.8% below twice v = 1.
   real*8, parameter :: e_vib_v1_eV = 0.515920d0
   real*8, parameter :: e_vib_v2_eV = 1.002661d0
   ! The two levels the H3+ dissociative-recombination product distribution
   ! peaks between (R6 of molecular_reaction_heat), from the same table:
   ! (v=5,J=0) lies 26605.762 K and (v=6,J=0) 30942.019 K above (v=0,J=0).
   real*8, parameter :: e_vib_v5_eV = 2.292710d0
   real*8, parameter :: e_vib_v6_eV = 2.666380d0

   ! ----- Radiative and collisional rates of v = 1 ----- !
   ! Hollenbach & McKee (1979), p. 583, section VI b: "A_10 = 8.3 x 10^-7
   ! s^-1, ... taken from Turner, Kirby-Docken, and Dalgarno (1977) as
   ! average A values for a J level in the given v states". This is the
   ! J-averaged TOTAL decay rate of v = 1, not the rate of one quadrupole
   ! line (those have "A values generally less than 3 x 10^-7 s^-1", same
   ! page and paragraph).
   ! Checked against the original: Turner et al.'s tables give the v = 1 ->
   ! v = 0 quadrupole branch sum as 8.54e-7 s^-1 at J = 0, 8.52e-7 at J = 1,
   ! 8.46e-7 at J = 2, 8.34e-7 at J = 3, 8.13e-7 at J = 4 and 6.10e-7 at
   ! J = 8, so Hollenbach & McKee's 8.3e-7 is that branch-summed rate at the
   ! low J levels that carry the population. Read the "generally less than
   ! 3 x 10^-7 s^-1" as a loose description of typical individual lines and
   ! not as a bound: the branch sum at J = 0 is 8.54e-7 s^-1, well above it.
   real*8, parameter :: a10_h2 = 8.3d-7

   ! ----- The all-level radiative rate of the cascade ----- !
   ! The largest total spontaneous decay rate over the 302 bound
   ! rovibrational levels of X^1 Sigma_g^+ [s^-1], reduced from the Roueff
   ! et al. (2019) line list by h2_vibrational_relaxation_init. This is the
   ! rate n_cr is built from, so that the fraction states the condition of
   ! the EVERY-LEVEL cascade and not of v = 1 alone (module header).
   !
   ! Until that initializer has run the value is Hollenbach & McKee's
   ! adopted v = 1 rate, which is what the fraction reduces to and which is
   ! the published number the consumers used before the reduction existed;
   ! it is never a sentinel, so a caller that has skipped the prologue gets
   ! the older published bound and not a meaningless one.
   ! h2_thermochemistry_init (mol_rates) calls the initializer in the serial
   ! prologue every molecular configuration passes through, before the first
   ! parallel region opens.
   real*8,  save :: a_cascade_max     = a10_h2
   logical, save :: a_cascade_ready   = .false.

   ! ----- Collisional de-excitation by H2 ----- !
   ! Le Bourlot, Pineau des Forets & Flower (1999), MNRAS 305, 802,
   ! state-to-state H2 + H2 rate coefficients, read from the tables the
   ! Cloudy c25.00 distribution carries (coll_rates_H2para_LeBourlot.dat,
   ! 627 de-excitation rows over 51 levels to (v = 3, j = 8), on the grid
   ! T = 100, 300, 500, 1000, 1500, 2000, 3000, 4500, 6000 K).
   !
   ! WHICH OF THE TWO DISTRIBUTED H2 FILES, AND WHY ONLY ONE.  Cloudy
   ! carries two, named for an ortho-H2 and for a para-H2 PERTURBER
   ! (parse_atom_h2.cpp fills two separate collider slots from them).  THE
   ! PAPER ITSELF SAYS THERE IS ONLY ONE SET.  Section 2.1, published text:
   ! Flower & Roueff (1998a) computed the rovibrational excitation of
   ! ortho- and para-H2 "by para-H2 in its rotational ground state
   ! (J = 0)"; test calculations with ground-state ortho-H2 as the
   ! perturber "showed that the rate of v = 1 -> 0 vibrational relaxation
   ! was insensitive to the rotational state of the perturber"; and "we
   ! have, therefore, adopted the rate coefficients calculated by Flower &
   ! Roueff (1998a) for the excitation of ortho- and para-H2 by para-H2
   ! (J = 0), and applied them also to the case of excitation by ortho-H2".
   ! The para file is that set: it carries the level list the paper states,
   ! "A total of 51 rovibrational energy levels ... (v = 0, J <= 16; v = 1,
   ! J <= 13; v = 2, J <= 10; v = 3, J <= 8)", while the file named for the
   ! ortho perturber holds 39 levels reaching only (1,13) and its rates
   ! stand a median 1.74 above the para file's on the 361 rows they share
   ! (MEASURED at 1000 K).  A perturber-resolved pair is therefore not what
   ! the paper published, and the para file alone is used for an H2
   ! collider of either form.
   !
   ! CHECKED AGAINST THE PAPER'S OWN TABLE, which is the only place it
   ! prints numbers: Table 1 gives the critical density of the 1-0 S(1)
   ! transition, n_cr = A/q, for each of the three perturbers at 500, 1000
   ! and 2000 K.  With A(1-0 S(1)) = 3.470e-07 s^-1 from the tree's own line
   ! list and q the (v = 1, j = 3) -> (v = 0, j = 1) row of the file,
   ! MEASURED: the para file gives 2.33e+11, 1.65e+10 and 5.80e+08 cm^-3
   ! against the published 3.7e+11, 1.8e+10 and 4.2e+08, i.e. ratios 1.59,
   ! 1.09 and 0.72, scattered about one; the file named for the ortho
   ! perturber gives 1.13e+11, 7.5e+09 and 2.94e+08, low by 3.3, 2.4 and
   ! 1.4 throughout.  The SAME procedure on the helium file of the same
   ! paper returns 5.59e+10, 1.56e+09 and 5.50e+07 against the published
   ! 4.2e+10, 1.5e+09 and 6.8e+07, ratios 0.75, 0.96 and 1.24, which is what
   ! establishes that the procedure and not the data is being tested.
   !
   ! The published paper is references/LeBourlot_1999_MNRAS_305_802.pdf;
   ! what the file itself was checked for is that every row is a
   ! DE-EXCITATION, upper (v', j') to lower (v, j), which is what its header
   ! states and what the low-temperature behavior confirms (a rate that is
   ! finite at 100 K has no threshold). The file distributes only that
   ! direction, so no detailed-balance identity can be formed from it; the
   ! helium set does distribute both.
   !
   ! WHAT IS TABULATED is the same quantity as for the other two colliders:
   ! the thermal v = 1 -> v' = 0 relaxation, summed over final j' and
   ! Boltzmann averaged over the initial j of v = 1 with the weights
   ! g_I (2j+1) exp(-E(1,j)/kT) of the code's own ladder.
   !
   ! WHAT IT CHANGES, AND IT IS LARGE.  Hollenbach & McKee's (1979)
   ! eq. (6.29) fit, gamma_10^H2 = 1.4e-12 T^1/2 exp[-12000/(T + 1200)],
   ! which this module carried before, is 500 times this calculation at
   ! 300 K, 147 at 808 K, 72 at 1000 K and 27 at 1527 K (MEASURED). This
   ! also settles the discrepancy the module used to record as unchecked:
   ! Burton, Hollenbach & Tielens (1990) Table 7 reprint the same
   ! coefficient with 18100 in the exponent where eq. (6.29) has 12000, and
   ! against the calculation that reading is 8.57, 7.04, 4.50 and 2.89 times
   ! too large at the same four temperatures (MEASURED). Burton et al. are
   ! therefore the closer of the two readings and both are too large; the
   ! question of which fit is right is retired, because neither is used.
   ! The fit is a 1979 V-T estimate and the calculation is a scattering
   ! result; H2
   ! is in fact a POOR vibrational quencher of H2, comparable to helium
   ! (7.0e-15 against 2.6e-15 cm^3 s^-1 at 1000 K) and three orders below
   ! atomic hydrogen, whose exchange channel is what makes it efficient.
   ! The ordering gamma_H >> gamma_H2 ~ gamma_He is therefore the one the
   ! three sources give, and the fit had the H2 collider three orders too
   ! close to the atomic one.
   !
   ! WHERE IT MATTERS.  Not in this layer: at the certified base the
   ! collider sum is 99.7 per cent atomic hydrogen, so the fraction is
   ! unchanged. It matters in a cold H2 base with no atomic hydrogen, where
   ! this coefficient alone sets n_cr, and there the fit was claiming
   ! fifty times too much thermalization.
   !
   ! RANGE. 100-6000 K, with the end values held outside rather than
   ! extrapolated; the molecular layer runs 200-1600 K, inside.
   !
   ! WHAT IS NOT USED, AND WHY. The Cloudy distribution also carries
   ! coll_rates_H2ortho_ORNL.dat and coll_rates_H2para_ORNL.dat (Wan et al.
   ! 2018, ApJ 862, 132), a newer H2-H2 set on a finer temperature grid.
   ! MEASURED: all 240 of their rows are PURELY ROTATIONAL, v = 0 to v = 0,
   ! so that set contains no vibrational transition and cannot give this
   ! coefficient.
   integer, parameter :: n_g10_h2 = 9
   real*8,  parameter :: t_g10_h2(n_g10_h2) = (/                         &
        100.0d0, 300.0d0, 500.0d0, 1000.0d0, 1500.0d0, 2000.0d0,         &
        3000.0d0, 4500.0d0, 6000.0d0 /)
   real*8,  parameter :: g10_h2_tab(n_g10_h2) = (/                       &
        5.9364d-18, 1.6258d-17, 7.9918d-17, 2.6310d-15, 2.3015d-14,      &
        9.5076d-14, 5.1397d-13, 1.8230d-12, 3.5736d-12 /)

   ! ----- Collisional de-excitation by atomic H, from the accurate
   !       quantum calculation rather than from a fit ----- !
   ! Lique (2015), MNRAS 453, 810, "Revisited study of the ro-vibrational
   ! excitation of H2 by H": nearly exact time-independent quantum scattering
   ! with the hydrogen exchange channels treated rigorously, state-to-state
   ! k(v, j -> v', j'; T) distributed as supplementary data to the article
   ! (references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat, 1431
   ! transitions, T = 100 to 5000 K in steps of 100 K, READ).
   !
   ! WHAT IS TABULATED BELOW is the THERMAL v = 1 -> v' = 0 vibrational
   ! relaxation: for each temperature of the file's own grid, the rates out
   ! of the fifteen initial levels j = 0 to 14 of v = 1 are summed over
   ! every final j' of v' = 0 and averaged over the initial j with the
   ! Boltzmann weights g_I (2j+1) exp(-E(1,j)/kT) of the code's own ladder.
   ! That is the quantity Lique's Table 2 quotes and it validates against
   ! it: 1.781e-13 cm^3 s^-1 at 300 K here (MEASURED) against his 1.8e-13,
   ! a measured 3.0 +/- 1.5e-13 (Heidner & Kasper 1972), 0.7e-13 for
   ! Wrathmall et al. (2007) and 0.8e-13 for Martin & Mandy (1995).
   !
   ! WHY IT REPLACES THE FIT.  Hollenbach & McKee's eq. (6.29),
   ! gamma_10^H = 1.0e-12 T^1/2 exp(-1000/T), is 1.14 times this at 200 K,
   ! 3.47 at 300 K, 1.89 at 800 K, 1.55 at 1000 K and 1.05 at 1600 K, and
   ! crosses 1 near 1750 K (MEASURED).  It therefore overstated the
   ! de-excitation, and so the heated share, everywhere in this layer.
   ! Atomic hydrogen is not a detail here: at the certified base (808 K,
   ! x2 = 0.355, He/H = 2.13) the collider sum below is 99.7 per cent
   ! atomic hydrogen and 0.3 per cent H2 (MEASURED), so this coefficient is
   ! what sets n_cr numerically.
   !
   ! RANGE AND WHAT IS OUTSIDE IT.  100-5000 K, with the end values held
   ! outside rather than extrapolated; the molecular layer runs 200-1600 K,
   ! inside.  The calculation covers internal energies below 22000 K, which
   ! is 55 of the code's 302 bound levels (v = 0 to 4), so it gives the
   ! v = 1 rate this module uses and says nothing about the
   ! near-dissociation levels a nascent molecule is born in; using the v = 1
   ! coefficient for those is the same approximation as before and still
   ! understates their de-excitation, which errs towards radiating rather
   ! than heating.
   integer, parameter :: n_g10_h = 50
   real*8,  parameter :: t_g10_h_lo = 100.0d0, t_g10_h_step = 100.0d0
   real*8,  parameter :: g10_h_tab(n_g10_h) = (/                         &
        6.7952d-14, 8.3919d-14, 1.7809d-13, 4.8006d-13, 1.0502d-12,      &
        1.8889d-12, 2.9765d-12, 4.2899d-12, 5.8071d-12, 7.5097d-12,      &
        9.3787d-12, 1.1398d-11, 1.3554d-11, 1.5830d-11, 1.8215d-11,      &
        2.0696d-11, 2.3260d-11, 2.5897d-11, 2.8595d-11, 3.1346d-11,      &
        3.4138d-11, 3.6963d-11, 3.9813d-11, 4.2680d-11, 4.5557d-11,      &
        4.8435d-11, 5.1313d-11, 5.4179d-11, 5.7034d-11, 5.9872d-11,      &
        6.2686d-11, 6.5476d-11, 6.8238d-11, 7.0965d-11, 7.3659d-11,      &
        7.6323d-11, 7.8947d-11, 8.1529d-11, 8.4077d-11, 8.6582d-11,      &
        8.9044d-11, 9.1463d-11, 9.3840d-11, 9.6175d-11, 9.8469d-11,      &
        1.0071d-10, 1.0292d-10, 1.0508d-10, 1.0720d-10, 1.0928d-10 /)

   ! ----- Collisional de-excitation by helium ----- !
   ! Jozwiak, Thibault, Viel, Wcislo & Lique (2024), A&A 685, A113,
   ! "Revisiting the rovibrational (de-)excitation of molecular hydrogen by
   ! helium": quantum scattering on a state-of-the-art surface,
   ! state-to-state rate coefficients for 1059 transitions between
   ! rovibrational levels of H2 with internal energies up to 15000 cm^-1,
   ! T = 20 to 8000 K (READ, their abstract and the VizieR ReadMe of
   ! references/Jozwiak_2024J_A+A_685_A113/, 26 ortho and 27 para levels,
   ! 520 + 539 converged transitions on a 43-point temperature grid).
   !
   ! WHAT IS TABULATED is the same quantity as for atomic hydrogen, reduced
   ! the same way: the thermal v = 1 -> v' = 0 relaxation, summed over final
   ! j' within each nuclear-spin symmetry (helium does not convert ortho to
   ! para) and Boltzmann averaged over the initial j of v = 1.  MEASURED
   ! here, 2.654e-15 cm^3 s^-1 at 808 K and 7.015e-15 at 1000 K.
   !
   ! WHAT IT CHANGES.  Helium used to be left out of the collider sum for
   ! want of a rate.  It is now in it, with its own coefficient, and the
   ! statement that used to stand here -- that leaving it out understated
   ! the de-excitation by an unknown amount -- is replaced by the size of
   ! the term: helium is 9.3e-04 of the atomic-hydrogen coefficient at
   ! 1000 K and 5.9e-04 at 800 K (MEASURED), so at the base helium density
   ! of this configuration it carries about a part in a thousand of the
   ! sum.  It is carried because it is physically present, not because it
   ! moves a number here.
   !
   ! RANGE.  20-8000 K, wider than this layer at both ends, with the end
   ! values held outside rather than extrapolated.  The set reaches v <= 3,
   ! so like the atomic-hydrogen one it gives the v = 1 rate and not the
   ! near-dissociation levels.
   integer, parameter :: n_g10_he = 43
   real*8,  parameter :: t_g10_he(n_g10_he) = (/                         &
        20.0d0, 30.0d0, 40.0d0, 50.0d0, 60.0d0, 70.0d0, 80.0d0, 90.0d0,  &
        100.0d0, 120.0d0, 140.0d0, 160.0d0, 180.0d0, 200.0d0, 250.0d0,   &
        300.0d0, 350.0d0, 400.0d0, 450.0d0, 500.0d0, 550.0d0, 600.0d0,   &
        650.0d0, 700.0d0, 750.0d0, 800.0d0, 850.0d0, 900.0d0, 950.0d0,   &
        1000.0d0, 1100.0d0, 1200.0d0, 1300.0d0, 1400.0d0, 1500.0d0,      &
        1750.0d0, 2000.0d0, 3000.0d0, 4000.0d0, 5000.0d0, 6000.0d0,      &
        7000.0d0, 8000.0d0 /)
   real*8,  parameter :: g10_he_tab(n_g10_he) = (/                       &
        6.4077d-20, 8.6066d-20, 1.2051d-19, 1.6964d-19, 2.3360d-19,      &
        3.1271d-19, 4.0773d-19, 5.1973d-19, 6.5002d-19, 9.7228d-19,      &
        1.3918d-18, 1.9327d-18, 2.6295d-18, 3.5332d-18, 7.3118d-18,      &
        1.5542d-17, 3.3338d-17, 6.8729d-17, 1.3232d-16, 2.3659d-16,      &
        3.9512d-16, 6.2209d-16, 9.3222d-16, 1.3410d-15, 1.8651d-15,      &
        2.5238d-15, 3.3391d-15, 4.3375d-15, 5.5502d-15, 7.0146d-15,      &
        1.0881d-14, 1.6372d-14, 2.4044d-14, 3.4561d-14, 4.8682d-14,      &
        1.0527d-13, 2.0308d-13, 1.2074d-12, 3.3057d-12, 6.2460d-12,      &
        9.6318d-12, 1.3111d-11, 1.6433d-11 /)

   ! Quadratic fit to the nine computed points (see the block above).
   real*8, parameter :: eb_c0 =  2.0247811011d0
   real*8, parameter :: eb_c1 =  5.6035354525d-05
   real*8, parameter :: eb_c2 = -7.4141671520d-09
   real*8, parameter :: eb_tmin = 700.0d0, eb_tmax = 3200.0d0
   real*8, parameter :: ev2erg = 1.602176634d-12

   contains

   ! Reduce the line list to the largest TOTAL spontaneous decay rate over
   ! the bound ladder: for each level, the sum of the A values of every line
   ! that leaves it downward; the maximum of those sums is what n_cr is
   ! built from. Idempotent, and serial by construction: called from
   ! h2_thermochemistry_init (mol_rates), the prologue that runs before any
   ! parallel region opens.
   !
   ! Lines are grouped by the UPPER TERM ENERGY the list carries with each
   ! of them, which is the level's own energy: two lines leave the same
   ! level exactly when their h2_line_Tu agree. The level ladder itself is
   ! not read here, so the one Boltzmann sum over it stays the one place
   ! that touches it. The closest pair of levels in this molecule is
   ! 0.199 K apart (MEASURED on the same table), so the 0.05 K window below
   ! separates every pair and is far above the rounding of the table.
   !
   ! Seventeen levels have no downward line and contribute zero: (0,0) and
   ! (0,1), which cannot decay, and fifteen near-dissociation levels of
   ! v = 10 to v = 14 the list does not reach. Those fifteen would trap
   ! population in a level-resolved cascade; they cannot raise this maximum,
   ! which belongs to a high-J level of low v.
   subroutine h2_vibrational_relaxation_init
      integer :: il, jl
      real*8  :: a_tot
      real*8, parameter :: dT_same_level = 0.05d0   ! [K]
      if (a_cascade_ready) return
      a_cascade_max = 0.0d0
      do il = 1, n_h2_line
         a_tot = 0.0d0
         do jl = 1, n_h2_line
            if (abs(h2_line_Tu(jl) - h2_line_Tu(il))                     &
                .le. dT_same_level) a_tot = a_tot + h2_line_A(jl)
         enddo
         if (a_tot .gt. a_cascade_max) a_cascade_max = a_tot
      enddo
      a_cascade_ready = .true.
   end subroutine h2_vibrational_relaxation_init

   ! The reduced rate [s^-1], for the consumers that state the bound and for
   ! the assertion that holds it to the ladder.
   double precision function h2_total_decay_rate_max() result(a_max)
      a_max = a_cascade_max
   end function h2_total_decay_rate_max

   ! Whether the reduction above has run.
   logical function h2_vibrational_relaxation_ready() result(is_ready)
      is_ready = a_cascade_ready
   end function h2_vibrational_relaxation_ready

   ! Mean internal energy of the X(v'', J'') level a Lyman-Werner
   ! fluorescence lands on, at gas temperature T [K], in eV. Per excitation
   ! that ends BOUND: the caller supplies the bound fraction.
   elemental function h2_energy_per_bound_fluorescence_eV(T) result(e)
      real*8, intent(in) :: T
      real*8 :: e, x
      x = min(max(T, eb_tmin), eb_tmax)
      e = eb_c0 + x*(eb_c1 + x*eb_c2)
   end function h2_energy_per_bound_fluorescence_eV

   ! The same in erg.
   elemental function h2_energy_per_bound_fluorescence_erg(T) result(e)
      real*8, intent(in) :: T
      real*8 :: e
      e = ev2erg*h2_energy_per_bound_fluorescence_eV(T)
   end function h2_energy_per_bound_fluorescence_erg

   ! Rate coefficient for thermal collisional de-excitation of H2 v = 1 into
   ! v' = 0 by atomic hydrogen, at temperature T [K]: interpolation of
   ! log k, linear in T, on Lique's (2015) own 100 K grid, with the end
   ! values held at 100 and 5000 K rather than extrapolated (see the table
   ! above for the source, the reduction and the range).
   !
   ! THE INTERPOLATION IS OF log k, NOT OF k, for all three colliders. A
   ! vibrational de-excitation rate is thermally activated and rises by
   ! orders across a grid interval at the low end; interpolating k itself
   ! between 500 and 1000 K on the H2 table, whose grid is the coarsest of
   ! the three, reads 2.3 times high at 808 K (MEASURED).
   elemental function gamma_10_H(T) result(gamma_deex)
      real*8, intent(in) :: T
      real*8  :: gamma_deex, x, frac
      integer :: i
      x = (min(max(T, t_g10_h_lo),                                       &
               t_g10_h_lo + t_g10_h_step*dble(n_g10_h - 1))              &
           - t_g10_h_lo)/t_g10_h_step
      i = min(int(x) + 1, n_g10_h - 1)
      frac = x - dble(i - 1)
      gamma_deex = exp( log(g10_h_tab(i))                                &
                      + frac*(log(g10_h_tab(i + 1))                      &
                              - log(g10_h_tab(i))) )
   end function gamma_10_H

   ! The same quantity for a helium collider, on Jozwiak et al.'s (2024)
   ! own 43-point grid, which is not uniform, so the interval is found by
   ! search; log k is interpolated, and the end values are held outside
   ! 20-8000 K rather than extrapolated.
   elemental function gamma_10_He(T) result(gamma_deex)
      real*8, intent(in) :: T
      real*8  :: gamma_deex, tc, frac
      integer :: i, k
      tc = min(max(T, t_g10_he(1)), t_g10_he(n_g10_he))
      i  = 1
      do k = 1, n_g10_he - 1
         if (tc .ge. t_g10_he(k)) i = k
      enddo
      frac = (tc - t_g10_he(i))/(t_g10_he(i + 1) - t_g10_he(i))
      gamma_deex = exp( log(g10_he_tab(i))                               &
                      + frac*(log(g10_he_tab(i + 1))                     &
                              - log(g10_he_tab(i))) )
   end function gamma_10_He

   ! The same quantity for an H2 collider, on Le Bourlot et al.'s (1999)
   ! own nine-point grid, which is not uniform, so the interval is found by
   ! search; log k is interpolated, and the end values are held outside
   ! 100-6000 K rather than extrapolated.
   elemental function gamma_10_H2(T) result(gamma_deex)
      real*8, intent(in) :: T
      real*8  :: gamma_deex, tc, frac
      integer :: i, k
      tc = min(max(T, t_g10_h2(1)), t_g10_h2(n_g10_h2))
      i  = 1
      do k = 1, n_g10_h2 - 1
         if (tc .ge. t_g10_h2(k)) i = k
      enddo
      frac = (tc - t_g10_h2(i))/(t_g10_h2(i + 1) - t_g10_h2(i))
      gamma_deex = exp( log(g10_h2_tab(i))                               &
                      + frac*(log(g10_h2_tab(i + 1))                     &
                              - log(g10_h2_tab(i))) )
   end function gamma_10_H2

   ! Critical density against a given collider mix, n_cr = A/gamma_deex
   ! (Hollenbach & McKee 1979, p. 576) [cm^-3], with A the all-level maximum
   ! of the cascade. Diagnostic: the heat fraction below does not go through
   ! it, so that it stays finite when the collider density vanishes.
   elemental function h2_vibrational_critical_density(T, f_HI) result(ncr)
      real*8, intent(in) :: T      ! [K]
      real*8, intent(in) :: f_HI   ! atomic-H share of the collider density
      real*8 :: ncr, gamma_deex
      gamma_deex = f_HI*gamma_10_H(T) + (1.0d0 - f_HI)*gamma_10_H2(T)
      ncr = a_cascade_max/max(gamma_deex, 1.0d-99)
   end function h2_vibrational_critical_density

   ! Fraction of the vibrational excitation energy that becomes HEAT rather
   ! than being radiated away: Hollenbach & McKee's (1 + n_cr/n)^-1, written
   ! with the collider densities explicit so it is finite and correct at
   ! either end,
   !
   !     f = (gamma_H n_HI + gamma_H2 n_H2) / (gamma_H n_HI + gamma_H2 n_H2
   !                                           + A_max) .
   !
   ! with the helium term carried whenever the caller has a helium density,
   !
   !     f = (gamma_H n_HI + gamma_H2 n_H2 + gamma_He n_He)
   !         / (same + A_max) .
   !
   ! WHAT THIS FORM IS: A FIRST-EVENT BRANCHING MODEL, and it is not a
   ! bound. It asks what becomes of ONE excited molecule at its FIRST
   ! disposal, with A_max the largest total spontaneous decay rate over the
   ! whole bound ladder and the collisional side the v = 1 coefficient of
   ! each collider. Two things follow, and both are stated because the
   ! comment that stood here asserted the opposite:
   !
   !  * A FIRST-EVENT FRACTION IS NOT THE ENERGY FRACTION OF A CASCADE. A
   !    molecule born at v = 10 reaches v = 0 through many steps, each with
   !    its own branching, and the share of its energy that ends as heat is
   !    a product over the path and not this single ratio.
   !  * THE ALL-LEVEL RADIATIVE MAXIMUM DOES NOT MAKE THE RATIO A LOWER
   !    BOUND. A_max bounds every level's radiative loss from above, but the
   !    v = 1 collisional coefficient does not bound every level's
   !    collisional loss from below: a level with C_u < C_1 would have a
   !    collisional branching fraction below this one even though its
   !    A_u <= A_max. The argument from closely spaced high levels that used
   !    to be made here is not a proof -- selection rules, the collider's
   !    identity, rotational redistribution, reactive destruction of the
   !    excited molecule and upward thermal transitions all enter, and
   !    Lique's (2015) section 3.2, which reports a modest increase with
   !    initial vibrational state for a fixed vibrational change in the
   !    H-collision transitions he studied, is supporting evidence inside
   !    that data set and not a minimum for every level, rotational state
   !    and collider a nascent molecule reaches.
   !
   ! WHAT BRACKETS IT. A reduced statistical-equilibrium model solves the
   ! levels the published collision data cover, carries the levels above
   ! them as one explicitly uncertain group, and evaluates the exact form
   ! below on its solution.
   !
   ! WHERE THE FORM COMES FROM: Burton, Hollenbach & Tielens (1990)
   ! eq. (A1), whose structure and counting unit are set out in
   ! md/p39_lw_cross_section_sources.md section 5.2 and are not restated
   ! here; this function is their quenched share n gamma/(A + n gamma).
   !
   ! WHAT THIS FRACTION APPROXIMATES, written down so that the
   ! approximation has a target and not only a caveat. The exact statement
   ! of the same physics is the NET COLLISIONAL HEAT of the ladder,
   !
   !     Q = sum over level pairs (u, l) and colliders M of
   !         [ n_u C_ul(M) - n_l C_lu(M) ] n_M (E_u - E_l) ,
   !
   ! with the upward coefficients from the downward ones by detailed
   ! balance and the populations n_u from a statistical equilibrium of the
   ! whole X ladder. That form carries no n_cr, no single level standing in
   ! for the cascade and no fraction at all, and it is exact in both the
   ! collisionally starved and the thermalized limit. It is what Cloudy
   ! c25.00 evaluates over its 303-level model of the same state
   ! (md/cloudy_h2_model_reference_20260916.md section 5.1).
   !
   ! IT IS NOT BUILDABLE HERE, and the reason is data, not effort: the
   ! H2-H rate coefficients stop at v = 3 in everything that distributes
   ! numbers, while the molecules this fraction is applied to are born at
   ! v = 10 to 14 (R15) and v = 5 to 6 (R6). Until that gap closes there
   ! are no level populations to put in the expression above.
   !
   ! HOW FAR THIS LAYER IS FROM THE DISTINCTION MATTERING: 1 - f = 7.0e-09
   ! at the certified base with the collider sum this module now carries
   ! (MEASURED). That
   ! number is a property of THIS model, not a bound on its distance from a
   ! level-population one; what makes the two agree at this base is that
   ! every level of the ladder is collisionally dominated there by orders,
   ! so both forms return the thermalized limit whatever their branching
   ! details.
   !
   ! n_He IS OPTIONAL, AND WHY.  The argument is the interface's and not the
   ! callers': a caller that has no helium density to hand can omit it and
   ! still get the two-collider sum. EVERY PRODUCTION CALL SITE PASSES ONE.
   ! The chemical-heat ledger passes the ground-singlet neutral
   ! (molecular_reaction_heat), and so do the three photoelectric and
   ! Lyman-Werner sites: ionization_equilibrium.f90 line 1262 passes
   ! he_ground_singlet_density(nhei, nheiTR) into f_vib_quench, and
   ! util_ion_eq.f90 lines 2113 and 2385 pass the same ground-singlet
   ! density into the Lyman-Werner fluorescence channel of the heating
   ! assembly and of the band ledger (READ). The helium collider carries
   ! about a part in a thousand of the collider sum at this base
   ! (MEASURED) and nothing at all in f where f is 1 to a part in 1e6,
   ! which is why the reading was the same before those sites passed it.
   elemental function h2_vibrational_heat_fraction(T, n_HI, n_H2, n_He)   &
                      result(f)
      real*8, intent(in) :: T          ! [K]
      real*8, intent(in) :: n_HI, n_H2 ! [cm^-3]
      real*8, intent(in), optional :: n_He ! [cm^-3]
      real*8 :: f, q
      q = gamma_10_H(T)*max(n_HI, 0.0d0) + gamma_10_H2(T)*max(n_H2, 0.0d0)
      if (present(n_He)) q = q + gamma_10_He(T)*max(n_He, 0.0d0)
      f = q/(q + a_cascade_max)
   end function h2_vibrational_heat_fraction

   ! End of module
   end module h2_vibrational_relaxation

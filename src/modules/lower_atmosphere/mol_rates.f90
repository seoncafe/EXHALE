      module mol_rates
      ! molecular (H2 / H2+ / H3+ / HeH+) reaction-rate coefficients,
      ! transcribed VERBATIM from Koskinen et al. (2022), ApJ 929:52, Table 1
      ! (page 19; verified against the PDF in references/).  All rates in
      ! cgs (cm^3 s^-1; the two three-body rates R13/R15 in cm^6 s^-1 --
      ! multiply by the total density n as in the table's "n" factor).
      ! T = heavy-particle temperature [K]; Te = electron temperature [K]
      ! (EXHALE uses a common T).
      !
      ! ONE ENTRY IS NOT TRANSCRIBED: R12, the thermal dissociation of H2, is
      ! built by detailed balance from the R15 recombination instead of from
      ! the Table-1 fit, because the two Table-1 entries are the same channel
      ! in opposite directions and their source measurements do not overlap
      ! in temperature.  The full argument and its verification are at R12.
      !
      ! Photo-rates P1-P5 (H, He, H2 photoionization; H2 dissociative and
      ! double photoionization) are "SC" in the paper -- computed from cross
      ! sections x stellar flux x column densities.  In EXHALE these follow
      ! the existing PH_heat routines (util_ion_eq) once H2 cross sections
      ! are added; they are NOT part of this module.
      !
      ! THIRD BODY M (R12, R13, R15).  Koskinen et al. (2022) write the
      ! three-body rates as a two-body coefficient times "n" and do not say
      ! which particles that n counts; their models are H2/H-dominated
      ! Jovian and Neptunian envelopes, where M is H2 and H.  EXHALE passes
      ! the TOTAL gas-particle density (electrons excluded; ion_cell_state
      ! %ntot from calc_ntot), so every heavy particle is treated as an
      ! equally efficient third body.  A monatomic third body has no internal
      ! modes to take up the released energy, so He is in general a less
      ! efficient third body than H2, and R12/R13/R15 with He as M are then
      ! upper bounds.  No He-specific efficiency factor is applied because
      ! none of the sources quoted for these three reactions supplies one;
      ! the caveat is recorded, with the size of the affected terms
      ! measured, in docs/molecular_chemistry_audit_he_rich.md.
      !
      ! HE-DOMINATED LIMIT.  Two further limits of this network are recorded
      ! there rather than patched here, because both would replace a
      ! Koskinen Table-1 entry with a rate from another compilation:
      !   * HeH+ formation. Table 1 forms HeH+ only through R20
      !     (He+ + H2, 4.2e-13). It has no H2+ + He -> HeH+ + H channel,
      !     which Garcia Munoz (2025) Table A.6 (Black 1978) gives as
      !     1.0e-11 at 2000 K rising to 1.5e-10 at 10^4 K -- a route that
      !     needs no He+ and is therefore the one that grows with the He
      !     fraction.  Nor does it carry the ~10% associative branch of
      !     He(2^3S) + H / + H2 (Garcia Munoz 2025, Appendix), which also
      !     ends in HeH+.
      !   * HeH+ destruction and He+ + H2. Against the Garcia Munoz (2025)
      !     Tables A.6/A.7 values, R16 is 3.4-8.6x smaller, R18 1.2x larger,
      !     R19 1.4-2.6x smaller, R20 14x larger and R17 larger by up to
      !     1.9e4 at 10^4 K.  The last two follow from the same choice:
      !     Garcia Munoz assigns the Schauer et al. (1989) total, 3e-14 and
      !     T-independent, entirely to the dissociative channel, where
      !     Table 1 gives that channel the Moses & Bass (2000) Jovian-
      !     ionosphere fit and the HeH+ channel 4.2e-13 on top of it.  R17
      !     is the dominant H2 sink and H+ source once He+ is present, so
      !     the choice matters most exactly in the He-dominated limit.
      !
      ! NOTE (Koskinen 2022 baseline): neutral H2 photodissociation through
      ! the Lyman-Werner bands is NOT part of their Table 1; their
      ! sensitivity test with a Backx et al. (1976) cross section and
      ! dissociation probability 0.125 changed Mdot by <= 1.4x.  EXHALE adds
      ! it separately and opt-in, in
      ! src/modules/lower_atmosphere/lyman_werner.f90 (Draine & Bertoldi
      ! 1996; key "Stellar LW flux"), so it is deliberately absent here.

      ! Physical constants and the hydrogen-atom mass, from the one place
      ! that owns them.  Nothing else of global_parameters is in scope here,
      ! so the local names of this module are checked against these five and
      ! against the module's own declarations only.
      use global_parameters, only: kb_erg, hp_erg, c_light, mu, pi

      implicit none
      private
      public :: rk_R1_Hp_rec,   rk_R2_Hep_rec,  rk_R3_H_cion,            &
                rk_R4_He_cion,  rk_R5_H2p_dr,   rk_R6_H3p_dr_H2,         &
                rk_R7_H3p_dr_3H, rk_R8_H2p_H2,  rk_R9_H2p_H,             &
                rk_R10_Hp_H2v4, rk_R11_H3p_H,   rk_R12_H2_thdis,         &
                rk_R13_Hp_H2_M, rk_R14_H2_edis, rk_R15_3body_H2,         &
                rk_R16_HeHp_dr, rk_R17_Hep_H2_diss, rk_R18_HeHp_H2,      &
                rk_R19_HeHp_H,  rk_R20_Hep_H2_HeHp, rk_R21_H_Hep_cx,     &
                rk_R22_Hp_He_cx, rk_R23_H2_Hep_cx
      ! H2 thermochemistry: the three-body recombination coefficient that R15
      ! and R12 share, the H + H <-> H2 equilibrium constant that turns one
      ! into the other, its rovibrational partition function, and the table
      ! builder.  See the H2 THERMOCHEMISTRY block below.
      public :: k3b_H_H_to_H2, keq_H_H_to_H2, q_rovib_H2,                 &
                h2_thermochemistry_init

      ! ------------------------------------------------------------------ !
      ! H2 THERMOCHEMISTRY
      !
      ! Spectroscopic constants of H2 X^1 Sigma_g^+ (Huber & Herzberg 1979,
      ! "Constants of Diatomic Molecules"; the standard tabulated set) and
      ! the dissociation energy from v = 0, J = 0.  All in cm^-1.
      !
      ! D0 IS CHECKED AGAINST THE THERMOCHEMISTRY OF THE SAME REACTION.
      ! 36118.11 cm^-1 = 432.05 kJ/mol; adding the 298.15 K enthalpy
      ! functions of two H atoms (2 x 5/2 RT) and removing that of H2
      ! (7/2 RT) gives dH(298.15 K) = 435.8 kJ/mol, against the 436 kJ/mol
      ! tabulated on the H2 (+M) <-> H + H (+M) data sheet of Baulch et al.
      ! (1992) J. Phys. Chem. Ref. Data 21, 411, p. 550.
      real*8, parameter :: D0_H2_cm    = 36118.11d0   ! D0(v=0,J=0)
      real*8, parameter :: we_H2       = 4401.213d0   ! omega_e
      real*8, parameter :: wexe_H2     = 121.336d0    ! omega_e x_e
      real*8, parameter :: Be_H2       = 60.853d0     ! B_e
      real*8, parameter :: alpha_e_H2  = 3.062d0      ! vib-rot coupling

      ! Tabulation of ln K_eq on a log-spaced temperature grid.  The
      ! rovibrational sum of q_rovib_H2 costs some 800 exponentials, far too
      ! many to evaluate at every cell inside the Newton solve of the molecular
      ! system (~1e9 over a run), so it is evaluated once here and
      ! interpolated linearly in (ln T, ln K_eq).
      !
      ! GRID.  100-20000 K, 2001 points, i.e. 2.65e-3 in ln T.  ln K_eq is
      ! dominated by D0 hc / kT, whose curvature in ln T is D0 hc / kT and
      ! largest at the cold end; the linear-interpolation error is therefore
      ! bounded by (1/8) (D0 hc/kT) (dlnT)^2 = 4.6e-4 in ln K_eq at 100 K,
      ! i.e. 0.046 % in K_eq, and falls as 1/T from there.  Measured against
      ! the direct sum at the midpoint of every interval, the largest error
      ! is 4.69e-4 (0.047 %) in the first interval and 1.54e-6 (1.5e-4 %) in
      ! the last.
      !
      ! ENDS.  T is clamped into [100, 20000] K before the lookup.  Below
      ! 100 K the clamp overestimates K_eq's fall and so overestimates the
      ! dissociation rate, but K_eq(100 K) = 1.2e202 cm^3 already makes that
      ! rate 1e-234 cm^3 s^-1 -- and 100 K is also as cold as K_eq may be
      ! evaluated at all before exp(D0 hc/kT) overflows double precision.
      ! Above 20000 K the clamp underestimates the dissociation rate, at a
      ! temperature where it is already ~1e-10 cm^3 s^-1 times the total
      ! density and H2 cannot survive a single time step either way.
      !
      ! The table is written once, by h2_thermochemistry_init, and only read
      ! afterwards, so it is NOT threadprivate: every thread of the parallel
      ! cell sweep shares one copy.  EXHALE_main calls the initializer before
      ! any parallel region opens, exactly as it does for the infrared band
      ! tables; the guard inside keq_H_H_to_H2 is a safety net for any other
      ! entry point and never fires in the cell sweep.
      integer, parameter :: n_keq    = 2001
      real*8,  parameter :: T_keq_lo = 1.0d2
      real*8,  parameter :: T_keq_hi = 2.0d4
      real*8,  save      :: keq_lnT(n_keq), keq_lnK(n_keq)
      logical, save      :: keq_table_ready = .false.
      ! ------------------------------------------------------------------ !

      contains

      ! R1: H+ + e -> H + hv                    (Storey & Hummer 1995)
      double precision function rk_R1_Hp_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.0d-12*(300.0d0/Te)**0.64d0
      end function

      ! R2: He+ + e -> He + hv                  (Storey & Hummer 1995)
      double precision function rk_R2_Hep_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.6d-12*(300.0d0/Te)**0.64d0
      end function

      ! R3: H + e -> H+ + 2e                    (Voronov 1997)
      ! U = 13.6 eV / E_e;  E_e = kB Te in eV.
      double precision function rk_R3_H_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 13.6d0/(8.6173d-5*Te)
      k = 2.91d-8*U**0.39d0*exp(-U)/(0.232d0 + U)
      end function

      ! R4: He + e -> He+ + 2e                  (Voronov 1997)
      double precision function rk_R4_He_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 24.6d0/(8.6173d-5*Te)
      k = 1.75d-8*U**0.35d0*exp(-U)/(0.180d0 + U)
      end function

      ! R5: H2+ + e -> H + H                    (Auerbach et al. 1977)
      double precision function rk_R5_H2p_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.3d-8*(300.0d0/Te)**0.4d0
      end function

      ! R6: H3+ + e -> H2 + H                   (Larsson et al. 2008)
      double precision function rk_R6_H3p_dr_H2(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.16d-8*(300.0d0/Te)**0.65d0
      end function

      ! R7: H3+ + e -> H + H + H                (Larsson et al. 2008)
      double precision function rk_R7_H3p_dr_3H(Te) result(k)
      real*8, intent(in) :: Te
      k = 5.04d-8*(300.0d0/Te)**0.65d0
      end function

      ! R8: H2+ + H2 -> H3+ + H                 (Theard & Huntress 1974)
      double precision function rk_R8_H2p_H2() result(k)
      k = 2.0d-9
      end function

      ! R9: H2+ + H -> H+ + H2                  (Karpas et al. 1979)
      double precision function rk_R9_H2p_H() result(k)
      k = 6.4d-10
      end function

      ! R10: H+ + H2(v>=4) -> H2+ + H           (Yelle 2004)
      double precision function rk_R10_Hp_H2v4(T) result(k)
      real*8, intent(in) :: T
      k = 1.0d-9*exp(-21900.0d0/T)
      end function

      ! R11: H3+ + H -> H2+ + H2                (Harada et al. 2010)
      double precision function rk_R11_H3p_H(T) result(k)
      real*8, intent(in) :: T
      k = 2.1d-9*exp(-20000.0d0/T)
      end function

      ! R12: H2 + M -> H + H + M   (detailed balance of R15; see below)
      ! Two-body-equivalent: multiply by the third-body density n_M (see the
      ! THIRD BODY M note in the module header).
      !
      ! WHAT THIS RATE IS.  Not an independent fit.  R12 and R15 are the same
      ! channel run in opposite directions, so detailed balance genuinely
      ! applies to the pair:
      !     k_diss(T) = k_rec(T) / K_eq(T),   K_eq(T) = n_H2 / n_H^2  [cm^3]
      ! with k_rec = k3b_H_H_to_H2 (cm^6 s^-1, the coefficient of R15 WITHOUT
      ! its third-body density factor) and K_eq = keq_H_H_to_H2.  cm^6 s^-1
      ! divided by cm^3 leaves the cm^3 s^-1 this function returns.
      !
      ! (2) WHY DETAILED BALANCE RATHER THAN THE PUBLISHED DISSOCIATION FIT.
      ! The primary measurement behind the Table-1 dissociation fit was
      ! itself interpreted through detailed balance.  Breshears & Bird (1973)
      ! analyze their shock-tube data with the rate law
      !     -d(A2)/dt = k_d (M) [ (A2) - (A)^2 / K_eq ]
      ! (their Eq. 2) under the relation k_d / k_r = K_eq (their Eq. 3),
      ! introduced with the words "is assumed".  Deriving one direction of
      ! this channel from the other through K_eq is therefore not an
      ! assumption added here; it is the relation the measured coefficient
      ! already carries.
      !
      ! What that leaves is which direction to anchor on, and the Koskinen
      ! Table-1 pair does not overlap in temperature:
      !   * recombination -- Ham, Trainor & Kaufman (1970) J. Chem. Phys. 53,
      !     4395: (8.3 +/- 0.4)e-33 cm^6 s^-1 at 298 K following T^-0.6,
      !     measured over 77-300 K;
      !   * dissociation -- Baulch et al. (1992) J. Phys. Chem. Ref. Data 21,
      !     411, Table 2 and the data sheet on p. 550: 1.5e-9 exp(-48350/T)
      !     cm^3 s^-1 recommended over 2500-8000 K.
      ! The molecular layer modelled here sits at ~1000-2000 K, BETWEEN the
      ! two ranges, and nothing makes the pair thermodynamically consistent
      ! where neither was measured: taken at face value they imply K_eq
      ! values differing from the thermochemical one by up to a factor 10 at
      ! 1000 K.
      !
      ! The recombination is the direction anchored on, because a
      ! recombination coefficient exists whose published range covers the
      ! layer while no dissociation fit does.  R15 therefore carries Cohen &
      ! Westberg's (1983) recommended k1(H2) = 2.8e-31 T^-0.6 over
      ! 50-5000 K rather than Ham's 77-300 K measurement, which it
      ! reproduces in shape and exceeds by 14 % in amplitude (see R15 for the
      ! full argument), and R12 is its thermodynamic reverse.  Reversing the
      ! recombination is also what both reference networks in the workspace
      ! do with the near-identical Baulch p. 495 coefficient (VULCAN reaction
      ! 531, Photochem zahnle_earth "H + H (+ M) <=> H2 (+ M)", both
      ! 2.70e-31 T^-0.6).
      !
      ! WHICH VIBRATIONAL STATE THE COEFFICIENTS REFER TO is answered by the
      ! source rather than inferred.  Breshears & Bird read their initial
      ! rate off the minimum of the postshock density-gradient "dip", i.e.
      ! "the rate corresponding to conditions of vibrational equilibrium in
      ! the thermal sense but with negligible extent of dissociation".  So
      ! the dissociation coefficient is a vibrationally equilibrated rate,
      ! which is the state K_eq assumes on both sides; and at Ham's 77-300 K
      ! H2 is in v = 0 in any case.  The same paper bounds where its rate law
      ! applies: the "observation of induction periods for diatomic
      ! dissociation clearly indicates that the assumptions implicit in Eqs.
      ! (2)-(4) are not valid during the very early course of the reaction;
      ! nevertheless, they may be applicable during the 'steady state' which
      ! follows the induction period."  A hydrodynamic model that resolves
      ! nothing faster than a vibrational relaxation time is in that steady
      ! state throughout.
      !
      ! (1) THE THIRD BODY IS H2.  The anchored coefficient is Cohen &
      ! Westberg's k1(H2), labelled M = H2 on their data sheet and separate
      ! from their k1(H); the studies behind it are H2-as-third-body flow
      ! tubes, Ham et al.'s among them.  Every fit compared against here is
      ! the M = H2 one too: Baulch's "H2 + H2 -> 2H + H2" row, not the
      ! "H2 + Ar" row printed next to it, and Breshears & Bird's k_d^H2.
      ! (Breshears & Bird did not shock pure H2: they measured 5, 10 and
      ! 20 % H2 in Ar and in Xe and extracted k(M = H2) through the linear
      ! composition rule k_d = sum_M X_M k_d^M of their Eq. 4, whose validity
      ! is one of the things their precision was meant to test.)  What is
      ! derived here is therefore k(M = H2).  EXHALE passes the TOTAL
      ! heavy-particle density as M (see the THIRD BODY M note in the
      ! module header), i.e. it
      ! assumes every collider is as efficient as H2.  That assumption is
      ! weakest ABOVE the H2 front, where atomic H is the dominant collider
      ! -- which is also where H2 is nearly gone, so the error it makes is
      ! carried by a vanishing H2 density.
      !
      ! (3) VERIFICATION.  Measured with this implementation:
      !   k_diss(detailed balance) / k_diss(Baulch 1992 fit) =
      !     1.05 at 8000 K, 1.18 at 5000 K, 0.97 at 3000 K, 0.82 at 2500 K,
      !     0.61 at 2000 K, 0.36 at 1500 K, 0.11 at 1000 K.
      !   The agreement INSIDE Baulch's stated 2500-8000 K range, and the
      !   departure only below it, is the evidence that the construction is
      !   right: where the fit was measured the two agree, where it is an
      !   extrapolation they part company.  Against the PRIMARY measurement
      !   rather than the evaluation -- Breshears & Bird's own M = H2
      !   coefficient 5.48e-9 exp(-52989/T), their abstract's 3.30e15
      !   exp(-105300/RT) cc/mole/s converted -- the same ratio is 1.11 at
      !   3500 K, the bottom of their measured 3500-8000 K range, 0.99 at
      !   4000 K, 0.82 at 5000 K and 0.51 at 8000 K.  Evaluation and
      !   measurement are not the same curve: BB73 divided by Baulch is 0.035
      !   at 1000 K, 0.36 at 2000 K, 0.78 at 3000 K, 1.15-1.45 over
      !   4000-5000 K and 2.05 at 8000 K, matching in the middle and parting
      !   at both ends as a synthesis of several data sets does.
      ! K_eq itself has three independent confirmations:
      !   * the chemical-equilibrium H2 fit this code already carries,
      !     q_h2_equilibrium (lower_column.f90, Koskinen et al. 2022 Eq. 11
      !     quoting Visscher et al. 2006), evaluated in its DILUTE limit
      !     where n_H2 << n_H and q -> 0.831856 / 10^u: that K_eq divided by
      !     this one is 0.993 at 2000 K, 1.029 at 2500 K, 1.082 at 3000 K.
      !   * the NIST-JANAF Shomate enthalpies and entropies of H and H2 that
      !     oxygen_rates.f90 already carries: K_eq from dG(T) divided by this
      !     one is 1.010 at 600 K, 1.018 at 1000 K, 1.037 at 2000 K, 1.056 at
      !     3000 K, 1.101 at 5000 K.  The slow rise is the truncation of the
      !     Morse level sum at D0, which drops the quasi-bound levels a real
      !     H2 molecule still has.
      !   * Cohen & Westberg (1983) recommend BOTH directions on the same
      !     data sheet, so their ratio k1(H2)/k-1(H2) is an equilibrium
      !     constant from an evaluation that used none of the above:
      !     1.958e-25 T^0.1 exp(52530/T) cm^3 over 600-5000 K, which divided
      !     by this K_eq is 1.019 at 600 K, 0.966 at 1000 K, 1.063 at
      !     2000 K, 1.118 at 3000 K and 1.108 at 5000 K.  (It is also the
      !     reason the k_diss above is not simply k-1(H2): the two agree to
      !     those same 12 %, and taking K_eq from the physics instead of from
      !     a two-parameter fit is what keeps R12 and R15 exactly reversible
      !     inside the solver.)
      !
      ! (4) A COMPARISON AGAINST ANOTHER NETWORK IS MEANINGLESS UNLESS THE
      ! THIRD BODY MATCHES.  Breshears & Bird measured all three colliders in
      ! one experiment with one analysis, so their DISSOCIATION coefficients
      ! are the clean statement of how much the third body matters:
      !     M = Ar  1.55e-10 exp(-44736/T)
      !     M = H2  5.48e-9  exp(-52989/T)
      !     M = H   3.52e-9  exp(-43881/T)   [cm^3 s^-1]
      ! Inside their measured 3500-8000 K range k(M = H)/k(M = H2) is 8.7 at
      ! 3500 K and 4.0 at 5000 K, and Cohen & Westberg's independent
      ! RECOMBINATION pair -- k1(H) = 8.8e-33, temperature-independent over
      ! 50-5000 K, against k1(H2) = 2.8e-31 T^-0.6 -- gives 4.2 and 5.2 at
      ! the same two temperatures.  Detailed balance forces the two ratios to
      ! be equal, and where both sources are inside their ranges they are.
      ! They are NOT equal where the sources are extrapolated: at 1000 K the
      ! Breshears ratio reaches 5.8e3 (its two colliders have activation
      ! temperatures 52989 and 43881 K, so the ratio runs away below the
      ! measured range) while Cohen & Westberg give 2.0.  Take 2 to 5 as the
      ! collider spread where H2 exists, and read the 5.8e3 as what
      ! extrapolating an Arrhenius fit costs.  Note also that C&W's M = H
      ! recommendation is itself uncertain by +/-0.5 in log throughout.
      ! A second axis has to match as well: rates quoted as the n -> 0
      ! collisional-dissociation limit, with H2 held in v = 0 (Glover &
      ! Jappsen 2007, ApJ 666, 1, R9 and R10), are ~4e9 below this rate at
      ! 1000 K for the SAME third body, because this one is the vibrationally
      ! equilibrated regime of item (2).
      !
      ! (5) SOURCES, all read from the published papers in references/:
      !   Ham, D. O., Trainor, D. W. & Kaufman, F. 1970, J. Chem. Phys. 53,
      !     4395 (Ham_1970JCP.pdf)
      !   Baulch, D. L. et al. 1992, J. Phys. Chem. Ref. Data 21, 411,
      !     p. 550 and p. 495 (Baulch_1992JPCRD.pdf)
      !   Breshears, W. D. & Bird, P. F. 1973, 14th Symposium (International)
      !     on Combustion, 211 (Breshears_1973.pdf) -- the shock-tube
      !     measurement behind the Baulch dissociation recommendation
      !   Cohen, N. & Westberg, K. R. 1983, J. Phys. Chem. Ref. Data 12, 531,
      !     p. 559 (Cohen_1983JPCRD_12_531.pdf) -- the recombination
      !     coefficient R15 carries, and both directions of the pair
      double precision function rk_R12_H2_thdis(T) result(k)
      real*8, intent(in) :: T
      k = k3b_H_H_to_H2(T)/keq_H_H_to_H2(T)
      end function

      ! R13: H+ + H2 + M -> H3+ + M             (Miller et al. 1968)
      ! Returns the 2-body-equivalent rate: 3.2e-29 * n  [cm^3 s^-1], with n
      ! the third-body density (see the THIRD BODY M note in the header).
      double precision function rk_R13_Hp_H2_M(n) result(k)
      real*8, intent(in) :: n
      k = 3.2d-29*n
      end function

      ! R14: H2 + e -> H + H + e                (Stibbe & Tennyson 1999)
      double precision function rk_R14_H2_edis(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.33d-6*(300.0d0/Te)**0.91d0*exp(-55800.0d0/Te)
      end function

      ! Three-body recombination coefficient of H + H + H2 -> H2 + H2
      ! [cm^6 s^-1], Cohen & Westberg (1983) J. Phys. Chem. Ref. Data 12,
      ! 531, p. 559, recommended k1(H2) over 50-5000 K.  Defined here once
      ! because BOTH directions of this channel read it: R15 multiplies it by
      ! the third-body density, R12 divides it by the equilibrium constant.
      ! Provenance and validity are on R15 below; the reason R12 is built
      ! from it is on R12 above.
      double precision function k3b_H_H_to_H2(T) result(k)
      real*8, intent(in) :: T
      k = 2.8d-31*T**(-0.6d0)
      end function

      ! R15: H + H + M -> H2 + M           (Cohen & Westberg 1983)
      ! Returns the 2-body-equivalent rate: 2.8e-31 T^-0.6 * n, with n the
      ! third-body density (see the THIRD BODY M note in the header).
      ! Table 1 prints the temperature of this fit as Te; the reaction is a
      ! neutral three-body recombination with no electron in it, so the
      ! heavy-particle T is used.  EXHALE carries a single T, so the two
      ! readings coincide numerically.
      !
      ! WHICH COEFFICIENT, AND WHY NOT KOSKINEN'S.  Koskinen et al. (2022)
      ! Table 1 gives this rate as 8e-33 (300/T)^0.6, the room-temperature
      ! measurement of Ham, Trainor & Kaufman (1970) J. Chem. Phys. 53, 4395:
      ! k_R = (8.3 +/- 0.4)e-33 cm^6 s^-1 at 298 K, following T^-0.6 "fairly
      ! closely", over 77-300 K.  Written as A T^-0.6 that is 2.45e-31, and
      ! the molecular layer here (~1000-2000 K) sits three to seven times
      ! ABOVE the measured range.  What is used instead is the recommended
      ! k1(H2) of Cohen & Westberg (1983) J. Phys. Chem. Ref. Data 12, 531,
      ! p. 559: 2.8e-31 T^-0.6 cm^6 s^-1 over 50-5000 K, k(298) = 9.0e-33.
      ! Its range covers the layer outright, and its exponent is Ham's own
      ! regime carried across that range -- "the low temperature (77-300 K)
      ! data for M = H2 are in reasonable agreement and suggest a temperature
      ! dependence of T^-0.6 for k1(H2)" -- with Ham et al. one of the four
      ! studies the recommendation rests on.  The change of source is 14 % in
      ! amplitude and nothing in shape.  A third, independent evaluation
      ! agrees: Baulch et al. (1992) J. Phys. Chem. Ref. Data 21, 411,
      ! Table 3 and p. 495 recommend ko = 2.7e-31 T^-0.6 over 100-5000 K,
      ! 3.6 % below this one, and it is the coefficient both reference
      ! networks in the workspace carry (VULCAN reaction 531, Photochem
      ! zahnle_earth, both 2.70e-31 T^-0.6).
      !
      ! Stated uncertainty on the value used: +/-0.2 in log k1(H2) at 300 K,
      ! rising to +/-0.4 at 5000 K.
      !
      ! WHY NO HIGH-PRESSURE LIMIT IS CARRIED.  A falloff term would be
      ! needed where k_0 n approaches the high-pressure limit k_inf. The two
      ! reference networks in the workspace disagree about k_inf by three and
      ! a half orders (VULCAN 3.31e-6 T^-1, Photochem 1.0e-12 cm^3 s^-1), so
      ! the onset density n_crit = k_inf/k_0 was evaluated with BOTH.
      ! Measured on the five molecular regression cases and both diagnostic
      ! arms, the largest molecular-layer particle density is 5.4e13 cm^-3
      ! (matrix) and 7.2e14 cm^-3 (He/H = 10 arm), against n_crit >=
      ! 4.4e19 cm^-3 from the STRICTER of the two at this k_0. Every
      ! configuration is therefore at least four and a half orders below the
      ! onset, and the two references' disagreement never enters. The
      ! low-pressure form is exact here; no k_inf is carried.
      double precision function rk_R15_3body_H2(T, n) result(k)
      real*8, intent(in) :: T, n
      k = k3b_H_H_to_H2(T)*n
      end function

      ! ------------------------------------------------------------------ !
      ! H + H <-> H2 equilibrium constant and its rovibrational partition
      ! function. See the H2 THERMOCHEMISTRY block above the contains
      ! statement for the constants, the table and its error, and R12 for
      ! what the equilibrium constant is used for.
      ! ------------------------------------------------------------------ !

      ! Internal (rovibrational) partition function of H2, measured from
      ! v = 0, J = 0, with the ortho:para nuclear-spin weights 3:1 included:
      !     G(v)   = we (v+1/2) - wexe (v+1/2)^2
      !     B_v    = Be - alpha_e (v+1/2)
      !     E(v,J) = [G(v) - G(0)] + B_v J(J+1)
      !     q      = sum_v sum_J g_ns(J) (2J+1) exp(-E/kT), g_ns = 3 odd J,
      !                                                            1 even J
      ! The sums stop where the level is no longer bound: the J loop at
      ! E > D0, the v loop at B_v <= 0 or G(v) - G(0) > D0.
      !
      ! NO CENTRIFUGAL DISTORTION TERM. Subtracting De [J(J+1)]^2 with
      ! De = 4.71e-2 cm^-1 drives E NEGATIVE at moderate J, so the J loop
      ! never reaches its D0 stop and the sum runs away: q(2000 K) comes out
      ! 798, a factor 16 too large. A distortion term would need its own
      ! bound on J (the top of the rotational barrier), which is not worth
      ! carrying for a correction below 1 % where the levels are actually
      ! populated. As truncated here q(2000 K) = 50.22, which the classical
      ! estimate 4 kT / (2 Be) x q_vib = 48.1 confirms; the agreement of
      ! K_eq with the three independent routes quoted on R12 is the
      ! quantitative test.
      double precision function q_rovib_H2(T) result(q)
      real*8, intent(in) :: T
      real*8  :: hc_k, g_v0, g_v, b_v, v_half, e_vj, w_ns
      integer :: iv, jr
      hc_k = hp_erg*c_light/kb_erg          ! K per cm^-1
      g_v0 = 0.5d0*we_H2 - 0.25d0*wexe_H2   ! G(v=0)
      q    = 0.0d0
      iv   = 0
      do
         v_half = dble(iv) + 0.5d0
         b_v    = Be_H2 - alpha_e_H2*v_half
         g_v    = we_H2*v_half - wexe_H2*v_half*v_half - g_v0
         if (b_v .le. 0.0d0 .or. g_v .gt. D0_H2_cm) exit
         jr = 0
         do
            e_vj = g_v + b_v*dble(jr)*dble(jr+1)
            if (e_vj .gt. D0_H2_cm) exit
            if (mod(jr,2) .eq. 1) then
               w_ns = 3.0d0                 ! ortho (odd J)
            else
               w_ns = 1.0d0                 ! para  (even J)
            endif
            q  = q + w_ns*dble(2*jr+1)*exp(-e_vj*hc_k/T)
            jr = jr + 1
         enddo
         iv = iv + 1
      enddo
      end function q_rovib_H2

      ! Build the ln K_eq table. Called from EXHALE_main before any parallel
      ! region opens; idempotent, so a second call costs nothing.
      !
      !     K_eq(T) = Lam(2 m_H) q_int(T) exp(D0/kT) / [ Lam(m_H) x 4 ]^2
      !     Lam(m)  = (2 pi m kB T / h^2)^(3/2)                    [cm^-3]
      !
      ! The 4 on the H atom is 2 (electronic doublet) x 2 (nuclear spin
      ! I = 1/2); the 3:1 weights inside q_rovib_H2 are the matching
      ! nuclear-spin counting on the molecule, so both sides of the
      ! equilibrium are on the same convention and the spin factors cancel.
      subroutine h2_thermochemistry_init
      integer :: i
      real*8  :: ln_lo, d_ln, t_k, lam_h, lam_h2, hc_k
      if (keq_table_ready) return
      hc_k  = hp_erg*c_light/kb_erg
      ln_lo = log(T_keq_lo)
      d_ln  = (log(T_keq_hi) - ln_lo)/dble(n_keq - 1)
      do i = 1, n_keq
         keq_lnT(i) = ln_lo + d_ln*dble(i - 1)
         t_k        = exp(keq_lnT(i))
         lam_h      = (2.0d0*pi*mu*kb_erg*t_k/hp_erg**2)**1.5d0
         lam_h2     = (2.0d0*pi*2.0d0*mu*kb_erg*t_k/hp_erg**2)**1.5d0
         keq_lnK(i) = log(lam_h2) + log(q_rovib_H2(t_k))                  &
                    + D0_H2_cm*hc_k/t_k - 2.0d0*log(4.0d0*lam_h)
      enddo
      keq_table_ready = .true.
      end subroutine h2_thermochemistry_init

      ! K_eq(T) = n_H2 / n_H^2 [cm^3] at chemical equilibrium, interpolated
      ! from the table. T is clamped into [100, 20000] K (see the ENDS note
      ! above the contains statement).
      double precision function keq_H_H_to_H2(T) result(keq)
      real*8, intent(in) :: T
      real*8  :: ln_t, frac
      integer :: i
      if (.not. keq_table_ready) then
!$omp critical (h2_thermochemistry_table)
         if (.not. keq_table_ready) call h2_thermochemistry_init
!$omp end critical (h2_thermochemistry_table)
      endif
      ln_t = log(min(max(T, T_keq_lo), T_keq_hi))
      i    = int((ln_t - keq_lnT(1))/(keq_lnT(2) - keq_lnT(1))) + 1
      i    = min(max(i, 1), n_keq - 1)
      frac = (ln_t - keq_lnT(i))/(keq_lnT(i+1) - keq_lnT(i))
      keq  = exp(keq_lnK(i) + frac*(keq_lnK(i+1) - keq_lnK(i)))
      end function keq_H_H_to_H2

      ! R16: HeH+ + e -> He + H                 (Yousif & Mitchell 1989)
      ! 3.4x below the Florescu-Mitchell & Mitchell (2006) values tabulated
      ! by Garcia Munoz (2025) Table A.6 at 500 K and 8.6x below at 10^4 K:
      ! both fall with T, this one faster.  See the HE-DOMINATED LIMIT note
      ! in the header.
      double precision function rk_R16_HeHp_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.0d-8*(300.0d0/Te)**0.6d0
      end function

      ! R17: He+ + H2 -> H+ + H + He            (Moses & Bass 2000)
      ! The largest cross-source discrepancy in this network: Garcia Munoz
      ! (2025) Table A.7 gives 3e-14, T-independent (the Schauer et al. 1989
      ! total assigned to this channel), against 5.8e-11 at 2000 K and
      ! 5.7e-10 at 10^4 K here.  See the HE-DOMINATED LIMIT note above.
      double precision function rk_R17_Hep_H2_diss(T) result(k)
      real*8, intent(in) :: T
      k = 1.0d-9*exp(-5700.0d0/T)
      end function

      ! R18: HeH+ + H2 -> H3+ + He              (Bohme et al. 1980)
      double precision function rk_R18_HeHp_H2() result(k)
      k = 1.5d-9
      end function

      ! R19: HeH+ + H -> H2+ + He               (Karpas et al. 1979)
      double precision function rk_R19_HeHp_H() result(k)
      k = 9.1d-10
      end function

      ! R20: He+ + H2 -> HeH+ + H               (Schauer et al. 1989)
      double precision function rk_R20_Hep_H2_HeHp() result(k)
      k = 4.2d-13
      end function

      ! R21 and R22 are transcribed for completeness but are NOT the rates
      ! the code uses: the H <-> He charge-exchange pair is applied in every
      ! system that carries He, from Huang et al. (2023) Table 4 rows B1/B2
      ! (charge_exchange::he_h_cx_rates, which is where the live definition
      ! and its provenance are), whose values agree with these to 4% (R21)
      ! and to the rounding of the 128,000 K barrier (R22).
      !
      ! These two are NOT each other's reverse, and it is a mistake to test
      ! them against detailed balance -- as this comment did before
      ! 2026-08-31. R21 is RADIATIVE charge transfer, He+ + H -> He + H+ +
      ! photon (Stancil, Lepp & Dalgarno 1998, Table 1 row 19, from
      ! Zygelman et al. 1989); R22 is the NON-RADIATIVE collisional channel
      ! (Kimura et al. 1993). A photon-emitting process has no collisional
      ! reverse, so k(H+ + He)/k(He+ + H) is under no obligation to equal
      ! 4 exp(-127,500 K / T), and the ~10^2-10^3 factor by which it differs
      ! measures the separation of two channels rather than an error in
      ! either fit. See docs/molecular_chemistry_audit_he_rich.md section
      ! 4.3 and the Group B block of charge_exchange.f90.
      !
      ! R21: H + He+ -> H+ + He                 (Stancil et al. 1998)
      double precision function rk_R21_H_Hep_cx(T) result(k)
      real*8, intent(in) :: T
      k = 1.2d-15*(300.0d0/T)**(-0.25d0)
      end function

      ! R22: H+ + He -> H + He+                 (Glover & Jappsen 2007)
      double precision function rk_R22_Hp_He_cx(T) result(k)
      real*8, intent(in) :: T
      k = 1.75d-11*(300.0d0/T)**0.75d0*exp(-128000.0d0/T)
      end function

      ! R23: H2 + He+ -> H2+ + He               (Barlow 1984)
      double precision function rk_R23_H2_Hep_cx() result(k)
      k = 7.2d-15
      end function

      ! End of module
      end module mol_rates

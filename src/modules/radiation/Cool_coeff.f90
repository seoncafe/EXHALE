   module Cooling_Coefficients
   ! Module containing rate coefficients and cooling rates 
   ! Included processes: bremsstrahlung, collisional ionization,
   !     collisional excitation, recombination
   ! 
   ! --------------------
   ! 
   ! - AC, 14/04/22: Added He metastable collisional strength 
   !	coefficients
   
   use global_parameters
   ! Koskinen et al. (2022, ApJ 929, 52) Table 1 R1 and R2, transcribed once
   ! in mol_rates and selected here by "Atomic rate set: Koskinen2022". R3
   ! and R4 of the same table are the Voronov (1997) fits this module owns
   ! (voronov_ci, ci_HI_new, ci_HeI_new).
   use mol_rates,     only: rk_R1_Hp_rec, rk_R2_Hep_rec
   ! The ground-state photoionization cross sections of the radiative
   ! transfer, for the Milne relation of the ground-state captures
   ! (ground_capture_milne). Private here, so that a module using this one
   ! does not see them a second time.
   use Cross_sections, only: sigma, sigma_HeI
   use species_table, only: n_mion, melem_A_u,                          &
                            mion_iscool, mion_stage, mion_z2, mion_ethr,   &
                            iel_C, iel_N, iel_O,                           &
        im_CI, im_CII, im_CIII, im_OI, im_OII, im_OIII,                &
        im_NI, im_NII, im_NIII, im_MgI, im_MgII, im_MgIII,             &
        im_SiI, im_SiII, im_SiIII, im_CaI, im_CaII, im_CaIII,          &
        im_NaI, im_NaII, im_KI, im_KII, im_S_I, im_S_II,               &
        im_FeI, im_FeII, im_FeIII

   implicit none

   private :: sigma, sigma_HeI

   !--- Badnell total recombination fits (RR + DR) ------------------!
   ! Modern replacement for the Aldrovandi & Pequignot (1973) power-law
   ! recombination used previously (RR only, no dielectronic term).
   ! Ported from ATES_extended (metals_solve.f90).
   !
   ! RR: Badnell 2006 (ApJS 167, 334), total ground-level fit
   !   alpha_RR = A / [ sqrt(T/T0) (1+sqrt(T/T0))^(1-b)
   !                                (1+sqrt(T/T1))^(1+b) ],
   !   b = B + C*exp(-T2/T)
   ! DR: Badnell adf09 total ground-level fit
   !   alpha_DR = T^(-3/2) * sum_i dr_c(i) exp(-dr_e(i)/T)
   ! Source tables (clist_K): https://amdpp.phys.strath.ac.uk/tamoc/{RR,DR}
   !
   ! Keyed by the recombined (daughter) ion: bad_CI is the rate for
   ! CII + e -> CI, bad_NI for NII + e -> NI, etc.
   type :: rec_fit
      real*8  :: rr_A, rr_B, rr_T0, rr_T1, rr_C, rr_T2
      integer :: dr_n
      real*8, dimension(9) :: dr_c, dr_e
   end type rec_fit

   type(rec_fit), parameter :: bad_CI  = rec_fit(                         &
        2.995d-9, 0.7849d0, 6.670d-3, 1.943d6, 0.1597d0, 4.955d4, 5,      &
        [6.346d-9,9.793d-9,1.634d-6,8.369d-4,3.355d-4,0d0,0d0,0d0,0d0],   &
        [1.217d1,7.380d1,1.523d4,1.207d5,2.144d5,0d0,0d0,0d0,0d0])

   type(rec_fit), parameter :: bad_CII = rec_fit(                         &
        2.067d-9, 0.8012d0, 1.643d-1, 2.172d6, 0.0427d0, 6.341d4, 6,      &
        [3.489d-6,2.222d-7,1.954d-5,4.212d-3,2.037d-4,2.936d-4,0d0,0d0,0d0], &
        [2.660d3,3.756d3,2.566d4,1.400d5,1.801d6,4.307d6,0d0,0d0,0d0])

   type(rec_fit), parameter :: bad_NI  = rec_fit(                         &
        6.387d-10, 0.7308d0, 9.467d-2, 2.954d6, 0.2440d0, 6.739d4, 6,     &
        [1.658d-8,2.760d-8,2.391d-9,7.585d-7,3.012d-4,7.132d-4,0d0,0d0,0d0], &
        [1.265d1,8.425d1,2.964d2,5.923d3,1.278d5,2.184d5,0d0,0d0,0d0])

   type(rec_fit), parameter :: bad_NII = rec_fit(                         &
        2.410d-9, 0.7948d0, 1.231d-1, 3.016d6, 0.0774d0, 1.016d5, 6,      &
        [7.712d-8,4.839d-8,2.218d-6,1.536d-3,3.647d-3,4.234d-5,0d0,0d0,0d0], &
        [7.113d1,2.765d2,1.439d4,1.347d5,2.496d5,2.204d6,0d0,0d0,0d0])

   type(rec_fit), parameter :: bad_OI  = rec_fit(                         &
        6.622d-11, 0.6109d0, 4.136d0, 4.214d6, 0.4093d0, 8.770d4, 4,      &
        [5.629d-8,2.550d-7,6.173d-4,1.627d-4,0d0,0d0,0d0,0d0,0d0],        &
        [5.395d3,1.770d4,1.671d5,2.687d5,0d0,0d0,0d0,0d0,0d0])

   type(rec_fit), parameter :: bad_OII = rec_fit(                         &
        2.096d-9, 0.7668d0, 1.602d-1, 4.377d6, 0.1070d0, 1.392d5, 6,      &
        [1.627d-7,1.262d-7,6.663d-7,3.925d-6,2.406d-3,1.146d-3,0d0,0d0,0d0], &
        [4.535d1,2.847d2,4.166d3,2.877d4,1.953d5,3.646d5,0d0,0d0,0d0])

   ! Mg II + e -> Mg I (Mg isoelectronic sequence, Z=12, N_el=11)
   type(rec_fit), parameter :: bad_MgI = rec_fit(                         &
        5.452d-11, 0.6845d0, 5.637d0, 1.551d6, 0.3945d0, 8.360d5, 4,      &
        [3.871d-8,4.732d-7,1.599d-3,2.628d-5,0d0,0d0,0d0,0d0,0d0],        &
        [8.415d3,1.682d4,5.000d4,2.759d5,0d0,0d0,0d0,0d0,0d0])

   ! Mg III + e -> Mg II (Z=12, N_el=10)
   type(rec_fit), parameter :: bad_MgII = rec_fit(                        &
        1.345d-11, 0.1074d0, 7.877d2, 7.925d7, 0.4631d0, 5.027d5, 3,      &
        [6.269d-6,9.181d-4,3.082d-4,0d0,0d0,0d0,0d0,0d0,0d0],             &
        [4.104d5,5.766d5,7.310d5,0d0,0d0,0d0,0d0,0d0,0d0])

   ! --- New trace metals (Si, Ca, Na, K, S): Badnell TAMOC clist_K ---
   ! Same convention: bad_X is the rate for the recombining ion that
   ! produces daughter X (e.g. bad_SiI = Si II + e -> Si I). Na/K/S are
   ! two-stage so only the X+ -> X0 daughter (NaI/KI/SI) is tabulated.
   ! Ca I (recombining Ca II, K-like) is NOT in Badnell -> handled by the
   ! Shull & Van Steenberg 1982 power law in alpha_rr_metal (no entry here).
   ! Nor is iron: the recombining Fe II and Fe III are Mn-like and Cr-like,
   ! and the tabulation stops at the Ar-like sequence, so Fe I and Fe II take
   ! the Huang et al. 2023 fits near the end of this module.

   ! Si II + e -> Si I (Z=14, N_el=13)
   type(rec_fit), parameter :: bad_SiI = rec_fit(                         &
        3.262d-11, 0.6270d0, 1.590d1, 4.237d7, 0.2333d0, 5.828d4, 6,      &
        [3.408d-8,1.913d-7,1.679d-7,7.523d-7,8.386d-5,4.083d-3,0d0,0d0,0d0], &
        [2.431d1,1.293d2,4.272d2,3.729d3,5.514d4,1.295d5,0d0,0d0,0d0])

   ! Si III + e -> Si II (Z=14, N_el=12)
   type(rec_fit), parameter :: bad_SiII = rec_fit(                        &
        1.964d-10, 0.6287d0, 7.712d0, 2.951d7, 0.1523d0, 4.804d5, 5,      &
        [2.930d-6,2.803d-6,9.023d-5,6.909d-3,2.582d-5,0d0,0d0,0d0,0d0],   &
        [1.162d2,5.721d3,3.477d4,1.176d5,3.505d6,0d0,0d0,0d0,0d0])

   ! Ca III + e -> Ca II (Z=20, N_el=18)
   type(rec_fit), parameter :: bad_CaII = rec_fit(                        &
        2.248d-10, 0.6605d0, 6.175d0, 6.032d6, 0.3158d0, 2.100d5, 3,      &
        [3.843d-4,8.040d-3,8.670d-3,0d0,0d0,0d0,0d0,0d0,0d0],             &
        [2.282d5,3.682d5,4.479d5,0d0,0d0,0d0,0d0,0d0,0d0])

   ! Na II + e -> Na I (Z=11, N_el=10)
   type(rec_fit), parameter :: bad_NaI = rec_fit(                         &
        5.095d-12, 0.0000d0, 3.546d2, 2.310d6, 0.9395d0, 4.297d5, 3,      &
        [2.673d-6,1.918d-4,1.491d-5,0d0,0d0,0d0,0d0,0d0,0d0],             &
        [3.027d5,3.732d5,4.787d5,0d0,0d0,0d0,0d0,0d0,0d0])

   ! K II + e -> K I (Z=19, N_el=18)
   type(rec_fit), parameter :: bad_KI = rec_fit(                          &
        4.528d-11, 0.4234d0, 5.931d0, 2.897d9, 0.3049d0, 1.645d5, 4,      &
        [6.292d-4,2.350d-3,1.165d-2,2.470d-3,0d0,0d0,0d0,0d0,0d0],        &
        [2.451d5,3.504d5,4.094d5,4.766d5,0d0,0d0,0d0,0d0,0d0])

   ! S II + e -> S I (Z=16, N_el=15)
   type(rec_fit), parameter :: bad_SI = rec_fit(                          &
        1.384d-10, 0.6886d0, 1.074d0, 7.159d5, 0.1845d0, 1.858d4, 7,      &
        [7.300d-8,2.577d-7,4.961d-8,9.520d-7,9.586d-7,6.849d-4,6.539d-4,0d0,0d0], &
        [5.077d2,6.007d2,2.342d3,7.269d3,2.190d4,1.483d5,1.906d5,0d0,0d0])

   !--- CHIANTI collisional metal line cooling (Huang 2023) ---!
   ! Effective cooling coefficient Lambda(T) per (n_e * n_ion) [erg cm^3 s^-1],
   ! optically thin / coronal limit (only the ground level significantly
   ! populated; every collisional excitation is followed by a radiative decay).
   !
   ! Source: cooling_data/chianti_cooling.py + export_cooling_tables.py
   !   (CHIANTI v11.0.2 .scups; the Burgess & Tully (1992) scaled
   !   collision strengths are descaled with the natural cubic spline of
   !   CHIANTI's DESCALE_SCUPS.PRO through the stored points; auditable
   !   bridge).
   !   Mg I  2853 A  (3s2 1S0 -> 3s3p 1P1),  ground -> level 5: the coronal
   !         resonance-line skeleton. The cooling takes Mg I from the
   !         density-resolved table cool_logLrem_MgI, every channel.
   !   Mg II every excitation out of 3s 2S (h&k 2796/2804 A and the levels
   !         above 3p) with the OBSERVED level energies; built by
   !         cooling_data/magnesium_ii_line_cooling.py (cool_MgII_func).
   !   Ca II H&K 3934/3969 A (4s 2S -> 4p 2P),  ground -> levels 4,5: the
   !         coronal resonance-line skeleton. The cooling takes Ca II from
   !         the density-resolved table cool_logLrem_CaII, every channel.
   !   Na I  D doublet 5890/5896 A (3s 2S -> 3p 2P): not in CHIANTI v11, so the
   !         rate uses the rigorous Van Regemorter form with the effective Gaunt
   !         factor gbar=0.2 calibrated against CHIANTI Mg II/Ca II (gbar=0.14-0.24
   !         back-solved over 5000-30000 K); Na I is isoelectronic with Mg II.
   !   Fe II: the (T, n_e) statistical-equilibrium table cool_logL_FeII_ne
   !         (cooling_data/fe2_cooling.py); its low-density edge is the
   !         ground-level coronal rate, validated against Huang (2023)
   !         Fig. 6 red Cloudy coronal curve to ~1.5x.
   !   Fe I  permitted (E1) multiplets: not in CHIANTI v11 (no Fe I model atom),
   !         so - as in Huang (2023) Sec. 2.5 / Fig. 5 - built from NIST ASD
   !         oscillator strengths + the Van Regemorter compact form (Huang
   !         Eq. 10), with the lower levels in Boltzmann equilibrium over the
   !         EVEN-parity metastable manifold (E < 19351 cm^-1, the first odd
   !         level z7D deg). Source: cooling_data/fe1_cooling.py +
   !         fetch_fe1_nist.py (fe1_nist_lines.tsv, fe1_nist_levels.tsv).
   !
   ! Mg II and Na I are evaluated by CLOSED-FORM analytic fits
   ! (cool_MgII_func, cool_NaI_func below; fit scripts
   ! cooling_data/magnesium_ii_line_cooling.py and fit_cooling_formulas.py,
   ! max errors 0.2% and 0.01% over 1e3-1e5 K). Tabulated and interpolated
   ! log-log at runtime on a uniform log10(T) grid (3.0..5.0, dlogT=0.05)
   ! are Fe I (1-D, below), the density-dependent Fe II statistical-
   ! equilibrium coefficient (2-D, cool_logL_FeII_ne) and the
   ! density-resolved C/N/O, Ca II and Mg I tables (2-D, cool_logLrem_*;
   ! see metal_cooling_above_ground_term). The Mg I and Ca II resonance-line
   ! fits remain below as the coronal skeletons cool_MgI_2853_coronal and
   ! cool_CaII_HK_coronal.
   !
   ! Fe II cooling: the density-dependent 2-D table cool_logL_FeII_ne, a
   ! multilevel statistical-equilibrium solve over the full a6D
   ! fine-structure + metastable manifold. A ground-level coronal rate
   ! (Lambda propto n_e, no level saturation) badly overestimates the
   ! cooling at the dense base (n_e ~ 1e8, far above the n_crit ~ 1e4-1e7
   ! of the forbidden a6D IR/metastable lines), where the levels are
   ! collisionally saturated (LTE); at low n_e the table reduces to that
   ! coronal rate. See cool_FeII_ne_value / cool_table_value_2d.
   integer, parameter :: NCOOLT = 41
   real*8,  parameter :: cool_dlogT = 0.05d0
   real*8, parameter :: cool_logT(NCOOLT) = [ &
            3.00000d0,     3.05000d0,     3.10000d0,     3.15000d0,     3.20000d0, &
            3.25000d0,     3.30000d0,     3.35000d0,     3.40000d0,     3.45000d0, &
            3.50000d0,     3.55000d0,     3.60000d0,     3.65000d0,     3.70000d0, &
            3.75000d0,     3.80000d0,     3.85000d0,     3.90000d0,     3.95000d0, &
            4.00000d0,     4.05000d0,     4.10000d0,     4.15000d0,     4.20000d0, &
            4.25000d0,     4.30000d0,     4.35000d0,     4.40000d0,     4.45000d0, &
            4.50000d0,     4.55000d0,     4.60000d0,     4.65000d0,     4.70000d0, &
            4.75000d0,     4.80000d0,     4.85000d0,     4.90000d0,     4.95000d0, &
            5.00000d0 ]
   ! Fe I line cooling, NIST f-values + Van Regemorter compact form
   ! (Huang 2023 Fig. 5 / Sec. 2.5). Lower levels Boltzmann-populated over the
   ! even-parity metastable manifold (17 levels below the first odd level,
   ! z7D deg 19351 cm^-1); 828 permitted lines. Built by cooling_data/
   ! fe1_cooling.py from fe1_nist_{lines,levels}.tsv. Per (n_e n_FeI).
   ! LTE LIMIT OF THE LOWER LEVELS. Boltzmann populations of the metastable
   ! manifold at every density are the n_e -> infinity limit, reached where
   ! collisions outpace the forbidden radiative decays of the metastable
   ! levels. In the n_e -> 0 limit only the a5D4 ground level is populated
   ! and the coefficient is lower by the factor 0.65 (2e3 K), 0.60 (3e3 K),
   ! 0.51 (5e3 K), 0.59 (1e4 K), 0.79 (2e4 K); the a5D term alone gives the
   ! same to 3 per cent (cooling_data/fe1_level_population_limits.py). Over
   ! 2e3-3e4 K the table is therefore the upper limit of the two, at most
   ! 2.0x above the low-density value. The crossover density is not
   ! computed: it needs the forbidden radiative rates and the collision
   ! strengths among the metastable levels, CHIANTI v11.0.2 has no Fe I
   ! model atom, and the Van Regemorter form applies to permitted lines only.
   real*8, parameter :: cool_logL_FeI(NCOOLT) = [ &
            -34.64477d0,     -33.18846d0,     -31.81450d0,     -30.52353d0,     -29.32655d0, &
            -28.23000d0,     -27.23234d0,     -26.32678d0,     -25.50430d0,     -24.75577d0, &
            -24.07318d0,     -23.45016d0,     -22.88187d0,     -22.36457d0,     -21.89509d0, &
            -21.47039d0,     -21.08739d0,     -20.74292d0,     -20.43376d0,     -20.15673d0, &
            -19.90873d0,     -19.68685d0,     -19.48837d0,     -19.31080d0,     -19.15187d0, &
            -19.00953d0,     -18.88195d0,     -18.76748d0,     -18.66467d0,     -18.57222d0, &
            -18.48897d0,     -18.41390d0,     -18.34610d0,     -18.28476d0,     -18.22917d0, &
            -18.17869d0,     -18.13275d0,     -18.09086d0,     -18.05257d0,     -18.01748d0, &
            -17.98524d0 ]
   ! Fe II density-dependent line cooling, multilevel statistical
   ! equilibrium (fe2_cooling.py: lambda_eff_table).  Lambda_eff(T,ne)
   ! [erg cm^3 s^-1] = (sum_u n_u A_ul dE_ul)/ne; coronal at low ne,
   ! ~1/ne (LTE-saturated) at high ne.  logT 3.0..5.0 (NCOOLT=41),
   ! log10 ne 0..14 dlog 0.5 (NCOOLNE=29), clamped at the edges.
   ! Collision strengths: fe_2.scups (5-point knot values) descaled with
   ! the natural cubic spline of CHIANTI's DESCALE_SCUPS.PRO.
   integer, parameter :: NCOOLNE = 29
   real*8, parameter :: cool_logne(NCOOLNE) = [ &
             0.0000d0,      0.5000d0,      1.0000d0,      1.5000d0,      2.0000d0,      2.5000d0, &
             3.0000d0,      3.5000d0,      4.0000d0,      4.5000d0,      5.0000d0,      5.5000d0, &
             6.0000d0,      6.5000d0,      7.0000d0,      7.5000d0,      8.0000d0,      8.5000d0, &
             9.0000d0,      9.5000d0,     10.0000d0,     10.5000d0,     11.0000d0,     11.5000d0, &
            12.0000d0,     12.5000d0,     13.0000d0,     13.5000d0,     14.0000d0 ]
   real*8, parameter :: cool_logL_FeII_ne(NCOOLT,NCOOLNE) = &
      reshape( [ &
            -19.77686d0,     -19.74430d0,     -19.71302d0,     -19.68265d0,     -19.65250d0, &
            -19.62162d0,     -19.58887d0,     -19.55327d0,     -19.51428d0,     -19.47201d0, &
            -19.42729d0,     -19.38140d0,     -19.33583d0,     -19.29190d0,     -19.25052d0, &
            -19.21188d0,     -19.17537d0,     -19.13928d0,     -19.10079d0,     -19.05622d0, &
            -19.00185d0,     -18.93521d0,     -18.85616d0,     -18.76705d0,     -18.67179d0, &
            -18.57458d0,     -18.47898d0,     -18.38750d0,     -18.30171d0,     -18.22243d0, &
            -18.14994d0,     -18.08420d0,     -18.02495d0,     -17.97184d0,     -17.92448d0, &
            -17.88244d0,     -17.84529d0,     -17.81264d0,     -17.78408d0,     -17.75926d0, &
            -17.73785d0,     -19.77717d0,     -19.74461d0,     -19.71334d0,     -19.68296d0, &
            -19.65281d0,     -19.62192d0,     -19.58917d0,     -19.55356d0,     -19.51455d0, &
            -19.47227d0,     -19.42753d0,     -19.38164d0,     -19.33605d0,     -19.29211d0, &
            -19.25071d0,     -19.21207d0,     -19.17555d0,     -19.13945d0,     -19.10095d0, &
            -19.05636d0,     -19.00198d0,     -18.93532d0,     -18.85626d0,     -18.76714d0, &
            -18.67188d0,     -18.57466d0,     -18.47906d0,     -18.38757d0,     -18.30178d0, &
            -18.22249d0,     -18.15000d0,     -18.08426d0,     -18.02500d0,     -17.97190d0, &
            -17.92453d0,     -17.88249d0,     -17.84534d0,     -17.81268d0,     -17.78412d0, &
            -17.75930d0,     -17.73789d0,     -19.77813d0,     -19.74558d0,     -19.71432d0, &
            -19.68394d0,     -19.65379d0,     -19.62288d0,     -19.59010d0,     -19.55446d0, &
            -19.51542d0,     -19.47309d0,     -19.42831d0,     -19.38237d0,     -19.33675d0, &
            -19.29278d0,     -19.25134d0,     -19.21266d0,     -19.17611d0,     -19.13997d0, &
            -19.10143d0,     -19.05680d0,     -19.00238d0,     -18.93568d0,     -18.85659d0, &
            -18.76744d0,     -18.67215d0,     -18.57492d0,     -18.47930d0,     -18.38780d0, &
            -18.30200d0,     -18.22270d0,     -18.15020d0,     -18.08445d0,     -18.02519d0, &
            -17.97207d0,     -17.92470d0,     -17.88265d0,     -17.84549d0,     -17.81282d0, &
            -17.78426d0,     -17.75943d0,     -17.73801d0,     -19.78111d0,     -19.74861d0, &
            -19.71737d0,     -19.68701d0,     -19.65683d0,     -19.62588d0,     -19.59303d0, &
            -19.55729d0,     -19.51812d0,     -19.47567d0,     -19.43075d0,     -19.38468d0, &
            -19.33894d0,     -19.29485d0,     -19.25331d0,     -19.21452d0,     -19.17785d0, &
            -19.14161d0,     -19.10294d0,     -19.05819d0,     -19.00363d0,     -18.93681d0, &
            -18.85761d0,     -18.76837d0,     -18.67301d0,     -18.57572d0,     -18.48005d0, &
            -18.38852d0,     -18.30268d0,     -18.22336d0,     -18.15083d0,     -18.08505d0, &
            -18.02576d0,     -17.97262d0,     -17.92522d0,     -17.88315d0,     -17.84596d0, &
            -17.81327d0,     -17.78468d0,     -17.75983d0,     -17.73839d0,     -19.79015d0, &
            -19.75781d0,     -19.72666d0,     -19.69633d0,     -19.66611d0,     -19.63503d0, &
            -19.60196d0,     -19.56592d0,     -19.52639d0,     -19.48354d0,     -19.43822d0, &
            -19.39176d0,     -19.34564d0,     -19.30121d0,     -19.25933d0,     -19.22023d0, &
            -19.18323d0,     -19.14663d0,     -19.10760d0,     -19.06245d0,     -19.00749d0, &
            -18.94030d0,     -18.86076d0,     -18.77124d0,     -18.67565d0,     -18.57819d0, &
            -18.48238d0,     -18.39073d0,     -18.30480d0,     -18.22538d0,     -18.15277d0, &
            -18.08691d0,     -18.02755d0,     -17.97432d0,     -17.92685d0,     -17.88469d0, &
            -17.84743d0,     -17.81466d0,     -17.78599d0,     -17.76107d0,     -17.73955d0, &
            -19.81621d0,     -19.78428d0,     -19.75340d0,     -19.72315d0,     -19.69281d0, &
            -19.66137d0,     -19.62769d0,     -19.59081d0,     -19.55026d0,     -19.50628d0, &
            -19.45979d0,     -19.41221d0,     -19.36504d0,     -19.31962d0,     -19.27681d0, &
            -19.23679d0,     -19.19886d0,     -19.16128d0,     -19.12118d0,     -19.07491d0, &
            -19.01880d0,     -18.95050d0,     -18.86999d0,     -18.77966d0,     -18.68343d0, &
            -18.58547d0,     -18.48927d0,     -18.39729d0,     -18.31107d0,     -18.23140d0, &
            -18.15855d0,     -18.09246d0,     -18.03286d0,     -17.97941d0,     -17.93170d0, &
            -17.88931d0,     -17.85182d0,     -17.81882d0,     -17.78992d0,     -17.76478d0, &
            -17.74304d0,     -19.88517d0,     -19.85379d0,     -19.82316d0,     -19.79275d0, &
            -19.76177d0,     -19.72913d0,     -19.69366d0,     -19.65443d0,     -19.61109d0, &
            -19.56409d0,     -19.51455d0,     -19.46404d0,     -19.41414d0,     -19.36621d0, &
            -19.32105d0,     -19.27875d0,     -19.23851d0,     -19.19850d0,     -19.15577d0, &
            -19.10668d0,     -19.04768d0,     -18.97662d0,     -18.89365d0,     -18.80129d0, &
            -18.70348d0,     -18.60429d0,     -18.50712d0,     -18.41437d0,     -18.32748d0, &
            -18.24720d0,     -18.17376d0,     -18.10710d0,     -18.04693d0,     -17.99291d0, &
            -17.94462d0,     -17.90164d0,     -17.86357d0,     -17.82998d0,     -17.80051d0, &
            -17.77479d0,     -17.75250d0,     -20.04350d0,     -20.01119d0,     -19.97902d0, &
            -19.94636d0,     -19.91227d0,     -19.87552d0,     -19.83485d0,     -19.78944d0, &
            -19.73919d0,     -19.68495d0,     -19.62824d0,     -19.57093d0,     -19.51478d0, &
            -19.46114d0,     -19.41069d0,     -19.36336d0,     -19.31812d0,     -19.27288d0, &
            -19.22455d0,     -19.16947d0,     -19.10434d0,     -19.02742d0,     -18.93932d0, &
            -18.84282d0,     -18.74186d0,     -18.64034d0,     -18.54145d0,     -18.44735d0, &
            -18.35935d0,     -18.27807d0,     -18.20369d0,     -18.13609d0,     -18.07498d0, &
            -18.01999d0,     -17.97071d0,     -17.92672d0,     -17.88760d0,     -17.85297d0, &
            -17.82244d0,     -17.79567d0,     -17.77234d0,     -20.33079d0,     -20.29447d0, &
            -20.25743d0,     -20.21896d0,     -20.17794d0,     -20.13292d0,     -20.08251d0, &
            -20.02594d0,     -19.96348d0,     -19.89646d0,     -19.82699d0,     -19.75733d0, &
            -19.68946d0,     -19.62474d0,     -19.56372d0,     -19.50603d0,     -19.45029d0, &
            -19.39403d0,     -19.33388d0,     -19.26627d0,     -19.18848d0,     -19.09968d0, &
            -19.00139d0,     -18.89688d0,     -18.79005d0,     -18.68444d0,     -18.58276d0, &
            -18.48676d0,     -18.39742d0,     -18.31512d0,     -18.23989d0,     -18.17153d0, &
            -18.10969d0,     -18.05396d0,     -18.00392d0,     -17.95913d0,     -17.91916d0, &
            -17.88363d0,     -17.85215d0,     -17.82438d0,     -17.80001d0,     -20.72909d0, &
            -20.68830d0,     -20.64612d0,     -20.60179d0,     -20.55400d0,     -20.50096d0, &
            -20.44101d0,     -20.37320d0,     -20.29790d0,     -20.21679d0,     -20.13236d0, &
            -20.04720d0,     -19.96347d0,     -19.88252d0,     -19.80470d0,     -19.72930d0, &
            -19.65441d0,     -19.57706d0,     -19.49380d0,     -19.40168d0,     -19.29950d0, &
            -19.18837d0,     -19.07136d0,     -18.95236d0,     -18.83503d0,     -18.72220d0, &
            -18.61575d0,     -18.51671d0,     -18.42551d0,     -18.34214d0,     -18.26637d0, &
            -18.19782d0,     -18.13603d0,     -18.08048d0,     -18.03070d0,     -17.98617d0, &
            -17.94645d0,     -17.91110d0,     -17.87972d0,     -17.85194d0,     -17.82746d0, &
            -21.18928d0,     -21.14598d0,     -21.10099d0,     -21.05352d0,     -21.00205d0, &
            -20.94452d0,     -20.87885d0,     -20.80373d0,     -20.71923d0,     -20.62691d0, &
            -20.52918d0,     -20.42862d0,     -20.32726d0,     -20.22624d0,     -20.12557d0, &
            -20.02399d0,     -19.91909d0,     -19.80782d0,     -19.68769d0,     -19.55803d0, &
            -19.42065d0,     -19.27923d0,     -19.13813d0,     -19.00113d0,     -18.87100d0, &
            -18.74940d0,     -18.63713d0,     -18.53438d0,     -18.44094d0,     -18.35638d0, &
            -18.28018d0,     -18.21172d0,     -18.15040d0,     -18.09562d0,     -18.04679d0, &
            -18.00336d0,     -17.96479d0,     -17.93061d0,     -17.90037d0,     -17.87368d0, &
            -17.85018d0,     -21.67551d0,     -21.63122d0,     -21.58515d0,     -21.53644d0, &
            -21.48349d0,     -21.42401d0,     -21.35555d0,     -21.27630d0,     -21.18569d0, &
            -21.08448d0,     -20.97425d0,     -20.85669d0,     -20.73312d0,     -20.60423d0, &
            -20.46989d0,     -20.32911d0,     -20.18043d0,     -20.02289d0,     -19.85719d0, &
            -19.68620d0,     -19.51425d0,     -19.34590d0,     -19.18490d0,     -19.03377d0, &
            -18.89383d0,     -18.76556d0,     -18.64884d0,     -18.54319d0,     -18.44795d0, &
            -18.36240d0,     -18.28578d0,     -18.21734d0,     -18.15637d0,     -18.10219d0, &
            -18.05416d0,     -18.01168d0,     -17.97418d0,     -17.94115d0,     -17.91211d0, &
            -17.88665d0,     -17.86439d0,     -22.17101d0,     -22.12639d0,     -22.07994d0, &
            -22.03081d0,     -21.97733d0,     -21.91706d0,     -21.84724d0,     -21.76542d0, &
            -21.66980d0,     -21.55913d0,     -21.43247d0,     -21.28926d0,     -21.13002d0, &
            -20.95665d0,     -20.77184d0,     -20.57815d0,     -20.37780d0,     -20.17320d0, &
            -19.96761d0,     -19.76503d0,     -19.56945d0,     -19.38413d0,     -19.21125d0, &
            -19.05189d0,     -18.90629d0,     -18.77413d0,     -18.65475d0,     -18.54730d0, &
            -18.45088d0,     -18.36459d0,     -18.28758d0,     -18.21902d0,     -18.15814d0, &
            -18.10421d0,     -18.05656d0,     -18.01456d0,     -17.97763d0,     -17.94522d0, &
            -17.91687d0,     -17.89213d0,     -17.87063d0,     -22.66958d0,     -22.62485d0, &
            -22.57828d0,     -22.52900d0,     -22.47533d0,     -22.41470d0,     -22.34399d0, &
            -22.25965d0,     -22.15736d0,     -22.03153d0,     -21.87633d0,     -21.68928d0, &
            -21.47452d0,     -21.24123d0,     -20.99885d0,     -20.75396d0,     -20.51042d0, &
            -20.27101d0,     -20.03850d0,     -19.81586d0,     -19.60573d0,     -19.41001d0, &
            -19.22971d0,     -19.06503d0,     -18.91557d0,     -18.78058d0,     -18.65910d0, &
            -18.55010d0,     -18.45257d0,     -18.36551d0,     -18.28799d0,     -18.21915d0, &
            -18.15816d0,     -18.10426d0,     -18.05675d0,     -18.01498d0,     -17.97833d0, &
            -17.94626d0,     -17.91827d0,     -17.89393d0,     -17.87282d0,     -23.16912d0, &
            -23.12436d0,     -23.07775d0,     -23.02843d0,     -22.97467d0,     -22.91382d0, &
            -22.84222d0,     -22.75468d0,     -22.64279d0,     -22.49428d0,     -22.29826d0, &
            -22.05550d0,     -21.78063d0,     -21.49172d0,     -21.20135d0,     -20.91573d0, &
            -20.63776d0,     -20.36967d0,     -20.11390d0,     -19.87291d0,     -19.64862d0, &
            -19.44209d0,     -19.25355d0,     -19.08255d0,     -18.92821d0,     -18.78942d0, &
            -18.66500d0,     -18.55375d0,     -18.45451d0,     -18.36622d0,     -18.28786d0, &
            -18.21849d0,     -18.15721d0,     -18.10322d0,     -18.05577d0,     -18.01414d0, &
            -17.97771d0,     -17.94591d0,     -17.91822d0,     -17.89418d0,     -17.87340d0, &
            -23.66898d0,     -23.62420d0,     -23.57758d0,     -23.52824d0,     -23.47445d0, &
            -23.41344d0,     -23.34112d0,     -23.25066d0,     -23.12959d0,     -22.95942d0, &
            -22.72655d0,     -22.43776d0,     -22.11601d0,     -21.78239d0,     -21.44883d0, &
            -21.12117d0,     -20.80356d0,     -20.50007d0,     -20.21445d0,     -19.94938d0, &
            -19.70618d0,     -19.48503d0,     -19.28523d0,     -19.10554d0,     -18.94451d0, &
            -18.80058d0,     -18.67224d0,     -18.55805d0,     -18.45665d0,     -18.36683d0, &
            -18.28743d0,     -18.21740d0,     -18.15577d0,     -18.10163d0,     -18.05417d0, &
            -18.01266d0,     -17.97641d0,     -17.94482d0,     -17.91738d0,     -17.89360d0, &
            -17.87309d0,     -24.16893d0,     -24.12415d0,     -24.07753d0,     -24.02819d0, &
            -23.97438d0,     -23.91329d0,     -23.84058d0,     -23.74847d0,     -23.62196d0, &
            -23.43851d0,     -23.18189d0,     -22.85983d0,     -22.49634d0,     -22.11282d0, &
            -21.72371d0,     -21.34037d0,     -20.97250d0,     -20.62733d0,     -20.30901d0, &
            -20.01908d0,     -19.75726d0,     -19.52218d0,     -19.31194d0,     -19.12444d0, &
            -18.95756d0,     -18.80929d0,     -18.67776d0,     -18.56125d0,     -18.45822d0, &
            -18.36727d0,     -18.28713d0,     -18.21665d0,     -18.15476d0,     -18.10051d0, &
            -18.05304d0,     -18.01158d0,     -17.97543d0,     -17.94397d0,     -17.91666d0, &
            -17.89303d0,     -17.87266d0,     -24.66892d0,     -24.62414d0,     -24.57751d0, &
            -24.52817d0,     -24.47435d0,     -24.41324d0,     -24.34037d0,     -24.24754d0, &
            -24.11834d0,     -23.92642d0,     -23.64832d0,     -23.28409d0,     -22.85699d0, &
            -22.39837d0,     -21.93708d0,     -21.49414d0,     -21.08160d0,     -20.70459d0, &
            -20.36399d0,     -20.05842d0,     -19.78550d0,     -19.54242d0,     -19.32637d0, &
            -19.13461d0,     -18.96462d0,     -18.81408d0,     -18.68092d0,     -18.56326d0, &
            -18.45944d0,     -18.36796d0,     -18.28747d0,     -18.21677d0,     -18.15477d0, &
            -18.10046d0,     -18.05298d0,     -18.01153d0,     -17.97541d0,     -17.94399d0, &
            -17.91673d0,     -17.89314d0,     -17.87282d0,     -25.16891d0,     -25.12413d0, &
            -25.07751d0,     -25.02816d0,     -24.97435d0,     -24.91322d0,     -24.84028d0, &
            -24.74705d0,     -24.61536d0,     -24.41076d0,     -24.09048d0,     -23.64356d0, &
            -23.11648d0,     -22.57183d0,     -22.04973d0,     -21.56763d0,     -21.13062d0, &
            -20.73825d0,     -20.38777d0,     -20.07566d0,     -19.79824d0,     -19.55198d0, &
            -19.33360d0,     -19.14014d0,     -18.96888d0,     -18.81743d0,     -18.68361d0, &
            -18.56549d0,     -18.46134d0,     -18.36966d0,     -18.28904d0,     -18.21827d0, &
            -18.15622d0,     -18.10189d0,     -18.05439d0,     -18.01293d0,     -17.97680d0, &
            -17.94536d0,     -17.91809d0,     -17.89448d0,     -17.87414d0,     -25.66891d0, &
            -25.62413d0,     -25.57750d0,     -25.52816d0,     -25.47434d0,     -25.41321d0, &
            -25.34022d0,     -25.24636d0,     -25.10892d0,     -24.87136d0,     -24.45523d0, &
            -23.87841d0,     -23.25163d0,     -22.64870d0,     -22.09515d0,     -21.59606d0, &
            -21.14962d0,     -20.75183d0,     -20.39812d0,     -20.08400d0,     -19.80524d0, &
            -19.55804d0,     -19.33896d0,     -19.14495d0,     -18.97327d0,     -18.82151d0, &
            -18.68747d0,     -18.56923d0,     -18.46503d0,     -18.37335d0,     -18.29278d0, &
            -18.22206d0,     -18.16008d0,     -18.10582d0,     -18.05837d0,     -18.01694d0, &
            -17.98082d0,     -17.94938d0,     -17.92207d0,     -17.89842d0,     -17.87802d0, &
            -26.16891d0,     -26.12413d0,     -26.07750d0,     -26.02816d0,     -25.97434d0, &
            -25.91320d0,     -25.84009d0,     -25.74446d0,     -25.59006d0,     -25.26786d0, &
            -24.69065d0,     -23.99130d0,     -23.30574d0,     -22.67688d0,     -22.11136d0, &
            -21.60641d0,     -21.15699d0,     -20.75768d0,     -20.40321d0,     -20.08872d0, &
            -19.80980d0,     -19.56253d0,     -19.34342d0,     -19.14940d0,     -18.97774d0, &
            -18.82604d0,     -18.69212d0,     -18.57406d0,     -18.47011d0,     -18.37873d0, &
            -18.29851d0,     -18.22818d0,     -18.16658d0,     -18.11269d0,     -18.06560d0, &
            -18.02449d0,     -17.98864d0,     -17.95742d0,     -17.93028d0,     -17.90676d0, &
            -17.88643d0,     -26.66891d0,     -26.62413d0,     -26.57750d0,     -26.52816d0, &
            -26.47434d0,     -26.41319d0,     -26.33971d0,     -26.23857d0,     -26.03555d0, &
            -25.54537d0,     -24.80303d0,     -24.03430d0,     -23.32453d0,     -22.68638d0, &
            -22.11686d0,     -21.61007d0,     -21.15978d0,     -20.76009d0,     -20.40551d0, &
            -20.09105d0,     -19.81223d0,     -19.56508d0,     -19.34610d0,     -19.15222d0, &
            -18.98073d0,     -18.82922d0,     -18.69556d0,     -18.57781d0,     -18.47426d0, &
            -18.38336d0,     -18.30369d0,     -18.23398d0,     -18.17306d0,     -18.11989d0, &
            -18.07354d0,     -18.03317d0,     -17.99805d0,     -17.96754d0,     -17.94108d0, &
            -17.91818d0,     -17.89843d0,     -27.16891d0,     -27.12413d0,     -27.07750d0, &
            -27.02816d0,     -26.97434d0,     -26.91314d0,     -26.83851d0,     -26.72051d0, &
            -26.39798d0,     -25.68987d0,     -24.84570d0,     -24.04891d0,     -23.33072d0, &
            -22.68951d0,     -22.11871d0,     -21.61133d0,     -21.16079d0,     -20.76100d0, &
            -20.40642d0,     -20.09199d0,     -19.81323d0,     -19.56616d0,     -19.34727d0, &
            -19.15348d0,     -18.98209d0,     -18.83072d0,     -18.69722d0,     -18.57968d0, &
            -18.47640d0,     -18.38583d0,     -18.30655d0,     -18.23731d0,     -18.17693d0, &
            -18.12436d0,     -18.07868d0,     -18.03903d0,     -18.00469d0,     -17.97499d0, &
            -17.94938d0,     -17.92735d0,     -17.90849d0,     -27.66891d0,     -27.62413d0, &
            -27.57750d0,     -27.52816d0,     -27.47434d0,     -27.41297d0,     -27.33473d0, &
            -27.16792d0,     -26.62895d0,     -25.74798d0,     -24.86030d0,     -24.05380d0, &
            -23.33285d0,     -22.69065d0,     -22.11943d0,     -21.61186d0,     -21.16123d0, &
            -20.76141d0,     -20.40682d0,     -20.09241d0,     -19.81367d0,     -19.56662d0, &
            -19.34776d0,     -19.15401d0,     -18.98266d0,     -18.83133d0,     -18.69790d0, &
            -18.58045d0,     -18.47729d0,     -18.38685d0,     -18.30776d0,     -18.23873d0, &
            -18.17860d0,     -18.12633d0,     -18.08098d0,     -18.04172d0,     -18.00781d0, &
            -17.97858d0,     -17.95349d0,     -17.93202d0,     -17.91377d0,     -28.16891d0, &
            -28.12413d0,     -28.07750d0,     -28.02816d0,     -27.97432d0,     -27.91246d0, &
            -27.82302d0,     -27.53446d0,     -26.73855d0,     -25.76870d0,     -24.86559d0, &
            -24.05589d0,     -23.33402d0,     -22.69147d0,     -22.12009d0,     -21.61243d0, &
            -21.16174d0,     -20.76189d0,     -20.40727d0,     -20.09284d0,     -19.81409d0, &
            -19.56703d0,     -19.34817d0,     -19.15441d0,     -18.98307d0,     -18.83175d0, &
            -18.69834d0,     -18.58092d0,     -18.47779d0,     -18.38741d0,     -18.30837d0, &
            -18.23942d0,     -18.17939d0,     -18.12724d0,     -18.08202d0,     -18.04291d0, &
            -18.00916d0,     -17.98013d0,     -17.95525d0,     -17.93402d0,     -17.91603d0, &
            -28.66891d0,     -28.62413d0,     -28.57750d0,     -28.52816d0,     -28.47428d0, &
            -28.41086d0,     -28.28822d0,     -27.77193d0,     -26.78188d0,     -25.77740d0, &
            -24.86908d0,     -24.05823d0,     -23.33595d0,     -22.69318d0,     -22.12164d0, &
            -21.61386d0,     -21.16307d0,     -20.76312d0,     -20.40842d0,     -20.09392d0, &
            -19.81510d0,     -19.56799d0,     -19.34907d0,     -19.15527d0,     -18.98390d0, &
            -18.83256d0,     -18.69913d0,     -18.58171d0,     -18.47859d0,     -18.38822d0, &
            -18.30922d0,     -18.24030d0,     -18.18031d0,     -18.12820d0,     -18.08305d0, &
            -18.04400d0,     -18.01032d0,     -17.98137d0,     -17.95657d0,     -17.93543d0, &
            -17.91754d0,     -29.16891d0,     -29.12413d0,     -29.07750d0,     -29.02816d0, &
            -28.97414d0,     -28.90598d0,     -28.69573d0,     -27.89158d0,     -26.80300d0, &
            -25.78627d0,     -24.87585d0,     -24.06422d0,     -23.34143d0,     -22.69824d0, &
            -22.12633d0,     -21.61822d0,     -21.16713d0,     -20.76690d0,     -20.41195d0, &
            -20.09721d0,     -19.81818d0,     -19.57088d0,     -19.35179d0,     -19.15786d0, &
            -18.98637d0,     -18.83494d0,     -18.70146d0,     -18.58400d0,     -18.48086d0, &
            -18.39050d0,     -18.31152d0,     -18.24264d0,     -18.18269d0,     -18.13064d0, &
            -18.08553d0,     -18.04654d0,     -18.01292d0,     -17.98402d0,     -17.95928d0, &
            -17.93820d0,     -17.92035d0,     -29.66891d0,     -29.62413d0,     -29.57750d0, &
            -29.52815d0,     -29.47375d0,     -29.39207d0,     -29.00230d0,     -27.95774d0, &
            -26.82971d0,     -25.80766d0,     -24.89528d0,     -24.08221d0,     -23.35814d0, &
            -22.71377d0,     -22.14079d0,     -21.63169d0,     -21.17968d0,     -20.77861d0, &
            -20.42287d0,     -20.10742d0,     -19.82774d0,     -19.57985d0,     -19.36025d0, &
            -19.16588d0,     -18.99403d0,     -18.84233d0,     -18.70865d0,     -18.59108d0, &
            -18.48789d0,     -18.39753d0,     -18.31860d0,     -18.24980d0,     -18.18996d0, &
            -18.13803d0,     -18.09306d0,     -18.05420d0,     -18.02071d0,     -17.99193d0, &
            -17.96729d0,     -17.94629d0,     -17.92852d0,     -30.16891d0,     -30.12413d0, &
            -30.07750d0,     -30.02813d0,     -29.97278d0,     -29.85917d0,     -29.22134d0, &
            -28.03937d0,     -26.89438d0,     -25.86719d0,     -24.95088d0,     -24.13421d0, &
            -23.40677d0,     -22.75926d0,     -22.18335d0,     -21.67152d0,     -21.21697d0, &
            -20.81353d0,     -20.45558d0,     -20.13808d0,     -19.85652d0,     -19.60693d0, &
            -19.38584d0,     -19.19019d0,     -19.01730d0,     -18.86480d0,     -18.73054d0, &
            -18.61261d0,     -18.50926d0,     -18.41891d0,     -18.34011d0,     -18.27155d0, &
            -18.21203d0,     -18.16046d0,     -18.11586d0,     -18.07737d0,     -18.04424d0, &
            -18.01578d0,     -17.99143d0,     -17.97068d0,     -17.95310d0 ], &
      [NCOOLT, NCOOLNE] )

   ! van Hoof et al. (2014, MNRAS 444, 420) thermally-averaged free-free Gaunt
   ! factor <g_ff>(gamma^2), gamma^2 = Z_ion^2 Ry/(k T_e) (their eq. 21), Ry
   ! the infinite-mass Rydberg unit of energy (Ry_over_kB below).  The
   ! frequency average <g_ff> = Int e^-u g_ff du / Int e^-u du depends only on
   ! gamma^2; gbar_ff() interpolates this table in log10(gamma^2).
   ! Generated by cooling_data/gauntff_thermal_avg.py.
   integer, parameter :: NGFF = 161
   ! Ry/k [K] of gamma^2. van Hoof et al. (2014), section 2: "Ry is the
   ! infinite-mass Rydberg unit of energy given by 1Ry = alpha^2 m_e c^2/2
   ! ~ 2.179 87 x 10^-18 J". It is h c R_inf, with R_inf the CODATA 2018
   ! Rydberg constant 109737.31568160 cm^-1 and h, c, k_B the exact SI
   ! values the code carries (global_parameters): Ry = 2.1798724e-11 erg,
   ! Ry/k = 157887.51 K.
   real*8, parameter :: R_inf_cm   = 109737.31568160d0
   real*8, parameter :: Ry_over_kB = hp_erg*c_light*R_inf_cm/kb_erg
   real*8, parameter :: gff_lg2_min  = -6.0000d0   ! log10(gamma^2) start
   real*8, parameter :: gff_lg2_step = 0.1000d0   ! step [dex]
   real*8, parameter :: gff_avg(NGFF) = [ &
        1.105130d0, 1.105283d0, 1.105437d0, 1.105630d0, 1.105823d0, &
        1.106067d0, 1.106310d0, 1.106618d0, 1.106926d0, 1.107312d0, &
        1.107699d0, 1.108186d0, 1.108672d0, 1.109286d0, 1.109899d0, &
        1.110672d0, 1.111446d0, 1.112421d0, 1.113397d0, 1.114626d0, &
        1.115855d0, 1.117404d0, 1.118953d0, 1.120903d0, 1.122854d0, &
        1.125310d0, 1.127766d0, 1.130853d0, 1.133941d0, 1.137811d0, &
        1.141682d0, 1.146519d0, 1.151356d0, 1.157366d0, 1.163376d0, &
        1.170784d0, 1.178193d0, 1.187223d0, 1.196253d0, 1.207091d0, &
        1.217929d0, 1.230670d0, 1.243410d0, 1.257973d0, 1.272535d0, &
        1.288586d0, 1.304637d0, 1.321487d0, 1.338336d0, 1.354915d0, &
        1.371494d0, 1.386395d0, 1.401296d0, 1.412960d0, 1.424623d0, &
        1.431646d0, 1.438670d0, 1.440157d0, 1.441645d0, 1.437481d0, &
        1.433317d0, 1.424201d0, 1.415086d0, 1.402301d0, 1.389517d0, &
        1.374538d0, 1.359560d0, 1.343691d0, 1.327821d0, 1.312005d0, &
        1.296188d0, 1.281006d0, 1.265823d0, 1.251586d0, 1.237348d0, &
        1.224198d0, 1.211048d0, 1.199033d0, 1.187017d0, 1.176126d0, &
        1.165235d0, 1.155425d0, 1.145616d0, 1.136824d0, 1.128033d0, &
        1.120188d0, 1.112344d0, 1.105371d0, 1.098398d0, 1.092218d0, &
        1.086038d0, 1.080578d0, 1.075118d0, 1.070306d0, 1.065493d0, &
        1.061261d0, 1.057028d0, 1.053314d0, 1.049599d0, 1.046345d0, &
        1.043091d0, 1.040246d0, 1.037401d0, 1.034916d0, 1.032432d0, &
        1.030265d0, 1.028098d0, 1.026213d0, 1.024327d0, 1.022686d0, &
        1.021045d0, 1.019619d0, 1.018194d0, 1.016957d0, 1.015720d0, &
        1.014647d0, 1.013574d0, 1.012644d0, 1.011715d0, 1.010909d0, &
        1.010104d0, 1.009408d0, 1.008713d0, 1.008110d0, 1.007506d0, &
        1.006987d0, 1.006467d0, 1.006019d0, 1.005571d0, 1.005182d0, &
        1.004793d0, 1.004459d0, 1.004125d0, 1.003837d0, 1.003549d0, &
        1.003300d0, 1.003051d0, 1.002837d0, 1.002622d0, 1.002439d0, &
        1.002255d0, 1.002096d0, 1.001937d0, 1.001800d0, 1.001664d0, &
        1.001545d0, 1.001427d0, 1.001326d0, 1.001225d0, 1.001138d0, &
        1.001051d0, 1.000975d0, 1.000900d0, 1.000836d0, 1.000772d0, &
        1.000717d0, 1.000661d0, 1.000614d0, 1.000567d0, 1.000526d0, &
        1.000485d0 ]

   !--- Density-resolved line cooling above the ground term ------------!
   ! log10 Lambda_rem(T, n_e), the statistical-equilibrium line cooling of
   ! one ion divided by n_e, summed over every transition whose upper
   ! level lies outside the ground term. The axes are cool_logT and
   ! cool_logne above; what the runtime does with it is written at
   ! metal_cooling_above_ground_term. Level populations: the multilevel
   ! statistical equilibrium of cooling_data/multilevel_statistical_
   ! equilibrium.py, collision strengths descaled with the natural cubic
   ! spline of CHIANTI's DESCALE_SCUPS.PRO. Valid over 1e3-1e5 K and
   ! 1-1e14 cm^-3 (the axes; clamped at the edges), electron collisions
   ! only, optically thin.
   ! GENERATED by cooling_data/metal_cooling_density_resolved.py
   ! from CHIANTI v11.0.2. Do not edit by hand.

   ! CI
   real*8, parameter :: cool_logLrem_CI(NCOOLT,NCOOLNE) = reshape( [ &
            -26.202357d0,     -25.483057d0,     -24.838716d0,     -24.261378d0,     -23.744003d0,     -23.280368d0, &
            -22.864983d0,     -22.493016d0,     -22.160231d0,     -21.862807d0,     -21.596983d0,     -21.359309d0, &
            -21.146670d0,     -20.956232d0,     -20.785385d0,     -20.631720d0,     -20.493063d0,     -20.367355d0, &
            -20.252612d0,     -20.146920d0,     -20.048443d0,     -19.955429d0,     -19.866237d0,     -19.779411d0, &
            -19.693756d0,     -19.608418d0,     -19.522958d0,     -19.437404d0,     -19.352196d0,     -19.267962d0, &
            -19.185428d0,     -19.105174d0,     -19.027648d0,     -18.953374d0,     -18.882842d0,     -18.816258d0, &
            -18.753556d0,     -18.694752d0,     -18.639934d0,     -18.589187d0,     -18.542553d0,     -26.199096d0, &
            -25.479962d0,     -24.835793d0,     -24.258632d0,     -23.741436d0,     -23.277973d0,     -22.862748d0, &
            -22.490923d0,     -22.158258d0,     -21.860930d0,     -21.595188d0,     -21.357587d0,     -21.145024d0, &
            -20.954677d0,     -20.783945d0,     -20.630426d0,     -20.491933d0,     -20.366392d0,     -20.251801d0, &
            -20.146230d0,     -20.047825d0,     -19.954823d0,     -19.865589d0,     -19.778687d0,     -19.692941d0, &
            -19.607518d0,     -19.521997d0,     -19.436415d0,     -19.351193d0,     -19.266939d0,     -19.184366d0, &
            -19.104067d0,     -19.026531d0,     -18.952296d0,     -18.881807d0,     -18.815239d0,     -18.752553d0, &
            -18.693791d0,     -18.639052d0,     -18.588420d0,     -18.541930d0,     -26.194806d0,     -25.476140d0, &
            -24.832406d0,     -24.255646d0,     -23.738808d0,     -23.275662d0,     -22.860710d0,     -22.489120d0, &
            -22.156653d0,     -21.859493d0,     -21.593894d0,     -21.356422d0,     -21.143978d0,     -20.953746d0, &
            -20.783133d0,     -20.629738d0,     -20.491369d0,     -20.365945d0,     -20.251453d0,     -20.145954d0, &
            -20.047583d0,     -19.954567d0,     -19.865271d0,     -19.778266d0,     -19.692391d0,     -19.606831d0, &
            -19.521184d0,     -19.435505d0,     -19.350223d0,     -19.265938d0,     -19.183351d0,     -19.103045d0, &
            -19.025508d0,     -18.951294d0,     -18.880836d0,     -18.814294d0,     -18.751645d0,     -18.692940d0, &
            -18.638283d0,     -18.587757d0,     -18.541390d0,     -26.191532d0,     -25.473466d0,     -24.830258d0, &
            -24.253952d0,     -23.737500d0,     -23.274680d0,     -22.860003d0,     -22.488644d0,     -22.156372d0, &
            -21.859377d0,     -21.593921d0,     -21.356576d0,     -21.144246d0,     -20.954121d0,     -20.783607d0, &
            -20.630308d0,     -20.492028d0,     -20.366679d0,     -20.252247d0,     -20.146785d0,     -20.048420d0, &
            -19.955373d0,     -19.866003d0,     -19.778883d0,     -19.692858d0,     -19.607125d0,     -19.521299d0, &
            -19.435453d0,     -19.350035d0,     -19.265653d0,     -19.183003d0,     -19.102652d0,     -19.025079d0, &
            -18.950840d0,     -18.880373d0,     -18.813838d0,     -18.751209d0,     -18.692541d0,     -18.637934d0, &
            -18.587470d0,     -18.541171d0,     -26.191009d0,     -25.473442d0,     -24.830668d0,     -24.254737d0, &
            -23.738614d0,     -23.276083d0,     -22.861664d0,     -22.490538d0,     -22.158478d0,     -21.861677d0, &
            -21.596400d0,     -21.359222d0,     -21.147051d0,     -20.957076d0,     -20.786707d0,     -20.633542d0, &
            -20.495383d0,     -20.370137d0,     -20.255780d0,     -20.150355d0,     -20.051981d0,     -19.958864d0, &
            -19.869358d0,     -19.782030d0,     -19.695730d0,     -19.609666d0,     -19.523477d0,     -19.437260d0, &
            -19.351491d0,     -19.266796d0,     -19.183882d0,     -19.103317d0,     -19.025576d0,     -18.951211d0, &
            -18.880657d0,     -18.814067d0,     -18.751409d0,     -18.692733d0,     -18.638135d0,     -18.587689d0, &
            -18.541415d0,     -26.195119d0,     -25.478046d0,     -24.835753d0,     -24.260293d0,     -23.744635d0, &
            -23.282568d0,     -22.868613d0,     -22.497946d0,     -22.166336d0,     -21.869971d0,     -21.605116d0, &
            -21.368347d0,     -21.156573d0,     -20.966982d0,     -20.796980d0,     -20.644160d0,     -20.506312d0, &
            -20.381323d0,     -20.267147d0,     -20.161801d0,     -20.063375d0,     -19.970046d0,     -19.880144d0, &
            -19.792226d0,     -19.705148d0,     -19.618151d0,     -19.530926d0,     -19.443639d0,     -19.356824d0, &
            -19.271164d0,     -19.187407d0,     -19.106146d0,     -19.027864d0,     -18.953102d0,     -18.882260d0, &
            -18.815467d0,     -18.752676d0,     -18.693921d0,     -18.639283d0,     -18.588822d0,     -18.542548d0, &
            -26.210170d0,     -25.494106d0,     -24.852876d0,     -24.278534d0,     -23.764046d0,     -23.303191d0, &
            -22.890479d0,     -22.521071d0,     -22.190711d0,     -21.895565d0,     -21.631894d0,     -21.396274d0, &
            -21.185613d0,     -20.997093d0,     -20.828110d0,     -20.676235d0,     -20.539217d0,     -20.414891d0, &
            -20.301142d0,     -20.195910d0,     -20.097197d0,     -20.003101d0,     -19.911890d0,     -19.822095d0, &
            -19.732609d0,     -19.642776d0,     -19.552453d0,     -19.462003d0,     -19.372135d0,     -19.283682d0, &
            -19.197500d0,     -19.114255d0,     -19.034446d0,     -18.958566d0,     -18.886916d0,     -18.819552d0, &
            -18.756384d0,     -18.697401d0,     -18.642639d0,     -18.592123d0,     -18.545834d0,     -26.255323d0, &
            -25.541862d0,     -24.903399d0,     -24.331977d0,     -23.820539d0,     -23.362835d0,     -22.953336d0, &
            -22.587154d0,     -22.259969d0,     -21.967889d0,     -21.707162d0,     -21.474360d0,     -21.266381d0, &
            -21.080388d0,     -20.913736d0,     -20.763926d0,     -20.628588d0,     -20.505399d0,     -20.392055d0, &
            -20.286282d0,     -20.185865d0,     -20.088708d0,     -19.992968d0,     -19.897197d0,     -19.800479d0, &
            -19.702537d0,     -19.603736d0,     -19.504967d0,     -19.407364d0,     -19.312068d0,     -19.220120d0, &
            -19.132261d0,     -19.048963d0,     -18.970560d0,     -18.897116d0,     -18.828508d0,     -18.764537d0, &
            -18.705084d0,     -18.650090d0,     -18.599498d0,     -18.553231d0,     -26.373525d0,     -25.665620d0, &
            -25.032965d0,     -24.467554d0,     -23.962274d0,     -23.510801d0,     -23.107519d0,     -22.747437d0, &
            -22.426117d0,     -22.139552d0,     -21.883962d0,     -21.655894d0,     -21.452199d0,     -21.269954d0, &
            -21.106376d0,     -20.958740d0,     -20.824353d0,     -20.700487d0,     -20.584393d0,     -20.473362d0, &
            -20.364823d0,     -20.256506d0,     -20.146702d0,     -20.034478d0,     -19.919798d0,     -19.803495d0, &
            -19.687059d0,     -19.572276d0,     -19.460791d0,     -19.353953d0,     -19.252771d0,     -19.157820d0, &
            -19.069309d0,     -18.987220d0,     -18.911218d0,     -18.840894d0,     -18.775859d0,     -18.715834d0, &
            -18.660617d0,     -18.610044d0,     -18.563951d0,     -26.617462d0,     -25.917067d0,     -25.292028d0, &
            -24.734260d0,     -24.236562d0,     -23.792518d0,     -23.396403d0,     -23.043100d0,     -22.728022d0, &
            -22.447008d0,     -22.196177d0,     -21.971944d0,     -21.770948d0,     -21.589950d0,     -21.425713d0, &
            -21.274907d0,     -21.134098d0,     -20.999821d0,     -20.868761d0,     -20.738015d0,     -20.605346d0, &
            -20.469410d0,     -20.329983d0,     -20.187914d0,     -20.044897d0,     -19.903121d0,     -19.764878d0, &
            -19.632186d0,     -19.506462d0,     -19.388611d0,     -19.279120d0,     -19.178065d0,     -19.085205d0, &
            -19.000110d0,     -18.922075d0,     -18.850430d0,     -18.784614d0,     -18.724209d0,     -18.668905d0, &
            -18.618450d0,     -18.572620d0,     -26.990193d0,     -26.295453d0,     -25.675974d0,     -25.123618d0, &
            -24.631122d0,     -24.191999d0,     -23.800425d0,     -23.451140d0,     -23.139331d0,     -22.860522d0, &
            -22.610446d0,     -22.384981d0,     -22.180055d0,     -21.991538d0,     -21.815159d0,     -21.646546d0, &
            -21.481555d0,     -21.316767d0,     -21.149955d0,     -20.980258d0,     -20.807964d0,     -20.634169d0, &
            -20.460571d0,     -20.289156d0,     -20.121981d0,     -19.961010d0,     -19.807946d0,     -19.664061d0, &
            -19.530019d0,     -19.406069d0,     -19.292163d0,     -19.187962d0,     -19.092905d0,     -19.006316d0, &
            -18.927285d0,     -18.855006d0,     -18.788826d0,     -18.728258d0,     -18.672936d0,     -18.622570d0, &
            -18.576909d0,     -27.440716d0,     -26.748634d0,     -26.131720d0,     -25.581803d0,     -25.091575d0, &
            -24.654461d0,     -24.264458d0,     -23.915966d0,     -23.603587d0,     -23.321957d0,     -23.065614d0, &
            -22.828994d0,     -22.606524d0,     -22.392839d0,     -22.183144d0,     -21.973761d0,     -21.762793d0, &
            -21.550285d0,     -21.337714d0,     -21.127091d0,     -20.920222d0,     -20.718477d0,     -20.523045d0, &
            -20.335004d0,     -20.155408d0,     -19.985290d0,     -19.825561d0,     -19.676860d0,     -19.539355d0, &
            -19.412927d0,     -19.297259d0,     -19.191817d0,     -19.095900d0,     -19.008727d0,     -18.929308d0, &
            -18.856782d0,     -18.790460d0,     -18.729827d0,     -18.674496d0,     -18.624165d0,     -18.578568d0, &
            -27.923819d0,     -27.232711d0,     -26.616719d0,     -26.067638d0,     -25.578073d0,     -25.141221d0, &
            -24.750582d0,     -24.399581d0,     -24.081188d0,     -23.787747d0,     -23.511234d0,     -23.244122d0, &
            -22.980583d0,     -22.717345d0,     -22.453711d0,     -22.190907d0,     -21.931374d0,     -21.677839d0, &
            -21.432547d0,     -21.196867d0,     -20.971278d0,     -20.755733d0,     -20.550196d0,     -20.354781d0, &
            -20.169815d0,     -19.995794d0,     -19.833237d0,     -19.682492d0,     -19.543514d0,     -19.416027d0, &
            -19.299600d0,     -19.193618d0,     -19.097320d0,     -19.009882d0,     -18.930281d0,     -18.857635d0, &
            -18.791238d0,     -18.730560d0,     -18.675208d0,     -18.624872d0,     -18.579283d0,     -28.418335d0, &
            -27.727544d0,     -27.111827d0,     -26.562901d0,     -26.073133d0,     -25.635116d0,     -25.241009d0, &
            -24.881752d0,     -24.546616d0,     -24.224041d0,     -23.904261d0,     -23.582406d0,     -23.259230d0, &
            -22.938739d0,     -22.625220d0,     -22.321841d0,     -22.030791d0,     -21.753417d0,     -21.490284d0, &
            -21.241219d0,     -21.005478d0,     -20.782128d0,     -20.570531d0,     -20.370390d0,     -20.181746d0, &
            -20.004883d0,     -19.840152d0,     -19.687760d0,     -19.547547d0,     -19.419141d0,     -19.302036d0, &
            -19.195555d0,     -19.098890d0,     -19.011181d0,     -18.931385d0,     -18.858601d0,     -18.792102d0, &
            -18.731350d0,     -18.675942d0,     -18.625562d0,     -18.579939d0,     -28.916583d0,     -28.225880d0, &
            -27.610184d0,     -27.061037d0,     -26.570286d0,     -26.129204d0,     -25.727196d0,     -25.350672d0, &
            -24.983971d0,     -24.614193d0,     -24.237136d0,     -23.857697d0,     -23.484116d0,     -23.122727d0, &
            -22.776701d0,     -22.447172d0,     -22.134656d0,     -21.839391d0,     -21.561216d0,     -21.299430d0, &
            -21.052891d0,     -20.820374d0,     -20.601046d0,     -20.394468d0,     -20.200558d0,     -20.019470d0, &
            -19.851415d0,     -19.696456d0,     -19.554286d0,     -19.424409d0,     -19.306209d0,     -19.198915d0, &
            -19.101644d0,     -19.013481d0,     -18.933352d0,     -18.860323d0,     -18.793641d0,     -18.732745d0, &
            -18.677221d0,     -18.626745d0,     -18.581039d0,     -29.416024d0,     -28.725329d0,     -28.109547d0, &
            -27.559964d0,     -27.067681d0,     -26.622161d0,     -26.209339d0,     -25.811043d0,     -25.409312d0, &
            -24.995390d0,     -24.573557d0,     -24.154315d0,     -23.746117d0,     -23.353116d0,     -22.976893d0, &
            -22.618433d0,     -22.278980d0,     -21.959551d0,     -21.660335d0,     -21.380552d0,     -21.118748d0, &
            -20.873355d0,     -20.643256d0,     -20.427777d0,     -20.226627d0,     -20.039758d0,     -19.867170d0, &
            -19.708714d0,     -19.563885d0,     -19.432011d0,     -19.312326d0,     -19.203930d0,     -19.105837d0, &
            -19.017057d0,     -18.936467d0,     -18.863098d0,     -18.796156d0,     -18.735055d0,     -18.679362d0, &
            -18.628743d0,     -18.582913d0,     -29.915844d0,     -29.225139d0,     -28.609268d0,     -28.059310d0, &
            -27.565767d0,     -27.116714d0,     -26.695628d0,     -26.281636d0,     -25.856299d0,     -25.413524d0, &
            -24.959925d0,     -24.505516d0,     -24.057647d0,     -23.621529d0,     -23.202023d0,     -22.803712d0, &
            -22.430121d0,     -22.082787d0,     -21.761219d0,     -21.463492d0,     -21.187052d0,     -20.929501d0, &
            -20.689188d0,     -20.465116d0,     -20.256803d0,     -20.064052d0,     -19.886722d0,     -19.724508d0, &
            -19.576754d0,     -19.442636d0,     -19.321251d0,     -19.211572d0,     -19.112502d0,     -19.022973d0, &
            -18.941812d0,     -18.868007d0,     -18.800725d0,     -18.739349d0,     -18.683424d0,     -18.632601d0, &
            -18.586586d0,     -30.415786d0,     -29.725074d0,     -29.109156d0,     -28.559009d0,     -28.064837d0, &
            -27.614024d0,     -27.188726d0,     -26.766066d0,     -26.325011d0,     -25.856278d0,     -25.363806d0, &
            -24.858615d0,     -24.354363d0,     -23.864694d0,     -23.400369d0,     -22.967699d0,     -22.568971d0, &
            -22.203285d0,     -21.867687d0,     -21.558295d0,     -21.271265d0,     -21.003559d0,     -20.753382d0, &
            -20.519917d0,     -20.302973d0,     -20.102601d0,     -19.918780d0,     -19.751209d0,     -19.599151d0, &
            -19.461651d0,     -19.337655d0,     -19.225973d0,     -19.125358d0,     -19.034627d0,     -18.952540d0, &
            -18.878028d0,     -18.810196d0,     -18.748376d0,     -18.692077d0,     -18.640922d0,     -18.594600d0, &
            -30.915768d0,     -30.225052d0,     -29.609117d0,     -29.058897d0,     -28.564479d0,     -28.112913d0, &
            -27.685379d0,     -27.256164d0,     -26.797951d0,     -26.294115d0,     -25.749041d0,     -25.185356d0, &
            -24.629958d0,     -24.102577d0,     -23.613271d0,     -23.164799d0,     -22.755544d0,     -22.381255d0, &
            -22.036417d0,     -21.715522d0,     -21.414204d0,     -21.129970d0,     -20.862237d0,     -20.611580d0, &
            -20.378929d0,     -20.164973d0,     -19.969886d0,     -19.793255d0,     -19.634088d0,     -19.491117d0, &
            -19.362970d0,     -19.248157d0,     -19.145167d0,     -19.052625d0,     -18.969183d0,     -18.893680d0, &
            -18.825119d0,     -18.762746d0,     -18.706007d0,     -18.654482d0,     -18.607825d0,     -31.415762d0, &
            -30.725045d0,     -30.109104d0,     -29.558858d0,     -29.064345d0,     -28.612405d0,     -28.183282d0, &
            -27.747831d0,     -27.271249d0,     -26.731793d0,     -26.142252d0,     -25.540024d0,     -24.958612d0, &
            -24.416126d0,     -23.918362d0,     -23.463808d0,     -23.046882d0,     -22.659881d0,     -22.295021d0, &
            -21.946665d0,     -21.612757d0,     -21.294538d0,     -20.994894d0,     -20.716584d0,     -20.461310d0, &
            -20.229560d0,     -20.020830d0,     -19.833910d0,     -19.667070d0,     -19.518424d0,     -19.386110d0, &
            -19.268242d0,     -19.162997d0,     -19.068785d0,     -18.984143d0,     -18.907818d0,     -18.838703d0, &
            -18.775959d0,     -18.718970d0,     -18.667268d0,     -18.620476d0,     -31.915760d0,     -31.225043d0, &
            -30.609100d0,     -30.058845d0,     -29.564294d0,     -29.112163d0,     -28.682048d0,     -28.242370d0, &
            -27.753537d0,     -27.192650d0,     -26.580455d0,     -25.961867d0,     -25.370859d0,     -24.822321d0, &
            -24.317827d0,     -23.850986d0,     -23.411105d0,     -22.987344d0,     -22.573637d0,     -22.171089d0, &
            -21.785667d0,     -21.423918d0,     -21.090318d0,     -20.786786d0,     -20.513270d0,     -20.268485d0, &
            -20.050473d0,     -19.856938d0,     -19.685372d0,     -19.533344d0,     -19.398611d0,     -19.279006d0, &
            -19.172504d0,     -19.077382d0,     -18.992105d0,     -18.915366d0,     -18.845994d0,     -18.783102d0, &
            -18.726038d0,     -18.674308d0,     -18.627518d0,     -32.415760d0,     -31.725043d0,     -31.109098d0, &
            -30.558841d0,     -30.064275d0,     -29.612064d0,     -29.181508d0,     -28.739915d0,     -28.245670d0, &
            -27.675888d0,     -27.054974d0,     -26.430085d0,     -25.833701d0,     -25.276540d0,     -24.753551d0, &
            -24.250546d0,     -23.753250d0,     -23.257702d0,     -22.772368d0,     -22.310129d0,     -21.880670d0, &
            -21.488611d0,     -21.134790d0,     -20.817866d0,     -20.535433d0,     -20.284644d0,     -20.062537d0, &
            -19.866173d0,     -19.692632d0,     -19.539215d0,     -19.403501d0,     -19.283199d0,     -19.176197d0, &
            -19.080717d0,     -18.995193d0,     -18.918298d0,     -18.848834d0,     -18.785894d0,     -18.728813d0, &
            -18.677086d0,     -18.630311d0,     -32.915759d0,     -32.225042d0,     -31.609098d0,     -31.058840d0, &
            -30.564269d0,     -30.112030d0,     -29.681314d0,     -29.239028d0,     -28.742841d0,     -28.169882d0, &
            -27.545520d0,     -26.916320d0,     -26.310586d0,     -25.730071d0,     -25.157912d0,     -24.576655d0, &
            -23.988272d0,     -23.412632d0,     -22.869761d0,     -22.370709d0,     -21.918827d0,     -21.513222d0, &
            -21.151111d0,     -20.829002d0,     -20.543244d0,     -20.290273d0,     -20.066704d0,     -19.869342d0, &
            -19.695110d0,     -19.541210d0,     -19.405156d0,     -19.284614d0,     -19.177440d0,     -19.081837d0, &
            -18.996230d0,     -18.919282d0,     -18.849788d0,     -18.786833d0,     -18.729747d0,     -18.678024d0, &
            -18.631257d0,     -33.415759d0,     -32.725042d0,     -32.109098d0,     -31.558839d0,     -31.064267d0, &
            -30.612019d0,     -30.181249d0,     -29.738733d0,     -29.241882d0,     -28.667621d0,     -28.040491d0, &
            -27.402836d0,     -26.771610d0,     -26.131323d0,     -25.461733d0,     -24.774502d0,     -24.104342d0, &
            -23.478062d0,     -22.906887d0,     -22.392409d0,     -21.932000d0,     -21.521532d0,     -21.156546d0, &
            -20.832674d0,     -20.545801d0,     -20.292103d0,     -20.068048d0,     -19.870355d0,     -19.695893d0, &
            -19.541832d0,     -19.405664d0,     -19.285040d0,     -19.177809d0,     -19.082164d0,     -18.996529d0, &
            -18.919563d0,     -18.850060d0,     -18.787101d0,     -18.730015d0,     -18.678295d0,     -18.631533d0, &
            -33.915759d0,     -33.225042d0,     -32.609098d0,     -32.058839d0,     -31.564267d0,     -31.112015d0, &
            -30.681229d0,     -30.238635d0,     -29.741492d0,     -29.165977d0,     -28.532785d0,     -27.871376d0, &
            -27.174611d0,     -26.420256d0,     -25.631752d0,     -24.863104d0,     -24.149139d0,     -23.501275d0, &
            -22.919490d0,     -22.399614d0,     -21.936324d0,     -21.524244d0,     -21.158309d0,     -20.833855d0, &
            -20.546610d0,     -20.292665d0,     -20.068440d0,     -19.870628d0,     -19.696081d0,     -19.541957d0, &
            -19.405743d0,     -19.285085d0,     -19.177828d0,     -19.082165d0,     -18.996517d0,     -18.919545d0, &
            -18.850040d0,     -18.787083d0,     -18.730001d0,     -18.678286d0,     -18.631532d0,     -34.415759d0, &
            -33.725042d0,     -33.109098d0,     -32.558839d0,     -32.064266d0,     -31.612014d0,     -31.181222d0, &
            -30.738591d0,     -30.241107d0,     -29.662556d0,     -29.011695d0,     -28.286999d0,     -27.463163d0, &
            -26.575241d0,     -25.703544d0,     -24.895636d0,     -24.164546d0,     -23.509075d0,     -22.923725d0, &
            -22.402068d0,     -21.937828d0,     -21.525203d0,     -21.158935d0,     -20.834261d0,     -20.546862d0, &
            -20.292801d0,     -20.068487d0,     -19.870602d0,     -19.695996d0,     -19.541823d0,     -19.405568d0, &
            -19.284879d0,     -19.177597d0,     -19.081916d0,     -18.996257d0,     -18.919282d0,     -18.849780d0, &
            -18.786831d0,     -18.729762d0,     -18.678063d0,     -18.631328d0,     -34.915759d0,     -34.225042d0, &
            -33.609098d0,     -33.058839d0,     -32.564266d0,     -32.112013d0,     -31.681219d0,     -31.238537d0, &
            -30.740159d0,     -30.152467d0,     -29.451780d0,     -28.595020d0,     -27.617601d0,     -26.639070d0, &
            -25.729494d0,     -24.906993d0,     -24.170076d0,     -23.512103d0,     -22.925581d0,     -22.403317d0, &
            -21.938722d0,     -21.525861d0,     -21.159415d0,     -20.834593d0,     -20.547062d0,     -20.292884d0, &
            -20.068462d0,     -19.870481d0,     -19.695789d0,     -19.541542d0,     -19.405227d0,     -19.284489d0, &
            -19.177172d0,     -19.081469d0,     -18.995801d0,     -18.918830d0,     -18.849344d0,     -18.786420d0, &
            -18.729384d0,     -18.677724d0,     -18.631032d0,     -35.415759d0,     -34.725042d0,     -34.109098d0, &
            -33.558839d0,     -33.064266d0,     -32.612013d0,     -32.181214d0,     -31.738394d0,     -31.237273d0, &
            -30.622298d0,     -29.804007d0,     -28.767770d0,     -27.681932d0,     -26.662827d0,     -25.739564d0, &
            -24.912225d0,     -24.173434d0,     -23.514644d0,     -22.927703d0,     -22.405170d0,     -21.940361d0, &
            -21.527302d0,     -21.160661d0,     -20.835645d0,     -20.547925d0,     -20.293564d0,     -20.068973d0, &
            -19.870840d0,     -19.696017d0,     -19.541664d0,     -19.405267d0,     -19.284474d0,     -19.177126d0, &
            -19.081417d0,     -18.995767d0,     -18.918834d0,     -18.849406d0,     -18.786555d0,     -18.729606d0, &
            -18.678043d0,     -18.631457d0,     -35.915759d0,     -35.225042d0,     -34.609098d0,     -34.058839d0, &
            -33.564266d0,     -33.112013d0,     -32.681203d0,     -32.237957d0,     -31.728402d0,     -31.039688d0, &
            -30.025153d0,     -28.844261d0,     -27.708473d0,     -26.674909d0,     -25.747258d0,     -24.918459d0, &
            -24.179109d0,     -23.520035d0,     -22.932872d0,     -22.410107d0,     -21.945032d0,     -21.531671d0, &
            -21.164701d0,     -20.839342d0,     -20.551276d0,     -20.296584d0,     -20.071689d0,     -19.873293d0, &
            -19.698258d0,     -19.543749d0,     -19.407257d0,     -19.286429d0,     -19.179109d0,     -19.083486d0, &
            -18.997975d0,     -18.921230d0,     -18.852030d0,     -18.789446d0,     -18.732794d0,     -18.681553d0, &
            -18.635306d0 &
        ], [NCOOLT,NCOOLNE] )

   ! CII
   real*8, parameter :: cool_logLrem_CII(NCOOLT,NCOOLNE) = reshape( [ &
            -44.854834d0,     -41.950813d0,     -39.369656d0,     -37.076996d0,     -35.041835d0,     -33.235718d0, &
            -31.631664d0,     -30.203820d0,     -28.931030d0,     -27.795547d0,     -26.782122d0,     -25.878016d0, &
            -25.072972d0,     -24.358614d0,     -23.725790d0,     -23.165527d0,     -22.669540d0,     -22.230029d0, &
            -21.839512d0,     -21.490872d0,     -21.178099d0,     -20.895894d0,     -20.639455d0,     -20.404551d0, &
            -20.187593d0,     -19.985686d0,     -19.796644d0,     -19.618901d0,     -19.451274d0,     -19.293051d0, &
            -19.143896d0,     -19.003678d0,     -18.872347d0,     -18.749841d0,     -18.636039d0,     -18.530736d0, &
            -18.433667d0,     -18.344530d0,     -18.262999d0,     -18.188716d0,     -18.121215d0,     -44.850042d0, &
            -41.946270d0,     -39.365490d0,     -37.073336d0,     -35.038796d0,     -33.233380d0,     -31.630023d0, &
            -30.202739d0,     -28.930326d0,     -27.795038d0,     -26.781647d0,     -25.877463d0,     -25.072307d0, &
            -24.357885d0,     -23.725056d0,     -23.164845d0,     -22.668954d0,     -22.229570d0,     -21.839183d0, &
            -21.490650d0,     -21.177954d0,     -20.895796d0,     -20.639379d0,     -20.404474d0,     -20.187497d0, &
            -19.985558d0,     -19.796476d0,     -19.618686d0,     -19.451008d0,     -19.292736d0,     -19.143534d0, &
            -19.003275d0,     -18.871910d0,     -18.749377d0,     -18.635555d0,     -18.530239d0,     -18.433165d0, &
            -18.344028d0,     -18.262504d0,     -18.188232d0,     -18.120747d0,     -44.842838d0,     -41.939288d0, &
            -39.358946d0,     -37.067457d0,     -35.033806d0,     -33.229458d0,     -31.627212d0,     -30.200850d0, &
            -28.929077d0,     -27.794120d0,     -26.780778d0,     -25.876440d0,     -25.071064d0,     -24.356504d0, &
            -23.723652d0,     -23.163524d0,     -22.667807d0,     -22.228658d0,     -21.838522d0,     -21.490198d0, &
            -21.177655d0,     -20.895593d0,     -20.639217d0,     -20.404306d0,     -20.187286d0,     -19.985274d0, &
            -19.796096d0,     -19.618192d0,     -19.450390d0,     -19.291992d0,     -19.142672d0,     -19.002308d0, &
            -18.870853d0,     -18.748246d0,     -18.634368d0,     -18.529017d0,     -18.431927d0,     -18.342792d0, &
            -18.261284d0,     -18.187041d0,     -18.119598d0,     -44.836775d0,     -41.933267d0,     -39.353159d0, &
            -37.062123d0,     -35.029158d0,     -33.225707d0,     -31.624454d0,     -30.198952d0,     -28.927792d0, &
            -27.793157d0,     -26.779850d0,     -25.875329d0,     -25.069693d0,     -24.354959d0,     -23.722056d0, &
            -23.161997d0,     -22.666458d0,     -22.227567d0,     -21.837716d0,     -21.489636d0,     -21.177275d0, &
            -20.895328d0,     -20.639001d0,     -20.404079d0,     -20.186991d0,     -19.984867d0,     -19.795537d0, &
            -19.617446d0,     -19.449434d0,     -19.290819d0,     -19.141286d0,     -19.000725d0,     -18.869095d0, &
            -18.746342d0,     -18.632353d0,     -18.526929d0,     -18.429802d0,     -18.340665d0,     -18.259186d0, &
            -18.184998d0,     -18.117633d0,     -44.833706d0,     -41.930166d0,     -39.350125d0,     -37.059274d0, &
            -35.026627d0,     -33.223625d0,     -31.622892d0,     -30.197857d0,     -28.927038d0,     -27.792582d0, &
            -26.779287d0,     -25.874647d0,     -25.068840d0,     -24.353984d0,     -23.721035d0,     -23.161007d0, &
            -22.665571d0,     -22.226837d0,     -21.837167d0,     -21.489247d0,     -21.177007d0,     -20.895137d0, &
            -20.638842d0,     -20.403906d0,     -20.186762d0,     -19.984541d0,     -19.795076d0,     -19.616815d0, &
            -19.448605d0,     -19.289775d0,     -19.140024d0,     -18.999254d0,     -18.867434d0,     -18.744516d0, &
            -18.630397d0,     -18.524885d0,     -18.427712d0,     -18.338569d0,     -18.257119d0,     -18.182991d0, &
            -18.115711d0,     -44.832542d0,     -41.928981d0,     -39.348956d0,     -37.058166d0,     -35.025634d0, &
            -33.222799d0,     -31.622267d0,     -30.197414d0,     -28.926730d0,     -27.792346d0,     -26.779054d0, &
            -25.874362d0,     -25.068481d0,     -24.353571d0,     -23.720599d0,     -23.160581d0,     -22.665185d0, &
            -22.226518d0,     -21.836925d0,     -21.489073d0,     -21.176885d0,     -20.895050d0,     -20.638768d0, &
            -20.403824d0,     -20.186651d0,     -19.984381d0,     -19.794847d0,     -19.616495d0,     -19.448176d0, &
            -19.289225d0,     -19.139349d0,     -18.998455d0,     -18.866519d0,     -18.743499d0,     -18.629299d0, &
            -18.523729d0,     -18.426526d0,     -18.337377d0,     -18.255945d0,     -18.181852d0,     -18.114625d0, &
            -44.832151d0,     -41.928582d0,     -39.348561d0,     -37.057790d0,     -35.025296d0,     -33.222517d0, &
            -31.622053d0,     -30.197262d0,     -28.926624d0,     -27.792264d0,     -26.778973d0,     -25.874263d0, &
            -25.068355d0,     -24.353426d0,     -23.720446d0,     -23.160431d0,     -22.665049d0,     -22.226404d0, &
            -21.836838d0,     -21.489010d0,     -21.176842d0,     -20.895018d0,     -20.638741d0,     -20.403794d0, &
            -20.186611d0,     -19.984322d0,     -19.794760d0,     -19.616373d0,     -19.448012d0,     -19.289013d0, &
            -19.139086d0,     -18.998142d0,     -18.866158d0,     -18.743095d0,     -18.628860d0,     -18.523266d0, &
            -18.426050d0,     -18.336899d0,     -18.255473d0,     -18.181396d0,     -18.114190d0,     -44.832027d0, &
            -41.928455d0,     -39.348435d0,     -37.057671d0,     -35.025188d0,     -33.222428d0,     -31.621985d0, &
            -30.197214d0,     -28.926591d0,     -27.792239d0,     -26.778948d0,     -25.874231d0,     -25.068315d0, &
            -24.353380d0,     -23.720396d0,     -23.160382d0,     -22.665005d0,     -22.226368d0,     -21.836811d0, &
            -21.488990d0,     -21.176828d0,     -20.895008d0,     -20.638733d0,     -20.403785d0,     -20.186598d0, &
            -19.984302d0,     -19.794732d0,     -19.616332d0,     -19.447957d0,     -19.288942d0,     -19.138997d0, &
            -18.998035d0,     -18.866034d0,     -18.742957d0,     -18.628710d0,     -18.523108d0,     -18.425887d0, &
            -18.336735d0,     -18.255312d0,     -18.181240d0,     -18.114041d0,     -44.831996d0,     -41.928423d0, &
            -39.348403d0,     -37.057639d0,     -35.025160d0,     -33.222405d0,     -31.621968d0,     -30.197203d0, &
            -28.926585d0,     -27.792235d0,     -26.778943d0,     -25.874225d0,     -25.068306d0,     -24.353368d0, &
            -23.720384d0,     -23.160370d0,     -22.664994d0,     -22.226359d0,     -21.836804d0,     -21.488986d0, &
            -21.176826d0,     -20.895007d0,     -20.638732d0,     -20.403783d0,     -20.186595d0,     -19.984297d0, &
            -19.794724d0,     -19.616320d0,     -19.447940d0,     -19.288919d0,     -19.138969d0,     -18.998001d0, &
            -18.865995d0,     -18.742913d0,     -18.628662d0,     -18.523057d0,     -18.425834d0,     -18.336682d0, &
            -18.255260d0,     -18.181189d0,     -18.113993d0,     -44.832012d0,     -41.928437d0,     -39.348416d0, &
            -37.057651d0,     -35.025171d0,     -33.222416d0,     -31.621979d0,     -30.197215d0,     -28.926597d0, &
            -27.792247d0,     -26.778954d0,     -25.874235d0,     -25.068314d0,     -24.353376d0,     -23.720390d0, &
            -23.160376d0,     -22.664999d0,     -22.226365d0,     -21.836810d0,     -21.488993d0,     -21.176832d0, &
            -20.895013d0,     -20.638737d0,     -20.403788d0,     -20.186598d0,     -19.984300d0,     -19.794725d0, &
            -19.616319d0,     -19.447937d0,     -19.288914d0,     -19.138962d0,     -18.997992d0,     -18.865984d0, &
            -18.742901d0,     -18.628648d0,     -18.523042d0,     -18.425819d0,     -18.336667d0,     -18.255245d0, &
            -18.181175d0,     -18.113980d0,     -44.832098d0,     -41.928520d0,     -39.348494d0,     -37.057724d0, &
            -35.025238d0,     -33.222476d0,     -31.622034d0,     -30.197265d0,     -28.926644d0,     -27.792291d0, &
            -26.778996d0,     -25.874275d0,     -25.068352d0,     -24.353412d0,     -23.720425d0,     -23.160408d0, &
            -22.665030d0,     -22.226394d0,     -21.836837d0,     -21.489018d0,     -21.176856d0,     -20.895035d0, &
            -20.638757d0,     -20.403806d0,     -20.186614d0,     -19.984313d0,     -19.794736d0,     -19.616328d0, &
            -19.447944d0,     -19.288920d0,     -19.138966d0,     -18.997995d0,     -18.865986d0,     -18.742902d0, &
            -18.628649d0,     -18.523042d0,     -18.425819d0,     -18.336667d0,     -18.255245d0,     -18.181175d0, &
            -18.113980d0,     -44.832382d0,     -41.928793d0,     -39.348752d0,     -37.057964d0,     -35.025459d0, &
            -33.222676d0,     -31.622215d0,     -30.197430d0,     -28.926795d0,     -27.792433d0,     -26.779131d0, &
            -25.874404d0,     -25.068477d0,     -24.353531d0,     -23.720538d0,     -23.160515d0,     -22.665131d0, &
            -22.226488d0,     -21.836926d0,     -21.489100d0,     -21.176932d0,     -20.895105d0,     -20.638821d0, &
            -20.403863d0,     -20.186664d0,     -19.984357d0,     -19.794774d0,     -19.616361d0,     -19.447972d0, &
            -19.288944d0,     -19.138987d0,     -18.998014d0,     -18.866004d0,     -18.742919d0,     -18.628665d0, &
            -18.523058d0,     -18.425835d0,     -18.336682d0,     -18.255260d0,     -18.181190d0,     -18.113995d0, &
            -44.833275d0,     -41.929650d0,     -39.349563d0,     -37.058721d0,     -35.026154d0,     -33.223307d0, &
            -31.622785d0,     -30.197948d0,     -28.927273d0,     -27.792880d0,     -26.779555d0,     -25.874811d0, &
            -25.068869d0,     -24.353907d0,     -23.720896d0,     -23.160854d0,     -22.665450d0,     -22.226787d0, &
            -21.837204d0,     -21.489360d0,     -21.177172d0,     -20.895326d0,     -20.639021d0,     -20.404043d0, &
            -20.186824d0,     -19.984496d0,     -19.794894d0,     -19.616464d0,     -19.448061d0,     -19.289022d0, &
            -19.139057d0,     -18.998078d0,     -18.866064d0,     -18.742976d0,     -18.628721d0,     -18.523112d0, &
            -18.425888d0,     -18.336735d0,     -18.255312d0,     -18.181242d0,     -18.114047d0,     -44.836015d0, &
            -41.932287d0,     -39.352063d0,     -37.061056d0,     -35.028302d0,     -33.225259d0,     -31.624550d0, &
            -30.199557d0,     -28.928758d0,     -27.794271d0,     -26.780876d0,     -25.876078d0,     -25.070090d0, &
            -24.355078d0,     -23.722013d0,     -23.161912d0,     -22.666447d0,     -22.227721d0,     -21.838076d0, &
            -21.490171d0,     -21.177923d0,     -20.896016d0,     -20.639649d0,     -20.404606d0,     -20.187322d0, &
            -19.984931d0,     -19.795271d0,     -19.616789d0,     -19.448342d0,     -19.289268d0,     -19.139277d0, &
            -18.998280d0,     -18.866253d0,     -18.743156d0,     -18.628896d0,     -18.523284d0,     -18.426057d0, &
            -18.336903d0,     -18.255479d0,     -18.181408d0,     -18.114212d0,     -44.843982d0,     -41.939979d0, &
            -39.359386d0,     -37.067923d0,     -35.034649d0,     -33.231054d0,     -31.629815d0,     -30.204373d0, &
            -28.933218d0,     -27.798459d0,     -26.784862d0,     -25.879910d0,     -25.073787d0,     -24.358631d0, &
            -23.725404d0,     -23.165129d0,     -22.669481d0,     -22.230568d0,     -21.840736d0,     -21.492649d0, &
            -21.180220d0,     -20.898128d0,     -20.641571d0,     -20.406331d0,     -20.188850d0,     -19.986266d0, &
            -19.796425d0,     -19.617783d0,     -19.449202d0,     -19.290022d0,     -19.139951d0,     -18.998898d0, &
            -18.866832d0,     -18.743711d0,     -18.629434d0,     -18.523811d0,     -18.426578d0,     -18.337418d0, &
            -18.255990d0,     -18.181916d0,     -18.114718d0,     -44.864739d0,     -41.960135d0,     -39.378703d0, &
            -37.086181d0,     -35.051668d0,     -33.246735d0,     -31.644187d0,     -30.217625d0,     -28.945574d0, &
            -27.810126d0,     -26.796015d0,     -25.890671d0,     -25.084205d0,     -24.368672d0,     -23.735021d0, &
            -23.174284d0,     -22.678143d0,     -22.238721d0,     -21.848377d0,     -21.499784d0,     -21.186846d0, &
            -20.904231d0,     -20.647128d0,     -20.411324d0,     -20.193272d0,     -19.990132d0,     -19.799770d0, &
            -19.620665d0,     -19.451695d0,     -19.292209d0,     -19.141909d0,     -19.000692d0,     -18.868516d0, &
            -18.745322d0,     -18.630998d0,     -18.525344d0,     -18.428090d0,     -18.338915d0,     -18.257476d0, &
            -18.183392d0,     -18.116186d0,     -44.912922d0,     -42.007064d0,     -39.423893d0,     -37.129172d0, &
            -35.092073d0,     -33.284319d0,     -31.678983d0,     -30.250004d0,     -28.975997d0,     -27.839031d0, &
            -26.823776d0,     -25.917556d0,     -25.110316d0,     -24.393933d0,     -23.759317d0,     -23.197516d0, &
            -22.700231d0,     -22.259609d0,     -21.868036d0,     -21.518205d0,     -21.204000d0,     -20.920055d0, &
            -20.661550d0,     -20.424282d0,     -20.204741d0,     -20.000144d0,     -19.808422d0,     -19.628110d0, &
            -19.458131d0,     -19.297848d0,     -19.146953d0,     -19.005314d0,     -18.872849d0,     -18.749464d0, &
            -18.635014d0,     -18.529276d0,     -18.431962d0,     -18.342742d0,     -18.261267d0,     -18.187152d0, &
            -18.119917d0,     -45.019212d0,     -42.110648d0,     -39.523868d0,     -37.224666d0,     -35.182338d0, &
            -33.368876d0,     -31.757854d0,     -30.323878d0,     -29.045759d0,     -27.905538d0,     -26.887779d0, &
            -25.979592d0,     -25.170607d0,     -24.452339d0,     -23.815618d0,     -23.251510d0,     -22.751745d0, &
            -22.308497d0,     -21.914187d0,     -21.561535d0,     -21.244363d0,     -20.957245d0,     -20.695341d0, &
            -20.454500d0,     -20.231324d0,     -20.023186d0,     -19.828178d0,     -19.644975d0,     -19.472594d0, &
            -19.310424d0,     -19.158120d0,     -19.015472d0,     -18.882308d0,     -18.758442d0,     -18.643658d0, &
            -18.537678d0,     -18.440177d0,     -18.350803d0,     -18.269194d0,     -18.194959d0,     -18.127616d0, &
            -45.230992d0,     -42.318046d0,     -39.725486d0,     -37.419090d0,     -35.368259d0,     -33.545333d0, &
            -31.924644d0,     -30.481925d0,     -29.196401d0,     -28.050156d0,     -27.027629d0,     -26.115605d0, &
            -25.303168d0,     -24.581189d0,     -23.940337d0,     -23.371683d0,     -22.866955d0,     -22.418317d0, &
            -22.018183d0,     -21.659251d0,     -21.335194d0,     -21.040465d0,     -20.770265d0,     -20.520662d0, &
            -20.288630d0,     -20.071985d0,     -19.869234d0,     -19.679348d0,     -19.501507d0,     -19.335094d0, &
            -19.179627d0,     -19.034690d0,     -18.899891d0,     -18.774849d0,     -18.659189d0,     -18.552522d0, &
            -18.454452d0,     -18.364585d0,     -18.282535d0,     -18.207901d0,     -18.140197d0,     -45.570824d0, &
            -42.653641d0,     -40.055452d0,     -37.741919d0,     -35.682420d0,     -33.849497d0,     -32.218173d0, &
            -30.765410d0,     -29.471028d0,     -28.317313d0,     -27.288682d0,     -26.371613d0,     -25.554543d0, &
            -24.827506d0,     -24.180887d0,     -23.605624d0,     -23.093234d0,     -22.635582d0,     -22.224745d0, &
            -21.853112d0,     -21.514018d0,     -21.201904d0,     -20.912487d0,     -20.642826d0,     -20.391119d0, &
            -20.156329d0,     -19.937778d0,     -19.734833d0,     -19.546706d0,     -19.372515d0,     -19.211339d0, &
            -19.062277d0,     -18.924494d0,     -18.797244d0,     -18.679883d0,     -18.571840d0,     -18.472606d0, &
            -18.381724d0,     -18.298773d0,     -18.223334d0,     -18.154916d0,     -46.003372d0,     -43.083796d0, &
            -40.482391d0,     -38.164703d0,     -36.100024d0,     -34.260931d0,     -32.622788d0,     -31.163322d0, &
            -29.862796d0,     -28.703716d0,     -27.670562d0,     -26.749639d0,     -25.928906d0,     -25.197667d0, &
            -24.545834d0,     -23.963760d0,     -23.442008d0,     -22.971148d0,     -22.541885d0,     -22.145710d0, &
            -21.776032d0,     -21.428831d0,     -21.102549d0,     -20.797259d0,     -20.513561d0,     -20.251769d0, &
            -20.011556d0,     -19.791936d0,     -19.591414d0,     -19.408251d0,     -19.240686d0,     -19.087089d0, &
            -18.946045d0,     -18.816381d0,     -18.697154d0,     -18.587601d0,     -18.487099d0,     -18.395124d0, &
            -18.311220d0,     -18.234947d0,     -18.165803d0,     -46.479465d0,     -43.558932d0,     -40.956234d0, &
            -38.636855d0,     -36.570038d0,     -34.728352d0,     -33.087282d0,     -31.624871d0,     -30.321583d0, &
            -29.160036d0,     -28.124750d0,     -27.201923d0,     -26.379167d0,     -25.645076d0,     -24.988512d0, &
            -24.397988d0,     -23.861234d0,     -23.365549d0,     -22.899586d0,     -22.456014d0,     -22.032996d0, &
            -21.632849d0,     -21.259183d0,     -20.914687d0,     -20.600331d0,     -20.315499d0,     -20.058458d0, &
            -19.826806d0,     -19.617827d0,     -19.428790d0,     -19.257163d0,     -19.100740d0,     -18.957688d0, &
            -18.826548d0,     -18.706188d0,     -18.595726d0,     -18.494467d0,     -18.401849d0,     -18.317390d0, &
            -18.240639d0,     -18.171084d0,     -46.971601d0,     -44.050740d0,     -41.447596d0,     -39.127634d0, &
            -37.060074d0,     -35.217479d0,     -33.575376d0,     -32.111914d0,     -30.807631d0,     -29.645179d0, &
            -28.609058d0,     -27.685296d0,     -26.860936d0,     -26.123149d0,     -25.457862d0,     -24.848588d0, &
            -24.277101d0,     -23.728006d0,     -23.195064d0,     -22.682310d0,     -22.198327d0,     -21.750194d0, &
            -21.341326d0,     -20.971983d0,     -20.640411d0,     -20.343750d0,     -20.078616d0,     -19.841445d0, &
            -19.628715d0,     -19.437135d0,     -19.263780d0,     -19.106167d0,     -18.962276d0,     -18.830523d0, &
            -18.709694d0,     -18.598857d0,     -18.497289d0,     -18.404408d0,     -18.319724d0,     -18.242781d0, &
            -18.173062d0,     -47.469082d0,     -44.548115d0,     -41.944826d0,     -39.624673d0,     -37.556871d0, &
            -35.713979d0,     -34.071537d0,     -32.607729d0,     -31.303112d0,     -30.140328d0,     -29.103758d0, &
            -28.178936d0,     -27.351333d0,     -26.604104d0,     -25.915376d0,     -25.259020d0,     -24.614587d0, &
            -23.980033d0,     -23.369138d0,     -22.797216d0,     -22.273004d0,     -21.798888d0,     -21.373469d0, &
            -20.993525d0,     -20.655077d0,     -20.353898d0,     -20.085763d0,     -19.846588d0,     -19.632517d0, &
            -19.440035d0,     -19.266070d0,     -19.108039d0,     -18.963854d0,     -18.831886d0,     -18.710891d0, &
            -18.599922d0,     -18.498243d0,     -18.405268d0,     -18.320504d0,     -18.243492d0,     -18.173714d0, &
            -47.968282d0,     -45.047281d0,     -42.443946d0,     -40.123732d0,     -38.055853d0,     -36.212866d0, &
            -34.570315d0,     -33.106395d0,     -31.801655d0,     -30.638661d0,     -29.601402d0,     -28.673881d0, &
            -27.836861d0,     -27.062707d0,     -26.314936d0,     -25.566274d0,     -24.821713d0,     -24.107441d0, &
            -23.444665d0,     -22.842126d0,     -22.300330d0,     -21.816009d0,     -21.384508d0,     -21.000826d0, &
            -20.660015d0,     -20.357308d0,     -20.088170d0,     -19.848328d0,     -19.633811d0,     -19.441030d0, &
            -19.266861d0,     -19.108688d0,     -18.964401d0,     -18.832355d0,     -18.711298d0,     -18.600276d0, &
            -18.498552d0,     -18.405537d0,     -18.320738d0,     -18.243695d0,     -18.173891d0,     -48.468029d0, &
            -45.547017d0,     -42.943667d0,     -40.623434d0,     -38.555530d0,     -36.712514d0,     -35.069927d0, &
            -33.605965d0,     -32.301138d0,     -31.137754d0,     -30.098561d0,     -29.162880d0,     -28.298429d0, &
            -27.456618d0,     -26.599210d0,     -25.738758d0,     -24.915941d0,     -24.157639d0,     -23.472055d0, &
            -22.857744d0,     -22.309677d0,     -21.821862d0,     -21.388324d0,     -21.003405d0,     -20.661814d0, &
            -20.358602d0,     -20.089130d0,     -19.849065d0,     -19.634394d0,     -19.441506d0,     -19.267258d0, &
            -19.109025d0,     -18.964687d0,     -18.832594d0,     -18.711493d0,     -18.600428d0,     -18.498662d0, &
            -18.405606d0,     -18.320768d0,     -18.243689d0,     -18.173854d0,     -48.967949d0,     -46.046933d0, &
            -43.443579d0,     -41.123340d0,     -39.055428d0,     -37.212402d0,     -35.569803d0,     -34.105810d0, &
            -32.800801d0,     -31.636277d0,     -30.591150d0,     -29.631323d0,     -28.696867d0,     -27.731321d0, &
            -26.750624d0,     -25.812643d0,     -24.951635d0,     -24.175806d0,     -23.482034d0,     -22.863694d0, &
            -22.313517d0,     -21.824525d0,     -21.390289d0,     -21.004933d0,     -20.663055d0,     -20.359648d0, &
            -20.090037d0,     -19.849868d0,     -19.635117d0,     -19.442160d0,     -19.267849d0,     -19.109551d0, &
            -18.965145d0,     -18.832980d0,     -18.711799d0,     -18.600648d0,     -18.498792d0,     -18.405647d0, &
            -18.320722d0,     -18.243562d0,     -18.173655d0,     -49.467924d0,     -46.546907d0,     -43.943551d0, &
            -41.623310d0,     -39.555396d0,     -37.712366d0,     -36.069758d0,     -34.605699d0,     -33.300158d0, &
            -32.132135d0,     -31.069272d0,     -30.045871d0,     -28.979483d0,     -27.876654d0,     -26.815485d0, &
            -25.842257d0,     -24.966797d0,     -24.184796d0,     -23.488168d0,     -22.868372d0,     -22.317370d0, &
            -21.827861d0,     -21.393266d0,     -21.007641d0,     -20.665550d0,     -20.361967d0,     -20.092205d0, &
            -19.851902d0,     -19.637026d0,     -19.443945d0,     -19.269506d0,     -19.111070d0,     -18.966514d0, &
            -18.834183d0,     -18.712821d0,     -18.601478d0,     -18.499422d0,     -18.406075d0,     -18.320955d0, &
            -18.243611d0,     -18.173537d0,     -49.967915d0,     -47.046898d0,     -44.443542d0,     -42.123301d0, &
            -40.055385d0,     -38.212354d0,     -36.569729d0,     -35.105482d0,     -33.798355d0,     -32.620043d0, &
            -31.509399d0,     -30.358449d0,     -29.138423d0,     -27.946111d0,     -26.849201d0,     -25.862649d0, &
            -24.981802d0,     -24.197279d0,     -23.499242d0,     -22.878508d0,     -22.326787d0,     -21.836669d0, &
            -21.401530d0,     -21.015405d0,     -20.672852d0,     -20.368839d0,     -20.098680d0,     -19.858004d0, &
            -19.642772d0,     -19.449344d0,     -19.274559d0,     -19.115769d0,     -18.970848d0,     -18.838139d0, &
            -18.716389d0,     -18.604650d0,     -18.502200d0,     -18.408468d0,     -18.322980d0,     -18.245296d0, &
            -18.174912d0 &
        ], [NCOOLT,NCOOLNE] )

   ! NI
   real*8, parameter :: cool_logLrem_NI(NCOOLT,NCOOLNE) = reshape( [ &
            -31.835492d0,     -30.496014d0,     -29.299104d0,     -28.229327d0,     -27.272958d0,     -26.417794d0, &
            -25.652978d0,     -24.968831d0,     -24.356648d0,     -23.808676d0,     -23.318029d0,     -22.878553d0, &
            -22.484751d0,     -22.131772d0,     -21.815339d0,     -21.531572d0,     -21.276757d0,     -21.047723d0, &
            -20.841843d0,     -20.656993d0,     -20.491364d0,     -20.343066d0,     -20.209730d0,     -20.089171d0, &
            -19.979378d0,     -19.878813d0,     -19.785827d0,     -19.697635d0,     -19.610235d0,     -19.521191d0, &
            -19.430054d0,     -19.338170d0,     -19.246617d0,     -19.155557d0,     -19.065135d0,     -18.975990d0, &
            -18.889572d0,     -18.807428d0,     -18.730735d0,     -18.660200d0,     -18.596114d0,     -31.835636d0, &
            -30.496169d0,     -29.299271d0,     -28.229506d0,     -27.273151d0,     -26.418000d0,     -25.653199d0, &
            -24.969067d0,     -24.356899d0,     -23.808942d0,     -23.318310d0,     -22.878848d0,     -22.485059d0, &
            -22.132092d0,     -21.815670d0,     -21.531914d0,     -21.277109d0,     -21.048088d0,     -20.842221d0, &
            -20.657388d0,     -20.491779d0,     -20.343503d0,     -20.210193d0,     -20.089662d0,     -19.979899d0, &
            -19.879357d0,     -19.786386d0,     -19.698195d0,     -19.610782d0,     -19.521714d0,     -19.430547d0, &
            -19.338629d0,     -19.247041d0,     -19.155946d0,     -19.065492d0,     -18.976320d0,     -18.889885d0, &
            -18.807729d0,     -18.731029d0,     -18.660490d0,     -18.596402d0,     -31.836092d0,     -30.496659d0, &
            -29.299797d0,     -28.230072d0,     -27.273757d0,     -26.418650d0,     -25.653895d0,     -24.969811d0, &
            -24.357691d0,     -23.809782d0,     -23.319196d0,     -22.879778d0,     -22.486030d0,     -22.133101d0, &
            -21.816714d0,     -21.532991d0,     -21.278220d0,     -21.049236d0,     -20.843413d0,     -20.658631d0, &
            -20.493083d0,     -20.344879d0,     -20.211650d0,     -20.091208d0,     -19.981533d0,     -19.881067d0, &
            -19.788141d0,     -19.699950d0,     -19.612495d0,     -19.523350d0,     -19.432088d0,     -19.340064d0, &
            -19.248362d0,     -19.157156d0,     -19.066602d0,     -18.977350d0,     -18.890857d0,     -18.808664d0, &
            -18.731942d0,     -18.661392d0,     -18.597296d0,     -31.837528d0,     -30.498203d0,     -29.301455d0, &
            -28.231852d0,     -27.275667d0,     -26.420696d0,     -25.656084d0,     -24.972148d0,     -24.360179d0, &
            -23.812419d0,     -23.321977d0,     -22.882695d0,     -22.489073d0,     -22.136261d0,     -21.819980d0, &
            -21.536359d0,     -21.281690d0,     -21.052819d0,     -20.847128d0,     -20.662504d0,     -20.497142d0, &
            -20.349153d0,     -20.216170d0,     -20.095996d0,     -19.986588d0,     -19.886346d0,     -19.793548d0, &
            -19.705348d0,     -19.617749d0,     -19.528354d0,     -19.436784d0,     -19.344422d0,     -19.252364d0, &
            -19.160811d0,     -19.069950d0,     -18.980452d0,     -18.893783d0,     -18.811478d0,     -18.734692d0, &
            -18.664106d0,     -18.599990d0,     -31.842018d0,     -30.503026d0,     -29.306632d0,     -28.237403d0, &
            -27.281615d0,     -26.427063d0,     -25.662888d0,     -24.979402d0,     -24.367889d0,     -23.820578d0, &
            -23.330567d0,     -22.891688d0,     -22.498440d0,     -22.145967d0,     -21.829994d0,     -21.546658d0, &
            -21.292276d0,     -21.063719d0,     -20.858396d0,     -20.674213d0,     -20.509373d0,     -20.361988d0, &
            -20.229694d0,     -20.110262d0,     -20.001576d0,     -19.901914d0,     -19.809399d0,     -19.721072d0, &
            -19.632940d0,     -19.542700d0,     -19.450124d0,     -19.356688d0,     -19.263534d0,     -19.170943d0, &
            -19.079181d0,     -18.988971d0,     -18.901802d0,     -18.819184d0,     -18.742219d0,     -18.671540d0, &
            -18.607377d0,     -31.855748d0,     -30.517737d0,     -29.322378d0,     -28.254243d0,     -27.299605d0, &
            -26.446254d0,     -25.683323d0,     -25.001107d0,     -24.390863d0,     -23.844789d0,     -23.355947d0, &
            -22.918144d0,     -22.525863d0,     -22.174243d0,     -21.859013d0,     -21.576336d0,     -21.322592d0, &
            -21.094725d0,     -20.890219d0,     -20.707027d0,     -20.543370d0,     -20.397360d0,     -20.266615d0, &
            -20.148799d0,     -20.041586d0,     -19.942921d0,     -19.850540d0,     -19.761213d0,     -19.671002d0, &
            -19.577904d0,     -19.482143d0,     -19.385494d0,     -19.289254d0,     -19.193894d0,     -19.099830d0, &
            -19.007861d0,     -18.919488d0,     -18.836131d0,     -18.758764d0,     -18.687901d0,     -18.623674d0, &
            -31.895234d0,     -30.559790d0,     -29.367113d0,     -28.301769d0,     -27.350023d0,     -26.499647d0, &
            -25.739738d0,     -25.060542d0,     -24.453250d0,     -23.909977d0,     -23.423702d0,     -22.988161d0, &
            -22.597794d0,     -22.247712d0,     -21.933650d0,     -21.651836d0,     -21.398802d0,     -21.171665d0, &
            -20.968073d0,     -20.786088d0,     -20.623961d0,     -20.479767d0,     -20.351008d0,     -20.235023d0, &
            -20.128954d0,     -20.030027d0,     -19.935268d0,     -19.841078d0,     -19.743873d0,     -19.642565d0, &
            -19.538509d0,     -19.434198d0,     -19.331221d0,     -19.230265d0,     -19.131832d0,     -19.036685d0, &
            -18.946209d0,     -18.861600d0,     -18.783589d0,     -18.712477d0,     -18.648243d0,     -31.994847d0, &
            -30.664774d0,     -29.477604d0,     -28.417885d0,     -27.471847d0,     -26.627204d0,     -25.872965d0, &
            -25.199257d0,     -24.597120d0,     -24.058496d0,     -23.576161d0,     -23.143669d0,     -22.755314d0, &
            -22.406120d0,     -22.091827d0,     -21.808803d0,     -21.553896d0,     -21.324586d0,     -21.118859d0, &
            -20.934998d0,     -20.771312d0,     -20.625751d0,     -20.495471d0,     -20.377126d0,     -20.266938d0, &
            -20.161186d0,     -20.056278d0,     -19.948759d0,     -19.836239d0,     -19.719499d0,     -19.601589d0, &
            -19.485766d0,     -19.373631d0,     -19.265688d0,     -19.162157d0,     -19.063491d0,     -18.970769d0, &
            -18.884864d0,     -18.806219d0,     -18.734913d0,     -18.670763d0,     -32.200536d0,     -30.878663d0, &
            -29.699691d0,     -28.648116d0,     -27.710090d0,     -26.873194d0,     -26.126237d0,     -25.459056d0, &
            -24.862329d0,     -24.327545d0,     -23.846986d0,     -23.413760d0,     -23.021860d0,     -22.666255d0, &
            -22.342925d0,     -22.048787d0,     -21.781507d0,     -21.539412d0,     -21.321220d0,     -21.125673d0, &
            -20.951194d0,     -20.795517d0,     -20.655208d0,     -20.526058d0,     -20.403400d0,     -20.283006d0, &
            -20.161428d0,     -20.036084d0,     -19.906109d0,     -19.773945d0,     -19.643618d0,     -19.518405d0, &
            -19.399396d0,     -19.286547d0,     -19.179621d0,     -19.078700d0,     -18.984578d0,     -18.897883d0, &
            -18.818867d0,     -18.747469d0,     -18.683410d0,     -32.534055d0,     -31.220067d0,     -30.048666d0, &
            -29.004281d0,     -28.072911d0,     -27.241856d0,     -26.499434d0,     -25.834729d0,     -25.237439d0, &
            -24.697928d0,     -24.207493d0,     -23.758753d0,     -23.345996d0,     -22.965312d0,     -22.614389d0, &
            -22.292038d0,     -21.997621d0,     -21.730674d0,     -21.490570d0,     -21.276150d0,     -21.085519d0, &
            -20.915804d0,     -20.762734d0,     -20.621263d0,     -20.486110d0,     -20.352904d0,     -20.218526d0, &
            -20.081027d0,     -19.940304d0,     -19.799425d0,     -19.662559d0,     -19.532673d0,     -19.410394d0, &
            -19.295290d0,     -19.186840d0,     -19.084924d0,     -18.990189d0,     -18.903152d0,     -18.823976d0, &
            -18.752541d0,     -18.688526d0,     -32.962658d0,     -31.653146d0,     -30.485830d0,     -29.444976d0, &
            -28.516189d0,     -27.685964d0,     -26.941228d0,     -26.269042d0,     -25.656790d0,     -25.093090d0, &
            -24.569200d0,     -24.079952d0,     -23.623454d0,     -23.199830d0,     -22.809806d0,     -22.453796d0, &
            -22.131555d0,     -21.842233d0,     -21.584446d0,     -21.356128d0,     -21.154479d0,     -20.975789d0, &
            -20.815012d0,     -20.666479d0,     -20.524509d0,     -20.384632d0,     -20.243870d0,     -20.100547d0, &
            -19.954875d0,     -19.810121d0,     -19.670423d0,     -19.538542d0,     -19.414864d0,     -19.298787d0, &
            -19.189668d0,     -19.087307d0,     -18.992292d0,     -18.905091d0,     -18.825833d0,     -18.754371d0, &
            -18.690365d0,     -33.437081d0,     -32.129278d0,     -30.963245d0,     -29.922729d0,     -28.992135d0, &
            -28.155646d0,     -27.396574d0,     -26.697867d0,     -26.044674d0,     -25.427953d0,     -24.845827d0, &
            -24.300997d0,     -23.796846d0,     -23.335292d0,     -22.916445d0,     -22.539088d0,     -22.201225d0, &
            -21.900514d0,     -21.634406d0,     -21.399941d0,     -21.193631d0,     -21.011207d0,     -20.847145d0, &
            -20.695382d0,     -20.549965d0,     -20.406331d0,     -20.261608d0,     -20.114421d0,     -19.965367d0, &
            -19.817932d0,     -19.676250d0,     -19.542922d0,     -19.418163d0,     -19.301275d0,     -19.191558d0, &
            -19.088770d0,     -18.993469d0,     -18.906088d0,     -18.826724d0,     -18.755208d0,     -18.691182d0, &
            -33.928561d0,     -32.621085d0,     -31.454549d0,     -30.411250d0,     -29.472412d0,     -28.616756d0, &
            -27.821338d0,     -27.067118d0,     -26.346354d0,     -25.663014d0,     -25.025345d0,     -24.439178d0, &
            -23.906207d0,     -23.425180d0,     -22.993316d0,     -22.607266d0,     -22.263589d0,     -21.958993d0, &
            -21.690284d0,     -21.454039d0,     -21.246371d0,     -21.062601d0,     -20.896765d0,     -20.742319d0, &
            -20.592940d0,     -20.443966d0,     -20.292855d0,     -20.139016d0,     -19.984021d0,     -19.831848d0, &
            -19.686646d0,     -19.550704d0,     -19.423933d0,     -19.305491d0,     -19.194600d0,     -19.090958d0, &
            -18.995071d0,     -18.907310d0,     -18.827714d0,     -18.756066d0,     -18.691974d0,     -34.425664d0, &
            -33.117724d0,     -31.948976d0,     -30.898581d0,     -29.941523d0,     -29.048500d0,     -28.194255d0, &
            -27.370693d0,     -26.587272d0,     -25.857489d0,     -25.189271d0,     -24.584141d0,     -24.039677d0, &
            -23.551772d0,     -23.115879d0,     -22.727590d0,     -22.382835d0,     -22.077908d0,     -21.809295d0, &
            -21.573242d0,     -21.365385d0,     -21.180295d0,     -21.011033d0,     -20.849964d0,     -20.690119d0, &
            -20.527141d0,     -20.359952d0,     -20.190225d0,     -20.021698d0,     -19.859153d0,     -19.706485d0, &
            -19.565165d0,     -19.434389d0,     -19.312950d0,     -19.199848d0,     -19.094628d0,     -18.997672d0, &
            -18.909222d0,     -18.829203d0,     -18.757311d0,     -18.693091d0,     -34.924526d0,     -33.615683d0, &
            -32.443543d0,     -31.383007d0,     -30.401891d0,     -29.466139d0,     -28.558080d0,     -27.685991d0, &
            -26.868527d0,     -26.118107d0,     -25.438198d0,     -24.826662d0,     -24.278852d0,     -23.789420d0, &
            -23.353127d0,     -22.965168d0,     -22.621223d0,     -22.317357d0,     -22.049703d0,     -21.813869d0, &
            -21.604242d0,     -21.413311d0,     -21.231568d0,     -21.049371d0,     -20.860204d0,     -20.663206d0, &
            -20.462294d0,     -20.263402d0,     -20.072506d0,     -19.894129d0,     -19.730747d0,     -19.582137d0, &
            -19.446243d0,     -19.321168d0,     -19.205496d0,     -19.098497d0,     -19.000362d0,     -18.911162d0, &
            -18.830686d0,     -18.758528d0,     -18.694165d0,     -35.424028d0,     -34.114422d0,     -32.939530d0, &
            -31.871040d0,     -30.872314d0,     -29.909042d0,     -28.971922d0,     -28.077937d0,     -27.247189d0, &
            -26.489544d0,     -25.805917d0,     -25.192594d0,     -24.644111d0,     -24.154686d0,     -23.718844d0, &
            -23.331625d0,     -22.988549d0,     -22.685347d0,     -22.417314d0,     -22.178181d0,     -21.958737d0, &
            -21.746277d0,     -21.527124d0,     -21.293070d0,     -21.046208d0,     -20.796031d0,     -20.552551d0, &
            -20.322799d0,     -20.111108d0,     -19.919311d0,     -19.747437d0,     -19.593375d0,     -19.453854d0, &
            -19.326316d0,     -19.208967d0,     -19.100837d0,     -19.001967d0,     -18.912305d0,     -18.831549d0, &
            -18.759228d0,     -18.694777d0,     -35.923831d0,     -34.613845d0,     -33.437600d0,     -32.365272d0, &
            -31.358403d0,     -30.383312d0,     -29.434805d0,     -28.532928d0,     -27.697684d0,     -26.937717d0, &
            -26.252965d0,     -25.639170d0,     -25.090592d0,     -24.601312d0,     -24.165765d0,     -23.778883d0, &
            -23.435900d0,     -23.131644d0,     -22.858952d0,     -22.606144d0,     -22.355357d0,     -22.086939d0, &
            -21.792057d0,     -21.480009d0,     -21.169000d0,     -20.873477d0,     -20.600513d0,     -20.352355d0, &
            -20.129438d0,     -19.930855d0,     -19.754878d0,     -19.598274d0,     -19.457113d0,     -19.328492d0, &
            -19.210418d0,     -19.101808d0,     -19.002627d0,     -18.912773d0,     -18.831900d0,     -18.759511d0, &
            -18.695023d0,     -36.423762d0,     -35.113633d0,     -33.936882d0,     -32.863130d0,     -31.853296d0, &
            -30.874031d0,     -29.921645d0,     -29.017160d0,     -28.180464d0,     -27.419764d0,     -26.734667d0, &
            -26.120737d0,     -25.572149d0,     -25.082931d0,     -24.647456d0,     -24.260440d0,     -23.916287d0, &
            -23.607092d0,     -23.318589d0,     -23.026215d0,     -22.701821d0,     -22.336836d0,     -21.951827d0, &
            -21.574929d0,     -21.223927d0,     -20.905322d0,     -20.619208d0,     -20.363494d0,     -20.136195d0, &
            -19.935045d0,     -19.757546d0,     -19.600014d0,     -19.458263d0,     -19.329255d0,     -19.210925d0, &
            -19.102146d0,     -19.002857d0,     -18.912935d0,     -18.832021d0,     -18.759609d0,     -18.695107d0, &
            -36.923739d0,     -35.613562d0,     -34.436642d0,     -33.362415d0,     -32.351597d0,     -31.370964d0, &
            -30.417320d0,     -29.511999d0,     -28.674842d0,     -27.913909d0,     -27.228705d0,     -26.614734d0, &
            -26.066142d0,     -25.576924d0,     -25.141321d0,     -24.753411d0,     -24.405035d0,     -24.080236d0, &
            -23.747426d0,     -23.366635d0,     -22.928283d0,     -22.466848d0,     -22.021253d0,     -21.611642d0, &
            -21.243776d0,     -20.916396d0,     -20.625568d0,     -20.367235d0,     -20.138445d0,     -19.936433d0, &
            -19.758426d0,     -19.600586d0,     -19.458640d0,     -19.329505d0,     -19.211091d0,     -19.102256d0, &
            -19.002932d0,     -18.912987d0,     -18.832061d0,     -18.759641d0,     -18.695135d0,     -37.423732d0, &
            -36.113540d0,     -34.936565d0,     -33.862185d0,     -32.851051d0,     -31.869980d0,     -30.915935d0, &
            -30.010349d0,     -29.173045d0,     -28.412039d0,     -27.726801d0,     -27.112817d0,     -26.564216d0, &
            -26.074928d0,     -25.638798d0,     -25.247919d0,     -24.886413d0,     -24.517667d0,     -24.086429d0, &
            -23.576129d0,     -23.036230d0,     -22.518521d0,     -22.046142d0,     -21.624125d0,     -21.250346d0, &
            -20.920010d0,     -20.627629d0,     -20.368441d0,     -20.139169d0,     -19.936878d0,     -19.758708d0, &
            -19.600769d0,     -19.458761d0,     -19.329585d0,     -19.211144d0,     -19.102292d0,     -19.002956d0, &
            -18.913005d0,     -18.832074d0,     -18.759652d0,     -18.695144d0,     -37.923730d0,     -36.613533d0, &
            -35.436541d0,     -34.362112d0,     -33.350878d0,     -32.369667d0,     -31.415496d0,     -30.509825d0, &
            -29.672475d0,     -28.911446d0,     -28.226197d0,     -27.612207d0,     -27.063578d0,     -26.574046d0, &
            -26.136214d0,     -25.736082d0,     -25.336696d0,     -24.867027d0,     -24.291291d0,     -23.671465d0, &
            -23.077048d0,     -22.536325d0,     -22.054362d0,     -21.628168d0,     -21.252455d0,     -20.921165d0, &
            -20.628286d0,     -20.368826d0,     -20.139400d0,     -19.937021d0,     -19.758799d0,     -19.600829d0, &
            -19.458800d0,     -19.329612d0,     -19.211162d0,     -19.102305d0,     -19.002965d0,     -18.913011d0, &
            -18.832079d0,     -18.759656d0,     -18.695148d0,     -38.423729d0,     -37.113530d0,     -35.936533d0, &
            -34.862089d0,     -33.850823d0,     -32.869568d0,     -31.915356d0,     -31.009659d0,     -30.172295d0, &
            -29.411258d0,     -28.726006d0,     -28.112006d0,     -27.563289d0,     -27.072979d0,     -26.629813d0, &
            -26.202145d0,     -25.709979d0,     -25.082670d0,     -24.383308d0,     -23.706622d0,     -23.090822d0, &
            -22.542121d0,     -22.056999d0,     -21.629458d0,     -21.253127d0,     -20.921534d0,     -20.628496d0, &
            -20.368950d0,     -20.139476d0,     -19.937069d0,     -19.758830d0,     -19.600851d0,     -19.458816d0, &
            -19.329624d0,     -19.211172d0,     -19.102312d0,     -19.002972d0,     -18.913017d0,     -18.832085d0, &
            -18.759661d0,     -18.695153d0,     -38.923729d0,     -37.613530d0,     -36.436530d0,     -35.362081d0, &
            -34.350806d0,     -33.369537d0,     -32.415312d0,     -31.509607d0,     -30.672238d0,     -29.911199d0, &
            -29.225944d0,     -28.611920d0,     -28.062925d0,     -27.570161d0,     -27.110676d0,     -26.609916d0, &
            -25.954832d0,     -25.181463d0,     -24.417015d0,     -23.718371d0,     -23.095277d0,     -22.543976d0, &
            -22.057843d0,     -21.629873d0,     -21.253346d0,     -20.921657d0,     -20.628569d0,     -20.368996d0, &
            -20.139507d0,     -19.937092d0,     -19.758850d0,     -19.600867d0,     -19.458832d0,     -19.329639d0, &
            -19.211187d0,     -19.102327d0,     -19.002986d0,     -18.913031d0,     -18.832098d0,     -18.759675d0, &
            -18.695167d0,     -39.423729d0,     -38.113529d0,     -36.936530d0,     -35.862079d0,     -34.850800d0, &
            -33.869527d0,     -32.915298d0,     -32.009590d0,     -31.172220d0,     -30.411180d0,     -29.725920d0, &
            -29.111821d0,     -28.561948d0,     -28.061537d0,     -27.555319d0,     -26.904514d0,     -26.073762d0, &
            -25.218061d0,     -24.428259d0,     -23.722169d0,     -23.096711d0,     -22.544580d0,     -22.058126d0, &
            -21.630021d0,     -21.253433d0,     -20.921715d0,     -20.628613d0,     -20.369033d0,     -20.139541d0, &
            -19.937125d0,     -19.758884d0,     -19.600904d0,     -19.458871d0,     -19.329680d0,     -19.211228d0, &
            -19.102369d0,     -19.003029d0,     -18.913074d0,     -18.832140d0,     -18.759716d0,     -18.695208d0, &
            -39.923729d0,     -38.613529d0,     -37.436529d0,     -36.362078d0,     -35.350798d0,     -34.369524d0, &
            -33.415294d0,     -32.509585d0,     -31.672214d0,     -30.911173d0,     -30.225898d0,     -29.611562d0, &
            -29.058928d0,     -28.535395d0,     -27.915978d0,     -27.063781d0,     -26.119378d0,     -25.230352d0, &
            -24.431922d0,     -23.723423d0,     -23.097213d0,     -22.544821d0,     -22.058266d0,     -21.630121d0, &
            -21.253517d0,     -20.921793d0,     -20.628691d0,     -20.369114d0,     -20.139626d0,     -19.937217d0, &
            -19.758984d0,     -19.601013d0,     -19.458988d0,     -19.329804d0,     -19.211357d0,     -19.102500d0, &
            -19.003160d0,     -18.913205d0,     -18.832270d0,     -18.759845d0,     -18.695335d0,     -40.423729d0, &
            -39.113529d0,     -37.936529d0,     -36.862078d0,     -35.850798d0,     -34.869523d0,     -33.915293d0, &
            -33.009584d0,     -32.172212d0,     -31.411169d0,     -30.725847d0,     -30.110761d0,     -29.549532d0, &
            -28.961812d0,     -28.144709d0,     -27.129451d0,     -26.134994d0,     -25.234451d0,     -24.433229d0, &
            -23.723966d0,     -23.097520d0,     -22.545048d0,     -22.058467d0,     -21.630315d0,     -21.253715d0, &
            -20.922000d0,     -20.628910d0,     -20.369346d0,     -20.139876d0,     -19.937487d0,     -19.759280d0, &
            -19.601337d0,     -19.459338d0,     -19.330174d0,     -19.211741d0,     -19.102890d0,     -19.003551d0, &
            -18.913595d0,     -18.832657d0,     -18.760228d0,     -18.695714d0,     -40.923729d0,     -39.613529d0, &
            -38.436529d0,     -37.362078d0,     -36.350798d0,     -35.369523d0,     -34.415292d0,     -33.509583d0, &
            -32.672211d0,     -31.911162d0,     -31.225691d0,     -30.608245d0,     -30.021128d0,     -29.287974d0, &
            -28.252528d0,     -27.152864d0,     -26.140466d0,     -25.236181d0,     -24.434075d0,     -23.724575d0, &
            -23.098061d0,     -22.545570d0,     -22.058990d0,     -21.630850d0,     -21.254269d0,     -20.922581d0, &
            -20.629523d0,     -20.369996d0,     -20.140568d0,     -19.938234d0,     -19.760094d0,     -19.602225d0, &
            -19.460296d0,     -19.331187d0,     -19.212789d0,     -19.103954d0,     -19.004619d0,     -18.914659d0, &
            -18.833715d0,     -18.761279d0,     -18.696758d0,     -41.423729d0,     -40.113529d0,     -38.936529d0, &
            -37.862078d0,     -36.850798d0,     -35.869522d0,     -34.915292d0,     -34.009583d0,     -33.172210d0, &
            -32.411140d0,     -31.725200d0,     -31.100411d0,     -30.442126d0,     -29.478502d0,     -28.294169d0, &
            -27.161730d0,     -26.143441d0,     -25.237979d0,     -24.435604d0,     -23.726033d0,     -23.099498d0, &
            -22.547001d0,     -22.060424d0,     -21.632297d0,     -21.255742d0,     -20.924093d0,     -20.631084d0, &
            -20.371612d0,     -20.142260d0,     -19.940028d0,     -19.762024d0,     -19.604311d0,     -19.462527d0, &
            -19.333535d0,     -19.215215d0,     -19.106423d0,     -19.007109d0,     -18.917163d0,     -18.836234d0, &
            -18.763818d0,     -18.699324d0,     -41.923729d0,     -40.613529d0,     -39.436529d0,     -38.362078d0, &
            -37.350798d0,     -36.369522d0,     -35.415292d0,     -34.509583d0,     -33.672208d0,     -32.911073d0, &
            -32.223666d0,     -31.576733d0,     -30.760019d0,     -29.564339d0,     -28.311350d0,     -27.167853d0, &
            -26.147741d0,     -25.241938d0,     -24.439478d0,     -23.729858d0,     -23.103268d0,     -22.550705d0, &
            -22.064054d0,     -21.635857d0,     -21.259242d0,     -20.927549d0,     -20.634509d0,     -20.375025d0, &
            -20.145699d0,     -19.943553d0,     -19.765705d0,     -19.608197d0,     -19.466628d0,     -19.337835d0, &
            -19.219691d0,     -19.111060d0,     -19.011916d0,     -18.922167d0,     -18.841468d0,     -18.769315d0, &
            -18.705111d0 &
        ], [NCOOLT,NCOOLNE] )

   ! NII
   real*8, parameter :: cool_logLrem_NII(NCOOLT,NCOOLNE) = reshape( [ &
            -28.233068d0,     -27.217061d0,     -26.314092d0,     -25.511834d0,     -24.799293d0,     -24.166660d0, &
            -23.605181d0,     -23.107038d0,     -22.665263d0,     -22.273691d0,     -21.926853d0,     -21.619898d0, &
            -21.348523d0,     -21.108914d0,     -20.897686d0,     -20.711823d0,     -20.548616d0,     -20.405461d0, &
            -20.279813d0,     -20.169288d0,     -20.071665d0,     -19.984861d0,     -19.906895d0,     -19.835805d0, &
            -19.769624d0,     -19.706221d0,     -19.643146d0,     -19.577794d0,     -19.507821d0,     -19.431564d0, &
            -19.348361d0,     -19.258655d0,     -19.163846d0,     -19.065968d0,     -18.967318d0,     -18.870077d0, &
            -18.776178d0,     -18.687047d0,     -18.603525d0,     -18.525998d0,     -18.454541d0,     -28.231130d0, &
            -27.215667d0,     -26.313164d0,     -25.511307d0,     -24.799113d0,     -24.166783d0,     -23.605574d0, &
            -23.107676d0,     -22.666130d0,     -22.274772d0,     -21.928137d0,     -21.621372d0,     -21.350172d0, &
            -21.110719d0,     -20.899622d0,     -20.713858d0,     -20.550712d0,     -20.407578d0,     -20.281918d0, &
            -20.171360d0,     -20.073692d0,     -19.986844d0,     -19.908841d0,     -19.837730d0,     -19.771536d0, &
            -19.708117d0,     -19.645009d0,     -19.579597d0,     -19.509524d0,     -19.433130d0,     -19.349764d0, &
            -19.259890d0,     -19.164927d0,     -19.066918d0,     -18.968162d0,     -18.870835d0,     -18.776861d0, &
            -18.687657d0,     -18.604061d0,     -18.526461d0,     -18.454936d0,     -28.224777d0,     -27.210644d0, &
            -26.309316d0,     -25.508501d0,     -24.797235d0,     -24.165740d0,     -23.605288d0,     -23.108091d0, &
            -22.667204d0,     -22.276472d0,     -21.930431d0,     -21.624224d0,     -21.353543d0,     -21.114558d0, &
            -20.903865d0,     -20.718424d0,     -20.555500d0,     -20.412483d0,     -20.286855d0,     -20.176267d0, &
            -20.078534d0,     -19.991614d0,     -19.913554d0,     -19.842419d0,     -19.776223d0,     -19.712789d0, &
            -19.649622d0,     -19.584075d0,     -19.513772d0,     -19.437049d0,     -19.353287d0,     -19.263001d0, &
            -19.167656d0,     -19.069324d0,     -18.970306d0,     -18.872769d0,     -18.778612d0,     -18.689225d0, &
            -18.605438d0,     -18.527649d0,     -18.455945d0,     -28.209746d0,     -27.197631d0,     -26.298165d0, &
            -25.499074d0,     -24.789411d0,     -24.159416d0,     -23.600379d0,     -23.104533d0,     -22.664950d0, &
            -22.275483d0,     -21.930667d0,     -21.625640d0,     -21.356079d0,     -21.118138d0,     -20.908386d0, &
            -20.723753d0,     -20.561462d0,     -20.418893d0,     -20.293554d0,     -20.183133d0,     -20.085485d0, &
            -19.998608d0,     -19.920595d0,     -19.849545d0,     -19.783457d0,     -19.720102d0,     -19.656932d0, &
            -19.591254d0,     -19.520653d0,     -19.443458d0,     -19.359096d0,     -19.268168d0,     -19.172222d0, &
            -19.073384d0,     -18.973967d0,     -18.876118d0,     -18.781684d0,     -18.692001d0,     -18.607885d0, &
            -18.529747d0,     -18.457690d0,     -28.191442d0,     -27.180706d0,     -26.282556d0,     -25.484727d0, &
            -24.776282d0,     -24.147471d0,     -23.589593d0,     -23.094890d0,     -22.656447d0,     -22.268120d0, &
            -21.924442d0,     -21.620546d0,     -21.352100d0,     -21.115241d0,     -20.906519d0,     -20.722839d0, &
            -20.561386d0,     -20.419509d0,     -20.294728d0,     -20.184747d0,     -20.087440d0,     -20.000829d0, &
            -19.923038d0,     -19.852210d0,     -19.786347d0,     -19.723189d0,     -19.660165d0,     -19.594561d0, &
            -19.523942d0,     -19.446618d0,     -19.362031d0,     -19.270829d0,     -19.174623d0,     -19.075593d0, &
            -18.976071d0,     -18.878185d0,     -18.783721d0,     -18.693939d0,     -18.609620d0,     -18.531180d0, &
            -18.458754d0,     -28.182020d0,     -27.171589d0,     -26.273730d0,     -25.476182d0,     -24.768012d0, &
            -24.139473d0,     -23.581868d0,     -23.087443d0,     -22.649284d0,     -22.261248d0,     -21.917869d0, &
            -21.614279d0,     -21.346143d0,     -21.109597d0,     -20.901189d0,     -20.717821d0,     -20.556667d0, &
            -20.415059d0,     -20.290508d0,     -20.180711d0,     -20.083535d0,     -19.996991d0,     -19.919208d0, &
            -19.848350d0,     -19.782439d0,     -19.719232d0,     -19.656187d0,     -19.590626d0,     -19.520137d0, &
            -19.443029d0,     -19.358708d0,     -19.267788d0,     -19.171874d0,     -19.073163d0,     -18.974003d0, &
            -18.876512d0,     -18.782429d0,     -18.692946d0,     -18.608782d0,     -18.530338d0,     -18.457768d0, &
            -28.186510d0,     -27.175654d0,     -26.277389d0,     -25.479452d0,     -24.770912d0,     -24.142021d0, &
            -23.584084d0,     -23.089344d0,     -22.650887d0,     -22.262569d0,     -21.918922d0,     -21.615077d0, &
            -21.346698d0,     -21.109921d0,     -20.901294d0,     -20.717722d0,     -20.556379d0,     -20.414590d0, &
            -20.289859d0,     -20.179873d0,     -20.082488d0,     -19.995704d0,     -19.917640d0,     -19.846466d0, &
            -19.780214d0,     -19.716663d0,     -19.653301d0,     -19.587486d0,     -19.516844d0,     -19.439698d0, &
            -19.355434d0,     -19.264629d0,     -19.168871d0,     -19.070357d0,     -18.971442d0,     -18.874241d0, &
            -18.780465d0,     -18.691261d0,     -18.607295d0,     -18.528938d0,     -18.456348d0,     -28.214381d0, &
            -27.201996d0,     -26.302275d0,     -25.502953d0,     -24.793098d0,     -24.162958d0,     -23.603838d0, &
            -23.107979d0,     -22.668462d0,     -22.279141d0,     -21.934546d0,     -21.629804d0,     -21.360576d0, &
            -21.122995d0,     -20.913606d0,     -20.729311d0,     -20.567280d0,     -20.424833d0,     -20.299467d0, &
            -20.188862d0,     -20.090863d0,     -20.003460d0,     -19.924763d0,     -19.852935d0,     -19.786004d0, &
            -19.721748d0,     -19.657659d0,     -19.591111d0,     -19.519752d0,     -19.441934d0,     -19.357061d0, &
            -19.265722d0,     -19.169516d0,     -19.070650d0,     -18.971483d0,     -18.874123d0,     -18.780261d0, &
            -18.691012d0,     -18.607008d0,     -18.528595d0,     -18.455923d0,     -28.296333d0,     -27.280129d0, &
            -26.376736d0,     -25.573888d0,     -24.860651d0,     -24.227274d0,     -23.665058d0,     -23.166244d0, &
            -22.723913d0,     -22.331915d0,     -21.984775d0,     -21.677614d0,     -21.406090d0,     -21.166327d0, &
            -20.954862d0,     -20.768584d0,     -20.604650d0,     -20.460368d0,     -20.333226d0,     -20.220887d0, &
            -20.121181d0,     -20.032084d0,     -19.951688d0,     -19.878133d0,     -19.809419d0,     -19.743292d0, &
            -19.677216d0,     -19.608558d0,     -19.534996d0,     -19.454952d0,     -19.367923d0,     -19.274591d0, &
            -19.176627d0,     -19.076276d0,     -18.975895d0,     -18.877570d0,     -18.782953d0,     -18.693116d0, &
            -18.608650d0,     -18.529866d0,     -18.456891d0,     -28.485771d0,     -27.463027d0,     -26.553197d0, &
            -25.744025d0,     -25.024587d0,     -24.385142d0,     -23.817004d0,     -23.312424d0,     -22.864493d0, &
            -22.467065d0,     -22.114668d0,     -21.802424d0,     -21.525985d0,     -21.281468d0,     -21.065393d0, &
            -20.874621d0,     -20.706277d0,     -20.557646d0,     -20.426187d0,     -20.309535d0,     -20.205490d0, &
            -20.112003d0,     -20.027135d0,     -19.948974d0,     -19.875449d0,     -19.804214d0,     -19.732657d0, &
            -19.658125d0,     -19.578393d0,     -19.492096d0,     -19.399016d0,     -19.300115d0,     -19.197260d0, &
            -19.092779d0,     -18.989012d0,     -18.887967d0,     -18.791193d0,     -18.699657d0,     -18.613856d0, &
            -18.534022d0,     -18.460222d0,     -28.811803d0,     -27.782655d0,     -26.866332d0,     -26.050588d0, &
            -25.324513d0,     -24.678381d0,     -24.103521d0,     -23.592203d0,     -23.137537d0,     -22.733388d0, &
            -22.374288d0,     -22.055360d0,     -21.772243d0,     -21.521027d0,     -21.298181d0,     -21.100491d0, &
            -20.924989d0,     -20.768878d0,     -20.629540d0,     -20.504541d0,     -20.391643d0,     -20.288787d0, &
            -20.194019d0,     -20.105372d0,     -20.020666d0,     -19.937409d0,     -19.852883d0,     -19.764515d0, &
            -19.670444d0,     -19.569910d0,     -19.463389d0,     -19.352428d0,     -19.239229d0,     -19.126182d0, &
            -19.015495d0,     -18.908949d0,     -18.807838d0,     -18.712898d0,     -18.624431d0,     -18.542513d0, &
            -18.467082d0,     -29.238263d0,     -28.205494d0,     -27.285411d0,     -26.465772d0,     -25.735665d0, &
            -25.085368d0,     -24.506214d0,     -23.990479d0,     -23.531272d0,     -23.122449d0,     -22.758517d0, &
            -22.434544d0,     -22.146080d0,     -21.889072d0,     -21.659784d0,     -21.454732d0,     -21.270635d0, &
            -21.104399d0,     -20.953171d0,     -20.814399d0,     -20.685886d0,     -20.565761d0,     -20.452309d0, &
            -20.343733d0,     -20.237911d0,     -20.132356d0,     -20.024487d0,     -19.912226d0,     -19.794615d0, &
            -19.671982d0,     -19.545737d0,     -19.417970d0,     -19.290978d0,     -19.166908d0,     -19.047556d0, &
            -18.934254d0,     -18.827891d0,     -18.728864d0,     -18.637215d0,     -18.552818d0,     -18.475452d0, &
            -29.712155d0,     -28.677937d0,     -27.756332d0,     -26.935096d0,     -26.203315d0,     -25.551264d0, &
            -24.970274d0,     -24.452611d0,     -23.991366d0,     -23.580340d0,     -23.213936d0,     -22.887034d0, &
            -22.594874d0,     -22.332937d0,     -22.096857d0,     -21.882385d0,     -21.685450d0,     -21.502357d0, &
            -21.330059d0,     -21.166363d0,     -21.009974d0,     -20.860204d0,     -20.716448d0,     -20.577678d0, &
            -20.442178d0,     -20.307690d0,     -20.171980d0,     -20.033593d0,     -19.892386d0,     -19.749351d0, &
            -19.606191d0,     -19.464885d0,     -19.327331d0,     -19.195144d0,     -19.069594d0,     -18.951564d0, &
            -18.841580d0,     -18.739765d0,     -18.645962d0,     -18.559893d0,     -18.481227d0,     -30.203562d0, &
            -29.168846d0,     -28.246716d0,     -27.424926d0,     -26.692561d0,     -26.039894d0,     -25.458247d0, &
            -24.939865d0,     -24.477778d0,     -24.065652d0,     -23.697610d0,     -23.368025d0,     -23.071321d0, &
            -22.801798d0,     -22.553611d0,     -22.321004d0,     -22.098914d0,     -21.883779d0,     -21.674112d0, &
            -21.470330d0,     -21.273957d0,     -21.086496d0,     -20.908497d0,     -20.739166d0,     -20.576490d0, &
            -20.417824d0,     -20.260724d0,     -20.103777d0,     -19.946955d0,     -19.791235d0,     -19.638099d0, &
            -19.489143d0,     -19.345822d0,     -19.209324d0,     -19.080558d0,     -18.960119d0,     -18.848320d0, &
            -18.745127d0,     -18.650270d0,     -18.563392d0,     -18.484101d0,     -30.700809d0,     -29.665931d0, &
            -28.743631d0,     -27.921660d0,     -27.189104d0,     -26.536231d0,     -25.954351d0,     -25.435667d0, &
            -24.973097d0,     -24.560039d0,     -24.190063d0,     -23.856539d0,     -23.552277d0,     -23.269379d0, &
            -22.999651d0,     -22.735882d0,     -22.473686d0,     -22.212724d0,     -21.956097d0,     -21.708295d0, &
            -21.473198d0,     -21.252884d0,     -21.047310d0,     -20.854614d0,     -20.671754d0,     -20.495374d0, &
            -20.322723d0,     -20.152372d0,     -19.984362d0,     -19.819631d0,     -19.659465d0,     -19.505167d0, &
            -19.357858d0,     -19.218418d0,     -19.087492d0,     -18.965468d0,     -18.852502d0,     -18.748444d0, &
            -18.652942d0,     -18.565581d0,     -18.485926d0,     -31.199934d0,     -30.165006d0,     -29.242650d0, &
            -28.420622d0,     -27.688004d0,     -27.035062d0,     -26.453091d0,     -25.934249d0,     -25.471325d0, &
            -25.057399d0,     -24.685352d0,     -24.347224d0,     -24.033579d0,     -23.733424d0,     -23.435648d0, &
            -23.132468d0,     -22.822860d0,     -22.512477d0,     -22.209834d0,     -21.922273d0,     -21.654104d0, &
            -21.406468d0,     -21.177918d0,     -20.965265d0,     -20.764499d0,     -20.571780d0,     -20.384325d0, &
            -20.200936d0,     -20.021894d0,     -19.848195d0,     -19.680983d0,     -19.521298d0,     -19.369952d0, &
            -19.227528d0,     -19.094416d0,     -18.970797d0,     -18.856668d0,     -18.751759d0,     -18.655633d0, &
            -18.567812d0,     -18.487815d0,     -31.699657d0,     -30.664712d0,     -29.742340d0,     -28.920293d0, &
            -28.187655d0,     -27.534690d0,     -26.952684d0,     -26.433765d0,     -25.970634d0,     -25.556123d0, &
            -25.182472d0,     -24.840172d0,     -24.516540d0,     -24.195363d0,     -23.860802d0,     -23.506096d0, &
            -23.138531d0,     -22.773404d0,     -22.424555d0,     -22.100147d0,     -21.803040d0,     -21.532346d0, &
            -21.284760d0,     -21.055668d0,     -20.840147d0,     -20.633966d0,     -20.434360d0,     -20.240345d0, &
            -20.052380d0,     -19.871470d0,     -19.698609d0,     -19.534603d0,     -19.380010d0,     -19.235176d0, &
            -19.100289d0,     -18.975370d0,     -18.860289d0,     -18.754683d0,     -18.658044d0,     -18.569844d0, &
            -18.489567d0,     -32.199570d0,     -31.164620d0,     -30.242242d0,     -29.420189d0,     -28.687545d0, &
            -28.034572d0,     -27.452553d0,     -26.933601d0,     -26.470355d0,     -26.055407d0,     -25.680096d0, &
            -25.332055d0,     -24.991643d0,     -24.632222d0,     -24.234183d0,     -23.803872d0,     -23.367474d0, &
            -22.949202d0,     -22.562202d0,     -22.210757d0,     -21.894065d0,     -21.608553d0,     -21.349084d0, &
            -21.109865d0,     -20.885305d0,     -20.670913d0,     -20.463949d0,     -20.263562d0,     -20.070302d0, &
            -19.885155d0,     -19.709007d0,     -19.542508d0,     -19.386057d0,     -19.239857d0,     -19.103975d0, &
            -18.978334d0,     -18.862731d0,     -18.756748d0,     -18.659836d0,     -18.571436d0,     -18.491012d0, &
            -32.699542d0,     -31.664590d0,     -30.742211d0,     -29.920156d0,     -29.187510d0,     -28.534534d0, &
            -27.952510d0,     -27.433538d0,     -26.970173d0,     -26.554562d0,     -26.176056d0,     -25.815563d0, &
            -25.438652d0,     -25.006529d0,     -24.516067d0,     -24.005924d0,     -23.514744d0,     -23.061839d0, &
            -22.652753d0,     -22.286469d0,     -21.959043d0,     -21.665062d0,     -21.398278d0,     -21.152252d0, &
            -20.921109d0,     -20.700382d0,     -20.487545d0,     -20.281998d0,     -20.084442d0,     -19.895888d0, &
            -19.717143d0,     -19.548732d0,     -19.390913d0,     -19.243757d0,     -19.107222d0,     -18.981145d0, &
            -18.865255d0,     -18.759086d0,     -18.662055d0,     -18.573579d0,     -18.493103d0,     -33.199533d0, &
            -32.164581d0,     -31.242201d0,     -30.420146d0,     -29.687499d0,     -29.034522d0,     -28.452495d0, &
            -27.933499d0,     -27.469951d0,     -27.053159d0,     -26.668614d0,     -26.284766d0,     -25.847126d0, &
            -25.322686d0,     -24.749700d0,     -24.186673d0,     -23.664777d0,     -23.193139d0,     -22.771303d0, &
            -22.395168d0,     -22.059104d0,     -21.756635d0,     -21.480816d0,     -21.224861d0,     -20.982996d0, &
            -20.751298d0,     -20.528028d0,     -20.313278d0,     -20.108124d0,     -19.913637d0,     -19.730476d0, &
            -19.558909d0,     -19.398917d0,     -19.250321d0,     -19.112869d0,     -18.986233d0,     -18.870026d0, &
            -18.763692d0,     -18.666590d0,     -18.578094d0,     -18.497623d0,     -33.699531d0,     -32.664578d0, &
            -31.742198d0,     -30.920142d0,     -30.187496d0,     -29.534518d0,     -28.952489d0,     -28.433469d0, &
            -27.969718d0,     -27.551552d0,     -27.159833d0,     -26.749249d0,     -26.252256d0,     -25.659237d0, &
            -25.038848d0,     -24.448567d0,     -23.909853d0,     -23.425955d0,     -22.993534d0,     -22.606810d0, &
            -22.258804d0,     -21.941819d0,     -21.648053d0,     -21.370691d0,     -21.105091d0,     -20.849379d0, &
            -20.603981d0,     -20.370460d0,     -20.150399d0,     -19.944679d0,     -19.753404d0,     -19.576187d0, &
            -19.412397d0,     -19.261335d0,     -19.122347d0,     -18.994803d0,     -18.878100d0,     -18.771534d0, &
            -18.674357d0,     -18.585878d0,     -18.505471d0,     -34.199530d0,     -33.164577d0,     -32.242197d0, &
            -31.420141d0,     -30.687494d0,     -30.034517d0,     -29.452486d0,     -28.933452d0,     -28.469570d0, &
            -28.050500d0,     -27.653986d0,     -27.226107d0,     -26.695081d0,     -26.068830d0,     -25.427535d0, &
            -24.824593d0,     -24.276597d0,     -23.784074d0,     -23.341741d0,     -22.941682d0,     -22.574405d0, &
            -22.230015d0,     -21.900279d0,     -21.580821d0,     -21.271724d0,     -20.975977d0,     -20.697095d0, &
            -20.437525d0,     -20.198210d0,     -19.978776d0,     -19.778010d0,     -19.594395d0,     -19.426402d0, &
            -19.272658d0,     -19.132018d0,     -19.003506d0,     -18.886279d0,     -18.779471d0,     -18.682229d0, &
            -18.593791d0,     -18.513486d0,     -34.699530d0,     -33.664577d0,     -32.742196d0,     -31.920141d0, &
            -31.187494d0,     -30.534517d0,     -29.952485d0,     -29.433445d0,     -28.969508d0,     -28.550043d0, &
            -28.151428d0,     -27.716106d0,     -27.171283d0,     -26.532450d0,     -25.883271d0,     -25.274945d0, &
            -24.721674d0,     -24.221612d0,     -23.766064d0,     -23.342412d0,     -22.936868d0,     -22.539347d0, &
            -22.147839d0,     -21.767680d0,     -21.406596d0,     -21.070535d0,     -20.762353d0,     -20.482304d0, &
            -20.228997d0,     -20.000144d0,     -19.793117d0,     -19.605399d0,     -19.434764d0,     -19.279358d0, &
            -19.137702d0,     -19.008596d0,     -18.891049d0,     -18.784095d0,     -18.686816d0,     -18.598411d0, &
            -18.518182d0,     -35.199529d0,     -34.164577d0,     -33.242196d0,     -32.420141d0,     -31.687494d0, &
            -31.034517d0,     -30.452485d0,     -29.933442d0,     -29.469486d0,     -29.049882d0,     -28.650519d0, &
            -28.212565d0,     -27.662966d0,     -27.019820d0,     -26.367639d0,     -25.756129d0,     -25.196762d0, &
            -24.682684d0,     -24.197381d0,     -23.721274d0,     -23.242823d0,     -22.766150d0,     -22.305146d0, &
            -21.872571d0,     -21.475394d0,     -21.115568d0,     -20.791992d0,     -20.501992d0,     -20.242233d0, &
            -20.009186d0,     -19.799438d0,     -19.609965d0,     -19.438212d0,     -19.282107d0,     -19.140027d0, &
            -19.010673d0,     -18.892993d0,     -18.785978d0,     -18.688685d0,     -18.600295d0,     -18.520101d0, &
            -35.699529d0,     -34.664577d0,     -33.742196d0,     -32.920141d0,     -32.187494d0,     -31.534516d0, &
            -30.952485d0,     -30.433441d0,     -29.969478d0,     -29.549828d0,     -29.150219d0,     -28.711401d0, &
            -28.160217d0,     -27.515458d0,     -26.861220d0,     -26.244514d0,     -25.670269d0,     -25.119869d0, &
            -24.566373d0,     -23.997600d0,     -23.427951d0,     -22.881694d0,     -22.375346d0,     -21.915325d0, &
            -21.501854d0,     -21.132274d0,     -20.802745d0,     -20.509037d0,     -20.246929d0,     -20.012378d0, &
            -19.801663d0,     -19.611570d0,     -19.439424d0,     -19.283074d0,     -19.140845d0,     -19.011405d0, &
            -18.893677d0,     -18.786642d0,     -18.689344d0,     -18.600961d0,     -18.520780d0,     -36.199529d0, &
            -35.164577d0,     -34.242196d0,     -33.420141d0,     -32.687494d0,     -32.034516d0,     -31.452485d0, &
            -30.933441d0,     -30.469476d0,     -30.049811d0,     -29.650123d0,     -29.211022d0,     -28.659241d0, &
            -28.013274d0,     -27.354971d0,     -26.723937d0,     -26.109116d0,     -25.476248d0,     -24.812192d0, &
            -24.144864d0,     -23.510014d0,     -22.926947d0,     -22.400888d0,     -21.930257d0,     -21.510900d0, &
            -21.137931d0,     -20.806378d0,     -20.511425d0,     -20.248533d0,     -20.013481d0,     -19.802444d0, &
            -19.612144d0,     -19.439866d0,     -19.283434d0,     -19.141154d0,     -19.011685d0,     -18.893943d0, &
            -18.786901d0,     -18.689604d0,     -18.601224d0,     -18.521050d0,     -36.699529d0,     -35.664577d0, &
            -34.742196d0,     -33.920141d0,     -33.187494d0,     -32.534516d0,     -31.952485d0,     -31.433441d0, &
            -30.969475d0,     -30.549806d0,     -30.150092d0,     -29.710882d0,     -29.158621d0,     -28.510092d0, &
            -27.840032d0,     -27.168788d0,     -26.461535d0,     -25.701779d0,     -24.933038d0,     -24.204968d0, &
            -23.540132d0,     -22.942749d0,     -22.409668d0,     -21.935421d0,     -21.514098d0,     -21.140004d0, &
            -20.807778d0,     -20.512406d0,     -20.249246d0,     -20.014018d0,     -19.802864d0,     -19.612486d0, &
            -19.440157d0,     -19.283691d0,     -19.141390d0,     -19.011910d0,     -18.894163d0,     -18.787122d0, &
            -18.689829d0,     -18.601456d0,     -18.521291d0,     -37.199529d0,     -36.164577d0,     -35.242196d0, &
            -34.420141d0,     -33.687494d0,     -33.034516d0,     -32.452485d0,     -31.933441d0,     -31.469475d0, &
            -31.049804d0,     -30.650081d0,     -30.210775d0,     -29.657454d0,     -29.001368d0,     -28.297348d0, &
            -27.531360d0,     -26.681766d0,     -25.808482d0,     -24.980803d0,     -24.227066d0,     -23.551248d0, &
            -22.948934d0,     -22.413477d0,     -21.937994d0,     -21.515977d0,     -21.141462d0,     -20.808963d0, &
            -20.513402d0,     -20.250104d0,     -20.014772d0,     -19.803537d0,     -19.613098d0,     -19.440721d0, &
            -19.284221d0,     -19.141899d0,     -19.012408d0,     -18.894660d0,     -18.787626d0,     -18.690349d0, &
            -18.601998d0,     -18.521860d0,     -37.699529d0,     -36.664577d0,     -35.742196d0,     -34.920141d0, &
            -34.187494d0,     -33.534516d0,     -32.952485d0,     -32.433441d0,     -31.969475d0,     -31.549804d0, &
            -31.150071d0,     -30.710547d0,     -30.154087d0,     -29.475655d0,     -28.686665d0,     -27.765668d0, &
            -26.787473d0,     -25.852571d0,     -25.001000d0,     -24.238037d0,     -23.558358d0,     -22.954231d0, &
            -22.417802d0,     -21.941718d0,     -21.519279d0,     -21.144438d0,     -20.811667d0,     -20.515870d0, &
            -20.252364d0,     -20.016847d0,     -19.805451d0,     -19.614874d0,     -19.442387d0,     -19.285804d0, &
            -19.143427d0,     -19.013910d0,     -18.896163d0,     -18.789155d0,     -18.691927d0,     -18.603643d0, &
            -18.523587d0 &
        ], [NCOOLT,NCOOLNE] )

   ! OI
   real*8, parameter :: cool_logLrem_OI(NCOOLT,NCOOLNE) = reshape( [ &
            -30.434912d0,     -29.340103d0,     -28.361459d0,     -27.486185d0,     -26.702860d0,     -26.001312d0, &
            -25.372611d0,     -24.808873d0,     -24.303102d0,     -23.849095d0,     -23.441352d0,     -23.074995d0, &
            -22.745687d0,     -22.449409d0,     -22.182426d0,     -21.941430d0,     -21.723512d0,     -21.526131d0, &
            -21.347130d0,     -21.184682d0,     -21.037128d0,     -20.902905d0,     -20.780475d0,     -20.668208d0, &
            -20.564562d0,     -20.468146d0,     -20.377748d0,     -20.292333d0,     -20.211106d0,     -20.133531d0, &
            -20.059304d0,     -19.988318d0,     -19.920578d0,     -19.856160d0,     -19.794868d0,     -19.736323d0, &
            -19.680292d0,     -19.626709d0,     -19.575652d0,     -19.527292d0,     -19.481832d0,     -30.434910d0, &
            -29.340102d0,     -28.361458d0,     -27.486184d0,     -26.702859d0,     -26.001311d0,     -25.372610d0, &
            -24.808872d0,     -24.303101d0,     -23.849095d0,     -23.441352d0,     -23.074994d0,     -22.745687d0, &
            -22.449409d0,     -22.182426d0,     -21.941430d0,     -21.723512d0,     -21.526131d0,     -21.347130d0, &
            -21.184683d0,     -21.037128d0,     -20.902906d0,     -20.780475d0,     -20.668209d0,     -20.564563d0, &
            -20.468147d0,     -20.377749d0,     -20.292334d0,     -20.211107d0,     -20.133531d0,     -20.059305d0, &
            -19.988318d0,     -19.920579d0,     -19.856161d0,     -19.794868d0,     -19.736324d0,     -19.680292d0, &
            -19.626710d0,     -19.575653d0,     -19.527293d0,     -19.481833d0,     -30.434905d0,     -29.340097d0, &
            -28.361454d0,     -27.486180d0,     -26.702855d0,     -26.001308d0,     -25.372608d0,     -24.808870d0, &
            -24.303099d0,     -23.849093d0,     -23.441350d0,     -23.074993d0,     -22.745686d0,     -22.449409d0, &
            -22.182426d0,     -21.941430d0,     -21.723513d0,     -21.526132d0,     -21.347131d0,     -21.184684d0, &
            -21.037129d0,     -20.902907d0,     -20.780477d0,     -20.668210d0,     -20.564565d0,     -20.468148d0, &
            -20.377750d0,     -20.292336d0,     -20.211109d0,     -20.133534d0,     -20.059307d0,     -19.988321d0, &
            -19.920581d0,     -19.856163d0,     -19.794871d0,     -19.736326d0,     -19.680294d0,     -19.626712d0, &
            -19.575655d0,     -19.527295d0,     -19.481835d0,     -30.434888d0,     -29.340083d0,     -28.361441d0, &
            -27.486168d0,     -26.702845d0,     -26.001300d0,     -25.372600d0,     -24.808863d0,     -24.303094d0, &
            -23.849089d0,     -23.441347d0,     -23.074991d0,     -22.745685d0,     -22.449408d0,     -22.182425d0, &
            -21.941430d0,     -21.723514d0,     -21.526134d0,     -21.347133d0,     -21.184687d0,     -21.037133d0, &
            -20.902912d0,     -20.780481d0,     -20.668215d0,     -20.564570d0,     -20.468154d0,     -20.377757d0, &
            -20.292342d0,     -20.211116d0,     -20.133540d0,     -20.059314d0,     -19.988327d0,     -19.920588d0, &
            -19.856170d0,     -19.794877d0,     -19.736332d0,     -19.680300d0,     -19.626718d0,     -19.575661d0, &
            -19.527301d0,     -19.481841d0,     -30.434836d0,     -29.340036d0,     -28.361400d0,     -27.486132d0, &
            -26.702813d0,     -26.001272d0,     -25.372576d0,     -24.808843d0,     -24.303077d0,     -23.849074d0, &
            -23.441336d0,     -23.074983d0,     -22.745679d0,     -22.449405d0,     -22.182425d0,     -21.941432d0, &
            -21.723518d0,     -21.526140d0,     -21.347142d0,     -21.184697d0,     -21.037145d0,     -20.902925d0, &
            -20.780496d0,     -20.668232d0,     -20.564588d0,     -20.468173d0,     -20.377777d0,     -20.292363d0, &
            -20.211137d0,     -20.133561d0,     -20.059335d0,     -19.988348d0,     -19.920608d0,     -19.856190d0, &
            -19.794897d0,     -19.736352d0,     -19.680320d0,     -19.626737d0,     -19.575680d0,     -19.527320d0, &
            -19.481860d0,     -30.434672d0,     -29.339891d0,     -28.361271d0,     -27.486019d0,     -26.702714d0, &
            -26.001185d0,     -25.372501d0,     -24.808778d0,     -24.303022d0,     -23.849030d0,     -23.441301d0, &
            -23.074957d0,     -22.745662d0,     -22.449396d0,     -22.182423d0,     -21.941438d0,     -21.723531d0, &
            -21.526159d0,     -21.347167d0,     -21.184728d0,     -21.037182d0,     -20.902968d0,     -20.780544d0, &
            -20.668284d0,     -20.564644d0,     -20.468233d0,     -20.377839d0,     -20.292428d0,     -20.211203d0, &
            -20.133629d0,     -20.059402d0,     -19.988415d0,     -19.920675d0,     -19.856255d0,     -19.794961d0, &
            -19.736415d0,     -19.680383d0,     -19.626799d0,     -19.575742d0,     -19.527380d0,     -19.481920d0, &
            -30.434168d0,     -29.339444d0,     -28.360877d0,     -27.485671d0,     -26.702408d0,     -26.000918d0, &
            -25.372270d0,     -24.808581d0,     -24.302858d0,     -23.848896d0,     -23.441196d0,     -23.074880d0, &
            -22.745612d0,     -22.449370d0,     -22.182422d0,     -21.941459d0,     -21.723574d0,     -21.526222d0, &
            -21.347250d0,     -21.184831d0,     -21.037302d0,     -20.903105d0,     -20.780697d0,     -20.668452d0, &
            -20.564825d0,     -20.468425d0,     -20.378040d0,     -20.292635d0,     -20.211414d0,     -20.133842d0, &
            -20.059616d0,     -19.988628d0,     -19.920885d0,     -19.856463d0,     -19.795166d0,     -19.736616d0, &
            -19.680581d0,     -19.626994d0,     -19.575934d0,     -19.527571d0,     -19.482109d0,     -30.432698d0, &
            -29.338146d0,     -28.359731d0,     -27.484663d0,     -26.701526d0,     -26.000152d0,     -25.371611d0, &
            -24.808024d0,     -24.302396d0,     -23.848525d0,     -23.440913d0,     -23.074680d0,     -22.745490d0, &
            -22.449324d0,     -22.182447d0,     -21.941552d0,     -21.723732d0,     -21.526444d0,     -21.347532d0, &
            -21.185171d0,     -21.037698d0,     -20.903553d0,     -20.781194d0,     -20.668994d0,     -20.565407d0, &
            -20.469042d0,     -20.378685d0,     -20.293300d0,     -20.212093d0,     -20.134528d0,     -20.060303d0, &
            -19.989311d0,     -19.921560d0,     -19.857127d0,     -19.795819d0,     -19.737258d0,     -19.681211d0, &
            -19.627613d0,     -19.576542d0,     -19.528168d0,     -19.482695d0,     -30.429033d0,     -29.334925d0, &
            -28.356912d0,     -27.482208d0,     -26.699403d0,     -25.998336d0,     -25.370082d0,     -24.806764d0, &
            -24.301391d0,     -23.847764d0,     -23.440384d0,     -23.074373d0,     -22.745395d0,     -22.449431d0, &
            -22.182749d0,     -21.942043d0,     -21.724406d0,     -21.527297d0,     -21.348561d0,     -21.186370d0, &
            -21.039061d0,     -20.905073d0,     -20.782861d0,     -20.670795d0,     -20.567328d0,     -20.471064d0, &
            -20.380788d0,     -20.295463d0,     -20.214294d0,     -20.136745d0,     -20.062518d0,     -19.991508d0, &
            -19.923728d0,     -19.859258d0,     -19.797909d0,     -19.739303d0,     -19.683209d0,     -19.629562d0, &
            -19.578442d0,     -19.530017d0,     -19.484493d0,     -30.422717d0,     -29.329497d0,     -28.352296d0, &
            -27.478340d0,     -26.696229d0,     -25.995810d0,     -25.368165d0,     -24.805424d0,     -24.300601d0, &
            -23.847501d0,     -23.440627d0,     -23.075105d0,     -22.746600d0,     -22.451097d0,     -22.184867d0, &
            -21.944611d0,     -21.727424d0,     -21.530767d0,     -21.352483d0,     -21.190740d0,     -21.043870d0, &
            -20.910304d0,     -20.788489d0,     -20.676790d0,     -20.573651d0,     -20.477667d0,     -20.387611d0, &
            -20.302441d0,     -20.221361d0,     -20.143839d0,     -20.069582d0,     -19.998493d0,     -19.930596d0, &
            -19.865984d0,     -19.804476d0,     -19.745698d0,     -19.689421d0,     -19.635585d0,     -19.584271d0, &
            -19.535652d0,     -19.489933d0,     -30.418070d0,     -29.326140d0,     -28.350152d0,     -27.477344d0, &
            -26.696327d0,     -25.996960d0,     -25.370335d0,     -24.808587d0,     -24.304738d0,     -23.852597d0, &
            -23.446672d0,     -23.082093d0,     -22.754527d0,     -22.459966d0,     -22.194694d0,     -21.955418d0, &
            -21.739242d0,     -21.543629d0,     -21.366418d0,     -21.205764d0,     -21.059983d0,     -20.927484d0, &
            -20.806689d0,     -20.695938d0,     -20.593645d0,     -20.498368d0,     -20.408847d0,     -20.324019d0, &
            -20.243083d0,     -20.165514d0,     -20.091039d0,     -20.019591d0,     -19.951230d0,     -19.886086d0, &
            -19.824003d0,     -19.764618d0,     -19.707712d0,     -19.653236d0,     -19.601282d0,     -19.552032d0, &
            -19.505693d0,     -30.426671d0,     -29.336393d0,     -28.362041d0,     -27.490865d0,     -26.711491d0, &
            -26.013795d0,     -25.388877d0,     -24.828880d0,     -24.326830d0,     -23.876539d0,     -23.472519d0, &
            -23.109897d0,     -22.784342d0,     -22.491852d0,     -22.228729d0,     -21.991695d0,     -21.777863d0, &
            -21.584696d0,     -21.410017d0,     -21.251946d0,     -21.108748d0,     -20.978762d0,     -20.860329d0, &
            -20.751699d0,     -20.651193d0,     -20.557267d0,     -20.468571d0,     -20.383997d0,     -20.302752d0, &
            -20.224361d0,     -20.148634d0,     -20.075609d0,     -20.005467d0,     -19.938450d0,     -19.874464d0, &
            -19.813172d0,     -19.754383d0,     -19.698074d0,     -19.644359d0,     -19.593430d0,     -19.545495d0, &
            -30.473382d0,     -29.385538d0,     -28.413746d0,     -27.545277d0,     -26.768780d0,     -26.074151d0, &
            -25.452500d0,     -24.895969d0,     -24.397580d0,     -23.951136d0,     -23.551135d0,     -23.192687d0, &
            -22.871435d0,     -22.583374d0,     -22.324817d0,     -22.092491d0,     -21.883491d0,     -21.695247d0, &
            -21.525507d0,     -21.372269d0,     -21.233644d0,     -21.107781d0,     -20.992801d0,     -20.886739d0, &
            -20.787700d0,     -20.693938d0,     -20.603980d0,     -20.516736d0,     -20.431585d0,     -20.348321d0, &
            -20.267081d0,     -20.188236d0,     -20.112277d0,     -20.039695d0,     -19.970495d0,     -19.904359d0, &
            -19.841106d0,     -19.780725d0,     -19.723331d0,     -19.669104d0,     -19.618219d0,     -30.604800d0, &
            -29.520930d0,     -28.553403d0,     -27.689510d0,     -26.917923d0,     -26.228564d0,     -25.612532d0, &
            -25.061938d0,     -24.569769d0,     -24.129785d0,     -23.736429d0,     -23.384742d0,     -23.070283d0, &
            -22.788980d0,     -22.537098d0,     -22.311269d0,     -22.108456d0,     -21.925913d0,     -21.761158d0, &
            -21.611898d0,     -21.475916d0,     -21.351003d0,     -21.234898d0,     -21.125320d0,     -21.020153d0, &
            -20.917596d0,     -20.816389d0,     -20.715993d0,     -20.616564d0,     -20.518714d0,     -20.423303d0, &
            -20.331244d0,     -20.243383d0,     -20.160374d0,     -20.082150d0,     -20.008220d0,     -19.938265d0, &
            -19.872162d0,     -19.809928d0,     -19.751645d0,     -19.697392d0,     -30.868714d0,     -29.789561d0, &
            -28.827033d0,     -27.968427d0,     -27.202423d0,     -26.518940d0,     -25.909039d0,     -25.364775d0, &
            -24.879067d0,     -24.445598d0,     -24.058717d0,     -23.713335d0,     -23.404831d0,     -23.128922d0, &
            -22.881612d0,     -22.659184d0,     -22.458189d0,     -22.275448d0,     -22.108047d0,     -21.953287d0, &
            -21.808613d0,     -21.671545d0,     -21.539675d0,     -21.410837d0,     -21.283396d0,     -21.156406d0, &
            -21.029803d0,     -20.904419d0,     -20.781614d0,     -20.662798d0,     -20.549203d0,     -20.441776d0, &
            -20.341177d0,     -20.247758d0,     -20.161034d0,     -20.080118d0,     -20.004408d0,     -19.933577d0, &
            -19.867482d0,     -19.806080d0,     -19.749341d0,     -31.256487d0,     -30.180582d0,     -29.221441d0, &
            -28.366352d0,     -27.603989d0,     -26.924263d0,     -26.318199d0,     -25.777801d0,     -25.295930d0, &
            -24.866189d0,     -24.482787d0,     -24.140405d0,     -23.834044d0,     -23.558881d0,     -23.310209d0, &
            -23.083494d0,     -22.874509d0,     -22.679529d0,     -22.495442d0,     -22.319706d0,     -22.150209d0, &
            -21.985086d0,     -21.822691d0,     -21.661875d0,     -21.502298d0,     -21.344414d0,     -21.189452d0, &
            -21.039159d0,     -20.895270d0,     -20.759078d0,     -20.631416d0,     -20.512722d0,     -20.403142d0, &
            -20.302565d0,     -20.210061d0,     -20.124394d0,     -20.044735d0,     -19.970605d0,     -19.901760d0, &
            -19.838075d0,     -19.779456d0,     -31.714019d0,     -30.639568d0,     -29.681926d0,     -28.828376d0, &
            -28.067590d0,     -27.389469d0,     -26.785021d0,     -26.246219d0,     -25.765868d0,     -25.337458d0, &
            -24.954972d0,     -24.612666d0,     -24.304836d0,     -24.025666d0,     -23.769275d0,     -23.530083d0, &
            -23.303366d0,     -23.085736d0,     -22.875139d0,     -22.670316d0,     -22.470154d0,     -22.273328d0, &
            -22.078564d0,     -21.885364d0,     -21.694447d0,     -21.507534d0,     -21.326900d0,     -21.154803d0, &
            -20.992934d0,     -20.842171d0,     -20.702803d0,     -20.574728d0,     -20.457612d0,     -20.350934d0, &
            -20.253391d0,     -20.163459d0,     -20.080139d0,     -20.002850d0,     -19.931278d0,     -19.865245d0, &
            -19.804615d0,     -32.199681d0,     -31.125752d0,     -30.168646d0,     -29.315644d0,     -28.555416d0, &
            -27.877860d0,     -27.273973d0,     -26.735702d0,     -26.255795d0,     -25.827604d0,     -25.444831d0, &
            -25.101215d0,     -24.790249d0,     -24.505090d0,     -24.238891d0,     -23.985619d0,     -23.740963d0, &
            -23.502604d0,     -23.269409d0,     -23.040064d0,     -22.812340d0,     -22.583676d0,     -22.352840d0, &
            -22.121267d0,     -21.892689d0,     -21.671579d0,     -21.461777d0,     -21.265861d0,     -21.085033d0, &
            -20.919343d0,     -20.768251d0,     -20.630945d0,     -20.506517d0,     -20.393994d0,     -20.291677d0, &
            -20.197747d0,     -20.111029d0,     -20.030837d0,     -19.956785d0,     -19.888637d0,     -19.826210d0, &
            -32.695047d0,     -31.621290d0,     -30.664361d0,     -29.811540d0,     -29.051494d0,     -28.374123d0, &
            -27.770414d0,     -27.232300d0,     -26.752485d0,     -26.324221d0,     -25.941009d0,     -25.596235d0, &
            -25.282884d0,     -24.993534d0,     -24.720907d0,     -24.458892d0,     -24.203288d0,     -23.951269d0, &
            -23.699569d0,     -23.443104d0,     -23.176364d0,     -22.897488d0,     -22.611068d0,     -22.326225d0, &
            -22.052021d0,     -21.794700d0,     -21.557428d0,     -21.341113d0,     -21.145216d0,     -20.968333d0, &
            -20.808837d0,     -20.665149d0,     -20.535815d0,     -20.419472d0,     -20.314110d0,     -20.217688d0, &
            -20.128900d0,     -20.046975d0,     -19.971471d0,     -19.902111d0,     -19.838675d0,     -33.193571d0, &
            -32.119870d0,     -31.162997d0,     -30.310233d0,     -29.550246d0,     -28.872934d0,     -28.269281d0, &
            -27.731212d0,     -27.251411d0,     -26.823084d0,     -26.439639d0,     -26.094308d0,     -25.779860d0, &
            -25.488627d0,     -25.213079d0,     -24.946599d0,     -24.683299d0,     -24.415879d0,     -24.133354d0, &
            -23.824177d0,     -23.486283d0,     -23.132115d0,     -22.780452d0,     -22.446227d0,     -22.137330d0, &
            -21.856387d0,     -21.603102d0,     -21.375826d0,     -21.172292d0,     -20.989968d0,     -20.826514d0, &
            -20.679888d0,     -20.548336d0,     -20.430294d0,     -20.323595d0,     -20.226093d0,     -20.136417d0, &
            -20.053757d0,     -19.977644d0,     -19.907782d0,     -19.843932d0,     -33.693104d0,     -32.619420d0, &
            -31.662565d0,     -30.809820d0,     -30.049851d0,     -29.372557d0,     -28.768921d0,     -28.230866d0, &
            -27.751067d0,     -27.322714d0,     -26.939181d0,     -26.593645d0,     -26.278789d0,     -25.986799d0, &
            -25.709666d0,     -25.438865d0,     -25.162486d0,     -24.860659d0,     -24.510387d0,     -24.108129d0, &
            -23.679516d0,     -23.256781d0,     -22.860335d0,     -22.498583d0,     -22.172896d0,     -21.881517d0, &
            -21.621550d0,     -21.389852d0,     -21.183294d0,     -20.998834d0,     -20.833828d0,     -20.686048d0, &
            -20.553620d0,     -20.434901d0,     -20.327667d0,     -20.229727d0,     -20.139689d0,     -20.056727d0, &
            -19.980363d0,     -19.910292d0,     -19.846268d0,     -34.192956d0,     -33.119277d0,     -32.162428d0, &
            -31.309689d0,     -30.549726d0,     -29.872438d0,     -29.268808d0,     -28.730756d0,     -28.250958d0, &
            -27.822596d0,     -27.439033d0,     -27.093426d0,     -26.778407d0,     -26.485913d0,     -26.206580d0, &
            -25.926143d0,     -25.616544d0,     -25.236811d0,     -24.774112d0,     -24.269955d0,     -23.774301d0, &
            -23.313680d0,     -22.896570d0,     -22.523274d0,     -22.190808d0,     -21.895212d0,     -21.632469d0, &
            -21.398847d0,     -21.190893d0,     -21.005378d0,     -20.839549d0,     -20.691112d0,     -20.558151d0, &
            -20.438995d0,     -20.331394d0,     -20.233139d0,     -20.142827d0,     -20.059628d0,     -19.983058d0, &
            -19.912810d0,     -19.848635d0,     -34.692909d0,     -33.619232d0,     -32.662385d0,     -31.809647d0, &
            -31.049686d0,     -30.372400d0,     -29.768772d0,     -29.230722d0,     -28.750923d0,     -28.322558d0, &
            -27.938985d0,     -27.593350d0,     -27.278201d0,     -26.984856d0,     -26.700262d0,     -26.395263d0, &
            -26.011031d0,     -25.506513d0,     -24.931444d0,     -24.360790d0,     -23.831718d0,     -23.354224d0, &
            -22.927888d0,     -22.549000d0,     -22.212776d0,     -21.914428d0,     -21.649534d0,     -21.414155d0, &
            -21.204719d0,     -21.017927d0,     -20.850984d0,     -20.701572d0,     -20.567761d0,     -20.447867d0, &
            -20.339615d0,     -20.240774d0,     -20.149935d0,     -20.066265d0,     -19.989280d0,     -19.918667d0, &
            -19.854173d0,     -35.192894d0,     -34.119218d0,     -33.162371d0,     -32.309634d0,     -31.549674d0, &
            -30.872388d0,     -30.268760d0,     -29.730711d0,     -29.250912d0,     -28.822546d0,     -28.438969d0, &
            -28.093311d0,     -27.777948d0,     -27.482738d0,     -27.185948d0,     -26.828158d0,     -26.323878d0, &
            -25.695245d0,     -25.049009d0,     -24.445301d0,     -23.899986d0,     -23.413119d0,     -22.980426d0, &
            -22.596671d0,     -22.256406d0,     -21.954532d0,     -21.686485d0,     -21.448254d0,     -21.236224d0, &
            -21.047059d0,     -20.877951d0,     -20.726579d0,     -20.591014d0,     -20.469568d0,     -20.359924d0, &
            -20.259810d0,     -20.167809d0,     -20.083091d0,     -20.005173d0,     -19.933740d0,     -19.868527d0, &
            -35.692889d0,     -34.619214d0,     -33.662367d0,     -32.809630d0,     -32.049670d0,     -31.372384d0, &
            -30.768757d0,     -30.230707d0,     -29.750909d0,     -29.322543d0,     -28.938963d0,     -28.593267d0, &
            -28.277447d0,     -27.978021d0,     -27.654340d0,     -27.201028d0,     -26.557994d0,     -25.842872d0, &
            -25.159951d0,     -24.540021d0,     -23.985736d0,     -23.492670d0,     -23.055019d0,     -22.666964d0, &
            -22.322771d0,     -22.017204d0,     -21.745644d0,     -21.504073d0,     -21.288868d0,     -21.096684d0, &
            -20.924732d0,     -20.770719d0,     -20.632749d0,     -20.509155d0,     -20.397558d0,     -20.295629d0, &
            -20.201946d0,     -20.115702d0,     -20.036429d0,     -19.963815d0,     -19.897583d0,     -36.192888d0, &
            -35.119212d0,     -34.162365d0,     -33.309629d0,     -32.549669d0,     -31.872383d0,     -31.268756d0, &
            -30.730706d0,     -30.250908d0,     -29.822541d0,     -29.438956d0,     -29.093165d0,     -28.776117d0, &
            -28.465444d0,     -28.077933d0,     -27.471013d0,     -26.703225d0,     -25.938569d0,     -25.238461d0, &
            -24.611339d0,     -24.053066d0,     -23.557193d0,     -23.117255d0,     -22.727182d0,     -22.381087d0, &
            -22.073649d0,     -21.800219d0,     -21.556779d0,     -21.339715d0,     -21.145683d0,     -20.971924d0, &
            -20.816196d0,     -20.676654d0,     -20.551673d0,     -20.438821d0,     -20.335711d0,     -20.240932d0, &
            -20.153712d0,     -20.073618d0,     -20.000354d0,     -19.933645d0,     -36.692888d0,     -35.619212d0, &
            -34.662365d0,     -33.809628d0,     -33.049668d0,     -32.372383d0,     -31.768755d0,     -31.230706d0, &
            -30.750907d0,     -30.322540d0,     -29.938939d0,     -29.592860d0,     -29.272077d0,     -28.929009d0, &
            -28.403326d0,     -27.616643d0,     -26.769942d0,     -25.982628d0,     -25.275615d0,     -24.645845d0, &
            -24.086239d0,     -23.589520d0,     -23.148954d0,     -22.758366d0,     -22.411797d0,     -22.103880d0, &
            -21.829945d0,     -21.585975d0,     -21.368356d0,     -21.173745d0,     -20.999402d0,     -20.843116d0, &
            -20.703082d0,     -20.577704d0,     -20.464520d0,     -20.361110d0,     -20.266074d0,     -20.178668d0, &
            -20.098489d0,     -20.025262d0,     -19.958718d0,     -37.192887d0,     -36.119212d0,     -35.162365d0, &
            -34.309628d0,     -33.549668d0,     -32.872383d0,     -32.268755d0,     -31.730706d0,     -31.250907d0, &
            -30.822538d0,     -30.438888d0,     -30.091903d0,     -29.759604d0,     -29.330967d0,     -28.593277d0, &
            -27.676428d0,     -26.794747d0,     -25.998966d0,     -25.289521d0,     -24.658863d0,     -24.098844d0, &
            -23.601887d0,     -23.161167d0,     -22.770468d0,     -22.423806d0,     -22.115795d0,     -21.841755d0, &
            -21.597669d0,     -21.379923d0,     -21.185174d0,     -21.010689d0,     -20.854272d0,     -20.714132d0, &
            -20.588689d0,     -20.475469d0,     -20.372037d0,     -20.277001d0,     -20.189629d0,     -20.109534d0, &
            -20.036451d0,     -19.970116d0,     -37.692887d0,     -36.619212d0,     -35.662365d0,     -34.809628d0, &
            -34.049668d0,     -33.372382d0,     -32.768755d0,     -32.230706d0,     -31.750907d0,     -31.322532d0, &
            -30.938728d0,     -30.588900d0,     -30.222472d0,     -29.616846d0,     -28.676740d0,     -27.697919d0, &
            -26.803610d0,     -26.004955d0,     -25.294685d0,     -24.663718d0,     -24.103549d0,     -23.606505d0, &
            -23.165726d0,     -22.774987d0,     -22.428292d0,     -22.120250d0,     -21.846178d0,     -21.602059d0, &
            -21.384279d0,     -21.189494d0,     -21.014974d0,     -20.858528d0,     -20.718370d0,     -20.592926d0, &
            -20.479717d0,     -20.376305d0,     -20.281300d0,     -20.193975d0,     -20.113948d0,     -20.040958d0, &
            -19.974743d0 &
        ], [NCOOLT,NCOOLNE] )

   ! OII
   real*8, parameter :: cool_logLrem_OII(NCOOLT,NCOOLNE) = reshape( [ &
            -35.073807d0,     -33.276263d0,     -31.676879d0,     -30.254098d0,     -28.988702d0,     -27.863563d0, &
            -26.863407d0,     -25.974614d0,     -25.185043d0,     -24.483870d0,     -23.861435d0,     -23.309121d0, &
            -22.819248d0,     -22.384963d0,     -22.000084d0,     -21.659040d0,     -21.356824d0,     -21.089019d0, &
            -20.851824d0,     -20.641923d0,     -20.456298d0,     -20.292230d0,     -20.147358d0,     -20.019559d0, &
            -19.906831d0,     -19.807217d0,     -19.718692d0,     -19.638990d0,     -19.565403d0,     -19.494868d0, &
            -19.424208d0,     -19.350669d0,     -19.272394d0,     -19.188555d0,     -19.099345d0,     -19.006047d0, &
            -18.910679d0,     -18.815402d0,     -18.722183d0,     -18.632608d0,     -18.547805d0,     -35.074379d0, &
            -33.276803d0,     -31.677389d0,     -30.254579d0,     -28.989156d0,     -27.863992d0,     -26.863812d0, &
            -25.974996d0,     -25.185404d0,     -24.484211d0,     -23.861757d0,     -23.309424d0,     -22.819534d0, &
            -22.385232d0,     -22.000336d0,     -21.659276d0,     -21.357045d0,     -21.089226d0,     -20.852017d0, &
            -20.642103d0,     -20.456467d0,     -20.292390d0,     -20.147510d0,     -20.019705d0,     -19.906973d0, &
            -19.807355d0,     -19.718828d0,     -19.639125d0,     -19.565536d0,     -19.494999d0,     -19.424336d0, &
            -19.350792d0,     -19.272510d0,     -19.188662d0,     -19.099443d0,     -19.006136d0,     -18.910758d0, &
            -18.815472d0,     -18.722245d0,     -18.632662d0,     -18.547852d0,     -35.076163d0,     -33.278489d0, &
            -31.678981d0,     -30.256083d0,     -28.990578d0,     -27.865335d0,     -26.865081d0,     -25.976196d0, &
            -25.186538d0,     -24.485282d0,     -23.862768d0,     -23.310379d0,     -22.820433d0,     -22.386079d0, &
            -22.001131d0,     -21.660021d0,     -21.357742d0,     -21.089876d0,     -20.852624d0,     -20.642671d0, &
            -20.457000d0,     -20.292893d0,     -20.147989d0,     -20.020164d0,     -19.907418d0,     -19.807791d0, &
            -19.719257d0,     -19.639549d0,     -19.565956d0,     -19.495412d0,     -19.424738d0,     -19.351179d0, &
            -19.272875d0,     -19.189002d0,     -19.099753d0,     -19.006416d0,     -18.911007d0,     -18.815693d0, &
            -18.722439d0,     -18.632833d0,     -18.548002d0,     -35.081590d0,     -33.283625d0,     -31.683844d0, &
            -30.260686d0,     -28.994935d0,     -27.869459d0,     -26.868985d0,     -25.979891d0,     -25.190035d0, &
            -24.488591d0,     -23.865898d0,     -23.313337d0,     -22.823225d0,     -22.388709d0,     -22.003604d0, &
            -21.662341d0,     -21.359914d0,     -21.091906d0,     -20.854520d0,     -20.644445d0,     -20.458666d0, &
            -20.294466d0,     -20.149486d0,     -20.021602d0,     -19.908812d0,     -19.809155d0,     -19.720601d0, &
            -19.640878d0,     -19.567269d0,     -19.496704d0,     -19.425997d0,     -19.352390d0,     -19.274021d0, &
            -19.190066d0,     -19.100725d0,     -19.007291d0,     -18.911787d0,     -18.816383d0,     -18.723047d0, &
            -18.633367d0,     -18.548472d0,     -35.097111d0,     -33.298380d0,     -31.697869d0,     -30.274017d0, &
            -29.007606d0,     -27.881503d0,     -26.880432d0,     -25.990769d0,     -25.200372d0,     -24.498411d0, &
            -23.875220d0,     -23.322179d0,     -22.831599d0,     -22.396624d0,     -22.011067d0,     -21.669359d0, &
            -21.366498d0,     -21.098072d0,     -20.860292d0,     -20.649854d0,     -20.463752d0,     -20.299274d0, &
            -20.154067d0,     -20.026006d0,     -19.913086d0,     -19.813339d0,     -19.724725d0,     -19.644958d0, &
            -19.571303d0,     -19.500673d0,     -19.429866d0,     -19.356111d0,     -19.277540d0,     -19.193332d0, &
            -19.103707d0,     -19.009977d0,     -18.914181d0,     -18.818500d0,     -18.724910d0,     -18.635004d0, &
            -18.549912d0,     -35.137567d0,     -33.337016d0,     -31.734770d0,     -30.309264d0,     -29.041276d0, &
            -27.913670d0,     -26.911167d0,     -26.020142d0,     -25.228443d0,     -24.525233d0,     -23.900838d0, &
            -23.346622d0,     -22.854884d0,     -22.418752d0,     -22.032037d0,     -21.689173d0,     -21.385169d0, &
            -21.115630d0,     -20.876789d0,     -20.665365d0,     -20.478378d0,     -20.313137d0,     -20.167300d0, &
            -20.038749d0,     -19.925471d0,     -19.825478d0,     -19.736699d0,     -19.656808d0,     -19.583023d0, &
            -19.512205d0,     -19.441105d0,     -19.366915d0,     -19.287748d0,     -19.202798d0,     -19.112339d0, &
            -19.017743d0,     -18.921094d0,     -18.824606d0,     -18.730280d0,     -18.639720d0,     -18.554056d0, &
            -35.235138d0,     -33.430536d0,     -31.824400d0,     -30.395163d0,     -29.123602d0,     -27.992577d0, &
            -26.986807d0,     -26.092661d0,     -25.297981d0,     -24.591910d0,     -23.964752d0,     -23.407836d0, &
            -22.913423d0,     -22.474599d0,     -22.085157d0,     -21.739539d0,     -21.432785d0,     -21.160547d0, &
            -20.919104d0,     -20.705238d0,     -20.516034d0,     -20.348866d0,     -20.201430d0,     -20.071621d0, &
            -19.957411d0,     -19.856765d0,     -19.767534d0,     -19.687284d0,     -19.613111d0,     -19.541746d0, &
            -19.469821d0,     -19.394430d0,     -19.313648d0,     -19.226717d0,     -19.134056d0,     -19.037196d0, &
            -18.938336d0,     -18.839778d0,     -18.743579d0,     -18.651366d0,     -18.564268d0,     -35.442152d0, &
            -33.630884d0,     -32.018185d0,     -30.582496d0,     -29.304602d0,     -28.167373d0,     -27.155531d0, &
            -26.255448d0,     -25.454951d0,     -24.743159d0,     -24.110329d0,     -23.547732d0,     -23.047549d0, &
            -22.602794d0,     -22.207225d0,     -21.855300d0,     -21.542136d0,     -21.263498d0,     -21.015796d0, &
            -20.795968d0,     -20.601275d0,     -20.429245d0,     -20.277678d0,     -20.144502d0,     -20.027652d0, &
            -19.924982d0,     -19.834165d0,     -19.752527d0,     -19.676895d0,     -19.603719d0,     -19.529385d0, &
            -19.450814d0,     -19.366036d0,     -19.274450d0,     -19.176812d0,     -19.074998d0,     -18.971441d0, &
            -18.868600d0,     -18.768610d0,     -18.673121d0,     -18.583231d0,     -35.782562d0,     -33.965086d0, &
            -32.346079d0,     -30.903989d0,     -29.619608d0,     -28.475809d0,     -27.457310d0,     -26.550459d0, &
            -25.743031d0,     -25.024051d0,     -24.383633d0,     -23.812867d0,     -23.303739d0,     -22.849092d0, &
            -22.442603d0,     -22.078779d0,     -21.752928d0,     -21.461129d0,     -21.200179d0,     -20.967444d0, &
            -20.760621d0,     -20.577605d0,     -20.416425d0,     -20.275093d0,     -20.151473d0,     -20.043204d0, &
            -19.947611d0,     -19.861580d0,     -19.781456d0,     -19.703225d0,     -19.622912d0,     -19.537250d0, &
            -19.444344d0,     -19.343986d0,     -19.237546d0,     -19.127434d0,     -19.016380d0,     -18.906989d0, &
            -18.801422d0,     -18.701267d0,     -18.607512d0,     -36.216577d0,     -34.395735d0,     -32.773228d0, &
            -31.327499d0,     -30.039332d0,     -28.891579d0,     -27.868906d0,     -26.957551d0,     -26.145084d0, &
            -25.420192d0,     -24.772521d0,     -24.192590d0,     -23.671823d0,     -23.202660d0,     -22.778696d0, &
            -22.394753d0,     -22.046833d0,     -21.731948d0,     -21.447905d0,     -21.192996d0,     -20.965663d0, &
            -20.764301d0,     -20.587174d0,     -20.432291d0,     -20.297302d0,     -20.179459d0,     -20.075565d0, &
            -19.981893d0,     -19.894181d0,     -19.807880d0,     -19.718675d0,     -19.623241d0,     -19.519946d0, &
            -19.409141d0,     -19.292856d0,     -19.173952d0,     -19.055338d0,     -18.939623d0,     -18.828870d0, &
            -18.724512d0,     -18.627372d0,     -36.693406d0,     -34.871246d0,     -33.247351d0,     -31.800157d0, &
            -30.510423d0,     -29.360937d0,     -28.336207d0,     -27.422129d0,     -26.605643d0,     -25.874441d0, &
            -25.216843d0,     -24.621980d0,     -24.080298d0,     -23.584203d0,     -23.128418d0,     -22.709801d0, &
            -22.326671d0,     -21.978062d0,     -21.663187d0,     -21.381091d0,     -21.130432d0,     -20.909473d0, &
            -20.716143d0,     -20.548002d0,     -20.402187d0,     -20.275399d0,     -20.163865d0,     -20.063291d0, &
            -19.968901d0,     -19.875752d0,     -19.779331d0,     -19.676352d0,     -19.565453d0,     -19.447392d0, &
            -19.324582d0,     -19.200092d0,     -19.076847d0,     -18.957383d0,     -18.843638d0,     -18.736913d0, &
            -18.637903d0,     -37.185810d0,     -35.363201d0,     -33.738827d0,     -32.291108d0,     -31.000742d0, &
            -29.850340d0,     -28.823959d0,     -27.906527d0,     -27.083236d0,     -26.339212d0,     -25.659974d0, &
            -25.032994d0,     -24.449691d0,     -23.906259d0,     -23.402452d0,     -22.939347d0,     -22.517629d0, &
            -22.136987d0,     -21.796232d0,     -21.493542d0,     -21.226608d0,     -20.992819d0,     -20.789375d0, &
            -20.613234d0,     -20.461025d0,     -20.329002d0,     -20.212988d0,     -20.108315d0,     -20.009896d0, &
            -19.912569d0,     -19.811744d0,     -19.704220d0,     -19.588848d0,     -19.466642d0,     -19.340217d0, &
            -19.212719d0,     -19.087058d0,     -18.965694d0,     -18.850476d0,     -18.742617d0,     -18.642733d0, &
            -37.683380d0,     -35.860623d0,     -34.236084d0,     -32.788144d0,     -31.497363d0,     -30.345923d0, &
            -29.316790d0,     -28.392575d0,     -27.554574d0,     -26.783280d0,     -26.061803d0,     -25.380733d0, &
            -24.739546d0,     -24.142655d0,     -23.594367d0,     -23.096503d0,     -22.648423d0,     -22.247939d0, &
            -21.892162d0,     -21.577953d0,     -21.302058d0,     -21.061198d0,     -20.852103d0,     -20.671372d0, &
            -20.515334d0,     -20.379958d0,     -20.260774d0,     -20.152824d0,     -20.050763d0,     -19.949277d0, &
            -19.843795d0,     -19.731335d0,     -19.611104d0,     -19.484468d0,     -19.354265d0,     -19.223711d0, &
            -19.095683d0,     -18.972538d0,     -18.856003d0,     -18.747177d0,     -18.646584d0,     -38.182608d0, &
            -36.359803d0,     -34.735197d0,     -33.287114d0,     -31.995891d0,     -30.843003d0,     -29.809583d0, &
            -28.874461d0,     -28.013364d0,     -27.202811d0,     -26.429106d0,     -25.693125d0,     -25.003949d0, &
            -24.369658d0,     -23.793636d0,     -23.275280d0,     -22.811712d0,     -22.399123d0,     -22.033548d0, &
            -21.711159d0,     -21.428289d0,     -21.181400d0,     -20.967026d0,     -20.781573d0,     -20.621126d0, &
            -20.481314d0,     -20.357221d0,     -20.243386d0,     -20.134037d0,     -20.023682d0,     -19.908017d0, &
            -19.784779d0,     -19.654118d0,     -19.518205d0,     -19.380289d0,     -19.243654d0,     -19.111032d0, &
            -18.984518d0,     -18.865552d0,     -18.754987d0,     -18.653153d0,     -38.682363d0,     -36.859541d0, &
            -35.234905d0,     -33.786728d0,     -32.495140d0,     -31.340954d0,     -30.303548d0,     -29.358237d0, &
            -28.476363d0,     -27.633623d0,     -26.822608d0,     -26.052649d0,     -25.336734d0,     -24.682415d0, &
            -24.091136d0,     -23.560590d0,     -23.086767d0,     -22.665173d0,     -22.291430d0,     -21.961478d0, &
            -21.671542d0,     -21.418023d0,     -21.197334d0,     -21.005622d0,     -20.838484d0,     -20.690762d0, &
            -20.556480d0,     -20.429081d0,     -20.302108d0,     -20.170338d0,     -20.030919d0,     -19.883847d0, &
            -19.731452d0,     -19.577261d0,     -19.424871d0,     -19.277257d0,     -19.136572d0,     -19.004253d0, &
            -18.881150d0,     -18.767647d0,     -18.663722d0,     -39.182286d0,     -37.359458d0,     -35.734810d0, &
            -34.286586d0,     -32.994809d0,     -31.839926d0,     -30.800341d0,     -29.849449d0,     -28.956491d0, &
            -28.097556d0,     -27.269121d0,     -26.484200d0,     -25.756691d0,     -25.093396d0,     -24.494795d0, &
            -23.957885d0,     -23.478246d0,     -23.051151d0,     -22.672087d0,     -22.336910d0,     -22.041784d0, &
            -21.782978d0,     -21.556523d0,     -21.357748d0,     -21.180799d0,     -21.018415d0,     -20.862317d0, &
            -20.704553d0,     -20.539504d0,     -20.365419d0,     -20.184293d0,     -20.000259d0,     -19.817784d0, &
            -19.640585d0,     -19.471280d0,     -19.311522d0,     -19.162243d0,     -19.023880d0,     -18.896525d0, &
            -18.780018d0,     -18.673955d0,     -39.682261d0,     -37.859432d0,     -36.234779d0,     -34.786538d0, &
            -33.494688d0,     -32.339533d0,     -31.299090d0,     -30.345999d0,     -29.448726d0,     -28.583644d0, &
            -27.748785d0,     -26.958434d0,     -26.226689d0,     -25.560024d0,     -24.958563d0,     -24.419064d0, &
            -23.936968d0,     -23.507471d0,     -23.125997d0,     -22.788316d0,     -22.490380d0,     -22.227888d0, &
            -21.995539d0,     -21.786092d0,     -21.589762d0,     -21.395083d0,     -21.192124d0,     -20.976597d0, &
            -20.751088d0,     -20.522192d0,     -20.296656d0,     -20.079377d0,     -19.873235d0,     -19.679649d0, &
            -19.499130d0,     -19.331713d0,     -19.177188d0,     -19.035207d0,     -18.905336d0,     -18.787058d0, &
            -18.679738d0,     -40.182254d0,     -38.359423d0,     -36.734769d0,     -35.286522d0,     -33.994648d0, &
            -32.839400d0,     -31.798664d0,     -30.844823d0,     -29.946084d0,     -29.078932d0,     -28.241929d0, &
            -27.449775d0,     -26.716617d0,     -26.048813d0,     -25.446372d0,     -24.905971d0,     -24.423006d0, &
            -23.992635d0,     -23.610214d0,     -23.271285d0,     -22.971101d0,     -22.703506d0,     -22.459186d0, &
            -22.224262d0,     -21.982168d0,     -21.721136d0,     -21.441707d0,     -21.154665d0,     -20.872265d0, &
            -20.602729d0,     -20.349896d0,     -20.114803d0,     -19.897135d0,     -19.696069d0,     -19.510649d0, &
            -19.339979d0,     -19.183266d0,     -19.039792d0,     -18.908888d0,     -18.789887d0,     -18.682052d0, &
            -40.682251d0,     -38.859421d0,     -37.234766d0,     -35.786517d0,     -34.494635d0,     -33.339357d0, &
            -32.298526d0,     -31.344442d0,     -30.445228d0,     -29.577408d0,     -28.739716d0,     -27.946981d0, &
            -27.213368d0,     -26.545197d0,     -25.942436d0,     -25.401740d0,     -24.918483d0,     -24.487777d0, &
            -24.104795d0,     -23.764382d0,     -23.459604d0,     -23.178812d0,     -22.902570d0,     -22.606873d0, &
            -22.278521d0,     -21.927054d0,     -21.574020d0,     -21.236482d0,     -20.922486d0,     -20.633890d0, &
            -20.369608d0,     -20.127563d0,     -19.905596d0,     -19.701820d0,     -19.514656d0,     -19.342842d0, &
            -19.185365d0,     -19.041372d0,     -18.910110d0,     -18.790858d0,     -18.682844d0,     -41.182251d0, &
            -39.359420d0,     -37.734765d0,     -36.286516d0,     -34.994631d0,     -33.839343d0,     -32.798482d0, &
            -31.844321d0,     -30.944955d0,     -30.076922d0,     -29.239011d0,     -28.446092d0,     -27.712334d0, &
            -27.044046d0,     -26.441183d0,     -25.900390d0,     -25.417016d0,     -24.986055d0,     -24.602107d0, &
            -24.257852d0,     -23.939838d0,     -23.622454d0,     -23.271404d0,     -22.871141d0,     -22.443378d0, &
            -22.021769d0,     -21.627135d0,     -21.266620d0,     -20.940064d0,     -20.644476d0,     -20.376189d0, &
            -20.131778d0,     -19.908374d0,     -19.703700d0,     -19.515961d0,     -19.343772d0,     -19.186044d0, &
            -19.041881d0,     -18.910501d0,     -18.791166d0,     -18.683093d0,     -41.682250d0,     -39.859420d0, &
            -38.234765d0,     -36.786515d0,     -35.494630d0,     -34.339339d0,     -33.298468d0,     -32.344282d0, &
            -31.444869d0,     -30.576768d0,     -29.738788d0,     -28.945810d0,     -28.212007d0,     -27.543681d0, &
            -26.940786d0,     -26.399953d0,     -25.916470d0,     -25.484955d0,     -25.098232d0,     -24.742416d0, &
            -24.386914d0,     -23.984281d0,     -23.512231d0,     -23.005781d0,     -22.512408d0,     -22.056844d0, &
            -21.645510d0,     -21.276675d0,     -20.945816d0,     -20.647903d0,     -20.378307d0,     -20.133129d0, &
            -19.909260d0,     -19.704296d0,     -19.516372d0,     -19.344060d0,     -19.186249d0,     -19.042028d0, &
            -18.910607d0,     -18.791242d0,     -18.683148d0,     -42.182250d0,     -40.359420d0,     -38.734765d0, &
            -37.286515d0,     -35.994629d0,     -34.839338d0,     -33.798464d0,     -32.844270d0,     -31.944841d0, &
            -31.076720d0,     -30.238717d0,     -29.445721d0,     -28.711903d0,     -28.043565d0,     -27.440657d0, &
            -26.899785d0,     -26.416041d0,     -25.982864d0,     -25.587607d0,     -25.198251d0,     -24.753478d0, &
            -24.214913d0,     -23.628461d0,     -23.058918d0,     -22.536789d0,     -22.068578d0,     -21.651500d0, &
            -21.279911d0,     -20.947656d0,     -20.648995d0,     -20.378979d0,     -20.133554d0,     -19.909535d0, &
            -19.704475d0,     -19.516486d0,     -19.344128d0,     -19.186282d0,     -19.042034d0,     -18.910590d0, &
            -18.791206d0,     -18.683095d0,     -42.682250d0,     -40.859420d0,     -39.234765d0,     -37.786515d0, &
            -36.494629d0,     -35.339337d0,     -34.298463d0,     -33.344266d0,     -32.444833d0,     -31.576704d0, &
            -30.738695d0,     -29.945693d0,     -29.211870d0,     -28.543528d0,     -27.940608d0,     -27.399639d0, &
            -26.915095d0,     -26.476742d0,     -26.056038d0,     -25.582705d0,     -24.989654d0,     -24.323774d0, &
            -23.672821d0,     -23.077184d0,     -22.544799d0,     -22.072359d0,     -21.653413d0,     -21.280941d0, &
            -20.948240d0,     -20.649340d0,     -20.379189d0,     -20.133684d0,     -19.909613d0,     -19.704514d0, &
            -19.516494d0,     -19.344108d0,     -19.186234d0,     -19.041957d0,     -18.910483d0,     -18.791070d0, &
            -18.682929d0,     -43.182250d0,     -41.359419d0,     -39.734765d0,     -38.286515d0,     -36.994629d0, &
            -35.839337d0,     -34.798462d0,     -33.844265d0,     -32.944830d0,     -32.076700d0,     -31.238688d0, &
            -30.445684d0,     -29.711860d0,     -29.043515d0,     -28.440567d0,     -27.899298d0,     -27.412245d0, &
            -26.958063d0,     -26.469177d0,     -25.842538d0,     -25.102341d0,     -24.364784d0,     -23.687853d0, &
            -23.083128d0,     -22.547366d0,     -22.073564d0,     -21.654023d0,     -21.281270d0,     -20.948426d0, &
            -20.649451d0,     -20.379256d0,     -20.133723d0,     -19.909630d0,     -19.704513d0,     -19.516472d0, &
            -19.344062d0,     -19.186159d0,     -19.041847d0,     -18.910332d0,     -18.790872d0,     -18.682681d0, &
            -43.682250d0,     -41.859419d0,     -40.234765d0,     -38.786515d0,     -37.494629d0,     -36.339337d0, &
            -35.298462d0,     -34.344265d0,     -33.444829d0,     -32.576698d0,     -31.738685d0,     -30.945681d0, &
            -30.211856d0,     -29.543506d0,     -28.940473d0,     -28.398259d0,     -27.903394d0,     -27.403821d0, &
            -26.772411d0,     -25.972731d0,     -25.145071d0,     -24.378613d0,     -23.692726d0,     -23.085033d0, &
            -22.548189d0,     -22.073954d0,     -21.654223d0,     -21.281381d0,     -20.948492d0,     -20.649492d0, &
            -20.379282d0,     -20.133738d0,     -19.909636d0,     -19.704508d0,     -19.516454d0,     -19.344027d0, &
            -19.186103d0,     -19.041764d0,     -18.910216d0,     -18.790716d0,     -18.682476d0,     -44.182250d0, &
            -42.359419d0,     -40.734765d0,     -39.286515d0,     -37.994629d0,     -36.839337d0,     -35.798462d0, &
            -34.844264d0,     -33.944829d0,     -33.076698d0,     -32.238685d0,     -31.445680d0,     -30.711854d0, &
            -30.043486d0,     -29.440185d0,     -28.895002d0,     -28.376553d0,     -27.766715d0,     -26.939578d0, &
            -26.023637d0,     -25.159539d0,     -24.383108d0,     -23.694307d0,     -23.085664d0,     -22.548476d0, &
            -22.074102d0,     -21.654310d0,     -21.281439d0,     -20.948535d0,     -20.649525d0,     -20.379308d0, &
            -20.133757d0,     -19.909646d0,     -19.704509d0,     -19.516444d0,     -19.344004d0,     -19.186065d0, &
            -19.041709d0,     -18.910142d0,     -18.790619d0,     -18.682355d0,     -44.682250d0,     -42.859419d0, &
            -41.234765d0,     -39.786515d0,     -38.494629d0,     -37.339337d0,     -36.298462d0,     -35.344264d0, &
            -34.444829d0,     -33.576697d0,     -32.738684d0,     -31.945680d0,     -31.211851d0,     -30.543428d0, &
            -29.939282d0,     -29.384869d0,     -28.801210d0,     -27.998195d0,     -27.009477d0,     -26.041162d0, &
            -25.164314d0,     -24.384633d0,     -23.694897d0,     -23.085949d0,     -22.548647d0,     -22.074227d0, &
            -21.654413d0,     -21.281530d0,     -20.948617d0,     -20.649599d0,     -20.379371d0,     -20.133807d0, &
            -19.909680d0,     -19.704524d0,     -19.516437d0,     -19.343976d0,     -19.186015d0,     -19.041640d0, &
            -18.910057d0,     -18.790521d0,     -18.682247d0,     -45.182250d0,     -43.359419d0,     -41.734765d0, &
            -40.286515d0,     -38.994629d0,     -37.839337d0,     -36.798462d0,     -35.844264d0,     -34.944829d0, &
            -34.076697d0,     -33.238684d0,     -32.445680d0,     -31.711843d0,     -31.043246d0,     -30.436444d0, &
            -29.854340d0,     -29.124252d0,     -28.107889d0,     -27.034483d0,     -26.047181d0,     -25.166146d0, &
            -24.385411d0,     -23.695364d0,     -23.086307d0,     -22.548958d0,     -22.074512d0,     -21.654681d0, &
            -21.281782d0,     -20.948855d0,     -20.649817d0,     -20.379565d0,     -20.133967d0,     -19.909798d0, &
            -19.704593d0,     -19.516454d0,     -19.343941d0,     -19.185933d0,     -19.041519d0,     -18.909910d0, &
            -18.790362d0,     -18.682092d0,     -45.682250d0,     -43.859419d0,     -42.234765d0,     -40.786515d0, &
            -39.494629d0,     -38.339337d0,     -37.298462d0,     -36.344264d0,     -35.444829d0,     -34.576697d0, &
            -33.738684d0,     -32.945679d0,     -32.211818d0,     -31.542674d0,     -30.927631d0,     -30.270218d0, &
            -29.311661d0,     -28.150360d0,     -27.043781d0,     -26.050132d0,     -25.167706d0,     -24.386589d0, &
            -23.696399d0,     -23.087265d0,     -22.549863d0,     -22.075373d0,     -21.655501d0,     -21.282563d0, &
            -20.949591d0,     -20.650501d0,     -20.380180d0,     -20.134498d0,     -19.910228d0,     -19.704908d0, &
            -19.516651d0,     -19.344023d0,     -19.185917d0,     -19.041429d0,     -18.909777d0,     -18.790222d0, &
            -18.681981d0 &
        ], [NCOOLT,NCOOLNE] )

   ! Ca II and Mg I. The single ground level is the whole ground term, so
   ! these two tables carry every channel of the ion; they replace the
   ! resonance-line fits in metal_line_cooling_coefficient (through
   ! cool_CaII_ne_func / cool_MgI_ne_func). Ca II adds the 4s-3d excitation
   ! at 1.7 eV (the metastable 3d 2D levels, [Ca II] 7291/7324 A and the
   ! 3d-4p 8498/8542/8662 A triplet they feed); Mg I adds the 3s3p 3P term
   ! at 2.7 eV (the 4571 A line, the metastable 3P0 and 3P2 levels, and the
   ! 3s4s 3S 5167-5184 A triplet above them). Critical densities at
   ! 5e3-1e4 K: Ca II 3d 2D 4e5-1e7 cm^-3; Mg I 3P2 6e4-9e4 cm^-3, 3P1
   ! 4e9-5e9 cm^-3; Mg I 3P0 has no radiative decay. Ratio to the resonance
   ! fits at n_e = 1e9 cm^-3: Ca II 4.0 (3e3 K), 3.1 (5e3 K), 1.7 (1e4 K),
   ! 1.2 (2e4 K); Mg I 481 (3e3 K), 34 (5e3 K), 5.8 (1e4 K), 2.3 (2e4 K).
   ! Validity: optically thin; electron collisions only (neutral-hydrogen
   ! impact is not in CHIANTI, so a base with n_HI >> n_e has more
   ! collisional coupling than these tables); CHIANTI v11.0.2 carries no
   ! collision strength among the Mg I 3P_J levels, and with a mixing
   ! Upsilon of 0.1 / 1 / 10 among all three the Mg I coefficient would be
   ! up to 1.1 / 1.7 / 2.5 times the table at n_e = 1e5-1e9 cm^-3 over
   ! 3e3-3e4 K (under 2 per cent at n_e >= 1e11, 3-17 per cent at
   ! n_e <= 1e4). Below 1e3 cm^-3 the Mg I column is flat to 1.1 per cent
   ! (the 3P0 population does not depend on n_e
   ! there), at 0.90 to 0.95 of the coronal sum of its excitation rates over
   ! 3e3-3e4 K, because an excitation to 3P0 is not radiated.
   ! Atomic data as cited in the CHIANTI v11.0.2 files: Ca II, Melendez,
   ! Bautista & Badnell (2007, A&A 469, 1203); Mg I collision strengths,
   ! Barklem et al. (2017, A&A 606, A11), R-matrix in LS coupling split into
   ! levels by CHIANTI, which flags them as meant for ionization-balance
   ! models rather than line diagnostics.
   ! Tables: CHIANTI v11.0.2, the statistical equilibrium of cooling_data/
   ! multilevel_statistical_equilibrium.py (built by metal_cooling_density_
   ! resolved.py), collision strengths descaled with the natural cubic
   ! spline of CHIANTI's DESCALE_SCUPS.PRO. With the descaling made equal,
   ! ChiantiPy 0.15.2 populate gives the same Ca II and Mg I values to
   ! 4e-5 over 1e3-1e5 K and 1-1e14 cm^-3 (MEASURED). The bilinear log-log
   ! interpolation of
   ! remainder_table_value departs from that solve by at most 5 per cent
   ! (Ca II) and 6 per cent (Mg I) at the cell midpoints over 2e3-3e4 K and
   ! 1e2-1e12 cm^-3, median under 1 per cent.

   ! CaII
   real*8, parameter :: cool_logLrem_CaII(NCOOLT,NCOOLNE) = reshape( [ &
            -25.881431d0,     -24.978962d0,     -24.177843d0,     -23.467088d0,     -22.836899d0,     -22.278522d0, &
            -21.784119d0,     -21.346636d0,     -20.959676d0,     -20.617384d0,     -20.314348d0,     -20.045543d0, &
            -19.806346d0,     -19.592552d0,     -19.400415d0,     -19.226709d0,     -19.068776d0,     -18.924515d0, &
            -18.792327d0,     -18.671009d0,     -18.559634d0,     -18.457407d0,     -18.363584d0,     -18.277473d0, &
            -18.198430d0,     -18.125856d0,     -18.059200d0,     -17.997971d0,     -17.941750d0,     -17.890192d0, &
            -17.843023d0,     -17.800033d0,     -17.761050d0,     -17.725901d0,     -17.694396d0,     -17.666350d0, &
            -17.641572d0,     -17.619863d0,     -17.601020d0,     -17.584837d0,     -17.571112d0,     -25.881431d0, &
            -24.978962d0,     -24.177843d0,     -23.467088d0,     -22.836899d0,     -22.278522d0,     -21.784120d0, &
            -21.346637d0,     -20.959677d0,     -20.617384d0,     -20.314348d0,     -20.045544d0,     -19.806347d0, &
            -19.592552d0,     -19.400415d0,     -19.226709d0,     -19.068776d0,     -18.924515d0,     -18.792327d0, &
            -18.671010d0,     -18.559634d0,     -18.457407d0,     -18.363584d0,     -18.277473d0,     -18.198430d0, &
            -18.125856d0,     -18.059200d0,     -17.997971d0,     -17.941750d0,     -17.890192d0,     -17.843024d0, &
            -17.800033d0,     -17.761051d0,     -17.725901d0,     -17.694396d0,     -17.666350d0,     -17.641572d0, &
            -17.619863d0,     -17.601020d0,     -17.584837d0,     -17.571112d0,     -25.881432d0,     -24.978963d0, &
            -24.177844d0,     -23.467089d0,     -22.836900d0,     -22.278523d0,     -21.784120d0,     -21.346637d0, &
            -20.959677d0,     -20.617385d0,     -20.314348d0,     -20.045544d0,     -19.806347d0,     -19.592553d0, &
            -19.400415d0,     -19.226710d0,     -19.068776d0,     -18.924515d0,     -18.792327d0,     -18.671010d0, &
            -18.559634d0,     -18.457407d0,     -18.363585d0,     -18.277473d0,     -18.198431d0,     -18.125856d0, &
            -18.059200d0,     -17.997971d0,     -17.941750d0,     -17.890192d0,     -17.843024d0,     -17.800033d0, &
            -17.761051d0,     -17.725901d0,     -17.694396d0,     -17.666350d0,     -17.641572d0,     -17.619863d0, &
            -17.601020d0,     -17.584837d0,     -17.571112d0,     -25.881436d0,     -24.978966d0,     -24.177847d0, &
            -23.467092d0,     -22.836902d0,     -22.278525d0,     -21.784123d0,     -21.346640d0,     -20.959679d0, &
            -20.617387d0,     -20.314350d0,     -20.045546d0,     -19.806349d0,     -19.592554d0,     -19.400416d0, &
            -19.226711d0,     -19.068777d0,     -18.924516d0,     -18.792328d0,     -18.671011d0,     -18.559635d0, &
            -18.457408d0,     -18.363585d0,     -18.277474d0,     -18.198431d0,     -18.125857d0,     -18.059201d0, &
            -17.997972d0,     -17.941751d0,     -17.890193d0,     -17.843024d0,     -17.800034d0,     -17.761051d0, &
            -17.725902d0,     -17.694397d0,     -17.666351d0,     -17.641573d0,     -17.619864d0,     -17.601021d0, &
            -17.584837d0,     -17.571113d0,     -25.881447d0,     -24.978977d0,     -24.177857d0,     -23.467101d0, &
            -22.836911d0,     -22.278534d0,     -21.784130d0,     -21.346647d0,     -20.959686d0,     -20.617393d0, &
            -20.314356d0,     -20.045551d0,     -19.806353d0,     -19.592558d0,     -19.400420d0,     -19.226714d0, &
            -19.068780d0,     -18.924519d0,     -18.792331d0,     -18.671013d0,     -18.559637d0,     -18.457410d0, &
            -18.363587d0,     -18.277476d0,     -18.198434d0,     -18.125859d0,     -18.059203d0,     -17.997974d0, &
            -17.941753d0,     -17.890195d0,     -17.843027d0,     -17.800036d0,     -17.761054d0,     -17.725904d0, &
            -17.694399d0,     -17.666353d0,     -17.641575d0,     -17.619866d0,     -17.601023d0,     -17.584840d0, &
            -17.571115d0,     -25.881482d0,     -24.979010d0,     -24.177888d0,     -23.467130d0,     -22.836938d0, &
            -22.278559d0,     -21.784154d0,     -21.346669d0,     -20.959707d0,     -20.617412d0,     -20.314374d0, &
            -20.045567d0,     -19.806368d0,     -19.592572d0,     -19.400433d0,     -19.226725d0,     -19.068790d0, &
            -18.924528d0,     -18.792339d0,     -18.671021d0,     -18.559645d0,     -18.457417d0,     -18.363594d0, &
            -18.277483d0,     -18.198441d0,     -18.125866d0,     -18.059210d0,     -17.997981d0,     -17.941760d0, &
            -17.890202d0,     -17.843034d0,     -17.800043d0,     -17.761061d0,     -17.725911d0,     -17.694406d0, &
            -17.666360d0,     -17.641582d0,     -17.619873d0,     -17.601029d0,     -17.584846d0,     -17.571121d0, &
            -25.881594d0,     -24.979114d0,     -24.177986d0,     -23.467222d0,     -22.837024d0,     -22.278640d0, &
            -21.784230d0,     -21.346739d0,     -20.959772d0,     -20.617473d0,     -20.314430d0,     -20.045620d0, &
            -19.806416d0,     -19.592616d0,     -19.400472d0,     -19.226760d0,     -19.068822d0,     -18.924556d0, &
            -18.792365d0,     -18.671045d0,     -18.559668d0,     -18.457440d0,     -18.363616d0,     -18.277505d0, &
            -18.198463d0,     -18.125889d0,     -18.059233d0,     -17.998004d0,     -17.941783d0,     -17.890225d0, &
            -17.843057d0,     -17.800066d0,     -17.761083d0,     -17.725933d0,     -17.694428d0,     -17.666382d0, &
            -17.641604d0,     -17.619894d0,     -17.601050d0,     -17.584867d0,     -17.571142d0,     -25.881946d0, &
            -24.979445d0,     -24.178296d0,     -23.467513d0,     -22.837297d0,     -22.278895d0,     -21.784468d0, &
            -21.346962d0,     -20.959980d0,     -20.617666d0,     -20.314609d0,     -20.045784d0,     -19.806567d0, &
            -19.592753d0,     -19.400595d0,     -19.226871d0,     -19.068921d0,     -18.924646d0,     -18.792447d0, &
            -18.671121d0,     -18.559740d0,     -18.457510d0,     -18.363686d0,     -18.277575d0,     -18.198533d0, &
            -18.125959d0,     -18.059304d0,     -17.998076d0,     -17.941855d0,     -17.890297d0,     -17.843128d0, &
            -17.800137d0,     -17.761154d0,     -17.726004d0,     -17.694498d0,     -17.666451d0,     -17.641672d0, &
            -17.619961d0,     -17.601117d0,     -17.584932d0,     -17.571206d0,     -25.883057d0,     -24.980488d0, &
            -24.179275d0,     -23.468431d0,     -22.838157d0,     -22.279700d0,     -21.785222d0,     -21.347666d0, &
            -20.960636d0,     -20.618277d0,     -20.315174d0,     -20.046305d0,     -19.807044d0,     -19.593186d0, &
            -19.400986d0,     -19.227222d0,     -19.069236d0,     -18.924930d0,     -18.792706d0,     -18.671361d0, &
            -18.559968d0,     -18.457732d0,     -18.363905d0,     -18.277794d0,     -18.198753d0,     -18.126182d0, &
            -18.059529d0,     -17.998302d0,     -17.942082d0,     -17.890523d0,     -17.843354d0,     -17.800362d0, &
            -17.761377d0,     -17.726225d0,     -17.694717d0,     -17.666668d0,     -17.641886d0,     -17.620174d0, &
            -17.601326d0,     -17.585138d0,     -17.571409d0,     -25.886554d0,     -24.983771d0,     -24.182356d0, &
            -23.471321d0,     -22.840866d0,     -22.282238d0,     -21.787596d0,     -21.349885d0,     -20.962706d0, &
            -20.620201d0,     -20.316957d0,     -20.047947d0,     -19.808546d0,     -19.594551d0,     -19.402218d0, &
            -19.228328d0,     -19.070227d0,     -18.925823d0,     -18.793520d0,     -18.672118d0,     -18.560687d0, &
            -18.458428d0,     -18.364592d0,     -18.278481d0,     -18.199445d0,     -18.126881d0,     -18.060233d0, &
            -17.999011d0,     -17.942792d0,     -17.891234d0,     -17.844062d0,     -17.801066d0,     -17.762076d0, &
            -17.726918d0,     -17.695404d0,     -17.667348d0,     -17.642559d0,     -17.620838d0,     -17.601982d0, &
            -17.585784d0,     -17.572044d0,     -25.897435d0,     -24.994000d0,     -24.191965d0,     -23.480342d0, &
            -22.849329d0,     -22.290170d0,     -21.795025d0,     -21.356831d0,     -20.969188d0,     -20.626231d0, &
            -20.322544d0,     -20.053095d0,     -19.813257d0,     -19.598828d0,     -19.406074d0,     -19.231785d0, &
            -19.073323d0,     -18.928609d0,     -18.796057d0,     -18.674471d0,     -18.562917d0,     -18.460588d0, &
            -18.366721d0,     -18.280607d0,     -18.201583d0,     -18.129036d0,     -18.062405d0,     -18.001193d0, &
            -17.944980d0,     -17.893419d0,     -17.846240d0,     -17.803232d0,     -17.764226d0,     -17.729050d0, &
            -17.697517d0,     -17.669440d0,     -17.644628d0,     -17.622883d0,     -17.604000d0,     -17.587773d0, &
            -17.574002d0,     -25.930164d0,     -25.024855d0,     -24.221032d0,     -23.507703d0,     -22.875063d0, &
            -22.314351d0,     -21.817719d0,     -21.378095d0,     -20.989065d0,     -20.644751d0,     -20.339717d0, &
            -20.068923d0,     -19.827735d0,     -19.611959d0,     -19.417885d0,     -19.242344d0,     -19.082745d0, &
            -18.937052d0,     -18.803711d0,     -18.681536d0,     -18.569581d0,     -18.467011d0,     -18.373026d0, &
            -18.286877d0,     -18.207867d0,     -18.135352d0,     -18.068753d0,     -18.007561d0,     -17.951349d0, &
            -17.899774d0,     -17.852566d0,     -17.809517d0,     -17.770465d0,     -17.735236d0,     -17.703648d0, &
            -17.675513d0,     -17.650640d0,     -17.628828d0,     -17.609873d0,     -17.593569d0,     -17.579714d0, &
            -26.019919d0,     -25.110108d0,     -24.301927d0,     -23.584382d0,     -22.947663d0,     -22.383003d0, &
            -21.882534d0,     -21.439159d0,     -21.046420d0,     -20.698388d0,     -20.389576d0,     -20.114904d0, &
            -19.869729d0,     -19.649898d0,     -19.451795d0,     -19.272397d0,     -19.109275d0,     -18.960541d0, &
            -18.824733d0,     -18.700686d0,     -18.587412d0,     -18.483986d0,     -18.389499d0,     -18.303091d0, &
            -18.223965d0,     -18.151403d0,     -18.084772d0,     -18.023535d0,     -17.967254d0,     -17.915584d0, &
            -17.868261d0,     -17.825087d0,     -17.785902d0,     -17.750542d0,     -17.718825d0,     -17.690566d0, &
            -17.665566d0,     -17.643624d0,     -17.624530d0,     -17.608076d0,     -17.594059d0,     -26.221437d0, &
            -25.304246d0,     -24.488722d0,     -23.763876d0,     -23.119896d0,     -22.547991d0,     -22.040236d0, &
            -21.589429d0,     -21.188951d0,     -20.832669d0,     -20.514891d0,     -20.230401d0,     -19.974572d0, &
            -19.743473d0,     -19.533936d0,     -19.343516d0,     -19.170366d0,     -19.013036d0,     -18.870283d0, &
            -18.740931d0,     -18.623806d0,     -18.517706d0,     -18.421425d0,     -18.333832d0,     -18.253907d0, &
            -18.180767d0,     -18.113672d0,     -18.052024d0,     -17.995356d0,     -17.943317d0,     -17.895649d0, &
            -17.852160d0,     -17.812697d0,     -17.777100d0,     -17.745189d0,     -17.716771d0,     -17.691641d0, &
            -17.669583d0,     -17.650381d0,     -17.633819d0,     -17.619685d0,     -26.558692d0,     -25.634638d0, &
            -24.812051d0,     -24.079925d0,     -23.428398d0,     -22.848548d0,     -22.332193d0,     -21.871667d0, &
            -21.459673d0,     -21.089244d0,     -20.753932d0,     -20.448192d0,     -20.167830d0,     -19.910221d0, &
            -19.674079d0,     -19.458887d0,     -19.264233d0,     -19.089369d0,     -18.933062d0,     -18.793681d0, &
            -18.669372d0,     -18.558231d0,     -18.458443d0,     -18.368387d0,     -18.286682d0,     -18.212198d0, &
            -18.144036d0,     -18.081503d0,     -18.024084d0,     -17.971405d0,     -17.923200d0,     -17.879273d0, &
            -17.839468d0,     -17.803623d0,     -17.771552d0,     -17.743052d0,     -17.717904d0,     -17.695877d0, &
            -17.676743d0,     -17.660275d0,     -17.646252d0,     -26.991287d0,     -26.063495d0,     -25.236943d0, &
            -24.500549d0,     -23.844248d0,     -23.258661d0,     -22.734675d0,     -22.263059d0,     -21.834338d0, &
            -21.439326d0,     -21.070446d0,     -20.723285d0,     -20.397190d0,     -20.094222d0,     -19.817178d0, &
            -19.567923d0,     -19.346664d0,     -19.152037d0,     -18.981594d0,     -18.832365d0,     -18.701300d0, &
            -18.585547d0,     -18.482589d0,     -18.390314d0,     -18.307009d0,     -18.231326d0,     -18.162230d0, &
            -18.098948d0,     -18.040918d0,     -17.987740d0,     -17.939136d0,     -17.894900d0,     -17.854870d0, &
            -17.818876d0,     -17.786728d0,     -17.758217d0,     -17.733109d0,     -17.711168d0,     -17.692157d0, &
            -17.675841d0,     -17.661994d0,     -27.467590d0,     -26.538311d0,     -25.710091d0,     -24.971627d0, &
            -24.312241d0,     -23.721119d0,     -23.186373d0,     -22.694538d0,     -22.231760d0,     -21.787466d0, &
            -21.358362d0,     -20.948604d0,     -20.565734d0,     -20.216333d0,     -19.903998d0,     -19.629242d0, &
            -19.390193d0,     -19.183459d0,     -19.004915d0,     -18.850317d0,     -18.715701d0,     -18.597582d0, &
            -18.493027d0,     -18.399647d0,     -18.315552d0,     -18.239282d0,     -18.169735d0,     -18.106097d0, &
            -18.047781d0,     -17.994374d0,     -17.945589d0,     -17.901217d0,     -17.861089d0,     -17.825034d0, &
            -17.792859d0,     -17.764349d0,     -17.739269d0,     -17.717377d0,     -17.698433d0,     -17.682200d0, &
            -17.668449d0,     -27.959805d0,     -27.029954d0,     -26.200867d0,     -25.460557d0,     -24.796466d0, &
            -24.193563d0,     -23.632843d0,     -23.093610d0,     -22.562127d0,     -22.039538d0,     -21.538237d0, &
            -21.071566d0,     -20.648041d0,     -20.271096d0,     -19.940552d0,     -19.653888d0,     -19.407096d0, &
            -19.195346d0,     -19.013561d0,     -18.856869d0,     -18.720892d0,     -18.601878d0,     -18.496723d0, &
            -18.402931d0,     -18.318542d0,     -18.242056d0,     -18.172343d0,     -18.108575d0,     -18.050154d0, &
            -17.996665d0,     -17.947816d0,     -17.903394d0,     -17.863233d0,     -17.827156d0,     -17.794972d0, &
            -17.766463d0,     -17.741394d0,     -17.719521d0,     -17.700603d0,     -17.684401d0,     -17.670687d0, &
            -28.457272d0,     -27.527029d0,     -26.696716d0,     -25.952238d0,     -25.275460d0,     -24.640637d0, &
            -24.017703d0,     -23.389696d0,     -22.766131d0,     -22.169935d0,     -21.618919d0,     -21.121373d0, &
            -20.679179d0,     -20.290908d0,     -19.953395d0,     -19.662378d0,     -19.412838d0,     -19.199345d0, &
            -19.016449d0,     -18.859045d0,     -18.722608d0,     -18.603294d0,     -18.497939d0,     -18.404008d0, &
            -18.319522d0,     -18.242963d0,     -18.173196d0,     -18.109384d0,     -18.050929d0,     -17.997412d0, &
            -17.948541d0,     -17.904104d0,     -17.863931d0,     -17.827847d0,     -17.795660d0,     -17.767153d0, &
            -17.742087d0,     -17.720221d0,     -17.701311d0,     -17.685119d0,     -17.671418d0,     -28.956333d0, &
            -28.025313d0,     -27.191643d0,     -26.434900d0,     -25.722527d0,     -25.012233d0,     -24.282040d0, &
            -23.552587d0,     -22.858731d0,     -22.221645d0,     -21.648338d0,     -21.138657d0,     -20.689677d0, &
            -20.297475d0,     -19.957608d0,     -19.665145d0,     -19.414703d0,     -19.200640d0,     -19.017384d0, &
            -18.859749d0,     -18.723164d0,     -18.603752d0,     -18.498332d0,     -18.404357d0,     -18.319839d0, &
            -18.243258d0,     -18.173472d0,     -18.109646d0,     -18.051181d0,     -17.997655d0,     -17.948777d0, &
            -17.904335d0,     -17.864158d0,     -17.828073d0,     -17.795885d0,     -17.767377d0,     -17.742313d0, &
            -17.720449d0,     -17.701542d0,     -17.685354d0,     -17.671656d0,     -29.455610d0,     -28.522290d0, &
            -27.678409d0,     -26.886526d0,     -26.089869d0,     -25.255791d0,     -24.416327d0,     -23.620374d0, &
            -22.892833d0,     -22.239473d0,     -21.658143d0,     -21.144322d0,     -20.693092d0,     -20.299606d0, &
            -19.958977d0,     -19.666048d0,     -19.415316d0,     -19.201071d0,     -19.017698d0,     -18.859989d0, &
            -18.723356d0,     -18.603913d0,     -18.498472d0,     -18.404483d0,     -18.319954d0,     -18.243365d0, &
            -18.173574d0,     -18.109744d0,     -18.051275d0,     -17.997746d0,     -17.948867d0,     -17.904423d0, &
            -17.864245d0,     -17.828159d0,     -17.795970d0,     -17.767463d0,     -17.742400d0,     -17.720536d0, &
            -17.701630d0,     -17.685443d0,     -17.671746d0,     -29.954038d0,     -29.013605d0,     -28.139755d0, &
            -27.262183d0,     -26.327248d0,     -25.373981d0,     -24.469344d0,     -23.644345d0,     -22.904310d0, &
            -22.245366d0,     -21.661383d0,     -21.146214d0,     -20.694257d0,     -20.300356d0,     -19.959480d0, &
            -19.666399d0,     -19.415570d0,     -19.201263d0,     -19.017851d0,     -18.860116d0,     -18.723466d0, &
            -18.604012d0,     -18.498563d0,     -18.404569d0,     -18.320037d0,     -18.243445d0,     -18.173652d0, &
            -18.109820d0,     -18.051349d0,     -17.997820d0,     -17.948940d0,     -17.904496d0,     -17.864318d0, &
            -17.828232d0,     -17.796044d0,     -17.767537d0,     -17.742474d0,     -17.720611d0,     -17.701705d0, &
            -17.685519d0,     -17.671822d0,     -30.449334d0,     -29.487489d0,     -28.536499d0,     -27.510337d0, &
            -26.441128d0,     -25.419569d0,     -24.487894d0,     -23.652535d0,     -22.908312d0,     -22.247540d0, &
            -21.662690d0,     -21.147077d0,     -20.694873d0,     -20.300826d0,     -19.959856d0,     -19.666712d0, &
            -19.415840d0,     -19.201502d0,     -19.018068d0,     -18.860318d0,     -18.723656d0,     -18.604195d0, &
            -18.498740d0,     -18.404742d0,     -18.320207d0,     -18.243614d0,     -18.173819d0,     -18.109986d0, &
            -18.051515d0,     -17.997986d0,     -17.949105d0,     -17.904662d0,     -17.864485d0,     -17.828399d0, &
            -17.796213d0,     -17.767707d0,     -17.742646d0,     -17.720784d0,     -17.701879d0,     -17.685693d0, &
            -17.671996d0,     -30.934924d0,     -29.914308d0,     -28.814814d0,     -27.632744d0,     -26.485530d0, &
            -25.436160d0,     -24.495001d0,     -23.656179d0,     -22.910558d0,     -22.249153d0,     -21.663980d0, &
            -21.148178d0,     -20.695848d0,     -20.301706d0,     -19.960659d0,     -19.667451d0,     -19.416524d0, &
            -19.202141d0,     -19.018670d0,     -18.860892d0,     -18.724209d0,     -18.604732d0,     -18.499266d0, &
            -18.405259d0,     -18.320718d0,     -18.244120d0,     -18.174323d0,     -18.110488d0,     -18.052016d0, &
            -17.998487d0,     -17.949608d0,     -17.905166d0,     -17.864992d0,     -17.828910d0,     -17.796727d0, &
            -17.768226d0,     -17.743168d0,     -17.721309d0,     -17.702407d0,     -17.686222d0,     -17.672525d0, &
            -31.392920d0,     -30.243135d0,     -28.963309d0,     -27.683681d0,     -26.504268d0,     -25.445089d0, &
            -24.500656d0,     -23.660558d0,     -22.914332d0,     -22.252572d0,     -21.667145d0,     -21.151132d0, &
            -20.698611d0,     -20.304288d0,     -19.963069d0,     -19.669699d0,     -19.418627d0,     -19.204118d0, &
            -19.020544d0,     -18.862683d0,     -18.725936d0,     -18.606412d0,     -18.500911d0,     -18.406880d0, &
            -18.322321d0,     -18.245710d0,     -18.175904d0,     -18.112064d0,     -18.053589d0,     -18.000061d0, &
            -17.951185d0,     -17.906750d0,     -17.866584d0,     -17.830514d0,     -17.798343d0,     -17.769855d0, &
            -17.744809d0,     -17.722961d0,     -17.704065d0,     -17.687883d0,     -17.674185d0,     -31.786493d0, &
            -30.445658d0,     -29.035245d0,     -27.713128d0,     -26.521928d0,     -25.458984d0,     -24.512965d0, &
            -23.671949d0,     -22.925043d0,     -22.262699d0,     -21.676727d0,     -21.160187d0,     -20.707141d0, &
            -20.312296d0,     -19.970565d0,     -19.676710d0,     -19.425195d0,     -19.210301d0,     -19.026406d0, &
            -18.868290d0,     -18.731348d0,     -18.611677d0,     -18.506069d0,     -18.411960d0,     -18.327345d0, &
            -18.250693d0,     -18.180859d0,     -18.117001d0,     -18.058518d0,     -18.004991d0,     -17.956126d0, &
            -17.911710d0,     -17.871572d0,     -17.835535d0,     -17.803404d0,     -17.774956d0,     -17.749949d0, &
            -17.728132d0,     -17.709258d0,     -17.693086d0,     -17.679385d0,     -32.083778d0,     -30.571511d0, &
            -29.097651d0,     -27.758732d0,     -26.562010d0,     -25.496270d0,     -24.548234d0,     -23.705484d0, &
            -22.956976d0,     -22.293100d0,     -21.705628d0,     -21.187590d0,     -20.733035d0,     -20.336672d0, &
            -19.993442d0,     -19.698153d0,     -19.445326d0,     -19.229285d0,     -19.044435d0,     -18.885555d0, &
            -18.748022d0,     -18.627907d0,     -18.521970d0,     -18.427619d0,     -18.342826d0,     -18.266045d0, &
            -18.196116d0,     -18.132195d0,     -18.073679d0,     -18.020148d0,     -17.971309d0,     -17.926948d0, &
            -17.886891d0,     -17.850957d0,     -17.818944d0,     -17.790621d0,     -17.765731d0,     -17.744012d0, &
            -17.725207d0,     -17.709066d0,     -17.695355d0,     -32.328429d0,     -30.717521d0,     -29.216296d0, &
            -27.868048d0,     -26.665874d0,     -25.595691d0,     -24.643589d0,     -23.796977d0,     -23.044737d0, &
            -22.377199d0,     -21.786079d0,     -21.264349d0,     -20.806022d0,     -20.405815d0,     -20.058739d0, &
            -19.759729d0,     -19.503452d0,     -19.284364d0,     -19.096944d0,     -18.935986d0,     -18.796827d0, &
            -18.675467d0,     -18.568585d0,     -18.473513d0,     -18.388164d0,     -18.310947d0,     -18.240679d0, &
            -18.176506d0,     -18.117825d0,     -18.064220d0,     -18.015405d0,     -17.971162d0,     -17.931311d0, &
            -17.895656d0,     -17.863978d0,     -17.836011d0,     -17.811461d0,     -17.790030d0,     -17.771430d0, &
            -17.755393d0,     -17.741670d0,     -32.624426d0,     -30.969114d0,     -29.453791d0,     -28.097564d0, &
            -26.888693d0,     -25.812181d0,     -24.853915d0,     -24.001237d0,     -23.242975d0,     -22.569389d0, &
            -21.972098d0,     -21.443976d0,     -20.978960d0,     -20.571767d0,     -20.217520d0,     -19.911385d0, &
            -19.648332d0,     -19.423099d0,     -19.230360d0,     -19.064970d0,     -18.922213d0,     -18.797955d0, &
            -18.688715d0,     -18.591671d0,     -18.504619d0,     -18.425902d0,     -18.354316d0,     -18.289026d0, &
            -18.229466d0,     -18.175264d0,     -18.126164d0,     -18.081957d0,     -18.042441d0,     -18.007376d0, &
            -17.976483d0,     -17.949402d0,     -17.925730d0,     -17.905076d0,     -17.887077d0,     -17.871411d0, &
            -17.857796d0 &
        ], [NCOOLT,NCOOLNE] )

   ! MgI
   real*8, parameter :: cool_logLrem_MgI(NCOOLT,NCOOLNE) = reshape( [ &
            -31.166016d0,     -29.685005d0,     -28.364513d0,     -27.187299d0,     -26.138216d0,     -25.204051d0, &
            -24.373408d0,     -23.636377d0,     -22.983435d0,     -22.405761d0,     -21.895282d0,     -21.444433d0, &
            -21.045901d0,     -20.692677d0,     -20.379346d0,     -20.101410d0,     -19.854882d0,     -19.636249d0, &
            -19.442294d0,     -19.269738d0,     -19.115484d0,     -18.976765d0,     -18.851385d0,     -18.737529d0, &
            -18.633680d0,     -18.538597d0,     -18.451287d0,     -18.370937d0,     -18.296739d0,     -18.227995d0, &
            -18.164162d0,     -18.104862d0,     -18.049857d0,     -17.998897d0,     -17.951319d0,     -17.906636d0, &
            -17.864676d0,     -17.825488d0,     -17.789003d0,     -17.754709d0,     -17.722453d0,     -31.166029d0, &
            -29.685018d0,     -28.364526d0,     -27.187311d0,     -26.138228d0,     -25.204064d0,     -24.373420d0, &
            -23.636389d0,     -22.983446d0,     -22.405773d0,     -21.895293d0,     -21.444443d0,     -21.045911d0, &
            -20.692686d0,     -20.379355d0,     -20.101418d0,     -19.854890d0,     -19.636256d0,     -19.442301d0, &
            -19.269744d0,     -19.115490d0,     -18.976771d0,     -18.851391d0,     -18.737534d0,     -18.633685d0, &
            -18.538602d0,     -18.451292d0,     -18.370942d0,     -18.296745d0,     -18.228000d0,     -18.164168d0, &
            -18.104867d0,     -18.049862d0,     -17.998902d0,     -17.951324d0,     -17.906641d0,     -17.864681d0, &
            -17.825493d0,     -17.789007d0,     -17.754713d0,     -17.722457d0,     -31.166070d0,     -29.685058d0, &
            -28.364565d0,     -27.187351d0,     -26.138267d0,     -25.204102d0,     -24.373458d0,     -23.636426d0, &
            -22.983483d0,     -22.405808d0,     -21.895327d0,     -21.444476d0,     -21.045941d0,     -20.692715d0, &
            -20.379382d0,     -20.101444d0,     -19.854915d0,     -19.636279d0,     -19.442322d0,     -19.269764d0, &
            -19.115509d0,     -18.976789d0,     -18.851408d0,     -18.737551d0,     -18.633702d0,     -18.538618d0, &
            -18.451308d0,     -18.370958d0,     -18.296761d0,     -18.228017d0,     -18.164185d0,     -18.104884d0, &
            -18.049879d0,     -17.998919d0,     -17.951340d0,     -17.906657d0,     -17.864696d0,     -17.825507d0, &
            -17.789021d0,     -17.754726d0,     -17.722469d0,     -31.166198d0,     -29.685185d0,     -28.364691d0, &
            -27.187475d0,     -26.138391d0,     -25.204225d0,     -24.373579d0,     -23.636545d0,     -22.983599d0, &
            -22.405920d0,     -21.895434d0,     -21.444578d0,     -21.046039d0,     -20.692807d0,     -20.379469d0, &
            -20.101526d0,     -19.854992d0,     -19.636352d0,     -19.442390d0,     -19.269828d0,     -19.115569d0, &
            -18.976846d0,     -18.851462d0,     -18.737603d0,     -18.633753d0,     -18.538669d0,     -18.451359d0, &
            -18.371009d0,     -18.296813d0,     -18.228069d0,     -18.164238d0,     -18.104937d0,     -18.049933d0, &
            -17.998971d0,     -17.951392d0,     -17.906708d0,     -17.864745d0,     -17.825554d0,     -17.789065d0, &
            -17.754768d0,     -17.722508d0,     -31.166603d0,     -29.685585d0,     -28.365087d0,     -27.187868d0, &
            -26.138781d0,     -25.204612d0,     -24.373962d0,     -23.636920d0,     -22.983964d0,     -22.406274d0, &
            -21.895774d0,     -21.444902d0,     -21.046346d0,     -20.693098d0,     -20.379744d0,     -20.101786d0, &
            -19.855236d0,     -19.636581d0,     -19.442605d0,     -19.270030d0,     -19.115759d0,     -18.977025d0, &
            -18.851634d0,     -18.737769d0,     -18.633915d0,     -18.538830d0,     -18.451520d0,     -18.371171d0, &
            -18.296976d0,     -18.228235d0,     -18.164405d0,     -18.105106d0,     -18.050101d0,     -17.999138d0, &
            -17.951556d0,     -17.906867d0,     -17.864899d0,     -17.825702d0,     -17.789206d0,     -17.754900d0, &
            -17.722631d0,     -31.167878d0,     -29.686844d0,     -28.366334d0,     -27.189105d0,     -26.140009d0, &
            -25.205829d0,     -24.375165d0,     -23.638102d0,     -22.985116d0,     -22.407388d0,     -21.896844d0, &
            -21.445923d0,     -21.047315d0,     -20.694016d0,     -20.380611d0,     -20.102603d0,     -19.856005d0, &
            -19.637303d0,     -19.443282d0,     -19.270665d0,     -19.116356d0,     -18.977591d0,     -18.852174d0, &
            -18.738292d0,     -18.634426d0,     -18.539335d0,     -18.452024d0,     -18.371679d0,     -18.297491d0, &
            -18.228756d0,     -18.164931d0,     -18.105635d0,     -18.050629d0,     -17.999662d0,     -17.952071d0, &
            -17.907369d0,     -17.865384d0,     -17.826166d0,     -17.789646d0,     -17.755314d0,     -17.723018d0, &
            -31.171853d0,     -29.690770d0,     -28.370222d0,     -27.192963d0,     -26.143839d0,     -25.209628d0, &
            -24.378919d0,     -23.641789d0,     -22.988712d0,     -22.410869d0,     -21.900188d0,     -21.449115d0, &
            -21.050346d0,     -20.696887d0,     -20.383325d0,     -20.105162d0,     -19.858412d0,     -19.639563d0, &
            -19.445400d0,     -19.272651d0,     -19.118224d0,     -18.979359d0,     -18.853862d0,     -18.739920d0, &
            -18.636017d0,     -18.540907d0,     -18.453595d0,     -18.373260d0,     -18.299089d0,     -18.230373d0, &
            -18.166564d0,     -18.107276d0,     -18.052268d0,     -18.001288d0,     -17.953670d0,     -17.908927d0, &
            -17.866888d0,     -17.827605d0,     -17.791012d0,     -17.756601d0,     -17.724221d0,     -31.183876d0, &
            -29.702656d0,     -28.382000d0,     -27.204655d0,     -26.155452d0,     -25.221150d0,     -24.390313d0, &
            -23.652989d0,     -22.999645d0,     -22.421465d0,     -21.910383d0,     -21.458859d0,     -21.059612d0, &
            -20.705675d0,     -20.391639d0,     -20.113006d0,     -19.865792d0,     -19.646488d0,     -19.451886d0, &
            -19.278723d0,     -19.123924d0,     -18.984741d0,     -18.858988d0,     -18.744854d0,     -18.640822d0, &
            -18.545643d0,     -18.458313d0,     -18.377999d0,     -18.303870d0,     -18.235203d0,     -18.171435d0, &
            -18.112168d0,     -18.057154d0,     -18.006132d0,     -17.958433d0,     -17.913571d0,     -17.871374d0, &
            -17.831903d0,     -17.795094d0,     -17.760447d0,     -17.727818d0,     -31.217207d0,     -29.735670d0, &
            -28.414768d0,     -27.237225d0,     -26.187841d0,     -25.253328d0,     -24.422185d0,     -23.684389d0, &
            -23.030387d0,     -22.451363d0,     -21.939265d0,     -21.486580d0,     -21.086079d0,     -20.730859d0, &
            -20.415524d0,     -20.135575d0,     -19.887033d0,     -19.666400d0,     -19.470487d0,     -19.296072d0, &
            -19.140124d0,     -18.999938d0,     -18.873353d0,     -18.758572d0,     -18.654078d0,     -18.558608d0, &
            -18.471135d0,     -18.390795d0,     -18.316713d0,     -18.248122d0,     -18.184421d0,     -18.125182d0, &
            -18.070134d0,     -18.018999d0,     -17.971093d0,     -17.925924d0,     -17.883327d0,     -17.843372d0, &
            -17.806011d0,     -17.770755d0,     -17.737478d0,     -31.291307d0,     -29.809406d0,     -28.488229d0, &
            -27.310471d0,     -26.260884d0,     -25.326112d0,     -24.494553d0,     -23.756064d0,     -23.101046d0, &
            -22.520661d0,     -22.006846d0,     -21.552098d0,     -21.149239d0,     -20.791450d0,     -20.473338d0, &
            -20.190390d0,     -19.938641d0,     -19.714633d0,     -19.515252d0,     -19.337403d0,     -19.178202d0, &
            -19.035085d0,     -18.905983d0,     -18.789147d0,     -18.683067d0,     -18.586455d0,     -18.498228d0, &
            -18.417446d0,     -18.343141d0,     -18.274453d0,     -18.210702d0,     -18.151395d0,     -18.096211d0, &
            -18.044834d0,     -17.996537d0,     -17.950807d0,     -17.907480d0,     -17.866639d0,     -17.828254d0, &
            -17.791855d0,     -17.757341d0,     -31.399040d0,     -29.917444d0,     -28.596545d0,     -27.419040d0, &
            -26.369674d0,     -25.435070d0,     -24.603589d0,     -23.865032d0,     -23.209740d0,     -22.628802d0, &
            -22.114051d0,     -21.657861d0,     -21.252930d0,     -20.892347d0,     -20.570602d0,     -20.283099d0, &
            -20.025878d0,     -19.795599d0,     -19.589384d0,     -19.404497d0,     -19.238459d0,     -19.089073d0, &
            -18.954516d0,     -18.833156d0,     -18.723488d0,     -18.624154d0,     -18.533955d0,     -18.451809d0, &
            -18.376591d0,     -18.307300d0,     -18.243137d0,     -18.183516d0,     -18.128043d0,     -18.076341d0, &
            -18.027612d0,     -17.981297d0,     -17.937212d0,     -17.895443d0,     -17.855971d0,     -17.818326d0, &
            -17.782431d0,     -31.488542d0,     -30.007972d0,     -28.687948d0,     -27.511201d0,     -26.462522d0, &
            -25.528597d0,     -24.697866d0,     -23.960207d0,     -23.305946d0,     -22.726091d0,     -22.212324d0, &
            -21.756777d0,     -21.351813d0,     -20.990131d0,     -20.665915d0,     -20.374372d0,     -20.111513d0, &
            -19.874203d0,     -19.659997d0,     -19.466750d0,     -19.292609d0,     -19.135890d0,     -18.995072d0, &
            -18.868619d0,     -18.754972d0,     -18.652634d0,     -18.560239d0,     -18.476535d0,     -18.400234d0, &
            -18.330197d0,     -18.265520d0,     -18.205535d0,     -18.149787d0,     -18.097845d0,     -18.048849d0, &
            -18.002194d0,     -17.957677d0,     -17.915374d0,     -17.875264d0,     -17.836868d0,     -17.800119d0, &
            -31.534204d0,     -30.054441d0,     -28.735104d0,     -27.558953d0,     -26.510820d0,     -25.577453d0, &
            -24.747375d0,     -24.010557d0,     -23.357338d0,     -22.778690d0,     -22.266200d0,     -21.811826d0, &
            -21.407650d0,     -21.046004d0,     -20.720794d0,     -20.427041d0,     -20.160725d0,     -19.918894d0, &
            -19.699479d0,     -19.500811d0,     -19.321501d0,     -19.160208d0,     -19.015584d0,     -18.886108d0, &
            -18.770149d0,     -18.666094d0,     -18.572455d0,     -18.487871d0,     -18.410956d0,     -18.340495d0, &
            -18.275528d0,     -18.215343d0,     -18.159452d0,     -18.107400d0,     -18.058294d0,     -18.011507d0, &
            -17.966820d0,     -17.924308d0,     -17.883944d0,     -17.845243d0,     -17.808139d0,     -31.551963d0, &
            -30.072556d0,     -28.753525d0,     -27.577638d0,     -26.529749d0,     -25.596632d0,     -24.766854d0, &
            -24.030427d0,     -23.377704d0,     -22.799642d0,     -22.287793d0,     -21.834038d0,     -21.430326d0, &
            -21.068814d0,     -20.743265d0,     -20.448608d0,     -20.180811d0,     -19.937014d0,     -19.715335d0, &
            -19.514331d0,     -19.332822d0,     -19.169612d0,     -19.023417d0,     -18.892711d0,     -18.775823d0, &
            -18.671084d0,     -18.576953d0,     -18.492024d0,     -18.414869d0,     -18.344245d0,     -18.279167d0, &
            -18.218907d0,     -18.162966d0,     -18.110878d0,     -18.061738d0,     -18.014910d0,     -17.970170d0, &
            -17.927592d0,     -17.887145d0,     -17.848343d0,     -17.811119d0,     -31.559002d0,     -30.079706d0, &
            -28.760770d0,     -27.584966d0,     -26.537155d0,     -25.604118d0,     -24.774433d0,     -24.038129d0, &
            -23.385560d0,     -22.807681d0,     -22.296030d0,     -21.842463d0,     -21.438880d0,     -21.077374d0, &
            -20.751657d0,     -20.456621d0,     -20.188233d0,     -19.943672d0,     -19.721126d0,     -19.519239d0, &
            -19.336911d0,     -19.172996d0,     -19.026228d0,     -18.895079d0,     -18.777859d0,     -18.672879d0, &
            -18.578577d0,     -18.493531d0,     -18.416299d0,     -18.345625d0,     -18.280518d0,     -18.220243d0, &
            -18.164298d0,     -18.112211d0,     -18.063073d0,     -18.016246d0,     -17.971505d0,     -17.928919d0, &
            -17.888459d0,     -17.849636d0,     -17.812382d0,     -31.564577d0,     -30.085265d0,     -28.766319d0, &
            -27.590509d0,     -26.542691d0,     -25.609646d0,     -24.779949d0,     -24.043622d0,     -23.391021d0, &
            -22.813096d0,     -22.301384d0,     -21.847732d0,     -21.444026d0,     -21.082343d0,     -20.756373d0, &
            -20.461001d0,     -20.192195d0,     -19.947158d0,     -19.724115d0,     -19.521752d0,     -19.339000d0, &
            -19.174733d0,     -19.027688d0,     -18.896331d0,     -18.778961d0,     -18.673879d0,     -18.579512d0, &
            -18.494430d0,     -18.417185d0,     -18.346515d0,     -18.281425d0,     -18.221178d0,     -18.165269d0, &
            -18.113223d0,     -18.064130d0,     -18.017349d0,     -17.972653d0,     -17.930109d0,     -17.889686d0, &
            -17.850893d0,     -17.813661d0,     -31.576614d0,     -30.097139d0,     -28.778062d0,     -27.602144d0, &
            -26.554228d0,     -25.621076d0,     -24.791239d0,     -24.054714d0,     -23.401844d0,     -22.823580d0, &
            -22.311455d0,     -21.857315d0,     -21.453052d0,     -21.090744d0,     -20.764074d0,     -20.467924d0, &
            -20.198281d0,     -19.952387d0,     -19.728520d0,     -19.525416d0,     -19.342041d0,     -19.177280d0, &
            -19.029863d0,     -18.898240d0,     -18.780690d0,     -18.675499d0,     -18.581080d0,     -18.495991d0, &
            -18.418776d0,     -18.348167d0,     -18.283163d0,     -18.223023d0,     -18.167239d0,     -18.115330d0, &
            -18.066385d0,     -18.019757d0,     -17.975215d0,     -17.932821d0,     -17.892538d0,     -17.853869d0, &
            -17.816740d0,     -31.611097d0,     -30.131123d0,     -28.811644d0,     -27.635392d0,     -26.587171d0, &
            -25.653690d0,     -24.823423d0,     -24.086284d0,     -23.432589d0,     -22.853275d0,     -22.339871d0, &
            -21.884223d0,     -21.478238d0,     -21.114024d0,     -20.785248d0,     -20.486802d0,     -20.214737d0, &
            -19.966414d0,     -19.740250d0,     -19.535117d0,     -19.350061d0,     -19.183987d0,     -19.035593d0, &
            -18.903276d0,     -18.785264d0,     -18.679798d0,     -18.585253d0,     -18.500160d0,     -18.423040d0, &
            -18.352607d0,     -18.287849d0,     -18.228013d0,     -18.172579d0,     -18.121060d0,     -18.072532d0, &
            -18.026342d0,     -17.982247d0,     -17.940292d0,     -17.900424d0,     -17.862135d0,     -17.825336d0, &
            -31.704612d0,     -30.223451d0,     -28.903010d0,     -27.725958d0,     -26.677005d0,     -25.742728d0, &
            -24.911418d0,     -24.172777d0,     -23.517034d0,     -22.935075d0,     -22.418374d0,     -21.958737d0, &
            -21.548073d0,     -21.178541d0,     -20.843774d0,     -20.538724d0,     -20.259683d0,     -20.004394d0, &
            -19.771707d0,     -19.560869d0,     -19.371134d0,     -19.201427d0,     -19.050336d0,     -18.916100d0, &
            -18.796794d0,     -18.690526d0,     -18.595570d0,     -18.510375d0,     -18.433400d0,     -18.363316d0, &
            -18.299077d0,     -18.239902d0,     -18.185250d0,     -18.134609d0,     -18.087040d0,     -18.041874d0, &
            -17.998844d0,     -17.957966d0,     -17.919153d0,     -17.881871d0,     -17.845999d0,     -31.912062d0, &
            -30.429015d0,     -29.107037d0,     -27.928698d0,     -26.878558d0,     -25.942975d0,     -25.109922d0, &
            -24.368710d0,     -23.709330d0,     -23.122434d0,     -22.599202d0,     -22.131119d0,     -21.709862d0, &
            -21.327539d0,     -20.977787d0,     -20.655969d0,     -20.359309d0,     -20.086770d0,     -19.838360d0, &
            -19.614148d0,     -19.413701d0,     -19.235828d0,     -19.078739d0,     -18.940239d0,     -18.818009d0, &
            -18.709839d0,     -18.613761d0,     -18.528045d0,     -18.451015d0,     -18.381245d0,     -18.317623d0, &
            -18.259320d0,     -18.205758d0,     -18.156395d0,     -18.110271d0,     -18.066699d0,     -18.025393d0, &
            -17.986340d0,     -17.949422d0,     -17.914075d0,     -17.880136d0,     -32.254661d0,     -30.769926d0, &
            -29.446561d0,     -28.267049d0,     -27.215808d0,     -26.278964d0,     -25.444113d0,     -24.700021d0, &
            -24.036104d0,     -23.442218d0,     -22.908447d0,     -22.425059d0,     -21.982797d0,     -21.573691d0, &
            -21.192336d0,     -20.836315d0,     -20.505789d0,     -20.202383d0,     -19.927794d0,     -19.682692d0, &
            -19.466361d0,     -19.276845d0,     -19.111441d0,     -18.967123d0,     -18.840905d0,     -18.730078d0, &
            -18.632316d0,     -18.545632d0,     -18.468168d0,     -18.398373d0,     -18.335052d0,     -18.277323d0, &
            -18.224570d0,     -18.176226d0,     -18.131314d0,     -18.089140d0,     -18.049417d0,     -18.012131d0, &
            -17.977160d0,     -17.943930d0,     -17.912267d0,     -32.690100d0,     -31.204472d0,     -29.880366d0, &
            -28.700210d0,     -27.648309d0,     -26.710549d0,     -25.874032d0,     -25.126587d0,     -24.456135d0, &
            -23.850260d0,     -23.296248d0,     -22.782071d0,     -22.298336d0,     -21.840138d0,     -21.407424d0, &
            -21.003265d0,     -20.631524d0,     -20.295198d0,     -19.995620d0,     -19.732244d0,     -19.502903d0, &
            -19.304303d0,     -19.132636d0,     -18.984041d0,     -18.854930d0,     -18.742176d0,     -18.643166d0, &
            -18.555718d0,     -18.477840d0,     -18.407892d0,     -18.344621d0,     -18.287108d0,     -18.234714d0, &
            -18.186861d0,     -18.142560d0,     -18.101117d0,     -18.062250d0,     -18.025955d0,     -17.992117d0, &
            -17.960170d0,     -17.929948d0,     -33.167511d0,     -31.681535d0,     -30.357129d0,     -29.176665d0, &
            -28.124274d0,     -27.185369d0,     -26.345897d0,     -25.591306d0,     -24.905510d0,     -24.271026d0, &
            -23.671715d0,     -23.097636d0,     -22.547708d0,     -22.027091d0,     -21.542367d0,     -21.098458d0, &
            -20.697950d0,     -20.341450d0,     -20.027964d0,     -19.755094d0,     -19.519315d0,     -19.316369d0, &
            -19.141779d0,     -18.991222d0,     -18.860799d0,     -18.747177d0,     -18.647603d0,     -18.559806d0, &
            -18.481730d0,     -18.411695d0,     -18.348423d0,     -18.290978d0,     -18.238714d0,     -18.191043d0, &
            -18.146975d0,     -18.105815d0,     -18.067286d0,     -18.031388d0,     -17.998012d0,     -17.966595d0, &
            -17.936981d0,     -33.660115d0,     -32.174014d0,     -30.849468d0,     -29.668712d0,     -28.615410d0, &
            -27.673599d0,     -26.825899d0,     -26.051405d0,     -25.325785d0,     -24.627731d0,     -23.949258d0, &
            -23.297214d0,     -22.683578d0,     -22.117030d0,     -21.601444d0,     -21.137434d0,     -20.723931d0, &
            -20.358993d0,     -20.039984d0,     -19.763470d0,     -19.525277d0,     -19.320727d0,     -19.145070d0, &
            -18.993802d0,     -18.862909d0,     -18.748978d0,     -18.649205d0,     -18.561286d0,     -18.483143d0, &
            -18.413081d0,     -18.349812d0,     -18.292396d0,     -18.240180d0,     -18.192576d0,     -18.148594d0, &
            -18.107537d0,     -18.069130d0,     -18.033375d0,     -18.000164d0,     -17.968939d0,     -17.939543d0, &
            -34.157747d0,     -32.671584d0,     -31.346869d0,     -30.165427d0,     -29.109476d0,     -28.158925d0, &
            -27.287353d0,     -26.461349d0,     -25.653091d0,     -24.859471d0,     -24.098537d0,     -23.388719d0, &
            -22.739047d0,     -22.151052d0,     -21.622778d0,     -21.151151d0,     -20.732972d0,     -20.365096d0, &
            -20.044202d0,     -19.766458d0,     -19.527453d0,     -19.322365d0,     -19.146350d0,     -18.994845d0, &
            -18.863797d0,     -18.749767d0,     -18.649935d0,     -18.561985d0,     -18.483831d0,     -18.413773d0, &
            -18.350519d0,     -18.293127d0,     -18.240943d0,     -18.193379d0,     -18.149442d0,     -18.108437d0, &
            -18.070089d0,     -18.034400d0,     -18.001264d0,     -17.970121d0,     -17.940819d0,     -34.656985d0, &
            -33.170737d0,     -31.845576d0,     -30.662049d0,     -29.597967d0,     -28.621800d0,     -27.689172d0, &
            -26.762889d0,     -25.846135d0,     -24.970317d0,     -24.159726d0,     -23.422833d0,     -22.758739d0, &
            -22.162941d0,     -21.630316d0,     -21.156179d0,     -20.736501d0,     -20.367699d0,     -20.046210d0, &
            -19.768072d0,     -19.528797d0,     -19.323521d0,     -19.147377d0,     -18.995786d0,     -18.864684d0, &
            -18.750627d0,     -18.650789d0,     -18.562848d0,     -18.484717d0,     -18.414693d0,     -18.351480d0, &
            -18.294137d0,     -18.242007d0,     -18.194502d0,     -18.150629d0,     -18.109694d0,     -18.071420d0, &
            -18.035810d0,     -18.002758d0,     -17.971706d0,     -17.942499d0,     -35.156712d0,     -33.670226d0, &
            -32.343695d0,     -31.153715d0,     -30.065545d0,     -29.023568d0,     -27.975806d0,     -26.929574d0, &
            -25.931701d0,     -25.013131d0,     -24.181969d0,     -23.435224d0,     -22.766280d0,     -22.168018d0, &
            -21.634110d0,     -21.159294d0,     -20.739246d0,     -20.370227d0,     -20.048592d0,     -19.770338d0, &
            -19.530963d0,     -19.325604d0,     -19.149397d0,     -18.997770d0,     -18.866664d0,     -18.752632d0, &
            -18.652850d0,     -18.564992d0,     -18.486966d0,     -18.417066d0,     -18.353995d0,     -18.296807d0, &
            -18.244845d0,     -18.197516d0,     -18.153831d0,     -18.113095d0,     -18.075033d0,     -18.039645d0, &
            -18.006823d0,     -17.976010d0,     -17.947053d0,     -35.656526d0,     -34.169309d0,     -32.838528d0, &
            -31.629219d0,     -30.477726d0,     -29.310405d0,     -28.130104d0,     -27.001415d0,     -25.965095d0, &
            -25.030095d0,     -24.191860d0,     -23.441986d0,     -22.771656d0,     -22.172822d0,     -21.638738d0, &
            -21.163930d0,     -20.743959d0,     -20.375022d0,     -20.053441d0,     -19.775211d0,     -19.535846d0, &
            -19.330507d0,     -19.154357d0,     -19.002845d0,     -18.871921d0,     -18.758147d0,     -18.658693d0, &
            -18.571229d0,     -18.493654d0,     -18.424258d0,     -18.361735d0,     -18.305134d0,     -18.253790d0, &
            -18.207103d0,     -18.164087d0,     -18.124053d0,     -18.086721d0,     -18.052085d0,     -18.020027d0, &
            -17.990001d0,     -17.961847d0,     -36.156164d0,     -34.666712d0,     -33.323216d0,     -32.061656d0, &
            -30.783657d0,     -29.469006d0,     -28.200087d0,     -27.034150d0,     -25.983597d0,     -25.043055d0, &
            -24.202601d0,     -23.451884d0,     -22.781338d0,     -22.182604d0,     -21.648778d0,     -21.174299d0, &
            -20.754672d0,     -20.386068d0,     -20.064806d0,     -19.786903d0,     -19.547912d0,     -19.343045d0, &
            -19.167506d0,     -19.016777d0,     -18.886826d0,     -18.774214d0,     -18.676101d0,     -18.590135d0, &
            -18.514201d0,     -18.446571d0,     -18.385924d0,     -18.331284d0,     -18.281961d0,     -18.237336d0, &
            -18.196433d0,     -18.158575d0,     -18.123470d0,     -18.091083d0,     -18.061276d0,     -18.033519d0, &
            -18.007645d0 &
        ], [NCOOLT,NCOOLNE] )

   !--- Ground-term fine-structure lines ------------------------------!
   ! The forbidden lines emitted inside the split ground term of the C/N/O
   ! coolants. Level energies E/k [K], statistical weights and transition
   ! probabilities are CHIANTI v11.0.2 (elvlc / wgfa); level indices run in
   ! energy order, so for the inverted O I term level 1 is 3P2.
   !   C I  2p2 3P: g = 1, 3, 5      [C I]  609.1 / 370.4 um
   !   C II 2p  2P: g = 2, 4         [C II] 157.7 um
   !   N II 2p2 3P: g = 1, 3, 5      [N II] 205.3 / 121.8 um
   !   O I  2p4 3P: g = 5, 3, 1      [O I]   63.2 / 145.5 / 44.1 um
   ! The 3P2-3P0 lines of C I and N II have no CHIANTI transition
   ! probability (magnetic quadrupole, A < 1e-12 s^-1) and are carried as
   ! collisional couplings only.
   real*8, parameter :: Ek_CI2  =  23.6204d0, Ek_CI3  =  62.4631d0
   real*8, parameter :: Ek_CII2 =  91.2141d0
   real*8, parameter :: Ek_NII2 =  70.0684d0, Ek_NII3 = 188.1920d0
   real*8, parameter :: Ek_OI2  = 227.7080d0, Ek_OI3  = 326.5693d0
   real*8, parameter :: A_CI609  = 7.960d-8, A_CI370  = 2.660d-7
   real*8, parameter :: A_CII158 = 2.290d-6
   real*8, parameter :: A_NII205 = 2.080d-6, A_NII122 = 7.460d-6
   real*8, parameter :: A_OI63   = 8.542d-5, A_OI145  = 1.643d-5
   real*8, parameter :: A_OI44   = 1.380d-10

   ! Atomic weights [m_H] for the Doppler width of these lines
   ! Standard atomic weights in u, read from the species table (one
   ! definition, species_table.f90 melem_A_u); the fine-structure line
   ! opacities below use them for the Doppler width.
   real*8, parameter :: amu_C = melem_A_u(iel_C), amu_N = melem_A_u(iel_N), &
                        amu_O = melem_A_u(iel_O)

   ! Escape-probability slots, one per line with a transition probability.
   ! fine_structure_line_transfer (util_ion_eq.f90) fills them; the cooling
   ! coefficients consume them. Order: C I, C II, N II, O I.
   integer, parameter :: n_fsline    = 8
   integer, parameter :: ifs_CI609   = 1, ifs_CI370  = 2,              &
                         ifs_CII158  = 3,                              &
                         ifs_NII205  = 4, ifs_NII122 = 5,              &
                         ifs_OI63    = 6, ifs_OI145  = 7,              &
                         ifs_OI44    = 8
   ! emitting ion of each slot (canonical mion index)
   integer, parameter :: fsline_ion(n_fsline) =                        &
        [ im_CI, im_CI, im_CII, im_NII, im_NII, im_OI, im_OI, im_OI ]
   ! transition energy E_ul/k [K] of each slot, in the same order. Note the
   ! O I term is inverted, so slot ifs_OI44 (3P0-3P2) carries the full
   ! Ek_OI3 while ifs_OI145 (3P0-3P1) carries the difference.
   real*8, parameter :: fsline_Ek(n_fsline) =                          &
        [ Ek_CI2, Ek_CI3 - Ek_CI2, Ek_CII2,                            &
          Ek_NII2, Ek_NII3 - Ek_NII2,                                  &
          Ek_OI2, Ek_OI3 - Ek_OI2, Ek_OI3 ]

   !--- He(2^3S) + neutral: Penning / associative branching -----------!
   ! Collisional ionization of the He 2^3S metastable by a neutral partner
   ! runs through two channels that share the same total cross section,
   !   Penning:     He(2^3S) + X  -> He(1^1S) + X^+ + e^-
   !   associative: He(2^3S) + X  -> HeH^+ (+ fragment) + e^-,
   ! and Garcia Munoz (2025), A&A 698, A199, splits the total between them
   ! with "an average 0.9:0.1" (Sect. 2, the paragraph introducing
   ! Table A.5). The published network file applies exactly that split
   ! through the amplitudes of its two rows, for the H and the H2 partner
   ! alike (references/garcia_munoz_2025_network/SI_networkfile.txt lines
   ! 198/199 and 202/203). ioniz_HeI23S_H and ioniz_HeI23S_H2 below return
   ! the TOTAL, which is what removes the metastable; multiply by
   ! f_penning_HeI23S where a lasting ion and free electron are produced.
   real*8, parameter :: f_penning_HeI23S = 0.9d0

   !--- Electron-impact excitation rate coefficients -------------------!
   ! For a transition l -> u of energy E_ul with effective collision
   ! strength Upsilon (Maxwell-averaged), the excitation and de-excitation
   ! rate coefficients [cm^3 s^-1, T in K] are
   !   q_lu = coll_rate_prefactor Upsilon exp(-E_ul/kT) / (g_l sqrt(T)),
   !   q_ul = coll_rate_prefactor Upsilon / (g_u sqrt(T)),
   ! so q_ul/q_lu = (g_l/g_u) exp(E_ul/kT) (detailed balance) holds by
   ! construction, and the de-excitation form carries no Boltzmann factor
   ! that could overflow at low T. The constant is the Maxwell average of
   ! the Bethe-Born cross-section scale,
   !   h^2/((2 pi m_e)^(3/2) k^(1/2)) = 8.629132e-6 cm^3 s^-1 K^(1/2),
   ! formed here from the CODATA 2018 values of h, m_e and k (the literal
   ! 8.629e-6 it replaces was 1.5e-5 low). Oklopcic & Hirata (2018, ApJL
   ! 855, L11, eq. 16) write the same form as 2.1e-8 sqrt(13.6 eV/kT)
   ! Upsilon/omega_i, whose constant is this one rounded down:
   ! 8.629e-6/sqrt(Ry/k) = 2.172e-8, so 2.1e-8 is 3.3% low. Every
   ! Upsilon-form rate and cooling coefficient of this module uses this one
   ! constant.
   real*8, parameter :: coll_rate_prefactor =                             &
                        hp_erg**2/((2.0d0*pi*m_e)**1.5d0*sqrt(kb_erg))

   !--- He I levels of the metastable network -------------------------!
   ! Excitation energies above the ground 1^1S0 [eV], from the NIST ASD
   ! level energies as distributed in CHIANTI v11 he_1.elvlc: 2^3S1
   ! 159856.078, 2^1S0 166277.542, 2^1P1 171135.000 cm^-1. The 2^3S energy
   ! is written as the difference of the two ionization potentials of the
   ! code, e_th_HeI - e_th_HeTR = 19.819614 eV (159856.078 cm^-1 gives
   ! 19.819628 eV), so that ionizing the metastable after exciting it costs
   ! exactly e_th_HeI. Weights: g(1^1S) = 1, g(2^3S) = 3.
   real*8, parameter :: E_HeI_23S_eV     = e_th_HeI - e_th_HeTR
   real*8, parameter :: E_HeI_21S_23S_eV = 0.796160d0
   real*8, parameter :: E_HeI_21P_23S_eV = 1.398408d0

   ! Effective collision strengths of the three electron transitions of
   ! the network, Bray, Burgess, Fursa & Tully (2000, A&AS 146, 481),
   ! Table 2 (1^1S - 2^3S) and Table 3 (2^3S - 2^1S, 2^3S - 2^1P), at
   ! log10 T = 3.75 (0.25) 5.75, as printed; and, for the temperatures
   ! outside that table, the Burgess & Tully (1992) data of the same
   ! transitions in CHIANTI v11.0.2 he_1.scups (levels 1-2, 2-3, 2-7):
   ! scaling type, dE [Ry], C, and the five scaled knot values
   ! (upsilon_HeI_bray_table).
   real*8, parameter :: bray_logT_first = 3.75d0, bray_dlogT = 0.25d0
   real*8, parameter :: bray_ups_11S_23S(9) = [ 6.198d-2, 6.458d-2,       &
        6.387d-2, 6.157d-2, 5.832d-2, 5.320d-2, 4.787d-2, 4.018d-2,        &
        3.167d-2 ]
   real*8, parameter :: bray_ups_23S_21S(9) = [ 2.389d0, 2.456d0,         &
        2.275d0, 1.916d0, 1.496d0, 1.111d0, 8.003d-1, 5.660d-1, 3.944d-1 ]
   real*8, parameter :: bray_ups_23S_21P(9) = [ 7.965d-1, 9.579d-1,       &
        1.042d0, 1.015d0, 8.950d-1, 7.265d-1, 5.516d-1, 3.948d-1,          &
        2.677d-1 ]
   integer, parameter :: bt_type_11S_23S = 3, bt_type_23S_21S = 2,       &
                         bt_type_23S_21P = 3
   real*8, parameter :: bt_dE_11S_23S = 1.457d0, bt_C_11S_23S = 0.06d0
   real*8, parameter :: bt_dE_23S_21S = 0.05858d0, bt_C_23S_21S = 1.6d0
   real*8, parameter :: bt_dE_23S_21P = 0.1029d0, bt_C_23S_21P = 1.1d0
   real*8, parameter :: bt_knots_11S_23S(5) = [ 0.02482d0, 0.06481d0,     &
        0.07228d0, 0.07636d0, 0.08904d0 ]
   real*8, parameter :: bt_knots_23S_21S(5) = [ 0.7875d0, 2.222d0,        &
        2.317d0, 1.686d0, 0.3501d0 ]
   real*8, parameter :: bt_knots_23S_21P(5) = [ 0.1801d0, 1.132d0,        &
        2.257d0, 4.299d0, 8.113d0 ]

   !--- Energy removed from the electrons by a capture -----------------!
   ! For any capture process whose rate coefficient is a Maxwell average,
   ! alpha(T) = <sigma(E) v>, the kinetic energy carried off by the
   ! captured electrons is
   !   beta(T) = <sigma(E) v E> = k T alpha(T) (3/2 + dln alpha/dln T)
   ! [erg cm^3 s^-1]: differentiate alpha = C T^(-3/2) int sigma v E^(1/2)
   ! exp(-E/kT) dE with respect to T. It holds for radiative recombination
   ! into any set of levels (case A, case B, one level) and for dielectronic
   ! recombination alike. CHECKED against Hummer (1994, MNRAS 268, 109,
   ! Table 1): beta/alpha from his tabulated alpha_B and beta_B, and
   ! 3/2 + dln alpha_B/dln T from his alpha_B alone, agree to 0.1% over
   ! 1.6e3-6.3e4 K (case A to 0.2%). dlnT_capture is the half step of
   ! the centred logarithmic difference used for dln alpha/dln T
   ! (capture_energy_loss_rate); its truncation error is of order
   ! dlnT_capture^2/6 relative, 2e-7.
   real*8, parameter :: dlnT_capture = 1.0d-3

   !--- H I n = 2 --------------------------------------------------------!
   ! The 1s-2s excitation energy, 82258.954 cm^-1 (NIST ASD, as in CHIANTI
   ! v11 h_1.elvlc) = 10.19881 eV; the two 2p levels lie within 0.37 cm^-1
   ! of 2s, so this one energy serves the 1s-2s and 1s-2p channels.
   real*8, parameter :: E_HI_n2_eV = 10.19881d0

   !--- Tables of the capture coefficients built from quadratures ------!
   ! ln alpha and dln alpha/dln T of the ground captures of H I, He I and
   ! He II (Milne relation, ground_capture_milne_moments) and of the direct
   ! capture into n = 2 of a hydrogenic ion in the scaled temperature T/Z^2
   ! (Seaton 1959, seaton_n2_moments), on log10 T = 0 (0.005) 8. They are
   ! filled ONCE, serially, by capture_coefficient_tables_init after the
   ! input is read (the He I cross section depends on "ATES
   ! photoionization rate"), and read by cubic Hermite interpolation in
   ! (ln T, ln alpha) with the tabulated slopes, whose error at this node
   ! spacing is below 1e-9 relative (MEASURED against the quadrature,
   ! physics probe ground_capture_and_cascade_data). Until they are filled
   ! (a test driver that does not call the init) every value is the
   ! quadrature itself.
   integer, parameter :: n_capture_tab = 1601
   real*8,  parameter :: capture_tab_dlog10T = 0.005d0
   real*8,  save :: capture_tab_lnalpha(n_capture_tab,4) = 0.0d0
   real*8,  save :: capture_tab_slope(n_capture_tab,4)   = 0.0d0
   logical, save :: capture_tables_ready = .false.

   !--- Channels of the radiative cooling of a cell -------------------!
   ! Columns of radiative_cooling_of_cell and of eval_cool's breakdown:
   ! 1-6 the atomic and ionic channels, 7-10 the molecular ones, 10+i the
   ! line cooling of metal ion i (see radiative_cooling_of_cell).
   integer, parameter :: n_cool_chan = 10 + n_mion

   contains

   !---------------------------------------------------!

   ! =============================================================== !
   !  Default H/He recombination and collisional-ionization rates.   !
   !  Used unless legacy_hhe_rates = .true. (then the legacy HG97 /   !
   !  Abel fits below are taken instead) or atomic_rate_set_k22 =     !
   !  .true. (then the Koskinen et al. 2022 Table 1 fits of           !
   !  mol_rates are taken; the accessors near the end of this module  !
   !  are where all three sets meet).                                 !
   ! =============================================================== !

   ! Badnell radiative recombination fit (2023 update of Badnell 2006,
   ! ApJS 167, 334; https://amdpp.phys.strath.ac.uk/tamoc/RR).
   elemental double precision function rr_badnell(T,A,B,T0,T1,Cc,T2)
   real*8, intent(in) :: T,A,B,T0,T1,Cc,T2
   real*8 :: bp,tt
   bp = B + Cc*exp(-T2/T)
   tt = sqrt(T/T0)
   rr_badnell = A/( tt*(1.0d0+tt)**(1.0d0-bp)*(1.0d0+sqrt(T/T1))**(1.0d0+bp) )
   end function rr_badnell

   ! Badnell dielectronic recombination He II -> He I (adf09 three-term sum;
   ! negligible below ~5e4 K).
   elemental double precision function dr_HeII_badnell(T)
   real*8, intent(in) :: T
   dr_HeII_badnell = T**(-1.5d0)*( 1.417d-3*exp(-4.633d5/T)   &
                                 + 2.235d-4*exp(-5.532d5/T)   &
                                 - 2.185d-5*exp(-8.887d5/T) )
   end function dr_HeII_badnell

   ! THE CAPTURE DIRECTLY INTO THE GROUND LEVEL, alpha_1 [cm^3 s^-1], from
   ! the Milne relation on the ground-state photoionization cross section
   ! the radiative transfer of this code absorbs with (Cross_sections:
   ! sigma for H I and He II, sigma_HeI for He I), so that these captures
   ! and the absorption of the photons they emit are the two directions of
   ! ONE process. Detailed balance between photoionization of a level of
   ! weight g_n and radiative recombination onto it from an ion of weight
   ! g_+ gives
   !   alpha_1(T) = (g_n/g_+) sqrt(2/pi) / (c^2 (m_e k T)^(3/2))
   !                x INT_0^inf (I + E)^2 sigma(I + E) exp(-E/kT) dE ,
   ! I the ionization potential and E the kinetic energy of the captured
   ! electron (the Milne relation, Milne 1924, in the form of Hummer 1994,
   ! MNRAS 268, 109, eq. 1, and Hummer & Storey 1998, MNRAS 297, 1073,
   ! eq. 7; DERIVED here with the weights written out). Weights:
   ! H I 1s 2S (2) from H+ (1); He I 1s2 1S (1) from He+ 1s 2S (2); He II
   ! 1s 2S (2) from He2+ (1).
   !
   ! QUADRATURE. With E = I (exp(t) - 1) the integrand is smooth in t at
   ! every temperature, and a 24-point Gauss-Legendre rule on
   ! 0 <= t <= ln(1 + 50 kT/I) (the rest of the Maxwellian, e^-50, is
   ! below double precision) reproduces an adaptive quadrature of the same
   ! integral to 1.4e-10 relative over 10 K - 3e7 K for all three ions
   ! (MEASURED). CHECKS (MEASURED): H I, 0.998 of Hummer (1994, Table 1)
   ! alpha_1 over 10 K - 1e7 K (the transfer's cross section carries the
   ! 6.3e-18 cm^2 threshold value and the measured, reduced-mass
   ! potential); He II, 0.999 of 2 alpha_1(T/4) of the same table; He I,
   ! 1.006-1.011 of the 1s2 1S coefficient of Hummer & Storey (1998,
   ! Table 5) over 10 K - 2.5e4 K (their R-matrix cross section against
   ! the Verner et al. 1996 fit near threshold). VALIDITY: that of the
   ! cross section, every temperature the code reaches.
   elemental double precision function ground_capture_milne(T,ion)
   real*8,  intent(in) :: T
   integer, intent(in) :: ion      ! 1 = H I, 2 = He I, 3 = He II
   real*8 :: a, dlna
   if (capture_tables_ready) then
      ground_capture_milne = capture_table_value(T, ion)
   else
      call ground_capture_milne_moments(T, ion, a, dlna)
      ground_capture_milne = a
   endif
   end function ground_capture_milne

   ! The quadrature itself: alpha_1 and its logarithmic slope,
   ! dln alpha_1/dln T = -3/2 + <E/kT>, the mean of E/kT under the same
   ! integrand, from one pass of the rule.
   pure subroutine ground_capture_milne_moments(T, ion, alpha, dlnalpha)
   real*8,  intent(in)  :: T
   integer, intent(in)  :: ion
   real*8,  intent(out) :: alpha, dlnalpha
   ! 24-point Gauss-Legendre rule on [-1, 1]: the 12 positive nodes and
   ! their weights (the rule is symmetric).
   real*8, parameter :: x_gl(12) = [ 6.40568928626056300d-02,             &
        1.91118867473616311d-01, 3.15042679696163397d-01,                 &
        4.33793507626045127d-01, 5.45421471388839563d-01,                 &
        6.48093651936975546d-01, 7.40124191578554358d-01,                 &
        8.20001985973902947d-01, 8.86415527004400960d-01,                 &
        9.38274552002732798d-01, 9.74728555971309474d-01,                 &
        9.95187219997021311d-01 ]
   real*8, parameter :: w_gl(12) = [ 1.27938195346752215d-01,             &
        1.25837456346828303d-01, 1.21670472927803419d-01,                 &
        1.15505668053725613d-01, 1.07444270115965607d-01,                 &
        9.76186521041140648d-02, 8.61901615319532882d-02,                 &
        7.33464814110804109d-02, 5.92985849154367417d-02,                 &
        4.42774388174195510d-02, 2.85313886289337432d-02,                 &
        1.23412297999870909d-02 ]
   real*8  :: I_eV, g_ratio, kT_eV, t_max, t_q, E_ph, term, acc0, acc1
   integer :: k, sgn
   select case (ion)
      case (1)
         I_eV = e_th_HI;   g_ratio = 2.0d0
      case (2)
         I_eV = e_th_HeI;  g_ratio = 0.5d0
      case default
         I_eV = e_th_HeII; g_ratio = 2.0d0
   end select
   kT_eV = kb_eV*max(T, 1.0d0)
   t_max = log(1.0d0 + 50.0d0*kT_eV/I_eV)
   acc0  = 0.0d0
   acc1  = 0.0d0
   do k = 1,12
      do sgn = -1,1,2
         t_q  = 0.5d0*t_max*(1.0d0 + dble(sgn)*x_gl(k))
         E_ph = I_eV*exp(t_q)
         ! (I + E)^2 sigma(I + E) exp(-E/kT) dE/dt, dE/dt = I e^t = E_ph
         term = w_gl(k)*E_ph**3*ground_cross_section(E_ph, ion)          &
                *exp(-(E_ph - I_eV)/kT_eV)
         acc0 = acc0 + term
         acc1 = acc1 + term*(E_ph - I_eV)/kT_eV
      enddo
   enddo
   ! eV^3 x 1e-18 cm^2 -> erg^3 cm^2, then the Milne prefactor
   alpha = g_ratio*sqrt(2.0d0/pi)                                         &
        /(c_light**2*(m_e*kT_eV/erg2eV)**1.5d0)                           &
        *0.5d0*t_max*acc0*1.0d-18/erg2eV**3
   dlnalpha = -1.5d0 + acc1/max(acc0, tiny(1.0d0))
   end subroutine ground_capture_milne_moments

   ! The ground-state photoionization cross section of the transfer
   ! [1e-18 cm^2] at photon energy E [eV], ion as in ground_capture_milne.
   pure double precision function ground_cross_section(E,ion)
   real*8,  intent(in) :: E
   integer, intent(in) :: ion
   select case (ion)
      case (1)
         ground_cross_section = sigma(E, 1.0d0, e_th_HI)
      case (2)
         ground_cross_section = sigma_HeI(E)
      case default
         ground_cross_section = sigma(E, 2.0d0, e_th_HeII)
   end select
   end function ground_cross_section

   !--- Hydrogenic ions: the n = 2 level after a case-B recombination ---!

   ! DIRECT CAPTURE INTO n = 2 of a hydrogenic ion of nuclear charge Z
   ! [cm^3 s^-1], Seaton (1959, MNRAS 119, 81, eqs. 5-10, READ from the
   ! published article):
   !   alpha_n = D Z lambda^(1/2) (x_n/n) S_n(lambda),
   !   D = 5.197e-14 cm^3 s^-1,  lambda = 157890 Z^2/T,  x_n = lambda/n^2,
   !   S_n = INT_0^inf g_II(n,u) exp(-x_n u)/(1 + u) du,
   !   g_II = 1 + 0.1728 n^(-2/3) (u+1)^(-2/3) (u-1)
   !          - 0.0496 n^(-4/3) (u+1)^(-4/3) (u^2 + 4u/3 + 1),
   ! the three-term asymptotic expansion of the bound-free Gaunt factor.
   ! With u = exp(t) - 1 the integrand is g_II exp(-x_n u) dt and the
   ! 24-point rule of ground_capture_milne on 0 <= t <= ln(1 + 50/x_n)
   ! integrates it. Its photon, of energy I_2 + E_k (I_2 = Z^2 Ry/4), is
   ! the only case-B recombination continuum of He II that can ionize
   ! H I: I_2(He II) = 13.6047 eV against 13.5984 eV. ACCURACY: Hummer
   ! (1994, section 4) finds Seaton's alpha_1 and alpha_B within 0.1-2.4%
   ! of the exact values for T/Z^2 < 1e4 K and within 1-5% for 1e4 <
   ! T/Z^2 < 1e6 K; for n = 1 this form gives 0.977 of Hummer's alpha_1
   ! at 1e4 K (MEASURED).
   elemental double precision function alpha_n2_hydrogenic_seaton(T,Z)
   real*8, intent(in) :: T, Z
   real*8 :: a, dlna
   if (capture_tables_ready) then
      alpha_n2_hydrogenic_seaton = Z*capture_table_value(T/(Z*Z), 4)
   else
      call seaton_n2_moments(T/(Z*Z), a, dlna)
      alpha_n2_hydrogenic_seaton = Z*a
   endif
   end function alpha_n2_hydrogenic_seaton

   ! The quadrature for Z = 1 at the scaled temperature Ts = T/Z^2 (alpha_n
   ! scales as Z alpha_n(T/Z^2) in Seaton's form, lambda = 157890 Z^2/T):
   ! alpha_2, and dln alpha_2/dln T = -3/2 + (x_n/S) INT g_II u
   ! exp(-x_n u)/(1 + u) du (lambda^1/2 x_n goes as T^-3/2, and
   ! dln S/dln T = -x_n dln S/dx_n), from one pass of the rule, in which
   ! du = (1 + u) dt cancels the 1/(1 + u).
   pure subroutine seaton_n2_moments(Ts, alpha, dlnalpha)
   real*8, intent(in)  :: Ts
   real*8, intent(out) :: alpha, dlnalpha
   real*8, parameter :: x_gl(12) = [ 6.40568928626056300d-02,             &
        1.91118867473616311d-01, 3.15042679696163397d-01,                 &
        4.33793507626045127d-01, 5.45421471388839563d-01,                 &
        6.48093651936975546d-01, 7.40124191578554358d-01,                 &
        8.20001985973902947d-01, 8.86415527004400960d-01,                 &
        9.38274552002732798d-01, 9.74728555971309474d-01,                 &
        9.95187219997021311d-01 ]
   real*8, parameter :: w_gl(12) = [ 1.27938195346752215d-01,             &
        1.25837456346828303d-01, 1.21670472927803419d-01,                 &
        1.15505668053725613d-01, 1.07444270115965607d-01,                 &
        9.76186521041140648d-02, 8.61901615319532882d-02,                 &
        7.33464814110804109d-02, 5.92985849154367417d-02,                 &
        4.42774388174195510d-02, 2.85313886289337432d-02,                 &
        1.23412297999870909d-02 ]
   real*8, parameter :: n_lev = 2.0d0
   real*8  :: lam, x_n, t_max, t_q, u, g_II, term, acc0, acc1
   integer :: k, sgn
   lam   = 157890.0d0/max(Ts, 1.0d-2)
   x_n   = lam/n_lev**2
   t_max = log(1.0d0 + 50.0d0/x_n)
   acc0  = 0.0d0
   acc1  = 0.0d0
   do k = 1,12
      do sgn = -1,1,2
         t_q  = 0.5d0*t_max*(1.0d0 + dble(sgn)*x_gl(k))
         u    = exp(t_q) - 1.0d0
         g_II = 1.0d0 + 0.1728d0*n_lev**(-2.0d0/3.0d0)                   &
                *(u + 1.0d0)**(-2.0d0/3.0d0)*(u - 1.0d0)                  &
              - 0.0496d0*n_lev**(-4.0d0/3.0d0)*(u + 1.0d0)**(-4.0d0/3.0d0) &
                *(u*u + 4.0d0*u/3.0d0 + 1.0d0)
         term = w_gl(k)*g_II*exp(-x_n*u)
         acc0 = acc0 + term
         acc1 = acc1 + term*u
      enddo
   enddo
   alpha    = 5.197d-14*sqrt(lam)*(x_n/n_lev)*0.5d0*t_max*acc0
   dlnalpha = -1.5d0 + x_n*acc1/max(acc0, tiny(1.0d0))
   end subroutine seaton_n2_moments

   ! Cubic Hermite interpolation of ln alpha in ln T on the capture tables,
   ! column k (1-3 the ground captures, 4 the Seaton n = 2 capture at
   ! T/Z^2), with the tabulated logarithmic slopes; linear continuation at
   ! the two ends.
   pure double precision function capture_table_value(T, k)
   real*8,  intent(in) :: T
   integer, intent(in) :: k
   real*8  :: pos, f, h, h00, h10, h01, h11, lnT0
   integer :: i
   pos = log10(max(T, 1.0d-30))/capture_tab_dlog10T + 1.0d0
   h   = capture_tab_dlog10T*log(10.0d0)
   if (pos .le. 1.0d0) then
      lnT0 = 0.0d0
      capture_table_value = exp(capture_tab_lnalpha(1,k)                  &
           + capture_tab_slope(1,k)*(log(max(T, 1.0d-30)) - lnT0))
   else if (pos .ge. dble(n_capture_tab)) then
      lnT0 = dble(n_capture_tab - 1)*h
      capture_table_value = exp(capture_tab_lnalpha(n_capture_tab,k)      &
           + capture_tab_slope(n_capture_tab,k)*(log(T) - lnT0))
   else
      i   = int(pos)
      f   = pos - dble(i)
      h00 = (1.0d0 + 2.0d0*f)*(1.0d0 - f)**2
      h10 = f*(1.0d0 - f)**2
      h01 = f*f*(3.0d0 - 2.0d0*f)
      h11 = f*f*(f - 1.0d0)
      capture_table_value = exp(h00*capture_tab_lnalpha(i,k)              &
           + h10*h*capture_tab_slope(i,k)                                 &
           + h01*capture_tab_lnalpha(i+1,k)                               &
           + h11*h*capture_tab_slope(i+1,k))
   endif
   end function capture_table_value

   ! Fill the capture tables (serially, once, after the input is read).
   subroutine capture_coefficient_tables_init()
   integer :: i, ion
   real*8  :: Tn, a, dlna
   do i = 1, n_capture_tab
      Tn = 10.0d0**(dble(i-1)*capture_tab_dlog10T)
      do ion = 1,3
         call ground_capture_milne_moments(Tn, ion, a, dlna)
         capture_tab_lnalpha(i,ion) = log(a)
         capture_tab_slope(i,ion)   = dlna
      enddo
      call seaton_n2_moments(Tn, a, dlna)
      capture_tab_lnalpha(i,4) = log(a)
      capture_tab_slope(i,4)   = dlna
   enddo
   capture_tables_ready = .true.
   end subroutine capture_coefficient_tables_init

   ! FRACTION OF THE CASE-B RECOMBINATIONS OF A HYDROGENIC ION THAT END IN
   ! 2s (the rest end in 2p), alpha_2s/(alpha_2s + alpha_2p) from the
   ! effective case-B coefficients of Pengelly (1964, MNRAS 127, 145,
   ! Table I, case B, READ from the published article), tabulated as
   ! alpha(t')/Z at t' = 1e-4 T/Z^2, log2 t' = -3 (1) 3, so one table
   ! serves H I (Z = 1, 1250-80000 K) and He II (Z = 2, 5000-320000 K);
   ! interpolated linearly in log t' of the logarithms of the two
   ! coefficients and held at the table ends. For H I it gives 0.284,
   ! 0.293, 0.322 and 0.354 at 4e3, 5e3, 1e4 and 2e4 K (MEASURED); Draine
   ! (2011, Table 14.3, from Brown & Mathews 1970) lists f(2s) = 0.285,
   ! 0.325 and 0.356 at 4e3, 1e4 and 2e4 K, and his Table 14.2 (from
   ! Hummer & Storey 1987, n_e = 1e3 cm^-3) 0.305 at 5e3 K. Draine gives
   ! no fitting formula; the quadratic 0.282 + 0.047 T4 - 0.006 T4^2 of
   ! Christie et al. (2013, Table 2, R8), attributed there to Draine,
   ! lies 3.8% above this table at 5e3 K and turns over above 4e4 K
   ! (0.274 at 8e4 K against 0.414). For He II, Osterbrock & Ferland
   ! (2006, Table 2.6) imply 0.28 at 1e4 K and 0.31 at 2e4 K (DERIVED),
   ! against 0.264 and 0.293 here.
   ! The effective 2s coefficient of Pengelly at t' = 0.25 (He II at
   ! 1e4 K) equals the Storey & Hummer (1995) value at n_e = 1e2 cm^-3 to
   ! 0.3% (4.08e-13 against 4.091e-13 cm^3 s^-1; the latter as tabulated in
   ! the MoCHII data file heii_caseB_n2.txt, used here only as a check).
   ! Low-density values: the l-mixing of 2s by ion impact that moves a
   ! 2s atom to 2p at higher density is applied by the caller
   ! (l_mixing_2s2p_pengelly_seaton).
   elemental double precision function case_b_2s_fraction_hydrogenic(T,Z)
   real*8, intent(in) :: T, Z
   real*8, parameter :: a2s(7) = [ 30.6d0, 20.4d0, 13.3d0, 8.37d0,        &
                                   5.07d0, 2.93d0, 1.61d0 ]
   real*8, parameter :: a2p(7) = [ 97.5d0, 56.8d0, 32.1d0, 17.6d0,        &
                                   9.27d0, 4.68d0, 2.28d0 ]
   real*8  :: pos, f, l2s, l2p
   integer :: k
   pos = log(1.0d-4*max(T, 1.0d0)/(Z*Z))/log(2.0d0) + 4.0d0
   if (.not. (pos > 1.0d0)) then
      k = 1;  f = 0.0d0
   else if (pos .ge. 7.0d0) then
      k = 6;  f = 1.0d0
   else
      k = int(pos); f = pos - dble(k)
   endif
   l2s = exp(log(a2s(k)) + f*(log(a2s(k+1)) - log(a2s(k))))
   l2p = exp(log(a2p(k)) + f*(log(a2p(k+1)) - log(a2p(k))))
   case_b_2s_fraction_hydrogenic = l2s/(l2s + l2p)
   end function case_b_2s_fraction_hydrogenic

   ! 2s -> 2p COLLISIONAL TRANSFER BY A CHARGED PARTICLE in a hydrogenic
   ! ion of nuclear charge z [cm^3 s^-1], Pengelly & Seaton (1964, MNRAS
   ! 127, 165, READ from the published article): the impact-parameter
   ! cross section of their eqs. (34)-(37),
   !   Q(v) = pi R1^2 [1/2 + ln(Rc/R1)],
   !   R1^2 = 6 (Z_c/z)^2 (e^2/(hbar v))^2 n^2 (n^2 - l^2 - l - 1) a_0^2
   ! (their eq. 41, n = 2, l = 0), with the cut-off at large impact
   ! parameter "whichever is the smallest" of their three (sect. 4):
   !   Rc = min(1.12 hbar v/dE, 0.72 v tau, R_D),
   ! dE the energy separation of initial and final state (their eq. 29),
   ! tau the radiative lifetime of 2s (eq. 30) and R_D = (kT/(4 pi n_e
   ! e^2))^(1/2) the Debye radius (6.90 cm at 1e4 K, 1e4 cm^-3, as they
   ! quote). THE CUT-OFF THAT APPLIES. For 2s -> 2p the Lamb shift and the
   ! fine structure separate the final states: 2s1/2 -> 2p1/2 by
   ! 0.035 cm^-1 and 2s1/2 -> 2p3/2 by 0.331 cm^-1 in H I, 0.468 and
   ! 5.389 cm^-1 in He II (NIST ASD energies as distributed in CHIANTI v11
   ! h_1.elvlc and he_2.elvlc), which share the line strength 1:2. With
   ! these, 1.12 hbar v/dE is the smallest cut-off up to n_e ~ 4e12 cm^-3
   ! (H I, protons, 1e4 K), and P&S say so for these levels ("for the lower
   ! states with l = 0 and l = 1, it would be rather more accurate to use
   ! (29) for the cut-off (Seaton 1955)", their footnote to sect. 8); the
   ! Debye cut-off alone would give 2.9e-3 (n_e = 1e2) to 1.4e-3 (n_e =
   ! 1e10) cm^3 s^-1 for protons at 1e4 K, against 4.79e-4 with all three
   ! (MEASURED). Where Rc < R1 the component contributes pi Rc^2/2 (P = 1/2
   ! inside Rc and nothing beyond; P&S do not treat that corner, which
   ! carries a Maxwellian weight below 1e-4 at 1e4 K). The Maxwell average
   ! of Q v is written in closed form (lower incomplete gamma functions and
   ! the exponential integral E1 over the three velocity ranges the
   ! minimum defines) and agrees with an adaptive quadrature to 5e-5
   ! (MEASURED). Results (MEASURED): H I by protons 5.16e-4 (3e3 K), 4.79e-4
   ! (1e4 K), 4.26e-4 (2e4 K) cm^3 s^-1; by He+ 5.04e-4 and by He2+
   ! 1.46e-3 at 1e4 K. Osterbrock & Ferland (2006, Table 4.10, p. 87,
   ! taken from Pengelly & Seaton 1964) list for H I 2s -> 2p1/2 and 2p3/2
   ! protons 2.51e-4 + 2.23e-4 = 4.74e-4 (1e4 K) and 2.08e-4 + 2.19e-4 =
   ! 4.27e-4 (2e4 K), electrons 0.22e-4 + 0.35e-4 and 0.17e-4 + 0.27e-4,
   ! and totals 5.31e-4 and 4.71e-4 cm^3 s^-1, the sums of the proton and
   ! electron rows (their eq. 4.27 applies them as n_p q^p + n_e q^e); the
   ! same numbers are Seaton's (1955, Proc. Phys. Soc. A 68, 457, Table 1,
   ! approximation V). THE CUT-OFF RULE AGAINST SEATON (1955): he treats
   ! "separately the 2s -> 2p1/2 and 2s -> 2p3/2 transitions", with
   ! "energy differences ... dE' = 0.0354 cm-1 and dE'' = 0.327 cm-1
   ! (Lamb 1951)" and Omega' = Omega(dE')/3, Omega'' = 2 Omega(dE'')/3
   ! (eq. 18), the separation entering the logarithm ln(2 M v^2/dE) of his
   ! eq. (55); that is the large-distance cut-off 1.12 hbar v/dE of P&S
   ! eq. (29) with the same two separations and weights used here. Seaton
   ! has no Debye cut-off; P&S's minimum over the three cut-offs adds it
   ! for dense gas. Z_c the collider charge, mu_g the
   ! reduced mass [g], dE1, dE2 the two separations [erg] with weights 1/3
   ! and 2/3, A_2q the 2s radiative rate [s^-1]. VALIDITY: P&S's
   ! impact-parameter theory, which needs R1 well outside the atom; P&S
   ! judge the rates "quite reliable" because Rc enters only
   ! logarithmically.
   pure double precision function l_mixing_2s2p_pengelly_seaton          &
                         (T, n_e, z, Z_c, mu_g, dE1, dE2, A_2q) result(q)
   real*8, intent(in) :: T, n_e, z, Z_c, mu_g, dE1, dE2, A_2q
   real*8, parameter :: e_esu = 1.602176634d-19*c_light/10.0d0
   real*8, parameter :: hbar  = hp_erg/(2.0d0*pi)
   real*8, parameter :: a_bohr = hbar**2/(m_e*e_esu**2)
   real*8  :: Tk, R_D, K_1, v_th, a_c, x_1, x_a, x_D, A0, B0, comp
   real*8  :: dEs(2), wts(2)
   integer :: jc
   Tk   = max(T, 1.0d0)
   R_D  = sqrt(kb_erg*Tk/(4.0d0*pi*max(n_e, 1.0d-30)*e_esu**2))
   K_1  = sqrt(6.0d0*(Z_c/z)**2*12.0d0)*e_esu**2/hbar*a_bohr
   v_th = sqrt(2.0d0*kb_erg*Tk/mu_g)
   dEs  = [dE1, dE2]
   wts  = [1.0d0/3.0d0, 2.0d0/3.0d0]
   q    = 0.0d0
   do jc = 1,2
      a_c = min(1.12d0*hbar/dEs(jc), 0.72d0/A_2q)
      x_1 = K_1/(a_c*v_th**2)
      x_a = R_D**2/(a_c*v_th)**2
      x_D = (K_1/(R_D*v_th))**2
      B0  = 0.5d0 + log(R_D*v_th/K_1)
      if (x_1 .lt. x_a) then
         A0   = 0.5d0 + log(a_c*v_th**2/K_1)
         comp = 0.5d0*pi*a_c**2*v_th**3*lower_gamma_3(x_1)               &
              + pi*K_1**2/v_th*( A0*(exp(-x_1) - exp(-x_a))              &
                + log(x_1)*exp(-x_1) - log(x_a)*exp(-x_a)                &
                + exponential_integral_E1(x_1)                           &
                - exponential_integral_E1(x_a) )                         &
              + pi*K_1**2/v_th*( B0*exp(-x_a)                            &
                + 0.5d0*(log(x_a)*exp(-x_a) + exponential_integral_E1(x_a)) )
      else
         comp = 0.5d0*pi*a_c**2*v_th**3*lower_gamma_3(x_a)               &
              + 0.5d0*pi*R_D**2*v_th                                     &
                *(lower_gamma_2(x_D) - lower_gamma_2(x_a))              &
              + pi*K_1**2/v_th*( B0*exp(-x_D)                            &
                + 0.5d0*(log(x_D)*exp(-x_D) + exponential_integral_E1(x_D)) )
      endif
      q = q + wts(jc)*comp
   enddo
   q = 2.0d0/sqrt(pi)*q
   end function l_mixing_2s2p_pengelly_seaton

   ! Lower incomplete gamma functions gamma(3,x) = 2 - e^-x (x^2 + 2x + 2)
   ! and gamma(2,x) = 1 - e^-x (1 + x), by their series below x = 1e-2
   ! (where the closed forms lose digits to cancellation).
   pure double precision function lower_gamma_3(x)
   real*8, intent(in) :: x
   if (x .lt. 1.0d-2) then
      lower_gamma_3 = x**3*(1.0d0/3.0d0 - x/4.0d0 + x*x/10.0d0)
   else
      lower_gamma_3 = 2.0d0 - exp(-x)*(x*x + 2.0d0*x + 2.0d0)
   endif
   end function lower_gamma_3

   pure double precision function lower_gamma_2(x)
   real*8, intent(in) :: x
   if (x .lt. 1.0d-2) then
      lower_gamma_2 = x**2*(0.5d0 - x/3.0d0 + x*x/8.0d0)
   else
      lower_gamma_2 = 1.0d0 - exp(-x)*(1.0d0 + x)
   endif
   end function lower_gamma_2

   ! Exponential integral E1(x) = INT_x^inf e^-t/t dt, x > 0: the power
   ! series -gamma - ln x - sum (-x)^k/(k k!) for x <= 1 and the continued
   ! fraction of Abramowitz & Stegun (5.1.22) evaluated by the modified
   ! Lentz method above, both to 1e-15 relative.
   pure double precision function exponential_integral_E1(x) result(e1)
   real*8, intent(in) :: x
   real*8, parameter :: euler_gamma = 0.5772156649015329d0
   real*8  :: term, b, c, d, h, del
   integer :: k
   if (x .le. 0.0d0) then
      e1 = huge(1.0d0)
   else if (x .le. 1.0d0) then
      term = 1.0d0
      e1   = -euler_gamma - log(x)
      do k = 1,60
         term = -term*x/dble(k)
         del  = -term/dble(k)
         e1   = e1 + del
         if (abs(del) .lt. 1.0d-16*abs(e1)) exit
      enddo
   else
      b = x + 1.0d0
      c = 1.0d0/1.0d-300
      d = 1.0d0/b
      h = d
      do k = 1,200
         term = -dble(k)*dble(k)
         b    = b + 2.0d0
         d    = 1.0d0/(term*d + b)
         c    = b + term/c
         del  = c*d
         h    = h*del
         if (abs(del - 1.0d0) .lt. 1.0d-15) exit
      enddo
      e1 = h*exp(-x)
   endif
   end function exponential_integral_E1

   ! Case-B coefficients: the total recombination coefficient alpha_A of
   ! Badnell (radiative, and dielectronic for He II -> He I) minus the
   ! ground capture of the Milne relation above. alpha_A and alpha_1 are
   ! therefore from two calculations, but alpha_1 is the one the
   ! absorption of the same photons implies, which is what the case-B
   ! subtraction asserts. MEASURED: alpha_B(He I) = 2.763e-13 at 1e4 K
   ! against 2.755e-13 of Hummer & Storey (1998, Table 5); alpha_B(H) and
   ! alpha_B(He II) against Hummer (1994, Table 1): see the physics probe
   ! recombination_coefficient_fits.
   elemental double precision function alphaB_HII_new(T)
   real*8, intent(in) :: T
   alphaB_HII_new =                                                          &
        rr_badnell(T, 8.318d-11, 0.7472d0, 2.965d0, 7.001d5, 0.0d0, 0.0d0)   &
      - ground_capture_milne(T, 1)
   end function alphaB_HII_new

   elemental double precision function alphaB_HeII_new(T)
   real*8, intent(in) :: T
   alphaB_HeII_new =                                                         &
        rr_badnell(T, 5.235d-11, 0.6988d0, 7.301d0, 4.475d6, 0.0829d0,       &
                   1.682d5)                                                  &
      + dr_HeII_badnell(T)                                                   &
      - ground_capture_milne(T, 2)
   end function alphaB_HeII_new

   elemental double precision function alphaB_HeIII_new(T)
   real*8, intent(in) :: T
   alphaB_HeIII_new =                                                        &
        rr_badnell(T, 1.818d-10, 0.7492d0, 1.017d1, 2.786d6, 0.0d0, 0.0d0)   &
      - ground_capture_milne(T, 3)
   end function alphaB_HeIII_new

   ! Voronov (1997, ADNDT 65, 1, eq. 1) collisional-ionization fit
   ! [cm^3 s^-1]:
   !   k = A (1 + P sqrt(U)) U^K exp(-U) / (X + U),  U = dE[eV]/T[eV],
   ! one row (dE, P, A, X, K) per ion in his Table I, which also gives
   ! "Tmin, Tmax  Electron temperature range over which the fit has been
   ! made, with the lower bound listed in eV and the upper in keV": 1 eV
   ! (11605 K) to 20 keV for every row the code carries except He II
   ! (3 eV), N II and O II (2 eV). The fit reproduces the recommended
   ! (Belfast) rates "to within a few percentage points" (abstract); "The
   ! deviation of the fit values from the recommended data, less than 10%
   ! in most cases, is far less than the estimated accuracy of the
   ! recommended data, which is 40 - 60%" (Conclusion); at Te = 1 eV the
   ! fit is 95% (H I), 97% (He I), 105% (C I) and 101% (C II) of them
   ! (Table II). BELOW Tmin, the photoionized winds' 2e3-1e4 K, the code
   ! evaluates the same expression: it keeps the exact exp(-U) threshold
   ! factor and continues the fitted prefactor U^(K-1) (P = 0) or
   ! U^(K-1/2) (P = 1) to larger U, which neither Voronov nor his input
   ! data constrain there; the rate is then exp(-U)-small (H I 6.9e-23
   ! cm^3 s^-1 at 5e3 K), and no accuracy is claimed for it.
   elemental double precision function voronov_ci(T,dE_eV,P,A,X,Kexp)
   real*8, intent(in) :: T,dE_eV,P,A,X,Kexp
   real*8 :: U
   U = dE_eV/(kb_eV*T)
   voronov_ci = A*(1.0d0+P*sqrt(U))*U**Kexp*exp(-U)/(X+U)
   end function voronov_ci

   ! dE below is Voronov's own Table I entry, a COEFFICIENT OF HIS FIT and not
   ! the measured ionization potential of global_parameters: it enters the
   ! functional form through U^K, exp(-U) and X + U, and A, X, K, P were
   ! fitted with it.  It is left as published, as the fit's other parameters
   ! are; the measured potentials differ from it by 0.01 (H I), 0.05 (He I)
   ! and 0.03 (He II) per cent, far inside the fit's own accuracy.  The energy
   ! REMOVED per ionization is a different quantity and is the measured
   ! potential e_th_*_erg (util_ion_eq, T_equation).
   elemental double precision function ci_HI_new(T)
   real*8, intent(in) :: T
   ci_HI_new = voronov_ci(T, 13.6d0, 0.0d0, 2.91d-8, 0.232d0, 0.39d0)
   end function ci_HI_new

   elemental double precision function ci_HeI_new(T)
   real*8, intent(in) :: T
   ci_HeI_new = voronov_ci(T, 24.6d0, 0.0d0, 1.75d-8, 0.180d0, 0.35d0)
   end function ci_HeI_new

   elemental double precision function ci_HeII_new(T)
   real*8, intent(in) :: T
   ci_HeII_new = voronov_ci(T, 54.4d0, 1.0d0, 2.05d-9, 0.265d0, 0.25d0)
   end function ci_HeII_new

   ! Thermally-averaged free-free Gaunt factor from the van Hoof et al. (2014,
   ! MNRAS 444, 420) table gff_avg, evaluated at ion net charge zion.
   !   x = log10(gamma^2), gamma^2 = zion^2 Ry/(k T) (their eq. 21), with
   !   Ry/k = Ry_over_kB, the infinite-mass Rydberg of their section 2.
   ! Clamped to the tabulated range, linearly interpolated in x.
   elemental double precision function gbar_ff(T,zion)
   real*8, intent(in) :: T,zion
   real*8 :: x,frac
   integer :: i
   x = log10(zion*zion*Ry_over_kB/T)
   x = min(max(x, gff_lg2_min), gff_lg2_min+(NGFF-1)*gff_lg2_step)
   frac = (x-gff_lg2_min)/gff_lg2_step
   i = min(int(frac)+1, NGFF-1)
   frac = frac - dble(i-1)
   gbar_ff = gff_avg(i)*(1.0d0-frac) + gff_avg(i+1)*frac
   end function gbar_ff

   !---------------------------------------------------!

   !--- Recombination ---!

   ! Case B recombination coefficient of HII
   subroutine rec_HII_B_range(T,coeff_rec_HII_B,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: coeff_rec_HII_B
   coeff_rec_HII_B(j_lo:j_hi) = alpha_rec_HII_B(T(j_lo:j_hi))
   end subroutine rec_HII_B_range


   !--------------!

   ! Case B recombination coefficient of HeII
   subroutine rec_HeII_B_range(T,coeff_rec_HeII_B,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: coeff_rec_HeII_B
   coeff_rec_HeII_B(j_lo:j_hi) = alpha_rec_HeII_B(T(j_lo:j_hi))
   end subroutine rec_HeII_B_range


   !--------------!

   ! Case B recombination coefficient of HeIII
   subroutine rec_HeIII_B_range(T,coeff_rec_HeIII_B,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: coeff_rec_HeIII_B
   coeff_rec_HeIII_B(j_lo:j_hi) = alpha_rec_HeIII_B(T(j_lo:j_hi))
   end subroutine rec_HeIII_B_range


   !---------------------------------------------------!

   !--- Collisional ionization ---!

   ! Collisional ionization rate for HI
   subroutine ion_coeff_HI_range(T,a_ion_coeff_HI,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: a_ion_coeff_HI
   a_ion_coeff_HI(j_lo:j_hi) = ci_rate_HI(T(j_lo:j_hi))
   end subroutine ion_coeff_HI_range


   !--------------!

   ! Collisional ionization rate for HeI
   subroutine ion_coeff_HeI_range(T,a_ion_coeff_HeI,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: a_ion_coeff_HeI
   a_ion_coeff_HeI(j_lo:j_hi) = ci_rate_HeI(T(j_lo:j_hi))
   end subroutine ion_coeff_HeI_range


   !--------------!

   ! Collisional ionization rate for HeII
   subroutine ion_coeff_HeII_range(T,a_ion_coeff_HeII,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: a_ion_coeff_HeII
   a_ion_coeff_HeII(j_lo:j_hi) = ci_rate_HeII(T(j_lo:j_hi))
   end subroutine ion_coeff_HeII_range

   
   !---------------------------------------------------!
   
   !--- Electron collisions among the He I levels of the network -------!
   ! The metastable network (Oklopcic & Hirata 2018, ApJL 855, L11,
   ! eqs. 14-16) couples the ground 1^1S (weight 1), the metastable 2^3S
   ! (weight 3) and, through it, the singlets 2^1S and 2^1P by electron
   ! impact:
   !   q13  1^1S -> 2^3S  excitation,   E = E_HeI_23S_eV     = 19.8196 eV
   !   q31g 2^3S -> 1^1S  de-excitation, the detailed-balance reverse of q13
   !   q31a 2^3S -> 2^1S  excitation,   E = E_HeI_21S_23S_eV = 0.7962 eV
   !   q31b 2^3S -> 2^1P  excitation,   E = E_HeI_21P_23S_eV = 1.3984 eV
   ! (2^1S and 2^1P then decay radiatively to 1^1S). The rate form is the
   ! one written at coll_rate_prefactor, with the exact Maxwellian constant
   ! and the level energies of the NIST ASD. Oklopcic & Hirata leave q31g
   ! out as "less probable than collisional transitions to the excited
   ! singlet states"; that holds at 1e4 K, where q31g/(q31a + q31b) = 0.055,
   ! but not in cooler gas: 0.153 at 5e3 K, 0.563 at 3e3 K, 2.71 at 2e3 K
   ! (MEASURED with the strengths below).
   !
   ! EFFECTIVE COLLISION STRENGTHS: Bray, Burgess, Fursa & Tully (2000,
   ! A&AS 146, 481), convergent close-coupling, Table 2 (1^1S - 2^3S) and
   ! Table 3 (2^3S - 2^1S, 2^3S - 2^1P), the collision strengths Oklopcic
   ! & Hirata (2018) and Lampon et al. (2020, A&A 636, A13, Table 2: "Upsilon
   ! ij are the effective collision strengths taken from Bray et al. (2000)
   ! at the corresponding temperatures") name. Lampon et al. print no fit.
   ! The table covers "3.75 <= log T <= 5.75", and the authors "limit the
   ! low temperature end of our tabulation to about 6 000 degrees" because
   ! of "the uncertainty that pseudo resonances introduce into our cross
   ! sections at energies close to threshold" (p. 483); inside it "the error
   ! due to this threshold effect does not exceed a few percent" (Sect. 6).
   ! upsilon_HeI_bray_table evaluates, with the data declared at
   ! bray_ups_11S_23S:
   !   5.62e3 - 5.62e5 K: the table, interpolated by a cubic Hermite in
   !     (log10 T, log10 Upsilon) whose interior node slopes are the
   !     Fritsch-Carlson (PCHIP) ones and whose two end node slopes are
   !     those of the CHIANTI curve below; it passes through every tabulated
   !     value to round-off (physics probe helium_bray_collision_strengths);
   !   outside the table: the temperature dependence of the CHIANTI v11.0.2
   !     he_1 strength of the same transition (Burgess & Tully scaled knots,
   !     descaled with the natural spline of CHIANTI's DESCALE_SCUPS.PRO,
   !     chianti_bt_upsilon_5knot), scaled to the tabulated value at the
   !     table end, so that the strength and its slope are continuous there.
   ! WHY THIS BELOW 5.62e3 K, where the LHS 1140 b He 10830 line forms
   ! (1.5e3-6e3 K). Bray et al. give nothing there. CHIANTI states it uses
   ! its he_1 strengths (Sawey & Berrington 1993, ADNDT 55, 81, R-matrix;
   ! Bray et al. 2000) over 3.0 < log T < 5.75, so below log T = 3.75 its
   ! curve rests on the R-matrix data and on the threshold value of the
   ! scaled fit (Sawey & Berrington was not read); its magnitude differs
   ! from Bray's table by +7% (1^1S - 2^3S), -4% (2^3S - 2^1S) and +2%
   ! (2^3S - 2^1P) at 5.62e3 K (MEASURED), which the scaling removes. The
   ! alternative, extrapolating a fit through the table, has no data
   ! behind it: the two-exponential fits a exp(bT) + c exp(dT) this code
   ! carried before (coefficients from the ATES code, which calls them a
   ! "fit of Bray (2000)"; they pass exactly through the four table values
   ! 5.62e3-3.16e4 K) are 1.16, 1.12 and 1.10 times the strengths used
   ! here at 1.5e3 K, 1.12, 1.10, 1.06 at 2e3 K, 1.01, 1.01, 1.00 at 5e3 K,
   ! and fall away from the table above 3.2e4 K (2^3S - 2^1S: 0.73 of it
   ! at 1e5 K) (MEASURED, cooling_data/fit_helium_triplet_cascade.py). The
   ! strengths stay positive and finite as T -> 0 (0.0232, 0.821, 0.177
   ! at 1 K, MEASURED), so the rates are evaluated unclamped; below 1e3 K,
   ! where CHIANTI's range ends, they carry no accuracy claim.
   elemental double precision function upsilon_HeI_11S_23S(T)
   real*8, intent(in) :: T
   upsilon_HeI_11S_23S = upsilon_HeI_bray_table(bray_ups_11S_23S,        &
        bt_type_11S_23S, bt_dE_11S_23S, bt_C_11S_23S, bt_knots_11S_23S, T)
   end function upsilon_HeI_11S_23S

   elemental double precision function upsilon_HeI_23S_21S(T)
   real*8, intent(in) :: T
   upsilon_HeI_23S_21S = upsilon_HeI_bray_table(bray_ups_23S_21S,        &
        bt_type_23S_21S, bt_dE_23S_21S, bt_C_23S_21S, bt_knots_23S_21S, T)
   end function upsilon_HeI_23S_21S

   elemental double precision function upsilon_HeI_23S_21P(T)
   real*8, intent(in) :: T
   upsilon_HeI_23S_21P = upsilon_HeI_bray_table(bray_ups_23S_21P,        &
        bt_type_23S_21P, bt_dE_23S_21P, bt_C_23S_21P, bt_knots_23S_21P, T)
   end function upsilon_HeI_23S_21P

   ! One strength of the network from the Bray et al. (2000) table
   ! (9 values at log10 T = bray_logT_first + (i-1) bray_dlogT) and, outside
   ! it, the CHIANTI curve (itype, dE_Ry, C_bt, knots) scaled to the table
   ! end (see upsilon_HeI_11S_23S). Node slopes are in d log10 Upsilon per
   ! table interval; the CHIANTI end slopes are centered differences over
   ! +-1e-3 dex of a spline that is smooth there (truncation 1e-7
   ! relative).
   pure double precision function upsilon_HeI_bray_table(tab, itype,     &
                                    dE_Ry, C_bt, knots, T) result(ups)
   real*8,  intent(in) :: tab(9), dE_Ry, C_bt, knots(5), T
   integer, intent(in) :: itype
   real*8, parameter :: hl = 1.0d-3
   real*8  :: L(9), pos, f, lt_end, m0, m1, sl, sr, h00, h10, h01, h11
   integer :: k, i, inode(2)
   real*8  :: mnode(2)
   L   = log10(tab)
   pos = (log10(max(T, 1.0d-300)) - bray_logT_first)/bray_dlogT
   if (.not. (pos .ge. 0.0d0)) then
      lt_end = bray_logT_first
      ups = tab(1)*chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots, T) &
            /chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots,          &
                                      10.0d0**lt_end)
   else if (pos .gt. 8.0d0) then
      lt_end = bray_logT_first + 8.0d0*bray_dlogT
      ups = tab(9)*chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots, T) &
            /chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots,          &
                                      10.0d0**lt_end)
   else
      k = min(int(pos), 7) + 1               ! interval [k, k+1]
      f = pos - dble(k - 1)
      inode = [k, k + 1]
      do i = 1,2
         if (inode(i) .eq. 1 .or. inode(i) .eq. 9) then
            lt_end = bray_logT_first + dble(inode(i) - 1)*bray_dlogT
            mnode(i) = bray_dlogT/(2.0d0*hl)*log10(                         &
                 chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots,       &
                                          10.0d0**(lt_end + hl))           &
                /chianti_bt_upsilon_5knot(itype, dE_Ry, C_bt, knots,       &
                                          10.0d0**(lt_end - hl)))
         else
            sl = L(inode(i)) - L(inode(i) - 1)
            sr = L(inode(i) + 1) - L(inode(i))
            if (sl*sr .le. 0.0d0) then
               mnode(i) = 0.0d0
            else
               mnode(i) = 2.0d0/(1.0d0/sl + 1.0d0/sr)
            endif
         endif
      enddo
      m0 = mnode(1)
      m1 = mnode(2)
      h00 = (1.0d0 + 2.0d0*f)*(1.0d0 - f)**2
      h10 = f*(1.0d0 - f)**2
      h01 = f**2*(3.0d0 - 2.0d0*f)
      h11 = f**2*(f - 1.0d0)
      ups = 10.0d0**(h00*L(k) + h10*m0 + h01*L(k+1) + h11*m1)
   endif
   end function upsilon_HeI_bray_table

   ! Effective collision strength of one CHIANTI .scups transition stored
   ! as five Burgess & Tully (1992, A&A 254, 436) scaled knot values at
   ! s = 0, 1/4, 1/2, 3/4, 1, scaling types 2 and 3 only:
   !   x = kT/dE, s = x/(x + C), Upsilon = S(s) (type 2) or S(s)/(x + 1)
   !   (type 3),
   ! S the natural cubic spline through the knots (second derivative zero
   ! at s = 0 and 1), the interpolation of CHIANTI's DESCALE_SCUPS.PRO
   ! (cooling_data/chianti_cooling.py: upsilon, scaled_upsilon_spline, which
   ! the physics probe helium_bray_collision_strengths reproduces). dE in
   ! Rydberg, CODATA 2018 (13.605693122994 eV).
   pure double precision function chianti_bt_upsilon_5knot(itype, dE_Ry, &
                                    C_bt, y, T) result(ups)
   integer, intent(in) :: itype
   real*8,  intent(in) :: dE_Ry, C_bt, y(5), T
   real*8, parameter :: Ry_eV = 13.605693122994d0, h = 0.25d0
   real*8  :: x, s, r(3), M(5), w2, w3, d2, d3, t0, t1, S_s
   integer :: k
   x = max(kb_eV*T, 0.0d0)/(dE_Ry*Ry_eV)
   s = min(x/(x + C_bt), 1.0d0)
   ! natural spline: M(1) = M(5) = 0, M(k-1) + 4 M(k) + M(k+1) =
   ! 6 (y(k-1) - 2 y(k) + y(k+1))/h^2 for k = 2, 3, 4 (Thomas algorithm)
   r  = 6.0d0*(y(1:3) - 2.0d0*y(2:4) + y(3:5))/h**2
   w2 = 4.0d0 - 0.25d0
   d2 = r(2) - 0.25d0*r(1)
   w3 = 4.0d0 - 1.0d0/w2
   d3 = r(3) - d2/w2
   M(1) = 0.0d0
   M(5) = 0.0d0
   M(4) = d3/w3
   M(3) = (d2 - M(4))/w2
   M(2) = (r(1) - M(3))/4.0d0
   k  = min(int(s/h), 3) + 1                 ! interval [k, k+1]
   t0 = dble(k - 1)*h
   t1 = t0 + h
   S_s = M(k)*(t1 - s)**3/(6.0d0*h) + M(k+1)*(s - t0)**3/(6.0d0*h)      &
       + (y(k)/h - M(k)*h/6.0d0)*(t1 - s)                                &
       + (y(k+1)/h - M(k+1)*h/6.0d0)*(s - t0)
   if (itype .eq. 3) then
      ups = S_s/(x + 1.0d0)
   else
      ups = S_s
   endif
   end function chianti_bt_upsilon_5knot

   ! q13 [cm^3 s^-1]
   elemental double precision function excitation_rate_HeI_11S_23S(T)
   real*8, intent(in) :: T
   excitation_rate_HeI_11S_23S = coll_rate_prefactor/sqrt(T)             &
                               *upsilon_HeI_11S_23S(T)                   &
                               *exp(-E_HeI_23S_eV/(kb_eV*T))
   end function excitation_rate_HeI_11S_23S

   ! q31g [cm^3 s^-1]: q13 (g_1/g_3) exp(E/kT) with the Boltzmann factor
   ! cancelled analytically.
   elemental double precision function deexcitation_rate_HeI_23S_11S(T)
   real*8, intent(in) :: T
   deexcitation_rate_HeI_23S_11S = coll_rate_prefactor/(3.0d0*sqrt(T))    &
                                 *upsilon_HeI_11S_23S(T)
   end function deexcitation_rate_HeI_23S_11S

   ! q31a [cm^3 s^-1]
   elemental double precision function excitation_rate_HeI_23S_21S(T)
   real*8, intent(in) :: T
   excitation_rate_HeI_23S_21S = coll_rate_prefactor/(3.0d0*sqrt(T))      &
                               *upsilon_HeI_23S_21S(T)                   &
                               *exp(-E_HeI_21S_23S_eV/(kb_eV*T))
   end function excitation_rate_HeI_23S_21S

   ! q31b [cm^3 s^-1]
   elemental double precision function excitation_rate_HeI_23S_21P(T)
   real*8, intent(in) :: T
   excitation_rate_HeI_23S_21P = coll_rate_prefactor/(3.0d0*sqrt(T))      &
                               *upsilon_HeI_23S_21P(T)                   &
                               *exp(-E_HeI_21P_23S_eV/(kb_eV*T))
   end function excitation_rate_HeI_23S_21P

   ! Grid forms of the four rates (HeITR_coeffs).
   subroutine coex_HeI_1S_23S(T,coeff_coex_HeI_1S_23S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_1S_23S
   coeff_coex_HeI_1S_23S = excitation_rate_HeI_11S_23S(T)
   end subroutine coex_HeI_1S_23S

   subroutine deexc_HeI_23S_1S(T,coeff_deexc_HeI_23S_1S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_deexc_HeI_23S_1S
   coeff_deexc_HeI_23S_1S = deexcitation_rate_HeI_23S_11S(T)
   end subroutine deexc_HeI_23S_1S

   subroutine coex_HeI_23S_21S(T,coeff_coex_HeI_23S_21S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_23S_21S
   coeff_coex_HeI_23S_21S = excitation_rate_HeI_23S_21S(T)
   end subroutine coex_HeI_23S_21S

   subroutine coex_HeI_23S_21P(T,coeff_coex_HeI_23S_21P)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_23S_21P
   coeff_coex_HeI_23S_21P = excitation_rate_HeI_23S_21P(T)
   end subroutine coex_HeI_23S_21P


   !--------------!

   ! Collisional-excitation cooling of the He 2^3S metastable within the
   ! triplet system [erg cm^3 s^-1 per (n_e n_2^3S)]: every excitation
   ! 2^3S -> n^3L (2^3P, 3^3S, 3^3P, 3^3D, ... 5^3G, the triplet terms
   ! n <= 5 of Bray, Burgess, Fursa & Tully 2000, A&AS 146, 481, Table 3,
   ! with the strengths of upsilon_HeI_11S_23S: the table inside
   ! 10^3.75-10^5.75 K, the CHIANTI v11 he_1 shape scaled to it outside,
   ! held at the table value for 2^3S - 4^3S, which CHIANTI lacks)
   ! cascades back to 2^3S inside the triplets (the 2^3P -> 2^3S step is
   ! the 10830 A line)
   ! and radiates E(n^3L) - E(2^3S), leaving the metastable population
   ! unchanged. Coronal sum with observed level energies,
   !   Lambda = coll_rate_prefactor/(3 sqrt(T)) E_P Upsilon_c(T)
   !            exp(-E_P/kT),
   ! E_P = E(2^3P) - E(2^3S) = 9230.94 cm^-1 = 13281.3 K (g-weighted over
   ! 2^3P_0,1,2, NIST ASD as distributed in he_1.elvlc) and Upsilon_c the
   ! energy-weighted effective strength sum_u Upsilon_2u (E_u2/E_P)
   ! exp(-(E_u2 - E_P)/kT), log10 a quartic in log10(T/1e4 K) fitted over
   ! 1e3-1e5 K to 0.976-1.051, 0.978-1.026 over 3e3-3e4 K
   ! (cooling_data/fit_helium_triplet_cascade.py;
   ! T held inside that range by upsilon_quartic_logT, and the Boltzmann
   ! factor carries the temperature outside it). The Black (1981, MNRAS
   ! 197, 553, Table 3) coefficient 1.16e-20 T^1/2
   ! exp(-13179/T), which is his n = 2, 3, 4 triplet cooling rate divided
   ! by his steady-state n(2^3S) (his eq. 11): its 13179 K is a fitted
   ! threshold 0.8% below the 2^3P energy, and the whole coefficient is
   ! 0.74 (3e3 K), 0.81 (5e3 K), 0.80 (1e4 K), 0.82 (2e4 K) of the sum
   ! (MEASURED).
   elemental double precision function cooling_HeI_23S_triplets(T)
   real*8, intent(in) :: T
   real*8, parameter :: E_P_K = 13281.257d0
   cooling_HeI_23S_triplets = coll_rate_prefactor/(3.0d0*sqrt(T))         &
        *E_P_K*kb_erg*upsilon_quartic_logT(1.43825681d0, 9.76614520d-01,   &
        2.17740740d-03, -2.43394109d-01, 6.42716758d-02, T)               &
        *exp(-E_P_K/T)
   end function cooling_HeI_23S_triplets

   ! THE COLLISIONAL FEED OF 2^3S THROUGH THE HIGHER TRIPLETS [cm^3 s^-1]:
   ! electron-impact excitation 1^1S -> n^3L, every triplet term n <= 5
   ! except 2^3S itself (Bray et al. 2000, Table 2, with the strengths of
   ! upsilon_HeI_11S_23S inside and outside the table). A triplet level
   ! cannot decay to 1^1S on any competing time scale, so each such
   ! excitation ends in 2^3S after its radiative cascade within the
   ! triplets, and it enters the triplet balance exactly as the direct
   ! 1^1S -> 2^3S excitation does (HeITR_coeffs adds it to q13). Coronal sum,
   !   q_feed = coll_rate_prefactor/sqrt(T) Upsilon_f(T) exp(-E_f/kT),
   ! E_f = E(2^3P) = 169087.01 cm^-1 = 243278.5 K above 1^1S (g-weighted,
   ! observed) and Upsilon_f = sum_u Upsilon_1u exp(-(E_u - E_f)/kT),
   ! log10 a quartic in log10(T/1e4 K) fitted over 1e3-1e5 K to
   ! 0.985-1.010 (cooling_data/fit_helium_triplet_cascade.py). Against
   ! the direct 1^1S -> 2^3S rate the feed is 0.0028 (3e3 K), 0.019
   ! (5e3 K), 0.106 (1e4 K), 0.334 (2e4 K), 0.544 (3e4 K) (MEASURED on the
   ! Bray sums; the fitted function below gives 0.019, 0.107 and 0.331 at
   ! 5e3, 1e4 and 2e4 K, inside its stated accuracy; with the CHIANTI he_1
   ! strengths throughout the ratios are 0.0026, 0.018, 0.101, 0.328,
   ! 0.531). ENERGY: nothing is added to the cooling for it. The
   ! He I sum lambda_coex_HeI_except_23S already charges each such
   ! excitation its full E(n^3L); of that, E(n^3L) - E(2^3S) is radiated
   ! in the cascade and E(2^3S) is carried by the metastable it makes,
   ! which the metastable's own exits return or radiate exactly as for a
   ! direct excitation. NO REVERSE IS ADDED: the detailed-balance reverse
   ! of each channel is n^3L -> 1^1S de-excitation from the upper level,
   ! whose population (radiative lifetimes 1e-7 - 1e-8 s) is negligible;
   ! the reverse of the direct channel, q31g, stays that of 1^1S -> 2^3S.
   elemental double precision function excitation_rate_HeI_11S_triplets(T)
   real*8, intent(in) :: T
   real*8, parameter :: E_f_K = 243278.485d0
   excitation_rate_HeI_11S_triplets = coll_rate_prefactor/sqrt(T)        &
        *upsilon_quartic_logT(-1.58458877d0, 6.67140463d-01,               &
        -7.34789211d-03, -1.17759228d-01, -5.49937576d-02, T)             &
        *exp(-E_f_K/T)
   end function excitation_rate_HeI_11S_triplets

   ! THE SPIN-CHANGING EXCITATION OF 2^3S TO THE SINGLETS ABOVE 2^1P
   ! (3^1S, 3^1P, 3^1D, 4^1L, ... 5^1G, Bray et al. 2000, Table 3, with
   ! the strengths of upsilon_HeI_11S_23S)
   ! [cm^3 s^-1], which q31a (2^1S) and q31b (2^1P) leave out; each such
   ! atom cascades to 1^1S, so the rate destroys the metastable
   ! (HeITR_coeffs adds it to q31b, whose exits it shares: the cascade
   ! ends in 2^1P -> 1^1S or in a direct n^1P -> 1^1S resonance line, both
   ! above the H I edge, and is given the 584 A photon of q31b). Coronal
   ! sum, observed energies,
   !   q = coll_rate_prefactor/(3 sqrt(T)) Upsilon_s(T) exp(-E_s/kT),
   ! E_s = E(3^1S) - E(2^3S) = 35982.2 K, the lowest of the levels, and
   ! Upsilon_s = sum_u Upsilon_2u exp(-(E_u2 - E_s)/kT), log10 a quartic in
   ! log10(T/1e4 K) fitted over 1e3-1e5 K to 0.987-1.025
   ! (cooling_data/fit_helium_triplet_cascade.py). Against q31a + q31b it
   ! is 0.0013 (5e3 K), 0.020 (1e4 K), 0.079 (2e4 K), 0.128 (3e4 K)
   ! (MEASURED; the CHIANTI he_1 strengths, 14-20% above the table for
   ! these terms at 5.6e3-3e4 K, give 0.0016, 0.023, 0.090, 0.14).
   elemental double precision function excitation_rate_HeI_23S_singlets_n3(T)
   real*8, intent(in) :: T
   real*8, parameter :: E_s_K = 35982.161d0
   excitation_rate_HeI_23S_singlets_n3 = coll_rate_prefactor              &
        /(3.0d0*sqrt(T))*upsilon_quartic_logT(-8.52139935d-02,             &
        1.88951687d-01, -2.50604990d-01, -1.30265154d-01, 5.78292939d-02, T)&
        *exp(-E_s_K/T)
   end function excitation_rate_HeI_23S_singlets_n3

   ! Its energy [erg cm^3 s^-1 per (n_e n_2^3S)]: the excitation energy
   ! E(n^1L) - E(2^3S) each such collision takes from the electrons, the
   ! stored E(2^3S) having been charged with q13; energy-weighted strength
   ! Upsilon_e = sum_u Upsilon_2u (E_u2/E_s) exp(-(E_u2 - E_s)/kT), fitted
   ! as above to 0.985-1.029.
   elemental double precision function cooling_HeI_23S_singlets_n3(T)
   real*8, intent(in) :: T
   real*8, parameter :: E_s_K = 35982.161d0
   cooling_HeI_23S_singlets_n3 = coll_rate_prefactor/(3.0d0*sqrt(T))      &
        *E_s_K*kb_erg*upsilon_quartic_logT(-5.36819985d-02,                 &
        2.27762762d-01, -2.62734620d-01, -1.45042104d-01, 6.68587628d-02, T)&
        *exp(-E_s_K/T)
   end function cooling_HeI_23S_singlets_n3


   !--------------!

   ! Recombination coefficient of HeII on 23S state
   subroutine rec_HeII_23S(T,coeff_rec_HeII_23S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HeII_23S

   coeff_rec_HeII_23S = alpha_rec_HeII_23S(T)

   end subroutine rec_HeII_23S

   !--------------!

   ! TOTAL ionization rate coefficient of the He 2^3S metastable in
   ! collisions with atomic hydrogen,
   !   He(2^3S) + H -> He(1^1S) + H^+ + e^-   (Penning, 90%)
   !                -> HeH^+ + e^-            (associative, 10%),
   ! as the Maxwell-Boltzmann average of the Movre & Meyer (1997) total
   ! ionization cross sections.  Closed form published by Garcia Munoz
   ! (2025), A&A 698, A199, Fig. 5, and carried by that paper's network file
   ! (references/garcia_munoz_2025_network/SI_networkfile.txt, rows 198 and
   ! 199, label 'movre97', amplitudes 0.9e-9 and 0.1e-9 summing to the 1e-9
   ! written here):
   !   k = 1e-9 exp(c/T + d1 lnT + d2 (lnT)^2 + d3 (lnT)^3)  [cm^3 s^-1].
   ! NB the header of that network file writes the third term as exp(c), but
   ! only the exp(c/T) of the Fig. 5 caption reproduces the paper's own
   ! Table A.5; read c as c/T.
   !
   ! Validity: Garcia Munoz quotes errors of <2% against the cross-section
   ! average "from 200 to 10 000 K", and the expression reproduces the four
   ! points of Table A.5 (1.02e-9, 1.32e-9, 1.35e-9, 1.27e-9 at 500, 2000,
   ! 5000 and 10 000 K) to 0.4%.  It is evaluated unclamped outside that
   ! range: it stays smooth, positive and monotonic on either side of its
   ! 3.8e3 K maximum (6.93e-10 at 200 K, 1.19e-9 at 1.5e4 K, 1.11e-9 at
   ! 2e4 K), so clamping would only add a kink, but the <2% accuracy is not
   ! claimed there.  Below ~200 K it falls off faster than the underlying
   ! calculation, which Garcia Munoz reports as ~7e-10 cm^3 s^-1 at 100 K
   ! against 4.1e-10 from this expression; that region is far below the
   ! temperatures a photoionized wind reaches.
   !
   ! Replaces the two-branch power law of Taylor et al. (2025), ApJ 989:68,
   ! Table 2 (1.9e-9 (300/T)^0.07 for T <= 4000 K, 9.1e-9 (300/T)^0.50
   ! above), which was transcribed correctly but is not usable: it jumps by
   ! a factor 1.5724 at 4000 K (1.5849e-9 -> 2.4921e-9), and a rate
   ! coefficient is a continuous function of temperature.  That fit also
   ! disagrees with its own paper's Figure 19, where the Maxwell-Boltzmann
   ! average of the Cohen & Lane (1971) and Morgner & Niehaus (1979) cross
   ! sections rises to 1.50e-9 near 1400 K and falls to 1.29e-9 at 4000 K,
   ! while the tabulated branch decreases monotonically from 1.83e-9 to
   ! 1.58e-9; and above 4000 K it exceeds every cross-section determination
   ! collected in Garcia Munoz (2025), Fig. 5.  Figure 19 itself agrees with
   ! the expression used here to 10-20%, so the two independent calculations
   ! are consistent and it is the published fit that was the outlier.
   ! Before Taylor et al. the code used the temperature-independent 5e-10
   ! cm^3 s^-1 of Roberge & Dalgarno (1982), as adopted by Oklopcic & Hirata
   ! (2018) and Lampon et al. (2020).
   ! Units: cm^3 s^-1; T in K (gas temperature: this is a neutral-neutral
   ! collision).
   elemental double precision function ioniz_HeI23S_H(T) result(k)
   real*8, intent(in) :: T
   real*8 :: lnT

   lnT = log(T)
   k = 1.0d-9*exp(-8.64804d1/T - 2.86766d-1*lnT                       &
                  + 8.68445d-2*lnT**2 - 5.73001d-3*lnT**3)

   end function ioniz_HeI23S_H

   !--------------!

   ! TOTAL ionization rate coefficient of the He 2^3S metastable in
   ! collisions with molecular hydrogen,
   !   He(2^3S) + H2 -> He(1^1S) + H2^+ + e^-   (Penning, 90%)
   !                 -> H + HeH^+ + e^-         (associative, 10%),
   ! as the Maxwell-Boltzmann average of the Cohen & Lane (1977) cross
   ! sections.  Coefficients taken from the Garcia Munoz (2025) network file
   ! (references/garcia_munoz_2025_network/SI_networkfile.txt, rows 202 and
   ! 203, label 'cohenl77'), whose two amplitudes 4.86740e-12 and
   ! 5.40822e-13 sum to the total written here and realise the same 0.9:0.1
   ! split as the atomic channel:
   !   k = a T^b exp(c/T)  [cm^3 s^-1].
   ! It reproduces that paper's Table A.5 (8.94e-11, 6.48e-10, 1.48e-9,
   ! 2.54e-9 cm^3 s^-1 at 500, 2000, 5000 and 10 000 K) to 0.13%.  Validity
   ! follows the tabulation, 500-10 000 K; it is evaluated unclamped outside
   ! that range, where it stays smooth and positive but carries no accuracy
   ! claim.
   !
   ! This replaces a fit made here to the four Table A.5 points, which was
   ! accurate but reproduced the TOTAL while the code then charged all of it
   ! to the Penning channel, overstating H2^+ production from this reaction
   ! by 1/0.9 = 1.111.
   ! Units: cm^3 s^-1; T in K (gas temperature).  Elemental so it serves
   ! both the scalar calls in System_HeH_mol::set_mol_coeffs and the array
   ! evaluation of the heating term in ionization_equilibrium.
   elemental double precision function ioniz_HeI23S_H2(T) result(k)
   real*8, intent(in) :: T
   k = 5.408222d-12 * T**6.75388d-1 * exp(-6.96275d2/T)
   end function ioniz_HeI23S_H2

   !--------------!

   ! Electron-impact collisional ionization of the He 2^3S metastable,
   ! He(2^3S) + e^- -> He^+ + 2e^-, Black (1981, MNRAS 197, 553, eq. 12,
   ! READ from the published article):
   !   k(T) = 8.38e-10 T^1/2 exp(-55338/T)   [cm^3 s^-1],
   ! "based upon the cross sections of Taylor, Kingston & Bell (1979)" and
   ! "applicable in the range 1e3 <= T <= 1e5 K" (his words). The code
   ! formerly divided Black's Table 3 cooling coefficient by the potential
   ! (6.41e-21/7.6388e-12 = 8.391e-10, 0.13% above his printed rate). The
   ! 55338 K of the exponent is his fit's threshold (4.7687 eV, against
   ! 55328 K for the potential e_th_HeTR). The energy removed per
   ! ionization is the potential e_th_HeTR_erg (radiative_cooling_of_cell).
   ! Evaluated unclamped outside 1e3-1e5 K, where it carries no accuracy
   ! claim. T in K.
   elemental double precision function ionization_rate_HeI_23S(T)
   real*8, intent(in) :: T
   ionization_rate_HeI_23S = 8.38d-10*sqrt(T)*exp(-55338.0d0/T)
   end function ionization_rate_HeI_23S

   subroutine ci_HeI23S_range(T,coeff_ci_HeI23S,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: coeff_ci_HeI23S
   coeff_ci_HeI23S(j_lo:j_hi) = ionization_rate_HeI_23S(T(j_lo:j_hi))
   end subroutine ci_HeI23S_range



   !---------------------------------------------------!

   !--- Metal cooling (ported from AIOLOS photochem.cpp:95-127) ---!
   ! Functional form: Lambda(T) [erg cm^3 / s], per ion per electron.
   ! Total volumetric cooling = n_ion * n_e * Lambda(T).
   ! Provenance: introduced in AIOLOS commit 63695ea (Schulik, 2022-03-09);
   ! no in-code citation. Likely fit to Cloudy or similar tables; see
   ! Schulik & Booth (2023) for context. C++ (Cpp), O++ (Opp), H3+
   ! cooling are zero in AIOLOS and are intentionally omitted here.

   !--- CHIANTI C/N/O line cooling (DEFAULT; metals.inp 'cno_cool 0' ---!
   !--- reverts to the legacy AIOLOS fits) ------------------------------!
   ! The coronal (low-density) CHIANTI v11.0.2 line-cooling curves of
   ! C I, C II, N I, N II, O I and O II, Lambda(T) = T^-1/2 sum_i A_i
   ! exp(-T_i/T) per (n_e n_ion) with the ground-term fine structure
   ! Boltzmann-populated (fit script cooling_data/fit_cno_formulas.py, max
   ! error 0.15-4.6% over 1e3-1e5 K), are NOT what the cooling assembly
   ! uses. Their ground-term fine-structure lines ([C I] 609/370um,
   ! [C II] 158um, [N II] 205/122um, [O I] 63/145/44um) have critical
   ! densities of 1e0-1e5 cm^-3 and their metastable upper levels (2p2
   ! 1D2/1S0 for C I and N II, 2p4 1D2/1S0 for O I, 2p3 2D*/2P* for N I and
   ! O II, 2s2p2 4P for C II) of 1e4-1e9 cm^-3, at or below the electron
   ! density of a wind base, where the coronal form overestimates them by
   ! up to seven decades. The coefficients below (cool_CI_ne_func ...
   ! cool_OII_ne_func) instead solve the ground term exactly and read the
   ! channels above it from the density-resolved CHIANTI table
   ! (metal_cooling_above_ground_term). All CHIANTI-derived coefficients
   ! carry the coronal_excitation_cutoff factor, which removes them smoothly
   ! below the 1e3 K floor of the data they were fitted to (see that
   ! function).

   !--- Ground-term fine-structure statistical equilibrium -------------!
   ! The coronal CHIANTI curves carry the ground-term fine-structure
   ! (FS) lines in the optically thin, LOW-DENSITY limit. Their
   ! critical densities are of order 1e0-1e5 cm^-3, decades below the base
   ! density of an irradiated atmosphere (n_e ~ 1e9, n_HI ~ 1e14 cm^-3), so
   ! the coronal form overestimates the FS cooling there by up to seven
   ! decades. For every C/N/O coolant whose ground term is split, the FS
   ! part is therefore replaced by the EXACT statistical-equilibrium (SE)
   ! solution of the ground term,
   !   solve  sum_j f_j R_ji = f_i sum_j R_ij,   sum_i f_i = 1,
   !   R_ul = C_ul + beta_ul A_ul,  R_lu = C_ul (g_u/g_l) exp(-E_ul/kT),
   !   C_ul = n_e k_e,ul(T) + n_HI k_H,ul(T),
   !   W_FS = sum_{u>l} f_u beta_ul A_ul k_B E_ul       [erg/s per ion],
   ! and the coronal curve is replaced by a refit (Lambda_rem) that keeps
   ! only the channels LEAVING the ground term. The coefficient returned is
   !   Lambda_eff = W_FS/n_e + Lambda_rem      [erg cm^3 s^-1],
   ! so the assembly prefactor n_e*n_ion recovers the H-collision part
   ! exactly (the n_e cancels); n_e is floored to avoid 0/0.
   !
   ! IONS TREATED -- this is the complete set:
   !   C I  2p2 3P_0,1,2    [C I]  609.1 / 370.4 um
   !   C II 2p  2P_1/2,3/2  [C II] 157.7 um
   !   N II 2p2 3P_0,1,2    [N II] 205.3 / 121.8 um
   !   O I  2p4 3P_2,1,0    [O I]   63.2 / 145.5 / 44.1 um (inverted term)
   ! N I and O II have a single-level 4S ground term; Mg I/II, Ca II and
   ! Na I have a single ground level; Fe I is built from permitted lines
   ! only (no forbidden a5D fine structure in its table); and Fe II is
   ! already density-dependent through its 2-D SE table (cool_FeII_ne).
   !
   ! LIMITS. n_e, n_HI -> 0 collapses the SE populations onto the ground
   ! level and reproduces the coronal within-term channel to the accuracy
   ! of the Upsilon fits (1.9-7.8%). n_e or n_HI -> infinity saturates at
   ! the exact multilevel LTE emission sum_u f_u^Boltz A_ul h nu_ul (checked
   ! to 1e-14 relative). The earlier two-level treatment of [O I] 63um and
   ! [C II] 158um multiplied a two-level solution -- already normalized to
   ! the TWO-LEVEL partition function -- by the ground-term Boltzmann
   ! fraction of the lower level; that double normalization left its LTE
   ! limit low by a factor (1 + (g_u/g_l) exp(-E/kT)), i.e. 2.7x for
   ! [C II] 158um at the 530 K base of HD 189733 b.
   !
   ! Lambda_rem keeps the ground-term Boltzmann weighting of its lower
   ! levels (the convention of the CHIANTI curves it is fitted to). That is
   ! the right weighting wherever it matters: the FS levels sit above their
   ! critical densities, hence in LTE, throughout these winds. Where they
   ! do not, the ground-only and Boltzmann weightings agree to 1.4-1.9%
   ! anyway, because Upsilon_lu is roughly proportional to g_l for the
   ! LS-coupled transitions that leave the ground term.
   !
   ! LINE TRAPPING. beta_ul is the escape probability of that FS line and
   ! enters as A_ul -> beta_ul A_ul INSIDE the SE solution, not as a factor
   ! on the result. That is the physically correct place: in the
   ! subcritical limit the cooling is set by the collisional excitation
   ! rate and is independent of beta (every excitation still ends as an
   ! escaped photon, just later), while in the saturated limit the escaping
   ! flux is proportional to beta. beta = 1 reproduces the optically thin
   ! result exactly. Lambda_rem is left optically thin: it collects
   ! transitions to higher terms whose lower levels are far less populated,
   ! so their opacity is orders of magnitude below the FS lines'.
   !
   ! ATOMIC DATA. cooling_data/fit_fs_saturation.py is the auditable source
   ! and prints every coefficient below.
   !   levels, A values, electron Upsilon: CHIANTI v11.0.2 elvlc/wgfa/scups
   !     (Upsilon descaled with the natural cubic spline of CHIANTI's
   !     DESCALE_SCUPS.PRO)
   !   k_H (H-atom de-excitation):
   !     C I, N II  Yan & Babb (2023), MNRAS 518, 6004, Tables 1 and 2,
   !                relaxation rates (tabulated 10-1e4 K; fitted and clamped
   !                to 50-1e4 K), every value checked against the published
   !                article
   !     O I        Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171
   !                (20-1000 K, as tabulated in LAMDA oatom.dat)
   !     C II       Barinovs, van Hemert, Krems & Dalgarno (2005),
   !                ApJ 620, 537 (20-2000 K, LAMDA c+.dat)
   ! These replace the k_H values of the earlier two-level branch (a
   ! constant 4.0e-11 for [C II] 158um + H and 4.2e-11 (T/100)^0.67 for
   ! [O I] 63um + H), which were flagged approximate and are 5-20x below
   ! the quantum-scattering results above. H collisions on the ION N II are
   ! included on the same footing as on the neutrals; they are not what
   ! saturates N II at the base (n_crit,e of [N II] 205um is ~1e2 cm^-3),
   ! they matter only in gas thin enough for the electron channel alone to
   ! be subcritical.

   ! Effective collision strength of one transition as a fit to the CHIANTI
   ! .scups values descaled with the natural cubic spline of CHIANTI's
   ! DESCALE_SCUPS.PRO: log10(Upsilon) a quartic in x = log10(T/1e4),
   ! fitted over 1e3-1e5 K (the ground-term fine-structure strengths below
   ! to 8.2% max, reached at the 1e3 K end on C I, where a cubic
   ! reaches 9.0% on the C I 3P0-3P2 strength; H I 1s-2s and 1s-2p
   ! to 2.2% and 0.8%). T is clamped to the fitted range, so
   ! outside it the strength is held at the endpoint value rather than
   ! extrapolated by a quartic.
   pure double precision function upsilon_quartic_logT                 &
                                    (c0,c1,c2,c3,c4,T)
   real*8, intent(in) :: c0,c1,c2,c3,c4,T
   real*8 :: x
   x = log10(min(max(T,1.0d3),1.0d5)/1.0d4)
   upsilon_quartic_logT = 10.0d0**(c0 + x*(c1 + x*(c2 + x*(c3 + x*c4))))
   end function upsilon_quartic_logT

   ! Electron-impact de-excitation rate coefficient [cm^3 s^-1] of a
   ! transition with upper-level weight g_u and effective collision
   ! strength Upsilon: k_ul = coll_rate_prefactor Upsilon / (g_u sqrt(T)).
   pure double precision function electron_impact_deexcitation(ups,g_u,T)
   real*8, intent(in) :: ups, g_u, T
   electron_impact_deexcitation = coll_rate_prefactor*ups/(g_u*sqrt(max(T,1.0d0)))
   end function electron_impact_deexcitation

   ! H-atom impact de-excitation rate coefficient [cm^3 s^-1], log10(k_H) a
   ! quadratic in u = log10(T/1e3) fitted to the tabulated quantum results
   ! at T >= 50 K (max error 4.9%, O I 3->2; 3.2% for the Yan & Babb
   ! rates). T is clamped to [Tlo,Thi], the fitted part of the tabulated
   ! range (C I, N II 50-1e4 K; O I 50-1e3 K; C II 60-2e3 K), so the rate
   ! is held at its endpoint value outside it.
   pure double precision function h_impact_deexcitation                  &
                                    (a0,a1,a2,Tlo,Thi,T)
   real*8, intent(in) :: a0,a1,a2,Tlo,Thi,T
   real*8 :: u
   u = log10(min(max(T,Tlo),Thi)/1.0d3)
   h_impact_deexcitation = 10.0d0**(a0 + u*(a1 + u*a2))
   end function h_impact_deexcitation

   ! Exact statistical equilibrium of a TWO-level ground term.
   ! f_2 = R_12/(R_12 + R_21); returns the NET escaping power per ion
   ! [erg/s], emission minus absorption of the incident field.
   ! nb21 is the photon occupation number of that incident field at the line
   ! frequency (0 = no incident radiation, the optically thin vacuum limit),
   ! so the radiative rates are b A (1 + nb) down and b A (g2/g1) nb up, and
   !   W = kB E b A [ f2 (1 + nb) - f1 (g2/g1) nb ].
   ! With nb = W_dil/(exp(E/T_rad) - 1) this vanishes when the level ratio
   ! reaches the equilibrium set by the field, which is the radiative
   ! equilibrium floor of the line; below it W is negative, i.e. the line
   ! heats the gas.
   pure double precision function fine_structure_cooling_2level          &
                                    (T,g1,g2,E2,A21,C21,b21,nb21) result(W)
   real*8, intent(in) :: T,g1,g2,E2,A21,C21,b21,nb21
   ! The emission factor is written first and multiplied by (1 + nb21) so
   ! that nb21 = 0 reproduces the vacuum expression bit for bit.
   real*8 :: R12, R21, f2
   R21 = C21 + b21*A21*(1.0d0 + nb21)
   R12 = C21*(g2/g1)*exp(-E2/T) + b21*A21*(g2/g1)*nb21
   f2  = R12/max(R12 + R21, 1.0d-300)
   W   = kb_erg*E2*b21*A21*R12/max(R12 + R21, 1.0d-300)*(1.0d0 + nb21) &
       - kb_erg*E2*b21*A21*(1.0d0 - f2)*(g2/g1)*nb21
   end function fine_structure_cooling_2level

   ! Exact statistical equilibrium of a THREE-level ground term. Levels are
   ! in energy order with E1 = 0; C_ul are the total collisional
   ! de-excitation rates [s^-1], b_ul the escape probabilities and nb_ul the
   ! photon occupation numbers of the incident field at each line frequency
   ! (0 = vacuum). The radiative rates are b A (1 + nb) down and
   ! b A (g_u/g_l) nb up (see the two-level routine above).
   ! Eliminating f_1 = 1 - f_2 - f_3 from the two level-balance equations
   ! leaves a 2x2 system, solved by Cramer's rule:
   !   f_2 (R12+R21+R23) + f_3 (R12-R32) = R12
   !   f_2 (R13-R23) + f_3 (R13+R31+R32) = R13
   ! This is the ONE definition of that equilibrium in the code. The cooling
   ! power below is assembled from what it returns, and the exported O I
   ! level field (oxygen_ground_term_populations) calls it with the same
   ! rates, so the populations a transit forward model reads are by
   ! construction the populations the cooling was computed from.
   ! f1 + f2 + f3 = 1 identically.
   pure subroutine fine_structure_populations_3level                     &
                                    (T,g1,g2,g3,E2,E3,A21,A31,A32,      &
                                     C21,C31,C32,b21,b31,b32,           &
                                     nb21,nb31,nb32, f1,f2,f3)
   real*8, intent(in)  :: T,g1,g2,g3,E2,E3,A21,A31,A32
   real*8, intent(in)  :: C21,C31,C32,b21,b31,b32,nb21,nb31,nb32
   real*8, intent(out) :: f1,f2,f3
   real*8 :: R12,R13,R23,R21,R31,R32, m11,m12,m21,m22, det
   R21 = C21 + b21*A21*(1.0d0 + nb21)
   R31 = C31 + b31*A31*(1.0d0 + nb31)
   R32 = C32 + b32*A32*(1.0d0 + nb32)
   R12 = C21*(g2/g1)*exp(-E2/T) + b21*A21*(g2/g1)*nb21
   R13 = C31*(g3/g1)*exp(-E3/T) + b31*A31*(g3/g1)*nb31
   R23 = C32*(g3/g2)*exp(-(E3 - E2)/T) + b32*A32*(g3/g2)*nb32
   m11 = R12 + R21 + R23
   m12 = R12 - R32
   m21 = R13 - R23
   m22 = R13 + R31 + R32
   det = m11*m22 - m12*m21
   ! det > 0 for any non-degenerate set of positive rates; the guard only
   ! catches the fully depopulated limit (every rate zero), where the
   ! excited fractions are zero too.
   if (det .gt. 0.0d0) then
      f2 = (R12*m22 - m12*R13)/det
      f3 = (m11*R13 - m21*R12)/det
   else
      f2 = 0.0d0
      f3 = 0.0d0
   endif
   f1 = 1.0d0 - f2 - f3
   end subroutine fine_structure_populations_3level

   ! NET escaping power per ion [erg/s] of the same three-level term:
   ! emission minus absorption of the incident field.
   pure double precision function fine_structure_cooling_3level          &
                                    (T,g1,g2,g3,E2,E3,A21,A31,A32,      &
                                     C21,C31,C32,b21,b31,b32,           &
                                     nb21,nb31,nb32) result(W)
   real*8, intent(in) :: T,g1,g2,g3,E2,E3,A21,A31,A32
   real*8, intent(in) :: C21,C31,C32,b21,b31,b32,nb21,nb31,nb32
   real*8 :: f1,f2,f3
   call fine_structure_populations_3level(T,g1,g2,g3,E2,E3,A21,A31,A32, &
                                          C21,C31,C32,b21,b31,b32,      &
                                          nb21,nb31,nb32, f1,f2,f3)
   ! Each line contributes f_u b A (1 + nb) - f_l b A (g_u/g_l) nb; the
   ! emission terms are written exactly as in the vacuum expression so that
   ! nb = 0 reproduces it bit for bit.
   W = kb_erg*( f2*b21*A21*E2*(1.0d0 + nb21)                         &
              + f3*b31*A31*E3*(1.0d0 + nb31)                         &
              + f3*b32*A32*(E3 - E2)*(1.0d0 + nb32) )                &
     - kb_erg*( f1*b21*A21*E2*(g2/g1)*nb21                           &
              + f1*b31*A31*E3*(g3/g1)*nb31                           &
              + f2*b32*A32*(E3 - E2)*(g3/g2)*nb32 )
   end function fine_structure_cooling_3level

   ! C I: [C I] 609.1um (3P1-3P0) and 370.4um (3P2-3P1). The 3P2-3P0
   ! channel has no transition probability but does couple the levels
   ! collisionally, so it enters C31 with A31 = 0.
   elemental double precision function cool_CI_ne_func                   &
                                         (T,ne,nHI,b609,b370,n609,n370)
   real*8, intent(in) :: T, ne, nHI, b609, b370, n609, n370
   real*8 :: Ts, C21, C31, C32, W, lam_rem
   Ts  = max(T, 1.0d0)
   C21 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -0.36997513d0, 0.2677744d0, -0.44730464d0,                 &
             0.32574003d0, -0.066939801d0, Ts), 3.0d0, Ts)              &
       + nHI*h_impact_deexcitation(-9.6487495d0, 0.27178113d0,          &
             0.110106d0, 5.0d1, 1.0d4, Ts)
   C31 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -0.58761698d0, 0.12443835d0, -0.29369201d0,                  &
             0.39395787d0, -0.13898642d0, Ts), 5.0d0, Ts)               &
       + nHI*h_impact_deexcitation(-9.7250154d0, 0.35123928d0,          &
             0.032556919d0, 5.0d1, 1.0d4, Ts)
   C32 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
             0.04798486d0, 0.19117637d0, -0.35985595d0,                   &
             0.36015847d0, -0.10999838d0, Ts), 5.0d0, Ts)                &
       + nHI*h_impact_deexcitation(-9.2461144d0, 0.34543127d0,          &
             0.053368443d0, 5.0d1, 1.0d4, Ts)
   W = fine_structure_cooling_3level(Ts, 1.0d0, 3.0d0, 5.0d0,           &
          Ek_CI2, Ek_CI3, A_CI609, 0.0d0, A_CI370,                      &
          C21, C31, C32, b609, 1.0d0, b370, n609, 0.0d0, n370)
   ! The channels that LEAVE the ground term. lam_rem is their coronal
   ! limit, the reference the density-resolved table replaces and the
   ! value it falls back to for an ion with no table.
   lam_rem = ( 1.63414350d-18*exp(-16371.1d0/Ts)                        &
             + 4.47668751d-18*exp(-22644.4d0/Ts)                        &
             + 4.27846098d-17*exp(-59981.7d0/Ts)                        &
             + 2.88107272d-16*exp(-156387.0d0/Ts) )/sqrt(Ts)
   cool_CI_ne_func = W/max(ne, 1.0d-30)                                 &
        + metal_cooling_above_ground_term(im_CI, Ts, ne, lam_rem) &
                 *coronal_excitation_cutoff(Ts)
   end function cool_CI_ne_func

   ! C II: [C II] 157.7um (2P3/2-2P1/2). The ground term has only two
   ! levels, so the two-level solution is exact.
   elemental double precision function cool_CII_ne_func(T,ne,nHI,b158,n158)
   real*8, intent(in) :: T, ne, nHI, b158, n158
   real*8 :: Ts, C21, W, lam_rem
   Ts  = max(T, 1.0d0)
   C21 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
             0.34358283d0, 0.11624658d0, -0.17963499d0,                 &
            -0.061985731d0, 0.10606931d0, Ts), 4.0d0, Ts)               &
       + nHI*h_impact_deexcitation(-8.9730551d0, 0.19760195d0,          &
             0.048910789d0, 6.0d1, 2.0d3, Ts)
   W = fine_structure_cooling_2level(Ts, 2.0d0, 4.0d0, Ek_CII2,         &
                                     A_CII158, C21, b158, n158)
   lam_rem = ( 2.62718417d-17*exp(-61536.5d0/Ts)                        &
             + 1.07377677d-17*exp(-69811.7d0/Ts)                        &
             + 3.29228566d-16*exp(-123223.0d0/Ts)                       &
             + 1.39114570d-15*exp(-247558.0d0/Ts) )/sqrt(Ts)
   cool_CII_ne_func = W/max(ne, 1.0d-30)                                &
        + metal_cooling_above_ground_term(im_CII, Ts, ne, lam_rem) &
                 *coronal_excitation_cutoff(Ts)
   end function cool_CII_ne_func

   ! N II: [N II] 205.3um (3P1-3P0) and 121.8um (3P2-3P1); the 3P2-3P0
   ! channel again couples collisionally only.
   elemental double precision function cool_NII_ne_func                  &
                                         (T,ne,nHI,b205,b122,n205,n122)
   real*8, intent(in) :: T, ne, nHI, b205, b122, n205, n122
   real*8 :: Ts, C21, C31, C32, W, lam_rem
   Ts  = max(T, 1.0d0)
   C21 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -0.42450591d0, 0.07672821d0, 0.015109843d0,              &
            -0.038977811d0, 0.031565297d0, Ts), 3.0d0, Ts)              &
       + nHI*h_impact_deexcitation(-9.5141594d0, 0.11153196d0,          &
             0.11112557d0, 5.0d1, 1.0d4, Ts)
   C31 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -0.7278854d0, 0.40068378d0, 0.012785755d0,                 &
            -0.17640487d0, -0.0075761571d0, Ts), 5.0d0, Ts)             &
       + nHI*h_impact_deexcitation(-9.4220729d0, 0.29532393d0,          &
             0.022835614d0, 5.0d1, 1.0d4, Ts)
   C32 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
             0.16134799d0, 0.068443207d0, -0.068085615d0,               &
            -0.031838155d0, 0.0297447d0, Ts), 5.0d0, Ts)              &
       + nHI*h_impact_deexcitation(-8.990019d0, 0.28444121d0,           &
             0.0409254d0, 5.0d1, 1.0d4, Ts)
   W = fine_structure_cooling_3level(Ts, 1.0d0, 3.0d0, 5.0d0,           &
          Ek_NII2, Ek_NII3, A_NII205, 0.0d0, A_NII122,                  &
          C21, C31, C32, b205, 1.0d0, b122, n205, 0.0d0, n122)
   lam_rem = ( 7.71532491d-18*exp(-23578.5d0/Ts)                        &
             + 6.44871941d-18*exp(-49880.2d0/Ts)                        &
             + 1.40104558d-16*exp(-125780.0d0/Ts)                       &
             + 7.13236468d-16*exp(-247915.0d0/Ts) )/sqrt(Ts)
   cool_NII_ne_func = W/max(ne, 1.0d-30)                                &
        + metal_cooling_above_ground_term(im_NII, Ts, ne, lam_rem) &
                 *coronal_excitation_cutoff(Ts)
   end function cool_NII_ne_func

   ! O I: the term is inverted (3P2 lowest), so level 1 is 3P2 and the
   ! lines are [O I] 63.2um (3P1-3P2), 145.5um (3P0-3P1) and the very weak
   ! 44.1um (3P0-3P2).
   !
   ! Total collisional de-excitation rates [s^-1] of the three couplings,
   ! electrons plus H atoms. Written once here because two consumers need
   ! exactly these rates: the cooling coefficient below and the exported
   ! level populations (oxygen_ground_term_populations). T is expected
   ! already clamped away from zero; the clamp is repeated for safety and
   ! is idempotent.
   pure subroutine oxygen_ground_term_collisions(T,ne,nHI,C21,C31,C32)
   real*8, intent(in)  :: T, ne, nHI
   real*8, intent(out) :: C21, C31, C32
   real*8 :: Ts
   Ts  = max(T, 1.0d0)
   C21 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -2.0708554d0, 0.18660187d0, -0.22731732d0,                   &
             0.061602392d0, 0.049940423d0, Ts), 3.0d0, Ts)              &
       + nHI*h_impact_deexcitation(-9.0268259d0, 0.45769694d0,          &
             0.050610554d0, 5.0d1, 1.0d3, Ts)
   C31 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -2.4228892d0, 0.18267434d0, -0.22162813d0,                  &
             0.063243772d0, 0.045163335d0, Ts), 1.0d0, Ts)                &
       + nHI*h_impact_deexcitation(-9.1067415d0, 0.43724931d0,          &
             0.05693361d0, 5.0d1, 1.0d3, Ts)
   C32 = ne*electron_impact_deexcitation(upsilon_quartic_logT(        &
            -3.9241026d0, 0.27572369d0, -0.30036307d0,                  &
             0.13841574d0, -0.0033242401d0, Ts), 1.0d0, Ts)               &
       + nHI*h_impact_deexcitation(-8.9724741d0, 0.034669741d0,         &
            -0.35264943d0, 5.0d1, 1.0d3, Ts)
   end subroutine oxygen_ground_term_collisions

   ! Fractional populations of the three O I ground-term levels, in the
   ! same statistical equilibrium the cooling solves and with the same
   ! rates: f_3P2 + f_3P1 + f_3P0 = 1. These are the lower levels of the
   ! O I 1302.168 / 1304.858 / 1306.029 A resonance triplet respectively,
   ! which is what makes them an output field and not just an internal of
   ! the cooling: a transit forward model that applied the total O I
   ! density to all three components would count the same atoms three
   ! times.
   elemental subroutine oxygen_ground_term_populations                   &
                           (T,ne,nHI,b63,b145,b44,n63,n145,n44,          &
                            f_3P2,f_3P1,f_3P0)
   real*8, intent(in)  :: T, ne, nHI, b63, b145, b44, n63, n145, n44
   real*8, intent(out) :: f_3P2, f_3P1, f_3P0
   real*8 :: Ts, C21, C31, C32
   Ts = max(T, 1.0d0)
   call oxygen_ground_term_collisions(Ts,ne,nHI,C21,C31,C32)
   call fine_structure_populations_3level(Ts, 5.0d0, 3.0d0, 1.0d0,      &
          Ek_OI2, Ek_OI3, A_OI63, A_OI44, A_OI145,                      &
          C21, C31, C32, b63, b44, b145, n63, n44, n145,                &
          f_3P2, f_3P1, f_3P0)
   end subroutine oxygen_ground_term_populations

   elemental double precision function cool_OI_ne_func                   &
                                         (T,ne,nHI,b63,b145,b44,        &
                                          n63,n145,n44)
   real*8, intent(in) :: T, ne, nHI, b63, b145, b44, n63, n145, n44
   real*8 :: Ts, C21, C31, C32, W, lam_rem
   Ts  = max(T, 1.0d0)
   call oxygen_ground_term_collisions(Ts,ne,nHI,C21,C31,C32)
   W = fine_structure_cooling_3level(Ts, 5.0d0, 3.0d0, 1.0d0,           &
          Ek_OI2, Ek_OI3, A_OI63, A_OI44, A_OI145,                      &
          C21, C31, C32, b63, b44, b145, n63, n44, n145)
   lam_rem = ( 2.63734982d-19*exp(-23770.4d0/Ts)                        &
             + 1.18242203d-18*exp(-30965.4d0/Ts)                        &
             + 4.96806093d-18*exp(-58929.4d0/Ts)                        &
             + 3.12665954d-17*exp(-160409.0d0/Ts) )/sqrt(Ts)
   cool_OI_ne_func = W/max(ne, 1.0d-30)                                 &
        + metal_cooling_above_ground_term(im_OI, Ts, ne, lam_rem) &
                 *coronal_excitation_cutoff(Ts)
   end function cool_OI_ne_func

   ! N I and O II have a single-level 4S* ground term, so there is no
   ! ground-term fine structure to solve; what needs the electron density
   ! is the metastable manifold above it, 2p3 2D* and 2P* in both ions.
   ! These two coefficients therefore carry only the saturation; lam_rem
   ! is the coronal CHIANTI fit of the whole ion. The legacy (cno_chianti
   ! off) branches are selected by metal_line_cooling_coefficient and never
   ! reach here.

   ! N I: the whole fit is a channel leaving the 4S* ground level, so the
   ! suppression multiplies all of it.
   elemental double precision function cool_NI_ne_func(T,ne)
   real*8, intent(in) :: T, ne
   real*8 :: Ts, lam_rem
   Ts = max(T, 1.0d0)
   lam_rem = ( 1.51497279d-18*exp(-29686.4d0/Ts)                        &
             + 4.85982268d-18*exp(-34688.5d0/Ts)                        &
             + 1.14706763d-17*exp(-50580.8d0/Ts)                        &
             + 1.49777997d-16*exp(-138115.0d0/Ts)                       &
             + 4.45186229d-16*exp(-265801.0d0/Ts) )/sqrt(Ts)
   cool_NI_ne_func =                                                    &
        metal_cooling_above_ground_term(im_NI, Ts, ne, lam_rem) &
               *coronal_excitation_cutoff(Ts)
   end function cool_NI_ne_func

   ! O II: as for N I, the ground level 4S*3/2 is single and the whole fit
   ! leaves it.
   elemental double precision function cool_OII_ne_func(T,ne)
   real*8, intent(in) :: T, ne
   real*8 :: Ts, lam_rem
   Ts = max(T, 1.0d0)
   lam_rem = ( 1.76706465d-17*exp(-44127.3d0/Ts)                        &
             + 8.05195645d-18*exp(-59651.1d0/Ts)                        &
             + 2.32124773d-16*exp(-178708.0d0/Ts)                       &
             + 8.78052704d-16*exp(-324831.0d0/Ts) )/sqrt(Ts)
   cool_OII_ne_func =                                                   &
        metal_cooling_above_ground_term(im_OII, Ts, ne, lam_rem) &
               *coronal_excitation_cutoff(Ts)
   end function cool_OII_ne_func

   ! Ca II and Mg I: the single ground level (4s 2S1/2, 3s2 1S0) is the
   ! whole ground term, so the density-resolved table carries every channel
   ! of the ion at the local electron density. The coronal resonance-line
   ! skeleton is the argument metal_cooling_above_ground_term returns for
   ! an ion with no table.
   elemental double precision function cool_CaII_ne_func(T,ne)
   real*8, intent(in) :: T, ne
   real*8 :: Ts
   Ts = max(T, 1.0d0)
   cool_CaII_ne_func =                                                  &
        metal_cooling_above_ground_term(im_CaII, Ts, ne,                &
                                        cool_CaII_HK_coronal(Ts))       &
               *coronal_excitation_cutoff(Ts)
   end function cool_CaII_ne_func

   elemental double precision function cool_MgI_ne_func(T,ne)
   real*8, intent(in) :: T, ne
   real*8 :: Ts
   Ts = max(T, 1.0d0)
   cool_MgI_ne_func =                                                   &
        metal_cooling_above_ground_term(im_MgI, Ts, ne,                 &
                                        cool_MgI_2853_coronal(Ts))      &
               *coronal_excitation_cutoff(Ts)
   end function cool_MgI_ne_func



   ! Grid version of the exported O I level fractions, fed by the same
   ! beta_fs / nbar_fs that fine_structure_line_transfer produced for the
   ! cooling (cool_OI_ne_func).
   subroutine oxygen_ground_term_levels(T,ne,nHI,beta_fs,nbar_fs,       &
                                        f_3P2,f_3P1,f_3P0)
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T, ne, nHI
   real*8, dimension(1-Ng:N+Ng,n_fsline), intent(in) :: beta_fs, nbar_fs
   real*8, dimension(1-Ng:N+Ng), intent(out) :: f_3P2, f_3P1, f_3P0
   call oxygen_ground_term_populations(T,ne,nHI,beta_fs(:,ifs_OI63),    &
                        beta_fs(:,ifs_OI145),beta_fs(:,ifs_OI44),       &
                        nbar_fs(:,ifs_OI63),nbar_fs(:,ifs_OI145),       &
                        nbar_fs(:,ifs_OI44),                            &
                        f_3P2,f_3P1,f_3P0)
   end subroutine oxygen_ground_term_levels

   !--------------!



   !--------------!


   ! Value of a CHIANTI cooling table at one temperature -- the ONE definition
   ! of this interpolation (metal_line_cooling_coefficient, for Fe I).
   !
   ! cool_logT is uniform (dlogT, 3.0..5.0 in log10 K); outside the table the
   ! nearest endpoint is held. Returns Lambda(T) per (n_e * n_ion)
   ! [erg cm^3 s^-1].
   !
   ! Monotonicity-preserving piecewise-cubic Hermite (PCHIP) on log10(Lambda)
   ! vs log10(T), not plain linear interpolation. Linear interpolation is only
   ! C0: dLambda/dT jumps at every table node, and the semi-implicit energy
   ! solver finite-differences the cooling to get dC/dT -- so those slope jumps
   ! inject non-smooth derivatives that can destabilize the energy update for
   ! the metal-rich wind (the analytic C/N/O cooling is C-infinity and does
   ! not). PCHIP is C1 (continuous slope) and introduces no new
   ! extrema/overshoot, so it removes the interpolation kinks while staying
   ! monotone where the table is.
   pure double precision function cool_table_value(logL,Ts)
   real*8, dimension(NCOOLT), intent(in) :: logL
   real*8, intent(in) :: Ts
   integer :: k
   real*8  :: lt, pos, frac, m0, m1, h00, h10, h01, h11
   lt  = log10(max(Ts, 1.0d0))
   pos = (lt - cool_logT(1))/cool_dlogT + 1.0d0
   ! NB: the first branch is written .not.(pos>1) rather than (pos<=1) so a
   ! non-finite pos (transient NaN temperature during relaxation, or a NaN
   ! trial T from the post-process root finder) lands on the table edge
   ! instead of reaching int(pos) -> huge negative index -> out-of-bounds
   ! read (caught by -fcheck=bounds on WASP-121b).
   if (.not. (pos .gt. 1.0d0)) then
      cool_table_value = 10.0d0**logL(1)
   else if (pos .ge. dble(NCOOLT)) then
      cool_table_value = 10.0d0**logL(NCOOLT)
   else
      k    = int(pos)
      frac = pos - dble(k)                 ! in [0,1), interval [k, k+1]
      ! PCHIP node slopes (d log10(Lambda) per unit index; grid spacing = 1).
      m0 = pchip_slope(logL, k)
      m1 = pchip_slope(logL, k+1)
      ! cubic Hermite basis on the unit interval
      h00 = (1.0d0 + 2.0d0*frac)*(1.0d0 - frac)**2
      h10 = frac*(1.0d0 - frac)**2
      h01 = frac**2*(3.0d0 - 2.0d0*frac)
      h11 = frac**2*(frac - 1.0d0)
      cool_table_value = 10.0d0**( h00*logL(k)   + h10*m0          &
                                 + h01*logL(k+1) + h11*m1 )
   endif
   ! Below the 1e3 K table edge the value above is the held endpoint,
   ! which is an unconstrained clamp rather than a cooling rate.
   cool_table_value = cool_table_value*coronal_excitation_cutoff(Ts)
   end function cool_table_value


   ! Monotonicity-preserving (PCHIP / Fritsch-Carlson) node slope at index i,
   ! in units of d(logL) per unit index (the grid is uniform in index). At a
   ! local extremum (secants of opposite sign) the slope is set to zero so the
   ! cubic neither overshoots nor introduces a new extremum.
   pure real*8 function pchip_slope(L, i)
   real*8, dimension(NCOOLT), intent(in) :: L
   integer, intent(in) :: i
   real*8 :: sL, sR
   if (i .le. 1) then
      pchip_slope = L(2) - L(1)
   else if (i .ge. NCOOLT) then
      pchip_slope = L(NCOOLT) - L(NCOOLT-1)
   else
      sL = L(i)   - L(i-1)
      sR = L(i+1) - L(i)
      if (sL*sR .le. 0.0d0) then
         pchip_slope = 0.0d0
      else
         pchip_slope = 2.0d0/(1.0d0/sL + 1.0d0/sR)   ! weighted harmonic mean
      endif
   endif
   end function pchip_slope

   !--------------!



   !--------------!

   ! Closed-form analytic fits to the CHIANTI v11 cooling curves
   ! (cooling_data/fit_cooling_formulas.py; cf. metal_cooling_chianti.txt).
   ! Resonance lines use the exact two-level skeleton
   !    Lambda(T) = coll_rate_prefactor/(g_l sqrt(T)) * Ups(T) * dE * exp(-dE/kT)
   ! with the Burgess-Tully type-1 effective collision strength
   !    Ups(T) = a + b*ln(1 + T/T0)
   ! fitted to the CHIANTI table (collision strengths descaled with the
   ! natural cubic spline of CHIANTI's DESCALE_SCUPS.PRO) with dE free.
   ! For Mg I and Ca II the CHIANTI .scups energies equal the observed
   ! ones to 0.02 per cent, and the fitted dE lies 0.2 per cent (Mg I) and
   ! 0.05 per cent (Ca II) from the observed transition energy. Max
   ! error vs the 201-point source table over 1e3-1e5 K: Mg I 2.8%, Ca II
   ! 1.3%, Na I 0.01% (Na I is
   ! exactly Van Regemorter with constant Ups by construction). Mg II has
   ! its own fit with the observed energies (below).

   ! Mg I 2853 A (3s2 1S0 -> 3s3p 1P1); dE = 4.3381 eV, g_l = 1. Resonance
   ! line alone, coronal limit, fitted over 1e3-1e5 K, no low-T cutoff:
   ! the skeleton cool_MgI_ne_func passes for an ion with no table. The
   ! cooling uses the density-resolved table cool_logLrem_MgI instead.
   elemental double precision function cool_MgI_2853_coronal(T)
   real*8, intent(in) :: T
   real*8 :: ups
   ups = 0.35452071d0 + 15.099537d0*log(1.0d0 + T/93933.324d0)
   cool_MgI_2853_coronal = coll_rate_prefactor/sqrt(T)*ups*6.94992895d-12      &
                           *exp(-50338.130d0/T)
   end function cool_MgI_2853_coronal

   ! Mg II line cooling per (n_e n_MgII): every electron-impact excitation
   ! out of the 3s 2S1/2 ground level into the 32 bound levels of the
   ! CHIANTI v11.0.2 model ion (its levels above the ionization energy
   ! autoionize), in the coronal limit, with the OBSERVED level energies in
   ! the Boltzmann factor and in the radiated energy (the .scups header
   ! energy of 3s-3p, 4.27 eV, is the theoretical one and belongs only to
   ! the descaling of the collision strength). h&k 2796/2803 A (3s-3p):
   ! dE = 4.4300 eV, the g-weighted observed 3p 2P energy, g_l = 2,
   ! Ups(T) = a + b ln(1 + T/T0) fitted (collision strengths descaled with
   ! the natural cubic spline of CHIANTI's DESCALE_SCUPS.PRO through the
   ! 15 stored points); one exponential for the doublet
   ! costs 0.2 per cent at 1e3 K. The excitations above 3p (4s 8.65 eV,
   ! 3d 8.86 eV, 4p and higher) are one effective term
   ! c exp(-T_x/T)/sqrt(T): 0.3 per cent of the total at 1e4 K, 4.3 per
   ! cent at 2e4 K, 28 per cent at 1e5 K. Fit error 0.2 per cent max
   ! over 1e3-1e5 K (cooling_data/magnesium_ii_line_cooling.py).
   ! Validity: optically thin; n_e <= 1e12 cm^-3, where the full
   ! statistical equilibrium of the bound levels equals this coronal sum to
   ! 0.2 per cent (2 per cent at 1e13, 17 per cent at 1e14); every excited
   ! bound level decays by a permitted line.
   ! The dominant metal-line coolant of the upper thermosphere in
   ! ultrahot Jupiters (Huang+2023 Fig. 4).
   elemental double precision function cool_MgII_func(T)
   real*8, intent(in) :: T
   real*8 :: ups
   ups = 15.789567d0 + 46.028142d0*log(1.0d0 + T/208367.52d0)
   cool_MgII_func = ( coll_rate_prefactor/(2.0d0*sqrt(T))*ups*7.09764369d-12     &
                      *exp(-51408.024d0/T)                            &
                    + 4.17321676d-16*exp(-132045.745d0/T)/sqrt(T) )     &
                    *coronal_excitation_cutoff(T)
   end function cool_MgII_func

   ! Ca II H&K 3934/3969 A (4s 2S -> 4p 2P); dE = 3.1438 eV, g_l = 2.
   ! Resonance line alone, coronal limit, fitted over 1e3-1e5 K, no low-T
   ! cutoff: the skeleton cool_CaII_ne_func passes for an ion with no
   ! table. The cooling uses the density-resolved table cool_logLrem_CaII
   ! instead.
   elemental double precision function cool_CaII_HK_coronal(T)
   real*8, intent(in) :: T
   real*8 :: ups
   ups = 14.653524d0 + 18.705877d0*log(1.0d0 + T/33882.599d0)
   cool_CaII_HK_coronal = coll_rate_prefactor/(2.0d0*sqrt(T))*ups*5.03631038d-12 &
                          *exp(-36477.848d0/T)
   end function cool_CaII_HK_coronal

   ! Na I D 5890/5896 A (3s 2S -> 3p 2P); dE = 2.1037 eV, g_l = 2,
   ! Ups = 36.07 (constant; rigorous Van Regemorter, gbar = 0.2)
   elemental double precision function cool_NaI_func(T)
   real*8, intent(in) :: T
   cool_NaI_func = coll_rate_prefactor/(2.0d0*sqrt(T))*36.074352d0             &
                   *3.37049304d-12*exp(-24412.382d0/T)               &
                   *coronal_excitation_cutoff(T)
   end function cool_NaI_func




   !---------------------------------------------------!

   !--- Validity floor of the coronal line-cooling fits ----------------!
   ! Every metal line-cooling coefficient in this module is a fit to, or
   ! an interpolation of, a CHIANTI optically thin curve tabulated over
   ! 1e3-1e5 K: the analytic C/N/O and Mg/Ca/Na/Fe forms, and the 1-D /
   ! 2-D tables (whose log10(T) axis starts exactly at 3.0 and which
   ! otherwise hold their edge value indefinitely below it). Below that
   ! floor the fitted exponentials are unconstrained by the data behind
   ! them, and the softest one dominates: for O I it is exp(-930.111/T),
   ! whose 930 K matches NO [O I] ground-term splitting (the splittings
   ! are 227.7 K and 326.6 K); for C I it is exp(-2351.38/T) against
   ! splittings of 23.6 K and 62.4 K. Evaluated at the ~240 K base of
   ! HD 189733 b those components supplied 99.7% of the O I cooling rate
   ! and drove the base to a fifth of T_eq.
   !
   ! Physically, at T far below 1e3 K the only metal transitions still
   ! collisionally excitable are the ground-term fine-structure lines,
   ! and those are carried EXPLICITLY by the statistical-equilibrium
   ! solutions (cool_CI_ne_func, cool_CII_ne_func, cool_NII_ne_func,
   ! cool_OI_ne_func), which also include the critical-density saturation
   ! the coronal curves lack. The coronal part is therefore switched off
   ! below the fit floor.
   !
   ! The switch is a Gaussian in the fractional temperature deficit,
   !   x = (T_floor/T - 1)/w,   cutoff = exp(-x^2)   for T < T_floor,
   ! and exactly 1 for T >= T_floor. Both the value and the dT-derivative
   ! are continuous at T_floor (x = 0 and dx/dT finite, so d(exp(-x^2))/dT
   ! = 0 there), which matters because the Brent energy solve and the
   ! semi-implicit update finite-difference the cooling in T. At and above
   ! 1e3 K every coefficient is BIT-IDENTICAL to the unguarded form.
   !
   ! VALIDITY / CHOICE: w = 0.1 is the default ("Coronal cutoff width: <w>"
   ! in input.inp changes it). w is a modeling choice, not a measured
   ! quantity, and it was bounded to 0.08-0.13 in
   ! md/coronal_cutoff_width.md by (i) confining the band in which the
   ! unsaturated coronal residual still dominated and (ii) keeping the
   ! guard's own log slope within an order of magnitude of the 8.3-8.7 the
   ! cooling function already has inside the fitted range.
   !
   ! WHAT THE GUARD NOW DOES. That bound was set while the ground-term
   ! fine-structure floors of C I, N II and the [O I] 145.5um channel were
   ! still coronal, and they dominated everything the guard removed. Now
   ! that all four split ground terms are solved in statistical
   ! equilibrium, no coronal term reaching below 1e3 K is left: the softest
   ! exponential surviving in any of the remainders is exp(-16400/T), which
   ! at 500 K is 1e-15 of its 1e4 K value. The guard is back to its stated
   ! job -- refusing to extrapolate a fit below the data it was made from --
   ! and the base is no longer sensitive to w: the local balance
   ! temperature of the base cell is identical to six digits over
   ! w = 0.02-1.2 for all four planet runs, against a 500-650 K spread
   ! before (md/coronal_cutoff_width.md section 7).
   ! The legacy AIOLOS branch (cno_chianti = .false.) is deliberately NOT
   ! guarded: its constant floors (1.0e-24 etc.) are crude stand-ins for
   ! fine-structure cooling, not extrapolated coronal fits.
   elemental double precision function coronal_excitation_cutoff(T)
   real*8, intent(in) :: T
   real*8, parameter :: T_fit_floor = 1.0d3   ! [K] lower edge of the fits
   real*8 :: x
   ! .not.(T < floor) so a non-finite T falls on the unguarded branch,
   ! matching how the table interpolators handle a transient NaN.
   if (.not. (T .lt. T_fit_floor)) then
      coronal_excitation_cutoff = 1.0d0
   else
      ! coronal_cutoff_width (global, default 0.1) is the roll-off width in
      ! fractional temperature deficit; "Coronal cutoff width: <w>".
      x = (T_fit_floor/max(T,1.0d0) - 1.0d0)/coronal_cutoff_width
      if (x .gt. 26.0d0) then
         coronal_excitation_cutoff = 0.0d0     ! exp(-676) underflows anyway
      else
         coronal_excitation_cutoff = exp(-x*x)
      endif
   endif
   end function coronal_excitation_cutoff

   !---------------------------------------------------!

   !--- Ground-term fine-structure line transfer -----------------------!
   ! HISTORY. AIOLOS (chemistry.cpp:1006) multiplies its gray escape
   ! optical depth by an arbitrary 1e8, which drives beta -> 0 and so
   ! switches metal-line cooling off entirely. EXHALE replaced that with
   ! beta = 1 (fully optically thin), which is the right limit for the
   ! thin wind but overestimates the cooling of the dense base, where the
   ! dominant coolant lines are measured to be optically thick (tau ~ 3,
   ! beta ~ 0.16, for [O I] 63um through the HD 189733 b base). Both were gray in the
   ! XUV continuum opacity, which has nothing to do with line trapping.
   ! What follows computes the LINE-CENTER optical depth of the ground-term
   ! fine-structure lines the cooling assembly solves explicitly, derives
   ! beta from it, and (with "Base IR field" on) the thermal infrared field
   ! those same lines absorb from the lower atmosphere.
   !
   ! SCOPE. Trapping, and the incident base infrared field that goes with
   ! it, are applied to the eight ground-term fine-structure lines of C I,
   ! C II, N II and O I -- exactly the lines for which this module carries
   ! an explicit statistical-equilibrium solution, so emission, opacity and
   ! absorption of the incident field use one set of atomic data. All other
   ! metal-line cooling keeps beta = 1 and no incident field; the validity
   ! note below is why beta = 1 is the right effective treatment for the
   ! permitted resonance lines, and what follows is why the incident field
   ! is left out of the other channels:
   !   - the coronal REMAINDERS added on top of each fine-structure ground
   !     term (the exp(-E/T) sums in cool_CI_ne_func and friends) and the
   !     coronal curves of the other ions are fits to CHIANTI sums over many
   !     transitions with no line list in the code, so the matching
   !     absorption integral int kappa_nu(T) B_nu(T_rad) dnu cannot be
   !     formed from them. Their lowest terms sit at E/k >~ 1.6e4 K, where
   !     B_nu(T0 ~ 1.2e3 K) is down by e^-13, and the heating they would add
   !     was estimated at ~0.4% of the total radiative losses in the
   !     molecular layer of the hot Uranus gate -- small, but an
   !     approximation, not zero.
   !   - Fe II is a precomputed 2-D statistical-equilibrium table
   !     (cool_FeII_ne), which cannot take a radiation field as an argument.
   !   - H3+ (h3p_cooling) is the Miller et al. (2013) fit to the TOTAL
   !     optically thin emission of an LTE molecule, so it has no line list
   !     either -- but its 3-4 um bands ARE within reach of a 1.2e3 K
   !     blackbody, and with metals off it carries 100% of the cooling of
   !     the molecular base, so leaving it emitting into vacuum is not a
   !     small omission. Under the same switch its absorption is taken from
   !     the emission fit itself, at the radiating temperature, which is
   !     what Kirchhoff's law makes of a total emission integral; the
   !     approximation and its range are written at h3p_net_cooling_rate.
   !
   ! VALIDITY for the PERMITTED RESONANCE lines (Mg I 2853, Mg II h&k,
   ! Ca II H&K, Na I D, the Fe II UV multiplets). Those lines DO reach
   ! large line-center depths -- tau0 ~ 1e3-1e5 through the WASP-121 b
   ! wind -- but a large tau0 alone does not suppress their cooling.
   ! Trapping lengthens the random walk; it does not destroy the photon.
   ! In two-level equilibrium with escape probability beta,
   !   Lambda = h nu q_lu ne * [ beta A_ul / (beta A_ul + ne q_ul) ],
   ! and the bracket -- not beta -- is the correction to the optically
   ! thin coronal fits used here. It stays at 1 unless the gas is dense
   ! enough to de-excite the ion before the trapped photon works its way
   ! out, i.e. unless
   !   ne >~ n_crit,eff = beta A_ul / q_ul,   q_ul = coll_rate_prefactor Ups/(gu sqrt(T)).
   ! For these lines A_ul ~ 1e8 s^-1, so even at beta ~ 1e-4 the escape
   ! rate beta A_ul ~ 1e4 s^-1 leaves n_crit,eff ~ 1e10-1e15 cm^-3, one to
   ! six decades above the ne these winds reach (max 4.3e9 cm^-3 in the
   ! WASP-121 b run). beta = 1 is therefore the correct EFFECTIVE
   ! treatment for them here, and the bracket was measured to remove only
   ! 0.02% of the total radiative losses of WASP-121 b, 0.002% of
   ! HD 209458 b and 0.0015% of HD 189733 b, confined to the innermost
   ! ~0.007 R_p. The forbidden fine-structure lines above are the opposite
   ! case (A_ul ~ 1e-5-1e-3 s^-1, n_crit,eff ~ 1e0-1e5 cm^-3), which is
   ! why they, and only they, need beta. The approximation would have to
   ! be revisited for a base an order of magnitude denser in ne, or for a
   ! line whose ne q_ul is comparable to beta A_ul.
   ! Measurement and derivation: md/resonance_line_trapping.md.

   ! Line-center absorption coefficient [cm^-1] of a two-level line whose
   ! lower and upper levels follow the Boltzmann ratio,
   !   kappa_0 = lambda^3/(8 pi^3/2) (g_u/g_l) A_ul n_low (1 - e^-Ek/T)/v_th
   ! for a Doppler core of width v_th = sqrt(2 k T/m) (no turbulent
   ! broadening, so this is an upper bound on kappa_0). The last factor is
   ! the stimulated-emission correction 1 - (g_l n_u)/(g_u n_l). The
   ! Boltzmann assumption holds wherever the result matters: the dense
   ! base has Cdex >> A_ul, so the ground-term levels are in LTE; in the
   ! thin wind kappa_0 is negligible either way. Ek [K] is the transition
   ! energy in temperature units, which also fixes lambda = hc/(k Ek).
   ! atomic_weight_u is the emitter's standard atomic weight in u; the
   ! atomic mass unit itself is the global amu of parameters.f90 (one
   ! definition), so no local copy of it is kept here.
   elemental double precision function line_center_opacity_lte           &
                                        (T,n_low,Ek,A_ul,gu_gl,atomic_weight_u)
   real*8, intent(in) :: T, n_low, Ek, A_ul, gu_gl, atomic_weight_u
   real*8, parameter :: hc_over_k = 1.43877736d0     ! [cm K]
   real*8 :: lam, vth, Ts
   ! T is floored as in the table interpolators, so a transient non-physical
   ! trial temperature cannot turn the Doppler width into a NaN.
   Ts  = max(T, 1.0d0)
   lam = hc_over_k/Ek
   vth = sqrt(2.0d0*kb_erg*Ts/(atomic_weight_u*amu))
   line_center_opacity_lte = lam**3/(8.0d0*pi*sqrt(pi))*gu_gl*A_ul     &
                             *n_low*(1.0d0 - exp(-Ek/Ts))/vth
   end function line_center_opacity_lte

   ! Line-center opacity [cm^-1] of one ground-term fine-structure line.
   ! The lower-level population is the ground-term Boltzmann fraction, the
   ! same partition sum the cooling coefficients emit with; A_ul, E_ul and
   ! g_u/g_l are the CHIANTI values declared at module scope. n_ion is the
   ! TOTAL density of the emitting ion.
   elemental double precision function fine_structure_line_opacity       &
                                         (iline,T,n_ion) result(kap)
   integer, intent(in) :: iline
   real*8,  intent(in) :: T, n_ion
   real*8 :: Ts, Z
   Ts = max(T, 1.0d0)
   select case (iline)
      case (ifs_CI609)     ! 3P1 - 3P0
         Z = 1.0d0 + 3.0d0*exp(-Ek_CI2/Ts) + 5.0d0*exp(-Ek_CI3/Ts)
         kap = line_center_opacity_lte(Ts, n_ion/Z, Ek_CI2,             &
                                       A_CI609, 3.0d0, amu_C)
      case (ifs_CI370)     ! 3P2 - 3P1
         Z = 1.0d0 + 3.0d0*exp(-Ek_CI2/Ts) + 5.0d0*exp(-Ek_CI3/Ts)
         kap = line_center_opacity_lte(Ts,                              &
                  3.0d0*exp(-Ek_CI2/Ts)*n_ion/Z, Ek_CI3 - Ek_CI2,       &
                  A_CI370, 5.0d0/3.0d0, amu_C)
      case (ifs_CII158)    ! 2P3/2 - 2P1/2
         Z = 2.0d0 + 4.0d0*exp(-Ek_CII2/Ts)
         kap = line_center_opacity_lte(Ts, 2.0d0*n_ion/Z, Ek_CII2,      &
                                       A_CII158, 2.0d0, amu_C)
      case (ifs_NII205)    ! 3P1 - 3P0
         Z = 1.0d0 + 3.0d0*exp(-Ek_NII2/Ts) + 5.0d0*exp(-Ek_NII3/Ts)
         kap = line_center_opacity_lte(Ts, n_ion/Z, Ek_NII2,            &
                                       A_NII205, 3.0d0, amu_N)
      case (ifs_NII122)    ! 3P2 - 3P1
         Z = 1.0d0 + 3.0d0*exp(-Ek_NII2/Ts) + 5.0d0*exp(-Ek_NII3/Ts)
         kap = line_center_opacity_lte(Ts,                              &
                  3.0d0*exp(-Ek_NII2/Ts)*n_ion/Z, Ek_NII3 - Ek_NII2,    &
                  A_NII122, 5.0d0/3.0d0, amu_N)
      case (ifs_OI63)      ! 3P1 - 3P2 (inverted term)
         Z = 5.0d0 + 3.0d0*exp(-Ek_OI2/Ts) + exp(-Ek_OI3/Ts)
         kap = line_center_opacity_lte(Ts, 5.0d0*n_ion/Z, Ek_OI2,       &
                                       A_OI63, 0.6d0, amu_O)
      case (ifs_OI145)     ! 3P0 - 3P1
         Z = 5.0d0 + 3.0d0*exp(-Ek_OI2/Ts) + exp(-Ek_OI3/Ts)
         kap = line_center_opacity_lte(Ts,                              &
                  3.0d0*exp(-Ek_OI2/Ts)*n_ion/Z, Ek_OI3 - Ek_OI2,       &
                  A_OI145, 1.0d0/3.0d0, amu_O)
      case (ifs_OI44)      ! 3P0 - 3P2
         Z = 5.0d0 + 3.0d0*exp(-Ek_OI2/Ts) + exp(-Ek_OI3/Ts)
         kap = line_center_opacity_lte(Ts, 5.0d0*n_ion/Z, Ek_OI3,       &
                                       A_OI44, 0.2d0, amu_O)
      case default
         kap = 0.0d0
   end select
   end function fine_structure_line_opacity

   ! Probability that a photon emitted in a static, thermally broadened
   ! plane-parallel layer escapes through the NEARER of its two boundaries.
   ! This is de Jong, Boland & Dalgarno (1980), A&A 91, 68, eq. (B-7),
   !   beta = (1 - e^-a tau)/(2 a tau)            tau < 7,
   !   beta = 1/(4 tau sqrt(ln(tau/sqrt(pi))))    tau >= 7,
   ! with a = 2.34: the single-flight, complete-redistribution escape chance
   ! of a Doppler line, "accurate to within 10% for small and intermediate
   ! tau and exact at very large tau". beta(0) = 1/2, half of the photons
   ! leaving through the near face, which is why the closure adds the two
   ! faces, beta = beta_1(tau_up) + beta_1(tau_down). The two branches meet
   ! where 2 sqrt(ln(tau/sqrt(pi))) = a, i.e. at tau_c = sqrt(pi) exp(a^2/4)
   ! = 6.9676 (the paper rounds it to 7), so the switch is continuous in
   ! value to the e^-a tau_c term the thick branch drops, 8e-8 relative.
   !
   ! ARGUMENT. (B-7) is derived from beta = (1/2) int dx phi(x) E_2[tau
   ! phi(x)] with a NORMALIZED profile, int phi dx = 1, so phi(0) =
   ! 1/sqrt(pi) and its tau is the frequency-integrated depth, sqrt(pi)
   ! times the LINE-CENTRE depth. Written in one convention it is the same
   ! function as Hollenbach & McKee (1979), ApJS 41, 555, eq. (5.10), which
   ! is quoted in line-centre depth: 2 beta(sqrt(pi) tau_0) and their
   ! eps(tau_0) agree to five digits above tau_0 = 1e3.
   !
   ! fine_structure_line_transfer below therefore multiplies the
   ! line-centre depth of each cell by sqrt(pi) before any argument is
   ! formed (since 2026-09-07; until then the LINE-CENTRE column was passed,
   ! and the beta returned was too large by a factor rising from 1 in the
   ! thin limit through 2.11 at the branch point to 1.77 asymptotically:
   ! 1.12 at the largest line-centre depth any shipped case reaches, 0.14
   ! in [O I] 63um at the base of the hot Uranus gate, below 0.1% of the
   ! total radiative losses there, a factor of two in any base thick in
   ! these lines). Measurement: md/resonance_line_trapping.md section 10.
   elemental double precision function line_escape_probability_one_face  &
                                        (tau)
   real*8, intent(in) :: tau
   real*8, parameter :: a = 2.34d0
   real*8 :: x, tau_c
   x = a*tau
   tau_c = sqrt(pi)*exp(0.25d0*a*a)
   if (.not. (x .gt. 1.0d-8)) then
      line_escape_probability_one_face                                  &
                     = 0.5d0*(1.0d0 - 0.5d0*max(x,0.0d0))  ! series limit
   else if (tau .lt. tau_c) then
      line_escape_probability_one_face = 0.5d0*((1.0d0 - exp(-x))/x)
   else
      line_escape_probability_one_face                                  &
                     = 0.5d0*(1.0d0/(2.0d0*tau*sqrt(log(tau/sqrt(pi)))))
   endif
   end function line_escape_probability_one_face

   ! Photon occupation number of a blackbody at T_rad at a line whose
   ! transition energy is Ek = h nu / k [K]:  n = 1/(exp(Ek/T_rad) - 1).
   ! The mean intensity follows as Jbar = (2 h nu^3/c^2) n, so n is the form
   ! the level-balance rates b A (1+n) and b A (g_u/g_l) n need directly.
   elemental double precision function planck_photon_occupation          &
                                        (Ek,T_rad)
   real*8, intent(in) :: Ek, T_rad
   real*8 :: x
   x = Ek/max(T_rad, 1.0d0)
   if (x .gt. 7.0d2) then
      planck_photon_occupation = 0.0d0            ! exp() would overflow
   else
      planck_photon_occupation = 1.0d0/(exp(x) - 1.0d0)
   endif
   end function planck_photon_occupation

   ! Fraction of the sky covered, at radius rr [R_p], by a sphere of radius
   ! r_base [R_p]:  f = 1 - sqrt(1 - (r_base/rr)^2), i.e. twice the usual
   ! dilution factor W. It is 1 at the surface (the lower atmosphere fills
   ! the whole lower hemisphere, W = 1/2) and falls off as
   ! (1/2)(r_base/rr)^2 far away, which is what removes the base infrared
   ! field from the outer wind.
   elemental double precision function base_sky_fraction(rr,r_base)
   real*8, intent(in) :: rr, r_base
   real*8 :: q
   q = min(r_base/max(rr, 1.0d-30), 1.0d0)
   base_sky_fraction = 1.0d0 - sqrt(max(1.0d0 - q*q, 0.0d0))
   end function base_sky_fraction

   ! Line transfer closure of the ground-term fine-structure lines on the
   ! grid: the escape probability beta_fs of each line and the photon
   ! occupation number nbar_fs of the radiation incident on it.
   !
   ! tau_up(j) is the line-center column from the CENTER of cell j to the top
   ! of the domain (half of the emitting cell plus every cell above it) and
   ! tau_dn(j) the column to the bottom of the domain. Being columns rather
   ! than single cell widths they are grid-independent and converge under
   ! refinement, unlike the cell-width gray depth they replace.
   !
   ! WITHOUT the base infrared field (base_ir_field = .false., the default)
   ! the lower atmosphere is treated as a cold, perfectly absorbing floor:
   ! everything emitted downward is lost from the modeled gas and nothing
   ! comes back, so beta = 2 beta_1(tau_up) (the downward hemisphere is
   ! assumed as transparent as the upward one) and nbar = 0.
   !
   ! WITH the field on, the gas below the base is what it physically is: an
   ! optically thick H2 atmosphere at T0. Measured on the converged hot
   ! Uranus molecular solution, one pressure scale height below the base
   ! already carries tau = 1.08 in [O I] 63um and 0.32 in [O I] 145um, and
   ! [C I] 609/370um reach tau = 1 within 2-3 scale heights, so the lower
   ! hemisphere is a blackbody at T0 in every line that matters. Then
   !   beta   = beta_1(tau_up) + beta_1(tau_dn)     (still a loss downward:
   !            the reservoir absorbs it), and
   !   nbar   = beta_1(tau_dn) f_sky(r) / (exp(Ek/T0) - 1),
   ! i.e. a blackbody at T0 covering the fraction f_sky of the sky, hence a
   ! dilution factor f_sky/2, damped on its way to the cell by the same
   ! line opacity (the attenuation is 2 beta_1(tau_dn), which is 1 at
   ! tau_dn = 0, so at the base nbar = (1/2)/(exp(Ek/T0) - 1), the half-sky
   ! value). Using the escape probability for the penetration of the
   ! incident field is the usual reciprocity approximation; it is exact
   ! when the two columns are equal and the reservoir is black, and both
   ! hold in the molecular layer this is written for.
   !
   ! VALIDITY. T0 is the fixed base temperature, so the field is imposed,
   ! not solved for: the layer can be heated by the base but cannot heat it
   ! back. The closure is one-dimensional and static (no velocity shift
   ! between the emitter and the reservoir), which the molecular layer
   ! satisfies (v of order 1-1e2 cm/s there, four decades below the ~1 km/s
   ! Doppler width these lines are broadened by),
   ! and it applies only to the eight fine-structure lines with an explicit
   ! statistical-equilibrium solution -- see the scope note below the H3+
   ! and remainder terms in eval_cool.
   subroutine fine_structure_line_transfer(T,nm,beta_fs,nbar_fs)
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(in)  :: nm
   real*8, dimension(1-Ng:N+Ng,n_fsline), intent(out) :: beta_fs, nbar_fs
   integer :: j, k
   real*8  :: dl, col(n_fsline), b_dn
   real*8  :: dtau(1-Ng:N+Ng,n_fsline), tau_dn(1-Ng:N+Ng,n_fsline)
   ! Frequency-integrated depth of each cell, the argument (B-7) is written
   ! in: sqrt(pi) times the line-centre depth of a Doppler profile (header).
   ! Floored at zero: the ionization solve can leave a trace species with a
   ! small negative density, and an optical depth cannot be negative.
   do j = 1-Ng,N+Ng
      dl = dr_j(j)*R0*sqrt(pi)
      do k = 1,n_fsline
         dtau(j,k) = max(fine_structure_line_opacity(k, T(j),           &
                                     nm(j,fsline_ion(k)))*dl, 0.0d0)
      enddo
   enddo
   ! Column downward, from the bottom of the domain up
   col = 0.0d0
   do j = 1-Ng,N+Ng
      do k = 1,n_fsline
         tau_dn(j,k) = col(k) + 0.5d0*dtau(j,k)
         col(k) = col(k) + dtau(j,k)
      enddo
   enddo
   ! Column upward, from the top of the domain down; assemble beta and nbar
   col = 0.0d0
   nbar_fs = 0.0d0
   do j = N+Ng,1-Ng,-1
      do k = 1,n_fsline
         if (base_ir_field) then
            b_dn = line_escape_probability_one_face(tau_dn(j,k))
            beta_fs(j,k) = line_escape_probability_one_face             &
                              (col(k) + 0.5d0*dtau(j,k)) + b_dn
            nbar_fs(j,k) = b_dn*base_sky_fraction(r(j), 1.0d0)          &
                           *planck_photon_occupation(fsline_Ek(k), T0)
         else
            beta_fs(j,k) = 2.0d0*line_escape_probability_one_face       &
                              (col(k) + 0.5d0*dtau(j,k))
         endif
         col(k) = col(k) + dtau(j,k)
      enddo
   enddo
   end subroutine fine_structure_line_transfer

   !---------------------------------------------------!

   !--- Metal recombination rates (Badnell 2006 RR + adf09 DR) ---!
   ! Total (RR+DR) recombination, replacing the Aldrovandi & Pequignot
   ! (1973) power-law fits (which had no dielectronic term). The fit
   ! constants bad_* are declared at module scope above.

   ! Look up the Badnell fit for a recombined (daughter) ion.
   pure subroutine lookup_rec(daughter,f,known)
   character(len=*), intent(in)  :: daughter
   type(rec_fit),    intent(out) :: f
   logical,          intent(out) :: known
   known = .true.
   select case (daughter)
      case ('CI');   f = bad_CI
      case ('CII');  f = bad_CII
      case ('NI');   f = bad_NI
      case ('NII');  f = bad_NII
      case ('OI');   f = bad_OI
      case ('OII');  f = bad_OII
      case ('MgI');  f = bad_MgI
      case ('MgII'); f = bad_MgII
      case ('SiI');  f = bad_SiI
      case ('SiII'); f = bad_SiII
      case ('CaII'); f = bad_CaII
      case ('NaI');  f = bad_NaI
      case ('KI');   f = bad_KI
      case ('SI');   f = bad_SI
      ! NOTE: 'CaI' (K-like Ca II + e) is intentionally absent; it is
      ! handled by the Shull & Van Steenberg 1982 power law in
      ! alpha_rr_metal, with no dielectronic term.
      case default;  known = .false.
   end select
   end subroutine lookup_rec

   ! Radiative recombination only [cm^3/s] (Badnell 2006 form)
   pure double precision function alpha_rr_metal(daughter,T)
   character(len=*), intent(in) :: daughter
   real*8,           intent(in) :: T
   type(rec_fit) :: f
   logical :: known
   real*8  :: sT0,sT1,b

   ! Ca II + e -> Ca I (K-like Ca+) is not in the Badnell tabulation
   ! (its most electron-rich recombining ion is Ar-like, eighteen
   ! electrons; there is no K-like row), so the Shull &
   ! Van Steenberg (1982, ApJS 48, 95) radiative-recombination power law
   ! is taken, in the machine form of Verner's rrfit.f (version 4, 1999),
   !
   !     alpha_RR = rrec(1,Z,N) * (T/1e4 K)^(-rrec(2,Z,N)) ,
   !
   ! whose second index N is the electron count of the RECOMBINED ion:
   ! rnew(:,1,1) is the H+ + e -> H I fit and rnew(:,2,2) the He+ + e ->
   ! He I one.  Ca I has twenty electrons, so the row is rrec(:,20,20) =
   ! (1.120e-13, 0.9000), not rrec(:,20,19), which forms Ca II.  The
   ! ordering is monotone in the charge of the recombining ion, as a
   ! radiative recombination coefficient is: 1.12e-13 forming Ca I,
   ! 6.78e-13 forming Ca II, 3.96e-12 forming Ca III at 1e4 K.
   ! No dielectronic fit exists for the K-like sequence in the Badnell
   ! adf09 set, so alpha_dr_metal('CaI') returns 0.
   if (trim(daughter) == 'CaI') then
      if (T .le. 0.0d0) then
         alpha_rr_metal = 0.0d0
      else
         alpha_rr_metal = 1.120d-13*(T/1.0d4)**(-0.900d0)
      endif
      return
   endif

   call lookup_rec(daughter,f,known)
   if (.not.known .or. T.le.0.0d0) then
      alpha_rr_metal = 0.0d0
      return
   endif
   sT0 = sqrt(T/f%rr_T0)
   sT1 = sqrt(T/f%rr_T1)
   b   = f%rr_B + f%rr_C*exp(-f%rr_T2/T)
   alpha_rr_metal = f%rr_A                                            &
        / ( sT0*(1.0d0+sT0)**(1.0d0-b)*(1.0d0+sT1)**(1.0d0+b) )
   end function alpha_rr_metal

   ! Dielectronic recombination only [cm^3/s] (Badnell adf09 form)
   pure double precision function alpha_dr_metal(daughter,T)
   character(len=*), intent(in) :: daughter
   real*8,           intent(in) :: T
   type(rec_fit) :: f
   logical :: known
   real*8  :: acc
   integer :: i
   call lookup_rec(daughter,f,known)
   if (.not.known .or. T.le.0.0d0 .or. f%dr_n.le.0) then
      alpha_dr_metal = 0.0d0
      return
   endif
   acc = 0.0d0
   do i = 1,f%dr_n
      acc = acc + f%dr_c(i)*exp(-f%dr_e(i)/T)
   enddo
   alpha_dr_metal = T**(-1.5d0)*acc
   end function alpha_dr_metal

   ! Badnell RR + DR of one tabulated fit [cm^3/s], from the fit itself
   ! (the same arithmetic as alpha_rr_metal + alpha_dr_metal, without the
   ! lookup by name).
   pure double precision function badnell_total(f,T)
   type(rec_fit), intent(in) :: f
   real*8,        intent(in) :: T
   real*8  :: sT0, sT1, b, acc
   integer :: i
   if (T .le. 0.0d0) then
      badnell_total = 0.0d0
      return
   endif
   sT0 = sqrt(T/f%rr_T0)
   sT1 = sqrt(T/f%rr_T1)
   b   = f%rr_B + f%rr_C*exp(-f%rr_T2/T)
   badnell_total = f%rr_A                                             &
        / ( sT0*(1.0d0+sT0)**(1.0d0-b)*(1.0d0+sT1)**(1.0d0+b) )
   acc = 0.0d0
   do i = 1,f%dr_n
      acc = acc + f%dr_c(i)*exp(-f%dr_e(i)/T)
   enddo
   badnell_total = badnell_total + T**(-1.5d0)*acc
   end function badnell_total

   ! Total recombination = RR + DR [cm^3/s]
   pure double precision function alpha_rec_metal(daughter,T)
   character(len=*), intent(in) :: daughter
   real*8,           intent(in) :: T
   alpha_rec_metal = alpha_rr_metal(daughter,T)                      &
                   + alpha_dr_metal(daughter,T)
   end function alpha_rec_metal

   !--------------!

   ! Iron recombination (RR+DR total), Huang+2023 (ApJ 951, 123) Eqs (5)-(6).
   ! Form:  alpha(T) = A*T^-1.5*exp(-T0/T)*(1 + B*exp(-T1/T))      [DR]
   !                 + C*(T/1e4)^-eta                              [RR]
   ! with T in K and alpha in cm^3/s. Coefficients read from the published
   ! paper, ApJ 951, 123, Eqs. (5) and (6). The rest of the metal grid uses
   ! Badnell RR+DR (alpha_rec_metal); iron cannot, and this is not a
   ! deferral. Fe I comes from recombining Fe II (Mn-like, 25 electrons) and
   ! Fe II from Fe III (Cr-like, 24), while the most electron-rich
   ! recombining ion in the Badnell RR and DR tabulations is Ar-like
   ! (18 electrons) -- there is no Badnell fit to switch to for either
   ! stage. Huang's analytic form is therefore the rate for iron, which is
   ! also what the Huang reproduction plan prescribes for the stages where
   ! Badnell coverage stops.

   ! Fe II + e -> Fe I  (rate producing the Fe I daughter).
   pure double precision function alpha_rec_FeI_Huang(T)
   real*8, intent(in) :: T
   alpha_rec_FeI_Huang = 2.833d-8*T**(-1.5d0)*exp(-5.731d4/T)            &
                         *(1.0d0 + 1.383d4*exp(-120.4d0/T))             &
                       + 1.248d-12*(T/1.0d4)**(-0.485d0)
   end function alpha_rec_FeI_Huang

   ! Fe III + e -> Fe II  (rate producing the Fe II daughter).
   pure double precision function alpha_rec_FeII_Huang(T)
   real*8, intent(in) :: T
   alpha_rec_FeII_Huang = 1.094d-5*T**(-1.5d0)*exp(-1.490d4/T)           &
                          *(1.0d0 + 36.74d0*exp(-1.153d5/T))            &
                        + 1.728d-12*(T/1.0d4)**(-0.618d0)
   end function alpha_rec_FeII_Huang

   ! RECOMBINATION COEFFICIENT OF METAL ION i [cm^3 s^-1], the recombining
   ! ion in the canonical species_table order (ion i + e -> the stage below
   ! it): Badnell radiative + dielectronic (badnell_total of the fit of the
   ! recombined ion, the arithmetic of alpha_rec_metal), the Shull & Van
   ! Steenberg (1982) power law for
   ! Ca II + e (alpha_rr_metal('CaI')), Huang et al. (2023) for iron. THE
   ! CASE LIST LIVES HERE ONLY: the ionization balance (through
   ! rec_coeff_by_ion_range) and the recombination cooling of
   ! radiative_cooling_of_cell both call it. Ions that do not recombine in
   ! the ladder (the neutrals) return 0.
   pure double precision function metal_recombination_coefficient(i,T)
   integer, intent(in) :: i
   real*8,  intent(in) :: T
   select case (i)
      case (im_CII);   metal_recombination_coefficient = badnell_total(bad_CI,  T)
      case (im_CIII);  metal_recombination_coefficient = badnell_total(bad_CII, T)
      case (im_NII);   metal_recombination_coefficient = badnell_total(bad_NI,  T)
      case (im_NIII);  metal_recombination_coefficient = badnell_total(bad_NII, T)
      case (im_OII);   metal_recombination_coefficient = badnell_total(bad_OI,  T)
      case (im_OIII);  metal_recombination_coefficient = badnell_total(bad_OII, T)
      case (im_MgII);  metal_recombination_coefficient = badnell_total(bad_MgI, T)
      case (im_MgIII); metal_recombination_coefficient = badnell_total(bad_MgII,T)
      case (im_SiII);  metal_recombination_coefficient = badnell_total(bad_SiI, T)
      case (im_SiIII); metal_recombination_coefficient = badnell_total(bad_SiII,T)
      case (im_CaII);  metal_recombination_coefficient = alpha_rr_metal('CaI', T)
      case (im_CaIII); metal_recombination_coefficient = badnell_total(bad_CaII,T)
      case (im_NaII);  metal_recombination_coefficient = badnell_total(bad_NaI, T)
      case (im_KII);   metal_recombination_coefficient = badnell_total(bad_KI,  T)
      case (im_S_II);  metal_recombination_coefficient = badnell_total(bad_SI,  T)
      case (im_FeII);  metal_recombination_coefficient = alpha_rec_FeI_Huang(T)
      case (im_FeIII); metal_recombination_coefficient = alpha_rec_FeII_Huang(T)
      case default;    metal_recombination_coefficient = 0.0d0
   end select
   end function metal_recombination_coefficient

   !---------------------------------------------------!

   !--- Metal collisional ionization rates (Voronov 1997) ---!
   ! Q(T) = A * (1 + P*sqrt(U)) / (X+U) * U^K * exp(-U)
   ! U    = dE / (kB * T)
   ! One row (dE, P, A, X, K) per ion from Voronov (1997, ADNDT 65, 1,
   ! Table I, READ from the published article; every row agrees with the
   ! author-distributed transcription, D. A. Verner's cfit.f, version 2,
   ! 24 March 1997). Validity 1 eV - 20 keV for every row here (Table I,
   ! Tmin and Tmax), and below 1 eV as written at voronov_ci.
   ! dE is a COEFFICIENT OF THE FIT, as for H and He (see ci_HI_new): it
   ! enters through U^K, exp(-U) and X + U, and A, X, K, P were fitted
   ! with it, so it is used as published, not replaced by the measured
   ! ionization potential (Mg II 15.2 eV in the fit against 15.035 eV
   ! measured; the substitution moved the Mg II rate by x1.47 at 5e3 K and
   ! x1.21 at 1e4 K, MEASURED). Each ionization removes the measured
   ! potential of the ion (species_table mion_ethr) from the electron gas
   ! (radiative_cooling_of_cell, channel 2), the same energy the ion's
   ! photoelectrons are charged against.

   ! THE KINETIC ENERGY REMOVED BY THE RECOMBINATIONS OF METAL ION i
   ! [erg cm^3 s^-1 per (n_e n_i)], k T alpha (3/2 + dln alpha/dln T) on
   ! metal_recombination_coefficient, with the logarithmic slope of each
   ! published form written analytically (so one evaluation of the
   ! exponentials serves the coefficient and its slope):
   !   Badnell RR, alpha = A/[s0 (1+s0)^(1-b) (1+s1)^(1+b)], s = (T/T_x)^1/2,
   !     b = B + C exp(-T2/T): dln/dln T = -1/2 - (1-b) s0/(2(1+s0))
   !     - (1+b) s1/(2(1+s1)) - (db/dln T) ln[(1+s0)/(1+s1)],
   !     db/dln T = C (T2/T) exp(-T2/T);
   !   DR, T^-3/2 sum c_i exp(-E_i/T): d/dln T = -3/2 alpha_DR
   !     + T^-3/2 sum c_i (E_i/T) exp(-E_i/T);
   !   power law C (T/1e4)^-eta: -eta.
   ! The physics probe metal_ionization_energy_ledger checks it against a
   ! centred difference of metal_recombination_coefficient.
   pure double precision function metal_recombination_energy_rate(i,T)    &
                                                  result(beta)
   integer, intent(in) :: i
   real*8,  intent(in) :: T
   real*8 :: a, dadlnT
   select case (i)
      case (im_CII);   call badnell_moments(bad_CI,  T, a, dadlnT)
      case (im_CIII);  call badnell_moments(bad_CII, T, a, dadlnT)
      case (im_NII);   call badnell_moments(bad_NI,  T, a, dadlnT)
      case (im_NIII);  call badnell_moments(bad_NII, T, a, dadlnT)
      case (im_OII);   call badnell_moments(bad_OI,  T, a, dadlnT)
      case (im_OIII);  call badnell_moments(bad_OII, T, a, dadlnT)
      case (im_MgII);  call badnell_moments(bad_MgI, T, a, dadlnT)
      case (im_MgIII); call badnell_moments(bad_MgII,T, a, dadlnT)
      case (im_SiII);  call badnell_moments(bad_SiI, T, a, dadlnT)
      case (im_SiIII); call badnell_moments(bad_SiII,T, a, dadlnT)
      case (im_CaII)
         a = alpha_rr_metal('CaI', T);  dadlnT = -0.900d0*a
      case (im_CaIII); call badnell_moments(bad_CaII,T, a, dadlnT)
      case (im_NaII);  call badnell_moments(bad_NaI, T, a, dadlnT)
      case (im_KII);   call badnell_moments(bad_KI,  T, a, dadlnT)
      case (im_S_II);  call badnell_moments(bad_SI,  T, a, dadlnT)
      case (im_FeII)
         call huang_iron_moments(2.833d-8, 5.731d4, 1.383d4, 120.4d0,      &
                                 1.248d-12, 0.485d0, T, a, dadlnT)
      case (im_FeIII)
         call huang_iron_moments(1.094d-5, 1.490d4, 36.74d0, 1.153d5,      &
                                 1.728d-12, 0.618d0, T, a, dadlnT)
      case default
         a = 0.0d0;  dadlnT = 0.0d0
   end select
   beta = kb_erg*T*(1.5d0*a + dadlnT)
   end function metal_recombination_energy_rate

   ! alpha and dalpha/dln T of one Badnell RR + DR fit.
   pure subroutine badnell_moments(f, T, alpha, dadlnT)
   type(rec_fit), intent(in)  :: f
   real*8,        intent(in)  :: T
   real*8,        intent(out) :: alpha, dadlnT
   real*8  :: sT0, sT1, e2, b, l0, l1, a_rr, slope_rr, ex, acc, acc1
   integer :: k
   if (T .le. 0.0d0) then
      alpha = 0.0d0;  dadlnT = 0.0d0
      return
   endif
   sT0  = sqrt(T/f%rr_T0)
   sT1  = sqrt(T/f%rr_T1)
   e2   = exp(-f%rr_T2/T)
   b    = f%rr_B + f%rr_C*e2
   l0   = log(1.0d0 + sT0)
   l1   = log(1.0d0 + sT1)
   a_rr = f%rr_A/(sT0*exp((1.0d0 - b)*l0 + (1.0d0 + b)*l1))
   slope_rr = -0.5d0 - (1.0d0 - b)*0.5d0*sT0/(1.0d0 + sT0)              &
              - (1.0d0 + b)*0.5d0*sT1/(1.0d0 + sT1)                       &
              - f%rr_C*(f%rr_T2/T)*e2*(l1 - l0)
   acc  = 0.0d0
   acc1 = 0.0d0
   do k = 1,f%dr_n
      ex   = f%dr_c(k)*exp(-f%dr_e(k)/T)
      acc  = acc + ex
      acc1 = acc1 + ex*f%dr_e(k)/T
   enddo
   alpha  = a_rr + T**(-1.5d0)*acc
   dadlnT = a_rr*slope_rr + T**(-1.5d0)*(acc1 - 1.5d0*acc)
   end subroutine badnell_moments

   ! alpha and dalpha/dln T of the Huang et al. (2023) iron form
   ! A T^-3/2 exp(-T0/T) (1 + B exp(-T1/T)) + C (T/1e4)^-eta.
   pure subroutine huang_iron_moments(A, T0h, B, T1h, C, eta, T, alpha,   &
                                      dadlnT)
   real*8, intent(in)  :: A, T0h, B, T1h, C, eta, T
   real*8, intent(out) :: alpha, dadlnT
   real*8 :: e0, e1, dr, rr
   e0 = exp(-T0h/T)
   e1 = exp(-T1h/T)
   dr = A*T**(-1.5d0)*e0*(1.0d0 + B*e1)
   rr = C*(T/1.0d4)**(-eta)
   alpha  = dr + rr
   dadlnT = dr*(-1.5d0 + T0h/T) + A*T**(-1.5d0)*e0*B*e1*T1h/T - eta*rr
   end subroutine huang_iron_moments

   ! COLLISIONAL IONIZATION COEFFICIENT OF METAL ION i [cm^3 s^-1] (ion i
   ! -> the stage above it), voronov_ci with Voronov's row (dE, P, A, X, K)
   ! for the ion (Table I of Voronov 1997, rows as named in the physics
   ! probe voronov_table_rows, which checks every row). P = 0 for C I, N I,
   ! N II, O I, Mg I, Mg II, Ca I, Ca II and Fe I; P = 1 for C II, O II,
   ! Si I, Si II, Na I, K I, S I and Fe II. THE CASE
   ! LIST LIVES HERE ONLY: the ionization balance (ion_coeff_by_ion_range)
   ! and the ionization-energy cooling of radiative_cooling_of_cell call it.
   ! The top stage of each ladder returns 0.
   pure double precision function metal_collisional_ionization_coefficient &
                                    (i,T) result(c)
   integer, intent(in) :: i
   real*8,  intent(in) :: T
   select case (i)
      case (im_CI)
         c = voronov_ci(T, 11.3d0, 0.0d0, 6.85d-8, 0.193d0, 0.25d0)
      case (im_CII)
         c = voronov_ci(T, 24.4d0, 1.0d0, 1.86d-8, 0.286d0, 0.24d0)
      case (im_OI)
         c = voronov_ci(T, 13.6d0, 0.0d0, 3.59d-8, 0.073d0, 0.34d0)
      case (im_OII)
         c = voronov_ci(T, 35.1d0, 1.0d0, 1.39d-8, 0.212d0, 0.22d0)
      case (im_NI)
         c = voronov_ci(T, 14.5d0, 0.0d0, 4.82d-8, 0.0652d0, 0.42d0)
      case (im_NII)
         c = voronov_ci(T, 29.6d0, 0.0d0, 2.98d-8, 0.31d0, 0.3d0)
      case (im_MgI)
         c = voronov_ci(T, 7.6d0, 0.0d0, 6.21d-7, 0.592d0, 0.39d0)
      case (im_MgII)
         c = voronov_ci(T, 15.2d0, 0.0d0, 1.92d-8, 0.0027d0, 0.85d0)
      case (im_SiI)
         c = voronov_ci(T, 8.2d0, 1.0d0, 1.88d-7, 0.376d0, 0.25d0)
      case (im_SiII)
         c = voronov_ci(T, 16.4d0, 1.0d0, 6.43d-8, 0.632d0, 0.2d0)
      case (im_CaI)
         c = voronov_ci(T, 6.1d0, 0.0d0, 4.40d-7, 0.848d0, 0.33d0)
      case (im_CaII)
         c = voronov_ci(T, 11.9d0, 0.0d0, 5.22d-8, 0.151d0, 0.34d0)
      case (im_NaI)
         c = voronov_ci(T, 5.1d0, 1.0d0, 1.01d-7, 0.275d0, 0.23d0)
      case (im_KI)
         c = voronov_ci(T, 4.3d0, 1.0d0, 2.02d-7, 0.272d0, 0.31d0)
      case (im_S_I)
         c = voronov_ci(T, 10.4d0, 1.0d0, 5.49d-8, 0.1d0, 0.25d0)
      case (im_FeI)
         c = voronov_ci(T, 7.9d0, 0.0d0, 2.52d-7, 0.701d0, 0.25d0)
      case (im_FeII)
         c = voronov_ci(T, 16.2d0, 1.0d0, 2.21d-8, 0.033d0, 0.45d0)
      case default
         c = 0.0d0
   end select
   end function metal_collisional_ionization_coefficient

   ! ---------------------------------------------------------------- !

   ! The coefficients below are the ONE definition of each rate: the grid
   ! routines above evaluate these same elemental functions over their cell
   ! range, and the cooling of a cell (radiative_cooling_of_cell, at the end
   ! of this module) calls them directly.

   ! Recombination coefficient of HII: case B, or the published set the
   ! "Atomic rate set:" key selects.
   elemental double precision function alpha_rec_HII_B(T)
   real*8, intent(in) :: T
   real*8 :: xl

   if (atomic_rate_set_k22) then
      ! Koskinen et al. 2022, Table 1 R1 (Storey & Hummer 1995 power law)
      alpha_rec_HII_B = rk_R1_Hp_rec(T)
   else if (legacy_hhe_rates) then
      ! Hui & Gnedin 1997, MNRAS 292, 27, case-B fit. The 157807 K is
      ! that paper's own H I threshold temperature, a constant of the fit
      ! (lambda = 2 T_HI/T), and is kept as published rather than replaced
      ! by Ry_over_kB.
      xl = 2.0d0*157807.0d0/T
      alpha_rec_HII_B = 2.753d-14*xl**1.5d0/(1.0d0+(xl/2.740d0)**0.407d0)**2.242d0
   else
      alpha_rec_HII_B = alphaB_HII_new(T)
   endif

   end function alpha_rec_HII_B

   !--------------!

   ! Case B recombination coefficient of HeII
   elemental double precision function alpha_rec_HeII_B(T)
   real*8, intent(in) :: T
   real*8 :: xl

   if (atomic_rate_set_k22) then
      ! Koskinen et al. 2022, Table 1 R2 (Storey & Hummer 1995 power law).
      ! With He_rec_coupling on, the ground capture added on top of this
      ! is the Milne alpha_1 of this code -- Table 1 carries no such split
      ! -- so the closest like-for-like run also sets "He_rec_coupling:
      ! False".
      alpha_rec_HeII_B = rk_R2_Hep_rec(T)
   else if (legacy_hhe_rates) then
      ! Hui & Gnedin 1997, MNRAS 292, 27, case-B fit
      xl         = 2.0d0*285335.0d0/T
      alpha_rec_HeII_B = 1.26d-14*xl**0.750d0
   else
      alpha_rec_HeII_B = alphaB_HeII_new(T)
   endif

   end function alpha_rec_HeII_B

   !--------------!

   ! Case B recombination coefficient of HeIII
   elemental double precision function alpha_rec_HeIII_B(T)
   real*8, intent(in) :: T
   real*8 :: xl

   if (legacy_hhe_rates) then
      ! Hui & Gnedin 1997, MNRAS 292, 27, case-B fit
      xl = 2.0d0*631515.0d0/T
      alpha_rec_HeIII_B = 2.0d0*2.753d-14*xl**1.5d0                       &
                        /(1.0d0+(xl/2.740d0)**0.407d0)**2.242d0
   else
      alpha_rec_HeIII_B = alphaB_HeIII_new(T)
   endif

   end function alpha_rec_HeIII_B

   !--------------!

   ! THE GROUND CAPTURES of the three ions [cm^3 s^-1], one definition
   ! each (ground_capture_milne), whichever "Atomic rate set:" gives
   ! alpha_B.
   elemental double precision function alpha_1_HI(T)
   real*8, intent(in) :: T
   alpha_1_HI = ground_capture_milne(T, 1)
   end function alpha_1_HI

   elemental double precision function alpha_1_HeI(T)
   real*8, intent(in) :: T
   alpha_1_HeI = ground_capture_milne(T, 2)
   end function alpha_1_HeI

   elemental double precision function alpha_1_HeII(T)
   real*8, intent(in) :: T
   alpha_1_HeII = ground_capture_milne(T, 3)
   end function alpha_1_HeII

   !--------------!

   ! He I CASE-B RECOMBINATION BY SPIN SYSTEM. A capture into a triplet
   ! level cascades into 2^3S (the triplets reach 1^1S only through the
   ! 1.27e-4 s^-1 magnetic-dipole decay of 2^3S), a capture into an
   ! excited singlet cascades into 1^1S. Hummer & Storey (1998, MNRAS 297,
   ! 1073, Table 5, READ from the published article) tabulate
   ! T^1/2 sum_{n>=2} alpha(T; nlS) for S = 1 and S = 3 separately over
   ! log T = 1.0 (0.2) 4.4; their sum is their alpha_B. The triplet share
   ! of the case-B captures is their ratio, interpolated linearly in log T
   ! and held at the table ends outside it: 0.755 at 10 K, 0.771 at 1e4 K,
   ! 0.781 at 2.5e4 K. Above 2.5e4 K the share is held; there the
   ! dielectronic term of alpha_B (Badnell) grows, whose autoionizing
   ! states are populated in proportion to their statistical weights, the
   ! same 3:1 that the held share is close to. It splits the ONE alpha_B of
   ! the balance, so singlets + triplets + ground = alpha_A (Badnell)
   ! exactly. This replaces the alpha_3 = 2.10e-13 (T/1e4)^-0.778 fit
   ! quoted with Oklopcic & Hirata (2018), which is 0.989 of the Hummer &
   ! Storey triplet coefficient at 1e4 K, and their alpha_1 = 1.54e-13
   ! (T/1e4)^-0.486: that is the capture into 1^1S alone (0.966 of Hummer
   ! & Storey's 1.595e-13 at 1e4 K; Oklopcic & Hirata, eq. 14, feed only
   ! it into the singlet), which left the excited-singlet captures
   ! (6.3e-14 at 1e4 K) out of the network when the coupling was off.
   elemental double precision function case_b_triplet_share_HeI(T)
   real*8, intent(in) :: T
   ! Hummer & Storey (1998) Table 5: T^1/2 sum_{n>=2} alpha [cm^3 s^-1
   ! K^1/2], singlets and triplets, log T = 1.0, 1.2, ..., 4.4
   real*8, parameter :: s1(18) = [ 2.276d-11, 2.167d-11, 2.056d-11,       &
        1.943d-11, 1.830d-11, 1.716d-11, 1.602d-11, 1.488d-11, 1.375d-11,    &
        1.263d-11, 1.152d-11, 1.042d-11, 9.343d-12, 8.292d-12, 7.276d-12,    &
        6.306d-12, 5.393d-12, 4.550d-12 ]
   real*8, parameter :: s3(18) = [ 7.008d-11, 6.680d-11, 6.347d-11,       &
        6.009d-11, 5.669d-11, 5.328d-11, 4.987d-11, 4.648d-11, 4.310d-11,    &
        3.976d-11, 3.646d-11, 3.322d-11, 3.005d-11, 2.699d-11, 2.404d-11,    &
        2.124d-11, 1.862d-11, 1.618d-11 ]
   real*8  :: pos, f, sh_lo, sh_hi
   integer :: k
   pos = (log10(max(T, 1.0d0)) - 1.0d0)/0.2d0 + 1.0d0
   if (.not. (pos > 1.0d0)) then
      k = 1;  f = 0.0d0
   else if (pos .ge. 18.0d0) then
      k = 17; f = 1.0d0
   else
      k = int(pos); f = pos - dble(k)
   endif
   sh_lo = s3(k)/(s1(k) + s3(k))
   sh_hi = s3(k+1)/(s1(k+1) + s3(k+1))
   case_b_triplet_share_HeI = sh_lo + f*(sh_hi - sh_lo)
   end function case_b_triplet_share_HeI

   ! Captures into the triplets (all ending in 2^3S) and into the excited
   ! singlets (all ending in 1^1S), the two parts of alpha_B [cm^3 s^-1].
   elemental double precision function alpha_rec_HeII_23S(T)
   real*8, intent(in) :: T
   alpha_rec_HeII_23S = case_b_triplet_share_HeI(T)*alpha_rec_HeII_B(T)
   end function alpha_rec_HeII_23S

   elemental double precision function alpha_rec_HeII_excited_singlets(T)
   real*8, intent(in) :: T
   alpha_rec_HeII_excited_singlets =                                      &
        (1.0d0 - case_b_triplet_share_HeI(T))*alpha_rec_HeII_B(T)
   end function alpha_rec_HeII_excited_singlets

   !--------------!

   ! THE RECOMBINATION COEFFICIENTS THE IONIZATION BALANCE REMOVES IONS
   ! WITH [cm^3 s^-1]. A capture into the ground level emits a photon
   ! above the ionization edge of the recombined species; where that
   ! photon is re-absorbed by the same species in the cell it is emitted
   ! in, the capture and the re-ionization cancel (case B on the spot), and
   ! where it leaves or is taken by another absorber the capture is a net
   ! recombination (case A). y_gnd (He II -> He I), y_HI (H II -> H I)
   ! and y_HeII (He III -> He II) are the fractions of those photons NOT
   ! re-absorbed by the recombined species in the cell (utils_ion_eq:
   ! ground_capture_escape_weights), functions of the absorber densities
   ! and the cell width, never of T; zero for a switch that is off (case B).
   !
   ! alpha_rec_HeII_into_singlets is what the balance writes as the
   ! singlet coefficient (rcheiiB):
   !   no metastable   alpha_B + y_gnd alpha_1
   !   metastable      y_gnd alpha_1 + (1 - share) alpha_B
   ! (share: case_b_triplet_share_HeI), and alpha_rec_HeII_net adds the
   ! capture into the triplets, share alpha_B, where the metastable is a
   ! level of the network. The two configurations therefore remove He+ at
   ! the same total rate, alpha_B + y_gnd alpha_1.
   elemental double precision function alpha_rec_HeII_into_singlets(T,y_gnd)
   real*8, intent(in) :: T, y_gnd

   if (thereis_HeITR) then
      alpha_rec_HeII_into_singlets = y_gnd*alpha_1_HeI(T)                 &
                                   + alpha_rec_HeII_excited_singlets(T)
   else
      alpha_rec_HeII_into_singlets = alpha_rec_HeII_B(T)                  &
                                   + y_gnd*alpha_1_HeI(T)
   endif

   end function alpha_rec_HeII_into_singlets

   elemental double precision function alpha_rec_HeII_net(T,y_gnd)
   real*8, intent(in) :: T, y_gnd

   alpha_rec_HeII_net = alpha_rec_HeII_into_singlets(T,y_gnd)
   if (thereis_HeITR) alpha_rec_HeII_net = alpha_rec_HeII_net                &
                                         + alpha_rec_HeII_23S(T)

   end function alpha_rec_HeII_net

   ! H II -> H I, the same construction: case B plus the ground captures
   ! whose photons are not re-absorbed by H I in the cell.
   elemental double precision function alpha_rec_HII_net(T,y_HI)
   real*8, intent(in) :: T, y_HI
   alpha_rec_HII_net = alpha_rec_HII_B(T) + y_HI*alpha_1_HI(T)
   end function alpha_rec_HII_net

   ! He III -> He II, the same construction.
   elemental double precision function alpha_rec_HeIII_net(T,y_HeII)
   real*8, intent(in) :: T, y_HeII
   alpha_rec_HeIII_net = alpha_rec_HeIII_B(T) + y_HeII*alpha_1_HeII(T)
   end function alpha_rec_HeIII_net

   !--------------!

   ! Kinetic energy removed from the electron gas by captures of rate
   ! coefficient alpha [erg cm^3 s^-1]: k T [(3/2) alpha + dalpha/dln T]
   ! (the relation and its check are at dlnT_capture). a_lo, a_0 and a_hi
   ! are alpha at T exp(-dlnT_capture), T and T exp(+dlnT_capture).
   elemental double precision function capture_energy_loss_rate           &
                                         (T,a_lo,a_0,a_hi)
   real*8, intent(in) :: T, a_lo, a_0, a_hi
   capture_energy_loss_rate = kb_erg*T*(1.5d0*a_0                          &
                            + (a_hi - a_lo)/(2.0d0*dlnT_capture))
   end function capture_energy_loss_rate

   !--------------!

   ! RECOMBINATION COOLING [erg cm^3 s^-1 per (n_e n_ion)]: the kinetic
   ! energy of the captured electrons, from the relation at dlnT_capture
   ! applied to the coefficient the ionization balance uses, so the
   ! electrons the balance removes are charged the energy they carry, for
   ! every "Atomic rate set:" and helium configuration. The escape weight
   ! (y_HI, y_gnd, y_HeII) is FROZEN: it multiplies the coefficient at T
   ! and at the two shifted temperatures alike, so dalpha/dln T is the
   ! slope of the capture coefficients alone, and a caller that
   ! finite-differences the cooling in T at fixed densities (the
   ! semi-implicit energy update) differentiates the same function. The
   ! kinetic energy of a ground capture whose photon is absorbed by
   ! another species is returned there with the photoelectron
   ! (utils_ion_eq: recombination_radiation_absorbed); one re-absorbed by
   ! the recombined species removes and returns the same energy and is in
   ! neither.
   !
   ! H II. The Hui & Gnedin (1997, MNRAS 292, 27, Appendix A) case-B
   ! cooling fit agrees with this relation applied to its own alpha_B to
   ! 1% but departs from the relation applied to the default balance
   ! coefficient by -3.3% to +3.7% over 2e3-1e5 K, beyond its stated 2%
   ! accuracy (MEASURED). VALIDITY: the accuracy of the balance
   ! coefficient and of its logarithmic slope; at y_HI = 0 and with the
   ! default coefficient the result agrees with Hummer (1994, Table 1)
   ! beta_B to -0.1% ... +2.0% over 1e3-4e4 K, -1.5% at 6.3e4 K and -4.4%
   ! at 1e5 K (MEASURED before the Milne alpha_1; the physics probe
   ! recombination_energy_loss holds it to 2.5%).
   elemental double precision function lambda_rec_HII(T,y_HI)
   real*8, intent(in) :: T, y_HI
   lambda_rec_HII = capture_energy_loss_rate(T,                            &
                       alpha_rec_HII_net(T*exp(-dlnT_capture),y_HI),        &
                       alpha_rec_HII_net(T,y_HI),                           &
                       alpha_rec_HII_net(T*exp(dlnT_capture),y_HI))
   end function lambda_rec_HII

   ! He II, with the net coefficient of alpha_rec_HeII_net. With the
   ! dielectronic term of alphaB_HeII_new (Badnell) the relation also
   ! charges the resonance energy each dielectronic capture removes, which
   ! dominates above ~4e4 K (1.9 kT per capture at 5e4 K, MEASURED).
   elemental double precision function lambda_rec_HeII(T,y_gnd)
   real*8, intent(in) :: T, y_gnd
   lambda_rec_HeII = capture_energy_loss_rate(T,                           &
                        alpha_rec_HeII_net(T*exp(-dlnT_capture),y_gnd),     &
                        alpha_rec_HeII_net(T,y_gnd),                        &
                        alpha_rec_HeII_net(T*exp(dlnT_capture),y_gnd))
   end function lambda_rec_HeII

   ! He III. NB the case-B He III recombination-cooling fit printed by Hui
   ! & Gnedin (1997, Appendix A, p. 41), 8 x 3.435e-30 T lambda^1.970/
   ! [1 + (lambda/2.250)^0.376]^3.720, multiplies the H II fit by
   ! Z^3 = 8 but keeps T where hydrogenic scaling, Lambda_Z(T) =
   ! Z^3 Lambda_H(T/Z^2) (Hummer 1994, eq. 10 for alpha and beta), replaces
   ! it by T/Z^2: the printed form is 4x too large and removes 3.0 kT per
   ! capture at 1e4 K, above the (3/2 + dln alpha/dln T) kT < (3/2) kT a
   ! purely radiative capture (dln alpha/dln T < 0) can remove. With the
   ! factor 2 = Z^3/Z^2 that fit agrees with the relation applied to the
   ! default balance coefficient to 2.8% over 2e3-1e5 K. VALIDITY: at
   ! y_HeII = 0 and with the default coefficient the result agrees with
   ! 2 k T beta_B(T/4) of Hummer (1994, Table 1) to 1.7% over 4e3-1.6e5 K
   ! and 2.0% at 2.5e5 K (MEASURED).
   elemental double precision function lambda_rec_HeIII(T,y_HeII)
   real*8, intent(in) :: T, y_HeII
   lambda_rec_HeIII = capture_energy_loss_rate(T,                          &
                         alpha_rec_HeIII_net(T*exp(-dlnT_capture),y_HeII),  &
                         alpha_rec_HeIII_net(T,y_HeII),                     &
                         alpha_rec_HeIII_net(T*exp(dlnT_capture),y_HeII))
   end function lambda_rec_HeIII

   !--------------!

   ! Collisional ionization rate for HI
   elemental double precision function ci_rate_HI(T)
   real*8, intent(in) :: T
   real*8 :: th

   if (atomic_rate_set_k22) then
      ! Koskinen et al. 2022, Table 1 R3, which is the Voronov (1997) fit
      ! with the parameters ci_HI_new evaluates by default; the branch is
      ! what pins it when legacy_hhe_rates would otherwise take the Abel et
      ! al. fit below.
      ci_rate_HI = ci_HI_new(T)
   else if (legacy_hhe_rates) then
      ! Abel, Anninos, Zhang & Norman 1997, New Astronomy 2, 181
      !   (fit to Janev et al. 1987)
      th = log(T*8.61733d-5)
      ci_rate_HI = exp(-3.271396786d1 + 1.35365560d1*th &
               -5.73932875d0*th**2.0d0    + 1.56315498d0*th**3.0d0     &
               -2.87705600d-1*th**4.0d0 + 3.48255977d-2*th**5.0d0    &
               -2.63197617d-3*th**6.0d0 + 1.11954395d-4*th**7.0d0    &
               -2.03914985d-6*th**8.0d0)
   else
      ci_rate_HI = ci_HI_new(T)
   endif

   end function ci_rate_HI

   !--------------!

   ! Collisional ionization rate for HeI
   elemental double precision function ci_rate_HeI(T)
   real*8, intent(in) :: T
   real*8 :: th

   if (atomic_rate_set_k22) then
      ! Koskinen et al. 2022, Table 1 R4, again the Voronov (1997) fit
      ! ci_HeI_new evaluates by default (see ci_rate_HI).
      ci_rate_HeI = ci_HeI_new(T)
   else if (legacy_hhe_rates) then
      ! Abel, Anninos, Zhang & Norman 1997, New Astronomy 2, 181
      !   (fit to Janev et al. 1987)
      th = log(T*8.61733d-5)
      ci_rate_HeI = exp(-4.409864886d1 + 2.391596563d1*th &
            -1.07532302d1*th**2.0d0 + 3.05803875d0*th**3.0d0         &
            -5.6851189d-1*th**4.0d0 + 6.79539123d-2*th**5.0d0        &
            -5.0090561d-3*th**6.0d0 + 2.06723616d-4*th**7.0d0        &
            -3.64916141d-6*th**8.0d0)
   else
      ci_rate_HeI = ci_HeI_new(T)
   endif

   end function ci_rate_HeI

   !--------------!

   ! Collisional ionization rate for HeII
   elemental double precision function ci_rate_HeII(T)
   real*8, intent(in) :: T
   real*8 :: xl

   if (legacy_hhe_rates) then
      ! Hui & Gnedin 1997, MNRAS 292, 27, collisional-ionization fit
      xl = 2.0d0*631515.0d0/T
      ci_rate_HeII = 19.95d0*exp(-xl/2.0d0)*T**(-1.5d0)*                &
                           xl**(-1.089d0)/(1.0d0+(xl/0.553d0)**0.735d0)**1.275d0
   else
      ci_rate_HeII = ci_HeII_new(T)
   endif

   end function ci_rate_HeII


   !--------------!

   !--- Electron-impact excitation of H I ------------------------------!
   ! ONE rate set for the 1s -> 2s and 1s -> 2p excitation, used by the
   ! H I cooling below and by the H(n=2) population model
   ! (hydrogen_n2_rates, which also feeds the Ly-alpha transfer).
   ! Effective collision strengths: CHIANTI v11.0.2 h_1.scups, whose header
   ! cites Anderson, Ballance, Badnell & Summers (2000, J. Phys. B 33,
   ! 1255, "revised 2002"; R-matrix with 15 physical terms and 24
   ! pseudostates), with the 2p fine-structure levels summed and the five
   ! scaled knot values of each transition descaled with the natural cubic
   ! spline of CHIANTI's DESCALE_SCUPS.PRO. The paper (Table 2, nl -> n'l',
   ! n <= 5) tabulates Upsilon at Te = 0.5, 1, 3, 5, 10, 15, 20 and 25 eV
   ! only, 5.80e3-2.90e5 K; below 5.80e3 K the CHIANTI value is the scaled
   ! spline running to its threshold knot, which the paper does not
   ! constrain. CHIANTI / Table 2 (MEASURED): 1s-2s 0.993-1.012 and 1s-2p
   ! 0.975-1.017 at 0.5-5 eV; above 5 eV CHIANTI is higher (1s-2p 1.04 at
   ! 10 eV, 1.31 at 25 eV), its dipole fits ending at the Bethe limit
   ! 4gf/dE while the table's 1s-np values turn over; the source of the
   ! 2002 revision is not identified in the paper or in CHIANTI. log10
   ! Upsilon is a quartic in log10(T/1e4 K) fitted to CHIANTI over 1e3-1e5 K
   ! to 2.2% (1s-2s) and 0.8% (1s-2p); quartic / Table 2 at 0.5-5 eV
   ! (MEASURED): 0.991-1.015 (1s-2s), 0.976-1.015 (1s-2p)
   ! (cooling_data/fit_hydrogen_helium_excitation.py). VALIDITY: 1e3-1e5 K;
   ! T is held inside that range (upsilon_quartic_logT), so above 1e5 K the
   ! held 1s-2p value is low (Table 2: 1.75 at 10 eV against the held
   ! 1.66). Energy: the observed 1s-2s interval, 82258.954 cm^-1 =
   ! 10.19881 eV (NIST ASD, as in CHIANTI
   ! h_1.elvlc; the 2p levels lie within 0.37 cm^-1 of it). Weights g(1s) = 2,
   ! g(2s) = 2, g(2p) = 6.
   elemental double precision function upsilon_HI_1s2s(T)
   real*8, intent(in) :: T
   upsilon_HI_1s2s = upsilon_quartic_logT(-5.41913497d-01, 1.45680880d-01, &
                        -5.14786365d-02, -4.89245376d-02, 4.41143616d-02, T)
   end function upsilon_HI_1s2s

   elemental double precision function upsilon_HI_1s2p(T)
   real*8, intent(in) :: T
   upsilon_HI_1s2p = upsilon_quartic_logT(-3.02352097d-01, 3.07054967d-01, &
                        2.18544412d-01, 2.79631531d-02, -3.12591467d-02, T)
   end function upsilon_HI_1s2p

   elemental double precision function excitation_rate_HI_1s2s(T)
   real*8, intent(in) :: T
   excitation_rate_HI_1s2s = coll_rate_prefactor/(2.0d0*sqrt(T))          &
                           *upsilon_HI_1s2s(T)*exp(-E_HI_n2_eV/(kb_eV*T))
   end function excitation_rate_HI_1s2s

   elemental double precision function excitation_rate_HI_1s2p(T)
   real*8, intent(in) :: T
   excitation_rate_HI_1s2p = coll_rate_prefactor/(2.0d0*sqrt(T))          &
                           *upsilon_HI_1s2p(T)*exp(-E_HI_n2_eV/(kb_eV*T))
   end function excitation_rate_HI_1s2p

   ! De-excitation 2s -> 1s and 2p -> 1s, the detailed-balance reverses of
   ! the two rates above with the Boltzmann factor cancelled analytically.
   elemental double precision function deexcitation_rate_HI_2s1s(T)
   real*8, intent(in) :: T
   deexcitation_rate_HI_2s1s = coll_rate_prefactor/(2.0d0*sqrt(T))        &
                             *upsilon_HI_1s2s(T)
   end function deexcitation_rate_HI_2s1s

   elemental double precision function deexcitation_rate_HI_2p1s(T)
   real*8, intent(in) :: T
   deexcitation_rate_HI_2p1s = coll_rate_prefactor/(6.0d0*sqrt(T))        &
                             *upsilon_HI_1s2p(T)
   end function deexcitation_rate_HI_2p1s

   ! Collisional-excitation cooling of H I [erg cm^3 s^-1 per (n_e n_HI)],
   ! coronal limit: 10.2 eV times the two n = 2 rates above, plus every
   ! upper level with n >= 3 (CHIANTI h_1, 3s ... 5g, observed energies),
   ! fitted as C T^-1/2 (a + ln(1 + T/T0)) exp(-Tex/T) over 3e3-1e5 K to
   ! 0.983-1.009 (cooling_data/fit_hydrogen_helium_excitation.py), the
   ! knot values descaled with the natural cubic spline of CHIANTI's
   ! DESCALE_SCUPS.PRO. The n >= 3 part is 0.6% of the total at 5e3 K,
   ! 5.4% at 1e4 K and 15% at 2e4 K; the n >= 3 energy-weighted CHIANTI
   ! sum is 1.003-1.012 of Anderson et al. (2000, Table 2) at 0.5-5 eV and
   ! 1.07 at 10 eV (MEASURED). The Black (1981, Table 3) form
   ! with the factor of Cen (1992), 7.5e-19/(1 + (T/1e5)^1/2)
   ! exp(-118348/T), is 0.90, 0.97 and 0.91 of this CHIANTI sum at
   ! 5e3, 1e4 and 2e4 K (MEASURED). Collisional de-excitation of n = 2
   ! returns 10.2 eV to the
   ! electrons; with the H(n=2) model on that return is its own heating
   ! channel (excited_hydrogen: Hdx), from the same rates.
   elemental double precision function lambda_coex_HI(T)
   real*8, intent(in) :: T
   lambda_coex_HI = E_HI_n2_eV/erg2eV                                     &
                    *(excitation_rate_HI_1s2s(T) + excitation_rate_HI_1s2p(T)) &
                  + 3.23738685d-17/sqrt(T)                                 &
                    *(5.70053703d-01 + log(1.0d0 + T/2.47570408d4))       &
                    *exp(-140671.201d0/T)
   end function lambda_coex_HI

   !--------------!

   ! Collisional-excitation cooling of He I out of the 1^1S ground level
   ! [erg cm^3 s^-1 per (n_e n_1^1S)], coronal limit (every excitation
   ! followed by a radiative decay), summed over the 28 upper terms n <= 5
   ! of Bray, Burgess, Fursa & Tully (2000, A&AS 146, 481, Table 2), with
   ! the strengths the metastable network uses (the table inside
   ! 10^3.75-10^5.75 K, the CHIANTI v11 he_1 temperature dependence scaled
   ! to it outside; see upsilon_HeI_11S_23S) and observed term energies,
   ! fitted as C T^-1/2 (a + ln(1 + T/T0)) exp(-Tex/T)
   ! (cooling_data/fit_hydrogen_helium_excitation.py). The CHIANTI he_1 sum
   ! (Sawey & Berrington 1993 and Bray et al. 2000 in five-knot fits) is
   ! 1.067, 1.062, 1.063 and 1.095 of it at 5e3, 1e4, 2e4 and 5e4 K
   ! (MEASURED): CHIANTI carries 1^1S - 2^3S 7% above the table and several
   ! n = 3, 4 singlets up to 2-3 times above it, the n = 4 excess being
   ! the one Bray et al. (Sect. 6) find against the R-matrix data of Sawey
   ! & Berrington.
   !
   ! lambda_coex_HeI_all_levels: every upper term; the cooling of ground
   ! state helium when the 2^3S metastable is not a level of the network.
   ! Fit / sum 0.992-1.005 over 3e3-1e5 K. The coefficient
   ! 1.1e-19 T^0.082 exp(-2.3e5/T) of the ATES code (Caldiroli et al. 2021,
   ! which gives no source for it; it is not in Black 1981, Table 3) is
   ! 0.79, 0.85 and 0.95 of this sum at 3e3, 5e3 and 1e4 K (MEASURED).
   !
   ! lambda_coex_HeI_except_23S: every upper term except 2^3S. With the
   ! metastable tracked, 1^1S -> 2^3S is the network's own excitation
   ! rate q13, and its energy is charged in radiative_cooling_of_cell with
   ! q13 and its detailed-balance reverse q31g; this sum leaves that one
   ! transition out so that the 19.82 eV of an excitation is charged once.
   ! Fit / sum 0.981-1.009 over 3e3-1e5 K. The 2^3S channel is 97.4% of
   ! the full sum at 3e3 K, 72.0% at 1e4 K and 50.2% at 2e4 K (MEASURED).
   ! Excitations into the higher triplet levels (2^3P ...) cascade into
   ! 2^3S: HeITR_coeffs (util_ion_eq) adds them to q13 through
   ! excitation_rate_HeI_11S_triplets, and their energy above 2^3S is the
   ! part of this sum that is radiated.
   elemental double precision function lambda_coex_HeI_all_levels(T)
   real*8, intent(in) :: T
   lambda_coex_HeI_all_levels = 8.18542388d-17/sqrt(T)                    &
                              *(2.04216400d-01 + log(1.0d0 + T/7.40245815d4)) &
                              *exp(-230931.918d0/T)
   end function lambda_coex_HeI_all_levels

   elemental double precision function lambda_coex_HeI_except_23S(T)
   real*8, intent(in) :: T
   lambda_coex_HeI_except_23S = 7.23126663d-17/sqrt(T)                    &
                              *(1.20064198d-01 + log(1.0d0 + T/6.16209136d4)) &
                              *exp(-240448.382d0/T)
   end function lambda_coex_HeI_except_23S

   !--------------!

   ! Collisional-excitation cooling of He II [erg cm^3 s^-1 per
   ! (n_e n_He+)], coronal limit, summed over the upper levels of the
   ! CHIANTI v11 he_2 model (Ballance 2003, R-matrix with pseudostates),
   ! observed energies, knot values descaled with the natural cubic spline
   ! of CHIANTI's DESCALE_SCUPS.PRO, fitted as C T^-1/2 (a + ln(1 + T/T0))
   ! exp(-Tex/T): 0.995-1.008 over 5e3-2e5 K (cooling_data/
   ! fit_hydrogen_helium_excitation.py). The Black (1981, Table 3) form
   ! with the factor of Cen (1992), 5.54e-17 T^-0.397/(1 + (T/1e5)^1/2)
   ! exp(-473638/T), is 0.77, 0.73 and 0.46 of this sum at 1e4, 2e4 and
   ! 1e5 K (MEASURED).
   elemental double precision function lambda_coex_HeII(T)
   real*8, intent(in) :: T
   lambda_coex_HeII = 2.20200522d-16/sqrt(T)                              &
                    *(5.87104891d-01 + log(1.0d0 + T/2.83620056d5))       &
                    *exp(-473383.427d0/T)
   end function lambda_coex_HeII

   !--------------!

   ! Legacy metal line cooling ("cno_cool 0" in metals.inp): the AIOLOS
   ! fits of C I, C II, O I and O II (see the provenance note at the top
   ! of the metal cooling section). The default CHIANTI coefficients are
   ! cool_CI_ne_func ... cool_OII_ne_func.
   elemental double precision function lambda_line_CI(T)
   real*8, intent(in) :: T
   lambda_line_CI = 1.0d-24 + 3.1d-20*exp(-15162.0d0/T)                   &
                          *(1.0d0 + (T/2.0d4)**1.5d0)
   end function lambda_line_CI

   elemental double precision function lambda_line_CII(T)
   real*8, intent(in) :: T
   lambda_line_CII = 1.5d-23 + 3.1d-20*exp(-45162.0d0/T)                  &
                           *(1.0d0 + (T/0.75d4)**1.5d0)
   end function lambda_line_CII

   elemental double precision function lambda_line_OI(T)
   real*8, intent(in) :: T
   lambda_line_OI = 5.5d-24 + 1.1d-20*exp(-30162.0d0/T)                   &
                          *(1.0d0 + (T/0.75d4)**0.5d0)
   end function lambda_line_OI

   elemental double precision function lambda_line_OII(T)
   real*8, intent(in) :: T
   lambda_line_OII = 5.1d-20*exp(-35162.0d0/T)                            &
                         *(1.0d0 + (T/0.75d4)**0.5d0)
   end function lambda_line_OII

   !--------------!

   ! Grid forms of the metal rates, by the ion's canonical species_table
   ! index (im_* constants); the case lists are those of
   ! metal_recombination_coefficient and
   ! metal_collisional_ionization_coefficient.

   ! Recombination rate of the (recombining) ion: ion + e -> ion-1.
   subroutine rec_coeff_by_ion_range(i,T,out,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   integer, intent(in) :: i
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: out
   integer :: j
   do j = j_lo,j_hi
      out(j) = metal_recombination_coefficient(i, T(j))
   enddo
   end subroutine rec_coeff_by_ion_range


   !--------------!

   ! Collisional ionization rate of the ion: ion -> ion+1.
   subroutine ion_coeff_by_ion_range(i,T,out,j_lo,j_hi)
   integer, intent(in) :: j_lo,j_hi
   integer, intent(in) :: i
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
   real*8, dimension(1-Ng:N+Ng), intent(inout) :: out
   integer :: j
   do j = j_lo,j_hi
      out(j) = metal_collisional_ionization_coefficient(i, T(j))
   enddo
   end subroutine ion_coeff_by_ion_range


   !--------------!

   ! Radiative line-cooling coefficient of one metal ion [erg cm^3 s^-1 per
   ! (n_e n_ion)] at the temperature Ts, electron density nes and H I
   ! density nHIs of the cell, with the escape probabilities beta and the
   ! incident-field photon occupation numbers nbar of the eight ground-term
   ! fine-structure lines (1 and 0 in the optically thin vacuum limit; see
   ! fine_structure_line_transfer). THE CASE LIST LIVES HERE ONLY; its one
   ! caller is radiative_cooling_of_cell. C/N/O take the statistical-
   ! equilibrium coefficients under cno_chianti (the default) and the
   ! legacy AIOLOS fits otherwise, with no nitrogen cooling in the legacy
   ! set; Mg II and Na I their coronal fits; Ca II and Mg I their (T, n_e)
   ! tables, which carry every channel of the ion; Fe I its table; Fe II
   ! its (T, n_e) table. Ions with no entry return 0.
   pure double precision function metal_line_cooling_coefficient          &
                                    (i,Ts,nes,nHIs,beta,nbar) result(c)
   integer, intent(in) :: i
   real*8, intent(in)  :: Ts, nes, nHIs
   real*8, intent(in)  :: beta(n_fsline), nbar(n_fsline)
   select case (i)
      case (im_CI)
         if (cno_chianti) then
            c = cool_CI_ne_func(Ts, nes, nHIs, beta(ifs_CI609),            &
                                beta(ifs_CI370), nbar(ifs_CI609),          &
                                nbar(ifs_CI370))
         else
            c = lambda_line_CI(Ts)
         endif
      case (im_CII)
         if (cno_chianti) then
            c = cool_CII_ne_func(Ts, nes, nHIs, beta(ifs_CII158),          &
                                 nbar(ifs_CII158))
         else
            c = lambda_line_CII(Ts)
         endif
      case (im_NI)
         if (cno_chianti) then
            c = cool_NI_ne_func(Ts, nes)
         else
            !To Be Checked/AIOLOS tuning?
            c = 0.0d0
         endif
      case (im_NII)
         if (cno_chianti) then
            c = cool_NII_ne_func(Ts, nes, nHIs, beta(ifs_NII205),          &
                                 beta(ifs_NII122), nbar(ifs_NII205),       &
                                 nbar(ifs_NII122))
         else
            !To Be Checked/AIOLOS tuning?
            c = 0.0d0
         endif
      case (im_OI)
         if (cno_chianti) then
            c = cool_OI_ne_func(Ts, nes, nHIs, beta(ifs_OI63),             &
                                beta(ifs_OI145), beta(ifs_OI44),           &
                                nbar(ifs_OI63), nbar(ifs_OI145),           &
                                nbar(ifs_OI44))
         else
            c = lambda_line_OI(Ts)
         endif
      case (im_OII)
         if (cno_chianti) then
            c = cool_OII_ne_func(Ts, nes)
         else
            c = lambda_line_OII(Ts)
         endif
      case (im_MgI);  c = cool_MgI_ne_func (Ts, nes)
      case (im_MgII); c = cool_MgII_func(Ts)
      case (im_CaII); c = cool_CaII_ne_func(Ts, nes)
      case (im_NaI);  c = cool_NaI_func (Ts)
      case (im_FeI);  c = cool_table_value(cool_logL_FeI, Ts)
      case (im_FeII); c = cool_FeII_ne_value(Ts, nes)
      case default;   c = 0.0d0
   end select
   end function metal_line_cooling_coefficient

   !--------------!

   ! THE RADIATIVE COOLING OF ONE CELL of atomic and ionized gas, channel
   ! by channel [erg cm^-3 s^-1]. This is the ONE assembly: eval_cool
   ! (utils_ion_eq) calls it for every cell of the grid, and the
   ! cell-by-cell temperature root of the advection post-process
   ! (T_equation) calls it for its trial temperature, so the temperature
   ! that root returns balances the cooling eval_cool reports at it.
   !
   ! Inputs: the temperature Ts [K], the electron density nes (the caller's
   ! charge sum), the densities of H I, H II, the He I ground singlet
   ! n(1^1S), the He 2^3S metastable (zero when it is not tracked), He II,
   ! He III and the metal ions nm [cm^-3], the frozen ground-capture escape
   ! weights of the three recombinations (y_HI for H II -> H I, y_gnd for
   ! He II -> He I, y_HeII for He III -> He II; alpha_rec_HII_net,
   ! alpha_rec_HeII_net, alpha_rec_HeIII_net), and the fine-structure line
   ! transfer of the cell (beta_fs, nbar_fs).
   !
   ! Channels (the columns of eval_cool's breakdown, n_cool_chan in all):
   !   1  recombination: the captured electrons' kinetic energy, of H II,
   !      He II, He III and of every recombining metal ion
   !   2  collisional ionization: the ionization potential of each level,
   !      e_th_*_erg for H and He and mion_ethr for the metal ions, the
   !      same energies the photoelectrons are charged against
   !   3  collisional excitation of H I
   !   4  collisional excitation of He I: the ground singlet, and with the
   !      metastable tracked its own channels -- the net 1^1S <-> 2^3S
   !      exchange E(2^3S) (q13 n(1^1S) - q31g n(2^3S)), negative where the
   !      metastable is overpopulated against its Boltzmann ratio to the
   !      ground (the superelastic collisions then heat the electrons), the
   !      excitation within the triplets that radiates (10830 A and the
   !      higher triplet lines), and the threshold energies of the
   !      2^3S -> 2^1S, 2^1P and 2^3S -> n^1L (n >= 3) excitations
   !   5  collisional excitation of He II
   !   6  free-free emission of every ion by its net charge
   !   7-10  the molecular channels, zero here (eval_cool adds them)
   !   10+i  line cooling of metal ion i
   !
   ! THE METAL IONIZATION ENERGY. A photoionization of a metal ion leaves
   ! h nu - I in the gas and stores I in the ionization; a collisional
   ! ionization takes I out of the electron gas (channel 2), and a
   ! recombination returns I in the photon it emits while the electron gas
   ! loses only the captured electron's kinetic energy (channel 1, the
   ! capture relation applied to metal_recombination_coefficient, the
   ! coefficient the metal balance uses, radiative and dielectronic alike:
   ! metal_recombination_energy_rate). With both terms the electron
   ! gas, the ionization energy and the radiation close: the physics probe
   ! metal_ionization_energy_ledger checks that ledger.
   pure subroutine radiative_cooling_of_cell(Ts, nes, nhi, nhii, nheiS,     &
                                             nheiTR, nheii, nheiii,         &
                                             y_HI, y_gnd, y_HeII,           &
                                             nm, beta_fs, nbar_fs, chan)
   real*8, intent(in)  :: Ts, nes, nhi, nhii, nheiS, nheiTR, nheii, nheiii
   real*8, intent(in)  :: y_HI, y_gnd, y_HeII
   real*8, intent(in)  :: nm(n_mion), beta_fs(n_fsline), nbar_fs(n_fsline)
   real*8, intent(out) :: chan(n_cool_chan)
   real*8  :: gz1, gz2, charge_sum, coex_HeI, rec_metal, ion_metal
   integer :: i

   chan = 0.0d0

   chan(1) = nes*( lambda_rec_HII(Ts,y_HI)*nhii                            &
                 + lambda_rec_HeII(Ts,y_gnd)*nheii                         &
                 + lambda_rec_HeIII(Ts,y_HeII)*nheiii )

   chan(2) = nes*( e_th_HI_erg*ci_rate_HI(Ts)*nhi                          &
                 + e_th_HeI_erg*ci_rate_HeI(Ts)*nheiS                      &
                 + e_th_HeII_erg*ci_rate_HeII(Ts)*nheii                    &
                 + e_th_HeTR_erg*ionization_rate_HeI_23S(Ts)*nheiTR )

   ! The metal recombination and collisional-ionization energies.
   if (thereis_metals) then
      rec_metal = 0.0d0
      ion_metal = 0.0d0
      do i = 1, n_mion
         if (nm(i) .eq. 0.0d0) cycle
         if (mion_stage(i) .gt. 0) rec_metal = rec_metal + nm(i)          &
              *metal_recombination_energy_rate(i, Ts)
         if (mion_ethr(i) .gt. 0.0d0) ion_metal = ion_metal + nm(i)       &
              *metal_collisional_ionization_coefficient(i, Ts)           &
              *mion_ethr(i)/erg2eV
      enddo
      chan(1) = chan(1) + nes*rec_metal
      chan(2) = chan(2) + nes*ion_metal
   endif

   chan(3) = nes*lambda_coex_HI(Ts)*nhi

   if (thereis_HeITR) then
      coex_HeI = lambda_coex_HeI_except_23S(Ts)*nheiS                      &
               + E_HeI_23S_eV/erg2eV                                       &
                 *( excitation_rate_HeI_11S_23S(Ts)*nheiS                  &
                  - deexcitation_rate_HeI_23S_11S(Ts)*nheiTR )             &
               + ( cooling_HeI_23S_triplets(Ts)                            &
                 + cooling_HeI_23S_singlets_n3(Ts)                         &
                 + ( E_HeI_21S_23S_eV*excitation_rate_HeI_23S_21S(Ts)      &
                   + E_HeI_21P_23S_eV*excitation_rate_HeI_23S_21P(Ts) )    &
                   /erg2eV )*nheiTR
   else
      coex_HeI = lambda_coex_HeI_all_levels(Ts)*nheiS
   endif
   chan(4) = nes*coex_HeI

   chan(5) = nes*lambda_coex_HeII(Ts)*nheii

   ! Free-free scales with the ion NET charge Z_ion: H II, He II and
   ! singly ionized metals are Z_ion = 1, He III and doubly ionized metals
   ! Z_ion = 2 (mion_z2 = stage^2, zero for a neutral), with the Gaunt
   ! factor at that charge. The molecular ions of the cold base are left
   ! out: this hot-plasma expression is an extrapolation there, and
   ! negligible against the infrared coolants.
   gz1 = gbar_ff(Ts, 1.0d0)
   gz2 = gbar_ff(Ts, 2.0d0)
   charge_sum = gz1*nhii + gz1*nheii + 4.0d0*gz2*nheiii
   if (thereis_metals) then
      do i = 1, n_mion
         if (mion_stage(i) .eq. 2) then
            charge_sum = charge_sum + dble(mion_z2(i))*gz2*nm(i)
         else
            charge_sum = charge_sum + dble(mion_z2(i))*gz1*nm(i)
         endif
      enddo
   endif
   chan(6) = nes*1.426d-27*sqrt(Ts)*charge_sum

   if (thereis_metals) then
      do i = 1, n_mion
         if (.not. mion_iscool(i)) cycle
         chan(10+i) = nes*nm(i)*metal_line_cooling_coefficient(i, Ts, nes, &
                                          nhi, beta_fs, nbar_fs)
      enddo
   endif

   end subroutine radiative_cooling_of_cell

   !--------------!


   !--------------!

   ! Bilinear value of a 2-D (T, ne) cooling table at one cell. This is the
   ! ONE definition (cool_FeII_ne_value). The (log T, log n_e) axes are
   ! cool_logT (3.0..5.0, dlogT) and cool_logne (0..14, dlogne = 0.5). T
   ! clamps to the nearest edge; below the lowest density the coefficient
   ! is held (the coronal limit, Lambda per n_e n_ion independent of n_e),
   ! and above the highest it is continued by density_saturation_log.
   pure double precision function cool_table_value_2d(logL2d,Ts,nes)
   real*8, dimension(NCOOLT,NCOOLNE), intent(in) :: logL2d
   real*8, intent(in) :: Ts, nes
   real*8, parameter :: dlogne = 0.5d0
   integer :: kt, ke
   real*8  :: lt, ln, post, posn, ft, fn
   real*8  :: f00, f10, f01, f11, l0, l1
   lt   = log10(max(Ts, 1.0d0))
   post = (lt - cool_logT(1))/cool_dlogT + 1.0d0
   if (.not. (post > 1.0d0)) then
      kt = 1;          ft = 0.0d0
   else if (post .ge. dble(NCOOLT)) then
      kt = NCOOLT - 1; ft = 1.0d0
   else
      kt = int(post);  ft = post - dble(kt)
   endif
   ln   = log10(max(nes, 1.0d0))
   posn = (ln - cool_logne(1))/dlogne + 1.0d0
   if (.not. (posn > 1.0d0)) then
      ke = 1;          fn = 0.0d0
   else if (posn .ge. dble(NCOOLNE)) then
      ke = NCOOLNE - 1; fn = 1.0d0
   else
      ke = int(posn);  fn = posn - dble(ke)
   endif
   f00 = logL2d(kt,   ke  )
   f10 = logL2d(kt+1, ke  )
   f01 = logL2d(kt,   ke+1)
   f11 = logL2d(kt+1, ke+1)
   l0  = f00 + ft*(f10 - f00)
   l1  = f01 + ft*(f11 - f01)
   if (posn > dble(NCOOLNE)) then
      cool_table_value_2d = 10.0d0**density_saturation_log(l0, l1,     &
                                        ln - cool_logne(NCOOLNE))       &
                                 *coronal_excitation_cutoff(Ts)
   else
      cool_table_value_2d = 10.0d0**( l0 + fn*(l1 - l0) )           &
                                 *coronal_excitation_cutoff(Ts)
   endif
   end function cool_table_value_2d

   !--------------!

   !--- Density-resolved line cooling above the ground term ------------!
   ! (C/N/O below; Ca II and Mg I, whose tables carry the whole ion, are
   ! described at cool_logLrem_CaII.)
   ! The closed-form fits for C I, C II, N I, N II, O I and O II carry the
   ! channels LEAVING the ground term in their coronal limit: every
   ! collisional excitation is taken to radiate. That is the n_e -> 0
   ! limit of statistical equilibrium, and it overstates any line whose
   ! critical density is at or below the electron density of the gas. The
   ! metastable levels of these six ions have n_cr = 1e4 to 1e9 cm^-3
   ! (2p2 1D2/1S0 for C I and N II, 2p4 1D2/1S0 for O I, 2p3 2D*/2P* for
   ! N I and O II, 2s2p2 4P for C II), at or below a wind base, so the
   ! coronal form is outside its density domain there.
   !
   ! What replaces it is the statistical equilibrium of the same channels
   ! at the gas's own electron density, tabulated in (log T, log n_e):
   !   Lambda_rem(T,n_e) = ( sum_{u outside the ground term, l}
   !                          f_u(T,n_e) A_ul E_ul ) / n_e
   ! with f_u the level populations of the full CHIANTI model of the ion.
   ! The ground term is EXCLUDED from the sum because the coefficients
   ! solve it exactly themselves, with the neutral-hydrogen de-excitation
   ! channel and the line escape probability that no such table can carry;
   ! including it here would double count it.
   !
   ! The floor column of the table is the coronal limit of the same
   ! channels and is therefore directly comparable with the fit it
   ! replaces. They do NOT agree at low temperature, and the reason is a
   ! defect in the fits: they were built from the .scups THEORETICAL
   ! transition energies while the populations use the observed level
   ! energies, and the two differ by up to 14 percent in dE for the
   ! O II 4S* - 2D* excitation, which at 2000 K is a factor 14 in
   ! exp(-dE/kT). The ratio is 0.96 to 1.6 for C I, 0.89 to 1.0 for C II,
   ! 1.0 to 2.3 for N I, 1.0 to 4.2 for N II, 0.91 to 1.0 for O I and 1.1
   ! to 220 for O II over 1e3 to 2e4 K (MEASURED). The fits are retained as the
   ! coronal reference and for the legacy (cno_cool 0) branch.
   !
   ! Ca II and Mg I have a single ground level, so their tables carry the
   ! whole ion: the resonance line and the channels through the Ca II 3d 2D
   ! and Mg I 3s3p 3P levels, whose critical densities (6e4 to 5e9 cm^-3)
   ! lie inside the range of a wind base. Mg II and Na I have no table:
   ! every excited level of Mg II decays by a permitted line, and its full
   ! statistical equilibrium equals the coronal sum to 0.2 per cent for
   ! n_e <= 1e12 cm^-3; Na I has no CHIANTI model atom. Fe I has no CHIANTI
   ! model atom either; Fe II carries its own two-dimensional table
   ! (cool_logL_FeII_ne above).
   !
   ! Taking the cooling from the level populations at the local electron
   ! density is the method of MoCHII (~/MoCHII/MoCHII_v1.00,
   ! md/COOLING_LOCAL_NE_PLAN.md, the same author's photoionization code),
   ! and the case for a precomputed (T,n_e) table rather than a live solve
   ! is that document's design point D1. Tables: CHIANTI v11.0.2, the
   ! statistical equilibrium of cooling_data/multilevel_statistical_
   ! equilibrium.py with collision strengths descaled by the natural cubic
   ! spline of CHIANTI's DESCALE_SCUPS.PRO, built by cooling_data/
   ! metal_cooling_density_resolved.py on the cool_logT / cool_logne axes.

   ! Bilinear interpolation of log10 Lambda_rem in (log10 T, log10 n_e) on
   ! the cool_logT / cool_logne axes, with edge clamping.
   pure double precision function remainder_table_value(logL,Ts,nes)
   real*8, dimension(NCOOLT,NCOOLNE), intent(in) :: logL
   real*8, intent(in) :: Ts, nes
   real*8, parameter :: dlogne = 0.5d0
   integer :: kt, ke
   real*8  :: post, posn, ft, fn, l0, l1
   post = (log10(max(Ts, 1.0d0)) - cool_logT(1))/cool_dlogT + 1.0d0
   if (.not. (post > 1.0d0)) then
      kt = 1;          ft = 0.0d0
   else if (post .ge. dble(NCOOLT)) then
      kt = NCOOLT - 1; ft = 1.0d0
   else
      kt = int(post);  ft = post - dble(kt)
   endif
   posn = (log10(max(nes, 1.0d0)) - cool_logne(1))/dlogne + 1.0d0
   if (.not. (posn > 1.0d0)) then
      ke = 1;           fn = 0.0d0
   else if (posn .ge. dble(NCOOLNE)) then
      ke = NCOOLNE - 1; fn = 1.0d0
   else
      ke = int(posn);   fn = posn - dble(ke)
   endif
   l0 = logL(kt,ke)   + ft*(logL(kt+1,ke)   - logL(kt,ke))
   l1 = logL(kt,ke+1) + ft*(logL(kt+1,ke+1) - logL(kt,ke+1))
   if (posn > dble(NCOOLNE)) then
      remainder_table_value = 10.0d0**density_saturation_log(l0, l1,   &
                   log10(max(nes, 1.0d0)) - cool_logne(NCOOLNE))
   else
      remainder_table_value = 10.0d0**( l0 + fn*(l1 - l0) )
   endif
   end function remainder_table_value

   ! log10 Lambda ABOVE THE TOP DENSITY OF A (T, n_e) TABLE, dln (n_e
   ! beyond the top column, in dex) past it, from the last two columns
   ! l_below, l_top (log10 Lambda at log n_e = 13.5 and 14). Lambda is per
   ! (n_e n_ion): a channel whose upper level is thermalized (n_e above its
   ! critical density) emits a fixed power per ion, so its Lambda falls as
   ! 1/n_e (slope -1), and a channel still in its coronal regime has a
   ! Lambda independent of n_e (slope 0). A sum of such channels has a
   ! slope in [-1, 0], and the slope of the last interval, clipped to that
   ! range, is carried on. This is exact where the top column is fully
   ! thermalized or fully coronal; for a mixture the fraction still
   ! coronal can only grow with n_e, so the continuation errs toward too
   ! little cooling at densities far beyond the table. Holding Lambda, the
   ! former clamp, made the volume emission of thermalized channels keep
   ! growing as n_e.
   pure double precision function density_saturation_log(l_below, l_top, dln)
   real*8, intent(in) :: l_below, l_top, dln
   real*8, parameter  :: dlogne = 0.5d0
   density_saturation_log = l_top                                       &
        + min(0.0d0, max(-1.0d0, (l_top - l_below)/dlogne))*dln
   end function density_saturation_log

   ! Line cooling of one ion above its ground term, at the local electron
   ! density [erg cm^3 s^-1 per (n_e n_ion)]. THE TABLE FOR EACH ION IS
   ! SELECTED HERE AND NOWHERE ELSE. lam_coronal is that ion's closed-form
   ! coronal remainder, returned unchanged for an ion with no table.
   pure double precision function metal_cooling_above_ground_term          &
                                    (i,Ts,nes,lam_coronal)
   integer, intent(in) :: i
   real*8, intent(in) :: Ts, nes, lam_coronal
   select case (i)
      case (im_CI)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_CI,  Ts, nes)
      case (im_CII)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_CII, Ts, nes)
      case (im_NI)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_NI,  Ts, nes)
      case (im_NII)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_NII, Ts, nes)
      case (im_OI)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_OI,  Ts, nes)
      case (im_OII)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_OII, Ts, nes)
      case (im_CaII)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_CaII, Ts, nes)
      case (im_MgI)
         metal_cooling_above_ground_term =                                &
              remainder_table_value(cool_logLrem_MgI, Ts, nes)
      case default
         metal_cooling_above_ground_term = lam_coronal
   end select
   end function metal_cooling_above_ground_term

   !--------------!

   ! Density-dependent Fe II coefficient Lambda_eff(T,ne) per (n_e n_FeII).
   ! The table is selected here and nowhere else (the one caller is
   ! metal_line_cooling_coefficient).
   pure double precision function cool_FeII_ne_value(Ts,nes)
   real*8, intent(in) :: Ts, nes
   cool_FeII_ne_value = cool_table_value_2d(cool_logL_FeII_ne, Ts, nes)
   end function cool_FeII_ne_value

   !--------------!


   ! End of module
   end module Cooling_Coefficients

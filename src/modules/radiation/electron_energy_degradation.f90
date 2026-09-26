   module electron_energy_degradation
   ! How the excess energy of a fast photoelectron is degraded in a partly
   ! ionized H / H2 / He gas: what fraction ends as heat, and what fraction
   ! ionizes each neutral target.
   !
   ! Two published determinations are used, and the division of labour
   ! between them is deliberate.
   !
   ! * Shull & van Steenberg (1985, ApJ 298, 268) supply the ENERGY BUDGET OF
   !   THE IONIZATION CHANNELS as a function of the ionized fraction x. Their
   !   fitting form C(1 - x^a)^b has the correct limit at both ends -- it
   !   vanishes at x = 1, where no neutral is left to ionize -- so it may be
   !   evaluated over the whole domain of a photoevaporative wind. Their gas
   !   is pure H + He: it has NO H2 channel, in either direction.
   !
   ! * Dalgarno, Yan & Liu (1999, ApJS 125, 237) supply THE HEATING FRACTION
   !   and everything that depends on molecular hydrogen: how the ionization
   !   energy divides between H and H2, and the energy returned as heat when
   !   electronically excited H2 dissociates. Their calculation resolves the
   !   primary energy E0 over 30-1000 eV, which Shull & van Steenberg's does
   !   not, and it does not send the heating fraction to zero in a neutral
   !   gas, which Shull & van Steenberg's fitted form does.
   !
   ! The H2 terms below are IDENTICALLY ZERO when n(H2) = 0, and the routine
   ! then reduces, expression by expression, to the two-target atomic case.
   !
   ! ---------------------------------------------------------------------
   !  VALIDITY, and where our own models sit relative to it
   ! ---------------------------------------------------------------------
   !
   !  * Ionized fraction. Dalgarno et al. state eq. (14) "for x <= 0.1" and
   !    their eq. (13) fits carry "They should not be used for x > 0.1".
   !    Shull & van Steenberg fit over 0.0001 < x < 1.0 and say so
   !    ("reproduced to ~2% over the range 0.0001 < x < 1.0"). Our molecular
   !    layer sits at x = 1e-16 to 5e-5: INSIDE Dalgarno's range and BELOW
   !    Shull & van Steenberg's. The wind above it reaches x -> 1, which is
   !    OUTSIDE Dalgarno's stated range, and the partition below
   !    extrapolates eq. (14) there deliberately -- see the validity note on
   !    dalgarno_heating_efficiency for what justifies it and what does not.
   !
   !  * Primary energy. Dalgarno et al. tabulate 30, 50, 100, 200, 500 and
   !    1000 eV; the interpolation below is linear in log10(E0) and CLAMPED
   !    at both ends. Clamping below 30 eV is inert in practice because the
   !    caller's secondary-ionization threshold is above it; clamping above
   !    1 keV freezes the coefficients at their most nearly asymptotic
   !    values, which is the behaviour the source describes at high energy.
   !    Shull & van Steenberg's Table 2 carries the footnote "These formulae
   !    are only accurate for E0 > 100 eV" and their Table 1 is a 3 keV
   !    calculation, so their fits are a single high-energy asymptote with no
   !    E0 axis at all.
   !
   !  * Helium content. Dalgarno et al. fix helium at
   !    n(He) = 0.1[n(H) + 2 n(H2)] and vary only n(H2)/n(H) (the paragraph
   !    of their eq. 9). Shull & van Steenberg assume n(He)/n(H) = 0.1 with
   !    equal H and He ionized fractions. Runs at He/H of order 1-10 are
   !    therefore an order of magnitude or more OUTSIDE the composition both
   !    sources sampled, in the direction of more helium. What the split
   !    below corrects for is the BRANCHING between targets at the cell's own
   !    composition; the total ionization budget and the heating fraction
   !    still carry the helium content of the two Monte Carlo calculations,
   !    and no source consulted here covers a helium-dominated gas.
   !
   !  * Metals. Neither source has any metal in its gas. The ionized fraction
   !    the caller hands in is the TOTAL free-electron density over the H and
   !    He nuclei, so in a metals-on run it carries electrons the two Monte
   !    Carlo calculations never had. That is deliberate -- what degrades the
   !    photoelectron is Coulomb scattering off free electrons, and an
   !    electron released by a metal does that work as well as one released
   !    by hydrogen -- but it is an EXTENSION beyond both compositions, not
   !    something either paper supports. It also lets the ratio exceed one,
   !    which is why every entry point clamps x into [0,1].
   !
   !  * Temperature. Dalgarno et al. require H2 to be predominantly in its
   !    ground vibrational level, "which means that the temperature of the
   !    cloud must be less than about 2000 K". Our molecular layer runs
   !    866-1320 K, inside it; a hotter molecular base would not be.
   !
   !  * WHAT THE HEATING EFFICIENCY DOES AND DOES NOT INCLUDE. Dalgarno et
   !    al. define it (their section 5.3) as "the fraction of the primary
   !    energy that is converted into heat from elastic scattering and
   !    rotational excitations": VIBRATIONAL excitation is outside it, and
   !    so is the energy the B and C states return to the ground state when
   !    they fluoresce. Both are added here as separate channels, each
   !    multiplied by the fraction of the excitation energy that is
   !    collisionally de-excited rather than radiated
   !    (h2_vibrational_relaxation.f90). Without them the molecular heat
   !    fraction would be low by a factor of about three at 1 keV.
   !
   ! ---------------------------------------------------------------------
   !  DISSOCIATIVE IONIZATION OF H2 (see docs/Update_EXHALE_stage1.pdf section 130)
   ! ---------------------------------------------------------------------
   !  H2 + e -> H + H+ + 2e is carried, through dissoc_ion_per_H2p below:
   !  the secondary electrons make one such proton for every 22 H2+ ions.
   !  The photon side of the same channel is the H2 photoabsorption cross
   !  section (sigma_H2 in cross_sec.f90: Backx et al. 1976 below 18 eV,
   !  Samson & Haddad 1994 from 18 to 300 eV, the Yan, Sadeghpour &
   !  Dalgarno 1998 sum-rule tail above) times the H+/(H+ + H2+) branching
   !  of Chung et al. (1993) Table II (frac_H2_dissociative_ionization in
   !  cross_sec.f90), and its photoelectron enters this module as its own
   !  absorber, iabs_H2_di.
   !
   !  All quantities cgs except energies, which are eV.

   use h2_vibrational_relaxation, only: e_vib_v1_eV, e_vib_v2_eV

   implicit none
   private

   public :: photoelectron_partition_t
   public :: photoelectron_energy_partition, photoelectron_shares
   public :: photoelectron_energy_grid
   public :: dalgarno_energy_bracket, dalgarno_energy_node_weights
   public :: svs85_fion_HI, svs85_fion_HeI
   public :: dalgarno_x_h2he, h2_ionization_share
   public :: heh_neutral_svs85, dal_h2_weight, n_dal_E, dal_E
   public :: sec_kE, sec_wE, n_abs_fixed
   public :: iabs_HI, iabs_HeI, iabs_HeII, iabs_HeTR, iabs_H2, iabs_H2_di
   public :: dissoc_ion_per_H2p

   ! ----- Composition the SvS85 secondary-ionization fits carry ----- !
   ! Shull & van Steenberg (1985), section II: "We assume that n(He)/n(H) =
   ! 0.1, and that the hydrogen and helium ionization fractions are equal:
   ! n(H+)/n(Htot) = n(He+)/n(Hetot)." Equal ionization fractions make the
   ! NEUTRAL ratio their Monte Carlo sampled equal to the elemental one at
   ! every x, so the branching between the H I and He I secondary channels
   ! that their coefficients carry is the one belonging to
   ! n(He I)/n(H I) = 0.1. That is the reference the split renormalizes from.
   real*8, parameter :: heh_neutral_svs85 = 0.1d0

   ! -- Collision weight of H2 against H, from Dalgarno's own coefficients --
   ! Dalgarno et al. eq. (9) raises the mean energy per H+ by
   ! [1 + 1.89 n(H2)/n(H)] and eq. (10) raises the mean energy per H2+ by
   ! [1 + 0.53 n(H)/n(H2)]. The product 1.89 x 0.53 = 1.0017, so the two
   ! statements are, to 0.17%, one statement: the hydrogenic ionization
   ! energy divides between the two targets in the ratio n(H) : 1.89 n(H2).
   !
   ! THAT IDENTIFICATION IS NOT MADE IN THE PAPER. It is read off here from
   ! the two published coefficients, and it is what lets one weight replace
   ! the two mean-energy expressions -- with the advantage that the weight is
   ! bounded in [0,1] and needs no division by a vanishing density. The 0.17%
   ! is the residual of treating the two coefficients as exactly reciprocal.
   real*8, parameter :: dal_h2_weight = 1.89d0

   ! ---- Dissociative ionization of H2 by the secondary electrons ---- !
   ! Dalgarno, Yan & Liu (1999), the paragraph after their eq. (10): "The
   ! production of H+ ions from dissociative ionization of H2 is given by a
   ! mean energy per ion of 22W_m(H2+). The harmonic mean of 22W_m(H2+) and
   ! W_m(H+) is the mean energy for the production of H+ ions in the gas
   ! mixture." A mean energy 22 times larger is a yield 22 times smaller, and
   ! the harmonic mean says the two H+ sources ADD, so the channel makes one
   ! proton (and one H atom) for every 22 H2+ ions the same electrons make.
   !
   ! WHAT IS NOT RE-BILLED. Dalgarno et al.'s W is the primary energy divided
   ! by the number of ions, so it already contains every loss of the
   ! degradation; the extra 4.5% of ionizations is not deducted from the
   ! Shull & van Steenberg amplitude f_ion this module borrows for the
   ! ionization budget, which would double-count it. The two accountings are
   ! consistent to about the 6% at which the sources' mean energies per H2+
   ! ion already disagree (39.4 eV here against their 41.9 eV).
   real*8, parameter :: dissoc_ion_per_H2p = 1.0d0/22.0d0

   ! ----- Energy grid of the Dalgarno tables ----- !
   integer, parameter :: n_dal_E = 6
   real*8, parameter :: dal_E(n_dal_E) =                                 &
        (/ 30.0d0, 50.0d0, 100.0d0, 200.0d0, 500.0d0, 1000.0d0 /)
   real*8, parameter :: dal_logE(n_dal_E) =                              &
        (/ 1.4771212547196624d0, 1.6989700043360187d0,                   &
           2.0000000000000000d0, 2.3010299956639813d0,                   &
           2.6989700043360187d0, 3.0000000000000000d0 /)

   ! -- Table 4: mean energy per ion pair, eq. (13) W = W0(1 + C x^a) [eV] --
   ! H and He gas mixture rows, (W0, a, C). These carry the ENERGY DEPENDENCE
   ! of the ionization channels, which Shull & van Steenberg's fits do not
   ! have: their Table 2 footnote restricts them to E0 > 100 eV and their
   ! Table 1 is a 3 keV calculation. Used below as the ratio
   ! W(x, 1 keV)/W(x, E0), which is 1 at the top of the table and so leaves
   ! the high-energy limit exactly at Shull & van Steenberg's amplitude.
   ! Transcribed from the published table and checked against the rendered
   ! page; W0(H+) at 1 keV, 39.8 eV, is the value Dalgarno et al. quote in
   ! their section 6 (as 39.5) against Shull & van Steenberg's 37.2.
   real*8, parameter :: t4_Hp_hhe(3,n_dal_E) = reshape( (/              &
        60.0d0, 1.06d0 , 221.0d0,                                       &
        48.9d0, 1.01d0 ,  87.9d0,                                       &
        42.6d0, 0.964d0,  41.0d0,                                       &
        40.3d0, 0.928d0,  24.5d0,                                       &
        39.8d0, 0.890d0,  15.5d0,                                       &
        39.8d0, 0.866d0,  12.2d0 /), (/3,n_dal_E/) )
   real*8, parameter :: t4_Hep_hhe(3,n_dal_E) = reshape( (/             &
        6740.0d0, 1.08d0 , 522.0d0,                                     &
        2150.0d0, 1.07d0 ,  98.4d0,                                     &
        1030.0d0, 1.03d0 ,  40.9d0,                                     &
         700.0d0, 1.02d0 ,  24.1d0,                                     &
         540.0d0, 1.00d0 ,  15.6d0,                                     &
         487.0d0, 0.994d0,  12.5d0 /), (/3,n_dal_E/) )

   ! - Table 7: heating efficiency, eq. (14) eta = 1 + (eta0-1)/(1 + C x^a) -
   ! Rows are (eta0, a, C) at the six tabulated energies; two-gas mixtures.
   ! Transcribed from the published table and checked against the rendered
   ! page. The H, He block reproduces their Figure 23d to 5% at x = 0.1.
   real*8, parameter :: t7_hhe(3,n_dal_E) = reshape( (/                  &
        0.259d0, 0.984d0, 181.0d0,                                       &
        0.186d0, 0.834d0,  53.4d0,                                       &
        0.148d0, 0.781d0,  25.8d0,                                       &
        0.132d0, 0.743d0,  15.8d0,                                       &
        0.122d0, 0.703d0,  10.1d0,                                       &
        0.117d0, 0.678d0,   7.95d0 /), (/3,n_dal_E/) )
   real*8, parameter :: t7_h2he(3,n_dal_E) = reshape( (/                 &
        0.100d0, 0.515d0,  17.7d0,                                       &
        0.083d0, 0.468d0,   8.89d0,                                      &
        0.069d0, 0.438d0,   5.17d0,                                      &
        0.062d0, 0.410d0,   3.54d0,                                      &
        0.057d0, 0.381d0,   2.54d0,                                      &
        0.055d0, 0.366d0,   2.17d0 /), (/3,n_dal_E/) )

   ! -- Table 5: mean energy per excitation, eq. (13) W = W0(1 + C x^a) [eV]
   ! H2 and He gas mixture rows, (W0, a, C).
   ! B1Sigma_u+ and C1Pi_u include cascading from higher singlet states, as
   ! the published table states. "Dissociation" is the TOTAL dissociation
   ! yield excluding dissociative excitation [H + H(2p)], per the table's own
   ! footnote b, so it contains the singlet predissociation as well as the
   ! triplet route and the two must be separated before a heat is assigned.
   real*8, parameter :: t5_B(3,n_dal_E) = reshape( (/                    &
        124.0d0, 1.12d0,  221.0d0,                                       &
        113.0d0, 0.978d0,  55.3d0,                                       &
        108.0d0, 0.925d0,  25.9d0,                                       &
        111.0d0, 0.872d0,  15.1d0,                                       &
        115.0d0, 0.813d0,   9.24d0,                                      &
        117.0d0, 0.779d0,   7.09d0 /), (/3,n_dal_E/) )
   real*8, parameter :: t5_C(3,n_dal_E) = reshape( (/                    &
        168.0d0, 1.12d0,  228.0d0,                                       &
        145.0d0, 0.988d0,  55.7d0,                                       &
        135.0d0, 0.936d0,  26.0d0,                                       &
        134.0d0, 0.887d0,  15.1d0,                                       &
        134.0d0, 0.833d0,   9.05d0,                                      &
        132.0d0, 0.802d0,   6.88d0 /), (/3,n_dal_E/) )
   real*8, parameter :: t5_diss(3,n_dal_E) = reshape( (/                 &
         34.7d0, 0.872d0,  126.0d0,                                      &
         49.6d0, 0.799d0,   69.7d0,                                      &
         59.2d0, 0.713d0,   38.4d0,                                      &
         80.5d0, 0.661d0,   34.7d0,                                      &
         89.5d0, 0.603d0,   25.9d0,                                      &
         92.6d0, 0.574d0,   22.0d0 /), (/3,n_dal_E/) )

   ! Table 5, H2 and He gas mixture, the DIRECT electron-impact excitation of
   ! v = 1 and v = 2 out of v = 0 (rows "v = 1" and "v = 2"). These are not
   ! the fluorescence cascade, which comes from the B and C rows above, so
   ! the two channels do not overlap. The fits for H2 vibrational states are
   ! the loosest in the paper: "The fits for excitation of He triplet states
   ! and H2 vibrational states are good to 15%".
   real*8, parameter :: t5_v1(3,n_dal_E) = reshape( (/                   &
          4.15d0, 0.845d0,  9650.0d0,                                    &
          5.12d0, 0.866d0, 11400.0d0,                                    &
          6.31d0, 0.890d0, 12400.0d0,                                    &
          6.96d0, 0.923d0, 18500.0d0,                                    &
          7.61d0, 0.950d0, 22400.0d0,                                    &
          7.81d0, 0.955d0, 23500.0d0 /), (/3,n_dal_E/) )
   real*8, parameter :: t5_v2(3,n_dal_E) = reshape( (/                   &
         58.0d0, 0.768d0,  3300.0d0,                                     &
         69.7d0, 0.797d0,  4300.0d0,                                     &
         85.5d0, 0.814d0,  5120.0d0,                                     &
         96.5d0, 0.865d0,  7710.0d0,                                     &
        108.0d0, 0.902d0,  9690.0d0,                                     &
        109.0d0, 0.907d0, 10700.0d0 /), (/3,n_dal_E/) )

   ! ----- Dissociation energetics ----- !
   ! Radiative dissociation probability of the two singlet states, Dalgarno
   ! et al. section 5.3: "with a fractional probability of 0.27 for the
   ! B1Sigma_u+ state and of 0.036 for the C1Pi_u state (Stephens & Dalgarno
   ! 1972; Abgrall et al. 1992, 1997)."
   real*8, parameter :: p_diss_B = 0.27d0
   real*8, parameter :: p_diss_C = 0.036d0
   ! Kinetic energy released to the fragments when H2 is dissociated through
   ! a TRIPLET state. Dalgarno et al. section 5.3: the energy distribution
   ! "cannot be much different from a mean of 5.4 eV corresponding to a
   ! vertical transition".
   real*8, parameter :: e_diss_triplet = 5.4d0
   ! Kinetic energy released when a singlet state predissociates into the
   ! vibrational continuum of the ground state: Black & Dalgarno (1977, ApJS
   ! 34, 405), p. 418, "the yield is about 0.4 eV per atom pair". The same
   ! constant drives the photon-side Lyman-Werner channel
   ! (lyman_werner.f90, e_lw_fragment_erg).
   real*8, parameter :: e_diss_singlet = 0.4d0

   ! Smallest denominator treated as nonzero.
   real*8, parameter :: tiny_den = 1.0d-99

   ! -------------------------------------------------------------------- !
   ! The cell's partition, held at the six tabulated primary energies. The
   ! expensive part -- the powers x^a -- is evaluated once per cell, and
   ! every absorber then costs one gather and one linear interpolation per
   ! spectral bin, on weights that were built once at startup.
   type :: photoelectron_partition_t
      real*8  :: c_ion_HI (n_dal_E) = 0.0d0   ! per H I atom    [cm^3]
      real*8  :: c_ion_HeI(n_dal_E) = 0.0d0   ! per He I atom   [cm^3]
      real*8  :: c_ion_H2 (n_dal_E) = 0.0d0   ! per H2 molecule [cm^3]
      real*8  :: f_heat   (n_dal_E) = 1.0d0   ! heat share of the primary
      logical :: molecular = .false.
   end type photoelectron_partition_t

   ! -------------------------------------------------------------------- !
   ! Where each absorber's photoelectron energy E0 = hv - E_th falls in the
   ! Dalgarno energy grid. E0 depends only on the spectral bin and on the
   ! absorber's own threshold, both fixed for a run, so the bracket index and
   ! the interpolation weight are built once by photoelectron_energy_grid and
   ! read afterwards. Absorbers are the six H/He/H2 edges below, then the
   ! metal ions in their species_table order at n_abs_fixed + i.
   integer, parameter :: iabs_HI   = 1
   integer, parameter :: iabs_HeI  = 2
   integer, parameter :: iabs_HeII = 3
   integer, parameter :: iabs_HeTR = 4
   integer, parameter :: iabs_H2   = 5
   ! The dissociative H2 photoionization channel is its own absorber: its
   ! photoelectron carries hv - 18.08 eV, not hv - e_th_H2, so it lands
   ! elsewhere in the Dalgarno energy grid and thermalizes a different share.
   integer, parameter :: iabs_H2_di = 6
   integer, parameter :: n_abs_fixed = 6
   integer, allocatable, save :: sec_kE(:,:)   ! (Nl, n_abs)
   real*8,  allocatable, save :: sec_wE(:,:)   ! (Nl, n_abs)

   contains

   ! ================================================================= !
   !  Startup: the energy-interpolation table
   ! ================================================================= !

   subroutine photoelectron_energy_grid(e_v, e_th)
   ! Bracket index and weight of every (spectral bin, absorber) pair in the
   ! Dalgarno energy grid, for the photoelectron energy E0 = e_v - e_th.
   !
   ! Linear in log10(E0) and CLAMPED at both ends. Clamping below 30 eV is
   ! inert where it matters, because the caller applies the partition only
   ! above E_sec_ion = 30 eV, which is the bottom of the table; clamping
   ! above 1 keV freezes the coefficients at their most nearly asymptotic
   ! values, which is the behaviour the source describes at high energy.

   real*8, intent(in) :: e_v(:)     ! spectral grid [eV]
   real*8, intent(in) :: e_th(:)    ! ionization threshold of each absorber
   integer :: nl_loc, na, i, a, k
   real*8  :: w

   nl_loc = size(e_v)
   na     = size(e_th)
   if (allocated(sec_kE)) deallocate(sec_kE)
   if (allocated(sec_wE)) deallocate(sec_wE)
   allocate(sec_kE(nl_loc,na), sec_wE(nl_loc,na))

   do a = 1, na
      do i = 1, nl_loc
         call dalgarno_energy_bracket(e_v(i) - e_th(a), k, w)
         sec_kE(i,a) = k
         sec_wE(i,a) = w
      enddo
   enddo

   end subroutine photoelectron_energy_grid

   pure subroutine dalgarno_energy_bracket(E0, k, w)
   ! Bracket index k and weight w of an electron energy E0 [eV] in the
   ! Dalgarno energy grid: a quantity tabulated at dal_E is read at E0 as
   ! (1 - w) q(k) + w q(k+1), linear in log10(E0) and clamped at both ends
   ! (the validity note of photoelectron_energy_grid). The one definition
   ! of that interpolation, for the photoelectrons of the spectral grid and
   ! for an electron of any other energy (the Auger electrons of a metal
   ! inner-shell absorption, dalgarno_energy_node_weights).
   real*8,  intent(in)  :: E0
   integer, intent(out) :: k
   real*8,  intent(out) :: w
   real*8 :: le
   le = log10(max(E0, tiny_den))
   if (le .le. dal_logE(1)) then
      k = 1;  w = 0.0d0
   else if (le .ge. dal_logE(n_dal_E)) then
      k = n_dal_E - 1;  w = 1.0d0
   else
      k = 1
      do while (k .lt. n_dal_E-1 .and. le .gt. dal_logE(k+1))
         k = k + 1
      enddo
      w = (le - dal_logE(k))/(dal_logE(k+1) - dal_logE(k))
   endif
   end subroutine dalgarno_energy_bracket

   pure subroutine dalgarno_energy_node_weights(E0, wt)
   ! The same interpolation written as weights on the tabulated energies:
   ! a partition coefficient q tabulated at dal_E is sum_k wt(k) q(k) at
   ! E0. A photoabsorption that releases several electrons of different
   ! energies (a photoelectron and Auger electrons) is the sum of their
   ! weights times their energies, so the partition of all of them is one
   ! contraction with the cell's coefficients (photoelectron_partition_t)
   ! whatever the number of electrons.
   real*8, intent(in)  :: E0
   real*8, intent(out) :: wt(n_dal_E)
   integer :: k
   real*8  :: w
   call dalgarno_energy_bracket(E0, k, w)
   wt      = 0.0d0
   wt(k)   = 1.0d0 - w
   wt(k+1) = w
   end subroutine dalgarno_energy_node_weights

   ! ================================================================= !
   !  Shull & van Steenberg (1985) fits: their equations (1) and (2)
   !  with the coefficients of their Table 2.
   ! ================================================================= !

   ! A possible future refinement is to couple the SvS85 Ly-alpha excitation
   ! channel f_exc,Lya = 0.4766*(1-x^0.2735)^1.5221 to the Ly-alpha field of
   ! the excited-H model; here that energy is assumed to escape as line
   ! radiation and is not put into any rate.
   elemental function svs85_fion_HI(x) result(f)
      real*8, intent(in) :: x
      real*8 :: f
      f = 0.3908d0*(1.0d0 - x**0.4092d0)**1.7592d0
   end function svs85_fion_HI

   elemental function svs85_fion_HeI(x) result(f)
      real*8, intent(in) :: x
      real*8 :: f
      f = 0.0554d0*(1.0d0 - x**0.4614d0)**1.6660d0
   end function svs85_fion_HeI

   ! ================================================================= !
   !  Dalgarno, Yan & Liu (1999)
   ! ================================================================= !

   ! Ionized fraction to be used with the H2-He rows, their eq. (15). The
   ! H-H2-He mixture counts electrons per H and He NUCLEUS, the H2-He mixture
   ! per H2 MOLECULE and He atom, and this is the conversion between them.
   elemental function dalgarno_x_h2he(x) result(xp)
      real*8, intent(in) :: x
      real*8 :: xp
      xp = 1.83d0*x/(1.0d0 + 0.83d0*x)
   end function dalgarno_x_h2he

   ! Share of the hydrogenic ionization energy taken by H2 rather than H,
   ! from Dalgarno et al. eqs. (9) and (10) written as the collision weight
   ! the two coefficients imply (dal_h2_weight above). Bounded in [0,1] and
   ! finite as either density vanishes.
   elemental function h2_ionization_share(n_HI, n_H2) result(u)
      real*8, intent(in) :: n_HI, n_H2
      real*8 :: u, d
      d = n_HI + dal_h2_weight*n_H2
      if (d .gt. 0.0d0) then
         u = dal_h2_weight*n_H2/d
      else
         u = 0.0d0
      endif
   end function h2_ionization_share

   ! Share of the H2-He mixture's VIBRATIONAL excitation yield that survives
   ! in a mixture that also contains atomic hydrogen: Dalgarno et al.
   ! eqs. (11), (16) and (17), which are one expression,
   !
   !     Y = Y_HeH2 * r/(r + a) ,   r = n(H2)/n(H) ,
   !
   ! with a = 0.5 for x >= 1e-4 (their eq. 16, since 2n(H2)/[n(H)+2n(H2)] is
   ! r/(r+0.5)), a = a(x) for 1e-7 <= x < 1e-4 (their eq. 17), and a = 0.177
   ! for smaller x (their eq. 11). It is separate from the ionization weight
   ! above because vibrational excitation is done by the slowest secondaries,
   ! which the free electrons compete for; the ionization is not.
   !
   ! A TYPOGRAPHICAL CORRECTION, MADE HERE AND STATED. The paper prints
   ! "a(x) = 0.5 x (x/10^4)^0.15". With the exponent as printed, eq. (17)
   ! gives a = 0.0316 at x = 1e-4 where eq. (16) requires 0.5, and a = 0.0112
   ! at x = 1e-7 where eq. (11) requires 0.177 -- two discontinuities, both
   ! by the SAME factor 15.85. With 10^-4 in place of 10^4 the expression
   ! joins eq. (16) exactly at x = 1e-4 (a = 0.5) and eq. (11) exactly at
   ! x = 1e-7 (a = 0.5 x 10^-0.45 = 0.1774 against their 0.177), which is
   ! what the three formulas are stated to do. 15.85 = (10^8)^0.15 is exactly
   ! the ratio between the two readings, so the sign of the exponent is the
   ! whole discrepancy. The minus sign is taken to have been lost in
   ! typesetting; the same page prints a multiplication sign that the text
   ! extraction of this paper renders as a plus.
   elemental function h2_vibrational_share(x, n_HI, n_H2) result(wv)
      real*8, intent(in) :: x, n_HI, n_H2
      real*8 :: wv, a, d
      if (x .ge. 1.0d-4) then
         a = 0.5d0
      else if (x .ge. 1.0d-7) then
         a = 0.5d0*(x*1.0d4)**0.15d0
      else
         a = 0.177d0
      endif
      d = n_H2 + a*n_HI
      if (d .gt. 0.0d0) then
         wv = n_H2/d
      else
         wv = 0.0d0
      endif
   end function h2_vibrational_share

   ! Heating efficiency at one tabulated energy, their eq. (14):
   ! the fraction of the primary energy that ends as heat, through elastic
   ! scattering, Coulomb collisions with the thermal electrons, and -- in a
   ! molecular gas -- rotational excitation of H2 ("the heating efficiency,
   ! defined as the fraction of the primary energy that is converted into
   ! heat from elastic scattering and rotational excitations", their section
   ! 5.3). Vibrational excitation is NOT in it, and is a separate channel.
   !
   ! VALIDITY. Dalgarno et al. state eq. (14) "for x <= 0.1 with the values
   ! of eta0, c, and a listed in Table 7", and fit it to within 15%. IT IS
   ! EVALUATED HERE AT EVERY x, INCLUDING x -> 1 IN THE IONIZED WIND, WHICH
   ! IS EXTRAPOLATION BEYOND THAT STATEMENT. The closure applied in
   ! photoelectron_energy_partition forces the correct value at x = 1, but
   ! forcing an endpoint is not the same thing as the fit being valid in
   ! between, and the extrapolation is accepted on a measurement rather than
   ! on the closure: judged against Shull & van Steenberg's own Monte Carlo
   ! (their Table 1, an 18-point 3 keV calculation with 1 sigma widths) over
   ! x > 0.1, eq. (14) closed as below deviates by at most 5.6% and 1.01
   ! sigma, against 3.9% and 1.05 sigma for Shull & van Steenberg's own fit
   ! -- while below x = 0.1, where our molecular layer lives, eq. (14) stays
   ! inside 1.40 sigma and their fit reaches 2.77 sigma at x = 1e-4 and
   ! collapses to zero below it. The full point-by-point table is in
   ! docs/Update_EXHALE_stage1.pdf section 118.
   pure function dalgarno_heating_efficiency(tab, k, x) result(eta)
      real*8, intent(in)  :: tab(3,n_dal_E)
      integer, intent(in) :: k
      real*8, intent(in)  :: x
      real*8 :: eta
      eta = 1.0d0 + (tab(1,k) - 1.0d0)/(1.0d0 + tab(3,k)*x**tab(2,k))
   end function dalgarno_heating_efficiency

   ! Mean energy per ion pair or per excitation at one tabulated energy,
   ! their eq. (13) [eV].
   pure function dalgarno_mean_energy_per_ion_pair(tab, k, x) result(W_eV)
      real*8, intent(in)  :: tab(3,n_dal_E)
      integer, intent(in) :: k
      real*8, intent(in)  :: x
      real*8 :: W_eV
      W_eV = tab(1,k)*(1.0d0 + tab(3,k)*x**tab(2,k))
   end function dalgarno_mean_energy_per_ion_pair

   ! ================================================================= !
   !  The cell partition
   ! ================================================================= !

   pure subroutine photoelectron_energy_partition(x, n_HI, n_HeI, n_H2,   &
                                                  f_vib_heat, e_vib_bound, &
                                                  p)

   ! Branching of a fast photoelectron's excess energy between heat and the
   ! ionization of each neutral target, on the cell's own composition.
   !
   ! The three ionization channels are returned PER TARGET PARTICLE [cm^3],
   ! so the caller multiplies them by the deposited-energy integrand and adds
   ! the result straight to the photoionization rate, with no division by a
   ! neutral density anywhere.
   !
   ! WHAT IS TAKEN FROM WHERE.
   !
   ! The energy spent on ionization is Shull & van Steenberg's amplitude with
   ! Dalgarno et al.'s energy dependence,
   !
   !    f_HI (x,E0) = phi_HI (x)  * W_H+ (x, 1 keV) / W_H+ (x, E0)
   !    f_HeI(x,E0) = phi_HeI(x)  * W_He+(x, 1 keV) / W_He+(x, E0)
   !
   ! with W from their Table 4 and eq. (13). THE RATIO IS 1 AT 1 keV, so the
   ! high-energy limit is exactly Shull & van Steenberg's -- which is what
   ! their fits are, an E0 > 100 eV asymptote fitted to a 3 keV calculation.
   ! Below that a photoelectron makes fewer ion pairs per unit energy, and
   ! for helium the effect is large: W_He+ rises from 487 eV at 1 keV to
   ! 6740 eV at 30 eV, so a 30 eV photoelectron puts 14 times less of its
   ! energy into He I ionization than the unresolved fit asserts. That is
   ! physical -- 30 eV is barely above helium's 24.59 eV threshold.
   !
   ! Only the DIVISION of that energy among the targets is recomputed, on the
   ! cell's own densities, at each tabulated energy:
   !
   !    w_tot     = (n_HI + 1.89 n_H2) * href * f_HI  +  n_HeI * f_HeI
   !    c_ion_HI  = f_ion *        href*f_HI / w_tot
   !    c_ion_H2  = f_ion * 1.89 * href*f_HI / w_tot
   !    c_ion_HeI = f_ion *             f_HeI / w_tot
   !
   ! with href = heh_neutral_svs85 and f_ion = f_HI + f_HeI. No
   ! electron-impact cross sections are needed: their ratio cancels against
   ! the reference composition the fits were measured at. The construction
   ! conserves energy identically --
   !
   !    n_HI*c_ion_HI + n_H2*c_ion_H2 + n_HeI*c_ion_HeI = f_ion
   !
   ! -- at every energy and for every composition, which is the property that
   ! makes it safe to move energy between the channels.
   !
   ! The collision weight 1.89 is NOT given an energy dependence: Dalgarno et
   ! al. state their eqs. (9) and (10) coefficients once, not per energy.
   !
   ! LIMITS, which are the requirements this form was written to meet.
   !  * n_H2 = 0 and n_HeI/n_HI = 0.1, the gas Shull & van Steenberg
   !    sampled: c_ion_HI = f_HI/n_HI and c_ion_HeI = f_HeI/n_HeI exactly.
   !  * n_H2 = 0 at any He/H: the expressions reduce to the two-target form
   !    the code used before H2 existed, TERM BY TERM.
   !  * n_HI -> 0 with H2 present: the H I rate per atom stays finite, the
   !    volumetric one vanishes, and the hydrogenic budget goes to H2.
   !  * n_H2 -> 0 with H present: the mirror image.
   !
   ! THE H2 IONIZATION CHANNEL, CHECKED AGAINST ITS SOURCE. In a neutral,
   ! fully molecular gas this construction gives a mean energy per H2+ ion of
   ! e_th_H2/f_HI(0) = 15.4259/0.3908 = 39.5 eV, against the 41.9 eV of Dalgarno
   ! et al. Table 4 (H2 and He mixture, 1 keV): 6%, the same size as the
   ! disagreement between the two sources' mean energy per H+ (Dalgarno et
   ! al. section 6 quote 37.2 eV for Shull & van Steenberg against their own
   ! 39.5 eV) and well inside Dalgarno et al.'s stated 15%.
   !
   ! THE HEATING FRACTION is Dalgarno et al.'s, at the cell's own H2 content
   ! and at each of their six tabulated primary energies,
   !
   !    f_heat(E0) = 1 - [1 - eta_mix(x,E0,r)]*(1 - x)
   !                 + f_diss + f_vib + f_fluor
   !
   ! with eta_mix their eq. (12), the H2/H mixing of the two two-gas heating
   ! efficiencies. The three molecular terms are the kinetic energy of the
   ! fragments of a dissociated H2, the direct vibrational excitation of
   ! v = 1 and 2, and the vibrational energy the B and C states hand back to
   ! the ground state when they fluoresce instead of dissociating. All three
   ! vanish with n(H2); the last two also carry the collisional-quench
   ! fraction f_vib_heat, since a radiated excitation is not heat.
   !
   ! THE (1 - x) FACTOR IS A CLOSURE AND IS OURS, not either source's. It
   ! states that the channels that compete with heat -- ionization and
   ! excitation -- need a NEUTRAL to act on, so their share must vanish in
   ! proportion to the neutral fraction of the H and He nuclei, which is
   ! 1 - x when x is the electron density over those nuclei. Dalgarno et
   ! al.'s eq. (14) has no such factor: at x = 1 it still leaves
   ! (1 - eta0)/(1 + c) of the primary energy unthermalized, 10% of it at
   ! 1 keV, which is unphysical in a fully ionized gas and is why Shull &
   ! van Steenberg chose fitting forms "because of their appropriate limits
   ! at x = 0 and x = 1". The factor is 1 to within 1e-4 over the whole
   ! molecular layer, so it does not touch the regime the H2 terms are for;
   ! what it fixes is the ionized wind above.
   !
   ! WHY NOT JOIN THE TWO SOURCES AT x = 0.1 INSTEAD. Because there is no
   ! seam to define. Shull & van Steenberg's Table 2 fits their Figure 4,
   ! captioned "Limiting values (E0 > 100 eV)", and their Table 1 is a 3 keV
   ! calculation: their prescription is a single high-energy asymptote with
   ! no E0 axis, while Dalgarno et al.'s resolves E0 and moves by 25% between
   ! 100 eV and 1 keV at x = 0.1. Matching a curve to a point that exists at
   ! only one energy would put a discontinuity in E0 wherever the join was
   ! placed. One expression over the whole domain has none.
   !
   ! f_diss is the kinetic energy of the fragments of an H2 molecule
   ! dissociated by the photoelectron, as a share of the primary energy. It
   ! carries NO collisional-quench factor and should not: the energy is
   ! translational from the moment it is released. The vibrational and
   ! fluorescence terms are internal excitation and do carry one.

   real*8, intent(in)  :: x, n_HI, n_HeI, n_H2
   ! Fraction of an H2 vibrational excitation that is collisionally
   ! de-excited and so becomes heat, rather than being radiated away in the
   ! infrared quadrupole lines (h2_vibrational_relaxation.f90). Pass zero to
   ! keep both vibrational channels out; they are zero in any case where
   ! n_H2 = 0.
   real*8, intent(in)  :: f_vib_heat
   ! Mean internal energy [eV] of the X(v'', J'') level a B or C
   ! fluorescence lands on, h2_energy_per_bound_fluorescence_eV at the cell's
   ! temperature. It reaches the answer only through f_fluor below, which
   ! f_vib_heat already gates to zero for a cell with no molecular hydrogen.
   real*8, intent(in)  :: e_vib_bound
   type(photoelectron_partition_t), intent(out) :: p

   real*8  :: xc, xp, f_HI, f_HeI, f_ion, w_tot, n_Hw, u, n_Hden, wv, fq
   real*8  :: phi_HI, phi_HeI, rE_HI, rE_HeI, W_top_HI, W_top_HeI
   real*8  :: eta_hhe, eta_h2he, eta_mix, f_diss, f_vib, f_fluor, fh
   real*8  :: Wd, WB, WC, Nd, NB, NC, Ns, Nt
   integer :: k

   xc      = min(max(x, 0.0d0), 1.0d0)
   phi_HI  = svs85_fion_HI(xc)
   phi_HeI = svs85_fion_HeI(xc)

   p%molecular = (n_H2 .gt. 0.0d0)

   xp      = dalgarno_x_h2he(xc)
   u       = h2_ionization_share(n_HI, n_H2)
   wv      = h2_vibrational_share(xc, n_HI, n_H2)
   fq      = min(max(f_vib_heat, 0.0d0), 1.0d0)
   n_Hw    = n_HI + dal_h2_weight*n_H2
   n_Hden  = 10.0d0*n_H2 + n_HI
   f_diss  = 0.0d0
   f_vib   = 0.0d0
   f_fluor = 0.0d0

   ! Mean energies at the top of the Dalgarno grid, the reference the energy
   ! ratios are taken against.
   W_top_HI  = dalgarno_mean_energy_per_ion_pair(t4_Hp_hhe , n_dal_E, xc)
   W_top_HeI = dalgarno_mean_energy_per_ion_pair(t4_Hep_hhe, n_dal_E, xc)

   do k = 1, n_dal_E

      ! ---- ionization channels at this primary energy ---- !
      rE_HI  = W_top_HI                                                  &
               /max(dalgarno_mean_energy_per_ion_pair(t4_Hp_hhe , k, xc),&
                    tiny_den)
      rE_HeI = W_top_HeI                                                 &
               /max(dalgarno_mean_energy_per_ion_pair(t4_Hep_hhe, k, xc),&
                    tiny_den)
      f_HI  = phi_HI *rE_HI
      f_HeI = phi_HeI*rE_HeI
      f_ion = f_HI + f_HeI

      w_tot = n_Hw*heh_neutral_svs85*f_HI + n_HeI*f_HeI
      if (w_tot .gt. 0.0d0) then
         p%c_ion_HI (k) = f_ion*heh_neutral_svs85*f_HI/w_tot
         p%c_ion_H2 (k) = dal_h2_weight*p%c_ion_HI(k)
         p%c_ion_HeI(k) = f_ion*f_HeI/w_tot
      else
         ! Either no neutral target is left, or x = 1 has sent both fits to
         ! zero: there is no secondary ionization to distribute.
         p%c_ion_HI (k) = 0.0d0
         p%c_ion_H2 (k) = 0.0d0
         p%c_ion_HeI(k) = 0.0d0
      endif

      ! ---- heating fraction at this primary energy ---- !
      eta_hhe = dalgarno_heating_efficiency(t7_hhe , k, xc)

      if (p%molecular) then
         eta_h2he = dalgarno_heating_efficiency(t7_h2he, k, xp)
         ! Dalgarno et al. eq. (12), multiplied through by n(H) so that it
         ! stays finite when the gas is fully molecular.
         eta_mix = (10.0d0*n_H2*eta_h2he + n_HI*eta_hhe)/n_Hden

         ! Dissociation heat. The published "Dissociation" yield is the
         ! total, so the singlet route -- which releases 0.4 eV, not 5.4 --
         ! is subtracted out before the triplet energy is applied.
         ! Reproduces the per-ion-pair dissociation heat Dalgarno et al.
         ! state in their sections 5.3 and 6 to 2% at 30 eV and 12% at
         ! 1 keV, against 10% and 44% for the same expression without the
         ! split.
         Wd = dalgarno_mean_energy_per_ion_pair(t5_diss, k, xp)
         WB = dalgarno_mean_energy_per_ion_pair(t5_B   , k, xp)
         WC = dalgarno_mean_energy_per_ion_pair(t5_C   , k, xp)
         Nd = dal_E(k)/max(Wd, tiny_den)
         NB = dal_E(k)/max(WB, tiny_den)
         NC = dal_E(k)/max(WC, tiny_den)
         Ns = p_diss_B*NB + p_diss_C*NC
         Nt = max(Nd - Ns, 0.0d0)
         f_diss = u*(e_diss_triplet*Nt + e_diss_singlet*Ns)/dal_E(k)

         ! DIRECT vibrational excitation of v = 1 and v = 2 out of v = 0,
         ! Table 5's own rows for those levels, mixed onto the cell's H/H2
         ! ratio by the rule of eqs. (11)/(16)/(17) and delivered as heat
         ! only in so far as the molecule is collisionally de-excited.
         f_vib = fq*wv                                                    &
                 *( e_vib_v1_eV                                           &
                    /max(dalgarno_mean_energy_per_ion_pair(t5_v1,k,xp),   &
                         tiny_den)                                        &
                  + e_vib_v2_eV                                           &
                    /max(dalgarno_mean_energy_per_ion_pair(t5_v2,k,xp),   &
                         tiny_den) )

         ! FLUORESCENCE return of the B and C excitations that do NOT
         ! dissociate -- 1 - 0.27 of B and 1 - 0.036 of C -- each leaving
         ! the molecule with h2_energy_per_bound_fluorescence_eV(T) of
         ! internal energy in the bound ground state. That mean landing
         ! energy is computed for PHOTON pumping; the upper states and the
         ! cascade are the same under electron impact but the distribution
         ! over the upper vibrational levels is not, so sharing it here is an
         ! assumption (h2_vibrational_relaxation.f90 header). The
         ! complementary, dissociating branches are the singlet part of
         ! f_diss above, so the two do not overlap, and neither overlaps
         ! f_vib, which is direct impact out of v = 0.
         f_fluor = fq*u*e_vib_bound                                       &
                   *( (1.0d0 - p_diss_B)/max(WB, tiny_den)                &
                    + (1.0d0 - p_diss_C)/max(WC, tiny_den) )
      else
         eta_mix = eta_hhe
      endif

      ! The heat and the ionization channels are shares of one primary
      ! energy and cannot together exceed it; the excitation channels the
      ! code discards as escaping line radiation are what shrink first.
      fh = 1.0d0 - (1.0d0 - eta_mix)*(1.0d0 - xc)                        &
           + f_diss + f_vib + f_fluor
      p%f_heat(k) = min(max(fh, 0.0d0), max(1.0d0 - f_ion, 0.0d0))
   enddo

   end subroutine photoelectron_energy_partition

   ! The partition of one cell, expanded onto the spectral grid for one
   ! absorber. iabs selects the absorber, whose photoelectron energy
   ! E0 = hv - E_th placed the interpolation weights of sec_kE/sec_wE at
   ! startup; every bin then costs one gather and one linear interpolation.
   pure subroutine photoelectron_shares(p, iabs, f_heat, c_HI, c_HeI, c_H2)
   type(photoelectron_partition_t), intent(in) :: p
   integer, intent(in) :: iabs
   real*8, intent(out) :: f_heat(:), c_HI(:), c_HeI(:), c_H2(:)

   real*8  :: w
   integer :: i, k

   do i = 1, size(f_heat)
      k = sec_kE(i,iabs)
      w = sec_wE(i,iabs)
      f_heat(i) = p%f_heat   (k) + w*(p%f_heat   (k+1) - p%f_heat   (k))
      c_HI  (i) = p%c_ion_HI (k) + w*(p%c_ion_HI (k+1) - p%c_ion_HI (k))
      c_HeI (i) = p%c_ion_HeI(k) + w*(p%c_ion_HeI(k+1) - p%c_ion_HeI(k))
      c_H2  (i) = p%c_ion_H2 (k) + w*(p%c_ion_H2 (k+1) - p%c_ion_H2 (k))
   enddo

   end subroutine photoelectron_shares

   ! End of module
   end module electron_energy_degradation

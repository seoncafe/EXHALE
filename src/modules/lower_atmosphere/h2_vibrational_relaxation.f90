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
   !  * COLLISION PARTNERS. Hollenbach & McKee's eq. (6.29) gives rate
   !    coefficients for atomic H and for H2 only, and Burton et al.'s
   !    Table 7 adds electrons; NEITHER SOURCE GIVES A HELIUM RATE. A
   !    helium-dominated gas is therefore outside both, and helium is left
   !    out of the collider sum below rather than given a borrowed rate.
   !    That understates the de-excitation rate, so it errs towards
   !    radiating the energy away rather than towards claiming it as heat.
   !    Electrons are also left out: their rate coefficient is about two
   !    orders above the H2 one (Burton et al. Table 7), but the molecular
   !    layer runs at an electron fraction of 1e-16 to 5e-5, so they
   !    contribute a part in 1e3 or less of the collider sum there.
   !  * TEMPERATURE. The fits of eq. (6.29) are stated "good to 20% for
   !    T >~ 500 K and to 50% for T >~ 300 K". The molecular layer runs
   !    900-1350 K, inside the tighter of the two.
   !  * A SINGLE LEVEL. n_cr is built from v = 1, the level the energy ends
   !    up in after the cascade. Burton et al. use a pseudolevel at v = 6
   !    with A_(6->0) = 2e-7 s^-1 and gamma taken as gamma_(1->0); using
   !    v = 1 with its own A of 8.3e-7 s^-1 gives a critical density four
   !    times higher, i.e. the more conservative of the two.
   !
   !  The distinction only matters where n is comparable to n_cr. At
   !  900-1350 K, n_cr(v=1) is 2e6-6e6 cm^-3 against H2 and 5e4-8e4 cm^-3
   !  against atomic H, while the molecular layer carries
   !  n(H2) = 4e9 to 5e13 cm^-3: the fraction below is 1 to within a part in
   !  1e3 everywhere in it. It is written out rather than assumed because a
   !  shallower base, or the top of a thinner layer, would not be.
   !
   !  Energies in eV where marked, otherwise cgs.

   implicit none
   private

   public :: h2_energy_per_bound_fluorescence_eV,                         &
             h2_energy_per_bound_fluorescence_erg
   public :: e_vib_v1_eV, e_vib_v2_eV
   public :: h2_vibrational_heat_fraction, h2_vibrational_critical_density

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
   ! Source and arithmetic: docs/p39_lw_cross_section_sources.md sec. 5.2.
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

   ! Collisional de-excitation rate coefficients, Hollenbach & McKee (1979)
   ! eq. (6.29) [cm^3 s^-1]:
   !    gamma_10^H  = 1.0e-12 T^1/2 exp[-(1000/T)]
   !    gamma_10^H2 = 1.4e-12 T^1/2 exp{-[12000/(T + 1200)]}
   !
   ! A DISCREPANCY BETWEEN THE TWO SOURCES, NOT RESOLVED HERE. Burton et al.'s
   ! Table 7 reprints the same H2-H2 coefficient with 18100 in the exponent
   ! where Hollenbach & McKee's eq. (6.29) has 12000; their three H-H2 rows
   ! match eq. (6.29) exactly, so only this one row differs. Hollenbach &
   ! McKee is followed here because it is the primary of the two and because
   ! its H rows are the ones Burton et al. reproduce. Which is right has NOT
   ! been checked against the original rate source (Chu 1977; Shull &
   ! Hollenbach 1978). It changes n_cr by a factor of 16 at 1000 K -- 4.4e6
   ! against 7e7 cm^-3 -- and neither reading comes within three decades of
   ! the layer's own density, so the heat fraction below is 1 either way.
   real*8, parameter :: g10_h_pre  = 1.0d-12
   real*8, parameter :: g10_h2_pre = 1.4d-12

   ! Quadratic fit to the nine computed points (see the block above).
   real*8, parameter :: eb_c0 =  2.0247811011d0
   real*8, parameter :: eb_c1 =  5.6035354525d-05
   real*8, parameter :: eb_c2 = -7.4141671520d-09
   real*8, parameter :: eb_tmin = 700.0d0, eb_tmax = 3200.0d0
   real*8, parameter :: ev2erg = 1.602176634d-12

   contains

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

   ! Rate coefficient for collisional de-excitation of H2 v = 1 by atomic
   ! hydrogen and by H2, at temperature T [K].
   elemental function gamma_10_H(T) result(gamma_deex)
      real*8, intent(in) :: T
      real*8 :: gamma_deex
      gamma_deex = g10_h_pre*sqrt(max(T, 1.0d0))*exp(-1000.0d0/max(T, 1.0d0))
   end function gamma_10_H

   elemental function gamma_10_H2(T) result(gamma_deex)
      real*8, intent(in) :: T
      real*8 :: gamma_deex
      gamma_deex = g10_h2_pre*sqrt(max(T, 1.0d0))                                 &
          *exp(-12000.0d0/(max(T, 1.0d0) + 1200.0d0))
   end function gamma_10_H2

   ! Critical density of H2 v = 1 against a given collider mix,
   ! n_cr = A/gamma_deex
   ! (Hollenbach & McKee 1979, p. 576) [cm^-3]. Diagnostic: the heat fraction
   ! below does not go through it, so that it stays finite when the collider
   ! density vanishes.
   elemental function h2_vibrational_critical_density(T, f_HI) result(ncr)
      real*8, intent(in) :: T      ! [K]
      real*8, intent(in) :: f_HI   ! atomic-H share of the collider density
      real*8 :: ncr, gamma_deex
      gamma_deex = f_HI*gamma_10_H(T) + (1.0d0 - f_HI)*gamma_10_H2(T)
      ncr = a10_h2/max(gamma_deex, 1.0d-99)
   end function h2_vibrational_critical_density

   ! Fraction of the vibrational excitation energy that becomes HEAT rather
   ! than being radiated away: Hollenbach & McKee's (1 + n_cr/n)^-1, written
   ! with the collider densities explicit so it is finite and correct at
   ! either end,
   !
   !     f = (gamma_H n_HI + gamma_H2 n_H2) / (gamma_H n_HI + gamma_H2 n_H2
   !                                           + A_10) .
   !
   ! Helium is deliberately absent: see the module header.
   elemental function h2_vibrational_heat_fraction(T, n_HI, n_H2) result(f)
      real*8, intent(in) :: T          ! [K]
      real*8, intent(in) :: n_HI, n_H2 ! [cm^-3]
      real*8 :: f, q
      q = gamma_10_H(T)*max(n_HI, 0.0d0) + gamma_10_H2(T)*max(n_H2, 0.0d0)
      f = q/(q + a10_h2)
   end function h2_vibrational_heat_fraction

   ! End of module
   end module h2_vibrational_relaxation

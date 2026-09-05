      module lyman_werner_photodissociation
      ! Photodissociation of H2 in the Lyman and Werner bands (the Solomon
      ! process):
      !
      !     H2(X) + hv(11.2-13.6 eV)  ->  H2(B,C)  ->  H2(X, continuum)
      !                               ->  H + H .
      !
      ! This is the process that destroys H2 above the level where it can
      ! shield itself, and it is the H2 sink the Koskinen et al. (2022)
      ! Table-1 network (mol_rates) leaves out -- their photo-rates start at
      ! the 15.4 eV photoionization edge, so without this term the only H2
      ! losses below 15.4 eV are thermal (R12), electron impact (R14) and
      ! ion chemistry.  Absent it, EXHALE's network keeps a base far more
      ! molecular than photochemical models of the same layer give
      ! (docs/base_composition_handoff_plan.md sec. 5).
      !
      ! ---------------------------------------------------------------
      ! 1. Band and unattenuated rate
      !
      ! Draine & Bertoldi (1996), ApJ 468, 269 (hereafter DB96) work with
      ! the 912-1110 A interval: 912 A is the H Lyman edge (shortward of it
      ! atomic H absorbs everything), and essentially all H2 pumping out of
      ! v = 0 happens longward of 1110 A only through very weak lines
      ! (their footnote 4).  We use the same interval.
      !
      ! DB96 characterize a radiation field by the photon flux in that band
      ! (their eq. [21], F = c n_phot for a beam) and tabulate F and the
      ! unshielded dissociation rate for several spectra.  For the flat-F_lam
      ! spectrum (u_nu ~ nu^-2, their eq. [24], T_color = 2.9e4 K) at chi = 1:
      !
      !     F = 1.208e7 photons cm^-2 s^-1        (their Table 1, the
      !                                            F/chi column)
      !     zeta_pump = 3.09e-10 s^-1, <p_diss> = 0.135   (their Table 2,
      !                                            the T_r = 100 K row)
      !     zeta_diss(0) = 4.17e-11 s^-1          (the annotation inside
      !                                            their Fig. 1 panel; the
      !                                            same number appears in the
      !                                            Fig. 7 panel carrying an
      !                                            explicit chi factor, and
      !                                            in NEITHER caption)
      !
      ! so the dissociation rate per band photon is an effective cross
      ! section
      !
      !     sigma_LW = 4.17e-11 / 1.208e7 = 3.452e-18 cm^2 .
      !
      ! The code takes the band-integrated ENERGY flux at the planet,
      ! F_LW_star [erg cm^-2 s^-1], and converts with the mean photon energy
      ! of a flat F_lam band, <hv> = 2hc/(912 A + 1110 A) = 12.2635 eV, so
      !
      !     k_LW,thin = sigma_LW * F_LW_star / <hv> .
      !
      ! APPROXIMATION AND ITS RANGE.  sigma_LW depends on the shape of the
      ! spectrum WITHIN the band, because the Lyman/Werner lines sample it
      ! unevenly.  Repeating the arithmetic above for the (much softer)
      ! Draine 1978 field of DB96 eq. (23) -- zeta_pump = 2.78e-10,
      ! <p_diss> = 0.119 (again the T_r = 100 K row), F/chi = 1.232e7 --
      ! gives 2.685e-18 cm^2,
      ! 22% lower.  The two spectra bracket color temperatures 1.3e4-2.9e4 K,
      ! and DB96 note (their sec. 3) that PDR properties are insensitive to
      ! the spectrum for T_color >~ 1e4 K.  The flat-F_lam value adopted here
      ! is therefore good to about +-25% for a stellar FUV band that is not
      ! strongly tilted; a spectrum concentrated at one end of the band
      ! (a single dominant emission line) is outside that range.
      !
      ! ---------------------------------------------------------------
      ! 2. Self-shielding: A TABLE, not a fit
      !
      ! The pumping lines saturate, so the rate falls far faster than any
      ! continuum opacity would give.  The suppression factor is taken from
      ! h2_self_shielding_table, which tabulates a LEVEL-RESOLVED CLOUDY
      ! c25.00 calculation on our own (T, n_H, N_H2) grid.  That module's
      ! header states what was calculated, over what range it is valid, and
      ! the one place where it is an upper bound rather than a value.
      !
      ! WHY NEITHER PUBLISHED FIT IS USED FOR THE RATE ANY MORE.  Both were
      ! measured against that calculation (docs/h2_self_shielding_cloudy.md,
      ! 2026-09-02).  Over 1e18 <= N_H2 <= 5e20 -- the columns where our H2
      ! sits, and where the CLOUDY run is still self-consistent -- Draine &
      ! Bertoldi (1996) eq. (37) runs from 1.85x high at 900 K to 2.3x low at
      ! 2700 K, and Richings, Schaye & Oppenheimer (2014) is low by 3.3-4.3x
      ! at every temperature of the layer.  Neither follows the SHAPE of the
      ! level-resolved curve either: it has a trough near 1e16-1e17 and a
      ! shoulder near 1e19 that no two-term algebraic form of this kind
      ! carries.  Replacing one fit by the other would swap one wrong curve
      ! for another, so the calculation itself is tabulated instead.
      !
      ! THIS REVERSES A DECISION.  Update_EXHALE sec. 116 moved the rate from
      ! DB96 to R14 on the strength of R14's own finding that DB96
      ! overestimates the shielding factor by about 3 at T = 100 K.  That
      ! finding is not in doubt -- it is reproduced from our own CLOUDY deck,
      ! DB96/CLOUDY = 2.27-3.38 for 1e17 <= N_H2 <= 5e20 at 100 K.  What the
      ! measurement adds is that the sign of the DB96 error reverses between
      ! 100 K and 900 K, so a fit built to remove the 100 K excess keeps
      ! removing it after it has gone.  Our layer runs at 865-2724 K.
      !
      ! 2a. THE TWO FITS ARE STILL HERE, and are deliberately not deleted:
      ! h2_self_shielding_richings and h2_self_shielding_draine_bertoldi are
      ! kept so that the comparison in docs/h2_self_shielding_cloudy.md can be
      ! reproduced from this tree, and because the SECOND of them is still
      ! live for a different quantity (sec. 3).  NEITHER IS CALLED BY THE
      ! DISSOCIATION RATE.
      !
      ! 2b. DOPPLER PARAMETER: THERMAL ONLY, and now only for the band share.
      ! h2_doppler_parameter returns b = (2kT/m_H2)^(1/2) with no turbulent
      ! term.  The table needs no b argument -- b is thermal, hence a
      ! function of T, which is already one of its axes -- so b is read only
      ! by h2_band_equivalent_width and by the two retained fits.  CLOUDY's
      ! own Doppler width is sqrt(2kT/m + v_turb^2) with v_turb = 0 for the
      ! decks the table was built from, i.e. the same definition, so the
      ! table and this function are consistent.
      !
      ! 2c. NO ELWERT/CLOUDY RADIATION-FIELD EXPONENT, and now measured
      ! rather than argued.  Dropping the incident band flux by four decades
      ! (G_0 from 1.25e6 to 125) moves the level-resolved factor by at most
      ! 9 per cent anywhere on the grid and by under 5 per cent below
      ! N_H2 = 1e20.  There is no radiation-field axis to carry.
      !
      ! 2d. LINE OVERLAP IS INSIDE THE TABLE.  It has to be: overlap is what
      ! DB96 sec. 5.2 say produces "the rapid falloff ... for
      ! N2 > 1e20 cm^-2", and at the columns of a molecular base it is the
      ! dominant term, not a correction -- letting the lines absorb each
      ! other's beam lowers the surviving pumping by 24x at N_H2 = 1e21 and
      ! 630x at 4.4e21 (1300 K, n_H = 1e13).  The table is therefore built
      ! from a line-by-line calculation that carries every Lyman and Werner
      ! transition on one frequency grid, checked against the Meudon PDR code
      ! (Le Petit et al. 2006) with its exact UV line transfer to 5-22 per
      ! cent over 1e19 <= N_H2 <= 1e21.  There is no column above which the
      ! tabulated value becomes an upper bound; the only limit is the top of
      ! the column axis, h2_shield_max_column(), above which the edge value
      ! is returned.  docs/p38_line_overlap_shielding.md is the measurement.
      !
      ! 2e. THE TRAPPING INSIDE THE TABLE IS A SLAB QUANTITY, AND THE PRICE
      ! IS MEASURED.  sigma_diss carries the re-absorption of the fluorescent
      ! decay photons, and that factor alone is computed in a plane-parallel
      ! slab (h2_self_shielding_table header).  A matched slab/sphere
      ! calculation on the same line list gives the size and the sign: the
      ! slab OVER-TRAPS, so this module's rate is high, by 0.2 per cent at
      ! r/H = 41, by 7-10 per cent through an H2-bearing layer, by up to 20
      ! per cent where H2 is already negligible, and by 6.5 per cent weighted
      ! by the dissociation rate.  It is one-signed and it is smaller than
      ! the spread between escape-probability methods in the same slab.
      ! docs/e2_lw_geometry.md.
      !
      ! 2f. WHICH COLUMN, AND WHICH FLUX.  N_H2 below is the column toward
      ! the star, and the caller supplies the RADIAL column
      ! (calc_column_dens_one) -- which for a parallel beam is the substellar
      ! ray exactly, in a sphere as in a slab, so the table is not
      ! approximated by that choice.  What the choice does mean is that the
      ! rate returned here is the SUBSTELLAR rate.  The run-wide dayside
      ! convention that turns it into a shell average is applied to F_LW by
      ! the caller, through fuv_band_flux and global_parameters
      ! dayside_dilution (Update_EXHALE section 150).  The shell average of a
      ! SHIELDED band is not 1/2 but 0.20-0.33 over a hot-Uranus molecular
      ! layer, because the slant columns away from the substellar point are
      ! longer; that refinement is the separate, default-off shell-average
      ! key.
      !
      ! ---------------------------------------------------------------
      ! 3. The DB96 fit, which still owns the band share
      !
      ! DB96 eq. (40) is zeta_diss = f_shield * exp(-tau_dust) *
      ! zeta_diss(0), and their eq. (39) is the closed form of the integral
      ! of eq. (37) -- the summed equivalent width of the pumping lines,
      ! which is what h2_band_equivalent_width returns and what the FUV
      ! band ledger of water_photolysis.f90 needs.  Its derivation
      ! conditions are therefore still live, for that quantity:
      !
      ! DB96 STATE NO UPPER BOUND ON EQ. (37).  Their sec. 5.2 says only that
      ! it "does an excellent job in reproducing the initial rapid decline in
      ! f_shield for 1e14 < N2 < 1e17 cm^-2, the 'square root' behavior for
      ! 1e17 < N2 < 1e20 cm^-2, and the rapid falloff due to line overlap for
      ! N2 > 1e20 cm^-2" (their Figs. 3-5).  What IS bounded is the range over
      ! which they showed their pumping rates against an exact calculation:
      ! "the agreement is excellent even out to the largest column densities
      ! considered, N2 = 3e21 cm^-2, where line overlap effects suppress the
      ! pumping rates by a factor ~10" (same section, of their eq. 30).  So
      ! 3e21 is the largest column the fit was DEMONSTRATED at, not a stated
      ! limit of the fit; a deeper column is undemonstrated rather than
      ! disallowed.  (The simpler power law of their eq. 36 is the one that
      ! carries an explicitly quoted range, and the paper gives two different
      ! ones for it -- 1e15 to 1e21 in the summary, 1e14 to 5e20 in sec. 5.2.)
      !
      ! WHERE THE DB96 FIT WAS DERIVED.  DB96 built eq. (37) -- and hence
      ! its integral eq. (39) -- for interstellar PDR gas: stationary
      ! photodissociation fronts run at n_H = 1e2 cm^-3 with the
      ! temperature parameter T0 = 200 K (their Figs. 3-6 and 8-11; Figs. 1
      ! and 2 add 100 K and 300 K LTE comparisons, and Fig. 7 -- the one this
      ! note leans on below -- is run at T = 1e2 K with no dust), in which the
      ! level distribution of the SHIELDED H2 is not thermal at all but is set
      ! in NON-LTE by UV pumping and by formation on grains.  Three axes of
      ! that derivation matter for a planetary base:
      !  - COLUMN.  Demonstrated to 3e21 cm^-2 (see above).  The hot-Uranus
      !    molecular base of the regression matrix settles at N_H2 = 4.4e21
      !    cm^-2 and passes through higher columns while relaxing, i.e. past
      !    the largest column the fit was checked at, so the deepest cells
      !    are undemonstrated; the run reports the column once when it
      !    leaves the checked range of either fit.
      !  - DOPPLER PARAMETER.  b enters only through the saturated Doppler
      !    core, the (1 + x/b5)^-2 term.  Every exact multiline calculation
      !    DB96 checked the fit against was run at b = 3 km s^-1, so the b
      !    axis of the fit is constructed rather than tested; our thermal
      !    b = 2.79-4.48 km s^-1 sits inside the 1-10 km s^-1 span the form
      !    was built to cover, and brackets their one validated value.
      !  - LEVEL POPULATIONS (T and n_H).  Our layer runs T = 945-2430 K at
      !    n_H up to 9e13 cm^-3, far outside their 200 K and 1e2 cm^-3.
      !    Neither enters the fit as an argument -- T only through b -- so
      !    what they really change is the rovibrational population of the
      !    absorbing H2.  That is exactly the ingredient R14 found DB96 gets
      !    wrong for warm gas (sec. 2a), and it is why the RATE has moved
      !    off this fit.  For the equivalent width DB96's own Fig. 7 is what
      !    bounds the effect: it compares a 100 K LTE distribution against
      !    two UV-pumped non-LTE ones (chi = 1 and 10), and its CAPTION
      !    reads "It is seen that the three self-shielding factors are
      !    essentially identical for N(H2) >~ 2 x 10^17 cm^-2".  Their BODY
      !    text is more careful: "for 1e14 <~ N2 <~ 1e18 cm^-2,
      !    f_shield does show a dependence on the details of the H2
      !    rovibrational distribution.  When there is higher excitation, the
      !    pumping occurs via a larger number of lines, so that the H2 remains
      !    optically thin a bit longer, and self-shielding is less effective
      !    until one enters the heavily damped regime at N2 >~ 1e18 cm^-2.
      !    Overall, however, the self-shielding function f_shield(N2) is
      !    relatively insensitive to the effects of UV pumping with
      !    chi <~ 1e3."  Note that this bounds UV PUMPING, not thermal
      !    rotational excitation, which is the axis R14 measured; the band
      !    share carries that unquantified, and it is the weakest link left
      !    in the oxygen-chemistry beam split.
      !
      ! ---------------------------------------------------------------
      ! 4. What shields the band besides the H2 lines
      !  - Dust: EXHALE's metals are atomic and trace, so there are no
      !    grains, and the dust part of DB96 eq. (40) is identically 1.
      !    Dust is not inside f_shield either -- DB96 keep it in the separate
      !    exp(-tau_d) factor -- so a dust-free use of the fit is the use it
      !    was written for.
      !  - H2O and OH, when the oxygen chemistry carries them: they absorb
      !    912-1110 A as a CONTINUUM, which is exactly the absorber
      !    exp(-tau) of DB96 eq. (40) stands for.  That term is therefore
      !    not zero with the option on, and it is passed in as tau_cont
      !    below.  It is zero for every run without the oxygen chemistry,
      !    which is why this module could ignore it before that option
      !    existed.  The photons the H2 lines take out of the band are
      !    removed from what H2O and OH see by the same token; the other
      !    half of that bookkeeping is water_photolysis.f90 sec. 3, and
      !    h2_band_equivalent_width below is what it needs from here.
      !  - Trace-metal continuum: the neutral low-IP metals (Mg I, Fe I,
      !    Si I, Ca I, Na I, K I) do photoionize inside the band, but at
      !    solar abundance and sigma ~ 1e-18 cm^2 their optical depth is
      !    ~4e-23 N_H, i.e. <~ 0.05 at the base column where f_shield is
      !    already below 1e-4.  Neglected.
      !  - H Lyman-series lines (Ly-beta 1025.7, Ly-gamma 972.5, ... all lie
      !    in the band): DB96 include them in the equivalent width their fit
      !    was built on and state (their sec. 4.3) that "absorption by the H
      !    Lyman lines has only a small effect on the H2 pumping rates".
      !    They are not treated separately here.  In a wind the H I / H2
      !    column ratio is much larger than in the PDRs DB96 fitted, so this
      !    is the least controlled of the three; it can only reduce the rate.
      !
      ! ---------------------------------------------------------------
      ! 5. Heating
      !
      ! The fragments of a Solomon-process dissociation leave with kinetic
      ! energy.  Black & Dalgarno (1977), ApJS 34, 405, p. 418: "Fluorescent
      ! dissociation of H2 gives rise to a pair of energetic hydrogen atoms
      ! (Milgrom, Panagia, and Salpeter 1973; Stephens and Dalgarno 1973);
      ! for a typical ultraviolet radiation field, the yield is about 0.4 eV
      ! per atom pair, but it varies slightly with depth."  We adopt the
      ! 0.4 eV, held constant with depth.
      !
      ! The 4.48 eV H-H bond energy is paid by the absorbed photon, not by
      ! the gas, so it is NOT a thermal sink of this channel.  (A thermal
      ! dissociation-energy sink for R12/R14 is a separate open item.)
      !
      ! ---------------------------------------------------------------
      ! References: h2_self_shielding_table.f90 and
      ! docs/h2_self_shielding_cloudy.md (the level-resolved calculation the
      ! rate now uses, and the measurement that put it there);
      ! Draine & Bertoldi (1996) ApJ 468, 269 (published
      ! version, eqs. 20, 21, 24, 37, 39, 40, Tables 1-2, Fig. 7 caption);
      ! Richings, Schaye & Oppenheimer (2014) MNRAS 442, 2780 (published
      ! version, sec. 3.2 eqs. 3.11-3.17, Appendix B and the Fig. B1
      ! caption); Wolcott-Green, Haiman & Bryan (2011) MNRAS 418, 838, as
      ! quoted by R14; Black & Dalgarno (1977) ApJS 34, 405, p. 418.

      use h2_self_shielding_table, only:                                 &
                h2_lw_dissociation_cross_section,                        &
                h2_lw_dissociation_per_pump,                             &
                h2_lw_dissociation_per_absorbed_photon,                  &
                h2_self_shielding_level_resolved,                        &
                h2_shield_max_column

      implicit none
      private
      ! Re-exported from h2_self_shielding_table so that a caller of the
      ! Lyman-Werner rate needs one use statement, not two.
      public :: h2_self_shielding_level_resolved,                        &
                h2_shield_max_column
      public :: h2_doppler_parameter,                                    &
                h2_self_shielding_richings,                              &
                h2_self_shielding_draine_bertoldi,                       &
                h2_band_equivalent_width,                                &
                lyman_werner_dissociation_rate, e_lw_fragment_erg,       &
                e_lw_photon_erg,                                         &
                h2_lw_dissociation_cross_section,                        &
                h2_lw_dissociation_per_pump,                             &
                h2_lw_dissociation_per_absorbed_photon

      ! Mean photon energy of a flat-F_lam 912-1110 A band [erg]:
      ! 2hc/(912+1110 A) = 12.2635 eV.
      real*8, parameter :: e_lw_photon_erg = 1.96483d-11

      ! Kinetic energy released to the H + H pair, 0.4 eV in erg
      ! (Black & Dalgarno 1977, p. 418).
      real*8, parameter :: e_lw_fragment_erg = 6.40871d-13

      ! Boltzmann constant and the H2 mass, in the same cgs values the rest
      ! of the code uses (global_parameters kb_erg, mu); repeated here so
      ! the module has no dependency on the parameter block. Any change to
      ! kb_erg or mu there must be mirrored here (kb_erg -> CODATA
      ! 1.380649e-16 on 2026-08-15, synced below).
      real*8, parameter :: kb_lw = 1.380649d-16
      real*8, parameter :: m_h2  = 2.0d0*1.67353284d-24

      contains

      ! Thermal Doppler parameter of H2, b = (2 k T / m_H2)^(1/2) [cm s^-1].
      ! Both self-shielding fits take b from here, so the thermal-only
      ! choice is made once.
      !
      ! No turbulent contribution.  Two reasons, and the second is the
      ! binding one.  (i) The model has no sub-grid velocity field to set
      ! b_turb from.  (ii) R14 fitted their function to CLOUDY for purely
      ! thermal broadening, and their Fig. B1 caption reads: "The agreement
      ! between our best-fitting self-shielding function and CLOUDY is
      ! poorer when turbulence is included, as it was fitted to the purely
      ! thermal Doppler broadening case."  Adding the b_turb = 7.1 km s^-1
      ! they use for interstellar gas (their eqs. 3.16-3.17) would move us
      ! OUT of the calibration of the fit we call, not into a better
      ! description of it.
      !
      ! b enters both fits only through the saturated line core, where a
      ! larger b would raise the rate.
      double precision function h2_doppler_parameter(T) result(b)
      real*8, intent(in) :: T
      b = sqrt(2.0d0*kb_lw*max(T, 1.0d0)/m_h2)
      end function h2_doppler_parameter

      ! Richings, Schaye & Oppenheimer (2014) eqs. (3.12)-(3.15): their fitted
      ! H2 self-shielding factor for a star-ward H2 column N_H2 [cm^-2] at gas
      ! temperature T [K] and Doppler parameter b [cm s^-1].
      !
      ! NOT CALLED BY THE RATE.  It is retained so that the comparison against
      ! the level-resolved calculation (docs/h2_self_shielding_cloudy.md) can
      ! be reproduced from this tree; header sec. 2 gives the measurement that
      ! took the rate off it.  Over 1e18 <= N_H2 <= 5e20 it is low by 3.3-4.3x
      ! against that calculation at 900-2700 K.
      !
      ! b is the caller's, so that the thermal-only choice (header sec. 2b)
      ! is made in one place, h2_doppler_parameter.
      !
      ! T enters TWICE and the two are physically distinct: through b, the
      ! Doppler width of one line, and through w/alpha/N_crit, the
      ! rovibrational level populations that decide how many lines there
      ! are to pump.  The second is what DB96 has no handle on and what R14
      ! fitted against CLOUDY.
      !
      ! Validity as published (R14 sec. 3.2): within 30 per cent at 100 K
      ! for N_H2 < 10^21 cm^-2, within 60 per cent at 5000 K for
      ! N_H2 < 10^20 cm^-2.  Our layer is outside both.
      double precision function h2_self_shielding_richings(N_H2, T, b)    &
                                result(S_self)
      real*8, intent(in) :: N_H2, T, b
      real*8 :: x, b5, s, w_h2, alpha_t, n_crit, t_gas, t_cut
      t_gas = max(T, 1.0d0)

      ! R14 eq. (3.13): the weight of the damping-wing term.  The
      ! exp[-(T/3900 K)^14.6] cutoff is a fitted shape, not an asymptotic
      ! form, so it is evaluated as written.  Above T = 6.1e3 K
      ! (T/3900 K > 1.57) that exponent already exceeds 700 and the term
      ! has gone; the branch keeps the 14.6th power from overflowing on the
      ! way to a zero it is going to reach anyway.
      t_cut = t_gas/3900.0d0
      if (t_cut .gt. 1.57d0) then
         w_h2 = 0.0d0
      else
         w_h2 = 0.013d0*(1.0d0 + (t_gas/2700.0d0)**1.3d0)**(1.0d0/1.3d0)  &
              * exp(-t_cut**14.6d0)
      endif

      ! R14 eq. (3.14).
      if (t_gas .lt. 3.0d3) then
         alpha_t = 1.4d0
      else if (t_gas .lt. 4.0d3) then
         alpha_t = (t_gas/4500.0d0)**(-0.8d0)
      else
         alpha_t = 1.1d0
      endif

      ! R14 eq. (3.15).  N_crit replaces the fixed 5e14 cm^-2 of DB96: it
      ! is the column at which the pumping lines saturate, and warmer gas
      ! spreads the pumping over more lines, so it rises with T below
      ! 3000 K.
      if (t_gas .lt. 3.0d3) then
         n_crit = 1.3d14*(1.0d0 + (t_gas/600.0d0)**0.8d0)
      else if (t_gas .lt. 4.0d3) then
         n_crit = 1.0d14*(t_gas/4760.0d0)**(-3.8d0)
      else
         n_crit = 2.0d14
      endif

      x  = max(N_H2, 0.0d0)/n_crit
      b5 = max(b, 1.0d0)/1.0d5
      s  = sqrt(1.0d0 + x)
      S_self = (1.0d0 - w_h2)/(1.0d0 + x/b5)**alpha_t                     &
             * exp(-5.0d-7*(1.0d0 + x))                                   &
             + w_h2/s*exp(-8.5d-4*s)
      end function h2_self_shielding_richings

      ! DB96 eq. (37): H2 self-shielding factor for a star-ward H2 column
      ! N_H2 [cm^-2] and Doppler parameter b [cm s^-1].  Demonstrated against
      ! exact calculations up to N_H2 = 3e21 cm^-2, under the conditions
      ! sec. 3 of the header sets out; it tends to 1 as N_H2 -> 0 and falls as
      ! N_H2^-3/4 in the saturated regime, steeper than N_H2^-1/2 because of
      ! line overlap.
      !
      ! THIS IS NOT THE RATE'S SHIELDING FUNCTION.  It is retained because
      ! h2_band_equivalent_width below is the closed-form integral of THIS
      ! expression and of no other, and because that integral's
      ! normalization is exact against DB96's own tables -- and is confirmed
      ! independently by the same level-resolved calculation the rate now
      ! uses: sigma_pump measured off it at 1300 K is 2.610e-17 cm^2 against
      ! the 2.557e-17 that DB96's own tables give (the normalization is set
      ! out in the h2_band_equivalent_width header below), a 2 per cent
      ! agreement.  The rate is on h2_lw_dissociation_cross_section
      ! (header sec. 2).
      double precision function h2_self_shielding_draine_bertoldi(N_H2, b)&
                                result(f_shield)
      real*8, intent(in) :: N_H2, b
      real*8 :: x, b5, s
      x  = max(N_H2, 0.0d0)/5.0d14
      b5 = max(b, 1.0d0)/1.0d5
      s  = sqrt(1.0d0 + x)
      f_shield = 0.965d0/(1.0d0 + x/b5)**2                                &
               + 0.035d0/s*exp(-8.5d-4*s)
      end function h2_self_shielding_draine_bertoldi

      ! Fraction of the 912-1110 A band that the H2 Lyman and Werner lines
      ! have taken out of the beam by the time it has crossed a star-ward H2
      ! column N_H2 [cm^-2] at Doppler parameter b [cm s^-1].
      !
      ! DB96 integrate their own f_shield analytically.  Their eq. (38)
      ! defines the total dimensionless equivalent width of the pumping
      ! lines, W(N2) = Dln_nu (zeta_pump(0)/F) int_0^N2 f_shield dN2', and
      ! their eq. (39) is the closed form
      !
      !   W(N2) = Dln_nu 1.05 [ 1 + 0.0117 x/(1 + x/b5)
      !                         - exp( -8.5e-4 ((1+x)^0.5 - 1) ) ] ,
      !   x = N2/5e14 cm^-2 ,   b5 = b/1e5 cm s^-1        (their eq. 37) .
      !
      ! Dln_nu = ln(1110/912) = 0.1964753 is the logarithmic width of the
      ! band, so W is a width in ln(nu) and W/Dln_nu -- what this function
      ! returns -- is the fraction of the band the lines occupy.
      !
      ! WHAT THE FIT IS NORMALIZED TO.  DB96 built eq. (39) on their own
      ! Table 1 and Table 2 at the flat-F_lam spectrum of their eq. (24),
      ! chi = 1, and the T_r = 100 K row: band photon flux
      ! F = 1.208e7 cm^-2 s^-1, pump rate zeta_pump = 3.09e-10 s^-1, and
      ! <p_diss> = 0.135, so that the dissociation cross section is
      ! 3.452e-18 cm^2 and the PUMP cross section -- the one this equivalent
      ! width is a statement about, since every pump removes a photon
      ! whether or not it dissociates -- is
      !
      !     sigma_pump = 3.452e-18/0.135 = 2.557e-17 cm^2 .
      !
      ! THESE NUMBERS ARE NOT COMPUTED WITH ANYWHERE.  They are the source
      ! of the fit evaluated below and the check on it: the level-resolved
      ! CLOUDY calculation gives 2.610e-17 cm^2 for the same quantity at
      ! 1300 K, 2 per cent away, and reproduces the 3.452e-18 to 0.5 per
      ! cent in its untrapped limit (docs/p39_lw_cross_section_sources.md).
      ! The photodissociation RATE and both branching ratios come from
      ! h2_self_shielding_table, not from here.  Because it
      ! is a summed equivalent width and not an optical depth, the beam's
      ! line transmission is 1 - A -- exactly, for non-overlapping lines
      ! across a flat band; exp(-A) would be wrong.
      !
      ! LINE OVERLAP IS ALREADY INSIDE IT.  f_shield was constructed on the
      ! overlap-corrected pumping rates of DB96 eq. (30), i.e. their eq. (29)
      ! with W_max = ln(1110/912) ~ 0.2.  Applying a further overlap factor
      ! to A would count the same suppression twice.
      !
      ! LIMITS.  The bracket vanishes exactly at x = 0, so A = 0 where there
      ! is no H2, and tends to 1 + 0.0117 b5 as x -> infinity, i.e.
      ! A -> 1.05(1 + 0.0117 b5) = 1.084 (b = 2.79 km s^-1) to 1.105
      ! (4.48 km s^-1) over the Doppler parameters of our layer.  An
      ! asymptote above 1 is slack in the authors' own fit, not a defect
      ! here: DB96 cap the overlap-corrected equivalent width at
      ! W_max = ln(1110/912) (their eq. 29 and sec. 4.3), i.e. at A = 1, and
      ! remark below their eq. (39) only that it "corresponds to an
      ! equivalent width W ~ Dln_nu ~ 0.2 in the limit N2 -> infinity".  So
      ! the deep layer, where A is within a few per cent of its asymptote,
      ! carries an 8-11 per cent overshoot of DB96's own ceiling.  That is a
      ! property of the fit and is therefore the same whether A comes from
      ! eq. (39) or from integrating f_shield on the grid.  The caller clamps
      ! A at 1 so the transmission it hands the continuum absorbers of the
      ! same interval cannot go negative.
      double precision function h2_band_equivalent_width(N_H2, b)         &
                                result(A)
      real*8, intent(in) :: N_H2, b
      real*8 :: x, b5, s
      x  = max(N_H2, 0.0d0)/5.0d14
      b5 = max(b, 1.0d0)/1.0d5
      s  = sqrt(1.0d0 + x)
      A  = 1.05d0*(1.0d0 + 0.0117d0*x/(1.0d0 + x/b5)                     &
                 - exp(-8.5d-4*(s - 1.0d0)))
      end function h2_band_equivalent_width

      ! Photodissociation rate of H2 [s^-1] for a band-integrated stellar
      ! energy flux F_LW [erg cm^-2 s^-1] at the planet, a star-ward H2
      ! column N_H2 [cm^-2], a gas temperature T [K], a total hydrogen
      ! nucleus density n_H [cm^-3] and the star-ward CONTINUUM optical depth
      ! tau_cont of the same 912-1110 A interval.
      !
      ! This is DB96 eq. (40) in structure -- a line self-shielding factor,
      ! exp(-tau) for the continuum, no dust term -- with the level-resolved
      ! CLOUDY table in place of a closed-form fit (header sec. 2).  n_H
      ! enters only through that table, which carries a density axis because
      ! it can; the measured density dependence is 3-7 per cent per two
      ! decades below N_H2 = 1e19 and 40 per cent at 1e21.
      !
      ! tau_cont is the H2O + OH continuum of the oxygen chemistry
      ! (water_photolysis.f90) and is identically zero for a run without it,
      ! so a molecular run with no oxygen chemistry sees only the shielding.
      double precision function lyman_werner_dissociation_rate(F_LW,     &
                                N_H2, T, n_H, tau_cont) result(k)
      real*8, intent(in) :: F_LW, N_H2, T, n_H, tau_cont
      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif
      ! The cross section is per INCIDENT band photon and already carries the
      ! self-shielding of the lines and the trapping of the fluorescent decay
      ! photons; the continuum of the same interval is the caller's tau_cont.
      k = F_LW/e_lw_photon_erg                                           &
        * h2_lw_dissociation_cross_section(N_H2, T, n_H)
      if (tau_cont .gt. 0.0d0) k = k*exp(-tau_cont)
      end function lyman_werner_dissociation_rate

      ! End of module
      end module lyman_werner_photodissociation

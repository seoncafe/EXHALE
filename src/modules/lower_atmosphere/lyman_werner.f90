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
      ! molecular than photochemical models of the same layer give.
      !
      ! ---------------------------------------------------------------
      ! 1. Band and unattenuated rate
      !
      ! THE BAND IS 912-1201 A.  912 A is the H Lyman edge, shortward of
      ! which atomic H absorbs everything; 1201 A is where band B2 of the
      ! FUV band list begins (oxygen_rates.f90), and it is also the red end
      ! of the Lyman and Werner line list the self-shielding table is built
      ! from.
      !
      ! WHY NOT DRAINE & BERTOLDI'S 912-1110 A.  Their footnote 4 says that
      ! essentially all H2 pumping out of v = 0 happens longward of 1110 A
      ! only through very weak lines, and in interstellar gas v = 0 is the
      ! only level populated.  A planetary molecular base runs at
      ! 700-3200 K, where the vibrationally excited levels are populated and
      ! pump in lines that lie longward of 1110 A; at saturation those lines
      ! carry about a third of the dissociations the table rates
      ! (md/p38_line_overlap_shielding.md sec. 4.4).  The band ran to
      ! 1110 A until 2026-09-06 while the table's line list ran past it, and
      ! the table then rated 45 per cent more line absorptions than the beam
      ! lost (sec. 3a).  One band, one line list, one normalization.
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
      !     sigma_LW = 4.17e-11 / 1.208e7 = 3.452e-18 cm^2
      !
      ! per photon of THEIR band.  Their numbers are not the rate: the rate
      ! is the level-resolved table of sec. 2, normalized per photon of
      ! 912-1201 A.  A flat F_lambda carries 1.52529 times as many photons
      ! in 912-1201 A as in 912-1110 A (the ratio of the two int lambda
      ! dlambda), so DB96's cross section restated per 912-1201 A photon is
      ! 2.263e-18 cm^2, which is what it has to be compared against.
      !
      ! The code takes the band-integrated ENERGY flux at the planet,
      ! F_LW_star [erg cm^-2 s^-1], and converts with the mean photon energy
      ! of a flat F_lam band, <hv> = 2hc/(912 A + 1201 A) = 11.7354 eV, so
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
      ! h2_self_shielding_table, which tabulates a LINE-BY-LINE calculation
      ! over every Lyman and Werner transition on one frequency grid, so
      ! that line overlap is inside it, on our own (T, n_H, N_H2) grid; the
      ! one ingredient it takes from level-resolved CLOUDY c25.00 runs is
      ! the trapping ratio p_eff/p_single of the fluorescent decay photons.
      ! That module's header states what was calculated, over what range it
      ! is valid, and which of its ingredients is computed in slab geometry.
      !
      ! WHY NEITHER PUBLISHED FIT IS USED FOR THE RATE.  Both were measured
      ! against the level-resolved CLOUDY c25.00 calculation on the same
      ! grid.  Over 1e18 <= N_H2 <= 5e20 -- the columns where our H2
      ! sits, and where the CLOUDY run is still self-consistent -- Draine &
      ! Bertoldi (1996) eq. (37) runs from 1.85x high at 900 K to 2.3x low at
      ! 2700 K, and Richings, Schaye & Oppenheimer (2014; R14 below, not
      ! the R14 electron-impact reaction of mol_rates) is low by 3.3-4.3x
      ! at every temperature of the layer.  Neither follows the SHAPE of the
      ! level-resolved curve either: it has a trough near 1e16-1e17 and a
      ! shoulder near 1e19 that no two-term algebraic form of this kind
      ! carries.  Replacing one fit by the other would swap one wrong curve
      ! for another, so the calculation itself is tabulated instead.
      !
      ! THIS REVERSES A DECISION.  Update_EXHALE_stage1 sec. 116 moved the rate from
      ! DB96 to R14 on the strength of R14's own finding that DB96
      ! overestimates the shielding factor by about 3 at T = 100 K.  That
      ! finding is not in doubt -- it is reproduced from our own CLOUDY deck,
      ! DB96/CLOUDY = 2.27-3.38 for 1e17 <= N_H2 <= 5e20 at 100 K.  What the
      ! measurement adds is that the sign of the DB96 error reverses between
      ! 100 K and 900 K, so a fit built to remove the 100 K excess keeps
      ! removing it after it has gone.  Our layer runs at 865-2724 K.
      !
      ! 2a. NEITHER FIT IS CARRIED.  The Draine & Bertoldi (1996) eq. (37)
      ! and Richings, Schaye & Oppenheimer (2014) eqs. (3.12)-(3.15)
      ! self-shielding functions are not in this module: the rate, the band
      ! share and the photon ledger all read the table.
      !
      ! 2b. DOPPLER PARAMETER: THERMAL ONLY.  The table needs no b argument:
      ! its line-by-line calculation takes b = (2kT/m_H2)^(1/2) with no
      ! turbulent term, so b is a function of T, which is already one of its
      ! axes.  CLOUDY's own Doppler width, which sets the trapping ratio, is
      ! sqrt(2kT/m + v_turb^2) with v_turb = 0 for the decks the table was
      ! built from, i.e. the same definition.
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
      ! is returned.  md/p38_line_overlap_shielding.md is the measurement.
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
      !
      ! 2f. WHICH COLUMN, AND WHICH FLUX.  N_H2 below is the column toward
      ! the star, and the caller supplies the RADIAL column
      ! (calc_column_dens_one) -- which for a parallel beam is the substellar
      ! ray exactly, in a sphere as in a slab, so the table is not
      ! approximated by that choice.  What the choice does mean is that the
      ! rate returned here is the SUBSTELLAR rate.  The run-wide dayside
      ! convention that turns it into a shell average is applied to F_LW by
      ! the caller, through fuv_band_flux and global_parameters
      ! dayside_dilution (Update_EXHALE_stage1 section 150).  The shell average of a
      ! SHIELDED band is not 1/2 but 0.20-0.33 over a hot-Uranus molecular
      ! layer, because the slant columns away from the substellar point are
      ! longer; that refinement is NOT IMPLEMENTED, and no key asks for it,
      ! so the dilution applied here is the optically thin one.
      !
      ! ---------------------------------------------------------------
      ! 2g. THE BEAM LOSES WHAT THE TABLE RATES, AND IT IS ONE NUMBER
      !
      ! The photons a beam loses to the lines over a star-ward column N are
      ! N_b A(N) with
      !
      !     A(N) = int_0^N sigma_pump dN' ,   sigma_pump = sigma_diss/p_eff ,
      !
      ! because p_eff is the dissociations per photon REMOVED FROM THE BEAM.
      ! The transmission of the beam past the lines is 1 - A and not
      ! exp(-A): the lines occupy a SHARE of the band, they do not attenuate
      ! it uniformly.  A is served from the same table as the rate,
      ! h2_lw_band_photon_fraction_absorbed, as the exact column integral of
      ! that table's own interpolant, so the share the continuum absorbers
      ! of the same interval are denied is exactly the share the
      ! dissociation rate spends.  The FUV band ledger of write_output.f90
      ! closes in this band for that reason and not by a tolerance.
      !
      ! A CANNOT PASS 1, and the table respects it: MEASURED over the whole
      ! (T, n_H) grid, A at the top of the column axis (5e21 cm^-2) is
      ! 0.6903 to 0.9974 (h2_self_shielding_table).  Deeper than that axis
      ! the table is clamped at its edge cross section and A goes on growing
      ! linearly, so the caller caps it at 1 (util_ion_eq.f90); a capped A
      ! means the lines have taken the band.
      !
      ! 3. WHY THE BAND SHARE IS NOT DB96 EQ. (39)
      !
      ! DB96's closed-form integral of their eq. (37), their eq. (39), is an
      ! equivalent width of 912-1110 A carrying interstellar level
      ! populations, and the share it returns disagrees with the table's own
      ! column integral by a factor 1.4 at N_H2 = 1e21 cm^-2 and 1300 K and
      ! 6.6 at 1e18 cm^-2 and 2700 K, with a sign change between them.  The
      ! band share A is therefore the table's column integral (sec. 2g).
      !
      ! 3a. WHAT WAS MEASURED, AND HOW IT WAS CLOSED.  The 45 per cent
      ! excess of the FUV band ledger was two WAVELENGTH BANDS, not two
      ! fits.  MEASURED (2026-09-06 with the table then shipped): the column
      ! integral
      ! of sigma_diss/p_eff reached 1.5154 at the top of the column axis,
      ! where a share of the band cannot pass 1, and 1.5154 is the photon
      ! content of 912-1200 A over that of 912-1110 A for the flat-F_lambda
      ! spectrum the table is built on, 1.5193.  The table sampled
      ! 911.75-1200 A while being normalized per photon of 912-1110 A, and
      ! the beam lost the DB96 eq. (39) equivalent width of 912-1110 A
      ! alone.  On the oxygen_chemistry state the ledger rated 1.4544 line
      ! photons for every one the beam lost, and 249.4 erg cm^-2 s^-1 of
      ! band energy against an incident 171.5.
      !
      ! Both halves were repaired at once and in one direction, which is the
      ! direction that keeps every absorber the line list has: the band is
      ! now 912-1201 A everywhere (the FUV band list, the incident flux key,
      ! the table's normalization and its line list), and the beam's loss is
      ! the table's own pump absorption instead of eq. (39).
      ! md/p38_line_overlap_shielding.md sec. 4.4 states what the wider
      ! band is worth to the RATE: f_shield at 1300 K is 1.03x larger at
      ! N_H2 = 1e19 cm^-2 and 3.7x larger at 4.4e21 than on 912-1110 A.
      !
      ! 3b. WHAT THE DB96 NORMALIZATION IS STILL GOOD FOR.  DB96 built their
      ! Tables 1 and 2 at the flat-F_lam spectrum of their eq. (24), chi = 1
      ! and the T_r = 100 K row: band photon flux F = 1.208e7 cm^-2 s^-1,
      ! pump rate zeta_pump = 3.09e-10 s^-1 and <p_diss> = 0.135, so their
      ! PUMP cross section is 3.452e-18/0.135 = 2.557e-17 cm^2 per
      ! 912-1110 A photon, i.e. 1.676e-17 cm^2 per 912-1201 A photon.  The
      ! level-resolved table gives 1.7325e-17 cm^2 at 1300 K at the bottom
      ! of its column axis, 3 per cent away.  THAT IS AN AGREEMENT OF THE
      ! THIN LIMIT ONLY, i.e. of the SLOPE of A at zero column, and it is
      ! the only place the two normalizations were ever comparable.
      !
      ! 3c. WHERE THE DB96 FIT WAS DERIVED, AND WHAT IT COULD AND COULD NOT
      ! CARRY.
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
      !    would be undemonstrated for it.  No rate or band share of the
      !    code is on this fit any more (sec. 3); the run reports once when
      !    the column passes the top of the column axis of the line-by-line
      !    table instead (ionization_equilibrium).
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
      !    exp(-tau_d) factor -- and the table carries none.
      !  - H2O and OH, when the oxygen chemistry carries them: they absorb
      !    912-1201 A as a CONTINUUM, which is exactly the absorber
      !    exp(-tau) of DB96 eq. (40) stands for.  That term is therefore
      !    not zero with the option on, and it is passed in as tau_cont
      !    below.  It is zero for every run without the oxygen chemistry,
      !    which is why this module could ignore it before that option
      !    existed.  The photons the H2 lines take out of the band are
      !    removed from what H2O and OH see by the same token; the other
      !    half of that bookkeeping is water_photolysis.f90 sec. 3, and
      !    h2_lw_band_photon_fraction_absorbed is what it needs.
      !  - Trace-metal continuum: the neutral low-IP metals (Mg I, Fe I,
      !    Si I, Ca I, Na I, K I) do photoionize inside the band, but at
      !    solar abundance and sigma ~ 1e-18 cm^2 their optical depth is
      !    ~4e-23 N_H, i.e. <~ 0.05 at the base column where f_shield is
      !    already below 1e-4.  Neglected.
      !  - H Lyman-series lines (Ly-beta 1025.7, Ly-gamma 972.5, ... all lie
      !    in the band): DB96 include them in the equivalent width their fit
      !    was built on and state (their sec. 4.3) that "absorption by the H
      !    Lyman lines has only a small effect on the H2 pumping rates".
      !    They are not treated separately here, and the line-by-line
      !    table carries the H2 lines alone, so they are in neither the rate
      !    nor the band share.  In a wind the H I / H2 column ratio is much
      !    larger than in the PDRs DB96 fitted, so this is the least
      !    controlled of the three; it can only reduce the rate.
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
      ! CONFIRMED against the level-resolved values: Abgrall, Roueff &
      ! Drira (2000), A&AS 141, 297, give the mean kinetic energy of the
      ! dissociating products for every rovibronic level of the B, C, B' and
      ! D states, and weighting those by this module's own rate weights
      ! gives 0.397 eV at 100 K, 0.406 at 1300 K and 0.429 at 2700 K
      ! (MEASURED, md/p39_lw_cross_section_sources.md section 4.1), so the
      ! constant is right to 2 per cent over the layer.
      !
      ! The 4.48 eV H-H bond energy is paid by the absorbed photon, not by
      ! the gas, so it is NOT a thermal sink of this channel.  (A thermal
      ! dissociation-energy sink for R12/R14 is a separate open item.)
      !
      ! ---------------------------------------------------------------
      ! References: h2_self_shielding_table.f90 (the line-by-line
      ! calculation the rate uses, with the trapping ratio from level-
      ! resolved CLOUDY c25.00 runs);
      ! Draine & Bertoldi (1996) ApJ 468, 269 (published
      ! version, eqs. 20, 21, 24, 37, 39, 40, Tables 1-2, Fig. 7 caption);
      ! Richings, Schaye & Oppenheimer (2014) MNRAS 442, 2780 (published
      ! version, sec. 3.2 eqs. 3.11-3.17 and Appendix B); Black & Dalgarno
      ! (1977) ApJS 34, 405, p. 418.

      use h2_self_shielding_table, only:                                 &
                h2_lw_dissociation_cross_section,                        &
                h2_lw_dissociation_per_pump,                             &
                h2_lw_dissociation_per_absorbed_photon,                  &
                h2_lw_pump_cross_section,                                &
                h2_lw_band_photon_fraction_absorbed,                     &
                h2_self_shielding_level_resolved,                        &
                h2_shield_max_column

      implicit none
      private
      ! Re-exported from h2_self_shielding_table so that a caller of the
      ! Lyman-Werner rate needs one use statement, not two.
      public :: h2_self_shielding_level_resolved,                        &
                h2_shield_max_column,                                    &
                h2_lw_pump_cross_section,                                &
                h2_lw_band_photon_fraction_absorbed
      public :: lyman_werner_dissociation_rate,                          &
                lyman_werner_dissociation_rate_cell_mean,                &
                lyman_werner_band_absorption_rate_cell_mean,             &
                e_lw_fragment_erg,                                       &
                e_lw_photon_erg,                                         &
                h2_lw_dissociation_cross_section,                        &
                h2_lw_dissociation_per_pump,                             &
                h2_lw_dissociation_per_absorbed_photon

      ! Mean photon energy of a flat-F_lam 912-1201 A band [erg]:
      ! 2hc/(912+1201 A) = 11.7354 eV.  It is the exact photon content of
      ! the band for a flat F_lambda, so F_LW/e_lw_photon_erg is the band
      ! photon flux; the self-shielding table is normalized per photon of
      ! the same band with the same number (src/utils/h2_shielding_lbl).
      real*8, parameter :: e_lw_photon_erg = 1.88021d-11

      ! Kinetic energy released to the H + H pair, 0.4 eV in erg
      ! (Black & Dalgarno 1977, p. 418; confirmed to 2 per cent against the
      ! level-resolved values of Abgrall, Roueff & Drira 2000, section 5 of
      ! the header).
      real*8, parameter :: e_lw_fragment_erg = 6.40871d-13

      ! The shape of every cross section h2_self_shielding_table serves:
      ! a function of the star-ward H2 column, the gas temperature and the
      ! hydrogen nucleus density.  The band rate and its cell mean are
      ! written once against this shape and are then given sigma_diss or
      ! sigma_pump by their callers.
      abstract interface
         double precision function h2_cross_section_of_state(N_H2, T, n_H)
         real*8, intent(in) :: N_H2, T, n_H
         end function h2_cross_section_of_state
      end interface

      contains

      ! LOCAL photodissociation rate of H2 [s^-1] at ONE POINT of the
      ! column: a band-integrated stellar energy flux F_LW
      ! [erg cm^-2 s^-1] at the planet, a star-ward H2 column N_H2 [cm^-2],
      ! a gas temperature T [K], a total hydrogen nucleus density n_H
      ! [cm^-3] and the star-ward CONTINUUM optical depth tau_cont of the
      ! same 912-1201 A interval, all at that point.  What a grid cell needs
      ! is the mean of this over the cell, which is the function below.
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
      k = band_rate_from_cross_section(h2_lw_dissociation_cross_section, &
                                       F_LW, N_H2, T, n_H, tau_cont)
      end function lyman_werner_dissociation_rate

      ! The same contraction for any of the table's cross sections of
      ! (N_H2, T, n_H): the incident band photon fluence F_LW/<hv> times the
      ! cross section, attenuated by the continuum of the same interval.
      ! sigma_diss gives the dissociation rate, sigma_pump the rate at which
      ! the lines take photons out of the beam.  One expression, so the two
      ! cannot drift apart.
      double precision function band_rate_from_cross_section(sigma_at,   &
                                F_LW, N_H2, T, n_H, tau_cont) result(k)
      procedure(h2_cross_section_of_state) :: sigma_at
      real*8, intent(in) :: F_LW, N_H2, T, n_H, tau_cont
      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif
      ! The cross section is per INCIDENT band photon and already carries the
      ! self-shielding of the lines; the continuum of the same interval is
      ! the caller's tau_cont.
      k = F_LW/e_lw_photon_erg*sigma_at(N_H2, T, n_H)
      if (tau_cont .gt. 0.0d0) k = k*exp(-tau_cont)
      end function band_rate_from_cross_section

      ! MEAN of that rate over one grid cell [s^-1], for the star-ward face
      ! values (N_H2_out, tau_out) and the inner face values (N_H2_in,
      ! tau_in) of the H2 column and of the 912-1201 A continuum depth.
      !
      ! WHY THE MEAN AND NOT A FACE VALUE.  The dissociation cross section
      ! falls by four decades across the self-shielding transition, and it
      ! falls fastest at the H2 front, exactly where one cell can carry a
      ! large fraction of a decade of H2 column.  A rate taken at one face
      ! and applied to the whole cell is then wrong in one direction
      ! everywhere in that cell.  The H2O and OH continua of the SAME beam
      ! already take the exact cell mean (water_photolysis.f90 sec. 3,
      ! where the one-point rule was measured 30 per cent low in the
      ! Ly-alpha band), so this is what makes the two absorbers of one beam
      ! discretized alike.
      !
      ! WHAT IS AVERAGED.  Within a cell the H2 density and the continuum
      ! absorber densities are uniform: that is the rectangle rule the
      ! column integration itself uses (calc_column_dens_one).  Both the H2
      ! column N and the continuum depth tau therefore run LINEARLY across
      ! the cell from their star-ward face values to their inner face
      ! values, and the mean of the rate over the cell's radial extent is
      !
      !   <k> = int_0^1 k(N(s), tau(s)) ds ,
      !   N(s) = N_out + s dN ,   tau(s) = tau_out + s dtau .
      !
      ! HOW IT IS INTEGRATED, AND WHY NOT IN CLOSED FORM.  The continuum
      ! factor on its own would give the (1 - e^-dtau)/dtau of the water
      ! bands, because it IS exponential in the column.  The line term is
      ! not: sigma_diss comes from the level-resolved table, which is
      ! linear in (log N, log sigma) between its column knots, i.e. a
      ! piecewise power law of a slope that reaches 3.3 and steps by up to
      ! 0.8 from one knot interval to the next.  The product has no closed
      ! form, so the mean is taken by composite three-point
      ! Gauss-Legendre on segments cut GEOMETRICALLY in the column, at most
      ! seg_dex decades of column and seg_dtau of continuum depth wide.
      ! Geometric segments are what a power law needs: the rule is then
      ! integrating a nearly constant number of decades per segment
      ! wherever the cell sits on the transition, and the slope steps at
      ! the table's knots are resolved by several segments each.
      !
      ! MEASURED ACCURACY (2026-09-05, 1500 random cells covering the whole
      ! table, T = 700-3200 K, n_H = 3e11-3e14 cm^-3, column spans 0.01-10
      ! decades including cells whose star-ward column is zero, continuum
      ! depths to 3): the largest relative departure from a 4000-segment
      ! reference of the same integrand is 2.5e-4 at seg_dex = 0.05, and
      ! 9e-5 at 0.02.
      !
      ! The head segment.  A cell whose star-ward column is zero (the
      ! outermost cell, and any cell at the top of the H2 distribution) has
      ! no geometric starting point, so the segments start at
      ! col_head_ratio of the inner face column and the remainder below
      ! that is one further segment.  The table returns its edge value
      ! below N_H2 = 1e12 cm^-2, so with col_head_ratio = 1e-10 and any
      ! column the table carries, sigma_diss is constant over that head and
      ! the single segment integrates it exactly.
      double precision function column_cell_mean_of_a_cross_section(     &
                                sigma_at, F_LW, N_H2_out, N_H2_in, T,     &
                                n_H, tau_out, tau_in) result(k)
      procedure(h2_cross_section_of_state) :: sigma_at
      real*8, intent(in) :: F_LW, N_H2_out, N_H2_in, T, n_H
      real*8, intent(in) :: tau_out, tau_in
      ! Three-point Gauss-Legendre on [0,1]: nodes (1 -+ sqrt(3/5))/2 and
      ! 1/2, weights 5/18, 8/18, 5/18.  Exact for polynomials of degree 5.
      integer, parameter :: n_gl = 3
      real*8, parameter :: gl_s(n_gl) = (/ 0.1127016653792583d0,         &
                                           0.5d0,                        &
                                           0.8872983346207417d0 /)
      real*8, parameter :: gl_w(n_gl) = (/ 5.0d0/18.0d0, 8.0d0/18.0d0,   &
                                           5.0d0/18.0d0 /)
      real*8, parameter  :: seg_dex   = 0.05d0
      real*8, parameter  :: seg_dtau  = 0.5d0
      real*8, parameter  :: col_head_ratio = 1.0d-10
      integer, parameter :: nseg_max  = 256
      real*8  :: se(0:nseg_max+1)
      real*8  :: Nlo, Nhi, dN, dtau, tau_lo, Nstart, ratio, Nx
      real*8  :: s0, s1, ds, ss, acc
      integer :: nseg, nedge, i, g

      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif

      Nlo    = max(N_H2_out, 0.0d0)
      Nhi    = max(N_H2_in, Nlo)
      dN     = Nhi - Nlo
      tau_lo = max(tau_out, 0.0d0)
      dtau   = max(tau_in - tau_lo, 0.0d0)

      ! Segment count: one criterion per varying factor, and the finer wins.
      nseg   = 1
      Nstart = Nhi
      if (dN .gt. 0.0d0) then
         Nstart = max(Nlo, Nhi*col_head_ratio)
         nseg   = max(nseg, int(log10(Nhi/Nstart)/seg_dex) + 1)
      endif
      nseg = max(nseg, int(dtau/seg_dtau) + 1)
      nseg = min(nseg, nseg_max)

      ! Segment edges, as positions s across the cell.
      se(0) = 0.0d0
      nedge = 0
      if (dN .gt. 0.0d0) then
         if (Nstart .gt. Nlo) then
            nedge = 1
            se(1) = (Nstart - Nlo)/dN
         endif
         ratio = (Nhi/Nstart)**(1.0d0/dble(nseg))
         Nx    = Nstart
         do i = 1,nseg-1
            Nx = Nx*ratio
            se(nedge+i) = (Nx - Nlo)/dN
         enddo
         nedge = nedge + nseg
      else
         do i = 1,nseg-1
            se(i) = dble(i)/dble(nseg)
         enddo
         nedge = nseg
      endif
      se(nedge) = 1.0d0

      acc = 0.0d0
      do i = 1,nedge
         s0 = se(i-1)
         s1 = se(i)
         ds = s1 - s0
         if (ds .le. 0.0d0) cycle
         do g = 1,n_gl
            ss  = s0 + ds*gl_s(g)
            acc = acc + ds*gl_w(g)                                       &
                * band_rate_from_cross_section(sigma_at, F_LW,           &
                                          Nlo + ss*dN, T, n_H,           &
                                          tau_lo + ss*dtau)
         enddo
      enddo
      k = acc
      end function column_cell_mean_of_a_cross_section

      ! MEAN over one grid cell of the H2 photodissociation rate [s^-1].
      ! The cross section is the dissociations per incident band photon;
      ! this is the rate the H2 network destroys H2 at and the rate the
      ! 0.4 eV fragment heating is charged on.
      double precision function lyman_werner_dissociation_rate_cell_mean( &
                                F_LW, N_H2_out, N_H2_in, T, n_H,          &
                                tau_out, tau_in) result(k)
      real*8, intent(in) :: F_LW, N_H2_out, N_H2_in, T, n_H
      real*8, intent(in) :: tau_out, tau_in
      k = column_cell_mean_of_a_cross_section(                            &
              h2_lw_dissociation_cross_section, F_LW, N_H2_out, N_H2_in,  &
              T, n_H, tau_out, tau_in)
      end function lyman_werner_dissociation_rate_cell_mean

      ! MEAN over one grid cell of the rate at which the Lyman and Werner
      ! lines REMOVE BAND PHOTONS from the beam, per H2 molecule [s^-1].
      ! Same quadrature, same cell, sigma_pump in place of sigma_diss.
      !
      ! WHY IT IS NOT THE DISSOCIATION RATE DIVIDED BY p_eff.  p_eff varies
      ! with column, and a cell of a molecular base can span decades of it,
      ! so dividing the cell MEAN of sigma_diss by p_eff at one column is
      ! not the cell mean of sigma_diss/p_eff.  The photon ledger needs the
      ! second, because it is what telescopes to the column integral of
      ! sigma_pump, i.e. to the beam's own loss
      ! (h2_lw_band_photon_fraction_absorbed).  MEASURED: on
      ! cells spanning 2.5 H2 scale heights the one-column division stands
      ! 8 per cent above the beam loss and this form closes on it.
      double precision function                                           &
                lyman_werner_band_absorption_rate_cell_mean(              &
                                F_LW, N_H2_out, N_H2_in, T, n_H,          &
                                tau_out, tau_in) result(k)
      real*8, intent(in) :: F_LW, N_H2_out, N_H2_in, T, n_H
      real*8, intent(in) :: tau_out, tau_in
      k = column_cell_mean_of_a_cross_section(                            &
              h2_lw_pump_cross_section, F_LW, N_H2_out, N_H2_in,          &
              T, n_H, tau_out, tau_in)
      end function lyman_werner_band_absorption_rate_cell_mean

      ! End of module
      end module lyman_werner_photodissociation

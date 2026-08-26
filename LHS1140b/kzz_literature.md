# Eddy diffusion coefficient K_zz at ~1 microbar: what the literature supports

Prepared 2026-08-25 for the `He_Kzz` [cm^2/s] parameter of the binary H/He
element-diffusion operator in EXHALE (`src/modules/functions/
binary_element_diffusion.f90`, key `He_Kzz:` in `input.inp`).

The target condition is the LHS 1140 b wind base as it is actually set in the
run directory: `Base BC: pressure 1.0` in `LHS1140b/exhale/heh0p55_diff_ctrl/
input.inp` (1.0 dyn/cm^2 = 1.0 microbar = 1.0e-6 bar), `T_0 = 226 K` (the
equilibrium temperature echoed by `EXHALE_setup.out`), and the derived base
density `n_0 = 3.2049e13 cm^-3` reported by `run.log`
("Base BC pressure mode: derived n0"). That density is exactly p/kT, so it is
the total number density, not a hydrogen-only density.

Planet parameters used below: R_p = 1.73 R_earth, M_p = 5.60 M_earth,
T_eq = 226 K, giving g = 1837 cm/s^2 = 18.4 m/s^2 (computed here).

---

## 1. Conclusion

For a hydrogen/helium sub-Neptune atmosphere at 1 microbar, the published
estimates that can be extrapolated to that level cluster between about
**1e9 and 1e11 cm^2/s**, and the one paper that states a value at exactly
1e-6 bar for an H/He escape model (Taylor et al. 2025) adopts **1e9 cm^2/s**.
A physically derived homopause scaling (Arfaux & Lavvas 2023) applied to
LHS 1140 b's gravity and temperature gives a much lower **3e7-1e8 cm^2/s**,
and the one temperate sub-Neptune study in the list (Blain et al. 2021,
K2-18 b) adopts a nominal **1e6 cm^2/s** with an explored range of
1e5-1e10 cm^2/s. Every one of these numbers comes from an atmosphere warmer
and more strongly irradiated than LHS 1140 b's, and the scalings all point
the same direction, so the LHS 1140 b value likely sits at or below the low
end. A defensible working range appears to be **1e7-1e9 cm^2/s**, with
1e9 cm^2/s as an upper-end choice inherited from hot-Jupiter practice and
1e7 cm^2/s as the value the homopause scaling actually predicts for this
planet. The spread is two orders of magnitude and none of it is measured for
this object; `He_Kzz` should be treated as a scanned parameter, not a
constant taken from a table.

A useful bracket: the molecular H-He binary diffusion coefficient at the
EXHALE base, using the Banks & Kockarts (1973) hard-sphere form that
`binary_element_diffusion.f90` itself implements
(D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2)/n), is

    D(H,He) = 1.52e18 * sqrt(1 + 1/4) * sqrt(226) / 3.2049e13
            = 8.0e5 cm^2/s   (computed here)

so **any K_zz above ~1e6 cm^2/s places the homopause above the wind base**,
i.e. the base is well mixed and the eddy term, not molecular diffusion,
controls the composition there. Every literature value surveyed is above
that threshold. The choice of `He_Kzz` within 1e7-1e9 therefore sets *how far
above the base* the homopause sits, not whether the base separates.

---

## 2. Table of literature estimates

`K_zz` values are given in cm^2/s. Where a source quotes m^2/s the original
figure is shown with the conversion in brackets (1 m^2/s = 1e4 cm^2/s).
Numbers marked *(computed here)* are my evaluation of the source's own
formula, not a number printed in the source.

| reference | object / atmosphere type | K_zz [cm^2/s] | pressure at which quoted | how obtained | verified? |
|---|---|---|---|---|---|
| Charnay, Meadows & Leconte (2015) | GJ 1214 b, warm sub-Neptune, H2/He at 1x, 10x, 100x solar and pure H2O | K_zz = K_zz0 * P_bar^(-0.4) m^2/s with K_zz0 = 7e2, 2.8e3, 3e3, 3e2 m^2/s [7e6, 2.8e7, 3e7, 3e6 cm^2/s] for 1x, 10x, 100x solar and pure water | fit tuned for P < 10 mbar; GCM spans 80 bar to 3 Pa (3e-5 bar) | 3-D GCM with passive settling tracers, effective 1-D K_zz from the flux/gradient ratio | yes (published PDF, ApJ 813:15) |
| Parmentier, Showman & Lian (2013) | HD 209458 b, hot Jupiter | K_zz = 5e4 / sqrt(P_bar) m^2/s = 5e8 / sqrt(P_bar) cm^2/s | stated valid "from ~1 bar to a few microbar" (abstract); conclusions say "between ~1 bar and ~1 microbar" | 3-D GCM tracer fields fitted by a 1-D diffusion model | partial - numerals read from arXiv:1301.4522; the sentences "which is valid over a pressure range from" and "can be used in 1D models of HD 209458b" are confirmed present in the ADS full-text index of the published A&A article, but the published PDF was not retrievable (aanda.org served a captcha page) |
| Moses et al. (2011) | HD 189733 b, HD 209458 b, hot Jupiters | deep value K_zz = 1e10; low-pressure plateau read off Fig. 1 at roughly 1e11 (HD 189733 b) and 4e11 (HD 209458 b) | deep value in the adiabatic region; plateau applies for P below ~1e-3 bar, where "K_zz's are assumed to be constant at altitudes above the top level of the GCMs" | mixing-length / free convection in the deep atmosphere; K_zz = w(z)L(z) with w from Showman et al. (2009) GCM and L = H above it | partial - Fig. 1 and the caption read from arXiv:1102.0063; the caption sentence "is an estimate of the eddy diffusion coefficient in the adiabatic region" is confirmed present in the ADS full-text index of the published ApJ article; the plateau values are my reading of the figure, not printed numbers |
| Koskinen et al. (2013b), Icarus 226, 1695 | HD 209458 b, hot Jupiter, lower atmosphere / cold trap | K_zz >~ 1e5 m^2/s [1e9] sufficient to prevent cloud settling; quoting Moses et al. (2011), "the high pressure value of Kzz = 10^6 m^2 s^-1" [1e10] "which implies that Kzz >~ 10^7 m^2 s^-1" [1e11] "at p <~ 10 mbar" | 1-10 mbar and below | argument from particle settling velocity vs. mixing; citation of Moses et al. (2011) | yes (published PDF, `references/Koskinen_2013Icarus_226_1695.pdf`) |
| Koskinen et al. (2013a), Icarus 226, 1678 | HD 209458 b, thermosphere | no value quoted; argues eddy mixing "is unlikely to mix the atmosphere up to 3Rp and beyond" and that "there is considerable uncertainty over the values of Kzz" | thermosphere, above ~1 microbar | discussion; the paper attributes the observed lack of diffusive separation to rapid escape rather than to eddy mixing | yes (published PDF, `references/Koskinen_2013Icarus_226_1678.pdf`) |
| Koskinen et al. (2022), ApJ 929, 52 | hot Jupiter / sub-Neptune escape model | K_zz appears only as a symbol in the diffusion equation (Appendix B, Lambda_s = K_zz/D_s); no numerical value appears in the extractable text | - | formulation only | yes (published PDF, `references/Koskinen_2022_ApJ_929_52.pdf`); absence of a number confirmed by grep over the full extracted text |
| Taylor et al. (2025), ApJ 989, 68 | HD 209458 b, hot Jupiter, H/He escape model with multispecies diffusion | best-fit model K_zz = 1e5 m^2/s [1e9]; sensitivity runs at 1e6 and 1e7 m^2/s [1e10, 1e11] plus a no-diffusion case | Model B's lower boundary is a pressure of 1e-6 bar, i.e. exactly the microbar level of interest; K_zz is a single constant for the upper-atmosphere model | adopted, "based on a recent assessment of eddy diffusion in exoplanet atmospheres (A. Arfaux & P. Lavvas 2023)" | yes (published PDF, `references/Taylor_2025_ApJ_989_68.pdf`) |
| Arfaux & Lavvas (2023), MNRAS 522, 2525 | giant planets / hot Jupiters; homopause value K_top | K_zz = K_zz^J (v/v_J)(g_J/g)(T/T_J) with K_zz^J = 1e7 cm^2/s (Jupiter homopause, after Koskinen et al. 2010); v/v_J = 5 retained for HD 209458 b. Applied to LHS 1140 b: 2.8e7 (v/v_J = 1) to 1.4e8 (v/v_J = 5) *(computed here)* | at and above the homopause, where K_zz is taken constant | physically derived parameterization: deep convective branch, gravity-wave middle branch, constant homopause plateau | partial - read from arXiv:2304.06314; the phrase "the homopause eddy coefficient of Jupiter" is confirmed present in the ADS full-text index of the published MNRAS article; the publisher PDF was not retrievable |
| Blain, Charnay & Bezard (2021), A&A 646, A15 | K2-18 b, temperate sub-Neptune, H2/He | "typical Kzz values for K2-18b range from 10^6 to 10^9 cm^2 s^-1, with the highest values found in the convective layers"; nominal value 1e6; model grid widened to 1e5-1e10 | height-independent constant across the modeled column; the model is the observable/deep atmosphere, not a thermosphere | mixing-length estimates taken from the Exo-REM model (Charnay et al. 2018) | partial - numerals read from arXiv:2011.10459; the phrases "values for K2-18b range from" and "with the highest values found in the convective layers" and "we enlarged this range to" are confirmed present in the ADS full-text index of the published A&A article |
| Zhang & Showman (2018a, 2018b), ApJ 866, 1 and 866, 2 | general theory, fast rotators and tidally locked planets | no single value; the result is that K_zz depends on tracer chemical lifetime, circulation strength and horizontal eddy mixing, that "different chemical species in a single atmosphere should in principle have different eddy diffusion profiles", and that non-diffusive transport can make an inferred K_zz negative | all | analytic theory validated against 2-D and 3-D tracer transport simulations | partial - read from arXiv:1803.09149 and arXiv:1808.05365; publisher PDFs not retrievable |
| Cherubim et al. (2026), Science | LHS 1140 b, helium-dominated outflow | none - no eddy diffusion, molecular diffusion, or homopause appears anywhere in the paper | - | not modeled; composition is a free parameter of a homogeneous Parker wind (see section 5) | yes (published PDF; `grep -ci "eddy\|diffus\|homopause\|Kzz"` returns 0 for both the main text and the supplement) |

Transcription note: quoted passages are copied from `pdftotext` output of the
PDFs listed in section 6. Superscripts, subscripts and the relational glyphs
">~" and "<~" (for the "greater/less than or approximately" signs) do not
survive text extraction and are rendered in ASCII above; the wording is
otherwise unaltered.

### 2.1. A unit discrepancy in Taylor et al. (2025)

Taylor et al. (2025), section 3.5.2, motivates its K_zz choices with:

> J. I. Moses et al. (2011) suggested values between 10^8 and 10^9 m^2 s^-1
> for hot Jupiter atmospheres, while V. Parmentier et al. (2013) proposed
> lower values around 10^6-10^7 m^2 s^-1.

(verbatim from the published ApJ 989:68 PDF; superscripts flattened.)

Checked against the originals:

- **Moses et al. (2011)** plots K_zz on an axis labeled
  "Eddy Diffusion Coefficient (cm^2 s^-1)" spanning 1e8 to 1e12, i.e.
  1e4 to 1e8 m^2/s. Taylor's "10^8 to 10^9 m^2 s^-1" is 1e12-1e13 cm^2/s,
  about four orders of magnitude above anything in Moses's figure. Taylor's
  numerals match Moses's plotted range only if read as cm^2/s. This appears
  to be a cm^2 / m^2 slip in Taylor et al. (2025), though I cannot rule out
  that the authors intended some other quantity.
- **Parmentier et al. (2013)** gives K_zz = 5e4 / sqrt(P_bar) m^2/s, which
  runs from 5e4 m^2/s at 1 bar to 5e7 m^2/s at 1 microbar. Taylor's
  "around 10^6-10^7 m^2 s^-1" falls inside that range (it corresponds to
  roughly 1e-3 to 1e-5 bar), so this second figure is **not** obviously a
  unit slip - it appears to be a mid-range readout of Parmentier's formula
  in the original units.

Net effect: the ordering Taylor states ("Parmentier proposed lower values")
is correct in spirit - Parmentier explicitly finds his GCM-derived K_zz to be
"two orders of magnitude smaller than what is obtained when multiplying the
vertical scale height by the root mean square of the vertical velocity", the
w_rms * H estimate that Moses et al. (2011) used - but the Moses figure as
printed cannot be right in m^2/s.

Note also that Taylor et al. cite Arfaux & Lavvas as **2022** at the point
where Model B's K_zz = 1e5 m^2/s is introduced and as **2023** in section
3.5.2 where the same choice is justified. The 2023 paper (MNRAS 522, 2525)
is the eddy parameterization; the 2022 paper (MNRAS 515, 4753) is the
haziness study.

---

## 3. Extrapolation to 1 microbar

Power-law fits, written in the source's own symbols and units, then evaluated
at P = 1e-6 bar. All evaluations in this section are *computed here*.

**Charnay, Meadows & Leconte (2015)**, GJ 1214 b, from the abstract and
section 4.4:

    K_zz = K_zz0 * P_bar^(-0.4)   [m^2/s, P in bar]
    K_zz0 = 7e2 (1x solar), 2.8e3 (10x), 3e3 (100x), 3e2 (pure H2O)

At P = 1e-6 bar, P^(-0.4) = 251.2:

| composition | K_zz (1 microbar) [m^2/s] | [cm^2/s] |
|---|---|---|
| 1x solar | 1.76e5 | **1.8e9** |
| 10x solar | 7.03e5 | **7.0e9** |
| 100x solar | 7.54e5 | **7.5e9** |
| pure H2O | 7.54e4 | **7.5e8** |

**This is an extrapolation outside the fit's stated range in two ways.**
The paper says "The parameters in our fits were chosen to primarily match
Kzz and the mean tracer abundance (Figure 11) for pressures lower than
10 mbar", so the *low*-pressure side is where the fit is meant to work; but
the GCM itself has "the first level at 80 bars and the top level at 3 Pa",
i.e. 3e-5 bar. Evaluating at 1e-6 bar is about 1.5 decades above the model
top, where no GCM information exists. The paper offers a physical reason the
slope should persist - it argues that planetary waves "generally grow as
P^-1/2" and expects "an eddy diffusion coefficient with an exponent between
-1/3 and -1/2 for the pressure dependence" - but that is an expectation, not
a validated extrapolation.

**Parmentier, Showman & Lian (2013)**, HD 209458 b, Eq. (22) and conclusions:

    K_zz = 5e4 / sqrt(P_bar)  m^2/s  =  5e8 / sqrt(P_bar)  cm^2/s

At P = 1e-6 bar: K_zz = 5e7 m^2/s = **5e11 cm^2/s**. Unlike the Charnay fit,
this one is *inside* its stated validity range - the paper says the
parameterization holds "between ~1 bar and ~1 microbar".

**Moses et al. (2011)**, HD 189733 b / HD 209458 b: no formula; K_zz is
constant above the GCM top, at roughly 1e11 and 4e11 cm^2/s respectively
(read from Fig. 1). No extrapolation is needed - the plateau already covers
1 microbar by construction. Koskinen et al. (2013b) summarizes the same
profiles as "Kzz >~ 10^7 m^2 s^-1 at p <~ 10 mbar", i.e. >~ 1e11 cm^2/s,
which agrees with my reading of the figure.

**Arfaux & Lavvas (2023)** homopause scaling, applied to LHS 1140 b
(this is the only entry here derived for *this planet* rather than
extrapolated from another one):

    K_top = K_zz^J * (v/v_J) * (g_J/g) * (T/T_J),  K_zz^J = 1e7 cm^2/s

with g_J = 24.79 m/s^2, T_J = 110 K (Jupiter's equilibrium temperature),
g = 18.4 m/s^2 and T = 226 K for LHS 1140 b:

| assumed v/v_J | K_top [cm^2/s] |
|---|---|
| 1 (Jupiter-like turbulence) | **2.8e7** |
| 5 (the value the paper retains for HD 209458 b) | **1.4e8** |

Sanity check of the scaling: the same formula with HD 209458 b's
g = 9.4 m/s^2, T = 1450 K and v/v_J = 5 gives 1.7e9 cm^2/s, within a factor
of 2 of the 1e9 cm^2/s that Taylor et al. (2025) actually adopt on the
strength of this parameterization. The scaling therefore reproduces the
hot-Jupiter practice it is supposed to, which makes its much lower LHS 1140 b
prediction worth taking seriously.

**Taylor et al. (2025)**: 1e9 cm^2/s, quoted at a model lower boundary of
1e-6 bar. No extrapolation at all - this is a directly comparable number, for
a hot Jupiter.

**Blain et al. (2021)**: 1e6 cm^2/s nominal, height-independent, for a
temperate sub-Neptune. Their model is not a thermosphere model and the
constant is not attached to a microbar level, so this is not directly
transferable; it is included because it is the only entry at a comparable
temperature.

### Summary of the extrapolation

| source | K_zz at 1 microbar [cm^2/s] | in range? |
|---|---|---|
| Charnay+2015 (GJ 1214 b, 1x-100x solar) | 1.8e9 - 7.5e9 | no - ~1.5 decades above the GCM top |
| Charnay+2015 (pure H2O) | 7.5e8 | no - same |
| Parmentier+2013 (HD 209458 b) | 5e11 | yes |
| Moses+2011 (hot Jupiters) | ~1e11 - 4e11 | yes (constant plateau by assumption) |
| Taylor+2025 (HD 209458 b) | 1e9 | yes (quoted at 1e-6 bar) |
| Arfaux & Lavvas 2023 scaling, applied to LHS 1140 b | 2.8e7 - 1.4e8 | yes (homopause plateau) |
| Blain+2021 (K2-18 b) | 1e6 nominal, 1e5-1e10 explored | not pressure-resolved |

Taken together: **1e9 to 1e11 cm^2/s** for hot, H2-rich atmospheres at
1 microbar, dropping to **1e7-1e8 cm^2/s** when the one available physical
scaling is evaluated with LHS 1140 b's own gravity and temperature.

---

## 4. Temperature-scale caveat

Every atmosphere in the table above is warmer than LHS 1140 b's base by a
large factor: HD 209458 b's thermosphere base sits near 1270 K
(Taylor et al. 2025, Model B) and GJ 1214 b is a warm sub-Neptune orbiting
at 0.014 AU (Charnay et al. 2015), against T_0 = 226 K here; only K2-18 b
in the table is at a comparable temperature. What follows is an
argument from the scalings the papers themselves write down; **it is an
expectation, not a result, and none of it has been checked against a model
of this planet.**

**Free-convection / mixing-length branch.** Charnay et al. (2015) Eq. (16),
the Gierasch & Conrath (1985) form used in the Ackerman & Marley (2001) cloud
model, is

    K_zz = (H/3) (L/H)^(4/3) (R F_c / (mu rho_a c_p))^(1/3)

so K_zz scales as F_c^(1/3) and linearly with the scale height H. Both push
down for LHS 1140 b. The convective flux F_c in this form is the heat flux
the convection actually carries -- the *internal* heat flux, since absorbed
starlight is deposited above the convective region and stabilizes the
stratification rather than driving it (the day/night insolation contrast
belongs to the wave-driving branch below, not here; this paragraph
originally compared the two planets by insolation, corrected 2026-08-25 --
see docs/eddy_diffusion_kzz.tex section 3.1). A 5.6 M_Earth planet of
several Gyr has a far smaller internal flux than an inflated hot Jupiter,
which is the sense in which this branch pushes down here. The scale height is
H = kT/(mu m_H g) = 25 km for a helium-dominated base (mu = 4, computed here)
against about 860 km at HD 209458 b's thermosphere base (T = 1270 K,
mu = 1.3, g = 9.4 m/s^2; computed here). Charnay et al. also warn that this
formula "strongly overestimates the convective heat flux for irradiated
planets" when F_c is set to sigma T_eff^4, so the naive version of this
scaling is an upper bound in any case.

**Breaking-gravity-wave branch.** Charnay et al. note that wave amplitude
"grows as P^-1/2 until they break", which fixes the *slope* of K_zz(P) but
not its magnitude. The magnitude is set by the wave source. On tidally
locked hot planets that source is the day/night thermal contrast, which is
what Parmentier et al. (2013) identify as the driver: "vertical mixing
results not from small-scale convection but from the large-scale circulation
driven by the day-night heating contrast". LHS 1140 b's day/night contrast in
absolute terms is far smaller, so the same mechanism should be weaker,
although the P^-0.4 to P^-1/2 slope would likely survive.

**The homopause scaling makes the same prediction.** Arfaux & Lavvas (2023)
carry the temperature dependence explicitly, K_top proportional to
v (T/g). Cooler means smaller directly through T, and the characteristic
turbulent velocity v is itself expected to fall with weaker forcing. That is
why the scaling lands at 3e7-1e8 cm^2/s for LHS 1140 b rather than the
~1e9 cm^2/s adopted for HD 209458 b.

**Direction of the bias.** All three arguments push the same way: the
hot-Jupiter and GJ 1214 b numbers **likely overestimate** K_zz for
LHS 1140 b, plausibly by one to two orders of magnitude. Nothing in the
surveyed literature argues the other way.

**One counter-consideration, also speculative.** LHS 1140 b's base is
helium-dominated (Cherubim et al. 2026 retrieve H:He < 1e-3), so mu ~ 4
rather than ~2.3. At fixed g and T that halves the scale height, which
reduces K_zz through the H prefactor in the mixing-length form but also
compresses the whole atmosphere, so the *number of scale heights* between the
base and the homopause is not simply reduced. Whether the compositional
gradient sharpens or flattens is not something this survey can settle.

---

## 5. Cherubim et al. (2026): how the LHS 1140 b modeling handles mixing

**It does not.** A case-insensitive grep for `eddy`, `diffus`, `homopause`
and `Kzz` over the extracted text of both `Cherubim_2026Science.pdf` and
`Cherubim_2026Science_Supplement.pdf` returns zero hits. There is no eddy
term, no molecular diffusion, no homopause, and no discussion of
compositional stratification anywhere in the paper. The composition is a
fitted constant.

What the paper does say, verbatim from the published PDF (Science first
release, 16 July 2026; page numbers not final in that release):

Main text, on the model:

> To determine the physical properties of the atmospheric outflow, we
> modeled the observed spectra of LHS 1140b using the p-winds code (28).
> This code models an escaping atmosphere as a one-dimensional, isothermal
> outflow (29) then forward models a predicted transmission spectrum using a
> radiative transfer model.

Supplementary materials, "Retrieving atmospheric properties":

> The p-winds model assumes a spherical, homogeneous outflow.

and, on how composition enters:

> We used the p-winds models to constrain the mass-loss rate, outflow
> temperature, and H:He ratio.

The H:He ratio is thus a single number applied to the whole outflow. The
paper's central compositional claim - the retrieved H:He ratio below about
1e-3, stated in the main-text section "Atmospheric composition" (the
inequality sign does not survive text extraction, so the numeral is
reproduced here rather than quoted) - is a property of a uniform wind, not of
a stratified one. The supplement records one related numerical
test, but it concerns the mean molecular weight rather than mixing:

> We tested the assumption in the p-winds code that the electrons from
> helium ionization do not contribute to the atmospheric mean molecular
> weight by fixing the mean molecular weight at 4 amu, the theoretical
> maximum for a pure, neutral helium wind. The transmission spectra were
> unchanged by this assumption.

(Note: the extracted text renders several ligatures and minus signs
imperfectly - "e!ects" for "effects", an arrow glyph for the minus sign in
exponents. The quotations above were chosen to avoid affected characters and
were checked against the rendered PDF text.)

Consequence for EXHALE: the Cherubim et al. H:He is a *well-mixed* value with
no altitude dependence, so it cannot be compared directly with an EXHALE
profile in which He/H varies with radius. If EXHALE's `He_diffusion` produces
a He/H that changes appreciably between the base and the He 10830 formation
region, the two models are answering different questions, and the
Cherubim et al. number is best read as an effective value weighted toward
wherever the 10830 line forms.

---

## 6. References

- Arfaux, A. & Lavvas, P. 2023, MNRAS, 522, 2525, "A physically derived eddy
  parametrization for giant planet atmospheres with application to
  exoplanets" - bibcode `2023MNRAS.522.2525A`;
  `references/Arfaux_2023_MNRAS_522_2525_arXiv.pdf` (arXiv:2304.06314).
- Blain, D., Charnay, B. & Bezard, B. 2021, A&A, 646, A15, "1D atmospheric
  study of the temperate sub-Neptune K2-18b" - bibcode
  `2021A&A...646A..15B`; `references/Blain_2021_A+A_646_A15_arXiv.pdf`
  (arXiv:2011.10459).
- Charnay, B., Meadows, V. & Leconte, J. 2015, ApJ, 813, 15, "3D Modeling of
  GJ1214b's Atmosphere: Vertical Mixing Driven by an Anti-Hadley
  Circulation" - bibcode `2015ApJ...813...15C`;
  `references/Charnay_2015_ApJ_813_15.pdf` (publisher PDF).
- Cherubim, C. et al. 2026, Science, "Helium escaping from the atmosphere of
  a nearby rocky exoplanet orbiting in a habitable zone" -
  `references/Cherubim_2026Science.pdf` and
  `references/Cherubim_2026Science_Supplement.pdf` (publisher PDF, first
  release 16 July 2026).
- Koskinen, T. T. et al. 2013a, Icarus, 226, 1678, "The escape of heavy atoms
  from the ionosphere of HD209458b. I." - bibcode `2013Icar..226.1678K`;
  `references/Koskinen_2013Icarus_226_1678.pdf`.
- Koskinen, T. T. et al. 2013b, Icarus, 226, 1695, "The escape of heavy atoms
  from the ionosphere of HD209458b. II." - bibcode `2013Icar..226.1695K`;
  `references/Koskinen_2013Icarus_226_1695.pdf`.
- Koskinen, T. T. et al. 2022, ApJ, 929, 52 - bibcode `2022ApJ...929...52K`;
  `references/Koskinen_2022_ApJ_929_52.pdf`.
- Moses, J. I. et al. 2011, ApJ, 737, 15, "Disequilibrium Carbon, Oxygen, and
  Nitrogen Chemistry in the Atmospheres of HD 189733b and HD 209458b" -
  bibcode `2011ApJ...737...15M`;
  `references/Moses_2011_ApJ_737_15_arXiv.pdf` (arXiv:1102.0063).
- Parmentier, V., Showman, A. P. & Lian, Y. 2013, A&A, 558, A91, "3D mixing
  in hot Jupiters atmospheres. I. Application to the day/night cold trap" -
  bibcode `2013A&A...558A..91P`;
  `references/Parmentier_2013_A+A_558_A91_arXiv.pdf` (arXiv:1301.4522).
- Taylor, J. et al. 2025, ApJ, 989, 68 - bibcode `2025ApJ...989...68T`;
  `references/Taylor_2025_ApJ_989_68.pdf` (publisher PDF).
- Zhang, X. & Showman, A. P. 2018a, ApJ, 866, 1, "Global-mean Vertical Tracer
  Mixing in Planetary Atmospheres. I." - bibcode `2018ApJ...866....1Z`;
  `references/Zhang_2018_ApJ_866_1_arXiv.pdf` (arXiv:1803.09149).
- Zhang, X. & Showman, A. P. 2018b, ApJ, 866, 2, "Global-mean Vertical Tracer
  Mixing in Planetary Atmospheres. II. Tidally Locked Planets" - bibcode
  `2018ApJ...866....2Z`; `references/Zhang_2018_ApJ_866_2_arXiv.pdf`
  (arXiv:1808.05365).

Also consulted and found not to carry a usable number: Tsai et al. 2021
(ApJ 923, 264; `2021ApJ...923..264T`) uses a uniform K_zz = 1e8 cm^2/s for
its Jupiter model and takes exoplanet K_zz profiles from other work;
Wogan et al. 2025 (PSJ 6, 256; `references/Wogan_2025_Planet._Sci._J._6_256.
pdf`) defines K_zz in its transport equations but adopts profiles from the
cited literature rather than deriving values.

### Access notes

Publisher PDFs were retrieved for Charnay et al. (2015), Taylor et al.
(2025), the two Koskinen et al. (2013) papers, Koskinen et al. (2022) and
Cherubim et al. (2026). IOP (ApJ) and OUP (MNRAS) returned HTML rather than
PDF through the ADS link gateway, and aanda.org (A&A) served a captcha page;
for those the arXiv PDFs were read instead, and key sentences were checked
for presence in the ADS full-text index of the published article using
`body:"..."` phrase queries. Where a numeral could not be confirmed against
the published text, the table says so explicitly.

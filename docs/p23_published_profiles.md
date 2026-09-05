# P23 source material: the published hot-Uranus and super-Earth profiles, read off the figures

**What this file is.** Gate P23 compares the three-layer structure of our
converged hot-Uranus solution against two published models. Up to now that
comparison (`supersonic_molecular_base.md` sections 13.5 and 13.6) rested only
on *sentences* from the two papers; no number was taken from their figures.
This file supplies the missing quantitative side: the model inputs, the method
each paper uses for H2, the relevant passages quoted from the published text,
and the radial profiles of `f(H2)`, `n(H3+)`, `T` and `n_e` digitized from the
published figures, tabulated beside our own converged state on a common radial
grid.

**What this file deliberately does not do.** It draws no conclusion from the
comparison. Item P43 puts the H2 transport term into the code first, and the
judgment is made after that, against these numbers. Everything below is
source material.

Both papers were read from the publisher PDFs in `references/`:

* Koskinen, T. T., Lavvas, P., Huang, C., Bergsten, G., Fernandes, R. B., &
  Yelle, R. V. 2022, ApJ, 929, 52, "Mass Loss by Atmospheric Escape from
  Extremely Close-in Planets" (`references/Koskinen_2022_ApJ_929_52.pdf`).
  The hot-Uranus run used here is **Model A at a = 0.05 au**.
* Frelikh, R., & Murray-Clay, R. A. 2026, ApJ, 996, 96, "Efficiency of
  Hydrodynamic Atmospheric Escape in Hot Jupiters and Super-Earths"
  (`references/Frelikh_2026_ApJ_996_96.pdf`). The run used here is the
  **fiducial super-Earth of section 4**, parameters in their Table 3.

Our side is the converged `He/H = 0.0793` rung of
`supersonic_molecular_base.md` section 13.2 (run directory `heh0p0793c_r1`,
`Molecular chemistry: True`, metals off, `info = 0`,
`log10 Mdot = 10.42` g/s as printed by the run itself).

---

## 1. Model inputs, side by side

Section numbers and table numbers below were checked to exist in the published
versions.

| | this work, He/H = 0.0793 | Koskinen+2022 Model A | Frelikh & Murray-Clay 2026 super-Earth |
|---|---|---|---|
| planet | hot Uranus | tidally locked Uranus-like | fiducial super-Earth |
| source of inputs | `input.inp`, `base.inp`, `EXHALE_resolved.out` | section 3 text; Model A defined in section 3.2 | Table 3 (plus Table 2 for what Table 3 does not repeat) |
| mass | 0.0457 M_J = 8.674e28 g = 14.52 M_E | 86.813e24 kg = 8.6813e28 g = 14.54 M_E (Uranus) | 7.2 M_E |
| reference radius R_p | 0.49 R_J = 3.4256e9 cm (code constant R_J = 6.9911e9 cm) | 25,559 km = 2.5559e9 cm, the **polar 1 bar** radius of Uranus | 1.9 R_E = 1.2108e9 cm |
| radius of the model's lower boundary | 1.0 R_p = 3.4256e9 cm | 1.34 R_p = 3.4249e9 cm | 1.0 R_p = 1.2108e9 cm |
| gravity at the lower boundary | 493 cm s^-2 | 494 cm s^-2 | 1958 cm s^-2 |
| orbital distance | 0.048 au | 0.05 au | 0.05 au |
| host star | 1.2 M_sun | Sun-like (M_s = M_sun, L_s = L_sun) | 1 M_sun (Table 2) |
| base temperature T0 | 1140 K | 1140 K (stated with the lower boundary) | 1000 K |
| base pressure p0 | 8.99e-6 bar (8.99 dyn cm^-2) | 1e-6 bar (Appendix B: p1 = 1 microbar) | 0.47 dyn cm^-2 = 4.7e-7 bar |
| XUV / EUV flux at the planet | 10^3.193 = 1.56e3 erg cm^-2 s^-1 (13.6 eV - 1.24 keV) | F_XUV = 1.6 W m^-2 = 1.6e3 erg cm^-2 s^-1 (L_XUV = 1.12e21 W) | 2.6e3 erg cm^-2 s^-1 at 0.05 au (top axis label of their Figure 23) |
| spectrum | power law, index -1 | mean solar spectrum of Koskinen et al. (2013a, 2013b) | Richards et al. (1994) solar EUV, 37 bins, scaled by (au/a)^2; range 12-248 eV (Table 2) |
| He/H by number | 0.0793 | 7.93e-2 (Lodders 2003), used in their Equations (11)-(13) | not stated in the text; the base of their Figure 21 gives n(He)/n(H nuclei) ~ 0.085 (digitized, see section 5.3) |
| elements | H, He only (metals off) | H, He only ("we do not include elements heavier than H and He") | H, He only |
| species carried | H I, H II, He I, He II, He III, He(2^3S), H2, H2+, H3+, HeH+ | H2, H, He, H2+, H3+, H+, He+, HeH+ | H2, H2+, H, H+, He, He+, H3+, HeH+ |
| reaction list | `System_HeH_mol.f90` | Table 1 (P1-P5, R1-R23) | Table 1 (k1-k22) |
| H2 at the base | H nuclei in H2 = 0.925 (`q_H2_base = 0.75`, a representative photochemical value, not a computed one); volume mixing ratio n(H2)/n_tot = 0.861 | q0(H2) = 0.84, q0(H) = 0.026 | not stated; digitized n(H2)/n_tot at the first readable columns is 0.82-0.85 |
| tidal / Roche gravity | Roche potential, truncated at the Hill/L1 radius | Model A uses spherical gravity; Model B replaces it by the substellar Roche gravity | tidal term included in their Equation (7) |
| mass-loss rate | log10 Mdot = 10.42 (g s^-1), i.e. 2.1e9 g s^-1 sr^-1 | F_M = rho v r^2 = 1.5e6 kg s^-1 sr^-1 = 1.5e9 g s^-1 sr^-1; Mdot = 1.9e7 kg s^-1 = 1.9e10 g s^-1 | 7.4e8 g s^-1 sr^-1 at 0.05 au (digitized from their Figure 23), i.e. 9.3e9 g s^-1 over 4 pi |
| outer boundary of the domain | 4.73 R_p | 9.7 R_p = 7.2 x the lower boundary | 5.0e9 cm = 4.13 R_p |

Two notes on the radius normalization, because the three models do not
normalize alike:

1. **Koskinen's R_p is the 1 bar polar radius of Uranus and their model starts
   at 1.34 R_p.** Our R_p *is* the lower boundary of the escape model. The two
   lower boundaries are almost the same absolute radius: 3.4256e9 cm here
   against 3.4249e9 cm there, a ratio of 1.0002. Every joint table and the
   figure below therefore use `r/r_base`, r_base being each model's own lower
   boundary; for Koskinen that is `(r/R_p)/1.34`.
2. Section 13.5 of `supersonic_molecular_base.md` quotes our radius as
   `0.49 R_J = 5.49 R_E` and the base gravity as 473 cm s^-2. Those follow from
   R_J = 7.1492e9 cm (the equatorial radius); the code's constant is
   R_J = 6.9911e9 cm (`parameters.f90` line 454), which gives 5.38 R_E and
   493 cm s^-2. The values in the table above are the ones the run actually
   used. (That section is left as it stands; this is only a note that the two
   numbers differ and why.)

---

## 2. How each model treats H2: transport, or a local balance

This is the axis P43 is about, so it is worth separating what each code solves.

### 2.1 Koskinen et al. (2022): species continuity with advection and diffusion

Their Appendix B states the structure of the model:

> "As indicated above, the model solves separate continuity equations for the
> different neutral and ion species coupled to common momentum and energy
> equations for the bulk atmosphere. The use of the diffusion approximation is
> therefore necessary to solve for the individual species velocity
> perturbation ws."

Their Equation (B1) is a time-dependent continuity equation for each species
carrying the divergence of the species mass flux, and Equations (B9)-(B11) give
the diffusion velocities `ws` (molecular, ion-neutral, electron-ion, plus eddy
diffusion through `Lambda_s = Kzz/Ds`). Transport is therefore present twice:
in the bulk advection and in the species diffusion velocity.

They also show the H2 budget explicitly. Their Figure 16 is captioned
"Production, loss, and transport terms (absolute values) in the continuity
equation for H2 in Model A at 0.05 au", and the accompanying Appendix B text
reads:

> "We find a similar result for H2, as illustrated by Figure 16, which shows
> the continuity equation terms for H2. Chemical loss of H2 (dissociation and
> dissociative ionization reactions) dominates over production but outflow
> brings fresh H2 up from deeper layers of the atmosphere and primarily
> balances the total loss rate."

and, in section 3.2.1:

> "The prevalence of neutral H and H2 at high radii in the hot Uranus model is
> due to a combination of a relatively low temperature and high escape flux
> that replenishes H2 and H at high altitudes where they would otherwise be
> predominantly dissociated and ionized (see Appendix B)."

### 2.2 Frelikh & Murray-Clay (2026): species continuity with advection, diffusion off

Their Equation (6) in section 2.1 is the species continuity equation
`d(rho_s)/dt + (1/r^2) d(rho_s u r^2)/dr = S_s`, i.e. advection by the bulk
velocity plus a chemical source. Immediately after it they write:

> "We note that, while diffusion is implemented in our code, we do not expect
> it to significantly affect the purely hydrogen and helium models presented
> here, so we turn it off and leave further discussion for future work."

Section 2.9 confirms the coupling: "The evolution of the concentration of the
seven atomic/molecular species (+electrons) is given by Equation (6), where the
22 reactions of our chemical network are encapsulated by the source terms,
Ss."; the stiff network is integrated with a semi-implicit method (their
Appendix B), not solved as an equilibrium.

They also plot the transport term for H2 directly. In section 3.1.2:

> "As an aside, before we focus on the photochemistry, we derive the advective
> term, ntotu ∂fs/∂r, that transports a species s with the volumetric number
> density ns and number fraction fs = ns/ntot, where ntot = Σsns is the total
> number density of all species."

and their Figure 7 is captioned "Selected reaction rates involving H2. The
advective term is plotted as derived in Equation (37)."

### 2.3 EXHALE, as the run tabulated here was made

The H2 abundance is obtained from a **local** balance. In
`src/modules/nonlinear_system_solver/System_HeH_mol.f90` the H2 row of the
coupled system is

```
	! (4) H2 balance
	fvec(4) = k6*n_e*n_h3p + k9*n_h2p*n_hi + k11*n_h3p*n_hi           &
	        + k15*n_hi*n_hi                                           &
	        - ( g_h2 + g_lw + (k10 + k13)*n_hii + k12*ntot + k14*n_e  &
	          + k8*n_h2p + (k17 + k20 + k23)*n_heii + k18*n_hehp      &
	          + k_ion_H2*n_heiTR )*n_h2
```

which is production minus loss set to zero in the cell: there is no
`d(n_H2 u r^2)/dr` term in the row, and no species diffusion velocity enters
it. This is the closure difference `supersonic_molecular_base.md` section 13.6
measured from the other direction (the omitted advection term is 4 to 500 times
the net chemical rate in the region where `f(H2)` falls).

So of the three, two solve H2 with transport and one solves it locally. That
is a statement about the codes, not yet about the profiles.

---

## 3. Passages on H2, H3+ and temperature, quoted from the published text

Koskinen et al. (2022), section 3.2 (composition of Model A at 0.05 au):

> "The composition is shown by Figure 8. The volume mixing ratio of H2 at the
> lower boundary is about q0 = 0.84 while H is a minor species with q0 = 0.026.
> H becomes the dominant species at around r ≈ 2 Rp, but the density of H2
> remains significant at all altitudes. The dominant ion is H+, with a density
> that is at least 2-3 orders of magnitude higher than the densities of other
> ions. Charge balance means that the electron density is equal to the proton
> density. The density of neutral H is higher than the proton density at all
> radii."

Koskinen et al. (2022), section 3.2.1 (comparison with hot Jupiter
simulations):

> "Models of HD209458b predict significant dissociation of H2 near the bottom
> of the thermosphere that removes practically all H2 molecules and related
> molecular ions from the outflow (Yelle 2004; Koskinen et al. 2013a). In
> contrast to this, the density of H2 is significant at all radii of the hot
> Uranus model. This leads to the appearance of the molecular ions H+3 and
> HeH+ at high altitudes, neither of which are significantly produced in hot
> Jupiter simulations."

Koskinen et al. (2022), section 3.2, on the temperature and velocity of
Model A (their Figure 7):

> "The lower boundary of the model is at r0 = 1.34Rp (p0 = 10−6 bar,
> T0 = 1, 140 K) where the boundary conditions are consistent with our models
> of the lower and middle atmosphere (see Section 2.2). Above the lower
> boundary, the temperature first slightly decreases with the radius, and then
> increases rapidly to about 4400 K around r ≈ 3 Rp. At higher radii, the
> temperature increases slowly with the radius to about 5540 K at the upper
> boundary of the model."

Frelikh & Murray-Clay (2026), section 4:

> "In a hot Jupiter, the primary mode of dissociation that determines the
> transition from molecular hydrogen to atomic, thus shutting off of H+3
> cooling, is thermal dissociation of H2. As the temperature threshold for
> significant thermal dissociation is never reached in a super-Earth, molecular
> hydrogen is present throughout, instead of being confined to a layer at the
> base. Nevertheless, the atmosphere is quite optically thin to ionizing
> photons, allowing H+3 to accumulate where H2 is present. Thus, H+3 is an
> important radiative coolant for super-Earth-sized planets."

Frelikh & Murray-Clay (2026), caption of Figure 21:

> "The wind is mostly neutral H at the sonic point. The H2 is not cut off
> sharply as in the hot Jupiter case; its gradual drop off allows H+3 to be
> present throughout."

Frelikh & Murray-Clay (2026), caption of Figure 20:

> "Density (panel (a)), velocity (panel (b)), pressure (panel (c)), and
> temperature (panel (d)) profiles for our fiducial super-Earth. The velocity
> reaches 5.8 km s−2 at 4.4 planetary radii. The temperature peaks at 4200 K,
> and is 4000 K at the sonic point."

Frelikh & Murray-Clay (2026), section 5.4, comparing their super-Earth with
Koskinen's hot Uranus:

> "Our super-Earth runs are qualitatively similar to the results of their hot
> Uranus model, in which instead of the sharp transition between the molecular
> hydrogen layer and the atomic layer, as seen in our hot Jupiter simulations
> (and also in R. V. Yelle 2004; and A. García-Muñoz 2007), molecules are
> present throughout the outflow. Particularly interesting is the presence of
> H+3 throughout the outflow, instead of being confined to a narrow region at
> the base. In agreement with T. T. Koskinen et al. (2022), this is largely due
> to the lower temperature in the wind. Our temperature range is 1000–3800 K,
> in comparison to their temperature prediction of 4000–5000 K for the hot
> Uranus, continuing the trend that a planet with lower gravity (and similar
> stellar and orbital parameters) will host a lower-temperature outflow."

For the hot Jupiter (not the super-Earth), the same paper gives sharp
transition radii in the caption of its Figure 4, which is a useful contrast:

> "The transition from H2 (green)→H (dark blue) occurs at 1.01RP, H (dark
> blue) →H+ (light blue) at 1.8RP, and He (black) →He+ (gray) at 1.4RP. The
> molecular layer drops off sharply at 1.04RP."

---

## 4. Digitization: how the points were read, and how the reading was checked

Each figure page was rendered with `pdftoppm -r 300 -png`. Inside each panel
the axis frame and the tick marks were located by their pixel rows and columns,
so the axis calibration comes from the printed ticks rather than from a
by-eye placement of the corners. Every curve in these figures is drawn in a
distinct RGB color, so each curve was extracted by an exact color match
(tolerance: sum of absolute RGB differences < 40) and reduced to one value per
pixel column, the mean of the matching rows in that column. Legend boxes were
masked out. A 7-point running median is then applied along the columns; it
removes isolated single-column artifacts (a stray gray glyph landing in the
same column as the gray electron curve, for instance) and changes nothing else
at the precision quoted here.

Resolution of the readings:

| figure | quantity axis | pixels per decade (or per unit) | one pixel is | line width |
|---|---|---|---|---|
| Koskinen Fig. 8 | log n [m^-3] | 55.7 px/dex | 0.018 dex | 4-5 px |
| Koskinen Fig. 8 | r/R_p | 127.1 px per unit | 0.008 R_p | -- |
| Koskinen Fig. 7 | T [K] | 119.4 px per 1000 K | 8 K | 3-4 px |
| Koskinen Fig. 13 | log n [m^-3] | 56.5 px/dex | 0.018 dex | 3-4 px |
| Frelikh Fig. 21 | log n [cm^-3] | 43.4 px/dex | 0.023 dex | 6-7 px |
| Frelikh Fig. 20(d) | log T [K] | 480.5 px/dex | 0.002 dex | 3-4 px |

A realistic uncertainty is therefore about **+/- 0.05 dex** on Koskinen's
densities, **+/- 0.1 dex** on Frelikh's, and a few tens of K on both
temperatures, away from steep parts of the curves.

**Where the reading is not reliable.** In the first one or two percent in
radius above each model's lower boundary the curves are essentially vertical,
and the column mean there is an average over several decades rather than a
value. Those entries are marked in the tables and excluded from the figure.

**Calibration checks.** Quantities the papers state in words were recomputed
from the digitized curves; agreement is the evidence that the axis calibration
and the unit conversions are right.

| check | stated in the paper | recovered from the digitization |
|---|---|---|
| Koskinen q0(H2) at the lower boundary | 0.84 | 0.839 |
| Koskinen q0(H) at the lower boundary | 0.026 | 0.028 |
| Koskinen He/H by number | 7.93e-2 | 7.84e-2 |
| Koskinen T at the lower boundary | 1140 K | 1130 K (first column of Fig. 7) |
| Koskinen T "about 4400 K around r ≈ 3 Rp" | 4400 K | 4379 K at 3.00 R_p |
| Koskinen T at the upper boundary | 5540 K | 5527 K at 9.60 R_p |
| Koskinen v at the upper boundary | 9.8 km s^-1 | 9.77 km s^-1 |
| Frelikh T at the base | 1000 K (T0, Table 3) | 1000 K at the first column |
| Frelikh "The temperature peaks at 4200 K" | 4200 K | 4170 K, at r = 2.54e9 cm |
| Frelikh "is 4000 K at the sonic point" | 4000 K | T = 4000 K at r = 2.0e9 cm |

One internal inconsistency in the source is worth recording, since it affects
how their Figure 20 is read: the caption says the velocity reaches
5.8 km s^-1 "at 4.4 planetary radii", while the red marker in their panel (b)
sits at r ≈ 4.4e9 cm, which is 3.6 R_p for R_p = 1.9 R_E. The two readings of
"4.4" cannot both hold; the marker position is what the panel shows. Their
section 5.4 likewise quotes their own temperature range as 1000-3800 K where
the Figure 20 caption and the panel both give a peak of about 4200 K.

---

## 5. Digitized profiles

### 5.1 Koskinen et al. (2022) Model A, a = 0.05 au

Densities read from **Figure 8** (single panel: number density in m^-3 against
r/R_p); the tabulated values are log10 n in **cm^-3**, i.e. the reading minus 6.
`Hp` is the light purple curve, which the caption says is "practically
identical to the electron (em) density profile" and which is hidden underneath
the gray `em` curve over most of the range; `em` is therefore the column to use
for both. T and v are read from **Figure 7** (single panel, T on the left axis,
V on the right axis; the solid curve is T and the dashed curve is V).

| r/R_p | r/r_base | H2 | H | He | H+ | e- | He+ | H3+ | H2+ | HeH+ | T [K] | v [km/s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.34* | 1.000 | 12.69 | 11.21 | 11.89 | 4.59 | 4.82 | 0.17 | 4.41 | -2.83 | -3.39 | 1067 | -- |
| 1.40* | 1.045 | 11.84 | 10.42 | 11.08 | 4.84 | 6.71 | 1.87 | 3.43 | -0.84 | -1.67 | 1096 | -- |
| 1.50 | 1.119 | 10.86 | 9.59 | 10.11 | -- | 6.72 | 2.32 | 3.86 | 0.70 | -1.47 | 1611 | -- |
| 1.60 | 1.194 | 10.12 | 9.22 | 9.36 | -- | 6.96 | 3.15 | 3.94 | 1.78 | -0.51 | 1884 | -- |
| 1.75 | 1.306 | 9.31 | 8.88 | 8.60 | -- | 7.01 | 3.80 | -- | 2.73 | 0.14 | 2286 | 0.09 |
| 2.00 | 1.493 | 8.38 | 8.38 | 7.77 | -- | 6.77 | 4.07 | 4.03 | 3.28 | 0.26 | 2855 | 0.30 |
| 2.30 | 1.716 | 7.82 | 7.94 | 7.24 | -- | 6.48 | 4.05 | 3.97 | 3.55 | 0.22 | 3517 | 0.93 |
| 2.70 | 2.015 | 7.31 | 7.56 | 6.81 | -- | 6.24 | 4.04 | 3.86 | 3.61 | 0.14 | 4141 | 1.91 |
| 3.20 | 2.388 | 6.92 | 7.26 | 6.45 | 6.03 | 6.07 | 4.04 | 3.61 | 3.55 | 0.09 | 4488 | 2.92 |
| 4.00 | 2.985 | 6.51 | 6.92 | 6.09 | 5.87 | 5.89 | 4.05 | 3.29 | 3.48 | 0.03 | 4740 | 4.36 |
| 5.00 | 3.731 | 6.13 | 6.61 | 5.79 | 5.69 | 5.73 | 4.04 | 3.00 | 3.36 | -- | 4924 | 5.69 |
| 6.00 | 4.478 | 5.86 | 6.39 | 5.53 | 5.56 | 5.61 | 4.01 | 2.77 | 3.27 | -0.11 | 5087 | 6.84 |
| 7.50 | 5.597 | 5.55 | 6.12 | 5.24 | 5.40 | 5.43 | 3.95 | 2.52 | 3.16 | -0.23 | 5292 | 8.21 |
| 9.00 | 6.716 | -- | 5.90 | 5.03 | 5.26 | 5.30 | 3.88 | 2.34 | 3.08 | -0.33 | 5477 | 9.31 |
| 9.60 | 7.164 | 5.19 | 5.82 | 4.96 | 5.20 | 5.25 | 3.87 | 2.27 | 3.05 | -0.37 | 5527 | 9.70 |

\* the two rows nearest the lower boundary fall in the near-vertical part of
every curve; treat them as order-of-magnitude only. Entries "--" are columns
where the curve is hidden under another curve of the same darkness or leaves
the plotted range.

Derived from the same reading:

* volume mixing ratio n(H2)/n_tot = 0.839 at the lower boundary (0.830 with
  the median filter of section 4; the spread is the base-column caveat), 0.76 at
  1.20 r_base, 0.47 at 1.50, 0.31 at 2.00, 0.23 at 3.00, 0.14 at the top of
  the domain (7.16 r_base). It never falls by a full decade inside the domain.
* the fraction of H nuclei in H2, `f(H2) = 2 n(H2) / n(H nuclei)`: 0.984 at the
  lower boundary, 0.94 at 1.20, 0.69 at 1.50, 0.51 at 2.00, 0.40 at 3.00,
  0.27 at 7.16.
* n(H3+) has a **broad maximum of 1.2e4 cm^-3**, within 10 percent of that
  value everywhere between 1.23 and 1.54 r_base (1.65-2.06 R_p), and is still
  1.06e4 at 1.50 r_base, 7.2e3 at 2.00, 1.9e3 at 3.00 and 1.9e2 at 7.16.
* n_e (= the `em` curve) peaks at 1.1e7 cm^-3 near 1.25 r_base (within
  10 percent of that between 1.21 and 1.32) and is 1.8e5 cm^-3 at the top of
  the domain.

### 5.2 Koskinen et al. (2022) Figure 13: the revised Model A

Their Figure 13 repeats the Model A composition with an approximate treatment
of H2 neutral dissociation by photons longward of 80.35 nm switched on (their
Appendix B; the caption says "The results are qualitatively similar to the
model that does not include neutral dissociation of H2"). It is tabulated here
because it bounds how much of their H2 survival depends on that omitted
channel. Same units and conventions as section 5.1.

| r/R_p | r/r_base | H2 | H | He | H+ | e- | He+ | H3+ | H2+ | HeH+ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.34* | 1.000 | 12.73 | 11.22 | 11.93 | 4.57 | 4.80 | 0.05 | 4.44 | -2.84 | -3.50 |
| 1.40* | 1.045 | 12.00 | 10.58 | 11.21 | 4.94 | 6.58 | 2.04 | 3.45 | -1.00 | -1.54 |
| 1.50 | 1.119 | 10.73 | 9.92 | 10.03 | -- | 6.56 | 2.33 | 3.93 | 0.81 | -1.27 |
| 1.60 | 1.194 | 10.07 | 9.65 | 9.35 | -- | 6.88 | 3.11 | 3.95 | 1.79 | -0.55 |
| 1.75 | 1.306 | 9.28 | 9.14 | 8.63 | -- | 7.03 | 3.78 | 3.96 | 2.70 | 0.06 |
| 2.00 | 1.493 | 8.50 | 8.61 | 7.92 | -- | 6.87 | 4.07 | 3.90 | 3.19 | 0.25 |
| 2.30 | 1.716 | 7.88 | 8.14 | 7.37 | -- | 6.61 | 4.11 | 3.85 | 3.45 | 0.21 |
| 2.70 | 2.015 | 7.36 | 7.75 | 6.92 | -- | 6.36 | 4.11 | 3.73 | 3.58 | 0.14 |
| 3.20 | 2.388 | 6.92 | 7.43 | 6.56 | -- | 6.18 | 4.12 | 3.50 | 3.57 | 0.07 |
| 4.00 | 2.985 | 6.48 | 7.10 | 6.19 | 5.96 | 6.00 | 4.14 | 3.12 | 3.44 | 0.00 |
| 5.00 | 3.731 | 6.07 | 6.80 | 5.88 | 5.80 | 5.83 | 4.14 | 2.76 | 3.29 | -0.09 |
| 6.00 | 4.478 | 5.78 | 6.57 | 5.62 | 5.65 | 5.70 | 4.12 | 2.46 | 3.16 | -0.18 |
| 7.50 | 5.597 | 5.42 | 6.30 | 5.35 | 5.50 | 5.54 | 4.08 | 2.12 | 2.99 | -0.32 |
| 9.00 | 6.716 | 5.15 | 6.09 | 5.11 | 5.37 | 5.42 | 4.02 | 1.86 | 2.87 | -0.46 |
| 9.60 | 7.164 | 5.04 | 6.02 | -- | 5.32 | 5.36 | 4.00 | 1.75 | 2.81 | -0.51 |

Read against section 5.1, n(H2) is lower by 0.1-0.2 dex over most of the
domain and n(H3+) by up to 0.5 dex at the top, consistent with their statement
that "The abundance of H2 is somewhat reduced in the revised model and the
mass-loss rate is slightly higher but the results are qualitatively similar and
the mass-loss rates agree to within a factor of 1.4."

### 5.3 Frelikh & Murray-Clay (2026), fiducial super-Earth

Densities read from **Figure 21** (single panel: number density in cm^-3
against radius in cm), tabulated as log10 n in cm^-3. `n_e` is not a curve in
the figure; it is the sum H+ + He+ + H2+ + H3+ + HeH+ formed from the same
reading. T is read from **Figure 20, panel (d)** (temperature against radius,
log axis). r_base = 1.9 R_E = 1.2108e9 cm.

| r [1e9 cm] | r/r_base | H2 | H | He | H+ | He+ | H3+ | H2+ | HeH+ | n_e | T [K] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.22* | 1.008 | 12.57 | 11.16 | 11.81 | 7.09 | 3.70 | 5.30 | 1.57 | -0.02 | 7.09 | 935 |
| 1.25* | 1.033 | 11.85 | 11.58 | 11.11 | 8.34 | 4.68 | 4.15 | 2.23 | 0.96 | 8.34 | 967 |
| 1.30 | 1.074 | 10.59 | 10.84 | 10.10 | 8.46 | 6.16 | 3.88 | 3.12 | 2.32 | 8.47 | 1666 |
| 1.35 | 1.115 | 9.95 | 10.46 | 9.59 | 8.27 | 6.84 | 4.14 | 3.86 | 2.80 | 8.28 | 2237 |
| 1.45 | 1.198 | 9.11 | 9.99 | 9.01 | 8.15 | 7.37 | 3.66 | 4.21 | 2.99 | 8.22 | 2769 |
| 1.60 | 1.322 | 8.22 | 9.53 | 8.45 | 7.91 | 7.43 | 3.02 | 4.19 | 2.72 | 8.04 | 3329 |
| 1.80 | 1.487 | 7.54 | 9.11 | 7.98 | 7.69 | 7.27 | 2.42 | 4.02 | 2.30 | 7.83 | 3762 |
| 2.10 | 1.735 | 6.88 | 8.69 | 7.52 | 7.42 | 7.01 | 1.67 | 3.75 | 1.81 | 7.56 | 4052 |
| 2.50 | 2.065 | 6.38 | 8.30 | 7.13 | 7.32 | 6.74 | 1.04 | 3.42 | 1.35 | 7.42 | 4161 |
| 3.00 | 2.478 | 5.94 | 7.97 | 6.76 | 7.14 | 6.48 | 0.48 | 3.10 | 0.97 | 7.23 | 4151 |
| 3.50 | 2.891 | 5.63 | 7.72 | 6.51 | 7.00 | 6.30 | 0.09 | 2.87 | 0.68 | 7.08 | 4111 |
| 4.00 | 3.304 | 5.40 | 7.51 | 6.30 | 6.89 | 6.14 | -0.21 | 2.69 | 0.47 | 6.96 | 4052 |
| 4.50 | 3.717 | 5.20 | 7.34 | 6.14 | 6.78 | 6.00 | -0.44 | 2.55 | 0.28 | 6.85 | 4014 |
| 5.00 | 4.131 | 5.03 | 7.19 | 6.00 | 6.69 | 5.89 | -0.63 | 2.43 | 0.13 | 6.75 | 3975 |

\* as before, the first two rows sit in the near-vertical base columns.

Derived from the same reading:

* n(H2)/n_tot = 0.82-0.85 at the first readable columns (0.85 raw, 0.82 after
  the median filter), 0.10 at 1.20 r_base, 0.023 at 1.50, 0.010 at 2.00,
  0.0061 at 3.00.
* `f(H2)` = 2 n(H2)/n(H nuclei): 0.98 at the first columns, 0.20 at 1.20,
  0.048 at 1.50, 0.022 at 2.00, 0.013 at 3.00 and 0.011 at the outer boundary
  (4.13 r_base). It falls by one decade at 1.30 r_base and does not reach two
  decades of fall inside the plotted domain.
* n(H3+) has a local maximum of **1.5e4 cm^-3 at 1.10 r_base** (the curve also
  spikes at the base itself, in the unreadable columns) and is 2.4e2 at
  1.50 r_base, 14 at 2.00, 1.0 at 3.00 and 0.24 at the outer boundary.
* n(He nuclei)/n(H nuclei) at the base of the readable range is 0.085. This
  is the only handle the paper gives on its helium abundance, and it is read
  from the steepest part of the figure, so it should be treated as
  "near solar" rather than as a number.

---

## 6. The three models on one radial grid

`r/r_base` is each model's own lower boundary set to 1 (section 1). "ours" is
`heh0p0793c_r1`, physical cells only. `f(H2)` is the fraction of H nuclei bound
in H2, `q(H2) = n(H2)/n_tot` the volume mixing ratio against all heavy
particles (electrons excluded), both computed the same way for all three.
Published columns are interpolated in log space from the digitized curves,
after a 7-point running median along the pixel columns that removes isolated
single-column artifacts. Entries at `r/r_base = 1.00` come from the
near-vertical base columns of the published figures and carry the caveat of
section 4; there the filtered and unfiltered readings differ by up to 0.4 dex,
which is itself the measure of how little that column can be trusted. From
1.05 outward the filter changes nothing at the printed precision.

| r/r_base | 1.00 | 1.05 | 1.10 | 1.15 | 1.20 | 1.30 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|---|
| **f(H2)**, ours | 0.999 | 0.966 | 0.516 | 1.42e-2 | 1.10e-3 | 7.2e-5 | 2.8e-5 | 6.4e-6 | 1.5e-6 |
| **f(H2)**, Koskinen Fig. 8 | 0.984 | 0.982 | 0.976 | 0.963 | 0.936 | 0.843 | 0.688 | 0.508 | 0.403 |
| **f(H2)**, Frelikh Fig. 21 | 0.976 | 0.606 | 0.438 | 0.269 | 0.200 | 0.096 | 0.048 | 0.022 | 0.013 |
| **q(H2)**, ours | 0.861 | 0.810 | 0.314 | 6.6e-3 | 5.1e-4 | 3.3e-5 | 1.3e-5 | 2.9e-6 | 7.0e-7 |
| **q(H2)**, Koskinen Fig. 8 | 0.830 | 0.822 | 0.825 | 0.801 | 0.758 | 0.642 | 0.467 | 0.308 | 0.230 |
| **q(H2)**, Frelikh Fig. 21 | 0.817 | 0.387 | 0.254 | 0.142 | 0.102 | 0.046 | 0.023 | 0.010 | 6.1e-3 |
| **n(H3+)** [cm^-3], ours | 2.9e5 | 2.0e5 | 1.1e3 | 13.9 | 0.18 | 6.4e-3 | 7.9e-4 | 1.7e-5 | 6.9e-7 |
| **n(H3+)**, Koskinen Fig. 8 | 1.9e4 | 3.1e3 | 5.0e3 | 7.0e3 | 9.1e3 | 1.2e4 | 1.1e4 | 7.2e3 | 1.9e3 |
| **n(H3+)**, Frelikh Fig. 21 | 1.5e5 | 7.7e3 | 1.5e4 | 8.8e3 | 4.5e3 | 1.3e3 | 2.4e2 | 14.2 | 0.99 |
| **T** [K], ours | 1100 | 795 | 775 | 502 | 759 | 1610 | 2470 | 2720 | 2250 |
| **T**, Koskinen Fig. 7 | 1082 | 1180 | 1570 | 1720 | 1900 | 2270 | 2870 | 4120 | 4740 |
| **T**, Frelikh Fig. 20(d) | 939 | 1140 | 2090 | 2470 | 2770 | 3250 | 3790 | 4150 | 4090 |
| **n_e** [cm^-3], ours | 2.9e5 | 5.5e5 | 3.3e7 | 6.9e7 | 7.7e7 | 7.1e7 | 6.8e7 | 4.6e7 | 2.0e7 |
| **n_e**, Koskinen Fig. 8 (`em`) | 1.6e5 | 5.0e6 | 4.4e6 | 7.2e6 | 9.6e6 | 1.0e7 | 5.7e6 | 1.8e6 | 7.8e5 |
| **n_e**, Frelikh Fig. 21 (ion sum) | 1.6e7 | 2.8e8 | 2.0e8 | 1.9e8 | 1.6e8 | 1.2e8 | 6.6e7 | 2.8e7 | 1.1e7 |

Landmarks of the same three profiles, in the same normalization:

| | ours | Koskinen Model A | Frelikh super-Earth |
|---|---|---|---|
| `f(H2)` down by 1 decade at | 1.132 r_base | not reached in the domain (7.16 r_base) | 1.296 r_base |
| `f(H2)` down by 2 decades at | 1.156 | not reached | not reached (4.13 r_base) |
| `f(H2)` down by 4 decades at | 1.273 | not reached | not reached |
| `f(H2)` at the outer boundary | 3.5e-7 at 4.73 | 0.269 at 7.16 | 0.011 at 4.13 |
| n(H3+) maximum | 3.09e5 cm^-3 at 1.023 | 1.2e4 cm^-3, flat over 1.23-1.54 | 1.5e4 cm^-3 at 1.10 |
| n(H3+) at 2.00 r_base | 1.7e-5 | 7.2e3 | 14.2 |
| n(H3+) at 3.00 r_base | 6.9e-7 | 1.9e3 | 0.99 |
| T at the lower boundary | 1140 (base state) | 1140 (stated), 1130 (digitized) | 1000 (T0), 1000 (digitized) |
| T maximum, and where | 2751 K at 1.815 | 5527 K at 7.16 (the outer boundary; still rising) | 4170 K at 2.10 |
| n_e maximum | 7.8e7 cm^-3 at 1.19 | 1.1e7 cm^-3 near 1.25 | 3.2e8 cm^-3 near 1.07 |

---

## 7. Budget terms: the transport term in the H2 continuity equation

Section 2 established that both published models carry transport in their
species equations and that ours does not. This section reads the terms
themselves off the three budget figures, so that the comparison with the
EXHALE budget of `supersonic_molecular_base.md` section 13.6 is numerical
rather than verbal.

### 7.1 Which curve is which, and how that was decided

**Koskinen Figures 15 and 16** carry the same four-entry legend: "H+ (or H2)
production" black solid, "H+ (or H2) loss" black dashed, "Advection+diffusion"
black dotted, "Loss+Transport" light blue solid. Only the blue is separable by
color. The assignment was made as follows.

* The blue "Loss+Transport" curve lies exactly on one of the black curves over
  the whole domain. The Appendix B text says that this coincidence is the
  steady-state check ("the near perfect alignment of the proton production rate
  (solid black line) with the sum of the proton loss and transport rates (solid
  blue line)"), so the black curve underneath the blue is the **production**
  curve. The remaining two black curves were separated top-down in each pixel
  column.
* Both figures plot **absolute values** (their captions say so), and the
  steady-state statement in the text is "the sum of the chemical loss and
  transport terms equals the production rate". For H2 the transport term is a
  supply, i.e. negative, so in the plotted absolute values the identity reads
  `loss = production + |transport|`; for H+ the transport term is a sink and it
  reads `production = loss + |transport|`. Both forms were used as the check.
* In **Figure 16** the two remaining black curves never cross: the upper one
  was taken to be the **loss** (dashed) and the lower one the
  **advection+diffusion** (dotted). The readings satisfy the identity to
  1.4-3 percent at every radius tested (for example at r = 2.01 R_p,
  10^3.474 = 2980 against 10^2.799 + 10^3.362 = 2930). Line style was confirmed
  visually on a 600 dpi crop: the upper curve has long dashes, the lower has
  short dots.
* In **Figure 15** the two cross. Below the crossing the dashed (loss) is
  above the dotted (transport); above it the order reverses, which was checked
  on a 600 dpi zoom of r = 2.8-6.0 R_p where the middle curve is visibly dotted
  and the lower one visibly dashed. The identity holds to 1.4-2 percent at
  r = 2.0, 2.7 and 3.2 R_p.

**Frelikh Figure 7** gives every channel its own color, and the legend was
rendered at 600 dpi and matched swatch by swatch:
blue = `H2 + hv -> H2+ + e`; black = `H2 + M -> H + H + M`;
orange = `H2 + hv(FUV) -> H + H`; red solid = `H2+ + H2 -> H3+ + H`;
red dotted = `H3+ + H -> H2+ + H2`; gray solid = `H+ + H2 -> H2+ + H`;
gray dotted = `H2+ + H -> H+ + H2`; steel blue = `H2 + e_ph -> H2+ + e'_ph + e`;
light blue dashed = `- advection of H2`. The caption fixes the two
interpretations that matter: "The gray solid and dotted lines represent the
process of H2 and H exchanging an electron: this reaction pair does not lead to
the destruction of H2", and, for the advection curve, "its sign follows the
convention established in Figure 5: as the sign is negative, the advective term
represents an addition of H2 from the boundary."

**Two things about Frelikh Figure 7 have to be said before its numbers are
used.** First, it belongs to their **fiducial hot Jupiter** (section 3.1.2,
Table 2), not to the super-Earth of section 4; the paper publishes no
budget figure for the super-Earth. Second, its title is "*Selected* reaction
rates involving H2": the H2-returning channels `H3+ + e -> H2 + H` (their k6)
and `H + H + M -> H2 + M` (k13) are not drawn, so a true net chemical rate
cannot be formed from the figure. What is written below as "sum of the plotted
sinks" is
`H2+hv` + `H2+e_ph` + `H2+hv(FUV)` + `H2+M` + `H2+ + H2 -> H3+ + H`
minus `H3+ + H -> H2+ + H2`, with the gray pair excluded on the caption's
instruction. It is an upper bound on the net, not the net.

### 7.2 Koskinen et al. (2022) Model A, H2 (their Figure 16)

`P` = H2 production, `L` = H2 gross loss, `|T|` = advection+diffusion, all in
cm^-3 s^-1, read at the radii of the EXHALE budget table and a few beyond.
`n(H2)` comes from their Figure 8 (section 5.1). Their figure spans
r/R_p = 1.4-10, i.e. `r/r_base` = 1.045-7.46.

| r/r_base | r/R_p | P | L | \|T\| | \|T\|/P | \|T\|/L | n(H2)/\|T\| [s] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.05 | 1.407 | 1.4e2 | 1.4e3 | 1.2e3 | 8.6 | 0.81 | 5.0e8 |
| 1.15 | 1.541 | 4.0e2 | 3.2e3 | 2.7e3 | 6.7 | 0.84 | 1.2e7 |
| 1.25 | 1.675 | 9.7e2 | 6.4e3 | 5.5e3 | 5.7 | 0.86 | 9.7e5 |
| 1.50 | 2.010 | 6.1e2 | 2.9e3 | 2.3e3 | 3.7 | 0.79 | 1.0e5 |
| 2.00 | 2.680 | 1.6e2 | 5.9e2 | 4.3e2 | 2.8 | 0.74 | 4.7e4 |
| 3.00 | 4.020 | 2.2e1 | 8.6e1 | 6.4e1 | 2.9 | 0.75 | 4.7e4 |

Two features of the reading, stated without interpretation:

* `|T|` exceeds `P` at **every radius the figure covers**, by a factor 2.8 to
  8.6, and `|T|/L` stays between 0.74 and 0.88 throughout. There is no crossing
  inside the plotted range: whatever radius the transport term first becomes
  comparable to the chemistry, it is at or below the left edge of their figure,
  r/R_p = 1.4 (`r/r_base` = 1.045).
* In a steady state the identity `|T| = L - P` holds by construction, so the
  ratio "transport over *net* chemical rate" is 1 everywhere in their model by
  definition and carries no information. The informative ratios are the two
  tabulated.

### 7.3 Koskinen et al. (2022) Model A, H+ (their Figure 15)

For protons the loss (radiative recombination) and the transport term cross.
From the pixel separation of the two black curves, the gap between them is a
minimum at **r = 2.95 R_p, i.e. `r/r_base` = 2.20**, where the production rate
is 2.2e2 cm^-3 s^-1 and each of loss and transport is therefore about
1.1e2 cm^-3 s^-1. Inside that radius recombination is the larger sink; outside
it the transport term is. This is the figure behind the sentence in
Appendix B, "Instead of recombination, ionization of H is balanced by upward
transport of neutral H from below by escaping gas."

Values at r <= 3.2 R_p, where the two curves are far enough apart to separate:

| r/R_p | r/r_base | production | loss | \|T\| |
|---:|---:|---:|---:|---:|
| 2.00 | 1.49 | 9.9e2 | 7.7e2 | 2.4e2 |
| 2.30 | 1.72 | 5.5e2 | 4.0e2 | 1.6e2 |
| 2.70 | 2.01 | 3.0e2 | 1.9e2 | 1.3e2 |
| 3.20 | 2.39 | 1.7e2 | 7.6e1 | 9.5e1 |

Beyond about 4 R_p the dashed and dotted curves run close together and the
column separation stops being reliable; those radii are not quoted.

### 7.4 Frelikh & Murray-Clay (2026), H2 in their fiducial hot Jupiter (their Figure 7)

Rates in cm^-3 s^-1; `r_base = R_p = 1e10 cm` for this model, so the abscissa
of their figure is already `r/r_base`. `n(H2)` comes from their Figure 4
panel (a), digitized here for this purpose (green curve, same color key as
their Figure 21).

| r/r_base | H2+hv | H2+e_ph | H2+hv(FUV) | H2+M | H2+ +H2 | sum of plotted sinks | \|adv\| | \|adv\|/sum | n(H2)/\|adv\| [s] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.005 | 9.1e2 | 9.6e3 | 8.3e2 | -- | 6.3e4 | 7.4e4 | 1.5e5 | 2.04 | 9.6e6 |
| 1.008 | 9.4e2 | 9.8e3 | 8.6e2 | -- | 5.0e4 | 6.2e4 | 7.0e4 | 1.13 | 6.6e6 |
| 1.010 | 8.1e2 | 8.1e3 | 6.6e2 | -- | 4.7e4 | 5.6e4 | 4.6e4 | 0.82 | 5.7e6 |
| 1.012 | 6.9e2 | 6.3e3 | 4.8e2 | 1.7 | 4.4e4 | 5.1e4 | 3.0e4 | 0.59 | 5.2e6 |
| 1.016 | 5.5e2 | 3.8e3 | 2.3e2 | 5.3 | 3.8e4 | 4.2e4 | 1.5e4 | 0.35 | 4.4e6 |
| 1.020 | 5.8e2 | 2.3e3 | 1.1e2 | 1.4e1 | 3.2e4 | 3.5e4 | 7.7e3 | 0.22 | 3.5e6 |
| 1.025 | 4.9e2 | 9.4e2 | 3.0e1 | 8.8e1 | 2.4e4 | 2.6e4 | 2.9e3 | 0.11 | 2.5e6 |
| 1.030 | 1.4e2 | 1.7e2 | 4.0e0 | 1.0e3 | 7.6e3 | 8.8e3 | 3.0e2 | 0.034 | 3.2e6 |

* **The advection term equals the sum of the plotted sinks at
  `r/r_base` = 1.0084** (interpolated ratio 1.00), is twice the sum at 1.005
  and a thirtieth of it at 1.030. Their molecular layer ends at about
  1.04 r_base, so on this reading the transport-dominated part of it is the
  inner ~0.01 R_p.
* The gray pair (the H2/H electron exchange, which the caption says does not
  destroy H2) is the largest plotted rate for `r/r_base` > 1.012, reaching
  about 2e5 cm^-3 s^-1 near 1.030. The two gray curves nearly coincide over
  much of that range and the column separation of solid from dotted is not
  reliable there, so only the upper envelope is quoted.
* Below about `r/r_base` = 1.005 every curve in their Figure 7 is drawn in a
  faded tint, the deemphasis their Figures 5 and 6 use for the bolometrically
  heated region. Those columns were not read.

### 7.5 The EXHALE budget beside them

From `supersonic_molecular_base.md` section 13.6, with `n(H2)` taken from the
same converged run:

| r/r_base | n(H2) [cm^-3] | net chemical H2 loss | \|advection term\| | ratio | n(H2)/\|adv\| [s] | n(H2)/net [s] |
|---:|---:|---:|---:|---:|---:|---:|
| 1.05 | 6.0e12 | 6.75e3 | 2.96e5 | 44 | 2.0e7 | 8.9e8 |
| 1.15 | 1.8e9 | 1.33e2 | 1.55e4 | 117 | 1.1e5 | 1.3e7 |
| 1.25 | 7.0e5 | 8.51e1 | 4.74e1 | 0.6 | 1.5e4 | 8.2e3 |

The ratio crosses unity between 1.15 and 1.25 `r_base` in our run.

**The denominators are not the same statistic in the three models, and that
has to be carried with any use of these numbers.** Koskinen's "H2 loss" is a
*gross* rate: it includes the reactions that give the H2 straight back, and his
production curve includes their return legs. Our "net chemical H2 loss" is
*cycle-resolved*: section 13.6 assigns each molecular ion its expected net H2
yield and charges each primary channel only the H2 it actually removes, which
is why the null cycle R20 -> R19 -> R9 does not appear in it. Frelikh's column
is a sum over the *plotted* sinks only, and is an upper bound on their net.
So `|T|/L` = 0.74-0.88 (Koskinen), `|adv|/sum` = 2.0 down to 0.03 (Frelikh) and
`|adv|/net` = 44, 117, 0.6 (ours) are three different ratios, and the table of
section 7.6 is the like-for-like one.

### 7.6 The one quantity all three define the same way

`n(H2)` divided by the transport term is the time in which the flow replaces
the H2 in a cell, and it needs no convention about which chemical channels
count. Beside it, `n(H2)` divided by the chemical rate is the time the
chemistry needs to change the H2 there.

| | 1.05 r_base | 1.15 r_base | 1.25 r_base |
|---|---:|---:|---:|
| Koskinen, n(H2)/\|T\| | 5.0e8 s | 1.2e7 s | 9.7e5 s |
| Koskinen, n(H2)/L (gross) | 4.1e8 s | 1.0e7 s | 8.4e5 s |
| EXHALE, n(H2)/\|advection\| | 2.0e7 s | 1.1e5 s | 1.5e4 s |
| EXHALE, n(H2)/net chemical | 8.9e8 s | 1.3e7 s | 8.2e3 s |
| EXHALE, `H_p/v` (section 13.6) | 1.8e7 s | 2.4e5 s | 2.5e4 s |

Frelikh's hot Jupiter does not reach these radii; over the 1.005-1.030 r_base
range its `n(H2)/|adv|` falls from 9.6e6 to 3.2e6 s and its
`n(H2)/(sum of plotted sinks)` from 1.9e7 to 1.0e5 s.

---

## 8. Figures

`docs/lower_atmosphere_figs/fig_p23_published_vs_exhale.png` -- four panels,
`f(H2)`, `n(H3+)`, `T` and `n_e` against `r/r_base`. The solid black curve is
our converged He/H = 0.0793 state; open circles are Koskinen Model A digitized
from their Figures 8 and 7; open squares are the Frelikh super-Earth digitized
from their Figures 21 and 20(d). Markers begin at 1.05 r_base, above the
near-vertical base columns, and the digitized series carry the same
running-median filter as the table of section 6.

`docs/lower_atmosphere_figs/fig_p23_budget_terms.png` -- the budget terms of
section 7. Panel (a): the three H2 continuity terms of Koskinen Figure 16
against `r/r_base`, with the EXHALE net chemical rate and the omitted advection
term of section 13.6 as points. Panel (b): the transport term divided by a
chemical rate, plotted against `r/r_base - 1` so that all three models fit on
one axis (Koskinen over `|T|`/production and `|T|`/gross loss, Frelikh's hot
Jupiter over the sum of its plotted sinks, ours over the cycle-resolved net) --
the three denominators differ, as section 7.5 sets out. Panel (c): the
individual channels of Frelikh Figure 7 on their own abscissa. Panel (d): the
H2 replacement timescales of section 7.6.

The plotting scripts and the digitized arrays live in the campaign scratch
directory beside the runs (`p23/make_p23_fig.py`, `p23/make_p23_budget_fig.py`,
`p23/fig*_res.npy`, `p23/f1[56]_raw.npy`, `p23/f7_raw.npy`, `p23/f4a_H2.npy`);
neither figure is regenerated by
`docs/lower_atmosphere_figs/make_figures.py`.

---

## 9. What was not read, and what stays open

**Figures not digitized.**

* Koskinen Figure 9 (heating and cooling rates for Model A). Only the two
  numbers stated in the text were taken: the stellar XUV heating peak at
  1.66 R_p (p = 2.1e-9 bar) with a peak volume heating rate of
  1.8e-8 W m^-3, and that "At r ≲ 2.3 Rp, heating is mostly balanced by H+3
  cooling."
* Koskinen Figures 15 and 16 and Frelikh Figure 7 **have now been read**;
  section 7 has them. What is still missing from those three is listed at the
  end of this section.
* Koskinen Figure 10 (T and v for Model B at 0.02 au) and Figures 2-6, 11-12,
  which concern the Roche-lobe-overflow transition and the planet population,
  not the layer structure.
* Frelikh Figures 3, 5, 6, 8-19 (energy balance, grid of models, H3+ spectrum,
  conduction, efficiency) and 22 (super-Earth heating and cooling). Of their
  Figure 4 only panel (a) was read, and only its H2 curve, for the n(H2) that
  section 7.4 needs; the rest of panel (a) and all of panel (b) were not read.
* Frelikh Figure 20 panels (a), (b), (c) -- density, velocity and pressure of
  the super-Earth. Only panel (d) was digitized; panel (b) would give their
  sonic point directly and settle the caption ambiguity noted in section 4.

**Normalization points that stay open.**

* Frelikh's Figure 20 caption places the sonic point at "4.4 planetary radii"
  while the marker in the panel sits near 4.4e9 cm = 3.6 R_p (section 4). Until
  panel (b) is read, the sonic-point radius of their super-Earth should be
  quoted as ~4.4e9 cm, not in R_p.
* Their section 5.4 quotes the super-Earth temperature range as 1000-3800 K
  where their own Figure 20 caption and panel give a peak near 4200 K. Both
  numbers are in the published version; the tables above use the panel.
* Frelikh do not state a helium abundance for the super-Earth. The 0.085
  in section 5.3 is inferred from their Figure 21 at the base and carries the
  digitization uncertainty of that steep region.
* Koskinen's R_p is the 1 bar polar radius while ours is the escape base, so
  every "R_p" in their text is 1.34 times our unit. All joint tables here use
  `r/r_base`; anything lifted from this file into another document should carry
  that conversion with it.
* Our `q_H2_base = 0.75` is a representative photochemical value rather than a
  computed one (the header of the run's `base.inp` says so), so our base
  composition row in section 1 is an input choice, not a measurement.

**Left unread inside the three budget figures (section 7).**

* Frelikh Figure 7 publishes only *selected* channels: `H3+ + e -> H2 + H` and
  `H + H + M -> H2 + M` return H2 and are not drawn, so no true net chemical
  rate can be formed from it. The "sum of the plotted sinks" in section 7.4 is
  an upper bound on the net, not the net.
* The two gray curves of Frelikh Figure 7 (the H2/H electron-exchange pair)
  nearly coincide for `r/r_base` > 1.02, and separating solid from dotted there
  by pixel column is not reliable; only their upper envelope is quoted.
* Columns below `r/r_base` ~ 1.005 in Frelikh Figure 7 are drawn in the faded
  tint their Figures 5 and 6 use for the bolometrically heated region, and were
  not read.
* In Koskinen Figure 15 the dashed (H+ loss) and dotted (advection+diffusion)
  curves run close together beyond about 4 R_p; the separation is quoted only
  out to 3.2 R_p, plus the crossing radius.
* Koskinen Figure 16 starts at r/R_p = 1.4 (`r/r_base` = 1.045). The radius at
  which their transport term first becomes comparable to their chemistry, if it
  is inside their domain at all, lies below the left edge of the published
  figure and cannot be read from it.
* **Frelikh publish no budget figure for their super-Earth.** The only H2
  budget in that paper is their Figure 7, which is the fiducial hot Jupiter, a
  different planet from the super-Earth tabulated in sections 1, 5.3 and 6.


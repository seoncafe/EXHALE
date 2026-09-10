# P23: the heating and cooling budget of the hot Uranus, ours beside Koskinen et al. (2022) Figure 9

**What this file is.** `docs/p23_published_profiles.md` put the two codes'
*profiles* on one radial axis; `supersonic_molecular_base.md` section 14.5
recorded that at a 1 microbar base our converged hot Uranus is 550, 257 and
118 K below Koskinen et al. (2022) at 1.15, 1.20 and 1.50 `r_base` and
2616 K below at 3.00. This file puts the two codes' *energy budgets* on that
same axis, so that the temperature difference can be attributed to a term
rather than guessed at. It also carries the two 1 microbar ladder anomalies of
section 14.4 as far as measurement takes them.

---

## The judgment first

**The temperature difference is very unlikely to be a deficit of heating and
very unlikely to be an excess of radiative cooling. Both codes absorb the same
stellar flux; the difference is what each does with it, and that follows from
whether H2 survives above the base.**

* **Our photoheating is not small.** Our peak volumetric heating is
  2.23e-7 erg cm^-3 s^-1 at 1.111 `r_base`, against Koskinen's 1.84e-7 at
  1.238 `r_base` (their Figure 9, digitized; the paper states 1.8e-8 W m^-3 at
  1.66 R_p = 1.238 `r_base`). Over 1.05 to 3.00 `r_base` the ratio of our
  photoheating to theirs runs 0.50 to 1.81, i.e. the two are the same to a
  factor 2 either way. Candidate (a) is not supported.
* **Our radiative cooling is 1 to 2 percent of theirs over the front**, not
  larger: 1.4e-9 against 8.7e-8 at 1.15 `r_base`. Their radiative loss is
  H3+ infrared (their own text: "At r <~ 2.3 Rp, heating is mostly balanced by
  H+3 cooling"); ours is H recombination plus bremsstrahlung, because our
  `f(H2)` has fallen to 0.31 at 1.05 `r_base` and 2.7e-3 at 1.10, taking H3+
  with it, where theirs is 0.98 and 0.98. Candidate (b) in the form
  "our H3+ cooling function is too efficient" is **refuted twice over**: our
  H3+ formula applied to *their* densities and temperatures reproduces the
  cooling curve they plot to within a factor 1.0 to 2.9 (section 5), and an
  A/B that removes our H3+ cooling entirely moves our temperature by at most
  56 K, in the first cell only (section 6).
* **Where our budget differs is the expansion and advection terms**: our
  adiabatic `p div v` is 1.6 to 2.9 times theirs at the same radius and our
  advection term is 3.8 to 7.9 times theirs over 1.10 to 1.30 `r_base`.
* **The likely reading.** Because we radiate away almost none of the absorbed
  energy and they radiate away most of it, the net radiative energy input per
  unit area is 0.226 W m^-2 here against 0.11 W m^-2 there (their Equation 24),
  a factor 2.05 -- and our mass flux is 3.14e9 g s^-1 sr^-1 against their
  1.5e9, a factor 2.09. (Both are `rho v r^2` per steradian: ours over 4 pi is
  3.95e10 g s^-1, of which the `log10 Mdot = 10.30` the run prints is the half
  the `Rate/2 + Mdot/2` convention takes, and theirs over 4 pi is the
  1.9e10 g s^-1 they quote.) The extra energy that our missing coolant leaves in the
  gas goes into *driving more wind*, not into raising the temperature. On top of
  that, our gas above the front is atomic where theirs is 40 to 95 percent H2,
  so our mean molecular weight is 1.6 to 1.9 times smaller and the same thermal
  energy per unit mass buys proportionally less temperature. Corrected for the
  mean molecular weight alone, our thermal energy per unit mass is 1.1 to 1.8
  times *theirs* over 1.19 to 2.02 `r_base` and only falls below it beyond
  2.4 `r_base` (section 7). **On this evidence the 1.15-2.0 `r_base` part of
  the temperature gap is a composition effect and not an energy-budget deficit,
  and the outer part (3 `r_base` and beyond) is both.**
* **Heat conduction is not the explanation either.** It is off in our runs
  (no `Conduction:` key), and in theirs it is a diagnostic that never exceeds
  2 percent of the local heating anywhere in 1.45 <= r/R_p <= 9.
* **The 100-450 K molecular-layer deficit of `TO_BE_DONE.md` item (X) is not
  the H3+ coolant.** The A/B of section 6 sets the whole H3+ contribution at
  <= 56 K at 1.001 `r_base`, <= 11 K at 1.01 and below 3 K at 1.03 and above.
  Where the trough actually sits (415 K at 1.065 `r_base`) the H3+ cooling is
  already 3e-5 of the local heating.

**The two ladder anomalies split.** On a base grid refined by a factor 2
(section 8), the He/H = 30 loss of three decades of H2 inside the first cell is
**unchanged** (`x2` at cell 1 4.563e-4 against 4.593e-4) and is therefore the
interior chemistry -- helium ions, 93.5 percent of it the single reaction
He+ + H2 -> H+ + H + He -- meeting a handoff that states a composition the
chemistry does not reproduce. The He/H = 0.3 cold first cell and alternating
temperature are **gone** on the refined grid and were a discretization
artifact; that rung's peak n(H3+) moved by a factor 5.4 with them, so the
artifact was not only cosmetic.

Two defects were found on the way and are reported separately in section 9:
the `Heating_breakdown.txt` diagnostic drops 49 percent of the molecular
photoheating in the first cells, and the written state's own energy budget does
not close below 1.08 `r_base`. **Both have since been dealt with** -- the first
by a fix in the tree, the second by the caloric EOS plus a state finished on
the face-flux measure -- and section 10 re-measures the whole channel table on a
P53 state. The verdict above is unchanged by that re-measurement.

---

## 1. Which figure this is, and what it contains

Koskinen et al. (2022, ApJ 929, 52), **Figure 9**, captioned

> "Heating (red) and cooling (blue) rates for a model of a Uranus-like planet
> orbiting a Sun-like star at the orbital distance of 0.05 au."

The panel holds **five entries, not a channel decomposition**: `Stellar
heating` (red solid), `Radiative cooling` (blue dashed), `Adiabatic` (blue
dash-dot), `Advection` (blue dotted) and `Conduction` (dots, red where it heats
and blue where it cools). The radiative cooling is a single aggregate curve;
the paper does not plot H3+, Lyman-alpha and recombination separately. What it
says about the split is in words, in section 3.2.1:

> "The stellar XUV heating peak is at a radius of 1.66 Rp
> (p = 2.1 x 10-9 bar), as demonstrated by Figure 9, which shows the prominent
> heating and cooling terms. The peak volume heating rate is 1.8 x 10-8 W m-3.
> At r <~ 2.3 Rp, heating is mostly balanced by H+3 cooling. At higher radii,
> the outflow velocity is significant, and adiabatic cooling due to the
> expansion of the atmosphere primarily balances the stellar heating rate, with
> a smaller contribution from H+3 cooling. Advection makes a small contribution
> to cooling and heat conduction is largely negligible. We note that the rate
> for heat conduction in Figure 9 is diagnostic only."

The terms are the ones in their energy equation (their Equation B3): `q` is the
radiative heating and cooling, `p (1/r^2) d(r^2 w)/dr` the adiabatic term,
`(1/r^2) d(r^2 rho u w)/dr` the advection term and
`(1/r^2) d(r^2 kappa dT/dr)/dr` the conduction term, with `u = cv T`.

---

## 2. Digitization, and how it was checked

The figure page (page 10 of the publisher PDF, `references/Koskinen_2022_ApJ_929_52.pdf`)
was rendered with `pdftoppm -r 600 -png`. The axis frame was located by its own
pixel rows and columns and the calibration taken from the printed ticks: the
left spine carries seven major ticks 307.75 px apart, the bottom of the frame
being 1e-13 W m^-3 and the top 1e-7; the bottom axis runs from `r/R_p = 1` at
column 1376.5 to `r/R_p = 10` at column 4102.5, i.e. 302.889 px per unit.
One pixel is therefore 0.0032 dex in the rate and 0.0033 R_p in radius, and the
curves are 8 to 10 px thick, so a reading is good to about +/- 0.015 dex
(4 percent) away from steep parts.

The curves are drawn in pure red (255,0,0) and pure blue (0,0,255) and the
legend uses black handles, so no masking was needed. The three blue curves
share one color and differ only in dash pattern, so they were separated by
tracking: starting at `r/R_p = 8.91`, where all three are visible and 0.4 to
0.5 dex apart, each was followed column by column leftward, the next point
predicted from a local linear fit to the last 70 accepted points and accepted
within 22 px, with dash gaps of up to 90 columns bridged. The tracking carries
the dashed (radiative) and dash-dot (adiabatic) curves correctly through their
crossing; the overlay of the tracked curves on the rendered panel was inspected
and is the acceptance test.

**Checks against numbers the paper states.**

| check | stated | recovered |
|---|---|---|
| radius of the stellar-heating peak | 1.66 R_p | 1.659 R_p |
| peak volume heating rate | 1.8e-8 W m^-3 | 1.837e-8 W m^-3 |
| where H3+ cooling stops dominating | r <~ 2.3 R_p | radiative and adiabatic curves cross at 2.31 R_p |
| `eps F_XUV/4 = int (H - C) dr` (their Eq. 24) | 0.11 W m^-2 | 0.104 W m^-2 (1.401 R_p to infinity; the unread 1.34-1.40 sliver would add of order 0.015) |

**A stronger check the paper does not state.** Their four plotted terms are the
whole of their Equation (B3) up to conduction, so in a steady state
`stellar heating = radiative + adiabatic + advection`. The digitized curves
satisfy that to **within 4 percent at every radius from 1.10 to 6.7 `r_base`**
(ratios 0.96, 1.00, 0.97, 1.00, 1.01, 1.01, 1.02, 1.01), which tests the axis
calibration, the unit conversion and the assignment of the three blue curves at
once. Below 1.05 `r_base` the curves are near vertical and the closure degrades
to 1.17, so **the first 5 percent in radius is read as order of magnitude
only**, the same caveat as section 4 of `docs/p23_published_profiles.md`.

### 2.1 The digitized figure

Rates in **erg cm^-3 s^-1** (the figure axis is W m^-3; multiply by 10).
`r_base = 1.34 R_p` is their lower boundary.

| r/r_base | 1.05 | 1.10 | 1.15 | 1.20 | 1.25 | 1.30 | 1.40 | 1.50 | 1.75 | 2.00 | 2.50 | 3.00 | 4.00 | 5.00 | 6.00 | 6.70 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r/R_p | 1.407 | 1.474 | 1.541 | 1.608 | 1.675 | 1.742 | 1.876 | 2.010 | 2.345 | 2.680 | 3.350 | 4.020 | 5.360 | 6.700 | 8.040 | 8.978 |
| stellar heating | 1.44e-07 | 1.20e-07 | 1.37e-07 | 1.70e-07 | 1.82e-07 | 1.61e-07 | 9.90e-08 | 6.01e-08 | 2.11e-08 | 1.01e-08 | 3.86e-09 | 1.99e-09 | 7.82e-10 | 3.99e-10 | 2.35e-10 | 1.71e-10 |
| radiative cooling | 6.03e-08 | 6.79e-08 | 8.66e-08 | 1.30e-07 | 1.42e-07 | 1.25e-07 | 7.11e-08 | 3.86e-08 | 8.12e-09 | 2.99e-09 | 1.05e-09 | 5.04e-10 | 1.86e-10 | 8.91e-11 | 5.12e-11 | 3.68e-11 |
| adiabatic | 4.94e-08 | 4.80e-08 | 4.18e-08 | 3.63e-08 | 3.24e-08 | 2.76e-08 | 2.08e-08 | 1.63e-08 | 9.49e-09 | 5.66e-09 | 2.48e-09 | 1.35e-09 | 5.50e-10 | 2.85e-10 | 1.69e-10 | 1.24e-10 |
| advection | 1.29e-08 | 9.27e-09 | 9.37e-09 | 8.56e-09 | 8.20e-09 | 6.72e-09 | 5.10e-09 | 4.61e-09 | 3.08e-09 | 1.27e-09 | 2.92e-10 | 1.24e-10 | 4.70e-11 | 2.76e-11 | 1.63e-11 | 1.08e-11 |

The conduction points, over 1.45 <= r/R_p <= 9.05, lie between 1.5e-13 and
4.4e-11 W m^-3 and never exceed **2.0 percent** of the local stellar heating.

---

## 3. Our side: which state, and how the terms were obtained

**The state.** The converged He/H = 0.0793 rung of `supersonic_molecular_base.md`
section 14 (`p_base = 1e-6 bar`, `q_H2_base = 0.84`, molecular chemistry on,
metals off, `Secondary_ionization: Immediate`), restarted from its own output
and re-converged: `info = 0`, `||R|| = 2.714e-6`, flux spread 2.499e-3,
`log10 Mdot = 10.30` g/s. It is therefore a further chain cycle of that rung in
the sense of section 14.3, and it differs from the section 14 cycle-1 state by
5.1 percent in T at most (at the H2 front, 1.13 `r_base`) and by 0.6 percent
median above 1.10 `r_base` -- inside the reproducibility that section records
for this rung.

**Caveat on the binary.** The run used a scratch build of the repository source
as it stood on 2026-09-03 at 09:45 KST, with the two diagnostic-only changes of
section 9.1 and the A/B knob of section 6 added. The tree was being edited by
other work at the time, so this binary is not bit-identical to the one section
14 used (`md5 b07ed10c...`); the agreement quoted in the previous paragraph is
what bounds the difference.

**The radiative channels** come from `Heating_breakdown.txt` and
`Cooling_breakdown.txt`, physical cells only (rows 3 to 502).

**The hydrodynamic terms are not written by the code** and were evaluated from
the profile as
`adiabatic = p (1/r^2) d(r^2 v)/dr` and
`advection = (1/r^2) d(r^2 rho u v)/dr` with `u = p/((gamma-1) rho)`,
`gamma = 5/3`, by centered differences on the cell centers. `EXHALE_RESIDUAL=1`
prints only the residual norm, not a term dump, so this is the available route.

**How far that is trustworthy.** In a steady state the four terms must satisfy
`heating = radiative + adiabatic + advection`. They do, to **0.4 percent
between 1.25 and 3.0 `r_base` and 4 percent at 1.10**. They do **not** below
1.08: the ratio falls to 0.71 at 1.05 and 0.52 at 1.03. The same failure
appears in the total-energy form and is not removed by taking logarithmic
derivatives, and the cell-center mass flux `rho v r^2` -- which a steady state
makes constant -- varies by 43 percent over 1.005 to 1.10 while being constant
to 0.7 percent above 1.15. **Everything quoted below 1.08 `r_base` is therefore
a term as evaluated and not a closed budget** (section 9.2).

---

## 4. The two budgets side by side

> **PRE-P53, kept as measured (marked 2026-09-03).** The "ours" rows below were
> measured before the composition-dependent caloric EOS and the molecular
> reaction heat existed, and before the section 9.1 fix, so their
> `heat_total` is an under-report and they carry no reaction-heat channel.
> Section 10 has the same table re-measured on a P53 state. The Koskinen rows
> are unchanged.

Rates in erg cm^-3 s^-1. Ours interpolated onto the same `r/r_base` grid. The
Koskinen temperature row is their Figure 7 as digitized in
`docs/p23_published_profiles.md` section 6 (no entry is tabulated there at
1.25).

| r/r_base | 1.05 | 1.10 | 1.15 | 1.20 | 1.25 | 1.30 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|---|
| **T [K], ours** | 472 | 646 | 1224 | 1703 | 2056 | 2312 | 2765 | 2666 | 2116 |
| **T [K], Koskinen** | 1180 | 1566 | 1717 | 1901 | -- | 2267 | 2871 | 4120 | 4742 |
| photoheating, ours | 7.14e-08 | 2.18e-07 | 1.92e-07 | 1.46e-07 | 1.13e-07 | 9.01e-08 | 4.23e-08 | 1.08e-08 | 1.84e-09 |
| &nbsp;&nbsp;H I | 1.79e-08 | 8.91e-08 | 1.03e-07 | 8.98e-08 | 7.55e-08 | 6.34e-08 | 3.29e-08 | 8.75e-09 | 1.46e-09 |
| &nbsp;&nbsp;He I | 4.31e-08 | 1.12e-07 | 7.79e-08 | 4.84e-08 | 3.20e-08 | 2.25e-08 | 8.01e-09 | 1.80e-09 | 3.23e-10 |
| &nbsp;&nbsp;H2 | 1.03e-08 | 3.40e-10 | 2.03e-11 | 7.81e-12 | 4.71e-12 | 3.18e-12 | 8.11e-13 | 7.45e-14 | 6.15e-15 |
| &nbsp;&nbsp;He(2^3S)+H Penning | 2.39e-11 | 1.31e-08 | 1.01e-08 | 7.22e-09 | 5.28e-09 | 3.88e-09 | 1.22e-09 | 1.82e-10 | 3.06e-11 |
| &nbsp;&nbsp;He recombination photons | 1.58e-11 | 3.89e-09 | 1.07e-09 | 4.34e-10 | 2.26e-10 | 1.35e-10 | 3.15e-11 | 3.08e-12 | 1.63e-13 |
| radiative cooling, ours | 5.09e-10 | 1.37e-09 | 1.38e-09 | 1.46e-09 | 1.54e-09 | 1.59e-09 | 1.45e-09 | 5.40e-10 | 8.39e-11 |
| &nbsp;&nbsp;H recombination | 3.79e-10 | 1.08e-09 | 1.05e-09 | 1.08e-09 | 1.10e-09 | 1.09e-09 | 8.92e-10 | 3.22e-10 | 5.10e-11 |
| &nbsp;&nbsp;bremsstrahlung | 9.57e-11 | 2.88e-10 | 3.27e-10 | 3.69e-10 | 4.02e-10 | 4.19e-10 | 3.72e-10 | 1.34e-10 | 1.93e-11 |
| &nbsp;&nbsp;H3+ infrared | 3.48e-11 | 1.08e-14 | 2.03e-16 | 2.62e-16 | 2.55e-16 | 2.07e-16 | 3.83e-17 | 1.31e-18 | 4.78e-20 |
| &nbsp;&nbsp;He I collisional excitation | 1.12e-22 | 1.91e-16 | 5.62e-13 | 9.35e-12 | 3.68e-11 | 7.87e-11 | 1.82e-10 | 8.36e-11 | 1.36e-11 |
| adiabatic, ours | 9.10e-08 | 1.35e-07 | 1.19e-07 | 9.53e-08 | 7.69e-08 | 6.32e-08 | 3.27e-08 | 9.93e-09 | 2.11e-09 |
| advection, ours | 9.18e-09 | 7.30e-08 | 6.79e-08 | 4.86e-08 | 3.48e-08 | 2.55e-08 | 8.40e-09 | 3.84e-10 | -3.37e-10 |
| closure, ours (heating/losses) | 0.71 | 1.04 | 1.02 | 1.01 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| stellar heating, K22 | 1.44e-07 | 1.20e-07 | 1.37e-07 | 1.70e-07 | 1.82e-07 | 1.61e-07 | 6.01e-08 | 1.01e-08 | 1.99e-09 |
| radiative cooling, K22 | 6.03e-08 | 6.79e-08 | 8.66e-08 | 1.30e-07 | 1.42e-07 | 1.25e-07 | 3.86e-08 | 2.99e-09 | 5.04e-10 |
| adiabatic, K22 | 4.94e-08 | 4.80e-08 | 4.18e-08 | 3.63e-08 | 3.24e-08 | 2.76e-08 | 1.63e-08 | 5.66e-09 | 1.35e-09 |
| advection, K22 | 1.29e-08 | 9.27e-09 | 9.37e-09 | 8.56e-09 | 8.20e-09 | 6.72e-09 | 4.61e-09 | 1.27e-09 | 1.24e-10 |
| closure, K22 | 1.17 | 0.96 | 1.00 | 0.97 | 1.00 | 1.01 | 1.01 | 1.02 | 1.01 |

Ratios, ours over theirs:

| r/r_base | 1.05 | 1.10 | 1.15 | 1.20 | 1.25 | 1.30 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|---|
| photoheating | 0.50 | 1.81 | 1.40 | 0.86 | 0.62 | 0.56 | 0.70 | 1.07 | 0.92 |
| radiative cooling | 0.008 | 0.020 | 0.016 | 0.011 | 0.011 | 0.013 | 0.038 | 0.180 | 0.166 |
| adiabatic | 1.84 | 2.81 | 2.86 | 2.62 | 2.38 | 2.29 | 2.01 | 1.76 | 1.56 |
| advection | 0.71 | 7.88 | 7.25 | 5.68 | 4.25 | 3.80 | 1.82 | 0.30 | -2.73 |

(The last entry is negative because our advection term changes sign at
2.1 `r_base`, beyond which it is a heat source rather than a loss.)

The figure is `docs/lower_atmosphere_figs/fig_p23_thermal_budget.png`:
panel (a) the four aggregate terms of both codes plus their conduction points,
panel (b) our channel decomposition against their total radiative cooling.

---

## 5. The four candidates, sized

**(a) Is our photoheating too small at the same radius?** No, not by any
factor that matters. Peak 2.23e-7 erg cm^-3 s^-1 at 1.111 `r_base` against
1.84e-7 at 1.238; the radius of the peak differs by 0.13 `r_base` and the value
by 21 percent, in our favor. The two runs carry the same incident flux
(`F_XUV/4 = 390` erg cm^-2 s^-1 here against 400 there) and different spectra
(power law of index -1 here, the mean solar spectrum of Koskinen et al. 2013
there). The peak sits 0.13 `r_base` deeper here, at a particle density of
4.1e10 cm^-3 against about 9e9 there, so the two models put the same unit
optical depth at very different densities -- which is what a thinner, cooler
column below the front does. That is a description, not an attribution: the
spectra were not compared (section 11). **Nothing in the heating rates supports
a heating deficit.**

**(b) Is our H3+ cooling too large?** The opposite: above 1.05 `r_base` our
H3+ channel is 5e-4 of our own local photoheating at 1.05 and 1e-9 of it at
1.15, and our total radiative cooling -- which by then is H recombination and
bremsstrahlung, not H3+ -- is 1 to 2 percent of theirs. The one place H3+ is
large is the first percent in radius, where it reaches 24 times the local
photoheating at cell 1. Whether the H3+ *rate* is right is a separate
question, tested two ways:

*The formula.* Ours is Miller et al. (2013) Table 5 LTE emission per molecule
times their Table 6 non-LTE departure factor. Theirs, in their own words, is
"optically thin infrared cooling to space by H+3. The per molecule cooling rate
is based on line lists (Neale et al. 1996; Miller et al. 2013) and a correction
factor for non-LTE conditions that is derived from detailed balance
calculations (Koskinen et al. 2009)." Applying **our** formula to **their**
digitized `T`, `n(H2)` and `n(H3+)` gives:

| r/R_p | 1.50 | 1.60 | 2.00 | 2.30 | 2.70 | 3.20 | 4.00 | 5.00 | 6.00 |
|---|---|---|---|---|---|---|---|---|---|
| T [K] | 1611 | 1884 | 2855 | 3517 | 4141 | 4488 | 4740 | 4924 | 5087 |
| our non-LTE factor s | 0.613 | 0.346 | 0.047 | 0.003 | 0.002 | 0.002 | 0.001 | 0.001 | 0.001 |
| our H3+ rate on their state | 1.44e-07 | 1.80e-07 | 1.13e-07 | 1.05e-08 | 7.14e-09 | 3.22e-09 | 1.12e-09 | 3.71e-10 | 1.80e-10 |
| their total radiative cooling | 7.10e-08 | 1.15e-07 | 3.95e-08 | 1.03e-08 | 2.87e-09 | 1.30e-09 | 5.26e-10 | 2.34e-10 | 1.27e-10 |
| ratio | 2.03 | 1.57 | 2.87 | 1.03 | 2.49 | 2.48 | 2.13 | 1.59 | 1.42 |

Our formula reproduces the curve they plot to a factor 1.0 to 2.9 (median 2.0)
on their own composition, which is about what this test can resolve given the
+/- 0.05 dex on their densities, the steepness of the non-LTE factor in
`n(H2)`, and the fact that their curve also contains recombination and
Lyman-alpha. **The 60- to 100-fold difference between the two codes' radiative
cooling is not in the rate coefficient; it is in `n(H3+)`, which differs by
5 orders of magnitude at 1.12 `r_base` and 10 at 4.5 because our H2 is gone.**

*The sensitivity.* Section 6.

**(c) Is our expansion cooling too large?** Yes, by a factor 1.6 to 2.9 at the
same radius, and the advection term by up to 7.9. This is the term that carries
the extra energy our missing coolant leaves in the gas.

**(d) Is conduction the difference?** No. `Conduction:` is not set in any of
these runs, so our conduction is identically zero; theirs is a diagnostic that
never exceeds 2 percent of the local heating anywhere between 1.45 and 9 R_p,
and they say so themselves ("heat conduction is largely negligible").

---

## 6. The A/B: scaling the H3+ coolant

The most direct test of candidate (b), and of `TO_BE_DONE.md` item (X), is to
multiply the H3+ cooling rate and re-converge. An environment knob
`EXHALE_H3P_SCALE` was added **to the scratch source tree only** (it multiplies
the return of `h3p_cooling_rate`, which both the direct and the infrared-field
paths go through) and four states were converged from the same starting state
with the same binary, each reaching `info = 0` on both gates.

**Scratch source changes, and which of them is meant for the repository.** This
knob is a measurement instrument and is not; the `write_heat_breakdown_eq` fix
of section 9.1 is, and its patch and verification are recorded there
(`.../scratchpad/p23heat/p52_heating_breakdown_fix.diff`, applied with `-p1`
from `EXHALE_v1.00/`). The `Photorates_diag.txt` dump used by section 8.1 is
also scratch-only. Nothing in this campaign was written to the repository
`src/`.

| `EXHALE_H3P_SCALE` | 0.0 | 0.5 | 1.0 | 3.0 |
|---|---|---|---|---|
| T at 1.001 `r_base` [K] | 1058 | 1029 | 1002 | 923 |
| T at 1.003 | 893 | 867 | 846 | 794 |
| T at 1.01 | 696 | 690 | 685 | 673 |
| T at 1.03 | 659 | 657 | 656 | 652 |
| T at 1.05 | 475 | 473 | 472 | 469 |
| T at 1.10 | 642 | 644 | 646 | 650 |
| T at 1.15 | 1219 | 1222 | 1224 | 1228 |
| T at 1.50 | 2764 | 2765 | 2765 | 2767 |
| T at 3.00 | 2116 | 2116 | 2116 | 2116 |
| minimum T [K] (at `r/r_base`) | 415 (1.0647) | 415 (1.0639) | 415 (1.0639) | 417 (1.0630) |
| H2 front `f = 0.5` [`r_base`] | 1.0431 | 1.0431 | 1.0425 | 1.0419 |
| `log10 Mdot` | 10.30 | 10.30 | 10.30 | 10.29 |

**Removing the H3+ coolant entirely buys 56 K at 1.001 `r_base`, 11 K at 1.01,
3 K at 1.03 and nothing at all from 1.10 outward, and it does not move the
temperature minimum.** On this state the deficit against Koskinen's Figure 7 is
359 K at `(r - r_base)/r_base = 1e-2` and 206 K at 3e-2 (item (X)'s own row,
on the different state noted below, gives 435 and 468 K); removing the H3+
coolant entirely recovers 11 K and 3 K of that, and a factor-2 error in the rate
would account for 6 K and 2 K.

*One caveat on which state item (X) tabulates.* Its "EXHALE, 1 microbar base"
row is the `Base BC: pressure 1.0` run of `docs/p44_base_sawtooth.md`
section 9.3 (`r_front = 1.03551`), not the ladder rung used here
(`p_base 1.000e-06` in `base.inp`, `r_front = 1.0425`). The two agree at the
top of the layer (1030 against 1002 K at 1e-3, 849 against 846 at 3e-3) and put
their trough at a different radius (394 K near 3e-2 there against 415 K at
6.5e-2 here), so the trough is the same feature at the same depth in a slightly
different place. The A/B above was run on the ladder rung; nothing here re-runs
it on the p44 configuration. **The suspicion recorded in item (X) -- "an H3+
rate or an escape treatment that is too efficient would show exactly this" --
is not supported.** What the trough does respond to is the expansion term: at
1.02 to 1.05 the adiabatic term is 8.6e-8 to 9.1e-8 erg cm^-3 s^-1 against a
photoheating of 4.9e-8 to 7.2e-8 and a radiative cooling of 1.5e-8 falling to
5e-10.

---

## 7. What the extra energy does, and why the temperature is still lower

Two ledgers, both measured.

**The energy ledger.** Integrating the net radiative heating radially,
`int (H - C) dr`, gives **0.226 W m^-2** here against the **0.11 W m^-2** their
Equation (24) states, from an incident `F_XUV/4` of 390 against 400
erg cm^-2 s^-1. In their notation that is an efficiency `eps = 0.58` here
against `eps = 0.28` there. The ratio of the two net inputs is 2.05, and the
ratio of the two mass fluxes, both read as `rho v r^2` per steradian, is
**3.14e9 / 1.5e9 = 2.09**. Within
the reading error of the digitization and of their quoted rate, **the energy
per unit mass driven off the planet is the same in the two models, and the
factor 2 in the escape rate is exactly the factor 2 in how much of the absorbed
flux survives radiative cooling.**

**The composition ledger.** The temperature at fixed thermal energy per unit
mass scales with the mean molecular weight. Theirs is 2.24 mH at 1.12 `r_base`
falling to 1.32 at 4.48; ours is 1.22 falling to 0.68, because our gas is
atomic and increasingly ionized where theirs stays 40 to 95 percent H2.

| r/r_base | 1.119 | 1.194 | 1.306 | 1.493 | 1.716 | 2.015 | 2.388 | 2.985 | 3.731 | 4.478 |
|---|---|---|---|---|---|---|---|---|---|---|
| T, Koskinen [K] | 1611 | 1884 | 2286 | 2855 | 3517 | 4141 | 4488 | 4740 | 4924 | 5087 |
| T, ours [K] | 869 | 1653 | 2338 | 2759 | 2802 | 2657 | 2434 | 2123 | 1826 | 1601 |
| mu, Koskinen [mH] | 2.24 | 2.17 | 2.00 | 1.75 | 1.65 | 1.57 | 1.49 | 1.43 | 1.38 | 1.32 |
| mu, ours [mH] | 1.22 | 1.20 | 1.16 | 1.07 | 0.98 | 0.89 | 0.82 | 0.75 | 0.71 | 0.68 |
| (T/mu) ours / (T/mu) theirs | 0.99 | 1.58 | 1.76 | 1.57 | 1.35 | 1.14 | 0.99 | 0.85 | 0.72 | 0.61 |
| our T if our gas had their mu [K] | 1603 | 2980 | 4027 | 4489 | 4740 | 4705 | 4447 | 4038 | 3558 | 3091 |

**Read the last two rows.** Over 1.19 to 2.02 `r_base` our thermal energy per
unit mass is 1.1 to 1.8 times *theirs*; the temperature is lower only because
it is spread over particles half as heavy. Give our gas their mean molecular
weight and our temperature would be 1.1 to 1.8 times theirs over that interval.
Beyond about 2.4 `r_base` the specific thermal energy falls below theirs as
well (0.85 at 3.0, 0.61 at 4.5), so the outer part of the gap -- the 2616 K at
3.00 `r_base` that section 14.5 records -- is roughly half composition and half
energy. **The tentative conclusion is that the 1.15-2.0 `r_base` temperature
gap is a consequence of the H2 abundance and not of the thermal budget, which
places it under item (P43) (H2 transport) rather than under a cooling rate.**

---

## 8. The two anomalies of section 14.4

### 8.1 He/H = 30: three decades of H2 lost inside the first cell

The handoff states `x2 = 0.985448` at the ghost; cell 1 carries
`x2 = 4.56e-4`, and `n(H3+)` never exceeds 0.36 cm^-3 anywhere. The chemistry
of cell 1, computed from the converged state (`heh30c_r1` re-converged with the
photoionization-rate dump of section 9.1):

| | ghost (row 1) | cell 1 | cell 2 | cell 3 |
|---|---|---|---|---|
| r [R_p] | 1.00000 | 1.00019 | 1.00039 | 1.00058 |
| T [K] | 1140.0 | 1120.2 | 1115.5 | 1111.2 |
| n(H2) [cm^-3] | 1.03e11 | 4.73e07 | 4.54e07 | 4.39e07 |
| n(He II) [cm^-3] | 4.75e03 | 4.00e06 | 4.13e06 | 4.25e06 |
| n(H I) [cm^-3] | 2.99e09 | 2.07e11 | 2.05e11 | 2.03e11 |
| n_e [cm^-3] | 4.21e07 | 5.33e07 | 5.34e07 | 5.35e07 |

The H2 destruction budget at cell 1, in cm^-3 s^-1:

| channel | rate | share | `n(H2)/rate` [s] |
|---|---|---|---|
| R17 He+ + H2 -> H+ + H + He | 1.166e+03 | 93.5 % | 4.05e4 |
| R20 He+ + H2 -> HeH+ + H | 7.94e+01 | 6.4 % | 5.95e5 |
| R23 He+ + H2 -> H2+ + He | 1.36e+00 | 0.11 % | 3.47e7 |
| R13 H+ + H2 + M -> H3+ + M | 4.80e-01 | 0.04 % | 9.86e7 |
| R8 H2+ + H2 -> H3+ + H | 7.46e-02 | 0.01 % | 6.33e8 |
| photoionization (`P_H2 = 1.26e-9` s^-1) | 5.95e-02 | 0.005 % | 7.95e8 |
| R18 HeH+ + H2 -> H3+ + He | 3.91e-02 | 0.003 % | 1.21e9 |
| R10 H+ + H2 -> H2+ + H | 7.54e-03 | 0.001 % | 6.27e9 |
| R12 H2 + M -> H + H + M (thermal) | 1.34e-08 | 0 % | 3.54e15 |
| R14 e + H2 -> H + H + e | 2.35e-13 | 0 % | 2.01e20 |
| **total** | **1.248e+03** | | **3.79e4** |

**Helium ions destroy 99.9 percent of the H2, and R17 alone 93.5 percent.**
The photoionization channel is four orders below R17. The proximate cause is
the He II density, which is **843 times larger in cell 1 than in the ghost**
(4.00e6 against 4.75e3) at essentially the same temperature: the handoff writes
a composition that the interior chemistry does not reproduce, and the interior
chemistry, given 30 helium atoms per hydrogen and `P_HeI = 5.65e-10` s^-1,
puts helium into He II and He II then eats the H2. This is hypothesis H2 of
`supersonic_molecular_base.md` section 5 -- the pinned base state and the
interior chemistry disagree about what the base gas is -- restated with the
reaction that carries it.

**The timescales say the local-equilibrium closure is the wrong one here.** The
H2 chemical lifetime at cell 1 is 3.8e4 s, and the residence time in that cell
(`dr/|v|`, with `v = -323` cm s^-1) is 2.0e3 s, nineteen times shorter. A
carrier that were transported rather than held at local equilibrium would
therefore cross the first cell substantially undestroyed. **The three-decade
drop is a property of the local-equilibrium closure, not a statement about how
much H2 survives.**

**The grid test says the drop is not a resolution artifact.** A cold start of
the same rung with the base cells halved and the count doubled
(`Base grid [dr,cells]: 1.0e-4 100`, `Grid cells: 550`) converged on its own --
WENO3 at step 32,197, the staged secondary ionization at 85,697, the JFNK at
87,697, `info = 0`, `||R|| = 5.343e-6`, flux spread 4.936e-3,
`log10 Mdot = 10.30`.

| | 500 cells, base `2e-4 x 50` | 550 cells, base `1e-4 x 100` |
|---|---|---|
| `x2` at the ghost | 0.985448 | 0.985448 |
| `x2` at cell 1 (`r`) | 4.563e-4 (1.000193) | **4.593e-4** (1.000096) |
| `x2` at cell 2 / cell 5 | 4.432e-4 / 4.120e-4 | 4.515e-4 / 4.347e-4 |
| n(He II), ghost -> cell 1 | 4.75e3 -> 4.00e6 (x842) | 4.77e3 -> 3.97e6 (**x833**) |
| n(H2) at cell 1 [cm^-3] | 4.729e7 | 4.771e7 |
| peak n(H3+) [cm^-3] (at `r`) | 0.36393 (1.00019) | 0.36806 (1.00010) |
| `v` at cell 1 [cm/s] | -324.3 | **-153.2** |
| T max [K] | 7104 | 7164 |
| ionization front [R_p] | 1.3720 | 1.3680 |
| `log10 Mdot` | 10.29 | 10.30 |

**Halving the base cell changes `x2` at cell 1 by 0.6 percent and the He II
jump by 1 percent, while it halves the velocity artifact** (-324 to -153 cm/s,
the first-order-in-`dr` behavior `docs/p44_base_sawtooth.md` section 7.2
describes). The three-decade drop is therefore a property of the chemistry the
first cell is asked to solve, not of how finely the first cell is drawn -- which
is what a local-equilibrium composition, having no memory of `dr`, has to do.

### 8.2 He/H = 0.3: a cold first cell, and an alternating temperature

Re-converging that rung (`heh0p3c_r1`) reproduces the anomaly section 14.4
records: cell 1 at 309 K (312 K there), `n(H2)` at cell 1 3.7 times the ghost
value (3.6 there), and T alternating over the first six cells. The term-by-term
state, rows 0 and 1 being the ghosts:

| row | r [R_p] | T [K] | v [cm/s] | n(H2) | n_e | heating | cooling (= H3+ IR) | `p div v` | advection |
|---|---|---|---|---|---|---|---|---|---|
| 1 (ghost) | 1.00000 | 1140 | 6.60e-01 | 3.88e12 | 1.37e06 | 4.45e-08 | 6.55e-07 | -2.69e-04 | -3.98e-04 |
| 2 (cell 1) | 1.00019 | **309** | -3.55e+02 | 1.42e13 | 4.71e05 | 1.62e-07 | 3.81e-10 | 2.09e-04 | 3.06e-04 |
| 3 | 1.00039 | 1126 | 2.80e+02 | 3.49e12 | 1.83e06 | 4.31e-08 | 4.58e-07 | 2.36e-04 | 3.62e-04 |
| 4 | 1.00058 | 1302 | -3.22e+01 | 2.91e12 | 1.87e06 | 3.68e-08 | 7.83e-07 | -1.72e-04 | -2.61e-04 |
| 5 | 1.00077 | 1428 | 4.23e+01 | 2.57e12 | 1.28e06 | 3.30e-08 | 1.62e-06 | 1.81e-05 | 2.75e-05 |
| 6 | 1.00097 | 1176 | -6.95e+00 | 3.23e12 | 2.01e06 | 4.11e-08 | 4.92e-07 | -3.25e-05 | -4.90e-05 |
| 7 | 1.00116 | 1040 | -3.28e+00 | 3.68e12 | 1.77e06 | 4.69e-08 | 3.44e-07 | 3.37e-06 | 5.15e-06 |

**The source terms cannot be what makes cell 1 cold.** It has the *largest*
heating of the seven rows (1.62e-7, four times its neighbors, because its H2
column is 3.7 times theirs and the H2 photoheating channel scales with it) and
the *smallest* cooling (3.81e-10, three orders below its neighbors, because the
Miller et al. (2013) LTE emission of one H3+ molecule falls by a factor 6.7e3
between 1100 K and 309 K). If the local source balance set the temperature, cell 1
would be the hottest cell, not the coldest.

**What alternates is the velocity, and the temperature follows it.** `v` runs
-355, +280, -32, +42, -7.0, -3.3 cm s^-1 over the first six cells -- the
odd-even pattern -- and the two hydrodynamic terms evaluated from the profile
alternate in sign at 1e-4 to 3e-4 erg cm^-3 s^-1, three orders above the source
terms. The gate rung (He/H = 0.0793) carries the same single-cell velocity
artifact at cell 1 (`v = -234` cm s^-1) but its T is monotone
(1098, 1074, 1050, 1027, 1006 K), so what distinguishes He/H = 0.3 is that the
oscillation has reached the temperature and not that its energy balance is
different in kind. **The likely reading is that this is the first-order-in-`dr`
base discretization error of `docs/p44_base_sawtooth.md` section 7.2, coupled
into T through `T = p / n_part` -- the same route section 138 identifies -- and
not a thermal effect.** The alternating terms in the last two columns are
centered differences of an odd-even mode and are quoted only to show the
magnitude; they do not close (section 9.2).

**It is not the section 138 limiter firing.** The ladder ledger of
`supersonic_molecular_base.md` section 14.2 records zero energy-floor hits and
zero positivity-limiter activations -- face states, ghost cells and face fluxes
alike -- on all six rungs, He/H = 0.3 included, and zero trial states refused
for an uncertified composition. Whatever sustains the mode, the guards that
section 138 installed are not switching on and off under it.

**The grid test removes it.** The same refinement
(`Base grid [dr,cells]: 1.0e-4 100`, `Grid cells: 550`), cold start, converged
on its own: WENO3 at step 33,772, staged secondary ionization at 117,856, the
JFNK at 119,856, four solves all `info = 0`, final `||R|| = 5.356e-6`, flux
spread 4.667e-3, `log10 Mdot = 10.26`.

| | 500 cells, base `2e-4 x 50` | 550 cells, base `1e-4 x 100` |
|---|---|---|
| T over the first six cells [K] | 316, 1136, 1318, 1425, 1172, 1039 | **1104, 1095, 1086, 1078, 1070, 1061** |
| `v` over the same cells [cm/s] | -346, +226, -26.5, +45.4, -8.2, -1.4 | **-125, +19.7, +6.2, +9.6, +8.1, +8.6** |
| n(H2) at cell 1, over the ghost | 3.58 | **0.957** |
| minimum T beyond cell 6 [K] (at `r`) | 471 (1.0590) | 341 (1.0418) |
| peak n(H3+) [cm^-3] (at `r`) | 4.037e5 (1.00019) | **7.535e4** (1.00221) |
| `f(H2) = 0.5` front [R_p] | 1.03339 | 1.03105 |
| T max [K] | 3684 | 3728 |
| ionization front [R_p] | 2.0788 | 2.0296 |
| `log10 Mdot` | 10.29 | 10.26 |

**The cold cell and the alternation are gone.** T is monotone through the first
ten cells, `v` is positive from cell 2 onward with only the single-cell artifact
left at cell 1 -- and that artifact halves with the cell, -346 to -125 cm/s,
the same first-order-in-`dr` scaling the He/H = 30 rung shows. n(H2) at cell 1
is no longer 3.6 times the ghost value but 0.96 of it. **This is a
grid-dependent artifact, and the reading above is confirmed.**

Two consequences worth recording, because they are not small. The rung's
molecular-layer numbers moved with the artifact: **peak n(H3+) falls by a factor
5.4** (4.04e5 to 7.54e4 cm^-3) and moves out from 1.0002 to 1.0022 `R_p`, and
the layer temperature minimum falls from 471 K at 1.0590 to 341 K at 1.0418.
`log10 Mdot` moves by 0.03 dex and the ionization front by 0.05 R_p. So on this
rung the coarse base grid was not only carrying a cosmetic first-cell artifact
-- it was setting the H3+ peak.

---

## 9. Two defects found on the way

### 9.1 `Heating_breakdown.txt` under-reports the molecular photoheating -- SINCE FIXED

> **Fixed in the tree (marked 2026-09-03).** `write_heat_breakdown_eq` now passes `f_vibq` and `e_vibq` to `PH_heat_HHe`. On the section 10 state the dumped `heat_total` agrees with the solver's own `heat` column to 1.6e-6, against the 1.4e-3 measured below. The measurement that follows is kept as the record of the defect.

`write_heat_breakdown_eq` (`src/modules/radiation/util_ion_eq.f90`) recomputes
the photoheating with the same `PH_heat_HHe` the solver uses, but on the
molecular branch it omits the two optional arguments `f_vibq` and `e_vibq`
that `ionization_equilibrium` supplies. The routine's own comment says what
that costs: "Absent => the two vibrational heat channels are off". The result
is that the dumped `heat_total` and its `heat_H2` column are too small wherever
H2 is abundant. The clean measure is within one run: the dumped total against
the `heat` column of that same run's `Hydro_ioniz.txt`, which the solver wrote.

| r/r_base | 1.001 | 1.005 | 1.01 | 1.02 | 1.03 | 1.05 | 1.07 | 1.10 | 1.15 | 1.30 |
|---|---|---|---|---|---|---|---|---|---|---|
| as written | 0.514 | 0.514 | 0.520 | 0.564 | 0.672 | 0.906 | 0.993 | 1.000 | 1.000 | 1.001 |
| with `f_vibq`/`e_vibq` passed | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.001 | 1.000 | 1.000 | 1.000 | 1.000 |

So the diagnostic dropped **49 percent** of the heating in the first cells and
is exact once the two arguments are passed.

**Handoff, and how it was resolved.** The repository `src/` was off limits for
this work, so the change first sat in the campaign scratch as a patch for
whoever carried it in. **It has since been carried in**, and the record of that
is `Update_EXHALE_stage1.md` section 140.5 ("A second diagnostic-only defect travels
with it"), which measures the dump's photoheating total rising by up to a factor
1.97 at `r = 1.0002 R_p` and by 0.49 per cent volume-integrated on
`mol_base_handoff`, with `Hydro_ioniz.txt` and `Ion_species.txt` untouched. In
the tree as it now stands `write_heat_breakdown_eq`
(`src/modules/radiation/util_ion_eq.f90`) declares `P_H2_di` and passes
`f_vibq = f_vib_quench, e_vibq = e_vib_bound`, so the paragraphs below are the
record of the defect and of the patch, not an open handoff.

**Identifier.** This item was labeled "P52" while it was open. That label is
taken: `TO_BE_DONE.md` uses **(P52)** for a different item -- "a run could not
say why it stopped, and 'accepted' had two definitions", closed under
`Update_EXHALE_stage1.md` section 140 -- so referring to the heating dump by it points
a reader at the wrong entry. Call it the **heating-dump fix
(`Update_EXHALE_stage1.md` section 140.5)**. The scratch paths below keep the old
`p52_` and `_p52` spellings because those are the file names on disk:

```
.../scratchpad/p23heat/p52_heating_breakdown_fix.diff      # apply with -p1 from EXHALE_v1.00/
.../scratchpad/p23heat/util_ion_eq_full.diff               # the same, plus the two scratch-only changes below
.../scratchpad/p23heat/tree_p52/                           # a clean copy of src/ + Makefile with only the fix applied
.../scratchpad/p23heat/verify_p52/                         # the run that verified it
```

Three hunks, all inside `write_heat_breakdown_eq` and its `use` statement: add
`h2_energy_per_bound_fluorescence_eV` to the `h2_vibrational_relaxation` import,
declare `P_H2_di`, `f_vib_quench` and `e_vib_bound`, and pass
`P_H2_di, f_vibq = f_vib_quench, e_vibq = e_vib_bound` in the molecular
`PH_heat_HHe` call. **Verified against the repository source as it stood on
2026-09-03 at 10:56 KST**: `patch -p1 --dry-run` clean, full rebuild from a
clean copy of `src/` clean, and a restart of the gate rung on that binary
reproduces the untouched state exactly (`info = 0`, `||R|| = 2.714e-06`, flux
spread 2.499e-03, `log10 Mdot = 10.30`, the same three numbers the unpatched
binary gives) while the dumped total now agrees with the solver's `heat` column
to 1.000 at every radius. So it touches no solver path and no golden should
move; the byte-identity check on the regression matrix was the confirmation to
run before carrying it in. **That check was run**: `Update_EXHALE_stage1.md` section
140.7 reports all ten regression cases byte-identical in `Hydro_ioniz.txt` and
`Ion_species.txt` against both the control build and `backup/regression/golden/`.

**Two further changes in `.../scratchpad/p23heat/tree/` are scratch-only and
must NOT be carried in**: the `Photorates_diag.txt` dump (section 8.1 needs the
photoionization rates) and the `EXHALE_H3P_SCALE` knob in
`lower_atmosphere/h3p_cooling.f90` (section 6). `p52_heating_breakdown_fix.diff`
already excludes both.

The tables of section 4 above 1.10 `r_base` are unaffected by the defect; the
molecular-layer numbers of sections 6 and 8 use the corrected dump.

A second, smaller item in the same routine, and it was carried in with the
first: `P_H2_di` was likewise not passed,
which is harmless because it is an output, but the call is no longer the same
call the solver makes and that is the thing worth keeping true.

### 9.2 The converged state's own energy budget does not close below 1.08 `r_base` -- SINCE CLOSED

> **Closed (marked 2026-09-03), and it took two changes.** Evaluating the advection term with the caloric `gamma_eff` the code now integrates, on a state finished under the P54 face-flux measure, the three ratios below become 0.99, 1.01 and 1.02. Neither change alone does it. Section 10 has the side-by-side. The diagnosis recorded here -- that the cell-center products are not steady in the layer -- is the half that survives, and section 14.7 of `docs/supersonic_molecular_base.md` measures it directly.

Evaluated from the written profile, `heating = radiative + adiabatic +
advection` holds to 0.4 percent from 1.25 to 3.0 `r_base` and to 4 percent at
1.10, and fails below: 0.71 at 1.05, 0.52 at 1.03, 0.56 at 1.02. The same
failure appears in the total-energy form
`div[r^2 rho v (u + p/rho + v^2/2)]/r^2 + rho v GM/r^2 = H - C` (ratios 2.1 at
1.02, 2.0 at 1.03, 1.4 at 1.05, 0.99 at 1.08), and in the integral form over
shells (ratio 2.19 over [1.005, 1.05], 1.04 over [1.05, 1.10], 1.005 over
[1.10, 1.30], 1.048 over the whole domain). It is not removed by using
logarithmic derivatives.

The diagnosis is not settled. What is measured is that the cell-center mass
flux `rho v r^2`, which a steady state makes constant, **varies by 43 percent
over 1.005 to 1.10** and is constant to 0.7 percent above 1.15, while the
solver's own residual on the same state is `||R|| = 2.7e-6` and its flux-spread
gate reads 2.5e-3. Two readings are consistent with that and are not separated
here: the finite-volume state is steady in its own discrete sense while the
cell-center products used above are not (the profile changes by a factor 2 per
cell over 1.02 to 1.07), or the gate that certifies the state does not see the
base layer. **Either way, no budget quoted below 1.08 `r_base` in this file is a
closed budget, and the flux-spread gate's coverage of the base layer is worth a
separate look.**

### 9.3 The adiabatic index in the molecular layer -- SINCE ACTED ON (P53)

**Superseded 2026-09-03.** When this was written, `parameters.f90` fixed
`gamma = 5/3` for the whole domain, including the layer that is 80 percent H2
by volume at 1000 K, where the rotational degrees of freedom make the effective
value nearer 7/5. It changed the specific heat by a factor 5/3, hence the
advection term and the sound speed, though not `p div v`. It was recorded here
as a stated approximation with no note of its range at the definition, was not
tested, and was not claimed as a contributor to anything above.

P53 replaced it: `src/modules/states/caloric_eos.f90` carries a
composition-dependent caloric EOS in which H2 gets its own rovibrational heat
capacity from the 302-level Roueff et al. (2019) ladder, so
`gamma_eff = 1 + (n_tot + n_e)/C_V_over_k` is a function of composition and
temperature. On the state of section 10 it runs 1.409 at the base to 5/3 by
1.10 `r_base`. It carries no input key: `caloric_state_from_composition`
returns on the legacy constant-gamma branch only when the run has no molecular
chemistry, so for a molecular run with H2 in any cell it is unconditional.
**The approximation mattered**: evaluating the advection term
with `gamma_eff` instead of 5/3 is one of the two changes that close the
section 9.2 budget failure. See section 10.

---

## 10. After P53: the channel table re-measured

> **RE-MEASURED UNDER BLOCK I; SEE SECTION 13 (marked 2026-09-04).** The state
> below predates block I (sections 148-155 of `docs/Update_EXHALE_stage1.md`). Section
> 13 re-finishes the same gate rung under the current tree binary and
> re-measures every table of this section. **Every channel and every ratio
> below reproduces there**, the largest change being the H2 photoheating at
> 1.05 `r_base` (7.84e-8 here, 7.29e-8 there, -7 percent) under section 151's
> measured cross sections; the temperature deficits against Koskinen move by at
> most 3 K. The one row that does NOT carry over is the `rho v r^2` caveat at
> the end of "The state" paragraph: block I's base condition removes the
> base-adjacent part of that non-flatness. Read this section for the P53
> attribution and section 13 for the numbers under the current tree.

> **PROVISIONAL, under the current residual norm (marked 2026-09-03).** The
> state below was accepted by the two gates of section 133 of
> `docs/Update_EXHALE_stage1.md` as that norm stands in the tree at 17:38 KST.
> Whether the norm should be volume weighted is under review; if it changes,
> this state is re-finished and the table is re-measured. The runs and the
> scripts are kept for that.

**Judgment first. The molecular reaction heat is the largest single heating
channel of the molecular layer -- 5.26e-8 erg cm^-3 s^-1 at 1.05 `r_base`,
55 percent of the photoheating there and 71 percent of the total at its peak --
and the caloric EOS closes the base-layer energy budget that section 9.2
recorded as failing. The two channels are now separated (section 10.3): the
reaction heat carries the whole +0.05 dex in `Mdot` and the whole +169 K in the
layer temperature minimum, and the caloric EOS on its own leaves the rate where
it was. Neither changes the section-5 verdict: our radiative
cooling is still 1.0 to 2.0 percent of Koskinen's over the front, our
photoheating is still the same as theirs to a factor 1.6 either way, and our
adiabatic and advection terms are still 1.8 to 3.0 and 2.3 to 6.9 times
theirs. The attribution of the temperature gap is unchanged.**

**The energy budget now closes below 1.08 `r_base`, and it took two changes to
do it.** Section 9.2 measured `heating = radiative + adiabatic + advection` at
0.56, 0.52 and 0.71 at 1.02, 1.03 and 1.05 `r_base`. Evaluating the advection
term with the caloric `gamma_eff` of the cell instead of a constant 5/3 --
which is what the code itself now integrates -- and reading it off a state
finished on the P54 face-flux measure, the same three ratios are **0.99, 1.01
and 1.02**. Neither change alone does it: on the section-4 state `gamma_eff`
alone gives 0.72, 0.71, 0.93, and on the P53 state before its re-finish it
overshoots to 1.09, 1.42, 0.75. The reading this supports is that section 9.2
was two errors at once -- the specific heat used to evaluate the term, and a
layer that was not steady in the state the term was read from.

| closure `heating/(radiative + adiabatic + advection)` | 1.02 | 1.03 | 1.05 | 1.08 | 1.10 | 1.15 | 1.30 | 1.50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| section 4 state, `gamma = 5/3` (section 9.2) | 0.56 | 0.52 | 0.71 | 1.01 | 1.04 | 1.02 | 1.00 | 1.00 |
| section 4 state, `gamma_eff` | 0.72 | 0.71 | 0.93 | 1.02 | 1.04 | 1.02 | 1.00 | 1.00 |
| P53 cold-converged, `gamma_eff` | 1.09 | 1.42 | 0.75 | 0.93 | 0.97 | 1.03 | 1.00 | 1.00 |
| **P54 re-finished, `gamma_eff`** | **0.99** | **1.01** | **1.02** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** |
| P54 re-finished, `gamma = 5/3` | 0.89 | 0.82 | 0.71 | 0.97 | 1.00 | 1.00 | 1.00 | 1.00 |

**Section 9.1 is fixed in the tree.** `write_heat_breakdown_eq` now passes
`f_vibq` and `e_vibq` to `PH_heat_HHe`, and on this state the dumped
`heat_total` agrees with the solver's own `heat` column to **1.6e-6** over
`r > 1.001`, against 1.4e-3 on the section-4 state. The channel table below is
therefore the run's own heating, not an under-report of it.

**The state.** The converged He/H = 0.0793 rung of
`supersonic_molecular_base.md` section 14.8 (`p_base = 1e-6 bar`, molecular
chemistry on, metals off, `Secondary_ionization: Immediate`): cold-started
under the P53 binary (md5 `ea9e7492`, block-G goldens 10/10) and re-finished
under the P54 binary (md5 `ec1f367c`, 10/10 identical). `info = 0`,
`||R|| = 9.198e-6`, gate flux spread 2.501e-4, `log10 Mdot = 10.34` g/s.
Section 14.8's caveat applies here too: both gates are gates on the wind, and
over `r >= 1.03` this state's `rho v r^2` spread is 3.4e-3, not flat.

The channels come from `Heating_breakdown.txt` and `Cooling_breakdown.txt`,
physical cells only. The two hydrodynamic terms are still evaluated from the
profile as section 3 sets out, and are given with both choices of gamma.

#### 10.1 The budget, ours beside theirs

Rates in erg cm^-3 s^-1. Ours interpolated onto the same `r/r_base` grid.

| r/r_base | 1.05 | 1.1 | 1.15 | 1.2 | 1.25 | 1.3 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| **T [K], ours** | 633 | 592 | 983 | 1448 | 1826 | 2115 | 2683 | 2687 | 2164 |
| T [K], Koskinen | 1180 | 1566 | 1717 | 1901 | 2084 | 2267 | 2871 | 4120 | 4742 |
| f(H2), ours | 0.393 | 0.00312 | 0.000266 | 7.88e-05 | 5.01e-05 | 3.96e-05 | 1.98e-05 | 5.19e-06 | 1.47e-06 |
| gamma_eff, ours | 1.5782 | 1.6660 | 1.6666 | 1.6666 | 1.6667 | 1.6667 | 1.6667 | 1.6667 | 1.6667 |
| photoheating, ours | 9.61e-08 | 1.53e-07 | 1.92e-07 | 1.61e-07 | 1.28e-07 | 1.03e-07 | 4.88e-08 | 1.26e-08 | 2.14e-09 |
| &nbsp;&nbsp;H I | 9.09e-09 | 5.45e-08 | 8.91e-08 | 8.98e-08 | 7.95e-08 | 6.84e-08 | 3.70e-08 | 1.01e-08 | 1.69e-09 |
| &nbsp;&nbsp;He I | 2.65e-08 | 7.62e-08 | 8.84e-08 | 6.06e-08 | 4.05e-08 | 2.82e-08 | 9.66e-09 | 2.09e-09 | 3.72e-10 |
| &nbsp;&nbsp;H2 | 7.84e-09 | 2.30e-10 | 3.48e-11 | 1.06e-11 | 6.00e-12 | 4.07e-12 | 1.12e-12 | 9.79e-14 | 7.49e-15 |
| &nbsp;&nbsp;He(2^3S)+H Penning | 4.63e-12 | 7.57e-09 | 1.06e-08 | 8.11e-09 | 6.16e-09 | 4.68e-09 | 1.57e-09 | 2.21e-10 | 3.48e-11 |
| &nbsp;&nbsp;He recombination photons | 3.30e-12 | 3.67e-09 | 2.01e-09 | 7.59e-10 | 3.70e-10 | 2.10e-10 | 4.49e-11 | 4.20e-12 | 2.19e-13 |
| &nbsp;&nbsp;H2 Lyman-Werner | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;FUV photolysis | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;molecular reaction heat | 5.26e-08 | 1.04e-08 | 2.06e-09 | 1.43e-09 | 1.11e-09 | 8.69e-10 | 3.14e-10 | 4.63e-11 | 7.33e-12 |
| radiative cooling, ours | 6.39e-10 | 1.17e-09 | 1.29e-09 | 1.31e-09 | 1.38e-09 | 1.44e-09 | 1.44e-09 | 6.09e-10 | 9.87e-11 |
| &nbsp;&nbsp;H recombination | 2.75e-10 | 9.28e-10 | 9.99e-10 | 9.83e-10 | 1.01e-09 | 1.02e-09 | 9.07e-10 | 3.64e-10 | 5.99e-11 |
| &nbsp;&nbsp;bremsstrahlung | 7.48e-11 | 2.44e-10 | 2.94e-10 | 3.20e-10 | 3.53e-10 | 3.78e-10 | 3.72e-10 | 1.52e-10 | 2.29e-11 |
| &nbsp;&nbsp;H3+ infrared | 2.90e-10 | 5.03e-14 | 4.52e-16 | 2.20e-16 | 2.53e-16 | 2.43e-16 | 6.13e-17 | 1.91e-18 | 6.20e-20 |
| &nbsp;&nbsp;H2 infrared | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;He I collisional excitation | 2.33e-21 | 1.35e-17 | 4.68e-14 | 2.22e-12 | 1.46e-11 | 4.20e-11 | 1.62e-10 | 9.32e-11 | 1.60e-11 |
| adiabatic, ours (gamma = 5/3) | 1.26e-07 | 1.20e-07 | 1.26e-07 | 1.06e-07 | 8.62e-08 | 7.11e-08 | 3.68e-08 | 1.12e-08 | 2.39e-09 |
| advection, ours (gamma = 5/3) | 8.94e-09 | 3.17e-08 | 6.50e-08 | 5.40e-08 | 4.04e-08 | 3.02e-08 | 1.06e-08 | 7.72e-10 | -3.32e-10 |
| advection, ours (gamma_eff) | -3.24e-08 | 3.15e-08 | 6.50e-08 | 5.40e-08 | 4.04e-08 | 3.02e-08 | 1.06e-08 | 7.72e-10 | -3.32e-10 |
| closure, ours (heating/losses) | 0.71 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 |
| closure with the gamma_eff advection | 1.02 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 |
| stellar heating, K22 | 1.44e-07 | 1.20e-07 | 1.37e-07 | 1.70e-07 | 1.82e-07 | 1.61e-07 | 6.01e-08 | 1.01e-08 | 1.99e-09 |
| radiative cooling, K22 | 6.03e-08 | 6.79e-08 | 8.66e-08 | 1.30e-07 | 1.42e-07 | 1.25e-07 | 3.86e-08 | 2.99e-09 | 5.04e-10 |
| adiabatic, K22 | 4.94e-08 | 4.80e-08 | 4.18e-08 | 3.63e-08 | 3.24e-08 | 2.76e-08 | 1.63e-08 | 5.66e-09 | 1.35e-09 |
| advection, K22 | 1.29e-08 | 9.27e-09 | 9.37e-09 | 8.56e-09 | 8.20e-09 | 6.72e-09 | 4.61e-09 | 1.27e-09 | 1.24e-10 |

Ratios, ours over theirs:

| r/r_base | 1.05 | 1.1 | 1.15 | 1.2 | 1.25 | 1.3 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| photoheating | 0.67 | 1.27 | 1.40 | 0.95 | 0.70 | 0.64 | 0.81 | 1.25 | 1.07 |
| radiative cooling | 0.011 | 0.017 | 0.015 | 0.010 | 0.010 | 0.012 | 0.037 | 0.203 | 0.196 |
| adiabatic | 2.54 | 2.50 | 3.02 | 2.91 | 2.66 | 2.57 | 2.26 | 1.98 | 1.77 |
| advection | 0.69 | 3.42 | 6.94 | 6.31 | 4.92 | 4.49 | 2.30 | 0.61 | -2.69 |

#### 10.2 What moved against the section 4 table

**The reaction heat is a new channel and it is large.** At 1.05 `r_base` it is
5.26e-8 erg cm^-3 s^-1, against a photoheating total of 9.61e-8 there, and it
peaks at 71 percent of the total heating near 1.02 `r_base`. Section 4 had no
such row: before P53 the ionization energy that H2 photoabsorption spends was
returned to the gas by dissociative recombination and the code kept none of it.

**Our photoheating moved toward theirs.** The peak is now
1.94e-7 erg cm^-3 s^-1 at 1.141 `r_base` against 2.23e-7 at 1.111 before, and
Koskinen's 1.84e-7 at 1.238. The ratio of ours to theirs over 1.05 to 3.00
`r_base` is now 0.64 to 1.40, against 0.50 to 1.81 in section 4. Candidate (a)
of the judgment is still not supported, and by a slightly smaller margin.

**Our radiative cooling did not move.** It is 1.0 to 2.0 percent of theirs over
the front (1.29e-9 against 8.66e-8 at 1.15 `r_base`), the same one-to-two
percent section 5 refuted candidate (b) on. The H3+ infrared channel is
2.90e-10 at 1.05 `r_base` -- an order larger than section 4's 3.48e-11, because
the layer is hotter there -- and still 5e-14 by 1.10.

**The hydrodynamic terms moved up, not down.** Our adiabatic term is now 1.8 to
3.0 times theirs (was 1.6 to 2.9) and our advection 2.3 to 6.9 over 1.10 to
1.30 `r_base` (was 3.8 to 7.9). Both are the same statement section 5 made.

**The temperature comparison.** Ours is now 633 K at 1.05 `r_base` where
section 4 had 472 K, and 983 K at 1.15 where section 4 had 1224 K. The deficit
against Koskinen at 1.05 falls from 708 K to 547 K and at 1.15 rises from 493 K
to 734 K. Beyond 1.5 `r_base` it is unchanged to within 40 K.
`supersonic_molecular_base.md` section 14.8.4 carries the full column.


#### 10.3 The two channels separated

**The reaction heat carries essentially all of the effect on this planet.** The
section 10 state was re-finished once more with `Molecular reaction heat: False`
and everything else identical, same binary, and closed in 7.9 s with `info = 0`
on both gates (`||R|| = 5.087e-6`, gate flux spread 2.468e-4). The caloric EOS
has no input key, so this is **caloric EOS alone (B)** against **caloric EOS
plus reaction heat (C)**, not P53 off against P53 on; `gamma_eff` reads 1.409
in both columns, which is the check that the EOS is on in B.

| | B: caloric EOS only | C: EOS + reaction heat | C - B |
|---|---:|---:|---:|
| **`log10 Mdot` [g/s]** | 10.29 | 10.34 | **+0.05** |
| `heat_mol_chem/heat_total` at 1.02 `r_base` | 0.000 | 0.706 | +0.706 |
| layer `T` minimum [K] (at `r_base`) | 372.1 (1.0522) | 541.4 (1.0759) | +169.3 |
| `T` at the `f(H2) = 0.5` front [K] | 428.1 | 665.1 | +237.0 |
| H2 front `f = 0.5` [`r_base`] | 1.0394 | 1.0449 | +0.0055 |
| `f = 1e-2` / `f = 1e-4` | 1.0698 / 1.1447 | 1.0838 / 1.1851 | +0.014 / +0.040 |
| peak `n(H3+)` [cm^-3] (at `r_base`) | 1.539e5 (1.0037) | 1.456e5 (1.0021) | -8.3e3 |
| ionization front [`r_base`] | 2.3762 | 2.4804 | +0.104 |
| `T` at 1.05 / 1.15 / 1.50 / 3.00 `r_base` [K] | 374 / 1384 / 2793 / 2108 | 633 / 983 / 2683 / 2164 | +259 / -401 / -110 / +56 |
| deficit `T`(K22) - ours at 1.05 / 1.15 / 1.30 [K] | 806 / 333 / -125 | 547 / 734 / 152 | -259 / +401 / +277 |

**The two channels pull in opposite directions over the front.** The reaction
heat closes the 1.05 `r_base` deficit by 259 K and opens the 1.15, 1.20 and
1.30 deficits by 401, 385 and 277 K, so on the Koskinen comparison it is a
redistribution and not a closing -- the same conclusion section 141.4 of
`docs/Update_EXHALE_stage1.md` reached, now measured on states finished under the
face-flux measure. The caloric EOS alone leaves `log10 Mdot` at the pre-P53
10.29.

**Against the section 141.4 table (pre-G23).** The C column here reproduces
its C closely -- `log10 Mdot` 10.34 in both, layer minimum 541 K against 525 K,
`f = 0.5` front 1.0449 against 1.0429, peak `n(H3+)` 1.456e5 against
1.469e5 cm^-3 -- while B differs by more (minimum 372 K against 417 K), so the
reaction heat's increment on the rate agrees to 0.01 dex (+0.05 here, +0.046
there) and its increment on the layer minimum does not (+169 K here, +108 K
there).

#### 10.4 What section 10 does not settle

* **The channel split rests on one rung.** Section 10.3 separates the two on
  the gate rung; no other rung was run with `Molecular reaction heat: False`.
* **The closure result rests on one state.** The `gamma_eff` closure of the
  table above was evaluated on the gate rung only; the other five rungs of
  `supersonic_molecular_base.md` section 14.8 were not put through it.
* **The `gamma_eff` used in the closure is reconstructed in post-process**, from
  `p = (n_tot + n_e) k T` and the run's own `n(H2)`, following
  `src/modules/states/caloric_eos.f90`. The code does not write `gamma_eff` to
  any output file, so this is the same formula evaluated independently and not
  the value the run integrated. The two should agree to the accuracy of the
  code's 4096-node Hermite table for `u_rv`; that was not verified against an
  instrumented run.
* **The residual norm is under review**, per the box at the head of this
  section.
* **Everything section 10 of the old numbering listed as unmeasured is still
  unmeasured**, except the H2 dissociation ledger item, which P53 closed and
  which is marked accordingly in section 11.

---

## 11. What was not measured

* **The grid test of section 8 is one refinement, not a convergence study.**
  Both refined runs converged (`info = 0` on both gates) and both answered the
  question they were set, but a single factor-2 refinement of the base cells
  shows that a quantity moves or does not move; it does not establish that
  either state is grid converged. The He/H = 0.3 rung in particular moved its
  H3+ peak by a factor 5.4 under that one refinement, so nothing on that rung
  should be quoted as converged in `dr` until a third grid is run.
* **Nothing here re-runs the other four rungs of the 1 microbar ladder on the
  refined base grid.** If the He/H = 0.3 first cell was distorting that rung's
  H3+ peak by a factor 5.4, the same question is open for the rest of
  `supersonic_molecular_base.md` section 14.2, and it is not answered here.
* **Nothing here tests the H2 abundance itself.** The claim of section 7 is
  conditional: *given* that our H2 is gone above 1.05 `r_base` and theirs is
  not, the temperature and mass-flux differences follow. Whether ours should be
  gone is item (P43) and the transported-closure question of
  `supersonic_molecular_base.md` section 13.7, neither of which is touched.
* **The spectrum is not matched.** Ours is a power law of index -1 over
  13.6 eV to 1.24 keV; theirs is the mean solar spectrum of Koskinen et al.
  (2013a, 2013b) over their own range. The integrated flux agrees to
  2.5 percent, the shape was not compared, and the inward shift of our heating
  peak by 0.13 `r_base` was not attributed.
* **Conduction was not switched on** in any run here, so how much of the
  outer-region gap a Watson et al. (1981) conductivity would close is not
  known. Their own diagnostic says it would be small.
* **The chemical energy of H2 dissociation** was an open ledger item here and
  is **CLOSED on our side by P53** (marked 2026-09-03): the collisional H2/He
  network now deposits its chemical energy, and on the section 10 state that
  channel carries 71 percent of the total heating at its peak, and section 10.3
  measures what the run does without it. Whether their code charges it was
  still not checked.
* **The digitized values below 1.05 `r_base`** are order of magnitude only, for
  the reason section 2 gives, so the first column of every table above carries
  that caveat on the Koskinen side and the section 9.2 caveat on ours.

---

## 12. Files and reproduction

Campaign directory (scratch): `.../scratchpad/p23heat/`.

* `dig/` -- the 600 dpi render, the tracking code (`track.py`, `run_track.py`,
  `run_track2.py`), the extracted curves (`k22fig9_curves.json`,
  `k22fig9_conduction.json`) and the overlay used to accept the tracking.
* `ab_x1/`, `ab_x0p5/`, `ab_x0/`, `ab_x3/` -- the H3+ A/B of section 6.
  `ab_x1` is also the state every "ours" number above is read from.
* `diag_heh30c_r1/`, `diag_heh0p3c_r1/`, `diag_heh0p0793c_r1/` -- the three
  rungs re-converged with the photoionization-rate dump, for section 8.
* `grid/heh30_fine/`, `grid/heh0p3_fine/` -- the base-grid refinement runs of
  section 8.
* `tree/` -- the scratch source tree and its binary; the three changes in it are
  the section 9.1 fix, the `Photorates_diag.txt` dump and `EXHALE_H3P_SCALE`.
* `p52_heating_breakdown_fix.diff` -- the section 9.1 fix ALONE, the one thing
  here meant for the repository; `util_ion_eq_full.diff` is the same plus the
  two scratch-only changes. `tree_p52/` is a clean copy of `src/` with only the
  fix applied and `verify_p52/` the run that verified it (section 9.1).
* `ours.py`, `full.py`, `tabs.py`, `mu.py`, `h2budget.py`, `make_fig.py` --
  the analysis and the figure.

The figure is `docs/lower_atmosphere_figs/fig_p23_thermal_budget.png`; it is
drawn by `make_fig.py` in the campaign directory and is not regenerated by
`docs/lower_atmosphere_figs/make_figures.py`.

---

## 13. After block I

*(Numbered 13, not 11: sections 11 and 12 were already taken, and renumbering
them would break the references other files carry.)*

**Judgment first: the P23 gate is NOT passed, and block I does not move it.
The temperature deficit against Koskinen et al. (2022) Model A is unchanged to
within 3 K at every radius of section 2's grid -- +550, +973, +731, +450, +150,
+187, +1434 and +2578 K at 1.05 through 3.00 `r_base`, against +547, +974,
+734, +453, +152, +188, +1433 and +2578 K before -- and our H2 is still gone
where theirs survives: `f(H2)` falls to 3.1e-3 at 1.10 `r_base` where Model A
still carries 0.976, and to 1.3e-6 at 3.00 where Model A carries 0.403. H3+ is
still inverted in shape: ours peaks at 1.4e5 cm^-3 in the first few cells and
falls below 1 cm^-3 by 1.10 `r_base`, while Model A rises to 1.2e4 cm^-3 at
1.30 and holds 1.9e3 at 3.00. What remains is what section 7 already
attributed: the layer is thermally starved because the H2 is not there to
absorb, and the outer gap is both composition and energy.**

**Why the two block-I channels do not decide it, and the two reasons are
different in kind.** The brief for this campaign expected section 150's halving
of the Lyman-Werner dissociation rate to push toward MORE surviving H2 and
section 151's measured cross sections (85 times the old fit at threshold) to
push toward LESS. Measured, section 150 cannot act on this ladder at all, and
section 151 does act, is properly resolved, and is simply small:

* **The Lyman-Werner channel is switched off on this run, so section 150 is
  multiplying zero.** The setup report says so in as many words:
  `WARNING H2 Lyman-Werner photodissociation OFF: no "Stellar LW flux" key and
  no spectrum to integrate one from`, and it goes on: `The band exists whenever
  the star does, so this is a MISSING CHANNEL, not a modelling choice: H2 keeps
  a sink it should not.` The ladder's spectrum is the analytic power law of
  index -1 from 13.6 eV, there is no `Stellar LW flux` key in `input.inp`, and
  `heat_H2_LW` is 0.00e+00 at every radius of section 10.1 both before and
  after. **This ladder therefore cannot test section 150 at all.** Testing it
  needs a rung that supplies a `Stellar LW flux` or reads a spectrum file
  covering the band, and that is a separate campaign, not a re-reading of these
  runs.
* **Section 151 IS resolved on this grid, and its effect is genuinely small.**
  The power-law branch of `src/modules/init/set_energy_vectors.f90` does not
  bin the spectrum into the three numbers `input.inp` names; those are the
  endpoints of a logarithmic grid. It lays `num_HI = 50` points between
  `e_th_HI` = 13.6 eV and `e_th_HeI` = 24.6 eV, and **13 of them fall inside
  the 15.4-18.0 eV window Backx et al. (1976) measured, at a spacing of 0.185
  to 0.211 eV** (15.494, 15.679, 15.866, ... 17.862 eV). The threshold
  structure is carried on the grid. With it carried, the `heat_H2` channel at
  1.05 `r_base` reads 7.29e-9 against 7.84e-9 before, **-7 percent** (the two
  new mutually exclusive channels, H2 double ionization above 51.4 eV and H2
  neutral dissociation over 33-41 eV, both appear in the setup report as
  RESOLVED). `heat_H2` is 7.6 percent of the photoheating there, so that is
  **-0.6 percent of the local photoheating**. The 85-times factor at threshold
  does not become a large change in the rate because the channel it corrects is
  a small part of the budget on this planet, not because the grid cannot see
  it.

So the verdict is not "one direction won". **Section 151 was measured on this
ladder and is small; section 150 was not measured at all, because the
configuration switches its channel off.** Either way the H2 survival problem is
untouched by block I.

**The state.** The He/H = 0.0793 rung of `supersonic_molecular_base.md`
section 14.9 (`runs/heh0p0793_r1` in `.../scratchpad/ladder4/`), the same
inputs as section 10's state, re-finished under the current tree binary
(md5 `063b2673`, built 2026-09-04 02:23 KST) with `Load IC? True`,
`Solver: Newton 100.0`, `du_th [PLM,WENO3]: 1.0e9 1.0e-3` and
`Secondary_ionization: Immediate`. Three JFNK attempts, then `info = 0`;
`log10 Mdot = 10.34` g/s, unchanged. Its two gate numbers are the
own-composition residual `||R|| = 1.485e-4` (section 155: the residual gate is
OFF in this input, `resid_th = -1`, so this number is reported and not tested)
and the face mass-flux spread 1.135e-13 over `r >= 1.2` against a tolerance of
2e-5. Section 10's `rho v r^2` caveat is partly retired: the face flux is flat
to 1.5e-10 over `r >= 1.03`, the cell-1 velocity is +14.4 cm/s instead of
-15.0, and the cell-centre flux hole at 1.0006 `r_base` (0.73 of the wind
value) is gone. The cell-centre product over `r >= 1.03` still reads 3.42e-3,
the same as before, so the non-flatness ABOVE 1.03 `r_base` in that functional
is not retired.

### 13.1 The comparison, re-measured

Section 2's radii; the Koskinen column is the same digitization, and the
sub-1.05 `r_base` entries carry section 2's order-of-magnitude caveat on their
side.

| r/r_base | 1 | 1.05 | 1.1 | 1.15 | 1.2 | 1.3 | 1.5 | 2 | 3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **T [K], block I** | 1117 | 630 | 593 | 986 | 1451 | 2117 | 2684 | 2686 | 2164 |
| T [K], section 10 (pre-block-I) | 1117 | 633 | 592 | 983 | 1448 | 2115 | 2683 | 2687 | 2164 |
| T [K], Koskinen Model A | 1082 | 1180 | 1566 | 1717 | 1901 | 2267 | 2871 | 4120 | 4742 |
| **deficit T(K22) - ours, block I** | -35 | +550 | +973 | +731 | +450 | +150 | +187 | +1434 | +2578 |
| deficit, section 10 | -35 | +547 | +974 | +734 | +453 | +152 | +188 | +1433 | +2578 |
| ours/Koskinen, block I | 1.03 | 0.53 | 0.38 | 0.57 | 0.76 | 0.93 | 0.93 | 0.65 | 0.46 |
| **f(H2), block I** | 0.958 | 0.393 | 0.0031 | 0.000262 | 7.84e-05 | 3.95e-05 | 1.98e-05 | 5.08e-06 | 1.33e-06 |
| f(H2), section 10 | 0.958 | 0.393 | 0.00312 | 0.000266 | 7.88e-05 | 3.96e-05 | 1.98e-05 | 5.19e-06 | 1.47e-06 |
| **f(H2), Koskinen Model A** | 0.984 | 0.982 | 0.976 | 0.963 | 0.936 | 0.843 | 0.688 | 0.508 | 0.403 |
| **n(H3+) [cm^-3], block I** | 1.40e+05 | 7.87e+02 | 1.21e+00 | 3.14e-02 | 8.08e-03 | 2.41e-03 | 3.26e-04 | 9.95e-06 | 5.20e-07 |
| n(H3+), section 10 | 1.42e+05 | 7.87e+02 | 1.22e+00 | 3.18e-02 | 8.13e-03 | 2.42e-03 | 3.28e-04 | 1.02e-05 | 5.73e-07 |
| **n(H3+), Koskinen Model A** | 1.90e+04 | 3.14e+03 | 4.99e+03 | 7.02e+03 | 9.10e+03 | 1.16e+04 | 1.06e+04 | 7.16e+03 | 1.88e+03 |

### 13.2 The channel table, re-measured

Section 10.1's table on the block-I state. `heat_total` and the solver's own
`heat` column now agree to **2.1e-12** over `r > 1.001` (section 10 reported
1.6e-6 on its state, and 1.4e-3 before section 9.1's fix), so the dump is the
run's own heating to round-off.

| r/r_base | 1.05 | 1.1 | 1.15 | 1.2 | 1.25 | 1.3 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| **T [K], ours** | 630 | 593 | 986 | 1451 | 1829 | 2117 | 2684 | 2686 | 2164 |
| T [K], Koskinen | 1180 | 1566 | 1717 | 1901 | 2084 | 2267 | 2871 | 4120 | 4742 |
| f(H2), ours | 0.393 | 0.0031 | 0.000262 | 7.84e-05 | 5e-05 | 3.95e-05 | 1.98e-05 | 5.08e-06 | 1.33e-06 |
| gamma_eff, ours | 1.5783 | 1.6660 | 1.6666 | 1.6666 | 1.6667 | 1.6667 | 1.6667 | 1.6667 | 1.6667 |
| photoheating, ours | 9.61e-08 | 1.53e-07 | 1.92e-07 | 1.61e-07 | 1.28e-07 | 1.03e-07 | 4.88e-08 | 1.25e-08 | 2.14e-09 |
| &nbsp;&nbsp;H I | 9.22e-09 | 5.48e-08 | 8.93e-08 | 8.98e-08 | 7.94e-08 | 6.83e-08 | 3.70e-08 | 1.01e-08 | 1.69e-09 |
| &nbsp;&nbsp;He I | 2.68e-08 | 7.66e-08 | 8.83e-08 | 6.04e-08 | 4.04e-08 | 2.82e-08 | 9.64e-09 | 2.09e-09 | 3.72e-10 |
| &nbsp;&nbsp;H2 | 7.29e-09 | 2.19e-10 | 3.41e-11 | 1.05e-11 | 5.98e-12 | 4.08e-12 | 1.14e-12 | 9.93e-14 | 7.10e-15 |
| &nbsp;&nbsp;He(2^3S)+H Penning | 4.74e-12 | 7.64e-09 | 1.06e-08 | 8.09e-09 | 6.15e-09 | 4.67e-09 | 1.56e-09 | 2.21e-10 | 3.48e-11 |
| &nbsp;&nbsp;He recombination photons | 3.36e-12 | 3.69e-09 | 1.99e-09 | 7.54e-10 | 3.68e-10 | 2.09e-10 | 4.47e-11 | 4.19e-12 | 2.18e-13 |
| &nbsp;&nbsp;H2 Lyman-Werner | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;FUV photolysis | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;molecular reaction heat | 5.27e-08 | 1.03e-08 | 2.05e-09 | 1.43e-09 | 1.11e-09 | 8.68e-10 | 3.13e-10 | 4.62e-11 | 7.32e-12 |
| radiative cooling, ours | 6.30e-10 | 1.18e-09 | 1.29e-09 | 1.31e-09 | 1.38e-09 | 1.45e-09 | 1.44e-09 | 6.09e-10 | 9.86e-11 |
| &nbsp;&nbsp;H recombination | 2.73e-10 | 9.31e-10 | 9.99e-10 | 9.84e-10 | 1.01e-09 | 1.02e-09 | 9.07e-10 | 3.64e-10 | 5.98e-11 |
| &nbsp;&nbsp;bremsstrahlung | 7.43e-11 | 2.45e-10 | 2.94e-10 | 3.20e-10 | 3.54e-10 | 3.79e-10 | 3.72e-10 | 1.51e-10 | 2.28e-11 |
| &nbsp;&nbsp;H3+ infrared | 2.83e-10 | 4.92e-14 | 4.42e-16 | 2.21e-16 | 2.53e-16 | 2.43e-16 | 6.10e-17 | 1.87e-18 | 5.62e-20 |
| &nbsp;&nbsp;H2 infrared | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| &nbsp;&nbsp;He I collisional excitation | 2.26e-21 | 1.40e-17 | 4.86e-14 | 2.26e-12 | 1.48e-11 | 4.24e-11 | 1.62e-10 | 9.31e-11 | 1.60e-11 |
| adiabatic, ours (gamma = 5/3) | 1.26e-07 | 1.20e-07 | 1.26e-07 | 1.06e-07 | 8.61e-08 | 7.10e-08 | 3.68e-08 | 1.12e-08 | 2.39e-09 |
| advection, ours (gamma = 5/3) | 9.10e-09 | 3.22e-08 | 6.51e-08 | 5.39e-08 | 4.03e-08 | 3.01e-08 | 1.06e-08 | 7.68e-10 | -3.32e-10 |
| advection, ours (gamma_eff) | -3.24e-08 | 3.20e-08 | 6.51e-08 | 5.39e-08 | 4.03e-08 | 3.01e-08 | 1.06e-08 | 7.68e-10 | -3.32e-10 |
| closure, ours (heating/losses) | 0.71 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 |
| closure with the gamma_eff advection | 1.02 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 |

Koskinen et al. (2022) Model A, the same digitization section 10.1 uses:

| r/r_base | 1.05 | 1.1 | 1.15 | 1.2 | 1.25 | 1.3 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| stellar heating, K22 | 1.44e-07 | 1.20e-07 | 1.37e-07 | 1.70e-07 | 1.82e-07 | 1.61e-07 | 6.01e-08 | 1.01e-08 | 1.99e-09 |
| radiative cooling, K22 | 6.03e-08 | 6.79e-08 | 8.66e-08 | 1.30e-07 | 1.42e-07 | 1.25e-07 | 3.86e-08 | 2.99e-09 | 5.04e-10 |
| adiabatic, K22 | 4.94e-08 | 4.80e-08 | 4.18e-08 | 3.63e-08 | 3.24e-08 | 2.76e-08 | 1.63e-08 | 5.66e-09 | 1.35e-09 |
| advection, K22 | 1.29e-08 | 9.27e-09 | 9.37e-09 | 8.56e-09 | 8.20e-09 | 6.72e-09 | 4.61e-09 | 1.27e-09 | 1.24e-10 |


Ratios, ours over theirs:

| r/r_base | 1.05 | 1.1 | 1.15 | 1.2 | 1.25 | 1.3 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| photoheating | 0.67 | 1.28 | 1.40 | 0.95 | 0.70 | 0.64 | 0.81 | 1.24 | 1.07 |
| radiative cooling | 0.010 | 0.017 | 0.015 | 0.010 | 0.010 | 0.012 | 0.037 | 0.203 | 0.196 |
| adiabatic | 2.54 | 2.51 | 3.02 | 2.91 | 2.66 | 2.57 | 2.26 | 1.98 | 1.77 |
| advection | 0.71 | 3.47 | 6.95 | 6.30 | 4.91 | 4.48 | 2.29 | 0.60 | -2.69 |

**What moved against section 10.1, and it is only two things.** The H2
photoheating at 1.05 `r_base` falls 7.84e-9 to 7.29e-9 (-7 percent) and the H3+
infrared cooling there falls 2.90e-10 to 2.83e-10 (-2 percent), both inside the
molecular layer. Every one of the four ratios against Koskinen is unchanged to
the two digits it is quoted to except the advection at 1.05 `r_base` (0.69 to
0.71) and at 1.10 (3.42 to 3.47). The closure rows are identical: 1.02 at 1.05
`r_base` with the `gamma_eff` advection and 1.00 from 1.10 outward, so section
10's closure result -- that section 9.2's failure was the specific heat plus an
unsteady layer -- survives block I unchanged. **The section 5 attribution is
therefore also unchanged**: our radiative cooling is 1.0 to 2.0 percent of
theirs over the front, our photoheating agrees to a factor 1.6 either way, and
our adiabatic and advection terms are 1.8 to 3.0 and 2.3 to 7.0 times theirs.

### 13.3 What section 13 does not settle

* **Everything section 10.4 lists is still unsettled**, and for the same
  reasons; nothing here re-runs the other rungs, adds a third grid, switches on
  conduction, or matches the spectrum.
* **Section 150 was not tested, and this ladder cannot test it.** Its channel
  is off in this configuration, so nothing here bears on how large the
  Lyman-Werner geometry factor is where the channel is on. A rung carrying a
  `Stellar LW flux` or a spectrum file is what would say, and it is a separate
  campaign.
* **Section 151 WAS tested here, on one planet and one state.** The 15.4-18 eV
  structure is resolved by 13 grid points at 0.185 to 0.211 eV, and the
  measured effect is -7 percent on `heat_H2` and -0.6 percent on the local
  photoheating at 1.05 `r_base`. That the effect is small on THIS planet, whose
  H2 is nearly gone by 1.05 `r_base` already, says nothing about a planet whose
  molecular layer is thicker.
* **`gamma_eff` is still reconstructed in post-process**, as in section 10.4.
* **Two rungs of the ladder this state belongs to needed repeated restarts**
  (He/H = 3 four rounds, He/H = 10 nine; section 14.9.3 of
  `supersonic_molecular_base.md`). All seven are accepted, and the gate rung
  read here closed on the first round, so nothing in this section rests on
  that -- but the ladder is not as cheap to re-finish as section 14.8's was.
* **`ntot_bc_per_H` rose 0.55 percent on this rung** (0.5434782 to 0.5464827)
  on an unchanged `q_H2_base`, so the base number density is 0.55 percent
  lower than section 10's state carried. Which block-I section does that was
  not traced. It is 0.002 dex on the rate and below what the run prints.

Campaign directory (scratch): `.../scratchpad/ladder4/`, with `table4.py`,
`cmp4.py`, `p23tab.py` and `make_fig4.py`. The channel table above is
`p23tab_blockI.md` there.

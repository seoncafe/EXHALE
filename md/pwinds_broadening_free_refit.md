# LHS 1140 b He 10830: what the p-winds retrieval returns when the broadening is free

Cherubim et al. (2026) retrieved the mass-loss rate, outflow temperature and
H:He ratio of LHS 1140 b from the WINERED He I 10830 line with p-winds, with
the line broadening set by the model rather than fitted: thermal broadening at
the retrieved temperature, plus the Lampon et al. (2020) turbulence term, which
is a fixed function of that same temperature. The width measurements find that the measured line is about
22 km/s FWHM wider than any 1-D steady-state model produces, and that the
measurement side accounts for at most 0.4 km/s of it.

This memo asks the obvious follow-up: **if the broadening is made a free
parameter, do the retrieved composition and mass-loss rate survive?**

Run directory: `../LHS1140b/pwinds_refit/` (scripts, grids, logs, figures).
Everything below was computed there; nothing is quoted from the earlier oracle
except where marked.

---

## Verdict

**The broadening goes to v_nt = 7.5-8 km/s (1-D sigma). At that point the
line width, depth and equivalent width are all reproduced; the temperature
leaves the published interval; the mass-loss rate rises by about an order of
magnitude; and the composition stops being determined, retaining only a lower
bound on He/H.** Fitting the
released spectrum with the authors' own model, data and likelihood, and
freeing only the non-thermal dispersion:

1. **The kernel is strongly preferred.** With the equivalent width held at the
   measurement, chi^2 falls from 218.5 to 19.7 on 129 pixels for one added
   parameter, and the fitted line goes from 15.5 to 23.0 km/s FWHM against the
   23.3 km/s measured. The unconstrained joint control finds the same thing
   independently (chi^2 214.7 -> 18.7, v_nt = 8 km/s). Even with the
   composition and mass-loss rate frozen at the published values, 5 km/s of
   dispersion improves the fit by delta chi^2 = 148.
2. **The composition stops being determined.** Profiled over temperature and
   `v_nt` with the equivalent width held, chi^2 is 19.7-22.0 at *every*
   He/H between 30 and 3e4 -- flat to within delta chi^2 = 2.3 across three
   orders of magnitude -- and only rises below He/H ~ 9. The unconstrained
   joint control agrees, 18.7-22.6 from He/H = 17 to 3e4. At the fixed C26
   broadening the same axis is not flat: it has a minimum near
   He/H ~ 3e2-3e3 and rises by 60-170 in chi^2 towards the hydrogen-rich end,
   which is the shape that supports a helium-dominated reading. **Stated as a
   fact about this experiment: the He 10830 line constrains the composition in
   this model only as long as the broadening is prescribed.** The published
   He/H ~ 1e3 remains an allowed solution; it stops being a selected one, and
   what survives is a lower bound, He/H greater than about 10-17.
3. **The mass-loss rate moves up by about an order of magnitude.** The best
   equivalent-width-held solution carries Mdot = 2.1e9 g/s against the
   published 2.03e8 (+0.67/-0.58)e8 g/s, i.e. +28 sigma in the units of their
   published error bar; the flat composition ridge spans 3e8-6e9 g/s, and the
   unconstrained control gives 5.6e8 with an error-rescaled interval of
   4.2e8-6.8e9. The
   upper part of that range meets the 2e9 g/s the paper identifies as the rate
   needed to drag atomic oxygen, and approaches the 5e9 g/s it adopts as an
   evolutionary ceiling.
4. **The temperature falls below their prior.** The stage-2 minimum is at
   4800 K and the joint control at 4500 K, against the published
   5160 +46/-50 K (-7.2 sigma and lower). This is the expected direction: at
   fixed broadening the temperature is the only knob that can widen the line,
   and removing that duty lets it fall.
5. **The two stages do not separate cleanly.** Because the line is saturated,
   pinning the equivalent width still leaves Mdot coupled to the broadening.
   At the published composition and temperature the Mdot that holds the
   equivalent width varies between 0.73x and 1.99x of its v_nt = 0 value over
   0-18 km/s, and the equivalent width also selects two mass-loss branches
   rather than one. The residual correlation is a factor of about 2, not
   negligible but far smaller than the order-of-magnitude motion of the
   best-fit Mdot itself.

**Consistency with the width memos: confirmed, independently.** They put the
missing broadening at 21.96 km/s FWHM (sigma = 9.33 km/s) on top of an EXHALE
line of sigma = 2.42 km/s. Here a p-winds line that already carries thermal
plus Lampon turbulence (sigma = 4.4 km/s) asks for 7.5-8 km/s more, giving a
total dispersion of 8.6-9.1 km/s against the 9.6 km/s the memo's requirement
implies in quadrature, and a fitted line width of 23.0 km/s FWHM against the
memo's 21.96 km/s requirement and the 23.3 km/s measured. The memo's value
lies inside the error-rescaled 1 sigma interval of the freely fitted
dispersion (6.8-9.3 km/s). A second, independently written retrieval, given
the freedom, asks for essentially the amount of extra broadening the memo says
is missing.

*Interpretation, stated as such:* what these numbers show is that the
published composition, temperature and mass-loss rate are conditional on the
broadening prescription, not that they are wrong. The fixed Lampon term is a
modelling choice the field makes generally, and the paper is explicit about
using it. What changes when it is relaxed is that width and column trade
against each other in a saturated line, and the composition -- which entered
the published result through the depth and shape of a line whose width was
pinned to the temperature -- no longer has an independent handle. Nothing here
identifies the physical origin of the 7.5-8 km/s.

---

## Method

### The forward model is the authors' own

The model, the data, the likelihood, the SED and its normalization, the system
parameters and the priors are taken from the script released with the paper
(`LHS1140b_zenodo/lhs1140b_2024B_pwinds.py`, Zenodo 15723779):

- p-winds 1.4.7, isothermal Parker wind, `relax_solution`, `exact_phi`;
- SED: const-res v23 GJ 1132 scaled by the catalog Lx ratio 0.59 and the
  distance factor;
- R_p = 0.154 R_J, M_p = 0.0176 M_J, R_p/R_* = 0.072092, b = 0.23;
- 5-phase transit averaging over [-0.5, +0.5], grid 100 px supersampled 5x;
- `wind_broadening_method='average'`, `turbulence_broadening=True`;
- bulk line-of-sight velocity v_wind = 2.26 km/s, held fixed as in their grid;
- data: the released in-transit spectrum
  `lhs1140b_trans_spec_GP_IT_vac.pickle` converted to air, with the GP
  component added back, 129 pixels;
- instrument profile: 4.4 km/s FWHM Gaussian, convolved on the data grid;
- likelihood: theirs, Gaussian with the `log(2 pi sigma^2)` term kept.

This is a stricter reproduction than `../LHS1140b/pwinds_oracle/`, which
predates the release of the authors' script and was missing the turbulence
switch and the phase averaging.

**Validation.** At nine nodes of the likelihood cube the authors released
(`lhs1140b_grid_likelihood_results_1132_v26e3.pickle`) the reproduction agrees
with their stored values to |delta ln L| < 0.03 out of ln L ~ 630, the
difference being the stochastic tail of the relaxation iteration
(`pwinds_refit/verify.log`). At their published parameter vector the line
metrics (red depth 1.376 %, FWHM 0.575 A) reproduce the independent oracle
run `pwinds_oracle/reproduce_authors_fig4.py` to the third decimal.

### The free parameter

`v_nt` is an isotropic non-thermal velocity dispersion (1-D Gaussian sigma)
added **in quadrature to the Doppler width of the opacity**, alongside the
thermal and turbulence terms:

```
sigma_v^2 = k T / m_He  +  v_wind_broad^2  +  (5/6) k T / m_He  +  v_nt^2
             thermal        Parker LOS        Lampon turbulence   free
```

In p-winds' `'average'` mode the temperature argument of
`radiative_transfer_2d` enters nowhere except this width, so passing
`turbulence_broadening=False` with an effective temperature
`T_rt = (11/6) T + m_He v_nt^2 / k_B` implements it exactly; the atmospheric
model always sees the true T. The `v_nt = 0` case is bit-identical to calling
p-winds with the authors' `turbulence_broadening=True`
(`pwinds_refit/verify.log`, max |diff| = 0).

**Applying it to the opacity rather than to the emergent spectrum matters,
because the line is saturated.** At the published parameter vector, 9 km/s of
in-opacity broadening raises the red-component equivalent width from 0.00839 to
0.01268 A (+51 %) and lowers the peak depth from 1.376 % to 1.084 %, while
convolving the unbroadened model with the same 9 km/s Gaussian conserves the
equivalent width at 0.00836 A. A post-hoc kernel is therefore not the same
experiment, and it fits worse (chi^2 207 against 166).

### Two-stage procedure

The main line of the analysis follows the logic of the width memos, in which
the equivalent width fixes the absorbing column and the width is a separate,
kinematic question:

- **Stage 1 -- the equivalent width.** The red component of the line (the
  blended 10833.217/10833.306 A pair) has a measured equivalent width. With the
  broadening held at the authors' value, (T, h_fraction) are placed on a grid
  and log Mdot is solved for so that the model reproduces that equivalent
  width. The equivalent width is one number, so this defines a surface, not a
  point; the surface is what stage 1 maps.
- **Stage 2 -- the shape.** On that surface `v_nt` is freed. Because the line
  is saturated the equivalent width itself grows with `v_nt`, so log Mdot is
  re-solved at every node: the equivalent width stays pinned to the
  measurement and only the line shape is being fitted. The amount by which
  Mdot has to move to keep the equivalent width fixed as `v_nt` grows is the
  residual coupling between the two stages, and it is reported.

Aperture: the red trough, 10832.566-10834.166 A in vacuum (the aperture the
oracle uses), 44 of the 129 pixels. The quoted equivalent-width uncertainty
propagates the errors of individual pixels only; the released spectrum carries a Gaussian
process whose correlations are not propagated, so it is a lower bound.

**The equivalent width does not determine Mdot uniquely.** It rises with Mdot,
peaks, and falls again as the denser wind loses its metastable fraction, so
the measured value is generally reached twice, and the two solutions fit the
line profile very differently. At T = 5160 K, He/H = 999, v_nt = 6 km/s the
crossings are at log Mdot = 8.2 and 8.8, with chi^2 = 73 and 30. Both roots
are therefore found at every node and the one with the better full-spectrum
chi^2 is kept. (A first pass that took only the lower root gave a systematically
worse and composition-dependent stage 2; the numbers below are from the
two-root version.)

Two grids are run in this construction:

| grid | T [K] | log h_fraction | note |
|---|---|---|---|
| `authors_prior` | 4500-7000 | -4.5 to -1.5, step 0.25 | inside the support of the authors' own grid search (He/H >= 30.6) |
| `extended` | 4000-6500 | -4.5 to -0.25, step 0.5 | composition and temperature carried past their prior |

with `v_nt` on 0-18 km/s in both. The authors' composition prior already
excludes He/H < 30.6 and their temperature prior stops at 4500 K; where a
result sits on one of those edges it is said so explicitly.

### Control: joint likelihood grid

As a control, the same likelihood is evaluated on a joint grid without the
equivalent-width constraint, in the style of the authors' own grid search:
log Mdot on 25 nodes over [7, 10], T on 11 nodes over [4500, 7500] K,
log h_fraction on 17 nodes over [-4.5, -0.5], once at v_nt = 0 and once with
v_nt on 12 nodes over [0, 20] km/s. Profile-likelihood intervals
(delta chi^2 = 1) along each axis follow from it. This replaces an earlier
MCMC attempt, which is not part of the record: profiling a grid is both the
authors' own procedure and much the cheaper way to get the comparison the
question needs.

---

## Results

Figures live with the run, not in `docs/figures/`:
`../LHS1140b/pwinds_refit/fig1_spectrum_fits.pdf` (the spectrum with the
published vector, the stage-1 and the stage-2 model, plus residuals) and
`fig2_broadening.pdf` (a: chi^2 against `v_nt`; b: the composition axis with
the broadening fixed and free; c: the residual coupling).

### The reference points

All chi^2 values are on 129 pixels with the released errors.

| model | chi^2 | red EW [A] | full EW [A] | red depth | FWHM |
|---|---|---|---|---|---|
| flat continuum | 2009.6 | 0 | 0 | 0 | - |
| C26 published vector | 212.0 | 0.008394 | 0.011273 | 1.376 % | 0.575 A = 15.9 km/s |
| observed | - | 0.011223 | 0.013023 | 1.254 % | 0.841 A = 23.3 km/s |

The published vector is (Mdot 2.03e8 g/s, T 5160 K, H:He 1.01e-3, i.e.
He/H = 989). Its chi^2 per degree of freedom is 1.68. Two things follow
directly, before any refitting:

- **the model line carries 25 % less red-component equivalent width than the
  measurement** (13 % less over the whole triplet window; our full-window
  observed value, 0.01302 A, reproduces the 0.01304 A obtained by direct
  integration);
- and it is 1.5x too narrow, 15.9 against 23.3 km/s FWHM, which is the
  discrepancy the width memos describe.

### Stage 1: the equivalent width alone

With the broadening at the authors' value, **no point in the (T, h_fraction)
region their prior allows reproduces the observed red equivalent width at any
mass-loss rate in [1e7, 1e10] g/s, except at the cold edge of their
temperature prior (T = 4500-4800 K).** At T = 5160 K the largest reachable value is 0.00871 A
against the measured 0.01122 A.

The family at T = 5160 K, each node at the Mdot that comes closest to the
measured equivalent width:

| He/H | log Mdot | EW [A] | EW - EW_obs | chi^2 (full spectrum) |
|---|---|---|---|---|
| 31622 | 8.450 | 0.006448 | -43 % | 315.6 |
| 9999 | 8.532 | 0.007326 | -35 % | 244.2 |
| 3161 | 8.586 | 0.008088 | -28 % | 218.5 |
| 999 (C26) | 8.612 | 0.008601 | -23 % | 227.2 |
| 315 | 8.632 | 0.008678 | -23 % | 250.4 |
| 99 | 8.633 | 0.008062 | -28 % | 273.6 |
| 30.6 | 8.588 | 0.006465 | -42 % | 357.4 |

Where the equivalent width *is* reachable at fixed broadening -- twelve nodes,
all at the two coldest temperatures of the grid, T = 4500 and 4800 K, against
a cold prior edge of 4500 K -- the spectral fit is poor: chi^2 = 471-518 at ten of them
and 851 and 974 at the other two, because the line that carries the right area
there is far too deep and too narrow.

**Stage 1 finding.** At the C26 broadening the equivalent width and the line
profile cannot be matched at the same time. The best full-spectrum chi^2 on
the equivalent-width-matched family is 218.5, no better than the published
vector's 212.0.

### Stage 2: the broadening freed, the equivalent width held

Freeing `v_nt` while re-solving Mdot at every node so that the equivalent
width stays pinned changes the picture completely. Both grids return the same
minimum, and it lies inside the authors' composition prior:

| | chi^2 | v_nt [km/s] | T [K] | He/H | Mdot [g/s] | red EW | FWHM |
|---|---|---|---|---|---|---|---|
| stage 1 | 218.5 | 0 | 5160 | 3161 | 3.9e8 | -28 % | 15.5 km/s |
| stage 2 | 19.7 | 7.5 | 4800 | 3161 | 2.1e9 | matched | 23.0 km/s |
| observed | - | - | - | - | - | 0.011223 A | 23.3 km/s |

`delta chi^2 = 198.8` between stage 1 and stage 2, for one added parameter.
At the stage-2 minimum the three independent observables all land together:
red equivalent width 0.011218 A against 0.011223 measured, full-window
0.013043 against 0.013023, peak depth 1.277 % against 1.254 %, FWHM 23.02
against 23.29 km/s. The residuals (Figure 1) are flat where the
fixed-broadening models miss the blue component and both wings by 2-3 sigma.

**The broadening lands at v_nt = 7.5 km/s at essentially every node of both
grids**, and the chi^2 curve against `v_nt` (Figure 2a) is a sharp minimum:
218 at 0 km/s, 20 at 7.5, back above 200 by 15-18 km/s.

### Where the C26 parameters go

**Composition -- the central result.** Profiled over T and `v_nt` with the
equivalent width held (`authors_prior` grid, the finer one):

| He/H | chi^2, v_nt = 0 | chi^2, v_nt free | v_nt at the minimum | T | log Mdot |
|---|---|---|---|---|---|
| 31622 | 234.6 | 20.2 | 7.5 | 4500 | 9.24 |
| 17782 | 239.1 | 20.1 | 7.5 | 4500 | 9.37 |
| 9999 | 244.2 | 21.3 | 7.5 | 4500 | 9.49 |
| 5622 | 225.8 | 20.0 | 7.5 | 4800 | 9.24 |
| 3161 | 218.5 | **19.7** | 7.5 | 4800 | 9.32 |
| 1777 | 220.3 | 20.2 | 7.5 | 4800 | 9.37 |
| **999 (C26)** | 227.2 | 22.0 | 9.0 | 4500 | 9.78 |
| 561 | 237.7 | 21.4 | 7.5 | 4800 | 9.41 |
| 315 | 250.4 | 21.8 | 7.5 | 4800 | 9.38 |
| 177 | 258.3 | 21.7 | 9.0 | 4500 | 9.71 |
| 99 | 273.6 | 21.5 | 7.5 | 4800 | 9.21 |
| 55 | 297.7 | 20.6 | 7.5 | 4800 | 9.01 |
| 30.6 | 357.4 | 21.9 | 7.5 | 4800 | 8.47 |
| 9 (extended grid) | 416.3 | 51.5 | 7.5 | 4300 | 8.81 |
| 2.2 (extended grid) | 1127.5 | 1114.7 | 3.0 | 4000 | 7.72 |

Plotted in Figure 2b. At the C26 broadening the composition axis has
structure: it falls from 357
at He/H = 30 to a minimum of 218 near He/H = 3e3 and rises again, which is the
information the published retrieval uses. **With the broadening free the axis
becomes flat: chi^2 = 19.7-22.0 at every He/H between 30 and 3e4, a spread of
2.3 across three orders of magnitude in composition**, with a rise only below
He/H ~ 9. The published He/H ~ 1e3 sits 2.3 above the flat floor -- an allowed
solution, but no longer a selected one.

**Mass-loss rate.** Every equivalent-width-held solution on that flat ridge
carries log Mdot = 8.5-9.8, i.e. 3e8-6e9 g/s, against the published
2.03e8 (+0.67/-0.58)e8. The minimum is at 2.1e9 g/s, +28 sigma in the units
of their published error bar and a factor 10 above their value. The reason is
visible in stage 1: at the published broadening the model cannot make the
observed area at all, so the fit was never asked to; once it can, it needs
much more absorbing helium.

**Temperature.** The stage-2 minimum sits at 4800 K, -7.2 sigma from the
published 5160 +46/-50 K, and the unconstrained control goes to 4500 K, the
edge of their prior. Direction as expected: at fixed broadening the
temperature is the only knob that widens the line.

### The joint control: the same answer without the equivalent-width constraint

The two-stage construction imposes the equivalent width by hand. The control
does not: it is the same likelihood on a joint grid, profiled, which is also
the authors' own procedure.

**With the broadening fixed, the control reproduces the published
retrieval.** Best chi^2 = 214.7 at Mdot = 1.33e8 g/s, T = 5100 K,
He/H = 561, against the published 2.03e8 (+0.67/-0.58)e8, 5160 +46/-50 K,
He/H = 989 (+/- 0.26-0.29 dex): within about 1 sigma on all three axes, at a
grid resolution of 0.125 dex in Mdot and 200 K in T. Error-rescaled
profile intervals are Mdot 1.24-1.59e8, T 5090-5116 K, He/H 461-796.

**With the broadening free, chi^2 = 18.7 -- an improvement of 196.0 for one
parameter -- and the solution moves the same way the two-stage one does:**

| axis | fixed broadening | free broadening |
|---|---|---|
| v_nt [km/s] | 0 by construction | **8.0**, rescaled 1 sigma [6.8, 9.3] |
| Mdot [g/s] | 1.33e8 [1.24, 1.59]e8 | 5.6e8, rescaled 1 sigma [4.2e8, 6.8e9] |
| T [K] | 5100 [5090, 5116] | 4500 (prior edge), rescaled 1 sigma < 5100 |
| He/H | 561 [461, 796] | 17, bounded only from below: He/H > 17 |
| chi^2 | 214.7 | 18.7 |

and the composition profile is flat over the whole He-rich half of the grid:

| He/H | 31622 | 9999 | 3161 | 999 | 315 | 99 | 30.6 | 16.8 | 9 | 4.6 |
|---|---|---|---|---|---|---|---|---|---|---|
| chi^2, v_nt = 0 | 238 | 230 | 216 | 216 | 216 | 234 | 279 | 322 | 389 | 681 |
| chi^2, v_nt free | 21.9 | 18.8 | 22.6 | 21.2 | 19.0 | 18.9 | 21.5 | 18.7 | 162 | 632 |

At fixed broadening the axis has a genuine minimum near He/H ~ 3e2-3e3 and
rises by 60-170 towards the hydrogen-rich end, which is the shape that
supports a helium-dominated conclusion. With the broadening free the same
axis is flat to within delta chi^2 = 3.9 from He/H = 17 to 3e4, and only the
region below He/H ~ 10 is excluded. The two routes -- equivalent width pinned
by hand, and nothing pinned -- agree: **the broadening lands at 7.5-8 km/s,
the temperature below 5100 K, the mass-loss rate an order of magnitude up,
and the composition ceases to be determined except for a lower bound near
He/H ~ 10-17.**

### The residual coupling between the two stages

The two stages do not separate cleanly, for two reasons, both consequences of
the line being saturated.

First, the equivalent width does not select a unique Mdot: it rises with Mdot,
peaks and falls, so the measured value is generally reached on two branches
(at T = 5160 K, He/H = 999, v_nt = 6 km/s the crossings are log Mdot = 8.2 and
8.8, with chi^2 = 73 and 30). The two-stage procedure keeps whichever fits the
profile better; a first pass that took only the lower branch gave a
systematically worse and composition-dependent stage 2, and that pass is not
part of the record.

Second, along the pinned surface the required Mdot still moves with `v_nt`
(Figure 2c). At the C26 composition and temperature:

| v_nt [km/s] | 0 | 3 | 4.5 | 6 | 7.5 | 10.5 | 18 |
|---|---|---|---|---|---|---|---|
| Mdot [g/s] | 4.09e8 | 3.62e8 | 3.00e8 | 4.51e8 | 6.95e8 | 8.13e8 | 5.48e8 |
| ratio | 1.00 | 0.88 | 0.73 | 1.10 | 1.70 | 1.99 | 1.34 |

(the step between 4.5 and 6 km/s is the branch switch above). **So fixing the
equivalent width constrains the mass-loss rate to about a factor 2 while the
broadening is free**, which is the residual correlation the two-stage
construction was meant to expose. It is much smaller than the factor 10 by
which the *best-fitting* Mdot moves, so that motion is not an artifact of the
coupling: it comes from the equivalent width being reachable at all only once
the line can be made wide.

### The published vector with only the kernel added

A narrower test, with the composition and mass-loss rate held exactly at the
published values and only `v_nt` free:

| v_nt [km/s] | 0 | 3 | 4 | **5** | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|
| chi^2 | 212.0 | 118.8 | 82.6 | **63.6** | 65.7 | 122.4 | 214.4 |
| red EW [A] | 0.00839 | 0.00942 | 0.01007 | 0.01076 | 0.01140 | 0.01237 | 0.01287 |
| FWHM [A] | 0.575 | 0.665 | 0.727 | 0.800 | 0.879 | 1.048 | 1.216 |

Even without touching their parameters, 5 km/s of non-thermal dispersion
improves the fit by delta chi^2 = 148 and brings the equivalent width from
25 % low to 4 % low. The data want the kernel whether or not the other
parameters are allowed to move.

### Consistency with the width memos

The memos state the requirement as **21.96 km/s FWHM (sigma = 9.33 km/s) of
broadening that no 1-D steady-state solution provides**, measured against an
EXHALE line whose intrinsic dispersion is 2.42 km/s.

Here the p-winds line already carries thermal plus Lampon turbulence
(sigma = 4.41 km/s at 5160 K, 4.25 at 4800 K), and the fit adds
**v_nt = 7.5 km/s** on top of it (8.0 km/s in the unconstrained control, with
an error-rescaled 1 sigma interval of 6.8-9.3 km/s), for a total
non-instrumental dispersion of 8.6-9.1 km/s. The memo's requirement implies
9.64 km/s when its 9.33 is added in quadrature to its 2.42 km/s starting
point -- 6-11 % above what is fitted here, and inside the control's interval.
The two numbers are not the same quantity: the memo's 9.33 is what has to be
added to an EXHALE line, this 7.5-8 is what has to be added to a p-winds line
that already carries a larger turbulence term. In the observable, which is
the same quantity on both sides, they agree closely: the stage-2 line has
23.02 km/s FWHM against the 23.29 km/s measured on the same spectrum with the
same three-Gaussian estimator, and against the memo's requirement of
21.96 km/s.

**This is an independent confirmation of the memo's width claim, obtained from
a different code, a different wind structure and a different fitting
procedure**: an independent retrieval, given the freedom, asks for essentially
the amount of extra broadening the memo says is missing, and the free
parameter settles at the same value (7.5 km/s here) from every starting
composition.

It also adds a handle the memos do not use. They argue from the *shape* of the
line; the equivalent width is a second, independent constraint, because in a
saturated line the missing width also costs area. At the published parameters
the model is 25 % short in red-component equivalent width and, at the
published broadening, **cannot** be brought up to the measured area anywhere
in the authors' (T, h_fraction, Mdot) box. The kernel that fixes the width
fixes the area at the same time.

It does not confirm the memo's mechanism, only its size. Nothing here says
what the 7.5 km/s is.
### Caveats

1. **The absolute chi^2 scale is not calibrated.** The best fits reach
   chi^2/dof = 0.16, which says the released errors on individual pixels are
   conservative for a spectrum whose correlated component has already been
   removed by the GP. Rescaling the errors so that the best fit has
   chi^2/dof = 1 divides every delta chi^2 quoted here by 6.3 (6.7 for the
   control): the stage 1 -> stage 2 improvement becomes 31, the joint
   fixed -> free improvement 29, the 5 km/s improvement at the published
   vector 23, and the composition ridge, already flat at 2.3, becomes flat at
   0.4. Every interval quoted as "error-rescaled" already carries this factor.
   The ordering of the results is unaffected; the significance attached to any
   one of them is.
2. **The grids are coarse.** Minima are reported at grid nodes (0.25-0.5 in
   log h_fraction, 200-500 K in T, 1.5 km/s in `v_nt`, 0.125 dex in
   log Mdot), not from a converged optimizer. Treat the quoted parameters as
   locations to about half a grid step. The `v_nt` minimum is the exception:
   it returns 7.5 km/s at essentially every node of both two-stage grids and
   8 km/s in the joint control, so it is robust to where one starts.
3. **The equivalent-width uncertainty** propagates the errors of individual
   pixels only; the GP correlations are not carried, so quoted sigma values on
   equivalent widths are upper bounds on the significance. The percentage
   deficits quoted alongside them do not depend on that.
4. **v_wind is held at 2.26 km/s** throughout, as in the authors' released
   grid. The paper's own MCMC also fits it (2260 +330/-300 m/s), and in
   `'average'` mode it acts as a rigid shift, so it trades against the line
   position rather than the width; it is not expected to absorb `v_nt`, but
   this was not tested here.
5. **The `extended` grid leaves the authors' priors** in temperature (down to
   4000 K) and composition (up to He/H ~ 2), and is labeled as such wherever
   it is used. The stage-2 minimum itself does *not* need that freedom: it
   sits at He/H = 3161 and T = 4800 K, both inside their prior, and the two
   grids return the same solution.
6. **What is fitted is a Gaussian dispersion, not a mechanism.** `v_nt` enters
   the opacity as an isotropic term added in quadrature, exactly like the
   Lampon turbulence it replaces. A real velocity field with structure -- a
   faster outer wind, a stellar-wind interaction region, a non-Gaussian
   distribution -- would not produce the same line even at the same second
   moment. The result says how much velocity spread the data want, not how it
   is arranged.

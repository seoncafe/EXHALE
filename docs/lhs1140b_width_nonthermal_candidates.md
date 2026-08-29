# LHS 1140 b He 10830: can a non-thermal component supply the missing width?

Candidate 3 of the width investigation. The measurement side is settled in
`lhs1140b_width_measurement_audit.md`: the like-for-like requirement is
**21.96 km/s FWHM (sigma = 9.33 km/s)** of broadening that no 1-D
steady-state solution provides, and at most 0.4 km/s of it can be blamed on
the observation. This memo asks whether physics absent from a 1-D wind --
stellar-wind interaction, charge exchange and energetic neutral atoms,
turbulence, magnetic fields, 3-D flow -- can supply the rest.

Everything numerical below was computed here from the released spectrum and
from existing EXHALE outputs; nothing was rebuilt, re-run, or re-converged.
The reference solution is `../LHS1140b/exhale/flux_closure/heh11p1/k01/`
(He/H = 11.1, `He_diffusion: True`, lower-atmosphere profile, Mdot =
2.38e7 g/s measured from `4 pi r^2 rho v` at 10/20/29 R_p). Key numbers were
repeated on four other solutions (`heh10p3/k01`, `hi/k06`,
`heh2p13_diff_kzz1e9`, `heh0p55`) and are quoted with that spread.

---

## Verdict

**The width is a kinematic problem, not an abundance problem, and every
route that accelerates the planet's own gas fails on energy by one to three
orders of magnitude.** The measurement admits essentially no narrow
component at all, which kills the natural form of the stellar-wind /
charge-exchange candidate: whatever moves has to be the *whole* metastable
population, including the dense inner part where an interaction front cannot
reach.

Three results carry the verdict.

1. **The model already has the right amount of metastable helium.** The
   equivalent width of the observed feature is 0.01276 +/- 0.00035 A; the
   EXHALE model line integrates to 0.01227 A, i.e. 1.4 sigma low. The model
   line is 2.9x too deep and 3.0x too narrow, and the product is unity. So
   the entire discrepancy is that the same helium is moving 3x too slowly.
   Nothing is to be gained by adding or removing metastable helium.

2. **The data allow no narrow core.** Fitting the released spectrum with a
   narrow component at the model's intrinsic width (sigma = 2.42 km/s) plus a
   free broad component, the narrow amplitude is consistent with zero and its
   3 sigma upper limit is a blended red depth of **0.32-0.34 %**, i.e. **<= 8 %
   of the line area** (the limit is unchanged if the narrow component is given
   its own Doppler shift). The EXHALE line is 3.60 % deep -- **10x** over that
   limit. Any mechanism that broadens only the outer, interacting gas leaves
   the inner core behind, and the inner core is excluded.

3. **A fast helium atom is not a metastable helium atom.** The 2^3S level is
   made by He+ recombination (the collisional channel from the singlet ground
   state contributes < 1e-11 of the source anywhere in these solutions) and
   the recombination time is 6e5-9e5 s at 2 R_p, against a metastable
   lifetime of 7-38 s there. In gas shocked to 1e5-1e6 K -- what a stellar
   wind at 100-400 km/s produces -- the equilibrium metastable fraction of
   He+ collapses to 1e-7-7e-9, two to four orders below the 1e-5-1e-4 of the
   cool wind.

Candidate by candidate:

| candidate | verdict | the number that decides it |
|---|---|---|
| Stellar-wind interaction + charge-exchange He ENAs | **rejected as the source of the bulk of the width** | the front cannot get inside the absorbing core without contradicting the observed depth, and shocked gas carries no metastables (f <= 2e-7 at 1e5 K) |
| Sustained turbulence at sigma = 9.33 km/s | **rejected on power** | 24-230x the XUV heating deposited in the same shells |
| Faster ordered radial outflow | **rejected on power** | needs the velocity field scaled ~60x; kinetic power scales as k^3, ~4e2x the XUV heating |
| Magnetic (Alfvenic) broadening | **not excluded by field strength; constrained by the power it must be fed, and it cannot reach 1-3 R_p** | v_A = 9.33 km/s needs a 0.34-0.73 G dipole (plasma beta 0.16-0.38); the wave flux is 12-870x the XUV heating, which only a solar-strength or stronger stellar wind could supply |
| 3-D day-night flow, Coriolis, rotation | **scale-excluded** | corotation gives 0.065 km/s at 2 R_p, 0.97 km/s at 30 R_p; a coherent line-of-sight flow is capped at 4.7 km/s by the measured Doppler shift |
| Radiation pressure on the 2^3S level | **too small inside 8 R_p, and the wrong sign** | a_rad = 140 cm/s^2 gives 0.02-1.5 km/s at 2-5 R_p, and what it does give is a blueshift |

**Interpretation** (stated as such): the width is not explained by a 1-D
steady-state wind, and none of the non-thermal mechanisms examined here
supplies it either -- but the arguments that reject them are energy and
atomic-physics scaling arguments applied to *our* density and velocity
structure. They exclude the candidates as modifications of this solution.
They do not exclude the possibility that the absorbing gas is not this
solution's gas at all. The measured pre-ingress feature (depth 1.01 %,
FWHM 0.923 A -- the Supplement prints 0.392 A, which the measurement audit
shows is sigma, not FWHM) is at least as broad as the in-transit feature and
nearly as deep, and the planet occults nothing at that phase. That is direct
evidence that a structure outside any 1-D radial framework carries a large
share of the signal.

**And the problem is not this planet's.** Section 8 places the measurement in
the published population: He 10830 widths that have been fitted freely cluster
at 19-33 km/s FWHM, and the 23.9 km/s measured for LHS 1140 b is an ordinary
member of that set. The anomaly is on the model side, it is shared by an
independently written self-consistent code at the same factor of 2-3
(Taylor et al. 2026, after testing the same axes tested here), and the one
broadening term the field supplies -- `v_turb = sqrt(5kT/3m)`, the monatomic
sound speed -- is Mach 1 by definition and so cannot deliver the Mach 2-4 the
data ask for. Three published analyses say as much in their own systems.

---

## 1. Restating the requirement: the column is right, the velocity is not

**Finding.** The observed and modeled lines contain the same number of
metastable helium atoms to within 4 %. The discrepancy is entirely in how
those atoms are distributed in velocity.

Equivalent widths, both measured with the same integration over the same
triplet, and converted to a disk-averaged column with
`int sigma dlambda = 0.02654 f lambda^2/c` and f_total = 0.5392:

| | EW [A] | disk-averaged N(2^3S) [cm^-2], optically thin |
|---|---|---|
| observation (three-Gaussian fit to the released spectrum) | 0.01276 +/- 0.00035 | 2.28e10 |
| observation (direct trapezoid over the released points) | 0.01304 | 2.33e10 |
| EXHALE `heh11p1/k01`, R = 68,000 | 0.01227 | 2.19e10 |
| EXHALE `heh10p3/k01`, R = 68,000 | 0.01195 | 2.13e10 |

Two accuracy notes on the right-hand column, neither of which changes the
comparison because the same conversion is applied to both rows and the two
lines are similarly saturated (measured red-to-blue depth ratios 7.04 in the
model and 6.35 in the data, against the optically thin 8). First, that
saturation means a thin-limit conversion *understates* the true column, by
14 % and 23 % respectively. Second, integrating the solution's `n(2^3S)`
profile geometrically over the stellar disk,
`int 2 pi b N(b) db / (pi R_star^2)` for b <= R_star = 13.63 R_p, gives
1.90e10 cm^-2, 13 % below what the model's own equivalent width implies --
the opposite sign to the saturation correction, and not traced further here.
Take the entries as accurate to about 20 %.

Consistency of depth, width and area:

```
depth ratio  model/obs = 3.600 / 1.254 = 2.87
width ratio  model/obs = 0.2825 / 0.841 = 0.336
product                                 = 0.964
```

**So there is exactly one thing wrong with the model line.** Every candidate
below is therefore judged on whether it moves helium, not on whether it makes
helium.

This agreement is not a fit. The `flux_closure` runs converge the elemental
flux between the photochemical lower atmosphere and the wind
(`heh11p1/closure.json`: "the converged He/H is an output, not this number");
nothing in them is tuned to the He 10830 measurement, and indeed their line
depth misses it by 2.9x. That the equivalent width lands within 4 % is an
independent result.

### 1.1 No narrow core is allowed

Fit to `Cherubim_2026/.../LHS1140b_He10833_Fig3B_reproduced.csv`
(`absorption_depth_positive_percent`, 129 points), three Gaussians at the
vacuum triplet wavelengths in the optically thin 1:3:5 ratio, each component
convolved with the R = 68,000 kernel:

| narrow sigma | narrow blended depth | 3 sigma upper limit | narrow share of line area (3 sigma UL) | broad sigma | chi^2 |
|---|---|---|---|---|---|
| 2.42 km/s (EXHALE intrinsic), shared shift | 0.015 % | 0.339 % | 0.082 | 9.79 km/s | 18.8 / 125 |
| 2.42 km/s, free shift | 0.040 % | 0.317 % | 0.078 | 9.80 km/s | 18.6 / 124 |
| 4.15 km/s, free shift | 0.067 % | 0.470 % | 0.153 | 9.93 km/s | 18.6 / 124 |
| 5.54 km/s, free shift | 0.140 % | 0.810 % | 0.291 | 10.17 km/s | 18.5 / 124 |

Repeating on the pre-GP column (`raw_excess_absorption_percent`) gives a
3 sigma limit of 0.370 % and a narrow area share of 0.093, so the result does
not depend on the correlated-noise correction. Repeating with the paper's own
parameterization for the broad component -- three free amplitudes instead of
the fixed 1:3:5 -- gives 0.338 %, so it does not depend on the optically thin
amplitude ratio either. A single broad Gaussian already fits the data on its
own (chi^2 = 18.8 over 126 dof; that this is far below the dof means the
released uncertainties are conservative, so the limits above are, if
anything, loose).

**Fact.** The model line at 3.60 % blended depth exceeds the 3 sigma narrow-core
limit by a factor of **9.7-10.6**.

### 1.2 Where the model puts the absorbing helium

Fraction of the optically thin absorption (weight `2 pi b N(b)`) inside an
impact parameter b:

| | 25 % inside | 50 % inside | 75 % inside | 90 % inside |
|---|---|---|---|---|
| EXHALE `heh11p1/k01` | 2.54 R_p | 4.44 R_p | 7.39 R_p | 11.95 R_p |
| EXHALE, four-run spread | -- | 4.39-5.26 R_p | -- | 11.8-13.8 R_p |
| p-winds, the paper's fiducial (`pwinds_oracle/profile_published_gj1132.txt`) | 1.28 R_p | 1.43 R_p | 1.62 R_p | 1.81 R_p |

Two things follow. The EXHALE distribution is far more extended than the
paper's p-winds model, which concentrates 90 % of its absorption inside
1.81 R_p and reaches its width by saturation instead (red/blue = 2.36 against
the thin limit of 8). And a third of the EXHALE absorption comes from inside
3 R_p, where the gas is collisional -- Kn(H I) = 0.0010, 0.016 and 0.053 at
1.2, 2.0 and 3.0 R_p, rising through 0.19 at 5 R_p to 0.53 at 8 R_p
(`src/utils/collisional_validity.py` on this run, `_adv` pair). The numbers
first written here -- 0.0008, 0.019, 0.076, 0.29, 0.77 -- were a separate
recomputation using `L = min(H_p, r)`, the structure scale that memo carried
before 2026-08-28; `collisional_validity.md` section 2 now builds `L` from the
density, temperature and velocity gradients instead, and the conclusion of
this paragraph is the same under either.

### 1.3 What the broad component has to be

Putting 1.1 and 1.2 together, the broad component that has to carry the line
is specified as follows.

* **Velocity width.** sigma = 9.6-10.2 km/s intrinsic across the fit variants
  on the GP-corrected spectrum (8.6 km/s on the pre-GP column), i.e.
  FWHM 20-24 km/s. This is not a free parameter -- it *is* the measurement.
* **Share.** >= 92 % of the line area at 3 sigma, and the best fits put it at
  99.6 %.
* **Column.** The projected number of metastable atoms is fixed by the
  equivalent width, not by the geometry: 1.65e31 in the optically thin limit,
  2.15e31 after the 23 % saturation correction implied by the measured
  red/blue ratio of 6.35. Splitting that between covering fraction f_cov of
  the stellar disk and local column N:

  | f_cov | equivalent radius | N(2^3S) required |
  |---|---|---|
  | 1.0 | 13.6 R_p | 3.0e10 cm^-2 |
  | 0.30 | 7.5 R_p | 1.0e11 cm^-2 |
  | 0.10 | 4.3 R_p | 3.0e11 cm^-2 |
  | 0.030 | 2.4 R_p | 1.0e12 cm^-2 |
  | 0.010 | 1.4 R_p | 3.0e12 cm^-2 |

  The depth adds one more relation, `depth = f_cov (1 - exp(-tau_red))`: with
  the 1.24 % excess depth and tau_red = 0.557 from the red/blue ratio, f_cov
  = 0.029, i.e. an absorbing region of effective radius **2.3 R_p** carrying
  **1.0e12 cm^-2**. **Caveat:** the paper's red/blue ratio is 6.7 +12.7/-3.1,
  so tau is barely constrained; at red/blue = 4 the same construction gives
  f_cov = 0.015 and N = 3.6e12 cm^-2. The column is therefore known to a
  factor of a few, the covering fraction to a factor of two, and only the
  product -- the equivalent width -- is measured well.

For comparison, the EXHALE solution's tangential metastable column peaks at
3.6e11 cm^-2 at b = 1.2 R_p and is 1.4e11 at 2 R_p: **an order of magnitude
below what a compact broad component would need, spread over an area an order
of magnitude larger.** That is the same statement as "right area, wrong
distribution", now in column-density units.

The paper's p-winds fiducial sits at the other end: 4.15e12 cm^-2 at
b = 1 R_p with 90 % of its absorption inside 1.81 R_p, i.e. in the compact,
saturated regime, which is how it reaches 0.458-0.575 A (measured in
`lhs1140b_width_measurement_audit.md`) without any velocity field to speak of (its outflow is 0.02 km/s at 2 R_p). It gets a third of the
way to the observed width by saturation, and its red/blue ratio of 2.36
against the measured 6.35 is the price.

### 1.4 The velocity scales available

```
v_esc at 1 / 2 / 5 / 10 R_p   =  20.13 / 14.23 / 9.00 / 6.36 km/s
sound speed at 1.5-5 R_p       =  4.0-5.3 km/s   (mean mass 2.6-3.5 amu)
wind speed in the solution     =  0.006 km/s at 1.5 R_p, 0.27 at 6, 1.37 at 29
required sigma                 =  9.33 km/s   (= v_esc at 4.7 R_p)
required FWHM                  =  21.96 km/s  (= 1.09 v_esc at the surface)
```

**Interpretation.** The requirement is the planet's escape speed. That is
suggestive -- gas on ballistic or barely-bound trajectories would have this
dispersion naturally -- and it is also why the energy arguments below are so
tight: delivering the escape speed to the whole absorbing column is, by
definition, the entire energy-limited mass-loss budget.

---

## 2. Stellar-wind interaction: where the front sits

**Finding.** For any stellar wind strong enough to matter, the contact
surface stands off *outside* the bulk of the metastable column, so the
interaction can only broaden the outer gas -- which section 1.1 forbids.
Pushing the front inward far enough to sweep up the core requires a wind
mass-loss rate of order tens of solar, and is separately excluded by the
observed transit depth.

Pressure balance uses the solution's own total pressure `n_tot k T + rho v^2`
(thermal dominates everywhere: the wind is subsonic through the whole
absorbing region, and reaches Mach 1 only in the last cell of the domain --
max Mach 1.003 at 29.05 R_p for `heh11p1/k01` and 1.005 at 29.52 R_p for
`heh10p3/k01`, both computed with the file's own `sqrt(gamma p / rho)`;
`hi/k06` peaks at 0.550) against `p_ram = Mdot_* v_sw / (4 pi a^2)` with
`4 pi a^2 = 2.517e25 cm^2` at a = 0.0946 AU. A solar-strength wind
(Mdot_sun = 1.26e12 g/s, 400 km/s) gives p_ram = 2.00e-6 dyn/cm^2 at this
orbit.

| p_ram [dyn/cm^2] | Mdot_* at 300 km/s [g/s] | in Mdot_sun | standoff [R_p] | share of the model's He 10830 absorption outside it |
|---|---|---|---|---|
| 1e-4 | 8.4e13 | 67 | 1.65 | 0.87 |
| 1e-5 | 8.4e12 | 6.7 | 2.39 | 0.76 |
| 1e-6 | 8.4e11 | 0.67 | 3.78 | 0.58 |
| 1e-7 | 8.4e10 | 0.067 | 6.51 | 0.30 |
| 1e-8 | 8.4e9 | 0.0067 | 13.08 | 0.082 |
| 1e-9 | 8.4e8 | 0.00067 | 26.99 | 0.002 |

**Two hard constraints cross here.**

*From the depth.* The observed excess depth of 1.24 % with the zero-continuum
convention requires an equivalent opaque radius of 1.82 R_p
(`R_eff^2 = R_p^2 + delta R_star^2`). Absorbing material must therefore reach
1.82 R_p, so the contact surface cannot be inside it. That caps p_ram at
about 6e-5 dyn/cm^2.

*From the width.* To leave less than 8 % of the line area in a narrow core,
the front would have to stand at roughly 1.4-1.5 R_p -- inside the radius the
depth demands. **The two requirements are incompatible**: there is no standoff
distance at which a stellar-wind-driven broadening both reaches deep enough to
avoid a narrow core and stays shallow enough to produce the observed depth.

**Caveat, stated as a limit of this argument.** The share-outside column uses
*our* metastable distribution. A solution whose metastable helium were much
more extended would relax it. The p-winds fiducial goes the other way (90 %
inside 1.81 R_p), so among the models on the table this is the favorable case,
not the unfavorable one.

---

## 3. The decisive constraint: making metastable helium where the gas is fast

**Finding.** He 2^3S is a locally maintained population with a lifetime of
seconds to hours, fed only by He+ recombination on a timescale of 1e5-1e9 s.
The population therefore reports the conditions where it sits, and shocked or
entrained gas -- hot, ionized, and passing through quickly -- cannot carry it.

### 3.1 Source and sink, measured in the solution

The 2^3S balance is `ion_residual_core.f90:72-86` (`tr_triplet_row`): source
`n_e (n_HeII alpha_23S + n_HeI q13)`, sinks radiative decay `A31`, Penning
ionization with H `n_HI Q31`, electron collisional transfer `q31a + q31b`,
electron-impact ionization `b_heiTR`, and photoionization `g_heiTR`. The rate
coefficients are `Cool_coeff.f90`: `rec_HeII_23S` = 2.10e-13 (T/1e4)^-0.778,
`penning_HeI_23S` = 1.9e-9 (300/T)^0.07 below 4000 K and 9.1e-9 (300/T)^0.50
above (Taylor et al. 2025, ApJ 989, 68, Table 2), `A31` = 1.272e-4 s^-1.

The collisional source from the singlet ground state is negligible: its share
of the total source is 1e-13 at 1.2 R_p and falls to 1e-33 by 5 R_p. **The
metastable population is a recombination product, without qualification.**

| r [R_p] | T [K] | n_e [cm^-3] | t_rec = 1/(n_e alpha_23S) [s] | t_Penning = 1/(n_HI Q31) [s] | t_life (total) [s] | t_rec / t_Penning |
|---|---|---|---|---|---|---|
| 1.20 | 5125 | 3.4e7 | 8.3e4 | 1.1 | 0.68 | 7.7e4 |
| 1.50 | 5979 | 1.2e7 | 2.6e5 | 7.8 | 2.6 | 3.3e4 |
| 2.00 | 5223 | 3.3e6 | 8.8e5 | 37.7 | 12.3 | 2.3e4 |
| 3.02 | 3812 | 6.4e5 | 3.5e6 | 337 | 112 | 1.0e4 |
| 4.98 | 2643 | 1.1e5 | 1.6e7 | 2340 | 1084 | 6.7e3 |
| 7.96 | 2078 | 2.4e4 | 5.8e7 | 1.2e4 | 4024 | 5.0e3 |
| 20.12 | 555 | 1.9e3 | 2.7e8 | 1.7e5 | 6559 | 1.6e3 |

Across the other four solutions, `t_rec(2 R_p)` = 5.9e5-8.8e5 s and
`t_Penning(2 R_p)` = 7.3-37.7 s, so the ratio is 2.3e4-8.0e4 everywhere.

The metastable fraction is the ratio of these,
`n(2^3S)/n(He+) = alpha_23S n_e t_life`, measured as **8e-6 to 7e-5** across
1.2-8 R_p, and `n(2^3S)/n(He)` = 5e-8 to 1.4e-5. These are read off the
advection-corrected profiles, the ones the transit calculation consumes; the
equilibrium profiles of the same run reach 8e-4 and 5e-4 respectively around
8 R_p, so the local-equilibrium ceiling is about an order of magnitude above
what the model line actually carries.

**How far a metastable atom can travel while it is still metastable.**
Multiplying the total lifetime by the required dispersion,
`L = 9.33 km/s x t_life`, gives the distance a fast metastable atom covers
before it is destroyed:

| r [R_p] | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 8.0 | 12 | 20 | 29 |
|---|---|---|---|---|---|---|---|---|---|
| t_life [s] | 1.09 | 7.81 | 37.5 | 324 | 1.80e3 | 4.69e3 | 6.59e3 | 7.51e3 | 7.73e3 |
| L [R_p] | 0.0009 | 0.0066 | 0.032 | 0.27 | 1.53 | 3.97 | 5.57 | 6.35 | 6.55 |

Inside 3 R_p, where a third of the absorption sits, an atom moving at the
required dispersion covers less than 0.3 R_p before Penning ionization
removes it. **The inner metastable population cannot have been made
anywhere else.** Beyond about 5 R_p the constraint relaxes, because Penning
ionization becomes slower than radiative decay and the lifetime saturates at
`1/A31` = 7.86e3 s -- but that is the outer gas that section 1.1 has already
shown cannot carry the line on its own.

### 3.2 The same fraction in shocked gas

A stellar wind at v_sw shocked at the contact surface reaches
`T = (3/16) mu m_H v_sw^2 / k` = 3.4e4 K at 50 km/s, 1.4e5 K at 100 km/s and
1.2e6 K at 300 km/s. Evaluating the same equilibrium there, with electron
collisional ionization of the metastable (`ci_HeI23S`, threshold 4.8 eV) as
the dominant sink and Penning switched off because the gas is ionized:

| T [K] | n_e = 1e3 | 1e4 | 1e5 |
|---|---|---|---|
| 1e4 | 1.6e-6 | 1.6e-5 | 1.3e-4 |
| 3e4 | 6.0e-7 | 2.5e-6 | 3.7e-6 |
| 1e5 | 1.3e-7 | 2.1e-7 | 2.3e-7 |
| 3e5 | 2.9e-8 | 3.8e-8 | 3.9e-8 |
| 1e6 | 6.4e-9 | 7.3e-9 | 7.4e-9 |

(entries are `n(2^3S)/n(He+)`; the photoionization sink is omitted, which
makes these upper limits).

**Fact.** Above 1e5 K the metastable fraction is 1e-7 or below -- 100 to
1000 times smaller than the cool wind achieves. **The shocked layer is not a
metastable helium reservoir.** Metastable helium can only exist in the cool,
partly ionized planetary gas, which is exactly the gas the interaction has
not yet touched.

### 3.3 Charge exchange: which channels could even apply

Stated from energetics computed here; the cross sections themselves are a gap
recorded in "what was not determined".

* `H+ + H -> H + H+` is resonant and is the channel behind every published
  hydrogen ENA calculation. It has no helium analogue at the same strength.
* `He+ + H(1s) -> He(1^1S) + H+` is exothermic by
  24.587 - 13.598 = 10.99 eV, i.e. a large energy defect, so it is
  non-resonant and slow at low collision energies. It produces **ground-state**
  helium, which is invisible at 10830 A.
* `He+ + H(1s) -> He(2^3S) + H+` is the only single-step channel that would
  make a fast *metastable* atom. The 2^3S level lies 19.82 eV above the
  ground state, so the transferred electron ends up bound by
  24.587 - 19.82 = 4.77 eV and the reaction is **endothermic by 8.83 eV**.
  With mu(He,H) = 0.805 amu the threshold relative speed is **46.0 km/s**.
  It is therefore closed in the planetary wind (relative speeds < 2 km/s) and
  open only in the stellar-wind interaction -- which is the region section 3.2
  has just shown cannot hold metastables anyway.
* `He+ + He -> He + He+` is resonant and fast, but again yields ground-state
  helium.
* `He++ + H -> He+ + H+` is quasi-resonant and fast, and produces an ion, not
  a neutral.

**Interpretation.** Even granting an efficient helium ENA population, the
product of every channel except the endothermic one is ground-state helium.
Turning it into 2^3S requires ionizing it again and recombining, i.e.
`t_rec` = 1e5-1e9 s, against a residence time in the interaction region of
`~5 R_p / 30 km/s ~ 2e3 s`. The gap is four to six orders of magnitude.

### 3.4 What *is* possible: the metastable inherits the ion's motion

The one route that survives this section is the one Taylor et al. (2026, ApJ
999, 214, Section 4) invoke for HD 189733 b: since the metastable is made by
recombination, it is born with whatever velocity the recombining He+ had.
This memo confirms that the inheritance is not erased by collisions. The
number of momentum-transfer collisions a metastable atom suffers within its
own lifetime, computed with the same rigid-sphere model, is

```
r/R_p     1.1   1.2   1.5   2.0   3.0   5.0   8.0   12    20
N_coll   11.2   5.5   2.0   1.3   1.4   1.6   1.0   0.4   0.1
```

so a metastable born moving does not thermalize before it absorbs.
**Thermalization is not the obstacle; producing a fast He+ population is.**
And He+ ions in this region are themselves tightly collisional -- with an
assumed resonant He+ + He cross section of 1e-14 cm^2, the ion-neutral
coupling time at 2 R_p is of order 0.01 s against the 12 s metastable
lifetime, i.e. some three orders of magnitude shorter -- so a suprathermal
*ion* population cannot be sustained either unless it is continuously
driven. That converts the question
into an energy question, which is section 4.

---

## 4. Energy: what it costs to move the absorbing gas at 9.33 km/s

**Finding.** Every route fails by 1-3 orders of magnitude. The reference is
the solution's own volumetric heating, integrated over 1-30 R_p:
**1.224e20 erg/s** (ghost cells excluded). As a check that the solution
absorbs what it should: integrating the SED file the run actually reads over
10-1300 A gives F_XUV = 33.4 erg/s/cm^2 at the orbit -- matching the
0.033 W m-2 the paper quotes as its fiducial -- and `F_XUV pi R_p^2` =
1.28e20 erg/s.

### 4.1 Sustained turbulence

Dissipation of a cascade with outer scale L = r, `eps = rho v^3 / (2 L)`,
integrated over shells, against the XUV heating deposited in the same shells:

| shell | gas mass [g] | turbulent kinetic energy (3-D, sigma = 9.33 km/s) [erg] | dissipation [erg/s] | XUV heating in shell [erg/s] | ratio |
|---|---|---|---|---|---|
| 1.2-30 R_p | 7.3e13 | 9.6e25 | 1.62e22 | 6.9e19 | 234 |
| 1.5-30 R_p | 3.1e13 | 4.0e25 | 4.10e21 | 3.2e19 | 128 |
| 2-30 R_p | 2.0e13 | 2.6e25 | 1.68e21 | 2.0e19 | 86 |
| 5-30 R_p | 9.6e12 | 1.2e25 | 3.32e20 | 8.9e18 | 37 |
| 10-30 R_p | 5.8e12 | 7.5e24 | 1.28e20 | 5.2e18 | 24 |

An outer scale ten times larger would divide these by ten and still leave
2-23. The specific energy alone is 5.2x thermal at 2 R_p, i.e. a 3-D turbulent
Mach number of 3.1.

**This is Murray-Clay, Chiang & Murray (2009, ApJ 693, 23, Section 4, p. 38)
applied to this system.** Their statement about Garcia Munoz's turbulent
broadening of Ly-alpha reads:

> "But to generate the large velocities observed, the energy in such
> turbulence would need to exceed the thermal and bulk kinetic energies in the
> mean flow by a factor of ∼100. Such energy requirements seem
> insurmountable."

Here the factor in *power* is 24-234 and in *specific energy* 5.2. The
objection is weaker than for Ly-alpha at 100 km/s but still decisive.

### 4.2 A faster ordered radial outflow

Recomputing the optically thin line profile from this solution with the
radial velocity field scaled by k (the same triplet, the same thermal widths,
the same R = 68,000 kernel, integrated over impact parameter and along the
line of sight):

| k | max radial speed [km/s] | fitted FWHM [A] | FWHM [km/s] |
|---|---|---|---|
| 1 | 1.37 | 0.199 | 5.5 |
| 5 | 6.87 | 0.234 | 6.5 |
| 10 | 13.7 | 0.296 | 8.2 |
| 20 | 27.5 | 0.412 | 11.4 |
| 40 | 55.0 | 0.631 | 17.5 |

(this thin-limit reconstruction gives 0.199 A at k = 1 where the full transit
tool gives 0.2825 A, the difference being saturation and rotation, so the
column is a lower bound on the width; the trend is what matters).

A radial field is an inefficient broadener because the gas at closest
approach, which dominates the column, has zero line-of-sight velocity. It
does not, however, produce a distinctively non-Gaussian profile: the k = 40
curve is still fitted by a Gaussian triplet to within 6.5 % of its peak, so
the line *shape* adds no constraint beyond the width.
Reaching the observed 0.84 A needs k ~ 60-80, i.e. peak speeds of 80-110 km/s,
four to five times the surface escape speed. Scaling the velocity at fixed
density scales the mass flux with it, so the kinetic power goes as k^3:
`0.5 Mdot v^2 = 2.25e17 k^3 erg/s`, which equals the 1.22e20 erg/s heating at
**k = 8** (FWHM 7.7 km/s) and would need **4e2 x** the heating at k = 60.

**Fact.** The absolute ceiling on an ordered radial outflow from this
solution is FWHM ~ 8 km/s. The requirement is 22 km/s.

### 4.3 The measured Doppler shift caps any coherent line-of-sight flow

During transit the star-planet axis is the line of sight, so a day-to-night
flow, or radiation-driven anti-stellar acceleration, projects onto it
directly and shifts the line rather than broadening it. From the released
MCMC chain the shared shift is +0.072 A (+1.99 km/s) with percentiles

```
   0.135 %      2.5 %       16 %       50 %       84 %      97.5 %
  -4.71      -2.07      -0.03      +1.99      +4.14      +7.65   km/s
```

so a coherent blueshifted flow of the whole absorbing column is limited to
**<= 4.7 km/s** at the 3 sigma-equivalent level. That is a direct observational
cap on the day-night and radiation-pressure routes, independent of any model.

---

## 5. Magnetic fields

**Finding.** The field strength required to make the Alfven speed equal the
required dispersion is *not* implausible -- a few tenths of a gauss -- but
such a field would be dynamically dominant (plasma beta < 0.4) and would
restructure the wind rather than decorate it, and the wave energy flux it
would have to carry is the same 1-2 orders of magnitude over budget as
section 4.1.

`v_A = B / sqrt(4 pi rho) = 9.33 km/s` requires:

| r [R_p] | rho [g/cm^3] | B required [G] | equivalent dipole surface field [G] | plasma beta at that B | nu_ni / omega |
|---|---|---|---|---|---|
| 1.20 | 1.63e-14 | 0.423 | 0.73 | 0.28 | 96 |
| 1.50 | 1.27e-15 | 0.118 | 0.40 | 0.38 | 44 |
| 2.00 | 1.66e-16 | 0.0426 | 0.34 | 0.38 | 15 |
| 3.02 | 2.24e-17 | 0.0157 | 0.43 | 0.30 | 4.5 |
| 4.98 | 3.12e-18 | 0.0058 | 0.72 | 0.22 | 1.3 |
| 7.96 | 6.22e-19 | 0.0026 | 1.31 | 0.17 | 0.45 |

`nu_ni / omega` compares the neutral-ion collision frequency
(`n_ion x 2e-9 cm^3/s`) with a wave frequency `v_A / r`: neutrals are carried
by the waves inside about 5 R_p and decouple beyond 8 R_p, so the outer
metastable helium would not share the motion even if the inner gas did.

The energy flux such waves must carry, `rho (3 sigma^2) v_A` through a sphere:

```
r = 1.5 / 2 / 3 / 5 / 8 R_p  ->  1.1e23 / 2.5e22 / 7.6e21 / 2.9e21 / 1.5e21 erg/s
                             ->  872 / 203 / 62 / 24 / 12  x the XUV heating
```

For comparison, the stellar wind's kinetic power intercepted by an obstacle
of radius R_ob is `0.5 p_ram v_sw pi R_ob^2`: at p_ram = 1e-6 dyn/cm^2 and
300 km/s this is 2.3e20 erg/s for R_ob = 2 R_p (1.9x the XUV heating) and
9.2e20 for R_ob = 4 R_p (7.5x). **So a star-planet interaction could in
principle power a wave field of the required strength at 5-8 R_p, but only
with a stellar wind at least of order solar strength at this orbit, and it
still cannot reach the inner 1-3 R_p where a third of the absorption sits.**

Taylor et al. (2026), Section 4, propose exactly this route for HD 189733 b,
and are explicit that it is speculative:

> "A particularly interesting idea is that HD189733b hosts a population of
> trapped ions on low-latitude magnetic field lines that form an energetic,
> turbulent plasma population around the planet. During periods of higher
> stellar activity, the turbulence is strengthened, similarly to the Earth's
> magnetosphere during high solar activity. ... While these ideas are
> speculative, we encourage further exploration of plasma effects and magnetic
> field interactions to characterize He I 10830 Å transit depths and light
> curves."

Their host is an active K dwarf with an inferred surface polar dipole of
20 +/- 7 G on the planet. LHS 1140 is a 131 d rotator with
L_X/L_bol = 6.5e-6, i.e. the opposite end of the activity range, which makes
the magnetically driven route *less* attractive here than there, not more.

---

## 6. Rotation, Coriolis and 3-D flow

**Finding.** Scale-excluded by an order of magnitude, with the added
constraint of section 4.3.

With P = 24.73723 d and tidal locking, Omega = 2.940e-6 s^-1:

```
corotation speed   0.065 km/s at 2 R_p ; 0.324 at 10 R_p ; 0.972 at 30 R_p
Hill radius        4.417e10 cm = 40.1 R_p   (the 30 R_p domain is inside it)
orbital speed      41.6 km/s
Rossby number      u / (2 Omega L) ~ 0.19 at 2 R_p for u = 0.024 km/s
```

The Rossby number below unity means the flow *is* rotationally influenced --
Coriolis will bend streamlines -- but the velocity that rotation itself can
impart is 0.07-1.0 km/s, two orders below the requirement. A pressure-driven
day-night flow is bounded above by the thermal energy available, i.e. by a
few times the 4-5 km/s sound speed in the most extreme complete-conversion
limit, and by section 4.3's 4.7 km/s cap on any coherent line-of-sight
component. Taylor et al. (2026), Section 4, reach the same conclusion for
HD 189733 b: "typical GCM wind speeds (a few kilometers per second) are
unlikely to supply the additional ∼12 km s−1 required to match the observed
width."

---

## 7. Radiation pressure on the metastable level

Not among the mechanisms usually considered for this line, checked here
because the host is an M dwarf whose flux peaks near 1 micron and the
planet's surface gravity is only 1837 cm/s^2.

From the SED actually used by these runs
(`../LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt`, 49.05 erg/s/cm^2/A at
10830 A), the radiation force on an atom in 2^3S is
`(pi e^2 / m_e c) f_total F_nu / c` = 9.34e-22 dyn, i.e.

```
a_rad = 140.5 cm/s^2
g_planet = 1837 / 459 / 204 / 74 / 29 cm/s^2 at 1 / 2 / 3 / 5 / 8 R_p
```

so **radiation pressure on the metastable level exceeds the planet's gravity
beyond about 3.5 R_p.** What it can actually do is limited by the metastable
lifetime, since a destroyed atom returns its momentum to the bulk (the bulk
acceleration is `a_rad x f_ms` ~ 2e-3 cm/s^2, negligible):

```
velocity gained in one metastable lifetime:
  2 R_p  0.017 km/s ;  3 R_p  0.16 ;  5 R_p  1.5 ;  8 R_p  5.7 ;  12-20 R_p  ~9
```

**Verdict: too small where the absorption is, and the wrong sign where it is
not.** It cannot touch the inner 1-3 R_p; in the outer region (14 % of the
absorption beyond 10 R_p) it is of the required magnitude, but it drives a
systematic blueshift of the metastable population relative to the bulk gas,
and section 4.3 caps any such coherent shift at 4.7 km/s. It is nonetheless a
real effect absent from `EXHALE_transit.py`, which Doppler-shifts the
absorption only by the hydrodynamic `v_LOS` (`EXHALE_transit.py:536-537`,
`v_x = x_LOS*v_LOS/r_LOS`).

---

## 8. Precedent: the width in the wider He 10830 literature

**Finding.** The measured LHS 1140 b width is ordinary. Published He 10830
transit widths cluster at 19-33 km/s FWHM, and the 23.9 km/s measured here
sits in the middle of that range. What is unusual is not the observation but
the model: a cold, weakly escaping super-Earth is being asked for a line as
broad as the ones hot Jupiters and hot Neptunes show, and the 1-D
steady-state solution delivers a quarter of it. The same shortfall, at the
same factor of 2-3, is reported by an independent code for HD 189733 b after
testing the same axes tested here.

Provenance: the entries below were taken from the published journal versions
and the quotations are reproduced verbatim from them. Of these, the Taylor
et al. (2026) quotations and the p-winds source-code statement in section 8.2
were re-checked directly against `../../references/Taylor_2026_ApJ_999_214.pdf`
and `../../p-winds/p_winds/transit.py` while writing this section; the rest
carry over from the census compiled for this memo and were not re-opened here.

### 8.1 The population of measured widths

Only a few analyses fit the He 10830 width freely rather than predicting it.
Where they do, the answer is 19-33 km/s FWHM. Wavelengths are converted with
`FWHM[km/s] = 2.998e5 x FWHM[A] / 10833`.

| system | instrument (R) | published width | FWHM [km/s] | source |
|---|---|---|---|---|
| **LHS 1140 b** | WINERED (68,000) | 0.865 +0.261/-0.155 A | **23.9 +7.2/-4.3** | this memo, on the released spectrum |
| WASP-107 b | GIANO-B (~50,000) | 0.1176 nm | 32.6 +4.2/-3.6 | Guilluy et al. 2024, A&A 686, A83, Table A.3 |
| HAT-P-11 b (night 1/2/3) | GIANO-B | 0.1059 / 0.0686 / 0.0829 nm | 29.3 / 19.0 / 22.9 | Guilluy et al. 2024, Table A.3 |
| WASP-69 b | GIANO-B | 0.0960 nm | 26.6 | Guilluy et al. 2024, Table A.3 |
| GJ 3470 b | GIANO-B | 0.2180 nm | 60.3 | Guilluy et al. 2024, Table A.3 -- see caveat |
| four young sub-Neptunes | Keck/NIRSPEC (32,000) | ~1.1 A, "within 10 %" for all four | 27-33 | Zhang et al. 2023, AJ 165, 62 |
| SPIRou detections (typical) | SPIRou | ~0.08 nm | ~22 | Masson et al. 2024, A&A 688, A179 |
| HD 189733 b | CARMENES | free `b` = 15.1 +/- 0.6 km/s | (b, not FWHM) | Salz et al. 2018, A&A 620, A97, Table 3 |
| HAT-P-32 b | CARMENES | free `vt,a` = 15.2 +/- 3.1 km/s | (dispersion, not FWHM) | Czesla et al. 2022, A&A 657, A6 |

The GIANO-B instrumental FWHM is 6.0 km/s and the NIRSPEC one about 9 km/s,
so deconvolving would lower the entries by a few per cent to about 10 %; the
table quotes them as published. The GJ 3470 b entry should **not** be taken
at face value -- Guilluy et al. write that "The larger FWHM for GJ 3470 b is
believed to be attributed to residuals from OH line contamination". Ninan et
al. (2020, ApJ 894, 97) separately report a blue wing on the same planet
reaching -36 km/s, "the largest reported in the literature to date".

Zhang et al. (2023) use the width as *evidence for* photoevaporation rather
than as a problem, and their reasoning is the one that fails here:

> "The width is the most consistent metric and is within 10% of 1.1 A for all
> planets"

> "For a Parker outflow, the sound, bulk outflow, and thermal speeds are all
> close to sqrt(2kB T/mu mH) = 9.4 km s-1 for T = 7000 K and mu = 1.3, leading
> to an FWHM of 2.355 sigma ~ 22 km s-1."

That argument needs a 7000 K, hydrogen-dominated, transonic outflow. The
LHS 1140 b solution has 5200 K at 2 R_p, a mean mass of 2.6-3.5 amu because
of the helium enrichment, and a wind that reaches 0.024 km/s there
(section 1.4). The same construction applied to this solution gives the
model's 2.4 km/s intrinsic dispersion, not 9.4.

### 8.2 The literature's "turbulence" term is Mach 1 by construction

There is only one turbulence formula in this literature:

```
v_turb = sqrt(5 k T / 3 m)          Lampon et al. 2020, A&A 636, A13, Eq. (16)
```

with `m` the mass of a helium atom -- i.e. the adiabatic sound speed of a
monatomic gas, a Mach number of exactly 1. p-winds adopts it unchanged as its
Eq. (21) (Dos Santos et al. 2022, A&A 659, A62). It is a fixed function of the
temperature, not a free parameter, and at realistic outflow temperatures it
returns only 3.2-7.0 km/s.

Evaluated on this solution it returns

| r [R_p] | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 8.0 | 12 |
|---|---|---|---|---|---|---|---|
| T [K] | 5125 | 5979 | 5223 | 3812 | 2643 | 2078 | 1891 |
| v_turb [km/s] | 4.20 | 4.53 | 4.24 | 3.62 | 3.01 | 2.67 | 2.55 |
| 9.33 / v_turb | 2.22 | 2.06 | 2.20 | 2.58 | 3.10 | 3.49 | 3.66 |

with a metastable-weighted mean of 3.8-4.0 km/s over 1-8 R_p, i.e. a shortfall
factor of **2.3-2.5** on the absorber-weighted average and 2.1-3.7 across the
absorbing region. **The requirement is Mach 2-4 and the term on offer is
Mach 1, by its own definition. It cannot absorb the requirement, whatever one
thinks of its physical justification.**

Three published analyses say the same thing in their own systems, in the
published text:

- **Palle et al. 2020, A&A 638, A61** (GJ 3470 b): "The inclusion of the
  broadening of the lines due to turbulence (vturb = 5kT/3m, where m is the
  mass of a helium atom), in addition to the standard Doppler broadening, was
  not enough to explain the measured broadening in the observations"
- **Lampon et al. 2021, A&A 647, A129** (HD 189733 b): "Including a turbulence
  broadening component (see Sect. 3.3), the profile broadens, but it is still
  narrower than the measured absorption. ... Hence, our hydrodynamic model
  alone is not able to explain the width of the absorption profile."
- **Czesla et al. 2022, A&A 657, A6** (HAT-P-32 b): "Thermal broadening alone
  remains insufficient to explain the observed line width, even if the
  turbulence term is included."

**Source-code fact, not a published statement.** In the copy of p-winds v1.4.7
in this workspace, `turbulence_broadening=False` is the default of both
transit entry points (`p-winds/p_winds/transit.py:126` and `:318`), and the
term is implemented only inside the `_method == 'average'` branch
(`transit.py:492`). Selecting `wind_broadening_method='formal'` therefore
drops it silently, with no warning and no error. The Cherubim et al. (2026)
released script sets `turbulence_broadening=True` explicitly and uses
`'average'`, so their run does carry the term; but a reader who reproduces a
p-winds fit in the `'formal'` mode gets a different line for a reason the
argument list does not advertise.

### 8.3 What a freely fitted width returns instead

The two analyses that let the width float, rather than tying it to the
temperature, land at two to three times the sonic term:

- **Salz et al. 2018, A&A 620, A97** (HD 189733 b, CARMENES). Free Doppler
  parameter `b`, defined in their Table 3 footnote as "(b) Doppler parameter
  given by sqrt(2) times the velocity dispersion of the Gaussians." Mid-transit
  values: **b = 15.1 +/- 0.6 km/s** for the extended model (f = 20 %) and
  **11.3 +/- 0.5 km/s** for the compact one. Worth recording explicitly:
  **this paper contains no turbulence term at all** -- a full-text search for
  "turbulen" over the published version returns nothing.
- **Czesla et al. 2022** (HAT-P-32 b). Free velocity dispersion `vt,a`:
  **15.2 +/- 3.1 km/s** (annulus), **11.2 +/- 3.1 km/s** (annulus plus
  rotation), **11.7 +/- 1.0 km/s** (up-orbit stream) at He 10833.

Our own requirement, sigma = 9.33 km/s, is at the low end of that set. The
gap between what a physical model produces and what a free fit wants is not
specific to this planet.

### 8.4 The same shortfall from an independent code

The closest published analogue is **Taylor et al. (2026, ApJ 999, 214)**, a
1-D self-consistent model with lower-atmosphere coupling, energy balance,
multispecies diffusion and non-LTE H(n = 2) -- the same class of model as
EXHALE -- applied to HD 189733 b. From the abstract:

> "HD189733b exhibits comparable He I depths, but the broadest reported
> profiles require ∼12 km s−1 of additional nonthermal broadening, whereas
> more recent measurements are narrower, consistent with our predictions."

and from Section 3.2:

> "Several prior analyses have matched the HD 189733b width by introducing ad
> hoc (micro)turbulent or Gaussian broadening terms; here, we show that this
> step remains necessary even in our self-consistent framework."

**They tested the same axes tested here, and all of them failed.** Section 4:

> "We tested variations in the free parameters (heating efficiency, Kzz,
> stellar activity level/XUV input, base composition, and lower-boundary
> placement): configurations that increase the width either deepen the core
> transit depth too much or violate the Halpha constraints, and none
> simultaneously reproduce both core depth and width."

with the size of the residual stated in the same paragraph:

> "this additional broadening remains insufficient, underestimating the
> observed line width by a factor of ∼2-3."

Their proposed mechanism, which they label speculative themselves:

> "Given that excited He is mostly produced by recombination of helium ions,
> this broadening can be driven by plasma turbulence. A particularly
> interesting idea is that HD189733b hosts a population of trapped ions on
> low-latitude magnetic field lines that form an energetic, turbulent plasma
> population around the planet. ... While these ideas are speculative, we
> encourage further exploration of plasma effects and magnetic field
> interactions"

The first clause is the same physics section 3 establishes independently
here: the metastable is a recombination product, so whatever moves it has to
move the recombining ions.

Rumenskikh et al. (2022), as reported by Taylor et al., reached the observed
width in 3-D only by prescribing equatorial jets "of the order of 10–20 km
s−1" together with an empirical reduction of H line cooling and He/H = 0.005.

### 8.5 How LHS 1140 b compares

HD 189733 b needs ~12 km/s of added dispersion on top of a model line that is
already thermally broad at ~10^4 K; LHS 1140 b needs 9.33 km/s on top of a
2.4 km/s model line, i.e. a *larger relative* correction (a factor 3.0 in
width against their 2-3) around a far quieter star. The mechanism Taylor et
al. favor -- activity-driven plasma turbulence on trapped field lines --
scales the wrong way for a 131 d rotator with L_X/L_bol = 6.5e-6.

**The conclusion this supports.** The shortfall is not an EXHALE defect and
not a peculiarity of this planet. It is a property of the model class: 1-D
steady-state photoevaporative outflows produce He 10830 lines two to three
times narrower than the ones that are measured, across hosts from an active K
dwarf to an inactive M4.5 dwarf and across two independently written codes.
The response of the field has been to add the missing width by hand -- as a
sonic turbulence term (Lampon, Salz 2016, Czesla, and hence the p-winds run in
Cherubim et al. 2026), as a free Gaussian (Taylor et al. 2026), or as
prescribed jets (Rumenskikh et al. 2022) -- and, where the added term is the
sonic one, to report in the published text that it is not enough (section 8.2).

**The conventional term does not close the gap here.** The Mach-1 convention
multiplies the Doppler parameter by sqrt(11/6) = 1.354, which is what
`EXHALE_TRANSIT_TURB=1` applies. On this model that takes the intrinsic sigma
from 2.42 to 3.28 km/s and the convolved FWHM from 0.261 to 0.322 A -- still
2.6x short of the observed 0.841 A. Switching turbulence on is therefore a
matter of comparing like with like against the paper's p-winds run, not a
route to the observed width.

`../LHS1140b/turbulence_broadening_literature.md` records the wider in-repo
survey of which published He 10830 analyses carry a turbulence term and which
argue against it.

The paper itself reads its pre-ingress signal as a leading tail and notes
that such tails have been "interpreted as resulting from stellar winds or
interactions between the magnetic fields of the star and planet" (Cherubim et
al. 2026, main text, citing Matsakos, Uribe & Konigl and McCann et al.).

### 8.6 One citation to handle with care

Stated as an observation, not as a correction to anyone's result. Lampon et
al. (2023) and Yan et al. (2024) both cite **Salz et al. (2018)** as a source
of the `v_turb = sqrt(5kT/3m)` convention. In the published Salz et al. (2018)
the string "turbulen" does not occur at all: that paper fits a free Doppler
parameter `b` (section 8.3), which is a different object -- unconstrained
rather than pinned to the local sound speed. The intended reference appears to
be to Salz et al. (2018) as the source of the *measured* width, and the
convention itself traces to Salz et al. (2016), Eq. (5), which does state it.
As written, the attribution invites the reader to believe that a Mach-1
turbulence term is what reproduced the HD 189733 b width, when in that paper
it was a free 15.1 km/s Doppler parameter.

---

## 9. What was not determined

* **The stellar wind of LHS 1140.** No measurement or published model
  mass-loss rate or wind speed for this star was located and verified for
  this memo. Section 2 is therefore parameterized in p_ram rather than
  anchored, and the statement "Mdot_* of tens of solar would be needed" is a
  requirement derived here, not a comparison against a published value. This
  is the single largest gap: an astrospheric or rotation-evolution estimate
  for a mid-M dwarf at P_rot = 131 d would convert the section 2 table into a
  verdict without the conditional.
* **Charge-exchange cross sections for the helium channels.** The energetics
  in section 3.3 are computed here from ionization potentials and the 2^3S
  excitation energy; the cross sections themselves (magnitude of
  `He+ + H -> He + H+` at 1-100 eV/amu, and whether the endothermic
  `-> He(2^3S)` branch has any appreciable cross section above its 46 km/s
  threshold) were not verified against published atomic data. The conclusion
  of section 3 does not rest on them -- section 3.2 disposes of the shocked
  layer regardless -- but a published number would close the argument
  cleanly.
* **Whether any published work computes helium ENAs at all**, as opposed to
  hydrogen ENAs (Holmstrom et al., Ekenback et al., Tremblin & Chiang and
  successors, all hydrogen). Not verified here.
* **A homogenized width comparison.** Section 8.1 now carries the survey, but
  the entries are quoted as published: different instruments, different line
  extractors, different treatments of stellar contamination, and no common
  deconvolution of the instrumental profile. The clustering at 19-33 km/s is
  therefore a statement about published numbers, not a uniformly reduced
  measurement. Two further limits: the sample is dominated by hot Jupiters and
  hot Neptunes, so there is no temperate super-Earth to compare LHS 1140 b
  against; and only Guilluy et al. (2024), Zhang et al. (2023), Salz et al.
  (2018) and Czesla et al. (2022) fit a width freely at all -- most He 10830
  papers report depth, not width.
* **Whether a displaced stellar residual blended with a narrow planetary line
  could mimic the observed profile.** Carried over unresolved from
  `lhs1140b_width_measurement_audit.md`; the released data do not include the
  stellar template.
* **The shape of the released window.** The spectrum spans -90.8 to
  +37.1 km/s about the red pair, so a component broader than about 40-60 km/s
  is degenerate with the continuum and is neither required nor excluded by
  these fits. A pickup-ion or ENA population would naturally have a
  dispersion of order the stellar wind speed, i.e. hundreds of km/s -- the
  wrong shape as well as the wrong magnitude, but not testable in this
  window.

---

## 10. What to compute next, and what each would decide

Ordered by how much they would move the verdict.

1. **Let the wind reach its critical point.** `collisional_validity.md`
   established that LHS 1140 b's wind never crosses its sonic point inside
   the 30 R_p domain (max Mach 0.54-0.57 for the three runs it measured); the
   two flux-closure runs used here touch Mach 1.00 in their outermost cell,
   which is the boundary, not a resolved critical point. Either way the mass
   flux is set at the outer boundary rather than at a critical point.
   Section 4.2 shows that even an ordered outflow scaled to the energy
   ceiling gives only FWHM ~8 km/s -- but that ceiling was computed with
   *this* density profile. Extending the domain past the
   critical point, with an outer boundary that does not pin the pressure,
   would test whether a genuinely transonic solution redistributes the
   metastable helium outward far enough to change both the width and the
   b-distribution of section 1.2. **This is the one calculation that could
   move the width without new physics**, and it is a well-defined run.

2. **Put radiation pressure on the 2^3S level into the transit calculation.**
   `a_rad` = 140 cm/s^2 against a planetary gravity of 74 cm/s^2 at 5 R_p is
   not a term one leaves out by inspection. The calculation is a drift
   velocity `a_rad x min(t_life, t_coll)` added to the metastable population's
   line-of-sight velocity, using quantities the solution already carries. It
   would settle whether the outer 25 % of the absorption is systematically
   blueshifted by several km/s, which is testable against the measured
   +2.0 km/s shift.

3. **Anchor the stellar wind.** Obtain a published Mdot_* and v_sw for a
   quiet mid-M dwarf at P_rot ~ 131 d, evaluate p_ram at 0.0946 AU, and read
   the standoff off the section 2 table. This converts section 2 from
   "incompatible for any wind strong enough to matter" to a specific
   statement about this star.

4. **Test the two-component hypothesis against a physical model rather than
   Gaussians.** Section 1.1 rules out a narrow core using free Gaussians. The
   sharper test is to take the EXHALE profile as the narrow component with its
   amplitude fixed by the model (not free) and ask what broad component the
   residual demands: with the model line 10x over the allowed narrow depth,
   the answer is expected to be "none exists", and confirming that closes the
   hybrid option explicitly.

5. **Check the He/H and base-density degeneracy against the b-distribution,
   not the depth.** All of section 2's conclusions used the model's
   `2 pi b N(b)` weighting. The four solutions tested agree
   (b50 = 4.4-5.3 R_p), but they share a lower boundary condition. A solution
   with a much more extended metastable layer -- if one exists within the
   observational constraints on depth and red/blue ratio -- would weaken the
   section 2 exclusion.

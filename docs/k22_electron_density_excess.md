# Why the electron density is 3-12 times Koskinen et al. (2022) Model A (2026-09-05)

**Judgement first.** The excess is entirely H+ (99.9 percent of n_e in every
cell). Both codes solve the same local balance for it, and in an H2-bearing
wind that balance is not photoionization against radiative recombination: the
proton is destroyed by charge exchange with hot H2 (their R10, our
`rk_R10_Hp_H2v4`, `1e-9 exp(-21900/T)`), and that channel is **70 times**
radiative recombination in their outer wind. Our n_e is high because R10 is
weak in our wind -- our H2 is 3.8 times lower and our gas hotter at 2 r_base --
and because our photoionization rate there is 2.6 times theirs. The single
input behind the latter is a geometric convention: **Model A divides the
stellar XUV flux by 4** ("uniform redistribution of energy around the planet",
Koskinen et al. 2022 section 3; defined in Koskinen et al. 2013a section 2 as
"we divided the incident stellar flux by a factor of 4"), while every hot-Uranus
run of ours used `Rate/2`, i.e. a factor 2. The H2 deficit and the temperature
excess are a feedback driven by that difference, not independent defects.

Everything below was measured on the Koskinen-composition run `k22_test`
(section 6 of the record; snapshot at 618,327 steps, outer ionization steady to
under 1 percent over the last 190,000 steps) and on 500-step restarts from it.

## 1. Which ion

| r/r_base | n_e ours | n_e Model A | ratio | H+ share of n_e (ours) |
|---|---|---|---|---|
| 1.119 | 4.1e7 | 5.2e6 | 7.8 | 99.98 % |
| 1.306 | 3.2e7 | 1.0e7 | 3.1 | 99.92 % |
| 2.015 | 1.0e7 | 1.7e6 | 5.9 | 99.45 % |
| 2.985 | 9.1e6 | 7.8e5 | 11.7 | 98 % |

The neutral H densities agree to a few percent (4.1e9 vs 3.9e9 at 1.119;
1.7e9 vs 1.7e9 at 1.194), so this is the ionization fraction, not the gas.

## 2. What the balance is, in both codes

In EXHALE the H/H+ partition is a local equilibrium recomputed every step (only
H2 is transported), so the proton row is exactly
`P n_H = k_R10 n_H2 n_H+ + alpha_B n_e n_H+`. In Model A the same terms plus a
transport term; their Figure 15 (digitized, `docs/p23_published_profiles.md`
section 7.3) gives at 2.7 R_p = 2.01 r_base: production 300, loss 190,
|transport| 130 cm^-3 s^-1. Their loss decomposes, with their own densities and
rates, into R10 = 1e-9 exp(-21900/4120) x 2.0e7 x 1.7e6 = 167 and radiative
7.5e-13 x (1.7e6)^2 = 2.2. **Their proton sink is charge exchange with H2, not
recombination.**

Ours at 2.015 r_base: production 660 (P = 2.2e-5 s^-1, calibrated from the
run's own heating column: heat_HI / n_HI reproduces the two-band power law at
`Rate/2` to 0-13 percent), R10 = 1.6e-11 x 5.3e6 x 1.0e7 = 850, radiative 46.
Same equation, same rate coefficients (all of Koskinen's R1-R23 are in
`mol_rates.f90`); different inputs.

## 3. The inputs that differ, and the ones that do not

| candidate | test | effect on n_e | verdict |
|---|---|---|---|
| secondary ionization (ours Immediate; Model A: "100 percent of this energy heats the atmosphere") | switched off, 500-step restart | 1.2-1.5x at 1.05-1.19; **1.00** beyond 1.3 | minor, deep layer only |
| spectral hardness (power law -1: 56 percent of energy above 100 eV; eps Eri SED: 24 percent) | index -2, 500-step restart | deep layer 7.8 -> 5.2x; **outer 6.1 -> 10.6x** | softer makes the outer excess WORSE: sigma ~ E^-3, so the same erg/cm^2/s ionizes 3.7x more |
| recombination coefficient | Koskinen R1 vs our Badnell case B, 800-5300 K | theirs 1.2-1.5x larger | 10-20 percent, in their favour |
| steadiness of the outer ionization | x(2.0), x(3.0), x(4.0) over 600k steps | 0.26->0.20, 0.77->0.74, 0.83->0.82; <1 percent over the last 190k | steady; not a transient |
| **dayside dilution 1/2 vs 1/4** | `Rate/4`, 500-step restart | **1.2-1.6x immediately**, H2 and T not yet moved | the convention Model A uses; grows with the feedback below |
| **n_H2 at 1.5-3 r_base** | measured | ours 3.8x lower at 2.0 | the largest single factor, and a feedback (section 4) |

## 4. The feedback

At 2.0 r_base our H2 is destroyed 7 times faster per molecule than theirs
(loss/n_H2 = 2.1e-4 vs 3.0e-5 s^-1), and the R10 part alone is **19 times**
faster: 3.3x from the temperature (5290 vs 4120 K in exp(-21900/T)) and 6x from
n_H+. Less H2 means less H3+ (ours 930 vs 7200 cm^-3 at 2.0 -- the infrared
coolant) and less R10 recombination, so the gas is hotter and more ionized,
which destroys more H2. The push on this loop is the photoionization input:
our P and photo-heating are ~2.6x theirs at 2 r_base, half of it the dilution
convention and the rest their heavier self-shielding by the denser neutral
envelope the loop itself preserves.

## 5. Two corrections made along the way

* The first reading of "uniform redistribution of energy around the planet"
  in the 2022 paper (line near equation 4) is about the **Bond albedo and
  effective temperature** of the lower atmosphere, not the XUV. That inference
  was withdrawn. The XUV statement is the one quoted above (2022 section 3,
  "we use the same mean solar XUV spectrum as Koskinen et al. (2013a, 2013b)
  and assume uniform redistribution of energy around the planet"), and 2013a
  defines it: "In general, we divided the incident stellar flux by a factor of
  4 to account for uniform redistribution of energy around the planet" and
  "In order to simulate a global average, we divided the flux by a factor of 4
  in the model."
* The `Mdot` agreement of the P23 gate (10.34-10.42 vs 10.28) was partly a
  coincidence of conventions: `Rate/2 + Mdot/2` halves the flux AND
  subtracts log 2 from the reported rate; Model A quarters the flux and reports
  the full sphere. The like-for-like configuration is `Rate/4 + Mdot`, which
  the code maps `Rate/4` to.

## 6. What follows

The matching configuration for the Koskinen comparison is `2D approximate
method: Rate/4 + Mdot` -- LW off, IR field off, carrier transport on, Koskinen
composition initial condition -- marched to steady state. The 500-step
restart says n_e drops 1.2-1.6x at once; how far the H2/T feedback then takes
it is what the long run measures. The run has not been started.

Runs: `scratchpad/p54ad/{ne_base,ne_nosec,L_base,L_nosec,L_pl2,L_nosec_pl2,L_rate4}`.

## 7. With the Model-A convention: what closed and what the remainder is (r4_cont, ~690k steps)

`Rate/4 + Mdot`, Koskinen composition, LW off, IR off, transport on; steady
to a few percent (T 7 percent, x_H2 0.6 percent over the last 450k steps):

| | ours | Model A |
|---|---|---|
| full-sphere Mdot | **2.1e10 g/s (log 10.32)**, flat 1.5-3.0 r_base | 1.9e10 (log 10.28) |
| T at 1.15 / 1.20 / 1.30 | 1752 / 1957 / 2285 K | 1720 / 1900 / 2270 K |
| f(H2) at 1.10 / 1.30 / 2.0 / 3.0 | 0.934 / 0.754 / 0.371 / 0.263 | 0.976 / 0.843 / 0.508 / 0.403 |
| n_e at 2.0 / 3.0 | 5.8e6 / 5.2e6 (3.2x / 6.7x) | 1.8e6 / 7.8e5 |
| T at 3.0 | 3245 K | 4740 K |

**Closed by the convention alone**: the mass-loss rate agrees like for like
(the earlier 10.34-vs-10.28 was `Rate/2` with a halved report), the layer's
temperature over 1.15-1.30 is within 60 K, and the outer n_e excess halved
(5.9 -> 3.2x at 2.0; 11.7 -> 6.7x at 3.0). At 2.0 r_base the photoionization
rates are now equal (9.3e-6 vs 8.3e-6 s^-1) and the R10 sink per proton is
within 1.6x.

**The remainder is a method difference, and it is measured.** In the outer
wind the ionization time exceeds the flow time:

| r/r_base | v [km/s] | 1/P [s] | r/v [s] | P r/v | x, local equilibrium (ours) | x, advected along our own flow | x, Model A |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.87 | 1.7e5 | 6.1e4 | 0.35 | 0.037 | 0.037 (start) | 0.008 |
| 2.0 | 3.41 | 1.1e5 | 2.1e4 | 0.19 | 0.129 | **0.057** | 0.020 |
| 2.4 | 5.41 | 9.8e4 | 1.6e4 | 0.16 | 0.327 | **0.078** | 0.055 |
| 3.0 | 7.64 | 9.5e4 | 1.4e4 | 0.15 | 0.585 | **0.106** | 0.080 |
| 4.0 | 9.71 | 9.3e4 | 1.4e4 | 0.16 | 0.703 | 0.142 | -- |

Local equilibrium overshoots Model A by 5.9x at 2.4 and 7.3x at 3.0; the
advected fraction lands within **1.4x and 1.3x** of it -- with our own rates,
velocities, temperatures and H2 densities, and starting from our own (already
4.6x too high) value at 1.5.

The gas leaves each shell three to seven times faster than it can be ionized,
so the ionization fraction is not the local equilibrium value: it is whatever
the parcel has accumulated on its way up. EXHALE recomputes the H/H+
partition from scratch every step as a local equilibrium (the design
statement is in `diffusive_photochemistry.f90`, section 1 of its header:
"f_sp is not advected by the hydro update, it is recomputed from scratch
each step by ioniz_eq"), so it sits on the equilibrium curve by construction.
Model A advects every species, and its Figure 15 shows the transport term
carrying 43 percent of the proton production at 2.0 r_base and 56 percent
at 2.4. Integrating dx/dt = P(1-x) - (k_R10 n_H2 + alpha n_e) x along OUR
flow, with OUR P, v, T and n_H2, from our x at 1.5 r_base, gives x = 0.057 at
2.0, 0.078 at 2.4 and 0.106 at 3.0 -- against Model A's 0.020, 0.055 and
0.080, and against our equilibrium values 0.129, 0.327 and 0.585. The
advected fraction reproduces Model A's outer ionization to within 30-40
percent; the local equilibrium misses it by a factor six to seven. The same over-ionization is why our outer wind is colder: gas that is
already ionized at 3 r_base cannot absorb the photons that are still heating
theirs (their proton production at 2.4 r_base is 170 cm^-3 s^-1 against our
74), and why H3+ is 100 times low there.

For the hot Jupiters the local-equilibrium treatment was built for, P r/v is
large and the equilibrium holds; for a light planet at 0.05 au with the flux
quartered, it does not above ~1.5 r_base. **Advecting the ionization state
(H, H+ as transported species, the same way H2 is now carried) is the method
change this comparison asks for.** Under way 2026-09-05 as the default-off key
`Ionization transport: True`.

**On the deep layer (1.05-1.10, our n_e 3.6-8x).** The digitized Model-A
n_e there is not usable: at 1.045 the table has e- = 5.1e6 against H+ = 6.9e4
and H3+ = 2.7e3, which violates the charge balance the paper itself states
(n_e = n_H+), and those rows were flagged uncertain when digitized. The solid
n_e comparisons are the ones at 2.4 and 3.0 r_base, where H+ and e- were both
read and agree.

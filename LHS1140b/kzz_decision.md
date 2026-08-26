# Choosing `He_Kzz` for LHS 1140 b — material for a user decision

**Status: settled. Sections 1-4 lay the two halves of the evidence side by
side, section 5 makes the recommendation, and the user adopted it on
2026-08-25 (section 0). Section 6 is what that value does to the composition
the He 10830 line implies.**

Prepared 2026-08-25, after Phase D of `../docs/binary_diffusion_design.md`
(binary H/He element diffusion, through decision D3). Companion documents:

- `kzz_literature.md` — the literature survey behind section 1, with the
  verification status of every number.
- `../docs/binary_diffusion_design.md` section 9.4 — one-paragraph summary
  and pointer.

---

## 0. Adopted 2026-08-25: `He_Kzz = 1.0e9` (user decision)

The recommendation of section 5 is the value the LHS 1140 b runs use:

```
He_diffusion: True
He_Kzz: 1.0e9
```

The consequence is not cosmetic. He/H = 0.55 was calibrated in the well-mixed
limit, and section 3 shows it cannot reproduce the measured line at any
`K_zz`. With the adopted value the composition the line implies is
**He/H = 2.09** instead, measured in section 6. Every composition statement
about this planet is therefore conditional on the eddy coefficient, the same
way it is already conditional on the assumed SED.

## 1. Why there is a decision to make

`He_Kzz` defaults to 0, i.e. molecular diffusion only. On LHS 1140 b that
default removes the helium the observation is about.

The wind base sits at **p = 1.0 microbar** (`Base BC: pressure 1.0`; the run
log reports the derived `n_0 = 3.2049e13 cm^-3`), at **T_0 = 226.0 K**, on a
planet with **g = 1837 cm/s^2** (computed here from M_p = 5.60 M_earth,
R_p = 1.730 R_earth). At that temperature and gravity helium settles faster
than this slow wind lifts it, so with `He_Kzz = 0` the homopause falls
*below* the base and the helium column is gone by 1.05 R_p. Measured in
`exhale/heh0p55_diff_ctrl`: `(He/H)/HeH` is 0.012 at 1.02 R_p and 5e-4 at
1.05 R_p, and the He 10830 red-pair equivalent width collapses to
2.5e-5 %A.

That matters because **the He/H = 0.55 of this planet's run was calibrated in
the well-mixed limit**. The diffusion-off case `exhale/heh0p55` gives a
red-pair EW of 1.109 %A against the measured 1.108 +/- 0.030 %A
(Cherubim et al. 2026, both put through `../he_line_metrics.py`); that
agreement is what fixed He/H = 0.55. Turning diffusion on with no eddy term
does not refine that result, it destroys it. The question is which eddy
coefficient the planet actually has, and what the model does with it.

### The structure that decides the answer

What competes with `He_Kzz` is the molecular coefficient `D_eff` that the
diffusion operator computes per cell (stage-resolved, memo section 2.6),
written by `EXHALE_DIFFUSION_CHECK=1` into `diffusion_faceflux.txt`. Measured
on the converged `K_zz = 1e10` state, with the pressure of each radius from
the control run:

| r [R_p] | p [microbar] | `D_eff` [cm^2/s] |
|---|---|---|
| 1.000 (base) | 1.0 | 4.6e5 |
| 1.010 | 0.21 | 1.1e7 |
| 1.020 | 0.080 | 4.4e7 |
| 1.050 | 0.017 | 7.7e8 |
| 1.100 | 5.9e-3 | 5.5e9 |
| 1.200 | 2.0e-3 | 2.5e10 |
| 1.500 | 3.5e-4 | 1.0e11 |
| 2.000 | 5.9e-5 | 1.4e11 |

The base value agrees with the Banks & Kockarts hard sphere evaluated by hand
at the base conditions, 8.0e5 cm^2/s (`kzz_literature.md` section 1); across
the whole scan the first-face value stays between 4.6e5 and 1.3e6.

Two readings follow, and together they are the difficulty:

1. **Any `K_zz` above ~1e6 already mixes the base cell.** The homopause
   leaves the base as soon as the eddy term exceeds ~5e5 cm^2/s.
2. **That is nowhere near enough.** `D_eff` rises by four decades between the
   base and 1.2 R_p, so a `K_zz` that mixes the base is irrelevant a
   hundredth of a radius higher.

Where the line is made matters too. In the well-mixed reference run the
metastable population peaks at 1.25 R_p but half of the `n(2^3S) r dr`
integral lies above 4.2 R_p and 90% below 10.5 R_p — the He 10830 line is
made from helium that has crossed the ionization front and been carried
outward. Whether any helium gets that far is settled near the front.

---

## 2. What the literature gives at 1 microbar

Full table, formulas, and verification status: `kzz_literature.md`. Values
extrapolated or evaluated at the base pressure of this run:

| source | object | `K_zz` at 1 microbar [cm^2/s] | inside the source's stated range? |
|---|---|---|---|
| Charnay, Meadows & Leconte (2015), `K_zz = K_zz0 P_bar^-0.4` | GJ 1214 b, warm sub-Neptune | 7.5e8 (pure H2O) to 7.5e9 (100x solar) | no — about 1.5 decades above the GCM top at 3e-5 bar |
| Parmentier, Showman & Lian (2013), `K_zz = 5e8 P_bar^-1/2` cm^2/s | HD 209458 b, hot Jupiter | 5e11 | yes — stated to ~1 microbar |
| Moses et al. (2011), constant plateau above the GCM top | hot Jupiters | ~1e11 to 4e11 | yes, by assumption |
| **Taylor et al. (2025)** | HD 209458 b, H/He escape model | **1e9** | yes — quoted at a model boundary of exactly 1e-6 bar |
| Arfaux & Lavvas (2023) homopause scaling, evaluated with **this planet's** g and T | **LHS 1140 b** | **2.8e7 (v/v_J = 1) to 1.4e8 (v/v_J = 5)** | yes — homopause plateau |
| Blain, Charnay & Bezard (2021) | K2-18 b, temperate sub-Neptune | 1e6 nominal, 1e5-1e10 explored | not pressure-resolved |

Verification: Charnay, Taylor, both Koskinen (2013) papers, Koskinen (2022)
and Cherubim were read from publisher PDFs; Parmentier, Moses,
Arfaux & Lavvas, Blain and Zhang & Showman were read from arXiv with key
sentences confirmed present in the ADS full-text index of the published
article, and are marked per row in `kzz_literature.md`. The Moses value that
Taylor et al. (2025) quote in m^2/s appears to be a unit slip; both readings
are recorded there.

Three things to carry from this table.

- **Taylor et al. (2025) is the closest analogue in configuration** — an H/He
  escape model with multispecies diffusion whose lower boundary is the same
  1e-6 bar — and it adopts 1e9 cm^2/s. That is also the value already cited
  in `parameters.f90` and used by `examples/14_diffusion`. But it is a hot
  Jupiter.
- **The only entry derived for LHS 1140 b itself gives 3e7 to 1.4e8**, one to
  two decades lower, and every temperature scaling in the survey points the
  same way: the hot-planet numbers likely overestimate `K_zz` here. That
  direction is argued from the sources' own scalings, not measured.
- **Cherubim et al. (2026) has no eddy term, no molecular diffusion and no
  homopause anywhere in the paper.** Their H:He is a single fitted constant
  of a homogeneous Parker wind ("The p-winds model assumes a spherical,
  homogeneous outflow"), so it is a well-mixed number with no altitude
  dependence. Our scan says when that assumption is safe: only where the row
  below has `(He/H)/HeH` near 1 over the line-forming region, which no row
  achieves.

---

## 3. The sensitivity scan

Configuration identical to `exhale/heh0p55_diff_ctrl` in every respect except
the `He_Kzz` line and the wind residual target: GJ 1132 SED, He/H = 0.55,
`He_diffusion: True`, `He_ambipolar` at its default on, `He_alphaT: 0.0`,
spherical domain to 30 R_p. Route: restart from a converged state,
direct-steady JFNK (`EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0`),
then one `Do only PP` pass and `EXHALE_transit.py`.

Scripts, all in `exhale/`: `run_kzz_case.sh` builds a case from the control
and solves it; `continue_kzz_case.sh` repeats an invocation until the
composition drift stops moving; `ladder_kzz.sh` restarts a case from the case
one decade below it; `kzz_scan_table.py` reads every number below out of the
run directories.

**How the rows were obtained, and what had to be relaxed.** Restarting a case
at 1e8 or above from the `K_zz = 0` control asks the wind solver to absorb
the whole helium column in one step; its line search stalls and the
composition outer loop stops on its solver-failure guard after a single pass.
Rebuilding one decade at a time fixes the seeding but not the residual: on
this wind JFNK settles at `||R||` between 2.7e-3 and 4.7e-3 rather than the
1.0e-3 the control reaches, and any target below that is read as a failure.
The four upper rows therefore state a looser `Resid tol` — 3.0e-3, except
1e9 which needs 5.0e-3 — and the achieved `||R||` is in the table so the
reader can judge each row. A run at 2.0e-2 was tried and rejected: there the
solver stops responding to the eddy term altogether (`||R||` reaches 1.8e-2
and `log10 Mdot` sticks at 7.68 for every case), so 3e-3 to 5e-3 is the
tightest target each rung can actually meet. Every row below has `info = 0`
and a composition drift under 1.3e-3. First attempts are kept in each case's
`attempt_seed_ctrl/`.

| case | `K_zz` [cm^2/s] | `Resid tol` | info | `\|\|R\|\|` | drift | homopause r [R_p] | `log10 Mdot` [g/s] | He 10830 red depth [%] | red EW [%A] |
|---|---|---|---|---|---|---|---|---|---|
| `heh0p55_diff_ctrl` | 0 | 1.0e-3 | 0 | 9.96e-4 | 1.7e-4 | below the base | 7.61 | 9.97e-5 | 2.51e-5 |
| `heh0p55_diff_kzz1e6` | 1e6 | 1.0e-3 | 0 | 9.58e-4 | 9.6e-4 | 1.0004 | 7.68 | 1.00e-2 | 2.27e-3 |
| `heh0p55_diff_kzz1e7` | 1e7 | 1.0e-3 | 0 | 9.96e-4 | 7.3e-4 | 1.0108 | 7.68 | 2.57e-2 | 5.78e-3 |
| `heh0p55_diff_kzz1e8` | 1e8 | 3.0e-3 | 0 | 2.70e-3 | 9.4e-4 | 1.0282 | 7.68 | 2.86e-1 | 6.54e-2 |
| `heh0p55_diff_kzz1e9` | 1e9 | 5.0e-3 | 0 | 4.71e-3 | 6.3e-4 | 1.0556 | 7.73 | 1.289 | 3.02e-1 |
| `heh0p55_diff_kzz1e10` | 1e10 | 3.0e-3 | 0 | 2.81e-3 | 9.9e-4 | 1.1286 | 7.84 | 2.514 | 6.05e-1 |
| `heh0p55_diff_kzz1e11` | 1e11 | 3.0e-3 | 0 | 3.00e-3 | 7.6e-4 | 1.4354 | 7.84 | 2.538 | 6.18e-1 |
| *reference*: `heh0p55`, diffusion off | — | 1.0e-3 | 0 | — | — | — | 7.76 | 4.285 | **1.109** |
| *measurement* (Cherubim et al. 2026) | — | — | — | — | — | — | — | 1.24 +0.22/-0.23 | **1.108 +/- 0.030** |

The homopause is the radius where the run's own `D_eff` equals its `K_zz`,
interpolated from `diffusion_faceflux.txt`.

Elemental ratio `(He/H)/HeH`, from the `_adv` profiles:

| `K_zz` [cm^2/s] | 1.02 | 1.05 | 1.10 | 1.20 | 2 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.012 | 0.0005 | 0.0001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 1e6 | 0.074 | 0.006 | 0.002 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| 1e7 | 0.418 | 0.083 | 0.027 | 0.011 | 0.002 | 0.002 | 0.002 | 0.002 |
| 1e8 | 0.848 | 0.485 | 0.295 | 0.184 | 0.039 | 0.030 | 0.030 | 0.030 |
| 1e9 | 0.978 | 0.845 | 0.674 | 0.535 | 0.206 | 0.160 | 0.160 | 0.160 |
| 1e10 | 0.997 | 0.974 | 0.901 | 0.785 | 0.394 | 0.318 | 0.317 | 0.317 |
| 1e11 | 1.000 | 0.996 | 0.981 | 0.926 | 0.530 | 0.432 | 0.430 | 0.430 |

### What the scan shows

- **The minimum `K_zz` that puts helium back at the reservoir ratio just
  above the base is 1e10.** At 1e10 the ratio is 0.997 at 1.02 R_p and 0.974
  at 1.05 R_p; at 1e9 it is 0.978 and 0.845; at 1e8 it is 0.848 and 0.485.
  Read against the `D_eff` table of section 1, this is exactly the crossing
  condition: `K_zz` has to exceed `D_eff` out to the radius one wants mixed.
- **No value of `K_zz` makes the outer wind well mixed.** Even at 1e11 the
  ratio settles at 0.43 above 5 R_p. `D_eff` reaches 1.4e11 cm^2/s at
  2 R_p, so the eddy term loses to molecular diffusion in the region where
  most of the He 10830 line is made, and the eddy coefficient would have to
  be raised past any published value to change that.
- **The line therefore saturates below the measurement.** The red EW rises
  steeply from 2.5e-5 %A to 0.605 %A between `K_zz` = 0 and 1e10 and then
  stops: 1e11 adds only 2%. The diffusion-off value at the same reservoir
  composition is 1.109 %A. **With the Phase-D operator active, He/H = 0.55
  cannot reproduce the measured line at any `K_zz`** — the best the scan
  reaches is 55% of it.
- **The wind itself is nearly indifferent.** `log10 Mdot` moves from 7.61 to
  7.84 across eleven decades of `K_zz`, +0.23 dex, and the change is not
  monotonic in the way the line is (7.68 for 1e6 through 1e8, then rising).
  The eddy term is a composition knob, not a mass-loss knob.

---

## 4. The two halves overlaid

**The literature range and the model's well-mixed range do not overlap.**

- The literature at 1 microbar spans 1e7-1e11 cm^2/s depending on which
  argument one accepts, and the one scaling evaluated for *this* planet gives
  3e7-1.4e8.
- The model needs `K_zz >= 1e10` before helium is at the reservoir ratio even
  at 1.05 R_p, and no value makes the outer wind well mixed.

So the planet-specific literature estimate (3e7-1.4e8) sits in the part of
the scan where helium is *strongly* separated: `(He/H)/HeH` around 0.5-0.8 at
1.05 R_p and 0.03-0.04 in the outer wind, giving a red EW near 0.07 %A, about
17 times below the measurement. Taylor's hot-Jupiter value of 1e9 gives
0.30 %A, still 3.7 times below. Only the extrapolated hot-Jupiter GCM values
(1e10-5e11) reach the plateau of 0.60-0.62 %A, and even that is 45% short.

Two readings of that gap are available and the scan cannot choose between
them:

1. **The reservoir composition is higher than 0.55.** He/H = 0.55 was fitted
   without diffusion. With diffusion on, matching 1.108 %A requires more
   helium at the base. Cherubim et al. retrieve H:He below 1e-3, i.e. He/H
   above 1000, so a much larger reservoir ratio is not excluded by their fit
   — the He/H = 0.55 in our runs came from our own EW crossing, not from
   theirs.
2. **Something in the outer-wind separation is too strong.** The ratio
   freezing at 0.43 above 5 R_p even at `K_zz = 1e11` is a property of the
   operator in a slow, cool, helium-rich wind, and section 9.3 of the design
   memo established the freeze-out mechanism but on a hot Jupiter, not here.

Reading 1 is now measured, in section 6: at the adopted `K_zz` = 1e9 the
equivalent width crosses the measurement at He/H = 2.09. That does not
dispose of reading 2 -- a reservoir ratio four times larger and an
outer-wind separation that is too strong would look the same in this one
observable -- and separating them is a piece of work this document does not
attempt. What it does establish is that **the `He_Kzz` choice changes the
He 10830 line by four orders of magnitude and the mass-loss rate by 0.23 dex,
so it cannot be left at the default and it cannot be inherited from a hot
Jupiter without saying so.**

---

## 5. Recommendation, and the decision that is not mine

### Recommended: `He_Kzz: 1.0e9`

At that value:

- **He 10830 red-pair EW 0.302 %A**, red depth 1.289%, from
  `exhale/heh0p55_diff_kzz1e9`;
- **`log10 Mdot` = 7.73** (against 7.61 at `K_zz = 0` and 7.76 with diffusion
  off);
- homopause at 1.056 R_p; `(He/H)/HeH` = 0.98 at 1.02 R_p, 0.85 at 1.05,
  0.67 at 1.10, 0.16 in the outer wind;
- solver: `info = 0`, `||R||` = 4.7e-3 at a stated `Resid tol` of 5.0e-3,
  composition drift 6.3e-4 over 20 outer passes.

The reasons, in order:

1. It is the only value in the range with a published justification for a
   model of the same kind at the same pressure: Taylor et al. (2025) adopt
   `K_zz = 1e5 m^2/s` at a lower boundary of 1e-6 bar for an H/He escape
   model with multispecies diffusion, verified against the published PDF.
2. It is what EXHALE already uses elsewhere — `examples/14_diffusion` and the
   HD 209458 b Phase-D baselines — so the LHS 1140 b runs stay comparable to
   the rest of the tree without a second convention.
3. It sits between the two literature clusters: an order of magnitude above
   the planet-specific Arfaux & Lavvas scaling (3e7-1.4e8) and one to two
   below the hot-Jupiter GCM extrapolations (1e10-5e11), so it does not
   commit to either.
4. It is on the steep part of the response, not the plateau, which keeps the
   value honest: a reader can see that the answer depends on it.

**The recommendation is weak, and the reason is worth stating.** Every
argument for 1e9 is an argument from hot-Jupiter practice. The one estimate
derived from LHS 1140 b's own gravity and temperature is one to two decades
lower, and all three temperature scalings in the survey say cooler and less
irradiated should mean smaller. If the planet-specific scaling is right, the
value should be 1e8 or below, and the He 10830 line the model predicts at
He/H = 0.55 is then 17 times weaker than measured.

### Alternatives, each defensible

| choice | case for it | what it gives |
|---|---|---|
| `1.0e8` | The Arfaux & Lavvas homopause scaling with this planet's `g` and `T` gives 2.8e7-1.4e8 — the only planet-specific number in the survey. | EW 0.065 %A, `log10 Mdot` 7.68, homopause 1.028 R_p |
| `1.0e10` | The lowest value at which the base region is genuinely well mixed (`(He/H)/HeH` = 0.97 at 1.05 R_p), and inside the Charnay GJ 1214 b extrapolation. | EW 0.605 %A, `log10 Mdot` 7.84, homopause 1.129 R_p |
| `0` (keep the default) | States plainly that no eddy mixing is being modelled. | EW 2.5e-5 %A — the helium column is gone, and the run cannot be compared to the measurement at all |
| scan rather than choose | The spread is two decades and none of it is measured for this object. | the table of section 3 is the deliverable |

If a single number has to go into the paper runs, the honest form is to quote
the recommended value **and** the 1e8-1e10 bracket with the line strengths
those give, since the observable moves by a factor of 9 across it.

### This is a user decision

`He_Kzz` is not a number EXHALE can derive: it stands for turbulence and wave
breaking that a 1-D wind model does not resolve, and no measurement of it
exists for this planet. The default of 0 is deliberate — the design memo
records it as "an eddy term is a property of the atmosphere being modelled,
so it is stated, not inherited". This document supplies the range the
literature supports at the base pressure of this run and what the model does
across that range. **The value to adopt for the LHS 1140 b runs is for the
user to choose** — and 1.0e9 is what the user chose on 2026-08-25
(section 0), with the composition that follows from it in section 6.

Phase E would replace the hand-set value: `read_base_inp` already accepts a
`Kzz_base` key that overrides `He_Kzz`, so a lower-atmosphere handoff can
carry `K_zz` at the matching pressure instead of it being set by hand here.

---

## 6. Composition with diffusion, at the adopted `He_Kzz = 1.0e9`

Section 3 left the model 3.7 times below the measurement at He/H = 0.55, and
section 4 named two readings of that gap. This section measures the first of
them: **what reservoir ratio reproduces the measured line once the diffusion
operator is active at the adopted eddy coefficient.**

Configuration: identical to `heh0p55_diff_kzz1e9` in every respect except the
`He/H number ratio` line -- GJ 1132 SED, `He_diffusion: True`,
`He_Kzz: 1.0e9`, `He_ambipolar` default on, `He_alphaT: 0.0`, spherical
domain to 30 R_p, `Resid tol: 5.0e-3` as in that case. Each case restarts
from the converged state of the case beside it (He/H = 1 from the
diffusion-off `heh1`, then each higher ratio from the one below), takes the
direct-steady JFNK route repeatedly until the composition drift stops moving,
then one `Do only PP` pass and `EXHALE_transit.py`. Script:
`exhale/run_heh_diff_case.sh`; table: `exhale/heh_diff_scan_table.py`, which
reads the same red-pair equivalent width as `make_memo_figures.py` (the
instrument-convolved curve normalized to its own continuum, on the vacuum
window 10832.60-10834.20 A).

**No case needed the `K_zz` ladder of section 3.** Restarting from a
diffusion-off solution of the *same* composition -- which is well mixed near
the base, where the adopted `K_zz` also mixes -- is a step the line search
takes directly: every case below reached `info = 0` in its first pass-set and
then stopped moving. The relaxations of section 3 that remain in force are
the looser `Resid tol` of 5.0e-3 and the achieved `||R||` quoted per row.

| case | He/H | info | `\|\|R\|\|` | drift | homopause r [R_p] | `log10 Mdot` [g/s] | red depth [%] | red EW [%A] |
|---|---|---|---|---|---|---|---|---|
| `heh0p55_diff_kzz1e9` | 0.55 | 0 | 4.71e-3 | 6.3e-4 | 1.0556 | 7.73 | 1.289 | 0.302 |
| `heh1_diff_kzz1e9` | 1.00 | 0 | 3.27e-3 | 6.4e-4 | 1.0507 | 7.80 | 2.406 | 0.579 |
| `heh2_diff_kzz1e9` | 2.00 | 0 | 3.17e-3 | 6.5e-4 | 1.0465 | 7.80 | 4.155 | **1.065** |
| `heh2p13_diff_kzz1e9` | 2.13 | 0 | 3.56e-3 | 7.7e-4 | 1.0461 | 7.80 | 4.396 | **1.129** |
| `heh4_diff_kzz1e9` | 4.00 | 0 | 1.85e-3 | 6.7e-4 | 1.0437 | 7.81 | 6.037 | 1.636 |
| *reference*: `heh0p55`, diffusion off | 0.55 | 0 | -- | -- | -- | 7.76 | 4.285 | 1.109 |
| *measurement* (Cherubim et al. 2026) | -- | -- | -- | -- | -- | -- | 1.24 +0.22/-0.23 | **1.108 +/- 0.030** |

Elemental ratio `(He/H)/HeH` from the `_adv` profiles:

| case | 1.05 | 2 | 5 | 10 | 20 |
|---|---|---|---|---|---|
| `heh0p55_diff_kzz1e9` | 0.845 | 0.206 | 0.160 | 0.160 | 0.160 |
| `heh1_diff_kzz1e9` | 0.836 | 0.233 | 0.191 | 0.191 | 0.191 |
| `heh2_diff_kzz1e9` | 0.834 | 0.255 | 0.208 | 0.207 | 0.207 |
| `heh2p13_diff_kzz1e9` | 0.835 | 0.268 | 0.221 | 0.220 | 0.220 |
| `heh4_diff_kzz1e9` | 0.832 | 0.281 | 0.235 | 0.234 | 0.234 |

### The crossing

The two runs at He/H = 2.00 and 2.13 straddle the measurement, so the
crossing is interpolated inside a 6.5% wide bracket rather than extrapolated:

- **central, EW = 1.108 %A: He/H = 2.09**;
- 1 sigma low, 1.078 %A: He/H = 2.03 (inside the same bracket);
- 1 sigma high, 1.138 %A: He/H = 2.16 (bracket 2.13-4.00).

so **He/H = 2.09, with 2.03-2.16 for the measurement's 1 sigma**, from a
log-log interpolation of the scanned points solved by Brent's method. The
quoted 1 sigma range is the propagated measurement error alone; it carries
none of the model uncertainty, which is dominated by `K_zz` itself (section
3: the line moves by a factor of 9 across 1e8-1e10) and by the SED.

The escape rate barely moves across the scan -- `log10 Mdot` = 7.73 to 7.81,
the same near-independence of composition that the diffusion-off scan showed
-- so the composition is fixed by the line and not by the wind.

### He 10830 profile at the crossing

`heh2p13_diff_kzz1e9`, all from `tpm_He10830_metrics.txt` (three-Gaussian fit
of the instrument-convolved curve, air frame):

| quantity | `heh2p13_diff_kzz1e9` | `heh0p55`, diffusion off | measured |
|---|---|---|---|
| red depth [%] | 4.396 | 4.285 | 1.254 |
| blue depth [%] | 0.643 | 0.624 | -- |
| red/blue | 6.84 | 6.86 | -- |
| FWHM [A] | 0.242 | 0.244 | 0.841 |
| shift [A] | -7.1e-4 | -5.2e-4 | -- |
| sigma [A] | 0.0916 | 0.0925 | -- |
| red EW [%A] | 1.129 | 1.109 | 1.108 +/- 0.030 |

The line the model makes at the crossing is the same line the diffusion-off
solution made, to a few percent in every metric. The width discrepancy
against the measurement -- FWHM 0.24 A against 0.84 A -- is unchanged and
unrelated to diffusion; it is the broadening question of
`../docs/lhs1140b_exhale_vs_pwinds.tex`.

### What absorbs, and why the reservoir has to be four times larger

Comparing `heh2p13_diff_kzz1e9` with the diffusion-off `heh0p55` at the same
observable, from the `_adv` profiles:

| r [R_p] | 1.05 | 1.25 | 2 | 5 | 10 |
|---|---|---|---|---|---|
| n(He) elemental, He/H = 0.55, no diffusion [cm^-3] | 3.09e10 | 6.22e8 | 3.00e7 | 6.70e5 | 6.24e4 |
| n(He) elemental, He/H = 2.13, `K_zz` = 1e9 [cm^-3] | 3.61e10 | 7.78e8 | 3.36e7 | 6.98e5 | 6.35e4 |
| n(2^3S), He/H = 0.55, no diffusion [cm^-3] | 1.94 | 21.8 | 19.8 | 4.47 | 0.562 |
| n(2^3S), He/H = 2.13, `K_zz` = 1e9 [cm^-3] | 7.76 | 34.2 | 16.6 | 4.41 | 0.572 |

The two solutions carry **the same helium in the outer wind** -- 6.98e5
against 6.70e5 cm^-3 at 5 R_p, 4% apart, and the metastable densities there
agree to 1.3%. That is the whole mechanism. Separation removes 78% of the
helium above 2 R_p (`(He/H)/HeH` = 0.22 there against 1.00 with diffusion
off), and raising the reservoir by 3.9x puts back what separation took, so
the absorbing column above the ionization front is restored to what it was.
The metastable integral `int n(2^3S) r^2 dr` confirms it: 1164 against 1129,
3% apart, matching the 2% difference in equivalent width.

Where the two differ is below 2 R_p, and it is a redistribution rather than a
gain: the metastable peak moves inward, from 1.25 to 1.17 R_p, and is 57%
higher there (34.2 against 21.8 cm^-3), while above 2 R_p the diffusion case
is slightly the weaker of the two (16.6 against 19.8 at 2 R_p). The
half-weight radius of the `n(2^3S) r dr` integral moves from 4.17 to
3.99 R_p and the 90% radius is unchanged at 10.5 R_p, so the line is still
made mostly out beyond the front in both cases -- diffusion has not moved the
line-forming region, it has changed how much helium reaches it.

### Implication

**The composition inferred for LHS 1140 b is conditional on diffusion and on
`K_zz`, exactly as it is already conditional on the assumed SED.** The three
numbers now on record for the same measurement, same planet, same line:

| assumption | He/H implied |
|---|---|
| GJ 1132 SED, no diffusion (well mixed) | 0.55 |
| GJ 699 SED, no diffusion (well mixed) | 0.060 |
| GJ 1132 SED, diffusion on at `K_zz` = 1e9 | **2.09** |

The SED changes the answer by a factor of 9 downward, the eddy coefficient by
a factor of 3.8 upward, and neither is measured for this system. A retrieval
that quotes a single H:He from this line -- as Cherubim et al. (2026) do,
from a homogeneous Parker wind with no diffusion term at all -- is quoting a
number whose value depends on assumptions their model does not carry.

### 6.1 The crossing at the two alternative eddy coefficients

The subsection above fixes the composition at the adopted `K_zz` = 1e9 only.
The two alternatives of section 5, 1e8 and 1e10, have now been rescanned the
same way, so the crossing carries a `K_zz` bracket instead of a caveat.

Method identical to the 1e9 scan in every respect -- same control input, same
`Resid tol: 5.0e-3`, same direct-steady JFNK route, same `Do only PP` pass and
`EXHALE_transit.py`, same red-pair equivalent width on the vacuum window
10832.60-10834.20 A. The only change is the `He_Kzz` line and the composition
bracket the crossing needed. `exhale/heh_diff_scan_table.py` now reads both
`He/H number ratio` and `He_Kzz` out of each case's `input.inp` rather than
assuming them, and takes a scan name (`1e8`, `1e9`, `1e10`) as its argument;
re-running the 1e9 group through it reproduces the numbers of section 6 to the
last digit, which is the check that the three scans are measured alike.

Runs, all `info = 0` with composition drift below 1.0e-3 in their first
pass-set, all seeded as noted:

| case | seed | He/H | `\|\|R\|\|` | drift | homopause r [R_p] | `log10 Mdot` [g/s] | red depth [%] | red EW [%A] |
|---|---|---|---|---|---|---|---|---|
| `heh0p55_diff_kzz1e8` | (section 3) | 0.55 | 2.70e-3 | 9.4e-4 | 1.0282 | 7.68 | 0.286 | 0.065 |
| `heh2p7_diff_kzz1e8` | `heh3_diff_kzz1e8` | 2.70 | 3.48e-3 | 9.8e-4 | 1.0213 | 7.80 | 4.119 | **1.050** |
| `heh3_diff_kzz1e8` | `heh4_diff_kzz1e8` | 3.00 | 4.55e-3 | 6.5e-4 | 1.0211 | 7.80 | 4.458 | **1.149** |
| `heh3p5_diff_kzz1e8` | `heh4_diff_kzz1e8` | 3.50 | 4.03e-3 | 9.5e-4 | 1.0207 | 7.80 | 4.914 | 1.285 |
| `heh4_diff_kzz1e8` | `heh4_diff_kzz1e9` | 4.00 | 3.35e-3 | 8.0e-4 | 1.0205 | 7.80 | 5.290 | 1.400 |
| `heh10_diff_kzz1e8` | `heh10` (diffusion off) | 10.0 | 4.79e-3 | 5.4e-4 | 1.0006 | 7.76 | 83.76 | 23.45 (not usable, below) |
| `heh0p55_diff_kzz1e10` | (section 3) | 0.55 | 2.81e-3 | 9.9e-4 | 1.1286 | 7.84 | 2.514 | 0.605 |
| `heh1_diff_kzz1e10` | `heh1` (diffusion off) | 1.00 | 1.75e-3 | 8.5e-4 | 1.1056 | 7.80 | 3.235 | 0.801 |
| `heh1p4_diff_kzz1e10` | `heh1p5_diff_kzz1e10` | 1.40 | 3.09e-4 | 9.8e-4 | 1.0997 | 7.81 | 4.157 | **1.058** |
| `heh1p5_diff_kzz1e10` | `heh1_diff_kzz1e9` | 1.50 | 3.24e-3 | 6.2e-4 | 1.0988 | 7.80 | 4.230 | **1.087** |
| `heh1p6_diff_kzz1e10` | `heh1p5_diff_kzz1e10` | 1.60 | 2.79e-3 | 7.5e-4 | 1.0978 | 7.80 | 4.455 | **1.146** |

As at 1e9, no case needed the `K_zz` ladder of section 3: seeding from a
converged solution of the same composition -- diffusion-off where one exists,
otherwise the neighbouring composition already solved at the same `K_zz` --
was a step the line search took directly. `log10 Mdot` is 7.80 +/- 0.01 in
every run that carries the line, so the composition is again fixed by the
line and not by the wind.

**The crossing exists below saturation at 1e8.** The concern was that a
helium-rich reservoir would saturate the line before it reached the
measurement. It does not: at 1e8 the equivalent width is still rising with a
local logarithmic slope near 0.73 across He/H = 2.7 to 4.0, and it passes the
measured 1.108 %A at He/H = 2.87 -- inside a 11% wide bracket, not
extrapolated.

Crossings, log-log interpolation of the scanned points solved by Brent's
method, against the measured 1.108 +/- 0.030 %A:

| `K_zz` [cm^2/s] | He/H at the crossing | 1 sigma range | bracketing runs | achieved `\|\|R\|\|` of those runs |
|---|---|---|---|---|
| 1e8 | **2.87** | 2.78 - 2.97 | 2.70 [EW 1.050] and 3.00 [EW 1.149] | 3.48e-3, 4.55e-3 |
| 1e9 | **2.09** | 2.03 - 2.16 | 2.00 [EW 1.065] and 2.13 [EW 1.129] | 3.17e-3, 3.56e-3 |
| 1e10 | **1.54** | 1.47 - 1.59 | 1.50 [EW 1.087] and 1.60 [EW 1.146] | 3.24e-3, 2.79e-3 |

The 1 sigma range is the propagated measurement error alone within one
`K_zz`; at 1e9 and 1e10 the 1 sigma-low end falls in the adjacent bracket
(2.00-2.13 and 1.40-1.50 respectively), which is measured too.

**The degeneracy is shallow and very nearly a power law.** Over the two
decades,

    d log(He/H) / d log K_zz  =  -0.136

(-0.139 across 1e8-1e9, -0.133 across 1e9-1e10), that is
He/H proportional to `K_zz`^(-0.14). A factor of 10 in the eddy coefficient
moves the inferred reservoir ratio by only 37%, and the full 1e8-1e10 range
of section 5 -- the whole span the literature leaves open at 1 microbar --
spans He/H = 1.5 to 2.9, a factor of 1.9.

That is a much weaker degeneracy than the line depth itself suggested:
section 3 measured a factor of **9** in equivalent width across the same two
decades at fixed composition. The two are consistent because the line
responds to composition much faster than it responds to `K_zz` -- the local
slope d log EW / d log(He/H) is 0.73 to 0.92 near the crossing against
d log EW / d log K_zz of about 0.5 at He/H = 0.55 -- so a large change in the
eddy coefficient is bought back by a small change in the reservoir. The
practical consequence is that **`K_zz` is not the dominant uncertainty on the
inferred composition**, even though it is the dominant uncertainty on the
predicted line at fixed composition. The SED remains the larger lever: it
moved the diffusion-off crossing by a factor of 9 (0.55 to 0.060), against
1.9 for the entire `K_zz` range.

Updating the table of section 6:

| assumption | He/H implied |
|---|---|
| GJ 1132 SED, no diffusion (well mixed) | 0.55 |
| GJ 699 SED, no diffusion (well mixed) | 0.060 |
| GJ 1132 SED, diffusion on at `K_zz` <= 1e5 | 4.68 |
| GJ 1132 SED, diffusion on at `K_zz` = 1e8 | 2.87 |
| GJ 1132 SED, diffusion on at `K_zz` = 1e9 (adopted) | **2.09** |
| GJ 1132 SED, diffusion on at `K_zz` = 1e10 | 1.54 |
| GJ 1132 SED, diffusion on at `K_zz` = 1e11 | 1.11 |

#### Every decade from `K_zz` = 0 to 1e11

The three decades above are the ones section 5 left open. The scan has since
been carried down to `K_zz` = 0 and up to 1e11, one decade at a time, by the
identical route -- same control input, same `Resid tol: 5.0e-3`, same
direct-steady JFNK, same `Do only PP` pass and `EXHALE_transit.py`, same
red-pair window 10832.60-10834.20 A -- and the numbers of the three original
rows are reproduced to the last digit by the same reader, which is the check
that the whole ladder is measured alike.

Two things had to be measured separately, because a crossing scan at every
decade would have been wasted work below the point where the eddy term stops
acting at all.

**Where the eddy term stops acting.** The operator adds `K_zz` to the
molecular coefficient, so the decade at which the eddy term disappears is
fixed by `D_eff` at the base, which these He-rich solutions carry at
1.24-1.32e6 cm^2/s (first face of `diffusion_faceflux.txt`, r = 1.0003 R_p).
Below that the eddy term is a small correction to a coefficient the wind
already has. Measured at a fixed reservoir He/H = 5, against the same run
with `He_Kzz` absent:

| `K_zz` [cm^2/s] | red EW [%A] | EW - EW(`K_zz` = 0) | in units of the 0.030 %A measurement error |
|---|---|---|---|
| 0 | 1.173690 | — | — |
| 1e1 | 1.173690 | +6e-7 | 0.00002 |
| 1e2 | 1.173696 | +6e-6 | 0.0002 |
| 1e3 | 1.173749 | +5.9e-5 | 0.002 |
| 1e4 | 1.173849 | +1.6e-4 | 0.005 |
| 1e5 | 1.178971 | +5.3e-3 | 0.18 |
| 1e6 | 1.226357 | +5.3e-2 | 1.76 |

**`K_zz` <= 1e4 is not distinguishable from `K_zz` = 0 by this line**, and
1e5 is only marginally so; the first decade that moves the equivalent width
by more than the measurement error is 1e6. The ordering is the ratio
`K_zz`/`D_eff` at the base: 1e6 raises the base coefficient by 77%, 1e5 by
8%, 1e4 by 0.8%. The crossing composition was therefore bracketed at
`K_zz` = 0 and at every decade from 1e5 upward, and the decades 1e1 to 1e4
are reported as what they are, indistinguishable from `K_zz` = 0.

**No decade runs into saturation before the crossing.** The concern of the
1e8 scan applies with more force at low `K_zz`, where the crossing sits at
He/H near 5: at 1e6 the equivalent width is still climbing at He/H = 10
(2.2465 %A, local logarithmic slope 0.86 through the crossing), so every
crossing below is interpolated inside a bracket of two solved runs, none
extrapolated and none against a saturated line.

Crossings, log-log interpolation of the scanned points solved by Brent's
method against the measured 1.108 +/- 0.030 %A:

| `K_zz` [cm^2/s] | He/H at the crossing | 1 sigma range | bracketing runs [EW in %A] |
|---|---|---|---|
| 0 (molecular diffusion only) | **4.68** | 4.54 - 4.83 | 4.6 [1.0923] and 4.8 [1.1341] |
| 1e1 | indistinguishable from `K_zz` = 0 | — | EW probe above, 0.00002 sigma |
| 1e2 | indistinguishable from `K_zz` = 0 | — | EW probe above, 0.0002 sigma |
| 1e3 | indistinguishable from `K_zz` = 0 | — | EW probe above, 0.002 sigma |
| 1e4 | indistinguishable from `K_zz` = 0 | — | EW probe above, 0.005 sigma |
| 1e5 | **4.68** | 4.55 - 4.83 | 4.6 [1.0917] and 4.8 [1.1325] |
| 1e6 | **4.45** | 4.33 - 4.59 | 4.4 [1.0965] and 4.7 [1.1606] |
| 1e7 | **3.78** | 3.66 - 3.96 | 3.6 [1.0627] and 3.8 [1.1127] |
| 1e8 | **2.87** | 2.78 - 2.97 | 2.70 [1.0505] and 3.00 [1.1492] |
| 1e9 (adopted) | **2.09** | 2.03 - 2.16 | 2.00 [1.0653] and 2.13 [1.1289] |
| 1e10 | **1.54** | 1.47 - 1.59 | 1.50 [1.0872] and 1.60 [1.1460] |
| 1e11 | **1.11** | 1.07 - 1.15 | 1.00 [1.0192] and 1.20 [1.1801] |
| *reference*: diffusion off, well mixed | 0.55 | — | section 6, figure `lhs1140b_ew_vs_heh.pdf` |

The measured crossing at 1e5 (4.6797) sits 0.1% from the `K_zz` = 0 one
(4.6748), which is the same statement as the equivalent-width probe, made
the other way.

**The power law of the section above holds only above 1e7.** The local
exponent `d log(He/H) / d log K_zz`, decade by decade:

| decades | exponent |
|---|---|
| 1e5 - 1e6 | -0.022 |
| 1e6 - 1e7 | -0.071 |
| 1e7 - 1e8 | -0.119 |
| 1e8 - 1e9 | -0.139 |
| 1e9 - 1e10 | -0.133 |
| 1e10 - 1e11 | -0.141 |

so -0.136 is the exponent of the 1e8-1e10 window it was fitted in, and it
carries up to 1e11 (-0.138 as a mean over 1e8-1e11) but not down: below 1e7
the curve rolls over onto the `K_zz` = 0 plateau at He/H = 4.68. The whole
landscape is therefore three regimes, not one power law -- a plateau at
He/H = 4.68 for `K_zz` <= 1e5, a roll-off through 1e6-1e7, and a
`K_zz`^(-0.14) decline from 1e8 up -- and none of them reaches the
well-mixed crossing of 0.55, which the operator approaches only
asymptotically: even at 1e11 the inferred reservoir is twice it.

The three regimes are drawn in `../docs/figures/lhs1140b_heh_vs_kzz.pdf`
(built by `make_memo_figures.py`, which re-solves every crossing from the run
directories through its own `red_ew`, against the equivalent width it reads
straight from the Cherubim et al. digitization, 1.1084 +/- 0.0295 %A; the
crossings agree with the table above to better than 0.1%).

Runs behind the added rows, all `info = 0` with composition drift below
1.0e-3, all with `log10 Mdot` = 7.80-7.82:

| case | seed | He/H | `\|\|R\|\|` | drift | homopause r [R_p] | red depth [%] | red EW [%A] |
|---|---|---|---|---|---|---|---|
| `heh4p5_diff_kzz0` | `heh4p6_diff_kzz0` | 4.5 | 2.375e-3 | 7.1e-4 | below the base | 4.179 | 1.0683 |
| `heh4p6_diff_kzz0` | `heh5_diff_kzz1e6` | 4.6 | 3.159e-3 | 5.4e-4 | below the base | 4.263 | 1.0923 |
| `heh4p8_diff_kzz0` | `heh5_diff_kzz1e6` | 4.8 | 3.110e-3 | 6.6e-4 | below the base | 4.406 | 1.1341 |
| `heh4p9_diff_kzz0` | `heh4p8_diff_kzz0` | 4.9 | 2.240e-3 | 8.4e-4 | below the base | 4.445 | 1.1462 |
| `heh5_diff_kzz0` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 6.5e-4 | below the base | 4.540 | 1.1737 |
| `heh5_diff_kzz1e1` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 6.5e-4 | below the base | 4.540 | 1.1737 |
| `heh5_diff_kzz1e2` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 6.5e-4 | below the base | 4.540 | 1.1737 |
| `heh5_diff_kzz1e3` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 6.5e-4 | below the base | 4.540 | 1.1737 |
| `heh5_diff_kzz1e4` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 6.5e-4 | below the base | 4.540 | 1.1738 |
| `heh4p5_diff_kzz1e5` | `heh4p5_diff_kzz0` | 4.5 | 1.762e-3 | 9.5e-4 | below the base | 4.177 | 1.0681 |
| `heh4p6_diff_kzz1e5` | `heh4p6_diff_kzz0` | 4.6 | 2.352e-3 | 5.6e-4 | below the base | 4.259 | 1.0917 |
| `heh4p8_diff_kzz1e5` | `heh4p8_diff_kzz0` | 4.8 | 2.234e-3 | 9.6e-4 | below the base | 4.399 | 1.1325 |
| `heh4p9_diff_kzz1e5` | `heh4p9_diff_kzz0` | 4.9 | 1.498e-3 | 6.8e-4 | below the base | 4.462 | 1.1514 |
| `heh5_diff_kzz1e5` | `heh5_diff_kzz1e6` | 5.0 | 3.053e-3 | 5.8e-4 | below the base | 4.558 | 1.1790 |
| `heh2_diff_kzz1e6` | `heh2_diff_kzz1e9` | 2.0 | 4.485e-3 | 9.2e-4 | below the base | 1.794 | 0.4217 |
| `heh4p3_diff_kzz1e6` | `heh4p4_diff_kzz1e6` | 4.3 | 2.363e-3 | 7.3e-4 | below the base | 4.191 | 1.0717 |
| `heh4p4_diff_kzz1e6` | `heh5_diff_kzz1e6` | 4.4 | 3.159e-3 | 5.5e-4 | below the base | 4.277 | 1.0965 |
| `heh4p7_diff_kzz1e6` | `heh5_diff_kzz1e6` | 4.7 | 3.138e-3 | 9.2e-4 | below the base | 4.495 | 1.1606 |
| `heh5_diff_kzz1e6` | `heh4_diff_kzz1e8` | 5.0 | 3.797e-3 | 7.8e-4 | below the base | 4.718 | 1.2264 |
| `heh10_diff_kzz1e6` | `heh10` (diffusion off) | 10.0 | 3.981e-3 | 8.7e-4 | below the base | 8.337 | 2.2465 |
| `heh3p4_diff_kzz1e7` | `heh4p3_diff_kzz1e6` | 3.4 | 1.868e-3 | 8.4e-4 | 1.0065 | 3.938 | 0.9991 |
| `heh3p6_diff_kzz1e7` | `heh4_diff_kzz1e8` | 3.6 | 3.969e-3 | 7.1e-4 | 1.0064 | 4.162 | 1.0627 |
| `heh3p8_diff_kzz1e7` | `heh4_diff_kzz1e8` | 3.8 | 3.985e-3 | 5.7e-4 | 1.0063 | 4.334 | 1.1127 |
| `heh3p9_diff_kzz1e7` | `heh3p8_diff_kzz1e7` | 3.9 | 3.325e-3 | 6.0e-4 | 1.0063 | 4.375 | 1.1251 |
| `heh4_diff_kzz1e7` | `heh3p9_diff_kzz1e7` | 4.0 | 2.225e-3 | 7.0e-4 | 1.0063 | 4.443 | 1.1458 |
| `heh4p2_diff_kzz1e7` | `heh4_diff_kzz1e8` | 4.2 | 4.015e-3 | 6.4e-4 | 1.0062 | 4.648 | 1.2055 |
| `heh1_diff_kzz1e11` | `heh1_diff_kzz1e10` | 1.0 | 3.056e-3 | 7.7e-4 | 1.3232 | 4.011 | 1.0192 |
| `heh1p2_diff_kzz1e11` | `heh1p4_diff_kzz1e11` | 1.2 | 4.399e-3 | 6.8e-4 | 1.3096 | 4.570 | 1.1801 |
| `heh1p4_diff_kzz1e11` | `heh1p4_diff_kzz1e10` | 1.4 | 8.623e-4 | 5.9e-4 | 1.2986 | 4.923 | 1.2844 |
| `heh1p5_diff_kzz1e11` | `heh1p5_diff_kzz1e10` | 1.5 | 4.454e-3 | 6.4e-4 | 1.2932 | 5.053 | 1.3258 |

The `He/H` = 5 probe rows are one wind solved at five values of `He_Kzz`;
their identical `\|\|R\|\|` is the eddy term being negligible, not a copied
file -- the equivalent widths differ in the sixth decimal place.

`heh_diff_scan_table.py` now takes `0`, `1e5`, `1e6`, `1e7` and `1e11` as
scan names alongside the three it already had, plus `flat` for the fixed
composition probe, whose crossing solve is skipped because it holds a single
He/H.

### 6.2 A limitation the He-rich end exposed, in the post-processing

`heh10_diff_kzz1e8` (He/H = 10) is in the table above but excluded from the
crossing fit, and the reason is a limitation of `post_process_adv.f90` rather
than of that run's wind solution.

The advection correction is gated cell by cell, and one of its three gates is
the **hydrogen** Damkohler number: where `Da = (dr/v)*(P_HI + alpha_HII*n_e)`
exceeds 100 the cell keeps its equilibrium ionization, on the stated grounds
that hydrogen relaxes to local equilibrium many times over while the gas
crosses the cell, so the equilibrium solution already is the solution of the
ODE. That reasoning is sound for hydrogen. But the correction the gate
controls carries the **whole** species vector, including the He 2^3S
population, whose relaxation time is far longer than hydrogen's; the gate
therefore leaves the metastable at its equilibrium value in cells where it is
not in equilibrium.

At He/H = 10 the equilibrium solution develops a sharp hydrogen ionization
front near r = 5.5 R_p -- n(HI) falls from 1.6e3 to 5.2 cm^-3 over half a
planetary radius, and x_HII saturates at 1 beyond it. Every cell above the
front trips the Damkohler gate: 116 of 503 cells keep the equilibrium
ionization against 13 to 32 in every other case of this section. The result
is a step in the `_adv` metastable profile at the gate radius, n(2^3S) going
from 18 cm^-3 just below the gate radius near 6.5 R_p to 206 cm^-3
just above it, an order of magnitude across the switch. The transit
synthesis reads the `_adv` profiles, so the run reports a 84% red-pair depth
and EW = 23.45 %A -- a number produced by the discontinuity, not by the wind.

**None of the crossing runs is affected.** In all ten other cases of this
section the uncorrected cells sit at r <= 1.008 R_p, in the base where gate
(i)/(iii) is meant to act, and not one cell above r = 1.2 R_p keeps its
equilibrium value. The crossings quoted in section 6.1 are therefore clean.

The limitation itself is pre-existing, is not specific to diffusion, and was
not fixed in this pass: gating a coupled species vector on a single species'
Damkohler number is wrong wherever the species' relaxation times differ, and
correcting it means either gating the correction species by species or
raising the gate to the slowest relaxation time in the system.

**Fixed, `../docs/Update_EXHALE.md` section 72** (2026-08-25): the gate now
uses the slowest relaxation rate of the solved species vector. Re-post-
processing `heh10_diff_kzz1e8` with the corrected binary leaves 15 of 503
cells at equilibrium, all at `r <= 1.0034 R_p`, removes the step in the
metastable profile, and gives a red depth of 21.5% and EW = 5.69 %A in place
of 83.8% and 23.45 %A. The other ten cases of this section are byte-identical
under the fix -- their pinned cells are held by gate (i)/(iii), which is
unchanged -- so every number of section 6.1 stands. The stored files of
`heh10_diff_kzz1e8` were left as they are and still carry the old gate; the
run is still outside the crossing fit, which is now a statement about the run,
not about the post-processing.

The crossings of section 6.1 have since been folded into
`../docs/lhs1140b_exhale_vs_pwinds.tex`: section 8.3 carries the full-decade
crossing table, the degeneracy figure
(`docs/figures/lhs1140b_heh_vs_kzz.pdf`) and a structure figure
(`docs/figures/lhs1140b_kzz_profiles.pdf`, built by `make_memo_figures.py`
from the four fixed-composition runs of section 3), and the composition item
of "Status against the goal" states the `K_zz` conditionality alongside the
SED one. The well-mixed crossing of 0.55 is still the memo's headline
number, now labelled as the well-mixed limit.

# Choosing `He_Kzz` for LHS 1140 b — material for a user decision

**2026-09-13: every `exhale/...` path in this file now lives under
`archive_20260830/exhale/...`. The solutions it describes predate the code
changes listed in `MODELS.md` section 5 and are superseded by `models/`;
read `MODELS.md` for the tree of record.**

**Status: settled. Sections 1-4 lay the two halves of the evidence side by
side, section 5 makes the recommendation, and the user adopted it on
2026-08-25 (section 0). Section 6 is what that value does to the composition
the He 10830 line implies, section 7 makes that composition an output of the
flux closure, and section 8 asks the two questions together.**

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
**He/H = 2.09** instead, measured in section 6 -- on a scalar H/He base. On
the photochemical profile the flux closure hands over, which carries C, N
and O and their cooling, the same line asks for **He/H = 11.73** (section 8). Every composition statement
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
written by `EXHALE_DIFFUSION_CHECK=1` into
`output/element_flux_profile.txt` (it was `./diffusion_faceflux.txt` in the
run root when this section was written). Measured
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
interpolated from `output/element_flux_profile.txt`.

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

> **Every He/H number in sections 6 and 6.1 is on the retired Taylor (2025)
> Penning coefficient and on the pre-section-96 secondary-ionization split,
> and none of it is current.** The whole ladder has been re-solved twice
> since: on the Garcia Munoz (2025) coefficient
> (`exhale/ladder_gm25/results.txt`) and then on the composition-renormalized
> SvS85 split of `docs/Update_EXHALE_stage1.md` section 96
> (`exhale/crossings_j96/results.txt`). The current values are
> He/H = 3.8075 on the `K_zz` <= 1e4 plateau, 1.6108 at the adopted
> `K_zz` = 1e9, 0.8368 at 1e11, 0.4132 well mixed on the GJ 1132 proxy and
> 0.0483 on the GJ 699 one; `docs/lhs1140b_exhale_vs_pwinds.pdf`
> Table 2 carries the full ladder. What did not change is the shape: three
> regimes, and the local exponent d log(He/H)/d log `K_zz` reproduced decade
> by decade to 0.003. The text below is kept as the record of the scan that
> was run.

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
1.24-1.32e6 cm^2/s (first face of `output/element_flux_profile.txt`,
r = 1.0003 R_p).
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

**Fixed, `../docs/Update_EXHALE_stage1.md` section 72** (2026-08-25): the gate now
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

## 7. The flux-closed composition (2026-08-27): the composition is now an output

Everything above treats `He/H` as an input to be chosen — the value at which
the modelled He 10830 line meets the measured one, conditional on `K_zz`.
Milestone E4 of `../docs/phase_e_flux_closure_design.md` removes that framing
on the lower-atmosphere side: **the matching-level composition is now solved
for, not stated.** The photochemical column is given the elemental escape
fluxes as an upper boundary condition, the wind is solved on the profile that
comes back, its own elemental fluxes are measured, and the two are iterated
to agreement (`../docs/Update_EXHALE_stage1.md` section 79;
`exhale/flux_closure/`).

The measured answer on this planet, from three starting fluxes spanning a
factor of ten:

| start | k at convergence | `He/H` at the match | `F_H` [g/s] | `F_He` [g/s] | `log10 Mdot` |
|---|---|---|---|---|---|
| `1.0 x` | 0 | 2.0923516 | 1.81433e7 | 1.32040e7 | 7.500 |
| `0.3 x` | 5 | 2.0923486 | 1.82031e7 | 1.34266e7 | 7.500 |
| `3.0 x` | 6 | 2.0923602 | 1.82119e7 | 1.34710e7 | 7.500 |

The three agree in `He/H` to 5.5e-6, so the closure is single-valued here.

**What this does and does not change about section 0.** It does not change the
adopted value or the number the line implies: the closure returns
`He/H = 2.0923` against the 2.09 the reservoir was started from, a change of
1.1e-3 and well inside the width of the crossing. What it changes is the
standing of that number. `He/H = 2.09` is no longer only the value that fits
the line; it is also the value the coupled lower atmosphere and wind settle
on, and the two agree.

**And it sharpens the conditionality section 0 already states.** The closure
says explicitly *why* the composition at the match is what it is, and the
reason is `K_zz` and not the escape. Of the 1.1e-3 change, only 2e-4 comes
from the escape flux itself: at `K_zz = 1e9 cm^2/s` eddy mixing homogenizes
the column all the way to 1 microbar, so a 1.8e7 g/s hydrogen escape imposed
at the photochemical model top moves `He/H` there by 1.8e-2 and moves it at
the match by 2e-4. The eddy coefficient sets the composition at the match;
the escape does not. Section 0's sentence that every composition statement
about this planet is conditional on the eddy coefficient is therefore not
weakened by the composition becoming an output — it is the mechanism behind
it, and the closure measured it.

The fractionation itself is real and lives above the match, where it always
did: the wind removes helium and hydrogen in the mass ratio 0.73 against the
base reservoir's 8.31, and `He/H` falls from 2.09 at the base to 0.183 at
30 R_p. What the closure adds is that this separation does not reach back
down through the eddy-mixed column to the matching level at the adopted
`K_zz`. At a much smaller `K_zz` it would, and that is the experiment section
6.1's decade scan would now be repeated as.

## 8. Composition under flux closure (2026-08-27): the reservoir the line needs

Section 7 closed the loop at one reservoir, `He/H = 2.09`, and found the
closure single-valued there. It did not ask what that solution's He 10830
line looks like. It does not reproduce the measurement: the red-pair
equivalent width of the flux-closed solution is **0.421 %A against the
measured 1.108 +/- 0.030**, a factor 2.63 low. The scalar-base run at the
same composition (`exhale/heh2p13_diff_kzz1e9`, section 6) gives 1.129 and
does reproduce it, so the deficit belongs to the profile the closure hands
over and not to the composition.

Where it comes from is the elemental C, N and O the Photochem profile
carries into the wind, and the cooling they add. At the same reservoir,
measured from the `_adv` profiles:

| | scalar base, `He/H` 2.13 | flux-closed profile, `He/H` 2.09 |
|---|---|---|
| max wind `T` [K] | 5330 | 3945 |
| `T` at 2 `R_p` [K] | 4140 | 2971 |
| max `n(2^3S)` [cm^-3] | 42.9 | 50.3 |
| `int n(2^3S) dr` [cm^-3 R_p] | 67.3 | 40.7 |
| `log10 Mdot` [g/s] | 7.80 | 7.40 |
| red EW [%A] | 1.129 | 0.421 |

The metastable peak is slightly *higher* in the closed solution but the
column above it is 1.7 times smaller, and the line is 2.6 times weaker: the
cooler wind keeps the metastable helium concentrated near the base instead
of carrying it out through the line-forming region.

### 8.1 The reservoir ladder under closure

So the question section 6 answered for the scalar base is reopened for the
closed system: **what reservoir `He/H` makes the flux-closed solution
reproduce the measured line?** The ladder below runs the closure driver
(`../src/utils/element_flux_closure.py`) at eight further reservoirs, each
started from the converged flux and converged wind of the point below it, so
every rung begins near its own fixed point. Everything else is the section-7
configuration: Photochem `clima` + chemistry, `K_zz = 1e9`, `p_match = 1e-6`
bar, the GJ 1132 proxy SED, `tol = 0.05`, `k_max = 8`. Runs:
`exhale/flux_closure/heh{3,5,8,9p7,10p3,10p7,11p1,12_ctl}/`, table
`exhale/flux_closure/closure_heh_table.py`. A ninth rung,
`exhale/flux_closure/heh12/`, reaches the same 12.01 reservoir in one 1.50x
jump from the 8.0 rung instead of one step from 11.11; it is a convergence
artifact and is recorded separately in section 8.3.

| reservoir `He/H` | k | `F_H` [g/s] | `F_He` [g/s] | `He/H` at match | `log10 Mdot` | red depth [%] | blue [%] | FWHM [A] | red EW [%A] | EW / measured | cold-trap `O/H` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.09 (sec. 7) | 0 | 1.814e7 | 1.320e7 | 2.0924 | 7.500 | 1.699 | 0.246 | 0.2326 | 0.4206 | 0.380 | 4.95e-7 |
| 3.0 | 3 | 1.369e7 | 1.491e7 | 3.0034 | 7.460 | 2.122 | 0.307 | 0.2386 | 0.5388 | 0.486 | 5.54e-7 |
| 5.0 | 4 | 8.799e6 | 1.692e7 | 5.0057 | 7.410 | 2.795 | 0.404 | 0.2497 | 0.7433 | 0.671 | 6.43e-7 |
| 8.0 | 4 | 5.722e6 | 1.894e7 | 8.0090 | 7.390 | 3.388 | 0.487 | 0.2608 | 0.9437 | 0.852 | 7.94e-7 |
| 9.7 | 3 | 4.768e6 | 2.043e7 | 9.7109 | 7.400 | 3.660 | 0.526 | 0.2664 | 1.0403 | 0.939 | 8.78e-7 |
| 10.3 | 1 | 4.498e6 | 2.069e7 | 10.3116 | 7.400 | 3.696 | 0.530 | 0.2674 | 1.0563 | 0.953 | 9.09e-7 |
| 10.7 | 1 | 4.338e6 | 2.095e7 | 10.7121 | 7.400 | 3.733 | 0.535 | 0.2684 | 1.0706 | 0.966 | 9.39e-7 |
| 11.1 | 1 | 4.192e6 | 2.124e7 | 11.1125 | 7.400 | 3.771 | 0.540 | 0.2694 | 1.0853 | 0.980 | 9.74e-7 |
| 12.0 | 2 | 3.903e6 | 2.200e7 | 12.0135 | 7.410 | 3.856 | 0.552 | 0.2714 | 1.1181 | 1.009 | 1.066e-6 |
| *12.0, discarded (sec. 8.3)* | 4 | 3.910e6 | 2.238e7 | 12.0135 | 7.420 | 4.630 | 0.663 | 0.2654 | 1.3126 | 1.185 | 1.066e-6 |
| *measurement* | | | | | | 1.254 | | 0.841 | **1.108 +/- 0.030** | 1 | |

Every rung returns a matching-level `He/H` equal to its reservoir to four
digits, the same result section 7 reported at 2.09 and for the same reason:
at `K_zz = 1e9` the eddy-mixed column reaches the microbar match, so the
escape flux does not set the composition there.

**The crossing.** The runs at 11.11 (EW 1.0853) and 12.01 (1.1181) straddle
the measured 1.108 +/- 0.030 %A directly, and both of them lie inside its
1 sigma band, so the answer is read inside a measured bracket:

- **central, EW = 1.108 %A: `He/H` = 11.73**;
- 1 sigma low, 1.078 %A: 10.92; 1 sigma high, 1.138 %A: 12.58

Two readings of the same ladder agree. Chord interpolation in log-log
between the two straddling rungs gives 11.731 (10.918-12.581); a quadratic
in log-log through the three upper rungs 10.71, 11.11 and 12.01 gives 11.735
(10.913-12.565). The centers are 0.005 apart. (While the discarded 12.01 of
section 8.3 was the top rung, the same two methods disagreed by 0.16.) The
quoted range carries the measurement error alone and none of the model
uncertainty, and only its upper edge continues the top chord past the
highest rung. At the crossing the escape rate is `log10 Mdot` = 7.40 to
7.41, i.e. 2.5-2.6e7 g/s, against the 2.03 +0.58/-0.67 e8 g/s of the
Cherubim et al. p-winds retrieval -- a factor 8 below it -- and the
elemental fluxes are `F_H` ~ 4.0e6 and `F_He` ~ 2.2e7 g/s. The escape rate
is nearly flat across the whole ladder (7.39 to 7.50 over a factor 5.7 in
reservoir), so here too the composition is fixed by the line and not by the
wind.

**So there is a composition that satisfies both conditions, and it is 5.6
times the section-6 value.** The flux-closure condition by itself does not
select a reservoir on this planet -- it is satisfied at every rung -- so the
line is what picks 11.73, and the closure's contribution is that the answer
is now consistent with a photochemical lower atmosphere carrying its own
C/N/O rather than with a bare H/He base. The price is a helium fraction of
92 per cent by number.

**How firmly the crossing is bracketed.** The center is still an
interpolation -- no rung was run at 11.73 -- but the bracket around it is
measured and the local slope no longer jumps across it. The equivalent width
grows sublinearly with the reservoir along the whole ladder, at
`d ln EW / d ln (He/H)` = 0.688 (2.09-3.00), 0.630, 0.508, 0.506, 0.254,
0.353, 0.372, 0.382 over the successive intervals up to 12.01. No interval
departs from its neighbours, so the flattening near 10 is curvature in the
relation and not one contaminated rung, and the upper end where the crossing
is solved has settled to 0.37-0.38. The sublinearity also settles that the
equivalent width does *not* saturate before the measurement is reached,
which is what the He-rich end of the scalar ladder (section 6) had made a
live possibility.

**Convergence quality, rung by rung.** The wind residual the JFNK finish can
reach degrades as the reservoir grows: the line search stalls on the base
contact mode and returns its best iterate rather than converging. The target
was therefore raised rung by rung, and the achieved values are recorded rather
than smoothed over:

| rung | wind `Resid tol` | achieved `\|\|R\|\|` | `info` | steady-window `F_H` spread |
|---|---|---|---|---|
| 2.09 | 1.0e-4 | 9.81e-5 | 0 | 0.44 % |
| 3.0 | 2.0e-4 | 1.40e-4 | 0 | 0.64 % |
| 5.0 | 2.0e-4 | 1.70e-4 | 0 | 0.80 % |
| 8.0 | 4.0e-4 | 3.35e-4 | 0 | 1.14 % |
| 9.7 | 4.0e-4 | 2.44e-4 | 0 | 2.14 % |
| 10.3 | 4.0e-4 | 2.56e-4 | 0 | 2.61 % |
| 10.7 | 4.0e-4 | 2.94e-4 | 0 | 3.03 % |
| 11.1 | 4.0e-4 | 3.22e-4 | 0 | 3.23 % |
| 12.0 | 4.0e-4 | 3.86e-4 | 0 | 3.30 % |
| *12.0, discarded* | 4.0e-4 | 2.53e-4 | 0 | 1.74 % |

The acceptance test the closure actually applies is the last column -- the
radial spread of the elemental flux over the window it is measured on -- and
it stays between 0.4 and 3.6 per cent (`F_He`, 3.58 % at the top rung)
against the 5 per cent tolerance at every rung. Section 8.3 is the reason
that test, and the residual norm beside it, are not sufficient on their
own. `Resid tol` became a configuration key of the closure driver in
the course of this ladder; the default is unchanged at 1.0e-4.

**What the crossing does not fix.** The line is still too narrow: FWHM 0.269
A at the crossing against the measured 0.841 A, unchanged from section 6 and
from the diffusion-off solutions, so the modelled red depth (3.8 per cent)
overshoots the measured 1.254 per cent by the same factor the width falls
short. The composition is inferred from the equivalent width for that
reason. And the wind these numbers come from is *unvalidated as a continuum
solution*: it is subsonic to the 30 `R_p` domain edge and the critical region
is transitional (`../docs/collisional_validity.md`).

### 8.2 The XUV grid against the 2025 non-detection

Cherubim et al. detect the line in 2024 and not in 2025, with a limit of
0.6 per cent in depth. The grid below asks what fraction of the fiducial XUV
reproduces that. The spectrum is the fiducial GJ 1132 proxy with its flux
column multiplied by the stated factor (`sed/lhs1140_sed_gj1132_at_b_xuv*.txt`;
`read_sed` recomputes `L_X` and `L_EUV` from the file, so the luminosity
lines of `input.inp` are provenance only). The lower atmosphere is held
fixed in every case -- the closure loop is not re-run, and neither is the
chemistry that made the profile -- so the grid isolates the wind's response.

It is run from both models that reproduce the 2024 equivalent width, because
neither is the other: the scalar-base solution of section 6 at `He/H` = 2.13
(EW 1.129), and the flux-closed solution of section 8.1 at `He/H` = 9.71
(EW 1.040). Runs `exhale/xuv*_heh2p13/` and `exhale/xuv*_closure9p7/`.

| `F_XUV` / fiducial | scalar `He/H` 2.13: red [%] | EW [%A] | `log10 Mdot` | flux-closed `He/H` 9.71: red [%] | EW [%A] | `log10 Mdot` |
|---|---|---|---|---|---|---|
| 1.00 | 4.396 | 1.129 | 7.80 | 3.660 | 1.040 | 7.40 |
| 0.33 | 0.778 | 0.180 | 7.38 | **0.602** | 0.154 | 6.74 |
| 0.30 | **0.591** | 0.136 | 7.35 | -- | -- | -- |
| 0.25 | 0.334 | 0.077 | 7.32 | 0.511 | 0.129 | 6.55 |
| 0.20 | 0.034 | 0.009 | 7.23 | 0.361 | 0.090 | 6.40 |
| 0.15 | 0.007 | 0.002 | 7.14 | 0.256 | 0.063 | 6.22 |
| 0.10 | 0.003 | 0.001 | 6.93 | 0.165 | 0.040 | 5.97 |
| 0.01 | 0.000 | 0.000 | 6.10 | 0.190 | 0.055 | 5.51 |

**The 0.6 per cent limit is met at about 0.3 of the fiducial XUV in both**,
which is the useful part of the result: the two models disagree by a factor
4.6 in reservoir and a factor 2.5 in escape rate, and still put the
non-detection at the same place. In the scalar model the depth crosses 0.6
between 0.30 (0.591 per cent) and 0.33 (0.778); in the flux-closed model it
sits at 0.602 per cent at 0.33 and 0.511 at 0.25, so the crossing is just
under 0.33. Taking the fiducial `F_XUV` = 33 erg/cm^2/s of the paper, the
2025 epoch would need about 10 erg/cm^2/s or less.

Two caveats belong with that number. The decline is far steeper in the
scalar model than in the closed one -- a factor 200 in depth between 0.30 and
0.20 against a factor 1.6 -- so the *shape* of the XUV dependence is a
property of the model and not a result; only the crossing agrees. And the
lowest closed-model point (0.01) is not monotonic with the 0.10 one
(0.190 against 0.165 per cent) at a wind that has dropped to
`log10 Mdot` = 5.5; that solution is not to be read as a limit.

### 8.3 The outer region a converged wind inherits (2026-08-27)

The 12.01 rung was first reached in one step from the 8.0 rung, a 1.50x jump
in reservoir (`exhale/flux_closure/heh12/`, converged at k = 4). The rung
kept in the table above reaches the same reservoir from 11.11, a 1.08x step
(`exhale/flux_closure/heh12_ctl/`, k = 2). The two agree on everything the
closure and the chemistry set, and disagree on the line:

| | 12.0 discarded (`heh12`) | 12.0 kept (`heh12_ctl`) |
|---|---|---|
| `He/H` at match | 12.013522 | 12.013520 |
| cold-trap `O/H` | 1.0656e-6 | 1.0660e-6 |
| `exhale_info` | 0 | 0 |
| achieved `\|\|R\|\|` | 2.529e-4 | 3.864e-4 |
| `log10 Mdot` | 7.420 | 7.410 |
| red depth [%] | 4.630 | 3.856 |
| red EW [%A] | 1.3126 | 1.1181 |

Same base, same matching-level composition to seven digits, both converged,
and a 17 per cent difference in equivalent width. The difference is in the
outer wind alone. Two measures that do not depend on where a spectral window
is placed separate them, and both are monotonic along the ladder except at
the discarded point (`r_drop` = the first radius outside 1.5 `R_p` where `T`
falls below half of `T(12 R_p)`; the last column is the share of the
He 2^3S radial column outside 10 `R_p`):

| reservoir | `r_drop` [`R_p`] | `T(12 R_p)` [K] | He 2^3S column outside 10 `R_p` |
|---|---|---|---|
| 2.09 | 30.00 | 893 | 0.0119 |
| 3.0 | 27.63 | 1035 | 0.0116 |
| 5.0 | 23.76 | 1310 | 0.0105 |
| 8.0 | 20.80 | 1632 | 0.0104 |
| 9.7 | 19.79 | 1789 | 0.0108 |
| 10.3 | 19.79 | 1834 | 0.0108 |
| 10.7 | 19.47 | 1863 | 0.0109 |
| 11.1 | 19.47 | 1891 | 0.0109 |
| 12.0 kept | 19.47 | 1951 | 0.0112 |
| **12.0 discarded** | **25.83** | 1957 | **0.0219** |

`r_drop` moves inward monotonically with the reservoir, 30.0 to 19.5 `R_p`,
and the outer share of the metastable column stays in 0.0104-0.0119 across
the ladder. The discarded solution departs from both -- 25.8 `R_p` where its
own reservoir gives 19.5, and twice the outer share of any other rung --
while its `T(12 R_p)` matches its control to 0.3 per cent. It is a hotter,
more extended outer atmosphere on an identical base.

*Read as the mechanism*: the outer state is inherited from the seed and the
JFNK finish does not re-solve it. The density out there contributes almost
nothing to the residual norm, so `||R||` does not separate the two states --
both report `info = 0`, and the one that is wrong reports the *smaller*
residual. He 10830 is optically thick at these columns, so a hot, extended
outer region fills the line wings and lifts the equivalent width even though
it holds about one per cent of the metastable column.

**The physical judgement.** The 20-30 `R_p` region in which the two
solutions differ is where the collisional-validity measurement puts the
exobase, 19.8-26.7 `R_p`, and where `Kn > 1` (`../docs/collisional_validity.md`).
*The fluid solution is not valid in exactly the region that makes the
difference.* The ground for discarding the original 12.01 solution is
therefore continuity of the trend, not a physical criterion; and by the same
token the absolute equivalent widths of the whole ladder carry a systematic
error that the closure tolerance -- a 5 per cent flux spread -- does not
catch. Its size is small, because the outer region is about one per cent of
the metastable column, but it is not zero.

**What was and was not tested.** Path independence was checked directly at
`He/H` = 3.0, from a seed 3.7x away in composition and in the downward
direction the ladder never takes: it returns the same `r_drop` to the digit,
with everything else inside the closure tolerance. The 5.0, 8.0 and 2.09
rungs were not checked this way. Separately, `exhale/flux_closure/heh10p3_ctl/`
re-runs the 10.31 rung from the same seed and reproduces
`Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` bitwise, which fixes
determinism but says nothing about seed dependence.

**The practical rule.** A ladder scan is climbed one step at a time; a large
jump in reservoir can leave an outer region the solver will not revisit and
no convergence test will flag. Record: `../docs/Update_EXHALE_stage1.md` section 83.

---

## 9. Is the pressure form of `K_zz` carried in the answer? (2026-08-28)

Sections 3 and 6.1 measure what the **magnitude** of `K_zz` is worth: a factor
9 in equivalent width between 1e8 and 1e10, and a `He/H` crossing that moves as
`K_zz^-0.14` wherever the eddy term acts at all. Everything above, and every
run of the profile branch (sections 7-8), states the eddy coefficient as a
**constant** — `--kzz-const 1.0e9` at all 102 levels of the handoff. That form
is a second assumption sitting on top of the value, and it had never been
tested. It is now.

### The test

`K_zz(p) = 1.0e9 (p / 1.1e-8 bar)^(-1/2)` — the saturated gravity-wave form,
`--kzz-power 0.5 --kzz-ref 1.0e9 --kzz-ref-bar 1.1e-8` — **anchored at the top
of the profile**, not at 1 bar. The anchor is the point: EXHALE gives every
cell above the shallowest tabulated level that level's value
(`eddy_diffusion_on_grid`), so anchoring at the top leaves the number the wind
inherits above the file unchanged and confines the change to the interval the
file describes. Implied values: **1.05e8 at the 1 microbar match**, 1.05e5 at
1 bar, 2.56e4 at the 16.73 bar deep boundary — a factor 9.5 under the constant
at the match and four decades under it at the bottom.

Two tests, because the profile hands over both the composition and the eddy
coefficient:

- **A, chemistry only.** The column is solved on `K_zz(p)`, then its `Kzz`
  column alone is written back to 1e9 before the handoff. The wind sees what
  it saw before; only the composition changes.
- **B, chemistry and wind.** The profile as produced.

The reference is the constant-`K_zz` profile re-solved with the same seed and
the same binary. Reservoirs: the `11.11` rung (`flux_closure/heh11p1`, k = 1)
and the `2.09` one on its high-trial-flux branch (`flux_closure/hi`, k = 6).
Runs and the measurement script: `exhale/kzz_power/`.

### Result

| run | red [%] | FWHM [A] | EW [%A] | EW/obs | He/H at 2 R_p | `log10 Mdot` |
|---|---|---|---|---|---|---|
| **reservoir 11.11** | | | | | | |
| recorded, other seed | 3.5998 | 0.2825 | 1.0853 | 0.980 | 1.5536 | 7.400 |
| reference, `K_zz` = 1e9 | 3.6524 | 0.2825 | 1.1019 | 0.994 | 1.5764 | 7.410 |
| A, chemistry only | 3.6630 | 0.2825 | 1.1050 | 0.997 | 1.5797 | 7.410 |
| B, chemistry and wind | 3.6591 | 0.2825 | 1.1031 | 0.996 | 1.5630 | 7.410 |
| **reservoir 2.09** | | | | | | |
| recorded, other seed | 1.6233 | 0.2467 | 0.4272 | 0.386 | 0.2301 | 7.500 |
| reference, `K_zz` = 1e9 | 1.6282 | 0.2467 | 0.4285 | 0.387 | 0.2307 | 7.500 |
| A, chemistry only | 1.6380 | 0.2467 | 0.4311 | 0.389 | 0.2320 | 7.500 |
| B, chemistry and wind | 1.6265 | 0.2467 | 0.4279 | 0.386 | 0.2299 | 7.500 |

At the 1 microbar match, the ratio `K_zz(p)` / constant:

| quantity | 11.11 rung | 2.09 rung |
|---|---|---|
| **consumed by the wind** | | |
| `X_He` | 0.9985 | 0.9985 |
| `X_C` | 0.9899 | 0.9903 |
| `X_N` | 0.9878 | 0.9875 |
| `X_O` | 0.9894 | 0.9821 |
| `T` | 1.0000 | 1.0000 |
| **carried, not consumed** | | |
| `q_CH4` | 0.443 | 0.665 |
| `q_CO` | 54.7 | 47.6 |
| `q_CO2` | 9.3e3 | 4.3e3 |
| `q_NH3` | 0.0038 | 0.0496 |
| `q_N2` | 28.7 | 28.1 |
| `q_HCN` | 0.566 | 2.95 |
| `q_H2O` | 0.714 | 0.804 |
| `q_H` | 4.99 | 5.55 |

**What did not move.** The climate is identical — deep boundary 411.4 K,
tropopause 2.636 bar at 185.5 K, cold-trap `f_H2O` = 5.590e-8, match
temperature 185.478 K, the same digits — and identical **structurally**: the
adapter's `--climate` path solves the radiative-convective column with no eddy
coefficient entering it and hands the photochemistry a (P, T) pair with the
`K_zz` slot empty. FWHM agrees to four decimals in all six runs; `log10 Mdot`
does not move at all; the elemental ratios agree to 1.2% or better.

**What did move: the carriers, by up to four decades.** Nitrogen changes
carrier outright — on the 11.11 rung `q_NH3` at the match falls 6.32e-6 ->
2.43e-8 while `q_N2` rises 1.21e-7 -> 3.46e-6, so an NH3-dominated column
becomes N2-dominated. Carbon partly follows (CH4 down 2.3x, CO up 55x, CO2 up
nearly four decades from a negligible base). Deeper it is worse: at 1 bar the
ratios are 1.3e5 for CO, 3.0e7 for CO2 and 8.8e3 for N2. Elemental N at the
match nevertheless moves by -1.2%, because a carrier swap conserves nuclei and
nuclei are what the handoff transmits.

**The crossing, against the path spread.** Converted with the local slope
`d ln EW / d ln(He/H)` = 0.382 over the 11.11-12.01 interval (section 8.1),
the four equivalent widths move the crossing 11.73 -> 11.64 / 11.70 (11.11
rung) and 11.54 / 11.77 (2.09 rung): **-0.19 to +0.05 in He/H, under a tenth of
the 10.92-12.58 the measurement's own 1 sigma band allows**. It is also
smaller than the workflow's own path spread: re-solving the *same* profile
from a different seed moves the equivalent width by 1.5% (11.11 rung) and 0.31%
(2.09 rung), against the 0.10-0.61% the eddy form produces.

### Two limits on what this tested

1. **The deep half of the law is wrong, and knowingly so.** `P^-1/2` is the
   saturated gravity-wave scaling of a *stratified* region (Lindzen 1981), and
   the climate solve puts the tropopause at 2.636 bar; below it the column is
   convective, where the Gierasch-Conrath mixing-length form applies and
   `K_zz` should flatten or grow with depth. The amplitude is not defensible
   there either: 1.05e5 at 1 bar is at or below the floor of the modelling
   literature — 1.5-2.5 decades under Charnay et al. (2015) 3e6-3e7 for
   GJ 1214 b, 3.7 decades under Parmentier et al. (2013) 5e8 for HD 209458 b,
   and equal to the 1e5 minimum Ackerman & Marley (2001) impose by hand. The
   correct construction is two-branch, a wave law above the
   radiative-convective boundary joined to a mixing-length law below it; this
   test ran one branch to 16.73 bar. *Read as interpretation*: that is very
   likely why the carriers moved as far as they did, since the over-quenched
   0.1-10 bar interval is the one that fixes where CH4 and NH3 quench. **The
   defect does not reach the verdict**: a form this wrong in the deep column
   still failed to move an elemental ratio at the match by more than 1.2%.
2. **Test B, chemistry and wind, is a bounded test of the wind side.** The column reaches only
   1.0091 R_p above the match, and every cell above that takes the top level's
   value, so the wind saw a changed coefficient in a shell nine thousandths of
   a radius deep — 9.5x lower at the base, equal at the top — while the
   homopause under the constant sits at 1.056 R_p (section 3), outside it.
   Test B is not a test of lowering `K_zz` where helium settles; the unbounded
   version of that is the magnitude scan of section 3. Anchoring the same law
   at 1 microbar instead would raise the inherited value by 9.5x over the whole
   domain — a change of magnitude, a separate test, not run.

**One point in the law's favor.** Its value at the match, 1.05e8, falls inside
the 2.8e7-1.4e8 that the Arfaux & Lavvas homopause scaling gives when it is
evaluated with *this* planet's `g` and `T` (section 2) — the only
planet-specific estimate in the survey, and the direction in which section 5
already records the adopted 1e9 as having no LHS 1140 b basis.

### Verdict

**The constant assumption was carrying a great deal of the lower atmosphere's
carrier chemistry — up to four decades — and none of the conclusions of this
document.** EXHALE consumes elemental ratios and thermodynamics from the
profile; both agree at the match to 1.2% under a form that differs from the
constant by 9.5x there and four decades at the bottom, because the match is
still well mixed under either (1.05e8 is 200x the molecular coefficient
~4.6e5 cm^2/s at the base, section 1). **What the results carry is the
magnitude of `K_zz` — the 1e8-1e10 bracket, across which the equivalent width
moves by a factor 9 (section 3) — not its pressure form.**

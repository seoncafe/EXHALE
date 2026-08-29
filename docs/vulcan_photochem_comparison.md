# VULCAN and Photochem as the EXHALE lower-atmosphere pre-step

- Date: 2026-08-09
- Planet: HD 189733 b
- Codes: VULCAN (Tsai et al. 2017, 2021) and Photochem (Wogan et al. 2023, 2024)
- Scripts: `docs/compare_vulcan_photochem.py` (`run` / `base` / `compare`)
- Question: EXHALE's Tier-3 pre-step currently gets its base composition from
  VULCAN. Would Photochem give a different answer, and would EXHALE notice?

## Summary

Photochem and VULCAN disagree about the base composition by a large factor —
the neutral-hydrogen fraction at 1 ubar differs by 7.2x — and most of that
disagreement is the reaction network rather than the code. But **none of it
reaches EXHALE**: the `base.inp` interface carries four numbers, and the
hydrogen dissociation state is not one of them. Written from either code, the
handoff files agree to 0.02%.

So the choice of photochemistry code is, at present, not a question about
EXHALE's answer. It is a question about the interface.

> [2026-08-27] The two sentences above and the "four numbers" of this summary
> are the 2026-08-09 state and no longer hold. The interface changed twice --
> `q_H2_base` and `p_base` became read keys on 2026-08-10 (see "Where this
> leaves the coupling"), the elemental reservoir keys `<El>_H_base` followed,
> and since 2026-08-27 the handoff can be a *profile* file rather than
> scalars (`docs/phase_e_flux_closure_design.md`,
> `src/modules/files_IO/lower_atmosphere_profile.f90`). The code choice was
> taken with it: see P1.7 below.

## Setup

Everything that can be held fixed between the two codes was:

| held identical | value / source |
|---|---|
| T(p), Kzz | `VULCAN/atm/atm_HD189_Kzz.txt` — the same file |
| stellar spectrum | `VULCAN/atm/stellar_flux/sflux-HD189_Moses11.txt` |
| gravity | R_p = 1.138 R_J, g_s = 2140 cm/s^2 (M_p back-solved so g matches) |
| elemental abundances | Lodders, from `fastchem_vulcan/input/solar_element_abundances.dat` |
| zenith angle | 48 deg |

Not removable: the reaction networks, the radiative transfer, and the
integration schemes.

Three runs:

| | code | network |
|---|---|---|
| A | VULCAN | `NCHO_photo_network` |
| B | Photochem | Zahnle set restricted to H/He/N/O/C |
| C | Photochem | `NCHO_photo_network`, converted with `photochem.utils.vulcan2yaml` |

C exists so that A vs C isolates the code and B vs C isolates the network.

Run A reproduces the earlier VULCAN result recorded in
`lower_atmosphere_coupling.md` (q_H2 = 0.63, q_H = 0.23, T_base = 863 K, r_base
1.168 R_J) to about 2%, although that run used the SNCHO network — a check that
this configuration is the same one.

## Composition at 1 ubar

| | T [K] | q_H2 | q_H | q_He |
|---|---|---|---|---|
| A VULCAN / NCHO | 863.4 | 0.618 | **0.238** | 0.143 |
| B Photochem / Zahnle | 863.6 | 0.806 | **0.033** | 0.160 |
| C Photochem / NCHO | 864.2 | 0.717 | **0.130** | 0.152 |

The three temperatures agree to 0.8 K, confirming all three read the same T(p).

Splitting the disagreement in q_H:

- **network** (B to C, same code): 0.033 to 0.130, a factor **4.0**
- **code** (C to A, same chemistry): 0.130 to 0.238, a factor **1.8**
- combined: **7.2**

The network dominates. Photochem's own set is built for rocky-planet
atmospheres and appears to carry fewer of the H2 dissociation paths that matter
at 860 K; that reading is inferred from the result, not from a reaction-level
audit.

Deeper down the codes agree well. At 1000 ubar all three give H2O and CO to
within 10% and H2/He to within 2%. The disagreement is confined to the H2
dissociation region above about 10 ubar — which is exactly where EXHALE takes
its base.

## Runtime

| | steps | wall time |
|---|---|---|
| A VULCAN | ~4300 | 1370 s |
| B Photochem / Zahnle | ~600 | < 60 s |
| C Photochem / NCHO | ~700 | < 90 s |

Roughly 15-20x. This is observed wall time on a shared machine, not a
controlled benchmark: VULCAN is Python and single-threaded, Photochem is
Fortran over CVODE.

## What the handoff actually carries

`read_base_inp` (in `src/modules/files_IO/input_read.f90`) read exactly four
keys at the time of this comparison:

```
T_base -> T0        r_base -> R0        HeH_base -> HeH        Kzz_base -> he_kzz
```

`q_H2` and `q_H` are written by `src/utils/vulcan_to_base.py` as **comments**.
They do not reach the code.

*[2026-08-15: `read_base_inp` now reads six keys — the four above plus
`q_H2_base` (the photochemical H2 volume mixing ratio at the base, which drives
the molecular-base particle count) and `p_base` (the pressure level the handoff
describes, default 1 ubar). So `q_H2` does reach the code now, under the name
`q_H2_base`; the "the interface discards it" conclusion below was true when
this memo was written and no longer is. The rest of the comparison — that
T_base, r_base and HeH_base agree to 0.02% across the three photochemistry
solutions because they are set by inputs rather than results — is unaffected.]* Running the same hypsometric integration on all
three solutions:

| handoff | T_base [K] | r_base [R_J] | HeH_base | (q_H, dropped) |
|---|---|---|---|---|
| VULCAN / NCHO | 863.45 | 1.16844 | 0.096939 | 0.238 |
| Photochem / Zahnle | 863.64 | 1.16821 | 0.096951 | 0.033 |
| Photochem / NCHO | 864.16 | 1.16829 | 0.096950 | 0.130 |

**The four numbers EXHALE reads agree to 0.02%.** They agree because all three
are set by things that are inputs rather than results: T_base by the prescribed
T(p), HeH_base by the elemental abundance, r_base by an integration over mu(p)
that the H2/H split barely moves (mu changes by 10% while q_H changes by 7.2x,
because dissociating H2 into 2 H conserves mass).

The one quantity the photochemistry actually determines — and the one the
Tier-1 caveat was about, that a hot Jupiter's atomic base is made by
photochemistry rather than by temperature — is the one the interface discards.

## Running EXHALE from each handoff

Both handoffs were run through EXHALE from an identical initial condition (the
converged HD 189733 b state), same `input.inp`, same executable, differing only
in `base.inp`.

This did not resolve anything, for a reason worth recording. HD 189733 b does
not reach a steady state: over 87000 steps `du` wandered between 0.77 and 1.86
and never approached 1e-3 — the marching-time base breathing that HD 189733 b
has shown since ATES. (This originally cited `TO_BE_DONE.md` item (A); note
that item is now closed, and its diagnosis of the *JFNK* residual floor as a
base-momentum wall was refuted — `docs/newton_scaling_and_base_wall.md`. The
marching oscillation described here is a separate observation and stands.)
Cutting both runs at a fixed 20000 steps
(`EXHALE_MAXSTEPS`) and comparing on common physical radii gives median
differences of 0.1-2% and maxima of tens of percent concentrated in the
innermost cells at r ~ 1.17 R_J.

Those maxima are not a measurement of the handoff. A 0.02% shift in `r_base`
changes the phase of the base oscillation, so at a matched step the two runs
sit at different points of the same cycle. The medians bound the bulk effect;
the maxima measure the sensitivity of an unconverged transient.

A conclusive EXHALE-level comparison needs a planet that converges. It was not
run here because it would need new VULCAN and Photochem runs for that planet's
T(p) and spectrum.

## Traps met along the way

Each of these cost a run, and none of them announces itself.

1. **VULCAN `atom_list` must match the network's atoms.** Leaving `S` in
   `atom_list` while running the sulfur-free `NCHO_photo_network` produced a run
   that did 8017 steps in 66 minutes while advancing the simulated time to
   1e-10 s: `dt` collapsed to 2e-14 and stuck, `longdy` pinned at 1.0, the
   reported worst cell `nz = 0` and species `OH`, and the log filling with
   "Element conservation is violated too large". The corrected run finished in
   1370 s with zero conservation warnings. **The step counter is not the health
   indicator — the simulated elapsed time is.**

2. **`use_solar = True` makes VULCAN ignore the `O_H`/`C_H`/`N_H`/`He_H` lines
   of `vulcan_cfg.py`** and read `fastchem_vulcan/input/solar_element_abundances.dat`
   instead. Taking the config values as the abundances put He/H at 0.0838
   against VULCAN's actual 0.0969, and that mismatch alone moved `HeH_base` by
   16% — larger than any real difference between the codes, and in the one
   quantity the handoff is sensitive to. Read the table, not the config.

3. **`vulcan2yaml` emits a degenerate `He <=> He` reaction** (A = 0) because
   VULCAN's networks list He as a non-reacting bulk species. Photochem rejects
   the mechanism with `IOError: This reaction is a duplicate`. Dropping the
   entry changes no rate; `compare_vulcan_photochem.py` does it automatically.

4. **Photochem's grid must extend above the supplied climate grid.** A T(p)
   profile reaching higher than Photochem's top of atmosphere raises
   `The photochemical grid needs to extend above the climate grid`. Truncating
   the profile is the easy fix but it moves the upper boundary, which is why the
   1.8x "code" factor above should be read as an upper bound: VULCAN's domain
   top sits two decades above the comparison point and Photochem's one decade.

5. **`Load IC? True` needs `output/*_IC.txt` to exist already.** The copy from
   `Hydro_ioniz.txt` to `Hydro_ioniz_IC.txt` is done by the run script, not by
   `EXHALE.x`. Supplying only `Hydro_ioniz.txt` makes the Fortran `open`
   create an empty `_IC` file, and the run dies in `load_IC.f90` with
   `End of file`.

6. The step-cap environment variable is **`EXHALE_MAXSTEPS`**, not
   `ATES_MAXSTEPS`.

## Where this leaves the coupling

Switching the pre-step to Photochem is not urgent: with the interface as it
stands, EXHALE would not see the change. Photochem is faster and carries a
gas-giant extension (`photochem.extensions.gasgiants`, with quench-based
initial conditions and automatic boundary management) plus a climate model
(`clima`), which would matter the day the prescribed Guillot T(p) in
`run_lower.py` is replaced by a self-consistent one.

The more useful change is to the handoff itself. If `base.inp` carried `q_H2`
and `q_H` as read keys rather than comments, the pre-step would deliver what it
actually computes, the Tier-1 caveat would be quantified inside the model
rather than in a document, and the difference between photochemistry codes
would become measurable in EXHALE output. That is a change to `read_base_inp`
and to whatever consumes the base composition, and it has not been scoped here.

Neither code provides atomic-metal release (Na/Mg/Ca/Fe), which EXHALE needs
for the metal lines; that part of the base stays user-supplied.

**2026-08-10 — the interface change described above is implemented.**
`read_base_inp` now reads `q_H2_base` (the photochemical H2 volume mixing ratio
at the base) and `p_base` (the level it refers to), and both converters write
them. With `Molecular base: True` the photochemical value replaces the
chemical-equilibrium fit in the base particle count, so the difference between
networks now reaches the wind solve; `base.inp` files without the key behave
exactly as before. Design, size of the effect and the validation gates:
`docs/base_composition_handoff_plan.md`. `q_H` and the molecular mixing ratios
stay comments — `q_H` is implied by `q_H2_base` and `HeH_base`, and the
molecules have nothing to act on in EXHALE's atomic metal set.

**2026-08-19 — where the code choice is decided.** Because `q_H2_base` is now
read, the 7.2x spread measured here propagates into the wind solve, so choosing
between the two codes is no longer neutral. The decision, together with the
reasons to prefer Photochem that have nothing to do with chemistry (speed, the
`gasgiants` extension, and `clima` as the only path away from a prescribed
T(p)), is in `docs/oxygen_chemistry_options.md` §2.4. The recommendation there
is to fix the reaction network before changing codes, since the network carries
4.0x of the 7.2x. The plan of record for acting on this, including the matched
network-plus-domain rerun (its phase P1), is `docs/oxygen_chemistry_new_plan.md`.

---

# 2026-08-26 — Phase P1: matched network, matched domain, and the H2 budget

This section is the record of phase **P1** of `docs/oxygen_chemistry_new_plan.md`
(reinstall Photochem, rerun the comparison with matched network *and* matched
vertical domain, add Photochem's gas-giant mechanism as a third arm, export the
reaction budget at the handoff level, and measure the EXHALE-level effect on a
planet that converges). Everything below was measured on 2026-08-26 in this
working copy. The 2026-08-09 sections above are left unedited; where a number
here supersedes one there, it is said so explicitly.

Script: `docs/p1_matched_comparison.py` (`run` / `compare` / `budget` / `base`),
which imports the shared helpers of `docs/compare_vulcan_photochem.py` so the
2026-08-09 procedure is not duplicated. Run directories:
`vulcan_work/pc_compare_p1/`.

## P1.0 The Photochem environment

The plan recorded that "Photochem is not importable from any Python on this
machine". As of 2026-08-26 that was true of the default interpreter and false
of the conda environment built for the 2026-08-09 comparison, which still
exists and still works. (It is no longer true of the default interpreter
either: the corrected photochem 0.9.0 build was installed there on 2026-08-29,
`README_photochem.md`. That is a *different* build from the 0.8.4 this section
is about, and the two are not interchangeable.) **No reinstall was needed**; the environment was
verified by importing the package, the `extensions.gasgiants` and `utils`
submodules, and by running seven photochemical models to steady state through
it.

| item | value |
|---|---|
| interpreter | a dedicated conda environment, Python 3.11.15 |
| `photochem` | 0.8.4 (conda-forge, `py311h57bc489_0`) |
| `photochem_clima_data` | 0.3.1 (`pyhd8ed1ab_0`) |
| data directory | `.../site-packages/photochem_clima_data/data` |
| `zahnle_earth.yaml` | 136018 bytes, sha256 `bb72b9ac9fc43104...` |
| other | numba 0.66.0, numpy 2.4.6, scipy 1.17.1 |

Reproduction from scratch, if the environment is ever lost:
`conda create -n <name> -c conda-forge photochem` (the data package comes
in as a dependency). The workspace also carries a source clone at
`photochem/` (tag `v0.9.0`), which is **not** what is imported and was used only
as a reference for the gas-giant workflow.

## P1.1 What changed relative to the 2026-08-09 comparison

1. **Domain.** Photochem's top of atmosphere is set to VULCAN's own
   `P_t = 1e-2` dyn/cm^2 (`gdat.TOA_pressure_avg`), instead of leaving it at the
   1e-7 bar default with the climate grid cut at 0.5 dyn/cm^2. `gasgiants`
   requires `3*TOA < P_climate_top`, so the supplied T(p) is cut at
   3.05*TOA rather than at a fixed pressure.
2. **A third arm**, `zahnle_s`: the same Zahnle set with sulfur, i.e. the
   H/He/N/O/C/S gas-giant mechanism that `photochem.extensions.gasgiants`
   generates.
3. **Reaction budgets** exported (`EvoAtmosphere.production_and_loss`) for H2,
   H2O, OH and H on the model grid.
4. **A second planet**, HD 209458 b, run through all three Photochem arms
   against the existing VULCAN solution, because HD 189733 b does not reach a
   steady state in EXHALE and therefore cannot carry the wind-level measurement.
5. **VULCAN reruns with a truncated domain**, so the domain sensitivity is
   measured in both codes rather than in one.

Held fixed exactly as in 2026-08-09: T(p), Kzz, the stellar spectrum and its
(R_star/a)^2 dilution, gravity, the Lodders elemental abundances VULCAN actually
reads, and the 48 deg zenith angle.

## P1.2 HD 189733 b at 1 microbar, all arms

| run | code | network | model top [dyn/cm^2] | T [K] | q_H2 | q_H |
|---|---|---|---|---|---|---|
| A | VULCAN | NCHO | 1.00e-2 | 863.4 | 0.6181 | **0.2382** |
| B | Photochem | Zahnle H/He/N/O/C | 1.25e-1 | 863.6 | 0.8063 | **0.0331** |
| C | Photochem | NCHO (vulcan2yaml) | 1.61e-1 | 864.2 | 0.7170 | **0.1304** |
| B' | Photochem | Zahnle H/He/N/O/C | 9.14e-3 | 863.4 | 0.8041 | **0.0354** |
| C' | Photochem | NCHO (vulcan2yaml) | 6.71e-3 | 863.4 | 0.7081 | **0.1401** |
| D' | Photochem | Zahnle H/He/N/O/C/S | 9.10e-3 | 863.4 | 0.7997 | **0.0403** |
| C'' | Photochem | NCHO | 4.65e-3 | 863.5 | 0.7016 | **0.1471** |

B and C are the 2026-08-09 runs, reused unchanged; the primed runs are the
matched-domain repeats. All seven temperatures agree to 0.8 K.

**The decomposition the P1 gate asks for**, at matched model top (~1e-2
dyn/cm^2):

| factor | ratio in q_H | measured as |
|---|---|---|
| domain truncation | **1.07** | C -> C' (same code, same network) |
| reaction network | **3.95** | B' -> C' (same code, same domain) |
| code | **1.70** | C' -> A (same network, same domain) |
| product | 7.20 | = A/B, the 2026-08-09 total |

This supersedes the 2026-08-09 split (network 4.0, code 1.8): the network share
is essentially unchanged, the code share drops from 1.8 to 1.70, and the
remaining 1.07 was the truncation artifact. **The artifact is small**, which was
not the expectation the plan was written with.

## P1.3 Domain sensitivity is itself code-dependent

The 2026-08-09 caveat said the 1.8x code factor "should be read as an upper
bound" because Photochem's domain was the shorter one. Measured in both codes,
the sensitivity is strongly asymmetric:

| code | top moved from | to | q_H at 1 ubar | ratio |
|---|---|---|---|---|
| Photochem / NCHO | 1.61e-1 | 6.71e-3 | 0.1304 -> 0.1401 | 1.07 |
| Photochem / NCHO | 6.71e-3 | 4.65e-3 | 0.1401 -> 0.1471 | 1.05 |
| VULCAN / NCHO | 1.00e-2 | 1.61e-1 | 0.2382 -> 0.2663 | 1.12 |
| VULCAN / NCHO | 1.00e-2 | 5.0e-1 | 0.2382 -> 0.4042 | 1.70 |

The two VULCAN rows are fresh runs of the same configuration with `P_t` changed
(`vulcan_work/pc_compare_p1/vulcan_hd189_trunc016/` and `.../vulcan_hd189_trunc/`,
4821 and 1985 steps to steady state, no element-conservation warnings). The third row truncates VULCAN
to the top the 2026-08-09 Photochem runs actually had.

Facts: for the *same* truncation (top at 1.61e-1 dyn/cm^2) the two codes respond
by 1.12 and 1.07, in opposite directions — VULCAN gives more atomic H when
truncated, Photochem slightly less. Pushing VULCAN's top down further, to
5.0e-1, raises q_H by 1.70. Interpretation, not demonstrated here: truncation
removes the shielding column above the comparison level, which raises the local
photolysis rate, and it also removes the atomic H produced in the hot top that
would mix downward; the two have opposite signs and the codes weight them
differently.

**Practical consequence.** Matched by truncating VULCAN to the old Photochem top
instead of by extending Photochem, the "code" factor would have read 2.04
(0.2663/0.1304) rather than 1.70. The domain statement must accompany any
quoted code-to-code ratio, and the 2026-08-09 reading that the code factor was
an upper bound was the wrong way round: extending the short domain lowered it,
truncating the long one would have raised it.

**A limit worth recording.** Photochem cannot be given VULCAN's full T(p): with
the top three levels kept (T up to 6000 K) `set_temperature` raises
`The temperature is not within the ranges given for the thermodynamic data for
NH2`. Run C' therefore holds T constant at 1454 K above 3.29e-2 dyn/cm^2 where
VULCAN continues to 6000 K, and C'' (top of the thermal profile at 1911 K)
shows that pushing that limit further moves q_H by another 5%. The residual
"code" factor of 1.70 contains this thermal-top difference; it is not purely
numerics.

## P1.4 The reaction budget at the handoff level

`production_and_loss` for H2 at 1 microbar, HD 189733 b. Rates are cm^-3 s^-1;
"O family" means every channel with O, O(1D), OH or H2O as a reactant or
product.

| arm | gross loss | gross production | net loss | O family, gross loss | O family, net |
|---|---|---|---|---|---|
| C' Photochem/NCHO | 1.470e7 | 6.087e6 | 8.617e6 | 99.1% | **99.0%** |
| B' Photochem/Zahnle | 4.235e6 | 1.347e6 | 2.888e6 | 98.8% | **99.6%** |
| D' Photochem/Zahnle+S | 7.271e6 | 4.419e6 | 2.851e6 | 59.0% | **99.3%** |

The single largest channel in every arm is `OH + H2 -> H2O + H` (90.2% of gross
loss in C', 95.9% in B', 57.2% in D'), and the water made that way is returned
to OH mostly by photolysis (`H2O + hv -> OH + H`, 54% of H2O loss in C') and by
`H2O + H -> OH + H2`. In D' the sulfur channel `S + H2 -> H + HS` carries 40% of
the gross loss but is almost exactly balanced by its reverse
(`H + HS -> S + H2`, 2.933e6 against 2.943e6), so it contributes 0.4% of the
net.

**Gate answer for HD 189733 b: yes.** The H2O/OH cycle dominates H2 destruction
at the handoff level — 99% of the net in all three arms, independent of which
network is used.

The same budget on **HD 209458 b** reads differently, and the difference
explains the rest of this section:

| arm, level | gross loss | net loss | O family, gross | O family, net |
|---|---|---|---|---|
| C' at 1 ubar (T = 2331 K) | 8.367e9 | 6.089e6 | 99.9% | **29.6%** |
| B' at 1 ubar | 7.606e9 | 5.889e6 | 99.9% | 25.4% |
| D' at 1 ubar | 1.053e10 | 6.507e6 | 95.3% | 20.3% |
| C' at 1e-4 bar (T = 1830 K) | 2.565e12 | 1.306e8 | 100.0% | **2.6%** |

At 2331 K the oxygen channels run three orders of magnitude faster than the net,
in near-exact balance both ways, and the net H2 destruction is carried by
`H2 + M -> H + H + M` (4.30e6 of the 6.09e6 net in C'). **On this planet the
base partition is thermal, not photochemical.**

Chemical timescales at the handoff, `n_H2 / (net H2 loss)`:

| case | tau_chem(H2) |
|---|---|
| HD 189733 b, 1 ubar, 864 K | 6.0e5 s |
| HD 209458 b, 1 ubar, 2331 K | 2.2e5 s |
| HD 209458 b, 1e-4 bar, 1830 K | 2.3e6 s |

Against EXHALE's own advection time at its base cell in the converged
HD 209458 b run below (`H/v` with `v = 0.52` cm/s, `H = 1.27e8` cm, i.e.
2.4e8 s), `tau_chem/tau_adv ~ 1e-3`. That is the P4 gate's second condition, and
it is met at the handoff level — it says nothing about the cells further out.

## P1.5 HD 209458 b: the three arms agree, and why

Same three arms, the VULCAN solution being the converged
`vulcan_work/hd209_vulcan/output/HD209.vul` of 2026-08-10 (no rerun).

| level | A VULCAN/NCHO | C' pc/NCHO | B' pc/Zahnle | D' pc/Zahnle+S |
|---|---|---|---|---|
| q_H2 at 1e-4 bar | 0.8020 | 0.8027 | 0.8028 | 0.7757 |
| q_H2 at 1e-6 bar | 0.4237 | 0.4536 | 0.4626 | 0.3765 |
| q_H at 1e-6 bar | 0.4499 | 0.4171 | 0.4073 | 0.5011 |

The network factor that reaches 4.0 on HD 189733 b is **1.02** here
(B' vs C' in q_H), and the code factor is **1.08** (C' vs A). The reading that
fits both planets, and that the budget of P1.4 supports: where the base is cool
enough for the partition to be kinetic (864 K), the reaction network decides it;
where the base is hot enough for collisional dissociation to run the net budget
(2331 K), all networks converge on the same answer. Sulfur is the exception —
it moves q_H2 by 17% at 1e-6 bar and 3% at 1e-4 bar, in the direction of *less*
H2.

## P1.6 The EXHALE-level effect on a converged planet

Run directories `vulcan_work/pc_compare_p1/exhale_hd209/p1e-6_seed/<arm>/`. All
five runs are the same `HD209458b/input.inp` with `Molecular base: True`,
the same `metals.inp`, the same initial condition (the converged state in
`vulcan_work/hd209_wind_response/photo_newton/output`), the same executable and
the same 4 OpenMP threads; they differ in `base.inp` only. Each arm's `base.inp`
is written from its own solution, so `T_base`, `r_base` and `HeH_base` are that
arm's own values as well (they agree to 0.06%, as in 2026-08-09).

| arm | q_H2_base | ntot_bc | JFNK | ghost T [K] | log10 Mdot [g/s] | vs VULCAN |
|---|---|---|---|---|---|---|
| VULCAN / NCHO | 0.4237 | 0.703 | info=0, \|\|R\|\| 9.7e-4 | 1638.7 | **9.8853** | — |
| Photochem / NCHO | 0.4536 | 0.689 | info=0, \|\|R\|\| 3.6e-4 | 1605.0 | **9.8788** | -1.5% |
| Photochem / Zahnle | 0.4626 | 0.684 | info=0, \|\|R\|\| 1.4e-4 | 1595.2 | **9.8792** | -1.4% |
| Photochem / Zahnle+S | 0.3765 | 0.727 | info=0, \|\|R\|\| 8.2e-4 | 1694.9 | **9.9130** | +6.6% |
| chem.-eq. fit (no `q_H2_base`) | 0.0063 | 0.994 | info=2 (failed) | — | not converged | — |

log10 Mdot is `exhale_io.mdot_log10` on the `_adv` profiles; the code's own
printed values are 9.94 / 9.93 / 9.93 / 9.97 and give the same differences.

**Measured: the choice of photochemistry code and network moves HD 209458 b's
mass-loss rate by at most 0.035 dex (8.3%) across the four arms, and by 1.5%
between the two codes on the same network** — all four JFNK-converged. The
deeper-handoff set below shows that this configuration's own JFNK-to-JFNK
spread is about 0.025 dex, so 0.035 dex is an upper bound on the effect rather
than a resolved signal. The ordering follows `ntot_bc` exactly, as §11.5 above
predicted: less H2 at the base means more particles, a hotter ghost cell and a
stronger wind.

For contrast, the equilibrium-fit reference is a much larger perturbation
(`ntot_bc` 0.994 against 0.68-0.73) and its run did not converge: JFNK returned
`info=2` and the resumed marching was still at `du = 3.2` after 7000 steps. The
"photochemical against equilibrium" difference recorded in
`base_composition_handoff_plan.md` §11 (-12.9% at a matched marching state)
therefore remains the un-Newton-finished number it was; **only the
arm-against-arm comparison here is JFNK-converged.**

### The same measurement at the deeper handoff, and what it says about the noise

`p1e-4_seed/`, identical in construction but at `p_base = 1e-4` bar (T = 1850 K),
where the arms are physically almost the same solution: q_H2 = 0.8020 (VULCAN),
0.8027 (Photochem/NCHO), 0.8028 (Photochem/Zahnle), 0.7757 (Photochem/Zahnle+S).
All five runs reached `JFNK info=0`.

| arm | q_H2_base | ntot_bc | log10 Mdot | vs VULCAN |
|---|---|---|---|---|
| VULCAN / NCHO | 0.8020 | 0.556 | 9.7517 | — |
| Photochem / NCHO | 0.8027 | 0.555 | 9.7764 | +5.9% |
| Photochem / Zahnle | 0.8028 | 0.555 | 9.7764 | +5.9% |
| Photochem / Zahnle+S | 0.7757 | 0.564 | 9.7778 | +6.2% |
| chem.-eq. fit | 0.7940 | 0.558 | 9.7649 | +3.1% |

**A 0.1% difference in the input produces 5.9% in Mdot here, and the ordering
does not follow `ntot_bc`.** Read as a noise measurement, this configuration
reproduces log10 Mdot to about **0.025 dex** between JFNK solutions that started
from the same IC and differ negligibly in their boundary data. The honest
conclusion for §P1.6 is therefore: the arm-to-arm effect at the 1 microbar
handoff (0.035 dex spread, ordered by `ntot_bc`) is only marginally above that
floor, and the safe statement is **"the choice of photochemistry code or network
moves HD 209458 b's mass loss by at most a few percent, comparable to the
spread between JFNK solutions of the same configuration"** — not a clean
signal of 8.3%.

Two further sets were launched from the cold isothermal IC (`p1e-6/`, `p1e-4/`)
as an independent path check. They were still marching down when the session
ended and were stopped: `du` was 0.40 (all four 1 microbar arms; 1.17 for the
equilibrium-fit reference) at 71500 steps, and 0.33 (all five) at 71000 steps at
1e-4 bar. Nothing above uses them. The equilibrium-fit reference at 1 microbar
was also stopped: after its JFNK failure it sat between `du` = 0.076 and 0.096
for 70000 further steps, i.e. it neither converged nor diverged.

## P1.7 Evidence for the code choice — no decision taken

The plan leaves the VULCAN-or-Photochem decision to the user. What P1 measured:

> [decision taken 2026-08-27] Photochem is the production chemistry
> (`src/utils/photochem_to_lower_profile.py`) and VULCAN the cross-check arm
> (`src/utils/vulcan_to_lower_profile.py`), on the last row of the table
> below: only Photochem carries a climate model, so only it can be the route
> away from a prescribed T(p). Record: `docs/Update_EXHALE.md` sections 77-78,
> `docs/oxygen_chemistry_new_plan.md` ("P1 and P2 on a real target").
> The evidence table is left as measured.

| criterion | VULCAN | Photochem |
|---|---|---|
| installation | in-tree copy (`EXHALE_v1.00/VULCAN`), fetched by `src/utils/setup_vulcan.sh`; pure Python, no build | one conda package plus its data package; not importable from the default interpreter, so a dedicated environment is required |
| reproducibility of a run | config file edited in place; `chem_funs.py` regenerated per network | mechanism and thermodynamic files written by the script; data package versioned |
| runtime (HD 189733 b, this session) | 629 s / 1985 steps and 1593 s / 4821 steps (the two truncated domains, VULCAN's own CPU-time line); 1370 s / ~4300 steps (2026-08-09, full domain) | under 2 minutes per arm, all reaching steady state (wall clock from the log timestamps on a shared machine, not a controlled benchmark) |
| network flexibility | its own NCHO/SNCHO sets, plain text | its own Zahnle sets by atom list, plus any VULCAN network through `vulcan2yaml` |
| domain | any `P_t`; took VULCAN's own 6000 K thermospheric top | thermodynamic data limit the top to T well below 6000 K (NH2 failed); the model top must sit below the climate grid top by 3x |
| q_H agreement on the same network | reference | 1.70x on HD 189733 b (864 K base), 1.08x on HD 209458 b (2331 K base) |
| effect on EXHALE Mdot (HD 209458 b, converged) | reference | -1.5% on the same network; +6.6% with the sulfur mechanism |
| climate model | prescribed T(p) only | `clima` available, the only path away from a prescribed T(p) |
| metal release (Na/Mg/Ca/Fe) | not provided | not provided |

## P1.8 Traps met in P1

1. **`gasgiants` needs `3*TOA_pressure_avg < P_climate_top`.** Setting the top
   pressure without also cutting the supplied T(p) raises
   `The photochemical grid needs to extend above the climate grid`.
2. **Photochem's thermodynamic data cap the temperature.** VULCAN's HD 189733 b
   profile ends at 6000 K and Photochem rejects it (`NH2`). The matched-domain
   runs therefore match the top *pressure*, not the top *temperature*.
3. **An EXHALE run directory must contain `output/` before the run starts.**
   Without it `set_IC.f90` dies at line 201 with
   `Cannot open file 'output/IC_dump.txt'`, after the setup report has already
   been written — which makes it look like a physics failure.
4. `EvoAtmosphereGasGiant` orders `gdat.gas.atoms_names` differently from run to
   run; the elemental abundances must be mapped by name, never by position.

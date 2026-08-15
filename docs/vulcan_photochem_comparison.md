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

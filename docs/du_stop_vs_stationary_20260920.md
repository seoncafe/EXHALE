# A flux-converged state that this code refuses, against the certified stationary solution

2026-09-20. Case `LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13`,
binary `LHS1140b/models/EXHALE_75d55d9d.x`, md5
`75d55d9d4fd0e748cd01d6e35713e34c`, on `lart4` at `OMP_NUM_THREADS=8`.
Every number below is MEASURED in this session unless it is marked READ.

## The question

`docs/lhs1140b_catalog_state_20260920.tex`, section 3, states that several
models of the catalog meet the flux-spread criterion the literature converges
on, `du < 1e-3` (ATES, Caldiroli et al. 2021; CETIMB, Koskinen et al. 2013a,
who ask that `F_c = rho v r^2` be constant with altitude), and are
nonetheless refused by this code's row-by-row certification. The question
this document answers is what such a state costs: on ONE CERTIFIED case, how
far does a flux-converged but refused state sit from the certified stationary
solution, as a function of radius, and does the difference reach a measured
quantity.

## What could not be done, and why

**The marching route does not reach `du < 1e-3` on this case.** The recipe of
`README_HOWTO.md` (`Reconstruction scheme: PLM+WENO3`, `du_th [PLM,WENO3]:
0.5 1.0e-3`, no Newton finish, no `Resid tol`) was run twice on scratch
copies of this case's configuration, with the run's own `du` (the spread of
`rho v r^2` over `r >= Escape radius = 2 R_p`, `mass_flux_spread` in
`src/EXHALE_main.f90`) read from its log:

| start | steps | `du` at the first step | smallest `du` reached | `du` at the last step | rate |
|---|---|---|---|---|---|
| cold, the code's own isothermal initial condition | 2106 | 9.4544 | 5.5972 (the last step) | 5.5972 | 2.07 steps/s |
| the certified state of `HeH1.70`, the nearest certified case of the same group, mapped onto this grid with its reservoir carried to He/H = 2.13 | 7627 | 2.7911E-03 (step 2) | 2.7911E-03 (step 2) | 45.85 | 8.31 steps/s |

The cold march is descending at the point it was stopped and is still three
and a half decades above the threshold. The seeded march starts one step away
from the threshold, turns around, and climbs to 218 before oscillating
between 45 and 190. This reproduces what was already READ in `LHS1140b/MODELS.md`
section 6 for this configuration: "`du` does not settle: 9.45 (step 1) to
8.35 (step 5212) from the cold start; 1.32 (step 499) to 10.7 (step 8147)
from the mapped well-mixed state". The mass-flux spread is taken over a
window that reaches to the 30 R_p outer boundary, where the wind is subsonic
and the relaxation time is that of the whole domain, so a marching run would
have to integrate a flow time of the outer domain to flatten it.

**Taking another certified case of the same group does not change this.**
Every case of `atomic_scalar_gj1132_kzz1e9` carries the same input but the
He/H ratio, so the marching behavior above is the group's, not this case's.

**The initial condition of the certified solve is not on disk.** The case's
`REPRODUCE.md` states that the chain ends at generation
`g0002_20260917T171710Z_4148d706`, which names no parent, and
`provenance/g0001.../provenance_recovered.json` and
`provenance/g0002.../provenance_recovered.json` both record
`confidence: "not established"`. No case of this group carries a
`seed_identity.json`. `models/pick_seed.py atomic_scalar_gj1132_kzz1e9/HeH2.13`
run today returns this case's own certified generation, which makes any march
from it vacuous.

## The two states

The flux-converged but refused state was therefore produced by the route this
catalog actually uses, stopped early: the partitioned stationary route
(`Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`) from the mapped
`HeH1.70` seed above, with the outer loop given ONE pass
(`EXHALE_OUTER_PASSES=1`). Its single pass took 123.64 s and returned
`hydro info=0`; the run then wrote the state with `info = 1`,
`certified=F cert_reason=no_stationary_claim`.

Both states were then measured the same way, with `Restart intent: stationary
evaluate`, which takes no step.

| | certified stationary (S) | flux-converged, refused (A) |
|---|---|---|
| generation | `g0005_20260919T212717Z_f9d878de` | none, written in this session |
| `du`, the run's own flux gate | 2.462E-11 (READ, the case log) | 1.8340E-12 |
| `du`, cell-centered spread of `rho v r^2` over `r >= 2` | 1.9397E-04 | 1.8696E-04 |
| the same over `r >= 1.2` | 2.4578E-04 | 2.4542E-04 |
| `\|\|R\|\|` at the re-evaluation | 1.6115E-08 | 4.9690E-08 |
| `\|\|R\|\|` recorded for the state | 1.700E-08 (READ, the case log) | none recorded |
| hydrodynamic mass row | 6.610E-09 at cell 1, tol 7.9E-09, within | 6.610E-09 at cell 1, tol 7.9E-09, within |
| hydrodynamic momentum row | 3.449E-14, tol 1.0E-08, within | 3.449E-14 at cell 494, tol 1.0E-08, within |
| hydrodynamic energy row | 1.612E-08, tol 1.0E-06, within | 4.969E-08 at cell 1, tol 1.0E-06, within |
| elemental transport He/H partition, gated at `r >= 1.2` | within | **5.209E-05 above 1.0E-05 at cell 288** |
| level balance He 2^3S | within | 2.153E-19, tol 1.0E-06, within |
| eliminated-species closure `System_HeH_TR` | within | 6.992E-17, tol 1.0E-06, within |
| verdict | **CERTIFIED** | **NOT CERTIFIED, 1 refusing entry** |

Both states meet the flux-spread criterion by four to eight decades. The one
row that separates them is the elemental transport of the He/H partition in
cell 288, `r = 1.6951 R_p`, standing a factor 5.2 above its tolerance. Its unrestricted maximum, 7.301E-01 at cell 2, sits below the gate
radius and is reported only, as it is in the certified state too (4.051E-01
at cell 2, READ).

## How far apart they are, as a function of radius

Relative difference `(A - S)/|S|` in per cent, over the 500 physical cells:

| `r/R_p` | `T` | `v` | `n_H` | `x(H II)` | `x(He II)` | `n(He 2^3S)` | `rho v r^2` |
|---|---|---|---|---|---|---|---|
| 1.050 | +0.037 | +0.126 | -0.071 | +0.238 | +0.299 | +0.482 | +0.055 |
| 1.201 | -0.233 | +0.073 | -0.018 | -0.656 | -0.149 | -0.734 | +0.055 |
| 1.499 | -0.556 | +0.058 | -0.002 | -1.078 | -0.354 | -0.956 | +0.055 |
| 2.003 | -0.927 | +0.041 | +0.014 | -1.281 | -0.498 | -1.136 | +0.055 |
| 3.015 | -0.809 | +0.028 | +0.028 | -0.938 | -0.384 | -0.965 | +0.055 |
| 4.978 | -0.674 | +0.036 | +0.019 | -0.990 | -0.251 | -1.537 | +0.055 |
| 7.046 | -0.632 | +0.058 | -0.002 | -0.950 | -0.195 | -1.869 | +0.055 |
| 10.028 | -0.620 | +0.114 | -0.058 | -0.599 | -0.167 | -1.365 | +0.055 |
| 20.104 | -0.579 | +0.239 | -0.183 | -0.252 | -0.038 | -0.584 | +0.055 |
| 29.031 | -0.544 | +0.269 | -0.214 | -0.186 | +0.034 | -0.444 | +0.055 |

The largest difference of each quantity over the whole domain:

| quantity | largest relative difference | at `r/R_p` | certified value there | refused value there |
|---|---|---|---|---|
| temperature | -9.397E-03 | 2.133 | 4.2389E+03 K | 4.1991E+03 K |
| velocity | +2.694E-03 | 29.031 | 1.1993E+05 cm/s | 1.2025E+05 cm/s |
| hydrogen density | -2.143E-03 | 29.031 | 2.7305E+04 cm^-3 | 2.7247E+04 cm^-3 |
| `x(H II)` | -1.300E-02 | 1.872 | 3.9028E-02 | 3.8521E-02 |
| `x(He II)` | -4.979E-03 | 1.986 | 6.2636E-02 | 6.2324E-02 |
| `n(He 2^3S)` | -1.878E-02 | 6.738 | 4.3574E+01 cm^-3 | 4.2756E+01 cm^-3 |
| `rho v r^2` | +1.427E-03 | 1.001 | 5.0045E+11 | 5.0116E+11 |

The differences are nowhere larger than 1.9 per cent, and they sit exactly
where the He I 10830 line is formed. The He 2^3S column of the certified
state has its median at `r = 6.39 R_p` and its 5 to 95 per cent range over
`r = 1.69` to `16.72 R_p`, and over that range `n(He 2^3S)` differs by 0.4 to
1.9 per cent, the temperature by 0.5 to 0.9 per cent and `x(H II)` by 0.2 to
1.3 per cent. The mass flux itself is flat to 0.055 per cent between the two,
which is why the flux-spread criterion cannot see the difference at all.

Figures: `docs/figures/du_stop_vs_stationary_profiles.pdf` (the seven
profiles overlaid) and `docs/figures/du_stop_vs_stationary_differences.pdf`
(the relative differences).

## What it costs in a measured quantity

`EXHALE_transit.py` through the WINERED HIRES-Y kernel of the measurement
(`LHS1140b/winered_hires_y.sh`, `R = 68000`), run on both states with each of
the two state selections. The `adv` selection is the catalog's own and
reproduces its published number for this case (EW 1.3749 %A, READ from
`docs/lhs1140b_catalog_state_20260920.tex`); the `solution` selection reads
the state files themselves. The comparison is like for like in both columns,
the same selection on both states, and the two selections disagree with each
other far more than the two states do, since the advection-corrected
composition carries a `n(He 2^3S)` up to 99 per cent away from the
equilibrium one.

| quantity | selection | certified | refused | difference |
|---|---|---|---|---|
| red-pair equivalent width [%A] | `adv` | 1.374915 | 1.356726 | -1.32 % |
| red-pair depth [%] | `adv` | 4.957660 | 4.902392 | -1.11 % |
| blue-pair depth [%] | `adv` | 0.7137931 | 0.7057957 | -1.12 % |
| FWHM [A] | `adv` | 0.2603135 | 0.2598090 | -0.19 % |
| red-pair equivalent width [%A] | `solution` | 17.040779 | 16.916147 | -0.73 % |
| red-pair depth [%] | `solution` | 63.89678 | 63.54451 | -0.55 % |
| blue-pair depth [%] | `solution` | 14.51181 | 14.36944 | -0.98 % |
| FWHM [A] | `solution` | 0.2517373 | 0.2512328 | -0.20 % |

`log10 Mdot [g/s] = 7.87` for both, as the binary prints it. Taken from the
states themselves the wind-window mean of `rho v r^2` differs by +0.055 per
cent, which is 2.4E-04 in the logarithm.

Figure: `docs/figures/du_stop_vs_stationary_he10830.pdf`, the two line
profiles overlaid under each selection with their difference below.

## The answer

A state of this case that satisfies the flux-spread criterion of the
literature and is refused by this code's certification differs from the
certified stationary solution by at most 1.9 per cent in any profile, and
that difference DOES reach a measured quantity: the He I 10830 red-pair
equivalent width moves by 1.3 per cent and its depth by 1.1 per cent, while
the mass-loss rate does not move at the two decimals the code reports.

Where it sits is the point. The flux-spread criterion is a statement about
`rho v r^2`, and `rho v r^2` is the one quantity that agrees between the two
states to 0.055 per cent everywhere. What differs is the composition and the
temperature of the emitting layer between 1.7 and 17 R_p, driven by the one
row the certification refuses, the He/H elemental transport in the wind. The
criterion is therefore blind by construction to the physics that sets the
line: a wind can be flux-converged to eight decades and still carry a helium
partition that is wrong by a per cent where the line forms.

Against the observational side this is small. It is not small against the
differences the catalog resolves between neighboring models: the He/H rungs
of this group stand 0.178 %A apart in the equivalent width here (READ,
`docs/lhs1140b_catalog_state_20260920.tex`, section 4: 1.1970 at He/H = 1.70
against 1.3749 at 2.13), and the shift measured above is 0.018 %A, a tenth of
that gap, so a crossing He/H read off such a ladder would move by about a
tenth of a rung. The certification buys a resolution the flux criterion does not.

## What was not measured

The state compared here was produced by the stationary route stopped at one
outer pass, not by a marching run stopped at `du < 1e-3`, because no marching
run of this configuration reaches that threshold. A state reached by a
different route, in particular one whose hydrodynamic rows are far from
stationary while its flux happens to be flat (the catalog's
`scalar_x0.01_kzz1e9/HeH9.7`, `du = 4.4E-12` with `||R|| = 1.0`, READ), is
not represented by this comparison and would be further away.

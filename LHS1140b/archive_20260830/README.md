# LHS 1140 b: the runs made before 2026-09-13

This directory holds every EXHALE run of LHS 1140 b made between 2026-08-23
and 2026-08-30, moved here unchanged on 2026-09-13. It was `LHS1140b/exhale/`
and `LHS1140b/heh_series/` until then; nothing inside was recomputed, no
output file was rewritten, and no number was edited. The move renamed
directories and nothing else.

**Why it is an archive.** The 939 solutions below are not results of the
current code. The changes that landed after 2026-08-30 reach the atomic He/H
wind directly, and `MODELS.md` section 5 lists them: the photon-grid
quadrature (a 9.5 percent overcount of the H I photoionization rate removed),
the cell width `dr_j`, the Jupiter radius (7.1492e9 cm, +2.26 percent in
`R_0`), the He mass, the CODATA constants, the spectrum below 13.6 eV that
the He(2^3S) population depends on, the competing absorbers of the He
recombination photons, the characteristic lower boundary, and the ghost-row
fix in the readers of the transit synthesis. The radius change also moves the
grid, so the current binary cannot reload any state stored here without
`src/utils/map_state_to_grid.py`. The tree of record on the current code is
`LHS1140b/models/`, described in `LHS1140b/MODELS.md`.

**The scripts in here are a record, not a way to run anything.** The absolute
paths in `run_*.sh`, `finish_case.sh`, `resolve_wind.sh`, `run_closure.sh`,
`continue_*.sh`, `ladder_*.sh` and in every `closure.json` are the paths of
2026-08-30: the binaries they name have been rebuilt, the conda environments
several of them name are gone, and the directory they sit in has moved. They
are kept so that each result says how it was made. Reading tools that take
their own directory as the reference point still work in place
(`exhale/kzz_scan_table.py`, `exhale/heh_diff_scan_table.py`,
`heh_series/audit_checks.py`).

**The transit curves here carry the ghost rows.** Every `tpm_*` file and
every figure drawn from one was synthesized before the 2026-09-03 reader fix,
over a chord grid two cells too tall. Line-center depths move by -3.81
percent (the largest single move measured) down to less than 0.2 percent on
the weak lines, and 4 A band depths by -2.5 to +0.5 percent, almost always
shallower. `STALE_P48.md` in this directory is the full account. Nothing was
regenerated.

## Counts

| | runs (`input.inp`) | with an `output/` |
|---|---|---|
| `exhale/` | 933 | 932 |
| `heh_series/` | 6 | 6 |
| **this directory** | **939** | **938** |

Counted with `find <dir> -name input.inp`. The one run without an `output/`
is `exhale/basemetals_gm25/f_closure11p1`. The whole `LHS1140b/` tree holds
942 `input.inp` files: these 939 plus the three solutions of record in
`LHS1140b/examples/`.

No run in this directory turns molecular chemistry on: `Molecular chemistry`
appears in none of the 939 input files, so the `molecular_*` groups of
`MODELS.md` section 3 have no counterpart here.

## What the suffixes mean

A directory name ending in a suffix says which generation of the code and
of the chemistry produced it.

| suffix | meaning |
|---|---|
| none | the original generation, 2026-08-23 to 2026-08-28 |
| `gm25` | after the He(2^3S) + H Penning coefficient was replaced by the continuous Garcia Munoz (2025) form, `docs/Update_EXHALE_stage1.md` section 87, together with the advection post-process fix of section 88. Binary of 2026-08-29 04:56 |
| `pc090` | the same binary, with the chemistry from the corrected Photochem 0.9.0 build of section 90 (Equilibrate element-relative mass closure, bounded/scaled clima solves) |
| `j96` | after the Shull & van Steenberg (1985) secondary-ionization branching was renormalized onto the cell's own neutral `n(He I)/n(H I)`, section 96. Binary of 2026-08-30 05:02; `refresh_j96` uses the 2026-08-30 10:47 binary, carrying sections through 105 |

"Lower boundary" in the tables below is one of: **scalar** (T = 226 K,
`R_0` = 0.157692 R_J, p = 1 microbar, metal-free), **scalar + CNO** (the same
scalar base with the C, N, O reservoirs of the photochemical column handed
over through `base.inp`), or **profile** (the Photochem column as
`Lower atmosphere profile:`, with `T(p)`, `K_zz(p)` and the element
reservoirs at `p_match` = 1 microbar).

## The case directories directly under `exhale/`

These are the 2026-08-23 to 2026-08-27 runs, one directory per case, before
the work was organized into named series. `attempt_*/` subdirectories inside
them are discarded first attempts, not results.

| directories | what it was for | spectrum | lower boundary | diffusion and `K_zz` | He/H | runs | generation | `MODELS.md` group |
|---|---|---|---|---|---|---|---|---|
| `solar`, `heh0p25`, `heh0p5`, `heh0p55`, `heh0p6`, `heh0p7`, `heh1`, `heh10`, `heh100`, `heh1000`, `heh10000` | the first He/H scan of this planet: the composition axis of the He I 10830 line, and the He-rich audit of Phase C (`exhale/audit_summary.md`) | GJ 1132 | scalar | off | 0.0833 to 10000, 11 values | 22 | original | `atomic_scalar_gj1132_wellmixed` |
| `solar_gj699`, `heh0p04_gj699`, `heh0p06_gj699`, `heh0p25_gj699`, `heh0p55_gj699`, `heh1_gj699`, `heh1000_gj699` | the same scan on the alternate spectrum, and the two bracket cases that put the equivalent-width crossing on it without extrapolating | GJ 699 | scalar | off | 0.04 to 1000, 7 values | 14 | original | `atomic_scalar_gj699_wellmixed` |
| `heh0p55_diff`, `heh1000_diff`, `heh0p06_gj699_diff`, `heh0p55_diff_ctrl`, `heh0p55_diff_d3`, `heh0p55_diff_fix` | the first runs of the binary H/He element-diffusion operator, at its default `He_Kzz` = 0: does the separation move the wind or the line | GJ 1132 and GJ 699 | scalar | on, `K_zz` = 0 | 0.06, 0.55, 1000 | 6 | original | none (superseded by the `kzz*` groups) |
| `heh0p55_diff_kzz1e6` ... `heh0p55_diff_kzz1e11` | the eddy-coefficient sensitivity scan at one fixed composition, the material of `kzz_decision.md` section 3 behind the adopted `He_Kzz` = 1.0e9 | GJ 1132 | scalar | on, 1e6 to 1e11 | 0.55 | 10 | original | the `He/H = 0.55` case of each `atomic_scalar_gj1132_kzz*` group |
| `heh5_diff_kzz0`, `heh5_diff_kzz1e1` ... `heh5_diff_kzz1e6`, `heh5_diff_kzz1e9` | the low decades probed at a fixed reservoir instead of rescanned in composition, to find where the line stops distinguishing `K_zz` from zero | GJ 1132 | scalar | on, 0 to 1e9 | 5 | 8 | original | none (a fixed-composition probe) |
| the remaining `heh*_diff_kzz*` (39 directories, `kzz0` through `kzz1e11`) | the composition brackets that put the equivalent-width crossing at each decade of `K_zz`, `kzz_decision.md` sections 6 and 6.1 | GJ 1132 | scalar | on, 0 to 1e11 | 1 to 10, 26 values | 39 | original | `atomic_scalar_gj1132_kzz{0,1e5,...,1e11}` |
| `heh1000_fcheck`, `heh1000_ptc`, `solar_ptc` | solver and trapping checks, not science: a JFNK finish under `-fcheck=bounds,do,mem`, and two pseudo-transient starts | GJ 1132 | scalar | off | 0.0833, 1000 | 3 | original | none |
| `xuv0p01_heh2p13` ... `xuv0p33_heh2p13`, `xuv0p10_heh9p7` ... `xuv0p33_heh9p7` | the XUV grid against the 2025 non-detection, from the scalar-base model that reproduces the 2024 equivalent width | GJ 1132 scaled to 0.01-0.33 | scalar | on, 1e9 | 2.13, 9.7 | 12 | original | `atomic_scalar_gj1132x*_kzz1e9` |
| `xuv0p01_closure9p7` ... `xuv0p33_closure9p7` | the same grid from the flux-closed model, so the inference does not rest on one lower boundary | GJ 1132 scaled to 0.01-0.33 | profile | on, from the profile | 9.7 | 6 | original | `atomic_photochem_gj1132x*_kzzprofile` |

## The named series under `exhale/`

| directory | what it was for | spectrum | lower boundary | diffusion and `K_zz` | He/H | runs | generation | `MODELS.md` group |
|---|---|---|---|---|---|---|---|---|
| `flux_closure` | the first elemental-flux closure: the reservoir is solved for rather than prescribed, one rung per reservoir, each seeded from the rung below (Phase E, `docs/phase_e_flux_closure_design.md`) | GJ 1132 | profile | on, from the profile (1e9) | 2.09 to 12, 9 rungs | 54 | original | `atomic_photochem_gj1132_kzzprofile` |
| `kzz_profile_scan` | the `K_zz` scan of the crossing repeated under the profile boundary condition, where the scalar-base scan of `kzz_decision.md` section 6.1 had been made on scalars | GJ 1132 | profile | on, `K_zz` set per decade, 1e3 to 1e9 plus 5.6e6, 1.33e7, 1.78e7, 3.16e7 | 1.6 to 18, 33 values | 162 | original | `atomic_photochem_gj1132_kzzprofile` |
| `branch_diagnosis` | are the two converged states the `K_zz` = 1e7 scan reports at one reservoir two physical states? Verdict: the hotter one is a numerical artifact, a band of cells emptied of C, N and O between two shells at the handoff ratio, and refilling them collapses it onto the cooler branch | GJ 1132 | profile | on, from the profile | 3.399 to 14.643, 6 values | 14 | original | none (a diagnosis) |
| `metastable_step_diagnosis` | where the steps in the He 2^3S profile come from; they sit at the 4000 K seam of the retired Taylor (2025) two-branch Penning fit and at the outer boundary | GJ 1132 | scalar and profile | both | reads the stored runs | 0 | original | none (analysis scripts only) |
| `bump_gm25` | the six solutions in which the Schulik & Owen (2025) adiabatic-cooling bump of the He 2^3S population is measured, put on one vintage: the Garcia Munoz (2025) Penning coefficient and the current advection post-process (`LHS1140b/bump_analysis.py` reads them) | GJ 1132 and GJ 699 | scalar | off | 0.06, 0.0833, 0.55, 1000 | 6 | `gm25` | `atomic_scalar_gj{1132,699}_wellmixed` |
| `thermostat_gm25` | the three solutions behind the memo's thermostat table, re-solved: what the elemental C, N, O of the profile do to the wind temperature at fixed composition | GJ 1132 | scalar and profile | on, 1e9 and from the profile | 2.09, 2.13, 10.3 | 4 | `gm25` | `atomic_scalar_gj1132_kzz1e9`, `atomic_photochem_gj1132_kzzprofile` |
| `basemetals_gm25` | the one-key-at-a-time control ladder: He/H, base temperature, base radius, the CNO reservoirs, then the full profile, each changed alone, so the difference between the scalar base and the profile is attributed key by key | GJ 1132 | scalar, scalar + CNO, profile | on, 1e9 and from the profile | 2.13, 11.1 | 7 (6 solved) | `gm25` | `atomic_scalar_gj1132_kzz1e9`, `atomic_scalarCNO_gj1132_kzz1e9`, `atomic_photochem_gj1132_kzzprofile` |
| `xuvkzz_gm25` | three blocks re-solved on the new Penning coefficient: the XUV grid from both models, the pressure form of `K_zz`, and the outer region a converged wind inherits from its seed | GJ 1132, and scaled to 0.01-0.33 | scalar and profile | on, 1e9 and from the profile | 2.09 to 12, 10 values | 31 | `gm25` | `atomic_scalar_gj1132x*_kzz1e9`, `atomic_photochem_gj1132x*_kzzprofile` |
| `misc_gm25` | the six remaining places in the memo still standing on the retired Penning coefficient, re-solved: the GJ 699 seven-row scan, four `K_zz` decades, and four closure rungs | GJ 1132 and GJ 699 | scalar and profile | off and on, 1e8 to 1e10 and from the profile | 0.04 to 1000, 9 values | 17 | `gm25` | several (one per site) |
| `crossings_gm25` | the three equivalent-width crossings re-measured on the new Penning coefficient: well mixed, `K_zz` = 1e9 over a scalar base, and the flux-closed reservoir | GJ 1132 | scalar and profile | off and on, 1e9 and from the profile | 0.40 to 11.5, 31 values | 88 | `gm25` | `atomic_scalar_gj1132_wellmixed`, `atomic_scalar_gj1132_kzz1e9`, `atomic_photochem_gj1132_kzzprofile` |
| `ladder_gm25` | the `K_zz` ladder across every decade, the GJ 699 crossing, and the reservoir band the profile adapter had been refusing, all re-measured on the same binary | GJ 1132 and GJ 699 | scalar and profile | off and on, 0 to 1e11 and from the profile | 0.042 to 8.75, 43 values | 85 | `gm25` | `atomic_scalar_gj1132_kzz*`, `atomic_scalar_gj699_wellmixed`, `atomic_photochem_gj1132_kzzprofile` |
| `confirm_gm25` | every one of the five crossings had been read between two solved rungs and never solved at its own value; this directory solves each of the five there and compares | GJ 1132 and GJ 699 | scalar and profile | off and on, 0, 1e9 and from the profile | 0.0479, 0.4245, 1.6261, 3.8103, 9.048 | 11 | `gm25` | all five crossing groups |
| `nh_refusal_diagnosis` | why the profile adapter refused reservoir He/H in 8.5-8.7; a 40-point adapter scan from He/H = 1 to 15 with the conservation check disarmed, plus the same scan on the corrected Photochem build in `scan_pc090/`. Chemistry and climate only, no wind step | GJ 1132 | profile (the column itself is what is measured) | `K_zz` = 1e9 in the column | 1 to 15 | 0 | `gm25` / `pc090` | none (a diagnosis) |
| `clima_bracket_diagnosis` | why clima's background-pressure bracket failed at He/H = 9.4 and 9.5: the outer solve hands it a surface temperature of 6.5e8 K, past the 6000 K ceiling of the thermodynamic polynomials, because its forward-difference Jacobian is noise at MINPACK's default step. Chemistry and climate only | GJ 1132 | profile (the column itself) | `K_zz` = 1e9 in the column | 9.4, 9.5 | 0 | `pc090` and the 0.8.4 reference build | none (a diagnosis) |
| `deep_temperature_sensitivity` | what a 2.1e-6 K difference in the deep boundary temperature does to the handoff column, the observation left open by `crossings_pc090/results.txt` section 7. Thirteen adapter runs under `runs/`, chemistry and climate only | GJ 1132 | profile (the column itself) | `K_zz` = 1e9 in the column | 9.0 | 0 | `pc090` | none (a diagnosis) |
| `crossings_pc090` | the elemental-flux-closure crossing re-measured on the corrected Photochem 0.9.0 build, the same EXHALE binary as `crossings_gm25` | GJ 1132 | profile | on, from the profile | 8.45 to 11, 11 values | 36 | `pc090` | `atomic_photochem_gj1132_kzzprofile` |
| `clima_epsfcn_installed` | the clima finite-difference-step repair measured on the installed environment rather than on a replication of the outer solve; two closure rungs and their wind re-solves | GJ 1132 | profile | on, from the profile | 9.0, 9.1 | 6 | `pc090` plus the repair | `atomic_photochem_gj1132_kzzprofile` |
| `crossings_j96` | the five crossings re-measured after the secondary-ionization branching was renormalized on the cell's own composition; the largest single series here | GJ 1132 and GJ 699 | scalar and profile | off and on, 0 to 1e11 and from the profile | 0.042 to 12, 89 values | 181 | `j96` | all five crossing groups |
| `adapter_grid_fix` | one closure rung re-run with the adapter pinning Photochem's grid at the stated model top, everything else held: what the grid repair does to the He I 10830 equivalent width | GJ 1132 | profile | on, from the profile | 8.2117 | 4 | `j96` plus the adapter repair | `atomic_photochem_gj1132_kzzprofile` |
| `kzz_power` | is the pressure form of `K_zz` carried in the answer? `K_zz(p)` = 1e9 (p / 1.1e-8 bar)^-1/2 against the constant 1e9, in the chemistry alone and in chemistry plus wind | GJ 1132 | profile | on, `K_zz(p)` against constant 1e9 | 2.09, 11.1 | 6 | `j96` | `atomic_photochem_gj1132_kzzprofile` |
| `refresh_j96` | the whole memo re-measured on one binary after section 96: 82 JFNK re-solves in seven series, three closures re-run end to end, one HD 209458 b control. Subdirectories `basemetals`, `broadening`, `bump`, `closure`, `fig4style`, `misc`, `scan`, `thermostat`, `validity`, `xuvkzz`, each with its own `results.txt` | GJ 1132, GJ 699, and scaled to 0.01-0.33 | scalar, scalar + CNO, profile | off and on, 1e8 to 1e10 and from the profile | 0.04 to 1000, 23 values | 101 | `j96` (binary of 2026-08-30 10:47) | every group of `MODELS.md` section 3 except the `molecular_*` ones |

## `heh_series/`

`heh_series/` is not LHS 1140 b. It is the first attempt at the Phase C
helium-rich audit, run on the tutorial planet (1.38 R_J, 0.69 M_J,
`T_eq` = 1320 K, a = 0.047 AU) with a power-law spectrum, with the He/H
number ratio as the only parameter changed: `heh_1`, `heh_10`, `heh_100`,
`heh_1000`, plus two `heh_100` variants that changed how the marching was
handed to the steady solver. `heh_1`, `heh_10` and `heh_1000` completed;
`heh_100` did not reach a usable state in any of the three variants tried,
and Phase C was rescoped on 2026-08-24 onto the converged LHS 1140 b
solutions in `exhale/` instead. `audit_checks.py` is the audit itself and
takes case directories as arguments, so it is the tool that produced
`exhale/audit_summary.md`; with no arguments it audits the four cases here.
Six runs, all with an `output/`.

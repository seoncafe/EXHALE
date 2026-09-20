# LHS 1140 b: the model catalog

Every EXHALE model of LHS 1140 b, named by its physics, with what each one
is for and where its result is. This file is the record of the 2026-09-13
reorganization; the results table at the end is filled by
`models/status.py` as the runs complete.

## 1. Layout

```
LHS1140b/
  MODELS.md                 this catalog
  models/                   the model tree of record, solved on the current code
    <group>/HeH<value>/     one case: input.inp (+ base.inp or
                            lower_atmosphere_profile.dat), output/, tpm_*.txt,
                            REPRODUCE.md (how that case was reached),
                            state_index.json and states/<generation_id>/
                            (which state a reader gets, section 9)
    <group>/HeH<value>/tpm_turb/
                            where it exists, the same solution synthesized a
                            second time with EXHALE_TRANSIT_TURB=1, so that a
                            figure can draw the line with and without the
                            turbulence term on one wind (the well-mixed rungs
                            up to He/H = 1 carry it)
    diffusion_check/        the post-processing pass of four solved K_zz cases
                            repeated with EXHALE_DIFFUSION_CHECK=1, which is
                            the only way output/element_flux_profile.txt is
                            written for a prescribed-composition case; its
                            README.md states what was measured
    make_models.py          writes the tree from the group table below
    run_case.sh             one prescribed-composition case, start to line
    run_closure.sh          one flux-closure rung (Photochem + EXHALE)
    pick_seed.py            the solved state a case is seeded from
    publish_state.py        the one publisher of a state generation, and the
                            resolver every reader takes a state through
                            (section 9)
    import_legacy_states.py enters the states this tree already held into
                            that contract, with the provenance they have
    legacy_state_map.txt    a state of another directory that belongs to a
                            case, one line each, with the reason
    current_grid_Hydro_ioniz.txt
                            the cell centers every case is solved on, the
                            target grid of src/utils/map_state_to_grid.py
    write_reproduce.py      writes a case's REPRODUCE.md, called by the runners
    compare_archive.py      one case overlaid on its archived solution
    status.py               reads every case and writes sections 7 and 8
    README.md               how to reproduce the whole tree, and what it needs
  archive_20260830/         every run made before 2026-09-13 (formerly
                            exhale/ and heh_series/), untouched; its own
                            README.md maps each old directory to a group here
  examples/                 the three solutions of record of 2026-08-30
                            (superseded by models/; kept until the new
                            solutions of record are chosen)
  sed/, pwinds_oracle/, pwinds_refit/, lower_profile/, Cherubim_2026/
                            inputs and reference material, unchanged
```

## 2. Naming

A group name states the physics that does not change along its composition
ladder, in a fixed order:

```
<chemistry>_<lower boundary>_<spectrum>_<mixing>
```

| field | values | meaning |
|---|---|---|
| chemistry | `atomic` | H, He, He(2^3S), electrons; no molecules (`Molecular chemistry` absent) |
| | `molecular` | `Molecular chemistry: True` and `Molecular carrier transport: True`: H2, H2+, H3+, HeH+ with the H2 carrier transported, the molecular base from `base.inp` (`q_H2_base`, `p_base`) or from the profile |
| lower boundary | `scalar` | prescribed scalar base: T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free |
| | `scalarCNO` | the same scalar base plus the C, N, O reservoirs of the photochemical column at the matching level, through `base.inp` (`C_H_base`, `N_H_base`, `O_H_base`) |
| | `photochem` | the Photochem column as `Lower atmosphere profile:` (T(p), K_zz(p), q_H2, the element reservoirs at p_match = 1 microbar); the reservoir He/H is the rung, and for `atomic` groups the rung is solved by the elemental-flux closure (`src/utils/element_flux_closure.py`) so the composition is an output |
| spectrum | `gj1132` | the GJ 1132 Mega-MUSCLES v23 proxy at the b orbit, the paper's fiducial (`sed/lhs1140_sed_gj1132_at_b.txt`) |
| | `gj699` | the GJ 699 proxy, the paper's alternate (`sed/lhs1140_sed_gj699_at_b.txt`) |
| | `gj1132x0.30` | the GJ 1132 proxy with its XUV scaled to 0.30 of the fiducial (`sed/lhs1140_sed_gj1132_at_b_xuv0p30.txt`); values 0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.33 |
| mixing | `wellmixed` | no element diffusion (`He_diffusion` absent): He/H uniform |
| | `kzz1e9` | binary H/He element diffusion with `He_Kzz: 1.0e9` cm^2/s (the adopted value, `kzz_decision.md` section 0); `kzz0`, `kzz1e5` ... `kzz1e11` likewise |
| | `kzzprofile` | diffusion with K_zz(p) taken from the profile (the profile's constant 1e9 today) |

The case name is the reservoir He/H by number: `HeH0.42`, `HeH2.13`,
`HeH1000`. For a closure rung it is the STARTING reservoir; the converged
value is in the rung's last `kNN/` iterate.

Every case carries: `input.inp` exactly as solved, `output/` (the solved
state, the `_adv` profiles, the heating and cooling breakdowns,
`element_flux_profile.txt` where the diffusion operator writes it),
`EXHALE_setup.out`, `EXHALE_resolved.out`, `tpm_*.txt` (He I 10830 at
R = 68 000 through `winered_hires_y.sh`, H-alpha, H-beta, Ly-alpha) and
`run.log`, `pp.log`, `transit.log`.

## 3. The groups and their ladders

Every ladder brackets the measured red-pair equivalent width,
1.108 +/- 0.030 %A (Cherubim et al. 2026), on the crossing the 2026-08-30
solutions found (`WORKPLAN.md`, superseded-numbers table), widened by one
step on each side because the code has moved since (section 5).

| group | He/H cases | purpose |
|---|---|---|
| `atomic_scalar_gj1132_wellmixed` | 0.083, 0.40, 0.42, 0.44, 0.55, 1, 10, 100, 1000 | the well-mixed crossing (was 0.4132) and the He-rich audit ladder to He/H = 1000 |
| `atomic_scalar_gj699_wellmixed` | 0.042, 0.046, 0.050, 0.083, 1, 1000 | the same on the alternate spectrum (crossing was 0.0479) |
| `atomic_scalar_gj1132_kzz0` | 0.55, 2.6, 3.0, 3.5, 3.7, 3.9 | crossing at K_zz = 0 (was 3.81; the 2026-09-14 solutions put it below 3.5, so 2.6 and 3.0 were added) |
| `atomic_scalar_gj1132_kzz1e5` | 2.6, 3.0, 3.35, 3.64, 3.93 | crossing at 1e5 (2.6 and 3.0 added 2026-09-14, the crossing lying below 3.35) |
| `atomic_scalar_gj1132_kzz1e6` | 0.55, 2.4, 2.8, 3.19, 3.46, 3.74 | crossing at 1e6 (2.4 and 2.8 added 2026-09-14, the crossing lying below 3.19) |
| `atomic_scalar_gj1132_kzz1e7` | 0.55, 2.70, 2.94, 3.18 | crossing at 1e7 |
| `atomic_scalar_gj1132_kzz1e8` | 0.55, 2.05, 2.23, 2.41 | crossing at 1e8 |
| `atomic_scalar_gj1132_kzz1e9` | 0.55, 1.50, 1.60, 1.70, 2.13, 4.0, 9.7 | the adopted K_zz: crossing (was 1.6108), the 2.13 reference of the metal ladder, the 9.7 companion of the closure rungs |
| `atomic_scalar_gj1132_kzz1e10` | 0.55, 1.06, 1.15, 1.29 | crossing at 1e10 |
| `atomic_scalar_gj1132_kzz1e11` | 0.55, 0.795, 0.865, 0.93 | crossing at 1e11 |
| `atomic_scalarCNO_gj1132_kzz1e9` | 2.13 | the metal effect at fixed composition (the C -> D rung of the memo's ladder) |
| `atomic_photochem_gj1132_kzzprofile` | 2.09, 3, 5, 7, 8, 9, 9.7, 10, 12 | the flux-closure reservoir ladder (crossing was 8.212); the 9.7 rung is the column the XUV photochem grid holds fixed |
| `atomic_scalar_gj1132x{0.01,0.10,0.15,0.20,0.25,0.30,0.33}_kzz1e9` | 2.13, 9.7 | the XUV grid against the 2025 non-detection, scalar base |
| `atomic_photochem_gj1132x{0.01,0.10,0.15,0.20,0.25,0.30,0.33}_kzzprofile` | 9.7 | the XUV grid above the photochemical column: the column of the FIDUCIAL He/H = 9.7 closure rung held fixed (its converged iterate's `lower_atmosphere_profile.dat`), the wind re-solved on the scaled spectrum, no closure (as the 2026-08-30 grid did: a Photochem climate solve on a spectrum scaled as a whole freezes the deep atmosphere, 136 K at 17 bar, and its chemistry fails elemental closure) |
| `molecular_scalar_gj1132_wellmixed` | 0.083, 0.55, 2.13 | NEW: the molecular layer solved in EXHALE above the scalar base (`base.inp` with `q_H2_base` from the photochemical column's H2 fraction, `p_base` 1 microbar) **None of the three is solved** (all three re-tried 2026-09-18, item L34 step (c)). **A well-mixed molecular case can be entered only through the atomic-to-molecular conversion of its own certified atomic pair**: the run rescales the whole column onto its own He/H, and that rescale is refused outright when the state carries HeH+, which holds one nucleus of each element (`load_IC.f90` lines 650 to 663). Neither the element-treatment step from a certified `kzz1e9` state nor a composition step inside this group is therefore admissible, and 2.13, which has no atomic rung at its composition, was given one solved for the purpose in `models/.L34/atomic_wm2.13`. From those seeds the wind converges (`hydro info = 0` from pass 3 to 5) and the H2 carrier balance does not: 0.55 reaches 3.76e-02 at pass 11 of 31 and rises again as its front walks outward, 0.083 reaches 4.13e-02 over forty passes and rises over forty more, 2.13 flattens at 4.3e-01 after seven. Every relaxation of every pass ends on the composition movement bound and none on its own residual (`docs/lhs1140b_stationary_L34c_20260918.md`). Sections 7 and 8 are written by `status.py` from the case directories and carry the verdict alone, so the reason is here and in each case's own `not_solved.md`. |
| `molecular_scalar_gj1132_kzz1e9` | 0.083, 0.55, 2.13, 9.7 | NEW: the same with element diffusion |
| `molecular_photochem_gj1132_kzzprofile` | 2.09, 9 | NEW: the molecular layer above the converged column of ITS OWN atomic closure rung (`atomic_photochem_gj1132_kzzprofile/HeH2.09/k<last>/lower_atmosphere_profile.dat` and `HeH9/k<last>/...`), held fixed, no closure. **Why the rung's column and not a stored one** (changed 2026-09-16): a molecular case of this group is solved from the certified wind of its atomic pair, and `load_IC` admits a seed only when the two states agree on the column and on the reservoir to a part in 1e6. The 2026-08-30 stored columns carry He/H = 2.092388551 and 9.010268395 at the matching level only by coincidence of their own history; measured on 2026-09-16, the stored 2.09 column and the re-run 2.09 rung differ by 1.2e-04 in He/H, twenty times the threshold, and the seed was refused. The case name is now the rung's name and the He/H line carries the value that column actually holds, so the pair cannot drift apart again. (The former `HeH9.05` directory named the stored column's value and has been removed; it held no result.) **Neither case is solved** (re-tried 2026-09-18 from the certified wind of its own atomic pair, item L34 step (c)): a molecular SCALAR state cannot stand in for that seed, because this group's cases carry C/H, N/H and O/H reservoirs at the matching level and a metal-free scalar state carries none, which `load_IC` refuses on the carbon ratio. From the conversion seed, 9 is refused at outer pass 1 with its hydrodynamic rows at 2.8e-01 and 3.2e-01 (the element relaxation finds no admissible advance and restores its entry composition, which ends the outer loop by construction) and 2.09 bands between 1.9e-01 and 4.4e-01 over twelve passes (`docs/lhs1140b_stationary_L34c_20260918.md`). |

The `He/H = 0.55` case of each `kzz*` group is the K_zz sensitivity scan
of `kzz_decision.md` section 3 at one composition.

**How long a case of each group is given** (item D9 of
`docs/PLAN_20260918_rev2.md`, 2026-09-18). `models/run_campaign.sh` gives a
case a pass ceiling and a wall-clock ceiling, and both are EMPIRICAL
ALLOWANCES read off recorded solves of the same configuration class, not a
rule derived from a residual or a front speed. The `atomic_*` groups are
`atomic_prescribed`, 40 passes and 30 m; a rung of the elemental-flux closure
is `atomic_closure_rung`, 20 passes for each EXHALE solve of the rung and 45 m
for the series; the `molecular_*` groups are `molecular_alternation`, 40
passes and 6 h; a case that turns `Ionization transport` on is
`transported_ionization`, 90 passes and 6 h. Where each number comes from is
written beside the table in `run_campaign.sh`, and the one number of it that
is not a measurement of its own class says so: no wall clock was recorded for
a transported-ionization solve, so its six hours are the limit carried over
from the molecular runs. `models/budget_overrides.txt` gives one case another
budget and must name the authorization in the same line;
`EXHALE_OUTER_PASSES` or `CAMPAIGN_WALL` in the environment replaces the table
for a whole run. **A case that reaches either ceiling is UNCERTIFIED and
INCOMPLETE and is never a solution**: `campaign_status.txt` records
`ceiling=passes` or `ceiling=wall` beside it, and the case's own
`not_solved.md` carries the residual profile, the binding row and cell, the
solver verdict and the seed.

## 4. The molecular base

### What the base level is

The base level is the level at which the lower-atmosphere model hands the
column over, and `p_base`, `T_base`, `q_H2_base` and the elemental ratios
are that model's statement ABOUT THE LEVEL, not only about gas that happens
to be moving upward through it. The lower boundary reads them that way: the
interior state and the reservoir are matched acoustically first, the
direction of the contact is read off the matched face velocity, and the
face then carries the reservoir's density, temperature and composition
wherever gas enters the domain AND at a pressure-balanced contact at rest.
Only a reverse flow carries the interior's own entropy and composition out
through the level, and the reservoir then states the single condition its
one entering characteristic allows, the pressure. The physical assumption
is that the lower atmosphere below the level is dense and radiatively
controlled, so on the time scales of a stationary solution it is a heat
bath and a column at rest above it takes the level's temperature; its
validity is a subsonic base at a level inside the radiatively controlled
lower atmosphere, which is what 1 microbar is for these planets.

Each enabled transport operator carries its own base condition beside that
one, and none of them follows the direction the contact is upwinded on:
thermal conduction holds the level's own temperature, the element diffusion
holds the reservoir's composition as a Dirichlet value and reports the
diffusive flux it drives, and the molecular carrier transport carries no
diffusive flux across the base face and upwinds its advective term on the
face mass flux.

States carry the boundary model they were solved under in their header
(`# boundary_model`). The model above is
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`; a state
written under `..._smoothstep_v1` differs in what the face carries at and
near zero flow, and its certificate is that model's and does not transfer.

### The H2 partition

`q_H2_base` is the H2 volume mixing ratio of the inflowing gas. The
photochemical column at 1 microbar (`lower_profile/lower_atmosphere_profile.dat`,
He/H = 2.092) has q_H2 = 0.19265 and q_H = 6.5e-6, i.e. 99.998 percent of
its hydrogen nuclei are in H2. The scalar molecular cases keep that
fraction f = 0.99998 at every He/H, so

    q_H2_base = (f/2) / ((1 - f) + f/2 + He/H)

which is 0.85761 at He/H = 0.083, 0.47618 at 0.55, 0.19011 at 2.13,
0.04902 at 9.7, each just below the ceiling 0.5/(0.5 + He/H) the reader
enforces. The base temperature and radius are the scalar family's
(226 K, 0.157692 R_J) so that the only change against the atomic scalar
cases is the chemistry.

The 2.0e-5 in relative terms that separates `q_H2_base` from that ceiling is
NOT what stops the cold march of a molecular case. MEASURED 2026-09-13 at
f = 0.9998, which leaves 2.1e-4 to 3.7e-4 of room instead:
`molecular_scalar_gj1132_kzz1e9/HeH2.13` stops the same way, on the same
chemical-equilibrium refusal, 197 steps later (section 6 and
`docs/lhs1140b_stationary_L4c_20260913.md` section 3). The column's own
value is kept.

## 5. Why everything is re-solved

The 939 solutions under `archive_20260830/` were made between 2026-08-23
and 2026-08-30. Since then the code has moved in ways that reach the atomic
He/H wind directly (`docs/code_status_20260910.md` section 5,
`docs/Update_EXHALE_stage2.md` section 6): the photon-grid quadrature (a
9.5 percent overcount of the H I photoionization rate removed), the cell
width `dr_j`, the Jupiter radius (7.1492e9 cm, +2.26 percent in R_0, which
also changes the grid so that no archived state can be reloaded without
`src/utils/map_state_to_grid.py`), the He mass, the CODATA constants, the
spectrum below 13.6 eV that the He(2^3S) population depends on, the He
recombination photons' competing absorbers, the characteristic lower
boundary, and the readers' ghost-row fix of the transit synthesis. None of
the archived numbers is a result of the current code.

### Why the 2026-09-15 re-run

The campaign was solved on 2026-09-14 and then solved again on 2026-09-15.
Plan item L21 (`docs/lhs1140b_stationary_L21_20260915.md`) found two defects
in the base boundary condition, both of them in how it chooses the state of
the ghost below the first cell:

1. the entropy branch was keyed on the CELL-1 VELOCITY. On LHS 1140 b that
   velocity is the odd-even artifact of item P44 and not the flow, so the
   boundary selected the interior isentrope and put it at the base instead of
   the reservoir `base.inp` states;
2. the handover window was 1e-6 in face Mach number, above this planet's own
   base Mach number of 4.6e-7, so the branch could not have been decided by
   the flow even with the velocity clean.

The branch is now keyed on the wind-window mass flux and the window is 1e-8.
MEASURED on the fiducial case, corrected against certified: the ghost the
boundary states moves from 413.7 K to 226.0 K, which is `base.inp`'s `T_base`
to the digit it is given in; the first cell from 417.6 K to 238.0 K and its
density by a factor 1.733; the profiles cross near cell 21 and agree to 0.2
per cent in the outer wind; Mdot moves from 7.5096e7 to 7.3767e7 g/s, -1.77
per cent, log10 7.8756 to 7.8679. The movement is far above the 0.1 per cent
below which a change is treated as identical, and the base layer is exactly
the layer the He I 10830 equivalent width is most sensitive to, so the line
products are part of the re-run and not separable from it. The other
regression bases are byte-identical under the correction.

The results of 2026-09-14 are kept, unchanged, in
`LHS1140b/models_20260914_preL21/` (83 case directories: the 74 prescribed
atomic cases that solved and the 9 flux-closure rungs, each with its own
`REPRODUCE.md`, and that tree's `README.md` states what they are). They are
superseded and are not to be quoted as results; they are the record behind
the memo `docs/lhs1140b_exhale_vs_pwinds.pdf` and the sections below as they
stood, and they are the SEED of the re-run: `pick_seed.py` takes a case's own
2026-09-14 state as its first tier (`tier0`), the equations above the first
cell being the ones that state already solves.

The re-run is on one binary throughout, `EXHALE.x` md5 `db87b88d1ce5` (built
2026-09-15 20:38), which carries L4h, L7e, L7f, L18, L19 and L21, on two
machines of the same NFS mount (`lart4`, `lart3`); each case's `REPRODUCE.md`
names the binary, its md5 and the host. 79 of the 83 cases certify.

**The catalog of record is the FINAL PASS on `EXHALE.x` md5 `c2e9c9990b9f`**
(built 2026-09-16 10:43, manifest
`models/BINARY_MANIFEST_c2e9c9990b9f.txt`), which carries the review fixes as
well: every one of those 79 was solved again on it, each from its own db87
state, and all 79 certified in one outer pass. The two generations agree far
inside the certification tolerances -- the red-pair depth to the last bit in
all 79 cases, the mass-loss rate to the last bit in 52 of 79 and to 8.7
decimal digits at worst, the base cell to 8.3 -- which is what the same
equations on two binaries have to give. `models/compare_trees.py` prints that
block. The db87 results are preserved in `LHS1140b/models_20260915_db87/`.

FOUR CASES DO NOT SOLVE and are not retried:
`atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` and the photochemical cases at
0.15, 0.20 and 0.25 of the fiducial spectrum, all held by the hydrodynamic
rows under a composition that is still travelling, at the raised pseudo-time
start as at the default (the memo `docs/lhs1140b_rerun_20260915.md` section 4
measures both). With the three at 0.01 of the fiducial XUV (item L17) they
are the seven cases sections 7 and 8 below show as unsolved. Each
case's `REPRODUCE.md` names it and its md5 again. The three cases at 0.01 of
the fiducial XUV were not attempted: they did not solve before (item L17,
diagnosed as the HLLC contact selection on a base subsonic to Mach 1e-8) and
nothing in L21 touches that. **All seven are solved with the Roe flux from
2026-09-17**, on the same grid and by the same recipe; the subsection at the
end of section 6 states that decision, its reason and the flux systematic,
and the `flux` column of section 7 marks the rows. The nine molecular cases are not part of this
re-run; they follow item L7g.

## 6. Run recipe: what was measured on 2026-09-13 before any case was run

The archived solutions were finished on the code of 2026-08-30 by a JFNK
solve from a marched snapshot (`archive_20260830/exhale/finish_case.sh`:
`EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0`), the marches
themselves having run for many hours (the `heh1` march of 2026-08-23 ran
from 20:08 to 10:52 the next day and was then terminated; its `du` never
reached the 1e-2 hand-off). On the current code, MEASURED on the scalar
K_zz = 1e9 crossing case (`examples/scalar_base_kzz1e9`, 8 threads):

| route | result |
|---|---|
| reload the archived state directly | refused: the grid construction changed on 2026-09-10 (centers differ by up to 8.4e-4), whatever the base radius (the grid in R_p is the same for every base radius) |
| map the archived state onto the current grid (`src/utils/map_state_to_grid.py`), then the JFNK solve (`EXHALE_PTC`, dtau0 = 1.0, 0.1, 0.01) | stagnates at the cellwise-maximum residual 1.6 to 2.0 for 56 iterations; the worst rows are the momentum and energy of cells 1 to 10 (r < 1.003) and the loaded state's residual is O(1) only below r = 1.03 (mass 1.9 at cell 4, momentum 0.6 at cell 1, energy 2.0 at cell 119) |
| the same for the well-mixed crossing case (`crossings_j96/wm_heh0p42`) | the same: 1.96 to 1.78 in 8 iterations, worst rows at r = 1.000 to 1.001 |
| map, then march with the Newton hand-off (`Solver: Newton`, `du_th 0.5 1.0e-3`) | the mass-flux spread `du` of the mapped state is 0.8 to 8 (the archived solutions themselves have a spread of 0.09 above 2 R_p, the wind being subsonic to the 30 R_p boundary); it never falls below the 1e-2 hand-off in 5600 steps; the staged secondary ionization switches on at the first hand-off attempt and re-relaxes the wind (`Secondary_ionization: Immediate` avoids that on a restart) |
| cold start, two-stage march | 0.4 to 0.6 s per step at 8 threads; `du` 9.45 at step 1 and 9.58 at step 954 |
| the regression fixture of the atomic element-diffusion Newton finish (`backup/regression/atomic_elem_newton`, HD 209458 b) | ends `info = 2`, residual 1.89, NOT CERTIFIED (11 entries): the current code does not certify an atomic wind with element diffusion on any planet, which is the state the stage-2 solver program was closed in (2026-09-10) |
| the no-diffusion Newton fixture (`wasp_full_newton`) | `info = 0`, residual 4e-9, CERTIFIED |
| the partitioned stationary route (`Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`) from the mapped archived state | outer pass 1: hydro `info = 2` (mass 1.47, energy 1.33), the composition relaxation finds no admissible advance, REFUSED; `info = 1` after 322 s |
| the same route from the cold march's snapshot at step 5000 (du 8.4) and from the mapped well-mixed march's snapshot at step 8000 (du 10.7) | the snapshots are not near stationarity anywhere (cellwise residual 1.4 to 2.0 at cells 135 to 397, r = 1.05 to 5.7); the hydro solve finds no descent direction after 18 and 145 iterations; NOT CERTIFIED (`info = 1`, `info = 2`) |
| the marches themselves | `du` does not settle: 9.45 (step 1) to 8.35 (step 5212) from the cold start; 1.32 (step 499) to 10.7 (step 8147) from the mapped well-mixed state |

The conclusion at that point was that a solution of record for LHS 1140 b on
the current code needed either a long march to a plateau followed by a JFNK
solve that returns an uncertified state (what the diffusion fixture does), or
a solver item on the base rows. `docs/PLAN_20260913_lhs_stationary.md` is that
item; its L5c closed it with no source change, and the recipe below is what it
left.

### Recipe adopted 2026-09-13

**Why the well-balanced option is in every input.** At the base of this
planet the flow's own Mach number is 3e-7, and the contact speed the HLLC
flux puts through a base face is then set by the WENO3 reconstruction's
third-order pressure truncation error rather than by the velocity -- one to
three decades above it -- so the energy row of every base cell measures its
own largest term against itself and no seed can lower it. `Well balanced:
True` (K46) builds the Riemann pressure jump out of the departure from each
cell's own hydrostatic equilibrium, which is the quantity that is small
there; the same seed, the same input and the same binary then converge
quadratically and certify (`docs/lhs1140b_stationary_L5c_20260913.md`).

**The four passes of `models/run_case.sh`.**

| pass | what it does |
|---|---|
| 0, the seed | `models/pick_seed.py` names the solved state of the same physics nearest in log He/H. COMPATIBILITY COMES FIRST, in every tier including the case's own directory (item D9, 2026-09-18): a candidate is put against the restart contract of `load_IC.f90` before any distance is taken, and what is compared is the spectrum and its normalization, `K_zz`, the boundary kind, the elemental reservoir INVENTORY (which elements the state's `# reservoir` line names, the ratios themselves being carried across by the mapper), the four option tokens that decide which species the state files carry (`metals`, `mol`, `oxychem`, `carrier`), `he_diff`, and `iontrans`. The candidate's configuration is read from the STATE, out of the restart metadata block the binary wrote into it, falling back to `EXHALE_resolved.out` beside it and then to the directory's `input.inp`; the case's own is read from its `input.inp`, which is what the next run asks for. A state written without `Ionization transport` may still seed a run with it, but only as an explicit MODEL-OPTION TRANSITION: the case must name `iontrans` on its own `Restart option change:` line, the seed line says so in as many words, the candidate ranks after every candidate of the case's own system whatever its distance, and no certificate of the source problem transfers to the target. `--why` names every candidate that was refused and why. The tiers: a CERTIFIED case of the same group first, then a certified case of another group with the same physics, then the archive, newest generation first. `src/utils/map_state_to_grid.py` interpolates it onto `models/current_grid_Hydro_ioniz.txt`, the cell centers of the current code, and writes it as the case's `*_IC.txt`; where the seed was solved at another He/H, `--reservoir He/H` carries its helium onto this case's before it is written, so that the base rows arrive at this composition and a diffused He/H profile keeps its shape. Skipped where the case states `Load IC? False`. Where nothing carries the case's XUV normalization (the 0.30-scaled closure rung is the one such point) the nearest scaling of the same star and the same spectral shape is taken instead, and the line says so. Which tier the seed came from is on the line, and in `REPRODUCE.md` |
| 1, the wind | the binary, with `EXHALE_PTC_DTAU0=1.0` and `OMP_NUM_THREADS` (8 by default). The verdict is `the stationary solve returned info = 0` and the last certification block; `info = 2` is a state written but not certified |
| 2, the profiles | the same binary on its own solved state through `Restart intent: stationary evaluate`, which evaluates the held state with NO time step and writes the `*_adv.txt` files, the heating and cooling breakdowns, the steady-state `Mdot` and the certification from that evaluation (item L9, 2026-09-17). The route this replaced, `Do only PP: True` with `CFL: 1.0e-12`, is what every case on disk today was post-processed with, and the block below records what was wrong with it |
| 3, the line | `EXHALE_transit.py` through the WINERED HIRES-Y kernel of the measurement (`LHS1140b/winered_hires_y.sh`) |

**What a case directory holds, and when a second solve is taken (2026-09-18,
item D8 of `docs/PLAN_20260918_rev2.md`).** A case directory publishes ONE
solve, in `output/`, and it is the LATEST COMPLETED SOLVE: the runner compares
no two solves and claims none is the best. `ENDING` beside it says in one line
how that solve ended, every solve the runner superseded is kept whole in
`solve_<n>/` with an `ENDING` of its own, every seed attempt in `attempt_<n>/`
likewise, and `REPRODUCE.md` lists them all. Before anything is decided the
runner CLASSIFIES the ending out of the log and the state on disk, as `solved`,
`composition_refusal` (every hydrodynamic row within its own tolerance, so what
refuses is the composition), `hydrodynamic_refusal`, `nonfinite_state`,
`state_missing`, `no_verdict` or `no_row_measures`, and records whether the
outer loop spent its pass budget. The pseudo-time continuation of item L4e is
taken for `hydrodynamic_refusal` ALONE, because that is the ending it
addresses; every other refusal ends the case with its reason in the case's own
`not_solved.md`, which is the `reason` column of section 7. Until 2026-09-18
the continuation was taken on any nonzero `info`, and on
`molecular_scalar_gj1132_wellmixed/HeH0.083` (item L34c) it took the gated
carrier row from 4.13e-02 back to 7.37e-02 and published that state over the
better one. The evaluate pass of pass 2 runs in `eval/`, with its own copy of
the case's inputs, so the case's `input.inp` is never opened for writing and is
always the file the solution was started from. `run_campaign.sh` records the
exit status of each case and its ending class in `campaign_status.txt`, adds
the budget that case was given and whether it reached a ceiling, and exits
nonzero when any case failed.

**The closure rungs follow the same rule** (item D9, 2026-09-18).
`src/utils/element_flux_closure.py` used to restart a refused solve at the
raised pseudo-time start whenever `info != 0` and a state was written, with no
classification at all, so a rung continued on a composition-only refusal just
as the runner did. It now SOURCES the policy block of `models/run_case.sh`
(`RUN_CASE_POLICY_ONLY=1` defines it and returns) and takes the continuation
for the one ending that block addresses, keeping the first solve whole in the
iteration's own `solve_<n>/` with its `ENDING`; the rule is written once and
read once. `RUN_CASE_POLICY` names the file where it is elsewhere, and a
checkout without it takes no continuation and says so.

The keys `make_models.py` writes for an atomic case, beyond the physics of
its group:

```
Reconstruction scheme: PLM        the route rebuilds the state under WENO3
Load IC? True
Solver: Newton                    "Restart intent: stationary" needs it (K43)
Restart intent: stationary
Secondary_ionization: Immediate
Well balanced: True
```

A `molecular_*` case has none of them but `Well balanced: True`: it keeps
`Load IC? False`, `du_th [PLM,WENO3]: 0.5 1.0e-3` and `Solver: Newton`,
because the archive holds no molecular state to start from. The closure
rungs carry the same four keys through `closure.json`'s `input_keys`, which
`src/utils/element_flux_closure.py` adds to each iteration's input, with
`exhale_env` reduced to `EXHALE_PTC_DTAU0` alone.

**Why the secondary-ionization coupling is applied from the first
evaluation.** Its STAGED default arms only when a march first converges,
and this route takes no time step, so a state solved under the default is
written with `sec_ion=F` in its own coupling header while the
post-processing pass -- which applies the coupling unconditionally and is
where every reported number comes from -- reports the same wind with it.
MEASURED on the `K_zz` = 1e9, He/H = 1.6261 case from the mapped seed, 16
threads:

| | STAGED (the default) | `Secondary_ionization: Immediate` |
|---|---|---|
| certified at | outer pass 8 | outer pass 7 |
| coupling header of the state | `sec_ion=F` | `sec_ion=T` |
| log10 `Mdot` [g/s] | 8.29 | **7.87** |
| archived 2026-08-30 value | 7.80 | 7.80 |

and between the two solutions `n(He 2^3S)` differs by 74 percent beyond
1.5 R_p, where the He I 10830 line forms. The danger the staging was
written for is a cold start (K14c); this route starts from a solved state.

**Why the post-processing pass ran at `CFL: 1.0e-12`, and why it no longer
does (superseded 2026-09-17 by item L9).** The measurement below stands as
taken; what it did not measure is the three things that route was doing
besides suppressing the step. It reconstructed with PLM where the solution was
reached with WENO3; it overwrote the solved state's own `Hydro_ioniz.txt` with
an uncertified relaxation snapshot (`certified=F
cert_reason=no_stationary_claim`), so **every
`LHS1140b/models/*/output/Hydro_ioniz.txt` in this tree is such a snapshot**
and the certified state survives only in the `_IC` copy; and it wrote a
certification-claiming `# coupling:` header onto the `_adv` files, which
describe a different composition. The measure the advection operator itself
weights by, `adv_mass_row`, is nine decades smaller at the base and four in
the wind on the evaluate route (1.5731e+00 against 1.7338e-09 at cell 1;
3.0840e-05 against 3.9957e-14 at cell 500), and the He I 10830 red-pair
equivalent width moves by 5.1e-06 relative, so the line numbers of this
catalog are not in question, only the states the products were made from.
Record: `docs/lhs1140b_stationary_L9_20260916.md`. The original measurement
follows. It takes one time step before it stops. MEASURED on the certified state above, the state the
pass writes differs from the state it reads by

| | default CFL | `CFL: 1.0e-12` |
|---|---|---|
| `rho` | 1.2e-6 | 3.6e-14 |
| `v` | 0.59 | 2.3e-12 |

(largest relative difference over the physical cells, PLM). What remains at
1e-12 is the composition, which the pass re-equilibrates under the coupling
above and not through any time step.

**The molecular groups still have no working recipe**, and the campaign
should not spend time on them until they do; `models/run_case.sh` therefore
refuses a case whose input states `Molecular chemistry: True` unless
`MOLECULAR=1` is set (the change is held as
`models/run_case.sh.molecular.patch` while the atomic campaign is reading
that file). MEASURED 2026-09-13 on
`molecular_scalar_gj1132_kzz1e9/HeH2.13` and
`molecular_scalar_gj1132_wellmixed/HeH0.083`, 8 threads
(`docs/lhs1140b_stationary_L4c_20260913.md`):

- **seeding an atomic solution is refused by the restart contract**, as K44
  says it must be: `mol` and `carrier` decide how many unknowns a state
  has, so `Restart option change: mol carrier` stops at startup and the
  state is not a restart of that run whatever is named. The same rule
  refuses `carrier_newton`, so a `Coupled carrier solve: True` state can
  only come from a run that marched under that key.
- **the cold march stops at step 1199 on the chemical equilibrium of the
  cells below the base face**, with `ERROR STOP ioniz_eq: persistent
  non-root chemical equilibrium` after 1000 consecutive non-root sweeps.
  The base H2 fraction standing 2.0e-5 below its ceiling is NOT the cause:
  at f = 0.9998, which leaves 2.1e-4 to 3.7e-4 of room (section 4), the same
  run stops the same way at step 1396, at ghost cell -1 (r = 0.999807),
  element violation 2.395e-4 against 1e-6. The violation of that cell
  crosses the tolerance at step 984 and grows monotonically from 6.5e-6
  to 2.4e-4 without a single root in between.
- **the well-mixed He/H = 0.083 case stops earlier and for another reason**:
  at step 980 the energy update finds no temperature for four cells, the
  first at 7.06e4 K, because a cell that `caloric_eos` still counts as
  molecular is held below the 5e4 K top of the H2 rovibrational table
  (`T_ceiling_mol_K`), and nine step halvings down to dt/2^8 change nothing.
  The run exits 2 and writes no state.
- **a snapshot taken below both stops is not near stationarity**: at step
  900 the mass-flux spread `du` is 9.58, and the partitioned stationary
  route from it ends `info = 1` in every variant tried -- carriers
  transported, carriers in the Newton unknown vector, or the molecular
  content in local equilibrium. The obstruction is the hydrodynamic rows of
  the wind, mass 1.0 to 1.4, momentum 1.7 to 2.0 and energy 1.5 to 2.0 at
  cells 330 to 420, not the carrier rows. This is the same behaviour the
  atomic cold-march snapshots show above, and the atomic cases escape it
  only by starting from an archived solved state, which the molecular groups
  do not have.

A molecular recipe therefore needs a near-stationary molecular seed, and the
restart contract says one can come only from a molecular run. The two ways
out are a code-side initializer that builds the molecular state from a
certified atomic one of the same He/H, and a base composition far enough
from the ceiling to let a march run (`q_H2_base` = 0.15 at He/H = 2.13, 21
percent below the ceiling, runs past step 10000 but has `du` 9.5 to 10.1 and
is a different lower boundary from the photochemical column's). **The first
of the two was adopted afterwards** (item L7, `EXHALE_MOLECULAR_SEED`): the
binary converts a certified atomic state of the same He/H, and `SEED_X2`
says what H2 the conversion carries -- by default, in every cell, the
smaller of the thermochemical fit and the root of each cell's own H2 carrier
row, that root solved for and not divided out
(`docs/input_schema.md` appendix D). The second is not adopted.

**Against the archived solution.** `models/compare_archive.py` overlays a
solved case on its 2026-08-30 counterpart; for the `K_zz` = 1e9,
He/H = 1.6261 case the figure is
`docs/figures/lhs1140b_wb_vs_archive_kzz1e9.pdf` (T, v, rho, the hydrogen
ionized fraction, `n(He 2^3S)` and the He I 10830 profile). MEASURED:

| | EW red pair [%A] | red depth [%] | FWHM [A] | log10 Mdot |
|---|---|---|---|---|
| archived 2026-08-30 | 1.11635 | 4.17649 | 0.2507 | 7.80 (READ) |
| this recipe | 1.17576 | 4.34388 | 0.2543 | 7.87 |

+5.3 percent in the equivalent width, +4.0 in the red depth, +1.4 in the
width and +0.07 dex in the mass-loss rate, after every change listed in
section 5 and with the base layer solved on a different discretization.
The depth and the width are the three-Gaussian fit of
`tpm_He10830_metrics.txt`, the equivalent width the integral over the
measurement's own window.
The archived curve carries the pre-2026-09-03 transit reader's ghost rows
(`archive_20260830/README.md`), which is part of that difference.

**What the recipe does not settle.** Three of the five items
`docs/lhs1140b_stationary_L5c_20260913.md` section 7 left open stand
unchanged, and the campaign's numbers carry them. MEASURED on the certified
He/H = 1.6261 state:

- the base velocity artifact survives: `v` = -0.111, -0.133, +0.087, +0.048,
  +0.075, +0.071 cm/s in cells 1 to 6, so the cell-centered `rho v r^2`
  reads -1.98, -2.24, +1.39, +0.72, +1.08, +0.97 of the far-wind value
  before settling to 1.000 farther out, while the face fluxes carry a
  constant flux (mass row 3.5e-10). This is P44 section 3;
- the base thermal jump is unchanged: cell 1 sits at 395 K against the
  226 K reservoir at r = 1 (P44 section 8);
- the grid convergence of the well-balanced solution has not been measured.

The `Mdot` question of that section is answered: the factor 3.1 against the
archived value was the missing secondary-ionization coupling, and with it
applied the difference is +0.07 dex.

**What one case costs.** MEASURED end to end on
`atomic_scalar_gj1132_kzz1e9/HeH1.60` at 8 threads, seed to line:
**15 min 15 s** -- 7 outer passes of 117 to 166 s each, then the
post-processing pass, the transit synthesis and the record. The case
certified (`info = 0`, log10 Mdot 7.87, red-pair EW 1.1702 %A, red depth
4.3260 %) and is left in the tree as the first case of the campaign.

**The closure rungs do not yet reach a certified state.** MEASURED on
`atomic_photochem_gj1132_kzzprofile/HeH9`, iteration 0 alone (`--kmax 1`,
16 threads): the Photochem column is produced in 70 s and the partitioned
route starts from the mapped seed as it does for a prescribed case, but only
its first outer pass converges. Passes 2 to 4 return `hydro info = 2` with
the energy row at 3.7e-2, 4.6e-2 and 4.9e-2 at cell 235, the element row
falling 5.55e-4, 2.64e-4, 1.95e-4, 1.70e-4, and the rung ends
`the stationary solve returned info = 1`, NOT CERTIFIED on two rows
(hydrodynamic energy 4.931e-2 against 1e-6; elemental transport 1.699e-4
against 1e-5 at a wind cell). The driver stops there, as it should. What is
established is that the route, the adapter and the seed all work over the
profile boundary; what is not is a converged rung, and the 15 closure cases
of section 3 wait on it.

**The runner writes its own record.** Every case that finishes carries
`REPRODUCE.md`, written by `models/write_reproduce.py` from the run that
produced it: the binary and its md5, the archived state the seed came from
and the mapping command, every command with its environment, the outer
passes and the certification as `run.log` states them, the line measures and
the clock. `models/status.py --write` collects one row per case from those
files into section 8 below, and `models/README.md` is the same account for
the tree as a whole.

**The turbulence twin of a case.** `docs/lhs1140b_exhale_vs_pwinds.tex` draws
the He I 10830 line with and without the turbulence term of
`exhale_transit_lib.py` (`v_turb = sqrt(5/6 kT/m)` added in quadrature to the
thermal width, after Lampon et al. 2020, which is what p-winds uses). Both
curves have to stand on one wind, so the second synthesis is run on the
solved case itself and written beside it:

```bash
cd models/atomic_scalar_gj1132_wellmixed/HeH0.42
mkdir -p tpm_turb
. ../../../winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$EXHALE EXHALE_TRANSIT_TURB=1 \
    EXHALE_TRANSIT_SAVE_PREFIX=tpm_turb/ OMP_NUM_THREADS=1 \
    python3 $EXHALE/EXHALE_transit.py > tpm_turb/transit.log 2>&1
```

It reads the case's own `output/` and `input.inp` and writes only into
`tpm_turb/`, so the canonical `tpm_*.txt` of the case are untouched. MEASURED
2026-09-14: 21 s per case at one thread. Regenerated 2026-09-16 on the
re-solved cases, together with `diffusion_check/`, since the L21 re-run
replaced every case directory. The well-mixed rungs up to
He/H = 1 carry it; no other case does.


### The numerical flux of the seven lowest-XUV atomic cases (decided 2026-09-17)

**Decision of the user, 2026-09-17.** The catalog is solved with
`Numerical flux: HLLC`, except the seven lowest-XUV atomic cases, which are
solved with `Numerical flux: ROE` on the same 500-cell grid. They are, and
the `flux` column of section 7 says so case by case:

```
atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7
atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7
atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7
atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7
atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13
atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7
atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7
```

**The reason.** The HLLC flux does not reach the stationary root on the
catalog grid at these XUV levels. The Roe flux does: with the catalog recipe
otherwise unchanged, three of the seven certify
(`x0.10_kzz1e9/HeH9.7` at outer pass 14, photochem `x0.20` at pass 6,
`x0.01_kzz1e9/HeH2.13` at pass 39) and photochem `x0.01/HeH9.7` reaches a
wind whose three hydrodynamic rows are all inside their tolerances and is
refused only by the gated elemental He/H partition at 10.85 R_p. That the
obstruction is the discretization and not the physics is settled by
resolution: the UNMODIFIED HLLC flux certifies the same `x0.10/HeH9.7` case
on a grid of twice the cells, and the state it reaches is the 500-cell Roe
state, agreeing to 1.0 per cent in Mdot, 0.4 per cent in the equivalent
width and 0.1 to 3 per cent in T and rho everywhere above 1.005 R_p. All of
this is MEASURED in `docs/lhs1140b_stationary_L25_20260916.md`, step 3,
sections 3.3 to 3.6.

**The systematic to state wherever a ROE case is quoted.** Measured on the
two controls that solve under both fluxes
(`atomic_scalar_gj1132_kzz1e9/HeH2.13` and
`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`), Roe against HLLC gives

| quantity | movement |
|---|---|
| Mdot | -0.31 and -0.62 per cent |
| He 10830 red-pair equivalent width | -0.31 and -0.44 per cent |
| first-cell temperature | -5.70 and -7.21 per cent |
| first-cell mass density | +5.86 and +7.52 per cent |
| mass-flux spread `du` over the wind window | unchanged to three digits |

The first cell moves in the two directions at nearly equal size because the
base pressure is held by the boundary, so a colder first cell is a denser
one. Every certified Roe state closes its element budget at the rounding of
the budget and keeps its critical point outside the 30 R_p domain, as the
HLLC states do.

**What this section does NOT claim.** The states of the seven cases are made
by the campaign re-run on the final binary; no output of the L25 measurement
directories (`models/.L25/`) is copied into a case directory. Four of the
seven were not certified under ROE when the decision was recorded, and each
of those carries its own `not_solved.md` with the row that holds it.

## 7. Results

Written by `models/status.py`; every number is read from the case directory.

The `archived claim on 59bfdb3f` column states what the no-step evaluate route answered when the stationary claim each ARCHIVED state file carried was re-measured with `EXHALE.x` md5 `59bfdb3fc4d0104fc2e9c3734596d2f6` (`models/CLAIMS_59bfdb3fc4d0.md`, PLAN_20260917 item L29): of its 119 states, 101 reproduced, 4 refused, 6 no_claim, 7 no_state, 1 not_evaluable. An archived `certified=T` header is a claim made by the executable that wrote it and is not a current acceptance. A row reads `superseded` where the state in the case directory has been written since that table was made, so that the archived verdict is about a file that is no longer there: the comparison is the `run=` stamp in the provenance header of the very file the table read, named in its own "claim file" column, against the modification time of `models/CLAIMS_59bfdb3fc4d0.md`, and that file's modification time where the header carries no stamp.

| group | He/H | flux | state | archived claim on 59bfdb3f | reason | log10 Mdot | red EW [%A] | red depth [%] | FWHM [A] | \|\|R\|\| |
|---|---|---|---|---|---|---|---|---|---|---|
| atomic_scalar_gj1132_wellmixed | 0.083 | HLLC | evaluated certified | superseded | -- | 7.840 | 0.3797 | 1.533 | 0.2331 | -- |
| atomic_scalar_gj1132_wellmixed | 0.40 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.1011 | 4.085 | 0.2533 | -- |
| atomic_scalar_gj1132_wellmixed | 0.42 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.1324 | 4.185 | 0.2543 | -- |
| atomic_scalar_gj1132_wellmixed | 0.44 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.1626 | 4.281 | 0.2553 | -- |
| atomic_scalar_gj1132_wellmixed | 0.55 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.3118 | 4.741 | 0.2598 | -- |
| atomic_scalar_gj1132_wellmixed | 1 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.6959 | 5.821 | 0.2734 | -- |
| atomic_scalar_gj1132_wellmixed | 10 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.9347 | 5.818 | 0.3113 | -- |
| atomic_scalar_gj1132_wellmixed | 100 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.9644 | 5.756 | 0.3193 | -- |
| atomic_scalar_gj1132_wellmixed | 1000 | HLLC | evaluated certified | superseded | -- | 7.850 | 1.9197 | 5.590 | 0.3214 | -- |
| **atomic_scalar_gj1132_wellmixed** | crossing He/H = 0.4044  (log-log chord between the rungs 0.4 and 0.42) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj699_wellmixed | 0.042 | HLLC | evaluated certified | superseded | -- | 8.570 | 0.9091 | 3.035 | 0.2815 | -- |
| atomic_scalar_gj699_wellmixed | 0.046 | HLLC | evaluated certified | superseded | -- | 8.570 | 0.9857 | 3.289 | 0.2815 | -- |
| atomic_scalar_gj699_wellmixed | 0.050 | HLLC | evaluated certified | superseded | -- | 8.570 | 1.0608 | 3.538 | 0.2820 | -- |
| atomic_scalar_gj699_wellmixed | 0.083 | HLLC | evaluated certified | superseded | -- | 8.570 | 1.6271 | 5.398 | 0.2835 | -- |
| atomic_scalar_gj699_wellmixed | 1 | HLLC | evaluated certified | superseded | -- | 8.550 | 4.5250 | 13.597 | 0.3133 | -- |
| atomic_scalar_gj699_wellmixed | 1000 | HLLC | evaluated certified | superseded | -- | 8.420 | 2.8725 | 7.510 | 0.3587 | -- |
| **atomic_scalar_gj699_wellmixed** | crossing He/H = 0.0526  (log-log chord between the rungs 0.05 and 0.083) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz0 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.840 | 0.0305 | 0.127 | 0.2270 | -- |
| atomic_scalar_gj1132_kzz0 | 2.6 | HLLC | evaluated certified | superseded | -- | 7.860 | 0.9549 | 3.611 | 0.2487 | -- |
| atomic_scalar_gj1132_kzz0 | 3.0 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.0798 | 4.026 | 0.2517 | -- |
| atomic_scalar_gj1132_kzz0 | 3.5 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2175 | 4.470 | 0.2558 | -- |
| atomic_scalar_gj1132_kzz0 | 3.7 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2698 | 4.635 | 0.2573 | -- |
| atomic_scalar_gj1132_kzz0 | 3.9 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3168 | 4.781 | 0.2583 | -- |
| **atomic_scalar_gj1132_kzz0** | crossing He/H = 3.1011  (log-log chord between the rungs 3 and 3.5) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e5 | 2.6 | HLLC | evaluated certified | superseded | -- | 7.870 | 0.9631 | 3.639 | 0.2487 | -- |
| atomic_scalar_gj1132_kzz1e5 | 3.0 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.0880 | 4.053 | 0.2522 | -- |
| atomic_scalar_gj1132_kzz1e5 | 3.35 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1862 | 4.371 | 0.2548 | -- |
| atomic_scalar_gj1132_kzz1e5 | 3.64 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2627 | 4.613 | 0.2568 | -- |
| atomic_scalar_gj1132_kzz1e5 | 3.93 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3314 | 4.826 | 0.2588 | -- |
| **atomic_scalar_gj1132_kzz1e5** | crossing He/H = 3.0706  (log-log chord between the rungs 3 and 3.35) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e6 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.840 | 0.0591 | 0.245 | 0.2270 | -- |
| atomic_scalar_gj1132_kzz1e6 | 2.4 | HLLC | evaluated certified | superseded | -- | 7.860 | 0.9465 | 3.583 | 0.2482 | -- |
| atomic_scalar_gj1132_kzz1e6 | 2.8 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.0804 | 4.028 | 0.2517 | -- |
| atomic_scalar_gj1132_kzz1e6 | 3.19 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1949 | 4.398 | 0.2553 | -- |
| atomic_scalar_gj1132_kzz1e6 | 3.46 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2669 | 4.626 | 0.2573 | -- |
| atomic_scalar_gj1132_kzz1e6 | 3.74 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3388 | 4.849 | 0.2593 | -- |
| **atomic_scalar_gj1132_kzz1e6** | crossing He/H = 2.8930  (log-log chord between the rungs 2.8 and 3.19) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e7 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.850 | 0.1585 | 0.652 | 0.2290 | -- |
| atomic_scalar_gj1132_kzz1e7 | 2.70 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1996 | 4.414 | 0.2553 | -- |
| atomic_scalar_gj1132_kzz1e7 | 2.94 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2724 | 4.643 | 0.2573 | -- |
| atomic_scalar_gj1132_kzz1e7 | 3.18 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3265 | 4.811 | 0.2588 | -- |
| **atomic_scalar_gj1132_kzz1e7** | crossing He/H = 2.5365  (log-log chord between the rungs 0.55 and 2.7) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e8 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.850 | 0.3089 | 1.252 | 0.2321 | -- |
| atomic_scalar_gj1132_kzz1e8 | 2.05 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1513 | 4.259 | 0.2538 | -- |
| atomic_scalar_gj1132_kzz1e8 | 2.23 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2206 | 4.480 | 0.2558 | -- |
| atomic_scalar_gj1132_kzz1e8 | 2.41 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.2826 | 4.675 | 0.2573 | -- |
| **atomic_scalar_gj1132_kzz1e8** | crossing He/H = 1.9729  (log-log chord between the rungs 0.55 and 2.05) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e9 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.850 | 0.4793 | 1.910 | 0.2361 | -- |
| atomic_scalar_gj1132_kzz1e9 | 1.50 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.0996 | 4.091 | 0.2522 | -- |
| atomic_scalar_gj1132_kzz1e9 | 1.60 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1511 | 4.257 | 0.2538 | -- |
| atomic_scalar_gj1132_kzz1e9 | 1.70 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.1970 | 4.404 | 0.2553 | -- |
| atomic_scalar_gj1132_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3749 | 4.958 | 0.2603 | -- |
| atomic_scalar_gj1132_kzz1e9 | 4.0 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.8229 | 6.209 | 0.2754 | -- |
| atomic_scalar_gj1132_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.860 | 2.0592 | 6.548 | 0.2946 | -- |
| **atomic_scalar_gj1132_kzz1e9** | crossing He/H = 1.5161  (log-log chord between the rungs 1.5 and 1.6) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e10 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.860 | 0.6463 | 2.529 | 0.2406 | -- |
| atomic_scalar_gj1132_kzz1e10 | 1.06 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.0456 | 3.912 | 0.2512 | -- |
| atomic_scalar_gj1132_kzz1e10 | 1.15 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.1032 | 4.100 | 0.2527 | -- |
| atomic_scalar_gj1132_kzz1e10 | 1.29 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.1867 | 4.369 | 0.2548 | -- |
| **atomic_scalar_gj1132_kzz1e10** | crossing He/H = 1.1580  (log-log chord between the rungs 1.15 and 1.29) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e11 | 0.55 | HLLC | evaluated certified | superseded | -- | 7.860 | 0.8065 | 3.099 | 0.2447 | -- |
| atomic_scalar_gj1132_kzz1e11 | 0.795 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.0302 | 3.859 | 0.2507 | -- |
| atomic_scalar_gj1132_kzz1e11 | 0.865 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.0866 | 4.044 | 0.2522 | -- |
| atomic_scalar_gj1132_kzz1e11 | 0.93 | HLLC | evaluated certified | superseded | -- | 7.860 | 1.1364 | 4.206 | 0.2538 | -- |
| **atomic_scalar_gj1132_kzz1e11** | crossing He/H = 0.8927  (log-log chord between the rungs 0.865 and 0.93) |  |  |  |  |  |  |  |  |  |
| atomic_scalarCNO_gj1132_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.870 | 1.3721 | 4.948 | 0.2603 | -- |
| **atomic_scalarCNO_gj1132_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132_kzzprofile | 2.09 -> 2.0924 | HLLC | evaluated certified | superseded | -- | 7.910 | 1.4765 | 5.321 | 0.2603 | -- |
| atomic_photochem_gj1132_kzzprofile | 3 -> 3.0034 | HLLC | evaluated certified | superseded | -- | 7.900 | 1.7640 | 6.156 | 0.2684 | -- |
| atomic_photochem_gj1132_kzzprofile | 5 -> 5.0057 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.0746 | 6.933 | 0.2805 | -- |
| atomic_photochem_gj1132_kzzprofile | 7 -> 7.0080 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.1768 | 7.102 | 0.2871 | -- |
| atomic_photochem_gj1132_kzzprofile | 8 -> 8.0091 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.2004 | 7.115 | 0.2896 | -- |
| atomic_photochem_gj1132_kzzprofile | 9 -> 9.0103 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.2137 | 7.111 | 0.2916 | -- |
| atomic_photochem_gj1132_kzzprofile | 9.7 -> 9.7111 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.2206 | 7.101 | 0.2931 | -- |
| atomic_photochem_gj1132_kzzprofile | 10 -> 10.0114 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.2224 | 7.094 | 0.2936 | -- |
| atomic_photochem_gj1132_kzzprofile | 12 -> 12.0137 | HLLC | evaluated certified | superseded | -- | 7.900 | 2.2289 | 7.037 | 0.2966 | -- |
| **atomic_photochem_gj1132_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.01_kzz1e9 | 2.13 | ROE | evaluated uncertified | superseded | stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 5 had completed; that pass had written no complete state of its own, so nothing new is published | 5.790 | 0.0000 | 0.000 | 0.2699 | -- |
| atomic_scalar_gj1132x0.01_kzz1e9 | 9.7 | ROE | evaluated uncertified | superseded | stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 4 had completed; that pass had written no complete state of its own, so nothing new is published | 5.780 | 0.0000 | 0.000 | 0.3461 | -- |
| **atomic_scalar_gj1132x0.01_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.10_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 6.820 | 0.0001 | 0.000 | 0.2815 | -- |
| atomic_scalar_gj1132x0.10_kzz1e9 | 9.7 | ROE | evaluated certified | superseded | -- | 6.810 | 0.2353 | 0.879 | 0.2487 | -- |
| **atomic_scalar_gj1132x0.10_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.15_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.000 | 0.0003 | 0.001 | 0.2865 | -- |
| atomic_scalar_gj1132x0.15_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.000 | 0.4514 | 1.671 | 0.2517 | -- |
| **atomic_scalar_gj1132x0.15_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.20_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.130 | 0.0022 | 0.007 | 0.2785 | -- |
| atomic_scalar_gj1132x0.20_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.130 | 0.6294 | 2.285 | 0.2573 | -- |
| **atomic_scalar_gj1132x0.20_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.25_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.240 | 0.0700 | 0.279 | 0.2351 | -- |
| atomic_scalar_gj1132x0.25_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.230 | 0.7913 | 2.823 | 0.2618 | -- |
| **atomic_scalar_gj1132x0.25_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.30_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.320 | 0.1894 | 0.760 | 0.2341 | -- |
| atomic_scalar_gj1132x0.30_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.310 | 0.9395 | 3.301 | 0.2664 | -- |
| **atomic_scalar_gj1132x0.30_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.33_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | -- | 7.360 | 0.2555 | 1.020 | 0.2351 | -- |
| atomic_scalar_gj1132x0.33_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | -- | 7.360 | 1.0231 | 3.566 | 0.2684 | -- |
| **atomic_scalar_gj1132x0.33_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.01_kzzprofile | 9.7 | ROE | info=1 uncertified hydrodynamic_refusal | superseded | stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 4 had completed; that pass had written no complete state of its own, so nothing new is published | 5.770 | 0.0000 | 0.000 | 0.2987 | 1.95E+00 |
| **atomic_photochem_gj1132x0.01_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.10_kzzprofile | 9.7 | HLLC | evaluated certified | superseded | -- | 6.840 | 0.2674 | 1.001 | 0.2487 | -- |
| **atomic_photochem_gj1132x0.10_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.15_kzzprofile | 9.7 | ROE | evaluated certified | superseded | -- | 7.020 | 0.4822 | 1.783 | 0.2527 | -- |
| **atomic_photochem_gj1132x0.15_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.20_kzzprofile | 9.7 | ROE | evaluated certified | superseded | -- | 7.160 | 0.6711 | 2.433 | 0.2578 | -- |
| **atomic_photochem_gj1132x0.20_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.25_kzzprofile | 9.7 | ROE | evaluated certified | superseded | -- | 7.260 | 0.8429 | 3.002 | 0.2623 | -- |
| **atomic_photochem_gj1132x0.25_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.30_kzzprofile | 9.7 | HLLC | evaluated certified | superseded | -- | 7.350 | 1.0097 | 3.542 | 0.2664 | -- |
| **atomic_photochem_gj1132x0.30_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.33_kzzprofile | 9.7 | HLLC | evaluated certified | superseded | -- | 7.390 | 1.0979 | 3.821 | 0.2689 | -- |
| **atomic_photochem_gj1132x0.33_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| molecular_scalar_gj1132_wellmixed | 0.083 | HLLC | evaluated uncertified | superseded | eighty outer passes reach a converged wind whose carrier H2 balance does not close, the row falling to 4.1e-02 over the first forty and rising again over the second forty, with the movement bound attained in every relaxation | 6.980 | 0.1164 | 0.450 | 0.2427 | -- |
| molecular_scalar_gj1132_wellmixed | 0.55 | HLLC | evaluated uncertified | superseded | thirty-one outer passes on a converged wind, the carrier H2 row falling to 3.8e-02, rising to 1.6e-01 and falling again to 6.5e-02 while its refusing cell walks outward, and the run was stopped at the six-hour ceiling | 7.410 | 0.6807 | 2.388 | 0.2674 | -- |
| molecular_scalar_gj1132_wellmixed | 2.13 | HLLC | not solved [no index] | superseded | the case was started for the first time, from an atomic well-mixed wind solved for it at this composition, and its carrier H2 row flattened at 4.3e-01 after seven outer passes; the run was stopped at the six-hour ceiling | -- | -- | -- | -- | 1.57E-08 |
| **molecular_scalar_gj1132_wellmixed** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |
| molecular_scalar_gj1132_kzz1e9 | 0.083 | HLLC | evaluated uncertified | superseded | one scalar movement bound over a column holding a slow H2 front and a far wind | 7.390 | 0.0003 | 0.001 | 0.3304 | -- |
| molecular_scalar_gj1132_kzz1e9 | 0.55 | HLLC | evaluated uncertified | refused | from the certified reference solution of its own group the carrier H2 row fell by a factor 12 over five outer passes and then rose again, with the hydrodynamic solve refusing every pass, and the run was stopped at the six-hour ceiling | 7.640 | 0.1591 | 0.631 | 0.2366 | -- |
| molecular_scalar_gj1132_kzz1e9 | 2.13 | HLLC | evaluated certified | superseded | stopped by the wall ceiling of 6h that models/run_campaign.sh applies (SIGTERM 21600 s after the campaign started the case), in the continuation pass after outer pass 2 had completed; that pass had written no complete state of its own, so nothing new is published | 7.900 | 1.5603 | 5.571 | 0.2628 | -- |
| molecular_scalar_gj1132_kzz1e9 | 9.7 | HLLC | evaluated certified | superseded | the stationary route ended info=1 with the three hydrodynamic rows inside their own tolerances, so what refuses is the composition: carrier balance H2: gated row measure  5.130E-04 above  1.0E-05 at cell 262 (a wind cell), and the outer loop spent its whole budget of 40 passes | 7.940 | 2.4061 | 7.629 | 0.2956 | -- |
| **molecular_scalar_gj1132_kzz1e9** | crossing He/H = 1.7386  (log-log chord between the rungs 0.55 and 2.13) |  |  |  |  |  |  |  |  |  |
| molecular_photochem_gj1132_kzzprofile | 2.09 | HLLC | not solved [no index] | superseded | from the certified wind of its own atomic pair the carrier H2 row bands between 1.9e-01 and 4.4e-01 over twelve outer passes without a trend, and the run was stopped at the six-hour ceiling | -- | -- | -- | -- | 3.27E-07 |
| molecular_photochem_gj1132_kzzprofile | 9 | HLLC | evaluated uncertified | superseded | the atomic-to-molecular conversion of this case's own certified atomic wind leaves the hydrodynamic rows at the base two decades outside their tolerances, and the element composition relaxation of outer pass 1 then found no admissible advance and restored its entry composition | 7.930 | 2.3332 | 7.505 | 0.2911 | -- |
| **molecular_photochem_gj1132_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |  |  |  |

## 8. How each model was reached

One row per case that has been run. The `record` column says what the published generation (`latest_complete`) is: a solve, or an evaluation of an earlier state (`evaluate-only`, with the wall clock of the evaluation). The seed, the outer passes, the verdict and the wall clock are those of the SOLVE behind the published state, the generation its parent chain reaches, READ from that solve's run record (`runs/<run_id>/run.json`), its log, or its own `REPRODUCE.md`; a fact none of them records reads `not recorded`. The wall clock of a solve is measured from the start of the run to the publication of that generation, so for a continuation it includes the solve it continued. The certification and the line measures are those of the published generation. The linked file carries the commands themselves.

| case | record | seed of the solve | outer passes of the solve | solve verdict | certification (published) | log10 Mdot | EW [%A] | wall clock of the solve | record | why not solved |
|---|---|---|---|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172348Z_563cd549, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.84 | 0.3797 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.083/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172358Z_3dc43708, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.1011 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.40/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172358Z_b138473d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.1324 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.42/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172401Z_6bb13986, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.1626 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.44/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172402Z_e564aecb, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.3118 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH1` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172425Z_86e2eb16, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.6959 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH1/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH10` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172428Z_2945aaa6, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.9347 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH10/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH100` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172429Z_fffea24a, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.9644 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH100/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172430Z_29caa5c3, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 1.9197 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH1000/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172522Z_09b39c54, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.57 | 0.9091 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.042/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172546Z_52e3280c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.57 | 0.9857 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.046/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172550Z_2e9c8c05, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.57 | 1.0608 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.050/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172552Z_583adbd7, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.57 | 1.6271 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.083/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH1` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172553Z_802e4c06, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.55 | 4.5250 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH1/REPRODUCE.md) | -- |
| `atomic_scalar_gj699_wellmixed/HeH1000` | evaluate-only, of `g0005` (0m24s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172555Z_a71277dd, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 8.42 | 2.8725 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH1000/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172023Z_dd1eb3a9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.84 | 0.0305 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172023Z_81e1e6bf, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 0.9549 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH2.6/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172024Z_12be6c0c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.0798 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.0/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172024Z_ecce5ecd, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2175 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.5/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172102Z_b4ccab06, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2698 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172102Z_86298bf1, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3168 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.9/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172144Z_b78ce7bd, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 0.9631 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH2.6/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172145Z_157b092e, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.0880 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.0/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172156Z_52c3eb10, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1862 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.35/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172157Z_6061f00d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2627 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.64/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172158Z_7e55e593, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3314 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.93/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172200Z_73a43ada, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.84 | 0.0591 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172222Z_554a5c53, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 0.9465 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH2.4/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172224Z_0631bee4, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.0804 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH2.8/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172224Z_b5780692, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1949 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.19/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172225Z_623a52f4, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2669 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.46/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172237Z_632e42c9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3388 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.74/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172238Z_2549c2a8, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 0.1585 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172239Z_5ef87d54, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1996 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH2.70/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172240Z_6d50b19c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2724 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH2.94/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172303Z_d319938d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3265 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH3.18/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172306Z_f3e22e12, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 0.3089 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172306Z_527a7f3f, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1513 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.05/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172307Z_303be6fe, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2206 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.23/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172317Z_3631af5d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.2826 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.41/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172319Z_5be1b90b, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.85 | 0.4793 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172320Z_703d936d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.0996 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.50/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172321Z_6f049297, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1511 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.60/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172345Z_f9303f17, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.1970 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.70/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T171710Z_4148d706, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3749 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172346Z_f6bc8d35, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.8229 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH4.0/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172347Z_affaa41f, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 2.0592 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172104Z_67190c7f, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 0.6463 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172104Z_f3251a6a, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.0456 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.06/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172114Z_2dc94d98, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.1032 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.15/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172116Z_9d0023c8, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.1867 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.29/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172116Z_d7233a7b, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 0.8065 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.55/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172119Z_44da648e, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.0302 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.795/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172143Z_a0cb6b6f, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.0866 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.865/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172143Z_3d4e7a3c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.86 | 1.1364 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.93/REPRODUCE.md) | -- |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m33s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172023Z_d06b2880, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.87 | 1.3721 | not recorded | [REPRODUCE.md](models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_85d9e44c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.91 | 1.4765 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH2.09/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH3` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171824Z_dceeca4c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 1.7640 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH3/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH5` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_e2f996ad, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.0746 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH5/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH7` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_4d42a042, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.1768 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH8` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_37190131, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.2004 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH8/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH9` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_10a40b7f, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.2137 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH9/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_22279af2, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.2206 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH10` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_5ba13fd3, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.2224 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH10/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132_kzzprofile/HeH12` | solve and its evaluation `g0005` | not recorded: the chain ends at the evaluate generation g0002_20260917T171925Z_7f4dd973, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 2.2289 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH12/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | solve and its evaluation `g0002` | not established; a candidate is in `provenance/g0001_20260917T172702Z_f98147ad/` | 1 | info=0 | uncertified | -- | -- | 0m31s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/REPRODUCE.md) | [not_solved.md](models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/not_solved.md) |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | solve and its evaluation `g0002` | not established; a candidate is in `provenance/g0001_20260917T172820Z_3a57917b/` | 1 | info=0 | uncertified | -- | -- | 0m33s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/REPRODUCE.md) | [not_solved.md](models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/not_solved.md) |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m20s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172438Z_42d68318, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 6.82 | 0.0001 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | evaluate-only, of `g0004` (0m21s) | not established; a candidate is in `provenance/g0001_20260917T172701Z_c89bc6ea/` | 1 | info=0 | certified | 6.81 | 0.2353 | 0m32s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m20s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172439Z_a27641f9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.00 | 0.0003 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172441Z_966e33f8, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.00 | 0.4514 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m20s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172442Z_c3f4d7d6, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.13 | 0.0022 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172505Z_b01f9ba9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.13 | 0.6294 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172511Z_cadde2d9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.24 | 0.0700 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172513Z_0178f44b, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.23 | 0.7913 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172513Z_db70ffa4, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.32 | 0.1894 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172516Z_be8e7118, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.31 | 0.9395 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | evaluate-only, of `g0005` (0m21s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172517Z_17179424, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.36 | 0.2555 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13/REPRODUCE.md) | -- |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | evaluate-only, of `g0005` (0m22s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172520Z_5b983ec7, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.36 | 1.0231 | not recorded | [REPRODUCE.md](models/atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | solve `g0004` | atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7 g0004 (tier3, He/H 9.71107), then a continuation | 14 | info=1 | uncertified | -- | -- | 160m29s to the publication of `g0004` at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/REPRODUCE.md) | [not_solved.md](models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved.md) |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | evaluate-only, of `g0005` (0m32s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172021Z_6d7b27c7, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 6.84 | 0.2674 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | evaluate-only, of `g0004` (0m33s) | not established; a candidate is in `provenance/g0001_20260917T172820Z_3203da7b/` | 1 | info=0 | certified | 7.02 | 0.4822 | 0m49s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | evaluate-only, of `g0004` (0m33s) | not established; a candidate is in `provenance/g0001_20260917T172703Z_4dbda7ca/` | 1 | info=0 | certified | 7.16 | 0.6711 | 0m47s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | evaluate-only, of `g0004` (0m33s) | not established; a candidate is in `provenance/g0001_20260917T172820Z_e33bb9ad/` | 1 | info=0 | certified | 7.26 | 0.8429 | 0m48s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | evaluate-only, of `g0005` (0m34s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172021Z_03b0446c, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.35 | 1.0097 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | evaluate-only, of `g0005` (0m34s) | not recorded: the chain ends at the evaluate generation g0002_20260917T172022Z_4f04b0e9, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.39 | 1.0979 | not recorded | [REPRODUCE.md](models/atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7/REPRODUCE.md) | -- |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | evaluate-only, of `g0004` (0m21s) | not established; a candidate is in `provenance/g0002_20260918T034603Z_258bc5d7/` | 40 | info=1 | uncertified | 6.98 | 0.1164 | 319m25s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH0.083/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_wellmixed/HeH0.083/not_solved.md) |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | evaluate-only, of `g0004` (0m23s) | not recorded: the chain ends at g0002_20260915T221824Z_d7347336, whose phase is 'unknown' and not a solve | not recorded | not recorded | uncertified | 7.41 | 0.6807 | not recorded | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH0.55/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_wellmixed/HeH0.55/not_solved.md) |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | no index | not recorded | not recorded | not recorded | uncertified | -- | -- | not recorded | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH2.13/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_wellmixed/HeH2.13/not_solved.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | evaluate-only, of `g0004` (0m19s) | recorded in the manifest: generation g0001_20260916T023329Z_32b34e57 of this same case | not recorded | info=1 | uncertified | 7.39 | 0.0003 | not recorded | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.083/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.083/not_solved.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | evaluate-only, of `g0005` (0m20s) | not recorded | not recorded | not recorded | uncertified | 7.64 | 0.1591 | not recorded | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.55/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.55/not_solved.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | evaluate-only, of `g0008` (0m23s) | not recorded: the chain ends at the evaluate generation g0002_20260917T174914Z_19d7ee4d, which names no parent generation (it was imported by models/import_legacy_states.py, and the run that solved the state it measured predates the index) | not recorded | not recorded | certified | 7.90 | 1.5603 | not recorded | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_kzz1e9/HeH2.13/not_solved.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | evaluate-only, of `g0008` (0m23s) | molecular_scalar_gj1132_kzz1e9/HeH9.7 g0003 (given) | 9 | info=0 | certified | 7.94 | 2.4061 | 6m45s to the publication of `g0004` at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH9.7/REPRODUCE.md) | [not_solved.md](models/molecular_scalar_gj1132_kzz1e9/HeH9.7/not_solved.md) |
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | no index | not recorded | not recorded | not recorded | uncertified | -- | -- | not recorded | [REPRODUCE.md](models/molecular_photochem_gj1132_kzzprofile/HeH2.09/REPRODUCE.md) | [not_solved.md](models/molecular_photochem_gj1132_kzzprofile/HeH2.09/not_solved.md) |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | evaluate-only, of `g0005` (0m34s) | not established; a candidate is in `provenance/g0003_20260917T222518Z_c55563e3/` | 1 | info=1 | uncertified | 7.93 | 2.3332 | 37m41s at 8 thread(s) | [REPRODUCE.md](models/molecular_photochem_gj1132_kzzprofile/HeH9/REPRODUCE.md) | [not_solved.md](models/molecular_photochem_gj1132_kzzprofile/HeH9/not_solved.md) |
## 9. What a case publishes, and which state a reader gets

Written 2026-09-18 with item D8 steps 3 to 7 of `docs/PLAN_20260918_rev2.md`,
in the shape approved that day (`docs/DECISION_D2a_D8_review.md` sections 3
and 4). Item D8a, section 6 above, is its first increment and is unchanged by
this: `output/` still holds the products, `solve_<n>/` still holds a
superseded solve whole, `ENDING` still says in one line how the published
solve ended.

### 9.1 The directories

A STATE DIRECTORY is one that carries an `input.inp` and an `output/`: a
case, a flux-closure iterate `k<NN>/`, or any other run directory of this
tree. It publishes its states like this:

```
<state directory>/
  output/                        the latest products, as before
  state_index.json               which generation a reader gets
  states/<generation_id>/        immutable once published
      Hydro_ioniz.txt  Ion_species.txt
      [Lyman_Werner.txt  OI_levels.txt  carrier_row_terms.txt]
      certification.txt          the certification block the certificate is
      manifest.json              section 9.2
  runs/<run_id>/                 one attempt: the `*.inp` files as the binary
                                 read them, `run.json`, which names the logs
                                 and the generations that attempt left, and
                                 `seed_identity.json`, the state that attempt
                                 was seeded with (section 9.12)
  provenance/<generation_id>/    what a LATER reading of the tree established
      provenance_recovered.json  about a generation whose manifest does not
      [provenance_recovered_v2.json ...]   carry it, beside the generation and
                                 never inside it (section 9.14)
```

`<generation_id>` is `g<NNNN>_<UTC of the moment the state was written>_<8
hex of the two halves>`: the counter orders the generations of one case, the
stamp is the state's own `# provenance: ... run=` field, and the hash makes
the name unique across retries and across two publishers that raced.

A published generation is immutable: its files and its directory are
read-only. NOTHING IS EVER DELETED, here or anywhere in this contract; a
retention policy is a separate decision and does not exist.

A generation holds the STATE, not the products derived from it. The
advection-corrected profiles, the heating and cooling breakdowns and the
transit curves are derived from a state and do not define one, so they stay
in `output/` and the manifest lists them by name and md5 under
`products_beside_the_state`. The logs stay where they are too: one run has
one log, and a copy of it beside every generation would be two records of one
thing. What the generation does carry of a log is the certification block,
copied into `certification.txt`, because that block IS the certificate.

The generations are copies and not hard links: the binary opens
`output/Hydro_ioniz.txt` with a truncating open, and a link to that file
would be truncated with it, so a link is not an immutable generation.
MEASURED 2026-09-18 after the import of section 9.6: the 134 `states/`
directories of the catalog take 132 MB together (`du -c`), against 565 MB in
the `output*` directories they were read from.

### 9.2 What the manifest records

Schema `exhale_state_generation/1`. Every field is READ from the state, from
the configuration beside it or from the log that wrote it, and nothing in it
is inferred:

| field | what it is |
|---|---|
| `generation_id`, `case_id`, `published_at`, `published_by` | the identity of this generation and of the publication |
| `state_written` | the `run=` stamp of the state's own header, which is part of the bytes and survives a copy while a file mtime does not |
| `parent` | the generation of THIS case the state descends from: the generation a continuation continued, or the generation the seed was taken from where the seed is one of this case's. Null where the seed came from another case or from no generation at all, and then the field says so and `seed` carries the identity |
| `seed` | what this generation was started from (section 9.12): which case and which generation, or which file pair where the directory publishes none, with the md5 of both halves; who chose it; and what was done to it before the binary read it. Never silently null: a publication whose seed could not be established carries `established: false` and the reason |
| `iteration_phase` | marching, stationary alternation, coupled block, evaluate, legacy import, or unknown; READ from the run's own `(input_read) Restart intent:` echo |
| `run_id` | the attempt that wrote the state, `runs/<run_id>/`, whose `run.json` carries its seed, its budget and the wall clock to each solve it published; null for a publication no runner made and for every generation published before 2026-09-19 (the publisher accepted the option and did not store it until then) |
| `source_identity` | the binary, its md5 and its `BINARY_MANIFEST_*` where one is beside it |
| `configuration_identity` | md5 of `input.inp`, `base.inp`, `metals.inp`, `opacity.inp`, the lower-atmosphere profile and the spectrum the input names |
| `grid` | rows, cells, ghost cells, `R0`, the radius range, the grid mode, and the width of the first physical cell MEASURED from this state |
| `units`, `species_schema` | the `# constants` line and the `# columns` and `# species_columns` lines of both halves |
| `model_identity` | the `# options` field (which is the equation of state and the reaction network as `load_IC` compares them), the `# reservoir` field, and `# boundary_model` and `# boundary_reservoir` where the state carries them (null for every state written before D5b-2, 2026-09-18) |
| `components` | name, md5, bytes, rows and columns of each file the generation holds |
| `storage_complete`, `storage_checks` | both halves present, one row count, one radius column to 1e-12, every value finite |
| `admissible` | what was checked of the state itself |
| `ending` | the ending class and its reason, from the `ENDING` file where there is one and otherwise from `classify_ending` of `models/run_case.sh`, with the source named and that file's md5 (`source_md5`) as it stood at the publication; and `solver`, the `info` and `||R||` the log states |
| `certification` | the verdict of the certification block, the entries that refuse, the log it was read from AND THAT LOG'S md5 AND BYTE COUNT AS THEY STOOD AT THE PUBLICATION (`source_md5`, `source_bytes`: a path is reused by the next run of the case and bytes are not, so the record is bound to the bytes it quotes), the md5 pair the certificate is attached to, how that attachment was established, and `stale_under` where a model change has made it historical |
| `state_claim` | `certified=`, `cert_reason=` and `mode=` of the state's own header |

### 9.3 The index, and completion against certification

`state_index.json`, schema `exhale_state_index/1`, carries two references and
they are NOT the same statement:

- `latest_complete`: the newest generation whose storage checks pass. It says
  the snapshot is whole, and nothing about the solve.
- `latest_certified`: the newest generation that ALSO carries a certificate
  for its own bytes. A CERTIFICATE IS A PROPERTY OF THE STATE, MEASURED:
  every active equation of the certification inventory evaluated on those
  bytes and within its tolerance, under a named binary (section 9.15). The
  index, the generations, the parentage and the ending classes RECORD that
  measurement and never confer it. A null reference means that no measurement
  of this case has been recorded, and says nothing about whether the state
  solves the equations. Two things must hold and both are about the published
  bytes: the state's own header claims `certified=T`, and a certification
  block says CERTIFIED. An uncertified publication moves `latest_complete`
  and never `latest_certified`, so a continuation that ends worse than the
  solve it continued cannot take the certified reference with it (the L34c
  finding). A certificate marked `stale_under` certifies nothing, and
  neither does a generation whose INDEX ENTRY carries `stale_under`: an
  evaluate pass of that binary refused it (section 9.8).

There is no `best` reference. Ranking two states needs a rule for which is
better, and that rule is not written yet.

### 9.4 One publisher, and what atomic means here

`models/publish_state.py` is the only thing that writes a generation or the
index. It publishes in this order: every component into an unpublished
directory, `fsync`ed; a re-measurement of what arrived against the source;
the manifest; the directory synchronized and renamed under its generation
name and made read-only; and only then a replacement index written to a
temporary file in the same directory and moved into place with `os.replace`.

`states/.publishing.lock` is held for the whole of it and a second publisher
on the same case is refused by name. ATOMIC VISIBILITY IS NOT CRASH
DURABILITY: a reader sees the old index or the new one, which is tested on
this filesystem (`models/tests/state_generations.sh`, check S1, stops the
publisher by PID between the temporary write and the rename); what a server
failure would leave is NOT tested and is not claimed.

### 9.5 The reader rule

Every reader resolves the index ONCE and takes every component of a state
from the one generation it names: `run_case.sh` (the evaluate pass reads the
generation the solve published, or `latest_certified` under `--evaluate`),
`pick_seed.py` (a candidate's metadata, its certification and the path it
prints, all from the `latest_certified` generation where the index names one,
since every tier ranks a case of this tree as certified; tier 0 labels its
line with the generation it handed over), `status.py`, `write_reproduce.py`,
`src/utils/map_state_to_grid.py` and `src/utils/element_flux_closure.py`.

An incomplete generation is refused. An older generation can be read only by
naming it (`EXHALE_STATE_GENERATION` for the mapper, `--generation` for the
publisher's own resolver), and the reader then says it was named.

A directory that publishes no index is read as it stands, and both halves
must carry the SAME suffix: `Hydro_ioniz.txt` with `Ion_species.txt`, or
`Hydro_ioniz_IC.txt` with `Ion_species_IC.txt`, never one of each. The
search for each half on its own that preceded this rule could map the product
of one generation beside the seed of another.

WHERE THE INDEX NAMES NO CERTIFIED GENERATION, a reader reads the absence of
a measurement and not a verdict on the state. `status.py` says `not measured
under binary <md5>` in the inventory's own column for such a case, and
`certified by measurement on <date>, binary <md5>` where a certificate is on
record; `pick_seed.py` says of such a case that the index records no
certifying measurement for it. `run_case.sh --evaluate` measures it: it takes
`latest_complete` where the index names no certified generation, and says so
in those words.

`status.py` reads the index and nothing else for a case that has one: the
solver verdict, the ending and the certification of its row are the published
generation's. An EVALUATION carries no solver verdict of its own, so such a
row reads `evaluated certified` rather than an `info` that belongs to the
solve. A case with no index says `[no index]` in the same column.

### 9.6 The warm restart this contract supports

A restart from a generation is a WARM restart: it preserves the physical
state and restarts the solver controls. It is not an exact algorithmic
continuation, and this contract does not offer one.

**Authoritative**: the conserved primitives and the species of the two data
files. The conserved energy is NOT stored separately, because D5a MEASURED
the physical column to round trip through the file at 1e-16 on `u(1)`,
`u(2)`, `u(3)`, the species fractions, `p` and `T`
(`docs/lhs1140b_stationary_D5b2_20260918.md` section 4), so the stored
primitives reconstruct it.

**Re-derived**: the two lower ghost rows, which `load_IC` no longer reads
under boundary model v1 (D5b-2): the ghost takes the reservoir's elemental
abundances with the first physical cell's partition, and its molecular
partition is solved. The face states, the caloric state and the derived
thermodynamics are rebuilt from the physical column and this run's own
reservoir.

**Restarted, not continued**: `dtau`, the trust bounds, the relaxation
parameters, the progress history and every counter. A continuation therefore
takes a different trajectory from an uninterrupted solve, which is one reason
the continuation of D8a is a decision about an ENDING and not about a number.

**What is a restart and what is a new initialization**, from the refusals of
`src/modules/files_IO/load_IC.f90` (the table of the D9 memo, section 2):
`metals`, `mol`, `oxychem` and `carrier` decide which species the files
carry, so a state whose columns are not this run's columns is a cold start
and never a restart; `he_diff` changes what the column's helium fraction
means; `iontrans` is a change of the equations that `load_IC` takes only when
the run names the token on a `Restart option change:` line, and a candidate
offered under it is a model-option transition whose certificate does not
transfer; the grid field is compared as text; and an element present in one
reservoir and absent from the other is a different composition whatever the
tolerance.

### 9.7 The states this tree already held

`models/import_legacy_states.py` entered them into the contract without
recomputing anything. MEASURED 2026-09-18: 136 state directories, 222
generations published, 134 indexes, 86 of them with a `latest_certified` and
48 with none.

A certificate is READ, never manufactured. No run of this tree recorded the
md5 of the pair it certified, so an imported certificate is attached by WRITE
ORDER and by agreement: the log that finished nearest the second the state
states it was written in, and whose verdict is the verdict the state's own
header carries, is the log of the pass that wrote it. That is a surrogate for
the md5 the contract asks for; every manifest that carries it says so in
`certification.identity_basis`, and where no log agrees the generation is
complete, uncertified, and says that no log of its directory states a verdict
about it. MEASURED, 87 of the 222 are in that position and 84 of those 87
are `output_pre_L34/` states, which the old post-processing route wrote with
one CFL step and `certified=F` over a solve that had certified; the log of
that pass was overwritten by the run that moved the state aside.

A state of another directory becomes a generation of a case only where
`models/legacy_state_map.txt` names it for that case, with the reason. A run
that was merely SEEDED from a case's state is not a generation of that case.

Two defects of the catalog the import exposed, both of them D8's own subject:

- `molecular_scalar_gj1132_kzz1e9/HeH0.083` published the state that the
  forty passes of item L33 STARTED from: its `output/Hydro_ioniz.txt` is byte
  for byte `models/.L22/i4_kz0083/output/Hydro_ioniz_IC.txt`, md5
  `1229b961323d7ac710146650da72454a`. The state those passes LEFT is on disk
  in that same directory, md5 `1ee5430a3c8560c40f017bc5eb618a8b`, and is now
  the case's second generation and its `latest_complete`, with the first as
  its parent; the manifest of the first says whose seed it is.
- `molecular_scalar_gj1132_wellmixed/HeH0.083` has no `ENDING` file: its
  ending was classified from its log by `classify_ending` of `run_case.sh`,
  and the manifest says that is where it came from. No case of the catalog
  carries an `ENDING` yet, because none has been run since item D8a; the
  field names the log for every one of them.

### 9.8 The evaluate-only entry, and a refused re-evaluation

Added 2026-09-19 after the D9 step 3 catalog refresh
(`docs/lhs1140b_catalog_refresh_20260919.md`), which had to reproduce the
runner's evaluate block in a scratch driver because the runner had no entry
for it.

    models/run_case.sh --evaluate <group>/<case>

solves nothing. It resolves the case's `latest_certified` generation (and
stops if the index names none), copies it as the `_IC` pair into
`runs/<run_id>/eval/output/`, writes `runs/<run_id>/eval/input.inp` from the
case's input with the relative paths made absolute, `Load IC? True` and
`Restart intent: stationary evaluate`, runs the binary there with its log in
`runs/<run_id>/pp.log`, and publishes the state the pass writes back as a
child `evaluate` generation whose parent is that `latest_certified`. Then the
products are installed in the case directory and the transit synthesis and
`REPRODUCE.md` follow as after a solve. `FORCE` is not needed.

BOTH ROUTES now run the evaluate pass in `runs/<run_id>/eval/`. The solve
route used `<case>/eval/`, which it emptied with `rm -rf` at every run; that
removal is gone, and a directory an earlier runner left at `<case>/eval/` is
no longer touched. Every case-level file the pass or the transit synthesis is
about to replace, and that differs from its replacement, is first copied into
`runs/<run_id>/superseded_case_products/` (the files of `output/` under
`output/`, and `pp.log`, `tpm_*.txt`, `transit.log`).

The `_IC` pair copied out of a generation is made writable (`chmod u+w`):
`cp` carries the read-only mode of a published generation, and the next seeded
run of the case opens `output/*_IC.txt` for writing in
`src/utils/map_state_to_grid.py`, which stopped with `PermissionError` in the
refresh. The same is done to the atomic pair the molecular seed conversion is
handed.

A REFUSED RE-EVALUATION TAKES THE CERTIFIED REFERENCE AWAY. When
`publish_state.py` publishes an `evaluate` generation whose certificate says
NOT CERTIFIED over a parent that is certified, it records on the PARENT'S
INDEX ENTRY `stale_under: <md5 of the binary that refused it>` and
`stale_evidence: <the evaluate generation>`; the parent's manifest is
immutable and is not touched. Where that parent is `latest_certified`, the
reference moves to the newest other generation certified under the same
binary, meaning an `evaluate` generation that binary certified or a
generation whose `evaluate` child it certified, and to null where no such
generation was evaluated. A solve's own certificate does not count as an
evaluation here. `publish_state.py reassess <case> <evaluate generation>`
applies the same rule, through the same locked index replacement, to an
`evaluate` generation published before the rule existed, and `verify` now
also refuses a `latest_certified` whose entry carries `stale_under`. The
symmetric half, an evaluation that CERTIFIES the state it was handed, is
section 9.15, and `reassess` applies that half too.

Applied 2026-09-19 with `reassess` to the three cases the refresh's binary
(`EXHALE_7670f310.x`, md5 `7670f31031fb4db91d27b44cb0da6f70`) refused; in
each the parent is the only certified generation, so `latest_certified` is
now null (MEASURED, `verify` passes on each):

| case | parent, now `stale_under` 7670f310... | refusing `evaluate` generation |
|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `g0001_20260917T172702Z_f98147ad` | `g0002_20260919T004834Z_d59ec9ff` |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `g0001_20260917T172820Z_3a57917b` | `g0002_20260919T004843Z_f7485b14` |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `g0001_20260917T172820Z_7e745859` | `g0002_20260919T004137Z_b2a43801` |

Sections 7 and 8 were regenerated by `models/status.py --write` afterwards and
came out byte for byte as before: their rows read the published generation,
`latest_complete`, which was already the refusing evaluation. What changed is
what a seed reader gets: `pick_seed.py` no longer offers these three as
certified states.

Tests: `models/tests/evaluate_entry.sh`, 18 checks on a synthetic tree with a
stand-in binary (E1 the writable `_IC` pair, E2 tier 0 and the certified
generation, E3 the evaluate-only entry and nothing removed, E4 the demotion
and `reassess`, E5 the record).

### 9.9 A run stopped from outside, the seed attempts, and what section 8 reads

Added 2026-09-19 (D9 step 3c, `docs/lhs1140b_catalog_refresh_20260919.md`,
section "Runner records and the HeH9.7 continuation").

A STOP FROM OUTSIDE IS RECORDED. `run_case.sh` runs every binary it starts in
the background and waits for it, and traps TERM, INT and HUP
(`stopped_from_outside`). The trap stops the binary by the PID it holds,
publishes a complete state the interrupted pass wrote itself (both halves,
newer than the pass, finite) and writes `ENDING`, `not_solved.md`,
`run.json` and `REPRODUCE.md` with the ending class
`stopped_by_wall_ceiling` when `run_campaign.sh` told it the ceiling
(`RUN_CASE_WALL`, `RUN_CASE_WALL_S`, `RUN_CASE_T0`) and the ceiling has been
reached, and `stopped_by_signal` otherwise. The campaign's own status line
uses the same class name. A record that replaces an earlier one keeps the
earlier one in `runs/<run_id>/superseded_case_products/REPRODUCE.md`; the
line measures of a record are quoted only from products the run itself
wrote.

ONE BUDGETED RE-SOLVE IS ONE SOLVE. `SEED_ATTEMPTS` (default 3) is a
ceiling, and an attempt after the first is taken only when the attempt
before it ended before its first outer pass (`another_seed_allowed`). The
campaign status line states `seeds=<used>/<allowed>` and its header states
the rule.

SECTION 8 IS ABOUT THE SOLVE. Its `record` column says whether the published
generation is a solve or an evaluation (`evaluate-only`, with the
evaluation's wall clock). Seed, outer passes, verdict and wall clock are those
of the solve generation the published generation's parent chain reaches
(`publish_state.solve_ancestor`), READ from that solve's `run.json`, from the
log whose certification block the generation carries verbatim, or from a
`REPRODUCE.md` of a solve whose run window holds the state's write time.
Where the chain ends at an imported evaluation with no parent, as for the
evaluated atomic catalog, the row says `not recorded` and why: the solve
records of those states were overwritten on 2026-09-18 before the index
existed.


### 9.10 How a run ended: four facts, one rule, two readers

Added 2026-09-19 (PLAN_20260919_rev1 item P4a, review section 6.1).

An ending is FOUR statements and any of them can hold while another fails:
the PROCESS ended, and by which of the binary's own end markers; the ROUTE
the run took stated its own outcome, and what it stated; a COMPLETE FINITE
state of that pass is on disk; and the state as written was CERTIFIED. They
are recorded separately, in `ending.class`, `ending.solver` and
`ending.evidence` of the manifest, and none of them is read off another.

THE ROUTE IS DETERMINED FIRST AND ONLY ITS OWN TERMINAL EVENT IS ITS VERDICT.
The binary echoes what it was asked to do (`Restart intent:`) and each route
states its outcome in its own line: the stationary route in
`the stationary solve returned info = N`, the evaluate route in
`the loaded state was measured ...` (it solves nothing and states no solver
verdict at all), and the marching route in its own stop line,
`-> converged:` (a stationary claim) or `-> stopped:` (no claim, a relaxation
snapshot). ONE STOP LINE COVERS SIX ENDINGS (`EXHALE_main.f90` near line
4192), so the binary's own reason is read with it: a run that asked for no
time integration at all is `no_integration`, a stop on a NaN in the conserved
state is `nonfinite_state`, and the rest are `marching_stop` with the reason
carried into the record. An inner ` (JFNK|PTC) done info=` line belongs to ONE PASS of the
route: the partitioned route writes one per outer pass and the marching route
one per hand-off attempt, marching on when the attempt fails. It is recorded
as `inner_info`, as evidence, and it is never the route's verdict. Reading it
as one let a run killed after a successful pass read as `solved`.

A STORED SOLVE IS NOT CALLED SUCCESSFUL WITHOUT A COMPLETE FINITE STATE. The
state on disk is tested before any verdict is honored, so `state_missing`,
`nonfinite_state` and `stale_state` are reached from `info = 0` as well as
from a refusal. `stale_state` is the state that is OLDER than the pass this
log records, which the runner can tell because it touches a mark when the
pass begins and hands `classify_ending` that mark; without a mark the
freshness is `unknown` and is never assumed.

CONFLICTING EVIDENCE IS REPORTED, NOT RECONCILED. `conflicting_evidence` is
the class of a run whose route claims a stationary solution and whose
certification of the state as written refuses it, `info = 0` with a refusing
certificate being that case; the conflict is stated in the reason and in
`ending.evidence.conflicts`. A certificate is taken as this state's only
where the block is the block of the state AS WRITTEN: a block a pass left
behind is a verdict on a state the run then moved away from.

ONE RULE, AND BOTH READERS CALL IT. `classify_ending` and `solver_verdict` of
`models/run_case.sh` are the rule; `models/publish_state.py` reaches them
through `report_ending` with `RUN_CASE_POLICY_ONLY=1` (`run_case_policy`,
`solver_from_log`, `ending_evidence`, `classify_with_run_case`) and carries no
reading of its own, and `models/status.py` reads the same answer through the
publisher. Until this item the publisher had a second reading that fell back
to the last inner line, so a manifest could carry a verdict the runner never
gave.

Tests: `models/tests/termination_classification.sh`, 14 checks on synthetic
logs built from the lines this catalog's logs carry, with no binary run: a
truncated log, a timeout after a successful inner pass, a missing output, a
nonfinite output, a stale output, `info = 0` with a refusing certificate, an
ordinary success, and three checks that the publisher answers what the runner
answers. MEASURED against the text before the item: 13 of the 14 fail.

WHAT THIS DOES NOT DO. No stored record is re-classified. MEASURED 2026-09-19
over the 398 generation records of the catalog, read only: 311 name a log
that is on disk, 306 of those are still the record's own log (the block the
generation carries verbatim is in it), and 44 of the 306 would be classified
differently, all of them `no_verdict` becoming `no_integration` (37 closure
rungs, 4 `diffusion_check` cases, 3 generations of
`molecular_scalar_gj1132_kzz1e9/HeH0.55`): every one of those runs printed
`-> stopped: "Do only PP" -- no time integration was requested`, so its state
is the state it was given and no march stands behind it. No record stored as
`solved` moves.
The remaining 5 comparable-looking records name a case-level `run.log` that a
later run has since overwritten, so they cannot be re-read at all: a record
that names a path and not a copy is bound to a file that is reused.

### 9.11 The case inventory: the list the counts are of

Added 2026-09-19 (PLAN_20260919_rev1 item P0).

    models/status.py inventory [--write]

writes `models/CASE_INVENTORY.md` and `models/CASE_INVENTORY.json`: the list
the catalog counts are OF, so a later count is reproduced by running the
command and not recalled. It reads the tree and touches no state.

THE SELECTION RULE. One row per state index: every `state_index.json` under
`models/`, found by walking the directory recursively, hidden directories
included. An index is what a reader resolves to get a state (section 9.5), so
the count of indexes is the count of states a reader can be handed. Nothing
is excluded: flux-closure rungs `k<NN>/`, `.L*` study directories,
`diffusion_check/`, preserved trees (`.stopped/`, `.ab/`) and archives
(`pre_*/`) each carry their own index and each is one row, with the `kind`
column saying which it is. A narrower count is taken by reading that column,
never by walking the tree again with another rule.

THE CRITERION FOR CERTIFIED. The index names a `latest_certified` generation
AND that generation's index entry carries no `stale_under` (nor
`certificate_stale_under`): a refused re-evaluation takes the certified
reference away from a state (section 9.8). A non-null reference is not a
fresh revalidation under every current option; it says the state certified
under the binary that evaluated it, which the inventory's `certified under`
column names.

MEASURED 2026-09-19 21:52 KST, binary of record `EXHALE_2c3b0acc.x`, md5
`2c3b0acc9983aec03bb4f844294fed18`: 137 state indexes, 83 certified. By kind:
84 ladder cases (74 certified), 46 closure rungs (9), 4 `diffusion_check` (0),
2 archives (0), 1 study directory (0). The plan's recursive count of the same
day reads 136 with 83 certified and the review's reads 137 with 83: the tree
gained one index during that day, `.L22/i3_alt`, published at
2026-09-19T09:53:59 and the only index of the catalog born on that date, so a
walk taken before that instant reads 136 and one taken after reads 137. The
certified count is the same in all three.

### 9.12 The seed a solve was given

Added 2026-09-19 (PLAN_20260919_rev1 item P4b). A state solved from another
state has a starting point, and that starting point is part of what the state
is: two solves of the same input from two different seeds are two different
histories, and a Newton solve reaches the fixed point of the basin it was
started in. Before this item the runner passed the publisher only
`FIRST_GENERATION`, the generation an internal continuation continued, so a
first solve from a chosen or a named seed published `parent.generation_id:
null` with "not recorded by the run that wrote this state". The solve
`g0004_20260919T103629Z_9b8394a3` of `molecular_scalar_gj1132_kzz1e9/HeH9.7`
is that case: its own run record names `g0003` of the same case as the state
it was handed, and its manifest names nothing.

THE SEED IS RESOLVED ONCE, BY THE RUNNER, AND BOTH RECORDS READ THAT ONE
RESOLUTION. `models/run_case.sh` writes it into
`runs/<run id>/seed_identity.json` (`record_seed_identity`, with
`seed_selection_note` for who chose it and `seed_transformation` for what was
done to it), and hands that same file to `models/publish_state.py`, which puts
it in the manifest of every generation the run publishes, and to
`models/write_reproduce.py`, which states it in `REPRODUCE.md`. The two
cannot name different seeds: there is one file, and neither tool forms a
reading of its own.

WHAT THE RECORD HOLDS, schema `exhale_seed_identity/1`:

| field | what it is |
|---|---|
| `established` | whether the seed is a fact this publication could read. False carries the reason in `note` and is never left a silent null |
| `kind` | a published generation, a file pair with no generation published for it, an internal continuation, a cold start, or not established |
| `case_id`, `case_directory` | the case the seed came from, which is not this case where the seed is an external one |
| `generation_id` | the generation, where the directory IS one. `generation_by_md5` instead where the pair matches a published generation of that case without being read from it, which is state identity and does not say which run wrote it |
| `path`, `components` | the directory the pair was read from and the md5 of both halves as they stood when the runner read them |
| `selected_by` | the tier of `models/pick_seed.py` that chose it, `SEED=` where the caller named it, `NOSEED=1`, or the runner continuing a generation |
| `transformation` | used as it stands, interpolated onto this case's cell centers by `src/utils/map_state_to_grid.py`, or converted from an atomic state to a molecular one by the binary. These are three different starting states, and the conversion is not an interpolation: it changes which species the state carries |
| `transformation_changed` | what that transformation changed: the grid header the pair was written onto, the elemental reservoir carried onto this case's composition, the base partition extension the molecular carriers were formed at |
| `initial_seed_of_the_run` | for an internal continuation, the seed the RUN itself was given, so the chain is not lost at the first continuation |

AND THE PARENT FOLLOWS FROM IT. Where the seed is a generation of THIS case,
that generation is the `parent`: a continuation records the generation it
continued, and a first solve from this case's own certified state records that
state. Where the seed is a generation of another case the parent stays null
and says why, because the parent chain is the chain this case's index can
resolve (section 9.8 follows it for a refused re-evaluation) and a foreign id
is not in it; the identity is in `seed`, with its case and its md5 pair.
Where the run was given no seed at all, a cold start is an established fact
and is recorded as one.

WHAT THE OLDER RECORDS HOLD, MEASURED 2026-09-19 and changed in nothing: of
398 published manifests, 226 carry a null parent. Six of them carry evidence
of their seed anywhere in the tree: one in its own run record (the `HeH9.7`
`g0004` above, an external seed of its own case, which the rule above would
have recorded as the parent) and five in the `REPRODUCE.md` of their case,
which speaks for the last run of that case (four molecular conversions and
one `SEED=`, each of another case). The remaining 220 name no run record at
all. Recovering them is item P4c and is not this item: an immutable manifest
is not rewritten to pretend the fact was there.

### 9.13 The seed of a flux-closure rung

Added 2026-09-20 (PLAN_20260919_rev1 item P4b-2). A rung `k<NN>/` is a state
directory (section 9.1), and a ladder is a chain: iteration 0 is started from
the state `models/run_closure.sh` mapped onto the current cell centers, and
every iteration after it from the iterate before it. Item P4b gave the
catalog runner the seed resolution of section 9.12 and left this one without
it.

THE RUNG IS PUBLISHED BY THE ONE PUBLISHER, WITH THE SEED IT WAS GIVEN.
After the closure driver returns, `models/run_closure.sh` publishes every
rung it left through `models/publish_state.py`, in the order the ladder ran
them, and resolves each rung's seed ONCE into
`k<NN>/runs/<run id>/seed_identity.json`, the same `exhale_seed_identity/1`
record section 9.12 defines:

| the rung | its seed | who chose it | what was done to it |
|---|---|---|---|
| `k00` | the directory this script mapped onto the current cell centers | the tier of `models/pick_seed.py` that named it, or `SEED=` | interpolated onto the cell centers of this tree by `src/utils/map_state_to_grid.py` |
| `k<NN>`, NN > 0 | the rung before it, as the published generation its index names | the runner, continuing the ladder | interpolated again where the driver carried the moved elemental reservoir or the moved R0 onto it (it leaves a `seed.log`, and `target_grid_Hydro_ioniz.txt` where R0 moved), used as it stands where neither had moved |

A RUNG'S PARENT STAYS NULL AND SAYS WHY. The rung before it is a state
directory of its own, so its generation is not in the chain this rung's index
resolves; the identity is in `seed`, with its case, its generation and the
md5 of both halves, by the same rule section 9.12 states for a seed of
another case.

The publication is idempotent: a rung whose index already names a generation
holding exactly these two halves is not published again, so a resumed or
re-run ladder adds a generation only where the bytes are new. Where
`models/run_case.sh` cannot be read, no rung is published and the runner says
so: the ending of a solve is classified by that one rule and by no other.

Tests: `models/tests/closure_seed_provenance.sh`, 22 checks on a synthetic
ladder with no binary and no closure driver run (`RUN_CLOSURE_PUBLISH_ONLY=1`
defines the publication and runs nothing, as `RUN_CASE_POLICY_ONLY=1` does for
the catalog runner).

### 9.14 Provenance recovered after the fact

Added 2026-09-20 (PLAN_20260919_rev1 item P4c). A manifest is IMMUTABLE and is
never rewritten to pretend a fact was in it. What a later reading of the tree
establishes about a generation published before the publisher carried a seed
goes BESIDE the generation:

    <state directory>/provenance/<generation_id>/provenance_recovered.json

written by `models/recover_provenance.py` and by nothing else. Nothing is
deleted: a later recovery that says something different writes
`provenance_recovered_v2.json` beside the first and names what it supersedes,
and the readers take the newest.

EVERY ATTACHMENT CARRIES: what was recovered, the file it was read from with
that file's md5 and byte count, the rule that was applied, the date, and one
of three confidences.

| confidence | when |
|---|---|
| `established` | the evidence names THIS generation and states the seed: the manifest's own `parent` (which is RECORDED there, and the attachment says so), or a run record of the directory whose `solve_generation`, `solve_generations` or `evaluate_generation` names this generation and whose `seed` field is not empty |
| `inferred from a single match` | exactly one `REPRODUCE.md` of the case describes a run whose window holds the moment this state was written, and it states a seed. A window is not a name |
| `not established` | no evidence, or several candidates that the run, the configuration and the time do not separate. THE AMBIGUITY IS KEPT: every candidate is listed, because an identical state pair is evidence of state identity and does not identify a run |

A READER MAY USE A RECOVERED FACT ONLY WHERE IT IS `established`, AND SAYS IT
IS RECOVERED. `publish_state.recovered_provenance` is the one entry and
returns nothing for a weaker confidence unless the caller is reporting ON the
recovery (`whatever_the_confidence=True`). `models/status.py` prints such a
seed as `recovered: ...` in section 8, or as `recorded in the manifest: ...`
where the attachment restates the manifest's own parent, and where the
attachment is weaker it
says `not established; a candidate is in provenance/<generation>/` instead of
stating a seed from a window match. `models/pick_seed.py
--seed-provenance` appends to each candidate what that state was itself
started from, marked `recorded:` where the manifest carries it and
`recovered (...)` where the attachment does; without the option its lines are
unchanged.

THE LOG A RECORD NAMES. A record that names a case-level `run.log` or
`pp.log` is bound to a path the next run of that directory reuses. Every
attachment states, MEASURED when it is written, whether the file at that path
still carries the certification block the generation holds verbatim, and says
in as many words where it does not. A publication made from now on records
the md5 of that log at the moment it publishes (section 9.2), so the
statement no longer depends on the file surviving.

MEASURED 2026-09-20 over the 398 published generations of 137 state
directories, and NOTHING WAS REWRITTEN: 176 established (172 of them the
manifest's own parent, and 4 a run record naming the generation), 10 inferred
from a single match, 212 not established. Of the 85 records whose named log
no longer carries their block, 5 name a case-level `run.log` and 80 name a
`pp.log`; each of the 85 says so in its attachment.

Tests: `models/tests/provenance_recovery.sh`, 34 checks, no binary run.


### 9.15 A state that passes the certificate is certified

Added 2026-09-20 after the user's rule of that day, over the three molecular
states of section 9.8's own subject: a state is certified because it
satisfies the equations of the model to the stated tolerances, MEASURED, and
never because a procedure was or was not followed.

THE STATEMENT. A certificate is a measurement of the state AS WRITTEN under a
named binary: every active equation of the certification inventory evaluated
on those bytes and within its tolerance. It is not a statement about a run.
The generations, the index, the parentage and the ending classes are
bookkeeping: they record the measurement, and they are never its criterion.

THE RULE HAS TWO HALVES and `publish_state.apply_evaluation_verdict` holds
both, beside the solve path of `update_index`:

- **The solve path.** A generation a solve published, whose own block
  certifies the bytes it wrote, is `latest_certified`, whatever the runner's
  ending class for that run was: a pass ceiling, a wall ceiling, a stop from
  outside, or a continuation refused afterwards. `certified_here` reads the
  certificate, the state's own `certified=` claim and the storage checks, and
  reads no ending. MEASURED 2026-09-20: this was already what the publisher
  did, and the check `C7` of `models/tests/certified_evaluation.sh` passes
  on the text as it stood before this section was written.
- **The evaluate path**, the symmetric half of the demotion rule of section
  9.8. An `evaluate` generation whose certificate says CERTIFIED, over a
  parent that is complete, becomes `latest_certified` itself. The parent's
  index entry records `certified_under: <md5 of the binary that measured it>`
  and `certified_evidence: <the evaluate generation>`; the child's entry
  records `certificate_from: an evaluation of <parent>` and
  `certificate_binary_md5`. THE PARENT'S OWN ENDING IS NOT CHANGED AND ITS
  MANIFEST IS NOT REWRITTEN: a manifest is immutable, and the index is where
  the statement is recorded, exactly as for a refusal.

THE ENDING AND THE CERTIFICATE STAY TWO FACTS (item P4a). The ending says how
the run ended, the certificate says what the state is; both are recorded and
neither overrides the other.

THE GUARD RAILS, each of them tested:

- the certificate must be of the bytes the generation carries
  (`certification.assessed_state`, the md5 pair published);
- the state must be complete and finite (`storage_complete`): an incomplete
  or nonfinite state is never promoted, whatever its block says;
- an `evaluate` generation must be an evaluation OF ITS PARENT. The pair the
  pass was handed is the `_IC` pair the runner leaves beside the state it
  wrote; its md5 pair is recorded in the manifest as `evaluated_state`, and
  it must be the parent's two halves BITWISE
  (`evaluated_pair_is_the_parent`);
- a refused evaluation still demotes, by section 9.8.

WHICH STATE `--evaluate` MEASURES. Until this section the entry resolved
`latest_certified` and stopped where the index named none, so a case whose
solve ended uncertified could not have its state measured through the runner
at all, and a measurement made outside the runner published nothing. That
was a defect of the procedure and not of the states. The entry now reads:

    models/run_case.sh --evaluate [--generation <id>] <group>/<case>

`--generation <id>` measures that generation of the case. With neither, the
entry takes `latest_certified` where the index names one and `latest_complete`
where it does not, and says in its output which reference it took and why.
The rest is unchanged: the pass runs in `runs/<run_id>/eval/`, the state it
writes back is published as an `evaluate` child of the generation it
measured, and the products, the transit synthesis and `REPRODUCE.md` follow.

Tests: `models/tests/certified_evaluation.sh`, 17 checks on a synthetic tree
with a stand-in binary (C1 the named generation, C2 the fallback and its
words, C3 the promotion and the untouched parent, C4 the evaluation-of-its-
parent guard, C5 the nonfinite state, C6 the refusal still demoting, C7 the
solve path under a wall ceiling, C8 `--generation` refused without
`--evaluate`). MEASURED 2026-09-20: 3 of 17 pass on the text before this
section, 17 of 17 after, and the eleven existing fixtures pass unchanged.

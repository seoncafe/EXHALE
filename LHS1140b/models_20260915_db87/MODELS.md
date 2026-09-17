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
                            REPRODUCE.md (how that case was reached)
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
| `molecular_scalar_gj1132_wellmixed` | 0.083, 0.55, 2.13 | NEW: the molecular layer solved in EXHALE above the scalar base (`base.inp` with `q_H2_base` from the photochemical column's H2 fraction, `p_base` 1 microbar) |
| `molecular_scalar_gj1132_kzz1e9` | 0.083, 0.55, 2.13, 9.7 | NEW: the same with element diffusion |
| `molecular_photochem_gj1132_kzzprofile` | 2.09, 9 | NEW: the molecular layer above the converged column of ITS OWN atomic closure rung (`atomic_photochem_gj1132_kzzprofile/HeH2.09/k<last>/lower_atmosphere_profile.dat` and `HeH9/k<last>/...`), held fixed, no closure. **Why the rung's column and not a stored one** (changed 2026-09-16): a molecular case of this group is solved from the certified wind of its atomic pair, and `load_IC` admits a seed only when the two states agree on the column and on the reservoir to a part in 1e6. The 2026-08-30 stored columns carry He/H = 2.092388551 and 9.010268395 at the matching level only by coincidence of their own history; measured on 2026-09-16, the stored 2.09 column and the re-run 2.09 rung differ by 1.2e-04 in He/H, twenty times the threshold, and the seed was refused. The case name is now the rung's name and the He/H line carries the value that column actually holds, so the pair cannot drift apart again. (The former `HeH9.05` directory named the stored column's value and has been removed; it held no result.) |

The `He/H = 0.55` case of each `kzz*` group is the K_zz sensitivity scan
of `kzz_decision.md` section 3 at one composition.

## 4. The molecular base

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
names the binary, its md5 and the host. 79 of the 83 cases certify. The four
that do not are still solving at the time sections 7 and 8 below were
written: `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` and the photochemical
cases at 0.15, 0.20 and 0.25 of the fiducial spectrum, all of them held by
the hydrodynamic rows under a composition that is still travelling, at the
raised pseudo-time start as at the default (the memo
`docs/lhs1140b_rerun_20260915.md` section 4 measures both). Their rows below
state what they stand at. Each
case's `REPRODUCE.md` names it and its md5 again. The three cases at 0.01 of
the fiducial XUV were not attempted: they did not solve before (item L17,
diagnosed as the HLLC contact selection on a base subsonic to Mach 1e-8) and
nothing in L21 touches that. The nine molecular cases are not part of this
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
| 0, the seed | `models/pick_seed.py` names the solved state of the same physics (spectrum, boundary kind, diffusion, `K_zz`) nearest in log He/H: a CERTIFIED case of the same group first, then a certified case of another group with the same physics, then the archive, newest generation first. `src/utils/map_state_to_grid.py` interpolates it onto `models/current_grid_Hydro_ioniz.txt`, the cell centers of the current code, and writes it as the case's `*_IC.txt`; where the seed was solved at another He/H, `--reservoir He/H` carries its helium onto this case's before it is written, so that the base rows arrive at this composition and a diffused He/H profile keeps its shape. Skipped where the case states `Load IC? False`. Where nothing carries the case's XUV normalization (the 0.30-scaled closure rung is the one such point) the nearest scaling of the same star and the same spectral shape is taken instead, and the line says so. Which tier the seed came from is on the line, and in `REPRODUCE.md` |
| 1, the wind | the binary, with `EXHALE_PTC_DTAU0=1.0` and `OMP_NUM_THREADS` (8 by default). The verdict is `the stationary solve returned info = 0` and the last certification block; `info = 2` is a state written but not certified |
| 2, the profiles | the same binary on its own solved state (`Load IC? True`, `Do only PP: True`, `CFL: 1.0e-12`, no `Restart intent` and no `Solver`), which writes the `*_adv.txt` files and the steady-state `Mdot` |
| 3, the line | `EXHALE_transit.py` through the WINERED HIRES-Y kernel of the measurement (`LHS1140b/winered_hires_y.sh`) |

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

**Why the post-processing pass runs at `CFL: 1.0e-12`.** It takes one time
step before it stops. MEASURED on the certified state above, the state the
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
2026-09-14: 21 s per case at one thread. The well-mixed rungs up to
He/H = 1 carry it; no other case does.


## 7. Results

Written by `models/status.py`; every number is read from the case directory.

| group | He/H | state | log10 Mdot | red EW [%A] | red depth [%] | FWHM [A] | \|\|R\|\| |
|---|---|---|---|---|---|---|---|
| atomic_scalar_gj1132_wellmixed | 0.083 | info=0 certified | 7.840 | 0.3797 | 1.533 | 0.2331 | 2.49E-08 |
| atomic_scalar_gj1132_wellmixed | 0.40 | info=0 certified | 7.850 | 1.1011 | 4.085 | 0.2533 | 1.75E-08 |
| atomic_scalar_gj1132_wellmixed | 0.42 | info=0 certified | 7.850 | 1.1324 | 4.185 | 0.2543 | 1.56E-08 |
| atomic_scalar_gj1132_wellmixed | 0.44 | info=0 certified | 7.850 | 1.1626 | 4.281 | 0.2553 | 1.50E-08 |
| atomic_scalar_gj1132_wellmixed | 0.55 | info=0 certified | 7.850 | 1.3118 | 4.742 | 0.2598 | 1.78E-08 |
| atomic_scalar_gj1132_wellmixed | 1 | info=0 certified | 7.860 | 1.6959 | 5.821 | 0.2734 | 1.12E-08 |
| atomic_scalar_gj1132_wellmixed | 10 | info=0 certified | 7.850 | 1.9347 | 5.818 | 0.3113 | 3.52E-08 |
| atomic_scalar_gj1132_wellmixed | 100 | info=0 certified | 7.860 | 1.9645 | 5.757 | 0.3193 | 1.13E-08 |
| atomic_scalar_gj1132_wellmixed | 1000 | info=0 certified | 7.850 | 1.9197 | 5.590 | 0.3214 | 1.47E-08 |
| **atomic_scalar_gj1132_wellmixed** | crossing He/H = 0.4044  (log-log chord between the rungs 0.4 and 0.42) |  |  |  |  |  |  |
| atomic_scalar_gj699_wellmixed | 0.042 | info=0 certified | 8.570 | 0.9091 | 3.035 | 0.2815 | 9.98E-08 |
| atomic_scalar_gj699_wellmixed | 0.046 | info=0 certified | 8.570 | 0.9858 | 3.289 | 0.2815 | 9.54E-08 |
| atomic_scalar_gj699_wellmixed | 0.050 | info=0 certified | 8.570 | 1.0609 | 3.538 | 0.2820 | 1.71E-07 |
| atomic_scalar_gj699_wellmixed | 0.083 | info=0 certified | 8.570 | 1.6271 | 5.398 | 0.2835 | 7.46E-08 |
| atomic_scalar_gj699_wellmixed | 1 | info=0 certified | 8.550 | 4.5250 | 13.597 | 0.3133 | 4.03E-08 |
| atomic_scalar_gj699_wellmixed | 1000 | info=0 certified | 8.420 | 2.8725 | 7.510 | 0.3587 | 7.99E-09 |
| **atomic_scalar_gj699_wellmixed** | crossing He/H = 0.0526  (log-log chord between the rungs 0.05 and 0.083) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz0 | 0.55 | info=0 certified | 7.840 | 0.0305 | 0.127 | 0.2270 | 1.71E-08 |
| atomic_scalar_gj1132_kzz0 | 2.6 | info=0 certified | 7.860 | 0.9549 | 3.611 | 0.2487 | 1.18E-08 |
| atomic_scalar_gj1132_kzz0 | 3.0 | info=0 certified | 7.870 | 1.0798 | 4.026 | 0.2517 | 1.79E-08 |
| atomic_scalar_gj1132_kzz0 | 3.5 | info=0 certified | 7.870 | 1.2175 | 4.470 | 0.2558 | 1.56E-08 |
| atomic_scalar_gj1132_kzz0 | 3.7 | info=0 certified | 7.870 | 1.2698 | 4.635 | 0.2573 | 1.40E-08 |
| atomic_scalar_gj1132_kzz0 | 3.9 | info=0 certified | 7.870 | 1.3168 | 4.781 | 0.2583 | 2.73E-08 |
| **atomic_scalar_gj1132_kzz0** | crossing He/H = 3.1010  (log-log chord between the rungs 3 and 3.5) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e5 | 2.6 | info=0 certified | 7.870 | 0.9631 | 3.639 | 0.2487 | 1.60E-08 |
| atomic_scalar_gj1132_kzz1e5 | 3.0 | info=0 certified | 7.870 | 1.0880 | 4.053 | 0.2522 | 1.12E-08 |
| atomic_scalar_gj1132_kzz1e5 | 3.35 | info=0 certified | 7.870 | 1.1862 | 4.371 | 0.2548 | 1.58E-08 |
| atomic_scalar_gj1132_kzz1e5 | 3.64 | info=0 certified | 7.870 | 1.2627 | 4.613 | 0.2568 | 2.12E-08 |
| atomic_scalar_gj1132_kzz1e5 | 3.93 | info=0 certified | 7.870 | 1.3314 | 4.826 | 0.2588 | 1.52E-08 |
| **atomic_scalar_gj1132_kzz1e5** | crossing He/H = 3.0706  (log-log chord between the rungs 3 and 3.35) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e6 | 0.55 | info=0 certified | 7.840 | 0.0591 | 0.245 | 0.2270 | 1.63E-08 |
| atomic_scalar_gj1132_kzz1e6 | 2.4 | info=0 certified | 7.860 | 0.9465 | 3.583 | 0.2482 | 1.38E-08 |
| atomic_scalar_gj1132_kzz1e6 | 2.8 | info=0 certified | 7.870 | 1.0804 | 4.028 | 0.2517 | 1.43E-08 |
| atomic_scalar_gj1132_kzz1e6 | 3.19 | info=0 certified | 7.870 | 1.1949 | 4.398 | 0.2553 | 2.53E-08 |
| atomic_scalar_gj1132_kzz1e6 | 3.46 | info=0 certified | 7.870 | 1.2669 | 4.626 | 0.2573 | 1.87E-08 |
| atomic_scalar_gj1132_kzz1e6 | 3.74 | info=0 certified | 7.870 | 1.3388 | 4.849 | 0.2593 | 1.39E-08 |
| **atomic_scalar_gj1132_kzz1e6** | crossing He/H = 2.8930  (log-log chord between the rungs 2.8 and 3.19) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e7 | 0.55 | info=0 certified | 7.850 | 0.1585 | 0.652 | 0.2290 | 2.11E-08 |
| atomic_scalar_gj1132_kzz1e7 | 2.70 | info=0 certified | 7.870 | 1.1996 | 4.414 | 0.2553 | 1.76E-08 |
| atomic_scalar_gj1132_kzz1e7 | 2.94 | info=0 certified | 7.870 | 1.2724 | 4.643 | 0.2573 | 1.76E-08 |
| atomic_scalar_gj1132_kzz1e7 | 3.18 | info=0 certified | 7.870 | 1.3265 | 4.811 | 0.2588 | 4.26E-08 |
| **atomic_scalar_gj1132_kzz1e7** | crossing He/H = 2.5365  (log-log chord between the rungs 0.55 and 2.7) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e8 | 0.55 | info=0 certified | 7.850 | 0.3089 | 1.252 | 0.2321 | 1.75E-08 |
| atomic_scalar_gj1132_kzz1e8 | 2.05 | info=0 certified | 7.870 | 1.1513 | 4.259 | 0.2538 | 2.40E-08 |
| atomic_scalar_gj1132_kzz1e8 | 2.23 | info=0 certified | 7.870 | 1.2206 | 4.480 | 0.2558 | 1.44E-08 |
| atomic_scalar_gj1132_kzz1e8 | 2.41 | info=0 certified | 7.870 | 1.2826 | 4.675 | 0.2573 | 1.40E-08 |
| **atomic_scalar_gj1132_kzz1e8** | crossing He/H = 1.9729  (log-log chord between the rungs 0.55 and 2.05) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e9 | 0.55 | info=0 certified | 7.850 | 0.4793 | 1.910 | 0.2361 | 1.56E-08 |
| atomic_scalar_gj1132_kzz1e9 | 1.50 | info=0 certified | 7.870 | 1.0996 | 4.091 | 0.2522 | 9.88E-09 |
| atomic_scalar_gj1132_kzz1e9 | 1.60 | info=0 certified | 7.870 | 1.1511 | 4.257 | 0.2538 | 1.59E-08 |
| atomic_scalar_gj1132_kzz1e9 | 1.70 | info=0 certified | 7.870 | 1.1970 | 4.404 | 0.2553 | 1.63E-08 |
| atomic_scalar_gj1132_kzz1e9 | 2.13 | info=0 certified | 7.870 | 1.3749 | 4.958 | 0.2603 | 1.15E-08 |
| atomic_scalar_gj1132_kzz1e9 | 4.0 | info=0 certified | 7.870 | 1.8229 | 6.209 | 0.2754 | 1.34E-08 |
| atomic_scalar_gj1132_kzz1e9 | 9.7 | info=0 certified | 7.860 | 2.0592 | 6.548 | 0.2946 | 1.04E-08 |
| **atomic_scalar_gj1132_kzz1e9** | crossing He/H = 1.5161  (log-log chord between the rungs 1.5 and 1.6) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e10 | 0.55 | info=0 certified | 7.860 | 0.6463 | 2.529 | 0.2406 | 1.47E-08 |
| atomic_scalar_gj1132_kzz1e10 | 1.06 | info=0 certified | 7.860 | 1.0456 | 3.912 | 0.2512 | 1.68E-08 |
| atomic_scalar_gj1132_kzz1e10 | 1.15 | info=0 certified | 7.860 | 1.1032 | 4.100 | 0.2527 | 1.77E-08 |
| atomic_scalar_gj1132_kzz1e10 | 1.29 | info=0 certified | 7.860 | 1.1867 | 4.369 | 0.2548 | 1.20E-08 |
| **atomic_scalar_gj1132_kzz1e10** | crossing He/H = 1.1579  (log-log chord between the rungs 1.15 and 1.29) |  |  |  |  |  |  |
| atomic_scalar_gj1132_kzz1e11 | 0.55 | info=0 certified | 7.860 | 0.8065 | 3.099 | 0.2447 | 1.88E-08 |
| atomic_scalar_gj1132_kzz1e11 | 0.795 | info=0 certified | 7.860 | 1.0302 | 3.859 | 0.2507 | 2.10E-08 |
| atomic_scalar_gj1132_kzz1e11 | 0.865 | info=0 certified | 7.860 | 1.0866 | 4.044 | 0.2522 | 1.93E-08 |
| atomic_scalar_gj1132_kzz1e11 | 0.93 | info=0 certified | 7.860 | 1.1364 | 4.206 | 0.2538 | 1.42E-08 |
| **atomic_scalar_gj1132_kzz1e11** | crossing He/H = 0.8927  (log-log chord between the rungs 0.865 and 0.93) |  |  |  |  |  |  |
| atomic_scalarCNO_gj1132_kzz1e9 | 2.13 | info=0 certified | 7.870 | 1.3721 | 4.948 | 0.2603 | 1.54E-08 |
| **atomic_scalarCNO_gj1132_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132_kzzprofile | 2.09 -> 2.0924 | info=0 certified | 7.910 | 1.4765 | 5.321 | 0.2603 | 1.52E-08 |
| atomic_photochem_gj1132_kzzprofile | 3 -> 3.0034 | info=0 certified | 7.900 | 1.7640 | 6.156 | 0.2684 | 1.08E-08 |
| atomic_photochem_gj1132_kzzprofile | 5 -> 5.0057 | info=0 certified | 7.900 | 2.0746 | 6.933 | 0.2805 | 2.21E-08 |
| atomic_photochem_gj1132_kzzprofile | 7 -> 7.0080 | info=0 certified | 7.900 | 2.1768 | 7.102 | 0.2871 | 2.15E-08 |
| atomic_photochem_gj1132_kzzprofile | 8 -> 8.0091 | info=0 certified | 7.900 | 2.2004 | 7.115 | 0.2896 | 1.44E-08 |
| atomic_photochem_gj1132_kzzprofile | 9 -> 9.0103 | info=0 certified | 7.900 | 2.2137 | 7.111 | 0.2916 | 3.93E-08 |
| atomic_photochem_gj1132_kzzprofile | 9.7 -> 9.7111 | info=0 certified | 7.900 | 2.2206 | 7.101 | 0.2931 | 1.46E-08 |
| atomic_photochem_gj1132_kzzprofile | 10 -> 10.0114 | info=0 certified | 7.900 | 2.2224 | 7.094 | 0.2936 | 3.03E-08 |
| atomic_photochem_gj1132_kzzprofile | 12 -> 12.0137 | info=0 certified | 7.900 | 2.2289 | 7.037 | 0.2966 | 3.21E-08 |
| **atomic_photochem_gj1132_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.01_kzz1e9 | 2.13 | none | -- | -- | -- | -- | -- |
| atomic_scalar_gj1132x0.01_kzz1e9 | 9.7 | none | -- | -- | -- | -- | -- |
| **atomic_scalar_gj1132x0.01_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.10_kzz1e9 | 2.13 | info=0 certified | 6.820 | 0.0001 | 0.000 | 0.2815 | 4.63E-07 |
| atomic_scalar_gj1132x0.10_kzz1e9 | 9.7 | info=2 uncertified | -- | -- | -- | -- | 2.91E-01 |
| **atomic_scalar_gj1132x0.10_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.15_kzz1e9 | 2.13 | info=0 certified | 7.000 | 0.0003 | 0.001 | 0.2865 | 8.92E-08 |
| atomic_scalar_gj1132x0.15_kzz1e9 | 9.7 | info=0 certified | 7.000 | 0.4514 | 1.671 | 0.2517 | 1.11E-07 |
| **atomic_scalar_gj1132x0.15_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.20_kzz1e9 | 2.13 | info=0 certified | 7.130 | 0.0022 | 0.007 | 0.2785 | 9.81E-08 |
| atomic_scalar_gj1132x0.20_kzz1e9 | 9.7 | info=0 certified | 7.130 | 0.6294 | 2.285 | 0.2573 | 1.04E-07 |
| **atomic_scalar_gj1132x0.20_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.25_kzz1e9 | 2.13 | info=0 certified | 7.240 | 0.0700 | 0.279 | 0.2351 | 1.57E-07 |
| atomic_scalar_gj1132x0.25_kzz1e9 | 9.7 | info=0 certified | 7.230 | 0.7913 | 2.823 | 0.2618 | 1.12E-07 |
| **atomic_scalar_gj1132x0.25_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.30_kzz1e9 | 2.13 | info=0 certified | 7.320 | 0.1894 | 0.760 | 0.2341 | 3.79E-08 |
| atomic_scalar_gj1132x0.30_kzz1e9 | 9.7 | info=0 certified | 7.310 | 0.9395 | 3.301 | 0.2664 | 5.56E-08 |
| **atomic_scalar_gj1132x0.30_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_scalar_gj1132x0.33_kzz1e9 | 2.13 | info=0 certified | 7.360 | 0.2555 | 1.020 | 0.2351 | 3.25E-08 |
| atomic_scalar_gj1132x0.33_kzz1e9 | 9.7 | info=0 certified | 7.360 | 1.0231 | 3.566 | 0.2684 | 8.60E-08 |
| **atomic_scalar_gj1132x0.33_kzz1e9** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.01_kzzprofile | 9.7 | none | -- | -- | -- | -- | -- |
| **atomic_photochem_gj1132x0.01_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.10_kzzprofile | 9.7 | info=0 certified | 6.840 | 0.2674 | 1.001 | 0.2487 | 3.05E-07 |
| **atomic_photochem_gj1132x0.10_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.15_kzzprofile | 9.7 | info=2 uncertified | -- | -- | -- | -- | 1.18E+00 |
| **atomic_photochem_gj1132x0.15_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.20_kzzprofile | 9.7 | info=2 uncertified | -- | -- | -- | -- | 6.69E-01 |
| **atomic_photochem_gj1132x0.20_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.25_kzzprofile | 9.7 | info=1 uncertified | -- | -- | -- | -- | 1.56E+00 |
| **atomic_photochem_gj1132x0.25_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.30_kzzprofile | 9.7 | info=0 certified | 7.350 | 1.0097 | 3.542 | 0.2664 | 3.45E-08 |
| **atomic_photochem_gj1132x0.30_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| atomic_photochem_gj1132x0.33_kzzprofile | 9.7 | info=0 certified | 7.390 | 1.0979 | 3.821 | 0.2689 | 3.22E-08 |
| **atomic_photochem_gj1132x0.33_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| molecular_scalar_gj1132_wellmixed | 0.083 | info=1 uncertified | -- | -- | -- | -- | 2.03E-08 |
| molecular_scalar_gj1132_wellmixed | 0.55 | info=1 uncertified | -- | -- | -- | -- | 3.48E-08 |
| molecular_scalar_gj1132_wellmixed | 2.13 | running | -- | -- | -- | -- | -- |
| **molecular_scalar_gj1132_wellmixed** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |
| molecular_scalar_gj1132_kzz1e9 | 0.083 | info=0 uncertified | -- | -- | -- | -- | 6.95E-09 |
| molecular_scalar_gj1132_kzz1e9 | 0.55 | info=0 certified | 7.640 | 0.1591 | 0.631 | 0.2366 | 1.15E-08 |
| molecular_scalar_gj1132_kzz1e9 | 2.13 | info=0 certified | 7.910 | 1.5758 | 5.628 | 0.2628 | 9.40E-09 |
| molecular_scalar_gj1132_kzz1e9 | 9.7 | info=0 certified | 7.950 | 2.4212 | 7.683 | 0.2951 | 2.04E-08 |
| **molecular_scalar_gj1132_kzz1e9** | crossing He/H = 1.7300  (log-log chord between the rungs 0.55 and 2.13) |  |  |  |  |  |  |
| molecular_photochem_gj1132_kzzprofile | 2.09 | running | -- | -- | -- | -- | -- |
| molecular_photochem_gj1132_kzzprofile | 9 | running | -- | -- | -- | -- | -- |
| **molecular_photochem_gj1132_kzzprofile** | crossing He/H = no bracket (EW = 1.108 %A not spanned) |  |  |  |  |  |  |

## 8. How each model was reached

One row per case that has been run, read from the `REPRODUCE.md` its own run wrote: which state seeded it, how many outer passes the stationary solve took, what the solver and the certification said, what came out, and how long it took. The linked file carries the commands themselves.

| case | seed | outer passes | verdict | certification | log10 Mdot | EW [%A] | wall clock | record |
|---|---|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | cold start | 0 | info=0 | certified | 7.84 | 0.3797 | 2m17s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.083/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | cold start | 0 | info=0 | certified | 7.85 | 1.1011 | 3m23s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.40/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | cold start | 0 | info=0 | certified | 7.85 | 1.1324 | 3m54s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.42/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | cold start | 0 | info=0 | certified | 7.85 | 1.1626 | 3m52s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.44/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | cold start | 0 | info=0 | certified | 7.85 | 1.3118 | 4m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH1` | cold start | 0 | info=0 | certified | 7.86 | 1.6959 | 6m25s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH1/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH10` | cold start | 0 | info=0 | certified | 7.85 | 1.9347 | 17m21s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH10/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH100` | cold start | 0 | info=0 | certified | 7.86 | 1.9645 | 16m36s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH100/REPRODUCE.md) |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | cold start | 0 | info=0 | certified | 7.85 | 1.9197 | 23m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_wellmixed/HeH1000/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | cold start | 0 | info=0 | certified | 8.57 | 0.9091 | 0m29s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.042/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | cold start | 0 | info=0 | certified | 8.57 | 0.9858 | 0m29s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.046/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | cold start | 0 | info=0 | certified | 8.57 | 1.0609 | 0m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.050/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | cold start | 0 | info=0 | certified | 8.57 | 1.6271 | 0m33s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH0.083/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH1` | cold start | 0 | info=0 | certified | 8.55 | 4.5250 | 1m35s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH1/REPRODUCE.md) |
| `atomic_scalar_gj699_wellmixed/HeH1000` | cold start | 0 | info=0 | certified | 8.42 | 2.8725 | 2m38s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj699_wellmixed/HeH1000/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | cold start | 13 | info=0 | certified | 7.84 | 0.0305 | 39m53s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | cold start | 7 | info=0 | certified | 7.86 | 0.9549 | 23m00s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH2.6/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | cold start | 7 | info=0 | certified | 7.87 | 1.0798 | 22m01s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.0/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | cold start | 7 | info=0 | certified | 7.87 | 1.2175 | 22m40s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.5/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | cold start | 6 | info=0 | certified | 7.87 | 1.2698 | 20m07s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.7/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | cold start | 6 | info=0 | certified | 7.87 | 1.3168 | 26m49s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz0/HeH3.9/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | cold start | 7 | info=0 | certified | 7.87 | 0.9631 | 28m48s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH2.6/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | cold start | 7 | info=0 | certified | 7.87 | 1.0880 | 30m29s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.0/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | cold start | 7 | info=0 | certified | 7.87 | 1.1862 | 26m56s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.35/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | cold start | 6 | info=0 | certified | 7.87 | 1.2627 | 24m01s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.64/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | cold start | 6 | info=0 | certified | 7.87 | 1.3314 | 24m29s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e5/HeH3.93/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | cold start | 12 | info=0 | certified | 7.84 | 0.0591 | 31m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | cold start | 7 | info=0 | certified | 7.86 | 0.9465 | 21m32s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH2.4/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | cold start | 6 | info=0 | certified | 7.87 | 1.0804 | 19m25s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH2.8/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | cold start | 6 | info=0 | certified | 7.87 | 1.1949 | 20m07s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.19/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | cold start | 6 | info=0 | certified | 7.87 | 1.2669 | 19m54s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.46/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | cold start | 5 | info=0 | certified | 7.87 | 1.3388 | 18m06s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e6/HeH3.74/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | cold start | 6 | info=0 | certified | 7.85 | 0.1585 | 15m08s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | cold start | 3 | info=0 | certified | 7.87 | 1.1996 | 11m51s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH2.70/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | cold start | 3 | info=0 | certified | 7.87 | 1.2724 | 11m43s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH2.94/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | cold start | 5 | info=0 | certified | 7.87 | 1.3265 | 16m58s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e7/HeH3.18/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | cold start | 5 | info=0 | certified | 7.85 | 0.3089 | 13m07s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | cold start | 3 | info=0 | certified | 7.87 | 1.1513 | 10m56s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.05/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | cold start | 3 | info=0 | certified | 7.87 | 1.2206 | 11m11s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.23/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | cold start | 3 | info=0 | certified | 7.87 | 1.2826 | 11m17s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e8/HeH2.41/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | cold start | 3 | info=0 | certified | 7.85 | 0.4793 | 8m33s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | cold start | 3 | info=0 | certified | 7.87 | 1.0996 | 9m31s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.50/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | cold start | 3 | info=0 | certified | 7.87 | 1.1511 | 9m22s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.60/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | cold start | 3 | info=0 | certified | 7.87 | 1.1970 | 9m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH1.70/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | cold start | 3 | info=0 | certified | 7.87 | 1.3749 | 10m54s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | cold start | 1 | info=0 | certified | 7.87 | 1.8229 | 267m47s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH4.0/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | cold start | 2 | info=0 | certified | 7.86 | 2.0592 | 17m44s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | cold start | 2 | info=0 | certified | 7.86 | 0.6463 | 7m06s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | cold start | 2 | info=0 | certified | 7.86 | 1.0456 | 9m14s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.06/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | cold start | 2 | info=0 | certified | 7.86 | 1.1032 | 13m43s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.15/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | cold start | 2 | info=0 | certified | 7.86 | 1.1867 | 10m44s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e10/HeH1.29/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | cold start | 3 | info=0 | certified | 7.86 | 0.8065 | 13m50s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.55/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | cold start | 3 | info=0 | certified | 7.86 | 1.0302 | 17m03s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.795/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | cold start | 3 | info=0 | certified | 7.86 | 1.0866 | 14m00s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.865/REPRODUCE.md) |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | cold start | 3 | info=0 | certified | 7.86 | 1.1364 | 16m41s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132_kzz1e11/HeH0.93/REPRODUCE.md) |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | cold start | 3 | info=0 | certified | 7.87 | 1.3721 | 12m23s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09` | cold start | 3 | info=0 | certified | 7.91 | 1.4765 | 35m40s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH2.09/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH3` | cold start | 3 | info=0 | certified | 7.90 | 1.7640 | 38m19s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH3/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH5` | cold start | 2 | info=0 | certified | 7.90 | 2.0746 | 37m06s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH5/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH7` | cold start | 3 | info=0 | certified | 7.90 | 2.1768 | 99m00s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH7/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH8` | cold start | 3 | info=0 | certified | 7.90 | 2.2004 | 96m48s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH8/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH9` | cold start | 2 | info=0 | certified | 7.90 | 2.2137 | 91m19s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH9/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7` | cold start | 2 | info=0 | certified | 7.90 | 2.2206 | 91m35s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH9.7/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH10` | cold start | 2 | info=0 | certified | 7.90 | 2.2224 | 88m37s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH10/REPRODUCE.md) |
| `atomic_photochem_gj1132_kzzprofile/HeH12` | cold start | 3 | info=0 | certified | 7.90 | 2.2289 | 91m03s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132_kzzprofile/HeH12/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | cold start | 17 | info=0 | certified | 6.82 | 0.0001 | 16m12s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | cold start | 19 | info=0 | certified | 7.00 | 0.0003 | 251m40s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | cold start | 6 | info=0 | certified | 7.00 | 0.4514 | 90m28s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | cold start | 34 | info=0 | certified | 7.13 | 0.0022 | 16m35s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | cold start | 6 | info=0 | certified | 7.13 | 0.6294 | 59m02s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | cold start | 11 | info=0 | certified | 7.24 | 0.0700 | 18m13s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | cold start | 5 | info=0 | certified | 7.23 | 0.7913 | 51m48s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | cold start | 8 | info=0 | certified | 7.32 | 0.1894 | 7m11s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | cold start | 5 | info=0 | certified | 7.31 | 0.9395 | 35m37s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | cold start | 8 | info=0 | certified | 7.36 | 0.2555 | 6m11s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13/REPRODUCE.md) |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | cold start | 4 | info=0 | certified | 7.36 | 1.0231 | 25m39s at 8 thread(s) | [REPRODUCE.md](models/atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7/REPRODUCE.md) |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | cold start | 8 | info=0 | certified | 6.84 | 0.2674 | 27m10s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/REPRODUCE.md) |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | cold start | 17 | info=1 | not certified | -- | -- | 168m30s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7/REPRODUCE.md) |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | cold start | 4 | info=0 | certified | 7.35 | 1.0097 | 13m21s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7/REPRODUCE.md) |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | cold start | 4 | info=0 | certified | 7.39 | 1.0979 | 15m16s at 8 thread(s) | [REPRODUCE.md](models/atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7/REPRODUCE.md) |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | cold start | 24 | info=1 | not certified | -- | -- | 318m09s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH0.083/REPRODUCE.md) |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | cold start | 40 | info=1 | not certified | -- | -- | 506m01s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH0.55/REPRODUCE.md) |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | cold start | 0 | -- | -- | -- | -- | 0m00s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_wellmixed/HeH2.13/REPRODUCE.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | cold start | 40 | info=1 | not certified | -- | -- | 600m04s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.083/REPRODUCE.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | cold start | 4 | info=0 | certified | 7.64 | 0.1591 | 214m57s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH0.55/REPRODUCE.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | cold start | 5 | info=0 | certified | 7.91 | 1.5758 | 314m42s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH2.13/REPRODUCE.md) |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | cold start | 13 | info=0 | certified | 7.95 | 2.4212 | 90m32s at 8 thread(s) | [REPRODUCE.md](models/molecular_scalar_gj1132_kzz1e9/HeH9.7/REPRODUCE.md) |
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | seed from its own atomic rung | -- | -- | -- | -- | -- | -- | -- |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | seed from its own atomic rung | -- | -- | -- | -- | -- | -- | -- |

# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH7`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The rung has no `input.inp` of its own: `input_template.inp` beside `closure.json` is what each iteration starts from, and the closure driver writes that iteration's own `input.inp` into `kNN/` from it, adding the route's keys (`input_keys`) and the profile the chemistry step just produced.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b9419ba7cd8e519ca86896b29c382762` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 6; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/L8p0/output  tier4:archive/refresh_j96  HeH=8  target=7  dlog10=0.0580  candidates=1614  (an archived state)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/L8p0/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed
```

## The commands

```bash
# the seed of iteration 0, on the current cell centers
mkdir -p seed
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/L8p0/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed

# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/"omp_num_threads": *[0-9]*/"omp_num_threads": 6/' \
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/element_flux_closure.py . \
    --config closure_run.json \
    --phi0-H 4.7682069882E+06 --phi0-He 2.0433648550E+07 \
    --tol 0.05 --kmax 8 --seed /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH7/seed > closure_stdout.log 2>&1

# the transit spectrum of the iterate the rung converged on
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
cd $(ls -d k[0-9][0-9] | sort | tail -n 1)
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

### Iteration `k00`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.86E-10 | 3.06E-14 | 2.22E-08 | 3.68E-04 of 1.0E-05 (elemental transport He/H partition) | 272 | 67.38 |
| 2 | 0 | 1.04E-09 | 3.73E-14 | 1.87E-08 | 1.80E-04 of 1.0E-05 (elemental transport He/H partition) | 275 | 72.97 |
| 3 | 0 | 7.85E-10 | 1.02E-14 | 1.49E-08 | 8.80E-05 of 1.0E-05 (elemental transport He/H partition) | 277 | 64.46 |
| 4 | 0 | 9.52E-10 | 4.61E-14 | 3.19E-08 | 4.32E-05 of 1.0E-05 (elemental transport He/H partition) | 280 | 67.38 |
| 5 | 0 | 6.59E-10 | 2.10E-14 | 1.99E-08 | 2.14E-05 of 1.0E-05 (elemental transport He/H partition) | 283 | 56.26 |
| 6 | 0 | 8.34E-10 | 3.64E-14 | 1.47E-08 | 1.06E-05 of 1.0E-05 (elemental transport He/H partition) | 286 | 37.68 |
| 7 | 0 | 2.55E-09 | 6.04E-14 | 5.27E-08 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 289 | 42.30 |

- solver verdict: **info = 0** (last `||R||` = 5.268E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:40:04 to 2026-09-14 11:21:59

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 6.81E-10 | 2.55E-13 | 1.71E-08 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 289 | 130.54 |

- solver verdict: **info = 0** (last `||R||` = 1.708E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:40:04 to 2026-09-14 11:21:59

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 6.99E-10 | 1.40E-13 | 2.49E-08 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 289 | 97.94 |

- solver verdict: **info = 0** (last `||R||` = 2.491E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:40:04 to 2026-09-14 11:21:59

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.38E-09 | 1.11E-12 | 4.41E-08 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 289 | 114.84 |

- solver verdict: **info = 0** (last `||R||` = 4.553E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:40:04 to 2026-09-14 11:21:59

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.97E-10 | 5.54E-13 | 1.93E-08 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 289 | 86.70 |

- solver verdict: **info = 0** (last `||R||` = 1.931E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.2208** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.2400** %, FWHM = 0.2876 A (three-Gaussian fit)
- wall clock **41m55s** at `OMP_NUM_THREADS=6`, 2026-09-14 10:40:04 to 2026-09-14 11:21:59

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   8.8694149548E+06   6.8756538085E+07   9.1427572222E+06   7.1978117213E+07  2.98971E-02  4.47578E-02   0.5000 overlap   2.48031E-03  4.97230E-03  7.0079514E+00    7.910 e452e286a7119f94b7020660ab60f2d95566bf29f8b48344a704c873da010e78 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=6`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH9.7`

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
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C1p00/output  tier4:archive/refresh_j96  HeH=9.7  target=9.7  dlog10=0.0000  candidates=1614  (an archived state)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C1p00/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed
```

## The commands

```bash
# the seed of iteration 0, on the current cell centers
mkdir -p seed
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C1p00/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed

# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/"omp_num_threads": *[0-9]*/"omp_num_threads": 8/' \
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/element_flux_closure.py . \
    --config closure_run.json \
    --phi0-H 4.7682069882E+06 --phi0-He 2.0433648550E+07 \
    --tol 0.05 --kmax 8 --seed /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9.7/seed > closure_stdout.log 2>&1

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
| 1 | 0 | 8.76E-10 | 1.23E-14 | 3.13E-08 | 2.93E-04 of 1.0E-05 (elemental transport He/H partition) | 271 | 32.70 |
| 2 | 0 | 1.13E-09 | 2.84E-14 | 3.03E-08 | 1.44E-04 of 1.0E-05 (elemental transport He/H partition) | 273 | 34.53 |
| 3 | 0 | 2.43E-09 | 6.02E-14 | 5.99E-08 | 7.06E-05 of 1.0E-05 (elemental transport He/H partition) | 276 | 30.83 |
| 4 | 0 | 3.57E-10 | 2.55E-14 | 1.42E-08 | 3.48E-05 of 1.0E-05 (elemental transport He/H partition) | 279 | 30.56 |
| 5 | 0 | 7.71E-10 | 1.35E-14 | 1.97E-08 | 1.73E-05 of 1.0E-05 (elemental transport He/H partition) | 281 | 29.73 |
| 6 | 0 | 2.20E-09 | 8.93E-14 | 6.95E-08 | 8.61E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 23.18 |

- solver verdict: **info = 0** (last `||R||` = 6.952E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8`, 2026-09-14 11:23:44 to 2026-09-14 11:43:21

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.12E-09 | 1.26E-12 | 1.87E-08 | 8.61E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 61.73 |

- solver verdict: **info = 0** (last `||R||` = 1.873E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8`, 2026-09-14 11:23:44 to 2026-09-14 11:43:21

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.80E-10 | 6.58E-13 | 3.02E-08 | 8.61E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 58.02 |

- solver verdict: **info = 0** (last `||R||` = 3.022E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8`, 2026-09-14 11:23:44 to 2026-09-14 11:43:21

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.89E-10 | 1.82E-12 | 2.45E-08 | 8.61E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 39.97 |

- solver verdict: **info = 0** (last `||R||` = 2.446E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8`, 2026-09-14 11:23:44 to 2026-09-14 11:43:21

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.69E-09 | 1.28E-12 | 1.89E-08 | 8.61E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 12.86 |

- solver verdict: **info = 0** (last `||R||` = 1.891E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.2639** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.2286** %, FWHM = 0.2936 A (three-Gaussian fit)
- wall clock **19m37s** at `OMP_NUM_THREADS=8`, 2026-09-14 11:23:44 to 2026-09-14 11:43:21

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   6.5921711992E+06   7.0813419900E+07   6.7137575891E+06   7.4172075811E+07  1.81100E-02  4.52819E-02   0.5000 overlap   3.71629E-03  4.85064E-03  9.7110262E+00    7.910 6270e5f4a036cc14d1f8dbec0435e0a1b2457d8d0058e000b2b98f13b96a1e69 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH9`

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
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH9` |

## The seed

None: this case starts from the code's own initial condition (`Load IC? False`). The archive holds no state of its physics.

## The commands

```bash
# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/"omp_num_threads": *[0-9]*/"omp_num_threads": 6/' \
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/element_flux_closure.py . \
    --config closure_run.json \
    --phi0-H 4.7682069882E+06 --phi0-He 2.0433648550E+07 \
    --tol 0.05 --kmax 8 --resume > closure_stdout.log 2>&1

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
| 1 | 0 | 1.34E-09 | 1.20E-14 | 2.41E-08 | 2.91E-04 of 1.0E-05 (elemental transport He/H partition) | 271 | 49.04 |
| 2 | 0 | 6.57E-10 | 1.37E-14 | 2.09E-08 | 1.42E-04 of 1.0E-05 (elemental transport He/H partition) | 273 | 67.33 |
| 3 | 0 | 5.47E-10 | 3.42E-14 | 1.76E-08 | 6.98E-05 of 1.0E-05 (elemental transport He/H partition) | 276 | 49.84 |
| 4 | 0 | 9.45E-10 | 1.10E-14 | 1.48E-08 | 3.43E-05 of 1.0E-05 (elemental transport He/H partition) | 279 | 49.18 |
| 5 | 0 | 5.99E-10 | 1.48E-14 | 1.54E-08 | 1.70E-05 of 1.0E-05 (elemental transport He/H partition) | 281 | 47.05 |
| 6 | 0 | 2.75E-09 | 1.27E-13 | 6.32E-08 | 8.43E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 32.04 |

- solver verdict: **info = 0** (last `||R||` = 6.324E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:57:09 to 2026-09-14 11:10:55

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.10E-09 | 2.13E-12 | 3.35E-08 | 8.43E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 131.64 |

- solver verdict: **info = 0** (last `||R||` = 3.347E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:57:09 to 2026-09-14 11:10:55

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.18E-09 | 1.06E-12 | 1.64E-08 | 8.43E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 128.47 |

- solver verdict: **info = 0** (last `||R||` = 1.641E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:57:09 to 2026-09-14 11:10:55

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.70E-10 | 4.95E-13 | 2.19E-08 | 8.43E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 109.89 |

- solver verdict: **info = 0** (last `||R||` = 2.191E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:57:09 to 2026-09-14 11:10:55

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.23E-10 | 1.48E-12 | 2.70E-08 | 8.43E-06 of 1.0E-05 (elemental transport He/H partition) | 284 | 85.91 |

- solver verdict: **info = 0** (last `||R||` = 1.788E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.2585** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.2444** %, FWHM = 0.2921 A (three-Gaussian fit)
- wall clock **13m46s** at `OMP_NUM_THREADS=6`, 2026-09-14 10:57:09 to 2026-09-14 11:10:55

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   7.0538701546E+06   7.0386331141E+07   7.2062201939E+06   7.3716527108E+07  2.11415E-02  4.51757E-02   0.5000 overlap   3.59807E-03  4.87359E-03  9.0102279E+00    7.910 37cdcb6fcbdbeea9295d0bef8bd41ec57103629112515729cf6f918485dc7ed9 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=6`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

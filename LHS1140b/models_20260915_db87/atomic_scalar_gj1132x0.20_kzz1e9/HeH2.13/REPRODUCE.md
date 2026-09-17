# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132x0.20`
- mixing: binary H/He element diffusion with He_Kzz = 1e9 cm^2/s (`kzz1e9`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## Why this run was made the way it was

The first solve of this re-run was taken at the default pseudo-time start (`EXHALE_PTC_DTAU0=1.0`) and refused: forty outer passes with the gated species row already at 1e-12 and the hydrodynamic ENERGY row of the outermost physical cell above its tolerance (item L4f), which is the ramp collapse of item L4h. It was stopped, and its logs are in `models/.stopped/atomic_scalar_gj1132x0.20_kzz1e9_HeH2.13_dtau0-1_20260916055251/`. This run is what `run_case.sh` prescribes for that situation (item L14): the same mapped seed at `EXHALE_PTC_DTAU0=1.0e8`.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `db87b88d1ce53facf1d61084fa535ca5` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.20_kzz1e9_HeH2.13` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/output  tier0:models_20260914_preL21/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13  HeH=2.13  target=2.13  dlog10=0.0000  candidates=97  (this same case's own most recent certified state)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the advection-corrected profiles: the solved state becomes the IC and one
# time step is taken at a CFL small enough to leave it where it is
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
    -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > pp.log 2>&1
mv -f input.inp.solved input.inp

# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 3.35E-08 | 2.79E-14 | 2.46E-07 | 2.17E-04 of 1.0E-05 (elemental transport He/H partition) | 217 | 235.36 |
| 2 | 0 | 1.75E-08 | 1.29E-14 | 1.38E-07 | 3.26E-04 of 1.0E-05 (elemental transport He/H partition) | 290 | 35.43 |
| 3 | 0 | 5.38E-09 | 1.11E-14 | 5.17E-08 | 5.26E-04 of 1.0E-05 (elemental transport He/H partition) | 290 | 35.77 |
| 4 | 0 | 1.93E-08 | 2.33E-14 | 1.63E-07 | 5.75E-04 of 1.0E-05 (elemental transport He/H partition) | 290 | 35.36 |
| 5 | 0 | 7.56E-09 | 1.03E-14 | 6.08E-08 | 5.99E-04 of 1.0E-05 (elemental transport He/H partition) | 290 | 32.01 |
| 6 | 0 | 6.73E-09 | 1.11E-14 | 5.70E-08 | 6.02E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 33.94 |
| 7 | 0 | 1.29E-08 | 2.41E-14 | 9.71E-08 | 5.90E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 26.56 |
| 8 | 0 | 5.56E-09 | 1.13E-14 | 6.07E-08 | 5.65E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 33.43 |
| 9 | 0 | 7.04E-09 | 8.71E-15 | 7.59E-08 | 5.31E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 27.88 |
| 10 | 0 | 4.07E-08 | 3.96E-14 | 2.86E-07 | 4.92E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 22.71 |
| 11 | 0 | 3.47E-09 | 9.21E-15 | 6.41E-08 | 4.50E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 28.57 |
| 12 | 0 | 6.08E-09 | 9.36E-15 | 7.81E-08 | 4.06E-04 of 1.0E-05 (elemental transport He/H partition) | 291 | 28.70 |
| 13 | 0 | 9.81E-09 | 2.23E-14 | 6.47E-08 | 3.64E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 22.83 |
| 14 | 0 | 4.22E-09 | 8.97E-15 | 7.28E-08 | 3.22E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 27.75 |
| 15 | 0 | 1.25E-08 | 9.19E-15 | 8.92E-08 | 2.84E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 28.19 |
| 16 | 0 | 1.81E-08 | 1.13E-14 | 1.42E-07 | 2.48E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 22.10 |
| 17 | 0 | 4.00E-08 | 3.77E-14 | 2.89E-07 | 2.15E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 22.57 |
| 18 | 0 | 6.36E-09 | 1.34E-14 | 7.56E-08 | 1.85E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 21.35 |
| 19 | 0 | 3.59E-08 | 2.33E-14 | 2.70E-07 | 1.59E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 20.28 |
| 20 | 0 | 8.48E-09 | 1.07E-14 | 7.48E-08 | 1.36E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 14.46 |
| 21 | 0 | 1.70E-08 | 2.01E-14 | 1.12E-07 | 1.15E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 15.03 |
| 22 | 0 | 9.12E-09 | 9.66E-15 | 6.88E-08 | 9.72E-05 of 1.0E-05 (elemental transport He/H partition) | 292 | 21.14 |
| 23 | 0 | 1.02E-08 | 1.12E-14 | 8.14E-08 | 8.18E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 21.92 |
| 24 | 0 | 3.10E-09 | 9.74E-15 | 8.51E-08 | 6.86E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 20.88 |
| 25 | 0 | 1.47E-08 | 9.23E-15 | 9.30E-08 | 5.73E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 20.97 |
| 26 | 0 | 3.51E-08 | 2.61E-14 | 2.68E-07 | 4.77E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 13.06 |
| 27 | 0 | 1.36E-08 | 1.33E-14 | 1.10E-07 | 3.96E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 18.65 |
| 28 | 0 | 7.05E-09 | 1.27E-14 | 6.54E-08 | 3.27E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 13.01 |
| 29 | 0 | 4.13E-09 | 1.02E-14 | 7.02E-08 | 2.70E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 18.09 |
| 30 | 0 | 3.06E-09 | 8.42E-15 | 6.17E-08 | 2.21E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 14.42 |
| 31 | 0 | 3.83E-09 | 9.12E-15 | 1.07E-07 | 1.81E-05 of 1.0E-05 (elemental transport He/H partition) | 294 | 14.24 |
| 32 | 0 | 3.71E-08 | 3.08E-14 | 2.96E-07 | 1.48E-05 of 1.0E-05 (elemental transport He/H partition) | 294 | 6.95 |
| 33 | 0 | 1.75E-08 | 8.53E-15 | 1.30E-07 | 1.21E-05 of 1.0E-05 (elemental transport He/H partition) | 294 | 7.10 |
| 34 | 0 | 1.47E-08 | 9.10E-15 | 9.66E-08 | 9.81E-06 of 1.0E-05 (elemental transport He/H partition) | 294 | 7.02 |

- solver verdict: **info = 0** (last `||R||` = 9.813E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.13**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0022** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0075** %, FWHM = 0.2785 A (three-Gaussian fit)
- wall clock **16m35s** at `OMP_NUM_THREADS=8` on `lart4`, 2026-09-16 05:53:42 to 2026-09-16 06:10:17

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.

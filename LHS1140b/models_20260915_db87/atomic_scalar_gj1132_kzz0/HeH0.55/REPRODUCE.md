# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132_kzz0/HeH0.55`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 0 cm^2/s (`kzz0`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `db87b88d1ce53facf1d61084fa535ca5` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132_kzz0_HeH0.55` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132_kzz0/HeH0.55/output  tier0:preL21/atomic_scalar_gj1132_kzz0/HeH0.55  HeH=0.55  target=0.55  dlog10=0.0000  candidates=256  (this same case's own state of 2026-09-14, before the base boundary of plan item L21 was corrected)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own state of 2026-09-14, in `models_20260914_preL21/`, which the re-run on the corrected base boundary of plan item L21 continues from, `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132_kzz0/HeH0.55/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_scalar_gj1132_kzz0/HeH0.55/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
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
| 1 | 0 | 7.39E-10 | 1.07E-13 | 1.61E-08 | 4.80E-05 of 1.0E-05 (elemental transport He/H partition) | 217 | 260.35 |
| 2 | 0 | 5.39E-10 | 5.99E-14 | 1.62E-08 | 3.45E-04 of 1.0E-05 (elemental transport He/H partition) | 283 | 175.00 |
| 3 | 0 | 1.07E-09 | 7.90E-14 | 1.54E-08 | 4.85E-04 of 1.0E-05 (elemental transport He/H partition) | 284 | 172.33 |
| 4 | 0 | 8.29E-10 | 4.25E-14 | 1.75E-08 | 4.83E-04 of 1.0E-05 (elemental transport He/H partition) | 285 | 169.91 |
| 5 | 0 | 8.19E-10 | 4.10E-14 | 1.94E-08 | 3.97E-04 of 1.0E-05 (elemental transport He/H partition) | 286 | 167.35 |
| 6 | 0 | 1.00E-09 | 1.74E-12 | 1.83E-08 | 2.85E-04 of 1.0E-05 (elemental transport He/H partition) | 286 | 157.91 |
| 7 | 0 | 1.14E-09 | 7.59E-13 | 1.81E-08 | 1.86E-04 of 1.0E-05 (elemental transport He/H partition) | 287 | 154.03 |
| 8 | 0 | 8.21E-10 | 3.07E-13 | 1.54E-08 | 1.13E-04 of 1.0E-05 (elemental transport He/H partition) | 287 | 122.59 |
| 9 | 0 | 5.90E-10 | 1.11E-13 | 1.15E-08 | 6.58E-05 of 1.0E-05 (elemental transport He/H partition) | 287 | 179.42 |
| 10 | 0 | 1.05E-09 | 2.14E-12 | 1.93E-08 | 3.70E-05 of 1.0E-05 (elemental transport He/H partition) | 287 | 178.73 |
| 11 | 0 | 1.07E-09 | 6.60E-13 | 2.24E-08 | 2.04E-05 of 1.0E-05 (elemental transport He/H partition) | 288 | 185.28 |
| 12 | 0 | 9.14E-10 | 2.30E-13 | 1.43E-08 | 1.11E-05 of 1.0E-05 (elemental transport He/H partition) | 288 | 197.57 |
| 13 | 0 | 1.20E-09 | 1.60E-12 | 1.71E-08 | 5.95E-06 of 1.0E-05 (elemental transport He/H partition) | 288 | 225.86 |

- solver verdict: **info = 0** (last `||R||` = 1.707E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.84**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0305** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.1267** %, FWHM = 0.2270 A (three-Gaussian fit)
- wall clock **39m53s** at `OMP_NUM_THREADS=8`, 2026-09-15 20:50:53 to 2026-09-15 21:30:46

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.

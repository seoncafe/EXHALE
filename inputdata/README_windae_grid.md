# `windae_grid/` — Wind-AE starting-solution grid

The directory `inputdata/windae_grid/` holds the **Wind-AE solution grid** from
Broome et al. (2025a): a database of ~1000 converged Wind-AE wind solutions
spanning a wide range of planet mass, radius, and irradiation. EXHALE uses it to
**bootstrap the Wind-AE warm-start initial condition** (`IC mode: windae`): the
routine `pick_nearest_seed` (`src/modules/wind_ae/wae_continuation.f90`) scans
`windae_grid/manifest.csv`, picks the converged solution nearest to the target
planet in (M_p, R_p, F_tot), and ramps it to the requested planet as the seed.

## Not stored in the repository

The grid is **~552 MB**, so it is **git-ignored** (`.gitignore`:
`inputdata/windae_grid/`) and must be downloaded separately. EXHALE runs that do
**not** use `IC mode: windae` (i.e. the default cold IC) do not need it at all.

## How to obtain it

Download the grid from either source:

- **Google Drive:**
  <https://drive.google.com/file/d/1Ptz3YUMD3mNktR7snUUIRXfU5WbTcVb0/view?usp=sharing>
- **wind-ae GitHub repository:**
  <https://github.com/mibroome/wind-ae/tree/main/Notebooks/Broome%20et%20al.%20(2025a)/data>

## Where to put it

Extract/copy the contents so that the layout under
`EXHALE/inputdata/windae_grid/` is:

```
windae_grid/
├── manifest.csv          # one "path,Mp[g],Rp[cm],Ftot[erg/s/cm^2]" row per solution
├── High_flux/            # high-irradiation solutions (*.csv windsoln files)
├── Low_flux/             # low-irradiation solutions
└── starting_solutions/   # additional seed solutions
```

`manifest.csv` lists each solution file with its planet parameters; the paths in
it are **relative to the EXHALE root** (e.g.
`inputdata/windae_grid/High_flux/hi_10.08Me_1.85Re.csv`). If you place the grid
elsewhere, regenerate the manifest with
`src/utils/make_windae_manifest.sh` so the paths match.

If the manifest is missing, `pick_nearest_seed` falls back (with a warning) to
the single seed `inputdata/windae_seed.csv`.

## Reference

Broome et al. (2025a), the Wind-AE atmospheric-escape solver and its
starting-solution grid (`mibroome/wind-ae`).

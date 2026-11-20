# Photochem, as EXHALE needs it

The Photochem route of the lower-atmosphere profile handoff
(`src/utils/photochem_to_lower_profile.py`, driven by
`src/utils/element_flux_closure.py`) runs on upstream Photochem `v0.9.3`
with the changes of `src/utils/photochem_exhale.patch`. Photochem and its
data tables are third-party and are not part of this repository.

| What | Where | Version | Directory |
|---|---|---|---|
| Photochem (with Clima and Equilibrate in its tree) | <https://github.com/Nicholaswogan/photochem> | tag `v0.9.3` (`ba3542c`) | `photochem/` |
| its data | <https://github.com/Nicholaswogan/photochem_clima_data> | tag `v0.3.4` | installed Python package |

Cite, as upstream asks, Wogan et al. (2023), doi:10.3847/PSJ/aced83, or
Wogan et al. (2024), doi:10.3847/2041-8213/ad2616.

All commands below run from the top of the EXHALE repository.

## 1. Source

```bash
src/utils/setup_photochem.sh                 # -> photochem/
src/utils/setup_photochem.sh <dir>           # -> <dir>
src/utils/setup_photochem.sh --url <url> <dir>
```

The script clones Photochem, checks out `ba3542c` (`v0.9.3`), applies
`src/utils/photochem_exhale.patch` and verifies each change in the patched
files (one line `ok` or `MISSING` each; nonzero exit on any miss). A tree at
another commit, or one that is neither pristine `v0.9.3` nor patched, is
refused. Run on a patched tree it verifies and changes nothing. It does not
build.

## 2. A Python environment for Photochem

Photochem, its data and the build tools go into one Python environment,
which is then the Python that runs the adapter. A plain virtual environment
on Python 3.11 does not see the user site-packages
(`~/.local/lib/python3.11/site-packages`), so a `photochem` or NumPy
installed there cannot shadow the ones built here (with
`--system-site-packages` the user site is searched again, also by the NumPy
probe of the CMake build):

```bash
python3 -m venv photochem_env
. photochem_env/bin/activate
python -m pip install "cmake>=3.20,<4" scikit-build-core ninja fypp cython \
    numpy scipy pyyaml numba h5py astropy threadpoolctl requests \
    pytest astroquery
```

(`pytest` and `astroquery` are needed only by Photochem's tests, section 5.)

- CMake must be 3.x (the `cmake` package from pip provides it inside the
  environment): CMake 4.x stops in the CVODE 5.7.0 subproject.
- With NumPy 2 the current Cython works. If the environment has NumPy 1.26,
  use `"cython>=3.0,<3.1"`: Cython 3.3 generates NumPy 2 C-API calls that
  the NumPy 1.26 headers do not have.
- A Fortran compiler accepted by upstream is needed (gfortran 14 or newer)
  and a C compiler; name them with `FC` and `CC` when the first `gfortran`
  on the `PATH` is older.

## 3. Data

`import photochem` imports `photochem_clima_data` (for `DATA_DIR`), so the
data package must be installed in the same environment. Use `v0.3.4`:
Photochem `v0.9.3` refuses a mechanism whose gas thermodynamic polynomials
are discontinuous at a temperature join, and the reaction mechanism of
`v0.3.1` is.

```bash
git clone https://github.com/Nicholaswogan/photochem_clima_data
git -C photochem_clima_data checkout v0.3.4
python -m pip install ./photochem_clima_data
```

The CMake build of Photochem downloads `photochem_clima_data_v0.3.4.tar.gz`
into the source tree and unpacks its `data/` there for Photochem's own
tests; that does not install the Python package.

## 4. Build and install

```bash
FC=gfortran CC=gcc python -m pip wheel ./photochem \
    --no-build-isolation --no-deps -w wheels \
    -Cbuild-dir=build_photochem -Ccmake.build-type=Release \
    -Cbuild.tool-args=-j8
python -m pip install --no-deps wheels/photochem-0.9.3-*.whl
```

The build directory `build_photochem/` is created inside `photochem/`.

`pyproject.toml` sets the Python bindings on, the Fortran executables off,
OpenMP on and position-independent code; `-j8` limits the parallel
compilation. The shared libraries record the BLAS/LAPACK they were linked
against; check them with `ldd` on `photochem/_photochem*.so` inside the
environment's site-packages.

## 5. Verify

1. `src/utils/setup_photochem.sh` on the tree: every line `ok`.
2. The environment imports the library built from that tree and the data
   package installed above:
   ```bash
   cd /tmp && python -c "import photochem, photochem_clima_data; print(photochem.__version__, photochem.__file__, photochem_clima_data.__file__)"
   ```
3. The tests of the patched tree (the CMake build has unpacked `data/`):
   ```bash
   python -m pytest -q photochem/tests/photochem/test_python.py \
       -k "equilibrate or gas_giant"
   ```
   The test module imports `photochem.utils.stars`, which needs
   `astroquery`.

## 6. Running the adapter

The adapter writes into the run directory it is given (the first
argument): `lower_atmosphere_profile.dat`, a provenance `base.inp`, with
`--climate` and a conducted heat flux also `conducted_heat_column.txt`, the
face record `wind_face_radiation.txt`, and the generated mechanism and
climate files in `photochem_work/`. A run directory must not contain a
`photochem_work/` written with other data: the adapter reuses the generated
mechanism files it finds there.

- `python src/utils/photochem_to_lower_profile.py <run_dir> ...` runs it
  on one OpenMP thread: Photochem and Clima set the OpenMP limit to one when
  they are imported.
- `OMP_NUM_THREADS=N python src/utils/photochem_parallel.py <run_dir> ...`
  runs the same adapter, with the same arguments, after raising that limit
  to `N`. The threads are used by the photolysis wavelength loop, Clima's
  radiative transfer and the exchange sum of the column energy balance.
- The closure driver (`src/utils/element_flux_closure.py`) runs whatever
  `"adapter"` and `"python"` its configuration names, with the environment
  it was started in.

A complete example that uses only files of this repository, the LHS 1140 b
column at deep He/H = 1e5 with the climate solve and a conducted heat flux:

```bash
mkdir -p runs/lhs1140b_column
OMP_NUM_THREADS=12 python src/utils/photochem_parallel.py runs/lhs1140b_column \
    --mp 0.0176220 --r-ref 0.157692 --p-ref 1.0 --p-match 1.0e-6 \
    --climate --climate-p-deep 20.0 --boa-pressure-factor 1.0 \
    --stellar-flux LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
    --flux-at-planet --wavelength-unit A \
    --atoms H,He,N,O,C --abundances He=100000 --kzz-const 1.0e9 \
    --column-energy-treatment nlte --column-band-transfer cooling_to_space \
    --top-pressure-passes 12 --face-geometry wind --wind-dilution 1.0 \
    --trial-flux-H 655.230196665 --trial-flux-He 78326545.6645 \
    --trial-flux-H2 162746.103945 --iteration 1 --conducted-heat-flux 1.81195745
```

The trial fluxes and the conducted heat flux are the closure's values of
one LHS 1140 b iteration; the closure driver supplies them in its own runs.
The run writes an 836-level profile from 16.7 bar to the matching level at
1e-6 bar, and a 1590-node `conducted_heat_column.txt`. Measured for this
example at 12 threads on a 36-core Xeon host: 3 min 21 s of wall time and
2.5 GB of peak resident memory, with three chemistry solves and three
column energy balances. Other planets, grids and thread counts take
different times.

## 7. What the patch changes

`src/utils/photochem_exhale.patch` changes nine files of the upstream tree
and adds one test file.

- **Equilibrate mass balance** (`src/equilibrate/equilibrate_cea.f90`): the
  elemental mass-balance test of both update paths compares every requested
  element with its own relative tolerance (`mass_tol` times its abundance);
  the convergence flag of both abundance updates is set to false on entry.
- **What Equilibrate reports as converged** (`equilibrate_cea.f90`): a solve
  returns `True` only if the iteration met its convergence tests and
  - no condensate stage was ended by the rule that counts a stage still
    holding a negative condensed amount after 30 iterations as converged;
  - the LU factorization of the Newton system met no exactly zero pivot;
  - in the final state no absent condensate whose thermodynamic data cover
    the temperature passes the solver's own inclusion test or has a
    non-finite affinity, no included condensed amount is negative and every
    amount is finite.

  Otherwise it returns `False`, and `ChemEquiAnalysis.nonconvergence_reason`
  says why (empty after a converged solve; reset by every solve).
  `ChemEquiAnalysis.smallest_scaled_pivot` reports the smallest LU pivot
  relative to the largest entry of its row, as a diagnostic value only. The
  iteration itself is the upstream one, so these checks refuse states; they
  do not make a non-converging case converge.
- **The returned state** (`src/equilibrate/equilibrate.f90`): with
  `use_prev_guess`, condensed amounts of the previous solution that the new
  solve does not include are returned as zero, and the total mole and mass
  fractions are normalized over the gas and the included condensates. The
  gas quantities (`molfracs_species_gas`, `molfracs_atoms_gas`, `mubar`) are
  those of the iteration.
- **Gas giants** (`photochem/extensions/gasgiants.py`):
  `EvoAtmosphereGasGiant` and `GasGiantData` take `equilibrium_mass_tol`
  (default 1e-12); `composition_at_metallicity` checks at every level, with
  or without a condensed phase, that the returned total elemental inventory
  equals the level's request within `max(1e-8, 100*mass_tol)`, and raises if
  no attempt closes.
- **Clima** (`src/clima/clima_adiabat.f90`): `make_profile_bg_gas` brackets
  and bisects the background partial pressure; `simple_solver` solves for
  `log10(T_surf - T_trop)` with scaled residuals and rejects a returned point
  whose normalized residual exceeds 1e-7.
- **Bindings and tests**: the two attributes above in
  `bindings/equilibrate/`; tests in `tests/photochem/test_python.py`,
  `tests/equilibrate/test_returned_state_guards.py`,
  `tests/equilibrate/test_python.py` and
  `tests/equilibrate/test_equilibrate.f90`.

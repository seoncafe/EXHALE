# Photochem, as EXHALE needs it

Date: 2026-08-29

The Photochem arm of EXHALE's lower-atmosphere profile handoff
(`src/utils/photochem_to_lower_profile.py`) and the elemental-flux closure
that drives it (`src/utils/element_flux_closure.py`) call Photochem as a
library. What they need is **not** the release as published and **not**
`conda install -c conda-forge photochem`: those give a different elemental
closure at the handoff level, so results are not interchangeable with them.
What EXHALE runs is upstream Photochem `v0.9.0` with the four changes of
section 4 and nothing else. Which build wrote a given handoff is recorded in
the `# source_version photochem <version>` header line of every
`lower_atmosphere_profile.dat`, and by nothing else.

This file is what to download, what has to be installed to build it, how to
build it, what differs from upstream, and how to check the result. The
physics of the two solver corrections, why they were needed, and the
measurements that closed them are in
`docs/photochem_solver_modification_investigation.md` and
`docs/photochem_solver_modification_implementation.md`.


## 1. What to download from GitHub

Two repositories, neither of which is part of the EXHALE repository — both
are excluded in `.gitignore` (`/photochem/`, `/photochem_clima_data/`), so do
not expect to clone EXHALE and find either. Both live at the top of the
EXHALE working tree.

| | repository | version | goes to |
|---|---|---|---|
| the code | <https://github.com/Nicholaswogan/photochem> | tag `v0.9.0`, commit `e1e872528b61e8dd1db891b84869738e219b4f97` (2026-08-07) | `photochem/` |
| its data | <https://github.com/Nicholaswogan/photochem_clima_data> | tag `v0.3.1` | `photochem_clima_data/` |

### 1.1 The code

```bash
src/utils/setup_photochem.sh              # -> EXHALE_v1.00/photochem
src/utils/setup_photochem.sh <dir>        # -> <dir>
```

The script clones the URL above, checks out the pinned commit, applies
`src/utils/photochem_exhale.patch` — which is the whole of section 4 in one
file — and then verifies eleven things in the files themselves rather than
inferring them from `git apply` having exited 0. It refuses, rather than
half-applying, a tree at any other commit: the patch is cut against that
commit. Re-running it on an already-patched tree verifies and changes
nothing. It produces the source only; it does not build.

By hand instead:

```bash
git clone https://github.com/Nicholaswogan/photochem <dir>
git -C <dir> checkout e1e872528b61e8dd1db891b84869738e219b4f97
git -C <dir> apply --whitespace=nowarn src/utils/photochem_exhale.patch
```

The version is pinned for two reasons. `v0.9.0` is a release and is where
`origin/main` pointed when this was written; `origin/dev` is not a release
and still carries the solver behavior section 4 corrects. And the patch is a
context diff against that tree — at another commit it would apply
approximately or not at all.

Nothing else has to be fetched by hand for the build: Clima `v0.7.5` and
Equilibrate `v0.2.2` are named in `photochem/src/dependencies/CMakeLists.txt`
and downloaded by CPM during configure, which also applies the two dependency
patches sitting beside that file. The build therefore needs network access.

### 1.2 The data

Photochem imports `photochem_clima_data` at run time (`DATA_DIR`: cross
sections, reaction mechanisms, k-distributions, Rayleigh and CIA tables), and
it is a separate repository. Version `0.3.1` is not a choice — it is the
value of `PHOTOCHEM_CLIMA_DATA_VERSION` in Photochem 0.9.0's own top-level
`CMakeLists.txt`, so it is the data release that version of the code was
built and tested against.

```bash
git clone https://github.com/Nicholaswogan/photochem_clima_data
git -C photochem_clima_data checkout v0.3.1
python3 -m pip install --user -e ./photochem_clima_data
```

An editable install is what keeps the tables (38 MB measured here) in one
place instead of copying them into `site-packages`. A plain
`pip install ./photochem_clima_data` works too and costs the extra copy.


## 2. Required packages

Split by who needs them. Nothing here is EXHALE-specific — it is what
Photochem 0.9.0 needs to compile on a Linux machine.

### 2.1 Compiler and system tools

| | requirement | where the requirement comes from |
|---|---|---|
| Fortran | **GNU Fortran >= 14.0** | Photochem's own `CMakeLists.txt` raises `FATAL_ERROR` below it |
| C | a GCC matching the Fortran compiler | `FortranCInterface_VERIFY()` links the two at configure time |
| CMake | **>= 3.14**; with **CMake >= 4**, see the policy trap in section 3.2 | `cmake_minimum_required(VERSION "3.14")` |
| pkg-config | any | not Photochem's own CMake — one of the packages CPM downloads; configure stops without it |
| HDF5 | a build with **Fortran bindings** | likewise a downloaded dependency; `HDF5_ROOT` is how to point at it |
| BLAS / LAPACK | any (reference, OpenBLAS, MKL) | the linear algebra under the Fortran dependencies |

The last three are requirements of the dependency tree CPM fetches (futils,
fortran-yaml-c, CVODE 5.7, Clima, Differentia, Equilibrate), not of
Photochem's own build files, which is why they do not appear in
`photochem/CMakeLists.txt`. On a machine that has none of the six, a conda
environment supplies all of them — `gfortran_linux-64>=14`, `cmake`,
`pkg-config`, `hdf5` and a BLAS implementation are all conda-forge packages.
Where the system toolchain is already new enough, none of it is needed.

### 2.2 Python, to build

```bash
python3 -m pip install --user scikit-build ninja cython fypp
```

Python **3.11** is what the build here was made and verified on.

### 2.3 Python, at run time

Photochem's `install_requires` — install these **before** the build, because
the build is run with `--no-deps` (section 3.2):

```bash
python3 -m pip install --user numpy scipy pyyaml h5py numba astropy threadpoolctl requests
```

`photochem_clima_data` (section 1.2) is a run-time import as well, and is
not in `install_requires`.


## 3. Build and install

### 3.1 The command

One `pip install` of the patched source tree. `NP` is read from the same
interpreter Photochem is being installed into, so that the NumPy headers CMake
compiles against and the NumPy that will import the result are the same one.

```bash
PY=python3                                   # the interpreter to install into
NP=$($PY -c 'import numpy; print(numpy.get_include())')
BLAS=<the BLAS/LAPACK link line>             # e.g. -lopenblas, or an MKL chain

env FC=<gfortran >= 14> \
    CMAKE_POLICY_VERSION_MINIMUM=3.5 \
    HDF5_ROOT=<hdf5 prefix> \
    CMAKE_ARGS="-DNumPy_INCLUDE_DIR=$NP -DHDF5_ROOT=<hdf5 prefix> \
                -DBLAS_LIBRARIES=$BLAS -DLAPACK_LIBRARIES=$BLAS" \
    $PY -m pip install --user --no-deps --force-reinstall \
        --no-build-isolation ./photochem
```

Drop `--user` where the interpreter's own `site-packages` is writable. The
build takes about seven minutes from cold, most of it CVODE and the Fortran
dependencies.

### 3.2 Traps

Every one of these was hit while building this. They are properties of the
build, not of any one machine.

- **`--no-deps` is not optional.** Without it pip resolves
  `install_requires` and will happily upgrade NumPy underneath you. An
  extension compiled against the NumPy 1.x headers then aborts on import, and
  every other package installed against the old NumPy breaks with it. Install
  the run-time dependencies first (section 2.3), then build with `--no-deps`.

- **CMake 4 needs `CMAKE_POLICY_VERSION_MINIMUM=3.5` in the environment.**
  libyaml and CVODE, which arrive through the dependency chain, still declare
  a `cmake_minimum_required` that CMake 4 refuses. Passing the setting on the
  command line is not enough — it does not reach the sub-configure of CVODE's
  `try_compile`, and only the environment variable does.

- **the NumPy headers.** scikit-build's `FindNumPy` will take a system
  `/usr/include/numpy` ahead of the interpreter's own NumPy, and such a path
  is often a broken symlink into an old Python tree. Pass
  `-DNumPy_INCLUDE_DIR` explicitly, computed from the interpreter being
  installed into.

- **which Fortran compiler.** Putting a directory first on `PATH` does not
  select its compiler: CMake's Fortran search tries `f95` before `gfortran`,
  so a system `/usr/bin/f95` older than 14.0 ends configure at
  `Photochem will only work with gfortran >= 14.0.0` even with a new enough
  `gfortran` ahead of it. Name it with `FC` or `-DCMAKE_Fortran_COMPILER`.

- **BLAS/LAPACK autodetection.** Where CMake's `FindBLAS`/`FindLAPACK` come
  up empty or pick a library the Fortran compiler cannot link against, give
  `-DBLAS_LIBRARIES` and `-DLAPACK_LIBRARIES` the link line directly, as a
  `;`-separated list of absolute paths and flags. For a GNU-threaded MKL
  under a conda prefix `P`, that is

  ```
  P/lib/libmkl_gf_lp64.so;P/lib/libmkl_gnu_thread.so;P/lib/libmkl_core.so;P/lib/libgomp.so;-lpthread;-lm;-ldl
  ```

- **do not import `photochem` from inside the source tree.** There
  `import photochem` finds the pure-Python sources without their compiled
  extension and tells you nothing about what was built. Verify from another
  directory (section 5).

- **build in one run.** A `pip` killed and restarted immediately can leave a
  surviving child clone of a dependency racing the new one, and the result is
  a binary in which the patch was applied *after* the compile — installed,
  importable, and unrepaired. If a build is interrupted, let it finish dying
  before starting the next one, and verify the **installed** library rather
  than the source tree.


## 4. What differs from upstream

Two of the four changes are corrections to Photochem's dependencies, carried
as patch files and applied by CPM at configure time; two are in Photochem's
own Python layer. Paths in this section are relative to `photochem/`.

### 4.1 `src/dependencies/patches/equilibrate-element-relative-mass-closure.patch`

Equilibrate 0.2.2, both CEA update paths in `src/equilibrate_cea.f90`. The
elemental mass-balance test was one absolute threshold for every element,
`MAXVAL(b_0)*mass_tol`, scaled by the *most abundant* element and skipped
entirely for any element below `1e-6`. A trace element could therefore be
several orders of magnitude out of balance and the solve still report
`converged = True`. Each positive requested element is now tested against its
own relative tolerance,

```fortran
mval_mass_good = abs(b_0(i_atom))*self%mass_tol
```

and the `1e-6` exclusion is gone. On the LHS 1140 b column this moves the
deepest-level N/H departure from `1e-4`–`2e-4` to `1e-14`.

### 4.2 `src/dependencies/patches/clima-bounded-scaled-solvers.patch`

Clima 0.7.5, `src/adiabat/clima_adiabat.f90`. Two MINPACK `hybrd1` solves
are replaced by solves that stay in their physical domain and differentiate
at a step that resolves their own residual:

- `AdiabatClimate_make_profile_bg_gas` — the background pressure is positive
  and cannot exceed the requested total surface pressure. The unconstrained
  local solve is replaced by a bracketing scan over that interval in
  log pressure followed by bisection, so an isolated invalid thermodynamic
  state is stepped over rather than terminating the solve.
- `AdiabatClimate_simple_solver` — the surface-temperature solve. This is one
  routine, and `surface_temperature`, `surface_temperature_column` and
  `surface_temperature_bg_gas` are three thin wrappers that hand it different
  flux functions, so a change here reaches all three. The unknown is changed
  to the positive contrast `T_surf - T_trop`, so every trial satisfies
  `T_surf > T_trop`; the flux residual is divided by a scale built from the
  bolometric flux, so the mW/m² flux row and the temperature row of the
  coupled residual are comparable; the returned solution is rejected if its
  normalized residual exceeds `1e-7` instead of being accepted on MINPACK's
  exit code alone; and `hybrd1` is replaced by `hybrd` so that the
  forward-difference step can be set (`epsfcn = 1e-8`).

#### Why the difference step is set

`hybrd1` fixes MINPACK's `epsfcn` at zero, which makes the step
`sqrt(machine epsilon)*|x_j|`, 1.5e-8 relative. The unknowns here are
`log10(T_surf - T_trop)` and `log10(T_trop)`, so at the LHS 1140 b start
(`T_surf = 400 K`, `T_trop = 120 K`) that is a temperature step of 2e-5 K —
over which the normalized flux residual moves by about as much as its own
numerical jitter, ~3e-8. The Jacobian is then noise. Measured at that start
on three neighboring compositions, the flux row differentiated at the
default step reads −0.98, +0.28 and −4.92 where a converged central
difference gives −1.21, −1.22 and −1.30, and the small cross entry comes out
with the wrong sign. The first Newton step goes several decades in
`log10 T_surf`, and the profile is then asked for at a temperature the
thermodynamic polynomials — which stop at 6000 K — cannot supply, so
`make_profile_bg_gas` finds no valid point to bracket and reports that it
could not bracket the background pressure. That is where
`Could not bracket background pressure` comes from: not the bracketing scan.

`epsfcn = 1e-8` differences at 1e-4 relative, a 0.16 K step. Measured on this
build, that is about 1e4 above the jitter and the forward difference matches
a converged central difference to 0.1–0.5 percent at every composition
tested. The value is not chosen by which compositions solve: `epsfcn` =
1e-10, 1e-8, 1e-6 and 1e-4 all solve the same 61-point He/H grid
(8.00–11.00 in steps of 0.05), while the shipped `hybrd1` refuses 9.40, 9.50
and 10.45 and Photochem 0.8.4 refuses 8.10, 9.45, 9.85 and 10.15. It is
chosen by the accuracy of the Jacobian, which the grid does not see: at
`epsfcn = 1e-4` (a 16 K step) the flux row is 9 percent off the derivative
and the cross entry 38 percent off, and at 1e-12 and below the cross entry is
back in the noise. 1e-8 is the middle of the window where neither error
source is active.

Where the old build solves, the new one returns the same answer: over the
61-point grid the deep temperature agrees with the shipped `hybrd1` build to
1e-4 K and with Photochem 0.8.4 to 1.3e-4 K, both of which are the resolution
the sweep prints; read at full precision at four compositions the difference
against `hybrd1` is 2e-6 K (424.9534068835 against 424.9534089414 K at
He/H = 9.0). The three compositions `hybrd1` refuses come out at 422.0081,
421.3037 and 415.1543 K, on the smooth trend through their neighbours and
equal to 0.8.4's own answers at 9.40 and 9.50 to 3e-5 K. Across
He/H 9.30–9.60 in steps of 0.01 the second difference of the deep
temperature does not exceed 2e-4 K anywhere, including at the two points
that used to fail.

`AdiabatClimate_make_column`, the third `hybrd1` call in this file, is left
alone. Nothing measured here concerns it.

### 4.3 `photochem/extensions/gasgiants.py`

- `EvoAtmosphereGasGiant` and `GasGiantData` take `equilibrium_mass_tol`
  (default `1e-12`) and set it on the equilibrium solver; a non-finite or
  non-positive value is rejected at construction.
- `composition_at_metallicity` verifies elemental closure itself after each
  accepted solve, against `max(1e-8, 100*gas.mass_tol)`, and raises rather
  than returning a column that silently carries a large trace-element
  residual. **The check is stated on the total elemental vector, gas plus
  condensate, and is skipped at any level whose solution has a condensed
  phase** — where a condensate is present neither reported vector can be
  compared with the request (`molfracs_atoms_gas` is renormalized over the gas
  alone; `molfracs_atoms` keeps the previous level's condensate when the solve
  starts from it). Without that skip the check fires on every level above the
  first condensing one and no column with a cold trap can be built at all.
  Every level deep enough to set an elemental handoff is condensate-free,
  which is where the residual this check exists to catch appears.
- `molfracs_atoms = gas.molfracs_atoms_sun` is copied before being scaled by
  the metallicity, so repeated calls no longer accumulate onto the solver's
  own solar vector.

### 4.4 `tests/test_python.py`

Two tests for the above: the default `mass_tol` is `1e-12` after gas-giant
construction, and zero, negative, infinite and NaN `equilibrium_mass_tol` are
all rejected.

### 4.5 Nothing else

Measured 2026-08-29: a fresh clone at the pinned commit, patched by
`setup_photochem.sh`, against the `photochem/` tree the results were produced
with, `diff -r` excluding `.git/` reports **no content difference in any
file**. The tree is upstream text kept verbatim, so it carries upstream's own
wording — `v1.0_roadmap.md` uses two words this repository does not use in its
own documents. They are left alone: an edit there would cost the property this
arrangement exists for, and the file is not part of this repository.


## 5. Verify

Three checks, in order. Run all of them from a directory that is **not** the
Photochem source tree.

**5.1 The source is in the expected state.** `setup_photochem.sh` checks
eleven things in the files themselves — the settable equilibrium tolerance,
the closure check and its condensate skip, the copied solar vector, the Clima
`v0.7.5` pin, both dependency patches wired into CPM and present on disk, the
`epsfcn = 1.0e-8_dp` line, and the tolerance-validation test. It prints one
`ok` or `MISSING` line each and exits non-zero on any miss.

**5.2 The installed library is the patched one.** This is the check that
matters: 5.1 reads the source, and a build can be installed from a source
that was patched too late (section 3.2).

```bash
cd /tmp && python3 -c "
import photochem, inspect
from photochem import equilibrate, clima
import photochem.extensions.gasgiants as gg
print(photochem.__version__, equilibrate.__version__, clima.__version__)
print(inspect.signature(gg.GasGiantData.__init__))
print('molfracs_atoms_condensate' in inspect.getsource(gg.composition_at_metallicity))
"
```

Expected: `0.9.0 0.2.2 0.7.5`; `equilibrium_mass_tol=1e-12` in the signature;
`True` for the condensate guard. The guard is not optional — without it the
closure check fires at every level above the first condensing one, and no
column with a cold trap can be built at all.

The Clima repair is in a compiled Fortran extension and cannot be read from
Python source. Either disassemble the installed library for its two
constants — `epsfcn = 1e-08` and `factor = 100.0` — or, more directly, solve
the compositions the shipped `hybrd1` refuses: He/H = 9.40, 9.45 and 9.50 all
succeed on a repaired build and fail on an unrepaired one.

**5.3 Elemental closure at the deep states a handoff is set from.** This is
what the corrections were made for:

```bash
cd /tmp && python3 \
  <EXHALE>/LHS1140b/exhale/nh_refusal_diagnosis/closure_saved_states.py \
  <EXHALE>/LHS1140b/exhale/nh_refusal_diagnosis/scan
```

Every state should be accepted with a residual of order `1e-14`. Measured
2026-08-29: 39 states, largest relative residual `3.54e-13`, median
`1.90e-14`. On an unpatched Photochem the same script reports `1e-4`–`2e-4`,
which is what the adapter refuses.

**5.4 End to end, against a stored result.** The strongest check available,
and the one that says the physics is unchanged: re-run a converged elemental
closure and compare its He 10830 equivalent width with the stored value.
Measured 2026-08-29 on the LHS 1140 b He/H = 9.0 and 9.1 rungs, on a build
made by exactly the procedure of section 3: EW = 1.1068539227 and
1.1092155563 %A, agreeing with the stored reference in every printed digit
(relative 5e-12 and 2e-11). The stored reference is
`LHS1140b/exhale/clima_epsfcn_installed/results.txt`; the verification run
itself was deleted once its numbers were recorded here.

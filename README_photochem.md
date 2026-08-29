# The Photochem build EXHALE runs

Date: 2026-08-29

`photochem/` at the top of the EXHALE working tree is the Photochem source
EXHALE's lower-atmosphere handoff is actually run with, and the environment
built from it lives at `env/photochem`. Anything else on this machine — the
conda package `photochem 0.8.4` in `~/.conda/envs/photochem_cmp`, or a
`conda install -c conda-forge photochem` — is a *different* code that gives
different elemental closure at the handoff level, so results are not
interchangeable between them. Which build wrote a given handoff is recorded
per file, in the `# source_version photochem <version>` header line of every
`lower_atmosphere_profile.dat`.

`photochem/` is upstream Photochem with the four changes of section 2 applied
and nothing else. **It is not tracked**: `.gitignore` excludes the whole
directory, exactly as it excludes `VULCAN/`, and `src/utils/setup_photochem.sh`
puts it back from the upstream URL. Keeping it a clone plus one patch is the
point of the arrangement — `diff -r` against a clean checkout of the tag in
section 1 shows exactly section 2 and nothing else, and this file is the only
place any of it is written down.

The physics of the two solver corrections, why they were needed, and the
measurements that closed them are in
`docs/photochem_solver_modification_investigation.md` and
`docs/photochem_solver_modification_implementation.md`. This file is
provenance and build instructions only.


## 1. Where it came from

| | |
|---|---|
| upstream | <https://github.com/Nicholaswogan/photochem> |
| tag | `v0.9.0` |
| commit | `e1e872528b61e8dd1db891b84869738e219b4f97` (2026-08-07) |
| put in place | 2026-08-29 |

`v0.9.0` is what `origin/main` pointed at on that date. `origin/dev` was
inspected and not taken: it is not a release, and it still has the solver
behavior section 2 corrects.

`src/utils/setup_photochem.sh <dir>` builds the tree from nothing but the
upstream URL: it clones, checks out the commit above, applies
`src/utils/photochem_exhale.patch` — which is the whole of section 2 in one
file — and verifies eleven things in the files themselves. With no `<dir>` it
acts on `photochem/`. It refuses a tree at any other commit and is safe to
re-run; on an already-patched tree it verifies and exits without touching a
file.

Measured 2026-08-29: a fresh clone at that commit, patched by that script,
against `photochem/` in place, `diff -r` excluding `.git/` reports **no
content difference in any file**. The working copy in place was made by
copying rather than cloning, so it lacks the upstream repository metadata that
has no function in the build — `.gitattributes` (GitHub language-statistics
directives), `.github/workflows/test.yaml`, and an empty `.gitmodules`. A tree
produced by the script keeps those, and `.git` with them. The upstream
`.gitignore` is present either way, and is the right one for a build done in
place: its patterns are `build`, `_skbuild`, `*.so`, `data`,
`photochem.egg-info`.

The upstream clone the working copy was made from, `photochem/` one level
above this repository in the `ExoAtmosphere` workspace, still exists with its
history and is untouched.

The tree is upstream text kept verbatim, so it carries upstream's own wording
— `v1.0_roadmap.md` uses two words this repository does not use in its own
documents. They are left alone: an edit there would cost the property this
arrangement exists for, and the file is not part of this repository.


## 2. What differs from upstream

Two of the four changes are corrections to Photochem's dependencies, carried
as patch files and applied by CPM at configure time; two are in Photochem's
own Python layer. Paths in this section are relative to `photochem/`.

### 2.1 `src/dependencies/patches/equilibrate-element-relative-mass-closure.patch`

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

### 2.2 `src/dependencies/patches/clima-bounded-scaled-solvers.patch`

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

### 2.3 `photochem/extensions/gasgiants.py`

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

### 2.4 `tests/test_python.py`

Two tests for the above: the default `mass_tol` is `1e-12` after gas-giant
construction, and zero, negative, infinite and NaN `equilibrium_mass_tol` are
all rejected.


## 3. The environment

`env/photochem` is a conda environment holding the wheel built from that
source. It is what `src/utils/element_flux_closure.py` runs the chemistry in
when a `closure.json` does not name an interpreter of its own, and it is
resolved relative to the repository, not by absolute path.

```
env/photochem/bin/python      photochem 0.9.0, Equilibrate 0.2.2, Clima 0.7.5
```

About 1.6 GB, so it is in `.gitignore` (`env/`) along with `build/` and the
other build products — as is the `photochem/` source it is built from. Both
are made locally; do not expect to clone the repository and find either.

### Rebuilding it

`src/utils/setup_photochem.sh <dir>` reproduces the source, as section 1
describes, and verifies the result. It does not build.

The full procedure — the script, the toolchain (CMake **3.31.8**; CMake 4 does
not work, CVODE 5.7 is not compatible with it), the wheel build, the
environment, and the two verification steps — is one section of
`README_HOWTO.md`, "The Photochem environment". It is not repeated here.
Three things from it are worth having in front of you before you start:

- `/usr/include/numpy` on this machine is a **broken symlink into a Python 2.7
  tree** and scikit-build finds it first, so the build needs
  `-DNumPy_INCLUDE_DIR=$(python -c 'import numpy; print(numpy.get_include())')`.
- a user-site NumPy 1.26.4 **shadows the conda NumPy of every environment on
  this machine**, so the build needs `PYTHONNOUSERSITE=1`. It shadows at run
  time too; Photochem imports and runs either way, and the stored LHS 1140 b
  results were made under exactly that condition, so runs are deliberately
  left alone.
- CMake tries `f95` before `gfortran`, and `/usr/bin/f95` here is GNU 13.1,
  below the 14.0 Photochem requires. Putting the toolchain environment first
  on `PATH` does not stop it; `-DCMAKE_Fortran_COMPILER` and
  `-DCMAKE_C_COMPILER` have to name the compilers.

The environment at `env/photochem` was first made by cloning one that
already had a wheel, which is much cheaper than building:

```sh
conda create -p <EXHALE>/env/photochem --clone photochem_090_fix
```

It now carries a wheel built on 2026-08-29 from that source in
`~/.conda/envs/photochem_build_tc`, which is the change of section 2.2 to the
surface-temperature solve. `photochem_090_fix` is untouched and still holds
the build before it, so re-cloning is the way back.

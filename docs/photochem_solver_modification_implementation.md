# Photochem solver modification implementation

Date: 2026-08-29

## The judgement first

The deep elemental residual diagnosed in `deep_level_elemental_check.md` is
removed at its source, and **it does not move the observable**. Re-measured on
the corrected build with 11 arms, the LHS 1140 b flux-closure crossing of the
He I 10830 equivalent width is at He/H = **9.0484** (chord) / **9.0461**
(quadratic) -- the same to every digit the photochem 0.8.4 ladder printed. On
seven arms that are like-for-like the equivalent width moves by at most
`1.2e-6` in relative terms and `log10 Mdot` does not move in its fourth
decimal. The residual was real, the correction is right, and the science answer
is unchanged. Record: `LHS1140b/exhale/crossings_pc090/results.txt`.

Two things must be read with that conclusion:

- The first version of the strict closure test **rejected every composition**,
  including ones that ran on 0.8.4. It was a blocking defect, it was found only
  when a complete column was run, and it is fixed. Section 3.
- The corrected clima refuses two compositions the old build solved,
  He/H = 9.4 and 9.5. That is **not a regression of this change**: the same
  solver defect is in 0.8.4, which refuses two others on the same grid.
  Diagnosed. Section 4.4, and `Update_EXHALE.md` section 91.

## 1. What was built

The current stable Photochem release, `v0.9.0` at commit
`e1e872528b61e8dd1db891b84869738e219b4f97`, was fetched and modified. The
official `origin/main` points to this same commit. The newer `origin/dev`
branch was inspected but was not selected as a stable release; it also retained
the solver behavior diagnosed in
`photochem_solver_modification_investigation.md`.

The modified source builds as a Python 3.11 wheel with GNU Fortran 14.4 and
CMake 3.31.8. The build selects clima `v0.7.5`, the latest stable clima release
available on the implementation date, and Equilibrate `v0.2.2`. Reproducible
patches for both dependencies are applied by Photochem's CMake configuration.

The wheel is installed into a dedicated 0.9.0 environment. The 0.8.4
environment of the earlier comparisons is untouched and remains the
reproduction path for everything measured before 2026-08-29 11:30.

## 2. Implemented changes

### 2.1 Equilibrate elemental convergence

`src/dependencies/patches/equilibrate-element-relative-mass-closure.patch`
changes both CEA update paths. Each positive requested elemental abundance is
now tested against its own relative tolerance:

```fortran
mval_mass_good = abs(b_0(i_atom))*self%mass_tol
```

This replaces the common absolute threshold derived from `MAXVAL(b_0)` and
removes the exclusion of elements at or below `1e-6`. Therefore,
`converged=True` requires the documented relative mass accuracy for every
element present in the requested composition.

### 2.2 Photochem gas-giant equilibrium enforcement

`photochem/extensions/gasgiants.py` now:

1. exposes `equilibrium_mass_tol` on `EvoAtmosphereGasGiant`, with a default of
   `1e-12` and finite-positive validation;
2. assigns that tolerance to `ChemEquiAnalysis.mass_tol` before equilibrium
   initialization;
3. copies `molfracs_atoms_sun` before metallicity and C/O transformations,
   preventing mutation of the stored reference composition;
4. tests the returned elemental vector against the requested one **at levels
   whose solution has no condensed phase**, and accepts only a largest relative
   residual no greater than `closure_rtol = max(1e-8, 100*mass_tol)`;
5. raises a level-specific `RuntimeError` after all five temperature
   perturbations fail, instead of silently retaining the last result; and
6. restores `use_prev_guess` in a `finally` block.

Item 4 is stated exactly as the code states it. `gas.molfracs_atoms_condensate`
is read first; if any entry is positive the level is accepted on the solver's
own `converged` flag and no closure residual is formed. Only a condensate-free
level is compared, and it is compared on the total elemental vector
`gas.molfracs_atoms`. Section 3 is why.

### 2.3 Clima pressure and temperature solves

`src/dependencies/patches/clima-bounded-scaled-solvers.patch` makes three
changes to clima `v0.7.5`.

First, the scalar background-pressure MINPACK solve is replaced with a scan
and bracketed bisection in the physical interval
`P_surface*1e-12 <= P_background <= P_surface`. Invalid thermodynamic trial
states are skipped during bracketing instead of terminating the complete
climate solve. The normalized surface-pressure tolerance is `1e-10`.

Second, a coupled temperature solve uses

```text
log10(T_surface - T_trop), log10(T_trop)
```

when `solve_for_T_trop` is enabled. Every trial therefore has positive
temperatures and satisfies `T_surface > T_trop`.

Third, the energy residual is normalized by the stellar-flux scale and the
skin-temperature residual is normalized by `max(T_trop, 1 K)`. After MINPACK
returns, both normalized physical residuals must be at most `1e-7`; otherwise
the solution is rejected.

`src/dependencies/CMakeLists.txt` now selects clima `v0.7.5` and applies both
dependency patches through CPM.

## 3. The blocking defect in the first closure test, and its fix

**It rejected every composition it was given, including compositions that ran
on 0.8.4.** The test as first written compared the requested elemental vector
with `gas.molfracs_atoms`, the total over gas and condensate, at every level.
That comparison is not available at a level that has a condensed phase, and the
LHS 1140 b column has one.

The mechanism, measured on the LHS 1140 b column at He/H = 9:

- `composition_at_metallicity` is called with `rainout_condensed_atoms=True`,
  so what it hands to the next level is `gas.molfracs_atoms_gas`, the
  **gas-only** elemental vector, renormalized over the gas alone.
- The test read `gas.molfracs_atoms`, the **total**. Above the first condensing
  level those two are different objects: solved from scratch the total closes
  to `5.7e-14` against the request, but with `use_prev_guess=True` the solve is
  started from the previous level and its reported total still carries that
  level's condensate, so the same quantity reads **`1.1e+2`**.
- Read the other way round the mismatch reverses: at the first condensing
  level the total closes to `2e-14` while the gas-only vector is off by
  `9.9e-1`; at every level above it the total is off by `1.1e+2` while the
  gas-only vector closes to `6e-14`.

**This is Equilibrate's existing reporting behavior, not something the Fortran
patch introduced.** The same numbers are obtained with the unpatched 0.8.4
build and with the patched build.

The fix, already applied in `photochem/photochem/extensions/gasgiants.py`,
forms the residual **only where the reported solution has no condensed phase**
(`molfracs_atoms_condensate` all zero); a level with a condensate is accepted
on the solver's `converged` flag alone. Every level deep enough to set an
elemental handoff is condensate-free, and that is also where the residual this
test exists to catch appears, so the test keeps its detection power exactly
where it is read out.

## 4. Measured validation

### 4.1 Build and Python interface

The repository was built through its Python wheel path without local source
overrides. This verifies that CPM downloaded the selected release tags and
applied the source-tree patches. The resulting artifact was:

```text
photochem-0.9.0-cp311-cp311-linux_x86_64.whl
```

Two focused Python tests passed:

- gas-giant construction uses `mass_tol = 1e-12` by default;
- zero, negative, infinite, and NaN equilibrium tolerances are rejected.

Python compilation and `git diff --check` also passed.

### 4.2 Deep elemental closure

The installed wheel was applied to all 38 saved deep states under
`LHS1140b/exhale/nh_refusal_diagnosis/scan/`. Each test used the saved
pressure, temperature, elemental composition, and thermodynamic file.

| Quantity | Measured value |
|---|---:|
| converged and independently accepted cases | 38/38 |
| largest relative elemental residual | `3.539391003e-13` |
| He/H at the largest residual | `11.0` |
| median largest residual | `1.915134717e-14` |
| the same measure on the 0.8.4 build | `2.018197773e-4` |

Every value in that table was reproduced independently.

**The scope of that test, stated plainly: it ran one level for each state and
never carried a column.** Each of the 38 cases is a single call at the saved
deep pressure and temperature. `composition_at_metallicity` was not exercised
over a full climate grid, `use_prev_guess` chaining was not exercised, and no
level with a condensed phase was reached. That is exactly why the defect of
section 3 survived this test and appeared only when a complete column was run.
A reading of "38/38 converged and independently accepted" as evidence that the
closure test is safe on a column is not supported by this measurement; the
column evidence is section 5.

### 4.3 Radiative-convective solutions

The wheel was tested with the saved EXHALE stellar spectra, pressure, planet
properties, composition partition, and original temperature guesses. The table
also includes the successful neighboring compositions used to test continuity.

| He/H | deep temperature [K] | tropopause temperature [K] | flux residual [mW/m2] | normalized flux residual |
|---:|---:|---:|---:|---:|
| 0.0969 | 720.483384 | 184.226352 | `+4.185577e-3` | `+7.682239e-9` |
| 8.0 | 433.261442 | 185.421721 | `-2.773013e-7` | `-5.089609e-13` |
| 8.1 | 432.371517 | 185.424045 | `-2.999713e-4` | `-5.505696e-10` |
| 8.2 | 431.500985 | 185.426329 | `+8.872694e-4` | `+1.628501e-9` |
| 15.0 | 393.126945 | 185.519683 | `+3.546445e-4` | `+6.509173e-10` |

The three deep temperatures at 0.0969, 8.1 and 15 were reproduced
independently as 720.5, 432.4 and 393.1 K. The He/H = 8.1 temperatures lie
between the 8.0 and 8.2 results. The final flux residuals are all below the
normalized acceptance limit by more than an order of magnitude.

All three clima temperature entry points affected by the shared solver were
also exercised with the upstream Earth test inputs:

| API | deep temperature [K] | tropopause temperature [K] |
|---|---:|---:|
| `surface_temperature_column` | 289.005470 | 218.549434 |
| `surface_temperature_bg_gas` | 291.703216 | 218.805395 |
| `surface_temperature` with dayside redistribution | 431.462187 | 224.974987 |

### 4.4 Two compositions this build refuses, and why it is not a regression

**He/H = 9.4 and 9.5 fail on the corrected build and solved on 0.8.4.** Both
raise

```text
_clima.ClimaException: hybrd1 root solve failed: Could not bracket background
pressure in make_profile_bg_gas.
```

The first version of this section read that message at face value -- an
isolated failure of the bracketing scan of section 2.3, cause not diagnosed --
and counted a net of three failures removed and one introduced. **The
measurement contradicts both.** It is recorded in
`LHS1140b/exhale/clima_bracket_diagnosis/` and written up as section 91 of
`docs/Update_EXHALE.md`; the three points that matter here are:

- **The bracketing scan is not the cause.** The root is inside the interval it
  scans and the residual crosses zero once, monotonically, in the last scan
  step. The scan reports failure because the outer
  `surface_temperature_bg_gas` solve hands it `T_surf = 6.5e8` K, at which all
  49 scan points are rejected with `Failed to compute heat capacity` -- the
  thermodynamic polynomials stop at 6000 K. The cause is the outer solve's
  forward-difference Jacobian, formed by `hybrd1` at `sqrt(epsmch)*|x|`
  (`2.4e-5` K at `T_surf = 400` K), a step over which the flux residual moves
  by `3.8e-8` against its own `~3e-8` jitter.
- **0.8.4 has the same defect at the same rate**, at different compositions.
  Measured on two grids:

| grid | corrected build | 0.8.4 |
|---|---|---|
| He/H 9.30--9.60, step 0.01 (31 points) | 9.40, 9.50 | 9.37, 9.45 |
| He/H 8.00--11.00, step 0.05 (61 points) | 9.40, 9.50, 10.45 | 8.10, 9.45, 9.85, 10.15 |

  The `log10(T_surf - T_trop)` parameterization of section 2.3 makes 0.8.4's
  `T_surf is less than T_trop` refusal unreachable, so the same bad step
  surfaces as a bracket message instead. It changed which compositions are
  unlucky, not how many.
- **The net accounting is withdrawn.** "Three removed, one introduced" compared
  this build against 0.8.4 columns that had not been run at the matching
  reservoir values; 9.45, recorded there as `not attempted` on 0.8.4, is a
  failure on 0.8.4.

Neither composition is a bracket arm of the crossing, and the solution is
smooth across them (deep temperature falls monotonically from 422.72470 K at
He/H = 9.30 to 420.61128 K at 9.60). Setting `epsfcn = 1e-4`, which requires
`hybrd` in place of `hybrd1`, takes the 61-point grid from 54/61 to 61/61 in
replication and does the same on 0.8.4's parameterization. That change is being
applied to the clima patch and verified separately; nothing here is measured on
a rebuilt build. Item (M) of `TO_BE_DONE.md`. A composition that refuses can be
recovered today by passing `--climate-t-deep-guess` anything other than its
default 400.0.

## 5. Downstream: what the correction does to the LHS 1140 b result

This is the part the earlier version of this document listed as not
regenerated. It has now been done, and it is the scientific conclusion of the
whole change.

The flux-closure crossing was re-determined end to end on the corrected build
with 11 arms, using the same EXHALE binary, the same scripts, the same single
parent seed and the same tolerances as the stored 0.8.4 ladder -- only the
interpreter differs. Record:
`LHS1140b/exhale/crossings_pc090/results.txt`.

| quantity | 0.8.4, 16 arms | corrected build, 11 arms |
|---|---:|---:|
| crossing, chord | 9.0484 | 9.0484 |
| crossing, quadratic | 9.0461 | 9.0461 |
| 1 sigma band | 8.4923 -- 10.7716 | 8.4649 -- 10.7357 |

**The crossing does not move to any printed digit.** On the seven arms seeded
identically in both builds the equivalent width moves by at most `1.2e-6` in
relative terms and `log10 Mdot` is unchanged in its fourth decimal. The handoff
profile itself does move -- at He/H = 9 the deepest-level N/H departure goes
from `2.1e-5` to `1.2e-14`, and trace species such as `q_H2O`, `q_CO` and
`q_HCN` move by up to `1e-2` relative -- but none of that reaches He I 10830,
which is set by H, He and the temperature.

**The 1 sigma band moves for a different reason than the solver.** Its low end
was previously interpolated across the band He/H = 8.5 -- 8.7, which the
handoff check refused and which no arm covered; the corrected build runs that
band (departures `1.4e-14` to `3.2e-14`), so real arms at 8.5 and 8.6 now carry
that end and it moves from 8.4923 to 8.4649. The high end moves from 10.7716 to
10.7357 because the 10.6 arm was not seeded the same way in the two ladders --
a seed difference of `1.1e-3` in equivalent width, the same size as the seed
spread the stored ladder already measured, not a solver difference.

## 6. The handoff tolerance

`src/utils/lower_profile_schema.py`: the `--abundance-tol` default is
**`1.0e-3` -> `1.0e-10`**, the value the check was designed at. The grounds are
measurement on the corrected build, not judgement:

- Over the **68 handoff writes** made by the corrected build under
  `LHS1140b/exhale/` (43 in `nh_refusal_diagnosis/scan_pc090/`, 25 in
  `crossings_pc090/`), the largest deepest-level departure is **`3.28e-13`**
  and the median **`1.73e-14`**. `1e-10` stands a factor 305 above the worst of
  those, and a factor 280 above the worst of the 38 deep states of section 4.2.
- Detection power is kept, not given up. Each element sits in one dominant
  carrier at this level, so a genuine miscount -- a carrier dropped from the
  sum -- is an O(1) error. The smallest conceivable one, dropping N2, is
  `1.0e-5` of `X_N` in the coldest band and `4.0e-3` where the deep column is
  hot enough to make N2; both are far above `1e-10`.
- **A 0.8.4 handoff is refused at this default, and that is the intended
  signal.** Over the 43 stored 0.8.4 handoffs of `crossings_gm25/` the
  departures run `3.48e-8` (best) to `8.50e-5` (worst), median `6.21e-6`: even
  the best 0.8.4 write misses `1e-10` by a factor 348. Reproducing a stored
  0.8.4 result therefore requires passing `--abundance-tol` explicitly, which
  is how a deliberate reproduction states that it accepts that residual. The
  scans under `nh_refusal_diagnosis/` already pass `--abundance-tol 1.0` on
  their own command lines and are unaffected.

The earlier `1.0e-3`, and the reasoning that set it, were correct **for the
0.8.4 solver**: with elemental residuals scaled by the largest elemental
abundance, the check could not be tightened without refusing converged columns.
That constraint is gone with the corrected build.

## 7. Validation scope and remaining work

Tested: dependency patch application, gas-giant equilibrium initialization,
all 38 saved deep equilibrium states, all three recorded climate failures,
neighboring climate states, all shared temperature-solve entry points, a
45-point reservoir scan carrying complete columns, and 11 complete flux-closure
arms through to synthesized He I 10830 equivalent widths.

Not tested, and not claimed:

- **The repair for the clima temperature solve on a rebuilt build.** Why
  He/H = 9.4 and 9.5 refuse is diagnosed (section 4.4), and `epsfcn = 1e-4`
  solves the whole grid in replication, but it has not been applied to the
  patch or re-measured in an installed environment.
- Condensate-rich and sulfur-network equilibrium sweeps. Note that after the
  section 3 fix a condensing level is accepted on the solver's own `converged`
  flag, so the closure test does not speak for such a level at all.
- The stored 0.8.4 ladder was not re-run; the seven like-for-like arms are the
  comparison. He/H = 9.5 cannot be repeated on this build without the
  `--climate-t-deep-guess` workaround of section 4.4.
- The seed systematic was not re-measured on the corrected build.
- No EXHALE Fortran source or binary changed, no golden was refreshed, and
  `make check` was not run: this path is not in the regression matrix.

The temporary validation environment required CMake 3.31.8 because CMake 4
removes compatibility behavior used by CVODE 5.7, an unchanged Photochem
dependency. This is a build-tool compatibility issue separate from the solver
patches.

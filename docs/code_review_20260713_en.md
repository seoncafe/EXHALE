# EXHALE Code Review Report (English)

- Review date: 2026-07-13
- Reviewed commit: `cfc33b1fe36143e6192f2d854805e506ddcda966` (`update`)
- Scope: main solver, Wind-AE initial-condition code, steady solver, molecular and metal chemistry, Python post-processing, build and run scripts
- Methods: source review, warning-enabled debug builds, Python/shell syntax checks, a valid-input smoke run, and reproduction of an invalid configuration

## 1. Executive summary

The current tree completes a warning-enabled debug build, and `tutorial_nometals` ran through 474 steps without a bounds-check or floating-point exception. Earlier fixes for interpolation NaNs, the local He/H value in post-processing, and `q_abs`/logarithm guards remain in place. The former OpenMP `critical` section in the radiation loop has also been removed, so one previously reported serialization bottleneck is no longer present.

The following items should receive the highest priority.

| Priority | Finding | Confidence | Main impact |
|---|---|---|---|
| P0 | Molecular-chemistry compatibility checks run before the options are parsed | Reproduced | A molecular+metal input passes validation and then crashes out of bounds |
| P0 | The same array is passed to `intent(in)` and `intent(out)` arguments | Fortran conformance defect | Results may depend on compiler and optimization level |
| P1 | Loaded-SED lower bound ignores low-ionization-potential metals | High-confidence physics/logic defect | Photoionization of Mg, Si, Ca, Na, K, and Fe can be underestimated |
| P1 | The monochromatic H-only branch treats helium inconsistently | High-confidence state inconsistency | Initial density, particle count, and chemistry assume different compositions |
| P1 | Transit auto-window column calculation omits the metre-to-centimetre conversion | Confirmed arithmetic defect | The column is 100 times too small, potentially clipping broad wings |
| P1 | The SED reader lacks safe EOF, row-count, and ordering checks | High-confidence input defect | It can loop indefinitely or access `e_v(2)` out of bounds |
| P1 | HLLC/PLM hot loops create noncontiguous row-section temporaries | Observed at runtime | Per-cell and per-interface copies add avoidable overhead |
| P2 | User-facing string selections lack `case default` validation | High-confidence input defect | A typo can propagate as undefined state instead of failing immediately |
| P2 | The build stamp does not track `FFLAGS` changes | Reproduced | Debug/release flag changes can silently reuse stale objects |
| P2 | Constants and external calls generate substantial precision/interface warnings | Build evidence | Reduced numerical accuracy and missed call-signature errors |

The recommended order is: fix both P0 items; correct the P1 physics and unit defects; harden SED handling; profile and remove hot-loop temporaries; then address build hygiene and structural duplication.

## 2. Review and verification results

### 2.1 Warning-enabled builds

The complete Fortran code was forcibly rebuilt in a separate `/tmp` object directory with options equivalent to:

```text
-O0 -g -fopenmp -Wall -Wextra -Wimplicit-interface
-Wconversion-extra -Wsurprising -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow
```

The build succeeded. It produced 4,440 warnings, dominated by the following categories.

| Warning category | Count | Interpretation |
|---|---:|---|
| Tab characters | 2,779 | Mostly formatting debt, but it obscures meaningful warnings |
| Precision conversions | 1,417 | Many default-real literals are assigned to `real*8` values |
| Implicit interfaces | 51 | Type and rank checking is incomplete for MINPACK/LAPACK calls |
| Real equality comparisons | 51 | Some are intentional sentinel checks; others merit cleanup |
| Unused variables | 40 | Candidates for removing old paths and duplicate implementations |
| Do-subscript warnings | 26 | The inspected main cases were endpoint-branch false positives |

The standalone Wind-AE build also succeeded. Every tracked Python file passed `py_compile`, and the tracked shell scripts passed `bash -n`.

### 2.2 Runtime verification

The debug executable ran `examples/tutorial_nometals/input.inp` for 20 seconds and reached 474 steps without a bounds or FPE failure. The runtime did, however, repeatedly report array temporaries for noncontiguous sections at:

- `src/modules/states/PLM_rec.f90:26`
- `src/modules/time_step/RK_rhs.f90:41,49,53`

An input with both molecular chemistry and metals enabled was not rejected by the parser. It entered time integration and then failed with:

```text
Fortran runtime error: Index '9' of dimension 1 of array 'sys_x'
above upper bound of 8
at src/modules/radiation/ionization_equilibrium.f90:409
```

Finding 3.1 is therefore a reproduced defect, not merely a static concern.

## 3. Correctness and robustness findings

### 3.1 [P0, reproduced] Molecular compatibility validation runs at the wrong time

Evidence:

- `src/modules/files_IO/input_read.f90:262-276` checks molecular chemistry against helium, helium diffusion, and metals.
- `Molecular chemistry` is not parsed until `:394-397`, and `He diffusion` until `:428-432`.
- Both flags therefore still hold their default values when the checks execute.
- The reproduced molecular+metal case subsequently accessed `sys_x(9)` at `ionization_equilibrium.f90:409`, while the upper bound was 8.

Impact:

- The documented invalid combinations `molecular + no He`, `molecular + He diffusion`, and `molecular + metals` can pass the parser.
- Failure occurs deep inside time integration rather than at configuration loading, obscuring the cause.

Recommended fix:

1. Call one `validate_configuration()` after every optional field and `base.inp` override have been read, but before allocation and initialization.
2. Reject invalid combinations with a descriptive, nonzero `error stop`.
3. Keep cross-option validation in this single routine instead of distributing it among parsing locations.

Required regressions:

- Molecular+metals, molecular+He diffusion, and molecular+HeH=0 must all fail during input validation.
- The supported molecular example must continue to initialize normally.

### 3.2 [P0] Aliasing between `intent(in)` and `intent(out)` arguments

Representative cases include:

- `Apply_BC(u,u)` throughout `EXHALE_main.f90`, `steady_newton.f90`, and `init.f90`.
- `Apply_BC` itself calls `Apply_BC_W(W,W)`.
- `ioniz_eq(T,rho,f_sp,rho,f_sp,...)` aliases its input `rho/f_sp` with `rho_out/f_sp_out`.
- `post_process.f90` passes the same scratch variable `dum_v` to several distinct `intent(out)` arguments in a single call.

Fortran restricts access through one dummy argument while the same actual object is being defined through another dummy. Even if the present implementation happens to copy inputs into locals before producing outputs, the interface contract is nonconforming and an optimizing compiler may assume the arguments do not alias.

Recommended fix:

- Provide a single explicit `intent(inout)` API for truly in-place boundary operations.
- If both copy and in-place behavior are needed, use separate-buffer and in-place wrappers.
- Call `ioniz_eq` with `rho_new` and `f_sp_new`, then assign them explicitly.
- Use separate locals for discarded outputs, or provide optional outputs/a dedicated rates API.
- Run golden regressions under gfortran `-O0/-O3` and, where available, Intel ifx.

### 3.3 [P1] Loaded SEDs discard photons needed by low-IP metals

The power-law path in `src/modules/radiation/set_energy_vectors.f90:31-53` lowers the energy-grid floor to active metal ionization thresholds when `thereis_lowIP_metal` is true. The loaded-SED path in `src/modules/radiation/sed_read.f90:32-33` accounts only for the He I triplet and does not apply the low-IP metal rule.

Impact:

- Photons below 13.6 eV may be discarded even when they exist in the SED file.
- Photoionization rates, electron density, cooling, and transit ion fractions for Mg, Si, Ca, Na, K, and Fe can all be biased.

Recommended fix:

- Derive the minimum required threshold from active `species_table` metadata and use the same policy for power-law and loaded spectra.
- Warn or fail explicitly when the supplied SED does not cover the required range.
- Add rate regressions with bins on both sides of each active low-IP threshold.

### 3.4 [P1] Inconsistent helium state in the monochromatic H-only branch

`input_read.f90:174-176` sets `thereis_He=.False.` when the monochromatic photon energy is below 24.6 eV. However, `HeH`, `mass_per_H=1+4*HeH`, and the boundary density still include helium mass, and `set_IC.f90:151-156` still writes helium fractions. Conversely, `calc_rho` and `calc_ntot` in `utilities.f90` omit helium entirely when `thereis_He=False`.

This mixes two models: neutral helium that is present but not photoactive, and an H-only fluid with no helium.

Recommended fix:

- Physically, retain neutral-helium mass and particle count below 24.6 eV while disabling only its photoionization channel.
- If a true H-only calculation is intended, consistently set `HeH=0` and update `mass_per_H`, initial abundances, and the EOS.
- Separate `composition_has_He` from `He_photoactive` to make the distinction explicit.

### 3.5 [P1] Unsafe EOF and structural handling in the SED reader

`sed_read.f90:46` reads with `iostat`, uses the value at `:49`, may `cycle` at `:53`, and checks `iostat` only at `:59`. If every row lies outside the selected upper bound, EOF can retain an old/undefined value and repeatedly take the `cycle` path. The code also references `e_v(2)` at `:90` even if selection produced zero or one row.

Positive wavelengths, strict ordering, duplicate rows, and malformed rows are not validated either.

Recommended fix:

- Handle `iostat` immediately after each read and distinguish EOF from malformed data.
- Require at least two selected rows.
- Require positive, strictly increasing wavelengths.
- Include the filename, source row, and required energy range in failures.

### 3.6 [P1, confirmed] Transit auto-window column is off by a factor of 100

At `EXHALE_transit.py:406-409`, `_dr_cm = np.gradient(r) * Rp` is used. `Rp` is in metres, while density is in cm\(^{-3}\); the variable name and comment assume centimetres, but the `1e2` conversion is missing. The actual LOS integration later in the same file, at `:472-475`, correctly uses `r * Rp * 1e2`.

Impact:

- The column used for automatic wavelength-window sizing is 100 times too small.
- Because a damping-wing estimate scales roughly as \(\sqrt{N}\), the selected window can be about 10 times too narrow before its cap is applied.
- The primary risk is clipping He/Lyα wings; the actual LOS optical-depth calculation uses the correct conversion.

Use `_dr_cm = np.abs(np.gradient(r)) * Rp * 1e2` and add an analytic uniform-atmosphere unit test.

### 3.7 [P2] User-facing enum strings lack default error cases

The following `select case` statements do not provide a `case default`:

- reconstruction in `src/modules/states/Reconstruction.f90`
- numerical flux in `src/modules/states/Num_Fluxes.f90`
- grid type in `src/modules/grid/define_grid.f90`
- spectrum type in `src/modules/files_IO/input_read.f90:151-178`

A misspelled option can proceed without initializing flags or output arrays. Validate every user-facing selection immediately and report the received value plus the allowed values through `case default; error stop`.

### 3.8 [P2] Composition state should be finalized after `base.inp` overrides

The lower-atmosphere path can override `HeH` through `base.inp`. Final composition flags, `mass_per_H`, and molecular compatibility checks should be recomputed once from the post-override values. Partially updating this state during parsing makes it easy for `HeH`, `thereis_He`, and species counts to disagree.

## 4. Performance opportunities

### 4.1 [First target] Remove noncontiguous row slices from HLLC/PLM hot loops

State arrays use `(cell, component)` ordering, but procedures receive `W(j,:)` and `u(k,:)`. In Fortran column-major layout these sections are noncontiguous. The debug runtime confirmed repeated temporaries at `PLM_rec.f90:26` and `RK_rhs.f90:41,49,53`.

Recommended sequence:

1. Profile a production-sized case to measure the time and allocation share of flux/reconstruction.
2. As a local change, use scalar arguments or explicit length-three locals to remove hidden allocation/copy behavior.
3. Longer term, transpose state to `(component, cell)` or a structure-of-arrays layout and perform whole-range conversions once.
4. Compare mass/energy conservation and wall time before and after the change.

The layout change has a broad blast radius and should be a separate change after P0/P1 fixes.

### 4.2 Nested loops and repeated line calculations in `EXHALE_transit.py`

The script executes Python loops over impact parameter, wavelength, and line of sight. Some doublets are evaluated after a similar single-component pass, and coordinate arrays are grown with repeated `np.append` calls.

Recommended improvements:

- Compute chord geometry, temperature, bulk velocity, and abundance once per impact parameter.
- Vectorize the wavelength dimension in bounded-memory chunks.
- Replace `np.append` loops with preallocation or vector expressions.
- Store line definitions in a table and use one singlet/doublet solver.
- Record both spectral relative error and peak memory during optimization.

### 4.3 Repeated chemistry and cooling rates

The H/He, triplet, metal, and molecular systems recompute similar rates and electron densities for the same cell state. This should be confirmed with a profiler before refactoring. If it is material, use a named per-cell rate cache/derived type shared by residual and Jacobian evaluation. Its lifetime and temperature/density dependencies must be explicit to prevent stale-rate defects.

### 4.4 OpenMP status

The previous `critical` serialization in the radiation loop is no longer present. Further directives should follow measurements of 1/2/4/8-thread scaling, scheduling imbalance, and chemistry-solver cost. CI should compare one-thread and multithread results within a documented tolerance.

## 5. Duplication and simplification opportunities

### 5.1 Share the common HLLC and ROE wave-speed calculation

`speed_estimate_HLLC.f90` and `speed_estimate_ROE.f90` substantially duplicate primitive-state unpacking, PVRS estimates, TRRS/TSRS branches, and star pressure/velocity calculations.

A common routine can return:

- left/right primitive state,
- `p_star` and `u_star`, and
- shock/rarefaction correction factors.

Each solver would retain only its final wave-speed combination. Establish golden tests over representative Riemann states before refactoring to prevent changes at numerical branch boundaries.

### 5.2 Reduce duplication among H/He/TR/metal/molecular systems

The `System_HeH*`, `System_HeH_TR*`, metal, and implicit-advection variants repeat H/He reaction logic and pass position-dependent `params(N)` arrays. A single displaced index is not type-checkable.

Recommended structure:

- Define named fields in `type(ion_cell_state)` and `type(ion_rates)`.
- Split species contributions into small pure procedures.
- Assemble residuals/Jacobians by iterating active species metadata.
- Migrate incrementally—H/He common code, then metals, then molecules—with golden tests at each step.

The current omission of collisional ionization in the triplet path is documented as intentional in source comments. It should not be labelled a defect without a quantitative study; instead, document the temperature range over which its effect is negligible.

### 5.3 Centralize composition calculations

`calc_rho`, `calc_ntot`, initialization, boundaries, and chemistry each branch separately on helium, metals, and molecules. The inconsistency in 3.4 is a consequence of this distributed state logic. Extend `species_table` to H/He/molecules and centralize mass, particle-count, and electron-count calculations around the metadata.

### 5.4 Simplify the `post_process` output API

Callers alias `dum_v` across unwanted `intent(out)` values only because the APIs always return every component. Provide purpose-specific entry points—heating only, cooling components, diagnostic rates—or optional outputs. This removes the aliasing and makes calls self-documenting.

### 5.5 Modularize `EXHALE_transit.py`

The roughly 1,400-line top-level script mixes input parsing, physics, convolution, plotting, and interactive input.

Suggested modules are:

- `main()` plus argparse/config,
- EXHALE profile/header reader,
- line metadata table,
- pure Voigt/LOS/transit solver,
- instrument convolution, and
- plotting.

Metal columns should also be selected from the output `# columns` header rather than hard-coded indices such as 17, 23, and 25. A changed species order can currently make the script silently read the wrong column.

### 5.6 Unify input parsing

The Fortran parser mixes positional mandatory lines at the beginning with key-search optional lines later, while Python utilities parse the same file independently. A single key-value schema with validation and explicit legacy aliases would remove duplication and ordering sensitivity.

## 6. Build, portability, and maintenance

### 6.1 [Reproduced] `FFLAGS` changes do not trigger recompilation

The Makefile compiler stamp records the compiler but does not track all compilation flags as dependencies. After building objects with debug flags, running `make -n` in the same `OBJDIR` with `FFLAGS='-O3 -fopenmp'` returned `Nothing to be done for 'all'`.

Recommended fix:

- Generate a configuration stamp whose content includes the compiler, `FFLAGS`, module flag, and preprocessing options.
- Make link-library changes trigger relinking.
- Give `debug`, `release`, and `check` separate object directories by default.

### 6.2 Floating-point literal precision

Locations such as `src/parameters.f90:230-265` assign default-precision literals to `real*8` variables. The literal can be rounded to real32 before promotion to real64, contributing to 1,417 conversion warnings.

Use `iso_fortran_env, only: real64`, `real(real64)`, and literals such as `1.0_real64`. Because low bits of numerical results can change, migrate incrementally under tolerance-based physical regressions.

### 6.3 Implicit interfaces and linear-algebra ABI

MINPACK and `dgbtrf/dgbtrs` calls generated 51 implicit-interface warnings. Module or explicit interfaces would let the compiler check kind, rank, and intent.

In the current environment, the linker also warned that the executable loads `libgfortran.so.5` while system LAPACK requires `libgfortran.so.4`; `ldd` confirmed both runtimes. Although the build succeeds, this is an ABI/runtime-state risk. Production should use BLAS/LAPACK built against the same compiler runtime.

### 6.4 Hard-coded paths and launch directory

`input_read.f90:818-831` contains a developer-specific absolute fallback path for the lower-atmosphere helper. Other clones fail unless `EXHALE_ROOT` is set, and unquoted shell-command paths fail when directory names contain spaces.

`run_EXHALE.sh` treats `pwd`, rather than the script location, as the project root. Derive the directory from `${BASH_SOURCE[0]}` and change into it explicitly.

### 6.5 Undocumented Python dependency

`EXHALE_transit.py:4` imports `astropy.convolution` unconditionally, but README requirements list only numpy, scipy, matplotlib, and tkinter. Add Astropy to the documented/locked requirements, or replace that operation with a SciPy convolution and remove the dependency.

### 6.6 Fatal exit behavior

Some input failures use plain `stop`, which can appear as successful status to a shell. Fatal configuration and I/O errors should use a common error routine or nonzero `error stop`, including filename, key, and received value in the message.

## 7. Proposed test strategy

Wind-AE contains several standalone test programs, but the top-level project has no integrated `make check` or CI target, and benchmark inputs/reference data are not connected to automated assertions. A small test pyramid would provide substantial protection.

1. Unit tests
   - EOS/composition for H-only, H+He, triplet, molecules, and metals
   - SED selection around low-IP thresholds, EOF, malformed/one-row/nonmonotonic files
   - Reconstruction/flux for constant state, shock tubes, and positivity
   - Transit against a uniform slab with an analytic column
2. Configuration tests
   - Invalid combinations, invalid enum strings, and missing files must fail early with nonzero status
3. Short integrations
   - Run limited steps for tutorial_nometals, one metal, molecular, and lower-atmosphere inputs
   - Assert mass/energy behavior, abundance bounds, and absence of NaN/Inf
4. Compiler matrix
   - gfortran debug with `-fcheck=all` and FPE traps
   - optimized gfortran
   - one ifx configuration where available, especially for alias/conformance differences
5. Performance regression
   - Record steps/s and radiation/chemistry/flux timings at fixed thread count and input
   - Compare performance only after numerical tolerances pass

Rather than making all 4,440 existing warnings fatal immediately, start with a `strict-check` target, prohibit new warnings, and then reduce the conversion and interface backlog by category.

## 8. Recommended implementation roadmap

### Phase 1: Immediate stabilization

- Add final `validate_configuration()` and the molecular-combination crash regressions.
- Remove aliasing in `Apply_BC`, `ioniz_eq`, and post-processing calls.
- Correct the transit centimetre conversion and add an analytic test.
- Harden SED EOF, minimum-row, and ordering handling.

### Phase 2: Physical consistency

- Unify low-IP threshold policy for loaded and power-law spectra.
- Separate composition presence from photoactive channels.
- Add fail-fast validation for every user-facing enum.

### Phase 3: Build and numerical hygiene

- Correct configuration-stamp dependencies.
- Migrate constants incrementally to `real64`.
- Add explicit MINPACK/LAPACK interfaces and a consistent runtime stack.
- Introduce `make check` and CI.

### Phase 4: Performance and structure

- Profile and remove hot-loop array temporaries through API/layout improvements.
- Refactor shared wave-speed and chemistry logic under golden tests.
- Split the transit script into line-table-driven modules and vectorize its core loops.

## 9. Scope limitations

This review is not a full physical validation against long production convergence runs or reference observables. The valid example was run only as a 20-second smoke test. Inspected `do-subscript` warnings guarded by endpoint branches were treated as false positives. Performance findings other than the runtime-confirmed array temporaries should be reprioritized using a production profiler.

Nevertheless, the array overrun in 3.1, the unit omission in 3.6, and the stale-object behavior in 6.1 were directly reproduced or arithmetically confirmed and should be treated as immediate fixes.

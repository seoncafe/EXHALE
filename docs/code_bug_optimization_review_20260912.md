# EXHALE v1.00: Source-Based Bug and Optimization Review

Date: 2026-09-12  
Source revision: `e8eb6e6ebbbd5a64d6e3c8ba5e9184f10868a264`  
Scope: current implementation, with emphasis on recent chemistry, stationary-solver, FUV-input, and restart-mapping changes.

## 1. Executive assessment

The current code still contains demonstrable defects. The most important findings are incorrect integration at clipped FUV band edges, suppression of a valid integral when one spectral segment covers a band, and success returned by thermochemical relaxation without satisfying its convergence or admissibility contract. These affect physical inputs or the meaning of an accepted state; a small numerical effect would not make them correct.

The new restart mapper also needs stronger validation and an explicit metadata policy. Its interpolation is intentionally nonconservative, which can be acceptable for constructing an initialization seed. It must not silently represent the result as an unchanged, certified physical trajectory.

This is a targeted audit, not a claim that every module or every physical approximation has been validated. No production source was changed. No atmospheric result was regenerated.

| ID | Priority | Finding | Evidence |
| --- | --- | --- | --- |
| B1 | High | FUV integration uses the original endpoint fluxes after clipping a segment | Unchanged production routine executed with an analytic spectrum |
| B2 | High | A single segment covering the requested band is discarded | Unchanged production routine executed |
| B3 | High | Five unsuccessful thermochemical cycles still return success | Production routine body executed with controlled dependencies |
| B4 | High | Thermochemical acceptance ignores the reported element-simplex failure | Production routine body executed with a controlled ledger |
| B5 | Medium | Restart mapping accepts inconsistent source radii and implicit extrapolation | Production Python utility executed on separate manufactured cases |
| B6 | Medium | Nonconservative mapping preserves trajectory and certification metadata | Production utility executed; loader and main-program consumers inspected |

Priority describes the affected contract, not measured frequency in production runs. In particular, the B3/B4 tests establish control-flow defects, not how often a real atmospheric state reaches those conditions.

## 2. Evidence and reproducibility

I read the recent portions of [Update_EXHALE_stage2.md](Update_EXHALE_stage2.md), then inspected implementation and consumers rather than treating the log as evidence of current behavior. Main source areas inspected were:

- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`: fixed-conserved-state chemistry, carrier relaxation, attempted substeps, and checkpoint state copying.
- `src/modules/radiation/ionization_equilibrium.f90`: the chemistry ledger and its actual assignment.
- `src/modules/radiation/sed_read.f90` and `src/modules/files_IO/input_read.f90`: band integration and assignment to molecular radiation inputs.
- `src/modules/time_step/steady_newton.f90`: the recent ghost-density refresh immediately before residual assembly.
- `src/utils/map_state_to_grid.py`, `src/modules/files_IO/load_IC.f90`, and `src/EXHALE_main.f90`: mapping, restart metadata, and the physical clock.

The retained diagnostic entry point is [run_probes.py](audit_20260905/bug_optimization_review_20260912/run_probes.py). Run from the project directory:

```bash
python docs/audit_20260905/bug_optimization_review_20260912/run_probes.py
```

It extracts the current Fortran routine bodies without rewriting their statements, supplies a minimal module context, and compiles with `gfortran -O0 -fcheck=all -ffree-line-length-none`. The SED reader and integrator use their actual production bodies. Only the chemistry dependencies are controlled substitutes. The Python mapper is invoked directly as a subprocess.

The generated Fortran source, executable, module file, manufactured inputs, mapped outputs, and [results.json](audit_20260905/bug_optimization_review_20260912/results.json) remain in the same audit directory. Rerunning the script updates these diagnostic products only. The script records observations; a successful process exit does **not** mean the production defects passed a correctness test.

No full EXHALE build, full regression suite, long stationary solve, or physical-time integration was performed in this audit. No speedup was measured. These limits are important: the tests are sufficient to establish the local arithmetic and acceptance failures below, but not their integrated atmospheric impact. This review does not independently revalidate published rate coefficients or the full molecular network against papers.

## 3. Confirmed defects and proposed corrections

### B1. Incorrect trapezoid after clipping an FUV segment

**Location:** `src/modules/radiation/sed_read.f90:508`, especially the accumulation near line 574 in `sed_band_integrated_flux`.

**Physical judgment:** the requested band must receive the integral of the stated spectral interpolation over that band, not an average formed at different wavelengths.

The routine clips a spectral interval to `wa = max(w_prev,w_lo)` and `wb = min(w,w_hi)` but accumulates

```fortran
F_band = F_band + 0.5d0*(f_prev + f)*(wb - wa)
```

The mean of the original endpoint fluxes is not generally the mean over the clipped interval. For the piecewise-linear spectrum that the trapezoid assumes, the correct expression is

```text
s = (f - f_prev)/(w - w_prev)
fa = f_prev + s*(wa - w_prev)
fb = f_prev + s*(wb - w_prev)
delta_F = 0.5*(fa + fb)*(wb - wa)
```

The executed manufactured spectrum was `F_lambda = lambda`, tabulated at 900, 1000, 1100, and 1200 Angstrom. Over [950,1100] Angstrom:

| Quantity | Value |
| --- | ---: |
| Exact integral | 153750 |
| Current production routine | 152500 |
| Signed relative error | -0.813008% |

The normalization is artificial; this is an integration test, not a stellar measurement. Its purpose is to expose an error even for a linear spectrum, which the stated quadrature should integrate exactly.

**Affected path:** `input_read.f90:1574` calls the routine for LW when molecular chemistry and a loaded SED are enabled and no explicit LW flux was stated. Lines 1583–1604 also use it for B3/B4 with oxygen chemistry when those fluxes are not explicitly supplied. These become radiation inputs to molecular photochemistry. Explicit band-flux overrides bypass this particular error. The LW thin-rate prescription in `lyman_werner.f90` depends on the supplied LW flux; therefore the error is not merely a reporting discrepancy.

**Correction:** interpolate flux at each clipped boundary before integrating. Keep the wavelength and flux units unchanged. Add exact constant- and linear-spectrum tests with both edges clipped, either edge clipped, and edges coincident with tabulated nodes. Also test unequal segment widths so cancellation cannot conceal an error.

### B2. A fully covered band can be assigned zero flux

**Location:** `src/modules/radiation/sed_read.f90:586`, `if (nin .lt. 2) F_band = 0.0d0`.

Here `nin` counts contributing **segments**, not supporting tabulated points. Two points can bracket the entire requested band and define one valid linear segment. That is sufficient for an integral under the chosen interpolation model, even if it provides limited spectral resolution.

Using only the points `(900,900)` and `(1200,1200)` for the same [950,1100] band, the executed function returned **0**, whereas the exact linear-spectrum integral remains **153750**. This failure is independent of B1: even a constant positive spectrum would be discarded.

**Correction:** do not equate the number of segments with wavelength coverage. Return a separate coverage/read-status result, and integrate a single valid covering segment. Whether sparse sampling is scientifically adequate should be an explicit resolution warning or rejection policy, not conversion of a positive field into zero.

**Related robustness issue, established by source inspection:** this function returns on open failure, breaks on malformed input, and does not verify complete band coverage. If enough earlier segments contributed, a partial integral can remain a positive result that `input_read` marks as originating from the spectrum. Zero, incomplete coverage, malformed data, and a deliberately stated zero flux are physically different cases. The principal `read_sed` path may reject some bad files independently, but that is not a complete status contract for the separate FUV reread.

Recommended interface information is `(flux, coverage, status)`. Require finite, nonnegative flux and strictly increasing positive wavelengths throughout the relevant support, with an explicit policy for missing coverage. Test this both in isolation and through input initialization.

### B3. Thermochemical iteration exhaustion is reported as success

**Location:** `src/modules/lower_atmosphere/diffusive_photochemistry.f90:4880–4916`, `equilibrate_chemistry_at_fixed_conserved_state`.

The routine initializes `ok = .true.`, performs at most five chemistry/EOS cycles, and returns early if its temperature increment is below `1e-6`. If all five cycles fail this criterion, execution reaches the end without setting `ok = .false.`.

The controlled execution increased the EOS temperature by one unit at each evaluation. It returned:

```text
cycle_exhaustion_ok=T, cycles=5, last_relative_increment=0.009434
```

The final increment is far above `1e-6`, yet the routine signals success. The real chemistry was not used for this test; the test isolates the exact production acceptance logic.

**Actual caller:** `relax_photochemical_composition` calls this routine near line 5036 and marks a trial as `refused_chemistry` only when `chem_ok` is false. A false success therefore bypasses the intended chemical rejection/rollback branch. `EXHALE_main.f90:6538` invokes this carrier relaxation in the partitioned stationary solve and treats the returned composition as already chemically refreshed. Later certification may still reject the state; this audit does not claim that B3 automatically produces a falsely certified final solution.

**Physical consequence:** pressure and temperature may be algebraically consistent with the returned composition at fixed conserved energy, while chemistry, heating, cooling, and radiation-dependent rates remain evaluated at an earlier temperature. Algebraic EOS consistency is not thermochemical equilibrium.

**Correction:** initialize the convergence result to false and set it true only after all required checks succeed. Return an explicit exhaustion reason and the final normalized increment. Ensure the existing caller restores the trial state on exhaustion. Do not simply increase five to a larger unvalidated number: that postpones the same logical failure.

**Acceptance tests:** converged first cycle, converged last allowed cycle, exhaustion, nonfinite temperature, and a trial whose rejection must restore composition, thermal quantities, and background data. A subsequent test with actual molecular chemistry should measure how often exhaustion occurs and whether damping is needed.

### B4. A finite chemical state is accepted despite a failed admissibility ledger

**Location:** the same routine, especially lines 4905–4910. Compare `ionization_equilibrium.f90:426` and the actual ledger assignment near line 2732.

The current guard checks `ledger%n_nonfinite` and finite values of `f_sp`, but does not check `ledger%n_offsimplex`. The chemistry implementation sets the latter from molecular clamping and atomic handback failures. This is an outcome signal, not merely the history of an intermediate invalid iterate.

The controlled execution returned a finite composition, constant temperature, and `n_offsimplex=1`. The unchanged production routine returned:

```text
off_simplex_ok=T, cycles=1
```

**Physical judgment:** finite numbers do not establish element-budget admissibility or reaction balance. Even if chemistry supplies a finite repaired/retained state after failure, that state must not be presented as a successfully closed chemical solution on this basis.

**Correction:** consume the final admissibility signal and establish a separate reaction-residual acceptance test. Do not reject on `viol_worst` alone: the ledger documentation and code distinguish unsuccessful final outcomes from intermediate violations that can precede a valid root. Similarly, do not use a solver termination flag alone as a physical equilibrium certificate.

A small temperature increment is insufficient by itself: individual ionic or molecular populations and reaction residuals can change without a comparable change in total particle number. Also check finite, physically admissible `p`, `T`, `heat`, `cool`, and `eta`, and explicitly distinguish accepted equilibrium from an admissible but incomplete relaxation state.

This is related to B3 but requires a separate fix. Correct exhaustion handling alone does not reject the one-cycle failed-ledger example.

### B5. Restart mapping silently accepts invalid grid relationships

**Location:** `src/utils/map_state_to_grid.py:69–102`.

The tool derives interpolation coordinates exclusively from the hydro radius column. It reads the species radius column but never compares it with the hydro radii. Matching row counts are not sufficient to establish that a species row belongs to the same cell.

It also calls `numpy.interp` without testing whether target radii lie inside the source support. Targets outside the support consequently receive endpoint values rather than an explicitly selected physical extrapolation.

Separate executions established both behaviors:

- Source hydro radii `[1,2]`, species radii `[1.2,2.2]`, and target radii `[1,2]` returned exit code 0 and no error. Species values were attached to the hydro-based output coordinates. The appended center-shift comment even reported zero because it compares hydro grids only.
- With matching source radii `[1,2]` and target radii `[1,3]`, execution again succeeded. The output density at radius 3 was 5, the source endpoint density at radius 2, without an extrapolation warning.

These are small utility fixtures, not valid full EXHALE restart files. They demonstrate absent validation in the mapper; they do not demonstrate successful loading of these particular two-row fixtures.

`load_IC.f90:1686` validates output radii against the constructed grid, but cannot recover the original inconsistency after the mapper replaces both radius columns. Other loader checks may still reject a mapped product for independent reasons.

**Correction:** validate positive finite, strictly increasing coordinates; source hydro/species row counts and radius agreement; required column identities; and nonzero hydrogen-nucleus density before interpolation. Reject target physical cells outside source support unless a separately requested physical extension is implemented. Handle ghost cells under an explicit boundary convention rather than allowing their treatment to decide physical-domain validity.

The tool already documents its nonconservative nature. For initialization seeds, preserving local composition ratios can be reasonable; this review does not classify nonconservation alone as a bug. For physical trajectory transfer, remap conserved quantities using cell volumes and rebuild the thermodynamic state from the mapped composition and energy.

### B6. Mapped data retain claims that belonged to the source state

**Location:** `map_state_to_grid.py:103–115`; consumers in `load_IC.f90:1555–1573` and `EXHALE_main.f90:1681–1700`.

The mapper changes numerical rows but copies all headers except the row-count line, then appends a descriptive mapping comment. In the executed interior-remap example it preserved:

```text
# coupling mode=phys t_phys=123 certified=T
```

The loader interprets these tokens as a physical run mode, elapsed physical time, and a stationarity claim. The main program continues the physical clock from `ic_t_phys` when both modes are physical. The appended mapping comment does not establish a separate trajectory-transfer policy in those consumers.

This does **not** mean certification checks are bypassed: the stationary claim is reevaluated, and changed grid/schema metadata can cause rejection. It does mean a nonconservative interpolation tool preserves state-specific claims without deciding whether they remain valid. The risk is conditional on the rest of the restart checks succeeding, for example after a center-construction change with otherwise compatible metadata.

**Correction:** make the default output an initialization seed, explicitly invalidate certification, and establish a new physical time origin when it is used to start physical integration. Preserve source time and source certification as provenance fields rather than active claims. A physical continuation option requires a separate conservative transfer and validation contract. If the target changes grid metadata, generate accurate target metadata or refuse an unsupported conversion; copying old metadata is not a general solution.

## 4. Optimization opportunities, not measured speedups

### O1. Reduce repeated checkpoint allocation without weakening rollback

The path `relax_photochemical_composition` → `photochemical_transport_step` → `carrier_checkpoint_take` copies substantial state. The relaxation saves `f_sp`, `bg_cell`, and thermal arrays. The transport interval takes an entry checkpoint, then takes another checkpoint on every attempted substep near line 2156.

`carrier_checkpoint_take` has `intent(out)` on a derived type containing allocatable arrays. `save_carrier_module_state`, also using `intent(out)`, assigns numerous allocated radiation, composition, scaling, and residual arrays. Repeated calls therefore incur allocation/finalization and deep-copy work; this is a source-identified cost, not a measured runtime fraction.

Measure time, calls, and copied bytes separately for accepted and rejected attempts. If material, reuse allocated storage of compatible shape with explicit ownership and validity flags. Separate quantities frozen over an interval from values that a substep can mutate, after tracing all writers. Preserve the ability to restore the original allocation state as well as values. Do not remove checkpoint entries solely because their present numerical effect seems small: rollback correctness is the acceptance criterion.

### O2. Measure the chemistry work at the actual stationary-solver boundary

Carrier relaxation can perform up to five full chemistry sweeps for each kept transport trial, and the stationary residual has its own chemistry-sweep loop. Rejected relaxation trials can repeat much of this work. The inspected loop returns `n_cycles`, but the local caller does not use it to summarize closure cost or exhaustion.

Add counts and elapsed times for chemistry calls, successful closures, exhausted closures, rejected trials, and retry depth. Use those measurements to decide whether a damped thermochemical iteration or a coupled solve is justified. The first change must be B3/B4 correctness, not loosening the tolerance to make timings look better.

Reusing radiation or chemistry data is safe only when its density, composition, temperature, and column dependencies are unchanged. Sharing an evaluation record can help make this explicit. Blindly skipping a refresh is not an acceptable optimization.

### O3. Parse the numerical spectrum once for all FUV bands

`input_read` can invoke `sed_band_integrated_flux` three times: LW, B3, and B4. Each call reopens and scans the spectrum file. A single validated spectrum representation could serve these band integrals and make coverage checks consistent.

This is primarily an input-validation and maintainability improvement; it is startup work, so it should not be advertised as a major solver acceleration without measurement. Retain the distinction between the linear wavelength interpolation used here and the ionizing-bin representation in `read_sed`. Combining storage must not silently change the radiation discretization.

### O4. Do not optimize away the recent ghost-density correction

`steady_newton.f90:3481` now performs `U_to_W` and `get_species_densities` after the last boundary refresh and before `assemble_residual`. This uses the final ghost conserved state for the particle count consumed by the residual. The call exists and is on the actual residual path; the old missing-refresh defect should not be repeated as a current finding.

Its full-grid recomputation might eventually be narrowed to ghost quantities, but only after proving that interior quantities are already current and checking downstream conduction and viscosity. This audit did not perform that proof or a timing comparison. Keep the current correction until a physically equivalent implementation is demonstrated.

## 5. Recommended correction and validation order

1. Correct B1/B2 and add exact spectral quadrature tests. Verify LW/B3/B4 automatic-input paths and that explicit overrides, including zero, retain their intended precedence.
2. Correct B3/B4 as one thermochemical acceptance contract. Add controlled failure tests and rollback tests, then a bounded molecular fixture using actual chemistry. Record reaction residuals and closure increments, not only solver status.
3. Harden the mapper before creating more restart fixtures. Test matching and mismatching source grids, out-of-support targets, finite/monotonic coordinates, identity mapping, element ratios, and output metadata. Test full loader behavior separately from utility behavior.
4. Profile checkpoint copying and chemistry sweeps on a representative molecular case. Only then select a performance redesign.
5. For fixes affecting shared SED or EOS/chemistry behavior, expand validation to the dependent atomic/molecular configurations. Do not regenerate unrelated products merely to obtain a broad pass count.

A production correction should distinguish three questions: whether a candidate is physically admissible, whether the discrete equations are satisfied, and whether the algorithm made progress. Neither finite values nor a small iteration increment answers all three.

## 6. Deliverables and limitations

Added by this audit:

- This report, `docs/code_bug_optimization_review_20260912.md`.
- Diagnostic source and generated evidence under `docs/audit_20260905/bug_optimization_review_20260912/`.

No production fixes were applied, existing atmospheric products were preserved, and no Git staging or history operations were performed. The working tree was clean at the start of the audit; the new report and diagnostic directory are its additions. The recent development log's regression counts were not rerun and are not presented here as measured evidence.

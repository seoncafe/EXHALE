# Partitioned transport–wind experiments

Date: September 11, 2026  
Production source revision: `a0f4292b39e68cb83c1578a63cf4d50bd081242f`

## 1. Result and recommendation

There is experimental support for developing a partitioned approach, but the existing alternation is not a validated stationary solver. No state in the experiments below satisfied the complete selected certification contract.

The clearest positive result is that the existing H2 transport routine, on a fixed background, reduced the wind residual from `7.380744e-2` to `2.449021e-13` on the existing 500-cell grid. Once the chemical/radiative background was refreshed, that residual increased to `7.989781e-2`. This isolates a real difficulty in the outer coupling: an accurate frozen-background species solution is not necessarily close to the joint solution.

The atomic experiment also showed substantial reductions of elemental residuals after the hydrodynamic background had improved. However, each subsequent hydrodynamic update changed the elemental residuals again. Directly relaxing composition on the original, strongly nonstationary atomic background produced an invalid composition, including a `17.6857%` mass-closure error before the chemical refresh.

Recommendation: repair the subsolve acceptance and failure contracts before expanding the iteration budget or adopting the method. Then evaluate a controlled, conservative outer update against the full residual. Do not weaken the `1e-5` wind species gates to label these results successful.

## 2. What was actually executed

These are production-module experiments, not the small algebraic examples in the preceding review.

- Fixtures: `backup/regression/atomic_elem_newton` and `backup/regression/carrier_elem_newton`.
- The atomic fixture contains H/He and seven nonzero trace-element reservoirs, with element diffusion and transport. The molecular fixture is the hot-Uranus H2-carrier case. It does not test transported H+, OH, H2O, or CO.
- Each run uses copies of the fixture's inputs and saved initial-state files in a newly created directory. The atomic input's `Load IC?` setting is changed to `True` in that copy only.
- No physical time steps precede the measurement. The driver uses the production initialization, WENO3 residual, radiation/chemical elimination, stationary solver, transport relaxation, and certification routines.
- A fresh private build used GNU Fortran 16.2.0, production flags `-O3 -fopenmp`, and the OpenBLAS library under `/opt/miniconda3/lib`. Bounds checking was additionally enabled when compiling the diagnostic driver, not retroactively on all production objects.
- The hydrodynamic and coupled calls use `dtau0 = 1`, GMRES budget 40, and the fixture's `Resid tol = 1e-8`. The independent certification tolerances remain unchanged: mass `3e-12`, momentum `1e-8`, energy `1e-6`, and wind species `1e-5` for `r >= 1.2 R_p`.
- Each run has an external wall-time limit of 180 s. Time-limited runs are explicitly incomplete, not converged.

The standalone driver is [transport_wind_experiment.f90](audit_20260905/partition_experiment_20260911/transport_wind_experiment.f90). It links the current production objects except `EXHALE_main.o`. It calls the component routines explicitly so that it can inspect the stages that the production outer loop does not expose separately.

### Measurement sequence

The partitioned experiment takes the following sequence:

1. Evaluate the complete residual at the initialized state and adopt the eliminated chemical composition returned by that evaluation.
2. Remove species from the Newton registry and call the hydrodynamic solver.
3. Reevaluate the full-system state, including the species balances.
4. Hold the hydrodynamic background fixed; measure its species residual, call the existing composition relaxation, and measure that residual again before refreshing the chemical/radiative background.
5. Restore the full registry, refresh the eliminated chemistry and radiation through the production residual, and evaluate all active certification entries again.

This experiment deliberately continues after a failed hydrodynamic subsolve to inspect the effect of the next composition update. It does not convert the failure into success. The production outer loop currently exits when `jfnk_info /= 0`; consequently, this diagnostic sequence is not a claim that the unmodified production outer loop performs all these passes successfully.

The full-state measurements reassemble the hydrodynamic rows on the adopted composition and boundary, rather than combining the residual of one state with the composition of another. The difference from the preceding production residual evaluation is logged as `hydro_reassembly_difference`. It was zero at both initial states and small but sometimes nonzero later. Returned-state tables below use the explicit reassembled state.

The driver also records the chemical-map composition change. The measurement is therefore an explicit chemical-elimination/update stage, not a claim to be a read-only diagnostic of an arbitrary unrelaxed chemical composition. The raw frozen-background measurements are recorded separately for this reason.

## 3. Molecular case: small composition advances

The final confirmation used eight OpenMP threads, three outer passes, a hydrodynamic cap of 40 iterations in each pass, and the existing carrier movement setting `trust = 0.01`. It completed in `16.737055 s`. The hydro calls took 13, 7, and 7 iterations; the carrier routine took one transport step in each pass. An earlier one-thread pilot completed in `80.291291 s`; its entire `rows.tsv` table was identical to the confirmation's table in a direct file comparison.

Source: [run log](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.g44qGR/run.log), [row table](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.g44qGR/rows.tsv), and [step table](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.g44qGR/steps.tsv).

All numbers in the following table are directly executed normalized residuals. The hydro columns are maxima over the physical column; the species column is the gated wind maximum.

| State | Mass | Momentum | Energy | H2, wind |
| --- | ---: | ---: | ---: | ---: |
| Initial | `2.284006e-1` | `1.954291` | `1.870595` | `1.165916e-1` |
| After hydro 1 | `5.954113e-12` | `1.019682e-9` | `2.540981e-9` | `7.380744e-2` |
| After species 1 and refresh | `1.786211e-1` | `5.109083e-2` | `1.007056` | `6.960871e-2` |
| After hydro 2 | `8.002587e-12` | `8.779793e-10` | `3.035255e-9` | `6.967684e-2` |
| After species 2 and refresh | `1.539697e-1` | `4.754165e-2` | `1.008247` | `6.533034e-2` |
| After hydro 3 | `4.991987e-12` | `1.007949e-9` | `2.548997e-9` | `6.541660e-2` |
| After species 3 and refresh | `1.539526e-1` | `4.450620e-2` | `1.008240` | `6.089371e-2` |

The H2 residual decreases, but the complete system is never simultaneously stationary. In particular, ending the sequence immediately after a species update leaves large hydrodynamic residuals. Ending after a hydrodynamic update leaves the H2 balance unsatisfied.

The dimensional mass closure remained within approximately `1.3e-14` relative to the prescribed density, and the recorded species populations remained nonnegative. These are useful invariants, not substitutes for the missing stationary balances.

### 3.1 The hydro success flag matters

The three hydrodynamic calls returned `info = 2`, not zero. They satisfied the internal stopping condition, but the independent mass residual still exceeded `3e-12` at the returned state. Momentum and energy being much smaller than their tolerances does not waive the mass equation.

The production solver's loop-top test at `steady_newton.f90`, around lines 15039–15044, uses `steady_gates_met`. The final return at approximately lines 16048–16066 additionally checks the certified rows of the returned state. This permits an early exit from iteration followed by a correct refusal at return. The production outer loop then exits at `EXHALE_main.f90:6106`.

This is not a false certification: the return flag refuses correctly. It is an operational mismatch that must be resolved if an outer iteration is expected to continue until the joint gates are met.

### 3.2 The carrier movement setting is not an enforced maximum

Although `trust = 0.01` was supplied, the reported composition drifts were:

| Pass | Transport steps | Returned normalized drift |
| --- | ---: | ---: |
| 1 | 1 | `2.812676e-2` |
| 2 | 1 | `2.780587e-2` |
| 3 | 1 | `2.748751e-2` |

The code in `relax_photochemical_composition` tests the movement threshold after taking a step and exits if it was exceeded. It does not shorten or reject the step that exceeded the threshold. Thus decreasing this setting alone does not guarantee a smaller first composition update.

For a controlled outer method, either define this parameter honestly as a stopping trigger, or implement an actual accepted-movement limit through a reduced step/retry or a conservative outer update. The experiment did not implement such a correction.

## 4. Molecular case: a full frozen-background relaxation

To distinguish a weak transport solver from a deliberately short pass, one additional run used the same molecular initial state and first hydro solve, but set `trust = 1e6`. This disables that movement stopping trigger for the measured trajectory; it does not change the transport equations, their internal tolerances, the 400-step maximum, or the final certification gates.

Source: [run log](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.TZ326I/run.log), [row table](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.TZ326I/rows.tsv), and [provenance](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_partition.TZ326I/provenance.log).

| Measurement | Whole-column H2 residual | Wind H2 residual |
| --- | ---: | ---: |
| Before frozen-background relaxation | `4.927495e-1` | `7.380744e-2` |
| Immediately after relaxation, background still frozen | `2.346741e-12` | `2.449021e-13` |
| After chemical/radiative refresh | `9.999997e-1` | `7.989781e-2` |

The relaxation took 38 transport steps and `58.169630 s`, with a returned composition drift of `0.9787245`. The entire run took `107.943121 s`. After the refresh, the mass, momentum, and energy residuals were approximately `1.302329`, `1.776463`, and `1.000015`, respectively. The state remained finite and nonnegative; its mass closure was approximately `7.56e-15`.

This is positive evidence that the fixed-background H2 transport problem can be solved accurately. It also directly demonstrates that a residual near `1e-13` is attainable for this discrete species subproblem on the existing grid; a quoted spatial truncation error is not its algebraic floor.

It is not a proof that the joint molecular atmosphere has been solved. The subsequent background refresh changes the species equation. The experiment identifies that combined feedback but does not separate the individual contributions from radiation, temperature, particle counts, and eliminated chemical stages. That attribution would require an additional controlled experiment.

The large frozen-background update is therefore useful as a component and diagnostic. Handing the whole update directly to the next block is not demonstrated to be a convergent outer method.

## 5. Atomic case

### 5.1 Two completed outer passes

The final eight-thread run with a 40-iteration hydro cap and element under-relaxation `omega = 0.5` completed two outer passes in `135.224967 s`. An earlier three-pass pilot completed the same first two passes before its 180 s limit interrupted the third hydrodynamic solve. The table below uses the completed two-pass confirmation, not an interrupted third-pass state.

Source: [row table](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_partition.ouOtBg/rows.tsv), [run log](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_partition.ouOtBg/run.log), and [exit status](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_partition.ouOtBg/exit_status.log).

| State | Mass | Momentum | Energy | Worst elemental wind row |
| --- | ---: | ---: | ---: | ---: |
| Initial | `1.508647` | `3.348139e-1` | `1.892405` | `5.839988e-1` |
| After hydro 1 | `2.314430e-9` | `2.830229e-10` | `2.399068e-8` | `2.612907e-3` |
| After species 1 and refresh | `2.314430e-9` | `2.830229e-10` | `6.310148e-1` | `2.781845e-4` |
| After hydro 2 | `9.109392e-10` | `9.450913e-13` | `1.842713e-8` | `3.811843e-3` |
| After species 2 and refresh | `1.151322e-9` | `9.450913e-13` | `5.642080e-1` | `7.565920e-5` |

Each element relaxation took 30 transport steps, approximately `0.11 s`. The much larger recorded costs were in the hydrodynamic calls, approximately `63.83 s` and `70.74 s`. These are run timings, not isolated throughput benchmarks.

The first hydro update made the composition relaxation much better behaved than relaxation on the original state. At the same time, the second hydro update increased the worst wind elemental residual from `2.78e-4` to `3.81e-3` before transport reduced it again. This is another direct observation of outer feedback, not evidence that the coupling has been removed.

The hydro calls returned `info = 1`; neither had delivered a certified hydrodynamic subsolution within its 40 iterations. Consequently, these are inexact exploratory block updates, not a test of the convergence rate of exact block Gauss–Seidel.

The existing element relaxation also damps helium but leaves the independently relaxed trace metals at their updated values. Its `omega` is not a uniform damping parameter for the entire composition vector. That distinction matters when designing a common outer update.

### 5.2 A deliberately insufficient hydro warm-up does not repair the problem

A separate eight-thread run with only 10 hydrodynamic iterations completed one outer pass in `16.123621 s`. Before composition relaxation, the mass and energy residuals were still approximately `1.45` and `1.57`. After relaxation, the worst wind elemental residual increased to `1.0`, although mass closure remained near `1.44e-14`.

Source: [row table](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_partition.gVmnY3/rows.tsv) and [log](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_partition.gVmnY3/run.log).

Thus more nearly stationary hydrodynamic backgrounds matter physically and numerically. A small number of arbitrary hydro iterations followed by a very long composition relaxation is not a reliable replacement for a coupled solve.

### 5.3 Relaxation on the original atomic background exposes a failure-contract defect

The species-only stress test used the original saved atomic background, no hydro correction, and `omega = 1`. The background initially had large mass and energy residuals; it was intentionally not a certified steady wind.

The element relaxation returned after 183 transport steps. Its log reported an inner composition Newton failure after 30 iterations, with relative residual approximately `0.1385`. The raw returned species were finite and nonnegative, but the production `calc_rho` routine measured

\[
\max_j\frac{|\rho_{\mathrm{species},j}-\rho_j|}{|\rho_j|}
=0.1768571764708571.
\]

This was measured before the next chemical refresh. That refresh subsequently produced a nonfinite state, and the diagnostic stopped with a nonzero exit status.

Source: [instrumented log](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_species.MlxOaF/run.log), [raw returned composition](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_species.MlxOaF/frozen_after_01.state), and [exit status](audit_20260905/partition_experiment_20260911/runs/atomic_elem_newton_species.MlxOaF/exit_status.log). Earlier repetitions without the raw mass diagnostic are retained in the archive; they are not independent physical cases.

The source explains why this outcome can be passed onward:

- `solve_mass_fraction` can accept its last halved trial even when the residual was not reduced. It reports failure through a diagnostic and warning, not a returned success flag.
- `element_diffusion_step` clips the helium fraction and continues the composition projection and other transport work.
- `relax_element_composition` grows its time step and uses composition movement as an exit test; it does not require a valid stationary residual or a successful status from every internal solve.

This experiment does not locate the first instruction that introduced the mass discrepancy. It does establish that the discrepancy exists at the transport routine's return, before chemistry, and that finite/nonnegative checks alone would miss it.

There is also a physical compatibility issue. With density held fixed, a strongly nonzero divergence of the prescribed mass flux means that the assumed background is not a stationary mass solution. A very long conservative composition advance on that background is not guaranteed to admit a valid normalized composition. The stress test does not prove that the same failure occurs on a certified wind or in the controlled physical-time driver.

Required handling is explicit rejection or a bounded, validated update—not silently treating the warning, clipping, or a small subsequent drift as convergence. This should precede any production use of an aggressive frozen-background atomic relaxation.

## 6. Coupled control and comparison limits

The eight-thread molecular coupled control completed 40 outer iterations using the same initialized configuration and production equations. Its returned-state residuals were:

| Mass | Momentum | Energy | H2, wind | H2, whole column |
| ---: | ---: | ---: | ---: | ---: |
| `4.723616e-3` | `1.322218e-1` | `1.938361e-1` | `7.343477e-2` | `4.790024e-1` |

It returned `info = 1`, was not certified, and took `107.327893 s` including initialization and diagnostic work. Source: [control row table](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_coupled.VaAVVl/rows.tsv) and [log](audit_20260905/partition_experiment_20260911/runs/carrier_elem_newton_coupled.VaAVVl/run.log).

The final molecular partitioned run and this control use the same version-3 executable and eight-thread setting. The comparison is therefore not based solely on the earlier pilot with different instrumentation:

| Measured state | Mass | Momentum | Energy | H2, wind |
| --- | ---: | ---: | ---: | ---: |
| Coupled control, returned state | `4.723616e-3` | `1.322218e-1` | `1.938361e-1` | `7.343477e-2` |
| Partitioned, after the third hydro solve | `4.991987e-12` | `1.007949e-9` | `2.548997e-9` | `6.541660e-2` |
| Partitioned, after the following species update and refresh | `1.539526e-1` | `4.450620e-2` | `1.008240` | `6.089371e-2` |

The intermediate hydro-stage state improves every listed residual relative to this bounded coupled control, but still fails the mass and H2 gates. The subsequent species update reduces H2 further while spoiling hydrodynamic stationarity. This is promising progress, not convergence of the complete system. The prescribed partitioned experiment took `16.74 s` versus `107.33 s` for the control, but neither reached the same successful endpoint, so these timings are not a demonstrated speedup to a certified solution.

The eight-thread atomic coupled control did not return a solution within 180 s. Its log reaches iteration 26; this is not a measurement of the best returned state at the requested 40-iteration cap. One-thread coupled controls also reached the wall-time limit. These incomplete runs remain in the archive and are not used to manufacture a speedup ratio or a convergence claim.

Some diagnostic jobs overlapped, and the pilot runs used different thread counts and instrumentation versions. Their wall times document the cost actually incurred; they do not constitute a controlled efficiency benchmark. A lower H2 residual at one partitioned stage also does not prove a better complete solution when that stage has larger hydrodynamic residuals.

## 7. Defects in the current implementation

The findings below apply to the production revision identified at the beginning of this report. They distinguish confirmed implementation defects from an observed invalid result whose first corrupting operation has not been located, and from nonconvergence that is not itself evidence of a coding error. No production correction has been implemented as part of this experiment.

| ID | Finding | Classification | Priority |
| --- | --- | --- | --- |
| D1 | An unsuccessful element composition solve can return an unchecked update to the relaxation caller | Confirmed failure-handling defect; invalid returned composition reproduced | Highest: prevent invalid state propagation |
| D2 | Carrier `trust` does not enforce the advertised maximum composition movement | Confirmed update-acceptance defect; overshoot reproduced | High: enforce the intended movement contract |
| D3 | Newton iteration termination is weaker than final returned-state certification | Confirmed termination inconsistency; premature stop followed by rejection reproduced | High: align stopping and acceptance |
| C1 | Element `omega` damps helium, not the full transported element vector | Confirmed implementation limitation, explicitly stated locally; not independently established as a conservation defect | Address when defining outer composition control |

### 7.1 D1: failed element solves can propagate invalid composition

**Physical judgment.** A returned composition that does not reconstruct the prescribed mass density is inadmissible for this fixed-density calculation. Finite and nonnegative species populations are necessary but insufficient. Clipping a failed nonlinear iterate into the helium bounds does not establish that the discrete transport equation or the complete mass constraint is satisfied.

**Implementation evidence.** In [binary_element_diffusion.f90](../src/modules/functions/binary_element_diffusion.f90):

- `solve_mass_fraction`, lines 1990–2005, leaves its line search when either the residual decreases or the last allowed halving is reached. It then assigns `Xhe = Xtry` unconditionally. Reaching the halving limit therefore does not reject an unsuccessful trial.
- Lines 2020–2026 print a warning when the Newton iteration limit is exhausted above the warning threshold. The routine has no explicit success/failure output for its caller.
- `element_diffusion_step`, lines 616–636, calls that routine, clips helium into `[0,1]`, and continues. Recording the excursions does not itself reject this update in the relaxation path tested here.
- `relax_element_composition`, lines 970–982, repeatedly calls the step routine, tests composition movement, and increases the step size. It does not use an inner success result to decide whether to retain or discard the update.

**Reproduction and impact.** The species-only atomic stress case in Section 5.3 returned after 183 transport steps despite a logged Newton residual of approximately `0.1385` after 30 Newton iterations. The raw returned composition was finite and nonnegative, but its maximum relative mass-closure error was `0.1768571764708571`. The next chemical refresh produced a nonfinite state. The diagnostic terminated with exit status 1; no existing simulation result was overwritten.

**Boundary of the conclusion.** The background in this test was strongly nonstationary. The experiment does not prove failure on an already certified wind, nor does it establish the same behavior in the physical-time driver. The diagnostic deliberately called the relaxation on a background that the production outer loop would reject through its hydro return flag. The confirmed defect is that the called relaxation can hand back an invalid update without a usable failure result. The exact operation that first creates the mass discrepancy remains unlocated; the measurement does not justify blaming clipping or `project_elements` alone.

**Proposed correction.** Introduce an explicit solve result and propagate it through `solve_mass_fraction`, `element_diffusion_step`, and `relax_element_composition`. Reject a line search that exhausts its trials without satisfying its acceptance rule. Preserve the last valid composition and restore it after a rejected update; reduce the step and retry within a bounded budget, or return failure. Validate finite values, physical bounds, species-derived mass density, applicable element constraints, and the equation actually being solved before accepting an update. A converged implicit transport step and a converged stationary relaxation are different outcomes and must be reported separately. A small composition change must not stand in for a stationary residual test.

**Focused validation after correction.** Repeat the archived atomic stress case and require either a valid update or an explicit failure with the entry state restored. A downstream chemical refresh must not receive the reproduced invalid composition. Add a deliberately exhausted line-search case to verify rejection. Trace mass closure immediately after the nonlinear solve and after each subsequent composition transformation to locate the first invalid stage. Also repeat the nearly stationary atomic case so that failure handling is tested alongside a valid transport advance. Because the step routine is shared, inspect and test the affected physical-time callers when changing its interface; the existing stationary experiment does not cover those callers.

### 7.2 D2: the carrier movement limit is checked only after acceptance

**Numerical judgment.** An advertised maximum composition movement must constrain the state returned to the caller. Exiting after an excessive step has already been applied does not enforce that maximum. It also prevents an outer controller from reliably reducing the update simply by decreasing the setting.

**Implementation evidence.** The header of `relax_photochemical_composition` in [diffusive_photochemistry.f90](../src/modules/lower_atmosphere/diffusive_photochemistry.f90), lines 4718–4721, describes an upper bound. The loop calls `photochemical_transport_step` before calculating the displacement; line 4831 then exits when `drift/x_ref > trust`. There is no restoration or reduced-step retry for the step that exceeded the threshold. The stopping test normalizes by the current maximum H2 mixing ratio, whereas the final reported drift uses the larger of the entry and current maxima. A corrected contract should use one explicitly defined measure for enforcement and reporting.

**Reproduction and impact.** In the completed molecular partitioned confirmation, `trust = 0.01` returned drifts of `0.02812676`, `0.02780587`, and `0.02748751`. Each pass took only one transport step. The overshoot therefore cannot be prevented merely by ending the loop sooner. This finding concerns update control, not an established error in the molecular reaction rates or transport flux formula.

**Proposed correction.** For a true bound, preserve the state before each trial, compute the displacement relative to the outer entry state, and reject and retry any step that would exceed the bound, using a consistent normalization. A shorter step must retain the physical conservation constraints. If the intended behavior is only a stopping trigger, rename and document that meaning and remove all upper-bound claims at callers; that clarification alone would not provide the controlled outer update recommended here.

**Focused validation after correction.** Repeat the first molecular pass with `trust = 0.01` and at least one smaller value. Require the returned movement to satisfy the selected bound within an explicitly stated numerical tolerance, while preserving mass closure, finite populations, and nonnegativity. Include a first-step overshoot and a later cumulative overshoot. Evaluate the full refreshed residual afterward; satisfying the movement limit alone does not certify the atmosphere.

### 7.3 D3: iteration stopping and final certification use different predicates

**Numerical judgment.** The solver must not stop for success while a required equation still fails its acceptance threshold. Refusing that state at return is correct, but it does not repair the premature termination that prevented further iterations.

**Implementation evidence.** In [steady_newton.f90](../src/modules/time_step/steady_newton.f90), lines 15039–15044, the loop exits when `steady_gates_met` is true. At lines 16048–16066, final acceptance additionally requires the certified rows of the returned state and agreement with the certification measure. An earlier success is changed to `info = 2` when those additional checks fail. The production outer loop in [EXHALE_main.f90](../src/EXHALE_main.f90), line 6106, exits on a nonzero `jfnk_info`.

**Reproduction and impact.** All three hydro calls in the completed molecular partitioned case stopped with `info = 2`. Their independently measured mass residuals were approximately `5.95e-12`, `8.00e-12`, and `4.99e-12`, above the required `3e-12`. The stricter final check therefore acted correctly. The diagnostic continued only to measure the following composition update; the production outer loop would not perform that continuation.

**Proposed correction.** Make the iteration success test require the same relevant returned-state conditions as final acceptance, evaluated on a consistent state. If chemical elimination changes the state at return, recheck that actual state before declaring success. Continue within the remaining budget when the internal test passes but a required final row fails; if the budget is exhausted, report failure. An outer algorithm may deliberately accept an inexact block update only through a separately defined, validated acceptance rule. Do not reinterpret `info = 2` as success, weaken the mass tolerance, or remove the final certification check.

**Focused validation after correction.** Reproduce the molecular first hydro call. When its mass row remains above `3e-12`, require continued iteration or an explicit budget/failure exit rather than a success-triggered early stop. Check that every `info = 0` returned state satisfies all selected hydro gates and that the production caller preserves failure semantics. Extending the budget alone is not evidence that the stopping predicates have been aligned.

### 7.4 Findings that should not be mislabeled as additional confirmed errors

- **Outer molecular feedback:** reducing the frozen H2 residual to `2.45e-13` and then observing approximately `7.99e-2` after chemical/radiative refresh demonstrates changed background dependence. It is not, by itself, proof that the frozen transport solver or the chemistry equations are incorrect.
- **Helium-only damping:** lines 988–1008 of `binary_element_diffusion.f90` explicitly limit `omega` to helium. An outer controller therefore cannot assume that reducing `omega` reduces every elemental update. A common conservative composition control is a development requirement if that behavior is intended, but the asymmetry alone does not establish another mass-conservation error.
- **Incomplete or unsuccessful solves:** `info = 1`, an external timeout, and a residual above tolerance establish that the requested solve was not certified. They do not independently identify an erroneous physical equation.
- **Unlocated origin of mass loss or gain:** the atomic stress measurement establishes a mass discrepancy at the routine boundary, not the location or unique cause of its first occurrence. That distinction must remain explicit until intermediate-state instrumentation identifies it.

## 8. Changes warranted by the experiment

No production changes were made. The following are proposed follow-up changes, in priority order:

1. **Give the element relaxation a real failure result (D1).** Check inner-solve success, finite values, composition bounds, mass closure, and the stationary residual. Preserve the entry state so a rejected update can be discarded. Do not advance chemistry on the invalid returned composition reproduced here.
2. **Define a controlled outer composition update (D2 and C1).** The carrier setting currently checks movement after the step, and the element damping currently acts only on helium. A new outer contract must state which independent transported quantities it controls and enforce conservation and admissibility for their joint update.
3. **Align subsolve termination and returned-state acceptance (D3).** Continue solving when an independent hydro row still fails, or return an explicit failure that an outer controller can handle. Do not label `info = 2` as success or bypass it silently.
4. **Judge every accepted outer state on the full equations.** Record the raw frozen-background gain, but accept or reject the outer state only after refreshing all quantities on which its residual depends.
5. **Keep both partitioned acceleration and coupled preconditioning available.** The accurately solved frozen H2 block is potentially useful in either role. These experiments do not establish that one architecture is intrinsically superior.

The next useful convergence experiment would compare accepted outer updates after these contracts are defined, using a fixed executable and thread count, common initial states, and a full-system stopping predicate. Merely increasing the present number of outer passes would mix known interface defects with the question of the method's convergence.

## 9. Validation scope and retained data

The experiment exercised the current molecular and atomic operator paths, their callers through the standalone driver, and the downstream residual/certification consumers. It measured normalized stationary residuals, local mass closure, finite/nonnegative populations, internal warnings, solver return flags, and run costs. Later instrumentation also saves the numerical mass flux at the physical faces.

It did not independently validate global elemental flux budgets, an integrated energy budget, physical-time stability, spatial or temporal convergence, the transported-ion configuration, or the oxygen/carbon molecular carriers. No claim that these unperformed checks passed is implied. No existing output, initial condition, fixture, production source, or production executable was overwritten. No full regression suite or unrelated data generation was run.

All diagnostic code and run products are retained under [partition_experiment_20260911](audit_20260905/partition_experiment_20260911/README.md), inside the previously requested `docs/audit_20260905/` directory. The private build is retained at `/tmp/exhale_partition_20260911.ZxXJEB` and can be rebuilt from the archived scripts if that temporary directory is later removed.

The main conclusion is narrower and more useful than “partitioning works” or “partitioning fails”: an accurate frozen-background molecular transport solution is available, but the outer background update can undo that accuracy; the atomic path additionally needs robust handling of incompatible backgrounds and failed internal solves. Those are concrete development targets, not a justification for accepting an unsatisfied physical balance.

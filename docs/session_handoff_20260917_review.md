# Review of the September 17 handoff, implementation, and blocked work

Date: 2026-09-17

## 1. Executive judgment

The handoff records substantial useful work, but several of its strongest conclusions are not supported by the measurements as interpreted. Two diagnostic errors can be demonstrated without another stationary solve. There is also a genuine inconsistency in transport geometry, and important molecular approximations are described as physical bounds without sufficient justification.

The most important corrections are:

1. **The reported base-flux factors of 32.38 and 336.55 are PLM evaluations of WENO3 stationary states.** On the same states, the production WENO3 operator gives essentially constant face mass flux. These measurements do not demonstrate a current certification loophole.
2. **The old 0.03 XUV state is a separate case.** Its stored certification claim is not reproduced by the current executable. The current evaluation explicitly rejects its mass, momentum, and energy rows. An old `certified=T` header must not be combined with a new residual measurement and presented as a state that passes today's gate.
3. **The reported `steady_selfconsistent_residual` failure joins records from different stationary solves.** Its ratio of `2.83103e-08` is reproducible as a parser error. Matching records within the last solve gives `0.842987`, inside the test's factor-of-two interval.
4. **Diffusion and advection use different spherical cell volumes.** This is an implementation-level conservation inconsistency, already recognized by L12b but still present in both elemental and molecular-carrier transport. Correcting it is necessary; preserving existing numerical results is not a reason to retain it.
5. **The molecular changes are not all established corrections with known one-sided errors.** The helium/argon substitution, universal H(n=2) recipient, scalar quenching bound, and H2 fragment energy require a more precise validity contract and uncertainty treatment.
6. **The coupled solver deserves a bounded linear-algebra investigation, not another unrestricted campaign.** Poor preconditioning is plausible, but the reported mean Jacobian errors do not exclude inaccurate directional products, nonlinear elimination noise, or omitted long-range coupling.

I recommend repairing diagnostic fidelity first, then resolving transport conservation and the molecular validity contract, and only then using new molecular runs as reference results. Do not weaken certification tolerances to make the current stalled cases pass.

## 2. Scope, provenance, and what was actually executed

### 2.1 Implementation inspected

The review followed the relevant paths through:

- `src/EXHALE_main.f90`: stationary entry, reconstruction selection, evaluation, and final certification.
- `src/modules/states/{base_boundary,Apply_BC,Reconstruction}.f90` and `src/modules/init/init.f90`.
- `src/modules/functions/binary_element_diffusion.f90` and `src/modules/lower_atmosphere/diffusive_photochemistry.f90`.
- `src/modules/time_step/steady_newton.f90`, including the trust-region call path and GMRES restart implementation.
- `src/modules/lower_atmosphere/{molecular_reaction_heat,h2_vibrational_relaxation,mol_rates}.f90` and the excited-hydrogen recipient path.
- The base-continuity diagnostic, the stationary log-reading test, associated L7g/L7h, L12b, L15, L22, L25, and L26 records, and archived run logs.

The inspected repository HEAD was `3c73905ca8a7fe2af92a5c2b014c225b2feedff3`. The executable MD5 was `59bfdb3fc4d0104fc2e9c3734596d2f6`, matching the handoff. `md5sum --check --quiet LHS1140b/models/BINARY_MANIFEST_59bfdb3fc4d0.txt` completed successfully: the manifest's source checksums match the inspected source. The executable's historical build provenance and the current repository HEAD are different identifiers, not interchangeable evidence.

`make -q` requested a rebuild; inspecting the proposed commands showed a build-stamp rebuild and its dependencies. No production rebuild was performed. The numerical probe linked the existing production objects, and the separate stationary evaluation used the named executable directly. Source-manifest agreement is the source identity check used here, not a claim that a fresh build was tested.

### 2.2 New, bounded checks

The following were executed during this review:

- A small Fortran program linked against the existing production objects, comparing PLM and WENO3 on four copied restart states, with and without refreshing the boundary ghosts.
- The existing `steady_selfconsistent_residual` log reader against `.L22/i3_alt/run.log`, followed by an independent record-scope check of that same log.
- `Restart intent: stationary evaluate` on a scratch copy of `.L14/x003_HeH2.13/output/`, using the current executable. It took **zero attempted and zero accepted steps**, rejected the state, wrote products only in the scratch directory, and exited with status 2.
- Source-checksum verification and direct reading of published papers already in `references/`, extracted with `pdftotext -layout`.

The first compilation of the small probe exposed an invalid reference to a main-program variable; that reference was removed before measurement. Two initial scratch runs could not find relative spectrum paths; only those scratch input paths were changed to the existing absolute spectrum paths. The retained flux logs contain the successful runs.

No new stationary campaign, full regression suite, molecular trajectory, parallel scaling study, or catalog regeneration was run. The review changes no production source, inputs, reference results, or catalog state. Newly written diagnostic programs and generated logs are retained under `docs/audit_20260905/`; see section 11.

## 3. Confirmed diagnostic defects and the corrected interpretation

### 3.1 Reconstruction mismatch in the base-continuity probe (high priority)

**Implementation evidence.** `src/tests/grid_and_gates/base_boundary_continuity_probe.f90:109` initializes through `input_read` and `init`. Its `local_flux_candidates` routine, around lines 636–660, converts `Wq`, calls `Reconstruct`, and calls `RK_rhs`. It does not select the stationary reconstruction. `input_read.f90:374–386` dispatches the input's reconstruction scheme, and these inputs request PLM. In contrast, the stationary and residual routes in `EXHALE_main.f90` explicitly select WENO3; examples are lines 1231, 1531, and 5464. `Reconstruction.f90:77` dispatches on `rec_method`.

Calling the same initialization routines is therefore **not** equivalent to calling the same stationary operator.

**New measurements.** Let `F0` be `wind_window_mass_flux`'s mean of the cell-centered quantity `rho*v*r^2`. Values below are dimensionless ratios; the WENO3 columns include `Apply_BC` immediately before reconstruction.

| Copied state | PLM base face, original probe convention | WENO3 base face, refreshed ghosts | WENO3 minimum and maximum over all faces |
|---|---:|---:|---:|
| `.L26/fid_resolve/output/` | 32.3762799477 | 0.999984373657 | 0.999984373624 to 0.999984374845 |
| Catalog `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/output/*_IC.txt` | 336.551905978 | 0.999974828179 | 0.999974777751 to 0.999974828179 |
| Old `.L14/x003_HeH2.13/output/` | -404.318081737 | 607.887544875 | 0.999971914706 to 612.748500726 |
| Old `.L14/x002_HeH2.13/output/` transient | -584.452904971 | 855.800267543 | 0.999966257403 to 886.655843599 |

The first two PLM numbers reproduce L26 R3 to the reported precision. Their WENO3 face-flux ranges have widths about `1.22e-9 F0` and `5.04e-8 F0`, respectively. A nearly constant face flux need not equal the mean **cell-centered** product exactly; the roughly `1.6e-5` and `2.5e-5` offsets from unity are not mass leakage.

**Boundary refresh is a second diagnostic issue.** `init.f90:199–200` updates the conservative array `u` through `Apply_BC`, but does not copy its updated ghosts back into the returned primitive array `W`. The original diagnostic subsequently executes `W_to_U(Wq,u)`, restoring the file's primitive ghost values. Calling `base_boundary_states` into local output arrays is not a replacement for installing the boundary state through `Apply_BC`. `Apply_BC.f90:45–59` documents the saved face-state contract consumed by reconstruction. This matters particularly for the two old `.L14` states, whose refreshed-ghost results differ substantially from the original probe convention.

**Required repair.** Give diagnostics an explicit operator identity and a single production evaluation entry point. It should install the intended composition/caloric state, choose the same reconstruction and flux as the stationary route, refresh ghosts, reconstruct, and evaluate fluxes. Print these choices, the input-state identity, and the certification identity beside the measurement. Add a test that a PLM-requesting restart input still evaluates the stationary WENO3 operator when the stationary diagnostic is requested.

This does not establish that PLM is defective: a state solving one discretization is not required to solve another. It establishes that PLM fluxes cannot be used to accuse the WENO3 certification of accepting a mass imbalance.

### 3.2 The old 0.03 claim is rejected today, not silently certified

The large WENO3 values in the last two rows must not be hidden by the first two successes. To resolve the certified 0.03 case, I ran the production stationary-evaluation route on the same copied pair, without taking a solve or a time step.

The current executable reported:

| Current certification entry | Measured maximum | Cell | Tolerance |
|---|---:|---:|---:|
| Hydrodynamic mass | `9.984e-01` | 2 | `3.2e-10` |
| Hydrodynamic momentum | `1.955e-03` | 1 | `1.0e-08` |
| Hydrodynamic energy | `1.000e+00` | 2 | `1.0e-06` |

It printed `original claim: NOT REPRODUCED`, wrote `certified=F cert_reason=failing_entries`, and exited with status 2. The raw log is retained.

The copied hydro restart had MD5 `6b21753dcb264b47bb4a322ea9dc00b5`, identical to the archived source. This is not a comparison against a newly solved replacement.

**Conclusion.** Section 4.3 and the corresponding decision in section 6.2 of the handoff need correction. The first two alleged failures use the wrong reconstruction; the old 0.03 state actually fails the current gate; the 0.02 state was already a transient. None of these observations establishes the claimed current certification loophole.

The historical reason why the old 0.03 file carries its original claim has not been isolated here. It requires its original executable, boundary configuration, and output sequence. The review does not attribute that history solely to L26 or prove that the older executable certified the same refreshed state. The current boundary/reload incompatibility is real and localized to the base; its historical cause remains open.

### 3.3 Stationary log reader joins different solves (confirmed)

`src/tests/steady_selfconsistent_residual/run.sh:65–75` selects the last certification and completion record. Around lines 122–129 it independently selects the last `returning best iterate` record **anywhere in the entire log**, using an iteration record only if no such line exists.

On `.L22/i3_alt/run.log`, the records actually selected are:

| Selection | Log line | Residual |
|---|---:|---:|
| Earlier best-iterate return | 368 | `4.267e-01` |
| Last iteration of the final solve | 2232 | `1.433e-08` |
| Final completion | 2284 | `1.208e-08` |

The existing test computes `1.208e-08 / 4.267e-01 = 2.831028826e-08`. Selecting within the last solve gives `1.208e-08 / 1.433e-08 = 0.8429867411`. The latter passes the existing interval `[0.5,2]`.

The test's accompanying assertion that composition elimination is lagged is therefore not warranted by this failure. Correct record selection does not prove every elimination path is correct; it removes this particular false accusation.

**Repair.** Attach a solve identifier to iteration, best-state, certification, and completion records, and parse one completed solve at a time. Until that exists, delimit the latest solve and never search earlier records for its best iterate. Test multiple outer passes, a best-state restore only in an earlier solve, no restore in the final solve, and a truncated final solve. Report incomplete evidence instead of combining unmatched records. Structured scalar records would be more reliable than several independent `grep | tail` selections.

## 4. Actual transport inconsistency and how to unblock L12b

### 4.1 Different volume factors remain in production

The physically conservative spherical finite-volume divergence is

```text
D(F)_j = [r_+^2 F_+ - r_-^2 F_-] / V_j,
V_j = (r_+^3 - r_-^3)/3.
```

All flux contributions to the same conserved quantity must use the same `V_j`.

The current implementation uses:

- Exact spherical volume in `binary_element_diffusion.f90:2270–2279`, `element_advective_face_coefficients`.
- `1/[r_j^2 (r_+ - r_-)]` in the elemental diffusion matrix and residual, lines 2395 and 2589.
- The same approximate volume in the trace-element residual and matrix, lines 3204–3205 and 3337–3338.
- Exact spherical volume in carrier advection, `diffusive_photochemistry.f90:3628–3632` and 3663–3665.
- The approximate volume in carrier diffusion, its residual, and matrix, lines 4600 and 4864.

The `max(width,1 cm)` safeguard is an additional modification if that floor activates; it is not an exact geometric identity.

For an arithmetic midpoint, the difference can be seen directly:

```text
V_j = r_j^2 * dr_j + dr_j^3/12.
```

For other center conventions, compare the actual `V_j` and `r_j^2*dr_j` directly. When diffusion residuals are multiplied by the actual volume, their face terms acquire the factor `V_j/(r_j^2*dr_j)`. If that factor differs between adjacent cells, an internal face does not cancel in the column sum. Advection does cancel under those weights. There is consequently no single physical volume weighting under which the present mixed operator is the divergence of the claimed total flux.

This is not a claim of literal creation of nuclei by a reaction or of a newly measured effect size. It is a discrete conservation inconsistency established by the implemented coefficients. Its magnitude on a particular case remains to be measured after defining the correct operator.

### 4.2 Recommended implementation shape

The L12b diagnosis is justified, but lack of an exposed total flux is a software dependency, not a data barrier.

1. Define physical face areas and exact cell volumes once. Use a cancellation-resistant expression such as `dr*(r_+^2+r_+*r_-+r_-^2)/3` where appropriate.
2. Make the elemental transport path return the actual advective and diffusive **face fluxes**, with explicit units and species/element weights. `element_face_flux` already exists at lines 2174–2209 and is invoked in production; the missing piece is an exposed, consistently assembled total flux, not an absent diffusion routine.
3. Use those same fluxes and geometry in residuals, time updates, matrices, and diagnostics. Changing a diagnostic denominator alone would not repair the operator.
4. Define the ionization-stage flux from that nucleus flux:

   ```text
   F_k = x_k_face * N_element_face
         - n_element_face * K_face * (x_k_right - x_k_left)/distance.
   ```

   Preserve `sum_k x_k_face = 1` at the reconstructed face, not just at cell centers. Then `sum_k F_k = N_element_face` is a discrete identity. Reusing the full elemental eddy term inside every stage and adding another abundance-gradient eddy term would count part of the transport twice.
5. State the boundary nucleus flux and stage fractions consistently with the reservoir. At an outflow, do not independently impose incompatible stage and elemental fluxes.

**Acceptance tests.** Sum the exact-volume-weighted residuals and compare against boundary fluxes and integrated sources; test nonuniform grids, a zero-source constant-flux state, a static equilibrium, and a manufactured smooth profile. Check the stage-flux sum to floating-point roundoff, and repeat after positivity limiting and boundary reconstruction. Update preconditioners together with residual coefficients and verify directional derivatives afterward.

The L12b paper-derived charge-stage drift concern remains a separate physical limitation. A sum-preserving stage flux does not validate the assumption that ions and neutrals move together; collision/charge-exchange times and ambipolar drift must be checked for the intended altitude range.

## 5. Molecular physics: established improvements versus unsupported bounds

### 5.1 Published material checked

I read the relevant methods and results in the local final journal PDFs, including Takagi (2002), Giusti-Suzor et al. (1983), Guberman (1994), Kokoouline et al. (2001), Strasser et al. (2001), Cohen and Westberg (1983), and Paolini et al. (2011). Lique (2015) was also inspected for the scope of its collision calculations. These are not conclusions drawn solely from the code comments.

An ADS API metadata request was attempted first, as required by the repository instructions, but failed because the host could not be resolved. Publisher DOI access through the web tool also failed. The literature conclusions below therefore rely on the local published PDFs, not on a successful new literature search. This review does not claim an exhaustive search for additional rate tables.

### 5.2 R5/R16: subtracting one n=2 excitation is a restricted model

`molecular_reaction_heat.f90:522–550` subtracts `10.19883 eV` from the ground-reference enthalpy for every R5 and R16 event. The associated excited-hydrogen source counts one n=2 atom for every event. The recipient implementation exists and is invoked; it is not missing.

The physical issue is the generality of that prescription:

- Takagi's product-state result depends on the initial vibrational state and collision energy. Its n=2-dominated limit is a low-vibration, low-energy limit; higher electronic states occur outside it. The code itself says that its H2+ population is vibrationally hot and outside that condition.
- Guberman's calculation is for ground-vibrational `3HeH+`, with electron energies from 0.001 to 0.33 eV. The paper identifies a dominant n=2 route and isotope sensitivity; it does not establish a universal unit n=2 branching ratio for arbitrary HeH+ populations and temperatures.
- The actual energy available includes initial internal excitation. Schematically,

  ```text
  translational energy released
      = Q_ground + E_internal,reactants - E_internal,products,
  ```

  with incident electron energy and rate weighting handled consistently in the thermal energy equation. Increasing the product principal quantum number alone does not prove that `Q_ground-E(n=2)` is an upper bound when the initial vibrational energy also changes.

This is a defect in the claimed validity bound, not a measurement of the sign or magnitude of the resulting wind error. In a reduced model that has already deposited nascent ion excitation into the gas, simply adding the excitation again at recombination could itself count energy twice. The full reaction sequence must be audited before changing a coefficient.

**Recommendation.** Treat the current prescription explicitly as a low-state recipient approximation, not a universal correction. Record a validity flag when its assumed reactant-state distribution is unavailable. Longer term, use population- and temperature-weighted branching probabilities and a consistent internal-energy reservoir. For the immediate reduced model, propagate a stated recipient uncertainty through a small set of representative cells and complete reaction cycles. When excited H is disabled, explicitly identify the assumption that the removed excitation energy escapes rather than being collisionally recovered or reabsorbed.

### 5.3 R6: a distribution peak is not its mean energy

`molecular_reaction_heat.f90:554–576` uses the midpoint of H2 v=5 and v=6 energies, approximately 2.4795 eV, as the internal energy of every two-body recombination product. The code acknowledges that this is a peak-based estimate.

Kokoouline's quoted distribution comes from its direct-pathway calculation, not a complete calculation of all indirect pathways. Strasser measures a broad distribution and discusses rotational excitation of the initial ions and products. Neither a peak location nor uncertainty in that location establishes the mean of the full rovibrational energy distribution.

Therefore the stated approximately 15% uncertainty is not a demonstrated uncertainty in the mean energy needed by the heating ledger. The relevant quantity is `sum_(v,J) P(v,J)*E(v,J)`, with a normalized, applicable product distribution.

**Recommendation.** Keep the value only as a labeled model estimate until an applicable distribution is obtained. If a figure must be digitized, include digitization and population assumptions explicitly. Test separate uncertainty in branching and internal energy. The fact that a scalar quench factor approaches one can reduce sensitivity within that model; it does not validate the chosen distribution.

### 5.4 The scalar quenching expression is not a proved cascade bound

`h2_vibrational_relaxation.f90:565–633` computes

```text
C1 = n_H*k10_H + n_H2*k10_H2 + n_He*k10_He,
f  = C1/(C1 + A_max).
```

Using the largest total radiative decay rate is a meaningful improvement over using one transition's radiative rate. It does **not** establish a bound for the entire excited ladder unless the adopted collision rate also bounds the relevant energy-removal rates throughout that ladder. The comment that closely spaced high levels necessarily have larger effective de-excitation rates is not a proof. Selection rules, collider dependence, rotational redistribution, destructive reactions, and upward thermal transitions all matter.

Lique's section 3.2 does report a modest increase with initial vibrational state for a fixed vibrational change in the H-collision transitions it studied. That is useful supporting evidence within that data set. It does not establish the required minimum energy-removal rate for every higher level, rotational state, and collider reached by nascent molecules.

For example, if one important level had `C_u < C1`, even `A_u <= A_max` would not force its collisional branching fraction to exceed the adopted scalar fraction. Moreover, a first-event branching fraction is not automatically the total energy fraction recovered through a multistep cascade.

The handoff's quoted `1-f = 7e-9` measures the implemented scalar model. It is not a bound on its difference from an uncomputed level-population model. The comments claiming an all-level lower bound and negligible discrepancy should be corrected.

**A practical way past the missing-data problem:** build a reduced statistical-equilibrium model only for levels and colliders supported by available data; keep an explicitly uncertain high-energy population group rather than inventing precise missing rates. Enforce detailed balance for upward/downward thermal rates. In a collision-only thermal-bath test with no pumping or radiative escape, the populations should approach the intended Boltzmann distribution and net collisional heating should vanish. Then inject a known formation/pumping energy with the physical loss channels restored and test its partition into heat, escaping radiation, chemical destruction, and stored internal energy. Missing high-level data become quantified model uncertainty instead of a claim that no useful model can be built.

The comment near line 613 that photoelectric/Lyman-Werner callers omit helium is also stale: the current call sites in `ionization_equilibrium.f90` and `util_ion_eq.f90` pass a helium density. That part of L7h is implemented.

### 5.5 Argon as helium: an approximation, not a certified upper bound

`mol_rates.f90:540–621` correctly distinguishes collider weights and uses the same effective collider density in the R12/R15 pair. Cohen and Westberg's published recommendation supports the H2 and H coefficients quoted there. Preserving the forward/reverse ratio through a common coefficient is a useful consistency property.

However, `k3b_H_H_to_H2_monatomic` evaluates an argon rate and applies it to helium while claiming `k(He) <= k(Ar)` from mass and polarizability. Those qualitative properties do not prove an ordering of thermally averaged quantum recombination rates over the implemented temperature range.

Paolini et al. explicitly analyze He and Ar with different interaction potentials, resonance/continuum contributions, and approximations. Their discussion includes potential-energy-surface sensitivity, the distinction between equilibrium and steady-state intermediate populations, and an additional exchange contribution important for Ar at low temperature. This supports treating the substitution cautiously, not assigning a rigorous global upper bound from collider mass alone.

**Recommendation.** Label the rate as an Ar-based estimate for He, with its temperature range and uncertainty; do not present it as measured helium data or a proved bound. Obtain applicable helium measurements/theory where possible, including controlled extraction of published numerical curves when needed. Vary the uncertain helium efficiency in bounded cell calculations before choosing new catalog reference results. Keep the same choice in both R12 and R15 so detailed balance is maintained.

An additional validity error appears in the justification for extrapolating the pair: disappearance of H2 above the front does not make **both** reactions vanish. R15 forms H2 from atomic H and scales with `n(H)^2` times the collider sum. Its validity cannot be dismissed using the preexisting H2 abundance alone. Check the actual association source at out-of-range temperatures, not only the molecular abundance. A common treatment of both directions can preserve detailed balance; it is not necessary to extrapolate indefinitely just to preserve their ratio.

### 5.6 Consequences for catalog decisions

Old molecular roots do not remain roots after changing their equations; the handoff is right to flag that. But replacing reference outputs after a regression change proves consistency with the new implementation, not the physical correctness of every new approximation.

Before a full new molecular campaign, require:

- Energy accounting for complete reaction sequences, including radiation and stored excitation, without double counting.
- Collider-sum agreement between chemistry, carrier balance, and heat evaluation.
- Detailed balance and thermal-limit tests.
- Explicit validity/uncertainty records for the three approximations above.
- A small representative set spanning molecular depth, transition front, and dilute upper atmosphere.

These are acceptance criteria even if a particular observable changes only slightly. Unknown effect size is not a reason to call an unsupported bound established physics.

## 6. A concrete route through the coupled-solver blockage

### 6.1 What the existing evidence does and does not establish

The handoff reports 13 exhausted 40-vector cycles, a relative linear residual of 0.2597 against a requested 0.1, and certification by the alternating route at pass 12. These are **quoted historical run results**, not new solves performed here.

They justify investigating the linear system. They do not uniquely identify a defective species preconditioner. Mean directional-Jacobian discrepancies of `3.692e-4` and `3.792e-3` do not bound the largest error in the stiff rows that stop convergence. A finite-difference comparison can also share a residual defect with the action it is checking.

Existing capability must be used accurately: `EXHALE_GM_M` is read at `steady_newton.f90:16610–16613`; `EXHALE_GM_CYCLES` and `EXHALE_GM_CYCLES_FROM` are read around lines 7951–7966. `pgmres_with_restarts` exists at lines 13376–13455 and is invoked by the trust-region step near line 15485. A restart implementation is not missing. The plain `pgmres` call sites elsewhere do not automatically inherit that wrapper's behavior.

### 6.2 First experiment: one frozen rejected Newton system

At a saved coupled iterate, hold the unknowns, composition seed, column/row scales, active bounds, pseudo-time shift, and reconstruction linearization fixed. Do not let trial chemistry or caches silently change the operator between samples.

Measure:

1. Repeated residual and Jv evaluations, including interleaving an unrelated state and returning to the checkpoint.
2. A perturbation ladder for Jv, with maximum and RMS errors reported separately for mass, momentum, energy, elemental, and carrier rows. Locate the worst cells. Respect positivity and shared abundance bounds when choosing perturbations.
3. Directional additivity/homogeneity and consistency of the Arnoldi image with an independently evaluated action on the completed direction.
4. The independently recomputed linear residual `b-A*delta`, not just the recurrence residual. `A` must include the same scaling and pseudo-time shift as the solve.
5. Basis orthogonality, preconditioned residual history, and whether changing the chemistry stopping tolerance changes the attainable linear residual floor.

If Jv accuracy limits the linear solve, widening the Krylov subspace will not cure the underlying problem. If the measured action is stable and a stronger preconditioner sharply improves convergence, the preconditioning hypothesis gains direct support.

### 6.3 Bounded candidate comparisons

Compare a 40-vector cycle, two 40-vector cycles, and an 80-vector cycle on the same checkpoint and operator. Record total residual evaluations, achieved **true** linear residual, memory, and accepted nonlinear progress. Equal iteration counts are not equal cost, and timings across differently loaded hosts are not a clean comparison.

For the stronger model, use the block structure

```text
J = [ J_hh  J_hs ]
    [ J_sh  J_ss ],
S = J_ss - J_sh * inverse(J_hh) * J_hs.
```

Here `h` denotes hydrodynamic unknowns and `s` transported composition unknowns. A useful approximation to the species Schur complement should retain the measured hydro–composition feedback and full transport stencil, not just local chemical diagonal terms. Conversely, eliminating stiff local chemistry can reduce the retained system if the local solve and its derivative are controlled. Radiation attenuation and the nonlocal base closure can introduce couplings outside a narrow spatial band; measure those contributions rather than assuming the two-cell reconstruction entries explain the entire difficulty.

A small-grid explicitly assembled Jacobian and direct solve provide a diagnostic reference. They are not a proposal to use a dense direct solve for the full catalog. A base-only coupled block should be considered only if the measured problematic modes are actually confined to that region; the documented outer-column stalls make a universal base-only cure unlikely.

### 6.4 Slow outer iteration and the unused handover

The handoff correctly states that the automatic handover was not exercised. It should not be counted as a validated rescue route.

Its progress estimate also needs correction: `0.0816 -> 0.0479` over 11 passes is a geometric factor of approximately 0.953, or **4.7% reduction each pass**, not 10%. Continuing that unverified constant factor from 0.0479 to `1e-5` would require roughly 175 additional passes, not 80. Neither extrapolation predicts the actual future path.

Use residual histories and movement of the worst cell to distinguish slow transport relaxation from a fixed stalled mode. A moving front may call for multilevel acceleration, controlled continuation, or acceleration of the outer fixed-point map. Any acceleration must be accepted using the full physical residual and admissibility checks. Removing movement limits alone is not a demonstrated solution.

## 7. Boundary physics, flux choice, and grid evidence

### 7.1 What remains valid about L26

The repaired window weighting is present in `base_boundary.f90:542–614` and 747–769. The variance is accumulated in a second pass about the mean, avoiding the cancellation-prone subtraction of two large moments. The amplitude factor tends to zero quadratically as the window mean flux vanishes. These are sensible improvements to residual smoothness.

However, smoothness is not physical locality. The boundary still depends instantaneously on a remote wind window. That is a stated nonlocal closure, not a derivation of a local characteristic boundary. Equality on five sampled states does not establish that the closure changes no physics on arbitrary transients or candidate states.

The negative local **cell-centered** flux readings still warn that naive extrapolation is unreliable at the base. What must be withdrawn is the claimed production-face-flux evidence used to exclude all local closures. A local construction using the same well-balanced reconstruction and a characteristic boundary/Riemann compatibility solve has not been disproved by PLM evaluation of WENO3 roots.

Recommended bounded study: keep the imposed reservoir pressure/thermodynamic data, formulate compatibility with outgoing characteristics for the intended subsonic regime, use the actual production face flux, and test hydrostatic equilibrium, weak inflow, reversal, and perturbations that change only the remote window. Do not select a branch merely to reproduce an already certified profile. If the remote closure is retained, document its range and test sensitivity to window location and domain extension.

### 7.2 L25 is useful evidence, not a completed convergence study

Roe success and continuation success support the user's limited seven-case choice. They do not establish Roe accuracy for the whole catalog. Similarly, HLLC success at 1000 cells demonstrates that the 500-cell difficulty can be relieved by a changed grid; it does not by itself uniquely identify the mechanism or prove grid convergence.

As the handoff itself notes, the grid change also moved the endpoints. Hold physical face boundaries fixed and compare HLLC and Roe at matched resolutions, with the same source and boundary models. Compare actual conserved face fluxes, integrated energy balance, first-cell thermodynamics, and observational outputs. A later third resolution is needed to assess a convergence trend.

The optional low-Mach velocity-jump modification should remain a controlled experiment. Failure on this stratified near-static base does not invalidate the published method in its intended regime, and success of another flux does not validate the base thermodynamics automatically.

## 8. Parallel execution and optimization

The L15 stage-2 record is stronger evidence than the earlier unlocalized thread-dependence report: it compares full iterative traces on the stated fixture. Those results were read, not reproduced in this review. They do not establish race freedom on all molecular or coupled paths.

The inspected carrier residual uses a maximum reduction and a subsequent deterministic selection of the worst row (`diffusive_photochemistry.f90:4567` onward). This is a reasonable design. Chemistry reached from these parallel loops depends on thread-local module state; future derivative/coloring parallelization must audit every called routine's mutable state, not merely make the top-level arrays private.

The first optimization target should be **fewer expensive residual evaluations for validated progress**, not more threads. The handoff's quoted coupled-iteration cost is already dominated by model construction and Krylov products. Reusing a preconditioner can help only while a measured quality criterion remains satisfied; stale radiative or chemical derivatives can erase the saving. Avoid simultaneous unrestricted OpenMP and BLAS threading, and report both policies in benchmarks. No new race or scaling result is claimed here.

## 9. Recommended decisions and order of work

| Priority | Action | Completion criterion |
|---|---|---|
| 1 | Repair operator identity in diagnostics and solve scoping in log tests | Production and diagnostic fluxes agree on the same installed state; repeated-solve logs cannot cross-match records |
| 2 | Reclassify archived certification claims using the current evaluation route in isolated copies | Each claim is labeled reproduced or refused under a named binary/operator; old headers are not presented as current acceptance |
| 3 | Unify spherical transport geometry and expose consistent elemental face fluxes | Exact-volume conservation identity, correct boundary budget, and derivative checks pass |
| 4 | Replace unsupported molecular bounds with explicit models and validity/uncertainty records | Reaction-cycle energy accounting, detailed balance, and thermal-limit tests pass |
| 5 | Diagnose one frozen coupled linear system and compare bounded preconditioner/Krylov candidates | True residual reduction and useful nonlinear progress at measured evaluation cost |
| 6 | Validate one new molecular reference solution, then expand | Full residual and physical budget checks plus uncertainty assessment; no tolerance relaxation |
| 7 | Complete matched-domain flux/grid comparisons and stage transport | Resolution trends and stage-sum conservation are demonstrated |

The current evidence does **not** justify opening a project whose premise is “hundreds of times the wind flux pass the current mass gate.” It does justify adding an independent face-flux budget to certification diagnostics, provided it is evaluated on the same state and operator. It also justifies retaining the limited Roe choice pending grid evidence and investigating species preconditioning after checking Jv quality.

## 10. Handoff statements that should be revised

- Remove the obsolete “Until that block is filled” binary paragraph now that the final binary block is filled.
- Replace “No existing default changed” with a scoped statement. Molecular default equations and the boundary shape threshold did change; optional solver switches remaining off is a different claim.
- Replace the base-flux certification accusation using sections 3.1–3.2 of this review.
- Correct the interpretation of the log-reader failure using section 3.3.
- Describe the molecular recipient and collider prescriptions as established data **plus explicit approximations**, not all as proved corrections with rigorous bounds.
- Correct the first well-mixed run's average progress rate and identify all extrapolations as conditional arithmetic.
- Limit “all 95 outputs” to files that actually exist; the same handoff says seven cases have no output state.
- Separate historical molecular states failing changed equations from inability to obtain any new molecular root. The reported I3 alternating success already shows that these are different claims.
- Do not treat refreshed reference outputs or bit-identical selected tests as proof of physical validity outside the tested regime.

## 11. Retained evidence and reproducibility

New review artifacts in `docs/audit_20260905/`:

- `handoff_20260917_flux_operator_probe.f90`: direct production flux comparison.
- `handoff_20260917_flux_fiducial.log`.
- `handoff_20260917_flux_x010_HeH2.13.log`.
- `handoff_20260917_flux_x003_HeH2.13.log`.
- `handoff_20260917_flux_x002_HeH2.13.log`.
- `handoff_20260917_log_scope_probe.py` and `handoff_20260917_log_scope.log`.
- `handoff_20260917_x003_evaluate.log`: complete production evaluation of the refused historical claim.

The Fortran diagnostic was compiled with `gfortran -O0 -g -fbacktrace -fopenmp`, module includes from `build/`, and all `build/*.o` except `EXHALE_main.o`, `*_tests.o`, and `*_probe.o`. It linked the prefix's OpenBLAS and `-ldl`. Runs used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Copy the stated input and restart pair into a scratch directory, put the pair at `output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`, and resolve spectrum paths without altering the spectrum data.

The parser check is reproducible with:

```bash
python3 docs/audit_20260905/handoff_20260917_log_scope_probe.py \
  LHS1140b/models/.L22/i3_alt/run.log
```

The production evaluation used the old 0.03 pair with only the scratch input's spectrum path resolved and `Restart intent: stationary evaluate`. It is expected to exit nonzero because certification is refused; that is the result being checked, not a numerical crash.

Published source files used for the molecular analysis are available locally:

- [Takagi (2002)](../../references/Takagi_2002_PhST_96_52.pdf): initial-state-dependent dissociative recombination and product excitation.
- [Giusti-Suzor et al. (1983)](../../references/Giusti-Suzor_1983PRA_28_682.pdf): low-state H2+ recombination limits.
- [Guberman (1994)](../../references/Guberman_1994PRA_49_R4277.pdf): ground-vibrational HeH+ calculation and dominant excited-H route.
- [Kokoouline et al. (2001)](../../references/Kokoouline_2001Nature_412_891.pdf) and [Strasser et al. (2001)](../../references/Strasser_2001PhysRevLett_86_779.pdf): H3+ fragmentation and product internal excitation.
- [Cohen and Westberg (1983)](../../references/Cohen_1983JPCRD_12_531.pdf): collider-specific recommended association rates, including the H/H2 sheet on journal page 559.
- [Paolini et al. (2011)](../../references/Paolini_2011PRA_83_042713.pdf): quantum treatment of He/Ar recombination and its physical approximations.
- [Lique (2015)](../../references/Lique_2015_MNRAS_453_810.pdf): H–H2 rovibrational collision calculations and data scope.

No production fix is claimed in this review. The deliverables are this analysis and the retained diagnostic evidence. The confirmed diagnostic defects, geometric inconsistency, unresolved model-validity claims, and untested solver proposals are deliberately kept distinct.

# Review of `PLAN_20260920_rev2.md`

Review date: 2026-09-20 (KST).

## Overall verdict

Revision 2 is substantially improved and is close to an executable investigation plan. It correctly adopts the Mode R/Mode C distinction, corrects the sweep-count interpretation, removes `On stall` from the experiment list until repaired, adds output and exit gates, and separates residual-route, closure, Jacobian, conservation, M1, and M2 questions.

I recommend accepting the structure with a small set of mandatory corrections before execution:

1. Correct the carrier row-term channel count: the current writer outputs 15 H2 channels, not 13.
2. Define how the independent hydrodynamic conservation audit obtains momentum and energy face fluxes and source terms. The current plan says "exported faces," but the existing output path visibly exports face mass flux and carrier terms, not a complete hydrodynamic face/source ledger.
3. Replace a few causal statements about the alternating M1 route with measured, conditional wording.
4. Record the effective parsed configuration, not only environment variables and input hashes; defaults and derived flags materially control the route.
5. Pin the residual-route equivalence experiment so that the three paths evaluate the same state without silently changing composition, caches, reconstruction mode, or boundary state.

No production code or case state was changed for this review. The findings below come from direct inspection of the current source and the stored records cited by rev2.

## 1. Changes that are correct

### 1.1 Mode R and Mode C

The two-mode distinction is correct and useful.

Mode R is a valid experiment when its claim is limited to applying the current executable and current boundary reconstruction to an old stored state. `load_IC` explicitly reports a boundary-model mismatch, rebuilds the boundary from the current physical column and reservoir, and warns that the residual is not the residual under which the restart was written. This is not a historical-certificate reproduction and must not be used as evidence that the current branch converges.

Mode C is required for current-branch convergence, continuation, or solver-performance claims. The plan correctly assigns the coupled trial and Stage F to Mode C.

One addition is required: "current operator" must include the effective parsed configuration, not merely the binary hash, boundary identity, input hash, and environment. Several defaults are derived after parsing. For example, carrier transport can be enabled by the oxygen-chemistry configuration when the input does not state it, and the coupled route is refused for incompatible ionization-transport combinations. Record a machine-readable effective-options line or the complete parser summary beside every result.

### 1.2 Sweep-count correction

Rev2 correctly corrects rev1. `n_eq_sweeps_model` starts at one but is assigned from `n_eq_sweeps_last` and can be increased with `max` during the run. A Jacobian evaluation uses a fixed count for that evaluation, but the count can change between outer iterations. The proposed two-dimensional experiment in probe length and fixed closure count is therefore appropriate.

The experiment should additionally record whether the fixed count was actually honored when an evaluation exits early for an inadmissible composition, a missing chemical root, or a discarded state. "Requested count" and "completed count" are not enough if an evaluation returns a refusal flag; record the refusal reason and whether `F` is admissible for derivative use.

### 1.3 `On stall` diagnosis

The reachability defect remains correctly identified. The present ordering resets `carrier_outcome`, skips the relaxation when `outer_ending` is already `outer_no_progress`, and then requires the reset outcome to be `carrier_relax_movement_bound` before returning. The conjunction cannot hold in the current pass.

The proposed explicit transition and immutable entry snapshot are the correct design direction. The snapshot must include the adopted composition, fixed-wind state, face mass flux, row vectors/scales, movement-bound reason, progress history, and boundary/reconstruction state. The five proposed tests are also appropriate.

Do not implement this as a simple "last outcome" scalar. The coupled block must recompute and verify a fresh residual from the exact entry snapshot before accepting the handover. Otherwise a stale outcome token can authorize a state that no longer satisfies the condition that triggered the transition.

### 1.4 Diagnostic gates and residual separation

The requirement that every diagnostic produce an expected file, valid header, row count, state identity, and exit class is correct. "No output" and "timeout" must never be interpreted as a numerical pass.

The three-way separation

```text
r_model = -F - J_action(s) - shift(s)
r_trial = F(Y+s)
r_state = F(Y_adopted)
```

is also correct. The first tests the supplied shifted linear model, the second tests nonlinear prediction, and the third is the state-level result.

## 2. Confirmed remaining issue: H2 channel count

Section 9.3 states that `EXHALE_CARRIER_ROW_TERMS=1` writes "the thirteen channels" for H2. The current source defines:

```text
n_h2chan_mol      = 13
n_h2chan_oxy_prod = 14
n_h2chan_oxy_loss = 15
n_h2chan          = 15
```

and `carrier_row_terms_write` emits `(rowdump_h2ch(j,k), k = 1, n_h2chan)`. Therefore the output has 15 channel values when the oxygen chemistry channels are present; the first 13 are the molecular channels and the last two are oxygen-cycle production and loss channels. The header also uses `n_h2chan` when constructing the channel names.

This is not merely editorial: a parser or row-count check based on 13 will misread the output columns. Replace "thirteen channels" with "up to 15 channels: 13 molecular channels plus two oxygen-chemistry channels when enabled," and require the header to determine the active count.

## 3. Conservation audit is not yet operationally specified

The plan correctly asks for an independent finite-volume conservation audit and correctly warns against summing the same residual implementation twice. However, section 13 still assumes that all required hydrodynamic faces and sources are already exported. Current source inspection shows a more limited situation:

- `steady_residual.f90` assembles `R_mass`, `R_momentum`, and `R_energy` from `dF`, `S`, heating/cooling, and optional viscous/conductive sources.
- `Num_Fluxes.f90` forms the hydrodynamic flux, including energy flux `v*(E+p)`.
- The loaded-state residual diagnostic writes a residual profile with lower and upper **mass** face fluxes.
- The carrier row-term output writes carrier faces, production/loss, and reaction channels.

The plan does not identify an existing file that exports hydrodynamic momentum and energy face fluxes, geometric source terms, viscous/conductive source terms, and cell volumes in a form suitable for an independent reconstruction. Consequently, "reconstruct from exported faces and sources" is currently a design requirement, not an immediately runnable instruction.

Before scheduling this phase, choose one of two explicit designs:

1. Add a diagnostic-only export at the hydrodynamic assembly boundary, containing cell volume, both face fluxes for all three conserved variables, geometric/gravity source, viscous source, conductive source, radiative heating, radiative cooling, and the exact signs/units; or
2. Build an independent postprocessor from an immutable checkpoint and a separately documented re-evaluation, explicitly acknowledging that this reuses production primitive/flux routines and is not an independent implementation of the flux arithmetic.

The first option is more auditable. In either case, test the identity on subdomains and report signed defect, sum of absolute local defects, and largest local defect. Also state whether the diagnostic includes ghost cells, spherical face areas, and the well-balanced pressure/gravity decomposition. Without these details the conservation phase can produce a plausible but dimensionally incorrect number.

## 4. Residual-route equivalence requires stricter controls

Section 5.1 is valuable, but the three routes do not automatically evaluate the same map:

- `EXHALE_RESIDUAL` is a loaded-state diagnostic path with its own reconstruction and composition-refresh policy.
- Newton `eval_residual` uses a fixed sweep count for its Jacobian map and maintains module workspaces/caches.
- The marching/RK path advances a state and uses the marching reconstruction and source sequence.

The plan says to compare them at one checkpoint, but it must specify a non-mutating protocol. For each route, record before and after:

- conserved state and composition hashes;
- ghost state and base boundary cache;
- particle counts, temperature, heating, and cooling;
- face fluxes and reconstruction/limiter mode;
- composition sweep count and seed hash;
- all relevant module-cache or branch identifiers.

If a route changes the state while preparing derived quantities, save the pre-route state and restore it before the next route. A route equivalence test that runs the three paths sequentially without restoration can mistake state mutation for operator disagreement. If the marching route cannot be evaluated without advancing time, use a single RHS evaluation at a frozen state and label it as an RHS comparison rather than a stationary residual comparison.

The result must be a term-by-term map, not only three norms. Otherwise the proposed first measurement can itself create the ambiguity it is intended to resolve.

## 5. M1 wording and interpretation

The implementation evidence supports the statement that the default route performs an alternating hydrodynamic/composition iteration and that carrier/element rows are actively relaxed. One sentence in section 9.1 remains too strong:

> "the refusing row is the row the relaxation drives."

The relaxation targets the carrier or element balance, but the hydrodynamic update changes density, velocity, pressure, temperature, rates, face fluxes, and therefore the next composition residual. The refusing row is the residual of the coupled map, not necessarily the direct residual of the most recent relaxation. Replace that sentence with: "the refusing carrier or elemental row is one component of the block iteration; its change reflects both the fixed-wind relaxation and the subsequent hydrodynamic refresh."

The planned full-state norms and two-cycle test are correct. Add a convergence metric for the **map defect**: evaluate the next alternating-map state from the same input and report `||G(Y_k)-Y_k||` separately from `||Y_{k+1}-Y_k||`. A small state displacement can result from a movement bound while the unconstrained map defect remains large.

For the coupled trial, "same current checkpoint" must include the same composition refresh stage. If the `False` control is started from an already refreshed state but the `True` route initializes its carrier unknowns from a different pre-refresh state, the comparison is not matched.

## 6. Stage B/C remaining requirements

### 6.1 Projection and limiter state

The active-set mask and projection displacement are correctly added. Also record the bound treatment separately for each endpoint. If one endpoint is projected and the other is only step-shortened, the central difference is not symmetric and should not be labeled a central derivative.

Record whether WENO smoothness weights, contact-mode dissipation, positivity repairs, and well-balanced data were frozen or recomputed. The plan mentions reconstruction and limiter state, but these specific flags can change the operator and should be part of the effective configuration record.

### 6.2 Scaling provenance

The requirement to record the state used to form `Drow` and active boxes is correct. Add the exact arrays' hashes or a deterministic summary (min/max/norm and worst-cell indices), because a textual state identity alone does not prove that the scaling arrays were not overwritten by a diagnostic evaluation.

### 6.3 Closure count and derivative validity

The two-dimensional `h x k` design is sound. Add one rule: do not use a closure count for which the inner map has unresolved chemical roots or a branch change as a Jacobian validation point. Such a point may be scientifically interesting, but its finite difference is not a derivative of the intended smooth map.

## 7. M2 physical audit

The M2 section is materially improved. Including the base ghost and face state is important because pressure-preserving interior conversion does not guarantee preservation of the boundary particle count or face state.

One additional requirement is needed: distinguish "charge conservation" from the explicitly implemented conversion invariants. The conversion source comments and checks emphasize H nuclei, He nuclei, mass, nonnegative fractions, and EOS closure. If charge is checked, document the species and ionization bookkeeping used to derive it and whether electrons are included in the conserved charge. Do not report charge conservation as confirmed merely because nucleus and mass budgets pass.

For the energy decomposition, include the sign convention and units of every term, and report whether transport is active. The plan correctly notes that viscous and conductive sources enter the momentum and energy residual; these terms cannot be omitted from a cell-level explanation of an energy failure even if their integrated contribution is small.

## 8. H3+ diagnostics

The endpoint-counting defect and cumulative-counter limitation are correctly retained. The new output gate is appropriate.

The proposed failing-cell analysis should also record whether the H3+ function was evaluated more than once for that cell during the selected residual route. Because the counters are call counts, a cell-level table should be generated from a controlled final-state evaluation, not inferred by dividing the global counters by the number of cells. Report the actual `T`, `n(H2)`, and `n(H3+)` branch for that evaluation.

## 9. Evidence and execution discipline

Section 12 is strong. Two small additions are recommended:

- Record compiler and linked-library identities for any diagnostic executable, not only the EXHALE binary and source manifest. This matters for OpenMP, BLAS, LAPACK, and floating-point behavior.
- Record the effective CPU/thread configuration and floating-point mode where available. `OMP_NUM_THREADS` and BLAS thread count are necessary but may not be sufficient to reproduce reduction order.

The claim that inherited H3+ evidence was not rebuilt is now clearly labeled. The earlier successful audit log should be treated as inherited, while any new source inspection should remain labeled `INSPECTED` rather than `MEASURED`.

## 10. Recommended execution order

The rev2 order is appropriate with the following small adjustment:

1. A0 identity and effective parsed configuration.
2. Non-mutating residual-route equivalence on one checkpoint, or a documented reason why the marching route cannot be included.
3. Full-vector repeatability with fixed compiler/runtime/thread metadata.
4. Closure count versus probe length.
5. Carrier and hydrodynamic conservation budgets, after the hydrodynamic export design is fixed.
6. M1 alternation map-defect/cycle audit.
7. M2 invariant, charge bookkeeping, and base-ghost audit.
8. One admissible M1 coupled trial in Mode C, with a matched `False` control.
9. P2 Stage B/C/D diagnostics.
10. Stage F continuation only after separate authorization.

This order avoids spending campaign time on a spectrum or continuation before establishing that the compared residual routes, closures, and conservation identities are coherent.

## 11. Final decision

`PLAN_20260920_rev2.md` should be accepted after these edits:

- correct the H2 channel count from 13 to "up to 15" and make parsers header-driven;
- specify and implement the actual hydrodynamic face/source export needed for the independent conservation audit;
- pin and restore state/caches for residual-route equivalence;
- weaken the M1 causal wording and add an unconstrained map-defect metric;
- include effective parsed options and compiler/runtime metadata in every artifact;
- add charge-bookkeeping and transport-source details to the M2 audit.

After those changes, rev2 provides a physically and numerically defensible sequence for distinguishing stale state identity, closure error, operator error, active constraints, boundary effects, and genuine current-branch behavior. It still cannot prove nonexistence of a stationary atmosphere from a failed solve or a finite continuation ladder.

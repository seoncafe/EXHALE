# Recommendations for D2a and D8

Date: 2026-09-18

Status: Recommended decisions for user approval. This document does not record approval or implementation.

## 1. Executive recommendation

1. **D2a: Adopt B for the advective boundary, with separate, explicit boundary conditions for species diffusion and thermal conduction.** The reservoir supplies the entropy and composition of incoming material. The interior supplies the thermodynamic trace during reverse flow and at a pressure-balanced stationary contact. Pressure and velocity must satisfy the chosen acoustic boundary relation.
2. **D8: Adopt immutable generation directories, a manifest, and one atomic publication index.** Readers select one generation and read all state components from it. Define restart as a warm restart that preserves the physical state while resetting solver controls according to a documented policy.

These recommendations concern physical and storage contracts. Neither a favorable convergence result nor agreement with an existing output is sufficient to establish physical correctness.

## 2. Scope and evidence

This recommendation follows inspection of the current D2a/D8 discussion in [PLAN_20260918_rev2.md](PLAN_20260918_rev2.md) and the following implementation paths:

- [base_boundary.f90](../src/modules/states/base_boundary.f90): hydrostatic continuation, entropy-branch weights, density mixing, and the outgoing acoustic velocity relation, particularly lines 650–789 in the inspected version.
- [Representative molecular base.inp](../LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13/base.inp): reservoir pressure, temperature, molecular abundance, elemental abundance, and eddy diffusivity. Its description identifies the molecular composition as a condition on incoming gas.
- [binary_element_diffusion.f90](../src/modules/functions/binary_element_diffusion.f90): construction and use of the reservoir helium mass fraction in the diffusion solve, particularly lines 768–790.
- [viscous_conduction.f90](../src/modules/time_step/viscous_conduction.f90): the thermal-conduction coefficients use the base ghost temperature, particularly lines 455–480.
- [load_IC.f90](../src/modules/files_IO/load_IC.f90): selection and loading of the hydro and species components, including lines 424–425.
- [write_output.f90](../src/modules/files_IO/write_output.f90): hydro output fields and restart metadata, including lines 272–311.

The published [Thompson (1987) paper](../../references/Thompson_1987JCP_68_1.pdf) was also consulted for the characteristic interpretation. Its one-dimensional Euler eigenvalues include the acoustic speeds and the entropy/contact speed in Equation (60). The relevant example sets dissipative terms to zero at the boundaries; it does not establish a boundary condition for the diffusion and conduction operators used here.

This was a code and document inspection. No new numerical experiments, bitwise comparisons, restart round-trip tests, or NFS interruption tests were performed for this recommendation. Source line numbers identify the inspected version and can move with subsequent edits.

## 3. D2a: Adopt B for advection

### 3.1 Physical interpretation

The input quantities `p_base`, `T_base`, `q_H2_base`, and elemental abundances define an exterior reservoir. They do not, by themselves, require the entropy of stationary interior gas or gas returning to the reservoir to equal the reservoir entropy.

For a subsonic lower boundary, define positive flow as flow from the reservoir into the computational domain:

| Local boundary state | Entropy and composition carried by advection | Acoustic condition |
| --- | --- | --- |
| Inflow into the domain | Reservoir values | Prescribed incoming condition matched to the outgoing interior acoustic information |
| Reverse flow into the reservoir | Interior values | The same consistent acoustic matching, with the appropriate characteristic count |
| Pressure-balanced stationary contact | Interior trace retained | Preserve zero velocity and pressure balance |

For a pressure-prescribed, stationary boundary, the interior trace is conceptually

\[
\rho_b = \rho_{\mathrm{EOS}}(p_{\mathrm{res}},s_i,Y_i),
\]

where `s_i` and `Y_i` are the interior entropy and composition continued to the boundary radius. The continuation must use the molecular EOS and a consistent treatment of gravity. This expression describes the pressure-prescribed variant; an incoming-acoustic prescription that produces a different boundary pressure must state that pressure explicitly.

A stationary contact can have equal pressure and velocity but different density and entropy on its two sides. The exterior reservoir and the interior trace therefore need not collapse into a single averaged thermodynamic state. Which trace a numerical flux uses must follow the boundary formulation.

B is not the only model permitted by characteristic theory at zero contact speed. It is recommended because it avoids introducing an additional thermal-contact assumption into an advective boundary. Supersonic cases require a separate characteristic count and must not inherit the subsonic pressure prescription without justification.

### 3.2 Why changing the midpoint alone is insufficient

The inspected implementation first constructs a branch weight from the interior and wind-window information, then evaluates

```fortran
rho_rev = isentropic_density_at_pressure(1, rho_i, p_i, T_i, p_b)
rho_b   = (1.0d0 - w_rev)*rho_res + w_rev*rho_rev
v_b     = v_i + (p_b - p_i)/(rho_i*c_i)
```

At zero discriminants, the existing smooth transition produces an intermediate density. This is an interpolation rule, not an independently specified law of thermal contact.

Moreover, an initially stationary interior state need not produce a stationary boundary: a pressure mismatch can produce nonzero `v_b`. The entropy/composition selection must be consistent with the resulting local boundary flow, not merely with the velocity before acoustic matching.

Accordingly:

1. Define the pressure/acoustic boundary problem and its thermodynamic traces.
2. Determine the local mass-flux or contact-wave direction consistently with that problem.
3. Select the advected entropy and composition using that direction.
4. Apply exactly the same physical boundary operator in marching and stationary residual evaluation.

The remote wind window may remain useful for diagnostics or a justified numerical treatment. It must not silently override a local reverse flow and inject reservoir entropy or composition without a physical argument.

Simply assigning a new value at the single point `M = 0` is not a complete implementation. A contact with different entropy on its two sides has direction-dependent traces. If differentiability is required, distinguish numerical regularization from physical exchange, consider regularizing the flux consistently, and test convergence as the regularization width decreases. A smooth density blend is not automatically a physical heat-transfer law.

### 3.3 Diffusion and conduction remain separate

The phrase “the reservoir supplies composition only during inflow” must refer to **advection**, not all transport.

The current element-diffusion implementation constructs `X_base` from the reservoir He/H ratio and passes it to the diffusion solve. The thermal-conduction operator uses a base ghost temperature. Both processes can transport material or energy when the bulk velocity is zero.

The boundary contract must therefore distinguish:

- **Advection:** choose reservoir or interior entropy and composition according to the consistent flow direction.
- **Species diffusion:** explicitly choose reservoir composition, specified flux, zero flux, or a justified finite exchange law. Zero bulk velocity does not imply zero diffusive flux.
- **Thermal conduction:** when enabled, explicitly prescribe temperature, heat flux, or finite thermal contact. Changing an advective ghost state must not accidentally change this thermal boundary condition.

The same principle applies to molecular and eddy transport. Each enabled operator needs a stated boundary condition and a compatible conservation budget.

### 3.4 When A or C would be appropriate

**A** is defensible if the modeled lower boundary is explicitly a surface whose thermodynamic state is maintained by a reservoir, including at rest. That is an additional physical assumption. Its implementation must still respect the characteristic count during reverse flow and avoid overconstraining the interior.

**C** is appropriate if finite thermal contact is part of the intended physics. It requires an energy-exchange law and a justified conductance, length, or timescale, with the corresponding energy flux included in the budget. The present half-density interpolation has neither such a law nor such a scale and should not be presented as an established thermal-contact model.

Introducing an arbitrary relaxation time only to improve convergence is not recommended as a physical closure.

### 3.5 Acceptance tests and sequencing

The focused acceptance tests should cover:

1. Hydrostatic equilibrium compatible with the molecular EOS and the discrete gravity treatment.
2. A pressure-balanced stationary contact with unequal reservoir and interior entropy, verifying that the chosen trace does not introduce spurious flow.
3. Small pressure perturbations that launch both inflow and reverse flow, verifying acoustic signs and the entropy/composition source.
4. Zero bulk velocity with enabled diffusion or conduction, checking species and energy flux budgets rather than assuming zero transport.
5. Consistency between marching and stationary residual evaluation, and convergence with any numerical transition width.

Do not claim that all flowing solutions remain unchanged. Algebraic equality can hold where both closures select the same saturated branch, but weak nonzero flow can differ, and a changed boundary operator can change the global steady solution. Bitwise identity is a measured property of specified cases, not a physical acceptance criterion.

D5b may first establish reproducibility of the current model, provided it is labeled as such. The subsequent D2a change should have its own boundary-model identity and validation results. Reproducing the old interpolation does not validate its physics.

## 4. D8: Adopt a generation-based publication contract

### 4.1 Recommended structure

The following names are proposed, not existing interface guarantees:

```text
case/
  runs/<run_id>/
    ... execution inputs, logs, and isolated evaluation work ...
  states/<generation_id>/
    Hydro_ioniz.txt
    Ion_species.txt
    manifest.json
  state_index.json
```

The existing numeric data-file formats can initially be retained if they satisfy the physical-state round-trip contract below. A new container format is not required solely to make publication atomic.

Published generation directories are immutable. A generation identifier must remain unique across retries and concurrent attempts; an iteration number alone is insufficient. Subsequent evaluations should produce separate artifacts referring to the immutable state identity rather than modifying a published state.

### 4.2 Distinguish completion from certification

Use one atomic index containing, at minimum:

- `latest_complete`: the latest fully stored generation that passes the declared storage and state-consistency checks.
- `latest_certified`: the latest generation that also passes the stated physical and numerical certification criteria.

Here, completion means completion of the snapshot, not convergence of the solve. A valid, fully stored checkpoint can remain uncertified. A failed or unfinished continuation must not replace the previous certified reference with an uncertified claim.

Keep the solve ending reason, admissibility, evaluability, and certification status explicit. Do not introduce a `best` reference until its ranking rule is defined. If added later, update all references in the same index transaction.

### 4.3 Manifest and reader contract

The manifest should record:

- Schema version, run/attempt/generation identifiers, parent state identity, and iteration/snapshot phase.
- Source and executable identities and the effective configuration identity.
- Grid, units, species schema, EOS/reaction-network identity, and relevant boundary-model identity.
- Component filenames, dimensions, and hashes.
- Storage completion, admissibility, solve ending reason, and certification status.
- The exact state identity assessed by any attached certificate, including the assessment configuration.

Every reader must resolve the index **once**, select one generation, and read all components from that immutable directory. It must not resolve a moving reference separately for each component or substitute a missing component from another generation.

An incomplete or inconsistent generation is rejected. Selection of an older complete checkpoint is allowed only as an explicit, reported policy. Readers may see an older complete generation under filesystem caching; that is different from constructing a mixed state.

### 4.4 Warm restart contract

A warm restart preserves the defined physical state and reconstructs derived quantities while initializing solver controls according to a documented policy. It does not promise the same iteration trajectory as an uninterrupted solve.

Specify:

- Which physical variables are authoritative and which are reconstructed.
- The treatment of independent species, algebraically constrained species, and molecular thermodynamics.
- The accuracy with which density, momentum, energy, and species inventories survive a write/read cycle.
- How `dtau`, trust bounds, relaxation parameters, counters, and progress history are initialized.
- Which configuration changes are accepted as a restart, and which constitute a new initialization or mapped seed.

The current hydro output records primitive thermodynamic quantities. With a molecular EOS, verify that the stored primitives and species reconstruct the original conserved energy to the declared tolerance. If they do not, store the conserved energy explicitly. Warm restart is not permission to alter the physical state silently.

When the executable, EOS, reaction network, grid, or boundary model changes, an old certificate remains historical provenance. Certification under the new operator requires a new evaluation. Exact algorithmic continuation is a separate future capability and is not part of the recommended initial contract.

### 4.5 Publication, concurrency, and NFS

Publish in this order:

1. Write all components into a new, unpublished generation directory.
2. Close the files and check dimensions, finite values, component consistency, and hashes.
3. Complete the manifest and the required data/directory synchronization.
4. Write a replacement index on the same filesystem and atomically replace the published index.

A same-filesystem `rename` supports atomic replacement of the publication name. Atomic visibility is not the same as crash durability: synchronizing a file does not by itself synchronize its containing directory. See the primary Linux interface documentation for [rename](https://man7.org/linux/man-pages/man2/rename.2.html) and [fsync](https://man7.org/linux/man-pages/man2/fsync.2.html).

For this NFS workspace, initially establish and test a bounded guarantee: interruption of the publishing process leaves readers with either the previous complete generation or the new complete generation. Do not claim protection against server failure without validating the relevant storage behavior.

An atomic index replacement does not prevent competing publishers from losing each other's updates. Initially enforce one publisher for each case, with ownership or locking that is tested on the actual filesystem. Keep old generations; automatic deletion requires a separate retention policy.

### 4.6 Integration and acceptance tests

D8a can proceed independently with preservation, explicit continuation conditions, and isolated evaluation. The generation interface should then be introduced together with its readers, including `load_IC`, the mapper, seed selection, status, and reproduction tools.

Preserve existing products. A legacy file pair can be checked and imported as one explicitly identified state without recomputing its numerical data. Legacy import must not manufacture a certificate or claim provenance that is unavailable.

Test interruption before any data are written, during each component write, between component writes, during manifest preparation, and around index publication. Also test:

- A reader overlapping publication never mixes generations.
- An uncertified new generation does not overwrite the certified reference.
- A certificate always identifies the exact evaluated state.
- Competing publishers are rejected or correctly serialized.
- A warm restart reconstructs the physical state within its declared tolerance.
- Evaluation cannot modify solve inputs or published state files.

The focused checkpoint writer should accept one synchronized snapshot and must not change subsequent solver behavior. Reusing the full output path is acceptable only if that property is established; the existence of an output routine is not sufficient evidence.

## 5. Proposed approval text

### D2a

Adopt B for the subsonic advective lower boundary. The reservoir provides the entropy and composition of incoming material. Reverse flow and a pressure-balanced stationary contact retain the interior thermodynamic trace. Determine the flow direction consistently with acoustic matching. Specify independent physical boundary conditions for species diffusion and thermal conduction. Treat numerical smoothing as a declared regularization, not an unnamed thermal-contact model.

### D8

Adopt immutable generation directories, manifests, and one atomic publication index. Readers select one generation and use it for all state components. Distinguish storage completion from physical certification. Define restart as a warm restart with validated physical-state reconstruction and a documented solver-control initialization policy. Enforce one publisher for each case and test interruption behavior on the actual NFS filesystem. Preserve existing products and keep evaluation separate from published states.

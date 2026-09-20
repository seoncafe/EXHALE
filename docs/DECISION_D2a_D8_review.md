# Review of `DECISION_D2a_D8_codex.md`: the D2a and D8 recommendations against the source and the record

Written 2026-09-18 (KST) at the user's request. This is the advisor's
review of the Codex recommendation document; it records what of that
document was verified in the source, where the reviewer agrees, and where
the record of this code argues for a different choice. It records no
approval and no implementation. Every number is READ from the file named
beside it.

## 1. Verdict in two lines

- **D8: adopt as recommended.** The generation directories, the manifest,
  one atomic index with `latest_complete` and `latest_certified` kept
  apart, readers resolving the index once, a warm-restart contract, one
  publisher for each case and a bounded NFS guarantee are the right shape;
  they are what `PLAN_20260918_rev2.md` D8 steps 3 to 7 asked to have
  agreed before implementation, made concrete. Two additions are noted in
  section 3.
- **D2a: the recommendation of B is defensible as characteristic theory
  but conflicts with the stated meaning of the reservoir in this code's
  record; the reviewer recommends A at the pressure-balanced stationary
  contact, with everything else of the Codex text kept** (B's treatment
  of inflow and reverse flow, the separate boundary conditions for species
  diffusion and thermal conduction, the flow direction decided after the
  acoustic matching, the regularization named as numerical). The choice
  between A and B at the stationary contact is the user's; section 2 sets
  out the two readings and their consequences.

## 2. D2a

### 2.1 What the Codex document says, verified

| statement | where checked | verdict |
|---|---|---|
| `rho_rev = isentropic_density_at_pressure(1, rho_i, p_i, T_i, p_b)`, `rho_b = (1 - w_rev) rho_res + w_rev rho_rev`, `v_b = v_i + (p_b - p_i)/(rho_i c_i)` (the C- relation), so at zero discriminants the blend gives an intermediate density that is an interpolation and not a contact law | `base_boundary.f90` lines 650 to 789 (the reservoir isentrope at line 670, `p_res = base_reservoir_nhat rho_res T_res` at 671, the `C-` pair at 685) | verified |
| a stationary interior does not give a stationary face: a pressure mismatch launches `v_b` through `C-` | the same relation; L35's rest fixture (base mass flux 1.2751e-5 with the interior at rest) | verified |
| `base.inp` describes the molecular composition as a condition on the INFLOWING gas | `LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13/base.inp` header: "The inflowing gas keeps the hydrogen-nucleus fraction bound into H2 that the photochemical column carries at 1 microbar" | verified |
| the element-diffusion operator takes its base composition from the reservoir He/H as a Dirichlet value, independent of the flow direction | `binary_element_diffusion.f90` lines 778 to 788 (`X_base = m_He HeH/(m_1 + m_He HeH)`, "base reservoir composition (Dirichlet)") | verified |
| the thermal-conduction operator's base face uses the anchored ghost temperature | `viscous_conduction.f90` line 470 ("left face (base face at j = 1 uses the anchored ghost)") | verified |
| Thompson (1987) gives the eigenvalues and does not fix a boundary condition for diffusion or conduction | the local PDF, as the first review already read it | accepted |

### 2.2 Where the reviewer differs: what `T_base` means in this code

The Codex argument for B is that `p_base`, `T_base`, `q_H2_base` and the
element reservoirs "define an exterior reservoir" and "do not, by
themselves, require the entropy of stationary interior gas ... to equal
the reservoir entropy". As a statement about characteristic theory that is
right: at a pressure-balanced contact at rest the two sides may carry two
entropies. The question is what THIS model's lower boundary is meant to
be, and the record answers it:

- The base level is the hand-over level of the lower-atmosphere model (1
  microbar; Koskinen et al. 2022 for the hot Uranus, the Photochem column
  for the LHS 1140 b `photochem` cases), and `T_base` is the temperature
  that model gives AT that level. It is a temperature of the level, not
  only of the gas that happens to be moving upward through it. The lower
  atmosphere below it is dense and radiatively controlled, and on the
  time scales of a stationary solution it is a heat bath: a column at rest
  above it would take the bath's temperature at the level by conduction and
  radiation, which is the physics behind writing `T_base` as a prescribed
  temperature at all.
- Item L21 (2026-09-15/16, `docs/lhs1140b_stationary_L21_20260915.md`)
  found that "the 226 K base was entering the state with a weight of
  0.026, and the first cell of a '226 K base' run sits at 417.6 K", and
  treated it as a defect to be corrected: the interior isentrope standing
  at the base instead of the stated reservoir was wrong. Every catalog
  state since then was re-solved on that correction. B at the stationary
  contact reintroduces exactly that reading for the one state that does
  not flow: the face takes the interior isentrope at `p_res`, and
  `T_base` is not enforced when nothing flows.
- Where the two choices differ in practice: only at rest and in the
  weak-flow neighborhood of rest. On every flowing catalog state the branch
  is inflow, the reservoir owns the contact under A and under B alike, and
  L35 measured the two closures bit-identical there. The states that
  differ are static tests, hydrostatic starts and the base of a column
  whose wind has stalled, which are exactly the states where a
  well-balanced scheme is judged, and there A gives the reservoir
  temperature the input states and B gives the interior's.

The reviewer therefore recommends **A at the pressure-balanced stationary
contact**: the reservoir maintains the thermodynamic state of the level
including at rest, which is the additional physical assumption the Codex
text names in its section 3.4, and it is the assumption this model has
made since it wrote `T_base` into `base.inp` and since L21. Everything
else of the Codex D2a text stands: during reverse flow the interior trace
leaves with the outgoing characteristics (A must "respect the
characteristic count during reverse flow and avoid overconstraining the
interior", the Codex caveat, which is exactly the L26 reversal branch);
the flow direction is decided after the acoustic matching and not from the
interior velocity; a numerical transition width, if one remains, is
declared as regularization and its convergence tested; the wind window
does not override a local reverse flow.

If the user prefers B, the reading of `T_base` changes to "the temperature
of incoming material only", the L21 correction is re-read as "the
reservoir owns the base while the base is inflow", and `docs/input_schema.md`
and `MODELS.md` section 4 say so. That is a legitimate model, but it is a
change of the boundary's stated meaning, not a clarification of it.

### 2.3 The point the Codex document adds that the plan lacked, adopted

Section 3.3 of the Codex text is right and was missing from
`PLAN_20260918_rev2.md` D2: the advective boundary is not the whole
boundary. The element-diffusion operator imposes the reservoir composition
as a Dirichlet value whatever the flow does, the thermal conduction
operator reads the anchored ghost temperature, and the molecular carrier
transport has its own ghost rule (`X_ghost = X_N`, no diffusive flux across
the base, L22 step 2c). Each enabled transport operator therefore needs
its own stated boundary condition (prescribed value, prescribed flux, zero
flux, or a named exchange law) and its own place in the conservation
budget, and a change of the advective ghost must not silently change
them. D2b's specification and tests are extended with the Codex section
3.5 list: hydrostatic equilibrium with the molecular equation of state and
the discrete gravity, a stationary contact with unequal entropies and no
spurious flow, acoustic perturbations of both signs, zero bulk velocity
with diffusion and conduction on (species and energy flux budgets, not an
assumption of zero transport), marching against stationary evaluation,
and convergence in any regularization width. "All flowing solutions
unchanged" is not claimed; it is measured on named states.

### 2.4 Sequencing, agreed

D5b first makes the CURRENT boundary reproducible and labels it as the
current model (the plan's D5b already says so); the D2a change then
carries its own boundary-model identity in the state's provenance and its
own validation. Reproducing the old interpolation validates its
reproducibility, not its physics.

## 3. D8

### 3.1 Adopted as recommended

The structure (`states/<generation_id>/` immutable with the two data files
and a manifest; `runs/<run_id>/` for inputs, logs and isolated evaluation;
one `state_index.json` with `latest_complete` and `latest_certified`), the
manifest content (schema version, identifiers, parent state, iteration
phase, source and executable and configuration identities, grid, units,
species schema, equation-of-state and network and boundary-model identity,
component hashes and dimensions, storage completion, admissibility, ending
reason, certification status, the exact state identity a certificate
assessed), the reader rule (resolve the index once, read every component
from one immutable generation, refuse an incomplete one, choose an older
complete one only as an explicit reported policy), the publication order
(write, close and validate, manifest and synchronization, atomic index
replacement on the same filesystem), one publisher for each case, no
automatic deletion, and the warm-restart contract (authoritative against
reconstructed variables, the treatment of constrained species and the
molecular thermodynamics, the round-trip accuracy of density, momentum,
energy and species, the initialization of `dtau`, trust bounds, relaxation
parameters and counters, which configuration changes are a restart and
which a new initialization) are adopted as the agreed shape of D8 steps 3
to 7. The existing numeric file formats are kept; a container format is
not required for atomic publication.

Two points the Codex text makes deserve to be written into D8's
acceptance explicitly:
- **Completion is not certification**: a fully stored checkpoint can be
  uncertified, and a failed or unfinished continuation must never replace
  `latest_certified`. This is the concrete form of the L34c finding
  (the continuation's worse state overwrote the first solve's better one).
- **The conserved energy round trip**: the hydro file stores primitives;
  with the molecular equation of state, whether `rho`, `v`, `p` and the
  species reconstruct the original `u(3)` to the declared tolerance is a
  measurement D5a is making now; if they do not, the conserved energy is
  stored explicitly. D8's writer follows D5a's answer.

### 3.2 Two additions

1. **Legacy import**: every existing catalog `output/` pair is imported as
   one explicitly identified legacy generation with the provenance it
   actually has (the binary its `REPRODUCE.md` names, the L29 verdict, the
   L34 re-evaluation), without a manufactured certificate; `status.py`
   reads the index and nothing else. This is how the 95 catalog cases
   enter the new contract without recomputation.
2. **The D8a increment now in execution** (preservation without a new
   format, the continuation policy, the evaluate pass in its own
   directory, the campaign's aggregate status) is the first step and stays
   compatible with the structure above: `output_solve1/` and the `ENDING`
   file it introduces become legacy generations when the index arrives.

### 3.3 What is not agreed by this review

Exact algorithmic continuation is not part of the initial contract (the
Codex text says so too); the `best` reference is not introduced until its
ranking rule is defined; crash durability across a server failure is not
claimed until the storage behavior is tested on this NFS.

## 4. Proposed approval text, as the reviewer would write it

**D2a.** Adopt the Codex D2a text with one change: at the pressure-balanced
stationary contact the reservoir owns the thermodynamic state of the
level (A), the interior trace leaving only with the outgoing
characteristics of a reverse flow; the flow direction is decided after
the acoustic matching; species diffusion and thermal conduction carry
their own stated boundary conditions and budget entries; any transition
width is a declared numerical regularization with a convergence test.
(Or, if the user prefers B: adopt the Codex text as written and re-state
`T_base` in `docs/input_schema.md` and `MODELS.md` section 4 as the
temperature of incoming material only.)

**D8.** Adopt the Codex D8 text as written, with legacy import of the
existing products as identified generations and with D8a as its first
increment.

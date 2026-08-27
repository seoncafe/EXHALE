# Staged refactor plan: §5.2 System dedup, §5.3 composition, §5.6 parser

These three are the large, high-risk items from the 2026-07-13 code review
(`docs/code_review_20260713_en.md` §5.2, §5.3, §5.6). Each is a quality
refactor, not a bug fix, and touches validated, paper-critical code. They are
**not** one-shot changes: each is broken into increments that are individually
small, and **every increment is gated by the byte-identical regression**
(`make check`, then the wasp_full + wasp_he23off cases; seven cases as of
2026-08-27) plus, where a molecular or
loaded-SED path is affected, a molecular reference run (`examples/15_molecular`,
bounded) and a loaded-SED run. Stop and re-validate physics if any gate needs a
golden re-snapshot.

Guiding rule: **operation order matters.** Replacing an explicit expression with
a loop or a helper can reorder floating-point adds and break byte-identity (as
the §4.1 assumed-shape attempt showed). Preserve the exact arithmetic order, or
treat the increment as a deliberate golden re-snapshot and flag it.

---

## §5.2 — Deduplicate the `System_*` ionization residual modules

**Current state.** `System_HeH`, `System_HeH_TR`, `System_HeH_metals`,
`System_HeH_TR_metals`, `System_HeH_mol`, `System_H`, plus the advection
(`System_implicit_adv_*`) and post-process variants each repeat H/He reaction
logic and communicate through a position-indexed `params(N)` array. A single
displaced index is not type-checkable and is a known fragility.

Two facts the increments below must respect:

- **Analytic Jacobians duplicate the residuals.** `System_H`, `System_HeH`, and
  `System_HeH_metals` carry hand-written `jac_system_*` routines (used by the
  Newton solver, `newton_solver.f90`) — the same reaction physics in derivative
  form. The TR, TR_metals, and mol systems have no analytic Jacobian (hybrd1
  finite differences only). Any change to a residual body must keep its paired
  Jacobian consistent.
- **OpenMP threadprivate state.** The ionization sweep over cells runs
  OpenMP-parallel; `System_HeH_metals` holds threadprivate module state
  (`met_*`), and the charge-transfer state (`cx_*`) in `charge_exchange` is
  threadprivate too. Refactors that move this state between modules or into
  derived types can silently break the parallel sweep.

**Increments (each byte-identical unless noted):**

- **Inc 0 — named params indices.** Replace the magic `params(20)`, `params(21)`,
  … with named `integer, parameter` constants. **One flat constant set is not
  possible:** the slot meanings differ between families — in the equilibrium
  systems `params(1)` is P_HI (and `params(15)` is q13 in the TR variants),
  while in `System_implicit_adv_*` `params(1)` is dr/v and `params(15)` is the
  local He/H ratio. Define one named set for each family (e.g. `IPE_*` equilibrium,
  `IPA_*` advection) in a shared module, and apply it at BOTH ends: the packing
  sites (`ionization_equilibrium.f90`, `post_process_adv.f90`) and the
  `System_*` unpacking. Pure renaming → byte-identical. Removes the
  displaced-index fragility and is a prerequisite for the rest.
  Gate: WASP + molecular byte-identical.
  These named constants (`params_idx`) were the Inc 0-3 stepping stone and were retired by Inc 4's named-field cell state (`ion_cell_state`), which replaced all `params` packing/unpacking; the module was removed once nothing referenced it.
- **Inc 1 — H/He common core.** The H/He rows are identical within PAIRS, not
  across all four non-mol systems (verified against the bodies): the standard
  three rows (with collisional ionization) are verbatim-shared by `System_HeH`
  and `System_HeH_metals`; the TR-form four rows (no collisional ionization,
  Oklopcic He form) are verbatim-shared by `System_HeH_TR` and
  `System_HeH_TR_metals`, and the triplet row is additionally verbatim in
  `System_HeH_mol` (its row 8). Inc 1 extracts the standard three rows (plus
  the matching Jacobian local terms shared by `jac_system_HeH` and
  `jac_system_HeH_metals`) into small pure subroutines that write into
  `fvec` in the SAME order as today; the TR-form rows move in Inc 3. Convert ONE `System_*` at a time to call
  them; regress after each. The advection systems (`System_implicit_adv_*`)
  are fraction-form with time terms and share no verbatim block — they stay
  out of the extraction. Where the converted system has a paired analytic
  Jacobian (`jac_system_H/HeH/HeH_metals`), update the Jacobian in the SAME
  increment — either extract a matching derivative helper or re-check it term
  by term. A residual/Jacobian drift changes Newton's iteration path, and the
  hybrd1 fallback can rescue convergence on the gate cases while the drift
  bites elsewhere — do not rely on the byte gate alone to catch it. Gate each
  conversion.
- **Inc 2 — metals block.** Extract the metal-loop residual contribution shared
  by `System_HeH_metals` and `System_HeH_TR_metals` into one pure subroutine.
  The two copies are verbatim except the metal row base (`ix = 4 + 2*(e-1)` vs
  `5 + 2*(e-1)`) — pass the base as an argument (mirroring `cx_metal_base`).
  The safest shape for the threadprivate state (`met_*`, `cx_*`): keep it
  declared where it is and have the CALLERS pass their thread's arrays as
  actual arguments, so the helper stays pure and stateless; run the OpenMP
  identity check (see the checklist). Gate.
- **Inc 3 — TR and molecular blocks.** Extract the TR-form four rows
  (verbatim-shared by `System_HeH_TR` and `System_HeH_TR_metals`) the same way,
  with the He 2^3S triplet row as its own helper so `System_HeH_mol` (whose
  row 8 is verbatim that row; its molecular rows 1-7 are unique and stay put)
  calls it too. No analytic Jacobian exists for the TR/mol systems, so there
  is no Jacobian counterpart to keep in step. Gate (molecular reference).
- **Inc 4 — derived types (bigger).** Introduce `type(ion_cell_state)` and
  `type(ion_rates)` with named fields to replace the `params` packing in
  `ionization_equilibrium.f90` and all residuals. This touches packing/unpacking
  everywhere; do it after Inc 1-3 so the residual bodies already read through
  helpers. Keep each thread's cell state local (subroutine variables) or
  threadprivate — a shared module variable of the new types breaks the parallel
  sweep silently. Gate byte-identical (including the OpenMP identity check);
  watch operation order.
- **Inc 5 — metadata-driven assembly (optional).** Assemble residual/Jacobian by
  iterating active-species metadata (`species_table`). Highest reordering risk;
  may require a golden re-snapshot — decide explicitly.

  **Inc 5 decision (2026-07-16): NOT implemented — goal already met.** Scoping
  after Inc 0-4 showed the generic assembler would add nothing and cost real
  structure:
  - The only block that GROWS with new species — the metal rows — is already
    metadata-driven (`met_*` coefficient arrays in canonical `species_table`
    order, `met_top` staging, `metal_rows`/`metal_fractions`/
    `metal_electron_sum` helpers, `cx_add_to_fvec` generic charge exchange).
    Adding a metal element is already "add table rows".
  - The remaining cores are deliberately DIFFERENT physics, not duplicated
    code: the standard H/He rows carry collisional ionization, the TR form
    omits it (Oklopcic triplet formulation) and adds the 2^3S kinetics, and
    the molecular system has its own reaction network (rows 1-7). A generic
    species-iterating assembler would need a special case for each of these,
    reintroducing in table form the branching it is meant to remove.
  - The standard family's hand-written analytic Jacobians (Newton solver)
    have no generic counterpart; assembly-driven residuals would fall back to
    fdjac1 finite differences (slower) or require a matching generic Jacobian
    assembler (new, unvalidated code on paper-critical physics).
  - Any assembly reordering forces another golden re-snapshot, for zero
    physics payoff.
  If a future need arises (e.g. a species whose rows do not fit the current
  three shapes), revisit with that concrete case in hand.

**Recommended stop point:** Inc 0-3 give most of the dedup value at low risk.
Inc 4-5 are optional and carry byte-identity risk. (Executed 2026-07-16:
Inc 0-4 done; Inc 5 explicitly declined — see above.)

---

## §5.3 — Centralize composition calculations

**Current state.** `calc_rho`, `calc_ntot`, `calc_ne` (`utilities.f90`),
initialization, boundaries, and chemistry each branch separately on He / metals /
molecules. The §3.4 monochromatic-He inconsistency was a symptom of this scattered
state logic.

**Increments:**

- **Inc 0 — extend `species_table`.** Add H, He, and molecular species to the
  metadata (mass, charge, electron count for each), which is currently
  metals-focused. Pure data addition → byte-identical.
- **Inc 1 — metadata-driven EOS helpers.** Rewrite `calc_rho`/`calc_ntot`/
  `calc_ne` to sum mass / particle / electron counts by iterating the species
  metadata instead of hardcoded He/metal/mol branches. **Byte-identity risk:** a
  metadata loop can reorder the FP adds versus the current explicit expressions.
  Either preserve the exact term order, or accept a deliberate golden
  re-snapshot — flag this as a decision point before doing it.
- **Inc 2 — single source of truth.** Route `mass_per_H`, `rho_bc`, `ntot_bc`,
  and the BC/IC composition through the same centralized helpers so composition
  cannot disagree between code paths (the root cause behind §3.4). Gate.

  **Inc 2 executed 2026-07-16 (bit-identical, no re-snapshot).** The base
  scalars flow through `comp_mass_per_H`/`comp_ntot_bc`/`comp_rho_bc`
  (composition.f90); the mono-He path sets `HeH = 0` before that block, so it
  uses the same source. Two deliberate exceptions stay: (1) `load_IC.f90`
  keeps the historical H/He-only mass formula so legacy IC reloads remain
  bit-identical — documented in-file, do not "fix"; (2) `calc_mmw`
  (utilities.f90, a post-process diagnostic) still hardcodes the He mass 4.0
  and could be moved onto `bsp_mass` in a follow-up. The wind-ae converters
  intentionally use H/He-only conversions matching their metal-free oracle.

**Note:** §5.3 Inc 1 is the increment most likely to require a re-snapshot. Treat
it as a decision, not an automatic byte-identical change.

**Inc 1 executed 2026-07-16 (re-snapshot accepted).** Measured record: the
reorder seed is 1-3 ulp (max rel diff 7.3e-16 across all columns at step 1);
the converged wasp states differ by median ~1e-14 / max ~2e-10 with identical
step counts and unchanged Mdot (13.30 both cases); OMP 1-vs-16 stays mutually
byte-identical. Caveat for future gates: the bounded 400-step molecular
checkpoint sits in a chaotic transient — the stiff H2 front amplifies any
1-ulp seed by roughly e per step, so the front decorrelates by 1-2 cells at
that step count; this is transient decorrelation, not a physics change.
Goldens and the bounded molecular/SED/OMP references were re-snapshotted.

---

## §5.6 — Unify input parsing

**Current state.** The Fortran parser mixes positional mandatory header lines
(read by line order) with keyword-scan optional keys later; the Python utilities
(`examples/exhale_io.py`, and `exhale_transit_lib.py`'s `get_word` used by
`EXHALE_transit.py`) parse the same file independently. Highest back-compatibility risk of the three, because a
parser change affects every run and every shipped `input.inp`.

**Increments:**

- **Inc 0 — document the schema.** Write one authoritative schema table (every
  key, positional-vs-keyword, type, default, legacy alias). No code change;
  prerequisite for everything else.
- **Inc 1 — Fortran positional → keyword, with back-compat.** Convert the
  positional header reads to keyword `index(line,'KEY')` matching (extension keys
  already work this way), while STILL accepting the legacy positional layout
  (dual-mode or a one-time migration of all shipped inputs). **Hard constraint:**
  every existing `input.inp` (positional) must parse identically. Gate: parse
  each `examples/01`-`15` and every planet-folder input at startup and confirm
  identical derived parameters, plus byte-identical output where feasible. This
  is the risky increment — an exhaustive input-corpus regression is mandatory.
- **Inc 2 — shared Python schema.** Make the Python loaders read the same
  key-value schema, removing duplicate parsing logic.
- **Inc 3 — centralized validation.** One validation pass with explicit legacy
  aliases and clear errors.

**Recommended:** Inc 0 (document) is free and valuable on its own. Inc 1 only
with the full input-corpus regression; back-compat with positional inputs is the
gating constraint.

[2026-08-27 status, read off the code: **Inc 0 done** -- the schema table is
`docs/input_schema.md`. **Inc 1 done** -- the core block of
`src/modules/files_IO/input_read.f90` is matched by anchored label
(`lbl_match`, 73 sites), so physical line order no longer matters while
legacy positional files parse unchanged because their lines are
self-labeling; the input-corpus regression the increment made mandatory is
`EXHALE_PARSE_DUMP` / `write_parse_dump` in `EXHALE_main.f90` driven by
`backup/regression/run_parse_corpus.sh`, which snapshots `parse_dump.txt`
for every `input.inp` in the tree. **Inc 2 not done** --
`examples/exhale_io.py:read_input` still parses `input.inp` on its own, and
`exhale_transit_lib.py` still has its own `get_word`. Inc 3 was not checked
here.]

---

## Cross-cutting gate checklist (run for every increment)

0. Before the FIRST increment: confirm `make check` already passes at the
   starting commit (goldens current), so any later failure is attributable to
   an increment and not to stale goldens.
1. `make check` — wasp_full + wasp_he23off byte-identical.
2. Molecular reference (`examples/15_molecular`, bounded) byte-identical if the
   molecular path is touched.
3. Loaded-SED run byte-identical if `sed_read`/energy-grid is touched.
4. `run_fcheck` (bounds/FPE) periodically to catch out-of-bounds from index
   changes.
5. OpenMP identity when `System_HeH_metals`, `charge_exchange`, or the sweep
   scratch is touched: one bounded run with `OMP_NUM_THREADS=1` vs. 16 must
   stay byte-identical. The regression harness itself is single-thread and
   does NOT cover this.
6. If an increment intentionally changes the last bits (operation-order or
   precision), it is a **golden re-snapshot decision**, not a silent change —
   surface it, and remember re-snapshotting does not restore reproducibility
   against already-published numbers.

---

## 2026-07-16 physics-correctness audit (pre-Inc 4)

Criterion: physical correctness only — impact size and backward compatibility
are never arguments (see the user-level instructions). Two systematic sweeps
(mass/density accounting; electron/particle/pressure accounting) found that
the eos_metals policy had been wired everywhere but the LATER molecular
network had not. Fixed (all byte-identical for atomic runs; molecular
reference re-snapshotted for the corrected physics):

- calc_mmw: metal mass + nuclei under the eos policy (the _adv temperature
  solve was using an H/He-only mmw while its ne already carried the metal
  electrons). _adv T shifts ~0.3% median on wasp_full.
- load_IC: rho reconstruction now uses calc_rho (the run's own mass policy:
  HeITR + metals + molecules); the old H/He-only formula left a ~1.1% mass
  discontinuity on a metals-on restart (now closed to 4e-16). Molecular f_sp
  columns are now restored; f_sp fully zero-initialized.
- set_IC: f_sp zeroed before BOTH IC branches (the cold-hydrostatic branch
  left the molecular columns undefined).
- energy_semi_implicit: molecular densities passed to calc_ne/calc_ntot —
  neutral H2 was missing from n_tot (first-order at a molecular base; the
  bounded molecular checkpoint changes mainly at the H2 front, r~1.03).
- EXHALE_transit.py: the free-electron density is now metal-aware (stage-
  weighted metal columns; +molecular ions when present) — ne was underestimated
  up to ~30x at the near-neutral base where low-IP metals dominate, biasing
  the n=2 Balmer rates. Metals-off files keep the legacy expression.
- Documented-approximation notes added where molecular ions are deliberately
  neglected as trace electron donors (eval_cool, excited_hydrogen, T_equation,
  dp_bc), and in lower_column (metal mass absent from its mu). post_process_adv
  now states its composition scope (H/He + metals; no molecular treatment).
- Output label corrected: Hydro_ioniz column 2 is rho in m_H/cm^3 (metals
  included), not a number density; header now says rho[mH/cm3].

Still-correct-as-is (judged, not deferred): wind-ae converters (their oracle
is genuinely metal-free); charge_exchange (depends on neutral H, not ne);
brem Z^2 weighting; the System_* internal n_e sums.

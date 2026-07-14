# Staged refactor plan: §5.2 System dedup, §5.3 composition, §5.6 parser

These three are the large, high-risk items from the 2026-07-13 code review
(`docs/code_review_20260713_en.md` §5.2, §5.3, §5.6). Each is a quality
refactor, not a bug fix, and touches validated, paper-critical code. They are
**not** one-shot changes: each is broken into increments that are individually
small, and **every increment is gated by the byte-identical regression**
(`make check`, the wasp_full + wasp_he23off cases) plus, where a molecular or
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

**Increments (each byte-identical unless noted):**

- **Inc 0 — named params indices.** Replace the magic `params(20)`, `params(21)`,
  … with named `integer, parameter` constants (e.g. `IP_T`, `IP_NTOT`, …) in a
  shared module, used across all `System_*`. Pure renaming → byte-identical.
  Removes the displaced-index fragility and is a prerequisite for the rest.
  Gate: WASP + molecular byte-identical.
- **Inc 1 — H/He common core.** Extract the H and He ionization/recombination
  residual contributions that are IDENTICAL across the four non-mol systems into
  small pure subroutines (`add_H_residual`, `add_He_residual`) that write into
  `fvec` in the SAME order as today. Convert ONE `System_*` at a time to call
  them; regress after each. Gate each conversion.
- **Inc 2 — metals block.** Extract the metal-loop residual contribution shared
  by `System_HeH_metals` and `System_HeH_TR_metals` into one pure subroutine.
  Gate.
- **Inc 3 — TR and molecular blocks.** Extract the He 2^3S triplet block and the
  molecular block (`System_HeH_mol`) the same way. Gate (molecular reference).
- **Inc 4 — derived types (bigger).** Introduce `type(ion_cell_state)` and
  `type(ion_rates)` with named fields to replace the `params` packing in
  `ionization_equilibrium.f90` and all residuals. This touches packing/unpacking
  everywhere; do it after Inc 1-3 so the residual bodies already read through
  helpers. Gate byte-identical; watch operation order.
- **Inc 5 — metadata-driven assembly (optional).** Assemble residual/Jacobian by
  iterating active-species metadata (`species_table`). Highest reordering risk;
  may require a golden re-snapshot — decide explicitly.

**Recommended stop point:** Inc 0-3 give most of the dedup value at low risk.
Inc 4-5 are optional and carry byte-identity risk.

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

**Note:** §5.3 Inc 1 is the increment most likely to require a re-snapshot. Treat
it as a decision, not an automatic byte-identical change.

---

## §5.6 — Unify input parsing

**Current state.** The Fortran parser mixes positional mandatory header lines
(read by line order) with keyword-scan optional keys later; the Python utilities
(`examples/exhale_io.py`, `EXHALE_transit.py`'s `get_word` reads) parse the same
file independently. Highest back-compatibility risk of the three, because a
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

---

## Cross-cutting gate checklist (run for every increment)

1. `make check` — wasp_full + wasp_he23off byte-identical.
2. Molecular reference (`examples/15_molecular`, bounded) byte-identical if the
   molecular path is touched.
3. Loaded-SED run byte-identical if `sed_read`/energy-grid is touched.
4. `run_fcheck` (bounds/FPE) periodically to catch out-of-bounds from index
   changes.
5. If an increment intentionally changes the last bits (operation-order or
   precision), it is a **golden re-snapshot decision**, not a silent change —
   surface it, and remember re-snapshotting does not restore reproducibility
   against already-published numbers.

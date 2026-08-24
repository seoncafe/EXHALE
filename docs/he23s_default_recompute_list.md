# He 2³S default switched on — what needs recomputing

**Written 2026-08-22. Nothing in this list has been recomputed; this file is
the work list, not a record of work done.**

## What changed

The metastable helium triplet is now **on by default**:

- `parameters.f90`: `thereis_HeITR = .true.` (was `.false.`).
- `input_read.f90`: `Include He23S` became an **optional** key (it was
  mandatory through `req`, so the old initializer never actually applied to
  a valid run). A file that omits the line now gets the triplet; `Include
  He23S? False` is the explicit opt-out.
- Tk interface (`src/utils/EXHALE_interface_functions.py`): the checkbutton
  starts checked, the reset button restores it instead of clearing it, and
  the file reader treats the key as optional with the Fortran semantics.
  Since a file may now omit the line, the sequential reader checks the label
  before consuming the line, so an older file without it no longer shifts
  `Load IC` / `Do only PP` / `Force start` by one.
- Verified on a one-step run: key absent -> the setup report prints
  "Including helium triplet chemistry"; `Include He23S? False` -> it does
  not.  The Tk reader parses the tutorial input with and without the line to
  the same values.

Because the triplet feeds the electron budget, the energy budget (its
Penning, recombination and 10830 Å channels) and the observable itself, a
case that flips from off to on changes physically. Stored results from
before the flip are stale.

## Deliberately off — do NOT change, do NOT recompute

These exist precisely to exercise the triplet-off branch. They set the key
explicitly, so the default change does not touch them and the goldens stay
valid.

| case | why it stays off |
|---|---|
| `backup/regression/wasp_he23off/` | the HeITR-off branch check in the default regression matrix |
| `backup/regression/wasp_he23off_newton/` | the same with the Newton finish |
| `examples/01_legacy_marching/` | ladder base (decision below) |
| `examples/02_two_stage/` | ladder rung, signature option is the reconstruction |
| `examples/03_newton/` | ladder base for the physics rungs, signature option is the solver |
| `examples/04_newton_from_state/` | ladder rung, signature option is the restart |
| `examples/05_metals/` | ladder rung, signature option is `metals.inp` |
| `examples/09_spherical/` | ladder rung, signature option is the domain |
| `examples/10_warm_seed_ic/` | ladder rung, signature option is the IC |
| `examples/11_windae_ic/` | ladder rung, signature option is the Wind-AE IC |
| `examples/12_windae_ic_hd209/` | the same on HD 209458 b |

**The examples ladder (decided 2026-08-22).** `examples/01`--`12` is a ladder:
each folder is its base plus exactly one line, and `diff` against the base is
how the example documents its option. Atomic helium stays the ladder's
baseline, so `Include He23S? False` is kept in both bases and in every rung
whose signature option is something else. That is what leaves `06_he23s` a
one-line difference against `03_newton` rather than a folder that differs in
nothing. `06`, `07` and `08` carry the triplet, as they did before. Their
stored output is consistent with their input, so nothing here is recomputed.
Outside the ladder the triplet is on.

## Needs recomputing — inputs already flipped today

Their `input.inp` now says `True` while their stored `output/` predates the
change (2026-08-09).

| run directory | stored output | note |
|---|---|---|
| `examples/tutorial/` | 2026-08-09 | the minimal worked example |
| `examples/tutorial_nometals/` | 2026-08-09 | same config without `metals.inp`; flipped too so the pair keeps differing only in metals |

## Needs a decision, then flipping and recomputing

Still `False` today, outside the ladder, and **not** deliberate: the setting
was inherited when the case was created, and the case is about something else
entirely. Under the new default these should be `True`, but they carry
production numbers, so flipping them means recomputing them.

| run directory | stored output files | what the case is really about |
|---|---|---|
| `benchmarks/wasp121/` | 8 | WASP-121 b benchmark |
| `WASP-121b/solver_validation/hybrd1/` | 8 | solver comparison |
| `WASP-121b/solver_validation/newton/` | 8 | solver comparison |

The paper's WASP-121 b figures come from `WASP-121b/` proper (already
`True`), not from these.

## Flipped, nothing stored to recompute

`examples/13_lower_atmosphere/` ran two of its four planets with the triplet
and two without --- one example, one physics setup, an accidental split, and
one that mattered, since the He 2^3S + H2 Penning channel is active only
when a molecular run also tracks the triplet. All four now carry `True`
(`HD209458b/` and `WASP-121b/` were the two changed). None of the four has
stored output.

## Archive — leave alone unless asked

`backup/example_HD209458b/` is a stored configuration with its own
`EXHALE_setup.out` from a past run and no output directory; flipping the key
would only make that report disagree with the file next to it.

`backup/HD209458b_test/` holds eleven `False` directories with stored output
(`*_He0_*` in the names: helium deliberately absent or the triplet off in a
2026-06-era parameter study). They are a dated record of that study, not
live examples, and recomputing them would destroy the record without
replacing it with anything the current code needs.

## Not affected

- **Goldens.** Every `backup/regression/*/input.inp` sets `Include He23S`
  explicitly, so the default change moves nothing there. `make check` was
  5/5 byte-identical before this change and the compared cases are
  unchanged by it.
- **The 97 cases already `True`.** Their inputs and results are consistent.

## Side effects worth remembering

- **Parse corpus.** It snapshots `parse_dump.txt` for every `input.inp` in
  the tree, so each flipped input shows a diff until the corpus is
  re-snapshotted. That is a deliberate end-of-series action.
- **Example figures and notebooks** that read the flipped run directories
  will show the old numbers until those runs are redone.
- **Documentation of the ladder decision.** `examples/README.md` and
  Sects. 6.2--6.3 of the user manual now state that the ladder's baseline is
  atomic helium on purpose. Sect. 6.2 describes the tutorial, which carries
  the triplet, so its single-line comparison is now `Include He23S? False`,
  the change that removes it.

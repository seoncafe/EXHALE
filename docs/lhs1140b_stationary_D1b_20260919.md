# D1b: historical inputs pinned, and the default base cell width made exact

Item D1b of `docs/PLAN_20260918_rev2.md` section D1, approved by the user on
2026-09-19 ("User decisions of 2026-09-19", point 2). Done 2026-09-19 (KST) on
`lart4` with `GNU Fortran (conda-forge gcc 16.2.0-5) 16.2.0` and the OpenBLAS
of its prefix. Every number below is MEASURED unless marked READ.

---

## 1. Verdict

The contract D1 named is closed. An `input.inp` without the
`Base grid [dr,cells]:` key and one stating `Base grid [dr,cells]: 2.0e-4 50`
now build one grid, bit for bit, because the default is the double `2.0d-4`.
Before the default moved, every input that had run on the old default and has
results beside it (906 files) was given the old width explicitly, so none of
them changes grid: on the four cases measured, the pinned inputs on the new
build give the same data rows as the unpinned inputs on the entry build. A
state written on the old grid and handed to a run on the new default is
refused, and the refusal now says which width the run used and which line
restores the old grid.

## 2. The pinned line

```
# Base grid of the stored results: the Mixed-grid default before 2026-09-19 (docs/lhs1140b_stationary_D1b_20260919.md)
Base grid [dr,cells]: 1.9999999494757503e-4 50
```

placed right after `Grid type: Mixed`. `input.inp` accepts whole-line `#`
comments (`input_read.f90`, the keyword loop skips a line whose first
non-blank character is `#`; `refuse_duplicate_keys` and the unknown-line scan
skip them too). A comment after the value is not possible, because word 5 of
the key line is read as the cell count.

Why the string is right, MEASURED with gfortran 16.2 list-directed reads:

| quantity | bits |
|---|---|
| default-real literal `2.0e-4` widened to double (the old default) | `3F2A36E2E0000000` |
| read of `1.9999999494757503e-4` | `3F2A36E2E0000000` |
| read of `1.9999999494757503e-4 50` (the key line's words 4 and 5) | `3F2A36E2E0000000` |
| `2.0d-4`, the new default, and the read of `2.0e-4` | `3F2A36E2EB1C432D` |

The width is also the shortest decimal Python's `repr` gives for that double
(`0.00019999999494757503`). The grid it builds equals the old default grid bit
for bit: `base_grid_pinned_width.f90` compares every coordinate of the
production `define_grid` on its bits, at `r_max` 30 R_p and 10 R_p (section 5).

The cell count is 50 in every pinned file, READ: `N_low_cells` is set to 50
by the reader's defaults (`input_read.f90`, now `N_low_cells_default`) and the
`Base grid` key is its only other writer, so a file without the key resolves
to 50 whatever else it states. No file's count had to be guessed.

## 3. What was pinned, and what was not

`src/utils/pin_base_grid.py` (dry run by default, `--execute` to write,
`--list` to list every file) considers `input.inp` and `input_template.inp`
(the file a flux-closure rung copies into the `input.inp` of each iterate,
`element_flux_closure.py`), outside every `states/` and `runs/` directory, that
state `Grid type: Mixed` and no non-comment `Base grid` line.

Count, MEASURED (the files outside `states/` and `runs/`: 2445 `input.inp`,
2442 of them `Mixed`, 2430 of those without the key, plus 185
`input_template.inp`):

| class | area | pinned / noted |
|---|---|---|
| LIVE | `LHS1140b/models` | 770 pinned: 95 case-level inputs (cases, `seed/`, `diffusion_check/`), 46 closure iterates `k<NN>/`, 9 closure templates, 620 in the `.L*` and `.ab` experiment directories |
| LIVE | `backup/regression` | 51 (every matrix and named case) |
| LIVE | `vulcan_work` | 32 |
| LIVE | `WASP-52b` / `HD209458b` / `WASP-121b` | 23 / 5 / 3 |
| LIVE | `benchmarks` | 12 |
| LIVE | `examples` | 10 (below) |
| LIVE | `HD189733b` | 0: its input already states `Base grid [dr,cells]: 1.0e-4 100` |
| FROZEN | `LHS1140b/archive_20260830` | 1097 (939 inputs, 158 templates) |
| FROZEN | `LHS1140b/models_20260914_preL21` | 155 |
| FROZEN | `LHS1140b/models_20260915_db87` | 129 |
| FROZEN | `LHS1140b/examples` | 3 (the solutions of record of 2026-08-30, READ `LHS1140b/MODELS.md` section 1) |
| FROZEN | `backup/HD209458b_test`, `backup/phase_d_baseline`, `backup/lart_runs`, `backup/example_HD209458b` | 18, 11, 2, 1 |
| FROZEN | `backup/regression/_quarantined` | 5 (cases withdrawn from the matrix, each the record of a retired key) |
| FROZEN | `docs/audit_20260905`, `docs/lower_atmosphere_figs`, `docs/version_compare` | 2, 3, 4 (the runs behind documents) |
| TEMPLATE | `examples` | 14, not pinned (below) |
| TRANSIENT | `build/` (test work copies), `scratchpad/` (earlier work directories) | 87 and 178, not pinned: rewritten by whatever makes them |

After `--execute` wrote 906 files, a second dry run reports 0 files left in
every LIVE area, and a second `--execute` writes 0. Every pinned file differs
from its previous text by exactly the two added lines (checked on all 906
against a copy taken before the write).

Each FROZEN root carries a `GRID_DEFAULT_NOTE.md` stating that its inputs ran
on `1.9999999494757503e-4` R_p and 50 cells, that from 2026-09-19 the default
is `2.0e-4`, and the line to add to a copy to reproduce one of them.

The examples. Pinned: those with a stored state beside the input (`03_newton`,
`11_windae_ic`, `12_windae_ic_hd209`, `14_diffusion`, `15_molecular`,
`18_oxygen_chemistry`, `tutorial`, `tutorial_nometals`) and the two that the
`examples/README.md` recipe seeds from one of them: `04_newton_from_state`
(from `03_newton/output/`) and `16_molecular_metals` (from
`15_molecular/output/`); a seed and the run it seeds must build one grid.
Not pinned, taking the new default: `01_legacy_marching`, `02_two_stage`,
`05_metals`, `06_he23s`, `07_balmer_lya`, `08_full`, `09_spherical`,
`10_warm_seed_ic`, `13_lower_atmosphere/{HD189733b,HD209458b,WASP-121b,WASP-52b}`,
`17_lower_profile`, `19_molecular_ir_bands` (their `output/` directories are
empty or absent).

## 4. The default

- `parameters.f90`: `real*8, parameter :: dr_base_default = 2.0d-4`, the
  comment above it stating the present contract (the value, why an absent key
  and `2.0e-4` agree, that the older width is carried by the pinned key, and
  where the notes are). The cell count got a named default as well,
  `N_low_cells_default = 50`, used by both `parameters.f90` and the reader's
  default block in `input_read.f90`, which already read `dr_base_default`.
- `load_IC.f90`, `verify_restart_radii`: on a Mixed grid the refusal now prints
  the run's width at round-trip precision (`ES26.17E3`), its cell count,
  whether the width came from the key or the default, and the pinned line.
  This file was not named in the brief; the test the brief asks for (the
  refusal naming the width) cannot pass without it.
- `docs/input_schema.md` K33b: the trap text removed; the new default, the
  former one, the pinning and the notes stated.
- Found stale and corrected: `docs/EXHALE_BC_and_IC.tex` (said spelling the
  default does not reproduce a no-key run; PDF rebuilt), `README_HOWTO.md`
  ("the default reproduces the historical hardcoded grid"),
  `LHS1140b/models/write_reproduce.py` (its docstring and the paragraph it
  writes into every `REPRODUCE.md` said the default is the single-precision
  neighbor), `src/tests/grid_and_gates/README.md` (the shell tests do honor
  `EXHALE_TEST_OUT`).

## 5. Tests

RED on a build of the entry text (snapshot of `src/` at 08:23 KST, md5
`5a146a4823e2f035c56383ac14dee34b`), GREEN on the same snapshot with the three
changed sources (md5 `5221dc1a49f0836add439bbd40e0c298`):

| row | test | entry build | D1b build |
|---|---|---|---|
| `dr_base_default_is_the_read_of_2.0e-4` | `base_cell_width` | FAIL (`3F2A36E2E0000000` against `3F2A36E2EB1C432D`) | PASS |
| `pinned_width_is_the_pre_20260919_default` | `base_cell_width` | PASS | PASS |
| `base_grid_spelled_default_is_the_default` (no key vs `2.0e-4 50`, `r` column) | `base_grid_key` | FAIL | PASS |
| `base_grid_pinned_width_keeps_its_own_grid` (no key vs pinned) | `base_grid_key` | FAIL (identical) | PASS, 3.818119e-09 at row 366 of 504 |
| `default_grid_equals_key_2.0e-4_grid` at r_max 30, 10 | `base_grid_pinned` | FAIL (502 coordinates differ) | PASS |
| `pinned_grid_equals_pre_20260919_default_grid` at r_max 30, 10 | `base_grid_pinned` | PASS | PASS |
| `default_grid_differs_from_pre_20260919` at r_max 30, 10 | `base_grid_pinned` | FAIL | PASS: 6.57934e-09 at cell 342 (30 R_p), 5.02042e-09 at cell 357 (10 R_p) |
| `base_grid_pinned_state_loads_on_pinned_input` | `base_grid_restart` | PASS | PASS |
| `base_grid_old_state_refused_on_default` | `base_grid_restart` | FAIL (exit 0, loaded) | PASS (exit 1) |
| `base_grid_refusal_names_the_width` | `base_grid_restart` | FAIL | PASS |

The 6.58e-09 at 30 R_p reproduces L35's catalog-grid number (READ, L35 memo
section 8) independently.

The refusal as printed by the D1b build:

```
 (load_IC) ERROR: output/Hydro_ioniz_IC.txt holds cell centers of a different grid from the one this run built.
   this run: "Grid type: Mixed" with 500 cells;
   base grid: width  2.00000000000000010E-004 R_p x 50 uniform cells, the default (no "Base grid" key);
   a state written before 2026-09-19 on an input without that key needs
   "Base grid [dr,cells]: 1.9999999494757503e-4 50" (src/utils/pin_base_grid.py).
   cell 364: the run has r =  1.6997117887646787E+00
                    the file has r =  1.6997117822749761E+00
   worst relative difference  3.82E-09, tolerance 1.0e-10.
```

Suites run on both builds (every suite whose drivers call `define_grid` and
therefore see the default): `grid_and_gates` 291 PASS / 9 FAIL on the entry
build (exactly the nine rows above) and 300 PASS / 0 FAIL on the D1b build;
`adv_static_limit` 53/0, `boundary_state` 35/0, `element_operator` 40/0,
`ionization_stage_flux` 45/0, `physics_probe` 1604/0 on both.

## 6. Byte identity

Entry build on the unpinned inputs against the D1b build on the pinned
inputs, single threaded, scratch copies:

| case | steps | result |
|---|---|---|
| `wasp_full` | 400 | every output file identical except the two `# provenance:` lines of `Hydro_ioniz.txt` and `Hydro_ioniz_adv.txt` (run time stamp and `ck_input`, the checksum of `input.inp`); data rows identical |
| `mol_carrier` | 12000 (`maxsteps`) | same |
| `lower_profile` | 12000 (`maxsteps`) | same (including `element_flux_profile.txt` and `OI_levels*.txt`, identical) |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13`, `Restart intent: stationary evaluate` on its `latest_certified` generation `g0002_20260917T171710Z_4148d706`, prepared as `run_case.sh` prepares `eval/` | 0 (evaluate) | same; both exit 0 |

The records that change, and why: `EXHALE_resolved.out` row
`base_cell_width_source` reads `key` instead of `default`, and
`EXHALE_setup.out` says `(from the key)` instead of `(the default)`, with the
same width `1.99999994947575033E-004`; the run log gains the key's echo line.
This is the provenance flag doing its job for a pinned input.

## 7. The parse corpus and the generation manifests

Parse corpus. `run_parse_corpus.sh check` was NOT run as written: it runs the
tree's `EXHALE.x` in place in every directory holding an `input.inp` (its
exclusion list names `LHS1140b/exhale/`, which no longer exists), so it would
write `parse_dump.txt` and overwrite `EXHALE_setup.out` inside
`LHS1140b/archive_20260830/` and inside `states/` and `runs/`
directories. Instead the 108 corpus entries of `parse_golden/` that are LIVE
pinned directories were parsed (`EXHALE_PARSE_DUMP=1`) on scratch mirrors with
copies of each directory's files: entry build on the unpinned inputs, D1b
build on the pinned ones. 78 produced a dump on both, 30 (`vulcan_work`, all in
`parse_skip.txt`) on neither. The 78 dumps are byte-identical between the two
(`dr_base` is dumped with 16 digits and `N_low_cells` as an integer, and the
provenance flag is not dumped). All 78 already differed from `parse_golden/`
before this item, on rows unrelated to the grid (`R0`, `Mp`, `rho_bc`, `t_s`,
added keys such as `well_balanced`); no `dr_base` or `N_low_cells` row differs.
Nothing was refreshed.

Manifests. READ: `publish_state.py` writes the md5 of `input.inp` into each
manifest's `configuration_identity`, and the state header carries the
binary's own `ck_input`. No reader compares either with the current file
(`pick_seed.py`, `status.py`, `run_case.sh`, `run_closure.sh`,
`write_reproduce.py`, `import_legacy_states.py`, `compare_trees.py` searched;
`configuration_identity` is only written, `ck_input` only written).
`publish_state.py` states that it "records both and claims neither". So
pinning makes the recorded md5 of every published generation of a pinned case
differ from the present `input.inp`, and no reader refuses on it. Manifests
were not touched.

## 8. Found outside the scope, not changed

1. `LHS1140b/models/make_models.py` writes `Grid type: Mixed` and no
   `Base grid` line. A case it writes from now on builds the new default grid,
   while `run_case.sh` and `run_closure.sh` map seeds onto
   `models/current_grid_Hydro_ioniz.txt`, the old grid; `load_IC` would refuse
   the mapped seed. Proposed: write the pinned line after `Grid type: Mixed`
   in `make_models.py` (keeps the catalog on one grid), or deliberately move
   the catalog grid (regenerate the target grid file; a new-grid decision).
2. `backup/regression/run_parse_corpus.sh`: the exclusion list is stale
   (section 7) and the in-place run writes into frozen records and published
   generations. Proposed: exclude `archive_20260830/`, `models_20260914_preL21/`,
   `models_20260915_db87/`, and every `states/` and `runs/` directory, and
   refresh `parse_golden/` (stale on all 78 LIVE entries measured).
3. Any scratch copy of an unpinned input made before the pinning (for example
   other workers' copies) builds the new default grid on a binary built after
   this change and will refuse a catalog state; copy the input again.
4. `docs/ISSUES_20260909.md`, `docs/Update_EXHALE_stage2.md`,
   `docs/session_handoff_20260919.md` state the old default as the present one
   in records of their date; left as records.

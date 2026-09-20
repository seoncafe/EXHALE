# D1a: the base cell width a run resolved, disclosed and reproducible

Item D1a of `docs/PLAN_20260918_rev2.md` section D1, done 2026-09-18 (KST) on
the host `lart4` with `GNU Fortran (conda-forge gcc 16.2.0-4) 16.2.0` and the
OpenBLAS of its own prefix. D1a changes no number: the stored default keeps the
value it had, and the two regression cases run for this item are byte-identical
in every data line between a build of the entry text and a build carrying this
change. D1b, which moves the default to `2.0d-4`, is a separate item and is not
started here.

Every number below is MEASURED unless it is marked READ.

---

## 1. Verdict

The defect D1 names is a disclosure defect, and it is closed for reading:
a run now states which base cell width built its grid, to the digits that
rebuild it, and whether that width came from the key or from the code. The
value itself is untouched, so no grid in the tree moves.

Three things were true of the entry text and are no longer:

1. The default width was written twice, as the default-real literal `2.0e-4`
   in `parameters.f90` and again in `input_read.f90`. It is now one named
   constant, `dr_base_default`.
2. Nothing recorded whether `dr_base` came from `Base grid [dr,cells]:` or
   from the default. A logical `dr_base_from_key`, written only by the key's
   own branch of the reader, now does; value equality cannot answer it,
   because a key may state the default's own digits.
3. No output stated the width to a precision that reproduces it.
   `EXHALE_resolved.out` now carries the grid type, the cell count, the outer
   radius, the base cell width at round-trip precision, the uniform base cell
   count and the provenance flag; `EXHALE_setup.out` states the width and
   says "the default" or "from the key" in one line beside the rounded echo.

MEASURED, the numbers the disclosure is about:

| quantity | value |
|---|---|
| the stored default, `dr_base_default` | `1.99999994947575033E-004`, bits `3F2A36E2E0000000` |
| the double `2.0e-4` a key produces | `2.00000000000000010E-004`, bits `3F2A36E2EB1C432D` |
| relative difference of the two widths | 2.52621e-08 |
| max relative displacement of the cell centers, `roundtrip` grid (hot Uranus, 500 cells) | 3.818119e-09 at row 366 of 504 |
| max relative displacement, 500-cell LHS 1140 b catalog grid | 6.58e-09 (READ, L35 memo section 8) |
| the tolerance `load_IC` allows between a state's centers and the run's | 1e-10 relative (READ, `load_IC.f90` lines 2221 to 2244) |

The displacement is smaller than the width difference and differs from grid to
grid, because the uniform base block is rescaled onto a fixed `r_max`. That is
why the new test reports its own grid's number as an observation and gates only
on "identical" or "not identical".

## 2. The round-trip serializer

`ES24.17E3`, the first candidate, was rejected: MEASURED, it prints
`************************` for a negative value, 24 asterisks, and asterisks
parse as nothing. The production serializer is `round_trip_decimal` in module
`setup_report`, a field of 26 characters with `ES26.17E3`: 18 significant
digits, at least the 17 binary64 requires, and room for the sign of the value
beside a three-digit exponent of either sign.

MEASURED by `src/tests/grid_and_gates/base_cell_width_provenance.f90`, which
compares the `integer*8` transfer of the value written with that of the value
a list-directed read of the printed string returns:

| value | printed | bits written = bits read |
|---|---|---|
| 0 | `0.00000000000000000E+000` | yes |
| 1, -1 | `1.00000000000000000E+000` | yes |
| 1e-4, -1e-4 | `1.00000000000000005E-004` | yes |
| 2e-4, -2e-4 | `2.00000000000000010E-004` | yes |
| 1e-300, -1e-300 | `1.00000000000000003E-300` | yes |
| 1e300, -1e300 | `1.00000000000000005E+300` | yes |
| `dr_base_default`, its negative | `1.99999994947575033E-004` | yes |

The same driver asserts that `dr_base_default`, now written in
`parameters.f90` as the exact decimal `1.9999999494757503d-4`, carries the
bits of the default-real literal `2.0e-4` and differs from the double `2.0d-4`.
The spelling change is the reason that assertion exists: the exact decimal is
the same value only if the compiler rounds it correctly, so it is checked and
not assumed. MEASURED: both bit patterns are `3F2A36E2E0000000`.

## 3. The grid rows

Three assertions, in
`src/tests/grid_and_gates/base_grid_key_reproduces_default.sh`, on three
one-step runs of copies of `backup/regression/roundtrip` (`Grid type: Mixed`,
no `Base grid` key as shipped):

| run | key | result |
|---|---|---|
| A | none | `base_cell_width_source default`, width `1.99999994947575033E-004` |
| B | `Base grid [dr,cells]: 1.99999994947575033E-004 50`, read out of A's record | `base_cell_width_source key`; the `r` column BYTE-IDENTICAL to A's |
| C | `Base grid [dr,cells]: 2.0e-4 50` | the `r` column differs from A's, max relative difference 3.818119e-09 at row 366 of 504 |

Row B is the one that makes the record a reproduction instruction, and it is
also the one that shows why provenance has to be a flag: B states the
default's own digits and is still reported as a key.

RED before the change, MEASURED against a build of the entry text: the
Fortran driver does not compile (neither `dr_base_default` nor
`round_trip_decimal` exists), and the shell test stops at
`base_grid_resolved_record_present`, the entry text's `EXHALE_resolved.out`
carrying no `base_cell_width_Rp` row to state. GREEN after: every assertion of
both passes.

## 4. What no longer changes a number, and the check that says so

Two builds of the same snapshot of `src/`, differing only in
`parameters.f90`, `input_read.f90` and `write_setup_report.f90`, built the
same way (`make OBJDIR=build_D1a`, gfortran 16.2.0, OpenBLAS of its prefix):

| build | md5 |
|---|---|
| entry text | `4e02c6807587f0d88cdacb944a65ed25` |
| with D1a | `45d0aa6a4a8196be47c5e3d0b608747d` |

Both were run single-threaded on copies of two regression cases, and the two
output files the regression compares were compared directly:

| case | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|
| `wasp_full` | 400 | every data line identical; the file differs only in the `# provenance:` run timestamp | identical, whole file |
| `hydrostatic_column` | 300 | every data line identical; same timestamp line | identical, whole file |

No golden was written. `parse_dump.txt`, which the parse corpus compares, is
written by `write_parse_dump` and was not touched.

## 5. Where the numbers now are

- `EXHALE_resolved.out`: `grid_type`, `grid_cells`, `outer_radius_Rp`,
  `base_cell_width_Rp`, `base_uniform_cells`, `base_cell_width_source`
  (`default` or `key`), `base_grid_in_effect` (`F` for the `Uniform` and
  `Stretched` grid types, which ignore the width and the count).
- `EXHALE_setup.out`: the existing rounded echo, then
  `dr_base = <round-trip value> (the default)` or `(from the key)`.
- `REPRODUCE.md`, written by `LHS1140b/models/write_reproduce.py`: a
  `## The grid` section built from the run's own resolved record, with the
  `Base grid [dr,cells]:` line to state verbatim to put a run back on its
  grid. A record without those rows produces no section, so an older case
  directory is described as before rather than assigned the present default.
- `docs/input_schema.md`, row K33b: the stored default to full precision and
  the trap.

## 6. What D1a does not do

- It does not move the default. `2.0e-4` in the key still builds a different
  grid from an absent key; the difference is now stated everywhere instead of
  being silent. D1b is where the default moves, and the plan's precondition
  stands: a global default cannot tell an old input from a new one by its
  age, so every historical omitted-key input has to be pinned first.
- It does not touch the base face geometry. The `Mixed` construction ties the
  inner physical face to the first cell width, so refining the width moves the
  numerical boundary relative to the reservoir level; that is the separate
  item D1 names and nothing here addresses it.
- It does not regenerate any product. No archived state was read, mapped or
  rewritten.

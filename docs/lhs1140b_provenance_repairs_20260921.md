# Two records repairs, 2026-09-21

`docs/PLAN_20260920_rev9.md` leaves two records defects open: section 10.6, a
manifest that names the log its certificate was read from while nothing holds
that log still, and the section 16 item that the seed mode of a stored state
cannot be recovered once a solve has been taken. Both are repaired here. Both
are records defects: no state, residual or certificate value changes, and
nothing in either repair is a physical statement.

## 1. The certificate is recorded by its bytes, not by the path

`LHS1140b/models/publish_state.py`.

**What a manifest records now.** Every publication of a certificate writes,
inside `certification`:

| field | what it is |
| --- | --- |
| `source` | the log the block was read from, as before |
| `source_md5`, `source_bytes` | of that file as it was read at publication |
| `block_lines` | `[first, last]`, 1-based and inclusive, where the block stood in it |
| `block_file` | the copy of the block published inside the generation, `certification.txt` |
| `block_md5` | the md5 of that copy |
| `source_reused_by_later_runs` | true for the case-level `run.log` and `pp.log` that `models/run_case.sh` writes again at the next run of the directory |

All six fields are set on every manifest, `build_manifest` filling them with
nulls where a publication carries no certificate, so a null is a statement and
not a field that was forgotten. `verify` re-measures the published block
against `block_md5` and refuses a generation whose certificate came from a log
a later run replaces and that carries no copy of the block.

**What it does not do.** It does not touch a manifest already published. A
generation is immutable (`docs/PLAN_20260919_rev1.md` item P4c), so the three
atomic manifests of the P2 cases keep their missing fields, and what a reader
of those has to do is unchanged: locate the block by its own text and say in
which file it was found, which `LHS1140b/models/identity_record.py` does.
`verify` deliberately says nothing about them, since
`source_reused_by_later_runs` is absent there rather than false; making it
complain would turn every generation published before today red.

If the existing nulls are ever to be filled, the way to do it is a SEPARATE
provenance attachment beside the generation (the `provenance/` attachment
`models/recover_provenance.py` already writes, with its confidence field),
recording the file the block was actually found in, its md5 now, and the line
span of the match. That was not done here.

**Measured.** Published on a scratch copy of
`atomic_scalar_gj1132_kzz0/HeH2.6` outside the repository: the manifest
recorded `source_md5 c0d17a803be2b55c7f04ba70360d2edd`, `source_bytes 22996`,
`block_lines [238, 282]`, `block_file certification.txt`, `block_md5
edc0ea3c4fabfe2d655cad6459e76fa9` and `source_reused_by_later_runs true`. The
log was then overwritten with another run's. The generation still carries the
block at the recorded md5, and `identity_record.py` reads the state as

```
  md5 MEASURED      b39f4048349fbafac31fca5f145cc2f8
  md5 READ          c0d17a803be2b55c7f04ba70360d2edd
  verdict           MISMATCH: the log on disk is not the log the certificate was read from
```

which is the statement that could not be made before, the manifest of
`atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` generation
`g0004_20260919T094145Z_259fe9c3` carrying no `source_md5` key at all.

## 2. A solved state says which conversion produced its seed

`src/modules/init/molecular_seed_from_atomic_state.f90`,
`src/modules/files_IO/load_IC.f90`.

**What a state file records now.** Both halves of every state carry one
`# molecular_seed` line:

- `converted from <dir>; invariant <p|T>; x2_source <handoff|local|stated>`
  (with `; x2 <value>` in the stated mode) when this run built the seed;
- the same sentence, unchanged, when this run was started from a state that
  carries one, so it holds along a chain of restarts and a SOLVED molecular
  state states the mode and the invariant its seed came out of;
- `none: this state was not produced by an atomic-to-molecular conversion`
  for a state that was not seeded by a conversion;
- `not stated by the state this run was started from, which was written
  before this line existed` for a restart of an older state, which claims
  neither of the two.

The detailed `# molecular-seed-from:`, `# molecular_partition:` and
`# molecular_partition_local:` lines are still written by the converting run
alone, so the run that did the conversion stays distinguishable from the
states that descend from it. `docs/input_schema.md` appendix D carries the
same description.

**What it does not do.** It carries a sentence, not a state identity: the
line names the directory the atomic pair was read from and not that pair's
md5, so a reader who needs the exact source state still goes to the manifest
`seed` record of the generation. It says nothing about states written before
today, which is the fourth case above.

**Measured.** `src/tests/molecular_seed/run.sh` on the private binary: 26
assertions, all PASS. Its three conversions write `x2_source handoff`,
`x2_source local` and `x2_source stated; x2  0.0000000E+00`, and the atomic
state they were built from writes the `none` sentence. The `reload` directory of that
suite, a run that loads a seed and solves, writes
`# molecular_seed converted from .../atomic/output; invariant p; x2_source
handoff` into the state it solved, which is exactly what phase 8 of
`docs/lhs1140b_m2_conversion_audit_20260921.md` had to recover from two
`seed.log` files and the provenance record.

## 3. The regression

The state files' data rows are unchanged: only a `#` header line is added. A
private run outside the repository, on a copy of the harness and of the four
cases, with the delivered binary:

```
REGRESSION_EXE=<private build> run_check.sh check hp_front roundtrip mol_base_handoff mol_carrier
==> REGRESSION PASS (all cases byte-identical)
```

Every one of the sixteen file comparisons reads `PASS <file> (data
identical)`. The four cases ran on a build of this source whose only later
change was comment text; `hp_front` and `roundtrip` were then re-run on the
delivered binary (md5 13e981692643) with the same verdict. The harness compares numeric content alone (`grep -v '^ *#'`), so
it reports nothing about the header line either way; that the line is written
is shown by the suite above and by the case outputs, where `hp_front` carries
the `not stated` sentence, `roundtrip` the `none` one.

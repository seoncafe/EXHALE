# D8 steps 3 to 7: the generations a case publishes, the index a reader resolves, and the states this tree already held

2026-09-18. Item D8 of `docs/PLAN_20260918_rev2.md`, steps 3 to 7, in the
shape approved that day (`docs/DECISION_D2a_D8_review.md` sections 3 and 4:
the Codex structure of `docs/DECISION_D2a_D8_codex.md` section 4 as reviewed,
with legacy import of the existing products as identified generations and
with D8a as the first increment). Item D8a
(`docs/lhs1140b_stationary_D8a_20260918.md`) is unchanged by this and its
fixture still passes.

Files: `LHS1140b/models/publish_state.py` and
`LHS1140b/models/import_legacy_states.py` (new),
`LHS1140b/models/legacy_state_map.txt` (new),
`LHS1140b/models/run_case.sh`, `run_campaign.sh`, `pick_seed.py`,
`status.py`, `write_reproduce.py`,
`LHS1140b/models/tests/state_generations.sh` (new),
`tests/seed_compatibility.sh` (its N check),
`src/utils/map_state_to_grid.py`, `src/utils/element_flux_closure.py`
(reader resolution only), `LHS1140b/MODELS.md` section 9 (the contract, and
the two new tools in the section 1 layout). No file under `src/modules/` was
changed, no Fortran was changed, no golden was touched and no state was
recomputed. `LHS1140b/models/` is outside the git remote, so of the files
above only `MODELS.md` and the two under `src/utils/` are tracked.

---

## 1. What was built

`MODELS.md` section 9 is the contract as implemented and is not restated
here. In one paragraph: a state directory (one `input.inp` and one
`output/`) publishes each of its states as an immutable
`states/<generation_id>/` holding the state files, the certification block
and a `manifest.json`; `state_index.json` beside it carries
`latest_complete` and `latest_certified` as separate fields; one publisher,
`models/publish_state.py`, writes both, under a lock, by writing the
generation into an unpublished directory, validating what arrived, making it
read-only, and then replacing the index with an atomic rename on the same
filesystem; `runs/<run_id>/` holds the inputs of one attempt and a `run.json`
naming its logs and the generations it left; and every reader resolves the
index ONCE and takes both halves of a state from the one generation it names.

Three choices are worth naming here because they are not in the approved text:

**The generations are copies, not hard links.** The binary opens
`output/Hydro_ioniz.txt` with a truncating open, so a link to that file would
be truncated with it: a link is not an immutable generation. MEASURED after
the import of section 3: the 134 `states/` directories take 132 MB together
(`du -c`), against 565 MB in the `output*` directories they were read from,
and the whole `models/` tree went from 1992 MB to 2124 MB.

**A generation holds the state, not what was derived from it.** The
advection-corrected profiles, the heating and cooling breakdowns and the
transit curves are derived from a state and do not define one; they stay in
`output/` and the manifest lists them by name and md5. The logs stay where
they are too, and what travels into the generation is the certification block
(`certification.txt`), which is what the certificate is.

**The state's own claim is the anchor of the certificate.**
`write_output.f90` writes `certified=T|F cert_reason=...` into the state
header, so the claim is part of the bytes the manifest hashes and cannot
drift from them. `latest_certified` moves only when the header claims
`certified=T` AND a certification block says CERTIFIED AND the attachment of
that block to these bytes is established. Completion is a separate field and
a separate reference.

## 2. What the readers do now

| reader | before | now |
|---|---|---|
| `map_state_to_grid.py` | `name.txt` if present, else `name_IC.txt`, decided for each half on its own | resolves the index once and takes both halves from one generation; with no index, both halves must carry the same suffix, and a `Hydro_ioniz` of one generation beside an `Ion_species` of another is refused |
| `element_flux_closure.py` | the previous iterate's `output/` as it stands | resolves that directory the same way, through the mapper's own resolver, and logs which generation it took |
| `pick_seed.py` | `certified()` from the `REPRODUCE.md` line, else the log; the path printed was `<case>/output` | the index's `latest_certified` where there is one; the path printed is the generation directory, so the mapper is handed one immutable state |
| `status.py` | the solver verdict and the certification from `run.log` | the index and nothing else for a case that has one: the published generation's solver verdict, ending and certification. A case with no index says `[no index]` |
| `write_reproduce.py` | nothing about generations | names `latest_complete`, `latest_certified` and every generation of the case |
| `run_case.sh` | the evaluate pass read `output/`; the molecular seed read `<atomic>/output/*_IC.txt` | the evaluate pass reads the generation the solve published; the molecular seed is the atomic case's `latest_certified` generation, staged under `runs/<run_id>/molecular_seed_source/` with the names the conversion reads |
| `run_campaign.sh` | the status line ended at the ceiling | it also carries `published=<generation> certified=<generation>`, two fields because completion is not certification |

An EVALUATION carries no solver verdict of its own, so a `status.py` row
whose published generation is an evaluate pass reads `evaluated certified`
rather than an `info` that belongs to the solve. That is what the catalog
looks like today: 81 of the 134 published states were written by an evaluate
pass, and the solve that produced the state they were made from is not on
disk for any of them.

## 3. The states this tree already held

`models/import_legacy_states.py` entered them without recomputing anything.
MEASURED 2026-09-18: 136 state directories (90 cases and 46 flux-closure
iterates), 222 generations published, 134 indexes (two directories hold no
state pair), 86 indexes with a `latest_certified` and 48 with none. A second
run publishes nothing and leaves every index byte for byte. Nothing outside
`states/` and `state_index.json` was written: no file of the catalog older
than this item has a modification time inside it.

**How a certificate is attached, and why it is a surrogate.** No run of this
tree recorded the md5 of the pair it certified, so the contract's "the exact
state identity the certificate assessed" cannot be established for an
imported state by the means the contract names. What is used instead, and
named as such in `certification.identity_basis` of every manifest: the state
carries the second it was written in, in its own header
(`# provenance: ... run=`, `date_and_time` at the moment the header is
written, `utilities.f90`), and the log of the pass that wrote it is the one
that lies nearest that second AND whose verdict is the verdict the state's
own header carries. Where no log of the directory satisfies both, the
generation is complete, uncertified, and says that no log there states a
verdict about it. MEASURED: 135 of the 222 have the attachment established
and 87 do not; 84 of those 87 are `output_pre_L34/` states, which the old
post-processing route wrote with one CFL step and `certified=F` over a solve
that had certified, and whose own log was overwritten by the run that moved
them aside.

The agreement test is not decoration. Without it the nearest log wins, and on
`molecular_scalar_gj1132_kzz1e9/HeH2.13` the nearest log is `run.log`, which
ended 0.16 s before the state was written and says CERTIFIED, while the state
carries `certified=F cert_reason=failing_entries` and was written by the
evaluate pass whose `pp.log` ended 1.9 s after it and says NOT CERTIFIED.

**Two catalog cases lose a certified status they appeared to have.**
MEASURED over the whole catalog: `REPRODUCE.md` states CERTIFIED for 88
cases, the index sets `latest_certified` for 86, and the two that differ are
`molecular_scalar_gj1132_kzz1e9/HeH2.13` and `HeH9.7`. Their `REPRODUCE.md`
states the certification of the SOLVE; the pair in their `output/` was
written over that solve by the L34 evaluate pass and carries
`certified=F cert_reason=failing_entries`. These are the two molecular
reference states that D5b-2 found refused under boundary model v1, reached
from an independent direction; their manifests also carry
`stale_under: boundary_model_v1` with the D5b-2 memo named. `pick_seed.py`
no longer offers either as a certified seed, which is the correction it
looks like and not a regression: the state on disk is not the state the
certificate was printed for.

**A state of another directory becomes a generation of a case only where
`legacy_state_map.txt` names it**, with the reason, one line each. A run
merely SEEDED from a case's state is not a generation of that case, and the
import found three such relations in this tree of which only one is a
continuation of the case's own solve.

## 4. The two defects of `docs/lhs1140b_stationary_D8bound_20260918.md` section 7

**`molecular_scalar_gj1132_kzz1e9/HeH0.083` published the state item L33's
forty passes started from.** The import measures it and says so: the case's
`output/Hydro_ioniz.txt` has md5 `1229b961323d7ac710146650da72454a`, which is
byte for byte `models/.L22/i4_kz0083/output/Hydro_ioniz_IC.txt`. That state
is now the case's first generation, `g0001_20260916T023329Z_32b34e57`, and
its manifest names the run whose seed it is. The state those forty passes
LEFT, md5 `1ee5430a3c8560c40f017bc5eb618a8b`, is the second generation,
`g0002_20260916T220756Z_f9fca548`, with the first as its parent and the map
line as its provenance; it is what `latest_complete` now names, so a reader
of that case gets the state the run left. Both are NOT CERTIFIED and
`latest_certified` is null: the contract does not invent a verdict for
either.

**`molecular_scalar_gj1132_wellmixed/HeH0.083` has no `ENDING` file.** Its
manifest states the ending class its log gives under `classify_ending` of
`run_case.sh` and says in the `ending.source` field that the class was
classified from the log because the directory holds no `ENDING`. No case of
the catalog carries one yet, since none has been run since item D8a landed;
the field names the log for all 222 generations.

## 5. A defect of the ending classifier, found and fixed

`run_case.sh`'s `certification_block` looked for the block headed
`final state, as written` and fell back to the LAST `(certification)` line of
the log. The block has two headings, one for each routine that writes it:
`final state, as written` (`EXHALE_main.f90`, both write sites) and
`FINAL STATE AS WRITTEN -- this run is NOT a certified stationary solution`
(`certification.f90`, `certification_stop_uncertified`, the route of a run
that exits 2). On a log of the second kind the fallback took the one-line
exit message `the state was written in full; the run exits with status 2`,
which measures no row, so `classify_ending` returned `no_row_measures` for a
log whose block names the rows that refuse it, and `no_row_measures` is an
ending the seed walk answers with another seed. Both headings are now
recognized, in `run_case.sh` and in the publisher that reads the same block.
MEASURED: no `run.log` of the catalog carries the second heading today (86
carry the first), but `pp.log` files do.

## 6. Tests

`models/tests/state_generations.sh`, MEASURED 2026-09-18: 27 checks, 27
passed. No solve is run; the states are small synthetic pairs carrying the
header fields the contract reads.

| check | what it establishes |
|---|---|
| S1a (RED) | an index written in place is left half written and a reader refuses it |
| S1b to S1e (GREEN) | a publisher stopped by PID between the temporary write and the rename leaves the PREVIOUS index, complete and readable, a reader still resolves the generation it names, and the temporary file is left behind and is not the index |
| S2a to S2d | a pair whose halves carry different row counts is not complete, and the resolver, the mapper, the seed chooser and the status reader all refuse it |
| S3a to S3d | an uncertified publication moves `latest_complete` and not `latest_certified`, and a CERTIFIED block over a state whose header claims `certified=F` certifies nothing |
| S4a, S4b | a second publisher on one case is refused and writes no generation |
| S5a to S5c | a published generation's files and directory are read-only and `verify` re-measures them against the manifest |
| S6a to S6c | the mapper reads the generation the index names; a `Hydro_ioniz` of one generation beside an `Ion_species` of another is refused (GREEN) where the rule it replaced, which searched for each half on its own, took exactly that pair (RED) |
| S7a to S7c | the legacy import is idempotent: the second run publishes nothing and leaves the index byte for byte |
| S8a, S8b | the two catalog defects of section 4, read back out of the catalog's own index and manifests |
| S9a | the runner publishes a generation and leaves the record of its attempt, exercised through `run_case.sh`'s own functions |

The fixtures of D8a and D9 after these changes, MEASURED 2026-09-18:
`tests/run_case_policy.sh` 16 of 16, `tests/campaign_budget.sh` 12 of 12,
`tests/closure_continuation.sh` 10 of 10, `tests/seed_compatibility.sh` 11 of
11. The last one is 11 and not 10 because its N check was rewritten: it now
compares the BYTES of the state each of five real cases chooses, against the
recorded choice, rather than the path string, since the tool now prints a
generation directory; and a new N2 states the sixth case's correction of
section 3 as what it is.

## 7. What is not claimed

- **No crash durability.** A rename is atomic within a filesystem and that
  is what is tested here, by stopping the publisher between the temporary
  write and the rename. What a server failure would leave on this NFS is not
  tested and is not claimed.
- **No exact algorithmic continuation.** The contract supports a warm
  restart, stated in `MODELS.md` section 9.6: the conserved primitives and
  the species are authoritative, the lower ghost rows are re-derived under
  boundary model v1, and every solver control is restarted.
- **No `best` reference.** Ranking two states needs a rule for which is
  better and that rule is not written.
- **No state was recomputed and no verdict was invented.** Every certificate
  in the 222 manifests was read out of a log that this tree already carried.

## 8. Noticed outside the scope, reported and not acted on

1. **The catalog's published states are evaluate-pass states.** For 81 of the
   134 indexed directories the published pair was written by the no-step
   evaluate route, and the solve that produced the state it was handed is not
   on disk: the old post-processing route wrote over it. That is why those
   rows of `MODELS.md` section 7 can no longer show the solver's `info` and
   `||R||` from the index. New runs publish the solve and the evaluation as
   two generations with the first as the parent of the second, so the chain
   exists from here on.
2. **The `_IC` pair of a case directory is not a stable reference.** It is
   the mapped seed for a run that stopped before its evaluate pass and the
   solved state for one that did not, and `molecular_seed_from_atomic_state`
   reads exactly `<directory>/<name>_IC.txt`. This item stopped the runner
   from depending on it (section 2), but any other reader of a `*_IC.txt`
   file in this tree is reading whichever of the two that directory happens
   to hold.
3. **`atomic_scalar_gj1132_kzz1e9/HeH2.13` was re-measured from a pair that
   is not its own solved state.** Its `output/Hydro_ioniz_IC.txt`
   (md5 `96e13cc1d27a39ec3a58351583a727f4`) is byte for byte
   `output_pre_L34/Hydro_ioniz_IC.txt`, the seed of the earlier run, and not
   `output_pre_L34/Hydro_ioniz.txt` (md5 `28d59f00995eed6e14ef2f23ad963546`),
   the state that run wrote. What the L34 evaluate pass certified is the pair
   it WROTE, which is what its generation carries, so nothing in the index is
   affected; what is worth checking in L34's own record is which state it was
   handed.
4. **222 generations, 84 of which are `output_pre_L34/` states with no log
   of their own.**
   They are preserved and complete, and they carry no verdict. Whether the
   catalog wants to keep them is a retention question, and a retention policy
   is a separate decision; nothing here deletes anything.

## Addendum 2026-09-19: four runner defects found by the D9 step 3 refresh

The catalog refresh (`docs/lhs1140b_catalog_refresh_20260919.md`) re-measured
86 cases under `LHS1140b/models/EXHALE_7670f310.x` (md5
`7670f31031fb4db91d27b44cb0da6f70`) and reported, outside its scope, four
defects of the runner built here. All four are fixed; the rules are written
in `LHS1140b/MODELS.md` sections 9.3, 9.5 and 9.8.

1. **Read-only `_IC` pairs.** `cp` carries the mode, so the `_IC` pair the
   evaluate pass copies out of a read-only generation stayed read-only in the
   case's `output/`, and the next seeded run stopped in
   `src/utils/map_state_to_grid.py` with `PermissionError`. `run_case.sh`
   now makes that pair writable, and the atomic pair it hands the molecular
   seed conversion likewise. Those two are the only places the runner copies
   out of `states/`; the continuation copies from `solve_<n>/output/`, which
   is not published and is writable.
2. **`pick_seed.py` tier 0.** `certified()` answered from `latest_certified`
   while the path printed was resolved from `latest_complete`, so a refused
   evaluation published after a certified solve was handed over as "own most
   recent certified state". Every candidate from this tree is now resolved
   from `latest_certified` where the index names one (metadata, the
   state-holding test and the printed path), and the tier-0 line says which
   generation it handed over.
3. **No evaluate-only entry, and `rm -rf eval`.** `run_case.sh --evaluate
   <case>` evaluates the index's `latest_certified` in `runs/<run_id>/eval/`
   and publishes the written state as its `evaluate` child. Both routes now
   run the pass there, the removal of `<case>/eval/` is gone, and every
   case-level file the pass or the transit synthesis replaces is first copied
   into `runs/<run_id>/superseded_case_products/`.
4. **A refused parent stayed `latest_certified`.** The publisher now records
   `stale_under: <binary md5>` and `stale_evidence` on the index entry of a
   certified parent that an `evaluate` child refuses, and moves
   `latest_certified` off it (to the newest other generation an evaluate pass
   of that binary certified, else null). `publish_state.py reassess` applied
   this to the three refused cases of the refresh,
   `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`, `.../HeH9.7` and
   `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`: `latest_certified` is
   null in all three (MEASURED; `verify` passes).

`write_reproduce.py` also stated "the run loaded the state already in
`output/`" for an evaluate pass given no seed; under `--evaluate` the record
now names the generation evaluated, the run directory and the child
published. The 74 `REPRODUCE.md` files of the refresh that carry the old
sentence are not rewritten here (MEASURED count, `grep -l`); the refresh memo
states (READ) that each carries a note naming the parent, run and generation.

A fifth defect was found while testing: `publish_generation` of
`run_case.sh` wrote its progress note to standard output inside the
`$(...)` that reads its value, so every generation id the runner passed on
(the `--parent` of the solve route's evaluate child, `solve_generation` and
`evaluate_generation` of `run.json`) carried the note in front of the id.
The note now goes to standard error. No catalog manifest carries such a
parent (MEASURED: none of the published manifests has a parent id containing
the note), because the refresh published through its own driver.

Tests: `LHS1140b/models/tests/evaluate_entry.sh`, 18 checks on a synthetic
tree with a stand-in binary: 18 of 18 pass on the tools as they stand, and 1
of 18 (E4e, `verify` on a demoted index, which the earlier text also passes)
on copies of the tools taken before these changes. The existing fixtures
after the change, MEASURED 2026-09-19: `state_generations.sh` 27 of 27,
`seed_compatibility.sh` 11 of 11, `campaign_budget.sh` 12 of 12,
`run_case_policy.sh` 16 of 16, `closure_continuation.sh` 10 of 10.

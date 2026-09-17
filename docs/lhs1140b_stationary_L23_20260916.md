# L23: the certification pair of a restart header, under a stated metadata contract

Item L23 of `docs/PLAN_20260916_rev3.md` section 1 and of
`docs/PLAN_20260913_lhs_stationary.md`, with rows R1, R39 and R51 of the rev3
review tables.

Binaries, both built from the tree of 2026-09-16 (git `43bc28c`, working copy)
and differing only in `src/modules/files_IO/load_IC.f90` and
`src/EXHALE_main.f90`:

| | md5 |
|---|---|
| control, the entry text of those two files | `c6fb9e8a0850a5a0684aff5525315580` (`EXHALE_L23ctl.x`) |
| measured | `f23da8d84fafc2f0587b610add56e714` (`EXHALE_L23.x`) |
| the tree binary of the catalog | `c2e9c9990b9f14f1be8cd77abca68945` (`EXHALE.x`) |

## 1. The defect

`write_coupling_state_header` (`utilities.f90` lines 167 to 252) writes the
certification of a state as a PAIR: `certified=T|F` from `state_is_certified`
and, when it is not empty, `cert_reason=<token>` from
`state_certification_reason`. Both are set together and only by
`set_state_certified` (`utilities.f90` lines 31 to 38), which
`certification.f90` lines 1357 to 1369 call with `.true.` plus
`wind_certification_token(rep)` (`certified_in_wind` when species rows were
judged in the wind alone, empty otherwise), or with `.false.` plus
`no_stationary_claim` or `failing_entries`.

`parse_coupling_header` (`load_IC.f90`) read `certified` and had no
`cert_reason` case at all. A state file reloaded and written back without a
new measurement therefore lost, in the control binary:

* the reason, which is what `backup/regression/roundtrip/roundtrip_check.sh`
  reports (MEASURED, control binary):

```
  FAIL coupling header not preserved:
    ref  coupling: sec_ion=F sec_ion_step=-1 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
    dump coupling: sec_ion=F sec_ion_step=-1 recon=PLM certified=F mode=init
```

* and, where the file claimed one, the certification itself. MEASURED with the
  control binary on the certified state of `backup/regression/wasp_full_newton/IC/`
  reloaded with `EXHALE_DUMP_IC=1`: the file states
  `certified=T`, the written file states `certified=F`. `ic_certified` was read
  into the loader and used only to tell the stationary evaluation what claim it
  was answering (`certification_note_stationarity_claim`, `EXHALE_main.f90`);
  nothing put it back into the state metadata, so a raw reload reported "no
  certification was ever made" of a state that had been certified.

## 2. The contract

The pair a file states is metadata OF THE IMPORTED STATE. It is held in
`ic_certified` and the new `ic_cert_reason` (`load_IC.f90`), kept distinct from
`state_is_certified` and `state_certification_reason`, which are what THIS
executable writes after it has evaluated a state.

* A reload that writes the state back WITHOUT evaluating it carries the
  imported pair through unchanged. The one such route is the `EXHALE_DUMP_IC=1`
  dump, which writes the state as it was loaded, before the first equilibrium
  sweep, and stops; it is where the pair is put into the state metadata, next
  to the line that does the same for the secondary-ionization coupling. Every
  other route that writes after a reload measures the state first: the
  marching stop, the "Do only PP" stop and the stationary solve all reach the
  final write through `certification_evaluate`, and `Restart intent: stationary
  evaluate` measures the state it just read.
* An evaluation writes the pair it determined. No success token survives a
  change of the state or of the model, because `set_state_certified` is called
  by the evaluation with what it measured.
* The parser does not call `set_state_certified` and sets no state variable;
  it records what one file says.
* The Boolean comes from `certified=` alone. A `cert_reason` token this
  version does not know is provenance text, carried and never promoted to a
  certification.

## 3. The policies

Read by `parse_certification_claim` into one record per file
(`file_certification_claim`: was a `# coupling:` header seen at all, was each
of the two fields stated, and what it said), and resolved by
`adopt_certification_claim` after both halves of the state have been read.

| situation | what happens |
|---|---|
| `cert_reason` longer than the destination (32 characters, `len(state_certification_reason)`) | the load is refused, with the offending line printed. The length is measured on the token BEFORE it is copied, so what is refused is the token the file carries and not a cut version of it |
| the same key stated twice in one file's header | refused, with the line printed |
| `certified` or `cert_reason` stated by both files with different values | refused as two claims about one state, never merged |
| a field one half states and the other omits, both halves carrying a header | the stated value is taken and a note says so |
| `Ion_species_IC.txt` carrying no `# coupling:` line at all | the ordinary layout, since the writer states the claim on `Hydro_ioniz.txt` alone; silent |
| a token this version does not know | kept as provenance text; `certified=` alone decides the Boolean |

The claim of `Ion_species_IC.txt` is read for the first time here. The writer
puts the coupling line on `Hydro_ioniz.txt` only, so for every pair this code
produces the second half states nothing and the resolution is the first half's
own; the reader exists so that a pair whose two halves DO state different
things is refused instead of silently taking one of them.

## 4. The tests

`src/tests/grid_and_gates/restart_intent_and_metadata.sh`, sixteen new rows
(the sixteenth is section 6a).
Each row hands the certified state of `backup/regression/wasp_full_newton/IC/`
back with its `# coupling:` header edited and reloads it with
`EXHALE_DUMP_IC=1`, except the last, which evaluates a state under a changed
stellar EUV luminosity.

MEASURED, control binary (RED) and measured binary (GREEN), same script:

| row | control | measured |
|---|---|---|
| `cert_pair_true_with_token` | FAIL, written `certified=F` with no reason | PASS |
| `cert_pair_true_without_token` | FAIL, written `certified=F` | PASS |
| `cert_pair_false_no_stationary_claim` | FAIL, reason dropped | PASS |
| `cert_pair_false_failing_entries` | FAIL, reason dropped | PASS |
| `cert_pair_false_without_reason` | PASS (the one pair the control could carry) | PASS |
| `cert_pair_independent_of_field_order` | FAIL | PASS |
| `cert_pair_reason_at_the_field_length` (32 characters) | FAIL | PASS |
| `cert_pair_unknown_token_kept` | FAIL | PASS |
| `cert_unknown_token_does_not_certify` | PASS | PASS |
| `cert_reason_longer_than_the_field_refused` (33 characters) | FAIL, exit 0 | PASS, exit 1 |
| `cert_duplicate_key_refused` | FAIL, exit 0 | PASS, exit 1 |
| `cert_conflicting_pair_refused` | FAIL, exit 0 | PASS, exit 1 |
| `cert_omitted_field_taken_from_the_other_half` | FAIL | PASS, with the note |
| `reevaluation_writes_its_own_pair` | PASS | PASS |

The three refusals, MEASURED with the measured binary:

```
  p_over32     (load_IC) ERROR: output/Hydro_ioniz_IC.txt: the "cert_reason" token is longer than the field.
  p_duplicate  (load_IC) ERROR: output/Hydro_ioniz_IC.txt: the key "certified" is stated twice.
  p_conflict   (load_IC) ERROR: the two restart files state different stationary claims about one state:
```

The last row is the one that fixes the direction of the contract: the file was
handed over stating `certified=T cert_reason=certified_in_wind` and the run
read it under a stellar EUV luminosity twice the one the state was solved at.
The evaluation answered with what it measured, `certified=F
cert_reason=failing_entries`, and not with the token it was given. It passes
with both binaries, which is the point: the imported pair is carried where
nothing measures the state and is overruled where something does.

Every other row of the script passes with both binaries, including
`stationary_evaluate_certification_unchanged`, which reads the same field on
the evaluation route. The whole `grid_and_gates` suite with the delivered
binary: 220 PASS and one FAIL,
`outer_iteration_ending_is_the_stagnation_one` of
`output_state_consistency.sh`, which fails identically with the control binary
and belongs to whatever moved the outer iteration, not to this item.

## 5. The round trip

`backup/regression/roundtrip/roundtrip_check.sh` on a scratch copy of the case,
stage A re-run with each binary (MEASURED):

```
control  FAIL coupling header not preserved:
           ref  coupling: sec_ion=F sec_ion_step=-1 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
           dump coupling: sec_ion=F sec_ion_step=-1 recon=PLM certified=F mode=init

measured ok   coupling header preserved: coupling: sec_ion=F sec_ion_step=-1 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
         restart round trip is the identity to 1e-12
```

The case golden is the stage-A output of a COLD start, which loads no state and
is not touched by this change; no golden is refreshed for it.

## 6. The impact

No output number moves, and the changed code is not even reached by a run that
starts cold: `backup/regression/wasp_he23off` has `Load IC? False`, so
`load_IC` is not called and the `EXHALE_DUMP_IC` branch is not entered.
MEASURED on scratch copies of that case, `OMP_NUM_THREADS=1`, the tree binary
`c2e9c9990b9f` against `EXHALE_L23.x`:

| what was compared | result |
|---|---|
| `output/Ion_species.txt`, `output/Hydro_ioniz_adv.txt`, `output/Ion_species_adv.txt` of a run bounded at `EXHALE_MAXSTEPS=200` | byte-identical, whole file |
| `output/Hydro_ioniz.txt` of the same run | every data row identical; the files differ in one header character, the `run=` wall-clock stamp of the `# provenance:` line (`19:15:23` against `19:15:22`), which every invocation writes afresh. The `# coupling:` line is identical, `certified=F cert_reason=no_stationary_claim mode=init` on both |
| the marching log of the unbounded run, line by line | identical over 1484 lines, 1410 accepted steps, every printed `du` and `dtu` digit |
| `output/Hydro_ioniz.txt` and `output/Ion_species.txt` written by the unbounded runs at the same step | byte-identical, headers included |

The `cert_reason` field of the header does NOT change for this case: a cold
start writes what its own certification determined, which is what it wrote
before. The field changes only in a file written by a reload that does not
measure the state, where the control binary wrote no reason at all, and where
`certified=` itself was wrong whenever the imported state carried a claim.


## 6a. The line that is not a coupling header

The selector that chose which comment line to parse was
`index(line,'coupling:') > 0`, a SUBSTRING test. `map_state_to_grid.py` writes
two lines into a mapped seed: the state's own

```
# coupling: mode=init t_phys=0 certified=F sec_ion=T sec_ion_step=0 recon=PLM cert_reason=no_stationary_claim
```

which says what the seed IS, an initialization state on the new grid, and

```
# mapped-from-coupling: sec_ion=T sec_ion_step=0 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
```

which says what the SOURCE state was produced under. The second is provenance
of the mapping, not a statement about the state in the file, and the substring
test matched it as well.

**A pre-existing defect, exposed by the duplicate-key policy.** With the
substring test the loader parsed BOTH lines, in file order, and the later one
won every key it carried: `sec_ion`, `sec_ion_step`, `recon`, `certified`,
`mode` and `t_phys`. A mapped seed was therefore restarted with the SOURCE
state's armed step, reconstruction label and run state, which is exactly what
the mapper's own line was written to prevent (it states `mode=init t_phys=0`
because a mapped state is a seed and not a continuation of a trajectory).
MEASURED with the entry-text binary on a state whose provenance line was given
`sec_ion_step=99` against the state's own `sec_ion_step=2429`: the reload wrote
`sec_ion_step=99`. A seed mapped from a `mode=phys` state would in the same way
have claimed to stand on a trajectory it was never on.

The duplicate-key refusal of section 3 made that silent adoption loud: with
both `certified=` fields visible in one file, the load stopped and printed the
line it stopped on. MEASURED on
`LHS1140b/models/.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output/`:

```
 (load_IC) ERROR: output/Hydro_ioniz_IC.txt: the key "certified" is stated twice.
   # mapped-from-coupling: sec_ion=T sec_ion_step=0 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
```

**The fix.** Both selectors now match the line's OWN label,
`comment_field_is(line, 'coupling:')`, which is the same first-token test the
reader already used for the `columns` line and for the same reason: a field of
a header is its first token and not a substring of the line.
`mapped-from-coupling:` is provenance and is never parsed as the state's
coupling header, so the state's own line is the only one read and the pair and
the fields it states are the ones adopted.

Test row `mapped_from_provenance_is_not_a_coupling_header`: a header carrying
both lines, the provenance one stating a DIFFERENT pair and a different armed
step. MEASURED, one script, three binaries:

| binary | result |
|---|---|
| entry text | loads, and writes `certified=F cert_reason=` (none) `sec_ion_step=99`: the provenance line's fields, the pre-existing defect |
| the substring selector with the duplicate-key policy | refuses the load ("the key \"certified\" is stated twice") |
| the label selector | loads and writes `certified=T cert_reason=` (none) `sec_ion_step=2429`: the state's own line |

The mapped seed above then loads and is evaluated end to end (MEASURED,
`Restart intent: stationary evaluate`, `OMP_NUM_THREADS=1`, exit 0, written
header `certified=F cert_reason=no_stationary_claim mode=init`).

## 6b. How far the pre-existing defect reached in the kept seeds (advisor survey, 2026-09-16)

MEASURED over every `Hydro_ioniz_IC.txt` still in the tree that carries a
`# mapped-from-coupling:` line (`LHS1140b/models/*/*/output/`, `models/.stopped/*/output/`,
`models/.L*/*/output/`, `models_20260915_db87/`, `models_20260914_preL21/`):
the six fields the old selector let the provenance line override
(`sec_ion`, `sec_ion_step`, `recon`, `certified`, `mode`, `t_phys`) were
compared between the seed's own `# coupling:` line and the provenance line.
The only field that ever differed is `certified`: eight seeds under
`models/.L4h/`, `models/.L14/` and `models/.L18/probe_seedmap/` state
`certified=F` on their own line and `certified=T` on the provenance line
(their source states were certified), so the old reader loaded those seeds
as certified states. `sec_ion_step`, `recon`, `mode` and `t_phys` agreed on
every kept seed, so no kept seed was restarted with a wrong armed step,
reconstruction label or run state. The seeds of the campaign passes
themselves were overwritten by the solved states and could not be surveyed;
their provenance lines were written by the same mapper from the same kind of
source, so the same statement is expected to hold for them, and it is
expected, not measured.

## 6c. The route is metadata of a state, not another equation set (2026-09-17, item L22 step 3 increment I3)

The contract of section 2 says that what a file states about a state is
metadata OF that state. The `carrier_newton` token of the `# options` line was
read as something else: as a statement of WHICH EQUATIONS the state solves. It
was listed among the tokens `Restart option change:` may never name, beside
`metals`, `mol`, `oxychem`, `carrier` and `iontrans`, so an alternation state
reloaded under `Coupled carrier solve: True` was refused with "The state in the
file is a state of another equation set" and the coupled block could not be
entered on a reload at all (MEASURED and reported by increment I1).

That reading is wrong on the code. The five other tokens decide how many
unknowns a state has and which columns its files carry: `metals` adds the metal
ionization stages, `mol` the four molecular carriers, `oxychem` the three
oxygen carriers, `carrier` gives each carrier a continuity equation and a
column, `iontrans` makes the hydrogen ionization state a transported row.
`carrier_newton` adds no row and no column. It says whether the transported
balances are unknowns of the Newton vector, solved together with the wind as
one block, or are relaxed at a held wind in alternation with it. The balances
whose residual must vanish are the same balances, the certification evaluates
the same rows against the same tolerances, and a state that is stationary is
stationary under either. It is a route, which is metadata, and the contract
makes it admissible.

So `opt_is_route` is added beside `opt_changes_layout` in `load_IC.f90`, true
for `carrier_newton` alone. A difference in a route token is admissible with
nothing named; the load prints `the restart changes the ROUTE and not the
equations` with the two values; and the change is written into the state the
run produces as a `# route_change` line of the same form as `# option_change`,
inherited by the rungs that follow, so a ladder states which route reached each
of its rungs.

The second half is what the state WRITES. The token was formed from
`carrier_in_newton` alone, which `On stall` leaves at what the input said, so a
state produced after a handover recorded `carrier_newton=F` although the block
produced it (found by I2 in a file it could not edit). The token is now
`carrier_in_newton .or. carrier_rows_entered_newton`, the second being set at
the handover in `EXHALE_main.f90`, so the written token states the route that
produced the state.

Five rows in `src/tests/grid_and_gates/restart_option_change.sh` on the
`backup/regression/carrier_model_a_newton` fixture, which is molecular with the
H2 carrier transported and `Coupled carrier solve: False`, hence an alternation
state. MEASURED: four RED against a control binary built from the same tree
snapshot with these two files at their entry text, twenty of twenty GREEN
after. The fifth row, an alternation state reloaded under `On stall`, passes on
both binaries and is stated as such: `On stall` leaves the token at `F` until a
handover fires, so before one there is no difference to admit.

## 7. Files changed

| file | what |
|---|---|
| `src/modules/files_IO/load_IC.f90` | `ic_cert_reason` and `cert_reason_len` beside `ic_certified`; the `file_certification_claim` record; `next_coupling_field` factored out of `parse_coupling_header`; `parse_certification_claim`, `refuse_coupling_header` and `adopt_certification_claim` added; the `certified` case moved out of `parse_coupling_header`; both halves of the state read for the pair; the coupling line selected by its own label (section 6a) |
| `src/EXHALE_main.f90` | one line in the `EXHALE_DUMP_IC` branch: the imported pair put into the state metadata for a reload that writes without measuring |
| `src/tests/grid_and_gates/restart_intent_and_metadata.sh` | sixteen rows for the pair and for the mapped seed's provenance line, and `run_it` records each run's exit status beside the run |
| `src/modules/files_IO/load_IC.f90` (2026-09-17, section 6c) | `opt_is_route` beside `opt_changes_layout`, `carrier_newton` moved out of the layout tokens, the route branch and its report in `compare_options_field`, `carrier_rows_entered_newton` in the written token, `route_change` lines collected and inherited |
| `src/EXHALE_main.f90` (2026-09-17, section 6c) | `carrier_rows_entered_newton` set at the handover of `Coupled carrier solve: On stall` |
| `src/tests/grid_and_gates/restart_option_change.sh` (2026-09-17, section 6c) | five rows for the route token on the `carrier_model_a_newton` fixture |

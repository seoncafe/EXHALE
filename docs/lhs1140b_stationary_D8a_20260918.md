# D8, first increment: what a case directory publishes, and when a second solve is taken

2026-09-18. Item D8 of `docs/PLAN_20260918_rev2.md`, steps 1, 2 and 7 only
(preservation without a new format, the continuation policy, the evaluate pass
in its own directory, the campaign's aggregate status). Steps 3 to 6 of that
item, the generations with their manifests, the atomic pointer that publishes
one of them, the focused checkpoint writer and the warm-restart contract, are a
design to be agreed and are NOT built here.

No source under `src/` was changed and no file format was changed. This
increment is the runner's policy and its directory discipline.

Files: `LHS1140b/models/run_case.sh`, `LHS1140b/models/run_campaign.sh`,
`LHS1140b/models/tests/run_case_policy.sh` (new),
`LHS1140b/MODELS.md` section 6. `LHS1140b/models/` is outside the git remote;
only `MODELS.md` is tracked.

---

## 1. What a case directory holds

A case directory publishes ONE solve, in `output/`, and the policy names which
one: **the latest completed solve.** The runner compares no two solves and
makes no claim that the published one is the best of them; a same-state measure
that would let it rank them is step 4 of D8 and does not exist yet.

| what | where | written by |
|---|---|---|
| the published solve | `output/` | the last solve the runner ran |
| how it ended, one line | `ENDING` | `run_case.sh` |
| a solve the runner superseded | `solve_<n>/` holding `output/`, `run.log`, `EXHALE_setup.out`, `EXHALE_resolved.out`, `ENDING` | `keep_solve` |
| a seed attempt | `attempt_<n>/` holding the whole attempt, its own `solve_*` included, and `ENDING` | the seed walk |
| the reason a case did not solve | `not_solved.md`, first line `# not solved: <one sentence>` | `write_not_solved` |
| every generation, listed | `REPRODUCE.md`, under "Why this run was made the way it was" | `write_reproduce.py` through the runner's `--note` |
| the evaluate pass's own inputs and resolved configuration | `eval/` | the evaluate pass |

`n` in `solve_<n>` counts up within the case directory in the order the solves
were kept, so `solve_1/` is the earliest kept and the published state is
younger than all of them. `next_solve_dir` takes the first free index, so
nothing already there is written again.

`ENDING` replaces the `OUTCOME` file the seed walk used to write, so that one
kind of record has one name. `not_solved.md` is now the runner's: a file
already there is kept under `not_solved_<YYYYmmddHHMMSS>.md`, both when a new
run does not solve and when a run solves a case that did not solve before, so
the `reason` column of MODELS.md section 7 never states a reason for a state
that no longer has one.

**What is still not preserved.** The state a PREVIOUS RUN left in `output/` is
still overwritten by the seed of a new run. Keeping it means copying the state
pair at the start of every run of every case: MEASURED, `output/` is 1.7 MB in
the `carrier_model_a_newton` fixture and 2.5 MB in
`models/atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13`, so a re-run of the catalog
would add of order 250 MB each time. The lifetime of the generations and the
pointer that publishes one of them are steps 3 and 4 of D8, and this is left
for that decision rather than settled here.

## 2. The ending is classified before anything is decided on it

`classify_ending <log> <state directory>` READS the ending and names it. What
it reads: the solver's own verdict line
(`the stationary solve returned info = N`, or `(JFNK|PTC) done info=`), the
certification block of the state as written (the block headed
`final state, as written`, its three hydrodynamic row verdicts and the entries
listed after `NOT CERTIFIED:`), the last `outer pass N` line, and the two
halves of the state on disk.

| class | what it means |
|---|---|
| `solved` | `info = 0` |
| `no_verdict` | the log states no solver verdict at all |
| `state_missing` | refused, and the two halves are not both on disk |
| `nonfinite_state` | refused, and the state carries a value that is not finite |
| `composition_refusal` | refused with every hydrodynamic row within its own tolerance: the wind is stationary and what refuses is the composition (L14) |
| `hydrodynamic_refusal` | refused with at least one hydrodynamic row above its tolerance |
| `no_row_measures` | refused, and the block measured no hydrodynamic row |

Pass-budget exhaustion is not a class. It is a separate fact about the same
ending, carried in `BUDGET_EXHAUSTED` and stated in the reason sentence, so
that the plan's "composition-only budget exhaustion" is what it is: a
composition refusal that spent its budget.

The nonfinite test reads the state files past their comment lines, because the
word `provenance` in a header contains the letters `nan`.

## 3. The continuation is taken for one ending

`continuation_addresses` is true for `hydrodynamic_refusal` and nothing else.

Item L4e opened the continuation for that ending: the ramp of the
pseudo-transient hydrodynamic solve is tied to the line-search merit, which is
nearly flat on a state that is already close, so from such a state `dtau` never
leaves `dtau0 = 1.0` and the hydrodynamic rows stay stuck, while the same state
restarted at `DTAU0_CONTINUATION` reaches its root in a few Newton iterations
and a raw seed at that value fails.

It addresses no other ending:

- a **composition refusal** already has the stationary wind the solve found,
  and a higher pseudo-time start cannot move the element relaxation. MEASURED
  in item L34c on `molecular_scalar_gj1132_wellmixed/HeH0.083`: the first forty
  passes reached a gated carrier row of 4.13e-02, the continuation's forty took
  it back to 7.37e-02, and the state published was the continuation's;
- a **nonfinite or absent state** is not a state to restart from;
- with **no verdict** nothing says the solve ended rather than died;
- with **no row measures** nothing says the hydrodynamic rows are what refuses.

Until this change the runner launched the continuation whenever `info != 0` and
both output files existed, and classified the hydrodynamic rows only afterward,
where the classification decided the seed walk alone. The order is now
classification, then decision, and the class is printed to the log at both the
first solve and the continuation.

The seed walk keeps the endings it had: it is taken for
`hydrodynamic_refusal`, `no_row_measures`, `nonfinite_state` and
`state_missing`, and not for `composition_refusal` (L14) or `solved`.

## 4. The evaluate pass runs in `eval/`

`Restart intent: stationary evaluate` needs two keys of the input changed and
the binary reads `./input.inp` and writes `./output/` relative to the directory
it is started in. The runner used to change the case's own `input.inp` and put
it back when the binary returned, so the solve input was transformed for the
whole of the pass and permanently if the runner was killed inside it.

The pass now has a directory of its own. The case's `*.inp` files are copied
into `eval/`, the solved state is written into `eval/output/` as the `*_IC.txt`
pair `load_IC` reads, the two keys are changed in `eval/input.inp`, and the
binary is run with `eval/` as its working directory. What it produces is moved
into the case's `output/` afterwards, so the layout every reader of this tree
expects is unchanged: `Hydro_ioniz.txt` and `Ion_species.txt` written back with
the certification of this measurement, `*_adv.txt`, the heating and cooling
breakdowns, and `pp.log` in the case directory.

A value an input key names is written into `eval/input.inp` as an absolute
path, because `Spectrum file: ../../../sed/...` is relative to the case
directory and `eval/` is one level below it. The rule is mechanical: a value
with no space that is not already absolute and that names something reachable
from the case directory becomes that path, normalized.

Two consequences worth stating. `eval/` is left in place as the record of how
the products were made. And the `EXHALE_setup.out` and `EXHALE_resolved.out` of
the case directory now describe the SOLVE they sit beside, where before the
evaluate pass overwrote them with its own.

## 5. The campaign answers with a status file and an exit code

`run_campaign.sh`:

- the shell that runs each case has `set -o pipefail`, so the status recorded
  is `run_case.sh`'s own and not `tail`'s;
- each case writes its own status line into a directory of its own and the
  lines are gathered when the run ends, because the tree is on NFS where an
  append from two hosts at once is not atomic;
- `campaign_status.txt` holds one line for each case: its name, the exit status
  of `run_case.sh`, the ending class of the solve that case published (read
  from the case's own `ENDING`), and the last line the runner printed. A status
  file already there is kept under a dated name;
- the campaign exits 0 only when every case exited 0, and with the number of
  cases that did not otherwise, capped at 125.

## 6. Tests

`LHS1140b/models/tests/run_case_policy.sh`, on scratch copies of
`backup/regression/carrier_model_a_newton` with
`LHS1140b/models/EXHALE_3146d11b.x` (md5 `3146d11b4090306dcea75bb9718edd22`),
`EXHALE_OUTER_PASSES=2`, `EXHALE_JFNK_MAXIT=5`, one thread. The whole script
takes about 30 s. It sources `run_case.sh` with `RUN_CASE_POLICY_ONLY=1`, which
defines the policy and returns, so the rule is tested where it is written and
nowhere restated.

MEASURED 2026-09-18: 16 checks, 16 passed.

- **B**, real: the bounded carrier reload ends `info = 1` with all three
  hydrodynamic rows above tolerance, is classified `hydrodynamic_refusal`, the
  continuation is taken, and afterwards `solve_1/` holds the first solve's log
  and state unchanged (md5 compared) with its `ENDING`, while `output/` holds
  the continuation's.
- **A**, synthetic log, labeled: no short fixture ends on the composition only,
  so B's real log is edited (the three row verdicts set to `within`, their
  three refusal entries removed, the carrier rows left refusing). The
  classification is `composition_refusal`, the continuation is not taken,
  `output/` is byte-identical afterwards, no generation is made, and
  `not_solved.md` states the reason.
- **C**, real, with a control: the superseded sequence is run and the shell
  running it is stopped by signal once its window is open (its `input.inp`
  carries the evaluate key); it leaves the case's `input.inp` transformed, which
  is the RED. The sequence used now is run and stopped the same way and leaves
  the case's `input.inp` byte-identical, which is the GREEN. A static check adds
  that no statement of `run_case.sh` writes the case's own `input.inp` at all.
  Every process is stopped by PID after its working directory has been checked.
- **D**: the campaign is run with a stub case runner that fails one case of
  two, so the check is of the campaign's own aggregation and costs no solve; it
  exits nonzero and both cases have a status line with their exit status and
  ending class.

## 7. Noticed and not fixed

`src/utils/element_flux_closure.py` (lines 654 to 676) carries the same
continuation rule the runner had, and the same defect: it continues at the
raised pseudo-time start whenever `info != 0` and a state was written, with no
classification of the ending. It also moves the first log to
`run_dtau0_first.log`, the name this runner no longer uses. The file is not in
this increment's ownership.

# D9, steps 1 and 2: the budget of a solve by configuration, and seed compatibility before distance

2026-09-18. Item D9 of `docs/PLAN_20260918_rev2.md`, steps 1 and 2 only. Step
3, the reference refresh, is acceptance work that follows D6, D7a and any
boundary or closure correction, and is NOT done here. No source under
`src/modules/` was changed, no golden was refreshed and no catalog case
directory was written.

Files: `LHS1140b/models/run_campaign.sh`, `LHS1140b/models/pick_seed.py`,
`src/utils/element_flux_closure.py` (the continuation block and the two
helpers it calls), `LHS1140b/models/tests/campaign_budget.sh` and
`tests/seed_compatibility.sh` and `tests/closure_continuation.sh` (new),
`LHS1140b/MODELS.md` sections 3 and 6. `LHS1140b/models/` is outside the git
remote (`git ls-files LHS1140b/models/` is empty), so the runner, the seed
chooser and the tests exist only in this working copy; of the files touched,
`MODELS.md` and `src/utils/element_flux_closure.py` are the tracked ones.

---

## 1. The budget by configuration

A solve is given a PASS ceiling and a WALL-CLOCK ceiling. Both are EMPIRICAL
ALLOWANCES read off recorded solves of the same configuration class. **No rule
here derives a budget from a residual, a cell count or a front speed**: the
second review withdrew the one such rule that had been written ("about 3.3
cells a pass, at least 40 passes before any gated row can begin to fall"),
because the binding cell is the argmax of a normalized residual, it can jump
between cells and species, and its movement is neither a measured front speed
nor a lower bound on any other configuration.

What decides a case's class is which system it solves, read from what the case
states: `closure.json` for a closure rung, then `Ionization transport:` and
`Molecular chemistry:` in its input.

| class | passes | wall | where the numbers come from |
|---|---:|---:|---|
| `atomic_prescribed` | 40 | 30 m | READ, L34a: the seven low-XUV Roe cases were ACCEPTED at outer pass 1. READ, the L14 measurements quoted in `run_case.sh`: `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` certified at pass 13 and `x0.10/HeH2.13` was ACCEPTED at pass 30 exactly, the largest recorded atomic pass count, which is why the ceiling is 40 and not the binary's 20. MEASURED 2026-09-18 over the 86 `REPRODUCE.md` records of the atomic catalog: wall clock 31 s minimum, 40 s median, 19 m 34 s maximum (`atomic_photochem_gj1132_kzzprofile/HeH2.09`) |
| `atomic_closure_rung` | 20 per EXHALE solve | 45 m for the rung | MEASURED 2026-09-18 over the nine closure logs of `atomic_photochem_gj1132_kzzprofile`: every EXHALE solve of every rung ended at outer pass 1, the longest single solve took 107 s, a rung is 5 or 6 solves, and the longest whole rung took 1132 s (18 m 52 s). 20 is the binary's own default, which the closure driver leaves in place |
| `molecular_alternation` | 40 | 6 h | MEASURED 2026-09-18 from the `REPRODUCE.md` of the molecular catalog, and READ from L34b: the reference `molecular_scalar_gj1132_kzz1e9/HeH2.13` certified at pass 12 in 40 m 35 s and `HeH9.7` at pass 10 in 25 m 35 s; the five cases that did not solve ran 37 m 41 s, 233 m 16 s, and three of them 360 m 32 s, 360 m 41 s, 360 m 47 s and 361 m 30 s, which is the six hours at which they were stopped by hand. The ceiling states that stop instead of leaving it to the operator |
| `transported_ionization` | 90 | 6 h | READ, L36e section 4: the atomic fiducial with the key on was ACCEPTED and CERTIFIED at outer pass 56 inside a 90-pass allowance, having been NOT certified at 25. That is evidence for that configuration, seed, algorithm and tolerances and for no other. **No wall clock was recorded for that solve**, so the six hours are the operating limit carried over from the molecular runs and are not a measurement of this class; a case that reaches it says so in its status line |

**Overriding one case.** `models/budget_overrides.txt`:

```
<group>/<case>  passes=<N>  wall=<duration>  authorized: <sentence>
```

The `authorized:` sentence is not optional. An extension of a budget is a
decision about what the campaign will accept, so it is named where it is taken
and the campaign copies it into that case's status line. A line without it is
refused and the campaign stops before running anything.
`EXHALE_OUTER_PASSES` or `CAMPAIGN_WALL` in the environment replaces the table
for a whole run, and every status line then reads `source=environment`.

**A case that reaches a ceiling is UNCERTIFIED and INCOMPLETE, never a
solution.** At the PASS ceiling `run_case.sh` classifies the ending itself
(D8a) and writes the reason into `ENDING` and the residual profile's refusing
entries, the binding row and cell, the solver verdict, the pass count and the
seed into `not_solved.md`; the campaign adds `ceiling=passes` and the words
`UNCERTIFIED INCOMPLETE` to the status line. At the WALL ceiling the runner is
stopped from outside the binary and leaves no ending of its own, so the
campaign records `ceiling=wall`, the class `wall_ceiling_reached` and the same
words. Neither is relabeled a solution because a front is suspected in what
remains.

The status line:

```
<case> exit=<status of run_case.sh> <ending class> budget=<class>:<passes>p/<wall>
   source=<table|environment|override> spent=<s> ceiling=<none|passes|wall>
   [UNCERTIFIED INCOMPLETE] [authorized=<sentence>] | <last line the runner printed>
```

## 2. Seed compatibility, in every tier, before any distance

`pick_seed.py::Physics.matches` returned true for two inputs differing only in
the ionization-transport key (the second review's section 5.4 construction),
and tier 0, the case's own directory, was taken without any comparison at all.
Compatibility is now decided first, in every tier including tier 0, and only
what survives it is ranked by distance.

**What is compared, and the loader line each test anticipates.** Every line
number is of `src/modules/files_IO/load_IC.f90`.

| compared | the refusal it anticipates |
|---|---|
| the spectrum, and its normalization where the tier does not relax it | no loader line: the radiation field is not in the restart metadata block, and a state solved in another field is not a state of this problem |
| `K_zz`, where the tier does not relax it | likewise |
| the lower boundary (scalar, scalar with the C/N/O reservoirs, profile) | the `grid` field, compared as text (lines 182 to 186, 1222 to 1231): a profile case builds its grid on the profile's radius at the matching level, so its R0 is its own |
| the elemental reservoir INVENTORY, which elements the `# reservoir` line names | `compare_reservoir_field` and `refuse_reservoir_element`, lines 1772 to 1830: an element present in one reservoir and absent from the other is a different composition whatever the tolerance. The RATIOS are carried onto the case by `map_state_to_grid.py --reservoir` and are distance, not compatibility |
| `metals`, `mol`, `oxychem`, `carrier` | `opt_changes_layout`, lines 245 to 272: these four decide which species the state files carry, and a restart may never be told to allow them to differ. "A state whose columns are not this run's columns is not this run's state, so such a change is a cold start and not a restart" |
| `he_diff` | lines 650 to 685: under `he_diffusion` only the base cells are set to the reservoir composition and every cell above keeps the helium fraction it was written with, so a diffused column carries a He/H that varies with radius while a well-mixed run carries one number. Line 663 also refuses the one-factor rescale outright when the state carries HeH+ |
| `iontrans` | lines 262 to 278 and 1646 to 1690: the three ionization stages have a column in every state file, so a state written without the key IS an admissible starting point for a run with it, but the difference is a change of the equations and `load_IC` takes it only when the run names the token on a `Restart option change:` line |

**A model-option transition is printed as one.** A candidate whose only
difference is `iontrans` is offered when, and only when, the case names
`iontrans` on its own `Restart option change:` line. The seed line then reads

```
compat=model-option transition: iontrans F->T, permitted by "Restart option
change"; the seed solves another system and its certificate does not transfer
```

and the candidate ranks after every candidate of the case's own system,
whatever the tier and whatever the distance. Where the case names no such
permission the candidate is excluded and `--why` gives the reason.

**No certificate transfers.** A certified state is preferred as a seed because
it is a fixed point of the equations the new case will be solved with, not
because its verdict says anything about the new case; where the seed solves
another system the line says so in as many words.

**Where a candidate's configuration is read from.** In three layers, each
overwriting the one before it token by token: the directory's `input.inp`,
then `EXHALE_resolved.out` beside the state (`carrier_transport`,
`oxygen_chemistry`, `ionization_transport`, `he_diffusion`), then the restart
metadata block of `output/Hydro_ioniz.txt` (`# options`, `# reservoir`), which
is what the generation that wrote the state was configured with and is what
`load_IC` compares against. The case's own configuration is read from its
`input.inp` alone: its earlier state is a candidate, not a statement of what
the next run asks for. An archived state carries no metadata block and no
resolved configuration, so its tokens come from its own input and the seed
line says the configuration was read from the input.

A first draft compared only what both sides happened to state and skipped a
token neither the resolved configuration nor the state named. That silently
readmitted 262 archived states to a molecular case (MEASURED: candidates went
from 2 to 264), which is the failure `load_IC` refuses a whole file for ("the
comparison would silently skip part of the field", lines 1634 to 1641). The
layering above is the repair.

## 3. The closure driver takes the continuation by the same rule

`src/utils/element_flux_closure.py` launched the raised-pseudo-time
continuation whenever `info != 0` and a state was written, with no
classification of the ending, and renamed the first log `run_dtau0_first.log`.
A closure rung therefore continued on a composition-only refusal exactly as
the catalog runner did before D8a, and L34c measured what that costs on
`molecular_scalar_gj1132_wellmixed/HeH0.083`: the gated carrier row went from
4.13e-02 back to 7.37e-02 and the worse state was published.

**The rule is sourced, not reimplemented.** The driver runs `bash` on the
policy block of `LHS1140b/models/run_case.sh` with `RUN_CASE_POLICY_ONLY=1`,
which defines `classify_ending`, `continuation_addresses` and `keep_solve` and
returns, and reads back the ending class, the reason and the verdict. The
alternative offered was to reimplement the classification in Python over the
same log lines; sourcing was chosen because two readings of the same log are
two rules and will drift, and the one thing D8a established is that the ending
decides the continuation, so the ending cannot be read two ways.

The cost of that choice is a dependency in the direction
`src/utils/` -> `LHS1140b/models/`, which is outside the git remote. It is
made explicit: `RUN_CASE_POLICY` names the file where it is elsewhere, and a
checkout that does not carry it takes NO continuation and says so in the
closure log. That is the safe direction, an unclassified continuation being
the defect this replaces.

What the driver does now: classify, print the class and the reason into the
closure log, and where the class is `hydrodynamic_refusal` move the first
solve whole into the iteration's own `solve_<n>/` (its `output/`, `run.log`,
`EXHALE_setup.out`, `EXHALE_resolved.out` and a one-line `ENDING`), seed a
fresh `output/` from the kept state's pair and solve again at the raised
start. `run_dtau0_first.log` is gone.

## 4. The tests

Three scripts under `LHS1140b/models/tests/`, beside `run_case_policy.sh`.
None writes in a catalog case directory or in `backup/regression/`; the one
real solve any of them makes is a bounded reload of a scratch copy of
`backup/regression/carrier_model_a_newton` with
`LHS1140b/models/EXHALE_3146d11b.x` (md5
`3146d11b4090306dcea75bb9718edd22`), `EXHALE_OUTER_PASSES=2`,
`EXHALE_JFNK_MAXIT=5`, one thread, about 12 s.

**`tests/campaign_budget.sh`, MEASURED 2026-09-18: 12 checks, 12 passed.** No
solve at all: the case runner is a stub that records its environment and exits
with the status the check asks for, so what is under test is the campaign's
own arithmetic and what it hands the runner.

| check | what it establishes |
|---|---|
| E1 to E4 | the four configurations are classified from what they state; each class carries the pass and wall allowance of the table; the duration arithmetic; each allowance names its source in the file, the withdrawn rule is named as withdrawn, and the missing transported wall clock is stated |
| F1 to F3 | each case is run at its class's allowance (`EXHALE_OUTER_PASSES` 40, 40, 90 for the atomic, molecular and transported stubs), the status line names the class, the allowance, its source and what was spent, and the campaign exits 0 when every case did |
| G1 (RED) | an override naming no `authorized:` sentence is refused, nothing runs, no status file is written |
| G2 (GREEN) | an authorized override changes that case alone (150 passes, 12 h), `source=override`, and the sentence is in the status line |
| G3 | a budget stated in the environment replaces the table for every case and says `source=environment` |
| H1 | a case stopped at its wall ceiling: `exit=124`, class `wall_ceiling_reached`, `ceiling=wall UNCERTIFIED INCOMPLETE`, campaign exit 3 |
| H2 | a case whose own `ENDING` says it spent its pass budget: `ceiling=passes UNCERTIFIED INCOMPLETE`, campaign exit 3 |

**`tests/seed_compatibility.sh`, MEASURED 2026-09-18: 10 checks, 10 passed.**
Checks J to M run on a synthetic tree whose state files carry the restart
metadata block and one data row, the block being what decides compatibility;
check N runs the real tool on six real catalog cases.

| check | what it establishes |
|---|---|
| J1, J1b | two inputs differing ONLY in the ionization key: the local-equilibrium state is NOT offered (exit 3, start cold) and the reason names the key and the missing permission. RED, MEASURED against a copy of the tool as it stood before this item: it offered that state as `tier1`, `dlog10=0.0000` |
| J2 | with `Restart option change: iontrans` named, it IS offered, as a transition, and the line says the certificate does not transfer |
| J3 | a seed of the case's own system 0.6584 decades away in log He/H outranks a transition candidate at distance zero |
| K1 | a metal-free scalar state is refused for a profile case, naming the boundary |
| K2 | a state whose reservoir names C/H, N/H and O/H is refused for a case whose reservoir is He/H alone, naming both inventories |
| L1 | every diffused state of the tree is refused for a well-mixed case, naming `he_diff` |
| M1 | tier 0 is checked like every other: the case's own certified state, written without the ionization key while its input now asks for it, is refused and the seed comes from another tier |
| M2 | a molecular state is refused for an atomic case, naming the species the files carry |
| N1 | the six real cases `atomic_scalar_gj1132_kzz1e9/HeH2.13`, `atomic_scalar_gj1132_wellmixed/HeH0.42`, `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7`, `atomic_photochem_gj1132_kzzprofile/HeH9.7`, `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` and `molecular_scalar_gj1132_kzz1e9/HeH2.13` choose the seeds they chose before, with the same candidate counts (383, 267, 105, 1622, 7, 2) |

**`tests/closure_continuation.sh`, MEASURED 2026-09-18: 10 checks, 10
passed.** One real bounded solve gives a real hydrodynamic refusal; the
composition-only log is that log edited (the three row verdicts set to
`within`, their three refusal entries removed) and is labeled SYNTHETIC
wherever it is used, no short fixture ending on the composition alone.

| check | what it establishes |
|---|---|
| P0 | the driver resolves the policy block of the catalog runner |
| P1, P2 | the real refused solve is classified `hydrodynamic_refusal` (`info=1`, the mass row 4.896E-03 against 3.0E-12 at cell 196) and the continuation addresses it |
| P3 | the first solve is kept whole in `solve_1/` with its `ENDING`, its log and state byte-identical to what it wrote (md5 compared), `output/` is fresh and empty, and no `run_dtau0_first.log` is made |
| Q1, Q2 | the SYNTHETIC composition-only log is classified `composition_refusal` and the continuation does NOT address it |
| Q3 | nothing is kept and the state pair beside that log is byte-identical afterwards |
| R1 to R3 | static checks of the driver text: the continuation is taken under the classification, the first solve is kept, `run_dtau0_first` is gone, and the rule is sourced through `RUN_CASE_POLICY_ONLY` and not restated |

`tests/run_case_policy.sh` of item D8a was re-run after the campaign change:
MEASURED 2026-09-18, 16 checks, 16 passed. Its two campaign checks read the
status line by the prefix `<case> exit=<n> <class> `, which the new fields
follow rather than displace.

## 5. What this does not do

- **Step 3, the reference refresh, is not done.** The four key-on regression
  cases stale against the 2026-09-17 goldens (`hp_front` 1.6e-1 in the helium
  columns, `hp_zero_seed` 3.9e-3, `hp_trace_seed` 3.3e-3,
  `carrier_model_a_newton` 1.2e-4) stay observations. No golden was touched.
- **No declared matrix of seeds, resolutions and regimes was run.** The
  budgets above are read off the solves the catalog and the L34 and L36e items
  already recorded. Establishing which configurations need larger budgets from
  a matrix of runs is the part of step 1 that costs solver time and is not
  done here.
- **No bounded repeat of the L36e fiducial** at the 90-pass allowance was run,
  so the reproducibility of its pass-56 acceptance is not measured here.
- **The transported-ionization wall ceiling is not a measurement of its own
  class**, and says so wherever it appears.

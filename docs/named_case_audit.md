# Audit of the named cases in `backup/regression/`

**Sixteen of the case directories named below were renamed on 2026-09-16.**
This document uses the NEW names throughout; section 6 is the old-to-new
mapping, so a citation of an old name elsewhere still resolves.

Date: 2026-09-03.  Read-only: nothing under
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00` was modified, and every
row below is read from the case directories as they stand.

Prompted by the P55 finding that `jfnk_hd189` carries no `Solver:` key and so
never runs the JFNK solver its name promises.  **That is not one case.  It is a
property of nineteen of them**, and it has a second form nobody has looked for:
two cases carry the key and it has since changed meaning underneath them (1.8).
The two things worth separating are:

* the case whose INPUT does not do what its NAME says (`Solver:` absent, or
  present with a different meaning);
* the case whose STORED OUTPUT was not produced by the input beside it
  (an environment variable the harness does not set, an input edited after the
  run, or a binary two months older than the code).

**Everything here is a proposal.  No case was renamed, edited, quarantined or
deleted, and no decision is implied.**

## 0. What is in the directory

46 case directories (a directory with an `input.inp`), plus
`_quarantined/molecular_base_no_chemistry`.  Ten of them are `DEFAULT_CASES`; the other 36 are named
cases that only run when passed as arguments.

Directories that are not cases: `golden/` and the ten dated `golden_*/`
snapshots, `parse_golden/` and `parse_golden_blockE_20260903/` (stored parser
outputs, not runs), `_quarantined/`, `valve_sens/` (two sub-runs `eps4`,
`eps5`, no `input.inp` of its own) and `windae_oracle/` (three planet
sub-directories plus `run_oracle.py` and `relax_gate.py`).  Those last two are
*campaign* directories with their own structure; `run_check.sh` cannot take
either as an argument, and neither is mentioned in its header.

## 1. The findings, ranked

### 1.1 No case in `DEFAULT_CASES` ever runs the steady solver

Measured from the stored logs:

| case | `Solver:` | what actually stopped the run |
|---|---|---|
| `wasp_full` | **absent** | flux gate, `final: count=11289 du=9.999E-04` |
| `wasp_he23off` | **absent** | flux gate |
| the eight `mol_*` / `lower_profile` cases | `Newton` | `reached EXHALE_MAXSTEPS = 12000` |

The only "JFNK" string in any of the eight molecular logs is `input_read`
echoing `Solver: Newton, JFNK hand-off at du < 1.00E-02` -- the arming message.
The hand-off never fires, because every one of them is a fixed-step snapshot
that stops at `maxsteps` long before `du` reaches the switch.  `grep -c '(PTC)'`
and `grep -c 'Newton finish'` are zero in all ten.

So `solve_steady_jfnk` and `solve_steady_ptc` -- the routines sections 126, 133,
138, 139, 142 and 144 rewrote -- are **not guarded by the regression matrix at
all**.  The two cases that do reach a JFNK finish, `wasp_full_newton` and
`wasp_he23off_newton` (`converged: JFNK steady solution (||R|| < res...)`), are
outside `DEFAULT_CASES` and have no `golden/` entry.

*Proposal (decision: user).*  Either add `wasp_full_newton` and
`wasp_he23off_newton` to `DEFAULT_CASES` and snapshot goldens for them, or state
in the `run_check.sh` header that the matrix guards the marching operator only
and that the steady path is covered elsewhere.  Adding them costs one JFNK solve
per `make check`.

### 1.2 Nineteen cases carry no `Solver:` key

`crit_cold`, `jfnk_cold`, `jfnk_hd189`, `jfnk_hd189_tight`, `ptc_warm`,
`resid_cold_13p707`, `resid_golden`, `resid_on`, `resid_warm_13p716`,
`roundtrip`, `tpm_hd189`, `tpm_wasp`, `wasp_cno_aiolos`, `wasp_cno_chianti`,
`wasp_full`, `wasp_he23off`, `wasp_hybrid_finish`, `wasp_localdt`,
`wasp_localdt_cont`.

For the twelve whose names say `jfnk_`, `ptc_`, `resid_`, `crit_` or `tpm_`
this is a name/behavior mismatch, because those routes are selected by
ENVIRONMENT VARIABLES that `run_check.sh` does not set: their own logs say so.

| case | log says | selected by |
|---|---|---|
| `jfnk_cold` | `(ATES_main) PTC dtau0 = 1.01E-04`, `(JFNK) start ...` | `EXHALE_PTC=1 EXHALE_PTC_JFNK=1` |
| `jfnk_hd189`, `jfnk_hd189_tight` | `(EXHALE_main) PTC dtau0 = 1.00E+00`, `(JFNK) start ...` | same |
| `ptc_warm` | `(PTC) start ||R|| = 2.141E-02` | `EXHALE_PTC=1` |
| `resid_cold_13p707`, `resid_golden`, `resid_warm_13p716` | no log at all | `EXHALE_RESIDUAL=1` (they hold only an `_IC` pair) |
| `resid_on`, `tpm_wasp`, `crit_cold` | plain time integration | nothing; the name refers to a stored quantity |

Re-run under `run_check.sh` today, every one of them marches and stops on `du`.

*Proposal (decision: user).*  Three shapes, and they are not interchangeable:
(a) add the key the name promises (`Solver: Newton <du>` for the `jfnk_*`
cases), which changes what the case computes and therefore needs a fresh
golden; (b) record the environment in the case directory -- a one-line `env`
file that `run_check.sh` would source, in the same spirit as `maxsteps` -- which
keeps the stored result meaningful; (c) rename the case for the quantity it
actually stores.  (b) is the only one of the three that does not change any
number.

### 1.3 Seven cases carry one of the four keys Phase C is expected to refuse

Currently all four keys are ACCEPTED (`input_read.f90:659, 664, 706, 779`).

| case | key | name says |
|---|---|---|
| `jfnk_cold`, `jfnk_hd189`, `jfnk_hd189_tight`, `tpm_hd189` | `Valve eps: 1.0e-4` | nothing about the valve |
| `hydrostatic_base` | `Hydrostatic base: True` | case D, variant 4b |
| `base_ghost_T_continuous` | `Base ghost temperature: continuous` | case D, variant 4a |
| `base_velocity_massflux` | `Base velocity: massflux` | case D, variant 3 |

`base_velocity_massflux`, `base_ghost_T_continuous` and `hydrostatic_base` are the *point* of those three cases: each
exists to exercise one base-boundary variant, and the key IS the variant.  When
Phase C replaces the lower boundary with a characteristic face solve, those
three cases stop being runnable in their present form.  The four `Valve eps`
cases are different: there the key is incidental to the name, and three of them
already adopt the valve implicitly (the JFNK hand-off sets `valve_eps = 1e-4`
when the input leaves it off), so dropping the line would leave them equivalent.

*Proposal (decision: user).*  For `base_velocity_massflux` / `base_ghost_T_continuous` / `hydrostatic_base`: quarantine
them together with a `NOTE.txt` when Phase C lands, the way
`_quarantined/molecular_base_no_chemistry` is handled, and record in the note which base variant
each one measured -- they are the only record of that comparison.  For the four
`Valve eps` cases: drop the line, or keep it and expect a startup refusal.
Neither should be decided before Phase C fixes what the replacement key is.

### 1.4 Fourteen cases restart, and not one of their `_IC` files carries a coupling header

`Load IC? True`: `heh_1_lw_newton`, `heh_1_newton`,
`heh_1_newton_bigstack`, `jfnk_hd189`, `jfnk_hd189_tight`, `ptc_warm`,
`resid_cold_13p707`, `resid_golden`, `resid_warm_13p716`, `roundtrip`,
`tpm_hd189`, `tpm_wasp`, `wasp_hybrid_finish`, `wasp_localdt_cont`.

Every one of their `output/*_IC.txt` files was written before the
`# coupling:` line existed, so by the rule of section 144.2 each of them
restarts with **secondary ionization STAGED**, whatever the run that produced
the file was doing.  Two of them are older still and carry no `# columns`
header either -- `jfnk_hd189` (`Hydro_ioniz_IC.txt` 92232 B) and `resid_golden`
-- so they take the legacy headerless read and their **metal ionization state is
rebuilt from the abundance rather than restored**.

`roundtrip` is the extreme case and is already item 12.4: its
`Hydro_ioniz_IC.txt` is 0 bytes and there is no `Ion_species_IC.txt`, so it
aborts in `load_IC` before anything runs.

*Proposal (decision: user).*  Regenerate the `_IC` pair of any case whose
stored result is still cited, from a state the current binary wrote; for the
rest, add the fact to the `run_check.sh` header so that "this case restarts
staged and without its metals" is documented rather than discovered.  A schema
version in the restart file -- the review's own recommendation for item 11 --
would turn all of this into a startup message.

### 1.5 The stored `Resid tol` no longer matches the stored log

| case | `input.inp` | its own log |
|---|---|---|
| `jfnk_cold` | `1.0e-5` | `tol = 1.00E-03` |
| `jfnk_hd189` | `1.0e-5` | `1.00E-03` |
| `jfnk_hd189_tight` | `1.0e-6` | `1.00E-04` |
| `resid_on` | `1.0e-4` | `1.00E-02` |
| `tpm_hd189`, `tpm_wasp` | `1.0e-5` | `1.00E-03` |

A clean factor 100 in all six, which is the residual renormalization of section
133 (the rows stopped being divided by `|u|`, and the default moved from 1e-3 to
1e-5).  So the inputs were bulk-edited to the new scale and the outputs beside
them were not regenerated: **the stored result of each of these six was produced
by a different input from the one in the directory.**

*Proposal (decision: user).*  Re-run the six, or mark their stored outputs as
pre-section-133 in the header.  Leaving an unmarked stale value in place is the
one option the standing rule excludes.

### 1.6 Log provenance: twelve cases were last run before the ATES -> EXHALE rename

Their logs still say `(ATES_main)`: `crit_cold`, `jfnk_cold`, `newton_rsw01`,
`newton_rsw05`, `resid_on`, `tpm_hd189`, `tpm_wasp`, `wasp_cno_aiolos`,
`wasp_cno_chianti`, `wasp_hybrid_finish`, `wasp_localdt`, `wasp_localdt_cont`.
The rename is 2026-07-06, so these results are at least two months and several
physics defaults old (CHIANTI C/N/O cooling, `eos_metals`, the Penning rate, the
He 2^3S set of section 116, the residual scale of 133, the row measure of 143).
Their file mtimes all read `2026-08-09 07:52`, which is the workspace-restore
date and not a run date -- the log banner is the only usable dating.

Only four case directories carry an `EXHALE.x` at all, so for the rest even the
binary cannot be identified:

| binary | md5 (8) | date | cases |
|---|---|---|---|
| current | `b85ff921` | 2026-09-03 21:05 | the ten `DEFAULT_CASES` |
| block-D era | `b9130a61` | 2026-09-02 05:01 | `wasp_full_newton`, `wasp_he23off_newton` |
| pre-rebuild | `49ee0692` | 2026-08-23 08:06 | `roundtrip` |

*Proposal (decision: user).*  The `# provenance:` header another worker added to
the output files on 2026-09-03 solves this going forward.  For the twelve
ATES-era cases the choice is to re-run them or to say in the header that their
stored numbers are historical.

### 1.7 Four cases have no runnable stored result

`heh_1_x2matched` has an empty `output/` and no `run.log`.
`resid_cold_13p707`, `resid_golden` and `resid_warm_13p716` hold only their
`_IC` pair -- no `run.log`, no `Hydro_ioniz.txt`.  The three `resid_*` are
inputs to an `EXHALE_RESIDUAL=1` measurement whose result was never stored
beside them; `heh_1_x2matched` has neither input state nor result.

*Proposal (decision: user).*  `heh_1_x2matched` is the natural candidate for
retirement (nothing distinguishes it from `heh_0p3`/`heh_3`/`heh_10`/`heh_30` except a name).  The three
`resid_*` are usable as they stand IF the environment they need is recorded
(1.2 proposal (b)).

### 1.8 `Solver: <x>` changed meaning, and two cases still carry the old one

`newton_rsw01` and `newton_rsw05` are the only two cases whose stored log reads

```
 (input_read) Solver: Newton, warm-up until ||R|| < 1.00E-01
```

Every other case that carries the key logs the current form,

```
 (input_read) Solver: Newton, JFNK hand-off at du < 1.00E-02
```

The second field of `Solver: Newton <x>` used to be a RESIDUAL threshold and is
now `newton_du_switch`, the flux-metric `du` at which the JFNK hand-off fires
(`input_read.f90:647-658`).  The number in those two `input.inp` files did not
change; what it instructs the code to do did.  `rsw` in the directory name is
"residual switch", so the NAME is accurate to the run that produced the stored
result and inaccurate to what the same input does today -- the sharpest form of
the mismatch this audit is about, because nothing in the directory looks wrong.

*Proposal (decision: user).*  Re-run both under the current meaning and rename
them for it (`newton_du0p1`, `newton_du0p5`), or keep them as the record of the
residual-switch experiment with a `NOTE.txt` saying that their `input.inp` no
longer reproduces them.  The second costs nothing and loses nothing; the first
is only worth doing if the `du`-switch sensitivity is a question anyone still
has.

### 1.9 Most `arm*` cases already carry a `NOTE.txt`; most others do not

20 directories carry a `NOTE.txt` or `README.md`: every `arm*` case, plus
`lower_profile`, `mol_diffusion`, `_quarantined/molecular_base_no_chemistry`, and the three that
record a deleted byte-identical duplicate (`resid_on` <- `ptc_wasp`,
`wasp_full_newton` <- `solver_newton_cold`, `wasp_hybrid_finish` <-
`crit_warm`).  Those notes are good and they already answer several of the
questions above -- `heh_1_newton/NOTE.txt` states both configuration changes
it makes and why, and `heh_1_newton_bigstack/NOTE.txt` states the environment
it needs and that it has no result to read.

The cases that carry no note are exactly the ones this audit had the most
trouble with: `crit_cold`, `jfnk_*`, `ptc_warm`, `resid_cold_13p707`,
`resid_golden`, `resid_warm_13p716`, `tpm_*`, `wasp_cno_*`, `wasp_localdt*`,
`roundtrip`, `heh_1_x2matched`.

*Proposal (decision: user).*  A one-paragraph `NOTE.txt` in each of those --
what it measures, what environment it needs, and whether its stored result is
current -- would remove most of sections 1.2, 1.5, 1.6 and 1.7 as open
questions without changing a single number.

## 2. The ten matrix cases

These are the cases `make check` runs.  All ten are internally consistent --
name, input and golden agree -- and all ten were re-run with the current binary
(`b85ff921`, 2026-09-03 21:05).  The one systematic remark is 1.1: none of them
reaches the steady solver.

| case | what the name says | what the input does | remarks |
|---|---|---|---|
| `wasp_full` | WASP-121b, full physics | He 2^3S on, metals on, `Solver:` **absent** | stops on the flux gate; the marching operator only |
| `wasp_he23off` | the He 2^3S-off branch | same, `Include He23S? False` | same |
| `mol_base_handoff` | molecular + `base.inp` handoff | hot Uranus, `Molecular chemistry: True`, `base.inp`, `maxsteps 12000` | fixed-step snapshot; `Solver: Newton` never fires |
| `mol_metals` | the molecular+metals coupled system | + `metals.inp` | same |
| `mol_lyman_werner` | the H2 Lyman-Werner branch | + `Stellar LW flux: 343.0` | same |
| `mol_diffusion` | binary H/He element diffusion | + `He_diffusion: True`, `He_Kzz 1e9`; planet `hotUranus_d` | same |
| `mol_ir_bands` | the molecular IR cooling channels | + `Base IR field`, `Molecular IR bands` | same |
| `mol_sec_ion` | the immediate secondary-ionization branch | + `Secondary_ionization: Immediate` | same; the only molecular case in which a photoelectron is partitioned |
| `mol_carrier` | the H2 carrier transport operator | + `Molecular carrier transport: True` | same; newest case -- present in `golden_blockG_20260903` only |
| `lower_profile` | the profile handoff instead of `base.inp` | HD 209458 b, `Lower atmosphere profile:`, metals-on with no `metals.inp` | same |

Golden coverage: all ten are in `golden/`.  Of the dated blocks, nine appear in
all ten `golden_*` snapshots; `mol_sec_ion` appears in six and `mol_carrier` in
one (`blockG_20260903`), which is simply their age.

## 3. The 36 named cases outside the matrix

The table has 37 rows: the 36 named cases plus `_quarantined/molecular_base_no_chemistry`.
"Name says" is read from the directory name; "input does" from `input.inp`;
"mismatch" lists only what the two disagree about, or a key at risk.  Blank
means the two agree.

| case | name says | input does | mismatch / at-risk | proposal |
|---|---|---|---|---|
| `heh_1_lw_40k` | case A, Lyman-Werner on | hotUranus_d, `Solver: Newton`, `maxsteps 40000` | | keep |
| `heh_1_40k` | case A, no LW | same without the LW flux | | keep |
| `heh_1_12k` | case D variant 2 | `Solver: Newton`, `maxsteps 12000` | | keep |
| `heh_1_lw_12k` | D2 + LW | same + LW flux | | keep |
| `heh_1_lw_newton` | D2 + LW, Newton finish | `Solver: Newton 100.`, `Load IC? True`, `maxsteps 2100` | `_IC` has no `# coupling:` (1.4) | regenerate the `_IC` pair |
| `heh_1_newton` | D2, Newton finish | `Solver: Newton 100.`, `Load IC? True`, no `maxsteps` | same | same |
| `heh_1_newton_bigstack` | the same, larger stack | `input.inp` byte-identical to `heh_1_newton`; the difference is the environment | none: `NOTE.txt` states it -- `OMP_NUM_THREADS=8`, `OMP_STACKSIZE=1G`, crash reproduction only, and its log ends at the segfault so there is deliberately no result | keep; it is the model for the `env` proposal of 1.2 |
| `base_velocity_massflux` | case D variant 3 | + **`Base velocity: massflux`** | Phase C is expected to refuse the key (1.3) | quarantine with a `NOTE.txt` when Phase C lands |
| `base_ghost_T_continuous` | case D variant 4a | + **`Base ghost temperature: continuous`** | same | same |
| `hydrostatic_base` | case D variant 4b | + **`Hydrostatic base: True`** | same | same |
| `heh_0p3` | the He/H = 0.3 rung | `He/H number ratio: 0.3` | | keep |
| `heh_3` | He/H = 3 | matches | | keep |
| `heh_10` | He/H = 10 | matches | | keep |
| `heh_30` | He/H = 30 | matches | | keep |
| `heh_1_x2matched` | He/H = 1, x2 matched | input present, **no output, no `run.log`** (1.7) | nothing distinguishes it from the `heh_0p3`/`heh_3`/`heh_10`/`heh_30` rungs | retire, or run it and say what it measures |
| `crit_cold` | critical point, cold start | WASP-121b, `Solver:` absent, plain marching, ATES-era log | the name refers to a stored quantity, not to a solver | rename for the quantity, or record the environment |
| `jfnk_cold` | JFNK from a cold start | `Solver:` **absent**; log shows `(JFNK)` under `EXHALE_PTC=1`; **`Valve eps`**; `Resid tol` 1.0e-5 against a logged 1.00E-03 | 1.2, 1.3, 1.5, 1.6 | record the environment (`env` file); re-run or mark the stored numbers |
| `jfnk_hd189` | JFNK on HD 189733 b | `Solver:` **absent**; `Load IC? True` with a **headerless** `_IC`; **`Valve eps`**; tol mismatch; `Reconstruction: PLM` | 1.2, 1.3, 1.4 (metals rebuilt from abundance), 1.5 | the P55 case; regenerate the `_IC` pair AND record the environment |
| `jfnk_hd189_tight` | the same at a tighter tol | same, `Resid tol 1.0e-6` against a logged 1.00E-04 | same | same |
| `newton_rsw01` | Newton, residual switch 0.1 | `Solver: Newton 0.1`; log reads `warm-up until \|\|R\|\| < 1.00E-01` | **the key changed meaning** (1.8): the same line now means "hand off at du < 0.1"; ATES-era log (1.6) | `NOTE.txt` saying the input no longer reproduces the result, or re-run and rename |
| `newton_rsw05` | the same at 0.5 | `Solver: Newton 0.5`, same old semantics | same | same |
| `ptc_warm` | PTC from a warm start | `Solver:` **absent**; log shows `(PTC) start`; `Load IC? True`, `_IC` without a coupling header | 1.2, 1.4 | record `EXHALE_PTC=1`; regenerate the `_IC` |
| `resid_cold_13p707` | residual at the cold `log Mdot = 13.707` state | `Load IC? True`, `_IC` pair only, **no result stored** | 1.2, 1.4, 1.7 | record `EXHALE_RESIDUAL=1` and store the measurement beside the input |
| `resid_golden` | the reference residual state | same, and its `_IC` is **headerless** | same, plus metals from abundance | same |
| `resid_warm_13p716` | residual at the warm 13.716 state | same as `resid_cold_13p707` | same | same |
| `resid_on` | residual-based convergence on | `Solver:` absent; `Resid tol 1.0e-4` against a logged 1.00E-02; ATES-era log | 1.2, 1.5, 1.6 | re-run, or mark the stored numbers pre-133 |
| `roundtrip` | the restart round trip | `Load IC? True`, `_IC` 0 bytes, no `Ion_species_IC.txt` -- **cannot start** | item 12.4 | rebuilt case ready in `$S/roundtrip_case` (section 149.5) |
| `tpm_hd189` | transit post-process, HD 189733 b | **`Do only PP: True`** -- no wind is solved; carries `Resid tol` and **`Valve eps`**, both inert on that route | 1.3, 1.5, 1.6 | drop the two inert keys; the name is accurate |
| `tpm_wasp` | transit post-process, WASP-121b | `Do only PP: True`; `Resid tol` inert; ATES-era log | 1.5, 1.6 | same |
| `wasp_cno_aiolos` | C/N/O cooling, AIOLOS form | WASP-121b, plain marching, ATES-era log | 1.6 | the pair is the record of the `cno_cool` comparison; re-run both or mark both |
| `wasp_cno_chianti` | C/N/O cooling, CHIANTI fits | same with the other setting | 1.6 | same |
| `wasp_full_newton` | `wasp_full` with the Newton finish | `Solver: Newton`; log `converged: JFNK steady solution` | **not in `DEFAULT_CASES`, no golden** (1.1) | promote to the matrix, or say in the header that it is the only JFNK guard |
| `wasp_he23off_newton` | the He 2^3S-off counterpart | same | same | same |
| `wasp_hybrid_finish` | the hybrid stop | `Load IC? True`, `_IC` without a coupling header; ATES-era log; stopped on `du plateau` | 1.4, 1.6 | regenerate the `_IC`; the name is accurate |
| `wasp_localdt` | local time stepping | `Time stepping:` key present, ATES-era log | 1.6 | keep; re-run or mark |
| `wasp_localdt_cont` | its continuation | `Load IC? True`, `_IC` without a coupling header; ATES-era log | 1.4, 1.6 | same |
| `_quarantined/molecular_base_no_chemistry` | case D variant 1 | already quarantined: refused at startup by the P35 check (`Molecular base: True` with `Molecular chemistry: False`) | | leave as is |

## 4. What this audit did not do

* Nothing was re-run.  Every statement about what a case *did* is read from the
  log it stores; every statement about what it *would* do now is read from
  `input.inp` and from the code that parses it.  Where the two disagree the row
  says so, and I did not attempt to decide which is right by running it.
* The four at-risk keys are read as "Phase C is expected to refuse them".  Today
  `input_read.f90` accepts all four (`:659`, `:664`, `:706`, `:779`); the
  refusal does not exist yet.
* `valve_sens/` and `windae_oracle/` were listed, not audited: neither is a
  `run_check.sh` case and each has its own structure.
* Whether the `_IC` files of section 1.4 still describe a physically current
  state was not measured -- only which headers they carry.

---

## 5. What was executed, 2026-09-04

The audit above proposes; this section records what was actually done, by whose
decision, and what was deliberately not touched. Nothing in it re-ran a case or
rewrote an input value.

**Finding 1.1 was treated as a defect, not a tidying job.** The default matrix
never ran the steady solver, so every rewrite of `solve_steady_jfnk` /
`solve_steady_ptc` and of the acceptance measures (sections 126 to 155) was
unguarded by a golden. `wasp_full_newton` is now the eleventh default case: it
is the one named case that both converges (`info = 0`) and reaches the solver,
so it guards that path and is a byte-identity test of it at the same time. The
reasoning is written into the `run_check.sh` header beside the case list.

**The four keys the characteristic base boundary retires** (section 152:
`Valve eps`, `Hydrostatic base`, `Base ghost temperature`,
`Base velocity: massflux`) were split by whether the key was incidental to the
case or was the case:

* `jfnk_cold`, `jfnk_hd189`, `jfnk_hd189_tight`, `tpm_hd189` carried
  `Valve eps: 1.0e-4` incidentally. The line is deleted and each `NOTE.txt`
  records that the deletion is equivalent either way, the JFNK hand-off having
  adopted that same value of its own accord.
* `base_velocity_massflux`, `base_ghost_T_continuous`, `hydrostatic_base` exist *to exercise* their key. They are
  moved to `_quarantined/` beside `molecular_base_no_chemistry`, with their inputs untouched and a
  `NOTE.txt` saying that their stored results are from a boundary condition the
  code no longer has, and that whether the case was testing the option or the
  physics the option stood in for is not inferable from the input. That is the
  campaign's intent to state, not the code's.

**Finding 1.8, the silent meaning change, was recorded and nothing was
rewritten.** `newton_rsw01` and `newton_rsw05` carry `Solver: Newton 0.1` and
`Solver: Newton 0.5`, whose second field was a residual warm-up threshold when
they were written and is now `newton_du_switch`. The number is unchanged, so
the case now means a `du` switch; the stored log describes the older meaning.
Both readings are legitimate cases, so each `NOTE.txt` states which is which
and the inputs stand.

**Finding 1.4, headerless restarts**: the 13 cases whose `_IC` predates the
restart schema of section 144 each gained a `NOTE.txt` paragraph saying that a
headerless restart starts with staged secondary ionization by design rather
than restoring the coupling state the writing run reached, so re-running will
not reproduce the stored log. Two of them (`jfnk_hd189`, `resid_golden`) also
predate the `# columns` header and are read by the legacy positional rule,
which reconstructs the metal ionization state from abundances rather than
restoring it; their notes say so.

**The `run_check.sh` header** now also lists what in this directory is not a
case -- `golden/` and the dated `golden_block*`, `parse_golden*/`,
`_quarantined/`, and the campaign directories `valve_sens/` and
`windae_oracle/` -- and the three byte-duplicates removed on 2026-09-03 with
the survivor each one duplicated.

**Not done, and left as the user's call:** findings 1.2 (19 cases whose
`jfnk_`/`ptc_`/`resid_` names came from environment variables the harness does
not set), 1.5 (`Resid tol` inputs rescaled by section 133 without regenerating
the outputs beside them), 1.6 (12 cases whose logs predate the ATES to EXHALE
rename), and 1.7 (four cases with no results at all). Each is a naming or
provenance question rather than a code defect, and renaming a case discards the
only record of what it once measured.

---

## 6. The renaming of 2026-09-16

Sixteen case directories carried the noun "arm", which the naming rule of
`~/.claude/CLAUDE.md` forbids, and none of the names said what the case is.
They were renamed on 2026-09-16 (PLAN_20260916_rev3 section 10, the user's
decision). **No number was regenerated and no golden was refreshed**: none of
the sixteen is in `DEFAULT_CASES` and none has an entry in `golden/` or in any
dated `golden_*/`, so the rename moved directories and citations only. The
`input.inp`, `base.inp`, `output/` and `run.log` of every case were left
untouched; only `NOTE.txt` was edited, to follow the sibling names and to
record the old name.

The twelve active directories are all the same planet and base state (the
Tier-2 gate hot Uranus, `hotUranus_diffusion`, molecular chemistry on,
`He_diffusion: True`, `p_base` 1 microbar, `T_base` 1140 K) and differ only in
the He/H number ratio, in whether the stellar Lyman-Werner field is on, and in
what the stored state is. The new names say exactly that: `heh_<ratio>` for the
ladder rung, `_lw_` for the Lyman-Werner field, then the stored state.

| old | new | what the case is |
|---|---|---|
| `armHeH_0p3` | `heh_0p3` | ladder rung He/H = 0.3, `q_H2_base` 0.610353658, 12000-step marching state |
| `armHeH_3` | `heh_3` | He/H = 3, `q_H2_base` 0.140486207, 12000-step state |
| `armHeH_10` | `heh_10` | He/H = 10, `q_H2_base` 0.046893592, 12000-step state |
| `armHeH_30` | `heh_30` | He/H = 30, `q_H2_base` 0.016151029, 12000-step state |
| `arm_heh1_x2matched` | `heh_1_x2matched` | the He/H = 1 rung, `q_H2_base` 0.326896922; input only, no stored output |
| `armD_D2` | `heh_1_12k` | the same input, 12000-step marching state |
| `armA_noLW` | `heh_1_40k` | the same input, 40000-step marching state |
| `armD_D2_LW` | `heh_1_lw_12k` | He/H = 1 with the stellar Lyman-Werner field, 12000-step state |
| `armA_LW` | `heh_1_lw_40k` | the same with the field, 40000-step state |
| `armD_D2_newton` | `heh_1_newton` | JFNK steady solve restarted from `heh_1_12k` (`Solver: Newton 100.`, PLM threshold 1e9) |
| `armD_D2_LW_newton` | `heh_1_lw_newton` | the same with the field, restarted from `heh_1_lw_12k`, 2100-step cap |
| `armD_D2_newton_bigstack` | `heh_1_newton_bigstack` | `heh_1_newton`'s input at `OMP_NUM_THREADS=8`, `OMP_STACKSIZE=1G`; segfault reproduction, no result |

The four quarantined directories are named for the key each one exercises, or
for the combination that refuses it:

| old | new | what the case is |
|---|---|---|
| `_quarantined/armD_D1` | `_quarantined/molecular_base_no_chemistry` | `Molecular base: True` with `Molecular chemistry: False`; refused at startup by the P35 check |
| `_quarantined/armD_D3` | `_quarantined/base_velocity_massflux` | exercises `Base velocity: massflux`, retired by section 152 |
| `_quarantined/armD_D4a` | `_quarantined/base_ghost_T_continuous` | exercises `Base ghost temperature: continuous`, retired by section 152 |
| `_quarantined/armD_D4b` | `_quarantined/hydrostatic_base` | exercises `Hydrostatic base: True`, retired by section 152 |

The brief of the item listed twelve directories plus `_quarantined/armD_D1`.
`armD_D3`, `armD_D4a` and `armD_D4b` moved to `_quarantined/` on 2026-09-04,
after section 3 of this document was written, and carry the same forbidden
noun, so they were renamed in the same change rather than left as the only
remaining occurrences.

Read from the directories while renaming them, and corrected in the `NOTE.txt`
it stands in rather than left stale:

* `heh_1_40k`'s note said the case was "run to Newton convergence (no
  EXHALE_MAXSTEPS)". Its stored `run.log` ends at the 40000-step cap of its own
  `maxsteps` file with `du = 3.8E+08` (READ), so the stored state is a marching
  state, not a converged one. The same holds for `heh_1_lw_40k` (`du =
  2.8E+08`, READ).
* The notes of `heh_1_lw_40k` and `heh_1_lw_12k` gave the stellar Lyman-Werner
  flux as 343.0 erg/cm2/s, which is what their stored logs echo (READ). Both
  `input.inp` files now carry 480.9, the value `mol_lyman_werner` was raised to
  (`Update_EXHALE_stage2.md`), so re-running will not reproduce the stored log.
  Each note now says so; the inputs were not touched.

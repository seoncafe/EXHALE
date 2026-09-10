# The restart contract (decision 15, option a): one intent, one metadata block

Advisor design, 2026-09-09, for items N10 and the metadata half of N9. The
user decided option (a) of decision 15 (PLAN_20260909_rev1 section 9).
Nothing here is implemented; N10 implements it after N11 (both write the
output headers).

## 1. What a restart is today (READ)

`Load IC? True` reads `output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`
(a renamed copy of a run's two state files), rebuilds the state
(`load_IC.f90`, the grid guard against `Grid cells:`, the He/H check
`heh_dev > 1e-6`, `equilibrate_loaded_composition`), and enters the marching
loop; the Newton follows only when the marching hands off (two CFL steps at
least, B5f's O(dt) kick). The files carry `# columns`, `# coupling`, and
`# provenance` (git hash, timestamp, `mode=init|phys`, `ck_input`), not the
reservoir abundance, the species schema version, the physical grid
parameters, the constant set or the option set.

## 2. The intent selector

One key in `input.inp`, read by `input_read.f90` beside `Load IC?`:

```
Restart intent: trajectory | relaxation | stationary
```

- `trajectory` (default when `Run mode: phys`): continue a physical
  trajectory; the loaded `t_phys` is the clock (today the header carries no
  `t_phys`, so a `phys` restart starts its clock at zero: the metadata block
  below adds it); every accepted step obeys the physical-mode contract.
- `relaxation` (default when `Run mode: init` or absent, i.e. today's
  behavior): continue relaxing toward stationarity with the marching, then
  the Newton hand-off; no clock.
- `stationary`: load, rebuild the derived quantities in the N16d order
  (composition, boundary, pressure map, temperature) without a step,
  evaluate the stationary residual and the certification on the state AS
  LOADED, and then either (i) return the certified state unchanged (written
  again with its certification, `EXHALE_STATIONARY_EVALUATE_ONLY=1` or a
  second word `evaluate`) or (ii) enter `solve_steady_jfnk` at once (the
  default of this intent). No CFL step, no `equilibrate_loaded_composition`
  before the evaluation (a loaded composition that is not the sweep's own
  fixed point is a FINDING the evaluation reports, not something to hide;
  the equilibration is offered as the second word `equilibrate` and is
  reported when taken).

Consistency: `Restart intent` without `Load IC? True` is an input error;
`trajectory` with `Run mode: init` is an input error; `stationary` with
`Solver:` other than `Newton` is an input error (there is nothing to enter).

## 3. The metadata block

Written by `write_output.f90` into every state file after `# provenance`,
one line per field, versioned:

```
# restart_schema 1
# reservoir He/H <value> [C/H <value> N/H ... for metals.inp elements]
# species_columns <n> <names as in # columns>
# grid N <cells> R0[cm] <value> r_min[Rp] <value> r_max[Rp] <value> mode <uniform|log|...>
# constants set <name or hash of parameters.f90 constants> RJ[cm] <value> kB[erg/K] <value>
# options <the resolved keys the setup report prints, one token each: He23S=T metals=T mol=F carriers=H2 ...>
# t_phys[s] <value or 0 with mode=init>
# source <git hash> <EXHALE.x md5>
```

Read by `load_IC.f90`: `restart_schema` absent means a legacy file, loaded
as today and marked `provenance unknown` in the setup report and in the
files the run writes; present means every field is compared with the run's
own (reservoir ratios within `heh_dev_tol`, grid exactly, constants exactly,
options exactly for the physics keys, the source only reported), a
mismatch is a refusal with the field named, except the source hash, which
is informational. The grid guard of today stays (it is the `grid` line's
exact comparison). A changed Jupiter radius is a changed physical grid and
is refused; the conservative remap is a separate authorized workflow.

## 4. What changes for the user

- New key `Restart intent` (optional; absent = today's behavior).
- State files gain eight comment lines; every reader skips `#`.
- A `phys` restart carries its clock.
- `stationary` gives the Stage C certification a way to re-evaluate a
  written state without a step, which is how "reproduced from an
  independently constructed nearby initial state" will be measured.

## 5. Tests (N10)

- intent parsing and the three consistency errors;
- `stationary evaluate` on a certified state (once one exists; until then
  on `wasp_full_newton`'s state, the three-unknown route): the written
  state and its certification identical to the writer's, no step taken;
- `stationary` on the atomic reload: the first JFNK line equals the one a
  marched hand-off would print for the same state (the residual is a state
  function, N5);
- legacy file (no schema) loaded and marked;
- schema mismatch in each field refused with the field named;
- round trip of the metadata (write, load, write: identical block).

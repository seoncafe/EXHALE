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

## 6. Which half of the pair states the density, and how closely a re-evaluation comes back (item L18, 2026-09-15)

The tests of section 5 ask for "the written state and its certification
identical to the writer's". What identical can mean was settled by measurement
in `docs/lhs1140b_stationary_L18_20260915.md`; this section is that result, so
the contract is read where the contract is.

### 6.1 The mass density is read, not rebuilt

The pair states the density twice: the `rho` column of `Hydro_ioniz`, and the
species densities of `Ion_species`, which weigh `sum_i m_i n_i` under the mass
policy of `calc_rho`. They are the same number only while the composition
closes its own mass, `sum_i f_i A_i = 1` -- the definition of `f_sp`, not a
tolerance. It does not close: the closure of a state written at the end of a
stationary solve drifts about 1e-14 per outer pass, monotonically, and reaches
5.1e-13 after nineteen passes and 6.2e-13 after twenty-four (a production state
of `LHS1140b/models/` closes to 1e-15).

**The conserved variable is the authority.** `rho` is what the hydrodynamics
advances and what the residual is a function of; the composition is an
eliminated variable the first sweep re-solves anyway. The loader therefore
reads the density from its own column and multiplies the loaded species, cell
by cell, by the single factor that puts them on it -- every element ratio,
every ionization split and every metal-to-hydrogen ratio survive a common
factor untouched, and the composition then closes the density it is a
composition of. The departure is reported on every restart.

**WHICH ONE IS THE AUTHORITY IS DECIDED BY INTENT, NOT BY THE SIZE OF THE
DISAGREEMENT** (corrected 2026-09-15 on the Codex review). The blocks of the
loader that change the composition on purpose -- the loaded H/He carried onto
the input's He/H, an element the file does not carry rebuilt at its abundance,
a reservoir the handoff moved, the oxygen carriers of a pre-oxygen-chemistry
file seeded -- each SAY SO, and where one of them acted the file's `rho`
describes a gas this run is not loading and the density follows the
composition. Where none did, the state is a restart of the same equations and
the conserved density is the authority. A magnitude cannot tell a deliberate
change from a damaged file, and reading it as one was the defect: a mismatched
pair would have been accepted as a deliberate change.

**A restart of the same equations is then allowed only rounding.** Below
`restart_density_rounding_tol = 1e-10` the two halves agree to what the writer
and the sweep's projection leave (item L19: 4e-16 on a state this code writes,
6e-13 on the longest solve before that projection existed); between that and
`restart_density_agreement_tol = 1e-8` the pair is loaded and the departure
reported; **above 1e-8 the pair is REFUSED** -- two halves of one state that
disagree by more than a part in 1e8, with no block of the loader having touched
either, is a damaged or mismatched pair, and loading it would be choosing
silently between two different gases.

**Since item L19 the same rule holds inside a run** (2026-09-15,
`docs/lhs1140b_stationary_L19_L20_20260915.md`): the composition an ionization
sweep returns is projected onto the density the sweep was given, by one factor
per cell, immediately after the sweep. So a state written, read back and
re-measured passes through ONE projection and not two different ones, and the
closure a file carries is the rounding of that projection (4e-16) rather than
a quantity that ratcheted with the number of sweeps (up to 6e-13 over
twenty-four outer passes). The certification report prints it,
`max_j |sum_i f_i A_i - 1|` with its worst cell, and gates nothing.

### 6.2 What is identical, and what is only close

- **Identical, bitwise**: the radius, velocity and pressure columns, which are
  carried in 17 significant figures and read as written, and the mass density,
  since 6.1. The state re-measured is the state certified.
- **Close to the file's own closure**: the composition, which the projection
  moves by the closure departure -- 1e-15 on a production state, up to 6e-13 on
  a long ladder solve. No tolerance in the certification inventory is within
  nine decades of that.
- **Not identical, and not reachable**: the row measures. The first equilibrium
  sweep of the re-entry moves the composition by 1e-14 to 1e-11 (one Picard
  step of a nonlocal coupling, from a state that is not exactly its own fixed
  point), and the flux assembly of a subsonic base is at its rounding floor:
  MEASURED on `.L14/x003_HeH2.13`, one unit in the last place of the density
  moves the cell-wise maximum of the energy row by about 11 per cent, because
  that row cancels its largest term by 6.2e+05.

Over the 120 certified states of `LHS1140b/models/`, re-evaluated through their
own files, the ratio of the re-measured cell-wise maximum to the in-run one is

| row | median | 10th | 90th |
|---|---|---|---|
| mass | 1.12 | 0.78 | 1.71 |
| energy | 1.22 | 0.89 | 1.93 |
| momentum | 9.2 | 1.80 | 134 |

the momentum row sitting four to six decades below its tolerance, where the
ratio of two draws from the rounding floor says nothing about the state.

**So a state certified at more than about half of its tolerance may be refused
on re-evaluation, and that refusal is a property of the cancellation and not a
defect of the file.** One of the 127 certified states re-evaluated for item L18
is refused on that ground. A certification is a statement about a state, made
once, by the run that produced it; a re-evaluation is a second measurement of
the same state with a different rounding, and the contract promises that the
state is the same, not that the measurement is.

## 7. What a certification is a certification OF (item R1 of the review of 2026-09-15)

**Certification is certification against the stationary operator.** A state
written by the stationary route carries `certified=T` because the rows of the
stationary residual `R_stat` were evaluated on it and were within their
tolerances. It is not a statement about the residual of the marched operator
`R_time` at the same state.

The two are not interchangeable while any part of the discretization differs
between them, and for a time one part did: the base boundary read the wind
window's mass flux on the stationary route and the first interior cell's
velocity while marching. MEASURED on the certified
`atomic_scalar_gj1132_kzz1e9/HeH2.13` state (the fiducial re-solved for item L21), the same state reads

| row | flux-keyed base face | local-face base |
|---|---|---|
| mass | 2.175E-09 | 1.306E+00 |
| momentum | 5.417E-13 | 1.702E-04 |
| energy | 2.126E-08 | 1.053E+00 |

-- nine decades apart, because the local face turns the base into an outflow
carrying -1.19 of the wind's own flux while the wind above it carries +1.000 of
it. `R_stat(U*) = 0` therefore said nothing about `R_time(U*)`.

That split is removed: there is now ONE base boundary, keyed on whether the
wind window carries a usable flux and not on which route is evaluating
(`base_boundary.f90`, `base_branch_on_wind_flux` and the weight `s`). The
label above is kept nonetheless, and stated here rather than left implicit,
because one boundary is a necessary and not a sufficient condition: the
marching and stationary paths still assemble their fluxes through different
code, and until a state is measured to be a fixed point of both, a
certification names the operator that made it.

**One boundary, and now one boundary that is a function of its argument
everywhere** (item L26, 2026-09-17). The label above says a certification names
the operator that made it, and that is only meaningful while the operator is a
function. The base boundary's entropy branch read the wind window through a
scale-free gate, and a scale-free gate has no limit at a zero window: the same
zero window was reached with two different face densities depending on the
direction it was approached from, MEASURED at 37.4 per cent of the face density
on the LHS 1140 b wind states. A residual with two values at a state is not a
residual, and no certification against it means anything there. The window's
say now carries an amplitude factor that vanishes quadratically with the window
flux, so the zero window is an ordinary point, and the range in the gate is
replaced by a standard deviation, so the residual is differentiable where the
extremal cell of the window changes as well. Neither change moves any state the
code meets: the two certified LHS 1140 b states re-evaluate to outputs that
differ from the control's in the provenance line alone. Memo:
`docs/lhs1140b_stationary_L26_20260916.md`, section "Repair".

### 7.1 Evaluating a state and producing a run's outputs from it (item L9, 2026-09-16)

`Restart intent: stationary evaluate` measures the loaded state and takes no
step. It also writes every product a run of the case is read for, and each
product says which of three states it describes.

| state | what is held | which product |
|---|---|---|
| **loaded** | everything: the conserved variables and the composition the two `_IC` files carry, copied at entry and not written to again | the file's own `certified=`/`cert_reason=` pair |
| **work** | the conserved density, momentum and total energy of the loaded state; the composition is the one one equilibrium sweep returns, and p and T follow that composition at that conserved energy through the caloric equation of state | `Hydro_ioniz.txt`, `Ion_species.txt`, the channel breakdowns, the residual breakdown, the certification, the mass-loss line |
| **advection-derived** | nothing of the above: it is another composition | `Hydro_ioniz_adv.txt`, `Ion_species_adv.txt` |

**The work state holds u, not p.** The pressure of a molecular cell is the
inverse of the caloric equation of state and therefore a function of the
composition at a given thermal energy. Holding p across a change of
composition assigns that composition a thermal energy that is not the file's,
which is a change of the conserved state and not a measurement of it. It is
the contract `pressure_and_temperature_at_fixed_conserved_state` states for
the molecular relaxation. Where no cell carries H2 the two rules are one map
and the state is unchanged bit for bit.

**"Refreshed" is not "closed".** One sweep is a single Picard step, so the
work state reports its own closure defect: the largest normalized reaction
residual the sweep accepted a cell state at, and the fractional distance
between the temperature each cell was solved at and the temperature the
returned composition has at the held energy.

**Two answers are kept apart.** Whether the file's own stationary claim
reproduces is a question about the FILE and is answered against the imported
pair (section 3 of `docs/lhs1140b_stationary_L23_20260916.md`); what verdict
the work state gets is a question about this run. A passing work state is
never a confirmation of a claim that did not reproduce.

**The derived product certifies nothing.** `Hydro_ioniz_adv.txt` carries no
`# coupling:` header. It states `# derived_from: certified=<T|F>
[cert_reason=...]`, the pair of the state its rows were built from, as
provenance, with the sentence that fixes it as provenance and not as a
certification of those rows.

This route replaced the `CFL 1.0e-12` step with `Do only PP: True` that the
LHS 1140 b runner used to take to reach the advection post-process. That step
was taken with the reconstruction `input.inp` names rather than the one the
solution was reached with, and it wrote a relaxation snapshot over the solved
state.

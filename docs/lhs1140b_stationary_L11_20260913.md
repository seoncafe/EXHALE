# L11: seeding one composition from the certified solution of another

Item L11 of `docs/PLAN_20260913_lhs_stationary.md` (row L4 of its table),
against `LHS1140b/MODELS.md` section 6. Written 2026-09-13.

## Verdict

**A case of the composition ladder is seeded from the certified solution of
its neighbour, and that is a better seed than the archive.** Three things
were needed and all three are in place:

1. `src/utils/map_state_to_grid.py --reservoir-HeH <value>` carries a solved
   state's helium onto another reservoir. Every helium column of every row
   is multiplied by the ratio of the two reservoirs, so the base rows -- the
   two inner ghosts and cell 1, where the element-diffusion operator holds
   its Dirichlet He/H -- arrive at the new composition exactly and the shape
   of the diffused He/H profile is left as it was. The `# reservoir`
   metadata line states the new value, which is what `load_IC` compares
   against the run's input and what it refused before.
2. `LHS1140b/models/pick_seed.py` prefers a CERTIFIED case of `models/` to
   any archived state, and says on its line which tier the seed came from.
3. `LHS1140b/models/run_case.sh` adds `--reservoir-HeH` by itself whenever
   the seed states a reservoir that is not the case's, and
   `write_reproduce.py` records the option and the factor in `REPRODUCE.md`.

MEASURED on the case the archive could not seed,
`atomic_scalar_gj1132_kzz1e9/HeH2.13`: from the archived state its
hydrodynamic rows never left `||R|| = 1.78` in five outer passes and it was
never certified; from the certified `HeH1.60` solution of the same group the
very first outer pass returns `info = 0` at `||R|| = 2.3E-08` and the case
is CERTIFIED at outer pass 4, `||R|| = 3.265E-08`, log10 Mdot = 7.88,
13m36s.

**One rung further is not certified, and not because of the seed.**
`HeH4.0` seeded from the certified `HeH2.13` also reaches `info = 0` at
`||R|| = 4.88E-08` in its first hydrodynamic solve, but the joint
alternation then refuses: the element row stalls at 3.2E-05 against its
1.0E-05 gate, the energy row of the OUTERMOST cell rises to 7.2E-03, the
under-relaxation is cut to 0.125 and outer pass 4 REFUSES. That is the
partitioned solver's behaviour at a composition step of 1.88 in helium
(`steady_newton.f90`, another worker's file); it is reported in section 4
and left open.

## 1. What the option does

`--reservoir-HeH <value>` is an INITIALIZATION choice and the file says so:
the `# coupling` line stays `mode=init t_phys=0 certified=F`, as it does for
every mapped state.

| quantity | what happens |
|---|---|
| He I, He II, He III, He 2^3S, HeH+ | multiplied by k = value / (source reservoir He/H) in every row, ghosts included |
| every hydrogen and metal column | untouched, bit for bit |
| pressure, velocity | untouched, bit for bit |
| mass density | the helium the rescaling adds, at the weights `calc_rho` uses: 3.9715259 m_H per free helium nucleus, 4.9715259 per HeH+ |
| temperature | `T = T_src n_src/(n_src + dn)`, `n_src = p/(k_B T_src)` read off the source state and `dn` the particles and electrons the added helium brings (He I 1, He II 2, He III 3, HeH+ 2). This is the relation the run itself carries, `T = p/((n_tot + n_e) k_B)` (`src/EXHALE_main.f90`, the ordering comment of the RK3 step), so the file stays one state |
| `# reservoir` line | its He/H replaced by the new value, every other element left alone; written in both output files |
| `# mapped:` line | records the two reservoirs, the factor, the He/H of the physical cells and of the base rows, and the temperature at the first physical cell |

### Why one factor over the whole column, and not the base rows alone

`load_IC.f90` has its own rule for a restart at another composition, and
with `He_diffusion` on it rescales the base rows only (`sH_l`, `sHe_l` set
to 1 for every cell above): the diffused He/H profile of a restart file is
the state being restarted and the operator produced it, so it is kept. That
is the right rule for RELOADING a state at its own composition. It is the
wrong one for CARRYING a state to another composition: it would leave the
whole column at the donor's helium and only the base at the new value.

Multiplying every row by the one factor does both things the rule is about.
The base rows land on the new reservoir exactly, because a converged
diffusion solution has them at its own reservoir exactly (MEASURED on
`atomic_scalar_gj1132_kzz1e9/HeH1.60`: rows 1, 2 and 3 of the file at He/H =
1.6000000000 and the profile falling to 0.38836 at the top). And
`He/H(r)/He/H(base)`, the shape the diffusion and the advection balance
produced, is unchanged by construction. So the rule adopted here is
`load_IC`'s rule with the column carried along: base rows at the input
value, the profile above it kept in shape.

### What it refuses

- a source with no `# reservoir` line: an archived state written before the
  metadata block. `load_IC` accepts such a file and rescales it by its own
  rule, so nothing is lost; forming a factor from a composition the file
  does not state would be inventing one.
- a value that is not finite and positive.
- a state whose base rows do not come out at the new He/H. That is the test
  that catches HeH+: it carries one nucleus of each element, so multiplying
  it moves the hydrogen count as well and no single factor sets both. It is
  the same reason `load_IC` refuses a molecular restart at a composition
  other than its own.

### The seed the runner wrote, measured

`HeH1.60/output` carried onto He/H = 2.13, k = 1.3312500000000000:

| row | r [R_p] | He/H source | He/H seed | T source [K] | T seed [K] | rho source [m_H/cm3] | rho seed |
|---|---|---|---|---|---|---|---|
| 1 (ghost) | 0.99981 | 1.6000 | 2.1300 | 395.01 | 328.13 | 5.298E+13 | 6.814E+13 |
| 2 (ghost) | 1.00000 | 1.6000 | 2.1300 | 389.76 | 323.76 | 5.192E+13 | 6.678E+13 |
| 3 (cell 1) | 1.00019 | 1.6000 | 2.1300 | 393.50 | 326.87 | 4.973E+13 | 6.397E+13 |
| 125 | 1.03792 | 1.5355 | 2.0441 | 1391.72 | 1159.21 | 8.383E+11 | 1.077E+12 |
| 301 | 1.84226 | 0.5449 | 0.7253 | 4454.64 | 3984.58 | 3.370E+08 | 4.134E+08 |
| 504 (ghost) | 30.00000 | 0.3884 | 0.5170 | 426.55 | 384.84 | 2.509E+04 | 3.013E+04 |

`He/H(seed)/He/H(source)` is 1.3312500000000000 to 1.3E-15 in every row, the
velocity column is unchanged to the last bit, and the pressure column is the
one the same mapping writes without the option (the 2.2E-15 against the
source is the mapping's own log10 round trip, not the option).

## 2. Which state seeds which case

`pick_seed.py` now searches five tiers in order and prints which one it took
from. A state THIS code certified is a better seed than one the code of
2026-08-30 wrote, whatever its composition: every archived state is the
answer to equations that have since changed (`MODELS.md` section 5), while a
certified case of `models/` is a fixed point of the binary the new case will
be solved with.

| tier | what it is |
|---|---|
| `tier1` | a CERTIFIED case of the SAME GROUP, nearest in \|log10 He/H\| |
| `tier2` | a CERTIFIED case of another group with the same physics: same spectrum, same kind of lower boundary, same diffusion, same `K_zz` |
| `tier3` | the same at another XUV normalization of the same star and the same spectral shape |
| `tier4` | an archived state of the same physics, newest code generation first |
| `tier5` | an archived state at another XUV normalization |

CERTIFIED means what the run itself said: the
`- certification of the state written: **CERTIFIED**` line of the case's
`REPRODUCE.md`, or, where that file is absent, a `run.log` whose last
stationary verdict is `info = 0` and whose last certification block is
`CERTIFIED:` and not `NOT CERTIFIED:`.

MEASURED, the lines the tool writes now:

```
$ models/pick_seed.py atomic_scalar_gj1132_kzz1e9/HeH2.13
.../models/atomic_scalar_gj1132_kzz1e9/HeH1.60/output  tier1:models/atomic_scalar_gj1132_kzz1e9  HeH=1.6  target=2.13  dlog10=0.1243  candidates=1  (a certified case of the same group)

$ models/pick_seed.py atomic_scalar_gj1132_kzz1e9/HeH1.60
.../archive_20260830/exhale/refresh_j96/basemetals/a_heh/output  tier4:archive/refresh_j96  HeH=2.13  target=1.6  dlog10=0.1243  candidates=50  (an archived state)

$ models/pick_seed.py molecular_scalar_gj1132_kzz1e9/HeH2.13
pick_seed: no certified case and no archived state carries this physics
(spectrum lhs1140_sed_gj1132_at_b.txt, scalar boundary, diffusion True,
K_zz 1000000000.0): start cold          [exit status 3]
```

A case never seeds itself, and a case whose group holds no certified state
still falls through to the archive, so any case can still be run on its own.
What the tiers change is that the ORDER of a ladder now matters for the cost
of the run, and `models/README.md` says so.

## 3. Why the archived seed failed at He/H = 2.13 and not at 1.60

Both runs used the same binary (md5 `97e10317a710b9ccc63addbedde3586a`), the
same seed (`archive_20260830/exhale/refresh_j96/basemetals/a_heh/output`, a
run at He/H = 2.13) and input files that differ in the planet name and the
`He/H number ratio` line and nothing else. The archived state carries no
metadata block, so it is loaded rather than refused, and `load_IC`'s
`He_diffusion` branch sets its base rows to the run's He/H and keeps the
diffused profile above (top cell He/H = 0.47077 in both runs).

The failure is not in the state as loaded: the two runs measure it at the
SAME `||R|| = 1.9478`, with the same worst row (energy of cell 124, r =
1.0386) and the same window breakdown. What the archived state is, is far
from anything the current code solves near the base -- MEASURED, its first
rows sit at T = 226.0 K (the prescribed base value carried up) with rho =
9.71E+13 and a two-cell velocity sawtooth, v = +0.07, +0.07, -643, +20 cm/s,
while the certified solution of the current code has T = 395, 390, 394, 402
K, rho = 5.3 to 4.7E+13 and |v| < 3.1 cm/s over the same rows. The archived
run also had no `Well balanced` key, so its base pressure jump is the one
the plain HLLC contact speed builds, which on this planet is set by the
reconstruction's truncation error rather than by the flow
(`docs/lhs1140b_stationary_L5c_20260913.md`).

So both solves start O(1) away from stationarity below r = 1.03, and whether
they get out of it is a basin question. At He/H = 1.60 the JFNK finds the
descent: `||Fs||2` falls 2.06E+02, 4.90, 6.69E-02, 1.85E-02, 2.78E-03 over
the first five iterations and outer pass 1 ends at `info = 0`, mass
1.06E-09, energy 1.84E-08. At He/H = 2.13 it does not: `||Fs||2` falls to
1.70E-02 and stops there, the worst row is the energy of cell 12 (r = 1.002)
for twenty-five consecutive iterations at full step length, `||R||` moves
only between 1.76 and 1.93, and the element relaxation ends each pass on its
own fixed point with a drift of 1.8E-01. Five outer passes ran that way
(section 4) before the case was stopped; the logs are in
`LHS1140b/models/.stopped/`.

The account is therefore: the archived state is not the current code's base
solution for EITHER composition, and its distance from one is measurable
(temperature high by 70 percent at the first rows, density low by a factor
1.8, a velocity sawtooth of 600 cm/s where the solution has 3). Which of the
two compositions the solver could still descend from it is not something the
seed tells; what the tiers above do is remove the question, by seeding each
rung from a state the same binary has already certified.

## 4. The pass tables

`atomic_scalar_gj1132_kzz1e9/HeH2.13`, from the archived state
(`models/.stopped/HeH2.13_archive_seed_run.log`), 8 threads:

| pass | hydro | mass | momentum | energy | worst element row | s |
|---|---|---|---|---|---|---|
| 1 | info=2 | 1.23E+00 | 6.42E-04 | 1.78E+00 | 2.98E-04 at cell 317 | 547 |
| 2 | info=2 | 1.69E+00 | 6.41E-04 | 1.78E+00 | 3.16E-04 at cell 306 | 256 |
| 3 | info=2 | 1.42E+00 | 6.35E-04 | 1.78E+00 | 2.57E-04 at cell 303 | 385 |
| 4 | info=2 | 1.70E+00 | 6.21E-04 | 1.78E+00 | 1.80E-04 at cell 301 | 358 |
| 5 | info=2 | 1.54E+00 | 6.17E-04 | 1.79E+00 | 1.16E-04 at cell 300 | 330 |

The same case seeded from the certified `HeH1.60` solution through
`--reservoir-HeH 2.13`:

| pass | hydro | mass | momentum | energy | worst element row | s |
|---|---|---|---|---|---|---|
| 1 | info=0 | 1.22E-09 | 1.25E-14 | 2.31E-08 | 6.53E-05 at cell 287 | 184 |
| 2 | info=0 | 1.03E-09 | 1.25E-13 | 2.20E-08 | 3.52E-05 at cell 290 | 193 |
| 3 | info=0 | 7.05E-10 | 9.91E-14 | 2.10E-08 | 1.88E-05 at cell 293 | 183 |
| 4 | info=0 | 1.49E-09 | 1.17E-12 | 3.27E-08 | 9.99E-06 at cell 295 | 233 |

Outer pass 4 ACCEPTED: every active equation within its own tolerance. The
run wrote `info = 0`, `||R|| = 3.265E-08`, CERTIFIED, log10 Mdot [g/s] =
7.88 and a He I 10830 red-pair equivalent width of 1.3934 %A over 10832.60
to 10834.20 A (air); 13m36s end to end at 8 threads. The seed is
`tier1:models/atomic_scalar_gj1132_kzz1e9`, `HeH=1.6`, carried onto 2.13 by
the factor 1.331250000000.

The hydrodynamic rows are at their tolerances from the FIRST pass, which is
what the seed buys: the four passes are the element relaxation walking its
He/H partition down 6.53E-05, 3.52E-05, 1.88E-05, 9.99E-06 against its
1.0E-05 gate, the hydrodynamic solve re-converging after each. From the
archived state, by contrast, the hydrodynamic rows never came within nine
decades of their tolerances at all.

`atomic_scalar_gj1132_kzz1e9/HeH1.60` itself, for comparison, from the
archive (`tier4`), 7 passes, all `info = 0`, the element row falling
2.37E-04, 1.37E-04, 8.02E-05, 4.71E-05, 2.76E-05, 1.61E-05, 9.32E-06
against its 1.0E-05 tolerance, 15m14s end to end.

### The next rung, He/H = 4.0 from the certified 2.13: NOT certified

The ladder was carried one step further, `HeH4.0` seeded from the `HeH2.13`
solution just certified (`tier1`, `dlog10 = 0.2737`, the helium carried by
the factor 1.877934...). **It did not certify**, and the reason is not the
seed: the hydrodynamic solve of outer pass 1 reached `info = 0` at
`||R|| = 4.88E-08` from the rescaled state, which is the same verdict the
2.13 rung got. What refused is the joint alternation afterwards.

| pass | hydro | mass | momentum | energy | worst element row | omega | s |
|---|---|---|---|---|---|---|---|
| 1 | info=0 | 2.21E-09 | 1.42E-13 | 5.34E-08 | 8.93E-05 at cell 292 | 0.500 | 685 |
| 2 | info=2 | 2.28E-09 | 1.39E-13 | 5.36E-03 | 4.80E-05 at cell 294 | 0.250 | 375 |
| 3 | info=2 | 2.21E-09 | 1.37E-13 | 6.73E-03 | 3.68E-05 at cell 295 | 0.125 | 376 |
| 4 | info=2 | 2.21E-09 | 1.31E-13 | 7.25E-03 | 3.25E-05 at cell 295 | 0.125 | 344 |

Outer pass 4 REFUSED: the joint distance of the state has not fallen in
three consecutive passes. The final verdict is `info = 1`,
`||R|| = 7.246E-03`, NOT CERTIFIED on two entries, 29m41s at 8 threads:

```
hydrodynamic energy row: row measure  7.246E-03 above  1.0E-06 at cell 500
elemental transport He/H partition: gated row measure  3.248E-05 above  1.0E-05 at cell 295 (a wind cell)
```

Two things are worth recording about it. The energy row that refuses is at
cell 500, r = 29.031, the OUTERMOST physical cell and not the base: the
signed terms there are `R = 6.5E-13` against a scale `8.9E-11` that is the
flux divergence, with heating 8.9E-11 and cooling 4.5E-13, so the row is
measuring a nearly exact cancellation between the flux divergence and the
heating at the outflow boundary. And the element row stalls at 3.2E-05
against its 1.0E-05 gate while the under-relaxation is cut 0.5, 0.25, 0.125,
which is the alternation refusing to advance rather than diverging.

This is the partitioned stationary solver's own behaviour at a larger
composition step, not a property of the seed, and `steady_newton.f90` is
another worker's file. It is reported here and left open. What the rung does
establish is that the rescaled seed puts the hydrodynamic rows at their
tolerances in one pass at He/H = 4.0 as well, a factor 1.88 in helium away
from the state it came from.

## 5. Tests

`src/tests/state_mapper/run.sh` (46 assertions, all PASS; the suite builds
nothing and needs no binary):

| assertion | measured |
|---|---|
| `reservoir_HeH_HeI/HeII/HeIII/HeITR_scaled_by_the_ratio` | 0.0 against the option-free mapping, tol 1E-15 |
| `reservoir_HeH_HI_untouched`, `..._HII_untouched` | 0.0, exact |
| `reservoir_HeH_pressure_unchanged`, `..._velocity_unchanged` | 0.0, exact |
| `reservoir_HeH_base_rows_at_the_new_reservoir` | 2.08E-16, tol 1E-14 |
| `reservoir_HeH_diffused_profile_shape_kept` | 4.44E-16, tol 1E-12 |
| `reservoir_HeH_profile_is_not_uniform` | 4.00 (the fixture's He/H falls by four over the column, so the shape test has something to hold) |
| `reservoir_HeH_line_rewritten_*`, `..._other_elements_kept_*` | the new value in both files, `C/H` untouched |
| `reservoir_HeH_temperature_is_p_over_the_particle_count` | 2.44E-15, tol 1E-13 |
| `reservoir_HeH_density_is_the_species_mass_sum` | 4.44E-16, tol 1E-13 |
| `no_option_data_rows_byte_identical_*` | the data rows of a state WITH a reservoir line, mapped without the option, are byte-identical to those of the same state written without one |
| `no_option_mapping_line_states_no_rescaling` | the `# mapped:` line carries no rescaling clause |
| `reservoir_HeH_without_a_reservoir_line_refused`, `..._nonpositive_refused`, `..._with_HeHp_refused` | exit status 2, nothing written |

RED before the change: the option did not exist, so the three
`--reservoir-HeH` invocations exited on the usage text (status 1) and the
scaling and header assertions had nothing to read. GREEN after: 0 failures.

The pre-existing assertions (identity mapping, the shifted-grid round trip,
the He/H invariance, the three refusals, the ghost extrapolation) still
pass, and `--extrapolate-beyond` was exercised by hand on the `HeH1.60`
state (anchor 20 R_p) and still writes its continuation note and takes the
`# grid` line from the target.

## 6. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00

# the tests
$EX/src/tests/state_mapper/run.sh

# which state seeds a case, and from which tier
python3 $EX/LHS1140b/models/pick_seed.py atomic_scalar_gj1132_kzz1e9/HeH2.13

# the rescaling by hand
python3 $EX/src/utils/map_state_to_grid.py \
    $EX/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH1.60/output \
    $EX/LHS1140b/models/current_grid_Hydro_ioniz.txt /tmp/seed --ic \
    --reservoir-HeH 2.13

# the case, end to end (the runner adds --reservoir-HeH by itself)
cd $EX/LHS1140b/models
OMP_NUM_THREADS=8 ./run_case.sh atomic_scalar_gj1132_kzz1e9/HeH2.13
OMP_NUM_THREADS=8 ./run_case.sh atomic_scalar_gj1132_kzz1e9/HeH4.0
```

Files changed by this item:

| file | what |
|---|---|
| `src/utils/map_state_to_grid.py` | the `--reservoir-HeH` option |
| `src/tests/state_mapper/state_mapper_tests.py`, `run.sh` | its tests |
| `LHS1140b/models/pick_seed.py` | the five tiers, the certification reader, the tier on the line |
| `LHS1140b/models/run_case.sh` | the option added when the seed reservoir is not the case's |
| `LHS1140b/models/write_reproduce.py` | the seed section: the tier fields and the rescaling sentence |
| `LHS1140b/models/status.py` | the seed regex, which matched only archive paths |
| `LHS1140b/models/README.md`, `LHS1140b/MODELS.md` | the seeding as it now is, and that the order of a ladder decides the cost |

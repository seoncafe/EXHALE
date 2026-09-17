# L7: the molecular seed built from a certified atomic solution

Item L7 of `docs/PLAN_20260913_lhs_stationary.md`, section 5, as the review of
2026-09-13 (`docs/PLAN_20260913_lhs_stationary_review.md`, sections 2 to 5)
corrected it, and approved by the user at 17:20 as one bounded pilot: the
initialization product, and the isolated reproduction of the two molecular
blockers L4c recorded.

Every number is MEASURED on this tree with the private build
`EXHALE_L7.x` unless it says READ.

## 1. Verdict

**The seed exists, is exact, and is accepted by the restart contract; it does
not by itself open the molecular route.**

1. The conversion is a mode of the binary (`EXHALE_MOLECULAR_SEED`), not a
   Python re-implementation: it calls the production census, equation of
   state, boundary and metadata routines, and writes
   `output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt` with the
   target run's own option, species and reservoir metadata and a
   `# molecular_partition:` line in both halves. Every acceptance test of
   plan section 5.4 that concerns the product passes: the conversion identity
   to 2.6e-15, the H-nucleus, helium and mass budgets to 4.3e-16, the
   equation-of-state closure to 7e-16, and the loader takes the pair with no
   option permitted to differ. On the entry-text binary the same pair is
   refused, by name, in `mol`, `molbase` and `carrier`.
2. **The seed is far from the molecular stationary state, and the distance is
   thermal, not numerical.** At `x2` = 0.99998 -- what the photochemical
   handoff states -- the partition takes every neutral hydrogen atom the wind
   has, out to 30 R_p, so the certified atomic wind's own heating channel is
   gone: the energy row of the seed stands at 1.995 of its own scale at
   r = 1.165 R_p with a flux divergence of 7.3e-5 and `heat - cool` =
   -1.06e-1, and the H2 carrier row at 0.906 in the wind. The partitioned
   stationary route from it does not descend.
3. **What stops it is neither blocker.** With `EXHALE_MOLECULAR_SEED_X2=0` or
   `=local` the hydrodynamic rows of the molecular problem certify in the
   first JFNK and hold through five to nine outer passes on both certified
   seeds; the one row that refuses, on both planets and at every pass, is the
   carrier balance of H2, at exactly 1.000, in the 1.2 to 5 R_p band the wind
   heats to 4500 to 7500 K and where H2 is physically absent. That row is a
   relative measure with no absolute floor, so a density of 1e-27 cm^-3
   gates the solution (section 7.3).
4. **The two blockers are different in kind, and neither is on this path.**
   The ghost reservoir failure is an over-prescription that closes in one line
   of arithmetic and is reproduced here with its numbers; the caloric ceiling
   is a bracket set by a boolean while the term it guards is at most 1.9e-4 of
   the energy of the cell it refuses. Sections 5 and 6 give both
   reproductions and the grounds. Neither was changed: both touch operators
   shared with the regression matrix, and the plan asks for the decision and
   its evidence first.

## 2. What was built

New file `src/modules/init/molecular_seed_from_atomic_state.f90` (module
`molecular_seed`), and four small changes that wire it in.

| file | change |
|---|---|
| `src/modules/init/molecular_seed_from_atomic_state.f90` | new. The partition, the thermodynamic invariant, the checks, the two record lines |
| `src/modules/files_IO/load_IC.f90` | the restart pair's path comes from `molecular_seed_state_file` (unchanged for an ordinary restart); the six molecular option tokens may differ during a seed conversion and are reported by name |
| `src/modules/files_IO/write_output.f90` | a seed run writes the pair under the names a restart reads; both halves carry the seed's record lines |
| `src/EXHALE_main.f90` | `molecular_seed_configure()` before `init`, and the conversion, the write and the stop after it |
| `Makefile` | the new source in `SRC` |
| `docs/input_schema.md` | appendix D.3, the whole contract of the mode |

How it is asked for, in full: `EXHALE_MOLECULAR_SEED=<atomic output
directory>` with the TARGET molecular `input.inp` and `Load IC? True`;
`EXHALE_MOLECULAR_SEED_INVARIANT` = `p` (default) or `T`;
`EXHALE_MOLECULAR_SEED_X2` states the hydrogen-nucleus fraction directly in
place of the handoff's.

**The partition.** `x2 = 2 q (1 + He/H)/(1 + q)` with `q` = `q_H2_base` of the
handoff -- `base_h2_nuclei_fraction()`, the single definition the base ghost
and the equation of state already share -- applied to each cell's own
hydrogen-nucleus budget from the production census
(`element_nuclei_and_charge`) and capped by that cell's neutral hydrogen:

```
delta = min(x2 n_H,nuclei , n_HI),   n(H2) += delta/2,   n(HI) -= delta
```

which is the form `set_IC` uses. It is carried out on the mass fractions
`f = n/(rho n0)`, so `2 x 0.5 delta = delta` makes the H nuclei and (with
`bsp_mass(H2) = 2`) the mass exact rather than round-off-exact. `H2+`, `H3+`
and `HeH+` are written at zero. A run with no handoff is refused: the
alternative source is a local chemical-equilibrium fit, which is the seed's
own estimate and not upstream information.

**The admissibility statements, as the review asked.** `q` is compared with
its ceiling `0.5/(0.5 + He/H)` (MEASURED for the pilot: 0.23809027 against
0.23809524), `x2` with unity, and the transfer with the hydrogen each cell
actually has. The mixing-ratio ceiling is never applied to `x2`, which runs to
1 over that range.

**The thermodynamic invariant, stated in the file.** Two H atoms into one H2
keep the nuclei and the mass but remove half a particle per pair, so `rho`,
`v`, `p` and `T` cannot all be kept. `invariant=p` keeps `rho, v, p` and takes
`T = p/(n_tot + n_e)` from the new census; `invariant=T` keeps `rho, v, T` and
takes `p = (n_tot + n_e) T`. Both close the energy with the production caloric
equation of state, re-apply the lower and outer boundary so the ghost rows are
the target run's own, and measure the pressure round trip through
`W_to_U`/`U_to_W`.

## 3. The tests

T-L7-1, T-L7-2 and T-L7-3 are automated as `src/tests/molecular_seed/run.sh`
(fixture: the hot Uranus of `backup/regression/mol_base_handoff`, whose atomic
counterpart is the same input without the molecular chemistry and without the
handoff's `q_H2_base`, so the two configurations differ in the molecular
tokens alone). MEASURED, all seventeen assertions PASS:

| test | assertion | measured | tolerance |
|---|---|---|---|
| T-L7-1 | `x2 = 0`: every shared column of both files, physical cells | 2.63e-15 | 1e-12 |
| T-L7-1 | the four molecular columns are exactly zero | 0 | 0 |
| T-L7-2 | equation-of-state closure at `x2 = 0` | 3.72e-16 | 1e-12 |
| T-L7-2 | H nuclei / He nuclei / mass at the handoff's `x2` | 2.21e-16 / 0 / 4.32e-16 | 1e-12 |
| T-L7-2 | equation-of-state closure at the handoff's `x2` | 6.96e-16 | 1e-12 |
| T-L7-2 | largest `x2` actually transferred | 0.985448 | <= 1 |
| T-L7-3 | the record lines are in both halves | present | -- |
| T-L7-3 | the pair reloads into the target molecular run | exit 0 | -- |
| T-L7-3 | species columns of the written pair | 38 | 38 |
| T-L7-3 | the coupling line | `certified=F cert_reason=molecular_seed mode=init` | -- |
| refusal | an atomic target with the seed variable set | refused | -- |
| refusal | a molecular target with no `q_H2_base` handoff | refused | -- |

RED, on the entry text: with the certified atomic pair placed in `output/`,
the campaign binary `EXHALE.x` (md5 `97e10317a710b9ccc63addbedde3586a`) stops
at startup with

```
(load_IC) ERROR: metadata field "options" differs between the restart files and this run in
  3 token(s) the input did not name as allowed to change:
  mol: F -> T; molbase: F -> T; carrier: F -> T;
```

which is the refusal L4c recorded and the reason the molecular groups had no
route.

On the LHS pilot itself (seed `atomic_scalar_gj1132_kzz1e9/HeH1.60/output`,
CERTIFIED 2026-09-13 14:23, into the target case built at the same He/H,
section 8), MEASURED over the physical cells at `x2 = 0`: rho 1.63e-14,
v 2.30e-16, p 2.65e-16, T 5.27e-16, and every species column at or below
6.94e-16. The ghost rows differ because the conversion re-applies the target
run's own base boundary: the base ghost velocity goes from -3.04 cm/s (the
atomic solution's) to -578 cm/s, which is a statement of the molecular base
boundary and not of the transfer.

## 4. The two invariants, on one seed

Seed `atomic_scalar_gj1132_kzz1e9/HeH1.60/output`, target
`models/.L7/molecular_kzz1e9_HeH1.60` (He/H = 1.60, `q_H2_base` = 0.23809027,
`x2` = 0.99998316, 473 of 504 cells capped by their own neutral hydrogen).

| | invariant `p` | invariant `T` |
|---|---|---|
| largest relative move of T | 3.116e-1 | 0 |
| largest relative move of p | 0 | 3.116e-1 |
| T of the physical column [K] | 463.5 to 7312 | 393.5 to 5333 (the atomic solution's, unchanged) |
| equation-of-state closure | 6.28e-16 | 8.05e-16 |

The 31 percent is the particle count: at He/H = 1.60 the hydrogen is 38
percent of the nuclei and binding essentially all of it into H2 removes 19
percent of the heavy particles, so at fixed pressure the temperature rises by
about the same fraction. Neither invariant is small, and neither preserves the
chemical or the energy equilibrium of the source, which is what plan section
5.2 item 3a says a seed does not promise.

## 5. Blocker (a): the caloric domain

**Judgement: the 5e4 K ceiling is a property of the ENERGY UPDATE's bracket,
not of the caloric equation of state, and it is selected by a boolean while
the model error it guards is proportional to an abundance.**

READ, from the source:

- `caloric_eos.f90` (`caloric_state_from_composition`, near line 215) marks a
  cell molecular when `n(H2) > 0` -- any positive value.
- `energy_semi_implicit.f90` gives a molecular cell the bracket ceiling
  `T_ceiling_mol_K` = 5.0e4 K and any other cell 1.0e7 K.
- The H2 rovibrational table runs from `T_utab_lo` = 1 K to `T_utab_hi` =
  5.0e4 K, 4096 nodes -- so the ceiling is the top of the table.
- But `h2_rovibrational_energy_and_heat_capacity` does NOT stop there: above
  `T_utab_hi` it continues LINEARLY in T with the end-point heat capacity,
  `u_rv = u(T_hi) + c(T_hi)(T - T_hi)`, `c_rv = c(T_hi)`. The pair is
  internally consistent (`du/dT = c` exactly), monotone in T -- the property
  the energy-to-temperature inverse relies on and the module's own comment
  names -- and it carries no dissociation energy, so extending it cannot
  double count one. `internal_energy_of_mixture` weights it by `x2 =
  n(H2)/(n_tot + n_e)`.

So a cell at 7e4 K holding a trace of H2 has a perfectly well defined,
invertible internal energy; what refuses it is the bracket. A cell with
`x2 = 1e-12` is thermodynamically atomic to twelve digits and is nevertheless
given a ceiling 200 times lower than the atomic one. That asymmetry is the
defect: the guard is a step function of `n(H2) > 0` while the quantity it
guards enters `u` multiplied by `x2`.

MEASURED on the pilot seed: the ceiling is NOT reached. The certified atomic
LHS 1140 b wind at He/H = 1.60 spans 393.5 to 5333 K, and the seed spans 463.5
to 7312 K at invariant `p`; no cell of either exceeds 5e4 K, so no cell is
both molecular and above the ceiling. The blocker therefore does not arise on
this seed and is not what stops the route (section 7).

### 5.1 The failure, reproduced and weighed

Where it does arise is the helium-poor end of the ladder, which has no
certified atomic solution to seed from. MEASURED with `EXHALE_L7.x`, the cold
march of `molecular_scalar_gj1132_wellmixed/HeH0.083`, 8 threads, reproduces
the stop of L4c at step 983 (L4c READ 980):

```
energy update FAILURE, step 983, 4 cell(s), first cell 228 (T  7.07E+04 K)
  reason: no bracket below the equation-of-state ceiling, |R|/scale  2.9227E-01,
          last iterate  5.0000E+04 K, missing source -9.7085E-03 erg cm^-3 s^-1
ATTEMPTED STEP EXHAUSTED at step 983 after 9 attempts.
```

Nine halvings to `dt/2**8` give the same reason, so it is not a step-size
failure: the balance wants a temperature above the ceiling and the bracket
has none.

**How much H2 those cells hold.** MEASURED on the last accepted state
(the same march stopped at step 982, which writes it):

| cell | r [R_p] | T [K] | n(H2) [cm^-3] | `x2` = n(H2)/(n_tot + n_e) | 2 n(H2)/n_H |
|---|---|---|---|---|---|
| 226 | 1.2350 | 5.263e4 | 1.515e3 | 3.906e-4 | 1.692e-3 |
| 227 | 1.2392 | 6.777e4 | 1.241e3 | 4.004e-4 | 1.753e-3 |
| 228 | 1.2434 | 7.186e4 | 1.212e3 | 3.970e-4 | 1.746e-3 |
| 229 | 1.2477 | 7.284e4 | 1.207e3 | 3.873e-4 | 1.706e-3 |
| 230 | 1.2521 | 7.036e4 | 1.181e3 | 3.698e-4 | 1.624e-3 |
| 231 | 1.2565 | 5.934e4 | 8.566e2 | 2.982e-4 | 1.296e-3 |

Six cells, in a 0.02 R_p band above the H2 front, carry four parts in ten
thousand of the particle count as H2 and nothing else about them is
molecular. The molecular term of the internal energy is `x2 u_rv/T0`
(`internal_energy_of_mixture`) against `1.5 T`, and the mean energy of a
BOUND rovibrational ladder cannot exceed its dissociation energy,
`D0/k` = 5.20e4 K, so that term is at most

```
x2 * u_rv / (1.5 T) <= 4.0e-4 * 5.20e4 / (1.5 * 7.19e4) = 1.9e-4
```

of the cell's internal energy. **The bracket that refuses the step, and with
it the run, is guarding a term worth less than two parts in ten thousand of
the energy of the cell it refuses.** That is the asymmetry stated as a
number.

### 5.2 The three options, and the grounds

The options of plan section 5.5, with the grounds these measurements give:

1. **A consistent high-temperature extension.** It already exists: the linear
   continuation above the table is exactly that, with a consistent `(u, c_v)`
   pair and no dissociation energy. Adopting it means setting the molecular
   ceiling to the atomic one and writing the validity statement at the code
   site: above 5e4 K the bound-ladder energy of H2 is an extrapolation whose
   weight in `u` is `x2`, and H2 does not survive at those temperatures, so
   what must remove it is the chemistry and not the bracket. Cheapest, and it
   removes a step function; it does not by itself make the hot cell's
   chemistry right.
2. **Change the coupled chemistry/energy solve** so that a cell heated past
   the domain dissociates its H2 within the same step instead of being
   refused. Physically the most nearly right and the largest change.
3. **Refuse explicitly, outside a documented validity domain.** This is what
   happens today, except that the refusal is phrased as a missing bracket and
   names no domain. If this is kept, the message should say that the state
   left the domain of the molecular caloric model and at which abundance.

An abundance cutoff is not among them: the review is right that it would put a
discontinuity into `u` and break the energy consistency.

**Not implemented here.** `energy_semi_implicit.f90` sets the accepted-state
domain of an operator every molecular case in the regression matrix runs
through, so the change needs the scoped impact measurement of the worker rules
(`mol_base_handoff`, `mol_metals`, `mol_lyman_werner`, `mol_diffusion`,
`mol_ir_bands`, `mol_sec_ion`, `mol_carrier`), which this item is not the
place to run. It is a one-constant change plus a comment once that measurement
is made, and option 1 is the recommendation.

## 6. Blocker (b): the ghost reservoir chemistry

**Judgement: the boundary-imposed composition is NOT what has no root. At the
certified base state the prescription is solvable, and it stays solvable when
the room it leaves is moved by two decades. What had no root in L4c is the
MARCHED ghost, whose ionization had grown to 25 times the room the
prescription leaves; the prescription carries a feasibility condition that
nothing in the code states or checks.**

READ, `ionization_equilibrium.f90` near line 1951: for the two lower ghosts
(`j <= 0`), when a handoff stated the partition,

```fortran
ieq_cell%x_h2_fixed = .true.
ieq_cell%x_h2_fix   = base_h2_nuclei_fraction()
```

and the ionization stages are solved against it. What is pinned is the H2
share of ALL the hydrogen nuclei, so what is left for H I, H+ and the H bound
in H2+, H3+ and HeH+ is `1 - x2` -- 1.68e-5 at the pilot's handoff. Nothing
compares the ionization the ghost's own balance produces with that room.

The arithmetic of the two states:

| ghost -1 | T [K] | x2 imposed | room `1 - x2` | x(H+) | ratio |
|---|---|---|---|---|---|
| the certified atomic wind of He/H = 1.60 (MEASURED) | 395 | 0.99998316 | 1.684e-5 | 4.402e-6 | 0.26 |
| L4c's marched state at its stop (READ) | 271 | 0.99998316 | 1.684e-5 | 4.29e-4 | 25 |

and the evidence in the failing state is what the solver returned there: L4c
READ `2 n(H2)/n_H` = 0.99956 against the imposed 0.99998, i.e. 4.2e-4 of the
hydrogen taken back out of H2 to pay for the ionization -- a departure from
the imposed constraint, not from the element identity, and the run reported it
as an element violation without naming the imposed value or the room.

### 6.1 The failure, reproduced and taken apart

MEASURED with `EXHALE_L7.x`, the cold march of the pilot's own molecular case
(`models/.L7/molecular_kzz1e9_HeH1.60`, He/H = 1.60, 8 threads): the stop of
L4c reproduces at the same cell, on the same statement, at step 1116.

```
(ioniz_eq) STOP: a cell has rested on a NON-ROOT chemical equilibrium beyond the relaxation amnesty
  cell -1  step 1116  consecutive non-root sweeps 1000
  r [R_p] .999807  T [K]  2.7075E+02  n_tot [cm^-3]  3.5795E+10  n_e [cm^-3]  3.0866E+06
  solver info 1  element violation  3.115E-04  normalized reaction residual  3.114E-04
  candidate stage fractions:
 3.2353E-04  1.9294E-04  9.3485E-08  9.9967E-01  2.2481E-08  4.6434E-06
 8.5760E-08  4.8106E-08
```

READ, `System_HeH_mol.f90` lines 7 to 12, the unknowns those eight numbers
are: `x(1)` = n(H+)/n_H, `x(2)` = n(He+)/n_He, `x(3)` = n(He++)/n_He,
`x(4)` = 2 n(H2)/n_H, `x(5)` = 2 n(H2+)/n_H, `x(6)` = 3 n(H3+)/n_H,
`x(7)` = n(HeH+)/n_H, `x(8)` = n(He 2^3S)/n_He; and line 576,
`fvec(4) = x(4) - ieq_cell%x_h2_fix`, is the pinned row itself. The
arithmetic of the failure is then closed:

| quantity | value |
|---|---|
| imposed `x(4)` | 0.9999831600354949 |
| room the pin leaves, `1 - x(4)`, for `x(1)`, `x(5)`, `x(6)`, `x(7)` | 1.684e-5 |
| `x(1)` the ionization balance asks for | 3.2353e-4 |
| ratio | 19.2 |
| `x(1) + x(4)` at the imposed value | 1.000307, i.e. 3.07e-4 beyond the hydrogen the cell has |
| `x(4)` the candidate returned | 0.99967 |
| its departure from the pinned value | 3.132e-4 |
| the violation the run reported | 3.115e-4 |

**The reported violation IS the amount by which the pinned row had to be
broken to make room for the ionization**, to half a percent. The candidate is
inside the simplex (`n_HI = (1 - x1 - x4 - x5 - x6 - x7) n_H` = 7.5e-7 n_H,
positive) precisely because `x(4)` came down; what it is not is a root of the
pinned row, which is why the sweep reports a non-root for a thousand steps and
the run stops. No tolerance and no amnesty can repair it: the system is
infeasible, not badly converged.

### 6.2 Varying only the room, at a stationary base

The prescription was held exactly as the operator holds it -- the same
`base_h2_nuclei_fraction()` at the same two ghosts, the same base state, the
same thermodynamics -- and only the room it leaves was varied, by building
three seeds from the one certified atomic solution at three `q_H2_base`
values and equilibrating each (`EXHALE_RESIDUAL=1`, which puts the loaded
state on its chemical root before it measures anything):

| `f` = x2 | `q_H2_base` | room `1 - x2` | sweeps | cells with no admissible root | sweeps in which every cell had one |
|---|---|---|---|---|---|
| 0.99990 | 0.23806576 | 1.000e-4 | 65 | 0 | 65 of 65 |
| 0.99998316 | 0.23809027 | 1.684e-5 | 65 | 0 | 65 of 65 |
| 0.999999 | 0.23809494 | 1.000e-6 | 65 | 0 | 65 of 65 |

MEASURED. Moving the room by two decades, through and below the ghost's own
ionization fraction, changes nothing: no cell of any of the three lacks a
root, and the number of cells whose molecular equilibrium root leaves the
physical simplex and is clamped onto the element budget is the same 63 in all
three (those are wind cells, where the seed's imposed H2 is far from its
balance, not the ghosts). The same holds through the whole stationary route
of section 7: `no admissible root at 0 cell(s)` in every one of its sweeps.

So the proximity of `q_H2_base` to its ceiling is refuted a second time -- L4c
refuted it for the march, this refutes it for the boundary-imposed state
itself -- and the ghost failure is a property of the state the march had
reached, not of the prescription applied to a stationary base.

What the prescription still needs, and what it does not have:

- it is feasible only while `x(H+) + x(H in H2+, H3+, HeH+) <= 1 - x2`, and
  nothing evaluates that condition. At `f` = 0.99998 the margin is 1.7e-5, so
  a hydrogen ionization fraction above 1.7e-5 anywhere in the reservoir makes
  the constraint unsatisfiable while every message still speaks of a
  non-root;
- the distinction the review asks for is therefore: the boundary MAY prescribe
  the molecular partition of the gas the lower atmosphere hands over -- the
  hydrogen that the wind's own field has not ionized -- and MUST solve the
  ionization of it, and hence the H2 share of the TOTAL hydrogen, which is
  `x2 (1 - x_H+ - ...)` and not `x2`. The two agree to 1e-5 in a cold neutral
  reservoir, which is why the constraint has worked where it has;
- either way, the refusal should name the cell, the imposed value, the room
  and the ionization that consumed it, instead of reporting only a streak.

Not implemented here: `ionization_equilibrium.f90` is the operator every
molecular case in the regression matrix runs through, and the change is a
change of what the lower boundary means. It is the next item, with its own
impact measurement.

## 7. What the seed does and does not open

The seed is accepted, equilibrates, and every cell has a chemical root; what
it does not do is put the wind near its molecular stationary state.

MEASURED, the certification of the state as loaded (invariant `p`,
`Restart intent: stationary equilibrate`, `Well balanced: True`,
`Secondary_ionization: Immediate`, `Solver: Newton`, `EXHALE_PTC_DTAU0=1.0`,
the carrier partitioned route, 8 threads):

| row | measure | tolerance | where |
|---|---|---|---|
| hydrodynamic mass | 1.000 | 3.0e-12 | cell 1 |
| hydrodynamic momentum | 2.121e-1 | 1.0e-8 | cell 1 |
| hydrodynamic energy | 1.995 | 1.0e-6 | cell 206, r = 1.1654 |
| carrier balance H2 (gated, r >= 1.2) | 9.057e-1 | 1.0e-5 | cell 226 |
| elemental transport He/H (gated) | 2.381e-2 | 1.0e-5 | cell 373 |
| level balance He 2^3S | 3.700e-12 | 1.0e-6 | within |
| eliminated-species closure `System_HeH_mol` | 6.746e-7 | 1.0e-6 | within |
| cells without a chemical root | 0 | -- | -- |

and the signed terms of the worst cell: `R` = 1.065e-1, flux divergence
7.25e-5, `heat` = -5.337e-2, `cool` = +5.303e-2, `heat - cool` = -1.064e-1.
**The energy row of the seed is its radiative imbalance, not its flux
divergence** -- a factor 1.5e3 between the two. That is what the partition
did: taking every neutral hydrogen atom out to 30 R_p removes the H I
photoionization heating the certified atomic wind is held up by, and the
temperature the invariant then assigns is not the temperature that balances
what is left.

MEASURED, T-L7-5 at the handoff's own `x2`, both invariants, 45 minutes each
at 8 threads: neither reaches a second outer pass. The first JFNK runs its
whole budget with `||R||` between 1.90 and 2.86 (invariant `p`: 1.995 at the
start, 15 judged iterates, last 2.01; invariant `T`: 1.901 at the start, 11
judged iterates, last 2.06), the hydrodynamic rows carrying all of it and the
element and carrier rows none. That is the same behaviour L4c measured from a
marched snapshot, and for the same reason: the state handed to the solve is
not near the stationary one. It is NOT the base-layer discretization L5c
removed -- `Well balanced: True` is on -- and it is not a chemical failure,
because every cell has a root throughout.

### 7.1 The other end of the partition, and why it is not the answer either

The conversion supports the opposite choice directly: `x2 = 0` writes the
certified atomic wind with the four molecular columns at zero and the
molecular base boundary applied. MEASURED on the same seed and route
(He/H = 1.60, 45 minutes, 8 threads), it is a different régime altogether:

| | `x2` = 0.99998 (handoff) | `x2` = 0 |
|---|---|---|
| certification as loaded: refusing entries | 5 | 4 |
| hydrodynamic energy row as loaded | 1.995 at cell 206 | 1.063 at cell 2 |
| elemental transport as loaded (gated) | 2.381e-2 | 4.848e-6, WITHIN |
| first JFNK | no descent in 45 min, `||R||` 1.90 to 2.86 | `info = 0` at `||R||` 4.0e-8 in 406 s |
| outer passes completed in 45 min | 0 | 6, every one `info = 0` |
| hydrodynamic rows at pass 6 | -- | mass 8.4e-7, momentum 2.5e-10, energy 3.0e-6 of their own scales |
| the one refusing row at every pass | -- | carrier balance H2, **1.000 exactly**, at every cell |

So the hydrodynamics of the molecular problem is reachable from a certified
atomic wind, and in one JFNK. What is not reachable is the H2 column, and the
reason is in the code:

READ, `ionization_equilibrium.f90` lines 1897 to 1902 -- with the carrier
transported AND the state loaded, the sweep does not solve the H2 row, it is
HANDED the loaded value:

```fortran
if (thereis_mol .and. carrier_transport                   &
    .and. (bg_ready .or. do_load_IC)) then
   ieq_cell%x_h2_fixed = .true.
   ieq_cell%x_h2_fix  = 2.0d0*nmol_eq(j,1)/nh(j)
```

with the comment that re-solving it locally "would throw the transported
partition away". A seed whose H2 column is identically zero therefore pins
`x_h2 = 0` in every cell for the whole run, the carrier balance row is its own
residual over its own (empty) terms -- 1.000 by construction, at every cell,
in every pass -- and no number of outer passes can move it. **Zero is a fixed
point of the loaded transported carrier**, which is the same trap `set_IC`
warns about for the cold start ("a zero-H2 start makes the base solve land on
/ fail into the spurious atomic root").

### 7.2 The same on the second certified seed

`atomic_scalar_gj1132_kzz1e9/HeH2.13` certified at 2026-09-13 18:00
(`info = 0`, `||R||` 3.3e-8, log Mdot 7.88), so T-L7-5 was repeated on the
case the plan names, `molecular_scalar_gj1132_kzz1e9/HeH2.13`, with all three
partitions. MEASURED:

| partition | seed: `x2` requested | largest T move | first JFNK | outer passes | the refusing row |
|---|---|---|---|---|---|
| handoff | 0.9999832, 480 of 504 cells capped | 27.6 % | no descent, `||R||` 1.93 to 2.73 in 14 min | 0 | -- |
| `local` | per cell, 234 cells clipped at the ceiling | 18.2 % | `info = 0` | 5, every one `info = 0`; REFUSED at pass 4, no joint progress | carrier balance H2, 1.000 at cell 340 |
| `x2 = 0` | 0 | 0 | `info = 0` | 9 in 60 min (the cap), every one `info = 0` | carrier balance H2, 1.000 at cell 340 |

The two planets agree in every respect, so the behaviour is the partition and
the route, not one case.

### 7.3 What actually stops it: the carrier row has no absolute floor

Cell 340 is at r = 2.72 R_p, T = 3.3e3 to 4.6e3 K. MEASURED, the H2 there:

| state | n(H2) [cm^-3] at cell 340 |
|---|---|
| the `local` and `x2 = 0` seeds | 0 |
| after five outer passes of the `local` route | 5.6e-27 |
| the seed at the handoff's `x2` | 7.7e6 |

The carrier solve does move H2 in the band the local partition leaves empty --
it grows from an exact zero to 1e-29 to 1e-27 cm^-3 -- and the row measure
stays exactly 1.000 the whole way, because **the carrier balance row is a
RELATIVE measure**: its scale is "the sum of the row's own terms" and a column
at 1e-27 cm^-3 that is far from its own balance reads the same as a column at
1e12 cm^-3 that is far from its own. Against a gas density of 1e6 cm^-3 at
that radius, an H2 density of 1e-27 is not a species; it is 33 decades of
round-off. The element rows carry an absolute floor for exactly this reason
(`floor 1e-20 rho X_base/(R0/v0)`, the certification prints it); the carrier
row's floor is stated as "the 1e-20 free-e" and does not reach this regime.

So the molecular route from a certified atomic seed is gated by a row whose
value is O(1) wherever the molecule is physically absent. The hot Uranus of
`backup/regression/carrier_model_a_newton` certifies under the same gate
because its H2 is present everywhere its rows are judged; the LHS 1140 b wind
has a 1.2 to 5 R_p band that is 4500 to 7500 K and simply has no H2 in it.

### 7.4 What the pilot hands on

Three items, in the order they bind:

1. **The carrier balance row needs an absolute floor**, on the same footing as
   the element rows' `1e-20 rho X_base`, below which a carrier row is reported
   and does not gate. Without it no molecular LHS 1140 b solution can certify,
   whatever the seed. This is the first item, and it is a certification item,
   not a seed item.
2. **The partition between the two ends.** The handoff's `x2` carried to
   30 R_p leaves the thermal state an order of unity away and the solve does
   not descend; `x2 = 0` puts the hydrodynamics exactly where it belongs and
   pins the chemistry at zero. `EXHALE_MOLECULAR_SEED_X2=local` is the
   physical middle and reaches five `info = 0` passes; which of the three the
   route finally takes is settled only after item 1. Note the limitation the
   fit carries: it is a chemical-equilibrium estimate, and the outer wind at
   500 K and 1e-2 cm^-3 is not in chemical equilibrium with the
   Lyman-Werner field, so the `x2` it states above about 10 R_p is an
   initialization value and nothing more.
3. The two blockers of sections 5 and 6, neither of which is on the path of a
   seeded LHS 1140 b case at He/H >= 1.6 (the caloric ceiling is never reached
   there, and the ghost prescription is feasible at a stationary base), and
   both of which are on the path of the He-poor and the cold-march cases.

## 8. How to reproduce

`EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00`, from
`$EX/LHS1140b/models`. The private build is
`make OBJDIR=build_L7 EXE=EXHALE_L7.x`.

Two targets were used. `molecular_scalar_gj1132_kzz1e9/HeH2.13` is the case
the plan names, seeded from `atomic_scalar_gj1132_kzz1e9/HeH2.13`, which
certified while this item ran. `models/.L7/molecular_kzz1e9_HeH1.60` was built
first, when `HeH1.60` was the one certified atomic solution of the tree; its
`input.inp` and `base.inp` are the `HeH2.13` pair at He/H = 1.60, with
`q_H2_base` = 0.23809027395326968 from the same photochemical column fraction
`f` = 0.9999831600354949. The two agree in every measurement.

```bash
SRC=$EX/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH1.60/output

# the seed, invariant p.  EXHALE_MOLECULAR_SEED_INVARIANT=T keeps T instead;
# EXHALE_MOLECULAR_SEED_X2 takes 0 (the conversion identity), a number, or
# "local" (each cell's own chemical equilibrium).
d=.L7/seed_p; mkdir -p $d/output
cp .L7/molecular_kzz1e9_HeH1.60/base.inp $d/
sed 's/^Load IC?.*/Load IC? True/' .L7/molecular_kzz1e9_HeH1.60/input.inp > $d/input.inp
( cd $d && OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=$SRC $EX/EXHALE_L7.x > seed.log 2>&1 )

# the stationary route from it
s=.L7/route_p; mkdir -p $s/output; cp $d/base.inp $s/
cp $d/output/Hydro_ioniz_IC.txt $d/output/Ion_species_IC.txt $s/output/
sed -e 's/^Reconstruction scheme:.*/Reconstruction scheme: PLM/' \
    -e 's/^Load IC?.*/Load IC? True/' -e '/^du_th /d' \
    .L7/molecular_kzz1e9_HeH1.60/input.inp > $s/input.inp
printf 'Restart intent: stationary equilibrate\nSecondary_ionization: Immediate\n' >> $s/input.inp
( cd $s && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 $EX/EXHALE_L7.x > run.log 2>&1 )

# the same on the case the plan names, from its own certified atomic solution
SRC=$EX/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output
d=.L7/seed213_local; mkdir -p $d/output
cp molecular_scalar_gj1132_kzz1e9/HeH2.13/base.inp $d/
sed 's/^Load IC?.*/Load IC? True/' molecular_scalar_gj1132_kzz1e9/HeH2.13/input.inp > $d/input.inp
( cd $d && OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=$SRC \
     EXHALE_MOLECULAR_SEED_X2=local $EX/EXHALE_L7.x > seed.log 2>&1 )

# the two blocker reproductions: the ghost stop and the caloric ceiling
m=.L7/march_HeH1.60; mkdir -p $m/output
cp .L7/molecular_kzz1e9_HeH1.60/{input.inp,base.inp} $m/
( cd $m && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=3000 $EX/EXHALE_L7.x > march.log 2>&1 )
w=.L7/march_wellmixed; mkdir -p $w/output
cp molecular_scalar_gj1132_wellmixed/HeH0.083/{input.inp,base.inp} $w/
( cd $w && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=1200 $EX/EXHALE_L7.x > march.log 2>&1 )

# the acceptance tests
EXHALE_EXE=$EX/EXHALE_L7.x EXHALE_TEST_OUT=$EX/build_L7/tests/molecular_seed \
    $EX/src/tests/molecular_seed/run.sh
```

## 9. What was changed in the tree

- the six files of section 2 (five source and build files, one schema
  appendix);
- `src/tests/molecular_seed/run.sh` and its `compare_columns.py`, new;
- `LHS1140b/models/run_case.sh.molecular.patch`, rewritten around the route
  the seed makes possible, against `run_case.sh` md5
  `4b3d45d51bdc92edab91acefada9195e` (the version carrying the L11 worker's
  `--reservoir-HeH` seed rescale). It is held as a patch and NOT applied, as
  L4c left it: the atomic campaign re-reads `run_case.sh` for every case;
- `LHS1140b/models/.L7/`, the scratch tree of this item.

No golden was refreshed and no regression case was run: nothing here changes
an atomic run, and `molecular_seed_on()` is false in every run that does not
set the variable. The three source files the conversion touches are entered
only through that test: `load_IC` takes its file path from
`molecular_seed_state_file`, which returns `output/<name>_IC.txt` unchanged
when no seed is being built; the option exception is guarded by the same
flag; and `write_output` opens the same two files it always did.

Noticed outside the scope of this item, reported and not fixed:

- the carrier balance H2 row has no absolute floor (section 7.3), which is
  what gates the molecular LHS 1140 b solution;
- `ionization_equilibrium`'s refusal of an infeasible ghost prescription
  reports a streak and a violation and names neither the imposed value nor
  the room it left (section 6);
- `energy_semi_implicit`'s "no bracket below the equation-of-state ceiling"
  names no validity domain and no abundance (section 5).

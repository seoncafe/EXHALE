# L7b: the three changes the L7 pilot named

Item L7b of `docs/PLAN_20260913_lhs_stationary.md`, the follow-up the pilot
report `docs/lhs1140b_stationary_L7_20260913.md` hands on in its section 7.4.
Three changes, each reproduction first: an absolute floor under the carrier
balance row, the molecular energy bracket, and the base handoff's H2
partition. Every number is MEASURED on this tree with the private build
`EXHALE_L7b.x` unless it says READ.

## 1. Verdict

**All three are defects, all three are fixed, and the two that stop a run
stop it no longer.**

1. **The carrier balance row had no absolute floor under the SPECIES**, only
   under the row's terms, so a species that is physically absent gated the
   solution. It now reports and does not gate below 1e-20 of the free
   reservoir of its own element, which is the constant the elemental
   transport rows already use. On the LHS 1140 b seed the refusing cell moves
   from 340, where n(H2) = 0 exactly, to 218, where n(H2)/n_H = 2.5e-6 and
   the row is a statement about a species that is there.
2. **The `ionization_equilibrium` pin of x_H2 under carrier transport is
   DESIGN, not a defect.** What it pins is not the loaded value: `nmol_eq` is
   rebuilt from `f_sp` at the top of every sweep, so the pin is the value the
   transport operator last wrote, and the H2 source of the transport row
   (`carrier_source`, `fv(4)` of `mol_heh_rows`) contains formation terms
   that do not vanish at n(H2) = 0. Zero is therefore not a fixed point of
   the transported carrier, which the pilot's own measurement shows
   (0 -> 1e-29 -> 1e-27 over five passes). Nothing was changed here.
3. **The molecular energy bracket refused whole steps for a term worth 1.9e-4
   of the energy of the cell it refused.** The 5e4 K ceiling was the top of
   the H2 rovibrational table, selected by the boolean n(H2) > 0, while the
   caloric equation of state continues consistently above it. One ceiling now
   serves every cell, with the validity statement at the code site. The
   He-poor cold march that stopped at step 983 runs its whole 1200 steps, and
   its final state has no cell above 3.94e4 K: what the bracket refused was a
   transient excursion the solve relaxes out of when it is allowed to find
   the root.
4. **The base handoff prescribed the H2 partition of every hydrogen nucleus,
   which is infeasible as soon as the reservoir ionizes.** It now prescribes
   the partition of the NON-IONIZED hydrogen, x(4) = x2 (1 - x_ion), which
   leaves (1 - x2)(1 - x_ion) + x_ion for the ions and is therefore always
   feasible. The cold march that stopped at step 1240 on a cell resting on a
   non-root runs its whole 3000 steps, and its ghost carries
   x(4)/(1 - x_ion) = 0.99998324 against the handoff's x2 = 0.99998316 --
   the prescription is met as stated, at a ghost ionization of 1.9e-2, which
   is 1100 times the room the old form left.

## 2. What changed

| file | change |
|---|---|
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | `carrier_absent_fraction = 1e-20` as a named module constant (it already carried the row-scale floor as a literal, twice); `carrier_residual` and `carrier_steady_residual` return an optional per-cell mark of where a carrier stands below that fraction of the free reservoir of its own element |
| `src/modules/time_step/certification.f90` | the carrier row's gate takes that mark: an absent cell is reported and does not gate, the count and the first such cell are printed, and a column with no cell carrying the species says so instead of reading as satisfied |
| `src/modules/time_step/energy_semi_implicit.f90` | one bracket ceiling for every cell (`T_ceiling_K`), with the validity of the caloric extrapolation above the H2 table stated at the site; `T_ceiling_mol_K`, `molecular_cell` and the branch are gone |
| `src/EXHALE_main.f90` | the same two-line bracket expression at the two sites of the extrapolation ratio test, which had its own copy of the removed branch |
| `src/modules/radiation/ionization_equilibrium.f90` | the lower-boundary reservoir imposes `x2 (1 - x_ion)` on the H2 row, with `x_ion` the ghost's ionized hydrogen-nucleus fraction (H+ and the hydrogen in H2+, H3+, HeH+), and refuses by name a prescription that leaves less room than the ionization it must make room for |

`base_h2_nuclei_fraction()` stays the single definition of the handoff
partition and is not touched; what changed is the quantity it is a fraction
OF, and that statement is made in exactly one place. READ, and it settles the
particle-count question: `ioniz_eq` already refreshes `ntot_bc` from the
ghost's own particle count after every sweep in a molecular run
(`if (thereis_mol) ntot_bc = n_tot(1-Ng)/n0`), so the pressure boundary
condition follows the species state and does not keep the startup constant
that binds x2 of every nucleus.

## 3. The control

Two pairs of binaries were built, both privately
(`make OBJDIR=... EXE=...`), and the report says which pair each number comes
from.

**Pair 1**, `OBJDIR=build_L7b`: the CONTROL is the live tree as this item
found it (md5 `3c33454a0d15870dd44df9ace548d8b6`, kept at
`/tmp/EXHALE_L7b_control.x`, linked 19:01) and the MEASURED build is the same
tree with this item's five files changed (linked 19:10). This is a working
tree several items share, so it moves: between those two links one file that
is not this item's changed, `Cool_coeff.f90` (another item's new
density-resolved C / N / O tables and the `metal_cooling_above_ground_term`
call). That code is unreachable in a run without metals, so pair 1 isolates
this item exactly on every case and fixture below except the two that carry
`metals.inp`.

**Pair 2**, built for those two: one frozen copy of the tree
(`/tmp/L7b_frozen`), compiled once as it stands (`EXHALE_M.x`) and once with
this item's five files restored to their entry text (`EXHALE_C.x`, verified
to carry none of the new strings). The two differ in this item's five files
and in nothing else at all.

A third binary, an ATTRIBUTION PROBE, is built in section 5.3: the measured
objects with one object recompiled, so that one of the three changes can be
switched off alone.

## 4. Change 1: the absolute floor under the carrier balance row

### 4.1 The defect

READ. The carrier row's scale already carries a floor, `1e-20 n_ref sigrate`
(`carrier_residual`), and the elemental transport rows carry the same
constant against their own reservoir, `1e-20 rho X_base/(R0/v0)`
(`element_transport_residual`). Both are floors on the ROW'S TERMS, and they
answer "is this imbalance resolvable". Neither answers "is there a species
here", and that is the question the LHS 1140 b wind asks: where H2 is absent
the row's terms are not round-off at all -- the formation rate from H + H is
a real rate -- so the residual is the whole source, the scale is the same
source, and the measure is exactly 1.000 however little H2 the cell holds.

### 4.2 RED and GREEN

The state is the pilot's T-L7-5 seed, `molecular_scalar_gj1132_kzz1e9/HeH2.13`
converted from the certified atomic solution at the same He/H with
`EXHALE_MOLECULAR_SEED_X2=local`, loaded under `Restart intent: stationary
equilibrate` (section 9). MEASURED, the certification of the state as loaded:

| | control | measured |
|---|---|---|
| carrier balance H2, whole column | 1.000 at cell 1 | 1.000 at cell 1 |
| the gated measure and its cell | 1.000 at cell 340 | 9.953e-1 at cell 218 |
| cells reported as absent | -- (no such statement) | 163, first cell 194 |
| refusing entries | 5 | 5 |

MEASURED, the H2 of that seed at the two cells and around them:

| cell | r [R_p] | T [K] | n(H2) [cm^-3] | n(H2)/n_H |
|---|---|---|---|---|
| 194 | 1.134 | 3.70e3 | 0 | 0 |
| 218 | 1.204 | 4.94e3 | 3.50e3 | 2.52e-6 |
| 300 | 1.857 | 4.87e3 | 0 | 0 |
| 340 | 2.723 | 3.31e3 | 0 | 0 |
| 400 | 5.904 | 1.52e3 | 1.70e3 | 2.03e-3 |
| 450 | 12.73 | 8.62e2 | 2.15e4 | 3.05e-1 |

So the gate moves off a cell holding exactly no H2 and onto one holding
2.5e-6 of the hydrogen -- fourteen decades above the floor -- where the row
is a statement about a species the state has. **The floor removes the
refusal the pilot identified; it does not by itself certify the molecular
LHS 1140 b wind**, because the seed's H2 at 1.2 R_p is far from its own
balance. The pass table of section 7 is what that costs.

### 4.3 The hot Uranus fixture

MEASURED, `backup/regression/carrier_model_a_newton` reloaded on scratch
copies (`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40`, 4 threads), control
against measured:

| | control | measured |
|---|---|---|
| outer passes completed | 12 | 12 |
| hydrodynamic mass / momentum / energy row | 4.525e-11 / 1.566e-12 / 3.518e-09 | 5.972e-11 / 2.225e-12 / 3.261e-09 |
| carrier balance H2, gated | 6.996e-04 at cell 363 | 6.996e-04 at cell 363 |
| carrier balance H+, gated | 7.425e-04 at cell 323 | 7.425e-04 at cell 323 |
| cells reported as absent | -- | none (no line printed) |
| verdict | NOT CERTIFIED, 2 entries | NOT CERTIFIED, the same 2 entries |
| `Hydro_ioniz.txt` / `Ion_species.txt` against the control | -- | 4.85e-07 / 1.30e-06 worst relative |

Two things follow. The fixture does not certify on this tree in EITHER
build -- it refuses on the two carrier rows at 7e-4 against 1e-5 -- so the
pilot's reading that it certifies is not reproduced here and is not
something this change took away. And no cell of that hot Uranus carries a
carrier below the absolute floor, so the new gate is inactive there: its H2
is present everywhere its rows are judged, exactly as the pilot said. The
movement of the hydrodynamic rows is the other two changes (this fixture
carries `q_H2_base = 0.84`), and it stays three decades inside the
tolerances. The fixture carries no `metals.inp`, so pair 1 isolates this
item here.

### 4.4 The x_H2 pin: design, not a defect

READ, three facts settle it.

- `nmol_eq` is rebuilt from `f_sp` at the top of every sweep
  (`ionization_equilibrium`, the `nmol_eq(:,1) = f_sp(:,isp_H2)*rho*n0`
  block), so `x_h2_fix = 2 nmol_eq(j,1)/nh(j)` is the value the CARRIER
  TRANSPORT operator last wrote, not the value the file was loaded with.
- The H2 row the transport operator solves takes its source from
  `carrier_source`, which returns `fv(4)` of `mol_heh_rows`, the full H2
  balance including the formation channels. Those do not vanish at
  n(H2) = 0, so an empty cell has a nonzero source and zero is not a fixed
  point of the transported carrier.
- The pilot measured exactly that: n(H2) at cell 340 grows 0 -> 1e-29 ->
  1e-27 over five outer passes (L7 section 7.3, READ).

The pilot's reading that "zero is a fixed point of the loaded transported
carrier" is therefore refuted by its own measurement; what held the row at
1.000 was the missing absolute floor, and that is change 1. Nothing was
changed in the pin.

## 5. Change 2: the molecular energy bracket

### 5.1 The defect, and the validity of the extension

READ. `energy_semi_implicit` gave a cell with n(H2) > 0 a bracket ceiling of
5e4 K -- the top of the H2 rovibrational table -- and every other cell 1e7 K.
But `h2_rovibrational_energy_and_heat_capacity` does not stop at the table:
above `T_utab_hi` it returns `c_rv = dudlnT_utab(n_utab)/T_utab_hi` and
`u_rv = u_utab(n_utab) + c_rv (T - T_utab_hi)`, so du_rv/dT = c_rv exactly,
u_rv is monotone in T (the property the energy-to-temperature inverse relies
on), and the continuation carries no dissociation energy. What it is, above
the table, is an extrapolation of the BOUND ladder into a region where H2
does not survive; its weight in the energy of a cell is bounded by the
dissociation energy of that ladder, D0/k = 51966 K (READ, `caloric_eos`), so

```
x2 u_rv/(1.5 T)  <=  x2 (D0/k)/(1.5 T) ,   x2 = n(H2)/(n_tot + n_e) ,
```

which at the cells that were refused (T = 7.2e4 K, x2 = 4.0e-4; L7 section
5.1, READ) is 1.9e-4 of the internal energy. The bound falls as 1/T, so the
extension cannot become more important further out of the table than it is
at its edge.

### 5.2 RED and GREEN

`molecular_scalar_gj1132_wellmixed/HeH0.083`, cold march, 8 threads,
`EXHALE_MAXSTEPS=1200`.

| | control | measured |
|---|---|---|
| energy update failures | first at step 931 (cell 211, T 5.08e4 K), then 932, ... | 0 |
| the run | ATTEMPTED STEP EXHAUSTED at step 983 after 9 attempts (cell 228, T 7.07e4 K, `no bracket below the equation-of-state ceiling`, dt down to dt/2**8), `STOP 2` | 1200 steps, `final: count=1200 du= 9.5987E+00`, log Mdot 9.90, exit 0, 5 m 57 s |
| largest T of the state written | -- (the run wrote its initial condition) | 3.941e4 K at r = 1.683, and NO cell above 5e4 K |

The last line is the physical content. The probe of section 5.3, which
carries the other two changes and keeps the 5e4 K ceiling, reaches the same
four cells at 7.07e4 K on this path, so the state the ceiling refuses does
arise; the measured build carries the march through it and ends with no cell
above 3.94e4 K. The excursion is transient, and the bracket was refusing the
step that would have taken the column out of it.

### 5.3 The change is what unblocks it, measured

The control differs from the measured build in three places, so the march
was run a third time with an ATTRIBUTION PROBE: the measured build's own
objects with one object swapped, `energy_semi_implicit` compiled from a copy
that restores the 5e4 K molecular ceiling and changes nothing else. Changes
1 and 3 are present in it; change 2 is not.

| build | first failure | the stop |
|---|---|---|
| control | step 931, cell 211, T 5.08e4 K, `\|R\|/scale` 1.6048e-02 | ATTEMPTED STEP EXHAUSTED at step 983, cell 228, T 7.07e4 K |
| probe (changes 1 and 3, no change 2) | step 931, cell 211, T 5.08e4 K, `\|R\|/scale` 1.6196e-02 | ATTEMPTED STEP EXHAUSTED at step 982, cell 228, T 7.07e4 K, 624 failures |
| measured | none | 1200 steps, exit 0 |

Same cells, same temperatures, same statement; the one-step difference in
where the march exhausts and the 0.9 percent in the residual are the ghost
partition of change 3, which moves the state in the twelfth digit from step 2
onward. **The bracket is what stopped the run, and removing it is what lets
it through.**

## 6. Change 3: the base handoff's H2 partition

### 6.1 The defect

READ. `ieq_cell%x_h2_fix = base_h2_nuclei_fraction()` pinned x(4) = 2 n(H2)/n_H
over ALL the hydrogen nuclei of the two lower ghosts, leaving 1 - x2 for H+
and for the hydrogen in H2+, H3+ and HeH+ together. At the Koskinen 2022
handoff, x2 = 0.9999831600354949, that room is 1.68e-5, and nothing compared
it with the ionization the ghost's own balance produces.

### 6.2 RED and GREEN

`models/.L7/molecular_kzz1e9_HeH1.60`, cold march, 8 threads,
`EXHALE_MAXSTEPS=3000`.

| | control | measured |
|---|---|---|
| the run | `ERROR STOP ioniz_eq: persistent non-root chemical equilibrium`, ghost cell -1, step 1240, 1000 consecutive non-root sweeps | 3000 steps, `final: count=3000 du= 9.6073E+00`, log Mdot 9.71, exit 0, 9 m 20 s |
| the candidate the solve returned at the stop | x(1) = 4.3512e-4, x(4) = 0.99956, x(1) + x(4) = 1.00000 | -- |
| element violation reported | 4.253e-4, which is the amount by which the pinned row had to be broken | -- |
| ghost -1 of the final state: x(H+) | -- | 1.8546e-2 |
| ghost -1: x(4) = 2 n(H2)/n_H | -- | 0.981350552 |
| ghost -1: x(4)/(1 - x_ion) | -- | 0.999983241 against the handoff's x2 = 0.999983160 |
| ghost -1: x(1) + x(4) | 1.000307 at the pilot's stop (READ) | 0.999896741 |

The last two lines are the statement of the change. The prescription is met
as the handoff states it -- the molecular partition of the non-ionized
hydrogen agrees with x2 to eight digits -- and the sum of the pinned
molecular hydrogen and the ionized hydrogen is now inside the hydrogen the
cell has. The ghost's ionization at the end is 1.9e-2, which is 1100 times
the room the old prescription left it; under the old form the system was
infeasible by that factor and the sweep could only report a streak of
non-roots.

The refusal that was added fires on no cell of any run made here, which is
what the arithmetic says it should do: with x2 <= 1 the condition
x2 (1 - x_ion) + x_ion <= 1 is (1 - x2)(1 - x_ion) >= 0 and holds
identically. It is there so that an impossible handoff is named -- cell,
imposed value, room, and the ionization that consumed it -- instead of
reaching the solver and being reported as a streak.

## 7. Item 4: T-L7-5 again, on the HeH2.13 seeds

`molecular_scalar_gj1132_kzz1e9/HeH2.13` seeded from the certified atomic
solution at the same He/H, `Restart intent: stationary equilibrate`,
`EXHALE_PTC_DTAU0=1.0`, 8 threads. The `local` partition, measured:

| pass | control | measured |
|---|---|---|
| 1 | info = 0; 1.00E+00 at cell 340; mass 7.20e-07, momentum 5.35e-11, energy 6.42e-06 | info = 0; 9.95e-1 at cell 218; mass 7.21e-07, momentum 5.34e-11, energy 6.42e-06 |
| 2 | info = 0; 1.00E+00 at cell 340; mass 7.20e-07, energy 6.43e-06 | info = 0; 9.96e-1 at cell 218; mass 7.20e-07, energy 6.43e-06 |
| 3 | info = 0; 1.00E+00 at cell 340; mass 7.21e-07, energy 6.43e-06 | info = 0; 9.96e-1 at cell 218; mass 7.21e-07, energy 6.43e-06 |
| 4 | info = 0; 1.00E+00 at cell 340; mass 7.22e-07, energy 6.43e-06 | info = 0; 9.96e-1 at cell 218; mass 5.88e-08, energy 5.20e-07 |
| verdict | REFUSED at pass 4, no joint progress; NOT CERTIFIED, 4 entries: mass 7.221e-07 (tol 5.0e-09), energy 6.431e-06 (tol 1.0e-06), carrier H2 1.000 at cell 340, elemental He/H 4.854e-05 | REFUSED at pass 4, no joint progress; NOT CERTIFIED, 3 entries: mass 5.724e-08 (tol 4.7e-09), carrier H2 9.959e-1 at cell 218, elemental He/H 4.854e-05 |

READ, the pilot's own run of the same route with its own binary
(`.L7/route213_local`) agrees with the control to three digits at every pass
(1.00E+00 at cell 340 throughout, mass 7.21e-07 to 7.22e-07, energy 6.42e-06
to 6.43e-06, REFUSED at pass 4, 4 refusing entries), so the control
reproduces the pilot and the comparison above is between two builds of one
tree.

So the alternation behaves the same way and stops at the same place, and the
change shows up in two places: the carrier row is now judged where H2 exists
(cell 218, 9.96e-1) instead of where it does not (cell 340, exactly 1.000),
and the hydrodynamic rows of the last pass are an order of magnitude smaller
(mass 5.88e-08 against 7.22e-07, energy 5.20e-07 against 6.43e-06). The
energy row of the measured build's final state is WITHIN its tolerance and
the control's is not, so the state ends with three refusing entries instead
of four. The hydrodynamic gain belongs to change 3: this run's ghost carries
the handoff partition, and the control's is the over-prescribed one.

MEASURED, the H2 of the state the route wrote:

| cell | r [R_p] | n(H2) [cm^-3] | n(H2)/n_H |
|---|---|---|---|
| 194 | 1.134 | 5.15e-32 | 1.43e-41 |
| 218 | 1.204 | 3.28e+03 | 2.52e-06 |
| 250 | 1.358 | 9.64e-27 | 2.41e-35 |
| 340 | 2.723 | 1.84e-27 | 1.26e-34 |
| 400 | 5.904 | 1.41e+03 | 2.03e-03 |

163 cells stand between 34 and 41 decades below their element and are
reported without gating; the row that refuses is the H2 front cell at
1.204 R_p, which holds 2.5e-6 of the hydrogen. **What blocks the molecular
LHS 1140 b solution is now a statement about H2 that is there**, and that is
the next item's subject rather than this one's.

### 7.1 The handoff partition

MEASURED with the measured build, seeded at the handoff's own
x2 = 0.9999832: the first JFNK runs 62 iterations without descending, `||R||`
between 1.769 and 2.894 against a start of 1.989, the worst row the energy of
the cells at r = 1.13 to 3.8, and no outer pass completes (the run was
stopped there). That is what the pilot measured for this partition with its
own binary (L7 section 7.2, READ: no descent, `||R||` 1.93 to 2.73, zero
outer passes). The partition, not the floor, is what that seed fails on:
taking every neutral hydrogen atom out to 30 R_p removes the H I
photoionization heating the wind is held up by, and no change made here
touches that.

## 8. The seven `mol_*` regression cases: the movement

Each case was run on a SCRATCH COPY of its directory (never in
`backup/regression/<case>/`, which another item was using), single-threaded
with its own `maxsteps`, once with the control binary and once with the
measured one, and compared with `backup/regression/compare_within_tolerance.py`.
No golden was refreshed. MEASURED, the worst relative difference over the
four compared files:

| case | measured against control | the file and cell that carries it | measured against the golden |
|---|---|---|---|
| `mol_base_handoff` | 2.539e-05 | `Hydro_ioniz_adv` row 3 col 10 | 2.549e-05, WITHIN |
| `mol_lyman_werner` | 2.182e-05 | `Hydro_ioniz_adv` row 2 col 10 | 2.179e-05, WITHIN |
| `mol_diffusion` | 2.545e-05 | `Hydro_ioniz_adv` row 3 col 10 | 4.200e-04, WITHIN |
| `mol_carrier` | 2.076e-05 | `Hydro_ioniz_adv` row 3 col 10 | 2.072e-05, WITHIN |
| `mol_sec_ion` | 1.161e-04 | `Hydro_ioniz` row 177 col 3 (v) | 8.866e-03, EXCEEDS |
| `mol_metals` | 3.022e-06 | `Hydro_ioniz_adv` row 213 col 10 | 7.315e-01, EXCEEDS |
| `mol_ir_bands` | 2.683e-04 | `Hydro_ioniz_adv` row 37 col 10 | 9.117e-01, EXCEEDS |

Every case moves by less than the 1e-3 the harness calls identical, and in
six of the seven the cell that carries the movement is a boundary row of the
base -- which is where the one change that can reach an accepted molecular
state acts. `mol_sec_ion` is the exception: its worst cell is 177, in the
wind, and 1.2e-04 there is the staged secondary ionization following a base
that moved. In `mol_base_handoff`, MEASURED at ghost 0,

| | x(H+) | x(4) = 2 n(H2)/n_H | x(4)/(1 - x_ion) |
|---|---|---|---|
| control | 1.492834e-08 | 0.985447826087 | 0.985447852482 |
| measured | 1.492834e-08 | 0.985447799692 | 0.985447826087 |

and `h2_mixing_ratio_base()` = 0.84 gives x2 = 0.985447826087 exactly. The
control imposes x2 on every nucleus; the measured build imposes it on the
non-ionized ones and reproduces x2 to twelve digits as the partition OF
those. The ghost of this hot Uranus is ionized to 2.7e-8, so the change is
2.7e-8 of the H2 and 1.8e-6 of the ghost's H I, and that is the whole
movement of the metal-free cases. `mol_metals` moves a thousand times less
again because its metal electrons hold the ghost's ionization at 2.3e-11.

**Two of the numbers above needed a second control.** `mol_metals` and
`mol_ir_bands` are the two cases that carry `metals.inp`, and between the
control build and the measured build another item changed the C / N / O
metal cooling of `Cool_coeff.f90` (new density-resolved tables and the
`metal_cooling_above_ground_term` call in `cool_CI_ne` and its siblings).
That change is unreachable without metals, so the five metal-free rows above
are this item's own movement; for the two metal rows a SECOND PAIR of
binaries was built from one frozen copy of the tree, differing in this
item's five files and nothing else, and the numbers in the table are that
pair's. Measured against the FIRST control, which lacks the new cooling
tables, the same two cases move by 3.755e-03 and 1.237e-02 -- the metal
cooling, not this item.

The golden column is the tree's distance from a reference taken before this
stage's work, not this item's: `mol_lyman_werner` moves 1.3e-07 under this
item and sits 1.3e-07 from its golden, while `mol_metals` moves 3.0e-06
under this item and sits 7.3e-01 from its golden. No golden was refreshed.

## 9. How to reproduce

`EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00`. The private build
is `make OBJDIR=build_L7b EXE=EXHALE_L7b.x`; every run below is on a scratch
copy under `$EX/LHS1140b/models/.L7b/`.

```bash
# the two cold marches (changes 2 and 3)
for c in wellmixed HeH1.60; do mkdir -p .L7b/march_$c/output; done
\cp -f molecular_scalar_gj1132_wellmixed/HeH0.083/{input.inp,base.inp} .L7b/march_wellmixed/
\cp -f .L7/molecular_kzz1e9_HeH1.60/{input.inp,base.inp}               .L7b/march_HeH1.60/
( cd .L7b/march_wellmixed && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=1200 $EX/EXHALE_L7b.x > march.log 2>&1 )
( cd .L7b/march_HeH1.60   && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=3000 $EX/EXHALE_L7b.x > march.log 2>&1 )

# T-L7-5 (change 1 and item 4): the route from the pilot's HeH2.13 seeds
for p in local handoff; do
  d=.L7b/route213_$p; mkdir -p $d/output
  \cp -f .L7/route213_local/{input.inp,base.inp} $d/
  \cp -f .L7/seed213_$p/output/{Hydro_ioniz_IC.txt,Ion_species_IC.txt} $d/output/
  ( cd $d && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 $EX/EXHALE_L7b.x > run.log 2>&1 )
done

# the hot Uranus carrier fixture, on a scratch copy
d=.L7b/fixture; mkdir -p $d/output
\cp -f $EX/backup/regression/carrier_model_a_newton/{input.inp,base.inp} $d/
\cp -f $EX/backup/regression/carrier_model_a_newton/IC/*_IC.txt $d/output/
( cd $d && EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=4 \
      $EX/EXHALE_L7b.x > run.log 2>&1 )

# the seven mol_* regression cases, on scratch copies, never in place
#   <case>/{input.inp,base.inp,metals.inp,maxsteps} copied to the scratch
#   directory, then  env EXHALE_MAXSTEPS=$(cat maxsteps) OMP_NUM_THREADS=1
#   and backup/regression/compare_within_tolerance.py against the control
#   run and against backup/regression/golden/<case>/.
```

The test suites were run with `EXHALE_OBJDIR=$EX/build_L7b`,
`EXHALE_TEST_OUT` and `EXHALE_TEST_OBJDIR` pointing into
`$EX/build_L7b/tests/<suite>`, and `EXHALE_EXE`/`EXHALE_RESID_EXE` pointing
at the private binary where the suite runs one.

## 10. The test suites

Every suite's driver links the whole production object set, so every suite
links the four modules changed here; all of them were run against
`build_L7b`, and the ones that run a binary were run against
`EXHALE_L7b.x`. MEASURED:

| suite | PASS | FAIL |
|---|---|---|
| certification | 84 | 0 |
| energy_update | 53 | 0 |
| carrier_reference_scales | 14 | 0 |
| carrier_returned_state_acceptance | 36 | 0 |
| carrier_constraint_attribution | 10 | 0 |
| carrier_retry | 142 | 0 |
| steady_species_rows | 195 | 0 |
| ionization_imposed_fractions | 31 | 0 |
| acceptance_classes | 22 | 0 |
| coupled_source_step | 30 | 0 |
| attempted_step | 70 | 0 |
| adv_static_limit | 53 | 0 |
| krylov_and_dogleg | 334 | 0 |
| element_operator | 28 | 0 |
| species_masses | 9 | 0 |
| constrained_network_layout | 8 | 0 |
| species_face_flux | 1 | 0 |
| molecular_seed | 17 | 0 |
| run_mode | 29 | 1 |
| residual_determinism | 5 | 1 |
| grid_and_gates | 195 | 4 |

The three suites with a FAIL were each re-run with the CONTROL binary, and
every one of those failures is on the tree before this item:

| failing assertion | control | measured |
|---|---|---|
| `run_mode` / `a_phys_header_without_its_clock_is_refused` | accepted (FAIL) | accepted (FAIL) |
| `residual_determinism` / `closure_spread_within_the_row_tolerance_atomic_elem_newton` | 1.713e+04 (FAIL) | 1.715e+04 (FAIL) |
| `grid_and_gates` / `stationary_evaluate_certification_unchanged`, `..._exit_status`, `..._rebuilt_within_sweep_budget`, `outer_iteration_ending_is_the_stagnation_one` | FAIL, rebuilt 6.427e-01 | FAIL, rebuilt 6.510e-01 |

None of the three is reachable by this item's changes: the first is a
restart-header check in `load_IC`, the second an ATOMIC element reload
(neither the carrier gate nor the handoff partition nor a molecular bracket
exists in it), and the third a WASP-121 b state with metals whose stored
certification no longer reproduces under the metal cooling this tree now
carries. They are reported here and not fixed: each belongs to another
item's file.

## 11. Noticed outside this item's scope, reported and not fixed

- `src/tests/residual_determinism/run.sh`, `run_mode/run.sh`,
  `adv_static_limit/run.sh`, `grid_and_gates/run.sh`, `molecular_seed/run.sh`
  and `fuv_band_ledger/run.sh` default to `$ROOT/EXHALE.x`, the SHARED
  binary. A worker running them while another worker rebuilds that binary
  measures someone else's build: `molecular_seed` reported 11 failures that
  way and 0 against the private binary. The default is worth changing to a
  refusal when `EXHALE_EXE` is not given.
- `src/tests/carrier_retry/run.sh` and the other suites that read
  `EXHALE_TEST_OBJDIR` (not `EXHALE_TEST_OUT`) write their driver into the
  SHARED `build/tests/<suite>/` unless that variable is set, which two
  workers cannot do at once.
- The three pre-existing suite failures of section 10.

## 12. What this item left in the tree

- the five source files of section 2;
- this report;
- `LHS1140b/models/.L7b/`, the scratch tree of the item: the two cold
  marches with each build, the attribution probe, the two T-L7-5 routes with
  each build, the carrier fixture with each build, and `reg/` with the four
  regression sets (`control`, `modified`, and `fcontrol`/`fmeasured` for the
  two metal cases). 84 MB.

Nothing else in the repository. The private build directory and its binary
were deleted, no golden was refreshed, and no regression case was run in
`backup/regression/<case>/` -- every one ran on a scratch copy.

One thing WAS written outside that, and is recorded here rather than
repaired: the first pass of the test suites was launched with
`EXHALE_TEST_OUT` but not `EXHALE_TEST_OBJDIR`, and the suites that read the
second variable (`attempted_step`, `carrier_constraint_attribution`,
`carrier_reference_scales`, `carrier_retry`,
`carrier_returned_state_acceptance`, `constrained_network_layout`,
`ionization_imposed_fractions`, `species_masses`) therefore built their
drivers into the shared `build/tests/<suite>/` at 19:23 to 19:24 before that
was noticed. Nothing of another item's was deleted; the drivers are
rebuilt by whoever runs those suites next. Every suite result quoted in
section 10 comes from the second pass, with both variables set.

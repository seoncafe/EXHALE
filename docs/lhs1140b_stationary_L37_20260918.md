# L37: the round trip of the base mass row on a molecular state

Item L37 of `docs/PLAN_20260917.md`, opened after L34b reported that the
no-step evaluate route refuses the certified molecular reference solution on
one entry, the hydrodynamic mass row of cell 1. 2026-09-18.

Binary for every measurement: `LHS1140b/models/EXHALE_3146d11b.x`, md5
`3146d11b4090306dcea75bb9718edd22`, the one the L34 states were made on.
Every run at `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, on scratch copies.
No tolerance, no anchor and no key default was changed by this item. One
source file was edited, `src/tests/physics_probe/molecular_energy_recipients.f90`
(section 6). Nothing under `backup/regression/` and no catalog case directory
was written to.

Every number is labelled MEASURED (this item ran it) or READ (from a file, a
log or a document).

## 1. Verdict

**The band.** On the six certified molecular states in hand the cell-1
hydrodynamic mass row comes back from a single no-step re-entry at **2.39 to
10.09 times** the value the run printed for the same cell, and at **16.0 to
31.2 times that cell's own rounding floor**, where the tolerance is ten
floors. On three certified ATOMIC states of the same planet, the same grid
and the same Kzz the ratio is **0.68 to 1.35** and the row stands at **0.21
to 1.52 floors**. The atomic figure is L18's 0.78 to 1.71 band met again on
three states; the molecular one is a decade above it.

**It is a transient of the re-entry, not a property of the state.** Feeding
the evaluation its own written state again, and again, the row settles: four
of the six molecular states reach a reading that repeats exactly by the third
or fourth re-entry, and at that reading cell 1 stands at **0.08 to 7.5
floors, inside the existing anchor on every one of the six**. The atomic
states are already at that reading on the FIRST re-entry (three of three
reproduce their own output exactly at cells 3 to 10, two of three at every
cell). The settled molecular reading also agrees with what the run itself
printed: 5.5734e-09 against the run's 5.213e-09 on the reference solution,
a factor 1.07.

**The mechanism is the base ghost cells and the caloric equation of state.**
The physical cells round-trip exactly: `rho`, `v` and `p` of cells 1 and 2
come back to all 17 digits and the composition of cell 1 to 1.7e-16. What
does not come back is the pair of inner GHOST rows, whose density and
composition the loader re-derives and projects; they move by parts in 1e-8
(density, cell 0) to 9e-6 (a trace level, cell -1). In the molecular base
that ghost is **99.998 per cent H2**, so the caloric equation of state makes
the base face flux a function of its H2 count, and MEASURED the amplification
into the cell-1 mass row is **0.15 to 0.33 of the row's own scale per
relative unit** of the ghost's H2 (or of its density, which the composition is
projected onto). In an atomic base the same perturbation, taken up to 1e-6,
moves the row by **nothing at the printed precision**, because there
`p = (gamma_ad - 1) rho e` carries no composition.

**The composition refresh of the physical cells is not the carrier.** No key
skips the one equilibrium sweep (READ, `stationary_state_of_the_loaded_restart`
and the `EXHALE_RESIDUAL` block of `src/EXHALE_main.f90`: both call
`ioniz_eq` unconditionally; `EXHALE_RELOAD_EQ=0` removes only the EXTRA
`equilibrate_loaded_composition` of the residual diagnostic, and the second
word `equilibrate` adds more sweeps, not fewer), so the sweep could not be
switched off. What can be measured instead is its size: the sweep moves the
physical particle count by 1.295e-13 (READ, the route's own line on the
reference solution), which at the amplification measured above is 4e-14 of
the row, five decades below the 6.9e-09 the re-entry moves it.

**The proposed anchor, and whether the rule settles it.** Following the rule
of `docs/certification_tolerance_anchoring_20260910.md` anchor (6) --
`c_round` is a round number above the largest measured ratio of the step the
row takes to its signal estimate -- the re-entry perturbation gives a largest
ratio of 31.2, and with P16's own margin of 3.7 that is 116, i.e. a round
**100** where the code carries 10. **The rule does not make the choice
unambiguous**, because it does not say which perturbation the ratio is
measured under: P16 measured one ulp of the density of every physical cell,
which the arithmetic forces, and got 1.59 to 2.68 (READ); the round trip is a
perturbation the arithmetic does NOT force, and two further re-entries of the
same state remove it. Raising `c_round` to 100 would widen a gate to cover a
restart that does not reproduce its own base ghost. **This item proposes no
change**; section 5 states the three ways out for the user.

## 2. The states, and the three readings

The certified molecular states in hand, from the two catalog trees and the
L22 record. Their claims are READ from the `# coupling:` line of each
`Hydro_ioniz_IC.txt`.

| name here | state | claim | binary it was made on |
|---|---|---|---|
| `heh213_new` | `molecular_scalar_gj1132_kzz1e9/HeH2.13/output` | `certified=T certified_in_wind` | `3146d11b` (L34b) |
| `heh213_pre` | the same case, `output_pre_L34` | `certified=T certified_in_wind` | `c2e9c999` (pre-L34) |
| `heh97_new` | `molecular_scalar_gj1132_kzz1e9/HeH9.7/output` | `certified=T certified_in_wind` | `3146d11b` |
| `heh97_pre` | the same case, `output_pre_L34` | `certified=T certified_in_wind` | `c2e9c999` |
| `heh055_new` | `molecular_scalar_gj1132_kzz1e9/HeH0.55/output`, copied 2026-09-18 06:2x | `certified=T certified_in_wind` | `c2e9c999` |
| `i3_alt` | `.L22/i3_alt/output` | `certified=T certified_in_wind` | the L22 binary |
| `at213`, `at97`, `at055` | `atomic_scalar_gj1132_kzz1e9/{HeH2.13, HeH9.7, HeH0.55}/output` | `certified=T certified_in_wind` | controls |

Two entries of the brief's list are not what they were taken for, and are
reported rather than used:

- **`heh055_new` was overwritten under this item.** The copy taken for this
  work carries `certified=T cert_reason=certified_in_wind` and reads the same
  ten rows as `output_pre_L34` to every printed digit; at 06:42 on 2026-09-18
  the file in the tree became `certified=F cert_reason=mapped_seed`, written
  by the L34 step (c) run of that case, which is still going. The
  measurement below is therefore of
  the archived certified state, under the name `heh055_new`, and the case
  directory was not touched.
- **No molecular regression fixture is a certified stationary state.** READ,
  `backup/regression/*/NOTE.txt` and the two READMEs: the eight `mol_*`
  DEFAULT_CASES are 12000-step marching snapshots pinned by their own
  `maxsteps`, and `carrier_elem_newton` and `carrier_model_a_newton` are
  reload entry points, the first of which the anchoring memo records at
  `||R||` 2.5e-1 with four entries refusing it (READ). None was re-run.

The three readings the brief names, on the same state and with no step taken
anywhere:

| tag | route | what it is |
|---|---|---|
| **E1** | `Restart intent: stationary evaluate` | the accepted route, the state as the file carries it |
| **R1** | `EXHALE_RESIDUAL=1` | the route L34b used for its third reading. Its `reload_equilibrate` defaults ON, so it puts the loaded composition on its own fixed point BEFORE measuring |
| **R0** | `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0` | the same without that extra equilibration |
| **E2, E3, E4** | `stationary evaluate` on the state the previous evaluation wrote | the re-entry ladder, which is how a one-shot reading is told from a settled one |

MEASURED, and it identifies L34b's third reading: on `heh213_new`, **R0 is E1
to every printed digit** (1.2475e-08) and **R1 is E2** (7.0724e-09, which is
the 7.072e-09 L34b reported). The third reading is not a second opinion about
the state; it is one settling pass further along the same ladder. The same
holds on `heh97_new` (R0 = E1 = 2.3398e-08, R1 = E2 = 2.3132e-09) and on
`at213` (R0 = R1 = E1, the atomic state having nothing to settle). On
`i3_alt` R0 reads 1.6332e-08 against E1's 1.4863e-08, a tenth, the two routes
differing in what they rebuild before the assembly.

## 3. The band, MEASURED

All rows at the cell the run's own log names for the mass row, `|R_1|/s_1`,
from `EXHALE_MASS_FLOOR_SCAN=1` (which reproduces the given state's mass rows
to 0.00E+00 on every run here, so the scan is not moving what it measures).
"in-run" is READ from the certification block of each state's own `run.log`;
"floor" and "tol" are the cell's own, `tol = 10 x floor`.

| state | cell | in-run | E1 (one re-entry) | E4 (settled) | floor | E1/in-run | E1/floor | E4/floor |
|---|---|---|---|---|---|---|---|---|
| `heh213_new` | 1 | 5.213e-09 | 1.2475e-08 | 5.5734e-09 | 7.8046e-10 | **2.39** | 15.98 | 7.14 |
| `heh213_pre` | 1 | 3.616e-09 | 1.7031e-08 | 3.0823e-09 | 7.6854e-10 | **4.71** | 22.16 | 4.01 |
| `heh97_new` | 1 | 2.319e-09 | 2.3398e-08 | 1.7066e-09 | 7.4894e-10 | **10.09** | 31.24 | 2.28 |
| `heh97_pre` | 1 | 2.330e-09 | 1.6456e-08 | 3.7884e-09 | 7.4264e-10 | **7.06** | 22.16 | 5.10 |
| `heh055_new` | 1 | 1.065e-08 | 2.5454e-08 | 5.0843e-09 | 1.2867e-09 | **2.39** | 19.78 | 3.95 |
| `i3_alt` | 1 | 3.249e-09 | 1.4863e-08 | 6.3734e-11 | 7.8046e-10 | **4.57** | 19.04 | 0.08 |
| `at213` | 1 | 1.189e-09 | 1.2016e-09 | 1.1852e-09 | 7.9044e-10 | **1.01** | 1.52 | 1.50 |
| `at97` | 2 | 1.481e-09 | 1.9924e-09 | 1.9924e-09 | 7.9345e-10 | **1.35** | 2.51 | 2.51 |
| `at055` | 1 | 7.414e-10 | 5.0212e-10 | 5.0212e-10 | 6.7147e-10 | **0.68** | 0.75 | 0.75 |

The in-run reading of `heh213_new` is the one the JFNK handed-back
certification printed at cell 1; the accepted-state certification of the same
conserved state printed its maximum at cell 3 (1.206e-09) with the same
verdict cell 188 and the same 2.611e-11 there. **Two certifications inside
one run, of one conserved state, already differ by a factor 4.3 at the base
cells and agree to four digits at the verdict cell** (READ, `run.log` lines
2352 and 2438): the re-entry is not the only place this row moves, it is
where it is measured.

### The same rows on cells 1 to 10

E1 against E4 and both in units of the cell's own floor. Cells 3 to 10 are
reproduced EXACTLY (ratio 1.000) by the re-entry in six of the nine states;
what moves is cells 1 and 2.

| state | cell | r | E1 | E4 | E1/E4 | E1/floor | E4/floor |
|---|---|---|---|---|---|---|---|
| `heh213_new` | 1 | 1.00019 | 1.2475e-08 | 5.5734e-09 | 2.24 | 15.98 | 7.14 |
| | 2 | 1.00039 | 2.7968e-09 | 3.9796e-11 | 70.3 | 7.30 | 0.10 |
| | 3 to 10 | | | | **1.000 each** | 0.68 to 3.90 | the same |
| `heh97_new` | 1 | 1.00019 | 2.3398e-08 | 1.7066e-09 | 13.7 | 31.24 | 2.28 |
| | 2 | 1.00039 | 5.3769e-09 | 1.6208e-09 | 3.32 | 11.75 | 3.54 |
| | 4 to 7 | | | | 0.55 to 2.39 | 0.41 to 1.45 | 0.17 to 0.87 |
| `i3_alt` | 1 | 1.00019 | 1.4863e-08 | 6.3734e-11 | 233 | 19.04 | 0.08 |
| | 2 | 1.00039 | 3.4221e-09 | 2.8248e-10 | 12.1 | 8.94 | 0.74 |
| `at213` | 1 | 1.00019 | 1.2016e-09 | 1.1852e-09 | 1.01 | 1.52 | 1.50 |
| | 2 | 1.00039 | 6.2266e-11 | 5.4141e-10 | 0.115 | 0.09 | 0.75 |
| | 3 to 10 | | | | **1.000 each** | 0.12 to 2.03 | the same |
| `at97` | 1 to 10 | | | | **1.000 each** | 0.08 to 2.51 | the same |
| `at055` | 1 to 10 | | | | **1.000 each** | 0.05 to 1.43 | the same |

**The ratio of the evaluate reading to the accepted one, taken cell by
cell, cannot be formed for cells 2 to 10**: a run prints the mass row's maximum, its cell and
the verdict cell, and no other cell's value, so there is one in-run number per
state and it is the one in the table above. E1 over E4 is what a round trip
can be measured as cell by cell, and it is what the second table carries.

### The base cell each of these rows belongs to

READ from row 3 of each `Hydro_ioniz_IC.txt`; the Mach number is MEASURED by
the scan.

| state | v [cm/s] | rho [m_H/cm3] | T [K] | Mach | 2 n(H2)/n_H at cell 1 |
|---|---|---|---|---|---|
| `heh213_new` | -5.62 | 2.834e+13 | 778.8 | 3.09e-05 | 0.3354 |
| `heh213_pre` | -5.73 | 2.741e+13 | 808.3 | 3.10e-05 | 0.3553 |
| `heh97_new` | -4.07 | 4.857e+13 | 522.2 | 2.93e-05 | 0.0676 |
| `heh97_pre` | -4.30 | 4.477e+13 | 567.5 | 2.97e-05 | 0.0780 |
| `heh055_new` | -5.92 | 1.925e+13 | 948.1 | 2.75e-05 | -- |
| `i3_alt` | -5.62 | 2.834e+13 | 778.8 | 3.09e-05 | 0.3354 |
| `at213` | -0.60 | 8.658e+13 | 238.0 | 5.74e-06 | atomic |
| `at97` | -1.02 | 1.013e+14 | 245.5 | 1.07e-05 | atomic |
| `at055` | -0.23 | 6.168e+13 | 231.4 | 1.84e-06 | atomic |

**The Mach number is not what separates the two groups.** The molecular base
cells run at 2.7e-05 to 3.1e-05 and the atomic ones at 1.8e-06 to 1.1e-05,
so `eps/Mach` is if anything SMALLER in the molecular base, and the rounding
floors agree to within a factor 2 across all nine states (6.7e-10 to
1.3e-09). The difference is in what a base face flux is a function of.

## 4. The mechanism, MEASURED

### 4.1 What does not round-trip

The state the evaluation writes back, against the state it was handed
(`heh213_new`, largest relative difference over the 504 rows):

| row | rho | v | p | T |
|---|---|---|---|---|
| ghost -1 | 9.07e-16 | **2.90e-11** | **1.82e-11** | **9.89e-12** |
| ghost 0 | 1.36e-16 | **2.90e-11** | 9.92e-14 | 5.50e-14 |
| every physical cell | 0.00 | 0.00 | 0.00 | <= 3.2e-15 |

and in the composition file, ghost -1 moves by 9.139e-06 (the He 2^3S
column), ghost 0 by 4.835e-08 (H I), and cell 1 by 4.7e-14 (HeH+). The
conserved state of the column is returned bit for bit; the two ghost rows
are not.

### 4.2 The ladder

Each evaluation handed the state the previous one wrote (`|R_1|/s_1` at cell
1):

| state | E1 | E2 | E3 | E4 |
|---|---|---|---|---|
| `heh213_new` | 1.2475e-08 | 7.0724e-09 | **5.5734e-09** | **5.5734e-09** |
| `heh97_new` | 2.3398e-08 | 2.3132e-09 | **1.7066e-09** | **1.7066e-09** |
| `heh97_pre` | 1.6456e-08 | 4.2485e-09 | **3.7884e-09** | **3.7884e-09** |
| `i3_alt` | 1.4863e-08 | 1.1231e-09 | **6.3734e-11** | **6.3734e-11** |
| `heh213_pre` | 1.7031e-08 | 1.6628e-09 | 4.3212e-09 | 3.0823e-09 |
| `heh055_new` | 2.5454e-08 | 4.6998e-09 | 5.0843e-09 | 5.0843e-09 (cells 7 and 8 still moving) |
| `at213` | 1.2016e-09 | 1.1852e-09 | 1.2016e-09 | 1.1852e-09 (a period-2 cycle of 1.4 per cent) |
| `at97`, `at055` | -- | **identical to E1 at every one of the ten cells** | | |

Four of the six molecular states reach a reading that repeats; two do not
settle at cells the verdict does not take. Two of the three atomic states are
at their settled reading immediately, and the third alternates between two
values 1.4 per cent apart.

### 4.3 Which half of the file carries it

The fixed-point state and the archived state differ only in their ghost rows,
so the two halves can be exchanged. MEASURED, cell 1 of `heh213_new`:

| state given to the evaluation | cell-1 row |
|---|---|
| the archived state (control) | 1.2475e-08 |
| the settled state (control) | 5.5734e-09 |
| the settled state with the ARCHIVED ghost composition rows | **1.2475e-08** |
| the archived state with the SETTLED ghost composition rows | 7.0724e-09 |
| the settled state with the ARCHIVED ghost hydrodynamic rows | 7.0724e-09 |
| the archived state with the SETTLED ghost hydrodynamic rows | 1.2475e-08 |

Giving the settled state the archived ghost composition reproduces the
archived reading exactly. The ghost rows are the carrier, and both files'
ghost rows act, because the composition is projected onto the density the
`Hydro_ioniz` ghost row carries.

### 4.4 The sensitivity, molecular against atomic

One column of the INNER GHOST row is scaled by `1 + eps` in the file and the
evaluation re-run from the settled state. Cell-1 row, MEASURED:

| perturbed | eps | `heh213_new` (molecular) | `at213` (atomic) | `at97` (atomic) |
|---|---|---|---|---|
| ghost n(H2) (molecular) / n(H I) (atomic) | 0 | 5.5734e-09 | 1.1852e-09 | 1.8517e-10 |
| | 1e-08 | 6.9972e-09 | 9.5494e-10 | 1.8526e-10 |
| | 1e-07 | 3.5562e-08 | 9.5494e-10 | 1.8526e-10 |
| | 1e-06 | 3.3171e-07 | 9.5494e-10 | 1.8526e-10 |
| ghost rho | 1e-08 | 7.0724e-09 | no change | no change |
| | 1e-06 | -- | no change | no change |
| ghost p | 1e-12, 1e-10, 1e-08 | no change | no change | -- |
| ghost v | 1e-10, 1e-08 | no change | -- | -- |

**The molecular response is linear with slope 0.30 to 0.33** of the row's own
scale per relative unit of the ghost H2 (3.3171e-07 at 1e-06 and 3.5562e-08
at 1e-07, against 5.5734e-09 unperturbed), i.e. **a ghost H2 count off by
1e-08 is four rounding floors of the cell-1 mass row**, and the measured
ghost drift of section 4.1 is of that size. **The atomic response is zero at
the printed precision up to eps = 1e-06**, a hundred times the molecular
perturbation that matters. The ghost pressure and velocity columns of the
file are inert in both, the base boundary rebuilding them.

The reason is in the pressure map. Where `caloric_mixture_active` is false
the map is `(gamma_ad - 1) rho e` and carries no composition (READ,
`docs/lhs1140b_stationary_L9_20260916.md` section 4, MEASURED there as the
two pressure rules agreeing to the last bit on an atomic fiducial and to
5.9e-10 on a molecular one); the inner ghost of these molecular columns is
**2 n(H2)/n_H = 0.99998**, so its pressure, its sound speed and therefore the
base face flux of the HLLC solve are functions of its H2 count. The mass row
of cell 1 is the difference of that face flux and the next, divided by a
scale which is the larger of the two, and at a base Mach number of 3e-05 the
difference is 3e-05 of each: a part in 1e8 of the face state is a part in
1e3 of the difference. That is the amplification the anchor already calls
`eps/Mach`, met here with a perturbation four decades above one ulp.

## 5. The anchor this measurement would give, and the contract question

### 5.1 What the anchoring rule says, and where it stops

`docs/certification_tolerance_anchoring_20260910.md`, anchor (6), fixes the
mass-row tolerance as `tol_mass(j) = max(3e-12, min(1, c_round x floor(j)))`
with `floor(j) = eps max_faces[rho(|v|+c_s)A] / max_faces|A F|`, and sets
`c_round = 10` as "a round number 3.7 times above the largest ratio measured
here" (READ), the ratios being the step each mass row takes when one ulp is
added to the density of every physical cell: 1.5931, 1.8857, 2.6806 and
1.7733 on the four states P16 measured (READ).

Applying the same rule with the re-entry as the perturbation, the ratios
MEASURED here are 15.98, 19.04, 19.78, 22.16, 22.16 and 31.24 on the six
molecular states, and 0.75, 1.52 and 2.51 on the three atomic ones. The
largest is 31.24, and P16's own margin of 3.7 above the largest gives 116:

> **the proposal the rule would produce: `c_round = 100` for a base cell of a
> molecular column, a decade above the value the code carries.**

**The rule does not make that choice unambiguous.** It fixes the margin above
a measured ratio but not which perturbation the ratio is measured under, and
the two candidates are not the same kind of thing:

- one ulp of the density is a perturbation the arithmetic forces on any
  evaluation of the state, and no re-evaluation removes it (P16 measured the
  state reproducing its own mass rows to 0.00E+00, and so does every run of
  this item);
- the restart round trip is not forced: **two further re-entries of the same
  state remove it**, and at the settled reading every one of the six
  molecular states stands at 0.08 to 7.5 floors, inside the tolerance the
  code carries today.

Anchoring on the second would widen a gate to admit a restart that does not
reproduce its own base ghost, which is the move the standing rule forbids.
The number is therefore reported and NOT proposed for adoption.

### 5.2 The three ways out, for the user

1. **Make the restart carry its base ghost.** The state pair already writes
   the ghost rows; what the loader does not reproduce is the ghost
   composition and the density it is projected onto, and the evaluation
   rebuilds them. If the re-entry landed on the settled reading, the existing
   `c_round = 10` would pass every molecular state measured here with a
   margin of 1.3 to 125.
2. **Make the certification's input the route's own fixed point**, i.e.
   re-enter until the reading repeats (two extra no-step evaluations, 7 s
   each on one thread here) and certify that. This is a contract change, and
   it certifies a state one ghost rebuild away from the one the run handed
   back.
3. **Re-derive `c_round` on the round-trip perturbation**, which gives 100
   for a molecular base. Not proposed.

### 5.3 The contract question the brief asks

**Should the evaluate route's composition refresh be the certification's
input at all for a molecular base?** On the measurement, the refresh of the
PHYSICAL cells is innocent: their composition returns to 1.7e-16, the sweep
moves the physical particle count by 1.295e-13, and at the amplification of
section 4.4 that is 4e-14 of the row, five decades below what the re-entry
moves. What the evaluation does that the run did not is rebuild the base
GHOST state from a file that does not carry it to the precision this row
needs. So the contract question is not about the sweep; it is:

> **what must a restart pair carry so that the base face flux of a molecular
> column is the same flux after a write and a read?**

That is the user's and the advisor's to answer. It is also the reason the
question is worth answering rather than tolerating: the same ghost carries
the base mass flux the mass-loss rate is read from, and section 4.4 says a
part in 1e8 of it is four rounding floors of the row that gates the state.

## 6. The probe can now be pointed at a state

`src/tests/physics_probe/molecular_energy_recipients.f90` held its three
representative cells as `real*8, parameter` arrays transcribed by hand, so
every re-solve of the molecular fiducial made the L31 uncertainty bracket
stale. It now reads them.

- `EXHALE_PROBE_STATE` names a directory; the pair read is
  `{Hydro_ioniz_IC.txt, Ion_species_IC.txt}` (the state a run hands back; the
  pair without `_IC` is the product written FROM a measurement of it and is a
  different state, section 4.1). Unset, the compiled values stand and nothing
  is read, so the suite runs with no state in the tree.
- The species columns are taken BY NAME from the file's own `# columns`
  header. The electron density and the particle count are formed as the
  code's mass policy forms them: a column whose name ends in III carries two
  charges and one ending in II carries one, H2p, H3p and HeHp carry one each,
  and HeITR is a level inside He I and enters neither sum. Both reproduce the
  transcribed values exactly.
- The three cells are a rule: **molecular depth** = physical cell 1;
  **H2 front** = the first cell where `2 n(H2)/n_H` has fallen to half its
  value at cell 1 (L31 section 5.4's own definition); **dilute upper column**
  = the cell nearest `r_dilute_upper_column = 1.8422600916216840` R_p. L31
  states no rule drawn from the profile for the third, and none exists: at
  that radius `2 n(H2)/n_H` is flat and rising outward, so no threshold picks
  the cell. The radius IS the choice, and it is recorded as the radius of the
  cell the record was first written on.

Reproduction, MEASURED:

| pointed at | cells | T of the three cells [K] | distance from the compiled values |
|---|---|---|---|
| unset | 1, 178, 299 | 808.32, 2128.37, 5083.44 | -- (the compiled values) |
| `HeH2.13/output_pre_L34` | 1, 178, 299 | 808.319, 2128.367, 5083.438 | **0.00000E+00** over all 33 numbers |
| `HeH2.13/output` (the L34b solution) | 1, **177**, 299 | 778.79, 2116.61, 5061.59 | 3.41525E-01 |

The archived reading reproduces L31 section 5.4 digit for digit (r 1.000193,
1.101040, 1.842260; `2 n(H2)/n_H` 0.35530, 0.17366, 3.6e-06; n(H I) 1.868e+12,
1.564e+10, 8.162e+07; n_e 9.000e+06, 2.333e+07, 6.901e+06). The new reading
reproduces L34b block G: the three temperatures 778.79 / 2116.61 / 5061.59 and
the helium third body at the depth **+42.05 / -21.07 per cent** of a nominal
heat of 1.04863e-06 (L34b: +42.0 / -21.1), the R5/R16 recipient at the front
**+4.37 per cent** of 8.60218e-09 (L34b: +4.4). The front moves one cell out
of 178 into 177, which is the movement L34b's comparison table reports.

Suite: `src/tests/physics_probe/run.sh` with the variable unset, built into a
scratch object directory, **1604 PASS / 0 FAIL, `physics_probe: PASSED`**,
this driver contributing the 81 assertions L34b also counted (READ). A
control run of the entry text was NOT taken: the file in the tree stands
ahead of the last commit by other items' work, so `HEAD` is not this
driver's entry text and no other copy of it exists. What is shown instead is
that the reader reads: pointed at the archived fiducial the distance from the
compiled values is 0.00000E+00 and pointed at the new solution it is
3.41525E-01, so the values are not the compiled ones coming back; pointed at
an atomic state the driver stops with `FAIL probe_state_has_no_column_H2`
and pointed at a directory with no state pair with
`FAIL probe_state_file_unreadable`, rather than reading a zero (MEASURED,
both exit 1). The driver alone, pointed at each of the two molecular states,
81 PASS / 0 FAIL.

## 7. How to repeat any of it

- the three readings of one state: copy `input.inp`, `base.inp` and the
  `_IC` pair into a scratch directory, set `Load IC? True` and
  `Restart intent: stationary evaluate`, and run the binary with
  `EXHALE_MASS_FLOOR_SCAN=1`, whose block prints `|R_1|/s_1`, both floor
  estimates and the Mach number of the first thirty cells. `EXHALE_RESIDUAL=1`
  and `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0` are the other two.
- the ladder: hand the next evaluation `output/Hydro_ioniz.txt` and
  `output/Ion_species.txt` of the previous one, renamed `_IC`.
- the ghost sensitivity: scale one column of row 2 (the inner ghost) of
  either file by `1 + eps` and evaluate again.
- the probe on a state: `EXHALE_PROBE_STATE=<dir>` with
  `src/tests/physics_probe/run.sh` or the driver alone.

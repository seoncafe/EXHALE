# LHS 1140 b, item L36d: the ionization stages of both elements carried in the run

Third increment of item L36 of `docs/PLAN_20260917.md` section 10, which
carries stages B to F of `docs/lhs1140b_stationary_L12a_design_20260913.md`
as `docs/lhs1140b_stationary_L12b_derivation_20260916.md` sequenced them. It
follows `docs/lhs1140b_stationary_L36_20260918.md` (the stage flux as a
module with no caller) and `docs/lhs1140b_stationary_L36c_20260918.md` (the
H II row solved in the run), and it closes that memo's section 2, "what the
helium stages and the atomic gas still need", on the helium side.

Every number is MEASURED (run or computed here) unless marked READ.

---

## 0. Verdict, before the detail

1. **The carried set is now x(H II) per hydrogen nucleus and x(He II),
   x(He III) per helium nucleus.** The three fractions of L12a section 0
   item 1 are transported rows of the transport-chemistry operator, each on
   its own element's nucleus face flux, and each is handed to every
   ionization sweep the run makes. The neutral stage of each element closes
   its simplex and is not a row.

2. **The stage sum identity holds for BOTH elements on a real solved state,
   to the rounding of the sums.** The certification carries one entry per
   element, `ionization stage nucleus sum H` and
   `ionization stage nucleus sum He`; on the solved states of the one
   fixture that runs the key through the partitioned stationary route they
   read 2.218e-16 (H) and 4.266e-16 (He), against 2.21e-16 and 3.35e-16 on
   the suite's synthetic column and against 1.5e-4 to 3.6e-4 when the
   construction is broken. Both entries REPORT and do not gate (section 5).

3. **The anchoring measurement L12a section 2.4 and L36c owed is made, and
   it re-derives the value the stage rows already carry.** A manufactured
   ionization column whose continuous stationary solution is known, refined
   250/500/1000 over one fixed domain, gives an observed order 2.886 and a
   row residual of 1.333e-7 at dr = 1.0e-3 R_p. Extrapolated to the spacing
   the carrier gate was anchored at (dr = 1.8447e-3 R_p, READ) it is
   7.80e-7, and a decade above it is 7.8e-6. `cert_tol_carrier_wind` stands
   at 1e-5 (READ), so the memo's rule gives that value for a stage row as
   well: it is anchored for these rows and not borrowed. **No tolerance was
   set or changed.** Section 4.

4. **The atomic gas is STILL refused, and the blocker is not the imposition
   gate L36c named.** The imposition is lifted -- the stage fractions are
   now imposed on the key alone and not on the molecular carrier
   configuration -- and that was necessary but is not sufficient: the
   transport-chemistry operator has two entry points, both in
   `src/EXHALE_main.f90`, and both are gated on
   `thereis_mol .and. carrier_transport`; the frozen background the operator
   reads (`bg_cell`, `bg_ready`) is filled only inside the molecular cell
   state; and the row source of every stage is `mol_heh_rows`, the molecular
   H/He network. Section 7 names every site. The refusal is kept and its
   message now states this reason. MEASURED: the certified atomic fiducial
   with `Ionization transport: True` aborts in `input_read`.

5. **The key-on stationary solve is NOT certified, and the refusing row is
   `carrier balance He+`.** On `carrier_model_a_newton` through the
   partitioned route, pass budget 25, the He+ row measure runs 9.99e-1 at
   pass 1 to 7.74e-1 at pass 25 without a monotone fall, while the control
   build (hydrogen alone carried) takes its H+ row from 2.17e-1 to 5.71e-4
   over the same 25 passes. Section 6 reads the row term by term: the He+
   row is carried entirely by its chemistry, whose net source stands 68
   times the stage flux divergence at the gated cell and 5.0e4 times it in
   the shielded layer. **The two builds' final compositions agree to 0.2 per
   cent** (He II at the gated cell, 1.2951e3 against 1.2977e3 cm^-3), and
   the same operator evaluating the same row on the CONTROL's own state
   reads **4.15e-2** instead of 7.74e-1: the row is not broken, the state it
   is given is 0.2 per cent from the one that would close it, and even the
   locally solved composition is three and a half decades above the row's
   1e-5 gate. No composition either route produces certifies that row on
   this fixture.

6. **Key-off byte identity: MEASURED** on `wasp_full` (400 steps),
   `mol_diffusion` and `mol_carrier` (12000 steps each): the data rows of
   `Hydro_ioniz.txt` and `Ion_species.txt` are byte-identical to the control
   build on all three (section 3). `hp_front` carries the key ON and moves: the
   largest movement is 1.56e-1 on n(He III) at the H2 front, which is the
   helium ionization partition changing from a local root to a transported
   one. No golden was refreshed.

7. **The `carrier_retry` row that L36c left RED is diagnosed and is not a
   defect.** The binding species, HeH+ at 1.3e-12 of the gas in the base
   cell, sits at the end of a chain whose links the further sweep moves by
   1.47e-8 (H2+), 1.05e-7 (H3+) and 1.04e-6 (HeH+) while every parent
   density moves at 1.3e-11 or below, and H3+ carries half the electrons of
   that cell, so the chain closes on itself. A THIRD sweep moves the
   composition by 1.01e-10, four decades below the second: the state
   SETTLES, so what the row measures is a seed that needs one more sweep and
   not a repeating arithmetic. The tolerance and the 1e-20 presence
   threshold were NOT touched. Section 8.

---

## 1. What the helium rows are, term by term

The unknowns of the two new rows are

```
    x(He II)(j)  = n(He II)(j)  / n_He(j) ,
    x(He III)(j) = n(He III)(j) / n_He(j) ,
```

`n_He` the helium NUCLEUS density of the cell, from the one stoichiometric
map of the species table (`element_nucleus_counts`: He I + He II + He III +
HeH+, with He 2^3S skipped as an excited level inside He I). It is frozen
with the rest of the background of a step (`cbg_nHenuc`), for the reason the
hydrogen count is: transport moves a nucleus from cell to cell, and which
stage it sits in is what the row solves.

Each face flux is equation (1) of `ionization_stage_transport` written on
the HELIUM nucleus face flux `N_He(f)` of the same evaluation of the element
operator that gives hydrogen its `N_H(f)`:

```
    F_k(f) = x_k(f) N_He(f)  -  n_He(f) K(f) [x_k(j+1) - x_k(j)]/dr(f) .
```

`F_k` depends on `x_k` alone -- the stages of an element couple only through
the closing stage, which is not a row -- so the two helium rows are assembled
one at a time, exactly as the hydrogen row is, and the coupling appears where
it belongs:

| where the two helium stages are ONE object | what it means |
|---|---|
| the closing face fraction | `x(He I)(f) = 1 - x(He II)(f) - x(He III)(f)`, so the reconstructed face fractions of the element sum to one |
| the closing eddy term | minus the sum of the two carried ones, which is condition C3 of the identity |
| the sum identity | `F_HeII + F_HeIII + F_close = N_He` at every face, reported per element |
| the admissible set | ONE simplex, `x(He II) + x(He III) <= 1`, projected together: taking each to [0,1] on its own would admit a cell whose two ionized stages hold more helium than the cell has (the RED row of section 2) |
| the write-back | the helium nuclei the two stages did not take go to He I; the helium of HeH+ is left where it is, because rescaling that molecule would move a hydrogen nucleus with it, and He 2^3S is rescaled with He I so that the level keeps its share of the neutral helium |

The chemistry source of each row is the corresponding row of
`mol_heh_rows`, taken whole: row (2), the complete He+ balance
(photoionization and collisional ionization of He I, the He 2^3S channels,
radiative and dielectronic recombination, the molecular sinks R17 and R23
and the He+ + CO channel), and row (3), the He++ balance. There is no helium
chemistry written in the transport operator that the local solve does not
already have; the difference is only where the row is evaluated. The trial
densities of both stages, the neutral helium that closes the element and the
electron density they contribute to are moved to the trial values inside the
solve, for the reason the proton's are: a row evaluated at the last sweep's
value has a zero Jacobian column for its own unknown.

`mol_heh_rows` hands out the production/loss SPLIT for rows (1) and (4)
alone, so the diagnostic record `EXHALE_CARRIER_ROW_TERMS` carries the
helium rows' NET on the side its sign falls and not the two magnitudes. The
row scale is built from the net in every case (`carrier_residual`), so
nothing in the solution or the acceptance reads the split; section 6 says
what that costs a fast stage.

## 2. Tests, RED and GREEN

`src/tests/ionization_stage_flux/run.sh`, private `OBJDIR=build_L36d`,
`OMP_NUM_THREADS=1`: **TOTAL 38 passed, 0 failed** (MEASURED). The 27 rows of
L36c are unchanged. Eleven are new: nine in one block that runs the
PRODUCTION operator (`carrier_stage_face_state`,
`carrier_stage_face_flux`, `ionization_stage_nucleus_sum` and
`carrier_ionization_stage_projection` are now public for it), and two on
the manufactured column of section 4:

| row | measured | tolerance |
|---|---|---|
| `operator_stage_rows_ride_on_their_own_element_flux` | 0.0000 (bitwise) | 1e-14 |
| `operator_stage_sum_is_the_hydrogen_nucleus_flux` | 2.2113e-16 | 1e-13 |
| `operator_stage_sum_is_the_helium_nucleus_flux` | 3.3525e-16 | 1e-13 |
| `charging_helium_to_the_hydrogen_flux_breaks_the_sum` (RED by construction) | 1.2614e+00 | must exceed 1e-3 |
| `helium_stages_are_returned_to_one_shared_simplex` | 0.0000 | 1e-15 |
| `the_projection_keeps_the_ratio_of_the_helium_stages` | 1.6653e-16 | 1e-14 |
| `the_hydrogen_stage_inside_its_simplex_is_untouched` | 0.0000 | 0 |
| `projecting_each_helium_stage_alone_leaves_the_sum` (RED by construction) | 4.0000e-01 | must exceed 1e-3 |
| `the_carried_set_is_x_HII_x_HeII_x_HeIII` | 0 | 0 |
| `the_manufactured_stage_row_converges_at_second_order_or_better` | 2.886 | >= 2 |
| `the_manufactured_row_falls_with_the_grid` | 1.3327e-07 against 7.2833e-06 | strict |

**RED before.** The nine operator rows cannot be written against the entry
text at all: `ic_HeII` and `ic_HeIII` do not exist there, so the suite does
not compile. What stands in the tree instead are the two RED rows above,
which measure the two mistakes the construction avoids and keep measuring
them: charging the helium stages to the hydrogen flux breaks the sum
identity by 1.26 of itself, and projecting each helium stage onto [0,1] on
its own leaves a cell holding 1.4 of its helium nuclei in ionized stages.

`src/tests/ionization_imposed_fractions/run.sh`: **48 assertions passed**
(31 on the entry text, READ from L36c). The new ones state the three-flag
contract of the substitution, each measured in every system that has the
rows:

| row | what it states |
|---|---|
| `routine_writes_the_constraint_in_row_two` / `..._three` | rows 2 and 3 are `x(2) - x_heii_fix` and `x(3) - x_heiii_fix`, exactly |
| `routine_writes_no_row_beyond_the_third` | nothing else is touched |
| `routine_leaves_the_unflagged_helium_rows` | the proton carried alone leaves the helium rows as the system wrote them |
| `routine_leaves_the_unflagged_proton_row` / `routine_writes_row_two_without_row_one` | and the other way round: three independent flags |
| `<system>_imposed_HeII_row_is_the_constraint` / `..._HeIII_...` | in each of the six systems that offer rows 2 and 3 |
| `imposed_helium_rows_agree_across_the_six` | one substitution, written once |

**RED before**: `x_heii_fixed` and `x_heiii_fixed` do not exist in the entry
text, so these rows do not compile against it either.

Every suite whose driver links a changed object, on the measured build:

| suite | result |
|---|---|
| `ionization_stage_flux` | 38 passed, 0 failed |
| `ionization_imposed_fractions` | 48 passed, 0 failed |
| `carrier_reference_scales` | 16 passed, 0 failed |
| `steady_species_rows` | 197 passed, 0 failed |
| `element_operator` | 40 passed, 0 failed |
| `species_face_flux` | 16 passed, 0 failed |
| `carrier_boundary_jacobian` | 12 passed, 0 failed |
| `carrier_outer_boundary` | 16 passed, 0 failed |
| `carrier_constraint_attribution` | 10 passed, 0 failed |
| `carrier_returned_state_acceptance` | 36 passed, 0 failed |
| `certification` | 84 passed, 0 failed |
| `attempted_step` | 70 passed, 0 failed |
| `constrained_network_layout` | 8 passed, 0 failed |
| `acceptance_classes` | 22 passed, 0 failed |
| `adv_static_limit` | 53 passed, 0 failed |
| `coupled_source_step` | 30 passed, 0 failed |
| `energy_update` | 53 passed, 0 failed |
| `fuv_band_ledger` | 24 passed, 0 failed |
| `h2_level_ladder` | 32 passed, 0 failed |
| `krylov_and_dogleg` | 334 passed, 0 failed |
| `low_mach_stress_energy` | 10 passed, 0 failed |
| `molecular_seed` | 25 passed, 0 failed |
| `physics_probe` | 1604 passed, 0 failed |
| `run_mode` | 31 passed, 0 failed |
| `species_masses` | 9 passed, 0 failed |
| `state_mapper` | 76 passed, 0 failed |
| `carrier_retry` | 151 passed, 1 RED (section 8) |
| `grid_and_gates` | 251 passed, 4 failed, all four on the entry text's own refusal (section 9 item 5) |
| `residual_determinism` | 5 passed, 1 failed, the same 1.715e+04 on the control build (section 9 item 6) |

## 3. Key-off byte identity, and the key-on movement

Control build `EXHALE_L36d_ctl.x`, md5 `91a0dca3daa3e72cf6ea5b61a16c52d6`,
a build of the ENTRY TEXT of every file of this item, which is the binary
the L36c increment delivered. Measured build `EXHALE_L36d_measured.x`, md5
`cf3d2941d1ebe7864f14bd411d5bdfb9`. Scratch case copies,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, identical step bounds on the two
builds.

| case | key | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|---|
| `wasp_full` | off | 400 | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_diffusion` | off | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_carrier` | off | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `hp_front` | **on** | 100 | 1.600e-02 | 1.561e-01 |

`mol_carrier` and `mol_diffusion` run the carrier transport operator itself,
which is the operator this item adds two rows to, so they are the cases that
reach the changed code with `carrier_solved(ic_HeII)` false; `wasp_full`
reaches the changed sweep through the He 2^3S and metal configuration with no
molecules at all.

**`hp_front` carries the key and moves, and the movement is the physics this
item adds.** Column by column, the largest relative difference over the 504
data rows:

| column | largest movement |
|---|---|
| He III | 1.561e-01 (cell 205, 2.471972e+00 to 2.086097e+00 cm^-3) |
| He I | 2.563e-02 |
| He II | 4.271e-02 |
| He 2^3S | 4.368e-02 |
| HeH+ | 4.368e-02 |
| H2+ | 4.025e-02 |
| H3+ | 3.201e-02 |
| cooling rate | 1.600e-02 |
| heating rate | 3.341e-03 |
| every other column, H I and H II among them | below 1e-3 |

Nothing of hydrogen moves above 1e-3: what changed is that the helium
ionization partition is now what the flow accumulated instead of the local
root of each cell, and the molecular ions that are functions of it follow.

**The state stays inside the simplex, MEASURED** on the final `hp_front`
state of the measured build: no species column is negative anywhere,
x(He II) spans 5.10e-12 to 0.6500, x(He III) 4.88e-16 to 0.02451 with
x(He II) + x(He III) at most 0.6745, and x(H II) spans 1.53e-08 to 0.7716.

**The helium write-back moves charge and not nuclei, MEASURED**: over the
500 physical cells of the final `hp_front` state the elemental ratio
n(He)/n(H), formed from the output columns, is 0.0793000000 with a spread of
**2.945e-12** on the measured build against 2.895e-12 on the control, so
transporting the two helium stages leaves the element total where the
element operator put it.

`hp_front`'s golden is therefore due a refresh together with the two goldens
L36c already reported as stale; **no golden was refreshed here**, `EXHALE.x`
was never touched and nothing was written under `backup/regression/`.

## 4. The anchoring measurement, and the anchor

L12a section 2.4 (READ) states what is owed: "a smooth manufactured
ionization column (a chosen `x_eq(r)` and the `P`, `L` that produce it),
refined by a factor two at the production spacing, giving the discretization
error of the row itself. The gate is then a decade above that, as the
carrier gate is a decade above its own anchor."

**The column** (`test_manufactured_ionization_column`, in the stage-flux
suite, so it is repeated whenever the suite runs). A spherical column over
[1, 2] whose faces are the same two radii at every resolution, with

```
   n_el(r) = n_0 exp[-a(r-1)] ,  N_el(r) = C/r^2 ,  K(r) = K_0 r^2 ,
   x(r)    = x_0 + dx (6 s^5 - 15 s^4 + 10 s^3) ,   s = r - 1 ,
```

`n_0 = 1e10 cm^-3`, `a = 4`, `C = 1e13 nuclei cm^-2 s^-1`,
`K_0 = 1e12 cm^2 s^-1`, `x_0 = 0.02`, `dx = 0.9`. `N_el` is divergence free,
so the column carries a steady nucleus flux; the quintic is monotone, so the
limiter never fires, and it has `x' = x'' = 0` at both ends, so the
operator's own boundary rule (no eddy flux through the two end faces) is
EXACT there. The continuous flux of equation (1) and its exact divergence
follow analytically, and that divergence is the source that makes `x(r)` an
exact stationary solution of the continuous row. What the DISCRETE row reads
on it is the truncation of the discretization and nothing else. The measure
is the row residual over the row's own terms, which is the measure the
certification gates: the two face contributions and the source.

MEASURED, interior cells (the two end cells carry the boundary rule and are
reported separately):

| cells | dr [R_p] | row measure | with the end cells |
|---|---|---|---|
| 250 | 4.000e-3 | 7.2833e-06 | 8.6047e-05 |
| 500 | 2.000e-3 | 1.0135e-06 | 2.1983e-05 |
| 1000 | 1.000e-3 | 1.3327e-07 | 5.5589e-06 |

observed order **2.845** (250 to 500), **2.927** (500 to 1000), **2.886**
over the pair. This is the order of the NORMALIZED measure and not of the
truncation: the scale it divides by contains the two face MAGNITUDES over
the cell volume, about 2|F|/dr for a smooth nonzero flux, so it grows as
1/h and the normalized order stands one above the dimensional one. MEASURED
on the same three grids (item D3a, `docs/lhs1140b_stationary_D3a_20260918.md`
section 5.1): the dimensional residual is 3.7493e-2, 1.0076e-2 and
2.6089e-3 cm^-3 s^-1, an observed order of **1.923**, which is the
conservative divergence's second order. The extrapolation below is an
extrapolation of the normalized measure by the normalized order, which is
the quantity the certification gates, so the anchored value is unchanged.
The advective and eddy halves of the source are within a factor 2.38 of each
other, so neither half carries the measurement alone.

**The anchor.** The carrier gate `cert_tol_carrier_wind = 1e-5` (READ) was
anchored on a smooth manufactured ELEMENT column reading 1.52e-6 at
dr = 1.8447e-3 R_p, the spacing of the atomic fixture's grid at the binding
radius (both READ, `docs/certification_tolerance_anchoring_20260910.md`
section 2 (3) and its decision 22 (a)). The stage column above extrapolated
to that same spacing with its own order reads **7.80e-7**, and a decade
above it is **7.8e-6**.

**So the value for a stage row is 1e-5, and it is the value the stage rows
already carry.** The 1e-5 is not borrowed from the transport rows against
L12a's warning: it is re-derived from a measurement of an ionization row.
**No tolerance was set or changed by this item.**

Two readings that belong with it:

- On the fixture's own Mixed grid the spacing at the radius from which a
  species row is gated (`cert_regime_wind_r = 1.20`) is 3.5215e-3 R_p, where
  the same column extrapolates to **5.04e-6**. Applied there, the memo's
  rule would give 5e-5, which is LOOSER than the gate in place.
  `PLAN_20260917.md` section 12 forbids weakening a certification tolerance,
  so the tighter of the two readings stands and nothing is changed; what the
  measurement adds is that 1e-5 is a factor 2 above the smooth-column
  truncation at the radius where it first gates, so it is a tight gate for
  these rows rather than a loose one.
- A smooth column is a LOWER bound on the production discretization, as the
  anchoring memo says of its own: a real ionization front is steeper than
  this quintic.

**The stage-sum entry takes no tolerance.** The identity is exact in exact
arithmetic, so what a tolerance has to clear is the rounding of the sums, and
the anchoring memo's five anchors are all measurements of an equation with a
truncation error -- none of them is a rule for an identity that has none.
MEASURED: 2.21e-16 and 3.35e-16 on the suite's synthetic 500-cell column,
2.218e-16 (H) and 4.266e-16 (He) on the solved state of section 6, against
1.48e-4 to 3.62e-4 when the construction is broken (READ, L36). The proposal
of L36c section 7, **1e-13**, stands as a proposal; the entry REPORTS and
does not gate.

## 5. The certification entries

One entry per element, `ionization stage nucleus sum H` and
`ionization stage nucleus sum He`, printed beside the carrier rows they were
formed with. The identity is a statement about ONE element's own nucleus
flux, so a single entry over both could not say which element's construction
it belongs to, which is the reason each trace element's transport balance is
its own entry. Both are `cert_evaluated` with `tol = 0` and
`within_tol = .true.`: they report and do not gate, with the reason printed
("reported, not gating: no anchored tolerance").

MEASURED on the solved state of section 6:

| entry | max | cell |
|---|---|---|
| `ionization stage nucleus sum H` | 2.218e-16 | face 460 |
| `ionization stage nucleus sum He` | 4.266e-16 | face 267 |

The helium entry is the one this item makes possible: with the helium stages
uncarried it reads `cert_unavailable` and says so.

## 6. The key-on stationary solve, and the row that refuses it

`backup/regression/carrier_model_a_newton` on a scratch copy: the one
fixture that carries `Ionization transport: True` with
`Restart intent: stationary`, `Solver: Newton 100.0` and
`Coupled carrier solve: False`, which is the partitioned route. Run with
`EXHALE_PTC_DTAU0=1.0`, `EXHALE_OUTER_PASSES=25`, `OMP_NUM_THREADS=8`, on
both builds from the same restart.

**The state is NOT certified on either build.** The refusing row and its
pass history:

| pass | control, worst gated species row | measured build, worst gated species row |
|---|---|---|
| 1 | 2.17e-01 `carrier balance H+` | 9.99e-01 `carrier balance He+` |
| 2 | 3.07e-02 | 9.95e-01 |
| 3 | 1.25e-02 | 9.85e-01 |
| 12 | 5.30e-03 | 7.77e-01 |
| 13 | 1.75e-03 | 6.67e-01 |
| 20 | 1.06e-03 | 5.09e-01 |
| 25 | 5.71e-04 | 7.74e-01 |

every pass at `hydro info=0`, every pass at cell 248, r = 1.2013, the first
cell the species rows are gated at. The control's H+ row falls by two and a
half decades over the 25 passes; the measured build's He+ row does not fall.

**The final certification, measured build:**

| entry | max | cell | gated value at r >= 1.20 | verdict |
|---|---|---|---|---|
| `carrier balance H2` | 5.334e-06 | 290 | 5.334e-06 | within |
| `carrier balance H+` | 1.172e-02 | 124 | 2.256e-04 at cell 248 | ABOVE |
| `carrier balance He+` | 1.000e+00 | 107 | 7.736e-01 at cell 248 | ABOVE |
| `carrier balance He++` | 5.767e-03 | 145 | 1.795e-04 at cell 248 | ABOVE |
| `ionization stage nucleus sum H` | 2.218e-16 | 460 | reported | within |
| `ionization stage nucleus sum He` | 4.266e-16 | 267 | reported | within |
| hydrodynamic mass / momentum / energy | 4.617e-11 / 5.442e-13 / 3.355e-09 | | | all within |

against the control's `carrier balance H+` 5.165e-02 (cell 124), gated
5.708e-04 at cell 248. **The proton row is BETTER on the measured build at
the gated cell** (2.256e-04 against 5.708e-04) and the He++ row is of the
same size; it is the He+ row alone that refuses.

**What the He+ row is made of**, read term by term from the row record
(`EXHALE_CARRIER_ROW_TERMS=1`) of the state that was handed back, every rate
a volumetric rate in cm^-3 s^-1:

| cell | r [R_p] | n(He+) | x(He II) | stage flux divergence | net chemical source | residual | measure |
|---|---|---|---|---|---|---|---|
| 107 | 1.0263 | 2.0884e+01 | 1.0087e-11 | 5.068e-06 | +2.5127e-01 | -2.5127e-01 | 0.99996 |
| 248 | 1.2013 | 1.2951e+03 | 6.7967e-08 | 1.4595e-02 | +9.9249e-01 | -9.7789e-01 | 0.9710 |

and, for comparison, the He++ row of the same cells, whose transport and
chemistry DO balance:

| cell | stage flux divergence | net chemical source | residual | measure |
|---|---|---|---|---|
| 107 | -1.3170e-10 | -1.3370e-10 | +2.001e-12 | 2.755e-03 |
| 248 | +5.6722e-05 | +5.6866e-05 | -1.441e-07 | 1.269e-03 |

The record is written by the certification call that produced it, and the
gated numbers of the table above come from the pass-25 certification, so the
two are two evaluations of the same row one composition refresh apart (0.971
against 0.774 at cell 248); the reading does not depend on which of them is
taken.

So the He+ row is carried entirely by its chemistry: the net source stands
68 times the stage flux divergence at the gated cell and 5.0e4 times it in
the shielded layer, and the stationary row therefore demands that the
transported x(He II) sit on the local root of the helium balance to that
many decades.

**And the two builds' compositions agree.** MEASURED at the same cells on
the final states:

| cell | He II, control | He II, measured | He III, control | He III, measured |
|---|---|---|---|---|
| 107 | 2.0932e+01 | 2.0884e+01 | 8.1262e-04 | 8.1973e-04 |
| 248 | 1.2977e+03 | 1.2951e+03 | 4.4634e+00 | 3.5714e+00 |

0.23 per cent in He II at cell 107 and 0.20 per cent at cell 248, against a
He+ row measure of 1.0 and 0.97. **The row refuses a state whose helium is
within a fifth of a per cent of the locally solved one**: the row measure
falls from 0.97 to 4.15e-2 when the composition moves by that 0.2 per cent
(the control below), so the He+ balance is stiff in its own unknown, its net
source changing by of order its whole size for a fifth of a per cent in
n(He II). The separate production and loss that would show this directly are
not available: `mol_heh_rows` returns rows (2) and (3) as one balance. Two things follow, and they are reported and not fixed here:

- The bounded alternation does not remove that last fraction of a per cent.
  The carrier relaxation's movement bound stood at `trust 2.5e-3` for the
  later passes, and 25 passes of at most that do not close a gap the sweep
  closes in one pass when the stage is not carried.
- The row's SCALE is the magnitude of the NET source (`carrier_residual`
  builds it as `|div F| + |src|`), and for a stage whose production and loss
  cancel to three decades the net is the same size as the residual. The
  module's own header says this of the diagnostic split -- "where the
  chemistry is fast the two cancel to many digits and |src| stands orders
  below either, so a row that reads as one term is not a row with one term"
  -- and the same statement applies to the scale the acceptance divides by.
  Whether a fast ionization stage should be measured against its production
  plus its loss instead is a change to a shared acceptance interface with its
  own anchoring, and it is NOT made here.

**A control on the row itself.** The same operator evaluated the same three
rows on the CONTROL's own final state, through `Restart intent: stationary
evaluate` with the MEASURED binary, so that only the STATE differs
(MEASURED, gated value at r >= 1.20):

| row | on the control's state (He+ the sweep's local root) | on the transported state |
|---|---|---|
| `carrier balance H+` | 5.708e-04 | 2.256e-04 |
| `carrier balance He+` | **4.151e-02** | 7.736e-01 |
| `carrier balance He++` | 7.765e-02 | 1.795e-04 |
| `ionization stage nucleus sum H` / `He` | 2.201e-16 / 4.336e-16 | 2.218e-16 / 4.266e-16 |

So the He+ row of this operator is not a broken row: given a composition
whose He+ IS the local root it reads 4.15e-02, a factor 19 below what the
transported state gives it. And it is above its 1e-5 gate by three and a half
decades even there, so on this fixture no composition either route produces
certifies that row.

**Which regime this fixture is.** L12a's measurement that decided the
carried set was taken on the LHS 1140 b wind, where the Damkohler numbers of
all three fractions fall below 0.3 above 4 R_p (READ, L12a section 0 item 1).
This fixture is a hot Uranus whose domain ends at 7.24 R_p and whose gated
cell is at 1.20 R_p, where the He+ chemistry is fast by the two ratios in the
table above. It is the only fixture that admits the key with the partitioned
route, and it is not the configuration the carried helium stages were
designed for.

## 7. The atomic gas: what the refusal is now, and what it takes

The imposition gate L36c section 2 item 3 named is lifted: the stage
fractions are imposed on `ionization_transport .and. (bg_ready .or.
do_load_IC)` alone, with no molecular carrier condition, so a sweep that runs
in any configuration is handed the transported partition. That was necessary
and it is not sufficient. MEASURED: the certified atomic fiducial
`LHS1140b/models/.L26/fid_resolve` with `Ionization transport: True`
appended aborts in `input_read` with the refusal below, and the refusal is
correct.

What an atomic run still needs, all of it READ from the source:

1. **The operator has no caller.** Both entry points of the
   transport-chemistry operator are in `src/EXHALE_main.f90` and both are
   gated on `thereis_mol .and. carrier_transport`: the operator-split step of
   the marching loop (the call at line 2499 through
   `photochemical_transport_step`, whose interval returns at once on
   `.not. carrier_transport`) and the fixed-wind relaxation of the
   partitioned stationary route (`relax_photochemical_composition` at line
   7210, in the block opened at line 7203). The same condition also sets
   `species_alternated` and the pass cap of that route (line 6763) and
   gates its diagnostics and its progress control at lines 4704, 5953, 6023, 6422, 7079, 7087, 7301, 7426, 7487,
   7501, 7526, 7555 and 7651. Fourteen sites in a file this item does not
   own, and the change is to the control flow of the stationary outer
   iteration, not a line.
2. **The operator has no frozen background.** `bg_cell(j) = ieq_cell` sits
   inside `if (thereis_mol)` (`ionization_equilibrium.f90` line 2049) and
   `bg_ready` is set only `if (thereis_mol)` (line 2784), while
   `carrier_state` reads `bg_cell(j)%ntot` and `%T_K` for every cell and
   `carrier_transport_interval` and `carrier_steady_residual` both return at
   once on `.not. bg_ready`.
3. **The row source is the molecular network.** `carrier_source` takes every
   stage's chemistry from `mol_heh_rows` (`System_HeH_mol`), which reduces to
   the atomic balances at zero molecular densities but is the network of a
   molecular run; whether an atomic run should evaluate it, or take its rows
   from `System_HeH`, is a design decision and not a mechanical edit.

The refusal's message now states (1) instead of the sweep, which was L36c's
reason and is no longer the binding one.

The two refusals that are KEPT for what the operator does not carry are
unchanged: a mixture with no helium (the element nucleus flux is the binary
H/He operator's and returns nothing for one element, and the carried set now
includes two helium stages), and `Coupled carrier solve: True` (its row
registry carries every carrier unknown as the species mass fraction of its
`f_sp` column, and a stage unknown is a fraction per element nucleus).

**What could not be measured because of this.** The He I 10830 equivalent
width from a transported composition: the atomic fiducial is refused and the
one fixture that runs the key has `Include He23S? False`. L36's factor 13.12
between the `_adv` and the solved composition on that fiducial (20.3950
against 1.5544 percent A, READ) stands as the size of the thing the carried
fractions move; it is not yet the answer they give. **Stage D, the 5 per cent
He 2^3S promotion rule of L12a section 3, still has no precondition**: the
rule asks what the post-process does to the He 2^3S column when the three
ionization fractions already agree to their own gate, and on the one fixture
that runs the key the He+ row is at 0.97 of its scale (section 6).

## 8. The `carrier_retry` row L36c left RED

`one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
reads **1.044879e-06** against `ieq_res_tol` = 1e-6 on the measured build,
against 1.034627e-06 on the entry text and 2.685237e-07 one increment before
it (both READ, L36c). The binding species is HeH+ at 1.33e-12 of the gas in
the base cell. **It is a seed that has not settled, not the imposition and
not the sweep's own noise.** Four measurements, all new diagnostics of the
fixture and none of them an assertion:

1. **The chain.** At the binding cell the further sweep moves, relative:
   H I 1.32e-11, H II 1.39e-16, He I 2.30e-16, He II 2.13e-16, H2 1.86e-16,
   **H2+ 1.47e-08, H3+ 1.05e-07, HeH+ 1.04e-06**, and the electron density
   by 5.53e-08. Every parent of the HeH+ balance moves at 1e-11 or below and
   the three molecular ions move in a chain with a factor 7 to 10 per link.
2. **The chain closes on itself.** The densities of that cell are
   n_e = 1.197e+05, H+ = 5.636e+04, H2+ = 3.450e+01, H3+ = 6.312e+04,
   HeH+ = 1.226e-02 cm^-3: **H3+ carries 53 per cent of the electrons** and
   the imposed proton the rest, so re-solving the trace molecular ions moves
   the electron density that their own balances read.
3. **It settles.** A THIRD sweep, from the state the second one left, moves
   the composition by **1.0085e-10**, four decades below the second, its
   worst row being H2+ at cell 12. A repeating arithmetic would repeat; this
   does not.
4. **And it is not the imposition.** The same further sweep with the
   transported fractions NOT imposed moves the composition by **1.5067e+08**,
   its worst being He III at cell 4: the local root of the helium ionization
   at that cell is eight decades from the transported value. The imposition
   is what holds the state, not what breaks it.

**The tolerance and the 1e-20 presence threshold were NOT touched**, and the
fixture decision -- raise the threshold above a species thirteen decades
below the gas, or accept the row as a statement about the first sweep after a
composition change -- remains the advisor's.

## 9. Noticed outside this item's scope

1. **Two `carrier_retry` rows about the superseded row scale needed a
   shorter step.** `the_full_term_scale_accepts_the_imbalance_at_a_short_step`
   and `the_full_term_measure_of_the_same_defect_falls_with_the_step`
   demonstrate that a measure carrying the time term lets a shortened step
   buy an acceptance. The demonstration ran at dt/16 and the acceptance reads
   the WORST row of the system; with the ionization stages carried that row
   is a helium stage, whose time term is a smaller share of its own scale, so
   a 1e-7 physical imbalance read 4.10e-8 at dt and 1.31e-8 at dt/16 and no
   longer crossed the fixed 1e-8 floor. The demonstration now runs at dt/64,
   where it reads 4.12e-9, and both rows are GREEN. **No tolerance was
   changed**: what moved is how short the short step of the demonstration is,
   and the paragraph beside it says why.
2. **`carrier_ionization_stage_projection` wrote into an unallocated cell
   record.** It marks `pct_cell_constrained`, which only
   `limit_to_element_budget` allocated, so a caller that projected on its own
   segfaulted. FIXED: the routine allocates it if it is not there. Found by
   the new acceptance row that calls it directly.
3. **`constrained_chemical_equilibrium.f90` had to learn the two helium
   constraints.** Its continuation judges the fraction system whose rows 2
   and 3 are now `x - x_fix`, so solving them as reaction rows would return
   the local helium ionization root and be measured against rows that were
   never asked for it, and its seed would start off its own constraint
   surface. `species_fixed(is_HeII)`, `species_fixed(is_HeIII)` and the seed
   are set the way the proton's already were, and the cell-state dump the
   file writes for debugging carries the two new flags on a record of their
   own. That file is not in this item's list; it is reported here because
   leaving it would have made the continuation solve a different system from
   the one it is judged against.
4. **`mol_heh_rows` hands out the production/loss split for rows (1) and (4)
   only.** The helium rows therefore have no separate production and loss for
   the `EXHALE_CARRIER_ROW_TERMS` record, and what it carries for them is the
   net on the side its sign falls, so that `production - loss` is the row's
   own source for every row. The two magnitudes of a helium row would need
   the split to be added at the one place the terms are written
   (`System_HeH_mol.f90`), which cannot be done without restating the
   summation order of rows (2) and (3) and moving every molecular golden;
   it is reported and not done.
5. **Four `grid_and_gates` rows fail on the ENTRY TEXT's refusal, not on
   this item.** `restart_option_change.sh`'s route_change block builds its
   rungs from `backup/regression/carrier_model_a_newton`, which carries
   `Ionization transport: True`, and its rung B forces
   `Coupled carrier solve: True` on it, which the refusal L36c added aborts.
   The four rows are `route_change_block_loads_an_alternation_state`,
   `route_change_line_written_in_both_halves`,
   `route_change_written_state_records_its_own_route` and
   `route_change_history_inherited_by_the_next_rung`; the rung's own log
   carries the refusal. `restart_option_change.sh` is unmodified in the
   working tree and predates L36c. NOT FIXED: the fix is to point that block
   at a fixture that admits the block route, which changes what the rows
   test, and `grid_and_gates` is being edited by another item.
6. **`residual_determinism` reads 1.715e+04 on
   `closure_spread_within_the_row_tolerance_atomic_elem_newton` on BOTH
   builds** (MEASURED here with `EXHALE_RESID_EXE` set to the control
   binary). Pre-existing; that fixture carries neither the ionization key
   nor the molecular network.
7. **`docs/input_schema.md` K15f is made false by this item in five places**
   -- the carried set, the identity entries, the admissible set, the
   refusal's reason and the sentence "x(He II) and x(He III) are NOT
   carried" -- and is REWRITTEN in those five places; the rest of the entry
   is left as it stands.
8. **`hp_front`'s golden is now stale by 1.6e-1 in the helium columns**, and
   joins the two L36c reported (`hp_zero_seed` and `hp_trace_seed`, already
   3.9e-3 and 3.3e-3 stale against the entry text). All three are key-on
   cases; the refresh they need is one refresh at the close of the series and
   is not taken here.

---

## 10. Where the L12 stages stand after this increment

| stage (L12b's sequence) | state |
|---|---|
| B, the stage flux on the shared element flux | landed with L36, and now carried for both elements |
| C, the rows solved in the run | landed with L36c for hydrogen and here for helium; the ATOMIC configuration is the one piece that is not built, and section 7 names every site it needs |
| D, the He 2^3S promotion rule | no precondition: it asks what the post-process does to the He 2^3S column when the three ionization fractions already agree to their own gate, and on the one fixture that runs the key the He+ row is at 0.97 of its scale. The side the level's own Damkohler measurement supports is unchanged (smallest 15.7 anywhere in the physical domain, READ, L12a section 1.1), and the systematic it would be judged against is L36's factor 6.15 in the He 2^3S radial column between the two closures (READ) |
| E, the certification entries and the tolerance anchoring | CLOSED. Each carried stage has its row entry, each element has its stage-sum entry, and the tolerance question of L12a section 2.4 has its measurement and its anchor (section 4). No tolerance was set or changed |
| F, the `_adv` retirement as the source of the carried fractions | landed with L36 on the transit side; the line result from a transported composition still waits on the atomic configuration (section 7) |

## 11. What this item does not do

- It does not make the atomic configuration run: the refusal is kept and
  section 7 names what it takes, including fourteen gate sites in a file
  this item does not own.
- It does not change any tolerance, and it does not gate the stage-sum
  entries.
- It does not refresh a golden, and it names the three key-on cases whose
  goldens the series owes.
- It does not change how a fast ionization stage's row is scaled, which is
  the one thing that decides whether `carrier balance He+` can certify
  (section 6); that is a shared acceptance interface with its own anchoring.
- It does not split the production and loss of `mol_heh_rows` rows (2) and
  (3), which would restate their summation order.

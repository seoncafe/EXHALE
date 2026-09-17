# L21: the base boundary decides its branch on a quantity that is not a flux

Item L21 of `docs/PLAN_20260913_lhs_stationary.md`, opened on the measurement
reported by item L7f and approved by the user on 2026-09-15.

Binary: `EXHALE_L18.x`, the tree of 2026-09-15 with items L18, L19 and L20 in
it. The old behaviour is reached from the same binary with
`EXHALE_BASE_BRANCH_ON_CELL1=1 EXHALE_BASE_MACH_BLEND=1.0e-6`, so every
comparison below is one binary against itself.

## 1. Verdict

**The base face of the certified LHS 1140 b states is an INFLOW face carrying
exactly the wind's own mass flux, and the boundary condition was running it 97
per cent in the REVERSAL branch.** The discriminant was the cell-centred
product `rho_1 v_1 r_1^2` of the first interior cell, and at this base that
product is not a flux: it is the collocated odd-even velocity mode of
`docs/p44_base_sawtooth.md`. Measured on
`atomic_scalar_gj1132_kzz1e9/HeH2.13`, the Riemann face mass flux is
**+1.00000 F_wind at every one of the 500 faces including the base face**
(largest departure 1.6e-9), while the cell-centred product reads **-2.00
F_wind at cell 1 and -2.27 at cell 2**.

**NOTE ADDED 2026-09-17 (item L26 section R9): the sentence above is no longer
true of the states this catalog now carries.** Re-measured on the states as
they stand after the re-run, the Riemann base face mass flux in units of the
window's mean flux reads +32.38 on the fiducial, +336.55 on the certified 0.10,
-404.32 on the certified 0.03, -584.45 on the 0.02 transient and -4735.68 on
the cold start, with the faces just above the base at -18.6, -203.7, -153.9 and
-232.5 before settling near +1 at faces 4 to 8. This record is left as taken,
on the state it was taken on; that the base face flux is not the wind's on the
present states is an open item, attributed to the base odd-even mode of
`docs/p44_base_sawtooth.md` and not to the boundary closure, and recorded in
`docs/lhs1140b_stationary_L26_20260916.md` section R9 and in
`docs/TO_BE_DONE.md`. What does NOT change is the reason this item gives for
keying the branch on a flux rather than on the cell-1 velocity.

The consequence is not small. On the reversal branch the face density is the
interior's isentrope instead of the reservoir's, and the two are not close:
**rho_rev/rho_res = 0.531**, which is exactly the temperature ratio
**T_i/T_res = 420.4 K / 223.2 K = 1.88** at the common face pressure. So
`base.inp`'s 226 K base was entering the state with a weight of 0.026, and the
first cell of a "226 K base" run sits at 417.6 K.

Two things were wrong and only one of them is the window.

- **The sign.** No choice of blend width repairs it: the discriminant is
  negative, so a narrower window takes the reversal branch MORE completely
  (`w_rev = 1.000` at a window of 1e-8) and a wider one tends to 0.5. Keying
  the branch on the mass flux the state carries is the only fix.
- **The width.** `base_face_mach_blend = 1e-6` was chosen with the note "the
  window must be far BELOW the physical operating point" and the hot-Uranus
  base Mach number of 6e-6 quoted for it. **LHS 1140 b's base operating Mach
  number is 4.6e-7 -- BELOW the window.** Even with the sign repaired the old
  width leaves `w_rev = 0.215`, so a fifth of the reservoir's entropy would
  still be blended away on a clean inflow.

Both are fixed. On the same state the branch now reads `M_branch = +4.02e-07`
and `w_rev = 0.000`: the reservoir states the entropy of the gas that enters,
which is what the boundary exists to do.

**The certified LHS 1140 b campaign was solved with this boundary and has to be
re-solved.** The corrected boundary asks the base face for 1.88 times the
density the certified states carry there, and the old certified state is not
near a solution of it: re-entered, its mass row reads **9.41e-01** against the
1.49e-09 it certifies at under the old boundary. Section 5 reports the
re-solve.

## 2. Sign convention, established from the code

`v > 0` is outward. At the base face, gas ENTERING the domain moves outward, so

- `v_face > 0` : inflow, the reservoir states (p, s), `w_rev = 0`;
- `v_face < 0` : reversal, the entropy at the face is the interior's advected
  out, the reservoir states p alone, `w_rev = 1`.

`characteristic_base_face_state` built its discriminant as

```
v_i = W1(1)*W1(2)*r1*r1/(rho_i*rb*rb)     ! cell 1's rho v r^2, mapped to the face
M_i = v_i/c_i
w_rev = characteristic_branch_weight(M_i/base_face_mach_blend)
```

so `sign(M_i) = sign(v(1))`. The branch logic is right; the quantity it is
given is not.

## 3. What the base of the certified state actually carries

`atomic_scalar_gj1132_kzz1e9/HeH2.13`, the certified state re-read with
`EXHALE_RESIDUAL=1` (the residual diagnostic now writes the stored Riemann face
mass flux beside the cell-centred state, so this no longer has to be
reconstructed by integrating the mass row as `p44_base_sawtooth.md` section 3
had to):

| cell | r | v [cm/s] | cell-centred `rho v r^2 / F_wind` | face flux below / F_wind | face flux above / F_wind |
|---|---|---|---|---|---|
| 1 | 1.000193 | **-0.1127** | **-2.004** | **+1.00000** | +1.00000 |
| 2 | 1.000387 | **-0.1348** | **-2.273** | +1.00000 | +1.00000 |
| 3 | 1.000580 | +0.0874 | +1.394 | +1.00000 | +1.00000 |
| 4 | 1.000773 | +0.0474 | +0.717 | +1.00000 | +1.00000 |
| 5 | 1.000967 | +0.0751 | +1.077 | +1.00000 | +1.00000 |
| 6 | 1.001160 | +0.0710 | +0.968 | +1.00000 | +1.00000 |
| 10 | 1.001933 | +0.0889 | +1.000 | +1.00000 | +1.00000 |
| 100 | 1.024957 | +1.392 | +1.000 | +1.00000 | +1.00000 |
| 400 | 5.904338 | +3.10e4 | +1.000 | +1.00000 | +1.00000 |

Two physical cells carry a negative cell-centred product; every face of the
grid carries the same positive flux to nine digits. This is the P44 signature
at a tenth of the hot-Uranus amplitude (-2 F_wind against -196), and the
answer to question (1) of the item is **(a): the odd-even artifact, not a
reversal**.

The boundary's own report on the same state, old behaviour:

```
   interior at the face: v_i =-1.115484E-01 cm/s,  M_i =-8.065711E-07
   the Mach the branch is decided by: M_branch =-8.065711E-07
   face state:           v_b =-3.100953E+00 cm/s,  M_b =-2.268070E-05
   branch weight w_rev = 9.737E-01  (window  1.000E-06 in face Mach)
   rho_res = 9.505841E+13  rho_rev = 5.046861E+13  [mH/cm3], ratio  0.53092
   T_res =   223.195 K   T_i =   420.397 K
```

and with the fix:

```
   the Mach the branch is decided by: M_branch = 4.024422E-07
   branch weight w_rev = 0.000E+00  (window  1.000E-08 in face Mach)
```

The module's own claim that the branch is confined to where "the two isentropes
agree to the boundary's own error" is false at this base: they differ by 1.88,
the ratio of the two temperatures, and the smoothstep was interpolating across
a factor of two.

## 4. What changed in the source

`src/modules/states/base_boundary.f90`.

1. **The discriminant.** `wind_window_mass_flux` returns the mean of
   `rho v r^2` over the wind window `r >= r_flux` -- the window
   `flux_spread_of_state` is already defined on, and where the cell-centred
   product IS the conserved flux -- and the branch weight is keyed on that flux
   mapped onto the face by the interior density,
   `M_branch = F_wind/(rho_i r_b^2)/c_i`, **where that Mach number is outside
   the handover window**. Inside it the wind is not selecting a branch and the
   first interior cell answers, as it did before: a cold start's wind window
   carries no wind, and section 7 measures what happens without that guard.
   The (C-) relation keeps `v_i`: the outgoing acoustic invariant is a genuine
   property of cell 1 and the item does not touch it.
   `EXHALE_BASE_BRANCH_ON_CELL1=1` restores the old discriminant.
2. **The width.** `base_face_mach_blend` 1e-6 -> **1e-8**, 46 times below the
   LHS 1140 b operating point and 600 times below the hot-Uranus one, which is
   what the parameter's own design note asks for and 1e-6 did not deliver.
   `EXHALE_BASE_MACH_BLEND=<value>` overrides it.
3. **The boundary can be asked what it did.** `report_base_face_state` prints
   the interior face state, the Mach the branch was decided by, the weight, the
   two densities it mixes and the two temperatures. Nothing printed them before;
   the counters existed and were never read.

`src/EXHALE_main.f90`: the `EXHALE_RESIDUAL=1` profile gains the two stored
Riemann face mass fluxes of each cell, and calls the report above.

**Why the face flux itself is not read.** It is produced by the Riemann solve
this boundary feeds, so reading it would be either the circular dependency the
module exists to break or history in a residual that has to stay a function of
its argument (`residual_determinism`). The wind window is the same conserved
number, read where the collocated mode is absent. The limitation is stated in
the code: during a transient the base can carry a flux the wind does not yet,
and the branch then follows the wind's sign; no converged state has the two
disagreeing.

**Why not a local estimator.** Every cell-centred estimator was checked against
both bases and none survives: the mean of `rho v r^2` over cells 1-2 is -2.14
F_wind here, over 1-4 it is -0.54, over 1-6 it is -0.02, and on the hot-Uranus
base of `p44_base_sawtooth.md` cell 1 alone is -196 F_wind, so even a mean over
the whole `r < 1.03` layer is negative there. A window that works on one base
does not work on the other, which is why the flux has to be read where it is a
flux.

## 5. The certified case re-solved, and what moves

`atomic_scalar_gj1132_kzz1e9/HeH2.13`, started from its own certified state, in
a `.L21/` copy, 8 threads, `Restart intent: stationary`. The control run is the
same binary with `EXHALE_BASE_BRANCH_ON_CELL1=1 EXHALE_BASE_MACH_BLEND=1.0e-6`
and reproduces the reference (mass 1.486e-09 against the 1.487e-09 of
`REPRODUCE.md`, energy 3.439e-08 against 3.265e-08, CERTIFIED, `info = 0`), so
the comparison is one binary against itself.

**The corrected boundary converges**: `info = 0`, CERTIFIED at outer pass 3,
233 s + 77 s + 79 s of solve.

| row | certified (old boundary) | corrected |
|---|---|---|
| hydrodynamic mass | 1.487e-09 at cell 4 | 2.175e-09 at cell 1 |
| hydrodynamic momentum | 1.170e-12 at cell 500 | 5.417e-13 at cell 500 |
| hydrodynamic energy | 3.265e-08 at cell 8 | 2.126e-08 at cell 1 |
| elemental transport He/H (reported, not gated) | 6.166e-01 at cell 2 | 4.051e-01 at cell 2 |
| verdict | CERTIFIED | CERTIFIED |

**The base thermal structure is a different one, and it is the stated one.**

| cell | r | T certified | T corrected | rho corrected / certified |
|---|---|---|---|---|
| ghost 0 | 1.000000 | 413.7 K | **226.0 K** | 1.857 |
| 1 | 1.000193 | 417.6 K | 238.0 K | 1.733 |
| 2 | 1.000387 | 426.6 K | 257.7 K | 1.597 |
| 5 | 1.000967 | 456.2 K | 316.3 K | 1.321 |
| 11 | 1.002127 | 513.0 K | 412.0 K | 1.077 |
| 21 | 1.004060 | 599.5 K | 537.2 K | 0.925 |
| 101 | 1.025419 | 1240.1 K | 1273.0 K | 0.794 |
| 201 | 1.151525 | 4085.0 K | 4224.8 K | 0.906 |
| 401 | 5.990596 | 1493.3 K | 1489.7 K | 0.984 |

The ghost the boundary states now sits at **226.0 K, which is `base.inp`'s
`T_base` to the digit it is given in**; it sat at 413.7 K. The first cell falls
from 417.6 K to 238.0 K, the profiles cross near cell 21 and agree to 0.2 per
cent in the outer wind.

**Mass-loss rate**, from the wind flux at cell 401:

| | Mdot [g/s] | log10 |
|---|---|---|
| certified | 7.5096e+07 | 7.8756 |
| corrected | 7.3767e+07 | 7.8679 |

**-1.77 per cent**, and `REPRODUCE.md` quotes log10 Mdot = 7.88 from the
post-processing pass, which the certified value reproduces.

**So the movement is far above 0.1 per cent** -- 1.8 per cent in Mdot and up to
46 per cent in the base temperature -- **and the LHS 1140 b campaign has to be
re-solved.** The scale: 120 certified states in `LHS1140b/models/` (74 cases
plus 46 closure iterations, counted for item L18), at three outer passes and
about 400 s of solve each on this case, plus the post-processing and transit
passes `run_case.sh` adds. The transit equivalent widths were not re-synthesized
here: the base layer they are most sensitive to is exactly the layer that moved,
so they are part of the re-run and not separable from it.

## 6. Tests

`EXHALE_OBJDIR=build_L18`, the same binary with and without the two control
variables.

| suite | control | corrected |
|---|---|---|
| `adv_static_limit` | -- | 53 PASS, 0 FAIL |
| `carrier_retry` | -- | 143 PASS, 0 FAIL |
| `grid_and_gates` | 196 PASS, 4 FAIL | 196 PASS, 4 FAIL |
| `certification` | 84 PASS, 0 FAIL | 84 PASS, 0 FAIL |

The two `grid_and_gates` runs give **the same verdict on every one of the 200
assertions, in the same order** (`diff` of the two verdict lists is empty). The
four that fail are the ones items L18 and L20 already account for: three read
the stale pinned fixture `backup/regression/wasp_full_newton/IC/`, and
`outer_iteration_ending_is_the_stagnation_one` is the one-bit outcome that
flips on any last-bit change.

`adv_static_limit` and `carrier_retry` are the two suites outside
`grid_and_gates` whose drivers link `base_boundary` (`grep -rl base_boundary
src/tests`), and both pass with the corrected boundary.

Inside `grid_and_gates` the drivers that exercise this module are
`hydrostatic_residual` and `free_outflow_boundary`. The well-balanced rows --
the ones that state the boundary is exact on a hydrostatic atmosphere at rest,
where the branch weight is evaluated at `M = 0` and therefore at the centre of
the window -- read at the rounding floor with the control settings
(momentum 4.9e-14 to 5.5e-14, mass 4.6e-17 to 1.1e-16, energy 7.0e-17 to
1.6e-16, all against 1e-13) in all four reconstruction/flux combinations.

## 7. Regression

From copies of `backup/regression/`, one thread, the same binary with and
without the two control variables, so the comparison isolates this item.

| case | steps | `Hydro_ioniz` | `Ion_species` | `Hydro_ioniz_adv` | `Ion_species_adv` | base ghost T |
|---|---|---|---|---|---|---|
| `hydrostatic_column` | 300 | **0** | **0** | **0** | **0** | 1140.00 K both |
| `hp_front` | 100 | **0** | **0** | **0** | **0** | 1140.00 K both |
| `wasp_full` | 300 | **0** | **0** | -- | -- | 2357.28 K both |

**Byte-identical.** The correction is inert on these bases, and the reason is
the one the design note gives: their operating point sits far ABOVE the old
window with the right sign, so `w_rev` was 0 and stays 0 when the window
narrows. The change acts only where the discriminant had the wrong sign, which
is what it was made for.

**How the third row became a zero, and what it cost.** Without the guard of
section 4, `wasp_full` moved by 5.1e-03 in pressure and 2.0e-02 in the
composition at the ionization front (r = 1.19, cell 349) while its base was
unchanged -- ghost temperature identical to the digit, cell-1 density equal to
1 part in 1e8. The cause is the one the code now states: a COLD START's wind
window carries no wind, the isothermal initial condition being at rest there,
so during the first steps the mean flux is whatever the transient puts into it
and its sign says nothing about the base. The flux is now consulted only where
the face Mach number it implies is outside the handover window and therefore
selects a branch instead of sitting in the blend; below that the first interior
cell answers, which is what this boundary did before. With the guard the
300-step cold start is bit for bit what it was, and the LHS 1140 b branch is
unchanged (`M_branch = +4.024422e-07`, `w_rev = 0.000`).

`wasp_full` was run capped at 300 steps under both settings rather than to
convergence: the question this item asks of a hot Jupiter is whether its base
branch moves at all, the branch is evaluated on every call to the boundary, and
a cold start's first three hundred steps are where the wind window is least
able to answer -- so the cap tests the change where it is most likely to act,
not least. `wasp_he23off` was not run: it differs from `wasp_full` in the
helium triplet and not in the base condition, so it exercises the same boundary
at the same base Mach number.

## 8. Reproduce

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd $EX && make OBJDIR=build_L18 EXE=EXHALE_L18.x
M=$EX/LHS1140b/models;  L=$M/.L21

# section 3: the face flux, the cell-centred product and the branch, on the
# certified state.  The diagnostic writes both faces of every cell and the
# boundary reports what it did.
mkdir -p $L/probe/output
C=$M/atomic_scalar_gj1132_kzz1e9/HeH2.13
cp $C/input.inp $L/probe/;  cp $C/metals.inp $L/probe/ 2>/dev/null
cp $C/output/Hydro_ioniz.txt $L/probe/output/Hydro_ioniz_IC.txt
cp $C/output/Ion_species.txt $L/probe/output/Ion_species_IC.txt
sed -i 's/^Restart intent:.*/Restart intent: relaxation/' $L/probe/input.inp
cd $L/probe && OMP_NUM_THREADS=1 EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 $EX/EXHALE_L18.x
# the old behaviour, same binary:
OMP_NUM_THREADS=1 EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 \
  EXHALE_BASE_BRANCH_ON_CELL1=1 EXHALE_BASE_MACH_BLEND=1.0e-6 $EX/EXHALE_L18.x

# section 5: the re-solve, from the case's own certified state
sed -i 's/^Restart intent:.*/Restart intent: stationary/' input.inp
OMP_NUM_THREADS=8 $EX/EXHALE_L18.x                       # corrected
OMP_NUM_THREADS=8 EXHALE_BASE_BRANCH_ON_CELL1=1 \
  EXHALE_BASE_MACH_BLEND=1.0e-6 $EX/EXHALE_L18.x         # control
```

## 9. What this leaves for the user

1. **The LHS 1140 b campaign has to be re-solved.** 120 certified states, about
   400 s of solve each on this planet plus the post-processing and transit
   passes. Nothing in `LHS1140b/models/` was rewritten by this item; every run
   above is in `LHS1140b/models/.L21/`.
2. **Every number the campaign has published about the BASE is affected**, and
   in the direction that matters: the base of those states was 417 K where
   `base.inp` states 226 K, and the transit equivalent widths are formed in
   exactly the layer that moves.
3. **The hot-Uranus and hot-Jupiter bases are untouched** (section 7), so the
   goldens and the regression matrix are unaffected and nothing needs
   refreshing.


---

## 10. The Codex review of 2026-09-15, applied

Three corrections, all verified to leave the LHS 1140 b answer bit for bit
where it was (section 10.4).

### 10.1 The handover is C1, and it blends the WEIGHTS

The first form of this item took the wind's discriminant through a threshold,
`if |M_wind| > blend then M_branch = M_wind else M_i`. That is a jump between
two numbers OF OPPOSITE SIGN -- on this base `M_i = -8.1e-07` and
`M_wind = +4.0e-07` -- so the entropy source of the face, and with it the ghost
density, stepped discontinuously and the residual carried a jump the Newton
solve cannot differentiate.

The handover is now the same cubic smoothstep the branch itself uses, read as a
function of `|M_wind|/blend`: 0 at 0, 1/2 at the blend width, 1 at twice it.
Its derivative vanishes at both ends, so the `|.|` puts no kink at
`M_wind = 0` either.

**And what is blended is the weight, not the Mach number.** Blending the Machs
was tried first and measured: `w_rev` inherits `(1 - s) M_i`, and `|M_i|` can be
hundreds of blend widths at a base inside a collocated mode, so the blended
Mach sweeps hundreds of widths while `s` moves by a per cent. Measured on the
synthetic sweep of the new test, with `M_i = -547` blend the blended Mach
crossed zero between two samples 3 per cent apart in the wind speed and the
face density jumped by 0.398 -- continuous, but with a derivative of order
`|M_i|/blend`, which for a Newton solve is nearly as bad as the jump it
replaced. Each discriminant's WEIGHT is bounded in [0,1] by construction, so
blending those is bounded as well and is still C1.

### 10.2 The wind window is read on the stationary route only -- WITHDRAWN

**This subsection records a decision that was made on 2026-09-15 and reversed
the same day. It is kept for the record; what the code does is section 11.**

Reading a flux from `r >= r_flux` is a STATIONARY statement: it is the flux
through the base face only because a steady wind carries one flux. Using it as a
boundary condition in physical time would let a distant cell set the base of a
transient whose base and wind need not carry the same flux. The marching path
therefore keeps the local face and the stationary residual, its solve and its
certification use the wind (`set_base_branch_stationary_route`, called at the
three entries of the stationary route and at the `EXHALE_RESIDUAL` diagnostic).

Two costs, recorded at the place in the code where the question arises rather
than left implicit:

- the marched operator and the stationary residual are then not the same
  discrete operator at the base, so a state marched to rest and a state the
  stationary route certifies are fixed points of two boundary conditions, not
  one. Every certified LHS 1140 b state is reached through the stationary
  route, so its states are states of the boundary that certified them; a
  marched state is not claimed to be;
- on the stationary route the residual of the base cells depends on cells at
  `r >= r_flux`, which the banded preconditioner does not contain. The Krylov
  products see it through the finite-difference action; the preconditioner does
  not, and that is one more term in the model error item L17 already measures
  at a factor 27.

### 10.3 Four new assertions

`src/tests/grid_and_gates/base_branch_discriminant.f90`, run by the suite as
`base_branch`. It builds a column and calls the production
`base_boundary_states` on it; nothing is re-implemented.

| assertion | measured |
|---|---|
| `base_branch_wind_wins_where_the_window_carries_flux` (renamed, section 11.2) | `w_rev` = 0 with a first cell moving inward at 1000 blend widths inside a wind of +100 |
| `base_branch_local_face_wins_where_it_does_not` (renamed, section 11.2) | `w_rev` = 1 on the same state with no flux in the window; the two face densities differ by exactly the factor 2 of the two isentropes |
| `base_branch_at_rest_face_velocity_is_zero` | 0.000e+00 |
| `base_branch_at_rest_is_reproducible` | 0.000e+00 |
| `base_branch_handover_is_continuous` | the largest neighbouring change of the face density over a sweep through both thresholds falls from 2.405e-02 at 401 samples to 1.204e-02 at 801, **ratio 0.501** -- what a continuous function gives, where the threshold form measured 1.000 |
| `base_branch_independent_of_the_flux_window` | 0.000e+00 between `j_flux = N/4` and `3N/4` on a state of uniform mass flux |

### 10.4 Nothing of the LHS 1140 b answer moved

The certified state of section 5, re-entered with `Restart intent: stationary
evaluate`, is **bit for bit identical** between the campaign binary
(`db87b88d1ce53facf1d61084fa535ca5`) and this build
(`e8188190f74a8067a7ebb0bcdbbde8d7`) in both state files, and re-certifies with
the same rows (mass 2.175e-09, momentum 5.417e-13, energy 2.126e-08). It has to
be: `|M_wind| = 40` blend widths there, so the weight of the wind is exactly 1
and the threshold form and the smooth form agree.


### 10.5 Regression after the review

The three cases the item names, run from copies with the review build against
the pre-L21 control, one thread, same step caps:

| case | steps | `Hydro_ioniz` | `Ion_species` | `Hydro_ioniz_adv` | `Ion_species_adv` |
|---|---|---|---|---|---|
| `wasp_full` | 300 | **0** | **0** | **0** | **0** |
| `hydrostatic_column` | 300 | **0** | **0** | **0** | **0** |
| `hp_front` | 100 | 1.102e-04 | 1.533e-07 | 5.552e-04 | 3.743e-06 |

`wasp_full` and `hydrostatic_column` are byte-identical, which is what gating
the wind window to the stationary route is for: the marching path is the
boundary it always was.

**`hp_front`'s 1.1e-04 is not this item's.** The review build in the EXACT
control configuration of the base boundary
(`EXHALE_BASE_BRANCH_ON_CELL1=1 EXHALE_BASE_MACH_BLEND=1.0e-6`, so the
boundary is bit for bit the old one) reads the same 1.102e-04, which puts the
difference outside `base_boundary.f90`; its loader line is identical in both
(4.862e-16, the conserved density taken from its column); and `hp_front` is the
only one of the three that is molecular (`Molecular chemistry`, carrier and
ionization transport), while `mol_rates.f90`, `molecular_reaction_heat.f90` and
`steady_newton.f90` were edited in this tree by another worker between the
campaign binary's build and this one. A clean control build of "the current
tree minus this item" was not made, because the tree is being edited
concurrently; that is stated here rather than left as an assumption.


### 10.6 The suites after the review

| suite | review build |
|---|---|
| `certification` | 84 PASS, 0 FAIL |
| `run_mode` | 31 PASS, 0 FAIL |
| `grid_and_gates` | **202 PASS, 4 FAIL** (196 + the six new `base_branch` rows) |

The four that fail are the ones items L18 and L20 already account for: three read
the stale pinned fixture `backup/regression/wasp_full_newton/IC`, and
`outer_iteration_ending_is_the_stagnation_one` is the one-bit outcome that flips
on any last-bit change.

**One row is worth more than its count.** `stationary_residual_equals_the
_diagnostic` compares the same state measured by three routes -- the evaluation
of the state as loaded, the `EXHALE_RESIDUAL` diagnostic, and the first rows the
stationary solve judges -- and it reads

```
  the evaluation of the state as loaded      ||R|| = 1.6737E+00
  EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0       ||R|| = 1.6737E+00
  the first rows the stationary solve judged ||R|| = 1.674E+00
```

so the three see ONE operator. That is the check the stationary-route gate of
10.2 could have broken and does not: the wind window is read at the evaluate, at
the diagnostic and inside the solve, and nowhere else.

## 11. The second Codex review of 2026-09-15: one boundary, not two

Section 10.2 above is **withdrawn**. Keying the base boundary on which route
was evaluating was wrong, and the measurement below is why.

### 11.1 The two boundaries measured on one certified state

`LHS1140b/models/.L21/resolve_new`, the certified
`atomic_scalar_gj1132_kzz1e9/HeH2.13` state of section 5, re-entered with
`Restart intent: stationary evaluate` under each boundary in turn. The state
is the same file in both columns; only the boundary differs.

| quantity at the base | flux-keyed face | local face |
|---|---|---|
| base face mass flux `Phi(0) r^2` | +6.306542E-07 (= +1.000 `F_wind`) | **-7.493857E-07 (= -1.188 `F_wind`)** |
| wind flux at cell 400 | +6.306542E-07 | +6.306542E-07 |
| face / ghost density | 9.505841E+13 (`rho_res`) | 8.810973E+13 (`rho_rev`), ratio 0.92690 |
| `T_res` / `T_i` at the face | 223.195 K / 240.807 K | the same |
| branch weight `w_rev` | 0.000 | 1.000 |
| `R_mass(1)` | 7.091301E-12 | -6.798026E-03 |
| `R_mom(1)` | -2.708833E-12 | -4.882242E-02 |
| `R_energy(1)` | 7.326027E-12 | -6.768184E-03 |
| certification rows, mass / momentum / energy | 2.174586E-09 / 5.417328E-13 / 2.125677E-08 | 1.305538E+00 / 1.702142E-04 / 1.053417E+00 |

`n(1)`, `v(1)` and `T(1)` are identical; only the two cells at `r <= 1.0004`
differ by more than 1e-6. **`R_stat(U*) = 0` said nothing whatever about
`R_time(U*)`**: nine decades separate them on one state, and the certified
state is not a steady state of the marched operator at all. A state is a
steady state of ONE operator, and two boundaries meant two operators.

### 11.2 What the code now does

`base_branch_stationary_route`, `set_base_branch_stationary_route` and its four
call sites are removed. The branch asks only whether the wind window carries a
flux that says anything:

```fortran
w_i    = characteristic_branch_weight(M_i/base_face_mach_blend)
M_wind = 0.0d0
s_wind = 0.0d0
w_rev  = w_i
if (base_branch_on_wind_flux .and. have_F) then
   M_wind = F_wind/(rho_i*rb*rb)/c_i
   s_wind = characteristic_branch_weight(1.0d0 - abs(M_wind)/base_face_mach_blend)
   w_rev  = s_wind*characteristic_branch_weight(M_wind/base_face_mach_blend)   &
          + (1.0d0 - s_wind)*w_i
endif
```

`have_F` is false where the window is empty or its mean is not a usable
number, and `s_wind` vanishes with `|M_wind|`; where the wind says nothing the
first interior cell answers, which is the cold start and the early transient.
`EXHALE_BASE_BRANCH_ON_CELL1=1` still restores the pre-L21 discriminant.

The two tests that named the routes are renamed for what they measure:
`base_branch_wind_wins_where_the_window_carries_flux` and
`base_branch_local_face_wins_where_it_does_not`. All six assertions of
`base_branch` pass.

The preconditioner's boundary block stays a PROPOSAL and is recorded as one in
the source: the residual of the base cells depends on cells at `r >= r_flux`
which the banded preconditioner does not contain, the Krylov products see that
dependence through the finite-difference action and the preconditioner does
not, and that is one more term in the model error item L17 measures at a factor
27. Nothing is built for it here.

The certification label is written into
`docs/restart_contract_design_20260909.md` section 7: certification is
certification against the stationary operator, and stays so until a state is
measured to be a fixed point of both paths.

### 11.3 The certified answer does not move, and the regression trio does

The fiducial of section 5, re-evaluated with the one-boundary build
(`6463f0fdf181bd1859d59da3226e7819`, and again with the comment added at the
blend width, `142ed07857266ea832f9389ecc210d59`) against the campaign binary
(`db87b88d1ce53facf1d61084fa535ca5`), is **bit for bit identical** in both
state files and re-certifies with the same rows (mass 2.175E-09, momentum
5.417E-13, energy 2.126E-08, CERTIFIED). It has to be: `|M_wind|` is 30 blend
widths there, so `s_wind` is exactly 1 and the wind decides under either
keying.

The regression trio, control (`922dec0fa183c2611d1bfb33c218a6cd`, the tree
before this item) against measured, single-threaded, `REGRESSION_REL_TOL=0`,
at the step counts the review asked for:

| case | steps | control vs one boundary | control vs SAME BINARY with the pre-L21 discriminant | the two measured runs |
|---|---|---|---|---|
| `hydrostatic_column` | 300 | 9.6e-08 `Hydro`, 1.0 `Ion` | 9.6e-08 `Hydro`, 1.0 `Ion` | **IDENTICAL** |
| `hp_front` | 100 | 2.5e-01 `Hydro`, 1.0 `Ion` | 2.5e-01 `Hydro`, 1.0 `Ion` | **IDENTICAL** |
| `wasp_full` | 300 | 3.954e-04 `Hydro`, 1.580e-03 `Ion` | 2.9e-11 `Hydro`, 7.9e-05 `Ion` | 3.954e-04 / 1.580e-03 |

So `hydrostatic_column` and `hp_front` do not move under this item AT ALL --
their difference from the control binary is other source changes in this
working tree, not the boundary -- and `wasp_full` does move, by 0.040 per cent
in the pressure and 0.158 per cent in an O I density of 1.8e-06 at
r = 1.194 R_p. **That movement is removed in section 12**; the row is kept
because the measurement is what identified its cause. The `Ion` entry of 1.0 for `hydrostatic_column` is an underflowed density
(1.16e-170 against 1.26e-180); for `hp_front` it is the H3+ density at
r = 1.074 R_p (6.49e-02 against 1.08e-05), a column the molecular chemistry of
item L7f is being changed under, and it is identical between the two boundary
settings of one binary.

### 11.4 Why `wasp_full` moves, and what it says about the blend width

The expectation was zero, on the argument that a cold start's wind window
carries no wind and the guard hands the branch to the local face there.
MEASURED, that is not quite true, and the reason is worth recording because it
is a defect of CALIBRATION and not of the construction.

Factoring the two settings apart on the same binary, one marching step of
`wasp_full`:

| base boundary | blend width | max rel (rho, v, p, T) after ONE step |
|---|---|---|
| local face | 1e-6 | reference |
| local face | 1e-8 | 0.000e+00 |
| wind window | 1e-6 | 4.451e-08 |
| wind window | 1e-8 | **1.829e-04** |

The mover is the BLEND WIDTH INSIDE THE WIND GUARD, not the discriminant.
Measured on that state: the cold start's wind window carries
`F_wind = 1.216e-14` in code units, which is `M_wind = +3.70e-09` at the base,
while the first cell is at `M_i = +1.15e-02`. With `base_face_mach_blend`
lowered to 1e-8 by this item, `|M_wind|` is 0.37 blend widths, so
`s_wind = 0.090` -- the window that carries no wind is given a 9 per cent say
-- and, being inside the blend window, its own reversal weight is 0.235, so
`w_rev = 0.021` where the local face gives 0. Two per cent of the reversal
datum at a face whose two isentropes differ by 8 per cent is the 1.8e-04, and
it is there from the first step.

**One symbol is doing two jobs.** `base_face_mach_blend` is (a) the width over
which the reversal handover of a face Mach number is smoothed and (b) the scale
below which the wind window is deemed silent. These are different physical
statements. It was lowered from 1e-6 to 1e-8 for (b) alone: at the LHS 1140 b
base of the state this item was diagnosed on, `M_wind = +4.0e-07` and
`M_i = -8.1e-07` (on the state re-solved under the fix they read +2.99e-07 and
-5.61e-06), so at 1e-6 the window would have had `s_wind = 0.10`
and the local face's `w_rev = 0.974` would have carried the branch -- the very
defect this item fixes. Job (a) needed no change there: `|M_i| = 8.1e-07` is
0.81 blend widths at 1e-6 and `w_i = 0.974`, which is the same verdict as the
1.000 it gives at 1e-8.

The separation available to a single number is only a factor 81 between the
certified LHS wind (2.99e-07) and this cold start (3.70e-09), which is not
much.
PROPOSED here, and **DONE in section 12** after the user's decision of
2026-09-16: give the wind-window guard its own statement, scale-free, as the
window's own flux spread. Section 12 carries the construction, the threshold's
measurement and the four results; with it `wasp_full` returns to identical.

### 11.5 What was NOT measured

`wasp_full` was run to 300 steps, which is the probe the review asked for and
not the converged case of the matrix (12247 steps). The converged case was not
re-run; its `golden/` entry is in any event stale by three decades and the
matrix is read control-against-measured, not against `golden/`.

## 12. The window's standing is read off the window, not off a Mach number (user decision, 2026-09-16)

Section 11.4 measured that `base_face_mach_blend` was serving two statements
at once and proposed separating them. The decision of 2026-09-16 separates
them and makes the second one SCALE-FREE.

### 12.1 What the branch now asks

The question put to the wind window is not "is its flux large" -- which needs
a scale and had been borrowing the reversal handover's -- but **"is its flux
ONE flux"**. A window that carries a wind carries the same `rho v r^2` at
every altitude in it, and that is a property of the state with no scale of its
own to borrow. So the window's weight in the branch is read off the window's
own relative spread,

```
du_window = ( max_j - min_j ) / |mean_j|  of  rho v r^2  over  j = j_flux..N ,
```

the same functional the convergence gate is defined on
(`flux_spread_of_state`), computed in `wind_window_mass_flux` on the same pass
that forms the mean, from the same cells, with no new state dependence and no
second traversal.

```fortran
s_wind = characteristic_branch_weight(                          &
            (2.0d0*du_window - 3.0d0*base_wind_window_spread)   &
            /base_wind_window_spread)
```

which is the same cubic smoothstep the branch itself uses, equal to **1 at or
below `base_wind_window_spread` and 0 at or above twice it**, with a vanishing
derivative at both ends. `base_face_mach_blend` now carries ONE statement, the
width of the reversal handover of a face Mach number, and keeps 1e-8.

### 12.2 The threshold is the measured separation

Chosen the way `flux_spread_th_default` was (`parameters.f90`): the geometric
middle of the gap between the two populations, to one digit, fitted to no
single state. MEASURED 2026-09-16 on this quantity, over this window, on the
physical cells:

| population | `du_window` |
|---|---|
| certified stationary states, 147 of them (`LHS1140b/models` and `models_20260914_preL21`) | 1.51e-04 to **2.95e-03**, median 2.67e-04 |
| a converged Newton state (`wasp_full_newton`'s pinned pair) | 8.93e-04 |
| *the gap, a factor 25 wide* | |
| a partly relaxed reload (`carrier_model_a_newton`'s IC) | **7.50e-02** |
| `hp_front` | 7.99e-01 |
| `carrier_elem_newton` | 1.02 |
| `atomic_elem_newton` | 3.51 |
| `wasp_full` cold start, 1 step to 300 steps | 4.35 to 11.5 |
| `hydrostatic_column`, 300 steps | 12.6 |

`sqrt(2.95e-03 * 7.50e-02) = 1.5e-02`, and **1e-2** to one digit: 3.4 times
above the loosest certified state and 7.5 times below the tightest state that
is not a wind. Every fixture in the matrix saturates the weight on one side or
the other; none sits inside the handover.

**One deviation from the decision as stated, and the measurement that forces
it.** The decision asked for `s_wind = 1` at `du_window <= 0`. The
cell-centred product is a reconstruction of the state, so its spread over the
window has a floor set by the reconstruction and not by the wind -- the
certified states sit at 1.51e-04 to 2.95e-03 and **no state reaches zero**. A
weight reaching one only at `du_window = 0` would therefore be below one on
every state, and a converged wind would carry a few parts in ten thousand of
the branch it does not belong to, which is this item's own defect in a smaller
size; the certified fiducial would not have re-evaluated bitwise. The full
weight is therefore reached at the threshold itself, and the zero side is
exactly as decided, at twice it.

(The decision's `du ~ 1e-10` for a stationary wind is the FACE-flux gate. That
number is not available where the boundary runs: `flux_spread_of_state` reads
the Riemann face fluxes cached by `assemble_residual`, and the boundary is
evaluated BEFORE those fluxes are formed -- it is what feeds them. The
cell-centred product over the same window is what the routine can read, it is
what it already averages, and its own floor is the 1.5e-04 above.)

### 12.3 The four measurements

Binary `EXHALE_L21b.x`, md5 `2966d426e782f5d36e6fbcde34ffb2c0` (the `22e5bfa0ca8cb7687e89fb7e406f0655` build
differs from it by one comment block and gives the same four results), one
thread.

| # | what | result |
|---|---|---|
| 1 | `wasp_full`, 300 steps, cold start, against the SAME binary carrying the pre-L21 boundary | **`Hydro_ioniz.txt` and `Ion_species.txt` IDENTICAL** (`du_window` 4.35 at step 1 and 11.5 at step 300, so `s_wind = 0` throughout and the local face answers, as it did before this item). The 3.954e-04 of section 11.3 is gone. |
| 2 | the certified fiducial, `Restart intent: stationary evaluate`, against the campaign binary `db87b88d1ce53facf1d61084fa535ca5` | **both files IDENTICAL**, CERTIFIED, rows 2.175E-09 / 5.417E-13 / 2.126E-08. The diagnostic reads `du_window = 2.458e-04`, `s_wind = 1.000000`, `w_rev = 0.000`. |
| 3 | `hydrostatic_column` 300 steps and `hp_front` 100 steps, same comparison | **IDENTICAL**, both files, both cases (`du_window` 12.6 and 0.80) |
| 4 | the L14 fixture `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13`, re-solved from its own seed at `EXHALE_PTC_DTAU0=1.0e8`, one thread, against the Mach-keyed build | the seed reads `du_window = 2.52e-03`, so `s_wind = 1` under the new keying and under the old one alike, and the branch cannot differ. MEASURED, the solve now complete: both reach **CERTIFIED at outer pass 35**, all 35 outer-pass lines are identical field for field, both written state files are **bitwise identical**, and the two `run.log`s differ in nothing but their wall times. |

The `base_branch` suite is restated on the new discriminant and all six
assertions pass:

| assertion | measured |
|---|---|
| `base_branch_wind_wins_where_the_window_is_one_wind` | `w_rev` = 0 on a window of uniform mass flux (`du_window` = 0) with the first cell moving inward at 1000 blend widths |
| `base_branch_local_face_wins_where_it_is_not` | `w_rev` = 1 on the SAME mean flux scattered to `du_window` = 3 thresholds; the two face densities differ by exactly the factor 2 of the two isentropes |
| `base_branch_at_rest_face_velocity_is_zero` | 0.000e+00 |
| `base_branch_at_rest_is_reproducible` | 0.000e+00 |
| `base_branch_handover_is_continuous` | sweeping `du_window` from half to two and a half thresholds, the largest neighbouring change of the face density falls from 3.746e-03 at 401 samples to 1.873e-03 at 801, **ratio 0.5000118**; the handover's midpoint is at `du_window` = 1.50e-02 where `s` = 0.4998 |
| `base_branch_independent_of_the_flux_window` | 0.000e+00 between `j_flux = N/4` and `3N/4` on a state of uniform mass flux |

The test column is now built with `v(j) = v_out (r_1/r_j)^2` so that the
window really does carry one mass flux, which is what the assertion about the
window's start had been asserting without providing.

Suites on this build: `certification` 84 PASS / 0 FAIL, `run_mode` 31 / 0,
`grid_and_gates` 202 PASS / 4 FAIL -- the four being the pre-existing failures
of the stale `wasp_full_newton/IC` fixture and the stagnation-ending row,
reproduced identically on `EXHALE_L14.x` of 2026-09-14 and recorded in
`lhs1140b_stationary_L7b_20260913.md` and `Update_EXHALE_stage2.md` section 11.

### 12.4 What is left of section 11.4

Nothing to act on. `base_face_mach_blend` keeps 1e-8 and one job; the wind
window's standing has its own scale-free statement and its own threshold with
its own measurement. `EXHALE_BASE_WIND_SAY_ON_MACH=1` restores the Mach-keyed
weight and `EXHALE_BASE_WIND_SPREAD=<value>` moves the threshold, both for
control experiments only.

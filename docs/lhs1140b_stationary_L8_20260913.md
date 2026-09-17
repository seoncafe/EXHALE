# LHS 1140 b, item L8: what the outer boundary of the subsonic wind states, and how much of the answer it sets

Item L8 of `docs/PLAN_20260913_lhs_stationary.md`, raised by section 6 of
`docs/PLAN_20260913_lhs_stationary_review.md`: the free-outflow continuation
of `states/Apply_BC.f90` is derived for an outflow with no incoming
characteristic, and the LHS 1140 b wind is subsonic at 30 R_p, so `v - c`
enters and the continuation IS the external condition. This is the external
condition, stated; the sensitivity of the mass-loss rate, the He 2^3S column
and the He I 10830 line to the outer radius, measured; and the verdict on
whether the certified state is a breeze selected by the environment or the
subsonic part of a transonic escape.

Every number is MEASURED with the binary of 2026-09-13 (md5
`97e10317a710b9ccc63addbedde3586a`, `EXHALE.x` as delivered, not rebuilt)
unless marked READ. **No file in `src/` other than
`src/utils/map_state_to_grid.py` was changed**, and that one only to let a
solution on a 30 R_p grid seed a run on a wider one.

---

## 1. Verdict

**The certified 30 R_p state is the subsonic part of a TRANSONIC escape, not
a breeze. Its critical point sits at 40.06 R_p, just outside the domain, so
the 30 R_p run cannot see it; the isothermal hydrostatic ghost is what stands
in its place, and it stands in well enough that the mass flux is unmoved.**

1. **What the continuation imposes.** At the outer face of the certified
   state the flow is outward at Mach 0.521, so of the three characteristic
   speeds `v-c = -2.30e5`, `v = 1.24e5`, `v+c = 4.78e5 cm/s`, one points
   inward and the boundary owes the domain exactly one condition. The ghost
   rule supplies it as a state of rest: `p_g/p_N = rho_g/rho_N =
   exp[-(phi_g - phi_N) rho_N/p_N]` at fixed `p/rho`, with `v_g = v_N`. That
   is the assertion **"above r_out the gas is a static isothermal atmosphere
   at the temperature of cell N"**, whose limit is a finite external
   pressure `p_inf = p_N exp(-GM rho_N/(r_N p_N)) = 0.134 p_N =
   2.04e-10 dyn/cm^2`. Equivalently, in the steady momentum equation it sets
   the ram term `rho v dv/dr` to zero at the face and leaves
   `dp/dr = -rho GM/r^2`.
2. **The interior state meets that condition rather than fighting it.** In
   the last physical cell the measured ram term is 0.22 percent of the
   weight and `dp/dr + rho g + rho v dv/dr` closes to 0.13 percent of the
   weight, so the solution has already gone hydrostatic where the ghost says
   it is hydrostatic. Taken alone that is consistent with a breeze, and the
   wind-equation numerator (section 2) does vanish at 29.11 R_p, one fifth
   of a cell outside the last cell center, while the flow is at Mach 0.52 --
   the textbook signature of a breeze. **That signature is manufactured by
   the boundary**: the two ghost rules together are the statement `N = 0` at
   the outer face.
3. **Moving the boundary settles it.** At `Outer radius 45` and `60 R_p`,
   from the same certified state as seed, the same binary and the same
   input, the solution crosses its critical point at **40.06 R_p** (r45) and
   **40.08 R_p** (r60) and leaves the domain supersonic at Mach 1.065 and
   1.292. `N = 0` and `v = c` at the same radius to 0.02 R_p, so it is a
   genuine critical point and not a turning point. The two extended
   solutions agree with each other to better than 0.1 percent at every
   radius: once the boundary is outside the sonic point the answer stops
   depending on where it is, which is what the free-outflow rule is built
   for.
4. **What that costs the campaign, and what it does not.** `log10 Mdot` is
   **7.8737 / 7.8740 / 7.8735** at 30 / 45 / 60 R_p -- a spread of
   5e-4 dex, 0.12 percent, over a factor two in domain size. The He 2^3S
   radial column moves by 0.23 percent and the He I 10830 red-pair
   equivalent width by up to 1.7 percent. **The profiles above 10 R_p are
   another matter**: at 28 R_p the 30 R_p solution is 25 percent too hot,
   29 percent too dense and 22 percent too slow against the converged
   answer (table T3).
5. **None of this is a validated continuum result above 7.6 R_p.** The
   exobase (`Kn_bulk = 1`) of the extended solutions is at **29.0 R_p**,
   BELOW the critical point at 40.06 R_p, and `Kn_bulk` passes 0.1 at
   7.65 R_p. The continuum model therefore fixes its own topology in gas it
   cannot describe. The mass-flux insensitivity of item 4 is a statement
   about the discrete problem, not a physical validation.
6. **There is no input key for an alternative outer prescription**, so the
   variant of the ghost rule at a fixed 30 R_p could not be run. `Apply_BC_W`
   calls `free_outflow_ghost` unconditionally (`Apply_BC.f90:182-185`), the
   routine branches only on its positivity fallback, and the `known_keys`
   list of `input_read.f90:125-164` carries no outer-boundary key:
   `Domain mode` and `Outer radius` place the boundary, and nothing in the
   input states what is beyond it.

---

## 2. The criterion: breeze or transonic escape

For steady spherical flow with `c^2 = gamma p/rho` and `q = (Gamma-Lambda)/rho`
the net heating per unit mass, continuity, momentum and energy combine into

```
(v^2 - c^2) (1/v) dv/dr = N(r),    N(r) = 2c^2/r - GM/r^2 - (gamma-1) q/v
```

with `gamma = 5/3` and `phi = -GM/r` (spherical mode: `grav_field.f90:16-19`).
The criterion this report uses:

- a **transonic escape** has `N = 0` and `v = c` at the SAME radius; there
  `dv/dr` is finite and the solution crosses to supersonic, and the mass flux
  is fixed by that crossing, not by anything downstream;
- a **breeze** has `N = 0` at a radius where `v < c`; there `dv/dr = 0` and
  beyond it `N > 0` with `v^2 - c^2 < 0` decelerates the flow, `v -> 0` and
  `p -> p_inf > 0`, so the mass flux is what the external pressure selects;
- a solution whose `N` is still negative at the edge of its domain has not
  decided, and what decides it is whatever is imposed there.

`N` is measured from the solved state's own columns, and the identity above
reproduces the measured `dv/dr` of the certified state to 0.3-1 percent over
the outer 5 R_p, which is how the evaluation is checked.

**T1. The certified 30 R_p state near its boundary** (raw solved state,
`output/Hydro_ioniz.txt`)

| r [R_p] | Mach | 2c^2/r | GM/r^2 | (g-1)q/v | N/(2c^2/r) | dv/dr measured | from N |
|---|---|---|---|---|---|---|---|
| 25.813 | 0.5131 | 3.922 | 2.636 | 1.551 | -0.0676 | 7.735e-7 | 7.726e-7 |
| 27.145 | 0.5184 | 3.455 | 2.325 | 1.348 | -0.0397 | 4.390e-7 | 4.381e-7 |
| 28.073 | 0.5204 | 3.578 | 2.229 | 1.425 | -0.0210 | 2.270e-7 | 2.257e-7 |
| 28.549 | 0.5209 | 3.516 | 2.155 | 1.402 | -0.0117 | 1.255e-7 | 1.239e-7 |
| 29.031 | 0.5212 | 3.457 | 2.084 | 1.381 | -0.0025 | 3.760e-8 | 2.608e-8 |

A straight line through the last ten cells puts `N = 0` at **29.11 R_p**,
0.08 R_p past the last cell center, with Mach 0.521 there. Read alone, that
is a breeze. It is the boundary condition: `v_g = v_N` makes the limited
velocity slope of cell N zero under PLM and the hydrostatic `p_g` makes the
pressure gradient balance the weight, which is `N = 0` written out.

---

## 3. The experiment

### 3.1 Recipe

The certified case
`LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH1.60` (its `REPRODUCE.md` is
the prototype) at `Outer radius [R_p]` 30, 45 and 60, everything else
unchanged, each solved by `models/run_case.sh` with the certified state as
seed. Run directories `LHS1140b/models/.L8/`, outside the campaign tree.

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
L8=$EX/LHS1140b/models/.L8
CERT=$EX/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH1.60

# 1. the target grid of each outer radius: a cold three-step run of the case
for R in 30 45 60; do
   mkdir -p $L8/grid_r$R
   sed -e 's/^Load IC?.*/Load IC? False/' -e '/^Restart intent:/d' -e '/^Solver:/d' \
       -e "s/^Outer radius \[R_p\]:.*/Outer radius [R_p]: ${R}.0/" \
       $CERT/input.inp > $L8/grid_r$R/input.inp
   ( cd $L8/grid_r$R && mkdir -p output && OMP_NUM_THREADS=2 EXHALE_MAXSTEPS=3 $EX/EXHALE.x > gridgen.log 2>&1 )
done

# 2. the seed: the certified 30 R_p state, continued above its last cell
for R in 30 45 60; do
   mkdir -p $L8/r$R/output
   sed -e "s/^Outer radius \[R_p\]:.*/Outer radius [R_p]: ${R}.0/" $CERT/input.inp > $L8/r$R/input.inp
   python3 $EX/src/utils/map_state_to_grid.py $CERT/output \
      $L8/grid_r$R/output/Hydro_ioniz.txt $L8/r$R/output --ic \
      --extrapolate-beyond 29.031224061195065
done

# 3. wind, advection-corrected profiles, transit spectrum, record
cd $EX/LHS1140b/models
for R in 30 45 60; do NOSEED=1 OMP_NUM_THREADS=8 ./run_case.sh .L8/r$R; done
```

The `grid_r30` grid is byte-identical to the campaign's
`models/current_grid_Hydro_ioniz.txt`, which is what makes step 1 a
trustworthy way to produce the other two.

### 3.2 The seed above 29.03 R_p

`src/utils/map_state_to_grid.py` refuses a physical target cell outside the
source's support, so the 45 and 60 R_p grids had no seed. The option
`--extrapolate-beyond <r/R_p>` was added: above `r` the state is not
interpolated but continued -- constant `T`, isothermal hydrostatic `rho` and
`p` (the same expression `free_outflow_ghost` writes into the ghosts, so the
seed and the boundary state one stratification), constant mass flux `v`, and
frozen composition ratios `n_s/n_H`. The gravity is the one the source state
carries: `g = -(1/rho) dp/dr - v dv/dr` from one-sided differences of its
last interior rows, scaled as `1/r^2`. The measured value,
`GM/R_p = 1.977931e12 (cm/s)^2 R_p`, is 0.13 percent from `G M_p / R_p`
computed from the input file's mass and radius -- the check that the source
is steady enough for its own momentum balance to name its gravity. The
`# mapped:` line of the seed records the anchor, that value and how many rows
were filled that way.

The `# grid` header line is now taken from the target file instead of being
copied from the source, because `load_IC` refuses a state whose grid field is
not the run's (`load_IC.f90:152-163`); without it the 45 and 60 R_p runs
stopped at `(load_IC) ERROR: metadata field "grid" differs`. A target on the
same grid leaves that line unchanged.

Tests: `src/tests/state_mapper/run.sh` PASSED before (17 assertions) and
PASSED after, with the same verdicts; and the campaign mapping
(`archive_20260830/.../a_heh/output` onto `current_grid_Hydro_ioniz.txt`)
is byte-identical between the entry text and the changed tool, so the
existing behavior did not move.

### 3.3 What the wider grids cost in resolution

`Grid cells: 500` is fixed, so a wider domain is also a coarser one, and this
is the confound of the experiment. Of the 500 cells, 499 / 477 / 463 lie
below 29.03 R_p at 30 / 45 / 60 R_p, and the cell width at 10 R_p grows by
5.6 and 9.7 percent. The interior discretization therefore changes by under
10 percent while the boundary moves by a factor of two, and the two extended
solutions -- which differ from each other in resolution by the same kind of
amount -- agree to 0.1 percent, so the resolution change is not what the
tables below are showing.

---

## 4. Results

**T2. The three solutions** (raw solved state for the boundary quantities,
`_adv` profiles for `Kn` and the line, as `collisional_validity.py` and
`EXHALE_transit.py` read them)

| `Outer radius` [R_p] | 30 | 45 | 60 |
|---|---|---|---|
| last physical cell [R_p] | 29.031 | 43.445 | 57.835 |
| solver, last residual norm | info = 0, 3.17e-8 | info = 0, 2.15e-8 | info = 0, 2.34e-8 |
| certification | CERTIFIED | CERTIFIED | CERTIFIED |
| wall clock, 8 threads | 0m25s | 13m06s | 8m58s |
| Mach at the outer face (raw) | 0.5212 | 1.0651 | 1.2918 |
| Mach at the outer face (`_adv`) | 0.3506 | 0.7854 | 0.6352 |
| `T` at the outer face (raw) [K] | 429.2 | 252.8 | 210.1 |
| critical point (raw) [R_p] | none in domain | 40.054 | 40.081 |
| `N = 0` (raw) [R_p] | 29.11 (extrapolated) | 40.071 | 40.064 |
| `Kn_bulk` at the outer face | 0.645 | 2.15 | 2.76 |
| `Kn(H I)` at the outer face | 0.916 | 3.06 | 4.01 |
| exobase, `Kn_bulk = 1` [R_p] | above the domain | 29.000 | 29.000 |
| `log10 Mdot` [g/s], outer cell | 7.8737 | 7.8740 | 7.8735 |
| `log10 Mdot`, median over r > 10 | 7.8737 | 7.8738 | 7.8739 |
| He 2^3S radial column [cm^-2] | 6.7556e10 | 6.7409e10 | 6.7562e10 |
| He I 10830 red-pair EW [%A] | 1.1702 | 1.1504 | 1.1537 |
| red-pair depth [%] | 4.3260 | 4.2595 | 4.2714 |
| FWHM [A] | 0.2538 | 0.2538 | 0.2538 |

Relative to 60 R_p: `log10 Mdot` +0.0002 / +0.0005 / 0 dex, EW +1.4 / -0.3 /
0 percent, depth +1.3 / -0.3 / 0 percent. The 45 and 60 R_p entries bracket
the 30 R_p one from the same side and differ from each other by less than
their difference from it, which is the sense in which the pair is converged
and the 30 R_p value is the outlier.

**T3. Where the 30 R_p boundary reaches into the domain** (raw state, 30 R_p
against 60 R_p at the same radius)

| r [R_p] | T [K] 30 / 60 | rho 30 / 60 [m_H/cm^3] | v 30 / 60 [cm/s] | dT | drho | dv |
|---|---|---|---|---|---|---|
| 1.05 | 1578.0 / 1577.6 | 4.447e11 / 4.440e11 | 5.708 / 5.721 | +0.03% | +0.17% | -0.22% |
| 2.00 | 4105.8 / 4100.1 | 2.216e8 / 2.214e8 | 3156 / 3160 | +0.14% | +0.09% | -0.13% |
| 4.00 | 2056.2 / 2046.7 | 1.111e7 / 1.106e7 | 1.574e4 / 1.582e4 | +0.46% | +0.47% | -0.49% |
| 8.00 | 1042.9 / 1025.5 | 9.110e5 / 8.941e5 | 4.800e4 / 4.891e4 | +1.70% | +1.89% | -1.87% |
| 12.00 | 727.4 / 701.2 | 2.553e5 / 2.452e5 | 7.610e4 / 7.929e4 | +3.74% | +4.13% | -4.02% |
| 16.00 | 579.8 / 543.0 | 1.123e5 / 1.045e5 | 9.728e4 / 1.047e5 | +6.76% | +7.52% | -7.05% |
| 20.00 | 499.9 / 449.8 | 6.250e4 / 5.555e4 | 1.119e5 / 1.259e5 | +11.1% | +12.5% | -11.1% |
| 24.00 | 455.5 / 388.7 | 4.032e4 / 3.375e4 | 1.205e5 / 1.440e5 | +17.2% | +19.5% | -16.3% |
| 28.00 | 432.6 / 345.5 | 2.882e4 / 2.240e4 | 1.238e5 / 1.594e5 | +25.2% | +28.7% | -22.3% |

The disturbance decays inward roughly as a power of the distance from the
boundary and is below 0.5 percent inside 4 R_p, which is why the He 2^3S
column and the 10830 line -- built where the metastable population lives,
1.05 to 4 R_p -- move by under 2 percent while the outer profiles move by
tens of percent. The He 2^3S density itself agrees to 0.3 percent at 1.05 to
4 R_p, 1.4 percent at 8 R_p and 8 percent at 16 R_p.

**Figure** `docs/figures/lhs1140b_L8_outer_boundary.pdf`: Mach number, the
normalized numerator `N/(2c^2/r)`, and temperature of the three solved
states, with the exobase marked. The 45 and 60 R_p curves lie on top of each
other; the 30 R_p curve leaves them above about 10 R_p and turns over at its
own boundary.

---

## 5. Why the extended solutions stop depending on the boundary

Above the sonic point every wave speed at the outer face is positive, HLLC
returns `Phys_flux(WL)` and never forms an expression in the ghost
(`Num_Fluxes.f90` clamps `SL = min(0, ...)`), so the ghost reaches cell N
only through that cell's own reconstruction stencil. The free-outflow rule is
then imposing nothing the flow did not already carry, which is the regime the
routine's own comment says it is built for. The measurement here is that the
LHS 1140 b wind reaches that regime at 40.06 R_p and that 45 R_p is far
enough past it: 45 and 60 R_p agree to 0.1 percent everywhere, while 30 R_p,
which stops 11 R_p short of the crossing, does not.

---

## 6. What must accompany a campaign number

The campaign of `LHS1140b/MODELS.md` runs at `Outer radius 30`. On the
evidence above:

1. **`log10 Mdot` may be quoted as it stands.** Over a factor two in outer
   radius it moves by 5e-4 dex, a tenth of the last digit the campaign
   prints. The sentence to attach is that
   the domain ends below the wind's critical point, so the value is the flux
   the interior sets and not a flux the boundary selected -- established by
   this test and not by the 30 R_p run itself.
2. **The He I 10830 metrics carry a 1.7 percent boundary term** (EW 1.1702
   at 30 R_p against 1.1504 and 1.1537 at 45 and 60), with the FWHM
   unmoved to four digits. Quote them with that as a systematic, not as a
   precision.
3. **No profile quantity above about 10 R_p may be quoted from a 30 R_p
   run**: the density, temperature and velocity there are set by the
   truncation (T3). Anything built by integrating the outer wind -- a Lyman
   alpha or H alpha absorption reaching to tens of R_p, a Roche-lobe
   argument, a comparison with an observed high-velocity tail -- needs a
   domain past 45 R_p.
4. **Nothing above 7.65 R_p is a validated continuum result**, at any outer
   radius. `Kn_bulk` passes 0.1 there, the exobase is at 29.0 R_p and the
   continuum critical point at 40.06 R_p is above it, so the model fixes its
   topology in collisionless gas. `collisional_validity.py` returns
   UNVALIDATED for all three runs and that verdict is unchanged by this item.
   The mass-flux insensitivity of point 1 is a property of the discrete
   problem and is not evidence that the continuum answer is right.

---

## 7. Noticed outside the scope of this item, not acted on

- **The advection-corrected temperature and the solved temperature diverge
  in the outer wind.** At 29.03 R_p the certified state carries `T = 429 K`
  and `Hydro_ioniz_adv.txt` carries `1513 K`, a factor 3.5; the two agree to
  2 percent at 2 R_p and separate outward. The Mach number of the outer face
  is 0.521 read from the solved state and 0.351 read from the `_adv` profile,
  and `collisional_validity.py`, which reads `_adv`, therefore reports a
  critical point at 32.5 R_p for the 45 R_p run where the solved state has it
  at 40.05. Which of the two is the wind's temperature above the exobase is
  worth deciding before either is quoted.
- **The `_adv` profiles of the 45 and 60 R_p runs disagree with each other
  above 30 R_p** although the solved states agree to 0.1 percent there: at
  35 R_p the advection-corrected temperature is 447 K in one and 1269 K in
  the other. Below 29 R_p they agree to 0.2 percent. The instability is in
  the post-processing of the collisionless outer region, not in the solve.

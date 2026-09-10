# A handoff column that comes back 2.4 % apart in pressure: the deep boundary temperature, and what actually moves

Written 2026-08-29. Primary measurement record:
`LHS1140b/exhale/deep_temperature_sensitivity/` (README, six scripts, three
tables, thirteen adapter runs under `runs/`). Every number quoted here was
measured for this document on the corrected Photochem 0.9.0 build
(`README_photochem.md`); nothing is carried over from an earlier run except
where it says so.

**Repaired 2026-08-30, section 8.** The diagnosis of sections 1-7 stands
unchanged and is left as it was written; what section 8 said was not repaired
now is, in the adapter, and the two discrete outcomes are gone. Sections 6 and
9 carry the numbers the repair replaces, and say so where they do.

---

## The judgement first

**Numerical, and not an amplification.** The 2.4 % is the extent of Photochem's
altitude grid at the instant the run is declared steady, and nothing else. The
grid is uniform in altitude between a pinned bottom and `top_atmos`
(`vertical_grid`, `photochem/src/photochem_eqns.f90`); `top_atmos` is re-pinned
to the requested top pressure only every `freq_update_TOA = 1000` internal
steps, and the exit test accepts any top pressure within a **factor 3** of it
(`TOA_in_range`, `photochem/extensions/gasgiants.py`). Two runs of the same
problem exit on different sides of that cycle, their grids differ by 2.5 % in
top pressure, and the converged composition differs with the grid.

A change of 2.1e-6 K in the deep boundary temperature does not produce that
difference: it only decides which side of the cycle the integration lands on.
The two outcomes are discrete -- every quantity of the written column takes one
of two values, labelled by the model top pressure -- and inside each of them
the column reproduces to 1e-6 or better across initial guesses spanning
300-450 K and across a factor 1e5 in stopping time.

**It does not limit the closure chain.** The matching pressure does not move at
all (it is an exact inserted node, `1.000000e-06` bar in every profile), and of
the base state EXHALE actually reads at that node the spread is 2.7e-09 in T,
4.8e-09 in r, 5.7e-08 in `q_H2` and **6.5e-08 in He/H**, the quantity the
closure iterates on. Against `element_flux_closure.py`'s `tol = 0.05` that is
six orders of margin. Section 6.

**Corrected in passing.** `insert_level` in `src/utils/lower_profile_schema.py`
interpolated `n_tot` and `rho` alongside `p` and `T`, which broke `n = p/(k_B
T)` at exactly one level of the table -- the matching level, the one a reader
takes its base state from -- by 7.7e-04. Section 7.

---

## 1. The observation, and what it was about

`LHS1140b/exhale/crossings_pc090/results.txt` section 7 recorded, on the
LHS 1140 b He/H = 9.0 arm:

> At He/H = 9.0 the two builds put the deep boundary 2.1e-6 K apart
> (424.9534068835 vs 424.9534089414 K, tropopause 1.2e-7 K apart), and out of
> that the k00 column comes back differing by up to 2.4 percent in pressure and
> up to 15 percent in one trace species, with the elemental O/H the profile
> hands over differing by 0.5 percent. [...] The size looks like the
> photochemical steady-state solver's own tolerance rather than an
> amplification of 2e-6 K, but which it is was not established here.

Three readings have to be fixed before anything can be measured.

- **"the deep boundary"** is `T_deep`, the surface temperature the climate
  solve returns (`radiative_convective_column.solve_radiative_convective_column`
  -> `AdiabatClimate.surface_temperature_bg_gas`), not the deep boundary
  *pressure*, which is the stated input `--climate-p-deep 20.0` bar.
- **"2.4 percent in pressure"** is a comparison at fixed **level index**, row i
  against row i. It is *not* the matching pressure. `p_match` is inserted as an
  exact node by `lower_profile_schema.insert_level` and reads
  `1.000000000000e-06` bar in all thirteen columns written for this document.
- **"the k00 column"** is the handoff profile of closure iteration 0, i.e. one
  `photochem_to_lower_profile.py --climate` run.

The two builds compared there are both photochem 0.9.0 -- one without and one
with the `clima` `epsfcn` repair of `docs/Update_EXHALE_stage1.md` section 91 -- and
the 2.1e-6 K is what that repair moved a composition that already solved. That
comparison turns out not to be needed: the whole spread is reproduced inside
one build, in section 3, on the environment that carries the repair.

The Fortran quoted below is `photochem/src/` in this repository, the tree that
reproduces the installed library (section 91 of `docs/Update_EXHALE_stage1.md`); the
Python is the installed copy, not the source tree.

## 2. Where the deep boundary temperature becomes a pressure

The chain has three pressure grids and only the last one moves.

1. **The climate grid.** `AdiabatClimate` is given `P_surf = 20` bar and
   `P_top = 1e-2` dyn/cm^2 and lays 60 layers between them; the spacing is
   fixed by those two inputs, so `T_deep` does not move it. Measured: the
   deepest climate level is `1.673098e+01` bar in every run here, to all
   printed digits.
2. **The photochemical grid.** `initialize_to_climate_equilibrium_PT` puts the
   model bottom at the deepest quench level divided by
   `BOA_pressure_factor` (1 here, so the column bottom itself) and then builds
   its own grid -- **in altitude**. `vertical_grid` sets
   `dz = (top_atmos - bottom)/nz` and places the nodes uniformly in z, so a
   single number, `top_atmos`, fixes the pressure of every level. This is where
   a temperature becomes a pressure: `p(z_i)` follows from hydrostatic
   integration through the solution's own T and mean molecular weight.
3. **The written table.** The adapter writes `p = sol['pressure']` as it comes
   back from the solver, then inserts one node at `p_match`.

`top_atmos` is a solved quantity, and the adapter already says so in a comment
at `photochem_to_lower_profile.py`: "the grid the solver settles on need not be
the grid it started from".

## 3. Reproducing it inside one build

The climate solve is the same physical problem whatever `--climate-t-deep-guess`
says: the guess is where MINPACK starts looking, not an input to the physics.
Ten guesses over 300-450 K
(`LHS1140b/exhale/deep_temperature_sensitivity/climate_guess_ladder_heh9.txt`):

| quantity | value | spread over the ladder |
|---|---|---|
| `T_deep` | 424.95340719 -- 424.95340876 K | **2.121e-06 K** (4.99e-09 relative) |
| `P_trop` | 2.40746907 -- 2.40746913 bar | 2.55e-08 relative |
| `T_trop` | 185.443185350 -- 185.443185858 K | 2.74e-09 relative |
| `f_H2O` above the tropopause | 6.0837806e-08 -- 6.0837810e-08 | 6.61e-08 relative |

The 2.1e-6 K of the observation is therefore not a difference between builds at
all: it is **the climate root solve's own reproducibility**, and it is reached
by changing the initial guess by 1e-6 K. Re-run in a fresh process the ladder
is bit-identical.

Seven full adapter runs at seven of those guesses
(`runs/g4*`, `spread_table.txt`):

| run | guess [K] | `T_deep` [K] | steady at t [s] | exit branch | p_top [bar] | O/H at the match |
|---|---:|---:|---:|---|---:|---:|
| `g400a` | 400.000000 | 424.95340688 | 1.288e+17 | `equilibrium_time` | 1.09345112e-08 | 8.846077187e-07 |
| `g400b` | 400.000000 | 424.95340688 | 1.288e+17 | `equilibrium_time` | 1.09345112e-08 | 8.846077187e-07 |
| `g450` | 450.000000 | 424.95340718 | 3.765e+11 | `conv_longdy` | 1.09347540e-08 | 8.846114524e-07 |
| `g400_1em4` | 400.000100 | 424.95340736 | 1.144e+11 | `conv_longdy` | 1.12074787e-08 | 8.892242582e-07 |
| `g401` | 401.000000 | 424.95340782 | 1.803e+17 | `equilibrium_time` | 1.09355075e-08 | 8.846229961e-07 |
| `g410` | 410.000000 | 424.95340876 | 1.142e+11 | `conv_longdy` | 1.12074709e-08 | 8.889333331e-07 |
| `g400_1em6` | 400.000001 | 424.95340900 | 1.169e+11 | `conv_longdy` | 1.12074781e-08 | 8.892181366e-07 |

Sorted by `T_deep`, as the table is, the model top alternates
1.093, 1.093, 1.093, 1.121, 1.094, 1.121, 1.121: **the two outcomes interleave
in the deep boundary temperature and are not separated by a threshold in it.**

`g400a` and `g400b` are the same command run twice and write **byte-identical**
files, so the code is deterministic and nothing below is run-to-run scatter.

The spread over the seven, three ways (`spread_table.txt`):

| quantity | at fixed level index | at fixed pressure | at the matching level |
|---|---:|---:|---:|
| p | 2.470e-02 | -- (the axis) | 0 (exact node) |
| T | 1.099e-03 | 5.061e-04 | 2.741e-09 |
| `n_tot`, `rho` | 2.470e-02 | 9.676e-04 | 4.724e-04 |
| r | 4.996e-05 | 3.103e-07 | 4.849e-09 |
| `q_H2` | 1.961e-04 | 1.762e-05 | 5.694e-08 |
| He/H | 2.081e-04 | 1.956e-05 | **6.484e-08** |
| C/H | 1.130e-03 | 1.078e-04 | 3.072e-07 |
| N/H | 1.187e-03 | 1.133e-04 | 3.206e-07 |
| O/H | 2.140e-02 | 1.345e-02 | **5.207e-03** |
| `q_OH` | 5.837e-01 | 5.765e-01 | 5.747e-03 |
| `q_CO2` | 9.789e-03 | 1.163e-02 | 1.115e-02 |

The first column reproduces every number the observation recorded: 2.4 % in
pressure (2.470e-02), a trace species in the tens of percent (`q_OH`), and
0.5 % in the elemental O/H handed over (5.207e-03 at the match).

## 4. The cause, in the code

**Everything is two-valued, and the label is the model top pressure.** Over all
eleven production and tightened runs the profile top is either
`1.0934-1.0936e-08` bar or `1.12075e-08` bar and never between; the O/H handed
over at the match is either `8.8461e-07` or `8.8892e-07`; `n_tot` there is
either `3.90877e+13` or `3.91059e+13`. Within a group the columns agree to
1e-6 or better; between groups they differ by the numbers of section 3.

**The two groups are two grids.** Removing the inserted node so the level
indices line up, the pressure difference between `g400a` and `g400_1em6` is
zero at the bottom -- both are at `1.6731e+01` bar -- and grows monotonically
and linearly in level index to 2.496e-02 at the top:

| level | 0 | 20 | 40 | 60 | 80 | 100 |
|---|---:|---:|---:|---:|---:|---:|
| p [bar] | 1.67e+01 | 4.96e-01 | 5.38e-03 | 6.31e-05 | 8.00e-07 | 1.09e-08 |
| dp/p | 0 | +5.06e-03 | +1.03e-02 | +1.53e-02 | +2.02e-02 | +2.50e-02 |

That is the signature of a uniform-in-altitude grid whose spacing differs by a
constant relative amount, which is exactly what `dz = (top_atmos - bottom)/nz`
gives for two different `top_atmos`.

**Nothing pins `top_atmos` at the exit.** In `robust_step`
(`photochem/extensions/gasgiants.py`) the accept test is

```python
TOA_in_range = gdat.TOA_pressure_avg/3 < TOA_pressure < gdat.TOA_pressure_avg*3
if chemistry_converged and TOA_in_range:
    reached_steady_state = True
```

and the re-pinning happens on a step counter,

```python
if not (self.wrk.nsteps % gdat.freq_update_TOA) or (chemistry_converged and not TOA_in_range):
    self.update_vertical_grid(TOA_pressure=gdat.TOA_pressure_avg)
```

with `freq_update_TOA = 1000`. A run is therefore accepted with whatever
`top_atmos` the last re-pinning left, anywhere in a factor-3 window. The
observed 2.5 % is two phases of that cycle, and the 2.1e-6 K is enough to
change which one the integration reaches.

## 5. What was excluded, and by what measurement

- **An amplification of 2.1e-6 K.** Excluded, because the outcome is not a
  function of the temperature. Sorted by `T_deep` the seven runs interleave:
  `g401` at 424.95340782 K carries the *lower* model top while `g400_1em4` at
  424.95340736 K carries the higher one. A pair 3.0e-07 K apart in the same
  group (`g400a`, `g450`) agrees to 4.2e-06 in O/H at the match; a pair
  4.6e-07 K apart across the groups (`g400_1em4`, `g401`) differs by 5.2e-03.
  The separation in temperature does not predict the separation in the column;
  the grid group does.
- **The photochemical convergence tolerance.** Excluded. Four runs repeated
  with `conv_longdy` 1e-2 -> 1e-4, `conv_longdydt` 1e-6 -> 1e-8 and
  `equilibrium_time` 1e17 -> 1e22 s (`runs/t1em4_*`,
  `spread_table_tightened.txt`) reproduce the same two groups and the same
  spreads -- 2.466e-02 in pressure at fixed index, 5.199e-03 in O/H at the
  match -- and each tightened run reproduces its loose-tolerance twin to 4e-07
  in O/H. A hundredfold tighter stopping test changes nothing.
- **A second chemical solution, or a branch change.** Excluded. Within a grid
  group the composition reproduces to 1e-06 across initial guesses spanning
  300-450 K and across a factor 1e5 in stopping time. The two groups differ
  smoothly and in the same direction in every species; there is no
  reordering of carriers.
- **A cold-trap cell switching on or off in the climate solve.** Excluded. The
  climate grid does not move (fixed `P_surf`, `P_top`, 60 layers), and over the
  guess ladder the tropopause reproduces to 2.6e-08 relative and the water
  above it to 6.6e-08. The oxygen that does move is the *photochemical* water
  near the tropopause, and it moves with the grid, not with the trap.
- **Round-off amplification, or an ill-conditioned solve.** Excluded. The code
  is bit-deterministic (`g400a` / `g400b` byte-identical; the climate ladder
  bit-identical across processes), and the outcome is two discrete states, not
  a continuum of noise.
- **The matching pressure moving.** Refuted outright: it is an inserted node
  and reads `1.000000000000e-06` bar in every profile. The 2.4 % is at fixed
  level index.

## 6. What it means for the closure chain

`element_flux_closure.py` converges on `eps = |F - Phi|/|F|` for H and He
against `--tol 0.05`, where `F` is the elemental flux the wind carries across a
radial window and `Phi` the trial flux the profile was solved with. The chain
inherits the profile in exactly three ways, and each is measured above:

1. **The base state at the matching level.** `input_read.f90` takes `T`, `r`,
   `q_H2`, He/H and the elemental ratios from the node at `p_match`; it takes
   `p_base` from the header. Spreads: T 2.7e-09, r 4.8e-09, `q_H2` 5.7e-08,
   He/H 6.5e-08, C/H 3.1e-07, N/H 3.2e-07, O/H 5.2e-03. `n_tot` and `rho` are
   carried by the file and explicitly **not** imposed -- the base density is
   `n0` of `input.inp`.
2. **`K_zz(p)` on the grid.** Constant `1e9` in this configuration, identical
   in every run.
3. **`lap_r_top_RJ`**, read as `value_at_pressure(r, p_top_bar)` and used in
   `binary_element_diffusion.f90` as the **upper edge of the elemental flux
   window**. This is the one place the 2.4 % reaches EXHALE, and it moves the
   edge by 8.2e-06 R_J = **5.2e-05 R_p**.

So the handoff's own reproducibility contributes at most 6.5e-08 to the He
residual and nothing to the H residual, against `tol = 0.05` and against the
0.084-0.085 flux-window spread the two LHS 1140 b arms actually carry
(recorded in `crossings_pc090/results.txt` section 7, not re-measured here).
The oxygen reservoir, the one quantity that does move by 0.5 %, does not enter
the H or He residual at all.

**The 2.4 % therefore does not limit the reproducibility of the closure chain.**
It limits the reproducibility of the *stored profile as a table over level
index*, which is not how any consumer reads it. The end-to-end confirmation on
record -- the same two builds giving a steady H flux that agrees to 9.5e-08
relative and a crossing that moves by 2.3e-06 in reservoir He/H -- is quoted
from `crossings_pc090/results.txt` and was not re-measured for this document.

One caveat is worth stating rather than hiding. The closure re-runs the adapter
in every iteration directory (`solve_photochemical_lower_profile`), so each
iterate lands in one grid group or the other independently. The oxygen the
profile hands over therefore carries a 0.5 % step between iterations of the
same ladder. Nothing in the He I 10830 chain reads it, but a result that did
read the oxygen reservoir would inherit that step. (Measured on one rung and
removed by the repair: 1.55e-03 over three iterates before, 8.04e-06 after.
Section 8.3.)

## 7. A defect found on the way: an ideal gas that is not one, at one level

`insert_level` (`src/utils/lower_profile_schema.py`) interpolated **every**
column linearly in log p. Two of them are not free: at every node the solution
was computed on, `n_tot = p/(k_B T)` and `rho = n_tot mu m_u` hold exactly, and
interpolating the densities alongside p and T breaks the identity.

Measured on the columns written here: the identity holds to double-precision
round-off at 101 of 102 levels and is violated at **exactly one -- the matching
level**, by 7.7e-04 (`g400a`) and 1.24e-03 (`g400_1em6`). That is the one level
a reader takes its base state from, and it is larger than the run-to-run spread
of everything else at that node.

`insert_level` now recomputes `n_tot` from the inserted pressure and
temperature, and `rho` from the interpolated mean molecular weight `rho/n_tot`,
which is the intensive quantity to interpolate. Verified by rerunning `g400a`
with the corrected schema (`runs/fix_g400`):

- exactly **two cells** of the file change, `n_tot` and `rho` at the matching
  level; every other cell is bit-identical and the header, including
  `solution_id`, is unchanged;
- `n_tot/(p/k_B T) - 1` over the whole table goes from 7.706e-04 to 2.2e-16;
- the mean molecular weight at the node, 3.898746334 amu, sits between its
  neighbours 3.898749187 and 3.898746241 as it must;
- the run-to-run spread of `n_tot` at the matching level falls from 4.646e-04
  to 4.682e-10, i.e. to the spread of the temperature it is now computed from
  (`runs/fix_g400` against `runs/fix_g400_1em6`).

`make check` does not exercise this path -- the `lower_profile` regression case
reads a stored `lower_atmosphere_profile.dat` and never runs the adapter -- so
no golden is affected and no EXHALE source or binary changed.

## 8. The repair (2026-08-30), and what it removes

Section 8 as written on 2026-08-29 recorded two possible repairs and made
neither. The second of them was made on 2026-08-30, in the adapter, and this
section is what it is and what it measures. The first -- re-pinning and
re-testing convergence inside Photochem, before `reached_steady_state` is set
-- is still the tidier place for it and still needs a rebuild of `photochem/`;
it was not done, and nothing in `photochem/` was touched.

### 8.1 What the adapter now does

`steady_state_at_stated_model_top` in
`src/utils/photochem_to_lower_profile.py` runs after
`photochemical_steady_state` has reached a steady state. It calls
`pc.update_vertical_grid(TOA_pressure=args.toa)`, which re-lays the altitude
grid so that `top_atmos` is the altitude at which the current solution's own
hydrostatic pressure equals the stated top, re-initializes the robust stepper
and runs the chemistry to steady state again.

The two steps are **iterated**, not done once. Re-pinning perturbs the
solution, the perturbed solution moves `top_atmos` again through the
hydrostatic integration, and what the column is written on is their common
fixed point: `top_atmos` is the altitude at which the steady state computed on
the grid it defines reaches the stated pressure. That is a property of the
problem, so the written column no longer records which phase of Photochem's
`freq_update_TOA` cycle the integration happened to exit on.

**The acceptance test is that two successive passes agree on the model top**,
to 1e-05 relative, and not that the model top equals the stated pressure. The
reported top is the top CELL CENTRE while `top_atmos` is the domain edge, and
the half cell between them is about 9 % in pressure at this resolution --
which is also why an unpinned run's top sits 9-12 % above the `--toa` it was
given, in both grid groups, and why the file's `p_top_bar` header has always
read about 1.09e-08 bar for `--toa 1.0e-2`. 1e-05 is 2500 times below the
2.5e-02 that separates the two grids, so what is left of the grid dependence
is about 2e-06 in the elemental O/H at the matching level.

**Failure is a refusal.** If the chemistry does not re-reach steady state on
the pinned grid, or the model top will not settle within six passes, or
`update_vertical_grid` raises, the adapter refuses through `sch.refuse` and
writes no file -- the same class as every other unmet condition it checks. The
alternative, falling back to the unpinned grid, would put back exactly the
dependence this removes, and put it back invisibly: nothing in the written file
records which grid the chemistry stopped on. A `PhotoException` out of
`update_vertical_grid` is reported with the library's own message rather than
re-raised bare, because the adapter can say what it was trying to do and has
nothing to write without it.

### 8.2 The two discrete outcomes are gone

The guess ladder of section 3 re-run through the pinned adapter, seven runs at
the same seven initial guesses (record
`LHS1140b/exhale/deep_temperature_sensitivity/runs/pin_g4*`,
`spread_table_pinned.txt`, written by the unmodified `run_guess.sh` and
`spread_table.py`):

| quantity | unpinned (section 3) | pinned | factor |
|---|---:|---:|---:|
| p, at fixed level index | 2.470e-02 | 2.149e-08 | 1.1e+06 |
| `n_tot`, `rho`, at fixed index | 2.470e-02 | 2.415e-08 | 1.0e+06 |
| `q_OH`, at fixed index | 5.837e-01 | 7.431e-05 | 7.9e+03 |
| T, at the matching level | 2.741e-09 | 2.741e-09 | 1 |
| He/H, at the matching level | 6.484e-08 | 4.791e-13 | 1.4e+05 |
| **O/H, at the matching level** | **5.207e-03** | **2.996e-08** | **1.7e+05** |
| `q_CO2`, at the matching level | 1.115e-02 | 8.384e-05 | 1.3e+02 |

All seven runs now stop on the same grid -- top cell 1.093596e-02 dyn/cm^2 in
every one of them, where the unpinned ladder gave 1.0934-1.0936e-08 bar or
1.12075e-08 bar and nothing between. The two-valuedness of section 4 is not
reduced, it is absent: every column of the profile is single-valued to 1e-07
or better except the trace species, whose worst case is 1.2e-04 (`q_CO2`) where
it was 1.1e-02.

The pinned fixed point is close to, but not identical with, the lower of the
two former groups: against `runs/fix_g400` it moves the top by 1.3e-04 and the
matching-level O/H by 2.5e-05; against the upper group (`runs/fix_g400_1em6`)
by 2.4e-02 and 5.2e-03. The pin cost two passes at every one of the seven
guesses -- the second pass moved the top by 4.2e-08 to 7.6e-06, all inside the
1e-05 test -- each pass a full re-convergence, and about ten per cent of the
adapter's wall time.

### 8.3 What it does to an observable, measured end to end

One elemental-flux closure rung was re-run through the pinned adapter: the
crossing arm of `LHS1140b/exhale/crossings_j96` section 6, reservoir 8.2117,
same binary, same seed, same `closure.json`, same trial fluxes, same tolerance,
same wind protocol -- the adapter is the only difference. Record
`LHS1140b/exhale/adapter_grid_fix/` (the stored ladder was read only).

The column moved by the full size of the effect. The stored rung's last
iterate had stopped on the upper grid (`p_top = 1.121744e-08` bar) and the
pinned one lands at `1.093866e-08`: 2.485e-02 in pressure at fixed level index,
1.543e-03 in the elemental O/H at the matching level. **The line moves by
+1.5e-06 relative:**

| | EW, red pair [%A] |
|---|---:|
| `crossings_j96/rw_cf8p2117`, unpinned adapter | 1.107774992 |
| `adapter_grid_fix/rw_cf8p2117`, pinned adapter | 1.107776648 |
| difference | +1.656e-06 (+1.495e-06 relative) |

and `log10 Mdot` 7.40, `||R||` 3.3e-04, the flux spread 8.63e-02, two
metastable step flags, no metal dropout, `r_drop` 20.799, `outer10` 0.0088 and
`T12` 1679 are the same in both to the precision they are quoted at. This is
the prediction of section 6 measured rather than argued: what the He I 10830
chain reads at the matching level -- T, r, `q_H2`, He/H -- moves by 0, 1.6e-08,
1.6e-07 and 1.6e-07 across the largest grid step there is.

**Against the +0.73 % of `crossings_j96` section 5.** That section split the
+1.96 % the reservoir-9.0 arm moved into +1.22 % from the section-96
secondary-ionization branching and +0.73 % from the adapter, and attributed the
adapter's share to the `insert_level` correction of section 7. The two arms
compared here differ only in the grid, and by the largest grid step there is,
and the line moves by +0.00015 %. The grid is 2.0e-04 of the +0.73 %, so the
+0.73 % is not the grid: the attribution in `crossings_j96` section 5 stands.

**The oxygen step between iterations of one ladder is gone too**, which is the
caveat section 6 ends on. Over the three iterates of this rung, at the matching
level:

| iterate | unpinned `p_top` [bar] | unpinned O/H | pinned `p_top` [bar] | pinned O/H |
|---|---:|---:|---:|---:|
| k00 | 1.093718e-08 | 8.1109417e-07 | 1.093843e-08 | 8.1108833e-07 |
| k01 | 1.093620e-08 | 8.1110377e-07 | 1.093858e-08 | 8.1109268e-07 |
| k02 | 1.121744e-08 | 8.0984548e-07 | 1.093866e-08 | 8.1109486e-07 |
| spread | | **1.552e-03** | | **8.044e-06** |

The He/H at the match marches by 1.2e-05 over the three iterates in both runs
alike -- that is the composition change the closure is making, and it is
untouched. The closure itself is unchanged: the same three iterations, the same
residuals to five digits, the same converged He/H.

### 8.4 What this breaks, and by how much

**Every future handoff column is a different column from the stored ones.**
The `crossings_gm25` ladder (photochem 0.8.4), `crossings_pc090` (0.9.0),
`flux_closure`, `ladder_gm25` and the arms `crossings_j96` re-measured were all
written by an adapter that accepted whatever grid the chemistry stopped on. A
column written now is on the pinned grid, so a stored arm and a re-run of the
same command are no longer comparable at the level of the file. The size of the
break is measured above and is bounded by the two grid groups:

- at most **2.5e-02** in pressure, `n_tot` and `rho` at fixed level index, and
  **1.5e-03 to 5.2e-03** in the elemental O/H handed over at the matching
  level -- everything else at the matching level moves by 1.6e-07 or less;
- **+1.5e-06** in the He I 10830 equivalent width of a closed rung, i.e.
  5.0e-05 %A against a measurement error of 0.030 %A, and, on the closure
  ladder's own slope of EW proportional to (He/H) to the 0.20, about 7.5e-06 in
  the reservoir He/H a crossing is read at -- 6e-05 on a crossing of 8.21.

So no stored result changes at its quoted precision, and none was re-run: the
break is in the file, not in what any published number was read from.

### 8.5 Not done

The Photochem-side repair (section 8, first bullet) -- nothing in `photochem/`
was touched and no rebuild was made. The same measurement at another reservoir
value, another planet, another `K_zz` or another `--toa`; whether the pinned
fixed point is unique, which the seven guesses are consistent with but do not
prove. No EXHALE Fortran source or binary changed and `make check` was not
re-run: the `lower_profile` regression case reads a stored
`lower_atmosphere_profile.dat` and never runs the adapter, so no golden can see
this.

## 9. To carry into `docs/lhs1140b_exhale_vs_pwinds.tex` (not edited here)

- The handoff column's reproducibility floor **before the repair of section 8**
  was 6.5e-08 in the He/H it hands over and 5.2e-03 in the oxygen reservoir,
  both at the 1 microbar matching level; the 2.4 % figure recorded earlier is a
  fixed-level-index comparison and is not what a consumer of the file sees.
- **After the repair** the same floor is 4.8e-13 in He/H and 3.0e-08 in the
  oxygen reservoir, and the column at fixed level index reproduces to 2.1e-08.
  Columns written from 2026-08-30 on are not comparable, at the level of the
  file, with the stored `crossings_gm25`, `crossings_pc090`, `flux_closure` and
  `crossings_j96` arms; the He I 10830 equivalent width of a closed rung moves
  by +1.5e-06 relative between the two, so nothing quoted from those ladders
  changes at its printed precision (section 8.4).
- The deep boundary temperature of the LHS 1140 b He/H = 9.0 climate solution
  is reproducible to 2.1e-06 K, its tropopause to 2.6e-08 relative.

## 10. Files

Measurement record, `LHS1140b/exhale/deep_temperature_sensitivity/`:

| file | what it holds |
|---|---|
| `README.md` | the record's own summary |
| `climate_guess_ladder.py`, `climate_guess_ladder_heh9.txt` | the climate solve alone, from ten initial guesses |
| `run_guess.sh`, `runs/g4*` | seven handoff columns, one per guess, `g400a`/`g400b` the same command twice |
| `run_tightened.py`, `run_tight.sh`, `runs/t1em4_*` | four of them with the stopping test tightened a hundredfold |
| `runs/fix_g400`, `runs/fix_g400_1em6` | two of them repeated after the `insert_level` correction |
| `compare_columns.py` | two columns, at fixed index, at fixed pressure and at the match |
| `spread_table.py`, `spread_table.txt`, `spread_table_tightened.txt` | the spread over a set of runs, the same three ways |
| `runs/pin_g4*`, `spread_table_pinned.txt` | the same seven guesses through the pinned adapter (2026-08-30, section 8.2), written by the same two scripts unchanged |

Second measurement record, `LHS1140b/exhale/adapter_grid_fix/` (2026-08-30,
section 8.3): one elemental-flux closure rung at reservoir 8.2117 re-run
through the pinned adapter, its wind re-solved, and the line measured --
`results.txt`, `reproduce.sh`, `cf_C8p2117/`, `rw_cf8p2117/`, with
`run_closure.sh`, `resolve_wind.sh` and `measure.py` copied unchanged from
`LHS1140b/exhale/crossings_j96/`.

Code: `src/utils/lower_profile_schema.py` (`insert_level`),
`src/utils/photochem_to_lower_profile.py`
(`steady_state_at_stated_model_top`),
`src/utils/radiative_convective_column.py`,
`src/modules/files_IO/lower_atmosphere_profile.f90`,
`src/modules/files_IO/input_read.f90`,
`src/modules/functions/binary_element_diffusion.f90`.
Photochem: `photochem/src/photochem_eqns.f90` (`vertical_grid`),
`photochem/photochem/extensions/gasgiants.py` (`robust_step`).
The observation this answers: `LHS1140b/exhale/crossings_pc090/results.txt`
section 7. The series it came out of: `docs/deep_level_elemental_check.md`,
`docs/Update_EXHALE_stage1.md` sections 89-91. The handoff design:
`docs/phase_e_flux_closure_design.md` sections 2 and 6.

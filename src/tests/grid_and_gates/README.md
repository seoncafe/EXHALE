# Grid and gate tests

Nineteen tests of the discretization and of the run's own record of itself. Six
were built for Phase 0 of `docs/development_plan_20260905_rev3.md` (section
10.4, "Phase 0 gains tests"); `grid_window` and `sed_coverage` came with
Phase 1 batch 2a, items 2a-GUARDS and 2a-SED; `momentum_row`,
`threshold_edges`, `hydrostatic_residual` and `free_outflow_boundary` came
later, with the items named in their entries.
Most of them target a defect the independent source audit of
2026-09-05 recorded in section 10.1 or 10.2, so **most of them are RED at HEAD
`35d9dd5` by design**: they are the statements the Phase 1 corrections have to
turn green, and they exist now so that the corrections can be judged by
something other than a golden diff. Which of them have turned green since is
the summary at the end.

Run them all:

```
src/tests/grid_and_gates/run.sh                 # every test
src/tests/grid_and_gates/run.sh grid_width      # one of them
```

Names: `grid_width`, `grid_window`, `photon_quadrature`,
`threshold_edges`, `hydrostatic_residual`, `free_outflow_boundary`,
`flux_spread`, `output_state`, `base_level`, `restart_grid`,
`restart_round_trip`, `restart_intent`, `restart_option_change`,
`direct_steady_setup`, `sed_coverage`, `coupled_carrier_h2`,
`carrier_transport_inert`, `momentum_row`, `fpe_traps`.

`run.sh` builds into `build/tests/grid_and_gates/`, prints one
`PASS|FAIL <name> measured=<v> reference=<r> tol=<t>` line per assertion, and
exits nonzero if any assertion fails. Lines beginning with two spaces, or
carrying `DIAGNOSTIC`, are context and carry no verdict.

Three variables point the suite at another build, which is what lets it run
while someone else is rebuilding `build/`:

```bash
make OBJDIR=build_x EXE=EXHALE_x.x
EXHALE_OBJDIR=$PWD/build_x EXHALE_EXE=$PWD/EXHALE_x.x EXHALE_TEST_OUT=/tmp/gg \
   src/tests/grid_and_gates/run.sh
```

`EXHALE_OBJDIR` is the object directory the Fortran drivers link, `EXHALE_EXE`
the binary the shell tests run, and `EXHALE_TEST_OUT` where the test
executables and the Fortran drivers' scratch directories are written. With
either of the first two set the `make -q` staleness check is skipped, because
the objects then need not match the default `build/`. With none of them set
the suite behaves exactly as described above, staleness check included. The
shell tests place their copies of the regression cases under
`build/tests/grid_and_gates/` whatever `EXHALE_TEST_OUT` says, since each of
them resolves that path itself.

Nothing here re-implements the quantity it tests. The Fortran drivers
link the **production objects** in `build/` (`run.sh` refuses to run if
`make -q` says they are behind the sources), so the grid they read is the one
`define_grid` builds and the photon grid is the one `set_energy_vectors`
builds, and the right-hand side is the one `RK_rhs` and `Source` assemble.
The shell tests run `EXHALE.x` itself on **copies** of
regression cases placed under `build/tests/grid_and_gates/`; no regression
case directory and no golden is written to. The remaining test reads the
baseline output files with Python.

## The tests

### 1. `grid_width_identity.f90` -> `grid_width`

| | |
|---|---|
| Origin | section 10.1 item 1 (`dr_j(j)` is the width of cell j+1) |
| Quantity | `dr_j(j)` and `r_edg(j) - r_edg(j-1)`, and `sum_{j=1..N} dr_j(j)` against `r_edg(N) - r_edg(0)` |
| Reference | the identity itself. `r_edg(j) = r_{j+1/2}` is the stated convention of `define_grid.f90` and the one `RK_rhs` uses (`rp = r_edg(j)`, `rm = r_edg(j-1)` for cell j), so the width of cell j is `r_edg(j) - r_edg(j-1)` exactly |
| Tolerance | 1e-12 relative (an exact identity of the discretization) |
| Configurations | Mixed at `r_max` 10 and 2, Uniform at 10, Stretched at 10, all with N = 500, `N_low_cells` = 50, `dr_base` = 2.0e-4 |
| Expected at HEAD | **RED** on Mixed and Stretched, GREEN on Uniform |

`define_grid.f90` 152 and 185 write
`dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)`, which is
`r_edg(j+1) - r_edg(j)`: the width of the **next** cell. Uniform passes
because every cell has the same width there. The ratio the test prints,
`dr_j(j)/(r_edg(j)-r_edg(j-1))`, is the local stretch factor: measured
1.0145082 (Mixed, `r_max` = 10), 1.0100000 (Mixed, `r_max` = 2), 1.0046033
(Stretched, `r_max` = 10).

### 2. `photon_grid_quadrature.f90` -> `photon_quadrature`

| | |
|---|---|
| Origin | section 10.1 item 2 (the bin centered on a threshold charges the cross section below the threshold) |
| Quantity | `P_HI`, `P_HeI`, `P_HeII` and the H I photoheating rate of one atom, at zero optical depth, formed exactly as `util_ion_eq.f90` 770-772, 1015-1017, 1043 and 1055 form them, on the arrays `set_energy_vectors` filled |
| Reference | the same integrand from the same production functions (`J_inc`, `photoion_sigma`, `photoelectron_share`), trapezoid rule on 400001 logarithmic points over the photon band of the run, `[grid floor, e_top]`. A second reference on 200001 points is printed: the fine integral is converged to 3e-5 |
| Tolerance | 1e-3 relative (the tolerance section 10.1 item 2 states for this gate) |
| Configurations | (i) `Include He23S? True`, grid floor 4.80 eV; (ii) no metastable, Na I active, floor 5.139 eV; (iii) neither, floor 13.6 eV. Spectrum from `backup/regression/wasp_full/input.inp`: Power-law, index -1, [13.60, 123.98, 1.24e3] eV, LX 29.46, LEUV 30.42, a = 0.02544 AU, `Rate/2 + Mdot/2` |
| Expected | **RED before item 2c-QUAD** on `P_HI` in (i) and (ii), on `P_HeI` and on `P_HeII` in all three; GREEN on the heating and on `P_HI` in (iii). **After it**: GREEN except `P_HeI`, which is 1.2e-3 low for a reason of its own, below |

Code/exact ratios before 2c-QUAD: `P_HI` 1.095077, 1.089024, 1.000338;
`P_HeI` 1.015747, 1.015731, 1.015756; `P_HeII` 1.029919, 1.029878, 1.029884;
heating 0.999840 in all three. `P_HeI` is not named in section 10.1 item 2 but
failed for the same reason: 24.6 eV was a grid point in every configuration,
so its bin was charged over its sub-threshold half whether or not the
sub-Lyman grid was present.

After 2c-QUAD, with the thresholds on bin edges: `P_HI` 0.999878, 0.999867,
0.999919; `P_HeII` 0.999828, 0.999821, 0.999799; heating 1.000063; `P_HeI`
0.998775, 0.998783, 0.998772. The remaining 1.2e-3 of `P_HeI` is not a
quadrature error: the He I band edge is the global constant
`e_th_HeI` = 24.6 eV while the run's He I cross section, the VFKY96 fit of
`cross_sec.f90` 66, turns on at its own `E_th` = 24.59 eV, so the band
[24.59, 24.60] eV carries a cross section no bin center of the grid sees. That
band is 0.115% of the exact `P_HeI`, which is the whole of the shortfall.
Closing it means writing the He I ionization potential once (24.5874 eV, NIST
ASD) instead of twice.

The comparison interval is the photon band `[grid floor, e_top]`, which the
bin partition covers exactly; `sum(de_v)` against that width is printed as a
DIAGNOSTIC (0 after 2c-QUAD; 4.5% short of it before, the grid having stopped
at the last X-ray POINT, 1184.19 eV, rather than at `e_top` = 1240 eV).

### 3. `mass_flux_spread_functional.py` -> `flux_spread`

| | |
|---|---|
| Origin | section 10.2 item 5, decision 15 (the stop measured the spread of `abs(rho v r^2)` over its minimum while the reported `flux_spread` measured the signed spread over the mean: two functionals of one window) |
| Quantity | the value of `mass_flux_spread` the run used at its last step, printed in the run summary with the window it was taken over, against the same functional recomputed from the `output/Hydro_ioniz.txt` that run wrote; and the functional on a synthetic window carrying `-F` over its inner half and `+3F` over its outer half |
| Reference | assertion 1: the number in the run summary. Assertion 2: 4, the value of `(max - min)/abs(mean)` on that window (the magnitude-first, minimum-normalized form reads 2 on it) |
| Tolerance | 1e-9 relative on assertion 1 (full-precision columns, a scale-free functional), 1e-12 on assertion 2 |
| Configuration | a copy of `backup/regression/mol_base_handoff` capped at 200 steps, in `build/tests/grid_and_gates/flux_spread` |
| Expected before item 2c-DU | **RED**: the run summary carries no `mass_flux_spread:` line, because neither functional was a named quantity |

The test no longer transcribes the production expression. `mass_flux_spread`
is a contained function of the main program and cannot be linked against from
a driver, so the run prints the value it used:

```
     mass_flux_spread: window cells 392..500  value=  6.961360897137848E-01
```

and the test recomputes it from the profile the same run wrote. The 200-step
cap matters: the run must end by marching, since a JFNK finish moves the state
after the last `du`.

Measured with item 2c-DU in place, `mol_base_handoff` at 200 steps: window
cells 392..500 (`r` = 2.003 to 4.725, no sign change), summary 2.826589937536,
recomputed 2.826589937536.

### 4. `output_state_consistency.sh` -> `output_state`

Two checks since 2026-09-11 (item P14): the heat-column consistency below, and
`outer_iteration_ending_hands_back_one_state`, which forces the stationary
outer iteration's stagnation ending on the hot-Uranus carrier reload
(`EXHALE_CARRIER_TRUST=1e-6`, `EXHALE_JFNK_MAXIT=40`, since the joint
progress measure of 2026-09-12) and asserts that the
written state re-evaluates to itself (relative temperature movement below
1e-12; MEASURED 6.5e-14 after P14 against 1.3e-8 before).

| | |
|---|---|
| Origin | section 10.2 item 2 (the final output claims a coupling its heat column was not produced under) |
| Quantity | `heat` of `output/Hydro_ioniz.txt` against `heat_total` of `output/Heating_breakdown.txt`, cell by cell over the physical cells, next to the `# coupling: sec_ion=` flag the same run wrote |
| Reference | 0 relative difference. The breakdown file's own header states the relation ("Channel sum reproduces the heat_total column (and the Hydro_ioniz.txt heat column up to convergence)"): on one state the two columns are one number |
| Tolerance | 1e-6 relative |
| What is run | a copy of `backup/regression/mol_base_handoff` (its `input.inp` and `base.inp`) capped at 200 steps |
| Expected at HEAD | **RED** |

Measured at 200 steps: the header says `sec_ion=F`, and `heat_total/heat` has
median 1.000000, minimum 0.999992, maximum 1.000012 over the 500 physical
cells; the largest relative difference is 1.19e-5.

What remains is a rate lag of one sweep, not a coupling mismatch. The `heat`
column is what the ionization sweep returned: the post-sweep composition
contracted with the pre-sweep rates
(`ionization_equilibrium.f90` 1983-1989). `Heating_breakdown.txt` is written
afterwards by `write_heat_breakdown_eq` (`util_ion_eq.f90` 2142), which
re-solves the radiation field and every rate on the state as written, with an
`excited_H_update` on that state in between (`EXHALE_main.f90` 1933). The
final write leaves `sec_ion_active` as the marching loop left it
(`EXHALE_main.f90` 1826-1848); the two activations are inside the loop
(`:1649`, `:1698`). The difference is a property of one sweep, so the shorter
run is as informative as the 12000-step snapshot.

### 5. `base_level_single_statement.sh` -> `base_level`

| | |
|---|---|
| Origin | section 10.2 item 4 (with a lower-atmosphere profile the base level is stated twice and never compared) |
| Quantity | the base pressure in `EXHALE_setup.out` ("Base level: ... -> p = ... bar", from the density key) against the profile's matching level in the run log ("profile: p_base -> ... bar") |
| Reference | the profile's `p_match_bar`: the pressure at which the wind starts is one number |
| Tolerance | 1% relative |
| What is run | a copy of `backup/regression/lower_profile` (`input.inp`, `base.inp`, `lower_atmosphere_profile.dat`) capped at 2 steps; both lines are written during startup |
| Expected at HEAD | **RED**, by a factor of about 20 |

Measured: 2.0035e-5 bar in the setup report against 1.00e-6 bar in the log, a
ratio of 20.035. `input_read.f90` 1652 skips the base-level consistency check
whenever a profile is in use, and `apply_lower_atmosphere_profile` (2326)
adopts `p_match_bar` as `p_base_bar` while `n0` still comes from the density
key, so nothing in the code compares them.

### 6. `restart_grid_guard.sh` -> `restart_grid`

| | |
|---|---|
| Origin | section 10.2 item 3 (a restart overwrites the cell centers with the file's radii and keeps the run's faces, widths and window indices) |
| Quantity | the `r` column of a run started with `Load IC? True` from state files written on a different grid |
| Reference | either a nonzero exit (the mismatch refused) or the run's own grid, measured by a second cold run of the same input |
| Tolerance | 1e-12 relative |
| What is run | three short runs from `backup/regression/roundtrip`: (A) as shipped, `Grid type: Mixed`, 40 steps, the state to restart from; (B) the same input with `Grid type: Uniform`, no restart, 1 step, whose `r` column is the run's own grid; (C) `Grid type: Uniform` with A's outputs handed back as `output/*_IC.txt`, 1 step |
| Expected at HEAD | **GREEN**, by the refusal branch |

Measured 2026-09-10: C stops in `load_IC` with "holds cell centers of a
different grid from the one this run built" and exit 1, so the row passes as
`restart_grid_mismatch_refused` and the `r`-column comparison is not reached.
The guard compares the centers of cells 1..N against the ones `define_grid`
built, to 1e-10 relative; the row's own 1e-12 tolerance applies only to the
second branch, where a run that accepted the load would have to have kept its
own grid.

### 7. `grid_window_indices.f90` -> `grid_window`

| | |
|---|---|
| Origin | section 10.3 (`minloc(..., dim=1)` used as a declared subscript in `set_IC.f90` 116-118; `j_min` has no lower clamp in `define_grid.f90` 193-197) |
| Quantity | the cell index `cell_nearest_radius(x)` returns, and the first indices `j_min`, `j_flux` of the constant-momentum and flux windows |
| Reference | an independent brute-force scan of `r` in the driver: the cell that minimizes \|r(j) - x\| over 1-Ng..N+Ng, and the first cell with r >= the window radius; and the first physical cell, 1, for the clamp |
| Tolerance | 0 (integer indices) |
| Configurations | Mixed grid, N = 500, `N_low_cells` = 50, `dr_base` = 2.0e-4, `r_max` = 10; windows at `r_esc` = 1.50 / `r_flux` = 1.20 (an ordinary case) and at `r_esc` = `r_flux` = 1.00 (escape radius at the base radius) |
| Expected at HEAD | **RED**: `cell_nearest_radius` does not exist at 35d9dd5 (the driver does not build), and `j_min` = 0 at `r_esc` = 1.00 |

`minloc` counts positions from 1 whatever the declared lower bound is, so on
the grid arrays (declared `1-Ng:N+Ng`) a bare `minloc` result names a cell Ng
below the one meant; `set_IC` read the mid-domain density through it. The
window search returns `j = 2-Ng = 0` when the escape radius is the base
radius, which opens the convergence window in the ghost cells.

### 8. `sed_coverage_stop.sh` -> `sed_coverage`

| | |
|---|---|
| Origin | section 10.5 decision 17 (a loaded SED that ends above the photon-grid floor an active absorber needs) |
| Quantity | the exit status of `EXHALE.x` and the text of its message, for a run whose `Spectrum type: Load` table stops above the grid floor |
| Reference | decision 17: the run stops, there is no key to continue, and the message names the file, the missing band in eV and in A, the absorbers that need it and the remedies |
| Tolerance | 0 (an exit status, the presence of named strings, and the remedy wavelength to 1 A) |
| Configurations | `benchmarks/wasp52`, one step each: (A) as shipped (`Include He23S? True`, both stellar lines) with `Spectrum file:` set to `inputdata/scaled_solar_hd189.sed` (last row 1899.5 A = 6.5272 eV), floor 3.400 eV = 3647.0 A; (B) the same table with `Include He23S? False`, no metal and the stellar lines deleted, floor 13.60 eV, which the table covers; (C) B plus a `metals.inp` carrying Na (5.139 eV = 2412.6 A) and K (4.341 eV = 2856.1 A), floor 4.341 eV; (D) `Include He23S? False` with the stellar lines kept and the eps Eri table cut at 2600 A, so that H(n=2) is the only absorber below the file's lowest photon; (E) the benchmark as shipped, with its own eps Eri table reaching 3700 A |
| Expected before the Phase 1 (2a) change | **RED** on the four assertions of A and C: the branch printed a warning and the run continued |
| Expected before item 2c-N2FLOOR | **RED** on the four H(n=2) assertions: D stopped in `excited_hydrogen` and named no file, and E stopped although its table reaches 3647 A |

Eleven assertions: A, C and D exit nonzero and their messages name the file,
the band and the absorbers that lose it, A's remedy wavelength is the n=2 edge
3647.0 A to 1 A, D names H(n=2) and not the metastable; B and E exit 0 and
write `output/Hydro_ioniz.txt`. Configuration C is the branch no shipped
`metals.inp` reaches, since none of them lists K, whose neutral threshold is
the lowest in `species_table.f90`. E is the WASP-52b configuration, which the
n=2 grid floor is what lets run at all.

### 9. `momentum_row_from_fluxes_only.sh` -> `momentum_row`

| | |
|---|---|
| Origin | section 10.2 item 1 (`u(2,:) = u(2,:) + 1.0e-16` on every step, to keep the next step's `dtu` denominator away from zero) |
| Quantity | assertion 1: the number of statements in `src/EXHALE_main.f90` that add a numeric constant to a whole row of the conserved state. Assertion 2: the number of steps of a 200-step run whose reported `dtu` is NaN or Infinity |
| Reference | 0 and 0 |
| Tolerance | 0 (a count) |
| Configuration | a copy of `backup/regression/mol_base_handoff` capped at 200 steps, in `build/tests/grid_and_gates/momentum_row` |
| Expected before item 2c-KICK | assertion 1 **RED** (one statement, `EXHALE_main.f90` 1357), assertion 2 green |

Assertion 1 is a source assertion, not a run measurement: an increment of
1.0e-16 of the code momentum unit is absorbed into a state that stays
self-consistent, `u`, `rho` and `v` being re-derived from each other every
step, so no single run's output files expose it. What a pair of runs shows is
a changed state, which is a movement measurement and not a criterion. The
explicit energy source `u(3,:) = u(3,:) + dt_loc*(heat - cool)` carries
`dt_loc` and is an operator of the split, so it does not match the pattern.

Assertion 2 guards the replacement: with the increment gone, a conserved row
that is exactly zero at the start of a step has no relative change defined,
and the run must report an unbounded change rather than a NaN that no
comparison can act on.

### 10. `photon_grid_threshold_edges.f90` -> `threshold_edges`

| | |
|---|---|
| Origin | item 2c-METALEDGE: item 2c-QUAD put the H and He thresholds on bin edges, but every OTHER threshold of an active absorber -- the He 2^3S metastable when a lower metal floors the grid, H2, and each active metal ion -- still fell inside a bin, which then charged that cross section over its sub-threshold part |
| Quantity | assertion (a): the relative distance from each active ionization threshold to the nearest bin edge of the partition `set_energy_vectors` built, the edges rebuilt from `photon_grid_floor_eV()` and the exact widths `de_v`. Assertion (b): `P_m` for C I, C II, O I, Mg I, Mg II, Ca I, Na I, Fe I at zero optical depth, formed exactly as `util_ion_eq.f90` 1150-1157 forms it. Plus `sum(de_v)` against the band width, and `Nl` against `Nl_fix + NlTR` |
| Reference | (a) zero. (b) the same integrand from the same production functions (`J_inc`, `metal_photoion_sigma`), trapezoid rule on 400001 logarithmic points over `[threshold, e_top]`; a 200001-point integral is printed, the fine one being converged to 8e-6 or better |
| Tolerance | (a) 1e-12 relative (an edge either is the threshold or is not; only the round-off of rebuilding the edges is admissible). (b) 1e-3 relative for every ion, the four neutral metals of the sub-Lyman band included since that band carries 89 bins (the driver's header carries the measured `num_TR` ladder that fixed the budget) |
| Configurations | `wasp_full` (He 2^3S on, stellar lines set, solar C/N/O/Mg/Ca/Na/Fe) and `mol_metals` (the same metals on the hot Uranus with `Molecular chemistry: True`, so H2 is an absorber too) |
| Expected before item 2c-METALEDGE | **RED**: both `threshold_on_bin_edge` assertions (worst Na I, 2.23e-2 of its own threshold) and six of the eight rates |

Code/exact before, `wasp_full`: C I 1.044961, C II 1.007546, O I 1.000425,
Mg I 1.000143, Mg II 0.991179, Ca I 0.933302, Na I 1.030389, Fe I 0.965337.
After: 0.999996, 0.999994, 1.000024, 0.997527, 1.000506, 0.997368, 0.997523,
0.994182. Mg I passed before by accident, its 7.646 eV threshold falling
6.7e-4 from a bin edge of the old partition.

### 11. `hydrostatic_residual.f90` -> `hydrostatic_residual`

| | |
|---|---|
| Origin | increment B4-5 of `docs/b4_spatial_operator_design_20260906.md` section 2.6 (the hydrostatic residual as an operator identity), D2 of `docs/development_plan_20260905_rev3.md` section 4.5 |
| Quantity | `R_mom(j) = dF(2,j) - S(2,j)`, the momentum one call of the production right-hand side adds to a state that should stay at rest, normalized by the local weight `rho_j Dphi(r_j)`; plus the velocity `dv = R_mom dt_CFL/rho` it would produce, over the sound speed |
| State | the **analytic** isothermal hydrostatic column of `backup/regression/hydrostatic_column` (0.0457 M_J, 0.49 R_J, `T_eq` 1140 K, He/H 0.0793, Spherical, `r_max` 3), `rho(r) = exp[-(phi(r) - phi(1))/nhat]`, `p = nhat rho`, `v = 0`, as cell volume averages by an 8-point Gauss-Legendre rule. The marched case has no discrete hydrostatic state to measure against, which is why the analytic one is the fixture |
| Reference | assertion 1: zero, exactly. With `b0 = 0` and a uniform state the pressure divergence and the geometric source are the same number formed two ways. Assertion 2: 1, the ratio of the interior maximum at N = 2000 to the one at N = 250. Assertion 3: 1, the lowest order in `1/N` a consistent outer boundary closure can have |
| Tolerance | 1e-14 relative on assertion 1 (round-off of that pair); 0 on assertions 2 and 3, which are inequalities |
| Configurations | `Grid cells:` 250, 500, 1000, 2000, each with PLM and ESWENO3 and with the ROE and HLLC fluxes: 12 assertions over 16 measurements |
| Expected at HEAD | **GREEN**, all three. Assertion 3 was RED under the zero-gradient outflow ghost item B4-4a replaced (order 0.006 in `1/N` under PLM, both Riemann solvers) and is green with the hydrostatic continuation (1.80 under PLM, 1.59 under ESWENO3) |

The driver carries three further groups of assertions on the momentum row's
reference scale, on the same ladder and the same four scheme/solver pairs, so
the whole program is 60 assertions.

**The well-balanced key on its own discrete equilibrium (12, item N37).** The state whose
two neighboring equilibrium extrapolations agree at every shared face,
integrated outward from the analytic density of cell 1, is built beside the
analytic column; with the key on, the momentum, mass and energy rows of the
interior cells 3..N-2 must fall to the rounding level there (bound 1e-13 of
the size of the terms each row is built from, the arithmetic bound
`epsilon p r^2/dV` being 1e-12 of the momentum weight at the base of this
grid). MEASURED: the momentum row falls from 7.0e-5 (PLM) and 2.5e-3 (WENO3)
of the local weight with the key off to 4.1e-14 and 3.9e-14 with it at
N = 250, and stays at 4.5e-14 to 5.5e-14 at N = 500, 1000, 2000, while the
base scheme falls at its design order. The analytic column is measured too
and reported: on it the key and the base scheme both carry the truncation error of the
discretization of gravity, which is what separates that error from the
floating-point assembly.

**The key's momentum row scale (16, item P4 of
`docs/PLAN_20260911_partitioned_solver.md`).** The quantity is the row
DIVIDED BY ITS OWN REFERENCE SCALE, `residual_row_scale(2,...)`, the one
expression the stationary solver's acceptance test reads. With the key on the
equilibrium pressure force and the gravitational source cancel in the algebra
before the row is formed, so a scale built from the row's remaining terms
alone is the row itself and reads one however small the imbalance becomes.
The four statements: on the discrete equilibrium the scaled row is at the
rounding level (bound 1e-10, twelve decades below what a degenerate scale
returns; MEASURED 4.5e-14 to 7.1e-14); on the same column with the pressure
multiplied by `1 + eps sin(...)` it is below one (MEASURED 5.3e-6 at
eps = 1e-6) and halves with the perturbation (MEASURED ratio 2.0000 to 1e-8,
the departure being homogeneous of degree one in eps); and with `b0 = 0` the
equilibrium pressure force is zero and the scale is the row's own terms
`max(|ram|, |dp/dr|)` exactly, each formed in the driver from the face data
the evaluation stored (MEASURED 0). The departure of the key's scale from the
BASE SCHEME's at zero gravity is reported and not asserted: the two are the
same expression up to the rounding of the two flux assemblies (MEASURED 1.7e-13
to 3.7e-13 over the four pairs). Expected on the entry text of P4: the first
three groups RED at exactly 1.000 (the degenerate scale), the zero-gravity
limit GREEN. Its fourth statement was RED on the entry text of P11 as well
(MEASURED 0.750 for all four pairs), because it read the scale against
`max(|dF_2|, |S_2|)`, which under the key is the SUM of the two terms and is
smaller than either wherever they cancel; 0.944 of what that statement
reported at zero gravity was the PLM geometric pressure term, which item P11
removed from the scale.

**The scale against the terms of the momentum equation (20, item P11 of
`docs/PLAN_20260911_partitioned_solver.md`).** The scale is the largest
PHYSICAL term of
`d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc`, and the
discretization's pieces `dF_2` and `S_2` are not those terms: under PLM the
pressure sits partly in the momentum flux (`Phys_flux` adds `p`) and partly in
the geometric source `(A+ - A-) p_c/dV`, each carrying an O(2 p/r) part that
cancels against the other. Four states whose terms are known in closed form,
the base scheme and the key, interior cells 3..N-2:

| statement | bound | entry text | after |
|---|---|---|---|
| `zero_gravity_uniform_state_momentum_scale_vanishes` (no gravity, uniform `p`, at rest; the scale over `(A+ - A-) p_j/dV`) | 1e-8 | **1.00000E+00 FAIL** under PLM, 3.3e-308 PASS under WENO3 | 3.3e-308 (the `tiny` floor of `momentum_row_scale`, which is where a fully zero row is divided) |
| `momentum_scale_is_not_below_the_weight_on_the_discrete_equilibrium` (`(w - s)/w`, with `w` the weight `source` forms from the two face densities) | 1e-12 | **2.79E-01 FAIL** under PLM, 0 PASS under WENO3 | 0, the weight being one of the terms the max runs over |
| `well_balanced_momentum_scale_is_the_equilibrium_pressure_force` (the same state under the key, `(|s - epf|- |R|)/epf`) | 1e-12 | 0 PASS | 0 |
| `momentum_scaled_row_tracks_the_imbalance` (the pressure perturbed by `eps` and `eps/2`, the ratio of the scaled imbalance the perturbation ADDS) | \|ratio - 2\| < 0.1 | 2.0000 PASS | 2.0000 |
| `supersonic_uniform_flow_momentum_scale_is_the_ram_divergence` (no gravity, uniform `p`, `v = 3 c_s`; against `(A+ - A-) rho v v/dV`) | 1e-8 | **6.66667E-02 FAIL** under PLM, which is `1/(gamma M M)`, the pressure the momentum flux still carries; 4.2e-13 PASS under WENO3 | 1.9e-14 (PLM), 4.2e-13 (WENO3) |

Ten of the sixty assertions were RED on the entry text of P11 and all sixty
are GREEN after. Two quantities are reported and not asserted, because on the
base scheme they are neither round-off nor a term of the equation: the
departure `|s - w|/w` on the discrete equilibrium (MEASURED 8.5e-4 at N = 250
falling to 2.3e-6 at N = 2000 under PLM with ROE, 2.2e-3 to 4.5e-5 under
WENO3 with HLLC), which carries the base scheme's own truncation error on a
state that is the equilibrium of the WELL-BALANCED reconstruction and, with HLLC, the
dissipative part of a momentum flux whose face pressure is one side's own
value rather than an average; and the zero-gravity scale itself.

The undifferenced form of the halving statement is why the perturbation
statement subtracts the unperturbed row: the base scheme's row on the
well-balanced discrete equilibrium is its own truncation error, which is larger than a
1e-6 perturbation's imbalance and does not halve with it (MEASURED on the
entry text: the undifferenced ratio is 1.0012 at N = 250 and 1.586 at
N = 2000, the truncation falling with the grid).

One grid per process: `define_grid`'s Newton convergence flag is
`real*8 :: tol = 1.0`, which is SAVEd, so a second grid built in the same
process keeps the initial stretch guess 1.01 (section 8.4 of the design note).
`run.sh` therefore invokes the driver once per `N`, each appending its records
to `hydrostatic_ladder.dat` in the working directory, and once more with no
argument to read them back, print the ladder and give the verdicts. Invoked as
`hydrostatic_residual.x <N> [dr_base N_low]` it measures a base-refinement
ladder instead, which is what section 8.2 of the design note reports.

MEASURED 2026-09-06 on the prescribed ladder, under the ghost B4-4a replaced:
in the deep interior `R` falls as `dr^1.96` (PLM) and `dr^1.98` (WENO3) in the
local width, so the pair is at PLM's design order and one below WENO3's; the
largest residual was at the outer free-outflow cell, 0.750 of its own weight
under PLM at every resolution (order 0.006 in `1/N`) and 6.4e-2 falling at
order 1.43 under WENO3; the base cell is 1.0e-4 (PLM) / 1.6e-4 (WENO3) and
first order in `dr_j(1)`; the largest interior residual sits at the
uniform/stretched grid junction `r = 1.0100`. Full tables and the verdicts:
design note section 8.

MEASURED 2026-09-06 with the isothermal hydrostatic outflow ghost of B4-4a:
the interior and the base cell are unchanged to the printed digit, and the
outer cell falls from 4.4e-3 to 1.0e-4 over the same ladder, order 1.80 in
`1/N` (PLM) and 1.59 (WENO3). The largest residual of the domain is then the
grid junction, not the boundary, and the largest velocity one CFL step of the
operator produces falls from 1.32e-4 `c_s` (47 cm/s, at the outermost cell,
grid independent) to 4.8e-6 `c_s` at N = 250 and 2.4e-7 at N = 2000 (1.7 and
0.084 cm/s), located at the junction and converging at order 1.44 in `1/N`.

### 12. `free_outflow_boundary.f90` -> `free_outflow_boundary`

| | |
|---|---|
| Origin | item B4-4a, the outer ghost of `Apply_BC.f90`; B4-4 of `docs/b4_spatial_operator_design_20260906.md` section 2.5 (the outer face imposes no condition and the state is checked to be outflowing) |
| Quantity | (a) the production `Num_flux` at one face, called twice with two different right states; (b) the ghost `Apply_BC` writes at the outer boundary; (c) the two face states of the outermost physical cell, which differ by that cell's limited slope |
| State | a smooth isothermal radial outflow on the same geometry as test 11, `v` ramping linearly from 0.2 to 2.5 `c_s` across the domain with `rho` from steady continuity, so that the outer face is supersonic (Mach 2.45) and the base face subsonic (Mach 0.20) |
| Reference | (a) bitwise equality at the supersonic face: HLLC clamps `SL = min(0, ...)`, so with every wave speed positive it returns `Phys_flux(WL)` and never forms an expression in the right state; a nonzero change at the subsonic face, which is the discriminating half; (b) `rho > 0`, `p > 0`, `v_ghost v_N >= 0` (a reversed ghost velocity would make the boundary an inflow, which needs a reservoir this closure has none of), and `p/rho` of the ghost equal to cell N's, so the ghost carries cell N's temperature; (c) a nonzero slope |
| Tolerance | 0 (bitwise) on the supersonic flux; 1e-14 relative on the ghost temperature (round-off of one division); 1e-12 relative on the two inequalities |
| Configurations | PLM and ESWENO3, HLLC: 18 assertions |
| Expected at HEAD | **GREEN**. `outer_cell_carries_a_pressure_slope[PLM]` was RED before B4-4a (measured 1.4e-15, the round-off of a reconstruction whose limited slope is exactly zero) and measures 3.6e-2 with the continuation |

The supersonic assertion is about the SOLVER and not about the stencil, and
the entry above says so because the difference is easy to lose: a supersonic
outflow carries no information inward through the Riemann solution, but the
ghost still reaches cell N through that cell's reconstruction stencil. A
change of the outflow ghost therefore moves a supersonic run, and by the size
of the reconstruction's own slope; what it cannot do is enter through the
upwinding.

### 13. `coupled_carrier_h2_row.sh` -> `coupled_carrier_h2`

| | |
|---|---|
| Origin | item B5i; the decision of `docs/To_be_determined_by_user_20260906.md` section 10, taken 2026-09-08 (option a), on the measurement of item B5h |
| Quantity | the exit status of `EXHALE.x` and the text it prints for three key combinations resolved in `input_read.f90` |
| Reference | `Coupled carrier solve: True` with `Molecular chemistry: True` and the carriers not transported is refused with status 1 and a message naming `Molecular carrier transport: True`; the same input with that key set is not refused; and `Coupled carrier solve: True` on an atomic `He_diffusion` run, whose Newton row is an element and not H2, is not refused either |
| Why | with the carriers eliminated the layer's H2 content is not an unknown and the local-equilibrium sweep cannot determine it (the fast chemistry conserves H2 nuclei in the shielded layer), so the stationary residual is not a function of its unknowns: MEASURED seed dependence 3.2e5 times `Resid tol` per unit relative seed change, against 9.7e-10 with the H2 row carried |
| Tolerance | 0 (an exit status and the presence of five strings) |
| Configurations | copies of `examples/15_molecular` (twice) and `examples/14_diffusion`, each capped at one step, in `ccA/`, `ccB/`, `ccC/` |
| Expected before B5i | **RED** on the two stage-A assertions: the eliminated-H2 coupled configuration ran to completion and printed nothing |

### 14. `carrier_transport_inert_report.sh` -> `carrier_transport_inert`

| | |
|---|---|
| Origin | item HYG-CHECKED; `docs/input_schema.md` K15d and the user manual both stated "Needs `Molecular chemistry: True`" while nothing checked it |
| Quantity | the exit status of `EXHALE.x` and the text it prints when `Molecular carrier transport: True` is set with and without `Molecular chemistry: True` |
| Reference | on an atomic run the key draws a report naming itself, that it has no effect, and `Molecular chemistry: True` as the remedy, and the run continues to write a profile; on a molecular run the report is absent |
| Why | the transported carriers are species of the molecular network and every consumer of `carrier_transport` is guarded by `thereis_mol`, so with the network off the key changes nothing and the state produced is the correct atomic one; what is wrong is the input file, which is a report and not a refusal (the treatment `Molecular IR bands` and `Stellar LW flux` already get) |
| Tolerance | 0 (an exit status and the presence of four strings) |
| Configurations | copies of `examples/14_diffusion` and `examples/15_molecular`, each with the key appended and capped at one step, in `ctiA/`, `ctiB/` |
| Expected before HYG-CHECKED | **RED** on the first two assertions: the atomic run with the key set ran to completion and printed nothing about it |

### 15. `restart_round_trip.sh` -> `restart_round_trip`

| | |
|---|---|
| Origin | PLAN_20260909_rev1 item N9, review F13 and ISSUES 3.6: a restart of a molecular state was refused with a message whose two He/H values printed identically, and no test separated serialization from reconstruction |
| Quantity | six statements about a write -> read -> write of the hot-Uranus molecular gate: the stored scalar fields `r`, `v`, `p`; the species columns and the inventories the loader rebuilds (H, He, O, C nuclei, free charge, He/H, H2 mixing ratio); the mass density against the mass policy of the species columns; the maximum of each row of the steady residual and the cell it sits in, before and after the trip; the second trip against the first; and the cell the He/H refusal names |
| Reference | `r`, `v`, `p` bitwise; state and inventories within a round-off budget; the loader's `rho` equal to the mass policy of the species it read; the residual maxima and their cells unchanged; the second trip no farther than the first; the refusal naming the cell whose He/H decided it |
| Why | the writer emits dimensional text and the loader rebuilds `rho` from the species with the run's mass policy, so bitwise identity of every column is the wrong gate: what must hold is that nothing but round-off enters the conserved state and the inventories, and that no chemistry, remap, boundary rebuild or clock advance happens on load. The writer's `rho` column is not read at all, and the distance between it and the mass policy of its own species columns is the writer state's mass-closure defect, printed as context (MEASURED 4.1e-14 on this fixture) |
| Tolerance | 0 for the stored fields and for the residual cells, 1e-14 (about forty-five units in the last place) for the reconstruction, 1e-6 for the residual maxima |
| Configurations | six copies of `backup/regression/roundtrip` in `restart_round_trip/A..F`: a 40-step cold start, two `EXHALE_DUMP_IC=1` reloads, two `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0` evaluations, and one reload of a state whose helium at physical cell 199 has been raised by half |
| Expected before N9 | **RED** on `restart_refusal_names_the_deviating_cell` only: the message printed the TOP cell's ratio, which on a state whose He/H varies with radius agrees with the input to four digits while a front disagrees by percent, so the message said the file was written at the composition the input asks for and then refused it |
| Diagnostic, not asserted | the window integrals of the residual. The assembled residual of a relaxation snapshot is not a Lipschitz function of the state in the odd-even band above the base: MEASURED, one unit in the last place on every species column of this gate moves the 1.03-1.10 mass window by 4 percent and flips the sign of the cell residual at r = 1.075, while every row maximum and its cell stay where they were |

### 16. `restart_intent_and_metadata.sh` -> `restart_intent`

| | |
|---|---|
| Origin | PLAN_20260909_rev1 item N10 and `docs/restart_contract_design_20260909.md` (decision 15, option a): a restart said nothing about what it was continuing, and the state files said nothing about the configuration their state solves |
| Quantity | nineteen statements: the three input refusals of the intent selector; that `stationary evaluate` takes no step, returns the stored fields and the mass density within their budgets, returns the rebuilt quantities within the chemistry solve's own floor, writes back the certification it measured and the step its coupling was armed at; that the metadata block round trips field for field; that a pair with no block is loaded and marked `provenance unknown` in the log and in the files the run writes; that a disagreeing `reservoir` (a ratio and an element set), `grid`, `constants` or `options` field refuses the load naming the field, and so does a pair whose halves carry different blocks; and that the stationary residual of one state is the same number by three routes |
| Reference | the refusals taken with the field or the key named; 0 accepted and 0 attempted steps; the residual of the evaluation equal to the `EXHALE_RESIDUAL` diagnostic's to every printed digit and to the first rows the stationary solve judges to the digits it prints |
| Why | a stationary state is a state the equations hold on, and until this item every path into the solver went through the marching loop, which takes at least two CFL steps before it hands off: what was certified was never the state the file carried. MEASURED on the atomic element-row reload: the state as loaded stands at `\|\|R\|\|` 1.8924, the marched hand-off reported 1.888 for the state two steps later (`backup/regression/atomic_elem_newton/README.md`) |
| Tolerance | 0 for every refusal, for the step count and for the block comparison; 1e-14 for the stored fields (`r`, `v`, `p`), 1e-12 for the mass density (the writer state's own mass closure), 1e-8 for the temperature, species, heating and cooling, which come from one equilibrium sweep whose cell solve stops at `xtol = sqrt(eps)` ~ 1.5e-8 |
| Configurations | thirteen short runs: `A` and `B` evaluate `backup/regression/wasp_full_newton`'s own written state and then `A`'s output (the block round trip); `reservoir`, `grid`, `constants`, `options`, `elements`, `pair` each change ONE field of that block; `e_noload`, `e_traj`, `e_nosolver` are the three contradictory inputs; `x_eval`, `x_resid`, `x_solve` measure `backup/regression/atomic_elem_newton/IC` by the three routes. Nothing in `backup/` is written to |
| Expected before N10 | **RED** throughout. MEASURED with the entry binary: `Restart intent` is an unrecognized line (a warning, no refusal, exit 0 or the certification's own 2), the restart marches instead of measuring the state (1 accepted of 1 attempted step under `EXHALE_MAXSTEPS=1`), and the files it writes carry none of the eight block lines, so the round trip, the marking and the refusals have nothing to compare |

### 17. `restart_option_change.sh` -> `restart_option_change`

| | |
|---|---|
| Origin | PLAN_20260909_rev1 item N10b and decision 21 of `docs/To_be_determined_by_user_20260906.md`, option a: the restart contract of N10 refuses any difference in the `options` field, which forbids the option ladder this code is converged with (converge without an option, restart with it on, converge again) |
| Quantity | eleven statements: a rung whose named option differs loads, says so, and writes the change into both halves of the state it produces as one `# option_change` line; the same difference with nothing named refuses the load and states the TOKEN that differs, not only the field; a named token that does not differ is reported and nothing else, and writes no line; the key stops the run on an unknown token, on a token that decides how many unknowns the state has (`mol`), on a word that names the grid instead of an option (`N`), and without `Load IC? True`; a rung inherits the `# option_change` lines of the rung it was restarted from; and `EXHALE_setup.out` states both what was known about the loaded state's configuration and what this run was allowed to change |
| Reference | the refusals taken with the token or the key named; one `# option_change` line after rung 1 and two after rung 2, the second naming `cond=F -> cond=T` |
| Why | a state written under one option set is not a solution of another, which is why N10 refuses the difference; but it is a legitimate STARTING POINT for it, and that is the whole ladder workflow. The key makes the change explicit and auditable instead of silent (option b) or impossible (option c) |
| Tolerance | 0: every row is a refusal taken or not taken, a line present or absent, or a count |
| Configurations | nine short runs on copies of `backup/regression/atomic_elem_newton`, whose `IC/` pair is pinned and carries no configuration block, so rung 0 is the provenance-unknown load that produces the block-carrying state the rest start from. Every run takes `Restart intent: stationary evaluate`, so none of them marches or solves; the one cold run (`noload`) is bounded to a single step with no solver. Nothing in `backup/` is written to |
| Expected before N10b | **RED**, all eleven rows, MEASURED with the entry binary `EXHALE_lwv.x`: `Restart option change` is an unrecognized input line (a warning, no refusal), so both the named rung and its unnamed twin are refused by the exact comparison of the whole field, which names the field and not the token; no `# option_change` line is written and neither setup report line exists |

### 18. `direct_steady_setup_report.sh` -> `direct_steady_setup`

| | |
|---|---|
| Origin | PLAN_20260909_rev1 item N28, reported by N10: the direct steady route (`EXHALE_PTC=1`) writes no setup report |
| Quantity | `EXHALE_setup.out` of a run taken on the direct steady route: that it is written at all, and that it carries the report's header line and its reconstruction-method line |
| Reference | a non-empty file carrying both lines |
| Why | the setup report is the record of the configuration a run resolved, and every reader of a run directory attributes its outputs to it. The route solves and stops inside its own branch and never reaches the marching route's call, so what it left was the empty file the unit is opened with |
| Tolerance | 0: the file is written or it is not, the line is there or it is not |
| Configurations | one short run on a copy of `backup/regression/roundtrip` with `EXHALE_PTC=1`, its JFNK solve bounded to one outer iteration (`EXHALE_PTC_JFNK=1`, `EXHALE_JFNK_MAXIT=1`). The report is written before the solve, so one iteration reaches the point under test. Nothing in `backup/` is written to |
| Expected before N28 | **RED**, MEASURED with the entry binary `EXHALE_lwv.x`: `EXHALE_setup.out` has 0 lines; with the call on this route it has 66 |

## Working directories

`run.sh` creates these under `build/tests/grid_and_gates/`, and every one of
them is removed and rebuilt on each run:

| directory | test |
|---|---|
| `photon_run/` | the photon quadrature driver runs here so that no `opacity.inp` is in reach |
| `flux_spread/` | the `mol_base_handoff` copy of test 3 |
| `mbh/` | the `mol_base_handoff` copy of test 4 |
| `lp/` | the `lower_profile` copy of test 5 |
| `rgA/`, `rgB/`, `rgC/` | the three runs of test 6 |
| `sedA/`, `sedB/`, `sedC/` | the three runs of test 8 |
| `momentum_row/` | the `mol_base_handoff` copy of test 9 |
| `ccA/`, `ccB/`, `ccC/` | the three runs of test 13 |
| `ctiA/`, `ctiB/` | the two runs of test 14 |
| `restart_round_trip/A` .. `F` | the six runs of test 15 |
| `restart_intent_and_metadata/A`, `B`, `reservoir`, `grid`, `constants`, `options`, `elements`, `pair`, `e_noload`, `e_traj`, `e_nosolver`, `x_eval`, `x_resid`, `x_solve` | the thirteen runs of test 16 |
| `restart_option_change/rung0`, `rung1`, `rung2`, `unnamed`, `inert`, `unknown`, `layout`, `gridtok`, `noload` | the nine runs of test 17 |
| `hydrostatic_run/` | the four grid invocations of test 11 and their `hydrostatic_ladder.dat` |
| `outflow_run/` | test 12, which reads no input file and writes none |
| `threshold_run/` | the threshold-edge driver runs here so that no `opacity.inp` is in reach |

## Summary

**Phase 0, at `35d9dd5`.** Six programs, 28 assertions: 9 pass, 19 fail
(13.4 s). Failing: 6 of the 8 of `grid_width` (the two `Uniform` ones pass,
every cell having the same width), 8 of the 12 of `photon_quadrature`, 2 of the
3 of `flux_spread`, and 1 each of `output_state`, `base_level` and
`restart_grid`. Every failure is a defect of section 10.1 or 10.2, and the
tests are meant to stay red until the item that owns each one lands.

**Phase 1 batch 2a.** Two tests were added, each green with its item:
`grid_window` (seven assertions, item 2a-GUARDS, the out-of-bounds `minloc` of
`cell_nearest_radius` and the `j_min >= 1` floor) and `sed_coverage` (six
assertions, item 2a-SED, the coverage stop of decision 17). The six Phase 0
programs were not re-run in that pass and no item of batch 2a targets any of
them: the shifted `dr_j` of `define_grid.f90` 152 and 185, the threshold-bin
quadrature (re-measured under item 2a-PLANCK at `P_HI` code/exact 1.0951 with
the sub-Lyman grid, 1.0003 without, `P_HeII` 1.0299), the `du` functional, the
forced `sec_ion_active`, the two base-level statements and the restart grid all
stand as Phase 0 measured them.

**B4-5 (2026-09-06).** `hydrostatic_residual` was added as a MEASUREMENT: the
momentum residual of the production right-hand side on the analytic
hydrostatic column, at four grids and for both reconstructions and both
Riemann solvers. Its eight assertions are green at HEAD and were written to be
(they state identities that hold), so the numbers it prints, not its verdict,
are its product; they are in section 8 of
`docs/b4_spatial_operator_design_20260906.md`. No production source was
changed for it.

**B4-4a (2026-09-06).** The outer free-outflow ghost became one rule for every
reconstruction, the isothermal hydrostatic continuation of the last physical
cell, in place of a zero-gradient copy under PLM and a linear extrapolation
under ESWENO3. `hydrostatic_residual` gained its third assertion, on the order
of the outermost cell's residual, and `free_outflow_boundary` was added.

**O0FPE (2026-09-06).** `fpe_traps` was added: it runs 50 steps of
`hydrostatic_column` with the binary at `EXHALE_EXE` and fails on a signal or
a backtrace, so a build carrying `-ffpe-trap=invalid,zero,overflow` turns it
into a statement that no statement of the run forms an undefined or
out-of-range value. It is green with the production binary and with a
trapping build at `-O3`; it is red with a trapping build at `-O0`, and not at
a trap: gfortran evaluates both operands of the `.and.` of the duty-cycle
test in `src/EXHALE_main.f90` at `-O0`, so `mod(count, n_err_every)` divides by
zero whenever the step-doubling duty cycle is off. That site and the stall
detector's `abs(du - du_prev)/max(du,1.0d-30)`, which overflows against the
`huge(1.0d0)` sentinel on the first step of each stage, are named in the test
header with their fixes.

**HYG-CHECKED (2026-09-08).** `carrier_transport_inert` was added with the
startup report it states. Every requirement sentence of `input_read.f90` and
of `docs/input_schema.md` was checked against the code that enforces it; this
was the one whose enforcement was missing and whose stated form (a
requirement) did not match what the code does or should do (a report, the
atomic run being correct). Its four assertions are green with the report and
the first two are red without it.

**P11 (2026-09-11).** `hydrostatic_residual` gained twenty assertions on the
momentum row's reference scale against the TERMS OF THE MOMENTUM EQUATION,
after the scale stopped being `max(|dF_2|, |S_2|)`: on a state with no force
at all, on the discrete hydrostatic equilibrium with the key off and on, on a
perturbation of it, and on supersonic uniform flow. Ten of the sixty
assertions of the program were RED on the entry text, all PLM (the geometric
pressure term `(A+ - A-) p_c/dV` entering the max separately from the flux
pressure it cancels against) and the four zero-gravity rows of item P4's
group (which read the scale against `|dF_2|`, the SUM of the two terms the
row holds), and all sixty are GREEN after. MEASURED with the whole suite:
188 PASS / 10 FAIL on the entry text, 198 PASS / 0 FAIL after.

**B5i (2026-09-08).** `coupled_carrier_h2` was added with the startup refusal
it states: a coupled steady solve on a molecular configuration may not
eliminate H2, because the local-equilibrium elimination does not determine the
shielded layer's H2 content and the stationary residual is then not a function
of its unknowns. Its five assertions are green with the refusal in
`input_read.f90` and the first two are red without it.

**L26 (2026-09-16).** `base_boundary_continuity_probe.f90` -> `base_continuity`
was added as the frozen-state DIAGNOSTIC of `docs/PLAN_20260916_rev3.md`
section 6. It calls `input_read` and `init` on the run directory named by
`EXHALE_L26_STATE` and then calls `base_boundary_states` on copies of that
state along the seven paths of the plan, printing `rho_b`, the two ghost
densities, `M_i`, `M_wind`, `w_i`, `s_wind`, `w_rev`, `d_window`, `have_F`,
`rho_res` and `rho_rev` at each sampled point, plus one-sided directional
derivatives at three step sizes. Those lines are measurements and carry no
verdict.

**L26 repair (2026-09-17).** The entry is no longer a diagnostic with one
invariant: it is the suite of the repaired boundary, **32 `PASS|FAIL` rows**.
Seven are zero-window limit rows, the face density at the smallest sampled
amplitude of a path against its value at amplitude zero, tolerance 1e-12 of
`rho_b`; sixteen are derivative rows, asserting that the difference quotient
does not GROW as the step shrinks, which is the one thing a jump does; one is
the extremal-cell row, asserting that the two one-sided derivatives at the
argmax switch agree to 1e-3; four are static rows; and four are the R48
identities the entry already had, that the boundary is a function of its
argument (ten repeated evaluations, an evaluation of another state in between,
a complete `W_to_U` / `Reconstruct` / `RK_rhs` evaluation in between, a changed
call order). It also prints the group `local readings of the base mass flux`,
which is a measurement and carries no verdict: every reading of the base flux
that is a function of the base cells alone, beside the Riemann flux through the
base face on the same state.

All 32 are green at HEAD on five states (the three of the item plus the
certified `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` and the fiducial
`atomic_scalar_gj1132_kzz1e9/HeH2.13`), and the limit and extremal-cell rows
are RED against the entry text of that item by 3.7e11 and 1.9e4 respectively.
The measurements are in `docs/lhs1140b_stationary_L26_20260916.md`, section
"Repair". Without `EXHALE_L26_STATE` the entry prints one DIAGNOSTIC line and
skips, since it has no state to probe.

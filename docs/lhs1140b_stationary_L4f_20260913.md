# L4f: why the energy row of the outer boundary cell refuses once the element row is inside

Item L4f of `docs/PLAN_20260913_lhs_stationary.md`, raised by section 10 of
`docs/lhs1140b_stationary_L4d_20260913.md`. This item is a DIAGNOSIS: no file
under `src/` was changed. Every number below is MEASURED on this tree unless
it is marked READ. The case is
`LHS1140b/models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`, the isolating
control of L4b, and the controls are the certified metal-free case at the same
He/H and the three certified states of L8.

---

## 1. Verdict

**It is a property of the state and of the step the solve takes, not of the
outer continuation: reading (i).** Four measurements, each on its own, refuse
the boundary explanation.

1. **Three CERTIFIED states satisfy the same row, at the same cell, under the
   same ghost rule.** The energy row of cell 500 reads 1.66e-12 on the
   certified metal-free HeH2.13 state at 30 R_p and 8.81e-12, 1.96e-11 and
   3.50e-11 on the three certified L8 states at 30, 45 and 60 R_p, against a
   tolerance of 1e-6. The 45 and 60 R_p states leave the domain SUPERSONIC (Mach 1.065 and
   1.292, READ from L8 table T2), so
   their ghost never enters the Riemann flux, and the 30 R_p states leave it at
   Mach 0.52, so theirs does; the row is the same order in both regimes. What
   the boundary states is not what stops the C/N/O case.
2. **The row is smooth four decades below its tolerance.** Perturbing the
   outermost cell's pressure and temperature by +/- 8e-6 relative, with the
   full residual including the composition sweep, moves R_3(N) along a straight
   line whose largest departure is 2.5e-20 in code units, which is 2.7e-10 of
   the row's own scale. The radiative terms move in their eighth digit. There
   is no kink, no limiter switch and no rootless chemistry at that cell.
3. **The row's own diagonal block is O(1), and the ghost REINFORCES it.** The
   finite-difference column dR_3(N)/du(k,N) is (-5.82, +6.56, +2.84) in code
   units with the ghosts recomputed from the perturbed cell, which is what the
   solve linearizes, and (-2.44, +2.41, +2.28) with them held fixed: ratios
   1.24 to 2.72. The row does not chase itself; moving cell 500 moves its own
   row by more, not less, once the ghost follows.
4. **Cell 500 is not where the refusal is.** At the entry of outer pass 8 the
   energy row stands outside 1e-6 in **499 of the 500 physical cells** and
   outside 1e-4 in 371 of them; the worst cell is 157 at r = 1.0697 (1.032e-3),
   and cell 500 carries 3.371e-4, the value of a smooth plateau that the outer
   twelve cells all share to 7 percent. "The energy row of cell 500" is the
   name the solve's worst-row line puts on a column-wide imbalance.

**What cell 500 does have is amplification, and it is what makes it the row the
ladder reports.** The energy row's scale there is the heating, 4.64e-11, while
the two face energy fluxes the row differences are 3.80e-9 each: the row is a
1.2 percent cancellation, and the ratio |A+F_3/dV| / s_3 rises outward from
10 to 27 in the layer below 1.2 R_p to 81 at 29 R_p and to 111 at 43 R_p in a
wider domain, with no local maximum inside the wind. Measured
directly, the row measure at cell 500 moves by 76 per unit relative change of
that cell's pressure. The Newton step of pass 8 moves cell 500 by 1e-5 to
2.6e-5 relative at every iteration, and 2e-5 x 76 = 1.5e-3 is the swing the
outer row shows. The outer cell is the loudest reading of a step that is not
reducing any judged row, not the reason the step fails.

**What the step is actually doing, measured.** Run in isolation from the pass-7
state with the iteration cap at 60, the solve returns `info = 1` with `||R||`
1.032e-3, the state it was given, after 60 iterations and 58 non-monotone
accepts, while the line-search merit ||F/D||_2 falls from 1.69e-4 to 2.11e-5, a
factor 8. At every one of those iterations the Krylov cycle reached the
tolerance it was asked for in ONE product, and the unpreconditioned linear
residual it left, measured over the largest cell of each row of F, was 0.13 to
0.20 of the mass row, 0.98 rising to 5.1e+02 of the momentum row, and 0.32 to
0.56 of the energy row. The direction the cycle returns lowers one scalar
2-norm and leaves the rows the certification reads where they were. This is
L4d section 5 measured again at pass 8, now on the energy row as well as on
mass and momentum, and it is item L4e.

**Reading (ii) is refused.** There is no case for judging the last cell against
"what the boundary can deliver", and no rule of that shape is warranted. The
base cell's continuity anchor exists because the ARITHMETIC there cannot
resolve the row (`cert_tol_mass_at`, the rounding floor of a flux difference
formed at eps/Mach); at the outer cell the arithmetic resolves the row to
2.7e-10, four decades inside the tolerance, and four certified states stand
between 1.7e-12 and 3.5e-11. A tolerance introduced at the last cell would be a
number chosen to let a state pass.

**Reading (iii) is refused as a cure.** Solved at `Outer radius 45`, the same
case does not lose the outer-cell row: it gains it. Outer pass 1, which at
30 R_p reaches `info = 0` in 226 s, stagnates at 45 R_p after 2440 s at
`||R||` 4.84e-2 with all three rows at cell 500 (mass 8.30e-3, momentum
4.05e-8, energy 4.84e-2), and cell 500 holds the worst row in 306 of the 343
iterations that reported one. The amplification factor at the new outermost
cell is larger, not smaller (111 to 138 against 81). Section 6 states the
caveats of that test.

---

## 2. The energy row at the outer cell, term by term

`assemble_residual` (steady_residual.f90) builds

```
R_3(j) = dF_3(j) - S_3(j) - (heat(j) - cool(j))
dF_3(j) = [ A+ F_3(j) - A- F_3(j-1) + dF3p ] / dV
dF3p    = A+ F_1(j) (phi_i(j) - phi_c(j)) - A- F_1(j-1) (phi_i(j-1) - phi_c(j))
```

with A+ = r_edg(j)^2, A- = r_edg(j-1)^2 (RK_rhs.f90). The three pieces below
are named `out` = A+F_3(j)/dV, `in` = A-F_3(j-1)/dV and `grav` = dF3p/dV: the
outgoing and incoming enthalpy-plus-kinetic energy flux and the work against
gravity that RK_rhs folds into the same row. `S_3` is identically zero here and
so is the conduction source: the case sets neither `Viscosity:` nor
`Conduction:`, so `transport_active()` is false and no outer-face conduction
flux exists (viscous_conduction.f90:314-328). The row's scale is
`max(|dF_3|, |S_3|, heat, cool)` (`energy_row_scale`).

**T1. Signed terms at cell 500, r = 29.0312, code units**

| state | out | in | out - in | grav | heat | cool | dF_3 | R_3 | s_3 | R_3/s_3 |
|---|---|---|---|---|---|---|---|---|---|---|
| metal free, CERTIFIED, 30 R_p | 7.706725e-9 | 7.709863e-9 | -3.138e-12 | 9.599535e-11 | 9.322192e-11 | 3.648918e-13 | 9.285703e-11 | 1.545e-22 | 9.322192e-11 | **1.658e-12** |
| C/N/O, campaign state (entry of the L4d reproduction) | 4.274632e-9 | 4.274559e-9 | +7.3e-14 | 5.286450e-11 | 4.987150e-11 | 1.403067e-13 | 5.293745e-11 | 3.206260e-12 | 5.293745e-11 | **6.057e-2** |
| C/N/O, entry of outer pass 8 | 3.800813e-9 | 3.802518e-9 | -1.705e-12 | 4.799008e-11 | 4.639988e-11 | 1.305006e-13 | 4.628502e-11 | 1.564139e-14 | 4.639988e-11 | **3.371e-4** |

**What the 3.371e-4 is a fraction of**, at the entry of pass 8:

| against | value |
|---|---|
| the heating, which is the row's scale | 0.0337 percent |
| the work against gravity | 0.0326 percent |
| the enthalpy-plus-kinetic flux divergence `out - in`, the term the row has to trim | **0.92 percent** |
| the outgoing face flux `out` itself | 4.1e-6 |

The last line is the demand the tolerance makes: two face energy fluxes of
3.80e-9 each must differ by the right amount to within 4.6e-17, a relative
precision of 1.2e-8 between them. The certified metal-free state reaches
2.0e-14 of `out`, so the demand is met with six decades to spare where the
state is right.

**The physics the row holds at that radius.** `grav` and `heat - cool` are
equal to within 4 percent and `out - in` is the small remainder: the wind's heating
pays almost exactly for the work against gravity, and what is left sets
dv/dr. That is the same statement as L8's numerator N vanishing near the outer
edge. The row is therefore the discrete form of the wind's marginal balance,
which is why it is the deepest cancellation of the column (right panel of the
figure).

**How the outer face flux is formed.** `free_outflow_ghost` (Apply_BC.f90:192)
writes the ghost cell averages as the isothermal hydrostatic continuation of
cell N, `rho_g = rho_N exp[-(phi_g - phi_N) rho_N/p_N]`, `p_g` by the same
factor, `v_g = v_N`; the run's own ghosts, MEASURED at the entry of pass 8:

| cell | r [R_p] | rho | v | p | p/rho |
|---|---|---|---|---|---|
| 500 | 29.0312 | 4.700305e-10 | 8.102555e-1 | 8.455398e-10 | 1.798904 |
| ghost 501 | 29.5156 | 4.546012e-10 | 8.102555e-1 | 8.177840e-10 | 1.798904 |
| ghost 502 | 30.0000 | 4.401526e-10 | 8.102555e-1 | 7.917923e-10 | 1.798904 |

so the ratio is 0.96718 per ghost cell and `p/rho` is held to every printed
digit, as the routine states. The reconstructed states meeting at the outer
face are

```
WL(500) = (4.622359e-10, 8.102983e-1, 8.312954e-10)
WR(500) = (4.622479e-10, 8.102555e-1, 8.313156e-10)
```

a jump of 2.6e-5 in density, 5.3e-5 in velocity and 2.4e-5 in pressure. The
HLLC dissipation that jump carries is part of the discrete operator and the
converged state is defined with it: the certified metal-free state has the same
kind of jump (3.7e-5 in density) and its row still closes to 1.66e-12.

**Does the row chase itself?** No. MEASURED finite-difference column of
R_3(N), heat and cool held, transport off, step 1e-6 relative:

**T2. dR_3(N)/du(k,j), code units**

| state | j | k | ghost follows | ghost frozen | ratio |
|---|---|---|---|---|---|
| C/N/O, entry of pass 8 | 499 | 1,2,3 | +6.625, -7.351, -3.356 | identical | 1.000 |
| C/N/O, entry of pass 8 | 500 | 1 | -5.824 | -2.442 | 2.385 |
| | 500 | 2 | +6.558 | +2.414 | 2.717 |
| | 500 | 3 | +2.836 | +2.281 | 1.243 |
| metal free, CERTIFIED | 500 | 1,2,3 | -6.456, +6.615, +3.019 | -2.835, +2.573, +2.350 | 2.277, 2.571, 1.285 |
| C/N/O, campaign state | 500 | 1,2,3 | -6.265, +6.682, +2.946 | -2.690, +2.529, +2.331 | 2.329, 2.642, 1.264 |

The ghost is a function of cell N alone, so the column at cell N-1 is the same
either way; at cell N the ghost's response has the SAME sign as the direct one
and adds 24 to 172 percent to it. Removing R_3(N) = 1.564e-14 by moving that
cell's total energy density alone would take 5.5e-15, a relative move of
3.9e-6.

---

## 3. Where the energy row actually stands, cell by cell

**T3. The energy row measure over the column** (`EXHALE_RESIDUAL=1` on each
state as written; the certified states are read from their own output files)

| r [R_p] | metal free, CERTIFIED, 30 R_p | metal free, CERTIFIED, 45 R_p | C/N/O, campaign state | C/N/O, entry of pass 8 |
|---|---|---|---|---|
| 1.05 | 4.56e-10 | 7.17e-10 | 4.63e-2 | 9.44e-4 |
| 1.20 | 4.70e-11 | 2.98e-11 | 9.82e-2 | 6.14e-4 |
| 2.00 | 1.74e-11 | 5.52e-12 | 1.11e-2 | 1.89e-5 |
| 4.01 | 3.38e-12 | 1.87e-12 | 2.35e-2 | 1.02e-4 |
| 7.95 | 5.08e-12 | 3.33e-12 | 4.31e-3 | 4.51e-5 |
| 11.94 | 4.79e-12 | 7.07e-12 | 2.88e-2 | 1.72e-4 |
| 15.97 | 8.64e-13 | 8.94e-12 | 4.23e-2 | 2.43e-4 |
| 20.10 | 8.32e-12 | 1.29e-11 | 5.08e-2 | 2.88e-4 |
| 24.14 | 4.66e-12 | 8.25e-12 | 5.63e-2 | 3.16e-4 |
| 28.07 | 1.39e-11 | 1.95e-11 | 5.99e-2 | 3.34e-4 |
| outermost cell | 1.66e-12 (29.03) | 1.96e-11 (43.44) | 6.06e-2 (29.03) | 3.37e-4 (29.03) |
| worst cell of the column | 3.56e-8 at cell 16 | 2.46e-8 at cell 46 | 1.12e-1 at cell 245 | 1.03e-3 at cell 157 |

At the entry of pass 8, 499 of 500 cells stand above 1e-6 and 371 above 1e-4,
over 1.0002 to 29.0312 R_p. The outer twelve cells read 3.16e-4 to 3.37e-4, a
smooth ramp; nothing distinguishes the last one.

**T4. The amplification of the energy row, |A+F_3/dV| / s_3**

| r [R_p] | metal free 30 R_p | metal free 45 R_p | C/N/O pass-8 entry |
|---|---|---|---|
| 1.20 | 40.9 | 38.6 | 9.9 |
| 4.01 | 68.4 | 65.0 | 68.4 |
| 11.94 | 79.2 | 74.4 | 77.7 |
| 20.10 | 83.2 | 82.0 | 81.1 |
| 28.07 | 82.2 | 90.6 | 81.3 |
| outermost cell | 82.7 | 111.1 | 81.3 |

**T5. The row measure's sensitivity to the cell's own pressure**, measured on
the certified metal-free 30 R_p state by perturbing p(j) and T(j) by 1e-6
relative and reading that cell's own row

| cell | r [R_p] | d(measure)/d(ln p_j) | |A+F_3/dV| / s_3 |
|---|---|---|---|
| 100 | 1.025 | 2.64e+06 | 52.0 |
| 200 | 1.149 | 4.25e+04 | 34.7 |
| 300 | 1.857 | 4.35e+03 | 65.2 |
| 400 | 5.904 | 2.51e+02 | 69.6 |
| 450 | 12.726 | 9.40e+01 | 80.0 |
| 480 | 20.782 | 6.64e+01 | 83.2 |
| 495 | 26.694 | 6.04e+01 | 82.6 |
| 500 | 29.031 | 7.59e+01 | 82.7 |

In the wind the sensitivity IS the amplification of T4 to within 10 percent; in
the base layer it is four decades larger, so the outer WIND is the least
sensitive part of the column, not the most, and the outermost cell is
unremarkable within it. At cell 500 the 1e-6 tolerance
means locating that cell's pressure to 1.3e-8 relative, and the pass-8 stall at
3.4e-4 means it stands 4.4e-6 away.

**T6. Smoothness of R_3(500)**, thirteen states differing only in p(500) and
T(500) by -8e-6 to +8e-6 relative, full residual including the composition
sweep

| state | slope dR_3/d(eps) | largest departure from the straight line | in row units |
|---|---|---|---|
| metal free, CERTIFIED | 7.077e-9 | 2.52e-20 | 2.7e-10 |
| C/N/O, campaign state | 3.991e-9 | 1.41e-20 | 2.7e-10 |

`heat` and `cool` move in their eighth digit across the whole scan.

---

## 4. What pass 8 is doing

Reproduction of passes 1 to 7 from this case's own seed (section 7) followed by
one isolated solve from the state they end on, with `EXHALE_LINEAR_ROWS=1`,
8 threads, iteration cap 60.

**The reproduction of passes 1 to 7 is the L4d ladder.** The gated species row
of passes 1 to 7 is 8.82e-4, 2.86e-4, 1.06e-4, 4.14e-5, 1.96e-5, 1.15e-5,
6.37e-6 at cells 250, 253, 309, 313, 217, 217, 283, and pass 7 ends
`info = 2` with the energy row 1.03e-3: every one of those is the number L4d
section 4 records, to every printed digit. The hydrodynamic rows differ at the
1e-8 level (pass 1 energy 4.62e-8 here against 2.72e-8 there), which is the
tree having moved to `43bc28c` since L4d built its binary; nothing in this item
turns on it.

**The isolated pass-8 solve.**

| entering | after 60 iterations |
|---|---|
| `\|\|R\|\|` 1.032e-3 | `\|\|R\|\|` 1.032e-3, `info = 1`, THE ENTRY STATE RETURNED |
| merit `\|\|F/D\|\|_2` 1.69e-4 | 2.11e-5, a factor 8 |
| | 46 damped Gauss-Newton escapes, 58 non-monotone accepts |

**T7. What the linear solve leaves, and what the step does to the outer cell**

| iteration | linear residual left over the row of F (mass, momentum, energy) | at cell N, over F of cell N | dY/u at cell N (mass, momentum, energy) | `\|\|R\|\|` after |
|---|---|---|---|---|
| 1 | 1.139e+3, 4.443e+3, 3.613e-1 | 4.8, 22.8, 1.1e-3 (mass and momentum at their floors, 6e-24) | -3.0e-7, -2.3e-6, -1.1e-5 | 5.799e-4 |
| 3 | 1.708e-1, 1.005e+0, 4.375e-1 | 2.6e-3, 2.6e-3, 2.7e-3 | -5.3e-6, +7.4e-6, -6.4e-6 | 3.618e-3 |
| 5 | 1.713e-1, 1.009e+0, 4.387e-1 | 2.9e-3, 2.9e-3, 2.8e-3 | -9.5e-6, +1.3e-5, -1.1e-5 | 3.393e-3 |
| 9 | 1.719e-1, 9.997e-1, 4.373e-1 | 3.5e-3, 3.4e-3, 3.4e-3 | -1.5e-5, +2.0e-5, -1.8e-5 | 3.067e-3 |
| 30 | 2.001e-1, 1.546e+2, 5.556e-1 | | | 3.489e-3 |
| 60 | 2.034e-1, 5.130e+2, 5.632e-1 | | | 2.720e-3 |

Read it in that order. The cycle reports its tolerance reached in one product
of the forty available at every iteration, and what it leaves in the rows the
certification reads is 20 percent of the mass row, 56 percent of the energy row
and, by iteration 60, five hundred times the momentum row. From iteration 3 on, at cell N itself
the cycle is two to three decades better than the domain maximum (0.26 to
0.35 percent of that cell's F), so nothing is wrong with the linear solve
THERE. The step
nevertheless moves cell N by 1e-5 to 2.6e-5 relative each iteration, and the
outer row measure swings between 3.0e-3 and 7.5e-3 in step with it, which is
the 76-fold sensitivity of T5 applied to a 2e-5 step (2e-5 x 76 = 1.5e-3). The
worst-row line then names cell 500 or cell 499 at almost every iteration, which
is how the outer boundary came to be suspected.

So the pass-8 grind is not a boundary condition that cannot be satisfied. It is
a Newton direction that satisfies a scalar 2-norm tolerance while leaving the
judged rows untouched, seen through the cell that magnifies a step the most in
the wind.

---

## 5. Is there a rule in the certification for the last cell, and should there be

READ from `src/modules/time_step/certification.f90`:

- The ONLY cell-dependent hydrodynamic tolerance is the continuity row's,
  `cert_tol_mass_at` (line 782): `tol(j) = max(3e-12, min(1, 10 x floor(j)))`
  with `floor(j)` the measured rounding floor of that cell's own flux
  difference. Its justification is arithmetic: in a nearly hydrostatic layer
  the mass row is a difference formed at eps/Mach of the flux, so a fixed
  3e-12 asks one state for a tenth of an ulp.
- The species rows carry two radius regimes (`cert_tol_element_at`,
  `cert_tol_carrier_at`): reported only below 1.20 R_p, gated at 1e-5 above it.
- `hydro_row_entry` (line 1407) gives EVERY physical cell the same
  `cert_tol_momentum = 1e-8` and `cert_tol_energy = 1e-6`. There is no rule for
  the base cell of those two rows and none for the last cell.

**No rule for the last cell is warranted by anything measured here.** The
analogy with the base cell fails at the only place it would have to hold: the
base anchor exists because the assembly cannot RESOLVE the row there, and at
the outer cell the assembly resolves it to 2.7e-10 (T6), four decades inside
the tolerance, while four certified states stand at 1.7e-12 to 3.5e-11 (T1,
section 1). A tolerance introduced at cell 500 would be a number chosen so that
a state that is 4.4e-6 from the root in that cell's pressure could be called
certified.

One caveat belongs beside that, and it is L8's, not this item's: a satisfied
energy row at 29 R_p says the discrete equation balances, not that the profile
there is physical. L8 measured the 30 R_p outer profile as 25 percent too hot
and 29 percent too dense at 28 R_p against the converged answer, and
`Kn_bulk` passes 0.1 at 7.65 R_p (both READ from L8), so nothing above that radius is a validated
continuum result at any outer radius.

---

## 6. The same case at 45 R_p

Method as L8 section 3.1: a cold three-step run of the case at
`Outer radius 45` builds the target grid, `map_state_to_grid.py --ic
--extrapolate-beyond 29.031224061195065` maps the campaign 30 R_p state onto
it, and the binary is run with `EXHALE_PTC_DTAU0=1.0` at 8 threads. Stopped at
51 minutes.

| | 30 R_p (L4d control, same seed) | 45 R_p |
|---|---|---|
| outer pass 1 | `info = 0`, 226 s (this reproduction; 185 s in L4d's) | **`info = 2`, 2440 s** |
| `\|\|R\|\|` at the end of pass 1 | 2.72e-8 (energy row) | 4.844e-2 |
| rows at the end of pass 1 | mass 1.3e-9, momentum 1.4e-14, energy 2.7e-8 | mass 8.30e-3, momentum 4.05e-8, energy 4.84e-2 |
| gated species row | 8.82e-4 at cell 250 | 2.05e-4 at cell 290 |
| worst row, over the reported iterations | cells 250 to 500 | **cell 500 in 306 of 343** (mass 210, energy 57, momentum 39) |

**The cell-500 row does not disappear; it becomes the whole obstruction.**
Pass 2 was under way at 114 iterations when the run was stopped.

Conditions this test does NOT hold fixed, stated because they limit what it
proves: the seed is an extrapolation above 29.03 R_p and not a solution, the
grid is the same 500 cells over a 1.5 times wider domain (477 of them below
29 R_p), and the seed's outer face is SUBSONIC (Mach 0.42 at 43.4 R_p, MEASURED
from the mapped file), so this run never entered the regime the converged
45 R_p metal-free state of L8 is in. What the test does establish is that a
wider domain is not a cure for the row: the row is there at the new outermost
cell from the first iteration, and the amplification factor at that cell is
larger (T4).

The independent measurement that closes the boundary question is the CERTIFIED
45 and 60 R_p states of L8 (metal free): both leave the domain supersonic, so
the ghost never enters the Riemann flux (`Num_Fluxes.f90` clamps
`SL = min(0, ...)`), and their energy rows at cell 500 read 1.96e-11 and
3.50e-11, the same order as the subsonic 30 R_p state's 8.81e-12.

---

## 7. Reproduction

The measurements were made from an ISOLATED copy of the tree, `git archive
HEAD` (43bc28c) plus the working tree's `steady_newton.f90`, which is L4d's
change and is not in HEAD. The reason is concurrency, not compilation: seven
source files were being edited by other workers while this item ran
(`Cool_coeff.f90`, `diffusive_photochemistry.f90`, `certification.f90`,
`energy_semi_implicit.f90`, `composition.f90`, `ion_residual_core.f90`, the
`System_*` family), and a binary built from the live tree would have carried
unfinished physics into every number above.

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
T=$EX/scratchpad/L4f/tree
mkdir -p $T && cd $EX && git archive HEAD | tar -x -C $T
\cp -f $EX/src/modules/time_step/steady_newton.f90 $T/src/modules/time_step/
cd $T && PATH=/usr/bin:$PATH make OBJDIR=build_L4f EXE=EXHALE_L4f.x -j8

# the diagnostic binary: apply scratchpad/L4f/L4f_diagnostic.patch to $T and
#   PATH=/usr/bin:$PATH make OBJDIR=build_L4f_diag EXE=EXHALE_L4f_diag.x -j4
# It adds, to the isolated tree only and behind EXHALE_L4F=1:
#   - the signed terms of the three rows of cells N-4..N with the ghost states
#     and the face states, printed by write_residual_breakdown;
#   - output/row_measures.txt, the same quantities for every physical cell;
#   - the finite-difference column dR_3(N)/du with the ghosts recomputed and
#     with them held;
#   - in steady_newton, the linear residual and the step of cell N alone,
#     beside the existing EXHALE_LINEAR_ROWS print.

# passes 1 to 7, the state pass 8 is entered at (8 threads, 62 min here)
D=$EX/LHS1140b/models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13
W=$EX/scratchpad/L4f/pass7;  mkdir -p $W/output
\cp -f $D/input.inp $D/base.inp $W/
sed -i "s|^Spectrum file: .*|Spectrum file: $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" $W/input.inp
\cp -f $D/output/Hydro_ioniz_IC.txt $D/output/Ion_species_IC.txt $W/output/
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=7 $T/EXHALE_L4f.x > run.log 2>&1

# the terms of that state (seconds)
W2=$EX/scratchpad/L4f/res_pass8; mkdir -p $W2/output
\cp -f $W/input.inp $W/base.inp $W2/
\cp -f $W/output/Hydro_ioniz.txt $W2/output/Hydro_ioniz_IC.txt
\cp -f $W/output/Ion_species.txt $W2/output/Ion_species_IC.txt
cd $W2 && OMP_NUM_THREADS=4 EXHALE_L4F=1 EXHALE_RESIDUAL=1 $T/EXHALE_L4f_diag.x > resid.log 2>&1

# pass 8 in isolation, with the row-resolved linear residual
cd $EX/scratchpad/L4f/pass8_iso && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=60 EXHALE_LINEAR_ROWS=1 \
   $T/EXHALE_L4f_diag.x > run.log 2>&1

# the smoothness scan of T6 and the sensitivity scan of T5
python3 $EX/scratchpad/L4f/scan_smooth.py res_cno scan_cno
python3 $EX/scratchpad/L4f/scan_cell.py  res_mf  scan_cells 100,200,300,400,450,480,495,500
```

> Note added 2026-09-14 15:20: the private builds of this memo were made with `PATH=/usr/bin:$PATH make`, which selects the Ubuntu gfortran 13.1 and its LAPACK, not the conda-forge gfortran 16.2 + OpenBLAS the tree binary `EXHALE.x` is built with (bare `make`). The control and the measured build of the memo share one toolchain, so every comparison here stands; reproduce with bare `make` (`OBJDIR=... EXE=...` as written) to obtain the tree's toolchain.

The certified controls are `EXHALE_RESIDUAL=1 EXHALE_L4F=1` on the states as
written, in
`LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output/` and
`LHS1140b/models/.L8/r{30,45,60}/output/`, with the `Spectrum file:` line made
absolute and nothing else changed.

Figure: `docs/figures/lhs1140b_L4f_outer_energy_row.pdf`. Left, the energy row
measure over the column for the two certified states and the two C/N/O states,
with the 1e-6 tolerance marked; right, the amplification |A+F_3/dV| / s_3 of
the same four.

The private build (`build_L4f/`, `build_L4f_diag/`, `EXHALE_L4f.x`,
`EXHALE_L4f_diag.x`, the isolated tree) was deleted at the end of the item; the
diagnostic is kept as `scratchpad/L4f/L4f_diagnostic.patch`. No process was
left running.

---

## 8. Minimum change proposed, NOT implemented

**P1, and it is already item L4e: give the linear solve a row scaling by the
size each row is JUDGED at**, `tol_k x s_k(j)`, so the quantity the Krylov
cycle minimizes is the certification distance and not a 2-norm in which rows
ten decades apart are one number. The measurement that asks for it is T7: the
cycle reaches the tolerance it is given in one product at every iteration of
pass 8 and still leaves 20 percent of the mass row, 56 percent of the energy
row and, by iteration 60, five hundred times the momentum row; the solve
returns the state it was given after 60 iterations while the merit falls by a
factor 8. Tightening the scalar tolerance cannot reach this, for the reason
L4d section 5 already measured.

**P2: let the line search descend on the functional the ledger ranks by.** The
acceptance test is `||F/D||_2` with D the state scales, and over the 60
iterations it fell by a factor 8 while no judged row improved and 58
non-monotone steps were accepted. The judged distance is already computed at
every iteration (it is printed on the "judged rows" line) and is already what
chooses the returned state; using it, or a smooth surrogate of it, in the
acceptance test would stop the search from paying for merit with rows. This is
smaller than P1 and independent of it.

**P3: nothing in the certification.** Section 5 gives the reasons. If the
question is reopened, the numbers to hold any proposed rule against are the
2.7e-10 non-smoothness floor of the row at that cell and the 1.7e-12 to
3.5e-11 that four certified states reach there.

**P4, a print and not physics: report how many cells stand outside each row's
tolerance beside the cell that is worst.** The pass line and the JFNK
worst-row line name one cell, and on the pass-8 entry state 499 of 500 cells
are outside the energy tolerance; the single name sent this item looking at the
boundary. `certification_row_measure` already sweeps the column, so the count
costs one integer.

**Not proposed: moving the campaign's outer radius.** Section 6 measures that
45 R_p does not remove the row, and L8 already priced the change at a
1.7 percent systematic on the He I 10830 equivalent width and 13 minutes a case
against 25 seconds.

---

## 9. Noticed outside this item, reported and not acted on

- `LHS1140b/models/make_models.py` line 64 sets `Resid tol` to 5.0e-5 for the
  `scalar` group and 4.0e-6 for `scalarCNO` and `photochem`, a factor 12.5
  tighter for exactly the group that stalls. It is not the cause here (the
  pass-8 state stands at 1.03e-3 and misses both, and the metal-free case
  certified at 3.27e-8, far inside either), but the two groups are being asked
  for different things and the file records no measurement behind the
  difference.
- The energy row's scale `max(|dF_3|, |S_3|, heat, cool)` reads the ALREADY
  DIFFERENCED flux divergence, where the momentum row's scale deliberately
  gathers the terms of the equation back out of `dF_2` and `S_2` because a
  differenced quantity understates them. For the energy row that is defensible
  as it stands, since `out - in` is a genuine term of the energy equation while
  `out` alone is not, and T1 shows `grav` and `heat` within 3 percent of the
  scale. It is written here because the two rows look inconsistent until that
  is checked.

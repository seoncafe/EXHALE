# HD 189733 b: origin of the stationary cell-to-cell velocity checkerboard at the base

Investigation date: 2026-08-11. All numbers below were measured in this
investigation unless explicitly marked as an estimate. No main-tree source file
was modified; every code experiment was run from an instrumented copy of the
tree in a scratch directory (paths in §8).

---

## 1. Summary

**The controlling parameter appears to be the base density scale height measured
in cells, `H/dr`.** Where `H/dr` is small the discretization supports a
stationary, essentially undamped 2*dr* entropy mode; where it is large the same
mode decays quickly. Two independent routes to a small `H/dr` were separated
experimentally:

* coarsening the base grid at a *warm* base (T ~ 1000 K) reproduces the full
  checkerboard, and
* the collapse of the HD 189733 b base temperature to ~236 K (a factor 5 below
  `T_eq = 1183 K`) shrinks `H` by the same factor 5 on the standard grid.

So the cold base is **not required** for the mode, but on HD 189733 b it is what
produces the small `H/dr`. HD 189733 b is the planet that shows the problem
because it starts with the least margin: at `T_eq` its base scale height spans
only ~22 cells, against 53-83 for HD 209458 b, WASP-121 b and WASP-52 b (§4.1).

**The mode is an entropy (contact) mode, not an acoustic one.** In the converged
HD 189733 b solution the alternating amplitude of `ln(p/rho^gamma)` is 0.0875
while that of `ln p` is 0.00197 - a factor 44. The momentum equation is
satisfied to 1e-3 - 1e-6 of the gravity term, with *smooth* interface pressures,
so the checkerboard is not a momentum-balance failure (§6).

**The measured decay rates say the mode is essentially undamped when
underresolved.** An imposed isobaric 2*dr* perturbation of amplitude 0.03 decayed
by a factor 50 in 2000 steps at `H/dr = 135`, by 21% at `H/dr = 5.3`, and by 3%
at `H/dr = 1.7` (§5.3).

**On the base temperature (Q1), the finding is a physical-correctness problem,
independent of the checkerboard.** The base cools because 90% of the cooling at
cell 1 is `[O I]` line cooling and 10% is `C I`, applied in the optically thin
limit (`beta_esc` is hardwired to 1.0 in `util_ion_eq.f90`), while the same lines
are measured to be optically thick through the cold layer (`tau ~ 7` for
`[O I] 63 um`, `tau ~ 14` for `[C I] 609 um`, i.e. escape probabilities 0.03 and
0.013). In addition, 99.7% of the `[O I]` rate at 236 K comes from a
multi-exponential CHIANTI refit component whose stated validity range is
1e3-1e5 K, evaluated 4.2x below its floor (§3). A temperature of 236 K is also
below `T_eq = 1183 K`, which the model has no mechanism to prevent because it
contains no stellar optical/IR absorption and no thermal background - once the
XUV is shielded out, nothing sets a floor.

> [2026-08-15: a thermal background now exists as an option. `Base IR field:
> True` lets the eight ground-term fine-structure lines and the H3+ bands see
> the atmosphere below the base as a black surface radiating `B_nu(T0)` over the
> sky fraction `1 - sqrt(1 - (R_p/r)^2)`, and returns the net rate, so each of
> those channels stops cooling at its own radiative-equilibrium temperature. It
> is **off by default** (the no-incident-field limit above is what an ordinary
> run still does), and it covers only those channels -- there is still no
> stellar optical/IR absorption.]

**A direct marching test of the cooling hypothesis is not feasible.** The base
radiative time is `u/(cool - heat) = 8.3e5 s ~ 9.6 d`, against a CFL timestep of
0.77 s: about 1.1e6 steps. The production run that produced the converged state
marched 398601 steps (3.1e5 s), which is the right order; the 3000-step
experiments below cover 0.3% of it and cannot move the base temperature. This is
stated explicitly wherever it limits a conclusion.

---

## 2. The mode as measured

Converged HD 189733 b (`HD189733b/output/`), cells 1-12 above the base:

| quantity | cell 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `v` [cm/s] | -2651 | +711 | -302 | +267 | -67 | +117 | +21 | +77 |
| `T` [K] | 236.1 | 258.1 | 239.1 | 249.6 | 243.5 | 247.4 | 247.8 | 252.0 |
| `rho` [mH/cm^3] | 5.60e14 | 3.74e14 | 3.42e14 | 2.66e14 | 2.27e14 | 1.82e14 | 1.50e14 | 1.21e14 |
| `p` [cgs] | 14.66 | 12.05 | 9.873 | 8.066 | 6.583 | 5.362 | 4.366 | 3.555 |

Alternating amplitude (projection of the discrete second difference onto
`(-1)^j`, cells 1-12):

| field | amplitude |
|---|---|
| `ln(p/rho^gamma)` (entropy) | 0.0875 |
| `ln rho` | 0.0537 |
| `ln T` | 0.0443 |
| `ln p` | **0.00197** |

The pressure is smooth; density and temperature alternate in anti-phase at
essentially constant pressure. That is the signature of a stationary entropy
(contact) mode, not of an acoustic/odd-even pressure mode.

A second measured feature of the HD 189733 b base is a **density inversion at
the boundary**: `rho(cell 1)/rho(ghost) = 4.50`. The ghost is pinned to the
input base density (`Log10 lower boundary number density: 14.00`) and to
`p = ntot_bc` (i.e. `T = T_eq`), while the first computational cell holds 4.5x
that mass density at 90% of the pressure. The reported base temperature is then
simply `p/rho`: `236.1/1183.0 = 0.1996` matches
`(rho_ghost/rho_1) x (p_1/p_ghost) = 0.1995`. The effective base density of the
model is therefore not the value requested in `input.inp`.

---

## 3. Q1 - what sets `T1 = 236 K`

### 3.1 Which channels (measured, `output/Cooling_breakdown.txt`, `Heating_breakdown.txt`)

Cell 1 (`r = 1.000193`, `T = 236.1 K`, `n_e = 9.14e8`):

| cooling channel | rate [erg/cm^3/s] | share |
|---|---|---|
| O I | 3.049e-05 | 89.6% |
| C I | 3.472e-06 | 10.2% |
| Fe II | 4.16e-08 | 0.12% |
| bremsstrahlung | 2.01e-08 | 0.06% |
| **total** | **3.402e-05** | |

| heating channel | rate [erg/cm^3/s] | share |
|---|---|---|
| metals (photoionization) | 7.345e-06 | 98.3% |
| He I | 7.72e-08 | 1.03% |
| H I | 3.34e-08 | 0.45% |
| **total** | **7.474e-06** | |

So (c) is answered: H I photoionization heating is shielded out at the base
(0.45% of a total that is itself 4.6x smaller than the cooling); what is left is
metal photoionization heating. The base is **not** in local radiative balance -
cooling exceeds heating by a factor 4.55 there, and by 3-4 through cells 2-6.

### 3.2 Validity of the cooling used (measured + code reading)

**(i) The escape probability is hardwired to 1.** `util_ion_eq.f90` computes a
`beta_esc(j)` from a local optical depth, then overwrites it:

```
	beta_esc = 1.0d0
```

with the comment "Here we instead assume 100% escape (optically-thin limit)".
The `[O I] 63 um` and `[C I] 609 um` line optical depths of the cold base layer,
computed from the run's own `Ion_species.txt` densities with a Doppler core and
no turbulence:

| line | tau per cell (cell 1) | tau summed over cells 1-28 | escape probability at that tau |
|---|---|---|---|
| `[O I] 63 um` | 2.04 | 7.32 | 0.029 |
| `[C I] 609 um` | 3.90 | 13.85 | 0.013 |

These are lower bounds (down-going photons see the whole lower atmosphere).
Applying the code's own `beta_esc` formula at those depths would remove ~97% of
the base metal-line cooling.

**(ii) The dominant term is a fit extrapolated below its stated range.** The
`[O I]` coefficient is `cool_OI_ne_func = W_FS/ne + Lambda_rem(T)`, where `W_FS`
is the exact two-level `[O I] 63 um` solution (which *does* handle the
critical-density saturation) and `Lambda_rem` is the multi-exponential CHIANTI
refit of the rest. Decomposed at the base conditions:

| term | rate [erg/cm^3/s] | share |
|---|---|---|
| saturated two-level `[O I] 63 um` (`W_FS * n_OI`) | 7.95e-08 | 0.26% |
| CHIANTI multi-exponential remainder (`n_e n_OI Lambda_rem`) | 3.042e-05 | **99.74%** |
| total | 3.050e-05 | (file: 3.049e-05) |

The header of `Cool_coeff.f90` states the fits' accuracy as "0.2-2.1% over
1e3-1e5 K". At 236 K the value is set by the softest exponential of the fit,
`exp(-930.111/T)`; 930 K corresponds to no `[O I]` ground-term splitting (the
splittings are 227.7 K and 326.6 K). The same holds for `C I`, whose leading
term is `exp(-2351.38/T)` while its ground-term splittings are 23.6 K and
62.4 K. Below 1000 K these components are unconstrained by the data they were
fitted to.

**(iii) No radiative-equilibrium floor exists in the model.** EXHALE contains no
stellar optical/near-IR absorption and no thermal IR background, so nothing
prevents a shielded layer from cooling below `T_eq`. The ghost-cell pin at `T0 =
T_eq` is the model's stand-in for the lower atmosphere; cell 1 is free of it.

### 3.3 Comparison with the other planets (measured)

| run | `T_ghost` | `T1` | `rho1/rho_ghost` | `H/dr` |
|---|---|---|---|---|
| HD 189733 b (converged) | 1183.0 | 236.1 | 4.50 | 2.49 |
| WASP-52 b | 1304.0 | 181.7 | 6.89 | 1.11 |
| HD 209458 b (metals off, 2026-08-10) | 1450.0 | 385.5 | 3.63 | 6.29 |
| HD 209458 b (metals on, run of 2026-08-11) | 1450.0 | 522.7 | 2.70 | 16.80 |
| WASP-121 b | 2358.0 | 2348.8 | 1.00 | 302.28 |

Only WASP-121 b, whose base is hot enough (2358 K) that the `[O I]`/`C I`
low-temperature terms are irrelevant, keeps `T1 = T_eq`. The cold-base problem
appears wherever the base falls into the few-hundred-K regime.

### 3.4 What could not be tested here

Whether removing the excess base cooling actually warms the base **was not
measured.** Two experiments were tried and both are inconclusive for a stated
reason:

* Disabling the `beta_esc = 1.0` override (scratch build) reduced the base
  cooling from 2.90e-05 to 1.35e-08 erg/cm^3/s - i.e. removed metal-line cooling
  at the base entirely - yet after 3000 warm-start steps `T1` moved only from
  239.4 K to 239.9 K. That is exactly what the timescales predict: with the
  cooling gone the base heating time is `u/heat = 2.9e6 s`, and 3000 steps is
  2.3e3 s, so the expected rise is ~0.8 K against ~0.5 K observed. The test
  simply did not run long enough (it would need ~1e6 steps).
* Handing off to the JFNK steady solver, which has no timescale, did not help:
  it stalls on this case (`lam = 9.5e-7`, no descent), which is the known
  `steady_newton` problem under separate repair.

A quantitative expectation can still be stated: applying `beta_esc ~ 0.03` to
the metal-line channels would drop the base cooling from 3.40e-05 to ~1.4e-06,
i.e. **below** the 7.47e-06 heating, so the base would be expected to settle
substantially warmer. This is an inference from the measured breakdown, not a
run.

---

## 4. Q2 - what maintains the checkerboard

### 4.1 Grid margin at `T_eq` (computed from `input.inp` values)

With the hardcoded base spacing `drc = 2.0e-4 R_p` (`define_grid.f90`, `N_low =
50`) and a neutral H/He mean molecular weight:

| planet | `g` [cm/s^2] | `H(T_eq)` [cm] | `dr` [cm] | `H/dr` at `T_eq` |
|---|---|---|---|---|
| HD 189733 b | 2154 | 3.68e7 | 1.71e6 | **21.6** |
| HD 209458 b | 909 | 1.07e8 | 2.00e6 | 53.4 |
| WASP-52 b | 707 | 1.44e8 | 1.82e6 | 79.1 |
| WASP-121 b | 601 | 2.62e8 | 3.16e6 | 83.0 |

HD 189733 b has 2.5-4x less grid margin at the base than the other three, purely
from its higher surface gravity (Jeans parameter 192 against 61-78). Dividing
that 21.6 by the factor 5 temperature collapse leaves ~4 cells; the measured
local value in the converged run is 2.49.

### 4.2 Boundary-condition experiments (warm start, 3000 steps, default grid)

Baseline reproduces the reference state, so the comparison is like for like.

| # | change | `T_ghost` | `T1` | `H/dr` | `A(ln rho)` | `A(ln p)` | `v` sign flips /11 |
|---|---|---|---|---|---|---|---|
| ref | converged reference | 1183.0 | 236.1 | 2.49 | 0.0537 | 0.00197 | 5 |
| E0 | baseline, 3000 steps | 1183.0 | 239.4 | 2.57 | 0.0518 | 0.00196 | 5 |
| E1 | `Hydrostatic base: True` | 314.1 | 237.4 | 2.32 | 0.0195 | 0.00117 | **0** |
| E2 | ghost `T` set to `T(cell 1)` (scratch patch) | 281.6 | 281.8 | 4.52 | **0.0029** | 0.00266 | 1 |
| E3 | `beta_esc` override disabled | 1183.0 | 239.9 | 2.58 | 0.0515 | 0.00196 | 5 |
| E4 | metals off (`metals.inp` removed) | 1183.0 | 268.3 | 4.05 | 0.0397 | 0.00187 | 3 |

Reading:

* **E2 is the strongest single intervention**: removing the ghost-to-cell-1
  temperature jump (ghost pressure set so `T_ghost = T_1` at the pinned density)
  cuts the alternating density amplitude by a factor 18 and removes the density
  inversion (`rho1/rho_ghost` 4.44 -> 0.92). `H/dr` rises from 2.57 to 4.52.
* **E1 removes the velocity alternation entirely** but replaces it with a
  uniform -3200 to -6200 cm/s inflow through the whole base - the base is
  draining. After 3000 steps this is still a transient, so E1 should not be read
  as a working fix.
* **E3 and E4 are inconclusive for the reason given in §3.4** (timescale), even
  though E3 verifiably removed the metal-line cooling (`Cooling_breakdown` total
  at cell 1: 2.90e-05 -> 1.35e-08).

### 4.3 Grid-resolution experiments (cold start, 8000 steps, `drc` patched)

These start from a cold hydrostatic IC, so after 8000 steps (6.2e3 s) the base
is still near `T_eq` on every grid - the base temperature is therefore *not* a
variable here, and the only thing that changes is `dr`.

| # | `drc` | `dr` [R_p] | `T1` [K] | `H/dr` | `A(ln rho)` | `A(ln p)` | flips /11 |
|---|---|---|---|---|---|---|---|
| G1 | 5.0e-5 (`N_low = 200`) | 4.66e-05 | 1178.8 | 171.7 | **0.0001** | 0.00010 | 0 |
| G0 | 2.0e-4 (default, `N_low = 50`) | 1.93e-04 | 1088.8 | 8.5 | 0.0046 | 0.00052 | 2 |
| G2 | 8.0e-4 (`N_low = 25`) | 7.83e-04 | 973.7 | 2.0 | **0.0410** | 0.00269 | 8 |

**This is the decisive experiment.** G2 has a warm base (974 K, no density
inversion: `rho1/rho_ghost = 1.09`) and nevertheless develops the full
checkerboard, with an amplitude comparable to the converged HD 189733 b value
(0.0410 vs 0.0537) and *more* sign flips (8 vs 5). Refining by 4x removes the
mode to the 1e-4 level. The cold base is therefore sufficient but not necessary;
small `H/dr` is what the mode tracks.

### 4.4 The mode across all runs, ordered by `H/dr`

| run | `H/dr` | `A(ln rho)` |
|---|---|---|
| WASP-121 b | 302 | 0.0002 |
| G1 fine grid | 172 | 0.0001 |
| HD 209458 b (metals on) | 16.8 | 0.0216 |
| G0 default grid, cold | 8.5 | 0.0046 |
| HD 209458 b (metals off) | 6.3 | 0.0335 |
| E2 ghost `T` fix | 4.5 | 0.0029 |
| E4 metals off | 4.1 | 0.0397 |
| E0 baseline | 2.6 | 0.0518 |
| HD 189733 b converged | 2.5 | 0.0537 |
| E1 hydrostatic base | 2.3 | 0.0195 |
| G2 coarse grid | 2.0 | 0.0410 |
| WASP-52 b | 1.1 | 0.0794 |

The trend is monotone in the aggregate but not cell-exact - several of these
runs are mid-transient, so the ordering within a factor ~2 should not be
over-read. The clean statement is: amplitudes are at the 1e-4 level for
`H/dr > 100`, in the 1e-2 range for `H/dr < 5`, and the transition is somewhere
around `H/dr ~ 8-17`.

---

## 5. Q2 continued - is the mode generated, or merely not damped?

### 5.1 It is not inherited from the IC

The cold hydrostatic IC of the reference run
(`HD189733b/output/Hydro_ioniz_IC.txt`) has a smooth base: `T = 451, 471, 490,
508, 528, ... K`, `rho` ratios 0.851, 0.858, 0.863, 0.865, ... The checkerboard
develops during marching.

### 5.2 It is not maintained by a growing instability

Over 3000 marching steps from the converged state the amplitude changes by <4%
(`A(ln rho)` 0.0537 -> 0.0518). The mode is stationary, consistent with the
earlier diagnosis.

### 5.3 Decay of an imposed perturbation (the direct test)

An isobaric 2*dr* perturbation was written into a restart (`rho_j ->
rho_j (1 + 0.03 (-1)^j)` over the first 16 cells, `p` unchanged, all species
scaled together so the ionization fractions are preserved), then marched 2000
steps on each grid with the matching binary. Controls are the same restarts
without the perturbation.

| grid | `H/dr` | `A(ln rho)` at start | after 2000 steps | ratio |
|---|---|---|---|---|
| fine, `dr = 4.66e-5` | 135 | 0.0295 | 0.00046 | **0.02** |
| fine, control | 137 | 0.00009 | 0.00039 | 4.5 |
| default, `dr = 1.93e-4` | 5.3 | 0.0340 | 0.0268 | **0.79** |
| default, control | 7.2 | 0.0046 | 0.0064 | 1.4 |
| coarse, `dr = 7.83e-4` | 1.7 | 0.0704 | 0.0681 | **0.97** |
| coarse, control | 1.9 | 0.0410 | 0.0516 | 1.3 |

For scale, 2000 steps at `dt = 0.77 s` is 1.5e3 s, roughly 280 sound-crossing
times of one base cell on the default grid. An acoustic disturbance would be
gone; this one is not.

Two things follow. First, **when the stratification is well resolved the scheme
damps the 2*dr* entropy mode strongly** (factor 50 in 2000 steps). Second, **when
it is not, the mode is essentially undamped** (3-21% in the same interval) *and*
weakly generated - the unperturbed controls drift upward by 26-40% toward their
own attractor amplitude.

---

## 6. Q3 - the discretization view

### 6.1 Direct measurement of the momentum terms

`RK_rhs.f90` and `Num_Fluxes.f90` were instrumented in a scratch build to dump,
per cell and per RK stage, the interface pressures `p_out`, which side of the
Riemann problem each was taken from, the non-conservative pressure gradient
`(pR - pL)/dr`, the gravity source `S(2)`, the momentum advection, and the
residual. Warm-started from the converged state, WENO3 stage:

| j | `p` at face `j-1` | `p` at face `j` | side | `(pR-pL)/dr` | gravity `S(2)` | residual | `|res|/|grav|` |
|---|---|---|---|---|---|---|---|
| 1 | 9.5771e-01 | 8.1316e-01 | 2->2 | -7.483e+02 | -9.985e+02 | +2.40e+00 | 2.4e-03 |
| 2 | 8.1316e-01 | 6.6759e-01 | 2->2 | -7.535e+02 | -7.182e+02 | -1.13e+00 | 1.6e-03 |
| 3 | 6.6759e-01 | 5.4578e-01 | 2->1 | -6.305e+02 | -6.431e+02 | -7.69e-01 | 1.2e-03 |
| 4 | 5.4578e-01 | 4.4599e-01 | 1->2 | -5.166e+02 | -5.111e+02 | +4.89e-02 | 9.6e-05 |
| 6 | 3.6351e-01 | 2.9605e-01 | 2->2 | -3.492e+02 | -3.479e+02 | +2.84e-03 | 8.2e-06 |
| 8 | 2.4101e-01 | 1.9644e-01 | 2->2 | -2.307e+02 | -2.305e+02 | +1.97e-04 | 8.6e-07 |

Observations:

* Pressure and gravity cancel to 3-6 significant digits, confirming the expected
  4-digit cancellation. The mode is **not** a residual momentum imbalance.
* The interface pressures are smooth and monotone; there is no 2*dr* structure in
  them. The alternation lives entirely in the cell-centered `rho` and `T`.
* The Riemann-solver side selection does flip at cells 3-4, where the velocity
  changes sign - see §6.3.

### 6.2 Why no term in the scheme sees a stationary 2*dr* entropy mode

From reading the code (this part is structural, not measured):

* **Gravity source, `Source.f90`.** `S(2) = -0.5 (rhoL + rhoR) (Gphi_i(j) -
  Gphi_i(j-1))/dr`, where `rhoL` and `rhoR` are cell `j`'s *own* reconstructed
  values at its two faces. The term is entirely local to cell `j`; it introduces
  no neighbor coupling and so cannot by itself resolve or damp a 2*dr* mode.
* **Pressure gradient, WENO3 branch.** `Phys_flux` omits `p` from the momentum
  flux when `use_weno3`, and `RK_rhs` adds `(p_out(j) - p_out(j-1))/dr`
  non-conservatively. Only the two *interface* pressures enter. A cell-centered
  checkerboard whose reconstruction maps to a smooth interface sequence is
  invisible to this operator - which is exactly what §6.1 measures.
* **HLLC contact resolution.** HLLC resolves the contact wave exactly. For a
  stationary contact (`v -> 0`, as in the base) the entropy mode has zero
  characteristic speed, so the Riemann solver adds no dissipation to it at any
  wavelength, including 2*dr*.
* **No short-wavelength sink.** The model has no conduction and no physical
  viscosity, and the radiative source term is evaluated cell by cell, so it is
  diagonal and cannot mix neighbors either. The Shapiro filter is the only
  existing operator that acts on 2*dr*, and it is off by default.

  > [2026-08-15: all three clauses have changed. Explicit viscosity and
  > conduction exist as `src/modules/time_step/viscous_conduction.f90` (keys
  > `Viscosity:` / `Conduction:`), a gated fourth-difference dissipation of
  > exactly this 2*dr* contact mode exists as
  > `src/modules/flux/low_mach_dissipation.f90` (key `Low-Mach damping`), and
  > the metal line cooling is no longer cell-local: §10 made the escape
  > probability a function of the column above each cell. All three are off by
  > default, so the reading above still describes a default-configuration run.]

Under this reading the only thing that damps the mode is the reconstruction
itself, through the limiter/weights recognizing the profile as smooth - which is
precisely what fails when `H/dr` approaches unity, and precisely the dependence
measured in §5.3.

### 6.3 A specific detail worth recording

In the HLLC branch of `Num_flux`, `p_out` is **not** the HLLC star pressure. It
is the raw reconstructed pressure of whichever side the contact speed selects:

```
         if(SL.ge.(0.0)) then
            call Phys_flux(WL,NF)
            p_out = pL
         elseif(SL.lt.(0.0).and.S_star.ge.(0.0)) then
            ...
            p_out = pL
         elseif(S_star.lt.(0.0).and.SR.ge.(0.0)) then
            ...
            p_out = pR
         else
            ...
            p_out = pR
         endif
```

In the base `S_star ~ v ~ 0` and alternates in sign, so this selection flips from
face to face (measured at cells 3-4 in §6.1). Each flip changes the selected
pressure by roughly one reconstruction jump, ~18% of `p` per cell here, in a term
that has to cancel gravity to four digits. This is a plausible amplification
path for the mode once it exists; it was **not** isolated experimentally, and is
recorded as a lead rather than a conclusion. Note that `LLF` and `ROE` instead
return `p_out = 0.5 (pL + pR)`, a centered average, which is the classical
2*dr*-blind form.

### 6.4 Scheme variants (confounded - read with care)

Cold start, 8000 steps, on the coarse base grid. These runs land on *different*
base states, so the scheme and `H/dr` are not separated; the `H/dr` column is
given so the reader can see the confound.

| variant | `T1` [K] | `H/dr` | `A(ln rho)` | flips /11 |
|---|---|---|---|---|
| `PLM+WENO3` (never switched; still PLM at step 8000) | 973.7 | 2.02 | 0.0410 | 8 |
| `PLM` single stage | 973.7 | 2.02 | 0.0410 | 8 |
| `WENO3` single stage | 1236.3 | 4.51 | 0.0024 | 1 |
| `LLF` flux | 1146.9 | 7.52 | 0.0011 | 1 |
| `ROE` flux | 420.2 | 0.30 | 0.1605 | 1 |

The `PLM+WENO3` and `PLM` runs are identical, confirming the two-stage run never
switched within 8000 steps. The `ROE` run degenerated (`T` = 7887 K at cell 2)
and should not be used. Every non-degenerate row is again consistent with `H/dr`
being the controlling variable rather than the flux or reconstruction choice.

---

## 7. Candidate remedies

### Measured effect

1. **Refine the base grid.** `drc = 5e-5` with `N_low = 200` removed the mode
   (`A(ln rho) = 1e-4`) and damped an imposed perturbation by 50x. `drc` and
   `N_low` are currently hardcoded local variables in
   `src/modules/init/define_grid.f90`; they are not exposed as input keys.
   [2026-08-15: they are, as `Base grid [dr,cells]: <dr_base> [<N_low_cells>]`;
   see section 11.] Cost:
   4x the base cells, i.e. more steps at the same CFL. This is the most direct
   fix and the one with the clearest evidence, but it treats the symptom - it
   does not address the base temperature.
2. **Remove the ghost/cell-1 temperature jump.** Setting the ghost pressure so
   that `T_ghost = T_1` at the pinned density cut `A(ln rho)` by 18x and removed
   the 4.5x density inversion (E2). This is close in spirit to the existing
   `Hydrostatic base: True` key but is not the same operator; `Hydrostatic base`
   (E1) removed the velocity alternation but started draining the base.

### Inferred, not yet measured

3. **Stop applying optically thin metal-line cooling at the base.** The code
   already computes a `beta_esc` and then discards it. Applying an escape
   probability appropriate to the *line* opacity (rather than the current gray
   XUV proxy, which gives `tau ~ 4e3` and `beta ~ 2e-5` at the base and is
   equally wrong in the other direction) would reduce base cooling by roughly the
   factor 30-80 that the measured `[O I]`/`[C I]` optical depths imply. Expected
   consequence: base cooling drops below base heating, the base warms toward
   `T_eq`, `H` grows by up to a factor 5, and the mode moves into the regime
   where §5.3 shows it decays. This chain has **not** been demonstrated end to
   end; see §3.4 for why (1e6 steps needed, or a working steady solver).
4. **Guard the low-temperature extrapolation.** The CHIANTI refits are stated
   valid over 1e3-1e5 K; the base evaluates them at 236 K, where 99.7% of the
   `[O I]` rate comes from a fit component with no corresponding atomic
   splitting. Some explicit treatment below the fit floor - a documented
   extrapolation policy, or the two-level solution taken alone - would remove an
   uncontrolled term.
5. **A radiative floor at the base.** The model has no way to keep a shielded
   layer at or above `T_eq`. The Tier-1 lower-atmosphere column
   (`src/modules/lower_atmosphere/lower_column.f90`, default off) is the existing
   place where such a constraint could live.

### Not recommended on the evidence

6. **Shapiro filter.** Previously measured to remove the alternation but to
   change the base temperature structure (`T1` 236 -> 427 K) and to make the JFNK
   fail with `info = 2`. It suppresses the symptom by adding dissipation at 2*dr*,
   without addressing either the resolution or the cooling.

---

## 8. Raw data and reproduction

Scratch tree (instrumented builds, all experiment directories, logs)
(session scratchpad, no longer present):

```
/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/
  31dad923-90fb-485a-bbf3-1cf82bd83337/scratchpad/hd189_diag/
    EXHALE.x                     baseline instrumented build (base_trace.txt writer)
    variants/ghostT/             ghost pressure set to give T_ghost = T_1
    variants/betaesc/            beta_esc = 1.0 override removed
    variants/gridfine/           drc = 5.0e-5, N_low = 200
    variants/gridcoarse/         drc = 8.0e-4, N_low = 25
    variants/rkdiag/             momentum-term dump, cell by cell (rk_momentum_terms.txt)
    exp/analyze.py               alternating-amplitude metric used in every table
    exp/perturb.py               writes the isobaric 2*dr* perturbed restart
    exp/base_warm.inp            HD189733b input, Newton finish removed, du stop disabled
    exp/base_cold.inp            same, cold start
    exp/E0 E1_hydrobase E2_ghostT E3_betaesc E4_metalsoff
    exp/G0_cold_drc2em4 G1_cold_drc5em5 G2_cold_drc8em4
    exp/S1_plm S2_weno3 S3_llf S4_roe
    exp/P0..P5                   perturbation-decay runs
    exp/D1_rkterms               momentum-term dump (§6.1)
    exp/batch1.sh batch2.sh batch3.sh
```

Reference outputs read (not modified): `HD189733b/output/`,
`HD209458b/output/`, `HD209458b/output_pre_metals_20260811/`,
`WASP-121b/output/`, `WASP-52b/output/`.

## 9. Scope

Checked: the base cells (1-20) of HD 189733 b and, for comparison, of
HD 209458 b, WASP-52 b and WASP-121 b; the cooling and heating channel
breakdowns; `define_grid.f90`, `Apply_BC.f90`, `Source.f90`, `RK_rhs.f90`,
`Num_Fluxes.f90`, `Cool_coeff.f90`, `util_ion_eq.f90`, `load_IC.f90`.

Not checked: the wind region above ~1.01 `R_p`; mass-loss rates (none of these
short runs is converged, and no quantitative `Mdot` should be taken from them);
the ionization solve itself; whether the checkerboard measurably affects the
transit observables. The scratch tree carries the pre-repair `steady_newton.f90`,
so nothing here should be read as a statement about the JFNK line search.

---

## 10. Fix applied (2026-08-11)

Both physical-correctness problems of §3.2 were repaired the same day, in
`src/modules/radiation/Cool_coeff.f90` and `src/modules/radiation/util_ion_eq.f90`
(with the matching post-process paths in `T_equation.f90` and
`post_process_adv.f90`). The changelog entry is `docs/Update_EXHALE.md` §42.
The checkerboard itself (§4-§6) is NOT addressed here; only the cooling is.

### 10.1 What was changed

**(a) Line trapping.** The `beta_esc = 1.0` override, and the gray
XUV-continuum depth over one cell width behind it, are gone. `beta` is now the
line-center escape probability of the emitting line, from the column between
the cell center and the top of the domain:

* `kappa_OI63` / `kappa_CII158` give the line-center absorption coefficient of
  `[O I] 63um` and `[C II] 158um` [2026-08-15: now the single elemental
  function `fine_structure_line_opacity`, built on `line_center_opacity_lte`
  and covering all eight lines], using the SAME `A_ul`, level energies and
  ground-term partition sums the two-level emission terms use, a Doppler core
  `v_th = sqrt(2kT/m)` (no turbulence) and the stimulated-emission correction
  `1 - exp(-E/kT)` for Boltzmann level ratios.
* `fine_structure_escape` [2026-08-15: now `fine_structure_line_transfer`]
  accumulates `tau(j)` downward from the top: half of
  the emitting cell plus every cell above it. Being a column, `tau` is
  grid-independent and converges under refinement, which the cell-width depth
  it replaces did not.
* `line_escape_probability(tau)` [2026-08-15: now
  `line_escape_probability_one_face(tau)`, summed over the two faces of the
  cell instead of renormalized by a factor 2] is the plane-parallel Doppler
  form of
  Hollenbach & McKee (1979) / de Jong, Boland & Dalgarno (1980),
  `(1-e^-a tau)/(a tau)` and `1/(2 tau sqrt(ln(tau/sqrt(pi))))` with
  `a = 2.34`, **renormalized by a factor 2 so that `beta(0) = 1` exactly**.
  The published form tends to 1/2 because it counts escape through one face of
  a slab; here a photon sent downward is absorbed by the lower atmosphere,
  which the model treats as a fixed reservoir, so it leaves the modeled gas
  either way and the thin limit must be full escape. That normalization is
  also what keeps the optically thin wind at the previous `beta = 1`
  behavior. The branches are switched where they cross,
  `tau_c = sqrt(pi) exp(a^2/4) = 6.967`, so the switch is continuous in value.
* `beta` enters as `A_ul -> beta A_ul` INSIDE the two-level solution. This is
  the physically correct place: in the subcritical limit the cooling equals
  the collisional excitation rate and must not depend on `beta` (every
  excitation still ends as an escaped photon), while in the saturated limit
  the escaping flux is proportional to `beta`. Multiplying the *result* by
  `beta`, as the old assembly did, is wrong by a factor `beta` in the first
  limit.
* Scope: trapping is applied to `[O I] 63um` and `[C II] 158um` only, the two
  lines for which the code carries an explicit two-level solution. All other
  metal-line cooling keeps `beta = 1`.
  [2026-08-15: the scope is now eight lines -- `[C I] 609/370um`,
  `[C II] 158um`, `[N II] 205/122um` and `[O I] 63/145/44um` (`n_fsline = 8` in
  `Cool_coeff.f90`) -- each with its own optical depth and escape probability,
  entering the ground-term statistical equilibrium.] That is correct in the wind, and at the
  HD 189733 b base the other coolants are together < 0.3% of the total once
  (b) is applied. The thick resonance lines of a metal-rich wind (Mg II h&k)
  were recorded here as untreated; they were measured afterwards and `beta = 1`
  turns out to be the correct effective treatment for them, because their
  `A_ul ~ 1e8 s^-1` keeps the trapped-photon escape rate `beta A_ul` far above
  the collisional de-excitation rate `ne q_ul`. See
  `docs/resonance_line_trapping.md`.

**(b) Coronal-fit validity floor.** `coronal_excitation_cutoff(T)` multiplies
every CHIANTI-derived coefficient — the analytic C/N/O and Mg/Ca/Na/Fe fits,
the 1-D and 2-D tables (whose `log10 T` axis starts exactly at 3.0 and which
otherwise hold their edge value indefinitely below it), and the coronal
remainder of `cool_OI_ne_func` / `cool_CII_ne_func`, but not their two-level
parts. It is exactly 1 for `T >= 1e3 K` and `exp(-((T_floor/T - 1)/w)^2)`
below, with `w = 0.5`. Value and `dT`-slope are continuous at the floor, which
matters because the Brent energy solve and the semi-implicit update
differentiate the cooling in `T`; every coefficient is bit-identical at and
above 1e3 K. The legacy AIOLOS branch (`cno_cool 0`) is deliberately not
guarded: its constant floors are crude fine-structure stand-ins, not
extrapolated coronal fits. `w = 0.5` is a modeling choice, not a measurement.

> [2026-08-15: `w` is neither 0.5 nor hardcoded any more. The default is
> `w = 0.1` and the value is set by the input key `Coronal cutoff width: <w>`
> (`docs/coronal_cutoff_width.md`). That memo's section 7.2 also records that
> the ground-term statistical equilibrium removed what the guard was
> suppressing, so the base result no longer depends on `w`.]

### 10.2 Optical depths and escape probabilities actually reached

Computed with the code's own formula from each run's `Ion_species.txt`
(replica evaluation, for reporting; the code applies the same expression):

| run | `T1` | `tau([O I] 63um)` at cell 1 | column total | `beta` at cell 1 |
|---|---|---|---|---|
| HD 189733 b, converged reference | 236.1 | 2.73 | 3.25 | **0.156** |
| HD 189733 b, after 30000 warm steps with the fix | 256.9 | 2.01 | 2.45 | 0.210 |
| WASP-121 b, converged | 2343.2 | 0.039 | 0.039 | 0.956 |

The `tau = 7.3` / `beta = 0.029` quoted in §3.2 was computed without the
stimulated-emission correction and with a different lower-level weighting; the
correct value for this line at these temperatures is `tau ~ 3`, `beta ~ 0.16`.
The correction matters: at 236 K the `3P1`-`3P2` splitting is only 227.7 K, so
`1 - exp(-E/kT) = 0.62` already, and for `[C I] 609um` (splitting 23.6 K) it
is 0.095, which is why that line is NOT optically thick at the base
(`tau ~ 0.25` over the whole cold layer) despite the §3.2 estimate.

### 10.3 Base energy budget, measured

HD 189733 b cell 1, `output/Cooling_breakdown.txt` + `Heating_breakdown.txt`:

| | `T1` [K] | cooling | heating | net |
|---|---|---|---|---|
| before (converged reference, §3.1) | 236.1 | 3.402e-05 | 7.474e-06 | **-2.65e-05** |
| after (warm start + 5000 steps) | 240.9 | 2.681e-08 | 6.997e-06 | **+6.97e-06** |

All rates in erg cm^-3 s^-1. **The sign flips**: the base changes from cooling
at 4.6x the heating to heating with the cooling 260x below it. The reduction is
carried almost entirely by (b) — the coronal remainder was 99.74% of the `[O I]`
rate (§3.2) — while (a) alone is a factor ~5 on the two-level term.

### 10.4 Where the base now balances (measured with the code's own cooling)

The marching timescale of §3.4 makes a direct relaxation test impractical, so
the balance temperature was measured instead. A driver program linked against
the compiled `Cooling_Coefficients` object (`scratchpad/fix/probe/`; session
scratchpad, no longer present) assembles
the total metal line cooling exactly as `T_equation` does, at the FROZEN cell-1
state of the converged reference (`n_e = 9.14e8`, `n_HI = 4.15e14`, the 27 metal
densities from `Ion_species.txt`), and sweeps `T`. Two binaries were used: the
current one, and a control built from the same sources with only the two new
factors forced to 1, so the control is the old physics exactly.

| `T` [K] | cooling, old | cooling, new | new/old |
|---|---|---|---|
| 237 | 3.462e-05 | 1.244e-08 | 3.6e-04 |
| 398 | 2.706e-04 | 4.405e-08 | 1.6e-04 |
| 501 | 6.306e-04 | 1.201e-05 | 1.9e-02 |
| 708 | 1.775e-03 | 8.986e-04 | 5.1e-01 |
| 891 | 2.988e-03 | 2.815e-03 | 9.4e-01 |
| 1000 | 3.7022e-03 | 3.7021e-03 | 1.000 |
| 1189 | 4.861e-03 | 4.861e-03 | 1.000 |
| 1995 | 9.800e-03 | 9.800e-03 | 1.000 |

Rates in erg cm^-3 s^-1. The curves are identical from 1e3 K up, as designed.
Against the local heating measured at the same cell, `7.474e-06`:

| | balance temperature | cooling at 236 K |
|---|---|---|
| old physics (`beta = 1`, no cutoff) | **168.0 K** | 3.400e-05 |
| (a) line trapping only | 168.2 K | 3.394e-05 |
| (b) coronal-fit cutoff only | 485.5 K | 7.942e-08 |
| (a)+(b), as shipped | **485.8 K** | 1.241e-08 |

**(b) does essentially all of the work at the base**, because the base cooling
was carried by the coronal remainder, which is not a trapped line. (a) is
required for physical correctness and does move the two-level `[O I] 63um`
term by a factor 6.4 at 236 K, but that term is small against the heating with
or without it. The two changes are independent and both are in.

The converged run sits at 236 K, between the two, because the base is not in
local radiative balance (it also exchanges energy with its neighbors). The
2.9x rise in the balance temperature is the size of the effect. Held fixed in
this measurement: the ionization state, `n_e`, and the heating rate, all taken
from the 236 K solution; a self-consistent solve would move them, and the metal
photoionization heating in particular should rise as the metals stay neutral
over a warmer, still shielded base. So 486 K is the balance point of the
cooling change alone, not a prediction of the converged base temperature.

### 10.5 What the marching test does and does not show

Warm-started from the converged reference with the `du` stop and the Newton
finish disabled (`exp/base_warm.inp`), the fixed binary and the control were
marched side by side:

| steps | `T1` new | `H/dr` new | `A(ln rho)` new | `T1` old | `H/dr` old | `A(ln rho)` old |
|---|---|---|---|---|---|---|
| 0 (reference) | 236.1 | 2.49 | 0.0537 | 236.1 | 2.49 | 0.0537 |
| 25000 | 255.9 | 2.82 | 0.0442 | 248.7 | 2.67 | 0.0466 |
| 50000 | 274.4 | 3.18 | 0.0391 | 297.8 | 4.91 | 0.0325 |
| 75000 | 295.7 | 3.72 | 0.0347 | 310.9 | 5.15 | 0.0309 |
| 100000 | 318.2 | 4.44 | 0.0310 | 321.9 | 5.21 | 0.0300 |

(The fixed run carries an extra 5000-step lead-in, so its step labels are
5000 short of its true totals.) The chains were stopped at 100000 steps.

**This test does not resolve the radiative difference.** The control, which
still cools the base 1000x harder, warms just as fast — faster over steps
25000-75000 — and by 100000 steps the two runs agree to 4 K in `T1` and to 3%
in the alternating amplitude. The reference state was produced with a Newton
finish and is not a fixed point of pure marching, so both runs are riding a
base transient (`v` at cell 1 is an inflow of -2.6e3 cm/s) whose adiabatic
heating dominates the radiative term on this timescale — as §3.4 predicted it
would, for 1e5 steps against a radiative time of order 1e6. Nothing about the
sign or size of the cooling change should be read off these columns; §10.4 is
the measurement that speaks to it. Both runs do leave the worst of the
checkerboard behind (`A(ln rho)` 0.054 -> 0.031, `H/dr` 2.5 -> 4.4-5.2, sign
flips 5 -> 1), but they do it together, so that is the transient too.

The steady route was also tried and abandoned: warm-started with the production
`input.inp` (`Solver: Newton`), `du` rose from 1.8e-2 at step 5500 to 1.16 by
step 11000 and never reached the 1e-2 hand-off, identically for both binaries.

### 10.6 Non-regression of the wind

* WASP-121 b, full converged run (base at 2358 K, so the cutoff is inactive
  and `beta = 0.956`): `log10 Mdot = 13.17` with and without the fix; outer
  mass flux `rho v r^2` agrees to `5.3e-7` relative; largest profile
  differences `6.5e-5` (velocity, at the base cell), `6.9e-6` (density),
  `9.1e-6` (temperature). The JFNK still reaches `info = 0`.
* `make check`: `mol_base_handoff` (metals off) is byte-identical, PASS on both
  files. The two metal cases fail by the intended physics, and the shift is
  small:

  | case | `log10 Mdot` | outer mass flux new/gold | largest profile difference |
  |---|---|---|---|
  | `wasp_full` | 13.22 (unchanged) | 1 + 2.0e-8 | 3.7e-07 (`rho`), 2.2e-06 (`cool`) |
  | `wasp_he23off` | 13.22 (unchanged) | 1 - 1.2e-5 | 1.6e-03 (`rho`, `T`) at `r = 1.012` |

  Step counts move by 2 and 6 (13486 -> 13488, 13488 -> 13482). `wasp_he23off`
  is the more sensitive of the two because its largest differences sit in the
  narrow band at `r = 1.012-1.014` where the run already carries **negative
  species densities in both the golden and the new output** (H I, O I, O II,
  S II; 17 such cells in the `wasp_full` golden, 3 in the `wasp_he23off`
  golden). Those cells are excluded from the numbers above. The negative
  densities are PRE-EXISTING — they are in the goldens — and are unrelated to
  this change, but they are a physical-correctness defect in their own right
  and are recorded here. **Goldens were not re-snapshotted.**
  [2026-08-15: they have been since, most recently on 2026-08-15. All five
  golden `Ion_species.txt` files under `backup/regression/golden/` now carry
  zero negative densities.]

### 10.7 What is still not fixed

1. **No radiative floor.** §3.2(iii) stands unchanged: the model has no
   stellar optical/IR absorption and no thermal background, so nothing
   prevents a shielded layer from settling below `T_eq`. The equilibrium the
   base now heads for is set by where the guarded cooling meets the metal
   photoionization heating, not by `T_eq`.
   [2026-08-15: partly addressed. `Base IR field: True` (default off) gives the
   eight fine-structure lines and the H3+ bands the thermal field of the layer
   below the base, so those channels reach radiative equilibrium instead of
   radiating into vacuum. Stellar optical/IR absorption is still absent.]
2. **Ions with no explicit two-level term have no cooling below ~700 K.** The
   guard removes their coronal fit and nothing replaces it. For C I this omits
   the real `[C I] 609/370um` lines; their LTE rate at the HD 189733 b base is
   ~1e-10 erg cm^-3 s^-1, i.e. 1e-4 of the local heating, so it is negligible
   there but would not be in a colder or more carbon-rich base.
   [2026-08-15: closed. `[C I] 609um` and `[C I] 370um` are slots 1 and 2 of the
   eight-line set, and the split ground terms of C I, C II, N II and O I are
   solved in statistical equilibrium at the local `(ne, nHI)`
   (`docs/coronal_cutoff_width.md` section 7; `TO_BE_DONE.md` item (C)).]
3. ~~**Thick resonance lines in the wind are still treated as thin** (Mg II h&k
   in an ultrahot Jupiter). Nothing in this change moves that either way.~~
   **Settled 2026-08-11, no change needed.** The lines are thick
   (`tau0 ~ 7.6e4` for Mg II k at the WASP-121 b base) but their cooling is
   not suppressed by it: the correction to an optically thin coronal fit is
   `beta A_ul/(beta A_ul + ne q_ul)`, not `beta`, and with `A_ul ~ 1e8 s^-1`
   the effective critical density `beta A_ul/q_ul` sits one to six decades
   above the `ne` these winds reach. Measured effect on the total radiative
   losses: 0.02% (WASP-121 b), 0.002% (HD 209458 b), 0.0015% (HD 189733 b).
   `docs/resonance_line_trapping.md`.
4. **The checkerboard.** The mode tracks `H/dr` (§4-§5), and a warmer base
   raises `H`. Whether the warming is enough to leave the damped regime is not
   settled by these runs; see §10.4.

---

## 11. Base grid resolution as an input key (2026-08-11)

§7 recommendation 1 - refine the base grid - required editing
`define_grid.f90`, because `drc` and `N_low` were hardcoded locals. They are
now an input key:

```
Base grid [dr,cells]:  2.0e-4 50     # the default = the grid used everywhere above
Base grid [dr,cells]:  5.0e-5 200    # the 4x refinement of §4.3 G1
```

The two numbers are on one line because their product is the extent of the
uniform region (0.01 R_p), which the refinement is meant to hold fixed. The key
applies to `Grid type: Mixed` only. `EXHALE_setup.out` now also reports
`Base scale-height resolution: H(T_eq)/dr = 1/(b0*dr_j(1))` - the §4.1 quantity,
computed from the code's own Jeans parameter and the actual first cell size
after the grid smoothing - and warns below 10 cells. Its values are 26.9 cells
for HD 189733 b and 102.3 for WASP-121 b on the default grid, i.e. the same
ordering as the hand computation of §4.1 (21.6 / 83.0), the difference being
that §4.1 used a hand-built `mu` and the nominal `drc` rather than the smoothed
`dr_j(1)`.

### 11.1 The key reproduces the patched builds exactly

The §4.3 experiments were run from scratch builds with `drc` and `N_low`
edited in the source. Re-running the same configurations from the key, on the
current tree (which now also carries the metal-cooling fix of §10), reproduces
them digit for digit at 8000 cold-start steps - as it must, since the base is
still near `T_eq` there and the cooling fix has not had time to act:

| run | `dr` [R_p] | `T1` [K] | `H/dr` | `A(ln rho)` | `A(ln p)` | flips /11 |
|---|---|---|---|---|---|---|
| key `2.0e-4 50` (default) | 1.93e-04 | 1088.8 | 8.47 | 0.0046 | 0.00052 | 2 |
| §4.3 G0 (patched source) | 1.93e-04 | 1088.8 | 8.5 | 0.0046 | 0.00052 | 2 |
| key `1.0e-4 100` | 9.55e-05 | 1160.2 | 166.7 | 0.0002 | 0.00024 | 0 |
| key `5.0e-5 200` | 4.66e-05 | 1178.8 | 171.7 | 0.0001 | 0.00010 | 0 |
| §4.3 G1 (patched source) | 4.66e-05 | 1178.8 | 171.7 | 0.0001 | 0.00010 | 0 |

So the key is the same experiment, not a new one, and the 2x step (`1.0e-4
100`, not tried in §4.3) already removes the mode at this stage.

### 11.2 The cost of refining at fixed `N`

`N = 500` is a compile-time constant, so cells given to the base come out of
the stretched region.
[2026-08-15: no longer. `N` is a runtime variable set by the input key
`Grid cells: <N>`, defaulting to the 500 that used to be compiled in, so cells
given to the base need not come out of the stretched region.] Computed from `define_grid` for the HD 189733 b domain
(`r_max = 4.43 R_p`):

| `Base grid` | stretch ratio | `dr` at base | `dr` at 1.5 R_p | `dr` at `r_max` |
|---|---|---|---|---|
| `2.0e-4 50` | 1.01189 | 1.93e-04 | 6.06e-03 | 3.90e-02 |
| `1.0e-4 100` | 1.01587 | 9.55e-05 | 7.91e-03 | 5.11e-02 |
| `5.0e-5 200` | 1.02515 | 4.66e-05 | 1.26e-02 | 7.82e-02 |

A 4x finer base makes the wind grid 2x coarser at 1.5 R_p, and the CFL step
falls with the smallest cell, so the same physical time costs ~4x more steps.
Refining the base is therefore not free for `Mdot` or for the transmission
spectrum, both of which are formed in the region being coarsened. Raising `N`
alongside would need the static `(1-Ng:N+Ng)` arrays to become allocatable.
[2026-08-15: done. The grid-sized arrays are allocated by
`allocate_grid_arrays` once `N` is read, and `Grid cells:` raises `N`, so the
trade-off above applies only at fixed `N`.]

### 11.3 The alternating-amplitude metric mixes two different things

The `A(ln rho)` numbers of §2 and §4 are taken over cells 1-12, which includes
the ghost-to-cell-1 pressure/density jump of §2 (the 4.5x density inversion).
That jump is a lower-boundary artifact, not the 2*dr* mode, and it does not go
away with grid refinement. Splitting the window separates them; measured on the
production runs of §11.4 at the same stage of their transient, against the
pre-fix converged reference:

| run | `T1` [K] | cells 1-12 | cells 3-12 | cells 5-14 | flips (3-12) /9 |
|---|---|---|---|---|---|
| pre-fix converged reference | 236.1 | 0.0537 | **0.00903** | 0.00292 | 3 |
| `2.0e-4 50` (default) | 600.9 | 0.0248 | **0.00157** | 0.00006 | 1 |
| `1.0e-4 100` | 603.8 | 0.0171 | **0.00009** | 0.00001 | 0 |
| `5.0e-5 200` | 607.3 | 0.0148 | **0.00005** | 0.00000 | 2 |

Read off the cells 3-12 column, which is the mode proper:

* the metal-cooling fix of §10 alone (default grid, base at 607 K instead of
  236 K) drops the interior mode by a factor 6, from 9.0e-3 to 1.5e-3, but does
  **not** remove it - one velocity sign flip survives;
* the base refinement drops it a further factor 17-31, to the 1e-4 -- 1e-5
  level, which is the "mode removed" regime of §4.4;
* 2x (`1.0e-4 100`) is already at that level, so the useful step is smaller
  than the 4x of §4.3. The residual sign flips at 4x are in a velocity field
  whose interior amplitude is 3e-5 of the base value, i.e. they are not the
  mode.

The cells 1-12 column stays at 1.3-2.4e-2 in every case because it is measuring
the boundary jump. Nothing here addresses that jump; §7 remedy 2 (removing the
ghost/cell-1 temperature jump) is still the open item for it.

Caveat: the three runs above were mid-transient, at the same point in their own
relaxation but not at the same physical time as each other, and they are not
converged. The ordering is robust (it is the same ordering as §4.4 and §5.3)
but the individual values should not be quoted as converged amplitudes.

### 11.4 Production re-convergence: not completed at the time this was written

> **Superseded 2026-08-11 (later the same day): it was completed.** See §12.
> `HD189733b/input.inp` now carries `Base grid [dr,cells]: 1.0e-4 100` and the
> directory holds a Newton-converged solution. The state described below is the
> mid-transient one; the ceiling argument at the end of the section still
> stands and is why 2x, not 4x, is what production uses.


Cold-start HD 189733 b re-convergences at `2.0e-4 50`, `1.0e-4 100` and
`5.0e-5 200`, everything else as in `HD189733b/input.inp` (`Solver: Newton`,
`du_th 0.5 1e-3`, `CFL 0.3`, metals, He 2^3S, He/metal diffusion), were still
marching when this section was written and had not reached the JFNK hand-off
(`du < 1e-2`). State at 5.4e5 steps (default), 5.3e5 (2x), 6.1e5 (4x):
`du = 3.2e-2 / 2.2e-2 / 1.2`, `||R|| = 1.1e-2 / 1.0e-2 / 3.3e-1`, all three
still descending. For scale, the
pre-fix production run needed 3.99e5 steps before its hand-off and floored at
`du = 1.5e-2`. **No converged `Mdot` is quoted from these runs, and
`HD189733b/` was not touched.**

There is a hard ceiling worth recording. `count_max = 1000000` is a
compile-time parameter (`parameters.f90`; the `EXHALE_MAXSTEPS` environment
hook only lowers the cap, it does not raise it), and the CFL step falls with
the smallest cell.
[2026-08-15: the cap is an input key now, `Max steps: <N>`; `EXHALE_MAXSTEPS`
still only lowers it.] A 4x-refined base therefore reaches at most 1/4 of the
physical time the default grid gets within the same cap - the pre-fix run used
its full 1e6 steps. On this planet 4x refinement is thus not reachable to
convergence without raising `count_max`, whereas 2x is, and §11.3 says 2x is
where the mode already goes. That is the configuration to try first.
[2026-08-15: "not reachable without raising `count_max`" is no longer a
limitation -- `Max steps:` raises it from `input.inp`.]

---

## 12. Production re-convergence, and what the converged state still carries (2026-08-11)

The 2x configuration of §11.4 was carried to convergence and adopted.
`HD189733b/input.inp` now reads `Base grid [dr,cells]: 1.0e-4 100`; everything
else is as before (`Solver: Newton`, `du_th 0.5 1e-3`, `CFL 0.3`, metals,
He 2^3S, He/metal diffusion).

### 12.1 What the production run reached (from the run logs)

| stage | log | result |
|---|---|---|
| cold start, 2x base grid, marched to the `count_max` cap | `run_20260811_coldstart_x2.log` | `log10 Mdot = 9.24`, `du`-stopped |
| Newton finish restarted from it (`Load IC`, hand-off at step 2002) | `run_20260811_newton_finish.log` | JFNK `info = 0`, `\|\|R\|\| = 5.369e-04`, **`log10 Mdot = 9.04`** |

The 0.20 dex between the two is the familiar point that a `du`-threshold stop is
not quantitative for `Mdot`; 9.04 is the Newton-grade number for the
configuration of this section, i.e. with the legacy `T_eq`-pinned ghost.
(**Superseded 2026-08-11, later the same day.** The folder was re-converged with
`Base ghost temperature: continuous` and the ionization-root validation, and
gave 9.05; see §14. **Superseded again 2026-08-13**: after the ground-term
fine-structure statistical equilibrium, the base-ghost composition fix of §15 and
the H(n=2) rate corrections, the folder re-converges to `log10 Mdot = 9.14`
(`run_20260812_lyafix.log`, JFNK `info = 0`, `||R|| = 8.78e-04`), and that is
what `paper/ms.tex` and `python/paper_data.py` carry.)

### 12.2 The mode in the converged state (measured here)

Measured on `HD189733b/output/Hydro_ioniz.txt` with the §11.3 metric (projection
of the discrete second difference onto `(-1)^j`, divided by 4). The metric was
checked first against the §11.3 reference row and reproduces it digit for digit,
so the numbers below are on the same scale as §2, §4 and §11.3.

| run | `T` cell 1 [K] | `A(ln rho)` 1-12 | 3-12 | 5-14 | `A(ln p)` 1-12 |
|---|---|---|---|---|---|
| pre-fix converged reference (§2, §11.3) | 236.1 | 0.0537 | 9.03e-3 | 2.92e-3 | 1.97e-3 |
| **production, 2x grid, Newton-finished** | **493.6** | **0.0264** | **4.41e-3** | **1.62e-3** | **5.5e-4** |

Readings:

* The base temperature moved 236 K -> 494 K, as the §10 cooling fix predicts
  (§10.4 swept the code's own cooling to a balance point of 486 K; the run
  settles 8 K above that, which is within what the inflow term can account for).
* The pressure field is now smooth to 5.5e-4, a factor 3.6 better than the
  reference. The entropy/density alternation is what remains.
* **The interior mode is reduced by a factor 2, not by the factor 17-31 that
  §11.3 measured.** This is the part worth flagging. §11.3 compared runs
  mid-transient, at 8000 cold-start steps, where the 2x grid put the interior
  amplitude at 9e-5. The converged production state sits at 4.4e-3, i.e. ~50x
  above what the mid-transient comparison suggested, and only 2x below the
  pre-fix reference. Grid refinement alone therefore does **not** carry the
  planet to the "mode removed" regime once the run is taken all the way to a
  Newton-converged steady state.
* The cells 1-12 window stays at 2.6e-2 for the same reason as before: it is
  dominated by the ghost-to-cell-1 boundary jump (§11.3), which no amount of
  base refinement touches. §7 remedy 2 (removing the ghost/cell-1 temperature
  jump) remains the open item for it, and on the evidence of this section it is
  now the more promising of the two remedies, not the secondary one.

### 12.3 Attributed but not reproducible from the artifacts in this tree

Three further findings about this mode were reported alongside the production
re-convergence and are recorded here as attributions, because the runs that
produced them are not in the tree and the numbers could not be re-measured:

* that a single JFNK finish can *excite* the mode rather than damp it, taking
  the interior amplitude from 2.6e-3 to 7.0e-2;
* that a frozen-base anchor does not converge on this planet;
* that alternating marching and damping cycles removes about 4% of the residual
  amplitude per cycle, leaving a residual `A = 0.0138`.

None of the three could be checked against the stored outputs: the value 0.0138
does not correspond to any window or field of the converged production run under
the §11.3 metric (the candidates are listed in §12.2), and the intermediate
states of the excitation and damping experiments are not stored. The
*qualitative* claim that the Newton finish does not simply damp this mode is
independently supported by §12.2 -- the converged state carries a larger interior
amplitude than the mid-transient 2x run -- but the specific amplitudes above
should be re-measured before being quoted.

### 12.4 Scope

This section changes no source file. It records the configuration adopted for
production and one measurement on the resulting output. The open items of §10.7
are unchanged, and item 4 there ("whether the warming is enough to leave the
damped regime") is now answered: on the evidence of §12.2, not by itself.

---

## 13. Remedy 2 implemented and measured: `Base ghost temperature: continuous` (2026-08-11)

§7 remedy 2 (remove the ghost/cell-1 temperature jump) is now an input key,
`Base ghost temperature: isothermal | continuous`, default `isothermal` (the
legacy `T_0` pin, so every existing run is unchanged; the byte-identical
regression matrix passes).

### 13.1 What the key does

The lower ghost keeps the base composition -- `ntot_bc` nuclei and `dp_bc`
electrons at the pinned `rho_bc`, exactly the particle count the isothermal pin
uses -- but carries the temperature of the first interior cell:

```
p_ghost = (ntot_bc + dp_bc) * T_1 ,   T_1 = p_1 / (n_tot + n_e)_1 ,
```

in the code's `T0` units. `(n_tot + n_e)_1` is taken from the composition solve
(`get_species_densities`, the single policy point for what counts as a
particle), so the ghost pressure stays a differentiable function of the interior
pressure -- what the JFNK line search needs -- while the ionization state it
divides by is lagged exactly like every other composition quantity in a hydro
step. This differs from the §4.2 E2 scratch patch (`p_ghost = rho_bc*p_1/rho_1`),
which assumed the ghost carries the cell-1 composition per unit mass; the two
agree to 0.07% in `T_ghost` (the E2 row shows `T_ghost = 281.6` against
`T_1 = 281.8`), and the form above closes that gap by construction.

The key sets the same quantity as `Hydrostatic base: True` (the ghost pressure),
so the two are mutually exclusive; `Hydrostatic base` is tested first and
input_read warns when both are given. With `Base BC: pressure` the microbar
target is imposed at `T0` when `n0` is derived, so under a continuous-`T` ghost
only the base density remains anchored; input_read notes this too.

### 13.2 Measurement (production configuration, HD 189733 b)

Both runs are the production `HD189733b/input.inp` (`Load IC`,
`Base grid [dr,cells]: 1.0e-4 100`, `Solver: Newton 5.0e-2`, metals, He 2^3S,
He/metal diffusion), restarted from the converged production state in
`HD189733b/output/` and Newton-finished, in a scratch copy; the planet folder
itself was not written to. The baseline restart reproduces the stored production
state digit for digit under the §11.3 metric, so the comparison is like for like.

Both columns are the state as it stood **before** the ionization-root validation
of `docs/ionization_root_validation.md`; §14 gives the production numbers after
it.

| quantity | `isothermal` (production at the time) | `continuous` |
|---|---|---|
| `T_ghost` [K] | 1183.0 | 539.9 |
| `T` cell 1 [K] | 493.6 | 540.2 |
| `rho_1/rho_ghost` | 2.33 | 0.95 |
| `A(ln rho)` cells 1-12 | 0.0264 | **0.0092** |
| `A(ln rho)` cells 3-12 | 4.41e-3 | 9.28e-3 |
| `A(ln rho)` cells 5-14 | 1.62e-3 | 8.95e-4 |
| `A(ln p)` cells 1-12 | 5.46e-4 | 3.80e-4 |
| JFNK | `info = 0`, `\|\|R\|\| = 5.369e-04` | `info = 0`, `\|\|R\|\| = 2.757e-04` |
| `log10 Mdot` | 9.04 | 9.05 |

Window sweep of `A(ln rho)` over 10 cells, by first cell of the window:

| first cell | 1 | 3 | 5 | 7 | 9 | 11 | 13 |
|---|---|---|---|---|---|---|---|
| `isothermal` | 5.28e-2 | 7.74e-3 | 4.21e-4 | 1.21e-3 | 4.84e-4 | 6.62e-5 | 3.84e-6 |
| `continuous` | 1.25e-2 | 1.16e-2 | 4.04e-3 | 9.23e-5 | 1.14e-6 | 1.87e-6 | 1.75e-6 |

Readings:

* **The boundary jump is gone.** `T_ghost` and `T_1` now agree to 0.06% and the
  density inversion is removed (2.33 -> 0.95). The cells 1-12 amplitude, which
  §11.3 and §12.2 attributed to that jump, falls by 2.9x.
* **The mode is localized, not removed.** Beyond cell 7 the alternation falls
  13x (7-16 window) to 420x (9-18), i.e. the extended tail the isothermal pin
  carried out to cell ~20 is gone. What replaces it is a stronger disturbance
  confined to cells 2-4 (`T` 556, 592, 524 K), which is why the 3-12 window
  reads 2.1x *worse* while 5-14 reads 1.8x better. The velocity alternation
  decays monotonically (-259, +228, -120, +52, -16, +9, -1 cm/s) instead of
  sitting at +-9 cm/s out to cell 10.
* **The base velocity field changes character.** The isothermal ghost holds a
  -1125 cm/s inflow at cell 1 against a valved ~0 ghost; the continuous ghost
  gives +1389 cm/s outflow through the ghost and cell 1.
* **The steady residual improves 1.9x** and `Mdot` moves by +0.01 dex, i.e. the
  wind solution is unaffected within the du-stop path dependence.
* **Repeated restart cycles do not reduce the mode further.** Three further
  restart+Newton cycles of each run reproduce their own output digit for digit:
  both states are exact fixed points of the production workflow, so no residual
  decay is available from cycling. (This also settles the §12.3 attribution
  about damping cycles for this configuration: none of the ~4% decay per cycle
  reported there is reproducible here.)

### 13.3 Scope

Checked: the two restarts above and their Newton finishes; the byte-identical
regression matrix (`wasp_full`, `wasp_he23off`, `mol_base_handoff`) with the key
absent; a full WASP-121 b run with the key absent, which reproduces the HEAD
binary byte for byte.

Not checked: a cold start under the continuous ghost (the production cold start
costs 6 h and 1e6 steps, so this measurement says how the converged isothermal
state responds to the new closure, not what a run started from scratch under it
would settle on); any other planet; the transit observables; whether the
remaining cells 2-4 disturbance responds to the base grid.

---

## 14. Production state after adopting the key and the root validation (2026-08-11)

> **Dated 2026-08-13.** The state described in this section is the 2026-08-11
> one and is no longer what the folder holds. The current production run
> (`run_20260812_lyafix.log`) carries three later changes -- the ground-term
> fine-structure statistical equilibrium of `docs/coronal_cutoff_width.md` §7,
> the base-ghost composition fix of §15 below, and the H(n=2) rate corrections
> of `docs/lya_destruction_channels.md` §11 -- and converges to
> `log10 Mdot = 9.14` with `||R|| = 8.78e-04`, at a base temperature of 551 K
> rather than the 529 K tabulated here. The ghost/cell-1 agreement (0.05%) and
> the removed density inversion (`rho_1/rho_ghost = 0.95`) still hold in that
> run; the alternating amplitudes below were not re-measured on it.

`HD189733b/input.inp` now carries `Base ghost temperature: continuous` alongside
`Base grid [dr,cells]: 1.0e-4 100`, and the folder was re-converged with the
ionization-root validation of `docs/ionization_root_validation.md` in place
(`run_20260811_rootfix.log`, restarted from the stored IC, JFNK hand-off at step
2002): `log10 Mdot = 9.05`, JFNK `info = 0`, `||R|| = 1.724e-04`.

Measured on `HD189733b/output/Hydro_ioniz.txt` with the §11.3 metric, against
the two §13.2 columns (both of which predate the root validation):

| run | `T_ghost` [K] | `T` cell 1 [K] | `rho_1/rho_ghost` | `A(ln rho)` 1-12 | 3-12 | 5-14 | `A(ln p)` 1-12 |
|---|---|---|---|---|---|---|---|
| `isothermal`, pre-root-fix (§13.2) | 1183.0 | 493.6 | 2.33 | 0.0264 | 4.41e-3 | 1.62e-3 | 5.46e-4 |
| `continuous`, pre-root-fix (§13.2) | 539.9 | 540.2 | 0.95 | 0.0092 | 9.28e-3 | 8.95e-4 | 3.80e-4 |
| **production now** (`continuous` + root fix) | **529.4** | **529.6** | **0.95** | **0.0019** | **7.30e-4** | **5.58e-4** | **3.77e-5** |

Readings, stated tentatively because the two changes were not separated in this
state (the §13.2 A/B isolated the ghost closure; the root-validation gates
isolated the root fix; this row carries both):

* The boundary jump stays removed -- `T_ghost` and `T_1` agree to 0.04% and the
  density inversion is still gone.
* Every alternating amplitude is now lower than in either §13.2 column: 1-12 by
  14x against the isothermal pin and 4.8x against the pre-root-fix continuous
  run, and the cells 2-4 disturbance that §13.2 flagged as what replaces the
  extended tail has come down with it (`T` cells 1-8: 529.6, 539.6, 540.0,
  541.0, 540.0, 543.9, 542.3, 542.4 K, against 540.2, 556.1, 592.3, 523.8 in
  §13.2). The base velocity alternation is a few tens of cm/s (-58, +19, -5,
  +14, -2, +11, +1, +4) against +1389, -259, +228 before.
* The pressure field is smooth to 3.8e-5, an order of magnitude below anything
  measured earlier in this document.
* `A(ln rho)` over cells 3-12, the mode proper by §11.3, reads 7.3e-4 -- between
  the mid-transient 2x value (9e-5) and the converged pre-fix reference
  (9.0e-3), i.e. the interior mode is reduced but not removed.

The likely reading is that the root validation removed part of what §13.2
measured as the residual cells 2-4 disturbance: five of those cells carried
negative species densities in the stored state, and the run rejects exactly
those five at step 1 (see `docs/ionization_root_validation.md`). Separating the
two contributions would need a `continuous`-with-old-binary rerun, which was not
done.

Not checked here: a cold start under this configuration, the other planets, and
whether the remaining cells 1-4 alternation responds to further base refinement.

## 15. The cell 1-2 steady residual was a stale composition in the ghost, not the ghost closure (2026-08-12)

The measurement that opened this section asked which term of the finite-volume
steady residual carries the cells 1-2 spike that a `Base ghost temperature:
continuous` run of HD 189733 b reports (mass 53-60, energy 106-120 per sound
crossing time, unchanged under 8x base refinement and unchanged across initial
conditions, and absent from an `isothermal` run). The answer is that the spike
is not a property of the closure. It is produced by the residual evaluation
itself, which built the ghost from a particle count belonging to a different
state.

### 15.1 Term split of the residual (measured)

Instrumentation: a scratch copy of the tree with one added routine
(`base_face_terms` in `steady_residual.f90`) that repeats the production
`Reconstruct` / `Num_flux` / `source` path and writes, for each of the first
cells, the inner-face flux, the outer-face flux and the source separately, plus
the reconstructed states on both sides of every face. Configuration: the
production HD 189733 b input (`Base ghost temperature: continuous`, `Base grid
[dr,cells]: 1.0e-4 100`, `N = 500`, WENO3), evaluated on the converged state of
the grid-response series. `Fc` below is the wind mass-flux constant
`<rho v r^2>` over `[j_min:N]`, `5.9718e-08` in code units.

For mass, `S = 0` identically, so `R(1,j)` is nothing but the difference of the
two face flows. Measured face flows `r_edg^2 (rho v)_HLLC`, in units of `Fc`:

| face | 0 (base) | 1 | 2 | 3 | 10 | 20 |
|---|---|---|---|---|---|---|
| as evaluated | 2.19e5 | 1.05e5 | 128 | 115 | 101 | 88 |
| composition refreshed first | 142 | 265 | 128 | 115 | 101 | 88 |

The two innermost faces carried five orders more mass than the wind while faces
2 and outward sat at ~1e2 `Fc`. Cell by cell (code units, per sound crossing
time), with `out` and `in` the outer- and inner-face contributions:

| cell | `R1` | mass out | mass in | `R3` | energy out | energy in | `heat-cool` |
|---|---|---|---|---|---|---|---|
| 1 | -71.19 | 65.73 | 136.92 | -69.64 | 70.63 | 142.11 | 4.77e-03 |
| 2 | -65.64 | 0.080 | 65.72 | -69.93 | 0.087 | 70.61 | 4.47e-03 |
| 3 | -8.11e-03 | 0.0719 | 0.0800 | -1.21e-02 | 0.0777 | 0.0868 | 4.31e-03 |

Readings: the energy spike is not radiative -- `heat - cool` is four orders
below `R3` -- it is the same acoustic flux carrying enthalpy. The momentum row
behaves the same way: at cell 1 the WENO3 pressure term and gravity cancel to
`-4.07` out of `-230.3` and `-226.2`, and the `+58.2` flux difference is what is
left. Every component therefore points at one place, the flux through the base
face.

### 15.2 What the base face was seeing

`BC_component_constrho` closes the continuous ghost as
`p_ghost = (ntot_bc + dp_bc) * p(1) / n_part_cell1`, i.e. the base particle
count times the interior temperature `T(1) = p(1)/n_part_cell1`.
`n_part_cell1` is a module variable written by `get_species_densities`
(the single composition policy point) and read by `Apply_BC`. `input_read`
initializes it to the placeholder `ntot_bc + dp_bc`.

In the state measured, `ntot_bc + dp_bc = 1.00084` and the true cell-1 particle
count is `0.958583`, a ratio of `1.0441`. The standalone residual diagnostic
calls `Apply_BC` before any composition solve, so the cell-average ghost was
built with the placeholder, i.e. with ratio exactly 1:

* `W(3,0) = W(3,1)` to every printed digit -- a zero pressure gradient, not the
  continuous-temperature ghost the key asks for;
* that zero left-hand jump made ESWENO3 flatten cell 1, giving
  `WL(3,1) = 0.4990916` against `WR(3,1) = 0.4884203` from cell 2, a 2.2%
  pressure discontinuity *at face 1*;
* `Rec_BC`, which runs inside `assemble_residual` after `ioniz_eq` has refreshed
  `n_part_cell1`, rebuilt the ghost with the correct ratio,
  `WL(3,0) = 0.5210929` against `WR(3,0) = 0.4990927`, a 4.4% discontinuity *at
  face 0*.

So a single residual evaluation carried two mutually inconsistent ghosts, and
HLLC read both jumps as contact discontinuities and returned the corresponding
acoustic mass flux. This also explains the three properties that made the spike
look structural: the jump is a ratio of two particle counts, hence independent
of `dr` and of the initial condition, and `n_part_cell1` is read only by the
continuous branch, hence absent from `isothermal` runs.

### 15.3 Fix

The rule applied is that the ghost is built from the composition of the state it
bounds. Two call sites needed it:

* `init.f90` -- evaluate `get_species_densities` on the initial state before the
  first `Apply_BC`. This is the only entry point where `n_part_cell1` can still
  be the placeholder; every `Apply_BC` in the marching loop is already preceded
  by a composition solve.
* `steady_newton.f90`, `eval_residual` -- same call after `unpack_U` and before
  `Apply_BC`. Without it `F(Y)` also depended on the *previous* `Y` through
  `n_part_cell1`, so a finite-difference Jacobian column mixed two states.

The mass-flux base velocity carries the same defect and is seeded in the same
place (section 15.8).

Runs that do not set `Base ghost temperature: continuous` never read
`n_part_cell1`, so nothing outside that key can change; `make check` is
byte-identical (section 15.5).

With the composition refreshed, `W(3,0) = 0.5210898 = 1.0441 p(1)` as the key
intends, and `WL(3,1) = 0.4884517` now agrees with `WR(3,1) = 0.4884203` to
6.4e-05, i.e. 0.006%.

### 15.4 Effect on the residual (measured)

Local fractional residual rate `|R_k| / |u_k|` per sound crossing time, same
state, WENO3:

| grid | | cell 1 | cell 2 | cell 3 | wind `r > 2` |
|---|---|---|---|---|---|
| `dr = 9.55e-05` | mass, before | 5.97e+01 | 5.92e+01 | 7.67e-03 | 8.59e-04 |
| | mass, after | 6.41e-02 | 7.70e-02 | 7.67e-03 | 8.59e-04 |
| | energy, before | 1.12e+02 | 1.20e+02 | 2.18e-02 | 6.22e-03 |
| | energy, after | 1.31e-01 | 1.54e-01 | 2.18e-02 | 6.22e-03 |
| `dr = 2.40e-05` | mass, before | 5.59e+01 | 5.70e+01 | 1.29e-02 | 6.99e-04 |
| | mass, after | 9.98e-01 | 2.06e-01 | 1.29e-02 | 6.99e-04 |
| | energy, before | 1.10e+02 | 1.15e+02 | 2.40e-02 | 3.97e-03 |
| | energy, after | 2.00e+00 | 4.12e-01 | 2.40e-02 | 3.97e-03 |

Cells 3 and outward and the wind window are untouched to every digit, which is
the expected signature of a two-face boundary effect. The reduction at cell 1 is
930x on the default grid and 56x at 4x refinement; what is left grows under
refinement rather than staying flat, so it is not the same object.

An `isothermal` run of the same planet reproduces bit for bit across the fix,
and against it the corrected continuous closure is now the better one at every
one of the first five cells (mass 6.4e-02, 7.7e-02, 7.7e-03, 2.1e-03, 6.8e-04
against 8.2e-02, 1.2e-01, 3.4e-02, 2.4e-02, 1.9e-02) -- which restores, rather
than contradicts, the section 13-14 reading.

### 15.5 Gates

* Re-convergence, HD 189733 b, restart from the stored state, 6000-step cap, JFNK
  finish. The marching phase is unchanged to four digits (`step 2000`: flux
  spread `3.601e-03`, `||R||(ref) 1.941e-03` in all three builds). The JFNK
  scaled merit `||Fs||_2`, which unlike `||R||` includes the base cells, reads
  `1.43` before, `5.02` with the `init` fix alone (the residual then re-imported
  its own lag at every call) and `0.315` with both -- 4.5x below the starting
  point. `info = 0` throughout; `||R||` `7.13e-04` / `9.90e-04` / `8.48e-04`;
  `log10 Mdot` `9.15` / `9.14` / `9.14`.
* Residual of each build on its own converged state: cells 1-2 mass
  `6.06e+01, 6.15e+01` before against `1.43e-01, 2.58e-02` after (425x, 2380x),
  energy `1.13e+02, 1.21e+02` against `2.78e-01, 5.07e-02`.
* `make check` byte-identical, the fix being unreachable for the golden cases.
  WASP-121 b is covered by `wasp_full`, and it does not set the key.

### 15.6 What this does not fix

After the correction the first faces still carry ~1.4e+2 to 2.7e+2 times the
wind mass flux, and the excess decays outward over `r - 1 ~ 0.05` (roughly 500
cells on this grid), so it is the broad near-base non-stationarity documented in
sections 12 and 14, not a boundary closure defect. The factor ~2 by which face 1
exceeds its neighbors corresponds to the sign alternation of the velocity
between cells 1 and 2 (-21 and +21 cm/s in the state measured), the checkerboard
remnant of section 14. The cell-1 residual is still ~1.7e+2 times the wind level
and grows under base refinement.

On this evidence the alternative ghost closures that were drafted before the
measurement -- anchoring `rho_bc` on the base face instead of the cell center,
or deriving the ghost velocity from the mass-flux constant -- were aimed at a
disturbance that was not there, and none was implemented. The reading is
tentative in one respect: what remains has been localized (faces, and the
cell 1-2 velocity alternation) but not attributed to a specific term.

### 15.7 A bounded inconsistency left in place

The in-loop residual monitor calls `Apply_BC` with the previous step's
`n_part_cell1` and then `Reconstruct` with the current one, so it carries the
same mismatch, of the size of one step's change in the cell-1 particle count
(~1e-08 relative near convergence, against the 4.4% above). Inside a marching
step the three RK stages are self-consistent -- their `Apply_BC` and their
`Rec_BC` read the same value, lagged by one step like every other composition
quantity in an operator-split step -- while the `Apply_BC` calls that follow the
composition solve read the value that solve just wrote, so the lag is not uniform
across the step either. Both are bounded by one step's change and both are left
as they are; changing them would alter continuous-key results for no measurable
gain.

### 15.8 The same defect in the mass-flux base velocity

`Base velocity: massflux` sets the ghost velocity to `F_c/(rho_bc r^2)`, with
`F_c` the wind mass-flux constant that the marching loop refreshes each step.
`F_c` starts at `-1`, and `Apply_BC` guards on `base_flux_const > 0`, so before
the first refresh the key silently falls back to the valve. Measured: with the
key appended to the HD 189733 b input, the residual diagnostic of the pre-fix
build reproduces the valve run bit for bit, i.e. the key had no effect on
anything that stops right after `init`. `init` now seeds `F_c` from the initial
state by the same average over `[j_min:N]` the loop takes, and the key does
change the residual (cells 1-2 mass `6.377e-02, 7.699e-02` against
`6.412e-02, 7.701e-02` for the valve). The guard stays, so a state with no
outward flux still falls back as before.

### 15.9 Reproduction

Scratch tree, diagnostic routine and run directories:
`scratchpad/base_close/` (session scratchpad, no longer present) --
`t_cont` / `t_cont_ord` / `t_fix` (default grid,
before / order forced / fixed), `t_x4` / `t_x4o` (`N = 916`), `t_iso` /
`t_iso_fix`, `march_pre` / `march_post` / `march_post2` (re-convergence),
`rz_pre` / `rz_post2` (residual of each converged state). The added routine
writes `output/base_terms.txt`; it is a scratch instrument and was not added to
the tree.

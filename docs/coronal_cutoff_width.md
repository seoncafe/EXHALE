# The coronal-fit cutoff width `w`: sensitivity and a physical bound

**2026-08-12.** `coronal_excitation_cutoff(T)` in
`src/modules/radiation/Cool_coeff.f90` multiplies every CHIANTI-derived metal
line-cooling coefficient by

```
cutoff(T) = 1                             T >= T_floor = 1e3 K
          = exp(-x^2),  x = (T_floor/T - 1)/w    T <  T_floor
```

`w` is exposed as `Coronal cutoff width: <w>` and was shipped at 0.5. It is a
modeling choice, not a measured quantity. This note measures what it controls,
proposes two criteria that can be evaluated numerically, and reports the value
they select.

---

## 0. Verdict

**`w = 0.5` is not supported by the physics; the measured justified window is
`0.08 <= w <= 0.13`, and the default has been changed to `w = 0.1`.** The
argument, in order:

> [2026-08-15: read section 7.2 before using the window. The ground-term
> statistical equilibrium of section 7 removed the coronal residual the window
> was bounding, and the base balance temperature is now identical to every
> printed digit over `w = 0.02-1.20`. The window is a record of how the default
> was chosen, not a live constraint; `w = 0.1` is kept because a fit outside its
> range still needs a guard.]

1. What the guard removes below `1e3` K is dominated by the **ground-term
   fine-structure floor of C I**, which the CHIANTI fits carry in the
   optically thin *coronal* (low-density) limit. At the base densities of
   these planets those lines are far above their critical densities, so the
   coronal value exceeds the LTE-saturated rate the same lines can actually
   radiate by roughly seven decades. Removing it faster therefore moves the
   cooling function **toward** the physically correct value, not away from it,
   over the whole range below the floor.
2. At the temperature the HD 189733 b base actually occupies (~500-600 K),
   `w <= 0.3` already reproduces the density-saturated cooling to the digits
   printed, while `w = 0.5` leaves a spurious residual ~60x larger than the
   true metal cooling there.
3. The counter-consideration is numerical: a narrow guard is steep. The
   measured steepest slope `max |dln(Lambda)/dln T|` grows as roughly
   `5/w + 13`, against the `8.3-8.7` the cooling curve already has inside the
   fitted range. Keeping the guard within an order of magnitude of that
   intrinsic steepness gives `w >= 0.08`.
4. Confining the band in which the unsaturated residual still dominates to
   within 30% of the fit floor gives `w <= 0.13`.

**What does not change with `w`.** Mass-loss rate and transmission spectra are
insensitive. Over `w = 0.05-3.0` (a factor 60) on HD 189733 b the mass-loss
rate moves by 1.7% and the transit line depths by <= 1.8%; only the innermost
`r < 1.02 R_p` moves appreciably.

**What no `w` fixes.** The guard is 1 at `T_floor` by construction (it is
continuous in value and in `dT`-slope there), so just below the floor the
unsaturated residual is `1.1e4-5.4e4` times the explicitly saturated
fine-structure cooling for *every* `w`. Above the floor the guard is inactive
by design, and the same unsaturated coronal floors are fully in play. The
substantive fix is to give C I, the `[O I] 145.5 um` channel and N II the same
critical-density saturation that `[O I] 63 um` and `[C II] 158 um` already
have in `cool_OI_ne_func` / `cool_CII_ne_func`. After that the base result
should become nearly `w`-independent. `w` is a stopgap, and this note only
bounds it.
[2026-08-15: that fix is implemented, the ground terms of C I, C II, N II and
O I are solved in statistical equilibrium at the local `(ne, nHI)`
(section 7; `TO_BE_DONE.md` item (C), RESOLVED 2026-08-12). The base result is
`w`-independent to every printed digit (section 7.2).]

---

## 1. How the numbers were obtained

**Local-equilibrium sweep.** A driver was linked against the compiled
`Cooling_Coefficients` object and assembles the total radiative cooling of one
cell exactly as `T_equation` / `eval_cool` do in the default branch
(`cno_chianti = .true.`, `use_2lev_cool = .false.`):

```
cool = ne*(brem + coex + reco + coio) + ne*sum_i n_i c_i(T)
```

with the `[O I] 63 um` / `[C II] 158 um` coefficients replaced by their
density-dependent two-level forms and Fe II by the 2-D statistical-equilibrium
table. Frozen while `T` is swept: the ionization state, `ne`, `nHI`, the escape
probabilities `beta` (a column, so not a function of this cell's trial
temperature) and the heating rate, all read from the converged run.

`T_eq,base(w)` below is the root of `cool(T) = heat` with everything else held
at the converged state. It is the balance point of the cooling change alone,
not a prediction of the converged base temperature, which also exchanges energy
with its neighbors.

The base cell is `j = 1`. Note that the output files carry the ghost cells
(`write_output` loops `j = 1-Ng .. N+Ng` with `Ng = 2`), so the physical base
cell is the **third** row, not the first. The first two rows are the inner
boundary condition; for HD 209458 b, WASP-52 b and WASP-121 b that ghost sits
at `T_eq` while the physical base cell sits far below it.

**Re-convergence.** HD 189733 b was re-run at `w = 0.05, 0.1, 0.2, 0.5, 1.2,
3.0` from the converged reference as initial condition (`Load IC? True`),
production `input.inp` otherwise unchanged, followed by
`EXHALE_transit.py`. Every run reached `info = 0` in the JFNK finish in the
same 2002 marching steps.

Base cell states (converged runs as of this date):

| planet | `r/Rp` | `T` [K] | `ne` [cm^-3] | `nHI` [cm^-3] | heat [erg cm^-3 s^-1] | `beta([O I] 63um)` |
|---|---|---|---|---|---|---|
| HD 189733 b | 1.000096 | 528.5 | 8.41e8 | 8.78e13 | 3.40e-6 | 0.732 |
| HD 209458 b | 1.000193 | 518.8 | 1.49e8 | 2.51e14 | 1.93e-7 | 0.103 |
| WASP-52 b | 1.000192 | 794.2 | 4.73e8 | 8.85e12 | 1.80e-6 | 0.737 |
| WASP-121 b | 1.000196 | 2406.9 | 1.88e9 | 2.80e12 | 8.59e-6 | 0.925 |

Only HD 189733 b lets its base float (`Base ghost temperature: continuous`);
the other three pin the inner ghost at `T_eq`. Three of the four nevertheless
have their physical base cell below the `1e3` K floor, so the guard is active
there. WASP-121 b is the exception.

---

## 2. What the guard is actually removing

Volumetric cooling at `T = 999` K, one line per line-cooling ion, split into
the surviving coronal part and the explicit two-level fine-structure part
(erg cm^-3 s^-1; ions below 1e-12 omitted):

| ion | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| C I (coronal) | 6.57e-4 | 3.35e-4 | 3.54e-5 | 3.21e-6 |
| O I (coronal remainder) | 5.82e-5 | 2.95e-5 | 3.31e-6 | 4.15e-6 |
| Fe II | 3.97e-8 | 4.40e-9 | 1.33e-8 | 5.99e-9 |
| C II (coronal remainder) | 2.44e-8 | 5.3e-16 | 1.56e-8 | 3.32e-7 |
| **coronal total** | **7.15e-4** | **3.64e-4** | **3.88e-5** | **7.71e-6** |
| **explicit two-level FS** | **1.68e-8** | **6.79e-9** | **1.71e-9** | **6.83e-10** |
| ratio | 4.25e4 | 5.36e4 | 2.26e4 | 1.13e4 |

C I carries 82-92% of the residual in the three cool-base runs. Its coefficient
there is dominated by the softest exponential of the fit,
`1.106e-20 exp(-2351.38/T)/sqrt(T)`, whose 2351 K matches no C I level: the
ground-term splittings are 23.6 K and 62.4 K and the next term (`1D2`) is at
1.47e4 K. The exponential is a shape parameter of the multi-exponential
representation, and what it is representing at the low end is the
`[C I] 609/370 um` ground-term fine-structure emission.

**That emission is not coronal at these densities.** Taking the standard
`[C I]` transition probabilities (`A(609 um) = 7.93e-8 s^-1`,
`A(370 um) = 2.65e-7 s^-1`) and the ground-term Boltzmann populations, the LTE
ceiling on C I cooling at the HD 189733 b base (`n(C I) = 2.35e10 cm^-3`) is

```
~2.0e-11 erg cm^-3 s^-1   (essentially flat over 500-1000 K)
```

i.e. `6e-6` of the local heating, and about `3e7` times *below* the coronal
value at 999 K. The critical densities behind that are of order 1e0-1e3 cm^-3
for these lines, against `ne = 8.4e8` and `nHI = 8.8e13` at the base. (This
estimate was computed outside the code, from standard atomic data, and is
quoted to an order of magnitude.)

So below the fit floor the correct metal cooling is essentially the two
explicitly saturated lines plus a negligible saturated C I term, and the
residual the guard removes is an artifact of applying a coronal curve at
`n >> n_crit`. **Removing it faster is more correct, not less.**

Measured, for the HD 189733 b base cell: total cooling `cool(T)` in
erg cm^-3 s^-1, against the local heating `3.40e-6`. The `w = 0.02` column is
the physically motivated reference (guard fully closed; it equals the
saturated `[O I] 63 um` + `[C II] 158 um` cooling plus the H/He terms).

| `T` [K] | reference | `w=0.1` | `w=0.2` | `w=0.3` | `w=0.5` | `w=1.2` |
|---|---|---|---|---|---|---|
| 302 | 3.33e-8 | 3.33e-8 | 3.33e-8 | 3.33e-8 | 3.33e-8 | 4.74e-7 |
| 502 | 4.12e-8 | 4.12e-8 | 4.12e-8 | 4.35e-8 | **2.46e-6** | 6.23e-5 |
| 694 | 4.67e-8 | 4.67e-8 | 2.52e-6 | 3.73e-5 | 1.49e-4 | 2.85e-4 |
| 796 | 4.92e-8 | 7.13e-7 | 8.91e-5 | 2.21e-4 | 3.51e-4 | 4.36e-4 |
| 894 | 5.14e-8 | 1.42e-4 | 4.09e-4 | 4.97e-4 | 5.50e-4 | 5.76e-4 |
| 1003 | 7.21e-4 | 7.21e-4 | 7.21e-4 | 7.21e-4 | 7.21e-4 | 7.21e-4 |

The bold entry is the one that matters: at the temperature this base actually
sits at, `w = 0.5` is 60x above the physically motivated value while
`w <= 0.3` is on it. The 1003 K row is above the floor, where the guard is 1
for every `w` including the reference column -- that row is the unsaturated
coronal value, and it is the residual defect item (1) of section 6 addresses.

---

## 3. Criteria and the window they select

### (i) Confine the band where the unsaturated residual still dominates

Below the floor the guard has to bring the residual under the explicitly
saturated fine-structure term. Writing `R` for the ratio of the two at the
floor (table in section 2, `1.1e4-5.4e4`), the crossing sits at
`x = sqrt(ln R) = 3.0-3.3`, i.e.

```
T_cross = T_floor / (1 + w*sqrt(ln R))  ~  T_floor / (1 + 3.2 w)
```

which reproduces the measured values to within 1% for `w <= 0.1`, 2% at
`w = 0.2` and 7% at `w = 0.5` (`R` drifts slightly with `T`). Measured
`T_cross` [K]:

| `w` | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| 0.05 | 861 | 860 | 865 | 868 |
| 0.10 | 759 | 756 | 764 | 768 |
| 0.20 | 617 | 615 | 625 | 627 |
| 0.30 | 525 | 522 | 533 | 532 |
| 0.50 | 410 | 407 | 418 | 411 |
| 1.20 | 247 | 244 | 254 | 236 |
| 3.00 | 139 | 138 | 140 | 116 |

Requiring `T_cross >= 0.7 T_floor` gives `w <= 0.13` (the binding planet is
HD 209458 b, with the largest `R`). `w = 0.5` puts the crossing at
`0.41 T_floor`: the artifact is the dominant metal coolant over a factor 2.4 in
temperature, and the HD 189733 b base sits inside that band.

### (ii) Do not make the guard the stiffest feature of the energy equation

`max |dln(Lambda)/dln T|` over 50-1000 K, and (last row) the steepest slope the
same cooling function already has *inside* the fitted range 1e3-1e5 K, where
the guard is exactly 1 and plays no part:

| `w` | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| 0.02 | 258 | 281 | 227 | 100 |
| 0.05 | 111 | 121 | 97 | 42 |
| 0.10 | 62.0 | 68.0 | 53.9 | 22.3 |
| 0.20 | 37.3 | 41.2 | 32.1 | 12.6 |
| 0.30 | 28.9 | 32.0 | 24.7 | 9.3 |
| 0.50 | 21.7 | 24.1 | 18.4 | 6.7 |
| 1.20 | 13.8 | 15.5 | 11.5 | 4.1 |
| 3.00 | 9.4 | 11.1 | 7.2 | 2.8 |
| *intrinsic (1e3-1e5 K)* | *8.3* | *8.3* | *8.3* | *8.7* |

The intrinsic maximum sits near 1.2e4 K and is set by the Lyman-alpha
collisional-excitation term, i.e. by physics the solvers already handle.
Requiring the guard to stay within an order of magnitude of it gives
`w >= 0.08`.

This criterion is numerical, not physical, and in these runs it was **not
observed to bind**: HD 189733 b converged with `info = 0` in an identical 2002
steps at every `w` tested down to 0.05 (slope 111). The final JFNK residual is
the only visible trace, and it is not monotone in `w`
(`9.6e-4`, `9.7e-4`, `8.7e-4`, `1.2e-4`, `2.2e-4`, `8.2e-4` for
`w = 0.05, 0.1, 0.2, 0.5, 1.2, 3.0`). The criterion is kept as a guard rail
rather than as evidence of a problem.

### Window

`0.08 <= w <= 0.13`. `w = 0.1` is the round value inside it and is the new
default.
[2026-08-15: superseded as a live constraint by section 7.2; `w = 0.1` remains
the default.]

---

## 4. Sensitivity map

### 4.1 Local balance temperature `T_eq,base(w)` [K]

| `w` | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| 0.02 | 956 | 948 | 966 | 1049 |
| 0.05 | 898 | 881 | 921 | 1049 |
| 0.10 | 818 | 789 | 857 | 1049 |
| 0.20 | 700 | 660 | 759 | 1049 |
| 0.30 | 618 | 572 | 688 | 1049 |
| 0.50 | 511 | 460 | 594 | 1049 |
| 0.80 | 421 | 367 | 513 | 1049 |
| 1.20 | 355 | 299 | 456 | 1049 |
| 2.00 | 290 | 231 | 404 | 1049 |
| 3.00 | 253 | 191 | 379 | 1049 |

WASP-121 b is flat because its balance point (1049 K) lies above the fit floor,
where the guard is inactive; its converged base (2407 K) is higher still and is
held by the inner boundary, so the guard never touches it. The other three all
have balance points below the floor and move by 500-700 K across this range of
`w`. The `w -> 0` limit is not unbounded: it saturates just under 1000 K,
because at the floor the guard is 1 and the C I residual is already ~200x the
local heating. **The whole family of cutoffs can only move the base within
roughly 200-1000 K, and the upper end is set by the position of the fit floor,
not by `w`.**

### 4.2 Re-converged HD 189733 b

Warm restart from the converged reference, full production settings, Newton
finish in every case.

| `w` | `T` base [K] | `T(r=1.01)` [K] | Mdot at `r=2` [g/s] | JFNK `\|\|R\|\|` | steps |
|---|---|---|---|---|---|
| 0.05 | 579.9 | 969 | 2.2823e9 | 9.56e-4 | 2002 |
| 0.10 | 567.3 | 1145 | 2.2598e9 | 9.70e-4 | 2002 |
| 0.20 | 612.4 | 1127 | 2.2583e9 | 8.74e-4 | 2002 |
| 0.50 | 517.7 | 1845 | 2.2490e9 | 1.15e-4 | 2002 |
| 1.20 | 500.1 | 2062 | 2.2464e9 | 2.19e-4 | 2002 |
| 3.00 | 486.8 | 2093 | 2.2434e9 | 8.17e-4 | 2002 |

Mdot spans 1.7% over a factor 60 in `w`. The converged base temperature is not
monotone in `w` (the base is not in local radiative balance; it exchanges
energy with its neighbors and with the boundary), but it stays in a ~110 K
band, consistent with the "~100 K level" already recorded in the code.

Radial extent of the difference, relative to `w = 0.5`:

| quantity | `r < 1.1` | `r >= 1.1` |
|---|---|---|
| `\|dT/T\|`, `w=0.2` | 0.39 | 6.1e-4 |
| `\|drho/rho\|`, `w=0.2` | 1.05 | 1.1e-2 |
| `\|dT/T\|`, `w=1.2` | 0.17 | 2.8e-4 |
| `\|drho/rho\|`, `w=1.2` | 0.18 | 2.2e-3 |

Everything outside `r ~ 1.05` agrees to a few percent or better; by `r = 1.1`
the profiles are identical to ~1e-3.

### 4.3 Transmission spectra (HD 189733 b)

Line depth `1 - min(T_rot+instr)`, relative to `w = 0.5`:

| line | depth at `w=0.5` | `w=0.05` | `w=0.1` | `w=0.2` | `w=1.2` | `w=3.0` |
|---|---|---|---|---|---|---|
| He I 10830 | 1.949e-2 | +0.95% | +0.38% | +0.34% | -0.06% | -0.14% |
| Ly-alpha | 2.762e-1 | +0.37% | +0.10% | +0.10% | -0.07% | -0.08% |
| H-alpha | 5.168e-3 | +1.16% | +0.51% | +0.46% | -0.10% | -0.16% |
| Mg II | 3.854e-3 | +0.60% | +0.34% | +0.30% | -0.09% | -0.13% |
| Ca II | 2.374e-3 | +0.35% | +0.20% | +0.21% | -0.08% | -0.12% |
| Na I D | 2.950e-3 | +1.62% | +1.03% | +0.97% | -0.33% | -0.48% |

Equivalent widths move more, because they pick up the wings that sample the
inner layer: up to +4.8% (Mg II), +4.6% (Na I D), +3.2% (Ly-alpha) at
`w = 0.05` relative to `w = 0.5`.

**Verdict: the observables are robust.** The choice of `w` is a statement about
the base structure, not about anything this code compares to data.

---

## 5. Change made

- `src/modules/init/parameters.f90`: `coronal_cutoff_width` default
  `0.5d0 -> 0.1d0`.
- Comment blocks in `Cool_coeff.f90` (`coronal_excitation_cutoff`),
  `input_read.f90`, `docs/input_schema.md` and
  `docs/EXHALE_user_manual.tex` updated to the new default and the window.
  [2026-08-15: `README.md` was listed here as well but never carried the key;
  it is not part of this change.]

Nothing else changed. Any run can restore the old behavior with
`Coronal cutoff width: 0.5` in `input.inp`.

Golden regression: `make check` after the change is **PASS, all three cases
byte-identical** (`wasp_full` 13483 steps, `wasp_he23off` 13478,
`mol_base_handoff` 12000). That is not a validation of the change: the two
WASP-121 b cases never put a cell below `1e3` K (converged minima 2163 K and
2157 K), and `mol_base_handoff`, which does go down to 860 K, runs metals-off
-- the guard multiplies only the CHIANTI metal coefficients. **The regression
matrix does not exercise this guard.** The goldens were not re-snapshotted
(there was nothing to re-snapshot).

---

## 6. What this left open

1. **The saturation, not the cutoff, is the physics.** C I, the
   `[O I] 145.5 um` channel and N II need the same critical-density treatment
   that `[O I] 63 um` and `[C II] 158 um` have. Until then the coronal floors
   remain overestimated *above* `1e3` K as well, where no cutoff acts. The
   three cool-base runs here sit in that region.
   **Done, 2026-08-12; section 7.**
2. **The guard is temperature-keyed, not density-keyed.** In genuinely
   low-density cold gas the coronal limit is the valid one and the guard would
   remove real cooling. In the runs examined this is moot: the only
   low-density cells below `1e3` K are the outermost four
   (`rho ~ 1.3e5 cm^-3` at `r = 4.3` in HD 189733 b) and they carry 0.00% of
   the integrated cooling. A density-aware saturation would remove the concern
   along with item 1. **Removed with item 1: the fine-structure lines that
   dominated below the floor are now solved in statistical equilibrium at the
   local density, and what the guard still switches off is negligible at any
   density (section 7.1).**
3. **The base has no radiative floor**. With the guard narrowed, the
   HD 189733 b base at ~570 K has a
   metal cooling of `~5e-8` against a heating of `3.4e-6`: it is not in local
   radiative balance at all, and its temperature is set by the boundary and by
   the flow. A wide `w` hides this behind an artificial coolant rather than
   fixing it. **Still open, and now visible without the artifact: section
   7.4.**

---

## 7. Follow-up (2026-08-12): the fine-structure floors are saturated, and
##    the base is no longer sensitive to `w`

Item 1 above is implemented. Every C/N/O coolant with a split ground term --
C I (`2p2 3P`), C II (`2p 2P`), N II (`2p2 3P`) and O I (`2p4 3P`) -- now has
its ground term solved in exact statistical equilibrium at the local
`(ne, nHI)`, with electron and H-atom collisions and with the escape
probability of each line entering as `A_ul -> beta_ul A_ul` inside the
solution. The coronal curve is kept only for the channels that leave the
ground term. Physics, atomic data and provenance:
`docs/cooling_formulas.tex` section "Ground-term fine-structure statistical
equilibrium"; generator `cooling_data/fit_fs_saturation.py`.

This set is complete: N I and O II have a single-level `4S` ground term,
Mg I/II, Ca II and Na I have a single ground level, Fe I is built from
permitted lines only (its `a5D` fine structure is not in the table at all),
and Fe II is already density-dependent through its 2-D statistical-equilibrium
table.

### 7.1 How large the overestimate was

Measured at the converged base cell of each run: the coronal within-term
channel as the CHIANTI fits carry it, against the exact statistical-equilibrium
value at the same `(T, ne, nHI)`, both volumetric [erg cm^-3 s^-1] with
`beta = 1`. The last column is the saturated value as a fraction of the local
heating rate.

| run | ion | coronal | SE (saturated) | decades | SE/heat |
|---|---|---|---|---|---|
| HD 189733 b, `T = 528` K | C I | 1.61e-4 | 2.00e-11 | 6.9 | 5.9e-6 |
| | C II | 1.14e-4 | 2.09e-12 | 7.7 | 6.2e-7 |
| | N II | 7.47e-10 | 5.97e-17 | 7.1 | 1.8e-11 |
| | O I | 1.47e-4 | 3.08e-8 | 3.7 | 9.1e-3 |
| HD 209458 b, `T = 519` K | C I | 7.77e-5 | 5.76e-11 | 6.1 | 3.0e-4 |
| | C II | 2.50e-12 | 2.56e-19 | 7.0 | 1.3e-12 |
| | N II | 1.28e-15 | 5.72e-22 | 6.3 | 3.0e-15 |
| | O I | 7.42e-5 | 8.77e-8 | 2.9 | 4.5e-1 |
| WASP-52 b, `T = 794` K | C I | 2.15e-5 | 1.94e-12 | 7.0 | 1.1e-6 |
| | C II | 6.14e-5 | 2.43e-12 | 7.4 | 1.4e-6 |
| | N II | 8.42e-10 | 1.48e-16 | 6.8 | 8.2e-11 |
| | O I | 8.92e-6 | 3.40e-9 | 3.4 | 1.9e-3 |
| WASP-121 b, `T = 2407` K | C I | 8.08e-6 | 4.48e-14 | 8.3 | 5.2e-9 |
| | C II | 8.19e-4 | 1.33e-11 | 7.8 | 1.6e-6 |
| | N II | 1.14e-8 | 8.49e-16 | 7.1 | 9.9e-11 |
| | O I | 1.26e-5 | 1.20e-9 | 4.0 | 1.4e-4 |

The `coronal` column is what the fits carry; it is not in every case what the
code was using, because C II 158 um and O I 63 um already had a two-level
saturation. Section 7.3 separates the two.

### 7.2 The `w` sensitivity is gone

Same construction as section 4.1 (local balance of the base cell, everything
but `T` frozen at the converged profile), run with the code before and after
the change. `T_eq,base` [K]:

| `w` | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| **before** 0.02 | 956.04 | 948.10 | 966.40 | 1046.94 |
| 0.10 | 817.67 | 789.46 | 856.90 | 1046.94 |
| 0.20 | 699.79 | 659.54 | 759.00 | 1046.94 |
| 0.30 | 617.63 | 571.57 | 688.32 | 1046.94 |
| 0.50 | 511.05 | 460.09 | 594.10 | 1046.94 |
| 1.20 | 354.77 | 298.81 | 455.61 | 1046.94 |
| spread | **601 K** | **649 K** | **511 K** | 0 |
| **after** 0.02-1.20 | 1304.69 | 1101.43 | 1607.05 | 2443.45 |
| spread | **0.00 K** | **0.00 K** | **0.00 K** | **0.00 K** |

The "after" rows are identical to every digit printed across a factor 60 in
`w`. The reason is structural rather than numerical: with the fine-structure
floors saturated, the softest exponential surviving in any of the coronal
remainders is `exp(-16400/T)` (C I), which at 500 K is `1e-15` of its 1e4 K
value. There is no longer any coronal cooling below `1e3` K for the guard to
remove, so the balance point moves above the fit floor -- where the guard is 1
by construction, for every `w`. The window `0.08 <= w <= 0.13` derived in
section 3 was a bound on an artifact that no longer exists; `w = 0.1` is kept
because the guard is still the right thing to do with a fit outside its range,
but nothing in these runs now depends on its value.

A side result: the WASP-121 b base cell balances at 2443 K against a converged
base of 2407 K, i.e. it is now close to local radiative balance. Before the
change its balance point was 1047 K, 1360 K below where the boundary holds it.

### 7.3 Effect on the converged solutions

> **Dated 2026-08-13.** The runs below isolate this change and are not the
> production set. The planet folders were re-converged again after the base-ghost
> composition fix and the H(n=2) rate corrections, and now give
> `log10 Mdot = 9.46` (HD 209458 b), `9.14` (HD 189733 b), `11.70` (WASP-52 b)
> and `13.20` (WASP-121 b); see each folder's `run_20260812_lyafix.log`.

**HD 189733 b**, full re-convergence at `w = 0.1` from the same warm restart,
production `input.inp` unchanged, JFNK finish reached `info = 0` in the same
2002 marching steps as the reference.

| quantity | before | after | change |
|---|---|---|---|
| base `T` [K] | 528.5 | 613.2 | +84.7 K |
| `T(r = 1.01)` [K] | 1927.7 | 1952.9 | +1.3% |
| max `T` [K] | 12236 | 12115 | -1.0% |
| `Mdot` at `r = 2` [g/s] | 2.3509e9 | 2.8498e9 | **+21.2%** |

What moved is the base layer. At `r = 1.0005` the C I coronal floor had been
`3.3e-6` erg cm^-3 s^-1 against a local heating of `3.4e-6` -- it was the
dominant coolant of the whole region below `r ~ 1.01`, and the `w = 0.1` guard
only cut it down to that value, not away. It is now `1.0e-11`. The layer warms,
the density at `1.01 R_p` rises by a factor 3.3, and the mass flux follows.
Outside `r ~ 1.1` the temperature profile moves by less than 4%.

Note that this is *not* the same defect the cutoff was addressing. The guard
acts below `1e3` K; most of the C I residual removed here sits at 500-1000 K
where the guard was already suppressing it by `1e-11` and it was *still*
comparable to the heating.

**WASP-121 b** (`backup/regression/wasp_full`, the golden case, cold start):
`Mdot` at `2 R_p` moves from 3.4815e13 to 3.5494e13 g/s, `log10 Mdot`
13.5418 -> 13.5502, i.e. **+2.0%**; base `T` +0.7%, `T(1.01)` +7.0%, and the
profile outside `r = 1.05` moves by <=3% in `T` and <=8% in density. The
smaller shift is expected: at 2400 K the C I and N II densities at the base are
`5e7` and `1e4` cm^-3, and C II 158 um / O I 63 um already had a two-level
saturation, so what is left is the newly saturated `[O I] 145.5/44 um`
channels, the corrected level normalization (below) and the updated H-atom
collision rates.

### 7.4 Two things fixed on the way, and one still open

- **Double normalization.** The previous two-level treatment of [O I] 63 um
  and [C II] 158 um multiplied a two-level solution -- already normalized to
  the two-level partition function -- by the ground-term Boltzmann fraction of
  its lower level. Its LTE limit was therefore low by
  `1 + (g_u/g_l) exp(-E/kT)`, i.e. 2.7x for [C II] 158 um at the 530 K base of
  HD 189733 b. The statistical-equilibrium solution has no such factor and
  reproduces the exact multilevel LTE emission to 1e-14 relative.
- **H-atom collision rates.** The two-level branch used a constant
  `4.0e-11 cm^3 s^-1` for [C II] 158 um + H and `4.2e-11 (T/100)^0.67` for
  [O I] 63 um + H, both flagged approximate in the source. They are replaced by
  fits to published quantum scattering calculations -- Yan & Babb (2023) for
  C I and N II, Abrahamsson, Krems & Dalgarno (2007) for O I, Barinovs et al.
  (2005) for C II -- which are 5-20x larger.
- **Still open: the base has no radiative floor.** With the artifact gone this
  is plain rather than hidden. At the HD 189733 b base the total radiative
  cooling is `3.5e-8` against a heating of `3.4e-6`, a factor 100 imbalance;
  the base temperature is set by the inner boundary and by the flow, not by
  local balance. That is item 3 of section 6 and is unchanged by this work.

### 7.5 Regression

`make check` after the change: `wasp_full` and `wasp_he23off` differ (metals
on -- the numbers above), `mol_base_handoff` is byte-identical (metals off, so
no metal line cooling is evaluated at all). The goldens were **not**
re-snapshotted here.

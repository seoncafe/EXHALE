# Where the post-process advection correction is valid

> **STALE since 2026-09-09 (N11, N11b).** This account describes the single
> `adv_status` column, the enthalpy-ratio refusal and the counts of that
> time. The current product has two status fields (`adv_T_status`,
> `adv_comp_status`), the measure column `adv_mass_row`, the row decision
> by the certification's face-flux mass operator against
> `adv_conditional_tol = 1e-2`, and the header block `# adv_schema 2`. The
> current description is in `README.md`, `README_HOWTO.md` and the user manual
> section 4. The measured
> numbers below are of the old product and are kept for the record.

2026-08-12. Files: `src/modules/post_process/post_process_adv.f90`,
`src/modules/nonlinear_system_solver/ion_cell_state.f90`,
`src/modules/nonlinear_system_solver/System_implicit_adv_{H,HeH,HeH_TR}.f90`.

## Verdict

Two defects, both in the advection-corrected post-process (`*_adv.txt` only;
the converged wind and the equilibrium `*.txt` files are produced before
`post_process_adv` runs and are untouched by anything here).

1. **The electron density inside the advection residuals was wrong.** It
   counted only the H and He electrons, `xe = x_HII + (He/H)(x_HeII +
   2 x_HeIII)`, while the equilibrium residual it is supposed to correct counts
   the metal electrons as well (`metal_electron_sum`, `ion_residual_core.f90`).
   In the shielded base the metals are the *dominant* electron donors --
   measured `n_e,metal/n_e = 0.48-1.00` in the base cells of all four paper
   planets -- so the recombination terms of the advection residual were low by
   up to six orders of magnitude there. This is a physics
   error independent of its size: the same routine already used the
   metal-inclusive `n_e` for its heating and cooling, so the two halves of the
   post-process disagreed about how many free electrons the gas has.

2. **The correction was being applied where it cannot be computed and where it
   is not needed.** In a shielded base the equilibrium ion fraction reaches
   `x_HII ~ 1e-13`, four to five orders of magnitude below the solver's
   absolute resolution, and those cells are simultaneously in local ionization
   equilibrium to within a factor `Da = 350-2900`. The advection correction
   there returns numerical noise of arbitrary sign in place of an answer that
   is already known.

Both are fixed. The `_adv` files of all four paper planets are now free of
negative densities and the affected base cells sit exactly on the converged
equilibrium solution.

## 1. The electron density of the advection residuals

`System_implicit_adv_{H,HeH,HeH_TR}` build the electron density of a cell from
the H/He ionized fractions alone. `post_process_adv` meanwhile computes `ne`
with `calc_ne(nhii,nheii,nheiii,ne,nm_w)`, i.e. including the metal electrons,
and hands *that* to the temperature solve and the cooling. The equilibrium
ionization systems (`System_HeH_metals`, `System_HeH_TR_metals`) also include
them, through `metal_electron_sum`.

Measured metal fraction of the electron density in the first 16 cells above the
base, from the converged equilibrium profiles of each planet (`Ion_species.txt`,
X+ counted once and X++ twice, exactly as `metal_electron_sum` does):

| planet | `n_e,metal/n_e`, base | `n_e,metal/n_e` at `r = 1.3 R_p` |
|---|---|---|
| HD 209458 b | 1.000 | 0.0010 |
| HD 189733 b | 0.997-0.998 | 0.0000 |
| WASP-52 b | 0.996-0.997 | - |
| WASP-121 b | 0.477-0.518 | 0.0013 |

The three cold-based planets have essentially no H/He electrons at the base at
all; WASP-121 b, whose base is hot (`T = 2407 K`, `x_HII,eq = 3e-4`), still gets
half of its base electrons from the metals. In the wind the metals are a
0.1% correction, which is why the change is confined to the base.

The consequence for the ODE the correction solves is direct. With `u = x_HII`
small, the steady advection-ionization balance is
`u ~ P_HI / (alpha_HII n_e)`, so an electron density low by `10^6` puts the
fixed point of the residual `10^6` above the equilibrium value. Measured at the
HD 209458 b base (cell j = 4, `r = 1.00077 R_p`), solving the residual *as
coded* to machine precision with Brent gives `u = 4.74e-8` against an
equilibrium `x_HII,eq = 7.46e-14`, and the H/He-only electron density used
there is `6.6e-7` of the true one.

The fix adds one field, `adv_cell%xe_metal`, the metal electrons per H nucleus,
set per cell in `post_process_adv` and added to `xe` in the three advection
residuals. It uses the same definition as `calc_ne` and `metal_electron_sum`
(stage 1 contributes one electron, stage 2 two). Unlike `calc_ne` it is *not*
conditioned on `eos_include_metals`: that switch governs whether the metals
enter the gas mass and particle budget, whereas the recombination terms need
the true free electron density. It is zero when metals are off or when the
post-process runs metal-free (`pp_metals 0`), since `nm_w` is zero there.

## 2. The validity range of the correction

The correction replaces the local ionization balance of a cell by the steady
advection-ionization ODE, integrated upwind across the cell:

```
x_j - x_{j-1} = (dr_j / v_{j-1}) * [ source - sink ](x_j)
```

Four conditions make that replacement carry no information. Where any of them
holds, the cell keeps the converged equilibrium ionization.

**(i) Inflow, `v <= 0` on either face -- physical.** The upwind
discretization takes the upstream state from the cell below, which is not the
upstream cell when the gas moves inward (the breathing base of a
Roche-filling planet). The residence time `dr/v` is negative as well. This
condition already existed but was gated on `pp_metal_on`, i.e. it was switched
on and off by the metal-cooling flag; the discretization does not care about
metals, so the gate is removed.

**(ii) `Da = (dr/v) min(nu) > 100` -- physical.** The Damkohler number compares
the time the gas spends in the cell with the relaxation time of the level
populations. `Da >> 1` means they relax to local equilibrium many times over
while the gas crosses the cell, so the equilibrium solution *is* the solution
of the ODE and the correction can only add integration error. `n_e` here is
the metal-inclusive electron density, the same one the equilibrium solve used.
Measured `Da` in the cells that motivated this: 350-2900 at the HD 209458 b
base, up to 2500 at HD 189733 b.

`min(nu)` is the slowest relaxation rate of the species the advection system
solves -- H I/H II, He I/He II, He II/He III and, with the triplet on,
He(2^3S) -- each being the total rate at which that population is destroyed
and re-formed. The systems solve the whole H/He vector at once, so a cell may
be pinned to equilibrium only when every population it solves is equilibrated.
This condition originally read `Da = (dr/v) (P_HI + alpha_HII n_e)`, the
hydrogen rate alone, which froze the far more slowly relaxing He(2^3S)
metastable at its equilibrium value in cells where it is advected.

**(iii) `x_HII,eq < 1e-6` -- numerical.** The residuals carry the *neutral*
fraction `x_HI` and the ion density is extracted as `(1 - x_HI) n_h`, so the
ion fraction inherits the solver's *absolute* resolution on `x_HI`:
`xtol = sqrt(eps) = 1.49e-8`, which is also the forward-difference step of the
MINPACK Jacobian. Below `1e-6` the extracted ion fraction is worse than 1%
relative. In a shielded base with `x_HII,eq ~ 1e-13` it is quantized at `1e-8`
with an arbitrary sign: at HD 209458 b cell j = 4 `hybrd1` returns
`x_HII = -4.96e-8` where the exact root of the same residual is `+4.74e-8` --
the right magnitude to within 5% and the wrong sign. A negative ion density in
one cell is then carried outward by the upwind cascade, which is how 42
negative entries appeared across 14 cells of the stored 2026-08-11 HD 209458 b
output (48 in the re-run of the same configuration).

**(iv) The mass flux through the cell is not stationary -- physical.** Added
2026-09-08; stated and measured in "The refusal of a cell whose flow is not
stationary" below, which is also where the reason it is not covered by (i) or
(ii) is written out. Every equation this post-process solves is a steady
equation integrated along the recorded flow, and the steady internal-energy
equation says when it applies in its own terms: where the enthalpy flux of the
mass-flux divergence exceeds both terms the balance keeps, the cell has no
steady solution the gas realizes, and neither the ionization nor the energy
correction is a correction there.

Conditions (i), (ii) and (iv) are statements about the flow and (iii) about the
representation of the unknown, so none of them depends on whether metal
cooling is switched on. The guard is therefore unconditional.

Neither threshold is a tuning knob fitted to an output. `Da = 100` is two
orders of magnitude into the equilibrium-dominated regime; `x_HII,eq = 1e-6` is
the point at which the solver's fixed absolute resolution leaves less than 1%
relative accuracy. Neither condition can trigger in a wind that the correction
is meant for: at `r > 1.05 R_p` no cell of any of the four planets satisfies
either.

Clipping the extracted densities at zero was tested and rejected: it hides the
noise instead of removing it, and it breaks the H nucleus budget, since
`n_HI + n_HII = n_h` is imposed by construction and a clipped `n_HII` no longer
satisfies it.

### The energy equation

**Updated 2026-09-08.** The temperature loop applies condition (i) and, since
the refusal below, condition (iv); it also has a thermal Damkohler condition of
its own -- with
its `pp_metal_on` gate likewise removed, for the same reason. Conditions (ii)
and (iii) are statements about the ionization balance and its representation,
so they do not carry over; the statement the energy equation needs is a
*thermal* Damkohler number, and it is now evaluated
(`thermal_damkohler_number`, `post_process_adv.f90`):

    Da = t_cross * |heating - cooling| / u_th ,     t_cross = dr/v

the number of times over the local net radiative rate could rewrite the gas
internal energy while the gas crosses the cell. Above one the cell keeps the
temperature the run's own energy equation converged to, which is the root of
that same local balance with every channel the run solved; the post-process
carries fewer channels, so its own root there would be worse. Unity is the
statement that one term of the equation is larger than the other, not a
threshold with a value to choose. The pair is formed at the state each pass
starts from, since the condition has to be decided before the temperature is
solved for.

The run reports the largest Damkohler number it reached and where, so the
margin of the condition on a given state is in the log whether or not any cell
crossed unity. MEASURED 2026-09-08 on the last pass, over the cells the
condition tests (both faces outflowing), by post-processing the recorded state
of each case with the enthalpy flux of the mass-flux divergence carried (the
subsection below):

| case | largest `Da` | where | cells above one |
|---|---|---|---|
| `lower_profile` | 1.08 | 2.2217 `R_p` | 1 of 502 |
| `wasp_full` | 1.07 | 1.0327 `R_p` | 1 of 502 |
| `mol_base_handoff` | 0.048 | 1.2732 `R_p` | 0 of 502 |
| `wasp_full_newton` | 0.035 | 1.2629 `R_p` | 0 of 502 |
| `hydrostatic_column` | 3.4e-8 | 1.2160 `R_p` | 0 of 502 |

So the condition is a guard at the edge of a real wind rather than a switch
that removes the correction: at the time of that measurement the slowest
outflowing cell of the `wasp_full` base and one cell of `lower_profile` crossed
it and everything else was below it, `wasp_full_newton` by a factor of 30 and
the mechanical column by eight orders of magnitude. The numbers ADV-STATIC
measured before the enthalpy flux was carried (0.951 for `wasp_full`, 0.534
for `lower_profile`, 0.055 for `mol_base_handoff`, 3.0e-19 for the column) are
of the same size on the same cells: `Da` is formed at the temperature each pass
starts from, so restoring the term moves it through `u_th` and through the
radiative rates.

**Re-measured 2026-09-08 with the refusal (iv) in place**, which is tested
before this condition and takes the cells it would have held:

| case | largest `Da` | where | cells above one |
|---|---|---|---|
| `lower_profile` | 0.069 | 1.5524 `R_p` | 0 of 502 |
| `mol_base_handoff` | 0.049 | 1.2732 `R_p` | 0 of 502 |
| `wasp_full` | 0.039 | 1.2925 `R_p` | 0 of 502 |
| `wasp_full_newton` | 0.035 | 1.2629 `R_p` | 0 of 502 |
| `hydrostatic_column` | 4.3e-94 | 1.0004 `R_p` | 0 of 502 |

No case in the matrix now reaches it. The two cells that used to cross it are
the two whose flow is not stationary, so they are refused by (iv) first and
the thermal condition never sees them; it remains armed for a cell whose mass
flux IS stationary and whose radiative rate is nevertheless fast, which is the
state it was written for.

### The enthalpy flux of the mass-flux divergence

**Restored 2026-09-08.** The equation the temperature loop solves is now the
full steady internal-energy equation of the supplied profile,

    div(u v) + p div(v)  =  heating - cooling
      = rho v de/dr - p v dln(rho)/dr + h div(rho v)                     (E)

upwind-differenced, with `h = e + p/rho` the enthalpy per unit mass; the two
forms are one equation, through
`div(rho e v) = e div(rho v) + rho v de/dr` and
`p div(v) = (p/rho) div(rho v) - p v dln(rho)/dr`. Until 2026-09-08 the loop
solved (E) without its third term, that is with `div(v)` replaced by
`-v dln(rho)/dr`, which is the same equation only where the mass flux
`rho v r^2` is stationary. No case in the matrix is stationary by the
cell-to-cell measure, so the term was not a small correction anywhere.

The third term reaches the residual (`T_equation.f90`) through a named field of
the cell record, `teq_cell%div_rhov`, as a separate additive term: dividing the
residual by `mup*mum*dr` puts it in the form above with `e = E(x_H2,T)/mu` and
`p/rho = T/mu`, so the term to add is
`mum*dr*(E(x_H2,x) + x)*div_rhov` in the caloric branch and
`gamma_ad*mum*dr*x*div_rhov` in the monatomic one, which carries the extra
factor `gamma_ad - 1` of its own scaling. It is proportional to the unknown,
because `h` is a function of `T`. None of `mum`, `mup`, `coeff` is redefined,
so each still means what its name says, and a zero field contributes an exact
zero.

`div(rho v)` is formed with the operator the mass row of the state uses,
`(A_p F_p - A_m F_m)/dV` with `F = rho v`, `A = r^2` and `dV = d(r^3)/3`
(`RK_rhs`), over the control volume bounded by the two points the upwind energy
difference is taken between, `r(j-1)` and `r(j)`. Those are the only two points
at which the residual evaluates the state, so the term is an exact zero
wherever the two carry the same `rho v r^2`, and a stationary wind is left
where it was. The stored face mass flux (`face_mass_flux_of_state`) is not read:
it belongs to the last state whose mass row was assembled, which in a
post-process-only run is no state at all, and its faces are not the two points
of this difference.

The run reports the size of the term against the other two, in the residual's
own variables, and that ratio is also the refusal criterion of the section
below. MEASURED 2026-09-08 on the state the post-process is handed, over every
cell of it:

| case | largest term ratio | where | above one in |
|---|---|---|---|
| `hydrostatic_column` | 1.9e3 | 2.0287 `R_p` | 481 of 503 |
| `lower_profile` | 9.2e2 | 1.9131 `R_p` | 321 of 503 |
| `wasp_full` | 4.4e2 | 1.0320 `R_p` | 297 of 503 |
| `mol_base_handoff` | 2.1e2 | 4.7083 `R_p` | 33 of 503 |
| `wasp_full_newton` | 24 | 1.5705 `R_p` | 2 of 503 |

The reference point matters for comparing these with the numbers ADV-STATIC and
ADV-ENERGY reported (67 in 161 cells for `wasp_full`, 2.0 in 418 for the
column): those were formed on the last pass at the corrected iterate and only
over the cells that reached the measurement, which excludes the inflow cells of
a breathing base; these are formed once, on the recorded state, over every cell.
The ordering is the same either way, and `wasp_full_newton` is the smallest in
both.

**What the restored term moves.** MEASURED 2026-09-08 by post-processing each
case's recorded golden state with the binary before and after the change ("Do
only PP", so both post-process the same state). `Hydro_ioniz.txt` and
`Ion_species.txt` are byte-identical in every case, and so are the `rho` and
`v` columns of `Hydro_ioniz_adv.txt` and every metal column of
`Ion_species_adv.txt`; what moves is the corrected temperature, the pressure
that follows it, and the H/He ionization that the next pass solves at that
temperature:

| case | largest `T_adv` movement | where | largest H/He column movement |
|---|---|---|---|
| `wasp_full_newton` | 1.4 per cent | 1.5744 `R_p` | 1.9 per cent (He 2^3S) |
| `wasp_full` | 16 per cent | 1.0883 `R_p` | 16 per cent (He 2^3S) |
| `mol_base_handoff` | 301 per cent | 4.2886 `R_p` | 93 per cent (He 2^3S) |
| `lower_profile` | 456 per cent | 1.0411 `R_p` | 23800 per cent (He 2^3S) |

Every movement is relative to the value the pre-change binary wrote.

The ordering is the point: `wasp_full_newton` is the only Newton-converged
state of the four, and it is the one the term barely moves, which is what a
term proportional to the divergence of the mass flux has to do. `wasp_full`
stops on the `du` threshold; `mol_base_handoff` and `lower_profile` are
relaxation snapshots pinned to a step count, and their mass flux is not
stationary anywhere.

**`hydrostatic_column` shows what the state itself is worth, either way.** It
is a 300-step mechanical column with the radiation switched off, whose mass
flux `rho v r^2` runs over seven orders of magnitude across the domain
(4.3e11 to 9.9e18 in the run's own units, i.e. 58 times the median) and falls
by about 8 per cent from one cell to the next through the inner half of it.
Its `Hydro_ioniz.txt` is isothermal at 1084 to 1182 K. Without the third term
its `Hydro_ioniz_adv.txt` fell to 0.78 K at 1.396 `R_p`, the adiabat `T`
proportional to `rho^(gamma-1)` of the supplied density. With the third term it
rises instead, to 1.7e7 K at 1.40 `R_p` and 1.6e9 K at 2.96 `R_p`
(MEASURED), because the exact equation reads a mass flux falling 8 per cent
per cell as a compression at that rate and heats the gas accordingly. That is
the analytic solution of (E) for this profile: with heating and cooling
negligible (E) integrates to `w` proportional to
`rho^(gamma-1) F^(-gamma)`, `F = rho v r^2`, and the computed profile follows
it in trend, departing by 2.3 at the top of the domain over 502 cells of a
first-order upwind difference whose factor from one cell to the next is 1.08 (MEASURED).
So the number is the solution of the right equation on a state that does not
satisfy continuity, and neither 0.78 K nor 1.7e7 K is a temperature the gas
has. What removes it is refusing the cell, and that is the section below.

ADV-STATIC predicted 569 K at 1.40 `R_p` and 121 K at 2.96 `R_p` for the
restored equation. **That prediction is corrected here**: it was measured with
`p div(v)` put back in place of `-p v dln(rho)/dr` and with the other half of
the term, `e div(rho v)`, still absent, which is neither (E) nor the equation
that was solved before.

### The refusal of a cell whose flow is not stationary

**Added 2026-09-08.** The exact steady equation is
correct on a stationary state and unbounded on any other, so the correction is
now **refused** wherever its own assumption fails, and that cell keeps the run's
own temperature and the ionization equilibrium at that temperature -- exactly
what the thermal Damkohler condition already does where radiation rather than
the flow sets the state. The `_adv` product therefore means: **the exact steady
advective correction wherever a steady correction exists, and the run's own
state elsewhere.**

**The criterion is the term ratio of the equation being solved.** With the
common `1/dr` divided out, the three terms of (E) on one cell are, in the
variables of the residual,

    q_adv  = |rho v (e_j - e_{j-1})|
    q_prs  = |w v (rho_j - rho_{j-1})|            w = p/rho
    q_enth = |h div(rho v) dr|                    h = e + w

and the cell is refused where `q_enth > max(q_adv, q_prs)`
(`enthalpy_flux_term_ratio`, `post_process_adv.f90`). Unity is not a tunable: it
is the statement that one term is larger than the other. The inequality is
STRICT, so a flux that is stationary to round-off is never refused.

**Why this ratio and not `|dln F|` alone.** How far from stationary the mass flux
is matters only through the size of the term it puts into the equation, and the
same `|dln F|` is negligible in one cell and dominant in another according to how
large the two kept terms are there. MEASURED on `wasp_full`: the median
cell-to-cell `|dln F|` is 5.5e-3 (ADV-STATIC) and the ratio exceeds one in 297 of
503 cells, reaching 438.

**The refusal covers the composition as well as the temperature.** The
ionization correction of a cell is the same kind of object, the steady
advection-ionization ODE along the same flow, so a cell whose flow is not
stationary has no steady advective ionization correction either. Conditions (i)
and (ii) do not cover these cells: (i) tests the SIGN of `v` and (ii) the ratio
of the residence time to the ionization relaxation time, and both are satisfied
by a fast, well-directed outflow whose mass flux is nowhere stationary. MEASURED
2026-09-08, cells refused by (iv) that (i) to (iii) do not already refuse: 160 of
297 in `wasp_full`, 308 of 321 in `lower_profile`, 33 of 33 in
`mol_base_handoff`, 2 of 2 in `wasp_full_newton` and 0 of 481 in
`hydrostatic_column` (whose every cell is already at `Da_ion > 100`).

**Formed once, on the state handed in.** `rho`, `v` and `r` are the recorded
profile and the post-process never changes them, so the refusal is a property of
that state; it is evaluated on the first pass, at the run's own temperature and
the mean molecular weight of the composition the equilibrium solve returned.
Recomputing it per pass would let the refused set move with the iterate, and
every pass feeds its upwind cascade from that set, so the converged `_adv`
product would depend on the iteration history. This also changes the reference
point of the reported ratio: ADV-STATIC and ADV-ENERGY measured it on the last
pass at the iterate and over the cells that reached the measurement (inflow
cells excluded), which is why `wasp_full` reported 67 in 161 cells there and 438
in 297 cells here.

**What the refusal refuses, per case.** MEASURED 2026-09-08 by post-processing
each case's recorded state with the binary before and after the change ("Do only
PP", so both post-process the same state); `hydrostatic_column` was run in full
for its recorded 300 steps. `Hydro_ioniz.txt` and `Ion_species.txt` are
byte-identical in every case, and so are the `rho` and `v` columns of
`Hydro_ioniz_adv.txt`.

| case | 0 corrected | 1 `Da_th` | 2 not stationary | 3 ionization | 4 solve |
|---|---|---|---|---|---|
| `wasp_full_newton` | 501 | 0 | 2 | 1 | 0 |
| `mol_base_handoff` | 470 | 0 | 33 | 1 | 0 |
| `wasp_full` | 205 | 0 | 297 | 2 | 0 |
| `lower_profile` | 182 | 0 | 321 | 1 | 0 |
| `hydrostatic_column` | 0 | 0 | 481 | 23 | 0 |

(504 rows: the 500 physical cells and four ghost rows.) The ordering is the
result that matters: **`wasp_full_newton` is the only Newton-converged state of
the five and it is the one the refusal barely touches**, which is what a
condition on the divergence of the mass flux has to do. `wasp_full` stops on the
`du` threshold and its breathing base is refused; `mol_base_handoff` and
`lower_profile` are relaxation snapshots pinned to a step count; the mechanical
column is refused throughout.

| case | max &#124;dT_adv&#124;/T vs ADV-ENERGY | median over refused rows | max `T_adv` | max `T_run` | max `T_adv` before |
|---|---|---|---|---|---|
| `wasp_full_newton` | 1.9 per cent | 1.9 per cent | 11588 K | 11365 K | 11588 K |
| `wasp_full` | 10 per cent | 2.2 per cent | 11413 K | 11198 K | 11413 K |
| `mol_base_handoff` | 88 per cent | 35 per cent | 2955 K | 2615 K | 10652 K |
| `lower_profile` | 5.1 | 27 per cent | 7000 K | 7027 K | 16928 K |
| `hydrostatic_column` | 1.0 | 0.98 | 1182.32 K | 1182.32 K | 2.6e10 K |

The last three columns are the point of the change: the corrected temperature of
a state that does not satisfy continuity was unbounded above and is now bounded
by the state it was made from. `hydrostatic_column` falls from 2.6e10 K to
exactly the run's own maximum, with no row above it (435 rows were above it
before); `lower_profile` from 16928 K to below the run's own maximum, 78 rows
above it before and none now; `mol_base_handoff` from 10652 K to 2955 K, 84 rows
above the run's maximum before and 8 now, those 8 being cells whose flow IS
stationary and whose correction is therefore a correction. The two `wasp_full`
cases barely move at the top of the column because their non-stationary cells are
at the base.

Refusal is exact, not approximate: MEASURED, every refused row carries the run's
own temperature to the BIT (299 of 299 rows in `wasp_full`, 322 of 322 in
`lower_profile`, 483 of 504 in `hydrostatic_column` -- the 21 exceptions there
are rows refused by (ii) or (iii) alone, whose temperature is still corrected
because those conditions are about the ionization balance), and every row refused
by (iv) carries the equilibrium density of every species exactly. The one
exception is the He I column of `lower_profile`, which departs by 1.8e-16
relative: `pin_cell_to_equilibrium` writes the two neutral-helium populations and
their sum, so the summed column is a floating-point reassociation of the
equilibrium value rather than the same bits. That is a property of carrying the
singlet and the metastable separately and predates this change.

**The ionization columns move with the temperature they are solved at.** MEASURED
against the ADV-ENERGY binary on the same state: the largest movement of any
species column is 21 per cent (He 2^3S, `wasp_full`, 1.1308 `R_p`), 5.4 per cent
(He III, `wasp_full_newton`), 886x (He 2^3S, `mol_base_handoff`, 4.2886 `R_p`),
7.0x (He III, `lower_profile`) and exactly zero for `hydrostatic_column`, whose
`Ion_species_adv.txt` is byte-identical because its ionization was already pinned
to equilibrium in every cell.

**What it does to a spectrum.** MEASURED 2026-09-08, `EXHALE_transit.py` on
`wasp_full` before and after, same recorded state: 295 of the 500 rows inside the
impact-parameter range `1 <= b <= 1.5666 R_p` are refused, and the transit
metrics move at the 1e-4 relative level or less -- He 10830 red 4.102 per cent
both, blue 2.843 -> 2.844 per cent, red/blue 1.44 both, FWHM 0.899 A both. The
refused rows are in the base, where the disk-integrated depth is dominated by the
extended wind above. The largest movement of any line in that run is the O I
1302 triplet with rotation and instrument convolution, 0.5357 -> 0.5386 per cent.
So the refusal changes what the `_adv` file MEANS considerably and what this
particular spectrum predicts very little; a case whose refused cells sit in the
line-forming region would move more, which is why the count is printed.

**Every `_adv` golden of the matrix moves by design.** None was refreshed and
none was compared; `make check` and `run_check.sh` were not run.

## Measured effect

Runs: the four paper planets, from their own directories, `OMP_NUM_THREADS=16`,
each re-run in full. Previous outputs kept as `output_pre_ppadv_20260812/` and
`tpm_pre_ppadv_20260812/`.

### Mass-loss rate, negative densities, guard extent

Isolated A/B: the pristine `HEAD` binary and the modified one, run from the same
configuration. Their `Hydro_ioniz.txt`, `Ion_species.txt`,
`Cooling_breakdown.txt` and `Excited_H.txt` come out byte-identical -- only the
`_adv` files differ, as the call structure requires.

> **Dated 2026-08-13.** The `Mdot` column below is the A/B of two binaries on the
> 2026-08-11 production configuration and is a record of that measurement, not
> the current production set. The planet folders have since been re-converged
> under the ground-term fine-structure statistical equilibrium, the base-ghost
> composition fix and the H(n=2) rate corrections, and now give
> `log10 Mdot = 9.46` (HD 209458 b), `9.14` (HD 189733 b), `11.70` (WASP-52 b)
> and `13.20` (WASP-121 b) -- see each folder's `run_20260812_lyafix.log`. The
> finding this table supports (the `_adv` files change, the equilibrium files do
> not) is unaffected.

| planet | Mdot [log g/s] | negative `_adv` entries | cells held at eq | newly held | outermost held | max change, `r > 1.05` |
|---|---|---|---|---|---|---|
| HD 209458 b | 9.31 -> 9.31 | 48 -> 0 | 103 -> 124 | 21 | 1.0316 | 1.4% (He III) |
| HD 189733 b | 9.05 -> 9.05 | 0 -> 0 | 135 -> 212 | 77 | 1.0380 | 0.37% (He III) |
| WASP-121 b | 13.17 -> 13.17 | 0 -> 0 | 0 -> 0 | 0 | - | 3.8% (He III) |
| WASP-52 b | 11.63 -> 11.63 | 0 -> 0 | 3 -> 57 | 54 | 1.0106 | 3.3% (He III) |

Every cell whose equilibrium ion fraction is below `1e-6` is now held at the
equilibrium value (96, 35, 0 and 57 cells), and no newly held cell lies above
`r = 1.038 R_p`. WASP-121 b never triggers the guard: its base is hot and its
smallest equilibrium ion fraction is `2.9e-4`, well above the representability
limit, and no cell reaches `Da = 100`.

The two changes separate cleanly. A build carrying the guard alone (metal
electrons forced to zero) removes all 48 HD 209458 b negative entries by itself
and changes the profiles above `r = 1.05 R_p` by *exactly* zero on HD 209458 b,
HD 189733 b and WASP-121 b, and by `4.2e-5` on WASP-52 b (whose held cells do
feed the upwind cascade at that level). The 0.4-3.8% change in the wind is
therefore the metal electrons -- the physics correction -- to three significant
figures, and essentially none of it is the guard.

### Transmission spectra

Peak absorption depth of the rotation- and instrument-convolved profile, same
A/B:

| line | HD 209458 b | HD 189733 b | WASP-121 b | WASP-52 b |
|---|---|---|---|---|
| He I 10830 | -0.070% | -0.001% | +0.003% | -0.027% |
| Ly-alpha | +0.006% | +0.000% | -0.000% | +0.000% |
| H-alpha | -0.034% | +0.000% | +0.017% | +0.023% |
| H-beta | +0.017% | -0.008% | +0.036% | +0.013% |
| Mg II | +0.005% | +0.008% | +0.003% | +0.001% |
| Ca II | +0.004% | +0.020% | +0.002% | +0.000% |
| Na I D | +0.004% | +0.028% | +0.011% | -0.001% |

The largest move anywhere is 0.07%. Equivalent widths behave the same way (at
most 0.08%, H-alpha of WASP-121 b).

### A caveat about the stored 2026-08-11 outputs

The re-run reproduces the stored `output/` of WASP-121 b and WASP-52 b
byte-identically, but not that of HD 209458 b and HD 189733 b. Those two carry
`Load IC? True` and read `output/*_IC.txt`, and those files had been refreshed
to the *converged* state of the 2026-08-11 run. Re-running therefore continues
the convergence from a better starting point instead of repeating it: HD 209458 b
starts at `du = 2.79e-3` instead of `1.21e-2` and finishes at `3.13e-3` instead
of `1.199e-2`. `Mdot` is unchanged (9.31, 9.05), but the profiles move by up to
9% in density. **This has nothing to do with the change documented here** -- the
pristine `HEAD` binary reproduces the same new trajectory -- and it is why the
comparison above is an A/B between two binaries on the same configuration
rather than a comparison against the stored files. Comparing the regenerated
transmission spectra against the stored ones instead mixes the two effects in:
1.3% on HD 209458 b (H-beta) and 0.6% on HD 189733 b (Ca II), against 0.04% and
0.03% for WASP-121 b and WASP-52 b, whose wind did not move.

## Scope of the verification

Checked: `make check` (the goldens compare `Hydro_ioniz.txt` and
`Ion_species.txt`, which the post-process does not write, so byte-identity is
the expected and observed result); the four paper planets, end to end, wind and
`_adv` profiles and transmission spectra.

Not checked: runs with helium off (the `adv_implicit_H` path), runs with metals
off, and the molecular post-process. The guard and the electron term apply to
them by construction -- with metals off `xe_metal` is identically zero -- but no
such run was re-measured.

## An unrelated defect found in passing, not fixed (since fixed)

`System_implicit_adv_HeH.f90` (helium on, He 2^3S off) writes its electron-
impact ionization terms as `-(ghi + ionhi) x_HI` and so on, without the
`xe * n_h` factor. Every other form of the same physics in the code carries it:
`adv_implicit_H`, `adv_implicit_HeH_TR`, and the equilibrium rows
(`heh_rows`: `n_hi * b_hi * n_e`). The coefficients are rate coefficients
[cm^3 s^-1], so as written those three terms are also dimensionally
inconsistent with the photoionization rates they are added to. The path is
taken only when He is on and He 2^3S is off, so none of the four paper planets
uses it; fixing it changes that configuration's `_adv` output and is left to a
separate change. Marked at the code site.

**Fixed 2026-08-12** (commit `b0d44bc`, "adv_implicit_HeH: electron density in
the collisional-ionization terms"). Verified in the code: `adv_implicit_HeH` in
`src/modules/nonlinear_system_solver/System_implicit_adv_HeH.f90` now writes the
three rows as `-(ghi + ionhi*xe*n_h)*xhi + ahii*xhii*xe*n_h` and the He I / He II
equivalents, so the electron-impact coefficients carry the `xe*n_h` factor and
are dimensionally consistent with the photoionization rates. The paragraph above
is left in place as the record of the state at the time of the diagnosis.

## Raw diagnostic data

The diagnosis was carried out with an instrumented scratch build (environment-
switched candidate cures, full-domain dumps of the Damkohler number and the
electron budget) under

```
<session scratchpad>/ppadv_diag/
    src/          instrumented tree (candidates V1 guard, V2 metal electrons,
                  V3 clipping)
    runs/         the four planet configurations, post-process only
    var/<planet>_<variant>/   one output tree per candidate
    analysis/     stiffprof.py (Damkohler / electron budget per cell)
                  root.py      (exact Brent root vs. what hybrd1 returns)
                  table2.py    (negative entries, wind-region change per variant)
                  transit2.py  (equivalent widths per variant)
```

That tree is scratch and is not part of the repository; the numbers quoted
above are reproduced by the production runs recorded in this document.

---

# Two more, 2026-08-29: a population solved as a difference, and an unread `info`

Same file, same post-process, found the same way (a one-cell step in an `_adv`
metastable profile).

## 3. The ground singlet was solved as a difference

`System_implicit_adv_HeH_TR` carried the *summed* He I fraction as `x(2)` and
the metastable as `x(4)`, and built the ground singlet -- which every rate in
rows 2 and 4 needs -- as `xheiS = x(2) - x(4)`. Nothing in the system keeps the
two apart. Where helium is heavily ionized and the little neutral helium left
sits mostly in the metastable, the difference loses every significant digit and
finally evaluates to exactly zero; row 2 then no longer contains the unknown it
determines, and `hybrd1` stalls at `info = 4` with a residual of order `1e-2`.

Fixed by solving for the two populations, `x(2) = n(1^1S)/n_He` and
`x(4) = n(2^3S)/n_He`, and forming the summed He I as their **sum**. Row 2
becomes the ground-singlet balance (the summed row minus the metastable row),
which is the same system in different variables. `post_process_adv` carries
`nheiS` and `nheiTR` as the profiles it solves for and forms
`nhei = nheiS + nheiTR` where a routine wants the total; the only subtraction
left is at the entry point, where the equilibrium solution hands over a summed
He I and a metastable, and it is well conditioned there.

Measured, PP-only re-runs at `OMP_NUM_THREADS=1` against a baseline built from
the same tree without the change: `heh0p55_diff_ctrl` (the He-settled control,
`K_zz = 0`) goes from **310 non-converged advection cell solves of 4660 to 0**,
and its one-cell steps in `n(2^3S)` from 9 to 1; `heh2p13_diff_kzz1e9` and
`flux_closure/heh11p1/k01` move by `2e-9` at most, which is the roundoff of a
change of variables.

The neighbours of the stalled cells returned `info = 1` while carrying
`fvec(2) = 1.2e-2`: `hybrd1`'s `info = 1` tests the increment between iterates,
not the residual. The guard of item 4 below would not have caught them.

## 4. Nothing read the solver return code

`post_process_adv` calls `hybrd1` three times per cell -- the advection
ionization system, the metal stage re-solve (`pp_metals = 2`), the energy
equation -- and wrote the returned iterate out in every case. A non-converged
iterate satisfies neither the advection balance it was asked to solve nor the
equilibrium balance it started from, so it is not a state of the gas; the
equilibrium solution of that cell is, and it is already what the three validity
conditions above fall back to. Each call now keeps the equilibrium state on
`info /= 1` and counts it, summed over the ten passes (a cell reverted in an
early pass feeds that pass's upwind cascade whether or not the last pass
converges).

On the four LHS 1140 b runs checked the guard is quiet: nothing in the
advection system once item 3 is in, and one cell of the energy equation in two
of the ten passes of `heh0p55`, worth `1e-4` to `1e-3` in that run's `_adv`
profiles. The metal guard is untested by them -- `pp_metals` defaults to the
frozen mode, so the re-solve block never runs.

## 5. Undefined helium arrays with helium off

With `thereis_He = .false.` the post-process left `nhei`, `nheii`, `nheiii`,
`nheiTR` -- and, once the singlet is a profile of its own, `nheiS` -- undefined
until the helium-free branch of the ionization loop zeroes them at its *end*.
The first pass reads them before that, in `nhe`, in the electron sum and in
`eval_cool`. They are now zeroed where the helium arrays are initialized;
helium-on runs are bit-identical across the change.

## The two that were left open, and are now closed

The control run that exposed both was still unphysical after them, for a
different reason: `P_HeI` in its post-process reached `7e15 s^-1`. The SvS85
secondary channel was added as
`P_HeI = P_HeI + R_secHeI/max(nheiS, 1e-99)` (`util_ion_eq.f90`), a volumetric
rate carrying SvS85's own helium abundance divided by the actual singlet
density: it diverged as the singlet disappeared, and the ten post-process passes
closed the loop. In a `Do only PP: True` run the equilibrium solve additionally
ran with the coupling *off* and the post-process with it *on*, because
`EXHALE_main.f90` flips `sec_ion_active` after the time loop.

Both are fixed (2026-08-30). The branching is now resolved on the cell's own
neutral He/H ratio and returned per target atom, so nothing divides by a
vanishing density (`svs85_secondary_branching`;
`TO_BE_DONE.md` item (J)), and `sec_ion_active` is armed before the loop when
`do_only_pp` is set, so the one equilibrium solve of a PP-only run uses the same
physics as the post-process beside it (`TO_BE_DONE.md` item (K)). On the control run the cells above 1.5 R_p with a
collapsed ground singlet go from 171 of 233 to none.

Fixed in passing, outside the post-process: `LHS1140b/make_memo_figures.py` and
`LHS1140b/exhale/kzz_scan_table.py` formed the elemental helium density as
`HeI + HeII + HeIII + HeITR`, counting the metastable twice -- it is a level of
He I and is already inside the `HeI` column (`species_table.f90` says so in
capitals). Both now sum the three ion stages only. The stored figures and
tables move by at most `2e-6` relative (measured on the four `K_zz` runs), so
nothing was regenerated.

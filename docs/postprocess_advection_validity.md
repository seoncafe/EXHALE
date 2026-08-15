# Where the post-process advection correction is valid

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

Three conditions make that replacement carry no information. Where any of them
holds, the cell now keeps the converged equilibrium ionization.

**(i) Inflow, `v <= 0` on either face -- physical.** The upwind
discretization takes the upstream state from the cell below, which is not the
upstream cell when the gas moves inward (the breathing base of a
Roche-filling planet). The residence time `dr/v` is negative as well. This
condition already existed but was gated on `pp_metal_on`, i.e. it was switched
on and off by the metal-cooling flag; the discretization does not care about
metals, so the gate is removed.

**(ii) `Da = (dr/v) (P_HI + alpha_HII n_e) > 100` -- physical.** The Damkohler
number compares the time the gas spends in the cell with the H
ionization/recombination time. `Da >> 1` means the ionization state relaxes to
local equilibrium many times over while the gas crosses the cell, so the
equilibrium solution *is* the solution of the ODE and the correction can only
add integration error. `n_e` here is the metal-inclusive electron density, the
same one the equilibrium solve used. Measured `Da` in the cells that motivated
this: 350-2900 at the HD 209458 b base, up to 2500 at HD 189733 b.

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

Conditions (i) and (ii) are statements about the flow and (iii) about the
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

### What the guard does *not* cover

The temperature loop of `post_process_adv` applies condition (i) only -- with
its `pp_metal_on` gate likewise removed, for the same reason. Conditions (ii)
and (iii) are statements about the ionization balance and its representation;
the corresponding statement for the energy equation would need a thermal
Damkohler number, which is not evaluated. Where the ionization is held at
equilibrium but the flow is an outflow, the energy equation is still solved
with advection.

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

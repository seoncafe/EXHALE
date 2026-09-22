# Photochem solver modification investigation

Date: 2026-08-29

> **Implementation status:** The recommended changes were subsequently
> implemented in the local Photochem `v0.9.0` source tree. Build details and
> measured validation results are recorded in
> `photochem_solver_modification_implementation.md`.
>
> Two corrections to what is written below, both established by measurement
> after this document was written:
>
> - **Recommendation 4 of section 2.5 is wrong as stated.** Comparing the
>   solver's *total* elemental vector `gas.molfracs_atoms` with the request at
>   every level rejects every composition, including ones that ran on 0.8.4.
>   The correction is marked at that item and its mechanism is section 3 of the
>   implementation record.
> - **The handoff tolerance did not stay at `1e-3`.** With the corrected solver
>   in production the `--abundance-tol` default went to its design value
>   `1e-10`, as section 2.4 anticipated it could. See section 6 of the
>   implementation record.

## Executive conclusion

The two recorded failures are real, but
they do not have the same owner.

1. **The deep elemental residual can be fixed through the Photochem path.**
   The immediate control already exists: `ChemEquiAnalysis.mass_tol`. Setting
   it to `1e-12` before `composition_at_metallicity` is called reduced the
   largest N/H residual over all 38 saved LHS 1140 b deep states from
   `2.018e-4` to `3.478e-8`, with all 38 equilibrium solves still reporting
   convergence. The proper upstream correction belongs in **Equilibrate**,
   because its convergence test scales every elemental residual by the largest
   elemental abundance instead of the abundance of the element being tested.
   Photochem should additionally stop accepting a result when that stricter
   test fails.

2. **The three climate failures are not chemical-equilibrium failures.** They
   occur in the **clima** library included by Photochem. The climate driver
   sends two unscaled residuals with very different units and magnitudes to an
   unconstrained MINPACK `hybrd1` solve. A failed trial evaluation aborts the
   complete solve. The same models have valid radiative-balance roots when the
   temperature solve is bracketed. The proper correction is a bounded,
   scaled climate solve in clima, preferably nested one-dimensional bracketed
   solves. Changing EXHALE's initial temperature guesses is not a reliable
   fix.

Therefore, modifying only the Photochem Python repository can provide a good
equilibrium workaround, but a complete upstream solution requires changes to
all three repositories:

| Repository | Required responsibility |
|---|---|
| Equilibrate | Define elemental convergence with an element-relative test |
| Photochem | Request strict equilibrium closure and reject failed or poorly closed results |
| clima | Replace or harden the unconstrained climate and background-pressure root solves |

No source code was changed for this investigation. The recommended changes
below are designs, not implemented patches.

## 1. Implementation and version audit

The production environment used by the saved diagnosis imports:

```text
photochem 0.8.4
equilibrate 0.2.2
clima 0.7.4
```

The local `photochem/` source tree is clean at commit
`e1e872528b61e8dd1db891b84869738e219b4f97` (`v0.9.0`). Its build configuration
still selects Equilibrate `v0.2.2` and clima `v0.7.4` in
`photochem/src/dependencies/CMakeLists.txt`.

The upstream repositories were checked on 2026-08-29:

- Photochem `origin/main` is still the local `v0.9.0` commit. Its current
  `origin/dev` retains the same `composition_at_metallicity` behavior and the
  comment `# Do not enforce convergence.`
- Equilibrate `main` is still version `0.2.2`; the relevant solver files are
  identical to the `v0.2.2` release.
- clima `main` identifies itself as `0.7.5`, but its relevant
  `src/adiabat/clima_adiabat.f90` is identical to `v0.7.4`.

Thus, upgrading to the current upstream branches does not remove either
problem.

Primary upstream sources:

- <https://github.com/Nicholaswogan/photochem>
- <https://github.com/Nicholaswogan/Equilibrate>
- <https://github.com/Nicholaswogan/clima>

## 2. Deep elemental closure

### 2.1 Exact call path

The EXHALE adapter follows this path:

```text
photochem_to_lower_profile.py
  -> EvoAtmosphereGasGiant(...)
  -> pc.gdat.gas.molfracs_atoms_sun = requested elemental vector
  -> pc.initialize_to_climate_equilibrium_PT(...)
  -> gasgiants.composition_at_metallicity(...)
  -> ChemEquiAnalysis.solve(...)
  -> Equilibrate CEA convergence test
```

`EvoAtmosphereGasGiant._initialize_atmosphere` then applies a fixed partial
pressure lower boundary to every initialized gas. This confirms the original
diagnosis: the deep Photochem cell preserves the equilibrium initializer. The
photochemical kinetics and transport solve cannot repair an elemental residual
already present there.

### 2.2 The convergence criterion that produces the residual

`ChemEquiAnalysis` exposes `mass_tol` to Python and initializes it to `1e-6` in
`Equilibrate/src/equilibrate.f90`. Both the short and long update paths in
`Equilibrate/src/equilibrate_cea.f90` calculate

```fortran
mval_mass_good = MAXVAL(b_0)*self%mass_tol
...
IF ((abs(b_0(i_atom) - mass) > mval_mass_good) .AND. &
    (b_0(i_atom) > 1d-6)) THEN
  mass_good = .FALSE.
END IF
```

This has two important consequences.

1. The allowed absolute error for N is scaled by the dominant element, not by
   N. A single `mass_tol` therefore implies very different relative accuracy
   for H, He, C, N, and O.
2. An element whose normalized abundance is at most `1e-6` is omitted from the
   mass-convergence decision entirely.

Consequently, `converged=True` does not mean that every requested element is
closed to `mass_tol` in relative terms. It means that the result satisfies the
current global absolute criterion. The observed N/H behavior is consistent
with the implementation; it is not a false Boolean returned after failure of
that implementation's own test.

The subsequent normalization of `molfracs_species_gas` and
`molfracs_atoms_gas` in `Equilibrate/src/equilibrate.f90` does not create the
problem. The elemental residual is already in the species amounts returned by
the CEA iteration, as the saved diagnosis measured.

### 2.3 Direct tolerance experiment

The installed solver was called at the state used by the existing standalone
probe, `P = 1.673098e7 dyn/cm2`, `T = 428.10 K`, and reservoir He/H = 8.6.
Only `gas.mass_tol` was changed:

| `mass_tol` | `converged` | relative N/H residual |
|---:|:---:|---:|
| `1e-4` | true | `+7.464195e-5` |
| `1e-6` | true | `+7.464195e-5` |
| `1e-8` | true | `+7.464195e-5` |
| `1e-10` | true | `+2.786325e-9` |
| `1e-12` | true | `+2.786325e-9` |
| `1e-14` | true | `-1.221245e-14` |

The discontinuous improvement is expected from a stopping test: the previous
iterate is accepted until the global threshold becomes smaller than its
elemental residual.

A second test used the actual deepest `P` and `T` read from each of the 38
saved profiles under
`LHS1140b/exhale/nh_refusal_diagnosis/scan/`. Each case used a new equilibrium
object and its own saved thermo file.

| `mass_tol` | converged cases | largest absolute N/H residual | median absolute N/H residual |
|---:|---:|---:|---:|
| `1e-6` | 38/38 | `2.017761e-4` | `1.008858e-5` |
| `1e-10` | 38/38 | `1.643215e-6` | `4.418635e-9` |
| `1e-12` | 38/38 | `3.477772e-8` | `1.368947e-10` |

This directly demonstrates that a strict `mass_tol` solves the EXHALE refusal
in the measured range. It does not yet prove that `1e-12` succeeds for every
thermochemical network, condensate transition, planet, or pressure-temperature
profile.

### 2.4 Recommended immediate EXHALE change

The least invasive production correction is to set the existing public
property before initialization:

```python
pc.gdat.gas.mass_tol = args.equilibrium_mass_tol
pc.initialize_to_climate_equilibrium_PT(Pc, Tc, Kc, 1.0, 1.0)
```

Recommended location:
`EXHALE_v1.00/src/utils/photochem_to_lower_profile.py`, immediately after
`EvoAtmosphereGasGiant` is constructed and before
`initialize_to_climate_equilibrium_PT` is called.

Recommended initial value: `1e-12`, exposed as an adapter option such as
`--equilibrium-mass-tol` and written to the profile provenance. The measured
38-case result supports this value for the current H/He/N/O/C mechanism.

This is preferable to relying only on the relaxed handoff
`--abundance-tol=1e-3`. The relaxed check is a defensible guard against the
current external solver residual, but it does not make the equilibrium state
elementally accurate. Once strict equilibrium closure is in production, the
handoff tolerance can be selected from the new measured residual rather than
from the old solver artifact.

That is what happened. On the corrected build the worst deepest-level
departure over 68 handoff writes is `3.28e-13` and the median `1.73e-14`, and
the `--abundance-tol` default went to its design value `1e-10`.

### 2.5 Recommended Photochem change

Photochem should make equilibrium accuracy explicit in the gas-giant API.
A suitable design is:

1. Add `equilibrium_mass_tol` to `GasGiantData` or to the
   `EvoAtmosphereGasGiant` constructor, with a strict gas-giant default.
2. Assign it to `self.gdat.gas.mass_tol` before any call to
   `composition_at_metallicity`.
3. In `composition_at_metallicity`, reject a level if all retries fail instead
   of unconditionally copying the last state.
4. After a reported success, independently compare the solver's total
   elemental vector, `gas.molfracs_atoms`, with the requested vector. Use the
   total vector for this test, not `molfracs_atoms_gas`, because a condensed
   phase may legitimately remove atoms from the gas at that level.

   **Corrected by measurement: this is right only where there is no condensed
   phase, and as written above it rejects every composition.** With
   `rainout_condensed_atoms=True` the vector handed to the next level is the
   gas-only one, while `gas.molfracs_atoms` started from the previous level
   through `use_prev_guess` still carries that level's condensate: measured on
   the LHS 1140 b column at He/H = 9, the total reads `1.1e+2` against the
   request where the same state solved from scratch closes to `5.7e-14`. The
   implemented test therefore forms the residual only at a level whose
   `molfracs_atoms_condensate` is entirely zero, which is every level deep
   enough to set an elemental handoff. Section 3 of
   `photochem_solver_modification_implementation.md`.
5. Reset `gas.use_prev_guess` in a `finally` block so an exception cannot leave
   the reusable solver in a different mode.

The independent validation is important even after Equilibrate is corrected:
it protects Photochem from future changes to a dependency's convergence
meaning.

### 2.6 Recommended Equilibrate correction

The physical quantity that must close is each nonzero input element. Replace
the single largest-abundance-scaled test with an element-relative test, with a
separate small absolute floor for numerical underflow. Conceptually:

```fortran
allowed_i = atom_mass_atol + atom_mass_rtol*abs(b_0(i_atom))
if (abs(b_0(i_atom) - mass) > allowed_i) then
  mass_good = .false.
endif
```

The exact API should expose relative and absolute tolerances separately. The
current `b_0(i_atom) > 1d-6` exclusion should not silently remove a requested
element from convergence. If very small abundances are intentionally ignored,
that cutoff should be a named option and the solver should report which
elements were excluded.

Both duplicated implementations must be changed:

- `ec_UPDATE_ABUNDS_SHORT` in `src/equilibrate_cea.f90`
- `ec_UPDATE_ABUNDS_LONG` in `src/equilibrate_cea.f90`

A shared elemental-closure routine would prevent the two paths from acquiring
different criteria later.

Required focused tests are:

- H/He/N/O/C at He/H from solar through the helium-rich LHS 1140 b range;
- the saved 300--2000 K temperature probe;
- elements on both sides of the current `1e-6` exclusion;
- gas-only and gas-plus-condensate cases;
- a check that `converged=True` implies the documented element-relative bound.

## 3. Radiative-convective root failures

### 3.1 Exact call path

The separate climate failure follows

```text
radiative_convective_column.solve_radiative_convective_column
  -> photochem.clima.AdiabatClimate.surface_temperature_bg_gas
  -> clima AdiabatClimate_simple_solver
  -> MINPACK hybrd1 on log10(T_deep) and log10(T_trop)
  -> TOA_fluxes_bg_gas
  -> make_profile_bg_gas
  -> a second MINPACK hybrd1 solve for background-gas pressure
```

The source is fetched from the clima repository during a Photochem build. It
is not implemented by the Photochem kinetic solver.

### 3.2 Weak points in the current clima implementation

`AdiabatClimate_simple_solver` in
`clima/src/adiabat/clima_adiabat.f90` gives `hybrd1` these residuals directly:

```fortran
fvec_(1) = ISR*rad_enhancement - OLR + self%surface_heat_flow
fvec_(2) = skin_temperature(...) - T_trop
```

The first is an energy flux and is typically of order `1e5 mW/m2` in the
failed LHS 1140 b cases. The second is a temperature difference of order
`1e2 K`. They are neither dimensionless nor similarly scaled. MINPACK therefore
sees a numerically ill-conditioned two-component residual whose norm is
dominated by radiative balance.

The variables are logarithms, which enforces positive temperatures but does
not enforce the physically required ordering or the valid range of every
thermodynamic polynomial. If any trial reaches a state that `make_profile` or
the nested background-pressure solve cannot evaluate, the callback sets a
negative `iflag`. MINPACK then terminates the complete solve instead of
rejecting that trial and remaining in the valid domain.

`make_profile_bg_gas` also uses unconstrained `hybrd1` for a scalar monotonic
problem in log background pressure. It tries only two initial scales, `1.0`
and `0.1`, and propagates an invalid profile evaluation as a fatal error.

These implementation properties explain both recorded exception forms:

- He/H = 8.1 and 15 end with the generic `hybrd1 root solve failed`, consistent
  with stagnation or lack of progress.
- Solar He/H ends through `make_profile_bg_gas` with `Failed to compute heat
  capacity`, showing that an intermediate trial left the thermodynamic domain.

The last point is an inference from the call path. The compiled library does
not report the trial temperature or species, so this investigation did not
identify the exact invalid trial.

### 3.3 A valid energy root exists

Changing the initial guesses over

```text
T_deep guess = 250, 300, 350, 400, 450, 500, 600, 800 K
T_trop guess = 80, 100, 120, 150, 180, 220 K
```

did not make the coupled He/H = 8.1 solve succeed. A multi-start wrapper is
therefore not an adequate correction.

However, holding `T_trop = 185 K` and solving only the radiative-balance
equation with a bracketed scalar method produced:

| reservoir He/H | bracket [K] | radiative-balance `T_deep` [K] | final absolute flux residual [mW/m2] |
|---:|---:|---:|---:|
| 0.0969 | 600--800 | 719.779874 | `3.21e-3` |
| 8.1 | 400--450 | 432.728611 | `1.80e-4` |
| 15 | 350--450 | 393.498939 | `6.38e-5` |

Nearby successful coupled models have `T_trop = 184.8--185.5 K`; for example,
He/H = 8.0 and 8.2 converge to deep temperatures of 433.3 and 431.5 K. The
bracketed He/H = 8.1 root lies smoothly between them.

This test proves that the recorded failure is not simply absence of a
radiative-balance solution. It does not constitute a complete replacement
climate solution because the skin-temperature equation was held fixed rather
than solved simultaneously.

### 3.4 Recommended clima redesign

The robust design is to use the problem's scalar structure instead of an
unbounded two-variable local solve.

1. For a trial `T_trop`, bracket `T_deep` over the common valid temperature
   interval of the active thermodynamic data and solve the radiative-balance
   equation with Brent's method.
2. Evaluate the skin-temperature residual at that radiative-balance solution.
3. Bracket and solve the remaining scalar `T_trop` equation.
4. Require `T_deep >= T_trop` in the admissible domain.
5. For `make_profile_bg_gas`, bracket the background partial pressure and use
   a scalar bracketed solve. Pressure positivity is then automatic and an
   isolated invalid trial does not destroy an otherwise valid bracket.
6. Report the bracket, final residuals, iteration counts, and the temperature
   and species when thermodynamic data cannot be evaluated.

This retains the public `surface_temperature_bg_gas` API while making the
physical domain and failure reason explicit.

If MINPACK must be retained temporarily, the minimum acceptable hardening is:

- nondimensionalize radiative balance by a representative stellar or thermal
  flux and the skin-temperature residual by a representative temperature;
- use variables that enforce `T_deep >= T_trop`;
- impose the intersection of all active thermodynamic temperature ranges;
- treat a failed trial as a rejected step rather than an immediate fatal exit;
- verify both physical residuals independently after `hybrd1` returns.

The bracketed formulation is preferred because it removes sensitivity to the
initial guesses seen in this scan.

### 3.5 How the fix enters Photochem and EXHALE

The permanent climate patch should be made and tested in the clima repository.
Photochem should then update the clima commit selected in
`src/dependencies/CMakeLists.txt`, rebuild its compiled extension, and expose
the improved diagnostics through `photochem.clima.ClimaException`.

EXHALE should not hide these failures by silently changing a composition or
accepting a nearby reservoir value. Its useful responsibilities are limited to:

- passing physical temperature bounds if clima exposes them;
- recording the solver and dependency versions;
- reporting the full diagnostic returned by clima;
- optionally offering a deliberately selected alternate initial guess for
  experiments, without treating it as the production fix.

## 4. Recommended implementation order

### Phase A: narrow EXHALE correction

1. Add `--equilibrium-mass-tol`, initially `1e-12` for the current mechanism.
2. Set `pc.gdat.gas.mass_tol` before equilibrium initialization.
3. Record the value in the handoff provenance.
4. Re-run only the equilibrium initialization and handoff checks for the 38
   saved states. A full EXHALE wind regression is not needed for this isolated
   initializer setting unless the resulting molecular partition changes enough
   to affect downstream inputs.

### Phase B: upstream equilibrium correction

1. Implement element-relative convergence in Equilibrate.
2. Add closure tests at trace abundances and condensation transitions.
3. Make Photochem enforce both the solver status and independent elemental
   closure.
4. Rebuild Photochem and repeat the targeted LHS 1140 b initialization scan.

### Phase C: upstream climate correction

1. Implement the bracketed background-pressure solve in clima.
2. Implement scaled, bounded climate temperature solves, preferably the nested
   bracketed formulation.
3. Add regression cases for He/H = 0.0969, 8.1, and 15, plus their successful
   neighbors.
4. Rebuild Photochem with the corrected clima commit.
5. Run the three climate columns and verify both radiative and skin-temperature
   residuals, profile validity, and smoothness against neighboring reservoir
   values.

## 5. Validation performed and not performed

Performed:

- inspected the EXHALE caller, Photochem gas-giant initialization, Equilibrate
  convergence and output construction, and clima root-solve callbacks;
- verified the installed and source dependency versions;
- fetched current upstream state and checked that the relevant behavior remains;
- measured equilibrium closure at six `mass_tol` values for one state;
- measured three tolerances over all 38 saved deep states;
- tested 48 coupled climate initial-guess pairs at He/H = 8.1;
- evaluated fixed-`T_trop` radiative residuals and found bracketed roots for all
  three failed compositions.

Not performed:

- no source patch was implemented or compiled;
- no complete photochemical steady-state column was regenerated with strict
  equilibrium tolerance;
- no corrected two-equation climate solution was produced;
- no full EXHALE wind or spectral calculation was run;
- no condensate-rich or sulfur-network tolerance sweep was run.

Those omissions are appropriate for the present request, which is an
implementation investigation. The measurements establish where to modify the
software and show that the proposed equilibrium control resolves the recorded
deep residual, while the proposed climate algorithm still requires
implementation and focused validation.

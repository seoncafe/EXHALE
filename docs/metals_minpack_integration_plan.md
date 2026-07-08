> **Status note (2026):** This plan has been **realized** in `EXHALE` (the renamed `EXHALE`): metals C/N/O are solved inside the coupled MINPACK 9-equation system, with Badnell RR+DR recombination and Kingdon&Ferland charge transfer. See `Update_EXHALE_early_phase`, Part II. The plan below is kept for historical reference.

# Plan: integrating trace metals into the MINPACK ionization system

## 1. Motivation

In `ATES_extended` the trace metals (C, N, O) are currently solved by a
**coronal-balance solver** (`metals_solve.f90::coronal_fractions`,
driven by `metals_drive.f90::solve_metals_post`) that is **decoupled**
from the H/He MINPACK system (`ion_system_H/HeH/HeH_TR`). The metals see
the H/He radiation field and electron density, but they do **not** feed
back into it. This is the standard "trace species" approximation and it
is adequate when metals are a negligible source of electrons.

That assumption breaks down when:

- metal abundance is high (super-solar) or the metals are strongly
  ionized, so that C/N/O contribute a non-negligible fraction of the
  free electrons;
- one wants the metal ionization to be *self-consistent* with the same
  electron density and radiation field that fixes H/He, rather than
  lagged/decoupled;
- one wants the metal ion fractions advected consistently in the
  time-loop (they are currently recomputed from coronal balance each
  step, not carried in `f_sp`).

Folding the metals into the MINPACK system removes these limitations: a
single nonlinear solve returns H, He, (He triplet) and all metal stages
simultaneously, with one shared electron density.

Note: the **physics already exists** in the sister tree
`EXHALE`, where a 7-equation `System_HeHCO` (HII, HeII, HeIII,
CII, CIII, OII, OIII) was built. This plan adapts that to
`ATES_extended` and additionally supports the He triplet and a variable
metal list.

## 2. Current vs. target architecture

| | Current (coronal) | Target (MINPACK) |
|---|---|---|
| Metal solve | independent `coronal_fractions` per ion | one coupled `hybrd1` system |
| Electron density | H/He only; metals lag | includes metal ions self-consistently |
| Metal ion fractions | recomputed each step, not stored | carried in `f_sp`, advected |
| He triplet + metals | independent (already works) | one combined system |
| `f_sp` width | 6 (H/He only) | 6 + 3·N_elem (or 2·N_elem if X2+ dropped) |

## 3. Variable layout of the combined system

Unknowns are ionization fractions relative to each element's total. With
the He triplet and `n_elem` metal elements (each carrying X+/X and
X2+/X):

```
x(1) = n_HII   / n_H
x(2) = n_HeII  / n_He
x(3) = n_HeIII / n_He
x(4) = n_HeITR / n_He         (only if thereis_HeITR)
then per metal element e (C, N, O):
x(k)   = n_(e^+)  / n_e_tot
x(k+1) = n_(e^2+) / n_e_tot
```

`N_eq = 3 + (1 if HeITR) + 2*n_elem`. The neutral fractions follow from
conservation. The free-electron density inside the residual is

```
n_e = n_HII + n_HeII + 2 n_HeIII
    + sum_e [ n_(e^+) + 2 n_(e^2+) ]
```

(He triplet is neutral, no electron contribution).

The residual equations reuse the existing balance template
`fvec = production - loss`:

- H, He, HeIII: as in `ion_system_HeH` / `ion_system_HeH_TR`.
- metal X0<->X+:  `n_X0*Gamma0 + (n_X0*C0 - alpha0*n_X+)*n_e`
- metal X+<->X2+: `n_X+*Gamma1 + (n_X+*C1 - alpha1*n_X2+)*n_e`

where `Gamma` are photoionization rates (already available:
`photoion_rate_metal`, `photoion_rate_metal_ion`), `alpha` the Badnell
recombination rates (`alpha_rec_metal`), and `C` optional collisional
ionization (currently neglected in coronal mode; can stay neglected
initially for a like-for-like comparison).

## 4. New / modified files

**New**
- `nonlinear_system_solver/System_HeHCO_TR.f90` — the combined residual
  routine `ion_system_HeHCO_TR(N_eq, x, fvec, iflag, params)`. Generalize
  `EXHALE/.../System_HeHCO.f90` to (a) read the active metal list
  and coefficients for each ion from `params`, and (b) optionally include the
  HeITR block.

**Modified**
- `init/parameters.f90` — none required for variables (params array
  already sized 40 in the metal work; confirm it is large enough for
  HeITR(7) + metals(up to 8 ions × {Gamma0,Gamma1,alpha0,alpha1,n}) and
  enlarge if needed).
- `files_IO/input_read.f90` — set `N_eq` for the combined case:
  `N_eq = 3 + merge(1,0,thereis_HeITR) + 2*n_elem`, and size the MINPACK
  work arrays accordingly. Requires `n_metals`/`n_elem` to be known
  before this point (move `read_metals_input` ahead of the allocation,
  or recompute `lwa` after it).
- `radiation/ionization_equilibrium.f90` — replace the metal cooling
  hook (which calls the coronal `solve_metals_post`) with extraction of
  metal densities from the MINPACK solution, and dispatch:
  ```
  if (thereis_metals) then
     call hybrd1(ion_system_HeHCO_TR, ...)   ! handles He, HeITR, metals
  else if (thereis_HeITR) then
     call hybrd1(ion_system_HeH_TR, ...)
  else
     call hybrd1(ion_system_HeH, ...)
  endif
  ```
  Fill `params` with the metal photoionization/recombination
  coefficients per cell; unpack `x` into `f_sp(:,7:)`.
- `init/set_IC.f90`, `files_IO/load_IC.f90`, `files_IO/write_output.f90`,
  `EXHALE_main.f90` — widen `f_sp` from 6 to `6 + 3*n_elem` and add the
  metal columns to `Ion_species.txt` (mirrors the EXHALE change).
- `radiation/metals_cool.f90` — unchanged; `eval_metal_cooling` now takes
  the metal densities from `f_sp`/the MINPACK solution instead of from
  `solve_metals_post`.
- `metals_drive.f90` — `solve_metals_post` becomes optional (kept for a
  post-processing-only / diagnostic mode, or removed once the MINPACK
  path is validated).

## 5. Implementation steps

1. **Port the residual.** Copy `System_HeHCO.f90` from `EXHALE`,
   rename to `System_HeHCO_TR.f90`, and (a) drive the metal block off the
   active `metal_list` so it works for any subset of C/N/O, (b) add the
   optional HeITR equation, (c) update the electron-density sum.
2. **Coefficient plumbing.** In `ionization_equilibrium.f90`, compute per
   cell the metal `Gamma0`, `Gamma1` (from `photoion_rate_metal*`),
   `alpha0`, `alpha1` (from `alpha_rec_metal`), and pack them into
   `params` in a documented layout.
3. **`f_sp` widening + N_eq.** Make `f_sp` width and `N_eq`/`lwa` depend
   on `n_elem`; thread the new columns through IC, output, and main.
4. **Dispatch + unpack.** Add the `ion_system_HeHCO_TR` branch; unpack
   the solution into H/He/HeITR and metal densities; keep the cooling
   call reading those densities.
5. **Build via `run_EXHALE.sh`** (the compile list must gain
   `System_HeHCO_TR.f90`).
6. **Compare to coronal.** Run the same planet through both paths and
   confirm they agree in the trace limit (low abundance), and document
   the difference at high abundance (where the coupling matters).

## 6. Risks and trade-offs

- **Convergence/stiffness.** A larger nonlinear system is harder for
  `hybrd1`; good initial guesses (propagate inward, as the current code
  does) and sensible scaling are essential. Watch the partially-ionized
  layer where Jacobians are stiff.
- **Cost.** `N_eq` grows from 3–4 to up to ~10; each cell solve is more
  expensive. The time-loop already dominates runtime, so profile before
  optimizing.
- **Backward compatibility.** Keep `n_metals == 0` reducing exactly to
  the present H/He(/HeITR) system (bit-for-bit), as the current code
  does. Guard every metal addition behind `thereis_metals`.
- **Two code paths during transition.** Until validated, keep the coronal
  `solve_metals_post` available (e.g. behind a flag) so results can be
  cross-checked.

## 7. Validation

1. **Reduction test:** `n_metals = 0` must reproduce the current H/He
   results bit-for-bit.
2. **Trace-limit agreement:** at low metal abundance the MINPACK metal
   fractions should match the coronal `solve_metals_post` output to the
   solver tolerance.
3. **Coupling test:** at high abundance, show the electron-density
   feedback changes the H/He ionization (the effect that motivates this
   work), and that the run still converges.
4. **He triplet + metals:** confirm the combined system reproduces the
   independent He-triplet result when metals are negligible.

## 8. Effort estimate

Roughly comparable to the EXHALE metal port: ~1 new module
(~150 LOC) plus edits to ~6 existing files for the `f_sp` widening and
dispatch. The coefficient and recombination data already exist
(`cross_sec_metals.f90`, Badnell rates in `metals_solve.f90`), so the
work is wiring and the nonlinear-solve robustness, not new physics.

## 9. Recommendation

Implement incrementally and keep the coronal path as a fallback:
1. start with C and O only, no HeITR, neglecting collisional ionization
   (closest to the current coronal physics) to validate the reduction
   and trace-limit tests;
2. add the HeITR equation;
3. add collisional ionization and N if needed.
This staging keeps each step independently testable and preserves the
verified behavior at every stage.

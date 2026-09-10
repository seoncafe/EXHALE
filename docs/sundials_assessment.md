# SUNDIALS adoption assessment for EXHALE

Date: 2026-08-10
Supersedes the recommendation in `ATES/ATES_sundials_solver_plan.md` (2026-06-08), which predates the in-house JFNK steady solver.

## Verdict

Do not adopt SUNDIALS (https://computing.llnl.gov/projects/sundials) now. The
place where it would have paid off most has since been filled by an in-house
implementation, and the remaining candidate slots gain nothing from it. It
becomes the leading candidate again the moment kinetic (time-dependent)
chemistry with species transport is adopted: see
`docs/kompot_comparison.pdf`, Sect. 5.2, item (d).

## Why the earlier plan is superseded

The 2026-06-08 plan ranked KINSOL applied to the *global steady-state
residual* as the top-priority use. That capability now exists in-house:
`src/modules/time_step/steady_newton.f90` implements pseudo-transient
continuation with a Jacobian-free Newton-Krylov solver (right-preconditioned
GMRES, SER ramp, non-monotone Grippo line search), plus a matrix-based PTC
Newton variant. It is the recommended path for a quantitative mass-loss rate
(marching warm-up followed by the Newton finish), selected per run by
`Solver: Newton`; the marching-only path remains the default. The plan's priority-1 item is therefore done
without SUNDIALS.

## Module-by-module check against the current code

| SUNDIALS module | Corresponding place in EXHALE | Assessment |
|---|---|---|
| KINSOL (nonlinear system, global) | `steady_newton.f90` JFNK | Already implemented in-house. |
| KINSOL (nonlinear system, cell-local) | Ionization equilibrium in each cell: analytic-Jacobian Newton with the in-tree MINPACK `hybrd1` fallback (`src/modules/nonlinear_system_solver/`) | No gain. The systems are small and dense (`N_eq = 4 + 2*n_melem` with `n_melem = 10`, i.e. 24 unknowns with the triplet and every metal on), so the current solver is close to optimal; the solve runs inside the OpenMP cell sweep, so KINSOL contexts would have to be created for every thread, adding overhead only. |
| CVODE (stiff ODE) | No corresponding place: EXHALE chemistry is local algebraic equilibrium, not time integration | Meaningless today; becomes the standard tool if a kinetic network is ever integrated in time. |
| ARKODE (IMEX time integration) | SSP-RK3 hyperbolic step + semi-implicit energy update | The goal is a steady state, so higher time accuracy has no practical value. |
| CVODES/IDAS (sensitivities) | None | Only relevant if sensitivity analysis becomes a goal. |

## Adoption cost (why "no gain" also means "net negative")

- EXHALE is currently fully self-contained: MINPACK is in-tree, the build is
  a plain Makefile, and three compilers are supported (gfortran/ifort/ifx).
  SUNDIALS is a CMake-built C library with Fortran 2003 interface `.mod`
  files, so it adds compiler-specific module compatibility and version
  management across every machine the code runs on.
- Replacing any solver changes the numerical path, which invalidates the
  byte-identical regression baseline (`make check` goldens would all need to
  be re-established).

## When SUNDIALS becomes the right answer

If the kinetic-chemistry path of `docs/kompot_comparison.pdf` Sect. 5.2 (d)
is ever pursued (time integration of a stiff reaction network with species
advection and diffusion), CVODE (BDF with a sparse Jacobian) is the standard
engine for the chemistry step: `photochem` is built on it. At that point the
comparison to make is CVODE versus a small in-house Rosenbrock integrator
(the route taken by Kompot with RODAS3 and by VULCAN). Until then, no part of
EXHALE poses a problem that SUNDIALS solves better than what is already in
the tree.

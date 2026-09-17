# `hydrostatic_column` -- a mechanical column, with the radiation switched off

## Purpose

The fixture for item D2 of `docs/development_plan_20260905_rev3.md`
(section 4.5, Phase 4): the discrete pressure-gravity residual and the
spurious velocity of the spatial operator, measured where nothing else can
produce a flow. Everything that puts energy into the gas or takes it out
across a cell face is off, so a departure from hydrostatic balance in this
run is a property of the discretization, not of the physics.

**The test itself is Phase 4 work and is not in this directory.** This case
pins the configuration and its output, and its verdict in the matrix today is
the ordinary bitwise one: the run reproduces its recorded profiles. What
would constitute a pass is stated below, so that the Phase 4 test can be
written against it rather than around it.

## What a pass would be

A converged mechanical column at cell count `N` gives, for every physical
cell,

    R_j = [ (p_{j+1/2} - p_{j-1/2}) / dr_j + rho_j g_j ]

and the case passes when

1. `max_j |R_j| / max_j |rho_j g_j|` falls at the design order of the
   reconstruction as `N` is doubled (`Grid cells:` 250 / 500 / 1000 / 2000 on
   this same input), and
2. the spurious velocity `max_j |v_j| / c_s(T_eq)` falls with it,

both measured on the relaxed state, not on a fixed step count. Neither number
is asserted here. Two ladders are needed, not one: with `Base grid [dr,cells]:
2.0e-4 50` held fixed, doubling `Grid cells:` refines the stretched region
(`dr_j` at r = 2 falls 18x from 250 to 2000) but leaves the base cell at
1.9e-4 to 2.0e-4 R_p, so the base term needs its own `Base grid` ladder
(measured 2026-09-06, `src/tests/grid_and_gates/hydrostatic_residual.f90`,
which also measures the pair on the ANALYTIC column without a relaxed state:
interior at design order, the O(1) term at the outer free-outflow cell under
PLM, see `docs/b4_spatial_operator_design_20260906.md` section 8).

## Configuration, and what each key switches off

| What | How | Note |
|---|---|---|
| Irradiation | `Use only EUV? True` and `Log10 of EUV luminosity [erg/s]: 10.0` | 1e10 erg/s at 0.048 AU is 1.5e-15 erg cm^-2 s^-1, i.e. 18 orders below the value the same geometry carries in `mol_base_handoff`. There is no key that removes the radiation field, so it is made negligible instead. |
| X-rays | `Use only EUV? True` | `LX` is then not read and the X-ray band is absent. |
| Conduction | `Conduction: False` | the default, stated so the case says what it holds fixed |
| Viscosity | `Viscosity: False` | the default, likewise |
| Molecular chemistry, metals, oxygen | keys absent, no `metals.inp`, no `base.inp` | atomic H and He only |
| He 2^3S | `Include He23S? False` | one fewer level to populate |
| Steady solver | no `Solver:` key | marching only; the case is a snapshot |

**Cooling cannot be switched off.** There is no input key for it: the
radiative losses are evaluated unconditionally in `Cool_coeff.f90` through
`eval_cool`, and `Energy solver` (K25 of `docs/input_schema.md`) selects the
explicit or the semi-implicit update of the energy equation, not whether the
sink exists. The case therefore holds the gas where the sink is negligible
instead: the hot-Uranus geometry at `T_eq = 1140 K`, where the gas is
essentially neutral and every atomic channel is far below its excitation
threshold. MEASURED (below), the column stays isothermal to 2 per cent over
300 steps, which is the check that this substitute worked.

The planet is the low-gravity hot Uranus of the Tier-2 gate (0.0457 M_J,
0.49 R_J, `T_eq` 1140 K, He/H 0.0793) and NOT a hot Jupiter, for a mechanical
reason: the isothermal scale height must be resolvable over a domain the code
can carry. Here `H/R_p = 0.044` and the escape parameter is
`lambda = GM mu m_H / (k T R_p) = 22.8`, so a hydrostatic column falls by
`exp(-lambda(1 - R_p/r))`, a factor 2.5e-7 out to the 3 R_p outer radius: from
1e14 cm^-3 at the base to 2.5e7 cm^-3 at the top, still above the density
floor. The same construction on HD 209458 b (`H/R_p = 0.010`, `lambda = 96`)
puts the whole domain above 1 R_p on the floor, and the floor value, not the
stratification, then sets what the operator is measured on -- MEASURED
2026-09-05 before this geometry was adopted: with the HD 209458 b parameters
the outer half of the domain sat at the floor 9.7e5 cm^-3 and carried
`|v|` up to 5.4e5 cm s^-1.

## What the state at 300 steps actually is (MEASURED 2026-09-05)

Single-threaded, gfortran 16.2, binary `e779298d0b482bfbc40ae2f9aa80d60a`,
4.9 s wall, `final: count=300  du= 4.1979E+00`.

| r [R_p] | rho [m_H cm^-3] | v [cm s^-1] | T [K] |
|---|---|---|---|
| 1.0002 | 1.217e+14 | 8.18e+04 | 1137.9 |
| 1.0218 | 6.137e+13 | 1.08e+05 | 1115.8 |
| 1.1367 | 1.227e+12 | 1.10e+05 | 1127.2 |
| 1.6777 | 1.606e+07 | 1.43e+05 | 1121.7 |
| 2.9600 | 9.736e+05 | 2.81e+05 | 1120.0 |

The temperature is flat to 2 per cent, which is what switching the radiation
off was for. The velocity is not zero: `max |v|` over the physical cells is
2.8e5 cm s^-1, against a sound speed of 3.5e5 cm s^-1 at 1140 K, so the column
is draining at a substantial fraction of `c_s` through the outflow boundary at
3 R_p. **This state is not a hydrostatic solution and the case does not claim
one.**

Marching longer does not produce one. MEASURED at 3000 steps on the same
input: `du` grows from 4.20 to 2.13e3, the layer between 1.02 and 1.2 R_p
cools to 524 K (adiabatic expansion, with no heating to balance it), the base
velocity reverses to -2.5e4 cm s^-1 and a 9.6e5 cm s^-1 spike stands at
1.68 R_p. The column drains rather than settles. 300 steps was chosen because
it is the state in which the intended isolation still holds; whether a
mechanical column has a discrete steady state at all under the present
boundary treatment is a Phase 4 question (D3, the boundary count) and is
recorded here as an open one.

## An observation outside this task's scope

With the flux made negligible the run prints

    He recombination coupling WARNING: 6132 cell-step(s) with
    dP_HI(He rec)/P_HI up to 1.31E+80 (threshold 1.0E+03)

`ionization_equilibrium.f90` 947 tests `dP_HI_hrc(j) > he_rec_dominant_ratio *
P_HI(j)` and 952 reports `dP_HI_hrc(j)/max(P_HI(j), 1.0d-99)`. With no
stellar photons `P_HI` is at the 1e-99 guard, so the test is true in every
cell and the reported ratio is the guard, not a measurement. It is a
diagnostic-only artifact of the zero-flux limit, it produces no NaN, and it
was not changed: production source is outside this task.

## Provenance

Built 2026-09-05 for Phase 0 item 4 of
`docs/development_plan_20260905_execution.md`, against
`docs/development_plan_20260905_rev3.md` section 3.3 ("Hydrostatic column
(mechanical)") and section 4.5 (D2). Since the 2026-09-07 series gate it IS
in `DEFAULT_CASES` (`run_check.sh`) and `golden/hydrostatic_column/` holds
its four files and `run.log`; the pre-series reference sits in
`baseline_post170_20260905/hydrostatic_column/`. Its `_adv` golden pins the
post-process's temperature, which item ADV-STATIC showed to be the discrete
adiabat of the dropped enthalpy-flux term (0.78 K at 1.40 R_p) until item
ADV-ENERGY restores that term.

# `carrier_reference_scales`

What the transported carriers of `diffusive_photochemistry.f90` are scaled
against, and whether their chemistry rows can be differentiated at all.

## What is asserted

| assertion | quantity | reference | tolerance |
|---|---|---|---|
| `carrier_reference_*` | `carrier_element_reference_density` of each of the five carriers | the free density of the element that carrier is made of (H2 and H+ hydrogen, OH and H2O oxygen, CO the scarcer of oxygen and carbon) | 1e-12 relative |
| `free_oxygen_density_is_zero`, `free_carbon_density_is_zero` | the column the derivatives are taken on | 0 | exact |
| `carrier_step_on_hydrogen_floor_*`, `carrier_step_follows_carrier_*` | `carrier_source_derivative_step` for the proton | 1e-12 of the free hydrogen density where the proton is smaller than that, 1e-6 of the proton where it is larger | 1e-15 relative |
| `proton_row_diagonal_derivative_*` | `d src(H+)/d n(H+)` by forward difference with the production step, at n(H+) = 0, 1e-12 n_H and 1e-3 n_H | a central difference with a step of 1e-7 n_H, chosen in the driver and read from no production routine | 1e-3 relative |

The tolerance of the last row is set by the cancellation floor of a forward
difference, `eps |f| / (dn |f'|)`, which is about 4e-5 relative at these
densities; MEASURED 3.7e-5 and 4.5e-5. A step that is not scaled to the
element cannot come near it: with the proton charged to the oxygen
reservoir, which is identically zero in a hydrogen and helium atmosphere,
the same two points return a derivative of exactly 0.

## How to run

    src/tests/carrier_reference_scales/run.sh

It builds into `build/tests/carrier_reference_scales/` and prints one
`PASS|FAIL <name> measured= reference= tol=` line per assertion, exiting
nonzero if any fails.

## Why it is a directory of its own

The driver reaches the frozen cell state `bg_cell`, which lives in
`ionization_equilibrium`. The source closure of that module reaches the
equilibrium solvers, and those call MINPACK (`hybrd1`, `hybrd`, `dpmpar`)
and LAPACK (`dgesvd`) as free subroutines. `source_closure.py` follows `use`
statements, so it cannot find them, and the `physics_probe` link line
carries no library. This suite's `run.sh` adds those sources and the LAPACK
the Makefile resolves; nothing else differs.

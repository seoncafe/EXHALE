# Species mass tests

One species, one mass. `bsp_mass` of `src/modules/init/species_table.f90` is
the definition of what every species weighs, in units of the hydrogen atom
`mu = 1.67353284e-24 g` (`parameters.f90`), and `calc_rho` forms the gas
density from it. Any module that restates a mass instead of reading that table
can drift away from it, and did: decision 14 of
`docs/development_plan_20260905_rev3.md` section 10.5 moved helium from 4.0 to
3.9715 `m_H` (the helium-4 atom over the hydrogen atom), and four modules kept
their own 4.0.

This suite holds the assertions that need more of the production tree than the
drivers of `src/tests/physics_probe` link. The other three copies of the helium
mass are asserted there, in `atomic_mass_and_radius_constants.f90`; the one
here reaches the ionization sweep through the carrier transport and therefore
needs the MINPACK subprograms and LAPACK, which that script does not link.

## Running

```bash
src/tests/species_masses/run.sh
```

Objects, module files and executables go to `build/tests/species_masses/`,
which the script creates. It needs nothing else built first. `FC` selects the
compiler (default `gfortran`); the LAPACK is the one in that compiler's own
prefix, by the same rule the `Makefile` follows, because a LAPACK built
against another `libgfortran` runtime resolves the symbols and then misbehaves.
The verdict format and the `assertion_report` module are shared with
`src/tests/physics_probe`.

## Files

| File | Contents |
|---|---|
| `run.sh` | Builds into `build/tests/species_masses/`, runs every driver, exits nonzero on any failure. |
| `carrier_background_masses.f90` | Group 1 below. |

## 1. Background collision partners of the carrier transport (`carrier_background_masses.f90`)

Origin: item 2c-MASS2, closing the copies decision 14 left behind.
Production data: `carrier_background_mass` of
`src/modules/lower_atmosphere/diffusive_photochemistry.f90` against `bsp_mass`.

The molecular carrier transport (`Molecular carrier transport: True`) moves a
carrier through a background of nine species by Blanc's law, and the
coefficient of each pair is the rigid-sphere coefficient of Banks & Kockarts
(1973), `D = 1.52e18 (1/A_s + 1/A_t)^(1/2) sqrt(T)/n`, whose `A_s` and `A_t`
are the two masses in units of the hydrogen atom. The same masses set the
gravitational settling drift. If they are not the masses `calc_rho` weighs
those species with, the friction and the density describe different gases.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `hydrogen_i_mass_of_background` | `carrier_background_mass(1)` | `bsp_mass(H I)` | exact |
| `hydrogen_ii_mass_of_background` | `carrier_background_mass(2)` | `bsp_mass(H II)` | exact |
| `h2_mass_of_background` | `carrier_background_mass(3)` | `bsp_mass(H2)` | exact |
| `h2_plus_mass_of_background` | `carrier_background_mass(4)` | `bsp_mass(H2+)` | exact |
| `h3_plus_mass_of_background` | `carrier_background_mass(5)` | `bsp_mass(H3+)` | exact |
| `helium_i_mass_of_background` | `carrier_background_mass(6)` | `bsp_mass(He I)` | exact |
| `helium_ii_mass_of_background` | `carrier_background_mass(7)` | `bsp_mass(He II)` | exact |
| `helium_iii_mass_of_background` | `carrier_background_mass(8)` | `bsp_mass(He III)` | exact |
| `heh_plus_mass_of_background` | `carrier_background_mass(9)` | `bsp_mass(HeH+)` | exact |

MEASURED, RED before item 2c-MASS2 on the last four: the module's own list gave
helium 4.0 and HeH+ 5.0 against the table's 3.9715 and 4.9715, 0.72 and 0.57
per cent high. The five hydrogen rows passed already. GREEN after all nine.

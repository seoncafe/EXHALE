# Constrained-network layout tests

What the constrained chemical-equilibrium continuation builds for one cell:
which species are unknowns of the Newton system, how many reaction and
conservation rows it has, and whether the two counts agree.

The driver links the **production** sources under `src/modules`, never a copy
of them. It calls `constrained_network_layout_of_cell`
(`src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90`),
which runs the production `set_molecular_network_layout` and
`set_rung_partition` on the cell standing in `ieq_cell` and reports what they
built.

This is a suite of its own rather than a driver of `src/tests/physics_probe/`
because its closure needs two things that suite's does not: the in-tree
MINPACK routines, which are plain external subroutines that no `use`
statement reaches, and LAPACK, for the opt-in singular-value diagnostic the
module carries.

## Running

```bash
src/tests/constrained_network_layout/run.sh
```

Objects, module files and the executable go to
`build/tests/constrained_network_layout/`, which the script creates. `FC`
selects the compiler (default `gfortran`) and `LAPACK_LIBS` the LAPACK, which
otherwise follows the same rule as the Makefile: the OpenBLAS of the
compiler's own prefix when there is one, `-llapack` otherwise. The verdict
lines are the ones `src/tests/physics_probe/assertion_report.f90` prints, and
that module is compiled from there.

Output is one line per assertion,

```
PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
```

and the exit status is nonzero if any assertion failed or any build failed.

## What the test asserts

`transported_proton_constraint_layout.f90`. The cell is the H/He block alone
(no metastable, no metals, no oxygen carriers): seven balance rows (H+, He+,
He++, H2, H2+, H3+, HeH+) and two conservation rows (H, He), at a composition
far above the trace-holdout threshold so that nothing is held out and the
cell layout alone is measured.

`ion_cell_state` carries three flags saying that a partition of the cell is
owned from outside it: `x_h2_fixed` (the lower-boundary reservoir),
`x_ox_fixed` (carrier transport of OH and H2O) and `x_hp_fixed` (the
transported ionization state, the input key `Ionization transport`). For each
of them the same statement has to hold: the species is handed its density, is
therefore not an unknown, and its own balance row is not a row of the system.
For H+ the reason is sharp: the fraction system the continuation's candidate
is judged against replaces the H+ balance row by the constraint
`x(1) = x_hp_fix` (`System_HeH_mol_metals.f90` 282-290), so a continuation
that solves the local H+ balance instead returns a composition measured
against a row it was never asked to satisfy.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `local_proton_is_an_unknown` | H+ is an unknown with `x_hp_fixed` false | 1 | exact |
| `transported_proton_is_not_an_unknown` | H+ is an unknown with `x_hp_fixed` true | 0 | exact |
| `reservoir_h2_is_not_an_unknown` | H2 is an unknown with `x_h2_fixed` true | 0 | exact |
| `transported_proton_unknown_count` | Unknowns with `x_hp_fixed` true | one fewer than with it false | exact |
| `transported_proton_reaction_row_count` | Reaction rows with `x_hp_fixed` true | one fewer than with it false | exact |
| `transported_proton_conservation_row_count` | Conservation rows with `x_hp_fixed` true | the same as with it false: pinning a species removes no element | exact |
| `local_proton_layout_is_square`, `transported_proton_layout_is_square` | reaction rows + conservation rows - unknowns | 0 | exact |

The tolerances are exact because the quantities are integer counts and flags
read as 0 or 1.

## Measured

**RED at `35d9dd5`**, three assertions:
`transported_proton_is_not_an_unknown` (1 against 0),
`transported_proton_unknown_count` (9 against 8) and
`transported_proton_reaction_row_count` (7 against 6).
`set_molecular_network_layout` read `x_h2_fixed` and `x_ox_fixed` and never
`x_hp_fixed`, so the transported proton stayed an unknown and its balance row
stayed a row.

**GREEN** with the correction of `docs/development_plan_20260905_rev3.md`
section 10.2 item 12: all eight assertions pass.

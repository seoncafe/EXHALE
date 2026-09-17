# ionization_imposed_fractions

What the assertions compare, and why each one is where the statement can
break.

## The statement

Where an ionization stage is transported, its fraction in a cell is not a
root of that cell's photoionization balance: the ionization time exceeds
the flow time, and the composition of the gas is the one the flow brought.
The balance row that would have computed the fraction is then replaced by
the carried value,

```
    fvec(1) = x(1) - x_hp_fix ,      x(1) = n(H II) / n(H nuclei) ,
```

an identity row whose Jacobian row is the identity. Every other row of the
system keeps its balance and is solved against it, which is what keeps the
remaining stages -- helium, the molecular ions, the metals, and the electron
density they all share -- consistent with the transported composition.

Rows 2 and 3 of the same systems are He II and He III per He nucleus. The
cell state states no imposed helium fraction, so nothing is written there
and the assertions cover row 1 alone.

## The routine

`impose_transported_ionization_fractions`
(`src/modules/nonlinear_system_solver/ion_residual_core.f90`) is the one
place the substitution is written. Seven systems can reach the rows it
writes, and all seven call it:

| system | unknowns reaching row 1 |
|---|---|
| `ion_system_H` | H+ alone |
| `ion_system_HeH` | H+, He+, He++ |
| `ion_system_HeH_TR` | those and the He 2^3S metastable |
| `ion_system_HeH_metals` | those and two stages of each metal element |
| `ion_system_HeH_TR_metals` | the metastable and the metals together |
| `ion_system_HeH_mol` | the molecular network above the H/He block |
| `ion_system_HeH_mol_metals` | that network with the metals above it |

Seven transcriptions of one substitution is the drift the `g2s`/`G2s`
shadowing of 2026-08-12 stands as the warning about, so the second RED below
is the one the suite exists for.

## The assertions

The cell is a warm, partly ionized, partly molecular gas with every rate the
seven residuals read set nonzero, so no row is trivially satisfied and the
balance of row 1 stands far from the constraint. The tolerances are exact:
an imposed row is an assignment, not an arithmetic result, and a row the
substitution does not touch is the same floating-point number it was.

| assertion | what it measures | tolerance |
|---|---|---|
| `routine_writes_nothing_when_unset` | largest change the routine makes to a filled residual with no fraction imposed | exact |
| `routine_writes_the_constraint_in_row_one` | `fvec(1) - (x(1) - x_hp_fix)` with the fraction imposed | exact |
| `routine_leaves_the_helium_rows_alone` | largest change to rows 2 and above | exact |
| `<system>_imposed_row_is_the_constraint` | the same difference, inside each of the seven systems | exact |
| `<system>_imposed_moves_no_other_row` | largest difference between the imposed and unimposed residual over rows 2 and above | exact |
| `<system>_unimposed_ignores_the_carried_value` | largest difference between the residuals at two different carried fractions, neither imposed | exact |
| `<system>_unimposed_row_is_a_balance` | how far row 1 stands from the constraint when nothing is imposed | strictly positive |
| `imposed_row_agrees_across_the_seven_systems` | largest difference between the imposed row 1 of any two systems | exact |

`<system>_unimposed_ignores_the_carried_value` is the one that says a cell
the transport operator does not own is solved by the rows it was solved by
before the substitution existed: the carried value is in the cell state and
is not read.

No row is FAIL by design.

## RED, MEASURED

Both were run in a copy of `src/` with the production text changed there and
nowhere else.

- The substitution removed from the shared routine: 9 assertions FAIL, one
  per system plus `routine_writes_the_constraint_in_row_one` and
  `imposed_row_agrees_across_the_seven_systems`.
- The call removed from `System_HeH_TR` alone, the routine intact: 2
  assertions FAIL, `HeH_TR_imposed_row_is_the_constraint` and
  `imposed_row_agrees_across_the_seven_systems`. This is the drift case, and
  it is caught in the one system that drifted.

## Running it

```
src/tests/ionization_imposed_fractions/run.sh
EXHALE_TEST_OBJDIR=<dir> src/tests/ionization_imposed_fractions/run.sh
```

The driver links the production sources under `src/modules`, never a copy:
`source_closure.py` follows its `use` statements. `EXHALE_TEST_OBJDIR`
selects another object directory so two builds can run side by side, and
`EXHALE_OBJDIR` names the build tree the `build_stamp.f90` is taken from.

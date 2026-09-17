# `oxygen_chemistry` -- the oxygen cycle in the regression matrix

## Why it exists

Item 10.3 of `docs/development_plan_20260905_rev3.md` records that no case of
the matrix ran the oxygen chemistry:

> The four A2/E1 test programs assert nothing (all `write`); no
> oxygen-chemistry case in the matrix -- `src/tests/a2_*`, `e1_h2`,
> `run_check.sh` 187 -- the whole oxygen path is outside the regression.

`mol_carrier`, the one case that runs the carrier transport operator, runs it
with `Oxygen chemistry` off, so `oxygen_rates.f90`, `water_photolysis.f90`,
the OH / H2O / CO carriers of `diffusive_photochemistry.f90` and the five FUV
photolysis bands had no golden of any kind. This case gives them one.

## Configuration

`input.inp` and `metals.inp` are byte copies of
`examples/18_oxygen_chemistry/`, the example written for this option:
HD 209458 b, `Molecular chemistry: True`, `Oxygen chemistry: True`, the C, N
and O abundances of `metals.inp` (2.69e-4, 6.76e-5, 4.90e-4 by number
relative to H), and the five band fluxes at the planet

    LW 3.430e+02   B1 1.379e+02   B2 (Lya) 4.809e+03
    B3 6.948e+02   B4 6.588e+05     erg cm^-2 s^-1

`Solver: Newton` is kept, as the example has it, and is inert here: the JFNK
hand-off fires at `du < 1.00E-02` and the run is stopped by `maxsteps` long
before that (`du = 9.03E+02` at step 1000), exactly as in the eight molecular
gates of `DEFAULT_CASES`. The case is a relaxation snapshot, not a converged
solution.

## Step count

`maxsteps` = 1000, chosen from a measurement rather than a guess. MEASURED
2026-09-05, single-threaded (`OMP_NUM_THREADS=1`), gfortran 16.2 build
`e779298d0b482bfbc40ae2f9aa80d60a`: **3 min 52 s wall for 1000 steps**, ending
at `final: count=1000  du= 9.0284E+02  dtu= 1.7485E-01`, Mdot 9.62
(log10 g/s). That wall time was taken while four other EXHALE runs of other
workers were on the machine, so it is an upper bound on the uncontended cost;
the requirement it was chosen against is "under 10 minutes single-threaded"
and 1000 steps meets it with the contention included. 2000 steps would sit at
about 7.7 min under the same contention, which is too close to the limit to
pin.

The oxygen network is active from the first step: the run log carries the
band fluxes above, the `Molecular base` particle count
(`q_H2(base, chem.eq. fit) = 0.831 -> ntot_bc = 0.547`) and the
self-shielding-table warning at the H2 column, so the state being pinned is
not an untouched initial condition.

## Provenance

Built 2026-09-05 for Phase 0 item 4 of
`docs/development_plan_20260905_execution.md`. Not in `DEFAULT_CASES`:
`golden/` has no entry for it and is not refreshed before the end of Phase 1.
Its reference output sits in `baseline_post170_20260905/oxygen_chemistry/`.

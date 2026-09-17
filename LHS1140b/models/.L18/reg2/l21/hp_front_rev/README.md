# `hp_front` -- the transported ionization state, restarted as it was written

## What this case is

The hot-Uranus molecular gate of `mol_base_handoff` (1 microbar base, the
Koskinen et al. 2022 handoff level, `q_H2_base = 0.84`, hydrogen and helium
only), restarted from that case's converged-at-12000-steps state and marched
100 further steps with the ionization state carried as a transported quantity:

    Molecular chemistry: True
    Molecular carrier transport: True
    Ionization transport: True

H+ is then the fifth carrier of the transport-chemistry solve
(`diffusive_photochemistry.f90`) instead of the local root of each cell's
ionization balance, and the sweep is handed the transported value.

`hp_front` is the unmodified member of the trio: its `IC/Ion_species_IC.txt`
is a byte copy of `mol_base_handoff/output/Ion_species.txt`, so the proton
profile the transport starts from is the ionization front the marching run
built. (`IC/Hydro_ioniz_IC.txt` differs from the same case's re-run only in
its `# provenance:` comment, which names the build the copy was taken from;
MEASURED 2026-09-05 the numeric data of the two is byte-identical.) `hp_zero_seed` and `hp_trace_seed` are the same configuration with the
H II column of that file replaced (see their READMEs).

## Why the initial condition is a restart, and not a key

The carrier constraint is imposed on the equilibrium sweep only when the
background is ready or an initial condition was loaded
(`src/modules/radiation/ionization_equilibrium.f90` 1203-1211):

    if (thereis_mol .and. carrier_transport                   &
        .and. (bg_ready .or. do_load_IC)) then

On a cold start the very first sweep is what INITIALIZES the carriers, from
the local equilibrium of the seeded state, so a proton seed written into a
cold start is overwritten before anything reads it. A restart is therefore
the only way a chosen H II profile reaches the solver, which is why all three
cases carry an IC pair rather than an input key.

## Where the IC pair lives

`load_IC.f90` reads `output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`
(lines 119, 161, 181), and `run_check.sh` clears `output/*.txt` before every
run. The pair is therefore kept in `IC/` and restored into `output/` by
`check_case` immediately after that clear. Do not put it in `output/`: it
would survive exactly one invocation.

## Configuration, and the two keys that are refused

`input.inp` is `mol_base_handoff/input.inp` with `Load IC? True`, the two
carrier keys above, and **`Solver: Newton` removed**. `input_read.f90`
1870-1881 refuses `Ionization transport: True` together with `Solver: Newton`
(the steady solve would hold the composition at its own local root and undo
on the last iteration exactly the departure from local equilibrium this
option computes) and 1853-1861 refuses it together with
`Coupled carrier solve: True`. Neither key is present.

`base.inp` is a byte copy of `mol_base_handoff/base.inp`.

## Step count

`maxsteps` = 100. MEASURED 2026-09-05, single-threaded, gfortran 16.2 build
`e779298d0b482bfbc40ae2f9aa80d60a`: 5.7 s wall, `final: count=100
du= 1.2828E+00`, Mdot 10.51 (log10 g/s). The case is a relaxation snapshot,
not a converged solution, and nothing here claims otherwise.

## What is guarded, and what is not

GUARDED (bitwise, like any matrix case): that the transported-proton path
runs at all, and that its 100-step state does not move. MEASURED 2026-09-05,
the three cases are mutually distinct states, so the seed is not being
discarded on the way in: against `hp_zero_seed`'s `Ion_species.txt` the worst
relative difference is 1.0 (a column that is zero in one and not in the
other), and `hp_zero_seed` against `hp_trace_seed` likewise.

NOT GUARDED: the assertion the plan asks for -- the carrier state read
immediately before `carrier_source` and before the Jacobian assembly. There is
no diagnostic hook that dumps the carrier fractions at that point:
`grep EXHALE_DUMP src/` finds `EXHALE_DUMP_IC` only, which writes the loaded
state before the first sweep, not the state the carrier solve is handed. The
hook belongs at the entry of `solve_carriers`
(`src/modules/lower_atmosphere/diffusive_photochemistry.f90` 1717) for the
state the solve starts from, and inside `carrier_residual` beside the
`carrier_source` call (line 1596) for the state each row is built from. Adding
it is production-source work and was out of this task's scope.

## Provenance

Built 2026-09-05 for Phase 0 item 4 of
`docs/development_plan_20260905_execution.md`, from the `mol_base_handoff`
outputs then in the tree (`git=6d07d48afd41 tree=dirty`,
`ck_input=2432225406772013242`, `ck_base=579698832458720949`). Not in
`DEFAULT_CASES`: `golden/` has no entry for it and is not refreshed before the
end of Phase 1. Its reference output sits in
`baseline_post170_20260905/hp_front/`.

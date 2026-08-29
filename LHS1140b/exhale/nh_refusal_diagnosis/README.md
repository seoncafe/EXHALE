# Why the profile adapter refuses reservoir He/H = 8.5-8.7 on LHS 1140 b

Measurement record, 2026-08-29, this directory.  Nothing outside it was
written; no code, golden or stored arm was touched.  Every adapter run here
is chemistry + climate only (`photochem_to_lower_profile.py`), with no wind
step and no closure iteration.

## What was measured

1. `scan_reservoir.sh` -> `scan/`, `scan_table.txt`
   The adapter, in the closure ladder's configuration verbatim, at 40
   reservoir values from He/H = 1 to 15, with the conservation check
   disarmed on the command line (`--abundance-tol 1.0`) so the profile is
   written and the departure can be read.  The default tolerance in
   `lower_profile_schema.py` was not changed.

2. `probe_initial_vs_steady.py` -> `probe/`, `probe_heh8p6_summary.txt`
   The same adapter call at He/H = 8.6, with Photochem's gas-giant class
   subclassed to snapshot the composition (a) as chemical equilibrium plus
   the quench overwrite hands it over on the climate grid, (b) as the
   photochemical model's initial state, (c) converged.

3. `equilibrium_closure.py` -> `equilibrium_closure_heh.txt`,
   `equilibrium_closure_temperature.txt`
   Photochem's equilibrium solver (`equilibrate.ChemEquiAnalysis`) alone, at
   the deepest level of the LHS 1140 b climate column: its elemental closure
   against the abundances it was given.

4. `equilibrium_deep_composition.py` -> `equilibrium_deep_composition.txt`
   The same solver, reporting the deep mixing ratios themselves, to separate
   "a different chemical branch" from "a larger solver residual".

## What the measurements say

**The departure is the residual of Photochem's equilibrium solver, not a
nitrogen leak and not a change of chemical branch.**

- The deepest-level departure is already present in the equilibrium
  composition, before any photochemistry: (a), (b) and (c) agree to every
  digit printed (13 significant figures), and the converged profile's deep
  `q_NH3` is 8.995632484e-06 against the initial 8.995632e-06.  Photochem
  pins that cell: `_initialize_atmosphere` sets `bc_type='press'` for every
  gas species at the model bottom, so the deepest level is a Dirichlet
  boundary held at the equilibrium composition and the photochemical solve
  cannot move it.
- `composition_at_metallicity` in `photochem/extensions/gasgiants.py` takes
  the equilibrium result whether or not it converged -- its own comment is
  `# Do not enforce convergence.` -- and `solve` reports `converged = True`
  at every point measured here while its elemental closure on N is off by up
  to 2.0e-4.
- The discrepancy is NH3-stoichiometric: one excess N per three excess H.
  `dev_N / dev_He = 4072` in every one of the 38 runs, against
  1/(3 N/H) = 1/(3 x 8.18465e-5) = 4073.  He, C and O share one common
  departure (they carry the H error alone); N carries the H error plus a
  relative error in the equilibrium NH3 mole fraction, and essentially all
  the nitrogen at the deep level is NH3 (N2 is 3e-10, HCN 4e-28).
  The check therefore tests the equilibrium solver's NH3 mole fraction 4072
  times more tightly than it tests anything else.
- The deep composition is smooth and monotone straight through the refused
  band (`equilibrium_deep_composition.txt`): H2, H2O, CH4, NH3, N2, He all
  vary by less than 1% across 8.45 -> 8.75 with no jump.  What jumps is the
  residual: 1.99e-4 at He/H = 8.71, 5.2e-6 at 8.72, a factor 38 across a
  0.1% change in composition.  Same at 9.15 -> 9.2 (6.6e-5 -> 4.5e-8) and
  10.25 -> 10.3 (2.0e-4 -> 3.5e-8).
- 8.5-8.7 is not a special band.  At tolerance 1e-4 the scan also refuses
  He/H = 6.0 (1.31e-4) and 10.20-10.25 (1.5-2.0e-4), and 10.1 sits at
  6.7e-5, just under.
- It is not a property of this planet either.  At solar composition
  (He/H = 0.0969) the same solver's N closure is 1.50e-4 at 800 K and
  1.14e-4 at 1600 K -- both above the tolerance.
- The adapter's own carrier accounting at the deepest level is complete:
  the profile's diagnostic columns (`q_NH3` + 2 `q_N2` + `q_HCN`) account
  for 100.00% of `X_N` there.  (At the matching level they account for
  93.9%, the rest being N, NH, NH2 and the like; that is a property of the
  diagnostic columns, not of the `X_N` column, and it is not what the check
  reads.)

## Unrelated failure seen while scanning

`clima`'s `surface_temperature_bg_gas` root solve fails outright at
He/H = 0.0969, 8.1 and 15 in this configuration (`hybrd1 root solve
failed`), so no column is produced at all.  That is a different failure from
the refusal and was not investigated here.

# Recipe: LHS 1140 b with the Photochem closure

The escape wind of LHS 1140 b joined at 1 microbar to a Photochem/clima
column below it, iterated to the fixed point at which the elemental fluxes
the chemistry is handed equal those the wind carries, and the He I 10830
equivalent width of the result. This directory holds no input files: the
recipe points at the closure configuration of record and its scripts.
Labels: MEASURED = run for this recipe on 2026-10-10 (lart4, the `EXHALE.x`
of the tree built 2026-10-10 15:02); READ = from a stored log or the source
text named.

## The pieces

| Piece | Path |
|---|---|
| closure driver | `src/utils/element_flux_closure.py` (`--print-config-template` prints a commented configuration) |
| chemistry/climate adapter | `src/utils/photochem_to_lower_profile.py` (needs Photochem; `README_photochem.md`) |
| configuration of record, He/H 2.09 | `LHS1140b/models/.P2/ew_fit/phase3/energy_coupling_nlte_structural_a76da494/`: `closure.json`, `input_template.inp`, `launch.sh` (the command), `closure.log`, `closure_history.txt`, iterations `k00` to `k05` |
| SED | `LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt` (at the planet) |
| He 10830 EW of a state | `LHS1140b/models/.P1/fiducial_conduction/ew_of_state.sh`, `LHS1140b/he10830_equivalent_width.py` |

The stored iterations were made with `LHS1140b/models/.P2/EXHALE_a76da494.x`,
before update log stage 3 sections 125-134. They are a record, not a current
result: the k05 state evaluated with the current binary (`Restart intent:
stationary evaluate`, MEASURED) exits 2, not certified, its largest row the
energy row at 0.63 of its scale at cell 170 (the H3+ cooling function changed
in section 132, among others).

## What the wind input sets (`input_template.inp`)

`Domain mode: Spherical`, `Outer radius [R_p]: 30.0` (outermost cell at
Mach 1.16 on k05, MEASURED); `2D approximate method: Mdot` (dayside dilution
1, the substellar beam: a physical choice of the model, report it with every
result); `He/H number ratio: 2.092388551`; `He_diffusion`, `Conduction`,
`Molecular chemistry`, `Molecular carrier transport`, `Ionization transport`,
`Well balanced`: all `True`; `Lower atmosphere profile:
lower_atmosphere_profile.dat` (written by the adapter in each iteration).
The configuration adds `Restart intent: stationary`, `Solver: Newton`,
`Secondary_ionization: Immediate` (the same as the default `True` since
2026-10-09) and the residual tolerance `Resid tol: 4.0e-6`; `exhale_env`
holds the outer-alternation settings (`EXHALE_DIFF_OMEGA 0.5`,
`EXHALE_CARRIER_TRUST_HOLD 0.04`, `EXHALE_OUTER_PASSES 60`, ...).

## What the chemistry is handed (adapter defaults since 2026-10-10)

- `--face-geometry wind` (default): Photochem and the climate column take the
  stellar energy the wind's input states, xi_w F with xi_w the wind's dayside
  dilution. The closure adds `--face-geometry wind --wind-dilution xi_w` to
  every chemistry call unless `adapter_args` names a geometry, with xi_w from
  the template's `2D approximate method` (`Rate/4` 0.25, `Rate/2` 0.5,
  otherwise 1). `per_model` reproduces runs made before that date (the stored
  k00-k05 among them). Check the wind's value in `EXHALE_setup.out`:
  `Dayside dilution applied to every stellar beam:  1.000`.
- `--eddy-heat-transport off` (default): the chemistry's `--kzz-const 1.0e9`
  mixes tracers and carries no heat. As a heat diffusivity no reference
  climate exists above K_zz of about 1.4e6 cm^2 s^-1 at xi_w = 1 (update log
  stage 3 sections 130 and 133); `on` is kept for sensitivity tests.
- The hydrogen flux handed to the chemistry is the hydrogen-NUCLEUS flux
  F_Hnuc (`lower_profile_F_Hnuc_*`, `steady_F_Hnuc_*` in
  `EXHALE_resolved.out`), not the component-1 mass flux F_H, which includes
  the trace metals and the C and O bound in molecules (section 134). MEASURED
  on k05 re-evaluated: steady F_Hnuc/F_H = 1.2945e7/1.3003e7 = 0.9956.

## Command

To run a new closure, copy the configuration directory's `closure.json` and
`input_template.inp` into a new case directory, set `exhale_bin` to the
binary to use and `input_template` to the copied template, and run, as
`launch.sh` does:

```sh
/opt/miniconda3/bin/python3 src/utils/element_flux_closure.py <case dir> \
    --config <case dir>/closure.json --seed <a converged output/ directory> \
    --phi0-H <g/s> --phi0-He <g/s> --fcond0 <erg cm^-2 s^-1> \
    --omega 0.5 --tol 0.05 --kmax 16
```

`--resume` continues from the last completed iteration. Each iteration `k<n>/`
holds `adapter.log` (chemistry), `lower_atmosphere_profile.dat`, `base.inp`,
`run.log` (wind) and `output/`; `closure.log` and `closure_history.txt` record
the fluxes, residuals, F_cond and T_match of every iteration.

Then the line, width-matched (project `CLAUDE.md`, "Line width of a model
transit line against an observation"):

```sh
EW_NONTHERMAL_KERNEL=width_match \
  LHS1140b/models/.P1/fiducial_conduction/ew_of_state.sh <case dir>/k<n> <new work dir> <binary>
```

`EW.txt` gives `EW=` (unbroadened window EW) and, for each product, `fwhm_kms`
(the Gaussian velocity FWHM f at which the three-Gaussian red-pair FWHM is
the observed 0.841 A) and `EW_kernel` (the width-matched window EW,
10832.60-10834.20 A vacuum, the number compared with the observed EW). The
kernel stands for a velocity field the 1-D model does not contain; report f
with every EW.

## Run time

READ (`closure.log` of the configuration of record, 12 threads,
binary a76da494): a chemistry call 4.4-6.5 min; a wind solve 5-35 min;
iterations k00-k03 took 77 min; the closure converged at k = 5 (largest
residual 0.0366 <= 0.05, T_match moved 1.43 K < 6 K). With the current
binary not measured. MEASURED: the evaluate of the stored k05 state with the
current binary, 6.9 s on 2 threads.

## Acceptance

- Wind exit status: 0 certified or no claim; 2 converged and written but not
  certified, a flag the closure carries (`wind: converged but NOT certified
  (exit 2)` in `closure.log`); 2 without the closing line `(certification) the
  state was written in full; the run exits with status 2 because it is not
  certified.` is a gfortran runtime error and stops the closure; 3 the
  stationary solve failed, which stops it. `ew_of_state.sh` treats an
  evaluate exit 2 the same way (products made, labeled uncertified).
- A closure is converged when the undamped flux residual is below `--tol`
  (and, with `energy_coupling`, T_match moves less than `T_match_tol` and the
  conducted-flux residual is below `--tol`).
- The result is accepted by the precision of its observable
  (`md/acceptance_by_observable_precision_20261009.md`): continuing changes
  the width-matched EW by less than 0.002 %A, the mass and energy terms
  balance over 1.5-12 R_p, and T, v, n_e, n(2^3S) change by less than 1%
  there. The certification is a diagnostic.
- Limitations to state with every result (plan rev4, step F): the exobase
  lies inside the sonic point (the continuum treatment beyond it), the 1-D
  line width needs the kernel above, K_zz has no observational basis (state
  it as a sensitivity, 1e8 to 1e10), and the column takes the substellar
  geometry.

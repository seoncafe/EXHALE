# Recipe: an atomic hot Jupiter, HD 209458 b

A complete escape wind of an atomic H/He atmosphere (H, He, He 2^3S, the
Ly-alpha escape-probability field), from a cold start to a Newton-finished
stationary state. Labels: MEASURED = run for this recipe on 2026-10-10 (lart4,
the `EXHALE.x` of the tree built 2026-10-10 15:02); READ = taken from a stored
log or from the source text named.

## Input files

`input.inp` (this directory) is the planet folder `HD209458b/input.inp` with
three changes: `Load IC? False` (cold start instead of a stored state), the
default switch of `Solver: Newton`, and no `Base grid [dr,cells]` line (the
grid of the current default). No `metals.inp`, `base.inp` or `opacity.inp`:
each is switched on by its presence, so this run is metals-off with the base
state of `input.inp`. To add the trace metals, copy `HD209458b/metals.inp`
here.

Keys that matter (each checked in `src/modules/files_IO/input_read.f90` and
`src/modules/init/parameters.f90`):

| Key | Value here | What it does |
|---|---|---|
| `Spectrum type` / `Spectrum file` | `Load from file..` / `../../inputdata/sed/hd209458b_solar_whi2008.txt` | stellar SED at the orbit; the path is relative to the run directory. The two `Log10 of ... luminosity` lines are still required but the file's integral replaces them |
| `2D approximate method` | `Rate/2 + Mdot/2` | dayside dilution 0.5 of every stellar beam (printed in `EXHALE_setup.out` as `Dayside dilution applied to every stellar beam`); the reported Mdot is corrected by log10(2) |
| `Reconstruction scheme` + `du_th [PLM,WENO3]` | `PLM+WENO3`, `0.5 1.0e-3` | two-stage marching: PLM until the cell-centered du < 0.5, then WENO3. The second number only makes the run two-stage (it must be below the first); it no longer stops the run |
| `Solver` | `Newton` | JFNK Newton finish. The hand-off fires when the FACE mass-flux spread falls below the switch, default 0.1 (`Solver: Newton <value>` changes it); before 2026-10-10 the value was a du threshold |
| `Secondary_ionization` | (omitted) | default `True`: the photoelectron coupling is applied from step 0. `Staged` gives the earlier behavior (switched on once the marching has settled), `False` turns it off |
| `Include He23S?` | `True` | the metastable He 2^3S level (He I 10830) |
| `Jlya escape-prob`, `Stellar Lya flux`, `Lya stellar halfwidth`, `Lya stellar boost`, `Stellar Teff`, `Stellar radius` | as shown | Ly-alpha field by the escape-probability closure and the H(n=2) Balmer continuum |
| `Outer radius [R_p]` | (omitted) | the domain ends at the Roche L1 radius, 4.048 R_p here (MEASURED, `outer_radius_Rp` in `EXHALE_resolved.out`) |

## Command

```sh
cd examples/recipe_hd209458b_atomic
OMP_NUM_THREADS=8 ../../EXHALE.x > run.log 2>&1; echo "exit status $?"
```

The binary reads `./input.inp` and writes `./output/`.

## Run time

- Not measured end to end with the current binary. MEASURED: 5 steps
  (`EXHALE_MAXSTEPS=5`, 2 threads) take 12.8 s including startup and the
  post-processing pass; that is not a step rate.
- READ, earlier binaries (different stop and hand-off rules, so only an order
  of magnitude): `benchmarks/hd209/run_20260815_kb.log`, the same planet with
  element diffusion, 60 975 steps and 23.7 min on 8 threads;
  `HD209458b/run_20260810.log`, a cold start of the planet folder under the
  old du hand-off, 924 705 steps and 5.2 h on 16 threads.

## How the marching stops, and the hand-off

- The marching stop is the face mass-flux spread below 4e-4
  (`march_flux_spread_th`) AND the largest relative state change over one
  monitor window of 500 steps below 1e-4 (`march_change_th`). Every monitor
  step prints a `[stop]` line with both numbers. The cell-centered `du` is a
  printed diagnostic (`[diag]` line) and drives only the PLM -> WENO3 switch.
  `EXHALE_MARCH_FLUX_TH`, `EXHALE_MARCH_CHANGE_TH` and `EXHALE_MONITOR_STEPS`
  replace the two thresholds and the window for measurement.
- With `Solver: Newton` the JFNK hand-off happens in the WENO3 stage at a
  monitor step where the face spread is below the switch (0.1), or where the
  spread has settled within 5 times the switch.

## Outputs

- `output/Hydro_ioniz.txt`: r [R_p], rho, v, p, T, heating, cooling;
  `output/Ion_species.txt`: the species densities; the `*_adv.txt` pair: the
  advection-corrected profiles that `EXHALE_transit.py` reads;
  `output/Heating_breakdown.txt`, `output/Cooling_breakdown.txt`: channels;
  `output/Excited_H.txt`: H(n=2); `output/conserved_state.txt`: the exact
  conserved state for a restart. Every file has a `# columns` header; load
  them with `examples/exhale_io.py` (`exhale_io.Run('output/')`).
- `EXHALE_setup.out`: the configuration as the run understood it, and at the
  end the Mdot and the Mach number of the outermost cell.
  `EXHALE_resolved.out`: the same configuration in machine-readable form.
- `run.log` ends with the `----- Results -----` block:
  `Log10 of steady-state Mdot` and `Mach number of the outermost cell`. Below
  Mach 1.1 a NOTE follows: the outer face is subsonic or near sonic, the outer
  ghost condition then supplies the incoming characteristic and can set Mdot.
  For a quantitative Mdot put the outer face past the sonic point with
  `Outer radius [R_p]: <value>` (measured on the WASP-121 b regression case
  `wasp_full_newton`: L1 at Mach 0.95 against 2 R_p at Mach 1.5 moves rho and v by up to
  15-17% near 1.31 R_p and log10 Mdot by 0.01, READ,
  `md/outer_boundary_subsonic_review_20261009.md` section 5).
- Transit spectra: `MPLBACKEND=Agg python3 ../../EXHALE_transit.py` in this
  directory writes `tpm_*.txt` (He I 10830, Ly-alpha, H-alpha, H-beta).

## Acceptance

Exit status (`src/modules/time_step/certification.f90`):

| Status | Meaning |
|---|---|
| 0 | the state is certified, or the run made no stationary claim (the marching stop, a step cap, a stall, a NaN): a refused state is written `cert_reason=no_stationary_claim` and the log says `relaxation snapshot ... no stationary claim` |
| 2 | the run claimed a stationary state (a Newton finish, or the residual gate met while marching; the marching stop is not a claim, update log stage 3 section 135) and the certification refused it; the state is written in full. A flag, not a refusal. A gfortran runtime error also exits 2: the binary's own status 2 is the one whose log ends with `(certification) the state was written in full; the run exits with status 2 because it is not certified.` |
| 3 | a stationary solve failed (info /= 0 on the stationary restart or the direct steady route); the state written is the solver's last iterate with `cert_reason=stationary_solve_failed` |
| 1 | a Fortran `error stop` (input refused, a state that cannot be loaded, ...) |

A result is accepted by the precision of the observable it is used for, not
by the certification (`md/acceptance_by_observable_precision_20261009.md`):
continuing the solve changes the observable by less than its budget (log10
Mdot: 0.3% in Mdot; profiles in the line-forming region: 1%), the mass and
energy terms balance over the region that forms the observable, and the
profiles are stable. The certification block in `run.log` (its failing rows
and their cells) is a diagnostic: a row far above its scale (a factor 100) is
looked at as a possible defect before the state is used.

Two cautions for this planet: quote Mdot only from the Newton-finished state,
never from a marching snapshot, and check the Mach line (the L1 face lies at
4.05 R_p; the stored state of the planet folder, written 2026-08-19 on an
older grid, has its last physical cell at 4.15 R_p and Mach 1.29, computed
here with gamma = 5/3 from its `output/Hydro_ioniz.txt`).

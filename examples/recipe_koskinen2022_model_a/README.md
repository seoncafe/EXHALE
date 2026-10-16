# Recipe: a molecular base, the Koskinen et al. (2022) Model A setup

The hot Uranus of Koskinen et al. (2022, ApJ 929, 52) Model A with a
1 microbar molecular base: H2/H2+/H3+/HeH+ chemistry, the H2 carrier and the
H/He ionization stages transported, heat conduction, no Lyman-Werner band,
H/He only. The input files are those of the matched comparison run
`benchmarks/koskinen2022_model_a/matched_hnu_minus_I/` (its condition table and
overlaid profiles: `docs/koskinen2022_model_a_comparison.pdf`); only the
spectrum path is made relative. Labels: MEASURED = run for this recipe on
2026-10-10 (lart4, the `EXHALE.x` of the tree built 2026-10-10 15:02); READ =
from a stored log or the source text named.

## Input files

- `input.inp`: the matched run's input.
- `base.inp`: the base state handed over, `T_base 1140`, `p_base 1e-6` bar
  (this fixes the base level; `input.inp` carries no density key),
  `q_H2_base 0.84` (Koskinen et al. 2022 section 3.2), He/H 0.0793.
- `output/Hydro_ioniz_IC.txt`, `output/Ion_species_IC.txt`: the starting
  state, made in step 1 below (not stored here): the matched run's final
  state, so the run continues where the published comparison stood (that run
  was itself started from an earlier state interpolated onto its grid, the
  benchmark's README).

Keys that matter (checked in `input_read.f90` / `parameters.f90`):

| Key | Value | What it does |
|---|---|---|
| `Domain mode` / `Outer radius [R_p]` | `Spherical` / `7.24` | planetary potential only, domain to Model A's top (9.6 R_p in base radii); the outermost cell is at Mach 1.26 (MEASURED, 5-step run from the mapped state) |
| `2D approximate method` | `Rate/4 + Mdot` | dayside dilution 0.25 (`EXHALE_setup.out`) |
| `Spectrum file` | `../../benchmarks/koskinen2022_model_a/solar_ribas2005_bands_0p05au.txt` | the Ribas et al. (2005) band spectrum at 0.05 au, `E_high` 2 keV |
| `Molecular chemistry` | `True` | the molecular network in the coupled ionization system, and the molecular base particle count |
| `Molecular carrier transport` | `True` | H2 transported by advection-diffusion |
| `Ionization transport` | `True` | x(H II), x(He II), x(He III) transported |
| `Conduction` | `True` | heat conduction, kappa of electrons, H I, He I and H2 |
| `Atomic rate set` | `Koskinen2022` | their Table 1 R1-R4 rates |
| `Stellar LW flux [erg/cm2/s]` | `0.0` | no H2 photodissociation (Model A has none) |
| `Secondary_ionization`, `He_rec_coupling`, `Molecular reaction heat` | `False` | the photoelectron accounting of Model A, all of h nu - I heats (the condition table of the comparison memo). The code default of `Secondary_ionization` is `True`, on from step 0 |
| `Include He23S?` | `False` | no metastable helium |
| `du_th [PLM,WENO3]` | `1.0e9 1.0e-3` | the PLM stage ends at once (WENO3 from the first step); the second number no longer stops the run |
| `Stall [tol,N]` | `0.0 2000` | `N` is the window of the stall net (settled state while the face flux has stopped falling); `tol` acts only in the PLM stage, which is empty here |
| `Load IC?` | `True` | start from `output/*_IC.txt` |

## Commands

From this directory:

```sh
# 1. the starting state: the matched run's state mapped onto this grid
#    (load_IC refuses a state whose cell centers differ by more than 1e-10;
#    the stored state predates the grid change of 2026-09-27)
python3 ../../src/utils/map_state_to_grid.py \
    ../../benchmarks/koskinen2022_model_a/matched_hnu_minus_I/output \
    ../../backup/regression/carrier_model_a_newton/IC/Hydro_ioniz_IC.txt \
    output --ic
# 2. the run, capped as the published comparison was
EXHALE_MAXSTEPS=200000 OMP_NUM_THREADS=16 ../../EXHALE.x > run.log 2>&1
echo "exit status $?"
```

The second argument of step 1 is only a grid: any `Hydro_ioniz*.txt` written
by a run of this `input.inp` serves (the regression file named is one, and
`backup/` exists only in this working copy).

## Run time

- READ: the matched run, 200 000 steps on 16 threads, 19 min
  (`benchmarks/koskinen2022_model_a/matched_hnu_minus_I/run.log`, 2026-09-05
  binary).
- MEASURED: step 1 and a 5-step run (`EXHALE_MAXSTEPS=5`, 2 threads) from the
  mapped state end in 0.7 s with exit status 0; this confirms that every key
  above parses and the mapped state loads, not a step rate.

## Outputs and what to compare

- `output/Hydro_ioniz.txt`, `output/Ion_species.txt` and the `*_adv.txt`
  pair; `output/Heating_breakdown.txt`, `output/Cooling_breakdown.txt` (the
  H3+ channel is the MiZATeP cooling function of Mizus et al. 2017 since
  2026-10-10, update log stage 3 section 132; the stored benchmark states used
  the Miller et al. 2013 fit).
- `EXHALE_setup.out`, `EXHALE_resolved.out`: the configuration in force.
- `run.log`, `----- Results -----`: log10 Mdot and the Mach number of the
  outermost cell. The stop line says why the marching ended: the face
  mass-flux spread < 4e-4 AND the state change over a 500-step window < 1e-4
  (each monitor step prints a `[stop]` line with both), the stall net, or the
  step cap.
- Compare with Model A by `python3 ../../docs/k22_model_a_figures.py output`
  (its first argument is the output directory, default the stored matched
  run; figures go to `docs/figures/k22_model_a/` unless `K22_FIG_DIR` is
  set), against the digitized
  `model_a_fig8_digitized.txt` and `model_a_fig9_digitized.txt`. Overlay the
  whole profiles, not one value.

## Acceptance

Exit status: 0 certified or no stationary claim (a run ending on
`EXHALE_MAXSTEPS` makes none, so the published-comparison run exits 0 and its
log says `relaxation snapshot`); 2 a claimed stationary state that the
certification refuses, written in full, a flag (a gfortran runtime error also
exits 2; the binary's own 2 has the closing line `(certification) the state was
written in full; the run exits with status 2 because it is not certified.`);
3 a stationary solve that failed; 1 an `error stop`.

The state is a marching end point, not a Newton root (the benchmark's
README). It is accepted for a comparison with Model A by the stability of the
compared profiles under continued marching (1% in the region compared) and
by its term balance, `md/acceptance_by_observable_precision_20261009.md`; the
certification block is reported as a diagnostic.

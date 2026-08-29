# Why `clima`'s background-pressure bracket fails at He/H = 9.4 and 9.5

Diagnosis of `TO_BE_DONE.md` item (M), measured 2026-08-29.  Every number
below was produced by the scripts in this directory against the two installed
builds; nothing in `photochem/` or in either environment was modified.

- corrected build: `EXHALE_v1.00/env/photochem/bin/python`, photochem 0.9.0,
  Clima 0.7.5, carrying `clima-bounded-scaled-solvers.patch`
- reference build: `/home/kiseon/.conda/envs/photochem_cmp/bin/python`,
  photochem 0.8.4

Both are driven with the LHS 1140 b closure ladder's configuration: 20 bar
deep boundary, TOA 1e-2 dyn/cm^2, 60 layers, `solve_for_T_trop` on,
T_deep guess 400 K, T_trop guess 120 K, the GJ 1132 SED at the planet, and
the Lodders elemental vector with He/H overridden.  `flux_photochem.txt`
here is byte-identical to the one the ladder wrote.

## Finding

The bracketing scan is not the cause.  It fails because the outer solve
hands it a surface temperature of 6.5e8 K, at which `make_profile` refuses
every one of the 49 scan points with "Failed to compute heat capacity" --
the thermodynamic polynomials stop at 6000 K (`T_ceiling.py`).  With no
valid point in the interval there is no sign change to find, so the scan
reports that it could not bracket.

The outer solve arrives at 6.5e8 K because its forward-difference Jacobian
is numerical noise.  MINPACK's default step is `sqrt(eps)*|x|`, which at
T_surf = 400 K, T_trop = 120 K is 2.4e-5 K.  Over that step the flux
residual moves by about the same amount as its own jitter, ~3e-8 in
normalized units (`jacobian_noise.py`, `noise_source.py`), so the Jacobian
entries come out with the wrong magnitude and often the wrong sign
(`first_step.py`).  The cross term df1/dx2 is tens of times the diagonal
df1/dx1, so an error in the small, noisy diagonal is amplified into a step
of several decades in log10 T_surf.

The 0.8.4 build has the same defect (`first_step_old.py`): at He/H = 9.37
its noisy Jacobian sends the first trial to T_surf = 0.10 K, which is
exactly the "T_surf is less than T_trop" refusal that build reports there.
The patch changed the parameterization and the scaling, which reshuffled
which compositions are unlucky; it did not create the trap.

The root itself is in the interval and always was.  At He/H = 9.4 the 0.8.4
solution has log10(P_bg) = 7.2785, inside the scanned
[-4.6990, 7.3010] (`old_build.py solve`).  The residual over that interval
is monotone with exactly one sign change, in the last scan step
(`scan_residual.py`); the grid cannot step over it.

## Failure sets, measured

He/H 9.30-9.60 in steps of 0.01, 31 points
(`sweep_new_9p30_9p60.txt`, `sweep_old_9p30_9p60.txt`):

| build | failures |
|---|---|
| 0.9.0 + patch | 9.40, 9.50 |
| 0.8.4 | 9.37, 9.45 |

He/H 8.00-11.00 in steps of 0.05, 61 points
(`sweep_new_wide.txt`, `sweep_old_wide.txt`):

| build | failures |
|---|---|
| 0.9.0 + patch | 9.40, 9.50, 10.45 |
| 0.8.4 | 8.10, 9.45, 9.85, 10.15 |

Isolated points in both, at the same rate.  The solution is a smooth
function of composition across them -- T_deep falls monotonically from
422.72 K at 9.30 to 420.61 K at 9.60 -- so nothing about the composition
distinguishes a failing point.  `tguess_scan.py` shows the same: at both
9.4 and 9.5 every deep-temperature guess tried except exactly 400 K solves,
and each returns the 0.8.4 answer to five decimals.

## Repair, tested in replication

`fix_variants.py` drives the patched residual through the same MINPACK
`hybrd` from Python, reproducing the installed build's converged deep
temperature to 5e-5 K where it solves.  Over the 61-point grid:

| variant | solved | failed or stalled |
|---|---|---|
| as patched | 54 | 9.40, 9.50, 9.80, 10.05, 10.25, 10.70, 10.80 |
| initial trust-region factor 1 | 54 | the same seven |
| trial T_surf confined to (T_trop, 6000 K) | 55 | 8.50, 9.55, 9.80, 9.85, 10.20, 10.35 |
| `epsfcn = 1e-4` | **61** | none |

`epsfcn` is not the difference step itself: MINPACK's `fdjac1` forms
`eps = sqrt(max(epsfcn, epsmch))` and steps `h = eps*|x|`, so `epsfcn = 1e-4`
is a relative step of `1e-2` in the log variable -- 16.2 K at the starting
point, against 2.4e-5 K at the default. (Corrected 2026-08-29: an earlier
version of this file read the value as the step and quoted 0.064 K.) The
Jacobian it produces is smooth and nearly composition-independent,
`[[-1.21, -0.55], [5.0e-4, -3.557]]` across the whole 9.30-9.55 range,
against entries that vary from +52 to -42 at the default step.

The same step applied to 0.8.4's own parameterization also solves all 61
points (`old_epsfcn.py`, `variant_old_epsfcn.txt`), so the repair does not
depend on the parameterization the patch changed.  Deep temperatures agree
with the 0.8.4 build to 1.5e-4 K wherever that build solves, and the three
compositions the corrected build refuses come out at 422.00813 (9.40),
421.30374 (9.50) and 415.15428 (10.45) K -- on the smooth trend, and equal
to what 0.8.4 returns at 9.40 and 9.50.

In `clima` this means replacing `hybrd1` with `hybrd` in
`AdiabatClimate_surface_temperature_bg_gas` (and in `surface_temperature`,
which has the same structure) so that `epsfcn` can be set; `hybrd1` fixes it
at machine epsilon.  Confining the trial temperature to the range of the
thermodynamic data is worth doing as well, but on its own it does not fix
the solve -- it only converts one failure mode into another, because the
tropopause variable is still unbounded.  Not implemented here.

## Scripts

| file | what it measures |
|---|---|
| `bracket_lib.py` | the shared climate object, and the patch's scan transcribed |
| `reproduce.py` | pass/fail of `surface_temperature_bg_gas` on the corrected build |
| `scan_residual.py` | the 49-point residual table inside one bracketing call |
| `bracket_vs_T.py` | bracket survival against trial surface temperature |
| `T_ceiling.py` | the temperature above which `make_profile` refuses |
| `outer_trace.py` | the trial (T_surf, T_trop) pairs the outer solve visits |
| `old_build.py` | 0.8.4 solutions, and whether their P_bg is inside the interval |
| `first_step.py`, `first_step_old.py` | the first Newton step at two difference steps |
| `jacobian_noise.py`, `noise_source.py` | the residual's jitter against step size |
| `tguess_scan.py` | sensitivity to the deep-temperature guess |
| `fix_variants.py`, `old_epsfcn.py` | the candidate repairs over the composition grid |

The `work_*` directories the scripts create hold only regenerated climate
species and settings files and were removed.

## Not established

The source of the ~4e-7 relative jitter in ISR - OLR was not isolated to a
particular part of the radiative transfer; the background pressure the inner
solve returns is bit-identical across the steps that produce it
(`noise_source.py`), so it is not the inner tolerance.  Whether the same trap
is reachable at other planets, other K_zz or other deep pressures was not
tested -- only LHS 1140 b's configuration was scanned.

# LHS 1140b with EXHALE — the He/H scan

Six compositions of the same planet, with the parameters of
`../system_parameters.md` and the stellar spectrum of `../sed/`, so that the
comparison against `../pwinds_oracle/` and against the 2024 measurement in
`../observed/` is one planet, one spectrum, one composition grid.

| directory | He/H | H:He |
|---|---|---|
| `solar/` | 0.0833 | 12 |
| `heh1/` | 1 | 1 |
| `heh10/` | 10 | 0.1 |
| `heh100/` | 100 | 0.01 |
| `heh1000/` | 1000 | 0.001 |
| `heh10000/` | 10000 | 0.0001 |

The p-winds scan covers the same set (`../pwinds_oracle/scan_hhe.txt`, listed
by H:He), so the two models can be read against each other row by row.

## Instrument: WINERED HIRES-Y, R = 68,000

The 2024 He 10830 transit of LHS 1140 b was taken with WINERED in HIRES-Y
mode, whose resolving power is 68,000 (Cherubim et al. 2026, Supplement).
`EXHALE_transit.py` carries the CARMENES 8e4 as its He I 10830 default,
which is right for the other planets in this repository and wrong for this
one, so every run script here sources `../winered_hires_y.sh` before calling
the transit tool. Run the tool by hand and you must do the same:

```bash
. $EX/LHS1140b/winered_hires_y.sh      # EXHALE_TRANSIT_RES_HETR=68000
MPLBACKEND=Agg python3 $EX/EXHALE_transit.py
```

`transit.log` echoes the value it used on the `resolving powers` line; it
must read `He 6.8e+04`. The convolution conserves equivalent width, so the
resolving power moves the modelled widths and peak depths but not the EW
crossings. Curves synthesized at the earlier 8e4 are kept beside the current
ones as `tpm_*_R80k.txt` / `transit_R80k.log`. Background:
`../../docs/lhs1140b_width_measurement_audit.md`.

## What differs between the two models, deliberately

p-winds imposes an isothermal Parker wind at the paper's retrieved
(Mdot = 2.03e8 g/s, T = 6100 K) and solves ionization and the metastable
population on that fixed structure. EXHALE solves the wind: mass, momentum
and energy with photoionization heating and radiative cooling, so Mdot and
T(r) are outputs. A difference in line depth is therefore a difference in the
wind as much as in the line.

`2D approximate method: Mdot` is set on purpose: it applies the full
substellar flux and does not halve the reported Mdot, which is what p-winds
does with `spectrum_at_planet`. The EXHALE default `Rate/2 + Mdot/2` would
halve both and make the two incomparable.

## Domain

The first attempt used `Outer radius [R_p]: 10.0` and is kept in each case's
`attempt_rmax10/`. It is not a result. p-winds puts this planet's sonic point
at **8.0-9.5 R_p** (`../pwinds_oracle/run_pwinds.log`), so that domain ended
at the sonic point: there was no supersonic region for the outflow boundary
to sit in. Every case stalled -- the PLM stage left on its own stall detector
rather than on du, and du then froze near 2-11 while dtu decayed to ~1e-7 --
and `heh1000` ended in a NaN at r = 1.0079.

The runs here use `Outer radius [R_p]: 30.0`, about three times the sonic
radius. The Hill radius is 401 R_p, so nothing tidal is being cut.

## How the cases actually converged (2026-08-24)

Cold marching never finished: after 800k steps the interior was relaxed
(||R|| ~ 0.07-0.10) but the outer region was still settling the cold
hydrostatic IC -- infall at r > 5 R_p -- and the du metric never dropped to
the PLM hand-off. The states the marching *had* reached were good enough for
the steady solver, so each case was finished by `finish_case.sh`:

1. take the case's latest marched snapshot (`output/` updates every 1000
   steps) as the IC;
2. JFNK directly on it (`EXHALE_PTC=1 EXHALE_PTC_JFNK=1`, the
   `examples/04_newton_from_state` route) -- every case converged in
   seconds, info=0, ||R|| = 4e-4 to 1e-3;
3. one `Do only PP` pass for the `_adv` profiles and the Mdot report;
4. `EXHALE_transit.py`.

`heh10`, whose cold march diverged, was seeded from `heh1`'s converged
solution instead of its own snapshot (`finish_case.sh heh10 heh1`); load_IC
carries the composition onto He/H = 10.

The converged winds are slow, cold breezes-to-winds (v reaching only
1.5-1.8 km/s at 30 R_p, T falling to 240-450 K) with
log10 Mdot = 7.69-7.80 across the whole composition axis -- a factor ~3.5
below the paper's imposed 2.03e8 g/s. The marching logs are kept as
`run_march.log`; `run.log` is the JFNK, `pp.log` the post-processing pass.

## GJ 699 proxy SED repeat (2026-08-24)

Five companion cases carry the same compositions on the GJ 699 proxy
spectrum: `solar_gj699`, `heh0p25_gj699`, `heh0p55_gj699`, `heh1_gj699`,
`heh1000_gj699`. Each `input.inp` is a copy of the matching GJ 1132 case with
exactly two lines changed -- `Planet name:` gains the `_gj699` suffix and
`Spectrum file:` points at `../../sed/lhs1140_sed_gj699_at_b.txt`. Nothing
else moves: `sed_read` recomputes the X-ray and EUV luminosities by
integrating the file, so the `Log10 of ... luminosity` keys in `input.inp` are
overridden and were left untouched.

Seed convention: `./finish_case.sh <tag>_gj699 <tag>` -- the IC is the
converged GJ 1132 solution of the same composition, and load_IC carries it
onto the new spectrum. The GJ 1132 directories are read-only seeds and were
not modified. The one exception is `heh0p55_gj699`: seeded from `heh0p55` the
JFNK line search stalled at ||R|| = 1.9e-3 (also at `EXHALE_PTC_DTAU0=0.1`,
1.8e-3), so it was re-seeded from the already converged `heh0p25_gj699`
instead -- same spectrum, neighboring composition -- and then reached info=0,
||R|| = 9.5e-4. The stalled attempt is kept in `attempt_seed_gj1132/`.

The GJ 699 spectrum carries 5.1x the F_XUV of the GJ 1132 one, and the winds
respond: log10 Mdot = 8.58-8.62 for He/H <= 1 (8.43 at He/H = 1000) against
7.74-7.80 on GJ 1132, a factor ~7. Each case also has a second transit under
`tpm_turb/`, produced with `EXHALE_TRANSIT_TURB=1`.

Two bracket cases, `heh0p04_gj699` and `heh0p06_gj699` (He/H = 0.04 and
0.06), were added to measure the equivalent-width crossing on this spectrum
instead of extrapolating it. Both are copies of `solar_gj699/input.inp` with
`Planet name:` and `He/H number ratio:` changed, seeded from `solar_gj699`,
and both reach info=0 (||R|| = 7.9e-4 and 5.0e-4; log10 Mdot = 8.58 and
8.59). Their red-pair EW, 0.761 and 1.104 %A, straddles the measured
1.108 +/- 0.030 %A, and log-log interpolation over the seven GJ 699 cases
puts the crossing at He/H = 0.060 (0.060 with turbulence on) against 0.550
on the GJ 1132 spectrum.

## Diffusive separation on three of the cases (2026-08-25)

`heh0p55_diff`, `heh1000_diff` and `heh0p06_gj699_diff` repeat three cases of
the scan with `He_diffusion: True`, the binary two-component element
diffusion operator of `docs/binary_diffusion_design.md` (milestone M3;
results and their reading in `docs/Update_EXHALE.md` section 69).  Each is a
copy of its seed's `input.inp` with `Planet name:` and that one key changed
-- `He_Kzz` stays at its default 0 (no eddy term) and `He_ambipolar` at its
default on -- restarted from the seed's converged output and finished with
`./finish_case.sh <case> <seed>`.  The seed directories were not modified.

All three converged, info=0 (||R|| = 4.4e-4, 7.9e-4, 8.5e-4).  The wind is
untouched by the separation: log10 Mdot moves by <= 5e-4 dex from the seed
(7.7642 -> 7.7645, 7.7442 -> 7.7437, 8.5885 -> 8.5883) and the He 10830 red
depth by less than 0.1% relative.  The separation itself is small on these
winds: between the base and 20 R_p the element ratio stays within ~1% of the
reservoir He/H (helium mildly enriched outward, most at He/H = 1000), and it
falls to 0.63-0.83 of the reservoir value only in the last cells below the
30 R_p outer boundary.  These are fast, low-gravity winds whose advection
time across the domain is short against the diffusion time, so the mixture
is carried out before it can separate.

Two caveats to carry: the composition outer loop left at its five-pass
ceiling in all three (final drift 4.7e-2, 5.2e-2, 2.2e-2, dominated by the
near-base cells), so these are not composition-converged to the loop's own
1e-3; and the equivalent-width crossings quoted above are from the
diffusion-off cases and were not recomputed.

## The eddy-diffusion scan (2026-08-25)

`heh0p55_diff_kzz1e6` ... `heh0p55_diff_kzz1e11` repeat `heh0p55_diff_ctrl`
with `He_Kzz` set to 1e6 through 1e11 cm^2/s and nothing else changed. The
question they answer is which eddy coefficient this planet needs: with
`He_Kzz` at its default 0 the homopause of a 226 K, `g = 1837 cm/s^2` column
sits below the 1 microbar wind base, helium is gone by 1.05 R_p and the
He 10830 line goes with it. The proposal built on these runs, with the
literature at 1 microbar beside them, is `../kzz_decision.md`; the survey it
rests on is `../kzz_literature.md`; the design memo carries a summary in
`../../docs/binary_diffusion_design.md` section 9.4.

Read the table with `./kzz_scan_table.py` (no arguments = every case). It
reports the solver outcome, the homopause radius where the run's own `D_eff`
equals its `K_zz`, `(He/H)/HeH` on a fixed radius list, `log10 Mdot`, and the
He 10830 red depth and equivalent width, all out of the run directories.

Two things about how the cases were built, both in the scripts:

- `run_kzz_case.sh <tag> <Kzz>` builds a case from the control and solves it;
  `continue_kzz_case.sh <tag> [n]` repeats the invocation until the
  composition drift stops moving; `ladder_kzz.sh <seed> <tag> <Kzz> [n]`
  restarts a case from the case one decade below it, which is what the rows
  at 1e8 and above needed -- from the `K_zz = 0` control the wind solver's
  line search stalls and the composition outer loop stops on its
  solver-failure guard after one pass. The stalled first attempts are kept in
  each case's `attempt_seed_ctrl/`.
- The rows at 1e8 and above also state a looser `Resid tol` (3.0e-3, and
  5.0e-3 at 1e9): on this wind JFNK settles at `||R||` of 2.7e-3 to 4.7e-3
  rather than the 1.0e-3 the control reaches. 2.0e-2 was tried and rejected
  -- there the solver stops responding to the eddy term at all. Every row
  ends `info = 0` with composition drift under 1.3e-3, and each row's
  achieved `||R||` is in the table.

## The composition scan under diffusion (2026-08-25)

`He_Kzz = 1.0e9` was adopted for this planet on 2026-08-25
(`../kzz_decision.md` section 0), so the composition had to be rescanned:
`heh1_diff_kzz1e9`, `heh2_diff_kzz1e9`, `heh2p13_diff_kzz1e9` and
`heh4_diff_kzz1e9` join `heh0p55_diff_kzz1e9` and differ from it only in the
`He/H number ratio` line. The He 10830 equivalent width crosses the measured
1.108 +/- 0.030 %A at **He/H = 2.09**, bracketed by the runs at 2.00 and
2.13, against 0.55 with diffusion off. Record: `../kzz_decision.md` section 6.

- `run_heh_diff_case.sh <seed> <tag> <He/H> [n]` builds and solves one point,
  restarting from the seed case named; `KZZ` and `RESID` override the
  defaults 1.0e9 and 5.0e-3. No case needed the `K_zz` ladder -- seeded from
  a diffusion-off solution of the same composition, each reached `info = 0`
  in its first pass-set.
- `./heh_diff_scan_table.py [0|1e5|1e6|1e7|1e8|1e9|1e10|1e11|flat|<case>...]`
  prints the table: it
  reuses the readers of `kzz_scan_table.py`, takes He/H and `He_Kzz` from
  each case's `input.inp` instead of assuming them, reports the full
  He 10830 metric set, and solves for the equivalent-width crossing by
  Brent's method on a log-log interpolation of the scanned points. No
  argument = the adopted 1e9 scan; `flat` is the fixed-composition probe
  below, whose crossing solve is skipped because it holds a single He/H.

## The same scan at the two alternative eddy coefficients (2026-08-25)

The composition scan above was repeated at the 1e8 and 1e10 alternatives of
`../kzz_decision.md` section 5, so the inferred composition now carries a
`K_zz` bracket. Nine new cases: `heh2p7_diff_kzz1e8`, `heh3_diff_kzz1e8`,
`heh3p5_diff_kzz1e8`, `heh4_diff_kzz1e8`, `heh10_diff_kzz1e8`, and
`heh1_diff_kzz1e10`, `heh1p4_diff_kzz1e10`, `heh1p5_diff_kzz1e10`,
`heh1p6_diff_kzz1e10`. Same `run_heh_diff_case.sh` with `KZZ` overridden,
same `Resid tol: 5.0e-3`; all reached `info = 0` with composition drift under
1.0e-3 in their first pass-set, and none needed the `K_zz` ladder. Seeds and
achieved `||R||` per case: `../kzz_decision.md` section 6.1.

The measured 1.108 +/- 0.030 %A is crossed at **He/H = 2.87 at `K_zz` = 1e8,
2.09 at 1e9 and 1.54 at 1e10** -- a shallow, nearly power-law degeneracy,
`d log(He/H)/d log K_zz = -0.136`. The crossing at 1e8 is reached well before
the line saturates: the equivalent width is still rising with a local slope
near 0.73 there.

`heh10_diff_kzz1e8` is kept for the record but is **not usable as a scan
point**; `heh_diff_scan_table.py` lists it in the table but excludes it from
the crossing fit through `EXCLUDE_FROM_CROSSING`. At He/H = 10
the equilibrium solution develops a sharp hydrogen ionization front near
5.5 R_p; every cell above it trips the advection correction's hydrogen
Damkohler gate (116 of 503 cells kept at equilibrium against 13-32 elsewhere),
and since that gate controls the whole species vector it also freezes the much
slower He 2^3S population, leaving an order-of-magnitude step in the `_adv`
metastable profile that the transit synthesis then reads. Written up in
`../kzz_decision.md` section 6.2. Fixed since, in
`../../docs/Update_EXHALE.md` section 72: the gate now uses the slowest
relaxation rate of the solved species vector. The stored files of that case
still carry the old gate and it stays outside the crossing fit.

## Every decade from `K_zz` = 0 to 1e11 (2026-08-26)

The ladder was carried down to `K_zz` = 0 and up to 1e11 by the same route.
Because the operator adds `K_zz` to the molecular coefficient, and these
He-rich solutions carry `D_eff` = 1.24-1.32e6 cm^2/s at the base, the low
decades were first probed at a **fixed** reservoir He/H = 5 -- one wind
solved at `He_Kzz` = 0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6 (`heh5_diff_kzz0`
through `heh5_diff_kzz1e6`) -- rather than rescanned in composition. The red
EW moves from the `K_zz` = 0 value by 0.00002, 0.0002, 0.002, 0.005, 0.18 and
1.76 times the 0.030 %A measurement error, so **`K_zz` <= 1e4 is not
distinguishable from `K_zz` = 0 by this line** and only 1e6 upward is.

The crossing was then bracketed at `K_zz` = 0 (`heh4p5_diff_kzz0` ...
`heh5_diff_kzz0`), 1e5 (`heh4p5_diff_kzz1e5` ... `heh5_diff_kzz1e5`), 1e6
(`heh2_diff_kzz1e6`, `heh4p3`, `heh4p4`, `heh4p7`, `heh5`, `heh10`), 1e7
(`heh3p4_diff_kzz1e7` ... `heh4p2_diff_kzz1e7`) and 1e11
(`heh1_diff_kzz1e11` ... `heh1p5_diff_kzz1e11`). All reached `info = 0` with
composition drift under 1.0e-3 in their first pass-set; none needed the
`K_zz` ladder, and no crossing runs into line saturation (at 1e6 the EW is
still climbing at He/H = 10). Crossings: **4.68 at `K_zz` <= 1e5, 4.45 at
1e6, 3.78 at 1e7, 2.87 at 1e8, 2.09 at 1e9, 1.54 at 1e10, 1.11 at 1e11**,
against 0.55 with the operator off. The `-0.136` power law is the 1e8-1e11
behaviour only; below 1e7 the curve rolls over onto the `K_zz` = 0 plateau.
Full tables, seeds and achieved `||R||`: `../kzz_decision.md` section 6.1;
figure: `../../docs/figures/lhs1140b_heh_vs_kzz.pdf`.

## Re-post-processed under the slowest-species Damkohler gate (2026-08-26)

Every diffusion-off case in this directory -- `solar`, `heh0p25`, `heh0p5`,
`heh0p55`, `heh0p6`, `heh0p7`, `heh1`, `heh10`, `heh100`, `heh1000`,
`heh10000` and the seven `_gj699` companions -- was re-post-processed with
the binary that forms the advection gate from the slowest relaxing population
rather than from hydrogen alone (`docs/Update_EXHALE.md` section 72). The
winds were **not** re-solved: each case's JFNK solution of record,
`output/*_IC.txt`, was left untouched and read as the input of a single
`Do only PP: True` pass, followed by `EXHALE_transit.py` in both the nominal
and the `EXHALE_TRANSIT_TURB=1` (`tpm_turb/`) configuration. Everything the
pass overwrote is preserved beside it in each case's `pre_gate72_adv/`,
together with the `input.inp` it was produced with.

Five cases were affected, all of them helium-dominated: `heh10`, `heh100`,
`heh1000`, `heh10000` and `heh1000_gj699`. Their red-pair depth and
equivalent width fall by 89-93%, e.g. `heh10` from 75.97% and 21.53 %A to
7.42% and 2.270 %A. The count of cells whose `_adv` metastable density is
still exactly its equilibrium value falls from 173-437 of 504 to 13-22,
and the He-rich end of the equivalent-width curve no longer stands an
order of magnitude above the He/H = 1 case.

The other thirteen cases move by at most 0.06% in depth and equivalent
width, and their equilibrium-valued cell counts do not move at all (44 and
45 of 504 for `heh0p55` and `solar`, before and after). That residual is a
binary-version drift in the equilibrium solve itself (<= 4e-4 relative,
dating from the ghost-cell primitive conversion of section 66 and the
changes after it, all of which postdate these outputs), not the gate. The
equivalent-width crossing does not move: He/H = 0.550 (0.541 with
turbulence) on the GJ 1132 proxy and 0.060 on the GJ 699 one, unchanged,
because both are interpolated inside He/H <= 1 where nothing changed.

One artifact survives the fix. In `heh100` the corrected `_adv` metastable
density jumps from 4.4e2 to 3.0e6 cm^-3 across r = 1.058 -> 1.059 R_p and
returns to 5.3e3 at 1.068 R_p -- a 3.8 dex step, the largest in any case
here. It sits on the hydrogen ionization front (`x_HI` collapses from
3.5e8 to 2.5e-74 cm^-3 over the same cell) where the base velocity is
still sign-changing at +-1e2 cm/s, so it is the breathing base rather than
the gate. The shell is 0.008 R_p thick against a stellar radius of
6.8 R_p, so even fully opaque it can darken at most 0.04% of the disk
against the case's 5.13% red depth; the other four affected cases keep
their largest adjacent-cell step under 0.16 dex.

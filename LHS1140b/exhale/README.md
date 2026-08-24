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

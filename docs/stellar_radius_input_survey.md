# `input.inp` survey: the stellar radius, and a truncated line

**Written 2026-08-23.** Inputs are fixed; **nothing has been recomputed**.
This is the work list for what still needs to be.

## The defect

`EXHALE_transit.py` needs a stellar radius. It normalizes every transit
depth by $\pi R_\star^2$ and caps the impact-parameter integration at
$R_\star/R_p$. `Stellar radius [R_sun]` is an optional key, so a run without
it gets a built-in fallback of **0.44 R_sun** -- an M-dwarf radius. The
script warns, but the warning is one line in a long startup block and the
run continues.

Because the radius enters twice -- once as the normalizing area, once as the
integration cap -- the error is not a single scale factor. Measured on the
tutorial planet (a 1.13 M_sun star, so R_star = 1.155 R_sun):

| quantity | with the 0.44 fallback | with 1.155 R_sun |
|---|---|---|
| He 10830 peak depth, He/H = 1 | 22.77 % | 3.96 % |
| blended red depth, He/H = 1 | 23.27 % | 3.91 % |
| blended red depth, He/H = 1000 | 40.91 % | 8.57 % |

The inflated numbers are what made the LHS 1140b notebook's helium series
look implausible, which is how this was found.

A second defect turned up in the same sweep: four files carry a lone `t` on
line 24, exactly where their siblings carry `Stellar Teff [K]: ...`. It is a
truncated line. The parser matches by label, so an unrecognized line is
ignored and the `t` changed no result -- but two of those four files also
ended without a trailing newline, so appending a key to them merged it onto
the `t`. Both defects predate this session and are in the committed history.

## Survey

131 `input.inp` under the repository. Every file that carries the key pairs
it with the right star, and no file pairs a star mass with another system's
radius:

| `Parent star mass` | `Stellar radius` | system | files |
|---|---|---|---|
| 0.87 | 0.79 | WASP-52 | 25 |
| 0.940 | 0.765 | HD 189733 | 9 |
| 1.13 / 1.200 | 1.155 | HD 209458 (and the tutorial planet, modeled on it) | 20 |
| 1.3521 | 1.458 | WASP-121 | 32 |

## Fixed

**`Stellar radius` added** (value taken from the system's own value already
used elsewhere in the repository, keyed on `Parent star mass`):

- `LHS1140b/heh_series/` -- all six case directories, 1.155
- `examples/tutorial/`, `examples/tutorial_nometals/` -- 1.155
- `examples/01`--`06`, `09`, `10`, `11` -- 0.765 (HD 189733 b)
- `examples/12`, `13_lower_atmosphere/HD209458b`, `14`, `15`, `16` -- 1.155
- `docs/version_compare/v1_nodiff/`, `v2_diff/` -- 1.155
- `backup/regression/tpm_hd189/` -- 0.765

`examples/07_balmer_lya` and `08_full` already carried it, as the Balmer
model requires. The `01`--`12` ladder keeps its one-line-diff property:
every rung and both bases got the same line.

**Truncated line removed** (`t` on line 24):
`docs/version_compare/v1_nodiff/`, `docs/version_compare/v2_diff/`,
`examples/13_lower_atmosphere/HD209458b/`, `backup/lart_runs/hd209_run/`.
The line is deleted, not reconstructed: what it was meant to say is not
recoverable from the file, and the position where siblings put
`Stellar Teff` is a guess, not evidence.

## The key does not touch the wind

In the Fortran, `R_star` is read by `input_read.f90` and used in exactly one
place, `excited_hydrogen.f90`, for the diluted-blackbody Balmer continuum
that photoionizes and heats H(n=2). That model is gated on **both**
`Stellar Teff` and `Stellar radius` being present, so adding the radius to a
file with no `Stellar Teff` leaves the wind untouched.

Verified, not inferred: 25-step runs with and without the line, on both
systems that received it, give byte-identical `Hydro_ioniz.txt` and
`Ion_species.txt`.

- tutorial configuration (1.155): identical
- `examples/06_he23s` configuration (0.765): identical

So only transit spectra need recomputing, never the winds.

## Needs recomputing

| directory | what exists | what to do |
|---|---|---|
| `docs/version_compare/v1_nodiff/` | `tpm_tpm_He10830.txt`, `tpm_tpm_Halpha.txt` | re-run `EXHALE_transit.py` in place; the wind output stands |
| `docs/version_compare/v2_diff/` | the same two | the same |

Those two are the He/H-diffusion version comparison, and their numbers are
quoted in `docs/version_compare.md` and
`docs/version_compare.tex`; those two documents need their figures
and quoted depths refreshed after the transit is redone.

Already done: `LHS1140b/heh_series/heh_1`, `heh_10`, `heh_1000` were
recomputed on 2026-08-23, wind untouched, superseded products kept in each
case's `rstar0p44_default/`.

`examples/tutorial/` and `examples/tutorial_nometals/` also need their
transit redone, but they are already on the recompute list for a
different reason (their stored output predates the triplet default).

## Left as they are, and why

| directory | n | why |
|---|---|---|
| `backup/HD209458b_test/` | 18 | a dated 2026-06 parameter study; no transit products, and editing an archive's inputs makes them disagree with the `EXHALE_setup.out` stored beside them |
| `backup/example_HD209458b/` | 1 | the same, a stored configuration with its own setup report |
| `backup/regression/jfnk_hd189{,_tight}/`, `mol_base_handoff/`, `mol_lyman_werner/`, `mol_metals/` | 5 | golden-pinned cases that synthesize no spectrum. The molecular three are a synthetic hot-Uranus gate (R_p = 0.49 R_J) with no real host star, so any radius would be invented |
| `docs/lower_atmosphere_figs/data_g{1a,1m,2}/` | 3 | figure inputs with no wind output and no transit |

## One open decision

`backup/lart_runs/hd209_run/` sets `Stellar Teff [K]: 6065.0`,
`Stellar Lya flux`, `Lya stellar halfwidth` and `Lya stellar boost` -- the
whole excited-hydrogen configuration -- but has no `Stellar radius`. The
model is therefore **off** in that run: it is gated on both keys. The run
was configured for Balmer physics it did not get.

Adding the radius there is not a documentation fix; it switches the physics
on and the wind changes, so the run has to be redone. Whether that is worth
doing depends on what the run is still for, which is why it is listed here
rather than done.

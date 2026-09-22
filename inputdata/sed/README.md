# Stellar spectra of the planet calculations of record

Decision 16: every
planet calculation of record uses `Spectrum type: Load` with the SED the
reference literature used, or assumed, for that planet; one file per planet
lives here, with the paper, the proxy star, the scaling to the planet's orbit,
the wavelength coverage and the unit conversion stated.

This directory holds the inventory and the files that could be built from data
already in the workspace. Three of them are now the spectrum of their
configuration: `hd209458b_solar_whi2008.txt`, `hd189733b_epseri_salz2016.txt`
and `wasp121b_composite_huang2023.txt` were switched in on 2026-09-06 (step
C3), each recorded in the README of the run
directory that names it. Every other file here is built and unused. A switch is
a configuration change on its own and never a rerun: the outputs standing in a
switched directory were produced under its previous spectrum and are labelled
as such.

## File format, and what the reader demands

`src/modules/radiation/sed_read.f90` (`read_sed`, `sed_next_row`):

* two whitespace-separated columns, **wavelength in A** and **flux at the
  planet in erg cm^-2 s^-1 A^-1** (`F_lambda`, not `F_nu`, not the stellar
  surface flux); the reader converts to `F_E` itself;
* blank lines and lines whose first non-blank character is `#` are skipped, so
  a provenance header is allowed. Every file here carries one;
* wavelengths must be positive and **strictly increasing**; a repeated or
  decreasing value stops the run;
* the reader keeps the rows between `e_top_read` and `e_low` and needs at least
  two of them. `e_top_read` is `e_top` (1.24e3 eV, i.e. 10.0 A, in every planet
  input here) with X-rays on, `e_mid` otherwise;
* `e_low` is 13.6 eV, lowered by `photon_grid_floor_eV()` to the lowest
  ionization threshold of an absorber the run actually carries. Three things
  can lower it, and which of them apply is a property of the **input file**,
  not of the table:
  * **3.399609 eV (3647.0 A)**, H(n=2) at the head of the Balmer continuum,
    whenever `Stellar Teff [K]:` and `Stellar radius [R_sun]:` are **both**
    positive. That pair, and nothing else, arms the excited-hydrogen coupling
    (`input_read.f90`, `use_excited_H`);
  * **4.767775 eV (2600.5 A)**, the He 2^3S metastable, whenever
    `Include He23S?` is True. The key is optional and the metastable defaults
    to on, so an input that omits the line still carries it;
  * the neutral stage of an active low-IP metal. Of the ten elements in
    `species_table.f90` the lowest is K I at 4.341 eV (2856.1 A); the
    `metals.inp` files of the planet folders carry C, N, O, Mg, Ca, Na, Fe,
    whose lowest is Na I at 5.139 eV (2412.6 A).

  Every planet input of record sets both stellar keys, so their binding floor
  is the 3647.0 A of H(n=2), not the 2600 A of the metastable. The Koskinen
  hot-Uranus input and the LHS 1140 b examples set at most one of them, so
  theirs is 2600.5 A;
* a file that ends above `e_low` stops the run: `sed_read.f90` prints the
  file, the band it does not cover in eV and in A, every absorber whose
  ionization threshold lies below the lowest photon the file carries, and the
  remedies, then `error stop 1`. There is no key to continue (decision 17).
  Coverage is therefore a statement about a **configuration and a table
  together**, and `build/sed_coverage_and_order_check.py` checks it that way,
  one entry for each input file of record.

The band flux that drives H2 photodissociation is integrated from the same
file over 912-1110 A when `Stellar LW flux` is absent
(`lyman_werner_band_flux_from_sed`), which also skips `#` lines.

## What a table value means

**The reader treats a row as a histogram value in photon energy.** From the
selected rows `E_1 < ... < E_N` (ascending in energy; the file is ascending in
wavelength), `read_sed` builds bins whose interior edges are the geometric
mean `sqrt(E_k E_k+1)` and whose two end edges are `E_1` and `E_N`, gives the
whole of bin `k` the value `F_E,k = F_k hc/E_k^2` read at the row, and inserts
every active ionization threshold inside the span as a further edge, both
halves keeping the row's value. `e_v` is the geometric centre of a bin and
`de_v` its exact width, so a grid integral is the rectangle rule
`sum(F_XUV de_v)`. A geometric mean in energy is a geometric mean in
wavelength, so the reader's bin of row `k` is
`[sqrt(l_k-1 l_k), sqrt(l_k l_k+1)]` in wavelength.

**Two consumers, two reconstructions.** The photoionization, heating and
setup-report integrals are formed on that grid. The Balmer-continuum rate and
heating of `excited_hydrogen.f90` are not: they read `J_incident`'s
`stellar_flux_eV`, which interpolates the table rows **log-log** and is exact
on a power-law segment, and they integrate it on their own 400-point rule
uniform in photon energy over 3.3996 to 13.598 eV, a step of 0.0255 eV.

**What the shipped files are.** All of them are bin-averaged spectra, except
one component that is a point sample:

| File | source of the values | convention, and where it is stated |
|---|---|---|
| `wasp52b_epseri_yan2022.txt`, `hd189733b_epseri_salz2016.txt` | MUSCLES v22 const-res SED of eps Eri | bin average, 1 A bins. The FITS carries `WAVELENGTH0`/`WAVELENGTH1` (bin edges), `WAVELENGTH` is their exact midpoint and the header describes `FLUX` as "average flux over the bin" (read off the product, 2026-09-06) |
| `lhs1140b_gj1132_cherubim2026.txt` | Mega-MUSCLES v23 const-res SED of GJ 1132 | the same product family and the same convention, 1 A bins |
| `wasp121b_wasp17_huang2023.txt` | MUSCLES v24 const-res SED of WASP-17 | the same, 1 A bins |
| `hd209458b_solar_whi2008.txt` | WHI 2008 solar reference spectrum | bin average, 1 A bins: the rows sit at 0.5, 1.5, 2.5 A, the centres of the 1 A intervals the reference spectrum is tabulated on |
| `hot_uranus_solar_koskinen2022.txt` | Ribas band fluxes on the Gueymard composite | bin average, 10 A bins: the Gueymard file is tabulated at 0.5, 1.5, 2.5 nm, the centres of its 1 nm intervals |
| `hd189733b_bourrier2020.txt` | the released Bourrier et al. 2020 DEM reconstruction, 0.016 A rows, with the eps Eri photosphere above 1600 A | the released `readme_MNRAS.txt` states only "the wavelength, in angstroms" and "the stellar flux, in erg/cm2/A/s" and gives no bin convention. At 0.016 A the two readings differ by less than 1e-9 of a row, far below anything the reader resolves; above 1600 A the file is the 1 A eps Eri bin average |
| `wasp121b_composite_huang2023.txt` | WHI 2008 below 1700 A; the MUSCLES v24 const-res SED of WASP-17 above it | bin average throughout, 1 A bins: both components are bin averages on a 1 A step, the WHI rows at 0.5, 1.5, ... A and the MUSCLES rows at 1700.0, 1701.0, ... A |
| `wasp121b_solar_huang2023.txt` | WHI 2008 below 1700 A; `pi B_lambda(6459 K) (R_star/a)^2` above it | bin average below 1700 A; **point sample** above it, the blackbody evaluated at 1700.5, 1701.5, ... A by `build_wasp121b_solar_huang.py`. For a curve this smooth the point value and the 1 A bin average differ by order `(dl/l)^2 ~ 3e-7` |

**The consequence, measured.** A bin-averaged table is exact for the band
energy of its own bin, so reading it as a histogram loses nothing there; what
is left is the difference between a histogram in energy on geometric-mean
edges and a histogram in wavelength on the declared bin, which is second order
in `dl/l`. A point-sampled table is not exact for band energy, and its error
is first order in the curvature of `F_lambda` over a bin, again `(dl/l)^2`
relative for a smooth continuum but unbounded across a line narrower than the
step.

Exactly this is measured
through the production consumers, on four synthetic tables built from an
analytic `F_lambda` by exact bin averaging (a narrow emission line, a sharp
step at the He I threshold, one bin ten times its neighbours, and a steep
`F_lambda ~ l^4` continuum), at the 1 A step of the shipped products and at a
half and a quarter of it. Measured on the 2026-09-06 build:

* the band energy, the photon number, `P_HI`, `P_HeI`, `P_HeII`, `P_HeTR` and
  the H I photoheating integral reproduce the declared histogram to between
  1.4e-7 and 8.2e-5 on every one of the four tables, and their deviation from
  the sampled spectrum falls at **order 2** when the table step is halved
  (band energy 5.7e-7, 1.4e-7, 3.0e-8 on the narrow-line table);
* the log-log interpolation the Balmer band reads is exact on a power-law
  segment: on the smooth tables the refined Balmer integrals reproduce the
  sampled spectrum to 2.6e-8;
* the production 400-point Balmer rule steps 0.0255 eV, so it cannot resolve
  an emission line narrower than that. On the narrow-line table its rate is
  3.5e-4 and its heating 1.4e-3 away from the same field integrated on 40001
  points, and the gap **grows** as the table is refined (rate 3.5e-4, 9.7e-4,
  2.2e-3), because a finer table makes the interpolated line narrower while
  the rule's step does not move. A measured stellar spectrum carries such a
  line at Ly-alpha. This is an open defect of `excited_hydrogen.f90`, recorded
  here and not fixed by that test.

## Source rule (user, 2026-09-06)

A missing SED is built from the MUSCLES / Mega-MUSCLES data held locally in
`~/Exoplanetary_Atmosphere/MUSCLES` (v22 to v25 FITS of 36 stars, WASP-17,
WASP-43, WASP-77A, WASP-127, HAT-P-12, HAT-P-26, HD 149026 among them) or
downloaded from the MUSCLES site, and then normalized and scaled **the way the
planet's reference paper processed its spectrum** (band-flux normalization,
distance scaling, the joins between components). The proxy star and every
scaling step are stated in the file header and in the inventory row.

## Inventory

| Planet | Reference for the SED | Proxy / construction | Data in workspace | Status | File |
|---|---|---|---|---|---|
| WASP-52 b | Yan et al. 2022, ApJ 936, 177 (`2022ApJ...936..177Y`) | eps Eri (MUSCLES v22), x (d/a)^2 | yes (MUSCLES v22 FITS, local holdings) | **BUILT** | `wasp52b_epseri_yan2022.txt` |
| hot Uranus (Koskinen Model A) | Koskinen et al. 2022, ApJ 929, 52 (`2022ApJ...929...52K`) | mean solar, Ribas 2005 band fluxes on a Gueymard composite, at 0.05 au | partly (the true SOLAR2000 + Woods and Rottman spectrum is not here) | **BUILT** (stand-in, coverage extended to 2995 A) | `hot_uranus_solar_koskinen2022.txt` |
| HD 209458 b | Oklopcic and Hirata 2018, ApJL 855, L11 (`2018ApJ...855L..11O`) | solar irradiance at minimum (WHI 2008), x (1 au/a)^2 | yes (WHI reference spectrum, p-winds reference data) | **IN USE** by `HD209458b/` and `benchmarks/hd209/` (switched 2026-09-06, step C3) | `hd209458b_solar_whi2008.txt` |
| LHS 1140 b | Cherubim et al. 2026, Science, doi:10.1126/science.aea9708 | GJ 1132 (Mega-MUSCLES v23) normalized to the measured X-ray flux, x (d/a)^2 | yes (already built under `LHS1140b/sed/`) | **BUILT** (copy with header) | `lhs1140b_gj1132_cherubim2026.txt` |
| HD 189733 b | Salz et al. 2016, A&A 586, A75 (`2016A&A...586A..75S`) | eps Eri (MUSCLES v22) for the shapes, each of Salz's four regions normalized to his Table 3 luminosities, x (d/a)^2 | yes (MUSCLES v22 FITS, local holdings) | **IN USE** by `HD189733b/` and `benchmarks/hd189/` as the run of record (switched 2026-09-06, step C3); the Bourrier file is the observed-star alternative | `hd189733b_epseri_salz2016.txt` |
| HD 189733 b (recommended literature SED) | Bourrier et al. 2020, MNRAS 493, 559 (`2020MNRAS.493..559B`), MOVES III | the authors' own DEM reconstruction, mean of Visits B-E, x (1 au/a)^2 at Salz's a = 0.031 AU, eps Eri photosphere above 1600 A | yes (the four released tables, `references/Bourrier2020_supplemental_files/`) | **BUILT** | `hd189733b_bourrier2020.txt` |
| WASP-121 b | Huang et al. 2023, ApJ 951, 123 (`2023ApJ...951..123H`) | solar SED (WHI 2008) below 1700 A rescaled to F_XUV = 1.6e6; 6459 K blackbody above 1700 A in place of LLmodels | yes (WHI reference spectrum) | **BUILT** (candidate, recommended) | `wasp121b_solar_huang2023.txt` |
| WASP-121 b (alternative) | the same paper | WASP-17 (MUSCLES v24, F6V) throughout, same three normalizations | yes (MUSCLES v24 FITS, local holdings) | **BUILT** (candidate) | `wasp121b_wasp17_huang2023.txt` |
| WASP-121 b (provisional composite) | the same paper | candidate A's solar (WHI 2008) XUV below 1700 A, candidate B's bolometrically scaled WASP-17 photosphere above it, neither renormalized to close the join | yes (both components already built) | **IN USE** by `WASP-121b/` and `benchmarks/wasp121/`, labelled provisional (switched 2026-09-06, step C3); candidates A and B stay as the sensitivity inputs | `wasp121b_composite_huang2023.txt` |

Measured band fluxes of the built files (trapezoid on the stored bins,
`build/sed_coverage_and_order_check.py`; erg cm^-2 s^-1):

| File | F(1-100 A) | F(100-912 A) | F(1-912 A) | F(10-912 A) | F(912-2600 A) | beta_m | last row |
|---|---|---|---|---|---|---|---|
| `wasp52b_epseri_yan2022.txt` | 7.471e3 | 2.680e4 | 3.431e4 | 3.396e4 | 2.633e5 | 0.218 | 3699.9 A |
| `hot_uranus_solar_koskinen2022.txt` | 3.986e2 | 1.191e3 | 1.600e3 | 1.541e3 | 1.154e6 | 0.249 | 2995.0 A |
| `hd209458b_solar_whi2008.txt` | 1.075e2 | 1.232e3 | 1.340e3 | 1.340e3 | 1.194e6 | 0.080 | 4999.5 A |
| `lhs1140b_gj1132_cherubim2026.txt` | 9.430 | 1.005e1 | 1.948e1 | 1.941e1 | 2.642e1 | 0.484 | 30000.0 A |
| `hd189733b_epseri_salz2016.txt` | 5.629e3 | 1.507e4 | 2.073e4 | 2.043e4 | 1.841e5 | 0.272 | 4999.9 A |
| `wasp121b_solar_huang2023.txt` | 1.283e5 | 1.470e6 | 1.600e6 | 1.599e6 | 1.833e8 | 0.080 | 4999.5 A |
| `wasp121b_wasp17_huang2023.txt` | 1.620e5 | 1.438e6 | 1.600e6 | 1.600e6 | 1.205e8 | 0.101 | 5000.0 A |
| `wasp121b_composite_huang2023.txt` | 1.283e5 | 1.470e6 | 1.600e6 | 1.599e6 | 1.028e8 | 0.080 | 5000.0 A |
| `hd189733b_bourrier2020.txt` | 1.008e4 | 2.070e4 | 3.078e4 | 2.748e4 | 1.846e5 | 0.327 | 4999.9 A |

`beta_m = F(1-100)/F(1-912)`, the hardness ratio Yan et al. define. The
`wasp121b_wasp17` file starts at 10.0 A, so its F(1-100 A) entry is
F(10-100 A) and its `beta_m` is a lower bound.

**F(10-912 A) is the band a run actually integrates.** Every input of record
states `E_high = 1.24e3 eV` with X-rays on, so `e_top_read` is 1240 eV and the
photon grid starts at 9.9987 A: nothing shortward of that enters the run,
whatever the file carries there. The two columns differ by 11 per cent for
`hd189733b_bourrier2020.txt`, whose released tables reach 0.016 A, and by less
than 1 per cent for the rest. Both are printed by
`build/sed_coverage_and_order_check.py`, and the normalization assertions state
which band their source normalized over (Salz 1-912 A, Huang solar 1-912 A,
Huang WASP-17 10-912 A).

**Coverage, per configuration of record** (measured 2026-09-06 by
`build/sed_coverage_and_order_check.py`, which reads each `input.inp` and its
`metals.inp` and derives the floor as `photon_grid_floor_eV` does). Three
entries do not pass:

| Configuration | H(n=2) armed | floor set by | required long endpoint | table it would take | verdict |
|---|---|---|---|---|---|
| `WASP-52b` | yes | H(n=2) | 3647.0 A | `wasp52b_epseri_yan2022.txt` (3699.9 A) | covers. Its `Spectrum file:` line still names `.../EXHALE/WASP-52b/eps_eri_sed_fxuv1p0.txt`, a path that no longer exists after the 2026-07-14 move to `EXHALE_v1.00/`, so the run as configured stops in `read_sed`'s `open` |
| `benchmarks/wasp52` | yes | H(n=2) | 3647.0 A | `../../WASP-52b/eps_eri_sed_fxuv1p0.txt` (3699.9 A) | covers |
| `HD189733b` | yes | H(n=2) | 3647.0 A | `hd189733b_bourrier2020.txt` or `hd189733b_epseri_salz2016.txt` (both 4999.9 A) | covers |
| `HD209458b` | yes | H(n=2) | 3647.0 A | `hd209458b_solar_whi2008.txt` (4999.5 A) | covers, since the rebuild of 2026-09-06 (`--w-max 5000.0`; the WHI source runs to 23999.5 A). The earlier cut at 3000 A stopped the run in `read_sed` naming H(n=2), measured on the production binary |
| `WASP-121b` | yes | H(n=2) | 3647.0 A | `wasp121b_solar_huang2023.txt` (4999.5 A), `wasp121b_wasp17_huang2023.txt` (5000.0 A) or `wasp121b_composite_huang2023.txt` (5000.0 A) | long end covers. The WASP-17 variant starts at exactly 10.0 A against the 9.9987 A the 1240 eV top asks for, so it states no field over the topmost 1.3e-3 A of the requested band; `read_sed` does not refuse that, it simply puts the grid top at 1239.84 eV. The composite starts at 0.5 A, since its XUV is the solar component, so it covers the short endpoint |
| `benchmarks/koskinen2022_model_a` | **no** (neither stellar key is stated) | He 2^3S | 2600.5 A | `hot_uranus_solar_koskinen2022.txt` (2995.0 A) | covers. This answers the open question about the stand-in: it is short of the Balmer continuum, but this configuration does not arm H(n=2), so the requirement never applies. Confirmed on the production binary, which reads the file without complaint |
| `benchmarks/koskinen2022_model_a/matched_*` (three) | no | `E_low` of the input, no sub-Lyman absorber (`Include He23S? False`) | 911.6 A | `solar_ribas2005_bands_0p05au.txt` (1995.0 A) | covers; `E_high = 2.0e3` eV there, so the short endpoint is 6.20 A and the 5.0 A first row covers it |
| `LHS1140b/examples/*` (three) | **no** (`Stellar radius` is set, `Stellar Teff` is not) | He 2^3S | 2600.5 A | `lhs1140_sed_gj1132_at_b.txt` = `lhs1140b_gj1132_cherubim2026.txt` (30000.0 A) | covers |

The two short files in the tree do not cover a metastable run and stop it:
`inputdata/scaled_solar_hd189.sed` leaves the band 4.77-6.53 eV uncovered, and
`benchmarks/koskinen2022_model_a/solar_ribas2005_bands_0p05au.txt` leaves
4.77-6.21 eV. The second is used only by the `matched_*` runs, which carry
`Include He23S? False` and so never ask for that band.

## Papers consulted, and what is not in `references/`

Read in the local publisher PDFs for this inventory (bibcodes confirmed through
ADS): Huang et al. 2023 (`Huang_2023_ApJ_951_123.pdf`), Yan et al. 2022
(`Yan 2022.pdf`), Salz et al. 2016 (`Salz_2016A&A_586_A75.pdf`), Oklopcic and
Hirata 2018 (`Oklopcic_2018_ApJL_855_L11.pdf`), Koskinen et al. 2013a and 2013b
(`Koskinen_2013Icarus_226_1678.pdf`, `..._1695.pdf`), Koskinen et al. 2022
(`Koskinen_2022_ApJ_929_52.pdf`), Cherubim et al. 2026
(`Cherubim_2026Science.pdf` and its Supplement), Salz et al. 2019
(`Salz_2019A&A_623_A57.pdf`), Bourrier et al. 2020
(`Bourrier_2020MNRAS_493_559.pdf`) and Salz et al. 2018
(`Salz_2018A&A_620_A97.pdf`).

Not in `references/`, and named where a build depends on them: Lampon et al.
2020 and 2021 (the He 10830 retrievals for HD 209458 b and HD 189733 b),
Woods and Rottman 2002 and Tobiska et al. 2000 (the solar spectra behind
Koskinen's average Sun and Salz's EUV shape), Sanz-Forcada et al. 2011 and
Hunsch et al. 1999 (the HD 189733 X-ray and EUV luminosities Salz et al. 2018
quotes), and Andretta et al. 2017 (the He I equivalent widths and the eps Eri
X-ray luminosity Salz et al. 2018 compares HD 189733 with).

**How the Bourrier et al. 2020 tables got here.** They are MNRAS supplementary
data of `doi:10.1093/mnras/staa256`, not a CDS catalogue, and they cannot be
fetched from this host: measured 2026-09-06, `academic.oup.com` and
`oup.silverchair-cdn.com` return HTTP 403 to every request (bot protection),
VizieR answers "Table or Catalog not found: J/MNRAS/493/559", the ADS record
carries no CDS data link, and `arXiv:2001.11048` has no ancillary files. The
user downloaded them from a browser session; they are in
`references/Bourrier2020_supplemental_files/` as
`HD189733_XUV_spectrum_Visit{B,C,D,E}.dat` plus `readme_MNRAS.txt`.

## Build and check

```
python3 build/build_wasp52b_epseri.py
python3 build/build_hot_uranus_solar.py
python3 build/build_hd209458b_solar_whi.py
python3 build/build_lhs1140b_gj1132.py
python3 build/build_hd189733b_epseri_salz.py
python3 build/build_wasp121b_solar_huang.py
python3 build/build_wasp121b_wasp17.py
python3 build/build_hd189733b_bourrier2020.py
python3 build/sed_coverage_and_order_check.py      # exits nonzero on any failure
```

`sed_coverage_and_order_check.py` asserts, one printed line per assertion:

* the `read_sed` file rules above for every file here (two positive strictly
  increasing columns, no negative flux, at least two rows);
* **coverage for each configuration of record**, not for each file: it reads
  the `input.inp` and the `metals.inp` of every entry of the table above,
  derives the long endpoint exactly as `photon_grid_floor_eV` does and the
  short endpoint from `E_high` with X-rays on, prints the endpoint used and
  the absorber that set it, and checks the table the input names now together
  with the table it would name after its C3 switch;
* the normalizations the built files were made to, over the band the source
  normalized, with F(1-912 A) and F(10-912 A) both printed;
* that the built files reproduce the tables they came from.

It exits nonzero on any failure, and it currently reports one: the WASP-17
file starting at 10.0 A against the 9.9987 A of a 1240 eV grid top (the
MUSCLES v24 table itself begins at 10.0 A, so this is the data's limit;
`read_sed` puts the grid top at 1239.84 eV). The two failures of the first
run were closed the same day: `WASP-52b/input.inp` now names its own
`eps_eri_sed_fxuv1p0.txt`, and the WHI file is built to 5000 A. The
remaining one is stated in the coverage table
above; none of them is fixed by the check itself.

The reconstruction the reader applies to a table, and what the production
consumers make of it, are tested separately.

---

## WASP-52 b: `wasp52b_epseri_yan2022.txt` (BUILT)

**Paper.** Yan, Seon, Guo, Chen and Li (2022), ApJ 936, 177,
`2022ApJ...936..177Y`.

**What the paper says.** Section 2.1: "The integrated flux in the XUV band is
an important input in the simulations. In the absence of direct observation of
the stellar XUV, we use the spectrum of eps Eri from the MUSCLES Treasury
Survey (France et al. 2016) to calculate the F_XUV received by the planet. The
F_XUV is about 34,168 erg cm-2 s-1 at the orbital distance of about 0.0272 au."
Section 2.3: "Here the incident ultraviolet spectrum (912-2600 A) that ionizes
He(23S) is also taken from the MUSCLES Treasure Survey, as shown in
Figure 1(b). The integrated flux in this far-ultraviolet (FUV) and NUV band is
about 2.6 x 10^5 erg cm-2 s-1." Figure 1 caption: "The stellar XUV and FUV/NUV
SED, which is the spectrum of eps Eri taken from the MUSCLES Treasury Survey
(France et al. 2016) and scaled at the orbital distance of WASP-52b."
(The Greek letters are as printed; the pdftotext rendering of the subscripts is
not reproduced.)

**Source data.** MUSCLES Treasury Survey v22, `v-eps-eri`, broadband
**const-res** SED, flux at Earth on a 1 A grid. Local holdings:
`~/Exoplanetary_Atmosphere/MUSCLES/v-eps-eri/`. Originator: France et al.
(2016), ApJ 820, 89 (`2016ApJ...820...89F`).

**Scaling.** x (d/a)^2 = 5.932852e14 with d = 3.212 pc and a = 0.0272 AU.
17 negative source rows (background subtraction) are set to zero. Window
5.939229-3699.939229 A.

**Relation to the existing files.** `WASP-52b/eps_eri_sed_fxuv1p0.txt` is
exactly this construction: measured 2026-09-05, the ratio of that file to the
const-res product on the same 3695 rows is one constant, 5.932852e14, on
every row to 1e-6, and the negative source rows are zero there. (The
`adapt-const-res` product matches only 93 per cent of the rows, so the parent
is the plain `const-res` one.) `eps_eri_sed_fxuv0p5.txt` and
`..._fxuv0p25.txt` are that file times 0.5 and 0.25 to 2.5e-7; they belong to
Yan's F_XUV parameter study (their best fit is 0.5 F_0, which is the 0.25x
file), not to the spectrum. Only the unscaled proxy is kept here.

**Check against the paper.** F(1-912 A) = 3.431e4 against Yan's 34,168
(+0.4 per cent); beta_m = 0.218 against their 0.22; F(912-2600 A) = 2.633e5
against their "about 2.6 x 10^5".

**Open at the switch.** `WASP-52b/input.inp` line 12 points at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE/WASP-52b/eps_eri_sed_fxuv1p0.txt`,
a path that no longer exists (the tree is `EXHALE_v1.00/`); only
`benchmarks/wasp52/input.inp`, which uses a relative path, still resolves.

---

## Hot Uranus of Koskinen et al. 2022: `hot_uranus_solar_koskinen2022.txt` (BUILT, stand-in)

**Paper.** Koskinen, Lavvas, Huang, Bergsten, Fernandes and Yelle (2022),
ApJ 929, 52, `2022ApJ...929...52K`.

**What the paper says.** Section 3.2 (Numerical Modeling, p. 8): "In both
cases, we use the same mean solar XUV spectrum as Koskinen et al. (2013a,
2013b) and assume uniform redistribution of energy around the planet". Section
2.5 (Upper Atmosphere Escape Model, p. 6): "For energy-limited escape, we use
the mean solar XUV flux of 1.6 W m-2 at 0.05 au integrated over the wavelength
range of 0.1-91.1 nm (Ribas et al. 2005)". Page 9: "With a stellar XUV
luminosity of L_XUV = 1.12 x 10^21 W, the flux incident on Model A at
a = 0.05 au is F_XUV = 1.6 W m-2."

The spectrum they point to is defined in Koskinen et al. (2013a), Icarus 226,
1678 (`2013Icar..226.1678K`), section 2.1: "We simulated heating and
photoionization self-consistently by using the model density profiles and the
UV spectrum of the average Sun. The spectrum covers wavelengths between 0.1 and
3000 A. The XUV spectrum between 0.1 and 1050 A was generated by the SOLAR2000
model (Tobiska et al., 2000). It includes strong emission lines separately and
weaker lines binned by 50 A. The Lyman a line was included with a wavelength
spacing of 0.5 A from Lemaire et al. (2005) and the rest of the spectrum was
taken from Woods and Rottman (2002)."

**What we have instead.** Neither SOLAR2000 nor the Woods and Rottman table is
in the workspace, so this file is **not** that spectrum. It is the assembly
already used by `benchmarks/koskinen2022_model_a`: the band fluxes of Ribas
et al. (2005), ApJ 622, 680 (`2005ApJ...622..680R`) Table 4 for the 4.56 Gyr
Sun over 1-20, 20-100, 100-360, 360-920 and 920-1180 A; the shape inside each
band from the Gueymard (2003/2004) synthetic composite
(`VULCAN/atm/stellar_flux/Gueymard_solar.txt`, a surface flux, /46250 gives
1 au); moved to 0.05 au and normalized so F(1-911 A) = 1.6e3 erg cm^-2 s^-1.
The provenance of that part is
`benchmarks/koskinen2022_model_a/solar_ribas2005_bands_0p05au.README`.

**What this file adds.** Coverage only. The benchmark file stops at 1995 A,
above the 2600 A the metastable needs; with `Include He23S? True` the reader
stops the run. Longward of 1180 A the benchmark
file is the Gueymard composite times one constant, measured over its 82 rows
there as 8.4453229e-3 with a relative spread of 7.5e-7 (the round-off of the
7-digit stored file), so the continuation from 2005 to 2995 A is the same
source times the same constant. **No flux is invented**, and the 200 data lines of
5-1995 A are byte-for-byte those of the benchmark file (`diff` on the data
lines, measured; the check asserts the numbers to 1e-6). Koskinen's own
spectrum runs to 3000 A, so the extension is in the direction of the reference,
not away from it.

**Check against the paper.** F(1-912 A) = 1.6000e3 erg cm^-2 s^-1, which is the
1.6 W m^-2 of section 2.5 to the digits given.

**Open at the switch.** `benchmarks/koskinen2022_model_a/input.inp` states
`Orbital distance [AU]: 0.0480` while the paper's Model A and this spectrum are
at 0.05 au. `a_orb` sets the tidal term (`atilde`), the Balmer dilution and the
reported luminosities, not the loaded flux, so the run's geometry sits at
0.048 au while its radiation is that of 0.05 au: 4.2 per cent in distance,
8.5 per cent in flux. The three `matched_*` subdirectories already load the
1995 A file and carry `Include He23S? False`, which is why their floor is
13.6 eV and the 1995 A table covers it.

---

## HD 209458 b: `hd209458b_solar_whi2008.txt` (BUILT)

**Paper.** Oklopcic and Hirata (2018), ApJL 855, L11, `2018ApJ...855L..11O`.

**What the paper says.** Section 3.2 (Steady-state Hydrogen Distribution,
p. 3): "For HD209458, a G0-type star, we use the SORCE solar spectral
irradiance data from the LASP Interactive Solar Irradiance Datacenter.
Because HD209458 is an inactive star (Czesla et al. 2017 and references
therein), we use the solar spectrum data recorded during a solar minimum. To
fill in a gap in the data in the wavelength range ~400-1150 A, we use the
scaling relations between the Lya flux and fluxes in EUV bands from Linsky
et al. (2014)."

**Source data.** The Solar Irradiance Reference Spectra for the 2008 Whole
Heliosphere Interval (Woods et al. 2009, GRL 36, L01101, `2009GeoRL..36.1101W`),
file `ref_solar_irradiance_whi-2008_ver2.dat` from LASP LISIRD, present in the
workspace as a p-winds reference spectrum
(`~/Exoplanetary_Atmosphere/p-winds_data/p-winds_reference_spectra/`). Column
"March 30 - April 4 (6-day average)". Grid: 0.1 nm intervals on 0.05 nm
centers, i.e. 1 A bins centered at 0.5, 1.5, ... A. Irradiance in W m^-2 nm^-1
at 1 au; x100 gives erg cm^-2 s^-1 A^-1.

**Scaling.** x (1 au / a)^2 = 451.348084 with a = 0.04707 AU. No stellar-radius
correction: the star is given the Sun's surface flux, as in the reference.

**Why this column and this scaling.** p-winds ships exactly this table as
`data/solar_spectrum_scaled_lambda.dat`, its HD 209458 b input; measured
2026-09-05, that file is this WHI column times 451.348084 on all 3000 rows to
1e-6, and the file built here reproduces it to 0.0e0 relative. The choice of
column is therefore the published one, not ours. The WHI release also offers an
"April 10-16" quiet-Sun average, closer to cycle minimum; switching to it is a
one-line change in the build script and would be a departure from the p-winds
construction.

**Departures from the reference, stated.** (1) Oklopcic and Hirata filled
400-1150 A with the Linsky et al. (2014) relations; the WHI reference spectrum
carries values there already, and its own source flag marks part of 912-1200 A
as filled rather than measured. The two fillings are not the same numbers.
(2) Other reference treatments of this planet differ: Koskinen et al. (2013a)
used the average Sun (SOLAR2000 + Woods and Rottman) at 0.047 au, and Salz
et al. (2016) reconstructed the SED from the measured luminosities (their
Table 3 for HD 209458: log L_X < 26.40, log L_Lya = 28.77, log L_EUV < 27.84,
log F_XUV < 3.06). This file gives F(1-912 A) = 1.341e3, which is 1.17x Salz's
upper limit; the two are not the same statement about the star, and the
difference is the point of the comparison, not an error to hide.

**Open at the switch.** `HD209458b/input.inp` states `Orbital distance [AU]:
0.0480` against the 0.04707 AU used here (4.0 per cent in flux); one of the two
has to move. The power law it runs today gives F_XUV = 1.558e3 erg cm^-2 s^-1
at 0.0480 AU from log L_X = 27.20 and log L_EUV = 27.93, against 1.340e3 for
this file at 0.04707 AU, and a much harder shape (beta_m 0.157 against 0.080).

**Coverage.** The file was first built to 3000 A and stopped the run for
`HD209458b/input.inp`, whose two stellar keys arm the H(n=2) coupling and
put the grid floor at 3647.0 A (measured 2026-09-06 on the production binary:
`read_sed` refused it naming H(n=2) and the band 2999.5 to 3647.0 A). Since
the same day the build cuts at 5000 A (`--w-max 5000.0`, the WHI source runs
to 23999.5 A) and the configuration inventory of the coverage script passes.

---

## LHS 1140 b: `lhs1140b_gj1132_cherubim2026.txt` (BUILT, copy)

**Paper.** Cherubim et al. (2026), Science, doi:10.1126/science.aea9708. ADS
carries the preprint as `2026arXiv260714326C`; the published bibcode is not in
ADS yet.

**What the paper says.** Supplementary Materials, "X-ray observations and
analysis": "We estimate the XUV SED of LHS 1140 by comparison to two M dwarf
stars - GJ 699 and GJ 1132 - that are of similar spectral type (M4 and M3.5,
respectively) and rotation period (145 and 125 days, respectively) to LHS 1140
... Both stars have measured X-ray fluxes and panchromatic SEDs from X-ray to
infrared wavelengths from the Mega-Measurements of the Ultraviolet Spectral
Characteristics of Low-mass Exoplanetary Systems (Mega-MUSCLES) Treasury Survey
... We normalize the panchromatic SEDs of GJ 699 and GJ 1132 so their
integrated X-ray flux equals that measured for LHS 1140. The multiplicative
normalization factor (A) for GJ 699 (A_699) and GJ 1132 (A_1132) are
A_699 = 0.056 and A_1132 = 0.59, respectively."

**Source data.** Already built, with a long build record, as
`LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt`: Mega-MUSCLES v23 const-res SED of
GJ 1132 x A = 0.5900 (the catalog X-ray flux ratio, which is what reproduces
the paper's quoted constants), carried to the orbit by (d/a)^2 = 1.063974e15
with d = 14.96 pc and a = 0.0946 AU, 978 negative rows clipped and 4 repeated
wavelengths merged. The reasoning, the version difference between Mega-MUSCLES
v23 and v25, and the GJ 699 alternative are in `LHS1140b/sed/README.md`.

**What this file is.** The same two columns with a provenance header, so that
every planet of record has its spectrum in one place; the numbers are
identical (asserted at tolerance 0 by the check). `LHS1140b/sed/` stays as the
build tree, since it holds `make_sed.py` and the Mega-MUSCLES inputs. When the
inputs are switched, `LHS1140b/examples/*/input.inp` should point here and the
duplicate should go.

**Check against the paper.** F(10-1300 A) = 33.4 erg cm^-2 s^-1 against the
paper's fiducial 0.033 W m^-2 = 33.

**Noticed while copying.** The last section of `LHS1140b/sed/README.md` states
that the shipped products are integral-matched with an `A_eff` in their
headers. They are not: the shipped `lhs1140_sed_gj1132_at_b.txt` header reads
`A = 0.5900 (catalog X-ray flux ratio; authors' released script)` and
`make_sed.py` lines 58-63 and 95-98 say the same. The README paragraph is
stale; the file and the script are what this copy follows.

---

## HD 189733 b: `hd189733b_epseri_salz2016.txt` (BUILT, candidate)

**Reference.** Salz, Czesla, Schneider and Schmitt (2016), A&A 586, A75,
`2016A&A...586A..75S`, the escape simulation this planet's runs are compared
with. (Lampon et al. 2020/2021 are not in `references/`. Salz et al. 2018,
the He 10830 detection for this planet, is now present but constrains no SED:
it runs no atmosphere model, and section 4.2 states "we do not fit a proper
atmosphere model here". See the Bourrier et al. 2020 subsection below.)

**What the paper says.** Section 2.2 lists the four ranges of the
reconstruction: "1. For the X-ray spectrum (0 to 100 A) we use a 2 or 4 x 10^6 K
plasma emission model from CHIANTI (Dere et al. 1997, 2009) for inactive or
active host stars respectively (active: L_X > 10^28 erg s-1). The SED is
normalized to the observed X-ray luminosities of the host stars." "2. For the
hydrogen Lya emission line of the host stars a Gaussian with a full width at
half maximum (FWHM) of 9.4 A is used. The irradiation strength is normalized
according to the host star's total Lya luminosity." "3. For the EUV range
(100 A to ~912 A) the luminosity of the host stars is predicted based on the
Lya luminosity given by step 2 (Linsky et al. 2014). ... The shape of the SED
is taken from an active or inactive solar spectrum (Woods & Rottman 2002),
depending on the activity of the host star defined in step 1. This shape is
continued to a connecting point with the photospheric blackbody, where we use a
visual best fit for the connection point between 1500 to 4000 A." "4. For the
remaining part of the spectrum up to 5 x 10^4 A we choose a blackbody according
to the host star's effective temperature."

**What the paper supplies.** Table 3 gives HD 189733 b log L_X = 28.18
(0.124-2.48 keV, i.e. 5.0-100.0 A), log L_Lya = 28.43, log L_EUV = 28.61
(100-912 A) and log F_XUV = 4.32 (<912 A at the planet); Table 2 gives
spectral type K0-2V, T_eff = 5040 K, d = 19 pc and a = 0.031 AU.

**What is still missing, and what stands in for it.** The CHIANTI plasma model
spectrum and the Woods and Rottman solar shape are not in the workspace, so
only Salz's *normalizations* are used, and the *shapes* come from a measured
K dwarf under the source rule: the MUSCLES Treasury Survey v22 const-res
broadband SED of eps Eri (K2V), whose X-ray part is an APEC plasma model and
whose photosphere is a PHOENIX model, standing in for the CHIANTI model of
step 1 and the blackbody of step 4.

**Why eps Eri and not HD 97658.** Measured from the released products,
eps Eri at d = 3.212 pc has log L_X(5-100 A) = 28.19 against Salz's 28.18 for
HD 189733, so step 1 is a 3 per cent correction (c_X = 0.974) rather than a
substitution, and the file's flux at 100 A moves by only 1.28 across the
step 1 / step 3 join. HD 97658 (K1V) at d = 21.0 pc has
log L_X(5-100 A) = 26.98, a factor 16 below the target: normalizing it to
Salz's L_X leaves a factor 5 discontinuity at 100 A (join 0.198) and 0.539 at
912 A. Both proxies reproduce log F_XUV to 0.005 dex by construction, but only
eps Eri does it without deforming its own X-ray spectrum. Spectral type also
favours eps Eri (K2V against Salz's K0-2V for HD 189733).

**Construction, measured.** Flux at Earth x (d/a)^2 = 4.567494e14 with
d = 3.212 pc (the proxy) and a = 0.031 AU, i.e. the same equal-surface-flux,
equal-radius convention as the WASP-52 b and HD 209458 b files. Then, in
Salz's own bands: lambda < 100 A x 0.9737515 to L_X/(4 pi a^2) = 5.600e3;
100-911.6 A x 0.7304595 to L_EUV/(4 pi a^2) = 1.507e4; the Lya line above its
local continuum in 1206-1226 A x 0.3487423 to L_Lya/(4 pi a^2) = 9.959e3;
longward of 911.6 A nothing but the geometric carry. 17 negative source rows
set to zero. Coverage 3.939229-4999.939229 A, 4997 rows, strictly increasing.

**Salz's 9.4 A Lya Gaussian is deliberately not reproduced.** The paper states
the width "is adapted to the resolution of our photoionization solver, which
is 6.1 Angstrom at the Lyman alpha line" and is "about a factor of 10 broader
than usual Lya stellar emission lines". It is a property of their solver, not
of the star, so the proxy's reconstructed intrinsic profile is kept and only
the line flux is normalized.

**Check against the paper.** F(1-912 A) = 2.0726e4 erg cm^-2 s^-1,
log = 4.3165, against Salz's log F_XUV = 4.32; the residual is the trapezoid
rule on the stored bins, since the two normalized bands sum to 4.3153 by
construction. The HD 97658 variant gives 4.3176.

**Open at the switch.** (1) The proxy's own reconstructed Lya luminosity,
measured, is log 28.89 against Salz's 28.43, so step 2 divides the eps Eri
line by 2.87; if the HD 189733 Lya reconstruction of Bourrier and Lecavelier
des Etangs (2013), which is Salz's source, is preferred over Salz's Table 3
rounding, that constant moves. (2) The file's
F(912-2600.5 A) = 1.841e5 erg cm^-2 s^-1, the band the He 2^3S metastable
photoionizes (4.767775 eV = 2600.5 A), rests entirely on the equal-radius
assumption and is not normalized to any measured quantity of HD 189733; the
HD 97658 variant gave a factor 1.42 more there when it was measured over
912-2583 A, which is the size of the proxy uncertainty in that band. (The
2583 A of that earlier measurement was a superseded value of the metastable
edge; 0.7 per cent of band width moves none of the ratios quoted here.)

The orbital distance is settled: 0.031 AU, Salz's Table 2 value, is the
distance of record for this planet, because its calculations are compared with
his simulation (user decision, 2026-09-06). `HD189733b/input.inp` carries it,
and so do both HD 189733 b files here.

### Bourrier et al. 2020 (MOVES III): the measurement-based alternative

**Paper.** Bourrier, Wheatley, Lecavelier des Etangs et al. (2020), MNRAS 493,
559, `2020MNRAS.493..559B`, read in `references/Bourrier_2020MNRAS_493_559.pdf`.

**What it does.** Section 4.4, p. 573: "Most of the stellar EUV spectrum is not
observable from Earth because of ISM absorption. We therefore reconstructed the
entire XUV spectrum up to 1600 Angstrom using the DEM retrieval technique
described in Louden, Wheatley & Briggs (2017). This reconstruction is based on
the quiescent X-ray flux for the coronal region (from the data shown in
Section 3), and on the flux derived for the intrinsic FUV stellar lines
(Section 4.2) for the chromosphere and transition region. The XUV spectra could
thus be reconstructed for Visits B to E. ... The generated XUV spectra for the
four visits are available online as machine readable tables." Section 6, p. 577:
"These synthetic spectra, which extend up to 1600 Angstrom, are available
online. ... We thus include our reconstructed profiles for the intrinsic Lyman
alpha line of HD 189733 in its synthetic spectra."

This is a reconstruction of *this star* from simultaneous X-ray (XMM-Newton
EPIC-pn, Swift) and HST/STIS FUV line measurements in five epochs, not a proxy
shape carrying scaling-relation luminosities. Table 3 (READ) gives, as total
fluxes at 1 au and with the bands stated in its own note ("in the synthetic EUV
spectra (62-912 A), and in the model X-ray spectra (0.2-2.4 keV = 5.2-62.0 A)"):

| Quantity, at 1 au | Visit B | Visit C | Visit D | Visit E | Average |
|---|---|---|---|---|---|
| Lyman alpha (erg cm^-2 s^-1) | 12.83 | 13.38 | 13.21 | 11.25 | 12.77 +/- 0.25 |
| X-ray 5.2-62.0 A | 4.86 | 6.54 | 7.08 | 6.92 | 6.35 |
| EUV 62-912 A | 18.14 | 24.70 | 22.53 | 20.76 | 21.53 |

(Visit A has no X-ray coverage. Fig. 19's caption says 62-920 A for the same
posterior; the Table 3 note's 62-912 A is used here.) Table 1 (READ):
D_star = 19.78 pc, R_star = 0.780 R_sun, a_p/R_star = 8.863, so
a = 0.032149 au (MEASURED from those two, the paper tabulates no a in au).

**Converted, MEASURED**, with L = F(1 au) x 4 pi (1 au)^2 = F x 2.81216e27:
log L_Lya = 28.555, log L_X(5.2-62 A) = 28.252, log L_EUV(62-912 A) = 28.782,
and F_XUV(5.2-912 A) at a = 0.032149 au = 2.697e4 erg cm^-2 s^-1,
log F_XUV = 4.431. The epoch-to-epoch spread of that last number is 0.13 dex
(Visit B 4.348 to Visit C 4.480), which is variability the file cannot carry
and which no single-epoch normalization can represent.

**Side by side with the two Salz papers** (log10, erg s^-1 for luminosities):

| Quantity | Salz et al. 2016 Table 3 | Bourrier et al. 2020 Table 3 avg | Salz et al. 2018 section 2 |
|---|---|---|---|
| log L_X | 28.18, band 0.124-2.48 keV = 5.0-100 A, from Sanz-Forcada et al. 2011 | 28.252, band 0.2-2.4 keV = 5.2-62.0 A, measured here | "approximately 2 x 10^28", i.e. 28.30, from Hunsch et al. 1999 |
| log L_EUV | 28.61, band 100-912 A, Linsky et al. 2014 relation from L_Lya | 28.782, band 62-912 A, DEM reconstruction | "3 x 10^28", i.e. 28.48, from Sanz-Forcada et al. 2011 |
| log L_Lya | 28.43, from Bourrier and Lecavelier des Etangs 2013 | 28.555, measured HST/STIS, five epochs | not quoted |
| log F_XUV at the planet | 4.32, <912 A, a = 0.031 au | 4.431, 5.2-912 A, a = 0.032149 au | not quoted |

Bourrier's L_X is 1.18x Salz 2016's *in a narrower band*, his EUV 1.48x, his
Lya 1.34x, and the resulting XUV at the planet 1.29x. Salz et al. 2018 quote
different literature values again, 1.3x and 0.74x Salz 2016's; the three papers
do not agree on this star to better than a factor 1.5 in any band.

**Salz et al. 2018 (A&A 620, A97) uses no SED at all.** Read in
`references/Salz_2018A&A_620_A97.pdf`: it is the He I 10830 detection paper and
runs no atmosphere model. Section 4.2, p. 8: "While we do not fit a proper
atmosphere model here, the derived He I triplet state column density is to be
understood as an effective value, which can be compared to theoretical models.
From the evaporation model of Oklopcic and Hirata (2018) for HD 209458 b, we
derived a weighted mean absorption height of 1.6 Rp with a column of about
7.9 x 10^11 cm^-2." Its only high-energy numbers are the two quoted literature
luminosities in the table above. So it constrains no SED, and it cannot be used
to choose one.

It does, however, carry a published endorsement of this file's proxy.
Section 4.1.1, p. 7: "Among the stellar sample studied by Andretta et al.
(2017), eps Eri shows properties comparable to HD 189733. In particular, this
active K dwarf shows values of 258 and 51 mAngstrom for the EWs of the stellar
infrared triplet components and 18.1 mAngstrom for the optical line along with
an X-ray luminosity of 2.1 x 10^28 erg s^-1." That is the same pairing this
file makes, made in the He I 10830 context, and at an X-ray luminosity within
0.1 dex of the "approximately 2 x 10^28" the same section gives for HD 189733.

### The released tables, and the file built from them

**Source.** `references/Bourrier2020_supplemental_files/`:
`HD189733_XUV_spectrum_Visit{B,C,D,E}.dat` and `readme_MNRAS.txt`. The readme
gives the two columns as "the wavelength, in angstroms" and "the stellar flux,
in erg/cm2/A/s", says the spectra "extend up to 1600 A, and are defined as a
function of wavelength in the stellar rest frame", and **names no reference
distance**. All four are on one identical grid, 0.016-1600.000 A in 0.016 A
steps, 99999 rows, no negative flux (measured).

**The reference distance, MEASURED rather than assumed.** Integrating each
released table over 62-912 A returns the EUV entry of the paper's Table 3,
which the table's own note defines as "the total flux at 1 au from the star":

| Visit | F(62-912 A) of the table | Table 3 | ratio | F(5.2-62 A) | Table 3 | ratio | F(1206-1226 A) | Table 3 Lya | ratio |
|---|---|---|---|---|---|---|---|---|---|
| B | 18.1417 | 18.14 | 1.000093 | 5.3362 | 4.86 | 1.0980 | 13.4005 | 12.83 | 1.0445 |
| C | 24.6999 | 24.70 | 0.999995 | 5.3657 | 6.54 | 0.8204 | 13.9414 | 13.38 | 1.0420 |
| D | 22.5320 | 22.53 | 1.000087 | 6.1010 | 7.08 | 0.8617 | 13.7729 | 13.21 | 1.0426 |
| E | 20.7596 | 20.76 | 0.999981 | 5.9961 | 6.92 | 0.8665 | 11.7941 | 11.25 | 1.0484 |
| mean | 21.5333 | 21.53 | 1.000152 | 5.6997 | 6.35 | 0.8976 | 13.2273 | 12.77 | 1.0358 |

(erg cm^-2 s^-1.) The EUV agreement to 1e-5 in every visit settles it: **the
tables are the flux at 1 au**, already in the units `read_sed` wants, so no
unit conversion is applied. The Lya column runs 4 per cent high because the
1206-1226 A window also contains Si III 1206.5 and O V 1218.3 A and the
continuum, not because the normalization differs. The X-ray column is the one
band that does not reproduce Table 3, by up to 18 per cent in either
direction: Table 3's X-ray entry is the APEC fit of section 4.3, while the
table carries the DEM output, and section 4.4 claims only that "the X-ray flux
was recovered with consistent values to those found in Section 4.3".

**The released tables reach further than Table 3 reports.** Measured on the
four-visit mean at 1 au, F(0.016-5.2 A) = 2.840 erg cm^-2 s^-1, 9 per cent of
F(0.016-912 A) = 30.074, peaking near 1.8 A. That band lies outside the
0.2-2.4 keV in which Table 3's X-ray flux is quoted, and it is also above
`e_top` = 1.24e3 eV, so `read_sed` discards everything below 10 A: the run
sees F(10-912 A) = 26.408 at 1 au, i.e. 2.7481e4 at the planet,
log = 4.4390 at a = 0.031 AU. The full F(1-912 A) of the file is 3.0782e4,
log = 4.4883.

**Visit versus average.** The shipped file is the arithmetic mean of the four
visits, which is the paper's own summary: Table 3's average EUV entry, 21.53,
is the mean of its four visit entries to four digits. (Its average Lya entry,
12.77 +/- 0.25, is over five visits; Visit A has no XUV spectrum, and the mean
of the four Lya entries that do is 12.6675.) The numbers of each visit are in the
table above; the epoch spread is 0.13 dex in log F_XUV, which is larger than
the gap between any two of the three papers' normalizations, and a single-epoch
file cannot represent it. Rebuilding on one visit is a change of one line in
the script.

**Scaling.** x (1 au / a)^2 = 1040.582726 with **a = 0.031 AU**. That is the
value of Salz et al. (2016) Table 2, adopted as the orbital distance of record
for HD 189733 b because this planet's calculations are compared with his
simulation (user decision, 2026-09-06; `HD189733b/input.inp` now carries
0.031). Bourrier's own Table 1 implies a = 0.032149 AU
(a_p/R_star = 8.863 times R_star = 0.780 R_sun; the paper tabulates no
semimajor axis in au), which would lower every flux in this file by 7.5 per
cent; `--a-au 0.032149` builds that variant.

**Longward of 1600 A.** The tables stop at 1600.000 A, far short of the 2600 A
the He 2^3S metastable needs, so 1600-5000 A carries the same measured K2V
photosphere as `hd189733b_epseri_salz2016.txt`: the MUSCLES v22 const-res SED
of eps Eri, flux at Earth x (d/a)^2 = 4.567494e14 with d = 3.212 pc at
a = 0.031 AU. Since both files now sit at that a, this is the *same* number
the Salz file uses, and measured, the 3400 photospheric rows of the two files
agree bit for bit (max relative difference 0.0e0 on every row above 1600 A).
**The two HD 189733 files differ only in the part Bourrier et al. measured**,
0.016-1600 A.

**The photospheric radius, and why no correction is applied.** Carrying the
proxy with (d/a)^2 gives HD 189733 eps Eri's surface flux *and* eps Eri's
radius. Putting that surface flux on the R_star = 0.780 R_sun of Bourrier's
Table 1 would multiply the photospheric part by (R_star/R_proxy)^2, and eps
Eri's radius is not a published number available here. It can be measured from
the MUSCLES product itself: the PHOENIX component's NORMFAC = 2.459320e-25 is
the factor carrying the BT-Settl model to the flux at Earth (verified here,
NORMFAC x model reproduces the released SED to 2 per cent at 8000 and
20000 A), and BT-Settl models are surface fluxes in erg s^-1 cm^-2 cm^-1, so
NORMFAC = (R/d)^2 x 1e-8 and R_proxy = 0.7065 R_sun at d = 3.212 pc, giving
(R_star/R_proxy)^2 = 1.219. **That chain rests on the unit convention of the
model grid, not on anything measured in this repository, so the default applies
no correction**; `--photo-scale 1.219` applies it. Applying it would also make
the two HD 189733 files differ in their photospheres, which is not what they
are for.

**The join at 1600 A is a real step and is not smoothed.** Measured,
f(first row above)/f(last row below) = 6.2614 and
F(1600-1650 A)/F(1550-1600 A) = 5.1315. Bourrier's reconstruction models the
corona, transition region and chromosphere and carries no photospheric
continuum, which near 1600 A is already the dominant term for a K dwarf; the
proxy's spectrum there is measured HST data and includes it. So the step is the
photosphere appearing, and it says the file understates the total flux over
roughly the last hundred Angstrom below the join. No absorber threshold lies in
that interval (the Lyman edge is at 912 A, the He 2^3S threshold at 2600 A), so
no photoionization rate is affected.

**No rebinning.** The released 0.016 A grid is kept as delivered, so no
flux-conservation question arises. The file is 103399 rows and 2.6 MB, the
largest here; `read_sed` reads it row by row and keeps the 100631 rows inside
[10.0, 2856.1] A.

**The shipped file against `hd189733b_epseri_salz2016.txt`** (erg cm^-2 s^-1,
measured on the two files):

both at a = 0.031 AU, so the comparison is of the spectra alone:

| Band | Salz file | Bourrier file | ratio |
|---|---|---|---|
| 1-100 A | 5.6286e3 | 1.0077e4 | 1.790 |
| 10-100 A (what the reader keeps) | 5.3318e3 | 6.7761e3 | 1.271 |
| 100-912 A | 1.5074e4 | 2.0704e4 | 1.374 |
| 1-912 A | 2.0726e4 | 3.0782e4 | 1.485 |
| 10-912 A | 2.0430e4 | 2.7481e4 | 1.345 |
| 5.2-62 A | 4.8802e3 | 5.9310e3 | 1.215 |
| 62-912 A | 1.5816e4 | 2.2407e4 | 1.417 |
| 912-1110 A (Lyman-Werner) | 2.6111e3 | 1.8413e3 | 0.705 |
| 1206-1226 A (Lya window) | 9.9968e3 | 1.3764e4 | 1.377 |
| 912-2600.5 A | 1.8409e5 | 1.8457e5 | 1.003 |
| 2600.5-5000 A | 6.7777e7 | 6.7777e7 | 1.000 |
| beta_m = F(1-100)/F(1-912) | 0.272 | 0.327 | |

The photosphere is identical, as it must be. The Bourrier file is harder and
brighter in the XUV, by 1.35 in the band the reader actually keeps, and
**fainter by 30 per cent in the 912-1110 A Lyman-Werner band**, which is the
band that drives H2 photodissociation in the molecular lower atmosphere; that
is the largest difference of the two in the direction of *less* flux and is
worth watching if this planet is ever run with molecules on.

**Recommendation.** `hd189733b_bourrier2020.txt` is the SED of record for
HD 189733 b: it is the authors' own reconstruction of this star from
simultaneous X-ray and FUV measurements, it carries the measured intrinsic Lya
profile rather than a Gaussian or a scaling relation, and it is the only one of
the three that quantifies the epoch-to-epoch variability.
`hd189733b_epseri_salz2016.txt` **stays as built and its normalization is
unchanged**, because the comparison of record for this planet is against Salz's
own simulated mass-loss rate (his Table 3, log Mdot_sim = 9.61), and moving
F_XUV by 1.35-1.49 would break that like-for-like. Two files, one per
comparison, on one orbital distance.

**An interim variant, measured but not shipped.** `build_hd189733b_epseri_salz.py`
now takes `--w-split` and `--w-xnorm-lo`, so the eps Eri shape can be
renormalized in Bourrier's bands with one command:

```
python3 build/build_hd189733b_epseri_salz.py --w-split 62.0 --w-xnorm-lo 5.2 \
    --log-lx 28.2519 --log-leuv 28.7822 --log-llya 28.5553 --a-au 0.032149 \
    --out <somewhere outside inputdata/sed>
```

MEASURED, that gives c_X = 1.3186, c_EUV = 1.0486, c_Lya = 0.4654,
F(1-912 A) = 2.705e4 (log F_XUV = 4.432), F(912-2583 A) = 1.593e5 (over the
superseded 2583 A metastable edge; the threshold of record is 2600.5 A),
joins 1.557
at 62 A and 0.656 at 912 A. The notable number is c_EUV = 1.049: eps Eri's own
EUV, carried to this orbit, reproduces Bourrier's DEM reconstruction of
HD 189733 to 5 per cent with no correction at all, which is a stronger
statement for the proxy than anything in the Salz normalization. Such a file
carries a `# NOTE:` block naming the luminosities actually used, so it cannot
be mistaken for the Salz build. **Defaults are unchanged and the shipped file
is byte-identical to the one built before these options were added (verified by
`diff`).**

**A rebinned copy of the Bourrier spectrum is already in the workspace, and the
earlier verdict on it was wrong.** `VULCAN/atm/stellar_flux/sflux-HD189_B2020.txt`
carries the header "binned and reconstructed from V. Bourrier et al. 2020
'MOVES III. Simultaneous X-ray and ultraviolet observations unveiling the
variable environment of the hot Jupiter HD 189733b'", 16.79-8000 A, in
erg cm^-2 s^-1 nm^-1. This README previously recorded that "its normalization
convention could not be established" and that "it is neither a surface flux nor
a flux at 1 au", on the strength of its 800 nm value against a blackbody. That
test was invalid: Bourrier's spectrum stops at 1600 A, so 800 nm is not his
data, and the file's own grid spacing changes from 1.6 A to about 0.5 A exactly
at 1599.2 / 1599.8 A (MEASURED), which is where it stops being his. With
Table 1's R_star = 0.780 R_sun in hand, the surface-flux hypothesis now checks
out: (1 au / R_star)^2 = 7.6001e4, and dividing the file by that gives, at
1 au, a Lya-window flux of 12.47 against Bourrier's 12.77 average (2.3 per cent
low) and an EUV(62-912 A) flux of 16.96 against his 21.53 average and 18.14 in
Visit B. Its X-ray band cannot be compared because the file starts at 16.79 A,
above the 5.2 A blue edge of Bourrier's band. **It is a stellar surface flux
carrying the Bourrier spectrum below 1600 A.** It is still not adopted, for a
different reason: it is a third-party rebinning, not the authors' released
table, its EUV sits 21 per cent below the published average with no stated
epoch, and above 1600 A it is a different spectrum altogether (its ratio to
pi B(5040 K) is 0.42, 0.30 and 0.24 at 4000, 6000 and 8000 A, MEASURED, so it
is a blanketed photosphere model, not a blackbody). Adopting it is a decision
to be taken deliberately, not by default.

**Note on the run.** `HD189733b/input.inp` runs the power law with
log L_X = 28.420 and log L_EUV = 28.680, which is 1.7x and 1.2x Salz's Table 3
values. At the 0.031 AU it now carries, that power law gives
F_XUV = 2.744e4 erg cm^-2 s^-1 against Salz's log F_XUV = 4.32 (2.09e4) and
this directory's Bourrier file's 3.078e4.

### `inputdata/scaled_solar_hd189.sed`: an orphan, not the SED of record

The file exists (1900 rows, 0.5-1899.5 A) and nothing in the tree references
it; `git log` puts it in one bulk commit, 586520e of 2026-06-14, with the
message "updates", alongside 23 unrelated files, and there is no build script
and no memo. What could be measured about it:

* its wavelength grid is the WHI convention exactly, 0.1 nm intervals on
  0.05 nm centers, so the parent is a solar reference spectrum of that family;
* it carries a solar Lyman alpha line (peak 4.15e4 erg cm^-2 s^-1 A^-1 at
  1215.5 A);
* **F(1-912 A) = 2.4224e4 erg cm^-2 s^-1 is the F_XUV that
  `HD189733b/input.inp` produced from its own `Log10 of X-ray luminosity`
  28.420 and `Log10 of EUV luminosity` 28.680 at the 0.0330 AU it carried until
  2026-09-06, 2.4217e4, to 0.03 per cent.** So it was normalized to that run's
  total XUV. That input now states 0.031 AU, so the match no longer holds;
* the X-ray/EUV split was **not** preserved: the file's
  beta_m = F(1-100)/F(1-912) = 0.089 against 0.355 for the same two
  luminosities. Its X-rays are four times too faint for its own normalization;
* it is not a constant multiple of any solar spectrum in the workspace. Against
  the WHI 2008 columns the ratio drifts from 4.7e3 near 5 A to 1.15e4 over
  1700-1900 A; against `Gueymard_solar.txt`, `VPL_solar.txt`,
  `sflux-HD189_Moses11.txt` and `sflux-HD189_B2020.txt` the scatter is a factor
  6 to 100. The parent spectrum is therefore **not identified**;
* it stops at 1899.5 A, so with `Include He23S? True` (which `HD189733b`,
  `benchmarks/hd189` and the regression cases all set) the reader stops the
  run, naming He 2^3S as the absorber left without a field. Measured, the two
  band lines of the message:

  ```
      band not covered     :  4.80000E+00 to  6.52720E+00 eV
                           =  1.89950E+03 to  2.58300E+03 A
  ```

**Verdict: orphan.** It is not the HD 189733 b SED of record: no paper is
attached to it, the shape's provenance cannot be recovered, its hardness
contradicts the luminosities it was normalized to, and it does not reach the
metastable band. It should not be adopted. Deleting it is a separate decision;
until then this section is its documentation.

---

## WASP-121 b: two candidates (BUILT)

**Paper.** Huang, Koskinen, Lavvas and Fossati (2023), ApJ 951, 123,
`2023ApJ...951..123H`.

**What the paper says.** Section 2.2: "In order to constrain the XUV spectrum
of the host star, we use the integrated XUV flux at the planet
F_XUV = 1.6 x 10^6 erg cm-2 s-1 from Salz et al. (2019), which is based on the
measured stellar rotation rate. The spectral energy distribution (SED) at
wavelengths shorter than 1700 A is generated based on the solar SED. The SED is
rescaled by the integrated XUV flux and the stellar radius considering that
WASP-121 is larger than the Sun. The shape of the spectrum is consistent with
observations of the star Procyon by the Extreme Ultraviolet Explorer (EUVE) but
the flux is about 10 times higher than that of Procyon in the same wavelength
range, which is expected given that Procyon is likely less active, while
WASP-121 is a fast rotator and therefore likely more active (Fossati et al.
2018). The SED at longer wavelengths is generated by LLmodels (Shulyak et al.
2004), which has been specifically designed for simulating the spectra of F-,
A-, and B-type stars. The spectra between 1210 and 1220 A have a 0.5 A
wavelength bin width, while the rest of the spectra have a 5 A bin width." And:
"The synthetic spectrum of WASP-121 is shown by the black line in Figure 1. The
average solar spectrum from Koskinen et al. (2013a), a 6459 K blackbody, and
the observed out-of-transit NUV spectrum (Sing et al. 2019) are shown for
comparison. We note that the Lyman continuum flux of the spectrum is
2.69 x 10^5 erg cm-2 s-1."

Salz et al. (2019) is `2019A&A...623A..57S`, "Swift UVOT near-UV transit
observations of WASP-121 b".

**What is missing.** Two pieces: the SOLAR2000 + Woods and Rottman average
solar spectrum Huang's sub-1700 A part is built on, and the LLmodels synthetic
spectrum for T_eff = 6459 K above 1700 A. Neither is in the workspace, so two
candidates are built, differing in what stands in for them, with the same
three normalizations taken from the paper.

**The three normalizations, all from section 2.2 quoted above.** (1) One
constant on everything shortward of 1700 A so that the integrated XUV flux at
the planet is F_XUV = 1.6e6 erg cm^-2 s^-1. (2) The Lya line above its local
continuum in the window 1206-1226 A scaled so that the line flux at the planet
is 1.0e5 erg cm^-2 s^-1. (3) One constant longward of 1700 A so that the
bolometric flux at the planet is sigma T_eff^4 (R_star/a)^2 = 7.003e9
erg cm^-2 s^-1 with T_eff = 6459 K (the blackbody of their Figure 1),
R_star = 1.4572 R_sun and a = 0.02544 au (their Table 1). Huang's "rescaled by
the integrated XUV flux and the stellar radius" is contained in (1) and (3);
no distance or radius of a proxy star enters anywhere.

**Candidate A, `wasp121b_solar_huang2023.txt` (recommended).** Shortward of
1700 A the solar SED is the WHI 2008 reference spectrum, the same table the
HD 209458 b file is built from, x 5.387291e5. Longward of 1700 A, LLmodels is
replaced by pi B_lambda(6459 K) x (R_star/a)^2 = 7.095733e-2, which is the
blackbody Huang et al. plot in their own Figure 1 as the comparison curve for
exactly this band, and which is also Salz's step 4 for the same purpose.
5000 rows, 0.5-4999.5 A. The Lya constant is 2.961511e-2.

**Candidate B, `wasp121b_wasp17_huang2023.txt`.** The MUSCLES v24 const-res
broadband SED of WASP-17 (F6V), the closest match in spectral type to
WASP-121 among the 36 MUSCLES / Mega-MUSCLES stars held locally, carries both
shapes: x 7.450159e20 shortward of 1700 A and x 1.206591e19 longward of it,
the latter from the BOLOFLUX = 5.803769e-10 erg cm^-2 s^-1 of its FITS primary
header. 4991 rows, 10.0-5000.0 A, 4 negative source rows set to zero. Its
value is that longward of 1700 A it is HST/STIS G230L and G430L data plus a
PHOENIX photosphere, so it carries the UV line blanketing a blackbody lacks.

**Which one reproduces the paper's numbers.** Measured, against the three band
fluxes Huang et al. state for their own synthetic spectrum (erg cm^-2 s^-1):

| Band | Huang et al. | Candidate A | ratio | Candidate B | ratio |
|---|---|---|---|---|---|
| 700-800 A | 1.4e4 | 6.840e4 | 4.89 | 1.298e5 | 9.27 |
| 800-912 A | 3.6e4 | 1.740e5 | 4.83 | 3.276e5 | 9.10 |
| 912-1170 A | 6.3e4 | 3.002e5 | 4.76 | 2.566e5 | 4.07 |

Candidate A is high by one constant factor, 4.8 in all three bands, so its
*shape* over 700-1170 A is Huang's to 3 per cent. Candidate C carries
candidate A's XUV unchanged, so its row of this table is A's: 6.840e4, 1.740e5
and 3.002e5. Candidate B is high by 9.2,
9.1 and 4.1, so its shape is not. Both reproduce F_XUV = 1.6e6 exactly, by
construction, and both put the Lya line flux at 1.0e5.

**A tension inside the paper, stated and not resolved here.** Section 2.2 says
both that the SED is normalized to F_XUV = 1.6e6 and that "the Lyman continuum
flux of the spectrum is 2.69 x 10^5 erg cm-2 s-1". Those two cannot both hold
for a solar-shaped SED: rebuilt with `--fxuv 2.69e5`, candidate A gives
1.150e4, 2.925e4 and 5.048e4 in the three bands, within 20-25 per cent of
Huang's 1.4e4, 3.6e4 and 6.3e4 and in the same direction, i.e. the quoted band
fluxes are the ones consistent with F(1-912 A) = 2.69e5, not with 1.6e6. The
factor between the two readings is 5.95. Pulling the other way, the join at
Huang's own 1700 A is continuous only under 1.6e6: measured, candidate A's
step there is f(above)/f(below) = 1.1794, and it would be 7.015 under 2.69e5.
Both files are therefore built to the paper's explicit normalization
instruction, 1.6e6, which is also the F_XUV that `WASP-121b/input.inp`
produces from its power law today (1.604e6 at 0.02544 AU from log L_X = 29.46
and log L_EUV = 30.42); `--fxuv 2.69e5` rebuilds either file under the other
reading. Settled 2026-09-06 with the published Salz et al. (2019), A&A 623,
A57 (`references/Salz_2019A&A_623_A57.pdf`): its introduction defines "XUV,
lambda < 912 A", and section 7 derives "log10 FXUV (erg cm-2 s-1) = 6.2" for
WASP-121 b from the canonical log10 LX = 30.0 of a 1 d rotator (Pizzolato et
al. 2003) through the X-ray-EUV relations of King et al. (2018). So the
1.6e6 IS the whole 1-912 A band and the normalization of both candidates is
the one Huang et al. prescribe; the 2.69e5 sentence describes some other
quantity of their spectrum and is not used.

**Candidate C, `wasp121b_composite_huang2023.txt` (provisional composite,
BUILT 2026-09-06 by `build/build_wasp121b_composite.py`).** Candidate A below
1700 A and candidate B above it, each component carried unchanged, row for
row: x 5.387291e5 on the WHI 2008 solar SED with the Lya constant 2.961511e-2
below the join, x 1.206591e19 on the WASP-17 MUSCLES v24 SED above it. 5001
rows, 0.5-5000.0 A, 1 A bin averages throughout. The build script imports the
reader and the band integral of `build_wasp121b_solar_huang.py`, and the
coverage script asserts that all 1700 rows below the join are candidate A's
and all 3301 rows above it are candidate B's, exactly.

*Neither component is renormalized to force continuity at 1700 A.* Each keeps
the published normalization of its own band, F_XUV at the planet below the
join and the bolometric flux at the planet above it; scaling either to close
the step would discard the normalization the paper prescribes for it. The step
is measured and stated instead: f(1700.0 A)/f(1699.5 A) = 0.1324, and
F(1700-1750 A)/F(1650-1700 A) = 0.4272 (4.8313e5 / 1.1310e6). Candidate A's
step at the same place is 1.1794 upward. The composite therefore has a real
discontinuity at the join, and it is the honest one: the solar FUV continuum
just below 1700 A, lifted to WASP-121's XUV, is 7.6x the F6V photosphere just
above it. Both readings cannot be right, and nothing in the workspace settles
which; the header says so.

*The header says what the file is:* provisional composite, solar XUV shape
below 1700 A at the F_XUV = 1.6e6 normalization Huang et al. prescribe (the
band definition is Salz et al. 2019), WASP-17 photosphere above 1700 A as a
proxy for WASP-121 whose T_eff, gravity and line blanketing are the proxy's
and not the target's, **not a reproduction of the Huang et al. spectrum**
(their LLmodels run above 1700 A is unavailable).

**Sensitivity of the production rates to the candidate input.** Measured
2026-09-06, loading each file through the production `set_energy_vectors` with the
configuration of `WASP-121b/input.inp` and its `metals.inp` (He 2^3S on,
Stellar Teff 6459 K and Stellar radius 1.458 R_sun, Stellar Lya flux 1.0e5,
solar C/N/O/Mg/Ca/Na/Fe, a = 0.02544 AU, [13.60, 123.98, 1.24e3] eV, X-rays
on, dayside dilution 0.5), on ONE prescribed unattenuated field at the orbit:
no atmosphere, no column, no optical depth. The four photoionization rates are
formed as `util_ion_eq` forms them; the Balmer rate and heating are
`gamma_n2_balmer` and `heat_n2_balmer` of `excited_hydrogen`; the band fluxes
are the prescription of `lyman_werner_band_flux_from_sed`, checked against
that function itself on all three files.

**These are differences between candidate inputs, not a bias of any
candidate.** A bias would need a justified reference field for WASP-121, and
the workspace has none.

| Quantity | A | B | C | C/A | C/B |
|---|---|---|---|---|---|
| `P_HI` [s^-1] | 3.308690e-2 | 5.294887e-2 | 3.308690e-2 | 1.0000 | 0.6249 |
| `P_HeI` [s^-1] | 2.081327e-2 | 1.631648e-2 | 2.081327e-2 | 1.0000 | 1.2756 |
| `P_HeII` [s^-1] | 1.708502e-3 | 1.153713e-3 | 1.708502e-3 | 1.0000 | 1.4809 |
| `P_HeTR` [s^-1] | 4.638721e1 | 2.710192e1 | 2.489033e1 | 0.5366 | 0.9184 |
| `P_HeTR`, part below the 1700 A join | 4.595538e1 | 2.445920e1 | 2.445967e1 | 0.5322 | 1.0000 |
| `P_HeTR`, part above the join | 4.318252e-1 | 2.642711 | 4.306575e-1 | 0.9973 | 0.1630 |
| Balmer rate [s^-1] | 6.143432e2 | 4.738348e2 | 4.729875e2 | 0.7699 | 0.9982 |
| Balmer heating [erg s^-1] | 4.839797e-10 | 3.465302e-10 | 3.406892e-10 | 0.7039 | 0.9831 |
| F(LW 912-1110 A) | 2.539392e5 | 1.967098e5 | 2.539392e5 | 1.0000 | 1.2909 |
| F(B1 1110-1201 A) | 1.387185e5 | 1.898337e5 | 1.387185e5 | 1.0000 | 0.7307 |
| F(B2 1206-1226 A) | 1.476257e5 | 1.567026e5 | 1.476257e5 | 1.0000 | 0.9421 |
| F(B3 1231-1450 A) | 5.326850e5 | 2.087076e6 | 5.326850e5 | 1.0000 | 0.2552 |
| F(B4 1451-2304 A) | 8.159480e7 | 7.588232e7 | 5.973456e7 | 0.7321 | 0.7872 |

What the table says, item by item.

* The three ground-state rates of C are A's to machine zero, because their
  thresholds are 13.598, 24.587 and 54.418 eV and the integrands live below
  912 A, where C carries A's rows. Against B they move by factors 0.62, 1.28
  and 1.48: the solar and the WASP-17 XUV shapes are not the same spectrum
  even at the same F_XUV.
* `P_HeTR` is **not** a photosphere-only quantity. Its threshold is
  4.767775 eV, so the integrand runs from 2600.5 A down to the 10 A top of the
  grid. Split at the join, C's part below it is B's to 1.9e-5 and C's part
  above it is A's to 2.7e-3; neither is exact, because the one bin straddling
  1700 A has different edges in the three files (the geometric-mean edge is
  built from the last row below the join and the first row above it, and those
  pairs are 1699.5/1700.5 in A, 1699.0/1700.0 in B and 1699.5/1700.0 in C).
  The full rate of C sits 8 per cent below B's, all of it in the XUV part.
* Between A and C the metastable rate moves by 0.54, which is the same factor
  as F(1700-2600.5 A): 1.7985e8 for A against 9.9033e7 for B and C. That is
  the whole size of the photosphere choice for the He 10830 problem.
* The Balmer continuum reads its own field over 3.400-13.598 eV, i.e.
  912-3647 A, which straddles the join. C's rate and heating are B's to 0.2
  and 1.7 per cent and A's to 0.77 and 0.70.
* The band fluxes are the spectrum integrated over each band. What the
  equilibrium solve uses is `fuv_band_flux(ib)`, which returns the input keys
  (`Stellar LW / Lya / B3 / B4 flux`), and only the LW band is taken
  from the spectrum file, and only in a molecular run. B3 and B4 are
  therefore identical for the three candidates as `WASP-121b/input.inp`
  stands, and the columns above are what a molecular run or a future
  band-from-spectrum route would see.

**Recommendation for the switch: candidate C**, pending the user's C3
decision. It takes the XUV from the component the paper literally prescribes
and the 1700-2600.5 A band from a measured F-star photosphere with its line
blanketing, which is the band the He 2^3S metastable photoionizes out of. It
is a provisional composite, not a reproduction of Huang's spectrum, and the
1700 A step is a stated property of it, not a defect that was hidden.

**Between the two single-source candidates, A**, on three counts: it is what
section 2.2 literally prescribes shortward of 1700 A; its 700-1170 A shape reproduces
Huang's to 3 per cent while candidate B's does not; and it joins the
photosphere at 1700 A within 18 per cent, while candidate B's XUV
normalization sits a factor 61 above its own photospheric level and leaves a
factor 61 step at 1700 A, so its 912-1700 A FUV continuum is not physical.

**What candidate A costs, and what candidate B is for.** Measured,
F(1700-2600.5 A) is 1.7943e8 for candidate A against 9.898e7 for candidate B,
so the blackbody is 1.81x brighter than the measured F6V photosphere in
exactly the band the He 2^3S metastable photoionizes (4.767775 eV = 2600.5 A);
over 2600.5-5000 A, where line blanketing is weak, the two agree to 0.02 per
cent (2.0866e9 against 2.0869e9),
which is the check that the bolometric normalization is doing its job. The
He 2^3S photoionization rate is 1.8x higher on candidate A than on candidate B
(4.639e1 against 2.710e1 s^-1, measured); that is the sensitivity of the rate
to the photosphere assumed, not a bias of either file, which would need a
reference field for WASP-121. The construction that takes the XUV from A and
the photosphere from B is candidate C above.

**What to ask for, still.** The synthetic WASP-121 spectrum of Huang et al.
(2023) Figure 1 (black line), from the authors; failing that, an LLmodels run
at T_eff = 6459 K with the log g and [Fe/H] of WASP-121 for the NUV, plus the
Koskinen et al. (2013a) average solar spectrum for the XUV, plus Salz et al.
(2019) for the F_XUV band definition.

**Note on the run.** The power law of `WASP-121b/input.inp` already puts
F_XUV = 1.604e6 erg cm^-2 s^-1 at 0.02544 AU (log L_X = 29.46,
log L_EUV = 30.42), i.e. Huang's normalization to three digits, and its orbital
distance already matches Huang's Table 1. What the switch
would change is the shape, not the total: an index -1 power law against a solar
SED, and the sub-13.6 eV band that decision 13 now defines as part of the same
power law. `benchmarks/wasp121/input.inp` carries `Include He23S? False` and so
would not need the 2600 A floor; `WASP-121b/input.inp` carries `True` and would.

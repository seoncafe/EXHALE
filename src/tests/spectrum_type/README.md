# spectrum_type

Tests of decision 13 of `docs/development_plan_20260905_rev3.md` section 10.5:
**one spectrum type builds every band of the photon grid.** `Spectrum type:`
in `input.inp` states the spectrum of the whole grid, the XUV and the band
below 13.6 eV alike, which is the band where the He 2^3S metastable
(4.768 eV), the low-IP metals and hydrogen in n = 2 (3.400 eV) absorb.
`Planck` is the photospheric-blackbody type added by that decision.

Run everything with

    src/tests/spectrum_type/run.sh [planck_field] [balmer_field] \
                                  [balmer_quad] [n2_floor] [sed_edges] \
                                  [sed_semantics] [wasp121_sed] \
                                  [spectrum_gate]

`EXHALE_OBJDIR` and `EXHALE_EXE` point the suite at another build (a private
`OBJDIR`/`EXE` pair); with either of them set the `make -q` staleness check is
skipped. Output is one `PASS|FAIL <name> measured=... reference=... tol=...`
line per assertion; the suite exits nonzero if any of them fails.

## planck_field (`planck_photon_field.f90`)

Links the production objects and calls the production `set_energy_vectors`
twice on the stellar and spectral input of `backup/regression/wasp_full`
(T_eff 6459 K, R_star 1.458 Rsun, a 0.02544 AU, power-law index -1 on
[13.60, 123.98, 1.24e3] eV, dilution 1/2, He 2^3S on): once as `Power-law`,
once as `Planck`.

| assertion | reference | tolerance |
|---|---|---|
| `planck_grid_identical_to_power_law` | the power law's own `e_v`, `de_v` | 0 (identity) |
| `planck_field_on_grid` | `xi pi B_nu(T_eff) (R_star/a)^2 dnu/dE`, re-evaluated in the driver | 1e-12 relative |
| `planck_flux_stefan_boltzmann` | `sigma T_eff^4 (R_star/a)^2`, sigma = 5.670374419e-5 (CODATA 2018) | 1e-6 relative |
| `planck_field_admissible` | every grid point finite and >= 0 (the 1.24 keV top of the grid is `h nu / k T = 2228`) | 0 |

It also prints, as a DIAGNOSTIC, the two types side by side at four energies.

## balmer_field (`balmer_continuum_field.f90`)

Links the production objects and calls the production `gamma_n2_balmer` and
`heat_n2_balmer`, the two scalar rates of the H(n=2) Balmer continuum
(3.40-13.6 eV), on the same stellar and spectral input, under each spectrum
type. Decision 13 requires them to read the run's own type; until item
2c-BALMER they were a photospheric blackbody in every run, power-law runs
included.

| assertion | reference | tolerance |
|---|---|---|
| `balmer_rate_power_law`, `balmer_heat_power_law` | the same two integrals of the production `J_inc`, on 40001 uniform points instead of the production 400 | 1e-3 relative |
| `balmer_rate_planck_frequency_form`, `balmer_heat_planck_frequency_form` | the frequency writing of the same band integral, `INT F_nu/(h nu) sigma_2 dnu`, on 400001 uniform points. `F_E dE = F_nu dnu`, so the two writings are one number and what is left is the reference's own quadrature error: 1e-6 relative, MEASURED 2.3e-10 | 1e-6 relative |
| `loaded_table_field_interpolated` | a synthetic power-law table on its own rows (`e_sed_node`/`F_sed_node`): the field returned between two rows is that power law | 1e-12 relative |
| `balmer_band_covered_power_law`, `balmer_band_covered_planck`, `balmer_band_uncovered_monochromatic`, `balmer_band_uncovered_short_table` | `spectrum_covers_eV` at the n=2 edge; the last two are the configurations in which the production functions stop the run | 0 |
| `balmer_band_covered_table_reaching_floor`, `loaded_table_field_at_the_grid_floor` | a table read down to the grid floor whose lowest selected row sits a fraction of a row above it: the floor is covered, and the field there is the table's own segment continued | 0, and 1e-12 relative |

It prints, as a DIAGNOSTIC, the n=2 photoionization rate under the two
analytic types side by side: 0.41857 s^-1 for the power law of
`backup/regression/wasp_full` against 1230.6 s^-1 for its blackbody, a factor
2940.

### State

RED before item 2c-BALMER (the driver does not build: `stellar_flux_eV` and
`spectrum_covers_eV` do not exist and the two production functions take a
temperature and a dilution factor), GREEN after it.

## balmer_quad (`balmer_band_quadrature.f90`)

Links the production objects and measures the QUADRATURE of the same two
functions: `gamma_n2_balmer` and `heat_n2_balmer` against the exact integral
of the very field they read, `stellar_flux_eV`. What field a loaded table
states is the question of `sed_semantics` below; here the field is taken as
given and only the rule that integrates it is under test.

The reference is 5-point Gauss-Legendre on 16 geometric sub-intervals of each
interval between consecutive NODES of the field (the table's own rows inside
the band, plus the two band edges; a geometric subdivision of the band for the
analytic types), which is smooth inside one such interval whatever the field
does across the band. Doubling the subdivision is asserted to leave it
unchanged, at 1e-12. The 400001-point uniform trapezoid is formed as well and
printed beside it, so the accuracy of that reference is on the record too
(MEASURED 7.4e-10 on the narrow-line table, 6.3e-7 on the eps Eri file).

| configuration | field |
|---|---|
| `narrow_line`, `narrow_line_h2`, `narrow_line_h4` | the bin-averaged synthetic table of `sed_semantics`, continuum plus two emission lines narrower than the table step (200.0 A and 1215.67 A, the second one inside this band), at the 1 A step of the shipped products and at a half and a quarter of it |
| `wasp52` | `inputdata/sed/wasp52b_epseri_yan2022.txt`, which carries a measured Ly-alpha line inside the band |
| `power_law` | `Spectrum type: Power-law`, index -1, the field of `backup/regression/wasp_full` |
| `planck` | `Spectrum type: Planck`, the 5000 K photosphere of `benchmarks/wasp52` at 0.0272 AU |

| assertion | reference | tolerance |
|---|---|---|
| `<config>_rate_quadrature`, `<config>_heat_quadrature` | the node-following reference above | 1e-6 relative |
| `<config>_reference_converged` | the reference against itself at twice the subdivision | 1e-12 relative |
| `narrow_line_quadrature_converges_rate`, `..._heat` | the deviation at a quarter of the table step is not above the deviation at the full step | 0 |

The deviation of the FIXED 400-point uniform rule from the same reference is
printed as a DIAGNOSTIC for every configuration. MEASURED: 3.5e-4 (rate) and
1.4e-3 (heating) on the narrow-line table, GROWING to 2.2e-3 and 2.0e-2 as
that table is refined, because the interpolated line narrows while the step
does not move; 2.9e-2 and 2.7e-2 on the eps Eri file; 9.5e-5 and 6.0e-5 on the
power law; 3.7e-4 and 3.6e-4 on the Planck field.

### State

RED before item BALMERQ of `docs/PLAN_20260906_rev2.md` (the production
functions were that 400-point rule: 14 of the 21 assertions fail, the three
narrow-line resolutions and both convergence statements among them), GREEN
after it, every deviation at 1e-14.

## n2_floor (`photon_grid_n2_floor.f90`)

Links the production objects and calls the production `set_energy_vectors`
twice on the same stellar and spectral input, He 2^3S on and solar Na
present: once with the excited-hydrogen coupling armed off, once with it on.
Item 2c-N2FLOOR makes the n = 2 threshold `e_th_HI_n2 = e_th_HI/4 =
3.400 eV` (3647 A) floor the photon grid whenever that coupling is armed, so
that a loaded SED is read down to it; this is the statement that doing so
costs the analytic types nothing, because the band it opens carries no cross
section the grid integrates.

| assertion | reference | tolerance |
|---|---|---|
| `n2_floor_grid_reaches_the_n2_edge` | `sum(de_v)` against `e_top - e_th_HI_n2` | 1e-12 relative |
| `n2_floor_added_bin_count` | `num_n2` = 20 of `set_energy_vectors` | 0 |
| `n2_floor_bins_above_the_absorber_edge_identical` | `e_v`, `de_v`, `F_XUV` of the armed-off grid | 0 (exact equality) |
| `n2_floor_cross_sections_vanish_below_the_absorber_edge` | 0, for H, He, He 2^3S, H2 and every ion of an element `metals.inp` carries | 0 |
| `n2_floor_photoionization_rates_unchanged`, `n2_floor_metal_rates_unchanged` | the same `sum(F sigma/E de_v)` on the armed-off grid | 0 (exact equality) |

The metal cross-section table is filled for all 17 photo-ionizable ions
whatever `metals.inp` carries, so K I (4.341 eV) is non-zero in it below the
metastable edge; it is inert unless K is present, and if K were present it
would be the lowest bulk-gas threshold and the band split would move down to
it. The assertions therefore run over the ions of the active elements.

### State

RED before item 2c-N2FLOOR on the first two assertions (the grid stopped at
the metastable edge, so it added no bin and its width was
`e_top - e_th_HeTR`), GREEN after it.

## sed_edges (`loaded_sed_threshold_edges.f90`)

Links the production objects and calls the production `set_energy_vectors` on
the WASP-52 b configuration of `benchmarks/wasp52/input.inp` with the eps Eri
table `inputdata/sed/wasp52b_epseri_yan2022.txt` (a 0.0272 AU,
[13.60, 123.98, 1.24e3] eV, X-rays on, He 2^3S on, stellar Teff and radius
set, so the grid floor is the H(n=2) edge 3.400 eV), once as shipped and once
with solar C, N, O, Mg, Ca, Na and Fe added so that the metal thresholds are
exercised too.

The grid a loaded table is put on: rows `E_1 < ... < E_N`, interior edges the
geometric mean `sqrt(E_k E_k+1)`, end edges `E_1` and `E_N`, each row's flux
over the whole of its bin, and every active ionization threshold strictly
inside the span inserted as a further edge. Splitting a bin at `t` and giving
both halves the same flux leaves `F(t-b) + F(b'-t) = F(b'-b)`, so the split is
exact for the integrated flux.

| assertion | reference | tolerance |
|---|---|---|
| `*_thresholds_are_bin_edges` | the nearest bin edge, rebuilt from `E_1` and the cumulative `de_v` | 1e-12 relative |
| `*_grid_spans_the_table` | `sum(de_v)` against `E_N - E_1` | 1e-12 relative |
| `*_integrated_flux_preserved` | `sum(F_XUV de_v)` against the histogram integral of the table over its own span, diluted | 1e-12 relative |
| `*_P_HI`, `*_P_HeI`, `*_P_HeII`, `*_P_HeTR` | the same integral of the same table over `[threshold, E_N]`: flux constant over each table bin, cross section on 1000 midpoints of each clipped bin | 1e-3 relative |
| `*_metal_rates` | the same, for every photo-ionizable ion of an element the run carries | 1e-3 relative |

It prints, as a DIAGNOSTIC, the deviation of each rate from that reference on
the PREVIOUS loaded grid (table rows as bin points, central-difference widths,
no threshold an edge, reproduced in the driver) beside the deviation of the
grid that replaces it: on the eps Eri table, P_HI 1.77e-3 -> 4.6e-7, P_HeI
1.10e-3 -> 1.1e-6, P_HeII 8.39e-3 -> 8.5e-6, P_HeTR 9.58e-5 -> 6.5e-8, and
over the thirteen metal ions the worst 3.22e-3 (O II) -> 5.5e-4 (Mg II). What
is left is the rectangle-rule error of the table's own 1 A rows, not a
straddled threshold.

### State

RED before item 2c-SEDGRID: the driver does not build (`read_sed` took no
argument and the table rows were not kept in `e_sed_node`/`F_sed_node`), and
the grid it asserts on did not exist. GREEN after it.

## sed_semantics (`loaded_sed_table_semantics.f90`)

What a ROW of a loaded table means, and what the production consumers of the
loaded field make of it. `sed_edges` above tests the layout of the grid;
this one tests the RECONSTRUCTION: is the table a point sample of `F_lambda`
or a bin average, does the reader's histogram-in-energy reproduce it, and do
the rates converge as the table is refined (review finding 8.4 of
`docs/To_be_determined_by_user_recommend_20260906.md`).

Four synthetic tables are built from an analytic `F_lambda` by exact bin
averaging, so their declared semantics is known, and each is read by the
production `set_energy_vectors` at the 1 A step of the shipped products and at
a half and a quarter of it:

| table | what it carries |
|---|---|
| `narrow_line` | two emission lines narrower than the table step, at 200.0 A and at 1215.67 A |
| `step_edge` | a continuum ten times brighter shortward of the He I threshold, a discontinuity inside a bin |
| `coarse_bin` | one bin ten times the width of its neighbours, across the H I edge |
| `steep_power` | `F_lambda ~ lambda^4`, i.e. `F_E ~ E^-6` |

Nine consumers are formed exactly as production forms them: the band energy
and photon number of the grid, `P_HI`, `P_HeI`, `P_HeII`, `P_HeTR` and the H I
photoheating integral of `util_ion_eq`, and the Balmer rate and heating of
`excited_hydrogen`, which read `stellar_flux_eV` instead of the grid. Two
references are formed for each: the HISTOGRAM the table declares, and the
EXACT integral of the analytic spectrum it was sampled from.

| assertion | reference | tolerance |
|---|---|---|
| `<table>_selected_span_is_table_rows`, `<table>_selected_row_count` | the rows integrated are the rows the reader kept | 1e-7, 0 |
| `<table>_band_energy` ... `<table>_heating` | the HISTOGRAM reference of the same span | 1e-3 relative |
| `<table>_balmer_rate_quadrature`, `<table>_balmer_heat_quadrature` | the same interpolated field on 40001 points instead of the production 400 | 1e-3 relative |
| `<table>_converges_<quantity>` | the deviation from the sampled spectrum at a quarter of the table step is below the deviation at the full step | strict |

Measured on the 2026-09-06 build: the seven grid-integrated consumers sit
between 1.4e-7 and 8.2e-5 of the declared histogram on all four tables, and
their deviation from the sampled spectrum falls at order 2 (band energy on the
narrow-line table 5.7e-7, 1.4e-7, 3.0e-8). The log-log interpolation the
Balmer band reads is exact on a power-law segment, so on the smooth tables the
refined Balmer integrals sit 2.6e-8 from the sampled spectrum.

**Two production findings, reported and not fixed here** (neither file is in
this item's scope):

* `excited_hydrogen.f90` integrates the Balmer band on 400 points uniform in
  energy, a step of 0.0255 eV, and cannot resolve a line narrower than that.
  On the narrow-line table the rate is 3.5e-4 and the heating 1.4e-3 from the
  same field on 40001 points, and the gap GROWS as the table is refined (rate
  3.5e-4, 9.7e-4, 2.2e-3): a finer table makes the interpolated line narrower
  while the rule's step does not move. A measured stellar spectrum carries
  such a line at Ly-alpha. The quadrature assertion is therefore made on the
  three tables the rule can resolve, and the narrow-line numbers are printed
  as a DIAGNOSTIC labelled as an open defect.
* `sed_read.f90` converts wavelength to energy with the SINGLE-PRECISION
  literal `1e-8` at its two `hp_eV*c_light/(w*1e-8)` sites, which shifts every
  photon energy of a loaded SED by 6.077e-9 relative. That is the residual the
  `_selected_span_is_table_rows` assertion measures.

### State

The driver is new with item C1 of `docs/PLAN_20260906_rev2.md`; before it
there was no measurement of the loaded-table reconstruction through the
production consumers, and the two findings above were unrecorded.

## wasp121_sed (`wasp121b_sed_sensitivity.f90`)

Item C2 of `docs/PLAN_20260906_rev2.md`. Loads the three WASP-121 b spectrum
candidates through the production `set_energy_vectors` with the configuration
of `WASP-121b/input.inp` and its `metals.inp` (He 2^3S on, Stellar Teff
6459 K, Stellar radius 1.458 R_sun, Stellar Lya flux 1.0e5, solar
C/N/O/Mg/Ca/Na/Fe, a = 0.02544 AU, [13.60, 123.98, 1.24e3] eV, X-rays on,
dilution 1/2), on ONE prescribed unattenuated field at the orbit: no
atmosphere, no column, no optical depth. Candidate A is
`wasp121b_solar_huang2023.txt`, B is `wasp121b_wasp17_huang2023.txt`, C is the
provisional composite `wasp121b_composite_huang2023.txt` (A's solar XUV below
1700 A, B's WASP-17 photosphere above it).

The four photoionization rates are formed as `util_ion_eq` forms them; the
Balmer rate and heating are `gamma_n2_balmer` and `heat_n2_balmer`; the band
fluxes use the prescription of `lyman_werner_band_flux_from_sed`. Everything
is printed side by side with the C/A and C/B ratios as a **sensitivity of the
production rates to the candidate input**. It is not a bias of any candidate:
a bias would need a justified reference field for WASP-121, which the
workspace does not have.

| assertion | reference | tolerance |
|---|---|---|
| `C_equals_A_P_HI`, `_P_HeI`, `_P_HeII` | candidate A's rates: the thresholds are 13.598, 24.587 and 54.418 eV, so the integrands live below 912 A, where C carries A's rows. MEASURED 0 | 1e-6 relative |
| `C_P_HeTR_is_B_below_join_plus_A_above` | B's part of `P_HeTR` below the 1700 A join plus A's part above it. `P_HeTR` is not photosphere-only: its threshold is 4.767775 eV, so the integrand reaches from 2600.5 A to the 10 A grid top. MEASURED 2.8e-5 | 1e-3 relative |
| `C_equals_B_P_HeTR_below_the_join` | B's part below the join. MEASURED 1.9e-5 | 1e-3 relative |
| `C_equals_A_P_HeTR_above_the_join` | A's part above the join. MEASURED 2.7e-3; the departure is the one bin straddling 1700 A, whose geometric-mean edges come from the row pairs 1699.5/1700.5 in A, 1699.0/1700.0 in B and 1699.5/1700.0 in C | 1e-2 relative |
| `C_F_1_912_is_the_stated_F_XUV` | 1.6e6 erg cm^-2 s^-1, the normalization the header states (Huang et al. 2023 section 2.2 with the band definition of Salz et al. 2019). MEASURED 1.8e-4 | 1e-3 relative |
| `C_F_10_912_equals_A` | A's F(10-912 A): the whole XUV component is A's. MEASURED 0 | 1e-12 relative |
| `C_F_1700_2600_equals_B` | B's F(1700-2600.5 A): the whole metastable band is B's. MEASURED 0 | 1e-12 relative |
| `band_integral_is_production_LW_A/B/C` | the production `lyman_werner_band_flux_from_sed()` on the same file, so the band fluxes printed are the production prescription and not a second one. The LW band is 912-1201 A on both sides since 2026-09-06. MEASURED 0 | 1e-12 relative |

The measured table is reproduced in `inputdata/sed/README.md`, WASP-121 b
section. The headline numbers: C's ground-state rates are A's exactly and move
by 0.62, 1.28 and 1.48 against B; C's `P_HeTR` is 0.54 of A's, which is the
size of the photosphere choice for the He 10830 problem, and 0.92 of B's, the
8 per cent difference sitting entirely in the XUV part of that rate.

### State

RED before the change of item C2: the composite file does not exist, so
`set_energy_vectors` stops in `read_sed`'s `open`. GREEN after it.

## spectrum_gate (`spectrum_type_gate.sh`)

Three one-step runs of the binary on copies of `backup/regression/wasp_full`
(the regression directory is not written to): the case as shipped, the same
case as `Spectrum type: Planck`, and that one with the two stellar lines
deleted. Asserts that the first two run, that the third is refused with a
message naming both keys, and that the setup report of every run states the
spectrum type, the source of the band below 13.6 eV, and the integrated grid
flux against the nominal `(10^LX + 10^LEUV)/(4 pi a^2)`.

## State

RED before the change of item 2a-PLANCK (there is no `Planck` type: the driver
does not build, and the gate's runs B and C stop with "unknown spectrum
type"), GREEN after it.

# Benchmark: WASP-52 b against Yan et al. (2022)

Purpose: compare EXHALE with the published 1-D escape model of WASP-52 b of
Yan et al. (2022), with the inputs matched one by one (plan rev4 step E). The
comparison quantities are the mass-loss rates of their Table 4 and the
profiles of their Figures 9, 12 and 23. This directory holds the setup, the
digitized published values and the six runs; the comparison itself is a later
phase.

Reference: Yan, D., Seon, K., Guo, J., Chen, G., and Li, L. 2022, ApJ, 936,
177, "Modeling the H-alpha and He 10830 Transmission Spectrum of WASP-52b",
`2022ApJ...936..177Y`, doi:10.3847/1538-4357/ac8793. Read in the publisher PDF
`references/Yan 2022.pdf` (`pdftotext -layout`; figure pages rendered with
`pdftoppm`). Their hydrodynamic code is that of Yang & Guo (2018), Chinese
Astron. Astrophys. 42, 81 (`references/Yang 2018.pdf`), also used by Yan et al.
(2021), ApJL 907, L47; where Yan et al. (2022) does not state an item, the
table says so and names what those papers state.

Runs use only the frozen binary `LHS1140b/models/.P2/EXHALE_f8a38866.x`
(plan rev4 A-D checkpoint, update log stage 3 section 136).

## Condition table

"READ" = quoted from the paper (section/table/figure given); "MEASURED" = run
or computed here; "not stated" = absent from Yan et al. (2022).

| Input | Published value (where READ) | EXHALE value | Matched? | Reason if not |
|---|---|---|---|---|
| Planet mass | 0.46 M_J (sec. 2.1) | `Planet mass [M_J]: 0.46` | yes | |
| Planet radius | 1.27 R_J (sec. 2.1) | `Planet radius [R_J]: 1.27` | yes | |
| Orbital distance | 0.0272 au (sec. 2.1); L = 44.861 R_p (Fig. 3 caption) | 0.0272 au (= 44.86 R_p with R_J = 7.1492e9 cm, MEASURED) | yes | |
| Stellar mass | 0.87 M_sun (sec. 2.1) | 0.87 | yes | |
| Stellar radius | 0.79 R_sun (sec. 2.1); R_S = 6.067 R_p (Fig. 3) | not used by the hydrodynamics (no `Stellar radius` key, so no Balmer continuum) | n/a | enters only Yan's Ly-alpha source geometry |
| SED source and shape | eps Eri, MUSCLES Treasury Survey, scaled to the orbit (sec. 2.1, Fig. 1); beta_m = 0.22 reference spectrum | `inputdata/sed/wasp52b_epseri_yan2022.txt` (MUSCLES v22 eps Eri x (d/a)^2); beta_m = 0.217 (MEASURED) | yes | Yan bins the XUV into 53 bins (sec. 2.1); EXHALE reads the 1 A table |
| F_XUV at the orbit | 34 168 erg cm^-2 s^-1 (sec. 2.1) | F(1-912 A) = 3.436e4, F(10-912 A) = 3.401e4 (MEASURED, trapezoid on the table); the run prints log F_XUV = 4.532, the 10-912 A band its photon grid covers | yes (0.6%) | |
| F0 and geometry | F0 = F_XUV/4, "uniform redistribution ... around the spherical planet" (sec. 2.1) | `2D approximate method: Rate/4`: dilution 0.250 of every stellar beam (printed), Mdot = 4 pi r^2 rho v (`Rate/4 + Mdot`) | yes | Yan does not state the Mdot formula. Fig. 9 settles it: 4 pi r^2 rho v from the digitized Fig. 9(b,c) curves (n_H1s + n_H+, He/H 0.087) is 2.62e11, 4.78e11, 8.48e11 g/s at 0.5, 1, 2 F0, flat to 2% over 2-8 R_p (MEASURED), against Table 4 2.73e11, 5.04e11, 9.08e11: the full-sphere rate, not 1/4 of it as in Yang & Guo (2018) Eq. 7 |
| F_XUV grid | 0.5, 1, 2 F0 (of 0.25-8 F0, Table 4) | 1 F0: the file above; 0.5 and 2 F0: `wasp52b_epseri_yan2022_xuv0p5.txt`, `..._xuv2.txt`, rows below 912 A times 0.5 or 2 (`inputdata/sed/build/build_wasp52b_epseri_xuv_scaled.py`); log F_XUV printed 4.230 at 0.5 F0 | yes | |
| FUV/NUV (912-2600 A) | 2.6e5 erg cm^-2 s^-1, one spectrum for all models (sec. 2.3, Fig. 1b), used for the He 2^3S photoionization in post-processing | 2.633e5 at the orbit (MEASURED), not scaled with F_XUV, diluted by 0.25 with every beam inside the run | partly | Whether Yan dilutes the 912-2600 A band by 4 is not stated. In Yan it enters only the post-processed He 2^3S balance; in EXHALE it photoionizes the coupled He 2^3S |
| Base radius | r = R_p (sec. 2.1: "temperature at the bottom boundary R_p") | `Base BC: pressure 1.0` at r = 1 R_p | yes | |
| Base pressure | 1 microbar (sec. 2.1) | 1.0e-6 bar, n0 = 5.554e12 cm^-3 H+He nuclei (printed) | yes | |
| Base temperature | 1304 K = T_eq (sec. 2.1) | `Equilibrium temperature [K]: 1304.0` | yes, as stated | CONFLICT inside the paper: the curves of Figs. 9(a), 12(a), 23(a) start near 2000-2300 K at r = 1-1.05 R_p (READ off the figures; digitized 1840-2100 K at 1.01-1.05 R_p) |
| Base composition | not stated in Yan 2022; the Yang & Guo (2018) code carries H2, H2+, H3+, HeH+ (their sec. 2). Yan's H(1s) is about 8e10 cm^-3 at 1.01 R_p (Fig. 9c, digitized) against 5.55e12 cm^-3 particles at 1 microbar and 1304 K, so either the base hydrogen is almost all molecular or the base density differs | atomic H and He, neutral | NO | molecules are off here (see "Molecules") |
| He/H | 92/8 and 98/2 (sec. 2.1, sec. 4, Table 4) | `He/H number ratio`: 0.0869565 (8/92) and 0.0204082 (2/98) | yes | |
| Metals | "solar abundance except [Fe/H] = 0.03" (sec. 2.1); the Yang & Guo (2018) code carries C, N, O, Si and their ions (photoionization and chemistry, sec. 2) | none (no `metals.inp`) | NO | Yan's cooling list (sec. 4.1: Ly-alpha, H3+, recombination, free-free) has no metal line cooling; EXHALE metals bring CHIANTI line cooling Yan lacks, so their photoabsorption cannot be matched without adding a cooling Yan does not have. Effect not measured (trace abundances) |
| Domain outer radius | 10 R_p (sec. 2.1; Fig. 3 R_A = 10 R_p) | `Outer radius [R_p]: 10.0` (printed "Upper boundary of the domain: 10.00 [R_p]") | yes | |
| Sonic point | not stated; from the digitized Fig. 9 profiles (adiabatic c_s, H ionization from Fig. 9c, He neutral) 1.99, 1.88, 1.81 R_p at 0.5, 1, 2 F0, Mach about 20 at 9.8 R_p (MEASURED) | outer face at 10 R_p, far past it; the outer-cell Mach number is read from each run's results | yes | both outer faces are supersonic, so no L1-domain variant is needed |
| Gravity | planet + "stellar tidal force" (sec. 2.1); Yang & Guo (2018) Eq. 4: gravity + stellar tide + revolution (centrifugal) | Roche mode (default): full tidal + centrifugal potential along the star-planet axis, extended past L1 to 10 R_p by `Outer radius` | yes | L1 at about 2.5 R_p (Hill radius 2.48 R_p, MEASURED arithmetic) lies inside both domains |
| Stellar heating, secondary ionization | Yang & Guo (2018) sec. 2.1: photoelectron heating with multiple (secondary) ionization for H, He, H2 | `Secondary_ionization` default True (Shull & van Steenberg 1985 budget) | yes (method not compared) | |
| Ly-alpha cooling | included (sec. 2.1; Fig. 10) | collisional-excitation Ly-alpha cooling, default on | yes | |
| Recombination cooling, free-free | included (sec. 4.1, Fig. 10) | default on | yes | |
| H3+ cooling | dominant below 1.1 R_p (sec. 4.1, Fig. 10) | off (no molecules) | NO | follows from "Molecules" |
| Molecules (H2, H2+, H3+, HeH+) | Yang & Guo (2018) code (sec. 2); H3+ cooling in Yan Fig. 10 | `Molecular chemistry` off | NO | the atomic run is the configuration the earlier memo compared with Table 4; the molecular option of the frozen binary is untested on this planet. A molecular variant is the advisor's decision |
| Charge exchange | "charge exchange" among the reactions (Yang & Guo 2018 sec. 2.1) | He+ + H charge exchange, default on | yes (rates not compared) | |
| Collisional ionization | included (Yang & Guo 2018 sec. 2.1) | default on | yes (rates not compared) | |
| Conduction, viscosity | not stated (nor in Yang & Guo 2018) | off (defaults) | unknown | |
| H(n=2) | post-processing only (sec. 2.2, LaRT Ly-alpha), not in the hydrodynamics | off: no `Stellar Teff` / `Stellar radius`, so no Balmer-continuum photoionization of H(n=2) in the run | yes for the hydrodynamics | H-alpha would be a post-processing comparison (EXHALE_transit.py or LaRT), not part of this setup |
| He 2^3S | post-processing only (sec. 2.3, Eq. 9, advection neglected; rates Table 2, Q31 = 5.0e-10 constant) | `Include He23S? True`: coupled in the run, with advection and the code's own rates (temperature-dependent Penning rate) | partly | kept on so that the He 2^3S profile of Fig. 23(d) can be compared; its feedback on the hydrodynamics is not in Yan's model (size not measured) |
| Diffusion | not stated | off | unknown | |
| Ly-alpha stellar line | double Gaussian, 49 km/s width at +-74 km/s, 37 300 erg cm^-2 s^-1 (EP Eri) (sec. 2.2), LaRT only | not used (H(n=2) off) | n/a | |

Numerical recipe (common brief): `Reconstruction scheme: PLM+WENO3`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Solver: Newton`, `IC mode: auto` (it chose
the transonic isothermal wind, cold sonic point 2.375 R_p, printed), Mixed
grid, 500 cells, HLLC.

## Published results to compare

- Table 4 (page 14), Mdot in 1e11 g/s, beta_m = 0.22 (READ):

  | F_XUV | H/He 92/8 | H/He 98/2 |
  |---|---|---|
  | 0.5 F0 | 2.73 | 2.59 |
  | 1 F0 | 5.04 | 4.92 |
  | 2 F0 | 9.08 | 9.01 |

  `published_yan2022_table4.txt` holds the whole table.
- Figure 9 (page 15), 92/8: T, v, n(H 1s), n(H+) against r/R_p for every
  F_XUV; digitized for 0.5, 1, 2 F0 in `published_yan2022_fig9_hhe92_8.txt`.
  Panels (d)-(f) are post-processed H(2) column densities (not digitized).
- Figure 12(c) (page 18), F0: n(H 1s), n(H+) for the four H/He; 98/2 in
  `published_yan2022_fig12c_hhe98_2.txt`.
- Figure 23 (page 29), 98/2: T, v, n(He), n(He+), n(He 2^3S) for every F_XUV;
  digitized for 0.5, 1, 2 F0 in `published_yan2022_fig23_hhe98_2.txt`.
- Figure 10 (heating and cooling rates, 92/8, F0) and Figure 22 (He, 98/2 to
  99.5/0.5 at F0) are not digitized.

Digitization: `digitize/` (color tracking on 300 dpi renders; reading
precision in each file header). Check against Table 4: the 4 pi r^2 rho v of
the digitized Fig. 9 curves lies 4-7% below Table 4 at all three fluxes (row
"F0 and geometry").

## Runs

| Directory | H/He | F_XUV | Spectrum file |
|---|---|---|---|
| `h92he8_f0p5` | 92/8 | 0.5 F0 | `wasp52b_epseri_yan2022_xuv0p5.txt` |
| `h92he8_f1` | 92/8 | 1 F0 | `wasp52b_epseri_yan2022.txt` |
| `h92he8_f2` | 92/8 | 2 F0 | `wasp52b_epseri_yan2022_xuv2.txt` |
| `h98he2_f0p5` | 98/2 | 0.5 F0 | `wasp52b_epseri_yan2022_xuv0p5.txt` |
| `h98he2_f1` | 98/2 | 1 F0 | `wasp52b_epseri_yan2022.txt` |
| `h98he2_f2` | 98/2 | 2 F0 | `wasp52b_epseri_yan2022_xuv2.txt` |

`run_table4_grid.sh` starts them with the frozen binary, at most three at once
with 8 threads each; start times and PIDs are in `run_table4_grid.log`.

## Earlier work

`WASP-52b/WASP-52b.tex` (2026-08-11, older code) used `Rate/2 + Mdot/2` with
the SED scaled by 0.5 for 1 F0: the same dilution (0.25 F_XUV) but half the
full-sphere Mdot. Its agreement with Table 4 (log Mdot 11.69 against 11.70 at
92/8, 1 F0) was therefore at a different Mdot convention from Yan's; the runs
here use the convention Fig. 9 establishes.

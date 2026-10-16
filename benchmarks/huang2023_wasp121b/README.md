# `benchmarks/huang2023_wasp121b`: WASP-121 b against Huang et al. (2023), cases A, B and D

## Purpose

Plan rev4 step E: benchmark EXHALE against a published escape model with the
inputs matched one by one. Reference: Huang, Koskinen, Lavvas & Fossati 2023,
*The Astrophysical Journal* 951, 123 (`2023ApJ...951..123H`,
doi:10.3847/1538-4357/accd5e), CETIMB coupled to a Ly-alpha Monte Carlo
transfer. Local PDF: `references/Huang_2023_ApJ_951_123.pdf` (read with
`pdftotext -layout`; figures rendered at 300 dpi). The text extract
`references/Huang_2023_ApJ_951_123_results.txt` agrees with the PDF where
compared (sections 2.7 and 3); every value below is READ from the PDF.

Cases run (all three published with Mdot and T, v profiles):

| case | what it is (paper) | profiles published | run directory |
|---|---|---|---|
| A | spherical potential, globally averaged (zenith 60 deg, flux / 2), section 3 | Fig. 9 (T, v), Fig. 10 (heating/cooling), Fig. 11 (proton rates), Fig. 12 (n of e, p, H(2p), metals) | `caseA/` |
| B | A plus the tidal (Roche) potential along the substellar line, vertical incidence, flux / 4, section 5 | Fig. 20 (T, v) | `caseB/` |
| D | the preferred model (Table 3 note e): B with M_p 5% lower, lower atmosphere +350 K, XUV x 0.74, Ly-alpha x 0.35, Fe x 4, section 6.1 | Fig. 25 (T, v), Figs. 26-29 | `caseD/` |

Case D was left out of the spring 2026 reproduction because the L1-truncated
domain was subsonic (update log stage 1 section 13). The current binary
extends a Roche domain past L1 (`Outer radius`, stage 1 section 20), and
Huang's case D sonic point (1.73 R_px) lies just past L1 (1.71 R_px), so the
outer face can be put past it; case D is therefore included. Cases C and E-K
are parameter variations without published profiles and are not run.

Binary: the frozen `LHS1140b/models/.P2/EXHALE_f8a38866.x` only.

## Spectrum: which normalization the paper's model used

Section 2.2 states both "F_XUV = 1.6 x 10^6 erg cm^-2 s^-1" (from Salz et al.
2019) and "the Lyman continuum flux of the spectrum is 2.69 x 10^5". No
solar-shaped SED satisfies both. Everything the paper computed with its
spectrum follows the second number:

* its stated band fluxes of the synthetic spectrum (1.4e4, 3.6e4, 6.3e4 over
  70-80, 80-91.2, 91.2-117 nm): a WHI 2008 solar shape scaled to
  F(1-912 A) = 1.6e6 gives 4.9x each (MEASURED), scaled to 2.69e5 it gives
  0.82, 0.81, 0.80 (MEASURED);
* Figure 1: F_lambda ~ 1e2-1e3 over 100-912 A at the orbit; its blackbody
  curve equals pi B_lambda(6459 K) (R_star/a)^2 at 3500 A (8.7e5, MEASURED),
  so the plot is undiluted at the orbit;
* Table 3 "M_XUV" = 0.16 M_p/Gyr (cases A, B) with their eq. (19),
  pi R_p^3 F / (G M_p): F = 2.69e5 gives 0.160, F = 1.6e6 gives 0.95
  (MEASURED arithmetic; case D with 0.74 x 2.69e5 gives 0.125 against 0.13).

So the runs use `inputdata/sed/wasp121b_composite_lyc_huang2023.txt` (cases A,
B) and `..._caseD.txt` (everything below 1700 A x 0.74), built by
`inputdata/sed/build/build_wasp121b_composite_lyc_huang.py`: WHI 2008 below
1700 A with F(1-912 A) = 2.69e5 and the Ly-alpha line 1.0e5 (D: 1.99e5 and
7.4e4), and the bolometrically scaled MUSCLES v24 WASP-17 (F6V) photosphere
above 1700 A. MEASURED band fluxes (1 A bins): F(10-912 A) = 2.691e5 (D
1.991e5); F(912-1170) = 5.08e4; F(1700-3647 A) = 6.62e8 erg cm^-2 s^-1;
SHA-256 `5b7ca2ec...` and `3d4d20eb...`. WASP-17 is chosen over the 6459 K
blackbody of the other candidate because Huang's Figure 1 synthetic and the
observed WASP-121 NUV spectrum both lie well below the blackbody in
1700-3000 A (line blanketing) and a measured F6V photosphere carries it; with
the 2.69e5 scaling the join at 1700 A closes to f(1700)/f(1699.5) = 0.79
(MEASURED). The old `wasp121b_composite_huang2023.txt` (1.6e6) is NOT used.

**Photosphere mismatch (unmatched):** LLmodels (Shulyak et al. 2004) for
WASP-121 is not available; WASP-17 (T_eff ~6550 K, its own gravity and
metallicity) stands in above 1700 A. This sets the Balmer-continuum
(3.4-13.6 eV) field that photoionizes H(n=2) and the field that ionizes the
neutral metals. Against Huang Fig. 1 by eye: 2000 A 8.4e4 vs ~5e4, 3000 A
5.1e5 vs ~5e5, 3500 A 6.2e5 vs ~7e5. Below 1700 A the shape is WHI 2008
(solar minimum), Huang's is "the solar SED" of Koskinen et al. 2013a (average
Sun); the three band ratios above (0.80-0.82) measure that shape difference.

## Condition table

R0 = EXHALE base radius ("Planet radius" key), the 1 microbar level.
"READ" = from the paper; "MEASURED" = computed or measured here.

| input | published (where READ) | EXHALE A / B / D | matched? | reason if not |
|---|---|---|---|---|
| planet mass | 1.1824 M_J (Tab. 1, Tab. 3 A, B); 1.1204 (Tab. 3 D) | 1.1824 / 1.1824 / 1.1204 | yes | |
| planet radius R_p (transit) | 1.766 R_J (Tab. 1) | not an input; used for unit conversion | n/a | EXHALE's radius key is the base level |
| base radius (1 microbar) | from log g @1 microbar = 2.84 / 2.70 / 2.62 (Tab. 3), g = G M_p / r^2 | 2.058 / 2.418 / 2.581 R_J (MEASURED from log g); = 1.165 / 1.369 / 1.461 R_p | yes | D agrees with the 1.46 R_p of section 7.1; base points of Figs. 9/20/25 at 1.160 R_p, 1.238 R_px (2.407 R_J), 1.311 R_px (2.568 R_J), within 0.5% |
| base pressure | 1 microbar (section 2.6; top axes 997-999 nbar) | `Base BC: pressure 1.0` (n0 = 1.83e12 / 2.69e12 / 2.37e12 cm^-3 printed) | yes | |
| base temperature | not tabulated; first point of Figs. 9 / 20 / 25 | 3960 / 2690 / 3050 K (READ from the digitized curves, +-60 K) | yes, to reading precision | Huang's base T comes from their photochemical lower atmosphere (Figs. 8, 21) |
| base composition | atomic, molecules depleted at 1 microbar (section 2.6); species densities from the photochemical model | atomic H/He + metals at solar abundance, the ionization equilibrium of the base cell | partly | Huang's base mixing ratios come from the photochemical model and are not tabulated; no molecules on either side |
| He/H | solar, Asplund et al. 2009 (section 2.4) | 0.0851 | yes | |
| metals | H, He, Mg, Fe, Si, O, C, N, S, Ca, Na, K, solar Asplund 2009 (section 2.4); D: Fe x 4 (Tab. 3) | the same ten metals, Asplund 2009 Table 1 values (`metals.inp`); D Fe 1.264e-4 | yes | |
| ionization stages | Mg, Fe, Si, Ca up to III; the rest neutral and II (section 2.4) | C, N, O, Mg, Si, Ca, Fe up to III; Na, K, S up to II | partly | EXHALE also carries C III, N III, O III; no key removes them |
| stellar mass, radius, a | 1.3521 M_sun, 1.4572 R_sun (Tab. 1); a = 0.02545 au (Tab. 3; Tab. 1 prints 0.02544) | 1.3521, 1.4572, 0.02545 | yes | 0.02545 au gives the Tab. 1 period 1.27484 d by Kepler (MEASURED); 0.02544 gives 1.27408 d |
| SED below 1700 A | solar SED scaled; Lyman continuum 2.69e5 at the orbit (section 2.2) | WHI 2008 shape, F(10-912 A) = 2.691e5 (`Log XUV flux` printed 5.430); D x 0.74 (printed 5.299) | integral yes, shape partly | the 1.6e6 of section 2.2 conflicts with the paper's own band fluxes, Fig. 1 and Tab. 3 (above); WHI 2008 vs Koskinen 2013a solar shape: bands 0.80-0.82 of Huang's |
| SED above 1700 A | LLmodels (section 2.2) | MUSCLES WASP-17 scaled to sigma T_eff^4 (R_star/a)^2, T_eff = 6459 K (the blackbody of Fig. 1) | no | LLmodels spectrum not available |
| stellar Ly-alpha | double-Voigt fit to HD 106516, line flux 1.0e5 at the orbit (section 2.2); D x 0.35 (Tab. 3) | `Stellar Lya flux` 1.0e5 / 1.0e5 / 3.5e4; flat-top line of half width 70 km/s | integral yes, profile no | EXHALE's beam is a top-hat closure; Huang's profile is not tabulated |
| geometry / dilution | A: zenith 60 deg and flux / 2 (section 2.2); B, D: vertical, flux / 4 (section 5.1) | A `Rate/2 + Mdot/2` (printed 0.500); B, D `Rate/4 + Mdot` (printed 0.250) | A partly, B, D yes | EXHALE has no slant path: A's column is vertical, Huang's is 1/cos 60 = 2 times longer. A's printed Mdot is halved by `Mdot/2`; compare 2x printed (= 4 pi r^2 rho v, Huang's definition, section 5.3) |
| Ly-alpha beam dilution | stellar Ly-alpha / 2 (section 2.7, stated for the method) | the same dilution as the XUV beam: 0.5 / 0.25 / 0.25 | A yes, B, D unclear | the paper does not say whether the Ly-alpha factor became 4 with the XUV in cases B, D |
| gravity | A: planet only; B, D: Roche potential eq. (18) along the substellar line | A `Domain mode: Spherical`; B, D Roche (planet + star + centrifugal, Kepler Omega) | yes | EXHALE potential (input_read.f90 / wae_ic_writer phiR) is eq. (18) with Omega^2 = G(M*+Mp)/a^3 |
| domain outer face | 31.7 R_J above the 1 microbar level, always past the sonic point (section 2.3): 16.40 / 14.11 / 13.28 R0 | 16.403 / 8.850 / 8.200 R0 | A yes, B, D no | at 14.11 and 12.0 R0 (B) and 13.28 R0 (D) this binary stops at the outermost cell on step 0-1 (suspected defect, below); B, D use the extent of Figs. 20 and 25 (11.0 R_px, 10.8 R_px), past the published sonic points (1.536, 1.313 R0) |
| grid | 580 stretched cells, dr0 = 10 km, stretch 1.014 | 580 cells, Mixed grid (50 uniform cells of 2e-4 R0 at the base) | partly | resolution choice, not physics |
| upper BC | outflow, constant-slope extrapolation (section 2.3) | EXHALE free outflow ghost | yes in kind | |
| photoelectron heating | constant efficiency 0.93 of the photoelectron energy (section 2.4.1) | `Photoelectron heating: excess` (1.00 of h nu - I), `Secondary_ionization: False` | no | no key for a fraction of h nu - I; heating is 7.5% higher than Huang's at the same absorption |
| H recombination | case B, Storey & Hummer fit eq. (3); cooling eq. (4) | `H_rec_escape: False` = case B everywhere; EXHALE's own case B coefficient | yes in kind | coefficient fits differ; not compared numerically |
| He recombination | Yelle 2004 rate (section 2.4.3) | `He_rec_coupling: False` (case B, cascade photons lost) | partly | EXHALE coefficients are not Yelle's |
| metal recombination | Shull & van Steenberg, UMIST, Verner & Ferland, Nahar fits (section 2.4.3) | EXHALE's metal rates (stage 1 sections 2-5) | not verified | not compared rate by rate |
| collisional ionization | Voronov 1997 for all species (section 2.4.1) | Voronov 1997 (H, He; metals as implemented) | yes for H, He | metals not verified |
| charge exchange | 65 reactions, Table 4 (section 2.4.2) | `cx_full 1` (Table 4 groups A, C, D), He/H pair on, `cx_O2p_H 0`, `cx_N2p_H 0` | yes | |
| Ly-alpha cooling | CHIANTI rates; Black (1981) is ~10x higher (section 2.5, Fig. 7) | CHIANTI H I collisional excitation (`Cool_coeff.f90 lambda_coex_HI`); its comment states Black/Cen is 0.90-0.97 of this CHIANTI sum | CONFLICT | Huang's CHIANTI Ly-alpha cooling is ~10x below Black 1981, EXHALE's CHIANTI value is within 10% of Black; the cause is not identified (paper conflict, below) |
| metal line cooling | Mg I, Mg II, Na I (two-level), Ca II, Fe II (CHIANTI), Fe I (Van Regemorter + NIST), free-free (section 2.5) | Mg, Ca, Na, Fe II, Fe I (Van Regemorter), C/N/O CHIANTI fits, bremsstrahlung | partly | EXHALE also cools by C, N, O lines (no off key; `cno_cool 0` swaps fits); Huang lists none for C, N, O |
| conduction | included (section 2.3) | `Conduction: True` | yes | conductivity forms not compared |
| viscosity | viscous drag and dissipation (section 2.3) | `Viscosity: True` | yes | |
| diffusion | multispecies diffusion (section 2.3) | off (`He_diffusion` default False) | no | EXHALE is single-fluid; its binary H/He element diffusion is not switched on here |
| secondary ionization | not included (constant heating efficiency) | `Secondary_ionization: False` | yes | |
| H(n=2) | 2s, 2p populations, Balmer-continuum ionization and heating, collisional de-excitation heating (section 2.7) | on (`Stellar Teff 6459`, `Stellar radius 1.4572`; grid floor 3.40 eV printed), `Deexc heat: True` | yes in kind | |
| Ly-alpha transfer | Monte Carlo, Hummer IIB partial redistribution, iterated with the hydro (section 2.7) | escape-probability closure (`Jlya escape-prob: True`), `Lya stellar boost 5` | no | EXHALE has no MC coupling; the boost was tuned to Huang Fig. 11 (`lya_rt.f90`), so an H(n=2) agreement is not independent |
| Ly-alpha lower boundary | not stated | EXHALE default (reflecting) | n/a | |
| He 2^3S | not modeled | `Include He23S? False` | yes | |
| molecules | none above 1 microbar (section 2.6) | none | yes | |

## Published results to compare

Mdot (Table 3, M_p/Gyr, 1 M_p/Gyr = 7.1e13 g/s): A 0.052 (3.72e12 g/s,
section 3), B 0.32, D 1.03. EXHALE: case A compare 2 x the printed Mdot
(`Mdot/2`), B and D the printed value (4 pi r^2 rho v). Huang's B and D
values are substellar-streamline rates (an upper limit, section 5.3).

Digitized profiles (headers give method and precision):

* `published_caseA_fig9_T_v.txt`: Fig. 9, T and v_z, 1.16-8.94 R_p;
* `published_caseB_fig20_T_v.txt`: Fig. 20, T 1.24-10.9 R_px, v_z to 8.05 R_px;
* `published_caseD_fig25_T_v.txt`: Fig. 25, T and v_z, 1.31-10.8 R_px;
* `published_caseA_fig12_densities.txt`: Fig. 12, n_e, H(2p), Mg I, Mg II,
  Fe I, Fe II, Ca II for r >= 1.17 R_p. Neutral H is not plotted, so the H
  ionization fraction is available only as the 50% markers;
* `published_values.txt`: Table 3 rows A, B, D, the figure markers (sonic
  point, L1, 50% ionization, tau = 1) and the domain statement.

Not digitized: Figs. 10, 11, 26-29 (rates, case D densities). The scripts,
300 dpi crops and trace overlays are in `digitization/`
(`digitize_fig9_fig20.py`, `digitize_fig12.py`, `write_published_profiles.py`).

## Runs (started 2026-10-11, lart4, 8 threads each)

| case | started | PID | outer face | IC (auto) |
|---|---|---|---|---|
| A | 00:39:15 | 58085 | 16.403 R0 | cold hydrostatic (no interior cold sonic point) |
| B | 00:41:50 | 62384 | 8.850 R0 | transonic isothermal (cold sonic point 1.379 R0) |
| D | 00:41:50 | 62385 | 8.200 R0 | transonic isothermal (cold sonic point 1.266 R0) |

The first B and D starts (00:39:15, outer face 14.110 and 13.282 R0) exited
with status 2 on step 0 and are kept in `caseB/failed_outer_face_31p7RJ/` and
`caseD/failed_outer_face_31p7RJ/`.

## Suspected code defect (reported, not fixed)

Roche domain extended far past L1: with the outer face at 12.0 or 14.11 R0
(case B) or 13.28 R0 (case D) the energy update of the outermost cell (582)
fails on step 0-1 ("non-finite temperature, step or source on entry",
ionization residual NaN; at 12.0 R0 "balance asks for a temperature below
the floor ... missing source 2.0e78 erg cm^-3 s^-1"). It does not depend on
viscosity, conduction or H(n=2) (each switched off in a 20-step probe,
MEASURED); at 8.85 R0 the same input runs. A source of 1e78 is not physical.
The star's surface is at 16.2 R0 (B), so the face is outside the planet's
reach but inside the stellar one. Probes: `~/.claude/jobs/1ff088bf/tmp/huang/probe/`.

Also: the setup report prints "Domain mode: Roche (full tidal potential,
truncated at the Hill/L1 radius)" for B and D although the upper boundary it
prints (8.85, 8.20 R_p) is past L1; the message text is stale for an
extended Roche domain.

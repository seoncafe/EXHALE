# Atomic Data in EXHALE versus MoCHII: Comparison and the 2026-07-17 Rate Update

Kwang-Il Seon — 2026-07-17

## Purpose and scope

This memo has two parts. Part 1 compares the atomic data that EXHALE and
MoCHII (a separate code, not part of this repository) use for the
microphysical processes the
two codes have in common: photoionization cross sections, recombination,
collisional ionization, charge exchange, cooling, and heating. The conclusion is
that, where the two codes overlap, they draw on the same lineage of published
data; the differences are in coverage (MoCHII spans more elements; EXHALE adds
the metastable He I 2^3S network and planet-atmosphere species) and in a few
modeling choices.

Part 2 documents the rate update applied on 2026-07-17. Before that update EXHALE
inherited its H/He recombination and collisional-ionization fits from ATES, and
several of those fits carried no source annotation. The update makes the EXHALE
default match the Badnell/Mao recombination and Voronov collisional ionization
that MoCHII uses, introduces the `legacy_hhe_rates` switch to reproduce the old
fits, replaces the free-free Gaunt factor with the van Hoof et al. (2014) thermal
average, and corrects a physical error in the free-free charge weighting.

## Part 1 — Process-by-process comparison

Both codes descend from the same published sources for each overlapping process.
EXHALE rows marked "[updated]" reflect the 2026-07-17 default (Part 2).

| Process | EXHALE | MoCHII |
|---|---|---|
| Photoionization cross sections | Verner et al. 1996 (VFKY96) with Verner & Yakovlev 1995 fallback; H I and He II from the hydrogenic analytic form. Adds He I 2^3S (VFKY96 wing bridge) and H2 (Yan, Sadeghpour & Dalgarno 1998). Tuned for planet atmospheres (Na, K). | Verner et al. 1996 with Verner & Yakovlev 1995 fallback; H I and He II hydrogenic. Wider element set (Ne, Cl, Ar). |
| Recombination (H, He) | [updated] Badnell (2023) radiative alpha_A minus Mao & Kaastra (2016) alpha_1, giving case B; He II includes Badnell dielectronic recombination. Identical coefficients to MoCHII. | Badnell (2023) radiative alpha_A and Mao & Kaastra (2016) alpha_1 for case B, with Badnell dielectronic recombination for He II. |
| Recombination (metals) | Badnell radiative/dielectronic recombination (Badnell 2006 and adf48 data), Ca I from Shull & Van Steenberg 1982, Fe from the Huang et al. 2023 fit
(neither the K-like nor the Mn-like/Cr-like sequences are covered by Badnell). | Badnell radiative/dielectronic recombination (`badnell_rr/dr.dat`) with a CHIANTI fallback. |
| Collisional ionization (H, He) | [updated] Voronov 1997. | Voronov 1997 (unified fits), with Dere 2007 as an option. |
| Collisional ionization (metals) | Voronov 1997. | Voronov 1997. |
| Charge exchange | Huang et al. 2023 (ApJ 951, 123) Table 4, based on Kingdon & Ferland 1996. 63 reactions (metal+H by default; He and metal-metal via the `cx_full` option). | Huang et al. 2023 Table 4 (Kingdon & Ferland 1996), metals only. |
| Recombination cooling | Hui & Gnedin 1997 case-B cooling fits. | Hui & Gnedin 1997 lineage, with a case-A/case-B choice. |
| Free-free cooling | [updated] Constant 1.426e-27; thermally-averaged Gaunt factor from van Hoof et al. 2014 (161-point table). | Constant 1.426e-27; Gaunt from Hummer 1988 / van Hoof et al. 2014. |
| Metal line cooling | CHIANTI v11.0.2 fit Lambda(T) = T^(-1/2) sum_i A_i exp(-T_i/T) (`cooling_data/fit_cno_formulas.py`), extended with density-dependent treatments (see below). | CHIANTI v11.0.2 fit of the same form (`chianti_cooling.py`), in the low-density limit; emission lines handled by a separate n-level solver. |
| Heating (photoelectric) | [updated] Shull & van Steenberg (1985) secondary ionization in the main loop: a photoelectron with E0 > 40 eV deposits only f_heat(x) of its excess as heat and drives H I / He I secondary ionizations (default on, `use_sec_ion`). Below 40 eV the photoelectron thermalizes fully. | Full photoelectron excess energy thermalized; no secondary ionization. |

### Photoionization cross sections
Both codes use Verner et al. (1996) fits with the Verner & Yakovlev (1995)
inner-shell data as a fallback, and both take H I and He II from the hydrogenic
analytic form. EXHALE adds two cross sections that MoCHII does not carry: the
metastable He I 2^3S level (bridged onto the VFKY96 wing) and molecular hydrogen
(Yan, Sadeghpour & Dalgarno 1998). MoCHII covers a broader element list
(including Ne, Cl, and Ar); EXHALE instead carries the alkali species (Na, K)
relevant to planet atmospheres.

### Recombination
Before 2026-07-17 the EXHALE H/He case-B coefficients were unattributed fits
inherited from ATES. A numerical check showed them to be exactly the Hui & Gnedin
(1997) case-B fits. MoCHII builds its default case-B coefficients from the Badnell
(2023) radiative alpha_A and the Mao & Kaastra (2016) alpha_1 subtraction, with
Badnell dielectronic recombination added for He II. The 2026-07-17 update (Part 2)
adopts the same construction as the EXHALE default, so the two codes now agree on
the H/He recombination coefficients. For metals, both codes already used Badnell
radiative and dielectronic recombination, with small differences in the fallback
data for a few ions.

### Collisional ionization
Before the update the EXHALE H I and He I collisional-ionization rates came from
Abel et al. (1997), which reproduces the Janev et al. (1987) fits, and He II from
a Hui & Gnedin (1997) fit; none of these carried a source annotation. MoCHII uses
the unified Voronov (1997) fits throughout, with Dere (2007) as an option. The
2026-07-17 update makes Voronov (1997) the EXHALE default for H and He as well;
the metal rates were already Voronov (1997) in both codes.

### Charge exchange
Both codes transcribe Table 4 of Huang et al. (2023, ApJ 951, 123), which rests
on Kingdon & Ferland (1996). EXHALE carries 63 of its 65 rows (23 in Group A, 2
in B, 6 in C, 32 in D; the two missing ones are in the metal-metal set). Two
switches select them. The 23 metal+H reactions of Group A are on by default and
`cx_full` adds the metal+He (C) and metal+metal (D) reactions. The He<->H pair
(B1, B2) is separate: it is not in the generic reaction list at all but is
applied by dedicated routines in every ionization system containing He, gated by
its own input key `He_H_charge_exchange` (default `True`) and independent of
`cx_full`. MoCHII carries the metal reactions only.

The printed Table 4 exchanges the reactant labels of its two oxygen rows: the
`exp(-227/T)` Boltzmann factor appears on `O+ + H0`, but `IP(O I) > IP(H I)`
makes `O0 + H+` the endothermic direction that must carry it. EXHALE assigns
the two rate coefficients to the physically correct rows; see
`HUANG2023_TABLE4_OXYGEN_ERRATUM.md` for the detailed-balance and Cloudy
cross-checks. Any code transcribing the table as printed inherits the error.

### Cooling
Recombination cooling in both codes follows the Hui & Gnedin (1997) case-B family
(MoCHII offers a case-A/case-B choice). Both use the same free-free constant,
1.426e-27; the Gaunt factor is discussed in Part 2. Metal line cooling in both
codes is a CHIANTI v11.0.2 fit of the form Lambda(T) = T^(-1/2) sum_i A_i
exp(-T_i/T) (EXHALE: `cooling_data/fit_cno_formulas.py`; MoCHII:
`chianti_cooling.py`). The MoCHII fit is the low-density limit, with the emission
lines themselves computed by a separate n-level solver.

### Heating
EXHALE now includes Shull & van Steenberg (1985) secondary ionization by fast
photoelectrons in the main loop (default on, `use_sec_ion`; see Part 2). MoCHII
thermalizes the full excess energy of the photoelectron and does not add
secondary ionization, which is appropriate in the nearly fully ionized H II
region it targets. The `wind_ae` initial-condition generator separately uses the
same SvS85 partition together with the Dere (2007) coefficients.

### Processes carried by only one code
EXHALE carries several networks that MoCHII does not: the metastable He I 2^3S
chemistry (Oklopcic & Hirata 2018; Bray 2000; Penning ionization from Taylor et
al. 2025), the H(n=2) 2s/2p populations (Christie, Arras & Li 2013; Draine 2011),
Lyman-alpha radiative transfer (Neufeld 1990), and the H2/H3+ molecular network
(Miller et al. 2013). Its cooling also adds density-dependent treatments absent
from the MoCHII low-density fits: a two-dimensional (T, n_e) statistical
equilibrium table for Fe II, the ground terms of C I, C II, N II and O I solved
in statistical equilibrium (which is where the [C II] 158 um and [O I] 63 um
saturation comes from), the
Mg II h&k, Ca II H&K, and Na I D two-level terms, and a Fe I table from NIST.

MoCHII in turn carries an n-level statistical-equilibrium emission-line solver
(CHIANTI v11), Storey & Hummer (1995) recombination lines, Porter et al. (2012,
2013) He I case-B emissivities, the nebular continuum (free-bound plus two-photon,
from Nussbaumer & Schmutz 1984 / Almog & Netzer 1989), and dust photoelectric
heating (Bakes & Tielens 1994).

One quantity crosses between them. MoCHII transports the He I 2^1S two-photon
continuum as sampled packets and so needs the full Drake, Victor & Dalgarno
(1969) shape; EXHALE never forms a spectrum, but its He II cascade coupling
(`he_rec_coupling`) needs two moments of that same shape above the H I edge --
how many ionizing photons a 2^1S decay yields and how much energy each of them
deposits. Both had been round numbers, 0.56 and 3.0 eV. Integrating the shape
gives 0.5564 and 2.512 eV: the photon count was right, but the deposited energy
was 20% high, being close to the 3.511 eV a uniform in-band distribution would
give rather than to the peaked one the pair actually has. MoCHII found the same
defect in its own sampler, which drew the in-band photons flat. EXHALE now
carries the two integrals as `f_2q_HeI` and `Ee_2q_HeI`.

## Part 2 — The 2026-07-17 rate update

The changes below are implemented in the radiation modules `Cool_coeff.f90` and
`util_ion_eq.f90`, in the temperature solver `T_equation.f90`, and in the input
parsing under `src/modules/files_IO/`. Items 1, 2, and 3 are gated by the new
`legacy_hhe_rates` switch; items 4 and 5 apply always, independent of the switch.

### 1. H/He recombination default
The default H/He case-B coefficients are now

    alpha_B = alpha_A(Badnell) - alpha_1(Mao),

where alpha_A(Badnell) is the Badnell radiative total (with the three-term
Badnell dielectronic recombination added for He II) and alpha_1(Mao) is the Mao &
Kaastra (2016, A&A 587, A84) ground-state coefficient. The implementation is in
`Cool_coeff.f90`: `rr_badnell`, `dr_HeII_badnell`, `rr_mao`, and the assembled
`alphaB_HII_new` / `alphaB_HeII_new` / `alphaB_HeIII_new`. The coefficients
match those in the MoCHII `src/recomb_mod.f90` (a file of the separate MoCHII
code, not of this repository).

### 2. H/He collisional ionization default
The default H/He collisional-ionization rates are now the Voronov (1997, ADNDT
65, 1) fits, implemented through `voronov_ci` (`Cool_coeff.f90`). The
adopted parameters are:

| Ion | dE [eV] | P | A | X | K |
|---|---|---|---|---|---|
| H I   | 13.6 | 0 | 2.91e-8 | 0.232 | 0.39 |
| He I  | 24.6 | 0 | 1.75e-8 | 0.180 | 0.35 |
| He II | 54.4 | 1 | 2.05e-9 | 0.265 | 0.25 |

### 3. The `legacy_hhe_rates` switch
A new input flag `legacy_hhe_rates` (default `False`) restores the previous fits:
Hui & Gnedin (1997) recombination and the Abel et al. (1997) / Hui & Gnedin (1997)
collisional ionization. It is parsed in `input_read.f90` (default set with the
other optional-block defaults, read in the keyword loop), declared in
`parameters.f90`, and echoed to `parse_dump.txt` by `write_setup_report.f90`.
It is registered as row
K14b in `docs/input_schema.md`. Setting it to `True` reverts items 1 and 2
exactly; the free-free changes (items 4 and 5) are not affected.

### 4. Free-free Gaunt factor
The free-free Gaunt factor is now the thermally-averaged non-relativistic
<g_ff>(gamma^2) of van Hoof et al. (2014, MNRAS 444, 420), with gamma^2 =
Z_ion^2 Ry/kT. The frequency average

    <g_ff>(gamma^2) = integral[ exp(-u) g_ff(gamma^2, u) du ] / integral[ exp(-u) du ]

is tabulated at 161 points in log10(gamma^2) = -6 to 10 in steps of 0.1 dex
(generated by `cooling_data/gauntff_thermal_avg.py`; the table is the `gff_avg`
parameter array in `Cool_coeff.f90`). This replaces the previous unattributed two-branch
logarithmic fit and applies always, regardless of `legacy_hhe_rates`. The
tabulated value at T = 1e4 K and Z = 1 is <g_ff> = 1.266.

### 5. Free-free charge-weighting correction
The previous code weighted He II (net charge +1) by Z^2 = 4 and evaluated the
metal Gaunt factor at the nuclear atomic number (Z = 8 for O, Z = 26 for Fe). This
is physically wrong: free-free emission scales with the ion *net* charge Z_ion,
not the nuclear number. The correction sets Z_ion = 1 for H II, He II, and
singly-ionized metals, and Z_ion = 2 for He III and doubly-ionized metals
(`util_ion_eq.f90` `eval_cool`, the `brem` accumulator; `T_equation.f90`, the
matching free-free term).
The most visible consequence is that the free-free cooling of a He II-dominated
layer drops by a factor of four. Like the Gaunt-table change, this correction
applies always.

### 6. Benchmark values
New default rates (alpha_B and collisional ionization ci; units cm^3 s^-1), with
legacy values at 1e4 K for comparison. The legacy alpha_B(H II) at 1e4 K is the
Hui & Gnedin (1997) fit, 2.592e-13.

| Coefficient | T = 1e4 K | T = 1e5 K | T = 1e6 K | legacy (1e4 K) |
|---|---|---|---|---|
| alpha_B(H II)   | 2.606e-13 | 3.070e-14 | 2.163e-15 | 2.592e-13 |
| alpha_B(He II)  | 2.807e-13 | 5.006e-13 | 1.015e-12 | 2.616e-13 |
| alpha_B(He III) | 1.536e-12 | 2.370e-13 | 2.244e-14 | --- |
| ci(H I)         | 7.448e-16 | 3.963e-9  | 3.103e-8  | 7.247e-16 |

The Fortran implementation reproduces an independent Python recomputation to all
printed digits, and the legacy branch reproduces the pre-update values exactly.
The new alpha_B(H II) at 1e4 K differs from the legacy Hui & Gnedin (1997) value
by about 1%; the He II and He III case-B coefficients differ more, mainly because
the new He II value carries the Badnell dielectronic contribution explicitly.

### 7. Secondary ionization in the main loop
A photoelectron ejected with kinetic energy E0 = e_v - E_th above about 40 eV
does not thermalize immediately: before it does, it collisionally ionizes further
H I and He I atoms. Shull & van Steenberg (1985) give asymptotic fits for the
partition of the photoelectron's energy as a function of the ionized fraction x
of the H+He nuclei,

  f_heat(x)    = 0.9971 (1 - (1 - x^0.2663)^1.3163),
  f_ion,HI(x)  = 0.3908 (1 - x^0.4092)^1.7592,
  f_ion,HeI(x) = 0.0554 (1 - x^0.4614)^1.6660.

As x -> 1 the heating fraction approaches unity and both ionization fractions go
to zero; as x -> 0 the heating fraction vanishes and the ionization fractions
reach their maxima (0.3908 and 0.0554). In neutral gas a large share of the
photoelectron energy therefore goes into secondary ionizations rather than heat.

A 40 eV threshold gates the partition: for E0 <= 40 eV the photoelectron
thermalizes fully (no secondary ionization), matching the wind_ae X-ray cutoff.
SvS85 is strictly an E0 >~ 100 eV asymptotic fit, so applying it down to 40 eV is
a deliberate approximation adopted here.

Only H I and He I receive secondary ionizations, the two channels SvS85 resolves.
The SvS85 Ly-alpha excitation channel is assumed to escape as line radiation and
is not put into any rate (a future refinement could couple it to the Ly-alpha
field of the excited-H model). Each absorbing species s (H I, He I, He II, the He I
triplet when active, H2 when present, and the photoionizable metal ions)
contributes photoelectrons of energy E0 = e_v - E_th,s that drive f_ion,HI
E0/E_th,HI secondary H I ionizations and f_ion,HeI E0/E_th,HeI secondary He I
ionizations. The He I triplet (2^3S, threshold 4.8 eV) previously contributed
opacity and its photoionization rate but no photoheating; its photoelectron
energy now enters the heating and absorbed-energy budgets and the secondary
source like every other absorber.

This differs from how wind_ae handles the same physics. wind_ae combines the
SvS85 f_heat / f_exc partition with the Dere (2007) coefficient tables to split
the secondary ionizations across species, but those Dere tables are tied to the
fixed wind_ae spectral grid and cannot be reused in the main loop, which runs on
a runtime SED and a variable energy grid. The main loop therefore applies the
species-resolved SvS85 fits (f_ion,HI, f_ion,HeI) directly.

The feature is gated by `use_sec_ion` (default True). With `use_sec_ion = False`
the photoelectron thermalizes fully and the heating and photoionization rates are
bit-identical to the legacy path. It is implemented in the `PH_heat_H` and
`PH_heat_HHe` routines of `util_ion_eq.f90`, wired through the callers in
`ionization_equilibrium.f90` and `post_process_adv.f90`; the flag is declared in
`parameters.f90`, parsed in `input_read.f90`, echoed by `write_setup_report.f90`,
and documented as key K14c in `input_schema.md`.

*[2026-08-15: since 2026-07 the coupling is staged rather than applied from step 0.
`use_sec_ion` still selects the physics, but a second runtime flag `sec_ion_active`
(`parameters.f90`, initially `.false.`) decides when it is actually applied: from a
cold IC the secondary-ionization base feedback amplifies the startup transient into
a runaway, so `EXHALE_main` flips it on only after the wind has first converged
without it, re-arms the marching stops, and lets the wind re-relax under the full
physics. The input override `Secondary_ionization: Immediate` (`sec_ion_immediate`,
default `.false.`) restores the pre-staging behavior for A/B tests.]*

### 8. Effect on existing results
Because the default rates change, results computed with the new default differ
from the pre-update golden and gate outputs. Setting `legacy_hhe_rates = True`
reproduces the old H/He rates and recovers the previous ionization balance, with
two exceptions that apply in all cases: the free-free charge-weighting correction
(item 5) and the van Hoof et al. (2014) Gaunt table (item 4). Those two changes
are physical corrections and are intentionally not tied to the switch. Secondary
ionization (item 7) is a further default change on top of the rate update;
`use_sec_ion = False` restores the full-thermalization heating and photoionization
rates bit-identically.

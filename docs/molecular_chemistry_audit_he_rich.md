# Molecular chemistry in the He-dominated limit: an audit

Phase F of `lhs1140b_lower_atmosphere_plan_new.md` asks for an audit of
R16-R20, R23, the shared H-He charge-exchange path, He 2<sup>3</sup>S + H<sub>2</sub> Penning
ionization, and the electron and H-nucleus closures "in the He-dominated
limit (baseline row 10)".  This is the record.

## Verdict

Three code errors, one of them changing results; four published-rate
caveats that are recorded rather than patched; and the closures hold.

* **The third body M of the three-body reactions was the wrong density.**
  `ionization_equilibrium` passed `n_in_dim` = &rho;/m<sub>H</sub> - a *mass*
  density in m<sub>H</sub> units - where R12/R13/R15 want the number density of
  third bodies, and where the interface it fills (`ion_cell_state%ntot`)
  says "total particle density".  It over-counts M by the mean particle
  mass: 2.3x at an H<sub>2</sub>-rich base and up to 4x in a He-dominated one, i.e.
  the error grows exactly along the axis this audit examines.  Fixed
  (`calc_ntot`).  Goldens move.
* **The retry seed's gas pressure had the same error.** The
  chemical-equilibrium restart of a failed molecular cell was seeded at
  p = (&rho;/m<sub>H</sub>) k<sub>B</sub> T instead of the EOS's own
  p = (n<sub>tot</sub> + n<sub>e</sub>) k<sub>B</sub> T. Fixed.
* **The admissibility test and the clamp ignored the He nucleus bound in
  HeH<sup>+</sup>.** The molecular systems close neutral helium as
  n<sub>He</sub>(1 - x<sub>2</sub> - x<sub>3</sub>) - n<sub>HeH+</sub>, so a
  root with x<sub>2</sub> + x<sub>3</sub> = 1 and x<sub>7</sub> > 0 gives a
  negative neutral He and was still accepted.  Fixed in
  `ionization_fractions_physical`, `element_budget_violation` and
  `clamp_fractions_to_element_budget`.  Unreachable in every case measured
  here (HeH<sup>+</sup>/He stays below 1e-10 over He/H = 0.079 to 1000), so it moves
  no result; it is a guard that was not guarding.

**Closures.** With molecules on and He/H from 0.079 to 1000, charge
neutrality &Sigma;Zn = n<sub>e</sub> closes to 4e-16, the H-nucleus/mass
closure to 4e-16, and `src/utils/element_budget.py` closes both elements to
6e-14 at He/H = 1000.  Nothing degrades with the He fraction.  What does
degrade is the *solve*: the number of cells whose molecular root leaves the
physical simplex rises steeply with He/H, and two of the four scans still
end on a NaN abort (below).

**Rates.** All 23 Koskinen et al. (2022) Table-1 coefficients are
transcribed correctly, verbatim, including R16-R20 and R23; so are the
Huang et al. (2023) Table-4 B1/B2 charge-exchange pair, the Taylor et al.
(2025) He(2<sup>3</sup>S)+H fit and the Garc&iacute;a Mu&ntilde;oz (2025) He(2<sup>3</sup>S)+H<sub>2</sub> fit.  No
coefficient was wrong.  Four physics caveats specific to a He-dominated gas
are recorded in section 4 and, in shorter form, at the code sites.

## 1. Sources

| Quantity | Source, verified against | Where |
|---|---|---|
| R1-R23 | Koskinen et al. (2022), ApJ 929:52, Table 1 (p. 19) | `references/Koskinen_2022_ApJ_929_52.pdf` |
| H &harr; He charge exchange (B1, B2) | Huang et al. (2023), ApJ 951:123, Table 4, group B | `references/Huang_2023_ApJ_951_123.pdf`, transcription `docs/charge_exchange_table4.md` |
| He(2<sup>3</sup>S) + H Penning | Taylor et al. (2025), ApJ 989:68, Table 2 and Sec. 2.4 | `references/Taylor_2025_ApJ_989_68.pdf` |
| He(2<sup>3</sup>S) + H<sub>2</sub> Penning, HeH<sup>+</sup> cross-check | Garc&iacute;a Mu&ntilde;oz (2025), A&A 698, A199, Tables A.5-A.7 | `references/GarciaMunoz_2025_A+A_698_A199.pdf` |

Text was extracted with `pdftotext -layout` from the publisher PDFs and
compared line by line with the code.

## 2. Reaction table

Code expressions are from `src/modules/lower_atmosphere/mol_rates.f90`;
T is the single temperature EXHALE carries (the table's T and T<sub>e</sub> coincide).
Rows R16-R23 and the three three-body rates are the audit's subject; R1-R15
are listed for completeness of the balance check.

| # | Reaction | Code | Table 1 | Verdict | He-dominated note |
|---|---|---|---|---|---|
| R5 | H<sub>2</sub><sup>+</sup> + e &rarr; H + H | 2.3e-8 (300/T)<sup>0.4</sup> | same | ok | |
| R6 | H<sub>3</sub><sup>+</sup> + e &rarr; H<sub>2</sub> + H | 2.16e-8 (300/T)<sup>0.65</sup> | same | ok | |
| R7 | H<sub>3</sub><sup>+</sup> + e &rarr; 3H | 5.04e-8 (300/T)<sup>0.65</sup> | same | ok | |
| R8 | H<sub>2</sub><sup>+</sup> + H<sub>2</sub> &rarr; H<sub>3</sub><sup>+</sup> + H | 2.0e-9 | same | ok | |
| R9 | H<sub>2</sub><sup>+</sup> + H &rarr; H<sup>+</sup> + H<sub>2</sub> | 6.4e-10 | same | ok | |
| R10 | H<sup>+</sup> + H<sub>2</sub>(v&ge;4) &rarr; H<sub>2</sub><sup>+</sup> + H | 1e-9 exp(-21900/T) | same | ok | Table 1 misprints the products as "H<sup>+</sup> + H<sub>2</sub>"; the code has the Yelle (2004) reaction |
| R11 | H<sub>3</sub><sup>+</sup> + H &rarr; H<sub>2</sub><sup>+</sup> + H<sub>2</sub> | 2.1e-9 exp(-20000/T) | same | ok | |
| R12 | H<sub>2</sub> + M &rarr; H + H + M | 1.5e-9 exp(-48350/T), &times; n<sub>M</sub> | **superseded 2026-09-01**: the code now builds this rate by detailed balance of R15, `k3b_H_H_to_H2(T)/keq_H_H_to_H2(T)` | ok | M identity: sec. 4.1. The Table-1 fit is valid only over 2500-8000 K, above the molecular layer; the argument, the sources and the verification numbers are at the R12 comment in `mol_rates.f90` |
| R13 | H<sup>+</sup> + H<sub>2</sub> + M &rarr; H<sub>3</sub><sup>+</sup> + M | 3.2e-29 n<sub>M</sub> | same | ok | M identity: sec. 4.1 |
| R14 | H<sub>2</sub> + e &rarr; H + H + e | 1.33e-6 (300/T)<sup>0.91</sup> exp(-55800/T) | same | ok | |
| R15 | H + H + M &rarr; H<sub>2</sub> + M | 8e-33 (300/T)<sup>0.6</sup> n<sub>M</sub> | **changed 2026-09-01**: 2.8e-31 T<sup>-0.6</sup> n<sub>M</sub>, the Cohen &amp; Westberg (1983) recommendation (50-5000 K), 14% above the Ham et al. value the table carries | ok | a neutral three-body recombination has no electron in it; the heavy-particle T is used, and the two readings coincide in EXHALE. Ham et al. measured 77-300 K only, so the Table-1 value is extrapolated across the molecular layer; see the R15 comment in `mol_rates.f90` |
| R16 | HeH<sup>+</sup> + e &rarr; He + H | 1e-8 (300/T)<sup>0.6</sup> | same | transcription ok | 3.4x below Garc&iacute;a Mu&ntilde;oz Table A.6 at 500 K and 8.6x below at 10<sup>4</sup> K (both fall with T, this one faster): sec. 4.2 |
| R17 | He<sup>+</sup> + H<sub>2</sub> &rarr; H<sup>+</sup> + H + He | 1e-9 exp(-5700/T) | same | transcription ok | up to 1.9e4 above Garc&iacute;a Mu&ntilde;oz Table A.7: sec. 4.2. Dominant H<sub>2</sub> sink and H<sup>+</sup> source once He<sup>+</sup> exists |
| R18 | HeH<sup>+</sup> + H<sub>2</sub> &rarr; H<sub>3</sub><sup>+</sup> + He | 1.5e-9 | same | ok | 1.2x above Table A.6 (Orient 1977) |
| R19 | HeH<sup>+</sup> + H &rarr; H<sub>2</sub><sup>+</sup> + He | 9.1e-10 | same | ok | 1.4-2.6x below Table A.6 (De Fazio 2014) |
| R20 | He<sup>+</sup> + H<sub>2</sub> &rarr; HeH<sup>+</sup> + H | 4.2e-13 | same | transcription ok | 14x above the Schauer et al. (1989) total Garc&iacute;a Mu&ntilde;oz adopts; the only HeH<sup>+</sup> source in the network: sec. 4.2 |
| R21 | H + He<sup>+</sup> &rarr; H<sup>+</sup> + He | 1.2e-15 (T/300)<sup>0.25</sup> | same | transcribed, **not used** | superseded by Huang B2, 4% away: sec. 4.3 |
| R22 | H<sup>+</sup> + He &rarr; H + He<sup>+</sup> | 1.75e-11 (300/T)<sup>0.75</sup> exp(-128000/T) | same | transcribed, **not used** | superseded by Huang B1: sec. 4.3 |
| R23 | H<sub>2</sub> + He<sup>+</sup> &rarr; H<sub>2</sub><sup>+</sup> + He | 7.2e-15 | same | ok | |

Balance check of `mol_heh_rows` (`System_HeH_mol.f90`): every reaction above
appears in each row it belongs to, with the right stoichiometric weight and
sign, and the H-nucleus and He-nucleus closures are the residual's own
(H I = 1 - x<sub>1</sub> - x<sub>4</sub> - ... - x<sub>7</sub> per H nucleus;
free neutral He = n<sub>He</sub>(1 - x<sub>2</sub> - x<sub>3</sub>) - n<sub>HeH+</sub>).
No missing or double-counted term was found in the eight rows.

### Shared H-He charge-exchange path

`System_HeH_mol` does not use R21/R22.  It calls `he_h_cx_fvec`, the
dedicated Huang Table-4 B1/B2 pair that every EXHALE system carrying He
shares, so the molecular-free limit reproduces the atomic systems exactly.
The two pairs agree: B2 = 1.25e-15 (T/300)<sup>0.25</sup> against R21's 1.2e-15 (4%),
B1 = 1.75e-11 (T/300)<sup>-0.75</sup> exp(-12.75/T<sub>4</sub>) against R22's identical
prefactor and exponent with the 127,500 K barrier rounded to 128,000 K.
The module header of `System_HeH_mol` previously said the H-He exchange was
"excluded"; it is not, and the comment is corrected.

### He(2<sup>3</sup>S) Penning channels

*Superseded 2026-08-28; see section 4.4 and `docs/Update_EXHALE_stage1.md` section 86.*
Both rate coefficients now come from Garc&iacute;a Mu&ntilde;oz (2025) and both are
branched:

| Channel | Code (total ionization) | Penning branch | Source |
|---|---|---|---|
| He(2<sup>3</sup>S) + H | 1e-9 exp(-86.4804/T - 0.286766 lnT + 0.0868445 (lnT)<sup>2</sup> - 0.00573001 (lnT)<sup>3</sup>) | 0.9 x total | Garc&iacute;a Mu&ntilde;oz (2025) Fig. 5, from the Movre &amp; Meyer (1997) cross sections; network rows 198/199 |
| He(2<sup>3</sup>S) + H<sub>2</sub> | 5.408222e-12 T<sup>0.675388</sup> exp(-696.275/T) | 0.9 x total | Garc&iacute;a Mu&ntilde;oz (2025) Table A.5, from the Cohen &amp; Lane (1977) cross sections; network rows 202/203 |

Both tabulated values are *totals*, and Garc&iacute;a Mu&ntilde;oz partitions them
"an average 0.9:0.1" between Penning and associative ionization.  The
published network file (`references/garcia_munoz_2025_network/`) applies
that split through the amplitudes of two separate rows per channel, which
settles what the paper text leaves implicit.  EXHALE now follows it: the
*total* removes the metastable, and `f_penning_HeI23S` = 0.9 scales the
terms that create a lasting proton or H<sub>2</sub><sup>+</sup> and that deposit the Penning
exothermicity.  The molecular system, which carries HeH<sup>+</sup> explicitly, takes
the remaining 0.1 into its HeH<sup>+</sup> row, so this is no longer an omitted
HeH<sup>+</sup> source there.  The atomic systems still drop it, on the ground that
HeH<sup>+</sup> dissociatively recombines back to He + H in an H<sub>2</sub>-poor gas.

Previously the code used Taylor et al. (2025) Table 2 for the H channel and
a fit made here to the four Table A.5 points for the H<sub>2</sub> channel, and
charged 100% of both totals to the Penning branch.

The Penning heating term
(n<sub>HeTR</sub> n<sub>H2</sub> k [E(2<sup>3</sup>S) - IP(H<sub>2</sub>)], with
E(2<sup>3</sup>S) = e_th_HeI - e_th_HeTR = 19.79 eV and IP(H<sub>2</sub>) = 15.43 eV,
so 4.36 eV to the electron) was checked and is right.

## 3. Closures measured

Configuration: `backup/regression/mol_diffusion` (hot Uranus, molecular
chemistry + He 2<sup>3</sup>S + binary H/He diffusion + a `base.inp` handoff), copied
to scratch with **only** `He/H number ratio` in `input.inp` and `HeH_base`
in `base.inp` changed together.  `OMP_NUM_THREADS=1`, `EXHALE_MAXSTEPS=12000`,
marching only, no Newton finish.  The He/H = 0.0793 run reproduces the
golden bitwise before the fix, which is what makes the scan a controlled
one.

These are relaxation snapshots, not converged winds (the flux spread `du`
ends between 2.6 and 717), so the `Mdot` column is a diagnostic of where the
run went, **not** a mass-loss rate.

Before the fix:

| He/H | steps run | NaN abort | H<sub>2</sub> fraction, cell 1 | front, f<sub>H2</sub> = 0.5 | front, f<sub>H2</sub> = 0.01 | max HeH<sup>+</sup>/He | max &vert;&Sigma;Zn - n<sub>e</sub>&vert;/n<sub>e</sub> | mass closure | clamped cells |
|---|---|---|---|---|---|---|---|---|---|
| 0.0793 | 12000 | no | 0.9997 | 1.160 | 1.193 | 8.6e-12 | 2.9e-16 | 4.5e-16 | 0 |
| 1 | 3402 | yes | 0.8768 | 1.224 | 1.277 | 9.3e-11 | 4.0e-16 | 4.2e-16 | 371 |
| 10 | 4855 | yes | 0.0000 | 1.263 | 1.767 | 4.2e-15 | 3.5e-16 | 3.9e-16 | 181 |
| 1000 | 12000 | no | 0.8776 | 1.003 | 1.020 | 2.4e-13 | 2.8e-16 | 4.0e-16 | 37 |

After the fix:

| He/H | steps run | NaN abort | H<sub>2</sub> fraction, cell 1 | front, f<sub>H2</sub> = 0.5 | front, f<sub>H2</sub> = 0.01 | max HeH<sup>+</sup>/He | max &vert;&Sigma;Zn - n<sub>e</sub>&vert;/n<sub>e</sub> | mass closure | clamped cells |
|---|---|---|---|---|---|---|---|---|---|
| 0.0793 | 12000 | no | 0.9998 | 1.160 | 1.193 | 8.7e-12 | 3.0e-16 | 4.1e-16 | 0 |
| 1 | 5538 | yes | 0.0000 | 1.260 | 1.273 | 7.4e-11 | 2.8e-16 | 4.1e-16 | 370 |
| 10 | 12000 | no | 0.0000 | 1.036 | 1.060 | 1.2e-15 | 4.2e-16 | 4.1e-16 | 169 |
| 1000 | 12000 | no | 0.0000 | 1.004 | 1.021 | 4.9e-13 | 3.8e-16 | 4.0e-16 | 4 |

`src/utils/element_budget.py` passes on all eight runs: the H and He
reservoirs close to 8.6e-16 at He/H = 10 and to 6.2e-14 at He/H = 1000
(base cell, which is where the check is made when `He_diffusion` is on).
No negative density appears in any output column of any run.

Reading:

* **The closures do hold to He/H = 1000**, molecules included, before and
  after the fix.  Charge neutrality is the residual's own definition
  (n<sub>e</sub> = n<sub>H+</sub> + n<sub>H2+</sub> + n<sub>H3+</sub> + n<sub>HeH+</sub> + n<sub>He+</sub> + 2n<sub>He++</sub>) and it survives the
  hydro round trip to machine precision.  Element conservation is not
  where the He-dominated limit hurts.
* **The molecular solve is where it hurts.**  Cells whose every molecular
  root leaves the physical simplex, and which are therefore clamped onto
  the element budget, go from 0 (He/H = 0.079) to hundreds (He/H = 1 and
  10).  The suspected mechanism is row scaling: rows 2 and 3 balance
  helium and carry terms of order n<sub>He</sub><sup>2</sup>, rows 4-7 balance the molecules
  and carry terms of order n<sub>H</sub><sup>2</sup>, so at He/H = 10<sup>3</sup> the two blocks of the
  same residual vector differ by 10<sup>6</sup>.  `hybrd1` is given no row scaling.
  This was not pursued further here.  **Pursued 2026-08-29 (section 7): the
  scaling is now applied, and what it removes is the solver failure, not the
  simplex failure.**
* **The M fix visibly improves the He-rich end**: He/H = 10 no longer
  aborts on a NaN and runs the full 12000 steps, and the clamped-cell count
  at He/H = 1000 falls from 37 to 4.  He/H = 1 still aborts, later (5538
  steps instead of 3402).  The abort is a controlled one - the marching
  loop's NaN detector stops the run and writes the state - and it happens
  in the hydro (p and T go non-finite at one cell while that cell's
  species densities stay finite), not in the network.
* At He/H = 1000 the base is molecular but **isolated cells fall into the
  atomic basin**: f<sub>H2</sub> reads 0.88, 0, 0, 0.885, 0.872, ... over the first
  cells.  This is the documented bistability of the molecular `hybrd1`
  solve showing up cell by cell, not a budget violation - those cells still
  close every element.

## 4. Caveats recorded, not patched

### 4.1 What is M?

Koskinen et al. write R12/R13/R15 as a two-body coefficient times "n" and
never say which particles n counts.  Their models are H<sub>2</sub>/H-dominated, where
M is H<sub>2</sub> and H.  EXHALE now passes the total gas-particle density
(electrons excluded), so every heavy particle is an equally efficient third
body.  A monatomic third body has no internal modes to absorb the released
energy, so He is in general *less* efficient than H<sub>2</sub>, and in a
He-dominated envelope R12/R13/R15 are upper bounds.  No He-specific
efficiency factor is applied, because none of the three quoted sources
(Baulch et al. 1992; Miller et al. 1968; Ham et al. 1970) supplies one that
we have verified; a factor would have to come from the primary literature
and is a change of physics, not of transcription.

Scale of the term: with the correct M, R15 and R13 are the reactions that
build H<sub>2</sub> and H<sub>3</sub><sup>+</sup> at the base, and the audit measured what the *wrong* M
did - see section 3, where the He-rich bases change basin.  An efficiency
factor of order 0.5 for He would move R12/R13/R15 by that factor in the
He-dominated base and nothing elsewhere.

### 4.2 HeH<sup>+</sup> in a He-dominated gas

Table 1 forms HeH<sup>+</sup> through R20 only, He<sup>+</sup> + H<sub>2</sub>.  That route needs He<sup>+</sup>,
which the shielded molecular base does not have.  The route that does not
need He<sup>+</sup> - H<sub>2</sub><sup>+</sup> + He &rarr; HeH<sup>+</sup> + H - is absent from Table 1, and
Garc&iacute;a Mu&ntilde;oz (2025) Table A.6 gives it (Black 1978) as 4.4e-16, 1.0e-11,
7.8e-11, 1.5e-10 cm<sup>3</sup> s<sup>-1</sup> at 500, 2000, 5000, 10<sup>4</sup> K.  It is the HeH<sup>+</sup>
source that grows with the He fraction, and it is the one the network does
not have.  The ~10% associative branches of both Penning channels (section 2)
end in HeH<sup>+</sup> too and are likewise not resolved.

Against the same Garc&iacute;a Mu&ntilde;oz tables, the HeH<sup>+</sup> destruction rates the
network does carry differ by factors of order unity to ten (R16 3.4-8.6x
low, R18 1.2x high, R19 1.4-2.6x low), and He<sup>+</sup> + H<sub>2</sub> differs far more: R17
is up to 1.9e4 above and R20 14x above the single T-independent 3e-14 that
Garc&iacute;a Mu&ntilde;oz adopts from the Schauer et al. (1989) experiment for the whole
He<sup>+</sup> + H<sub>2</sub> process.

Both compilations are published and internally consistent; picking between
them is a modelling decision, not a bug fix, and it would replace a
Koskinen Table-1 entry with a rate from another network.  The measured
HeH<sup>+</sup> abundance in this configuration is at most 1e-10 of the helium at any
He/H tested, so nothing observable rests on it today.  **If HeH<sup>+</sup> or the
He<sup>+</sup> + H<sub>2</sub> branching is ever quoted as a result, this paragraph is the
thing to resolve first.**

One HeH<sup>+</sup> formation path the network does not carry, recorded so that its
absence is a choice rather than an oversight: **radiative association**,
He<sup>+</sup> + H &rarr; HeH<sup>+</sup> + &nu;, which is reaction 7 of Courtney, Forrey, McArdle,
Stancil &amp; Babb (2021), ApJ 919, 70. Read from their Figure 3, its rate
coefficient is about 2.5e-16 cm<sup>3</sup> s<sup>-1</sup> near 10<sup>3</sup> K - roughly 7x below the
radiative charge-transfer channel in the same figure, and emphatically not the
vanishing number it is sometimes assumed to be. It is nonetheless negligible
*here*, for a specific reason: against R20 (He<sup>+</sup> + H<sub>2</sub> &rarr; HeH<sup>+</sup> + H,
4.2e-13 cm<sup>3</sup> s<sup>-1</sup>) it contributes only where n(H I)/n(H<sub>2</sub>) exceeds about 1700,
i.e. only where H<sub>2</sub> has already gone - and there HeH<sup>+</sup> itself is negligible.
Not adopted.

### 4.3 The two H &harr; He charge-exchange rates are different channels (revised 2026-08-31)

**This section previously reported the pair as violating detailed balance and
weighed "correcting" it by overriding one published rate with the
detailed-balance image of the other. That framing was wrong, and so was its
title.** The two rates describe two different physical processes, between
which no detailed-balance relation holds. The measured ratios below are
unchanged and correct; only their interpretation is.

Tracing each direction to its primary source:

* **B1**, He + H<sup>+</sup> &rarr; He<sup>+</sup> + H, is **non-radiative** collisional charge
  transfer, from Kimura, Lane, Dalgarno &amp; Dixson (1993), ApJ 405, 801,
  whose Table 3 tabulates it from 6000 K to 10<sup>5</sup> K only.
* **B2**, He<sup>+</sup> + H &rarr; He + H<sup>+</sup> + &nu;, is **radiative** charge transfer,
  printed with the photon in the exit channel as row (19) of Table 1 of
  Stancil, Lepp &amp; Dalgarno (1998), ApJ 509, 1, from Zygelman et al. (1989)
  and already carrying a 0.25 approach-probability factor.

The photon is the whole point: a photon-emitting channel has no collisional
reverse, so the equilibrium relation fixed by the statistical weights and the
10.98 eV (127,500 K) endothermicity,

    k(H+ + He) / k(He+ + H) = (g_He+ g_H)/(g_He g_H+) exp(-127500/T) = 4 exp(-127500/T),

was never binding on this pair. Both coded fits share the same exponential, so
the departure from it is exactly `1.05e6/T`: measured, **939x at 1140 K, 211x
at 5000 K and 105x at 10<sup>4</sup> K** (3500x at 300 K, and the Koskinen R22/R21
pair, from the same underlying fits, gives 104x, 198x and 689x). Those factors
measure the separation of two channels, not an error in either fit. The
label in Huang et al. (2023) Table 4, which both this audit and the code
followed, had dropped the `+ &nu;` and with it the reason.

The endothermic direction remains negligible below 10<sup>4</sup> K either way - at
8000 K, k(H<sup>+</sup> + He) = 1.8e-19 against k(He<sup>+</sup> + H) = 2.8e-15 - and B1 falls to
1.7e-60 at 1140 K, so in the shielded molecular base the pair acts in one
direction only.

**Independent corroboration that this pairing is normal practice.** Ziegler, U.
(2018), A&amp;A 620, A81, a 121-species / 426-reaction chemistry and cooling module
for NIRVANA, lists in its Table A.1

    7   He+ + H -> He + H+     ref 5
    8   He + H+ -> He+ + H     ref 6

with its reference list giving `5: Zygelman et al. (1989), 6: Kimura et al.
(1993)` - exactly the two sources behind EXHALE's B2 and B1. An independent
network built five years before ours therefore pairs the same two calculations
for the same two directions. That is outside confirmation that carrying two
independently computed rates here reflects the two channels, and is not a
defect anyone has been overlooking.

Ziegler carries He<sup>+</sup> + H as a **single** reaction row from the Zygelman
calculation rather than splitting it into radiative and non-radiative parts,
which is the same shape as the treatment adopted in section 4.3.1 below (one
summed removal rate). The comparison stops there: his numerical coefficients
are not printed in the paper, so **whether he actually included both channels
cannot be decided from it**, and nothing further is claimed. What can be said
is that EXHALE is now the more explicit of the two about what went into the
sum.

His section 2.1 also supplies a published precedent for bounding a Kingdon &amp;
Ferland fit outside its stated range - see section 4.3.1.

What the old framing did obscure is a real omission, now closed in section
4.3.1: the **non-radiative** He<sup>+</sup> + H channel was absent from the network
altogether.

#### 4.3.1 The missing non-radiative He<sup>+</sup> + H channel (adopted 2026-08-31)

Zygelman, Dalgarno, Kimura &amp; Lane (1989), Phys. Rev. A 40, 2340 - the common
source behind both coded rates, read here from the published pages - computes
the two channels separately and keeps them apart. Its Table I is captioned
"Rate coefficients for direct radiative charge transfer and for total
radiative decay in units of 10<sup>-15</sup> cm<sup>3</sup> sec<sup>-1</sup>" and tabulates only the
**radiative** processes, over T = 1 to 1000 K:

| T (K) | 1 | 10 | 100 | 200 | 400 | 1000 |
|---|---|---|---|---|---|---|
| Direct | 5.36 | 4.21 | 4.41 | 4.50 | 4.83 | 5.99 |
| Total | 15.4 | 11.7 | 8.34 | 7.49 | 7.16 | 7.71 |

The "Direct" column is their process (2), radiative charge transfer, and it is
what the code's B2 carries: Stancil, Lepp &amp; Dalgarno (1998) multiply it by the
0.25 approach-probability factor, and interpolating Direct to 300 K gives
4.67e-15, of which 0.25 is 1.17e-15 - the 1.25e-15 of the coded fit, to the
rounding. The 0.25 applies to **Direct**, not to "Total", which additionally
contains radiative association to HeH<sup>+</sup>.

The **non-radiative** channel appears in that paper only as cross sections
over 1-100 eV (their Fig. 8, "Comparison between nonradiative (direct) charge
transfer cross sections (circles) and radiative charge transfer cross sections
(triangles)"), which is the collision-energy range that maps onto the
6e3-1e5 K validity range of the Kingdon &amp; Ferland (1996) fit. The two are
therefore different parts of the same calculation and are **additive**; the
code now sums them. Zygelman et al. also supply the physical reason it matters
in a wind: the radiative processes dominate at low collision energy, while
above a few eV the direct process becomes the faster removal mechanism, the
two cross sections being comparable at about 1e-20 cm<sup>2</sup> in the 5-8 eV range.

Measured sizes of the added channel relative to the radiative one: **0.36x at
1140 K, 1.16x at 3020 K, 2.99x at 10<sup>4</sup> K**, so total He<sup>+</sup> + H removal rises by
1.36x, 2.16x and 3.99x. Note this is *not* uniformly the "4x" that the 10<sup>4</sup> K
figure alone suggests - at the molecular base the correction is under 40%.

**Extrapolation policy, and its precedent.** The Kingdon &amp; Ferland fit is valid
over 6e3-1e5 K, while the molecular base of interest sits near 1140 K. The code
evaluates it as fitted below the range - it decays as t<sub>4</sub><sup>2.06</sup> there, so
downward extrapolation under-weights the channel and cannot make it spuriously
large - and caps the temperature at the published ceiling above it, where the
fit grows without bound. Clamping to the 6000 K value instead was rejected: it
would freeze the rate at 5.1e-15 and overstate the channel at the base by about
a factor 8. Ziegler (2018) section 2.1 does the same kind of thing to a
coefficient from the same Kingdon &amp; Ferland table, requiring that rate
coefficients not show "unphysical behavior or unboundedness ... in the
asymptotic limits" and, for their C<sup>+</sup> + H rate, setting it to zero below the
5390 K at which it turns negative and imposing an upper floor at 1e9 K. The
precedent is for the practice, not the specific action: that fit goes negative
and is zeroed, whereas this one stays positive and is evaluated.

### 4.4 The 4000 K step in the He(2<sup>3</sup>S) + H rate (fixed 2026-08-28)

Taylor et al.'s two-branch fit is discontinuous at its own break point: the
low branch gives 1.5849e-9 at 4000 K and the high branch 2.4921e-9, a
factor 1.5724.  The code reproduced the published fit exactly, so this was
never a transcription error, but a rate coefficient is a continuous
function of temperature and the step is not physical.  It was also visible:
every cell whose temperature crossed 4000 K carried a discontinuity in the
He 2<sup>3</sup>S density, and those were the steps seen in the LHS 1140 b memo
figures.

Two further problems with that fit came out of checking it against its own
paper.  Figure 19 of Taylor et al. plots their Maxwell-Boltzmann average of
the Cohen &amp; Lane (1971) and Morgner &amp; Niehaus (1979) cross sections: it
rises to 1.50e-9 near 1400 K and falls to 1.29e-9 at 4000 K, while the
tabulated low branch decreases monotonically from 1.83e-9 to 1.58e-9 and
cannot reproduce a maximum at all.  And above 4000 K the high branch
exceeds every cross-section determination collected in Garc&iacute;a Mu&ntilde;oz
(2025) Fig. 5.  Figure 19 itself agrees with the Garc&iacute;a Mu&ntilde;oz curve
to 10-20%, so the two independent calculations are consistent and it was
the published fit that was the outlier.

The rate is now the Garc&iacute;a Mu&ntilde;oz (2025) closed form, continuous
everywhere and stated valid to better than 2% from 200 to 10<sup>4</sup> K.

The earlier claim here that the step entered the Jacobian was wrong and is
withdrawn: the coefficient is constant within a cell, so it does not enter
the density Jacobian, and `T_equation` holds the heating fixed while
varying T, so it does not enter the temperature root-find either.  The step
appeared only as a jump in the loss rate between neighbouring cells.

## 4b. Does this move the LHS 1140 b results?

No, and not by a small margin: **no LHS 1140 b run in the repository has
molecular chemistry on.**  None of the 170 `input.inp` files under
`LHS1140b/` sets `Molecular chemistry`, and none of the 193
`Ion_species.txt` files carries an H2 / H2+ / H3+ / HeH+ column.  Every
statement corrected above lives inside a `thereis_mol` branch, so the He/H =
2.09 solution, the `heh{3,5,8,9p7,10p3,12}` closure ladder and the
`heh2p13_diff_kzz1e9` run are untouched, and no re-convergence is needed.

What can be measured on those stored solutions is the size of the factor
that *would* have been wrong had they been molecular.  The mis-supplied
density was &rho;/m<sub>H</sub>; the right one is the EOS particle count
p/(k<sub>B</sub>T) = n<sub>tot</sub> + n<sub>e</sub>.  Their ratio at the base cell,
read from the stored `Hydro_ioniz.txt` of each run:

| stored solution | He/H | (&rho;/m<sub>H</sub>)/(n<sub>tot</sub>+n<sub>e</sub>) at the base |
|---|---|---|
| `LHS1140b/exhale/flux_closure/hi/k06` | 2.09 | 3.03 |
| `LHS1140b/exhale/heh2p13_diff_kzz1e9` | 2.13 | 3.04 |
| `LHS1140b/exhale/flux_closure/heh9p7/k03` | 9.7 | 3.72 |
| `backup/regression/golden/mol_diffusion` | 0.0793 | 2.27 |

So a molecular LHS 1140 b base would have had its three-body reactions
R12/R13/R15 evaluated with a third-body density 3.0-3.7 times too large -
R13 and R15 linearly, R12 linearly - against 2.3 times at the hot-Uranus
base that the regression case covers.  That is the quantitative form of the
statement that the error grew along the composition axis Phase C opened.

## 4c. Golden movement

`make check` splits exactly along the `thereis_mol` branch: `wasp_full`,
`wasp_he23off` and `lower_profile` PASS byte-identical; the four molecular
cases move.  Largest relative difference over the column, new output against
the previous golden, on the 12000-step snapshots:

| case | rho | p | T | heat | cool | H2 front (f = 0.5) | log10 Mdot |
|---|---|---|---|---|---|---|---|
| `mol_base_handoff` | 1.4% | 2.6% | 1.2% | 1.6% | 3.0% | 1.1597 &rarr; 1.1597 | 10.58 &rarr; 10.58 |
| `mol_metals` | 3.0% | 5.0% | 2.2% | 3.5% | 6.7% | 1.1617 &rarr; 1.1617 | 10.58 &rarr; 10.58 |
| `mol_lyman_werner` | 0.34% | 0.53% | 0.21% | 0.44% | 0.71% | 1.1304 &rarr; 1.1304 | 10.58 &rarr; 10.58 |
| `mol_diffusion` | 1.4% | 2.6% | 1.2% | 1.6% | 3.0% | 1.1597 &rarr; 1.1597 | 10.58 &rarr; 10.58 |

The velocity column moves by more (2.5-8x), but only in the base sound-wave
layer of these unconverged snapshots, where v is small and oscillating.

The H2 front does not move at all and the base H2 fraction shifts in the
fourth decimal.  That is not a small effect hiding: these are H-rich
hot-Uranus bases already at f_H2 ~ 1, and a fraction pinned at unity cannot
respond to a 2.3x change in the third-body density -- it responds in rho, p
and T instead.  It is the same statement the He-rich scan of section 3 makes
from the other side, where the base *is* free to move and does.

Goldens were re-snapshotted for the four moved cases at the end of the
change series (`run_check.sh golden mol_base_handoff mol_metals
mol_lyman_werner mol_diffusion`), and `make check` re-run afterwards.

## 5. Files changed

| File | Change |
|---|---|
| `src/modules/radiation/ionization_equilibrium.f90` | `n_tot` from `calc_ntot` replaces `n_in_dim` as the third body M and as the seed pressure (*the R13/R15 half of this, the argument of `set_mol_coeffs`, was missed and corrected on 2026-08-29; see section 7*); HeH<sup>+</sup> nucleus added to the He row of `ionization_fractions_physical`, `element_budget_violation`, `clamp_fractions_to_element_budget` |
| `src/modules/nonlinear_system_solver/ion_cell_state.f90` | `ntot` documented as the electron-free gas-particle density, not &rho;/m<sub>H</sub> |
| `src/modules/nonlinear_system_solver/System_HeH_mol.f90` | header: the H-He exchange is carried (Huang B1/B2), R21/R22 unused; params slot 21 described correctly |
| `src/modules/lower_atmosphere/mol_rates.f90` | header notes on the third body M and the He-dominated limit; per-reaction notes on R12/R13/R15, R16, R17, R21/R22 |
| `src/modules/radiation/Cool_coeff.f90` | both Penning fits: the tabulated value is the Penning + associative total, the 0.9:0.1 split is from the Appendix (not Fig. 4), and the 4000 K step is the published fit. *Superseded 2026-08-28*: both rates replaced by the Garc&iacute;a Mu&ntilde;oz (2025) continuous forms, renamed `ioniz_HeI23S_H` / `ioniz_HeI23S_H2` because they return the total, and the 0.9:0.1 split applied through `f_penning_HeI23S` |

## 6. What was not checked

* The photo-rates P1-P5 of Table 1 (they are "SC" in the paper and follow
  EXHALE's own cross sections, audited elsewhere).
* R1-R4, which duplicate rates EXHALE takes from its own atomic modules;
  the molecular system uses the EXHALE arrays, not `mol_rates`, for them.
* Convergence.  Every run here is a 12000-step relaxation snapshot with no
  Newton finish, chosen so that the He/H = 0.0793 run reproduces a golden
  bitwise.  No `Mdot` in this document is a converged mass-loss rate.
* ~~Whether the row scaling of the molecular system is what drives the
  simplex failures at high He/H.~~  **Answered 2026-08-29**, see section 7.
* ~~The He/H = 1 NaN abort.~~  **Answered 2026-08-29**, see section 7.

## 7. The two open items, closed (2026-08-29)

Both items section 6 left open were taken up on 2026-08-29 and are recorded in
full in `docs/Update_EXHALE_stage1.md` section 93. In short:

**They are not the same problem.** The row scaling is a solver defect and the
NaN abort is a hydrodynamic one; the abort happens at the same step and the
same face with the scaling on and off.

**Row scaling is now applied**, and it is applied per row to the row's own
turnover rate -- the rate at which the species that row balances is produced or
destroyed in the cell, built from the cell's `n_H`, `n_He`, `n_e` and `n_tot`
and the rate coefficients, and held constant across the cell's solve so the
finite-difference Jacobian still sees a smooth function. The fraction of
`hybrd1` attempts ending on `info = 4` falls from 54.6% to 0.01% at He/H = 1000,
41.3% to 0.01% at 100, and 36.3% to 0.26% at 10. What section 3 suspected -- a
correlation between the row spread and the *clamped cells* -- is not what the
scaling removes: the clamp count is a much rarer population (4e2 against 6e6
cell solves) and it does not fall with the solver failures. The clamp counts a
root that no starting point could put inside the simplex, and with the rows
scaled `hybrd1` reaches its tolerance on nearly every attempt, so a converged
root sitting on a face of the simplex is now rejected where a non-converged one
inside it used to be kept.

**The M fix of section 3 had reached only R12.** R13 and R15 are written in
`mol_rates` as two-body-equivalent coefficients that fold in the third-body
density, so their M enters through the argument of `set_mol_coeffs` -- which
`ionization_equilibrium` was still calling with `n_in_dim`. Corrected. Every
"after the fix" row of the section-3 tables therefore describes a state in which
the H3+ and H2 three-body sources still carried the wrong M; the tables are left
as the record of what was measured, and the scan of section 93 supersedes them.

**The He/H = 1 NaN abort is a negative pressure at face j = 232.** It reappears
once R13/R15 have the right third body (it is absent from the tree as found).
Caught with `-ffpe-trap=invalid,zero,overflow`, the first invalid operation is
`aR = sqrt(g*pR/rhoR)` at `src/modules/flux/Num_Fluxes.f90:46`, at face j = 232,
r = 1.1398 R_p, where the right state has p = -3.6e-2 in code units. Below that
face the run carries a cold hypersonic shell -- 650 K, 1e7 cm/s -- against a
dense hot wall at r = 1.144.

**Corrected 2026-08-30** (`Update_EXHALE_stage1.md` section 95, `TO_BE_DONE.md` items
(O) and (P)): the negative pressure is **not** made by the reconstruction. The
cell average of cell 233 already carries the identical `-3.64985e-2`, its
pressure slope limited to zero, so both of its face states are its own average;
the conservative update takes that cell from p = +5.58 to negative in one step,
in a base region running at Mach 38 to 224. The same run at `CFL: 0.2` completes
12000 steps with no NaN. The repair is a positivity test on each RK stage with a
dt bisection when it fails, plus a positivity guard on the reconstructed face
states -- which was a real and separate gap, firing 1164 times in `wasp_full`
and feeding negative pressures to HLLC, whose NaN sound speed the wave-speed
`min`/`max` was silently discarding.

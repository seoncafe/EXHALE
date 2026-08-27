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
| R12 | H<sub>2</sub> + M &rarr; H + H + M | 1.5e-9 exp(-48350/T), &times; n<sub>M</sub> | same | ok | M identity: sec. 4.1 |
| R13 | H<sup>+</sup> + H<sub>2</sub> + M &rarr; H<sub>3</sub><sup>+</sup> + M | 3.2e-29 n<sub>M</sub> | same | ok | M identity: sec. 4.1 |
| R14 | H<sub>2</sub> + e &rarr; H + H + e | 1.33e-6 (300/T)<sup>0.91</sup> exp(-55800/T) | same | ok | |
| R15 | H + H + M &rarr; H<sub>2</sub> + M | 8e-33 (300/T)<sup>0.6</sup> n<sub>M</sub> | same, but printed with T<sub>e</sub> | ok | a neutral three-body recombination has no electron in it; the heavy-particle T is used, and the two readings coincide in EXHALE |
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

| Channel | Code | Source | Verdict |
|---|---|---|---|
| He(2<sup>3</sup>S) + H | 1.9e-9 (300/T)<sup>0.07</sup> for T &le; 4000 K, 9.1e-9 (300/T)<sup>0.5</sup> above | Taylor et al. (2025) Table 2 | exact transcription; the 1.57x step at 4000 K is the published fit |
| He(2<sup>3</sup>S) + H<sub>2</sub> | 5.3791e-12 T<sup>0.676</sup> exp(-695.21/T) | fit to Garc&iacute;a Mu&ntilde;oz (2025) Table A.5 | reproduces the four tabulated points to 0.13% (re-measured here) |

Both are *totals*: Taylor's Section 2.4 states the fitted cross sections are
the sum of Penning and associative ionization, and Garc&iacute;a Mu&ntilde;oz partitions
the same totals "an average 0.9:0.1" between the two.  EXHALE assigns 100%
to the Penning branch in both cases, so the ~10% that physically ends in
HeH<sup>+</sup> is booked as a proton (H channel) or as H<sub>2</sub><sup>+</sup> (H<sub>2</sub> channel).  In a
He-dominated envelope this is the largest HeH<sup>+</sup> source the network omits
after the one in section 4.2.  The code comment attributing the 0.9:0.1
split to "Fig. 4" is corrected: it is in the Appendix paragraph that
introduces Table A.5.

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
  This was not pursued further here.
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

### 4.3 The H &harr; He charge-exchange pair does not satisfy detailed balance

For H<sup>+</sup> + He &harr; He<sup>+</sup> + H the equilibrium constant is fixed by the
statistical weights and the 10.98 eV (127,500 K) endothermicity:

    k(H+ + He) / k(He+ + H) = (g_He+ g_H)/(g_He g_H+) exp(-127500/T) = 4 exp(-127500/T).

The Huang B1/B2 pair gives 105x that at 10<sup>4</sup> K, 210x at 5000 K and 3500x at
300 K; the Koskinen R22/R21 pair, from the same underlying fits, gives
104x, 198x and 689x.  The two published fits come from different
calculations covering different temperature ranges, and neither pair was
constructed to be mutually consistent.  The endothermic direction is
negligible either way below 10<sup>4</sup> K - at 8000 K, k(H<sup>+</sup> + He) = 1.8e-19
against k(He<sup>+</sup> + H) = 2.8e-15 - but the forward term carries n<sub>He</sub> and the
reverse carries n<sub>H</sub>, so its weight in the He<sup>+</sup> balance scales as
n<sub>He</sub>/n<sub>H</sub>: a 10<sup>2</sup> inconsistency at 10<sup>4</sup> K is amplified by 10<sup>3</sup> at He/H = 10<sup>3</sup>.
It is recorded here and at the code site rather than "corrected", because
correcting it means overriding one published rate with the detailed-balance
image of the other, which is a decision for whoever needs the number.

### 4.4 The 4000 K step in the He(2<sup>3</sup>S) + H Penning rate

Taylor et al.'s two-branch fit is discontinuous at its own break point: the
low branch gives 1.585e-9 at 4000 K and the high branch 2.492e-9, a factor
1.57.  The code reproduces the published fit exactly.  Any cell that
crosses 4000 K sees that step, which is a small non-smoothness in the
He 2<sup>3</sup>S balance and in the Jacobian; it is not a transcription error.

## 4b. Does this move the LHS 1140 b results?

No, and not by a small margin: **no LHS 1140 b run in the repository has
molecular chemistry on.**  None of the 170 `input.inp` files under
`LHS1140b/` sets `Molecular chemistry`, and none of the 193
`Ion_species.txt` files carries an H2 / H2+ / H3+ / HeH+ column.  Every
statement corrected above lives inside a `thereis_mol` branch, so the He/H =
2.09 solution, the `heh{3,5,8,9p7,10p3,12}` closure ladder and the
`heh2p13_diff_kzz1e9` arm are untouched, and no re-convergence is needed.

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
| `src/modules/radiation/ionization_equilibrium.f90` | `n_tot` from `calc_ntot` replaces `n_in_dim` as the third body M and as the seed pressure; HeH<sup>+</sup> nucleus added to the He row of `ionization_fractions_physical`, `element_budget_violation`, `clamp_fractions_to_element_budget` |
| `src/modules/nonlinear_system_solver/ion_cell_state.f90` | `ntot` documented as the electron-free gas-particle density, not &rho;/m<sub>H</sub> |
| `src/modules/nonlinear_system_solver/System_HeH_mol.f90` | header: the H-He exchange is carried (Huang B1/B2), R21/R22 unused; params slot 21 described correctly |
| `src/modules/lower_atmosphere/mol_rates.f90` | header notes on the third body M and the He-dominated limit; per-reaction notes on R12/R13/R15, R16, R17, R21/R22 |
| `src/modules/radiation/Cool_coeff.f90` | both Penning fits: the tabulated value is the Penning + associative total, the 0.9:0.1 split is from the Appendix (not Fig. 4), and the 4000 K step is the published fit |

## 6. What was not checked

* The photo-rates P1-P5 of Table 1 (they are "SC" in the paper and follow
  EXHALE's own cross sections, audited elsewhere).
* R1-R4, which duplicate rates EXHALE takes from its own atomic modules;
  the molecular system uses the EXHALE arrays, not `mol_rates`, for them.
* Convergence.  Every run here is a 12000-step relaxation snapshot with no
  Newton finish, chosen so that the He/H = 0.0793 arm reproduces a golden
  bitwise.  No `Mdot` in this document is a converged mass-loss rate.
* Whether the row scaling of the molecular system is what drives the
  simplex failures at high He/H.  The correlation is measured; the cause
  is not.
* The He/H = 1 NaN abort.  It survives the fix and is not diagnosed here.

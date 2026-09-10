# A2 reaction audit (milestone M1)

**Status: complete. Nothing wired into the code.** This is the M1 deliverable
of `docs/a2_oxygen_option_design.md` section 7: every reaction of that
document's section 2.4 traced to a publication that was read, its validity
range and quoted uncertainty recorded, element and charge conservation checked
reaction by reaction, an explicit verdict on the reactions the design left
open, and a standalone driver that reproduces a published rate curve for each.
Decision D2 of section 8 also puts the thermodynamic data for the reverse
rates here.

What M1 produced:

- `src/modules/lower_atmosphere/oxygen_rates.f90`: the coefficients, the
  Shomate thermodynamic table, and the detailed-balance routine. It is
  **deliberately absent from `SRC` in the `Makefile`** and no existing module
  was touched; M2 connects it.
- `src/tests/a2_m1/a2_m1_rate_check.f90` and its saved output
  `a2_m1_rate_check.out`: the Fortran driver.
- `src/tests/a2_m1/a2_m1_thermo_xcheck.py` and `a2_m1_thermo_xcheck.out`: the
  thermodynamic cross-check against a second, independent thermodynamic table.

Neither test is part of any build. To regenerate them, from `src/tests/a2_m1/`:

```bash
gfortran -O2 -o a2_m1_rate_check.x \
    ../../modules/lower_atmosphere/oxygen_rates.f90 a2_m1_rate_check.f90
./a2_m1_rate_check.x > a2_m1_rate_check.out
python3 a2_m1_thermo_xcheck.py > a2_m1_thermo_xcheck.out
```

Every number in this document was measured by running those two programs or
read out of a publication opened for this audit. Nothing is quoted from
memory or from a rate database without saying so.

---

## 1. Verdict first

Eleven reactions were audited (O1-O12 of the design, with O11 and O12 resolved
as described in section 6). **Seven are in the minimal set the measured budget
supports: three rates and four photolysis channels; four are excluded, each
with a measured reason, and O11 turned out not to be a separate reaction at
all.** The network conserves H, O, C and charge reaction by reaction: the
driver checks this mechanically and prints `CONSERVATION: PASS`.

Eight source discrepancies were found. Four of them change a rate that this
option will actually use:

1. **Three rates come from a later published evaluation than the design
   assumed.** The design cites Baulch et al. (1992), JPCRD **21**, 411
   for O1, O2 and O10, following the `Ba92` key of `zahnle_earth.yaml`. The
   same body re-evaluated all three in Baulch et al. (2005), JPCRD **34**, 757,
   and the later recommendations are adopted here. The change is 11% for O1,
   **a factor 1.84 for O2 at the HD 189733 b base temperature**, and a factor
   2.4 for O10.
2. **`O(1D) + H2` (O6) had no traceable source, and the adopted value is now
   2.6x smaller than the one the design's measured budget was computed with.**
   `zahnle_earth.yaml` labels its 1.5e-10 with the key `Ba92`; Baulch et al.
   (1992), Baulch et al. (2005) and Tsang & Hampson (1986) are combustion
   evaluations and **none of them contains any O(1D) chemistry at all**. The
   value adopted here is the IUPAC evaluation, 1.1e-10. VULCAN's 2.87e-10,
   cited to Tully (1975), is the 300 K output of a statistical phase-space
   model, read in the published paper, and disagrees with four independent
   room-temperature measurements by a factor 2.6.
3. **The apparent two-decade disagreement on O9 is a third-body-efficiency
   difference, not a disagreement about the reaction.** VULCAN's value *is*
   Baulch's H2O-collider recommendation; Photochem's is a single shock-tube
   measurement valid only above 2790 K. Neither is appropriate to an H2/He
   bath. Resolved in section 4.
4. **The design's `O10` source attribution `Ba92, Li91` does not hold.**
   Lifshitz & Michael (1991) studied the reverse direction and is not in the
   1992 data sheet's reference list; it is cited by the 2005 evaluation.

Four more affect the design document's text rather than a number in use, and
are listed in section 8 with the corrections made to `a2_oxygen_option_design.md`.

**The detailed-balance implementation is validated against a published
reverse rate, not only against itself.** Baulch et al. (1992) recommends the
reverse of O1 as a rate in its own right. Reversing their forward rate with
the NIST-JANAF Shomate data of this module reproduces their reverse rate to
**0.3-6.1% over 300-2500 K** (section 5). That is the strongest test in this
audit and it passes.

---

## 2. The audit table

Units are cgs: two-body `k` in cm^3 molecule^-1 s^-1, three-body `k0` in
cm^6 molecule^-2 s^-1, `T` in K. "Δlog k" is the source's own quoted
reliability. "Share" is the gross forward rate as a fraction of O1's, measured
at the HD 189733 b base (section 3).

| id | reaction | adopted expression | range | Δlog k | source, read for this audit | share | in set? |
|---|---|---|---|---|---|---|---|
| **O1** | `OH + H2 -> H2O + H` | `3.6e-16 T^1.52 exp(-1740/T)` | 250-2500 K | ±0.1 at 250 K rising to ±0.3 at 2500 K | Baulch et al. (2005), JPCRD 34, 757, Table 4.1 and data sheet **p. 1029** | 1 (reference) | **yes** |
| **O2** | `O + H2 -> OH + H` | `6.34e-12 exp(-4000/T) + 1.46e-9 exp(-9650/T)` | 298-3300 K | ±0.2 over the whole range | Baulch et al. (2005), Table 4.1 and data sheet **p. 804** | 2.6e-2 | **yes** |
| **O3** | `H2O + hv -> OH + H` | quantum yield 0.89 / 0.78 / 0.89 / 1.00 on bands B1-B4 | see section 7 | see section 7 | Stief, Payne & Klemm (1975) and Slanger & Black (1982), by way of JPL Publication 19-5 entry **B2, pp. 4-37 to 4-42** | dominant H2O sink | **yes** |
| **O4** | `H2O + hv -> H2 + O(1D)` | 0.11 / 0.10 / 0.11 / 0.00 | as O3 | as O3 | as O3 | - | **yes** |
| **O5** | `H2O + hv -> O + H + H` | 0.00 / 0.12 / 0.00 / 0.00 | as O3 | as O3 | as O3 | - | **yes** |
| **O6** | `O(1D) + H2 -> OH + H` | `1.1e-10`, temperature-independent | 200-350 K, used beyond it | ±0.1 at 298 K | Atkinson et al. (2004), Atmos. Chem. Phys. 4, 1461, data sheet **I.A2.18, pp. 1508-1509** | 1.9e-2 | **yes** |
| **O7** | `OH + hv -> O + H` | single merged channel, yield 1 | 0.06-264.9 nm grid | not evaluated | Leiden default; Heays, Bosman & van Dishoeck (2017), A&A 602, A105, **sec. 3.1 and Table 1** | - | **yes, with a caveat** |
| **O8** | `O + H + M -> OH + M` | `k0 = 1.3e-29 T^-1` | none stated | uncertainty factor **10** | Tsang & Hampson (1986), JPCRD 15, 1087, entry 5,4, data sheet **p. 1111**, summary **p. 1091** | 7.8e-9 | no |
| **O9** | `H + OH + M -> H2O + M` | `k0(N2) = 6.1e-26 T^-2.0` (Ar 2.3e-26, H2O 3.9e-25) | 300-3000 K | ±0.5 (N2), ±0.3 (Ar), ±0.5 (H2O) | Baulch et al. (1992) Table 3 **p. 428**, data sheets **pp. 496-498**; unchanged in Baulch et al. (2005) **p. 913** | 9.7e-8 | no |
| **O10** | `OH + OH -> H2O + O` | `5.56e-20 T^2.42 exp(+970/T)` | 250-2400 K | ±0.15 | Baulch et al. (2005), Table 4.1 and data sheet **p. 1032** | 5.5e-7 | no |
| **O11** | `O(1D) + M -> O + M` | - | - | - | resolved inside O6's coefficient; see section 6 | ≤5% of O6 | no |
| **O12** | `O(1D) + H2O -> OH + OH` | `2.2e-10`, temperature-independent | 200-350 K | ±0.1 at 298 K | Atkinson et al. (2004), data sheet **I.A2.19, p. 1510** | 2.8e-5 | no |

**Reverse rates are not transcribed.** O1, O2, O8, O9 and O10 run in near
cancellation on a hot base, so an independently transcribed reverse produces an
arbitrary net rather than a small error. All five reverses come from
`rate_from_detailed_balance` (section 5), which is the design's decision D2 and
what both reference networks do.

### Conservation

The driver builds the (H, O, C, charge) count of every species and balances
every reaction. Its section 0 output:

```
    reaction                            dH   dO   dC   dq
    O1   OH + H2   -> H2O + H             0    0    0    0
    O2   O  + H2   -> OH  + H             0    0    0    0
    O3   H2O + hv  -> OH  + H             0    0    0    0
    O4   H2O + hv  -> H2  + O(1D)         0    0    0    0
    O5   H2O + hv  -> O   + H + H         0    0    0    0
    O6   O(1D)+ H2 -> OH  + H             0    0    0    0
    O7   OH  + hv  -> O   + H             0    0    0    0
    O8   O + H + M -> OH  + M             0    0    0    0
    O9   H + OH+ M -> H2O + M             0    0    0    0
    O10  OH + OH   -> H2O + O             0    0    0    0
    CONSERVATION: PASS
    O12  O(1D) + H2O -> OH + OH           0    0    0    0
```

CO appears in no reaction: decision D4 carries it as an unreactive oxygen
reservoir, so its C and O move only with the gas. The third body M cancels
between the two sides of O8 and O9. `hv` carries no nuclei and no charge. No
reaction in the set changes any species' charge, which is the statement the
design's section 2.4 makes when it excludes the water-family ions.

---

## 3. Where the include/exclude verdicts come from

The reference state is the HD 189733 b base of the P1 Photochem run with the
NCHO network, `vulcan_work/pc_compare_p1/hd189_toa1e-2/pc_ncho_solution.pkl`,
level 121: P = 0.9386 dyn cm^-2, T = 863.911 K, n = 7.869e12 cm^-3. The
abundances were read out of that file for this audit:

| species | mixing ratio | n [cm^-3] |
|---|---|---|
| H | 1.4741e-01 | 1.160e12 |
| H2 | 7.0132e-01 | 5.519e12 |
| He | 1.5027e-01 | 1.183e12 |
| H2O | 5.1164e-04 | 4.026e09 |
| CO | 4.2353e-04 | 3.333e09 |
| OH | 2.4664e-07 | 1.941e06 |
| O | 1.0848e-07 | 8.537e05 |
| O(1D) | 6.0386e-11 | 4.752e02 |

Gross forward rates with the adopted coefficients (driver section 5):

| channel | rate [cm^-3 s^-1] | / O1 |
|---|---|---|
| O1 `OH + H2 -> H2O + H` | 1.4957e7 | 1 |
| O1 reverse `H2O + H -> OH + H2` | 4.4756e6 | 2.99e-1 |
| O2 `O + H2 -> OH + H` | 3.8822e5 | 2.60e-2 |
| O2 reverse `OH + H -> O + H2` | 2.4812e5 | 1.66e-2 |
| O6 `O(1D) + H2 -> OH + H` | 2.8846e5 | 1.93e-2 |
| O12 `O(1D) + H2O -> OH + OH` | 4.2089e2 | 2.81e-5 |
| O10 `OH + OH -> H2O + O` | 8.221 | 5.50e-7 |
| O9 `H + OH + M -> H2O + M` | 1.448 | 9.68e-8 |
| O8 `O + H + M -> OH + M` | 0.1173 | 7.84e-9 |

O8, O9, O10 and O12 are seven decades or more below O1 and are excluded. They
are transcribed in the module anyway, so that the record is a value with a
verdict rather than an omission: the same treatment `mol_rates.f90` gives
R21/R22.

**Both three-body reactions are deeply in the low-pressure limit here.** The
measured reduced pressure `Pr = k0 n / k_inf` is 1.2e-8 for O8 and 2.6e-9 for
O9, so the falloff form and the high-pressure limits are numerically
irrelevant at the base. That disposes of the fact that neither reference
network's `k_inf` for these reactions comes from the cited source (section 4).

### What the change of evaluation does to the measured budget

The design's section 2.2 reports the HD 189733 b budget as measured from a
Photochem run of the VULCAN NCHO network. Those percentages are a correct
report of that run. With the coefficients adopted here they move:

| channel | with the network rates | with the adopted rates |
|---|---|---|
| O2 `O + H2` gross, as a fraction of O1 | 5.26e-2 | 2.60e-2 |
| O6 `O(1D) + H2` gross, as a fraction of O1 | 5.52e-2 | 1.93e-2 |

**The design's finding that O(1D) beats ground-state O at this level does not
survive the change of O6's source.** With VULCAN's Tully value O6 is 5.5% of
O1 against O2's 5.3%; with the IUPAC value O6 is 1.9% against O2's 2.6%. The
qualitative point the design draws from it (that O(1D) is not lumpable into
ground-state O, being three decades less abundant and still comparable in
effect) stands and is if anything strengthened. The specific claim that it is
"twice as effective" is a property of the rate set used, not of the physics.
Section 8 records the correction made to the design.

---

## 4. The two-network comparison, reaction by reaction

Both reference networks are on this machine and both were read:
`photochem/data/reaction_mechanisms/zahnle_earth.yaml` (`k = A T^b exp(-Ea/T)`,
`Ea` in K) and `EXHALE_v1.00/VULCAN/thermo/NCHO_photo_network.txt`
(`k = A T^B exp(-C/T)`, odd ids forward, even ids the thermodynamic reverse).
Ratios below are from the driver's section 1.

### O1 `OH + H2 -> H2O + H`

| source | expression | range |
|---|---|---|
| **adopted**: Baulch et al. (2005) p. 1029 | `3.6e-16 T^1.52 exp(-1740/T)` | 250-2500 K |
| Baulch et al. (1992) p. 552 = `zahnle_earth.yaml` (as 1.740708e-16) | `1.7e-16 T^1.6 exp(-1660/T)` | 300-2500 K |
| VULCAN id 1, Oldenborg & Loge (1992) | `3.57e-16 T^1.52 exp(-1740/T)` | 250-2580 K |

Measured ratios to the adopted value: Baulch 1992 is 0.973 at 300 K, 0.890 at
864 K, 0.912 at 2500 K; VULCAN is 0.9917 at every temperature.

**The two reference networks are not two readings of one evaluation.** They
are the two successive evaluations of the same reaction by the same body, and
VULCAN's is the later one to 0.8%. The design's remark that "the two `A`
differ by about a factor of 2 with similar effective rates" is correct
arithmetic (2.05x in `A`, under 9% in `k`) but the reason is not a
transcription difference.

The 1992 data sheet states its provenance: the recommendation is Zellner's
(J. Phys. Chem. **83**, 18, 1979).

### O2 `O + H2 -> OH + H`

| source | expression | range | Δlog k |
|---|---|---|---|
| **adopted**: Baulch et al. (2005) p. 804 | `6.34e-12 exp(-4000/T) + 1.46e-9 exp(-9650/T)` | 298-3300 K | ±0.2 |
| Baulch et al. (1992) p. 430 | `8.5e-20 T^2.67 exp(-3163/T)` | 300-2500 K | ±0.5 at 300 K falling to ±0.2 above 500 K |
| both networks | `8.5e-20 T^2.67 exp(-3160/T)` | - | - |

Measured ratios to the adopted value: the 1992 form is 0.897 at 300 K,
**1.836 at 864 K**, 1.750 at 1000 K, 0.885 at 2500 K.

Two findings. First, **both networks carry `exp(-3160/T)` where the published
value is `exp(-3163/T)`**: a 3 K slip worth 0.1% at 1000 K, harmless but a
transcription error in both. Second, the 2005 re-evaluation changed the
functional form and the rate at the base temperature by a factor 1.84, which is
larger than either evaluation's own uncertainty band in the 700-1200 K window.
This is the largest single change this audit makes to a rate in the minimal
set.

The two networks agree with each other to 6.6e-4 in `A` and exactly in `b` and
`Ea`, so the transcription of the 1992 form had an independent check, which is
exactly why neither network caught that the recommendation had been superseded.

### O6 `O(1D) + H2 -> OH + H`

| source | value | range | what it is |
|---|---|---|---|
| **adopted**: Atkinson et al. (2004) I.A2.18 | `1.1e-10` | 200-350 K | mean of four independent measurements, IUPAC evaluated |
| VULCAN id 615, cited to Tully (1975) | `2.87e-10` | 300 K, model computed over 100-2100 K | statistical phase-space model, Table II "Theory (present results)" |
| `zahnle_earth.yaml`, key `Ba92` | `1.5e-10` | none stated | **source not traceable** |

The IUPAC preferred value and its reliability were read from the data sheet:
"k = 1.1×10−10 cm3 molecule−1 s−1 , independent of temperature over the range
200-350 K", Δlog k = ±0.1 at 298 K, and the comment that the recommendation is
the mean of Davidson et al. (1976, 1977), Wine and Ravishankara (1981), Force
and Wiesenfeld (1981) and Talukdar and Ravishankara (1996), "all of which are
in excellent agreement."

Three separate problems, all recorded at the code site:

- **Range.** The evaluation covers 200-350 K; the base is 560-2400 K. The value
  is used outside its range because the same evaluation finds the rate
  temperature-independent where it was measured, which is what a barrierless
  insertion should be. That is an extrapolation and the module says so.
- **The `Ba92` key is wrong.** Baulch et al. (1992), Baulch et al. (2005) and
  Tsang & Hampson (1986) were all searched for O(1D) chemistry and contain
  none; they are combustion evaluations. Where Photochem's 1.5e-10 comes from
  could not be established.
- **Tully's number is a model, not a measurement.** The paper was obtained and
  read. His abstract describes "A statistical model of chemical reaction ...
  applied to collisions of O(1D) atoms with the molecules H2, N2, CO, CO2, N2O,
  O3, and H2O. Rate constants for reaction and deactivation are computed over
  the temperature range 100-2100 K"; the phase-space form is that of his
  Refs. 4 (Robinson & Holbrook) and 5 (Pechukas, Light & Rankin). The 2.87 is
  the "Theory (present results)" column of his Table II, p. 1896, whose caption
  reads "Rate constants for removal of O(1D) at 300 K", in units of
  1e-10 cm^3/sec. So the VULCAN transcription of the digits holds, but two
  things it carries do not. It is a **300 K** value: Tully's Fig. 2 has the H2
  removal rate rising with temperature, so a flat 2.87e-10 to 2100 K flattens
  his own temperature dependence as well as leaving his computed range.
  And the "Experiment" it agrees with, 2.9e-10, is his Ref. 10, the preferred
  values of Hampson et al. (1973), J. Phys. Chem. Ref. Data 2, 267, drawn from
  Cvetanovic's compilation of **relative** rate constants: a 1973 evaluation,
  not one of the four independent determinations of 1976-1996 that the IUPAC
  evaluation averages. Adopting Tully's value would mean preferring a 1975
  statistical model, anchored to a 1973 relative-rate evaluation, to four later
  direct measurements.
- **What Tully does support is the channel.** For hydrogen his model "predicts
  that essentially all (> 99.9%) of the O(1D) removal rate is due to reaction
  producing OH + H. This is in complete agreement with experiment." That is the
  same conclusion as the IUPAC comment used above, from an independent
  direction, and it is why the evaluated total is the right number for the
  reactive channel.

**This is the largest remaining uncertainty in the A2 rate set**: a factor 2.6
on a channel that is 2-6% of the net H2 loss, i.e. a few percent on the answer
the option exists to compute. M2 should carry it as an explicit sensitivity.

### O7 `OH + hv -> O + H`

Not a rate but a branching, and the branching is a merge rather than a
measurement. Heays, Bosman & van Dishoeck (2017) section 3.1 states plainly:

> We generally neglected further division of the photoabsorption cross section
> into decay channels leading to distinct dissociation products … As an
> exception, in Sect. 8.6 we undertake to characterise the photodissociation
> branching of H2O into OH and H products, and NH3 into NH2 and NH.

OH is not among the exceptions. Their Table 1 gives OH a dissociation
threshold of 279 nm and an ab initio cross section from van Dishoeck &
Dalgarno; the 264.9 nm at which the data file ends is a file endpoint, not a
statement in the paper. VULCAN's `OH_branch.csv` says in its own header that
it combines the O(1D) branch below 150 nm with the O(3P) branch above it, and
the PHIDRATES tables put the O(1D) branch at 0.88 at Ly-alpha.

**This matters for A2 specifically.** O(1D) produced by OH photolysis feeds
straight back into O6, and a single merged channel with unit yield erases the
distinction. M2 must either split the OH channel or state that it does not.

### O8 `O + H + M -> OH + M`

The Tsang & Hampson data sheet is unusually explicit about what the number is
worth:

> There are no definitive measurements on the rate of this process. Numbers
> given in the literature are rough estimates. We suggest using
> k(O + H + M) = 1.3x10⁻²⁹T⁻¹cm⁶molecule⁻²s⁻¹. Baulch et al. have summarized
> other estimates. They range from 10⁻³³ to 10⁻³⁰cm⁶molecule⁻²s⁻¹. The
> uncertainty is a factor of 10. This reaction is not very important under
> combustion conditions.

**No temperature range is attached to the recommendation.** Photochem carries
1.29e-29 (0.8% low); VULCAN's id 657 carries 1.30e-29 with a range 300-2500 K
that the source does not give. Their high-pressure limits differ by a factor
of 10 (1.0e-12 against 1.0e-11) and **neither is in Tsang & Hampson**;
Photochem's 1.0e-12 is the same placeholder it puts on `H + H (+M)` and
`O + O (+M)`, so it carries no information. Since `Pr` is 1.2e-8 at the base,
the choice does not matter there.

### O9 `H + OH + M -> H2O + M`: the two-decade disagreement, resolved

Baulch et al. (1992) p. 496-498, carried unchanged into Baulch et al. (2005)
p. 913, recommends three collider-resolved low-pressure limits over
300-3000 K:

| collider | `k0` | Δlog k |
|---|---|---|
| Ar | `2.3e-26 T^-2.0` | ±0.3 |
| N2 | `6.1e-26 T^-2.0` | ±0.5 |
| H2O | `3.9e-25 T^-2.0` | ±0.5 |

**VULCAN's id 659, `k0 = 3.89e-25 T^-2`, is Baulch's H2O-collider value.**

Photochem's `1.050748e-26 T^-2.1` is Javoy et al. (2003), section 3.4, which
has now been read. The abstract lists the recombination as reaction (5),
`H + OH + Ar -> H2O + Ar (2790-3200 K)`; section 3.4 numbers the measured
dissociation (5) and the recombination (-5). Section 3.4 reports

> "The rate constant of the water decomposition reaction has been studied in
> the temperature range 2790-3200 K at total pressure of about 250 kPa and
> using mixtures containing 1200-4500 ppm of H2O diluted in Ar. The rate
> constant of the termolecular recombination reaction between H, OH and Ar as
> collision partner was calculated from H2O dissociation rate constant
> measurements and the equilibrium constant."

with `k_-5 = 3.75e21 T(K)^-2.1 cm^6 mol^-2 s^-1`, quoted at ±25%, and a
water-collider value `k_-6 = 6.75e22 T(K)^-2.1 cm^6 mol^-2 s^-1` at ±30%
obtained from it by assuming a relative collision efficiency of 18 for H2O
with respect to argon. Dividing by `N_A^2 = 3.6266e47` puts the argon value at
`1.034e-26 T^-2.1 cm^6 molecule^-2 s^-1`; the coefficient Photochem carries is
1.6% above that, far inside the paper's own ±25%. So the NIST transcription
(`2003JAV/NAU371-377:10`) is confirmed in coefficient, exponent and the
2790-3200 K range, and it adds one thing the transcription in this document
did not carry: **the collider is argon**, not a generic third body.

That changes the comparison. The two networks are quoting different colliders,
VULCAN H2O and Photochem Ar, on top of a source disagreement. Like for like
against Baulch's own argon value, Javoy is 3.9x low at 300 K and 5.0x low at
3000 K; against the N2 value adopted here it is 10.4x low at 300 K and 13.1x
low at 3000 K, which is what the earlier "about 10 and 13" said. The 66-81x
between the two networks over 300-2500 K is therefore **about 17x of collider
times 3.9-4.9x of source** (Baulch's own H2O/Ar ratio is 17.0 and Javoy's
assumed efficiency is 18) and neither network is using a value appropriate to
an H2/He bath.

This module adopts Baulch's N2 value, N2 being the diatomic collider of the
three and the closest published analogue to H2; Ar is the monatomic analogue
for He, a factor 2.7 below it. Reading Javoy narrows the source disagreement
from a decade to a factor of 4-5 but does not move that choice: the measurement
is argon, and its 2790-3200 K validity range begins well above the molecular
layer this module is written for.

Photochem's high-pressure limit, `2.7e-10 exp(-75/T)`, is Cobos & Troe (1985),
Table I entry (14) `H + OH -> H2O`, which has also now been read. It is one of
26 systems in that survey, so the paper does apply to this reaction, but it
prints no Arrhenius form: it tabulates `k_rec,inf = 2.1e-10` at 300 K and
`2.6e-10` at 2100 K in `cm^3 molecule^-1 s^-1`, and `2.7e-10 exp(-75/T)` is
simply the two-point fit through them (2.10e-10 and 2.61e-10). That is where
the "300-2100 K" of the NIST record (`1985COB/TRO1010-1015:15`) comes from,
the two tabulated temperatures, not a stated validity range. Two details the
transcription loses: the method is the simplified statistical adiabatic channel
model of the paper's Part I, with the looseness parameter fitted (`alpha = 1.0
A^-1`, `beta = 2.1 A^-1`, `alpha/beta = 0.48`), not transition-state theory as
the NIST record labels it; and the 300 K experimental anchor cited in the
table's footnote 14 is the isotope exchange `OH + D -> OD + H` (Margitan,
Kaufman & Anderson 1975; Howard & Smith 1982), not a direct measurement of
`H + OH` association. Baulch recommends no `k_inf` at all.

The reaction stays out of the minimal set (it is 9.7e-8 of O1), so none of
this changes a number the option uses. It is recorded because a promotion of
O9 must not inherit either network's value.

### O10 `OH + OH -> H2O + O`

| source | expression | range | Δlog k |
|---|---|---|---|
| **adopted**: Baulch et al. (2005) p. 1032 | `5.56e-20 T^2.42 exp(+970/T)` | 250-2400 K | ±0.15 |
| Baulch et al. (1992) p. 555 = `zahnle_earth.yaml` (as 2.549944e-15) | `2.5e-15 T^1.14 exp(-50/T)` | 250-2500 K | ±0.2 |

Measured, the 1992 form is 1.013x the 2005 one at 300 K, **2.41x at 864 K** and
1.34x at 2500 K.

**The `Li91` half of the yaml's `ref: Ba92, Li91` does not hold for the value
it labels.** Lifshitz & Michael (1991) is a study of the reverse direction,
`O + H2O -> OH + OH`, and it is not in the 1992 data sheet's reference list; it
is cited by the 2005 evaluation, which is the one that revised the
recommendation the yaml does *not* carry. The 1992 recommendation is taken
there from Ernst, Wagner & Zellner.

### O12 `O(1D) + H2O -> OH + OH`

Photochem's 2.2e-10 is **the IUPAC evaluated value** (Atkinson et al. 2004,
data sheet I.A2.19, p. 1510: `k(298) = 2.2e-10`, Δlog k = ±0.1, 200-350 K),
which the audit traced even though the yaml's `Sa03` key does not resolve in
`photochem/data/bib.bib`. VULCAN's id 619, `1.62e-10 exp(+65/T)` over
235-370 K, is 21% below the evaluation at 300 K.

### A cross-check on the existing network

While the H-atom three-body rates were open, `mol_rates.f90`'s R15
(`H + H + M -> H2 + M`, `8e-33 (300/T)^0.6`, from Ham et al. 1970 by way of
Koskinen et al. 2022 Table 1) was compared with Baulch et al. (1992) Table 3
p. 428, which recommends `k0 = 2.7e-31 T^-0.6` for M = H2 over 100-5000 K with
Δlog k = ±0.5. **The two agree to exactly 1.10 at every temperature.** No
change is proposed; the agreement is recorded because it is a free check on a
coefficient the code already uses.

---

## 5. Reverse rates: the thermodynamic data and its validation

Decision D2 puts every reverse rate on

```
k_rev = k_fwd / K_c ,    K_c = exp(-dG/RT) (P0 / kB T)^dn ,   P0 = 1 bar,
```

with `dn = n_products - n_reactants`. This is `rate_from_detailed_balance` in
`oxygen_rates.f90`, and it is what both reference codes do
(`photochem/src/photochem_common.f90` lines 161-186 and
`photochem_eqns.f90:45`; `VULCAN/thermo/gibbs_text.txt`).

### The table

Shomate coefficients A-G for H, H2, O, OH, H2O and CO, transcribed from the
`species:` block of `zahnle_earth.yaml` and **checked digit by digit against
the NIST Chemistry WebBook (SRD 69) gas-phase thermochemistry pages**, whose
source is Chase (1998), NIST-JANAF Thermochemical Tables, 4th ed., JPCRD
Monograph 9. Every coefficient of H, H2, OH, H2O and CO in the NIST
temperature intervals matches the WebBook exactly. NIST's eighth coefficient
`H` is not carried, because NIST sets `H = dfH(298.15 K)`, so dropping it makes
the Shomate enthalpy the absolute enthalpy a reaction `dG` needs.

Atomic oxygen has **no Shomate table on the WebBook**, so its coefficients were
checked against the WebBook's CODATA entries instead (Cox, Wagman & Medvedev
1984): the module returns `dfH(298.15) = 249.173 kJ/mol` against the tabulated
249.18 ± 0.10, and `S(298.15) = 161.053 J/mol/K` against 161.059 ± 0.003.

### Round-trip residual

`k_fwd -> k_rev -> k_fwd` returns the input to 7.1e-15 or better at every
temperature for O1, O2 and O9 (driver section 3). That only tests the
arithmetic.

### The real test: a published reverse rate

Baulch et al. (1992) Table 1 p. 418 recommends the reverse of O1 as a rate in
its own right: `H + H2O -> OH + H2 = 7.5e-16 T^1.6 exp(-9270/T)`, 300-2500 K,
Δlog k = ±0.2, data sheet p. 504. Reversing their own forward rate with the
Shomate table must reproduce it:

| T [K] | reversed from the 1992 forward | published reverse | ratio |
|---|---|---|---|
| 300 | 2.630e-25 | 2.623e-25 | 1.003 |
| 500 | 1.364e-19 | 1.385e-19 | 0.985 |
| 864 | 8.538e-16 | 8.201e-16 | 1.041 |
| 1000 | 4.697e-15 | 4.458e-15 | 1.054 |
| 1500 | 1.989e-13 | 1.874e-13 | 1.061 |
| 2000 | 1.446e-12 | 1.392e-12 | 1.038 |
| 2500 | 5.059e-12 | 5.028e-12 | 1.006 |

**0.3-6.1% over the whole range.** Reversing the adopted 2005 forward rate
instead gives 1.03-1.19, the extra being the change in the forward rate itself.

The same data sheet prints its own thermodynamic data for the reaction, taken
there from the Sandia Chemkin compilation (Kee, Rupley & Miller 1987):
`dH(298) = -62.9 kJ/mol`, `dS(298) = -10.9 J/K/mol`. The Shomate table gives
**-62.815 kJ/mol** and **-10.837 J/K/mol**: 0.13% and 0.6%. Its
three-constant fit `Kp = 0.113 T^0.0639 exp(7680/T)` is reproduced to
1.10-1.19 over 300-2500 K, the residual being the fit's own coarseness; the
1992 evaluation warns in its section 2.6 that a rate obtained by reversal is
very sensitive to the thermodynamic data, particularly to `dH`, and this is the
size of that sensitivity.

### Which thermodynamic table, and the OH enthalpy

`a2_m1_thermo_xcheck.py` compares the adopted NIST-JANAF Shomate data with
Burcat's NASA-9 polynomials, which VULCAN uses
(`VULCAN/thermo/NASA9/`). Measured:

- H, H2, O, CO agree in enthalpy to better than 0.02 kJ/mol at every
  temperature tested; H2O to 0.003 kJ/mol up to 1500 K, drifting to 0.25 at
  2000 K and 0.62 at 2500 K where the two fits part company in their top
  interval.
- **OH differs by a constant +1.709 kJ/mol.** The NIST-JANAF enthalpy of
  formation is 38.99 kJ/mol; Burcat's is 37.28. The entropies agree to
  0.03 J/mol/K, so it is purely an enthalpy-of-formation difference.

That offset multiplies any reverse rate by `exp(1.709 kJ/mol / RT)` per OH in
the reaction: 1.99 at 300 K, 1.23 at 1000 K, 1.09 at 2000 K, and the square of
that for a reaction with two OH.

**The NIST-JANAF value is kept, and the reason is measured, not assumed.** The
rate evaluations being reversed were themselves referred to JANAF-derived
thermochemistry (Baulch's own `dH(298) = -62.9` for O1 is reproduced by the
JANAF table to 0.13%), so reversing them with a different enthalpy of formation
would introduce an inconsistency that is not in either the rate or the
thermodynamics separately. A run that needs the modern Active Thermochemical
Tables value must change the OH `F` coefficient by -1.71 and say so.

For completeness the O10 pair was also tested: reversing `OH + OH -> H2O + O`
against VULCAN's separately listed `O + H2O -> OH + OH` (id 5). The 1992
forward reverses to within 0.82-1.36 of it and the 2005 forward to within
0.44-1.00. **That test turned out to be weaker than it looks**: the much better
agreement of the 1992 pairing suggests VULCAN's id 5 was itself obtained by
reversing the 1992 recommendation rather than transcribed independently (the
network file gives it no reference). The O1 test above is the one that carries
weight.

---

## 6. Explicit verdicts on O10, O11 and O12

Section 7 of the design asks M1 for "an explicit verdict on O10/O11/O12", and
section 2.4 defines only O1-O11. O12 did not exist. The gap is closed as
follows and the design is corrected (section 8).

**O10 `OH + OH -> H2O + O`: EXCLUDED from the minimal set, transcribed.**
Measured at 5.5e-7 of O1 at the HD 189733 b base, seven decades below the
dominant channel. It is kept in the module because it is the only reaction
whose opposite direction appears explicitly in the other network, which made it
the first candidate for a detailed-balance test; section 5 explains why that
test proved weak. Its source attribution in `zahnle_earth.yaml` is partly
wrong and its recommendation was superseded in 2005, both recorded above.

**O11 `O(1D) + M -> O + M`: EXCLUDED, and the reason is now published rather
than assumed.** The design left this as "to be identified at M1". The IUPAC
data sheet that supplies O6 settles it: the quoted `1.1e-10` is the **total**
`k1 + k2`, where channel 2 is exactly this quenching, `O(1D) + H2 -> O(3P) +
H2`, and the evaluation's comment is that channel 1 "appears to be the dominant
pathway (>95%) for the reaction (Wine and Ravishankara, 1982)". So quenching by
H2 is inside O6's own coefficient, is under 5% of it, and has no separately
evaluated rate to transcribe. Neither reference network carries a quenching
channel for H2 or for He, the two colliders that are 85% of this gas. And O6
alone gives O(1D) a chemical lifetime of **1.6e-3 s** at the HD 189733 b base,
which is shorter than every other time scale in the problem by many decades, so
the local steady state the design proposes for O(1D) is the accurate treatment
and not a shortcut.

**O12: the design's M1 row names a reaction its section 2.4 never defined.**
The audit fills the slot with the reaction that most deserves it,
`O(1D) + H2O -> OH + OH`: it is the second O(1D) sink after O6, it is carried
by both reference networks, and it is the one channel that could have upset the
O(1D) steady state, because H2O is 43-55% of the oxygen and O(1D) is made
inside the water. **EXCLUDED, measured.** At the HD 189733 b base it is
2.8e-5 of O1 and it gives O(1D) a lifetime of 1.13 s against 1.6e-3 s from O6:
700x slower, so it shortens the O(1D) lifetime by 0.15%. Its adopted value
is the IUPAC evaluation, as above.

---

## 7. Photolysis: what the branching table actually is

The design's section 2.6 attributes the H2O branching ratios to "JPL
Evaluation 19 (Burkholder et al. 2020) with Stief et al. (1975) and Slanger &
Black (1982)". The entry was read (JPL Publication 19-5, Section 4B,
**B2. H2O (water)**, pp. 4-37 to 4-42) and the attribution needs restating on
four counts. All four are now written into the module header.

1. **JPL 19-5 does not recommend these quantum yields.** Within entry B2 the
   word "recommend" attaches only to the absorption cross section. The yields
   appear in the "Photolysis Quantum Yield and Product Studies" paragraph as a
   literature summary, and at Ly-alpha the entry sets out two *disagreeing*
   sets: Slanger & Black's 0.78 / 0.10 / 0.12 and Mordaunt et al.'s
   0.64 / 0.11 / 0.11 with a fourth channel `OH(A 2Sigma+) + H` at 0.14,
   without choosing between them.
2. **The Ly-alpha 0.78 is `Phi1 + Phi2`** in Slanger & Black's own numbering,
   read in the published paper: their Eqs. (1)-(4) are `OH(X 2Pi) + H`,
   `OH(A 2Sigma+) + H`, `O(1D) + H2(X 1Sigma_g+)` and `O(3P) + 2H`, and their
   p. 2435 reads "We may thus set yields of 78% for processes 1 + 2, 10% for
   process 3, and 12% for process 4". So the 0.78 is the ground-state channel
   plus the electronically excited one, and writing it as a single `OH + H`
   channel buries the OH(A) production. Slanger & Black put the OH(A) fraction
   inside the 0.78 at 8% ("the known 8% yield for process 2 at 1216 A",
   p. 2436, from their reference 28); Mordaunt et al. split off 0.14. Their 8%
   is a v = 0 fraction only, because OH(A 2Sigma+) is predissociated for
   v >= 1, and they therefore read part of the 12% three-body branch as
   H + OH(A, v >= 1) that then falls apart: "the yields into processes 2 and 4
   should perhaps be summed". The final products are `O(3P) + 2H` on either
   reading, so 0.12 stands as the three-body yield.
   **The three yields are not independent measurements.** Slanger & Black
   adopt the O(1D) yield of Stief et al. as the standard the oxygen signal is
   scaled against ("earlier work had indicated an O(1D) yield of 11% in this
   wavelength region"), obtain the three-body branch relative to it ("If 10% is
   taken as the O(1D) quantum yield, then the yield for process 4 is 12%"), and
   the 78% is what normalization leaves. A revision of the O(1D) yield moves
   all three.
3. **The 1.00 / 0.00 of the long-wavelength band is a rounding** of Stief,
   Payne & Klemm's `>= 0.99` and `<= 0.01` for 145-185 nm. Their 105-145 nm
   values, 0.89 / 0.11, are band averages over a lithium-fluoride flash lamp;
   **they give no single-wavelength Ly-alpha value at all.**
4. **VULCAN's `H2O_branch.csv` is not the same table as a function of
   wavelength**, although its four ratio triplets are the same numbers. Its
   Ly-alpha node sits at 121.0 nm, 0.567 nm short of the line, and VULCAN
   interpolates linearly between nodes, so at the actual line center it uses
   **0.837 / 0.105 / 0.058**: the three-body branch is 0.058 there, not 0.12.
   It also has no flat 1231-1450 A interval, ramping instead from 0.89 to 1.00.
   A third source, the PHIDRATES tables behind Huebner & Mukherjee (2015),
   gives 0.75 / 0.10 / 0.15 flat over 984-1304 A and attributes it to the same
   Slanger & Black paper. **The spread on the Ly-alpha three-body branch across
   these three sources is 0.058 to 0.15.**

The four-interval table the design gives *is* an exact reading of
`photochem/data/xsections/H2O.h5`, whose quantum-yield nodes sit at 105,
120.1, 120.2, 123.0, 123.1, 145.0, 145.1 and 185 nm and are constantly
extrapolated outside them, and whose 105 / 145 / 185 nm nodes line up with
Stief's intervals. That is what the module carries, with the four caveats
above at the code site.

Cross sections, for M2: the H2O grid of `H2O.h5` runs 0.1-230.413 nm with a
photodissociation peak of 2.75515e-17 cm^2 at 111.5 nm; the OH grid runs
0.06-264.9 nm with a peak of 1.36397e-17 cm^2 at 100 nm. VULCAN's
`H2O_cross.csv` covers only 6.2-230.413 nm and peaks 5.7% lower, at
2.598e-17 cm^2. The wavelength ranges quoted in the design are Photochem's
concatenation choices, not statements of any one paper: Heays et al.'s own
Table 7 runs to 193.9 nm, not 192 nm, and Ranjan et al. (2020) measured
H2O at 292 K over 186-230 nm, which is where the 230.413 nm endpoint comes
from.

The channel threshold energies the module computes from its own Shomate table,
for the photolysis heating term of M2:

| channel | dH(298.15 K) | erg | eV |
|---|---|---|---|
| `H2O -> OH + H` | 498.81 kJ/mol | 8.2829e-12 | 5.1698 |
| `H2O -> H2 + O(3P)` | 491.00 kJ/mol | 8.1532e-12 | 5.0888 |
| `H2O -> O + H + H` | 926.99 kJ/mol | 1.5393e-11 | 9.6076 |
| `OH -> O + H` | 428.18 kJ/mol | 7.1102e-12 | 4.4378 |
| O(1D) excitation above O(3P) | 189.50 kJ/mol | 3.1468e-12 | 1.9641 |

These are reaction enthalpies at 298.15 K, not 0 K dissociation energies; the
difference is the 0-298 K enthalpy content of the fragments minus the parent,
of order a few kJ/mol out of 428-927, i.e. under 1%. The Shomate fits do not
reach 0 K, so the approximation is used and stated at the code site. The
O(1D) excitation energy is the difference of the `F` coefficients that
`zahnle_earth.yaml` carries for `O1D` and `O`, and that file's own note on the
entry is "Estimated from thermodynamic data at 298 K and species O": it is a
thermal value referred to the statistically averaged 3P term, not a
spectroscopic `1D2 - 3P2` term difference.

---

## 8. Corrections made to `docs/a2_oxygen_option_design.md`

Eight discrepancies between section 2 of the design and what the sources
actually say. All were corrected in place.

| # | where | what the design said | what the sources say |
|---|---|---|---|
| 1 | §2.4 O1, O2, O10 source column | Baulch et al. (1992) | superseded for all three by Baulch et al. (2005); adopted, sizes in section 4 |
| 2 | §2.4 O6 source column | "VULCAN's `615`, cited to Tully (1975)": with no note that Photochem gives 1.5e-10 | Tully's is an RRKM extrapolation contradicted by four measurements; Photochem's `Ba92` key is untraceable; the IUPAC evaluation is adopted |
| 3 | §2.4 O8 "why it is in" | "the three-body closure of O2 at the dense base" | O8 is `O + H + M -> OH + M` and has nothing to do with O2; also, Tsang & Hampson give it no temperature range and call it a factor-10 estimate |
| 4 | §2.4 O9 | "the one place the two networks disagree by two decades in the low-pressure limit", presented as unresolved | the disagreement is 6.4x of third-body efficiency times ~11x of source; resolved in section 4 |
| 5 | §2.4 O10 source | "Baulch et al. (1992) + Lifshitz & Michael (1991)" | Lifshitz & Michael is the reverse direction and is not in the 1992 data sheet; it is a 2005 reference |
| 6 | §2.4 O3-O5 source, §2.6 | "Branching ratios: JPL Evaluation 19 (Burkholder et al. 2020) with …" | JPL 19-5 does not recommend H2O quantum yields; see section 7 |
| 7 | §2.6 | "`VULCAN/thermo/photo_cross/H2O/H2O_branch.csv`, whose four rows are the same numbers" | the ratio triplets are the same but the wavelength dependence is not, and the difference is largest at Ly-alpha; see section 7 |
| 8 | §2.5 | "five species … of order 30 numbers" | five species over one to four Shomate intervals each is 91 coefficients, 112 with CO |

Two further corrections that are not source discrepancies:

- §7, the M1 row, asks for "an explicit verdict on O10/O11/O12" while §2.4
  defines only O1-O11. §2.4 now carries an O12 row and §7 names O9-O12.
- §2.2's measured O(1D) and ground-O shares are a correct report of a run made
  with the VULCAN NCHO rate set; a note now records that the O(1D)/O ordering
  it draws from them depends on which O6 rate is used (section 3).

---

## 9. What M1 does not settle

- **Every publication this audit rests on has now been read.** Tully (1975),
  for O6's alternative value, and Slanger & Black (1982), for the Ly-alpha
  water photolysis yields, were obtained after the first pass and are quoted
  from the published papers in sections 4 and 7; both database transcriptions
  of the digits held, and reading them settled what the digits are (Tully: a
  computed 300 K value with its own temperature dependence, compared against a
  1973 relative-rate evaluation; Slanger & Black: 0.78 is Phi1 + Phi2 in their
  numbering, and only one of the three yields is independent of Stief et al.).
  Javoy et al. (2003) and Cobos & Troe (1985), for O9's Photochem low- and
  high-pressure limits, were likewise obtained after the first pass and are now
  quoted from the published papers in section 3; both NIST transcriptions held,
  and reading them added the argon collider of the Javoy measurement and the
  tabulated, rather than fitted, form of the Cobos & Troe value.
- **O6's factor 2.6 is the open physics.** It is 2-6% of the net H2 loss on the
  one planet the option exists for, and it should be an M2 sensitivity rather
  than a fixed number.
- **The OH photolysis channel is merged** (section 4, O7). Splitting it changes
  where O(1D) comes from, and O(1D) is why O6 is in the set at all.
- **The FUV band flux keys themselves are M2's subject**, not M1's. This
  document fixes the band edges and the yields on them; the three-key split of
  decision D1 and the energy bookkeeping of gate G4 are not exercised here.
- **The Shomate table is not exercised above 2500 K.** The H2O fit drifts from
  Burcat's by 0.62 kJ/mol at 2500 K and the two were not compared above it.

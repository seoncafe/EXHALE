# P39: where the two Lyman-Werner band constants come from, and what the
# 37 per cent against CLOUDY actually is

*2026-09-02. The subject is `sigma_lw = 3.4520e-18 cm^2` and
`p_diss_lw = 0.135` in `src/modules/lower_atmosphere/lyman_werner.f90`, and
the CLOUDY c25.00 comparison recorded in `docs/h2_self_shielding_cloudy.md`
sec. 8, Table 4. Sections 1-4 are the source check. Sections 5-8 are what a
patched CLOUDY then measured directly, which retracts one number of the first
draft of this note and settles the choice of fix.*

---

## 1. Verdict, stated first

**The two H2 line data sets agree. They are not the source of the gap.**
Rebuilt from CLOUDY's own line files, the band-averaged dissociation cross
section for a thermal (LTE) H2 population at 100 K in a flat-`F_lambda` field
is `3.479e-18 cm^2`, against the `3.4520e-18 cm^2` the module takes from
Draine & Bertoldi (1996). That is **0.8 per cent**, not 37. The pump cross
section agrees even better: `2.5565e-17` against the module's `2.5570e-17`,
**0.02 per cent**. Draine & Bertoldi's line data (Abgrall & Roueff 1989 and
Abgrall et al. 1992, 1993a, 1993b, in the form Roueff provided in 1992) and
CLOUDY's (Abgrall, Roueff & Drira 2000) are the same calculations by the same
group, and for this band average they give the same answer.

**What the 37 per cent is, most likely: line trapping of the fluorescent
decay photons, carried by CLOUDY's converged pass and absent from DB96's
unshielded Table 2.** CLOUDY's dissociation probability per pump is
`H2_dissprob / H2_rad_rate_out`, and `H2_rad_rate_out` accumulates
`A_ul * Ploss` rather than `A_ul` (`source/mole_h2.cpp:1961-1980, 1993`;
`source/mole_h2_etc.cpp:52`). `Ploss` is the line loss probability, one when
the line is thin. On CLOUDY's **first** pass, when no optical depths are yet
known, the face value of the same cross section the note tabulates is
`3.469e-18 cm^2` at 100 K -- 0.5 per cent from the module, and 0.3 per cent
from the reconstruction above. On the **converged** pass it is
`4.733e-18 cm^2`. The whole factor of 1.37 appears between the two passes of
the same run, with the same atomic data.

Two further facts point the same way. At fixed `T = 1300 K` the converged
face value moves by 13 per cent between `n_H = 1e12` and `1e14 cm^-3`
(`5.061`, `4.619`, `4.483 e-18`) while the first-pass value is identical to
five digits (`3.8192e-18`) -- **atomic data cannot depend on density, an
optical depth can**. And the excess is largest exactly where the H2 is
coldest and its absorption is concentrated in fewest lines (1.36 at 100 K,
1.15 at 900 K, 1.35 at 2700 K where the population spreads again).

**So the module's constants are not wrong against the literature, and the
CLOUDY numbers are not wrong either. They are two different quantities.**
The module's `sigma_lw` is the *unshielded, untrapped* cross section; the
face value the note reads off CLOUDY already carries the trapping enhancement
appropriate to a slab with `N(H2) = 5e21 cm^-2` behind it.

**The consequence, and it is the actionable part.** The tabulated shielding
factor is a ratio taken inside each run,
`f_shield(N) = zeta(N)/zeta(0)`, so its normalization is the *trapped* face
value. Multiplying an untrapped `sigma_lw` by it gives
`zeta_code(N) = (sigma_lw / sigma_CLOUDY(0)) zeta_CLOUDY(N)`, i.e. a rate low
by a constant 1/1.37 to 1/1.15 at every column. **The mixture is the
inconsistency, not either constant.** Section 6 lists the ways out.

---

## 2. The DB96 arithmetic, line by line against the published tables

Read from the published paper, `references/Draine_1996ApJ_468_269.pdf`
(ApJ 468, 269; the ADS scan, tables read from the page images).

**Table 1**, "Interstellar Ultraviolet Radiation Fields (lambda > 912 A)",
gives `chi`, `F/chi`, `dln u_nu/dln nu` at 1000 A and `T_color`. Its rows:

| DB96 Table 1 row | chi | F/chi [cm^-2 s^-1] | dln u/dln nu | T_color [1e4 K] |
|---|---|---|---|---|
| `u_nu ~ nu^-2` (eq. [24]) | -- | `1.208e7` | -2 | 2.90 |
| Draine 1978 (eq. [23]) | 1.71 | `1.232e7` | -8.119 | 1.29 |

`F` is defined by their eq. (21) as the photon flux in the 1110-912 A
interval; the band choice is their footnote 4, which reads in part:
"Photons shortward of 912 A are of course absorbed by atomic hydrogen. H2 in
the v = 0 levels has only weak absorptions longward of 1110 A: the longest
wavelength absorption out of the J = 0 and 1 levels is the Lyman 0-0 P(1)
line at 1110.066 A, with an oscillator strength f_lu = 0.00058 ... it is fair
to assume that the bulk of UV pumping of H2 is due to photons in the
1110-912 A interval."

Integrating their eq. (24), `lambda u_lambda = 4e-14 chi lambda_3 erg cm^-3`,
over 912-1110 A and multiplying by `c` gives `F = 1.2084e7` photons cm^-2
s^-1 at `chi = 1`; the same for their eq. (23) gives `1.2313e7`. Both
reproduce the Table 1 entries to better than 0.1 per cent, so the module's
reading of that column is right and `F = c n_phot` is the right sense of it.

**Table 2**, "Pumping and Dissociation Rates for Unshielded H2", is tabulated
per `H2` level `(v, J)` for `u_lambda` from eq. (24) and from eq. (23), both
at `chi = 1`, and closes with three rows labeled `T_r = 50, 100, 200`. Those
rows are the thermal case: sec. 5.1 describes models "where the distribution
of the absorbing H2 over rotation-vibration levels is assumed to be a thermal
distribution with a specified excitation temperature `T_ex` = 100 K", called
"LTE" models there. The `T_r = 100` row reads

| | zeta_pump [s^-1] | `<p_diss>` |
|---|---|---|
| eq. (24), `chi = 1` | `3.09e-10` | 0.135 |
| eq. (23), `chi = 1` | `2.78e-10` | 0.119 |

Every number the module header quotes at `lyman_werner.f90:18-60` is in the
published tables, and the arithmetic checks:

| module step | value | check |
|---|---|---|
| `sigma_lw = 4.17e-11 / 1.208e7` | `3.4520e-18 cm^2` | exact |
| the `4.17e-11` itself (Fig. 1 annotation) | | `zeta_pump x <p_diss> = 3.09e-10 x 0.135 = 4.1715e-11`, so the figure annotation and Table 2 are the same number |
| Draine-field variant `2.78e-10 x 0.119 / 1.232e7` | `2.6852e-18 cm^2` | module says `2.685e-18`; ratio to the flat-`F_lam` value 0.778, i.e. the "22 per cent lower" of the header |
| `sigma_lw_pump = sigma_lw / 0.135` | `2.5570e-17 cm^2` | exact |

**No error was found in the module's transcription of DB96.**

---

## 3. What each side's line list is, from the sources

| | states | data | provenance as printed |
|---|---|---|---|
| DB96 | `X`, and `B 1Sigma_u+`, `C+ 1Pi_u`, `C- 1Pi_u` -- sec. 2.2: "the first three electronically excited states"; 299 bound `X` states with `J <= 29` | "Energy levels E_i, transition probabilities A_ij, and dissociation probabilities p_diss,i have been published by Abgrall & Roueff (1989) and Abgrall et al. (1992, 1993a, 1993b). We used data generously provided by Roueff (1992), covering levels up to J = 29." | DB96 sec. 2.2 |
| CLOUDY c25.00 | `X`, `B`, `B' 1Sigma_u+`, `C+`, `C-`, `D+ 1Pi_u`, `D-` | `data/h2/transprob_*.dat` and `dissprob_*.dat`, six of each | every file carries `#>>refer H2 As Abgrall, H., Roueff, E., & Drira, I. 2000, A&AS, 141, 297-300` |

CLOUDY carries two electronic states DB96 does not, `B'` and `D`. **They do
not enter this band.** Out of `X(v = 0, J = 0)` the `B'-X` system spans
844.8-905.2 A and `D-X` spans 747.9-885.5 A, entirely shortward of the 912 A
Lyman edge; the `B'`/`D` lines that do land in 912-1110 A start on
vibrationally excited `X` levels, which hold nothing at our temperatures.
Dropping them from the reconstruction of sec. 4 changes the dissociation
cross section by 0.00 per cent at 100 K and 0.01 per cent at 2700 K, and the
pump cross section by 0.5 per cent at 2700 K.

The one difference DB96 themselves flag is internal to the `C` state: their
sec. 2.2 notes that the data they used "include dissociation probabilities
p_diss,i for the `C+(v, J)` levels that are considerably larger than the
dissociation probabilities of Stephens & Dalgarno (1972)". That change is
already in DB96's numbers, and Abgrall, Roueff & Drira (2000) is the
published form of the same group's calculation.

---

## 4. The reconstruction: CLOUDY's line data, DB96's definition

`.../scratchpad/p39/lw_from_cloudy_data.py` reads CLOUDY's
`energy_*.dat`, `transprob_*.dat` and `dissprob_*.dat` directly and evaluates,
for a flat-`F_lambda` band between 912 and 1110 A normalized to DB96's
`chi = 1` (`F = 1.208e7` band photons cm^-2 s^-1, i.e. `2.3735e-4` erg cm^-2
s^-1):

    zeta_pump = sum_lines  n_l  (pi e^2 / m_e c) f_lu  N_nu(lambda)
    zeta_diss = sum_lines  ... x  D_u / (D_u + sum_l A_ul)

with `n_l` the LTE population of `X(v, J)` including the ortho/para
statistical weight, `f_lu` from the tabulated `A_ul`, and `D_u` the
tabulated continuum dissociation rate. There is **no line transfer in it**:
every decay is allowed to escape, which is DB96's unshielded definition and
CLOUDY's first pass. The `X` levels are restricted to `E <= 4100 cm^-1`,
CLOUDY's `ENERGY_H2_STAR` threshold (`source/h2.cpp:10`), so that the result
is the same population average as CLOUDY's `Solomon_dissoc_rate_g`.

| T [K] | sigma_pump [cm^2] | sigma_diss [cm^2] | `<p_diss>` | CLOUDY pass 1, face | ratio |
|---|---|---|---|---|---|
| 100 | 2.5565e-17 | 3.4791e-18 | 0.1361 | 3.4694e-18 | 1.003 |
| 300 | 2.5642e-17 | 3.5259e-18 | 0.1375 | 3.5051e-18 | 1.006 |
| 900 | 2.5987e-17 | 3.7289e-18 | 0.1435 | 3.6979e-18 | 1.008 |
| 1000 | 2.6045e-17 | 3.7610e-18 | 0.1444 | 3.7336e-18 | 1.007 |
| 1300 | 2.6211e-17 | 3.8497e-18 | 0.1469 | 3.8192e-18 | 1.008 |
| 1800 | 2.6438e-17 | 3.9681e-18 | 0.1501 | 3.9406e-18 | 1.007 |
| 2700 | 2.6706e-17 | 4.1062e-18 | 0.1538 | 4.0762e-18 | 1.007 |

(The CLOUDY column is `Shield(H2) x G(TH85)` in the first zone of the **first**
iteration of the runs of `docs/h2_self_shielding_cloudy.md`, divided by the
same `1.751e13` band photon flux that note used. The last column is a uniform
0.7 per cent offset, which is the precision of the three-digit `save h2 rates`
output plus the small band-edge and non-LTE differences.)

Against DB96's Table 2 at `T_r = 100 K`, the reconstruction gives
`zeta_pump = 3.088e-10 s^-1` where DB96 print `3.09e-10`, and
`<p_diss> = 0.1361` where DB96 print `0.135`.

**Both sides therefore agree, at the same definition, to under one per
cent.** The DB96 constant in the module is the 100 K row of this same table;
what the module lacks is only the mild temperature dependence, `+11 per cent`
in `sigma_diss` and `+13 per cent` in `<p_diss>` between 100 and 2700 K.

### 4.1 A byproduct: the 0.4 eV fragment energy is confirmed

`e_lw_fragment_erg = 6.40871e-13` (0.4 eV, Black & Dalgarno 1977) is the
kinetic energy given to the `H + H` pair. Weighting the mean dissociation
kinetic energy tabulated alongside each dissociation probability in CLOUDY's
`dissprob_*.dat` by the same rate weights gives **0.397 eV at 100 K, 0.406 eV
at 1300 K and 0.429 eV at 2700 K**. The module's constant is right to 2 per
cent over the whole range of the layer, and its temperature dependence is
weaker than that of `sigma_lw` itself.

---

## 5. The depth dependence, measured rather than inferred

**CLOUDY does not write the single-pump branching out**, because its Solomon
rate carries the escape probability inside it. Four columns were added to
`save h2 rates` in a scratch build of c25.00 -- the fresh pump rate of the H2g
population, the dissociation rate that pump gives with `D/(D + sum A_ul)`, the
rate it gives with CLOUDY's `D/(D + sum A_ul Ploss)`, and the band field, all
at full precision. The addition is inert: the 27 decks re-run with it
reproduce the unpatched `f_shield` curve with a maximum relative difference of
**exactly zero**, and the table built from them is byte-identical to the one
in the tree. The patch is `Solomon_pump_rate_g` / `Solomon_dissoc_single_g` in
`H2_Solomon_rate` (`mole_h2_etc.cpp`) plus two `fprintf` fields in
`mole_h2_io.cpp`.

**The attribution of section 1 is now a measurement, not an inference.** At
the illuminated face of the 1300 K run the single-pump branching is **0.1466
on both passes** -- it cannot move, there is no escape probability in it --
while the effective one goes from 0.1466 on the first pass to 0.1771 on the
converged one. The 1.37 of section 1 is that step and nothing else.

The two, through the 1300 K, `n_H = 1e13 cm^-3` run:

| N(H2) [cm^-2] | p_single | p_eff | p_single / 0.135 |
|---|---|---|---|
| face (2.4e10) | 0.147 | 0.177 | 1.09 |
| 1e14 | 0.151 | 0.184 | 1.12 |
| 1e16 | 0.219 | 0.278 | 1.62 |
| 1e17 | 0.226 | 0.293 | 1.67 |
| 1e18 | 0.206 | 0.275 | 1.53 |
| 1e19 | 0.158 | 0.219 | 1.17 |
| 1e20 | 0.144 | 0.209 | 1.07 |
| 5e20 | 0.148 | 0.220 | 1.09 |
| 1e21 | 0.149 | 0.226 | 1.10 |
| 5e21 | 0.151 | 0.215 | 1.12 |

Two readings.

* **The depth dependence of `p_single` is real and it is not trapping.**
  Self-shielding removes the strongest pumping lines first and the lines that
  survive to depth have a different branching; the excursion peaks near
  `N(H2) = 1e17` and comes back. `p_eff` carries that shape plus the trapping
  on top of it.
* **`docs/h2_self_shielding_cloudy.md` sec. 8 understates the excursion of
  `p_eff`.** It reports "0.177 at the face to 0.217 at `N(H2) = 5e21`", the
  two endpoints; the profile peaks near **0.29 at `N(H2) ~ 1e17`**.

### 5.1 Which p the heat return needs -- a retraction

**The first draft of this note said the module overstates the vibrational heat
return by a factor 1.8-2.5 by using 0.135 where CLOUDY gives 0.28. That is
wrong and is withdrawn.** The comparison was against the wrong quantity.

The heat return counts fluorescent decays per dissociation. Take one absorption:
it dissociates with probability `p_single` and otherwise fluoresces. When the
layer is thick the fluorescent photon can be re-absorbed, and the chain
continues; per fresh pump the chain ends in dissociation with probability
`p_eff = D/(D + A Ploss)` and contains `(D + A)/(D + A Ploss)` pump events. So

    pumps per dissociation = [(D + A)/(D + A Ploss)] / [D/(D + A Ploss)]
                           = (D + A)/D = 1/p_single ,

**independent of Ploss**, and the decays per dissociation are
`(1 - p_single)/p_single`. The trapping cancels, and it must: a trapped decay
plus the re-absorption that follows it moves one molecule down and another up,
depositing nothing. CLOUDY's own bookkeeping says the same thing -- the return
flux into `X(v'', J'')` is `Pop_up x A_ul x Ploss` (`mole_h2.cpp:2035-2039`)
and the excited-state population is `rate_in/(D + sum A_ul Ploss)`
(`:1993-1999`), so the arrivals in vibrationally excited X per dissociation
are `A_+/D`, with no escape probability in them.

So `p_diss_lw = 0.135` is compared against `p_single`, not against `p_eff`,
and the correction is the last column of the table above: the heat-return
prefactor `(1 - p)/p` falls from 6.41 to between 3.85 (at `N(H2) = 1e18`) and
5.8 (at the face and in the deepest cells) -- **a factor 1.1 to 1.7, not 1.8
to 2.5**.

**The photon ledger is the other count.** Only fresh pumps take photons out of
the stellar beam, so the number of band photons removed per dissociation is
`1/p_eff` and not `1/p_single`. `write_output.f90` builds `ph_lw` from that
ratio, so it needs `p_eff` where the energy equation needs `p_single`. The two
differ by 20-40 per cent through the layer.

### 5.2 The 2.0 eV, from the two papers and from the line data

**What Burton, Hollenbach & Tielens (1990) actually define.** Their Appendix A
reads: "Vibrational heating occurs by the collisional deexcitation of excited
levels of the H2 molecule, populated by radiative decay of the UV excited
Lyman and Werner bands. The process is approximated as per Paper I, except
that the effective heating per pump, `E_*`, has been taken as 2.0 eV. An
excited pseudolevel at `v = 6`, `n*_H2`, has been taken to represent the
population of UV-excited molecules." So the counting unit is **the pump, not
the dissociation**.

**What the 3.2e-12 erg multiplies.** Neither the pump rate nor the
dissociation rate. Their eq. (A1) is

    H_vib = 3.2e-12 [gamma_(1->0),H n_H + gamma_(1->0),H2 n_H2] n*_H2 ,

i.e. `E_*` times the **collisional de-excitation rate of the UV-excited
pseudolevel population**. Their eq. (2) puts that population in a steady state
against radiative decay `A`, re-pumping `P` and dissociation, and the text
takes the limit explicitly: "when `n >> n_crit` and `n/G_0 >> 1`, collisional
deexcitation dominates and nearly every UV pump into excited `v` states
(typically ~2 eV above the ground state) leads to the transformation of the
vibrational energy into heat." That is the structure `ionization_equilibrium`
has: a rate times `E_*` times `h2_vibrational_heat_fraction`, the last factor
being the quenched share their `n gamma/(A + n gamma)` supplies.

**So `E_*` is per pump that lands in a bound excited level**, with 2.0 eV the
mean energy of the level reached -- not an average over all pumps that counts
the dissociating ones as zeros. Black & van Dishoeck (1987), p. 412-413, fix
that split in the same terms: "The initial absorptions are followed by
fluorescence to the vibrational continuum of the ground state with a typical
average probability of 0.10 ... The remaining 90% of the fluorescent
transitions populate various bound excited vibration-rotation levels of the
ground state."

**Which makes `(1 - p_single)/p_single x 2.0 eV` the right form, and the line
data say so to 2-6 per cent.** Summing CLOUDY's Abgrall, Roueff & Drira (2000)
`A_ul` over every decay of every in-band pump, weighted by the `X(v, J)`
energy the decay lands on, gives the internal energy left behind per pump
directly:

| T [K] | p_single | `<E_int>` per pump [eV] | `2.0 x (1 - p_single)` [eV] | share of decays landing on `v'' = 0` |
|---|---|---|---|---|
| 100 | 0.136 | 1.764 | 1.728 | 0.162 |
| 900 | 0.144 | 1.772 | 1.713 | 0.163 |
| 1300 | 0.147 | 1.779 | 1.706 | 0.163 |
| 2700 | 0.154 | 1.795 | 1.693 | 0.165 |

The module's form is **2 per cent low at 100 K and 6 per cent low at 2700 K**
-- the direction matters: the sixth of the decays that land back on `v'' = 0`
and leave only rotational energy is *already inside* the measured 1.76-1.80 eV,
so it must NOT be subtracted a second time. **No correction factor is added or
removed.**

### 5.3 The adopted 2.0 eV is replaced by the computed value

Burton et al.'s 2.0 eV is a number they adopt in one sentence of an appendix,
without a derivation; the sum above is the same quantity computed from the
line data. So the constant goes and the computation stays, as
`h2_energy_per_bound_fluorescence_eV(T)` in
`src/modules/lower_atmosphere/h2_vibrational_relaxation.f90`.

**Which of the two forms is stored, and why.** Per PUMP the sum is
1.769-1.799 eV over 700-3200 K; per BOUND FLUORESCENCE -- dividing by
`1 - p_single` -- it is 2.060-2.128 eV. The second is what the module carries,
because it is the quantity that belongs to the fluorescence cascade alone and
not to the branching of whatever excited the molecule: a Lyman-Werner photon
brings its own `p_single`, and a photoelectron brings a different, state
resolved one (0.27 for B, 0.036 for C, `electron_energy_degradation.f90`).
Both callers therefore keep their own `(1 - p_diss)` factor and share the
energy, which is how the adopted constant was used and leaves no room for the
factor to be applied twice. For the Lyman-Werner path the product is identical
to using the per-pump form with the pump rate:
`k_diss (1 - p_single)/p_single x E_bound = (k_diss/p_single) x E_pump`.

| T [K] | 700 | 900 | 1300 | 1800 | 2700 | 3200 |
|---|---|---|---|---|---|---|
| `E_bound` [eV] | 2.060 | 2.069 | 2.085 | 2.102 | 2.121 | 2.129 |
| per pump [eV] | 1.769 | 1.772 | 1.779 | 1.786 | 1.795 | 1.799 |
| against 2.0 eV | +3.0% | +3.4% | +4.3% | +5.1% | +6.1% | +6.4% |

A quadratic in `T` reproduces the nine points to 0.03 per cent (a straight
line would do to 0.31 per cent), clamped at both ends of 700-3200 K.

**What the number does not carry.** The weights are those of an unattenuated
band. With the shielded weights of the line model of section 4 the same sum
gives 2.17-2.25 eV at `N(H2) = 1e18-1e21 cm^-2` -- **+4 to +9 per cent**,
because self-shielding removes the strongest lines first and the survivors
land higher. The face value is kept, so the term is a lower bound inside the
shielded layer by about as much as the correction it applies; the two wing
treatments of that model bracket the shifted value within 4 per cent, which is
why it is quoted as a range and not adopted.


---

## 6. What was adopted

**Option C: the table carries the absolute cross section and both branchings.**
The constraint every option had to respect is the band share `A(N)` the H2O
and OH continua see (`util_ion_eq.f90:190, 246`). `A(N)` is evaluated from
Draine & Bertoldi's eq. (39) in closed form (`h2_band_equivalent_width`), and
`sigma_lw_pump = sigma_lw/p_diss_lw` is the differential normalization that
statement rests on rather than an input to it -- so leaving the fit alone
leaves the share alone, and the constant is the check on it. It is the *pump*
cross section, the one quantity here that no line transfer touches: `2.5570e-17 cm^2` in the module against
`2.5565e-17` reconstructed at 100 K and `2.62e-17` at 1300 K. **Its agreement
with CLOUDY is not a cancellation of two errors** -- it holds because
`sigma_lw` and `p_diss_lw` are *both* untrapped DB96 values and their ratio is
trapping-free. It is therefore left exactly as it was, and so is the 0.4 eV
fragment energy, which section 4.1 confirms to 2 per cent.

One consequence to note: with the rate and both branchings on the table,
`sigma_lw`, `p_diss_lw` and `sigma_lw_pump` no longer enter any number the
code produces. They stay in `lyman_werner.f90` as the stated normalization of
the eq. (39) fit that `A(N)` does use, and the module header says that is all
they are.

What changed instead:

| quantity | before | after |
|---|---|---|
| dissociation rate | `sigma_lw x f_shield(N, T, n_H)`, an untrapped constant times a face-normalized factor | `sigma_diss(N, T, n_H)`, the measured cross section per incident band photon, self-shielding and trapping both inside it |
| fluorescence heat return | `(1 - 0.135)/0.135` | `(1 - p_single(N, T, n_H))/p_single(...)` |
| band photon ledger | `k_diss/0.135` | `k_diss/p_eff(N, T, n_H)` |
| band share `A(N)` | `sigma_lw_pump` x DB96 eq. (37) | unchanged |
| fragment energy | 0.4 eV | unchanged |

The mixture was the defect: `f_shield` is a ratio taken inside each run, so its
normalization is the *trapped* face value, and multiplying an untrapped
`sigma_lw` by it gave a rate low by a constant 1/1.37 to 1/1.15 at every
column. Tabulating the absolute cross section removes the question rather than
answering it.

**What the adopted table inherits, and it is stated in its header.** The
trapping inside `sigma_diss` and `p_eff` is computed from escape probabilities
in a plane-parallel slab illuminated on one face and closed on the other. Our
layer has the same thick-inward side, which is the side that matters, but it
is spherical and open outward, and the size of the enhancement depends on the
column behind the point -- a property of the run, not of H2. The measured
sensitivity: at fixed `T = 1300 K` the converged face cross section moves 13
per cent between `n_H = 1e12` and `1e14` while the untrapped value is identical
to five digits. `p_single` is free of this. **This is item (P47) of
`TO_BE_DONE.md`**: the trapping is the one ingredient of the table still read
from the CLOUDY runs, and it is the one carrying a geometry our layer does not
have. (The sentence that stood here, that everything above
`h2_shielding_overlap_column` remains an upper bound, is no longer true:
`Update_EXHALE_stage1.md` section 135 rebuilt the table from a calculation in which
the lines absorb each other's beam and removed that function.)

---

## 7. What was not established

* **`Ploss` was still not read out of a run.** The identification of the
  pass-1 to pass-2 step with line trapping now rests on a direct measurement
  of the two branchings (section 5) rather than on the code path alone, but
  the escape probabilities themselves were not printed. What is measured is
  that the difference between the two branchings is exactly the factor CLOUDY
  puts on `A_ul`, and that the single-pump value is identical on both passes.
* ~~Abgrall, Roueff & Drira (2000), A&AS 141, 297 was not read.~~ **Read
  2026-09-02** (`references/Abgrall_2000AAS_141_297.pdf`); section 3 above is
  still written from the CLOUDY file headers and is to be raised to the
  paper's own section and table numbers when this change reaches the tree.
* **Abgrall & Roueff (1989) and Abgrall et al. (1992, 1993a, 1993b)** were
  likewise not read; DB96 sec. 2.2 is quoted for what they used, and the
  "data generously provided by Roueff (1992)" is by construction not
  available to check.
* ~~Burton, Hollenbach & Tielens (1990) Appendix A was not re-read.~~
  **Closed, section 5.2.**
* The reconstruction of section 4 assumes **LTE populations**. That is what
  DB96's `T_r` rows are and what the CLOUDY decks should be at
  `n_H >= 1e12 cm^-3` (`docs/h2_self_shielding_cloudy.md` sec. 9), but it is
  not what a pumping-dominated PDR has, and none of the numbers here apply to
  that case.
* **The line-by-line model of section 4 is not a substitute for the runs.**
  Extended to shielded columns with Voigt profiles it brackets, rather than
  reproduces, CLOUDY's shielding -- 1.5-3.7 times high with full damping
  wings and about 8 times low with a Doppler core only. Its unshielded limit
  is exact (section 4) and that is all it is used for here.
